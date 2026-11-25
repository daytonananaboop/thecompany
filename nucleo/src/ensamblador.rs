use crate::ensamblaje::{MedidaN50, medir_n50};
use crate::formatos::{LectorFastq, Lectura, abrir_fastq, escribir_fasta};
use crate::secuencia::complemento_inverso;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::hash::{BuildHasher, Hasher};
use std::io::{self, Write};
use std::path::Path;

pub const K_MINIMO: usize = 5;
pub const K_MAXIMO: usize = 31;
pub const COBERTURA_MAXIMA_PEDIDA: usize = 1_000_000;
pub const ULTIMA_BARRA: usize = 255;
pub const MINIMO_SIN_VALLE: u32 = 2;
pub const MAXIMAS_RONDAS: usize = 10;
pub const RAMA_MAXIMA_EN_K: usize = 3;
pub const DIFERENCIA_MINIMA_DE_RAMAS: usize = 3;
pub const MAXIMOS_UNITIGS_EN_DESVIO: usize = 20;
pub const MAXIMAS_VISITAS_POR_RAMA: usize = 500;
pub const MAXIMO_DE_BASES: usize = 300_000_000;
pub const MAXIMO_DE_KMEROS_DISTINTOS: usize = 100_000_000;
pub const LARGO_MAXIMO_DE_LECTURA: usize = 20_000_000;
const BITS_DE_PARTE: u32 = 6;
const PARTES: usize = 1 << BITS_DE_PARTE;
const LECTURAS_POR_BLOQUE: usize = 1024;
const BASES_POR_BLOQUE: usize = 1_000_000;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Kmero {
    pub directo: u64,
    pub reverso: u64,
}

impl Kmero {
    pub fn canonico(self) -> u64 {
        self.directo.min(self.reverso)
    }

    pub fn girado(self) -> Kmero {
        Kmero {
            directo: self.reverso,
            reverso: self.directo,
        }
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct Mezcla;

pub struct Revuelto {
    valor: u64,
}

impl BuildHasher for Mezcla {
    type Hasher = Revuelto;

    fn build_hasher(&self) -> Revuelto {
        Revuelto { valor: 0 }
    }
}

impl Hasher for Revuelto {
    fn write(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            self.valor = mezclar(self.valor ^ byte as u64);
        }
    }

    fn write_u64(&mut self, numero: u64) {
        self.valor = mezclar(self.valor ^ numero);
    }

    fn finish(&self) -> u64 {
        self.valor
    }
}

fn mezclar(numero: u64) -> u64 {
    let mut valor = numero.wrapping_add(0x9e37_79b9_7f4a_7c15);
    valor = (valor ^ (valor >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    valor = (valor ^ (valor >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    valor ^ (valor >> 31)
}

pub fn parte_de(canonico: u64) -> usize {
    (canonico.wrapping_mul(0x9e37_79b9_7f4a_7c15) >> (64 - BITS_DE_PARTE)) as usize
}

pub type Tabla = HashMap<u64, u32, Mezcla>;

pub struct Kmeros {
    pub k: usize,
    pub partes: Vec<Tabla>,
}

impl Kmeros {
    pub fn nuevo(k: usize) -> Kmeros {
        Kmeros {
            k,
            partes: (0..PARTES).map(|_| Tabla::default()).collect(),
        }
    }

    pub fn cantidad(&self) -> usize {
        self.partes.iter().map(|tabla| tabla.len()).sum()
    }

    pub fn cuenta(&self, canonico: u64) -> u32 {
        self.partes[parte_de(canonico)]
            .get(&canonico)
            .copied()
            .unwrap_or(0)
    }

    pub fn contiene(&self, canonico: u64) -> bool {
        self.partes[parte_de(canonico)].contains_key(&canonico)
    }

    pub fn agregar(&mut self, canonico: u64, veces: u32) {
        let cuenta = self.partes[parte_de(canonico)].entry(canonico).or_insert(0);
        *cuenta = cuenta.saturating_add(veces);
    }

    pub fn quitar(&mut self, canonico: u64) {
        self.partes[parte_de(canonico)].remove(&canonico);
    }

    pub fn retener(&mut self, minimo: u32) {
        self.partes
            .par_iter_mut()
            .for_each(|tabla| tabla.retain(|_, cuenta| *cuenta >= minimo));
    }

    pub fn siguientes(&self, kmero: Kmero) -> impl Iterator<Item = Kmero> + '_ {
        (0..4)
            .map(move |base| avanzar(kmero, base, self.k))
            .filter(|siguiente| self.contiene(siguiente.canonico()))
    }

    pub fn anteriores(&self, kmero: Kmero) -> impl Iterator<Item = Kmero> + '_ {
        self.siguientes(kmero.girado()).map(Kmero::girado)
    }
}

pub fn codigo_de_base(letra: u8) -> Option<u64> {
    match letra {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

pub fn letra_de_base(codigo: u64) -> u8 {
    b"ACGT"[(codigo & 3) as usize]
}

pub fn avanzar(kmero: Kmero, base: u64, k: usize) -> Kmero {
    let mascara = (1u64 << (2 * k)) - 1;
    Kmero {
        directo: ((kmero.directo << 2) | base) & mascara,
        reverso: (kmero.reverso >> 2) | ((3 - base) << (2 * (k - 1))),
    }
}

pub fn desde_canonico(canonico: u64, k: usize) -> Kmero {
    let mut kmero = Kmero {
        directo: 0,
        reverso: 0,
    };
    for posicion in (0..k).rev() {
        kmero = avanzar(kmero, (canonico >> (2 * posicion)) & 3, k);
    }
    kmero
}

pub fn letras_de(codigo: u64, k: usize) -> Vec<u8> {
    (0..k)
        .rev()
        .map(|posicion| letra_de_base(codigo >> (2 * posicion)))
        .collect()
}

pub fn recorrer_kmeros(secuencia: &[u8], k: usize, mut usar: impl FnMut(usize, Kmero)) {
    if k == 0 || k > K_MAXIMO {
        return;
    }
    let mut kmero = Kmero {
        directo: 0,
        reverso: 0,
    };
    let mut seguidas = 0;
    for (posicion, &letra) in secuencia.iter().enumerate() {
        let Some(base) = codigo_de_base(letra) else {
            seguidas = 0;
            continue;
        };
        kmero = avanzar(kmero, base, k);
        seguidas += 1;
        if seguidas >= k {
            usar(posicion + 1 - k, kmero);
        }
    }
}

pub fn validar_k(k: usize) -> Result<(), String> {
    if !(K_MINIMO..=K_MAXIMO).contains(&k) || k.is_multiple_of(2) {
        return Err(format!(
            "k debe ser impar y estar entre {K_MINIMO} y {K_MAXIMO}: impar para que ningún k-mero sea su propio complemento inverso, y hasta {K_MAXIMO} para que quepa en 64 bits."
        ));
    }
    Ok(())
}

pub struct Conteo {
    pub kmeros: Kmeros,
    pub lecturas: usize,
    pub bases: usize,
    pub kmeros_leidos: usize,
}

pub fn contar_kmeros(ruta: &Path, k: usize) -> Result<Conteo, String> {
    validar_k(k)?;
    let mut lector = abrir_fastq(ruta)?;
    let mut conteo = Conteo {
        kmeros: Kmeros::nuevo(k),
        lecturas: 0,
        bases: 0,
        kmeros_leidos: 0,
    };
    let mut tanda = leer_tanda(&mut lector)?;
    while !tanda.is_empty() {
        let (siguiente, cubetas) =
            rayon::join(|| leer_tanda(&mut lector), || repartir_tanda(&tanda, k));
        sumar_cubetas(&mut conteo.kmeros, &cubetas);
        for bloque in &tanda {
            conteo.lecturas += bloque.len();
            conteo.bases += bloque
                .iter()
                .map(|lectura| lectura.secuencia.len())
                .sum::<usize>();
        }
        for cubeta in cubetas.iter().flatten() {
            conteo.kmeros_leidos += cubeta.len();
        }
        if conteo.bases > MAXIMO_DE_BASES {
            return Err(format!(
                "El archivo pasa de {} millones de bases; este ensamblador está pensado para genomas bacterianos de hasta unas 60× de cobertura.",
                MAXIMO_DE_BASES / 1_000_000
            ));
        }
        revisar_kmeros_distintos(conteo.kmeros.cantidad())?;
        tanda = siguiente?;
    }
    if conteo.lecturas == 0 {
        return Err("El archivo no tiene ninguna lectura FASTQ.".to_string());
    }
    if conteo.kmeros.cantidad() == 0 {
        return Err(format!(
            "Ninguna lectura tiene {k} letras seguidas de A, C, G o T: no hay k-meros que contar. Prueba con un k menor."
        ));
    }
    Ok(conteo)
}

pub fn revisar_kmeros_distintos(cantidad: usize) -> Result<(), String> {
    if cantidad > MAXIMO_DE_KMEROS_DISTINTOS {
        return Err(format!(
            "Hay más de {} millones de k-meros distintos y no caben en la memoria de un computador de casa: el archivo tiene demasiados errores o no son lecturas de un solo genoma pequeño. Recorta las lecturas por calidad o usa un subconjunto.",
            MAXIMO_DE_KMEROS_DISTINTOS / 1_000_000
        ));
    }
    Ok(())
}

fn leer_tanda(lector: &mut LectorFastq) -> Result<Vec<Vec<Lectura>>, String> {
    let bloques_por_tanda = rayon::current_num_threads().max(1) * 2;
    let mut tanda = Vec::with_capacity(bloques_por_tanda);
    while tanda.len() < bloques_por_tanda {
        let mut bloque = Vec::with_capacity(LECTURAS_POR_BLOQUE);
        let mut bases = 0;
        while bloque.len() < LECTURAS_POR_BLOQUE && bases < BASES_POR_BLOQUE {
            match lector.next() {
                Some(lectura) => {
                    let lectura = lectura?;
                    if lectura.secuencia.len() > LARGO_MAXIMO_DE_LECTURA {
                        return Err(format!(
                            "Hay una lectura de más de {} millones de letras: este ensamblador es para lecturas cortas, no para genomas enteros ni lecturas larguísimas.",
                            LARGO_MAXIMO_DE_LECTURA / 1_000_000
                        ));
                    }
                    bases += lectura.secuencia.len();
                    bloque.push(lectura);
                }
                None => break,
            }
        }
        if bloque.is_empty() {
            break;
        }
        tanda.push(bloque);
    }
    Ok(tanda)
}

fn repartir_tanda(tanda: &[Vec<Lectura>], k: usize) -> Vec<Vec<Vec<u64>>> {
    tanda
        .par_iter()
        .map(|bloque| repartir_bloque(bloque, k))
        .collect()
}

fn repartir_bloque(bloque: &[Lectura], k: usize) -> Vec<Vec<u64>> {
    let mut cubetas = vec![Vec::new(); PARTES];
    for lectura in bloque {
        recorrer_kmeros(&lectura.secuencia, k, |_, kmero| {
            let canonico = kmero.canonico();
            cubetas[parte_de(canonico)].push(canonico);
        });
    }
    cubetas
}

fn sumar_cubetas(kmeros: &mut Kmeros, cubetas: &[Vec<Vec<u64>>]) {
    kmeros
        .partes
        .par_iter_mut()
        .enumerate()
        .for_each(|(parte, tabla)| {
            for bloque in cubetas {
                for &canonico in &bloque[parte] {
                    let cuenta = tabla.entry(canonico).or_insert(0);
                    *cuenta = cuenta.saturating_add(1);
                }
            }
        });
}

pub fn histograma(kmeros: &Kmeros) -> Vec<usize> {
    let mut barras = vec![0; ULTIMA_BARRA + 1];
    for tabla in &kmeros.partes {
        for &cuenta in tabla.values() {
            barras[(cuenta as usize).min(ULTIMA_BARRA)] += 1;
        }
    }
    barras
}

pub fn primer_valle(histograma: &[usize]) -> Option<usize> {
    let ultima = histograma.len().min(ULTIMA_BARRA).saturating_sub(1);
    (1..ultima).find(|&cobertura| histograma[cobertura + 1] > histograma[cobertura])
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Umbral {
    pub minimo: u32,
    pub valle: Option<usize>,
    pub automatico: bool,
}

pub fn elegir_umbral(histograma: &[usize], pedido: usize) -> Umbral {
    let valle = primer_valle(histograma);
    if pedido > 0 {
        return Umbral {
            minimo: pedido.min(COBERTURA_MAXIMA_PEDIDA) as u32,
            valle,
            automatico: false,
        };
    }
    Umbral {
        minimo: valle.map_or(MINIMO_SIN_VALLE, |cobertura| cobertura as u32),
        valle,
        automatico: true,
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Lado {
    pub unitig: usize,
    pub al_reves: bool,
}

impl Lado {
    pub fn opuesto(self) -> Lado {
        Lado {
            unitig: self.unitig,
            al_reves: !self.al_reves,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Unitig {
    pub secuencia: Vec<u8>,
    pub cobertura: f64,
    pub primero: Kmero,
    pub ultimo: Kmero,
    pub al_final: Vec<Lado>,
    pub al_inicio: Vec<Lado>,
}

pub struct GrafoDeUnitigs {
    pub k: usize,
    pub unitigs: Vec<Unitig>,
}

impl GrafoDeUnitigs {
    pub fn salidas(&self, lado: Lado) -> &[Lado] {
        let unitig = &self.unitigs[lado.unitig];
        if lado.al_reves {
            &unitig.al_inicio
        } else {
            &unitig.al_final
        }
    }
}

pub fn compactar(kmeros: &Kmeros) -> GrafoDeUnitigs {
    let mut visitados: HashSet<u64, Mezcla> =
        HashSet::with_capacity_and_hasher(kmeros.cantidad(), Mezcla);
    let mut unitigs = Vec::new();
    for tabla in &kmeros.partes {
        for &semilla in tabla.keys() {
            if visitados.insert(semilla) {
                let camino =
                    caminar_unitig(kmeros, desde_canonico(semilla, kmeros.k), &mut visitados);
                unitigs.push(armar_unitig(kmeros, &camino));
            }
        }
    }
    ordenar_y_conectar(unitigs, kmeros)
}

fn ordenar_y_conectar(mut unitigs: Vec<Unitig>, kmeros: &Kmeros) -> GrafoDeUnitigs {
    unitigs.par_sort_unstable_by(|uno, otro| {
        otro.secuencia
            .len()
            .cmp(&uno.secuencia.len())
            .then_with(|| uno.secuencia.cmp(&otro.secuencia))
    });
    conectar(&mut unitigs, kmeros);
    GrafoDeUnitigs {
        k: kmeros.k,
        unitigs,
    }
}

pub fn recompactar(
    anterior: GrafoDeUnitigs,
    kmeros: &Kmeros,
    quitados: &[usize],
) -> GrafoDeUnitigs {
    let mut fuera = vec![false; anterior.unitigs.len()];
    for &numero in quitados {
        fuera[numero] = true;
    }
    let mut tocados = Vec::new();
    for &numero in quitados {
        for al_reves in [false, true] {
            for vecino in anterior.salidas(Lado {
                unitig: numero,
                al_reves,
            }) {
                if !fuera[vecino.unitig] {
                    tocados.push(vecino.unitig);
                }
            }
        }
    }
    let mut extremos: HashMap<u64, usize, Mezcla> =
        HashMap::with_capacity_and_hasher(2 * anterior.unitigs.len(), Mezcla);
    for (numero, unitig) in anterior.unitigs.iter().enumerate() {
        extremos.insert(unitig.primero.canonico(), numero);
        extremos.insert(unitig.ultimo.canonico(), numero);
    }
    let mut visitados: HashSet<u64, Mezcla> = HashSet::default();
    let mut nuevos = Vec::new();
    for tocado in tocados {
        let semilla = anterior.unitigs[tocado].primero;
        if fuera[tocado] || !visitados.insert(semilla.canonico()) {
            continue;
        }
        let camino = caminar_unitig(kmeros, semilla, &mut visitados);
        for kmero in &camino {
            if let Some(&absorbido) = extremos.get(&kmero.canonico()) {
                fuera[absorbido] = true;
            }
        }
        nuevos.push(armar_unitig(kmeros, &camino));
    }
    let mut unitigs: Vec<Unitig> = anterior
        .unitigs
        .into_iter()
        .zip(fuera)
        .filter(|(_, afuera)| !afuera)
        .map(|(unitig, _)| unitig)
        .collect();
    unitigs.extend(nuevos);
    ordenar_y_conectar(unitigs, kmeros)
}

pub fn caminar_unitig(
    kmeros: &Kmeros,
    semilla: Kmero,
    visitados: &mut HashSet<u64, Mezcla>,
) -> Vec<Kmero> {
    let hacia_adelante = extender(kmeros, semilla, visitados);
    let hacia_atras = extender(kmeros, semilla.girado(), visitados);
    let mut camino: Vec<Kmero> = hacia_atras.into_iter().rev().map(Kmero::girado).collect();
    camino.push(semilla);
    camino.extend(hacia_adelante);
    camino
}

pub fn extender(
    kmeros: &Kmeros,
    inicio: Kmero,
    visitados: &mut HashSet<u64, Mezcla>,
) -> Vec<Kmero> {
    let mut camino = Vec::new();
    let mut actual = inicio;
    while let Some(siguiente) = unico(kmeros.siguientes(actual)) {
        if unico(kmeros.anteriores(siguiente)).is_none() {
            break;
        }
        if !visitados.insert(siguiente.canonico()) {
            break;
        }
        camino.push(siguiente);
        actual = siguiente;
    }
    camino
}

pub fn unico(mut vecinos: impl Iterator<Item = Kmero>) -> Option<Kmero> {
    let primero = vecinos.next()?;
    match vecinos.next() {
        Some(_) => None,
        None => Some(primero),
    }
}

fn armar_unitig(kmeros: &Kmeros, camino: &[Kmero]) -> Unitig {
    let mut secuencia = letras_de(camino[0].directo, kmeros.k);
    let mut suma: u64 = 0;
    for kmero in camino {
        suma += kmeros.cuenta(kmero.canonico()) as u64;
    }
    for kmero in &camino[1..] {
        secuencia.push(letra_de_base(kmero.directo));
    }
    let mut unitig = Unitig {
        secuencia,
        cobertura: suma as f64 / camino.len() as f64,
        primero: camino[0],
        ultimo: camino[camino.len() - 1],
        al_final: Vec::new(),
        al_inicio: Vec::new(),
    };
    orientar(&mut unitig);
    unitig
}

fn orientar(unitig: &mut Unitig) {
    let inverso = complemento_inverso(&unitig.secuencia);
    if inverso < unitig.secuencia {
        unitig.secuencia = inverso;
        let primero = unitig.primero;
        unitig.primero = unitig.ultimo.girado();
        unitig.ultimo = primero.girado();
    }
}

fn conectar(unitigs: &mut [Unitig], kmeros: &Kmeros) {
    let mut extremos: HashMap<u64, Lado, Mezcla> =
        HashMap::with_capacity_and_hasher(2 * unitigs.len(), Mezcla);
    for (numero, unitig) in unitigs.iter().enumerate() {
        extremos.insert(
            unitig.primero.directo,
            Lado {
                unitig: numero,
                al_reves: false,
            },
        );
        extremos.insert(
            unitig.ultimo.reverso,
            Lado {
                unitig: numero,
                al_reves: true,
            },
        );
    }
    unitigs.par_iter_mut().for_each(|unitig| {
        unitig.al_final = kmeros
            .siguientes(unitig.ultimo)
            .filter_map(|vecino| extremos.get(&vecino.directo).copied())
            .collect();
        unitig.al_inicio = kmeros
            .siguientes(unitig.primero.girado())
            .filter_map(|vecino| extremos.get(&vecino.directo).copied())
            .collect();
    });
}

#[derive(Debug, Clone, PartialEq)]
pub struct Punta {
    pub unitig: usize,
    pub secuencia: Vec<u8>,
    pub cobertura: f64,
    pub cobertura_rival: f64,
}

pub fn buscar_puntas(grafo: &GrafoDeUnitigs) -> Vec<Punta> {
    let mut puntas = Vec::new();
    for (numero, unitig) in grafo.unitigs.iter().enumerate() {
        if kmeros_del_unitig(grafo, numero) >= 2 * grafo.k {
            continue;
        }
        let punta = Lado {
            unitig: numero,
            al_reves: !unitig.al_final.is_empty(),
        };
        if !grafo.salidas(punta).is_empty() || grafo.salidas(punta.opuesto()).is_empty() {
            continue;
        }
        let rival = cobertura_rival(grafo, punta);
        if rival > unitig.cobertura {
            puntas.push(Punta {
                unitig: numero,
                secuencia: unitig.secuencia.clone(),
                cobertura: unitig.cobertura,
                cobertura_rival: rival,
            });
        }
    }
    puntas
}

fn cobertura_rival(grafo: &GrafoDeUnitigs, punta: Lado) -> f64 {
    let mut mejor: f64 = 0.0;
    for &atras in grafo.salidas(punta.opuesto()) {
        for &hermano in grafo.salidas(atras.opuesto()) {
            if hermano.unitig != punta.unitig {
                mejor = mejor.max(grafo.unitigs[hermano.unitig].cobertura);
            }
        }
    }
    mejor
}

#[derive(Debug, Clone, PartialEq)]
pub struct Burbuja {
    pub unitig_quitado: usize,
    pub quitada: Vec<u8>,
    pub cobertura_quitada: f64,
    pub conservada: Vec<u8>,
    pub cobertura_conservada: f64,
    pub unitigs_conservados: usize,
}

pub fn buscar_burbujas(grafo: &GrafoDeUnitigs) -> Vec<Burbuja> {
    let mut burbujas = Vec::new();
    for numero in 0..grafo.unitigs.len() {
        let rama = Lado {
            unitig: numero,
            al_reves: false,
        };
        let largo = kmeros_del_unitig(grafo, numero);
        let entradas = grafo.salidas(rama.opuesto());
        let destinos = grafo.salidas(rama);
        let con_lazo = entradas
            .iter()
            .chain(destinos)
            .any(|lado| lado.unitig == numero);
        if largo > RAMA_MAXIMA_EN_K * grafo.k
            || entradas.len() != 1
            || destinos.len() != 1
            || con_lazo
        {
            continue;
        }
        let origen = entradas[0].opuesto();
        let mut busqueda = Busqueda {
            destinos,
            evitar: [numero, origen.unitig],
            largo_buscado: largo,
            visitas: 0,
            mejor: None,
        };
        busqueda.explorar(grafo, origen, &mut Vec::new(), 0, f64::INFINITY);
        let Some((camino, cobertura)) = busqueda.mejor else {
            continue;
        };
        let rama_debil = &grafo.unitigs[numero];
        if cobertura > rama_debil.cobertura {
            burbujas.push(Burbuja {
                unitig_quitado: numero,
                quitada: rama_debil.secuencia.clone(),
                cobertura_quitada: rama_debil.cobertura,
                conservada: deletrear(grafo, &camino),
                cobertura_conservada: cobertura,
                unitigs_conservados: camino.len(),
            });
        }
    }
    burbujas
}

fn kmeros_del_unitig(grafo: &GrafoDeUnitigs, numero: usize) -> usize {
    grafo.unitigs[numero].secuencia.len() + 1 - grafo.k
}

struct Busqueda<'a> {
    destinos: &'a [Lado],
    evitar: [usize; 2],
    largo_buscado: usize,
    visitas: usize,
    mejor: Option<(Vec<Lado>, f64)>,
}

impl Busqueda<'_> {
    fn explorar(
        &mut self,
        grafo: &GrafoDeUnitigs,
        desde: Lado,
        camino: &mut Vec<Lado>,
        largo: usize,
        cobertura: f64,
    ) {
        for &siguiente in grafo.salidas(desde) {
            self.visitas += 1;
            if self.visitas > MAXIMAS_VISITAS_POR_RAMA {
                return;
            }
            if self.destinos.contains(&siguiente) {
                self.anotar(camino, largo, cobertura);
                continue;
            }
            let repetido = camino.iter().any(|lado| lado.unitig == siguiente.unitig);
            if repetido
                || self.evitar.contains(&siguiente.unitig)
                || camino.len() >= MAXIMOS_UNITIGS_EN_DESVIO
            {
                continue;
            }
            let nuevo_largo = largo + kmeros_del_unitig(grafo, siguiente.unitig);
            if nuevo_largo > self.largo_buscado
                && !largos_parecidos(nuevo_largo, self.largo_buscado)
            {
                continue;
            }
            let nueva_cobertura = cobertura.min(grafo.unitigs[siguiente.unitig].cobertura);
            camino.push(siguiente);
            self.explorar(grafo, siguiente, camino, nuevo_largo, nueva_cobertura);
            camino.pop();
        }
    }

    fn anotar(&mut self, camino: &[Lado], largo: usize, cobertura: f64) {
        if camino.is_empty() || !largos_parecidos(largo, self.largo_buscado) {
            return;
        }
        let mejora = match &self.mejor {
            Some((_, anterior)) => cobertura > *anterior,
            None => true,
        };
        if mejora {
            self.mejor = Some((camino.to_vec(), cobertura));
        }
    }
}

fn deletrear(grafo: &GrafoDeUnitigs, camino: &[Lado]) -> Vec<u8> {
    let mut secuencia: Vec<u8> = Vec::new();
    for lado in camino {
        let unitig = &grafo.unitigs[lado.unitig].secuencia;
        let orientado = if lado.al_reves {
            complemento_inverso(unitig)
        } else {
            unitig.clone()
        };
        let solapado = if secuencia.is_empty() { 0 } else { grafo.k - 1 };
        secuencia.extend_from_slice(&orientado[solapado..]);
    }
    secuencia
}

fn tolerancia_de_largo(largo: usize) -> usize {
    DIFERENCIA_MINIMA_DE_RAMAS.max(largo / 10)
}

pub fn largos_parecidos(uno: usize, otro: usize) -> bool {
    uno.abs_diff(otro) <= tolerancia_de_largo(uno.max(otro))
}

pub fn quitar_unitig(kmeros: &mut Kmeros, unitig: &Unitig) {
    let k = kmeros.k;
    recorrer_kmeros(&unitig.secuencia, k, |_, kmero| {
        kmeros.quitar(kmero.canonico())
    });
}

pub struct Limpieza {
    pub grafo: GrafoDeUnitigs,
    pub unitigs_iniciales: usize,
    pub puntas: Vec<Punta>,
    pub burbujas: Vec<Burbuja>,
    pub rondas: usize,
}

pub fn limpiar(kmeros: &mut Kmeros) -> Limpieza {
    let mut grafo = compactar(kmeros);
    let unitigs_iniciales = grafo.unitigs.len();
    let mut puntas = Vec::new();
    let mut burbujas = Vec::new();
    let mut rondas = 0;
    while rondas < MAXIMAS_RONDAS {
        let puntas_de_la_ronda = buscar_puntas(&grafo);
        let burbujas_de_la_ronda = buscar_burbujas(&grafo);
        if puntas_de_la_ronda.is_empty() && burbujas_de_la_ronda.is_empty() {
            break;
        }
        rondas += 1;
        let mut quitados = Vec::new();
        for punta in &puntas_de_la_ronda {
            quitar_unitig(kmeros, &grafo.unitigs[punta.unitig]);
            quitados.push(punta.unitig);
        }
        for burbuja in &burbujas_de_la_ronda {
            quitar_unitig(kmeros, &grafo.unitigs[burbuja.unitig_quitado]);
            quitados.push(burbuja.unitig_quitado);
        }
        grafo = recompactar(grafo, kmeros, &quitados);
        puntas.extend(puntas_de_la_ronda);
        burbujas.extend(burbujas_de_la_ronda);
    }
    Limpieza {
        grafo,
        unitigs_iniciales,
        puntas,
        burbujas,
        rondas,
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Contig {
    pub secuencia: Vec<u8>,
    pub cobertura: f64,
}

pub fn elegir_contigs(grafo: GrafoDeUnitigs) -> (Vec<Contig>, usize) {
    let largo_minimo = 2 * grafo.k;
    let mut contigs = Vec::new();
    let mut cortos = 0;
    for unitig in grafo.unitigs {
        if unitig.secuencia.len() >= largo_minimo {
            contigs.push(Contig {
                secuencia: unitig.secuencia,
                cobertura: unitig.cobertura,
            });
        } else {
            cortos += 1;
        }
    }
    (contigs, cortos)
}

pub struct Ensamblaje {
    pub k: usize,
    pub lecturas: usize,
    pub bases: usize,
    pub kmeros_leidos: usize,
    pub kmeros_distintos: usize,
    pub histograma: Vec<usize>,
    pub umbral: Umbral,
    pub kmeros_solidos: usize,
    pub unitigs_iniciales: usize,
    pub puntas: Vec<Punta>,
    pub burbujas: Vec<Burbuja>,
    pub rondas: usize,
    pub contigs: Vec<Contig>,
    pub cortos: usize,
    pub medida: MedidaN50,
}

pub fn ensamblar(ruta: &Path, k: usize, cobertura_minima: usize) -> Result<Ensamblaje, String> {
    if cobertura_minima > COBERTURA_MAXIMA_PEDIDA {
        return Err(format!(
            "La cobertura mínima debe estar entre 0 (automática) y {COBERTURA_MAXIMA_PEDIDA}."
        ));
    }
    let mut conteo = contar_kmeros(ruta, k)?;
    let kmeros_distintos = conteo.kmeros.cantidad();
    let histograma = histograma(&conteo.kmeros);
    let umbral = elegir_umbral(&histograma, cobertura_minima);
    conteo.kmeros.retener(umbral.minimo);
    let kmeros_solidos = conteo.kmeros.cantidad();
    let limpieza = limpiar(&mut conteo.kmeros);
    let (contigs, cortos) = elegir_contigs(limpieza.grafo);
    let largos: Vec<usize> = contigs
        .iter()
        .map(|contig| contig.secuencia.len())
        .collect();
    let medida = medir_n50(&largos);
    Ok(Ensamblaje {
        k,
        lecturas: conteo.lecturas,
        bases: conteo.bases,
        kmeros_leidos: conteo.kmeros_leidos,
        kmeros_distintos,
        histograma,
        umbral,
        kmeros_solidos,
        unitigs_iniciales: limpieza.unitigs_iniciales,
        puntas: limpieza.puntas,
        burbujas: limpieza.burbujas,
        rondas: limpieza.rondas,
        contigs,
        cortos,
        medida,
    })
}

pub fn largos_sin_limpiar(ruta: &Path, k: usize) -> Result<Vec<usize>, String> {
    let conteo = contar_kmeros(ruta, k)?;
    let grafo = compactar(&conteo.kmeros);
    Ok(grafo
        .unitigs
        .iter()
        .map(|unitig| unitig.secuencia.len())
        .collect())
}

pub fn escribir_contigs(destino: &mut impl Write, contigs: &[Contig]) -> io::Result<()> {
    for (numero, contig) in contigs.iter().enumerate() {
        let nombre = format!(
            "contig_{} largo={} cobertura_kmeros={:.1}",
            numero + 1,
            contig.secuencia.len(),
            contig.cobertura
        );
        escribir_fasta(destino, &nombre, &contig.secuencia)?;
    }
    Ok(())
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct IndiceDeReferencia {
    pub k: usize,
    pub largo: usize,
    pub regiones: Vec<(usize, usize)>,
    pub kmeros: Vec<(u64, u32)>,
}

impl IndiceDeReferencia {
    pub fn posiciones(&self, canonico: u64) -> &[(u64, u32)] {
        let desde = self.kmeros.partition_point(|&(otro, _)| otro < canonico);
        let hasta = self.kmeros.partition_point(|&(otro, _)| otro <= canonico);
        &self.kmeros[desde..hasta]
    }

    pub fn ajenos(&self, secuencia: &[u8]) -> usize {
        let mut ajenos = 0;
        recorrer_kmeros(secuencia, self.k, |_, kmero| {
            if self.posiciones(kmero.canonico()).is_empty() {
                ajenos += 1;
            }
        });
        ajenos
    }
}

pub fn indexar_referencia(
    referencia: &[u8],
    regiones: &[(usize, usize)],
    k: usize,
) -> IndiceDeReferencia {
    let mut indice = IndiceDeReferencia {
        k,
        largo: referencia.len(),
        regiones: Vec::new(),
        kmeros: Vec::new(),
    };
    if validar_k(k).is_err() {
        return indice;
    }
    for &(inicio, fin) in regiones {
        let fin = fin.min(referencia.len());
        if inicio < fin {
            indice.regiones.push((inicio, fin));
            recorrer_kmeros(&referencia[inicio..fin], k, |posicion, kmero| {
                indice
                    .kmeros
                    .push((kmero.canonico(), (inicio + posicion) as u32));
            });
        }
    }
    indice.kmeros.par_sort_unstable();
    indice
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Comparacion {
    pub largo_referencia: usize,
    pub cubiertas: usize,
    pub cubiertas_por_largos: usize,
    pub huecos: Vec<(usize, usize)>,
    pub kmeros_de_contigs: usize,
    pub kmeros_ajenos: usize,
    pub ubicaciones: Vec<Option<(usize, usize)>>,
}

pub fn comparar_con_referencia(
    contigs: &[Contig],
    indice: &IndiceDeReferencia,
    largo_de_contig_largo: usize,
) -> Comparacion {
    let mut comparacion = Comparacion::default();
    if validar_k(indice.k).is_err() {
        return comparacion;
    }
    let k = indice.k;
    let mut cubierto = vec![false; indice.largo];
    let mut cubierto_por_largos = vec![false; indice.largo];
    for contig in contigs {
        let es_largo = contig.secuencia.len() >= largo_de_contig_largo;
        let mut ubicacion: Option<(usize, usize)> = None;
        recorrer_kmeros(&contig.secuencia, k, |_, kmero| {
            let posiciones = indice.posiciones(kmero.canonico());
            comparacion.kmeros_de_contigs += 1;
            if posiciones.is_empty() {
                comparacion.kmeros_ajenos += 1;
            }
            for &(_, posicion) in posiciones {
                let posicion = posicion as usize;
                cubierto[posicion..posicion + k].fill(true);
                if es_largo {
                    cubierto_por_largos[posicion..posicion + k].fill(true);
                }
            }
            if let [(_, unica)] = posiciones {
                let posicion = *unica as usize;
                ubicacion = Some(match ubicacion {
                    Some((primera, ultima)) => (primera.min(posicion), ultima.max(posicion + k)),
                    None => (posicion, posicion + k),
                });
            }
        });
        comparacion.ubicaciones.push(ubicacion);
    }
    for &(inicio, fin) in &indice.regiones {
        comparacion.largo_referencia += fin - inicio;
        let mut hueco_desde: Option<usize> = None;
        for posicion in inicio..fin {
            if cubierto[posicion] {
                comparacion.cubiertas += 1;
                if let Some(desde) = hueco_desde.take() {
                    comparacion.huecos.push((desde, posicion));
                }
            } else if hueco_desde.is_none() {
                hueco_desde = Some(posicion);
            }
            if cubierto_por_largos[posicion] {
                comparacion.cubiertas_por_largos += 1;
            }
        }
        if let Some(desde) = hueco_desde {
            comparacion.huecos.push((desde, fin));
        }
    }
    comparacion
}
