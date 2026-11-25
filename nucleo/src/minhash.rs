use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::BTreeSet;

pub const K_POR_DEFECTO: usize = 21;
pub const CUPO_POR_DEFECTO: usize = 1000;
pub const K_MAXIMO: usize = 32;
pub const CUPO_MAXIMO: usize = 100_000;
pub const GENOMAS_MAXIMOS: usize = 100;
pub const HASHES_MAXIMOS_EN_TOTAL: usize = 1_000_000;

const SEMILLA: u64 = 0x9e37_79b9_7f4a_7c15;
const NO_ES_BASE: u8 = 4;
const CODIGO_DE_BASE: [u8; 256] = tabla_de_codigos();

const fn tabla_de_codigos() -> [u8; 256] {
    let mut tabla = [NO_ES_BASE; 256];
    tabla[b'A' as usize] = 0;
    tabla[b'a' as usize] = 0;
    tabla[b'C' as usize] = 1;
    tabla[b'c' as usize] = 1;
    tabla[b'G' as usize] = 2;
    tabla[b'g' as usize] = 2;
    tabla[b'T' as usize] = 3;
    tabla[b't' as usize] = 3;
    tabla[b'U' as usize] = 3;
    tabla[b'u' as usize] = 3;
    tabla
}

pub fn mezclar(kmero: u64) -> u64 {
    let mut x = kmero ^ SEMILLA;
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51_afd7_ed55_8ccd);
    x ^= x >> 33;
    x = x.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    x ^= x >> 33;
    x
}

#[derive(Debug, Clone)]
pub struct Ventana {
    pub k: usize,
    pub adelante: u64,
    pub atras: u64,
    pub validas: usize,
    mascara: u64,
    corrimiento: u32,
}

impl Ventana {
    pub fn nueva(k: usize) -> Ventana {
        let k = k.clamp(1, K_MAXIMO);
        Ventana {
            k,
            adelante: 0,
            atras: 0,
            validas: 0,
            mascara: u64::MAX >> (64 - 2 * k),
            corrimiento: 2 * (k as u32 - 1),
        }
    }

    pub fn reiniciar(&mut self) {
        self.adelante = 0;
        self.atras = 0;
        self.validas = 0;
    }

    pub fn agregar(&mut self, base: u8) -> Option<u64> {
        let codigo = CODIGO_DE_BASE[base as usize];
        if codigo == NO_ES_BASE {
            self.validas = 0;
            return None;
        }
        let codigo = codigo as u64;
        self.adelante = ((self.adelante << 2) | codigo) & self.mascara;
        self.atras = (self.atras >> 2) | ((3 - codigo) << self.corrimiento);
        self.validas += 1;
        if self.validas < self.k {
            return None;
        }
        Some(mezclar(self.adelante.min(self.atras)))
    }
}

#[derive(Debug, Clone)]
pub struct Menores {
    pub cupo: usize,
    pub conjunto: BTreeSet<u64>,
    pub mayor_guardado: u64,
}

impl Menores {
    pub fn nuevo(cupo: usize) -> Menores {
        Menores {
            cupo,
            conjunto: BTreeSet::new(),
            mayor_guardado: u64::MAX,
        }
    }

    pub fn ofrecer(&mut self, hash: u64) -> bool {
        let lleno = self.conjunto.len() >= self.cupo;
        if lleno && hash >= self.mayor_guardado {
            return false;
        }
        let es_nuevo = self.conjunto.insert(hash);
        if es_nuevo && self.conjunto.len() > self.cupo {
            self.conjunto.pop_last();
        }
        if self.conjunto.len() == self.cupo
            && let Some(&mayor) = self.conjunto.last()
        {
            self.mayor_guardado = mayor;
        }
        es_nuevo
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Bosquejo {
    pub k: usize,
    pub cupo: usize,
    pub hashes: Vec<u64>,
    pub longitud: usize,
    pub kmeros: usize,
}

pub fn bosquejar(piezas: &[Vec<u8>], k: usize, cupo: usize) -> Bosquejo {
    let mut ventana = Ventana::nueva(k);
    let mut menores = Menores::nuevo(cupo);
    let mut longitud = 0;
    let mut kmeros = 0;
    for pieza in piezas {
        ventana.reiniciar();
        longitud += pieza.len();
        for &base in pieza {
            if let Some(hash) = ventana.agregar(base) {
                kmeros += 1;
                menores.ofrecer(hash);
            }
        }
    }
    Bosquejo {
        k,
        cupo,
        hashes: menores.conjunto.into_iter().collect(),
        longitud,
        kmeros,
    }
}

pub fn en_fraccion(hash: u64) -> f64 {
    hash as f64 / 2f64.powi(64)
}

pub fn estimar_kmeros_distintos(bosquejo: &Bosquejo) -> f64 {
    let guardados = bosquejo.hashes.len();
    match bosquejo.hashes.last() {
        Some(&mayor) if guardados >= bosquejo.cupo && guardados > 1 => {
            (guardados - 1) as f64 / en_fraccion(mayor)
        }
        _ => guardados as f64,
    }
}

pub fn contar_compartidos(a: &[u64], b: &[u64], cupo: usize) -> (usize, usize) {
    let mut i = 0;
    let mut j = 0;
    let mut compartidos = 0;
    let mut vistos = 0;
    while vistos < cupo && i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            Ordering::Equal => {
                compartidos += 1;
                i += 1;
                j += 1;
            }
            Ordering::Less => i += 1,
            Ordering::Greater => j += 1,
        }
        vistos += 1;
    }
    let sobrantes = (a.len() - i) + (b.len() - j);
    vistos += sobrantes.min(cupo - vistos);
    (compartidos, vistos)
}

pub fn fraccion(parte: usize, total: usize) -> f64 {
    if total == 0 {
        return 0.0;
    }
    parte as f64 / total as f64
}

pub fn distancia_mash(jaccard: f64, k: usize) -> f64 {
    if jaccard <= 0.0 {
        return 1.0;
    }
    let distancia = -(2.0 * jaccard / (1.0 + jaccard)).ln() / k as f64;
    distancia.clamp(0.0, 1.0)
}

pub fn probabilidad_de_kmero_al_azar(longitud: usize, k: usize) -> f64 {
    let espacio = 4f64.powi(k as i32);
    1.0 / (1.0 + espacio / longitud as f64)
}

pub fn jaccard_al_azar(longitud_a: usize, longitud_b: usize, k: usize) -> f64 {
    let pa = probabilidad_de_kmero_al_azar(longitud_a, k);
    let pb = probabilidad_de_kmero_al_azar(longitud_b, k);
    let en_alguno = pa + pb - pa * pb;
    if en_alguno <= 0.0 {
        return 0.0;
    }
    pa * pb / en_alguno
}

pub fn log10_valor_p(compartidos: usize, vistos: usize, azar: f64) -> f64 {
    if compartidos == 0 || compartidos > vistos || azar >= 1.0 {
        return 0.0;
    }
    if azar <= 0.0 {
        return f64::NEG_INFINITY;
    }
    let total = vistos as f64;
    let mut termino = ln_combinaciones(vistos, compartidos)
        + compartidos as f64 * azar.ln()
        + (total - compartidos as f64) * (-azar).ln_1p();
    let mut cola = termino;
    let razon = azar.ln() - (-azar).ln_1p();
    for exitos in compartidos..vistos {
        termino += ((vistos - exitos) as f64 / (exitos + 1) as f64).ln() + razon;
        cola = sumar_en_logaritmos(cola, termino);
        if termino < cola - 40.0 && exitos as f64 > total * azar {
            break;
        }
    }
    (cola / std::f64::consts::LN_10).min(0.0)
}

fn ln_combinaciones(total: usize, elegidos: usize) -> f64 {
    let mut suma = 0.0;
    for i in 1..=elegidos {
        suma += ((total - elegidos + i) as f64 / i as f64).ln();
    }
    suma
}

fn sumar_en_logaritmos(a: f64, b: f64) -> f64 {
    let mayor = a.max(b);
    if mayor == f64::NEG_INFINITY {
        return mayor;
    }
    mayor + ((a - mayor).exp() + (b - mayor).exp()).ln()
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Medida {
    pub compartidos: usize,
    pub vistos: usize,
    pub jaccard: f64,
    pub distancia: f64,
    pub identidad: f64,
    pub log10_valor_p: f64,
}

pub fn medir(a: &Bosquejo, b: &Bosquejo) -> Medida {
    let (compartidos, vistos) = contar_compartidos(&a.hashes, &b.hashes, a.cupo.min(b.cupo));
    let jaccard = fraccion(compartidos, vistos);
    let distancia = distancia_mash(jaccard, a.k);
    let azar = jaccard_al_azar(a.longitud, b.longitud, a.k);
    Medida {
        compartidos,
        vistos,
        jaccard,
        distancia,
        identidad: 1.0 - distancia,
        log10_valor_p: log10_valor_p(compartidos, vistos, azar),
    }
}

pub fn identidad_de_contencion(contencion: f64, k: usize) -> f64 {
    if contencion <= 0.0 {
        return 0.0;
    }
    contencion.powf(1.0 / k as f64)
}

fn hashes_de_todos(bosquejos: &[Bosquejo]) -> Vec<u64> {
    let mut todos: Vec<u64> = bosquejos
        .iter()
        .flat_map(|bosquejo| bosquejo.hashes.iter().copied())
        .collect();
    todos.sort_unstable();
    todos.dedup();
    todos
}

fn marcar_presentes(todos: &[u64], piezas: &[Vec<u8>], k: usize) -> Vec<bool> {
    let mut presentes = vec![false; todos.len()];
    let Some(&tope) = todos.last() else {
        return presentes;
    };
    let mut ventana = Ventana::nueva(k);
    for pieza in piezas {
        ventana.reiniciar();
        for &base in pieza {
            if let Some(hash) = ventana.agregar(base)
                && hash <= tope
                && let Ok(posicion) = todos.binary_search(&hash)
            {
                presentes[posicion] = true;
            }
        }
    }
    presentes
}

fn posiciones_en(todos: &[u64], bosquejo: &Bosquejo) -> Vec<usize> {
    bosquejo
        .hashes
        .iter()
        .map(|&hash| todos.partition_point(|&otro| otro < hash))
        .collect()
}

fn fraccion_presente(posiciones: &[usize], presentes: &[bool]) -> f64 {
    let hallados = posiciones
        .iter()
        .filter(|&&posicion| presentes[posicion])
        .count();
    fraccion(hallados, posiciones.len())
}

pub fn contenciones(bosquejos: &[Bosquejo], genomas: &[Vec<Vec<u8>>], k: usize) -> Vec<Vec<f64>> {
    let todos = hashes_de_todos(bosquejos);
    let presentes: Vec<Vec<bool>> = genomas
        .par_iter()
        .map(|piezas| marcar_presentes(&todos, piezas, k))
        .collect();
    let posiciones: Vec<Vec<usize>> = bosquejos
        .par_iter()
        .map(|bosquejo| posiciones_en(&todos, bosquejo))
        .collect();
    posiciones
        .par_iter()
        .map(|chico| {
            presentes
                .iter()
                .map(|en_grande| fraccion_presente(chico, en_grande))
                .collect()
        })
        .collect()
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Par {
    pub primero: usize,
    pub segundo: usize,
    pub medida: Medida,
}

fn medir_todos_los_pares(bosquejos: &[Bosquejo]) -> Vec<Par> {
    let mut indices = Vec::new();
    for primero in 0..bosquejos.len() {
        for segundo in primero + 1..bosquejos.len() {
            indices.push((primero, segundo));
        }
    }
    indices
        .par_iter()
        .map(|&(primero, segundo)| Par {
            primero,
            segundo,
            medida: medir(&bosquejos[primero], &bosquejos[segundo]),
        })
        .collect()
}

#[derive(Debug, Clone, PartialEq)]
pub struct Comparacion {
    pub bosquejos: Vec<Bosquejo>,
    pub pares: Vec<Par>,
    pub contencion: Vec<Vec<f64>>,
}

pub fn validar(cantidad: usize, k: usize, cupo: usize) -> Result<(), String> {
    if k == 0 || k > K_MAXIMO {
        return Err(format!(
            "k debe estar entre 1 y {K_MAXIMO}: cada k-mero se guarda en 64 bits, 2 por base. Lo usual es 21."
        ));
    }
    if cupo == 0 || cupo > CUPO_MAXIMO {
        return Err(
            "Cada bosquejo debe guardar entre 1 y 100 000 hashes. Lo usual es 1 000.".to_string(),
        );
    }
    if cantidad < 2 {
        return Err("Hacen falta al menos dos secuencias para compararlas.".to_string());
    }
    if cantidad > GENOMAS_MAXIMOS {
        return Err(format!(
            "Se pueden comparar hasta {GENOMAS_MAXIMOS} secuencias a la vez y hay {cantidad}."
        ));
    }
    if cantidad * cupo > HASHES_MAXIMOS_EN_TOTAL {
        return Err(format!(
            "{cantidad} secuencias con {cupo} hashes cada una son {} hashes en total y el máximo es {}: baja los hashes por bosquejo o compara menos secuencias.",
            cantidad * cupo,
            HASHES_MAXIMOS_EN_TOTAL
        ));
    }
    Ok(())
}

pub fn comparar(genomas: &[Vec<Vec<u8>>], k: usize, cupo: usize) -> Result<Comparacion, String> {
    validar(genomas.len(), k, cupo)?;
    let bosquejos: Vec<Bosquejo> = genomas
        .par_iter()
        .map(|piezas| bosquejar(piezas, k, cupo))
        .collect();
    if let Some(vacio) = bosquejos
        .iter()
        .position(|bosquejo| bosquejo.hashes.is_empty())
    {
        return Err(format!(
            "La secuencia número {} no tiene ningún tramo de {k} letras A, C, G o T seguidas.",
            vacio + 1
        ));
    }
    let pares = medir_todos_los_pares(&bosquejos);
    let contencion = contenciones(&bosquejos, genomas, k);
    Ok(Comparacion {
        bosquejos,
        pares,
        contencion,
    })
}

pub fn matriz_de_distancias(comparacion: &Comparacion) -> Vec<Vec<f64>> {
    let cantidad = comparacion.bosquejos.len();
    let mut matriz = vec![vec![0.0; cantidad]; cantidad];
    for par in &comparacion.pares {
        matriz[par.primero][par.segundo] = par.medida.distancia;
        matriz[par.segundo][par.primero] = par.medida.distancia;
    }
    matriz
}

pub fn todos_los_kmeros(piezas: &[Vec<u8>], k: usize) -> Vec<u64> {
    let mut ventana = Ventana::nueva(k);
    let mut hashes = Vec::new();
    for pieza in piezas {
        ventana.reiniciar();
        for &base in pieza {
            if let Some(hash) = ventana.agregar(base) {
                hashes.push(hash);
            }
        }
    }
    hashes.par_sort_unstable();
    hashes.dedup();
    hashes
}

pub fn contar_comunes(a: &[u64], b: &[u64]) -> usize {
    let mut i = 0;
    let mut j = 0;
    let mut comunes = 0;
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            Ordering::Equal => {
                comunes += 1;
                i += 1;
                j += 1;
            }
            Ordering::Less => i += 1,
            Ordering::Greater => j += 1,
        }
    }
    comunes
}

pub fn jaccard_exacto(a: &[u64], b: &[u64]) -> f64 {
    let comunes = contar_comunes(a, b);
    fraccion(comunes, a.len() + b.len() - comunes)
}

pub fn contencion_exacta(chico: &[u64], grande: &[u64]) -> f64 {
    fraccion(contar_comunes(chico, grande), chico.len())
}
