use crate::formatos::{LectorFastq, Lectura, abrir_fastq};
use flate2::Compression;
use flate2::write::GzEncoder;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub const CALIDAD_MAXIMA: usize = 41;
pub const MAXIMO_POSICIONES: usize = 1000;
pub const LECTURAS_EN_MUESTRA: usize = 200_000;
pub const LECTURAS_PARA_ADAPTADOR: usize = 100_000;
pub const LONGITUD_SEMILLA: usize = 10;
pub const SOLAPAMIENTO_MINIMO: usize = 3;
pub const BASES_POR_ERROR: usize = 8;
const MAXIMO_ADAPTADOR: usize = 50;
const VECES_MAS_AL_FINAL: u32 = 2;
const DESVIOS_SOBRE_EL_AZAR: f64 = 5.0;
const ACUERDO_PARA_EXTENDER: f64 = 0.7;
const MINIMO_DE_VOTOS: usize = 10;
const MAXIMO_SOBRERREPRESENTADAS: usize = 20;
const LECTURAS_POR_BLOQUE: usize = 8192;
const BASES_POR_BLOQUE: usize = 4_000_000;
const INDICE_N: usize = 4;

#[derive(Debug, Clone, PartialEq)]
pub struct CalidadEnPosicion {
    pub lecturas: usize,
    pub media: f64,
    pub cuartil_inferior: u8,
    pub mediana: u8,
    pub cuartil_superior: u8,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Sobrerrepresentada {
    pub secuencia: Vec<u8>,
    pub veces: usize,
    pub fraccion: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AdaptadorDetectado {
    pub secuencia: Vec<u8>,
    pub lecturas: usize,
    pub revisadas: usize,
    pub fraccion: f64,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct InformeCalidad {
    pub lecturas: usize,
    pub bases: usize,
    pub longitud_minima: usize,
    pub longitud_maxima: usize,
    pub lecturas_por_longitud: Vec<usize>,
    pub calidad_por_posicion: Vec<CalidadEnPosicion>,
    pub composicion_por_posicion: Vec<[usize; 5]>,
    pub lecturas_por_gc: Vec<usize>,
    pub bases_q20: usize,
    pub bases_q30: usize,
    pub bases_n: usize,
    pub bases_gc: usize,
    pub lecturas_muestreadas: usize,
    pub duplicacion: f64,
    pub sobrerrepresentadas: Vec<Sobrerrepresentada>,
    pub adaptador: Option<AdaptadorDetectado>,
}

impl InformeCalidad {
    pub fn fraccion_q20(&self) -> f64 {
        cociente(self.bases_q20, self.bases)
    }

    pub fn fraccion_q30(&self) -> f64 {
        cociente(self.bases_q30, self.bases)
    }

    pub fn fraccion_n(&self) -> f64 {
        cociente(self.bases_n, self.bases)
    }

    pub fn gc_medio(&self) -> f64 {
        cociente(self.bases_gc, self.bases.saturating_sub(self.bases_n))
    }
}

fn cociente(numerador: usize, denominador: usize) -> f64 {
    if denominador == 0 {
        0.0
    } else {
        numerador as f64 / denominador as f64
    }
}

pub fn phred(simbolo: u8) -> usize {
    simbolo.saturating_sub(33) as usize
}

fn indice_de_base(base: u8) -> usize {
    const INDICES: [usize; 256] = {
        let mut tabla = [INDICE_N; 256];
        tabla[b'A' as usize] = 0;
        tabla[b'C' as usize] = 1;
        tabla[b'G' as usize] = 2;
        tabla[b'T' as usize] = 3;
        tabla[b'a' as usize] = 0;
        tabla[b'c' as usize] = 1;
        tabla[b'g' as usize] = 2;
        tabla[b't' as usize] = 3;
        tabla
    };
    INDICES[base as usize]
}

pub fn analizar_fastq(ruta: &Path) -> Result<InformeCalidad, String> {
    let mut lector = abrir_fastq(ruta)?;
    let mut conteos = Conteos::nuevo();
    let mut muestra: Vec<Vec<u8>> = Vec::new();
    let mut tanda = leer_tanda(&mut lector)?;
    while !tanda.is_empty() {
        let (siguiente, parciales) =
            rayon::join(|| leer_tanda(&mut lector), || contar_tanda(&tanda));
        for parcial in &parciales {
            conteos.sumar(parcial);
        }
        guardar_en_muestra(&mut muestra, tanda);
        tanda = siguiente?;
    }
    Ok(armar_informe(conteos, &muestra))
}

fn leer_tanda(lector: &mut LectorFastq) -> Result<Vec<Vec<Lectura>>, String> {
    let bloques_por_tanda = rayon::current_num_threads().max(1) * 2;
    let mut tanda = Vec::with_capacity(bloques_por_tanda);
    while tanda.len() < bloques_por_tanda {
        let bloque = leer_bloque(lector)?;
        if bloque.is_empty() {
            break;
        }
        tanda.push(bloque);
    }
    Ok(tanda)
}

fn leer_bloque(lector: &mut LectorFastq) -> Result<Vec<Lectura>, String> {
    let mut bloque = Vec::with_capacity(LECTURAS_POR_BLOQUE);
    let mut bases = 0;
    while bloque.len() < LECTURAS_POR_BLOQUE && bases < BASES_POR_BLOQUE {
        match lector.next() {
            Some(lectura) => {
                let lectura = lectura?;
                bases += lectura.secuencia.len();
                bloque.push(lectura);
            }
            None => break,
        }
    }
    Ok(bloque)
}

fn contar_tanda(tanda: &[Vec<Lectura>]) -> Vec<Conteos> {
    tanda
        .par_iter()
        .map(|bloque| {
            let mut conteos = Conteos::nuevo();
            for lectura in bloque {
                conteos.contar_lectura(lectura);
            }
            conteos
        })
        .collect()
}

fn guardar_en_muestra(muestra: &mut Vec<Vec<u8>>, tanda: Vec<Vec<Lectura>>) {
    for bloque in tanda {
        for lectura in bloque {
            if muestra.len() >= LECTURAS_EN_MUESTRA {
                return;
            }
            muestra.push(lectura.secuencia);
        }
    }
}

struct Conteos {
    lecturas: usize,
    bases: usize,
    lecturas_por_longitud: Vec<usize>,
    calidades_por_posicion: Vec<[usize; CALIDAD_MAXIMA + 1]>,
    suma_de_calidad_por_posicion: Vec<usize>,
    composicion_por_posicion: Vec<[usize; 5]>,
    lecturas_por_gc: Vec<usize>,
    bases_q20: usize,
    bases_q30: usize,
    bases_n: usize,
    bases_gc: usize,
}

impl Conteos {
    fn nuevo() -> Conteos {
        Conteos {
            lecturas: 0,
            bases: 0,
            lecturas_por_longitud: Vec::new(),
            calidades_por_posicion: Vec::new(),
            suma_de_calidad_por_posicion: Vec::new(),
            composicion_por_posicion: Vec::new(),
            lecturas_por_gc: vec![0; 101],
            bases_q20: 0,
            bases_q30: 0,
            bases_n: 0,
            bases_gc: 0,
        }
    }

    fn hacer_lugar(&mut self, posiciones: usize) {
        if self.calidades_por_posicion.len() < posiciones {
            self.calidades_por_posicion
                .resize(posiciones, [0; CALIDAD_MAXIMA + 1]);
            self.suma_de_calidad_por_posicion.resize(posiciones, 0);
            self.composicion_por_posicion.resize(posiciones, [0; 5]);
        }
    }

    fn contar_lectura(&mut self, lectura: &Lectura) {
        let longitud = lectura.secuencia.len();
        self.lecturas += 1;
        self.bases += longitud;
        if self.lecturas_por_longitud.len() <= longitud {
            self.lecturas_por_longitud.resize(longitud + 1, 0);
        }
        self.lecturas_por_longitud[longitud] += 1;
        self.hacer_lugar(longitud.min(MAXIMO_POSICIONES));
        let mut letras = [0; 5];
        for (posicion, (&base, &simbolo)) in
            lectura.secuencia.iter().zip(&lectura.calidad).enumerate()
        {
            let calidad = phred(simbolo);
            let letra = indice_de_base(base);
            letras[letra] += 1;
            self.bases_q20 += usize::from(calidad >= 20);
            self.bases_q30 += usize::from(calidad >= 30);
            if posicion < MAXIMO_POSICIONES {
                self.calidades_por_posicion[posicion][calidad.min(CALIDAD_MAXIMA)] += 1;
                self.suma_de_calidad_por_posicion[posicion] += calidad;
                self.composicion_por_posicion[posicion][letra] += 1;
            }
        }
        let gc = letras[1] + letras[2];
        let acgt = letras[0] + gc + letras[3];
        self.bases_n += letras[INDICE_N];
        self.bases_gc += gc;
        if let Some(porcentaje_gc) = (100 * gc + acgt / 2).checked_div(acgt) {
            self.lecturas_por_gc[porcentaje_gc] += 1;
        }
    }

    fn sumar(&mut self, otro: &Conteos) {
        self.lecturas += otro.lecturas;
        self.bases += otro.bases;
        if self.lecturas_por_longitud.len() < otro.lecturas_por_longitud.len() {
            self.lecturas_por_longitud
                .resize(otro.lecturas_por_longitud.len(), 0);
        }
        for (mio, suyo) in self
            .lecturas_por_longitud
            .iter_mut()
            .zip(&otro.lecturas_por_longitud)
        {
            *mio += suyo;
        }
        self.hacer_lugar(otro.calidades_por_posicion.len());
        for posicion in 0..otro.calidades_por_posicion.len() {
            for calidad in 0..=CALIDAD_MAXIMA {
                self.calidades_por_posicion[posicion][calidad] +=
                    otro.calidades_por_posicion[posicion][calidad];
            }
            self.suma_de_calidad_por_posicion[posicion] +=
                otro.suma_de_calidad_por_posicion[posicion];
            for letra in 0..5 {
                self.composicion_por_posicion[posicion][letra] +=
                    otro.composicion_por_posicion[posicion][letra];
            }
        }
        for (mio, suyo) in self.lecturas_por_gc.iter_mut().zip(&otro.lecturas_por_gc) {
            *mio += suyo;
        }
        self.bases_q20 += otro.bases_q20;
        self.bases_q30 += otro.bases_q30;
        self.bases_n += otro.bases_n;
        self.bases_gc += otro.bases_gc;
    }
}

fn armar_informe(conteos: Conteos, muestra: &[Vec<u8>]) -> InformeCalidad {
    let (adaptador, (duplicacion, sobrerrepresentadas)) =
        rayon::join(|| detectar_adaptador(muestra), || repeticiones(muestra));
    let calidad_por_posicion = conteos
        .calidades_por_posicion
        .iter()
        .zip(&conteos.suma_de_calidad_por_posicion)
        .map(|(histograma, &suma)| resumir_posicion(histograma, suma))
        .collect();
    let longitudes = &conteos.lecturas_por_longitud;
    InformeCalidad {
        lecturas: conteos.lecturas,
        bases: conteos.bases,
        longitud_minima: longitudes.iter().position(|&veces| veces > 0).unwrap_or(0),
        longitud_maxima: longitudes.iter().rposition(|&veces| veces > 0).unwrap_or(0),
        lecturas_por_longitud: conteos.lecturas_por_longitud,
        calidad_por_posicion,
        composicion_por_posicion: conteos.composicion_por_posicion,
        lecturas_por_gc: conteos.lecturas_por_gc,
        bases_q20: conteos.bases_q20,
        bases_q30: conteos.bases_q30,
        bases_n: conteos.bases_n,
        bases_gc: conteos.bases_gc,
        lecturas_muestreadas: muestra.len(),
        duplicacion,
        sobrerrepresentadas,
        adaptador,
    }
}

fn resumir_posicion(histograma: &[usize], suma_de_calidad: usize) -> CalidadEnPosicion {
    let lecturas: usize = histograma.iter().sum();
    CalidadEnPosicion {
        lecturas,
        media: cociente(suma_de_calidad, lecturas),
        cuartil_inferior: percentil_de_histograma(histograma, 0.25),
        mediana: percentil_de_histograma(histograma, 0.5),
        cuartil_superior: percentil_de_histograma(histograma, 0.75),
    }
}

pub fn percentil_de_histograma(histograma: &[usize], fraccion: f64) -> u8 {
    let total: usize = histograma.iter().sum();
    if total == 0 {
        return 0;
    }
    let puesto = ((fraccion * total as f64).ceil() as usize).max(1);
    let mut acumuladas = 0;
    for (valor, &veces) in histograma.iter().enumerate() {
        acumuladas += veces;
        if acumuladas >= puesto {
            return valor as u8;
        }
    }
    (histograma.len() - 1) as u8
}

fn repeticiones(muestra: &[Vec<u8>]) -> (f64, Vec<Sobrerrepresentada>) {
    let veces = contar_secuencias(muestra);
    (
        tasa_de_duplicacion(&veces, muestra.len()),
        secuencias_sobrerrepresentadas(&veces, muestra.len()),
    )
}

pub fn contar_secuencias(muestra: &[Vec<u8>]) -> HashMap<&[u8], usize> {
    let mut veces: HashMap<&[u8], usize> = HashMap::with_capacity(muestra.len());
    for secuencia in muestra {
        *veces.entry(secuencia.as_slice()).or_insert(0) += 1;
    }
    veces
}

pub fn tasa_de_duplicacion(veces: &HashMap<&[u8], usize>, total: usize) -> f64 {
    cociente(total.saturating_sub(veces.len()), total)
}

pub fn secuencias_sobrerrepresentadas(
    veces: &HashMap<&[u8], usize>,
    total: usize,
) -> Vec<Sobrerrepresentada> {
    let mut frecuentes: Vec<Sobrerrepresentada> = Vec::new();
    for (&secuencia, &cantidad) in veces {
        if secuencia.is_empty() || cantidad < 2 || cantidad * 1000 <= total {
            continue;
        }
        frecuentes.push(Sobrerrepresentada {
            secuencia: secuencia.to_vec(),
            veces: cantidad,
            fraccion: cociente(cantidad, total),
        });
    }
    frecuentes.sort_by(|una, otra| {
        otra.veces
            .cmp(&una.veces)
            .then_with(|| una.secuencia.cmp(&otra.secuencia))
    });
    frecuentes.truncate(MAXIMO_SOBRERREPRESENTADAS);
    frecuentes
}

pub fn detectar_adaptador(muestra: &[Vec<u8>]) -> Option<AdaptadorDetectado> {
    let muestra = &muestra[..muestra.len().min(LECTURAS_PARA_ADAPTADOR)];
    if muestra.is_empty() {
        return None;
    }
    let semilla = elegir_semilla(muestra)?;
    let secuencia = extender_semilla(muestra, &semilla);
    let lecturas = muestra
        .par_iter()
        .filter(|lectura| buscar_adaptador(lectura, &secuencia, LONGITUD_SEMILLA).is_some())
        .count();
    if lecturas == 0 {
        return None;
    }
    Some(AdaptadorDetectado {
        secuencia,
        lecturas,
        revisadas: muestra.len(),
        fraccion: cociente(lecturas, muestra.len()),
    })
}

struct ConteoDeSemillas {
    al_inicio: Vec<u32>,
    al_final: Vec<u32>,
}

impl ConteoDeSemillas {
    fn en_cero() -> ConteoDeSemillas {
        let posibles = 1 << (2 * LONGITUD_SEMILLA);
        ConteoDeSemillas {
            al_inicio: vec![0; posibles],
            al_final: vec![0; posibles],
        }
    }

    fn anotar_semillas(&mut self, secuencia: &[u8]) {
        let mitad = secuencia.len() / 2;
        let mascara = (1 << (2 * LONGITUD_SEMILLA)) - 1;
        let mut codigo = 0;
        let mut seguidas = 0;
        for (posicion, &base) in secuencia.iter().enumerate() {
            let valor = indice_de_base(base);
            if valor == INDICE_N {
                seguidas = 0;
                continue;
            }
            codigo = ((codigo << 2) | valor) & mascara;
            seguidas += 1;
            if seguidas < LONGITUD_SEMILLA {
                continue;
            }
            let inicio = posicion + 1 - LONGITUD_SEMILLA;
            if inicio < mitad {
                self.al_inicio[codigo] += 1;
            } else {
                self.al_final[codigo] += 1;
            }
        }
    }

    fn unir(mut self, otro: ConteoDeSemillas) -> ConteoDeSemillas {
        for (mio, suyo) in self.al_inicio.iter_mut().zip(&otro.al_inicio) {
            *mio += suyo;
        }
        for (mio, suyo) in self.al_final.iter_mut().zip(&otro.al_final) {
            *mio += suyo;
        }
        self
    }
}

fn contar_semillas(muestra: &[Vec<u8>]) -> ConteoDeSemillas {
    let hilos = rayon::current_num_threads().max(1);
    let lecturas_por_trozo = muestra.len().div_ceil(hilos).max(1);
    muestra
        .par_chunks(lecturas_por_trozo)
        .map(|trozo| {
            let mut conteo = ConteoDeSemillas::en_cero();
            for secuencia in trozo {
                conteo.anotar_semillas(secuencia);
            }
            conteo
        })
        .reduce_with(ConteoDeSemillas::unir)
        .unwrap_or_else(ConteoDeSemillas::en_cero)
}

fn elegir_semilla(muestra: &[Vec<u8>]) -> Option<Vec<u8>> {
    let conteo = contar_semillas(muestra);
    let minimo = (muestra.len() / 1000).max(20) as u32;
    let mut mejor: Option<usize> = None;
    for codigo in 0..conteo.al_final.len() {
        let al_final = conteo.al_final[codigo];
        if al_final < minimo || !enriquecida_al_final(al_final, conteo.al_inicio[codigo]) {
            continue;
        }
        if mejor.is_some_and(|elegido| conteo.al_final[elegido] >= al_final) {
            continue;
        }
        if es_poco_variada(&semilla_desde_codigo(codigo)) {
            continue;
        }
        mejor = Some(codigo);
    }
    mejor.map(semilla_desde_codigo)
}

fn enriquecida_al_final(al_final: u32, al_inicio: u32) -> bool {
    let diferencia = al_final as f64 - al_inicio as f64;
    let ruido = ((al_final + al_inicio) as f64).sqrt();
    al_final >= VECES_MAS_AL_FINAL * al_inicio && diferencia >= DESVIOS_SOBRE_EL_AZAR * ruido
}

fn semilla_desde_codigo(codigo: usize) -> Vec<u8> {
    let mut semilla = vec![b'A'; LONGITUD_SEMILLA];
    let mut resto = codigo;
    for posicion in (0..LONGITUD_SEMILLA).rev() {
        semilla[posicion] = b"ACGT"[resto & 3];
        resto >>= 2;
    }
    semilla
}

fn es_poco_variada(semilla: &[u8]) -> bool {
    let mut vistas = [false; 5];
    for &base in semilla {
        vistas[indice_de_base(base)] = true;
    }
    vistas[..4].iter().filter(|&&vista| vista).count() < 3
}

fn extender_semilla(muestra: &[Vec<u8>], semilla: &[u8]) -> Vec<u8> {
    let apariciones: Vec<(usize, usize)> = muestra
        .par_iter()
        .enumerate()
        .filter_map(|(indice, lectura)| {
            buscar_exacto(lectura, semilla).map(|inicio| (indice, inicio))
        })
        .collect();
    let minimo_de_votos = (apariciones.len() / 5).max(MINIMO_DE_VOTOS);
    let mut hacia_atras: Vec<u8> = Vec::new();
    while semilla.len() + hacia_atras.len() < MAXIMO_ADAPTADOR {
        let distancia = hacia_atras.len() + 1;
        let mut votos = [0; 5];
        for &(indice, inicio) in &apariciones {
            if inicio >= distancia {
                votar(&mut votos, muestra[indice][inicio - distancia]);
            }
        }
        match consenso(&votos, minimo_de_votos) {
            Some(base) => hacia_atras.push(base),
            None => break,
        }
    }
    let mut adaptador: Vec<u8> = hacia_atras.into_iter().rev().collect();
    adaptador.extend_from_slice(semilla);
    let mut agregadas = 0;
    while adaptador.len() < MAXIMO_ADAPTADOR {
        let mut votos = [0; 5];
        for &(indice, inicio) in &apariciones {
            let despues = inicio + semilla.len() + agregadas;
            if let Some(&base) = muestra[indice].get(despues) {
                votar(&mut votos, base);
            }
        }
        match consenso(&votos, minimo_de_votos) {
            Some(base) => {
                adaptador.push(base);
                agregadas += 1;
            }
            None => break,
        }
    }
    adaptador
}

fn buscar_exacto(lectura: &[u8], semilla: &[u8]) -> Option<usize> {
    lectura
        .windows(semilla.len())
        .position(|ventana| ventana == semilla)
}

fn votar(votos: &mut [usize], base: u8) {
    votos[indice_de_base(base)] += 1;
}

fn consenso(votos: &[usize], minimo_de_votos: usize) -> Option<u8> {
    let acgt = &votos[..4];
    let total: usize = acgt.iter().sum();
    if total < minimo_de_votos {
        return None;
    }
    let (mejor, &cantidad) = acgt.iter().enumerate().max_by_key(|&(_, &veces)| veces)?;
    if cantidad as f64 >= ACUERDO_PARA_EXTENDER * total as f64 {
        Some(b"ACGT"[mejor])
    } else {
        None
    }
}

pub fn buscar_adaptador(
    secuencia: &[u8],
    adaptador: &[u8],
    solapamiento_minimo: usize,
) -> Option<usize> {
    if adaptador.is_empty() {
        return None;
    }
    let minimo = solapamiento_minimo.max(1);
    for inicio in 0..secuencia.len() {
        let solapamiento = adaptador.len().min(secuencia.len() - inicio);
        if solapamiento < minimo {
            return None;
        }
        let permitidos = solapamiento / BASES_POR_ERROR;
        let tramo = &secuencia[inicio..inicio + solapamiento];
        if coincide_con_errores(tramo, &adaptador[..solapamiento], permitidos) {
            return Some(inicio);
        }
    }
    None
}

fn coincide_con_errores(tramo: &[u8], adaptador: &[u8], permitidos: usize) -> bool {
    let errores = tramo
        .iter()
        .zip(adaptador)
        .filter(|(base, esperada)| base != esperada)
        .count();
    errores <= permitidos
}

pub fn corte_por_calidad(calidad: &[u8], ventana: usize, umbral: usize) -> usize {
    let ventana = ventana.max(1);
    let mut fin = calidad.len();
    let mut inicio = fin.saturating_sub(ventana);
    let mut suma: usize = calidad[inicio..fin]
        .iter()
        .map(|&simbolo| phred(simbolo))
        .sum();
    while fin > 0 && suma < umbral.saturating_mul(fin - inicio) {
        fin -= 1;
        suma -= phred(calidad[fin]);
        if inicio > 0 {
            inicio -= 1;
            suma += phred(calidad[inicio]);
        }
    }
    fin
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ModoAdaptador {
    Detectar,
    Ninguno,
    Conocido(Vec<u8>),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OpcionesRecorte {
    pub adaptador: ModoAdaptador,
    pub ventana: usize,
    pub calidad_minima: usize,
    pub maximo_n: usize,
    pub longitud_minima: usize,
}

impl Default for OpcionesRecorte {
    fn default() -> OpcionesRecorte {
        OpcionesRecorte {
            adaptador: ModoAdaptador::Detectar,
            ventana: 4,
            calidad_minima: 20,
            maximo_n: 5,
            longitud_minima: 50,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Descarte {
    DemasiadasN,
    MuyCorta,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Recorte {
    pub inicio_adaptador: Option<usize>,
    pub fin: usize,
    pub descarte: Option<Descarte>,
}

pub fn recortar_lectura(
    lectura: &Lectura,
    adaptador: &[u8],
    opciones: &OpcionesRecorte,
) -> Recorte {
    let inicio_adaptador = buscar_adaptador(&lectura.secuencia, adaptador, SOLAPAMIENTO_MINIMO);
    let sin_adaptador = inicio_adaptador
        .unwrap_or(lectura.secuencia.len())
        .min(lectura.calidad.len());
    let fin = corte_por_calidad(
        &lectura.calidad[..sin_adaptador],
        opciones.ventana,
        opciones.calidad_minima,
    );
    let descarte = if contar_n(&lectura.secuencia[..fin]) > opciones.maximo_n {
        Some(Descarte::DemasiadasN)
    } else if fin < opciones.longitud_minima {
        Some(Descarte::MuyCorta)
    } else {
        None
    };
    Recorte {
        inicio_adaptador,
        fin,
        descarte,
    }
}

fn contar_n(secuencia: &[u8]) -> usize {
    secuencia
        .iter()
        .filter(|&&base| indice_de_base(base) == INDICE_N)
        .count()
}

fn contar_q30(calidad: &[u8]) -> usize {
    calidad
        .iter()
        .filter(|&&simbolo| phred(simbolo) >= 30)
        .count()
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct ResumenRecorte {
    pub adaptador: Vec<u8>,
    pub lecturas_entrada: usize,
    pub lecturas_salida: usize,
    pub bases_entrada: usize,
    pub bases_salida: usize,
    pub lecturas_con_adaptador: usize,
    pub lecturas_con_adaptador_corto: usize,
    pub bases_de_adaptador: usize,
    pub lecturas_recortadas_por_calidad: usize,
    pub bases_por_calidad: usize,
    pub descartadas_por_n: usize,
    pub descartadas_por_longitud: usize,
    pub bases_q30_entrada: usize,
    pub bases_q30_salida: usize,
}

impl ResumenRecorte {
    fn anotar(&mut self, lectura: &Lectura, recorte: &Recorte) {
        let longitud = lectura.secuencia.len();
        let sin_adaptador = recorte.inicio_adaptador.unwrap_or(longitud);
        self.lecturas_entrada += 1;
        self.bases_entrada += longitud;
        self.bases_q30_entrada += contar_q30(&lectura.calidad);
        if recorte.inicio_adaptador.is_some() {
            let cortadas = longitud - sin_adaptador;
            self.lecturas_con_adaptador += 1;
            self.bases_de_adaptador += cortadas;
            if cortadas < BASES_POR_ERROR {
                self.lecturas_con_adaptador_corto += 1;
            }
        }
        if recorte.fin < sin_adaptador {
            self.lecturas_recortadas_por_calidad += 1;
            self.bases_por_calidad += sin_adaptador - recorte.fin;
        }
        match recorte.descarte {
            Some(Descarte::DemasiadasN) => self.descartadas_por_n += 1,
            Some(Descarte::MuyCorta) => self.descartadas_por_longitud += 1,
            None => {
                self.lecturas_salida += 1;
                self.bases_salida += recorte.fin;
                self.bases_q30_salida += contar_q30(&lectura.calidad[..recorte.fin]);
            }
        }
    }

    fn acumular(&mut self, otro: &ResumenRecorte) {
        self.lecturas_entrada += otro.lecturas_entrada;
        self.lecturas_salida += otro.lecturas_salida;
        self.bases_entrada += otro.bases_entrada;
        self.bases_salida += otro.bases_salida;
        self.lecturas_con_adaptador += otro.lecturas_con_adaptador;
        self.lecturas_con_adaptador_corto += otro.lecturas_con_adaptador_corto;
        self.bases_de_adaptador += otro.bases_de_adaptador;
        self.lecturas_recortadas_por_calidad += otro.lecturas_recortadas_por_calidad;
        self.bases_por_calidad += otro.bases_por_calidad;
        self.descartadas_por_n += otro.descartadas_por_n;
        self.descartadas_por_longitud += otro.descartadas_por_longitud;
        self.bases_q30_entrada += otro.bases_q30_entrada;
        self.bases_q30_salida += otro.bases_q30_salida;
    }
}

pub fn recortar_fastq(
    entrada: &Path,
    salida: &Path,
    opciones: &OpcionesRecorte,
) -> Result<ResumenRecorte, String> {
    revisar_opciones(entrada, salida, opciones)?;
    let adaptador = elegir_adaptador(entrada, &opciones.adaptador)?;
    let mut lector = abrir_fastq(entrada)?;
    let comprimir = salida
        .extension()
        .is_some_and(|extension| extension == "gz");
    let archivo = File::create(salida)
        .map_err(|error| format!("No se pudo crear {}: {error}", salida.display()))?;
    let mut destino = BufWriter::with_capacity(1 << 20, archivo);
    let mut resumen = ResumenRecorte {
        adaptador: adaptador.clone(),
        ..ResumenRecorte::default()
    };
    let mut escritos = 0;
    let mut tanda = leer_tanda(&mut lector)?;
    while !tanda.is_empty() {
        let (siguiente, recortados) = rayon::join(
            || leer_tanda(&mut lector),
            || recortar_tanda(&tanda, &adaptador, opciones, comprimir),
        );
        for recortado in &recortados? {
            escritos += recortado.texto.len();
            destino
                .write_all(&recortado.texto)
                .map_err(|error| format!("No se pudo escribir {}: {error}", salida.display()))?;
            resumen.acumular(&recortado.resumen);
        }
        tanda = siguiente?;
    }
    if comprimir && escritos == 0 {
        destino
            .write_all(&comprimir_con_gzip(&[])?)
            .map_err(|error| format!("No se pudo escribir {}: {error}", salida.display()))?;
    }
    destino
        .flush()
        .map_err(|error| format!("No se pudo escribir {}: {error}", salida.display()))?;
    Ok(resumen)
}

fn revisar_opciones(
    entrada: &Path,
    salida: &Path,
    opciones: &OpcionesRecorte,
) -> Result<(), String> {
    if opciones.ventana == 0 {
        return Err("La ventana de calidad debe tener al menos 1 base.".to_string());
    }
    if let (Ok(original), Ok(destino)) = (entrada.canonicalize(), salida.canonicalize())
        && original == destino
    {
        return Err("El archivo de salida no puede ser el mismo que el de entrada.".to_string());
    }
    Ok(())
}

fn elegir_adaptador(entrada: &Path, modo: &ModoAdaptador) -> Result<Vec<u8>, String> {
    match modo {
        ModoAdaptador::Ninguno => Ok(Vec::new()),
        ModoAdaptador::Conocido(secuencia) => Ok(secuencia.clone()),
        ModoAdaptador::Detectar => {
            let muestra = leer_muestra(entrada)?;
            Ok(detectar_adaptador(&muestra)
                .map(|detectado| detectado.secuencia)
                .unwrap_or_default())
        }
    }
}

fn leer_muestra(ruta: &Path) -> Result<Vec<Vec<u8>>, String> {
    let mut muestra = Vec::new();
    for lectura in abrir_fastq(ruta)?.take(LECTURAS_PARA_ADAPTADOR) {
        muestra.push(lectura?.secuencia);
    }
    Ok(muestra)
}

struct BloqueRecortado {
    texto: Vec<u8>,
    resumen: ResumenRecorte,
}

fn recortar_tanda(
    tanda: &[Vec<Lectura>],
    adaptador: &[u8],
    opciones: &OpcionesRecorte,
    comprimir: bool,
) -> Result<Vec<BloqueRecortado>, String> {
    tanda
        .par_iter()
        .map(|bloque| recortar_bloque(bloque, adaptador, opciones, comprimir))
        .collect()
}

fn recortar_bloque(
    bloque: &[Lectura],
    adaptador: &[u8],
    opciones: &OpcionesRecorte,
    comprimir: bool,
) -> Result<BloqueRecortado, String> {
    let mut texto = Vec::new();
    let mut resumen = ResumenRecorte::default();
    for lectura in bloque {
        let recorte = recortar_lectura(lectura, adaptador, opciones);
        resumen.anotar(lectura, &recorte);
        if recorte.descarte.is_none() {
            escribir_recortada(&mut texto, lectura, recorte.fin);
        }
    }
    if comprimir && !texto.is_empty() {
        texto = comprimir_con_gzip(&texto)?;
    }
    Ok(BloqueRecortado { texto, resumen })
}

fn escribir_recortada(texto: &mut Vec<u8>, lectura: &Lectura, fin: usize) {
    texto.push(b'@');
    texto.extend_from_slice(&lectura.nombre);
    texto.push(b'\n');
    texto.extend_from_slice(&lectura.secuencia[..fin]);
    texto.extend_from_slice(b"\n+\n");
    texto.extend_from_slice(&lectura.calidad[..fin]);
    texto.push(b'\n');
}

fn comprimir_con_gzip(texto: &[u8]) -> Result<Vec<u8>, String> {
    let mut compresor = GzEncoder::new(Vec::with_capacity(texto.len() / 3), Compression::fast());
    compresor
        .write_all(texto)
        .map_err(|error| format!("No se pudo comprimir la salida: {error}"))?;
    compresor
        .finish()
        .map_err(|error| format!("No se pudo comprimir la salida: {error}"))
}
