use crate::ecoli::Gen;
use crate::formatos::{self, RegistroSam, Variante};
use crate::secuencia::{aminoacido, complemento, complemento_inverso, indice_base};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::io::{self, BufRead, Write};
use std::path::Path;

const BASES: &[u8] = b"ACGT";
const SECUNDARIA: u16 = 0x100;
const NO_PASA_EL_CONTROL: u16 = 0x200;
const DUPLICADA: u16 = 0x400;
const SUPLEMENTARIA: u16 = 0x800;
const CALIDAD_SI_FALTA: u8 = 30;
const LINEAS_POR_TANDA: usize = 50_000;
const LONGITUD_MAXIMA_DE_REFERENCIA: usize = 16_000_000;
const SIN_COLUMNA: u8 = 4;
const SITIOS_DE_RESISTENCIA: &[(&str, usize, usize, &str)] = &[
    ("rpoB", 507, 533, "rifampicina"),
    ("gyrA", 83, 83, "quinolonas"),
    ("gyrA", 87, 87, "quinolonas"),
    ("parC", 80, 80, "quinolonas"),
    ("parC", 84, 84, "quinolonas"),
];

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Umbrales {
    pub calidad_mapeo_minima: u8,
    pub calidad_base_minima: u8,
    pub profundidad_minima: u32,
    pub fraccion_minima: f64,
    pub fraccion_minima_indel: f64,
}

impl Default for Umbrales {
    fn default() -> Umbrales {
        Umbrales {
            calidad_mapeo_minima: 20,
            calidad_base_minima: 20,
            profundidad_minima: 5,
            fraccion_minima: 0.8,
            fraccion_minima_indel: 0.5,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Cambio {
    pub posicion: usize,
    pub referencia: Vec<u8>,
    pub alternativa: Vec<u8>,
}

impl Cambio {
    pub fn es_snp(&self) -> bool {
        self.referencia.len() == 1 && self.alternativa.len() == 1
    }

    pub fn es_insercion(&self) -> bool {
        self.alternativa.len() > self.referencia.len()
    }

    pub fn es_delecion(&self) -> bool {
        self.referencia.len() > self.alternativa.len()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Llamada {
    pub cambio: Cambio,
    pub profundidad: u32,
    pub a_favor: u32,
    pub fraccion: f64,
    pub calidad: f64,
    pub anotaciones: Vec<Anotacion>,
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ResumenDeProfundidad {
    pub con_lecturas: usize,
    pub suficientes: usize,
    pub media: f64,
    pub sin_confirmar: usize,
}

#[derive(Debug, Clone)]
pub struct Informe {
    pub nombre_referencia: String,
    pub llamadas: Vec<Llamada>,
    pub lecturas_usadas: usize,
    pub lecturas_descartadas: usize,
    pub profundidad: ResumenDeProfundidad,
}

pub fn llamar_variantes(
    ruta_sam: &Path,
    referencia: &[u8],
    genes: &[Gen],
    umbrales: &Umbrales,
) -> Result<Informe, String> {
    revisar_umbrales(referencia, umbrales)?;
    let referencia = &referencia.to_ascii_uppercase();
    let apilamiento = apilar(formatos::abrir(ruta_sam)?, referencia, umbrales)?;
    let mut llamadas = llamar_snps(&apilamiento, referencia, umbrales);
    llamadas.extend(llamar_indels(&apilamiento, umbrales));
    llamadas.sort_by(|una, otra| una.cambio.cmp(&otra.cambio));
    let indice = IndiceDeGenes::nuevo(genes, referencia.len());
    for llamada in &mut llamadas {
        llamada.anotaciones = anotar(&llamada.cambio, referencia, &indice);
    }
    Ok(Informe {
        profundidad: resumir_profundidad(&apilamiento, referencia, umbrales),
        nombre_referencia: apilamiento.nombre_referencia,
        lecturas_usadas: apilamiento.lecturas_usadas,
        lecturas_descartadas: apilamiento.lecturas_descartadas,
        llamadas,
    })
}

fn revisar_umbrales(referencia: &[u8], umbrales: &Umbrales) -> Result<(), String> {
    if referencia.is_empty() {
        return Err("La referencia está vacía.".to_string());
    }
    if umbrales.profundidad_minima == 0 {
        return Err("La profundidad mínima debe ser al menos 1.".to_string());
    }
    for fraccion in [umbrales.fraccion_minima, umbrales.fraccion_minima_indel] {
        let valida = fraccion > 0.0 && fraccion <= 1.0;
        if !valida {
            return Err("Las fracciones mínimas deben estar entre 0 y 1.".to_string());
        }
    }
    Ok(())
}

#[derive(Debug, Clone, Default)]
pub struct Apilamiento {
    pub nombre_referencia: String,
    pub bases: Vec<[u32; 4]>,
    pub borrados: Vec<u32>,
    pub cruzan: Vec<u32>,
    pub suma_error: Vec<f32>,
    pub indels: HashMap<Cambio, u32>,
    pub lecturas_usadas: usize,
    pub lecturas_descartadas: usize,
}

impl Apilamiento {
    fn vacio(longitud: usize) -> Apilamiento {
        Apilamiento {
            nombre_referencia: String::new(),
            bases: vec![[0; 4]; longitud],
            borrados: vec![0; longitud],
            cruzan: vec![0; longitud],
            suma_error: vec![0.0; longitud],
            indels: HashMap::new(),
            lecturas_usadas: 0,
            lecturas_descartadas: 0,
        }
    }

    pub fn profundidad(&self, posicion: usize) -> u32 {
        self.bases[posicion].iter().sum::<u32>() + self.borrados[posicion]
    }
}

pub fn apilar(
    mut fuente: impl BufRead,
    referencia: &[u8],
    umbrales: &Umbrales,
) -> Result<Apilamiento, String> {
    if referencia.len() > LONGITUD_MAXIMA_DE_REFERENCIA {
        return Err(format!(
            "La referencia mide {} pb y el máximo es 16 millones de pb (el tamaño de un genoma bacteriano grande).",
            referencia.len()
        ));
    }
    let errores = tabla_de_errores();
    let mut apilamiento = Apilamiento::vacio(referencia.len());
    let mut cabecera: Vec<(String, usize)> = Vec::new();
    let mut tanda = Tanda::default();
    loop {
        let quedan_lineas = tanda.llenar(&mut fuente, &mut cabecera)?;
        for (numero, mut lectura) in tanda.leer()? {
            if apilamiento.nombre_referencia.is_empty() && lectura.mapeado() {
                apilamiento.nombre_referencia =
                    elegir_referencia(&cabecera, &lectura.referencia, referencia.len())?;
            }
            if !sirve(
                &lectura,
                &apilamiento.nombre_referencia,
                referencia.len(),
                umbrales,
            ) {
                apilamiento.lecturas_descartadas += 1;
                continue;
            }
            if lectura.calidad.is_empty() {
                lectura.calidad = vec![CALIDAD_SI_FALTA + 33; lectura.secuencia.len()];
            }
            agregar_lectura(&mut apilamiento, &lectura, referencia, umbrales, &errores)
                .map_err(|error| format!("Línea {numero} del SAM: {error}"))?;
            apilamiento.lecturas_usadas += 1;
        }
        if !quedan_lineas {
            break;
        }
    }
    if apilamiento.nombre_referencia.is_empty() {
        apilamiento.nombre_referencia = cabecera
            .iter()
            .find(|(_, largo)| *largo == referencia.len())
            .map(|(nombre, _)| nombre.clone())
            .unwrap_or_else(|| "referencia".to_string());
    }
    Ok(apilamiento)
}

#[derive(Default)]
struct Tanda {
    texto: Vec<u8>,
    lineas: Vec<(usize, usize, usize)>,
    lineas_leidas: usize,
}

impl Tanda {
    fn llenar(
        &mut self,
        fuente: &mut impl BufRead,
        cabecera: &mut Vec<(String, usize)>,
    ) -> Result<bool, String> {
        self.texto.clear();
        self.lineas.clear();
        while self.lineas.len() < LINEAS_POR_TANDA {
            let inicio = self.texto.len();
            let leidos = fuente
                .read_until(b'\n', &mut self.texto)
                .map_err(|error| format!("No se pudo leer el SAM: {error}"))?;
            if leidos == 0 {
                return Ok(false);
            }
            self.lineas_leidas += 1;
            let linea = &self.texto[inicio..];
            if linea.starts_with(b"@SQ") {
                cabecera.extend(leer_secuencia_de_cabecera(linea));
            }
            if linea.starts_with(b"@") {
                self.texto.truncate(inicio);
            } else {
                self.lineas
                    .push((self.lineas_leidas, inicio, self.texto.len()));
            }
        }
        Ok(true)
    }

    fn leer(&self) -> Result<Vec<(usize, RegistroSam)>, String> {
        let leidas: Vec<Result<Option<RegistroSam>, String>> = self
            .lineas
            .par_iter()
            .map(|&(_, inicio, fin)| formatos::leer_sam(&self.texto[inicio..fin]))
            .collect();
        let mut lecturas = Vec::with_capacity(leidas.len());
        for (&(numero, _, _), leida) in self.lineas.iter().zip(leidas) {
            match leida {
                Ok(Some(lectura)) => lecturas.push((numero, lectura)),
                Ok(None) => {}
                Err(error) => return Err(format!("Línea {numero} del SAM: {error}")),
            }
        }
        Ok(lecturas)
    }
}

fn leer_secuencia_de_cabecera(linea: &[u8]) -> Option<(String, usize)> {
    let texto = String::from_utf8_lossy(linea);
    let mut nombre = None;
    let mut largo = None;
    for campo in texto.trim_end().split('\t') {
        if let Some(valor) = campo.strip_prefix("SN:") {
            nombre = Some(valor.to_string());
        }
        if let Some(valor) = campo.strip_prefix("LN:") {
            largo = valor.parse().ok();
        }
    }
    Some((nombre?, largo?))
}

fn elegir_referencia(
    cabecera: &[(String, usize)],
    nombre_en_la_lectura: &str,
    longitud: usize,
) -> Result<String, String> {
    if cabecera.is_empty() {
        return Ok(nombre_en_la_lectura.to_string());
    }
    if let Some((nombre, _)) = cabecera.iter().find(|(_, largo)| *largo == longitud) {
        return Ok(nombre.clone());
    }
    let descritas: Vec<String> = cabecera
        .iter()
        .take(3)
        .map(|(nombre, largo)| format!("{nombre} ({largo} pb)"))
        .collect();
    Err(format!(
        "El SAM se alineó contra {}, pero la referencia elegida mide {longitud} pb. Usa la misma referencia que en el mapeo.",
        descritas.join(", ")
    ))
}

pub fn sirve(
    lectura: &RegistroSam,
    nombre_referencia: &str,
    longitud_referencia: usize,
    umbrales: &Umbrales,
) -> bool {
    let descartables = SECUNDARIA | NO_PASA_EL_CONTROL | DUPLICADA | SUPLEMENTARIA;
    lectura.mapeado()
        && lectura.posicion <= longitud_referencia
        && lectura.bandera & descartables == 0
        && lectura.calidad_mapeo >= umbrales.calidad_mapeo_minima
        && lectura.referencia == nombre_referencia
        && !lectura.secuencia.is_empty()
        && !lectura.cigar.is_empty()
}

pub fn probabilidad_de_error(calidad: u8) -> f64 {
    10f64.powf(-(calidad as f64) / 10.0)
}

fn columna_de_letra(letra: u8) -> u8 {
    const COLUMNAS: [u8; 256] = {
        let mut tabla = [SIN_COLUMNA; 256];
        tabla[b'A' as usize] = 0;
        tabla[b'C' as usize] = 1;
        tabla[b'G' as usize] = 2;
        tabla[b'T' as usize] = 3;
        tabla
    };
    COLUMNAS[letra as usize]
}

fn tabla_de_errores() -> Vec<f32> {
    (0..=255)
        .map(|calidad| probabilidad_de_error(calidad) as f32)
        .collect()
}

fn agregar_lectura(
    apilamiento: &mut Apilamiento,
    lectura: &RegistroSam,
    referencia: &[u8],
    umbrales: &Umbrales,
    errores: &[f32],
) -> Result<(), String> {
    revisar_longitudes(lectura)?;
    let inicio = lectura.posicion - 1;
    let mut en_referencia = inicio;
    let mut en_lectura = 0;
    for &(largo, operacion) in &lectura.cigar {
        let largo = largo as usize;
        match operacion {
            b'M' | b'=' | b'X' => {
                contar_alineadas(
                    apilamiento,
                    en_referencia,
                    &lectura.secuencia[en_lectura..en_lectura + largo],
                    &lectura.calidad[en_lectura..en_lectura + largo],
                    umbrales,
                    errores,
                );
                en_referencia += largo;
                en_lectura += largo;
            }
            b'I' => {
                if en_referencia > inicio {
                    let insertadas = &lectura.secuencia[en_lectura..en_lectura + largo];
                    registrar_insercion(apilamiento, referencia, en_referencia - 1, insertadas);
                }
                en_lectura += largo;
            }
            b'D' => {
                if en_referencia > inicio {
                    registrar_delecion(apilamiento, referencia, en_referencia - 1, largo);
                }
                let hasta = (en_referencia + largo).min(referencia.len());
                if en_referencia < hasta {
                    for borrado in &mut apilamiento.borrados[en_referencia..hasta] {
                        *borrado += 1;
                    }
                }
                en_referencia += largo;
            }
            b'N' => en_referencia += largo,
            b'S' => en_lectura += largo,
            _ => {}
        }
    }
    let fin = en_referencia.min(referencia.len());
    if inicio + 1 < fin {
        for cuenta in &mut apilamiento.cruzan[inicio..fin - 1] {
            *cuenta += 1;
        }
    }
    Ok(())
}

fn revisar_longitudes(lectura: &RegistroSam) -> Result<(), String> {
    let mut largo_en_lectura = 0;
    for &(largo, operacion) in &lectura.cigar {
        if matches!(operacion, b'M' | b'I' | b'S' | b'=' | b'X') {
            largo_en_lectura += largo as usize;
        }
    }
    if largo_en_lectura != lectura.secuencia.len() {
        return Err("el CIGAR no coincide con la longitud de la secuencia.".to_string());
    }
    if !lectura.calidad.is_empty() && lectura.calidad.len() != lectura.secuencia.len() {
        return Err("la calidad y la secuencia no tienen la misma longitud.".to_string());
    }
    Ok(())
}

fn contar_alineadas(
    apilamiento: &mut Apilamiento,
    desde: usize,
    bases: &[u8],
    calidades: &[u8],
    umbrales: &Umbrales,
    errores: &[f32],
) {
    let hasta = (desde + bases.len()).min(apilamiento.bases.len());
    if desde >= hasta {
        return;
    }
    let calidad_minima = umbrales.calidad_base_minima;
    let posiciones = apilamiento.bases[desde..hasta]
        .iter_mut()
        .zip(&mut apilamiento.suma_error[desde..hasta]);
    let observaciones = bases.iter().zip(calidades);
    for ((columna, error_sumado), (&base, &simbolo)) in posiciones.zip(observaciones) {
        let calidad = simbolo.saturating_sub(33);
        let indice = columna_de_letra(base);
        if calidad < calidad_minima || indice == SIN_COLUMNA {
            continue;
        }
        columna[indice as usize] += 1;
        *error_sumado += errores[calidad as usize];
    }
}

fn registrar_insercion(
    apilamiento: &mut Apilamiento,
    referencia: &[u8],
    ancla: usize,
    insertadas: &[u8],
) {
    let Some(&base_ancla) = referencia.get(ancla) else {
        return;
    };
    let mut alternativa = vec![base_ancla];
    alternativa.extend_from_slice(insertadas);
    let cambio = Cambio {
        posicion: ancla,
        referencia: vec![base_ancla],
        alternativa,
    };
    sumar_indel(apilamiento, &cambio, referencia);
}

fn registrar_delecion(
    apilamiento: &mut Apilamiento,
    referencia: &[u8],
    ancla: usize,
    largo: usize,
) {
    let fin = ancla + largo + 1;
    if fin > referencia.len() {
        return;
    }
    let cambio = Cambio {
        posicion: ancla,
        referencia: referencia[ancla..fin].to_vec(),
        alternativa: vec![referencia[ancla]],
    };
    sumar_indel(apilamiento, &cambio, referencia);
}

fn sumar_indel(apilamiento: &mut Apilamiento, cambio: &Cambio, referencia: &[u8]) {
    if cambio.referencia == cambio.alternativa {
        return;
    }
    *apilamiento
        .indels
        .entry(normalizar(cambio, referencia))
        .or_insert(0) += 1;
}

pub fn normalizar(cambio: &Cambio, referencia: &[u8]) -> Cambio {
    let fin = cambio.posicion + cambio.referencia.len();
    if cambio.referencia == cambio.alternativa
        || fin > referencia.len()
        || referencia[cambio.posicion..fin] != cambio.referencia[..]
    {
        return cambio.clone();
    }
    let mut posicion = cambio.posicion;
    let mut alelo_referencia = cambio.referencia.clone();
    let mut alelo_alternativo = cambio.alternativa.clone();
    loop {
        if alelo_referencia.is_empty() || alelo_alternativo.is_empty() {
            if posicion == 0 {
                return anclar_a_la_derecha(
                    cambio,
                    alelo_referencia,
                    alelo_alternativo,
                    referencia,
                );
            }
            posicion -= 1;
            alelo_referencia.insert(0, referencia[posicion]);
            alelo_alternativo.insert(0, referencia[posicion]);
        } else if alelo_referencia.last() == alelo_alternativo.last() {
            alelo_referencia.pop();
            alelo_alternativo.pop();
        } else {
            break;
        }
    }
    while alelo_referencia.len() > 1
        && alelo_alternativo.len() > 1
        && alelo_referencia[0] == alelo_alternativo[0]
    {
        alelo_referencia.remove(0);
        alelo_alternativo.remove(0);
        posicion += 1;
    }
    Cambio {
        posicion,
        referencia: alelo_referencia,
        alternativa: alelo_alternativo,
    }
}

fn anclar_a_la_derecha(
    original: &Cambio,
    mut alelo_referencia: Vec<u8>,
    mut alelo_alternativo: Vec<u8>,
    referencia: &[u8],
) -> Cambio {
    let Some(&siguiente) = referencia.get(alelo_referencia.len()) else {
        return original.clone();
    };
    alelo_referencia.push(siguiente);
    alelo_alternativo.push(siguiente);
    Cambio {
        posicion: 0,
        referencia: alelo_referencia,
        alternativa: alelo_alternativo,
    }
}

pub fn llamar_snps(
    apilamiento: &Apilamiento,
    referencia: &[u8],
    umbrales: &Umbrales,
) -> Vec<Llamada> {
    let mut llamadas = Vec::new();
    for (posicion, &base_referencia) in referencia.iter().enumerate().take(apilamiento.bases.len())
    {
        let profundidad = apilamiento.profundidad(posicion);
        if profundidad < umbrales.profundidad_minima {
            continue;
        }
        let Some(indice_referencia) = indice_base(base_referencia) else {
            continue;
        };
        let (indice_alternativo, a_favor) =
            alternativa_mas_frecuente(&apilamiento.bases[posicion], indice_referencia);
        if a_favor == 0 {
            continue;
        }
        let fraccion = a_favor as f64 / profundidad as f64;
        if fraccion < umbrales.fraccion_minima {
            continue;
        }
        let error_hacia_esa_base = error_medio(apilamiento, posicion, umbrales) / 3.0;
        llamadas.push(Llamada {
            cambio: Cambio {
                posicion,
                referencia: vec![BASES[indice_referencia]],
                alternativa: vec![BASES[indice_alternativo]],
            },
            profundidad,
            a_favor,
            fraccion,
            calidad: calidad_phred(a_favor, profundidad, error_hacia_esa_base),
            anotaciones: Vec::new(),
        });
    }
    llamadas
}

pub fn alternativa_mas_frecuente(conteos: &[u32], indice_referencia: usize) -> (usize, u32) {
    let mut mejor = (indice_referencia, 0);
    for (indice, &conteo) in conteos.iter().enumerate() {
        if indice != indice_referencia && conteo > mejor.1 {
            mejor = (indice, conteo);
        }
    }
    mejor
}

pub fn error_medio(apilamiento: &Apilamiento, posicion: usize, umbrales: &Umbrales) -> f64 {
    let observadas: u32 = apilamiento.bases[posicion].iter().sum();
    if observadas == 0 {
        return probabilidad_de_error(umbrales.calidad_base_minima);
    }
    apilamiento.suma_error[posicion] as f64 / observadas as f64
}

pub fn llamar_indels(apilamiento: &Apilamiento, umbrales: &Umbrales) -> Vec<Llamada> {
    let mut llamadas = Vec::new();
    for (cambio, &a_favor) in &apilamiento.indels {
        let Some(&cruzan) = apilamiento.cruzan.get(cambio.posicion) else {
            continue;
        };
        let profundidad = cruzan.max(a_favor);
        let fraccion = a_favor as f64 / profundidad as f64;
        if profundidad < umbrales.profundidad_minima || fraccion < umbrales.fraccion_minima_indel {
            continue;
        }
        let error = error_medio(apilamiento, cambio.posicion, umbrales);
        llamadas.push(Llamada {
            cambio: cambio.clone(),
            profundidad,
            a_favor,
            fraccion,
            calidad: calidad_phred(a_favor, profundidad, error),
            anotaciones: Vec::new(),
        });
    }
    llamadas.sort_by(|una, otra| {
        una.cambio
            .posicion
            .cmp(&otra.cambio.posicion)
            .then(otra.a_favor.cmp(&una.a_favor))
            .then(una.cambio.cmp(&otra.cambio))
    });
    llamadas.dedup_by_key(|llamada| llamada.cambio.posicion);
    llamadas
}

pub fn calidad_phred(a_favor: u32, profundidad: u32, probabilidad_error: f64) -> f64 {
    let posible = probabilidad_error < 1.0;
    if a_favor == 0 || !posible {
        return 0.0;
    }
    let p = probabilidad_error.max(1e-12);
    let n = profundidad.max(a_favor);
    let razon = (p / (1.0 - p)).ln();
    let mut logaritmo_termino = logaritmo_combinaciones(n, a_favor)
        + a_favor as f64 * p.ln()
        + (n - a_favor) as f64 * (1.0 - p).ln();
    let mut logaritmo_cola = logaritmo_termino;
    for exitos in a_favor..n {
        logaritmo_termino += ((n - exitos) as f64 / (exitos + 1) as f64).ln() + razon;
        logaritmo_cola = sumar_en_logaritmos(logaritmo_cola, logaritmo_termino);
    }
    (-10.0 * logaritmo_cola / std::f64::consts::LN_10).max(0.0)
}

fn logaritmo_combinaciones(n: u32, k: u32) -> f64 {
    let menor = k.min(n - k);
    let mut suma = 0.0;
    for j in 1..=menor {
        suma += ((n - menor + j) as f64 / j as f64).ln();
    }
    suma
}

fn sumar_en_logaritmos(uno: f64, otro: f64) -> f64 {
    let mayor = uno.max(otro);
    let menor = uno.min(otro);
    mayor + (menor - mayor).exp().ln_1p()
}

fn resumir_profundidad(
    apilamiento: &Apilamiento,
    referencia: &[u8],
    umbrales: &Umbrales,
) -> ResumenDeProfundidad {
    let mut resumen = ResumenDeProfundidad::default();
    let mut suma: u64 = 0;
    for (posicion, &base_referencia) in referencia.iter().enumerate().take(apilamiento.bases.len())
    {
        let profundidad = apilamiento.profundidad(posicion);
        if profundidad == 0 {
            continue;
        }
        resumen.con_lecturas += 1;
        suma += profundidad as u64;
        if profundidad >= umbrales.profundidad_minima {
            resumen.suficientes += 1;
            continue;
        }
        let Some(indice_referencia) = indice_base(base_referencia) else {
            continue;
        };
        let (_, a_favor) =
            alternativa_mas_frecuente(&apilamiento.bases[posicion], indice_referencia);
        if a_favor >= 2 && a_favor as f64 / profundidad as f64 >= umbrales.fraccion_minima {
            resumen.sin_confirmar += 1;
        }
    }
    if resumen.con_lecturas > 0 {
        resumen.media = suma as f64 / resumen.con_lecturas as f64;
    }
    resumen
}

#[derive(Debug, Clone, PartialEq)]
pub struct Vecino {
    pub nombre: String,
    pub distancia: usize,
    pub corriente_arriba: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Efecto {
    Intergenica {
        izquierda: Option<Vecino>,
        derecha: Option<Vecino>,
    },
    NoCodificante,
    Sinonima,
    CambioDeAminoacido,
    ParadaPrematura,
    PerdidaDeParada,
    PerdidaDeInicio,
    CorrimientoDeMarco,
    InsercionEnMarco,
    DelecionEnMarco,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Anotacion {
    pub nombre_gen: String,
    pub efecto: Efecto,
    pub proteina: String,
    pub resistencia: Option<&'static str>,
}

pub struct IndiceDeGenes<'a> {
    genes: Vec<&'a Gen>,
    largo_maximo: usize,
}

impl<'a> IndiceDeGenes<'a> {
    pub fn nuevo(genes: &'a [Gen], longitud_referencia: usize) -> IndiceDeGenes<'a> {
        let mut ordenados: Vec<&Gen> = genes
            .iter()
            .filter(|candidato| {
                candidato.inicio >= 1
                    && candidato.inicio <= candidato.fin
                    && candidato.fin <= longitud_referencia
            })
            .collect();
        ordenados.sort_by_key(|candidato| candidato.inicio);
        let largo_maximo = ordenados
            .iter()
            .map(|candidato| candidato.fin - candidato.inicio + 1)
            .max()
            .unwrap_or(0);
        IndiceDeGenes {
            genes: ordenados,
            largo_maximo,
        }
    }

    pub fn esta_vacio(&self) -> bool {
        self.genes.is_empty()
    }

    fn que_tocan(&self, desde: usize, hasta: usize) -> Vec<&'a Gen> {
        let empiezan_antes = self
            .genes
            .partition_point(|candidato| candidato.inicio - 1 <= hasta);
        let mut tocados = Vec::new();
        for &candidato in self.genes[..empiezan_antes].iter().rev() {
            let ultima_posible = candidato.inicio - 1 + self.largo_maximo - 1;
            if ultima_posible < desde {
                break;
            }
            let ultima = candidato.fin - 1;
            if ultima >= desde {
                tocados.push(candidato);
            }
        }
        tocados.reverse();
        tocados
    }

    fn vecinos(&self, desde: usize, hasta: usize) -> (Option<&'a Gen>, Option<&'a Gen>) {
        let empiezan_antes = self
            .genes
            .partition_point(|candidato| candidato.inicio - 1 <= hasta);
        let derecho = self.genes.get(empiezan_antes).copied();
        let mut izquierdo: Option<&'a Gen> = None;
        for &candidato in self.genes[..empiezan_antes].iter().rev() {
            let fin_mas_tardio_posible = candidato.inicio + self.largo_maximo - 1;
            if izquierdo.is_some_and(|mejor| fin_mas_tardio_posible < mejor.fin) {
                break;
            }
            let ultima = candidato.fin - 1;
            if ultima < desde && izquierdo.is_none_or(|mejor| candidato.fin > mejor.fin) {
                izquierdo = Some(candidato);
            }
        }
        (izquierdo, derecho)
    }
}

pub fn anotar(cambio: &Cambio, referencia: &[u8], indice: &IndiceDeGenes) -> Vec<Anotacion> {
    if indice.esta_vacio() || cambio.referencia.is_empty() || cambio.alternativa.is_empty() {
        return Vec::new();
    }
    let (desde, hasta) = huella(cambio);
    let tocados = indice.que_tocan(desde, hasta);
    if tocados.is_empty() {
        return vec![anotar_intergenica(desde, hasta, indice)];
    }
    tocados
        .into_iter()
        .map(|elegido| anotar_en_gen(cambio, desde, hasta, referencia, elegido))
        .collect()
}

fn huella(cambio: &Cambio) -> (usize, usize) {
    let largo_referencia = cambio.referencia.len();
    let largo_alternativa = cambio.alternativa.len();
    let ultima = cambio.posicion + largo_referencia - 1;
    if largo_alternativa > largo_referencia {
        (ultima + 1, ultima)
    } else if largo_referencia > largo_alternativa {
        (cambio.posicion + largo_alternativa, ultima)
    } else {
        (cambio.posicion, ultima)
    }
}

fn anotar_intergenica(desde: usize, hasta: usize, indice: &IndiceDeGenes) -> Anotacion {
    let (izquierdo, derecho) = indice.vecinos(desde, hasta);
    let izquierda = izquierdo.map(|elegido| Vecino {
        nombre: elegido.nombre.clone(),
        distancia: desde - (elegido.fin - 1),
        corriente_arriba: elegido.hebra == b'-',
    });
    let derecha = derecho.map(|elegido| Vecino {
        nombre: elegido.nombre.clone(),
        distancia: (elegido.inicio - 1) - hasta,
        corriente_arriba: elegido.hebra == b'+',
    });
    Anotacion {
        nombre_gen: String::new(),
        efecto: Efecto::Intergenica { izquierda, derecha },
        proteina: String::new(),
        resistencia: None,
    }
}

struct Consecuencia {
    efecto: Efecto,
    primer_codon: usize,
    ultimo_codon: usize,
    proteina: String,
}

fn anotar_en_gen(
    cambio: &Cambio,
    desde: usize,
    hasta: usize,
    referencia: &[u8],
    elegido: &Gen,
) -> Anotacion {
    let consecuencia = if elegido.proteina.is_empty() {
        Consecuencia {
            efecto: Efecto::NoCodificante,
            primer_codon: 0,
            ultimo_codon: 0,
            proteina: String::new(),
        }
    } else if cambio.referencia.len() == cambio.alternativa.len() {
        efecto_de_sustitucion(cambio, desde, referencia, elegido)
    } else {
        efecto_de_indel(cambio, desde, hasta, referencia, elegido)
    };
    Anotacion {
        nombre_gen: elegido.nombre.clone(),
        resistencia: resistencia(&elegido.nombre, &consecuencia),
        efecto: consecuencia.efecto,
        proteina: consecuencia.proteina,
    }
}

fn indice_en_el_gen(elegido: &Gen, posicion: usize) -> usize {
    if elegido.hebra == b'-' {
        (elegido.fin - 1) - posicion
    } else {
        posicion - (elegido.inicio - 1)
    }
}

fn posicion_en_genoma(elegido: &Gen, indice: usize) -> usize {
    if elegido.hebra == b'-' {
        (elegido.fin - 1) - indice
    } else {
        (elegido.inicio - 1) + indice
    }
}

pub fn codon(
    referencia: &[u8],
    elegido: &Gen,
    numero: usize,
    cambio: Option<&Cambio>,
) -> Option<Vec<u8>> {
    let largo_gen = elegido.fin - elegido.inicio + 1;
    if numero == 0 || numero * 3 > largo_gen {
        return None;
    }
    let mut bases = Vec::with_capacity(3);
    for desplazamiento in 0..3 {
        let posicion = posicion_en_genoma(elegido, (numero - 1) * 3 + desplazamiento);
        let mut letra = referencia.get(posicion).copied().unwrap_or(b'N');
        if let Some(cambio) = cambio
            && posicion >= cambio.posicion
            && posicion < cambio.posicion + cambio.alternativa.len()
        {
            letra = cambio.alternativa[posicion - cambio.posicion];
        }
        if elegido.hebra == b'-' {
            bases.push(complemento(letra));
        } else {
            bases.push(letra);
        }
    }
    Some(bases)
}

fn es_codon_de_inicio(bases: &[u8]) -> bool {
    matches!(bases, b"ATG" | b"GTG" | b"TTG")
}

fn aminoacido_en(referencia: &[u8], elegido: &Gen, numero: usize) -> u8 {
    match codon(referencia, elegido, numero, None) {
        Some(bases) if numero == 1 && es_codon_de_inicio(&bases) => b'M',
        Some(bases) => aminoacido(&bases),
        None => b'X',
    }
}

fn residuo(referencia: &[u8], elegido: &Gen, numero: usize) -> String {
    format!(
        "{}{numero}",
        aminoacido_en(referencia, elegido, numero) as char
    )
}

fn efecto_de_sustitucion(
    cambio: &Cambio,
    desde: usize,
    referencia: &[u8],
    elegido: &Gen,
) -> Consecuencia {
    let dentro = desde.clamp(elegido.inicio - 1, elegido.fin - 1);
    let numero = indice_en_el_gen(elegido, dentro) / 3 + 1;
    let (Some(antes), Some(despues)) = (
        codon(referencia, elegido, numero, None),
        codon(referencia, elegido, numero, Some(cambio)),
    ) else {
        return Consecuencia {
            efecto: Efecto::NoCodificante,
            primer_codon: numero,
            ultimo_codon: numero,
            proteina: String::new(),
        };
    };
    let era_inicio = numero == 1 && es_codon_de_inicio(&antes);
    let sigue_siendo_inicio = era_inicio && es_codon_de_inicio(&despues);
    let aminoacido_antes = aminoacido_en(referencia, elegido, numero);
    let aminoacido_despues = if sigue_siendo_inicio {
        b'M'
    } else {
        aminoacido(&despues)
    };
    let efecto = if era_inicio && !sigue_siendo_inicio {
        Efecto::PerdidaDeInicio
    } else if aminoacido_antes == aminoacido_despues {
        Efecto::Sinonima
    } else if aminoacido_despues == b'*' {
        Efecto::ParadaPrematura
    } else if aminoacido_antes == b'*' {
        Efecto::PerdidaDeParada
    } else {
        Efecto::CambioDeAminoacido
    };
    Consecuencia {
        efecto,
        primer_codon: numero,
        ultimo_codon: numero,
        proteina: format!(
            "{}{numero}{}",
            aminoacido_antes as char, aminoacido_despues as char
        ),
    }
}

fn efecto_de_indel(
    cambio: &Cambio,
    desde: usize,
    hasta: usize,
    referencia: &[u8],
    elegido: &Gen,
) -> Consecuencia {
    let (primera_movida, primer_codon, ultimo_codon) = if cambio.es_insercion() {
        let despues_de = indice_en_el_gen(elegido, desde).min(indice_en_el_gen(elegido, hasta));
        (despues_de + 1, despues_de / 3 + 1, despues_de / 3 + 2)
    } else {
        let uno = indice_en_el_gen(elegido, desde.max(elegido.inicio - 1));
        let otro = indice_en_el_gen(elegido, hasta.min(elegido.fin - 1));
        (uno.min(otro), uno.min(otro) / 3 + 1, uno.max(otro) / 3 + 1)
    };
    let primera_movida = if primera_movida < 3 && conserva_el_inicio(cambio, referencia, elegido) {
        3
    } else {
        primera_movida
    };
    let diferencia = cambio.alternativa.len().abs_diff(cambio.referencia.len());
    let (efecto, proteina) = if primera_movida < 3 {
        (Efecto::PerdidaDeInicio, "M1?".to_string())
    } else if !diferencia.is_multiple_of(3) {
        let numero = primera_movida / 3 + 1;
        (
            Efecto::CorrimientoDeMarco,
            format!("{}fs", residuo(referencia, elegido, numero)),
        )
    } else if cambio.es_insercion() {
        (
            Efecto::InsercionEnMarco,
            format!(
                "{}_{}ins",
                residuo(referencia, elegido, primer_codon),
                residuo(referencia, elegido, ultimo_codon)
            ),
        )
    } else if primer_codon == ultimo_codon {
        (
            Efecto::DelecionEnMarco,
            format!("{}del", residuo(referencia, elegido, primer_codon)),
        )
    } else {
        (
            Efecto::DelecionEnMarco,
            format!(
                "{}_{}del",
                residuo(referencia, elegido, primer_codon),
                residuo(referencia, elegido, ultimo_codon)
            ),
        )
    };
    Consecuencia {
        efecto,
        primer_codon,
        ultimo_codon,
        proteina,
    }
}

fn conserva_el_inicio(cambio: &Cambio, referencia: &[u8], elegido: &Gen) -> bool {
    let posicion = cambio.posicion;
    let largo_referencia = cambio.referencia.len();
    let largo_alternativa = cambio.alternativa.len();
    if posicion + largo_referencia > referencia.len() {
        return false;
    }
    let primera = elegido.inicio - 1;
    let desde = primera.min(posicion);
    let hasta = elegido.fin.max(posicion + largo_referencia);
    let mut ventana = referencia[desde..posicion].to_vec();
    ventana.extend_from_slice(&cambio.alternativa);
    ventana
        .extend_from_slice(&referencia[posicion + largo_referencia..hasta.min(referencia.len())]);
    let corrida = |lugar: usize, ultimo_conservado: usize| {
        if lugar <= ultimo_conservado {
            lugar - desde
        } else {
            let desplazado = (lugar + largo_alternativa).saturating_sub(largo_referencia + desde);
            desplazado.max(posicion - desde + largo_alternativa)
        }
    };
    let inicio_mutado = corrida(primera, posicion);
    let fin_mutado = corrida(elegido.fin, posicion + 1);
    if inicio_mutado + 3 > fin_mutado || fin_mutado > ventana.len() {
        return false;
    }
    let gen_mutado = &ventana[inicio_mutado..fin_mutado];
    let primer_codon = if elegido.hebra == b'-' {
        complemento_inverso(&gen_mutado[gen_mutado.len() - 3..])
    } else {
        gen_mutado[..3].to_vec()
    };
    es_codon_de_inicio(&primer_codon)
}

fn resistencia(nombre_gen: &str, consecuencia: &Consecuencia) -> Option<&'static str> {
    let cambia_la_proteina = matches!(
        consecuencia.efecto,
        Efecto::CambioDeAminoacido | Efecto::InsercionEnMarco | Efecto::DelecionEnMarco
    );
    if !cambia_la_proteina {
        return None;
    }
    for &(gen_resistente, desde, hasta, antibiotico) in SITIOS_DE_RESISTENCIA {
        let toca_el_sitio =
            consecuencia.primer_codon <= hasta && consecuencia.ultimo_codon >= desde;
        if nombre_gen == gen_resistente && toca_el_sitio {
            return Some(antibiotico);
        }
    }
    None
}

pub fn es_transicion(referencia: u8, alternativa: u8) -> bool {
    matches!(
        (referencia, alternativa),
        (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C')
    )
}

pub fn escribir_llamadas_vcf(
    destino: &mut impl Write,
    cromosoma: &str,
    longitud: usize,
    llamadas: &[Llamada],
) -> io::Result<()> {
    formatos::escribir_cabecera_vcf(destino, cromosoma, longitud)?;
    for llamada in llamadas {
        let variante = Variante {
            cromosoma: cromosoma.to_string(),
            posicion: llamada.cambio.posicion + 1,
            referencia: llamada.cambio.referencia.clone(),
            alternativa: llamada.cambio.alternativa.clone(),
            calidad: llamada.calidad,
            profundidad: llamada.profundidad,
            frecuencia: llamada.fraccion,
        };
        formatos::escribir_vcf(destino, &variante)?;
    }
    Ok(())
}

pub fn leer_cambios_vcf(fuente: impl BufRead) -> Result<Vec<Cambio>, String> {
    let mut cambios = Vec::new();
    for (indice, linea) in fuente.split(b'\n').enumerate() {
        let linea = linea.map_err(|error| format!("No se pudo leer el VCF: {error}"))?;
        let numero = indice + 1;
        let Some(variante) = formatos::leer_vcf(&linea)
            .map_err(|error| format!("Línea {numero} del VCF: {error}"))?
        else {
            continue;
        };
        if variante.posicion == 0 || !son_bases(&variante.referencia) {
            return Err(format!(
                "Línea {numero} del VCF: la posición o la base de referencia no son válidas."
            ));
        }
        for alternativa in variante.alternativa.split(|&letra| letra == b',') {
            if son_acgt(alternativa) && alternativa != &variante.referencia[..] {
                cambios.push(Cambio {
                    posicion: variante.posicion - 1,
                    referencia: variante.referencia.clone(),
                    alternativa: alternativa.to_vec(),
                });
            }
        }
    }
    Ok(cambios)
}

fn son_acgt(alelo: &[u8]) -> bool {
    !alelo.is_empty() && alelo.iter().all(|letra| BASES.contains(letra))
}

fn son_bases(alelo: &[u8]) -> bool {
    !alelo.is_empty() && alelo.iter().all(|letra| letra.is_ascii_alphabetic())
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Comparacion {
    pub verdaderos_positivos: Vec<Cambio>,
    pub falsos_positivos: Vec<Cambio>,
    pub falsos_negativos: Vec<Cambio>,
}

impl Comparacion {
    pub fn sensibilidad(&self) -> Option<f64> {
        let reales = self.verdaderos_positivos.len() + self.falsos_negativos.len();
        if reales == 0 {
            return None;
        }
        Some(self.verdaderos_positivos.len() as f64 / reales as f64)
    }

    pub fn precision(&self) -> Option<f64> {
        let llamadas = self.verdaderos_positivos.len() + self.falsos_positivos.len();
        if llamadas == 0 {
            return None;
        }
        Some(self.verdaderos_positivos.len() as f64 / llamadas as f64)
    }
}

pub fn comparar_con_verdad(
    llamadas: &[Cambio],
    verdad: &[Cambio],
    referencia: &[u8],
) -> Comparacion {
    let llamados: HashSet<Cambio> = llamadas
        .iter()
        .map(|cambio| normalizar(cambio, referencia))
        .collect();
    let verdaderos: HashSet<Cambio> = verdad
        .iter()
        .map(|cambio| normalizar(cambio, referencia))
        .collect();
    let mut comparacion = Comparacion::default();
    for cambio in &llamados {
        if verdaderos.contains(cambio) {
            comparacion.verdaderos_positivos.push(cambio.clone());
        } else {
            comparacion.falsos_positivos.push(cambio.clone());
        }
    }
    for cambio in &verdaderos {
        if !llamados.contains(cambio) {
            comparacion.falsos_negativos.push(cambio.clone());
        }
    }
    comparacion.verdaderos_positivos.sort();
    comparacion.falsos_positivos.sort();
    comparacion.falsos_negativos.sort();
    comparacion
}
