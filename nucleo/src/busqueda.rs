use crate::ecoli;
use crate::formatos::Registro;
use crate::matrices::{AMINOACIDOS, blosum62};
use crate::secuencia::{complemento_inverso, indice_base, traducir};
use rayon::prelude::*;
use std::collections::HashMap;

pub const LONGITUD_MAXIMA_PROTEINA: usize = 10_000;
pub const LONGITUD_MAXIMA_ADN: usize = 100_000;
pub const LONGITUD_MAXIMA_BLASTX: usize = 30_000;
pub const LETRAS_MAXIMAS_BASE: usize = 100_000_000;
pub const LIMITE_DE_GOLPES: u64 = 200_000_000;
pub const LIMITE_DE_GOLPES_ADN: u64 = 2_000_000;
pub const LIMITE_DE_TRAMOS: usize = 1_000;
pub const LONGITUD_MAXIMA_SUJETO_PROTEINA: usize = 2_000_000;

pub const PAREJA: u8 = b'M';
pub const HUECO_EN_SUJETO: u8 = b'I';
pub const HUECO_EN_CONSULTA: u8 = b'D';

const LONGITUD_PALABRA_PROTEINA: usize = 3;
const LONGITUD_PALABRA_ADN: usize = 11;
const PALABRAS_DE_PROTEINA: usize = 20 * 20 * 20;
const PALABRAS_DE_ADN: usize = 1 << (2 * LONGITUD_PALABRA_ADN);
const SIN_PALABRA: u16 = u16::MAX;
const SIN_CODIGO: u32 = u32::MAX;
const FUERA_DEL_ALFABETO: u8 = u8::MAX;
const NUNCA: i64 = i64::MIN / 2;
const MENOS_INFINITO: i32 = i32::MIN / 4;
const VENTANA_DE_SEMILLA: usize = 11;
const LIMITE_DE_CELDAS: usize = 50_000_000;

const ORIGEN: u8 = 3;
const DESDE_DIAGONAL: u8 = 0;
const DESDE_HORIZONTAL: u8 = 1;
const DESDE_VERTICAL: u8 = 2;
const SIGUE_HORIZONTAL: u8 = 4;
const SIGUE_VERTICAL: u8 = 8;

#[derive(Debug, Clone)]
pub struct Puntuacion {
    filas: Vec<[i32; 256]>,
}

impl Puntuacion {
    pub fn blosum62() -> Puntuacion {
        let mut filas = vec![[-4; 256]; 256];
        for primera in b'A'..=b'Z' {
            for segunda in b'A'..=b'Z' {
                filas[primera as usize][segunda as usize] = blosum62(primera, segunda);
            }
        }
        Puntuacion { filas }
    }

    pub fn adn(premio: i32, castigo: i32) -> Puntuacion {
        let mut filas = vec![[castigo; 256]; 256];
        for &base in b"ACGT" {
            filas[base as usize][base as usize] = premio;
        }
        Puntuacion { filas }
    }

    pub fn puntaje(&self, primera: u8, segunda: u8) -> i32 {
        self.filas[primera as usize][segunda as usize]
    }

    fn fila(&self, primera: u8) -> &[i32; 256] {
        &self.filas[primera as usize]
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Estadistica {
    pub lambda: f64,
    pub k: f64,
    pub h: f64,
}

#[derive(Debug, Clone)]
pub struct Parametros {
    pub adn: bool,
    pub longitud_palabra: usize,
    pub umbral_vecinas: i32,
    pub ventana_dos_golpes: usize,
    pub caida_sin_huecos: i32,
    pub disparo_con_huecos: i32,
    pub caida_preliminar: i32,
    pub caida_con_huecos: i32,
    pub apertura: i32,
    pub extension: i32,
    pub puntuacion: Puntuacion,
    pub estadistica: Estadistica,
}

pub fn parametros_blastp() -> Parametros {
    Parametros {
        adn: false,
        longitud_palabra: LONGITUD_PALABRA_PROTEINA,
        umbral_vecinas: 11,
        ventana_dos_golpes: 40,
        caida_sin_huecos: 16,
        disparo_con_huecos: 41,
        caida_preliminar: 38,
        caida_con_huecos: 64,
        apertura: 11,
        extension: 1,
        puntuacion: Puntuacion::blosum62(),
        estadistica: Estadistica {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
        },
    }
}

pub fn parametros_blastn() -> Parametros {
    Parametros {
        adn: true,
        longitud_palabra: LONGITUD_PALABRA_ADN,
        umbral_vecinas: 0,
        ventana_dos_golpes: 0,
        caida_sin_huecos: 22,
        disparo_con_huecos: 28,
        caida_preliminar: 33,
        caida_con_huecos: 110,
        apertura: 5,
        extension: 2,
        puntuacion: Puntuacion::adn(2, -3),
        estadistica: Estadistica {
            lambda: 0.625,
            k: 0.41,
            h: 0.78,
        },
    }
}

pub fn bits(puntaje: i32, estadistica: &Estadistica) -> f64 {
    (estadistica.lambda * puntaje as f64 - estadistica.k.ln()) / std::f64::consts::LN_2
}

pub fn evalor(bits: f64, espacio: f64) -> f64 {
    espacio * 2f64.powf(-bits)
}

pub fn puntaje_necesario(evalor_maximo: f64, espacio: f64, estadistica: &Estadistica) -> i32 {
    ((estadistica.k * espacio / evalor_maximo).ln() / estadistica.lambda).ceil() as i32
}

pub fn ajuste_de_longitud(
    longitud_consulta: f64,
    longitud_base: f64,
    secuencias: f64,
    estadistica: &Estadistica,
) -> f64 {
    let ajuste = (estadistica.k * longitud_consulta * longitud_base).ln() / estadistica.h;
    let tope_consulta = longitud_consulta - 1.0 / estadistica.k;
    let tope_base = (longitud_base - 1.0 / estadistica.k) / secuencias;
    ajuste.min(tope_consulta).min(tope_base).max(0.0)
}

pub fn espacio_de_busqueda(
    longitud_consulta: f64,
    longitud_base: f64,
    secuencias: f64,
    estadistica: &Estadistica,
) -> f64 {
    let ajuste = ajuste_de_longitud(longitud_consulta, longitud_base, secuencias, estadistica);
    (longitud_consulta - ajuste) * (longitud_base - secuencias * ajuste)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Tramo {
    pub inicio_consulta: usize,
    pub inicio_sujeto: usize,
    pub longitud: usize,
    pub puntaje: i32,
}

pub fn extender_sin_huecos(
    consulta: &[u8],
    sujeto: &[u8],
    posicion_consulta: usize,
    posicion_sujeto: usize,
    longitud_palabra: usize,
    puntuacion: &Puntuacion,
    caida: i32,
) -> Tramo {
    if posicion_consulta + longitud_palabra > consulta.len()
        || posicion_sujeto + longitud_palabra > sujeto.len()
    {
        return Tramo {
            inicio_consulta: posicion_consulta,
            inicio_sujeto: posicion_sujeto,
            longitud: 0,
            puntaje: 0,
        };
    }
    let mut puntaje_palabra = 0;
    for paso in 0..longitud_palabra {
        puntaje_palabra += puntuacion.puntaje(
            consulta[posicion_consulta + paso],
            sujeto[posicion_sujeto + paso],
        );
    }

    let derecha_consulta = posicion_consulta + longitud_palabra;
    let derecha_sujeto = posicion_sujeto + longitud_palabra;
    let mut suma = 0;
    let mut mejor_derecha = 0;
    let mut largo_derecha = 0;
    let mut paso = 0;
    while derecha_consulta + paso < consulta.len() && derecha_sujeto + paso < sujeto.len() {
        suma += puntuacion.puntaje(
            consulta[derecha_consulta + paso],
            sujeto[derecha_sujeto + paso],
        );
        paso += 1;
        if suma > mejor_derecha {
            mejor_derecha = suma;
            largo_derecha = paso;
        } else if suma < mejor_derecha - caida {
            break;
        }
    }

    let mut suma = 0;
    let mut mejor_izquierda = 0;
    let mut largo_izquierda = 0;
    let mut paso = 0;
    while paso < posicion_consulta && paso < posicion_sujeto {
        paso += 1;
        suma += puntuacion.puntaje(
            consulta[posicion_consulta - paso],
            sujeto[posicion_sujeto - paso],
        );
        if suma > mejor_izquierda {
            mejor_izquierda = suma;
            largo_izquierda = paso;
        } else if suma < mejor_izquierda - caida {
            break;
        }
    }

    Tramo {
        inicio_consulta: posicion_consulta - largo_izquierda,
        inicio_sujeto: posicion_sujeto - largo_izquierda,
        longitud: largo_izquierda + longitud_palabra + largo_derecha,
        puntaje: puntaje_palabra + mejor_izquierda + mejor_derecha,
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlineamientoLocal {
    pub inicio_consulta: usize,
    pub fin_consulta: usize,
    pub inicio_sujeto: usize,
    pub fin_sujeto: usize,
    pub puntaje: i32,
    pub operaciones: Vec<u8>,
}

pub fn alinear_con_huecos(
    consulta: &[u8],
    sujeto: &[u8],
    semilla_consulta: usize,
    semilla_sujeto: usize,
    parametros: &Parametros,
    caida: i32,
) -> AlineamientoLocal {
    let semilla_consulta = semilla_consulta.min(consulta.len());
    let semilla_sujeto = semilla_sujeto.min(sujeto.len());
    let izquierda = extender_con_huecos(
        &consulta[..semilla_consulta],
        &sujeto[..semilla_sujeto],
        true,
        parametros,
        caida,
    );
    let derecha = extender_con_huecos(
        &consulta[semilla_consulta..],
        &sujeto[semilla_sujeto..],
        false,
        parametros,
        caida,
    );
    let mut operaciones = izquierda.operaciones;
    operaciones.reverse();
    operaciones.extend(derecha.operaciones);
    AlineamientoLocal {
        inicio_consulta: semilla_consulta - izquierda.usadas_consulta,
        fin_consulta: semilla_consulta + derecha.usadas_consulta,
        inicio_sujeto: semilla_sujeto - izquierda.usadas_sujeto,
        fin_sujeto: semilla_sujeto + derecha.usadas_sujeto,
        puntaje: izquierda.puntaje + derecha.puntaje,
        operaciones,
    }
}

struct Extension {
    puntaje: i32,
    usadas_consulta: usize,
    usadas_sujeto: usize,
    operaciones: Vec<u8>,
}

#[derive(Clone, Copy, PartialEq)]
enum Capa {
    Principal,
    Horizontal,
    Vertical,
}

fn letra(secuencia: &[u8], indice: usize, hacia_atras: bool) -> u8 {
    if hacia_atras {
        secuencia[secuencia.len() - 1 - indice]
    } else {
        secuencia[indice]
    }
}

fn elegir_mejor(diagonal: i32, horizontal: i32, vertical: i32) -> (i32, u8) {
    if diagonal >= horizontal && diagonal >= vertical {
        (diagonal, DESDE_DIAGONAL)
    } else if horizontal >= vertical {
        (horizontal, DESDE_HORIZONTAL)
    } else {
        (vertical, DESDE_VERTICAL)
    }
}

fn extender_con_huecos(
    consulta: &[u8],
    sujeto: &[u8],
    hacia_atras: bool,
    parametros: &Parametros,
    caida: i32,
) -> Extension {
    let abrir = parametros.apertura + parametros.extension;
    let alargar = parametros.extension;
    let mut mejor = 0;
    let mut mejor_fila = 0;
    let mut mejor_columna = 0;

    let mut rastros: Vec<u8> = vec![DESDE_DIAGONAL];
    let mut primera_columna: Vec<usize> = vec![0];
    let mut comienzo_de_fila: Vec<usize> = vec![0];
    let mut principal_anterior: Vec<i32> = vec![0];
    let mut vertical_anterior: Vec<i32> = vec![MENOS_INFINITO];
    let mut columna = 1;
    while columna <= sujeto.len() {
        let puntaje = -(abrir + alargar * (columna as i32 - 1));
        if puntaje < -caida {
            break;
        }
        principal_anterior.push(puntaje);
        vertical_anterior.push(MENOS_INFINITO);
        if columna == 1 {
            rastros.push(DESDE_HORIZONTAL);
        } else {
            rastros.push(DESDE_HORIZONTAL | SIGUE_HORIZONTAL);
        }
        columna += 1;
    }

    let mut inicio_anterior = 0;
    let mut principal_actual: Vec<i32> = Vec::new();
    let mut vertical_actual: Vec<i32> = Vec::new();
    for fila in 1..=consulta.len() {
        if rastros.len() > LIMITE_DE_CELDAS {
            break;
        }
        let puntajes = parametros
            .puntuacion
            .fila(letra(consulta, fila - 1, hacia_atras));
        let fin_anterior = inicio_anterior + principal_anterior.len();
        primera_columna.push(inicio_anterior);
        comienzo_de_fila.push(rastros.len());
        principal_actual.clear();
        vertical_actual.clear();
        let mut principal_izquierda = MENOS_INFINITO;
        let mut horizontal = MENOS_INFINITO;
        let mut columna = inicio_anterior;
        while columna <= sujeto.len() {
            let mut banderas = 0;
            let diagonal = if columna > inicio_anterior && columna <= fin_anterior {
                let letra_sujeto = letra(sujeto, columna - 1, hacia_atras);
                principal_anterior[columna - 1 - inicio_anterior] + puntajes[letra_sujeto as usize]
            } else {
                MENOS_INFINITO
            };
            let mut vertical = MENOS_INFINITO;
            if columna < fin_anterior {
                let alargando = vertical_anterior[columna - inicio_anterior] - alargar;
                let abriendo = principal_anterior[columna - inicio_anterior] - abrir;
                if alargando >= abriendo {
                    banderas |= SIGUE_VERTICAL;
                    vertical = alargando;
                } else {
                    vertical = abriendo;
                }
            }
            let alargando = horizontal - alargar;
            let abriendo = principal_izquierda - abrir;
            if alargando >= abriendo {
                banderas |= SIGUE_HORIZONTAL;
                horizontal = alargando;
            } else {
                horizontal = abriendo;
            }

            let (mut principal, origen) = elegir_mejor(diagonal, horizontal, vertical);
            if principal < mejor - caida {
                principal = MENOS_INFINITO;
                horizontal = MENOS_INFINITO;
                vertical = MENOS_INFINITO;
            } else if principal > mejor {
                mejor = principal;
                mejor_fila = fila;
                mejor_columna = columna;
            }
            rastros.push(origen | banderas);
            principal_actual.push(principal);
            vertical_actual.push(vertical);
            principal_izquierda = principal;
            if columna >= fin_anterior && principal == MENOS_INFINITO {
                break;
            }
            columna += 1;
        }

        let Some(primera_viva) = principal_actual
            .iter()
            .position(|&valor| valor > MENOS_INFINITO)
        else {
            break;
        };
        let ultima_viva = principal_actual
            .iter()
            .rposition(|&valor| valor > MENOS_INFINITO)
            .unwrap_or(primera_viva);
        principal_anterior.clear();
        principal_anterior.extend_from_slice(&principal_actual[primera_viva..=ultima_viva]);
        vertical_anterior.clear();
        vertical_anterior.extend_from_slice(&vertical_actual[primera_viva..=ultima_viva]);
        inicio_anterior += primera_viva;
    }

    let mut operaciones = Vec::new();
    let mut fila = mejor_fila;
    let mut columna = mejor_columna;
    let mut capa = Capa::Principal;
    while fila > 0 || columna > 0 {
        let codigo = rastros[comienzo_de_fila[fila] + columna - primera_columna[fila]];
        match capa {
            Capa::Principal => match codigo & ORIGEN {
                DESDE_HORIZONTAL => capa = Capa::Horizontal,
                DESDE_VERTICAL => capa = Capa::Vertical,
                _ => {
                    operaciones.push(PAREJA);
                    fila -= 1;
                    columna -= 1;
                }
            },
            Capa::Horizontal => {
                operaciones.push(HUECO_EN_CONSULTA);
                if codigo & SIGUE_HORIZONTAL == 0 {
                    capa = Capa::Principal;
                }
                columna -= 1;
            }
            Capa::Vertical => {
                operaciones.push(HUECO_EN_SUJETO);
                if codigo & SIGUE_VERTICAL == 0 {
                    capa = Capa::Principal;
                }
                fila -= 1;
            }
        }
    }
    operaciones.reverse();

    Extension {
        puntaje: mejor,
        usadas_consulta: mejor_fila,
        usadas_sujeto: mejor_columna,
        operaciones,
    }
}

struct Indice {
    inicio_por_palabra: Vec<u32>,
    posiciones: Vec<u32>,
}

impl Indice {
    fn nuevo(pares: &[(u32, u32)], cantidad_de_palabras: usize) -> Indice {
        let mut inicio_por_palabra = vec![0u32; cantidad_de_palabras + 1];
        for &(palabra, _) in pares {
            inicio_por_palabra[palabra as usize] += 1;
        }
        let mut acumulado = 0;
        for inicio in inicio_por_palabra.iter_mut() {
            let cuantas = *inicio;
            *inicio = acumulado;
            acumulado += cuantas;
        }
        let mut siguiente = inicio_por_palabra.clone();
        let mut posiciones = vec![0u32; pares.len()];
        for &(palabra, posicion) in pares {
            posiciones[siguiente[palabra as usize] as usize] = posicion;
            siguiente[palabra as usize] += 1;
        }
        Indice {
            inicio_por_palabra,
            posiciones,
        }
    }

    fn posiciones_de(&self, palabra: u32) -> &[u32] {
        let desde = self.inicio_por_palabra[palabra as usize] as usize;
        let hasta = self.inicio_por_palabra[palabra as usize + 1] as usize;
        &self.posiciones[desde..hasta]
    }
}

fn indices_de_aminoacidos() -> [u8; 256] {
    let mut indices = [FUERA_DEL_ALFABETO; 256];
    for (indice, &aminoacido) in AMINOACIDOS.iter().enumerate() {
        indices[aminoacido as usize] = indice as u8;
    }
    indices
}

fn matriz_por_indices() -> [[i32; 20]; 20] {
    let mut matriz = [[0; 20]; 20];
    for (fila, &primera) in AMINOACIDOS.iter().enumerate() {
        for (columna, &segunda) in AMINOACIDOS.iter().enumerate() {
            matriz[fila][columna] = blosum62(primera, segunda);
        }
    }
    matriz
}

fn codigo_con_indices(palabra: &[u8], indices: &[u8; 256]) -> u16 {
    let mut codigo: u16 = 0;
    for &aminoacido in palabra {
        let indice = indices[aminoacido as usize];
        if indice == FUERA_DEL_ALFABETO {
            return SIN_PALABRA;
        }
        codigo = codigo * 20 + indice as u16;
    }
    codigo
}

pub fn codigo_de_palabra(palabra: &[u8]) -> Option<u16> {
    if palabra.len() != LONGITUD_PALABRA_PROTEINA {
        return None;
    }
    match codigo_con_indices(palabra, &indices_de_aminoacidos()) {
        SIN_PALABRA => None,
        codigo => Some(codigo),
    }
}

fn agregar_vecinas(
    letras: [usize; 3],
    umbral: i32,
    matriz: &[[i32; 20]; 20],
    vecinas: &mut Vec<u16>,
) {
    let [primera, segunda, tercera] = letras;
    let mejor_segunda = matriz[segunda].iter().copied().max().unwrap_or(0);
    let mejor_tercera = matriz[tercera].iter().copied().max().unwrap_or(0);
    for (x, &con_una) in matriz[primera].iter().enumerate() {
        if con_una + mejor_segunda + mejor_tercera < umbral {
            continue;
        }
        for (y, &puntaje_y) in matriz[segunda].iter().enumerate() {
            let con_dos = con_una + puntaje_y;
            if con_dos + mejor_tercera < umbral {
                continue;
            }
            for (z, &puntaje_z) in matriz[tercera].iter().enumerate() {
                if con_dos + puntaje_z >= umbral {
                    vecinas.push((x * 400 + y * 20 + z) as u16);
                }
            }
        }
    }
}

fn letras_de_palabra(palabra: &[u8], indices: &[u8; 256]) -> Option<[usize; 3]> {
    if codigo_con_indices(palabra, indices) == SIN_PALABRA {
        return None;
    }
    Some([
        indices[palabra[0] as usize] as usize,
        indices[palabra[1] as usize] as usize,
        indices[palabra[2] as usize] as usize,
    ])
}

pub fn palabras_vecinas(palabra: &[u8], umbral: i32) -> Vec<u16> {
    let mut vecinas = Vec::new();
    if palabra.len() != LONGITUD_PALABRA_PROTEINA {
        return vecinas;
    }
    if let Some(letras) = letras_de_palabra(palabra, &indices_de_aminoacidos()) {
        agregar_vecinas(letras, umbral, &matriz_por_indices(), &mut vecinas);
    }
    vecinas
}

fn tabla_de_vecinas(consulta: &[u8], umbral: i32) -> Indice {
    let indices = indices_de_aminoacidos();
    let matriz = matriz_por_indices();
    let mut pares: Vec<(u32, u32)> = Vec::new();
    let mut vecinas = Vec::new();
    for inicio in 0..consulta.len().saturating_sub(LONGITUD_PALABRA_PROTEINA - 1) {
        let palabra = &consulta[inicio..inicio + LONGITUD_PALABRA_PROTEINA];
        let Some(letras) = letras_de_palabra(palabra, &indices) else {
            continue;
        };
        vecinas.clear();
        agregar_vecinas(letras, umbral, &matriz, &mut vecinas);
        for &vecina in &vecinas {
            pares.push((vecina as u32, inicio as u32));
        }
    }
    Indice::nuevo(&pares, PALABRAS_DE_PROTEINA)
}

fn palabras_de_proteina(secuencia: &[u8], indices: &[u8; 256]) -> Vec<u16> {
    let mut palabras = Vec::new();
    for inicio in 0..secuencia
        .len()
        .saturating_sub(LONGITUD_PALABRA_PROTEINA - 1)
    {
        palabras.push(codigo_con_indices(
            &secuencia[inicio..inicio + LONGITUD_PALABRA_PROTEINA],
            indices,
        ));
    }
    palabras
}

pub struct BaseDeProteinas {
    pub registros: Vec<Registro>,
    pub total_letras: usize,
    palabras: Vec<Vec<u16>>,
    veces_por_palabra: Vec<u64>,
}

pub fn indexar_proteinas(registros: Vec<Registro>) -> Result<BaseDeProteinas, String> {
    let registros: Vec<Registro> = registros
        .into_iter()
        .filter(|registro| !registro.secuencia.is_empty())
        .collect();
    validar_base(&registros)?;
    if let Some(larga) = registros
        .iter()
        .find(|registro| registro.secuencia.len() > LONGITUD_MAXIMA_SUJETO_PROTEINA)
    {
        return Err(format!(
            "La proteína «{}» tiene {} letras y el máximo por proteína es {LONGITUD_MAXIMA_SUJETO_PROTEINA}. ¿Es una base de ADN?",
            larga.nombre,
            larga.secuencia.len()
        ));
    }
    Ok(construir_base_de_proteinas(registros))
}

fn validar_base(registros: &[Registro]) -> Result<(), String> {
    if registros.is_empty() {
        return Err("La base de datos no tiene ninguna secuencia.".to_string());
    }
    let total: usize = registros
        .iter()
        .map(|registro| registro.secuencia.len())
        .sum();
    if total > LETRAS_MAXIMAS_BASE {
        return Err(format!(
            "La base de datos tiene {total} letras y el máximo es {LETRAS_MAXIMAS_BASE}."
        ));
    }
    Ok(())
}

fn construir_base_de_proteinas(mut registros: Vec<Registro>) -> BaseDeProteinas {
    let indices = indices_de_aminoacidos();
    let mut veces_por_palabra = vec![0u64; PALABRAS_DE_PROTEINA];
    let mut palabras = Vec::with_capacity(registros.len());
    let mut total_letras = 0;
    for registro in &mut registros {
        registro.secuencia.make_ascii_uppercase();
        total_letras += registro.secuencia.len();
        let codigos = palabras_de_proteina(&registro.secuencia, &indices);
        for &codigo in &codigos {
            if codigo != SIN_PALABRA {
                veces_por_palabra[codigo as usize] += 1;
            }
        }
        palabras.push(codigos);
    }
    BaseDeProteinas {
        registros,
        total_letras,
        palabras,
        veces_por_palabra,
    }
}

pub fn proteinas_de_ecoli() -> BaseDeProteinas {
    let registros = ecoli::genes()
        .iter()
        .map(|anotado| Registro {
            nombre: format!(
                "{} ({}) {}",
                anotado.nombre, anotado.locus, anotado.producto
            ),
            secuencia: anotado.proteina.clone(),
        })
        .collect();
    construir_base_de_proteinas(registros)
}

pub fn codigos_de_adn(secuencia: &[u8]) -> Vec<u32> {
    let mut codigos = Vec::with_capacity(secuencia.len());
    let mascara = (PALABRAS_DE_ADN - 1) as u32;
    let mut codigo: u32 = 0;
    let mut seguidas = 0;
    for (posicion, &base) in secuencia.iter().enumerate() {
        match indice_base(base) {
            Some(indice) => {
                codigo = ((codigo << 2) | indice as u32) & mascara;
                seguidas += 1;
            }
            None => seguidas = 0,
        }
        if posicion + 1 >= LONGITUD_PALABRA_ADN {
            if seguidas >= LONGITUD_PALABRA_ADN {
                codigos.push(codigo);
            } else {
                codigos.push(SIN_CODIGO);
            }
        }
    }
    codigos
}

pub struct BaseDeAdn {
    pub nombres: Vec<String>,
    pub total_letras: usize,
    inicios: Vec<usize>,
    longitudes: Vec<usize>,
    secuencia: Vec<u8>,
    indice: Indice,
}

impl BaseDeAdn {
    pub fn cantidad(&self) -> usize {
        self.nombres.len()
    }

    pub fn registro(&self, numero: usize) -> &[u8] {
        match (self.inicios.get(numero), self.longitudes.get(numero)) {
            (Some(&inicio), Some(&longitud)) => &self.secuencia[inicio..inicio + longitud],
            _ => &[],
        }
    }

    fn registro_en(&self, posicion: usize) -> usize {
        self.inicios.partition_point(|&inicio| inicio <= posicion) - 1
    }
}

pub fn indexar_adn(registros: Vec<Registro>) -> Result<BaseDeAdn, String> {
    let registros: Vec<Registro> = registros
        .into_iter()
        .filter(|registro| !registro.secuencia.is_empty())
        .collect();
    validar_base(&registros)?;
    Ok(construir_base_de_adn(registros))
}

fn construir_base_de_adn(registros: Vec<Registro>) -> BaseDeAdn {
    let mut nombres = Vec::with_capacity(registros.len());
    let mut inicios = Vec::with_capacity(registros.len());
    let mut longitudes = Vec::with_capacity(registros.len());
    let mut secuencia = Vec::new();
    for registro in registros {
        if !secuencia.is_empty() {
            secuencia.push(b'N');
        }
        inicios.push(secuencia.len());
        longitudes.push(registro.secuencia.len());
        secuencia.extend(
            registro
                .secuencia
                .iter()
                .map(|base| base.to_ascii_uppercase()),
        );
        nombres.push(registro.nombre);
    }
    let pares: Vec<(u32, u32)> = codigos_de_adn(&secuencia)
        .into_iter()
        .enumerate()
        .filter(|&(_, codigo)| codigo != SIN_CODIGO)
        .map(|(posicion, codigo)| (codigo, posicion as u32))
        .collect();
    BaseDeAdn {
        total_letras: longitudes.iter().sum(),
        nombres,
        inicios,
        longitudes,
        secuencia,
        indice: Indice::nuevo(&pares, PALABRAS_DE_ADN),
    }
}

pub fn genoma_de_ecoli() -> BaseDeAdn {
    construir_base_de_adn(vec![Registro {
        nombre: ecoli::IDENTIFICADOR.to_string(),
        secuencia: ecoli::genoma().to_vec(),
    }])
}

#[derive(Debug, Clone, PartialEq)]
pub struct Acierto {
    pub sujeto: usize,
    pub marco: i8,
    pub inicio_consulta: usize,
    pub fin_consulta: usize,
    pub inicio_sujeto: usize,
    pub fin_sujeto: usize,
    pub puntaje: i32,
    pub bits: f64,
    pub evalor: f64,
    pub menos_log10_evalor: f64,
    pub longitud: usize,
    pub identidades: usize,
    pub positivos: usize,
    pub desapareamientos: usize,
    pub aperturas: usize,
    pub huecos: usize,
    pub consulta_alineada: Vec<u8>,
    pub linea_media: Vec<u8>,
    pub sujeto_alineado: Vec<u8>,
}

pub fn armar_acierto(
    consulta: &[u8],
    sujeto: &[u8],
    numero: usize,
    alineamiento: &AlineamientoLocal,
    parametros: &Parametros,
    espacio: f64,
) -> Acierto {
    let bits_del_acierto = bits(alineamiento.puntaje, &parametros.estadistica);
    let mut acierto = Acierto {
        sujeto: numero,
        marco: 0,
        inicio_consulta: alineamiento.inicio_consulta,
        fin_consulta: alineamiento.fin_consulta,
        inicio_sujeto: alineamiento.inicio_sujeto,
        fin_sujeto: alineamiento.fin_sujeto,
        puntaje: alineamiento.puntaje,
        bits: bits_del_acierto,
        evalor: evalor(bits_del_acierto, espacio),
        menos_log10_evalor: bits_del_acierto * std::f64::consts::LOG10_2 - espacio.log10(),
        longitud: alineamiento.operaciones.len(),
        identidades: 0,
        positivos: 0,
        desapareamientos: 0,
        aperturas: 0,
        huecos: 0,
        consulta_alineada: Vec::with_capacity(alineamiento.operaciones.len()),
        linea_media: Vec::with_capacity(alineamiento.operaciones.len()),
        sujeto_alineado: Vec::with_capacity(alineamiento.operaciones.len()),
    };
    let mut en_consulta = alineamiento.inicio_consulta;
    let mut en_sujeto = alineamiento.inicio_sujeto;
    let mut anterior = PAREJA;
    for &operacion in &alineamiento.operaciones {
        match operacion {
            HUECO_EN_SUJETO => {
                acierto.consulta_alineada.push(consulta[en_consulta]);
                acierto.sujeto_alineado.push(b'-');
                acierto.linea_media.push(b' ');
                en_consulta += 1;
            }
            HUECO_EN_CONSULTA => {
                acierto.consulta_alineada.push(b'-');
                acierto.sujeto_alineado.push(sujeto[en_sujeto]);
                acierto.linea_media.push(b' ');
                en_sujeto += 1;
            }
            _ => {
                let de_consulta = consulta[en_consulta];
                let de_sujeto = sujeto[en_sujeto];
                let puntaje = parametros.puntuacion.puntaje(de_consulta, de_sujeto);
                acierto.consulta_alineada.push(de_consulta);
                acierto.sujeto_alineado.push(de_sujeto);
                if de_consulta == de_sujeto && puntaje > 0 {
                    acierto.identidades += 1;
                    acierto.positivos += 1;
                    acierto
                        .linea_media
                        .push(if parametros.adn { b'|' } else { de_consulta });
                } else {
                    acierto.desapareamientos += 1;
                    if puntaje > 0 {
                        acierto.positivos += 1;
                        acierto.linea_media.push(b'+');
                    } else {
                        acierto.linea_media.push(b' ');
                    }
                }
                en_consulta += 1;
                en_sujeto += 1;
            }
        }
        if operacion != PAREJA {
            acierto.huecos += 1;
            if operacion != anterior {
                acierto.aperturas += 1;
            }
        }
        anterior = operacion;
    }
    acierto
}

pub fn punto_de_partida(
    consulta: &[u8],
    sujeto: &[u8],
    tramo: &Tramo,
    puntuacion: &Puntuacion,
) -> (usize, usize) {
    if tramo.longitud <= VENTANA_DE_SEMILLA {
        return (
            tramo.inicio_consulta + tramo.longitud / 2,
            tramo.inicio_sujeto + tramo.longitud / 2,
        );
    }
    let par = |desplazamiento: usize| {
        puntuacion.puntaje(
            consulta[tramo.inicio_consulta + desplazamiento],
            sujeto[tramo.inicio_sujeto + desplazamiento],
        )
    };
    let mut suma: i32 = (0..VENTANA_DE_SEMILLA).map(par).sum();
    let mut mejor = suma;
    let mut mejor_inicio = 0;
    for inicio in 1..=tramo.longitud - VENTANA_DE_SEMILLA {
        suma += par(inicio + VENTANA_DE_SEMILLA - 1) - par(inicio - 1);
        if suma > mejor {
            mejor = suma;
            mejor_inicio = inicio;
        }
    }
    (
        tramo.inicio_consulta + mejor_inicio + VENTANA_DE_SEMILLA / 2,
        tramo.inicio_sujeto + mejor_inicio + VENTANA_DE_SEMILLA / 2,
    )
}

fn contiene(alineamiento: &AlineamientoLocal, tramo: &Tramo) -> bool {
    tramo.inicio_consulta >= alineamiento.inicio_consulta
        && tramo.inicio_consulta + tramo.longitud <= alineamiento.fin_consulta
        && tramo.inicio_sujeto >= alineamiento.inicio_sujeto
        && tramo.inicio_sujeto + tramo.longitud <= alineamiento.fin_sujeto
}

fn comparten_extremo(uno: &AlineamientoLocal, otro: &AlineamientoLocal) -> bool {
    let mismo_inicio =
        uno.inicio_consulta == otro.inicio_consulta && uno.inicio_sujeto == otro.inicio_sujeto;
    let mismo_fin = uno.fin_consulta == otro.fin_consulta && uno.fin_sujeto == otro.fin_sujeto;
    mismo_inicio || mismo_fin
}

pub fn alinear_tramos(
    consulta: &[u8],
    sujeto: &[u8],
    mut tramos: Vec<Tramo>,
    parametros: &Parametros,
    minimo: i32,
) -> Vec<AlineamientoLocal> {
    tramos.sort_by(|uno, otro| {
        otro.puntaje
            .cmp(&uno.puntaje)
            .then(uno.inicio_consulta.cmp(&otro.inicio_consulta))
            .then(uno.inicio_sujeto.cmp(&otro.inicio_sujeto))
    });
    tramos.truncate(LIMITE_DE_TRAMOS);
    let mut alineamientos: Vec<AlineamientoLocal> = Vec::new();
    for tramo in &tramos {
        if alineamientos.iter().any(|hecho| contiene(hecho, tramo)) {
            continue;
        }
        let (semilla_consulta, semilla_sujeto) =
            punto_de_partida(consulta, sujeto, tramo, &parametros.puntuacion);
        let preliminar = alinear_con_huecos(
            consulta,
            sujeto,
            semilla_consulta,
            semilla_sujeto,
            parametros,
            parametros.caida_preliminar,
        );
        if preliminar.puntaje < minimo {
            continue;
        }
        let nuevo = alinear_con_huecos(
            consulta,
            sujeto,
            semilla_consulta,
            semilla_sujeto,
            parametros,
            parametros.caida_con_huecos,
        );
        if alineamientos
            .iter()
            .any(|hecho| comparten_extremo(hecho, &nuevo))
        {
            continue;
        }
        alineamientos.push(nuevo);
    }
    alineamientos
}

#[derive(Clone, Copy)]
struct Diagonal {
    sujeto: usize,
    ultimo_golpe: i64,
    fin_extension: usize,
}

const DIAGONAL_SIN_USAR: Diagonal = Diagonal {
    sujeto: usize::MAX,
    ultimo_golpe: NUNCA,
    fin_extension: 0,
};

fn tramos_de_dos_golpes(
    consulta: &[u8],
    sujeto: &[u8],
    numero: usize,
    palabras_del_sujeto: &[u16],
    tabla: &Indice,
    parametros: &Parametros,
    diagonales: &mut Vec<Diagonal>,
) -> Vec<Tramo> {
    let mut tramos = Vec::new();
    let necesarias = consulta.len() + sujeto.len();
    if diagonales.len() < necesarias {
        diagonales.resize(necesarias, DIAGONAL_SIN_USAR);
    }
    let ventana = parametros.ventana_dos_golpes as i64;
    let palabra = parametros.longitud_palabra as i64;
    for (posicion_sujeto, &codigo) in palabras_del_sujeto.iter().enumerate() {
        if codigo == SIN_PALABRA {
            continue;
        }
        for &posicion_consulta in tabla.posiciones_de(codigo as u32) {
            let posicion_consulta = posicion_consulta as usize;
            let diagonal = &mut diagonales[posicion_sujeto + consulta.len() - posicion_consulta];
            if diagonal.sujeto != numero {
                *diagonal = Diagonal {
                    sujeto: numero,
                    ..DIAGONAL_SIN_USAR
                };
            }
            if posicion_sujeto < diagonal.fin_extension {
                continue;
            }
            let distancia = posicion_sujeto as i64 - diagonal.ultimo_golpe;
            if distancia >= ventana {
                diagonal.ultimo_golpe = posicion_sujeto as i64;
                continue;
            }
            if distancia < palabra {
                continue;
            }
            let tramo = extender_sin_huecos(
                consulta,
                sujeto,
                posicion_consulta,
                posicion_sujeto,
                parametros.longitud_palabra,
                &parametros.puntuacion,
                parametros.caida_sin_huecos,
            );
            diagonal.ultimo_golpe = posicion_sujeto as i64;
            diagonal.fin_extension = tramo.inicio_sujeto + tramo.longitud;
            if tramo.puntaje >= parametros.disparo_con_huecos {
                tramos.push(tramo);
            }
        }
    }
    tramos
}

fn demasiados_golpes(golpes: u64, limite: u64) -> String {
    format!(
        "La consulta produce {golpes} coincidencias de palabras con la base y el máximo es {limite}: suele pasar con secuencias repetitivas o de baja complejidad. Acorta la consulta o quita las repeticiones."
    )
}

pub fn buscar_proteinas(
    consulta: &[u8],
    base: &BaseDeProteinas,
    parametros: &Parametros,
    evalor_maximo: f64,
) -> Result<Vec<Acierto>, String> {
    if base.registros.is_empty() {
        return Err("La base de datos no tiene ninguna secuencia.".to_string());
    }
    let tabla = tabla_de_vecinas(consulta, parametros.umbral_vecinas);
    let mut golpes: u64 = 0;
    for (palabra, &veces) in base.veces_por_palabra[..PALABRAS_DE_PROTEINA]
        .iter()
        .enumerate()
    {
        let en_consulta = tabla.posiciones_de(palabra as u32).len() as u64;
        golpes += en_consulta * veces;
    }
    if golpes > LIMITE_DE_GOLPES {
        return Err(demasiados_golpes(golpes, LIMITE_DE_GOLPES));
    }
    let espacio = espacio_de_busqueda(
        consulta.len() as f64,
        base.total_letras as f64,
        base.registros.len() as f64,
        &parametros.estadistica,
    );
    let minimo = puntaje_necesario(evalor_maximo, espacio, &parametros.estadistica);
    let aciertos = (0..base.registros.len())
        .into_par_iter()
        .map_init(Vec::new, |diagonales, numero| {
            let sujeto = &base.registros[numero].secuencia;
            let tramos = tramos_de_dos_golpes(
                consulta,
                sujeto,
                numero,
                &base.palabras[numero],
                &tabla,
                parametros,
                diagonales,
            );
            alinear_tramos(consulta, sujeto, tramos, parametros, minimo)
                .iter()
                .map(|alineamiento| {
                    armar_acierto(consulta, sujeto, numero, alineamiento, parametros, espacio)
                })
                .filter(|acierto| acierto.evalor <= evalor_maximo)
                .collect::<Vec<Acierto>>()
        })
        .flatten()
        .collect();
    Ok(aciertos)
}

pub fn ordenar_y_recortar(aciertos: &mut Vec<Acierto>, maximo: usize) {
    aciertos.sort_by(|uno, otro| {
        otro.menos_log10_evalor
            .total_cmp(&uno.menos_log10_evalor)
            .then(otro.puntaje.cmp(&uno.puntaje))
            .then(uno.sujeto.cmp(&otro.sujeto))
            .then(uno.inicio_consulta.cmp(&otro.inicio_consulta))
            .then(uno.inicio_sujeto.cmp(&otro.inicio_sujeto))
    });
    aciertos.truncate(maximo);
}

fn preparar_consulta(
    consulta: &[u8],
    minimo: usize,
    maximo: usize,
    unidad: &str,
) -> Result<Vec<u8>, String> {
    if consulta.len() < minimo {
        return Err(format!(
            "La consulta debe tener al menos {minimo} {unidad}."
        ));
    }
    if consulta.len() > maximo {
        return Err(format!(
            "La consulta tiene {} {unidad} y el máximo es {maximo}.",
            consulta.len()
        ));
    }
    Ok(consulta.to_ascii_uppercase())
}

pub fn blastp(
    consulta: &[u8],
    base: &BaseDeProteinas,
    evalor_maximo: f64,
    maximo: usize,
) -> Result<Vec<Acierto>, String> {
    let consulta = preparar_consulta(
        consulta,
        LONGITUD_PALABRA_PROTEINA,
        LONGITUD_MAXIMA_PROTEINA,
        "aminoácidos",
    )?;
    let parametros = parametros_blastp();
    let mut aciertos = buscar_proteinas(&consulta, base, &parametros, evalor_maximo)?;
    ordenar_y_recortar(&mut aciertos, maximo);
    Ok(aciertos)
}

pub fn blastx(
    consulta: &[u8],
    base: &BaseDeProteinas,
    evalor_maximo: f64,
    maximo: usize,
) -> Result<Vec<Acierto>, String> {
    let consulta = preparar_consulta(
        consulta,
        3 * LONGITUD_PALABRA_PROTEINA,
        LONGITUD_MAXIMA_BLASTX,
        "nucleótidos",
    )?;
    let parametros = parametros_blastp();
    let inversa = complemento_inverso(&consulta);
    let mut aciertos = Vec::new();
    for marco in [1i8, 2, 3, -1, -2, -3] {
        let hebra = if marco > 0 { &consulta } else { &inversa };
        let desplazamiento = marco.unsigned_abs() as usize - 1;
        let proteina = traducir(&hebra[desplazamiento..]);
        if proteina.len() < LONGITUD_PALABRA_PROTEINA {
            continue;
        }
        for mut acierto in buscar_proteinas(&proteina, base, &parametros, evalor_maximo)? {
            let inicio = desplazamiento + 3 * acierto.inicio_consulta;
            let fin = desplazamiento + 3 * acierto.fin_consulta;
            if marco > 0 {
                acierto.inicio_consulta = inicio;
                acierto.fin_consulta = fin;
            } else {
                acierto.inicio_consulta = consulta.len() - fin;
                acierto.fin_consulta = consulta.len() - inicio;
            }
            acierto.marco = marco;
            aciertos.push(acierto);
        }
    }
    ordenar_y_recortar(&mut aciertos, maximo);
    Ok(aciertos)
}

pub fn golpes_de_adn(consulta: &[u8], base: &BaseDeAdn) -> u64 {
    codigos_de_adn(consulta)
        .into_iter()
        .filter(|&codigo| codigo != SIN_CODIGO)
        .map(|codigo| base.indice.posiciones_de(codigo).len() as u64)
        .sum()
}

pub fn tramos_de_semillas_exactas(
    consulta: &[u8],
    base: &BaseDeAdn,
    parametros: &Parametros,
) -> Vec<(usize, Tramo)> {
    let mut tramos = Vec::new();
    let mut fin_por_diagonal: HashMap<i64, usize> = HashMap::new();
    for (posicion_consulta, codigo) in codigos_de_adn(consulta).into_iter().enumerate() {
        if codigo == SIN_CODIGO {
            continue;
        }
        for &posicion in base.indice.posiciones_de(codigo) {
            let posicion = posicion as usize;
            let diagonal = posicion as i64 - posicion_consulta as i64;
            if fin_por_diagonal
                .get(&diagonal)
                .is_some_and(|&fin| posicion < fin)
            {
                continue;
            }
            let numero = base.registro_en(posicion);
            let inicio_del_registro = base.inicios[numero];
            let tramo = extender_sin_huecos(
                consulta,
                base.registro(numero),
                posicion_consulta,
                posicion - inicio_del_registro,
                parametros.longitud_palabra,
                &parametros.puntuacion,
                parametros.caida_sin_huecos,
            );
            fin_por_diagonal.insert(
                diagonal,
                inicio_del_registro + tramo.inicio_sujeto + tramo.longitud,
            );
            if tramo.puntaje >= parametros.disparo_con_huecos {
                tramos.push((numero, tramo));
            }
        }
    }
    tramos
}

pub fn voltear_a_la_hebra_directa(acierto: &mut Acierto, longitud_consulta: usize) {
    let inicio = longitud_consulta - acierto.fin_consulta;
    acierto.fin_consulta = longitud_consulta - acierto.inicio_consulta;
    acierto.inicio_consulta = inicio;
    acierto.consulta_alineada = complemento_inverso(&acierto.consulta_alineada);
    acierto.sujeto_alineado = complemento_inverso(&acierto.sujeto_alineado);
    acierto.linea_media.reverse();
}

pub fn blastn(
    consulta: &[u8],
    base: &BaseDeAdn,
    evalor_maximo: f64,
    maximo: usize,
) -> Result<Vec<Acierto>, String> {
    let consulta = preparar_consulta(
        consulta,
        LONGITUD_PALABRA_ADN,
        LONGITUD_MAXIMA_ADN,
        "nucleótidos",
    )?;
    let parametros = parametros_blastn();
    let inversa = complemento_inverso(&consulta);
    let golpes = golpes_de_adn(&consulta, base) + golpes_de_adn(&inversa, base);
    if golpes > LIMITE_DE_GOLPES_ADN {
        return Err(demasiados_golpes(golpes, LIMITE_DE_GOLPES_ADN));
    }
    let espacio = espacio_de_busqueda(
        consulta.len() as f64,
        base.total_letras as f64,
        base.cantidad() as f64,
        &parametros.estadistica,
    );
    let minimo = puntaje_necesario(evalor_maximo, espacio, &parametros.estadistica);
    let mut aciertos = Vec::new();
    for (hebra, secuencia) in [(1i8, &consulta), (-1i8, &inversa)] {
        let mut tramos = tramos_de_semillas_exactas(secuencia, base, &parametros);
        tramos.sort_by_key(|&(numero, _)| numero);
        for grupo in tramos.chunk_by(|uno, otro| uno.0 == otro.0) {
            let numero = grupo[0].0;
            let sujeto = base.registro(numero);
            let del_registro: Vec<Tramo> = grupo.iter().map(|&(_, tramo)| tramo).collect();
            for alineamiento in alinear_tramos(secuencia, sujeto, del_registro, &parametros, minimo)
            {
                let mut acierto = armar_acierto(
                    secuencia,
                    sujeto,
                    numero,
                    &alineamiento,
                    &parametros,
                    espacio,
                );
                if acierto.evalor > evalor_maximo {
                    continue;
                }
                acierto.marco = hebra;
                if hebra < 0 {
                    voltear_a_la_hebra_directa(&mut acierto, consulta.len());
                }
                aciertos.push(acierto);
            }
        }
    }
    ordenar_y_recortar(&mut aciertos, maximo);
    Ok(aciertos)
}
