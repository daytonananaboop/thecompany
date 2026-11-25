use crate::formatos::{
    LectorFastq, Lectura, RegistroSam, abrir_fastq, escribir_cabecera_sam, escribir_sam,
};
use crate::secuencia::{complemento_inverso, indice_base};
use rayon::prelude::*;
use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

pub const K: usize = 15;
pub const W: usize = 10;
pub const FRECUENCIA_MAXIMA: usize = 50;
pub const LONGITUD_MAXIMA_DE_LECTURA: usize = 1000;
pub const LARGO_VENTANA: usize = 10_000;
pub const CALIDAD_MAXIMA: u8 = 60;

pub const MAXIMO_SALTO: usize = 500;
pub const BANDA_DE_CADENA: usize = 100;
const PREDECESORES_A_REVISAR: usize = 50;
pub const MINIMO_DE_ANCLAS: usize = 2;
pub const PUNTAJE_MINIMO_DE_CADENA: f64 = 20.0;
const MAXIMO_DE_CANDIDATAS: usize = 20;
pub const PUNTAJE_MINIMO_DE_ALINEAMIENTO: i32 = 30;

pub const COINCIDENCIA: i32 = 1;
pub const DIFERENCIA: i32 = 4;
pub const ABRIR_HUECO: i32 = 6;
pub const EXTENDER_HUECO: i32 = 1;
pub const PREMIO_POR_LLEGAR_AL_EXTREMO: i32 = 5;
pub const MARGEN_DE_BANDA: i64 = 10;
const MAXIMO_DE_CELDAS: usize = 4_000_000;
const MENOS_INFINITO: i32 = i32::MIN / 4;

pub const LECTURAS_POR_TANDA: usize = 16_384;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Minimizador {
    pub valor: u64,
    pub posicion: usize,
    pub reverso: bool,
}

pub struct Indice {
    pub k: usize,
    pub w: usize,
    pub semillas: Vec<Minimizador>,
    pub cajones: Vec<usize>,
    pub bits_de_cajon: u32,
    pub longitud_referencia: usize,
    pub frecuencia_maxima: usize,
    pub descartados: usize,
    pub apariciones_descartadas: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Ancla {
    pub reverso: bool,
    pub en_referencia: usize,
    pub en_lectura: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Cadena {
    pub reverso: bool,
    pub anclas: Vec<Ancla>,
    pub puntaje: f64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alineamiento {
    pub posicion: usize,
    pub cigar: Vec<(u32, u8)>,
    pub puntaje: i32,
    pub coincidencias: usize,
    pub diferencias: usize,
    pub insertadas: usize,
    pub borradas: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Mapeo {
    pub reverso: bool,
    pub alineamiento: Alineamiento,
    pub calidad_mapeo: u8,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Tramo {
    pub inicio: usize,
    pub fin: usize,
    pub profundidad_media: f64,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct ResumenMapeo {
    pub lecturas: usize,
    pub mapeadas: usize,
    pub en_hebra_reversa: usize,
    pub con_recorte: usize,
    pub demasiado_largas: usize,
    pub calidades: Vec<usize>,
    pub coincidencias: usize,
    pub diferencias: usize,
    pub insertadas: usize,
    pub borradas: usize,
    pub longitud_referencia: usize,
    pub largo_ventana: usize,
    pub profundidad_por_ventana: Vec<f64>,
    pub bases_cubiertas: usize,
    pub profundidad_media: f64,
    pub tramos: Vec<Tramo>,
    pub minimizadores_en_indice: usize,
    pub minimizadores_descartados: usize,
}

pub fn barajar(codigo: u64) -> u64 {
    let mut valor = codigo.wrapping_add(0x9E37_79B9_7F4A_7C15);
    valor = (valor ^ (valor >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    valor = (valor ^ (valor >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    valor ^ (valor >> 31)
}

pub fn minimizadores(secuencia: &[u8], k: usize, w: usize) -> Vec<Minimizador> {
    let mut elegidos = Vec::new();
    if k == 0 || k > 31 || w == 0 || secuencia.len() < k {
        return elegidos;
    }
    let mascara = (1u64 << (2 * k)) - 1;
    let lugar_de_la_primera = 2 * (k as u64 - 1);
    let mut directo = 0u64;
    let mut inverso = 0u64;
    let mut bases_seguidas = 0;
    let mut ventana: VecDeque<Minimizador> = VecDeque::new();
    for (indice, &base) in secuencia.iter().enumerate() {
        match indice_base(base) {
            Some(numero) => {
                let numero = numero as u64;
                directo = ((directo << 2) | numero) & mascara;
                inverso = (inverso >> 2) | ((3 - numero) << lugar_de_la_primera);
                bases_seguidas += 1;
            }
            None => bases_seguidas = 0,
        }
        if indice + 1 < k {
            continue;
        }
        let inicio = indice + 1 - k;
        if bases_seguidas >= k && directo != inverso {
            let candidato = Minimizador {
                valor: barajar(directo.min(inverso)),
                posicion: inicio,
                reverso: inverso < directo,
            };
            while ventana
                .back()
                .is_some_and(|ultimo| ultimo.valor > candidato.valor)
            {
                ventana.pop_back();
            }
            ventana.push_back(candidato);
        }
        while ventana
            .front()
            .is_some_and(|primero| primero.posicion + w <= inicio)
        {
            ventana.pop_front();
        }
        if inicio + 1 >= w
            && let Some(&menor) = ventana.front()
            && elegidos.last() != Some(&menor)
        {
            elegidos.push(menor);
        }
    }
    elegidos
}

pub fn construir_indice(
    referencia: &[u8],
    k: usize,
    w: usize,
    frecuencia_maxima: usize,
) -> Result<Indice, String> {
    if k == 0 || k > 31 {
        return Err("El largo k de los minimizadores debe estar entre 1 y 31.".to_string());
    }
    if w == 0 {
        return Err("La ventana w debe abarcar al menos un k-mero.".to_string());
    }
    if frecuencia_maxima == 0 {
        return Err("La frecuencia máxima de un minimizador debe ser al menos 1.".to_string());
    }
    if referencia.len() < k + w - 1 {
        return Err(format!(
            "La referencia debe medir al menos {} bases.",
            k + w - 1
        ));
    }
    let mut semillas = minimizadores(referencia, k, w);
    if semillas.is_empty() {
        return Err(
            "La referencia no tiene ningún tramo de A, C, G y T donde buscar minimizadores."
                .to_string(),
        );
    }
    semillas.par_sort_unstable_by_key(|semilla| (semilla.valor, semilla.posicion));
    let mut indice = Indice {
        k,
        w,
        semillas: Vec::with_capacity(semillas.len()),
        cajones: Vec::new(),
        bits_de_cajon: 0,
        longitud_referencia: referencia.len(),
        frecuencia_maxima,
        descartados: 0,
        apariciones_descartadas: 0,
    };
    guardar_las_poco_frecuentes(&mut indice, &semillas);
    indice.bits_de_cajon = (usize::BITS - indice.semillas.len().leading_zeros()).clamp(1, 24);
    indice.cajones = repartir_en_cajones(&indice.semillas, indice.bits_de_cajon);
    Ok(indice)
}

fn guardar_las_poco_frecuentes(indice: &mut Indice, ordenadas: &[Minimizador]) {
    let mut inicio = 0;
    while inicio < ordenadas.len() {
        let mut fin = inicio + 1;
        while fin < ordenadas.len() && ordenadas[fin].valor == ordenadas[inicio].valor {
            fin += 1;
        }
        let apariciones = fin - inicio;
        if apariciones <= indice.frecuencia_maxima {
            indice.semillas.extend_from_slice(&ordenadas[inicio..fin]);
        } else {
            indice.descartados += 1;
            indice.apariciones_descartadas += apariciones;
        }
        inicio = fin;
    }
}

fn cajon_de(valor: u64, bits_de_cajon: u32) -> usize {
    (valor >> (64 - bits_de_cajon)) as usize
}

fn repartir_en_cajones(semillas: &[Minimizador], bits_de_cajon: u32) -> Vec<usize> {
    let cantidad = 1usize << bits_de_cajon;
    let mut cajones = Vec::with_capacity(cantidad + 1);
    let mut siguiente = 0;
    for cajon in 0..=cantidad {
        while siguiente < semillas.len()
            && cajon_de(semillas[siguiente].valor, bits_de_cajon) < cajon
        {
            siguiente += 1;
        }
        cajones.push(siguiente);
    }
    cajones
}

impl Indice {
    pub fn buscar(&self, valor: u64) -> &[Minimizador] {
        let cajon = cajon_de(valor, self.bits_de_cajon);
        let del_cajon = &self.semillas[self.cajones[cajon]..self.cajones[cajon + 1]];
        let inicio = del_cajon.partition_point(|semilla| semilla.valor < valor);
        let mut fin = inicio;
        while fin < del_cajon.len() && del_cajon[fin].valor == valor {
            fin += 1;
        }
        &del_cajon[inicio..fin]
    }
}

pub fn buscar_anclas(lectura: &[u8], indice: &Indice) -> Vec<Ancla> {
    let mut anclas = Vec::new();
    for minimizador in minimizadores(lectura, indice.k, indice.w) {
        for semilla in indice.buscar(minimizador.valor) {
            let reverso = minimizador.reverso != semilla.reverso;
            let en_lectura = if reverso {
                lectura.len() - indice.k - minimizador.posicion
            } else {
                minimizador.posicion
            };
            anclas.push(Ancla {
                reverso,
                en_referencia: semilla.posicion,
                en_lectura,
            });
        }
    }
    anclas
}

pub fn encadenar(anclas: &[Ancla], k: usize) -> Vec<Cadena> {
    let mut anclas = anclas.to_vec();
    anclas.sort_unstable_by_key(|ancla| (ancla.reverso, ancla.en_referencia, ancla.en_lectura));
    let cantidad = anclas.len();
    let mut puntaje = vec![0.0; cantidad];
    let mut anterior: Vec<Option<usize>> = vec![None; cantidad];
    for actual in 0..cantidad {
        puntaje[actual] = k as f64;
        let desde = actual.saturating_sub(PREDECESORES_A_REVISAR);
        for previa in (desde..actual).rev() {
            if anclas[previa].reverso != anclas[actual].reverso
                || anclas[actual].en_referencia - anclas[previa].en_referencia > MAXIMO_SALTO
            {
                break;
            }
            if let Some(ganancia) = ganancia_de_enlace(anclas[previa], anclas[actual], k) {
                let propuesta = puntaje[previa] + ganancia;
                if propuesta > puntaje[actual] {
                    puntaje[actual] = propuesta;
                    anterior[actual] = Some(previa);
                }
            }
        }
    }
    separar_cadenas(&anclas, &puntaje, &anterior)
}

pub fn ganancia_de_enlace(previa: Ancla, actual: Ancla, k: usize) -> Option<f64> {
    if actual.en_referencia <= previa.en_referencia || actual.en_lectura <= previa.en_lectura {
        return None;
    }
    let avance_en_referencia = actual.en_referencia - previa.en_referencia;
    let avance_en_lectura = actual.en_lectura - previa.en_lectura;
    let desvio = avance_en_referencia.abs_diff(avance_en_lectura);
    if avance_en_referencia > MAXIMO_SALTO
        || avance_en_lectura > MAXIMO_SALTO
        || desvio > BANDA_DE_CADENA
    {
        return None;
    }
    let bases_nuevas = avance_en_referencia.min(avance_en_lectura).min(k);
    Some(bases_nuevas as f64 - costo_de_desvio(desvio, k))
}

fn costo_de_desvio(desvio: usize, k: usize) -> f64 {
    if desvio == 0 {
        return 0.0;
    }
    0.01 * k as f64 * desvio as f64 + 0.5 * (desvio as f64).log2()
}

fn separar_cadenas(anclas: &[Ancla], puntaje: &[f64], anterior: &[Option<usize>]) -> Vec<Cadena> {
    let mut finales: Vec<usize> = (0..anclas.len()).collect();
    finales.sort_by(|&una, &otra| puntaje[otra].total_cmp(&puntaje[una]).then(una.cmp(&otra)));
    let mut usada = vec![false; anclas.len()];
    let mut cadenas = Vec::new();
    for fin in finales {
        if usada[fin] {
            continue;
        }
        let mut miembros = Vec::new();
        let mut ya_contado = 0.0;
        let mut actual = Some(fin);
        while let Some(indice) = actual {
            if usada[indice] {
                ya_contado = puntaje[indice];
                break;
            }
            usada[indice] = true;
            miembros.push(anclas[indice]);
            actual = anterior[indice];
        }
        miembros.reverse();
        cadenas.push(Cadena {
            reverso: anclas[fin].reverso,
            anclas: miembros,
            puntaje: puntaje[fin] - ya_contado,
        });
    }
    cadenas.sort_by(|una, otra| otra.puntaje.total_cmp(&una.puntaje));
    cadenas
}

impl Cadena {
    pub fn diagonales(&self) -> (i64, i64) {
        let mut menor = i64::MAX;
        let mut mayor = i64::MIN;
        for ancla in &self.anclas {
            let diagonal = ancla.en_referencia as i64 - ancla.en_lectura as i64;
            menor = menor.min(diagonal);
            mayor = mayor.max(diagonal);
        }
        (menor, mayor)
    }

    fn se_superpone_con(&self, otra: &Cadena) -> bool {
        let (Some(primera), Some(ultima)) = (self.anclas.first(), self.anclas.last()) else {
            return false;
        };
        let (Some(primera_otra), Some(ultima_otra)) = (otra.anclas.first(), otra.anclas.last())
        else {
            return false;
        };
        let (menor, mayor) = self.diagonales();
        let (menor_otra, mayor_otra) = otra.diagonales();
        self.reverso == otra.reverso
            && primera.en_referencia <= ultima_otra.en_referencia
            && primera_otra.en_referencia <= ultima.en_referencia
            && menor <= mayor_otra + MARGEN_DE_BANDA
            && menor_otra <= mayor + MARGEN_DE_BANDA
    }
}

pub fn elegir_candidatas(cadenas: Vec<Cadena>) -> Vec<Cadena> {
    let mut candidatas: Vec<Cadena> = Vec::new();
    for cadena in cadenas {
        if cadena.anclas.len() < MINIMO_DE_ANCLAS || cadena.puntaje < PUNTAJE_MINIMO_DE_CADENA {
            continue;
        }
        if let Some(mejor) = candidatas.first()
            && cadena.puntaje < mejor.puntaje / 2.0
        {
            break;
        }
        if candidatas.iter().any(|otra| otra.se_superpone_con(&cadena)) {
            continue;
        }
        candidatas.push(cadena);
        if candidatas.len() == MAXIMO_DE_CANDIDATAS {
            break;
        }
    }
    candidatas
}

#[derive(Clone, Copy, PartialEq)]
enum Origen {
    Recorte,
    Diagonal,
    Borrado,
    Insercion,
}

enum Capa {
    Principal,
    Borrado,
    Insercion,
}

#[derive(Clone, Copy)]
struct Puntajes {
    mejor: i32,
    borrado: i32,
    insercion: i32,
}

const PUNTAJES_VACIOS: Puntajes = Puntajes {
    mejor: MENOS_INFINITO,
    borrado: MENOS_INFINITO,
    insercion: MENOS_INFINITO,
};

#[derive(Clone, Copy)]
struct Rastro {
    origen: Origen,
    borrado_abierto: bool,
    insercion_abierta: bool,
}

const RASTRO_VACIO: Rastro = Rastro {
    origen: Origen::Recorte,
    borrado_abierto: false,
    insercion_abierta: false,
};

struct Tabla {
    ancho: usize,
    primera_diagonal: i64,
    rastros: Vec<Rastro>,
    mejor_final: Option<(usize, usize, i32)>,
}

impl Tabla {
    fn fin_en_referencia(&self, fila: usize, columna: usize) -> i64 {
        fila as i64 + self.primera_diagonal + columna as i64
    }
}

fn puntaje_de_par(base_lectura: u8, base_referencia: u8) -> i32 {
    if base_lectura == base_referencia {
        COINCIDENCIA
    } else {
        -DIFERENCIA
    }
}

pub fn alinear_con_recorte(
    lectura: &[u8],
    referencia: &[u8],
    diagonal_menor: i64,
    diagonal_mayor: i64,
) -> Option<Alineamiento> {
    if lectura.is_empty() || referencia.is_empty() || diagonal_menor > diagonal_mayor {
        return None;
    }
    let ancho = diagonal_mayor
        .saturating_sub(diagonal_menor)
        .saturating_add(2 * MARGEN_DE_BANDA + 1);
    let ancho = usize::try_from(ancho).ok()?;
    let filas = lectura.len() + 1;
    if filas.saturating_mul(ancho) > MAXIMO_DE_CELDAS
        || diagonal_menor > referencia.len() as i64
        || diagonal_mayor < -(filas as i64)
    {
        return None;
    }
    if diagonal_menor == diagonal_mayor
        && let Some(directo) = alinear_sin_huecos(lectura, referencia, diagonal_menor)
    {
        return Some(directo);
    }
    let mut tabla = Tabla {
        ancho,
        primera_diagonal: diagonal_menor - MARGEN_DE_BANDA,
        rastros: vec![RASTRO_VACIO; filas * ancho],
        mejor_final: None,
    };
    llenar_tabla(&mut tabla, lectura, referencia);
    reconstruir(&tabla, lectura, referencia)
}

fn alinear_sin_huecos(lectura: &[u8], referencia: &[u8], diagonal: i64) -> Option<Alineamiento> {
    let inicio = usize::try_from(diagonal).ok()?;
    let tramo = referencia.get(inicio..inicio + lectura.len())?;
    let mut diferencias = 0;
    for (base_lectura, base_referencia) in lectura.iter().zip(tramo) {
        if base_lectura != base_referencia {
            diferencias += 1;
            if diferencias > 1 {
                return None;
            }
        }
    }
    let coincidencias = lectura.len() - diferencias;
    let puntaje = coincidencias as i32 * COINCIDENCIA - diferencias as i32 * DIFERENCIA;
    if puntaje <= 0 {
        return None;
    }
    Some(Alineamiento {
        posicion: inicio,
        cigar: vec![(lectura.len() as u32, b'M')],
        puntaje,
        coincidencias,
        diferencias,
        insertadas: 0,
        borradas: 0,
    })
}

fn llenar_tabla(tabla: &mut Tabla, lectura: &[u8], referencia: &[u8]) {
    let ancho = tabla.ancho;
    let largo_referencia = referencia.len() as i64;
    let mut anterior = vec![PUNTAJES_VACIOS; ancho];
    let mut actual = vec![PUNTAJES_VACIOS; ancho];
    for (columna, puntajes) in anterior.iter_mut().enumerate() {
        let fin = tabla.fin_en_referencia(0, columna);
        if fin >= 0 && fin <= largo_referencia {
            puntajes.mejor = PREMIO_POR_LLEGAR_AL_EXTREMO;
        }
    }
    let mut mejor_total = 0;
    for fila in 1..=lectura.len() {
        let base_lectura = lectura[fila - 1];
        let premio_final = if fila == lectura.len() {
            PREMIO_POR_LLEGAR_AL_EXTREMO
        } else {
            0
        };
        for columna in 0..ancho {
            actual[columna] = PUNTAJES_VACIOS;
            let fin = tabla.fin_en_referencia(fila, columna);
            if fin < 0 || fin > largo_referencia {
                continue;
            }
            let mut puntajes = PUNTAJES_VACIOS;
            let mut rastro = RASTRO_VACIO;
            if columna > 0 {
                let izquierda = actual[columna - 1];
                let abrir = izquierda.mejor - ABRIR_HUECO - EXTENDER_HUECO;
                let seguir = izquierda.borrado - EXTENDER_HUECO;
                puntajes.borrado = abrir.max(seguir);
                rastro.borrado_abierto = abrir >= seguir;
            }
            if columna + 1 < ancho {
                let arriba = anterior[columna + 1];
                let abrir = arriba.mejor - ABRIR_HUECO - EXTENDER_HUECO;
                let seguir = arriba.insercion - EXTENDER_HUECO;
                puntajes.insercion = abrir.max(seguir);
                rastro.insercion_abierta = abrir >= seguir;
            }
            puntajes.mejor = 0;
            if fin >= 1 {
                let base_referencia = referencia[(fin - 1) as usize];
                let diagonal =
                    anterior[columna].mejor + puntaje_de_par(base_lectura, base_referencia);
                if diagonal >= puntajes.mejor {
                    puntajes.mejor = diagonal;
                    rastro.origen = Origen::Diagonal;
                }
            }
            if puntajes.borrado > puntajes.mejor {
                puntajes.mejor = puntajes.borrado;
                rastro.origen = Origen::Borrado;
            }
            if puntajes.insercion > puntajes.mejor {
                puntajes.mejor = puntajes.insercion;
                rastro.origen = Origen::Insercion;
            }
            actual[columna] = puntajes;
            tabla.rastros[fila * ancho + columna] = rastro;
            if rastro.origen != Origen::Recorte && puntajes.mejor + premio_final > mejor_total {
                mejor_total = puntajes.mejor + premio_final;
                tabla.mejor_final = Some((fila, columna, puntajes.mejor));
            }
        }
        std::mem::swap(&mut anterior, &mut actual);
    }
}

fn reconstruir(tabla: &Tabla, lectura: &[u8], referencia: &[u8]) -> Option<Alineamiento> {
    let (fila_final, columna_final, mut puntaje) = tabla.mejor_final?;
    let mut operaciones = vec![b'S'; lectura.len() - fila_final];
    let mut fila = fila_final;
    let mut columna = columna_final;
    let mut capa = Capa::Principal;
    loop {
        let rastro = tabla.rastros[fila * tabla.ancho + columna];
        match capa {
            Capa::Principal => {
                if fila == 0 {
                    puntaje -= PREMIO_POR_LLEGAR_AL_EXTREMO;
                    break;
                }
                match rastro.origen {
                    Origen::Recorte => {
                        operaciones.extend(std::iter::repeat_n(b'S', fila));
                        break;
                    }
                    Origen::Diagonal => {
                        operaciones.push(b'M');
                        fila -= 1;
                    }
                    Origen::Borrado => capa = Capa::Borrado,
                    Origen::Insercion => capa = Capa::Insercion,
                }
            }
            Capa::Borrado => {
                operaciones.push(b'D');
                if rastro.borrado_abierto {
                    capa = Capa::Principal;
                }
                columna -= 1;
            }
            Capa::Insercion => {
                operaciones.push(b'I');
                if rastro.insercion_abierta {
                    capa = Capa::Principal;
                }
                fila -= 1;
                columna += 1;
            }
        }
    }
    operaciones.reverse();
    let inicio = tabla.fin_en_referencia(fila, columna) as usize;
    let alineamiento = describir(&operaciones, lectura, referencia, inicio, puntaje);
    if alineamiento.coincidencias == 0 || alineamiento.puntaje <= 0 {
        return None;
    }
    Some(alineamiento)
}

fn describir(
    operaciones: &[u8],
    lectura: &[u8],
    referencia: &[u8],
    inicio: usize,
    puntaje: i32,
) -> Alineamiento {
    let mut alineamiento = Alineamiento {
        posicion: inicio,
        cigar: comprimir(operaciones),
        puntaje,
        coincidencias: 0,
        diferencias: 0,
        insertadas: 0,
        borradas: 0,
    };
    let mut en_lectura = 0;
    let mut en_referencia = inicio;
    for &operacion in operaciones {
        match operacion {
            b'M' => {
                if lectura[en_lectura] == referencia[en_referencia] {
                    alineamiento.coincidencias += 1;
                } else {
                    alineamiento.diferencias += 1;
                }
                en_lectura += 1;
                en_referencia += 1;
            }
            b'I' => {
                alineamiento.insertadas += 1;
                en_lectura += 1;
            }
            b'D' => {
                alineamiento.borradas += 1;
                en_referencia += 1;
            }
            _ => en_lectura += 1,
        }
    }
    alineamiento
}

fn comprimir(operaciones: &[u8]) -> Vec<(u32, u8)> {
    let mut cigar: Vec<(u32, u8)> = Vec::new();
    for &operacion in operaciones {
        match cigar.last_mut() {
            Some((cantidad, ultima)) if *ultima == operacion => *cantidad += 1,
            _ => cigar.push((1, operacion)),
        }
    }
    cigar
}

impl Alineamiento {
    pub fn largo_en_referencia(&self) -> usize {
        self.coincidencias + self.diferencias + self.borradas
    }

    pub fn tiene_recorte(&self) -> bool {
        self.cigar.iter().any(|&(_, operacion)| operacion == b'S')
    }

    fn recorte_inicial(&self) -> usize {
        match self.cigar.first() {
            Some(&(cantidad, b'S')) => cantidad as usize,
            _ => 0,
        }
    }

    fn diagonal(&self) -> i64 {
        self.posicion as i64 - self.recorte_inicial() as i64
    }

    pub fn se_superpone_con(&self, otro: &Alineamiento) -> bool {
        self.posicion < otro.posicion + otro.largo_en_referencia()
            && otro.posicion < self.posicion + self.largo_en_referencia()
            && (self.diagonal() - otro.diagonal()).abs() <= MARGEN_DE_BANDA
    }
}

pub fn calidad_de_mapeo(mejor: i32, segundo: i32, anclas: usize) -> u8 {
    if mejor <= 0 {
        return 0;
    }
    let segundo = segundo.clamp(0, mejor);
    let unicidad = 1.0 - segundo as f64 / mejor as f64;
    let apoyo = (anclas as f64 / 10.0).min(1.0);
    let calidad = 40.0 * unicidad * apoyo * (mejor as f64).ln();
    calidad.clamp(0.0, CALIDAD_MAXIMA as f64) as u8
}

struct Candidato {
    reverso: bool,
    alineamiento: Alineamiento,
    anclas: usize,
}

fn agregar_sin_repetir(candidatos: &mut Vec<Candidato>, nuevo: Candidato) {
    for existente in candidatos.iter_mut() {
        if existente.reverso == nuevo.reverso
            && existente.alineamiento.se_superpone_con(&nuevo.alineamiento)
        {
            if nuevo.alineamiento.puntaje > existente.alineamiento.puntaje {
                *existente = nuevo;
            }
            return;
        }
    }
    candidatos.push(nuevo);
}

pub fn mapear_lectura(
    lectura: &[u8],
    indice: &Indice,
    referencia: &[u8],
    desempate: u64,
) -> Option<Mapeo> {
    if lectura.len() > LONGITUD_MAXIMA_DE_LECTURA {
        return None;
    }
    let anclas = buscar_anclas(lectura, indice);
    let cadenas = elegir_candidatas(encadenar(&anclas, indice.k));
    if cadenas.is_empty() {
        return None;
    }
    let inversa = complemento_inverso(lectura);
    let mut candidatos: Vec<Candidato> = Vec::new();
    for cadena in &cadenas {
        let orientada = if cadena.reverso { &inversa } else { lectura };
        let (menor, mayor) = cadena.diagonales();
        let Some(alineamiento) = alinear_con_recorte(orientada, referencia, menor, mayor) else {
            continue;
        };
        if alineamiento.puntaje >= PUNTAJE_MINIMO_DE_ALINEAMIENTO {
            let nuevo = Candidato {
                reverso: cadena.reverso,
                alineamiento,
                anclas: cadena.anclas.len(),
            };
            agregar_sin_repetir(&mut candidatos, nuevo);
        }
    }
    let mejor = candidatos
        .iter()
        .map(|candidato| candidato.alineamiento.puntaje)
        .max()?;
    let empatados: Vec<usize> = (0..candidatos.len())
        .filter(|&numero| candidatos[numero].alineamiento.puntaje == mejor)
        .collect();
    let elegido = empatados[(desempate % empatados.len() as u64) as usize];
    let mut segundo = 0;
    for (numero, candidato) in candidatos.iter().enumerate() {
        if numero != elegido {
            segundo = segundo.max(candidato.alineamiento.puntaje);
        }
    }
    let ganador = candidatos.swap_remove(elegido);
    Some(Mapeo {
        reverso: ganador.reverso,
        calidad_mapeo: calidad_de_mapeo(mejor, segundo, ganador.anclas),
        alineamiento: ganador.alineamiento,
    })
}

pub fn mapear_lecturas(
    ruta_lecturas: &Path,
    referencia: &[u8],
    nombre_referencia: &str,
    ruta_sam: &Path,
) -> Result<ResumenMapeo, String> {
    let indice = construir_indice(referencia, K, W, FRECUENCIA_MAXIMA)?;
    mapear_fastq(
        ruta_lecturas,
        &indice,
        referencia,
        nombre_referencia,
        ruta_sam,
    )
}

pub fn desempate(nombre: &[u8]) -> u64 {
    let mut valor = 0u64;
    for &letra in nombre {
        valor = barajar(valor ^ letra as u64);
    }
    valor
}

pub fn mapear_fastq(
    ruta_lecturas: &Path,
    indice: &Indice,
    referencia: &[u8],
    nombre_referencia: &str,
    ruta_sam: &Path,
) -> Result<ResumenMapeo, String> {
    if indice.longitud_referencia != referencia.len() {
        return Err("El índice no corresponde a esta referencia: miden distinto.".to_string());
    }
    if es_el_mismo_archivo(ruta_lecturas, ruta_sam) {
        return Err("El SAM de salida no puede ser el mismo archivo que las lecturas.".to_string());
    }
    let mut lector = abrir_fastq(ruta_lecturas)?;
    let archivo = File::create(ruta_sam)
        .map_err(|error| format!("No se pudo crear {}: {error}", ruta_sam.display()))?;
    let error_de_escritura =
        |error: io::Error| format!("No se pudo escribir {}: {error}", ruta_sam.display());
    let mut destino = BufWriter::with_capacity(1 << 20, archivo);
    escribir_cabecera_sam(&mut destino, nombre_referencia, referencia.len())
        .map_err(error_de_escritura)?;
    let mut resumen = ResumenMapeo {
        calidades: vec![0; CALIDAD_MAXIMA as usize + 1],
        longitud_referencia: referencia.len(),
        largo_ventana: largo_de_ventana(referencia.len()),
        minimizadores_en_indice: indice.semillas.len(),
        minimizadores_descartados: indice.descartados,
        ..ResumenMapeo::default()
    };
    let mut cambios_de_profundidad = vec![0i64; referencia.len() + 1];
    loop {
        let tanda = leer_tanda(&mut lector)?;
        if tanda.is_empty() {
            break;
        }
        let mapeos = mapear_en_paralelo(&tanda, indice, referencia);
        for (lectura, mapeo) in tanda.iter().zip(&mapeos) {
            anotar(
                &mut resumen,
                &mut cambios_de_profundidad,
                lectura,
                mapeo.as_ref(),
            );
            let registro = registro_sam(lectura, mapeo.as_ref(), nombre_referencia);
            escribir_sam(&mut destino, &registro).map_err(error_de_escritura)?;
        }
    }
    destino.flush().map_err(error_de_escritura)?;
    medir_profundidad(&mut resumen, &cambios_de_profundidad);
    Ok(resumen)
}

fn es_el_mismo_archivo(una: &Path, otra: &Path) -> bool {
    match (una.canonicalize(), otra.canonicalize()) {
        (Ok(primera), Ok(segunda)) => primera == segunda,
        _ => false,
    }
}

pub fn largo_de_ventana(longitud_referencia: usize) -> usize {
    if longitud_referencia >= 100 * LARGO_VENTANA {
        LARGO_VENTANA
    } else {
        (longitud_referencia / 100).max(1)
    }
}

fn leer_tanda(lector: &mut LectorFastq) -> Result<Vec<Lectura>, String> {
    let mut tanda = Vec::with_capacity(LECTURAS_POR_TANDA);
    while tanda.len() < LECTURAS_POR_TANDA {
        match lector.next() {
            Some(lectura) => tanda.push(lectura?),
            None => break,
        }
    }
    Ok(tanda)
}

pub fn mapear_en_paralelo(
    tanda: &[Lectura],
    indice: &Indice,
    referencia: &[u8],
) -> Vec<Option<Mapeo>> {
    tanda
        .par_iter()
        .map(|lectura| {
            mapear_lectura(
                &lectura.secuencia,
                indice,
                referencia,
                desempate(&lectura.nombre),
            )
        })
        .collect()
}

pub fn anotar(
    resumen: &mut ResumenMapeo,
    cambios_de_profundidad: &mut [i64],
    lectura: &Lectura,
    mapeo: Option<&Mapeo>,
) {
    resumen.lecturas += 1;
    if lectura.secuencia.len() > LONGITUD_MAXIMA_DE_LECTURA {
        resumen.demasiado_largas += 1;
    }
    let Some(mapeo) = mapeo else {
        return;
    };
    let alineamiento = &mapeo.alineamiento;
    resumen.mapeadas += 1;
    if mapeo.reverso {
        resumen.en_hebra_reversa += 1;
    }
    if alineamiento.tiene_recorte() {
        resumen.con_recorte += 1;
    }
    resumen.calidades[mapeo.calidad_mapeo as usize] += 1;
    resumen.coincidencias += alineamiento.coincidencias;
    resumen.diferencias += alineamiento.diferencias;
    resumen.insertadas += alineamiento.insertadas;
    resumen.borradas += alineamiento.borradas;
    cambios_de_profundidad[alineamiento.posicion] += 1;
    cambios_de_profundidad[alineamiento.posicion + alineamiento.largo_en_referencia()] -= 1;
}

fn registro_sam(lectura: &Lectura, mapeo: Option<&Mapeo>, nombre_referencia: &str) -> RegistroSam {
    let nombre = nombre_sin_comentario(&lectura.nombre);
    let Some(mapeo) = mapeo else {
        return RegistroSam {
            nombre,
            bandera: 4,
            referencia: String::new(),
            posicion: 0,
            calidad_mapeo: 0,
            cigar: Vec::new(),
            secuencia: lectura.secuencia.clone(),
            calidad: lectura.calidad.clone(),
        };
    };
    let (bandera, secuencia, calidad) = if mapeo.reverso {
        (
            16,
            complemento_inverso(&lectura.secuencia),
            lectura.calidad.iter().rev().copied().collect(),
        )
    } else {
        (0, lectura.secuencia.clone(), lectura.calidad.clone())
    };
    RegistroSam {
        nombre,
        bandera,
        referencia: nombre_referencia.to_string(),
        posicion: mapeo.alineamiento.posicion + 1,
        calidad_mapeo: mapeo.calidad_mapeo,
        cigar: mapeo.alineamiento.cigar.clone(),
        secuencia,
        calidad,
    }
}

fn nombre_sin_comentario(nombre: &[u8]) -> Vec<u8> {
    let corto: Vec<u8> = nombre
        .iter()
        .copied()
        .take_while(|letra| !letra.is_ascii_whitespace())
        .collect();
    if corto.is_empty() {
        b"*".to_vec()
    } else {
        corto
    }
}

pub fn medir_profundidad(resumen: &mut ResumenMapeo, cambios_de_profundidad: &[i64]) {
    let profundidad = profundidad_por_base(cambios_de_profundidad, resumen.longitud_referencia);
    let mut suma: u64 = 0;
    for &valor in &profundidad {
        if valor > 0 {
            resumen.bases_cubiertas += 1;
            suma += valor as u64;
        }
    }
    if resumen.bases_cubiertas > 0 {
        resumen.profundidad_media = suma as f64 / resumen.bases_cubiertas as f64;
    }
    resumen.profundidad_por_ventana = promedios_por_ventana(&profundidad, resumen.largo_ventana);
    resumen.tramos = tramos_cubiertos(&profundidad);
}

fn profundidad_por_base(cambios: &[i64], longitud: usize) -> Vec<u32> {
    let mut profundidad = Vec::with_capacity(longitud);
    let mut actual: i64 = 0;
    for &cambio in &cambios[..longitud] {
        actual += cambio;
        profundidad.push(actual as u32);
    }
    profundidad
}

fn promedios_por_ventana(profundidad: &[u32], largo_ventana: usize) -> Vec<f64> {
    profundidad
        .chunks(largo_ventana)
        .map(|ventana| {
            let suma: u64 = ventana.iter().map(|&valor| valor as u64).sum();
            suma as f64 / ventana.len() as f64
        })
        .collect()
}

fn tramos_cubiertos(profundidad: &[u32]) -> Vec<Tramo> {
    let mut tramos = Vec::new();
    let mut inicio = 0;
    while inicio < profundidad.len() {
        if profundidad[inicio] == 0 {
            inicio += 1;
            continue;
        }
        let mut fin = inicio;
        let mut suma: u64 = 0;
        while fin < profundidad.len() && profundidad[fin] > 0 {
            suma += profundidad[fin] as u64;
            fin += 1;
        }
        tramos.push(Tramo {
            inicio,
            fin,
            profundidad_media: suma as f64 / (fin - inicio) as f64,
        });
        inicio = fin;
    }
    tramos
}
