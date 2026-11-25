use crate::secuencia::complemento_inverso;
use flate2::read::GzDecoder;
use std::io::Read;
use std::sync::OnceLock;

pub const NOMBRE: &str = "Escherichia coli K-12 MG1655 (NC_000913.3)";
pub const IDENTIFICADOR: &str = "NC_000913.3";

const GENOMA_EMPACADO: &[u8] = include_bytes!("../datos/ecoli.genoma");
const GENES_COMPRIMIDOS: &[u8] = include_bytes!("../datos/ecoli.genes.gz");

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Gen {
    pub locus: String,
    pub nombre: String,
    pub inicio: usize,
    pub fin: usize,
    pub hebra: u8,
    pub producto: String,
    pub proteina: Vec<u8>,
}

pub fn genoma() -> &'static [u8] {
    static GENOMA: OnceLock<Vec<u8>> = OnceLock::new();
    GENOMA.get_or_init(|| desempacar(GENOMA_EMPACADO))
}

pub fn desempacar(empacado: &[u8]) -> Vec<u8> {
    const BASES: &[u8; 4] = b"ACGT";
    let mut longitud_bytes = [0u8; 8];
    longitud_bytes.copy_from_slice(&empacado[..8]);
    let longitud = u64::from_le_bytes(longitud_bytes) as usize;
    let mut genoma = Vec::with_capacity(longitud);
    for &byte in &empacado[8..] {
        for desplazamiento in [0, 2, 4, 6] {
            if genoma.len() == longitud {
                break;
            }
            genoma.push(BASES[((byte >> desplazamiento) & 0b11) as usize]);
        }
    }
    genoma
}

pub fn genes() -> &'static [Gen] {
    static GENES: OnceLock<Vec<Gen>> = OnceLock::new();
    GENES.get_or_init(|| {
        let mut texto = String::new();
        GzDecoder::new(GENES_COMPRIMIDOS)
            .read_to_string(&mut texto)
            .expect("la tabla de genes incluida es válida");
        texto.lines().filter_map(leer_gen).collect()
    })
}

fn leer_gen(linea: &str) -> Option<Gen> {
    let columnas: Vec<&str> = linea.split('\t').collect();
    if columnas.len() != 7 {
        return None;
    }
    Some(Gen {
        locus: columnas[0].to_string(),
        nombre: columnas[1].to_string(),
        inicio: columnas[2].parse().ok()?,
        fin: columnas[3].parse().ok()?,
        hebra: columnas[4].bytes().next()?,
        producto: columnas[5].to_string(),
        proteina: columnas[6].as_bytes().to_vec(),
    })
}

pub fn buscar_gen(nombre: &str) -> Option<&'static Gen> {
    genes().iter().find(|candidato| {
        candidato.nombre.eq_ignore_ascii_case(nombre)
            || candidato.locus.eq_ignore_ascii_case(nombre)
    })
}

pub fn region(inicio: usize, fin: usize) -> Option<&'static [u8]> {
    let genoma = genoma();
    if inicio == 0 || inicio > fin || fin > genoma.len() {
        return None;
    }
    Some(&genoma[inicio - 1..fin])
}

pub fn secuencia_gen(elegido: &Gen) -> Vec<u8> {
    let tramo = &genoma()[elegido.inicio - 1..elegido.fin];
    if elegido.hebra == b'-' {
        complemento_inverso(tramo)
    } else {
        tramo.to_vec()
    }
}
