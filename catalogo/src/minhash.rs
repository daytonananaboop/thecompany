use crate::formato::{decimal, miles, porcentaje};
use crate::pizarra::{Marca, Narracion, Pizarra};
use crate::{Algoritmo, Bloque, Campo, Entrada, Paso, Salida, Tipo};
use nucleo::filogenia;
use nucleo::formatos;
use nucleo::minhash::{self, Bosquejo, Comparacion, Menores, Par, Ventana};
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::HashSet;
use std::io::{BufRead, Read};
use std::path::Path;
use std::time::{Duration, Instant};

const EJEMPLO: &str = "@ecoli:355001-375000
@ecoli:365001-385000
@gen:lacZ
@ecoli:2000001-2020000
>lacZ_1
ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTT
TCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTTTGCCTGGTTTCCGGCAGC
AGAAGCGGTGCCGGAAAGCTGGCTGGAGTGCGATCTTCCTGAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCACGGTTACCATGCGCCCATC
TACACCAACGTGACCTATCCCATTACGGTCAATCCGCCGTTTGTTCCCACGGAGAATCCGACGTGTTGTTACTCGCTCACATTTAATGTTGATGAAAGCT
GGCTACAGGAAGGCCAGACGCGAATTATTTTTGATGGCGTTAACTCGGCGTTTCATCTGTGGTGCAACGGGCGCTGGGTCGGTTACGGCCAGGACAGTCG
TTTGCCGTCTGAATTTGACCTGAGCGCATTTTTACGCGCCGGAGAAAACCGCCTCGCGGTGATGGTGCTGCGCTGGAGTGACGGCAGTTATCTGGAAGAT
CAGGATATGTGGCGGATGAGCGGCATTTTCCGTGACGTCTCGTTGCTGCATAAACCGACGACACAAATCAGCGATTTCCATGTTGCCACTCGCTTTAATG
ATGATTTCAGCCGCGCTGTACTGGAGGCTGAAGTTCAGATGTGCGGCGAGTTGCGTGACTACCTACGGGTAACAGTTTCTTTATAGCAGGGCGAAACACA
GGTCGCCAGCGGCACCGCGCCTTTCGGCGGTGAAATTATCGATGAGCGTGGTGGTTATGCCGATCGCGTGACACTACGGCTGAACGTCGAAAACCCGAAA
CTGTGGAGTGCCGAAATCCCGAATCTCTATCGTGCGGTGGTTGAACTGCACACCGCCGACGGCACGCTGATTGAAGCACAAGCCTGCGATGTCGGTTTCC
GCGAGGTGCGGATTGAAAATGGTCTGCTGCTGCTGAACGGCAAGCAGTTGCTGATTCGAGGCGTTAACCGTCACGAGCATCATCCTCTGCATGGTCAGGT
CATGGATGAGCAGACGATGGTGCAGGATATCCTGCTGATGAAGCAGAACAACTTTAACGCCGTGCGCTGTTCGCATTATCCCAACCATCCGCTGTGGTAC
ACGCTGTGCGACCGCTACGGCCTGTATGTGGTGGATGAAGCCAATATTGAAACCCACGGCATGGTGCCAATGAATCGTCTGACCGATGATCCGCGCTGGC
TACCGGCGATGAGCGAACGCGTAACGCGAATGGTGAAGCGCGATCGTGATCACCCGAGTGTGATCATCTGGTCGCTGGGGAATGAATCAGGCCACGGCGC
TAATCACGACGCGCTGTATCGCTGGATCAAATCTGTCGATCCTTCCCGCCCGGTGCAGTATGAAGGCGGGGGAGCCGACACCACGGCCACCGATATTATT
TCCCCGATGTACGCGCGCGTGGATGAAGACCAGCCCTTCCCGGCTGTGCCGAAATGGTCCATCAAAAAATGGCTTTCGCTACCTGGAGAGACGCGCCCGC
TGATCCTTTGCGAATACGCCCACGCGATGGGTAACAGTCTTGGCGGTTTCGCTAGATACTGGCAGGCGTTTCGTCAGTATCCCCGTTTACAGGGCGGCTT
CGTCTGGGACTGGGTGGATCAGTCGCTGATTAAATATGATGAAAACGGCAACCCGTGGTCGGCTTACGGCGGTGATTTTGGCGATACGCCGAAAGATCGC
CAGTTCTGTATGAACGGTCTGGTGTTTGCCGACCGCACGCCGCATCCAGCGCTGACGGAAGCAAAACACCAGCAGCAGTTTTTCCAGTTCCGTTTATCCG
GGCAAACCATCGAGGTGACCAGCGAATACCTGGTCCGTCATAGCGATAACGAGCTCCTGCACTGGATGGTGGCGCTGGATGGTAAGCCTCTGGCAAGCGG
TGAAGTGCCTCTGGATGTCGCTCCACAAGGTAAACAGTTGATTGAACTGCCTGAACTACCGCCGCCGGAGAGCGCCGGGCAACTCTGGCTCACAGTACGC
GTAGTGCAACCGAACGCGACCGCATGGTCAGAAGCCGGGCACATCAGCGCCTGGCAGCAGTGGCGTCTGGCGGAAAACCTCAGTGTGACGCTCCCCGCCG
AGTCCCACACCATCCCGCATCTGACCACCAGCGAAATGGATTTTTGCATCGAGCTGGGTAATAAGCGTTGGCAATTTAACCGCCAGTCAGGCTTTCTTTC
ACAGATGTGGATTGGCGATAAAAAACAACTGCTGACGCCGCTGCGCGATCAGTTCACCCGTGCACCGCTGGATAACGACACTGGCGTAAGTGAAGCGACC
CGCATTGACCCTAACGCCTGGGTCGAACGCTGGAAGGCGGCGGGCCATTACCAGGCCGAAGCAGCTTTGTTGCAGTGCACGGCAGATACACTTGCTGATG
CGGTGCTGATTACGACCGCTCACGCGTGGCAGCATCAGGGGAAAACCTTATTTATCAGCCGGAAAACCTACCGGATTGATGGTAGTGGTCAAATGGCGAT
TACCGTTGATATTGAAGTGGCGAGCGATACACCGCATCCGGCGCGGATTGGCCTGAACTGCCAGCTGGCGCAGGTAGCAGAGCGGGTAAACTGGCTCGGA
TTAGGGCCGCAAGAAAAATATCCCGACCGCCTTACTGCCGCCTGTTTTGACCGCTGGGATCTGCCATTGTCAGACATGTATACCCCGTACGTCTTCCCGA
GCGAAAACGGTCTGCGCTGCGGGACGCGCGAATTGAGTTATGGCCCACACCAGTGGCGCGGCGACTTCCAGTTCAACATCAGCCGCTACAGTCAACAGCA
ACTGATGGAAACCAGCCATCGCCATCTGCTGCACGCGGAAGAAGGCACATGGCTGAATATCGACGGTTTCCATATGGGGATTGGTGGCGACGACTCCTGG
AGCCCGTCAGTATCGGCGGAATTCCAGCTGAGCGCCGGTCGCTACCATTACCAGTTGGTCTGGTGTCAAAAATAA
>lacZ_5
ATTACCATGATTACGGATTCAATGGCCGTCGTCTTACAACGTAGCGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCATTGCGGCACATCCCCCTT
TCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGTTTTGCCTGGTTTCCGGCACC
AGAAGCGGTGCCGGAAATCTGGCTGGAGTGGGATCTTCCCGAGGCCGGTACTGTCGTCGTCCCATCAAACTGGCAGATGCATGGTTACGATGAGCACATC
TACACCAACGTGACCTATCCCATTACGGTCAACCCGTCGTTTGTTCCCACGGAGAATCCGACGGGTTGTTACTCGCTCACATTTAGTGTTGATGAAAGCT
GGCTACCGGAAGGCCAGACGCGAATTATTTTTGATGGCGTTAACTCGGCGTTTCATCTTTTGTGCAACGGGCGCTGGGTCGGCTACCGCCAGGACAGCCC
TTTGCCGTCTGAATTTGACCTGAGCGCATTTTTACGCCCCGGATAAAACCGCCTCGCGGTGATGGTGCTGCGCTGCAGTGACGGCAGTCATCTGGAAGAT
CACGATATGTGGTGGATGAGCGGCATTTTCCGTGACGTCTCGTGGCTGCATAAACCGACTACACAAATCAGCGATTTCCATGTTGCCACTCGCTTTAATC
ATGATTTCAGCCGCGTTGTACTGGAGGCTGAAGTTCAGGTTTGCGGCGAGTTGCGTGAGTACTTACGGGTAACAGTTCCTTTATGGCAGGGTGAAACGCA
GGTCGCCAGCGGCACCGCGCCTTTCGGCGGTGAAATTATCGATGAGCGTGGCGGTTATGCCGATCGCGTCACCCTACGTCTGAACGTCGAAAACCCGAAA
CTGTGGAGCGCCGAAATCCCGAATCTCTATCGTGCCGTGGTTGAACTGCACACCGCCGACGGCACGCTGATTGAAGCAGAAGCCTGCGATGTCGGTTTCC
GCGAGGCGCGTATTGAAAATCGTCAGCTGCTGCTGAACGGCAAGCCGTTGCTCATTCGAGGCGTTAACCGTGATGAGCATCATCCTCTGCATGGTCAGGT
CATGGATGATCAGACGATGGTGCAGGATATCCGGCTGATGAAGAAGAAAAACTTTAACCCCGTGCGCTGTTCGCATTATCCGAACCATCCGCGGTGGTAC
ACGATGTGCGACGGCTACGGCCTGTACATGGTGGATGAAGCCAATATTGAAACCCACGGCATGGTGCCAATGAATCGACTGACCGATGATCCTCGCTGGC
TACCGGCGATGAGCGAACGCGTAACGCGAATGGTGCAGCGCGATCGTAAGCACGCAAGTGTGATCGTCTGGTCACTGGGGAATTCATCACGCCACGGCGC
TAATCACGACGCGCTGTATCGCTGGATCAAATCTCTCGATCCTTCCCGCCCGGTGCAGTATGAAGGCGGCGGAGCCGACACCACGGCCACCGATATTAAT
TGCCCGATGTACGCGCGCGTGGATGAATACCAGCCCTACCCGACTGTGCCGAAATGGTCCATCACAGAATGGCTTTCGCTACCTGGCGAGACGCGCCCCC
TGATCCATTTCGAATACGCCCACGCGATGGGTAACAGTCTTGGCGGTTACGCTAAACACTGGCAGGGGTTTCGTCAGTATCCCCGTTTACAGGGCGGCTT
CGTCTGGGACTGGGTAGATCAGTCGCTGATTAAATATGATGAAAACGGCAACCGGTGGTCGGCTTACCGCGGTGATTTTGGCGATACGCCGAACGATCGC
CAGCTCTGTATGATCGGTCTGGTCTTTGCCGACGGCACGCCGCATCCAGCGCTGACGGAAGCAAAACACCTACAGCAGTTTTTCCAGTTCCGTTTATCCG
GGCAAACCATCGAAGTGACCAGCGAACACCTGTTCCGTCATAGCGATAACGCGCTCCTGCACTGGATGGTGGCACTGGATGGCAAGCCGCGGGCAAGCGG
TGAAGTGCCTCTGGATGTCGCTCCACAAGGTAAACAATTGATGCAACTGCCTGAACTACCGCAGCGGGAGATCGCCGGGCAAGTCTGGCTCACAGTACGC
GTAGTGCAACCGAACGCCGCCGCATGGTCAGAAGCCGGGAACATCAGCGCCTGACAGCAGTGGCGTCTGGCGGAAAACCTCAGAGGAACGCTCACCGCCG
CGTCCCACGCCATCCCGCGTCTGACAACCAACGATATGGATTTTTGCATTGAGCTGGGTCATAAGCATTGGCCATTTAACCGCCAGACAGGCTTTCTTTC
ACAGATGTGGATTGGCGATAAAAAACAACTGCTGGCGCCGCTGCGCGATCAGTTCACCCGTGCACCGCTGGATAACGACTTTGGCGTAAGTGAAGCGACC
CGCATTGACCCTAACGCCTGGGTCGAACGCTTGAAGTCGGCGGGCCATTACCAGGCCGAAGCAGCGTTGTTGCAGTGCACGGCAGATACACTTGCTGGTG
CGCTGCTGATTACGACCGCTCACGCGTGGCAGCATCAGCGGAAAACCTTATTTATCACCCGGAAAACCTACTGGATTGATGGTAGTGGTCAAATGGCGAT
TACCGTTGATGTTGCCGTGACGAGCGTTACACCGCATCCGTCGAGGATTGGCCTGAAGTTCCAGCTGGCGCAGGTAGCAGAGCGGGTAACCTGGCTCGGA
TTAGGGCCGCAAGAAAACTATCCCGACCGCCTTACTGCCGCCTGGTTTGACCGCTGGGATCTGCCATTGACAGACTTGTATACCCCGTACGTCTTCCCGA
GTGAAAACGGTCTGCGCTGCGGGACGCGCGAATTGAATTATGTCCCACACCAGTGGCGCGGCGACTTCCGGTTAAACATCAGCCGCGACAGTCAACAGCA
ACTGATGGAAACCAGCCATCGCCATCTGCTGCACGCGGAAGAAGGCACATGGCTGACTATAGTCAGTTTCCATATGGGGATTGGTGGCGACGGCTCCTGG
AGCCCGTTAGTATCGGCGGAATTCCAGCTGAGCGCCGGTAGCTACCATTACCAGTTGTTCTGGTGTCAAAAATAA";

const LETRAS_MAXIMAS: usize = 300_000_000;
const BYTES_MAXIMOS_DE_ARCHIVO: u64 = 300_000_000;
const BYTES_MAXIMOS_DESCOMPRIMIDOS: u64 = 330_000_000;
const LETRAS_PARA_EXACTO: usize = 12_000_000;
const FRONTERA_DE_ESPECIE: f64 = 0.95;
const ANI_CONFIABLE: f64 = 0.90;
const ALFA: f64 = 0.05;
const CASI_ENTERA: f64 = 0.9;
const TRAMO_MINIMO: f64 = 0.2;
const AZAR_TOLERADO: f64 = 0.01;
const FILAS_EN_PIZARRA: usize = 5;
const KMEROS_DETALLADOS: usize = 6;
const RECHAZOS_DETALLADOS: usize = 3;
const SALIDAS_DETALLADAS: usize = 6;
const LETRAS_RARAS_DETALLADAS: usize = 2;
const REPETIDOS_DETALLADOS: usize = 2;
const VUELTAS_DETALLADAS: usize = 8;
const COMPARTIDOS_DETALLADOS: usize = 2;
const PARES_EN_PIZARRA: usize = 12;

const L_COMPARAR: &str = "pub fn comparar(";
const L_VALIDAR: &str = "validar(genomas.len(), k, cupo)?;";
const L_MENORES: &str = "let mut menores = Menores::nuevo(cupo);";
const L_CODIGO: &str = "let codigo = CODIGO_DE_BASE[base as usize];";
const L_NO_ES_BASE: &str = "if codigo == NO_ES_BASE {";
const L_ADELANTE: &str = "self.adelante = ((self.adelante << 2) | codigo) & self.mascara;";
const L_ATRAS: &str = "self.atras = (self.atras >> 2) | ((3 - codigo) << self.corrimiento);";
const L_INCOMPLETO: &str = "if self.validas < self.k {";
const L_CANONICO: &str = "Some(mezclar(self.adelante.min(self.atras)))";
const L_MEZCLAR: &str = "let mut x = kmero ^ SEMILLA;";
const L_OFRECER: &str = "menores.ofrecer(hash);";
const L_RECHAZO: &str = "if lleno && hash >= self.mayor_guardado {";
const L_INSERTAR: &str = "let es_nuevo = self.conjunto.insert(hash);";
const L_SACAR: &str = "self.conjunto.pop_last();";
const L_UMBRAL: &str = "self.mayor_guardado = mayor;";
const L_FIN_BOSQUEJO: &str = "hashes: menores.conjunto.into_iter().collect(),";
const L_BOSQUEJAR_TODOS: &str = ".map(|piezas| bosquejar(piezas, k, cupo))";
const L_PARES: &str = "let pares = medir_todos_los_pares(&bosquejos);";
const L_CONTAR: &str = "contar_compartidos(&a.hashes, &b.hashes";
const L_MIENTRAS: &str = "while vistos < cupo && i < a.len() && j < b.len() {";
const L_IGUAL: &str = "compartidos += 1;";
const L_MENOR: &str = "Ordering::Less => i += 1,";
const L_MAYOR: &str = "Ordering::Greater => j += 1,";
const L_VISTOS: &str = "vistos += 1;";
const L_SOBRANTES: &str = "vistos += sobrantes.min(cupo - vistos);";
const L_JACCARD: &str = "let jaccard = fraccion(compartidos, vistos);";
const L_DISTANCIA: &str = "let distancia = -(2.0 * jaccard / (1.0 + jaccard)).ln() / k as f64;";
const L_IDENTIDAD: &str = "identidad: 1.0 - distancia,";
const L_PROBABILIDAD: &str = "1.0 / (1.0 + espacio / longitud as f64)";
const L_AZAR: &str = "let azar = jaccard_al_azar(a.longitud, b.longitud, a.k);";
const L_VALOR_P: &str = "log10_valor_p: log10_valor_p(compartidos, vistos, azar),";
const L_MEDIR: &str = "medida: medir(&bosquejos[primero], &bosquejos[segundo]),";
const L_CONTENCION: &str = "let contencion = contenciones(&bosquejos, genomas, k);";
const L_TODOS: &str = "let todos = hashes_de_todos(bosquejos);";
const L_PRESENTE: &str = "presentes[posicion] = true;";
const L_FRACCION_PRESENTE: &str = ".map(|en_grande| fraccion_presente(chico, en_grande))";
const L_FIN: &str = "Ok(Comparacion {";

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![Algoritmo {
        id: "puerto.minhash",
        capitulo: "puerto",
        nombre: "MinHash: comparar genomas sin alinear (Mash)",
        resumen: "Resume cada genoma en un bosquejo de 1 000 hashes y con eso estima el Jaccard, la distancia de Mash, el ANI y la contención en milisegundos; a diferencia de alinear o de las distancias p, nunca compara base por base.",
        campos: &[
            Campo {
                id: "secuencias",
                etiqueta: "Secuencias: una por línea (@ecoli:INICIO-FIN, @gen:NOMBRE) o en FASTA",
                tipo: Tipo::Secuencias,
                ejemplo: EJEMPLO,
            },
            Campo {
                id: "archivo",
                etiqueta: "Genomas en archivos FASTA (opcional; varios separados por ;)",
                tipo: Tipo::Archivo,
                ejemplo: "",
            },
            Campo {
                id: "k",
                etiqueta: "Largo de k-mero (k)",
                tipo: Tipo::Entero,
                ejemplo: "21",
            },
            Campo {
                id: "cupo",
                etiqueta: "Hashes por bosquejo",
                tipo: Tipo::Entero,
                ejemplo: "1000",
            },
        ],
        funcion: "minhash::comparar",
        ejecutar: ejecutar_minhash,
        pasos: Some(pasos_minhash),
    }]
}

struct Pedido {
    nombres: Vec<String>,
    genomas: Vec<Vec<Vec<u8>>>,
    k: usize,
    cupo: usize,
}

impl Pedido {
    fn bases(&self, indice: usize) -> usize {
        self.genomas[indice].iter().map(Vec::len).sum()
    }
}

fn es_separador(letra: char) -> bool {
    letra.is_whitespace() || letra == ',' || letra == ';'
}

fn solo_letras(texto: &str) -> Vec<u8> {
    texto
        .bytes()
        .filter(u8::is_ascii_alphabetic)
        .map(|letra| letra.to_ascii_uppercase())
        .collect()
}

fn nombre_de_referencia(pieza: &str) -> String {
    let sin_arroba = pieza.trim_start_matches('@');
    match sin_arroba.strip_prefix("gen:") {
        Some(nombre) => nombre.trim().to_string(),
        None => sin_arroba.to_string(),
    }
}

fn nombre_de_cabecera(cabecera: &str, cuantos: usize) -> String {
    cabecera
        .split_whitespace()
        .next()
        .map(str::to_string)
        .unwrap_or_else(|| format!("secuencia {}", cuantos + 1))
}

fn agregar_pieza(
    pieza: &str,
    nombres: &mut Vec<String>,
    genomas: &mut Vec<Vec<Vec<u8>>>,
) -> Result<(), String> {
    if pieza.starts_with('@') {
        let secuencia = Entrada::nueva(&[("secuencias", pieza)]).secuencia("secuencias")?;
        nombres.push(nombre_de_referencia(pieza));
        genomas.push(vec![secuencia]);
        return Ok(());
    }
    let secuencia = solo_letras(pieza);
    if secuencia.is_empty() {
        return Err(format!(
            "«{pieza}» no es una secuencia: escribe letras A, C, G y T, una referencia como @gen:lacZ o un FASTA."
        ));
    }
    nombres.push(format!("s{}", nombres.len() + 1));
    genomas.push(vec![secuencia]);
    Ok(())
}

fn leer_secuencias(
    texto: &str,
    nombres: &mut Vec<String>,
    genomas: &mut Vec<Vec<Vec<u8>>>,
) -> Result<(), String> {
    let mut en_fasta = false;
    for linea in texto.lines() {
        let linea = linea.trim();
        if linea.is_empty() {
            continue;
        }
        if let Some(cabecera) = linea.strip_prefix('>') {
            nombres.push(nombre_de_cabecera(cabecera, nombres.len()));
            genomas.push(vec![Vec::new()]);
            en_fasta = true;
            continue;
        }
        if en_fasta
            && !linea.starts_with('@')
            && let Some(ultimo) = genomas.last_mut()
        {
            ultimo[0].extend(solo_letras(linea));
            continue;
        }
        en_fasta = false;
        for pieza in linea.split(es_separador).filter(|pieza| !pieza.is_empty()) {
            agregar_pieza(pieza, nombres, genomas)?;
        }
    }
    Ok(())
}

fn nombre_de_archivo(ruta: &Path) -> String {
    let mut nombre = ruta
        .file_name()
        .map(|nombre| nombre.to_string_lossy().to_string())
        .unwrap_or_else(|| "archivo".to_string());
    for terminacion in [".gz", ".fasta", ".fa", ".fna", ".fas", ".ffn"] {
        if let Some(recortado) = nombre.strip_suffix(terminacion) {
            nombre = recortado.to_string();
        }
    }
    nombre
}

fn revisar_archivo(ruta: &Path) -> Result<(), String> {
    let datos = std::fs::metadata(ruta)
        .map_err(|_| format!("No existe el archivo «{}».", ruta.display()))?;
    if !datos.is_file() {
        return Err(format!("«{}» no es un archivo.", ruta.display()));
    }
    if datos.len() > BYTES_MAXIMOS_DE_ARCHIVO {
        return Err(format!(
            "El archivo «{}» pesa {} bytes y el máximo es {}.",
            ruta.display(),
            miles(datos.len() as usize),
            miles(BYTES_MAXIMOS_DE_ARCHIVO as usize)
        ));
    }
    let mut lector = formatos::abrir(ruta)?;
    let inicio = lector
        .fill_buf()
        .map_err(|error| format!("No se pudo leer {}: {error}", ruta.display()))?;
    if inicio.iter().find(|letra| !letra.is_ascii_whitespace()) != Some(&b'>') {
        return Err(format!(
            "«{}» no empieza con «>»: aquí van archivos FASTA de ADN (no FASTQ ni otros formatos).",
            ruta.display()
        ));
    }
    Ok(())
}

fn leer_piezas_de_archivo(ruta: &Path, maximo: u64) -> Result<Vec<Vec<u8>>, String> {
    let mut contenido = Vec::new();
    formatos::abrir(ruta)?
        .take(maximo + 1)
        .read_to_end(&mut contenido)
        .map_err(|error| format!("No se pudo leer {}: {error}", ruta.display()))?;
    if contenido.len() as u64 > maximo {
        return Err(format!(
            "«{}» ocupa más de {} bytes al descomprimirse y ese es el máximo.",
            ruta.display(),
            miles(maximo as usize)
        ));
    }
    Ok(formatos::leer_fasta(&contenido)
        .into_iter()
        .map(|registro| registro.secuencia)
        .filter(|secuencia| !secuencia.is_empty())
        .collect())
}

fn letras_en_total(genomas: &[Vec<Vec<u8>>]) -> usize {
    genomas.iter().flatten().map(Vec::len).sum()
}

fn revisar_letras(genomas: &[Vec<Vec<u8>>]) -> Result<(), String> {
    let letras = letras_en_total(genomas);
    if letras > LETRAS_MAXIMAS {
        return Err(format!(
            "Entre todas las secuencias hay {} letras y el máximo es {}.",
            miles(letras),
            miles(LETRAS_MAXIMAS)
        ));
    }
    Ok(())
}

fn leer_archivos(
    texto: &str,
    nombres: &mut Vec<String>,
    genomas: &mut Vec<Vec<Vec<u8>>>,
) -> Result<(), String> {
    for pieza in texto.split([';', '\n']) {
        let pieza = pieza.trim();
        if pieza.is_empty() {
            continue;
        }
        let ruta = Path::new(pieza);
        revisar_archivo(ruta)?;
        let piezas = leer_piezas_de_archivo(ruta, BYTES_MAXIMOS_DESCOMPRIMIDOS)?;
        if piezas.is_empty() {
            return Err(format!("«{}» no tiene secuencias.", ruta.display()));
        }
        nombres.push(nombre_de_archivo(ruta));
        genomas.push(piezas);
        revisar_letras(genomas)?;
    }
    Ok(())
}

fn hacer_unicos(nombres: &mut [String]) {
    let mut usados: HashSet<String> = HashSet::new();
    for nombre in nombres.iter_mut() {
        let base = nombre.clone();
        let mut numero = 2;
        while usados.contains(nombre.as_str()) {
            *nombre = format!("{base}_{numero}");
            numero += 1;
        }
        usados.insert(nombre.clone());
    }
}

fn tiene_un_kmero(piezas: &[Vec<u8>], k: usize) -> bool {
    piezas.iter().any(|pieza| {
        let mut seguidas = 0;
        for &letra in pieza {
            if matches!(letra, b'A' | b'C' | b'G' | b'T' | b'U') {
                seguidas += 1;
                if seguidas >= k {
                    return true;
                }
            } else {
                seguidas = 0;
            }
        }
        false
    })
}

fn leer_pedido(entrada: &Entrada) -> Result<Pedido, String> {
    let k = entrada.natural("k")?;
    let cupo = entrada.natural("cupo")?;
    let mut nombres = Vec::new();
    let mut genomas = Vec::new();
    if entrada.tiene("secuencias") {
        leer_secuencias(&entrada.texto("secuencias")?, &mut nombres, &mut genomas)?;
    }
    if entrada.tiene("archivo") {
        leer_archivos(&entrada.texto("archivo")?, &mut nombres, &mut genomas)?;
    }
    if genomas.is_empty() {
        return Err("Escribe al menos dos secuencias o elige archivos FASTA.".to_string());
    }
    minhash::validar(genomas.len(), k, cupo)?;
    revisar_letras(&genomas)?;
    for (nombre, piezas) in nombres.iter().zip(&genomas) {
        if !tiene_un_kmero(piezas, k) {
            return Err(format!(
                "«{nombre}» no tiene ningún tramo de {k} letras A, C, G o T seguidas, así que no tiene k-meros."
            ));
        }
    }
    hacer_unicos(&mut nombres);
    Ok(Pedido {
        nombres,
        genomas,
        k,
        cupo,
    })
}

struct Exacto {
    distintos: Vec<usize>,
    jaccard: Vec<Vec<f64>>,
    tiempo: Duration,
}

fn calcular_exacto(pedido: &Pedido) -> Option<Exacto> {
    if letras_en_total(&pedido.genomas) > LETRAS_PARA_EXACTO {
        return None;
    }
    let inicio = Instant::now();
    let conjuntos: Vec<Vec<u64>> = pedido
        .genomas
        .par_iter()
        .map(|piezas| minhash::todos_los_kmeros(piezas, pedido.k))
        .collect();
    let genomas_en_total = conjuntos.len();
    let mut jaccard = vec![vec![1.0; genomas_en_total]; genomas_en_total];
    for primero in 0..genomas_en_total {
        for segundo in primero + 1..genomas_en_total {
            let valor = minhash::jaccard_exacto(&conjuntos[primero], &conjuntos[segundo]);
            jaccard[primero][segundo] = valor;
            jaccard[segundo][primero] = valor;
        }
    }
    Some(Exacto {
        distintos: conjuntos.iter().map(Vec::len).collect(),
        jaccard,
        tiempo: inicio.elapsed(),
    })
}

fn texto_de_valor_p(log10: f64) -> String {
    if log10 >= 0.0 {
        return "1".to_string();
    }
    if log10 == f64::NEG_INFINITY {
        return "0".to_string();
    }
    if log10 > -3.0 {
        return decimal(10f64.powf(log10), 4);
    }
    let mut exponente = log10.floor();
    let mut mantisa = 10f64.powf(log10 - exponente);
    if mantisa >= 9.95 {
        mantisa = 1.0;
        exponente += 1.0;
    }
    format!("{}e{}", decimal(mantisa, 1), exponente as i64)
}

fn cantidad(numero: usize, singular: &str, plural: &str) -> String {
    let palabra = if numero == 1 { singular } else { plural };
    format!("{} {palabra}", miles(numero))
}

fn texto_de_tiempo(tiempo: Duration) -> String {
    let milisegundos = tiempo.as_secs_f64() * 1000.0;
    if milisegundos >= 1000.0 {
        format!("{} s", decimal(milisegundos / 1000.0, 2))
    } else {
        format!("{} ms", decimal(milisegundos, 1))
    }
}

fn hash_en_texto(hash: u64, cifras: usize) -> String {
    format!("{:.*}", cifras, minhash::en_fraccion(hash))
}

fn cifras_para(comparacion: &Comparacion) -> usize {
    let menor = comparacion
        .bosquejos
        .iter()
        .filter_map(|bosquejo| bosquejo.hashes.last())
        .map(|&mayor| minhash::en_fraccion(mayor))
        .fold(1.0, f64::min)
        .max(1e-12);
    ((-menor.log10()).ceil() as usize + 3).clamp(4, 12)
}

fn ani_en_texto(par: &Par) -> String {
    if par.medida.compartidos == 0 {
        "—".to_string()
    } else {
        porcentaje(par.medida.identidad, 2)
    }
}

fn pares_ordenados(comparacion: &Comparacion) -> Vec<Par> {
    let mut pares = comparacion.pares.clone();
    pares.sort_by(|uno, otro| {
        uno.medida
            .distancia
            .total_cmp(&otro.medida.distancia)
            .then(otro.medida.compartidos.cmp(&uno.medida.compartidos))
    });
    pares
}

fn tabla_de_pares(pedido: &Pedido, comparacion: &Comparacion, exacto: Option<&Exacto>) -> Bloque {
    let mut columnas = vec!["A", "B", "Hashes compartidos", "Jaccard estimado"];
    if exacto.is_some() {
        columnas.push("Jaccard exacto");
    }
    columnas.extend([
        "Distancia de Mash",
        "ANI estimado",
        "Valor p",
        "A dentro de B",
        "B dentro de A",
    ]);
    let filas = pares_ordenados(comparacion)
        .iter()
        .map(|par| {
            let (a, b) = (par.primero, par.segundo);
            let mut fila = vec![
                pedido.nombres[a].clone(),
                pedido.nombres[b].clone(),
                format!(
                    "{} de {}",
                    miles(par.medida.compartidos),
                    miles(par.medida.vistos)
                ),
                decimal(par.medida.jaccard, 4),
            ];
            if let Some(exacto) = exacto {
                fila.push(decimal(exacto.jaccard[a][b], 4));
            }
            fila.extend([
                decimal(par.medida.distancia, 4),
                ani_en_texto(par),
                texto_de_valor_p(par.medida.log10_valor_p),
                decimal(comparacion.contencion[a][b], 3),
                decimal(comparacion.contencion[b][a], 3),
            ]);
            fila
        })
        .collect();
    Bloque::tabla("Pares, del más parecido al menos", &columnas, filas)
}

fn tabla_de_bosquejos(
    pedido: &Pedido,
    comparacion: &Comparacion,
    exacto: Option<&Exacto>,
) -> Bloque {
    let mut columnas = vec![
        "Secuencia",
        "Bases",
        "k-meros leídos",
        "Hashes guardados",
        "Mayor hash guardado (de 0 a 1)",
        "k-meros distintos (estimado)",
    ];
    if exacto.is_some() {
        columnas.push("k-meros distintos (exacto)");
    }
    let filas = comparacion
        .bosquejos
        .iter()
        .enumerate()
        .map(|(indice, bosquejo)| {
            let mayor = bosquejo
                .hashes
                .last()
                .map(|&hash| decimal(minhash::en_fraccion(hash), 6))
                .unwrap_or_default();
            let mut fila = vec![
                pedido.nombres[indice].clone(),
                miles(bosquejo.longitud),
                miles(bosquejo.kmeros),
                miles(bosquejo.hashes.len()),
                mayor,
                miles(minhash::estimar_kmeros_distintos(bosquejo).round() as usize),
            ];
            if let Some(exacto) = exacto {
                fila.push(miles(exacto.distintos[indice]));
            }
            fila
        })
        .collect();
    Bloque::tabla("Bosquejos", &columnas, filas)
}

fn nombres_para_arbol(nombres: &[String]) -> Vec<String> {
    let mut limpios: Vec<String> = nombres
        .iter()
        .map(|nombre| {
            nombre
                .chars()
                .map(|letra| {
                    if letra.is_alphanumeric() || matches!(letra, '_' | '.' | '-') {
                        letra
                    } else {
                        '_'
                    }
                })
                .collect()
        })
        .collect();
    hacer_unicos(&mut limpios);
    limpios
}

fn arbol_de_vecinos(pedido: &Pedido, matriz: &[Vec<f64>]) -> Result<String, String> {
    let hojas = nombres_para_arbol(&pedido.nombres);
    let mut arbol = filogenia::union_de_vecinos(matriz, &hojas)?;
    for nodo in hojas.len()..arbol.nombres.len() {
        arbol.nombres[nodo] = String::new();
    }
    for vecinos in arbol.vecinos.iter_mut() {
        for (_, longitud) in vecinos.iter_mut() {
            *longitud = longitud.max(0.0);
        }
    }
    let raiz = (0..arbol.nombres.len())
        .rev()
        .find(|&nodo| arbol.grado(nodo) >= 2)
        .unwrap_or(0);
    Ok(filogenia::escribir_newick(&arbol, raiz))
}

fn ejecutar_minhash(entrada: &Entrada) -> Result<Salida, String> {
    let pedido = leer_pedido(entrada)?;
    let inicio = Instant::now();
    let comparacion = minhash::comparar(&pedido.genomas, pedido.k, pedido.cupo)?;
    let tiempo = inicio.elapsed();
    let exacto = calcular_exacto(&pedido);
    let matriz = minhash::matriz_de_distancias(&comparacion);
    let mut salida = Salida::default();
    salida.agregar(tabla_de_pares(&pedido, &comparacion, exacto.as_ref()));
    salida.agregar(Bloque::calor(
        "Distancia de Mash entre cada par (0 = iguales, 1 = nada en común)",
        pedido.nombres.clone(),
        pedido.nombres.clone(),
        matriz.clone(),
    ));
    salida.agregar(Bloque::arbol(
        "Árbol de unión de vecinos con las distancias de Mash (sin raíz)",
        arbol_de_vecinos(&pedido, &matriz)?,
    ));
    salida.agregar(tabla_de_bosquejos(&pedido, &comparacion, exacto.as_ref()));
    salida.agregar(Bloque::valor(
        "Parámetros",
        format!(
            "k = {} · {} hashes por bosquejo · {} pares",
            pedido.k,
            miles(pedido.cupo),
            miles(comparacion.pares.len())
        ),
    ));
    let mut tiempos = format!("MinHash: {}", texto_de_tiempo(tiempo));
    if let Some(exacto) = &exacto {
        tiempos.push_str(&format!(
            " · todos los k-meros (exacto): {}",
            texto_de_tiempo(exacto.tiempo)
        ));
    }
    salida.agregar(Bloque::valor("Tiempo", tiempos));
    interpretar(&mut salida, &pedido, &comparacion, exacto.as_ref(), tiempo);
    Ok(salida)
}

fn k_recomendado(longitud: usize) -> usize {
    let necesario = longitud as f64 * (1.0 - AZAR_TOLERADO) / AZAR_TOLERADO;
    (necesario.ln() / 4f64.ln()).ceil().max(1.0) as usize
}

fn frase_del_mejor(pedido: &Pedido, mejor: &Par, comparacion: &Comparacion) -> String {
    let (a, b) = (
        &pedido.nombres[mejor.primero],
        &pedido.nombres[mejor.segundo],
    );
    let medida = mejor.medida;
    if medida.compartidos == 0 {
        return format!(
            "Ningún par comparte hashes: con k = {} no hay parecido medible entre estas secuencias (todas las distancias valen 1).",
            pedido.k
        );
    }
    if medida.distancia == 0.0 {
        return format!(
            "«{a}» y «{b}» comparten todos los hashes de la muestra (distancia 0): con k = {} no se distinguen; pueden ser la misma secuencia o una leída en la otra hebra.",
            pedido.k
        );
    }
    let frontera = if medida.identidad >= FRONTERA_DE_ESPECIE {
        "por encima del 95 %, la frontera usual entre especies bacterianas (Jain y colegas, 2018)"
    } else {
        "por debajo del 95 % que suele separar especies bacterianas (Jain y colegas, 2018)"
    };
    let sobre_la_frontera = comparacion
        .pares
        .iter()
        .filter(|par| par.medida.compartidos > 0 && par.medida.identidad >= FRONTERA_DE_ESPECIE)
        .count();
    format!(
        "El par más parecido es «{a}» y «{b}»: comparten {} de {} hashes, distancia de Mash {} y ANI estimado {}, {frontera}; {} de {} pares quedan sobre esa frontera.",
        miles(medida.compartidos),
        miles(medida.vistos),
        decimal(medida.distancia, 4),
        porcentaje(medida.identidad, 2),
        miles(sobre_la_frontera),
        miles(comparacion.pares.len())
    )
}

fn frase_del_valor_p(pedido: &Pedido, comparacion: &Comparacion) -> String {
    let pares = comparacion.pares.len();
    let umbral = ALFA / pares as f64;
    let es_significativo =
        |par: &&Par| par.medida.compartidos > 0 && par.medida.log10_valor_p < umbral.log10();
    let significativos = comparacion.pares.iter().filter(es_significativo).count();
    let sin_nada = comparacion
        .pares
        .iter()
        .filter(|par| par.medida.compartidos == 0)
        .count();
    let dudosos = pares - significativos - sin_nada;
    let ruidosos = comparacion
        .pares
        .iter()
        .filter(es_significativo)
        .filter(|par| par.medida.identidad < ANI_CONFIABLE)
        .count();
    let mut frase = format!(
        "Con corrección de Bonferroni (0.05 entre {}, valor p < {}), {} más hashes de los que daría el azar",
        cantidad(pares, "par", "pares"),
        texto_de_valor_p(umbral.log10()),
        cantidad(significativos, "par comparte", "pares comparten")
    );
    if sin_nada > 0 {
        frase.push_str(&format!(
            "; {} ninguno (distancia 1: sin parecido medible con k = {})",
            cantidad(sin_nada, "no comparte", "no comparten"),
            pedido.k
        ));
    }
    if dudosos > 0 {
        frase.push_str(&format!(
            "; {} tan pocos que podría ser azar",
            cantidad(dudosos, "comparte", "comparten")
        ));
    }
    if ruidosos > 0 {
        frase.push_str(&format!(
            "; en {} el ANI estimado baja de 90 %, donde quedan pocos hashes en común y Mash se vuelve impreciso",
            miles(ruidosos)
        ));
    }
    frase.push('.');
    frase
}

fn jaccard_del_par(comparacion: &Comparacion, a: usize, b: usize) -> f64 {
    comparacion
        .pares
        .iter()
        .find(|par| (par.primero, par.segundo) == (a.min(b), a.max(b)))
        .map(|par| par.medida.jaccard)
        .unwrap_or(0.0)
}

fn par_contenido(pedido: &Pedido, comparacion: &Comparacion) -> Option<(usize, usize)> {
    let genomas_en_total = comparacion.bosquejos.len();
    let mut mejor: Option<(usize, usize, f64, f64)> = None;
    for chico in 0..genomas_en_total {
        for grande in 0..genomas_en_total {
            let contencion = comparacion.contencion[chico][grande];
            if chico == grande
                || pedido.bases(chico) >= pedido.bases(grande)
                || contencion < CASI_ENTERA
                || comparacion.contencion[grande][chico] >= CASI_ENTERA
            {
                continue;
            }
            let jaccard = jaccard_del_par(comparacion, chico, grande);
            let es_mejor = mejor.is_none_or(|(_, _, valor, otro_jaccard)| {
                contencion > valor || (contencion == valor && jaccard < otro_jaccard)
            });
            if es_mejor {
                mejor = Some((chico, grande, contencion, jaccard));
            }
        }
    }
    mejor.map(|(chico, grande, _, _)| (chico, grande))
}

fn par_que_se_encima(comparacion: &Comparacion) -> Option<&Par> {
    comparacion
        .pares
        .iter()
        .filter(|par| {
            let ida = comparacion.contencion[par.primero][par.segundo];
            let vuelta = comparacion.contencion[par.segundo][par.primero];
            par.medida.identidad < 0.99
                && ida.min(vuelta) >= TRAMO_MINIMO
                && ida.max(vuelta) < CASI_ENTERA
        })
        .max_by(|uno, otro| {
            let menor = |par: &Par| {
                comparacion.contencion[par.primero][par.segundo]
                    .min(comparacion.contencion[par.segundo][par.primero])
            };
            menor(uno).total_cmp(&menor(otro))
        })
}

fn frase_de_precision(
    pedido: &Pedido,
    comparacion: &Comparacion,
    ordenados: &[Par],
    exacto: Option<&Exacto>,
) -> String {
    let mut frase = format!(
        "Con {} por bosquejo",
        cantidad(pedido.cupo, "hash", "hashes")
    );
    let medido = ordenados
        .iter()
        .find(|par| par.medida.compartidos > 0 && par.medida.compartidos < par.medida.vistos);
    if let Some(par) = medido {
        let jaccard = par.medida.jaccard;
        let error_jaccard = (jaccard * (1.0 - jaccard) / par.medida.vistos as f64).sqrt();
        let error_ani = error_jaccard / (pedido.k as f64 * jaccard * (1.0 + jaccard));
        frase.push_str(&format!(
            ", el error típico del Jaccard de «{}» y «{}» es √(J(1−J)/s) ≈ {}, unos ±{} puntos de ANI",
            pedido.nombres[par.primero],
            pedido.nombres[par.segundo],
            decimal(error_jaccard, 4),
            decimal(error_ani * 100.0, 2)
        ));
    }
    frase.push_str("; más hashes dan más precisión, a cambio de más memoria");
    match exacto {
        Some(exacto) => {
            let diferencia = comparacion
                .pares
                .iter()
                .map(|par| {
                    (par.medida.jaccard - exacto.jaccard[par.primero][par.segundo]).abs()
                })
                .fold(0.0, f64::max);
            frase.push_str(&format!(
                ", y contra el Jaccard exacto con todos los k-meros la mayor diferencia fue {}",
                decimal(diferencia, 4)
            ));
        }
        None => frase.push_str(
            " (el Jaccard exacto no se calculó porque las secuencias pasan de 12 millones de letras)",
        ),
    }
    let mayor = (0..pedido.genomas.len())
        .map(|indice| pedido.bases(indice))
        .max()
        .unwrap_or(0);
    let recomendado = k_recomendado(mayor);
    if pedido.k < recomendado {
        frase.push_str(&format!(
            "; ojo: con k = {} y {} bases, un k-mero cualquiera aparece por azar con más de 1 % de probabilidad, así que conviene k ≥ {recomendado}",
            pedido.k,
            miles(mayor)
        ));
    }
    frase.push('.');
    frase
}

fn interpretar(
    salida: &mut Salida,
    pedido: &Pedido,
    comparacion: &Comparacion,
    exacto: Option<&Exacto>,
    tiempo: Duration,
) {
    let ordenados = pares_ordenados(comparacion);
    let Some(mejor) = ordenados.first() else {
        return;
    };
    salida.explicar(frase_del_mejor(pedido, mejor, comparacion));
    salida.explicar(frase_del_valor_p(pedido, comparacion));
    if let Some((chico, grande)) = par_contenido(pedido, comparacion) {
        let contencion = comparacion.contencion[chico][grande];
        salida.explicar(format!(
            "«{}» está casi entera dentro de «{}» (contención {}, identidad ≈ {}), aunque su Jaccard es solo {}: la distancia de Mash supone dos genomas de tamaño parecido y castiga la diferencia de tamaño, así que para un trozo dentro de un genoma lee la contención (Mash Screen).",
            pedido.nombres[chico],
            pedido.nombres[grande],
            decimal(contencion, 3),
            porcentaje(minhash::identidad_de_contencion(contencion, pedido.k), 2),
            decimal(jaccard_del_par(comparacion, chico, grande), 3)
        ));
    }
    if let Some(par) = par_que_se_encima(comparacion) {
        let ida = comparacion.contencion[par.primero][par.segundo];
        let vuelta = comparacion.contencion[par.segundo][par.primero];
        salida.explicar(format!(
            "«{}» y «{}» dan un ANI estimado de {}, pero ese número supone que divergieron por mutaciones a lo largo de toda su longitud; si solo comparten un tramo (aquí {} de una está en la otra y {} al revés), ese tramo puede ser idéntico y la contención lo dice mejor.",
            pedido.nombres[par.primero],
            pedido.nombres[par.segundo],
            porcentaje(par.medida.identidad, 2),
            porcentaje(ida, 0),
            porcentaje(vuelta, 0)
        ));
    }
    salida.explicar(frase_de_precision(pedido, comparacion, &ordenados, exacto));
    let mut frase = format!(
        "Frente a «Distancias p y Jukes–Cantor» o «Alineamiento global», que comparan base por base secuencias alineadas, MinHash no alinea nada: guarda {} hashes (unos {} KB) por genoma y compara cada par en a lo más {} pasos; aquí todo tardó {}",
        miles(pedido.cupo),
        miles((pedido.cupo * 8).div_ceil(1024)),
        miles(2 * pedido.cupo),
        texto_de_tiempo(tiempo)
    );
    if let Some(exacto) = exacto {
        frase.push_str(&format!(
            ", y contar todos los k-meros para el valor exacto, {}",
            texto_de_tiempo(exacto.tiempo)
        ));
    }
    frase.push('.');
    salida.explicar(frase);
}

struct Contexto {
    k: usize,
    cupo: usize,
    cifras: usize,
}

#[derive(Default)]
struct Recuento {
    kmeros: usize,
    entraron: usize,
    salieron: usize,
    rechazados: usize,
    repetidos: usize,
    letras_raras: usize,
}

fn decodificar(codigo: u64, largo: usize) -> String {
    let mut texto = String::with_capacity(largo);
    for posicion in (0..largo).rev() {
        texto.push(b"ACGT"[((codigo >> (2 * posicion)) & 3) as usize] as char);
    }
    texto
}

fn adelante_en_texto(ventana: &Ventana) -> String {
    decodificar(ventana.adelante, ventana.validas.min(ventana.k))
}

fn atras_en_texto(ventana: &Ventana) -> String {
    let largo = ventana.validas.min(ventana.k);
    if largo == 0 {
        return String::new();
    }
    decodificar(ventana.atras >> (2 * (ventana.k - largo)), largo)
}

fn umbral_en_texto(menores: &Menores, cifras: usize) -> String {
    if menores.conjunto.len() < menores.cupo {
        "sin tope (aún no se llena)".to_string()
    } else {
        hash_en_texto(menores.mayor_guardado, cifras)
    }
}

fn vista_de_menores(
    menores: &Menores,
    resaltado: Option<u64>,
    cifras: usize,
) -> (Vec<String>, Vec<usize>) {
    let ordenados: Vec<u64> = menores.conjunto.iter().copied().collect();
    let lugar = resaltado.and_then(|hash| ordenados.binary_search(&hash).ok());
    let mut indices: Vec<usize> = Vec::new();
    let total = ordenados.len();
    if total <= 12 {
        indices.extend(0..total);
    } else {
        indices.extend(0..3);
        if let Some(lugar) = lugar {
            indices.extend(lugar.saturating_sub(2)..(lugar + 3).min(total));
        }
        indices.extend(total - 3..total);
        indices.sort_unstable();
        indices.dedup();
    }
    let mut elementos = Vec::new();
    let mut marcados = Vec::new();
    let mut anterior: Option<usize> = None;
    for indice in indices {
        if anterior.is_some_and(|antes| indice > antes + 1) || (anterior.is_none() && indice > 0) {
            elementos.push("…".to_string());
        }
        if Some(indice) == lugar {
            marcados.push(elementos.len());
        }
        elementos.push(format!(
            "{}) {}",
            miles(indice + 1),
            hash_en_texto(ordenados[indice], cifras)
        ));
        anterior = Some(indice);
    }
    (elementos, marcados)
}

fn pizarra_de_kmero(
    ventana: &Ventana,
    menores: &Menores,
    hash: Option<u64>,
    recuento: &Recuento,
    contexto: &Contexto,
) -> Pizarra {
    Pizarra::nueva()
        .variable("ventana.validas", ventana.validas)
        .variable("ventana.adelante", adelante_en_texto(ventana))
        .variable("ventana.atras", atras_en_texto(ventana))
        .variable(
            "hash",
            hash.map(|valor| hash_en_texto(valor, contexto.cifras))
                .unwrap_or_else(|| "—".to_string()),
        )
        .variable("kmeros", miles(recuento.kmeros))
        .variable("menores.conjunto.len()", miles(menores.conjunto.len()))
        .variable(
            "menores.mayor_guardado",
            umbral_en_texto(menores, contexto.cifras),
        )
}

fn con_lista(pizarra: Pizarra, menores: &Menores, hash: Option<u64>, cifras: usize) -> Pizarra {
    let (elementos, marcados) = vista_de_menores(menores, hash, cifras);
    pizarra.lista(
        &format!(
            "Los hashes más chicos guardados ({} de {})",
            miles(menores.conjunto.len()),
            miles(menores.cupo)
        ),
        elementos,
        &marcados,
    )
}

fn narrar_primera_base(
    narracion: &mut Narracion,
    ventana: &Ventana,
    base: u8,
    pieza: &[u8],
    contexto: &Contexto,
) {
    let letra = base as char;
    let codigo = match base.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        _ => 3,
    };
    let marcas = [(0, 1, Marca::Actual)];
    let complemento = b"TGCA"[codigo] as char;
    narracion.agregar(
        Pizarra::nueva()
            .variable("base", letra)
            .variable("codigo", codigo)
            .secuencia(pieza, &marcas)
            .paso(
                L_CODIGO,
                format!(
                    "Lees la primera base, «{letra}». La tabla CODIGO_DE_BASE la vuelve un número de 2 bits: A = 0, C = 1, G = 2, T = 3; cualquier otra letra rompe la ventana."
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("base", letra)
            .variable("codigo", codigo)
            .variable("ventana.adelante", adelante_en_texto(ventana))
            .secuencia(pieza, &marcas)
            .paso(
                L_ADELANTE,
                "adelante corre dos bits a la izquierda y mete el código por la derecha: guarda el k-mero leído de izquierda a derecha, 2 bits por base.",
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("base", letra)
            .variable("codigo", codigo)
            .variable("ventana.adelante", adelante_en_texto(ventana))
            .variable("ventana.atras", atras_en_texto(ventana))
            .secuencia(pieza, &marcas)
            .paso(
                L_ATRAS,
                format!(
                    "atras guarda al mismo tiempo el complemento inverso: el complemento de {letra} ({complemento}) entra por la izquierda. Así nunca hay que voltear la secuencia."
                ),
            ),
    );
    if contexto.k < 2 {
        return;
    }
    let siguientes = if contexto.k > 2 {
        format!(
            " Lo mismo pasa con las siguientes {} bases.",
            contexto.k - 2
        )
    } else {
        String::new()
    };
    narracion.agregar(
        Pizarra::nueva()
            .variable("ventana.validas", ventana.validas)
            .variable("k", contexto.k)
            .secuencia(pieza, &marcas)
            .paso(
                L_INCOMPLETO,
                format!(
                    "Llevas {} base válida y k = {}: todavía no hay un k-mero completo, así que no sale ningún hash.{siguientes}",
                    ventana.validas,
                    contexto.k
                ),
            ),
    );
}

fn narrar_primer_kmero(
    narracion: &mut Narracion,
    ventana: &Ventana,
    hash: u64,
    pieza: &[u8],
    posicion: usize,
    contexto: &Contexto,
) {
    let inicio = posicion + 1 - contexto.k;
    let marcas = [(inicio, posicion + 1, Marca::Actual)];
    let adelante = adelante_en_texto(ventana);
    let atras = atras_en_texto(ventana);
    let elegido = if ventana.adelante <= ventana.atras {
        "adelante"
    } else {
        "atras"
    };
    narracion.agregar(
        Pizarra::nueva()
            .variable("ventana.validas", ventana.validas)
            .variable("ventana.adelante", &adelante)
            .variable("ventana.atras", &atras)
            .secuencia(pieza, &marcas)
            .paso(
                L_CANONICO,
                format!(
                    "Con {} ya hay k-mero: {adelante}, y su complemento inverso {atras}. El canónico es el menor de los dos como número (aquí {elegido}); así la misma secuencia leída en la otra hebra da el mismo hash.",
                    cantidad(contexto.k, "base", "bases")
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("ventana.adelante", &adelante)
            .variable("ventana.atras", &atras)
            .variable("hash", hash_en_texto(hash, contexto.cifras))
            .secuencia(pieza, &marcas)
            .paso(
                L_MEZCLAR,
                format!(
                    "mezclar revuelve los 64 bits (el paso final de MurmurHash3) para que el número parezca tirado al azar. Te lo muestro dividido entre 2^64: este k-mero cae en {} de la recta de 0 a 1. Es una biyección, así que dos k-meros distintos nunca chocan.",
                    hash_en_texto(hash, contexto.cifras)
                ),
            ),
    );
}

struct Evento {
    hash: u64,
    es_nuevo: bool,
    lleno_antes: bool,
    umbral_antes: u64,
    largo_antes: usize,
}

fn nota_de_evento(
    evento: &Evento,
    menores: &Menores,
    contexto: &Contexto,
) -> Option<(&'static str, String)> {
    let hash = hash_en_texto(evento.hash, contexto.cifras);
    if evento.lleno_antes && evento.hash >= evento.umbral_antes {
        return Some((
            L_RECHAZO,
            format!(
                "Hash {hash} ≥ umbral {}: se descarta con una sola comparación. Eso hace rápido a MinHash: casi todos los k-meros terminan aquí.",
                hash_en_texto(evento.umbral_antes, contexto.cifras)
            ),
        ));
    }
    if !evento.es_nuevo {
        return Some((
            L_INSERTAR,
            "Este k-mero ya había salido antes: insert avisa que no es nuevo y el conjunto no cambia (el bosquejo cuenta k-meros distintos).".to_string(),
        ));
    }
    let lugar = menores.conjunto.range(..evento.hash).count() + 1;
    if evento.largo_antes >= contexto.cupo {
        return Some((
            L_SACAR,
            format!(
                "Hash {hash} < umbral {}: entra en el lugar {} y, como ya había {}, sale el mayor. El nuevo umbral es {}.",
                hash_en_texto(evento.umbral_antes, contexto.cifras),
                miles(lugar),
                miles(contexto.cupo),
                hash_en_texto(menores.mayor_guardado, contexto.cifras)
            ),
        ));
    }
    if menores.conjunto.len() == contexto.cupo {
        return Some((
            L_UMBRAL,
            format!(
                "Con este ya hay {} hashes guardados: el mayor, {}, es ahora el umbral. Desde aquí un hash solo entra si es menor que él.",
                miles(contexto.cupo),
                hash_en_texto(menores.mayor_guardado, contexto.cifras)
            ),
        ));
    }
    Some((
        L_INSERTAR,
        format!(
            "Hash {hash}. Todavía hay menos de {} guardados, así que entra sin competir, en el lugar {} de {} de la lista ordenada.",
            miles(contexto.cupo),
            miles(lugar),
            miles(menores.conjunto.len())
        ),
    ))
}

fn debe_narrar(evento: &Evento, recuento: &Recuento, menores: &Menores, cupo: usize) -> bool {
    if evento.lleno_antes && evento.hash >= evento.umbral_antes {
        return recuento.rechazados <= RECHAZOS_DETALLADOS;
    }
    if !evento.es_nuevo {
        return recuento.repetidos <= REPETIDOS_DETALLADOS;
    }
    if evento.largo_antes >= cupo {
        return recuento.salieron <= SALIDAS_DETALLADAS;
    }
    menores.conjunto.len() == cupo || recuento.kmeros <= KMEROS_DETALLADOS
}

fn anotar(evento: &Evento, recuento: &mut Recuento, cupo: usize) {
    if evento.lleno_antes && evento.hash >= evento.umbral_antes {
        recuento.rechazados += 1;
    } else if !evento.es_nuevo {
        recuento.repetidos += 1;
    } else {
        recuento.entraron += 1;
        if evento.largo_antes >= cupo {
            recuento.salieron += 1;
        }
    }
}

fn narrar_control(
    narracion: &mut Narracion,
    leidas: usize,
    total: usize,
    recuento: &Recuento,
    menores: &Menores,
    contexto: &Contexto,
) {
    narracion.agregar(
        Pizarra::nueva()
            .variable("kmeros", miles(recuento.kmeros))
            .variable("menores.conjunto.len()", miles(menores.conjunto.len()))
            .variable(
                "menores.mayor_guardado",
                umbral_en_texto(menores, contexto.cifras),
            )
            .paso(
                L_OFRECER,
                format!(
                    "Vas en la base {} de {}: {} k-meros leídos; {} hashes entraron alguna vez, {} volvieron a salir y {} se descartaron con una comparación. El umbral solo baja: cada vez es más difícil entrar.",
                    miles(leidas),
                    miles(total),
                    miles(recuento.kmeros),
                    miles(recuento.entraron),
                    miles(recuento.salieron),
                    miles(recuento.rechazados)
                ),
            ),
    );
}

fn narrar_primer_bosquejo(
    narracion: &mut Narracion,
    nombre: &str,
    piezas: &[Vec<u8>],
    contexto: &Contexto,
) -> Vec<u64> {
    let mut ventana = Ventana::nueva(contexto.k);
    let mut menores = Menores::nuevo(contexto.cupo);
    let total: usize = piezas.iter().map(Vec::len).sum();
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", contexto.k)
            .variable("cupo", miles(contexto.cupo))
            .paso(
                L_MENORES,
                format!(
                    "Empiezas con «{nombre}» ({} bases). Vas a pasar por cada k-mero, convertirlo en un hash y quedarte solo con los {} hashes más chicos: ese es el bosquejo.",
                    miles(total),
                    miles(contexto.cupo)
                ),
            ),
    );
    let mut recuento = Recuento::default();
    let mut leidas = 0;
    let mut control = 1;
    for (numero, pieza) in piezas.iter().enumerate() {
        ventana.reiniciar();
        for (posicion, &base) in pieza.iter().enumerate() {
            leidas += 1;
            let salida = ventana.agregar(base);
            if numero == 0 && posicion == 0 && ventana.validas == 1 {
                narrar_primera_base(narracion, &ventana, base, pieza, contexto);
            }
            if control <= 3 && total >= 400 && leidas * 4 >= total * control {
                control += 1;
                narrar_control(narracion, leidas, total, &recuento, &menores, contexto);
            }
            let Some(hash) = salida else {
                if ventana.validas == 0 {
                    recuento.letras_raras += 1;
                    if recuento.letras_raras <= LETRAS_RARAS_DETALLADAS {
                        narracion.agregar(
                            pizarra_de_kmero(&ventana, &menores, None, &recuento, contexto)
                                .secuencia(pieza, &[(posicion, posicion + 1, Marca::Actual)])
                                .paso(
                                    L_NO_ES_BASE,
                                    format!(
                                        "La letra «{}» no es A, C, G ni T: validas vuelve a 0 y hay que juntar {} bases nuevas antes del siguiente k-mero.",
                                        base as char,
                                        contexto.k
                                    ),
                                ),
                        );
                    }
                }
                continue;
            };
            recuento.kmeros += 1;
            if recuento.kmeros == 1 {
                narrar_primer_kmero(narracion, &ventana, hash, pieza, posicion, contexto);
            }
            let evento = Evento {
                hash,
                es_nuevo: false,
                lleno_antes: menores.conjunto.len() >= contexto.cupo,
                umbral_antes: menores.mayor_guardado,
                largo_antes: menores.conjunto.len(),
            };
            let es_nuevo = menores.ofrecer(hash);
            let evento = Evento { es_nuevo, ..evento };
            anotar(&evento, &mut recuento, contexto.cupo);
            if debe_narrar(&evento, &recuento, &menores, contexto.cupo)
                && let Some((linea, nota)) = nota_de_evento(&evento, &menores, contexto)
            {
                let inicio = posicion + 1 - contexto.k;
                narracion.agregar(
                    con_lista(
                        pizarra_de_kmero(&ventana, &menores, Some(hash), &recuento, contexto),
                        &menores,
                        Some(hash),
                        contexto.cifras,
                    )
                    .secuencia(pieza, &[(inicio, posicion + 1, Marca::Actual)])
                    .paso(linea, nota),
                );
            }
        }
    }
    let hashes: Vec<u64> = menores.conjunto.iter().copied().collect();
    let estimado = minhash::estimar_kmeros_distintos(&Bosquejo {
        k: contexto.k,
        cupo: contexto.cupo,
        hashes: hashes.clone(),
        longitud: total,
        kmeros: recuento.kmeros,
    });
    let mayor = hashes
        .last()
        .map(|&hash| hash_en_texto(hash, contexto.cifras))
        .unwrap_or_default();
    let lleno = hashes.len() >= contexto.cupo;
    let cierre = if lleno {
        format!(
            "Como el hash reparte los k-meros al azar, los más chicos son una muestra pareja de todos: con el mayor guardado ({mayor}) estimas unos {} k-meros distintos sin haberlos guardado.",
            miles(estimado.round() as usize)
        )
    } else {
        format!(
            "No llegó a {} k-meros distintos, así que el bosquejo los guarda todos: aquí la comparación será exacta.",
            miles(contexto.cupo)
        )
    };
    narracion.agregar(
        con_lista(
            Pizarra::nueva()
                .variable("kmeros", miles(recuento.kmeros))
                .variable("menores.conjunto.len()", miles(hashes.len())),
            &menores,
            None,
            contexto.cifras,
        )
        .paso(
            L_FIN_BOSQUEJO,
            format!(
                "Bosquejo de «{nombre}» listo: {} hashes de {} k-meros leídos ({} repetidos). {cierre}",
                miles(hashes.len()),
                miles(recuento.kmeros),
                miles(recuento.repetidos)
            ),
        ),
    );
    hashes
}

fn filas_de_bosquejos(
    pedido: &Pedido,
    comparacion: &Comparacion,
    cifras: usize,
) -> Vec<Vec<String>> {
    comparacion
        .bosquejos
        .iter()
        .enumerate()
        .map(|(indice, bosquejo)| {
            vec![
                pedido.nombres[indice].clone(),
                miles(bosquejo.longitud),
                miles(bosquejo.hashes.len()),
                bosquejo
                    .hashes
                    .last()
                    .map(|&hash| hash_en_texto(hash, cifras))
                    .unwrap_or_default(),
            ]
        })
        .collect()
}

fn par_para_mezclar(comparacion: &Comparacion) -> (usize, usize) {
    comparacion
        .pares
        .iter()
        .filter(|par| par.medida.compartidos < par.medida.vistos)
        .max_by_key(|par| par.medida.compartidos)
        .or(comparacion.pares.first())
        .map(|par| (par.primero, par.segundo))
        .unwrap_or((0, 1))
}

fn ventana_de_hashes(lista: &[u64], desde: usize, cifras: usize, fila: usize) -> String {
    lista
        .get(desde + fila)
        .map(|&hash| hash_en_texto(hash, cifras))
        .unwrap_or_default()
}

fn pizarra_de_mezcla(
    a: &[u64],
    b: &[u64],
    posiciones: (usize, usize),
    cuentas: (usize, usize),
    cifras: usize,
) -> Pizarra {
    let (i, j) = posiciones;
    let (compartidos, vistos) = cuentas;
    let filas = (0..FILAS_EN_PIZARRA)
        .map(|fila| {
            vec![
                ventana_de_hashes(a, i, cifras, fila),
                ventana_de_hashes(b, j, cifras, fila),
            ]
        })
        .collect();
    Pizarra::nueva()
        .variable("i", i)
        .variable("j", j)
        .variable("a[i]", ventana_de_hashes(a, i, cifras, 0))
        .variable("b[j]", ventana_de_hashes(b, j, cifras, 0))
        .variable("compartidos", miles(compartidos))
        .variable("vistos", miles(vistos))
        .tabla(
            "Lo que falta de cada bosquejo (desde i y desde j)",
            &["a (bosquejo de A)", "b (bosquejo de B)"],
            filas,
            &[(0, 0), (0, 1)],
        )
}

fn narrar_mezcla(
    narracion: &mut Narracion,
    a: &Bosquejo,
    b: &Bosquejo,
    nombres: (&str, &str),
    contexto: &Contexto,
) -> (usize, usize) {
    let cupo = a.cupo.min(b.cupo);
    let (x, y) = (&a.hashes, &b.hashes);
    narracion.agregar(
        pizarra_de_mezcla(x, y, (0, 0), (0, 0), contexto.cifras).paso(
            L_CONTAR,
            format!(
                "Comparas «{}» (A) con «{}» (B). contar_compartidos recorre los dos bosquejos ordenados a la vez, como al juntar dos mazos ordenados, hasta ver {} hashes distintos de la unión.",
                nombres.0,
                nombres.1,
                miles(cupo)
            ),
        ),
    );
    let mut i = 0;
    let mut j = 0;
    let mut compartidos = 0;
    let mut vistos = 0;
    let mut vueltas = 0;
    let mut compartidos_narrados = 0;
    let mut control = 1;
    while vistos < cupo && i < x.len() && j < y.len() {
        let (antes_i, antes_j) = (i, j);
        let (izquierda, derecha) = (
            hash_en_texto(x[i], contexto.cifras),
            hash_en_texto(y[j], contexto.cifras),
        );
        let orden = x[i].cmp(&y[j]);
        match orden {
            Ordering::Equal => {
                compartidos += 1;
                i += 1;
                j += 1;
            }
            Ordering::Less => i += 1,
            Ordering::Greater => j += 1,
        }
        vistos += 1;
        vueltas += 1;
        let es_igual = orden == Ordering::Equal;
        let narrar = vueltas <= VUELTAS_DETALLADAS
            || (es_igual && compartidos_narrados < COMPARTIDOS_DETALLADOS);
        if narrar {
            if es_igual && vueltas > VUELTAS_DETALLADAS {
                compartidos_narrados += 1;
            }
            let (linea, nota) = match orden {
                Ordering::Equal => (
                    L_IGUAL,
                    format!(
                        "a[i] = b[j] = {izquierda}: este hash está en los dos bosquejos. compartidos sube a {} y avanzas en ambos.",
                        miles(compartidos)
                    ),
                ),
                Ordering::Less => (
                    L_MENOR,
                    format!(
                        "a[i] = {izquierda} es menor que b[j] = {derecha}: ese hash solo está en A. Avanzas i; cuenta para la unión, pero no como compartido."
                    ),
                ),
                Ordering::Greater => (
                    L_MAYOR,
                    format!(
                        "b[j] = {derecha} es menor que a[i] = {izquierda}: ese hash solo está en B. Avanzas j; cuenta para la unión, pero no como compartido."
                    ),
                ),
            };
            narracion.agregar(
                pizarra_de_mezcla(
                    x,
                    y,
                    (antes_i, antes_j),
                    (compartidos, vistos),
                    contexto.cifras,
                )
                .paso(linea, nota),
            );
        }
        if control <= 3 && cupo >= 40 && vistos * 4 >= cupo * control {
            control += 1;
            narracion.agregar(
                pizarra_de_mezcla(x, y, (i, j), (compartidos, vistos), contexto.cifras).paso(
                    L_VISTOS,
                    format!(
                        "Ya viste {} hashes de la unión y {} estaban en los dos. Por ahora la fracción compartida es {}.",
                        miles(vistos),
                        miles(compartidos),
                        decimal(minhash::fraccion(compartidos, vistos), 3)
                    ),
                ),
            );
        }
    }
    let sobrantes = (x.len() - i) + (y.len() - j);
    let agregados = sobrantes.min(cupo - vistos);
    vistos += agregados;
    let (linea, nota) = if agregados > 0 {
        (
            L_SOBRANTES,
            format!(
                "Un bosquejo se acabó antes de llegar a {}: los {} hashes que quedan en el otro también son de la unión (y no compartidos), así que suman a vistos.",
                miles(cupo),
                miles(agregados)
            ),
        )
    } else {
        (
            L_MIENTRAS,
            format!(
                "vistos llegó a {}: ya tienes los {} hashes más chicos de la unión, una muestra al azar de todos los k-meros de A y B juntos.",
                miles(vistos),
                miles(vistos)
            ),
        )
    };
    narracion.agregar(
        pizarra_de_mezcla(x, y, (i, j), (compartidos, vistos), contexto.cifras).paso(linea, nota),
    );
    (compartidos, vistos)
}

fn narrar_medidas(
    narracion: &mut Narracion,
    a: &Bosquejo,
    b: &Bosquejo,
    cuentas: (usize, usize),
    nombres: (&str, &str),
) {
    let (compartidos, vistos) = cuentas;
    let jaccard = minhash::fraccion(compartidos, vistos);
    let distancia = minhash::distancia_mash(jaccard, a.k);
    let pa = minhash::probabilidad_de_kmero_al_azar(a.longitud, a.k);
    let azar = minhash::jaccard_al_azar(a.longitud, b.longitud, a.k);
    let log10 = minhash::log10_valor_p(compartidos, vistos, azar);
    let base = || {
        Pizarra::nueva()
            .variable("compartidos", miles(compartidos))
            .variable("vistos", miles(vistos))
    };
    narracion.agregar(base().variable("jaccard", decimal(jaccard, 4)).paso(
        L_JACCARD,
        format!(
            "Jaccard estimado = {} / {} = {}: la fracción de la muestra de la unión que está en las dos secuencias.",
            miles(compartidos),
            miles(vistos),
            decimal(jaccard, 4)
        ),
    ));
    let nota_distancia = if compartidos == 0 {
        "No comparten ningún hash: la distancia de Mash vale 1, sin parecido medible con este k."
            .to_string()
    } else {
        format!(
            "D = −ln(2J / (1 + J)) / k = {}. Si cada base cambia con probabilidad d, un k-mero sobrevive con (1 − d)^k; despejando d de esa idea sale esta fórmula.",
            decimal(distancia, 4)
        )
    };
    narracion.agregar(
        base()
            .variable("jaccard", decimal(jaccard, 4))
            .variable("distancia", decimal(distancia, 4))
            .paso(L_DISTANCIA, nota_distancia),
    );
    narracion.agregar(
        base()
            .variable("jaccard", decimal(jaccard, 4))
            .variable("distancia", decimal(distancia, 4))
            .variable("identidad", porcentaje(1.0 - distancia, 2))
            .paso(
                L_IDENTIDAD,
                format!(
                    "La identidad promedio de nucleótidos (ANI) estimada es 1 − D = {}.",
                    porcentaje(1.0 - distancia, 2)
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("longitud", miles(a.longitud))
            .variable("k", a.k)
            .paso(
                L_PROBABILIDAD,
                format!(
                    "Ahora el valor p. Un k-mero cualquiera aparece por azar en «{}» ({} bases) con probabilidad 1 / (1 + 4^{} / {}) = {}.",
                    nombres.0,
                    miles(a.longitud),
                    a.k,
                    miles(a.longitud),
                    texto_de_valor_p(pa.log10())
                ),
            ),
    );
    narracion.agregar(
        base().variable("azar", texto_de_valor_p(azar.log10())).paso(
            L_AZAR,
            format!(
                "Haciendo lo mismo con «{}», dos secuencias al azar de esos tamaños tendrían un Jaccard de r = {}.",
                nombres.1,
                texto_de_valor_p(azar.log10())
            ),
        ),
    );
    let veredicto = if log10 < ALFA.log10() {
        "Es tan chico que este parecido no es casualidad."
    } else {
        "No es chico: este parecido podría ser casualidad."
    };
    narracion.agregar(
        base()
            .variable("azar", texto_de_valor_p(azar.log10()))
            .variable("log10_valor_p", decimal(log10, 1))
            .paso(
                L_VALOR_P,
                format!(
                    "Valor p = probabilidad de ver {} o más compartidos de {} si cada uno coincidiera solo con probabilidad r (cola binomial): {}. {veredicto}",
                    miles(compartidos),
                    miles(vistos),
                    texto_de_valor_p(log10)
                ),
            ),
    );
}

fn narrar_pares(narracion: &mut Narracion, pedido: &Pedido, comparacion: &Comparacion) {
    let ordenados = pares_ordenados(comparacion);
    let filas: Vec<Vec<String>> = ordenados
        .iter()
        .take(PARES_EN_PIZARRA)
        .map(|par| {
            vec![
                format!(
                    "{} · {}",
                    pedido.nombres[par.primero], pedido.nombres[par.segundo]
                ),
                format!(
                    "{} de {}",
                    miles(par.medida.compartidos),
                    miles(par.medida.vistos)
                ),
                decimal(par.medida.distancia, 4),
                ani_en_texto(par),
            ]
        })
        .collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("pares", miles(comparacion.pares.len()))
            .tabla(
                "Pares, del más parecido al menos",
                &["par", "compartidos", "distancia", "ANI"],
                filas,
                &[(0, 2)],
            )
            .paso(
                L_MEDIR,
                format!(
                    "Lo mismo con los {} pares, en paralelo: cada uno cuesta a lo más {} pasos, sin importar el largo de los genomas.",
                    miles(comparacion.pares.len()),
                    miles(2 * pedido.cupo)
                ),
            ),
    );
}

fn par_de_contencion(pedido: &Pedido, comparacion: &Comparacion) -> (usize, usize) {
    par_contenido(pedido, comparacion).unwrap_or_else(|| {
        let mut mejor = (0, 1);
        let mut valor = -1.0;
        for chico in 0..comparacion.bosquejos.len() {
            for grande in 0..comparacion.bosquejos.len() {
                if chico != grande && comparacion.contencion[chico][grande] > valor {
                    valor = comparacion.contencion[chico][grande];
                    mejor = (chico, grande);
                }
            }
        }
        mejor
    })
}

fn narrar_contencion(narracion: &mut Narracion, pedido: &Pedido, comparacion: &Comparacion) {
    let mut todos: Vec<u64> = comparacion
        .bosquejos
        .iter()
        .flat_map(|bosquejo| bosquejo.hashes.iter().copied())
        .collect();
    todos.sort_unstable();
    todos.dedup();
    narracion.agregar(
        Pizarra::nueva().variable("todos.len()", miles(todos.len())).paso(
            L_TODOS,
            format!(
                "Falta la contención: qué fracción de una secuencia está dentro de otra. Juntas los hashes de todos los bosquejos en una lista ordenada: {} distintos.",
                miles(todos.len())
            ),
        ),
    );
    let (chico, grande) = par_de_contencion(pedido, comparacion);
    let contencion = comparacion.contencion[chico][grande];
    let bosquejo = &comparacion.bosquejos[chico];
    let hallados = (contencion * bosquejo.hashes.len() as f64).round() as usize;
    narracion.agregar(
        Pizarra::nueva()
            .variable("todos.len()", miles(todos.len()))
            .variable("hallados", miles(hallados))
            .paso(
                L_PRESENTE,
                format!(
                    "Recorres cada k-mero de «{}» (sus {} bases, no su bosquejo) y marcas los hashes de la lista que aparecen. De los {} hashes del bosquejo de «{}», {} quedaron marcados.",
                    pedido.nombres[grande],
                    miles(pedido.bases(grande)),
                    miles(bosquejo.hashes.len()),
                    pedido.nombres[chico],
                    miles(hallados)
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("contencion", decimal(contencion, 3))
            .paso(
                L_FRACCION_PRESENTE,
                format!(
                    "Contención de «{}» en «{}» = {} / {} = {}: como el bosquejo es una muestra al azar de sus k-meros, esa fracción estima cuánto de «{}» está en «{}» (Mash Screen).",
                    pedido.nombres[chico],
                    pedido.nombres[grande],
                    miles(hallados),
                    miles(bosquejo.hashes.len()),
                    decimal(contencion, 3),
                    pedido.nombres[chico],
                    pedido.nombres[grande]
                ),
            ),
    );
}

fn trazar(pedido: &Pedido) -> Result<Vec<Paso>, String> {
    let comparacion = minhash::comparar(&pedido.genomas, pedido.k, pedido.cupo)?;
    let contexto = Contexto {
        k: pedido.k,
        cupo: pedido.cupo,
        cifras: cifras_para(&comparacion),
    };
    let mut narracion = Narracion::nueva(400);
    let resumen: Vec<String> = pedido
        .nombres
        .iter()
        .enumerate()
        .map(|(indice, nombre)| format!("{nombre}: {} bases", miles(pedido.bases(indice))))
        .collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("genomas.len()", pedido.genomas.len())
            .variable("k", pedido.k)
            .variable("cupo", miles(pedido.cupo))
            .lista("Secuencias", resumen, &[])
            .paso(
                L_COMPARAR,
                format!(
                    "Vas a comparar {} secuencias sin alinearlas. Cada una se resume en un bosquejo: los {} hashes más chicos de sus k-meros de largo {}. Después se comparan los bosquejos, no las secuencias.",
                    pedido.genomas.len(),
                    miles(pedido.cupo),
                    pedido.k
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", pedido.k)
            .variable("cupo", miles(pedido.cupo))
            .paso(
                L_VALIDAR,
                format!(
                    "validar revisa que k esté entre 1 y {} (cada k-mero cabe en 64 bits, 2 por base), que el cupo esté entre 1 y {} y que entre todos los bosquejos no pasen de {} hashes.",
                    minhash::K_MAXIMO,
                    miles(minhash::CUPO_MAXIMO),
                    miles(minhash::HASHES_MAXIMOS_EN_TOTAL)
                ),
            ),
    );
    let hashes_del_primero = narrar_primer_bosquejo(
        &mut narracion,
        &pedido.nombres[0],
        &pedido.genomas[0],
        &contexto,
    );
    narracion.agregar(
        Pizarra::nueva()
            .tabla(
                "Bosquejos",
                &["secuencia", "bases", "hashes", "mayor hash"],
                filas_de_bosquejos(pedido, &comparacion, contexto.cifras),
                &[],
            )
            .paso(
                L_BOSQUEJAR_TODOS,
                format!(
                    "Lo mismo con las otras {} secuencias, en paralelo (par_iter reparte una secuencia por núcleo). Cada bosquejo pesa {} KB aunque la secuencia tenga millones de bases.",
                    pedido.genomas.len() - 1,
                    miles((pedido.cupo * 8).div_ceil(1024))
                ),
            ),
    );
    let (primero, segundo) = par_para_mezclar(&comparacion);
    narracion.agregar(
        Pizarra::nueva()
            .variable("pares", miles(comparacion.pares.len()))
            .paso(
                L_PARES,
                format!(
                    "Con {} secuencias hay {} pares. Te muestro uno de cerca: «{}» con «{}», el que comparte más hashes sin ser idéntico.",
                    pedido.genomas.len(),
                    miles(comparacion.pares.len()),
                    pedido.nombres[primero],
                    pedido.nombres[segundo]
                ),
            ),
    );
    let a = &comparacion.bosquejos[primero];
    let b = &comparacion.bosquejos[segundo];
    let nombres = (
        pedido.nombres[primero].as_str(),
        pedido.nombres[segundo].as_str(),
    );
    let (compartidos, vistos) = narrar_mezcla(&mut narracion, a, b, nombres, &contexto);
    let medida = minhash::medir(a, b);
    if hashes_del_primero != comparacion.bosquejos[0].hashes
        || (compartidos, vistos) != (medida.compartidos, medida.vistos)
    {
        return Err("El paso a paso no coincide con el cálculo del núcleo.".to_string());
    }
    narrar_medidas(&mut narracion, a, b, (compartidos, vistos), nombres);
    narrar_pares(&mut narracion, pedido, &comparacion);
    narracion.agregar(Pizarra::nueva().paso(
        L_CONTENCION,
        "La distancia de Mash supone dos genomas de tamaño parecido. Para saber si una secuencia chica está dentro de una grande se mide la contención.",
    ));
    narrar_contencion(&mut narracion, pedido, &comparacion);
    let ultimo = Pizarra::nueva()
        .variable("pares", miles(comparacion.pares.len()))
        .paso(
            L_FIN,
            "Listo: bosquejos, distancias, ANI, valores p y contenciones. Con la matriz de distancias se arma además el árbol de unión de vecinos. Mira «Qué significa» para leer los números.",
        );
    Ok(narracion.terminar(ultimo))
}

fn pasos_minhash(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let pedido = leer_pedido(entrada)?;
    trazar(&pedido)
}

#[cfg(test)]
mod pruebas {
    use super::*;

    fn ejemplo() -> Entrada {
        Entrada::nueva(&[
            ("secuencias", EJEMPLO),
            ("archivo", ""),
            ("k", "21"),
            ("cupo", "1000"),
        ])
    }

    #[test]
    fn el_ejemplo_lee_referencias_y_fasta() {
        let pedido = leer_pedido(&ejemplo()).unwrap();
        assert_eq!(
            pedido.nombres,
            vec![
                "ecoli:355001-375000",
                "ecoli:365001-385000",
                "lacZ",
                "ecoli:2000001-2020000",
                "lacZ_1",
                "lacZ_5"
            ]
        );
        assert_eq!(pedido.bases(0), 20_000);
        assert_eq!(pedido.bases(4), 3075);
        assert_eq!(pedido.bases(5), 3075);
    }

    #[test]
    fn el_paso_a_paso_da_lo_mismo_que_el_nucleo() {
        let pedido = leer_pedido(&ejemplo()).unwrap();
        let pasos = trazar(&pedido).unwrap();
        assert!(pasos.len() >= 30 && pasos.len() <= 400);
        let codigo: String = crate::codigo("puerto.minhash")
            .unwrap()
            .iter()
            .map(|fragmento| fragmento.fuente.clone())
            .collect::<Vec<_>>()
            .join("\n");
        for paso in &pasos {
            assert!(
                codigo
                    .lines()
                    .any(|linea| linea.contains(paso.linea.as_str())),
                "{}",
                paso.linea
            );
        }
    }

    #[test]
    fn el_paso_a_paso_coincide_con_entradas_al_azar() {
        let mut azar = nucleo::azar::Azar::nuevo(17);
        for _ in 0..25 {
            let genomas_en_total = 2 + azar.entero(4);
            let base: Vec<u8> = (0..300 + azar.entero(6000))
                .map(|_| b"ACGTN"[azar.entero(5).min(azar.entero(60).min(4))])
                .collect();
            let mut texto = String::new();
            for numero in 0..genomas_en_total {
                let copia: Vec<u8> = base
                    .iter()
                    .map(|&letra| {
                        if azar.entero(50) == 0 {
                            b"ACGT"[azar.entero(4)]
                        } else {
                            letra
                        }
                    })
                    .collect();
                texto.push_str(&format!(
                    ">c{numero}\n{}\n",
                    String::from_utf8(copia).unwrap()
                ));
            }
            let k = (5 + azar.entero(28)).to_string();
            let cupo = (1 + azar.entero(1500)).to_string();
            let entrada = Entrada::nueva(&[
                ("secuencias", texto.as_str()),
                ("k", k.as_str()),
                ("cupo", cupo.as_str()),
            ]);
            if let Ok(pedido) = leer_pedido(&entrada) {
                let pasos = trazar(&pedido).unwrap();
                assert!(!pasos.is_empty() && pasos.len() <= 400);
            }
        }
    }

    #[test]
    fn el_ejemplo_muestra_los_ani_esperados() {
        let salida = ejecutar_minhash(&ejemplo()).unwrap();
        let Bloque::Tabla {
            filas, columnas, ..
        } = &salida.bloques[0]
        else {
            panic!("la primera salida debe ser la tabla de pares");
        };
        let ani = columnas.iter().position(|c| c == "ANI estimado").unwrap();
        let buscar = |a: &str, b: &str| {
            filas
                .iter()
                .find(|fila| fila[0] == a && fila[1] == b || fila[0] == b && fila[1] == a)
                .unwrap()
                .clone()
        };
        let valor =
            |fila: Vec<String>| -> f64 { fila[ani].replace("\u{a0}%", "").parse().unwrap() };
        let uno = valor(buscar("lacZ", "lacZ_1"));
        let cinco = valor(buscar("lacZ", "lacZ_5"));
        assert!((uno - 99.0).abs() < 0.3, "{uno}");
        assert!((cinco - 95.0).abs() < 0.8, "{cinco}");
        assert!(salida.interpretacion.len() >= 2 && salida.interpretacion.len() <= 6);
    }

    #[test]
    fn el_genoma_completo_se_narra_rapido() {
        let entrada = Entrada::nueva(&[
            ("secuencias", "@ecoli @ecoli:1-50000"),
            ("k", "21"),
            ("cupo", "1000"),
        ]);
        let pasos = pasos_minhash(&entrada).unwrap();
        assert!(pasos.len() <= 400);
        let salida = ejecutar_minhash(&entrada).unwrap();
        assert!(!salida.interpretacion.is_empty());
    }

    #[test]
    fn entradas_malas_dan_errores_claros() {
        let con = |secuencias: &str, k: &str, cupo: &str| {
            leer_pedido(&Entrada::nueva(&[
                ("secuencias", secuencias),
                ("k", k),
                ("cupo", cupo),
            ]))
        };
        assert!(con("ACGTACGT", "21", "1000").is_err());
        assert!(con("@ecoli:1-100 @ecoli:50-150", "0", "1000").is_err());
        assert!(con("@ecoli:1-100 @ecoli:50-150", "33", "1000").is_err());
        assert!(con("@ecoli:1-100 @ecoli:50-150", "21", "0").is_err());
        assert!(con("@ecoli:1-100 NNNNNNNNNNNNNNNNNNNNNNNNNN", "21", "100").is_err());
        assert!(con("1 2 3", "21", "100").is_err());
        assert!(con("@nada @ecoli:1-100", "21", "100").is_err());
        let bien = con(
            "@ecoli:1-100 @ecoli:50-150, acgtacgtacgtacgtacgtacgt",
            "21",
            "100",
        )
        .unwrap();
        assert_eq!(bien.nombres, vec!["ecoli:1-100", "ecoli:50-150", "s3"]);
        let repetidos = con("@gen:lacZ @gen:lacZ", "21", "100").unwrap();
        assert_eq!(repetidos.nombres, vec!["lacZ", "lacZ_2"]);
    }

    #[test]
    fn lee_archivos_fasta() {
        let carpeta = std::env::temp_dir().join("la-compania-minhash-prueba");
        std::fs::create_dir_all(&carpeta).unwrap();
        let ruta = carpeta.join("dos.fasta");
        std::fs::write(
            &ruta,
            ">uno\nACGTACGGTACCATGACCGTTAGC\n>dos\nTTGACCATGGCATTACGGATCCA\n",
        )
        .unwrap();
        let mala = carpeta.join("lecturas.fastq");
        std::fs::write(&mala, "@r1\nACGT\n+\nIIII\n").unwrap();
        let texto = ruta.display().to_string();
        let pedido = leer_pedido(&Entrada::nueva(&[
            ("secuencias", "@gen:lacZ"),
            ("archivo", texto.as_str()),
            ("k", "15"),
            ("cupo", "100"),
        ]))
        .unwrap();
        assert_eq!(pedido.nombres, vec!["lacZ", "dos"]);
        assert_eq!(pedido.genomas[1].len(), 2);
        let texto_malo = mala.display().to_string();
        assert!(
            leer_pedido(&Entrada::nueva(&[
                ("secuencias", "@gen:lacZ"),
                ("archivo", texto_malo.as_str()),
                ("k", "15"),
                ("cupo", "100"),
            ]))
            .is_err()
        );
    }

    #[test]
    fn un_archivo_que_se_descomprime_de_mas_se_rechaza() {
        let carpeta = std::env::temp_dir().join("la-compania-minhash-tope");
        std::fs::create_dir_all(&carpeta).unwrap();
        let ruta = carpeta.join("largo.fasta");
        std::fs::write(&ruta, format!(">uno\n{}\n", "ACGT".repeat(100))).unwrap();
        assert!(leer_piezas_de_archivo(&ruta, 100).is_err());
        assert_eq!(leer_piezas_de_archivo(&ruta, 10_000).unwrap()[0].len(), 400);
    }

    #[test]
    fn demasiados_hashes_en_total_se_rechazan() {
        let texto = (0..100).map(|_| "@gen:lacZ").collect::<Vec<_>>().join(" ");
        let pedir = |cupo: &str| {
            leer_pedido(&Entrada::nueva(&[
                ("secuencias", texto.as_str()),
                ("k", "21"),
                ("cupo", cupo),
            ]))
        };
        assert!(pedir("100000").is_err());
        assert!(pedir("10000").is_ok());
    }

    #[test]
    fn valor_p_en_texto() {
        assert_eq!(texto_de_valor_p(0.0), "1");
        assert_eq!(texto_de_valor_p(-1.0), "0.1000");
        assert_eq!(texto_de_valor_p(-1043.5), "3.2e-1044");
        assert_eq!(texto_de_valor_p(f64::NEG_INFINITY), "0");
    }
}
