use crate::formato::{cientifico, decimal, miles, porcentaje};
use crate::pizarra::{Marca, Narracion, Pizarra};
use crate::{Algoritmo, Bloque, Campo, Entrada, FilaAlineada, Paso, Salida, Tipo};
use nucleo::busqueda::{
    self, Acierto, BaseDeAdn, BaseDeProteinas, LIMITE_DE_GOLPES_ADN, LONGITUD_MAXIMA_ADN,
    LONGITUD_MAXIMA_BLASTX, LONGITUD_MAXIMA_PROTEINA, Parametros, Tramo,
};
use nucleo::ecoli;
use nucleo::formatos::{Registro, leer_fasta, leer_fasta_archivo};
use nucleo::matrices::AMINOACIDOS;
use nucleo::secuencia::{complemento_inverso, traducir};
use std::collections::BTreeSet;
use std::path::Path;
use std::sync::OnceLock;

const EJEMPLO_BLASTN: &str = "CGTTGAATACGACACCCCGACCCGTCACTACGCACACGTAGACTGCCCGGGGCACGCCGACTATGTTAAAAACATGATCACCGGTGCTGCTCATATGGACGGTGCGATCCTGGTAGTTGCTGCGATTGACGGCCCTATGCCGCAGACTCGTGAGCTCATCCTGCTGGGTCGTCAGGTAGGCCTTCCGTACATCATCGTGTTCCTGAACACATACGACATGGTTGATGACGAAGAGCTGCTTGAACTGGTTGAAATAGAAGTTCGTGAACCTCCGTCTCAGTACGACTTCCCGGGCGGCGACACTCCGATCGTTCGTGGTTCTGCTTTGAAAGCGCTGGAAGGCGACGCAGAGTGGGAAGCGAAAATACTGTACCTGTCTGGCTTCCTTTATTCTTATATTCCGGAACCAGAGCGTGAGATTGACAAGCCGTTCCTGCTGCCGATCGAAGACGTACTCTCCATCTCCGGTCGTGGTACCGTTGTTACCGGTCGTGTAGAAC";

const EJEMPLO_BLASTX: &str = "GCGCTTGGGGCAGGTGGTCTGCCGATGGGCCGTATCGTCGAAATCTACGGACCGGAATCTTCCGGTAAAACCACGCTGACGCTGCAGGTGATCGCCGCAGCGCAGCGTGAAGGTAAAACCTGTGCGTTTATCGATGCTGAACACGCGCTGGACCCAATCTACGCACGTAAACTGGGCGTCGATATCGACAACCTGCTGTGCTCCCAGCCGGACACCGGCGAGCAGGCACTGGAAATCTGTGACGCCCTGGCGCGTTCTGGCGCAGTAGACGTTATCGTCGTTGACTCCGTGGCGGCACTGCGCCGAAAGCGGAAATCGAAGGCGAAATCGGCGACTCTCACATGGGCCTTGCGGCACGTATGATGAGCCAGGCGATGCGTAAGCTGGCGGGTAACCTGAAGCAGTCCAACACGCTGCTGATCTTCATCAACCAGATCCGTATGAAAATTGGTGTGATGTTCGGTAACCCGGAAACCACTACCGGTGGTAACGCGCTGAAATTCTACGCCTCTGTTCGTCTCGACATCCGTCGTATCGGCGCGGTGAAAGAGGGCGAAAACGTGGTGGGTAGCGAAACCCGCGTGAAAGTGGTGAAGAAC";

const EVALOR_TOPE: f64 = 1000.0;
const ACIERTOS_TOPE: usize = 500;
const ALINEAMIENTOS_A_MOSTRAR: usize = 3;
const COLUMNAS_A_MOSTRAR: usize = 2000;
const BARRAS_A_MOSTRAR: usize = 15;
const BYTES_MAXIMOS_DE_ARCHIVO: u64 = 500_000_000;
const SIGNIFICATIVO: f64 = 1e-5;
const CASI_SEGURO: f64 = 1e-50;

const CAMPO_EVALOR: Campo = Campo {
    id: "evalor",
    etiqueta: "Valor E máximo",
    tipo: Tipo::Decimal,
    ejemplo: "0.001",
};

const CAMPO_MAXIMO: Campo = Campo {
    id: "maximo",
    etiqueta: "Máximo de aciertos",
    tipo: Tipo::Entero,
    ejemplo: "20",
};

const BASE_DE_PROTEINAS: Campo = Campo {
    id: "base",
    etiqueta: "Base de datos: @ecoli (proteínas de E. coli incluidas) o la ruta de un FASTA de proteínas",
    tipo: Tipo::Base,
    ejemplo: "@ecoli",
};

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![
        Algoritmo {
            id: "puerto.blastp",
            capitulo: "puerto",
            nombre: "BLASTP: buscar proteínas parecidas",
            resumen: "Busca en una base de proteínas las que se parecen a tu proteína y dice cuánto se puede confiar en cada parecido.",
            campos: &[
                Campo {
                    id: "consulta",
                    etiqueta: "Proteína de consulta",
                    tipo: Tipo::Proteina,
                    ejemplo: "@proteina:parC",
                },
                BASE_DE_PROTEINAS,
                CAMPO_EVALOR,
                CAMPO_MAXIMO,
            ],
            funcion: "busqueda::blastp",
            ejecutar: ejecutar_blastp,
            pasos: Some(pasos_blastp),
        },
        Algoritmo {
            id: "puerto.blastn",
            capitulo: "puerto",
            nombre: "BLASTN: buscar ADN parecido",
            resumen: "Busca tu secuencia de ADN en las dos hebras de un genoma (o de tu propia base) y la ubica con su valor E.",
            campos: &[
                Campo {
                    id: "consulta",
                    etiqueta: "ADN de consulta",
                    tipo: Tipo::Secuencia,
                    ejemplo: EJEMPLO_BLASTN,
                },
                Campo {
                    id: "base",
                    etiqueta: "Base de datos: @ecoli (genoma de E. coli incluido) o la ruta de un FASTA de ADN",
                    tipo: Tipo::Base,
                    ejemplo: "@ecoli",
                },
                CAMPO_EVALOR,
                CAMPO_MAXIMO,
            ],
            funcion: "busqueda::blastn",
            ejecutar: ejecutar_blastn,
            pasos: Some(pasos_blastn),
        },
        Algoritmo {
            id: "puerto.blastx",
            capitulo: "puerto",
            nombre: "BLASTX: ADN traducido contra proteínas",
            resumen: "Traduce tu ADN en los seis marcos de lectura y busca proteínas parecidas a cada traducción.",
            campos: &[
                Campo {
                    id: "consulta",
                    etiqueta: "ADN de consulta (se traduce en 6 marcos)",
                    tipo: Tipo::Secuencia,
                    ejemplo: EJEMPLO_BLASTX,
                },
                BASE_DE_PROTEINAS,
                CAMPO_EVALOR,
                CAMPO_MAXIMO,
            ],
            funcion: "busqueda::blastx",
            ejecutar: ejecutar_blastx,
            pasos: Some(pasos_blastx),
        },
    ]
}

fn proteinas_de_ecoli() -> &'static BaseDeProteinas {
    static BASE: OnceLock<BaseDeProteinas> = OnceLock::new();
    BASE.get_or_init(busqueda::proteinas_de_ecoli)
}

fn genoma_de_ecoli() -> &'static BaseDeAdn {
    static BASE: OnceLock<BaseDeAdn> = OnceLock::new();
    BASE.get_or_init(busqueda::genoma_de_ecoli)
}

pub fn preparar() {
    proteinas_de_ecoli();
    genoma_de_ecoli();
}

#[derive(Clone, Copy, PartialEq)]
enum Programa {
    Blastp,
    Blastn,
    Blastx,
}

impl Programa {
    fn nombre(self) -> &'static str {
        match self {
            Programa::Blastp => "blastp",
            Programa::Blastn => "blastn",
            Programa::Blastx => "blastx",
        }
    }

    fn unidad(self) -> &'static str {
        match self {
            Programa::Blastp => "aminoácidos",
            Programa::Blastn | Programa::Blastx => "nucleótidos",
        }
    }
}

struct Fila<'a> {
    acierto: &'a Acierto,
    sujeto: String,
    descripcion: String,
    corto: String,
}

struct Busqueda<'a> {
    programa: Programa,
    longitud_consulta: usize,
    evalor_maximo: f64,
    en_ecoli: bool,
    secuencias_en_base: usize,
    letras_en_base: usize,
    aviso: Option<String>,
    filas: Vec<Fila<'a>>,
}

fn leer_limites(entrada: &Entrada) -> Result<(f64, usize), String> {
    let evalor = entrada.decimal("evalor")?;
    if evalor <= 0.0 || evalor > EVALOR_TOPE {
        return Err(format!(
            "El valor E máximo debe ser mayor que 0 y no pasar de {}. Lo usual es 0.001.",
            decimal(EVALOR_TOPE, 0)
        ));
    }
    let maximo = entrada.natural("maximo")?;
    if maximo == 0 || maximo > ACIERTOS_TOPE {
        return Err(format!(
            "El máximo de aciertos debe estar entre 1 y {ACIERTOS_TOPE}."
        ));
    }
    Ok((evalor, maximo))
}

fn revisar_longitud(
    longitud: usize,
    minimo: usize,
    maximo: usize,
    programa: Programa,
) -> Result<(), String> {
    if longitud < minimo {
        return Err(format!(
            "La consulta de {} debe tener al menos {minimo} {}.",
            programa.nombre(),
            programa.unidad()
        ));
    }
    if longitud > maximo {
        let consejo = match programa {
            Programa::Blastp => "Usa una proteína más corta o solo el dominio que te interesa.",
            Programa::Blastn | Programa::Blastx => {
                "Usa un tramo más corto, por ejemplo @ecoli:1-5000."
            }
        };
        return Err(format!(
            "La consulta tiene {} letras y {} acepta hasta {} {}. {consejo}",
            miles(longitud),
            programa.nombre(),
            miles(maximo),
            programa.unidad()
        ));
    }
    Ok(())
}

fn leer_proteina(entrada: &Entrada) -> Result<Vec<u8>, String> {
    let proteina = entrada.proteina("consulta")?;
    revisar_longitud(
        proteina.len(),
        3,
        LONGITUD_MAXIMA_PROTEINA,
        Programa::Blastp,
    )?;
    Ok(proteina)
}

fn leer_adn(
    entrada: &Entrada,
    minimo: usize,
    maximo: usize,
    programa: Programa,
) -> Result<Vec<u8>, String> {
    let mut adn = entrada.secuencia("consulta")?;
    for base in adn.iter_mut() {
        if *base == b'U' {
            *base = b'T';
        }
    }
    if let Some(&rara) = adn
        .iter()
        .find(|base| !matches!(base, b'A' | b'C' | b'G' | b'T' | b'N'))
    {
        return Err(format!(
            "La consulta de {} solo admite A, C, G, T (o U) y N, pero contiene «{}».",
            programa.nombre(),
            rara as char
        ));
    }
    revisar_longitud(adn.len(), minimo, maximo, programa)?;
    Ok(adn)
}

fn fraccion_de_adn<'a>(secuencias: impl Iterator<Item = &'a [u8]>) -> f64 {
    let mut letras = 0usize;
    let mut de_adn = 0usize;
    for secuencia in secuencias {
        letras += secuencia.len();
        de_adn += secuencia
            .iter()
            .filter(|letra| matches!(letra, b'A' | b'C' | b'G' | b'T' | b'U' | b'N'))
            .count();
    }
    if letras == 0 {
        return 0.0;
    }
    de_adn as f64 / letras as f64
}

fn registros_de_la_base(entrada: &Entrada) -> Result<Option<Vec<Registro>>, String> {
    let texto = entrada.texto("base")?;
    let limpio = texto.trim();
    if limpio == "@ecoli" {
        return Ok(None);
    }
    let registros = if limpio.starts_with('>') {
        leer_fasta(limpio.as_bytes())
    } else if limpio.starts_with('@') {
        return Err(format!(
            "No conozco la base «{limpio}». Escribe @ecoli para usar E. coli o la ruta de un archivo FASTA."
        ));
    } else {
        leer_archivo_de_base(Path::new(limpio))?
    };
    let registros: Vec<Registro> = registros
        .into_iter()
        .filter(|registro| !registro.secuencia.is_empty())
        .collect();
    if registros.is_empty() {
        return Err("La base de datos no tiene ninguna secuencia.".to_string());
    }
    Ok(Some(registros))
}

fn leer_archivo_de_base(ruta: &Path) -> Result<Vec<Registro>, String> {
    let datos = std::fs::metadata(ruta).map_err(|_| {
        format!(
            "No existe el archivo «{}». Escribe @ecoli o la ruta de un archivo FASTA.",
            ruta.display()
        )
    })?;
    if !datos.is_file() {
        return Err(format!("«{}» no es un archivo.", ruta.display()));
    }
    if datos.len() > BYTES_MAXIMOS_DE_ARCHIVO {
        return Err(format!(
            "El archivo pesa {} bytes y el máximo es {}.",
            miles(datos.len() as usize),
            miles(BYTES_MAXIMOS_DE_ARCHIVO as usize)
        ));
    }
    leer_fasta_archivo(ruta)
}

fn base_de_proteinas_propia(entrada: &Entrada) -> Result<Option<BaseDeProteinas>, String> {
    let Some(registros) = registros_de_la_base(entrada)? else {
        return Ok(None);
    };
    if fraccion_de_adn(
        registros
            .iter()
            .map(|registro| registro.secuencia.as_slice()),
    ) >= 0.9
    {
        return Err(
            "Esa base parece de ADN (casi solo A, C, G y T). Para buscar en ADN usa blastn."
                .to_string(),
        );
    }
    busqueda::indexar_proteinas(registros).map(Some)
}

fn base_de_adn_propia(entrada: &Entrada) -> Result<Option<BaseDeAdn>, String> {
    let Some(registros) = registros_de_la_base(entrada)? else {
        return Ok(None);
    };
    if fraccion_de_adn(
        registros
            .iter()
            .map(|registro| registro.secuencia.as_slice()),
    ) < 0.9
    {
        return Err(
            "Esa base no parece de ADN: tiene muchas letras que no son A, C, G, T o N. Para proteínas usa blastp o blastx."
                .to_string(),
        );
    }
    busqueda::indexar_adn(registros).map(Some)
}

fn partir_nombre(nombre: &str) -> (String, String) {
    let limpio = nombre.trim();
    match limpio.split_once(char::is_whitespace) {
        Some((primero, resto)) => (primero.to_string(), resto.trim().to_string()),
        None => (limpio.to_string(), String::new()),
    }
}

fn marco_en_texto(marco: i8) -> String {
    if marco > 0 {
        format!("+{marco}")
    } else {
        format!("−{}", marco.unsigned_abs())
    }
}

fn filas_de_proteinas<'a>(
    aciertos: &'a [Acierto],
    base: &BaseDeProteinas,
    programa: Programa,
) -> Vec<Fila<'a>> {
    aciertos
        .iter()
        .map(|acierto| {
            let (sujeto, descripcion) = partir_nombre(&base.registros[acierto.sujeto].nombre);
            let corto = if programa == Programa::Blastx {
                format!("{sujeto} ({})", marco_en_texto(acierto.marco))
            } else {
                sujeto.clone()
            };
            Fila {
                acierto,
                sujeto,
                descripcion,
                corto,
            }
        })
        .collect()
}

fn genes_que_cruza(acierto: &Acierto) -> Vec<&'static ecoli::Gen> {
    ecoli::genes()
        .iter()
        .filter(|anotado| {
            anotado.inicio <= acierto.fin_sujeto && anotado.fin > acierto.inicio_sujeto
        })
        .collect()
}

fn describir_genes(genes: &[&ecoli::Gen]) -> String {
    match genes {
        [] => "fuera de genes que codifican proteínas".to_string(),
        [unico] => format!("{} ({}) {}", unico.nombre, unico.locus, unico.producto),
        varios => {
            let nombres: Vec<&str> = varios
                .iter()
                .take(3)
                .map(|anotado| anotado.nombre.as_str())
                .collect();
            if varios.len() > 3 {
                format!("{} y {} genes más", nombres.join(", "), varios.len() - 3)
            } else {
                nombres.join(", ")
            }
        }
    }
}

fn filas_de_adn<'a>(aciertos: &'a [Acierto], base: &BaseDeAdn, en_ecoli: bool) -> Vec<Fila<'a>> {
    aciertos
        .iter()
        .map(|acierto| {
            if en_ecoli {
                let genes = genes_que_cruza(acierto);
                let corto = match genes.first() {
                    Some(anotado) => anotado.nombre.clone(),
                    None => format!("pos. {}", miles(acierto.inicio_sujeto + 1)),
                };
                Fila {
                    acierto,
                    sujeto: ecoli::IDENTIFICADOR.to_string(),
                    descripcion: describir_genes(&genes),
                    corto,
                }
            } else {
                let (sujeto, descripcion) = partir_nombre(&base.nombres[acierto.sujeto]);
                Fila {
                    acierto,
                    corto: sujeto.clone(),
                    sujeto,
                    descripcion,
                }
            }
        })
        .collect()
}

fn parece_adn(secuencia: &[u8]) -> bool {
    fraccion_de_adn(std::iter::once(secuencia)) >= 0.9
}

fn ejecutar_blastp(entrada: &Entrada) -> Result<Salida, String> {
    let consulta = leer_proteina(entrada)?;
    let (evalor_maximo, maximo) = leer_limites(entrada)?;
    let propia = base_de_proteinas_propia(entrada)?;
    let base = propia.as_ref().unwrap_or_else(|| proteinas_de_ecoli());
    let aciertos = busqueda::blastp(&consulta, base, evalor_maximo, maximo)?;
    let aviso = parece_adn(&consulta).then(|| {
        "Tu consulta parece ADN (casi solo A, C, G y T), pero blastp lee cada letra como un aminoácido. Para ADN usa blastn, o blastx si quieres compararlo con proteínas."
            .to_string()
    });
    Ok(armar_salida(&Busqueda {
        programa: Programa::Blastp,
        longitud_consulta: consulta.len(),
        evalor_maximo,
        en_ecoli: propia.is_none(),
        secuencias_en_base: base.registros.len(),
        letras_en_base: base.total_letras,
        aviso,
        filas: filas_de_proteinas(&aciertos, base, Programa::Blastp),
    }))
}

fn ejecutar_blastn(entrada: &Entrada) -> Result<Salida, String> {
    let consulta = leer_adn(entrada, 11, LONGITUD_MAXIMA_ADN, Programa::Blastn)?;
    let (evalor_maximo, maximo) = leer_limites(entrada)?;
    let propia = base_de_adn_propia(entrada)?;
    let base = propia.as_ref().unwrap_or_else(|| genoma_de_ecoli());
    let aciertos = busqueda::blastn(&consulta, base, evalor_maximo, maximo)?;
    Ok(armar_salida(&Busqueda {
        programa: Programa::Blastn,
        longitud_consulta: consulta.len(),
        evalor_maximo,
        en_ecoli: propia.is_none(),
        secuencias_en_base: base.cantidad(),
        letras_en_base: base.total_letras,
        aviso: None,
        filas: filas_de_adn(&aciertos, base, propia.is_none()),
    }))
}

fn ejecutar_blastx(entrada: &Entrada) -> Result<Salida, String> {
    let consulta = leer_adn(entrada, 9, LONGITUD_MAXIMA_BLASTX, Programa::Blastx)?;
    let (evalor_maximo, maximo) = leer_limites(entrada)?;
    let propia = base_de_proteinas_propia(entrada)?;
    let base = propia.as_ref().unwrap_or_else(|| proteinas_de_ecoli());
    let aciertos = busqueda::blastx(&consulta, base, evalor_maximo, maximo)?;
    Ok(armar_salida(&Busqueda {
        programa: Programa::Blastx,
        longitud_consulta: consulta.len(),
        evalor_maximo,
        en_ecoli: propia.is_none(),
        secuencias_en_base: base.registros.len(),
        letras_en_base: base.total_letras,
        aviso: None,
        filas: filas_de_proteinas(&aciertos, base, Programa::Blastx),
    }))
}

fn coordenadas_de_consulta(programa: Programa, acierto: &Acierto) -> (usize, usize) {
    if programa == Programa::Blastx && acierto.marco < 0 {
        (acierto.fin_consulta, acierto.inicio_consulta + 1)
    } else {
        (acierto.inicio_consulta + 1, acierto.fin_consulta)
    }
}

fn coordenadas_de_sujeto(programa: Programa, acierto: &Acierto) -> (usize, usize) {
    if programa == Programa::Blastn && acierto.marco < 0 {
        (acierto.fin_sujeto, acierto.inicio_sujeto + 1)
    } else {
        (acierto.inicio_sujeto + 1, acierto.fin_sujeto)
    }
}

fn identidad(acierto: &Acierto) -> f64 {
    if acierto.longitud == 0 {
        return 0.0;
    }
    acierto.identidades as f64 / acierto.longitud as f64
}

fn cobertura(acierto: &Acierto, longitud_consulta: usize) -> f64 {
    if longitud_consulta == 0 {
        return 0.0;
    }
    (acierto.fin_consulta - acierto.inicio_consulta) as f64 / longitud_consulta as f64
}

fn descripcion_de_la_base(busqueda: &Busqueda) -> String {
    let unidad = match busqueda.programa {
        Programa::Blastn => "pb",
        Programa::Blastp | Programa::Blastx => "aminoácidos",
    };
    let medida = format!(
        "{} secuencias, {} {unidad}",
        miles(busqueda.secuencias_en_base),
        miles(busqueda.letras_en_base)
    );
    match (busqueda.en_ecoli, busqueda.programa) {
        (true, Programa::Blastn) => format!(
            "Genoma de {} ({} pb), las dos hebras",
            ecoli::NOMBRE,
            miles(busqueda.letras_en_base)
        ),
        (true, _) => format!("Proteínas de {}: {medida}", ecoli::NOMBRE),
        (false, _) => format!("Tu base: {medida}"),
    }
}

fn parametros_en_texto(programa: Programa) -> String {
    match programa {
        Programa::Blastn => "Semillas exactas de 11 nt en las dos hebras; premio +2, castigo −3; huecos: abrir 5, alargar 2; estadística de Karlin y Altschul con λ = 0.625, K = 0.41, H = 0.78.".to_string(),
        Programa::Blastp | Programa::Blastx => "Palabras de 3 aminoácidos con vecinas de puntaje ≥ 11 (BLOSUM62); dos golpes a menos de 40 posiciones en la misma diagonal; huecos: abrir 11, alargar 1; estadística de Karlin y Altschul con λ = 0.267, K = 0.041, H = 0.14.".to_string(),
    }
}

fn fila_de_tabla(programa: Programa, fila: &Fila) -> Vec<String> {
    let acierto = fila.acierto;
    let (inicio_consulta, fin_consulta) = coordenadas_de_consulta(programa, acierto);
    let (inicio_sujeto, fin_sujeto) = coordenadas_de_sujeto(programa, acierto);
    let mut celdas = vec![
        fila.sujeto.clone(),
        fila.descripcion.clone(),
        porcentaje(identidad(acierto), 1),
        miles(acierto.longitud),
        miles(acierto.desapareamientos),
        miles(acierto.aperturas),
        miles(inicio_consulta),
        miles(fin_consulta),
        miles(inicio_sujeto),
        miles(fin_sujeto),
    ];
    match programa {
        Programa::Blastn => celdas.push(if acierto.marco < 0 { "−" } else { "+" }.to_string()),
        Programa::Blastx => celdas.push(marco_en_texto(acierto.marco)),
        Programa::Blastp => {}
    }
    celdas.push(cientifico(acierto.evalor));
    celdas.push(decimal(acierto.bits, 1));
    celdas
}

fn bloque_de_tabla(busqueda: &Busqueda) -> Bloque {
    let mut columnas = vec![
        "Sujeto",
        "Descripción",
        "% identidad",
        "Longitud",
        "Desapareamientos",
        "Aperturas de hueco",
        "Inicio consulta",
        "Fin consulta",
        "Inicio sujeto",
        "Fin sujeto",
    ];
    match busqueda.programa {
        Programa::Blastn => columnas.push("Hebra"),
        Programa::Blastx => columnas.push("Marco"),
        Programa::Blastp => {}
    }
    columnas.push("Valor E");
    columnas.push("Bits");
    let filas = busqueda
        .filas
        .iter()
        .map(|fila| fila_de_tabla(busqueda.programa, fila))
        .collect();
    Bloque::tabla(
        "Aciertos (posiciones contando desde 1, como la tabla de BLAST)",
        &columnas,
        filas,
    )
}

fn texto_visible(letras: &[u8]) -> String {
    String::from_utf8_lossy(&letras[..letras.len().min(COLUMNAS_A_MOSTRAR)]).into_owned()
}

fn bloque_de_alineamiento(numero: usize, busqueda: &Busqueda, fila: &Fila) -> Bloque {
    let acierto = fila.acierto;
    let (inicio_consulta, fin_consulta) = coordenadas_de_consulta(busqueda.programa, acierto);
    let (inicio_sujeto, fin_sujeto) = coordenadas_de_sujeto(busqueda.programa, acierto);
    let consulta = match busqueda.programa {
        Programa::Blastx => format!(
            "Consulta traducida, marco {} (nt {}–{})",
            marco_en_texto(acierto.marco),
            miles(inicio_consulta),
            miles(fin_consulta)
        ),
        _ => format!(
            "Consulta {}–{}",
            miles(inicio_consulta),
            miles(fin_consulta)
        ),
    };
    let recorte = if acierto.longitud > COLUMNAS_A_MOSTRAR {
        format!(
            " (primeras {} de {} columnas)",
            miles(COLUMNAS_A_MOSTRAR),
            miles(acierto.longitud)
        )
    } else {
        String::new()
    };
    Bloque::alineamiento(
        &format!(
            "Alineamiento {numero}: {} — {} de identidad, E = {}{recorte}",
            fila.corto,
            porcentaje(identidad(acierto), 1),
            cientifico(acierto.evalor)
        ),
        vec![
            FilaAlineada {
                nombre: consulta,
                secuencia: texto_visible(&acierto.consulta_alineada),
            },
            FilaAlineada {
                nombre: String::new(),
                secuencia: texto_visible(&acierto.linea_media),
            },
            FilaAlineada {
                nombre: format!(
                    "{} {}–{}",
                    fila.sujeto,
                    miles(inicio_sujeto),
                    miles(fin_sujeto)
                ),
                secuencia: texto_visible(&acierto.sujeto_alineado),
            },
        ],
    )
}

fn bloque_de_barras(busqueda: &Busqueda) -> Bloque {
    let mostradas: Vec<&Fila> = busqueda.filas.iter().take(BARRAS_A_MOSTRAR).collect();
    Bloque::barras(
        "Qué tan improbable es cada acierto por azar: −log10(valor E), más alto es más seguro",
        "Acierto",
        "−log10(E)",
        mostradas
            .iter()
            .enumerate()
            .map(|(posicion, fila)| format!("{}. {}", posicion + 1, fila.corto))
            .collect(),
        mostradas
            .iter()
            .map(|fila| fila.acierto.menos_log10_evalor)
            .collect(),
    )
}

fn armar_salida(busqueda: &Busqueda) -> Salida {
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        "Aciertos",
        format!(
            "{} con E ≤ {}",
            miles(busqueda.filas.len()),
            cientifico(busqueda.evalor_maximo)
        ),
    ));
    salida.agregar(Bloque::valor(
        "Base de datos",
        descripcion_de_la_base(busqueda),
    ));
    salida.agregar(Bloque::texto(
        "Parámetros",
        parametros_en_texto(busqueda.programa),
    ));
    if !busqueda.filas.is_empty() {
        salida.agregar(bloque_de_tabla(busqueda));
        for (posicion, fila) in busqueda
            .filas
            .iter()
            .take(ALINEAMIENTOS_A_MOSTRAR)
            .enumerate()
        {
            salida.agregar(bloque_de_alineamiento(posicion + 1, busqueda, fila));
        }
        salida.agregar(bloque_de_barras(busqueda));
    }
    interpretar(&mut salida, busqueda);
    salida
}

fn interpretar(salida: &mut Salida, busqueda: &Busqueda) {
    if let Some(aviso) = &busqueda.aviso {
        salida.explicar(aviso.clone());
    }
    let Some(mejor) = busqueda.filas.first() else {
        salida.explicar(frase_sin_aciertos(busqueda));
        return;
    };
    salida.explicar(frase_del_evalor(mejor));
    if mejor.acierto.evalor <= SIGNIFICATIVO {
        salida.explicar(frase_de_identidad(busqueda, mejor));
        if let Some(frase) = frase_de_familia(busqueda) {
            salida.explicar(frase);
        }
        for frase in frases_del_programa(busqueda, mejor) {
            salida.explicar(frase);
        }
    }
    salida.explicar(format!(
        "El valor E cuenta cuántos aciertos así de buenos esperarías por azar en una base de este tamaño; los bits ({} en el mejor) miden el parecido sin depender del tamaño de la base, por eso sirven para comparar búsquedas distintas.",
        decimal(mejor.acierto.bits, 1)
    ));
}

fn frase_sin_aciertos(busqueda: &Busqueda) -> String {
    let consejo = match busqueda.programa {
        Programa::Blastn => {
            "Si buscas parientes lejanos, prueba blastx: las proteínas cambian más despacio que el ADN porque el código genético es redundante."
        }
        Programa::Blastp | Programa::Blastx => {
            "Puedes subir el umbral (por ejemplo a 10) para ver parecidos débiles, pero trátalos con desconfianza."
        }
    };
    format!(
        "Ningún acierto llega a E ≤ {}: nada en la base se parece a tu consulta más de lo que esperarías por azar. {consejo}",
        cientifico(busqueda.evalor_maximo)
    )
}

fn frase_del_evalor(mejor: &Fila) -> String {
    let evalor = mejor.acierto.evalor;
    let texto_evalor = if evalor == 0.0 {
        "0 (tan pequeño que la computadora lo redondea a cero)".to_string()
    } else {
        cientifico(evalor)
    };
    let inicio = format!(
        "El mejor acierto, {}, tiene E = {texto_evalor}",
        mejor.corto
    );
    if evalor <= CASI_SEGURO {
        format!(
            "{inicio}: con un valor E tan pequeño (≤ 10⁻⁵⁰) la homología es prácticamente segura."
        )
    } else if evalor <= SIGNIFICATIVO {
        format!(
            "{inicio}: es significativo (E ≤ 10⁻⁵), así que probablemente comparte un ancestro con tu consulta."
        )
    } else if evalor <= 1.0 {
        format!(
            "{inicio}: es débil (entre 10⁻⁵ y 1). Puede ser un pariente lejano o una casualidad; antes de creerlo busca más pruebas, como que el alineamiento cubra un dominio conocido o que aparezca también con otra base."
        )
    } else {
        format!(
            "{inicio}: con E > 1 esperarías al menos un acierto así de bueno por puro azar, así que no es prueba de parentesco."
        )
    }
}

fn frase_de_identidad(busqueda: &Busqueda, mejor: &Fila) -> String {
    let acierto = mejor.acierto;
    let parecido = identidad(acierto);
    let cubierto = cobertura(acierto, busqueda.longitud_consulta);
    let texto_identidad = porcentaje(parecido, 1);
    let texto_cobertura = porcentaje(cubierto, 0);
    if parecido >= 0.95 && cubierto >= 0.90 {
        let origen = if busqueda.en_ecoli {
            "de E. coli o de una cepa muy cercana"
        } else {
            "quizá de un organismo muy cercano"
        };
        return format!(
            "Tiene {texto_identidad} de identidad sobre {texto_cobertura} de tu consulta: lo más probable es que sea el mismo gen ({origen})."
        );
    }
    if parecido >= 0.95 {
        return format!(
            "Es casi idéntico ({texto_identidad}), pero solo en {texto_cobertura} de tu consulta: puede ser el mismo gen con tu secuencia más larga, o un tramo que comparten."
        );
    }
    let parcial = if cubierto < 0.5 {
        format!(
            " Ojo: solo {texto_cobertura} de tu consulta queda alineada, así que pueden compartir solo una parte (por ejemplo un dominio)."
        )
    } else {
        String::new()
    };
    if busqueda.programa == Programa::Blastn || parecido >= 0.30 {
        format!(
            "Con {texto_identidad} de identidad es un homólogo: comparten ancestro pero ya divergieron, así que puede ser el mismo gen en otra especie o un parálogo.{parcial}"
        )
    } else {
        format!(
            "Con solo {texto_identidad} de identidad está en la «zona crepuscular» (menos de 30 %): ahí la identidad sola no prueba parentesco, así que fíjate en el valor E y en que el alineamiento sea largo.{parcial}"
        )
    }
}

fn frase_de_familia(busqueda: &Busqueda) -> Option<String> {
    let significativas: Vec<&Fila> = busqueda
        .filas
        .iter()
        .filter(|fila| fila.acierto.evalor <= SIGNIFICATIVO)
        .collect();
    if busqueda.programa == Programa::Blastn {
        if significativas.len() < 2 {
            return None;
        }
        let donde = if busqueda.en_ecoli {
            "el genoma"
        } else {
            "la base"
        };
        return Some(format!(
            "Tu secuencia aparece {} veces en {donde} con E ≤ 10⁻⁵: hay copias repetidas, como genes duplicados, operones de ARN ribosómico o elementos móviles.",
            significativas.len()
        ));
    }
    let distintas: BTreeSet<usize> = significativas
        .iter()
        .map(|fila| fila.acierto.sujeto)
        .collect();
    if distintas.len() < 2 {
        return None;
    }
    let mismo_genoma = if busqueda.en_ecoli {
        " Como todas son de E. coli, son parálogas entre sí: copias de una duplicación antigua que luego divergieron."
    } else {
        " Si vienen del mismo genoma son parálogas; si vienen de especies distintas pueden ser ortólogas."
    };
    Some(format!(
        "Hay {} proteínas distintas con E ≤ 10⁻⁵: es la huella de una familia de genes.{mismo_genoma}",
        distintas.len()
    ))
}

fn frases_del_programa(busqueda: &Busqueda, mejor: &Fila) -> Vec<String> {
    let mut frases = Vec::new();
    let acierto = mejor.acierto;
    match busqueda.programa {
        Programa::Blastp => {}
        Programa::Blastn => {
            if busqueda.en_ecoli {
                frases.push(format!(
                    "Está en las posiciones {}–{} de {} (contando desde 1, como la anotación del genoma), en la hebra {}: {}.",
                    miles(acierto.inicio_sujeto + 1),
                    miles(acierto.fin_sujeto),
                    ecoli::IDENTIFICADOR,
                    if acierto.marco < 0 { "−" } else { "+" },
                    ubicacion_en_genes(acierto)
                ));
            }
            if acierto.marco < 0 {
                frases.push("Hebra − quiere decir que tu secuencia es el complemento inverso de esa región; por eso, en la tabla, el inicio del sujeto es mayor que el fin.".to_string());
            }
        }
        Programa::Blastx => {
            if let Some(frase) = frase_de_cambio_de_marco(busqueda) {
                frases.push(frase);
            } else if acierto.marco < 0 {
                frases.push(format!(
                    "El mejor acierto se lee en el marco {}: la proteína está en la hebra complementaria, al revés respecto a la secuencia que escribiste.",
                    marco_en_texto(acierto.marco)
                ));
            } else {
                frases.push(format!(
                    "El mejor acierto se lee en el marco {}: la proteína está en la misma hebra que escribiste, con el primer codón empezando en la base {}.",
                    marco_en_texto(acierto.marco),
                    acierto.marco
                ));
            }
            let paradas = acierto
                .consulta_alineada
                .iter()
                .filter(|&&letra| letra == b'*')
                .count();
            if paradas > 0 {
                frases.push(format!(
                    "La traducción alineada tiene {paradas} codón(es) de parada en medio: un gen que funciona no debería tenerlos, así que sospecha de un error de secuencia o de un pseudogén."
                ));
            }
        }
    }
    frases
}

fn ubicacion_en_genes(acierto: &Acierto) -> String {
    let genes = genes_que_cruza(acierto);
    if genes.is_empty() {
        "fuera de los genes que codifican proteínas; puede ser una región entre genes o un gen de ARN".to_string()
    } else {
        format!("dentro de {}", describir_genes(&genes))
    }
}

fn frase_de_cambio_de_marco(busqueda: &Busqueda) -> Option<String> {
    for fila in &busqueda.filas {
        if fila.acierto.evalor > SIGNIFICATIVO {
            continue;
        }
        let marcos: BTreeSet<i8> = busqueda
            .filas
            .iter()
            .filter(|otra| {
                otra.acierto.sujeto == fila.acierto.sujeto
                    && otra.acierto.evalor <= SIGNIFICATIVO
                    && otra.acierto.marco.signum() == fila.acierto.marco.signum()
            })
            .map(|otra| otra.acierto.marco)
            .collect();
        if marcos.len() >= 2 {
            let lista: Vec<String> = marcos.iter().map(|&marco| marco_en_texto(marco)).collect();
            return Some(format!(
                "{} aparece en los marcos {}: casi seguro hay un cambio de marco, una base de más o de menos en tu secuencia (por un error de secuenciación o porque es un pseudogén).",
                fila.sujeto,
                lista.join(" y ")
            ));
        }
    }
    None
}

const LETRAS_MAXIMAS_NARRADAS: usize = 50_000_000;
const PALABRAS_POSIBLES: usize = 8000;
const PALABRAS_A_DETALLAR: usize = 3;
const FILAS_EN_PIZARRA: usize = 6;
const ITERACIONES_A_DETALLAR: usize = 3;
const LETRAS_EN_PIZARRA: usize = 28;
const LETRAS_DE_ALINEAMIENTO: usize = 60;
const GOLPES_A_GUARDAR: usize = 2;
const NUNCA: i64 = i64::MIN / 2;
const VENTANA_DE_SEMILLA: usize = 11;

type Rastro = (Vec<Paso>, Vec<Acierto>, Vec<(usize, i32)>);

struct Libreta {
    narracion: Narracion,
    variables: Vec<(String, String)>,
}

impl Libreta {
    fn nueva() -> Libreta {
        Libreta {
            narracion: Narracion::nueva(400),
            variables: Vec::new(),
        }
    }

    fn poner(&mut self, nombre: &str, valor: impl ToString) {
        let valor = valor.to_string();
        match self
            .variables
            .iter_mut()
            .find(|(existente, _)| existente == nombre)
        {
            Some(par) => par.1 = valor,
            None => self.variables.push((nombre.to_string(), valor)),
        }
    }

    fn pizarra(&self) -> Pizarra {
        self.variables
            .iter()
            .fold(Pizarra::nueva(), |pizarra, (nombre, valor)| {
                pizarra.variable(nombre, valor)
            })
    }

    fn agregar(&mut self, pizarra: Pizarra, linea: &str, nota: impl Into<String>) {
        self.narracion.agregar(pizarra.paso(linea, nota));
    }

    fn anotar(&mut self, linea: &str, nota: impl Into<String>) {
        let pizarra = self.pizarra();
        self.agregar(pizarra, linea, nota);
    }

    fn terminar(self, pizarra: Pizarra, linea: &str, nota: impl Into<String>) -> Vec<Paso> {
        self.narracion.terminar(pizarra.paso(linea, nota))
    }
}

struct Datos<'a> {
    consulta: &'a [u8],
    parametros: &'a Parametros,
    espacio: f64,
    minimo: i32,
    evalor_maximo: f64,
    filtro: &'static str,
}

fn abreviar(letras: &[u8]) -> String {
    let inicio: String = String::from_utf8_lossy(letras)
        .chars()
        .take(LETRAS_EN_PIZARRA)
        .collect();
    if letras.len() > LETRAS_EN_PIZARRA {
        format!("{inicio}… ({} letras)", miles(letras.len()))
    } else {
        inicio
    }
}

fn recortar(letras: &[u8]) -> String {
    let visibles: String = String::from_utf8_lossy(letras)
        .chars()
        .take(LETRAS_DE_ALINEAMIENTO)
        .collect();
    if letras.len() > LETRAS_DE_ALINEAMIENTO {
        format!("{visibles}…")
    } else {
        visibles
    }
}

fn revisar_tamano_de_base(letras: usize) -> Result<(), String> {
    if letras > LETRAS_MAXIMAS_NARRADAS {
        return Err(format!(
            "Para el paso a paso la base puede tener como mucho {} letras y la tuya tiene {}. Usa una base más pequeña.",
            miles(LETRAS_MAXIMAS_NARRADAS),
            miles(letras)
        ));
    }
    Ok(())
}

fn base_en_pizarra(en_ecoli: bool, secuencias: usize, letras: usize, unidad: &str) -> String {
    let origen = if en_ecoli { "E. coli" } else { "tu base" };
    format!(
        "{origen}: {} secuencias, {} {unidad}",
        miles(secuencias),
        miles(letras)
    )
}

fn tabla_de_indices() -> [u8; 256] {
    let mut indices = [u8::MAX; 256];
    for (indice, &aminoacido) in AMINOACIDOS.iter().enumerate() {
        indices[aminoacido as usize] = indice as u8;
    }
    indices
}

fn codigo_de_palabra(ventana: &[u8], indices: &[u8; 256]) -> Option<usize> {
    let mut codigo = 0;
    for &letra in ventana {
        let indice = indices[letra as usize];
        if indice == u8::MAX {
            return None;
        }
        codigo = codigo * 20 + indice as usize;
    }
    Some(codigo)
}

fn palabra_en_texto(codigo: usize) -> String {
    [codigo / 400, codigo / 20 % 20, codigo % 20]
        .iter()
        .map(|&indice| AMINOACIDOS[indice] as char)
        .collect()
}

fn palabras_de_la_consulta(consulta: &[u8], umbral: i32) -> Vec<Vec<u32>> {
    let mut por_palabra: Vec<Vec<u32>> = vec![Vec::new(); PALABRAS_POSIBLES];
    for inicio in 0..consulta.len().saturating_sub(2) {
        for vecina in busqueda::palabras_vecinas(&consulta[inicio..inicio + 3], umbral) {
            por_palabra[vecina as usize].push(inicio as u32);
        }
    }
    por_palabra
}

fn narrar_palabras(
    libreta: &mut Libreta,
    consulta: &[u8],
    parametros: &Parametros,
) -> Vec<Vec<u32>> {
    let umbral = parametros.umbral_vecinas;
    let por_palabra = palabras_de_la_consulta(consulta, umbral);
    let total = consulta.len().saturating_sub(2);
    libreta.anotar(
        "let tabla = tabla_de_vecinas(consulta, parametros.umbral_vecinas);",
        format!(
            "No se busca solo tu palabra exacta, también las parecidas. Tu consulta tiene {} palabras de 3 aminoácidos (una por posición) y de cada una se buscan sus «vecinas».",
            miles(total)
        ),
    );
    libreta.poner("vecinas.len()", 0);
    let indices = tabla_de_indices();
    for inicio in 0..total.min(PALABRAS_A_DETALLAR) {
        let palabra = &consulta[inicio..inicio + 3];
        libreta.poner("inicio", inicio);
        libreta.poner("palabra", String::from_utf8_lossy(palabra));
        let pizarra = libreta
            .pizarra()
            .secuencia(consulta, &[(inicio, inicio + 3, Marca::Actual)]);
        libreta.agregar(
            pizarra,
            "let palabra = &consulta[inicio..inicio + LONGITUD_PALABRA_PROTEINA];",
            format!(
                "Toma la palabra de 3 letras que empieza en la posición {inicio} (se cuenta desde 0): {}.",
                String::from_utf8_lossy(palabra)
            ),
        );
        if codigo_de_palabra(palabra, &indices).is_none() {
            libreta.anotar(
                "let Some(letras) = letras_de_palabra(palabra, &indices) else {",
                "Esta palabra tiene una letra que no es uno de los 20 aminoácidos normales, así que no tiene vecinas y se salta.",
            );
            continue;
        }
        let vecinas = busqueda::palabras_vecinas(palabra, umbral);
        let mut puntuadas: Vec<(i32, usize)> = vecinas
            .iter()
            .map(|&codigo| {
                let texto = palabra_en_texto(codigo as usize);
                let puntaje: i32 = texto
                    .bytes()
                    .zip(palabra)
                    .map(|(vecina, &propia)| parametros.puntuacion.puntaje(propia, vecina))
                    .sum();
                (puntaje, codigo as usize)
            })
            .collect();
        puntuadas.sort_by(|uno, otro| otro.0.cmp(&uno.0).then(uno.1.cmp(&otro.1)));
        let filas: Vec<Vec<String>> = puntuadas
            .iter()
            .take(FILAS_EN_PIZARRA)
            .map(|&(puntaje, codigo)| vec![palabra_en_texto(codigo), puntaje.to_string()])
            .collect();
        libreta.poner("vecinas.len()", vecinas.len());
        let propia = filas
            .iter()
            .position(|fila| fila[0].as_bytes() == palabra)
            .map(|fila| vec![(fila, 1)])
            .unwrap_or_default();
        let pizarra = libreta.pizarra().tabla(
            "Las mejores vecinas de esta palabra",
            &["palabra vecina", "puntaje contra la tuya"],
            filas,
            &propia,
        );
        libreta.agregar(
            pizarra,
            "agregar_vecinas(letras, umbral, &matriz, &mut vecinas);",
            format!(
                "Cambia las letras por otras parecidas (según la matriz BLOSUM62) y se queda con las palabras que, comparadas letra a letra con la tuya, suman al menos {umbral}. Para {} salen {}; la tabla muestra las mejores.",
                String::from_utf8_lossy(palabra),
                miles(vecinas.len())
            ),
        );
    }
    let pares: usize = por_palabra.iter().map(Vec::len).sum();
    let distintas = por_palabra
        .iter()
        .filter(|posiciones| !posiciones.is_empty())
        .count();
    libreta.poner("pares.len()", pares);
    libreta.anotar(
        "pares.push((vecina as u32, inicio as u32));",
        format!(
            "Lo mismo se hace con las {} palabras. Cada vecina se anota junto con la posición de tu consulta que la generó: {} pares en total, con {} palabras distintas.",
            miles(total),
            miles(pares),
            miles(distintas)
        ),
    );
    libreta.anotar(
        "Indice::nuevo(&pares, PALABRAS_DE_PROTEINA)",
        format!(
            "Con los pares se arma un índice, como el de un libro: para cada una de las {} palabras posibles guarda en qué posiciones de tu consulta aparece (o aparece una vecina). Así, cuando la base muestre una palabra, se encuentra al instante dónde coincide con tu consulta.",
            miles(PALABRAS_POSIBLES)
        ),
    );
    por_palabra
}

struct Conteo {
    golpes: u64,
    por_sujeto: Vec<u64>,
    veces: Vec<u64>,
}

fn contar_golpes(base: &BaseDeProteinas, por_palabra: &[Vec<u32>], indices: &[u8; 256]) -> Conteo {
    let mut conteo = Conteo {
        golpes: 0,
        por_sujeto: vec![0; base.registros.len()],
        veces: vec![0; PALABRAS_POSIBLES],
    };
    for (numero, registro) in base.registros.iter().enumerate() {
        for ventana in registro.secuencia.windows(3) {
            if let Some(codigo) = codigo_de_palabra(ventana, indices) {
                conteo.veces[codigo] += 1;
                let en_consulta = por_palabra[codigo].len() as u64;
                conteo.por_sujeto[numero] += en_consulta;
                conteo.golpes += en_consulta;
            }
        }
    }
    conteo
}

fn narrar_golpes(libreta: &mut Libreta, por_palabra: &[Vec<u32>], conteo: &Conteo) {
    libreta.poner("golpes", 0);
    libreta.anotar(
        "let mut golpes: u64 = 0;",
        "Ahora se cuentan los «golpes»: las veces que una palabra de la base coincide con una palabra (o vecina) de tu consulta.",
    );
    libreta.poner("golpes", miles(conteo.golpes as usize));
    let mut aportes: Vec<(u64, usize)> = (0..PALABRAS_POSIBLES)
        .map(|palabra| {
            (
                por_palabra[palabra].len() as u64 * conteo.veces[palabra],
                palabra,
            )
        })
        .filter(|&(aporte, _)| aporte > 0)
        .collect();
    aportes.sort_by(|uno, otro| otro.0.cmp(&uno.0).then(uno.1.cmp(&otro.1)));
    let elementos: Vec<String> = aportes
        .iter()
        .take(3)
        .map(|&(aporte, palabra)| {
            format!(
                "{}: {} en la consulta × {} en la base = {}",
                palabra_en_texto(palabra),
                miles(por_palabra[palabra].len()),
                miles(conteo.veces[palabra] as usize),
                miles(aporte as usize)
            )
        })
        .collect();
    let pizarra = libreta
        .pizarra()
        .lista("Las palabras que más golpes dan", elementos, &[0]);
    libreta.agregar(
        pizarra,
        "golpes += en_consulta * veces;",
        format!(
            "Para cada palabra posible multiplica cuántas veces aparece en tu consulta por cuántas veces aparece en la base y lo suma. Resultado: {} golpes.",
            miles(conteo.golpes as usize)
        ),
    );
    libreta.anotar(
        "if golpes > LIMITE_DE_GOLPES {",
        format!(
            "{} golpes no pasan del límite de {}, así que la búsqueda sigue. Con una secuencia repetitiva habría tantos que se detendría y te pediría acortarla.",
            miles(conteo.golpes as usize),
            miles(busqueda::LIMITE_DE_GOLPES as usize)
        ),
    );
}

fn narrar_estadistica(
    libreta: &mut Libreta,
    longitud_consulta: usize,
    longitud_base: usize,
    secuencias: usize,
    parametros: &Parametros,
    evalor_maximo: f64,
) -> (f64, i32) {
    let estadistica = &parametros.estadistica;
    libreta.poner("longitud_consulta", longitud_consulta);
    libreta.poner("longitud_base", longitud_base);
    libreta.poner("secuencias", secuencias);
    libreta.anotar(
        "let espacio = espacio_de_busqueda(",
        format!(
            "Para saber si un parecido puede ser casualidad hay que contar cuántos lugares hay donde podría aparecer por azar: el «espacio de búsqueda», cerca de {} × {}. Pero un alineamiento no cabe pegado a los bordes, así que se corrige.",
            miles(longitud_consulta),
            miles(longitud_base)
        ),
    );
    let ajuste = busqueda::ajuste_de_longitud(
        longitud_consulta as f64,
        longitud_base as f64,
        secuencias as f64,
        estadistica,
    );
    libreta.poner("ajuste", decimal(ajuste, 2));
    libreta.anotar(
        "let ajuste = (estadistica.k * longitud_consulta * longitud_base).ln() / estadistica.h;",
        format!(
            "Un alineamiento real «gasta» unas {} letras en cada secuencia, así que no cuentan como lugares posibles de inicio. Se calcula con ln(K · {} · {}) / H, con un tope para que nunca pase del largo de las secuencias.",
            decimal(ajuste, 2),
            miles(longitud_consulta),
            miles(longitud_base)
        ),
    );
    let espacio = busqueda::espacio_de_busqueda(
        longitud_consulta as f64,
        longitud_base as f64,
        secuencias as f64,
        estadistica,
    );
    libreta.poner("espacio", cientifico(espacio));
    libreta.anotar(
        "(longitud_consulta - ajuste) * (longitud_base - secuencias * ajuste)",
        format!(
            "Resta ese ajuste a cada lado y multiplica: el espacio de búsqueda es {}.",
            cientifico(espacio)
        ),
    );
    let minimo = busqueda::puntaje_necesario(evalor_maximo, espacio, estadistica);
    libreta.poner("evalor_maximo", cientifico(evalor_maximo));
    libreta.poner("minimo", minimo);
    libreta.anotar(
        "let minimo = puntaje_necesario(evalor_maximo, espacio, &parametros.estadistica);",
        format!(
            "Con ese espacio, ¿cuántos puntos hacen falta para que el valor E sea de {} o menos? Se despeja: {minimo} puntos. Todo alineamiento con menos puntos se descarta sin más cuentas.",
            cientifico(evalor_maximo)
        ),
    );
    let cociente = estadistica.k * espacio / evalor_maximo;
    libreta.anotar(
        "((estadistica.k * espacio / evalor_maximo).ln() / estadistica.lambda).ceil() as i32",
        format!(
            "La cuenta: K · espacio / E = {}; su logaritmo natural es {}; dividido entre λ ({}) da {}, y se redondea hacia arriba a {minimo}.",
            cientifico(cociente),
            decimal(cociente.ln(), 2),
            decimal(estadistica.lambda, 3),
            decimal(cociente.ln() / estadistica.lambda, 2)
        ),
    );
    (espacio, minimo)
}

#[derive(Clone)]
struct EventoDeGolpe {
    posicion_consulta: usize,
    posicion_sujeto: usize,
    ultimo_golpe: i64,
    tramo: Option<Tramo>,
}

#[derive(Default)]
struct Recorrido {
    tramos: Vec<Tramo>,
    primeros: Vec<EventoDeGolpe>,
    extensiones: Vec<EventoDeGolpe>,
    mejor: Option<EventoDeGolpe>,
    golpes: u64,
    abiertos: u64,
    solapados: u64,
    cubiertos: u64,
    extendidos: u64,
}

fn recorrer_sujeto(
    consulta: &[u8],
    sujeto: &[u8],
    indices: &[u8; 256],
    por_palabra: &[Vec<u32>],
    parametros: &Parametros,
) -> Recorrido {
    let mut recorrido = Recorrido::default();
    let diagonales = consulta.len() + sujeto.len();
    let mut ultimo_golpe = vec![NUNCA; diagonales];
    let mut fin_extension = vec![0usize; diagonales];
    let ventana = parametros.ventana_dos_golpes as i64;
    let palabra = parametros.longitud_palabra as i64;
    for (posicion_sujeto, letras) in sujeto.windows(3).enumerate() {
        let Some(codigo) = codigo_de_palabra(letras, indices) else {
            continue;
        };
        for &posicion in &por_palabra[codigo] {
            let posicion_consulta = posicion as usize;
            recorrido.golpes += 1;
            let diagonal = posicion_sujeto + consulta.len() - posicion_consulta;
            if posicion_sujeto < fin_extension[diagonal] {
                recorrido.cubiertos += 1;
                continue;
            }
            let anterior = ultimo_golpe[diagonal];
            let distancia = posicion_sujeto as i64 - anterior;
            let mut evento = EventoDeGolpe {
                posicion_consulta,
                posicion_sujeto,
                ultimo_golpe: anterior,
                tramo: None,
            };
            if distancia >= ventana {
                recorrido.abiertos += 1;
                if recorrido.primeros.len() < GOLPES_A_GUARDAR {
                    recorrido.primeros.push(evento);
                }
                ultimo_golpe[diagonal] = posicion_sujeto as i64;
                continue;
            }
            if distancia < palabra {
                recorrido.solapados += 1;
                continue;
            }
            let tramo = busqueda::extender_sin_huecos(
                consulta,
                sujeto,
                posicion_consulta,
                posicion_sujeto,
                parametros.longitud_palabra,
                &parametros.puntuacion,
                parametros.caida_sin_huecos,
            );
            ultimo_golpe[diagonal] = posicion_sujeto as i64;
            fin_extension[diagonal] = tramo.inicio_sujeto + tramo.longitud;
            recorrido.extendidos += 1;
            evento.tramo = Some(tramo);
            if recorrido.extensiones.len() < GOLPES_A_GUARDAR {
                recorrido.extensiones.push(evento.clone());
            }
            let supera = recorrido
                .mejor
                .as_ref()
                .and_then(|mejor| mejor.tramo)
                .is_none_or(|mejor| tramo.puntaje > mejor.puntaje);
            if supera {
                recorrido.mejor = Some(evento);
            }
            if tramo.puntaje >= parametros.disparo_con_huecos {
                recorrido.tramos.push(tramo);
            }
        }
    }
    recorrido
}

fn puntuacion_en_texto(parametros: &Parametros) -> &'static str {
    if parametros.adn {
        "+2 si las bases son iguales, −3 si no"
    } else {
        "la matriz BLOSUM62"
    }
}

struct Lado {
    nombre: &'static str,
    mejor: &'static str,
    largo: &'static str,
    sube: &'static str,
    corta: &'static str,
    termina: &'static str,
}

const DERECHA: Lado = Lado {
    nombre: "derecha",
    mejor: "mejor_derecha",
    largo: "largo_derecha",
    sube: "if suma > mejor_derecha {",
    corta: "else if suma < mejor_derecha - caida {",
    termina: "while derecha_consulta + paso < consulta.len() && derecha_sujeto + paso < sujeto.len() {",
};

const IZQUIERDA: Lado = Lado {
    nombre: "izquierda",
    mejor: "mejor_izquierda",
    largo: "largo_izquierda",
    sube: "if suma > mejor_izquierda {",
    corta: "else if suma < mejor_izquierda - caida {",
    termina: "while paso < posicion_consulta && paso < posicion_sujeto {",
};

struct Cruce<'a> {
    largo_derecha: usize,
    consulta: &'a [u8],
    sujeto: &'a [u8],
    posicion_consulta: usize,
    posicion_sujeto: usize,
    parametros: &'a Parametros,
}

fn narrar_lado(libreta: &mut Libreta, cruce: &Cruce, lado: &Lado) -> (i32, usize) {
    let derecha = lado.nombre == "derecha";
    let palabra = cruce.parametros.longitud_palabra;
    let caida = cruce.parametros.caida_sin_huecos;
    let puntuacion = &cruce.parametros.puntuacion;
    let (limite, desde_consulta, desde_sujeto) = if derecha {
        (
            (cruce.consulta.len() - cruce.posicion_consulta - palabra)
                .min(cruce.sujeto.len() - cruce.posicion_sujeto - palabra),
            cruce.posicion_consulta + palabra,
            cruce.posicion_sujeto + palabra,
        )
    } else {
        (
            cruce.posicion_consulta.min(cruce.posicion_sujeto),
            cruce.posicion_consulta,
            cruce.posicion_sujeto,
        )
    };
    libreta.poner("suma", 0);
    libreta.poner(lado.mejor, 0);
    libreta.poner(lado.largo, 0);
    libreta.poner("paso", 0);
    let mut suma = 0;
    let mut mejor = 0;
    let mut largo = 0;
    let mut paso = 0;
    let mut filas: Vec<Vec<String>> = Vec::new();
    let mut corto = false;
    while paso < limite {
        let (letra_consulta, letra_sujeto) = if derecha {
            (
                cruce.consulta[desde_consulta + paso],
                cruce.sujeto[desde_sujeto + paso],
            )
        } else {
            (
                cruce.consulta[desde_consulta - paso - 1],
                cruce.sujeto[desde_sujeto - paso - 1],
            )
        };
        let puntaje = puntuacion.puntaje(letra_consulta, letra_sujeto);
        suma += puntaje;
        paso += 1;
        let mejora = suma > mejor;
        if mejora {
            mejor = suma;
            largo = paso;
        }
        corto = !mejora && suma < mejor - caida;
        filas.push(vec![
            paso.to_string(),
            (letra_consulta as char).to_string(),
            (letra_sujeto as char).to_string(),
            puntaje.to_string(),
            suma.to_string(),
        ]);
        if paso <= ITERACIONES_A_DETALLAR {
            libreta.poner("suma", suma);
            libreta.poner(lado.mejor, mejor);
            libreta.poner(lado.largo, largo);
            libreta.poner("paso", paso);
            let hallado_desde = if derecha {
                cruce.posicion_sujeto
            } else {
                cruce.posicion_sujeto - largo
            };
            let hallado_hasta = if derecha {
                cruce.posicion_sujeto + palabra + largo
            } else {
                cruce.posicion_sujeto + palabra + cruce.largo_derecha
            };
            let actual = if derecha {
                desde_sujeto + paso - 1
            } else {
                desde_sujeto - paso
            };
            let pizarra = libreta
                .pizarra()
                .secuencia(
                    cruce.sujeto,
                    &[
                        (hallado_desde, hallado_hasta, Marca::Hallado),
                        (actual, actual + 1, Marca::Actual),
                    ],
                )
                .tabla(
                    &format!("Extensión hacia la {}", lado.nombre),
                    &["paso", "consulta", "sujeto", "puntaje", "suma"],
                    filas.clone(),
                    &[(filas.len() - 1, 4)],
                );
            let nota = if mejora {
                format!(
                    "Compara {} con {}: puntaje {puntaje}, y la suma sube a {suma}. Es la mejor hasta ahora, así que guarda {} = {largo}.",
                    letra_consulta as char, letra_sujeto as char, lado.largo
                )
            } else if corto {
                format!(
                    "Compara {} con {}: puntaje {puntaje}. La suma ({suma}) ya quedó más de {caida} puntos por debajo de la mejor ({mejor}): no vale la pena seguir y deja de extender.",
                    letra_consulta as char, letra_sujeto as char
                )
            } else {
                format!(
                    "Compara {} con {}: puntaje {puntaje}, y la suma queda en {suma}. No supera a la mejor ({mejor}), pero tampoco cayó más de {caida} puntos bajo ella, así que sigue por si se recupera.",
                    letra_consulta as char, letra_sujeto as char
                )
            };
            libreta.agregar(pizarra, if mejora { lado.sube } else { lado.corta }, nota);
        }
        if corto {
            break;
        }
    }
    libreta.poner("suma", suma);
    libreta.poner(lado.mejor, mejor);
    libreta.poner(lado.largo, largo);
    libreta.poner("paso", paso);
    let motivo = if corto {
        "porque la suma cayó demasiado bajo la mejor"
    } else {
        "porque una de las dos secuencias se acabó"
    };
    libreta.anotar(
        lado.termina,
        format!(
            "La extensión hacia la {} terminó tras {paso} comparaciones, {motivo}. Se queda con el mejor punto: {largo} letras más, con suma {mejor}.",
            lado.nombre
        ),
    );
    (mejor, largo)
}

fn narrar_extension(
    libreta: &mut Libreta,
    consulta: &[u8],
    sujeto: &[u8],
    posicion_consulta: usize,
    posicion_sujeto: usize,
    parametros: &Parametros,
) -> Tramo {
    let real = busqueda::extender_sin_huecos(
        consulta,
        sujeto,
        posicion_consulta,
        posicion_sujeto,
        parametros.longitud_palabra,
        &parametros.puntuacion,
        parametros.caida_sin_huecos,
    );
    let palabra = parametros.longitud_palabra;
    if posicion_consulta + palabra > consulta.len() || posicion_sujeto + palabra > sujeto.len() {
        return real;
    }
    libreta.poner("posicion_consulta", posicion_consulta);
    libreta.poner("posicion_sujeto", posicion_sujeto);
    let mut puntaje_palabra = 0;
    let mut filas: Vec<Vec<String>> = Vec::new();
    for paso in 0..palabra {
        let letra_consulta = consulta[posicion_consulta + paso];
        let letra_sujeto = sujeto[posicion_sujeto + paso];
        let puntaje = parametros.puntuacion.puntaje(letra_consulta, letra_sujeto);
        puntaje_palabra += puntaje;
        filas.push(vec![
            (paso + 1).to_string(),
            (letra_consulta as char).to_string(),
            (letra_sujeto as char).to_string(),
            puntaje.to_string(),
        ]);
    }
    libreta.poner("puntaje_palabra", puntaje_palabra);
    let pizarra = libreta
        .pizarra()
        .secuencia(
            sujeto,
            &[(posicion_sujeto, posicion_sujeto + palabra, Marca::Hallado)],
        )
        .tabla(
            "La palabra que coincidió, letra contra letra",
            &["letra", "consulta", "sujeto", "puntaje"],
            filas,
            &[],
        );
    libreta.agregar(
        pizarra,
        "puntaje_palabra += puntuacion.puntaje(",
        format!(
            "Primero se califica la palabra que coincidió, letra contra letra, con {}: suma {puntaje_palabra}.",
            puntuacion_en_texto(parametros)
        ),
    );
    let mut cruce = Cruce {
        largo_derecha: 0,
        consulta,
        sujeto,
        posicion_consulta,
        posicion_sujeto,
        parametros,
    };
    let (mejor_derecha, largo_derecha) = narrar_lado(libreta, &cruce, &DERECHA);
    cruce.largo_derecha = largo_derecha;
    let (mejor_izquierda, largo_izquierda) = narrar_lado(libreta, &cruce, &IZQUIERDA);
    let tramo = Tramo {
        inicio_consulta: posicion_consulta - largo_izquierda,
        inicio_sujeto: posicion_sujeto - largo_izquierda,
        longitud: largo_izquierda + palabra + largo_derecha,
        puntaje: puntaje_palabra + mejor_izquierda + mejor_derecha,
    };
    debug_assert_eq!(tramo, real);
    libreta.poner("tramo.puntaje", tramo.puntaje);
    let pizarra = libreta.pizarra().secuencia(
        sujeto,
        &[(
            tramo.inicio_sujeto,
            tramo.inicio_sujeto + tramo.longitud,
            Marca::Hallado,
        )],
    );
    libreta.agregar(
        pizarra,
        "puntaje: puntaje_palabra + mejor_izquierda + mejor_derecha,",
        format!(
            "El tramo suma {puntaje_palabra} (la palabra) + {mejor_izquierda} (izquierda) + {mejor_derecha} (derecha) = {} puntos. Cubre {} letras: de la posición {} a la {} de la consulta y de la {} a la {} del sujeto (se cuenta desde 0, sin incluir el final).",
            tramo.puntaje,
            miles(tramo.longitud),
            tramo.inicio_consulta,
            tramo.inicio_consulta + tramo.longitud,
            tramo.inicio_sujeto,
            tramo.inicio_sujeto + tramo.longitud
        ),
    );
    real
}

fn evalor_en_texto(evalor: f64) -> String {
    if evalor == 0.0 {
        "0 (tan pequeño que la computadora lo redondea a cero)".to_string()
    } else {
        cientifico(evalor)
    }
}

fn descripcion_de_tramo(tramo: &Tramo) -> String {
    format!(
        "consulta {}–{}, sujeto {}–{} · {} puntos",
        tramo.inicio_consulta,
        tramo.inicio_consulta + tramo.longitud,
        tramo.inicio_sujeto,
        tramo.inicio_sujeto + tramo.longitud,
        tramo.puntaje
    )
}

fn narrar_alineamientos(
    libreta: &mut Libreta,
    datos: &Datos,
    sujeto: &[u8],
    numero: usize,
    tramos: Vec<Tramo>,
) -> Vec<(usize, i32)> {
    let parametros = datos.parametros;
    let consulta = datos.consulta;
    let mut ordenados = tramos.clone();
    ordenados.sort_by(|uno, otro| {
        otro.puntaje
            .cmp(&uno.puntaje)
            .then(uno.inicio_consulta.cmp(&otro.inicio_consulta))
            .then(uno.inicio_sujeto.cmp(&otro.inicio_sujeto))
    });
    libreta.poner("tramos", format!("{} tramos", ordenados.len()));
    let elementos: Vec<String> = ordenados
        .iter()
        .take(FILAS_EN_PIZARRA)
        .map(descripcion_de_tramo)
        .collect();
    let pizarra = libreta
        .pizarra()
        .lista("Tramos, del más fuerte al más débil", elementos, &[0]);
    libreta.agregar(
        pizarra,
        "tramos.sort_by(|uno, otro| {",
        "Los tramos con puntaje suficiente se ordenan del más fuerte al más débil: se empieza por el mejor.",
    );
    let Some(&primero) = ordenados.first() else {
        return Vec::new();
    };
    let (semilla_consulta, semilla_sujeto) =
        busqueda::punto_de_partida(consulta, sujeto, &primero, &parametros.puntuacion);
    libreta.poner("semilla_consulta", semilla_consulta);
    libreta.poner("semilla_sujeto", semilla_sujeto);
    let pizarra = libreta.pizarra().secuencia(
        sujeto,
        &[
            (
                primero.inicio_sujeto,
                primero.inicio_sujeto + primero.longitud,
                Marca::Comparado,
            ),
            (semilla_sujeto, semilla_sujeto + 1, Marca::Actual),
        ],
    );
    libreta.agregar(
        pizarra,
        "punto_de_partida(consulta, sujeto, tramo, &parametros.puntuacion)",
        format!(
            "Dentro del tramo se elige el punto de arranque: el centro de la ventana de {} letras que más puntos suma. Es la posición {semilla_consulta} de la consulta contra la {semilla_sujeto} del sujeto.",
            VENTANA_DE_SEMILLA
        ),
    );
    let preliminar = busqueda::alinear_con_huecos(
        consulta,
        sujeto,
        semilla_consulta,
        semilla_sujeto,
        parametros,
        parametros.caida_preliminar,
    );
    libreta.poner("preliminar.puntaje", preliminar.puntaje);
    libreta.anotar(
        "let preliminar = alinear_con_huecos(",
        format!(
            "Desde ese punto se alinea hacia los dos lados permitiendo huecos (insertar o borrar una letra), con una poda estricta: se abandona cualquier camino que quede más de {} puntos bajo el mejor. Es una prueba rápida y da {} puntos.",
            parametros.caida_preliminar, preliminar.puntaje
        ),
    );
    let pasa = preliminar.puntaje >= datos.minimo;
    libreta.poner("minimo", datos.minimo);
    libreta.anotar(
        "if preliminar.puntaje < minimo {",
        if pasa {
            format!(
                "{} puntos son al menos los {} que se necesitan, así que este tramo merece el alineamiento completo.",
                preliminar.puntaje, datos.minimo
            )
        } else {
            format!(
                "{} puntos no llegan a los {} que se necesitan: este tramo se descarta y se sigue con el siguiente.",
                preliminar.puntaje, datos.minimo
            )
        },
    );
    if pasa {
        let nuevo = busqueda::alinear_con_huecos(
            consulta,
            sujeto,
            semilla_consulta,
            semilla_sujeto,
            parametros,
            parametros.caida_con_huecos,
        );
        libreta.poner("nuevo.puntaje", nuevo.puntaje);
        let huecos = nuevo
            .operaciones
            .iter()
            .filter(|&&operacion| operacion != busqueda::PAREJA)
            .count();
        let pizarra = libreta.pizarra().secuencia(
            sujeto,
            &[
                (
                    nuevo.inicio_sujeto,
                    nuevo.fin_sujeto.max(nuevo.inicio_sujeto),
                    Marca::Hallado,
                ),
                (semilla_sujeto, semilla_sujeto + 1, Marca::Actual),
            ],
        );
        libreta.agregar(
            pizarra,
            "let nuevo = alinear_con_huecos(",
            format!(
                "Ahora se repite con una poda más amplia ({} puntos de caída) para que el alineamiento pueda crecer todo lo que valga la pena. Resultado: {} puntos, {} columnas y {huecos} hueco(s).",
                parametros.caida_con_huecos,
                nuevo.puntaje,
                miles(nuevo.operaciones.len())
            ),
        );
        libreta.anotar(
            "let izquierda = extender_con_huecos(",
            format!(
                "Hacia la izquierda del punto de arranque se alinearon {} letras de la consulta contra {} del sujeto. En cada casilla de la tabla se elige la mejor de tres jugadas: emparejar las dos letras, o abrir o alargar un hueco.",
                semilla_consulta - nuevo.inicio_consulta,
                semilla_sujeto - nuevo.inicio_sujeto
            ),
        );
        libreta.anotar(
            "let derecha = extender_con_huecos(",
            format!(
                "Hacia la derecha: {} letras de la consulta contra {} del sujeto. Un hueco cuesta {} puntos al abrirse y {} por cada letra extra; las casillas que caen demasiado bajo el mejor puntaje se podan.",
                nuevo.fin_consulta - semilla_consulta,
                nuevo.fin_sujeto - semilla_sujeto,
                parametros.apertura + parametros.extension,
                parametros.extension
            ),
        );
        libreta.anotar(
            "alineamientos.push(nuevo);",
            "Si no comparte un extremo con uno ya hecho, el alineamiento se guarda.",
        );
    }
    let hechos = busqueda::alinear_tramos(consulta, sujeto, tramos, parametros, datos.minimo);
    libreta.poner("alineamientos", format!("{} alineamientos", hechos.len()));
    let puntajes: Vec<String> = hechos
        .iter()
        .take(FILAS_EN_PIZARRA)
        .map(|hecho| format!("{} puntos", hecho.puntaje))
        .collect();
    let pizarra = libreta
        .pizarra()
        .lista("Alineamientos de esta secuencia", puntajes, &[]);
    libreta.agregar(
        pizarra,
        "for tramo in &tramos {",
        format!(
            "Los demás tramos se tratan igual, uno por uno; los que caen dentro de un alineamiento ya hecho se saltan. En total quedan {} alineamiento(s) en esta secuencia.",
            hechos.len()
        ),
    );
    let Some(mejor) = hechos.iter().max_by_key(|hecho| hecho.puntaje) else {
        return Vec::new();
    };
    let acierto =
        busqueda::armar_acierto(consulta, sujeto, numero, mejor, parametros, datos.espacio);
    libreta.poner("alineamiento.puntaje", mejor.puntaje);
    let bits = busqueda::bits(mejor.puntaje, &parametros.estadistica);
    libreta.poner("bits_del_acierto", decimal(bits, 1));
    libreta.anotar(
        "let bits_del_acierto = bits(alineamiento.puntaje, &parametros.estadistica);",
        format!(
            "Los puntos dependen de la matriz que se use; los bits no. Se pasan con bits = (λ · puntos − ln K) / ln 2 = ({} · {} − ln {}) / ln 2 = {}.",
            decimal(parametros.estadistica.lambda, 3),
            mejor.puntaje,
            decimal(parametros.estadistica.k, 3),
            decimal(bits, 1)
        ),
    );
    libreta.poner("evalor", cientifico(acierto.evalor));
    libreta.anotar(
        "evalor: evalor(bits_del_acierto, espacio),",
        format!(
            "El valor E es cuántos alineamientos así de buenos esperarías por puro azar: espacio × 2^(−bits) = {} × 2^(−{}) = {}.",
            cientifico(datos.espacio),
            decimal(bits, 1),
            evalor_en_texto(acierto.evalor)
        ),
    );
    libreta.poner("acierto.identidades", acierto.identidades);
    libreta.poner("acierto.desapareamientos", acierto.desapareamientos);
    libreta.poner("acierto.huecos", acierto.huecos);
    let lineas = vec![
        recortar(&acierto.consulta_alineada),
        recortar(&acierto.linea_media),
        recortar(&acierto.sujeto_alineado),
    ];
    let pizarra = libreta.pizarra().lista(
        "El alineamiento (consulta, coincidencias, sujeto)",
        lineas,
        &[1],
    );
    libreta.agregar(
        pizarra,
        "if de_consulta == de_sujeto && puntaje > 0 {",
        format!(
            "Se recorre el alineamiento columna por columna: {} columnas con letras iguales (identidades), {} distintas y {} con hueco. Eso da {} de identidad.",
            miles(acierto.identidades),
            miles(acierto.desapareamientos),
            miles(acierto.huecos),
            porcentaje(identidad(&acierto), 1)
        ),
    );
    let se_queda = acierto.evalor <= datos.evalor_maximo;
    libreta.anotar(
        datos.filtro,
        format!(
            "E = {} {} {}: {}.",
            evalor_en_texto(acierto.evalor),
            if se_queda { "≤" } else { ">" },
            cientifico(datos.evalor_maximo),
            if se_queda {
                "el acierto se queda en la lista"
            } else {
                "el acierto se descarta"
            }
        ),
    );
    hechos
        .iter()
        .map(|hecho| {
            busqueda::armar_acierto(consulta, sujeto, numero, hecho, parametros, datos.espacio)
        })
        .filter(|acierto| acierto.evalor <= datos.evalor_maximo)
        .map(|acierto| (numero, acierto.puntaje))
        .collect()
}

fn distancia_en_texto(golpe: &EventoDeGolpe) -> String {
    if golpe.ultimo_golpe == NUNCA {
        "ninguno".to_string()
    } else {
        golpe.ultimo_golpe.to_string()
    }
}

fn narrar_golpe(
    libreta: &mut Libreta,
    datos: &Datos,
    sujeto: &[u8],
    golpe: &EventoDeGolpe,
    es_primero: bool,
) {
    let parametros = datos.parametros;
    let palabra = parametros.longitud_palabra;
    let en_sujeto =
        String::from_utf8_lossy(&sujeto[golpe.posicion_sujeto..golpe.posicion_sujeto + palabra])
            .into_owned();
    let en_consulta = String::from_utf8_lossy(
        &datos.consulta[golpe.posicion_consulta..golpe.posicion_consulta + palabra],
    )
    .into_owned();
    libreta.poner("posicion_sujeto", golpe.posicion_sujeto);
    libreta.poner("posicion_consulta", golpe.posicion_consulta);
    let pizarra = libreta.pizarra().secuencia(
        sujeto,
        &[(
            golpe.posicion_sujeto,
            golpe.posicion_sujeto + palabra,
            Marca::Actual,
        )],
    );
    libreta.agregar(
        pizarra,
        "for &posicion_consulta in tabla.posiciones_de(codigo as u32) {",
        format!(
            "En la posición {} de la proteína está {en_sujeto}, que coincide (como palabra vecina) con {en_consulta}, la palabra de tu consulta que empieza en la posición {}. Eso es un golpe. Su diagonal es posición en la base − posición en la consulta = {}.",
            golpe.posicion_sujeto,
            golpe.posicion_consulta,
            golpe.posicion_sujeto as i64 - golpe.posicion_consulta as i64
        ),
    );
    libreta.poner("diagonal.ultimo_golpe", distancia_en_texto(golpe));
    let ventana = parametros.ventana_dos_golpes;
    if golpe.ultimo_golpe == NUNCA {
        libreta.poner("distancia", "sin golpe anterior");
    } else {
        libreta.poner(
            "distancia",
            golpe.posicion_sujeto as i64 - golpe.ultimo_golpe,
        );
    }
    let nota = if es_primero {
        if golpe.ultimo_golpe == NUNCA {
            format!(
                "Es el primer golpe de esa diagonal. Solo se anota (ultimo_golpe = {}): un golpe aislado no basta, hace falta otro cerca.",
                golpe.posicion_sujeto
            )
        } else {
            format!(
                "El golpe anterior de esa diagonal quedó a {} posiciones, igual o más lejos que la ventana de {ventana}. Se trata como uno nuevo: solo se anota (ultimo_golpe = {}).",
                golpe.posicion_sujeto as i64 - golpe.ultimo_golpe,
                golpe.posicion_sujeto
            )
        }
    } else {
        format!(
            "Ya había otro golpe en esa diagonal, a {} posiciones (menos que la ventana de {ventana}) y sin solaparse: son dos golpes cercanos, señal de un parecido real. Ahora sí se extiende.",
            golpe.posicion_sujeto as i64 - golpe.ultimo_golpe
        )
    };
    libreta.anotar("if distancia >= ventana {", nota);
}

fn narrar_decision_de_tramo(libreta: &mut Libreta, datos: &Datos, tramo: &Tramo) {
    libreta.poner("tramo.puntaje", tramo.puntaje);
    let disparo = datos.parametros.disparo_con_huecos;
    libreta.anotar(
        "if tramo.puntaje >= parametros.disparo_con_huecos {",
        if tramo.puntaje >= disparo {
            format!(
                "{} puntos llegan al mínimo de {disparo}: el tramo se guarda para alinearlo con huecos más tarde.",
                tramo.puntaje
            )
        } else {
            format!(
                "{} puntos no llegan a los {disparo} necesarios: el tramo se descarta.",
                tramo.puntaje
            )
        },
    );
}

fn narrar_sujeto_de_proteina(
    libreta: &mut Libreta,
    datos: &Datos,
    base: &BaseDeProteinas,
    numero: usize,
    por_palabra: &[Vec<u32>],
) -> Vec<(usize, i32)> {
    let sujeto = &base.registros[numero].secuencia;
    let (nombre, _) = partir_nombre(&base.registros[numero].nombre);
    let parametros = datos.parametros;
    libreta.poner("numero", numero);
    libreta.anotar(
        "let aciertos = (0..base.registros.len())",
        format!(
            "Ahora se revisa cada una de las {} proteínas de la base (la computadora las reparte entre sus núcleos para ir más rápido). Veamos con detalle una: {nombre}.",
            miles(base.registros.len())
        ),
    );
    libreta.poner("sujeto", abreviar(sujeto));
    let pizarra = libreta.pizarra().secuencia(sujeto, &[]);
    libreta.agregar(
        pizarra,
        "let sujeto = &base.registros[numero].secuencia;",
        format!(
            "Es la proteína número {numero}, con {} aminoácidos.",
            miles(sujeto.len())
        ),
    );
    let indices = tabla_de_indices();
    let recorrido = recorrer_sujeto(datos.consulta, sujeto, &indices, por_palabra, parametros);
    libreta.poner("tramos", format!("{} tramos", recorrido.tramos.len()));
    if recorrido.golpes == 0 {
        libreta.anotar(
            "let tramos = tramos_de_dos_golpes(",
            "Se recorre la proteína palabra por palabra, pero ninguna coincide con las de tu consulta: sin golpes no hay nada que extender.",
        );
        return Vec::new();
    }
    libreta.anotar(
        "let tramos = tramos_de_dos_golpes(",
        format!(
            "Se recorre la proteína palabra por palabra y se anotan los golpes. Aquí hay {}: {} abrieron una diagonal, {} tuvieron un segundo golpe cercano y se extendieron, y el resto cayó dentro de un tramo ya extendido o solapado.",
            miles(recorrido.golpes as usize),
            miles(recorrido.abiertos as usize),
            miles(recorrido.extendidos as usize)
        ),
    );
    for golpe in &recorrido.primeros {
        narrar_golpe(libreta, datos, sujeto, golpe, true);
    }
    let mejor = recorrido.mejor.clone();
    let mismo = |uno: &EventoDeGolpe, otro: &EventoDeGolpe| {
        uno.posicion_consulta == otro.posicion_consulta
            && uno.posicion_sujeto == otro.posicion_sujeto
    };
    if let Some(primera) = recorrido.extensiones.first() {
        narrar_golpe(libreta, datos, sujeto, primera, false);
        let es_la_mejor = mejor.as_ref().is_some_and(|mejor| mismo(primera, mejor));
        if es_la_mejor {
            if let Some(tramo) = primera.tramo {
                narrar_extension(
                    libreta,
                    datos.consulta,
                    sujeto,
                    primera.posicion_consulta,
                    primera.posicion_sujeto,
                    parametros,
                );
                narrar_decision_de_tramo(libreta, datos, &tramo);
            }
        } else {
            if let Some(tramo) = primera.tramo {
                libreta.poner("tramo.puntaje", tramo.puntaje);
                libreta.anotar(
                    "let tramo = extender_sin_huecos(",
                    format!(
                        "Extendido sin huecos desde ese golpe queda un tramo de {} letras con {} puntos (más adelante se ve una extensión completa).",
                        miles(tramo.longitud),
                        tramo.puntaje
                    ),
                );
                narrar_decision_de_tramo(libreta, datos, &tramo);
            }
            if let Some(mejor) = &mejor {
                libreta.anotar(
                    "let tramo = extender_sin_huecos(",
                    "Entre todas las extensiones de esta proteína, veamos la que más puntos dio.",
                );
                narrar_golpe(libreta, datos, sujeto, mejor, false);
                if let Some(tramo) = mejor.tramo {
                    narrar_extension(
                        libreta,
                        datos.consulta,
                        sujeto,
                        mejor.posicion_consulta,
                        mejor.posicion_sujeto,
                        parametros,
                    );
                    narrar_decision_de_tramo(libreta, datos, &tramo);
                }
            }
        }
    } else {
        libreta.anotar(
            "if distancia < palabra {",
            "En esta proteína ninguna diagonal llegó a tener dos golpes cercanos, así que no se extiende nada.",
        );
    }
    libreta.anotar(
        "tramos.push(tramo);",
        format!(
            "Al terminar la proteína quedan {} tramo(s) con al menos {} puntos.",
            recorrido.tramos.len(),
            parametros.disparo_con_huecos
        ),
    );
    if recorrido.tramos.is_empty() {
        libreta.anotar(
            "alinear_tramos(consulta, sujeto, tramos, parametros, minimo)",
            "Sin tramos suficientes, de esta proteína no sale ningún alineamiento.",
        );
        return Vec::new();
    }
    narrar_alineamientos(libreta, datos, sujeto, numero, recorrido.tramos)
}

fn elegir_sujeto(aciertos: &[Acierto], por_sujeto: &[u64]) -> Option<usize> {
    let mut copia = aciertos.to_vec();
    busqueda::ordenar_y_recortar(&mut copia, 1);
    if let Some(mejor) = copia.first() {
        return Some(mejor.sujeto);
    }
    por_sujeto
        .iter()
        .enumerate()
        .filter(|(_, golpes)| **golpes > 0)
        .max_by_key(|(numero, golpes)| (**golpes, std::cmp::Reverse(*numero)))
        .map(|(numero, _)| numero)
}

fn narrar_buscar_proteinas(
    libreta: &mut Libreta,
    consulta: &[u8],
    base: &BaseDeProteinas,
    parametros: &Parametros,
    evalor_maximo: f64,
    aciertos: &[Acierto],
) -> Vec<(usize, i32)> {
    libreta.poner("consulta", abreviar(consulta));
    let por_palabra = narrar_palabras(libreta, consulta, parametros);
    let indices = tabla_de_indices();
    let conteo = contar_golpes(base, &por_palabra, &indices);
    narrar_golpes(libreta, &por_palabra, &conteo);
    let (espacio, minimo) = narrar_estadistica(
        libreta,
        consulta.len(),
        base.total_letras,
        base.registros.len(),
        parametros,
        evalor_maximo,
    );
    let datos = Datos {
        consulta,
        parametros,
        espacio,
        minimo,
        evalor_maximo,
        filtro: ".filter(|acierto| acierto.evalor <= evalor_maximo)",
    };
    let detalle = match elegir_sujeto(aciertos, &conteo.por_sujeto) {
        Some(numero) => narrar_sujeto_de_proteina(libreta, &datos, base, numero, &por_palabra),
        None => {
            libreta.anotar(
                "let aciertos = (0..base.registros.len())",
                "Ninguna palabra de la base coincide con las de tu consulta, así que no hay nada que extender ni alinear en ninguna proteína.",
            );
            Vec::new()
        }
    };
    libreta.poner("aciertos", format!("{} aciertos", aciertos.len()));
    libreta.anotar(
        ".flatten()",
        format!(
            "Con todas las proteínas hecho lo mismo, se juntan los resultados: {} acierto(s) con E ≤ {}.",
            miles(aciertos.len()),
            cientifico(evalor_maximo)
        ),
    );
    detalle
}

fn lista_de_mejores(filas: &[Fila]) -> Vec<String> {
    filas
        .iter()
        .take(FILAS_EN_PIZARRA)
        .enumerate()
        .map(|(posicion, fila)| {
            format!(
                "{}. {} · {} bits · E = {}",
                posicion + 1,
                fila.corto,
                decimal(fila.acierto.bits, 1),
                cientifico(fila.acierto.evalor)
            )
        })
        .collect()
}

fn nota_final(filas: &[Fila], evalor_maximo: f64) -> String {
    match filas.first() {
        Some(mejor) => format!(
            "Quedan {} acierto(s). El mejor es {}: {} de identidad y E = {}. Pasa a «Qué significa» para ver qué tan confiable es.",
            miles(filas.len()),
            mejor.corto,
            porcentaje(identidad(mejor.acierto), 1),
            evalor_en_texto(mejor.acierto.evalor)
        ),
        None => format!(
            "Ningún acierto llegó a E ≤ {}. También eso es un resultado: pasa a «Qué significa» para ver qué probar.",
            cientifico(evalor_maximo)
        ),
    }
}

fn rastrear_blastp(entrada: &Entrada) -> Result<Rastro, String> {
    let consulta = leer_proteina(entrada)?;
    let (evalor_maximo, maximo) = leer_limites(entrada)?;
    let propia = base_de_proteinas_propia(entrada)?;
    let base = propia.as_ref().unwrap_or_else(|| proteinas_de_ecoli());
    revisar_tamano_de_base(base.total_letras)?;
    let consulta = consulta.to_ascii_uppercase();
    let parametros = busqueda::parametros_blastp();
    let mut aciertos = busqueda::buscar_proteinas(&consulta, base, &parametros, evalor_maximo)?;
    let todos = aciertos.clone();
    busqueda::ordenar_y_recortar(&mut aciertos, maximo);

    let mut libreta = Libreta::nueva();
    libreta.poner("consulta", abreviar(&consulta));
    libreta.poner("evalor_maximo", cientifico(evalor_maximo));
    libreta.poner("maximo", maximo);
    libreta.poner(
        "base",
        base_en_pizarra(
            propia.is_none(),
            base.registros.len(),
            base.total_letras,
            "aminoácidos",
        ),
    );
    let pizarra = libreta.pizarra().secuencia(&consulta, &[]);
    libreta.agregar(
        pizarra,
        "let consulta = preparar_consulta(",
        format!(
            "Quieres encontrar proteínas parecidas a tu consulta de {} aminoácidos en la base elegida, y quedarte con los aciertos de valor E ≤ {} (hasta {maximo}). Primero se pasa la consulta a mayúsculas y se revisa que tenga entre 3 y {} letras.",
            miles(consulta.len()),
            cientifico(evalor_maximo),
            miles(LONGITUD_MAXIMA_PROTEINA)
        ),
    );
    libreta.poner("parametros.longitud_palabra", parametros.longitud_palabra);
    libreta.poner("parametros.umbral_vecinas", parametros.umbral_vecinas);
    libreta.anotar(
        "let parametros = parametros_blastp();",
        format!(
            "Se cargan los ajustes de BLASTP: palabras de {} aminoácidos, vecinas con puntaje de al menos {} en la matriz BLOSUM62, y la estadística (λ, K, H) para convertir puntos en valor E.",
            parametros.longitud_palabra, parametros.umbral_vecinas
        ),
    );
    libreta.anotar(
        "buscar_proteinas(&consulta, base, &parametros, evalor_maximo)?",
        "Toda la búsqueda ocurre dentro de buscar_proteinas. Vamos a sus pasos.",
    );
    let detalle = narrar_buscar_proteinas(
        &mut libreta,
        &consulta,
        base,
        &parametros,
        evalor_maximo,
        &todos,
    );
    let filas = filas_de_proteinas(&aciertos, base, Programa::Blastp);
    libreta.poner("aciertos", format!("{} aciertos", aciertos.len()));
    let pizarra = libreta
        .pizarra()
        .lista("Los mejores aciertos", lista_de_mejores(&filas), &[0]);
    libreta.agregar(
        pizarra,
        "ordenar_y_recortar(&mut aciertos, maximo);",
        format!(
            "Los aciertos se ordenan del más improbable por azar al menos improbable (−log10 del valor E) y se recorta la lista a {maximo} como máximo."
        ),
    );
    let final_nota = nota_final(&filas, evalor_maximo);
    let pizarra = libreta
        .pizarra()
        .lista("Los mejores aciertos", lista_de_mejores(&filas), &[0]);
    let pasos = libreta.terminar(pizarra, "Ok(aciertos)", final_nota);
    Ok((pasos, aciertos, detalle))
}

fn pasos_blastp(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    rastrear_blastp(entrada).map(|(pasos, _, _)| pasos)
}

struct MarcoTraducido {
    marco: i8,
    desplazamiento: usize,
    proteina: Vec<u8>,
    encontrados: Vec<Acierto>,
}

fn convertir_a_nucleotidos(
    acierto: &mut Acierto,
    marco: i8,
    desplazamiento: usize,
    longitud_consulta: usize,
) {
    let inicio = desplazamiento + 3 * acierto.inicio_consulta;
    let fin = desplazamiento + 3 * acierto.fin_consulta;
    if marco > 0 {
        acierto.inicio_consulta = inicio;
        acierto.fin_consulta = fin;
    } else {
        acierto.inicio_consulta = longitud_consulta - fin;
        acierto.fin_consulta = longitud_consulta - inicio;
    }
    acierto.marco = marco;
}

fn mejor_de(aciertos: &[Acierto]) -> Option<Acierto> {
    let mut copia = aciertos.to_vec();
    busqueda::ordenar_y_recortar(&mut copia, 1);
    copia.into_iter().next()
}

fn rastrear_blastx(entrada: &Entrada) -> Result<Rastro, String> {
    let consulta = leer_adn(entrada, 9, LONGITUD_MAXIMA_BLASTX, Programa::Blastx)?;
    let (evalor_maximo, maximo) = leer_limites(entrada)?;
    let propia = base_de_proteinas_propia(entrada)?;
    let base = propia.as_ref().unwrap_or_else(|| proteinas_de_ecoli());
    revisar_tamano_de_base(base.total_letras)?;
    let consulta = consulta.to_ascii_uppercase();
    let parametros = busqueda::parametros_blastp();
    let inversa = complemento_inverso(&consulta);
    let mut aciertos = Vec::new();
    let mut marcos: Vec<MarcoTraducido> = Vec::new();
    for marco in [1i8, 2, 3, -1, -2, -3] {
        let hebra = if marco > 0 { &consulta } else { &inversa };
        let desplazamiento = marco.unsigned_abs() as usize - 1;
        let proteina = traducir(&hebra[desplazamiento..]);
        let mut encontrados = Vec::new();
        if proteina.len() >= 3 {
            encontrados = busqueda::buscar_proteinas(&proteina, base, &parametros, evalor_maximo)?;
            for acierto in &encontrados {
                let mut convertido = acierto.clone();
                convertir_a_nucleotidos(&mut convertido, marco, desplazamiento, consulta.len());
                aciertos.push(convertido);
            }
        }
        marcos.push(MarcoTraducido {
            marco,
            desplazamiento,
            proteina,
            encontrados,
        });
    }
    busqueda::ordenar_y_recortar(&mut aciertos, maximo);

    let elegido = marcos
        .iter()
        .enumerate()
        .filter_map(|(posicion, marco)| mejor_de(&marco.encontrados).map(|mejor| (posicion, mejor)))
        .max_by(|uno, otro| {
            uno.1
                .menos_log10_evalor
                .total_cmp(&otro.1.menos_log10_evalor)
                .then(otro.0.cmp(&uno.0))
        })
        .map(|(posicion, _)| posicion)
        .or_else(|| marcos.iter().position(|marco| marco.proteina.len() >= 3));

    let mut libreta = Libreta::nueva();
    libreta.poner("consulta", abreviar(&consulta));
    libreta.poner("evalor_maximo", cientifico(evalor_maximo));
    libreta.poner("maximo", maximo);
    libreta.poner(
        "base",
        base_en_pizarra(
            propia.is_none(),
            base.registros.len(),
            base.total_letras,
            "aminoácidos",
        ),
    );
    let pizarra = libreta.pizarra().secuencia(&consulta, &[]);
    libreta.agregar(
        pizarra,
        "let consulta = preparar_consulta(",
        format!(
            "Tienes {} nucleótidos de ADN y quieres encontrar proteínas parecidas a lo que ese ADN podría codificar, con E ≤ {} (hasta {maximo} aciertos). Se pasa a mayúsculas y se revisa que tenga entre 9 y {} nucleótidos.",
            miles(consulta.len()),
            cientifico(evalor_maximo),
            miles(LONGITUD_MAXIMA_BLASTX)
        ),
    );
    libreta.poner("parametros.longitud_palabra", parametros.longitud_palabra);
    libreta.poner("parametros.umbral_vecinas", parametros.umbral_vecinas);
    libreta.anotar(
        "let parametros = parametros_blastp();",
        "BLASTX compara proteínas, así que usa los mismos ajustes que BLASTP: palabras de 3 aminoácidos, vecinas según BLOSUM62.",
    );
    libreta.poner("inversa", abreviar(&inversa));
    let pizarra = libreta.pizarra().secuencia(&inversa, &[]);
    libreta.agregar(
        pizarra,
        "let inversa = complemento_inverso(&consulta);",
        "El ADN tiene dos hebras y el gen puede estar en cualquiera. La otra hebra se lee al revés y cambiando cada base por su complementaria (A↔T, C↔G): es el complemento inverso.",
    );
    let filas: Vec<Vec<String>> = marcos
        .iter()
        .map(|marco| {
            let mejor = mejor_de(&marco.encontrados);
            vec![
                marco_en_texto(marco.marco),
                miles(marco.proteina.len()),
                miles(marco.encontrados.len()),
                mejor
                    .map(|acierto| cientifico(acierto.evalor))
                    .unwrap_or_else(|| "—".to_string()),
            ]
        })
        .collect();
    let resaltar: Vec<(usize, usize)> = elegido.map(|posicion| (posicion, 2)).into_iter().collect();
    let pizarra = libreta.pizarra().tabla(
        "Los seis marcos de lectura",
        &["marco", "aminoácidos", "aciertos", "mejor valor E"],
        filas,
        &resaltar,
    );
    libreta.agregar(
        pizarra,
        "for marco in [1i8, 2, 3, -1, -2, -3] {",
        "Un mismo ADN se puede leer de seis formas: en cada hebra hay tres marcos, según empiece el primer codón en la base 1, 2 o 3. Cada marco se traduce a proteína y se busca por separado, igual que en BLASTP. La tabla resume los seis.",
    );
    let mut detalle = Vec::new();
    if let Some(posicion) = elegido {
        let marco = &marcos[posicion];
        libreta.poner("marco", marco_en_texto(marco.marco));
        libreta.anotar(
            "let hebra = if marco > 0 { &consulta } else { &inversa };",
            format!(
                "Veamos con detalle el marco {}, el del mejor acierto (o el primero si no hubo ninguno). Su signo dice qué hebra se lee: {}.",
                marco_en_texto(marco.marco),
                if marco.marco > 0 {
                    "la que escribiste"
                } else {
                    "el complemento inverso"
                }
            ),
        );
        libreta.poner("desplazamiento", marco.desplazamiento);
        libreta.anotar(
            "let desplazamiento = marco.unsigned_abs() as usize - 1;",
            format!(
                "El número del marco indica en qué base empieza el primer codón: se salta {} base(s).",
                marco.desplazamiento
            ),
        );
        let hebra = if marco.marco > 0 { &consulta } else { &inversa };
        let codones: Vec<String> = hebra[marco.desplazamiento..]
            .as_chunks::<3>()
            .0
            .iter()
            .take(FILAS_EN_PIZARRA)
            .zip(&marco.proteina)
            .map(|(codon, &aminoacido)| {
                format!(
                    "{} → {}",
                    String::from_utf8_lossy(codon),
                    aminoacido as char
                )
            })
            .collect();
        libreta.poner("proteina", abreviar(&marco.proteina));
        let pizarra = libreta
            .pizarra()
            .secuencia(
                hebra,
                &[(
                    marco.desplazamiento,
                    (marco.desplazamiento + 3).min(hebra.len()),
                    Marca::Actual,
                )],
            )
            .lista("Los primeros codones", codones, &[0]);
        libreta.agregar(
            pizarra,
            "let proteina = traducir(&hebra[desplazamiento..]);",
            format!(
                "Se lee el ADN de tres en tres (cada trío es un codón) y cada codón se cambia por su aminoácido con el código genético. Salen {} aminoácidos; un «*» sería un codón de parada.",
                miles(marco.proteina.len())
            ),
        );
        libreta.anotar(
            "if proteina.len() < LONGITUD_PALABRA_PROTEINA {",
            "La proteína tiene al menos 3 aminoácidos, así que se puede buscar. Si fuera más corta, el marco se saltaría.",
        );
        libreta.anotar(
            "for mut acierto in buscar_proteinas(&proteina, base, &parametros, evalor_maximo)? {",
            "Esta proteína traducida se busca igual que en BLASTP: dentro de buscar_proteinas, donde ahora «consulta» es la proteína traducida.",
        );
        detalle = narrar_buscar_proteinas(
            &mut libreta,
            &marco.proteina,
            base,
            &parametros,
            evalor_maximo,
            &marco.encontrados,
        );
        libreta.poner("consulta", abreviar(&consulta));
        if let Some(mut mejor) = mejor_de(&marco.encontrados) {
            libreta.poner("acierto.inicio_consulta", mejor.inicio_consulta);
            libreta.poner("acierto.fin_consulta", mejor.fin_consulta);
            let inicio = marco.desplazamiento + 3 * mejor.inicio_consulta;
            let fin = marco.desplazamiento + 3 * mejor.fin_consulta;
            libreta.poner("inicio", inicio);
            libreta.poner("fin", fin);
            libreta.anotar(
                "let inicio = desplazamiento + 3 * acierto.inicio_consulta;",
                format!(
                    "El mejor acierto de este marco está en los aminoácidos {}–{} de la proteína traducida. Cada aminoácido son 3 nucleótidos, así que en el ADN empieza en {} + 3 × {} = {inicio}.",
                    mejor.inicio_consulta,
                    mejor.fin_consulta,
                    marco.desplazamiento,
                    mejor.inicio_consulta
                ),
            );
            convertir_a_nucleotidos(
                &mut mejor,
                marco.marco,
                marco.desplazamiento,
                consulta.len(),
            );
            let pizarra = libreta.pizarra().secuencia(
                &consulta,
                &[(mejor.inicio_consulta, mejor.fin_consulta, Marca::Hallado)],
            );
            libreta.agregar(
                pizarra,
                if marco.marco > 0 {
                    "acierto.inicio_consulta = inicio;"
                } else {
                    "acierto.inicio_consulta = consulta.len() - fin;"
                },
                if marco.marco > 0 {
                    format!(
                        "En la hebra + las posiciones se quedan tal cual: el acierto cubre los nucleótidos {}–{} de tu ADN.",
                        mejor.inicio_consulta, mejor.fin_consulta
                    )
                } else {
                    format!(
                        "En la hebra − se contó sobre el complemento inverso, así que se invierte: {} − {fin} = {}. El acierto cubre los nucleótidos {}–{} de tu ADN.",
                        consulta.len(),
                        mejor.inicio_consulta,
                        mejor.inicio_consulta,
                        mejor.fin_consulta
                    )
                },
            );
            libreta.poner("acierto.marco", marco_en_texto(mejor.marco));
            libreta.anotar(
                "acierto.marco = marco;",
                format!(
                    "Se anota en qué marco se encontró: {}.",
                    marco_en_texto(mejor.marco)
                ),
            );
        }
    }
    let filas = filas_de_proteinas(&aciertos, base, Programa::Blastx);
    libreta.poner("aciertos", format!("{} aciertos", aciertos.len()));
    let pizarra = libreta
        .pizarra()
        .lista("Los mejores aciertos", lista_de_mejores(&filas), &[0]);
    libreta.agregar(
        pizarra,
        "ordenar_y_recortar(&mut aciertos, maximo);",
        format!(
            "Con los aciertos de los seis marcos juntos, se ordenan del más improbable por azar al menos improbable y se recorta la lista a {maximo} como máximo."
        ),
    );
    let final_nota = nota_final(&filas, evalor_maximo);
    let pizarra = libreta
        .pizarra()
        .lista("Los mejores aciertos", lista_de_mejores(&filas), &[0]);
    let pasos = libreta.terminar(pizarra, "Ok(aciertos)", final_nota);
    Ok((pasos, aciertos, detalle))
}

fn pasos_blastx(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    rastrear_blastx(entrada).map(|(pasos, _, _)| pasos)
}

struct HebraBuscada<'a> {
    hebra: i8,
    secuencia: &'a [u8],
    tramos: Vec<(usize, Tramo)>,
}

fn semilla_del_tramo(
    secuencia: &[u8],
    sujeto: &[u8],
    tramo: &Tramo,
    parametros: &Parametros,
) -> (usize, usize) {
    let palabra = parametros.longitud_palabra;
    let mut primera = None;
    let mut intentos = 0;
    for desplazamiento in 0..=tramo.longitud.saturating_sub(palabra) {
        let posicion_consulta = tramo.inicio_consulta + desplazamiento;
        let posicion_sujeto = tramo.inicio_sujeto + desplazamiento;
        if posicion_consulta + palabra > secuencia.len() || posicion_sujeto + palabra > sujeto.len()
        {
            break;
        }
        let de_consulta = &secuencia[posicion_consulta..posicion_consulta + palabra];
        let de_sujeto = &sujeto[posicion_sujeto..posicion_sujeto + palabra];
        if de_consulta != de_sujeto || !de_consulta.iter().all(|base| b"ACGT".contains(base)) {
            continue;
        }
        primera.get_or_insert((posicion_consulta, posicion_sujeto));
        let extendido = busqueda::extender_sin_huecos(
            secuencia,
            sujeto,
            posicion_consulta,
            posicion_sujeto,
            palabra,
            &parametros.puntuacion,
            parametros.caida_sin_huecos,
        );
        if extendido == *tramo {
            return (posicion_consulta, posicion_sujeto);
        }
        intentos += 1;
        if intentos >= 50 {
            break;
        }
    }
    primera.unwrap_or((tramo.inicio_consulta, tramo.inicio_sujeto))
}

fn inicio_del_registro(base: &BaseDeAdn, numero: usize) -> usize {
    (0..numero)
        .map(|anterior| base.registro(anterior).len() + 1)
        .sum()
}

fn describir_tramo_de_adn(base: &BaseDeAdn, numero: usize, tramo: &Tramo) -> String {
    format!(
        "{} · {}",
        partir_nombre(&base.nombres[numero]).0,
        descripcion_de_tramo(tramo)
    )
}

fn rastrear_blastn(entrada: &Entrada) -> Result<Rastro, String> {
    let consulta = leer_adn(entrada, 11, LONGITUD_MAXIMA_ADN, Programa::Blastn)?;
    let (evalor_maximo, maximo) = leer_limites(entrada)?;
    let propia = base_de_adn_propia(entrada)?;
    let base = propia.as_ref().unwrap_or_else(|| genoma_de_ecoli());
    revisar_tamano_de_base(base.total_letras)?;
    let consulta = consulta.to_ascii_uppercase();
    let parametros = busqueda::parametros_blastn();
    let inversa = complemento_inverso(&consulta);
    let golpes_directa = busqueda::golpes_de_adn(&consulta, base);
    let golpes_inversa = busqueda::golpes_de_adn(&inversa, base);
    let golpes = golpes_directa + golpes_inversa;
    if golpes > LIMITE_DE_GOLPES_ADN {
        return Err(format!(
            "La consulta produce {golpes} coincidencias de palabras con la base y el máximo es {LIMITE_DE_GOLPES_ADN}: suele pasar con secuencias repetitivas o de baja complejidad. Acorta la consulta o quita las repeticiones."
        ));
    }
    let espacio = busqueda::espacio_de_busqueda(
        consulta.len() as f64,
        base.total_letras as f64,
        base.cantidad() as f64,
        &parametros.estadistica,
    );
    let minimo = busqueda::puntaje_necesario(evalor_maximo, espacio, &parametros.estadistica);
    let mut aciertos = Vec::new();
    let mut buscadas: Vec<HebraBuscada> = Vec::new();
    for (hebra, secuencia) in [(1i8, &consulta), (-1i8, &inversa)] {
        let mut tramos = busqueda::tramos_de_semillas_exactas(secuencia, base, &parametros);
        buscadas.push(HebraBuscada {
            hebra,
            secuencia,
            tramos: tramos.clone(),
        });
        tramos.sort_by_key(|&(numero, _)| numero);
        for grupo in tramos.chunk_by(|uno, otro| uno.0 == otro.0) {
            let numero = grupo[0].0;
            let sujeto = base.registro(numero);
            let del_registro: Vec<Tramo> = grupo.iter().map(|&(_, tramo)| tramo).collect();
            for alineamiento in
                busqueda::alinear_tramos(secuencia, sujeto, del_registro, &parametros, minimo)
            {
                let mut acierto = busqueda::armar_acierto(
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
                    busqueda::voltear_a_la_hebra_directa(&mut acierto, consulta.len());
                }
                aciertos.push(acierto);
            }
        }
    }
    busqueda::ordenar_y_recortar(&mut aciertos, maximo);

    let elegido: Option<(i8, usize)> = match aciertos.first() {
        Some(mejor) => Some((mejor.marco, mejor.sujeto)),
        None => buscadas
            .iter()
            .flat_map(|buscada| {
                buscada
                    .tramos
                    .iter()
                    .map(move |(numero, tramo)| (tramo.puntaje, buscada.hebra, *numero))
            })
            .max_by_key(|&(puntaje, hebra, numero)| (puntaje, hebra, std::cmp::Reverse(numero)))
            .map(|(_, hebra, numero)| (hebra, numero)),
    };

    let mut libreta = Libreta::nueva();
    libreta.poner("consulta", abreviar(&consulta));
    libreta.poner("evalor_maximo", cientifico(evalor_maximo));
    libreta.poner("maximo", maximo);
    libreta.poner(
        "base",
        base_en_pizarra(propia.is_none(), base.cantidad(), base.total_letras, "pb"),
    );
    let pizarra = libreta.pizarra().secuencia(&consulta, &[]);
    libreta.agregar(
        pizarra,
        "let consulta = preparar_consulta(",
        format!(
            "Quieres ubicar tu secuencia de {} nucleótidos dentro de la base, en cualquiera de las dos hebras, y quedarte con los aciertos de valor E ≤ {} (hasta {maximo}). Primero se pasa a mayúsculas y se revisa que tenga entre 11 y {} nucleótidos.",
            miles(consulta.len()),
            cientifico(evalor_maximo),
            miles(LONGITUD_MAXIMA_ADN)
        ),
    );
    libreta.poner("parametros.longitud_palabra", parametros.longitud_palabra);
    libreta.anotar(
        "let parametros = parametros_blastn();",
        "Se cargan los ajustes de BLASTN: semillas exactas de 11 letras, +2 por base igual y −3 por base distinta, hueco: 5 para abrirlo y 2 por letra, y la estadística para convertir puntos en valor E.",
    );
    libreta.poner("inversa", abreviar(&inversa));
    let pizarra = libreta.pizarra().secuencia(&inversa, &[]);
    libreta.agregar(
        pizarra,
        "let inversa = complemento_inverso(&consulta);",
        "Tu secuencia puede estar en cualquiera de las dos hebras del ADN. La otra hebra se lee al revés y cambiando cada base por su complementaria (A↔T, C↔G): es el complemento inverso, y también se buscará.",
    );
    let codigos = busqueda::codigos_de_adn(&consulta);
    let validos: Vec<usize> = codigos
        .iter()
        .enumerate()
        .filter(|(_, codigo)| **codigo != u32::MAX)
        .map(|(inicio, _)| inicio)
        .take(PALABRAS_A_DETALLAR)
        .collect();
    let palabra_de_adn = parametros.longitud_palabra;
    if !validos.is_empty() {
        let filas: Vec<Vec<String>> = validos
            .iter()
            .map(|&inicio| {
                vec![
                    inicio.to_string(),
                    String::from_utf8_lossy(&consulta[inicio..inicio + palabra_de_adn])
                        .into_owned(),
                    codigos[inicio].to_string(),
                ]
            })
            .collect();
        let pizarra = libreta.pizarra().tabla(
            "Las primeras palabras de 11 letras",
            &["inicio", "palabra", "número"],
            filas,
            &[],
        );
        libreta.agregar(
            pizarra,
            "codigo = ((codigo << 2) | indice as u32) & mascara;",
            "Cada palabra de 11 bases se guarda como un número: cada base ocupa 2 bits (A=0, C=1, G=2, T=3). Al avanzar una base, el número se desplaza 2 bits y entra la base nueva; así se obtienen todas las palabras de la consulta sin releer cada una.",
        );
        let filas: Vec<Vec<String>> = validos
            .iter()
            .map(|&inicio| {
                let palabra = &consulta[inicio..inicio + palabra_de_adn];
                vec![
                    String::from_utf8_lossy(palabra).into_owned(),
                    miles(busqueda::golpes_de_adn(palabra, base) as usize),
                ]
            })
            .collect();
        let pizarra = libreta.pizarra().tabla(
            "Cuántas veces aparece cada palabra en la base",
            &["palabra", "veces en la base"],
            filas,
            &[],
        );
        libreta.agregar(
            pizarra,
            ".map(|codigo| base.indice.posiciones_de(codigo).len() as u64)",
            "La base tiene un índice, como el de un libro: dado el número de una palabra dice al instante en qué posiciones aparece. Aquí se cuenta cuántas veces aparece cada palabra de tu consulta.",
        );
    }
    libreta.poner("golpes", miles(golpes as usize));
    libreta.anotar(
        "let golpes = golpes_de_adn(&consulta, base) + golpes_de_adn(&inversa, base);",
        format!(
            "Se suman las coincidencias de todas las palabras: {} en tu hebra y {} en la inversa, {} golpes en total.",
            miles(golpes_directa as usize),
            miles(golpes_inversa as usize),
            miles(golpes as usize)
        ),
    );
    libreta.anotar(
        "if golpes > LIMITE_DE_GOLPES_ADN {",
        format!(
            "{} golpes no pasan del límite de {}. Si una secuencia repetitiva produjera millones, la búsqueda se detendría y te pediría acortarla.",
            miles(golpes as usize),
            miles(LIMITE_DE_GOLPES_ADN as usize)
        ),
    );
    narrar_estadistica(
        &mut libreta,
        consulta.len(),
        base.total_letras,
        base.cantidad(),
        &parametros,
        evalor_maximo,
    );
    libreta.poner("aciertos", "[]");
    let mut detalle = Vec::new();
    libreta.anotar(
        "for (hebra, secuencia) in [(1i8, &consulta), (-1i8, &inversa)] {",
        "Se busca en dos pasadas: primero con tu secuencia (hebra +) y luego con su complemento inverso (hebra −).",
    );
    for buscada in &buscadas {
        libreta.poner("hebra", if buscada.hebra > 0 { "+" } else { "−" });
        libreta.poner("tramos", format!("{} tramos", buscada.tramos.len()));
        let mut mejores = buscada.tramos.clone();
        mejores.sort_by(|uno, otro| {
            otro.1
                .puntaje
                .cmp(&uno.1.puntaje)
                .then(uno.1.inicio_consulta.cmp(&otro.1.inicio_consulta))
        });
        let elementos: Vec<String> = mejores
            .iter()
            .take(FILAS_EN_PIZARRA)
            .map(|(numero, tramo)| describir_tramo_de_adn(base, *numero, tramo))
            .collect();
        let pizarra =
            libreta
                .pizarra()
                .lista("Los tramos más fuertes de esta hebra", elementos, &[0]);
        libreta.agregar(
            pizarra,
            "let mut tramos = tramos_de_semillas_exactas(secuencia, base, &parametros);",
            format!(
                "En la hebra {} cada palabra de 11 letras que aparece exacta en la base es una semilla; cada semilla se extiende sin huecos hacia los dos lados. Quedan {} tramo(s) con al menos {} puntos.",
                if buscada.hebra > 0 { "+" } else { "−" },
                miles(buscada.tramos.len()),
                parametros.disparo_con_huecos
            ),
        );
        let Some((hebra_elegida, numero)) = elegido else {
            continue;
        };
        if hebra_elegida != buscada.hebra {
            continue;
        }
        let del_registro: Vec<Tramo> = buscada
            .tramos
            .iter()
            .filter(|(otro, _)| *otro == numero)
            .map(|&(_, tramo)| tramo)
            .collect();
        let Some(mejor_tramo) = del_registro.iter().copied().max_by(|uno, otro| {
            uno.puntaje
                .cmp(&otro.puntaje)
                .then(otro.inicio_consulta.cmp(&uno.inicio_consulta))
                .then(otro.inicio_sujeto.cmp(&uno.inicio_sujeto))
        }) else {
            continue;
        };
        let sujeto = base.registro(numero);
        let (posicion_consulta, posicion_sujeto) =
            semilla_del_tramo(buscada.secuencia, sujeto, &mejor_tramo, &parametros);
        let desfase = inicio_del_registro(base, numero);
        libreta.poner("posicion_consulta", posicion_consulta);
        libreta.poner("posicion", posicion_sujeto + desfase);
        libreta.poner(
            "diagonal",
            (posicion_sujeto + desfase) as i64 - posicion_consulta as i64,
        );
        let pizarra = libreta.pizarra().secuencia(
            sujeto,
            &[(
                posicion_sujeto,
                posicion_sujeto + palabra_de_adn,
                Marca::Actual,
            )],
        );
        libreta.agregar(
            pizarra,
            "let diagonal = posicion as i64 - posicion_consulta as i64;",
            format!(
                "Veamos un tramo con detalle, el de la hebra {}. La palabra que empieza en la posición {posicion_consulta} de la consulta aparece exacta en la base, en la posición {}. La diagonal (posición en la base − posición en la consulta) identifica los alineamientos que avanzan en paralelo: palabras de una misma diagonal ya cubiertas por una extensión no se vuelven a extender.",
                if buscada.hebra > 0 { "+" } else { "−" },
                posicion_sujeto + desfase
            ),
        );
        libreta.anotar(
            "let tramo = extender_sin_huecos(",
            "Desde esa semilla se extiende sin huecos, comparando letra a letra hacia los dos lados.",
        );
        let tramo = narrar_extension(
            &mut libreta,
            buscada.secuencia,
            sujeto,
            posicion_consulta,
            posicion_sujeto,
            &parametros,
        );
        libreta.poner("fin_por_diagonal", tramo.inicio_sujeto + tramo.longitud);
        libreta.anotar(
            "fin_por_diagonal.insert(",
            format!(
                "Se anota hasta dónde llegó esta diagonal (posición {}): las palabras que caigan antes de ese punto ya están cubiertas y se saltan.",
                tramo.inicio_sujeto + tramo.longitud
            ),
        );
        let datos = Datos {
            consulta: buscada.secuencia,
            parametros: &parametros,
            espacio,
            minimo,
            evalor_maximo,
            filtro: "if acierto.evalor > evalor_maximo {",
        };
        narrar_decision_de_tramo(&mut libreta, &datos, &tramo);
        libreta.anotar(
            "tramos.sort_by_key(|&(numero, _)| numero);",
            format!(
                "Los tramos se agrupan por secuencia de la base (aquí {}), para alinear cada grupo por separado.",
                partir_nombre(&base.nombres[numero]).0
            ),
        );
        libreta.poner("numero", numero);
        detalle = narrar_alineamientos(&mut libreta, &datos, sujeto, numero, del_registro);
        if let Some(acierto) = aciertos.first() {
            libreta.poner("acierto.marco", hebra_elegida);
            libreta.anotar(
                "acierto.marco = hebra;",
                format!(
                    "Se anota en qué hebra se encontró: {}.",
                    if hebra_elegida > 0 { "+" } else { "−" }
                ),
            );
            if hebra_elegida < 0 {
                let mut antes = acierto.clone();
                busqueda::voltear_a_la_hebra_directa(&mut antes, consulta.len());
                libreta.anotar(
                    "voltear_a_la_hebra_directa(&mut acierto, consulta.len());",
                    format!(
                        "En la hebra − las posiciones se contaron sobre el complemento inverso. Para que coincidan con tu consulta se invierten: lo que era {}–{} pasa a {}–{} (longitud − fin, longitud − inicio), y las letras del alineamiento también se pasan a la hebra directa.{}",
                        antes.inicio_consulta,
                        antes.fin_consulta,
                        acierto.inicio_consulta,
                        acierto.fin_consulta,
                        if antes.inicio_consulta == acierto.inicio_consulta
                            && antes.fin_consulta == acierto.fin_consulta
                        {
                            " Aquí el alineamiento cubre toda tu consulta, por eso los números no cambian."
                        } else {
                            ""
                        }
                    ),
                );
            }
        }
    }
    let filas = filas_de_adn(&aciertos, base, propia.is_none());
    libreta.poner("aciertos", format!("{} aciertos", aciertos.len()));
    let pizarra = libreta
        .pizarra()
        .lista("Los mejores aciertos", lista_de_mejores(&filas), &[0]);
    libreta.agregar(
        pizarra,
        "ordenar_y_recortar(&mut aciertos, maximo);",
        format!(
            "Con los aciertos de las dos hebras juntos, se ordenan del más improbable por azar al menos improbable y se recorta la lista a {maximo} como máximo."
        ),
    );
    let final_nota = nota_final(&filas, evalor_maximo);
    let pizarra = libreta
        .pizarra()
        .lista("Los mejores aciertos", lista_de_mejores(&filas), &[0]);
    let pasos = libreta.terminar(pizarra, "Ok(aciertos)", final_nota);
    Ok((pasos, aciertos, detalle))
}

fn pasos_blastn(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    rastrear_blastn(entrada).map(|(pasos, _, _)| pasos)
}

#[cfg(test)]
mod pruebas {
    use super::*;

    fn acierto(sujeto: usize, marco: i8) -> Acierto {
        Acierto {
            sujeto,
            marco,
            inicio_consulta: 0,
            fin_consulta: 30,
            inicio_sujeto: 0,
            fin_sujeto: 10,
            puntaje: 100,
            bits: 50.0,
            evalor: 1e-20,
            menos_log10_evalor: 20.0,
            longitud: 10,
            identidades: 10,
            positivos: 10,
            desapareamientos: 0,
            aperturas: 0,
            huecos: 0,
            consulta_alineada: b"AAAAAAAAAA".to_vec(),
            linea_media: b"AAAAAAAAAA".to_vec(),
            sujeto_alineado: b"AAAAAAAAAA".to_vec(),
        }
    }

    fn busqueda_con<'a>(aciertos: &'a [Acierto]) -> Busqueda<'a> {
        Busqueda {
            programa: Programa::Blastx,
            longitud_consulta: 30,
            evalor_maximo: 1e-3,
            en_ecoli: false,
            secuencias_en_base: 1,
            letras_en_base: 10,
            aviso: None,
            filas: aciertos
                .iter()
                .map(|acierto| Fila {
                    acierto,
                    sujeto: "p".to_string(),
                    descripcion: String::new(),
                    corto: "p".to_string(),
                })
                .collect(),
        }
    }

    #[test]
    fn cambio_de_marco_pide_la_misma_hebra() {
        let opuestos = [acierto(0, 1), acierto(0, -2)];
        assert!(frase_de_cambio_de_marco(&busqueda_con(&opuestos)).is_none());
        let misma_hebra = [acierto(0, 1), acierto(0, 3)];
        assert!(frase_de_cambio_de_marco(&busqueda_con(&misma_hebra)).is_some());
        let negativos = [acierto(0, -1), acierto(0, -3)];
        assert!(frase_de_cambio_de_marco(&busqueda_con(&negativos)).is_some());
    }

    fn entrada_de_ejemplo(id: &str) -> Entrada {
        Entrada::desde_json(crate::ejemplo(id).unwrap()).unwrap()
    }

    fn adn_de_prueba(largo: usize, semilla: u64) -> String {
        let mut estado = semilla;
        (0..largo)
            .map(|_| {
                estado = estado
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(estado >> 33) as usize % 4] as char
            })
            .collect()
    }

    fn proteina_de_prueba(largo: usize, semilla: u64) -> String {
        let mut estado = semilla;
        (0..largo)
            .map(|_| {
                estado = estado
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                nucleo::matrices::AMINOACIDOS[(estado >> 33) as usize % 20] as char
            })
            .collect()
    }

    fn verificar_blastp(entrada: &Entrada) -> Rastro {
        let rastro = rastrear_blastp(entrada).unwrap();
        let consulta = leer_proteina(entrada).unwrap();
        let (evalor, maximo) = leer_limites(entrada).unwrap();
        let propia = base_de_proteinas_propia(entrada).unwrap();
        let base = propia.as_ref().unwrap_or_else(|| proteinas_de_ecoli());
        let real = busqueda::blastp(&consulta, base, evalor, maximo).unwrap();
        assert_eq!(rastro.1, real);
        assert!(!rastro.0.is_empty() && rastro.0.len() <= 400);
        if let Some(mejor) = real.first() {
            assert!(rastro.2.contains(&(mejor.sujeto, mejor.puntaje)));
        }
        rastro
    }

    fn verificar_blastx(entrada: &Entrada) -> Rastro {
        let rastro = rastrear_blastx(entrada).unwrap();
        let consulta = leer_adn(entrada, 9, LONGITUD_MAXIMA_BLASTX, Programa::Blastx).unwrap();
        let (evalor, maximo) = leer_limites(entrada).unwrap();
        let propia = base_de_proteinas_propia(entrada).unwrap();
        let base = propia.as_ref().unwrap_or_else(|| proteinas_de_ecoli());
        let real = busqueda::blastx(&consulta, base, evalor, maximo).unwrap();
        assert_eq!(rastro.1, real);
        assert!(!rastro.0.is_empty() && rastro.0.len() <= 400);
        if let Some(mejor) = real.first() {
            assert!(rastro.2.contains(&(mejor.sujeto, mejor.puntaje)));
        }
        rastro
    }

    fn verificar_blastn(entrada: &Entrada) -> Rastro {
        let rastro = rastrear_blastn(entrada).unwrap();
        let consulta = leer_adn(entrada, 11, LONGITUD_MAXIMA_ADN, Programa::Blastn).unwrap();
        let (evalor, maximo) = leer_limites(entrada).unwrap();
        let propia = base_de_adn_propia(entrada).unwrap();
        let base = propia.as_ref().unwrap_or_else(|| genoma_de_ecoli());
        let real = busqueda::blastn(&consulta, base, evalor, maximo).unwrap();
        assert_eq!(rastro.1, real);
        assert!(!rastro.0.is_empty() && rastro.0.len() <= 400);
        if let Some(mejor) = real.first() {
            assert!(rastro.2.contains(&(mejor.sujeto, mejor.puntaje)));
        }
        rastro
    }

    #[test]
    fn el_paso_a_paso_de_blastp_llega_al_mismo_resultado() {
        let rastro = verificar_blastp(&entrada_de_ejemplo("puerto.blastp"));
        assert!(!rastro.1.is_empty());
    }

    #[test]
    fn el_paso_a_paso_de_blastx_llega_al_mismo_resultado() {
        let rastro = verificar_blastx(&entrada_de_ejemplo("puerto.blastx"));
        assert!(!rastro.1.is_empty());
    }

    #[test]
    fn el_paso_a_paso_de_blastn_llega_al_mismo_resultado() {
        let rastro = verificar_blastn(&entrada_de_ejemplo("puerto.blastn"));
        assert!(!rastro.1.is_empty());
    }

    #[test]
    fn blastp_narra_una_base_propia_con_y_sin_aciertos() {
        let familiar = proteina_de_prueba(120, 7);
        let mut variante = familiar.clone().into_bytes();
        for posicion in (5..variante.len()).step_by(9) {
            variante[posicion] = b'A';
        }
        let base = format!(
            ">uno parecida\n{}\n>dos lejana\n{}\n",
            String::from_utf8(variante).unwrap(),
            proteina_de_prueba(150, 99)
        );
        let con_aciertos = Entrada::nueva(&[
            ("consulta", familiar.as_str()),
            ("base", base.as_str()),
            ("evalor", "0.001"),
            ("maximo", "5"),
        ]);
        let rastro = verificar_blastp(&con_aciertos);
        assert!(!rastro.1.is_empty());
        let sin_aciertos = Entrada::nueva(&[
            ("consulta", "WWWWCCCCHHHHWWWW"),
            ("base", base.as_str()),
            ("evalor", "0.001"),
            ("maximo", "5"),
        ]);
        let rastro = verificar_blastp(&sin_aciertos);
        assert!(rastro.1.is_empty());
    }

    #[test]
    fn blastx_narra_una_base_propia_en_la_hebra_negativa() {
        let proteina = proteina_de_prueba(80, 3);
        let codones = |letra: char| match letra {
            'A' => "GCT",
            'R' => "CGT",
            'N' => "AAT",
            'D' => "GAT",
            'C' => "TGT",
            'Q' => "CAA",
            'E' => "GAA",
            'G' => "GGT",
            'H' => "CAT",
            'I' => "ATT",
            'L' => "CTG",
            'K' => "AAA",
            'M' => "ATG",
            'F' => "TTT",
            'P' => "CCG",
            'S' => "TCT",
            'T' => "ACC",
            'W' => "TGG",
            'Y' => "TAT",
            _ => "GTT",
        };
        let adn: String = proteina.chars().map(codones).collect();
        let inversa = String::from_utf8(complemento_inverso(adn.as_bytes())).unwrap();
        let base = format!(
            ">objetivo\n{proteina}\n>otra\n{}\n",
            proteina_de_prueba(90, 41)
        );
        let entrada = Entrada::nueva(&[
            ("consulta", inversa.as_str()),
            ("base", base.as_str()),
            ("evalor", "0.001"),
            ("maximo", "5"),
        ]);
        let rastro = verificar_blastx(&entrada);
        assert!(rastro.1.first().is_some_and(|mejor| mejor.marco < 0));
    }

    #[test]
    fn blastn_narra_una_base_propia_con_dos_registros_y_ambas_hebras() {
        let objetivo = adn_de_prueba(300, 5);
        let inversa = String::from_utf8(complemento_inverso(objetivo.as_bytes())).unwrap();
        let base = format!(
            ">uno\n{}\n>dos\n{}{}{}\n",
            adn_de_prueba(400, 11),
            adn_de_prueba(200, 12),
            inversa,
            adn_de_prueba(200, 13)
        );
        let entrada = Entrada::nueva(&[
            ("consulta", objetivo.as_str()),
            ("base", base.as_str()),
            ("evalor", "0.001"),
            ("maximo", "5"),
        ]);
        let rastro = verificar_blastn(&entrada);
        assert!(rastro.1.first().is_some_and(|mejor| mejor.marco < 0));
        let sin_aciertos = Entrada::nueva(&[
            ("consulta", adn_de_prueba(60, 77).as_str()),
            ("base", base.as_str()),
            ("evalor", "0.001"),
            ("maximo", "5"),
        ]);
        let rastro = verificar_blastn(&sin_aciertos);
        assert!(rastro.1.is_empty());
    }
}
