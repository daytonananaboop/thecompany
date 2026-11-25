use crate::formato::{decimal, miles, porcentaje};
use crate::pizarra::{Marca as Resalte, Narracion, Pizarra};
use crate::{Algoritmo, Bloque, Campo, Entrada, Marca, Paso, Salida, Tipo};
use nucleo::ecoli;
use nucleo::formatos::{Lectura, Registro, abrir_fastq, leer_fasta_archivo};
use nucleo::mapeo::{
    self, Alineamiento, Ancla, CALIDAD_MAXIMA, Cadena, FRECUENCIA_MAXIMA, Indice, K,
    LONGITUD_MAXIMA_DE_LECTURA, Mapeo, Minimizador, ResumenMapeo, Tramo, W,
};
use nucleo::secuencia::{complemento_inverso, indice_base};
use std::collections::{HashMap, VecDeque};
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;
use std::sync::atomic::{AtomicU64, Ordering};

static CONTADOR: AtomicU64 = AtomicU64::new(0);

const BYTES_MAXIMOS_DE_REFERENCIA: u64 = 100_000_000;
const LONGITUD_MAXIMA_DE_REFERENCIA: usize = 20_000_000;
const FRACCION_DE_ADN_MINIMA: f64 = 0.9;
const TRAMO_LARGO: usize = 10_000;
const HUECO_MINIMO: usize = 100;
const PROFUNDIDAD_MINIMA: f64 = 10.0;
const VENTANAS_MINIMAS: usize = 10;
const MAXIMO_DE_FRASES: usize = 6;
const MAXIMO_DE_MARCAS: usize = 10;
const MAXIMO_DE_FILAS: usize = 2000;
const REGIONES_A_NOMBRAR: usize = 4;
const GENES_A_NOMBRAR: usize = 4;

const GRUPOS_DE_MAPQ: [(usize, usize, &str); 8] = [
    (0, 0, "0"),
    (1, 9, "1–9"),
    (10, 19, "10–19"),
    (20, 29, "20–29"),
    (30, 39, "30–39"),
    (40, 49, "40–49"),
    (50, 59, "50–59"),
    (60, 60, "60"),
];

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![Algoritmo {
        id: "puerto.mapear_lecturas",
        capitulo: "puerto",
        nombre: "Mapear lecturas a un genoma",
        resumen: "Ubica cada lectura de un FASTQ en el genoma con minimizadores, cadenas de anclas y un alineamiento con recorte (como minimap2) y guarda un SAM.",
        campos: &[
            Campo {
                id: "lecturas",
                etiqueta: "Lecturas cortas en FASTQ (.fastq o .fastq.gz)",
                tipo: Tipo::Archivo,
                ejemplo: "@ejemplo:lecturas",
            },
            Campo {
                id: "referencia",
                etiqueta: "Referencia: @ecoli (genoma incluido) o la ruta de un FASTA con un solo cromosoma",
                tipo: Tipo::Base,
                ejemplo: "@ecoli",
            },
        ],
        funcion: "mapeo::mapear_lecturas",
        ejecutar: ejecutar_mapeo,
        pasos: Some(pasos_mapeo),
    }]
}

fn indice_de_ecoli() -> Result<&'static Indice, String> {
    static INDICE: OnceLock<Result<Indice, String>> = OnceLock::new();
    INDICE
        .get_or_init(|| mapeo::construir_indice(ecoli::genoma(), K, W, FRECUENCIA_MAXIMA))
        .as_ref()
        .map_err(|error| error.clone())
}

pub fn preparar() {
    let _ = indice_de_ecoli();
}

enum Referencia {
    Ecoli,
    Propia { nombre: String, secuencia: Vec<u8> },
}

impl Referencia {
    fn es_ecoli(&self) -> bool {
        matches!(self, Referencia::Ecoli)
    }
}

fn ejecutar_mapeo(entrada: &Entrada) -> Result<Salida, String> {
    let lecturas = entrada.ruta("lecturas")?;
    let referencia = leer_referencia(entrada)?;
    let destino = crate::ruta_de_salida(&lecturas, "mapeo", "sam");
    let resumen = mapear_y_guardar(&lecturas, &referencia, &destino)?;
    let cobertura = Cobertura::medir(&resumen);
    let mut salida = Salida::default();
    agregar_valores(&mut salida, &resumen, &cobertura);
    salida.agregar(Bloque::archivo("Alineamientos (SAM)", &destino));
    agregar_graficos(&mut salida, &resumen, &cobertura);
    agregar_regiones(&mut salida, &cobertura, &referencia);
    explicar(&mut salida, &resumen, &cobertura, &referencia);
    Ok(salida)
}

fn leer_referencia(entrada: &Entrada) -> Result<Referencia, String> {
    let texto = entrada.texto("referencia")?;
    let pedido = texto.trim();
    if pedido == "@ecoli" {
        return Ok(Referencia::Ecoli);
    }
    if pedido.starts_with('@') {
        return Err(format!(
            "No conozco la referencia «{pedido}». Escribe @ecoli (el genoma incluido) o la ruta de un archivo FASTA."
        ));
    }
    let ruta = Path::new(pedido);
    let datos = fs::metadata(ruta).map_err(|_| {
        format!(
            "No existe el archivo «{}». Escribe @ecoli o la ruta de un archivo FASTA.",
            ruta.display()
        )
    })?;
    if !datos.is_file() {
        return Err(format!("«{}» no es un archivo.", ruta.display()));
    }
    if datos.len() > BYTES_MAXIMOS_DE_REFERENCIA {
        return Err(format!(
            "El FASTA de referencia pesa {} bytes y el máximo es {}.",
            miles(datos.len() as usize),
            miles(BYTES_MAXIMOS_DE_REFERENCIA as usize)
        ));
    }
    let registros: Vec<Registro> = leer_fasta_archivo(ruta)?
        .into_iter()
        .filter(|registro| !registro.secuencia.is_empty())
        .collect();
    revisar_referencia(registros)
}

fn revisar_referencia(mut registros: Vec<Registro>) -> Result<Referencia, String> {
    if registros.is_empty() {
        return Err("El FASTA de referencia no tiene ninguna secuencia.".to_string());
    }
    if registros.len() > 1 {
        return Err(format!(
            "El FASTA tiene {} secuencias. Este mapeo usa un solo cromosoma: deja en el archivo solo esa secuencia.",
            miles(registros.len())
        ));
    }
    let registro = registros.remove(0);
    let longitud = registro.secuencia.len();
    if longitud > LONGITUD_MAXIMA_DE_REFERENCIA {
        return Err(format!(
            "La referencia mide {} pb y el máximo es {} pb, más que cualquier genoma bacteriano.",
            miles(longitud),
            miles(LONGITUD_MAXIMA_DE_REFERENCIA)
        ));
    }
    let de_adn = registro
        .secuencia
        .iter()
        .filter(|base| matches!(base, b'A' | b'C' | b'G' | b'T' | b'N'))
        .count();
    if (de_adn as f64) < FRACCION_DE_ADN_MINIMA * longitud as f64 {
        return Err(
            "La referencia no parece ADN: menos del 90 % de sus letras son A, C, G, T o N."
                .to_string(),
        );
    }
    let nombre = registro
        .nombre
        .split_whitespace()
        .next()
        .unwrap_or("referencia")
        .to_string();
    Ok(Referencia::Propia {
        nombre,
        secuencia: registro.secuencia,
    })
}

fn mapear_y_guardar(
    lecturas: &Path,
    referencia: &Referencia,
    destino: &Path,
) -> Result<ResumenMapeo, String> {
    let temporal = ruta_temporal(destino);
    let resultado = mapear(lecturas, referencia, &temporal).and_then(|resumen| {
        if resumen.lecturas == 0 {
            return Err("El archivo no tiene ninguna lectura FASTQ.".to_string());
        }
        fs::rename(&temporal, destino)
            .map_err(|error| format!("No se pudo guardar {}: {error}", destino.display()))?;
        Ok(resumen)
    });
    if resultado.is_err() {
        let _ = fs::remove_file(&temporal);
    }
    resultado
}

fn mapear(lecturas: &Path, referencia: &Referencia, sam: &Path) -> Result<ResumenMapeo, String> {
    match referencia {
        Referencia::Ecoli => mapeo::mapear_fastq(
            lecturas,
            indice_de_ecoli()?,
            ecoli::genoma(),
            ecoli::IDENTIFICADOR,
            sam,
        ),
        Referencia::Propia { nombre, secuencia } => {
            mapeo::mapear_lecturas(lecturas, secuencia, nombre, sam)
        }
    }
}

fn ruta_temporal(destino: &Path) -> PathBuf {
    let numero = CONTADOR.fetch_add(1, Ordering::Relaxed);
    let nombre = destino
        .file_name()
        .map(|nombre| nombre.to_string_lossy().to_string())
        .unwrap_or_else(|| "mapeo.sam".to_string());
    destino.with_file_name(format!(".{}-{numero}-{nombre}", std::process::id()))
}

fn fraccion(parte: usize, total: usize) -> f64 {
    if total == 0 {
        0.0
    } else {
        parte as f64 / total as f64
    }
}

fn contar(cantidad: usize, singular: &str, plural: &str) -> String {
    if cantidad == 1 {
        format!("1 {singular}")
    } else {
        format!("{} {plural}", miles(cantidad))
    }
}

fn desde_uno(inicio: usize, fin: usize) -> String {
    format!("{}–{}", miles(inicio + 1), miles(fin))
}

struct Region {
    inicio: usize,
    fin: usize,
    bases_con_lecturas: usize,
    suma_de_profundidad: f64,
}

impl Region {
    fn largo(&self) -> usize {
        self.fin - self.inicio
    }

    fn es_larga(&self) -> bool {
        self.largo() >= TRAMO_LARGO
    }

    fn profundidad_media(&self) -> f64 {
        if self.bases_con_lecturas == 0 {
            0.0
        } else {
            self.suma_de_profundidad / self.bases_con_lecturas as f64
        }
    }
}

fn unir_tramos(tramos: &[Tramo]) -> Vec<Region> {
    let mut regiones: Vec<Region> = Vec::new();
    for tramo in tramos {
        let largo = tramo.fin - tramo.inicio;
        let suma = tramo.profundidad_media * largo as f64;
        match regiones.last_mut() {
            Some(ultima) if tramo.inicio - ultima.fin < HUECO_MINIMO => {
                ultima.fin = tramo.fin;
                ultima.bases_con_lecturas += largo;
                ultima.suma_de_profundidad += suma;
            }
            _ => regiones.push(Region {
                inicio: tramo.inicio,
                fin: tramo.fin,
                bases_con_lecturas: largo,
                suma_de_profundidad: suma,
            }),
        }
    }
    regiones
}

struct Cobertura {
    fraccion: f64,
    regiones: Vec<Region>,
    dirigida: bool,
    huecos: Vec<(usize, usize)>,
}

impl Cobertura {
    fn medir(resumen: &ResumenMapeo) -> Cobertura {
        let cubierta = fraccion(resumen.bases_cubiertas, resumen.longitud_referencia);
        let regiones = unir_tramos(&resumen.tramos);
        let en_largas: usize = regiones
            .iter()
            .filter(|region| region.es_larga())
            .map(|region| region.bases_con_lecturas)
            .sum();
        let dirigida = cubierta < 0.5
            && en_largas > 0
            && en_largas as f64 >= 0.9 * resumen.bases_cubiertas as f64;
        let huecos = if cubierta >= 0.5 || dirigida {
            posibles_deleciones(&regiones, dirigida, resumen.profundidad_media)
        } else {
            Vec::new()
        };
        Cobertura {
            fraccion: cubierta,
            regiones,
            dirigida,
            huecos,
        }
    }

    fn largas(&self) -> Vec<&Region> {
        self.regiones
            .iter()
            .filter(|region| region.es_larga())
            .collect()
    }

    fn islas(&self) -> Vec<&Region> {
        self.regiones
            .iter()
            .filter(|region| !region.es_larga())
            .collect()
    }
}

fn posibles_deleciones(
    regiones: &[Region],
    dirigida: bool,
    profundidad_media: f64,
) -> Vec<(usize, usize)> {
    let mut huecos = Vec::new();
    if profundidad_media < PROFUNDIDAD_MINIMA {
        return huecos;
    }
    for par in regiones.windows(2) {
        let antes = &par[0];
        let despues = &par[1];
        let sin_lecturas = despues.inicio - antes.fin;
        if dirigida && (sin_lecturas >= TRAMO_LARGO || !antes.es_larga() || !despues.es_larga()) {
            continue;
        }
        huecos.push((antes.fin, despues.inicio));
    }
    huecos
}

fn agregar_valores(salida: &mut Salida, resumen: &ResumenMapeo, cobertura: &Cobertura) {
    let mapeadas = resumen.mapeadas;
    salida.agregar(Bloque::valor(
        "Lecturas mapeadas",
        format!(
            "{} de {} ({})",
            miles(mapeadas),
            miles(resumen.lecturas),
            porcentaje(fraccion(mapeadas, resumen.lecturas), 1)
        ),
    ));
    salida.agregar(Bloque::valor(
        "Profundidad media donde hay lecturas",
        format!("{}×", decimal(resumen.profundidad_media, 1)),
    ));
    salida.agregar(Bloque::valor(
        "Referencia con al menos una lectura",
        format!(
            "{} ({} de {} pb)",
            porcentaje(cobertura.fraccion, 1),
            miles(resumen.bases_cubiertas),
            miles(resumen.longitud_referencia)
        ),
    ));
    salida.agregar(Bloque::valor(
        "Lecturas con MAPQ 0 (posición ambigua)",
        format!(
            "{} ({} de las mapeadas)",
            miles(resumen.calidades[0]),
            porcentaje(fraccion(resumen.calidades[0], mapeadas), 1)
        ),
    ));
    salida.agregar(Bloque::valor(
        "Lecturas mapeadas con un extremo recortado (S)",
        format!(
            "{} ({})",
            miles(resumen.con_recorte),
            porcentaje(fraccion(resumen.con_recorte, mapeadas), 1)
        ),
    ));
    salida.agregar(Bloque::valor(
        "Diferencias con la referencia",
        format!(
            "{} por cada 100 bases alineadas ({} cambios, {} bases insertadas, {} borradas)",
            decimal(100.0 * tasa_de_diferencias(resumen), 2),
            miles(resumen.diferencias),
            miles(resumen.insertadas),
            miles(resumen.borradas)
        ),
    ));
    salida.agregar(Bloque::valor(
        "Índice de minimizadores",
        format!(
            "{} minimizadores (k = {K}, w = {W}); {} descartados por aparecer más de {FRECUENCIA_MAXIMA} veces",
            miles(resumen.minimizadores_en_indice),
            miles(resumen.minimizadores_descartados)
        ),
    ));
}

fn tasa_de_diferencias(resumen: &ResumenMapeo) -> f64 {
    let alineadas = resumen.coincidencias + resumen.diferencias + resumen.insertadas;
    let cambios = resumen.diferencias + resumen.insertadas + resumen.borradas;
    fraccion(cambios, alineadas)
}

fn ventana_en_texto(largo_ventana: usize) -> String {
    if largo_ventana >= 1000 && largo_ventana.is_multiple_of(1000) {
        format!("{} kb", miles(largo_ventana / 1000))
    } else {
        format!("{} pb", miles(largo_ventana))
    }
}

fn agregar_graficos(salida: &mut Salida, resumen: &ResumenMapeo, cobertura: &Cobertura) {
    let (escala, unidad) = if resumen.longitud_referencia >= 1_000_000 {
        (1e6, "Mb")
    } else {
        (1e3, "kb")
    };
    let mut puntos = Vec::with_capacity(resumen.profundidad_por_ventana.len());
    for (numero, &profundidad) in resumen.profundidad_por_ventana.iter().enumerate() {
        let centro = (numero as f64 + 0.5) * resumen.largo_ventana as f64;
        puntos.push([centro / escala, profundidad]);
    }
    let marcas = cobertura
        .huecos
        .iter()
        .take(MAXIMO_DE_MARCAS)
        .map(|&(inicio, fin)| Marca {
            x: (inicio + fin) as f64 / 2.0 / escala,
            texto: format!("Sin lecturas {}: ¿deleción?", desde_uno(inicio, fin)),
        })
        .collect();
    salida.agregar(Bloque::linea(
        &format!(
            "Profundidad a lo largo de la referencia (promedio en ventanas de {})",
            ventana_en_texto(resumen.largo_ventana)
        ),
        &format!("Posición en la referencia ({unidad})"),
        "Lecturas por base",
        &puntos,
        marcas,
    ));
    let mut categorias = Vec::new();
    let mut valores = Vec::new();
    for (desde, hasta, nombre) in GRUPOS_DE_MAPQ {
        let cantidad: usize = (desde..=hasta)
            .filter_map(|calidad| resumen.calidades.get(calidad))
            .sum();
        categorias.push(nombre.to_string());
        valores.push(cantidad as f64);
    }
    salida.agregar(Bloque::barras(
        "Calidad de mapeo (MAPQ) de las lecturas mapeadas",
        "MAPQ (0 = posición ambigua, 60 = única)",
        "Lecturas",
        categorias,
        valores,
    ));
}

fn genes_de_la_region(inicio: usize, fin: usize) -> Vec<&'static str> {
    ecoli::genes()
        .iter()
        .filter(|elegido| elegido.inicio <= fin && elegido.fin > inicio)
        .map(|elegido| {
            if elegido.nombre.is_empty() {
                elegido.locus.as_str()
            } else {
                elegido.nombre.as_str()
            }
        })
        .collect()
}

fn describir_genes(nombres: &[&str]) -> String {
    match nombres {
        [] => "—".to_string(),
        [primero, .., ultimo] if nombres.len() > GENES_A_NOMBRAR => {
            format!("{} genes, de {primero} a {ultimo}", miles(nombres.len()))
        }
        _ => nombres.join(", "),
    }
}

fn agregar_regiones(salida: &mut Salida, cobertura: &Cobertura, referencia: &Referencia) {
    let mut columnas = vec![
        "Inicio (contando desde 1)".to_string(),
        "Fin".to_string(),
        "Longitud (pb)".to_string(),
        "Profundidad media".to_string(),
    ];
    if referencia.es_ecoli() {
        columnas.push("Genes de proteínas de E. coli".to_string());
    }
    let mut filas = Vec::new();
    for region in cobertura.regiones.iter().take(MAXIMO_DE_FILAS) {
        let mut fila = vec![
            miles(region.inicio + 1),
            miles(region.fin),
            miles(region.largo()),
            format!("{}×", decimal(region.profundidad_media(), 1)),
        ];
        if referencia.es_ecoli() {
            fila.push(describir_genes(&genes_de_la_region(
                region.inicio,
                region.fin,
            )));
        }
        filas.push(fila);
    }
    salida.agregar(Bloque::Tabla {
        etiqueta: "Regiones de la referencia con lecturas (se unen los tramos separados por menos de 100 pb sin lecturas)".to_string(),
        columnas,
        filas,
        total: cobertura.regiones.len(),
    });
}

fn explicar(
    salida: &mut Salida,
    resumen: &ResumenMapeo,
    cobertura: &Cobertura,
    referencia: &Referencia,
) {
    let mut frases = vec![frase_de_mapeadas(resumen)];
    if resumen.demasiado_largas > 0 {
        frases.push(format!(
            "Hay {} de más de {} bases sin mapear: este mapeador es para lecturas cortas, como las de Illumina.",
            contar(resumen.demasiado_largas, "lectura", "lecturas"),
            miles(LONGITUD_MAXIMA_DE_LECTURA)
        ));
    }
    if resumen.mapeadas == 0 {
        frases.push(
            "Ninguna lectura se pudo alinear con la referencia, así que el SAM solo tiene lecturas sin mapear (bandera 4) y no hay profundidad que medir."
                .to_string(),
        );
    } else {
        frases.push(frase_de_cobertura(resumen, cobertura));
        if cobertura.dirigida && !cobertura.huecos.is_empty() {
            frases.push(frase_de_huecos_en_regiones(cobertura));
        }
        if cobertura.dirigida && !cobertura.islas().is_empty() && resumen.calidades[0] > 0 {
            frases.push(frase_de_islas(cobertura, referencia));
        }
        frases.push(frase_de_ambiguas(resumen, referencia));
        frases.push(frase_de_profundidad(resumen));
        if fraccion(resumen.con_recorte, resumen.mapeadas) > 0.05 {
            frases.push(format!(
                "El {} de las lecturas mapeadas tiene un extremo recortado (S en el CIGAR) porque no se parece al genoma: son restos de adaptador o colas de mala calidad. Pasa del 5 %, así que conviene limpiarlas antes con «Recortar lecturas».",
                porcentaje(fraccion(resumen.con_recorte, resumen.mapeadas), 1)
            ));
        }
        if tasa_de_diferencias(resumen) > 0.01 {
            frases.push(format!(
                "Hay {} diferencias por cada 100 bases alineadas, más que el 1 % que dejan los errores de secuenciación con calidad Q20: la muestra puede ser otra cepa o tener baja calidad.",
                decimal(100.0 * tasa_de_diferencias(resumen), 1)
            ));
        }
        if cobertura.fraccion >= 0.5
            && let Some(frase) = frase_de_uniformidad(resumen, cobertura)
        {
            frases.push(frase);
        }
    }
    frases.truncate(MAXIMO_DE_FRASES);
    for frase in frases {
        salida.explicar(frase);
    }
}

fn frase_de_mapeadas(resumen: &ResumenMapeo) -> String {
    let mapeadas = fraccion(resumen.mapeadas, resumen.lecturas);
    let cuantas = format!(
        "el {} de las lecturas ({} de {})",
        porcentaje(mapeadas, 1),
        miles(resumen.mapeadas),
        miles(resumen.lecturas)
    );
    if mapeadas >= 0.95 {
        format!("Se ubicó {cuantas}: la muestra es de esta referencia o de una cepa muy cercana.")
    } else if mapeadas >= 0.70 {
        format!(
            "Se ubicó {cuantas}, menos del 95 % que se espera si la muestra fuera este mismo genoma: puede ser otra cepa, con genes que la referencia no tiene, o una contaminación parcial. Busca con BLASTN qué son las lecturas sin mapear (bandera 4 en el SAM)."
        )
    } else {
        format!(
            "Solo se ubicó {cuantas}: probablemente vienen de otro organismo o hay mucha contaminación. Revisa que la referencia sea la correcta."
        )
    }
}

fn lista_de_regiones(regiones: &[&Region]) -> String {
    let mut nombradas: Vec<String> = regiones
        .iter()
        .take(REGIONES_A_NOMBRAR)
        .map(|region| desde_uno(region.inicio, region.fin))
        .collect();
    if regiones.len() > REGIONES_A_NOMBRAR {
        nombradas.push(format!(
            "{} más",
            miles(regiones.len() - REGIONES_A_NOMBRAR)
        ));
    }
    match nombradas.split_last() {
        Some((ultima, resto)) if !resto.is_empty() => format!("{} y {ultima}", resto.join(", ")),
        _ => nombradas.join(""),
    }
}

fn frase_de_cobertura(resumen: &ResumenMapeo, cobertura: &Cobertura) -> String {
    let cubierta = porcentaje(cobertura.fraccion, 1);
    if cobertura.dirigida {
        let largas = cobertura.largas();
        let regiones = contar(largas.len(), "región", "regiones");
        return format!(
            "Solo el {cubierta} de la referencia tiene lecturas, y casi todas caen en {regiones}: {} (contando desde 1). Las lecturas vienen de regiones elegidas (secuenciación dirigida o enriquecida), no del genoma completo, así que las zonas sin lecturas no son deleciones.",
            lista_de_regiones(&largas)
        );
    }
    if cobertura.fraccion < 0.5 {
        return format!(
            "Solo el {cubierta} de la referencia tiene lecturas, repartidas en {}: hay pocas lecturas para una referencia de este tamaño, y con más lecturas se cubriría el resto.",
            contar(cobertura.regiones.len(), "región", "regiones")
        );
    }
    if resumen.profundidad_media < PROFUNDIDAD_MINIMA {
        return format!(
            "El {cubierta} de la referencia tiene lecturas, pero con {}× de profundidad media es normal que queden huecos por azar: para buscar deleciones hacen falta al menos 10×.",
            decimal(resumen.profundidad_media, 1)
        );
    }
    if cobertura.huecos.is_empty() {
        return format!(
            "El {cubierta} de la referencia tiene lecturas y no hay tramos de 100 pb o más sin ninguna lectura entre zonas cubiertas: no se ven deleciones grandes."
        );
    }
    format!(
        "Hay {} entre zonas cubiertas. Con {}× de profundidad media eso no pasa por azar: {} de la muestra respecto de la referencia.",
        describir_huecos(&cobertura.huecos),
        decimal(resumen.profundidad_media, 1),
        conclusion_de_huecos(&cobertura.huecos)
    )
}

fn describir_huecos(huecos: &[(usize, usize)]) -> String {
    let (inicio, fin) = mayor_hueco(huecos).unwrap_or_default();
    let tramo = format!(
        "{}, {} pb, contando desde 1",
        desde_uno(inicio, fin),
        miles(fin - inicio)
    );
    if huecos.len() == 1 {
        format!("1 tramo de 100 pb o más sin ninguna lectura ({tramo})")
    } else {
        format!(
            "{} tramos de 100 pb o más sin ninguna lectura (el mayor: {tramo})",
            miles(huecos.len())
        )
    }
}

fn conclusion_de_huecos(huecos: &[(usize, usize)]) -> &'static str {
    if huecos.len() == 1 {
        "es una posible deleción"
    } else {
        "cada uno es una posible deleción"
    }
}

fn mayor_hueco(huecos: &[(usize, usize)]) -> Option<(usize, usize)> {
    let mut mayor: Option<(usize, usize)> = None;
    for &(inicio, fin) in huecos {
        if mayor.is_none_or(|(otro_inicio, otro_fin)| fin - inicio > otro_fin - otro_inicio) {
            mayor = Some((inicio, fin));
        }
    }
    mayor
}

fn frase_de_huecos_en_regiones(cobertura: &Cobertura) -> String {
    format!(
        "Dentro de esas regiones hay {}: {}.",
        describir_huecos(&cobertura.huecos),
        conclusion_de_huecos(&cobertura.huecos)
    )
}

fn frase_de_islas(cobertura: &Cobertura, referencia: &Referencia) -> String {
    let islas = cobertura.islas();
    let bases: usize = islas.iter().map(|isla| isla.bases_con_lecturas).sum();
    let mut ejemplos: Vec<&str> = Vec::new();
    if referencia.es_ecoli() {
        for isla in &islas {
            if let Some(&primero) = genes_de_la_region(isla.inicio, isla.fin).first()
                && !ejemplos.contains(&primero)
            {
                ejemplos.push(primero);
            }
            if ejemplos.len() == GENES_A_NOMBRAR {
                break;
            }
        }
    }
    let donde = if ejemplos.is_empty() {
        String::new()
    } else {
        format!(", en genes como {}", ejemplos.join(", "))
    };
    if islas.len() == 1 {
        return format!(
            "La única zona corta con lecturas fuera de esas regiones ({} pb{donde}) probablemente es la copia de una repetición donde cayeron, al azar, lecturas ambiguas (MAPQ 0).",
            miles(bases)
        );
    }
    format!(
        "Las {} zonas cortas con lecturas fuera de esas regiones ({} pb en total{donde}) probablemente son copias de repeticiones donde cayeron, al azar, lecturas ambiguas (MAPQ 0).",
        miles(islas.len()),
        miles(bases)
    )
}

fn frase_de_ambiguas(resumen: &ResumenMapeo, referencia: &Referencia) -> String {
    let ambiguas = fraccion(resumen.calidades[0], resumen.mapeadas);
    if ambiguas >= 0.01 {
        let ejemplo = if referencia.es_ecoli() {
            " (en E. coli, sobre todo los 7 operones de ARN ribosomal, los elementos IS y los genes gemelos tufA y tufB)"
        } else {
            ""
        };
        format!(
            "El {} de las lecturas mapeadas tiene MAPQ 0: su secuencia está igual en otra parte de la referencia{ejemplo}, así que su posición se eligió al azar entre las copias. No las uses para buscar variantes.",
            porcentaje(ambiguas, 1)
        )
    } else {
        format!(
            "Solo el {} de las lecturas mapeadas tiene MAPQ 0 (posición ambigua): casi todas tienen un único lugar posible en la referencia.",
            porcentaje(ambiguas, 1)
        )
    }
}

fn frase_de_profundidad(resumen: &ResumenMapeo) -> String {
    let profundidad = decimal(resumen.profundidad_media, 1);
    if resumen.profundidad_media < PROFUNDIDAD_MINIMA {
        format!(
            "Donde hay lecturas, cada base está cubierta en promedio por {profundidad} lecturas: para buscar variantes con confianza se suelen pedir al menos 10."
        )
    } else {
        format!(
            "Donde hay lecturas, cada base está cubierta en promedio por {profundidad} lecturas: alcanza para buscar variantes (se suelen pedir al menos 10). El siguiente paso es «Llamado de variantes» con este SAM."
        )
    }
}

fn ventanas_completas(resumen: &ResumenMapeo, cobertura: &Cobertura) -> Vec<f64> {
    let largo_ventana = resumen.largo_ventana.max(1);
    let mut profundidades = Vec::new();
    for region in &cobertura.regiones {
        let primera = region.inicio.div_ceil(largo_ventana);
        let ultima = region.fin / largo_ventana;
        for numero in primera..ultima {
            if let Some(&profundidad) = resumen.profundidad_por_ventana.get(numero) {
                profundidades.push(profundidad);
            }
        }
    }
    profundidades
}

fn frase_de_uniformidad(resumen: &ResumenMapeo, cobertura: &Cobertura) -> Option<String> {
    let mut profundidades = ventanas_completas(resumen, cobertura);
    if profundidades.len() < VENTANAS_MINIMAS {
        return None;
    }
    profundidades.sort_by(f64::total_cmp);
    let mediana = profundidades[profundidades.len() / 2];
    let altas = profundidades
        .iter()
        .filter(|&&profundidad| profundidad > 2.0 * mediana)
        .count();
    let bajas = profundidades
        .iter()
        .filter(|&&profundidad| profundidad < 0.5 * mediana)
        .count();
    let ventana = ventana_en_texto(resumen.largo_ventana);
    let mediana = decimal(mediana, 1);
    if fraccion(altas + bajas, profundidades.len()) > 0.05 {
        Some(format!(
            "La profundidad es despareja entre las ventanas de {ventana}: {} más del doble de la mediana ({mediana}×), lo que sugiere regiones duplicadas, y {} menos de la mitad, típico de zonas difíciles de secuenciar, como las muy ricas en GC.",
            contar(altas, "tiene", "tienen"),
            contar(bajas, "tiene", "tienen")
        ))
    } else {
        Some(format!(
            "La profundidad es pareja: el {} de las ventanas de {ventana} está entre la mitad y el doble de la mediana ({mediana}×).",
            porcentaje(1.0 - fraccion(altas + bajas, profundidades.len()), 1)
        ))
    }
}

const MAXIMO_DE_LECTURAS_ANIMABLES: usize = 200_000;
const MAXIMO_DE_BASES_ANIMABLES: usize = 50_000_000;
const BASES_DE_LA_DEMOSTRACION: usize = 60;
const ANCLAS_A_MOSTRAR: usize = 8;
const ENLACES_A_MOSTRAR: usize = 3;
const CADENAS_A_MOSTRAR: usize = 5;
const CADENAS_DETALLADAS: usize = 3;
const LECTURAS_UNA_A_UNA: usize = 4;
const CORTES_DEL_RECUENTO: usize = 9;
const LETRAS_A_MOSTRAR: usize = 24;
const TRAMOS_A_MOSTRAR: usize = 5;
const POSICIONES_A_MOSTRAR: usize = 6;

struct Contexto<'a> {
    genoma: &'a [u8],
    nombre: &'a str,
    indice: &'a Indice,
}

#[derive(Clone)]
struct Candidato {
    reverso: bool,
    alineamiento: Alineamiento,
    anclas: usize,
}

#[derive(Clone, Copy, PartialEq)]
enum Efecto {
    Nuevo,
    Reemplaza,
    Repetido,
}

struct Intento {
    cadena: Cadena,
    diagonales: (i64, i64),
    alineamiento: Option<Alineamiento>,
    aceptado: bool,
    efecto: Option<Efecto>,
}

#[derive(Default)]
struct Analisis {
    demasiado_larga: bool,
    minimizadores: Vec<Minimizador>,
    anclas: Vec<Ancla>,
    todas: Vec<Cadena>,
    candidatas: Vec<Cadena>,
    intentos: Vec<Intento>,
    candidatos: Vec<Candidato>,
    mejor: i32,
    elegido: usize,
    empatados: usize,
    segundo: i32,
    mapeo: Option<Mapeo>,
}

struct Instante {
    inicio: usize,
    kmero: Vec<u8>,
    valido: bool,
    valor: u64,
    quitados: usize,
    caducados: usize,
    ventana: Vec<Minimizador>,
    menor: Option<Minimizador>,
    nuevo: bool,
}

fn corto(texto: &[u8]) -> String {
    let visibles = String::from_utf8_lossy(&texto[..texto.len().min(LETRAS_A_MOSTRAR)]);
    if texto.len() > LETRAS_A_MOSTRAR {
        format!("{visibles}…")
    } else {
        visibles.to_string()
    }
}

fn huella(valor: u64) -> String {
    format!("{:06x}", valor >> 40)
}

fn hebra(reverso: bool) -> &'static str {
    if reverso { "reversa" } else { "directa" }
}

fn cigar_en_texto(cigar: &[(u32, u8)]) -> String {
    cigar
        .iter()
        .map(|&(cantidad, operacion)| format!("{cantidad}{}", operacion as char))
        .collect()
}

fn nombre_de(lectura: &Lectura) -> String {
    String::from_utf8_lossy(&lectura.nombre)
        .split_whitespace()
        .next()
        .unwrap_or("*")
        .to_string()
}

fn diagonal_de(ancla: &Ancla) -> i64 {
    ancla.en_referencia as i64 - ancla.en_lectura as i64
}

fn seguir_minimizadores(secuencia: &[u8], k: usize, w: usize) -> (Vec<Instante>, Vec<Minimizador>) {
    let mut instantes = Vec::new();
    let mut elegidos: Vec<Minimizador> = Vec::new();
    if k == 0 || k > 31 || w == 0 || secuencia.len() < k {
        return (instantes, elegidos);
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
        let mut valido = false;
        let mut valor = 0;
        let mut quitados = 0;
        if bases_seguidas >= k && directo != inverso {
            let candidato = Minimizador {
                valor: mapeo::barajar(directo.min(inverso)),
                posicion: inicio,
                reverso: inverso < directo,
            };
            valido = true;
            valor = candidato.valor;
            while ventana
                .back()
                .is_some_and(|ultimo| ultimo.valor > candidato.valor)
            {
                ventana.pop_back();
                quitados += 1;
            }
            ventana.push_back(candidato);
        }
        let mut caducados = 0;
        while ventana
            .front()
            .is_some_and(|primero| primero.posicion + w <= inicio)
        {
            ventana.pop_front();
            caducados += 1;
        }
        let mut menor = None;
        let mut nuevo = false;
        if inicio + 1 >= w
            && let Some(&primero) = ventana.front()
        {
            menor = Some(primero);
            if elegidos.last() != Some(&primero) {
                elegidos.push(primero);
                nuevo = true;
            }
        }
        instantes.push(Instante {
            inicio,
            kmero: secuencia[inicio..inicio + k].to_vec(),
            valido,
            valor,
            quitados,
            caducados,
            ventana: ventana.iter().copied().collect(),
            menor,
            nuevo,
        });
    }
    (instantes, elegidos)
}

fn agregar_sin_repetir(candidatos: &mut Vec<Candidato>, nuevo: Candidato) -> Efecto {
    for existente in candidatos.iter_mut() {
        if existente.reverso == nuevo.reverso
            && existente.alineamiento.se_superpone_con(&nuevo.alineamiento)
        {
            if nuevo.alineamiento.puntaje > existente.alineamiento.puntaje {
                *existente = nuevo;
                return Efecto::Reemplaza;
            }
            return Efecto::Repetido;
        }
    }
    candidatos.push(nuevo);
    Efecto::Nuevo
}

fn analizar(lectura: &Lectura, indice: &Indice, genoma: &[u8]) -> Analisis {
    let mut analisis = Analisis::default();
    let secuencia = &lectura.secuencia;
    if secuencia.len() > LONGITUD_MAXIMA_DE_LECTURA {
        analisis.demasiado_larga = true;
        return analisis;
    }
    analisis.minimizadores = mapeo::minimizadores(secuencia, indice.k, indice.w);
    analisis.anclas = mapeo::buscar_anclas(secuencia, indice);
    analisis.todas = mapeo::encadenar(&analisis.anclas, indice.k);
    analisis.candidatas = mapeo::elegir_candidatas(analisis.todas.clone());
    if analisis.candidatas.is_empty() {
        return analisis;
    }
    let inversa = complemento_inverso(secuencia);
    let mut candidatos: Vec<Candidato> = Vec::new();
    for cadena in &analisis.candidatas {
        let orientada = if cadena.reverso { &inversa } else { secuencia };
        let diagonales = cadena.diagonales();
        let alineamiento =
            mapeo::alinear_con_recorte(orientada, genoma, diagonales.0, diagonales.1);
        let mut intento = Intento {
            cadena: cadena.clone(),
            diagonales,
            alineamiento: alineamiento.clone(),
            aceptado: false,
            efecto: None,
        };
        if let Some(alineamiento) = alineamiento
            && alineamiento.puntaje >= mapeo::PUNTAJE_MINIMO_DE_ALINEAMIENTO
        {
            intento.aceptado = true;
            let nuevo = Candidato {
                reverso: cadena.reverso,
                alineamiento,
                anclas: cadena.anclas.len(),
            };
            intento.efecto = Some(agregar_sin_repetir(&mut candidatos, nuevo));
        }
        analisis.intentos.push(intento);
    }
    analisis.candidatos = candidatos.clone();
    let Some(mejor) = candidatos
        .iter()
        .map(|candidato| candidato.alineamiento.puntaje)
        .max()
    else {
        return analisis;
    };
    let empatados: Vec<usize> = (0..candidatos.len())
        .filter(|&numero| candidatos[numero].alineamiento.puntaje == mejor)
        .collect();
    let elegido = empatados[(mapeo::desempate(&lectura.nombre) % empatados.len() as u64) as usize];
    let mut segundo = 0;
    for (numero, candidato) in candidatos.iter().enumerate() {
        if numero != elegido {
            segundo = segundo.max(candidato.alineamiento.puntaje);
        }
    }
    let ganador = candidatos.swap_remove(elegido);
    analisis.mejor = mejor;
    analisis.elegido = elegido;
    analisis.empatados = empatados.len();
    analisis.segundo = segundo;
    analisis.mapeo = Some(Mapeo {
        reverso: ganador.reverso,
        calidad_mapeo: mapeo::calidad_de_mapeo(mejor, segundo, ganador.anclas),
        alineamiento: ganador.alineamiento,
    });
    analisis
}

fn leer_lecturas_a_narrar(ruta: &Path) -> Result<Vec<Lectura>, String> {
    let mut lecturas = Vec::new();
    let mut bases = 0;
    for lectura in abrir_fastq(ruta)? {
        let lectura = lectura?;
        bases += lectura.secuencia.len();
        if lecturas.len() >= MAXIMO_DE_LECTURAS_ANIMABLES || bases > MAXIMO_DE_BASES_ANIMABLES {
            return Err(format!(
                "El paso a paso admite hasta {} lecturas o {} bases y este archivo tiene más. Prueba con un archivo más pequeño.",
                miles(MAXIMO_DE_LECTURAS_ANIMABLES),
                miles(MAXIMO_DE_BASES_ANIMABLES)
            ));
        }
        lecturas.push(lectura);
    }
    if lecturas.is_empty() {
        return Err("El archivo no tiene ninguna lectura FASTQ.".to_string());
    }
    Ok(lecturas)
}

fn pasos_mapeo(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let ruta = entrada.ruta("lecturas")?;
    let referencia = leer_referencia(entrada)?;
    let (genoma, nombre): (&[u8], &str) = match &referencia {
        Referencia::Ecoli => (ecoli::genoma(), ecoli::IDENTIFICADOR),
        Referencia::Propia { nombre, secuencia } => (secuencia.as_slice(), nombre.as_str()),
    };
    let propio;
    let indice: &Indice = match &referencia {
        Referencia::Ecoli => indice_de_ecoli()?,
        Referencia::Propia { .. } => {
            propio = mapeo::construir_indice(genoma, K, W, FRECUENCIA_MAXIMA)?;
            &propio
        }
    };
    let lecturas = leer_lecturas_a_narrar(&ruta)?;
    let mapeos = mapeo::mapear_en_paralelo(&lecturas, indice, genoma);
    let contexto = Contexto {
        genoma,
        nombre,
        indice,
    };
    let mut narracion = Narracion::nueva(400);
    narrar_entrada(&contexto, &lecturas, &mut narracion);
    narrar_demostracion(&contexto, &mut narracion);
    narrar_indice(&contexto, &mut narracion);
    narrar_lectura_del_archivo(&contexto, &lecturas, &mapeos, &mut narracion);
    let resumen = narrar_recuento(&contexto, &lecturas, &mapeos, &mut narracion);
    Ok(narracion.terminar(paso_final(&contexto, &resumen)))
}

fn narrar_entrada(c: &Contexto, lecturas: &[Lectura], narracion: &mut Narracion) {
    narracion.agregar(
        Pizarra::nueva()
            .variable(
                "referencia",
                format!("{} ({} pb)", c.nombre, miles(c.genoma.len())),
            )
            .variable("k", K)
            .variable("w", W)
            .variable("frecuencia_maxima", FRECUENCIA_MAXIMA)
            .secuencia(c.genoma, &[(0, K, Resalte::Actual)])
            .paso(
                "let indice = construir_indice(referencia, K, W, FRECUENCIA_MAXIMA)?;",
                format!(
                    "Tienes {} lecturas cortas y un genoma de referencia ({}, {} pb). Mapear es decidir en qué lugar del genoma cae cada lectura. Para no comparar cada una con todo el genoma, primero se arma un índice de «palabras» de k = {K} letras, una por ventana de w = {W}.",
                    miles(lecturas.len()),
                    c.nombre,
                    miles(c.genoma.len())
                ),
            ),
    );
}

fn paso_de_instante(tramo: &[u8], instantes: &[Instante], numero: usize) -> Paso {
    let instante = &instantes[numero];
    let elegidos = instantes[..=numero]
        .iter()
        .filter(|otro| otro.nuevo)
        .count();
    let inicio_de_ventana = (instante.inicio + 1).saturating_sub(W);
    let mut marcas = vec![
        (inicio_de_ventana, instante.inicio + K, Resalte::Comparado),
        (instante.inicio, instante.inicio + K, Resalte::Actual),
    ];
    if let Some(menor) = instante.menor {
        marcas.push((menor.posicion, menor.posicion + K, Resalte::Hallado));
    }
    let valor = if instante.valido {
        huella(instante.valor)
    } else {
        "—".to_string()
    };
    let (linea, nota) = if !instante.valido {
        (
            "if bases_seguidas >= k && directo != inverso {",
            format!(
                "El k-mero {} tiene una letra que no es A, C, G ni T, o es igual a su complemento inverso: no se usa.",
                corto(&instante.kmero)
            ),
        )
    } else if numero < 2 {
        (
            "valor: barajar(directo.min(inverso)),",
            format!(
                "El k-mero {} ocupa las posiciones {}–{}. Su valor es un número revuelto, la «huella» {valor}: se saca del menor entre el k-mero y su complemento inverso, así sale igual lo leas en una hebra o en la otra.",
                corto(&instante.kmero),
                instante.inicio + 1,
                instante.inicio + K
            ),
        )
    } else if numero == 2 {
        let quitados = if instante.quitados == 0 {
            "Nadie en la ventana tenía una huella mayor, así que solo se añade al final."
                .to_string()
        } else {
            format!(
                "Se quitan {} candidato(s) con huella mayor: mientras este k-mero siga en la ventana, ya no pueden ser el menor.",
                instante.quitados
            )
        };
        (
            "ventana.push_back(candidato);",
            format!(
                "El k-mero {} (huella {valor}) entra a la ventana de candidatos. {quitados}",
                corto(&instante.kmero)
            ),
        )
    } else if numero == W - 1 {
        let menor = instante.menor.map(|menor| menor.posicion + 1).unwrap_or(0);
        (
            "elegidos.push(menor);",
            format!(
                "Ya hay {W} k-meros: la primera ventana está completa. El primero de la lista tiene la huella menor y empieza en la posición {menor}: ese k-mero es el primer minimizador."
            ),
        )
    } else if instante.nuevo {
        let menor = instante.menor.map(|menor| menor.posicion + 1).unwrap_or(0);
        let salida = if instante.caducados > 0 {
            " El k-mero más viejo salió de la ventana (ya no cabe en w)."
        } else {
            ""
        };
        (
            "elegidos.push(menor);",
            format!(
                "La ventana avanzó una base.{salida} Ahora el menor empieza en la posición {menor}: es un minimizador nuevo."
            ),
        )
    } else {
        let menor = instante.menor.map(|menor| menor.posicion + 1).unwrap_or(0);
        (
            "&& elegidos.last() != Some(&menor)",
            format!(
                "La ventana avanzó una base y el menor sigue siendo el de la posición {menor}. Ya estaba elegido, así que no se repite."
            ),
        )
    };
    Pizarra::nueva()
        .variable("k", K)
        .variable("w", W)
        .variable("inicio", instante.inicio)
        .variable("valor", valor)
        .variable("elegidos", elegidos)
        .secuencia(tramo, &marcas)
        .lista(
            "Ventana de candidatos (el primero es el menor)",
            instante
                .ventana
                .iter()
                .map(|candidato| {
                    format!(
                        "posición {}: huella {}",
                        candidato.posicion + 1,
                        huella(candidato.valor)
                    )
                })
                .collect(),
            &[0],
        )
        .paso(linea, nota)
}

fn narrar_demostracion(c: &Contexto, narracion: &mut Narracion) {
    let tramo = &c.genoma[..c.genoma.len().min(BASES_DE_LA_DEMOSTRACION)];
    let (instantes, elegidos) = seguir_minimizadores(tramo, K, W);
    if instantes.len() < W {
        return;
    }
    let mut cuales = vec![0, 1, 2, W - 1, W, W + 1];
    cuales.retain(|&numero| numero < instantes.len());
    cuales.dedup();
    for numero in cuales {
        narracion.agregar(paso_de_instante(tramo, &instantes, numero));
    }
    let marcas: Vec<(usize, usize, Resalte)> = elegidos
        .iter()
        .map(|elegido| (elegido.posicion, elegido.posicion + K, Resalte::Hallado))
        .collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", K)
            .variable("w", W)
            .variable("inicio", instantes.len() - 1)
            .variable("valor", "—")
            .variable("elegidos", elegidos.len())
            .secuencia(tramo, &marcas)
            .paso(
                "elegidos.push(menor);",
                format!(
                    "En estas {} bases hay {} k-meros, pero solo {} minimizadores (marcados). Con el genoma completo se hace lo mismo, ventana tras ventana: así se guarda una fracción de las palabras sin perder ninguna región.",
                    tramo.len(),
                    miles(instantes.len()),
                    elegidos.len()
                ),
            ),
    );
}

fn narrar_indice(c: &Contexto, narracion: &mut Narracion) {
    let indice = c.indice;
    let largo = c.genoma.len();
    let guardados = indice.semillas.len();
    let total = guardados + indice.apariciones_descartadas;
    let referencia = format!("{} ({} pb)", c.nombre, miles(largo));
    let cada = largo as f64 / total.max(1) as f64;
    narracion.agregar(
        Pizarra::nueva()
            .variable("referencia", &referencia)
            .variable("k", indice.k)
            .variable("w", indice.w)
            .variable("semillas", miles(total))
            .secuencia(c.genoma, &[(0, indice.k + indice.w - 1, Resalte::Comparado)])
            .paso(
                "let mut semillas = minimizadores(referencia, k, w);",
                format!(
                    "Se hace con todo el genoma: salen {} minimizadores, uno cada {} bases más o menos, en lugar de los {} k-meros posibles.",
                    miles(total),
                    decimal(cada, 1),
                    miles(largo.saturating_sub(indice.k - 1))
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("referencia", &referencia)
            .variable("k", indice.k)
            .variable("w", indice.w)
            .variable("semillas", miles(total))
            .paso(
                "semillas.par_sort_unstable_by_key(|semilla| (semilla.valor, semilla.posicion));",
                "Se ordenan por huella. Así las huellas iguales quedan juntas y se pueden contar y buscar rápido, como las palabras de un diccionario.",
            ),
    );
    let descarte = if indice.descartados == 0 {
        "En esta referencia ninguna huella se repite tanto, así que no se tira ninguna.".to_string()
    } else {
        format!(
            "Aquí se descartaron {} huellas distintas ({} apariciones): son repeticiones del genoma.",
            miles(indice.descartados),
            miles(indice.apariciones_descartadas)
        )
    };
    narracion.agregar(
        Pizarra::nueva()
            .variable("frecuencia_maxima", indice.frecuencia_maxima)
            .variable("indice.descartados", indice.descartados)
            .variable("indice.semillas", miles(guardados))
            .paso(
                "if apariciones <= indice.frecuencia_maxima {",
                format!(
                    "Una huella que aparece más de {} veces no sirve para ubicar una lectura: cabría en demasiados lugares. {descarte} Quedan {} minimizadores en el índice.",
                    indice.frecuencia_maxima,
                    miles(guardados)
                ),
            ),
    );
    let cajones = 1usize << indice.bits_de_cajon;
    let por_cajon = if guardados / cajones.max(1) == 0 {
        "menos de una semilla por cajón en promedio".to_string()
    } else {
        format!("unas {} semillas por cajón", miles(guardados / cajones))
    };
    narracion.agregar(
        Pizarra::nueva()
            .variable("indice.semillas", miles(guardados))
            .variable("indice.bits_de_cajon", indice.bits_de_cajon)
            .variable("indice.cajones", miles(indice.cajones.len()))
            .paso(
                "indice.bits_de_cajon = (usize::BITS - indice.semillas.len().leading_zeros()).clamp(1, 24);",
                format!(
                    "Para no buscar entre {} semillas, se reparten en {} cajones según los primeros {} bits de su huella: así tocan {por_cajon} y buscar es casi inmediato.",
                    miles(guardados),
                    miles(cajones),
                    indice.bits_de_cajon
                ),
            ),
    );
}

fn base_de_lectura(lectura: &Lectura) -> Pizarra {
    Pizarra::nueva().variable(
        "lectura",
        format!(
            "{} ({} pb)",
            nombre_de(lectura),
            miles(lectura.secuencia.len())
        ),
    )
}

fn tabla_de_anclas(
    anclas: &[Ancla],
    diagonal_comun: Option<i64>,
) -> (Vec<Vec<String>>, Vec<(usize, usize)>) {
    let mut filas = Vec::new();
    let mut resaltar = Vec::new();
    for (fila, ancla) in anclas.iter().take(ANCLAS_A_MOSTRAR).enumerate() {
        filas.push(vec![
            hebra(ancla.reverso).to_string(),
            miles(ancla.en_lectura + 1),
            miles(ancla.en_referencia + 1),
            diagonal_de(ancla).to_string(),
        ]);
        if Some(diagonal_de(ancla)) == diagonal_comun {
            resaltar.push((fila, 3));
        }
    }
    (filas, resaltar)
}

fn diagonal_mas_comun(anclas: &[Ancla]) -> Option<(i64, usize)> {
    let mut cuentas: HashMap<i64, usize> = HashMap::new();
    for ancla in anclas {
        *cuentas.entry(diagonal_de(ancla)).or_insert(0) += 1;
    }
    cuentas
        .into_iter()
        .max_by_key(|&(diagonal, veces)| (veces, std::cmp::Reverse(diagonal)))
}

fn narrar_anclas(c: &Contexto, lectura: &Lectura, a: &Analisis, narracion: &mut Narracion) {
    let indice = c.indice;
    let secuencia = &lectura.secuencia;
    let marcas: Vec<(usize, usize, Resalte)> = a
        .minimizadores
        .iter()
        .map(|minimizador| {
            (
                minimizador.posicion,
                minimizador.posicion + indice.k,
                Resalte::Hallado,
            )
        })
        .collect();
    let nota = if a.minimizadores.is_empty() {
        "La lectura es tan corta o tiene tantas N que no da ningún minimizador, así que no habrá anclas.".to_string()
    } else {
        format!(
            "Para ubicar la lectura no se usan todas sus {} letras: se sacan sus {} minimizadores (marcados), con el mismo método de antes.",
            miles(secuencia.len()),
            miles(a.minimizadores.len())
        )
    };
    narracion.agregar(
        base_de_lectura(lectura)
            .variable("minimizadores", a.minimizadores.len())
            .secuencia(secuencia, &marcas)
            .lista(
                "Primeros minimizadores de la lectura",
                a.minimizadores
                    .iter()
                    .take(POSICIONES_A_MOSTRAR)
                    .map(|minimizador| {
                        format!(
                            "posición {}: huella {}",
                            minimizador.posicion + 1,
                            huella(minimizador.valor)
                        )
                    })
                    .collect(),
                &[],
            )
            .paso(
                "for minimizador in minimizadores(lectura, indice.k, indice.w) {",
                nota,
            ),
    );
    let elegido = a
        .minimizadores
        .iter()
        .find(|minimizador| !indice.buscar(minimizador.valor).is_empty())
        .or(a.minimizadores.first());
    if let Some(&minimizador) = elegido {
        let encontradas = indice.buscar(minimizador.valor);
        let cajon = (minimizador.valor >> (64 - indice.bits_de_cajon)) as usize;
        let del_cajon = indice.cajones[cajon + 1] - indice.cajones[cajon];
        let mut paso = base_de_lectura(lectura)
            .variable("minimizador.valor", huella(minimizador.valor))
            .variable("cajon", cajon)
            .variable("del_cajon.len()", del_cajon)
            .variable("fin - inicio", encontradas.len());
        let nota = if encontradas.is_empty() {
            format!(
                "Se toma la huella {} del minimizador de la posición {}. Cae en el cajón {cajon} y no hay ninguna coincidencia: esa palabra no está en el genoma o era una repetición descartada.",
                huella(minimizador.valor),
                minimizador.posicion + 1
            )
        } else {
            paso = paso.secuencia(
                c.genoma,
                &[(
                    encontradas[0].posicion,
                    encontradas[0].posicion + indice.k,
                    Resalte::Hallado,
                )],
            );
            format!(
                "Se busca la huella {} del minimizador de la posición {}. Su cajón es el {cajon}, con {} semillas; una búsqueda binaria halla {} exacta(s) en el genoma. Cada una dará una ancla candidata.",
                huella(minimizador.valor),
                minimizador.posicion + 1,
                miles(del_cajon),
                contar(encontradas.len(), "coincidencia", "coincidencias")
            )
        };
        narracion.agregar(
            paso.lista(
                "Dónde está en el genoma",
                encontradas
                    .iter()
                    .take(POSICIONES_A_MOSTRAR)
                    .map(|semilla| format!("posición {}", miles(semilla.posicion + 1)))
                    .collect(),
                &[0],
            )
            .paso("for semilla in indice.buscar(minimizador.valor) {", nota),
        );
        if let Some(semilla) = encontradas.first() {
            let reverso = minimizador.reverso != semilla.reverso;
            let en_lectura = if reverso {
                secuencia.len() - indice.k - minimizador.posicion
            } else {
                minimizador.posicion
            };
            let explicacion = if reverso {
                format!(
                    "Se leyeron en direcciones distintas, así que la lectura viene de la hebra reversa y su posición se cuenta desde el otro extremo: {} − {} − {} = {}.",
                    miles(secuencia.len()),
                    indice.k,
                    minimizador.posicion,
                    en_lectura
                )
            } else {
                "Se leyeron en la misma dirección: la lectura viene de la hebra directa y la posición queda igual.".to_string()
            };
            narracion.agregar(
                base_de_lectura(lectura)
                    .variable("minimizador.reverso", minimizador.reverso)
                    .variable("semilla.reverso", semilla.reverso)
                    .variable("reverso", reverso)
                    .paso(
                        "let reverso = minimizador.reverso != semilla.reverso;",
                        format!(
                            "El minimizador de la lectura y la semilla del genoma tienen la misma huella; ahora se comprueba en qué dirección estaba escrita cada una. {explicacion}"
                        ),
                    ),
            );
            narracion.agregar(
                base_de_lectura(lectura)
                    .variable("reverso", reverso)
                    .variable("en_referencia", semilla.posicion)
                    .variable("en_lectura", en_lectura)
                    .secuencia(
                        c.genoma,
                        &[(semilla.posicion, semilla.posicion + indice.k, Resalte::Hallado)],
                    )
                    .paso(
                        "anclas.push(Ancla {",
                        format!(
                            "Una ancla es una pareja: la posición {} de la lectura equivale a la posición {} del genoma. Se guarda y se sigue con el siguiente minimizador.",
                            en_lectura + 1,
                            miles(semilla.posicion + 1)
                        ),
                    ),
            );
        }
    }
    let comun = diagonal_mas_comun(&a.anclas);
    let (filas, resaltar) = tabla_de_anclas(&a.anclas, comun.map(|(diagonal, _)| diagonal));
    let nota = match comun {
        Some((diagonal, veces)) => format!(
            "Repetido con todos los minimizadores, salen {} anclas. {} comparten la misma diagonal ({}), es decir, la misma resta posición en el genoma − posición en la lectura: ahí es donde cae la lectura.",
            miles(a.anclas.len()),
            if veces == a.anclas.len() {
                "Todas".to_string()
            } else {
                format!("{} de ellas", miles(veces))
            },
            miles_con_signo(diagonal)
        ),
        None => "Ningún minimizador de la lectura está en el índice: no hay anclas.".to_string(),
    };
    narracion.agregar(
        base_de_lectura(lectura)
            .variable("anclas", miles(a.anclas.len()))
            .tabla(
                "Primeras anclas",
                &["hebra", "en lectura", "en genoma", "diagonal"],
                filas,
                &resaltar,
            )
            .paso("let anclas = buscar_anclas(lectura, indice);", nota),
    );
    if a.anclas.is_empty() {
        return;
    }
    let mut ordenadas = a.anclas.clone();
    ordenadas.sort_unstable_by_key(|ancla| (ancla.reverso, ancla.en_referencia, ancla.en_lectura));
    let (filas, resaltar) = tabla_de_anclas(&ordenadas, comun.map(|(diagonal, _)| diagonal));
    narracion.agregar(
        base_de_lectura(lectura)
            .variable("anclas", miles(a.anclas.len()))
            .tabla(
                "Anclas ordenadas",
                &["hebra", "en lectura", "en genoma", "diagonal"],
                filas,
                &resaltar,
            )
            .paso(
                "anclas.sort_unstable_by_key(|ancla| (ancla.reverso, ancla.en_referencia, ancla.en_lectura));",
                "Para encadenarlas hay que ponerlas en orden: primero por hebra y luego por su posición en el genoma. Las que vienen de la misma lectura quedan en fila.",
            ),
    );
}

fn miles_con_signo(numero: i64) -> String {
    if numero < 0 {
        format!("−{}", miles(numero.unsigned_abs() as usize))
    } else {
        miles(numero as usize)
    }
}

fn narrar_encadenado(c: &Contexto, lectura: &Lectura, a: &Analisis, narracion: &mut Narracion) {
    let k = c.indice.k;
    if let Some(cadena) = a.todas.first()
        && cadena.anclas.len() >= 2
    {
        let mut propuesta = k as f64;
        for (numero, par) in cadena.anclas.windows(2).take(ENLACES_A_MOSTRAR).enumerate() {
            let (previa, actual) = (par[0], par[1]);
            let Some(ganancia) = mapeo::ganancia_de_enlace(previa, actual, k) else {
                break;
            };
            propuesta += ganancia;
            let avance_en_referencia = actual.en_referencia - previa.en_referencia;
            let avance_en_lectura = actual.en_lectura - previa.en_lectura;
            let desvio = avance_en_referencia.abs_diff(avance_en_lectura);
            let nuevas = avance_en_referencia.min(avance_en_lectura).min(k);
            let (linea, nota) = if numero == 0 {
                (
                    "let desvio = avance_en_referencia.abs_diff(avance_en_lectura);",
                    format!(
                        "Dos anclas seguidas de la mejor cadena: la lectura avanzó {avance_en_lectura} letras y el genoma {avance_en_referencia}, así que el desvío es {desvio}. La ganancia es {nuevas} letras nuevas menos {} de costo por el desvío: {}. Como cada ancla ya vale k = {k}, la cadena suma {}.",
                        decimal(nuevas as f64 - ganancia, 2),
                        decimal(ganancia, 2),
                        decimal(propuesta, 1)
                    ),
                )
            } else {
                (
                    "let propuesta = puntaje[previa] + ganancia;",
                    format!(
                        "Siguiente enlace: la lectura avanza {avance_en_lectura} y el genoma {avance_en_referencia} (desvío {desvio}). Gana {} y el puntaje acumulado sube a {}.",
                        decimal(ganancia, 2),
                        decimal(propuesta, 1)
                    ),
                )
            };
            narracion.agregar(
                base_de_lectura(lectura)
                    .variable("avance_en_referencia", avance_en_referencia)
                    .variable("avance_en_lectura", avance_en_lectura)
                    .variable("desvio", desvio)
                    .variable("ganancia", decimal(ganancia, 2))
                    .variable("propuesta", decimal(propuesta, 1))
                    .secuencia(
                        c.genoma,
                        &[
                            (
                                previa.en_referencia,
                                previa.en_referencia + k,
                                Resalte::Comparado,
                            ),
                            (
                                actual.en_referencia,
                                actual.en_referencia + k,
                                Resalte::Actual,
                            ),
                        ],
                    )
                    .paso(linea, nota),
            );
        }
    }
    if a.todas.is_empty() {
        return;
    }
    let mejor = &a.todas[0];
    let mut filas = Vec::new();
    for cadena in a.todas.iter().take(CADENAS_A_MOSTRAR) {
        filas.push(vec![
            hebra(cadena.reverso).to_string(),
            miles(cadena.anclas.len()),
            decimal(cadena.puntaje, 1),
            descripcion_de_cadena(cadena),
        ]);
    }
    narracion.agregar(
        base_de_lectura(lectura)
            .variable("anclas", miles(a.anclas.len()))
            .variable("cadenas", miles(a.todas.len()))
            .tabla(
                "Cadenas armadas (de mayor a menor puntaje)",
                &["hebra", "anclas", "puntaje", "posiciones en el genoma"],
                filas,
                &[(0, 2)],
            )
            .paso(
                "separar_cadenas(&anclas, &puntaje, &anterior)",
                format!(
                    "Cada ancla recuerda cuál es su mejor vecina anterior. Siguiendo esas flechas hacia atrás desde la ancla de mayor puntaje se arma la mejor cadena, y las anclas ya usadas no se repiten. Salen {}; la mejor tiene {} anclas y {} puntos.",
                    contar(a.todas.len(), "cadena", "cadenas"),
                    miles(mejor.anclas.len()),
                    decimal(mejor.puntaje, 1)
                ),
            ),
    );
}

fn descripcion_de_cadena(cadena: &Cadena) -> String {
    match (cadena.anclas.first(), cadena.anclas.last()) {
        (Some(primera), Some(ultima)) => format!(
            "{}–{}",
            miles(primera.en_referencia + 1),
            miles(ultima.en_referencia + 1)
        ),
        _ => "—".to_string(),
    }
}

fn narrar_candidatas(lectura: &Lectura, a: &Analisis, narracion: &mut Narracion) {
    let filas: Vec<Vec<String>> = a
        .candidatas
        .iter()
        .take(CADENAS_A_MOSTRAR)
        .map(|cadena| {
            vec![
                hebra(cadena.reverso).to_string(),
                miles(cadena.anclas.len()),
                decimal(cadena.puntaje, 1),
                descripcion_de_cadena(cadena),
            ]
        })
        .collect();
    let nota = if a.candidatas.is_empty() {
        format!(
            "Ninguna cadena cumple: hace falta que tenga al menos {} anclas y {} puntos. Sin candidatas no hay dónde alinear.",
            mapeo::MINIMO_DE_ANCLAS,
            mapeo::PUNTAJE_MINIMO_DE_CADENA
        )
    } else {
        format!(
            "Se descartan las cadenas con menos de {} anclas o menos de {} puntos, las que valen menos de la mitad de la mejor y las que se superponen con una ya elegida. Quedan {}.",
            mapeo::MINIMO_DE_ANCLAS,
            mapeo::PUNTAJE_MINIMO_DE_CADENA,
            contar(a.candidatas.len(), "candidata", "candidatas")
        )
    };
    narracion.agregar(
        base_de_lectura(lectura)
            .variable("anclas", miles(a.anclas.len()))
            .variable("cadenas", a.candidatas.len())
            .tabla(
                "Cadenas candidatas",
                &["hebra", "anclas", "puntaje", "posiciones en el genoma"],
                filas,
                &[],
            )
            .paso(
                "let cadenas = elegir_candidatas(encadenar(&anclas, indice.k));",
                nota,
            ),
    );
    if a.candidatas.is_empty() {
        narracion.agregar(base_de_lectura(lectura).variable("cadenas", 0).paso(
            "if cadenas.is_empty() {",
            "No hay cadenas candidatas: la lectura termina aquí, sin mapear.",
        ));
    }
}

fn tabla_de_operaciones(cigar: &[(u32, u8)]) -> Vec<Vec<String>> {
    cigar
        .iter()
        .take(6)
        .map(|&(cantidad, operacion)| {
            let significado = match operacion {
                b'M' => "alineadas con el genoma (coincidan o no)",
                b'S' => "recortadas: no se parecen al genoma",
                b'I' => "sobran en la lectura",
                _ => "faltan en la lectura",
            };
            vec![
                (operacion as char).to_string(),
                miles(cantidad as usize),
                significado.to_string(),
            ]
        })
        .collect()
}

fn narrar_alineamientos(c: &Contexto, lectura: &Lectura, a: &Analisis, narracion: &mut Narracion) {
    let cantidad = a.intentos.len();
    for (numero, intento) in a.intentos.iter().take(CADENAS_DETALLADAS).enumerate() {
        let cadena = &intento.cadena;
        let etiqueta = format!("{} de {}", numero + 1, cantidad);
        let (menor, mayor) = intento.diagonales;
        let largo = lectura.secuencia.len();
        let marcas_de_cadena: Vec<(usize, usize, Resalte)> =
            match (cadena.anclas.first(), cadena.anclas.last()) {
                (Some(primera), Some(ultima)) => vec![
                    (
                        primera.en_referencia,
                        primera.en_referencia + c.indice.k,
                        Resalte::Hallado,
                    ),
                    (
                        ultima.en_referencia,
                        ultima.en_referencia + c.indice.k,
                        Resalte::Hallado,
                    ),
                ],
                _ => Vec::new(),
            };
        if cadena.reverso {
            narracion.agregar(
                base_de_lectura(lectura)
                    .variable("cadena", &etiqueta)
                    .variable("cadena.reverso", true)
                    .secuencia(&complemento_inverso(&lectura.secuencia), &[])
                    .paso(
                        "let orientada = if cadena.reverso { &inversa } else { lectura };",
                        "Las anclas de esta cadena están en la hebra reversa: la lectura es el complemento inverso del genoma. Se le da la vuelta (se invierte el orden y A↔T, C↔G) y esa es la que se alinea.",
                    ),
            );
        }
        let ancho = (mayor - menor + 2 * mapeo::MARGEN_DE_BANDA + 1).max(0) as usize;
        narracion.agregar(
            base_de_lectura(lectura)
                .variable("cadena", &etiqueta)
                .variable("menor", miles_con_signo(menor))
                .variable("mayor", miles_con_signo(mayor))
                .secuencia(c.genoma, &marcas_de_cadena)
                .paso(
                    "let (menor, mayor) = cadena.diagonales();",
                    format!(
                        "Las anclas de la cadena dicen que la lectura empieza en el genoma entre las diagonales {} y {}. Solo se alinea alrededor de esa franja (con {} de margen a cada lado, {} columnas en total) en vez de contra todo el genoma.",
                        miles_con_signo(menor),
                        miles_con_signo(mayor),
                        mapeo::MARGEN_DE_BANDA,
                        ancho
                    ),
                ),
        );
        let Some(alineamiento) = &intento.alineamiento else {
            narracion.agregar(
                base_de_lectura(lectura)
                    .variable("cadena", &etiqueta)
                    .paso(
                        "let Some(alineamiento) = alinear_con_recorte(orientada, referencia, menor, mayor) else {",
                        "No salió un alineamiento válido en esa franja (la lectura se sale del genoma o no se parece lo bastante): esta cadena se ignora.",
                    ),
            );
            continue;
        };
        let marcas_de_alineamiento = [(
            alineamiento.posicion,
            alineamiento.posicion + alineamiento.largo_en_referencia(),
            Resalte::Hallado,
        )];
        let directo = menor == mayor
            && alineamiento.cigar == [(largo as u32, b'M')]
            && alineamiento.diferencias <= 1;
        if directo {
            narracion.agregar(
                base_de_lectura(lectura)
                    .variable("cadena", &etiqueta)
                    .variable("alineamiento.puntaje", alineamiento.puntaje)
                    .variable("alineamiento.diferencias", alineamiento.diferencias)
                    .secuencia(c.genoma, &marcas_de_alineamiento)
                    .paso(
                        "return Some(directo);",
                        format!(
                            "Hay una sola diagonal y la lectura casi no cambia, así que basta comparar base contra base: {} coincidencias y {} diferencia(s). Puntaje = {} × {} − {} × {} = {}.",
                            miles(alineamiento.coincidencias),
                            alineamiento.diferencias,
                            miles(alineamiento.coincidencias),
                            mapeo::COINCIDENCIA,
                            alineamiento.diferencias,
                            mapeo::DIFERENCIA,
                            alineamiento.puntaje
                        ),
                    ),
            );
        } else {
            narracion.agregar(
                base_de_lectura(lectura)
                    .variable("cadena", &etiqueta)
                    .variable("filas", largo + 1)
                    .variable("ancho", ancho)
                    .secuencia(c.genoma, &marcas_de_cadena)
                    .paso(
                        "llenar_tabla(&mut tabla, lectura, referencia);",
                        format!(
                            "Se arma una tabla de {} filas (una por base de la lectura) y {ancho} columnas (la franja). Cada celda guarda el mejor puntaje para llegar ahí: +{} por coincidencia, −{} por diferencia, −{} por abrir un hueco y −{} por cada base más del hueco.",
                            miles(largo + 1),
                            mapeo::COINCIDENCIA,
                            mapeo::DIFERENCIA,
                            mapeo::ABRIR_HUECO,
                            mapeo::EXTENDER_HUECO
                        ),
                    ),
            );
            narracion.agregar(
                base_de_lectura(lectura)
                    .variable("cadena", &etiqueta)
                    .variable("alineamiento.puntaje", alineamiento.puntaje)
                    .secuencia(c.genoma, &marcas_de_alineamiento)
                    .tabla(
                        "Operaciones del alineamiento (CIGAR)",
                        &["operación", "bases", "significado"],
                        tabla_de_operaciones(&alineamiento.cigar),
                        &[],
                    )
                    .paso(
                        "reconstruir(&tabla, lectura, referencia)",
                        format!(
                            "Se parte de la mejor celda y se vuelve hacia atrás siguiendo de dónde vino cada puntaje. El camino, leído al derecho, es el CIGAR {}: empieza en la posición {} del genoma y suma {} puntos (los extremos que no se parecen se recortan, S, y llegar a un extremo da +{}).",
                            cigar_en_texto(&alineamiento.cigar),
                            miles(alineamiento.posicion + 1),
                            alineamiento.puntaje,
                            mapeo::PREMIO_POR_LLEGAR_AL_EXTREMO
                        ),
                    ),
            );
        }
        let nota = if intento.aceptado {
            format!(
                "El puntaje {} es al menos {}: este alineamiento pasa a ser candidato.",
                alineamiento.puntaje,
                mapeo::PUNTAJE_MINIMO_DE_ALINEAMIENTO
            )
        } else {
            format!(
                "El puntaje {} no llega a {}: es demasiado flojo y se descarta.",
                alineamiento.puntaje,
                mapeo::PUNTAJE_MINIMO_DE_ALINEAMIENTO
            )
        };
        narracion.agregar(
            base_de_lectura(lectura)
                .variable("cadena", &etiqueta)
                .variable("alineamiento.puntaje", alineamiento.puntaje)
                .variable("alineamiento.coincidencias", alineamiento.coincidencias)
                .variable("alineamiento.diferencias", alineamiento.diferencias)
                .variable("alineamiento.insertadas", alineamiento.insertadas)
                .variable("alineamiento.borradas", alineamiento.borradas)
                .paso(
                    "if alineamiento.puntaje >= PUNTAJE_MINIMO_DE_ALINEAMIENTO {",
                    nota,
                ),
        );
        if let Some(efecto) = intento.efecto {
            let texto = match efecto {
                Efecto::Nuevo => "Es un lugar nuevo: se añade a los candidatos.",
                Efecto::Reemplaza => {
                    "Cae sobre un candidato que ya había, pero con mejor puntaje: lo reemplaza."
                }
                Efecto::Repetido => {
                    "Cae sobre un candidato que ya había y no lo mejora: no se repite."
                }
            };
            narracion.agregar(
                base_de_lectura(lectura)
                    .variable("cadena", &etiqueta)
                    .variable("alineamiento.puntaje", alineamiento.puntaje)
                    .paso("agregar_sin_repetir(&mut candidatos, nuevo);", texto),
            );
        }
    }
}

fn tabla_de_candidatos(a: &Analisis) -> Vec<Vec<String>> {
    a.candidatos
        .iter()
        .take(CADENAS_A_MOSTRAR)
        .map(|candidato| {
            vec![
                hebra(candidato.reverso).to_string(),
                miles(candidato.alineamiento.posicion + 1),
                cigar_en_texto(&candidato.alineamiento.cigar),
                candidato.alineamiento.puntaje.to_string(),
                miles(candidato.anclas),
            ]
        })
        .collect()
}

fn narrar_eleccion(lectura: &Lectura, a: &Analisis, detallada: bool, narracion: &mut Narracion) {
    let titulo = "Candidatos a posición";
    let columnas = ["hebra", "posición", "CIGAR", "puntaje", "anclas"];
    let omitidas = a.intentos.len().saturating_sub(CADENAS_DETALLADAS);
    let Some(mapeo) = &a.mapeo else {
        let nota = if a.candidatos.is_empty() {
            format!(
                "Ningún alineamiento llegó a {} puntos, así que no hay candidatos y la lectura queda sin mapear.",
                mapeo::PUNTAJE_MINIMO_DE_ALINEAMIENTO
            )
        } else {
            "No quedó ninguna posición posible: la lectura queda sin mapear.".to_string()
        };
        narracion.agregar(
            base_de_lectura(lectura)
                .variable("candidatos", a.candidatos.len())
                .paso("let mejor = candidatos", nota),
        );
        narracion.agregar(
            base_de_lectura(lectura)
                .variable("candidatos", a.candidatos.len())
                .paso(
                    "bandera: 4,",
                    "En el SAM se escribe con bandera 4 (sin mapear) y sin posición. Más adelante solo suma a «lecturas».",
                ),
        );
        return;
    };
    let resto = if detallada && omitidas > 0 {
        format!(" Las otras {omitidas} cadenas se procesaron igual.")
    } else {
        String::new()
    };
    narracion.agregar(
        base_de_lectura(lectura)
            .variable("candidatos", a.candidatos.len())
            .variable("mejor", a.mejor)
            .tabla(
                titulo,
                &columnas,
                tabla_de_candidatos(a),
                &[(a.elegido.min(CADENAS_A_MOSTRAR - 1), 3)],
            )
            .paso(
                "let mejor = candidatos",
                format!(
                    "Hay {}; el mejor puntaje es {}.{resto}",
                    contar(a.candidatos.len(), "candidato", "candidatos"),
                    a.mejor
                ),
            ),
    );
    if detallada && a.empatados > 1 {
        narracion.agregar(
            base_de_lectura(lectura)
                .variable("mejor", a.mejor)
                .variable("empatados.len()", a.empatados)
                .paso(
                    "let elegido = empatados[(desempate % empatados.len() as u64) as usize];",
                    format!(
                        "Hay {} candidatos empatados con {} puntos. Se elige uno con un número sacado del nombre de la lectura: no es azar de verdad, así el resultado siempre es el mismo.",
                        a.empatados, a.mejor
                    ),
                ),
        );
    }
    let ganador_anclas = a
        .candidatos
        .get(a.elegido)
        .map(|candidato| candidato.anclas)
        .unwrap_or(0);
    let segundo = a.segundo.clamp(0, a.mejor);
    let unicidad = 1.0 - segundo as f64 / a.mejor as f64;
    let apoyo = (ganador_anclas as f64 / 10.0).min(1.0);
    let bruta = 40.0 * unicidad * apoyo * (a.mejor as f64).ln();
    if detallada {
        narracion.agregar(
            base_de_lectura(lectura)
                .variable("mejor", a.mejor)
                .variable("segundo", a.segundo)
                .paso(
                    "segundo = segundo.max(candidato.alineamiento.puntaje);",
                    if a.candidatos.len() > 1 {
                        format!(
                            "El mejor de los demás candidatos tiene {} puntos. Si se parece al mejor, la lectura cabe en dos lugares y es ambigua.",
                            a.segundo
                        )
                    } else {
                        "No hay otros candidatos, así que el segundo puntaje es 0: no hay un segundo lugar posible.".to_string()
                    },
                ),
        );
        narracion.agregar(
            base_de_lectura(lectura)
                .variable("mejor", a.mejor)
                .variable("segundo", a.segundo)
                .variable("unicidad", decimal(unicidad, 2))
                .paso(
                    "let unicidad = 1.0 - segundo as f64 / mejor as f64;",
                    format!(
                        "Unicidad = 1 − {segundo} / {} = {}. Vale 1 si no hay rival y 0 si hay otro lugar igual de bueno.",
                        a.mejor,
                        decimal(unicidad, 2)
                    ),
                ),
        );
        narracion.agregar(
            base_de_lectura(lectura)
                .variable("anclas", ganador_anclas)
                .variable("unicidad", decimal(unicidad, 2))
                .variable("apoyo", decimal(apoyo, 2))
                .paso(
                    "let apoyo = (anclas as f64 / 10.0).min(1.0);",
                    format!(
                        "Apoyo = {ganador_anclas} anclas / 10, con tope en 1: {}. Más anclas en la cadena ganadora dan más confianza.",
                        decimal(apoyo, 2)
                    ),
                ),
        );
    }
    narracion.agregar(
        base_de_lectura(lectura)
            .variable("mejor", a.mejor)
            .variable("segundo", a.segundo)
            .variable("anclas", ganador_anclas)
            .variable("unicidad", decimal(unicidad, 2))
            .variable("apoyo", decimal(apoyo, 2))
            .variable("calidad", mapeo.calidad_mapeo)
            .paso(
                if detallada {
                    "let calidad = 40.0 * unicidad * apoyo * (mejor as f64).ln();"
                } else {
                    "calidad_mapeo: calidad_de_mapeo(mejor, segundo, ganador.anclas),"
                },
                format!(
                    "MAPQ = 40 × {} × {} × ln({}) = {}; se limita a 0–{} y se redondea hacia abajo: {}. {}",
                    decimal(unicidad, 2),
                    decimal(apoyo, 2),
                    a.mejor,
                    decimal(bruta, 1),
                    CALIDAD_MAXIMA,
                    mapeo.calidad_mapeo,
                    if mapeo.calidad_mapeo == 0 {
                        "MAPQ 0 significa posición ambigua: no la uses para buscar variantes."
                    } else if mapeo.calidad_mapeo >= 40 {
                        "Es una posición muy segura."
                    } else {
                        "Es una posición razonable, pero no del todo segura."
                    }
                ),
            ),
    );
    narracion.agregar(
        base_de_lectura(lectura)
            .variable("mapeo.reverso", mapeo.reverso)
            .variable("posicion", mapeo.alineamiento.posicion + 1)
            .variable("calidad_mapeo", mapeo.calidad_mapeo)
            .variable("cigar", cigar_en_texto(&mapeo.alineamiento.cigar))
            .paso(
                "posicion: mapeo.alineamiento.posicion + 1,",
                format!(
                    "En el SAM queda: bandera {} ({}), posición {} (contando desde 1), MAPQ {} y CIGAR {}.",
                    if mapeo.reverso { 16 } else { 0 },
                    hebra(mapeo.reverso),
                    miles(mapeo.alineamiento.posicion + 1),
                    mapeo.calidad_mapeo,
                    cigar_en_texto(&mapeo.alineamiento.cigar)
                ),
            ),
    );
}

type CasoDeLectura = (&'static str, fn(&Option<Mapeo>) -> bool, Detalle);

#[derive(Clone, Copy, PartialEq)]
enum Detalle {
    Completo,
    ConAlineamiento,
    Resumen,
}

fn narrar_lectura(
    c: &Contexto,
    lectura: &Lectura,
    introduccion: &str,
    detalle: Detalle,
    narracion: &mut Narracion,
) {
    let detallada = detalle == Detalle::Completo;
    let analisis = analizar(lectura, c.indice, c.genoma);
    let largo = lectura.secuencia.len();
    if detallada || analisis.demasiado_larga {
        let nota = if analisis.demasiado_larga {
            format!(
                "{introduccion} Mide {} bases y este mapeador es para lecturas cortas (máximo {}): se deja sin mapear.",
                miles(largo),
                miles(LONGITUD_MAXIMA_DE_LECTURA)
            )
        } else {
            format!(
                "{introduccion} Mide {} bases; primero se revisa que no pase de {} (este mapeador es para lecturas cortas).",
                miles(largo),
                miles(LONGITUD_MAXIMA_DE_LECTURA)
            )
        };
        narracion.agregar(
            base_de_lectura(lectura)
                .secuencia(&lectura.secuencia, &[])
                .paso("if lectura.len() > LONGITUD_MAXIMA_DE_LECTURA {", nota),
        );
    }
    if analisis.demasiado_larga {
        narracion.agregar(base_de_lectura(lectura).paso(
            "bandera: 4,",
            "Se escribe en el SAM con bandera 4 (sin mapear).",
        ));
        return;
    }
    if detallada {
        narrar_anclas(c, lectura, &analisis, narracion);
        narrar_encadenado(c, lectura, &analisis, narracion);
        narrar_candidatas(lectura, &analisis, narracion);
        narrar_alineamientos(c, lectura, &analisis, narracion);
    } else {
        narracion.agregar(
            base_de_lectura(lectura)
                .variable("anclas", miles(analisis.anclas.len()))
                .variable("cadenas", analisis.candidatas.len())
                .paso(
                    "let cadenas = elegir_candidatas(encadenar(&anclas, indice.k));",
                    format!(
                        "{introduccion} Sus minimizadores dan {} anclas, que se encadenan y dejan {}.",
                        miles(analisis.anclas.len()),
                        contar(
                            analisis.candidatas.len(),
                            "cadena candidata",
                            "cadenas candidatas"
                        )
                    ),
                ),
        );
        if analisis.candidatas.is_empty() {
            narracion.agregar(
                base_de_lectura(lectura)
                    .variable("cadenas", 0)
                    .paso(
                        "if cadenas.is_empty() {",
                        "No quedó ninguna cadena: la lectura no se parece a nada del genoma y queda sin mapear.",
                    ),
            );
            narracion.agregar(base_de_lectura(lectura).paso(
                "bandera: 4,",
                "En el SAM se escribe con bandera 4 (sin mapear) y sin posición.",
            ));
            return;
        }
    }
    if analisis.candidatas.is_empty() {
        narracion.agregar(base_de_lectura(lectura).paso(
            "bandera: 4,",
            "En el SAM se escribe con bandera 4 (sin mapear) y sin posición.",
        ));
        return;
    }
    if detalle == Detalle::ConAlineamiento {
        narrar_alineamientos(c, lectura, &analisis, narracion);
    }
    narrar_eleccion(lectura, &analisis, detallada, narracion);
}

fn narrar_lectura_del_archivo(
    c: &Contexto,
    lecturas: &[Lectura],
    mapeos: &[Option<Mapeo>],
    narracion: &mut Narracion,
) {
    let bases: usize = lecturas.iter().map(|lectura| lectura.secuencia.len()).sum();
    let menor = lecturas
        .iter()
        .map(|lectura| lectura.secuencia.len())
        .min()
        .unwrap_or(0);
    let mayor = lecturas
        .iter()
        .map(|lectura| lectura.secuencia.len())
        .max()
        .unwrap_or(0);
    let medidas = if menor == mayor {
        format!("todas de {} bases", miles(menor))
    } else {
        format!(
            "de entre {} y {} bases (en promedio {})",
            miles(menor),
            miles(mayor),
            miles(bases / lecturas.len())
        )
    };
    narracion.agregar(
        Pizarra::nueva()
            .variable("tanda", format!("{} lecturas", miles(lecturas.len())))
            .secuencia(&lecturas[0].secuencia, &[])
            .paso(
                "let tanda = leer_tanda(&mut lector)?;",
                format!(
                    "Se lee el FASTQ, en tandas de {} lecturas. Hay {} lecturas, {medidas}. La primera se llama {}.",
                    miles(mapeo::LECTURAS_POR_TANDA),
                    miles(lecturas.len()),
                    nombre_de(&lecturas[0])
                ),
            ),
    );
    let mapeadas = mapeos.iter().flatten().count();
    narracion.agregar(
        Pizarra::nueva()
            .variable("tanda", format!("{} lecturas", miles(lecturas.len())))
            .variable("mapeos", format!("{} con posición", miles(mapeadas)))
            .paso(
                "let mapeos = mapear_en_paralelo(&tanda, indice, referencia);",
                format!(
                    "Cada lectura pasa por mapear_lectura, y todas a la vez, repartidas entre los {} núcleos de tu computadora. Salen {} con posición y {} sin mapear. Veamos la primera con todo detalle.",
                    crate::nucleos(),
                    miles(mapeadas),
                    miles(lecturas.len() - mapeadas)
                ),
            ),
    );
    narrar_lectura(
        c,
        &lecturas[0],
        "Esta es la primera lectura del archivo.",
        Detalle::Completo,
        narracion,
    );
    let mut usadas = vec![0];
    let casos: [CasoDeLectura; 4] = [
        (
            "Esta es la primera lectura con un extremo recortado o con un hueco, así que no basta comparar base contra base.",
            |mapeo| {
                mapeo.as_ref().is_some_and(|mapeo| {
                    let alineamiento = &mapeo.alineamiento;
                    alineamiento.tiene_recorte()
                        || alineamiento.insertadas + alineamiento.borradas > 0
                })
            },
            Detalle::ConAlineamiento,
        ),
        (
            "Esta es la primera lectura que cayó en la hebra reversa.",
            |mapeo| mapeo.as_ref().is_some_and(|mapeo| mapeo.reverso),
            Detalle::Resumen,
        ),
        (
            "Esta es la primera lectura ambigua (MAPQ 0).",
            |mapeo| mapeo.as_ref().is_some_and(|mapeo| mapeo.calidad_mapeo == 0),
            Detalle::Resumen,
        ),
        (
            "Esta es la primera lectura que no se pudo mapear.",
            |mapeo| mapeo.is_none(),
            Detalle::Resumen,
        ),
    ];
    for (introduccion, cumple, detalle) in casos {
        let Some(numero) = (0..lecturas.len())
            .find(|&numero| !usadas.contains(&numero) && cumple(&mapeos[numero]))
        else {
            continue;
        };
        usadas.push(numero);
        let texto = format!("{introduccion} Es la lectura {}.", miles(numero + 1));
        narrar_lectura(c, &lecturas[numero], &texto, detalle, narracion);
    }
}

fn variables_del_recuento(pizarra: Pizarra, resumen: &ResumenMapeo) -> Pizarra {
    pizarra
        .variable("resumen.lecturas", miles(resumen.lecturas))
        .variable("resumen.mapeadas", miles(resumen.mapeadas))
        .variable("resumen.en_hebra_reversa", miles(resumen.en_hebra_reversa))
        .variable("resumen.con_recorte", miles(resumen.con_recorte))
        .variable("resumen.calidades[0]", miles(resumen.calidades[0]))
}

fn narrar_recuento(
    c: &Contexto,
    lecturas: &[Lectura],
    mapeos: &[Option<Mapeo>],
    narracion: &mut Narracion,
) -> ResumenMapeo {
    let mut resumen = ResumenMapeo {
        calidades: vec![0; CALIDAD_MAXIMA as usize + 1],
        longitud_referencia: c.genoma.len(),
        largo_ventana: mapeo::largo_de_ventana(c.genoma.len()),
        minimizadores_en_indice: c.indice.semillas.len(),
        minimizadores_descartados: c.indice.descartados,
        ..ResumenMapeo::default()
    };
    let mut cambios_de_profundidad = vec![0i64; c.genoma.len() + 1];
    let total = lecturas.len();
    let cortes: Vec<usize> = (1..=CORTES_DEL_RECUENTO)
        .map(|numero| numero * total / (CORTES_DEL_RECUENTO + 1))
        .collect();
    let mut primera_mapeada: Option<(usize, usize)> = None;
    for (numero, (lectura, resultado)) in lecturas.iter().zip(mapeos).enumerate() {
        mapeo::anotar(
            &mut resumen,
            &mut cambios_de_profundidad,
            lectura,
            resultado.as_ref(),
        );
        if let Some(mapeo) = resultado
            && primera_mapeada.is_none()
        {
            primera_mapeada = Some((numero, mapeo.alineamiento.posicion));
        }
        if numero < LECTURAS_UNA_A_UNA {
            narracion.agregar(paso_de_recuento(
                c,
                numero,
                lectura,
                resultado.as_ref(),
                &resumen,
            ));
            if let Some(mapeo) = resultado
                && primera_mapeada.map(|(primera, _)| primera) == Some(numero)
            {
                narracion.agregar(paso_de_extremos(c, numero, lectura, mapeo, &resumen));
                narracion.agregar(paso_de_registro(c, numero, lectura, mapeo, &resumen));
            }
        } else if cortes.contains(&(numero + 1)) {
            narracion.agregar(variables_del_recuento(Pizarra::nueva(), &resumen).paso(
                "for (lectura, mapeo) in tanda.iter().zip(&mapeos) {",
                format!(
                    "Van {} de {} lecturas: {} mapeadas ({}), {} por la hebra reversa, {} con un extremo recortado y {} con MAPQ 0.",
                    miles(numero + 1),
                    miles(total),
                    miles(resumen.mapeadas),
                    porcentaje(fraccion(resumen.mapeadas, resumen.lecturas), 1),
                    miles(resumen.en_hebra_reversa),
                    miles(resumen.con_recorte),
                    miles(resumen.calidades[0])
                ),
            ));
        }
    }
    narrar_profundidad(
        c,
        &mut resumen,
        &cambios_de_profundidad,
        primera_mapeada,
        narracion,
    );
    resumen
}

fn paso_de_recuento(
    c: &Contexto,
    numero: usize,
    lectura: &Lectura,
    resultado: Option<&Mapeo>,
    resumen: &ResumenMapeo,
) -> Paso {
    let pizarra = variables_del_recuento(Pizarra::nueva(), resumen);
    let Some(mapeo) = resultado else {
        return pizarra.paso(
            "let Some(mapeo) = mapeo else {",
            format!(
                "Lectura {} ({}): sin mapear. Solo suma 1 a resumen.lecturas y termina aquí.",
                numero + 1,
                nombre_de(lectura)
            ),
        );
    };
    let alineamiento = &mapeo.alineamiento;
    let mut cambios = vec!["suma 1 a mapeadas".to_string()];
    if mapeo.reverso {
        cambios.push("suma 1 a en_hebra_reversa".to_string());
    }
    if alineamiento.tiene_recorte() {
        cambios.push("suma 1 a con_recorte (tiene S en el CIGAR)".to_string());
    }
    cambios.push(format!("suma 1 a calidades[{}]", mapeo.calidad_mapeo));
    pizarra
        .secuencia(
            c.genoma,
            &[(
                alineamiento.posicion,
                alineamiento.posicion + alineamiento.largo_en_referencia(),
                Resalte::Hallado,
            )],
        )
        .paso(
            "resumen.mapeadas += 1;",
            format!(
                "Lectura {} ({}): cayó en la posición {} (hebra {}) con MAPQ {} y CIGAR {}. Por eso {}; además anota sus {} coincidencias y {} diferencias.",
                numero + 1,
                nombre_de(lectura),
                miles(alineamiento.posicion + 1),
                hebra(mapeo.reverso),
                mapeo.calidad_mapeo,
                cigar_en_texto(&alineamiento.cigar),
                cambios.join(", "),
                miles(alineamiento.coincidencias),
                miles(alineamiento.diferencias)
            ),
        )
}

fn paso_de_registro(
    c: &Contexto,
    numero: usize,
    lectura: &Lectura,
    mapeo: &Mapeo,
    resumen: &ResumenMapeo,
) -> Paso {
    let alineamiento = &mapeo.alineamiento;
    let sentido = if mapeo.reverso {
        "Como se mapeó en la hebra reversa, bandera 16, se escribe su complemento inverso."
    } else {
        "Bandera 0: hebra directa."
    };
    variables_del_recuento(Pizarra::nueva(), resumen)
        .lista(
            "Línea del SAM",
            vec![
                format!("nombre: {}", nombre_de(lectura)),
                format!("bandera: {}", if mapeo.reverso { 16 } else { 0 }),
                format!("referencia: {}", c.nombre),
                format!("posición: {}", miles(alineamiento.posicion + 1)),
                format!("MAPQ: {}", mapeo.calidad_mapeo),
                format!("CIGAR: {}", cigar_en_texto(&alineamiento.cigar)),
                format!("secuencia: {}", corto(&lectura.secuencia)),
            ],
            &[],
        )
        .paso(
            "escribir_sam(&mut destino, &registro).map_err(error_de_escritura)?;",
            format!(
                "La lectura {} se escribe como una línea del archivo SAM, con sus columnas en este orden. {sentido} Todas las lecturas se guardan así, una línea cada una.",
                numero + 1
            ),
        )
}

fn paso_de_extremos(
    c: &Contexto,
    numero: usize,
    lectura: &Lectura,
    mapeo: &Mapeo,
    resumen: &ResumenMapeo,
) -> Paso {
    let alineamiento = &mapeo.alineamiento;
    let fin = alineamiento.posicion + alineamiento.largo_en_referencia();
    variables_del_recuento(Pizarra::nueva(), resumen)
        .variable("alineamiento.posicion", miles(alineamiento.posicion))
        .variable("fin", miles(fin))
        .secuencia(
            c.genoma,
            &[
                (alineamiento.posicion, alineamiento.posicion + 1, Resalte::Actual),
                (fin - 1, fin, Resalte::Comparado),
            ],
        )
        .paso(
            "cambios_de_profundidad[alineamiento.posicion] += 1;",
            format!(
                "Para medir cuántas lecturas cubren cada base sin repasar cada lectura base por base, solo se anotan sus extremos: la lectura {} ({}) suma 1 donde empieza (posición {}) y resta 1 justo después de donde termina (posición {}).",
                numero + 1,
                nombre_de(lectura),
                miles(alineamiento.posicion + 1),
                miles(fin + 1)
            ),
        )
}

fn narrar_profundidad(
    c: &Contexto,
    resumen: &mut ResumenMapeo,
    cambios: &[i64],
    primera_mapeada: Option<(usize, usize)>,
    narracion: &mut Narracion,
) {
    mapeo::medir_profundidad(resumen, cambios);
    if let Some((_, inicio)) = primera_mapeada {
        let desde = inicio.saturating_sub(2);
        let hasta = (inicio + 5).min(c.genoma.len());
        let mut actual: i64 = cambios[..desde].iter().sum();
        let antes = actual;
        let mut en_el_inicio = 0;
        let mut elementos = Vec::new();
        for (posicion, &cambio) in cambios.iter().enumerate().take(hasta).skip(desde) {
            actual += cambio;
            if posicion == inicio {
                en_el_inicio = actual;
            }
            elementos.push(format!(
                "posición {}: cambio {cambio:+} → profundidad {actual}",
                miles(posicion + 1)
            ));
        }
        narracion.agregar(
            Pizarra::nueva()
                .variable("actual", en_el_inicio)
                .variable("resumen.longitud_referencia", miles(resumen.longitud_referencia))
                .secuencia(
                    c.genoma,
                    &[(inicio, inicio + 1, Resalte::Actual)],
                )
                .lista("Profundidad base por base", elementos, &[(inicio - desde).min(hasta - desde - 1)])
                .paso(
                    "actual += cambio;",
                    format!(
                        "Se recorre el genoma de izquierda a derecha acumulando los cambios: donde empieza una lectura la suma sube y donde termina baja. Hasta la posición {} la suma lleva {antes}; al pasar por los cambios de la posición {} (donde empieza la lectura que vimos) queda en {en_el_inicio}: esas son las lecturas que cubren esa base.",
                        miles(desde),
                        miles(inicio + 1)
                    ),
                ),
        );
    }
    narracion.agregar(
        Pizarra::nueva()
            .variable("resumen.bases_cubiertas", miles(resumen.bases_cubiertas))
            .variable("resumen.longitud_referencia", miles(resumen.longitud_referencia))
            .paso(
                "resumen.bases_cubiertas += 1;",
                format!(
                    "Cada base con profundidad mayor que 0 cuenta como cubierta. En total hay {} de {} bases ({}).",
                    miles(resumen.bases_cubiertas),
                    miles(resumen.longitud_referencia),
                    porcentaje(
                        fraccion(resumen.bases_cubiertas, resumen.longitud_referencia),
                        1
                    )
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("resumen.bases_cubiertas", miles(resumen.bases_cubiertas))
            .variable("resumen.profundidad_media", decimal(resumen.profundidad_media, 2))
            .paso(
                "resumen.profundidad_media = suma as f64 / resumen.bases_cubiertas as f64;",
                format!(
                    "La profundidad media es la suma de la profundidad de todas las bases cubiertas dividida entre cuántas son: {}×. Solo cuenta donde hay lecturas.",
                    decimal(resumen.profundidad_media, 1)
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("resumen.largo_ventana", miles(resumen.largo_ventana))
            .variable(
                "resumen.profundidad_por_ventana",
                format!("{} ventanas", miles(resumen.profundidad_por_ventana.len())),
            )
            .paso(
                "resumen.profundidad_por_ventana = promedios_por_ventana(&profundidad, resumen.largo_ventana);",
                format!(
                    "Para la gráfica de profundidad se promedia cada ventana de {} pb: salen {} puntos.",
                    miles(resumen.largo_ventana),
                    miles(resumen.profundidad_por_ventana.len())
                ),
            ),
    );
    let filas: Vec<Vec<String>> = resumen
        .tramos
        .iter()
        .take(TRAMOS_A_MOSTRAR)
        .map(|tramo| {
            vec![
                miles(tramo.inicio + 1),
                miles(tramo.fin),
                miles(tramo.fin - tramo.inicio),
                format!("{}×", decimal(tramo.profundidad_media, 1)),
            ]
        })
        .collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("resumen.tramos", format!("{} tramos", miles(resumen.tramos.len())))
            .tabla(
                "Primeros tramos con lecturas",
                &["inicio", "fin", "largo (pb)", "profundidad media"],
                filas,
                &[],
            )
            .paso(
                "resumen.tramos = tramos_cubiertos(&profundidad);",
                format!(
                    "Un tramo es una racha seguida de bases con profundidad mayor que 0. Hay {}; los huecos entre tramos son zonas sin ninguna lectura.",
                    miles(resumen.tramos.len())
                ),
            ),
    );
}

fn paso_final(c: &Contexto, resumen: &ResumenMapeo) -> Paso {
    variables_del_recuento(Pizarra::nueva(), resumen)
        .variable("resumen.bases_cubiertas", miles(resumen.bases_cubiertas))
        .variable("resumen.profundidad_media", decimal(resumen.profundidad_media, 2))
        .paso(
            "Ok(resumen)",
            format!(
                "Resultado: {} de {} lecturas se ubicaron en {} ({}); la profundidad media donde hay lecturas es {}× y el {} de la referencia quedó cubierto. Abajo, en «Qué significa», verás qué dicen esos números sobre tu muestra.",
                miles(resumen.mapeadas),
                miles(resumen.lecturas),
                c.nombre,
                porcentaje(fraccion(resumen.mapeadas, resumen.lecturas), 1),
                decimal(resumen.profundidad_media, 1),
                porcentaje(fraccion(resumen.bases_cubiertas, resumen.longitud_referencia), 1)
            ),
        )
}

#[cfg(test)]
mod pruebas {
    use super::*;
    use nucleo::azar::Azar;
    use nucleo::formatos::{Lectura, RegistroSam, abrir, escribir_fastq, leer_sam};
    use nucleo::secuencia::complemento_inverso;
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::{BufRead, BufWriter, Write};

    fn carpeta(nombre: &str) -> PathBuf {
        let carpeta =
            std::env::temp_dir().join(format!("mapeo-catalogo-{}-{nombre}", std::process::id()));
        fs::create_dir_all(&carpeta).unwrap();
        carpeta
    }

    fn leer_registros(ruta: &Path) -> Vec<RegistroSam> {
        abrir(ruta)
            .unwrap()
            .lines()
            .filter_map(|linea| leer_sam(linea.unwrap().as_bytes()).unwrap())
            .collect()
    }

    fn ruta_del_sam(salida: &Salida) -> PathBuf {
        salida
            .bloques
            .iter()
            .find_map(|bloque| match bloque {
                Bloque::Archivo { ruta, .. } => Some(PathBuf::from(ruta)),
                _ => None,
            })
            .unwrap()
    }

    fn filas_de_la_tabla(salida: &Salida) -> Vec<Vec<String>> {
        salida
            .bloques
            .iter()
            .find_map(|bloque| match bloque {
                Bloque::Tabla { filas, .. } => Some(filas.clone()),
                _ => None,
            })
            .unwrap()
    }

    fn numero(texto: &str) -> usize {
        texto.replace('\u{a0}', "").parse().unwrap()
    }

    fn escribir_lecturas(ruta: &Path, secuencias: &[Vec<u8>]) {
        let mut destino = BufWriter::new(File::create(ruta).unwrap());
        for (numero, secuencia) in secuencias.iter().enumerate() {
            let lectura = Lectura {
                nombre: format!("l{numero}").into_bytes(),
                secuencia: secuencia.clone(),
                calidad: vec![b'I'; secuencia.len()],
            };
            escribir_fastq(&mut destino, &lectura).unwrap();
        }
        destino.flush().unwrap();
    }

    fn lecturas_cada(muestra: &[u8], paso: usize) -> Vec<Vec<u8>> {
        let mut lecturas = Vec::new();
        let mut inicio = 0;
        while inicio + 150 <= muestra.len() {
            let tramo = &muestra[inicio..inicio + 150];
            if lecturas.len() % 2 == 0 {
                lecturas.push(tramo.to_vec());
            } else {
                lecturas.push(complemento_inverso(tramo));
            }
            inicio += paso;
        }
        lecturas
    }

    #[test]
    fn el_ejemplo_coincide_con_los_alineamientos_verdaderos() {
        let id = "puerto.mapear_lecturas";
        let salida = crate::ejecutar(id, crate::ejemplo(id).unwrap()).unwrap();
        let verdad = Entrada::nueva(&[("verdad", "@ejemplo:alineamientos")])
            .ruta("verdad")
            .unwrap();
        let verdaderos: HashMap<Vec<u8>, RegistroSam> = leer_registros(&verdad)
            .into_iter()
            .map(|registro| (registro.nombre.clone(), registro))
            .collect();
        let obtenidos = leer_registros(&ruta_del_sam(&salida));
        assert_eq!(obtenidos.len(), verdaderos.len());
        let mut correctas = 0;
        let mut ambiguas = 0;
        let mut erradas = 0;
        for obtenido in &obtenidos {
            let esperado = &verdaderos[&obtenido.nombre];
            if obtenido.mapeado()
                && obtenido.hebra_reversa() == esperado.hebra_reversa()
                && obtenido.posicion.abs_diff(esperado.posicion) <= 5
            {
                correctas += 1;
            } else if obtenido.mapeado() && obtenido.calidad_mapeo == 0 {
                ambiguas += 1;
            } else {
                erradas += 1;
            }
        }
        let total = obtenidos.len() as f64;
        assert!(correctas as f64 >= 0.96 * total);
        assert!((correctas + ambiguas) as f64 >= 0.99 * total);
        assert!(erradas as f64 <= 0.001 * total);
        let texto = salida.interpretacion.join(" ");
        assert!(salida.interpretacion.len() >= 2 && salida.interpretacion.len() <= 6);
        assert!(texto.contains("la muestra es de esta referencia"));
        assert!(texto.contains("casi todas caen en 2 regiones"));
        assert!(texto.contains("regiones elegidas"));
        assert!(!texto.contains("posible deleción"));
        let largas: Vec<(usize, usize)> = filas_de_la_tabla(&salida)
            .iter()
            .map(|fila| (numero(&fila[0]), numero(&fila[1])))
            .filter(|(inicio, fin)| fin - inicio >= 10_000)
            .collect();
        assert_eq!(largas.len(), 2);
        assert!(largas[0].0.abs_diff(2_300_001) <= 50 && largas[0].1.abs_diff(2_400_000) <= 50);
        assert!(largas[1].0.abs_diff(3_900_001) <= 50 && largas[1].1.abs_diff(4_200_000) <= 50);
        assert!(texto.contains(&miles(largas[0].0)) && texto.contains(&miles(largas[1].1)));
        assert!(texto.contains("MAPQ 0"));
        assert!(texto.contains("Recortar lecturas"));
    }

    #[test]
    fn una_referencia_propia_muestra_su_delecion() {
        let carpeta = carpeta("propia");
        let referencia = ecoli::region(100_001, 160_000).unwrap().to_vec();
        let ruta_referencia = carpeta.join("referencia.fasta");
        fs::write(
            &ruta_referencia,
            format!(
                ">mi_cromosoma una prueba\n{}\n",
                String::from_utf8_lossy(&referencia)
            ),
        )
        .unwrap();
        let mut muestra = referencia[..20_000].to_vec();
        muestra.extend_from_slice(&referencia[22_000..]);
        let mut lecturas = lecturas_cada(&muestra, 5);
        let mapeables = lecturas.len();
        let mut azar = Azar::nuevo(7);
        for _ in 0..mapeables / 10 {
            lecturas.push((0..150).map(|_| b"ACGT"[azar.entero(4)]).collect());
        }
        let ruta_lecturas = carpeta.join("lecturas.fastq");
        escribir_lecturas(&ruta_lecturas, &lecturas);
        let entrada = Entrada::nueva(&[
            ("lecturas", ruta_lecturas.to_str().unwrap()),
            ("referencia", ruta_referencia.to_str().unwrap()),
        ]);
        let salida = ejecutar_mapeo(&entrada).unwrap();
        let ruta_sam = ruta_del_sam(&salida);
        assert_eq!(ruta_sam, carpeta.join("lecturas.mapeo.sam"));
        let cabecera: Vec<String> = abrir(&ruta_sam)
            .unwrap()
            .lines()
            .take(2)
            .map(|linea| linea.unwrap())
            .collect();
        assert_eq!(cabecera[1], "@SQ\tSN:mi_cromosoma\tLN:60000");
        let registros = leer_registros(&ruta_sam);
        assert_eq!(registros.len(), lecturas.len());
        let mapeadas = registros.iter().filter(|r| r.mapeado()).count();
        assert_eq!(mapeadas, mapeables);
        let texto = salida.interpretacion.join(" ");
        assert!(texto.contains("otra cepa"));
        assert!(texto.contains("1 tramo de 100 pb o más"));
        assert!(texto.contains("posible deleción"));
        assert!(!texto.contains("regiones elegidas"));
        let marcas = salida
            .bloques
            .iter()
            .find_map(|bloque| match bloque {
                Bloque::Linea { marcas, .. } => Some(marcas.clone()),
                _ => None,
            })
            .unwrap();
        assert_eq!(marcas.len(), 1);
        assert!((marcas[0].x - 21.0).abs() < 0.05);
        let regiones = filas_de_la_tabla(&salida);
        assert_eq!(regiones.len(), 2);
        assert_eq!(numero(&regiones[0][0]), 1);
        assert_eq!(numero(&regiones[1][1]), 60_000);
        assert!(numero(&regiones[0][1]).abs_diff(20_000) <= 10);
        assert!(numero(&regiones[1][0]).abs_diff(22_001) <= 10);
        fs::remove_dir_all(&carpeta).unwrap();
    }

    #[test]
    fn lecturas_ajenas_quedan_sin_mapear() {
        let carpeta = carpeta("ajenas");
        let referencia = carpeta.join("referencia.fa");
        let tramo = ecoli::region(500_001, 530_000).unwrap();
        fs::write(
            &referencia,
            format!(">trozo\n{}\n", String::from_utf8_lossy(tramo)),
        )
        .unwrap();
        let mut azar = Azar::nuevo(8);
        let ajenas: Vec<Vec<u8>> = (0..300)
            .map(|_| (0..150).map(|_| b"ACGT"[azar.entero(4)]).collect())
            .collect();
        let lecturas = carpeta.join("ajenas.fq");
        escribir_lecturas(&lecturas, &ajenas);
        let entrada = Entrada::nueva(&[
            ("lecturas", lecturas.to_str().unwrap()),
            ("referencia", referencia.to_str().unwrap()),
        ]);
        let salida = ejecutar_mapeo(&entrada).unwrap();
        assert_eq!(salida.interpretacion.len(), 2);
        assert!(
            salida.interpretacion[0].contains(&format!("Solo se ubicó el {}", porcentaje(0.0, 1)))
        );
        assert!(salida.interpretacion[1].contains("bandera 4"));
        let registros = leer_registros(&ruta_del_sam(&salida));
        assert_eq!(registros.len(), 300);
        assert!(registros.iter().all(|registro| registro.bandera == 4));
        fs::remove_dir_all(&carpeta).unwrap();
    }

    #[test]
    fn referencias_malas_dan_error_en_castellano() {
        let carpeta = carpeta("malas");
        let lecturas = carpeta.join("lecturas.fastq");
        escribir_lecturas(
            &lecturas,
            &lecturas_cada(ecoli::region(1, 2_000).unwrap(), 100),
        );
        let dos = carpeta.join("dos.fasta");
        fs::write(
            &dos,
            ">a\nACGTACGTACGTACGTACGTACGTACGT\n>b\nACGTACGTACGTACGTACGTAC\n",
        )
        .unwrap();
        let proteina = carpeta.join("proteina.fasta");
        fs::write(&proteina, ">p\nMKLVQEDFRSWYHHMKLVQEDFRSWYHH\n").unwrap();
        let corta = carpeta.join("corta.fasta");
        fs::write(&corta, ">c\nACGTACGT\n").unwrap();
        let vacia = carpeta.join("vacia.fasta");
        fs::write(&vacia, "").unwrap();
        let casos = [
            ("@nada", "No conozco"),
            ("@ecoli:1-100", "No conozco"),
            ("xyz", "No existe"),
            (carpeta.to_str().unwrap(), "no es un archivo"),
            (dos.to_str().unwrap(), "2 secuencias"),
            (proteina.to_str().unwrap(), "no parece ADN"),
            (corta.to_str().unwrap(), "al menos"),
            (vacia.to_str().unwrap(), "ninguna secuencia"),
        ];
        for (referencia, mensaje) in casos {
            let entrada = Entrada::nueva(&[
                ("lecturas", lecturas.to_str().unwrap()),
                ("referencia", referencia),
            ]);
            let error = ejecutar_mapeo(&entrada).err().unwrap();
            assert!(error.contains(mensaje), "{referencia}: {error}");
        }
        let sin_lecturas = carpeta.join("sin_lecturas.fastq");
        fs::write(&sin_lecturas, "").unwrap();
        let entrada = Entrada::nueva(&[
            ("lecturas", sin_lecturas.to_str().unwrap()),
            ("referencia", "@ecoli"),
        ]);
        assert!(
            ejecutar_mapeo(&entrada)
                .err()
                .unwrap()
                .contains("ninguna lectura")
        );
        assert!(!carpeta.join("sin_lecturas.mapeo.sam").exists());
        let sobrantes = fs::read_dir(&carpeta)
            .unwrap()
            .filter(|archivo| {
                archivo
                    .as_ref()
                    .unwrap()
                    .file_name()
                    .to_string_lossy()
                    .starts_with('.')
            })
            .count();
        assert_eq!(sobrantes, 0);
        fs::remove_dir_all(&carpeta).unwrap();
    }

    #[test]
    fn muestra_el_codigo_del_indice_y_del_mapeo() {
        let fragmentos = crate::codigo("puerto.mapear_lecturas").unwrap();
        assert_eq!(fragmentos[0].nombre, "mapeo::mapear_lecturas");
        let nombres: Vec<&str> = fragmentos.iter().map(|f| f.nombre.as_str()).collect();
        for esperado in [
            "mapeo::construir_indice",
            "mapeo::minimizadores",
            "mapeo::mapear_fastq",
            "mapeo::encadenar",
            "mapeo::alinear_con_recorte",
            "mapeo::calidad_de_mapeo",
        ] {
            assert!(nombres.contains(&esperado), "falta {esperado}");
        }
    }

    fn lecturas_del_ejemplo() -> (PathBuf, Vec<Lectura>) {
        let ruta = Entrada::nueva(&[("lecturas", "@ejemplo:lecturas")])
            .ruta("lecturas")
            .unwrap();
        let lecturas = leer_lecturas_a_narrar(&ruta).unwrap();
        (ruta, lecturas)
    }

    #[test]
    fn el_paso_a_paso_llega_al_mismo_resultado_que_el_mapeo() {
        let (ruta, lecturas) = lecturas_del_ejemplo();
        let indice = indice_de_ecoli().unwrap();
        let genoma = ecoli::genoma();
        let mapeos = mapeo::mapear_en_paralelo(&lecturas, indice, genoma);
        for (lectura, esperado) in lecturas.iter().zip(&mapeos).take(3000) {
            assert_eq!(&analizar(lectura, indice, genoma).mapeo, esperado);
        }
        let contexto = Contexto {
            genoma,
            nombre: ecoli::IDENTIFICADOR,
            indice,
        };
        let mut narracion = Narracion::nueva(400);
        let resumen = narrar_recuento(&contexto, &lecturas, &mapeos, &mut narracion);
        let carpeta = carpeta("resumen");
        let sam = carpeta.join("real.sam");
        let real = mapeo::mapear_fastq(&ruta, indice, genoma, ecoli::IDENTIFICADOR, &sam).unwrap();
        assert_eq!(resumen, real);
        fs::remove_dir_all(&carpeta).unwrap();
    }

    #[test]
    fn la_demostracion_de_minimizadores_repite_al_nucleo() {
        let genoma = ecoli::genoma();
        for inicio in [0, 1_000, 250_000, 4_000_000] {
            let tramo = &genoma[inicio..inicio + 400];
            let (instantes, elegidos) = seguir_minimizadores(tramo, K, W);
            assert_eq!(elegidos, mapeo::minimizadores(tramo, K, W));
            assert_eq!(instantes.len(), tramo.len() - K + 1);
        }
        let mut con_n = ecoli::region(5_001, 5_300).unwrap().to_vec();
        con_n[100] = b'N';
        con_n[180] = b'N';
        let (_, elegidos) = seguir_minimizadores(&con_n, K, W);
        assert_eq!(elegidos, mapeo::minimizadores(&con_n, K, W));
    }

    #[test]
    fn el_paso_a_paso_del_ejemplo_cuenta_la_historia_completa() {
        let id = "puerto.mapear_lecturas";
        let lista = crate::pasos(id, crate::ejemplo(id).unwrap()).unwrap();
        assert!(lista.len() >= 40 && lista.len() <= 400, "{}", lista.len());
        let texto: Vec<&str> = lista.iter().map(|paso| paso.linea.as_str()).collect();
        for esperado in [
            "let mut semillas = minimizadores(referencia, k, w);",
            "let anclas = buscar_anclas(lectura, indice);",
            "let cadenas = elegir_candidatas(encadenar(&anclas, indice.k));",
            "let mejor = candidatos",
            "resumen.mapeadas += 1;",
            "resumen.tramos = tramos_cubiertos(&profundidad);",
            "Ok(resumen)",
        ] {
            assert!(texto.contains(&esperado), "falta {esperado}");
        }
        let ultimo = lista.last().unwrap();
        assert!(ultimo.nota.starts_with("Resultado:"));
        assert!(ultimo.nota.contains("Qué significa"));
        let tamano = serde_json::to_string(&lista).unwrap().len();
        assert!(tamano < 2_000_000, "{tamano}");
    }

    #[test]
    fn el_paso_a_paso_acepta_referencias_propias_y_lecturas_ajenas() {
        let carpeta = carpeta("pasos_propia");
        let referencia = carpeta.join("referencia.fa");
        let tramo = ecoli::region(500_001, 530_000).unwrap();
        fs::write(
            &referencia,
            format!(">trozo\n{}\n", String::from_utf8_lossy(tramo)),
        )
        .unwrap();
        let mut azar = Azar::nuevo(9);
        let mut lecturas = lecturas_cada(&tramo[..3_000], 150);
        let propias = lecturas.len();
        for _ in 0..20 {
            lecturas.push((0..150).map(|_| b"ACGT"[azar.entero(4)]).collect());
        }
        let ruta_lecturas = carpeta.join("lecturas.fq");
        escribir_lecturas(&ruta_lecturas, &lecturas);
        let entrada = Entrada::nueva(&[
            ("lecturas", ruta_lecturas.to_str().unwrap()),
            ("referencia", referencia.to_str().unwrap()),
        ]);
        let lista = pasos_mapeo(&entrada).unwrap();
        assert!(lista.len() <= 400);
        let ultimo = lista.last().unwrap();
        assert!(
            ultimo
                .nota
                .contains(&format!("{} de {}", propias, propias + 20))
        );
        assert!(ultimo.nota.contains("trozo"));
        assert!(lista.iter().any(|paso| paso.linea == "bandera: 4,"));
        fs::remove_dir_all(&carpeta).unwrap();
    }

    #[test]
    fn el_paso_a_paso_da_los_mismos_errores_que_el_mapeo() {
        let carpeta = carpeta("pasos_malas");
        let lecturas = carpeta.join("lecturas.fastq");
        escribir_lecturas(
            &lecturas,
            &lecturas_cada(ecoli::region(1, 2_000).unwrap(), 100),
        );
        let dos = carpeta.join("dos.fasta");
        fs::write(
            &dos,
            ">a\nACGTACGTACGTACGTACGTACGTACGT\n>b\nACGTACGTACGTACGTACGTAC\n",
        )
        .unwrap();
        let corta = carpeta.join("corta.fasta");
        fs::write(&corta, ">c\nACGTACGT\n").unwrap();
        for referencia in [
            "@nada",
            "xyz",
            "",
            "0",
            "-5",
            dos.to_str().unwrap(),
            corta.to_str().unwrap(),
        ] {
            let entrada = Entrada::nueva(&[
                ("lecturas", lecturas.to_str().unwrap()),
                ("referencia", referencia),
            ]);
            let esperado = ejecutar_mapeo(&entrada).err();
            let obtenido = pasos_mapeo(&entrada).err();
            assert_eq!(obtenido, esperado, "{referencia}");
        }
        let vacio = carpeta.join("vacio.fastq");
        fs::write(&vacio, "").unwrap();
        let entrada = Entrada::nueva(&[
            ("lecturas", vacio.to_str().unwrap()),
            ("referencia", "@ecoli"),
        ]);
        assert!(
            pasos_mapeo(&entrada)
                .err()
                .unwrap()
                .contains("ninguna lectura")
        );
        fs::remove_dir_all(&carpeta).unwrap();
    }
}
