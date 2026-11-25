use crate::formato::{decimal, miles, porcentaje};
use crate::pizarra::{Marca, Narracion, Pizarra};
use crate::{Algoritmo, Bloque, Campo, Entrada, Paso, Salida, Tipo};
use nucleo::ecoli;
use nucleo::ensamblador::{
    self, Burbuja, Comparacion, Contig, Ensamblaje, GrafoDeUnitigs, Kmero, Kmeros, Punta,
    ULTIMA_BARRA, Umbral,
};
use nucleo::ensamblaje::medir_n50;
use nucleo::formatos::{abrir_fastq, leer_vcf};
use nucleo::secuencia::complemento_inverso;
use nucleo::simulacion::{OpcionesSimulacion, REGIONES_DE_EJEMPLO};
use std::collections::HashSet;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static CONTADOR: AtomicU64 = AtomicU64::new(0);

const MAXIMO_PASOS: usize = 400;
const BARRAS_DE_CONTIGS: usize = 60;
const FILAS_DE_CONTIGS: usize = 50;
const FILAS_DE_HUECOS: usize = 200;
const CONTIG_LARGO: usize = 1000;
const BYTES_PARA_COMPARAR_SIN_LIMPIAR: u64 = 60_000_000;
const VECES_QUE_CRECE_UN_GZ: u64 = 4;
const KMEROS_DETALLADOS: usize = 3;
const PASOS_DE_CAMINO_DETALLADOS: usize = 5;
const PASOS_DEL_VALLE_DETALLADOS: usize = 12;
const LECTURAS_PARA_BUSCAR_UN_ERROR: usize = 2000;
const VISIBLE: usize = 40;
const COBERTURA_BAJA: usize = 8;
const VECES_PARA_REPETICION: f64 = 1.7;
const CONTEXTO_DE_MUTACION: usize = 20;
const LETRAS_DEL_ADAPTADOR: usize = 15;
const MARGEN_POR_MUTACION: usize = 3;
const NOMBRE_DEL_EJEMPLO: &str = "ejemplo_v1_lecturas.fastq";

const CAMPOS: &[Campo] = &[
    Campo {
        id: "lecturas",
        etiqueta: "Lecturas cortas (FASTQ o FASTQ.gz)",
        tipo: Tipo::Archivo,
        ejemplo: "@ejemplo:lecturas",
    },
    Campo {
        id: "k",
        etiqueta: "k (impar, de 5 a 31)",
        tipo: Tipo::Entero,
        ejemplo: "31",
    },
    Campo {
        id: "cobertura_minima",
        etiqueta: "Cobertura mínima de un k-mero (0 = automática: el primer valle del histograma)",
        tipo: Tipo::Entero,
        ejemplo: "0",
    },
];

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![Algoritmo {
        id: "ensamblaje.ensamblador",
        capitulo: "ensamblaje",
        nombre: "Ensamblador de lecturas reales (De Bruijn tolerante a errores)",
        resumen: "Arma contigs desde un FASTQ con errores: a diferencia de «Contigs», que supone k-meros perfectos, descarta los k-meros débiles con el valle de su histograma y limpia el grafo de De Bruijn quitando puntas y reventando burbujas, como Velvet y SPAdes.",
        campos: CAMPOS,
        funcion: "ensamblador::ensamblar",
        ejecutar: ensamblar,
        pasos: Some(pasos_ensamblar),
    }]
}

struct Pedido {
    ruta: PathBuf,
    k: usize,
    minimo: usize,
    es_ejemplo: bool,
}

fn leer_pedido(entrada: &Entrada) -> Result<Pedido, String> {
    let k = entrada.natural("k")?;
    ensamblador::validar_k(k)?;
    let minimo = entrada.natural("cobertura_minima")?;
    if minimo > ensamblador::COBERTURA_MAXIMA_PEDIDA {
        return Err(format!(
            "La cobertura mínima debe estar entre 0 (automática) y {}.",
            miles(ensamblador::COBERTURA_MAXIMA_PEDIDA)
        ));
    }
    let pide_el_ejemplo = entrada.texto("lecturas")?.trim() == "@ejemplo:lecturas";
    let ruta = entrada.ruta("lecturas")?;
    let es_ejemplo = pide_el_ejemplo || es_el_archivo_de_ejemplo(&ruta);
    Ok(Pedido {
        ruta,
        k,
        minimo,
        es_ejemplo,
    })
}

fn es_el_archivo_de_ejemplo(ruta: &Path) -> bool {
    let carpeta = std::env::temp_dir().join("la-compania");
    ruta.file_name()
        .is_some_and(|nombre| nombre == NOMBRE_DEL_EJEMPLO)
        && ruta.parent().is_some_and(|padre| padre == carpeta)
}

fn ensamblar(entrada: &Entrada) -> Result<Salida, String> {
    let pedido = leer_pedido(entrada)?;
    let (resultado, sin_limpiar) = rayon::join(
        || ensamblador::ensamblar(&pedido.ruta, pedido.k, pedido.minimo),
        || base_sin_limpiar(&pedido.ruta, pedido.k),
    );
    let resultado = resultado?;
    let destino = guardar_contigs(&pedido.ruta, &resultado.contigs)?;
    let contraste = if pedido.es_ejemplo {
        Some(contrastar_con_ecoli(&resultado))
    } else {
        None
    };
    let mut salida = Salida::default();
    agregar_valores(&mut salida, &resultado, sin_limpiar.as_deref());
    agregar_histograma(&mut salida, &resultado);
    agregar_contigs(&mut salida, &resultado, contraste.as_ref());
    if let Some(contraste) = &contraste {
        agregar_contraste(&mut salida, contraste);
    }
    salida.agregar(Bloque::archivo("Contigs (FASTA)", &destino));
    explicar(
        &mut salida,
        &resultado,
        sin_limpiar.as_deref(),
        contraste.as_ref(),
    );
    Ok(salida)
}

fn base_sin_limpiar(ruta: &Path, k: usize) -> Option<Vec<usize>> {
    let mut bytes = fs::metadata(ruta).ok()?.len();
    if ruta.extension().is_some_and(|extension| extension == "gz") {
        bytes = bytes.saturating_mul(VECES_QUE_CRECE_UN_GZ);
    }
    if bytes > BYTES_PARA_COMPARAR_SIN_LIMPIAR {
        return None;
    }
    ensamblador::largos_sin_limpiar(ruta, k).ok()
}

fn guardar_contigs(ruta: &Path, contigs: &[Contig]) -> Result<PathBuf, String> {
    let destino = crate::ruta_de_salida(ruta, "contigs", "fasta");
    let numero = CONTADOR.fetch_add(1, Ordering::Relaxed);
    let nombre = destino
        .file_name()
        .map(|nombre| nombre.to_string_lossy().to_string())
        .unwrap_or_else(|| "contigs.fasta".to_string());
    let temporal = destino.with_file_name(format!(".{}-{numero}-{nombre}", std::process::id()));
    let escrito = File::create(&temporal)
        .map_err(|error| format!("No se pudo crear {}: {error}", temporal.display()))
        .and_then(|archivo| {
            let mut escritor = BufWriter::new(archivo);
            ensamblador::escribir_contigs(&mut escritor, contigs)
                .and_then(|_| escritor.flush())
                .map_err(|error| format!("No se pudo escribir {}: {error}", temporal.display()))
        })
        .and_then(|_| {
            fs::rename(&temporal, &destino)
                .map_err(|error| format!("No se pudo guardar {}: {error}", destino.display()))
        });
    if let Err(error) = escrito {
        let _ = fs::remove_file(&temporal);
        return Err(error);
    }
    Ok(destino)
}

struct Contraste {
    comparacion: Comparacion,
    mutaciones: Vec<(usize, usize)>,
    mutaciones_en_contigs: usize,
    huecos_en_mutaciones: usize,
    puntas_con_adaptador: usize,
    kmeros_de_burbujas: usize,
    kmeros_de_burbujas_ajenos: usize,
}

fn contrastar_con_ecoli(resultado: &Ensamblaje) -> Contraste {
    let genoma = ecoli::genoma();
    let indice = ensamblador::indexar_referencia(genoma, &REGIONES_DE_EJEMPLO, resultado.k);
    let comparacion =
        ensamblador::comparar_con_referencia(&resultado.contigs, &indice, CONTIG_LARGO);
    let mut kmeros_de_burbujas = 0;
    let mut kmeros_de_burbujas_ajenos = 0;
    for burbuja in &resultado.burbujas {
        kmeros_de_burbujas += burbuja.quitada.len() + 1 - resultado.k;
        kmeros_de_burbujas_ajenos += indice.ajenos(&burbuja.quitada);
    }
    let mutaciones = mutaciones_del_ejemplo();
    let huecos_en_mutaciones = comparacion
        .huecos
        .iter()
        .filter(|hueco| {
            mutaciones
                .iter()
                .any(|mutacion| se_tocan(**hueco, (mutacion.0, mutacion.1)))
        })
        .count();
    let mutaciones_en_contigs = mutaciones
        .iter()
        .filter(|mutacion| mutacion_en_contigs(mutacion, &resultado.contigs, &comparacion))
        .count();
    let puntas_con_adaptador = resultado
        .puntas
        .iter()
        .filter(|punta| lleva_el_adaptador(&punta.secuencia))
        .count();
    Contraste {
        mutaciones: mutaciones
            .iter()
            .map(|mutacion| (mutacion.0, mutacion.1))
            .collect(),
        comparacion,
        mutaciones_en_contigs,
        huecos_en_mutaciones,
        puntas_con_adaptador,
        kmeros_de_burbujas,
        kmeros_de_burbujas_ajenos,
    }
}

fn se_tocan(uno: (usize, usize), otro: (usize, usize)) -> bool {
    uno.0 < otro.1 && otro.0 < uno.1
}

fn mutaciones_del_ejemplo() -> Vec<(usize, usize, Vec<u8>)> {
    let Ok(ruta) = crate::ejemplos::generar("mutaciones") else {
        return Vec::new();
    };
    let Ok(texto) = fs::read(&ruta) else {
        return Vec::new();
    };
    let mut mutaciones = Vec::new();
    for linea in texto.split(|&byte| byte == b'\n') {
        if let Ok(Some(variante)) = leer_vcf(linea)
            && variante.posicion > 0
        {
            let inicio = variante.posicion - 1;
            mutaciones.push((
                inicio,
                inicio + variante.referencia.len().max(1),
                variante.alternativa,
            ));
        }
    }
    mutaciones
}

fn mutacion_en_contigs(
    mutacion: &(usize, usize, Vec<u8>),
    contigs: &[Contig],
    comparacion: &Comparacion,
) -> bool {
    let genoma = ecoli::genoma();
    let (inicio, fin, alternativa) = mutacion;
    if *inicio < CONTEXTO_DE_MUTACION || fin + CONTEXTO_DE_MUTACION > genoma.len() {
        return false;
    }
    let mut contexto = genoma[inicio - CONTEXTO_DE_MUTACION..*inicio].to_vec();
    contexto.extend_from_slice(alternativa);
    contexto.extend_from_slice(&genoma[*fin..fin + CONTEXTO_DE_MUTACION]);
    let inverso = complemento_inverso(&contexto);
    let lo_tiene = |contig: &Contig| {
        contiene(&contig.secuencia, &contexto) || contiene(&contig.secuencia, &inverso)
    };
    let candidatos: Vec<&Contig> = contigs
        .iter()
        .zip(&comparacion.ubicaciones)
        .filter(|(_, ubicacion)| {
            ubicacion.is_some_and(|(desde, hasta)| desde <= *inicio && *fin <= hasta)
        })
        .map(|(contig, _)| contig)
        .collect();
    if candidatos.iter().any(|contig| lo_tiene(contig)) {
        return true;
    }
    contigs.iter().any(lo_tiene)
}

fn contiene(texto: &[u8], pieza: &[u8]) -> bool {
    pieza.len() <= texto.len() && texto.windows(pieza.len()).any(|ventana| ventana == pieza)
}

fn lleva_el_adaptador(secuencia: &[u8]) -> bool {
    let adaptador = OpcionesSimulacion::default().adaptador;
    let pieza = &adaptador[..LETRAS_DEL_ADAPTADOR.min(adaptador.len())];
    contiene(secuencia, pieza) || contiene(secuencia, &complemento_inverso(pieza))
}

fn pico(histograma: &[usize], desde: usize) -> Option<usize> {
    let mut mejor: Option<usize> = None;
    for cobertura in desde.max(1)..histograma.len().min(ULTIMA_BARRA) {
        if mejor.is_none_or(|actual| histograma[cobertura] > histograma[actual]) {
            mejor = Some(cobertura);
        }
    }
    mejor.filter(|&cobertura| {
        histograma[cobertura] > 0
            && histograma
                .get(ULTIMA_BARRA)
                .is_none_or(|&arriba| arriba <= histograma[cobertura])
    })
}

fn largo_medio(resultado: &Ensamblaje) -> f64 {
    resultado.bases as f64 / resultado.lecturas.max(1) as f64
}

fn cobertura_por_base(resultado: &Ensamblaje, pico_kmeros: usize) -> f64 {
    let largo = largo_medio(resultado);
    let kmeros_por_lectura = largo + 1.0 - resultado.k as f64;
    if kmeros_por_lectura <= 0.0 {
        return 0.0;
    }
    pico_kmeros as f64 * largo / kmeros_por_lectura
}

fn describir_umbral(umbral: &Umbral) -> String {
    match (umbral.automatico, umbral.valle) {
        (true, Some(_)) => format!(
            "{} (automática: primer valle del histograma)",
            umbral.minimo
        ),
        (true, None) => format!(
            "{} (automática: el histograma no tiene valle)",
            umbral.minimo
        ),
        (false, Some(valle)) => {
            format!("{} (la elegiste; el valle está en {valle})", umbral.minimo)
        }
        (false, None) => format!(
            "{} (la elegiste; el histograma no tiene valle)",
            umbral.minimo
        ),
    }
}

fn agregar_valores(salida: &mut Salida, resultado: &Ensamblaje, sin_limpiar: Option<&[usize]>) {
    let k = resultado.k;
    salida.agregar(Bloque::valor(
        &format!("Contigs (de 2k = {} letras o más)", 2 * k),
        miles(resultado.contigs.len()),
    ));
    salida.agregar(Bloque::valor(
        "N50",
        format!("{} letras", miles(resultado.medida.n50)),
    ));
    salida.agregar(Bloque::valor(
        "Contig más largo",
        format!(
            "{} letras",
            miles(
                resultado
                    .contigs
                    .first()
                    .map_or(0, |contig| contig.secuencia.len())
            )
        ),
    ));
    salida.agregar(Bloque::valor(
        "Largo total",
        format!("{} letras", miles(resultado.medida.total)),
    ));
    salida.agregar(Bloque::valor(
        "Lecturas",
        format!(
            "{} ({} bases)",
            miles(resultado.lecturas),
            miles(resultado.bases)
        ),
    ));
    salida.agregar(Bloque::valor(
        "k-meros",
        format!(
            "{} leídos, {} distintos, {} sólidos",
            miles(resultado.kmeros_leidos),
            miles(resultado.kmeros_distintos),
            miles(resultado.kmeros_solidos)
        ),
    ));
    salida.agregar(Bloque::valor(
        "Cobertura mínima usada",
        describir_umbral(&resultado.umbral),
    ));
    if let Some(cima) = pico(&resultado.histograma, resultado.umbral.minimo as usize) {
        salida.agregar(Bloque::valor(
            "Pico de cobertura de k-meros",
            format!(
                "{} (≈ {}× por base)",
                cima,
                decimal(cobertura_por_base(resultado, cima), 1)
            ),
        ));
    }
    salida.agregar(Bloque::valor(
        "Limpieza del grafo",
        format!(
            "{}; {} y {} en {}",
            cantidad(resultado.unitigs_iniciales, "unitig", "unitigs"),
            cantidad(resultado.puntas.len(), "punta quitada", "puntas quitadas"),
            cantidad(
                resultado.burbujas.len(),
                "burbuja reventada",
                "burbujas reventadas"
            ),
            cantidad(resultado.rondas, "ronda", "rondas")
        ),
    ));
    if let Some(largos) = sin_limpiar {
        let medida = medir_n50(largos);
        salida.agregar(Bloque::valor(
            "Sin filtro ni limpieza (como «Contigs»)",
            format!(
                "{} unitigs, N50 = {} letras",
                miles(largos.len()),
                miles(medida.n50)
            ),
        ));
    }
}

fn agregar_histograma(salida: &mut Salida, resultado: &Ensamblaje) {
    let minimo = resultado.umbral.minimo as usize;
    let cima = pico(&resultado.histograma, minimo).unwrap_or(minimo);
    let limite = (3 * cima).max(minimo + 5).clamp(10, 120).min(ULTIMA_BARRA);
    let categorias = (1..=limite)
        .map(|cobertura| {
            if cobertura == minimo {
                format!("{cobertura} (mínimo)")
            } else {
                cobertura.to_string()
            }
        })
        .collect();
    let valores = (1..=limite)
        .map(|cobertura| resultado.histograma[cobertura] as f64)
        .collect();
    salida.agregar(Bloque::barras(
        &format!(
            "Histograma de cobertura de k-meros: se descartan los que aparecen menos de {minimo} veces"
        ),
        "Veces que aparece un k-mero",
        "k-meros distintos",
        categorias,
        valores,
    ));
}

fn agregar_contigs(salida: &mut Salida, resultado: &Ensamblaje, contraste: Option<&Contraste>) {
    if resultado.contigs.is_empty() {
        return;
    }
    let mostrados = resultado.contigs.len().min(BARRAS_DE_CONTIGS);
    salida.agregar(Bloque::barras(
        &format!(
            "Largo de los {} contigs más largos (N50 = {})",
            mostrados,
            miles(resultado.medida.n50)
        ),
        "Contig (del más largo al más corto)",
        "Letras",
        (1..=mostrados).map(|numero| numero.to_string()).collect(),
        resultado.contigs[..mostrados]
            .iter()
            .map(|contig| contig.secuencia.len() as f64)
            .collect(),
    ));
    let mut columnas = vec!["Contig", "Largo", "Cobertura de k-meros"];
    if contraste.is_some() {
        columnas.push("Dónde cae en E. coli");
    }
    columnas.push("Comienzo");
    let filas = resultado
        .contigs
        .iter()
        .take(FILAS_DE_CONTIGS)
        .enumerate()
        .map(|(numero, contig)| {
            let mut fila = vec![
                (numero + 1).to_string(),
                miles(contig.secuencia.len()),
                decimal(contig.cobertura, 1),
            ];
            if let Some(contraste) = contraste {
                fila.push(describir_ubicacion(
                    contraste
                        .comparacion
                        .ubicaciones
                        .get(numero)
                        .copied()
                        .flatten(),
                ));
            }
            fila.push(recortar(&contig.secuencia, VISIBLE));
            fila
        })
        .collect();
    salida.agregar(Bloque::tabla_con_total(
        "Contigs más largos",
        &columnas,
        filas,
        resultado.contigs.len(),
    ));
}

fn describir_ubicacion(ubicacion: Option<(usize, usize)>) -> String {
    match ubicacion {
        Some((inicio, fin)) => format!("{}–{}", miles(inicio + 1), miles(fin)),
        None => "varias copias o fuera de las ventanas".to_string(),
    }
}

fn agregar_contraste(salida: &mut Salida, contraste: &Contraste) {
    let comparacion = &contraste.comparacion;
    salida.agregar(Bloque::valor(
        "Ventanas verdaderas de E. coli cubiertas",
        format!(
            "{} de {} letras ({})",
            miles(comparacion.cubiertas),
            miles(comparacion.largo_referencia),
            porcentaje(
                fraccion(comparacion.cubiertas, comparacion.largo_referencia),
                1
            )
        ),
    ));
    salida.agregar(Bloque::valor(
        &format!(
            "Cubierto por contigs de {} letras o más",
            miles(CONTIG_LARGO)
        ),
        porcentaje(
            fraccion(
                comparacion.cubiertas_por_largos,
                comparacion.largo_referencia,
            ),
            1,
        ),
    ));
    salida.agregar(Bloque::valor(
        "Mutaciones conocidas dentro de un contig",
        format!(
            "{} de {}",
            contraste.mutaciones_en_contigs,
            contraste.mutaciones.len()
        ),
    ));
    salida.agregar(Bloque::valor(
        "k-meros de los contigs que no están en E. coli",
        format!(
            "{} de {}",
            miles(comparacion.kmeros_ajenos),
            miles(comparacion.kmeros_de_contigs)
        ),
    ));
    let filas = comparacion
        .huecos
        .iter()
        .take(FILAS_DE_HUECOS)
        .map(|&(inicio, fin)| {
            let motivo = if contraste
                .mutaciones
                .iter()
                .any(|&mutacion| se_tocan((inicio, fin), mutacion))
            {
                "mutación de la muestra"
            } else {
                "repetición o poca cobertura"
            };
            vec![miles(inicio + 1), miles(fin - inicio), motivo.to_string()]
        })
        .collect();
    salida.agregar(Bloque::tabla_con_total(
        "Huecos frente al genoma de E. coli",
        &["Desde", "Letras", "Motivo probable"],
        filas,
        comparacion.huecos.len(),
    ));
}

fn fraccion(parte: usize, total: usize) -> f64 {
    if total == 0 {
        0.0
    } else {
        parte as f64 / total as f64
    }
}

fn recortar(secuencia: &[u8], maximo: usize) -> String {
    let texto = String::from_utf8_lossy(&secuencia[..secuencia.len().min(maximo)]).to_string();
    if secuencia.len() > maximo {
        format!("{texto}…")
    } else {
        texto
    }
}

fn mediana_de_cobertura(contigs: &[Contig]) -> f64 {
    let mut pares: Vec<(f64, usize)> = contigs
        .iter()
        .map(|contig| (contig.cobertura, contig.secuencia.len()))
        .collect();
    pares.sort_by(|uno, otro| uno.0.total_cmp(&otro.0));
    let total: usize = pares.iter().map(|par| par.1).sum();
    let mut acumulado = 0;
    for (cobertura, largo) in pares {
        acumulado += largo;
        if 2 * acumulado >= total {
            return cobertura;
        }
    }
    0.0
}

fn cantidad(numero: usize, singular: &str, plural: &str) -> String {
    if numero == 1 {
        format!("1 {singular}")
    } else {
        format!("{} {plural}", miles(numero))
    }
}

fn explicar(
    salida: &mut Salida,
    resultado: &Ensamblaje,
    sin_limpiar: Option<&[usize]>,
    contraste: Option<&Contraste>,
) {
    let k = resultado.k;
    let umbral = resultado.umbral;
    if resultado.contigs.is_empty() {
        salida.explicar(format!(
            "No salió ningún contig de 2k = {} letras: con cobertura mínima {} casi no quedaron k-meros sólidos; prueba con 0 para que el mínimo salga del histograma.",
            2 * k,
            umbral.minimo
        ));
    } else {
        salida.explicar(format!(
            "Salen {} de al menos 2k = {} letras; el mayor mide {} y el N50 es {}: la mitad de las {} letras ensambladas está en contigs de ese largo o más.",
            cantidad(resultado.contigs.len(), "contig", "contigs"),
            2 * k,
            miles(resultado.contigs[0].secuencia.len()),
            miles(resultado.medida.n50),
            miles(resultado.medida.total)
        ));
    }
    let descartados = resultado.kmeros_distintos - resultado.kmeros_solidos;
    let cima = pico(&resultado.histograma, umbral.minimo as usize);
    let aviso_de_cobertura = match cima {
        Some(cima) if cima < COBERTURA_BAJA => {
            " Es una cobertura baja: espera huecos y contigs cortos, porque nadie recupera lo que no se leyó."
        }
        _ => "",
    };
    match (umbral.automatico, umbral.valle, cima) {
        (true, Some(valle), Some(cima)) => salida.explicar(format!(
            "El histograma baja desde los k-meros vistos una vez (casi todos errores de lectura) hasta un valle en {valle} y sube a un pico en {cima}, la cobertura real (≈ {}× por base): se tiraron {} de {} k-meros distintos ({}).{aviso_de_cobertura}",
            decimal(cobertura_por_base(resultado, cima), 1),
            miles(descartados),
            miles(resultado.kmeros_distintos),
            porcentaje(fraccion(descartados, resultado.kmeros_distintos), 1)
        )),
        (true, Some(valle), None) => salida.explicar(format!(
            "El histograma baja hasta un valle en {valle}, pero el pico de los k-meros reales pasa de {ULTIMA_BARRA} copias (la última barra junta todo lo que está por encima): es una cobertura altísima, típica de un genoma muy pequeño como un plásmido o un fago, y el valle sigue siendo válido."
        )),
        (true, _, _) => salida.explicar(format!(
            "El histograma no tiene valle: con tan poca cobertura, o tantos errores, errores y genoma se mezclan; se usó {} por omisión y conviene desconfiar de los contigs.",
            umbral.minimo
        )),
        (false, valle, _) => salida.explicar(format!(
            "Elegiste cobertura mínima {}{}: muy por encima del valle pierdes k-meros reales y se abren huecos, y por debajo entran errores que fragmentan el grafo.{aviso_de_cobertura}",
            umbral.minimo,
            valle.map_or(String::new(), |valle| format!(
                " y el valle automático está en {valle}"
            ))
        )),
    }
    let comparado = match sin_limpiar.filter(|largos| !largos.is_empty()) {
        Some(largos) => format!(
            "; sin filtro ni limpieza, como supone «Contigs» del capítulo 3 al tratar cada k-mero como real, saldrían {} unitigs con N50 de {} letras",
            miles(largos.len()),
            miles(medir_n50(largos).n50)
        ),
        None => String::new(),
    };
    salida.explicar(format!(
        "La limpieza quitó {} (callejones de menos de 2k k-meros con menos cobertura que su rama hermana, como en Velvet) y reventó {} (ramas con otro camino de largo parecido y más cubierto entre los mismos puntos, como en SPAdes){comparado}.",
        cantidad(resultado.puntas.len(), "punta", "puntas"),
        cantidad(resultado.burbujas.len(), "burbuja", "burbujas")
    ));
    let mediana = mediana_de_cobertura(&resultado.contigs);
    let repetidos = resultado
        .contigs
        .iter()
        .filter(|contig| contig.cobertura >= VECES_PARA_REPETICION * mediana)
        .count();
    if repetidos > 0 && mediana > 0.0 {
        salida.explicar(format!(
            "{} {} o más veces la cobertura típica ({}): son repeticiones que el grafo juntó en una sola copia, y en sus bordes se cortan los contigs; confía más en los contigs largos de cobertura típica.",
            cantidad(repetidos, "contig tiene", "contigs tienen"),
            decimal(VECES_PARA_REPETICION, 1),
            decimal(mediana, 1)
        ));
    }
    match contraste {
        Some(contraste) => explicar_contraste(salida, resultado, contraste),
        None => salida.explicar(
            "Con lecturas cortas sin pares, ningún grafo de De Bruijn cruza repeticiones más largas que las lecturas: para unir estos contigs hacen falta lecturas pareadas o largas.",
        ),
    }
}

fn explicar_contraste(salida: &mut Salida, resultado: &Ensamblaje, contraste: &Contraste) {
    let comparacion = &contraste.comparacion;
    let total_mutaciones = contraste.mutaciones.len();
    let sobre_mutaciones = if total_mutaciones == 0 {
        String::new()
    } else if contraste.mutaciones_en_contigs == total_mutaciones {
        format!(
            " y las {total_mutaciones} mutaciones de la muestra aparecen dentro de un contig con su letra cambiada: no lo cortan, porque todas las lecturas las llevan, aunque {} de los {} huecos frente a la referencia caen sobre ellas",
            contraste.huecos_en_mutaciones,
            miles(comparacion.huecos.len())
        )
    } else {
        format!(
            " y solo {} de las {total_mutaciones} mutaciones de la muestra aparecen dentro de un contig con su letra cambiada: las otras quedaron en huecos o en contigs descartados, y {} de los {} huecos frente a la referencia caen sobre ellas",
            contraste.mutaciones_en_contigs,
            contraste.huecos_en_mutaciones,
            miles(comparacion.huecos.len())
        )
    };
    salida.explicar(format!(
        "Contra el genoma de E. coli, los contigs cubren {} de las dos ventanas verdaderas ({} en contigs de {} letras o más){sobre_mutaciones}; los demás huecos y los cortes suelen estar en repeticiones, como los operones de ARN ribosomal.",
        porcentaje(fraccion(comparacion.cubiertas, comparacion.largo_referencia), 1),
        porcentaje(
            fraccion(comparacion.cubiertas_por_largos, comparacion.largo_referencia),
            1
        ),
        miles(CONTIG_LARGO)
    ));
    let adaptador = if contraste.puntas_con_adaptador > 0 {
        format!(
            "{} de las {} quitadas llevan el adaptador de la librería (no son genoma)",
            contraste.puntas_con_adaptador,
            cantidad(resultado.puntas.len(), "punta", "puntas")
        )
    } else {
        String::new()
    };
    let ramas = if resultado.burbujas.is_empty() {
        String::new()
    } else if contraste.kmeros_de_burbujas_ajenos == 0 {
        "las dos ramas de cada burbuja existen en el genoma: son copias de una repetición que difieren en pocas letras, así que el contig resume todas las copias".to_string()
    } else {
        format!(
            "{} de los {} k-meros de ramas quitadas en burbujas no existen en el genoma: eran errores",
            miles(contraste.kmeros_de_burbujas_ajenos),
            miles(contraste.kmeros_de_burbujas)
        )
    };
    let explicables = contraste.mutaciones.len() * (resultado.k + MARGEN_POR_MUTACION);
    let ajenos = if comparacion.kmeros_ajenos <= explicables {
        format!(
            "solo {} de los {} k-meros de los contigs no existen en la referencia, no más de los que cruzan las mutaciones: los contigs no copiaron errores de lectura",
            miles(comparacion.kmeros_ajenos),
            miles(comparacion.kmeros_de_contigs)
        )
    } else {
        format!(
            "{} de los {} k-meros de los contigs no existen en la referencia, más de los que explican las mutaciones: desconfía de los contigs cortos",
            miles(comparacion.kmeros_ajenos),
            miles(comparacion.kmeros_de_contigs)
        )
    };
    let partes: Vec<String> = [adaptador, ramas, ajenos]
        .into_iter()
        .filter(|parte| !parte.is_empty())
        .collect();
    if !partes.is_empty() {
        let frase = partes.join("; ");
        let mut letras = frase.chars();
        let primera = letras
            .next()
            .map(|letra| letra.to_uppercase().to_string())
            .unwrap_or_default();
        salida.explicar(format!("{primera}{}.", letras.as_str()));
    }
}

fn pasos_ensamblar(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let pedido = leer_pedido(entrada)?;
    narrar(&pedido).map(|(pasos, _)| pasos)
}

#[derive(Debug, Clone, PartialEq)]
struct Final {
    contigs: usize,
    n50: usize,
    total: usize,
    mayor: usize,
    puntas: usize,
    burbujas: usize,
    rondas: usize,
}

fn letras(codigo: u64, k: usize) -> String {
    String::from_utf8_lossy(&ensamblador::letras_de(codigo, k)).to_string()
}

fn nombre_de_archivo(ruta: &Path) -> String {
    ruta.file_name()
        .map(|nombre| nombre.to_string_lossy().to_string())
        .unwrap_or_else(|| ruta.display().to_string())
}

fn narrar(pedido: &Pedido) -> Result<(Vec<Paso>, Final), String> {
    let k = pedido.k;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        Pizarra::nueva()
            .variable("ruta", nombre_de_archivo(&pedido.ruta))
            .variable("k", k)
            .variable("cobertura_minima", pedido.minimo)
            .paso(
                "let mut conteo = contar_kmeros(ruta, k)?;",
                format!(
                    "Vas a ensamblar «{}» con k = {k}. El plan: contar los k-meros de todas las lecturas, tirar los que aparecen muy pocas veces (errores), unir los que quedan en unitigs y limpiar el grafo de puntas y burbujas.",
                    nombre_de_archivo(&pedido.ruta)
                ),
            ),
    );
    narrar_primera_lectura(&mut narracion, &pedido.ruta, k)?;
    let mut conteo = ensamblador::contar_kmeros(&pedido.ruta, k)?;
    let kmeros_distintos = conteo.kmeros.cantidad();
    narracion.agregar(
        Pizarra::nueva()
            .variable("conteo.lecturas", miles(conteo.lecturas))
            .variable("conteo.bases", miles(conteo.bases))
            .variable("conteo.kmeros_leidos", miles(conteo.kmeros_leidos))
            .variable("distintos", miles(kmeros_distintos))
            .paso(
                "conteo.kmeros_leidos += cubeta.len();",
                format!(
                    "Así con todas: leyó {} lecturas ({} bases) por tandas y repartió {} k-meros en 64 cubetas; {} núcleos suman las cubetas a la vez. Hay {} k-meros distintos.",
                    miles(conteo.lecturas),
                    miles(conteo.bases),
                    miles(conteo.kmeros_leidos),
                    crate::nucleos(),
                    miles(kmeros_distintos)
                ),
            ),
    );
    let histograma = ensamblador::histograma(&conteo.kmeros);
    narrar_histograma(&mut narracion, &histograma);
    let umbral = ensamblador::elegir_umbral(&histograma, pedido.minimo);
    narrar_umbral(&mut narracion, &umbral, &histograma);
    narrar_lectura_con_error(&mut narracion, &pedido.ruta, &conteo.kmeros, umbral.minimo)?;
    conteo.kmeros.retener(umbral.minimo);
    let solidos = conteo.kmeros.cantidad();
    narracion.agregar(
        Pizarra::nueva()
            .variable("umbral.minimo", umbral.minimo)
            .variable("distintos", miles(kmeros_distintos))
            .variable("solidos", miles(solidos))
            .paso(
                "conteo.kmeros.retener(umbral.minimo);",
                format!(
                    "Tira de las 64 tablas, en paralelo, los k-meros con menos de {} copias: quedan {} k-meros sólidos de {} ({} se fueron).",
                    umbral.minimo,
                    miles(solidos),
                    miles(kmeros_distintos),
                    porcentaje(fraccion(kmeros_distintos - solidos, kmeros_distintos), 1)
                ),
            ),
    );
    let mut kmeros = conteo.kmeros;
    narrar_primer_unitig(&mut narracion, &kmeros);
    let mut grafo = ensamblador::compactar(&kmeros);
    narracion.agregar(
        Pizarra::nueva()
            .variable("grafo.unitigs.len()", miles(grafo.unitigs.len()))
            .paso(
                "let mut grafo = compactar(kmeros);",
                format!(
                    "Hace lo mismo desde cada k-mero sin visitar: salen {} unitigs{}. Los ordena del más largo al más corto y conecta los extremos que se solapan en k − 1 = {} letras.",
                    miles(grafo.unitigs.len()),
                    grafo
                        .unitigs
                        .first()
                        .map_or(String::new(), |unitig| format!(
                            " (el más largo mide {} letras)",
                            miles(unitig.secuencia.len())
                        )),
                    k - 1
                ),
            ),
    );
    let mut rondas = 0;
    let mut puntas = 0;
    let mut burbujas = 0;
    let mut punta_contada = false;
    let mut burbuja_contada = false;
    while rondas < ensamblador::MAXIMAS_RONDAS {
        let puntas_de_la_ronda = ensamblador::buscar_puntas(&grafo);
        let burbujas_de_la_ronda = ensamblador::buscar_burbujas(&grafo);
        if puntas_de_la_ronda.is_empty() && burbujas_de_la_ronda.is_empty() {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("rondas", rondas)
                    .paso(
                        "if puntas_de_la_ronda.is_empty() && burbujas_de_la_ronda.is_empty() {",
                        if rondas == 0 {
                            "No hay ninguna punta ni burbuja que quitar: el grafo ya está limpio.".to_string()
                        } else {
                            format!("En la ronda {} ya no aparece ninguna punta ni burbuja: el grafo quedó limpio.", rondas + 1)
                        },
                    ),
            );
            break;
        }
        rondas += 1;
        if !punta_contada && let Some(punta) = puntas_de_la_ronda.first() {
            narrar_punta(&mut narracion, &grafo, punta, pedido.es_ejemplo);
            punta_contada = true;
        }
        if !burbuja_contada && let Some(burbuja) = burbujas_de_la_ronda.first() {
            narrar_burbuja(&mut narracion, &grafo, burbuja, pedido.es_ejemplo);
            burbuja_contada = true;
        }
        let mut quitados = Vec::new();
        for punta in &puntas_de_la_ronda {
            ensamblador::quitar_unitig(&mut kmeros, &grafo.unitigs[punta.unitig]);
            quitados.push(punta.unitig);
        }
        for burbuja in &burbujas_de_la_ronda {
            ensamblador::quitar_unitig(&mut kmeros, &grafo.unitigs[burbuja.unitig_quitado]);
            quitados.push(burbuja.unitig_quitado);
        }
        grafo = ensamblador::recompactar(grafo, &kmeros, &quitados);
        puntas += puntas_de_la_ronda.len();
        burbujas += burbujas_de_la_ronda.len();
        narracion.agregar(
            Pizarra::nueva()
                .variable("rondas", rondas)
                .variable("puntas_de_la_ronda.len()", miles(puntas_de_la_ronda.len()))
                .variable("burbujas_de_la_ronda.len()", miles(burbujas_de_la_ronda.len()))
                .variable("grafo.unitigs.len()", miles(grafo.unitigs.len()))
                .paso(
                    "grafo = recompactar(grafo, kmeros, &quitados);",
                    format!(
                        "Ronda {rondas}: quitó {} y reventó {}. Vuelve a caminar solo los unitigs que tocaban a los quitados: los tramos que quedaron sin rival se unen y ahora hay {} unitigs.",
                        cantidad(puntas_de_la_ronda.len(), "punta", "puntas"),
                        cantidad(burbujas_de_la_ronda.len(), "burbuja", "burbujas"),
                        miles(grafo.unitigs.len())
                    ),
                ),
        );
    }
    let unitigs_finales = grafo.unitigs.len();
    let (contigs, cortos) = ensamblador::elegir_contigs(grafo);
    let largos: Vec<usize> = contigs
        .iter()
        .map(|contig| contig.secuencia.len())
        .collect();
    let medida = medir_n50(&largos);
    let mayor = largos.first().copied().unwrap_or(0);
    narracion.agregar(
        Pizarra::nueva()
            .variable("largo_minimo", 2 * k)
            .variable("contigs.len()", miles(contigs.len()))
            .variable("cortos", miles(cortos))
            .lista(
                "Los contigs más largos",
                largos
                    .iter()
                    .take(10)
                    .enumerate()
                    .map(|(numero, largo)| format!("contig {}: {} letras", numero + 1, miles(*largo)))
                    .collect(),
                &[],
            )
            .paso(
                "if unitig.secuencia.len() >= largo_minimo {",
                format!(
                    "De los {} unitigs finales se queda con los de al menos 2k = {} letras: {} contigs. Los {} más cortos son restos de repeticiones o de bordes y no se guardan.",
                    miles(unitigs_finales),
                    2 * k,
                    miles(contigs.len()),
                    miles(cortos)
                ),
            ),
    );
    let ultimo = Pizarra::nueva()
        .variable("medida.n50", miles(medida.n50))
        .variable("medida.total", miles(medida.total))
        .variable("mayor", miles(mayor))
        .paso(
            "let medida = medir_n50(&largos);",
            format!(
                "Resultado: {} contigs que suman {} letras; el mayor mide {} y el N50 es {}: la mitad del ensamblaje está en contigs de ese largo o más. Mira «Qué significa».",
                miles(contigs.len()),
                miles(medida.total),
                miles(mayor),
                miles(medida.n50)
            ),
        );
    let resumen = Final {
        contigs: contigs.len(),
        n50: medida.n50,
        total: medida.total,
        mayor,
        puntas,
        burbujas,
        rondas,
    };
    Ok((narracion.terminar(ultimo), resumen))
}

fn narrar_primera_lectura(narracion: &mut Narracion, ruta: &Path, k: usize) -> Result<(), String> {
    let Some(lectura) = abrir_fastq(ruta)?.next() else {
        return Ok(());
    };
    let lectura = lectura?;
    let nombre = String::from_utf8_lossy(&lectura.nombre).to_string();
    let mut encontrados: Vec<(usize, Kmero)> = Vec::new();
    ensamblador::recorrer_kmeros(&lectura.secuencia, k, |posicion, kmero| {
        encontrados.push((posicion, kmero))
    });
    narracion.agregar(
        Pizarra::nueva()
            .variable("lectura.nombre", &nombre)
            .variable("lectura.secuencia.len()", lectura.secuencia.len())
            .variable("k", k)
            .secuencia(&lectura.secuencia, &[])
            .paso(
                "recorrer_kmeros(&lectura.secuencia, k, |_, kmero| {",
                format!(
                    "Empieza por la primera lectura, «{nombre}», de {} letras. La recorre letra por letra guardando cada letra en 2 bits (A = 00, C = 01, G = 10, T = 11).",
                    lectura.secuencia.len()
                ),
            ),
    );
    if let Some(corte) = lectura
        .secuencia
        .iter()
        .position(|&letra| ensamblador::codigo_de_base(letra).is_none())
    {
        narracion.agregar(
            Pizarra::nueva()
                .variable("posicion", corte + 1)
                .secuencia(&lectura.secuencia, &[(corte, corte + 1, Marca::Comparado)])
                .paso(
                    "seguidas = 0;",
                    format!(
                        "La letra {} es «{}», que no es A, C, G ni T: corta la cuenta y ningún k-mero la incluye.",
                        corte + 1,
                        lectura.secuencia[corte] as char
                    ),
                ),
        );
    }
    let mut anterior: Option<usize> = None;
    for &(posicion, kmero) in encontrados.iter().take(KMEROS_DETALLADOS) {
        let marcas = [(posicion, posicion + k, Marca::Actual)];
        let tablero = || {
            Pizarra::nueva()
                .variable("posicion", posicion + 1)
                .variable("kmero.directo", letras(kmero.directo, k))
                .variable("kmero.reverso", letras(kmero.reverso, k))
                .secuencia(&lectura.secuencia, &marcas)
        };
        match anterior.filter(|&antes| antes + 1 == posicion) {
            None if anterior.is_some() => narracion.agregar(tablero().paso(
                "if seguidas >= k {",
                format!(
                    "Tras la letra que no es A, C, G ni T, vuelve a juntar k = {k} letras seguidas: este k-mero va de la letra {} a la {}, y se calcula desde cero."
                    , posicion + 1, posicion + k
                ),
            )),
            None => narracion.agregar(tablero().paso(
                "if seguidas >= k {",
                format!(
                    "Ya juntó k = {k} letras seguidas: el primer k-mero va de la letra {} a la {}. Lleva dos versiones a la vez: el directo, como se leyó, y el reverso, su complemento inverso (lo que diría la otra hebra).",
                    posicion + 1,
                    posicion + k
                ),
            )),
            Some(_) => {
                let nueva = lectura.secuencia[posicion + k - 1] as char;
                narracion.agregar(tablero().paso(
                    "directo: ((kmero.directo << 2) | base) & mascara,",
                    format!(
                        "Llega la letra {} («{nueva}»): el directo corre dos bits a la izquierda, pierde su primera letra y suma «{nueva}» al final, sin releer las otras {}.",
                        posicion + k,
                        k - 1
                    ),
                ));
                narracion.agregar(tablero().paso(
                    "reverso: (kmero.reverso >> 2) | ((3 - base) << (2 * (k - 1))),",
                    format!(
                        "El reverso hace lo contrario: corre a la derecha y pone el complemento de «{nueva}» al principio. Así las dos hebras se actualizan con una cuenta cada una."
                    ),
                ));
            }
        }
        let canonico = kmero.canonico();
        let cual = if canonico == kmero.directo {
            "el directo"
        } else {
            "el reverso"
        };
        narracion.agregar(
            tablero()
                .variable("canonico", letras(canonico, k))
                .paso(
                    "self.directo.min(self.reverso)",
                    format!(
                        "El canónico es el menor de los dos en orden alfabético, aquí {cual}. Así una lectura y otra de la hebra opuesta cuentan el mismo k-mero."
                    ),
                ),
        );
        if anterior.is_none() {
            narracion.agregar(
                tablero()
                    .variable("canonico", letras(canonico, k))
                    .variable("parte_de(canonico)", ensamblador::parte_de(canonico))
                    .paso(
                        "cubetas[parte_de(canonico)].push(canonico);",
                        format!(
                            "Lo echa en la cubeta {} de 64, elegida por una mezcla de sus bits. Cada cubeta tendrá su propia tabla, así varios núcleos cuentan a la vez sin pisarse.",
                            ensamblador::parte_de(canonico)
                        ),
                    ),
            );
        }
        anterior = Some(posicion);
    }
    let directos = encontrados
        .iter()
        .filter(|(_, kmero)| kmero.canonico() == kmero.directo)
        .count();
    narracion.agregar(
        Pizarra::nueva()
            .variable("lectura.secuencia.len()", lectura.secuencia.len())
            .variable("k-meros de la lectura", encontrados.len())
            .secuencia(&lectura.secuencia, &[])
            .paso(
                "usar(posicion + 1 - k, kmero);",
                format!(
                    "Hasta el final de la lectura salen {} k-meros: {} se guardan como se leyeron y {} como su reverso.",
                    encontrados.len(),
                    directos,
                    encontrados.len() - directos
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("cuenta", "+1")
            .paso(
                "*cuenta = cuenta.saturating_add(1);",
                "En cada tabla, la primera vez que llega un k-mero se anota con 1 y las siguientes se le suma 1.",
            ),
    );
    Ok(())
}

fn narrar_histograma(narracion: &mut Narracion, histograma: &[usize]) {
    let valle = ensamblador::primer_valle(histograma);
    let cima = pico(histograma, valle.unwrap_or(1));
    let filas_visibles = valle
        .map_or(12, |valle| valle + 3)
        .max(cima.map_or(0, |cima| cima + 1))
        .clamp(5, 30)
        .min(histograma.len().saturating_sub(1));
    let filas: Vec<Vec<String>> = (1..=filas_visibles)
        .map(|cobertura| vec![cobertura.to_string(), miles(histograma[cobertura])])
        .collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("histograma[1]", miles(histograma.get(1).copied().unwrap_or(0)))
            .tabla(
                "histograma",
                &["veces", "k-meros distintos"],
                filas.clone(),
                &[(0, 1)],
            )
            .paso(
                "barras[(cuenta as usize).min(ULTIMA_BARRA)] += 1;",
                format!(
                    "Arma el histograma: cuántos k-meros distintos aparecen 1 vez, 2 veces, 3… Un error de lectura crea k-meros que casi nadie más tiene, por eso {} aparecen una sola vez.",
                    miles(histograma.get(1).copied().unwrap_or(0))
                ),
            ),
    );
    let ultima = histograma.len().min(ULTIMA_BARRA).saturating_sub(1);
    let mut detallados = 0;
    for cobertura in 1..ultima {
        let sube = histograma[cobertura + 1] > histograma[cobertura];
        if detallados < PASOS_DEL_VALLE_DETALLADOS || sube {
            let resaltar = if cobertura < filas.len() {
                vec![(cobertura - 1, 1), (cobertura, 1)]
            } else {
                Vec::new()
            };
            let nota = if sube {
                format!(
                    "{} > {}: aquí el histograma vuelve a subir. El primer valle está en {cobertura}: a la izquierda dominan los errores, a la derecha el genoma.",
                    miles(histograma[cobertura + 1]),
                    miles(histograma[cobertura])
                )
            } else {
                format!(
                    "Compara {} veces con {}: {} ≤ {}, todavía baja.",
                    cobertura + 1,
                    cobertura,
                    miles(histograma[cobertura + 1]),
                    miles(histograma[cobertura])
                )
            };
            narracion.agregar(
                Pizarra::nueva()
                    .variable("cobertura", cobertura)
                    .variable("histograma[cobertura]", miles(histograma[cobertura]))
                    .variable("histograma[cobertura + 1]", miles(histograma[cobertura + 1]))
                    .tabla("histograma", &["veces", "k-meros distintos"], filas.clone(), &resaltar)
                    .paso(
                        "(1..ultima).find(|&cobertura| histograma[cobertura + 1] > histograma[cobertura])",
                        nota,
                    ),
            );
            detallados += 1;
        }
        if sube {
            return;
        }
    }
    narracion.agregar(Pizarra::nueva().variable("cobertura", ultima).paso(
        "(1..ultima).find(|&cobertura| histograma[cobertura + 1] > histograma[cobertura])",
        "El histograma nunca vuelve a subir: no hay valle que separe errores de genoma.",
    ));
}

fn narrar_umbral(narracion: &mut Narracion, umbral: &Umbral, histograma: &[usize]) {
    let cima = pico(histograma, umbral.minimo as usize);
    let sobre_el_pico = cima.map_or(String::new(), |cima| {
        format!(" El pico está en {cima}: esa es la cobertura típica de un k-mero real.")
    });
    let (linea, nota) = match (umbral.automatico, umbral.valle) {
        (true, Some(valle)) => (
            "minimo: valle.map_or(MINIMO_SIN_VALLE, |cobertura| cobertura as u32),",
            format!(
                "Como la cobertura mínima es 0 (automática), usa el valle: se quedan los k-meros con {valle} copias o más.{sobre_el_pico}"
            ),
        ),
        (true, None) => (
            "minimo: valle.map_or(MINIMO_SIN_VALLE, |cobertura| cobertura as u32),",
            format!(
                "Sin valle, usa {} por omisión: tira solo los k-meros vistos una vez.{sobre_el_pico}",
                umbral.minimo
            ),
        ),
        (false, _) => (
            "minimo: pedido.min(COBERTURA_MAXIMA_PEDIDA) as u32,",
            format!(
                "Pediste cobertura mínima {}, así que usa ese número en vez del valle{}.{sobre_el_pico}",
                umbral.minimo,
                umbral
                    .valle
                    .map_or(String::new(), |valle| format!(" (que está en {valle})"))
            ),
        ),
    };
    narracion.agregar(
        Pizarra::nueva()
            .variable(
                "umbral.valle",
                umbral
                    .valle
                    .map_or("ninguno".to_string(), |valle| valle.to_string()),
            )
            .variable("umbral.minimo", umbral.minimo)
            .paso(linea, nota),
    );
}

fn narrar_lectura_con_error(
    narracion: &mut Narracion,
    ruta: &Path,
    kmeros: &Kmeros,
    minimo: u32,
) -> Result<(), String> {
    let k = kmeros.k;
    for (numero, lectura) in abrir_fastq(ruta)?
        .take(LECTURAS_PARA_BUSCAR_UN_ERROR)
        .enumerate()
    {
        let lectura = lectura?;
        let mut debiles: Vec<(usize, u32)> = Vec::new();
        let mut total = 0;
        ensamblador::recorrer_kmeros(&lectura.secuencia, k, |posicion, kmero| {
            total += 1;
            let cuenta = kmeros.cuenta(kmero.canonico());
            if cuenta < minimo {
                debiles.push((posicion, cuenta));
            }
        });
        let (Some(primero), Some(ultimo)) = (debiles.first(), debiles.last()) else {
            continue;
        };
        if debiles.len() == total {
            continue;
        }
        let comun_desde = ultimo.0;
        let comun_hasta = primero.0 + k;
        let sospecha = if comun_desde < comun_hasta && comun_hasta - comun_desde <= 3 {
            if comun_hasta - comun_desde == 1 {
                format!(
                    " Todos comparten la letra {}: casi seguro, un error de lectura justo ahí.",
                    comun_desde + 1
                )
            } else {
                format!(
                    " Todos comparten las letras {} a {}: casi seguro, errores de lectura muy juntos.",
                    comun_desde + 1,
                    comun_hasta
                )
            }
        } else {
            String::new()
        };
        narracion.agregar(
            Pizarra::nueva()
                .variable("lectura", numero + 1)
                .variable("k-meros débiles", debiles.len())
                .variable("cuenta", debiles.iter().map(|par| par.1).max().unwrap_or(0))
                .secuencia(
                    &lectura.secuencia,
                    &[(primero.0, ultimo.0 + k, Marca::Comparado)],
                )
                .paso(
                    "tabla.retain(|_, cuenta| *cuenta >= minimo)",
                    format!(
                        "Mira la lectura {}: {} de sus {} k-meros aparecen menos de {minimo} veces en todo el archivo (entre las letras {} y {}).{sospecha} Esos k-meros se van a tirar.",
                        numero + 1,
                        debiles.len(),
                        total,
                        primero.0 + 1,
                        ultimo.0 + k
                    ),
                ),
        );
        return Ok(());
    }
    Ok(())
}

fn narrar_primer_unitig(narracion: &mut Narracion, kmeros: &Kmeros) {
    let k = kmeros.k;
    let Some(&canonico) = kmeros.partes.iter().flat_map(|tabla| tabla.keys()).next() else {
        narracion.agregar(Pizarra::nueva().paso(
            "if visitados.insert(semilla) {",
            "No quedó ningún k-mero sólido: con esa cobertura mínima no hay nada que unir.",
        ));
        return;
    };
    let semilla = ensamblador::desde_canonico(canonico, k);
    let mut visitados: HashSet<u64> = HashSet::new();
    visitados.insert(canonico);
    let mut secuencia = ensamblador::letras_de(semilla.directo, k);
    narracion.agregar(
        Pizarra::nueva()
            .variable("semilla", letras(semilla.directo, k))
            .variable("camino.len()", 1)
            .secuencia(&secuencia, &[(0, k, Marca::Actual)])
            .paso(
                "if visitados.insert(semilla) {",
                "Ahora une los k-meros sólidos en unitigs: tramos sin bifurcaciones. Toma un k-mero que nadie ha visitado, la semilla, y camina hacia adelante.",
            ),
    );
    let mut actual = semilla;
    let mut avances = 0;
    let motivo = loop {
        let siguientes: Vec<Kmero> = kmeros.siguientes(actual).collect();
        if siguientes.len() != 1 {
            break if siguientes.is_empty() {
                "ningún k-mero sólido continúa al último (un callejón: ahí faltan k-meros sólidos)"
                    .to_string()
            } else {
                format!(
                    "el último tiene {} sucesores posibles (una bifurcación)",
                    siguientes.len()
                )
            };
        }
        let siguiente = siguientes[0];
        let anteriores = kmeros.anteriores(siguiente).count();
        if anteriores != 1 {
            break format!("al siguiente llegan {anteriores} caminos (una confluencia)");
        }
        if !visitados.insert(siguiente.canonico()) {
            break "el siguiente ya estaba en el camino (un ciclo)".to_string();
        }
        secuencia.push(ensamblador::letra_de_base(siguiente.directo));
        avances += 1;
        if avances <= PASOS_DE_CAMINO_DETALLADOS {
            let inicio = secuencia.len() - k;
            let tablero = || {
                Pizarra::nueva()
                    .variable("actual", letras(actual.directo, k))
                    .variable("siguiente", letras(siguiente.directo, k))
                    .variable("camino.len()", avances + 1)
                    .secuencia(
                        &secuencia,
                        &[
                            (inicio - 1, inicio - 1 + k, Marca::Comparado),
                            (inicio, inicio + k, Marca::Actual),
                        ],
                    )
            };
            narracion.agregar(tablero().paso(
                "while let Some(siguiente) = unico(kmeros.siguientes(actual)) {",
                format!(
                    "De las 4 letras que podrían seguir, solo con «{}» sale un k-mero sólido: hay un único sucesor.",
                    *secuencia.last().unwrap_or(&b'N') as char
                ),
            ));
            if avances <= 2 {
                narracion.agregar(tablero().paso(
                    "if unico(kmeros.anteriores(siguiente)).is_none() {",
                    "Y a ese sucesor solo se llega desde aquí (un único antecesor): no es una confluencia, así que siguen en el mismo unitig.",
                ));
            }
            narracion.agregar(tablero().paso(
                "camino.push(siguiente);",
                format!("Lo agrega: el unitig ya mide {} letras.", secuencia.len()),
            ));
        }
        actual = siguiente;
    };
    let final_adelante = secuencia.len().saturating_sub(k);
    narracion.agregar(
        Pizarra::nueva()
            .variable("actual", letras(actual.directo, k))
            .variable("camino.len()", avances + 1)
            .secuencia(
                &secuencia,
                &[(final_adelante, final_adelante + k, Marca::Actual)],
            )
            .paso(
                "while let Some(siguiente) = unico(kmeros.siguientes(actual)) {",
                format!(
                    "Avanzó {} k-meros{} y se detiene porque {motivo}.",
                    miles(avances),
                    if avances > PASOS_DE_CAMINO_DETALLADOS {
                        " (los demás, igual que los primeros)"
                    } else {
                        ""
                    }
                ),
            ),
    );
    let mut atras = 0;
    let mut actual = semilla.girado();
    while let Some(siguiente) = ensamblador::unico(kmeros.siguientes(actual)) {
        if ensamblador::unico(kmeros.anteriores(siguiente)).is_none()
            || !visitados.insert(siguiente.canonico())
        {
            break;
        }
        atras += 1;
        actual = siguiente;
    }
    let largo = k + avances + atras;
    narracion.agregar(
        Pizarra::nueva()
            .variable("semilla.girado()", letras(semilla.reverso, k))
            .variable("hacia_atras.len()", atras)
            .variable("camino.len()", avances + atras + 1)
            .paso(
                "let hacia_atras = extender(kmeros, semilla.girado(), visitados);",
                format!(
                    "Luego camina hacia atrás: avanza desde el complemento inverso de la semilla, con las mismas reglas, y suma {} k-meros más. El primer unitig mide {} letras ({} k-meros).",
                    miles(atras),
                    miles(largo),
                    miles(avances + atras + 1)
                ),
            ),
    );
}

fn tablero_de_unitig(secuencia: &[u8], marcas: &[(usize, usize, Marca)]) -> Pizarra {
    Pizarra::nueva().secuencia(secuencia, marcas)
}

fn narrar_punta(
    narracion: &mut Narracion,
    grafo: &GrafoDeUnitigs,
    punta: &Punta,
    es_ejemplo: bool,
) {
    let k = grafo.k;
    let kmeros = punta.secuencia.len() + 1 - k;
    let marcas = [(0, punta.secuencia.len(), Marca::Actual)];
    let base = || {
        tablero_de_unitig(&punta.secuencia, &marcas)
            .variable("numero", punta.unitig)
            .variable("k-meros", kmeros)
            .variable("2 * grafo.k", 2 * k)
            .variable("unitig.cobertura", decimal(punta.cobertura, 1))
    };
    narracion.agregar(base().paso(
        "if kmeros_del_unitig(grafo, numero) >= 2 * grafo.k {",
        format!(
            "Primera limpieza: busca puntas. El unitig {} tiene {} k-meros ({} letras), menos que 2k = {}: es lo bastante corto para ser una punta.",
            punta.unitig,
            kmeros,
            punta.secuencia.len(),
            2 * k
        ),
    ));
    narracion.agregar(base().paso(
        "if !grafo.salidas(punta).is_empty() || grafo.salidas(punta.opuesto()).is_empty() {",
        "Por un extremo no lleva a ningún lado y por el otro se pega al grafo: es un callejón sin salida.",
    ));
    narracion.agregar(base().variable("rival", decimal(punta.cobertura_rival, 1)).paso(
        "let rival = cobertura_rival(grafo, punta);",
        format!(
            "En la bifurcación donde se pega, la rama hermana tiene cobertura {}; la punta, {}.",
            decimal(punta.cobertura_rival, 1),
            decimal(punta.cobertura, 1)
        ),
    ));
    let adaptador = if es_ejemplo && lleva_el_adaptador(&punta.secuencia) {
        " Además lleva el adaptador de la librería: viene de lecturas que lo leyeron al final."
    } else {
        ""
    };
    narracion.agregar(base().variable("rival", decimal(punta.cobertura_rival, 1)).paso(
        "if rival > unitig.cobertura {",
        format!(
            "{} > {}: la punta pierde. Corta y menos cubierta que su rival, casi seguro es un error al final de una lectura.{adaptador}",
            decimal(punta.cobertura_rival, 1),
            decimal(punta.cobertura, 1)
        ),
    ));
    narracion.agregar(base().paso(
        "quitar_unitig(kmeros, &grafo.unitigs[punta.unitig]);",
        format!("Borra sus {kmeros} k-meros de la tabla."),
    ));
}

fn narrar_burbuja(
    narracion: &mut Narracion,
    grafo: &GrafoDeUnitigs,
    burbuja: &Burbuja,
    es_ejemplo: bool,
) {
    let k = grafo.k;
    let largo = burbuja.quitada.len() + 1 - k;
    let diferencias: Vec<(usize, usize, Marca)> =
        if burbuja.quitada.len() == burbuja.conservada.len() {
            burbuja
                .quitada
                .iter()
                .zip(&burbuja.conservada)
                .enumerate()
                .filter(|(_, (una, otra))| una != otra)
                .map(|(posicion, _)| (posicion, posicion + 1, Marca::Comparado))
                .collect()
        } else {
            Vec::new()
        };
    let tabla = vec![
        vec![
            "rama débil".to_string(),
            miles(burbuja.quitada.len()),
            decimal(burbuja.cobertura_quitada, 1),
            recortar(&burbuja.quitada, VISIBLE),
        ],
        vec![
            format!(
                "otro camino ({})",
                cantidad(burbuja.unitigs_conservados, "unitig", "unitigs")
            ),
            miles(burbuja.conservada.len()),
            decimal(burbuja.cobertura_conservada, 1),
            recortar(&burbuja.conservada, VISIBLE),
        ],
    ];
    let columnas = ["rama", "letras", "cobertura", "secuencia"];
    let base = || {
        tablero_de_unitig(&burbuja.quitada, &diferencias)
            .variable("numero", burbuja.unitig_quitado)
            .variable("largo", largo)
            .variable(
                "rama_debil.cobertura",
                decimal(burbuja.cobertura_quitada, 1),
            )
    };
    narracion.agregar(base().paso(
        "let largo = kmeros_del_unitig(grafo, numero);",
        format!(
            "Segunda limpieza: burbujas. El unitig {} ({} k-meros, cobertura {}) sale de un punto del grafo y llega a otro: puede ser una rama de burbuja.",
            burbuja.unitig_quitado,
            largo,
            decimal(burbuja.cobertura_quitada, 1)
        ),
    ));
    narracion.agregar(base().paso(
        "busqueda.explorar(grafo, origen, &mut Vec::new(), 0, f64::INFINITY);",
        "Busca, con una búsqueda en profundidad acotada, otro camino desde el mismo punto de salida hasta el mismo punto de llegada, de largo parecido (±10 %).",
    ));
    narracion.agregar(
        base()
            .variable("cobertura", decimal(burbuja.cobertura_conservada, 1))
            .tabla("las dos ramas", &columnas, tabla.clone(), &[(1, 2)])
            .paso(
                "if self.destinos.contains(&siguiente) {",
                format!(
                    "Lo encuentra: {} unitig{} con {} letras y cobertura mínima {}.{}",
                    burbuja.unitigs_conservados,
                    if burbuja.unitigs_conservados == 1 {
                        ""
                    } else {
                        "s"
                    },
                    burbuja.conservada.len(),
                    decimal(burbuja.cobertura_conservada, 1),
                    if diferencias.is_empty() {
                        String::new()
                    } else {
                        format!(
                            " Las dos ramas difieren en {} letra{} (marcadas).",
                            diferencias.len(),
                            if diferencias.len() == 1 { "" } else { "s" }
                        )
                    }
                ),
            ),
    );
    let origen = if es_ejemplo {
        " En el ejemplo, ambas ramas suelen existir en el genoma: son copias de una repetición que difieren en pocas letras."
    } else {
        ""
    };
    narracion.agregar(
        base()
            .variable("cobertura", decimal(burbuja.cobertura_conservada, 1))
            .tabla("las dos ramas", &columnas, tabla, &[(1, 2), (0, 2)])
            .paso(
                "if cobertura > rama_debil.cobertura {",
                format!(
                    "{} > {}: gana el camino más cubierto y la rama débil se marca para quitar.{origen}",
                    decimal(burbuja.cobertura_conservada, 1),
                    decimal(burbuja.cobertura_quitada, 1)
                ),
            ),
    );
    narracion.agregar(base().paso(
        "quitar_unitig(kmeros, &grafo.unitigs[burbuja.unitig_quitado]);",
        format!("Borra los {largo} k-meros de la rama débil."),
    ));
}

#[cfg(test)]
mod pruebas {
    use super::*;
    use nucleo::formatos;

    fn pedido_de_ejemplo() -> Pedido {
        leer_pedido(&Entrada::nueva(&[
            ("lecturas", "@ejemplo:lecturas"),
            ("k", "31"),
            ("cobertura_minima", "0"),
        ]))
        .unwrap()
    }

    #[test]
    fn el_paso_a_paso_termina_igual_que_nucleo() {
        let pedido = pedido_de_ejemplo();
        let (pasos, resumen) = narrar(&pedido).unwrap();
        let real = ensamblador::ensamblar(&pedido.ruta, pedido.k, pedido.minimo).unwrap();
        assert_eq!(
            resumen,
            Final {
                contigs: real.contigs.len(),
                n50: real.medida.n50,
                total: real.medida.total,
                mayor: real.contigs[0].secuencia.len(),
                puntas: real.puntas.len(),
                burbujas: real.burbujas.len(),
                rondas: real.rondas,
            }
        );
        assert!(pasos.len() >= 20 && pasos.len() <= MAXIMO_PASOS);
        let texto = serde_json::to_string(&pasos).unwrap();
        assert!(texto.len() < 2_000_000);
        assert!(pasos.iter().any(|paso| paso.linea.starts_with("let rival")));
        assert!(
            pasos
                .iter()
                .any(|paso| paso.linea.starts_with("busqueda.explorar"))
        );
    }

    #[test]
    fn el_primer_unitig_narrado_esta_en_el_grafo() {
        let pedido = pedido_de_ejemplo();
        let mut conteo = ensamblador::contar_kmeros(&pedido.ruta, pedido.k).unwrap();
        let histograma = ensamblador::histograma(&conteo.kmeros);
        let umbral = ensamblador::elegir_umbral(&histograma, 0);
        conteo.kmeros.retener(umbral.minimo);
        let mut narracion = Narracion::nueva(MAXIMO_PASOS);
        narrar_primer_unitig(&mut narracion, &conteo.kmeros);
        let pasos = narracion.terminar(Pizarra::nueva().paso("x", "fin"));
        let nota = &pasos[pasos.len() - 2].nota;
        let grafo = ensamblador::compactar(&conteo.kmeros);
        let largos: HashSet<String> = grafo
            .unitigs
            .iter()
            .map(|unitig| miles(unitig.secuencia.len()))
            .collect();
        let dicho = nota
            .split("El primer unitig mide ")
            .nth(1)
            .and_then(|resto| resto.split(" letras").next())
            .unwrap()
            .to_string();
        assert!(largos.contains(&dicho), "{dicho}");
    }

    #[test]
    fn el_ejemplo_cubre_las_ventanas_verdaderas() {
        let salida = ensamblar(&Entrada::nueva(&[
            ("lecturas", "@ejemplo:lecturas"),
            ("k", "31"),
            ("cobertura_minima", "0"),
        ]))
        .unwrap();
        let valor = |buscada: &str| {
            salida
                .bloques
                .iter()
                .find_map(|bloque| match bloque {
                    Bloque::Valor { etiqueta, valor } if etiqueta.starts_with(buscada) => {
                        Some(valor.clone())
                    }
                    _ => None,
                })
                .unwrap()
        };
        let cubiertas = valor("Ventanas verdaderas");
        let porcentaje: f64 = cubiertas
            .split('(')
            .nth(1)
            .and_then(|resto| resto.split('\u{a0}').next())
            .unwrap()
            .parse()
            .unwrap();
        assert!(porcentaje > 95.0, "{cubiertas}");
        assert!(valor("Mutaciones conocidas").starts_with("2"));
        assert!(salida.interpretacion.len() >= 4);
        assert!(
            salida
                .bloques
                .iter()
                .any(|bloque| matches!(bloque, Bloque::Archivo { .. }))
        );
    }

    #[test]
    fn narra_archivos_pequenos_sin_puntas_ni_burbujas() {
        let ruta =
            std::env::temp_dir().join(format!("ensamblador-catalogo-{}.fastq", std::process::id()));
        let mut texto = String::new();
        let genoma = "ACGTTGCAAGGCTTAACCGGATCGATCGGCTAGCTTAGGCATCGATCGACTGACTAGCTAGCATCG";
        for inicio in 0..20 {
            let lectura = &genoma[inicio..inicio + 40];
            for _ in 0..3 {
                texto.push_str(&format!("@l{inicio}\n{lectura}\n+\n{}\n", "I".repeat(40)));
            }
        }
        fs::write(&ruta, texto).unwrap();
        let entrada = Entrada::nueva(&[
            ("lecturas", ruta.to_str().unwrap()),
            ("k", "15"),
            ("cobertura_minima", "2"),
        ]);
        let pasos = pasos_ensamblar(&entrada).unwrap();
        let salida = ensamblar(&entrada).unwrap();
        let destino = crate::ruta_de_salida(&ruta, "contigs", "fasta");
        let contigs = formatos::leer_fasta_archivo(&destino).unwrap();
        fs::remove_file(&ruta).ok();
        fs::remove_file(&destino).ok();
        assert_eq!(contigs.len(), 1);
        assert_eq!(contigs[0].secuencia.len(), 59);
        assert!(pasos.len() > 5);
        assert!(!salida.interpretacion.is_empty());
    }

    fn lecturas_pequenas(nombre: &str) -> PathBuf {
        let ruta = std::env::temp_dir().join(format!(
            "ensamblador-catalogo-{}-{nombre}.fastq",
            std::process::id()
        ));
        let genoma = "ACGTTGCAAGGCTTAACCGGATCGATCGGCTAGCTTAGGCATCGATCGACTGACTAGCTAGCATCG";
        let mut texto = String::new();
        for inicio in 0..20 {
            let lectura = &genoma[inicio..inicio + 40];
            for _ in 0..3 {
                texto.push_str(&format!("@l{inicio}\n{lectura}\n+\n{}\n", "I".repeat(40)));
            }
        }
        fs::write(&ruta, texto).unwrap();
        ruta
    }

    #[test]
    fn el_pico_no_se_inventa_si_la_ultima_barra_manda() {
        let mut histograma = vec![0usize; ULTIMA_BARRA + 1];
        histograma[1] = 100;
        histograma[2] = 10;
        histograma[3] = 20;
        histograma[ULTIMA_BARRA] = 500;
        assert_eq!(pico(&histograma, 3), None);
        histograma[ULTIMA_BARRA] = 5;
        assert_eq!(pico(&histograma, 3), Some(3));
    }

    #[test]
    fn la_narracion_de_la_primera_lectura_dice_donde_empieza_cada_kmero() {
        let ruta = std::env::temp_dir().join(format!(
            "ensamblador-catalogo-{}-enes.fastq",
            std::process::id()
        ));
        let lectura = "NACGTTGCAAGGNTTAACCGGATCGATCGG";
        fs::write(
            &ruta,
            format!("@a\n{lectura}\n+\n{}\n", "I".repeat(lectura.len())),
        )
        .unwrap();
        let mut narracion = Narracion::nueva(MAXIMO_PASOS);
        narrar_primera_lectura(&mut narracion, &ruta, 11).unwrap();
        fs::remove_file(&ruta).ok();
        let pasos = narracion.terminar(Pizarra::nueva().paso("x", "fin"));
        let notas: Vec<&str> = pasos.iter().map(|paso| paso.nota.as_str()).collect();
        assert!(
            notas
                .iter()
                .any(|nota| nota.contains("de la letra 2 a la 12"))
        );
        assert!(notas.iter().any(|nota| nota.contains("vuelve a juntar")));
    }

    #[test]
    fn no_asegura_que_las_mutaciones_no_cortan_si_no_es_cierto() {
        let ruta = lecturas_pequenas("contraste");
        let resultado = ensamblador::ensamblar(&ruta, 15, 2).unwrap();
        fs::remove_file(&ruta).ok();
        let contraste = Contraste {
            comparacion: Comparacion::default(),
            mutaciones: vec![(1, 2); 4],
            mutaciones_en_contigs: 1,
            huecos_en_mutaciones: 0,
            puntas_con_adaptador: 0,
            kmeros_de_burbujas: 0,
            kmeros_de_burbujas_ajenos: 0,
        };
        let mut salida = Salida::default();
        explicar_contraste(&mut salida, &resultado, &contraste);
        let texto = salida.interpretacion.join(" ");
        assert!(texto.contains("solo 1 de las 4 mutaciones"), "{texto}");
        assert!(!texto.contains("no lo cortan"), "{texto}");
    }
}
