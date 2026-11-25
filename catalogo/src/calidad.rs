use crate::formato::{decimal, miles, porcentaje};
use crate::pizarra::{Marca as Resalte, Narracion, Pizarra};
use crate::{Algoritmo, Bloque, Campo, Entrada, Marca, Paso, Salida, Tipo};
use nucleo::calidad::{
    self, Descarte, InformeCalidad, ModoAdaptador, OpcionesRecorte, ResumenRecorte,
};
use nucleo::formatos::{Lectura, abrir_fastq};
use nucleo::{ecoli, secuencia};
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static CONTADOR: AtomicU64 = AtomicU64::new(0);

const BASES: [&str; 5] = ["A", "C", "G", "T", "N"];
const COLUMNAS_DE_COMPOSICION: usize = 150;
const BARRAS_DE_LONGITUD: usize = 60;

const ARCHIVO: Campo = Campo {
    id: "archivo",
    etiqueta: "Archivo FASTQ (.fastq o .fastq.gz)",
    tipo: Tipo::Archivo,
    ejemplo: "@ejemplo:lecturas",
};

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![
        Algoritmo {
            id: "puerto.calidad_lecturas",
            capitulo: "puerto",
            nombre: "Calidad de lecturas",
            resumen: "Revisa un FASTQ: calidad por posición, composición, GC, duplicados, secuencias repetidas y adaptadores.",
            campos: &[ARCHIVO],
            funcion: "calidad::analizar_fastq",
            ejecutar: calidad_de_lecturas,
            pasos: Some(pasos_calidad_de_lecturas),
        },
        Algoritmo {
            id: "puerto.recortar_lecturas",
            capitulo: "puerto",
            nombre: "Recortar lecturas",
            resumen: "Quita adaptadores y extremos de baja calidad, descarta lecturas malas y guarda un FASTQ limpio.",
            campos: &[
                ARCHIVO,
                Campo {
                    id: "adaptador",
                    etiqueta: "Adaptador («auto» lo detecta, «no» no lo busca, o escribe su secuencia)",
                    tipo: Tipo::Texto,
                    ejemplo: "auto",
                },
                Campo {
                    id: "calidad",
                    etiqueta: "Calidad media mínima en la ventana (Phred)",
                    tipo: Tipo::Entero,
                    ejemplo: "20",
                },
                Campo {
                    id: "ventana",
                    etiqueta: "Bases de la ventana deslizante",
                    tipo: Tipo::Entero,
                    ejemplo: "4",
                },
                Campo {
                    id: "longitud",
                    etiqueta: "Longitud mínima para conservar una lectura (pb)",
                    tipo: Tipo::Entero,
                    ejemplo: "50",
                },
                Campo {
                    id: "maximo_n",
                    etiqueta: "Máximo de bases N por lectura",
                    tipo: Tipo::Entero,
                    ejemplo: "5",
                },
            ],
            funcion: "calidad::recortar_fastq",
            ejecutar: recortar_lecturas,
            pasos: Some(pasos_recortar_lecturas),
        },
    ]
}

fn calidad_de_lecturas(entrada: &Entrada) -> Result<Salida, String> {
    let ruta = entrada.ruta("archivo")?;
    let informe = calidad::analizar_fastq(&ruta)?;
    if informe.lecturas == 0 {
        return Err("El archivo no tiene ninguna lectura FASTQ.".to_string());
    }
    let mut salida = Salida::default();
    agregar_totales(&mut salida, &informe);
    agregar_graficos(&mut salida, &informe, "");
    agregar_tablas(&mut salida, &informe);
    explicar_calidad(&mut salida, &informe);
    Ok(salida)
}

fn agregar_totales(salida: &mut Salida, informe: &InformeCalidad) {
    salida.agregar(Bloque::valor("Lecturas", miles(informe.lecturas)));
    salida.agregar(Bloque::valor("Bases", miles(informe.bases)));
    salida.agregar(Bloque::valor(
        "Longitud de las lecturas",
        rango_de_longitudes(informe),
    ));
    salida.agregar(Bloque::valor(
        "Bases con calidad Q30 o más",
        porcentaje(informe.fraccion_q30(), 1),
    ));
    salida.agregar(Bloque::valor(
        "Bases con calidad Q20 o más",
        porcentaje(informe.fraccion_q20(), 1),
    ));
    salida.agregar(Bloque::valor("GC medio", porcentaje(informe.gc_medio(), 1)));
    salida.agregar(Bloque::valor(
        "Bases N (no leídas)",
        porcentaje(informe.fraccion_n(), 2),
    ));
    salida.agregar(Bloque::valor(
        &format!(
            "Lecturas duplicadas (en las primeras {})",
            miles(informe.lecturas_muestreadas)
        ),
        porcentaje(informe.duplicacion, 1),
    ));
}

fn rango_de_longitudes(informe: &InformeCalidad) -> String {
    if informe.longitud_minima == informe.longitud_maxima {
        format!("{} pb", miles(informe.longitud_maxima))
    } else {
        format!(
            "{}–{} pb",
            miles(informe.longitud_minima),
            miles(informe.longitud_maxima)
        )
    }
}

fn agregar_graficos(salida: &mut Salida, informe: &InformeCalidad, cuando: &str) {
    let medias: Vec<[f64; 2]> = informe
        .calidad_por_posicion
        .iter()
        .enumerate()
        .filter(|(_, calidad)| calidad.lecturas > 0)
        .map(|(posicion, calidad)| [(posicion + 1) as f64, calidad.media])
        .collect();
    let mut marcas = Vec::new();
    if let Some(posicion) = primera_posicion_bajo(informe, 30.0) {
        marcas.push(Marca {
            x: (posicion + 1) as f64,
            texto: "Media bajo Q30".to_string(),
        });
    }
    if let Some(posicion) = inicio_de_cola_baja(informe, 20.0) {
        marcas.push(Marca {
            x: (posicion + 1) as f64,
            texto: "Media bajo Q20 hasta el final".to_string(),
        });
    }
    salida.agregar(Bloque::linea(
        &format!("Calidad media por posición{cuando}"),
        "Posición en la lectura (contando desde 1)",
        "Calidad media (Phred)",
        &medias,
        marcas,
    ));
    let cuartiles: Vec<[f64; 2]> = informe
        .calidad_por_posicion
        .iter()
        .enumerate()
        .filter(|(_, calidad)| calidad.lecturas > 0)
        .map(|(posicion, calidad)| [(posicion + 1) as f64, calidad.cuartil_inferior as f64])
        .collect();
    salida.agregar(Bloque::linea(
        &format!("Cuartil inferior de calidad por posición{cuando} (una de cada cuatro bases está por debajo)"),
        "Posición en la lectura (contando desde 1)",
        "Calidad (Phred)",
        &cuartiles,
        Vec::new(),
    ));
    agregar_composicion(salida, informe, cuando);
    salida.agregar(Bloque::barras(
        &format!("Lecturas según su contenido GC{cuando}"),
        "GC de la lectura (%)",
        "Lecturas (%)",
        (0..informe.lecturas_por_gc.len())
            .map(|gc| gc.to_string())
            .collect(),
        informe
            .lecturas_por_gc
            .iter()
            .map(|&veces| 100.0 * veces as f64 / informe.lecturas as f64)
            .collect(),
    ));
    agregar_longitudes(salida, informe, cuando);
}

fn agregar_composicion(salida: &mut Salida, informe: &InformeCalidad, cuando: &str) {
    let posiciones = informe.composicion_por_posicion.len();
    let ancho = posiciones.div_ceil(COLUMNAS_DE_COMPOSICION).max(1);
    let mut columnas = Vec::new();
    let mut valores = vec![Vec::new(); BASES.len()];
    for inicio in (0..posiciones).step_by(ancho) {
        let fin = (inicio + ancho).min(posiciones);
        columnas.push(if fin - inicio == 1 {
            (inicio + 1).to_string()
        } else {
            format!("{}–{}", inicio + 1, fin)
        });
        let mut por_base = [0usize; 5];
        for conteo in &informe.composicion_por_posicion[inicio..fin] {
            for (letra, veces) in conteo.iter().enumerate() {
                por_base[letra] += veces;
            }
        }
        let total: usize = por_base.iter().sum();
        for (letra, fila) in valores.iter_mut().enumerate() {
            fila.push(if total == 0 {
                0.0
            } else {
                100.0 * por_base[letra] as f64 / total as f64
            });
        }
    }
    salida.agregar(Bloque::calor(
        &format!("Composición de bases por posición{cuando} (% de las lecturas, posiciones contando desde 1)"),
        BASES.iter().map(|base| base.to_string()).collect(),
        columnas,
        valores,
    ));
}

fn agregar_longitudes(salida: &mut Salida, informe: &InformeCalidad, cuando: &str) {
    let minima = informe.longitud_minima;
    let maxima = informe.longitud_maxima;
    let ancho = (maxima - minima + 1).div_ceil(BARRAS_DE_LONGITUD).max(1);
    let mut categorias = Vec::new();
    let mut valores = Vec::new();
    for inicio in (minima..=maxima).step_by(ancho) {
        let fin = (inicio + ancho - 1).min(maxima);
        let lecturas: usize = informe.lecturas_por_longitud[inicio..=fin].iter().sum();
        categorias.push(if inicio == fin {
            inicio.to_string()
        } else {
            format!("{inicio}–{fin}")
        });
        valores.push(100.0 * lecturas as f64 / informe.lecturas as f64);
    }
    salida.agregar(Bloque::barras(
        &format!("Lecturas por longitud{cuando}"),
        "Longitud (pb)",
        "Lecturas (%)",
        categorias,
        valores,
    ));
}

fn agregar_tablas(salida: &mut Salida, informe: &InformeCalidad) {
    match &informe.adaptador {
        Some(adaptador) => salida.agregar(Bloque::tabla(
            "Adaptador encontrado en los datos",
            &[
                "Secuencia",
                "Longitud (pb)",
                "Lecturas con 10 o más bases del adaptador",
                "Lecturas revisadas",
                "% de las revisadas",
            ],
            vec![vec![
                texto(&adaptador.secuencia),
                adaptador.secuencia.len().to_string(),
                miles(adaptador.lecturas),
                miles(adaptador.revisadas),
                porcentaje(adaptador.fraccion, 2),
            ]],
        )),
        None => salida.agregar(Bloque::valor(
            "Adaptador encontrado en los datos",
            "ninguno",
        )),
    }
    if informe.sobrerrepresentadas.is_empty() {
        salida.agregar(Bloque::valor(
            "Secuencias sobrerrepresentadas (más del 0.1 % de la muestra)",
            "ninguna",
        ));
        return;
    }
    let filas = informe
        .sobrerrepresentadas
        .iter()
        .map(|repetida| {
            vec![
                texto(&repetida.secuencia),
                miles(repetida.veces),
                porcentaje(repetida.fraccion, 2),
                posible_origen(&repetida.secuencia, informe).to_string(),
            ]
        })
        .collect();
    salida.agregar(Bloque::tabla(
        &format!(
            "Secuencias sobrerrepresentadas (más del 0.1 % de las primeras {} lecturas)",
            miles(informe.lecturas_muestreadas)
        ),
        &["Secuencia", "Veces", "% de la muestra", "Posible origen"],
        filas,
    ));
}

fn posible_origen(repetida: &[u8], informe: &InformeCalidad) -> &'static str {
    if let Some(adaptador) = &informe.adaptador
        && calidad::buscar_adaptador(repetida, &adaptador.secuencia, calidad::LONGITUD_SEMILLA)
            .is_some()
    {
        return "Contiene el adaptador";
    }
    let mut letras = [0usize; 256];
    for &base in repetida {
        letras[base as usize] += 1;
    }
    let mas_comun = letras.iter().max().copied().unwrap_or(0);
    if mas_comun * 10 >= repetida.len() * 9 {
        return "Poco variada (por ejemplo poli-G: el secuenciador no leyó señal)";
    }
    "Sin identificar (contaminación, ARN ribosómico o duplicados de PCR)"
}

fn texto(secuencia: &[u8]) -> String {
    String::from_utf8_lossy(secuencia).to_string()
}

fn posiciones_confiables(informe: &InformeCalidad) -> impl Iterator<Item = (usize, f64)> + '_ {
    let minimo = (informe.lecturas / 100).max(1);
    informe
        .calidad_por_posicion
        .iter()
        .enumerate()
        .filter(move |(_, calidad)| calidad.lecturas >= minimo)
        .map(|(posicion, calidad)| (posicion, calidad.media))
}

fn primera_posicion_bajo(informe: &InformeCalidad, umbral: f64) -> Option<usize> {
    posiciones_confiables(informe)
        .find(|&(_, media)| media < umbral)
        .map(|(posicion, _)| posicion)
}

fn inicio_de_cola_baja(informe: &InformeCalidad, umbral: f64) -> Option<usize> {
    let confiables: Vec<(usize, f64)> = posiciones_confiables(informe).collect();
    let mut inicio = None;
    for &(posicion, media) in confiables.iter().rev() {
        if media >= umbral {
            break;
        }
        inicio = Some(posicion);
    }
    inicio
}

fn explicar_calidad(salida: &mut Salida, informe: &InformeCalidad) {
    let q30 = informe.fraccion_q30();
    if q30 >= 0.8 {
        salida.explicar(format!(
            "El {} de las bases tiene calidad Q30 o más (a lo sumo 1 error cada 1000 bases): es una buena corrida para datos tipo Illumina, donde se espera 80 % o más.",
            porcentaje(q30, 1)
        ));
    } else if q30 >= 0.7 {
        salida.explicar(format!(
            "El {} de las bases tiene calidad Q30 o más: es aceptable (entre 70 % y 80 %), pero mira en qué posiciones cae la calidad.",
            porcentaje(q30, 1)
        ));
    } else {
        salida.explicar(format!(
            "Solo el {} de las bases tiene calidad Q30 o más (1 error cada 1000 bases): es bajo para datos tipo Illumina (lo bueno es 80 % o más), así que habrá más errores de lectura.",
            porcentaje(q30, 1)
        ));
    }
    explicar_caida(salida, informe);
    if informe.longitud_maxima > calidad::MAXIMO_POSICIONES {
        salida.explicar(format!(
            "Hay lecturas de más de {} pb: las gráficas por posición cubren solo las primeras {} bases de cada lectura.",
            miles(calidad::MAXIMO_POSICIONES),
            miles(calidad::MAXIMO_POSICIONES)
        ));
    }
    match &informe.adaptador {
        Some(adaptador) if adaptador.fraccion > 0.01 => salida.explicar(format!(
            "El adaptador {} aparece en el {} de las lecturas (más de 1 %): son bases artificiales que no están en el genoma. Recórtalo con «Recortar lecturas» antes de mapear o ensamblar.",
            texto(&adaptador.secuencia),
            porcentaje(adaptador.fraccion, 1)
        )),
        Some(adaptador) => salida.explicar(format!(
            "Hay adaptador en solo el {} de las lecturas (1 % o menos): es poco, aunque recortarlo no hace daño.",
            porcentaje(adaptador.fraccion, 2)
        )),
        None => salida.explicar(format!(
            "No se encontró ningún adaptador en las primeras {} lecturas: casi todos los fragmentos de ADN eran más largos que las lecturas, así que no hace falta recortar adaptadores.",
            miles(informe.lecturas.min(calidad::LECTURAS_PARA_ADAPTADOR))
        )),
    }
    if informe.fraccion_n() > 0.01 {
        salida.explicar(format!(
            "El {} de las bases son N (más de 1 %): el secuenciador no pudo decidir muchas bases, lo que suele indicar un problema en la corrida.",
            porcentaje(informe.fraccion_n(), 1)
        ));
    }
    explicar_duplicacion(salida, informe);
    explicar_gc(salida, informe);
}

fn explicar_caida(salida: &mut Salida, informe: &InformeCalidad) {
    let Some(primera) = informe
        .calidad_por_posicion
        .iter()
        .find(|calidad| calidad.lecturas > 0)
    else {
        return;
    };
    let ultima = posiciones_confiables(informe)
        .last()
        .map(|(_, media)| media)
        .unwrap_or(primera.media);
    if let Some(posicion) = inicio_de_cola_baja(informe, 20.0) {
        salida.explicar(format!(
            "La calidad media cae por debajo de Q20 (1 error cada 100 bases) desde la posición {} (contando desde 1) hasta el final: recorta el extremo 3' con «Recortar lecturas».",
            posicion + 1
        ));
    } else if let Some(posicion) = primera_posicion_bajo(informe, 20.0) {
        salida.explicar(format!(
            "La calidad media baja de Q20 en medio de la lectura (posición {}, contando desde 1) y luego se recupera: puede ser una burbuja o un problema puntual del instrumento.",
            posicion + 1
        ));
    } else if let Some(posicion) = primera_posicion_bajo(informe, 30.0) {
        salida.explicar(format!(
            "La calidad media nunca baja de Q20: va de Q{} al inicio a Q{} al final y queda bajo Q30 desde la posición {} (contando desde 1). Es la caída normal hacia el extremo 3'.",
            decimal(primera.media, 0),
            decimal(ultima, 0),
            posicion + 1
        ));
    } else {
        salida.explicar("La calidad media se mantiene en Q30 o más en toda la lectura.");
    }
}

fn explicar_duplicacion(salida: &mut Salida, informe: &InformeCalidad) {
    let repetidas = if informe.sobrerrepresentadas.is_empty() {
        String::new()
    } else {
        format!(
            " Además, {} secuencia(s) aparecen en más del 0.1 % de las lecturas (mira la tabla).",
            informe.sobrerrepresentadas.len()
        )
    };
    let muestra = miles(informe.lecturas_muestreadas);
    if informe.duplicacion > 0.2 {
        salida.explicar(format!(
            "El {} de las primeras {muestra} lecturas son copias exactas de otra: para ADN genómico es mucho (más de 20 %), probablemente por duplicados de PCR o poco ADN de partida.{repetidas}",
            porcentaje(informe.duplicacion, 1)
        ));
    } else {
        salida.explicar(format!(
            "Solo el {} de las primeras {muestra} lecturas son copias exactas de otra (menos de 20 %): la librería es diversa.{repetidas}",
            porcentaje(informe.duplicacion, 1)
        ));
    }
}

fn explicar_gc(salida: &mut Salida, informe: &InformeCalidad) {
    let gc = informe.gc_medio();
    let gc_ecoli = secuencia::contenido_gc(ecoli::genoma());
    if let Some((primero, segundo)) = dos_picos(&informe.lecturas_por_gc) {
        salida.explicar(format!(
            "El GC por lectura tiene dos picos (cerca de {primero} % y {segundo} %): un solo organismo da un solo pico, así que puede haber contaminación con otro."
        ));
    } else if (gc - gc_ecoli).abs() <= 0.05 {
        salida.explicar(format!(
            "El GC medio ({}) está cerca del de E. coli ({}) y la distribución tiene un solo pico: no hay señales de contaminación.",
            porcentaje(gc, 1),
            porcentaje(gc_ecoli, 1)
        ));
    } else {
        salida.explicar(format!(
            "El GC medio ({}) se aleja más de 5 puntos del de E. coli ({}): si tu muestra es de E. coli, sospecha contaminación; si es de otro organismo, compáralo con su GC esperado.",
            porcentaje(gc, 1),
            porcentaje(gc_ecoli, 1)
        ));
    }
}

fn dos_picos(lecturas_por_gc: &[usize]) -> Option<(usize, usize)> {
    let suavizada: Vec<f64> = (0..lecturas_por_gc.len())
        .map(|centro| {
            let desde = centro.saturating_sub(3);
            let hasta = (centro + 3).min(lecturas_por_gc.len() - 1);
            let suma: usize = lecturas_por_gc[desde..=hasta].iter().sum();
            suma as f64 / (hasta - desde + 1) as f64
        })
        .collect();
    let (principal, &altura_principal) = suavizada
        .iter()
        .enumerate()
        .max_by(|una, otra| una.1.total_cmp(otra.1))?;
    let mut mejor: Option<(usize, f64)> = None;
    for candidato in 1..suavizada.len().saturating_sub(1) {
        let altura = suavizada[candidato];
        let es_pico = altura > suavizada[candidato - 1] && altura >= suavizada[candidato + 1];
        if !es_pico || altura < 0.1 * altura_principal || candidato.abs_diff(principal) < 10 {
            continue;
        }
        let desde = candidato.min(principal);
        let hasta = candidato.max(principal);
        let valle = suavizada[desde..=hasta]
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min);
        if valle < 0.7 * altura && mejor.is_none_or(|(_, otra)| altura > otra) {
            mejor = Some((candidato, altura));
        }
    }
    let (secundario, _) = mejor?;
    Some((principal.min(secundario), principal.max(secundario)))
}

fn leer_adaptador(entrada: &Entrada) -> Result<ModoAdaptador, String> {
    let texto = entrada.texto("adaptador")?;
    let limpio = texto.trim();
    if limpio.eq_ignore_ascii_case("auto") {
        return Ok(ModoAdaptador::Detectar);
    }
    if limpio.eq_ignore_ascii_case("no") || limpio.eq_ignore_ascii_case("ninguno") {
        return Ok(ModoAdaptador::Ninguno);
    }
    let secuencia: Vec<u8> = limpio
        .bytes()
        .filter(|letra| !letra.is_ascii_whitespace())
        .map(|letra| letra.to_ascii_uppercase())
        .collect();
    if let Some(&letra) = secuencia
        .iter()
        .find(|letra| !matches!(letra, b'A' | b'C' | b'G' | b'T'))
    {
        return Err(format!(
            "El adaptador solo admite A, C, G y T (o «auto» para detectarlo, «no» para no buscarlo), pero contiene «{}».",
            letra as char
        ));
    }
    if !(6..=200).contains(&secuencia.len()) {
        return Err("El adaptador debe tener entre 6 y 200 bases.".to_string());
    }
    Ok(ModoAdaptador::Conocido(secuencia))
}

fn leer_opciones(entrada: &Entrada) -> Result<OpcionesRecorte, String> {
    let adaptador = leer_adaptador(entrada)?;
    let calidad_minima = entrada.natural("calidad")?;
    let ventana = entrada.natural("ventana")?;
    let longitud_minima = entrada.natural("longitud")?;
    let maximo_n = entrada.natural("maximo_n")?;
    if calidad_minima > 60 {
        return Err("La calidad mínima debe estar entre 0 y 60 (escala Phred).".to_string());
    }
    if !(1..=100).contains(&ventana) {
        return Err("La ventana debe tener entre 1 y 100 bases.".to_string());
    }
    if !(1..=100_000).contains(&longitud_minima) {
        return Err("La longitud mínima debe estar entre 1 y 100 000 pb.".to_string());
    }
    Ok(OpcionesRecorte {
        adaptador,
        ventana,
        calidad_minima,
        maximo_n,
        longitud_minima,
    })
}

fn ruta_temporal(destino: &Path) -> PathBuf {
    let numero = CONTADOR.fetch_add(1, Ordering::Relaxed);
    let nombre = destino
        .file_name()
        .map(|nombre| nombre.to_string_lossy().to_string())
        .unwrap_or_else(|| "recortadas.fastq".to_string());
    destino.with_file_name(format!(".{}-{numero}-{nombre}", std::process::id()))
}

fn recortar_y_revisar(
    ruta: &Path,
    temporal: &Path,
    opciones: &OpcionesRecorte,
) -> Result<(ResumenRecorte, InformeCalidad), String> {
    let resumen = calidad::recortar_fastq(ruta, temporal, opciones)?;
    if resumen.lecturas_entrada == 0 {
        return Err("El archivo no tiene ninguna lectura FASTQ.".to_string());
    }
    let informe = calidad::analizar_fastq(temporal)?;
    Ok((resumen, informe))
}

fn recortar_lecturas(entrada: &Entrada) -> Result<Salida, String> {
    let ruta = entrada.ruta("archivo")?;
    let opciones = leer_opciones(entrada)?;
    let extension = if ruta.to_string_lossy().ends_with(".gz") {
        "fastq.gz"
    } else {
        "fastq"
    };
    let destino = crate::ruta_de_salida(&ruta, "recortadas", extension);
    let temporal = ruta_temporal(&destino);
    let resultado = recortar_y_revisar(&ruta, &temporal, &opciones).and_then(|revisado| {
        fs::rename(&temporal, &destino)
            .map_err(|error| format!("No se pudo guardar {}: {error}", destino.display()))?;
        Ok(revisado)
    });
    let (resumen, informe) = match resultado {
        Ok(revisado) => revisado,
        Err(error) => {
            let _ = fs::remove_file(&temporal);
            return Err(error);
        }
    };
    let mut salida = Salida::default();
    salida.agregar(Bloque::archivo("Lecturas recortadas (FASTQ)", &destino));
    agregar_resumen_de_recorte(&mut salida, &resumen, &opciones);
    if informe.lecturas > 0 {
        agregar_graficos(&mut salida, &informe, " después de recortar");
    }
    explicar_recorte(&mut salida, &resumen, &informe, &opciones);
    Ok(salida)
}

fn de_un_total(parte: usize, total: usize) -> String {
    let fraccion = if total == 0 {
        0.0
    } else {
        parte as f64 / total as f64
    };
    format!(
        "{} de {} ({})",
        miles(parte),
        miles(total),
        porcentaje(fraccion, 1)
    )
}

fn agregar_resumen_de_recorte(
    salida: &mut Salida,
    resumen: &ResumenRecorte,
    opciones: &OpcionesRecorte,
) {
    salida.agregar(Bloque::valor(
        "Lecturas conservadas",
        de_un_total(resumen.lecturas_salida, resumen.lecturas_entrada),
    ));
    salida.agregar(Bloque::valor(
        "Bases conservadas",
        de_un_total(resumen.bases_salida, resumen.bases_entrada),
    ));
    salida.agregar(Bloque::valor(
        "Adaptador recortado",
        if resumen.adaptador.is_empty() {
            "ninguno".to_string()
        } else {
            texto(&resumen.adaptador)
        },
    ));
    let despues = if resumen.bases_salida == 0 {
        "—".to_string()
    } else {
        porcentaje(cociente(resumen.bases_q30_salida, resumen.bases_salida), 1)
    };
    salida.agregar(Bloque::valor(
        "Bases con calidad Q30 o más (antes → después)",
        format!(
            "{} → {despues}",
            porcentaje(
                cociente(resumen.bases_q30_entrada, resumen.bases_entrada),
                1
            )
        ),
    ));
    salida.agregar(Bloque::tabla(
        "Qué se quitó",
        &["Motivo", "Lecturas", "Bases"],
        vec![
            vec![
                "Adaptador (desde donde empieza hasta el final de la lectura)".to_string(),
                miles(resumen.lecturas_con_adaptador),
                miles(resumen.bases_de_adaptador),
            ],
            vec![
                format!(
                    "… de ellas, con solo {} a {} bases de adaptador al final",
                    calidad::SOLAPAMIENTO_MINIMO,
                    calidad::BASES_POR_ERROR - 1
                ),
                miles(resumen.lecturas_con_adaptador_corto),
                "—".to_string(),
            ],
            vec![
                format!(
                    "Extremo 3' con calidad media menor que Q{} en ventanas de {} bases",
                    opciones.calidad_minima, opciones.ventana
                ),
                miles(resumen.lecturas_recortadas_por_calidad),
                miles(resumen.bases_por_calidad),
            ],
            vec![
                format!(
                    "Lecturas descartadas por tener más de {} bases N",
                    opciones.maximo_n
                ),
                miles(resumen.descartadas_por_n),
                "—".to_string(),
            ],
            vec![
                format!(
                    "Lecturas descartadas por quedar con menos de {} pb",
                    opciones.longitud_minima
                ),
                miles(resumen.descartadas_por_longitud),
                "—".to_string(),
            ],
        ],
    ));
}

fn cociente(parte: usize, total: usize) -> f64 {
    if total == 0 {
        0.0
    } else {
        parte as f64 / total as f64
    }
}

fn explicar_recorte(
    salida: &mut Salida,
    resumen: &ResumenRecorte,
    informe: &InformeCalidad,
    opciones: &OpcionesRecorte,
) {
    let lecturas = cociente(resumen.lecturas_salida, resumen.lecturas_entrada);
    let bases = cociente(resumen.bases_salida, resumen.bases_entrada);
    if lecturas < 0.7 || bases < 0.7 {
        salida.explicar(format!(
            "Solo se conservó el {} de las lecturas y el {} de las bases (menos de 70 % en alguno de los dos): los parámetros pueden ser demasiado estrictos o la corrida es mala. Revisa primero «Calidad de lecturas».",
            porcentaje(lecturas, 1),
            porcentaje(bases, 1)
        ));
    } else {
        salida.explicar(format!(
            "Se conservó el {} de las lecturas y el {} de las bases: la mayoría de los datos sigue disponible para el análisis.",
            porcentaje(lecturas, 1),
            porcentaje(bases, 1)
        ));
    }
    let con_adaptador = cociente(resumen.lecturas_con_adaptador, resumen.lecturas_entrada);
    match &opciones.adaptador {
        ModoAdaptador::Detectar if resumen.adaptador.is_empty() => salida.explicar(
            "No se detectó ningún adaptador en los datos, así que solo se recortó por calidad.",
        ),
        ModoAdaptador::Detectar => salida.explicar(format!(
            "El adaptador {} se descubrió en los propios datos (semillas de {} bases que abundan cerca del extremo 3') y se quitó de {} lecturas ({}).",
            texto(&resumen.adaptador),
            calidad::LONGITUD_SEMILLA,
            miles(resumen.lecturas_con_adaptador),
            porcentaje(con_adaptador, 1)
        )),
        ModoAdaptador::Conocido(_) => salida.explicar(format!(
            "Se buscó el adaptador que escribiste (permitiendo 1 diferencia cada {} bases) y se quitó de {} lecturas ({}).",
            calidad::BASES_POR_ERROR,
            miles(resumen.lecturas_con_adaptador),
            porcentaje(con_adaptador, 1)
        )),
        ModoAdaptador::Ninguno => salida.explicar("No se buscaron adaptadores porque elegiste «no»."),
    }
    if resumen.lecturas_con_adaptador_corto > 0 {
        salida.explicar(format!(
            "De esas, {} lecturas tenían solo {} a {} bases de adaptador al final. Coincidencias tan cortas a veces son azar (se esperan en cerca del 2 % de las lecturas), pero cortarlas cuesta muy poco.",
            miles(resumen.lecturas_con_adaptador_corto),
            calidad::SOLAPAMIENTO_MINIMO,
            calidad::BASES_POR_ERROR - 1
        ));
    }
    if resumen.bases_por_calidad > 0 {
        salida.explicar(format!(
            "Por calidad se quitaron {} bases de {} lecturas: sus extremos 3' tenían ventanas de {} bases con media menor que Q{}.",
            miles(resumen.bases_por_calidad),
            miles(resumen.lecturas_recortadas_por_calidad),
            opciones.ventana,
            opciones.calidad_minima
        ));
    } else {
        salida.explicar(format!(
            "Ninguna lectura tenía su extremo 3' por debajo de Q{} en ventanas de {} bases, así que no hubo recorte por calidad.",
            opciones.calidad_minima, opciones.ventana
        ));
    }
    let descartadas = resumen.descartadas_por_n + resumen.descartadas_por_longitud;
    if descartadas > 0 {
        let motivo = if resumen.descartadas_por_longitud >= resumen.descartadas_por_n {
            format!(
                "sobre todo por quedar con menos de {} pb; si son muchas, baja la longitud mínima",
                opciones.longitud_minima
            )
        } else {
            format!("sobre todo por tener más de {} bases N", opciones.maximo_n)
        };
        salida.explicar(format!(
            "Se descartaron {} lecturas, {motivo}.",
            miles(descartadas)
        ));
    }
    if let Some(adaptador) = &informe.adaptador
        && adaptador.fraccion > 0.01
    {
        salida.explicar(format!(
            "Aún queda un adaptador ({}) en el {} de las lecturas recortadas: escríbelo en el campo «adaptador» y vuelve a recortar.",
            texto(&adaptador.secuencia),
            porcentaje(adaptador.fraccion, 1)
        ));
    } else if resumen.bases_salida > 0 {
        salida.explicar(format!(
            "Las bases con Q30 o más pasaron de {} a {}. Usa el archivo recortado en «Mapeo de lecturas».",
            porcentaje(
                cociente(resumen.bases_q30_entrada, resumen.bases_entrada),
                1
            ),
            porcentaje(cociente(resumen.bases_q30_salida, resumen.bases_salida), 1)
        ));
    }
}

const LECTURAS_A_REVISAR_PARA_EJEMPLOS: usize = 20_000;
const BASES_EN_DETALLE: usize = 3;

fn corto(secuencia: &[u8]) -> String {
    let limite = 28;
    let mut resultado: String = secuencia
        .iter()
        .take(limite)
        .map(|&base| base as char)
        .collect();
    if secuencia.len() > limite {
        resultado.push('…');
    }
    resultado
}

fn nombre_del_archivo(ruta: &Path) -> String {
    ruta.file_name()
        .map(|nombre| nombre.to_string_lossy().to_string())
        .unwrap_or_else(|| ruta.display().to_string())
}

fn valor_o_guion<T: ToString>(valor: Option<T>) -> String {
    valor
        .map(|dato| dato.to_string())
        .unwrap_or_else(|| "–".to_string())
}

fn leer_primeras(ruta: &Path, cantidad: usize) -> Result<Vec<Lectura>, String> {
    let mut lecturas = Vec::new();
    for lectura in abrir_fastq(ruta)?.take(cantidad) {
        lecturas.push(lectura?);
    }
    Ok(lecturas)
}

fn letra_de(base: u8) -> usize {
    match base.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

#[derive(Default)]
struct Contadores {
    lecturas: usize,
    bases: usize,
    posicion: Option<usize>,
    base: Option<u8>,
    simbolo: Option<u8>,
    calidad: Option<usize>,
    letras: [usize; 5],
    bases_q20: usize,
    bases_q30: usize,
    gc: Option<usize>,
    bases_n: usize,
    bases_gc: usize,
    porcentaje_gc: Option<usize>,
}

impl Contadores {
    fn pizarra(&self) -> Pizarra {
        Pizarra::nueva()
            .variable("lecturas", self.lecturas)
            .variable("bases", self.bases)
            .variable("posicion", valor_o_guion(self.posicion))
            .variable("base", valor_o_guion(self.base.map(|base| base as char)))
            .variable(
                "simbolo",
                valor_o_guion(self.simbolo.map(|simbolo| simbolo as char)),
            )
            .variable("calidad", valor_o_guion(self.calidad))
            .variable(
                "letras",
                format!(
                    "A:{} C:{} G:{} T:{} N:{}",
                    self.letras[0], self.letras[1], self.letras[2], self.letras[3], self.letras[4]
                ),
            )
            .variable("bases_q20", self.bases_q20)
            .variable("bases_q30", self.bases_q30)
            .variable("gc", valor_o_guion(self.gc))
            .variable("bases_n", self.bases_n)
            .variable("bases_gc", self.bases_gc)
            .variable("porcentaje_gc", valor_o_guion(self.porcentaje_gc))
    }
}

fn contar_con_detalle(
    contadores: &mut Contadores,
    lectura: &Lectura,
    narracion: Option<&mut Narracion>,
) {
    let mut narracion = narracion;
    let largo = lectura.secuencia.len();
    contadores.lecturas += 1;
    contadores.bases += largo;
    contadores.letras = [0; 5];
    contadores.posicion = None;
    contadores.base = None;
    contadores.simbolo = None;
    contadores.calidad = None;
    contadores.gc = None;
    contadores.porcentaje_gc = None;
    if let Some(narracion) = narracion.as_deref_mut() {
        narracion.agregar(
            contadores
                .pizarra()
                .secuencia(&lectura.secuencia, &[(0, largo, Resalte::Actual)])
                .paso(
                    "self.bases += longitud;",
                    format!(
                        "La primera lectura mide {largo} bases. Sumamos esas bases al total y contamos una lectura más. Ahora vamos letra por letra."
                    ),
                ),
        );
    }
    let mut filas: Vec<Vec<String>> = Vec::new();
    for (posicion, (&base, &simbolo)) in lectura.secuencia.iter().zip(&lectura.calidad).enumerate()
    {
        let calidad = calidad::phred(simbolo);
        let letra = letra_de(base);
        contadores.posicion = Some(posicion);
        contadores.base = Some(base);
        contadores.simbolo = Some(simbolo);
        contadores.calidad = Some(calidad);
        if posicion < BASES_EN_DETALLE {
            filas.push(vec![
                (posicion + 1).to_string(),
                (base as char).to_string(),
                (simbolo as char).to_string(),
                calidad.to_string(),
            ]);
        }
        if posicion < BASES_EN_DETALLE
            && let Some(narracion) = narracion.as_deref_mut()
        {
            narracion.agregar(
                contadores
                    .pizarra()
                    .secuencia(
                        &lectura.secuencia,
                        &[(posicion, posicion + 1, Resalte::Actual)],
                    )
                    .tabla(
                        "Esta lectura, base por base",
                        &["base n.º", "letra", "símbolo", "Phred"],
                        filas.clone(),
                        &[(posicion, 3)],
                    )
                    .paso(
                        "let calidad = phred(simbolo);",
                        format!(
                            "El símbolo «{}» del FASTQ es un código: su número ASCII menos 33 es la calidad Phred, aquí {calidad}. Un Phred de {calidad} significa una probabilidad de error de 1 en {}.",
                            simbolo as char,
                            miles(10f64.powf(calidad as f64 / 10.0).round() as usize)
                        ),
                    ),
            );
        }
        contadores.letras[letra] += 1;
        contadores.bases_q20 += usize::from(calidad >= 20);
        contadores.bases_q30 += usize::from(calidad >= 30);
        if posicion < BASES_EN_DETALLE
            && let Some(narracion) = narracion.as_deref_mut()
        {
            narracion.agregar(
                contadores
                    .pizarra()
                    .secuencia(
                        &lectura.secuencia,
                        &[(posicion, posicion + 1, Resalte::Actual)],
                    )
                    .tabla(
                        "Esta lectura, base por base",
                        &["base n.º", "letra", "símbolo", "Phred"],
                        filas.clone(),
                        &[(posicion, 1)],
                    )
                    .paso(
                        "self.bases_q20 += usize::from(calidad >= 20);",
                        format!(
                            "Apuntamos la letra {} en su casilla y miramos la calidad: {calidad} {} Q20 y {} Q30, así que los contadores quedan en {} y {}.",
                            base as char,
                            if calidad >= 20 { "alcanza" } else { "no alcanza" },
                            if calidad >= 30 { "sí alcanza" } else { "no alcanza" },
                            contadores.bases_q20,
                            contadores.bases_q30
                        ),
                    ),
            );
        }
    }
    contadores.posicion = None;
    contadores.base = None;
    contadores.simbolo = None;
    contadores.calidad = None;
    let gc = contadores.letras[1] + contadores.letras[2];
    let acgt = contadores.letras[0] + gc + contadores.letras[3];
    contadores.gc = Some(gc);
    if let Some(narracion) = narracion.as_deref_mut() {
        narracion.agregar(
            contadores
                .pizarra()
                .secuencia(&lectura.secuencia, &[(0, largo, Resalte::Actual)])
                .paso(
                    "let gc = letras[1] + letras[2];",
                    format!(
                        "Las demás {} bases se cuentan igual. Al terminar la lectura sabemos cuántas hay de cada letra; las C más las G suman {gc} (la fracción GC).",
                        largo.saturating_sub(BASES_EN_DETALLE)
                    ),
                ),
        );
    }
    contadores.bases_n += contadores.letras[4];
    contadores.bases_gc += gc;
    if let Some(narracion) = narracion.as_deref_mut() {
        narracion.agregar(
            contadores.pizarra().paso(
                "self.bases_gc += gc;",
                format!(
                    "Sumamos las {gc} C y G de esta lectura al total del archivo y las {} bases N (letras que la máquina no pudo leer) a su propio total.",
                    contadores.letras[4]
                ),
            ),
        );
    }
    contadores.porcentaje_gc = (100 * gc + acgt / 2).checked_div(acgt);
    if let (Some(narracion), Some(porcentaje_gc)) = (narracion, contadores.porcentaje_gc) {
        narracion.agregar(
            contadores.pizarra().paso(
                "self.lecturas_por_gc[porcentaje_gc] += 1;",
                format!(
                    "Esta lectura tiene {porcentaje_gc} % de GC. Anotamos una lectura más en la casilla de {porcentaje_gc} %: así se dibuja después el histograma de GC."
                ),
            ),
        );
    }
}

fn pasos_calidad_de_lecturas(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let ruta = entrada.ruta("archivo")?;
    let informe = calidad::analizar_fastq(&ruta)?;
    if informe.lecturas == 0 {
        return Err("El archivo no tiene ninguna lectura FASTQ.".to_string());
    }
    let primeras = leer_primeras(&ruta, 3)?;
    let primera = primeras
        .first()
        .ok_or_else(|| "El archivo no tiene ninguna lectura FASTQ.".to_string())?;
    let nombre = nombre_del_archivo(&ruta);
    let mut narracion = Narracion::nueva(400);
    narracion.agregar(
        Pizarra::nueva()
            .variable("ruta", &nombre)
            .secuencia(
                &primera.secuencia,
                &[(0, primera.secuencia.len(), Resalte::Actual)],
            )
            .paso(
                "let mut lector = abrir_fastq(ruta)?;",
                format!(
                    "Vamos a revisar «{nombre}», que tiene {} lecturas. En un FASTQ cada lectura ocupa 4 líneas: su nombre, sus letras (arriba, la primera), una línea «+» y un símbolo de calidad por cada letra. Abrimos el archivo para leerlo de a poco.",
                    miles(informe.lecturas)
                ),
            ),
    );
    let mut contadores = Contadores::default();
    narracion.agregar(contadores.pizarra().paso(
        "let mut conteos = Conteos::nuevo();",
        "Preparamos las libretas de cuentas, todas en cero: lecturas, bases, cuántas bases pasan Q20 y Q30, cuántas son N y cuántas son G o C.",
    ));
    let bloques = informe
        .lecturas
        .div_ceil(8192)
        .max(informe.bases.div_ceil(4_000_000));
    narracion.agregar(contadores.pizarra().paso(
        "let mut tanda = leer_tanda(&mut lector)?;",
        format!(
            "El archivo se lee en bloques (hasta 8 192 lecturas o 4 millones de bases cada uno): aquí son al menos {}. Mientras un bloque se cuenta, el siguiente ya se está leyendo, y varios hilos del procesador cuentan a la vez.",
            miles(bloques)
        ),
    ));
    let mut resumen_de_lecturas: Vec<Vec<String>> = Vec::new();
    for (indice, lectura) in primeras.iter().enumerate() {
        if indice == 0 {
            contar_con_detalle(&mut contadores, lectura, Some(&mut narracion));
        } else {
            contar_con_detalle(&mut contadores, lectura, None);
        }
        let gc = contadores.porcentaje_gc.unwrap_or(0);
        resumen_de_lecturas.push(vec![
            (indice + 1).to_string(),
            lectura.secuencia.len().to_string(),
            format!("{gc} %"),
        ]);
        if indice > 0 {
            narracion.agregar(
                contadores
                    .pizarra()
                    .secuencia(
                        &lectura.secuencia,
                        &[(0, lectura.secuencia.len(), Resalte::Actual)],
                    )
                    .tabla(
                        "Lecturas contadas hasta ahora",
                        &["lectura", "bases", "GC"],
                        resumen_de_lecturas.clone(),
                        &[(indice, 2)],
                    )
                    .paso(
                        "conteos.contar_lectura(lectura);",
                        format!(
                            "La lectura {} se cuenta igual que la primera. Los totales ya van en {} lecturas y {} bases.",
                            indice + 1,
                            contadores.lecturas,
                            miles(contadores.bases)
                        ),
                    ),
            );
        }
    }
    let totales = Contadores {
        lecturas: informe.lecturas,
        bases: informe.bases,
        bases_q20: informe.bases_q20,
        bases_q30: informe.bases_q30,
        bases_n: informe.bases_n,
        bases_gc: informe.bases_gc,
        ..Contadores::default()
    };
    narracion.agregar(totales.pizarra().paso(
        "conteos.sumar(parcial);",
        format!(
            "Esto se repite para las {} lecturas. Cada hilo devuelve sus cuentas parciales y aquí se suman todas: {} bases, de las cuales {} tienen Q30 o más.",
            miles(informe.lecturas),
            miles(informe.bases),
            porcentaje(informe.fraccion_q30(), 1)
        ),
    ));
    narracion.agregar(totales.pizarra().variable("muestra", informe.lecturas_muestreadas).paso(
        "guardar_en_muestra(&mut muestra, tanda);",
        format!(
            "Además guardamos las letras de las primeras {} lecturas (el máximo es 200 000): con esa muestra se buscan las duplicadas y el adaptador, que necesitan comparar lecturas entre sí.",
            miles(informe.lecturas_muestreadas)
        ),
    ));
    narracion.agregar(totales.pizarra().variable("muestra", informe.lecturas_muestreadas).paso(
        "rayon::join(|| detectar_adaptador(muestra), || repeticiones(muestra));",
        "Con las cuentas y la muestra se arma el informe. Buscar el adaptador y contar repeticiones son trabajos independientes, así que se hacen a la vez en dos hilos.",
    ));
    let posiciones = &informe.calidad_por_posicion;
    let mut filas: Vec<Vec<String>> = Vec::new();
    let mut elegidas: Vec<usize> = (0..posiciones.len().min(4)).collect();
    if posiciones.len() > 4 {
        elegidas.push(posiciones.len() - 1);
    }
    for &posicion in &elegidas {
        let dato = &posiciones[posicion];
        filas.push(vec![
            (posicion + 1).to_string(),
            decimal(dato.media, 1),
            dato.mediana.to_string(),
            format!("{}–{}", dato.cuartil_inferior, dato.cuartil_superior),
        ]);
    }
    let ultima = filas.len().saturating_sub(1);
    narracion.agregar(
        totales
            .pizarra()
            .variable("calidad_por_posicion", format!("{} posiciones", posiciones.len()))
            .tabla(
                "Calidad por posición en la lectura",
                &["posición", "media", "mediana", "cuartiles"],
                filas,
                &[(ultima, 1)],
            )
            .paso(
                ".map(|(histograma, &suma)| resumir_posicion(histograma, suma))",
                "Para cada posición de la lectura (la 1.ª letra, la 2.ª…) se juntan las calidades de todas las lecturas y se resumen: la media, la mediana y los cuartiles. Así se ve si la calidad cae hacia el final.",
            ),
    );
    let unicas = informe.lecturas_muestreadas.saturating_sub(
        (informe.duplicacion * informe.lecturas_muestreadas as f64).round() as usize,
    );
    narracion.agregar(
        totales
            .pizarra()
            .variable("veces.len()", unicas)
            .paso(
                "let veces = contar_secuencias(muestra);",
                format!(
                    "Contamos cuántas veces aparece cada secuencia distinta en la muestra de {} lecturas. Hay {} distintas.",
                    miles(informe.lecturas_muestreadas),
                    miles(unicas)
                ),
            ),
    );
    narracion.agregar(
        totales
            .pizarra()
            .variable("total", informe.lecturas_muestreadas)
            .variable("veces.len()", unicas)
            .paso(
                "cociente(total.saturating_sub(veces.len()), total)",
                format!(
                    "Duplicadas = lecturas totales menos secuencias distintas, entre lecturas totales: {} de {} sobran, o sea {}.",
                    miles(informe.lecturas_muestreadas.saturating_sub(unicas)),
                    miles(informe.lecturas_muestreadas),
                    porcentaje(informe.duplicacion, 1)
                ),
            ),
    );
    let frecuentes: Vec<String> = informe
        .sobrerrepresentadas
        .iter()
        .take(5)
        .map(|repetida| {
            format!(
                "{} — {} veces ({})",
                corto(&repetida.secuencia),
                miles(repetida.veces),
                porcentaje(repetida.fraccion, 2)
            )
        })
        .collect();
    let aviso = if frecuentes.is_empty() {
        "Ninguna secuencia pasa el filtro: para destacar, una secuencia debe repetirse al menos 2 veces y ser más del 0,1 % de las lecturas.".to_string()
    } else {
        format!(
            "Solo cuentan las secuencias que se repiten al menos 2 veces y son más del 0,1 % de las lecturas. Quedaron {}; las más frecuentes están en la lista. Una secuencia muy repetida suele ser un adaptador u otro contaminante.",
            informe.sobrerrepresentadas.len()
        )
    };
    narracion.agregar(
        totales
            .pizarra()
            .variable("total", informe.lecturas_muestreadas)
            .lista("Secuencias sobrerrepresentadas", frecuentes, &[0])
            .paso(
                "if secuencia.is_empty() || cantidad < 2 || cantidad * 1000 <= total {",
                aviso,
            ),
    );
    match &informe.adaptador {
        Some(adaptador) => {
            narracion.agregar(totales.pizarra().variable("muestra", informe.lecturas_muestreadas).paso(
                "let semilla = elegir_semilla(muestra)?;",
                "Para buscar un adaptador se cuentan todos los trozos de 10 letras al inicio y al final de las lecturas. Un adaptador aparece mucho más al final (donde la lectura se pasó del ADN real y entró en él) que al principio; el trozo más enriquecido al final es la semilla.",
            ));
            narracion.agregar(
                totales
                    .pizarra()
                    .variable("secuencia", corto(&adaptador.secuencia))
                    .paso(
                        "let secuencia = extender_semilla(muestra, &semilla);",
                        format!(
                            "La semilla se alarga base por base hacia los dos lados mientras al menos el 70 % de las lecturas que la contienen estén de acuerdo en la siguiente letra. Resultó {} ({} letras).",
                            corto(&adaptador.secuencia),
                            adaptador.secuencia.len()
                        ),
                    ),
            );
            narracion.agregar(
                totales
                    .pizarra()
                    .variable("secuencia", corto(&adaptador.secuencia))
                    .paso(
                        ".filter(|lectura| buscar_adaptador(lectura, &secuencia, LONGITUD_SEMILLA).is_some())",
                        format!(
                            "Se revisa cuántas de las {} lecturas revisadas contienen ese adaptador (permitiendo un error cada 8 bases): {}, o sea {}.",
                            miles(adaptador.revisadas),
                            miles(adaptador.lecturas),
                            porcentaje(adaptador.fraccion, 1)
                        ),
                    ),
            );
        }
        None => {
            narracion.agregar(totales.pizarra().variable("muestra", informe.lecturas_muestreadas).paso(
                "let semilla = elegir_semilla(muestra)?;",
                "Se cuentan los trozos de 10 letras al inicio y al final de las lecturas buscando uno que abunde mucho más al final. Ninguno lo hace con la fuerza suficiente: no se detecta ningún adaptador.",
            ));
        }
    }
    let longitudes = if informe.longitud_minima == informe.longitud_maxima {
        format!("{} pb", miles(informe.longitud_maxima))
    } else {
        format!(
            "{}–{} pb",
            miles(informe.longitud_minima),
            miles(informe.longitud_maxima)
        )
    };
    let adaptador_final = informe
        .adaptador
        .as_ref()
        .map(|adaptador| {
            format!(
                "{} ({})",
                corto(&adaptador.secuencia),
                porcentaje(adaptador.fraccion, 1)
            )
        })
        .unwrap_or_else(|| "ninguno".to_string());
    Ok(narracion.terminar(
        totales
            .pizarra()
            .tabla(
                "Informe terminado",
                &["dato", "valor"],
                vec![
                    vec!["Lecturas".into(), miles(informe.lecturas)],
                    vec!["Bases".into(), miles(informe.bases)],
                    vec!["Longitud".into(), longitudes],
                    vec!["Q30 o más".into(), porcentaje(informe.fraccion_q30(), 1)],
                    vec!["GC medio".into(), porcentaje(informe.gc_medio(), 1)],
                    vec!["Duplicadas".into(), porcentaje(informe.duplicacion, 1)],
                    vec!["Adaptador".into(), adaptador_final],
                ],
                &[(2, 1)],
            )
            .paso(
                "Ok(armar_informe(conteos, &muestra))",
                "Listo: con todas las cuentas se arma el informe. Mira «Qué significa» para saber si estas cifras son buenas o malas para tu archivo.",
            ),
    ))
}

struct Vista {
    lectura: String,
    inicio_adaptador: String,
    sin_adaptador: String,
    inicio: String,
    fin: String,
    suma: String,
    umbral: String,
    descarte: String,
}

impl Vista {
    fn nueva(lectura: String) -> Vista {
        let guion = || "–".to_string();
        Vista {
            lectura,
            inicio_adaptador: guion(),
            sin_adaptador: guion(),
            inicio: guion(),
            fin: guion(),
            suma: guion(),
            umbral: guion(),
            descarte: guion(),
        }
    }

    fn pizarra(&self, opciones: &OpcionesRecorte) -> Pizarra {
        Pizarra::nueva()
            .variable("ventana", opciones.ventana)
            .variable("calidad_minima", opciones.calidad_minima)
            .variable("maximo_n", opciones.maximo_n)
            .variable("longitud_minima", opciones.longitud_minima)
            .variable("lectura", &self.lectura)
            .variable("inicio_adaptador", &self.inicio_adaptador)
            .variable("sin_adaptador", &self.sin_adaptador)
            .variable("inicio", &self.inicio)
            .variable("fin", &self.fin)
            .variable("suma", &self.suma)
            .variable("umbral", &self.umbral)
            .variable("descarte", &self.descarte)
    }
}

struct RecorridoDeCorte {
    estados: Vec<(usize, usize, usize)>,
    ultimo: (usize, usize, usize),
    vueltas: usize,
}

fn recorrer_corte(calidad: &[u8], ventana: usize, umbral: usize) -> RecorridoDeCorte {
    let ventana = ventana.max(1);
    let mut fin = calidad.len();
    let mut inicio = fin.saturating_sub(ventana);
    let mut suma: usize = calidad[inicio..fin]
        .iter()
        .map(|&simbolo| calidad::phred(simbolo))
        .sum();
    let mut recorrido = RecorridoDeCorte {
        estados: vec![(fin, inicio, suma)],
        ultimo: (fin, inicio, suma),
        vueltas: 0,
    };
    while fin > 0 && suma < umbral.saturating_mul(fin - inicio) {
        fin -= 1;
        suma -= calidad::phred(calidad[fin]);
        if inicio > 0 {
            inicio -= 1;
            suma += calidad::phred(calidad[inicio]);
        }
        recorrido.vueltas += 1;
        recorrido.ultimo = (fin, inicio, suma);
        if recorrido.estados.len() < 6 {
            recorrido.estados.push((fin, inicio, suma));
        }
    }
    recorrido
}

fn bases_n(secuencia: &[u8]) -> usize {
    secuencia
        .iter()
        .filter(|&&base| letra_de(base) == 4)
        .count()
}

fn tabla_de_calidades(calidad: &[u8], hasta: usize) -> Vec<Vec<String>> {
    let desde = hasta.saturating_sub(8);
    (desde..hasta)
        .map(|posicion| {
            let simbolo = calidad[posicion];
            vec![
                (posicion + 1).to_string(),
                (simbolo as char).to_string(),
                calidad::phred(simbolo).to_string(),
            ]
        })
        .collect()
}

fn narrar_lectura(
    narracion: &mut Narracion,
    lectura: &Lectura,
    adaptador: &[u8],
    opciones: &OpcionesRecorte,
    posicion_en_archivo: usize,
    motivo: &str,
) {
    let recorte = calidad::recortar_lectura(lectura, adaptador, opciones);
    let largo = lectura.secuencia.len();
    let nombre = corto(&lectura.nombre);
    let mut vista = Vista::nueva(nombre.clone());
    narracion.agregar(
        vista
            .pizarra(opciones)
            .secuencia(&lectura.secuencia, &[(0, largo, Resalte::Actual)])
            .paso(
                "let recorte = recortar_lectura(lectura, adaptador, opciones);",
                format!(
                    "Tomamos la lectura n.º {} ({motivo}). Mide {largo} bases. recortar_lectura decide dónde cortarla o si hay que descartarla.",
                    miles(posicion_en_archivo + 1)
                ),
            ),
    );
    vista.inicio_adaptador = match recorte.inicio_adaptador {
        Some(inicio) => (inicio + 1).to_string(),
        None => "ninguno".to_string(),
    };
    let marcas = match recorte.inicio_adaptador {
        Some(inicio) => vec![
            (0, inicio, Resalte::Actual),
            (inicio, largo, Resalte::Hallado),
        ],
        None => vec![(0, largo, Resalte::Actual)],
    };
    let nota = if adaptador.is_empty() {
        "No hay adaptador que buscar en esta corrida, así que no se corta nada por este lado."
            .to_string()
    } else if let Some(inicio) = recorte.inicio_adaptador {
        format!(
            "Se compara el adaptador contra cada posición de la lectura, aceptando 1 error por cada 8 bases. Encaja desde la base {}: de ahí al final ({} bases, marcadas) es adaptador y sobra.",
            inicio + 1,
            largo - inicio
        )
    } else {
        "Se compara el adaptador contra cada posición de la lectura, aceptando 1 error por cada 8 bases. No encaja en ninguna parte: esta lectura no lo trae.".to_string()
    };
    narracion.agregar(
        vista
            .pizarra(opciones)
            .secuencia(&lectura.secuencia, &marcas)
            .paso(
                "let inicio_adaptador = buscar_adaptador(&lectura.secuencia, adaptador, SOLAPAMIENTO_MINIMO);",
                nota,
            ),
    );
    let sin_adaptador = recorte
        .inicio_adaptador
        .unwrap_or(largo)
        .min(lectura.calidad.len());
    vista.sin_adaptador = sin_adaptador.to_string();
    narracion.agregar(
        vista
            .pizarra(opciones)
            .secuencia(&lectura.secuencia, &[(0, sin_adaptador, Resalte::Actual)])
            .paso(
                "let sin_adaptador = inicio_adaptador",
                format!(
                    "Nos quedamos solo con las primeras {sin_adaptador} bases: todo lo que viene después del adaptador se ignora antes de mirar la calidad."
                ),
            ),
    );
    let calidades = &lectura.calidad[..sin_adaptador];
    let recorrido = recorrer_corte(calidades, opciones.ventana, opciones.calidad_minima);
    let (fin0, inicio0, suma0) = recorrido.estados[0];
    let umbral0 = opciones.calidad_minima.saturating_mul(fin0 - inicio0);
    vista.inicio = (inicio0 + 1).to_string();
    vista.fin = fin0.to_string();
    vista.suma = suma0.to_string();
    vista.umbral = umbral0.to_string();
    narracion.agregar(
        vista
            .pizarra(opciones)
            .secuencia(&lectura.secuencia, &[(inicio0, fin0, Resalte::Actual)])
            .tabla(
                "Calidades del final de la lectura",
                &["base n.º", "símbolo", "Phred"],
                tabla_de_calidades(calidades, fin0),
                &[],
            )
            .paso(
                "let mut suma: usize = calidad[inicio..fin]",
                format!(
                    "Se mira una ventana de {} bases pegada al final. Sus calidades suman {suma0}; para pasar necesitan sumar al menos {} × {} = {umbral0}.",
                    fin0 - inicio0,
                    opciones.calidad_minima,
                    fin0 - inicio0
                ),
            ),
    );
    for (turno, par) in recorrido.estados.windows(2).enumerate() {
        let (fin_antes, inicio_antes, suma_antes) = par[0];
        let (fin_despues, inicio_despues, suma_despues) = par[1];
        let necesita = opciones
            .calidad_minima
            .saturating_mul(fin_antes - inicio_antes);
        vista.inicio = (inicio_despues + 1).to_string();
        vista.fin = fin_despues.to_string();
        vista.suma = suma_despues.to_string();
        vista.umbral = opciones
            .calidad_minima
            .saturating_mul(fin_despues - inicio_despues)
            .to_string();
        narracion.agregar(
            vista
                .pizarra(opciones)
                .secuencia(
                    &lectura.secuencia,
                    &[(inicio_despues, fin_despues, Resalte::Actual)],
                )
                .paso(
                    "fin -= 1;",
                    format!(
                        "Vuelta {}: la ventana suma {suma_antes} y necesita {necesita}, así que es mala. Se descarta la última base y la ventana se corre una base hacia el inicio: ahora suma {suma_despues}.",
                        turno + 1
                    ),
                ),
        );
    }
    let (fin, inicio, suma) = recorrido.ultimo;
    vista.inicio = (inicio + 1).to_string();
    vista.fin = fin.to_string();
    vista.suma = suma.to_string();
    vista.umbral = opciones
        .calidad_minima
        .saturating_mul(fin - inicio)
        .to_string();
    let nota = if recorrido.vueltas == 0 {
        "La ventana final ya cumple: la lectura no se corta por calidad.".to_string()
    } else if fin == 0 {
        format!(
            "Tras {} vueltas se gastó toda la lectura: ninguna ventana alcanza la calidad pedida.",
            recorrido.vueltas
        )
    } else {
        format!(
            "Tras {} {} la ventana suma {suma} y ya alcanza lo pedido: se corta aquí. La lectura se queda con {fin} bases.{}",
            recorrido.vueltas,
            if recorrido.vueltas == 1 {
                "vuelta"
            } else {
                "vueltas"
            },
            if recorrido.vueltas > 5 {
                " Las vueltas intermedias se hacen igual que las anteriores."
            } else {
                ""
            }
        )
    };
    narracion.agregar(
        vista
            .pizarra(opciones)
            .secuencia(
                &lectura.secuencia,
                &[(0, fin, Resalte::Hallado), (inicio, fin, Resalte::Actual)],
            )
            .paso(
                "while fin > 0 && suma < umbral.saturating_mul(fin - inicio) {",
                nota,
            ),
    );
    let fin = recorte.fin;
    vista.fin = fin.to_string();
    let enes = bases_n(&lectura.secuencia[..fin]);
    vista.descarte = match recorte.descarte {
        Some(Descarte::DemasiadasN) => "Some(DemasiadasN)".to_string(),
        Some(Descarte::MuyCorta) => "Some(MuyCorta)".to_string(),
        None => "None".to_string(),
    };
    narracion.agregar(
        vista
            .pizarra(opciones)
            .secuencia(&lectura.secuencia, &[(0, fin, Resalte::Hallado)])
            .paso(
                "let descarte = if contar_n(&lectura.secuencia[..fin]) > opciones.maximo_n {",
                if recorte.descarte == Some(Descarte::DemasiadasN) {
                    format!(
                        "En lo que queda hay {enes} bases N (no leídas) y el máximo es {}: demasiadas, la lectura se descarta.",
                        opciones.maximo_n
                    )
                } else {
                    format!(
                        "En lo que queda hay {enes} bases N y el máximo permitido es {}: está bien, la lectura sigue.",
                        opciones.maximo_n
                    )
                },
            ),
    );
    if recorte.descarte != Some(Descarte::DemasiadasN) {
        narracion.agregar(
            vista
                .pizarra(opciones)
                .secuencia(&lectura.secuencia, &[(0, fin, Resalte::Hallado)])
                .paso(
                    "} else if fin < opciones.longitud_minima {",
                    if recorte.descarte == Some(Descarte::MuyCorta) {
                        format!(
                            "Quedaron {fin} bases y se piden al menos {}: es demasiado corta para servir, se descarta.",
                            opciones.longitud_minima
                        )
                    } else {
                        format!(
                            "Quedaron {fin} bases y se piden al menos {}: alcanza, la lectura se conserva.",
                            opciones.longitud_minima
                        )
                    },
                ),
        );
    }
    narracion.agregar(
        vista
            .pizarra(opciones)
            .secuencia(&lectura.secuencia, &[(0, fin, Resalte::Hallado)])
            .paso(
                "resumen.anotar(lectura, &recorte);",
                if recorte.descarte.is_some() {
                    "Se anota en el resumen como descartada y no se escribe en el archivo limpio.".to_string()
                } else {
                    format!(
                        "Se anota en el resumen y se escribe en el archivo limpio con {fin} de sus {largo} bases ({} {}).",
                        largo - fin,
                        if largo - fin == 1 { "quitada" } else { "quitadas" }
                    )
                },
            ),
    );
}

fn elegir_lecturas_de_ejemplo(
    ruta: &Path,
    adaptador: &[u8],
    opciones: &OpcionesRecorte,
) -> Result<Vec<(usize, Lectura, &'static str)>, String> {
    let mut con_adaptador = None;
    let mut por_calidad = None;
    let mut descartada = None;
    let mut primera = None;
    for (indice, lectura) in abrir_fastq(ruta)?
        .take(LECTURAS_A_REVISAR_PARA_EJEMPLOS)
        .enumerate()
    {
        let lectura = lectura?;
        let recorte = calidad::recortar_lectura(&lectura, adaptador, opciones);
        let sin_adaptador = recorte
            .inicio_adaptador
            .unwrap_or(lectura.secuencia.len())
            .min(lectura.calidad.len());
        if primera.is_none() {
            primera = Some((indice, lectura.clone(), "la primera del archivo"));
        }
        if con_adaptador.is_none() && recorte.inicio_adaptador.is_some() {
            con_adaptador = Some((indice, lectura.clone(), "la primera que trae adaptador"));
        } else if por_calidad.is_none()
            && recorte.inicio_adaptador.is_none()
            && recorte.fin < sin_adaptador
        {
            por_calidad = Some((
                indice,
                lectura.clone(),
                "la primera con el final de mala calidad",
            ));
        } else if descartada.is_none() && recorte.descarte.is_some() {
            descartada = Some((indice, lectura, "la primera que se descarta"));
        }
        if con_adaptador.is_some() && por_calidad.is_some() && descartada.is_some() {
            break;
        }
    }
    let mut elegidas: Vec<(usize, Lectura, &'static str)> =
        [con_adaptador, por_calidad, descartada]
            .into_iter()
            .flatten()
            .collect();
    if elegidas.is_empty() {
        elegidas.extend(primera);
    }
    elegidas.sort_by_key(|(indice, _, _)| *indice);
    Ok(elegidas)
}

fn pasos_recortar_lecturas(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let ruta = entrada.ruta("archivo")?;
    let opciones = leer_opciones(entrada)?;
    let temporal = ruta_temporal(&crate::carpeta_de_ejemplos()?.join("narracion.fastq"));
    let resultado = calidad::recortar_fastq(&ruta, &temporal, &opciones);
    let _ = fs::remove_file(&temporal);
    let resumen = resultado?;
    if resumen.lecturas_entrada == 0 {
        return Err("El archivo no tiene ninguna lectura FASTQ.".to_string());
    }
    let nombre = nombre_del_archivo(&ruta);
    let adaptador = resumen.adaptador.clone();
    let ejemplos = elegir_lecturas_de_ejemplo(&ruta, &adaptador, &opciones)?;
    let mut narracion = Narracion::nueva(400);
    let vista = Vista::nueva("–".to_string());
    let modo = match &opciones.adaptador {
        ModoAdaptador::Detectar => "se detecta solo",
        ModoAdaptador::Ninguno => "no se busca",
        ModoAdaptador::Conocido(_) => "lo escribiste tú",
    };
    narracion.agregar(vista.pizarra(&opciones).paso(
        "revisar_opciones(entrada, salida, opciones)?;",
        format!(
            "Vamos a limpiar «{nombre}» ({} lecturas). Cada lectura pasa por cuatro filtros: quitar el adaptador ({modo}), cortar el final si su calidad media en una ventana de {} bases baja de {}, descartar si quedan más de {} bases N y descartar si queda más corta que {} bases.",
            miles(resumen.lecturas_entrada),
            opciones.ventana,
            opciones.calidad_minima,
            opciones.maximo_n,
            opciones.longitud_minima
        ),
    ));
    let nota_del_adaptador = match &opciones.adaptador {
        ModoAdaptador::Ninguno => "Elegiste no buscar adaptador: la lista queda vacía.".to_string(),
        ModoAdaptador::Conocido(_) => format!(
            "Usamos el adaptador que escribiste ({} letras).",
            adaptador.len()
        ),
        ModoAdaptador::Detectar if adaptador.is_empty() => {
            "Se leyeron hasta 100 000 lecturas buscando un trozo que abunde al final de las lecturas, y no apareció ninguno: no se quitará adaptador.".to_string()
        }
        ModoAdaptador::Detectar => format!(
            "Se leyeron hasta 100 000 lecturas buscando un trozo que abunde al final de las lecturas, y se extendió hasta formar el adaptador de {} letras.",
            adaptador.len()
        ),
    };
    narracion.agregar(
        vista
            .pizarra(&opciones)
            .variable("adaptador", corto(&adaptador))
            .paso(
                "let adaptador = elegir_adaptador(entrada, &opciones.adaptador)?;",
                nota_del_adaptador,
            ),
    );
    narracion.agregar(vista.pizarra(&opciones).variable("adaptador", corto(&adaptador)).paso(
        "let mut tanda = leer_tanda(&mut lector)?;",
        "Las lecturas se leen en bloques y cada bloque lo recorta un hilo distinto, todos a la vez. Las lecturas no dependen unas de otras, por eso se pueden repartir. Veamos con detalle algunas lecturas representativas.",
    ));
    for (indice, lectura, motivo) in &ejemplos {
        narrar_lectura(
            &mut narracion,
            lectura,
            &adaptador,
            &opciones,
            *indice,
            motivo,
        );
    }
    let tabla = vec![
        vec![
            "Lecturas que entraron".into(),
            miles(resumen.lecturas_entrada),
        ],
        vec![
            "Con adaptador".into(),
            de_un_total(resumen.lecturas_con_adaptador, resumen.lecturas_entrada),
        ],
        vec![
            "Recortadas por calidad".into(),
            de_un_total(
                resumen.lecturas_recortadas_por_calidad,
                resumen.lecturas_entrada,
            ),
        ],
        vec![
            "Descartadas por N".into(),
            de_un_total(resumen.descartadas_por_n, resumen.lecturas_entrada),
        ],
        vec![
            "Descartadas por cortas".into(),
            de_un_total(resumen.descartadas_por_longitud, resumen.lecturas_entrada),
        ],
        vec![
            "Lecturas que quedaron".into(),
            de_un_total(resumen.lecturas_salida, resumen.lecturas_entrada),
        ],
    ];
    narracion.agregar(
        vista
            .pizarra(&opciones)
            .variable("adaptador", corto(&adaptador))
            .variable("lecturas_entrada", resumen.lecturas_entrada)
            .variable("lecturas_salida", resumen.lecturas_salida)
            .tabla("Cuentas de todo el archivo", &["qué", "cuántas"], tabla.clone(), &[(5, 1)])
            .paso(
                "resumen.acumular(&recortado.resumen);",
                format!(
                    "Lo mismo se hizo con las otras {} lecturas, y cada hilo entregó su resumen parcial; aquí se suman todos.",
                    miles(resumen.lecturas_entrada.saturating_sub(ejemplos.len()))
                ),
            ),
    );
    Ok(narracion.terminar(
        vista
            .pizarra(&opciones)
            .variable("adaptador", corto(&adaptador))
            .variable("lecturas_entrada", resumen.lecturas_entrada)
            .variable("lecturas_salida", resumen.lecturas_salida)
            .tabla(
                "Bases antes y después",
                &["qué", "bases"],
                vec![
                    vec!["Bases de entrada".into(), miles(resumen.bases_entrada)],
                    vec!["Bases que quedaron".into(), miles(resumen.bases_salida)],
                    vec![
                        "Bases de adaptador quitadas".into(),
                        miles(resumen.bases_de_adaptador),
                    ],
                    vec![
                        "Bases quitadas por calidad".into(),
                        miles(resumen.bases_por_calidad),
                    ],
                ],
                &[(1, 1)],
            )
            .paso(
                "Ok(resumen)",
                format!(
                    "Terminamos: quedaron {} de {} lecturas. En la ventana de resultados, el archivo limpio se guarda junto al original; aquí solo seguimos los cálculos. «Qué significa» te dice si el recorte fue suave o agresivo.",
                    miles(resumen.lecturas_salida),
                    miles(resumen.lecturas_entrada)
                ),
            ),
    ))
}

#[cfg(test)]
mod pruebas {
    use super::*;

    fn entrada_de_ejemplo(id: &str) -> Entrada {
        let valores = crate::ejemplo(id).unwrap();
        Entrada::desde_json(valores).unwrap()
    }

    fn variable(paso: &Paso, nombre: &str) -> String {
        paso.estado["variables"]
            .as_array()
            .unwrap()
            .iter()
            .find(|par| par[0] == nombre)
            .map(|par| par[1].as_str().unwrap().to_string())
            .unwrap_or_default()
    }

    #[test]
    fn la_narracion_de_calidad_cuadra_con_el_analisis() {
        let entrada = entrada_de_ejemplo("puerto.calidad_lecturas");
        let pasos = pasos_calidad_de_lecturas(&entrada).unwrap();
        assert!(pasos.len() > 10 && pasos.len() <= 400);
        let informe = calidad::analizar_fastq(&entrada.ruta("archivo").unwrap()).unwrap();
        let sumar = pasos
            .iter()
            .find(|paso| paso.linea == "conteos.sumar(parcial);")
            .unwrap();
        assert_eq!(variable(sumar, "lecturas"), informe.lecturas.to_string());
        assert_eq!(variable(sumar, "bases_q30"), informe.bases_q30.to_string());
        let ultimo = pasos.last().unwrap();
        assert!(ultimo.nota.contains("Qué significa"));
        assert_eq!(
            ultimo.estado["tabla"]["filas"][0][1],
            miles(informe.lecturas)
        );
    }

    #[test]
    fn la_narracion_de_recorte_cuadra_con_el_recorte() {
        let entrada = entrada_de_ejemplo("puerto.recortar_lecturas");
        let pasos = pasos_recortar_lecturas(&entrada).unwrap();
        assert!(pasos.len() > 10 && pasos.len() <= 400);
        let ruta = entrada.ruta("archivo").unwrap();
        let opciones = leer_opciones(&entrada).unwrap();
        let temporal = ruta_temporal(&crate::carpeta_de_ejemplos().unwrap().join("prueba.fastq"));
        let resumen = calidad::recortar_fastq(&ruta, &temporal, &opciones).unwrap();
        let _ = fs::remove_file(&temporal);
        let ultimo = pasos.last().unwrap();
        assert_eq!(
            variable(ultimo, "lecturas_salida"),
            resumen.lecturas_salida.to_string()
        );
        assert_eq!(
            variable(ultimo, "lecturas_entrada"),
            resumen.lecturas_entrada.to_string()
        );
    }

    #[test]
    fn el_corte_narrado_coincide_con_el_del_nucleo() {
        let ruta = entrada_de_ejemplo("puerto.recortar_lecturas")
            .ruta("archivo")
            .unwrap();
        for lectura in abrir_fastq(&ruta).unwrap().take(3000) {
            let lectura = lectura.unwrap();
            for (ventana, umbral) in [(4, 20), (1, 30), (10, 25)] {
                let narrado = recorrer_corte(&lectura.calidad, ventana, umbral);
                assert_eq!(
                    narrado.ultimo.0,
                    calidad::corte_por_calidad(&lectura.calidad, ventana, umbral)
                );
            }
        }
    }

    #[test]
    fn el_recorte_narrado_usa_el_mismo_adaptador() {
        let entrada = entrada_de_ejemplo("puerto.recortar_lecturas");
        let pasos = pasos_recortar_lecturas(&entrada).unwrap();
        let paso = pasos
            .iter()
            .find(|paso| paso.linea.starts_with("let adaptador = elegir_adaptador"))
            .unwrap();
        assert!(!variable(paso, "adaptador").is_empty());
    }
}
