use crate::formato::{decimal, miles, porcentaje};
use crate::pizarra::{Marca, Narracion, Pizarra};
use crate::{Algoritmo, Bloque, Campo, Entrada, Paso, Salida, Tipo};
use nucleo::ecoli::{self, Gen};
use nucleo::genes::{
    self, Anotado, CODONES_DE_INICIO, Comparacion, GenPredicho, Hebras, LONGITUD_MINIMA,
    LONGITUD_PARA_ENTRENAR, Marco, Modelo, Prediccion, RONDAS_DE_ENTRENAMIENTO,
    SECUENCIA_MINIMA_PARA_ENTRENAR, buscar_rbs, candidatos_del_marco, codon_de_inicio,
    elegir_genes, es_codon_de_parada, hexameros_de_fondo, hexameros_en_genes, marcos_abiertos,
    marcos_mas_largos, preparar_hebras, razon_logaritmica,
};
use nucleo::secuencia::{complemento_inverso, contenido_gc};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;
use std::sync::atomic::{AtomicU64, Ordering};

const LONGITUD_MAXIMA: usize = 20_000_000;
const LONGITUD_MEDIA_ECOLI: f64 = 950.0;
const USO_DE_INICIOS_ECOLI: [f64; 3] = [0.83, 0.14, 0.03];
const MINIMO_DE_GENES_PARA_OPINAR: usize = 10;
const LONGITUD_PARA_OPINAR_DEL_GC: usize = 1000;
const LETRAS_DE_ADN: &[u8] = b"ACGTNRYKMSWBDHV";
const CLASES_DE_LONGITUD: &[(usize, usize)] = &[
    (90, 300),
    (300, 600),
    (600, 900),
    (900, 1200),
    (1200, 1500),
    (1500, 2000),
    (2000, 3000),
    (3000, usize::MAX),
];

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![Algoritmo {
        id: "puerto.predecir_genes",
        capitulo: "puerto",
        nombre: "Predicción de genes bacterianos",
        resumen: "Encuentra los genes que codifican proteínas en ADN bacteriano con un modelo que aprende de la propia secuencia, al estilo de Prodigal.",
        campos: &[Campo {
            id: "adn",
            etiqueta: "ADN bacteriano (genoma o región de al menos 90 pb)",
            tipo: Tipo::Secuencia,
            ejemplo: "@ecoli:1-150000",
        }],
        funcion: "genes::predecir_genes",
        ejecutar: predecir,
        pasos: Some(trazar_genes),
    }]
}

fn modelo_ecoli() -> &'static Modelo {
    static MODELO: OnceLock<Modelo> = OnceLock::new();
    MODELO.get_or_init(|| genes::entrenar_modelo(ecoli::genoma()))
}

pub fn preparar() {
    modelo_ecoli();
}

fn leer_adn(entrada: &Entrada) -> Result<Vec<u8>, String> {
    let adn = entrada.secuencia("adn")?;
    if let Some(letra) = adn.iter().find(|letra| !LETRAS_DE_ADN.contains(letra)) {
        return Err(format!(
            "La secuencia solo puede tener bases de ADN (A, C, G, T, y N u otro código IUPAC para las dudosas), pero contiene «{}».",
            *letra as char
        ));
    }
    if adn.len() < LONGITUD_MINIMA {
        return Err(format!(
            "La secuencia mide {} pb y un gen necesita al menos {LONGITUD_MINIMA} pb.",
            adn.len()
        ));
    }
    if adn.len() > LONGITUD_MAXIMA {
        return Err(format!(
            "La secuencia mide {} pb; el límite es {} pb, más que el genoma bacteriano más grande conocido.",
            miles(adn.len()),
            miles(LONGITUD_MAXIMA)
        ));
    }
    let conocidas = adn
        .iter()
        .filter(|letra| matches!(letra, b'A' | b'C' | b'G' | b'T'))
        .count();
    if conocidas * 2 < adn.len() {
        return Err(
            "Más de la mitad de la secuencia son bases desconocidas; así no se pueden reconocer genes."
                .to_string(),
        );
    }
    Ok(adn)
}

struct Ubicacion {
    desplazamiento: usize,
    longitud: usize,
    inversa: bool,
}

impl Ubicacion {
    fn en_ecoli(&self, desde: usize, hasta: usize) -> (usize, usize) {
        if self.inversa {
            let final_region = self.desplazamiento + self.longitud;
            (final_region - hasta + 1, final_region - desde)
        } else {
            (self.desplazamiento + desde + 1, self.desplazamiento + hasta)
        }
    }

    fn nombre(&self) -> String {
        let primera = self.desplazamiento + 1;
        let ultima = self.desplazamiento + self.longitud;
        let completo = self.longitud == ecoli::genoma().len();
        match (completo, self.inversa) {
            (true, false) => ecoli::IDENTIFICADOR.to_string(),
            (true, true) => format!("{}_complementaria", ecoli::IDENTIFICADOR),
            (false, false) => format!("{}_{primera}-{ultima}", ecoli::IDENTIFICADOR),
            (false, true) => format!("{}_{primera}-{ultima}_complementaria", ecoli::IDENTIFICADOR),
        }
    }
}

fn buscar(genoma: &[u8], adn: &[u8]) -> Option<usize> {
    if adn.len() > genoma.len() {
        return None;
    }
    genoma.windows(adn.len()).position(|ventana| ventana == adn)
}

fn ubicar_en_ecoli(adn: &[u8]) -> Option<Ubicacion> {
    let genoma = ecoli::genoma();
    if let Some(desplazamiento) = buscar(genoma, adn) {
        return Some(Ubicacion {
            desplazamiento,
            longitud: adn.len(),
            inversa: false,
        });
    }
    let inversa = complemento_inverso(adn);
    buscar(genoma, &inversa).map(|desplazamiento| Ubicacion {
        desplazamiento,
        longitud: adn.len(),
        inversa: true,
    })
}

fn anotados_en_region(ubicacion: &Ubicacion) -> Vec<(&'static Gen, Anotado)> {
    let primera = ubicacion.desplazamiento;
    let ultima = ubicacion.desplazamiento + ubicacion.longitud;
    let mut anotados = Vec::new();
    for anotado in ecoli::genes() {
        let desde = anotado.inicio - 1;
        if desde < primera || anotado.fin > ultima {
            continue;
        }
        let en_la_entrada = if ubicacion.inversa {
            Anotado {
                desde: ultima - anotado.fin,
                hasta: ultima - desde,
                hebra: if anotado.hebra == b'-' { b'+' } else { b'-' },
            }
        } else {
            Anotado {
                desde: desde - primera,
                hasta: anotado.fin - primera,
                hebra: anotado.hebra,
            }
        };
        anotados.push((anotado, en_la_entrada));
    }
    anotados
}

fn nombre_de_archivo(adn: &[u8], ubicacion: Option<&Ubicacion>) -> String {
    match ubicacion {
        Some(lugar) => format!("genes_{}.gff3", lugar.nombre().replace('.', "_")),
        None => {
            let mut huella: u64 = 0xcbf2_9ce4_8422_2325;
            for &base in adn {
                huella = (huella ^ base as u64).wrapping_mul(0x0100_0000_01b3);
            }
            format!("genes_secuencia_{}pb_{:08x}.gff3", adn.len(), huella as u32)
        }
    }
}

fn guardar_gff(
    ruta: &Path,
    nombre: &str,
    longitud: usize,
    predichos: &[GenPredicho],
) -> Result<(), String> {
    static CONTADOR: AtomicU64 = AtomicU64::new(0);
    let numero = CONTADOR.fetch_add(1, Ordering::Relaxed);
    let temporal = ruta.with_extension(format!("tmp{}-{numero}", std::process::id()));
    let archivo = File::create(&temporal)
        .map_err(|error| format!("No se pudo crear {}: {error}", temporal.display()))?;
    let mut destino = BufWriter::new(archivo);
    genes::escribir_gff(&mut destino, nombre, longitud, predichos)
        .and_then(|_| destino.flush())
        .map_err(|error| format!("No se pudo escribir {}: {error}", temporal.display()))?;
    drop(destino);
    fs::rename(&temporal, ruta)
        .map_err(|error| format!("No se pudo guardar {}: {error}", ruta.display()))
}

fn densidad_codificante(predichos: &[GenPredicho], longitud: usize) -> f64 {
    let mut cubiertas = 0;
    let mut hasta_donde = 0;
    for predicho in predichos {
        let desde = predicho.desde.max(hasta_donde);
        if predicho.hasta > desde {
            cubiertas += predicho.hasta - desde;
        }
        hasta_donde = hasta_donde.max(predicho.hasta);
    }
    cubiertas as f64 / longitud.max(1) as f64
}

fn uso_de_inicios(predichos: &[GenPredicho]) -> ([usize; 3], usize) {
    let mut cuentas = [0; 3];
    for predicho in predichos {
        if let Some(indice) = CODONES_DE_INICIO
            .iter()
            .position(|codon| *codon == predicho.codon_de_inicio)
        {
            cuentas[indice] += 1;
        }
    }
    let total = cuentas.iter().sum();
    (cuentas, total)
}

fn texto_parcial(predicho: &GenPredicho) -> &'static str {
    match (predicho.parcial_al_inicio, predicho.parcial_al_final) {
        (false, false) => "no",
        (true, false) => "le falta el inicio",
        (false, true) => "le falta la parada",
        (true, true) => "le faltan ambos",
    }
}

fn inicio_biologico(desde: usize, hasta: usize, hebra: u8) -> usize {
    if hebra == b'-' { hasta } else { desde }
}

fn fila_del_gen(
    numero: usize,
    predicho: &GenPredicho,
    referencia: Option<(&Ubicacion, Option<(&Gen, Anotado)>)>,
) -> Vec<String> {
    let codones = predicho.longitud() / 3;
    let aminoacidos = if predicho.parcial_al_final {
        codones
    } else {
        codones.saturating_sub(1)
    };
    let (motivo, espaciador) = match predicho.rbs {
        Some(rbs) => (rbs.motivo.to_string(), rbs.espaciador.to_string()),
        None => ("—".to_string(), "—".to_string()),
    };
    let codon = if predicho.parcial_al_inicio {
        "—".to_string()
    } else {
        predicho.codon_de_inicio.to_string()
    };
    let mut fila = vec![
        miles(numero),
        miles(predicho.desde + 1),
        miles(predicho.hasta),
        (predicho.hebra as char).to_string(),
        miles(predicho.longitud()),
        miles(aminoacidos),
        codon,
        motivo,
        espaciador,
        decimal(predicho.puntaje, 1),
        texto_parcial(predicho).to_string(),
    ];
    if let Some((ubicacion, pareja)) = referencia {
        let (desde, hasta) = ubicacion.en_ecoli(predicho.desde, predicho.hasta);
        fila.push(format!("{}–{}", miles(desde), miles(hasta)));
        match pareja {
            Some((anotado, en_la_entrada)) => {
                fila.push(format!("{} ({})", anotado.nombre, anotado.locus));
                let exacto =
                    en_la_entrada.desde == predicho.desde && en_la_entrada.hasta == predicho.hasta;
                if exacto {
                    fila.push("inicio exacto".to_string());
                } else {
                    let predicho_inicio =
                        inicio_biologico(predicho.desde, predicho.hasta, predicho.hebra);
                    let anotado_inicio = inicio_biologico(
                        en_la_entrada.desde,
                        en_la_entrada.hasta,
                        en_la_entrada.hebra,
                    );
                    fila.push(format!(
                        "misma parada; el inicio difiere en {} pb",
                        miles(predicho_inicio.abs_diff(anotado_inicio))
                    ));
                }
            }
            None => {
                fila.push("—".to_string());
                if predicho.completo() {
                    fila.push("sin gen anotado".to_string());
                } else {
                    fila.push("cortado por el borde".to_string());
                }
            }
        }
    }
    fila
}

fn fraccion(parte: usize, total: usize) -> f64 {
    parte as f64 / total.max(1) as f64
}

fn tabla_de_comparacion(comparacion: &Comparacion) -> Bloque {
    let filas = vec![
        vec![
            "Genes anotados completos en la región".to_string(),
            miles(comparacion.anotados),
            "—".to_string(),
        ],
        vec![
            "Anotados de ≥ 300 pb hallados por su codón de parada".to_string(),
            miles(comparacion.largos_encontrados),
            porcentaje(
                fraccion(comparacion.largos_encontrados, comparacion.anotados_largos),
                1,
            ),
        ],
        vec![
            "Anotados de ≥ 300 pb con el inicio exacto".to_string(),
            miles(comparacion.largos_con_inicio_exacto),
            porcentaje(
                fraccion(
                    comparacion.largos_con_inicio_exacto,
                    comparacion.anotados_largos,
                ),
                1,
            ),
        ],
        vec![
            "Todos los anotados hallados por su codón de parada".to_string(),
            miles(comparacion.encontrados),
            porcentaje(fraccion(comparacion.encontrados, comparacion.anotados), 1),
        ],
        vec![
            "Todos los anotados con el inicio exacto".to_string(),
            miles(comparacion.con_inicio_exacto),
            porcentaje(
                fraccion(comparacion.con_inicio_exacto, comparacion.anotados),
                1,
            ),
        ],
        vec![
            "Predicciones completas sin gen anotado".to_string(),
            miles(comparacion.predichos_sin_anotacion),
            porcentaje(
                fraccion(
                    comparacion.predichos_sin_anotacion,
                    comparacion.predichos_completos,
                ),
                1,
            ),
        ],
    ];
    Bloque::tabla(
        "Comparación con la anotación de RefSeq (NC_000913.3)",
        &["Medida", "Genes", "Porcentaje"],
        filas,
    )
}

fn tabla_de_no_encontrados(
    anotados: &[(&'static Gen, Anotado)],
    comparacion: &Comparacion,
) -> Bloque {
    let filas: Vec<Vec<String>> = anotados
        .iter()
        .zip(&comparacion.hallado)
        .filter(|(_, encontrado)| !**encontrado)
        .map(|((anotado, _), _)| {
            vec![
                anotado.nombre.clone(),
                anotado.locus.clone(),
                miles(anotado.inicio),
                miles(anotado.fin),
                (anotado.hebra as char).to_string(),
                miles(anotado.fin + 1 - anotado.inicio),
                anotado.producto.clone(),
            ]
        })
        .collect();
    Bloque::tabla(
        "Genes anotados que no se encontraron (posiciones en E. coli, contando desde 1)",
        &[
            "Gen",
            "Locus",
            "Extremo izquierdo",
            "Extremo derecho",
            "Hebra",
            "Longitud (pb)",
            "Producto",
        ],
        filas,
    )
}

fn barras_de_longitud(predichos: &[GenPredicho]) -> Bloque {
    let mut categorias = Vec::new();
    let mut valores = Vec::new();
    for &(menor, mayor) in CLASES_DE_LONGITUD {
        let cuantos = predichos
            .iter()
            .filter(|predicho| (menor..mayor).contains(&predicho.longitud()))
            .count();
        if mayor == usize::MAX {
            categorias.push(format!("≥ {}", miles(menor)));
        } else {
            categorias.push(format!("{}–{}", miles(menor), miles(mayor - 1)));
        }
        valores.push(cuantos as f64);
    }
    Bloque::barras(
        "Longitud de los genes predichos",
        "Longitud (pb)",
        "Genes",
        categorias,
        valores,
    )
}

fn barras_de_inicios(cuentas: [usize; 3], total: usize) -> Bloque {
    Bloque::barras(
        "Codón de inicio de los genes con inicio visible",
        "Codón",
        "% de los genes",
        CODONES_DE_INICIO
            .iter()
            .map(|codon| codon.to_string())
            .collect(),
        cuentas
            .iter()
            .map(|&cuenta| (1000.0 * fraccion(cuenta, total)).round() / 10.0)
            .collect(),
    )
}

fn predecir(entrada: &Entrada) -> Result<Salida, String> {
    let adn = leer_adn(entrada)?;
    let prediccion = genes::predecir_genes(&adn, modelo_ecoli());
    let ubicacion = ubicar_en_ecoli(&adn);
    let nombre = match &ubicacion {
        Some(lugar) => lugar.nombre(),
        None => "secuencia".to_string(),
    };
    let ruta: PathBuf =
        crate::carpeta_de_ejemplos()?.join(nombre_de_archivo(&adn, ubicacion.as_ref()));
    guardar_gff(&ruta, &nombre, adn.len(), &prediccion.genes)?;
    let predichos = &prediccion.genes;
    let gc = contenido_gc(&adn);
    let densidad = densidad_codificante(predichos, adn.len());
    let longitud_media = predichos.iter().map(|p| p.longitud()).sum::<usize>() as f64
        / predichos.len().max(1) as f64;
    let (cuentas, con_inicio) = uso_de_inicios(predichos);

    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Genes predichos", miles(predichos.len())));
    salida.agregar(Bloque::valor(
        "Densidad codificante",
        porcentaje(densidad, 1),
    ));
    salida.agregar(Bloque::valor(
        "Longitud media",
        format!("{} pb", miles(longitud_media.round() as usize)),
    ));
    salida.agregar(Bloque::valor("GC de la secuencia", porcentaje(gc, 1)));
    let modelo = if prediccion.entrenado_con_la_entrada {
        format!(
            "entrenado con tu secuencia ({} de ≥ 300 pb)",
            plural(prediccion.modelo.genes_de_entrenamiento, "gen", "genes")
        )
    } else {
        "el de E. coli K-12 (tu secuencia mide menos de 20 kb)".to_string()
    };
    salida.agregar(Bloque::valor("Modelo", modelo));

    let anotados = ubicacion
        .as_ref()
        .map(anotados_en_region)
        .unwrap_or_default();
    let solo_anotados: Vec<Anotado> = anotados.iter().map(|(_, anotado)| *anotado).collect();
    let comparacion = ubicacion
        .as_ref()
        .map(|_| genes::comparar_con_anotacion(predichos, &solo_anotados));

    let mut columnas = vec![
        "N.º",
        "Extremo izquierdo (contando desde 1)",
        "Extremo derecho",
        "Hebra",
        "Longitud (pb)",
        "Aminoácidos",
        "Codón de inicio",
        "RBS",
        "Espaciador (pb)",
        "Puntaje",
        "¿Parcial?",
    ];
    if ubicacion.is_some() {
        columnas.extend(["En E. coli (desde 1)", "Gen de RefSeq", "Comparación"]);
    }
    let filas: Vec<Vec<String>> = predichos
        .iter()
        .enumerate()
        .map(|(indice, predicho)| {
            let referencia = match (&ubicacion, &comparacion) {
                (Some(lugar), Some(comparada)) => {
                    let pareja = comparada.pareja[indice].map(|j| anotados[j]);
                    Some((lugar, pareja))
                }
                _ => None,
            };
            fila_del_gen(indice + 1, predicho, referencia)
        })
        .collect();
    salida.agregar(Bloque::tabla("Genes predichos", &columnas, filas));
    if let Some(comparada) = &comparacion {
        salida.agregar(tabla_de_comparacion(comparada));
        if comparada.encontrados < comparada.anotados {
            salida.agregar(tabla_de_no_encontrados(&anotados, comparada));
        }
    }
    salida.agregar(barras_de_longitud(predichos));
    salida.agregar(barras_de_inicios(cuentas, con_inicio));
    salida.agregar(Bloque::archivo("Genes en formato GFF3", &ruta));

    let medidas = Medidas {
        gc,
        densidad,
        longitud_media,
        cuentas,
        con_inicio,
        corta: adn.len() < SECUENCIA_MINIMA_PARA_ENTRENAR,
        longitud: adn.len(),
    };
    explicar(&mut salida, &prediccion, &medidas, comparacion.as_ref());
    Ok(salida)
}

struct Medidas {
    gc: f64,
    densidad: f64,
    longitud_media: f64,
    cuentas: [usize; 3],
    con_inicio: usize,
    corta: bool,
    longitud: usize,
}

fn plural(cantidad: usize, singular: &str, plural: &str) -> String {
    if cantidad == 1 {
        format!("1 {singular}")
    } else {
        format!("{} {plural}", miles(cantidad))
    }
}

fn explicar(
    salida: &mut Salida,
    prediccion: &Prediccion,
    medidas: &Medidas,
    comparacion: Option<&Comparacion>,
) {
    let gc = medidas.gc;
    let gc_ecoli = modelo_ecoli().gc;
    let distinto_de_ecoli =
        medidas.longitud >= LONGITUD_PARA_OPINAR_DEL_GC && (gc - gc_ecoli).abs() > 0.05;
    let entrenamiento = prediccion.modelo.genes_de_entrenamiento;
    if prediccion.entrenado_con_la_entrada {
        let mut texto = format!(
            "El modelo aprendió de tu propia secuencia: usó {} de ≥ {LONGITUD_PARA_ENTRENAR} pb para saber qué hexámeros y qué sitios de unión al ribosoma (RBS) tienen sus genes.",
            plural(entrenamiento, "gen", "genes")
        );
        if entrenamiento < MINIMO_DE_GENES_PARA_OPINAR {
            texto.push_str(" Son muy pocos: el modelo queda débil, y una bacteria de este tamaño tendría muchos más; quizá la secuencia no es bacteriana.");
        }
        if distinto_de_ecoli {
            texto.push_str(&format!(
                " Su GC ({}) se aleja más de 5 puntos del de E. coli ({}): es otro tipo de genoma, o una región que llegó por transferencia horizontal.",
                porcentaje(gc, 1),
                porcentaje(gc_ecoli, 1)
            ));
        }
        salida.explicar(texto);
    } else {
        let mut texto = format!(
            "Tu secuencia mide menos de {} kb, muy poco para entrenar un modelo propio, así que se usó el de E. coli K-12.",
            SECUENCIA_MINIMA_PARA_ENTRENAR / 1000
        );
        if distinto_de_ecoli {
            texto.push_str(&format!(
                " Ojo: su GC ({}) es muy distinto del de E. coli ({}), así que el modelo puede equivocarse; si puedes, analiza el genoma completo.",
                porcentaje(gc, 1),
                porcentaje(gc_ecoli, 1)
            ));
        }
        salida.explicar(texto);
    }

    let predichos = &prediccion.genes;
    if predichos.is_empty() {
        salida.explicar(format!(
            "No apareció ningún marco abierto de al menos {LONGITUD_MINIMA} pb con aspecto de gen. Puede ser ADN no codificante, de un eucariota (con intrones) o una secuencia demasiado corta."
        ));
        return;
    }

    let densidad = medidas.densidad;
    let densidad_texto = porcentaje(densidad, 1);
    if medidas.corta {
        salida.explicar(format!(
            "{} el {densidad_texto} de la secuencia. En menos de {} kb este porcentaje depende de dónde empieza y termina el tramo, así que no lo compares con el 85–90 % típico de un genoma bacteriano completo.",
            if predichos.len() == 1 {
                "1 gen cubre".to_string()
            } else {
                format!("{} genes cubren", miles(predichos.len()))
            },
            SECUENCIA_MINIMA_PARA_ENTRENAR / 1000
        ));
    } else if (0.85..=0.90).contains(&densidad) {
        salida.explicar(format!(
            "La densidad codificante es {densidad_texto}: lo típico de una bacteria (85–90 %), cuyos genomas son compactos y casi no tienen ADN entre genes."
        ));
    } else if densidad > 0.90 {
        salida.explicar(format!(
            "La densidad codificante ({densidad_texto}) supera el 85–90 % típico de las bacterias: es una región muy apretada (como un operón) o hay genes falsos solapados."
        ));
    } else if densidad >= 0.5 {
        salida.explicar(format!(
            "La densidad codificante ({densidad_texto}) queda bajo el 85–90 % típico de las bacterias: puede haber genes de ARN (ribosomal o de transferencia), repeticiones o genes que el modelo no reconoció."
        ));
    } else {
        salida.explicar(format!(
            "Menos de la mitad de la secuencia codifica proteínas ({densidad_texto}), muy lejos del 85–90 % de las bacterias: quizá no es ADN bacteriano o es una región con pocos genes."
        ));
    }

    if predichos.len() >= MINIMO_DE_GENES_PARA_OPINAR {
        let longitud_media = medidas.longitud_media;
        let media = miles(longitud_media.round() as usize);
        if (800.0..=1100.0).contains(&longitud_media) {
            salida.explicar(format!(
                "Los genes miden {media} pb en promedio, cerca de los ~{LONGITUD_MEDIA_ECOLI} pb de E. coli."
            ));
        } else if longitud_media < 800.0 {
            salida.explicar(format!(
                "Los genes miden {media} pb en promedio, menos que los ~{LONGITUD_MEDIA_ECOLI} pb de E. coli: sobran genes cortos o hay genes cortados por los bordes."
            ));
        } else {
            salida.explicar(format!(
                "Los genes miden {media} pb en promedio, más que los ~{LONGITUD_MEDIA_ECOLI} pb de E. coli: puede haber inicios elegidos demasiado arriba o pocos genes cortos."
            ));
        }
    }

    let con_inicio = medidas.con_inicio;
    if con_inicio >= MINIMO_DE_GENES_PARA_OPINAR {
        let atg = fraccion(medidas.cuentas[0], con_inicio);
        let resumen = format!(
            "Usan ATG en {}, GTG en {} y TTG en {} de los casos (E. coli: ≈ {} ATG, {} GTG, {} TTG)",
            porcentaje(atg, 0),
            porcentaje(fraccion(medidas.cuentas[1], con_inicio), 0),
            porcentaje(fraccion(medidas.cuentas[2], con_inicio), 0),
            porcentaje(USO_DE_INICIOS_ECOLI[0], 0),
            porcentaje(USO_DE_INICIOS_ECOLI[1], 0),
            porcentaje(USO_DE_INICIOS_ECOLI[2], 0),
        );
        if atg >= 0.75 {
            salida.explicar(format!(
                "{resumen}. ATG domina, como en casi todas las bacterias."
            ));
        } else if gc > 0.60 {
            salida.explicar(format!(
                "{resumen}. Tantos GTG y TTG son comunes en genomas ricos en GC como el tuyo."
            ));
        } else {
            salida.explicar(format!(
                "{resumen}. Con un GC moderado, tantos GTG y TTG sugieren inicios mal elegidos o ADN que en realidad no codifica."
            ));
        }
    }

    let cortos_sin_rbs = predichos
        .iter()
        .filter(|p| p.completo() && p.longitud() < LONGITUD_PARA_ENTRENAR && p.rbs.is_none())
        .count();
    if cortos_sin_rbs > 0 {
        let cuantos = if cortos_sin_rbs == 1 {
            "1 gen mide menos de 300 pb y no tiene RBS: es la predicción menos fiable; confírmala"
                .to_string()
        } else {
            format!(
                "{} genes miden menos de 300 pb y no tienen RBS: son las predicciones menos fiables; confírmalas",
                miles(cortos_sin_rbs)
            )
        };
        salida.explicar(format!("{cuantos} buscando su proteína con BLAST."));
    }

    if let Some(comparada) = comparacion {
        if comparada.anotados_largos == 0 {
            salida.explicar(
                "Esta región de E. coli no tiene genes anotados de 300 pb o más con los que comparar.",
            );
        } else {
            let con_parada = fraccion(comparada.largos_encontrados, comparada.anotados_largos);
            let con_inicio = fraccion(
                comparada.largos_con_inicio_exacto,
                comparada.anotados_largos,
            );
            let cierre = if con_parada < 0.9 {
                "Se escaparon más de 1 de cada 10 genes largos: en una región corta el modelo tiene pocos datos y los genes que tocan el borde quedan cortados."
            } else if con_inicio + 0.05 < con_parada {
                "La parada casi siempre se acierta; lo difícil es el inicio, porque varios ATG, GTG o TTG compiten en el mismo marco."
            } else {
                "Casi todos los genes que encontraste tienen también el inicio correcto."
            };
            salida.explicar(format!(
                "Frente a RefSeq hallaste {} de los genes anotados de ≥ 300 pb por su codón de parada y {} con el inicio exacto; {} de tus predicciones completas no tienen gen anotado. {cierre}",
                porcentaje(con_parada, 1),
                porcentaje(con_inicio, 1),
                porcentaje(
                    fraccion(
                        comparada.predichos_sin_anotacion,
                        comparada.predichos_completos
                    ),
                    1
                ),
            ));
        }
    }
}

const MAXIMO_DE_PASOS: usize = 400;
const CODONES_DETALLADOS: usize = 3;
const CODONES_A_RECORRER: usize = 400;
const PASOS_DE_CODONES: usize = 14;
const FILAS_DE_TABLA: usize = 5;
const GENES_EN_LISTA: usize = 8;

struct Variables(Vec<(&'static str, String)>);

impl Variables {
    fn nueva() -> Variables {
        Variables(Vec::new())
    }

    fn poner(&mut self, nombre: &'static str, valor: impl ToString) {
        let texto = valor.to_string();
        match self.0.iter_mut().find(|(clave, _)| *clave == nombre) {
            Some(existente) => existente.1 = texto,
            None => self.0.push((nombre, texto)),
        }
    }

    fn pizarra(&self) -> Pizarra {
        let mut pizarra = Pizarra::nueva();
        for (nombre, valor) in &self.0 {
            pizarra = pizarra.variable(nombre, valor);
        }
        pizarra
    }
}

fn palabra_de_hexamero(numero: usize) -> String {
    (0..6)
        .rev()
        .map(|lugar| b"ACGT"[(numero >> (2 * lugar)) & 3] as char)
        .collect()
}

fn texto_de_inicios(inicios: &[(usize, &str)]) -> String {
    let mut partes: Vec<String> = inicios
        .iter()
        .take(4)
        .map(|(posicion, codon)| {
            let nombre = if codon.is_empty() { "borde" } else { codon };
            format!("{nombre}@{}", miles(posicion + 1))
        })
        .collect();
    if inicios.len() > 4 {
        partes.push("…".to_string());
    }
    format!("[{}]", partes.join(", "))
}

fn texto_del_gen(predicho: &GenPredicho) -> String {
    let codon = if predicho.codon_de_inicio.is_empty() {
        "borde"
    } else {
        predicho.codon_de_inicio
    };
    format!(
        "{}–{} ({}) {codon}, puntaje {}",
        miles(predicho.desde + 1),
        miles(predicho.hasta),
        predicho.hebra as char,
        decimal(predicho.puntaje, 1)
    )
}

fn lista_de_genes(genes: &[GenPredicho]) -> Vec<String> {
    genes
        .iter()
        .take(GENES_EN_LISTA)
        .map(texto_del_gen)
        .collect()
}

fn marcas_de_genes(genes: &[GenPredicho]) -> Vec<(usize, usize, Marca)> {
    genes
        .iter()
        .take(GENES_EN_LISTA)
        .enumerate()
        .map(|(posicion, predicho)| {
            let tipo = if posicion == 0 {
                Marca::Actual
            } else {
                Marca::Hallado
            };
            (predicho.desde, predicho.hasta, tipo)
        })
        .collect()
}

fn mejores_hexameros(puntajes: &[f64]) -> Vec<usize> {
    let mut orden: Vec<usize> = (0..puntajes.len()).collect();
    orden.sort_by(|a, b| puntajes[*b].total_cmp(&puntajes[*a]));
    orden.truncate(FILAS_DE_TABLA);
    orden
}

fn resumen_de_ronda(entrenamiento: &[GenPredicho]) -> String {
    plural(entrenamiento.len(), "gen", "genes")
}

fn trazar_genes(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let adn = leer_adn(entrada)?;
    Ok(narrar_genes(&adn, modelo_ecoli()).0)
}

fn narrar_genes(adn: &[u8], respaldo: &Modelo) -> (Vec<Paso>, Prediccion) {
    let mut narracion = Narracion::nueva(MAXIMO_DE_PASOS);
    let mut v = Variables::nueva();
    let hebras = preparar_hebras(adn);
    let longitud = hebras.longitud();
    let bases = &hebras.directa.bases;
    let entrenado_con_la_entrada = adn.len() >= SECUENCIA_MINIMA_PARA_ENTRENAR;

    v.poner("adn", format!("{} pb", miles(adn.len())));
    narracion.agregar(
        v.pizarra()
            .secuencia(bases, &[(0, 12.min(longitud), Marca::Actual)])
            .paso(
                "pub fn predecir_genes(adn: &[u8], respaldo: &Modelo) -> Prediccion {",
                format!(
                    "Recibes un ADN de {} pb. Un gen que fabrica una proteína es un tramo que empieza con un codón de inicio y sigue de tres en tres letras hasta un codón de parada. Vamos a buscar esos tramos, a puntuar cuáles suenan a gen de verdad y a quedarnos con los mejores.",
                    miles(adn.len())
                ),
            ),
    );

    v.poner(
        "hebras",
        format!("directa + inversa ({} pb cada una)", miles(longitud)),
    );
    narracion.agregar(v.pizarra().secuencia(bases, &[(0, 12.min(longitud), Marca::Actual)]).paso(
        "let hebras = preparar_hebras(adn);",
        "El ADN tiene dos hebras y un gen puede estar escrito en cualquiera de las dos. Por eso se prepara la hebra directa (tu secuencia en mayúsculas) y la inversa, que es la otra hebra leída en sentido contrario.",
    ));

    let inversa = &hebras.inversa.bases;
    let cola = String::from_utf8_lossy(&bases[longitud.saturating_sub(9)..]).to_string();
    let cabeza = String::from_utf8_lossy(&inversa[..9.min(longitud)]).to_string();
    narracion.agregar(
        v.pizarra()
            .secuencia(inversa, &[(0, 9.min(longitud), Marca::Actual)])
            .paso(
                "let inversa = complemento_inverso(&directa);",
                format!(
                    "La hebra inversa se arma al revés y cambiando cada letra por su pareja (A con T, C con G). Las últimas 9 bases de tu ADN, «{cola}», se vuelven las primeras de la inversa: «{cabeza}». Aquí ves la hebra inversa."
                ),
            ),
    );

    let primeros: Vec<usize> = hebras
        .directa
        .hexameros
        .iter()
        .enumerate()
        .filter_map(|(posicion, numero)| numero.map(|_| posicion))
        .take(4)
        .collect();
    let con_numero = hebras
        .directa
        .hexameros
        .iter()
        .filter(|numero| numero.is_some())
        .count();
    v.poner("hexameros", format!("{} números", miles(con_numero)));
    let filas: Vec<Vec<String>> = primeros
        .iter()
        .map(|&posicion| {
            let numero = hebras.directa.hexameros[posicion].unwrap_or(0);
            vec![
                miles(posicion + 1),
                String::from_utf8_lossy(&bases[posicion..posicion + 6]).to_string(),
                miles(numero as usize),
            ]
        })
        .collect();
    let marcas: Vec<(usize, usize, Marca)> = primeros
        .iter()
        .enumerate()
        .map(|(orden, &posicion)| {
            let tipo = if orden == 0 {
                Marca::Actual
            } else {
                Marca::Comparado
            };
            (posicion, posicion + 6, tipo)
        })
        .collect();
    narracion.agregar(
        v.pizarra()
            .secuencia(bases, &marcas)
            .tabla(
                "Primeras palabras de 6 letras",
                &["Posición", "Palabra", "Número (0–4 095)"],
                filas,
                &[(0, 2)],
            )
            .paso(
                "hexameros[posicion - 5] = Some(numero as u16);",
                format!(
                    "Cada palabra de 6 letras (hay 4 096 posibles) recibe un número: A=0, C=1, G=2, T=3 leídos como un número en base 4. La ventana avanza una letra a la vez y apunta el número de la palabra que empieza en cada posición; en total {} palabras. Más adelante el modelo dirá qué números suenan a gen.",
                    miles(con_numero)
                ),
            ),
    );

    let marcos_directa = marcos_abiertos(&hebras.directa);
    let marcos_inversa = marcos_abiertos(&hebras.inversa);
    v.poner("marcos", "(se buscan ahora)");
    narracion.agregar(
        v.pizarra()
            .secuencia(bases, &[(0, 12.min(longitud), Marca::Actual)])
            .paso(
                "let marcos = marcos_de_ambas_hebras(&hebras);",
                "Ahora buscamos los marcos abiertos: tramos que se leen de tres en tres sin toparse con una parada. Se prueban tres puntos de partida (desfases 0, 1 y 2) en cada una de las dos hebras: seis formas de leer el mismo ADN. Veamos cómo se hace en la primera.",
            ),
    );

    v.poner("desfase", 0);
    let mut inicios: Vec<(usize, &'static str)> = Vec::new();
    let mut posicion = 0;
    let mut detalle = 0;
    let mut recorridos = 0;
    let mut rbs_mostrado = false;
    let mut parada_hallada = false;
    v.poner("posicion", 0);
    v.poner("codon", "—");
    v.poner("inicios", "[]");
    while posicion + 3 <= bases.len()
        && recorridos < CODONES_A_RECORRER
        && detalle < PASOS_DE_CODONES
    {
        let codon = &bases[posicion..posicion + 3];
        let texto = String::from_utf8_lossy(codon).to_string();
        v.poner("posicion", miles(posicion));
        v.poner("codon", format!("«{texto}»"));
        let mut marcas: Vec<(usize, usize, Marca)> = inicios
            .iter()
            .filter(|(donde, _)| *donde != posicion)
            .map(|(donde, _)| (*donde, donde + 3, Marca::Hallado))
            .collect();
        marcas.insert(0, (posicion, posicion + 3, Marca::Actual));
        if es_codon_de_parada(codon) {
            narracion.agregar(v.pizarra().secuencia(bases, &marcas).paso(
                "if es_codon_de_parada(codon) {",
                format!(
                    "«{texto}» es un codón de parada (TAA, TAG o TGA): ahí termina la proteína. El tramo abierto acaba aquí."
                ),
            ));
            let fin = posicion + 3;
            let sirven = inicios
                .iter()
                .take_while(|(donde, _)| fin - donde >= LONGITUD_MINIMA)
                .count();
            let nota = if inicios.is_empty() {
                "No había ningún inicio antes de esta parada, así que no se guarda nada y se empieza de nuevo.".to_string()
            } else if sirven == 0 {
                format!(
                    "Entre el primer inicio y la parada hay {} pb, menos de {LONGITUD_MINIMA}: es demasiado corto para ser un gen y se descarta.",
                    fin - inicios[0].0
                )
            } else {
                format!(
                    "Se guarda un marco con {} de {} inicios posibles, de {} pb desde el primero. Esos inicios son candidatos a empezar el mismo gen: aún no sabemos cuál es el verdadero.",
                    sirven,
                    inicios.len(),
                    fin - inicios[0].0
                )
            };
            let mut marcas_guardado = vec![(posicion, fin, Marca::Actual)];
            marcas_guardado.extend(
                inicios
                    .iter()
                    .map(|(donde, _)| (*donde, donde + 3, Marca::Hallado)),
            );
            narracion.agregar(v.pizarra().secuencia(bases, &marcas_guardado).paso(
                "guardar_marco(&mut marcos, hebra.signo, &inicios, posicion + 3, true);",
                nota,
            ));
            parada_hallada = true;
            break;
        } else if let Some(nombre) = codon_de_inicio(codon) {
            inicios.push((posicion, nombre));
            v.poner("inicios", texto_de_inicios(&inicios));
            narracion.agregar(v.pizarra().secuencia(bases, &marcas).paso(
                "else if let Some(nombre) = codon_de_inicio(codon)",
                format!(
                    "«{texto}» es un codón de inicio posible (ATG, GTG o TTG). Se apunta su posición: un gen podría empezar aquí."
                ),
            ));
            detalle += 1;
            if !rbs_mostrado {
                rbs_mostrado = true;
                let rbs = buscar_rbs(bases, posicion);
                let mut marcas_rbs = vec![
                    (posicion, posicion + 3, Marca::Actual),
                    (posicion.saturating_sub(21), posicion, Marca::Comparado),
                ];
                let nota = match rbs {
                    Some(encontrado) => {
                        let desde = posicion - encontrado.espaciador - encontrado.motivo.len();
                        marcas_rbs.push((desde, desde + encontrado.motivo.len(), Marca::Hallado));
                        format!(
                            "Antes de cada inicio se busca el RBS, una señal parecida a AGGAGG que el ribosoma usa para agarrarse. Aquí aparece «{}» a {} bases del inicio: eso hace más creíble que sea un inicio de verdad.",
                            encontrado.motivo, encontrado.espaciador
                        )
                    }
                    None => "Antes de cada inicio se busca el RBS, una señal parecida a AGGAGG que el ribosoma usa para agarrarse. Aquí no hay ninguna entre 5 y 15 bases antes: este inicio queda sin RBS y valdrá menos.".to_string(),
                };
                narracion.agregar(
                    v.pizarra()
                        .secuencia(bases, &marcas_rbs)
                        .paso("rbs: buscar_rbs(bases, posicion),", nota),
                );
                detalle += 1;
            }
        } else if posicion == 0 {
            inicios.push((0, ""));
            v.poner("inicios", texto_de_inicios(&inicios));
            narracion.agregar(v.pizarra().secuencia(bases, &marcas).paso(
                "else if posicion == desfase {",
                format!(
                    "«{texto}» no es inicio ni parada, pero es el primer codón de la secuencia. El gen pudo empezar antes del borde de lo que tienes, así que se apunta como «inicio de borde»."
                ),
            ));
            detalle += 1;
        } else if recorridos < CODONES_DETALLADOS {
            narracion.agregar(v.pizarra().secuencia(bases, &marcas).paso(
                "posicion += 3;",
                format!(
                    "«{texto}» no es inicio ni parada, así que no se apunta nada. Se avanzan tres letras y se mira el siguiente codón."
                ),
            ));
            detalle += 1;
        }
        posicion += 3;
        recorridos += 1;
    }
    if !parada_hallada {
        narracion.agregar(v.pizarra().secuencia(bases, &[(posicion.min(longitud - 3), posicion.min(longitud - 3) + 3, Marca::Actual)]).paso(
            "while posicion + 3 <= bases.len() {",
            format!(
                "Así se recorre todo el ADN, codón por codón. Aquí dejamos de detallar tras {} codones sin llegar a una parada; el programa sigue igual hasta el final.",
                miles(recorridos)
            ),
        ));
    }
    narracion.agregar(v.pizarra().secuencia(bases, &[(0, 12.min(longitud), Marca::Actual)]).paso(
        "for desfase in 0..3 {",
        "Lo mismo se repite empezando en la segunda y en la tercera letra (desfases 1 y 2), y luego en la hebra inversa. Así no se escapa ningún gen, esté en la hebra que esté.",
    ));

    let mut marcos: Vec<Marco> = marcos_directa.clone();
    marcos.extend(marcos_inversa.iter().cloned());
    v.poner(
        "marcos",
        format!(
            "{} ({} en +, {} en −)",
            miles(marcos.len()),
            miles(marcos_directa.len()),
            miles(marcos_inversa.len())
        ),
    );
    let mut por_largo: Vec<&Marco> = marcos
        .iter()
        .filter(|marco| !marco.inicios.is_empty())
        .collect();
    por_largo.sort_by_key(|marco| std::cmp::Reverse(marco.fin - marco.inicios[0].posicion));
    let filas: Vec<Vec<String>> = por_largo
        .iter()
        .take(FILAS_DE_TABLA)
        .map(|marco| {
            vec![
                (marco.hebra as char).to_string(),
                miles(marco.inicios[0].posicion + 1),
                miles(marco.fin),
                miles(marco.fin - marco.inicios[0].posicion),
                miles(marco.inicios.len()),
                if marco.tiene_parada { "sí" } else { "no" }.to_string(),
            ]
        })
        .collect();
    let con_inicios: usize = marcos.iter().map(|marco| marco.inicios.len()).sum();
    let (texto_secuencia, marcas): (&[u8], Vec<(usize, usize, Marca)>) = match por_largo.first() {
        Some(marco) => {
            let mut marcas = vec![(
                marco.inicios[0].posicion,
                marco.inicios[0].posicion + 3,
                Marca::Actual,
            )];
            marcas.extend(
                marco
                    .inicios
                    .iter()
                    .skip(1)
                    .take(10)
                    .map(|inicio| (inicio.posicion, inicio.posicion + 3, Marca::Hallado)),
            );
            (&hebras.de(marco.hebra).bases, marcas)
        }
        None => (bases, Vec::new()),
    };
    narracion.agregar(
        v.pizarra()
            .secuencia(texto_secuencia, &marcas)
            .tabla(
                "Los marcos abiertos más largos",
                &[
                    "Hebra",
                    "Primer inicio",
                    "Fin",
                    "Largo (pb)",
                    "Inicios posibles",
                    "¿Con parada?",
                ],
                filas,
                &[(0, 3)],
            )
            .paso(
                "marcos.extend(marcos_abiertos(&hebras.inversa));",
                format!(
                    "Resultado de la búsqueda: {} marcos abiertos con al menos {LONGITUD_MINIMA} pb ({} en la hebra directa y {} en la inversa), con {} inicios posibles en total. La mayoría no son genes: en ADN al azar aparecen tramos abiertos por casualidad. Falta decidir cuáles son reales.",
                    miles(marcos.len()),
                    miles(marcos_directa.len()),
                    miles(marcos_inversa.len()),
                    miles(con_inicios)
                ),
            ),
    );

    v.poner("entrenado_con_la_entrada", entrenado_con_la_entrada);
    narracion.agregar(v.pizarra().secuencia(bases, &[(0, 12.min(longitud), Marca::Actual)]).paso(
        "let entrenado_con_la_entrada = adn.len() >= SECUENCIA_MINIMA_PARA_ENTRENAR;",
        if entrenado_con_la_entrada {
            format!(
                "Tu ADN mide {} pb, al menos {} pb: hay suficiente para que el modelo aprenda de tu propia secuencia. Así se adapta a cualquier bacteria, no solo a E. coli.",
                miles(adn.len()),
                miles(SECUENCIA_MINIMA_PARA_ENTRENAR)
            )
        } else {
            format!(
                "Tu ADN mide {} pb, menos de {} pb: es poco para aprender de él. Se usará el modelo que ya se entrenó con E. coli K-12.",
                miles(adn.len()),
                miles(SECUENCIA_MINIMA_PARA_ENTRENAR)
            )
        },
    ));

    let modelo = if entrenado_con_la_entrada {
        v.poner("modelo", "(se entrena ahora)");
        narracion.agregar(v.pizarra().secuencia(bases, &[(0, 12.min(longitud), Marca::Actual)]).paso(
            "let modelo = if entrenado_con_la_entrada {",
            "Entrenar es descubrir qué hace distinto al ADN de los genes del resto. Se hace en tres rondas: primero con los marcos más largos, que casi seguro son genes, y luego se repite con los genes que el propio modelo va eligiendo.",
        ));
        entrenar_narrado(&hebras, &marcos, &mut narracion, &mut v)
    } else {
        v.poner(
            "modelo",
            format!(
                "el de E. coli ({})",
                plural(respaldo.genes_de_entrenamiento, "gen", "genes")
            ),
        );
        narracion.agregar(v.pizarra().secuencia(bases, &[(0, 12.min(longitud), Marca::Actual)]).paso(
            "respaldo.clone()",
            format!(
                "Se toma una copia del modelo de E. coli, que aprendió de {} de su genoma. Su GC es {}. Cuanto más se parezca tu ADN al de E. coli, mejor funcionará.",
                plural(respaldo.genes_de_entrenamiento, "gen", "genes"),
                porcentaje(respaldo.gc, 1)
            ),
        ));
        respaldo.clone()
    };

    v.poner("candidatos", "(se puntúan ahora)");
    narracion.agregar(v.pizarra().secuencia(bases, &[(0, 12.min(longitud), Marca::Actual)]).paso(
        "let genes = genes_con_modelo(&hebras, &marcos, &modelo);",
        "Con el modelo listo se puntúa cada inicio posible de cada marco. El puntaje dice qué tan bien suena el tramo a gen y qué tan creíble es ese inicio.",
    ));

    let mut candidatos: Vec<GenPredicho> = Vec::new();
    let mut rangos = Vec::with_capacity(marcos.len());
    for marco in &marcos {
        let antes = candidatos.len();
        candidatos_del_marco(marco, hebras.de(marco.hebra), &modelo, &mut candidatos);
        rangos.push((antes, candidatos.len()));
    }
    let ejemplo = rangos
        .iter()
        .position(|(antes, despues)| despues - antes >= 2)
        .or_else(|| rangos.iter().position(|(antes, despues)| despues > antes));
    let elegidos = elegir_genes(&candidatos);
    let mut genes: Vec<GenPredicho> = elegidos.iter().map(|&indice| candidatos[indice]).collect();
    genes.sort_by_key(|predicho| (predicho.desde, predicho.hasta));

    if let Some(indice) = ejemplo {
        let marco = &marcos[indice];
        let (antes, despues) = rangos[indice];
        let sus = &candidatos[antes..despues];
        let hebra_del_marco = hebras.de(marco.hebra);
        let posicion_de = |candidato: &GenPredicho| {
            if candidato.hebra == b'-' {
                longitud - candidato.hasta
            } else {
                candidato.desde
            }
        };
        let mut marcas: Vec<(usize, usize, Marca)> = marco
            .inicios
            .iter()
            .take(12)
            .enumerate()
            .map(|(orden, inicio)| {
                let tipo = if orden == 0 {
                    Marca::Actual
                } else {
                    Marca::Comparado
                };
                (inicio.posicion, inicio.posicion + 3, tipo)
            })
            .collect();
        if marco.tiene_parada {
            marcas.push((marco.fin - 3, marco.fin, Marca::Hallado));
        }
        v.poner("candidatos", "[]");
        narracion.agregar(v.pizarra().secuencia(&hebra_del_marco.bases, &marcas).paso(
            "candidatos_del_marco(marco, hebras.de(marco.hebra), modelo, &mut candidatos);",
            format!(
                "Tomemos un marco de ejemplo en la hebra {}: tiene {} inicios posibles ({}) que comparten la misma parada. Cada inicio daría un gen de distinto largo, así que se puntúa cada uno por separado.",
                marco.hebra as char,
                miles(marco.inicios.len()),
                CODONES_DE_INICIO.join(", ")
            ),
        ));

        let filas_de_codificacion: Vec<Vec<String>> = sus
            .iter()
            .map(|candidato| {
                vec![
                    miles(posicion_de(candidato) + 1),
                    candidato.codon_de_inicio.to_string(),
                    decimal(candidato.puntaje_codificacion, 1),
                ]
            })
            .collect();
        v.poner("codificacion", "se suma");
        narracion.agregar(
            v.pizarra()
                .secuencia(&hebra_del_marco.bases, &marcas)
                .tabla(
                    "Puntaje de codificación por inicio",
                    &["Inicio (posición)", "Codón", "Codificación"],
                    filas_de_codificacion,
                    &[(0, 2)],
                )
                .paso(
                    "codificacion += modelo.codificacion[numero as usize];",
                    "Desde la parada hacia atrás se leen las palabras de 6 letras de tres en tres y se suma el puntaje que el modelo da a cada una. Un inicio más lejano suma más palabras: si esas palabras suenan a gen, el puntaje sube.",
                ),
        );

        let filas_de_inicio: Vec<Vec<String>> = sus
            .iter()
            .map(|candidato| {
                let rbs = match candidato.rbs {
                    Some(encontrado) => {
                        format!("{} a {} pb", encontrado.motivo, encontrado.espaciador)
                    }
                    None => "sin RBS".to_string(),
                };
                vec![
                    miles(posicion_de(candidato) + 1),
                    candidato.codon_de_inicio.to_string(),
                    rbs,
                    decimal(candidato.puntaje_inicio, 2),
                ]
            })
            .collect();
        v.poner("puntaje_inicio", "por inicio");
        narracion.agregar(
            v.pizarra()
                .secuencia(&hebra_del_marco.bases, &marcas)
                .tabla(
                    "Qué tan creíble es cada inicio",
                    &["Inicio (posición)", "Codón", "RBS", "Puntaje del inicio"],
                    filas_de_inicio,
                    &[(0, 3)],
                )
                .paso(
                    "let puntaje_inicio = puntaje_del_inicio(modelo, inicio.codon, inicio.rbs);",
                    "Cada inicio recibe además un puntaje propio: el modelo suma lo que aprendió sobre su codón (ATG es el más usado) y sobre su RBS. Un RBS fuerte y un ATG suben el puntaje.",
                ),
        );

        let filas_de_total: Vec<Vec<String>> = sus
            .iter()
            .map(|candidato| {
                vec![
                    miles(posicion_de(candidato) + 1),
                    decimal(candidato.puntaje_codificacion, 1),
                    decimal(genes::PESO_DEL_INICIO * candidato.puntaje_inicio, 1),
                    decimal(candidato.puntaje, 1),
                ]
            })
            .collect();
        v.poner("puntaje", "por inicio");
        narracion.agregar(
            v.pizarra()
                .secuencia(&hebra_del_marco.bases, &marcas)
                .tabla(
                    "Puntaje total de cada candidato",
                    &["Inicio (posición)", "Codificación", "3 × inicio", "Total"],
                    filas_de_total,
                    &[(0, 3)],
                )
                .paso(
                    "let puntaje = codificacion + PESO_DEL_INICIO * puntaje_inicio;",
                    format!(
                        "El puntaje total es lo que suena a gen más {} veces lo creíble que es el inicio. Un puntaje alto significa un candidato fuerte.",
                        decimal(genes::PESO_DEL_INICIO, 0)
                    ),
                ),
        );
        let descartados = marco.inicios.len() - sus.len();
        if descartados > 0 {
            narracion.agregar(v.pizarra().secuencia(&hebra_del_marco.bases, &marcas).paso(
                "if codificacion <= 0.0 || puntaje <= 0.0 {",
                format!(
                    "{} de los {} inicios de este marco tienen puntaje de codificación o total cero o negativo: suenan más a ADN cualquiera que a gen, y se descartan.",
                    miles(descartados),
                    miles(marco.inicios.len())
                ),
            ));
        }
    }

    let total_de_inicios: usize = marcos.iter().map(|marco| marco.inicios.len()).sum();
    v.poner(
        "candidatos",
        plural(candidatos.len(), "candidato", "candidatos"),
    );
    narracion.agregar(
        v.pizarra()
            .secuencia(bases, &[(0, 12.min(longitud), Marca::Actual)])
            .paso(
                "for marco in marcos {",
                format!(
                    "Se repite lo mismo en los {} marcos. De {} inicios posibles quedaron {} candidatos con puntaje positivo; el resto se descartó.",
                    miles(marcos.len()),
                    miles(total_de_inicios),
                    miles(candidatos.len())
                ),
            ),
    );

    let en_conflicto = ejemplo
        .map(|indice| {
            let (antes, despues) = rangos[indice];
            let filas: Vec<Vec<String>> = candidatos[antes..despues]
                .iter()
                .enumerate()
                .map(|(orden, candidato)| {
                    let elegido = elegidos.contains(&(antes + orden));
                    vec![
                        texto_del_gen(candidato),
                        if elegido { "sí" } else { "no" }.to_string(),
                    ]
                })
                .collect();
            let resaltar: Vec<(usize, usize)> = filas
                .iter()
                .enumerate()
                .filter(|(_, fila)| fila[1] == "sí")
                .map(|(orden, _)| (orden, 1))
                .collect();
            (filas, resaltar)
        })
        .unwrap_or_default();
    let mut tabla_de_eleccion = v.pizarra().secuencia(bases, &marcas_de_genes(&genes));
    if !en_conflicto.0.is_empty() {
        tabla_de_eleccion = tabla_de_eleccion.tabla(
            "Los candidatos del marco de ejemplo",
            &["Candidato", "¿Elegido?"],
            en_conflicto.0,
            &en_conflicto.1,
        );
    }
    narracion.agregar(tabla_de_eleccion.paso(
        "let mut genes: Vec<GenPredicho> = elegir_genes(&candidatos)",
        "Muchos candidatos se pisan: los de un mismo marco comparten parada y casi todos se solapan. Hay que escoger el conjunto que no choque y que sume el mayor puntaje total, como armar el mejor equipo sin que dos jugadores ocupen el mismo puesto.",
    ));
    v.poner("elegidos", plural(elegidos.len(), "gen", "genes"));
    narracion.agregar(
        v.pizarra()
            .secuencia(bases, &marcas_de_genes(&genes))
            .paso(
                "cadena[actual] = suma + candidato.puntaje;",
                format!(
                    "Los candidatos se ordenan por dónde terminan. Cada uno suma su puntaje al mejor encadenamiento de candidatos anteriores con los que no choca (solo se permite un pequeño solapamiento). Al final de ese recorrido se siguen las cadenas hacia atrás y quedan {} de los {} candidatos.",
                    miles(elegidos.len()),
                    miles(candidatos.len())
                ),
            ),
    );
    v.poner("genes", plural(genes.len(), "gen", "genes"));
    narracion.agregar(
        v.pizarra()
            .secuencia(bases, &marcas_de_genes(&genes))
            .lista("Primeros genes elegidos", lista_de_genes(&genes), &[0])
            .paso(
                "genes.sort_by_key(|predicho| (predicho.desde, predicho.hasta));",
                "Los genes elegidos se ordenan por su posición en el ADN, de izquierda a derecha. Esta es la lista final, con las posiciones contadas desde 1.",
            ),
    );

    let predichos = &genes;
    let densidad = densidad_codificante(predichos, adn.len());
    let longitud_media = predichos.iter().map(|p| p.longitud()).sum::<usize>() as f64
        / predichos.len().max(1) as f64;
    let parciales = predichos.iter().filter(|p| !p.completo()).count();
    let prediccion = Prediccion {
        genes,
        modelo,
        entrenado_con_la_entrada,
    };
    let ultimo = v
        .pizarra()
        .secuencia(bases, &marcas_de_genes(&prediccion.genes))
        .lista(
            "Primeros genes",
            lista_de_genes(&prediccion.genes),
            &[0],
        )
        .paso(
            "Prediccion {",
            format!(
                "Resultado: {} predichos, que cubren el {} del ADN, con {} pb de longitud media y {} cortados por el borde. La tabla de resultados los muestra todos; en «Qué significa» ves si esos números son normales para una bacteria.",
                plural(prediccion.genes.len(), "gen", "genes"),
                porcentaje(densidad, 1),
                miles(longitud_media.round() as usize),
                miles(parciales)
            ),
        );
    (narracion.terminar(ultimo), prediccion)
}

fn entrenar_narrado(
    hebras: &Hebras,
    marcos: &[Marco],
    narracion: &mut Narracion,
    v: &mut Variables,
) -> Modelo {
    let bases = &hebras.directa.bases;
    let inicio_de_vista = [(0, 12.min(bases.len()), Marca::Actual)];
    let fondo = hexameros_de_fondo(hebras);
    let total_de_fondo: f64 = fondo.iter().sum();
    v.poner(
        "fondo",
        format!("{} palabras contadas", miles(total_de_fondo as usize)),
    );
    let mejores: Vec<usize> = {
        let mut orden: Vec<usize> = (0..fondo.len()).collect();
        orden.sort_by(|a, b| fondo[*b].total_cmp(&fondo[*a]));
        orden.truncate(FILAS_DE_TABLA);
        orden
    };
    let filas: Vec<Vec<String>> = mejores
        .iter()
        .map(|&numero| {
            vec![
                palabra_de_hexamero(numero),
                miles(fondo[numero] as usize),
                porcentaje(fondo[numero] / total_de_fondo.max(1.0), 2),
            ]
        })
        .collect();
    narracion.agregar(
        v.pizarra()
            .secuencia(bases, &inicio_de_vista)
            .tabla(
                "Las palabras más comunes en todo el ADN",
                &["Palabra", "Veces", "Fracción"],
                filas,
                &[(0, 1)],
            )
            .paso(
                "let fondo = hexameros_de_fondo(hebras);",
                format!(
                    "Primero se cuentan todas las palabras de 6 letras de las dos hebras: {} en total. Eso es el «fondo»: cómo se ve el ADN en general, con genes y sin ellos.",
                    miles(total_de_fondo as usize)
                ),
            ),
    );

    let mut modelo = Modelo::neutro(nucleo::secuencia::contenido_gc(bases));
    v.poner("modelo", format!("vacío (GC {})", porcentaje(modelo.gc, 1)));
    narracion.agregar(
        v.pizarra().secuencia(bases, &inicio_de_vista).paso(
            "let mut modelo = Modelo::neutro(contenido_gc(&hebras.directa.bases));",
            format!(
                "El modelo empieza vacío: todos sus puntajes valen cero. Solo se anota el GC de tu ADN, {}, la proporción de letras G y C.",
                porcentaje(modelo.gc, 1)
            ),
        ),
    );

    let mut entrenamiento = marcos_mas_largos(hebras, marcos);
    v.poner("entrenamiento", resumen_de_ronda(&entrenamiento));
    narracion.agregar(
        v.pizarra()
            .secuencia(bases, &marcas_de_genes(&entrenamiento))
            .lista(
                "Marcos largos del primer ejemplo",
                lista_de_genes(&entrenamiento),
                &[0],
            )
            .paso(
                "let mut entrenamiento = marcos_mas_largos(hebras, marcos);",
                format!(
                    "Para empezar a aprender hacen falta genes seguros. Un marco abierto de {LONGITUD_PARA_ENTRENAR} pb o más rara vez ocurre por azar, así que se toman {} como si fueran genes.",
                    plural(entrenamiento.len(), "marco largo", "marcos largos")
                ),
            ),
    );

    for ronda in 1..=RONDAS_DE_ENTRENAMIENTO {
        v.poner("ronda", ronda);
        narracion.agregar(
            v.pizarra().secuencia(bases, &inicio_de_vista).paso(
                "for ronda in 1..=RONDAS_DE_ENTRENAMIENTO {",
                if ronda == 1 {
                    format!(
                        "Ronda 1 de {RONDAS_DE_ENTRENAMIENTO}: se aprende con los marcos largos."
                    )
                } else {
                    format!(
                        "Ronda {ronda} de {RONDAS_DE_ENTRENAMIENTO}: ahora se aprende con los genes que el modelo anterior eligió, que son mejores ejemplos que los marcos largos."
                    )
                },
            ),
        );
        if ronda > 1 {
            entrenamiento = genes::genes_con_modelo(hebras, marcos, &modelo);
            v.poner("entrenamiento", resumen_de_ronda(&entrenamiento));
            narracion.agregar(
                v.pizarra()
                    .secuencia(bases, &marcas_de_genes(&entrenamiento))
                    .paso(
                        "entrenamiento = genes_con_modelo(hebras, marcos, &modelo);",
                        format!(
                            "Con lo aprendido hasta ahora se vuelven a predecir los genes de todo el ADN: salen {}.",
                            plural(entrenamiento.len(), "gen", "genes")
                        ),
                    ),
            );
            entrenamiento.retain(|predicho| {
                predicho.completo() && predicho.longitud() >= LONGITUD_PARA_ENTRENAR
            });
            v.poner("entrenamiento", resumen_de_ronda(&entrenamiento));
            narracion.agregar(
                v.pizarra()
                    .secuencia(bases, &marcas_de_genes(&entrenamiento))
                    .paso(
                        "entrenamiento.retain(|predicho| {",
                        format!(
                            "Solo se conservan los genes completos de {LONGITUD_PARA_ENTRENAR} pb o más, para no aprender de ejemplos dudosos: quedan {}.",
                            plural(entrenamiento.len(), "gen", "genes")
                        ),
                    ),
            );
        }
        let en_genes = hexameros_en_genes(hebras, &entrenamiento);
        let total_en_genes: f64 = en_genes.iter().sum();
        v.poner(
            "en_genes",
            format!("{} palabras contadas", miles(total_en_genes as usize)),
        );
        narracion.agregar(
            v.pizarra()
                .secuencia(bases, &marcas_de_genes(&entrenamiento))
                .paso(
                    "let en_genes = hexameros_en_genes(hebras, &entrenamiento);",
                    format!(
                        "Ahora se cuentan las palabras de 6 letras pero solo dentro de esos genes, leídas de tres en tres como codones: {} palabras. Si una palabra aparece más aquí que en el fondo, probablemente es típica de los genes.",
                        miles(total_en_genes as usize)
                    ),
                ),
        );
        modelo.codificacion = razon_logaritmica(&en_genes, &fondo);
        let mejores = mejores_hexameros(&modelo.codificacion);
        let filas: Vec<Vec<String>> = mejores
            .iter()
            .map(|&numero| {
                vec![
                    palabra_de_hexamero(numero),
                    miles(en_genes[numero] as usize),
                    miles(fondo[numero] as usize),
                    decimal(modelo.codificacion[numero], 2),
                ]
            })
            .collect();
        narracion.agregar(
            v.pizarra()
                .secuencia(bases, &inicio_de_vista)
                .tabla(
                    "Las palabras que más suenan a gen",
                    &["Palabra", "En genes", "En todo el ADN", "Puntaje"],
                    filas,
                    &[(0, 3)],
                )
                .paso(
                    "modelo.codificacion = razon_logaritmica(&en_genes, &fondo);",
                    "A cada palabra se le da un puntaje: el logaritmo de cuánto más (o menos) frecuente es dentro de los genes que en todo el ADN. Positivo significa «suena a gen»; negativo, «suena a ADN que no codifica». Aquí están las cinco mejores.",
                ),
        );
        genes::entrenar_inicios(&mut modelo, marcos, &entrenamiento);
        modelo.genes_de_entrenamiento = entrenamiento.len();
        let filas: Vec<Vec<String>> = CODONES_DE_INICIO
            .iter()
            .zip(modelo.codones_de_inicio.iter())
            .map(|(codon, puntaje)| vec![codon.to_string(), decimal(*puntaje, 2)])
            .collect();
        v.poner(
            "modelo",
            format!(
                "ronda {ronda}: {}",
                plural(modelo.genes_de_entrenamiento, "gen", "genes")
            ),
        );
        narracion.agregar(
            v.pizarra()
                .secuencia(bases, &inicio_de_vista)
                .tabla(
                    "Puntaje de cada codón de inicio",
                    &["Codón", "Puntaje"],
                    filas,
                    &[(0, 1)],
                )
                .paso(
                    "entrenar_inicios(&mut modelo, marcos, &entrenamiento);",
                    format!(
                        "Lo mismo se hace con los inicios: se compara qué codones y qué RBS usan los genes con los de todos los inicios posibles. Si ATG gana puntos, es el favorito de esta bacteria. El modelo aprendió de {}.",
                        plural(modelo.genes_de_entrenamiento, "gen", "genes")
                    ),
                ),
        );
    }
    modelo
}

#[cfg(test)]
mod pruebas {
    use super::*;

    fn comprobar(adn: &[u8]) -> Vec<Paso> {
        let (pasos, prediccion) = narrar_genes(adn, modelo_ecoli());
        let real = genes::predecir_genes(adn, modelo_ecoli());
        assert_eq!(prediccion, real);
        assert!(!pasos.is_empty() && pasos.len() <= MAXIMO_DE_PASOS);
        let ultimo = pasos.last().unwrap();
        let genes_vistos = ultimo.estado["variables"]
            .as_array()
            .unwrap()
            .iter()
            .find(|par| par[0] == "genes")
            .map(|par| par[1].as_str().unwrap().to_string())
            .unwrap();
        assert_eq!(genes_vistos, plural(real.genes.len(), "gen", "genes"));
        pasos
    }

    #[test]
    fn el_ejemplo_da_el_mismo_resultado_que_el_nucleo() {
        let entrada = Entrada::nueva(&[("adn", "@ecoli:1-150000")]);
        let adn = leer_adn(&entrada).unwrap();
        let pasos = comprobar(&adn);
        assert!(pasos.len() > 30);
        let traza = trazar_genes(&entrada).unwrap();
        assert_eq!(traza.len(), pasos.len());
    }

    #[test]
    fn una_region_corta_usa_el_modelo_de_respaldo() {
        let genoma = ecoli::genoma();
        let pasos = comprobar(&genoma[10_000..13_000]);
        assert!(
            pasos
                .iter()
                .any(|paso| paso.linea.contains("respaldo.clone()"))
        );
    }

    #[test]
    fn una_secuencia_sin_genes_no_rompe() {
        let adn = vec![b'A'; 200];
        comprobar(&adn);
        let mezcla: Vec<u8> = (0..300).map(|i| b"ACGTN"[(i * i + i / 7) % 5]).collect();
        comprobar(&mezcla);
    }
}
