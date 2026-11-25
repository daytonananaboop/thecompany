use crate::ejemplos::escribir_simulacion;
use crate::formato::{decimal, miles, porcentaje};
use crate::pizarra::{Marca, Narracion, Pizarra};
use crate::{Algoritmo, Bloque, Campo, Entrada, Paso, Salida, Tipo};
use nucleo::azar::Azar;
use nucleo::ecoli;
use nucleo::formatos::{Lectura, RegistroSam};
use nucleo::secuencia::complemento_inverso;
use nucleo::simulacion::{
    Mutacion, OpcionesSimulacion, Simulacion, agregar_errores, aplicar_mutaciones, calidades,
    cigar_desde_origen, elegir_region, mutaciones_al_azar, mutaciones_de_resistencia,
    regiones_en_mutado, simular_lecturas,
};

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![Algoritmo {
        id: "puerto.simular_lecturas",
        capitulo: "puerto",
        nombre: "Simular lecturas de secuenciación",
        resumen: "Crea un FASTQ realista de E. coli con mutaciones conocidas, para practicar con las demás herramientas.",
        campos: &[
            Campo {
                id: "region",
                etiqueta: "Región de E. coli (INICIO-FIN, desde 1)",
                tipo: Tipo::Texto,
                ejemplo: "3900001-4200000",
            },
            Campo {
                id: "longitud",
                etiqueta: "Longitud de cada lectura (pb)",
                tipo: Tipo::Entero,
                ejemplo: "150",
            },
            Campo {
                id: "cobertura",
                etiqueta: "Cobertura (veces que se lee cada base)",
                tipo: Tipo::Decimal,
                ejemplo: "15",
            },
            Campo {
                id: "snps",
                etiqueta: "SNPs al azar",
                tipo: Tipo::Entero,
                ejemplo: "20",
            },
            Campo {
                id: "indels",
                etiqueta: "Inserciones y deleciones al azar",
                tipo: Tipo::Entero,
                ejemplo: "4",
            },
            Campo {
                id: "resistencia",
                etiqueta: "Agregar mutaciones de resistencia (rpoB S531F, gyrA S83L)",
                tipo: Tipo::Opcion(&["sí", "no"]),
                ejemplo: "sí",
            },
            Campo {
                id: "adaptador",
                etiqueta: "Fracción de lecturas con adaptador",
                tipo: Tipo::Decimal,
                ejemplo: "0.08",
            },
            Campo {
                id: "semilla",
                etiqueta: "Semilla",
                tipo: Tipo::Entero,
                ejemplo: "1",
            },
        ],
        funcion: "simulacion::simular_lecturas",
        ejecutar: simular,
        pasos: Some(pasos),
    }]
}

fn leer_region(texto: &str) -> Result<(usize, usize), String> {
    let longitud = ecoli::genoma().len();
    let limpio = texto.trim();
    if limpio.eq_ignore_ascii_case("todo") {
        return Ok((0, longitud));
    }
    let (inicio, fin) = limpio.split_once('-').ok_or_else(|| {
        "Escribe la región como INICIO-FIN, por ejemplo 3900001-4200000, o «todo».".to_string()
    })?;
    let inicio: usize = inicio
        .trim()
        .parse()
        .map_err(|_| "El inicio de la región debe ser un número.".to_string())?;
    let fin: usize = fin
        .trim()
        .parse()
        .map_err(|_| "El fin de la región debe ser un número.".to_string())?;
    if inicio == 0 || inicio > fin || fin > longitud {
        return Err(format!(
            "La región debe estar entre 1 y {} y el inicio no puede pasar del fin.",
            miles(longitud)
        ));
    }
    if fin - inicio + 1 < 1000 {
        return Err("La región debe medir al menos 1000 pb.".to_string());
    }
    Ok((inicio - 1, fin))
}

struct Preparacion {
    inicio: usize,
    fin: usize,
    cobertura: f64,
    resistencia: bool,
    fraccion: f64,
    semilla: u64,
    fuera: usize,
    mutaciones: Vec<Mutacion>,
    opciones: OpcionesSimulacion,
}

fn preparar(entrada: &Entrada) -> Result<Preparacion, String> {
    let (inicio, fin) = leer_region(&entrada.texto("region")?)?;
    let longitud = entrada.natural("longitud")?;
    let cobertura = entrada.decimal("cobertura")?;
    let snps = entrada.natural("snps")?;
    let indels = entrada.natural("indels")?;
    let resistencia = entrada.opcion("resistencia", &["sí", "no"])? == "sí";
    let fraccion = entrada.decimal("adaptador")?;
    let semilla = entrada.natural("semilla")? as u64;
    if snps + indels > 1000 {
        return Err("Pide como máximo 1000 mutaciones al azar en total.".to_string());
    }
    if !(0.0..=0.9).contains(&fraccion) {
        return Err("La fracción con adaptador debe estar entre 0 y 0.9.".to_string());
    }
    if (fin - inicio) as f64 * cobertura > 150_000_000.0 {
        return Err(
            "La simulación pasaría de 150 millones de bases; baja la cobertura o achica la región."
                .to_string(),
        );
    }
    let genoma = ecoli::genoma();
    let mut conocidas = if resistencia {
        mutaciones_de_resistencia()?
    } else {
        Vec::new()
    };
    let fuera = conocidas
        .iter()
        .filter(|mutacion| mutacion.posicion < inicio || mutacion.posicion >= fin)
        .count();
    conocidas.retain(|mutacion| mutacion.posicion >= inicio && mutacion.posicion < fin);
    let ocupadas: Vec<usize> = conocidas.iter().map(|mutacion| mutacion.posicion).collect();
    let mut azar = Azar::nuevo(semilla ^ 0x5EED);
    let mut mutaciones =
        mutaciones_al_azar(genoma, &[(inicio, fin)], snps, indels, &ocupadas, &mut azar);
    mutaciones.extend(conocidas);
    mutaciones.sort_by_key(|mutacion| mutacion.posicion);
    let opciones = OpcionesSimulacion {
        longitud_lectura: longitud,
        cobertura,
        fraccion_con_adaptador: fraccion,
        regiones: vec![(inicio, fin)],
        semilla,
        ..OpcionesSimulacion::default()
    };
    Ok(Preparacion {
        inicio,
        fin,
        cobertura,
        resistencia,
        fraccion,
        semilla,
        fuera,
        mutaciones,
        opciones,
    })
}

fn simular(entrada: &Entrada) -> Result<Salida, String> {
    let Preparacion {
        inicio,
        fin,
        cobertura,
        resistencia,
        fraccion,
        semilla,
        fuera,
        mutaciones,
        opciones,
    } = preparar(entrada)?;
    let longitud = opciones.longitud_lectura;
    let genoma = ecoli::genoma();
    let muestra = simular_lecturas(genoma, ecoli::IDENTIFICADOR, &mutaciones, &opciones)?;
    let carpeta = crate::carpeta_de_ejemplos()?;
    let base = format!("simulacion_{}_{}_{}", inicio + 1, fin, semilla);
    let ruta_lecturas = carpeta.join(format!("{base}.fastq"));
    let ruta_alineamientos = carpeta.join(format!("{base}.verdad.sam"));
    let ruta_mutaciones = carpeta.join(format!("{base}.mutaciones.vcf"));
    escribir_simulacion(
        &muestra,
        ecoli::IDENTIFICADOR,
        genoma.len(),
        &ruta_lecturas,
        &ruta_alineamientos,
        &ruta_mutaciones,
    )?;

    let total_bases: usize = muestra.lecturas.iter().map(|l| l.secuencia.len()).sum();
    let mut calidad_por_posicion = vec![0.0; longitud];
    for lectura in &muestra.lecturas {
        for (posicion, q) in lectura.calidad.iter().enumerate() {
            calidad_por_posicion[posicion] += (q - 33) as f64;
        }
    }
    let lecturas = muestra.lecturas.len().max(1) as f64;
    let puntos: Vec<[f64; 2]> = calidad_por_posicion
        .iter()
        .enumerate()
        .map(|(posicion, suma)| [(posicion + 1) as f64, suma / lecturas])
        .collect();
    let calidad_media = puntos.iter().map(|p| p[1]).sum::<f64>() / puntos.len().max(1) as f64;
    let filas: Vec<Vec<String>> = muestra
        .mutaciones
        .iter()
        .map(|mutacion| {
            vec![
                miles(mutacion.posicion + 1),
                String::from_utf8_lossy(&mutacion.referencia).to_string(),
                String::from_utf8_lossy(&mutacion.alternativa).to_string(),
                mutacion.nombre.clone(),
            ]
        })
        .collect();

    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Lecturas", miles(muestra.lecturas.len())));
    salida.agregar(Bloque::valor("Bases secuenciadas", miles(total_bases)));
    salida.agregar(Bloque::archivo("Lecturas (FASTQ)", &ruta_lecturas));
    salida.agregar(Bloque::archivo(
        "Posición verdadera de cada lectura (SAM)",
        &ruta_alineamientos,
    ));
    salida.agregar(Bloque::archivo(
        "Mutaciones verdaderas (VCF)",
        &ruta_mutaciones,
    ));
    salida.agregar(Bloque::tabla(
        "Mutaciones sembradas",
        &["Posición (desde 1)", "Referencia", "Alternativa", "Tipo"],
        filas,
    ));
    salida.agregar(Bloque::linea(
        "Calidad media por posición de la lectura",
        "Posición en la lectura",
        "Calidad (Phred)",
        &puntos,
        Vec::new(),
    ));

    let sin_cubrir = (-cobertura).exp();
    if cobertura >= 10.0 {
        salida.explicar(format!(
            "Cobertura de {}×: cada base se lee unas {} veces. Para un genoma haploide como el de E. coli, 10× o más basta para detectar SNPs con confianza.",
            decimal(cobertura, 1),
            decimal(cobertura, 0)
        ));
    } else {
        salida.explicar(format!(
            "Cobertura de {}×. Según el modelo de Lander y Waterman, se espera que una fracción e^(−c) ≈ {} de la región no tenga ninguna lectura; con menos de 10× muchas variantes quedarán con muy pocas lecturas para confirmarlas.",
            decimal(cobertura, 1),
            porcentaje(sin_cubrir, 2)
        ));
    }
    salida.explicar(format!(
        "La calidad baja de Q{} al inicio de la lectura a Q{} al final, como en un secuenciador real. Q30 significa 1 error cada 1000 bases; Q20, 1 cada 100. Calidad media: Q{}.",
        decimal(puntos.first().map(|p| p[1]).unwrap_or(0.0), 0),
        decimal(puntos.last().map(|p| p[1]).unwrap_or(0.0), 0),
        decimal(calidad_media, 1)
    ));
    if fraccion > 0.0 {
        salida.explicar(format!(
            "Alrededor de {} de las lecturas vienen de fragmentos más cortos que la lectura, así que terminan en adaptador. «Calidad de lecturas» debería detectarlo y recortarlo.",
            porcentaje(fraccion, 0)
        ));
    }
    if resistencia && fuera == 0 {
        salida.explicar(
            "Se sembraron rpoB S531F (resistencia a rifampicina) y gyrA S83L (resistencia a quinolonas). Mapea las lecturas y llama variantes para ver si las encuentras.",
        );
    } else if resistencia {
        salida.explicar(format!(
            "{fuera} mutación(es) de resistencia quedan fuera de la región elegida y no se sembraron. rpoB está en 4 181 245–4 185 273 y gyrA en 2 336 793–2 339 420."
        ));
    }
    salida.explicar(
        "Siguiente paso: «Calidad de lecturas» con el FASTQ, luego «Mapeo de lecturas» y «Llamado de variantes». Compara tus variantes con el VCF verdadero.",
    );
    Ok(salida)
}

const LIMITE_NARRADO: f64 = 20_000_000.0;

#[derive(Default, Clone)]
struct Variables(Vec<(&'static str, String)>);

impl Variables {
    fn poner(&mut self, nombre: &'static str, valor: impl ToString) {
        let valor = valor.to_string();
        match self.0.iter_mut().find(|(clave, _)| *clave == nombre) {
            Some(existente) => existente.1 = valor,
            None => self.0.push((nombre, valor)),
        }
    }

    fn cortar(&mut self, cuantas: usize) {
        self.0.truncate(cuantas);
    }

    fn pizarra(&self) -> Pizarra {
        self.0
            .iter()
            .fold(Pizarra::nueva(), |pizarra, (nombre, valor)| {
                pizarra.variable(nombre, valor)
            })
    }
}

fn corto(texto: &[u8], maximo: usize) -> String {
    let visible = String::from_utf8_lossy(&texto[..texto.len().min(maximo)]).to_string();
    if texto.len() > maximo {
        format!("{visible}…")
    } else {
        visible
    }
}

fn numeros_cortos(valores: &[u8], maximo: usize) -> String {
    let mut texto = valores
        .iter()
        .take(maximo)
        .map(|valor| valor.to_string())
        .collect::<Vec<_>>()
        .join(", ");
    if valores.len() > maximo {
        texto.push_str(", …");
    }
    texto
}

fn texto_cigar(cigar: &[(u32, u8)]) -> String {
    cigar
        .iter()
        .map(|(cuantas, operacion)| format!("{cuantas}{}", *operacion as char))
        .collect()
}

fn indices_mutados(ordenadas: &[Mutacion]) -> Vec<usize> {
    let mut desfase: isize = 0;
    let mut indices = Vec::with_capacity(ordenadas.len());
    for mutacion in ordenadas {
        indices.push((mutacion.posicion as isize + desfase) as usize);
        desfase += mutacion.alternativa.len() as isize - mutacion.referencia.len() as isize;
    }
    indices
}

fn describir_mutacion(mutacion: &Mutacion) -> String {
    format!(
        "{} → {}",
        String::from_utf8_lossy(&mutacion.referencia),
        String::from_utf8_lossy(&mutacion.alternativa)
    )
}

fn cuenta(filas: Vec<(&str, usize)>) -> Vec<Vec<String>> {
    filas
        .iter()
        .map(|(nombre, valor)| vec![nombre.to_string(), miles(*valor)])
        .collect()
}

fn pasos(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    trazar(entrada).map(|(pasos, _)| pasos)
}

fn trazar(entrada: &Entrada) -> Result<(Vec<Paso>, Simulacion), String> {
    let preparacion = preparar(entrada)?;
    if (preparacion.fin - preparacion.inicio) as f64 * preparacion.cobertura > LIMITE_NARRADO {
        return Err(
            "Para contarlo paso a paso, pide como máximo 20 millones de bases (región × cobertura); baja la cobertura o achica la región."
                .to_string(),
        );
    }
    let referencia = ecoli::genoma();
    let opciones = &preparacion.opciones;
    let mutaciones = &preparacion.mutaciones;
    let mut narracion = Narracion::nueva(400);
    let mut v = Variables::default();

    let snps_dados = mutaciones
        .iter()
        .filter(|m| m.referencia.len() == 1 && m.alternativa.len() == 1)
        .count();
    let filas_mutaciones: Vec<Vec<String>> = mutaciones
        .iter()
        .take(8)
        .map(|m| {
            vec![
                miles(m.posicion + 1),
                String::from_utf8_lossy(&m.referencia).to_string(),
                String::from_utf8_lossy(&m.alternativa).to_string(),
                m.nombre.clone(),
            ]
        })
        .collect();
    v.poner("mutaciones", mutaciones.len());
    v.poner("opciones.longitud_lectura", opciones.longitud_lectura);
    v.poner("opciones.cobertura", decimal(opciones.cobertura, 1));
    v.poner(
        "opciones.fraccion_con_adaptador",
        decimal(opciones.fraccion_con_adaptador, 2),
    );
    v.poner("opciones.semilla", opciones.semilla);
    narracion.agregar(
        v.pizarra()
            .tabla(
                "Mutaciones recibidas (primeras 8)",
                &["Posición", "Referencia", "Alternativa", "Tipo"],
                filas_mutaciones,
                &[],
            )
            .paso(
                "pub fn simular_lecturas(",
                format!(
                    "Vas a fabricar un FASTQ. Parto de la referencia de E. coli ({} pb) y de {} mutaciones ya elegidas ({} cambios de una base y {} inserciones o deleciones). Voy a «secuenciar» fragmentos al azar hasta cubrir la región {}× con lecturas de {} bases.",
                    miles(referencia.len()),
                    mutaciones.len(),
                    snps_dados,
                    mutaciones.len() - snps_dados,
                    decimal(opciones.cobertura, 1),
                    opciones.longitud_lectura
                ),
            ),
    );

    let longitud = opciones.longitud_lectura;
    v.poner("longitud", longitud);
    narracion.agregar(v.pizarra().paso(
        "let longitud = opciones.longitud_lectura;",
        format!("Guardo el largo de lectura en una variable más corta: longitud = {longitud}."),
    ));
    if !(30..=1000).contains(&longitud) {
        return Err("La longitud de lectura debe estar entre 30 y 1000 bases.".to_string());
    }
    narracion.agregar(v.pizarra().paso(
        "if !(30..=1000).contains(&longitud) {",
        format!("Compruebo que {longitud} esté entre 30 y 1000. Lo está, así que sigo; si no, el programa se detendría con un aviso."),
    ));
    if !(0.1..=200.0).contains(&opciones.cobertura) {
        return Err("La cobertura debe estar entre 0.1 y 200.".to_string());
    }
    narracion.agregar(v.pizarra().paso(
        "if !(0.1..=200.0).contains(&opciones.cobertura) {",
        format!(
            "Lo mismo con la cobertura: {}× está entre 0.1 y 200.",
            decimal(opciones.cobertura, 1)
        ),
    ));
    if referencia.len() < longitud * 2 {
        return Err(format!(
            "La referencia debe medir al menos {} bases para lecturas de {longitud}.",
            longitud * 2
        ));
    }
    narracion.agregar(v.pizarra().paso(
        "if referencia.len() < longitud * 2 {",
        format!(
            "La referencia ({} pb) es mucho más larga que dos lecturas ({}), así que alcanza.",
            miles(referencia.len()),
            miles(longitud * 2)
        ),
    ));

    let mut ordenadas = mutaciones.to_vec();
    ordenadas.sort_by_key(|mutacion| mutacion.posicion);
    for par in ordenadas.windows(2) {
        if par[0].posicion + par[0].referencia.len() > par[1].posicion {
            return Err("Hay mutaciones que se superponen.".to_string());
        }
    }
    let lista_ordenadas: Vec<String> = ordenadas
        .iter()
        .take(12)
        .map(|m| format!("{}: {}", miles(m.posicion + 1), describir_mutacion(m)))
        .collect();
    v.poner("ordenadas", format!("{} mutaciones", ordenadas.len()));
    narracion.agregar(
        v.pizarra()
            .lista("ordenadas (primeras 12)", lista_ordenadas, &[])
            .paso(
                "ordenadas.sort_by_key(|mutacion| mutacion.posicion);",
                "Pongo las mutaciones en orden de posición a lo largo del genoma, porque las voy a aplicar de izquierda a derecha.",
            ),
    );
    narracion.agregar(v.pizarra().paso(
        "for par in ordenadas.windows(2) {",
        "Reviso cada pareja de vecinas: ninguna mutación debe empezar antes de que termine la anterior. Aquí no se pisan.",
    ));

    let genoma = aplicar_mutaciones(referencia, &ordenadas);
    let indices = indices_mutados(&ordenadas);
    v.poner("genoma.bases.len()", miles(genoma.bases.len()));
    narracion.agregar(
        v.pizarra().paso(
            "let genoma = aplicar_mutaciones(referencia, &ordenadas);",
            format!(
                "Construyo el genoma de la muestra: una copia de la referencia con las mutaciones puestas. Mide {} bases ({} {} que la referencia). Además guardo en origen, base por base, de qué posición de la referencia viene (o nada si es una inserción).",
                miles(genoma.bases.len()),
                miles(genoma.bases.len().abs_diff(referencia.len())),
                if genoma.bases.len() >= referencia.len() {
                    "más"
                } else {
                    "menos"
                }
            ),
        ),
    );
    let primera_indel = ordenadas
        .iter()
        .position(|m| m.referencia.len() != m.alternativa.len());
    let mut ejemplos: Vec<usize> = Vec::new();
    if !ordenadas.is_empty() {
        ejemplos.push(0);
    }
    if let Some(indice) = primera_indel
        && !ejemplos.contains(&indice)
    {
        ejemplos.push(indice);
    }
    for indice in ejemplos {
        let mutacion = &ordenadas[indice];
        let en_mutado = indices[indice];
        v.poner("posicion", mutacion.posicion);
        narracion.agregar(
            v.pizarra()
                .secuencia(
                    referencia,
                    &[(
                        mutacion.posicion,
                        mutacion.posicion + mutacion.referencia.len(),
                        Marca::Actual,
                    )],
                )
                .paso(
                    "while posicion < mutacion.posicion {",
                    format!(
                        "Copio la referencia sin cambios hasta llegar a la posición {} ({}).",
                        miles(mutacion.posicion + 1),
                        mutacion.nombre
                    ),
                ),
        );
        let comunes = mutacion.referencia.len().min(mutacion.alternativa.len());
        v.poner("comunes", comunes);
        narracion.agregar(
            v.pizarra()
                .secuencia(
                    &genoma.bases,
                    &[(
                        en_mutado,
                        en_mutado + mutacion.alternativa.len(),
                        Marca::Hallado,
                    )],
                )
                .paso(
                    "bases.push(mutacion.alternativa[indice]);",
                    format!(
                        "Aquí la muestra cambia: donde la referencia tenía {} la muestra tiene {}. Esa base marcada es la que luego tendrán que encontrar las herramientas de variantes.",
                        String::from_utf8_lossy(&mutacion.referencia),
                        String::from_utf8_lossy(&mutacion.alternativa)
                    ),
                ),
        );
        if mutacion.alternativa.len() > comunes {
            narracion.agregar(v.pizarra().paso(
                "for &insertada in &mutacion.alternativa[comunes..] {",
                format!(
                    "Es una inserción: sobran {} base(s) que no existen en la referencia, y en origen quedan sin posición.",
                    mutacion.alternativa.len() - comunes
                ),
            ));
        } else if mutacion.referencia.len() > comunes {
            narracion.agregar(v.pizarra().paso(
                "posicion += mutacion.referencia.len();",
                format!(
                    "Es una deleción: salto {} base(s) de la referencia, así que la muestra queda más corta.",
                    mutacion.referencia.len() - comunes
                ),
            ));
        }
    }
    v.poner("posicion", referencia.len());
    narracion.agregar(v.pizarra().paso(
        "while posicion < referencia.len() {",
        "Después de la última mutación copio el resto de la referencia sin cambios. El genoma de la muestra está listo.",
    ));

    let regiones = regiones_en_mutado(&genoma, &opciones.regiones, referencia.len());
    let largo_total: usize = regiones.iter().map(|(inicio, fin)| fin - inicio).sum();
    let texto_regiones = regiones
        .iter()
        .map(|(inicio, fin)| format!("[{}, {})", miles(*inicio), miles(*fin)))
        .collect::<Vec<_>>()
        .join(" ");
    v.poner("regiones", texto_regiones);
    narracion.agregar(v.pizarra().paso(
        "let regiones = regiones_en_mutado(&genoma, &opciones.regiones, referencia.len());",
        "La región que pediste estaba medida en la referencia. Como las inserciones y deleciones mueven las posiciones, la traduzco a posiciones del genoma de la muestra.",
    ));
    if largo_total < longitud {
        return Err("Las regiones son más cortas que una lectura.".to_string());
    }
    v.poner("largo_total", miles(largo_total));
    narracion.agregar(v.pizarra().paso(
        "let largo_total: usize",
        format!(
            "Sumo el largo de las regiones: {} bases en total. Es el territorio que se va a secuenciar.",
            miles(largo_total)
        ),
    ));
    let cantidad = ((largo_total as f64 * opciones.cobertura) / longitud as f64).round() as usize;
    v.poner("cantidad", miles(cantidad));
    narracion.agregar(v.pizarra().paso(
        "let cantidad = ",
        format!(
            "Cuántas lecturas hacen falta: largo × cobertura ÷ longitud = {} × {} ÷ {} ≈ {} lecturas. Así cada base queda leída unas {} veces.",
            miles(largo_total),
            decimal(opciones.cobertura, 1),
            longitud,
            miles(cantidad),
            decimal(opciones.cobertura, 0)
        ),
    ));
    if cantidad * longitud > 150_000_000 {
        return Err(
            "La simulación pasaría de 150 millones de bases; baja la cobertura.".to_string(),
        );
    }
    let mut azar = Azar::nuevo(opciones.semilla);
    narracion.agregar(v.pizarra().paso(
        "let mut azar = Azar::nuevo(opciones.semilla);",
        format!(
            "Preparo el generador de azar con la semilla {}. Es azar «de mentira» repetible: con la misma semilla saldrán exactamente las mismas lecturas.",
            opciones.semilla
        ),
    ));
    let mut lecturas = Vec::with_capacity(cantidad);
    let mut verdad = Vec::with_capacity(cantidad);
    narracion.agregar(v.pizarra().paso(
        "let mut lecturas = Vec::with_capacity(cantidad);",
        "Preparo dos listas vacías: lecturas (lo que escribiría el secuenciador) y verdad (dónde estaba realmente cada lectura).",
    ));

    let base = v.0.len();
    let mut con_adaptador_total = 0usize;
    let mut reversas = 0usize;
    let mut bases_cambiadas = 0usize;
    let mut cubren = 0usize;
    let mut omitidas = 0usize;
    let mut visto_adaptador = false;
    let mut visto_error = false;
    let mut visto_cubre = false;
    for numero in 0..cantidad {
        let sorteo = azar.decimal();
        let con_adaptador = sorteo < opciones.fraccion_con_adaptador;
        let nivel = if numero == 0 || (con_adaptador && !visto_adaptador) {
            2
        } else if numero == 1 {
            1
        } else {
            0
        };
        if con_adaptador {
            visto_adaptador = true;
            con_adaptador_total += 1;
        }
        if nivel >= 1 {
            v.cortar(base);
            v.poner("numero", numero);
            if nivel == 2 {
                narracion.agregar(v.pizarra().paso(
                    "for numero in 0..cantidad {",
                    format!(
                        "Empieza la lectura número {} de {}. Cada vuelta del bucle fabrica una lectura completa.",
                        numero + 1,
                        miles(cantidad)
                    ),
                ));
                v.poner("con_adaptador", con_adaptador);
                narracion.agregar(v.pizarra().paso(
                    "let con_adaptador = azar.decimal() < opciones.fraccion_con_adaptador;",
                    if con_adaptador {
                        format!(
                            "Saco un número al azar entre 0 y 1: {}. Es menor que {}, así que este fragmento será corto y la lectura terminará en adaptador.",
                            decimal(sorteo, 3),
                            decimal(opciones.fraccion_con_adaptador, 2)
                        )
                    } else {
                        format!(
                            "Saco un número al azar entre 0 y 1: {}. No es menor que {}, así que este fragmento será de largo completo.",
                            decimal(sorteo, 3),
                            decimal(opciones.fraccion_con_adaptador, 2)
                        )
                    },
                ));
            }
        }
        let largo_fragmento = if con_adaptador {
            longitud / 3 + azar.entero(longitud - longitud / 3 - 10)
        } else {
            longitud
        };
        if nivel == 2 {
            v.poner("largo_fragmento", largo_fragmento);
            narracion.agregar(v.pizarra().paso(
                "let largo_fragmento = if con_adaptador {",
                if con_adaptador {
                    format!(
                        "El fragmento mide {largo_fragmento} bases, menos que la lectura ({longitud}). Faltarán {} bases y el secuenciador las llenará con adaptador.",
                        longitud - largo_fragmento
                    )
                } else {
                    format!("El fragmento mide {largo_fragmento}, justo lo que mide una lectura.")
                },
            ));
        }
        let (inicio_region, fin_region) = elegir_region(&regiones, &mut azar);
        if nivel == 2 {
            v.poner("inicio_region", miles(inicio_region));
            v.poner("fin_region", miles(fin_region));
            narracion.agregar(v.pizarra().paso(
                "let (inicio_region, fin_region) = elegir_region(&regiones, &mut azar);",
                "Elijo en qué región cae el fragmento (si hubiera varias, las más largas tienen más probabilidad). Aquí solo hay una.",
            ));
        }
        if fin_region - inicio_region < largo_fragmento {
            omitidas += 1;
            continue;
        }
        let inicio = inicio_region + azar.entero(fin_region - inicio_region - largo_fragmento + 1);
        let fragmento = &genoma.bases[inicio..inicio + largo_fragmento];
        let origen = &genoma.origen[inicio..inicio + largo_fragmento];
        if nivel >= 1 {
            v.poner("inicio", miles(inicio));
            narracion.agregar(
                v.pizarra()
                    .secuencia(
                        &genoma.bases,
                        &[(inicio, inicio + largo_fragmento, Marca::Actual)],
                    )
                    .paso(
                        "let inicio = inicio_region + azar.entero(",
                        format!(
                            "Elijo al azar dónde empieza el fragmento: posición {}. Todas las posiciones de la región son igual de probables, y por eso cada base termina cubierta por varias lecturas.",
                            miles(inicio + 1)
                        ),
                    ),
            );
        }
        let reversa = azar.entero(2) == 1;
        if reversa {
            reversas += 1;
        }
        if nivel >= 1 {
            v.poner("reversa", reversa);
            narracion.agregar(v.pizarra().paso(
                "let reversa = azar.entero(2) == 1;",
                if reversa {
                    "Lanzo una moneda: salió «reversa». Esta lectura viene de la hebra opuesta del ADN."
                } else {
                    "Lanzo una moneda: salió «directa». Esta lectura viene de la hebra de la referencia."
                },
            ));
        }
        let mut bases = if reversa {
            complemento_inverso(fragmento)
        } else {
            fragmento.to_vec()
        };
        let mut origen_lectura: Vec<Option<usize>> = if reversa {
            origen.iter().rev().copied().collect()
        } else {
            origen.to_vec()
        };
        if nivel >= 1 {
            v.poner("bases", corto(&bases, 24));
            narracion.agregar(
                v.pizarra()
                    .secuencia(&bases, &[(0, bases.len(), Marca::Actual)])
                    .paso(
                        "let mut bases = if reversa {",
                        if reversa {
                            "Como es de la hebra opuesta, tomo el complemento inverso del fragmento: se lee al revés y cambiando A↔T, C↔G."
                        } else {
                            "Como es directa, las bases son una copia del fragmento tal cual."
                        },
                    ),
            );
        }
        let mut relleno = 0;
        while bases.len() < longitud {
            let siguiente = if relleno < opciones.adaptador.len() {
                opciones.adaptador[relleno]
            } else {
                b"ACGT"[azar.entero(4)]
            };
            bases.push(siguiente);
            origen_lectura.push(None);
            relleno += 1;
        }
        if nivel == 2 && relleno > 0 {
            v.poner("relleno", relleno);
            narracion.agregar(
                v.pizarra()
                    .secuencia(
                        &bases,
                        &[
                            (0, largo_fragmento, Marca::Actual),
                            (largo_fragmento, bases.len(), Marca::Comparado),
                        ],
                    )
                    .paso(
                        "while bases.len() < longitud {",
                        format!(
                            "El fragmento era corto, así que completo la lectura con {relleno} bases: primero el adaptador (ACGGTAGC…) y, si sobra espacio, bases al azar. Esa parte no pertenece al genoma."
                        ),
                    ),
            );
        }
        let calidad_numerica = calidades(opciones, longitud, &mut azar);
        if nivel == 2 {
            v.poner("calidad_numerica", numeros_cortos(&calidad_numerica, 8));
            let media = calidad_numerica.iter().map(|&q| f64::from(q)).sum::<f64>()
                / calidad_numerica.len() as f64;
            narracion.agregar(v.pizarra().paso(
                "let valor = (media + azar.normal() * 2.5).round().clamp(2.0, 41.0);",
                format!(
                    "Cada base recibe una calidad Q. Empieza cerca de Q{} y baja hasta cerca de Q{} al final de la lectura, con algo de ruido. En esta lectura va de Q{} a Q{} (media Q{}).",
                    decimal(opciones.calidad_inicial, 0),
                    decimal(opciones.calidad_final, 0),
                    calidad_numerica[0],
                    calidad_numerica[longitud - 1],
                    decimal(media, 1)
                ),
            ));
            let q = calidad_numerica[0];
            let probabilidad_error = 10f64.powf(-f64::from(q) / 10.0);
            v.poner("q", q);
            v.poner("probabilidad_error", format!("{probabilidad_error:.5}"));
            narracion.agregar(v.pizarra().paso(
                "let probabilidad_error = 10f64.powf(-(q as f64) / 10.0);",
                format!(
                    "La calidad se vuelve probabilidad: la primera base tiene Q{q}, o sea {} de error, 1 equivocación cada {} bases.",
                    porcentaje(probabilidad_error, 3),
                    miles((1.0 / probabilidad_error).round() as usize)
                ),
            ));
        }
        let antes = bases.clone();
        agregar_errores(&mut bases, &calidad_numerica, &mut azar);
        let cambiadas: Vec<usize> = bases
            .iter()
            .zip(&antes)
            .enumerate()
            .filter(|(_, (nueva, vieja))| nueva != vieja)
            .map(|(posicion, _)| posicion)
            .collect();
        bases_cambiadas += cambiadas.len();
        let marcas_errores: Vec<(usize, usize, Marca)> = cambiadas
            .iter()
            .take(20)
            .map(|&posicion| (posicion, posicion + 1, Marca::Actual))
            .collect();
        if nivel == 2 || (nivel == 0 && !visto_error && !cambiadas.is_empty()) {
            let mut w = v.clone();
            if nivel == 0 {
                w.cortar(base);
                w.poner("numero", numero);
                w.poner("con_adaptador", con_adaptador);
                w.poner("reversa", reversa);
            }
            w.poner("bases", corto(&bases, 24));
            let nota = if cambiadas.is_empty() {
                "Para cada base saco un número al azar y, si es menor que su probabilidad de error, la cambio por otra. En esta lectura ninguna base se cambió: las calidades altas casi nunca fallan.".to_string()
            } else {
                format!(
                    "Para cada base saco un número al azar y, si es menor que su probabilidad de error, la cambio por otra. Aquí cambiaron {} base(s), marcadas abajo (primera: posición {} de la lectura). Son errores de secuenciación, no mutaciones de verdad.",
                    cambiadas.len(),
                    cambiadas[0] + 1
                )
            };
            narracion.agregar(
                w.pizarra()
                    .secuencia(&bases, &marcas_errores)
                    .paso("if azar.decimal() < probabilidad_error {", nota),
            );
        }
        if !cambiadas.is_empty() {
            visto_error = true;
        }
        let calidad: Vec<u8> = calidad_numerica.iter().map(|q| q + 33).collect();
        if nivel == 2 {
            v.poner("bases", corto(&bases, 24));
            v.poner("calidad_numerica", numeros_cortos(&calidad_numerica, 8));
            narracion.agregar(v.pizarra().paso(
                "let calidad: Vec<u8> = calidad_numerica.iter().map(|q| q + 33).collect();",
                format!(
                    "El FASTQ guarda cada calidad como una letra: sumo 33 a Q y uso el carácter con ese código. Las primeras: {}.",
                    corto(&calidad, 16)
                ),
            ));
        }
        let nombre = format!(
            "lectura{}_{}_{}",
            numero + 1,
            origen.iter().flatten().next().map(|p| p + 1).unwrap_or(0),
            if reversa { "menos" } else { "mas" }
        )
        .into_bytes();
        if nivel == 2 {
            v.poner("nombre", String::from_utf8_lossy(&nombre));
            narracion.agregar(v.pizarra().paso(
                "let nombre = format!(",
                "El nombre de la lectura dice su número, la posición de referencia donde empieza su fragmento y la hebra. Así luego puedes comprobar si el mapeo la puso en su sitio.",
            ));
        }
        let (secuencia_sam, calidad_sam, origen_sam) = if reversa {
            (
                complemento_inverso(&bases),
                calidad.iter().rev().copied().collect(),
                origen_lectura.iter().rev().copied().collect::<Vec<_>>(),
            )
        } else {
            (bases.clone(), calidad.clone(), origen_lectura.clone())
        };
        let (posicion, cigar) = cigar_desde_origen(&origen_sam);
        if nivel == 2 {
            v.poner("posicion", miles(posicion + 1));
            v.poner("cigar", texto_cigar(&cigar));
            narracion.agregar(v.pizarra().paso(
                "let (posicion, cigar) = cigar_desde_origen(&origen_sam);",
                format!(
                    "Calculo la respuesta correcta para el archivo SAM: la lectura empieza en la posición {} de la referencia y su CIGAR es {}. M son bases alineadas, S son bases que sobran (adaptador) e I/D son inserciones y deleciones.{}",
                    miles(posicion + 1),
                    texto_cigar(&cigar),
                    if reversa {
                        " Como la lectura es de la hebra opuesta, el SAM la guarda ya vuelta a invertir."
                    } else {
                        ""
                    }
                ),
            ));
        }
        verdad.push(RegistroSam {
            nombre: nombre.clone(),
            bandera: if reversa { 16 } else { 0 },
            referencia: ecoli::IDENTIFICADOR.to_string(),
            posicion: posicion + 1,
            calidad_mapeo: 60,
            cigar,
            secuencia: secuencia_sam,
            calidad: calidad_sam,
        });
        let primera_mutacion = indices.partition_point(|&indice| indice < inicio);
        let cubre = primera_mutacion < indices.len()
            && indices[primera_mutacion] < inicio + largo_fragmento;
        if cubre {
            cubren += 1;
        }
        if nivel == 2 {
            let registro = vec![
                format!("@{}", String::from_utf8_lossy(&nombre)),
                corto(&bases, 40),
                "+".to_string(),
                corto(&calidad, 40),
            ];
            narracion.agregar(
                v.pizarra()
                    .lista("registro FASTQ de la lectura", registro, &[0])
                    .paso(
                        "lecturas.push(Lectura {",
                        "Guardo la lectura: nombre, bases y calidades, las cuatro líneas de un registro FASTQ. La respuesta correcta queda aparte, en la lista verdad.",
                    ),
            );
        } else if nivel == 1 {
            narracion.agregar(v.pizarra().paso(
                "lecturas.push(Lectura {",
                "La segunda lectura sigue exactamente los mismos pasos, con otros números al azar. La guardo y paso a la siguiente.",
            ));
        }
        if cubre && !visto_cubre {
            visto_cubre = true;
            let mutacion = &ordenadas[primera_mutacion];
            let en_mutado = indices[primera_mutacion];
            let mut w = v.clone();
            w.cortar(base);
            w.poner("numero", numero);
            w.poner("inicio", miles(inicio));
            narracion.agregar(
                w.pizarra()
                    .secuencia(
                        &genoma.bases,
                        &[
                            (inicio, inicio + largo_fragmento, Marca::Actual),
                            (
                                en_mutado,
                                en_mutado + mutacion.alternativa.len(),
                                Marca::Hallado,
                            ),
                        ],
                    )
                    .paso(
                        "let fragmento = &genoma.bases[inicio..inicio + largo_fragmento];",
                        format!(
                            "La lectura {} es la primera que cae sobre una mutación: {} en la posición {} ({}). Su fragmento la incluye, así que esta lectura «dice» la base nueva. Cuando muchas lecturas coinciden en ella, «Llamado de variantes» la reconoce.",
                            numero + 1,
                            mutacion.nombre,
                            miles(mutacion.posicion + 1),
                            describir_mutacion(mutacion)
                        ),
                    ),
            );
        }
        lecturas.push(Lectura {
            nombre,
            secuencia: bases,
            calidad,
        });
        let siguiente_marca = (numero + 1) * 8 / cantidad;
        if nivel < 2 && numero + 1 < cantidad && siguiente_marca > numero * 8 / cantidad {
            let tabla = cuenta(vec![
                ("Lecturas guardadas", lecturas.len()),
                ("Con adaptador", con_adaptador_total),
                ("De la hebra opuesta", reversas),
                ("Bases cambiadas por error", bases_cambiadas),
                ("Lecturas sobre una mutación", cubren),
            ]);
            let mut w = v.clone();
            w.cortar(base);
            w.poner("numero", numero);
            narracion.agregar(
                w.pizarra()
                    .tabla("Cuenta hasta ahora", &["Qué", "Cuántas"], tabla, &[])
                    .paso(
                        "lecturas.push(Lectura {",
                        format!(
                            "Las lecturas siguientes repiten la misma receta con otros números al azar. Llevo {} de {} ({}).",
                            miles(numero + 1),
                            miles(cantidad),
                            porcentaje((numero + 1) as f64 / cantidad as f64, 0)
                        ),
                    ),
            );
        }
    }
    let total_bases: usize = lecturas.iter().map(|l| l.secuencia.len()).sum();
    let tabla = cuenta(vec![
        ("Lecturas", lecturas.len()),
        ("Bases secuenciadas", total_bases),
        ("Con adaptador", con_adaptador_total),
        ("De la hebra opuesta", reversas),
        ("Bases cambiadas por error", bases_cambiadas),
        ("Lecturas sobre una mutación", cubren),
        ("Fragmentos descartados", omitidas),
    ]);
    let mut fin = Variables::default();
    fin.poner("lecturas.len()", miles(lecturas.len()));
    fin.poner("verdad.len()", miles(verdad.len()));
    fin.poner("ordenadas.len()", ordenadas.len());
    let ultimo = fin
        .pizarra()
        .tabla("Resultado", &["Qué", "Cuántas"], tabla, &[])
        .paso(
            "Ok(Simulacion {",
            format!(
                "Listo: {} lecturas ({} bases), cada una con su respuesta correcta en verdad, y las {} mutaciones verdaderas. {} lecturas cayeron sobre una mutación y {} bases se cambiaron solo por errores de secuenciación. Mira «Qué significa» para interpretar la cobertura y la calidad; el botón Ejecutar guarda estos mismos datos en archivos FASTQ, SAM y VCF.",
                miles(lecturas.len()),
                miles(total_bases),
                ordenadas.len(),
                miles(cubren),
                miles(bases_cambiadas)
            ),
        );
    let pasos = narracion.terminar(ultimo);
    Ok((
        pasos,
        Simulacion {
            lecturas,
            verdad,
            mutaciones: ordenadas,
        },
    ))
}

#[cfg(test)]
mod pruebas {
    use super::*;

    fn entrada_de_ejemplo() -> Entrada {
        Entrada::desde_json(crate::ejemplo("puerto.simular_lecturas").unwrap()).unwrap()
    }

    #[test]
    fn el_paso_a_paso_da_el_mismo_resultado_que_el_algoritmo() {
        let entrada = entrada_de_ejemplo();
        let preparacion = preparar(&entrada).unwrap();
        let real = simular_lecturas(
            ecoli::genoma(),
            ecoli::IDENTIFICADOR,
            &preparacion.mutaciones,
            &preparacion.opciones,
        )
        .unwrap();
        let (pasos, narrada) = trazar(&entrada).unwrap();
        assert!(!pasos.is_empty() && pasos.len() <= 400);
        assert_eq!(real.lecturas.len(), narrada.lecturas.len());
        assert_eq!(real.verdad.len(), narrada.verdad.len());
        assert_eq!(real.mutaciones, narrada.mutaciones);
        for (a, b) in real.lecturas.iter().zip(&narrada.lecturas) {
            assert_eq!(a, b);
        }
        for (a, b) in real.verdad.iter().zip(&narrada.verdad) {
            assert_eq!(a, b);
        }
        println!("pasos del ejemplo: {}", pasos.len());
    }

    #[test]
    fn los_indices_mutados_apuntan_a_la_base_nueva() {
        let preparacion = preparar(&entrada_de_ejemplo()).unwrap();
        let mut ordenadas = preparacion.mutaciones.clone();
        ordenadas.sort_by_key(|m| m.posicion);
        let genoma = aplicar_mutaciones(ecoli::genoma(), &ordenadas);
        for (mutacion, indice) in ordenadas.iter().zip(indices_mutados(&ordenadas)) {
            assert_eq!(genoma.bases[indice], mutacion.alternativa[0]);
        }
    }

    #[test]
    fn rechaza_lo_que_no_se_puede_narrar() {
        let entrada = Entrada::nueva(&[
            ("region", "todo"),
            ("longitud", "150"),
            ("cobertura", "15"),
            ("snps", "0"),
            ("indels", "0"),
            ("resistencia", "no"),
            ("adaptador", "0"),
            ("semilla", "1"),
        ]);
        assert!(trazar(&entrada).is_err());
    }
}
