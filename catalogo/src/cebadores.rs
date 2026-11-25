use crate::formato::{decimal, miles, porcentaje};
use crate::pizarra::{Marca as Resalte, Narracion, Pizarra};
use crate::{Algoritmo, Bloque, Campo, Entrada, FilaAlineada, Marca, Paso, Salida, Tipo};
use nucleo::cebadores::{
    self, Apareamiento, BASES_DE_LA_PINZA, BASES_DEL_EXTREMO_3, Cebador, Condiciones,
    CorreccionSal, Diseno, Horquilla, OpcionesDiseno, PINZA_MAXIMA, PINZA_MINIMA, PLANTILLA_MAXIMA,
    PRODUCTO_MAXIMO, ParDeCebadores, UMBRAL_DIMERO, UMBRAL_DIMERO_EXTREMO_3, UMBRAL_HORQUILLA,
};
use nucleo::ecoli;
use nucleo::secuencia::complemento_inverso;

const CORRECCIONES: &[&str] = &["SantaLucia 1998 (Primer3)", "Owczarzy 2008"];
const TM_BUENA_MINIMA: f64 = 55.0;
const TM_BUENA_MAXIMA: f64 = 65.0;
const DIFERENCIA_IDEAL: f64 = 2.0;
const DIFERENCIA_PROBLEMATICA: f64 = 5.0;
const GC_MINIMO: f64 = 0.40;
const GC_MAXIMO: f64 = 0.60;
const QPCR_MINIMO: usize = 70;
const QPCR_MAXIMO: usize = 200;
const PRODUCTO_MINIMO_PERMITIDO: usize = 40;
const MARGEN_DE_TM: f64 = 3.0;
const VENTANA_DEL_PAISAJE: usize = 20;
const MEJORA_MINIMA_DEL_MAGNESIO: f64 = 0.1;
const LONGITUD_MINIMA_DEL_OLIGO: usize = 8;

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![
        Algoritmo {
            id: "puerto.disenar_cebadores",
            capitulo: "puerto",
            nombre: "Diseño de cebadores para PCR",
            resumen: "Busca los mejores pares de cebadores para amplificar una región, con las reglas de Primer3: Tm, GC, pinza GC, dímeros y especificidad.",
            campos: &[
                Campo {
                    id: "plantilla",
                    etiqueta: "Plantilla de ADN (hasta 20 000 pb)",
                    tipo: Tipo::Secuencia,
                    ejemplo: "@gen:lacZ",
                },
                Campo {
                    id: "producto",
                    etiqueta: "Tamaño del producto en pb (MÍNIMO-MÁXIMO)",
                    tipo: Tipo::Texto,
                    ejemplo: "100-300",
                },
                Campo {
                    id: "objetivo",
                    etiqueta: "Región que el producto debe contener, opcional (INICIO,LONGITUD; inicio desde 1)",
                    tipo: Tipo::Texto,
                    ejemplo: "1601,50",
                },
                Campo {
                    id: "tm_optima",
                    etiqueta: "Tm óptima (°C)",
                    tipo: Tipo::Decimal,
                    ejemplo: "60",
                },
            ],
            funcion: "cebadores::disenar_cebadores",
            ejecutar: disenar,
            pasos: Some(pasos_diseno),
        },
        Algoritmo {
            id: "puerto.temperatura_fusion",
            capitulo: "puerto",
            nombre: "Temperatura de fusión de cebadores",
            resumen: "Calcula la Tm de un cebador con vecinos cercanos y corrección por sales, y revisa si forma dímeros u horquillas.",
            campos: &[
                Campo {
                    id: "cebador",
                    etiqueta: "Cebador (5'→3')",
                    tipo: Tipo::Secuencia,
                    ejemplo: "TGTGCCGAAATGGTCCATCA",
                },
                Campo {
                    id: "pareja",
                    etiqueta: "Su pareja, opcional (5'→3')",
                    tipo: Tipo::Secuencia,
                    ejemplo: "GGCGTATCGCCAAAATCACC",
                },
                Campo {
                    id: "monovalentes",
                    etiqueta: "Sodio + potasio (mM)",
                    tipo: Tipo::Decimal,
                    ejemplo: "50",
                },
                Campo {
                    id: "magnesio",
                    etiqueta: "Magnesio, Mg²⁺ (mM)",
                    tipo: Tipo::Decimal,
                    ejemplo: "1.5",
                },
                Campo {
                    id: "dntp",
                    etiqueta: "dNTP (mM)",
                    tipo: Tipo::Decimal,
                    ejemplo: "0.6",
                },
                Campo {
                    id: "oligo",
                    etiqueta: "Concentración de cada cebador (nM)",
                    tipo: Tipo::Decimal,
                    ejemplo: "50",
                },
                Campo {
                    id: "correccion",
                    etiqueta: "Corrección por sales",
                    tipo: Tipo::Opcion(CORRECCIONES),
                    ejemplo: "SantaLucia 1998 (Primer3)",
                },
            ],
            funcion: "cebadores::temperatura_fusion",
            ejecutar: temperatura,
            pasos: Some(pasos_temperatura),
        },
    ]
}

fn grados(valor: f64) -> String {
    format!("{}\u{a0}°C", decimal(valor, 1))
}

fn energia(valor: f64) -> String {
    format!("{}\u{a0}kcal/mol", decimal(valor, 2).replacen('-', "−", 1))
}

fn energia_o_ninguno(union: &Apareamiento) -> String {
    if union.pares == 0 {
        "ninguno estable".to_string()
    } else {
        energia(union.energia)
    }
}

fn celda_de_energia(valor: f64) -> String {
    if valor == 0.0 {
        "—".to_string()
    } else {
        decimal(valor, 2)
    }
}

fn con_extremos(secuencia: &[u8]) -> String {
    format!("5'-{}-3'", String::from_utf8_lossy(secuencia))
}

fn texto(secuencia: &[u8]) -> String {
    String::from_utf8_lossy(secuencia).into_owned()
}

fn gc_en_la_pinza(secuencia: &[u8]) -> usize {
    let desde = secuencia.len().saturating_sub(BASES_DE_LA_PINZA);
    secuencia[desde..]
        .iter()
        .filter(|&&base| base == b'G' || base == b'C')
        .count()
}

fn leer_concentracion(entrada: &Entrada, campo: &str, maximo: f64) -> Result<f64, String> {
    let valor = entrada.decimal(campo)?;
    if !(0.0..=maximo).contains(&valor) {
        return Err(format!(
            "El campo «{campo}» debe estar entre 0 y {}.",
            miles(maximo as usize)
        ));
    }
    Ok(valor)
}

fn leer_condiciones(entrada: &Entrada) -> Result<Condiciones, String> {
    let correccion = if entrada.opcion("correccion", CORRECCIONES)? == CORRECCIONES[1] {
        CorreccionSal::Owczarzy
    } else {
        CorreccionSal::SantaLucia
    };
    let oligo = leer_concentracion(entrada, "oligo", 100_000.0)?;
    if oligo == 0.0 {
        return Err("La concentración del cebador debe ser mayor que cero.".to_string());
    }
    Ok(Condiciones {
        monovalentes_mm: leer_concentracion(entrada, "monovalentes", 1000.0)?,
        magnesio_mm: leer_concentracion(entrada, "magnesio", 100.0)?,
        dntp_mm: leer_concentracion(entrada, "dntp", 100.0)?,
        oligo_nm: oligo,
        correccion,
    })
}

fn filas_de_dimero(
    a: &[u8],
    nombre_a: &str,
    b: &[u8],
    nombre_b: &str,
    union: &Apareamiento,
) -> Vec<FilaAlineada> {
    let b_de_3_a_5: Vec<u8> = b.iter().rev().copied().collect();
    let enfrente = b.len() - 1 - union.en_b;
    let (sangria_a, sangria_b) = if enfrente >= union.en_a {
        (enfrente - union.en_a, 0)
    } else {
        (0, union.en_a - enfrente)
    };
    let ancho = (sangria_a + a.len()).max(sangria_b + b.len());
    let mut medio = vec![b' '; ancho];
    for paso in 0..union.pares {
        medio[sangria_a + union.en_a + paso] = b'|';
    }
    vec![
        FilaAlineada {
            nombre: nombre_a.to_string(),
            secuencia: rellenar(sangria_a, a, ancho),
        },
        FilaAlineada {
            nombre: String::new(),
            secuencia: texto(&medio),
        },
        FilaAlineada {
            nombre: nombre_b.to_string(),
            secuencia: rellenar(sangria_b, &b_de_3_a_5, ancho),
        },
    ]
}

fn rellenar(sangria: usize, letras: &[u8], ancho: usize) -> String {
    let mut fila = "-".repeat(sangria);
    fila.push_str(&texto(letras));
    while fila.len() < ancho {
        fila.push('-');
    }
    fila
}

fn filas_de_horquilla(oligo: &[u8], pliegue: &Horquilla) -> Vec<FilaAlineada> {
    let brazo_5 = &oligo[pliegue.inicio..pliegue.inicio + pliegue.pares];
    let brazo_3: Vec<u8> = oligo[pliegue.fin + 1 - pliegue.pares..=pliegue.fin]
        .iter()
        .rev()
        .copied()
        .collect();
    vec![
        FilaAlineada {
            nombre: "Brazo 5'→3'".to_string(),
            secuencia: texto(brazo_5),
        },
        FilaAlineada {
            nombre: String::new(),
            secuencia: "|".repeat(pliegue.pares),
        },
        FilaAlineada {
            nombre: "Brazo 3'→5'".to_string(),
            secuencia: texto(&brazo_3),
        },
    ]
}

fn preocupa(dimero: &cebadores::Dimero) -> bool {
    dimero.extremo_3.energia < UMBRAL_DIMERO_EXTREMO_3 || dimero.mas_estable.energia < UMBRAL_DIMERO
}

fn gravedad(dimero: &cebadores::Dimero) -> f64 {
    (UMBRAL_DIMERO_EXTREMO_3 - dimero.extremo_3.energia)
        .max(UMBRAL_DIMERO - dimero.mas_estable.energia)
}

fn union_a_mostrar(dimero: &cebadores::Dimero) -> Apareamiento {
    if dimero.extremo_3.energia < UMBRAL_DIMERO_EXTREMO_3 {
        dimero.extremo_3
    } else {
        dimero.mas_estable
    }
}

struct Interaccion<'a> {
    titulo: &'static str,
    a: &'a [u8],
    nombre_a: &'static str,
    b: &'a [u8],
    nombre_b: &'static str,
    dimero: cebadores::Dimero,
}

fn mostrar_el_peor_dimero(salida: &mut Salida, interacciones: &[Interaccion]) {
    let peor = interacciones
        .iter()
        .filter(|interaccion| preocupa(&interaccion.dimero))
        .max_by(|una, otra| gravedad(&una.dimero).total_cmp(&gravedad(&otra.dimero)));
    if let Some(peor) = peor {
        let union = union_a_mostrar(&peor.dimero);
        salida.agregar(Bloque::alineamiento(
            &format!("{} (ΔG = {})", peor.titulo, energia(union.energia)),
            filas_de_dimero(peor.a, peor.nombre_a, peor.b, peor.nombre_b, &union),
        ));
    }
}

fn explicar_tm_individual(salida: &mut Salida, tm: f64, sugerir_alineamiento: bool) {
    if (TM_BUENA_MINIMA..=TM_BUENA_MAXIMA).contains(&tm) {
        let sugerencia = if sugerir_alineamiento {
            format!(
                " Para empezar, usa una temperatura de alineamiento unos 5 °C por debajo, ≈ {}.",
                grados(tm - 5.0)
            )
        } else {
            String::new()
        };
        salida.explicar(format!(
            "Tm del cebador = {}: dentro del rango habitual para PCR (55–65 °C).{sugerencia}",
            grados(tm)
        ));
    } else if tm < TM_BUENA_MINIMA {
        salida.explicar(format!(
            "Tm del cebador = {}: baja para una PCR típica (55–65 °C). Alárgalo o muévelo a una zona con más G y C; si no, tendrás que alinear a una temperatura tan baja que también se pegará a sitios parecidos.",
            grados(tm)
        ));
    } else {
        salida.explicar(format!(
            "Tm del cebador = {}: alta (más de 65 °C). Suele indicar un cebador largo o muy rico en GC, que tolera desajustes y puede unirse a sitios parecidos; acórtalo o pásalo a una zona más equilibrada.",
            grados(tm)
        ));
    }
}

fn explicar_gc(salida: &mut Salida, gc: f64) {
    if gc < GC_MINIMO {
        salida.explicar(format!(
            "Su contenido GC es {}, menos del 40 % recomendado: las uniones A–T son más débiles, así que necesita ser más largo para alcanzar una buena Tm.",
            porcentaje(gc, 0)
        ));
    } else if gc > GC_MAXIMO {
        salida.explicar(format!(
            "Su contenido GC es {}, más del 60 % recomendado: se une con tanta fuerza que puede pegarse también donde no coincide del todo.",
            porcentaje(gc, 0)
        ));
    }
}

struct Oligo<'a> {
    nombre: &'static str,
    de_quien: &'static str,
    etiqueta: &'static str,
    secuencia: &'a [u8],
    tm: f64,
    duplex: cebadores::Duplex,
    propio: cebadores::Dimero,
    pliegue: Horquilla,
}

fn analizar_oligo<'a>(
    nombres: (&'static str, &'static str, &'static str),
    secuencia: &'a [u8],
    condiciones: &Condiciones,
) -> Result<Oligo<'a>, String> {
    let (nombre, de_quien, etiqueta) = nombres;
    Ok(Oligo {
        nombre,
        de_quien,
        etiqueta,
        secuencia,
        tm: cebadores::temperatura_fusion(secuencia, condiciones)?,
        duplex: cebadores::duplex(secuencia),
        propio: cebadores::dimero(secuencia, secuencia, condiciones)?,
        pliegue: cebadores::horquilla(secuencia, condiciones)?,
    })
}

fn fila_de_oligo(oligo: &Oligo) -> Vec<String> {
    vec![
        format!("{}: {}", oligo.etiqueta, texto(oligo.secuencia)),
        miles(oligo.secuencia.len()),
        porcentaje(
            oligo.duplex.bases_gc as f64 / oligo.secuencia.len() as f64,
            0,
        ),
        decimal(oligo.tm, 1),
        decimal(oligo.duplex.entalpia, 1),
        decimal(oligo.duplex.entropia, 1),
        celda_de_energia(oligo.propio.mas_estable.energia),
        celda_de_energia(oligo.propio.extremo_3.energia),
        celda_de_energia(oligo.pliegue.energia),
    ]
}

fn problemas_propios(oligo: &Oligo) -> Vec<String> {
    let mut problemas = Vec::new();
    if oligo.propio.extremo_3.energia < UMBRAL_DIMERO_EXTREMO_3 {
        problemas.push(format!(
            "{} forma un autodímero que toca su extremo 3' (ΔG = {})",
            oligo.nombre,
            energia(oligo.propio.extremo_3.energia)
        ));
    }
    if oligo.propio.mas_estable.energia < UMBRAL_DIMERO {
        problemas.push(format!(
            "{} forma un autodímero fuerte (ΔG = {})",
            oligo.nombre,
            energia(oligo.propio.mas_estable.energia)
        ));
    }
    if oligo.pliegue.energia < UMBRAL_HORQUILLA {
        problemas.push(format!(
            "{} se pliega en una horquilla de {} pb con un lazo de {} bases (ΔG = {})",
            oligo.nombre,
            oligo.pliegue.pares,
            oligo.pliegue.lazo,
            energia(oligo.pliegue.energia)
        ));
    }
    problemas
}

fn problemas_cruzados(cruzado: &cebadores::Dimero) -> Vec<String> {
    let mut problemas = Vec::new();
    if cruzado.extremo_3.energia < UMBRAL_DIMERO_EXTREMO_3 {
        problemas.push(format!(
            "el cebador y su pareja se aparean tocando un extremo 3' (ΔG = {})",
            energia(cruzado.extremo_3.energia)
        ));
    }
    if cruzado.mas_estable.energia < UMBRAL_DIMERO {
        problemas.push(format!(
            "el cebador y su pareja forman un dímero fuerte (ΔG = {})",
            energia(cruzado.mas_estable.energia)
        ));
    }
    problemas
}

fn explicar_estructuras(
    salida: &mut Salida,
    oligos: &[&Oligo],
    cruzado: Option<&cebadores::Dimero>,
) {
    let mut problemas = Vec::new();
    let mut dimeros: Vec<cebadores::Dimero> = Vec::new();
    for oligo in oligos {
        problemas.extend(problemas_propios(oligo));
        dimeros.push(oligo.propio);
    }
    if let Some(cruzado) = cruzado {
        problemas.extend(problemas_cruzados(cruzado));
        dimeros.push(*cruzado);
    }
    if problemas.is_empty() {
        let mas_estable = dimeros
            .iter()
            .map(|dimero| dimero.mas_estable.energia)
            .fold(0.0, f64::min);
        let en_3 = dimeros
            .iter()
            .map(|dimero| dimero.extremo_3.energia)
            .fold(0.0, f64::min);
        let extremo = if en_3 < 0.0 {
            format!(
                "el peor que toca un extremo 3', {} (hasta −5)",
                energia(en_3)
            )
        } else {
            "ninguno estable toca un extremo 3'".to_string()
        };
        salida.explicar(format!(
            "Sin dímeros ni horquillas preocupantes: el apareamiento más estable tiene ΔG = {} (se tolera hasta −9 kcal/mol) y {extremo}.",
            energia(mas_estable)
        ));
    } else {
        salida.explicar(format!(
            "Ojo: {}. Un oligo pegado a otro o plegado no se une al molde y, si su extremo 3' queda apareado, la polimerasa lo alarga y aparecen dímeros de cebador (bandas cortas en el gel).",
            problemas.join("; ")
        ));
    }
}

fn explicar_diferencia_de_tm(salida: &mut Salida, tm: f64, tm_otro: f64) {
    let diferencia = (tm - tm_otro).abs();
    let alineamiento = grados(tm.min(tm_otro) - 5.0);
    if diferencia <= DIFERENCIA_IDEAL {
        salida.explicar(format!(
            "Las Tm del par difieren en {}: ideal (2 °C o menos), los dos se unen bien a la misma temperatura de alineamiento, ≈ {alineamiento}.",
            grados(diferencia)
        ));
    } else if diferencia <= DIFERENCIA_PROBLEMATICA {
        salida.explicar(format!(
            "Las Tm del par difieren en {}: aceptable, pero alinea pensando en el más débil (≈ {alineamiento}); a esa temperatura el otro puede unirse a sitios parecidos.",
            grados(diferencia)
        ));
    } else {
        salida.explicar(format!(
            "Las Tm del par difieren en {}: más de 5 °C es problemático. A una temperatura buena para uno, el otro o no se une o se une en sitios equivocados; rediseña uno de los dos.",
            grados(diferencia)
        ));
    }
}

fn exigir_longitud_minima(oligo: &[u8], nombre: &str) -> Result<(), String> {
    if oligo.len() < LONGITUD_MINIMA_DEL_OLIGO {
        return Err(format!(
            "{nombre} mide {} bases y necesita al menos {LONGITUD_MINIMA_DEL_OLIGO}: con menos, el modelo de vecinos cercanos no da una Tm fiable.",
            oligo.len()
        ));
    }
    Ok(())
}

fn temperatura(entrada: &Entrada) -> Result<Salida, String> {
    let cebador = entrada.adn("cebador")?;
    exigir_longitud_minima(&cebador, "El cebador")?;
    let pareja = if entrada.tiene("pareja") {
        Some(entrada.adn("pareja")?)
    } else {
        None
    };
    if let Some(otro) = &pareja {
        exigir_longitud_minima(otro, "La pareja")?;
    }
    let condiciones = leer_condiciones(entrada)?;
    let principal = analizar_oligo(
        ("el cebador", "del cebador", "Cebador"),
        &cebador,
        &condiciones,
    )?;
    let segundo = match &pareja {
        Some(otro) => Some(analizar_oligo(
            ("la pareja", "de la pareja", "Pareja"),
            otro,
            &condiciones,
        )?),
        None => None,
    };
    let cruzado = match &pareja {
        Some(otro) => Some(cebadores::dimero(&cebador, otro, &condiciones)?),
        None => None,
    };
    let curva = cebadores::curva_de_fusion(&cebador, &condiciones)?;

    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Tm del cebador", grados(principal.tm)));
    let mut oligos = vec![&principal];
    let mut interacciones = vec![Interaccion {
        titulo: "Autodímero del cebador",
        a: &cebador,
        nombre_a: "Cebador 5'→3'",
        b: &cebador,
        nombre_b: "Cebador 3'→5'",
        dimero: principal.propio,
    }];
    if let (Some(otro), Some(cruzado)) = (&segundo, &cruzado) {
        salida.agregar(Bloque::valor("Tm de la pareja", grados(otro.tm)));
        salida.agregar(Bloque::valor(
            "Diferencia de Tm",
            grados((principal.tm - otro.tm).abs()),
        ));
        salida.agregar(Bloque::valor(
            "Dímero más estable entre ambos",
            energia_o_ninguno(&cruzado.mas_estable),
        ));
        salida.agregar(Bloque::valor(
            "Dímero entre ambos que toca un extremo 3'",
            energia_o_ninguno(&cruzado.extremo_3),
        ));
        oligos.push(otro);
        interacciones.push(Interaccion {
            titulo: "Autodímero de la pareja",
            a: otro.secuencia,
            nombre_a: "Pareja 5'→3'",
            b: otro.secuencia,
            nombre_b: "Pareja 3'→5'",
            dimero: otro.propio,
        });
        interacciones.push(Interaccion {
            titulo: "Dímero entre cebador y pareja",
            a: &cebador,
            nombre_a: "Cebador 5'→3'",
            b: otro.secuencia,
            nombre_b: "Pareja 3'→5'",
            dimero: *cruzado,
        });
    }
    salida.agregar(Bloque::tabla(
        "Propiedades de cada oligo",
        &[
            "Oligo",
            "Longitud (nt)",
            "GC",
            "Tm (°C)",
            "ΔH (kcal/mol)",
            "ΔS (cal/K·mol)",
            "Autodímero ΔG",
            "Autodímero en 3' ΔG",
            "Horquilla ΔG",
        ],
        oligos.iter().map(|oligo| fila_de_oligo(oligo)).collect(),
    ));
    let puntos: Vec<[f64; 2]> = curva
        .iter()
        .map(|(temperatura, fraccion)| [*temperatura, 100.0 * fraccion])
        .collect();
    salida.agregar(Bloque::linea(
        "Curva de fusión del cebador",
        "Temperatura (°C)",
        "Cebador en doble hélice (%)",
        &puntos,
        vec![Marca {
            x: principal.tm,
            texto: "Tm".to_string(),
        }],
    ));
    mostrar_el_peor_dimero(&mut salida, &interacciones);
    for oligo in &oligos {
        if oligo.pliegue.energia < UMBRAL_HORQUILLA {
            salida.agregar(Bloque::alineamiento(
                &format!(
                    "Horquilla {} (lazo de {} bases, ΔG = {})",
                    oligo.de_quien,
                    oligo.pliegue.lazo,
                    energia(oligo.pliegue.energia)
                ),
                filas_de_horquilla(oligo.secuencia, &oligo.pliegue),
            ));
        }
    }

    explicar_tm_individual(&mut salida, principal.tm, segundo.is_none());
    explicar_gc(
        &mut salida,
        principal.duplex.bases_gc as f64 / cebador.len() as f64,
    );
    explicar_estructuras(&mut salida, &oligos, cruzado.as_ref());
    match &segundo {
        Some(otro) => explicar_diferencia_de_tm(&mut salida, principal.tm, otro.tm),
        None => explicar_magnesio(&mut salida, &cebador, principal.tm, &condiciones),
    }
    if let Some((_, fraccion)) = curva.get(40) {
        salida.explicar(format!(
            "En la curva, a la Tm la mitad de las moléculas de cebador está en doble hélice; 5 °C por debajo ya lo está el {}. Por eso se alinea un poco por debajo de la Tm.",
            porcentaje(*fraccion, 0)
        ));
    }
    Ok(salida)
}

fn explicar_magnesio(salida: &mut Salida, cebador: &[u8], tm: f64, condiciones: &Condiciones) {
    if condiciones.magnesio_mm <= condiciones.dntp_mm {
        if condiciones.dntp_mm > 0.0 {
            salida.explicar(format!(
                "Con {} mM de Mg²⁺ y {} mM de dNTP no queda magnesio libre: los dNTP lo atrapan todo. La polimerasa necesita Mg²⁺ libre para copiar, así que en una PCR pon más Mg²⁺ que dNTP.",
                decimal(condiciones.magnesio_mm, 1),
                decimal(condiciones.dntp_mm, 1)
            ));
        }
        return;
    }
    let sin_magnesio = Condiciones {
        magnesio_mm: 0.0,
        dntp_mm: 0.0,
        ..*condiciones
    };
    if let Ok(tm_sin_magnesio) = cebadores::temperatura_fusion(cebador, &sin_magnesio) {
        if tm - tm_sin_magnesio < MEJORA_MINIMA_DEL_MAGNESIO {
            salida.explicar(format!(
                "Con {} mM de Mg²⁺ y {} mM de dNTP queda tan poco magnesio libre que casi no cambia la Tm: sin él sería {}.",
                decimal(condiciones.magnesio_mm, 1),
                decimal(condiciones.dntp_mm, 1),
                grados(tm_sin_magnesio)
            ));
            return;
        }
        salida.explicar(format!(
            "Con {} mM de Mg²⁺ y {} mM de dNTP queda magnesio libre, que estabiliza la doble hélice: sin él la Tm sería {}, es decir {} menos. Los dNTP atrapan parte del Mg²⁺, por eso importa la diferencia entre ambos.",
            decimal(condiciones.magnesio_mm, 1),
            decimal(condiciones.dntp_mm, 1),
            grados(tm_sin_magnesio),
            grados(tm - tm_sin_magnesio)
        ));
    }
}

struct Origen {
    primera: usize,
    inversa: bool,
}

impl Origen {
    fn en_ecoli(&self, posicion: usize) -> usize {
        if self.inversa {
            self.primera + 1 - posicion
        } else {
            self.primera + posicion - 1
        }
    }
}

fn buscar_en(genoma: &[u8], plantilla: &[u8]) -> Option<usize> {
    if plantilla.is_empty() || plantilla.len() > genoma.len() {
        return None;
    }
    genoma
        .windows(plantilla.len())
        .position(|ventana| ventana == plantilla)
}

fn origen_en_ecoli(escrito: &str, plantilla: &[u8]) -> Option<Origen> {
    let escrito = escrito.trim();
    if let Some(nombre) = escrito.strip_prefix("@gen:") {
        let elegido = ecoli::buscar_gen(nombre.trim())?;
        return Some(if elegido.hebra == b'-' {
            Origen {
                primera: elegido.fin,
                inversa: true,
            }
        } else {
            Origen {
                primera: elegido.inicio,
                inversa: false,
            }
        });
    }
    if let Some(rango) = escrito.strip_prefix("@ecoli:") {
        let (inicio, _) = rango.split_once('-')?;
        return Some(Origen {
            primera: inicio.trim().parse().ok()?,
            inversa: false,
        });
    }
    let genoma = ecoli::genoma();
    if let Some(lugar) = buscar_en(genoma, plantilla) {
        return Some(Origen {
            primera: lugar + 1,
            inversa: false,
        });
    }
    buscar_en(genoma, &complemento_inverso(plantilla)).map(|lugar| Origen {
        primera: lugar + plantilla.len(),
        inversa: true,
    })
}

fn leer_par_de_numeros(texto: &str, campo: &str, ejemplo: &str) -> Result<(usize, usize), String> {
    let numeros: Vec<&str> = texto
        .split(|c: char| c == '-' || c == ',' || c == ';' || c == '–' || c.is_whitespace())
        .filter(|pieza| !pieza.is_empty())
        .collect();
    let mensaje = format!("Escribe el campo «{campo}» como {ejemplo}.");
    if numeros.len() != 2 {
        return Err(mensaje);
    }
    let primero = numeros[0].parse::<usize>().map_err(|_| mensaje.clone())?;
    let segundo = numeros[1].parse::<usize>().map_err(|_| mensaje.clone())?;
    Ok((primero, segundo))
}

fn leer_opciones(entrada: &Entrada, plantilla: &[u8]) -> Result<OpcionesDiseno, String> {
    let (producto_minimo, producto_maximo) = leer_par_de_numeros(
        &entrada.texto("producto")?,
        "producto",
        "MÍNIMO-MÁXIMO, por ejemplo 100-300",
    )?;
    if producto_minimo < PRODUCTO_MINIMO_PERMITIDO
        || producto_minimo > producto_maximo
        || producto_maximo > PRODUCTO_MAXIMO
    {
        return Err(format!(
            "El producto debe cumplir {PRODUCTO_MINIMO_PERMITIDO} ≤ mínimo ≤ máximo ≤ {} pb.",
            miles(PRODUCTO_MAXIMO)
        ));
    }
    let tm_optima = entrada.decimal("tm_optima")?;
    if !(45.0..=75.0).contains(&tm_optima) {
        return Err("La Tm óptima debe estar entre 45 y 75 °C.".to_string());
    }
    let objetivo = if entrada.tiene("objetivo") {
        let (inicio, longitud) = leer_par_de_numeros(
            &entrada.texto("objetivo")?,
            "objetivo",
            "INICIO,LONGITUD, por ejemplo 1601,50",
        )?;
        if inicio == 0 || longitud == 0 || (inicio - 1).saturating_add(longitud) > plantilla.len() {
            return Err(format!(
                "La región objetivo debe empezar en 1 o más, medir al menos 1 base y terminar dentro de la plantilla ({} pb).",
                miles(plantilla.len())
            ));
        }
        Some((inicio - 1, longitud))
    } else {
        None
    };
    Ok(OpcionesDiseno {
        tm_minima: tm_optima - MARGEN_DE_TM,
        tm_optima,
        tm_maxima: tm_optima + MARGEN_DE_TM,
        producto_minimo,
        producto_maximo,
        objetivo,
        ..OpcionesDiseno::recomendadas()
    })
}

fn leer_plantilla(entrada: &Entrada) -> Result<Vec<u8>, String> {
    let plantilla = entrada.secuencia("plantilla")?;
    if let Some(letra) = plantilla
        .iter()
        .find(|letra| !matches!(letra, b'A' | b'C' | b'G' | b'T' | b'N'))
    {
        return Err(format!(
            "La plantilla solo puede tener A, C, G, T y N, pero contiene «{}».",
            *letra as char
        ));
    }
    if plantilla.len() > PLANTILLA_MAXIMA {
        return Err(format!(
            "La plantilla mide {} pb y el límite es {} pb: recorta la región que quieres amplificar dejando unos cientos de bases a cada lado.",
            miles(plantilla.len()),
            miles(PLANTILLA_MAXIMA)
        ));
    }
    Ok(plantilla)
}

fn tramo(cebador: &Cebador) -> String {
    format!(
        "{}–{}",
        miles(cebador.inicio + 1),
        miles(cebador.inicio + cebador.longitud)
    )
}

fn sitios(par: &ParDeCebadores) -> String {
    match (par.izquierdo.sitios, par.derecho.sitios) {
        (Some(izquierdo), Some(derecho)) => format!("{izquierdo} y {derecho}"),
        _ => "—".to_string(),
    }
}

fn peor_dimero(par: &ParDeCebadores) -> f64 {
    par.dimero
        .min(par.izquierdo.autodimero)
        .min(par.derecho.autodimero)
}

fn peor_dimero_en_3(par: &ParDeCebadores) -> f64 {
    par.dimero_extremo_3
        .min(par.izquierdo.autodimero_extremo_3)
        .min(par.derecho.autodimero_extremo_3)
}

fn fila_de_par(numero: usize, par: &ParDeCebadores) -> Vec<String> {
    vec![
        miles(numero),
        texto(&par.izquierdo.secuencia),
        tramo(&par.izquierdo),
        decimal(par.izquierdo.tm, 1),
        porcentaje(par.izquierdo.gc, 0),
        texto(&par.derecho.secuencia),
        tramo(&par.derecho),
        decimal(par.derecho.tm, 1),
        porcentaje(par.derecho.gc, 0),
        miles(par.producto),
        celda_de_energia(peor_dimero(par)),
        celda_de_energia(peor_dimero_en_3(par)),
        celda_de_energia(par.izquierdo.horquilla.min(par.derecho.horquilla)),
        sitios(par),
        decimal(par.penalizacion, 2),
    ]
}

fn paisaje_de_tm(plantilla: &[u8], condiciones: &Condiciones) -> Vec<[f64; 2]> {
    let mut puntos = Vec::new();
    for inicio in 0..plantilla.len().saturating_sub(VENTANA_DEL_PAISAJE - 1) {
        let ventana = &plantilla[inicio..inicio + VENTANA_DEL_PAISAJE];
        if let Ok(tm) = cebadores::temperatura_fusion(ventana, condiciones) {
            puntos.push([(inicio + 1) as f64, tm]);
        }
    }
    puntos
}

fn marcas_del_par(par: Option<&ParDeCebadores>, opciones: &OpcionesDiseno) -> Vec<Marca> {
    let mut marcas = Vec::new();
    if let Some(par) = par {
        marcas.push(Marca {
            x: (par.izquierdo.inicio + 1) as f64,
            texto: "cebador izquierdo".to_string(),
        });
        marcas.push(Marca {
            x: (par.derecho.inicio + par.derecho.longitud) as f64,
            texto: "cebador derecho".to_string(),
        });
    }
    if let Some((inicio, longitud)) = opciones.objetivo {
        marcas.push(Marca {
            x: (inicio + 1) as f64,
            texto: "inicio del objetivo".to_string(),
        });
        marcas.push(Marca {
            x: (inicio + longitud) as f64,
            texto: "fin del objetivo".to_string(),
        });
    }
    marcas
}

fn disenar(entrada: &Entrada) -> Result<Salida, String> {
    let plantilla = leer_plantilla(entrada)?;
    let opciones = leer_opciones(entrada, &plantilla)?;
    let origen = origen_en_ecoli(&entrada.texto("plantilla")?, &plantilla);
    let genoma = origen.as_ref().map(|_| ecoli::genoma());
    let diseno = cebadores::disenar_cebadores(&plantilla, &opciones, genoma)?;

    let mut salida = Salida::default();
    let mejor = diseno.pares.first();
    match mejor {
        Some(par) => {
            salida.agregar(Bloque::valor(
                "Cebador izquierdo",
                con_extremos(&par.izquierdo.secuencia),
            ));
            salida.agregar(Bloque::valor(
                "Cebador derecho",
                con_extremos(&par.derecho.secuencia),
            ));
            salida.agregar(Bloque::valor(
                "Producto",
                format!("{} pb", miles(par.producto)),
            ));
            salida.agregar(Bloque::valor(
                "Tm del par",
                format!("{} y {}", grados(par.izquierdo.tm), grados(par.derecho.tm)),
            ));
        }
        None => salida.agregar(Bloque::valor("Pares encontrados", "0")),
    }
    if !diseno.pares.is_empty() {
        let filas: Vec<Vec<String>> = diseno
            .pares
            .iter()
            .enumerate()
            .map(|(lugar, par)| fila_de_par(lugar + 1, par))
            .collect();
        salida.agregar(Bloque::tabla(
            &if diseno.pares.len() == 1 {
                "El mejor par (posiciones desde 1 en la plantilla)".to_string()
            } else {
                format!(
                    "Los {} mejores pares (posiciones desde 1 en la plantilla)",
                    diseno.pares.len()
                )
            },
            &[
                "#",
                "Izquierdo (5'→3')",
                "Posición izq.",
                "Tm izq. (°C)",
                "GC izq.",
                "Derecho (5'→3')",
                "Posición der.",
                "Tm der. (°C)",
                "GC der.",
                "Producto (pb)",
                "Peor dímero ΔG",
                "Peor dímero en 3' ΔG",
                "Peor horquilla ΔG",
                "Sitios en E. coli",
                "Penalización",
            ],
            filas,
        ));
    }
    salida.agregar(Bloque::linea(
        "Tm de cada tramo de 20 bases de la plantilla",
        "Posición en la plantilla (desde 1)",
        "Tm (°C)",
        &paisaje_de_tm(&plantilla, &opciones.condiciones),
        marcas_del_par(mejor, &opciones),
    ));
    agregar_descartes(&mut salida, &diseno, &opciones);
    if let Some(par) = mejor {
        let condiciones = &opciones.condiciones;
        let izquierdo = &par.izquierdo.secuencia;
        let derecho = &par.derecho.secuencia;
        let interacciones = [
            Interaccion {
                titulo: "Dímero entre los cebadores del mejor par",
                a: izquierdo,
                nombre_a: "Izquierdo 5'→3'",
                b: derecho,
                nombre_b: "Derecho 3'→5'",
                dimero: cebadores::dimero(izquierdo, derecho, condiciones)?,
            },
            Interaccion {
                titulo: "Autodímero del cebador izquierdo",
                a: izquierdo,
                nombre_a: "Izquierdo 5'→3'",
                b: izquierdo,
                nombre_b: "Izquierdo 3'→5'",
                dimero: cebadores::dimero(izquierdo, izquierdo, condiciones)?,
            },
            Interaccion {
                titulo: "Autodímero del cebador derecho",
                a: derecho,
                nombre_a: "Derecho 5'→3'",
                b: derecho,
                nombre_b: "Derecho 3'→5'",
                dimero: cebadores::dimero(derecho, derecho, condiciones)?,
            },
        ];
        mostrar_el_peor_dimero(&mut salida, &interacciones);
        explicar_mejor_par(&mut salida, par, &opciones, origen.as_ref());
    } else {
        explicar_sin_pares(&mut salida, &diseno, &opciones);
    }
    Ok(salida)
}

fn agregar_descartes(salida: &mut Salida, diseno: &Diseno, opciones: &OpcionesDiseno) {
    let descartes = &diseno.descartes;
    salida.agregar(Bloque::barras(
        "Qué pasó con cada candidato (las dos hebras, primera regla que falla)",
        "Regla",
        "Candidatos",
        vec![
            "Tiene N".to_string(),
            "Racha de más de 4 bases iguales".to_string(),
            "GC fuera de 40–60 %".to_string(),
            "Pinza GC: no 1–3 G/C en las últimas 5".to_string(),
            format!(
                "Tm fuera de {}–{} °C",
                decimal(opciones.tm_minima, 0),
                decimal(opciones.tm_maxima, 0)
            ),
            "Aceptados".to_string(),
        ],
        vec![
            descartes.bases_ambiguas as f64,
            descartes.repeticiones as f64,
            descartes.contenido_gc as f64,
            descartes.pinza_gc as f64,
            descartes.temperatura as f64,
            (diseno.izquierdos + diseno.derechos) as f64,
        ],
    ));
}

fn explicar_sin_pares(salida: &mut Salida, diseno: &Diseno, opciones: &OpcionesDiseno) {
    let descartes = &diseno.descartes;
    let motivos = [
        (descartes.bases_ambiguas, "tienen bases N"),
        (
            descartes.repeticiones,
            "tienen rachas de más de 4 bases iguales",
        ),
        (descartes.contenido_gc, "tienen un GC fuera de 40–60 %"),
        (
            descartes.pinza_gc,
            "no tienen una buena pinza GC en su extremo 3'",
        ),
        (descartes.temperatura, "tienen la Tm fuera del rango pedido"),
    ];
    if diseno.izquierdos == 0 || diseno.derechos == 0 {
        let principal = motivos
            .iter()
            .max_by_key(|(cuantos, _)| *cuantos)
            .map(|(cuantos, motivo)| format!("{} candidatos {motivo}", miles(*cuantos)))
            .unwrap_or_default();
        let cual = match (diseno.izquierdos, diseno.derechos) {
            (0, 0) => "ni izquierdo ni derecho",
            (0, _) => "izquierdo",
            _ => "derecho",
        };
        salida.explicar(format!(
            "No hay ningún par: ningún tramo de la plantilla sirve como cebador {cual}. El motivo más común: {principal}."
        ));
        if descartes.temperatura >= descartes.contenido_gc {
            salida.explicar(
                "Prueba con otra Tm óptima: una plantilla rica en A y T da cebadores de Tm baja, y una rica en G y C, de Tm alta.",
            );
        } else {
            salida.explicar(
                "La composición de la plantilla está lejos del 40–60 % de GC; amplía la región para que entren zonas más equilibradas.",
            );
        }
    } else {
        let objetivo = if opciones.objetivo.is_some() {
            ", que rodee la región objetivo"
        } else {
            ""
        };
        salida.explicar(format!(
            "Hay {} candidatos izquierdos y {} derechos, pero ninguna combinación da un producto de {}–{} pb{objetivo} con Tm que difieran 3 °C o menos.",
            miles(diseno.izquierdos),
            miles(diseno.derechos),
            miles(opciones.producto_minimo),
            miles(opciones.producto_maximo)
        ));
        salida.explicar(
            "Amplía el rango del producto o mueve la región objetivo: los cebadores necesitan sitio a ambos lados.",
        );
    }
}

fn explicar_mejor_par(
    salida: &mut Salida,
    par: &ParDeCebadores,
    opciones: &OpcionesDiseno,
    origen: Option<&Origen>,
) {
    let tm_izquierdo = par.izquierdo.tm;
    let tm_derecho = par.derecho.tm;
    let ambas_buenas = [tm_izquierdo, tm_derecho]
        .iter()
        .all(|tm| (TM_BUENA_MINIMA..=TM_BUENA_MAXIMA).contains(tm));
    let rango = if ambas_buenas {
        "dentro del rango habitual de 55–65 °C"
    } else {
        "fuera del rango habitual de 55–65 °C, así que ajusta bien la temperatura de alineamiento"
    };
    salida.explicar(format!(
        "El mejor par tiene Tm de {} y {}, {rango}. Para empezar, alinea a unos 5 °C por debajo de la más baja, ≈ {}.",
        grados(tm_izquierdo),
        grados(tm_derecho),
        grados(tm_izquierdo.min(tm_derecho) - 5.0)
    ));
    let diferencia = par.diferencia_tm;
    if diferencia <= DIFERENCIA_IDEAL {
        salida.explicar(format!(
            "Sus Tm difieren en {}: ideal (2 °C o menos), así ambos se unen igual de bien a la misma temperatura.",
            grados(diferencia)
        ));
    } else {
        salida.explicar(format!(
            "Sus Tm difieren en {}: aceptable, aunque lo ideal son 2 °C o menos; a la temperatura del más débil, el otro puede unirse a sitios parecidos.",
            grados(diferencia)
        ));
    }
    let pinzas = [
        gc_en_la_pinza(&par.izquierdo.secuencia),
        gc_en_la_pinza(&par.derecho.secuencia),
    ];
    if pinzas.contains(&PINZA_MAXIMA) {
        let cuales = match pinzas {
            [PINZA_MAXIMA, PINZA_MAXIMA] => "Los dos llevan",
            [PINZA_MAXIMA, _] => "El izquierdo lleva",
            _ => "El derecho lleva",
        };
        salida.explicar(format!(
            "{cuales} 3 G/C en sus últimas 5 bases: es el máximo aceptable. La pinza GC ayuda a que el extremo 3' quede bien unido, pero con más de 3 se pega con fuerza también a sitios equivocados."
        ));
    } else {
        salida.explicar(format!(
            "Terminan con {} y {} G/C en sus últimas 5 bases: una pinza GC moderada que ayuda a que el extremo 3', donde empieza a copiar la polimerasa, quede bien unido sin volverse pegajoso.",
            pinzas[0], pinzas[1]
        ));
    }
    let peor_en_3 = peor_dimero_en_3(par);
    let peor = peor_dimero(par);
    if peor_en_3 < UMBRAL_DIMERO_EXTREMO_3 || peor < UMBRAL_DIMERO {
        salida.explicar(format!(
            "Riesgo de dímeros de cebador: el dímero más estable tiene ΔG = {} y el que toca un extremo 3', {} (umbrales −9 y −5). Si en el gel aparece una banda por debajo de 100 pb, son dímeros; prueba otro par de la tabla.",
            energia(peor),
            energia(peor_en_3)
        ));
    } else {
        salida.explicar(format!(
            "Ningún dímero preocupante: el más estable tiene ΔG = {} (se tolera hasta −9 kcal/mol) y el que toca un extremo 3', {} (hasta −5).",
            energia(peor),
            energia(peor_en_3)
        ));
    }
    match (par.izquierdo.sitios, par.derecho.sitios, origen) {
        (Some(1), Some(1), _) => salida.explicar(format!(
            "Los últimos {BASES_DEL_EXTREMO_3} nt de cada cebador aparecen una sola vez en el genoma de E. coli, contando ambas hebras: el producto debería ser único."
        )),
        (Some(izquierdo), Some(derecho), _) => salida.explicar(format!(
            "Los últimos {BASES_DEL_EXTREMO_3} nt del cebador izquierdo aparecen {izquierdo} veces en E. coli y los del derecho {derecho}: donde hay más de un sitio puede salir un producto inespecífico."
        )),
        _ => salida.explicar(
            "La plantilla no viene de E. coli, así que no se revisó la especificidad: busca tus cebadores con BLASTN en el genoma del organismo que vas a amplificar.",
        ),
    }
    let uso_del_producto = if (QPCR_MINIMO..=QPCR_MAXIMO).contains(&par.producto) {
        "sirve también para PCR en tiempo real (qPCR), que prefiere productos de 70–200 pb"
    } else {
        "es para PCR convencional: lo verás como una banda en un gel de agarosa"
    };
    let ubicacion = match origen {
        Some(lugar) => {
            let uno = lugar.en_ecoli(par.izquierdo.inicio + 1);
            let otro = lugar.en_ecoli(par.derecho.inicio + par.derecho.longitud);
            format!(
                " En el genoma de E. coli (contando desde 1) ocupa de {} a {}.",
                miles(uno.min(otro)),
                miles(uno.max(otro))
            )
        }
        None => String::new(),
    };
    let objetivo = match opciones.objetivo {
        Some((inicio, longitud)) => format!(
            " Contiene la región objetivo ({}–{} de la plantilla).",
            miles(inicio + 1),
            miles(inicio + longitud)
        ),
        None => String::new(),
    };
    salida.explicar(format!(
        "El producto mide {} pb y {uso_del_producto}.{objetivo}{ubicacion}",
        miles(par.producto)
    ));
}

const PARES_DETALLADOS: usize = 4;
const BASES_MOSTRADAS: usize = 28;

fn numero(valor: f64, cifras: usize) -> String {
    decimal(valor, cifras).replacen('-', "−", 1)
}

fn corto(secuencia: &[u8]) -> String {
    if secuencia.len() <= BASES_MOSTRADAS {
        texto(secuencia)
    } else {
        format!("{}…", texto(&secuencia[..BASES_MOSTRADAS]))
    }
}

fn concentracion_en_molar(duplex: &cebadores::Duplex, condiciones: &Condiciones) -> f64 {
    let oligo = condiciones.oligo_nm * 1e-9;
    if duplex.autocomplementario {
        oligo
    } else {
        oligo / 4.0
    }
}

fn tm_a_un_molar(duplex: &cebadores::Duplex, condiciones: &Condiciones) -> f64 {
    let concentracion = concentracion_en_molar(duplex, condiciones);
    1000.0 * duplex.entalpia
        / (duplex.entropia + cebadores::CONSTANTE_DE_LOS_GASES * concentracion.ln())
}

fn pizarra_de_condiciones(pizarra: Pizarra, condiciones: &Condiciones) -> Pizarra {
    pizarra
        .variable(
            "condiciones.monovalentes_mm",
            decimal(condiciones.monovalentes_mm, 1),
        )
        .variable(
            "condiciones.magnesio_mm",
            decimal(condiciones.magnesio_mm, 1),
        )
        .variable("condiciones.dntp_mm", decimal(condiciones.dntp_mm, 1))
        .variable("condiciones.oligo_nm", decimal(condiciones.oligo_nm, 1))
}

fn pasos_temperatura(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let cebador = entrada.adn("cebador")?;
    exigir_longitud_minima(&cebador, "El cebador")?;
    let pareja = if entrada.tiene("pareja") {
        Some(entrada.adn("pareja")?)
    } else {
        None
    };
    if let Some(otro) = &pareja {
        exigir_longitud_minima(otro, "La pareja")?;
    }
    let condiciones = leer_condiciones(entrada)?;
    let tm = cebadores::temperatura_fusion(&cebador, &condiciones)?;
    let tm_de_la_pareja = match &pareja {
        Some(otro) => Some(cebadores::temperatura_fusion(otro, &condiciones)?),
        None => None,
    };
    let longitud = cebador.len();
    let todo = [(0, longitud, Resalte::Actual)];
    let owczarzy = condiciones.correccion == CorreccionSal::Owczarzy;
    let mut narracion = Narracion::nueva(400);

    narracion.agregar(
        pizarra_de_condiciones(
            Pizarra::nueva()
                .variable("oligo", corto(&cebador))
                .variable("oligo.len()", longitud),
            &condiciones,
        )
        .secuencia(&cebador, &todo)
        .paso(
            "validar_oligo(oligo)?;",
            format!(
                "Tienes un cebador de {longitud} bases, {}, y quieres su Tm: la temperatura a la que la mitad de las moléculas está pegada a su complemento y la otra mitad suelta. Primero el código comprueba que el oligo tenga entre 2 y 60 bases y solo A, C, G o T. Todo en regla.",
                corto(&cebador)
            ),
        ),
    );
    narracion.agregar(
        pizarra_de_condiciones(Pizarra::nueva(), &condiciones).paso(
            "validar_condiciones(condiciones)?;",
            format!(
                "Ahora revisa tus condiciones de reacción: {} mM de sodio + potasio, {} mM de Mg²⁺, {} mM de dNTP y {} nM de cada cebador. Ninguna es negativa y hay sal suficiente para que el ADN forme doble hélice.",
                decimal(condiciones.monovalentes_mm, 1),
                decimal(condiciones.magnesio_mm, 1),
                decimal(condiciones.dntp_mm, 1),
                decimal(condiciones.oligo_nm, 1)
            ),
        ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("entalpia", "0.0")
            .variable("entropia", "0.0")
            .secuencia(&cebador, &todo)
            .paso(
                "let mut entalpia = 0.0;",
                "La Tm sale de dos números: la entalpía ΔH (el calor que se libera al formarse la hélice) y la entropía ΔS (el orden que se pierde al juntarse las dos hebras). Empiezan en cero y vamos a sumar lo que aporta cada pareja de bases vecinas.",
            ),
    );

    let mut entalpia = 0.0;
    let mut entropia = 0.0;
    let mut filas: Vec<Vec<String>> = Vec::new();
    for (indice, par) in cebador.windows(2).enumerate() {
        let (entalpia_del_par, entropia_del_par) = cebadores::apilamiento(par[0], par[1]);
        entalpia += entalpia_del_par;
        entropia += entropia_del_par;
        filas.push(vec![
            texto(par),
            numero(entalpia_del_par, 1),
            numero(entropia_del_par, 1),
        ]);
        if indice < PARES_DETALLADOS {
            let ultima = filas.len() - 1;
            narracion.agregar(
                Pizarra::nueva()
                    .variable("par", texto(par))
                    .variable("entalpia_del_par", numero(entalpia_del_par, 1))
                    .variable("entropia_del_par", numero(entropia_del_par, 1))
                    .variable("entalpia", numero(entalpia, 1))
                    .variable("entropia", numero(entropia, 1))
                    .secuencia(&cebador, &[(indice, indice + 2, Resalte::Actual)])
                    .tabla(
                        "Pares de bases vecinas (ΔH en kcal/mol, ΔS en cal/K·mol)",
                        &["par", "ΔH", "ΔS"],
                        filas.clone(),
                        &[(ultima, 1), (ultima, 2)],
                    )
                    .paso(
                        "entalpia += entalpia_del_par;",
                        format!(
                            "El par {} aporta ΔH = {} y ΔS = {} (son valores medidos en el laboratorio para cada una de las 10 parejas posibles). Llevamos ΔH = {} y ΔS = {}.",
                            texto(par),
                            numero(entalpia_del_par, 1),
                            numero(entropia_del_par, 1),
                            numero(entalpia, 1),
                            numero(entropia, 1)
                        ),
                    ),
            );
        }
    }
    let pares = longitud - 1;
    narracion.agregar(
        Pizarra::nueva()
            .variable("entalpia", numero(entalpia, 1))
            .variable("entropia", numero(entropia, 1))
            .secuencia(&cebador, &todo)
            .tabla(
                "Pares de bases vecinas (ΔH en kcal/mol, ΔS en cal/K·mol)",
                &["par", "ΔH", "ΔS"],
                filas,
                &[],
            )
            .paso(
                "for par in oligo.windows(2) {",
                format!(
                    "El bucle pasó por los {pares} pares de vecinas del cebador y sumó todos los aportes: ΔH = {} kcal/mol y ΔS = {} cal/K·mol. Lo que sigue corrige los extremos.",
                    numero(entalpia, 1),
                    numero(entropia, 1)
                ),
            ),
    );

    let primera = cebador[0];
    let ultima = cebador[longitud - 1];
    let (entalpia_inicio, entropia_inicio) = cebadores::iniciacion(primera);
    let (entalpia_final, entropia_final) = cebadores::iniciacion(ultima);
    entalpia += entalpia_inicio + entalpia_final;
    entropia += entropia_inicio + entropia_final;
    narracion.agregar(
        Pizarra::nueva()
            .variable("primera", primera as char)
            .variable("ultima", ultima as char)
            .variable("entalpia_inicio", numero(entalpia_inicio, 1))
            .variable("entropia_inicio", numero(entropia_inicio, 1))
            .variable("entalpia_final", numero(entalpia_final, 1))
            .variable("entropia_final", numero(entropia_final, 1))
            .variable("entalpia", numero(entalpia, 1))
            .variable("entropia", numero(entropia, 1))
            .secuencia(
                &cebador,
                &[
                    (0, 1, Resalte::Comparado),
                    (longitud - 1, longitud, Resalte::Comparado),
                ],
            )
            .paso(
                "entalpia += entalpia_inicio + entalpia_final;",
                format!(
                    "Los dos extremos de la hélice cuestan algo según su base: la primera ({}) aporta ΔH = {} y ΔS = {}, y la última ({}) aporta ΔH = {} y ΔS = {}. Las G y C se agarran más fuerte que las A y T. Ahora ΔH = {} y ΔS = {}.",
                    primera as char,
                    numero(entalpia_inicio, 1),
                    numero(entropia_inicio, 1),
                    ultima as char,
                    numero(entalpia_final, 1),
                    numero(entropia_final, 1),
                    numero(entalpia, 1),
                    numero(entropia, 1)
                ),
            ),
    );

    let autocomplementario = cebadores::es_autocomplementario(&cebador);
    narracion.agregar(
        Pizarra::nueva()
            .variable("autocomplementario", autocomplementario)
            .secuencia(&cebador, &todo)
            .paso(
                "let autocomplementario = es_autocomplementario(oligo);",
                if autocomplementario {
                    "Este cebador es su propio complemento inverso: dos copias suyas se pueden pegar entre sí formando una hélice perfecta. Eso cambia la cuenta, porque las dos hebras son iguales."
                } else {
                    "Se revisa si el cebador es su propio complemento inverso (un palíndromo de ADN). No lo es, así que el cálculo sigue con dos hebras distintas."
                },
            ),
    );
    let duplex = cebadores::duplex(&cebador);
    if autocomplementario {
        let simetria = duplex.entropia - entropia;
        entropia += simetria;
        narracion.agregar(
            Pizarra::nueva()
                .variable("entropia", numero(entropia, 1))
                .paso(
                    "entropia += ENTROPIA_DE_SIMETRIA;",
                    format!(
                        "Por ser autocomplementario se suma una pequeña corrección de simetría ({} cal/K·mol). Ahora ΔS = {}.",
                        numero(simetria, 1),
                        numero(entropia, 1)
                    ),
                ),
        );
    }
    let marcas_gc: Vec<(usize, usize, Resalte)> = cebador
        .iter()
        .enumerate()
        .filter(|(_, base)| **base == b'G' || **base == b'C')
        .map(|(lugar, _)| (lugar, lugar + 1, Resalte::Hallado))
        .collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("bases_gc", duplex.bases_gc)
            .variable("longitud", duplex.longitud)
            .secuencia(&cebador, &marcas_gc)
            .paso(
                "bases_gc: contar_gc(oligo),",
                format!(
                    "Se cuentan las G y C: {} de {} bases ({}). Aparecen marcadas. Más G y C significa una hélice más difícil de separar y una Tm más alta.",
                    duplex.bases_gc,
                    duplex.longitud,
                    porcentaje(duplex.bases_gc as f64 / duplex.longitud as f64, 0)
                ),
            ),
    );

    let sodio_equivalente = cebadores::sodio_equivalente(&condiciones);
    let magnesio_libre = cebadores::magnesio_libre(&condiciones);
    let concentracion = concentracion_en_molar(&duplex, &condiciones);
    narracion.agregar(
        Pizarra::nueva()
            .variable("duplex.entalpia", numero(duplex.entalpia, 1))
            .variable("duplex.entropia", numero(duplex.entropia, 1))
            .variable(
                "condiciones.correccion",
                if owczarzy {
                    "Owczarzy"
                } else {
                    "SantaLucia"
                },
            )
            .paso(
                "match condiciones.correccion {",
                format!(
                    "Con ΔH = {} y ΔS = {} ya casi está, pero la Tm también depende de la sal y de cuánto cebador hay. Elegiste la corrección {}, así que el código sigue ese camino.",
                    numero(duplex.entalpia, 1),
                    numero(duplex.entropia, 1),
                    if owczarzy {
                        "de Owczarzy (2008)"
                    } else {
                        "de SantaLucia (1998)"
                    }
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("condiciones.magnesio_mm", decimal(condiciones.magnesio_mm, 1))
            .variable("condiciones.dntp_mm", decimal(condiciones.dntp_mm, 1))
            .variable("magnesio_libre", decimal(magnesio_libre, 2))
            .paso(
                "(condiciones.magnesio_mm - condiciones.dntp_mm).max(0.0)",
                format!(
                    "Los dNTP se agarran al Mg²⁺ y lo dejan sin usar: de {} mM de magnesio solo quedan libres {} mM ({} − {}, y nunca menos de cero).",
                    decimal(condiciones.magnesio_mm, 1),
                    decimal(magnesio_libre, 2),
                    decimal(condiciones.magnesio_mm, 1),
                    decimal(condiciones.dntp_mm, 1)
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("condiciones.monovalentes_mm", decimal(condiciones.monovalentes_mm, 1))
            .variable("magnesio_libre", decimal(magnesio_libre, 2))
            .variable("sodio_equivalente", decimal(sodio_equivalente, 1))
            .paso(
                "condiciones.monovalentes_mm + 120.0 * magnesio_libre(condiciones).sqrt()",
                format!(
                    "El magnesio libre estabiliza la hélice como si hubiera más sodio: {} + 120 × √{} = {} mM de «sodio equivalente».",
                    decimal(condiciones.monovalentes_mm, 1),
                    decimal(magnesio_libre, 2),
                    decimal(sodio_equivalente, 1)
                ),
            ),
    );

    if !owczarzy {
        let sodio = sodio_equivalente / 1000.0;
        let fosfatos = (duplex.longitud.saturating_sub(1)) as f64;
        let entropia_corregida =
            duplex.entropia + fosfatos * cebadores::correccion_de_entropia(&condiciones);
        narracion.agregar(Pizarra::nueva().variable("sodio", decimal(sodio, 4)).paso(
            "let sodio = sodio_equivalente(condiciones) / 1000.0;",
            format!(
                "Se pasa el sodio equivalente de milimolar a molar: {} mM = {} M.",
                decimal(sodio_equivalente, 1),
                decimal(sodio, 4)
            ),
        ));
        narracion.agregar(
            Pizarra::nueva()
                .variable("fosfatos", fosfatos)
                .paso(
                    "let fosfatos = duplex.longitud.saturating_sub(1) as f64;",
                    format!(
                        "Cada unión entre bases es un fosfato con carga negativa, y la sal los apantalla. Un oligo de {} bases tiene {} fosfatos.",
                        duplex.longitud, fosfatos
                    ),
                ),
        );
        narracion.agregar(
            Pizarra::nueva()
                .variable("duplex.entropia", numero(duplex.entropia, 1))
                .variable("fosfatos", fosfatos)
                .variable("sodio", decimal(sodio, 4))
                .variable("entropia", numero(entropia_corregida, 1))
                .paso(
                    "let entropia = duplex.entropia + ENTROPIA_POR_FOSFATO * fosfatos * sodio.ln();",
                    format!(
                        "La sal corrige la entropía: {} + 0.368 × {} × ln({}) = {} cal/K·mol. Con menos sal, la hélice es más inestable y la Tm baja.",
                        numero(duplex.entropia, 1),
                        fosfatos,
                        decimal(sodio, 4),
                        numero(entropia_corregida, 1)
                    ),
                ),
        );
        narracion.agregar(
            Pizarra::nueva()
                .variable("concentracion", format!("{} nM", decimal(concentracion * 1e9, 2)))
                .paso(
                    "let concentracion = concentracion_efectiva(duplex, condiciones);",
                    if duplex.autocomplementario {
                        format!(
                            "Como las dos hebras son iguales, cuenta la concentración completa: {} nM.",
                            decimal(concentracion * 1e9, 2)
                        )
                    } else {
                        format!(
                            "Para dos hebras distintas se usa la cuarta parte de la concentración del cebador: {} nM ÷ 4 = {} nM. Más cebador significa que se encuentran más fácil y la Tm sube un poco.",
                            decimal(condiciones.oligo_nm, 1),
                            decimal(concentracion * 1e9, 2)
                        )
                    },
                ),
        );
        let kelvin = 1000.0 * duplex.entalpia
            / (entropia_corregida + cebadores::CONSTANTE_DE_LOS_GASES * concentracion.ln());
        narracion.agregar(
            Pizarra::nueva()
                .variable("duplex.entalpia", numero(duplex.entalpia, 1))
                .variable("entropia", numero(entropia_corregida, 1))
                .variable("concentracion", format!("{} nM", decimal(concentracion * 1e9, 2)))
                .paso(
                    "1000.0 * duplex.entalpia / (entropia + CONSTANTE_DE_LOS_GASES * concentracion.ln())",
                    format!(
                        "La fórmula es Tm = ΔH ÷ (ΔS + R·ln C), con ΔH en cal y R = 1.987: 1000 × {} ÷ ({} + 1.987 × ln({})) = {} K. Restando 273.15 para pasar a grados: {}.",
                        numero(duplex.entalpia, 1),
                        numero(entropia_corregida, 1),
                        decimal(concentracion * 1e9, 2) + "e-9",
                        decimal(kelvin, 1),
                        grados(kelvin - cebadores::CERO_ABSOLUTO)
                    ),
                ),
        );
    } else {
        let monovalentes = condiciones.monovalentes_mm / 1000.0;
        let magnesio = magnesio_libre / 1000.0;
        let proporcion = if monovalentes > 0.0 {
            magnesio.sqrt() / monovalentes
        } else {
            f64::INFINITY
        };
        let a_un_molar = tm_a_un_molar(&duplex, &condiciones);
        narracion.agregar(
            Pizarra::nueva()
                .variable("duplex.entalpia", numero(duplex.entalpia, 1))
                .variable("duplex.entropia", numero(duplex.entropia, 1))
                .variable("concentracion", format!("{} nM", decimal(concentracion * 1e9, 2)))
                .paso(
                    "1000.0 * duplex.entalpia / (duplex.entropia + CONSTANTE_DE_LOS_GASES * concentracion.ln())",
                    format!(
                        "Primero se calcula la Tm como si la sal fuera 1 molar, la referencia de las tablas: {} K ({}). Después se corrige según tu sal.",
                        decimal(a_un_molar, 1),
                        grados(a_un_molar - cebadores::CERO_ABSOLUTO)
                    ),
                ),
        );
        let (linea, nota) = if proporcion < 0.22 {
            (
                "return correccion_monovalente(fraccion_gc, monovalentes);",
                format!(
                    "Se compara el efecto del magnesio con el del sodio y potasio: √Mg ÷ monovalentes = {}. Es menos de 0.22, así que manda la sal monovalente y la corrección solo la usa a ella.",
                    decimal(proporcion, 2)
                ),
            )
        } else if proporcion < 6.0 {
            (
                "if proporcion < 6.0 {",
                format!(
                    "Se compara el efecto del magnesio con el del sodio y potasio: √Mg ÷ monovalentes = {}. Está entre 0.22 y 6: los dos pesan, y se aplican ajustes extra a tres coeficientes de la fórmula.",
                    decimal(proporcion, 2)
                ),
            )
        } else {
            (
                "a + b * ln_mg",
                format!(
                    "Se compara el efecto del magnesio con el del sodio y potasio: √Mg ÷ monovalentes = {}. Es más de 6: manda el magnesio y se usa su fórmula directa.",
                    if proporcion.is_finite() {
                        decimal(proporcion, 2)
                    } else {
                        "infinito".to_string()
                    }
                ),
            )
        };
        narracion.agregar(
            Pizarra::nueva()
                .variable("monovalentes", decimal(monovalentes, 4))
                .variable("magnesio", decimal(magnesio, 5))
                .variable(
                    "proporcion",
                    if proporcion.is_finite() {
                        decimal(proporcion, 2)
                    } else {
                        "∞".to_string()
                    },
                )
                .paso(linea, nota),
        );
        let inversa = 1.0 / (tm + cebadores::CERO_ABSOLUTO);
        let correccion = inversa - 1.0 / a_un_molar;
        narracion.agregar(
            Pizarra::nueva()
                .variable("correccion_owczarzy", format!("{correccion:.3e}"))
                .variable("inversa", format!("{inversa:.6e}"))
                .paso(
                    "1.0 / tm_en_un_molar(duplex, condiciones) + correccion_owczarzy(duplex, condiciones);",
                    format!(
                        "Owczarzy corrige el inverso de la temperatura en kelvin: 1 ÷ {} K + {:.3e} = {:.6e}. La corrección sale del GC ({}) y del largo del cebador.",
                        decimal(a_un_molar, 1),
                        correccion,
                        inversa,
                        porcentaje(duplex.bases_gc as f64 / duplex.longitud as f64, 0)
                    ),
                ),
        );
        narracion.agregar(
            Pizarra::nueva()
                .variable("inversa", format!("{inversa:.6e}"))
                .paso(
                    "1.0 / inversa - CERO_ABSOLUTO",
                    format!(
                        "Se invierte de nuevo y se restan 273.15 para pasar de kelvin a grados: 1 ÷ {:.6e} − 273.15 = {}.",
                        inversa,
                        grados(tm)
                    ),
                ),
        );
    }

    narracion.agregar(
        Pizarra::nueva()
            .variable("tm", numero(tm, 1))
            .paso(
                "if !tm.is_finite() || tm <= -CERO_ABSOLUTO {",
                format!(
                    "Última revisión: la Tm tiene que ser un número real y mayor que el cero absoluto (−273.15 °C). Con {} pasa sin problema.",
                    grados(tm)
                ),
            ),
    );
    if let (Some(otro), Some(tm_otro)) = (&pareja, tm_de_la_pareja) {
        narracion.agregar(
            Pizarra::nueva()
                .variable("oligo", corto(otro))
                .variable("tm", numero(tm_otro, 1))
                .secuencia(otro, &[(0, otro.len(), Resalte::Comparado)])
                .paso(
                    "let tm = tm_del_duplex(&duplex(oligo), condiciones);",
                    format!(
                        "Para su pareja, {}, se repite exactamente lo mismo con sus propias G, C, A y T. Sale Tm = {}, que difiere {} de la del cebador. Lo ideal es que los dos se parezcan para que se peguen a la misma temperatura.",
                        corto(otro),
                        grados(tm_otro),
                        grados((tm - tm_otro).abs())
                    ),
                ),
        );
    }
    Ok(narracion.terminar(
        Pizarra::nueva()
            .variable("tm", numero(tm, 1))
            .paso(
                "Ok(tm)",
                format!(
                    "Resultado: la Tm del cebador es {}. Es el primer valor de la pantalla. Con las mismas condiciones, el resultado además busca dímeros y horquillas y dibuja la curva de fusión: mira «Qué significa» para saber si esta Tm sirve para tu PCR.",
                    grados(tm)
                ),
            ),
    ))
}

const PESO_DEL_SITIO_EXTRA: f64 = 5.0;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Veredicto {
    Ambigua,
    Repeticion,
    ContenidoGc,
    Pinza,
    Temperatura,
    Aceptada,
}

#[derive(Clone, Copy, Debug)]
struct Revision {
    veredicto: Veredicto,
    gc: Option<f64>,
    pinza: Option<usize>,
    tm: Option<f64>,
}

type Candidato = (usize, usize, Revision);

struct Escaneo {
    descartes: [usize; 5],
    aceptadas: usize,
    candidatos: usize,
    primera: Option<Candidato>,
    ejemplos: [Option<Candidato>; 6],
}

fn cuenta_gc(tramo: &[u8]) -> usize {
    tramo
        .iter()
        .filter(|&&base| base == b'G' || base == b'C')
        .count()
}

fn racha_mas_larga(tramo: &[u8]) -> usize {
    let mut mejor = 0;
    let mut actual = 0;
    for (lugar, base) in tramo.iter().enumerate() {
        if lugar > 0 && tramo[lugar - 1] == *base {
            actual += 1;
        } else {
            actual = 1;
        }
        mejor = mejor.max(actual);
    }
    mejor
}

fn revisar_candidato(
    plantilla: &[u8],
    inicio: usize,
    longitud: usize,
    inversa: bool,
    opciones: &OpcionesDiseno,
) -> Revision {
    let tramo = &plantilla[inicio..inicio + longitud];
    let mut revision = Revision {
        veredicto: Veredicto::Ambigua,
        gc: None,
        pinza: None,
        tm: None,
    };
    if tramo
        .iter()
        .any(|base| !matches!(base, b'A' | b'C' | b'G' | b'T'))
    {
        return revision;
    }
    revision.veredicto = Veredicto::Repeticion;
    if racha_mas_larga(tramo) > opciones.repeticion_maxima {
        return revision;
    }
    let gc = cuenta_gc(tramo) as f64 / longitud as f64;
    revision.gc = Some(gc);
    revision.veredicto = Veredicto::ContenidoGc;
    if gc < opciones.gc_minimo || gc > opciones.gc_maximo {
        return revision;
    }
    let pinza = if inversa {
        cuenta_gc(&tramo[..BASES_DE_LA_PINZA])
    } else {
        cuenta_gc(&tramo[longitud - BASES_DE_LA_PINZA..])
    };
    revision.pinza = Some(pinza);
    revision.veredicto = Veredicto::Pinza;
    if !(PINZA_MINIMA..=PINZA_MAXIMA).contains(&pinza) {
        return revision;
    }
    let tm = cebadores::tm_del_duplex(&cebadores::duplex(tramo), &opciones.condiciones);
    revision.tm = Some(tm);
    revision.veredicto = Veredicto::Temperatura;
    if !(opciones.tm_minima..=opciones.tm_maxima).contains(&tm) {
        return revision;
    }
    revision.veredicto = Veredicto::Aceptada;
    revision
}

fn escanear(plantilla: &[u8], inversa: bool, opciones: &OpcionesDiseno) -> Escaneo {
    let mut escaneo = Escaneo {
        descartes: [0; 5],
        aceptadas: 0,
        candidatos: 0,
        primera: None,
        ejemplos: [None; 6],
    };
    for inicio in 0..plantilla.len() {
        for longitud in opciones.longitud_minima..=opciones.longitud_maxima {
            if inicio + longitud > plantilla.len() {
                break;
            }
            let revision = revisar_candidato(plantilla, inicio, longitud, inversa, opciones);
            let lugar = revision.veredicto as usize;
            escaneo.candidatos += 1;
            if lugar < 5 {
                escaneo.descartes[lugar] += 1;
            } else {
                escaneo.aceptadas += 1;
            }
            if escaneo.primera.is_none() {
                escaneo.primera = Some((inicio, longitud, revision));
            }
            if escaneo.ejemplos[lugar].is_none() {
                escaneo.ejemplos[lugar] = Some((inicio, longitud, revision));
            }
        }
    }
    escaneo
}

fn marcas_del_candidato(
    inicio: usize,
    longitud: usize,
    inversa: bool,
    con_pinza: bool,
) -> Vec<(usize, usize, Resalte)> {
    let fin = inicio + longitud;
    if !con_pinza {
        vec![(inicio, fin, Resalte::Actual)]
    } else if inversa {
        vec![
            (inicio, inicio + BASES_DE_LA_PINZA, Resalte::Comparado),
            (inicio + BASES_DE_LA_PINZA, fin, Resalte::Actual),
        ]
    } else {
        vec![
            (inicio, fin - BASES_DE_LA_PINZA, Resalte::Actual),
            (fin - BASES_DE_LA_PINZA, fin, Resalte::Comparado),
        ]
    }
}

fn describir_candidato(
    plantilla: &[u8],
    candidato: &Candidato,
    inversa: bool,
    opciones: &OpcionesDiseno,
) -> (&'static str, String) {
    let (inicio, longitud, revision) = candidato;
    let tramo = &plantilla[*inicio..inicio + longitud];
    let lado = if inversa { "derecho" } else { "izquierdo" };
    let presentacion = format!(
        "Candidato a cebador {lado}: las bases {}–{} de la plantilla ({longitud} bases), {}.",
        miles(inicio + 1),
        miles(inicio + longitud),
        corto(tramo)
    );
    let extremo = if inversa {
        "las primeras 5 bases de la ventana (serán las últimas del cebador derecho, que se lee al revés)"
    } else {
        "las últimas 5 bases"
    };
    match revision.veredicto {
        Veredicto::Ambigua => (
            "return Err(Motivo::BasesAmbiguas);",
            format!(
                "{presentacion} Tiene una base N (desconocida): no se puede calcular su Tm ni saber si se pegaría bien, así que se descarta al instante y ya no se miran las otras reglas."
            ),
        ),
        Veredicto::Repeticion => (
            "return Err(Motivo::Repeticiones);",
            format!(
                "{presentacion} Contiene una racha de {} bases iguales seguidas y el máximo es {}. Las rachas hacen que el cebador resbale sobre la plantilla, así que se descarta.",
                racha_mas_larga(tramo),
                opciones.repeticion_maxima
            ),
        ),
        Veredicto::ContenidoGc => {
            let gc = revision.gc.unwrap_or(0.0);
            let motivo = if gc < opciones.gc_minimo {
                "pocas G y C se pegan flojo"
            } else {
                "tantas G y C se pegan tan fuerte que se agarra donde no debe"
            };
            (
                "return Err(Motivo::ContenidoGc);",
                format!(
                    "{presentacion} Tiene {} G o C de {longitud}: un GC de {}. Se pide entre {} y {}: {motivo}. Se descarta.",
                    cuenta_gc(tramo),
                    porcentaje(gc, 0),
                    porcentaje(opciones.gc_minimo, 0),
                    porcentaje(opciones.gc_maximo, 0)
                ),
            )
        }
        Veredicto::Pinza => (
            "return Err(Motivo::PinzaGc);",
            format!(
                "{presentacion} Su GC es {}, bien. Pero en {extremo} hay {} G o C, y se piden entre {PINZA_MINIMA} y {PINZA_MAXIMA}: con ninguna el extremo 3' se despega fácil, con demasiadas se pega aunque no coincida. Se descarta.",
                porcentaje(revision.gc.unwrap_or(0.0), 0),
                revision.pinza.unwrap_or(0)
            ),
        ),
        Veredicto::Temperatura => (
            "return Err(Motivo::Temperatura);",
            format!(
                "{presentacion} Pasa GC y pinza, pero su Tm es {} y se pide entre {} y {}. Se descarta.",
                grados(revision.tm.unwrap_or(0.0)),
                grados(opciones.tm_minima),
                grados(opciones.tm_maxima)
            ),
        ),
        Veredicto::Aceptada => (
            "Ok(tm) => aceptadas.push(Ventana {",
            format!(
                "{presentacion} Pasa todo: sin N, sin rachas, GC de {}, {} G o C en {extremo} y Tm de {}. Se guarda en `aceptadas`.",
                porcentaje(revision.gc.unwrap_or(0.0), 0),
                revision.pinza.unwrap_or(0),
                grados(revision.tm.unwrap_or(0.0))
            ),
        ),
    }
}

fn paso_de_candidato(
    plantilla: &[u8],
    candidato: &Candidato,
    inversa: bool,
    opciones: &OpcionesDiseno,
) -> Paso {
    let (inicio, longitud, revision) = candidato;
    let (linea, nota) = describir_candidato(plantilla, candidato, inversa, opciones);
    let mut pizarra = Pizarra::nueva()
        .variable("inicio", inicio)
        .variable("longitud", longitud);
    if let Some(gc) = revision.gc {
        pizarra = pizarra.variable("gc", porcentaje(gc, 0));
    }
    if let Some(pinza) = revision.pinza {
        pizarra = pizarra.variable("pinza", pinza);
    }
    if let Some(tm) = revision.tm {
        pizarra = pizarra.variable("tm", numero(tm, 1));
    }
    pizarra
        .secuencia(
            plantilla,
            &marcas_del_candidato(*inicio, *longitud, inversa, revision.pinza.is_some()),
        )
        .paso(linea, nota)
}

fn filas_de_descartes(
    descartes: &[usize; 5],
    aceptadas: usize,
    opciones: &OpcionesDiseno,
) -> Vec<Vec<String>> {
    let reglas = [
        "Tiene N".to_string(),
        format!("Racha de más de {} iguales", opciones.repeticion_maxima),
        format!(
            "GC fuera de {}–{}",
            porcentaje(opciones.gc_minimo, 0),
            porcentaje(opciones.gc_maximo, 0)
        ),
        format!("Pinza sin {PINZA_MINIMA}–{PINZA_MAXIMA} G/C"),
        format!(
            "Tm fuera de {}–{} °C",
            decimal(opciones.tm_minima, 0),
            decimal(opciones.tm_maxima, 0)
        ),
    ];
    let mut filas: Vec<Vec<String>> = reglas
        .iter()
        .zip(descartes)
        .map(|(regla, cuantos)| vec![regla.clone(), miles(*cuantos)])
        .collect();
    filas.push(vec!["Aceptados".to_string(), miles(aceptadas)]);
    filas
}

fn exceso_de_energia(energia: f64, umbral: f64) -> f64 {
    (umbral - energia).max(0.0)
}

fn terminos_de_penalizacion(cebador: &Cebador, opciones: &OpcionesDiseno) -> [f64; 4] {
    [
        (cebador.tm - opciones.tm_optima).abs(),
        cebador.longitud.abs_diff(opciones.longitud_optima) as f64,
        exceso_de_energia(cebador.autodimero, UMBRAL_DIMERO)
            + exceso_de_energia(cebador.autodimero_extremo_3, UMBRAL_DIMERO_EXTREMO_3)
            + exceso_de_energia(cebador.horquilla, UMBRAL_HORQUILLA),
        PESO_DEL_SITIO_EXTRA * cebador.sitios.unwrap_or(1).saturating_sub(1) as f64,
    ]
}

fn energia_o_cero(valor: f64) -> String {
    if valor == 0.0 {
        "ninguna estable".to_string()
    } else {
        energia(valor)
    }
}

fn texto_de_par(lugar: usize, par: &ParDeCebadores) -> String {
    format!(
        "{lugar}. {} + {} · {} pb · penalización {}",
        corto(&par.izquierdo.secuencia),
        corto(&par.derecho.secuencia),
        miles(par.producto),
        decimal(par.penalizacion, 2)
    )
}

fn redondear_penalizacion(valor: f64) -> f64 {
    (valor * 10_000.0).round() / 10_000.0
}

fn pasos_diseno(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let plantilla = leer_plantilla(entrada)?;
    let opciones = leer_opciones(entrada, &plantilla)?;
    let origen = origen_en_ecoli(&entrada.texto("plantilla")?, &plantilla);
    let genoma = origen.as_ref().map(|_| ecoli::genoma());
    let diseno = cebadores::disenar_cebadores(&plantilla, &opciones, genoma)?;
    let izquierda = escanear(&plantilla, false, &opciones);
    let derecha = escanear(&plantilla, true, &opciones);
    let mut narracion = Narracion::nueva(400);
    let objetivo = match opciones.objetivo {
        Some((inicio, longitud)) => {
            format!("bases {}–{}", miles(inicio + 1), miles(inicio + longitud))
        }
        None => "ninguno".to_string(),
    };
    let marcas_del_objetivo: Vec<(usize, usize, Resalte)> = match opciones.objetivo {
        Some((inicio, longitud)) => vec![(inicio, inicio + longitud, Resalte::Hallado)],
        None => Vec::new(),
    };

    narracion.agregar(
        Pizarra::nueva()
            .variable("plantilla", format!("{} pb", miles(plantilla.len())))
            .variable("opciones.tm_minima", decimal(opciones.tm_minima, 1))
            .variable("opciones.tm_optima", decimal(opciones.tm_optima, 1))
            .variable("opciones.tm_maxima", decimal(opciones.tm_maxima, 1))
            .variable("opciones.producto_minimo", opciones.producto_minimo)
            .variable("opciones.producto_maximo", opciones.producto_maximo)
            .variable("opciones.objetivo", &objetivo)
            .secuencia(&plantilla, &marcas_del_objetivo)
            .paso(
                "validar_diseno(plantilla, opciones)?;",
                format!(
                    "Tienes una plantilla de {} pb y quieres amplificar un producto de {} a {} pb con cebadores de Tm cercana a {} (se aceptan de {} a {}). Región objetivo: {objetivo}. Primero el código revisa que la petición tenga sentido: la plantilla cabe (máximo {} pb), el producto cabe (máximo {} pb) y las Tm y el GC están en orden. Todo en regla.",
                    miles(plantilla.len()),
                    miles(opciones.producto_minimo),
                    miles(opciones.producto_maximo),
                    grados(opciones.tm_optima),
                    grados(opciones.tm_minima),
                    grados(opciones.tm_maxima),
                    miles(PLANTILLA_MAXIMA),
                    miles(PRODUCTO_MAXIMO)
                ),
            ),
    );

    let mut suma_de_gc = vec![0usize; plantilla.len() + 1];
    for (lugar, base) in plantilla.iter().enumerate() {
        suma_de_gc[lugar + 1] = suma_de_gc[lugar] + usize::from(*base == b'G' || *base == b'C');
    }
    let filas_de_sumas: Vec<Vec<String>> = plantilla
        .iter()
        .take(8)
        .enumerate()
        .map(|(lugar, base)| {
            vec![
                miles(lugar + 1),
                (*base as char).to_string(),
                suma_de_gc[lugar + 1].to_string(),
            ]
        })
        .collect();
    let mut pizarra_de_sumas = Pizarra::nueva();
    let mut nota_de_sumas = "Para revisar miles de candidatos sin contar base por base cada vez, el código prepara «sumas acumuladas»: en cada posición guarda cuántas G o C hay desde el principio (y lo mismo para las N, las rachas y la energía de cada par de vecinas). Así el GC de cualquier ventana sale restando dos números.".to_string();
    if let Some((inicio, longitud, _)) = &izquierda.primera {
        let fin = inicio + longitud;
        pizarra_de_sumas = pizarra_de_sumas
            .variable("inicio", inicio)
            .variable("fin", fin)
            .variable("sumas.gc[inicio]", suma_de_gc[*inicio])
            .variable("sumas.gc[fin]", suma_de_gc[fin]);
        nota_de_sumas.push_str(&format!(
            " Por ejemplo, la ventana {}–{}: {} − {} = {} G o C.",
            miles(inicio + 1),
            miles(fin),
            suma_de_gc[fin],
            suma_de_gc[*inicio],
            suma_de_gc[fin] - suma_de_gc[*inicio]
        ));
    }
    narracion.agregar(
        pizarra_de_sumas
            .secuencia(&plantilla, &[(0, 8.min(plantilla.len()), Resalte::Actual)])
            .tabla(
                "sumas.gc: cuántas G o C hay hasta cada base (primeras 8)",
                &["posición", "base", "sumas.gc"],
                filas_de_sumas,
                &[],
            )
            .paso(
                "let sumas = sumas_acumuladas(plantilla, opciones.repeticion_maxima);",
                nota_de_sumas,
            ),
    );

    narracion.agregar(
        Pizarra::nueva()
            .variable("opciones.longitud_minima", opciones.longitud_minima)
            .variable("opciones.longitud_maxima", opciones.longitud_maxima)
            .variable("candidatos", miles(izquierda.candidatos))
            .paso(
                "buscar_ventanas(plantilla, &sumas, Hebra::Directa, opciones, &mut descartes);",
                format!(
                    "Primero los candidatos a cebador izquierdo, que se pega a la hebra tal como está escrita. El código prueba cada posición de inicio con cada longitud de {} a {} bases: {} candidatos en total. A cada uno le pasa cuatro reglas en orden y se descarta en la primera que falla.",
                    opciones.longitud_minima,
                    opciones.longitud_maxima,
                    miles(izquierda.candidatos)
                ),
            ),
    );
    let mut mostrados: Vec<(usize, usize)> = Vec::new();
    let mut ejemplos: Vec<Candidato> = Vec::new();
    if let Some(primera) = izquierda.primera {
        mostrados.push((primera.0, primera.1));
        ejemplos.push(primera);
    }
    for ejemplo in izquierda.ejemplos.iter().flatten() {
        if !mostrados.contains(&(ejemplo.0, ejemplo.1)) {
            mostrados.push((ejemplo.0, ejemplo.1));
            ejemplos.push(*ejemplo);
        }
    }
    for ejemplo in &ejemplos {
        narracion.agregar(paso_de_candidato(&plantilla, ejemplo, false, &opciones));
    }
    narracion.agregar(
        Pizarra::nueva()
            .variable("ventanas_izquierdas.len()", miles(izquierda.aceptadas))
            .tabla(
                "Candidatos izquierdos: qué pasó con cada uno",
                &["regla", "candidatos"],
                filas_de_descartes(&izquierda.descartes, izquierda.aceptadas, &opciones),
                &[(5, 0), (5, 1)],
            )
            .paso(
                "let ventanas_izquierdas =",
                format!(
                    "Terminó el recorrido de la hebra directa. De {} candidatos sobrevivieron {}: esos son los cebadores izquierdos posibles.",
                    miles(izquierda.candidatos),
                    miles(izquierda.aceptadas)
                ),
            ),
    );

    narracion.agregar(
        Pizarra::nueva()
            .variable("candidatos", miles(derecha.candidatos))
            .paso(
                "buscar_ventanas(plantilla, &sumas, Hebra::Inversa, opciones, &mut descartes);",
                "Ahora lo mismo para el cebador derecho, que se pega a la hebra de enfrente y se escribe como el complemento inverso del tramo. Se recorren los mismos tramos con las mismas reglas, con una diferencia: la pinza GC se revisa en las primeras 5 bases de la ventana, porque esas serán el extremo 3' del cebador derecho.",
            ),
    );
    let ejemplo_derecho = derecha.ejemplos[Veredicto::Aceptada as usize].or(derecha.primera);
    if let Some(candidato) = &ejemplo_derecho {
        narracion.agregar(paso_de_candidato(&plantilla, candidato, true, &opciones));
    }
    narracion.agregar(
        Pizarra::nueva()
            .variable("ventanas_derechas.len()", miles(derecha.aceptadas))
            .tabla(
                "Candidatos derechos: qué pasó con cada uno",
                &["regla", "candidatos"],
                filas_de_descartes(&derecha.descartes, derecha.aceptadas, &opciones),
                &[(5, 0), (5, 1)],
            )
            .paso(
                "let ventanas_derechas =",
                format!(
                    "Terminó la hebra inversa. De {} candidatos sobrevivieron {} cebadores derechos posibles. En total, las dos hebras sumaron {} descartes y {} aceptados.",
                    miles(derecha.candidatos),
                    miles(derecha.aceptadas),
                    miles(diseno.descartes.bases_ambiguas
                        + diseno.descartes.repeticiones
                        + diseno.descartes.contenido_gc
                        + diseno.descartes.pinza_gc
                        + diseno.descartes.temperatura),
                    miles(diseno.izquierdos + diseno.derechos)
                ),
            ),
    );

    let mejor = diseno.pares.first();
    if let Some(par) = mejor {
        let izquierdo = &par.izquierdo;
        let derecho = &par.derecho;
        let terminos = terminos_de_penalizacion(izquierdo, &opciones);
        let sin_sitios = terminos[0] + terminos[1] + terminos[2];
        let marcas_izquierdo = [(
            izquierdo.inicio,
            izquierdo.inicio + izquierdo.longitud,
            Resalte::Actual,
        )];
        narracion.agregar(
            Pizarra::nueva()
                .variable("ventanas_izquierdas.len()", miles(izquierda.aceptadas))
                .secuencia(&plantilla, &marcas_izquierdo)
                .paso(
                    "construir_cebadores(plantilla, &ventanas_izquierdas, Hebra::Directa, opciones);",
                    format!(
                        "Cada ventana aceptada se convierte en un cebador completo, al que se le mide el autodímero, la horquilla y una penalización total. Hay miles, así que seguimos al cebador izquierdo del mejor par: bases {}–{}, {}, Tm {}, GC {}.",
                        miles(izquierdo.inicio + 1),
                        miles(izquierdo.inicio + izquierdo.longitud),
                        con_extremos(&izquierdo.secuencia),
                        grados(izquierdo.tm),
                        porcentaje(izquierdo.gc, 0)
                    ),
                ),
        );
        narracion.agregar(
            Pizarra::nueva()
                .variable("secuencia", corto(&izquierdo.secuencia))
                .variable("propio.mas_estable.energia", numero(izquierdo.autodimero, 2))
                .variable(
                    "propio.extremo_3.energia",
                    numero(izquierdo.autodimero_extremo_3, 2),
                )
                .paso(
                    "let propio = buscar_dimero(&secuencia, &secuencia, sal);",
                    format!(
                        "Se acerca el cebador a una copia de sí mismo y se busca cómo se pegarían entre ellas (un autodímero, que gasta cebador sin amplificar nada). El apareamiento más estable: {}; el que toca un extremo 3': {}. Se toleran hasta {} y {} kcal/mol.",
                        energia_o_cero(izquierdo.autodimero),
                        energia_o_cero(izquierdo.autodimero_extremo_3),
                        numero(UMBRAL_DIMERO, 0),
                        numero(UMBRAL_DIMERO_EXTREMO_3, 0)
                    ),
                ),
        );
        narracion.agregar(
            Pizarra::nueva()
                .variable("pliegue.energia", numero(izquierdo.horquilla, 2))
                .paso(
                    "let pliegue = buscar_horquilla(&secuencia, sal);",
                    format!(
                        "Ahora se busca si el cebador puede doblarse sobre sí mismo como una horquilla, con un tramo pegado a otro del mismo cebador. Mejor horquilla: {}. Se tolera hasta {} kcal/mol.",
                        energia_o_cero(izquierdo.horquilla),
                        numero(UMBRAL_HORQUILLA, 0)
                    ),
                ),
        );
        narracion.agregar(
            Pizarra::nueva()
                .variable("cebador.penalizacion", decimal(sin_sitios, 2))
                .tabla(
                    "De qué se compone la penalización (menos es mejor)",
                    &["término", "valor"],
                    vec![
                        vec![
                            format!(
                                "Tm lejos de la óptima: |{} − {}|",
                                decimal(izquierdo.tm, 1),
                                decimal(opciones.tm_optima, 1)
                            ),
                            decimal(terminos[0], 2),
                        ],
                        vec![
                            format!(
                                "Longitud lejos de la óptima: |{} − {}|",
                                izquierdo.longitud, opciones.longitud_optima
                            ),
                            decimal(terminos[1], 2),
                        ],
                        vec![
                            "Estructura: cuánto pasa de los umbrales".to_string(),
                            decimal(terminos[2], 2),
                        ],
                        vec!["Total por ahora".to_string(), decimal(sin_sitios, 2)],
                    ],
                    &[(3, 1)],
                )
                .paso(
                    "cebador.penalizacion = penalizar_cebador(&cebador, opciones);",
                    format!(
                        "La penalización suma lo que se aleja el cebador de lo ideal: Tm, longitud y estructura. Este lleva {}. Cuanto menor, mejor candidato.",
                        decimal(sin_sitios, 2)
                    ),
                ),
        );
        let tramo_derecho = &plantilla[derecho.inicio..derecho.inicio + derecho.longitud];
        narracion.agregar(
            Pizarra::nueva()
                .variable("tramo", corto(tramo_derecho))
                .variable("secuencia", corto(&derecho.secuencia))
                .lista(
                    "El cebador derecho se escribe al revés y complementado",
                    vec![
                        format!("tramo de la plantilla: {}", corto(tramo_derecho)),
                        format!("cebador derecho: {}", con_extremos(&derecho.secuencia)),
                    ],
                    &[1],
                )
                .secuencia(
                    &plantilla,
                    &[(derecho.inicio, derecho.inicio + derecho.longitud, Resalte::Actual)],
                )
                .paso(
                    "Hebra::Inversa => complemento_inverso(tramo),",
                    format!(
                        "Con los derechos se hace lo mismo, pero antes se invierte y complementa el tramo: la polimerasa solo copia de 5' a 3', así que el cebador de enfrente se escribe leyendo la otra hebra. El mejor par usa el derecho de las bases {}–{}, Tm {}.",
                        miles(derecho.inicio + 1),
                        miles(derecho.inicio + derecho.longitud),
                        grados(derecho.tm)
                    ),
                ),
        );
        if let (Some(sitios_izquierdo), Some(sitios_derecho)) = (izquierdo.sitios, derecho.sitios) {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("extremos.len()", miles(diseno.izquierdos + diseno.derechos))
                    .variable("sitios (izquierdo)", sitios_izquierdo)
                    .variable("sitios (derecho)", sitios_derecho)
                    .paso(
                        "let sitios = contar_sitios(genoma, &extremos);",
                        format!(
                            "Como la plantilla viene de E. coli, se comprueba la especificidad: se toman los últimos {BASES_DEL_EXTREMO_3} nt de cada uno de los {} cebadores y se cuentan sus apariciones en las dos hebras del genoma ({} pb). El izquierdo del mejor par aparece {} {}, el derecho {} {}. Un solo sitio significa que solo amplificará donde quieres.",
                            miles(diseno.izquierdos + diseno.derechos),
                            miles(ecoli::genoma().len()),
                            sitios_izquierdo,
                            if sitios_izquierdo == 1 { "vez" } else { "veces" },
                            sitios_derecho,
                            if sitios_derecho == 1 { "vez" } else { "veces" }
                        ),
                    ),
            );
            let total = terminos[0] + terminos[1] + terminos[2] + terminos[3];
            narracion.agregar(
                Pizarra::nueva()
                    .variable("cebador.penalizacion", decimal(total, 2))
                    .tabla(
                        "Penalización del cebador izquierdo, ahora con los sitios",
                        &["término", "valor"],
                        vec![
                            vec!["Antes (Tm, longitud, estructura)".to_string(), decimal(sin_sitios, 2)],
                            vec![
                                format!("Sitios de más: 5 × ({sitios_izquierdo} − 1)"),
                                decimal(terminos[3], 2),
                            ],
                            vec!["Total".to_string(), decimal(total, 2)],
                        ],
                        &[(2, 1)],
                    )
                    .paso(
                        "cebador.penalizacion = penalizar_cebador(cebador, opciones);",
                        if terminos[3] > 0.0 {
                            format!(
                                "Cada sitio de más suma 5 puntos de penalización, porque podría amplificar otra región. La penalización de este cebador sube de {} a {}.",
                                decimal(sin_sitios, 2),
                                decimal(total, 2)
                            )
                        } else {
                            format!(
                                "Con un solo sitio no se suma nada: la penalización sigue en {}.",
                                decimal(total, 2)
                            )
                        },
                    ),
            );
        }

        let fin_izquierdo = izquierdo.inicio + izquierdo.longitud;
        let producto = par.producto;
        let mut condiciones_del_par = vec![
            vec![
                format!(
                    "El derecho empieza ({}) después de que termina el izquierdo ({})",
                    miles(derecho.inicio),
                    miles(fin_izquierdo)
                ),
                "sí".to_string(),
            ],
            vec![
                format!(
                    "Producto = {} + {} − {}, entre {} y {} pb",
                    miles(derecho.inicio),
                    derecho.longitud,
                    miles(izquierdo.inicio),
                    miles(opciones.producto_minimo),
                    miles(opciones.producto_maximo)
                ),
                format!("{} pb", miles(producto)),
            ],
        ];
        if let Some((inicio, longitud)) = opciones.objetivo {
            condiciones_del_par.push(vec![
                format!(
                    "Rodea la región objetivo {}–{}",
                    miles(inicio + 1),
                    miles(inicio + longitud)
                ),
                "sí".to_string(),
            ]);
        }
        condiciones_del_par.push(vec![
            format!(
                "Tm parecidas: |{} − {}| ≤ {}",
                decimal(izquierdo.tm, 1),
                decimal(derecho.tm, 1),
                decimal(opciones.diferencia_tm_maxima, 1)
            ),
            grados(par.diferencia_tm),
        ]);
        let filas_condiciones = condiciones_del_par.len();
        narracion.agregar(
            Pizarra::nueva()
                .variable("izquierdos", miles(diseno.izquierdos))
                .variable("derechos", miles(diseno.derechos))
                .variable("opciones.pares", opciones.pares)
                .paso(
                    "let pares = emparejar(&izquierdos, &derechos, opciones)?;",
                    format!(
                        "Ahora hay que formar parejas. Con {} izquierdos y {} derechos habría {} combinaciones: demasiadas para revisar una por una. El código ordena los izquierdos de menor a mayor penalización, solo mira los derechos que quedan a distancia de producto razonable y mantiene una lista con los {} mejores pares.",
                        miles(diseno.izquierdos),
                        miles(diseno.derechos),
                        miles((diseno.izquierdos as f64 * diseno.derechos as f64) as usize),
                        opciones.pares
                    ),
                ),
        );
        narracion.agregar(
            Pizarra::nueva()
                .variable("fin_izquierdo", fin_izquierdo)
                .variable("producto", producto)
                .variable("izquierdo.tm", numero(izquierdo.tm, 1))
                .variable("derecho.tm", numero(derecho.tm, 1))
                .secuencia(
                    &plantilla,
                    &[
                        (izquierdo.inicio, fin_izquierdo, Resalte::Actual),
                        (
                            derecho.inicio,
                            derecho.inicio + derecho.longitud,
                            Resalte::Comparado,
                        ),
                    ],
                )
                .tabla(
                    "Condiciones para que dos cebadores formen un par válido",
                    &["condición", "valor"],
                    condiciones_del_par,
                    &(0..filas_condiciones).map(|fila| (fila, 1)).collect::<Vec<_>>(),
                )
                .paso(
                    "let Some(producto) = producto_valido(izquierdo, derecho, opciones) else {",
                    format!(
                        "Mirando al mejor par: el izquierdo copia hacia la derecha y el derecho hacia la izquierda, y el producto es todo lo que queda entre el inicio del primero y el final del segundo: {} pb. Cumple las condiciones, así que se evalúa a fondo.",
                        miles(producto)
                    ),
                ),
        );
        narracion.agregar(
            Pizarra::nueva()
                .variable("cruzado.mas_estable.energia", numero(par.dimero, 2))
                .variable("cruzado.extremo_3.energia", numero(par.dimero_extremo_3, 2))
                .paso(
                    "let cruzado = buscar_dimero(&izquierdo.secuencia, &derecho.secuencia, sal);",
                    format!(
                        "Se acercan los dos cebadores entre sí para ver si se pegarían (un dímero de cebadores). Apareamiento más estable: {}; el que toca un extremo 3': {}. Si se pegan por el extremo 3', la polimerasa los alarga y salen bandas cortas inútiles en el gel.",
                        energia_o_cero(par.dimero),
                        energia_o_cero(par.dimero_extremo_3)
                    ),
                ),
        );
        let exceso_cruzado = exceso_de_energia(par.dimero, UMBRAL_DIMERO)
            + exceso_de_energia(par.dimero_extremo_3, UMBRAL_DIMERO_EXTREMO_3);
        let suma_de_par = izquierdo.penalizacion + derecho.penalizacion + par.diferencia_tm;
        narracion.agregar(
            Pizarra::nueva()
                .variable("diferencia_tm", decimal(par.diferencia_tm, 2))
                .variable("par.penalizacion", decimal(par.penalizacion, 2))
                .tabla(
                    "Penalización del par",
                    &["término", "valor"],
                    vec![
                        vec!["Cebador izquierdo".to_string(), decimal(izquierdo.penalizacion, 2)],
                        vec!["Cebador derecho".to_string(), decimal(derecho.penalizacion, 2)],
                        vec!["Diferencia de Tm".to_string(), decimal(par.diferencia_tm, 2)],
                        vec!["Dímero entre los dos".to_string(), decimal(exceso_cruzado, 2)],
                        vec!["Total".to_string(), decimal(par.penalizacion, 2)],
                    ],
                    &[(4, 1)],
                )
                .paso(
                    "penalizacion_sin_dimero(izquierdo, derecho) + PESO_ESTRUCTURA * estructura,",
                    format!(
                        "La penalización del par suma la de cada cebador ({} y {}), la diferencia de Tm ({}) y lo que el dímero entre ambos pase de los umbrales ({}): {}.",
                        decimal(izquierdo.penalizacion, 2),
                        decimal(derecho.penalizacion, 2),
                        decimal(par.diferencia_tm, 2),
                        decimal(exceso_cruzado, 2),
                        decimal(redondear_penalizacion(suma_de_par + exceso_cruzado), 2)
                    ),
                ),
        );
        let textos: Vec<String> = diseno
            .pares
            .iter()
            .enumerate()
            .map(|(lugar, otro)| texto_de_par(lugar + 1, otro))
            .collect();
        narracion.agregar(
            Pizarra::nueva()
                .variable("mejores.len()", diseno.pares.len())
                .variable("opciones.pares", opciones.pares)
                .lista("mejores: los pares con menor penalización", textos, &[0])
                .paso(
                    "guardar_si_es_de_los_mejores(&mut mejores, par, opciones.pares);",
                    format!(
                        "El par entra en la lista `mejores`, que se mantiene ordenada por penalización y nunca pasa de {} pares: si llega uno mejor, el peor sale. Así queda la lista al terminar la búsqueda.",
                        opciones.pares
                    ),
                ),
        );
        let peor = diseno
            .pares
            .last()
            .map_or(f64::INFINITY, |otro| otro.penalizacion);
        narracion.agregar(
            Pizarra::nueva()
                .variable(
                    "peor_penalizacion",
                    if diseno.pares.len() < opciones.pares {
                        "∞".to_string()
                    } else {
                        decimal(peor, 2)
                    },
                )
                .paso(
                    "if mejor_posible > peor_penalizacion(&mejores, opciones.pares) {",
                    if diseno.pares.len() < opciones.pares {
                        "Atajo del código: cuando la lista está llena y el mejor par que aún podría salir de un izquierdo ya es peor que el último guardado, se deja de buscar. Aquí la lista no llegó a llenarse, así que se revisaron todos los pares posibles.".to_string()
                    } else {
                        format!(
                            "Atajo del código: la lista ya está llena y el último par guardado tiene penalización {}. Cuando el mejor par que aún podría salir de un izquierdo pasa de ese valor, se deja de buscar: todos los izquierdos que faltan son peores. Por eso no se revisan millones de parejas.",
                            decimal(peor, 2)
                        )
                    },
                ),
        );
    } else {
        narracion.agregar(
            Pizarra::nueva()
                .variable("izquierdos", miles(diseno.izquierdos))
                .variable("derechos", miles(diseno.derechos))
                .paso(
                    "let pares = emparejar(&izquierdos, &derechos, opciones)?;",
                    format!(
                        "Ahora se intentan formar parejas con {} izquierdos y {} derechos que sobrevivieron. Ninguna combinación da un producto de {}–{} pb{} con Tm que difieran {} o menos.",
                        miles(diseno.izquierdos),
                        miles(diseno.derechos),
                        miles(opciones.producto_minimo),
                        miles(opciones.producto_maximo),
                        if opciones.objetivo.is_some() {
                            " que rodee la región objetivo"
                        } else {
                            ""
                        },
                        grados(opciones.diferencia_tm_maxima)
                    ),
                ),
        );
    }

    let filas_finales: Vec<Vec<String>> = diseno
        .pares
        .iter()
        .enumerate()
        .map(|(lugar, par)| {
            vec![
                (lugar + 1).to_string(),
                texto(&par.izquierdo.secuencia),
                texto(&par.derecho.secuencia),
                miles(par.producto),
                decimal(par.penalizacion, 2),
            ]
        })
        .collect();
    let marcas_finales: Vec<(usize, usize, Resalte)> = match mejor {
        Some(par) => vec![
            (
                par.izquierdo.inicio,
                par.izquierdo.inicio + par.izquierdo.longitud,
                Resalte::Hallado,
            ),
            (
                par.derecho.inicio,
                par.derecho.inicio + par.derecho.longitud,
                Resalte::Hallado,
            ),
        ],
        None => Vec::new(),
    };
    let nota_final = match mejor {
        Some(par) => format!(
            "Resultado: el mejor par es {} + {}, con un producto de {} pb y penalización {}. Salieron {} pares en total. Las Tm y los dímeros ya están revisados: en «Qué significa» te digo si este par sirve para tu PCR.",
            con_extremos(&par.izquierdo.secuencia),
            con_extremos(&par.derecho.secuencia),
            miles(par.producto),
            decimal(par.penalizacion, 2),
            diseno.pares.len()
        ),
        None => "Resultado: no hay ningún par que cumpla todo. Mira «Qué significa» para saber qué regla descartó más candidatos y qué cambiar (la Tm, el rango del producto o la región objetivo).".to_string(),
    };
    Ok(narracion.terminar(
        Pizarra::nueva()
            .variable("pares.len()", diseno.pares.len())
            .secuencia(&plantilla, &marcas_finales)
            .tabla(
                "Los mejores pares de cebadores",
                &["#", "izquierdo", "derecho", "producto (pb)", "penalización"],
                filas_finales,
                &if mejor.is_some() {
                    vec![(0, 1), (0, 2)]
                } else {
                    Vec::new()
                },
            )
            .paso("Ok(Diseno {", nota_final),
    ))
}

#[cfg(test)]
mod pruebas {
    use super::*;
    use crate::{codigo, ejemplo, pasos};
    use serde_json::json;

    fn entrada_de(id: &str) -> serde_json::Value {
        ejemplo(id).unwrap()
    }

    fn comprobar_lineas(id: &str, valor: serde_json::Value) -> Vec<Paso> {
        let lista = pasos(id, valor).unwrap();
        assert!(!lista.is_empty() && lista.len() <= 400);
        let fuente: String = codigo(id)
            .unwrap()
            .iter()
            .map(|fragmento| fragmento.fuente.as_str())
            .collect::<Vec<_>>()
            .join("\n");
        for paso in &lista {
            assert!(
                fuente.lines().any(|linea| linea.contains(&paso.linea)),
                "{}: línea ausente: {}",
                id,
                paso.linea
            );
            assert!(!paso.nota.is_empty());
        }
        lista
    }

    fn con_campos(id: &str, cambios: &[(&str, &str)]) -> serde_json::Value {
        let mut valor = entrada_de(id);
        for (campo, nuevo) in cambios {
            valor[*campo] = json!(nuevo);
        }
        valor
    }

    #[test]
    fn la_temperatura_narrada_coincide_con_la_del_nucleo() {
        let id = "puerto.temperatura_fusion";
        let casos = [
            con_campos(id, &[]),
            con_campos(id, &[("correccion", CORRECCIONES[1])]),
            con_campos(id, &[("pareja", "")]),
            con_campos(id, &[("cebador", "CGCGAATTCGCG"), ("pareja", "")]),
            con_campos(
                id,
                &[
                    ("correccion", CORRECCIONES[1]),
                    ("cebador", "CGCGAATTCGCG"),
                    ("magnesio", "5"),
                    ("monovalentes", "10"),
                ],
            ),
            con_campos(id, &[("correccion", CORRECCIONES[1]), ("magnesio", "0.2")]),
            con_campos(
                id,
                &[("correccion", CORRECCIONES[1]), ("monovalentes", "0")],
            ),
        ];
        for caso in casos {
            let lista = comprobar_lineas(id, caso.clone());
            let entrada = Entrada::desde_json(caso).unwrap();
            let cebador = entrada.adn("cebador").unwrap();
            let condiciones = leer_condiciones(&entrada).unwrap();
            let tm = cebadores::temperatura_fusion(&cebador, &condiciones).unwrap();
            let ultimo = lista.last().unwrap();
            assert_eq!(ultimo.linea, "Ok(tm)");
            assert!(ultimo.nota.contains(&grados(tm)));
            assert_eq!(ultimo.estado["variables"][0][1], numero(tm, 1));
        }
    }

    #[test]
    fn los_pasos_de_energia_suman_lo_mismo_que_el_nucleo() {
        let cebador = b"TGTGCCGAAATGGTCCATCA";
        let duplex = cebadores::duplex(cebador);
        let mut entalpia = 0.0;
        let mut entropia = 0.0;
        for par in cebador.windows(2) {
            let (h, s) = cebadores::apilamiento(par[0], par[1]);
            entalpia += h;
            entropia += s;
        }
        let (hi, si) = cebadores::iniciacion(cebador[0]);
        let (hf, sf) = cebadores::iniciacion(cebador[cebador.len() - 1]);
        assert!((entalpia + hi + hf - duplex.entalpia).abs() < 1e-9);
        assert!((entropia + si + sf - duplex.entropia).abs() < 1e-9);
        let condiciones = Condiciones::default();
        let sodio = cebadores::sodio_equivalente(&condiciones) / 1000.0;
        let corregida = duplex.entropia
            + (cebador.len() - 1) as f64 * cebadores::correccion_de_entropia(&condiciones);
        let concentracion = condiciones.oligo_nm * 1e-9 / 4.0;
        let kelvin = 1000.0 * duplex.entalpia
            / (corregida + cebadores::CONSTANTE_DE_LOS_GASES * concentracion.ln());
        let tm = cebadores::temperatura_fusion(cebador, &condiciones).unwrap();
        assert!(sodio > 0.0);
        assert!((kelvin - cebadores::CERO_ABSOLUTO - tm).abs() < 1e-9);
    }

    #[test]
    fn la_correccion_de_owczarzy_se_deduce_del_resultado() {
        let cebador = b"TGTGCCGAAATGGTCCATCA";
        let condiciones = Condiciones {
            correccion: CorreccionSal::Owczarzy,
            ..Condiciones::default()
        };
        let duplex = cebadores::duplex(cebador);
        let tm = cebadores::temperatura_fusion(cebador, &condiciones).unwrap();
        let a_un_molar = tm_a_un_molar(&duplex, &condiciones);
        let correccion = 1.0 / (tm + cebadores::CERO_ABSOLUTO) - 1.0 / a_un_molar;
        let reconstruida = 1.0 / (1.0 / a_un_molar + correccion) - cebadores::CERO_ABSOLUTO;
        assert!((reconstruida - tm).abs() < 1e-9);
        assert!(correccion.abs() < 1e-3);
    }

    fn diseno_de(entrada: &Entrada) -> (Vec<u8>, OpcionesDiseno, Diseno) {
        let plantilla = leer_plantilla(entrada).unwrap();
        let opciones = leer_opciones(entrada, &plantilla).unwrap();
        let origen = origen_en_ecoli(&entrada.texto("plantilla").unwrap(), &plantilla);
        let genoma = origen.as_ref().map(|_| ecoli::genoma());
        let diseno = cebadores::disenar_cebadores(&plantilla, &opciones, genoma).unwrap();
        (plantilla, opciones, diseno)
    }

    #[test]
    fn el_recorrido_de_ventanas_repite_los_descartes_del_nucleo() {
        let id = "puerto.disenar_cebadores";
        let casos = [
            entrada_de(id),
            con_campos(id, &[("objetivo", "")]),
            con_campos(id, &[("plantilla", "@gen:araC"), ("objetivo", "")]),
        ];
        for caso in casos {
            let entrada = Entrada::desde_json(caso).unwrap();
            let (plantilla, opciones, diseno) = diseno_de(&entrada);
            let izquierda = escanear(&plantilla, false, &opciones);
            let derecha = escanear(&plantilla, true, &opciones);
            assert_eq!(izquierda.aceptadas, diseno.izquierdos);
            assert_eq!(derecha.aceptadas, diseno.derechos);
            let sumados: Vec<usize> = (0..5)
                .map(|lugar| izquierda.descartes[lugar] + derecha.descartes[lugar])
                .collect();
            assert_eq!(
                sumados,
                vec![
                    diseno.descartes.bases_ambiguas,
                    diseno.descartes.repeticiones,
                    diseno.descartes.contenido_gc,
                    diseno.descartes.pinza_gc,
                    diseno.descartes.temperatura
                ]
            );
        }
    }

    #[test]
    fn la_penalizacion_desglosada_suma_la_del_nucleo() {
        let id = "puerto.disenar_cebadores";
        let entrada = Entrada::desde_json(entrada_de(id)).unwrap();
        let (_, opciones, diseno) = diseno_de(&entrada);
        let par = diseno.pares.first().unwrap();
        for cebador in [&par.izquierdo, &par.derecho] {
            let total: f64 = terminos_de_penalizacion(cebador, &opciones).iter().sum();
            assert!((total - cebador.penalizacion).abs() < 1e-9);
        }
        let exceso_cruzado = exceso_de_energia(par.dimero, UMBRAL_DIMERO)
            + exceso_de_energia(par.dimero_extremo_3, UMBRAL_DIMERO_EXTREMO_3);
        let suma = par.izquierdo.penalizacion
            + par.derecho.penalizacion
            + par.diferencia_tm
            + exceso_cruzado;
        assert!((redondear_penalizacion(suma) - par.penalizacion).abs() < 1e-9);
    }

    #[test]
    fn el_diseno_narrado_termina_en_el_mejor_par_del_nucleo() {
        let id = "puerto.disenar_cebadores";
        let casos = [
            entrada_de(id),
            con_campos(id, &[("objetivo", "")]),
            con_campos(id, &[("plantilla", "@gen:araC"), ("objetivo", "")]),
            con_campos(
                id,
                &[
                    (
                        "plantilla",
                        "ATGCGTACGTTAGCCTAGGCTAATCGGATCCGATTACAGGCTTAAGCCGGATAACGTTCAGGCATCGATCGGCTAAGCTTAGGCCATGCAAGTCCGATTGCAGCTAGGCTTAACGGATCCGTAAGCTTGGCATCCAGTAACGGTCAAGCTTGCCATGGATCCAAGTCGGCAATTCGGCGTACGATCCGATGCAAGCTTAGC",
                    ),
                    ("objetivo", ""),
                    ("producto", "80-150"),
                ],
            ),
            con_campos(
                id,
                &[
                    (
                        "plantilla",
                        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                    ),
                    ("objetivo", ""),
                    ("producto", "40-90"),
                ],
            ),
            con_campos(
                id,
                &[
                    (
                        "plantilla",
                        "ACGTNNACGTNNACGTNNACGTNNACGTNNACGTNNACGTNNACGTNNACGTNNACGTNNACGTNNACGTNNACGTNNACGTNNACGTNNACGTNN",
                    ),
                    ("objetivo", ""),
                    ("producto", "40-90"),
                ],
            ),
        ];
        for caso in casos {
            let lista = comprobar_lineas(id, caso.clone());
            let entrada = Entrada::desde_json(caso).unwrap();
            let (_, _, diseno) = diseno_de(&entrada);
            let ultimo = lista.last().unwrap();
            assert_eq!(ultimo.linea, "Ok(Diseno {");
            match diseno.pares.first() {
                Some(par) => {
                    assert!(
                        ultimo.nota.contains(&texto(&par.izquierdo.secuencia))
                            && ultimo.nota.contains(&texto(&par.derecho.secuencia))
                    );
                    let primera_fila = &ultimo.estado["tabla"]["filas"][0];
                    assert_eq!(primera_fila[1], texto(&par.izquierdo.secuencia));
                    assert_eq!(primera_fila[2], texto(&par.derecho.secuencia));
                    assert_eq!(
                        ultimo.estado["tabla"]["filas"].as_array().unwrap().len(),
                        diseno.pares.len()
                    );
                }
                None => assert!(ultimo.nota.contains("no hay ningún par")),
            }
        }
    }
}
