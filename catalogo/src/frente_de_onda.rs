use crate::formato::{decimal, miles, porcentaje};
use crate::pizarra::{Marca, Narracion, Pizarra};
use crate::{Algoritmo, Bloque, Campo, Entrada, FilaAlineada, Paso, Salida, Tipo};
use nucleo::alineamiento;
use nucleo::frente_de_onda::{
    self, AlineamientoWfa, BITS_POR_PALABRA, CASTIGOS_DEL_ARTICULO, Castigos, DistanciaMyers,
    Frente, FrentesDeOnda, Origen, Palabra, TablaDeIguales,
};
use std::cell::Cell;
use std::time::{Duration, Instant};

const CAPITULO: &str = "alineamiento";
const LIMITE_CLASICO_AFIN: usize = 25_000_000;
const LIMITE_CLASICO_EDICION: usize = 100_000_000;
const LIMITE_MYERS: usize = 400_000_000;
const MAXIMO_LETRAS_VISIBLES: usize = 20_000;
const MAXIMO_CIGAR: usize = 4000;
const MAXIMO_PASOS: usize = 400;
const DIAGONALES_A_CADA_LADO: i32 = 7;
const FRENTES_EN_DETALLE: usize = 6;
const FRENTES_RESUMIDOS: usize = 60;
const COLUMNAS_EN_DETALLE: usize = 8;
const COLUMNAS_RESUMIDAS: usize = 60;
const FILAS_VISIBLES: usize = 16;
const LETRAS_EN_PIZARRA: usize = 60;

const RECA_MUTADO: &str = ">recA con 3 % de cambios\nATGGCTATCGACGAAAACAAACAGAAAGCGTTGGCGGCAGCACTGGGCCAGATTGAGAAACAATTTGGTAAAGGCTCCATCAAGCGCCTGGGTGAAGACCGTTCCAATGTGGAAACCATCTCTACCGGTTTGCTTTCACTGGATCTCGCGCTTGGGGCAGGTGGTCTGCCGATGGGCCGTATCGTCGAAATCTACGGACCGGAACCTTCCGGTAAAACCACGCTGACGCTGCAGGTGATCGCCGCAGCGCAACGTGAAGGTAAAACCTGTGCGTTTATCGATGCTGAACACGCGCTGGACCCAATCTACGCACGTAAACTGGGCGTCGATATCGACAACCTGCTGTGCTCCCAGCCGGACACCGGCGAGCAGGCACTGGAAATCTGTGACGCCCAGGCGCGTTCTGGCGCAGTAGGCGTTATCGTCGTTGACTCCGTGGCGGCACTGACGCCGAAAGCGGAAATCGAAGGCGAAATCGGCGACTCTCACATGGGCCTTGCGGCACGTATGATGAGCCAGGATGCGTAAGCTGGCGGGTAACCTGAAGCAGTCCAACACGCTGCTGATCTTCATCAACCAGAACCGTATGAAAATCGGTGTGATGTTCGGTAACCCGGAAACCACTACCGGCGGTAACGCGCTAAAATTCTACGCCTCTGTTCTCGACATCCGTCGTACCCGCGCGGTGACCGAAGAGGGCGAAAACGTGGTGGGTAGCGAAACCCGCGTGAACGTGGTGAAGAATTCCAAAATCGCTGCGCCGTTTAAACAGGCTGTATTCCAGATCCTCTACGGCGAACGTATCAACTTCTACGGCGAACTGGTTGACCTGGGCGTAAAAGAGAAGCTGATCGAGAAAGCAGGCGCGTGGTACAGCTACAAAAGTGTGAAGATCGGTCAGGGTACAGCGAATGCGACTGCCTGGCTGAAAGATAACCCGGAAACCGCGACAGAGATCGAGAAGAAAGTACGTGAACGAGTTGCTGCTGAGCAACCCGAACTCAACGCCGGATTGCTCTGTAGATGATAGCGAAGGCGTAGCAGAAACTAGCGAAGATTTTTAA";

thread_local! {
    static CASTIGO_DE_DIFERENCIA: Cell<i32> = const { Cell::new(4) };
}

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![
        Algoritmo {
            id: "alineamiento.wfa",
            capitulo: CAPITULO,
            nombre: "Alineamiento por frentes de onda (WFA)",
            resumen: "Alineamiento global exacto con huecos afines que avanza de puntaje en puntaje por las diagonales en vez de llenar toda la tabla: con secuencias parecidas calcula una fracción mínima de las celdas del método clásico.",
            campos: &[
                Campo {
                    id: "v",
                    etiqueta: "Secuencia v (ADN o proteína)",
                    tipo: Tipo::Proteina,
                    ejemplo: "@gen:recA",
                },
                Campo {
                    id: "w",
                    etiqueta: "Secuencia w (ADN o proteína)",
                    tipo: Tipo::Proteina,
                    ejemplo: RECA_MUTADO,
                },
                Campo {
                    id: "diferencia",
                    etiqueta: "Castigo por letra distinta (x)",
                    tipo: Tipo::Entero,
                    ejemplo: "4",
                },
                Campo {
                    id: "apertura",
                    etiqueta: "Castigo por abrir un hueco (o)",
                    tipo: Tipo::Entero,
                    ejemplo: "6",
                },
                Campo {
                    id: "extension",
                    etiqueta: "Castigo por cada letra del hueco (e)",
                    tipo: Tipo::Entero,
                    ejemplo: "2",
                },
            ],
            funcion: "frente_de_onda::alinear_wfa",
            ejecutar: wfa,
            pasos: Some(pasos_wfa),
        },
        Algoritmo {
            id: "alineamiento.myers",
            capitulo: CAPITULO,
            nombre: "Distancia de edición con vectores de bits (Myers)",
            resumen: "La misma distancia de edición que la tabla clásica, pero cada columna se calcula de 64 filas a la vez con operaciones de bits y solo dentro de una banda de diagonales: miles de veces menos trabajo con secuencias parecidas.",
            campos: &[
                Campo {
                    id: "v",
                    etiqueta: "Secuencia v (ADN o proteína)",
                    tipo: Tipo::Proteina,
                    ejemplo: "@gen:cadA",
                },
                Campo {
                    id: "w",
                    etiqueta: "Secuencia w (ADN o proteína)",
                    tipo: Tipo::Proteina,
                    ejemplo: "@gen:ldcC",
                },
            ],
            funcion: "frente_de_onda::distancia_myers",
            ejecutar: myers,
            pasos: Some(pasos_myers),
        },
    ]
}

struct Pareja {
    v: Vec<u8>,
    w: Vec<u8>,
    nombre_v: String,
    nombre_w: String,
}

fn leer_pareja(entrada: &Entrada) -> Result<Pareja, String> {
    Ok(Pareja {
        v: entrada.proteina("v")?,
        w: entrada.proteina("w")?,
        nombre_v: nombre_de(entrada, "v"),
        nombre_w: nombre_de(entrada, "w"),
    })
}

fn nombre_de(entrada: &Entrada, campo: &str) -> String {
    let texto = entrada.texto(campo).unwrap_or_default();
    let limpio = texto.trim();
    if let Some(nombre) = limpio
        .strip_prefix("@gen:")
        .or_else(|| limpio.strip_prefix("@proteina:"))
    {
        return nombre.trim().to_string();
    }
    if let Some(rango) = limpio.strip_prefix("@ecoli:") {
        return format!("E. coli {}", rango.trim().replace('-', "–"));
    }
    if limpio == "@ecoli" {
        return "E. coli".to_string();
    }
    if let Some(cabecera) = limpio.strip_prefix('>') {
        let nombre: String = cabecera
            .lines()
            .next()
            .unwrap_or_default()
            .trim()
            .chars()
            .take(40)
            .collect();
        if !nombre.is_empty() {
            return nombre;
        }
    }
    campo.to_string()
}

fn leer_castigos(entrada: &Entrada) -> Result<Castigos, String> {
    Ok(Castigos {
        diferencia: leer_castigo(entrada, "diferencia")?,
        apertura: leer_castigo(entrada, "apertura")?,
        extension: leer_castigo(entrada, "extension")?,
    })
}

fn leer_castigo(entrada: &Entrada, campo: &str) -> Result<usize, String> {
    let valor = entrada.entero(campo)?;
    if !(0..=1000).contains(&valor) {
        return Err(format!(
            "El castigo «{campo}» debe estar entre 0 y 1000: es lo que suma a la penalización."
        ));
    }
    Ok(valor as usize)
}

fn leer_wfa(entrada: &Entrada) -> Result<(Pareja, Castigos), String> {
    let pareja = leer_pareja(entrada)?;
    let castigos = leer_castigos(entrada)?;
    frente_de_onda::revisar_entrada(&pareja.v, &pareja.w, castigos)?;
    Ok((pareja, castigos))
}

fn leer_myers(entrada: &Entrada) -> Result<Pareja, String> {
    let pareja = leer_pareja(entrada)?;
    let corta = pareja.v.len().min(pareja.w.len());
    let larga = pareja.v.len().max(pareja.w.len());
    let trabajo = corta.div_ceil(BITS_POR_PALABRA).saturating_mul(larga);
    if trabajo > LIMITE_MYERS {
        return Err(format!(
            "En el peor caso harían falta {} operaciones de palabra ({} × {} letras) y el límite es {}. Recorta las secuencias: hasta unos 100 000 × 250 000 letras.",
            miles(trabajo),
            miles(corta),
            miles(larga),
            miles(LIMITE_MYERS)
        ));
    }
    Ok(pareja)
}

fn parece_adn(secuencia: &[u8]) -> bool {
    !secuencia.is_empty()
        && secuencia
            .iter()
            .all(|letra| matches!(letra, b'A' | b'C' | b'G' | b'T' | b'U' | b'N'))
}

fn fraccion(parte: usize, total: usize) -> f64 {
    if total == 0 {
        0.0
    } else {
        parte as f64 / total as f64
    }
}

fn tiempo(duracion: Duration) -> String {
    let segundos = duracion.as_secs_f64();
    if segundos >= 1.0 {
        format!("{} s", decimal(segundos, 2))
    } else if segundos >= 0.01 {
        format!("{} ms", decimal(segundos * 1000.0, 1))
    } else {
        format!("{} ms", decimal(segundos * 1000.0, 3))
    }
}

fn veces(lento: Duration, rapido: Duration) -> f64 {
    lento.as_secs_f64() / rapido.as_secs_f64().max(1e-9)
}

fn texto_de_veces(razon: f64) -> String {
    if razon >= 100.0 {
        miles(razon.round() as usize)
    } else {
        decimal(razon, 1)
    }
}

fn recortado(texto: &[u8], maximo: usize) -> String {
    if texto.len() <= maximo {
        return String::from_utf8_lossy(texto).into_owned();
    }
    format!("{}…", String::from_utf8_lossy(&texto[..maximo]))
}

fn fila_visible(nombre: &str, fila: &[u8]) -> FilaAlineada {
    FilaAlineada {
        nombre: nombre.to_string(),
        secuencia: String::from_utf8_lossy(&fila[..fila.len().min(MAXIMO_LETRAS_VISIBLES)])
            .into_owned(),
    }
}

struct Cuenta {
    columnas: usize,
    iguales: usize,
    distintas: usize,
    letras_en_huecos: usize,
    huecos: usize,
}

fn contar(operaciones: &[u8]) -> Cuenta {
    let mut cuenta = Cuenta {
        columnas: operaciones.len(),
        iguales: 0,
        distintas: 0,
        letras_en_huecos: 0,
        huecos: 0,
    };
    let mut anterior = b'=';
    for &operacion in operaciones {
        match operacion {
            b'=' => cuenta.iguales += 1,
            b'X' => cuenta.distintas += 1,
            _ => {
                cuenta.letras_en_huecos += 1;
                if operacion != anterior {
                    cuenta.huecos += 1;
                }
            }
        }
        anterior = operacion;
    }
    cuenta
}

fn diferencia_con_castigo(a: u8, b: u8) -> i32 {
    if a == b {
        0
    } else {
        -CASTIGO_DE_DIFERENCIA.with(Cell::get)
    }
}

fn comparar_con_gotoh(v: &[u8], w: &[u8], castigos: Castigos) -> Option<(usize, Duration)> {
    if v.len().saturating_mul(w.len()) > LIMITE_CLASICO_AFIN {
        return None;
    }
    CASTIGO_DE_DIFERENCIA.with(|castigo| castigo.set(castigos.diferencia as i32));
    let inicio = Instant::now();
    let clasico = alineamiento::alineamiento_afin(
        v,
        w,
        diferencia_con_castigo,
        (castigos.apertura + castigos.extension) as i64,
        castigos.extension as i64,
    );
    let duracion = inicio.elapsed();
    Some((clasico.puntaje.unsigned_abs() as usize, duracion))
}

fn wfa(entrada: &Entrada) -> Result<Salida, String> {
    let (pareja, castigos) = leer_wfa(entrada)?;
    let Pareja {
        v,
        w,
        nombre_v,
        nombre_w,
    } = &pareja;
    let inicio = Instant::now();
    let resultado = frente_de_onda::alinear_wfa(v, w, castigos)?;
    let tiempo_wfa = inicio.elapsed();
    let clasico = comparar_con_gotoh(v, w, castigos);
    let celdas_clasicas = 3usize
        .saturating_mul(v.len() + 1)
        .saturating_mul(w.len() + 1);
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        "Penalización óptima",
        miles(resultado.penalizacion),
    ));
    let cuenta = resultado.operaciones.as_deref().map(contar);
    if let (Some(operaciones), Some(cuenta)) = (&resultado.operaciones, &cuenta) {
        salida.agregar(Bloque::valor(
            "Columnas del alineamiento",
            format!(
                "{}: {} iguales, {} distintas y {} en {}",
                miles(cuenta.columnas),
                miles(cuenta.iguales),
                miles(cuenta.distintas),
                cuantas(cuenta.letras_en_huecos, "letra", "letras"),
                cuantas(cuenta.huecos, "hueco", "huecos")
            ),
        ));
        salida.agregar(Bloque::valor(
            "Identidad",
            porcentaje(fraccion(cuenta.iguales, cuenta.columnas), 2),
        ));
        let (fila_v, fila_w) = frente_de_onda::filas_alineadas(v, w, operaciones);
        salida.agregar(Bloque::alineamiento(
            "Alineamiento",
            vec![
                fila_visible(nombre_v, &fila_v),
                fila_visible(nombre_w, &fila_w),
            ],
        ));
        salida.agregar(Bloque::texto(
            "CIGAR (= iguales, X distintas, I letras solo en w, D letras solo en v)",
            recortado(frente_de_onda::cigar(operaciones).as_bytes(), MAXIMO_CIGAR),
        ));
    }
    salida.agregar(Bloque::valor(
        "Frentes de onda calculados",
        format!(
            "{} (puntajes del 0 al {}; los demás no tienen ningún alineamiento con ese costo exacto)",
            miles(resultado.frentes),
            miles(resultado.penalizacion)
        ),
    ));
    salida.agregar(Bloque::valor(
        "Celdas calculadas por WFA (M̃, Ĩ y D̃)",
        miles(resultado.celdas),
    ));
    salida.agregar(Bloque::valor(
        "Celdas de la programación dinámica (3 tablas de (n+1)·(m+1))",
        miles(celdas_clasicas),
    ));
    salida.agregar(Bloque::valor(
        "Parte de la tabla que calculó WFA",
        porcentaje(fraccion(resultado.celdas, celdas_clasicas), 3),
    ));
    salida.agregar(Bloque::barras(
        "Celdas calculadas con la misma entrada",
        "Método",
        "Celdas",
        vec![
            "WFA".to_string(),
            "Programación dinámica (Gotoh)".to_string(),
        ],
        vec![resultado.celdas as f64, celdas_clasicas as f64],
    ));
    salida.agregar(Bloque::valor("Tiempo de WFA", tiempo(tiempo_wfa)));
    if let Some((_, tiempo_clasico)) = clasico {
        salida.agregar(Bloque::valor(
            "Tiempo de «Alineamiento con huecos afines» (Gotoh)",
            tiempo(tiempo_clasico),
        ));
        salida.agregar(Bloque::barras(
            "Tiempo con la misma entrada",
            "Método",
            "Milisegundos",
            vec![
                "WFA".to_string(),
                "Alineamiento con huecos afines".to_string(),
            ],
            vec![
                tiempo_wfa.as_secs_f64() * 1000.0,
                tiempo_clasico.as_secs_f64() * 1000.0,
            ],
        ));
    }
    explicar_wfa(
        &mut salida,
        &pareja,
        castigos,
        &resultado,
        cuenta.as_ref(),
        (tiempo_wfa, clasico, celdas_clasicas),
    );
    Ok(salida)
}

fn explicar_wfa(
    salida: &mut Salida,
    pareja: &Pareja,
    castigos: Castigos,
    resultado: &AlineamientoWfa,
    cuenta: Option<&Cuenta>,
    medidas: (Duration, Option<(usize, Duration)>, usize),
) {
    let (tiempo_wfa, clasico, celdas_clasicas) = medidas;
    let Castigos {
        diferencia,
        apertura,
        extension,
    } = castigos;
    match cuenta {
        Some(cuenta) => salida.explicar(format!(
            "La penalización óptima es {}: {} por {diferencia} y {} con {} en total, donde un hueco de L letras cuesta {apertura} + L × {extension}. Ningún alineamiento global de estas dos secuencias cuesta menos.",
            miles(resultado.penalizacion),
            cuantas(cuenta.distintas, "letra distinta", "letras distintas"),
            cuantas(cuenta.huecos, "hueco", "huecos"),
            cuantas(cuenta.letras_en_huecos, "letra", "letras")
        )),
        None => salida.explicar(format!(
            "La penalización óptima es {}, con un castigo de {diferencia} por letra distinta y {apertura} + L × {extension} por cada hueco de L letras. Guardar todos los frentes superaba la memoria permitida, así que WFA siguió solo con los últimos: el puntaje es exacto pero no hay alineamiento. La versión bidireccional (BiWFA, 2023) lo recupera con memoria proporcional a la penalización.",
            miles(resultado.penalizacion)
        )),
    }
    let parte = fraccion(resultado.celdas, celdas_clasicas);
    if parte < 0.01 {
        salida.explicar(format!(
            "WFA solo calculó el {} de las celdas que llena la programación dinámica: su trabajo crece con el largo por la penalización (n·s) y no con n·m, así que brilla cuando las secuencias se parecen.",
            porcentaje(parte, 3)
        ));
    } else if parte < 0.5 {
        salida.explicar(format!(
            "WFA calculó el {} de las celdas de la tabla clásica: todavía se ahorra trabajo, pero cuanto más distintas son las secuencias, mayor es la penalización y más diagonales recorre cada frente.",
            porcentaje(parte, 1)
        ));
    } else {
        salida.explicar(
            "WFA calculó casi tantas celdas como la tabla clásica, o más: con secuencias tan distintas su trabajo crece con el cuadrado de la penalización y pierde la ventaja. Para estas usa «Distancia de edición con vectores de bits (Myers)» o el método clásico.",
        );
    }
    match clasico {
        Some((penalizacion_clasica, tiempo_clasico))
            if penalizacion_clasica == resultado.penalizacion =>
        {
            let razon = veces(tiempo_clasico, tiempo_wfa);
            let comparacion = if razon >= 1.0 {
                format!("WFA fue {} veces más rápido", texto_de_veces(razon))
            } else {
                "aquí el método clásico fue más rápido".to_string()
            };
            salida.explicar(format!(
                "En esta computadora WFA tardó {} y «Alineamiento con huecos afines» del capítulo 5 (Gotoh, con apertura {} y extensión {extension}) tardó {}: {comparacion}, con la misma penalización óptima.",
                tiempo(tiempo_wfa),
                apertura + extension,
                tiempo(tiempo_clasico)
            ));
        }
        Some((penalizacion_clasica, _)) => salida.explicar(format!(
            "Atención: el método clásico dio una penalización de {} y WFA de {}. Deberían ser iguales; avisa de este caso.",
            miles(penalizacion_clasica),
            miles(resultado.penalizacion)
        )),
        None => salida.explicar(format!(
            "La tabla clásica tendría {} celdas, demasiadas para «Alineamiento con huecos afines» en esta app; WFA lo resolvió en {}.",
            miles(celdas_clasicas),
            tiempo(tiempo_wfa)
        )),
    }
    if let Some(cuenta) = cuenta
        && cuenta.columnas > 0
        && fraccion(cuenta.iguales, cuenta.columnas) >= 0.95
    {
        salida.explicar(format!(
            "Con {} de identidad las secuencias son casi idénticas: el caso ideal de WFA, como una lectura larga contra su referencia o dos copias recientes de un gen. Si hay empates, WFA devuelve uno de los alineamientos óptimos; otro con la misma penalización sería igual de válido.",
            porcentaje(fraccion(cuenta.iguales, cuenta.columnas), 1)
        ));
    }
    if !parece_adn(&pareja.v) || !parece_adn(&pareja.w) {
        salida.explicar(
            "Con proteínas WFA usa un solo castigo para cualquier cambio de letra: no premia los cambios entre aminoácidos parecidos como BLOSUM62 en «Alineamiento con huecos afines».",
        );
    } else if castigos != CASTIGOS_DEL_ARTICULO {
        salida.explicar(
            "Los castigos del artículo original son x = 4, o = 6 y e = 2; con otros valores cambia qué alineamiento es el óptimo y cuántos frentes hacen falta.",
        );
    }
}

fn myers(entrada: &Entrada) -> Result<Salida, String> {
    let pareja = leer_myers(entrada)?;
    let Pareja {
        v,
        w,
        nombre_v,
        nombre_w,
    } = &pareja;
    let inicio = Instant::now();
    let resultado = frente_de_onda::distancia_myers(v, w);
    let tiempo_myers = inicio.elapsed();
    let celdas_clasicas = v.len().saturating_mul(w.len());
    let clasico = if celdas_clasicas <= LIMITE_CLASICO_EDICION {
        let inicio = Instant::now();
        let distancia = alineamiento::distancia_de_edicion(v, w);
        Some((distancia, inicio.elapsed()))
    } else {
        None
    };
    let larga = v.len().max(w.len());
    let (nombre_patron, largo_patron) = if resultado.patron_es_v {
        (nombre_v, v.len())
    } else {
        (nombre_w, w.len())
    };
    let palabras_calculadas: usize = resultado
        .pasadas
        .iter()
        .map(|pasada| pasada.palabras_calculadas)
        .sum();
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        "Distancia de edición",
        miles(resultado.distancia),
    ));
    salida.agregar(Bloque::valor(
        "Cambios por cada 100 letras",
        decimal(100.0 * fraccion(resultado.distancia, larga), 2),
    ));
    salida.agregar(Bloque::valor(
        "Patrón empacado en bits",
        format!(
            "{nombre_patron}: {} letras en {} palabras de 64 bits por columna",
            miles(largo_patron),
            miles(resultado.palabras_por_columna)
        ),
    ));
    let filas = resultado
        .pasadas
        .iter()
        .map(|pasada| {
            vec![
                if pasada.limite >= v.len() + w.len() {
                    "toda la tabla".to_string()
                } else {
                    miles(pasada.limite)
                },
                match pasada.distancia {
                    Some(distancia) => format!("{} (exacta)", miles(distancia)),
                    None => format!("más de {}: se duplica el límite", miles(pasada.limite)),
                },
                miles(pasada.palabras_calculadas),
            ]
        })
        .collect();
    salida.agregar(Bloque::tabla(
        "Pasadas con banda",
        &["Límite k", "Distancia", "Palabras calculadas"],
        filas,
    ));
    salida.agregar(Bloque::valor(
        "Operaciones de palabra (64 celdas cada una)",
        miles(palabras_calculadas),
    ));
    salida.agregar(Bloque::valor(
        "Celdas de la tabla clásica (n·m)",
        miles(celdas_clasicas),
    ));
    salida.agregar(Bloque::valor("Tiempo de Myers", tiempo(tiempo_myers)));
    if let Some((_, tiempo_clasico)) = clasico {
        salida.agregar(Bloque::valor(
            "Tiempo de «Distancia de edición» (tabla clásica)",
            tiempo(tiempo_clasico),
        ));
        salida.agregar(Bloque::barras(
            "Tiempo con la misma entrada",
            "Método",
            "Milisegundos",
            vec![
                "Myers".to_string(),
                "Distancia de edición clásica".to_string(),
            ],
            vec![
                tiempo_myers.as_secs_f64() * 1000.0,
                tiempo_clasico.as_secs_f64() * 1000.0,
            ],
        ));
    }
    explicar_myers(
        &mut salida,
        &pareja,
        &resultado,
        palabras_calculadas,
        (tiempo_myers, clasico),
    );
    Ok(salida)
}

fn explicar_myers(
    salida: &mut Salida,
    pareja: &Pareja,
    resultado: &DistanciaMyers,
    palabras_calculadas: usize,
    medidas: (Duration, Option<(usize, Duration)>),
) {
    let (tiempo_myers, clasico) = medidas;
    let (v, w) = (&pareja.v, &pareja.w);
    let larga = v.len().max(w.len());
    let por_cien = 100.0 * fraccion(resultado.distancia, larga);
    salida.explicar(format!(
        "La distancia de edición es {}: el menor número de sustituciones, inserciones y borrados de una letra que convierten v en w, {} por cada 100 letras de la secuencia más larga.",
        miles(resultado.distancia),
        decimal(por_cien, 1)
    ));
    match clasico {
        Some((distancia_clasica, tiempo_clasico)) if distancia_clasica == resultado.distancia => {
            let razon = veces(tiempo_clasico, tiempo_myers);
            let comparacion = if razon >= 1.0 {
                format!("Myers fue {} veces más rápido", texto_de_veces(razon))
            } else {
                "aquí la tabla clásica fue más rápida".to_string()
            };
            salida.explicar(format!(
                "Coincide con «Distancia de edición» del capítulo 5, que llenó {} celdas en {}; Myers tardó {}: {comparacion}.",
                miles(v.len() * w.len()),
                tiempo(tiempo_clasico),
                tiempo(tiempo_myers)
            ));
        }
        Some((distancia_clasica, _)) => salida.explicar(format!(
            "Atención: la tabla clásica dio {} y Myers {}. Deberían ser iguales; avisa de este caso.",
            miles(distancia_clasica),
            miles(resultado.distancia)
        )),
        None => salida.explicar(format!(
            "La tabla clásica tendría {} celdas, demasiadas para compararla aquí; Myers la resolvió en {}.",
            miles(v.len().saturating_mul(w.len())),
            tiempo(tiempo_myers)
        )),
    }
    let completa = resultado
        .pasadas
        .last()
        .is_some_and(|pasada| pasada.limite >= v.len() + w.len());
    if completa && resultado.pasadas.len() == 1 {
        salida.explicar(format!(
            "Con un patrón tan corto, o con largos tan distintos, una banda no ahorraba nada: se calculó la tabla entera de una vez, 64 filas por operación ({} en total).",
            cuantas(palabras_calculadas, "operación de palabra", "operaciones de palabra")
        ));
    } else if completa {
        salida.explicar(format!(
            "Las secuencias son tan distintas que las bandas angostas no alcanzaron y al final se calculó la tabla entera, aun así de 64 filas por operación ({} contando las {} pasadas).",
            cuantas(palabras_calculadas, "operación de palabra", "operaciones de palabra"),
            resultado.pasadas.len()
        ));
    } else if let Some(ultima) = resultado.pasadas.last() {
        let pasadas = resultado.pasadas.len();
        let repeticion = if pasadas > 1 {
            format!(
                " Hicieron falta {pasadas} pasadas, duplicando el límite cada vez que la distancia lo superaba."
            )
        } else {
            String::new()
        };
        salida.explicar(format!(
            "Como la distancia no pasa del límite k = {}, bastó mirar una banda de unas {} diagonales alrededor de la principal (Ukkonen): un camino que se saliera de ella costaría más de k.{repeticion}",
            miles(ultima.limite),
            miles(ultima.limite + 1)
        ));
    }
    if parece_adn(v) && parece_adn(w) {
        let corta = v.len().min(w.len());
        let diferencia_de_largo = larga - corta;
        if fraccion(corta, larga) < 0.8 {
            let por_largo = 100.0 * fraccion(diferencia_de_largo, larga);
            salida.explicar(format!(
                "Las secuencias tienen largos muy distintos: {} letras de diferencia, {} por cada 100 letras de la más larga. Una parte grande de la distancia es solo esa diferencia de largo y no dice cuánto se parecen las letras; para saber si una secuencia está contenida en la otra usa un alineamiento local o semiglobal.",
                miles(diferencia_de_largo),
                decimal(por_largo, 1)
            ));
        } else if por_cien >= 45.0 {
            salida.explicar(
                "Dos secuencias de ADN al azar del mismo largo quedan a unos 52 cambios por cada 100 letras: por encima de 45 el parecido es casi el del azar y conviene buscar un alineamiento local.",
            );
        } else if por_cien <= 5.0 {
            salida.explicar(
                "Menos de 5 cambios por cada 100 letras: copias casi idénticas, como alelos o una duplicación reciente.",
            );
        } else {
            salida.explicar(
                "Entre 5 y 45 cambios por cada 100 letras: parecido claro con mucha divergencia, lo típico de genes parálogos, copias que se separaron hace mucho tiempo.",
            );
        }
    }
    salida.explicar(
        "Myers es exacto, no una aproximación: da el mismo número que la tabla clásica, pero solo la distancia, sin el alineamiento. Para un alineamiento con huecos afines usa «Alineamiento por frentes de onda (WFA)».",
    );
}

fn nulo_o_numero(valor: i32) -> String {
    if valor < 0 {
        "·".to_string()
    } else {
        miles(valor as usize)
    }
}

fn en(frente: Option<&Frente>, diagonal: i32) -> i32 {
    frente.map(|frente| frente.en(diagonal)).unwrap_or(-1)
}

fn mas_adelantada(frente: &Frente) -> (i32, i32) {
    let mut mejor = (frente.bajo, -1);
    for (posicion, &columna) in frente.desplazamientos.iter().enumerate() {
        if columna > mejor.1 {
            mejor = (frente.bajo + posicion as i32, columna);
        }
    }
    mejor
}

struct Instantanea {
    filas: Vec<Vec<String>>,
    resaltar: Vec<(usize, usize)>,
    bajo: i32,
    alto: i32,
}

fn instantanea(ondas: &FrentesDeOnda, antes: Option<&Frente>) -> Instantanea {
    let puntaje = ondas.puntaje;
    let parejas = ondas.parejas[puntaje].as_ref();
    let inserciones = ondas.inserciones[puntaje].as_ref();
    let borrados = ondas.borrados[puntaje].as_ref();
    let presentes = [parejas, inserciones, borrados];
    let bajo = presentes
        .iter()
        .flatten()
        .map(|f| f.bajo)
        .min()
        .unwrap_or(0);
    let alto = presentes
        .iter()
        .flatten()
        .map(|f| f.alto)
        .max()
        .unwrap_or(0);
    let centro = parejas.map(|f| mas_adelantada(f).0).unwrap_or(0);
    let hasta = (centro + DIAGONALES_A_CADA_LADO).min(alto);
    let desde = (hasta - 2 * DIAGONALES_A_CADA_LADO).max(bajo);
    let hasta = (desde + 2 * DIAGONALES_A_CADA_LADO).min(alto);
    let mut filas = Vec::new();
    let mut resaltar = Vec::new();
    for diagonal in desde..=hasta {
        let columna = en(parejas, diagonal);
        let celda = if columna < 0 {
            "·".to_string()
        } else {
            format!(
                "({}, {})",
                miles((columna - diagonal) as usize),
                miles(columna as usize)
            )
        };
        let cambio = match antes {
            Some(anterior) => anterior.en(diagonal) != columna,
            None => columna >= 0,
        };
        if cambio {
            resaltar.push((filas.len(), 1));
        }
        filas.push(vec![
            diagonal.to_string(),
            nulo_o_numero(columna),
            nulo_o_numero(en(inserciones, diagonal)),
            nulo_o_numero(en(borrados, diagonal)),
            celda,
        ]);
    }
    Instantanea {
        filas,
        resaltar,
        bajo,
        alto,
    }
}

const COLUMNAS_DEL_FRENTE: &[&str] = &[
    "diagonal k",
    "M̃ (columna j)",
    "Ĩ (inserción)",
    "D̃ (borrado)",
    "celda (i, j) de M̃",
];

fn pizarra_wfa(ondas: &FrentesDeOnda, antes: Option<&Frente>) -> Pizarra {
    let foto = instantanea(ondas, antes);
    Pizarra::nueva()
        .variable("puntaje", ondas.puntaje)
        .variable("diagonal_final", ondas.diagonal_final())
        .variable("frentes", miles(ondas.frentes))
        .variable("celdas", miles(ondas.celdas))
        .tabla(
            &format!(
                "Frente de onda con puntaje {} (diagonales {} a {})",
                ondas.puntaje, foto.bajo, foto.alto
            ),
            COLUMNAS_DEL_FRENTE,
            foto.filas,
            &foto.resaltar,
        )
}

struct Extension {
    marcas: Vec<(usize, usize, Marca)>,
    avanzaron: usize,
    salto: (i32, usize, usize),
    adelantada: (i32, usize),
}

fn extension(antes: &Frente, despues: &Frente, largo_w: usize) -> Extension {
    let mut avanzaron = 0;
    let mut salto = (despues.bajo, 0, 0);
    for diagonal in despues.bajo..=despues.alto {
        let (desde, hasta) = (antes.en(diagonal), despues.en(diagonal));
        if desde < 0 || hasta <= desde {
            continue;
        }
        avanzaron += 1;
        if (hasta - desde) as usize > salto.2 - salto.1 {
            salto = (diagonal, desde as usize, hasta as usize);
        }
    }
    let (diagonal, columna) = mas_adelantada(despues);
    let adelantada = (diagonal, columna.max(0) as usize);
    let (_, desde, hasta) = if avanzaron > 0 {
        salto
    } else {
        (diagonal, adelantada.1, adelantada.1)
    };
    let mut marcas = Vec::new();
    if hasta > desde {
        marcas.push((desde, hasta, Marca::Hallado));
    }
    if hasta < largo_w {
        marcas.push((hasta, hasta + 1, Marca::Actual));
    }
    Extension {
        marcas,
        avanzaron,
        salto,
        adelantada,
    }
}

fn nota_de_fuentes(ondas: &FrentesDeOnda, castigos: Castigos, vacios: usize) -> String {
    let puntaje = ondas.puntaje;
    let abrir = castigos.apertura + castigos.extension;
    let saltados = match vacios {
        0 => String::new(),
        1 => "El puntaje anterior no tiene frente: ningún alineamiento cuesta exactamente eso. "
            .to_string(),
        _ => format!(
            "Los {vacios} puntajes anteriores no tienen frente: ningún alineamiento cuesta exactamente eso. "
        ),
    };
    let huecos = if puntaje >= abrir {
        format!(
            "Ĩ y D̃ abren un hueco desde M̃ del puntaje {} (s − o − e) o alargan uno desde el puntaje {} (s − e)",
            puntaje - abrir,
            puntaje - castigos.extension
        )
    } else if puntaje >= castigos.extension {
        format!(
            "Ĩ y D̃ todavía no pueden abrir huecos (s − o − e sería negativo), solo alargar los del puntaje {}",
            puntaje - castigos.extension
        )
    } else {
        "Ĩ y D̃ todavía no pueden abrir ni alargar huecos".to_string()
    };
    let cambio = if puntaje >= castigos.diferencia {
        format!(
            "M̃ cambia una letra desde M̃ del puntaje {} (s − x)",
            puntaje - castigos.diferencia
        )
    } else {
        "M̃ todavía no puede venir de una letra distinta".to_string()
    };
    format!(
        "{saltados}Puntaje {puntaje}: {huecos}; {cambio} o cierra un hueco tomando Ĩ o D̃. Cada diagonal guarda solo la columna más lejana a la que se llega costando exactamente {puntaje}."
    )
}

type RelatoWfa = (Vec<Paso>, usize, Option<Vec<u8>>);

fn narrar_wfa(pareja: &Pareja, castigos: Castigos) -> Result<RelatoWfa, String> {
    let (v, w) = (&pareja.v, &pareja.w);
    let previsto = frente_de_onda::alinear_wfa(v, w, castigos)?;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let mut ondas = FrentesDeOnda::nuevos(v, w, castigos);
    narracion.agregar(
        pizarra_wfa(&ondas, None)
            .variable("v", recortado(v, LETRAS_EN_PIZARRA))
            .variable("w", recortado(w, LETRAS_EN_PIZARRA))
            .paso(
                "let mut ondas = FrentesDeOnda::nuevos(v, w, castigos);",
                format!(
                    "Vas a alinear {} ({} letras) con {} ({} letras) castigando {} cada letra distinta y {} + L × {} cada hueco de L letras. En vez de una tabla de {} × {}, WFA guarda por cada puntaje s y cada diagonal k = j − i la columna más lejana a la que se llega costando s. Empieza con s = 0 en la diagonal 0, columna 0.",
                    pareja.nombre_v,
                    miles(v.len()),
                    pareja.nombre_w,
                    miles(w.len()),
                    castigos.diferencia,
                    castigos.apertura,
                    castigos.extension,
                    miles(v.len() + 1),
                    miles(w.len() + 1)
                ),
            ),
    );
    let antes = ondas.parejas[0].clone();
    ondas.extender();
    if let (Some(antes), Some(despues)) = (&antes, &ondas.parejas[0]) {
        let avance = extension(antes, despues, w.len());
        narracion.agregar(
            pizarra_wfa(&ondas, Some(antes))
                .secuencia(w, &avance.marcas)
                .paso(
                    "*columna += letras_iguales(&v[i..], &w[j..]) as i32;",
                    format!(
                        "Extender: por la diagonal 0 las letras de v y w son iguales hasta la columna {}, así que el frente avanza {} letras sin pagar nada (compara de 8 en 8 bytes). La meta es la diagonal final {} en la columna {}.",
                        miles(avance.adelantada.1),
                        miles(avance.adelantada.1),
                        ondas.diagonal_final(),
                        miles(w.len())
                    ),
                ),
        );
    }
    let resto = previsto.frentes.saturating_sub(FRENTES_EN_DETALLE);
    let cada = resto.div_ceil(FRENTES_RESUMIDOS).max(1);
    let mut numero = 1;
    let mut vacios = 0;
    while !ondas.llego_al_final() {
        ondas.siguiente()?;
        let puntaje = ondas.puntaje;
        if ondas.parejas[puntaje].is_none() {
            vacios += 1;
            continue;
        }
        numero += 1;
        let en_detalle = numero <= FRENTES_EN_DETALLE;
        let narrado = en_detalle || numero % cada == 0 || numero == previsto.frentes;
        if !narrado {
            ondas.extender();
            vacios = 0;
            continue;
        }
        let antes = ondas.parejas[puntaje].clone();
        if en_detalle {
            narracion.agregar(pizarra_wfa(&ondas, None).paso(
                "ondas.siguiente()?;",
                nota_de_fuentes(&ondas, castigos, vacios),
            ));
        }
        ondas.extender();
        vacios = 0;
        let (Some(antes), Some(despues)) = (antes.as_ref(), ondas.parejas[puntaje].as_ref()) else {
            continue;
        };
        let avance = extension(antes, despues, w.len());
        let (diagonal, columna) = avance.adelantada;
        let (k_salto, desde, hasta) = avance.salto;
        let nota = if en_detalle && avance.avanzaron == 0 {
            format!(
                "Extender: en ninguna de las {} del frente la letra siguiente de v es igual a la de w, así que nada avanza gratis. La más adelantada sigue siendo la diagonal {diagonal}, en la columna {}.",
                cuantas(despues.ancho(), "diagonal", "diagonales"),
                miles(columna)
            )
        } else if en_detalle && despues.ancho() == 1 {
            format!(
                "Extender: en la única diagonal del frente, la {k_salto}, las letras siguen iguales de la columna {} a la {} ({} letras), así que el frente avanza sin pagar nada hasta donde v y w vuelven a diferir.",
                miles(desde),
                miles(hasta),
                miles(hasta - desde)
            )
        } else if en_detalle {
            format!(
                "Extender: {} de {} diagonales avanzaron por letras iguales. El salto más largo fue en la diagonal {k_salto}: de la columna {} a la {} ({} letras iguales seguidas), hasta donde v y w vuelven a diferir. La más adelantada es la diagonal {diagonal}, en la columna {}.",
                miles(avance.avanzaron),
                miles(despues.ancho()),
                miles(desde),
                miles(hasta),
                miles(hasta - desde),
                miles(columna)
            )
        } else {
            format!(
                "Puntaje {puntaje} (frente número {numero}): el frente ya cubre {} diagonales; la más adelantada es la {diagonal}, en la columna {} de {}. Se saltan frentes parecidos.",
                miles(despues.ancho()),
                miles(columna),
                miles(w.len())
            )
        };
        narracion.agregar(
            pizarra_wfa(&ondas, Some(antes))
                .secuencia(w, &avance.marcas)
                .paso("*columna += letras_iguales(&v[i..], &w[j..]) as i32;", nota),
        );
    }
    let penalizacion = ondas.puntaje;
    let operaciones = ondas.rastrear();
    let diagonal_final = ondas.diagonal_final();
    let llegada = format!(
        "En el puntaje {penalizacion}, M̃ en la diagonal final {diagonal_final} llega a la columna {}: ya se alinearon v y w completas. Como los puntajes se recorren de menor a mayor, {penalizacion} es la penalización óptima.",
        miles(w.len())
    );
    let mut ultimo_paso = pizarra_wfa(&ondas, None).paso(
        "Some(frente) => frente.en(self.diagonal_final()) == self.largo_w,",
        llegada,
    );
    let Some(lista) = &operaciones else {
        ultimo_paso.nota.push_str(" Las secuencias necesitaban demasiada memoria para guardar todos los frentes: se obtuvo la penalización exacta sin el alineamiento.");
        return Ok((narracion.terminar(ultimo_paso), penalizacion, None));
    };
    narracion.agregar(ultimo_paso);
    let (antes, origen) = ondas.origen_de_pareja(penalizacion, diagonal_final);
    let iguales_al_final = w.len() as i32 - antes;
    let de_donde = match origen {
        Origen::Inicio => "del inicio: no hubo ningún cambio".to_string(),
        Origen::Diferencia => format!(
            "de una letra distinta, desde M̃ del puntaje {}",
            penalizacion - castigos.diferencia
        ),
        Origen::Insercion => "de cerrar un hueco con letras solo de w (Ĩ)".to_string(),
        Origen::Borrado => "de cerrar un hueco con letras solo de v (D̃)".to_string(),
    };
    narracion.agregar(
        pizarra_wfa(&ondas, None)
            .variable("operaciones", miles(lista.len()))
            .paso(
                "let (antes, origen) = self.origen_de_pareja(puntaje, diagonal);",
                format!(
                    "Ahora rastrear vuelve hacia atrás desde la esquina final. Las últimas {} letras fueron iguales (las agregó extender); antes de eso, el valor de M̃ venía {de_donde}. Así, frente por frente, recupera cada operación sin haber guardado flechas.",
                    miles(iguales_al_final.max(0) as usize)
                ),
            ),
    );
    let forma = frente_de_onda::cigar(lista);
    let (fila_v, fila_w) = frente_de_onda::filas_alineadas(v, w, lista);
    let cuenta = contar(lista);
    let ultimo = pizarra_wfa(&ondas, None)
        .variable("operaciones", miles(lista.len()))
        .lista(
            "Alineamiento (primeras columnas)",
            vec![
                format!("v: {}", recortado(&fila_v, LETRAS_EN_PIZARRA)),
                format!("w: {}", recortado(&fila_w, LETRAS_EN_PIZARRA)),
                format!("CIGAR: {}", recortado(forma.as_bytes(), LETRAS_EN_PIZARRA)),
            ],
            &[2],
        )
        .paso(
            "operaciones.reverse();",
            format!(
                "Las operaciones se juntaron de atrás hacia adelante; al darlas vuelta queda el alineamiento: {} columnas, {} iguales, {} distintas y {} huecos, con penalización {penalizacion}. Mira «Qué significa».",
                miles(cuenta.columnas),
                miles(cuenta.iguales),
                miles(cuenta.distintas),
                miles(cuenta.huecos)
            ),
        );
    Ok((narracion.terminar(ultimo), penalizacion, operaciones))
}

fn pasos_wfa(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (pareja, castigos) = leer_wfa(entrada)?;
    let (pasos, _, _) = narrar_wfa(&pareja, castigos)?;
    Ok(pasos)
}

fn bits(valor: u64, filas: usize) -> String {
    (0..filas.min(BITS_POR_PALABRA))
        .map(|fila| if (valor >> fila) & 1 == 1 { '1' } else { '0' })
        .collect()
}

fn cuantas(cantidad: usize, singular: &str, plural: &str) -> String {
    if cantidad == 1 {
        format!("1 {singular}")
    } else {
        format!("{} {plural}", miles(cantidad))
    }
}

fn cambio_de_la_fila(horizontal: isize) -> &'static str {
    match horizontal {
        1 => "subió 1",
        -1 => "bajó 1",
        _ => "no cambió",
    }
}

#[derive(Debug, Clone, Copy, Default)]
struct Tablero {
    limite: Option<usize>,
    primera: Option<usize>,
    ultima: Option<usize>,
    fondo: Option<usize>,
    letra: Option<u8>,
    iguales: Option<u64>,
    vertical_mas: Option<u64>,
    vertical_menos: Option<u64>,
    x_vertical: Option<u64>,
    x_horizontal: Option<u64>,
    horizontal_mas: Option<u64>,
    horizontal_menos: Option<u64>,
}

impl Tablero {
    fn pizarra(&self, filas: usize) -> Pizarra {
        let numero = |valor: Option<usize>| valor.map(miles).unwrap_or_default();
        let cadena = |valor: Option<u64>| valor.map(|v| bits(v, filas)).unwrap_or_default();
        Pizarra::nueva()
            .variable("limite", numero(self.limite))
            .variable("primera", numero(self.primera))
            .variable("ultima", numero(self.ultima))
            .variable("fondo", numero(self.fondo))
            .variable(
                "letra",
                self.letra
                    .map(|letra| (letra as char).to_string())
                    .unwrap_or_default(),
            )
            .variable("iguales", cadena(self.iguales))
            .variable("vertical_mas", cadena(self.vertical_mas))
            .variable("vertical_menos", cadena(self.vertical_menos))
            .variable("x_vertical", cadena(self.x_vertical))
            .variable("x_horizontal", cadena(self.x_horizontal))
            .variable("horizontal_mas", cadena(self.horizontal_mas))
            .variable("horizontal_menos", cadena(self.horizontal_menos))
    }

    fn con_palabra(mut self, palabra: Palabra) -> Tablero {
        self.vertical_mas = Some(palabra.vertical_mas);
        self.vertical_menos = Some(palabra.vertical_menos);
        self
    }
}

fn tabla_de_columna(patron: &[u8], palabra: Palabra, j: usize) -> (Vec<Vec<String>>, usize) {
    let filas = patron.len().min(FILAS_VISIBLES);
    let mut valor = j as i64;
    let mut tabla = Vec::with_capacity(filas);
    for (fila, &letra) in patron.iter().take(filas).enumerate() {
        let mas = (palabra.vertical_mas >> fila) & 1;
        let menos = (palabra.vertical_menos >> fila) & 1;
        valor += mas as i64 - menos as i64;
        tabla.push(vec![
            (fila + 1).to_string(),
            (letra as char).to_string(),
            mas.to_string(),
            menos.to_string(),
            valor.to_string(),
        ]);
    }
    (tabla, filas)
}

struct Columna {
    j: usize,
    antes: Palabra,
    despues: Palabra,
    horizontal: isize,
}

fn narrar_columna(
    narracion: &mut Narracion,
    patron: &[u8],
    base: Tablero,
    columna: &Columna,
    filas: usize,
) {
    let Palabra {
        vertical_mas,
        vertical_menos,
    } = columna.antes;
    let iguales = base.iguales.unwrap_or(0);
    let x_vertical = iguales | vertical_menos;
    let x_horizontal =
        ((iguales & vertical_mas).wrapping_add(vertical_mas) ^ vertical_mas) | iguales;
    let horizontal_mas = vertical_menos | !(x_horizontal | vertical_mas);
    let horizontal_menos = vertical_mas & x_horizontal;
    let visibles = |valor: u64| (valor & mascara_de_filas(filas)).count_ones();
    let j = columna.j;
    let letra = base.letra.map(|letra| letra as char).unwrap_or(' ');
    let primero = base.con_palabra(columna.antes);
    let segundo = Tablero {
        x_vertical: Some(x_vertical),
        x_horizontal: Some(x_horizontal),
        ..primero
    };
    let tercero = Tablero {
        horizontal_mas: Some(horizontal_mas),
        horizontal_menos: Some(horizontal_menos),
        ..segundo
    };
    let fondo = base
        .fondo
        .unwrap_or(0)
        .saturating_add_signed(columna.horizontal);
    let cuarto = Tablero {
        fondo: Some(fondo),
        ..tercero
    }
    .con_palabra(columna.despues);
    if j <= 2 {
        narracion.agregar(primero.pizarra(filas).paso(
            "let iguales = &tabla.de(letra)[primera..=ultima];",
            format!(
                "Columna {j}: la letra del texto es «{letra}». iguales marca con 1 las filas del patrón que tienen esa letra ({}). Es la única vez que se mira la letra en toda la columna.",
                cuantas(visibles(iguales) as usize, "fila", "filas")
            ),
        ));
        narracion.agregar(segundo.pizarra(filas).paso(
            "((iguales & vertical_mas).wrapping_add(vertical_mas) ^ vertical_mas) | iguales;",
            format!(
                "x_horizontal marca {} donde llegar en diagonal desde la columna anterior no cuesta nada: por letra igual, o porque la suma arrastra el acarreo por una corrida de +1 de vertical_mas, como la tabla clásica arrastra un empate hacia abajo. Es un solo paso para 64 filas.",
                cuantas(visibles(x_horizontal) as usize, "fila", "filas")
            ),
        ));
        narracion.agregar(tercero.pizarra(filas).paso(
            "let horizontal_mas = vertical_menos | !(x_horizontal | vertical_mas);",
            format!(
                "horizontal_mas y horizontal_menos dicen si cada fila subió o bajó respecto de la columna anterior: subieron {} y bajaron {}. Solo OR, AND y NOT, sin recorrer filas.",
                cuantas(visibles(horizontal_mas) as usize, "fila", "filas"),
                cuantas(visibles(horizontal_menos) as usize, "fila", "filas")
            ),
        ));
    }
    let (tabla, mostradas) = tabla_de_columna(patron, columna.despues, j);
    let nota = if j <= 2 {
        format!(
            "Corriendo horizontal_mas y horizontal_menos una fila salen los nuevos vertical_mas y vertical_menos: la columna {j} entera de la tabla clásica, guardada como diferencias +1 y −1 entre filas vecinas (la tabla de la derecha las suma). La fila de abajo {}: ahora vale {}.",
            cambio_de_la_fila(columna.horizontal),
            miles(fondo)
        )
    } else {
        format!(
            "Columna {j} (letra «{letra}»): las mismas operaciones de bits dan la columna nueva. La fila de abajo {} y vale {}.",
            cambio_de_la_fila(columna.horizontal),
            miles(fondo)
        )
    };
    narracion.agregar(
        cuarto
            .pizarra(filas)
            .tabla(
                &format!("Columna {j} de la tabla clásica (primeras {mostradas} filas)"),
                &[
                    "fila",
                    "letra del patrón",
                    "vertical_mas",
                    "vertical_menos",
                    "valor",
                ],
                tabla,
                &[(mostradas - 1, 4)],
            )
            .paso(
                "palabra.vertical_mas = horizontal_menos | !(x_vertical | horizontal_mas);",
                nota,
            ),
    );
}

fn mascara_de_filas(filas: usize) -> u64 {
    if filas >= BITS_POR_PALABRA {
        u64::MAX
    } else {
        (1 << filas) - 1
    }
}

fn narrar_primera_pasada(
    narracion: &mut Narracion,
    patron: &[u8],
    texto: &[u8],
    tabla: &TablaDeIguales,
    limite: usize,
) -> (Option<usize>, Tablero) {
    let resta = texto.len() - patron.len();
    let filas_arriba = (limite + resta) / 2;
    let filas_abajo = (limite - resta) / 2;
    let ultima_palabra = tabla.palabras - 1;
    let bit_de_la_ultima_fila = 1u64 << ((patron.len() - 1) % BITS_POR_PALABRA);
    let filas = patron.len().min(BITS_POR_PALABRA);
    let mut palabras = vec![Palabra::inicial(); tabla.palabras];
    let (mut primera, mut ultima) = (0, 0);
    let mut fondo = filas;
    let mut tablero = Tablero {
        limite: Some(limite),
        ..Tablero::default()
    };
    let cada = texto
        .len()
        .saturating_sub(COLUMNAS_EN_DETALLE)
        .div_ceil(COLUMNAS_RESUMIDAS)
        .max(1);
    for (indice, &letra) in texto.iter().enumerate() {
        let j = indice + 1;
        let fila_mas_baja = (j + filas_abajo).min(patron.len());
        while ultima < (fila_mas_baja - 1) / BITS_POR_PALABRA {
            ultima += 1;
            fondo += (patron.len() - ultima * BITS_POR_PALABRA).min(BITS_POR_PALABRA);
        }
        let fila_mas_alta = j.saturating_sub(filas_arriba).max(1);
        primera = primera.max((fila_mas_alta - 1) / BITS_POR_PALABRA);
        let mascaras = tabla.de(letra);
        let antes = palabras[0];
        let base = Tablero {
            limite: Some(limite),
            primera: Some(primera),
            ultima: Some(ultima),
            fondo: Some(fondo),
            letra: Some(letra),
            iguales: (primera == 0).then_some(mascaras[0]),
            ..Tablero::default()
        };
        let mut horizontal = 1;
        for numero in primera..=ultima {
            let bit_alto = if numero == ultima_palabra {
                bit_de_la_ultima_fila
            } else {
                1 << 63
            };
            horizontal = frente_de_onda::avanzar_palabra(
                &mut palabras[numero],
                mascaras[numero],
                horizontal,
                bit_alto,
            );
        }
        fondo = fondo.saturating_add_signed(horizontal);
        let palabra_cero = primera == 0;
        tablero = Tablero {
            fondo: Some(fondo),
            ..base
        };
        if palabra_cero {
            tablero = tablero.con_palabra(palabras[0]);
        }
        if j <= COLUMNAS_EN_DETALLE && palabra_cero {
            let columna = Columna {
                j,
                antes,
                despues: palabras[0],
                horizontal,
            };
            narrar_columna(narracion, patron, base, &columna, filas);
        } else if j % cada == 0 || j == texto.len() {
            let nota = if palabra_cero {
                format!(
                    "Columna {} (letra «{}»): la fila de abajo {} y vale {}. Se saltan columnas que hacen lo mismo.",
                    miles(j),
                    letra as char,
                    cambio_de_la_fila(horizontal),
                    miles(fondo)
                )
            } else {
                format!(
                    "Columna {} (letra «{}»): la banda ya dejó atrás la palabra 0; solo se calculan las palabras {} a {}. La fila más baja de la banda vale {}.",
                    miles(j),
                    letra as char,
                    primera,
                    ultima,
                    miles(fondo)
                )
            };
            narracion.agregar(
                tablero
                    .pizarra(filas)
                    .paso("fondo = fondo.saturating_add_signed(horizontal);", nota),
            );
        }
    }
    ((fondo <= limite).then_some(fondo), tablero)
}

fn narrar_myers(pareja: &Pareja) -> (Vec<Paso>, usize, Option<Option<usize>>) {
    let (v, w) = (&pareja.v, &pareja.w);
    let resultado = frente_de_onda::distancia_myers(v, w);
    let (patron, texto, nombre_patron, nombre_texto) = if resultado.patron_es_v {
        (v, w, &pareja.nombre_v, &pareja.nombre_w)
    } else {
        (w, v, &pareja.nombre_w, &pareja.nombre_v)
    };
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let filas = patron.len().min(BITS_POR_PALABRA);
    let inicio = Tablero::default().con_palabra(Palabra::inicial());
    narracion.agregar(inicio.pizarra(filas).paso(
        "let (patron, texto) = if patron_es_v { (v, w) } else { (w, v) };",
        format!(
            "Vas a medir la distancia de edición entre {} ({} letras) y {} ({} letras). La más corta, {nombre_patron}, es el patrón: sus filas se empacan en palabras de 64 bits. La otra, {nombre_texto}, se recorre columna por columna. La tabla clásica tendría {} celdas.",
            pareja.nombre_v,
            miles(v.len()),
            pareja.nombre_w,
            miles(w.len()),
            miles(v.len().saturating_mul(w.len()))
        ),
    ));
    if patron.is_empty() {
        let ultimo = Pizarra::nueva().variable("distancia", texto.len()).paso(
            "distancia: Some(texto.len()),",
            format!(
                "Una de las secuencias está vacía: hay que insertar las {} letras de la otra. La distancia es {}.",
                miles(texto.len()),
                miles(texto.len())
            ),
        );
        return (narracion.terminar(ultimo), resultado.distancia, None);
    }
    let tabla = TablaDeIguales::para(patron);
    let mut distintas: Vec<u8> = patron.to_vec();
    distintas.sort_unstable();
    distintas.dedup();
    let filas_de_letras: Vec<Vec<String>> = distintas
        .iter()
        .take(24)
        .map(|&letra| vec![(letra as char).to_string(), bits(tabla.de(letra)[0], filas)])
        .collect();
    narracion.agregar(
        inicio
            .pizarra(filas)
            .tabla(
                "iguales de cada letra (primeras 64 filas del patrón)",
                &["letra", "un 1 donde el patrón tiene esa letra"],
                filas_de_letras,
                &[],
            )
            .paso(
                "let tabla = TablaDeIguales::para(patron);",
                format!(
                    "Antes de empezar se arma una máscara por letra, de {} de 64 bits, con un 1 en las filas donde aparece. vertical_mas empieza todo en 1 y vertical_menos en 0: en la columna 0 cada fila vale una más que la de arriba (0, 1, 2, …). El primer carácter de cada cadena de bits es la fila 1.",
                    cuantas(tabla.palabras, "palabra", "palabras")
                ),
            ),
    );
    let mut replica = None;
    for (numero, pasada) in resultado.pasadas.iter().enumerate() {
        let completa = pasada.limite >= v.len() + w.len();
        let (linea, nota) = if numero == 0 && completa {
            (
                "limite = patron.len() + texto.len();",
                "El patrón es corto o las secuencias difieren mucho en largo: una banda no ahorraría casi nada, así que se calcula la tabla entera, palabra por palabra.".to_string(),
            )
        } else if numero == 0 {
            (
                "let mut limite = (texto.len() - patron.len()).max(BITS_POR_PALABRA);",
                format!(
                    "Primera apuesta: la distancia no pasa de k = {}. Si es así, el camino óptimo nunca se aleja más de k diagonales de la principal (Ukkonen), y basta calcular las palabras que tocan esa banda.",
                    miles(pasada.limite)
                ),
            )
        } else if completa {
            (
                "limite = patron.len() + texto.len();",
                "La banda ya cubriría buena parte de la tabla, así que esta vez se calcula entera."
                    .to_string(),
            )
        } else {
            (
                "limite *= 2;",
                format!("Nueva apuesta con el doble: k = {}.", miles(pasada.limite)),
            )
        };
        let apuesta = Tablero {
            limite: Some(pasada.limite),
            ..inicio
        };
        narracion.agregar(apuesta.pizarra(filas).paso(linea, nota));
        let mut al_final = Tablero {
            limite: Some(pasada.limite),
            fondo: pasada.distancia,
            ..Tablero::default()
        };
        if numero == 0 {
            let (distancia, tablero) =
                narrar_primera_pasada(&mut narracion, patron, texto, &tabla, pasada.limite);
            replica = Some(distancia);
            al_final = tablero;
        }
        let cierre = match pasada.distancia {
            Some(distancia) => format!(
                "La pasada terminó: la última fila de la última columna vale {}, que no pasa del límite {}, así que es la distancia exacta. Se calcularon {} de 64 bits.",
                miles(distancia),
                miles(pasada.limite),
                cuantas(pasada.palabras_calculadas, "palabra", "palabras")
            ),
            None => format!(
                "La pasada terminó con un valor mayor que el límite {}: el camino óptimo pudo salirse de la banda, así que ese número no sirve y hay que repetir con un límite mayor ({} calculadas).",
                miles(pasada.limite),
                cuantas(pasada.palabras_calculadas, "palabra", "palabras")
            ),
        };
        narracion.agregar(al_final.pizarra(filas).paso(
            "let distancia = (fondo <= limite).then_some(fondo);",
            cierre,
        ));
    }
    let total: usize = resultado
        .pasadas
        .iter()
        .map(|pasada| pasada.palabras_calculadas)
        .sum();
    let ultimo = Pizarra::nueva()
        .variable("distancia", miles(resultado.distancia))
        .variable("palabras_calculadas", miles(total))
        .paso(
            "if let Some(distancia) = pasada.distancia {",
            format!(
                "La distancia de edición es {}. En total se hicieron {} de hasta 64 celdas cada una, en lugar de {} celdas de la tabla clásica. Mira «Qué significa».",
                miles(resultado.distancia),
                cuantas(total, "operación de palabra", "operaciones de palabra"),
                miles(v.len().saturating_mul(w.len()))
            ),
        );
    (narracion.terminar(ultimo), resultado.distancia, replica)
}

fn pasos_myers(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let pareja = leer_myers(entrada)?;
    let (pasos, _, _) = narrar_myers(&pareja);
    Ok(pasos)
}

#[cfg(test)]
mod pruebas {
    use super::*;

    fn entrada_de_ejemplo(id: &str) -> Entrada {
        let algoritmo = algoritmos()
            .into_iter()
            .find(|algoritmo| algoritmo.id == id)
            .unwrap();
        let pares: Vec<(&str, &str)> = algoritmo
            .campos
            .iter()
            .map(|campo| (campo.id, campo.ejemplo))
            .collect();
        Entrada::nueva(&pares)
    }

    #[test]
    fn el_paso_a_paso_de_wfa_llega_al_mismo_alineamiento() {
        let entrada = entrada_de_ejemplo("alineamiento.wfa");
        let (pareja, castigos) = leer_wfa(&entrada).unwrap();
        let (pasos, penalizacion, operaciones) = narrar_wfa(&pareja, castigos).unwrap();
        let directo = frente_de_onda::alinear_wfa(&pareja.v, &pareja.w, castigos).unwrap();
        assert_eq!(penalizacion, directo.penalizacion);
        assert_eq!(operaciones, directo.operaciones);
        assert!(pasos.len() > 10 && pasos.len() <= MAXIMO_PASOS);
    }

    #[test]
    fn el_paso_a_paso_de_myers_llega_a_la_misma_distancia() {
        let entrada = entrada_de_ejemplo("alineamiento.myers");
        let pareja = leer_myers(&entrada).unwrap();
        let directo = frente_de_onda::distancia_myers(&pareja.v, &pareja.w);
        let (pasos, distancia, replica) = narrar_myers(&pareja);
        assert_eq!(distancia, directo.distancia);
        assert_eq!(replica, Some(directo.pasadas[0].distancia));
        assert_eq!(
            alineamiento::distancia_de_edicion(&pareja.v, &pareja.w),
            directo.distancia
        );
        assert!(pasos.len() > 10 && pasos.len() <= MAXIMO_PASOS);
    }

    #[test]
    fn el_paso_a_paso_de_myers_con_patron_corto_reconstruye_la_columna() {
        let pareja = Pareja {
            v: b"GATTACA".to_vec(),
            w: b"GCATGCTA".to_vec(),
            nombre_v: "v".to_string(),
            nombre_w: "w".to_string(),
        };
        let (pasos, distancia, replica) = narrar_myers(&pareja);
        assert_eq!(
            distancia,
            alineamiento::distancia_de_edicion(&pareja.v, &pareja.w)
        );
        assert_eq!(replica, Some(Some(distancia)));
        let ultima_tabla = pasos
            .iter()
            .rev()
            .find_map(|paso| paso.estado.get("tabla"))
            .unwrap();
        let filas = ultima_tabla["filas"].as_array().unwrap();
        let ultima_fila = filas.last().unwrap().as_array().unwrap();
        assert_eq!(ultima_fila[4].as_str().unwrap(), distancia.to_string());
    }

    #[test]
    fn el_ejemplo_de_wfa_es_un_reca_con_cambios() {
        let entrada = entrada_de_ejemplo("alineamiento.wfa");
        let (pareja, castigos) = leer_wfa(&entrada).unwrap();
        assert_eq!(pareja.nombre_w, "recA con 3 % de cambios");
        let resultado = frente_de_onda::alinear_wfa(&pareja.v, &pareja.w, castigos).unwrap();
        let cuenta = contar(resultado.operaciones.as_deref().unwrap());
        let identidad = fraccion(cuenta.iguales, cuenta.columnas);
        assert!(identidad > 0.95 && identidad < 0.99, "{identidad}");
        let salida = wfa(&entrada).unwrap();
        assert!((2..=6).contains(&salida.interpretacion.len()));
        let salida = myers(&entrada_de_ejemplo("alineamiento.myers")).unwrap();
        assert!((2..=6).contains(&salida.interpretacion.len()));
    }

    #[test]
    fn los_castigos_imposibles_se_rechazan_con_un_mensaje() {
        for (campo, valor) in [("diferencia", "0"), ("extension", "0"), ("apertura", "-1")] {
            let mut pares = vec![
                ("v", "ACGTACGT"),
                ("w", "ACGAACGT"),
                ("diferencia", "4"),
                ("apertura", "6"),
                ("extension", "2"),
            ];
            for par in pares.iter_mut() {
                if par.0 == campo {
                    par.1 = valor;
                }
            }
            let entrada = Entrada::nueva(&pares);
            assert!(wfa(&entrada).is_err());
            assert!(pasos_wfa(&entrada).is_err());
        }
    }

    #[test]
    fn con_largos_muy_distintos_myers_no_habla_de_azar() {
        let corta = "ACGTTGCAAC".repeat(10);
        let larga = "TTGACCAGTA".repeat(100);
        let entrada = Entrada::nueva(&[("v", corta.as_str()), ("w", larga.as_str())]);
        let texto = myers(&entrada).unwrap().interpretacion.join(" ");
        assert!(!texto.contains("al azar"));
        assert!(texto.contains("largos muy distintos"));
        let parecida = "ACGTTGCAAC".repeat(10);
        let otra = "TTGACCAGTA".repeat(10);
        let entrada = Entrada::nueva(&[("v", parecida.as_str()), ("w", otra.as_str())]);
        let texto = myers(&entrada).unwrap().interpretacion.join(" ");
        assert!(!texto.contains("largos muy distintos"));
    }
}
