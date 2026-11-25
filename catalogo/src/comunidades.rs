use crate::formato::{decimal, miles};
use crate::pizarra::{Narracion, Pizarra};
use crate::{Algoritmo, Bloque, Campo, Entrada, Marca, Paso, Punto, Salida, Tipo};
use nucleo::agrupamiento::{self, Inicio, Medida};
use nucleo::azar::Azar;
use nucleo::comunidades::{
    self, CADA_CUANTO_MEDIR, Comunidades, Grafo, INTENTOS_DE_BETA, ITERACIONES_EXAGERADAS, Louvain,
    MAXIMO_PASADAS, MOMENTO_FINAL, MOMENTO_INICIAL, OpcionesTsne, PESO_MINIMO_SNN, Particion,
    TOLERANCIA_DE_ENTROPIA, Tsne,
};
use std::collections::BTreeMap;
use std::time::Instant;

const MAXIMO_FILAS: usize = 5000;
const MAXIMO_FILAS_TSNE: usize = 2000;
const MAXIMO_COLUMNAS: usize = 200;
const MAXIMO_VALOR: f64 = 1e12;
const MAXIMO_VECINOS: usize = 100;
const MAXIMO_RESOLUCION: f64 = 10.0;
const MAXIMO_EXAGERACION: f64 = 50.0;
const MAXIMO_ITERACIONES: usize = 5000;
const MAXIMO_CUENTAS_TSNE: f64 = 4e9;
const MAXIMO_K_MEDIAS: usize = 50;
const ITERACIONES_K_MEDIAS: usize = 100;
const MAXIMO_FILAS_DE_CALOR: usize = 300;
const MAXIMO_BARRAS: usize = 50;
const MAXIMO_CRUCE: usize = 10;
const NOMBRES_POR_COMUNIDAD: usize = 12;
const NOMBRES_A_LA_VISTA: usize = 5;
const MODULARIDAD_CLARA: f64 = 0.3;
const MARGEN_SOBRE_AZAR: f64 = 0.1;
const SEMILLA_REVUELTA: u64 = 0x5EED_0123_4567_89AB;
const ACUERDO_ALTO: f64 = 0.9;
const ACUERDO_MEDIO: f64 = 0.5;
const PERPLEJIDAD_BAJA: f64 = 5.0;
const PERPLEJIDAD_ALTA: f64 = 50.0;
const CAIDA_FINAL_NOTABLE: f64 = 0.01;
const ISLAS_CLARAS: f64 = 0.9;
const VENTAJA_CLARA: f64 = 0.1;
const NODOS_DETALLADOS: usize = 40;
const DISTANCIAS_DETALLADAS: usize = 6;
const FILAS_DE_DISTANCIAS: usize = 12;
const INTENTOS_DETALLADOS: usize = 10;
const ITERACIONES_DETALLADAS: usize = 5;
const CANDIDATAS_A_LA_VISTA: usize = 6;
const VECINOS_A_LA_VISTA: usize = 8;
const MAXIMO_PASOS: usize = 400;

const MEDIDA_PEARSON: &str = "Pearson (1 − r)";
const MEDIDA_EUCLIDEA: &str = "euclídea";
const MEDIDAS: &[&str] = &[MEDIDA_PEARSON, MEDIDA_EUCLIDEA];

const NOMBRE_LLOYD: &str = "k-means (algoritmo de Lloyd)";
const NOMBRE_COMPONENTES: &str = "Componentes principales (vista en 2D)";

const ETIQUETA_TABLA: &str = "Tabla de expresión: un gen por línea, su nombre y sus números separados por espacios, tabulaciones o comas (punto decimal)";

const TABLA_EJEMPLO: &str = "gen	0h	2h	4h	6h	8h	10h	12h	14h
gen01	2.45	1.04	1.57	-0.03	-0.32	-0.44	-1.74	-1.72
gen02	-0.75	-0.81	-0.44	0.02	0.57	0.74	1.56	1.50
gen03	-1.14	-0.75	-0.30	-0.21	1.48	0.92	1.21	1.20
gen04	2.13	1.93	1.20	0.46	-0.67	-1.11	-1.50	-1.67
gen05	1.18	0.84	0.43	-0.45	-0.54	-0.24	-0.65	-1.15
gen06	-1.36	0.05	0.63	2.00	2.20	0.95	-0.08	-1.32
gen07	-0.23	1.55	2.52	0.34	-1.43	-0.93	-1.03	-0.47
gen08	-0.52	-0.49	0.72	1.31	2.10	0.82	-0.83	-0.66
gen09	0.38	2.25	2.00	0.81	-1.21	-1.25	-0.66	0.46
gen10	1.66	1.25	0.59	0.11	-0.18	-1.00	-1.05	-2.19
gen11	0.88	1.43	1.99	0.73	-0.26	-0.98	-0.50	0.26
gen12	-0.06	1.04	1.20	-0.15	-0.96	-0.90	-0.46	-0.43
gen13	-0.39	-0.27	0.88	1.99	2.13	0.48	-0.60	-0.18
gen14	0.09	1.67	2.43	0.32	-1.47	-0.99	-0.45	0.48
gen15	-0.57	-0.71	1.09	2.27	2.61	0.80	-0.24	-1.03
gen16	-1.40	-1.39	-0.81	0.39	1.02	1.57	1.82	2.14
gen17	-1.84	-0.80	-0.67	0.20	0.48	0.99	1.81	2.05
gen18	-1.32	-0.91	-0.99	0.30	0.39	0.24	1.35	1.63
gen19	-1.71	-0.84	0.20	-0.18	0.77	0.85	1.95	2.30
gen20	1.75	1.01	1.29	0.41	-0.26	-1.04	-1.74	-1.39
gen21	0.08	1.78	2.09	-0.21	-1.00	-1.13	-0.98	-0.34
gen22	0.04	2.27	2.64	0.17	-0.94	-2.01	-0.95	0.20
gen23	1.45	1.10	0.60	-0.17	-0.81	-0.79	-1.73	-2.69
gen24	1.17	0.80	-0.04	0.32	0.05	-0.45	-1.18	-1.26
gen25	-0.29	-0.51	0.34	1.54	2.11	0.76	-0.28	-0.30
gen26	-0.07	1.30	1.52	-0.09	-1.02	-1.39	-0.28	-0.12
gen27	-0.71	-0.52	-0.48	-0.33	0.91	0.93	0.99	1.70
gen28	-1.92	-1.28	-0.42	-0.30	0.52	0.49	1.11	1.46
gen29	-1.40	-0.88	0.57	2.01	1.90	0.29	-0.75	-1.86
gen30	-1.17	-0.83	-0.12	-0.01	0.50	0.93	1.34	2.14
gen31	0.07	1.61	2.27	0.33	-0.89	-0.45	0.15	0.37
gen32	1.65	1.69	0.29	0.45	-0.80	-1.06	-1.31	-1.87
gen33	-2.47	-1.46	-0.96	-0.06	0.57	0.23	1.14	1.05
gen34	2.10	1.25	1.38	-0.20	-0.81	-1.13	-0.88	-2.11
gen35	0.24	2.29	2.53	0.10	-1.16	-1.44	-0.58	0.10
gen36	-1.29	-1.10	-0.22	-0.38	0.43	1.12	1.17	1.88
gen37	-1.25	-0.27	0.64	1.92	1.81	0.50	-0.29	-1.28
gen38	-1.60	-1.15	-0.98	-0.08	0.70	1.31	2.17	2.59
gen39	0.14	2.36	2.14	-0.41	-1.56	-0.65	-0.66	0.18
gen40	-0.44	-0.90	-0.05	1.72	1.82	0.87	0.29	-0.68
gen41	-0.97	-0.29	-0.01	1.66	1.59	0.41	-0.23	-1.29
gen42	-0.44	-0.22	1.41	2.41	2.57	1.45	0.25	-0.50
gen43	1.69	1.03	0.29	0.57	-0.21	-0.34	-1.42	-1.87
gen44	-0.44	-0.33	0.46	1.57	2.12	0.70	-0.51	-0.46
gen45	-1.72	-0.66	-0.59	-0.30	0.21	0.83	1.17	1.97
gen46	1.48	1.34	0.68	-0.01	0.31	-0.79	-1.13	-1.25
gen47	1.76	0.63	0.42	0.01	-1.01	-0.78	-1.74	-1.95
gen48	-0.98	-0.60	0.81	2.32	2.79	0.62	-0.20	-0.96
gen49	-0.42	1.97	2.61	0.26	-1.37	-2.16	-0.90	0.16
gen50	-0.77	-0.40	0.13	2.06	2.51	0.95	-0.62	-1.53
gen51	1.15	0.83	0.09	0.29	-0.40	-0.48	-0.81	-1.77
gen52	1.05	2.12	2.04	0.84	-1.18	-0.97	0.03	0.32
gen53	-2.20	-1.18	-0.22	0.27	0.60	1.84	2.54	2.77
gen54	-1.63	-0.82	-0.00	-0.25	0.27	1.16	1.26	1.67
gen55	2.12	1.47	0.79	0.54	-0.44	-0.59	-1.15	-1.86
gen56	-1.06	0.17	0.28	1.31	1.83	0.63	0.15	-1.16
gen57	0.02	1.93	1.71	0.31	-1.08	-0.84	-0.70	0.31
gen58	-0.87	-0.03	0.76	1.68	1.85	0.84	0.22	-0.82
gen59	0.45	2.47	2.57	0.42	-0.95	-1.66	-0.84	0.06
gen60	1.03	0.87	0.76	0.01	0.23	-0.75	-0.07	-0.94";

const DATOS: Campo = Campo {
    id: "datos",
    etiqueta: ETIQUETA_TABLA,
    tipo: Tipo::Texto,
    ejemplo: TABLA_EJEMPLO,
};

const DISTANCIA: Campo = Campo {
    id: "distancia",
    etiqueta: "Distancia entre genes",
    tipo: Tipo::Opcion(MEDIDAS),
    ejemplo: MEDIDA_PEARSON,
};

const VECINOS: Campo = Campo {
    id: "k",
    etiqueta: "k (vecinos más cercanos de cada gen; 10 es lo usual)",
    tipo: Tipo::Entero,
    ejemplo: "10",
};

const RESOLUCION: Campo = Campo {
    id: "resolucion",
    etiqueta: "Resolución (más alta = comunidades más pequeñas; 1 es la modularidad clásica)",
    tipo: Tipo::Decimal,
    ejemplo: "1",
};

const SEMILLA: Campo = Campo {
    id: "semilla",
    etiqueta: "Semilla del azar (el mismo número repite el mismo resultado)",
    tipo: Tipo::Entero,
    ejemplo: "1",
};

const PERPLEJIDAD: Campo = Campo {
    id: "perplejidad",
    etiqueta: "Perplejidad (vecinos efectivos de cada gen; entre 5 y 50 es lo usual)",
    tipo: Tipo::Decimal,
    ejemplo: "15",
};

const EXAGERACION: Campo = Campo {
    id: "exageracion",
    etiqueta: "Exageración temprana (durante las primeras 250 iteraciones)",
    tipo: Tipo::Decimal,
    ejemplo: "12",
};

const ITERACIONES: Campo = Campo {
    id: "iteraciones",
    etiqueta: "Iteraciones del descenso de gradiente (máximo 5000)",
    tipo: Tipo::Entero,
    ejemplo: "1000",
};

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![
        Algoritmo {
            id: "agrupamiento.louvain",
            capitulo: "agrupamiento",
            nombre: "Comunidades de Louvain en un grafo de vecinos",
            resumen: "Une cada gen con sus k vecinos más cercanos, pesa cada unión por los vecinos que comparten y busca con Louvain las comunidades de mayor modularidad: a diferencia de k-means, no hay que decir cuántos grupos hay ni suponer nubes redondas.",
            campos: &[DATOS, DISTANCIA, VECINOS, RESOLUCION, SEMILLA],
            funcion: "comunidades::agrupar_por_comunidades",
            ejecutar: louvain,
            pasos: Some(pasos_louvain),
        },
        Algoritmo {
            id: "agrupamiento.tsne",
            capitulo: "agrupamiento",
            nombre: "Mapa t-SNE (exacto)",
            resumen: "Dibuja los genes en un plano cuidando quién es vecino de quién, no las distancias grandes: a diferencia de componentes principales, que conserva la varianza en líneas rectas, separa grupos que no caben en un plano; los colores son las comunidades de Louvain.",
            campos: &[
                DATOS,
                DISTANCIA,
                PERPLEJIDAD,
                EXAGERACION,
                ITERACIONES,
                VECINOS,
                RESOLUCION,
                SEMILLA,
            ],
            funcion: "comunidades::tsne",
            ejecutar: tsne,
            pasos: Some(pasos_tsne),
        },
    ]
}

struct Expresion {
    nombres: Vec<String>,
    columnas: Vec<String>,
    filas: Vec<Vec<f64>>,
}

impl Expresion {
    fn cantidad(&self) -> usize {
        self.filas.len()
    }

    fn dimension(&self) -> usize {
        self.columnas.len()
    }
}

struct Ajustes {
    medida: Medida,
    k: usize,
    resolucion: f64,
    semilla: u64,
}

fn leer_expresion(entrada: &Entrada, maximo_filas: usize) -> Result<Expresion, String> {
    let texto = entrada.texto("datos")?;
    let mut nombres = Vec::new();
    let mut filas: Vec<Vec<f64>> = Vec::new();
    let mut encabezado: Option<Vec<String>> = None;
    for (numero, linea) in texto.lines().enumerate() {
        let piezas: Vec<&str> = linea
            .split(|letra: char| letra.is_whitespace() || letra == ',' || letra == ';')
            .filter(|pieza| !pieza.is_empty())
            .collect();
        if piezas.is_empty() {
            continue;
        }
        let es_numero = |pieza: &str| pieza.parse::<f64>().is_ok();
        if filas.is_empty()
            && encabezado.is_none()
            && piezas.len() > 1
            && piezas.iter().all(|pieza| !es_numero(pieza))
        {
            encabezado = Some(piezas.iter().map(|pieza| pieza.to_string()).collect());
            continue;
        }
        let (nombre, valores) = if es_numero(piezas[0]) {
            (format!("gen{}", filas.len() + 1), &piezas[..])
        } else {
            (piezas[0].to_string(), &piezas[1..])
        };
        if valores.is_empty() {
            return Err(format!(
                "La línea {} («{}») no tiene números.",
                numero + 1,
                recortar(linea.trim(), 30)
            ));
        }
        if valores.len() > MAXIMO_COLUMNAS {
            return Err(format!(
                "La línea {} tiene {} números; el máximo es {MAXIMO_COLUMNAS} columnas.",
                numero + 1,
                valores.len()
            ));
        }
        let mut fila = Vec::with_capacity(valores.len());
        for pieza in valores {
            let valor: f64 = pieza.parse().map_err(|_| {
                format!(
                    "Línea {}: «{}» no es un número.",
                    numero + 1,
                    recortar(pieza, 30)
                )
            })?;
            if !valor.is_finite() || valor.abs() > MAXIMO_VALOR {
                return Err(format!(
                    "Línea {}: «{}» no es un número finito entre −10¹² y 10¹².",
                    numero + 1,
                    recortar(pieza, 30)
                ));
            }
            fila.push(valor);
        }
        if let Some(primera) = filas.first()
            && primera.len() != fila.len()
        {
            return Err(format!(
                "La línea {} ({nombre}) tiene {} números y la primera fila tiene {}: todas las filas deben tener la misma cantidad.",
                numero + 1,
                fila.len(),
                primera.len()
            ));
        }
        filas.push(fila);
        nombres.push(nombre);
        if filas.len() > maximo_filas {
            return Err(format!(
                "La tabla tiene más de {} filas, el máximo que acepta este método.",
                miles(maximo_filas)
            ));
        }
    }
    if filas.is_empty() {
        return Err("La tabla no tiene ninguna fila con números.".to_string());
    }
    let dimension = filas[0].len();
    let columnas = match encabezado {
        Some(titulos) if titulos.len() == dimension + 1 => titulos[1..].to_vec(),
        Some(titulos) if titulos.len() == dimension => titulos,
        _ => (1..=dimension)
            .map(|columna| format!("c{columna}"))
            .collect(),
    };
    Ok(Expresion {
        nombres,
        columnas,
        filas,
    })
}

fn recortar(texto: &str, maximo: usize) -> String {
    if texto.chars().count() <= maximo {
        return texto.to_string();
    }
    let corto: String = texto.chars().take(maximo).collect();
    format!("{corto}…")
}

fn leer_ajustes(entrada: &Entrada, tabla: &Expresion) -> Result<Ajustes, String> {
    let medida = if entrada.opcion("distancia", MEDIDAS)? == MEDIDA_EUCLIDEA {
        Medida::Euclidea
    } else {
        Medida::Pearson
    };
    if tabla.cantidad() < 3 {
        return Err(format!(
            "Hacen falta al menos 3 genes y la tabla tiene {}.",
            tabla.cantidad()
        ));
    }
    let k = entrada.natural("k")?;
    let tope = MAXIMO_VECINOS.min(tabla.cantidad() - 1);
    if k == 0 || k > tope {
        return Err(format!(
            "k debe estar entre 1 y {tope} (como mucho {MAXIMO_VECINOS} y menos que la cantidad de genes)."
        ));
    }
    let resolucion = entrada.decimal("resolucion")?;
    if resolucion <= 0.0 || resolucion > MAXIMO_RESOLUCION {
        return Err(format!(
            "La resolución debe ser mayor que 0 y como mucho {MAXIMO_RESOLUCION}."
        ));
    }
    let semilla = entrada.natural("semilla")? as u64;
    Ok(Ajustes {
        medida,
        k,
        resolucion,
        semilla,
    })
}

fn leer_opciones_tsne(
    entrada: &Entrada,
    tabla: &Expresion,
    semilla: u64,
) -> Result<OpcionesTsne, String> {
    let n = tabla.cantidad();
    if n < 4 {
        return Err(format!(
            "t-SNE necesita al menos 4 genes y la tabla tiene {n}."
        ));
    }
    let perplejidad = entrada.decimal("perplejidad")?;
    let tope = (n - 1) as f64 / 3.0;
    if !(1.0..=tope).contains(&perplejidad) {
        return Err(format!(
            "Con {} genes la perplejidad debe estar entre 1 y {}: tres veces la perplejidad no puede pasar de los genes menos uno.",
            miles(n),
            decimal(tope, 1)
        ));
    }
    let exageracion = entrada.decimal("exageracion")?;
    if !(1.0..=MAXIMO_EXAGERACION).contains(&exageracion) {
        return Err(format!(
            "La exageración temprana debe estar entre 1 y {MAXIMO_EXAGERACION}."
        ));
    }
    let iteraciones = entrada.natural("iteraciones")?;
    if iteraciones == 0 || iteraciones > MAXIMO_ITERACIONES {
        return Err(format!(
            "Las iteraciones deben estar entre 1 y {}.",
            miles(MAXIMO_ITERACIONES)
        ));
    }
    let cuentas = (n * n) as f64 * iteraciones as f64;
    if cuentas > MAXIMO_CUENTAS_TSNE {
        return Err(format!(
            "t-SNE exacto compara todos los pares de genes en cada iteración: {} genes × {} iteraciones son demasiadas cuentas. Baja las iteraciones o la cantidad de genes.",
            miles(n),
            miles(iteraciones)
        ));
    }
    Ok(OpcionesTsne {
        perplejidad,
        exageracion,
        iteraciones,
        semilla,
    })
}

fn segun_cantidad<'a>(cantidad: usize, singular: &'a str, plural: &'a str) -> &'a str {
    if cantidad == 1 { singular } else { plural }
}

fn unir(nombres: &[String]) -> String {
    match nombres.len() {
        0 => String::new(),
        1 => nombres[0].clone(),
        cantidad => format!(
            "{} y {}",
            nombres[..cantidad - 1].join(", "),
            nombres[cantidad - 1]
        ),
    }
}

fn algunos(nombres: &[String], maximo: usize) -> String {
    if nombres.len() <= maximo {
        return unir(nombres);
    }
    format!("{}…", nombres[..maximo].join(", "))
}

fn muestra_pareja(total: usize, maximo: usize) -> Vec<usize> {
    if total <= maximo || maximo < 2 {
        return (0..total).collect();
    }
    (0..maximo)
        .map(|paso| paso * (total - 1) / (maximo - 1))
        .collect()
}

fn tamanos(comunidad: &[usize], cantidad: usize) -> Vec<usize> {
    let mut cuentas = vec![0; cantidad];
    for &c in comunidad {
        cuentas[c] += 1;
    }
    cuentas
}

fn miembros(tabla: &Expresion, comunidad: &[usize], elegida: usize) -> Vec<String> {
    comunidad
        .iter()
        .enumerate()
        .filter(|&(_, &c)| c == elegida)
        .map(|(gen_, _)| tabla.nombres[gen_].clone())
        .collect()
}

fn contar_aristas(grafo: &Grafo) -> usize {
    grafo.vecinos.iter().map(Vec::len).sum::<usize>() / 2
}

fn tabla_de_comunidades(tabla: &Expresion, louvain: &Louvain) -> Bloque {
    let cuentas = tamanos(&louvain.comunidad, louvain.cantidad);
    let filas = (0..louvain.cantidad)
        .map(|c| {
            vec![
                format!("{}", c + 1),
                miles(cuentas[c]),
                algunos(
                    &miembros(tabla, &louvain.comunidad, c),
                    NOMBRES_POR_COMUNIDAD,
                ),
            ]
        })
        .collect();
    Bloque::tabla("Comunidades", &["Comunidad", "Genes", "Quiénes"], filas)
}

fn calor_por_comunidad(tabla: &Expresion, comunidad: &[usize]) -> Bloque {
    let mut orden: Vec<usize> = (0..tabla.cantidad()).collect();
    orden.sort_by_key(|&gen_| (comunidad[gen_], gen_));
    let elegidos: Vec<usize> = muestra_pareja(orden.len(), MAXIMO_FILAS_DE_CALOR)
        .into_iter()
        .map(|posicion| orden[posicion])
        .collect();
    let filas = elegidos
        .iter()
        .map(|&gen_| {
            format!(
                "{} · comunidad {}",
                tabla.nombres[gen_],
                comunidad[gen_] + 1
            )
        })
        .collect();
    let valores = elegidos
        .iter()
        .map(|&gen_| tabla.filas[gen_].clone())
        .collect();
    Bloque::calor(
        "Expresión ordenada por comunidad",
        filas,
        tabla.columnas.clone(),
        valores,
    )
}

fn tabla_de_niveles(louvain: &Louvain) -> Bloque {
    let filas = louvain
        .niveles
        .iter()
        .enumerate()
        .map(|(numero, nivel)| {
            vec![
                format!("{}", numero + 1),
                miles(nivel.nodos),
                miles(nivel.movimientos),
                miles(nivel.comunidades),
                decimal(nivel.modularidad, 4),
            ]
        })
        .collect();
    Bloque::tabla(
        "Niveles de Louvain (mover nodos y luego juntar cada comunidad en un nodo)",
        &[
            "Nivel",
            "Nodos",
            "Movimientos",
            "Comunidades",
            "Modularidad",
        ],
        filas,
    )
}

struct Comparacion {
    rand: f64,
    grupos: Vec<usize>,
    k: usize,
    milisegundos: f64,
}

fn comparar_con_k_medias(
    tabla: &Expresion,
    louvain: &Louvain,
    semilla: u64,
) -> Result<Option<Comparacion>, String> {
    let k = louvain.cantidad;
    if !(2..=MAXIMO_K_MEDIAS).contains(&k) {
        return Ok(None);
    }
    let inicio = Instant::now();
    let resultado = agrupamiento::k_medias(
        &tabla.filas,
        k,
        Inicio::KMediasMasMas(semilla),
        ITERACIONES_K_MEDIAS,
    )?;
    let milisegundos = inicio.elapsed().as_secs_f64() * 1000.0;
    Ok(Some(Comparacion {
        rand: comunidades::indice_rand_ajustado(&louvain.comunidad, &resultado.grupos),
        grupos: resultado.grupos,
        k,
        milisegundos,
    }))
}

fn tabla_cruzada(louvain: &Louvain, comparacion: &Comparacion) -> Bloque {
    let mut cruce = vec![vec![0usize; comparacion.k]; louvain.cantidad];
    for (&c, &g) in louvain.comunidad.iter().zip(&comparacion.grupos) {
        cruce[c][g] += 1;
    }
    let mut columnas = vec!["Comunidad".to_string()];
    columnas.extend((1..=comparacion.k).map(|g| format!("k-means {g}")));
    let referencias: Vec<&str> = columnas.iter().map(String::as_str).collect();
    let filas = cruce
        .iter()
        .enumerate()
        .map(|(c, fila)| {
            let mut celdas = vec![format!("{}", c + 1)];
            celdas.extend(fila.iter().map(|&cuenta| miles(cuenta)));
            celdas
        })
        .collect();
    Bloque::tabla(
        "Genes de cada comunidad en cada grupo de k-means",
        &referencias,
        filas,
    )
}

fn tabla_revuelta(tabla: &Expresion, semilla: u64) -> Vec<Vec<f64>> {
    let mut azar = Azar::nuevo(semilla.wrapping_add(SEMILLA_REVUELTA));
    let columnas: Vec<Vec<f64>> = (0..tabla.dimension())
        .map(|columna| {
            comunidades::barajar(tabla.cantidad(), &mut azar)
                .iter()
                .map(|&origen| tabla.filas[origen][columna])
                .collect()
        })
        .collect();
    (0..tabla.cantidad())
        .map(|fila| columnas.iter().map(|columna| columna[fila]).collect())
        .collect()
}

fn modularidad_revuelta(tabla: &Expresion, ajustes: &Ajustes) -> Option<f64> {
    let filas = tabla_revuelta(tabla, ajustes.semilla);
    comunidades::agrupar_por_comunidades(
        &filas,
        ajustes.medida,
        ajustes.k,
        ajustes.resolucion,
        ajustes.semilla,
    )
    .ok()
    .map(|resultado| resultado.louvain.modularidad)
}

fn supera_al_azar(modularidad: f64, revuelta: Option<f64>) -> bool {
    revuelta.is_none_or(|nula| modularidad - nula > MARGEN_SOBRE_AZAR)
}

fn nombre_de_medida(medida: Medida) -> &'static str {
    match medida {
        Medida::Pearson => "Pearson (1 − r)",
        Medida::Euclidea => "euclídea",
    }
}

fn louvain(entrada: &Entrada) -> Result<Salida, String> {
    let tabla = leer_expresion(entrada, MAXIMO_FILAS)?;
    let ajustes = leer_ajustes(entrada, &tabla)?;
    let inicio = Instant::now();
    let resultado = comunidades::agrupar_por_comunidades(
        &tabla.filas,
        ajustes.medida,
        ajustes.k,
        ajustes.resolucion,
        ajustes.semilla,
    )?;
    let milisegundos = inicio.elapsed().as_secs_f64() * 1000.0;
    let comparacion = comparar_con_k_medias(&tabla, &resultado.louvain, ajustes.semilla)?;
    let revuelta = modularidad_revuelta(&tabla, &ajustes);
    let louvain = &resultado.louvain;
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Comunidades", miles(louvain.cantidad)));
    salida.agregar(Bloque::valor(
        "Modularidad Q",
        decimal(louvain.modularidad, 3),
    ));
    if let Some(nula) = revuelta {
        salida.agregar(Bloque::valor(
            "Q con los valores revueltos (sin estructura)",
            decimal(nula, 3),
        ));
    }
    salida.agregar(Bloque::valor(
        "Grafo de vecinos compartidos",
        format!(
            "{} genes, {} uniones, peso total {}",
            miles(tabla.cantidad()),
            miles(contar_aristas(&resultado.grafo)),
            decimal(resultado.grafo.peso_total, 2)
        ),
    ));
    let cuentas = tamanos(&louvain.comunidad, louvain.cantidad);
    let mostradas = louvain.cantidad.min(MAXIMO_BARRAS);
    salida.agregar(Bloque::barras(
        "Tamaño de cada comunidad",
        "Comunidad",
        "Genes",
        (1..=mostradas).map(|c| c.to_string()).collect(),
        cuentas[..mostradas].iter().map(|&c| c as f64).collect(),
    ));
    salida.agregar(tabla_de_comunidades(&tabla, louvain));
    salida.agregar(calor_por_comunidad(&tabla, &louvain.comunidad));
    salida.agregar(tabla_de_niveles(louvain));
    if let Some(comparacion) = &comparacion {
        salida.agregar(Bloque::valor(
            &format!("Rand ajustado contra k-means con k = {}", comparacion.k),
            decimal(comparacion.rand, 3),
        ));
        if louvain.cantidad <= MAXIMO_CRUCE {
            salida.agregar(tabla_cruzada(louvain, comparacion));
        }
    }
    let tiempos = match &comparacion {
        Some(comparacion) => format!(
            "vecinos + grafo + Louvain {} ms · k-means {} ms",
            decimal(milisegundos, 1),
            decimal(comparacion.milisegundos, 1)
        ),
        None => format!("vecinos + grafo + Louvain {} ms", decimal(milisegundos, 1)),
    };
    salida.agregar(Bloque::valor("Tiempo", tiempos));
    explicar_louvain(
        &mut salida,
        &tabla,
        &ajustes,
        &resultado,
        comparacion.as_ref(),
        revuelta,
    );
    Ok(salida)
}

fn explicar_louvain(
    salida: &mut Salida,
    tabla: &Expresion,
    ajustes: &Ajustes,
    resultado: &Comunidades,
    comparacion: Option<&Comparacion>,
    revuelta: Option<f64>,
) {
    let louvain = &resultado.louvain;
    if resultado.grafo.peso_total <= 0.0 {
        salida.explicar(format!(
            "Ningún par de genes comparte suficientes vecinos (se descartan uniones con un índice de Jaccard menor que 1/15), así que no hay grafo que partir y cada gen quedó solo. Sube k (ahora {}).",
            ajustes.k
        ));
        return;
    }
    let q = louvain.modularidad;
    if q > MODULARIDAD_CLARA && supera_al_azar(q, revuelta) {
        let contraste = revuelta.map_or(String::new(), |nula| {
            format!(
                " y bastante más que el {} que da este mismo método con los mismos valores revueltos (sin estructura)",
                decimal(nula, 3)
            )
        });
        salida.explicar(format!(
            "La modularidad es {}: mayor que 0.3, el umbral usual para decir que hay estructura de comunidades clara{contraste}. Dentro de cada comunidad hay mucho más peso de uniones que el que tendría un grafo al azar con los mismos grados.",
            decimal(q, 3)
        ));
    } else if q > MODULARIDAD_CLARA {
        salida.explicar(format!(
            "La modularidad es {}: pasa de 0.3, pero con los mismos valores revueltos (sin ninguna estructura) este método también da {}. Un grafo de vecinos siempre tiene comunidades aparentes, aunque los datos sean ruido, así que este número no prueba que los grupos sean reales.",
            decimal(q, 3),
            decimal(revuelta.unwrap_or(0.0), 3)
        ));
    } else {
        salida.explicar(format!(
            "La modularidad es {}: no pasa de 0.3, así que la estructura es débil; los genes no forman grupos bien separados y las fronteras entre comunidades son poco confiables.",
            decimal(q, 3)
        ));
    }
    if (ajustes.resolucion - 1.0).abs() > f64::EPSILON {
        salida.explicar(format!(
            "Con resolución {} la modularidad es la generalizada, que no se compara directamente con el 0.3 de la resolución 1.",
            decimal(ajustes.resolucion, 2)
        ));
    }
    let cuentas = tamanos(&louvain.comunidad, louvain.cantidad);
    let pequenas = cuentas.iter().filter(|&&cuenta| cuenta <= 2).count();
    let mut texto = format!(
        "Con resolución {} salieron {} {}. Subir la resolución pide comunidades más pequeñas y numerosas; bajarla, más grandes: no hay un número correcto, conviene probar varias y quedarse con las comunidades que no cambian.",
        decimal(ajustes.resolucion, 2),
        miles(louvain.cantidad),
        segun_cantidad(louvain.cantidad, "comunidad", "comunidades")
    );
    if pequenas > 0 {
        texto.push_str(&format!(
            " {} de 1 o 2 genes: suelen ser genes atípicos o con pocos vecinos compartidos.",
            match pequenas {
                1 => "Hay 1 comunidad".to_string(),
                _ => format!("Hay {} comunidades", miles(pequenas)),
            }
        ));
    }
    salida.explicar(texto);
    match comparacion {
        Some(comparacion) if comparacion.rand >= ACUERDO_ALTO => salida.explicar(format!(
            "El índice de Rand ajustado contra «{NOMBRE_LLOYD}» con k = {} es {}: prácticamente la misma partición. Cuando coinciden, los grupos son nubes compactas; la ventaja de Louvain es que encontró solo el número de grupos, que a k-means hubo que dárselo.",
            comparacion.k,
            decimal(comparacion.rand, 3)
        )),
        Some(comparacion) if comparacion.rand >= ACUERDO_MEDIO => salida.explicar(format!(
            "El índice de Rand ajustado contra «{NOMBRE_LLOYD}» con k = {} es {}: coinciden en buena parte. Mira la tabla cruzada para ver qué comunidades reparte k-means de otra forma.",
            comparacion.k,
            decimal(comparacion.rand, 3)
        )),
        Some(comparacion) => salida.explicar(format!(
            "El índice de Rand ajustado contra «{NOMBRE_LLOYD}» con k = {} es {}: las particiones difieren mucho. k-means busca nubes redondas alrededor de centros y Louvain sigue la conectividad entre vecinos, así que detecta grupos alargados o de densidad distinta que k-means parte o mezcla.",
            comparacion.k,
            decimal(comparacion.rand, 3)
        )),
        None => {}
    }
    if let Some(comparacion) = comparacion
        && ajustes.medida == Medida::Pearson
        && comparacion.rand < ACUERDO_ALTO
    {
        salida.explicar(
            "Ojo: k-means usó los valores tal cual, con distancia euclídea, mientras que el grafo midió la distancia de Pearson (solo la forma del perfil): parte de la diferencia viene de medir distinto y no solo del método.",
        );
    }
    salida.explicar(format!(
        "El grafo usa k = {} vecinos con distancia {}. Un k chico da un grafo ralo y comunidades pequeñas; uno grande une grupos vecinos. Louvain recorre los nodos en un orden al azar: si al cambiar la semilla las comunidades cambian mucho, no te fíes de esas fronteras.",
        ajustes.k,
        nombre_de_medida(ajustes.medida)
    ));
    salida.explicar(
        "Louvain (Blondel, Guillaume, Lambiotte y Lefebvre, 2008) puede dejar comunidades mal conectadas por dentro; Leiden (Traag, Waltman y van Eck, 2019) lo refina con un paso extra que lo garantiza, y es el que recomienda hoy Scanpy.",
    );
    if tabla.cantidad() < 3 * ajustes.k {
        salida.explicar(format!(
            "Hay pocos genes ({}) para k = {}: cada gen ve a casi todos los demás y el grafo pierde detalle. Prueba un k menor.",
            miles(tabla.cantidad()),
            ajustes.k
        ));
    }
}

struct Mapa {
    tsne: Tsne,
    louvain: Louvain,
    conservacion: f64,
    conservacion_pca: Option<f64>,
    islas: f64,
    revuelta: Option<f64>,
    milisegundos: f64,
}

fn calcular_mapa(
    tabla: &Expresion,
    ajustes: &Ajustes,
    opciones: OpcionesTsne,
) -> Result<Mapa, String> {
    let inicio = Instant::now();
    let mapa = comunidades::tsne(&tabla.filas, ajustes.medida, opciones)?;
    let milisegundos = inicio.elapsed().as_secs_f64() * 1000.0;
    let grupos = comunidades::agrupar_por_comunidades(
        &tabla.filas,
        ajustes.medida,
        ajustes.k,
        ajustes.resolucion,
        ajustes.semilla,
    )?;
    let filas = comunidades::filas_para_medir(&tabla.filas, ajustes.medida);
    let conservacion = comunidades::conservacion_de_vecinos(&filas, &mapa.coordenadas, ajustes.k);
    let conservacion_pca = if tabla.dimension() >= 2 {
        agrupamiento::componentes_principales(&filas)
            .ok()
            .map(|pca| comunidades::conservacion_de_vecinos(&filas, &pca.coordenadas, ajustes.k))
    } else {
        None
    };
    let plano: Vec<Vec<f64>> = mapa
        .coordenadas
        .iter()
        .map(|punto| punto.to_vec())
        .collect();
    let cercanos = comunidades::vecinos_mas_cercanos(&plano, 1);
    let iguales = cercanos
        .iter()
        .enumerate()
        .filter(|(gen_, lista)| {
            grupos.louvain.comunidad[lista[0]] == grupos.louvain.comunidad[*gen_]
        })
        .count();
    Ok(Mapa {
        islas: iguales as f64 / tabla.cantidad() as f64,
        revuelta: modularidad_revuelta(tabla, ajustes),
        tsne: mapa,
        louvain: grupos.louvain,
        conservacion,
        conservacion_pca,
        milisegundos,
    })
}

fn tsne(entrada: &Entrada) -> Result<Salida, String> {
    let tabla = leer_expresion(entrada, MAXIMO_FILAS_TSNE)?;
    let ajustes = leer_ajustes(entrada, &tabla)?;
    let opciones = leer_opciones_tsne(entrada, &tabla, ajustes.semilla)?;
    let mapa = calcular_mapa(&tabla, &ajustes, opciones)?;
    let mut salida = Salida::default();
    let puntos = mapa
        .tsne
        .coordenadas
        .iter()
        .enumerate()
        .map(|(gen_, punto)| Punto {
            x: punto[0],
            y: punto[1],
            grupo: mapa.louvain.comunidad[gen_],
            nombre: tabla.nombres[gen_].clone(),
        })
        .collect();
    salida.agregar(Bloque::dispersion(
        "Mapa t-SNE (color = comunidad de Louvain)",
        "t-SNE 1",
        "t-SNE 2",
        puntos,
    ));
    salida.agregar(Bloque::valor(
        "Divergencia KL final",
        decimal(mapa.tsne.divergencia, 4),
    ));
    let curva: Vec<[f64; 2]> = mapa
        .tsne
        .historia
        .iter()
        .map(|&(iteracion, valor)| [iteracion as f64, valor])
        .collect();
    let mut marcas = Vec::new();
    if opciones.iteraciones > ITERACIONES_EXAGERADAS {
        marcas.push(Marca {
            x: ITERACIONES_EXAGERADAS as f64,
            texto: "fin de la exageración".to_string(),
        });
    }
    salida.agregar(Bloque::linea(
        "Divergencia KL durante el descenso",
        "Iteración",
        "KL",
        &curva,
        marcas,
    ));
    salida.agregar(Bloque::valor(
        &format!("Vecinos conservados (de {} por gen)", ajustes.k),
        conservados(&mapa),
    ));
    salida.agregar(Bloque::valor(
        "Comunidades de Louvain (colores)",
        format!(
            "{} · modularidad {}",
            miles(mapa.louvain.cantidad),
            decimal(mapa.louvain.modularidad, 3)
        ),
    ));
    let mut sigmas = mapa.tsne.sigmas.clone();
    sigmas.sort_by(f64::total_cmp);
    salida.agregar(Bloque::valor(
        "Ancho σ de cada gaussiana (mínimo · mediana · máximo)",
        format!(
            "{} · {} · {}",
            decimal(sigmas[0], 3),
            decimal(sigmas[sigmas.len() / 2], 3),
            decimal(sigmas[sigmas.len() - 1], 3)
        ),
    ));
    salida.agregar(Bloque::valor(
        "Tasa de aprendizaje y tiempo",
        format!(
            "{} · {} ms",
            decimal(mapa.tsne.tasa, 1),
            decimal(mapa.milisegundos, 0)
        ),
    ));
    explicar_tsne(&mut salida, &tabla, &ajustes, &opciones, &mapa);
    Ok(salida)
}

fn conservados(mapa: &Mapa) -> String {
    match mapa.conservacion_pca {
        Some(pca) => format!(
            "t-SNE {} · componentes principales {}",
            decimal(mapa.conservacion * 100.0, 1) + " %",
            decimal(pca * 100.0, 1) + " %"
        ),
        None => format!("t-SNE {} %", decimal(mapa.conservacion * 100.0, 1)),
    }
}

fn explicar_tsne(
    salida: &mut Salida,
    tabla: &Expresion,
    ajustes: &Ajustes,
    opciones: &OpcionesTsne,
    mapa: &Mapa,
) {
    salida.explicar(format!(
        "La divergencia KL final es {}: mide cuánta información sobre quién es vecino de quién se pierde al pasar al plano. Solo sirve para comparar corridas con los mismos datos y la misma perplejidad; más baja es mejor.",
        decimal(mapa.tsne.divergencia, 4)
    ));
    salida.explicar(
        "En t-SNE las distancias entre grupos y el tamaño de cada grupo no significan nada: dos islas lejanas no son más distintas que dos cercanas. Lo que sí vale es qué genes quedan juntos.",
    );
    match mapa.conservacion_pca {
        Some(pca) if mapa.conservacion > pca + VENTAJA_CLARA => salida.explicar(format!(
            "De los {} vecinos más cercanos de cada gen, el mapa conserva en promedio el {} %, y la vista de «{NOMBRE_COMPONENTES}» solo el {} %: t-SNE muestra mucho mejor la vecindad local.",
            ajustes.k,
            decimal(mapa.conservacion * 100.0, 1),
            decimal(pca * 100.0, 1)
        )),
        Some(pca) => salida.explicar(format!(
            "De los {} vecinos más cercanos de cada gen, el mapa conserva en promedio el {} % y la vista de «{NOMBRE_COMPONENTES}» el {} %: con estos datos las dos vistas cuentan casi lo mismo, y componentes principales además conserva las distancias grandes.",
            ajustes.k,
            decimal(mapa.conservacion * 100.0, 1),
            decimal(pca * 100.0, 1)
        )),
        None => salida.explicar(format!(
            "De los {} vecinos más cercanos de cada gen, el mapa conserva en promedio el {} %.",
            ajustes.k,
            decimal(mapa.conservacion * 100.0, 1)
        )),
    }
    if mapa.islas >= ISLAS_CLARAS && !supera_al_azar(mapa.louvain.modularidad, mapa.revuelta) {
        salida.explicar(format!(
            "El {} % de los genes tiene como vecino más cercano en el mapa a uno de su misma comunidad, pero la modularidad ({}) no supera a la de los mismos valores revueltos ({}): t-SNE y Louvain dibujan islas hasta con ruido, así que no son prueba de grupos reales.",
            decimal(mapa.islas * 100.0, 0),
            decimal(mapa.louvain.modularidad, 3),
            decimal(mapa.revuelta.unwrap_or(0.0), 3)
        ));
    } else if mapa.islas >= ISLAS_CLARAS {
        salida.explicar(format!(
            "El {} % de los genes tiene como vecino más cercano en el mapa a uno de su misma comunidad de Louvain: las comunidades se ven como islas, y dos métodos independientes coinciden.",
            decimal(mapa.islas * 100.0, 0)
        ));
    } else {
        salida.explicar(format!(
            "Solo el {} % de los genes tiene como vecino más cercano en el mapa a uno de su misma comunidad: los colores se mezclan, señal de que las comunidades no están bien separadas o de que conviene cambiar la perplejidad.",
            decimal(mapa.islas * 100.0, 0)
        ));
    }
    let mut perplejidad = format!(
        "La perplejidad {} es el número efectivo de vecinos que mira cada gen. Lo usual va de 5 a 50 y debe ser menor que un tercio de los genes.",
        decimal(opciones.perplejidad, 1)
    );
    if opciones.perplejidad < PERPLEJIDAD_BAJA {
        perplejidad.push_str(
            " La tuya es baja: el mapa tiende a partirse en grupitos artificiales; súbela y mira qué islas sobreviven.",
        );
    } else if opciones.perplejidad > PERPLEJIDAD_ALTA {
        perplejidad.push_str(
            " La tuya es alta: el mapa mira estructura más global y puede juntar grupos pequeños.",
        );
    } else if opciones.perplejidad > tabla.cantidad() as f64 / 6.0 {
        perplejidad.push_str(
            " Para tan pocos genes ya es alta: cada gen mira a una parte grande de la tabla.",
        );
    }
    salida.explicar(perplejidad);
    if opciones.iteraciones <= ITERACIONES_EXAGERADAS {
        salida.explicar(format!(
            "Todas las iteraciones fueron con exageración temprana (las primeras {ITERACIONES_EXAGERADAS}): el mapa no tuvo tiempo de relajarse. Usa 1000."
        ));
    } else if let Some((caida, desde)) = caida_final(&mapa.tsne.historia, opciones.iteraciones) {
        if caida > CAIDA_FINAL_NOTABLE {
            salida.explicar(format!(
                "Desde la iteración {} hasta la última la KL todavía bajó un {} %: el descenso no terminó; sube las iteraciones.",
                miles(desde),
                decimal(caida * 100.0, 1)
            ));
        } else {
            salida.explicar(format!(
                "Desde la iteración {} hasta la última la KL casi no cambió: el descenso se estabilizó.",
                miles(desde)
            ));
        }
    }
    salida.explicar(
        "Es t-SNE exacto (van der Maaten y Hinton, 2008): compara todos los pares, así que es fiel pero crece con el cuadrado de los genes. Para decenas de miles de puntos se usan aproximaciones (Barnes-Hut, FIt-SNE) o UMAP; con la semilla fija el mapa se repite igual.",
    );
}

fn caida_final(historia: &[(usize, f64)], iteraciones: usize) -> Option<(f64, usize)> {
    let final_ = historia.last()?.1;
    let corte = iteraciones - iteraciones / 10;
    let &(desde, antes) = historia
        .iter()
        .rev()
        .find(|&&(iteracion, _)| iteracion <= corte)?;
    if antes <= 0.0 {
        return None;
    }
    Some(((antes - final_) / antes, desde))
}

fn corto(valor: f64) -> String {
    if valor != 0.0 && (valor.abs() < 1e-3 || valor.abs() >= 1e5) {
        format!("{valor:.3e}")
    } else {
        decimal(valor, 4)
    }
}

fn perfil(fila: &[f64]) -> String {
    let partes: Vec<String> = fila.iter().take(8).map(|v| decimal(*v, 2)).collect();
    if fila.len() > 8 {
        format!("{} …", partes.join(" "))
    } else {
        partes.join(" ")
    }
}

fn candidatas(
    grafo: &Grafo,
    particion: &Particion,
    nodo: usize,
    resolucion: f64,
) -> Vec<(usize, f64, f64)> {
    let m = grafo.peso_total;
    let grado = grafo.grados[nodo];
    let desde = particion.comunidad[nodo];
    let mut pesos: BTreeMap<usize, f64> = BTreeMap::new();
    pesos.insert(desde, 0.0);
    for &(otro, peso) in &grafo.vecinos[nodo] {
        *pesos.entry(particion.comunidad[otro]).or_insert(0.0) += peso;
    }
    let grado_sin = |c: usize| {
        if c == desde {
            particion.grado_de_comunidad[c] - grado
        } else {
            particion.grado_de_comunidad[c]
        }
    };
    let quedarse = pesos[&desde] - resolucion * grado_sin(desde) * grado / (2.0 * m);
    let mut lista: Vec<(usize, f64, f64)> = pesos
        .iter()
        .map(|(&c, &peso)| {
            let valor = peso - resolucion * grado_sin(c) * grado / (2.0 * m);
            (c, peso, (valor - quedarse) / m)
        })
        .collect();
    lista.sort_by(|a, b| b.2.total_cmp(&a.2).then(a.0.cmp(&b.0)));
    lista
}

fn pasos_louvain(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let tabla = leer_expresion(entrada, MAXIMO_FILAS)?;
    let ajustes = leer_ajustes(entrada, &tabla)?;
    let n = tabla.cantidad();
    let k = ajustes.k;
    let nombres = &tabla.nombres;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        Pizarra::nueva()
            .variable("puntos", format!("{} genes × {} columnas", miles(n), tabla.dimension()))
            .variable("medida", nombre_de_medida(ajustes.medida))
            .variable("k", k)
            .variable("resolucion", decimal(ajustes.resolucion, 2))
            .paso(
                "revisar_vecinos(puntos, k)?;",
                format!(
                    "Tienes {} genes con {} mediciones cada uno. La meta: unir cada gen con sus {k} vecinos más parecidos y buscar en ese grafo comunidades, grupos con muchas uniones por dentro y pocas hacia afuera.",
                    miles(n),
                    tabla.dimension()
                ),
            ),
    );
    let filas = comunidades::filas_para_medir(&tabla.filas, ajustes.medida);
    let nota_filas = match ajustes.medida {
        Medida::Pearson => format!(
            "Con Pearson, a cada gen se le resta su promedio y se divide por su largo: {} pasa de «{}» a «{}». Así, la distancia al cuadrado entre dos genes es 2 · (1 − r) y solo cuenta la forma del perfil, no su altura.",
            nombres[0],
            perfil(&tabla.filas[0]),
            perfil(&filas[0])
        ),
        Medida::Euclidea => format!(
            "Con distancia euclídea los genes se usan tal cual: {} es «{}». Cuenta la forma y también la altura de cada perfil.",
            nombres[0],
            perfil(&tabla.filas[0])
        ),
    };
    narracion.agregar(
        Pizarra::nueva()
            .variable("medida", nombre_de_medida(ajustes.medida))
            .variable("filas", format!("{} genes", miles(n)))
            .paso("let filas = filas_para_medir(puntos, medida);", nota_filas),
    );
    let punto = 0;
    let mut distancias = comunidades::distancias_desde(&filas, punto);
    let mut vistas: Vec<Vec<String>> = Vec::new();
    for (posicion, &(cuadrado, otro)) in distancias.iter().enumerate() {
        if vistas.len() < FILAS_DE_DISTANCIAS {
            vistas.push(vec![
                nombres[otro].clone(),
                decimal(comunidades::distancia_medida(cuadrado, ajustes.medida), 4),
            ]);
        }
        if posicion < DISTANCIAS_DETALLADAS {
            let ultima = vistas.len() - 1;
            narracion.agregar(
                Pizarra::nueva()
                    .variable("punto", &nombres[punto])
                    .variable("otro", &nombres[otro])
                    .variable("distancia", decimal(comunidades::distancia_medida(cuadrado, ajustes.medida), 4))
                    .tabla(
                        &format!("Distancias desde {}", nombres[punto]),
                        &["otro", "distancia"],
                        vistas.clone(),
                        &[(ultima, 1)],
                    )
                    .paso(
                        "distancias.push((distancia_al_cuadrado(&filas[punto], fila), otro));",
                        if posicion == 0 {
                            format!(
                                "Empezamos por {}: medimos su distancia a cada otro gen. A {} le toca {}.",
                                nombres[punto],
                                nombres[otro],
                                decimal(comunidades::distancia_medida(cuadrado, ajustes.medida), 4)
                            )
                        } else {
                            format!(
                                "Distancia de {} a {}: {}.",
                                nombres[punto],
                                nombres[otro],
                                decimal(comunidades::distancia_medida(cuadrado, ajustes.medida), 4)
                            )
                        },
                    ),
            );
        }
    }
    narracion.agregar(
        Pizarra::nueva()
            .variable("punto", &nombres[punto])
            .variable("otro", "…")
            .variable("distancia", "…")
            .tabla(
                &format!("Distancias desde {}", nombres[punto]),
                &["otro", "distancia"],
                vistas.clone(),
                &[],
            )
            .paso(
                "for (otro, fila) in filas.iter().enumerate() {",
                format!(
                    "Hacemos lo mismo con los {} genes restantes: {} tiene ya {} distancias, una por cada otro gen.",
                    miles(n.saturating_sub(1 + DISTANCIAS_DETALLADAS)),
                    nombres[punto],
                    miles(n - 1)
                ),
            ),
    );
    let lista = comunidades::quedarse_con_los_k(&mut distancias, k);
    let filas_vecinos: Vec<Vec<String>> = distancias
        .iter()
        .take(FILAS_DE_DISTANCIAS)
        .map(|&(cuadrado, otro)| {
            vec![
                nombres[otro].clone(),
                decimal(comunidades::distancia_medida(cuadrado, ajustes.medida), 4),
            ]
        })
        .collect();
    let resaltadas: Vec<(usize, usize)> = (0..filas_vecinos.len()).map(|fila| (fila, 0)).collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("punto", &nombres[punto])
            .variable("k", k)
            .tabla(
                &format!("Los {k} más cercanos a {}", nombres[punto]),
                &["otro", "distancia"],
                filas_vecinos,
                &resaltadas,
            )
            .paso(
                "distancias.select_nth_unstable_by(k - 1, mas_cerca);",
                format!(
                    "No hace falta ordenar las {} distancias: select_nth_unstable_by deja las {k} menores al principio en tiempo lineal, y solo esas se ordenan.",
                    miles(n - 1)
                ),
            ),
    );
    let nombres_lista: Vec<String> = lista.iter().map(|&otro| nombres[otro].clone()).collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("punto", &nombres[punto])
            .variable("k", k)
            .lista(
                &format!("Vecinos de {}", nombres[punto]),
                nombres_lista.clone(),
                &(0..nombres_lista.len()).collect::<Vec<_>>(),
            )
            .paso(
                "distancias.iter().map(|&(_, otro)| otro).collect()",
                format!(
                    "La lista de vecinos de {} es {}.",
                    nombres[punto],
                    algunos(&nombres_lista, VECINOS_A_LA_VISTA)
                ),
            ),
    );
    let vecinos = comunidades::vecinos_mas_cercanos(&filas, k);
    let resumen_listas: Vec<String> = vecinos
        .iter()
        .take(VECINOS_A_LA_VISTA)
        .enumerate()
        .map(|(gen_, lista)| {
            let cercanos: Vec<String> = lista.iter().map(|&otro| nombres[otro].clone()).collect();
            format!(
                "{}: {}",
                nombres[gen_],
                algunos(&cercanos, NOMBRES_A_LA_VISTA)
            )
        })
        .collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", k)
            .variable("vecinos", format!("{} listas de {k}", miles(n)))
            .lista("Listas de vecinos", resumen_listas, &[0])
            .paso(
                "let vecinos = vecinos_mas_cercanos(&filas, k);",
                format!(
                    "Lo mismo para cada gen, repartido entre los núcleos de tu computadora: {} listas de {k} vecinos.",
                    miles(n)
                ),
            ),
    );
    let vecindarios: Vec<Vec<usize>> = vecinos
        .iter()
        .enumerate()
        .map(|(gen_, lista)| comunidades::vecindario_con_el_mismo(lista, gen_))
        .collect();
    let mut filas_pesos: Vec<Vec<String>> = Vec::new();
    for &otro in &vecinos[punto] {
        let comunes = comunidades::compartidos(&vecindarios[punto], &vecindarios[otro]);
        let peso = comunidades::peso_jaccard(&vecindarios[punto], &vecindarios[otro]);
        let queda = peso >= PESO_MINIMO_SNN;
        filas_pesos.push(vec![
            nombres[otro].clone(),
            miles(comunes),
            decimal(peso, 3),
            if queda { "queda" } else { "se descarta" }.to_string(),
        ]);
        let ultima = filas_pesos.len() - 1;
        narracion.agregar(
            Pizarra::nueva()
                .variable("punto", &nombres[punto])
                .variable("otro", &nombres[otro])
                .variable("comunes", comunes)
                .variable("peso", decimal(peso, 3))
                .tabla(
                    &format!("Uniones de {}", nombres[punto]),
                    &["otro", "comparten", "peso", "¿queda?"],
                    filas_pesos.clone(),
                    &[(ultima, 2)],
                )
                .paso(
                    "let peso = peso_jaccard(&vecindarios[punto], &vecindarios[otro]);",
                    if queda {
                        format!(
                            "{} y {} comparten {comunes} genes en sus vecindarios de {} (cada vecindario incluye al propio gen): el peso de Jaccard, compartidos / juntos, es {comunes}/{} = {}. Pasa el mínimo de 1/15 y la unión queda.",
                            nombres[punto],
                            nombres[otro],
                            k + 1,
                            2 * (k + 1) - comunes,
                            decimal(peso, 3)
                        )
                    } else {
                        format!(
                            "{} y {} solo comparten {comunes}: peso {}, menor que 1/15. Esa unión se descarta, así el ruido no pega grupos distintos.",
                            nombres[punto],
                            nombres[otro],
                            decimal(peso, 3)
                        )
                    },
                ),
        );
    }
    let grafo = comunidades::grafo_snn(&vecinos);
    narracion.agregar(
        Pizarra::nueva()
            .variable("aristas", miles(contar_aristas(&grafo)))
            .variable("peso_total", decimal(grafo.peso_total, 3))
            .paso(
                "let grafo = grafo_snn(&vecinos);",
                format!(
                    "Con todas las parejas de vecinos, el grafo tiene {} uniones y un peso total m = {}. Dos genes quedan muy unidos si, además de ser vecinos, tienen los mismos vecinos.",
                    miles(contar_aristas(&grafo)),
                    decimal(grafo.peso_total, 3)
                ),
            ),
    );
    let comunidad = narrar_louvain(&mut narracion, &grafo, &tabla, &ajustes);
    let (_, cantidad) = comunidades::renumerar(&comunidad);
    let q = comunidades::modularidad(&grafo, &comunidad, ajustes.resolucion);
    let cuentas: Vec<String> = tamanos(&comunidad, cantidad)
        .iter()
        .map(|&cuenta| miles(cuenta))
        .collect();
    Ok(narracion.terminar(
        Pizarra::nueva()
            .variable("cantidad", cantidad)
            .variable("modularidad", decimal(q, 4))
            .lista(
                "Tamaños de las comunidades",
                cuentas.clone(),
                &[],
            )
            .paso(
                "let (comunidad, cantidad) = ordenar_por_tamano(&asignacion);",
                format!(
                    "Resultado: {} {} de {} genes, con modularidad {}. Se numeran de la más grande a la más chica. Mira «Qué significa» para leer el número.",
                    cantidad,
                    segun_cantidad(cantidad, "comunidad", "comunidades"),
                    unir(&cuentas),
                    decimal(q, 4)
                ),
            ),
    ))
}

fn narrar_louvain(
    narracion: &mut Narracion,
    grafo: &Grafo,
    tabla: &Expresion,
    ajustes: &Ajustes,
) -> Vec<usize> {
    let nombres = &tabla.nombres;
    let resolucion = ajustes.resolucion;
    let mut azar = Azar::nuevo(ajustes.semilla);
    let mut asignacion: Vec<usize> = (0..grafo.vecinos.len()).collect();
    let mut agregado: Option<Grafo> = None;
    let mut nivel = 0;
    if grafo.peso_total <= 0.0 {
        narracion.agregar(Pizarra::nueva().variable("peso_total", "0").paso(
            "if grafo.peso_total > 0.0 {",
            "El grafo no tiene uniones: no hay nada que agrupar y cada gen queda solo.",
        ));
        let (comunidad, _) = comunidades::ordenar_por_tamano(&asignacion);
        return comunidad;
    }
    loop {
        let actual = agregado.as_ref().unwrap_or(grafo);
        let mut particion = comunidades::particion_inicial(actual);
        let orden = comunidades::barajar(actual.vecinos.len(), &mut azar);
        let nombre_de = |nodo: usize| -> String {
            if nivel == 0 {
                nombres[nodo].clone()
            } else {
                format!("bloque {}", nodo + 1)
            }
        };
        if nivel == 0 {
            let primeros: Vec<String> = orden
                .iter()
                .take(VECINOS_A_LA_VISTA)
                .map(|&nodo| nombre_de(nodo))
                .collect();
            narracion.agregar(
                Pizarra::nueva()
                    .variable("nivel", nivel + 1)
                    .variable("orden", algunos(&primeros, NOMBRES_A_LA_VISTA))
                    .variable(
                        "modularidad",
                        decimal(comunidades::modularidad(grafo, &asignacion, resolucion), 4),
                    )
                    .paso(
                        "let mut particion = particion_inicial(actual);",
                        format!(
                            "Louvain empieza con cada gen en su propia comunidad (modularidad {}) y los visita en un orden al azar fijado por la semilla {}.",
                            decimal(comunidades::modularidad(grafo, &asignacion, resolucion), 4),
                            ajustes.semilla
                        ),
                    ),
            );
        }
        let mut movimientos = 0;
        for pasada in 1..=MAXIMO_PASADAS {
            let movidos = if nivel == 0 && pasada == 1 {
                narrar_primera_pasada(narracion, actual, &mut particion, &orden, tabla, resolucion)
            } else {
                comunidades::pasada_local(actual, &mut particion, &orden, resolucion)
            };
            movimientos += movidos;
            let compuesta: Vec<usize> = asignacion
                .iter()
                .map(|&nodo| particion.comunidad[nodo])
                .collect();
            let q = comunidades::modularidad(grafo, &compuesta, resolucion);
            narracion.agregar(
                Pizarra::nueva()
                    .variable("nivel", nivel + 1)
                    .variable("pasada", pasada)
                    .variable("movidos", movidos)
                    .variable("modularidad", decimal(q, 4))
                    .paso(
                        "let movidos = pasada_local(actual, &mut particion, &orden, resolucion);",
                        if movidos == 0 {
                            format!(
                                "Pasada {pasada} del nivel {}: ningún nodo mejora cambiándose de comunidad. Este nivel terminó con modularidad {}.",
                                nivel + 1,
                                decimal(q, 4)
                            )
                        } else {
                            format!(
                                "Pasada {pasada} del nivel {}: {} {} de comunidad y la modularidad queda en {}. Repetimos hasta que nadie se mueva.",
                                nivel + 1,
                                miles(movidos),
                                segun_cantidad(movidos, "nodo cambió", "nodos cambiaron"),
                                decimal(q, 4)
                            )
                        },
                    ),
            );
            if movidos == 0 {
                break;
            }
        }
        if movimientos == 0 {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("nivel", nivel + 1)
                    .variable("movimientos", 0)
                    .paso(
                        "if movimientos == 0 {",
                        if nivel == 0 {
                            "Ningún gen mejora la modularidad cambiándose de comunidad: con esta resolución juntar genes no la sube, así que cada gen queda solo y Louvain se detiene.".to_string()
                        } else {
                            format!(
                                "En el nivel {} ningún bloque quiso juntarse con otro: juntar dos comunidades ya no sube la modularidad, así que Louvain se detiene.",
                                nivel + 1
                            )
                        },
                    ),
            );
            break;
        }
        let (comunidad, cantidad) = comunidades::renumerar(&particion.comunidad);
        for nodo in asignacion.iter_mut() {
            *nodo = comunidad[*nodo];
        }
        narracion.agregar(
            Pizarra::nueva()
                .variable("nivel", nivel + 1)
                .variable("movimientos", movimientos)
                .variable("cantidad", cantidad)
                .variable(
                    "modularidad",
                    decimal(comunidades::modularidad(grafo, &asignacion, resolucion), 4),
                )
                .paso(
                    "let siguiente = agregar_comunidades(actual, &comunidad, cantidad);",
                    format!(
                        "Fin del nivel {}: {} nodos quedaron en {} comunidades. Cada comunidad se convierte en un solo nodo (un bloque); las uniones de adentro pasan a ser un lazo y las de afuera se suman. El siguiente nivel mueve bloques enteros.",
                        nivel + 1,
                        miles(actual.vecinos.len()),
                        miles(cantidad)
                    ),
                ),
        );
        agregado = Some(comunidades::agregar_comunidades(
            actual, &comunidad, cantidad,
        ));
        nivel += 1;
    }
    let (comunidad, _) = comunidades::ordenar_por_tamano(&asignacion);
    comunidad
}

fn narrar_primera_pasada(
    narracion: &mut Narracion,
    grafo: &Grafo,
    particion: &mut Particion,
    orden: &[usize],
    tabla: &Expresion,
    resolucion: f64,
) -> usize {
    let nombres = &tabla.nombres;
    let mut peso_hacia = vec![0.0; grafo.vecinos.len()];
    let mut tocadas = Vec::new();
    let mut movidos = 0;
    for (posicion, &nodo) in orden.iter().enumerate() {
        let opciones = if posicion < NODOS_DETALLADOS {
            candidatas(grafo, particion, nodo, resolucion)
        } else {
            Vec::new()
        };
        let representantes = if posicion < NODOS_DETALLADOS {
            representantes_de(&particion.comunidad)
        } else {
            Vec::new()
        };
        let movimiento = comunidades::mover_un_nodo(
            grafo,
            particion,
            nodo,
            resolucion,
            &mut peso_hacia,
            &mut tocadas,
        );
        if movimiento.hacia != movimiento.desde {
            movidos += 1;
        }
        if posicion >= NODOS_DETALLADOS {
            continue;
        }
        let comunidad_de = |c: usize| format!("la de {}", nombres[representantes[c]]);
        let filas: Vec<Vec<String>> = opciones
            .iter()
            .take(CANDIDATAS_A_LA_VISTA)
            .map(|&(c, peso, ganancia)| {
                vec![
                    if c == movimiento.desde {
                        format!("{} (la suya)", comunidad_de(c))
                    } else {
                        comunidad_de(c)
                    },
                    decimal(peso, 3),
                    format!("{:+.5}", ganancia),
                ]
            })
            .collect();
        let elegida: Vec<(usize, usize)> = opciones
            .iter()
            .take(CANDIDATAS_A_LA_VISTA)
            .position(|&(c, _, _)| c == movimiento.hacia)
            .map(|fila| vec![(fila, 0), (fila, 2)])
            .unwrap_or_default();
        let q = comunidades::modularidad(grafo, &particion.comunidad, resolucion);
        let se_movio = movimiento.hacia != movimiento.desde;
        narracion.agregar(
            Pizarra::nueva()
                .variable("nodo", &nombres[nodo])
                .variable("desde", comunidad_de(movimiento.desde))
                .variable("hacia", comunidad_de(movimiento.hacia))
                .variable("ganancia", corto(movimiento.ganancia))
                .variable("modularidad", decimal(q, 4))
                .tabla(
                    &format!("Comunidades vecinas de {}", nombres[nodo]),
                    &["comunidad", "peso hacia ella", "ΔQ"],
                    filas,
                    &elegida,
                )
                .paso(
                    if se_movio {
                        "particion.comunidad[nodo] = hacia;"
                    } else {
                        "if posible > mejor + MEJORA_MINIMA {"
                    },
                    if se_movio && posicion == 0 {
                        format!(
                            "{} prueba cada comunidad de sus vecinos y calcula cuánto subiría la modularidad si se mudara: ΔQ = peso hacia ella / m − resolución · grado de ella · grado propio / 2m², comparado con quedarse. La mejor es {} con ΔQ = {}, así que se muda; la modularidad sube a {}.",
                            nombres[nodo],
                            comunidad_de(movimiento.hacia),
                            corto(movimiento.ganancia),
                            decimal(q, 4)
                        )
                    } else if se_movio {
                        format!(
                            "{}: la mejor opción es {} (ΔQ = {}); se muda y la modularidad queda en {}.",
                            nombres[nodo],
                            comunidad_de(movimiento.hacia),
                            corto(movimiento.ganancia),
                            decimal(q, 4)
                        )
                    } else {
                        format!(
                            "Para {} ninguna comunidad vecina sube la modularidad: se queda en {}.",
                            nombres[nodo],
                            comunidad_de(movimiento.desde)
                        )
                    },
                ),
        );
    }
    if orden.len() > NODOS_DETALLADOS {
        narracion.agregar(
            Pizarra::nueva()
                .variable("nodo", "…")
                .variable("movidos", movidos)
                .paso(
                    "for &nodo in orden {",
                    format!(
                        "Los {} genes restantes de esta pasada siguen la misma regla sin mostrar cada uno.",
                        miles(orden.len() - NODOS_DETALLADOS)
                    ),
                ),
        );
    }
    movidos
}

fn representantes_de(comunidad: &[usize]) -> Vec<usize> {
    let mut representantes = vec![usize::MAX; comunidad.len()];
    for (gen_, &c) in comunidad.iter().enumerate() {
        if representantes[c] == usize::MAX {
            representantes[c] = gen_;
        }
    }
    representantes
}

fn nubosidad(coordenadas: &[[f64; 2]]) -> f64 {
    coordenadas
        .iter()
        .map(|punto| punto[0].abs().max(punto[1].abs()))
        .fold(0.0, f64::max)
}

fn pasos_tsne(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let tabla = leer_expresion(entrada, MAXIMO_FILAS_TSNE)?;
    let ajustes = leer_ajustes(entrada, &tabla)?;
    let opciones = leer_opciones_tsne(entrada, &tabla, ajustes.semilla)?;
    let Traza {
        narracion,
        coordenadas,
        historia,
        divergencia_inicial,
    } = trazar_tsne(&tabla, &ajustes, opciones)?;
    let divergencia = historia.last().map_or(0.0, |&(_, valor)| valor);
    Ok(narracion.terminar(
        Pizarra::nueva()
            .variable("iteracion", opciones.iteraciones)
            .variable("divergencia", decimal(divergencia, 4))
            .variable("coordenadas", format!("radio {}", corto(nubosidad(&coordenadas))))
            .paso(
                "let divergencia = historia.last().map_or(0.0, |&(_, valor)| valor);",
                format!(
                    "Listo: tras {} iteraciones la KL pasó de {} a {}. Cada gen tiene su punto en el plano; mira «Qué significa» para saber qué leer en el mapa y qué no.",
                    miles(opciones.iteraciones),
                    decimal(divergencia_inicial, 4),
                    decimal(divergencia, 4)
                ),
            ),
    ))
}

struct Traza {
    narracion: Narracion,
    coordenadas: Vec<[f64; 2]>,
    historia: Vec<(usize, f64)>,
    divergencia_inicial: f64,
}

fn trazar_tsne(
    tabla: &Expresion,
    ajustes: &Ajustes,
    opciones: OpcionesTsne,
) -> Result<Traza, String> {
    let n = tabla.cantidad();
    let nombres = &tabla.nombres;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        Pizarra::nueva()
            .variable("puntos", format!("{} genes × {} columnas", miles(n), tabla.dimension()))
            .variable("perplejidad", decimal(opciones.perplejidad, 1))
            .variable("exageracion", decimal(opciones.exageracion, 1))
            .variable("iteraciones", opciones.iteraciones)
            .paso(
                "revisar_tsne(puntos, &opciones)?;",
                format!(
                    "Tienes {} genes con {} mediciones. La meta: un punto en el plano por gen, de modo que los genes vecinos en la tabla queden vecinos en el dibujo.",
                    miles(n),
                    tabla.dimension()
                ),
            ),
    );
    let filas = comunidades::filas_para_medir(&tabla.filas, ajustes.medida);
    narracion.agregar(
        Pizarra::nueva()
            .variable("medida", nombre_de_medida(ajustes.medida))
            .paso(
                "let filas = filas_para_medir(puntos, medida);",
                match ajustes.medida {
                    Medida::Pearson => "Con Pearson, cada gen se centra y se divide por su largo: la distancia al cuadrado entre dos genes queda 2 · (1 − r).".to_string(),
                    Medida::Euclidea => "Con distancia euclídea los genes se usan tal cual.".to_string(),
                },
            ),
    );
    let propio = 0;
    let distancias: Vec<f64> = filas
        .iter()
        .map(|otra| agrupamiento::distancia_al_cuadrado(&filas[propio], otra))
        .collect();
    let mut probabilidades = vec![0.0; n];
    let objetivo = opciones.perplejidad.ln();
    let mut beta = comunidades::beta_inicial(&distancias, propio);
    let mut beta_minima = 0.0;
    let mut beta_maxima = f64::INFINITY;
    let mut entropia =
        comunidades::entropia_con_beta(&distancias, propio, beta, &mut probabilidades);
    let mut intentos = 1;
    let pizarra_beta = |beta: f64, entropia: f64, beta_minima: f64, beta_maxima: f64| {
        Pizarra::nueva()
            .variable("propio", &nombres[propio])
            .variable("objetivo", decimal(objetivo, 4))
            .variable("beta", corto(beta))
            .variable("entropia", decimal(entropia, 4))
            .variable("beta_minima", corto(beta_minima))
            .variable(
                "beta_maxima",
                if beta_maxima.is_finite() {
                    corto(beta_maxima)
                } else {
                    "∞".to_string()
                },
            )
    };
    narracion.agregar(pizarra_beta(beta, entropia, beta_minima, beta_maxima).paso(
        "let mut entropia = entropia_con_beta(distancias, propio, beta, probabilidades);",
        format!(
            "Para {} convertimos distancias en probabilidades con una campana: p = e^(−β · d²). Queremos que la entropía valga ln({}) = {}, o sea, que {} mire unos {} vecinos efectivos. Arrancamos con β = {} (uno entre la distancia² media a los demás) y la entropía es {}.",
            nombres[propio],
            decimal(opciones.perplejidad, 1),
            decimal(objetivo, 4),
            nombres[propio],
            decimal(opciones.perplejidad, 1),
            corto(beta),
            decimal(entropia, 4)
        ),
    ));
    while (entropia - objetivo).abs() > TOLERANCIA_DE_ENTROPIA && intentos < INTENTOS_DE_BETA {
        let mas_ancha = entropia > objetivo;
        if mas_ancha {
            beta_minima = beta;
            beta = if beta_maxima.is_finite() {
                (beta + beta_maxima) / 2.0
            } else {
                beta * 2.0
            };
        } else {
            beta_maxima = beta;
            beta = (beta + beta_minima) / 2.0;
        }
        entropia = comunidades::entropia_con_beta(&distancias, propio, beta, &mut probabilidades);
        intentos += 1;
        if intentos <= INTENTOS_DETALLADOS {
            narracion.agregar(pizarra_beta(beta, entropia, beta_minima, beta_maxima).paso(
                if mas_ancha {
                    "beta_minima = beta;"
                } else {
                    "beta_maxima = beta;"
                },
                if mas_ancha {
                    format!(
                        "La entropía era mayor que la meta: la campana es demasiado ancha. Subimos β a {} (σ = {}) y la entropía queda en {}.",
                        corto(beta),
                        decimal((1.0 / (2.0 * beta)).sqrt(), 4),
                        decimal(entropia, 4)
                    )
                } else {
                    format!(
                        "La entropía era menor que la meta: la campana es demasiado angosta. Bajamos β a {} (σ = {}) y la entropía queda en {}.",
                        corto(beta),
                        decimal((1.0 / (2.0 * beta)).sqrt(), 4),
                        decimal(entropia, 4)
                    )
                },
            ));
        }
    }
    let mut cercanos: Vec<usize> = (0..n).filter(|&otro| otro != propio).collect();
    cercanos.sort_by(|&a, &b| {
        probabilidades[b]
            .total_cmp(&probabilidades[a])
            .then(a.cmp(&b))
    });
    let filas_p: Vec<Vec<String>> = cercanos
        .iter()
        .take(VECINOS_A_LA_VISTA)
        .map(|&otro| vec![nombres[otro].clone(), decimal(probabilidades[otro], 4)])
        .collect();
    narracion.agregar(
        pizarra_beta(beta, entropia, beta_minima, beta_maxima)
            .tabla(
                &format!("Probabilidad de que {} elija a cada vecino", nombres[propio]),
                &["otro", "p"],
                filas_p,
                &[(0, 1)],
            )
            .paso(
                "while (entropia - objetivo).abs() > TOLERANCIA_DE_ENTROPIA && intentos < INTENTOS_DE_BETA {",
                format!(
                    "Tras {intentos} intentos de bisección, β = {} (σ = {}) da perplejidad e^{} = {}. Arriba, los vecinos que más probabilidad reciben.",
                    corto(beta),
                    decimal((1.0 / (2.0 * beta)).sqrt(), 4),
                    decimal(entropia, 4),
                    decimal(entropia.exp(), 3)
                ),
            ),
    );
    let (conjuntas, sigmas) = comunidades::afinidades(&filas, opciones.perplejidad);
    let mut ordenados = sigmas.clone();
    ordenados.sort_by(f64::total_cmp);
    let mejor = cercanos[0];
    narracion.agregar(
        Pizarra::nueva()
            .variable("sigmas", format!("{} · {} · {}", decimal(ordenados[0], 3), decimal(ordenados[n / 2], 3), decimal(ordenados[n - 1], 3)))
            .variable("conjuntas", format!("{n} × {n}"))
            .paso(
                "let conjunta = (condicionales[i * n + j] + condicionales[j * n + i]) / total;",
                format!(
                    "Cada gen calibra su propio σ (mínimo · mediana · máximo arriba): en zonas densas la campana es angosta y en zonas ralas es ancha. Luego se promedia p(j|i) con p(i|j) y se divide por 2n; por ejemplo, la afinidad entre {} y {} es {}.",
                    nombres[propio],
                    nombres[mejor],
                    corto(conjuntas[propio * n + mejor])
                ),
            ),
    );
    let tasa = comunidades::tasa_automatica(n, opciones.exageracion);
    let mut coordenadas = comunidades::coordenadas_iniciales(n, opciones.semilla);
    let mut velocidad = vec![[0.0, 0.0]; n];
    let mut ganancias = vec![[1.0, 1.0]; n];
    let mut historia = Vec::new();
    let divergencia_inicial = comunidades::divergencia_kl(&conjuntas, &coordenadas);
    let pizarra_iteracion = |iteracion: usize,
                             exageracion: f64,
                             momento: f64,
                             divergencia: f64,
                             coordenadas: &[[f64; 2]]| {
        Pizarra::nueva()
            .variable("iteracion", iteracion)
            .variable("exageracion", decimal(exageracion, 1))
            .variable("momento", decimal(momento, 1))
            .variable("tasa", decimal(tasa, 1))
            .variable("divergencia", decimal(divergencia, 4))
            .variable(
                "coordenadas",
                format!("radio {}", corto(nubosidad(coordenadas))),
            )
    };
    narracion.agregar(
        pizarra_iteracion(0, opciones.exageracion, MOMENTO_INICIAL, divergencia_inicial, &coordenadas).paso(
            "let mut coordenadas = coordenadas_iniciales(n, opciones.semilla);",
            format!(
                "Los {} puntos empiezan en una nube diminuta al azar (semilla {}). La divergencia KL entre las afinidades de la tabla y las del plano es {}: muy alta, porque en el plano todos están igual de cerca.",
                miles(n),
                opciones.semilla,
                decimal(divergencia_inicial, 4)
            ),
        ),
    );
    for iteracion in 1..=opciones.iteraciones {
        let exagerada = iteracion <= ITERACIONES_EXAGERADAS;
        let exageracion = if exagerada { opciones.exageracion } else { 1.0 };
        let momento = if exagerada {
            MOMENTO_INICIAL
        } else {
            MOMENTO_FINAL
        };
        if iteracion == ITERACIONES_EXAGERADAS + 1 {
            velocidad.fill([0.0, 0.0]);
            ganancias.fill([1.0, 1.0]);
            narracion.agregar(
                pizarra_iteracion(iteracion, exageracion, momento, historia.last().map_or(divergencia_inicial, |&(_, valor)| valor), &coordenadas).paso(
                    "velocidad.fill([0.0, 0.0]);",
                    format!(
                        "Terminan las {ITERACIONES_EXAGERADAS} iteraciones de exageración: las atracciones vuelven a su tamaño real, el momento sube a {} y la velocidad y las ganancias vuelven a empezar.",
                        decimal(MOMENTO_FINAL, 1)
                    ),
                ),
            );
        }
        let empuje = comunidades::gradiente(&conjuntas, &coordenadas, exageracion);
        comunidades::descender(
            &mut coordenadas,
            &mut velocidad,
            &mut ganancias,
            &empuje,
            momento,
            tasa,
        );
        let medir = iteracion % CADA_CUANTO_MEDIR == 0 || iteracion == opciones.iteraciones;
        if medir {
            historia.push((
                iteracion,
                comunidades::divergencia_kl(&conjuntas, &coordenadas),
            ));
        }
        if iteracion <= ITERACIONES_DETALLADAS {
            let divergencia = comunidades::divergencia_kl(&conjuntas, &coordenadas);
            let fuerza = empuje.iter().map(|g| g[0].hypot(g[1])).fold(0.0, f64::max);
            narracion.agregar(
                pizarra_iteracion(iteracion, exageracion, momento, divergencia, &coordenadas).paso(
                    "let empuje = gradiente(&conjuntas, &coordenadas, exageracion);",
                    if iteracion == 1 {
                        format!(
                            "Iteración 1: cada punto recibe un tirón hacia sus vecinos de la tabla (multiplicado por la exageración {}) y un empujón lejos de todos los demás; el mayor empuje mide {}. Luego cada punto avanza con tasa {} y momento {}; la KL queda en {}.",
                            decimal(exageracion, 1),
                            corto(fuerza),
                            decimal(tasa, 1),
                            decimal(momento, 1),
                            decimal(divergencia, 4)
                        )
                    } else {
                        format!(
                            "Iteración {iteracion}: mayor empuje {}, la nube mide {} de radio y la KL es {}.",
                            corto(fuerza),
                            corto(nubosidad(&coordenadas)),
                            decimal(divergencia, 4)
                        )
                    },
                ),
            );
        } else if medir {
            let divergencia = historia.last().map_or(0.0, |&(_, valor)| valor);
            let anterior = historia
                .len()
                .checked_sub(2)
                .map(|posicion| historia[posicion].1);
            let tendencia = match anterior {
                Some(antes) if divergencia < antes => "bajó",
                Some(antes) if divergencia > antes => "subió",
                _ => "no cambió",
            };
            narracion.agregar(
                pizarra_iteracion(iteracion, exageracion, momento, divergencia, &coordenadas).paso(
                    "historia.push((iteracion, divergencia_kl(&conjuntas, &coordenadas)));",
                    if exagerada && anterior.is_none() {
                        format!(
                            "Iteración {iteracion}: la KL vale {}. Mientras dura la exageración, los vecinos se atraen con fuerza de más: los grupos se juntan rápido en islas apretadas y la KL real sube y baja sin llegar a su mínimo.",
                            decimal(divergencia, 4)
                        )
                    } else if exagerada {
                        format!(
                            "Iteración {iteracion}: la KL {tendencia} a {}; seguimos con exageración.",
                            decimal(divergencia, 4)
                        )
                    } else {
                        format!(
                            "Iteración {iteracion}: la KL {tendencia} a {}; la nube mide {} de radio.",
                            decimal(divergencia, 4),
                            corto(nubosidad(&coordenadas))
                        )
                    },
                ),
            );
        }
    }
    Ok(Traza {
        narracion,
        coordenadas,
        historia,
        divergencia_inicial,
    })
}

#[cfg(test)]
mod pruebas {
    use super::*;
    use crate::ejemplo;

    fn entrada_de(id: &str, cambios: &[(&str, &str)]) -> Entrada {
        let mut valores = ejemplo(id).unwrap();
        for (campo, valor) in cambios {
            valores[*campo] = serde_json::Value::String(valor.to_string());
        }
        Entrada::desde_json(valores).unwrap()
    }

    fn louvain_del_ejemplo(entrada: &Entrada) -> Comunidades {
        let tabla = leer_expresion(entrada, MAXIMO_FILAS).unwrap();
        let ajustes = leer_ajustes(entrada, &tabla).unwrap();
        comunidades::agrupar_por_comunidades(
            &tabla.filas,
            ajustes.medida,
            ajustes.k,
            ajustes.resolucion,
            ajustes.semilla,
        )
        .unwrap()
    }

    #[test]
    fn el_ejemplo_tiene_cuatro_patrones_claros() {
        let entrada = entrada_de("agrupamiento.louvain", &[]);
        let resultado = louvain_del_ejemplo(&entrada);
        assert_eq!(resultado.louvain.cantidad, 4);
        assert!(resultado.louvain.modularidad > MODULARIDAD_CLARA);
        let tabla = leer_expresion(&entrada, MAXIMO_FILAS).unwrap();
        assert_eq!(tabla.cantidad(), 60);
        assert_eq!(tabla.dimension(), 8);
        let salida = louvain(&entrada).unwrap();
        assert!(salida.interpretacion.len() >= 4);
    }

    #[test]
    fn el_paso_a_paso_de_louvain_llega_al_mismo_resultado() {
        for cambios in [
            vec![],
            vec![("distancia", MEDIDA_EUCLIDEA), ("resolucion", "2.5")],
            vec![("k", "3"), ("semilla", "7")],
        ] {
            let entrada = entrada_de("agrupamiento.louvain", &cambios);
            let resultado = louvain_del_ejemplo(&entrada);
            let tabla = leer_expresion(&entrada, MAXIMO_FILAS).unwrap();
            let ajustes = leer_ajustes(&entrada, &tabla).unwrap();
            let filas = comunidades::filas_para_medir(&tabla.filas, ajustes.medida);
            let grafo =
                comunidades::grafo_snn(&comunidades::vecinos_mas_cercanos(&filas, ajustes.k));
            let mut narracion = Narracion::nueva(MAXIMO_PASOS);
            let comunidad = narrar_louvain(&mut narracion, &grafo, &tabla, &ajustes);
            assert_eq!(comunidad, resultado.louvain.comunidad);
            let pasos = pasos_louvain(&entrada).unwrap();
            assert!(pasos.len() <= MAXIMO_PASOS && pasos.len() > 20);
            let ultimo = pasos.last().unwrap();
            assert!(
                ultimo
                    .nota
                    .contains(&decimal(resultado.louvain.modularidad, 4))
            );
        }
    }

    #[test]
    fn la_narracion_de_un_movimiento_coincide_con_la_ganancia() {
        let entrada = entrada_de("agrupamiento.louvain", &[]);
        let resultado = louvain_del_ejemplo(&entrada);
        let grafo = &resultado.grafo;
        let mut particion = comunidades::particion_inicial(grafo);
        let mut peso_hacia = vec![0.0; grafo.vecinos.len()];
        let mut tocadas = Vec::new();
        for nodo in 0..grafo.vecinos.len() {
            let opciones = candidatas(grafo, &particion, nodo, 1.0);
            let movimiento = comunidades::mover_un_nodo(
                grafo,
                &mut particion,
                nodo,
                1.0,
                &mut peso_hacia,
                &mut tocadas,
            );
            let mejor = opciones[0];
            assert!((mejor.2 - movimiento.ganancia).abs() < 1e-12);
            if movimiento.ganancia > 0.0 {
                assert_eq!(mejor.0, movimiento.hacia);
            }
        }
    }

    #[test]
    fn el_paso_a_paso_de_tsne_llega_al_mismo_mapa() {
        let entrada = entrada_de("agrupamiento.tsne", &[("iteraciones", "400")]);
        let tabla = leer_expresion(&entrada, MAXIMO_FILAS_TSNE).unwrap();
        let ajustes = leer_ajustes(&entrada, &tabla).unwrap();
        let opciones = leer_opciones_tsne(&entrada, &tabla, ajustes.semilla).unwrap();
        let directo = comunidades::tsne(&tabla.filas, ajustes.medida, opciones).unwrap();
        let traza = trazar_tsne(&tabla, &ajustes, opciones).unwrap();
        assert_eq!(traza.coordenadas, directo.coordenadas);
        assert_eq!(traza.historia, directo.historia);
        let pasos = pasos_tsne(&entrada).unwrap();
        assert!(pasos.len() <= MAXIMO_PASOS);
        let mut probabilidades = vec![0.0; tabla.cantidad()];
        let filas = comunidades::filas_para_medir(&tabla.filas, ajustes.medida);
        let distancias: Vec<f64> = filas
            .iter()
            .map(|otra| agrupamiento::distancia_al_cuadrado(&filas[0], otra))
            .collect();
        let calibracion = comunidades::calibrar_un_punto(
            &distancias,
            0,
            opciones.perplejidad,
            &mut probabilidades,
        );
        assert!(pasos.iter().any(|paso| {
            paso.nota
                .contains(&format!("Tras {} intentos", calibracion.intentos))
        }));
    }

    #[test]
    fn el_ejemplo_de_tsne_separa_las_comunidades() {
        let entrada = entrada_de("agrupamiento.tsne", &[]);
        let tabla = leer_expresion(&entrada, MAXIMO_FILAS_TSNE).unwrap();
        let ajustes = leer_ajustes(&entrada, &tabla).unwrap();
        let opciones = leer_opciones_tsne(&entrada, &tabla, ajustes.semilla).unwrap();
        let mapa = calcular_mapa(&tabla, &ajustes, opciones).unwrap();
        assert!(mapa.islas >= ISLAS_CLARAS);
        assert!(mapa.conservacion > 0.5);
        assert!(mapa.tsne.divergencia < 1.0);
        let pasos = pasos_tsne(&entrada).unwrap();
        assert!(pasos.len() <= MAXIMO_PASOS && pasos.len() > 20);
    }

    #[test]
    fn las_entradas_malas_dan_errores_en_espanol() {
        for (campo, valor) in [
            ("k", "0"),
            ("k", "60"),
            ("resolucion", "0"),
            ("perplejidad", "30"),
            ("iteraciones", "0"),
            ("exageracion", "0.5"),
            ("datos", "1 2\n3"),
            ("datos", "a 1\nb 2\nc x"),
        ] {
            let entrada = entrada_de("agrupamiento.tsne", &[(campo, valor)]);
            assert!(tsne(&entrada).is_err(), "{campo} = {valor}");
            assert!(pasos_tsne(&entrada).is_err(), "{campo} = {valor}");
        }
        let entrada = entrada_de("agrupamiento.louvain", &[("datos", "a 1\nb 2")]);
        assert!(louvain(&entrada).is_err());
        let solo_numeros = entrada_de("agrupamiento.louvain", &[("datos", "1 2\n2 3\n3 4\n9 9")]);
        let salida = louvain(&entrada_de(
            "agrupamiento.louvain",
            &[("datos", "1 2\n2 3\n3 4\n9 9"), ("k", "2")],
        ))
        .unwrap();
        assert!(!salida.bloques.is_empty());
        assert!(louvain(&solo_numeros).is_err());
    }

    #[test]
    fn el_ruido_no_se_presenta_como_estructura_clara() {
        let mut azar = Azar::nuevo(5);
        let texto: Vec<String> = (1..=80)
            .map(|i| {
                let valores: Vec<String> =
                    (0..8).map(|_| format!("{:.3}", azar.normal())).collect();
                format!("g{i} {}", valores.join(" "))
            })
            .collect();
        let entrada = entrada_de("agrupamiento.louvain", &[("datos", &texto.join("\n"))]);
        let salida = louvain(&entrada).unwrap();
        let interpretacion = salida.interpretacion.join(" ");
        assert!(!interpretacion.contains("estructura de comunidades clara"));
        assert!(interpretacion.contains("revueltos"));
        let ejemplo = louvain(&entrada_de("agrupamiento.louvain", &[])).unwrap();
        assert!(
            ejemplo
                .interpretacion
                .join(" ")
                .contains("estructura de comunidades clara")
        );
    }

    #[test]
    fn cada_comunidad_se_nombra_con_uno_de_sus_miembros() {
        let comunidad = vec![3, 3, 0, 5, 0, 3];
        let representantes = representantes_de(&comunidad);
        for c in [0, 3, 5] {
            assert_eq!(comunidad[representantes[c]], c);
        }
    }

    #[test]
    fn la_caida_final_dice_desde_donde_se_mide() {
        let historia = [(225, 2.0), (250, 1.0), (275, 0.75), (300, 0.5)];
        let (caida, desde) = caida_final(&historia, 300).unwrap();
        assert_eq!(desde, 250);
        assert!((caida - 0.5).abs() < 1e-12);
    }

    #[test]
    fn con_resolucion_enorme_la_narracion_dice_que_nadie_se_mueve() {
        let entrada = entrada_de("agrupamiento.louvain", &[]);
        let tabla = leer_expresion(&entrada, MAXIMO_FILAS).unwrap();
        let mut ajustes = leer_ajustes(&entrada, &tabla).unwrap();
        ajustes.resolucion = 1000.0;
        let filas = comunidades::filas_para_medir(&tabla.filas, ajustes.medida);
        let grafo = comunidades::grafo_snn(&comunidades::vecinos_mas_cercanos(&filas, ajustes.k));
        let mut narracion = Narracion::nueva(MAXIMO_PASOS);
        let comunidad = narrar_louvain(&mut narracion, &grafo, &tabla, &ajustes);
        let (_, cantidad) = comunidades::renumerar(&comunidad);
        assert_eq!(cantidad, 60);
        let pasos = narracion.terminar(Pizarra::nueva().paso("fin", "fin"));
        assert!(
            pasos
                .iter()
                .any(|paso| paso.nota.contains("Ningún gen mejora"))
        );
        assert!(
            pasos
                .iter()
                .all(|paso| !paso.nota.contains("ningún bloque"))
        );
    }
}
