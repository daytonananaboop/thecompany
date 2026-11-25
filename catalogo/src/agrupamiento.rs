use crate::formato::{cientifico, decimal, miles, porcentaje};
use crate::pizarra::{Narracion, Pizarra};
use crate::{Algoritmo, Bloque, Campo, Entrada, Marca, Paso, Punto, Salida, Tipo};
use nucleo::agrupamiento::{self, Componentes, Inicio, KMediasSuaves, Medida, Union};
use nucleo::azar::Azar;
use serde_json::json;

const MAXIMO_PUNTOS: usize = 5000;
const MAXIMO_COLUMNAS: usize = 200;
const MAXIMO_K: usize = 50;
const MAXIMO_ITERACIONES: usize = 1000;
const MAXIMO_MATRIZ: usize = 2000;
const MAXIMO_K_DEL_CODO: usize = 15;
const MAXIMO_CUENTAS: f64 = 5e9;
const MAXIMO_VALOR: f64 = 1e12;
const MAXIMO_ANIMABLE: usize = 1000;
const MAXIMO_PASOS: usize = 400;
const MAXIMO_FILAS_DE_CALOR: usize = 300;
const MAXIMO_PUNTOS_DIBUJADOS: usize = 2000;
const MAXIMO_CALOR_DE_DISTANCIAS: usize = 100;
const MAXIMO_COLUMNAS_DE_GRUPOS: usize = 10;
const MAXIMO_GRUPOS_DEL_CORTE: usize = 10;
const NOMBRES_A_LA_VISTA: usize = 5;
const NOMBRES_POR_GRUPO: usize = 12;
const CORRELACION_FUERTE: f64 = 0.8;
const RESPONSABILIDAD_CLARA: f64 = 0.6;
const RESPONSABILIDAD_DURA: f64 = 0.99;
const MAYORIA: f64 = 0.9;
const MEJORA_NOTABLE: f64 = 0.05;
const CAIDA_CLARA: f64 = 0.5;
const VARIANZA_FIEL: f64 = 0.7;
const VARIANZA_POBRE: f64 = 0.5;
const MINIMO_TIEMPOS: usize = 4;
const TOLERANCIA: f64 = 1e-9;

const PRIMEROS_PUNTOS: &str = "los primeros k puntos";
const K_MEDIAS_MAS_MAS: &str = "k-means++";
const INICIOS: &[&str] = &[PRIMEROS_PUNTOS, K_MEDIAS_MAS_MAS];
const MEDIDA_PEARSON: &str = "Pearson (1 − r)";
const MEDIDA_EUCLIDEA: &str = "euclídea";
const MEDIDA_MATRIZ: &str = "la tabla ya es una matriz de distancias";
const MEDIDAS: &[&str] = &[MEDIDA_PEARSON, MEDIDA_EUCLIDEA, MEDIDA_MATRIZ];

const LINEA_PRIMEROS: &str = "Inicio::PrimerosPuntos => (0..k.min(puntos.len())).collect(),";
const LINEA_MAS_MAS: &str =
    "Inicio::KMediasMasMas(semilla) => k_medias_mas_mas(puntos, k, semilla),";
const LINEA_PRIMERA_ASIGNACION: &str = "let mut grupos = de_centros_a_grupos(puntos, &centros);";
const LINEA_MOVER: &str = "centros = de_grupos_a_centros(puntos, &grupos, &centros);";
const LINEA_ASIGNAR: &str = "let grupos_nuevos = de_centros_a_grupos(puntos, &centros);";
const LINEA_CONVERGIO: &str = "convergio = true;";
const LINEA_MAXIMO: &str = "while iteraciones < maximo_iteraciones {";

const ETIQUETA_TABLA: &str = "Tabla de expresión: un gen por línea, su nombre y sus números separados por espacios, tabulaciones o comas (punto decimal)";
const ETIQUETA_PUNTOS: &str = "Puntos: uno por línea, con sus coordenadas separadas por espacios (nombre opcional al principio)";

const TABLA_EJEMPLO: &str = "gen\t0h\t2h\t4h\t6h\t8h\t10h\t12h
gen01\t-0.17\t0.56\t0.99\t1.70\t2.17\t3.19\t3.54
gen02\t-0.43\t0.07\t1.26\t1.44\t1.87\t3.00\t2.52
gen03\t-0.04\t0.70\t0.91\t1.32\t2.43\t2.87\t2.69
gen04\t-0.45\t-0.38\t-0.54\t-1.16\t-1.05\t-2.02\t-2.32
gen05\t0.11\t-0.80\t-1.44\t-1.63\t-2.11\t-3.07\t-3.42
gen06\t0.10\t0.18\t1.38\t2.83\t1.64\t-0.12\t0.25
gen07\t0.17\t-0.85\t-0.70\t-1.05\t-1.26\t-1.55\t-1.95
gen08\t0.42\t0.21\t0.48\t1.05\t1.71\t2.12\t2.95
gen09\t-0.23\t-0.59\t-0.92\t-1.62\t-2.15\t-3.13\t-3.69
gen10\t0.17\t-0.01\t0.98\t1.42\t0.54\t0.11\t0.00
gen11\t0.42\t0.16\t0.75\t1.92\t0.81\t-0.05\t-0.28
gen12\t0.20\t0.10\t0.45\t1.09\t1.38\t1.40\t1.74
gen13\t-0.58\t-0.49\t-0.55\t-0.83\t-1.11\t-1.58\t-1.72
gen14\t0.10\t0.47\t1.67\t3.11\t1.55\t0.12\t-0.38
gen15\t-0.14\t0.39\t0.79\t1.50\t0.73\t0.04\t0.44
gen16\t0.13\t-0.55\t-1.24\t-1.39\t-2.36\t-2.73\t-3.43
gen17\t0.31\t0.32\t1.58\t2.44\t1.46\t0.57\t0.04
gen18\t0.23\t0.56\t1.45\t2.10\t1.84\t0.18\t-0.20
gen19\t0.14\t-0.29\t-1.00\t-0.72\t-2.06\t-1.77\t-2.02
gen20\t-0.07\t0.16\t1.80\t2.13\t2.49\t2.52\t3.45
gen21\t-0.19\t0.25\t1.50\t3.23\t1.80\t-0.01\t0.21
gen22\t0.06\t0.97\t0.57\t1.63\t2.24\t2.97\t3.25
gen23\t-0.11\t0.50\t1.96\t3.25\t1.61\t0.44\t-0.07
gen24\t-0.15\t-0.55\t-0.86\t-0.80\t-1.32\t-1.63\t-1.97
gen25\t-0.12\t-0.03\t1.10\t1.37\t1.84\t1.70\t1.92
gen26\t-0.01\t0.83\t1.30\t2.59\t1.02\t0.45\t0.12
gen27\t-0.18\t-0.30\t-0.52\t-0.35\t-0.78\t-1.10\t-1.36
gen28\t0.47\t-0.12\t-1.26\t-1.53\t-1.88\t-2.16\t-3.06
gen29\t0.06\t0.14\t-0.03\t0.71\t0.57\t1.08\t1.70
gen30\t-0.02\t0.36\t0.92\t1.66\t2.21\t3.04\t3.32";

const PERFILES_EJEMPLO: &str = "gen\t0h\t3h\t6h\t9h\t12h
gen01\t0.1\t0.9\t2.1\t2.9\t4.2
gen02\t0.0\t0.3\t0.6\t0.9\t1.2
gen03\t0.2\t0.1\t0.3\t0.2\t0.1
gen04\t3.9\t3.1\t2.0\t1.1\t0.1
gen05\t0.1\t1.8\t3.0\t1.7\t0.2
gen06\t1.3\t0.9\t0.6\t0.4\t0.0";

const NUBES_EJEMPLO: &str = "2.8 1.4
6.3 8.6
5.4 8.3
6.0 8.4
5.3 8.3
3.7 9.1
8.3 4.5
5.1 8.0
8.5 2.4
8.6 1.9
4.8 9.5
2.6 1.4
2.4 1.3
7.8 3.3
8.1 3.3
2.0 1.6
5.3 8.6
8.5 2.4
1.1 0.8
1.1 1.7
1.0 1.7
7.2 3.2
3.2 2.1
7.4 1.7";

const NUBES_CON_ATIPICO: &str = "2.8 1.4
6.3 8.6
5.4 8.3
6.0 8.4
5.3 8.3
3.7 9.1
8.3 4.5
5.1 8.0
8.5 2.4
8.6 1.9
4.8 9.5
2.6 1.4
14.5 13.0
2.4 1.3
7.8 3.3
8.1 3.3
2.0 1.6
5.3 8.6
8.5 2.4
1.1 0.8
1.1 1.7
1.0 1.7
7.2 3.2
3.2 2.1
7.4 1.7";

const CENTROS_EJEMPLO: &str = "2 2
8 3
5 8.5";

const DATOS_EXPRESION: Campo = Campo {
    id: "datos",
    etiqueta: ETIQUETA_TABLA,
    tipo: Tipo::Texto,
    ejemplo: TABLA_EJEMPLO,
};

const K_TRES: Campo = Campo {
    id: "k",
    etiqueta: "k (cantidad de grupos)",
    tipo: Tipo::Entero,
    ejemplo: "3",
};

const SEMILLA: Campo = Campo {
    id: "semilla",
    etiqueta: "Semilla del azar (el mismo número repite el mismo resultado)",
    tipo: Tipo::Entero,
    ejemplo: "1",
};

const ITERACIONES: Campo = Campo {
    id: "iteraciones",
    etiqueta: "Iteraciones (máximo 1000)",
    tipo: Tipo::Entero,
    ejemplo: "100",
};

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![
        Algoritmo {
            id: "agrupamiento.distancias",
            capitulo: "agrupamiento",
            nombre: "Distancia euclídea y de Pearson",
            resumen: "Qué tan distintos son los perfiles de expresión de cada par de genes: la distancia euclídea mira cuánto cambian, la de Pearson solo la forma del cambio.",
            campos: &[Campo {
                id: "datos",
                etiqueta: ETIQUETA_TABLA,
                tipo: Tipo::Texto,
                ejemplo: PERFILES_EJEMPLO,
            }],
            funcion: "agrupamiento::matriz_de_distancias",
            ejecutar: distancias,
            pasos: Some(pasos_distancias),
        },
        Algoritmo {
            id: "agrupamiento.mas_lejano_primero",
            capitulo: "agrupamiento",
            nombre: "k-centros: el más lejano primero",
            resumen: "Elige k centros entre los puntos: empieza por el primero y agrega siempre el punto más alejado de los centros ya elegidos.",
            campos: &[
                Campo {
                    id: "datos",
                    etiqueta: ETIQUETA_PUNTOS,
                    tipo: Tipo::Texto,
                    ejemplo: NUBES_CON_ATIPICO,
                },
                K_TRES,
            ],
            funcion: "agrupamiento::mas_lejano_primero",
            ejecutar: mas_lejano_primero,
            pasos: Some(pasos_mas_lejano_primero),
        },
        Algoritmo {
            id: "agrupamiento.distorsion",
            capitulo: "agrupamiento",
            nombre: "Distorsión de error cuadrático",
            resumen: "El promedio de la distancia al cuadrado de cada punto a su centro más cercano: el número que k-means intenta achicar.",
            campos: &[
                Campo {
                    id: "datos",
                    etiqueta: ETIQUETA_PUNTOS,
                    tipo: Tipo::Texto,
                    ejemplo: NUBES_EJEMPLO,
                },
                Campo {
                    id: "centros",
                    etiqueta: "Centros: uno por línea, con tantas coordenadas como los puntos",
                    tipo: Tipo::Texto,
                    ejemplo: CENTROS_EJEMPLO,
                },
            ],
            funcion: "agrupamiento::distorsion_cuadratica",
            ejecutar: distorsion,
            pasos: Some(pasos_distorsion),
        },
        Algoritmo {
            id: "agrupamiento.lloyd",
            capitulo: "agrupamiento",
            nombre: "k-means (algoritmo de Lloyd)",
            resumen: "Repite dos pasos hasta que nada cambia: cada punto va al centro más cercano y cada centro se mueve al centro de gravedad de sus puntos.",
            campos: &[
                DATOS_EXPRESION,
                K_TRES,
                Campo {
                    id: "inicio",
                    etiqueta: "Centros iniciales",
                    tipo: Tipo::Opcion(INICIOS),
                    ejemplo: PRIMEROS_PUNTOS,
                },
                SEMILLA,
                Campo {
                    id: "iteraciones",
                    etiqueta: "Iteraciones máximas (hasta 1000)",
                    tipo: Tipo::Entero,
                    ejemplo: "100",
                },
            ],
            funcion: "agrupamiento::k_medias",
            ejecutar: lloyd,
            pasos: Some(pasos_lloyd),
        },
        Algoritmo {
            id: "agrupamiento.codo",
            capitulo: "agrupamiento",
            nombre: "Elegir k con la curva del codo",
            resumen: "Corre k-means para k = 1, 2, 3… y grafica la distorsión: donde la curva deja de caer fuerte (el codo) suele estar la cantidad natural de grupos.",
            campos: &[
                DATOS_EXPRESION,
                Campo {
                    id: "k_maximo",
                    etiqueta: "k máximo a probar (hasta 15)",
                    tipo: Tipo::Entero,
                    ejemplo: "10",
                },
                SEMILLA,
            ],
            funcion: "agrupamiento::curva_del_codo",
            ejecutar: codo,
            pasos: Some(pasos_codo),
        },
        Algoritmo {
            id: "agrupamiento.k_medias_suaves",
            capitulo: "agrupamiento",
            nombre: "k-means suave",
            resumen: "Cada punto pertenece a todos los grupos a la vez, con una responsabilidad que depende de su distancia a cada centro y de la rigidez β.",
            campos: &[
                DATOS_EXPRESION,
                K_TRES,
                Campo {
                    id: "rigidez",
                    etiqueta: "Rigidez β (más grande, más parecido a k-means duro)",
                    tipo: Tipo::Decimal,
                    ejemplo: "1",
                },
                ITERACIONES,
            ],
            funcion: "agrupamiento::k_medias_suaves",
            ejecutar: k_medias_suaves,
            pasos: Some(pasos_k_medias_suaves),
        },
        Algoritmo {
            id: "agrupamiento.jerarquico",
            capitulo: "agrupamiento",
            nombre: "Agrupamiento jerárquico",
            resumen: "Une una y otra vez los dos grupos más cercanos, con la distancia promedio entre sus elementos, y dibuja el árbol de uniones.",
            campos: &[
                DATOS_EXPRESION,
                Campo {
                    id: "distancia",
                    etiqueta: "Distancia entre filas",
                    tipo: Tipo::Opcion(MEDIDAS),
                    ejemplo: MEDIDA_PEARSON,
                },
            ],
            funcion: "agrupamiento::agrupamiento_jerarquico",
            ejecutar: jerarquico,
            pasos: Some(pasos_jerarquico),
        },
        Algoritmo {
            id: "agrupamiento.cast",
            capitulo: "agrupamiento",
            nombre: "CAST: grupos por afinidad",
            resumen: "Arma grupos de genes cuya correlación promedio con el grupo llega a un umbral θ, sin decidir de antemano cuántos grupos habrá.",
            campos: &[
                DATOS_EXPRESION,
                Campo {
                    id: "umbral",
                    etiqueta: "Umbral θ de correlación (entre 0 y 1)",
                    tipo: Tipo::Decimal,
                    ejemplo: "0.7",
                },
            ],
            funcion: "agrupamiento::cast",
            ejecutar: cast,
            pasos: Some(pasos_cast),
        },
        Algoritmo {
            id: "agrupamiento.componentes_principales",
            capitulo: "agrupamiento",
            nombre: "Componentes principales (vista en 2D)",
            resumen: "Encuentra, con el método de la potencia, las dos direcciones en que más varían los datos y proyecta cada gen sobre ellas para verlos en un plano.",
            campos: &[DATOS_EXPRESION],
            funcion: "agrupamiento::componentes_principales",
            ejecutar: componentes_principales,
            pasos: Some(pasos_componentes),
        },
    ]
}

struct Tabla {
    nombres: Vec<String>,
    columnas: Vec<String>,
    filas: Vec<Vec<f64>>,
}

impl Tabla {
    fn cantidad(&self) -> usize {
        self.filas.len()
    }

    fn dimension(&self) -> usize {
        self.columnas.len()
    }
}

struct Plano {
    coordenadas: Vec<[f64; 2]>,
    ejes: [String; 2],
    componentes: Option<Componentes>,
}

#[derive(Clone, Copy, PartialEq)]
enum Patron {
    Sube,
    Baja,
    Pico,
    Valle,
    Plano,
}

fn leer_tabla(entrada: &Entrada, campo: &str, prefijo: &str) -> Result<Tabla, String> {
    let texto = entrada.texto(campo)?;
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
            (format!("{prefijo}{}", filas.len() + 1), &piezas[..])
        } else {
            (piezas[0].to_string(), &piezas[1..])
        };
        if valores.is_empty() {
            return Err(format!(
                "En «{campo}», la línea {} («{}») no tiene números.",
                numero + 1,
                recortar(linea.trim(), 30)
            ));
        }
        let mut fila = Vec::with_capacity(valores.len());
        for pieza in valores {
            let valor: f64 = pieza.parse().map_err(|_| {
                format!(
                    "En «{campo}», línea {}: «{}» no es un número.",
                    numero + 1,
                    recortar(pieza, 30)
                )
            })?;
            if !valor.is_finite() {
                return Err(format!(
                    "En «{campo}», línea {}: «{pieza}» no es un número finito.",
                    numero + 1
                ));
            }
            if valor.abs() > MAXIMO_VALOR {
                return Err(format!(
                    "En «{campo}», línea {}: «{pieza}» se sale del rango que acepta este capítulo (de −10¹² a 10¹²).",
                    numero + 1
                ));
            }
            fila.push(valor);
        }
        if fila.len() > MAXIMO_COLUMNAS {
            return Err(format!(
                "En «{campo}», la línea {} tiene {} números; el máximo es {MAXIMO_COLUMNAS} columnas.",
                numero + 1,
                fila.len()
            ));
        }
        if let Some(primera) = filas.first()
            && primera.len() != fila.len()
        {
            return Err(format!(
                "En «{campo}», la línea {} ({nombre}) tiene {} números y la primera fila tiene {}: todas las filas deben tener la misma cantidad.",
                numero + 1,
                fila.len(),
                primera.len()
            ));
        }
        filas.push(fila);
        nombres.push(nombre);
        if filas.len() > MAXIMO_PUNTOS {
            return Err(format!(
                "«{campo}» tiene más de {} filas, el máximo que acepta este capítulo.",
                miles(MAXIMO_PUNTOS)
            ));
        }
    }
    if filas.is_empty() {
        return Err(format!("«{campo}» no tiene ninguna fila con números."));
    }
    let dimension = filas[0].len();
    let columnas = match encabezado {
        Some(titulos) if titulos.len() == dimension + 1 => titulos[1..].to_vec(),
        Some(titulos) if titulos.len() == dimension => titulos,
        _ => (1..=dimension)
            .map(|columna| format!("c{columna}"))
            .collect(),
    };
    Ok(Tabla {
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

fn leer_k(entrada: &Entrada, campo: &str, cantidad: usize) -> Result<usize, String> {
    let k = entrada.natural(campo)?;
    if k == 0 || k > MAXIMO_K {
        return Err(format!("k debe estar entre 1 y {MAXIMO_K}."));
    }
    if k > cantidad {
        return Err(format!(
            "k = {k} es mayor que la cantidad de filas ({cantidad}): no se pueden formar más grupos que puntos."
        ));
    }
    Ok(k)
}

fn leer_iteraciones(entrada: &Entrada, campo: &str) -> Result<usize, String> {
    let iteraciones = entrada.natural(campo)?;
    if iteraciones == 0 || iteraciones > MAXIMO_ITERACIONES {
        return Err(format!(
            "Las iteraciones deben estar entre 1 y {}.",
            miles(MAXIMO_ITERACIONES)
        ));
    }
    Ok(iteraciones)
}

fn leer_semilla(entrada: &Entrada) -> Result<u64, String> {
    Ok(entrada.natural("semilla")? as u64)
}

fn leer_inicio(entrada: &Entrada) -> Result<Inicio, String> {
    let eleccion = entrada.opcion("inicio", INICIOS)?;
    let semilla = leer_semilla(entrada)?;
    if eleccion == K_MEDIAS_MAS_MAS {
        Ok(Inicio::KMediasMasMas(semilla))
    } else {
        Ok(Inicio::PrimerosPuntos)
    }
}

fn exigir_filas(tabla: &Tabla, minimo: usize) -> Result<(), String> {
    if tabla.cantidad() < minimo {
        return Err(format!(
            "Hacen falta al menos {minimo} filas y la tabla tiene {}.",
            tabla.cantidad()
        ));
    }
    Ok(())
}

fn limitar_matriz(tabla: &Tabla) -> Result<(), String> {
    if tabla.cantidad() > MAXIMO_MATRIZ {
        return Err(format!(
            "Este método compara todas las filas entre sí y admite hasta {} filas; la tabla tiene {}.",
            miles(MAXIMO_MATRIZ),
            miles(tabla.cantidad())
        ));
    }
    Ok(())
}

fn limitar_cuentas(cuentas: f64) -> Result<(), String> {
    if cuentas > MAXIMO_CUENTAS {
        return Err(format!(
            "Esta combinación pide unos {} millones de operaciones y el máximo es {} millones: baja k, las iteraciones o la cantidad de filas.",
            miles((cuentas / 1e6) as usize),
            miles((MAXIMO_CUENTAS / 1e6) as usize)
        ));
    }
    Ok(())
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

fn redondear(valor: f64) -> f64 {
    (valor * 10_000.0).round() / 10_000.0
}

fn plano_de_dibujo(tabla: &Tabla) -> Result<Plano, String> {
    if tabla.dimension() == 1 {
        return Ok(Plano {
            coordenadas: tabla.filas.iter().map(|fila| [fila[0], 0.0]).collect(),
            ejes: [tabla.columnas[0].clone(), String::new()],
            componentes: None,
        });
    }
    if tabla.dimension() == 2 {
        return Ok(Plano {
            coordenadas: tabla.filas.iter().map(|fila| [fila[0], fila[1]]).collect(),
            ejes: [tabla.columnas[0].clone(), tabla.columnas[1].clone()],
            componentes: None,
        });
    }
    let componentes = agrupamiento::componentes_principales(&tabla.filas)?;
    let total = componentes.varianza_total;
    let parte = |varianza: f64| {
        if total > 0.0 {
            porcentaje(varianza / total, 0)
        } else {
            porcentaje(0.0, 0)
        }
    };
    Ok(Plano {
        coordenadas: componentes.coordenadas.clone(),
        ejes: [
            format!(
                "Componente principal 1 ({} de la varianza)",
                parte(componentes.varianza_primera)
            ),
            format!(
                "Componente principal 2 ({} de la varianza)",
                parte(componentes.varianza_segunda)
            ),
        ],
        componentes: Some(componentes),
    })
}

fn ubicar(plano: &Plano, punto: &[f64]) -> [f64; 2] {
    match &plano.componentes {
        Some(componentes) => agrupamiento::proyectar(punto, componentes),
        None => [punto[0], punto.get(1).copied().unwrap_or(0.0)],
    }
}

fn dispersion(etiqueta: &str, plano: &Plano, nombres: &[String], grupos: &[usize]) -> Bloque {
    let puntos = muestra_pareja(nombres.len(), MAXIMO_PUNTOS_DIBUJADOS)
        .into_iter()
        .map(|indice| Punto {
            x: plano.coordenadas[indice][0],
            y: plano.coordenadas[indice][1],
            grupo: grupos[indice],
            nombre: nombres[indice].clone(),
        })
        .collect();
    Bloque::dispersion(etiqueta, &plano.ejes[0], &plano.ejes[1], puntos)
}

fn orden_por_grupos(grupos: &[usize]) -> Vec<usize> {
    let mut orden: Vec<usize> = (0..grupos.len()).collect();
    orden.sort_by_key(|&indice| (grupos[indice], indice));
    orden
}

fn calor_ordenado(
    etiqueta: &str,
    tabla: &Tabla,
    orden: &[usize],
    grupos: Option<&[usize]>,
) -> Bloque {
    let elegidos: Vec<usize> = muestra_pareja(orden.len(), MAXIMO_FILAS_DE_CALOR)
        .into_iter()
        .map(|posicion| orden[posicion])
        .collect();
    let filas = elegidos
        .iter()
        .map(|&indice| match grupos {
            Some(grupos) => format!("{} · grupo {}", tabla.nombres[indice], grupos[indice] + 1),
            None => tabla.nombres[indice].clone(),
        })
        .collect();
    let valores = elegidos
        .iter()
        .map(|&indice| tabla.filas[indice].clone())
        .collect();
    Bloque::calor(etiqueta, filas, tabla.columnas.clone(), valores)
}

fn clasificar(perfil: &[f64]) -> Patron {
    if perfil.len() < 2 {
        return Patron::Plano;
    }
    let tiempos: Vec<f64> = (0..perfil.len()).map(|tiempo| tiempo as f64).collect();
    let tendencia = agrupamiento::correlacion_pearson(perfil, &tiempos);
    if tendencia >= CORRELACION_FUERTE {
        return Patron::Sube;
    }
    if tendencia <= -CORRELACION_FUERTE {
        return Patron::Baja;
    }
    let (minimo, donde_minimo) = extremo(perfil, false);
    let (maximo, donde_maximo) = extremo(perfil, true);
    let rango = maximo - minimo;
    let ultimo = perfil.len() - 1;
    let bordes_altos = perfil[0].max(perfil[ultimo]);
    let bordes_bajos = perfil[0].min(perfil[ultimo]);
    let en_el_medio = |donde: usize| donde > 0 && donde < ultimo;
    if rango > 0.0 && en_el_medio(donde_maximo) && maximo - bordes_altos >= rango / 2.0 {
        return Patron::Pico;
    }
    if rango > 0.0 && en_el_medio(donde_minimo) && bordes_bajos - minimo >= rango / 2.0 {
        return Patron::Valle;
    }
    Patron::Plano
}

fn extremo(valores: &[f64], buscar_maximo: bool) -> (f64, usize) {
    let mut mejor = 0;
    for (indice, &valor) in valores.iter().enumerate() {
        let gana = if buscar_maximo {
            valor > valores[mejor]
        } else {
            valor < valores[mejor]
        };
        if gana {
            mejor = indice;
        }
    }
    (valores[mejor], mejor)
}

fn describir_perfil(perfil: &[f64], columnas: &[String]) -> String {
    let primero = una_cifra(perfil[0]);
    let ultimo = una_cifra(perfil[perfil.len() - 1]);
    let (minimo, donde_minimo) = extremo(perfil, false);
    let (maximo, donde_maximo) = extremo(perfil, true);
    match clasificar(perfil) {
        Patron::Sube => format!("sube (de {primero} a {ultimo})"),
        Patron::Baja => format!("baja (de {primero} a {ultimo})"),
        Patron::Pico => format!(
            "tiene un pico pasajero (llega a {} en {})",
            una_cifra(maximo),
            columnas[donde_maximo]
        ),
        Patron::Valle => format!(
            "tiene una caída pasajera (baja a {} en {})",
            una_cifra(minimo),
            columnas[donde_minimo]
        ),
        Patron::Plano => format!(
            "no tiene una tendencia clara (se mueve entre {} y {})",
            una_cifra(minimo),
            una_cifra(maximo)
        ),
    }
}

fn una_cifra(valor: f64) -> String {
    let texto = decimal(valor, 1);
    if texto == "-0.0" {
        "0.0".to_string()
    } else {
        texto
    }
}

fn cantidades_por_grupo(grupos: &[usize], cuantos: usize) -> Vec<usize> {
    let mut cantidades = vec![0; cuantos];
    for &grupo in grupos {
        cantidades[grupo] += 1;
    }
    cantidades
}

fn miembros_de(tabla: &Tabla, grupos: &[usize], grupo: usize) -> Vec<String> {
    grupos
        .iter()
        .enumerate()
        .filter(|(_, elegido)| **elegido == grupo)
        .map(|(indice, _)| tabla.nombres[indice].clone())
        .collect()
}

fn tabla_de_grupos(tabla: &Tabla, grupos: &[usize], centros: &[Vec<f64>]) -> Bloque {
    let con_tiempos = tabla.dimension() >= MINIMO_TIEMPOS;
    let filas = centros
        .iter()
        .enumerate()
        .map(|(grupo, centro)| {
            let miembros = miembros_de(tabla, grupos, grupo);
            let patron = if con_tiempos {
                describir_perfil(centro, &tabla.columnas)
            } else {
                "—".to_string()
            };
            vec![
                format!("{}", grupo + 1),
                miles(miembros.len()),
                patron,
                algunos(&miembros, NOMBRES_POR_GRUPO),
            ]
        })
        .collect();
    Bloque::tabla(
        "Grupos",
        &["grupo", "filas", "perfil promedio", "miembros"],
        filas,
    )
}

fn calor_de_centros(tabla: &Tabla, centros: &[Vec<f64>]) -> Bloque {
    Bloque::calor(
        "Perfil promedio de cada grupo (su centro)",
        (1..=centros.len())
            .map(|grupo| format!("grupo {grupo}"))
            .collect(),
        tabla.columnas.clone(),
        centros.to_vec(),
    )
}

fn centros_de(tabla: &Tabla, grupos: &[usize], cuantos: usize) -> Vec<Vec<f64>> {
    let vacios = vec![vec![0.0; tabla.dimension()]; cuantos];
    agrupamiento::de_grupos_a_centros(&tabla.filas, grupos, &vacios)
}

fn explicar_grupos(salida: &mut Salida, tabla: &Tabla, grupos: &[usize], centros: &[Vec<f64>]) {
    let cantidades = cantidades_por_grupo(grupos, centros.len());
    let cuentas: Vec<String> = cantidades.iter().map(|&cantidad| miles(cantidad)).collect();
    let mut texto = format!(
        "Las {} filas quedaron en {} {} de {} elementos.",
        miles(tabla.cantidad()),
        centros.len(),
        segun_cantidad(centros.len(), "grupo", "grupos"),
        unir(&cuentas)
    );
    let vacios: Vec<String> = cantidades
        .iter()
        .enumerate()
        .filter(|(_, cantidad)| **cantidad == 0)
        .map(|(grupo, _)| format!("{}", grupo + 1))
        .collect();
    if !vacios.is_empty() {
        texto.push_str(&format!(
            " {} {} {}: ninguna fila tiene ese centro como el más cercano; prueba con un k menor.",
            segun_cantidad(vacios.len(), "El grupo", "Los grupos"),
            unir(&vacios),
            segun_cantidad(vacios.len(), "quedó vacío", "quedaron vacíos")
        ));
    }
    let solitarios: Vec<String> = cantidades
        .iter()
        .enumerate()
        .filter(|(_, cantidad)| **cantidad == 1)
        .flat_map(|(grupo, _)| miembros_de(tabla, grupos, grupo))
        .collect();
    if !solitarios.is_empty() && tabla.cantidad() > centros.len() {
        texto.push_str(&format!(
            " {} {}: puede ser una fila atípica o una señal de que k es demasiado grande.",
            algunos(&solitarios, NOMBRES_A_LA_VISTA),
            segun_cantidad(
                solitarios.len(),
                "quedó solo en su grupo",
                "quedaron solos en sus grupos"
            )
        ));
    }
    salida.explicar(texto);
    explicar_patrones(salida, tabla, grupos, centros);
}

fn explicar_patrones(salida: &mut Salida, tabla: &Tabla, grupos: &[usize], centros: &[Vec<f64>]) {
    if tabla.dimension() < MINIMO_TIEMPOS {
        return;
    }
    let cantidades = cantidades_por_grupo(grupos, centros.len());
    let mut partes = Vec::new();
    let mut con_patron = 0;
    for (grupo, centro) in centros.iter().enumerate() {
        if cantidades[grupo] == 0 {
            continue;
        }
        if clasificar(centro) != Patron::Plano && cantidades[grupo] >= 2 {
            con_patron += 1;
        }
        partes.push(format!(
            "el grupo {} {}",
            grupo + 1,
            describir_perfil(centro, &tabla.columnas)
        ));
    }
    salida.explicar(format!(
        "Si las columnas son tiempos en orden, el perfil promedio dice que {}. Uso la correlación con el tiempo (≥ 0.8 sube, ≤ −0.8 baja) y si el máximo o el mínimo cae en el medio.",
        unir(&partes)
    ));
    if con_patron > 0 {
        salida.explicar(
            "Genes del mismo grupo cambian juntos: en biología eso sugiere que responden a la misma señal o al mismo regulador, pero es una hipótesis para comprobar, no una prueba.",
        );
    }
}

fn distancias(entrada: &Entrada) -> Result<Salida, String> {
    let tabla = leer_tabla(entrada, "datos", "p")?;
    exigir_filas(&tabla, 2)?;
    limitar_matriz(&tabla)?;
    let euclidea = agrupamiento::matriz_de_distancias(&tabla.filas, Medida::Euclidea);
    let pearson = agrupamiento::matriz_de_distancias(&tabla.filas, Medida::Pearson);
    let cantidad = tabla.cantidad();
    let mut filas = Vec::with_capacity(cantidad);
    let mut distintos = Vec::new();
    for fila in 0..cantidad {
        let cercano_euclideo = mas_cercano(&euclidea[fila], fila);
        let cercano_pearson = mas_cercano(&pearson[fila], fila);
        filas.push(vec![
            tabla.nombres[fila].clone(),
            tabla.nombres[cercano_euclideo].clone(),
            decimal(euclidea[fila][cercano_euclideo], 3),
            tabla.nombres[cercano_pearson].clone(),
            decimal(pearson[fila][cercano_pearson], 3),
            decimal(1.0 - pearson[fila][cercano_pearson], 3),
        ]);
        if cercano_euclideo != cercano_pearson {
            distintos.push((fila, cercano_euclideo, cercano_pearson));
        }
    }
    let mut salida = Salida::default();
    salida.agregar(Bloque::tabla(
        "Vecino más cercano de cada fila según cada distancia",
        &[
            "fila",
            "más cercana (euclídea)",
            "distancia euclídea",
            "más cercana (Pearson)",
            "distancia de Pearson (1 − r)",
            "correlación r",
        ],
        filas,
    ));
    if cantidad <= MAXIMO_CALOR_DE_DISTANCIAS {
        salida.agregar(Bloque::calor(
            "Distancia euclídea entre cada par",
            tabla.nombres.clone(),
            tabla.nombres.clone(),
            euclidea.clone(),
        ));
        salida.agregar(Bloque::calor(
            "Distancia de Pearson (1 − r) entre cada par",
            tabla.nombres.clone(),
            tabla.nombres.clone(),
            pearson.clone(),
        ));
    }
    if tabla.dimension() <= 2 {
        salida.explicar(format!(
            "Cada fila tiene solo {} valores: con tan pocos, la correlación de Pearson da casi siempre +1 o −1 y no dice nada. Hacen falta al menos 3 condiciones.",
            tabla.dimension()
        ));
    }
    if let Some(&(fila, euclideo, por_forma)) = distintos.first() {
        salida.explicar(format!(
            "En {} de {} filas el vecino más cercano cambia según la distancia. Por ejemplo, para {} la euclídea elige {} (valores parecidos) y Pearson elige {} (la misma forma, aunque cambie más o menos). La euclídea mira la magnitud; Pearson, solo la forma del perfil.",
            miles(distintos.len()),
            miles(cantidad),
            tabla.nombres[fila],
            tabla.nombres[euclideo],
            tabla.nombres[por_forma]
        ));
    } else {
        salida.explicar(
            "Las dos distancias eligen el mismo vecino más cercano para cada fila: aquí la magnitud y la forma cuentan la misma historia.",
        );
    }
    let (opuesto_a, opuesto_b, mas_lejos) = par_extremo(&pearson, true);
    if 1.0 - mas_lejos <= -CORRELACION_FUERTE {
        salida.explicar(format!(
            "{} y {} hacen lo contrario (r = {}): cuando uno sube, el otro baja. Pearson los pone lejos (cerca de 2), aunque podrían depender del mismo regulador, que activa a uno y reprime al otro.",
            tabla.nombres[opuesto_a],
            tabla.nombres[opuesto_b],
            decimal(1.0 - mas_lejos, 2)
        ));
    }
    let (parecido_a, parecido_b, mas_cerca) = par_extremo(&pearson, false);
    if 1.0 - mas_cerca >= CORRELACION_FUERTE {
        salida.explicar(format!(
            "El par con la forma más parecida es {} y {} (r = {}); una correlación de 0.8 o más se considera fuerte y los vuelve candidatos a estar regulados juntos.",
            tabla.nombres[parecido_a],
            tabla.nombres[parecido_b],
            decimal(1.0 - mas_cerca, 2)
        ));
    }
    let constantes: Vec<String> = tabla
        .filas
        .iter()
        .enumerate()
        .filter(|(_, fila)| fila.iter().all(|valor| *valor == fila[0]))
        .map(|(indice, _)| tabla.nombres[indice].clone())
        .collect();
    if !constantes.is_empty() {
        salida.explicar(format!(
            "{} no cambia en ninguna condición: su correlación con las demás no está definida y aquí vale 0 (distancia de Pearson 1).",
            algunos(&constantes, NOMBRES_A_LA_VISTA)
        ));
    }
    Ok(salida)
}

fn mas_cercano(distancias: &[f64], propio: usize) -> usize {
    let mut mejor = if propio == 0 { 1 } else { 0 };
    for (otro, &lejania) in distancias.iter().enumerate() {
        if otro != propio && lejania < distancias[mejor] {
            mejor = otro;
        }
    }
    mejor
}

fn par_extremo(distancias: &[Vec<f64>], buscar_maximo: bool) -> (usize, usize, f64) {
    let mut mejor = (0, 1, distancias[0][1]);
    for (fila, valores) in distancias.iter().enumerate() {
        for (columna, &valor) in valores.iter().enumerate().skip(fila + 1) {
            let gana = if buscar_maximo {
                valor > mejor.2
            } else {
                valor < mejor.2
            };
            if gana {
                mejor = (fila, columna, valor);
            }
        }
    }
    mejor
}

fn mas_lejano_primero(entrada: &Entrada) -> Result<Salida, String> {
    let tabla = leer_tabla(entrada, "datos", "p")?;
    let k = leer_k(entrada, "k", tabla.cantidad())?;
    let elegidos = agrupamiento::mas_lejano_primero(&tabla.filas, k)?;
    let centros: Vec<Vec<f64>> = elegidos
        .iter()
        .map(|&indice| tabla.filas[indice].clone())
        .collect();
    let grupos = agrupamiento::de_centros_a_grupos(&tabla.filas, &centros);
    let radio = agrupamiento::radio_de_los_centros(&tabla.filas, &centros);
    let cantidades = cantidades_por_grupo(&grupos, k);
    let mut filas = Vec::with_capacity(k);
    for (orden, &indice) in elegidos.iter().enumerate() {
        let al_elegirlo = if orden == 0 {
            "— (es el primer punto)".to_string()
        } else {
            let mut menor = f64::INFINITY;
            for anterior in &centros[..orden] {
                menor = menor.min(agrupamiento::distancia_euclidea(
                    &tabla.filas[indice],
                    anterior,
                ));
            }
            decimal(menor, 3)
        };
        let mut radio_del_grupo = 0.0f64;
        for (fila, &grupo) in grupos.iter().enumerate() {
            if grupo == orden {
                radio_del_grupo = radio_del_grupo.max(agrupamiento::distancia_euclidea(
                    &tabla.filas[fila],
                    &centros[orden],
                ));
            }
        }
        filas.push(vec![
            format!("{}", orden + 1),
            tabla.nombres[indice].clone(),
            al_elegirlo,
            miles(cantidades[orden]),
            decimal(radio_del_grupo, 3),
        ]);
    }
    let mut nombres = tabla.nombres.clone();
    for (orden, &indice) in elegidos.iter().enumerate() {
        nombres[indice] = format!("{} (centro {})", tabla.nombres[indice], orden + 1);
    }
    let plano = plano_de_dibujo(&tabla)?;
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        "Radio: la mayor distancia de un punto a su centro",
        decimal(radio, 3),
    ));
    salida.agregar(Bloque::tabla(
        "Centros en el orden en que se eligieron",
        &[
            "orden",
            "centro",
            "distancia a los centros anteriores",
            "puntos en su grupo",
            "radio del grupo",
        ],
        filas,
    ));
    salida.agregar(dispersion(
        "Cada punto, con el color de su centro más cercano",
        &plano,
        &nombres,
        &grupos,
    ));
    if tabla.dimension() >= 3 {
        salida.agregar(calor_ordenado(
            "Filas ordenadas por grupo",
            &tabla,
            &orden_por_grupos(&grupos),
            Some(&grupos),
        ));
    }
    salida.explicar(format!(
        "El primer centro es siempre la primera fila ({}); cada centro nuevo es el punto más lejano de los ya elegidos, así que el resultado depende del orden de las filas.",
        tabla.nombres[0]
    ));
    salida.explicar(format!(
        "Ningún punto está a más de {} de su centro: ese radio es lo que k-centros intenta achicar.",
        decimal(radio, 3)
    ));
    let atipicos: Vec<String> = elegidos
        .iter()
        .enumerate()
        .filter(|(orden, _)| cantidades[*orden] == 1)
        .map(|(_, &indice)| tabla.nombres[indice].clone())
        .collect();
    if !atipicos.is_empty() && tabla.cantidad() > k {
        salida.explicar(format!(
            "{} {}: {} lejos de todos. k-centros elige estos puntos porque mira la peor distancia; k-means, que promedia, se deja llevar menos por un solo punto raro.",
            unir(&atipicos),
            segun_cantidad(
                atipicos.len(),
                "quedó solo en su grupo",
                "quedaron solos en sus grupos"
            ),
            segun_cantidad(atipicos.len(), "es un punto atípico,", "son puntos atípicos,")
        ));
        if k < MAXIMO_K {
            salida.explicar(format!(
                "Como un centro se gastó en el atípico, otras nubes pueden haber quedado juntas en un mismo grupo: prueba con k = {}.",
                k + 1
            ));
        }
    }
    Ok(salida)
}

fn distorsion(entrada: &Entrada) -> Result<Salida, String> {
    let tabla = leer_tabla(entrada, "datos", "p")?;
    let centros = leer_tabla(entrada, "centros", "centro ")?;
    if centros.dimension() != tabla.dimension() {
        return Err(format!(
            "Los puntos tienen {} coordenadas y los centros {}: deben tener las mismas.",
            tabla.dimension(),
            centros.dimension()
        ));
    }
    if centros.cantidad() > MAXIMO_K {
        return Err(format!(
            "Hay {} centros; el máximo es {MAXIMO_K}.",
            centros.cantidad()
        ));
    }
    let valor = agrupamiento::distorsion_cuadratica(&tabla.filas, &centros.filas);
    let grupos = agrupamiento::de_centros_a_grupos(&tabla.filas, &centros.filas);
    let cantidades = cantidades_por_grupo(&grupos, centros.cantidad());
    let mut aportes = vec![0.0; centros.cantidad()];
    for (punto, &grupo) in tabla.filas.iter().zip(&grupos) {
        aportes[grupo] += agrupamiento::distancia_al_cuadrado(punto, &centros.filas[grupo]);
    }
    let total = tabla.cantidad() as f64;
    let filas = (0..centros.cantidad())
        .map(|centro| {
            let promedio = if cantidades[centro] > 0 {
                cientifico(aportes[centro] / cantidades[centro] as f64)
            } else {
                "—".to_string()
            };
            vec![
                centros.nombres[centro].clone(),
                miles(cantidades[centro]),
                cientifico(aportes[centro] / total),
                promedio,
            ]
        })
        .collect();
    let movidos = agrupamiento::de_grupos_a_centros(&tabla.filas, &grupos, &centros.filas);
    let mejor = agrupamiento::distorsion_cuadratica(&tabla.filas, &movidos);
    let plano = plano_de_dibujo(&tabla)?;
    let mut nombres = tabla.nombres.clone();
    let mut grupos_dibujados = grupos.clone();
    let mut dibujo = Plano {
        coordenadas: plano.coordenadas.clone(),
        ejes: plano.ejes.clone(),
        componentes: None,
    };
    for (indice, centro) in centros.filas.iter().enumerate() {
        dibujo.coordenadas.push(ubicar(&plano, centro));
        nombres.push(format!("{} (centro)", centros.nombres[indice]));
        grupos_dibujados.push(indice);
    }
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Distorsión", cientifico(valor)));
    salida.agregar(Bloque::valor(
        "Distancia típica a su centro (raíz de la distorsión)",
        cientifico(valor.sqrt()),
    ));
    salida.agregar(Bloque::tabla(
        "Cada centro y los puntos que tiene más cerca",
        &[
            "centro",
            "puntos",
            "aporte a la distorsión",
            "distancia al cuadrado promedio en el grupo",
        ],
        filas,
    ));
    salida.agregar(dispersion(
        "Puntos y centros (cada punto con el color de su centro más cercano)",
        &dibujo,
        &nombres,
        &grupos_dibujados,
    ));
    salida.explicar(format!(
        "La distorsión es {}: en promedio, cada punto está a una distancia al cuadrado de {} de su centro más cercano, o sea a unos {} en las unidades de tus datos.",
        cientifico(valor),
        cientifico(valor),
        cientifico(valor.sqrt())
    ));
    let vacios: Vec<String> = (0..centros.cantidad())
        .filter(|&centro| cantidades[centro] == 0)
        .map(|centro| centros.nombres[centro].clone())
        .collect();
    if !vacios.is_empty() {
        salida.explicar(format!(
            "{} {}",
            unir(&vacios),
            segun_cantidad(
                vacios.len(),
                "no es el centro más cercano de ningún punto: no aporta nada y podrías quitarlo.",
                "no son el centro más cercano de ningún punto: no aportan nada y podrías quitarlos."
            )
        ));
    }
    if valor > 0.0 && mejor < valor * (1.0 - TOLERANCIA) {
        salida.explicar(format!(
            "Estos centros no son los mejores para sus grupos: si mueves cada uno al centro de gravedad de sus puntos (un paso de Lloyd), la distorsión baja de {} a {} ({} menos).",
            cientifico(valor),
            cientifico(mejor),
            porcentaje((valor - mejor) / valor, 0)
        ));
    } else {
        salida.explicar(
            "Cada centro ya está en el centro de gravedad de sus puntos: un paso de Lloyd no los movería.",
        );
    }
    if centros.cantidad() >= 2 && valor > 0.0 {
        let (mayor, donde) = extremo(&aportes, true);
        let parte = mayor / total / valor;
        if parte > 0.5 {
            salida.explicar(format!(
                "{} aporta {} de la distorsión: su grupo es el más disperso, y es el primero que conviene revisar o partir en dos.",
                centros.nombres[donde],
                porcentaje(parte, 0)
            ));
        }
    }
    Ok(salida)
}

fn lloyd(entrada: &Entrada) -> Result<Salida, String> {
    let tabla = leer_tabla(entrada, "datos", "p")?;
    exigir_filas(&tabla, 2)?;
    let k = leer_k(entrada, "k", tabla.cantidad())?;
    let inicio = leer_inicio(entrada)?;
    let maximo = leer_iteraciones(entrada, "iteraciones")?;
    let resultado = agrupamiento::k_medias(&tabla.filas, k, inicio, maximo)?;
    let plano = plano_de_dibujo(&tabla)?;
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        "Distorsión final",
        cientifico(resultado.distorsion),
    ));
    salida.agregar(Bloque::valor(
        "Iteraciones",
        if resultado.convergio {
            format!("{} (se estabilizó)", miles(resultado.iteraciones))
        } else {
            format!("{} (llegó al máximo)", miles(resultado.iteraciones))
        },
    ));
    salida.agregar(tabla_de_grupos(
        &tabla,
        &resultado.grupos,
        &resultado.centros,
    ));
    salida.agregar(dispersion(
        "Cada fila con el color de su grupo",
        &plano,
        &tabla.nombres,
        &resultado.grupos,
    ));
    salida.agregar(calor_ordenado(
        "Filas ordenadas por grupo",
        &tabla,
        &orden_por_grupos(&resultado.grupos),
        Some(&resultado.grupos),
    ));
    salida.agregar(calor_de_centros(&tabla, &resultado.centros));
    explicar_grupos(&mut salida, &tabla, &resultado.grupos, &resultado.centros);
    if resultado.convergio {
        salida.explicar(format!(
            "Lloyd se detuvo tras {} {} porque ningún punto cambió de grupo; la distorsión final es {}.",
            miles(resultado.iteraciones),
            segun_cantidad(resultado.iteraciones, "iteración", "iteraciones"),
            cientifico(resultado.distorsion)
        ));
    } else {
        salida.explicar(format!(
            "Lloyd llegó al máximo de {} {} y los grupos todavía cambiaban: sube las iteraciones.",
            miles(maximo),
            segun_cantidad(maximo, "iteración", "iteraciones")
        ));
    }
    let semilla = leer_semilla(entrada)?;
    let (otra, nombre_otra) = match inicio {
        Inicio::PrimerosPuntos => (
            Inicio::KMediasMasMas(semilla),
            format!("k-means++ con semilla {semilla}"),
        ),
        Inicio::KMediasMasMas(_) => (Inicio::PrimerosPuntos, PRIMEROS_PUNTOS.to_string()),
    };
    let alternativa = agrupamiento::k_medias(&tabla.filas, k, otra, maximo)?;
    if alternativa.distorsion < resultado.distorsion * (1.0 - MEJORA_NOTABLE) {
        salida.explicar(format!(
            "Lloyd solo garantiza un mínimo local: empezando con {nombre_otra} la distorsión baja a {} ({} menos). Prueba otros centros iniciales.",
            cientifico(alternativa.distorsion),
            porcentaje(
                (resultado.distorsion - alternativa.distorsion) / resultado.distorsion,
                0
            )
        ));
    }
    Ok(salida)
}

struct Momento {
    fase: &'static str,
    linea: &'static str,
    nota: String,
    iteracion: usize,
    grupos: Option<Vec<usize>>,
    centros: Vec<[f64; 2]>,
    anteriores: Option<Vec<[f64; 2]>>,
    cambiaron: Option<usize>,
    distorsion: f64,
    terminado: bool,
}

fn proyectar_centros(plano: &Plano, centros: &[Vec<f64>]) -> Vec<[f64; 2]> {
    centros
        .iter()
        .map(|centro| {
            let lugar = ubicar(plano, centro);
            [redondear(lugar[0]), redondear(lugar[1])]
        })
        .collect()
}

fn distorsion_de_los_grupos(puntos: &[Vec<f64>], grupos: &[usize], centros: &[Vec<f64>]) -> f64 {
    if puntos.is_empty() {
        return 0.0;
    }
    let mut suma = 0.0;
    for (punto, &grupo) in puntos.iter().zip(grupos) {
        suma += agrupamiento::distancia_al_cuadrado(punto, &centros[grupo]);
    }
    suma / puntos.len() as f64
}

fn pasos_lloyd(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let tabla = leer_tabla(entrada, "datos", "p")?;
    exigir_filas(&tabla, 2)?;
    let k = leer_k(entrada, "k", tabla.cantidad())?;
    let inicio = leer_inicio(entrada)?;
    let maximo = leer_iteraciones(entrada, "iteraciones")?;
    if tabla.cantidad() > MAXIMO_ANIMABLE {
        return Err(format!(
            "El paso a paso dibuja hasta {} filas y la tabla tiene {}; para tablas más grandes usa «Ejecutar».",
            miles(MAXIMO_ANIMABLE),
            miles(tabla.cantidad())
        ));
    }
    let puntos = &tabla.filas;
    let plano = plano_de_dibujo(&tabla)?;
    let mut momentos = Vec::new();
    let mut centros = agrupamiento::centros_iniciales(puntos, k, inicio);
    let (linea_inicial, nota_inicial) = match inicio {
        Inicio::PrimerosPuntos => (
            LINEA_PRIMEROS,
            format!(
                "Los centros iniciales son las primeras {k} filas: {}.",
                algunos(&tabla.nombres[..k], NOMBRES_A_LA_VISTA)
            ),
        ),
        Inicio::KMediasMasMas(semilla) => {
            let elegidos: Vec<String> = agrupamiento::k_medias_mas_mas(puntos, k, semilla)
                .into_iter()
                .map(|indice| tabla.nombres[indice].clone())
                .collect();
            (
                LINEA_MAS_MAS,
                format!(
                    "k-means++ (semilla {semilla}) elige centros al azar, prefiriendo filas lejos de los ya elegidos: {}.",
                    algunos(&elegidos, NOMBRES_A_LA_VISTA)
                ),
            )
        }
    };
    momentos.push(Momento {
        fase: "inicio",
        linea: linea_inicial,
        nota: nota_inicial,
        iteracion: 0,
        grupos: None,
        centros: proyectar_centros(&plano, &centros),
        anteriores: None,
        cambiaron: None,
        distorsion: agrupamiento::distorsion_cuadratica(puntos, &centros),
        terminado: false,
    });
    let mut grupos = agrupamiento::de_centros_a_grupos(puntos, &centros);
    let mut distorsion = distorsion_de_los_grupos(puntos, &grupos, &centros);
    let cuentas: Vec<String> = cantidades_por_grupo(&grupos, k)
        .iter()
        .map(|&cantidad| miles(cantidad))
        .collect();
    momentos.push(Momento {
        fase: "asignar",
        linea: LINEA_PRIMERA_ASIGNACION,
        nota: format!(
            "Cada fila va al centro más cercano: quedan grupos de {} filas y la distorsión es {}.",
            unir(&cuentas),
            cientifico(distorsion)
        ),
        iteracion: 0,
        grupos: Some(grupos.clone()),
        centros: proyectar_centros(&plano, &centros),
        anteriores: None,
        cambiaron: None,
        distorsion,
        terminado: false,
    });
    let mut iteracion = 0;
    let mut convergio = false;
    while iteracion < maximo {
        iteracion += 1;
        let anteriores = proyectar_centros(&plano, &centros);
        centros = agrupamiento::de_grupos_a_centros(puntos, &grupos, &centros);
        let antes_de_mover = distorsion;
        distorsion = distorsion_de_los_grupos(puntos, &grupos, &centros);
        momentos.push(Momento {
            fase: "mover",
            linea: LINEA_MOVER,
            nota: format!(
                "Iteración {iteracion}: cada centro se mueve al centro de gravedad de sus filas y la distorsión pasa de {} a {}.",
                cientifico(antes_de_mover),
                cientifico(distorsion)
            ),
            iteracion,
            grupos: Some(grupos.clone()),
            centros: proyectar_centros(&plano, &centros),
            anteriores: Some(anteriores),
            cambiaron: None,
            distorsion,
            terminado: false,
        });
        let nuevos = agrupamiento::de_centros_a_grupos(puntos, &centros);
        let cambiaron = nuevos
            .iter()
            .zip(&grupos)
            .filter(|(nuevo, viejo)| nuevo != viejo)
            .count();
        if cambiaron == 0 {
            convergio = true;
            momentos.push(Momento {
                fase: "fin",
                linea: LINEA_CONVERGIO,
                nota: format!(
                    "Ninguna fila cambió de grupo: Lloyd terminó en {iteracion} {} con distorsión {}.",
                    segun_cantidad(iteracion, "iteración", "iteraciones"),
                    cientifico(distorsion)
                ),
                iteracion,
                grupos: Some(grupos.clone()),
                centros: proyectar_centros(&plano, &centros),
                anteriores: None,
                cambiaron: Some(0),
                distorsion,
                terminado: true,
            });
            break;
        }
        grupos = nuevos;
        let antes_de_asignar = distorsion;
        distorsion = distorsion_de_los_grupos(puntos, &grupos, &centros);
        momentos.push(Momento {
            fase: "asignar",
            linea: LINEA_ASIGNAR,
            nota: format!(
                "Con los centros nuevos, {} {} al centro que ahora tiene más cerca y la distorsión pasa de {} a {}.",
                miles(cambiaron),
                segun_cantidad(cambiaron, "fila se cambia", "filas se cambian"),
                cientifico(antes_de_asignar),
                cientifico(distorsion)
            ),
            iteracion,
            grupos: Some(grupos.clone()),
            centros: proyectar_centros(&plano, &centros),
            anteriores: None,
            cambiaron: Some(cambiaron),
            distorsion,
            terminado: false,
        });
    }
    if !convergio {
        momentos.push(Momento {
            fase: "fin",
            linea: LINEA_MAXIMO,
            nota: format!(
                "Se alcanzó el máximo de {} {} y los grupos todavía cambiaban.",
                miles(maximo),
                segun_cantidad(maximo, "iteración", "iteraciones")
            ),
            iteracion,
            grupos: Some(grupos.clone()),
            centros: proyectar_centros(&plano, &centros),
            anteriores: None,
            cambiaron: None,
            distorsion,
            terminado: true,
        });
    }
    let dibujo: Vec<[f64; 2]> = plano
        .coordenadas
        .iter()
        .map(|lugar| [redondear(lugar[0]), redondear(lugar[1])])
        .collect();
    let puntos_json = json!(dibujo);
    let nombres_json = json!(tabla.nombres);
    let ejes_json = json!(plano.ejes);
    let pasos = muestra_pareja(momentos.len(), MAXIMO_PASOS)
        .into_iter()
        .map(|indice| {
            let momento = &momentos[indice];
            let cantidades = momento
                .grupos
                .as_ref()
                .map(|grupos| cantidades_por_grupo(grupos, k));
            let estado = json!({
                "fase": momento.fase,
                "iteracion": momento.iteracion,
                "k": k,
                "ejes": ejes_json.clone(),
                "nombres": nombres_json.clone(),
                "puntos": puntos_json.clone(),
                "grupos": momento.grupos,
                "centros": momento.centros,
                "centros_anteriores": momento.anteriores,
                "cantidades": cantidades,
                "cambiaron": momento.cambiaron,
                "distorsion": redondear(momento.distorsion),
                "terminado": momento.terminado,
            });
            Paso::nuevo(momento.linea, momento.nota.clone(), estado)
        })
        .collect();
    Ok(pasos)
}

fn codo(entrada: &Entrada) -> Result<Salida, String> {
    let tabla = leer_tabla(entrada, "datos", "p")?;
    exigir_filas(&tabla, 2)?;
    let k_maximo = entrada.natural("k_maximo")?;
    if !(2..=MAXIMO_K_DEL_CODO).contains(&k_maximo) {
        return Err(format!(
            "El k máximo debe estar entre 2 y {MAXIMO_K_DEL_CODO}."
        ));
    }
    let semilla = leer_semilla(entrada)?;
    let hasta = k_maximo.min(tabla.cantidad());
    limitar_cuentas((tabla.cantidad() * tabla.dimension() * hasta * hasta) as f64 * 150.0)?;
    let curva = agrupamiento::curva_del_codo(&tabla.filas, k_maximo, semilla)?;
    let elegido = agrupamiento::k_del_codo(&curva);
    let puntos: Vec<[f64; 2]> = curva
        .iter()
        .enumerate()
        .map(|(indice, &valor)| [(indice + 1) as f64, valor])
        .collect();
    let marcas = match elegido {
        Some(k) => vec![Marca {
            x: k as f64,
            texto: format!("codo: k = {k}"),
        }],
        None => Vec::new(),
    };
    let inicial = curva[0];
    let filas = curva
        .iter()
        .enumerate()
        .map(|(indice, &valor)| {
            let caida = if indice == 0 || curva[indice - 1] <= 0.0 {
                "—".to_string()
            } else {
                porcentaje((curva[indice - 1] - valor) / curva[indice - 1], 0)
            };
            let explicado = if inicial > 0.0 {
                porcentaje(1.0 - valor / inicial, 0)
            } else {
                "—".to_string()
            };
            vec![
                format!("{}", indice + 1),
                cientifico(valor),
                caida,
                explicado,
            ]
        })
        .collect();
    let mut salida = Salida::default();
    salida.agregar(Bloque::linea(
        "Distorsión según la cantidad de grupos",
        "k (cantidad de grupos)",
        "distorsión",
        &puntos,
        marcas,
    ));
    salida.agregar(Bloque::tabla(
        "La curva, número por número",
        &[
            "k",
            "distorsión",
            "caída respecto de k − 1",
            "variación explicada",
        ],
        filas,
    ));
    match elegido {
        Some(k) => {
            let al_llegar = agrupamiento::caida_relativa(curva[k - 2], curva[k - 1]);
            let al_seguir = agrupamiento::caida_relativa(curva[k - 1], curva[k]);
            salida.explicar(format!(
                "El codo está en k = {k}: pasar de {} a {k} grupos quita {} de la distorsión que quedaba, pero pasar a {} solo quita {}. Es el punto donde la caída más se frena: desde ahí, cada grupo extra ayuda mucho menos.",
                k - 1,
                porcentaje(al_llegar, 0),
                k + 1,
                porcentaje(al_seguir, 0)
            ));
            salida.explicar(format!(
                "Con {k} grupos, la distancia de cada fila a su centro explica {} de la variación total (1 − distorsión con k = {k} / distorsión con k = 1).",
                porcentaje(1.0 - curva[k - 1] / inicial, 0)
            ));
            if al_llegar < CAIDA_CLARA {
                salida.explicar(
                    "Aun así, ni siquiera ese paso quita la mitad de lo que quedaba: el codo es débil y los datos quizá no forman grupos bien separados.",
                );
            }
        }
        None if inicial <= 0.0 => {
            salida.explicar(
                "La distorsión es 0 desde el principio: todas las filas son iguales y no hay grupos que buscar.",
            );
        }
        None if curva.len() < 3 => {
            salida.explicar(format!(
                "Con solo {} valores de k no se puede ver un codo: hacen falta al menos 3 para comparar la caída de antes con la de después.",
                curva.len()
            ));
        }
        None => {
            salida.explicar(
                "La caída nunca se frena de golpe: la curva no tiene un codo claro y los datos quizá no forman grupos bien separados.",
            );
        }
    }
    salida.explicar(
        "La distorsión siempre baja al agregar grupos (con un grupo por fila sería 0), así que no elijas el k más grande sino el del codo, donde agregar otro grupo deja de ayudar.",
    );
    let subidas: Vec<String> = curva
        .windows(2)
        .enumerate()
        .filter(|(_, par)| par[1] > par[0] * (1.0 + TOLERANCIA))
        .map(|(indice, _)| format!("{}", indice + 2))
        .collect();
    if !subidas.is_empty() {
        salida.explicar(format!(
            "En k = {} la distorsión subió: k-means quedó en un mínimo local. Cambia la semilla para probar otros comienzos.",
            unir(&subidas)
        ));
    }
    if hasta < k_maximo {
        salida.explicar(format!(
            "La curva llega solo hasta k = {hasta} porque no puede haber más grupos que filas."
        ));
    }
    Ok(salida)
}

fn k_medias_suaves(entrada: &Entrada) -> Result<Salida, String> {
    let tabla = leer_tabla(entrada, "datos", "p")?;
    let k = leer_k(entrada, "k", tabla.cantidad())?;
    let rigidez = entrada.decimal("rigidez")?;
    if rigidez <= 0.0 {
        return Err("La rigidez β debe ser mayor que 0.".to_string());
    }
    let iteraciones = leer_iteraciones(entrada, "iteraciones")?;
    limitar_cuentas((tabla.cantidad() * k * tabla.dimension() * iteraciones) as f64 * 2.0)?;
    let resultado = agrupamiento::k_medias_suaves(&tabla.filas, k, rigidez, iteraciones)?;
    let mut grupos = Vec::with_capacity(tabla.cantidad());
    let mut maximos = Vec::with_capacity(tabla.cantidad());
    for fila in &resultado.responsabilidades {
        let (maximo, donde) = extremo(fila, true);
        grupos.push(donde);
        maximos.push(maximo);
    }
    let mut columnas: Vec<String> = vec![
        "fila".to_string(),
        "grupo más probable".to_string(),
        "responsabilidad máxima".to_string(),
    ];
    let mostrar = k.min(MAXIMO_COLUMNAS_DE_GRUPOS);
    for grupo in 1..=mostrar {
        columnas.push(format!("grupo {grupo}"));
    }
    let filas = resultado
        .responsabilidades
        .iter()
        .enumerate()
        .map(|(indice, fila)| {
            let mut celdas = vec![
                tabla.nombres[indice].clone(),
                format!("{}", grupos[indice] + 1),
                decimal(maximos[indice], 2),
            ];
            celdas.extend(fila.iter().take(mostrar).map(|valor| decimal(*valor, 2)));
            celdas
        })
        .collect();
    let encabezados: Vec<&str> = columnas.iter().map(|columna| columna.as_str()).collect();
    let plano = plano_de_dibujo(&tabla)?;
    let orden = orden_por_grupos(&grupos);
    let elegidos = muestra_pareja(orden.len(), MAXIMO_FILAS_DE_CALOR);
    let mut salida = Salida::default();
    salida.agregar(Bloque::tabla(
        "Responsabilidad de cada grupo sobre cada fila",
        &encabezados,
        filas,
    ));
    salida.agregar(Bloque::calor(
        "Responsabilidades (filas ordenadas por su grupo más probable)",
        elegidos
            .iter()
            .map(|&posicion| tabla.nombres[orden[posicion]].clone())
            .collect(),
        (1..=k).map(|grupo| format!("grupo {grupo}")).collect(),
        elegidos
            .iter()
            .map(|&posicion| resultado.responsabilidades[orden[posicion]].clone())
            .collect(),
    ));
    salida.agregar(dispersion(
        "Cada fila con el color de su grupo más probable",
        &plano,
        &tabla.nombres,
        &grupos,
    ));
    salida.agregar(calor_ordenado(
        "Filas ordenadas por su grupo más probable",
        &tabla,
        &orden,
        Some(&grupos),
    ));
    salida.agregar(calor_de_centros(&tabla, &resultado.centros));
    let ambiguos: Vec<String> = maximos
        .iter()
        .enumerate()
        .filter(|(_, maximo)| **maximo < RESPONSABILIDAD_CLARA)
        .map(|(indice, _)| tabla.nombres[indice].clone())
        .collect();
    let cantidad = tabla.cantidad() as f64;
    let duras = maximos
        .iter()
        .filter(|maximo| **maximo > RESPONSABILIDAD_DURA)
        .count();
    let parejas = maximos
        .iter()
        .filter(|maximo| **maximo < 1.0 / k as f64 + 0.1)
        .count();
    if ambiguos.is_empty() {
        salida.explicar(format!(
            "Todas las filas pertenecen claramente a un grupo: su responsabilidad máxima es de {} o más.",
            decimal(RESPONSABILIDAD_CLARA, 1)
        ));
    } else {
        salida.explicar(format!(
            "{} de {} filas {} una responsabilidad máxima menor que {}: no {} claramente a ningún grupo ({}).",
            miles(ambiguos.len()),
            miles(tabla.cantidad()),
            segun_cantidad(ambiguos.len(), "tiene", "tienen"),
            decimal(RESPONSABILIDAD_CLARA, 1),
            segun_cantidad(ambiguos.len(), "pertenece", "pertenecen"),
            algunos(&ambiguos, NOMBRES_A_LA_VISTA)
        ));
        salida.explicar(
            "Un gen ambiguo puede mezclar dos patrones (por ejemplo, subir y además tener un pico) o responder a más de una señal; k-means duro lo habría metido a la fuerza en un solo grupo.",
        );
    }
    if k > 1 && duras as f64 >= MAYORIA * cantidad {
        salida.explicar(format!(
            "Con β = {} casi todas las responsabilidades son 0 o 1: el resultado es como el de k-means duro. Baja β para ver la incertidumbre.",
            decimal(rigidez, 2)
        ));
    }
    if k > 1 && parejas as f64 >= MAYORIA * cantidad {
        salida.explicar(format!(
            "Con β = {} casi todas las filas se reparten por igual entre los grupos: β es demasiado chico para separarlos. Súbelo.",
            decimal(rigidez, 2)
        ));
    }
    explicar_patrones(&mut salida, &tabla, &grupos, &resultado.centros);
    Ok(salida)
}

fn revisar_matriz_de_distancias(tabla: &Tabla) -> Result<(), String> {
    let cantidad = tabla.cantidad();
    if tabla.dimension() != cantidad {
        return Err(format!(
            "Una matriz de distancias es cuadrada: hay {cantidad} filas y {} columnas.",
            tabla.dimension()
        ));
    }
    for (fila, valores) in tabla.filas.iter().enumerate() {
        if valores[fila].abs() > TOLERANCIA {
            return Err(format!(
                "La distancia de {} consigo mismo debe ser 0 y es {}.",
                tabla.nombres[fila], valores[fila]
            ));
        }
        for (columna, &valor) in valores.iter().enumerate() {
            if valor < 0.0 {
                return Err(format!(
                    "Las distancias no pueden ser negativas (fila {}, columna {}).",
                    tabla.nombres[fila],
                    columna + 1
                ));
            }
            let espejo = tabla.filas[columna][fila];
            if (valor - espejo).abs() > TOLERANCIA * (1.0 + valor.abs()) {
                return Err(format!(
                    "La matriz debe ser simétrica: de {} a {} hay {} pero de vuelta hay {}.",
                    tabla.nombres[fila], tabla.nombres[columna], valor, espejo
                ));
            }
        }
    }
    Ok(())
}

struct Corte {
    grupos: Vec<usize>,
    cuantos: usize,
    abajo: f64,
    arriba: f64,
}

fn cortar_arbol(uniones: &[Union], cantidad: usize, orden: &[usize]) -> Option<Corte> {
    let mut mejor: Option<(usize, f64)> = None;
    for aplicadas in 1..uniones.len() {
        let quedan = cantidad - aplicadas;
        if !(2..=MAXIMO_GRUPOS_DEL_CORTE).contains(&quedan) {
            continue;
        }
        let salto = uniones[aplicadas].distancia - uniones[aplicadas - 1].distancia;
        if salto > mejor.map_or(0.0, |(_, anterior)| anterior) {
            mejor = Some((aplicadas, salto));
        }
    }
    let (aplicadas, _) = mejor?;
    let mut etiqueta: Vec<usize> = (0..cantidad).collect();
    for (paso, union) in uniones[..aplicadas].iter().enumerate() {
        for &miembro in &union.miembros {
            etiqueta[miembro] = cantidad + paso;
        }
    }
    let mut nuevas: Vec<Option<usize>> = vec![None; cantidad + aplicadas];
    let mut grupos = vec![0; cantidad];
    let mut siguiente = 0;
    for &hoja in orden {
        let vieja = etiqueta[hoja];
        let nueva = match nuevas[vieja] {
            Some(nueva) => nueva,
            None => {
                nuevas[vieja] = Some(siguiente);
                siguiente += 1;
                siguiente - 1
            }
        };
        grupos[hoja] = nueva;
    }
    Some(Corte {
        grupos,
        cuantos: siguiente,
        abajo: uniones[aplicadas - 1].distancia,
        arriba: uniones[aplicadas].distancia,
    })
}

fn jerarquico(entrada: &Entrada) -> Result<Salida, String> {
    let tabla = leer_tabla(entrada, "datos", "p")?;
    let eleccion = entrada.opcion("distancia", MEDIDAS)?;
    exigir_filas(&tabla, 2)?;
    limitar_matriz(&tabla)?;
    let es_matriz = eleccion == MEDIDA_MATRIZ;
    let distancias = if es_matriz {
        revisar_matriz_de_distancias(&tabla)?;
        tabla.filas.clone()
    } else if eleccion == MEDIDA_PEARSON {
        agrupamiento::matriz_de_distancias(&tabla.filas, Medida::Pearson)
    } else {
        agrupamiento::matriz_de_distancias(&tabla.filas, Medida::Euclidea)
    };
    let uniones = agrupamiento::agrupamiento_jerarquico(&distancias)?;
    let newick = agrupamiento::arbol_newick(&uniones, &tabla.nombres);
    let cantidad = tabla.cantidad();
    let orden: Vec<usize> = uniones
        .last()
        .map(|union| union.miembros.clone())
        .unwrap_or_else(|| (0..cantidad).collect());
    let filas = uniones
        .iter()
        .enumerate()
        .map(|(paso, union)| {
            let nombres: Vec<String> = union
                .miembros
                .iter()
                .map(|&miembro| tabla.nombres[miembro].clone())
                .collect();
            vec![
                format!("{}", paso + 1),
                algunos(&nombres, NOMBRES_POR_GRUPO),
                miles(union.miembros.len()),
                decimal(union.distancia, 3),
            ]
        })
        .collect();
    let mut salida = Salida::default();
    salida.agregar(Bloque::arbol("Árbol de uniones", newick));
    salida.agregar(Bloque::tabla(
        "Uniones en el orden en que ocurrieron",
        &["paso", "grupo nuevo", "elementos", "distancia promedio"],
        filas,
    ));
    salida.agregar(Bloque::lista(
        "Elementos de cada grupo nuevo (número de fila, contando desde 1)",
        uniones.iter().map(|union| {
            union
                .miembros
                .iter()
                .map(|miembro| (miembro + 1).to_string())
                .collect::<Vec<_>>()
                .join(" ")
        }),
    ));
    let corte = cortar_arbol(&uniones, cantidad, &orden);
    if es_matriz {
        let elegidos: Vec<usize> = muestra_pareja(orden.len(), MAXIMO_FILAS_DE_CALOR)
            .into_iter()
            .map(|posicion| orden[posicion])
            .collect();
        salida.agregar(Bloque::calor(
            "Matriz de distancias en el orden del árbol",
            elegidos
                .iter()
                .map(|&fila| tabla.nombres[fila].clone())
                .collect(),
            elegidos
                .iter()
                .map(|&fila| tabla.nombres[fila].clone())
                .collect(),
            elegidos
                .iter()
                .map(|&fila| {
                    elegidos
                        .iter()
                        .map(|&columna| distancias[fila][columna])
                        .collect()
                })
                .collect(),
        ));
    } else {
        let grupos_para_color = corte
            .as_ref()
            .map(|elegido| elegido.grupos.clone())
            .unwrap_or_else(|| vec![0; cantidad]);
        salida.agregar(calor_ordenado(
            "Filas en el orden de las hojas del árbol",
            &tabla,
            &orden,
            corte.as_ref().map(|elegido| elegido.grupos.as_slice()),
        ));
        let plano = plano_de_dibujo(&tabla)?;
        salida.agregar(dispersion(
            "Cada fila con el color del grupo que deja el corte sugerido",
            &plano,
            &tabla.nombres,
            &grupos_para_color,
        ));
    }
    let primera = &uniones[0];
    let (a, b) = (primera.miembros[0], primera.miembros[1]);
    if eleccion == MEDIDA_PEARSON {
        salida.explicar(format!(
            "Los primeros en unirse son {} y {}, a distancia {} (correlación r = {}): son los perfiles con la forma más parecida.",
            tabla.nombres[a],
            tabla.nombres[b],
            decimal(primera.distancia, 3),
            decimal(1.0 - primera.distancia, 2)
        ));
        salida.explicar(
            "Con Pearson dos genes se unen pronto si suben y bajan a la vez, aunque uno cambie mucho más que el otro; con la euclídea también importa cuánto cambian.",
        );
    } else {
        salida.explicar(format!(
            "Los primeros en unirse son {} y {}, a distancia {}: son el par más cercano de toda la tabla.",
            tabla.nombres[a],
            tabla.nombres[b],
            decimal(primera.distancia, 3)
        ));
    }
    salida.explicar(
        "La altura de cada unión es la distancia promedio entre todos los pares de elementos de los dos grupos; en el árbol, la suma de ramas entre dos hojas es esa distancia.",
    );
    match &corte {
        Some(elegido) => {
            salida.explicar(format!(
                "El mayor salto entre una unión y la siguiente va de {} a {}: si cortas el árbol ahí quedan {} grupos bien separados.",
                decimal(elegido.abajo, 3),
                decimal(elegido.arriba, 3),
                elegido.cuantos
            ));
            if !es_matriz {
                let centros = centros_de(&tabla, &elegido.grupos, elegido.cuantos);
                explicar_patrones(&mut salida, &tabla, &elegido.grupos, &centros);
            }
        }
        None => {
            if cantidad >= 3 {
                salida.explicar(
                    "Las alturas de las uniones crecen de a poco, sin un salto claro: el árbol no sugiere una cantidad natural de grupos.",
                );
            }
        }
    }
    Ok(salida)
}

fn cast(entrada: &Entrada) -> Result<Salida, String> {
    let tabla = leer_tabla(entrada, "datos", "p")?;
    let umbral = entrada.decimal("umbral")?;
    if umbral <= 0.0 || umbral > 1.0 {
        return Err(
            "El umbral θ es una correlación: debe estar entre 0 (sin incluir) y 1.".to_string(),
        );
    }
    exigir_filas(&tabla, 2)?;
    limitar_matriz(&tabla)?;
    let similitudes = agrupamiento::matriz_de_correlaciones(&tabla.filas);
    let encontrados = agrupamiento::cast(&similitudes, umbral)?;
    let cuantos = encontrados.len();
    let mut grupos = vec![0; tabla.cantidad()];
    for (grupo, miembros) in encontrados.iter().enumerate() {
        for &miembro in miembros {
            grupos[miembro] = grupo;
        }
    }
    let centros = centros_de(&tabla, &grupos, cuantos);
    let con_tiempos = tabla.dimension() >= MINIMO_TIEMPOS;
    let filas = encontrados
        .iter()
        .enumerate()
        .map(|(grupo, miembros)| {
            let nombres: Vec<String> = miembros
                .iter()
                .map(|&miembro| tabla.nombres[miembro].clone())
                .collect();
            vec![
                format!("{}", grupo + 1),
                miles(miembros.len()),
                similitud_interna(&similitudes, miembros),
                if con_tiempos {
                    describir_perfil(&centros[grupo], &tabla.columnas)
                } else {
                    "—".to_string()
                },
                algunos(&nombres, NOMBRES_POR_GRUPO),
            ]
        })
        .collect();
    let plano = plano_de_dibujo(&tabla)?;
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Grupos encontrados", miles(cuantos)));
    salida.agregar(Bloque::tabla(
        "Grupos en el orden en que CAST los armó",
        &[
            "grupo",
            "filas",
            "correlación promedio dentro del grupo",
            "perfil promedio",
            "miembros",
        ],
        filas,
    ));
    salida.agregar(dispersion(
        "Cada fila con el color de su grupo",
        &plano,
        &tabla.nombres,
        &grupos,
    ));
    salida.agregar(calor_ordenado(
        "Filas ordenadas por grupo",
        &tabla,
        &orden_por_grupos(&grupos),
        Some(&grupos),
    ));
    salida.explicar(format!(
        "CAST no necesita que le digas k: armó {cuantos} {} en {} cada gen tiene, en promedio, una correlación de al menos {} con su grupo.",
        segun_cantidad(cuantos, "grupo", "grupos"),
        segun_cantidad(cuantos, "el que", "los que"),
        decimal(umbral, 2)
    ));
    let solos: Vec<String> = encontrados
        .iter()
        .filter(|miembros| miembros.len() == 1)
        .map(|miembros| tabla.nombres[miembros[0]].clone())
        .collect();
    if cuantos == tabla.cantidad() {
        salida.explicar(
            "Cada fila quedó sola: θ es tan alto que ningún par se parece lo suficiente. Bájalo.",
        );
    } else if cuantos == 1 {
        salida.explicar(
            "Todas las filas quedaron en un solo grupo: θ es tan bajo que todo se parece a todo. Súbelo.",
        );
    } else if !solos.is_empty() {
        salida.explicar(format!(
            "{} {} aparte: no se parece lo suficiente a ningún grupo. Puede ser un gen con un patrón propio o con mucho ruido.",
            algunos(&solos, NOMBRES_A_LA_VISTA),
            segun_cantidad(solos.len(), "quedó", "quedaron")
        ));
    }
    if tabla.dimension() <= 2 {
        salida.explicar(
            "Con 2 columnas o menos la correlación de Pearson da casi siempre +1 o −1: CAST necesita al menos 3 condiciones para agrupar con sentido.",
        );
    }
    if cuantos > 1 && cuantos < tabla.cantidad() {
        explicar_patrones(&mut salida, &tabla, &grupos, &centros);
    }
    Ok(salida)
}

fn similitud_interna(similitudes: &[Vec<f64>], miembros: &[usize]) -> String {
    if miembros.len() < 2 {
        return "—".to_string();
    }
    let mut suma = 0.0;
    let mut pares = 0usize;
    for (posicion, &uno) in miembros.iter().enumerate() {
        for &otro in &miembros[posicion + 1..] {
            suma += similitudes[uno][otro];
            pares += 1;
        }
    }
    decimal(suma / pares as f64, 2)
}

fn componentes_principales(entrada: &Entrada) -> Result<Salida, String> {
    let tabla = leer_tabla(entrada, "datos", "p")?;
    exigir_filas(&tabla, 2)?;
    if tabla.dimension() < 2 {
        return Err("La tabla tiene una sola columna: no hay nada que proyectar.".to_string());
    }
    let componentes = agrupamiento::componentes_principales(&tabla.filas)?;
    let total = componentes.varianza_total;
    let parte_primera = if total > 0.0 {
        componentes.varianza_primera / total
    } else {
        0.0
    };
    let parte_segunda = if total > 0.0 {
        componentes.varianza_segunda / total
    } else {
        0.0
    };
    let mut salida = Salida::default();
    salida.agregar(dispersion(
        "Cada fila proyectada sobre las dos componentes principales",
        &Plano {
            coordenadas: componentes.coordenadas.clone(),
            ejes: [
                format!("Componente principal 1 ({})", porcentaje(parte_primera, 0)),
                format!("Componente principal 2 ({})", porcentaje(parte_segunda, 0)),
            ],
            componentes: None,
        },
        &tabla.nombres,
        &vec![0; tabla.cantidad()],
    ));
    salida.agregar(Bloque::valor(
        "Varianza explicada por las dos componentes",
        porcentaje(parte_primera + parte_segunda, 1),
    ));
    salida.agregar(Bloque::barras(
        "Peso de cada columna en la componente 1",
        "columna",
        "peso",
        tabla.columnas.clone(),
        componentes.primera.clone(),
    ));
    salida.agregar(Bloque::barras(
        "Peso de cada columna en la componente 2",
        "columna",
        "peso",
        tabla.columnas.clone(),
        componentes.segunda.clone(),
    ));
    salida.agregar(Bloque::tabla(
        "Varianza de cada componente",
        &["componente", "varianza", "parte del total"],
        vec![
            vec![
                "1".to_string(),
                cientifico(componentes.varianza_primera),
                porcentaje(parte_primera, 1),
            ],
            vec![
                "2".to_string(),
                cientifico(componentes.varianza_segunda),
                porcentaje(parte_segunda, 1),
            ],
            vec!["todas".to_string(), cientifico(total), porcentaje(1.0, 0)],
        ],
    ));
    let juntas = parte_primera + parte_segunda;
    if total <= 0.0 {
        salida.explicar("Todas las filas son iguales: no hay variación que mostrar.");
        return Ok(salida);
    }
    salida.explicar(format!(
        "Las dos primeras componentes juntan {} de la variación total ({} la primera, {} la segunda).",
        porcentaje(juntas, 0),
        porcentaje(parte_primera, 0),
        porcentaje(parte_segunda, 0)
    ));
    if juntas >= VARIANZA_FIEL {
        salida.explicar(
            "Con 70 % o más, el dibujo en 2D es fiel: filas cercanas en el plano tienen perfiles parecidos.",
        );
    } else if juntas < VARIANZA_POBRE {
        salida.explicar(
            "Con menos de 50 %, buena parte de la variación queda fuera del dibujo: dos filas pueden verse juntas en el plano y ser distintas en las otras direcciones.",
        );
    }
    let (_, donde) = extremo(
        &componentes
            .primera
            .iter()
            .map(|peso| peso.abs())
            .collect::<Vec<_>>(),
        true,
    );
    salida.explicar(format!(
        "La columna que más pesa en la componente 1 es {}: las filas se separan de izquierda a derecha sobre todo por su valor ahí.",
        tabla.columnas[donde]
    ));
    if tabla.dimension() >= MINIMO_TIEMPOS {
        let tiempos: Vec<f64> = (0..tabla.dimension()).map(|tiempo| tiempo as f64).collect();
        let tendencia = agrupamiento::correlacion_pearson(&componentes.primera, &tiempos);
        if tendencia.abs() >= CORRELACION_FUERTE {
            salida.explicar(
                "Los pesos de la componente 1 crecen (o decrecen) con el orden de las columnas: si son tiempos, esa componente separa a los genes que suben de los que bajan.",
            );
        }
    }
    if tabla.dimension() == 2 {
        salida
            .explicar("Tus datos ya tienen 2 columnas: la proyección solo los gira y los centra.");
    }
    Ok(salida)
}

#[cfg(test)]
mod pruebas {
    use super::*;

    fn entrada_de_lloyd(datos: &str, k: &str) -> Entrada {
        Entrada::nueva(&[
            ("datos", datos),
            ("k", k),
            ("inicio", PRIMEROS_PUNTOS),
            ("semilla", "1"),
            ("iteraciones", "100"),
        ])
    }

    fn distorsiones(pasos: &[Paso]) -> Vec<f64> {
        pasos
            .iter()
            .map(|paso| paso.estado["distorsion"].as_f64().unwrap())
            .collect()
    }

    #[test]
    fn el_ejemplo_de_lloyd_necesita_varias_vueltas() {
        let pasos = pasos_lloyd(&entrada_de_lloyd(TABLA_EJEMPLO, "3")).unwrap();
        let ultimo = pasos.last().unwrap();
        assert!(ultimo.estado["iteracion"].as_u64().unwrap() >= 4);
        assert_eq!(ultimo.estado["fase"], "fin");
        assert_eq!(ultimo.estado["terminado"], true);
        assert_eq!(ultimo.estado["cambiaron"], 0);
    }

    #[test]
    fn la_distorsion_baja_en_cada_medio_paso() {
        for datos in [TABLA_EJEMPLO, NUBES_EJEMPLO] {
            let pasos = pasos_lloyd(&entrada_de_lloyd(datos, "3")).unwrap();
            let valores = distorsiones(&pasos);
            assert!(valores.windows(2).all(|par| par[1] <= par[0] + 1e-4));
            let bajan_al_mover = pasos.windows(2).any(|par| {
                par[1].estado["fase"] == "mover"
                    && par[1].estado["distorsion"].as_f64() < par[0].estado["distorsion"].as_f64()
            });
            let bajan_al_asignar = pasos.windows(2).any(|par| {
                par[1].estado["fase"] == "asignar"
                    && par[1].estado["distorsion"].as_f64() < par[0].estado["distorsion"].as_f64()
            });
            assert!(bajan_al_mover && bajan_al_asignar);
        }
    }

    #[test]
    fn cada_estado_trae_todo_para_dibujarse_solo() {
        let pasos = pasos_lloyd(&entrada_de_lloyd(NUBES_EJEMPLO, "3")).unwrap();
        let claves = [
            "fase",
            "iteracion",
            "k",
            "ejes",
            "nombres",
            "puntos",
            "grupos",
            "centros",
            "centros_anteriores",
            "cantidades",
            "cambiaron",
            "distorsion",
            "terminado",
        ];
        for paso in &pasos {
            let estado = paso.estado.as_object().unwrap();
            assert_eq!(estado.len(), claves.len());
            for clave in claves {
                assert!(estado.contains_key(clave), "falta {clave}");
            }
            assert_eq!(estado["puntos"].as_array().unwrap().len(), 24);
            assert_eq!(estado["centros"].as_array().unwrap().len(), 3);
            assert_eq!(estado["puntos"][0], json!([2.8, 1.4]));
        }
        assert_eq!(pasos[0].estado["fase"], "inicio");
        assert!(pasos[0].estado["grupos"].is_null());
        assert_eq!(pasos[1].estado["fase"], "asignar");
        assert_eq!(pasos[2].estado["fase"], "mover");
        assert!(pasos[2].estado["centros_anteriores"].is_array());
    }

    #[test]
    fn el_codo_del_ejemplo_esta_en_tres() {
        let salida = codo(&Entrada::nueva(&[
            ("datos", TABLA_EJEMPLO),
            ("k_maximo", "10"),
            ("semilla", "1"),
        ]))
        .unwrap();
        let Bloque::Linea { marcas, .. } = &salida.bloques[0] else {
            panic!("el primer bloque no es una línea");
        };
        assert_eq!(marcas[0].x, 3.0);
        assert!(salida.interpretacion[0].starts_with("El codo está en k = 3"));
    }

    #[test]
    fn el_ejemplo_muestra_los_tres_patrones() {
        let salida = lloyd(&entrada_de_lloyd(TABLA_EJEMPLO, "3")).unwrap();
        let texto = salida.interpretacion.join(" ");
        assert!(texto.contains(" sube "));
        assert!(texto.contains(" baja "));
        assert!(texto.contains("pico pasajero"));
        let arbol = jerarquico(&Entrada::nueva(&[
            ("datos", TABLA_EJEMPLO),
            ("distancia", MEDIDA_PEARSON),
        ]))
        .unwrap();
        assert!(arbol.interpretacion.join(" ").contains("quedan 3 grupos"));
    }

    #[test]
    fn una_corrida_larga_se_reparte_en_400_pasos() {
        let elegidos = muestra_pareja(2003, MAXIMO_PASOS);
        assert_eq!(elegidos.len(), MAXIMO_PASOS);
        assert_eq!(elegidos[0], 0);
        assert_eq!(elegidos[MAXIMO_PASOS - 1], 2002);
        assert!(elegidos.windows(2).all(|par| par[0] < par[1]));
        assert_eq!(muestra_pareja(7, MAXIMO_PASOS), (0..7).collect::<Vec<_>>());
    }

    #[test]
    fn el_paso_a_paso_rechaza_tablas_demasiado_grandes() {
        let filas: Vec<String> = (0..MAXIMO_ANIMABLE + 1)
            .map(|fila| format!("{fila} {}", fila % 7))
            .collect();
        let error = pasos_lloyd(&entrada_de_lloyd(&filas.join("\n"), "3")).unwrap_err();
        assert!(error.contains("hasta 1\u{a0}000 filas"), "{error}");
    }

    #[test]
    fn los_valores_enormes_se_rechazan() {
        let error = lloyd(&entrada_de_lloyd("a 1 2\nb 1e300 3\nc 4 5", "2")).unwrap_err();
        assert!(error.contains("10¹²"), "{error}");
        assert!(lloyd(&entrada_de_lloyd("a 1 2\nb 1e11 3\nc 4 5", "2")).is_ok());
    }

    #[test]
    fn una_sola_iteracion_se_dice_en_singular() {
        let salida = lloyd(&entrada_de_lloyd("0\n10\n0.1\n10.1", "2")).unwrap();
        let texto = salida.interpretacion.join(" ");
        assert!(texto.contains("tras 1 iteración "), "{texto}");
        let uno = lloyd(&entrada_de_lloyd("0\n10\n0.1\n10.1", "1")).unwrap();
        assert!(uno.interpretacion[0].contains("en 1 grupo de"));
    }
}

const VALORES_A_LA_VISTA: usize = 6;
const FILAS_A_LA_VISTA: usize = 8;
const ELEMENTOS_A_LA_VISTA: usize = 24;
const PASOS_DETALLADOS: usize = 3;
const PASOS_DESEADOS: usize = 12;
const INTENTOS_DEL_CODO: u64 = 3;
const ITERACIONES_DEL_CODO: usize = 100;
const VUELTAS_DE_POTENCIA: usize = 1000;
const CAMBIO_DESPRECIABLE: f64 = 1e-24;
const PARTE_DESPRECIABLE: f64 = 1e-12;
const SEMILLA_DE_POTENCIA: u64 = 1;
const YA_ES_CENTRO: f64 = -1.0;
const EVENTOS_DEL_PRIMER_GRUPO: usize = 10;
const EVENTOS_DE_LOS_SIGUIENTES: usize = 4;

fn numero(valor: f64, cifras: usize) -> String {
    let texto = decimal(valor, cifras);
    match texto.strip_prefix('-') {
        Some(resto) if resto.chars().all(|letra| letra == '0' || letra == '.') => resto.to_string(),
        _ => texto,
    }
}

fn vector_corto(valores: &[f64]) -> String {
    let mostrados: Vec<String> = valores
        .iter()
        .take(VALORES_A_LA_VISTA)
        .map(|valor| numero(*valor, 2))
        .collect();
    if valores.len() > VALORES_A_LA_VISTA {
        format!("[{}, …]", mostrados.join(", "))
    } else {
        format!("[{}]", mostrados.join(", "))
    }
}

fn texto_de_distancia(valor: f64) -> String {
    if valor.is_infinite() {
        "∞".to_string()
    } else if valor < 0.0 {
        "centro".to_string()
    } else {
        numero(valor, 2)
    }
}

fn a_la_vista(total: usize, maximo: usize, destacados: &[usize]) -> Vec<usize> {
    let mut vista: Vec<usize> = (0..total.min(maximo)).collect();
    for &destacado in destacados {
        if destacado < total && !vista.contains(&destacado) {
            vista.push(destacado);
        }
    }
    vista.sort_unstable();
    vista
}

fn posiciones(vista: &[usize], destacados: &[usize]) -> Vec<usize> {
    destacados
        .iter()
        .filter_map(|destacado| vista.iter().position(|visto| visto == destacado))
        .collect()
}

fn salto_de(total: usize) -> usize {
    (total / PASOS_DESEADOS).max(1)
}

fn toca_narrar(indice: usize, salto: usize) -> bool {
    indice < PASOS_DETALLADOS || indice.is_multiple_of(salto)
}

fn lista_de(pizarra: Pizarra, titulo: &str, textos: &[String], destacados: &[usize]) -> Pizarra {
    let vista = a_la_vista(textos.len(), ELEMENTOS_A_LA_VISTA, destacados);
    let titulo = if vista.len() < textos.len() {
        format!("{titulo} (una parte)")
    } else {
        titulo.to_string()
    };
    let elementos = vista.iter().map(|&indice| textos[indice].clone()).collect();
    pizarra.lista(&titulo, elementos, &posiciones(&vista, destacados))
}

fn tabla_de_numeros(
    pizarra: Pizarra,
    titulo: &str,
    nombres: &[String],
    columnas: &[String],
    valores: &[Vec<f64>],
    cifras: usize,
    resaltar: &[(usize, usize)],
) -> Pizarra {
    let filas_pedidas: Vec<usize> = resaltar.iter().map(|marca| marca.0).collect();
    let columnas_pedidas: Vec<usize> = resaltar.iter().map(|marca| marca.1).collect();
    let filas_vistas = a_la_vista(valores.len(), FILAS_A_LA_VISTA, &filas_pedidas);
    let columnas_vistas = a_la_vista(columnas.len(), VALORES_A_LA_VISTA, &columnas_pedidas);
    let titulo = if filas_vistas.len() < valores.len() || columnas_vistas.len() < columnas.len() {
        format!("{titulo} (una parte)")
    } else {
        titulo.to_string()
    };
    let mut encabezados = vec!["fila".to_string()];
    encabezados.extend(
        columnas_vistas
            .iter()
            .map(|&columna| columnas[columna].clone()),
    );
    let referencias: Vec<&str> = encabezados.iter().map(String::as_str).collect();
    let cuerpo = filas_vistas
        .iter()
        .map(|&fila| {
            let mut celdas = vec![nombres[fila].clone()];
            celdas.extend(
                columnas_vistas
                    .iter()
                    .map(|&columna| numero(valores[fila][columna], cifras)),
            );
            celdas
        })
        .collect();
    let marcas: Vec<(usize, usize)> = resaltar
        .iter()
        .filter_map(|&(fila, columna)| {
            let dentro_fila = filas_vistas.iter().position(|&visto| visto == fila)?;
            let dentro_columna = columnas_vistas.iter().position(|&visto| visto == columna)?;
            Some((dentro_fila, dentro_columna + 1))
        })
        .collect();
    pizarra.tabla(&titulo, &referencias, cuerpo, &marcas)
}

fn tabla_de_textos(
    pizarra: Pizarra,
    titulo: &str,
    columnas: &[&str],
    filas: &[Vec<String>],
    destacada: Option<usize>,
) -> Pizarra {
    let pedidas: Vec<usize> = destacada.into_iter().collect();
    let vista = a_la_vista(filas.len(), FILAS_A_LA_VISTA, &pedidas);
    let titulo = if vista.len() < filas.len() {
        format!("{titulo} (una parte)")
    } else {
        titulo.to_string()
    };
    let cuerpo = vista.iter().map(|&fila| filas[fila].clone()).collect();
    let marcas: Vec<(usize, usize)> = posiciones(&vista, &pedidas)
        .into_iter()
        .flat_map(|fila| (0..columnas.len()).map(move |columna| (fila, columna)))
        .collect();
    pizarra.tabla(&titulo, columnas, cuerpo, &marcas)
}

fn etiqueta_de_grupo(tabla: &Tabla, miembros: &[usize]) -> String {
    let nombres: Vec<String> = miembros
        .iter()
        .map(|&miembro| tabla.nombres[miembro].clone())
        .collect();
    if nombres.len() == 1 {
        nombres[0].clone()
    } else {
        format!("{{{}}}", algunos(&nombres, 4))
    }
}

fn nombres_de(tabla: &Tabla, indices: &[usize]) -> String {
    let nombres: Vec<String> = indices
        .iter()
        .map(|&indice| tabla.nombres[indice].clone())
        .collect();
    algunos(&nombres, NOMBRES_POR_GRUPO)
}

fn pasos_distancias(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let tabla = leer_tabla(entrada, "datos", "p")?;
    exigir_filas(&tabla, 2)?;
    limitar_matriz(&tabla)?;
    Ok(narrar_distancias(&tabla).0)
}

fn narrar_distancias(tabla: &Tabla) -> (Vec<Paso>, [f64; 2]) {
    let puntos = &tabla.filas;
    let cantidad = puntos.len();
    let dimension = tabla.dimension();
    let pares = cantidad * (cantidad - 1) / 2;
    let nombre_a = &tabla.nombres[0];
    let nombre_b = &tabla.nombres[1];
    let (a, b) = (&puntos[0], &puntos[1]);
    let euclidea = agrupamiento::matriz_de_distancias(puntos, Medida::Euclidea);
    let pearson = agrupamiento::matriz_de_distancias(puntos, Medida::Pearson);
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        tabla_de_numeros(
            Pizarra::nueva().variable("cantidad", cantidad),
            "Tus datos: una fila por elemento",
            &tabla.nombres,
            &tabla.columnas,
            puntos,
            2,
            &[],
        )
        .paso(
            "let cantidad = puntos.len();",
            format!(
                "Tienes {} filas con {} números cada una, así que hay {} pares de filas. Vamos a medir qué tan separado está cada par de dos maneras: la distancia euclídea (cuánto se alejan los números) y la de Pearson (cuánto se parecen las formas).",
                miles(cantidad),
                miles(dimension),
                miles(pares)
            ),
        ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("cantidad", cantidad)
            .variable("fila", 0)
            .variable("columna", 1)
            .paso(
                "for columna in fila + 1..cantidad {",
                format!(
                    "Empezamos con la euclídea. Cada par se calcula una sola vez (columna siempre mayor que fila). El primer par es {nombre_a} contra {nombre_b}."
                ),
            ),
    );
    let mut suma = 0.0;
    let mut cuentas: Vec<Vec<String>> = Vec::new();
    for (indice, (x, y)) in a.iter().zip(b).enumerate() {
        let diferencia = x - y;
        suma += diferencia * diferencia;
        cuentas.push(vec![
            tabla.columnas[indice].clone(),
            numero(*x, 2),
            numero(*y, 2),
            numero(diferencia, 2),
            numero(diferencia * diferencia, 3),
        ]);
        if indice < VALORES_A_LA_VISTA {
            let nota = if indice == 0 {
                format!(
                    "En la primera columna {nombre_a} vale {} y {nombre_b} vale {}. Restamos y elevamos al cuadrado ({}) para que el signo no importe, y lo sumamos: suma empieza en 0 y pasa a {}.",
                    numero(*x, 2),
                    numero(*y, 2),
                    numero(diferencia * diferencia, 3),
                    numero(suma, 3)
                )
            } else {
                format!(
                    "Columna {}: diferencia {}, al cuadrado {}. La suma pasa a {}.",
                    tabla.columnas[indice],
                    numero(diferencia, 2),
                    numero(diferencia * diferencia, 3),
                    numero(suma, 3)
                )
            };
            narracion.agregar(
                tabla_de_textos(
                    Pizarra::nueva()
                        .variable("suma", numero(suma, 3))
                        .variable("x", numero(*x, 2))
                        .variable("y", numero(*y, 2))
                        .variable("diferencia", numero(diferencia, 2)),
                    "Cuenta de este par, columna por columna",
                    &[
                        "columna",
                        nombre_a.as_str(),
                        nombre_b.as_str(),
                        "diferencia",
                        "al cuadrado",
                    ],
                    &cuentas,
                    Some(indice),
                )
                .paso("suma += diferencia * diferencia;", nota),
            );
        }
    }
    if dimension > VALORES_A_LA_VISTA {
        narracion.agregar(Pizarra::nueva().variable("suma", numero(suma, 3)).paso(
            "suma += diferencia * diferencia;",
            format!(
                "Las otras {} columnas hacen exactamente lo mismo. Al terminar, la suma vale {}.",
                dimension - VALORES_A_LA_VISTA,
                numero(suma, 3)
            ),
        ));
    }
    let primera_distancia = suma.sqrt();
    narracion.agregar(
        Pizarra::nueva()
            .variable("suma", numero(suma, 3))
            .variable("valores[columna]", numero(primera_distancia, 3))
            .paso(
                "distancia_al_cuadrado(a, b).sqrt()",
                format!(
                    "Sacamos la raíz cuadrada de la suma: {} es la distancia en línea recta entre {nombre_a} y {nombre_b}, como el teorema de Pitágoras pero con {} lados.",
                    numero(primera_distancia, 3),
                    miles(dimension)
                ),
            ),
    );
    let mut vistas = vec![format!(
        "{nombre_a} ↔ {nombre_b}: {}",
        numero(primera_distancia, 3)
    )];
    for (columna, otro) in puntos.iter().enumerate().take(5).skip(2) {
        let lejania = agrupamiento::distancia_euclidea(a, otro);
        vistas.push(format!(
            "{nombre_a} ↔ {}: {}",
            tabla.nombres[columna],
            numero(lejania, 3)
        ));
        narracion.agregar(
            lista_de(
                Pizarra::nueva()
                    .variable("cantidad", cantidad)
                    .variable("fila", 0)
                    .variable("columna", columna)
                    .variable("valores[columna]", numero(lejania, 3)),
                "Distancias euclídeas de la primera fila a las demás",
                &vistas,
                &[vistas.len() - 1],
            )
            .paso(
                "valores[columna] = distancia_euclidea(&puntos[fila], &puntos[columna]);",
                format!(
                    "Mismo cálculo con {}: distancia {}.",
                    tabla.nombres[columna],
                    numero(lejania, 3)
                ),
            ),
        );
    }
    if pares > vistas.len() {
        narracion.agregar(
            Pizarra::nueva()
                .variable("cantidad", cantidad)
                .paso(
                    ".into_par_iter()",
                    format!(
                        "Los otros {} pares se calculan igual. Para no esperar, el programa reparte las filas entre todos los núcleos de tu computadora y los hace a la vez.",
                        miles(pares - vistas.len())
                    ),
                ),
        );
    }
    narracion.agregar(Pizarra::nueva().variable("cantidad", cantidad).paso(
        "fila[diagonal] = 0.0;",
        "La distancia de cada fila consigo misma es 0, así que la diagonal se rellena con ceros.",
    ));
    narracion.agregar(
        tabla_de_numeros(
            Pizarra::nueva().variable("cantidad", cantidad),
            "Distancia euclídea entre cada par",
            &tabla.nombres,
            &tabla.nombres,
            &euclidea,
            2,
            &[(0, 1), (1, 0)],
        )
        .paso(
            "copiar_al_otro_lado(&mut distancias);",
            format!(
                "Solo calculamos el triángulo de arriba. Como la distancia de {nombre_a} a {nombre_b} es la misma que de {nombre_b} a {nombre_a}, se copia al triángulo de abajo y la matriz queda completa."
            ),
        ),
    );
    narracion.agregar(
        Pizarra::nueva().variable("cantidad", cantidad).paso(
            "Medida::Pearson => matriz_de_correlaciones(puntos)",
            "Ahora la distancia de Pearson. Ella ignora cuánto valen los números y mira solo su forma: si dos filas suben y bajan juntas, están cerca aunque una valga mucho más que la otra. Para eso primero se «estandariza» cada fila.",
        ),
    );
    let promedio = a.iter().sum::<f64>() / dimension as f64;
    let centrados: Vec<f64> = a.iter().map(|valor| valor - promedio).collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("promedio", numero(promedio, 3))
            .paso(
                "let promedio = valores.iter().sum::<f64>() / valores.len() as f64;",
                format!(
                    "Promedio de {nombre_a}: {}. Es el nivel medio de la fila.",
                    numero(promedio, 3)
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("promedio", numero(promedio, 3))
            .variable("centrados", vector_corto(&centrados))
            .paso(
                "let mut centrados: Vec<f64> = valores.iter().map(|valor| valor - promedio).collect();",
                "Le restamos el promedio a cada número: ahora la fila se mueve alrededor de 0 y solo importa cuánto sube o baja respecto de su nivel medio.",
            ),
    );
    let magnitud = a.iter().fold(0.0f64, |mayor, valor| mayor.max(valor.abs()));
    let desvio = centrados
        .iter()
        .fold(0.0f64, |mayor, valor| mayor.max(valor.abs()));
    let estandar_a = agrupamiento::estandarizar(a);
    let estandar_b = agrupamiento::estandarizar(b);
    if desvio <= magnitud * PARTE_DESPRECIABLE {
        narracion.agregar(
            Pizarra::nueva()
                .variable("desvio", numero(desvio, 3))
                .paso(
                    "return vec![0.0; valores.len()];",
                    format!(
                        "{nombre_a} no cambia en ninguna columna, así que no tiene forma que comparar: su fila estandarizada queda en ceros y su correlación con las demás valdrá 0."
                    ),
                ),
        );
    } else {
        let escalados: Vec<f64> = centrados.iter().map(|valor| valor / desvio).collect();
        let largo = agrupamiento::producto_punto(&escalados, &escalados).sqrt();
        narracion.agregar(
            Pizarra::nueva()
                .variable("desvio", numero(desvio, 3))
                .paso(
                    "let desvio = centrados",
                    format!(
                        "desvio es el mayor valor absoluto de los centrados ({}). Dividimos entre él solo para evitar números enormes; no cambia la forma.",
                        numero(desvio, 3)
                    ),
                ),
        );
        narracion.agregar(
            Pizarra::nueva()
                .variable("desvio", numero(desvio, 3))
                .variable("largo", numero(largo, 3))
                .paso(
                    "let largo = producto_punto(&centrados, &centrados).sqrt();",
                    format!(
                        "El largo de la fila es la raíz de la suma de sus números al cuadrado: {}.",
                        numero(largo, 3)
                    ),
                ),
        );
        narracion.agregar(
            Pizarra::nueva()
                .variable("largo", numero(largo, 3))
                .variable("valores", vector_corto(&estandar_a))
                .paso(
                    "*valor /= largo;",
                    format!(
                        "Dividimos entre el largo: {nombre_a} queda con promedio 0 y largo 1. Ya solo conserva su forma, sin importar si sus números eran grandes o pequeños."
                    ),
                ),
        );
    }
    let pocas: Vec<Vec<f64>> = puntos
        .iter()
        .take(FILAS_A_LA_VISTA)
        .map(|fila| agrupamiento::estandarizar(fila))
        .collect();
    narracion.agregar(
        tabla_de_numeros(
            Pizarra::nueva().variable("cantidad", cantidad),
            "Filas estandarizadas",
            &tabla.nombres,
            &tabla.columnas,
            &pocas,
            2,
            &[],
        )
        .paso(
            "puntos.iter().map(|punto| estandarizar(punto)).collect()",
            "Lo mismo se hace con todas las filas. Todas quedan con promedio 0 y largo 1, y así se pueden comparar solo por su forma.",
        ),
    );
    let correlacion = agrupamiento::producto_punto(&estandar_a, &estandar_b).clamp(-1.0, 1.0);
    let mut productos: Vec<Vec<String>> = (0..dimension.min(VALORES_A_LA_VISTA))
        .map(|columna| {
            vec![
                tabla.columnas[columna].clone(),
                numero(estandar_a[columna], 3),
                numero(estandar_b[columna], 3),
                numero(estandar_a[columna] * estandar_b[columna], 3),
            ]
        })
        .collect();
    productos.push(vec![
        "suma".to_string(),
        String::new(),
        String::new(),
        numero(correlacion, 3),
    ]);
    let ultima = productos.len() - 1;
    narracion.agregar(
        tabla_de_textos(
            Pizarra::nueva()
                .variable("fila", 0)
                .variable("columna", 1)
                .variable("correlacion", numero(correlacion, 3)),
            "Producto de las dos filas estandarizadas",
            &["columna", nombre_a.as_str(), nombre_b.as_str(), "producto"],
            &productos,
            Some(ultima),
        )
        .paso(
            "let correlacion = producto_punto(&estandarizados[fila], &estandarizados[columna]);",
            format!(
                "Multiplicamos las dos filas columna por columna y sumamos: la correlación r entre {nombre_a} y {nombre_b} es {}. Vale 1 si tienen la misma forma, 0 si no se parecen y −1 si una sube justo cuando la otra baja.",
                numero(correlacion, 3)
            ),
        ),
    );
    let segunda_distancia = 1.0 - correlacion;
    narracion.agregar(
        Pizarra::nueva()
            .variable("correlacion", numero(correlacion, 3))
            .variable("1.0 - correlacion", numero(segunda_distancia, 3))
            .paso(
                ".map(|correlacion| 1.0 - correlacion)",
                format!(
                    "Para convertir la correlación en una distancia se le resta a 1: 0 significa «misma forma» y un valor cerca de 2 significa «formas opuestas». Aquí la distancia de Pearson es {}.",
                    numero(segunda_distancia, 3)
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva().variable("cantidad", cantidad).paso(
            "for columna in fila..cantidad {",
            format!(
                "Los otros pares de filas se calculan igual, también repartidos entre los núcleos de tu computadora: en total {} correlaciones.",
                miles(pares)
            ),
        ),
    );
    narracion.agregar(
        tabla_de_numeros(
            Pizarra::nueva().variable("cantidad", cantidad),
            "Distancia de Pearson (1 − r) entre cada par",
            &tabla.nombres,
            &tabla.nombres,
            &pearson,
            2,
            &[(0, 1), (1, 0)],
        )
        .paso(
            "copiar_al_otro_lado(&mut distancias);",
            "También aquí se completa el triángulo de abajo y la diagonal queda en 0. Las dos matrices están listas.",
        ),
    );
    let mut vecinos = Vec::new();
    let mut resaltar = Vec::new();
    let mut distintos = 0;
    for fila in 0..cantidad {
        let cercano_euclideo = mas_cercano(&euclidea[fila], fila);
        let cercano_pearson = mas_cercano(&pearson[fila], fila);
        if cercano_euclideo != cercano_pearson {
            distintos += 1;
        }
        if fila < FILAS_A_LA_VISTA {
            vecinos.push(vec![
                tabla.nombres[fila].clone(),
                tabla.nombres[cercano_euclideo].clone(),
                numero(euclidea[fila][cercano_euclideo], 2),
                tabla.nombres[cercano_pearson].clone(),
                numero(pearson[fila][cercano_pearson], 2),
            ]);
            if cercano_euclideo != cercano_pearson {
                resaltar.push(fila);
            }
        }
    }
    let ultimo = tabla_de_textos(
        Pizarra::nueva().variable("cantidad", cantidad),
        "Vecino más cercano de cada fila",
        &[
            "fila",
            "más cercana (euclídea)",
            "distancia",
            "más cercana (Pearson)",
            "distancia",
        ],
        &vecinos,
        resaltar.first().copied(),
    )
    .paso(
        "copiar_al_otro_lado(&mut distancias);",
        format!(
            "Resultado: dos matrices de {} × {}. En {} de {} filas el vecino más cercano cambia según la distancia (se marca la primera): la euclídea mira la magnitud y la de Pearson, la forma. Abajo, «Qué significa» te dice cómo leerlas.",
            miles(cantidad),
            miles(cantidad),
            miles(distintos),
            miles(cantidad)
        ),
    );
    (
        narracion.terminar(ultimo),
        [primera_distancia, segunda_distancia],
    )
}

fn pasos_mas_lejano_primero(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let tabla = leer_tabla(entrada, "datos", "p")?;
    let k = leer_k(entrada, "k", tabla.cantidad())?;
    Ok(narrar_mas_lejano(&tabla, k).0)
}

fn narrar_mas_lejano(tabla: &Tabla, k: usize) -> (Vec<Paso>, Vec<usize>) {
    let puntos = &tabla.filas;
    let cantidad = puntos.len();
    let textos = |distancias: &[f64]| -> Vec<String> {
        distancias
            .iter()
            .enumerate()
            .map(|(indice, &lejania)| {
                format!("{}: {}", tabla.nombres[indice], texto_de_distancia(lejania))
            })
            .collect()
    };
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let mut centros = vec![0];
    narracion.agregar(
        tabla_de_numeros(
            Pizarra::nueva()
                .variable("k", k)
                .variable("centros", nombres_de(tabla, &centros)),
            "Tus puntos",
            &tabla.nombres,
            &tabla.columnas,
            puntos,
            2,
            &[(0, 0)],
        )
        .paso(
            "let mut centros = vec![0];",
            format!(
                "Tienes {} puntos y quieres {k} {}. La regla: el primer centro es siempre el primer punto ({}); cada centro nuevo es el punto más lejano de todos los centros ya elegidos.",
                miles(cantidad),
                segun_cantidad(k, "centro", "centros"),
                tabla.nombres[0]
            ),
        ),
    );
    let mut distancia_al_centro = vec![f64::INFINITY; cantidad];
    narracion.agregar(
        lista_de(
            Pizarra::nueva().variable("k", k).variable("centros", nombres_de(tabla, &centros)),
            "Distancia de cada punto a su centro más cercano",
            &textos(&distancia_al_centro),
            &[],
        )
        .paso(
            "let mut distancia_al_centro = vec![f64::INFINITY; puntos.len()];",
            "Para cada punto guardamos qué tan lejos queda de su centro más cercano. Al principio no hay ninguna distancia conocida, así que todas valen infinito (∞).",
        ),
    );
    agrupamiento::acercar_al_nuevo_centro(&mut distancia_al_centro, puntos, &puntos[0]);
    let primer_lejano = agrupamiento::posicion_del_maximo(&distancia_al_centro);
    narracion.agregar(
        lista_de(
            Pizarra::nueva().variable("k", k).variable("centros", nombres_de(tabla, &centros)),
            "Distancia de cada punto a su centro más cercano",
            &textos(&distancia_al_centro),
            &[primer_lejano],
        )
        .paso(
            "acercar_al_nuevo_centro(&mut distancia_al_centro, puntos, &puntos[0]);",
            format!(
                "Medimos la distancia de cada punto al centro {}. El que queda más lejos es {} (a {}).",
                tabla.nombres[0],
                tabla.nombres[primer_lejano],
                numero(distancia_al_centro[primer_lejano], 2)
            ),
        ),
    );
    distancia_al_centro[0] = YA_ES_CENTRO;
    narracion.agregar(
        lista_de(
            Pizarra::nueva().variable("k", k).variable("centros", nombres_de(tabla, &centros)),
            "Distancia de cada punto a su centro más cercano",
            &textos(&distancia_al_centro),
            &[0],
        )
        .paso(
            "distancia_al_centro[0] = YA_ES_CENTRO;",
            format!(
                "Marcamos a {} como «centro» (con −1) para que nunca vuelva a ser elegido como el más lejano.",
                tabla.nombres[0]
            ),
        ),
    );
    let salto = salto_de(k.saturating_sub(1).max(1));
    while centros.len() < k {
        let vuelta = centros.len() - 1;
        let detallada = toca_narrar(vuelta, salto);
        let completa = vuelta < PASOS_DETALLADOS;
        if detallada && completa {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("k", k)
                    .variable("centros", nombres_de(tabla, &centros))
                    .paso(
                        "while centros.len() < k {",
                        format!(
                            "Tenemos {} de {k} centros, así que falta elegir otro.",
                            centros.len()
                        ),
                    ),
            );
        }
        let mas_lejano = agrupamiento::posicion_del_maximo(&distancia_al_centro);
        if detallada && completa {
            narracion.agregar(
                lista_de(
                    Pizarra::nueva()
                        .variable("k", k)
                        .variable("centros", nombres_de(tabla, &centros))
                        .variable("mas_lejano", mas_lejano),
                    "Distancia de cada punto a su centro más cercano",
                    &textos(&distancia_al_centro),
                    &[mas_lejano],
                )
                .paso(
                    "let mas_lejano = posicion_del_maximo(&distancia_al_centro);",
                    format!(
                        "Buscamos la distancia más grande de la lista: {} está a {} de su centro más cercano, más lejos que cualquier otro punto.",
                        tabla.nombres[mas_lejano],
                        numero(distancia_al_centro[mas_lejano], 2)
                    ),
                ),
            );
        }
        centros.push(mas_lejano);
        if detallada {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("k", k)
                    .variable("centros", nombres_de(tabla, &centros))
                    .variable("mas_lejano", mas_lejano)
                    .paso(
                        "centros.push(mas_lejano);",
                        format!(
                            "{} es el centro número {}. Los centros hasta ahora: {}.",
                            tabla.nombres[mas_lejano],
                            centros.len(),
                            nombres_de(tabla, &centros)
                        ),
                    ),
            );
        }
        let antes = distancia_al_centro.clone();
        agrupamiento::acercar_al_nuevo_centro(
            &mut distancia_al_centro,
            puntos,
            &puntos[mas_lejano],
        );
        let cambiados: Vec<usize> = (0..cantidad)
            .filter(|&indice| distancia_al_centro[indice] < antes[indice])
            .collect();
        if detallada && completa {
            let visibles: Vec<usize> = cambiados
                .iter()
                .copied()
                .filter(|&indice| indice < ELEMENTOS_A_LA_VISTA)
                .collect();
            narracion.agregar(
                lista_de(
                    Pizarra::nueva()
                        .variable("k", k)
                        .variable("centros", nombres_de(tabla, &centros))
                        .variable("mas_lejano", mas_lejano),
                    "Distancia de cada punto a su centro más cercano",
                    &textos(&distancia_al_centro),
                    &visibles,
                )
                .paso(
                    "acercar_al_nuevo_centro(&mut distancia_al_centro, puntos, &puntos[mas_lejano]);",
                    format!(
                        "Cada punto compara su distancia vieja con su distancia al centro nuevo y se queda con la menor. {} {} se acercaron al centro nuevo.",
                        miles(cambiados.len()),
                        segun_cantidad(cambiados.len(), "punto", "puntos")
                    ),
                ),
            );
        }
        distancia_al_centro[mas_lejano] = YA_ES_CENTRO;
        if detallada && completa {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("k", k)
                    .variable("centros", nombres_de(tabla, &centros))
                    .variable("mas_lejano", mas_lejano)
                    .paso(
                        "distancia_al_centro[mas_lejano] = YA_ES_CENTRO;",
                        format!(
                            "{} también se marca como «centro». Se repite hasta tener {k}.",
                            tabla.nombres[mas_lejano]
                        ),
                    ),
            );
        }
    }
    let radio = distancia_al_centro
        .iter()
        .fold(0.0f64, |mayor, &lejania| mayor.max(lejania));
    let ultimo = lista_de(
        Pizarra::nueva()
            .variable("k", k)
            .variable("centros", nombres_de(tabla, &centros)),
        "Centros en el orden en que se eligieron",
        &centros
            .iter()
            .enumerate()
            .map(|(orden, &indice)| format!("{}. {}", orden + 1, tabla.nombres[indice]))
            .collect::<Vec<String>>(),
        &[],
    )
    .paso(
        "Ok(centros)",
        format!(
            "Listo: los {k} centros son {}. Ningún punto queda a más de {} de su centro. En «Qué significa» verás cómo leer ese radio.",
            nombres_de(tabla, &centros),
            numero(radio, 3)
        ),
    );
    (narracion.terminar(ultimo), centros)
}

fn pasos_distorsion(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let tabla = leer_tabla(entrada, "datos", "p")?;
    let centros = leer_tabla(entrada, "centros", "centro ")?;
    if centros.dimension() != tabla.dimension() {
        return Err(format!(
            "Los puntos tienen {} coordenadas y los centros {}: deben tener las mismas.",
            tabla.dimension(),
            centros.dimension()
        ));
    }
    if centros.cantidad() > MAXIMO_K {
        return Err(format!(
            "Hay {} centros; el máximo es {MAXIMO_K}.",
            centros.cantidad()
        ));
    }
    Ok(narrar_distorsion(&tabla, &centros).0)
}

fn cuenta_al_cuadrado(punto: &[f64], centro: &[f64]) -> String {
    let mut partes: Vec<String> = punto
        .iter()
        .zip(centro)
        .take(3)
        .map(|(x, y)| format!("({} − {})²", numero(*x, 1), numero(*y, 1)))
        .collect();
    if punto.len() > 3 {
        partes.push("…".to_string());
    }
    partes.join(" + ")
}

fn narrar_distorsion(tabla: &Tabla, centros: &Tabla) -> (Vec<Paso>, f64) {
    let puntos = &tabla.filas;
    let cantidad = puntos.len();
    let cuantos = centros.cantidad();
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let nombres_de_centros: Vec<String> = (0..cuantos)
        .map(|centro| {
            format!(
                "{}: {}",
                centros.nombres[centro],
                vector_corto(&centros.filas[centro])
            )
        })
        .collect();
    narracion.agregar(
        lista_de(
            tabla_de_numeros(
                Pizarra::nueva().variable("suma", 0),
                "Tus puntos",
                &tabla.nombres,
                &tabla.columnas,
                puntos,
                2,
                &[],
            ),
            "Tus centros",
            &nombres_de_centros,
            &[],
        )
        .paso(
            "let mut suma = 0.0;",
            format!(
                "Tienes {} puntos y {} {}. Para cada punto buscamos su centro más cercano, elevamos esa distancia al cuadrado y lo vamos sumando; al final se promedia. Empezamos con la suma en 0.",
                miles(cantidad),
                miles(cuantos),
                segun_cantidad(cuantos, "centro", "centros")
            ),
        ),
    );
    let salto = salto_de(cantidad);
    let mut suma = 0.0;
    for (indice, punto) in puntos.iter().enumerate() {
        let mut menor = f64::INFINITY;
        let mut mejor = 0;
        let mut cuentas = Vec::with_capacity(cuantos);
        for (numero_de_centro, centro) in centros.filas.iter().enumerate() {
            let cuadrado = agrupamiento::distancia_al_cuadrado(punto, centro);
            if cuadrado < menor {
                mejor = numero_de_centro;
            }
            menor = menor.min(cuadrado);
            cuentas.push(cuadrado);
        }
        let antes_de_sumar = suma;
        suma += menor;
        let completa = indice < PASOS_DETALLADOS;
        if completa {
            let filas: Vec<Vec<String>> = (0..cuantos)
                .map(|centro| {
                    vec![
                        centros.nombres[centro].clone(),
                        cuenta_al_cuadrado(punto, &centros.filas[centro]),
                        numero(cuentas[centro], 3),
                    ]
                })
                .collect();
            narracion.agregar(
                tabla_de_textos(
                    Pizarra::nueva()
                        .variable("suma", numero(antes_de_sumar, 3))
                        .variable("menor", numero(menor, 3)),
                    &format!("Distancia al cuadrado de {} a cada centro", tabla.nombres[indice]),
                    &["centro", "cuenta", "distancia²"],
                    &filas,
                    Some(mejor),
                )
                .paso(
                    "menor = menor.min(distancia_al_cuadrado(punto, centro));",
                    format!(
                        "{} está en {}. Su distancia al cuadrado a cada centro está en la tabla y la menor es {}, la de {}: ese es su centro más cercano.",
                        tabla.nombres[indice],
                        vector_corto(punto),
                        numero(menor, 3),
                        centros.nombres[mejor]
                    ),
                ),
            );
        }
        if toca_narrar(indice, salto) {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("suma", numero(suma, 3))
                    .variable("menor", numero(menor, 3))
                    .paso(
                        "suma += menor_distancia_al_cuadrado(punto, centros);",
                        format!(
                            "Sumamos {} (de {} con {}) a la suma, que ahora vale {} tras {} {}.",
                            numero(menor, 3),
                            tabla.nombres[indice],
                            centros.nombres[mejor],
                            numero(suma, 3),
                            indice + 1,
                            segun_cantidad(indice + 1, "punto", "puntos")
                        ),
                    ),
            );
        }
    }
    let valor = suma / cantidad as f64;
    let ultimo = Pizarra::nueva()
        .variable("suma", numero(suma, 3))
        .variable("puntos.len()", cantidad)
        .variable("distorsión", cientifico(valor))
        .paso(
            "suma / puntos.len() as f64",
            format!(
                "Dividimos la suma ({}) entre los {} puntos: la distorsión es {}. Su raíz, {}, es la distancia típica de un punto a su centro en las unidades de tus datos. Mientras más baja, mejor ubicados están los centros: «Qué significa» lo interpreta.",
                numero(suma, 3),
                miles(cantidad),
                cientifico(valor),
                cientifico(valor.sqrt())
            ),
        );
    (narracion.terminar(ultimo), valor)
}

fn pasos_codo(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let tabla = leer_tabla(entrada, "datos", "p")?;
    exigir_filas(&tabla, 2)?;
    let k_maximo = entrada.natural("k_maximo")?;
    if !(2..=MAXIMO_K_DEL_CODO).contains(&k_maximo) {
        return Err(format!(
            "El k máximo debe estar entre 2 y {MAXIMO_K_DEL_CODO}."
        ));
    }
    let semilla = leer_semilla(entrada)?;
    let hasta = k_maximo.min(tabla.cantidad());
    limitar_cuentas((tabla.cantidad() * tabla.dimension() * hasta * hasta) as f64 * 150.0)?;
    Ok(narrar_codo(&tabla, k_maximo, semilla)?.0)
}

fn narrar_codo(
    tabla: &Tabla,
    k_maximo: usize,
    semilla: u64,
) -> Result<(Vec<Paso>, Vec<f64>), String> {
    let puntos = &tabla.filas;
    let hasta = k_maximo.min(puntos.len());
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        tabla_de_numeros(
            Pizarra::nueva()
                .variable("k_maximo", k_maximo)
                .variable("semilla", semilla),
            "Tus datos",
            &tabla.nombres,
            &tabla.columnas,
            puntos,
            2,
            &[],
        )
        .paso(
            "let mut curva = Vec::new();",
            format!(
                "Vamos a correr k-means con k = 1, 2, 3… hasta {hasta} y anotar la distorsión de cada uno (qué tan lejos queda cada punto de su centro). Esa lista de números es la curva. Como k-means depende del azar del inicio, para cada k probamos {INTENTOS_DEL_CODO} inicios distintos y nos quedamos con el mejor."
            ),
        ),
    );
    let mut curva: Vec<f64> = Vec::new();
    let mut historia: Vec<Vec<String>> = Vec::new();
    for k in 1..=hasta {
        let mut mejor = f64::INFINITY;
        let mut intentos: Vec<String> = Vec::new();
        let mut mejor_intento = 0;
        for intento in 0..INTENTOS_DEL_CODO {
            let semilla_usada = semilla.wrapping_add(intento);
            let inicio = Inicio::KMediasMasMas(semilla_usada);
            let resultado = agrupamiento::k_medias(puntos, k, inicio, ITERACIONES_DEL_CODO)?;
            if resultado.distorsion < mejor {
                mejor_intento = intentos.len();
            }
            mejor = mejor.min(resultado.distorsion);
            intentos.push(format!(
                "intento {} (semilla {semilla_usada}): distorsión {} en {} {}",
                intento + 1,
                cientifico(resultado.distorsion),
                miles(resultado.iteraciones),
                segun_cantidad(resultado.iteraciones, "iteración", "iteraciones")
            ));
            if k < PASOS_DETALLADOS {
                let explicacion = if k == 1 {
                    "Con k = 1 hay un solo grupo: su centro es el promedio de todos los puntos y la distorsión mide cuánto se dispersan los datos en total.".to_string()
                } else {
                    format!(
                        "Probamos k = {k} con k-means++ (semilla {semilla_usada}); {} tras {} {}. La distorsión es {}.",
                        if resultado.convergio {
                            "los grupos se estabilizaron"
                        } else {
                            "llegó al máximo de iteraciones"
                        },
                        miles(resultado.iteraciones),
                        segun_cantidad(resultado.iteraciones, "iteración", "iteraciones"),
                        cientifico(resultado.distorsion)
                    )
                };
                narracion.agregar(
                    lista_de(
                        Pizarra::nueva()
                            .variable("k", k)
                            .variable("mejor", cientifico(mejor))
                            .variable("intento", intento)
                            .variable("resultado.distorsion", cientifico(resultado.distorsion)),
                        "Inicios probados con este k",
                        &intentos,
                        &[intentos.len() - 1],
                    )
                    .paso(
                        "let resultado = k_medias(puntos, k, inicio, ITERACIONES_DEL_CODO)?;",
                        explicacion,
                    ),
                );
            }
        }
        curva.push(mejor);
        let caida = if k == 1 || curva[k - 2] <= 0.0 {
            "—".to_string()
        } else {
            porcentaje(agrupamiento::caida_relativa(curva[k - 2], mejor), 0)
        };
        historia.push(vec![format!("{k}"), cientifico(mejor), caida]);
        if k < PASOS_DETALLADOS || k == hasta || k.is_multiple_of(2) {
            narracion.agregar(
                tabla_de_textos(
                    lista_de(
                        Pizarra::nueva().variable("k", k).variable("mejor", cientifico(mejor)),
                        "Inicios probados con este k",
                        &intentos,
                        &[mejor_intento],
                    ),
                    "La curva hasta ahora",
                    &["k", "distorsión", "caída respecto de k − 1"],
                    &historia,
                    Some(k - 1),
                )
                .paso(
                    "curva.push(mejor);",
                    format!(
                        "Con k = {k}, el mejor de los {INTENTOS_DEL_CODO} inicios deja una distorsión de {}; esa es la que anotamos en la curva.",
                        cientifico(mejor)
                    ),
                ),
            );
        }
    }
    let elegido = agrupamiento::k_del_codo(&curva);
    let filas: Vec<Vec<String>> = (0..curva.len())
        .map(|indice| {
            let frenada = if indice >= 1 && indice + 1 < curva.len() {
                let al_llegar = agrupamiento::caida_relativa(curva[indice - 1], curva[indice]);
                let al_seguir = agrupamiento::caida_relativa(curva[indice], curva[indice + 1]);
                porcentaje(al_llegar - al_seguir, 0)
            } else {
                "—".to_string()
            };
            vec![
                format!("{}", indice + 1),
                cientifico(curva[indice]),
                historia[indice][2].clone(),
                frenada,
            ]
        })
        .collect();
    let nota = match elegido {
        Some(k) => format!(
            "La curva está lista. Para elegir el codo se mira, en cada k, cuánto se frena la caída: la caída al llegar a ese k menos la caída al seguir al siguiente. La mayor frenada está en k = {k}: ahí agregar un grupo más deja de ayudar mucho. Abajo, «Qué significa» lo explica."
        ),
        None => "La curva está lista. Para elegir el codo se mira, en cada k, cuánto se frena la caída, pero aquí ninguna frena de golpe (o hay muy pocos valores de k): no hay un codo claro. Abajo, «Qué significa» lo explica.".to_string(),
    };
    let ultimo = tabla_de_textos(
        Pizarra::nueva().variable("k_maximo", k_maximo),
        "La curva, número por número",
        &["k", "distorsión", "caída respecto de k − 1", "frenada"],
        &filas,
        elegido.map(|k| k - 1),
    )
    .paso("Ok(curva)", nota);
    Ok((narracion.terminar(ultimo), curva))
}

fn pasos_k_medias_suaves(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let tabla = leer_tabla(entrada, "datos", "p")?;
    let k = leer_k(entrada, "k", tabla.cantidad())?;
    let rigidez = entrada.decimal("rigidez")?;
    if rigidez <= 0.0 {
        return Err("La rigidez β debe ser mayor que 0.".to_string());
    }
    let iteraciones = leer_iteraciones(entrada, "iteraciones")?;
    limitar_cuentas((tabla.cantidad() * k * tabla.dimension() * iteraciones) as f64 * 2.0)?;
    Ok(narrar_suaves(&tabla, k, rigidez, iteraciones).0)
}

fn centros_con_nombres(k: usize) -> Vec<String> {
    (1..=k).map(|grupo| format!("grupo {grupo}")).collect()
}

fn narrar_suaves(
    tabla: &Tabla,
    k: usize,
    rigidez: f64,
    iteraciones: usize,
) -> (Vec<Paso>, KMediasSuaves) {
    let puntos = &tabla.filas;
    let cantidad = puntos.len();
    let grupos = centros_con_nombres(k);
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let mut centros: Vec<Vec<f64>> = puntos[..k].to_vec();
    narracion.agregar(
        tabla_de_numeros(
            Pizarra::nueva()
                .variable("k", k)
                .variable("rigidez", numero(rigidez, 2))
                .variable("iteraciones", iteraciones),
            "Los centros iniciales (los primeros k puntos)",
            &grupos,
            &tabla.columnas,
            &centros,
            2,
            &[],
        )
        .paso(
            "let mut centros: Vec<Vec<f64>> = puntos[..k].to_vec();",
            format!(
                "Empezamos con {k} centros: los primeros {k} puntos ({}). A diferencia de k-means normal, aquí cada punto no se va a un solo centro: pertenece a todos, un poco a cada uno. La rigidez β = {} decide qué tan fuerte es esa preferencia (más grande, más se parece a k-means normal).",
                nombres_de(tabla, &(0..k).collect::<Vec<usize>>()),
                numero(rigidez, 2)
            ),
        ),
    );
    let salto = salto_de(iteraciones);
    for vuelta in 0..iteraciones {
        if vuelta == 0 {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("k", k)
                    .variable("rigidez", numero(rigidez, 2))
                    .variable("iteraciones", iteraciones)
                    .paso(
                        "for _ in 0..iteraciones {",
                        format!(
                            "Cada vuelta tiene dos mitades: primero cada punto calcula cuánto le pertenece cada grupo (responsabilidades) y luego cada centro se mueve según esas responsabilidades. Haremos hasta {} {}.",
                            miles(iteraciones),
                            segun_cantidad(iteraciones, "vuelta", "vueltas")
                        ),
                    ),
            );
            let punto = &puntos[0];
            let distancias: Vec<f64> = centros
                .iter()
                .map(|centro| agrupamiento::distancia_euclidea(punto, centro))
                .collect();
            narracion.agregar(
                tabla_de_numeros(
                    Pizarra::nueva().variable("distancias", vector_corto(&distancias)),
                    &format!("Distancia de {} a cada centro", tabla.nombres[0]),
                    &grupos,
                    &["distancia".to_string()],
                    &distancias.iter().map(|&valor| vec![valor]).collect::<Vec<_>>(),
                    2,
                    &[],
                )
                .paso(
                    ".map(|centro| distancia_euclidea(punto, centro))",
                    format!(
                        "Tomamos el primer punto, {}, y medimos su distancia a cada uno de los {k} centros.",
                        tabla.nombres[0]
                    ),
                ),
            );
            let mas_cerca = distancias.iter().copied().fold(f64::INFINITY, f64::min);
            narracion.agregar(
                Pizarra::nueva()
                    .variable("distancias", vector_corto(&distancias))
                    .variable("mas_cerca", numero(mas_cerca, 3))
                    .paso(
                        "let mas_cerca = distancias.iter().copied().fold(f64::INFINITY, f64::min);",
                        format!(
                            "La distancia más corta es {}: la usamos como referencia para que los números de abajo no se hagan enormes ni diminutos.",
                            numero(mas_cerca, 3)
                        ),
                    ),
            );
            let atracciones: Vec<f64> = distancias
                .iter()
                .map(|lejania| (-rigidez * (lejania - mas_cerca)).exp())
                .collect();
            let filas_de_atraccion: Vec<Vec<String>> = (0..k)
                .map(|grupo| {
                    vec![
                        grupos[grupo].clone(),
                        numero(distancias[grupo], 3),
                        numero(atracciones[grupo], 3),
                    ]
                })
                .collect();
            narracion.agregar(
                tabla_de_textos(
                    Pizarra::nueva()
                        .variable("rigidez", numero(rigidez, 2))
                        .variable("mas_cerca", numero(mas_cerca, 3)),
                    "Atracción de cada centro sobre el primer punto",
                    &["centro", "distancia", "atracción"],
                    &filas_de_atraccion,
                    None,
                )
                .paso(
                    ".map(|lejania| (-rigidez * (lejania - mas_cerca)).exp())",
                    "La atracción de un centro es e elevado a −β por cuánto más lejos está que el más cercano. El más cercano siempre vale 1; los demás valen menos mientras más lejos estén (y caen más rápido con una rigidez grande).",
                ),
            );
            let total: f64 = atracciones.iter().sum();
            narracion.agregar(Pizarra::nueva().variable("total", numero(total, 3)).paso(
                "let total: f64 = atracciones.iter().sum();",
                format!("Sumamos todas las atracciones: {}.", numero(total, 3)),
            ));
            let responsabilidades_del_primero: Vec<f64> = atracciones
                .iter()
                .map(|atraccion| atraccion / total)
                .collect();
            let filas_de_responsabilidad: Vec<Vec<String>> = (0..k)
                .map(|grupo| {
                    vec![
                        grupos[grupo].clone(),
                        numero(atracciones[grupo], 3),
                        numero(responsabilidades_del_primero[grupo], 3),
                    ]
                })
                .collect();
            let (mayor, donde) = extremo(&responsabilidades_del_primero, true);
            narracion.agregar(
                tabla_de_textos(
                    Pizarra::nueva().variable("total", numero(total, 3)),
                    "Responsabilidad de cada grupo sobre el primer punto",
                    &["grupo", "atracción", "responsabilidad"],
                    &filas_de_responsabilidad,
                    Some(donde),
                )
                .paso(
                    "*atraccion /= total;",
                    format!(
                        "Dividimos cada atracción entre el total para que las responsabilidades sumen 1. {} pertenece un {} al {} y el resto se reparte entre los demás grupos.",
                        tabla.nombres[0],
                        porcentaje(mayor, 0),
                        grupos[donde]
                    ),
                ),
            );
        }
        let responsabilidades = agrupamiento::de_centros_a_grupos_suaves(puntos, &centros, rigidez);
        if vuelta == 0 {
            narracion.agregar(
                tabla_de_numeros(
                    Pizarra::nueva().variable("rigidez", numero(rigidez, 2)),
                    "Responsabilidades de cada punto",
                    &tabla.nombres,
                    &grupos,
                    &responsabilidades,
                    2,
                    &[(0, 0)],
                )
                .paso(
                    "let responsabilidades = de_centros_a_grupos_suaves(puntos, &centros, rigidez);",
                    format!(
                        "Lo mismo se hace con los {} puntos (en paralelo). Cada fila de la tabla suma 1: es cómo se reparte ese punto entre los grupos.",
                        miles(cantidad)
                    ),
                ),
            );
        }
        let nuevos = agrupamiento::de_grupos_suaves_a_centros(puntos, &responsabilidades, &centros);
        if vuelta == 0 {
            let mut suma = vec![0.0; tabla.dimension()];
            let mut peso_total = 0.0;
            let mut pesos = Vec::new();
            for (punto, pesos_del_punto) in puntos.iter().zip(&responsabilidades) {
                let peso = pesos_del_punto[0];
                peso_total += peso;
                pesos.push(vec![peso]);
                for (acumulado, valor) in suma.iter_mut().zip(punto) {
                    *acumulado += peso * valor;
                }
            }
            narracion.agregar(
                tabla_de_numeros(
                    Pizarra::nueva()
                        .variable("grupo", 0)
                        .variable("peso_total", numero(peso_total, 3)),
                    "Peso de cada punto en el grupo 1",
                    &tabla.nombres,
                    &["peso".to_string()],
                    &pesos,
                    2,
                    &[],
                )
                .paso(
                    "peso_total += peso;",
                    format!(
                        "Ahora los centros. Para el grupo 1 sumamos la responsabilidad que tiene sobre cada punto: {}. Es como «cuántos puntos» le pertenecen, pero no tiene que ser un entero.",
                        numero(peso_total, 3)
                    ),
                ),
            );
            narracion.agregar(
                Pizarra::nueva()
                    .variable("peso_total", numero(peso_total, 3))
                    .variable("suma", vector_corto(&suma))
                    .paso(
                        "*acumulado += peso * valor;",
                        "Cada punto aporta sus coordenadas multiplicadas por su peso: un punto que casi no pertenece al grupo casi no cuenta, y uno que le pertenece del todo cuenta completo.",
                    ),
            );
            let nuevo_a_mano: Vec<f64> = if peso_total > 0.0 {
                suma.iter().map(|valor| valor / peso_total).collect()
            } else {
                centros[0].clone()
            };
            narracion.agregar(
                Pizarra::nueva()
                    .variable("peso_total", numero(peso_total, 3))
                    .variable("centro nuevo", vector_corto(&nuevo_a_mano))
                    .paso(
                        "suma.into_iter().map(|valor| valor / peso_total).collect()",
                        format!(
                            "Dividimos entre el peso total: es un promedio ponderado de los puntos. El centro del grupo 1 pasa de {} a {}.",
                            vector_corto(&centros[0]),
                            vector_corto(&nuevo_a_mano)
                        ),
                    ),
            );
        }
        let movimiento = centros
            .iter()
            .zip(&nuevos)
            .map(|(antes, despues)| agrupamiento::distancia_euclidea(antes, despues))
            .fold(0.0f64, f64::max);
        centros = nuevos;
        if toca_narrar(vuelta, salto) {
            let nota = if vuelta == 0 {
                format!(
                    "Los {k} centros se mueven así. El que más se movió recorrió {}. Con esto termina la primera vuelta.",
                    numero(movimiento, 3)
                )
            } else if movimiento < TOLERANCIA {
                format!(
                    "Vuelta {}: los centros ya casi no se mueven (como máximo {}). Las vueltas que faltan darán lo mismo.",
                    miles(vuelta + 1),
                    cientifico(movimiento)
                )
            } else {
                format!(
                    "Vuelta {}: los centros se mueven otra vez, como máximo {}. Con centros nuevos, las responsabilidades cambian y se repite.",
                    miles(vuelta + 1),
                    numero(movimiento, 3)
                )
            };
            narracion.agregar(
                tabla_de_numeros(
                    Pizarra::nueva()
                        .variable("k", k)
                        .variable("rigidez", numero(rigidez, 2))
                        .variable("iteraciones", iteraciones),
                    "Los centros ahora",
                    &grupos,
                    &tabla.columnas,
                    &centros,
                    2,
                    &[],
                )
                .paso(
                    "centros = de_grupos_suaves_a_centros(puntos, &responsabilidades, &centros);",
                    nota,
                ),
            );
        }
    }
    let responsabilidades = agrupamiento::de_centros_a_grupos_suaves(puntos, &centros, rigidez);
    let firmes = responsabilidades
        .iter()
        .filter(|fila| extremo(fila, true).0 >= MAYORIA)
        .count();
    let marcados: Vec<(usize, usize)> = responsabilidades
        .iter()
        .enumerate()
        .take(FILAS_A_LA_VISTA)
        .map(|(fila, valores)| (fila, extremo(valores, true).1))
        .collect();
    let ultimo = tabla_de_numeros(
        Pizarra::nueva()
            .variable("k", k)
            .variable("rigidez", numero(rigidez, 2)),
        "Responsabilidades finales",
        &tabla.nombres,
        &grupos,
        &responsabilidades,
        2,
        &marcados,
    )
    .paso(
        "Ok(KMediasSuaves {",
        format!(
            "Con los centros finales se calculan las responsabilidades una última vez. {} de {} filas pertenecen al menos un {} a un solo grupo; las demás están en la frontera entre grupos. Abajo, «Qué significa» te dice cómo leerlo.",
            miles(firmes),
            miles(cantidad),
            porcentaje(MAYORIA, 0)
        ),
    );
    (
        narracion.terminar(ultimo),
        KMediasSuaves {
            centros,
            responsabilidades,
        },
    )
}

fn pasos_jerarquico(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let tabla = leer_tabla(entrada, "datos", "p")?;
    let eleccion = entrada.opcion("distancia", MEDIDAS)?;
    exigir_filas(&tabla, 2)?;
    limitar_matriz(&tabla)?;
    let es_matriz = eleccion == MEDIDA_MATRIZ;
    let distancias = if es_matriz {
        revisar_matriz_de_distancias(&tabla)?;
        tabla.filas.clone()
    } else if eleccion == MEDIDA_PEARSON {
        agrupamiento::matriz_de_distancias(&tabla.filas, Medida::Pearson)
    } else {
        agrupamiento::matriz_de_distancias(&tabla.filas, Medida::Euclidea)
    };
    Ok(narrar_jerarquico(&tabla, &distancias, &eleccion).0)
}

fn narrar_jerarquico(
    tabla: &Tabla,
    distancias: &[Vec<f64>],
    medida: &str,
) -> (Vec<Paso>, Vec<Union>) {
    let cantidad = distancias.len();
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        tabla_de_numeros(
            Pizarra::nueva().variable("cantidad", cantidad),
            &format!("Distancias entre filas ({medida})"),
            &tabla.nombres,
            &tabla.nombres,
            distancias,
            2,
            &[],
        )
        .paso(
            "let cantidad = distancias.len();",
            format!(
                "Tienes {} filas y la tabla de distancias entre cada par. Al principio cada fila es su propio grupo. La meta: unir una y otra vez los dos grupos más cercanos hasta que quede uno solo, y anotar cada unión. Eso es el árbol.",
                miles(cantidad)
            ),
        ),
    );
    let mut entre_grupos = distancias.to_vec();
    let mut miembros: Vec<Vec<usize>> = (0..cantidad).map(|punto| vec![punto]).collect();
    let mut nodo: Vec<usize> = (0..cantidad).collect();
    let mut activo = vec![true; cantidad];
    let mut vecino: Vec<usize> = (0..cantidad)
        .map(|grupo| agrupamiento::vecino_mas_cercano(grupo, &entre_grupos, &activo))
        .collect();
    let vecinos: Vec<String> = (0..cantidad)
        .map(|grupo| {
            format!(
                "{} → {} ({})",
                tabla.nombres[grupo],
                tabla.nombres[vecino[grupo]],
                numero(entre_grupos[grupo][vecino[grupo]], 2)
            )
        })
        .collect();
    narracion.agregar(
        lista_de(
            Pizarra::nueva().variable("cantidad", cantidad),
            "Vecino más cercano de cada grupo",
            &vecinos,
            &[],
        )
        .paso(
            "let mut vecino: Vec<usize> = (0..cantidad)",
            "Para no buscar en toda la tabla cada vez, cada grupo anota cuál es su vecino más cercano. El par más cercano de todos saldrá de esta lista.",
        ),
    );
    let mut uniones: Vec<Union> = Vec::with_capacity(cantidad - 1);
    let mut registro: Vec<String> = Vec::new();
    let salto = salto_de(cantidad - 1);
    for paso in 0..cantidad - 1 {
        let (a, b) = agrupamiento::par_mas_cercano(&entre_grupos, &vecino, &activo);
        let distancia = entre_grupos[a][b];
        let peso_a = miembros[a].len() as f64;
        let peso_b = miembros[b].len() as f64;
        let etiqueta_a = etiqueta_de_grupo(tabla, &miembros[a]);
        let etiqueta_b = etiqueta_de_grupo(tabla, &miembros[b]);
        let completa = paso < PASOS_DETALLADOS;
        if completa {
            let vivos: Vec<usize> = (0..cantidad).filter(|&grupo| activo[grupo]).collect();
            let textos: Vec<String> = vivos
                .iter()
                .map(|&grupo| etiqueta_de_grupo(tabla, &miembros[grupo]))
                .collect();
            narracion.agregar(
                lista_de(
                    Pizarra::nueva()
                        .variable("paso", paso)
                        .variable("a", a)
                        .variable("b", b)
                        .variable("distancia", numero(distancia, 3)),
                    "Grupos que quedan",
                    &textos,
                    &posiciones(&vivos, &[a, b]),
                )
                .paso(
                    "let (a, b) = par_mas_cercano(&entre_grupos, &vecino, &activo);",
                    format!(
                        "Unión {}: de todos los grupos que quedan, los dos más cercanos son {etiqueta_a} y {etiqueta_b}, a distancia {}.",
                        paso + 1,
                        numero(distancia, 3)
                    ),
                ),
            );
            let ejemplo = (0..cantidad).find(|&otro| activo[otro] && otro != a && otro != b);
            if let Some(otro) = ejemplo {
                let promedio = (peso_a * entre_grupos[a][otro] + peso_b * entre_grupos[b][otro])
                    / (peso_a + peso_b);
                narracion.agregar(
                    Pizarra::nueva()
                        .variable("peso_a", peso_a)
                        .variable("peso_b", peso_b)
                        .variable("otro", otro)
                        .variable("promedio", numero(promedio, 3))
                        .paso(
                            "let promedio = (peso_a * entre_grupos[a][otro] + peso_b * entre_grupos[b][otro])",
                            format!(
                                "El grupo nuevo necesita su distancia a cada otro grupo. Para {}: ({} × {} + {} × {}) / ({} + {}) = {}. Es el promedio de las dos distancias, contando más al grupo que tiene más filas.",
                                etiqueta_de_grupo(tabla, &miembros[otro]),
                                peso_a,
                                numero(entre_grupos[a][otro], 3),
                                peso_b,
                                numero(entre_grupos[b][otro], 3),
                                peso_a,
                                peso_b,
                                numero(promedio, 3)
                            ),
                        ),
                );
            }
        }
        for otro in 0..cantidad {
            if activo[otro] && otro != a && otro != b {
                let promedio = (peso_a * entre_grupos[a][otro] + peso_b * entre_grupos[b][otro])
                    / (peso_a + peso_b);
                entre_grupos[a][otro] = promedio;
                entre_grupos[otro][a] = promedio;
            }
        }
        activo[b] = false;
        let absorbidos = std::mem::take(&mut miembros[b]);
        miembros[a].extend(absorbidos);
        if completa {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("paso", paso)
                    .variable("a", a)
                    .variable("b", b)
                    .paso(
                        "activo[b] = false;",
                        format!(
                            "{etiqueta_b} se apaga: ya no es un grupo aparte, sus filas pasan a {etiqueta_a}. Ahora ese grupo tiene {} filas.",
                            miembros[a].len()
                        ),
                    ),
            );
        }
        uniones.push(Union {
            primero: nodo[a],
            segundo: nodo[b],
            distancia,
            miembros: miembros[a].clone(),
        });
        registro.push(format!(
            "{}. {etiqueta_a} + {etiqueta_b} a {}",
            paso + 1,
            numero(distancia, 3)
        ));
        nodo[a] = cantidad + paso;
        if toca_narrar(paso, salto) {
            narracion.agregar(
                lista_de(
                    Pizarra::nueva()
                        .variable("paso", paso)
                        .variable("a", a)
                        .variable("b", b)
                        .variable("distancia", numero(distancia, 3)),
                    "Uniones hechas",
                    &registro,
                    &[registro.len() - 1],
                )
                .paso(
                    "uniones.push(Union {",
                    format!(
                        "Anotamos la unión {}: {etiqueta_a} con {etiqueta_b}, a distancia {}. {}",
                        paso + 1,
                        numero(distancia, 3),
                        if paso == 0 {
                            "Cada unión es una rama del árbol y su distancia es la altura a la que se juntan."
                        } else {
                            "Quedan menos grupos en cada vuelta."
                        }
                    ),
                ),
            );
        }
        for otro in 0..cantidad {
            if activo[otro] && (otro == a || vecino[otro] == a || vecino[otro] == b) {
                vecino[otro] = agrupamiento::vecino_mas_cercano(otro, &entre_grupos, &activo);
            }
        }
        if paso == 0 {
            narracion.agregar(
                Pizarra::nueva().variable("paso", paso).paso(
                    "vecino[otro] = vecino_mas_cercano(otro, &entre_grupos, &activo);",
                    "Solo los grupos cuyo vecino era alguno de los dos que se unieron necesitan buscar un vecino nuevo; los demás siguen igual. Así cada vuelta es rápida.",
                ),
            );
        }
    }
    let ultima = uniones
        .last()
        .map(|union| union.distancia)
        .unwrap_or_default();
    let ultimo = lista_de(
        Pizarra::nueva().variable("cantidad", cantidad),
        "Uniones hechas",
        &registro,
        &[registro.len() - 1],
    )
    .paso(
        "Ok(uniones)",
        format!(
            "Después de {} uniones todo quedó en un solo grupo. La última unión ocurrió a distancia {}: es la altura de la raíz del árbol. «Qué significa» te enseña a leer el árbol y a cortarlo en grupos.",
            miles(cantidad - 1),
            numero(ultima, 3)
        ),
    );
    (narracion.terminar(ultimo), uniones)
}

fn pasos_cast(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let tabla = leer_tabla(entrada, "datos", "p")?;
    let umbral = entrada.decimal("umbral")?;
    if umbral <= 0.0 || umbral > 1.0 {
        return Err(
            "El umbral θ es una correlación: debe estar entre 0 (sin incluir) y 1.".to_string(),
        );
    }
    exigir_filas(&tabla, 2)?;
    limitar_matriz(&tabla)?;
    Ok(narrar_cast(&tabla, umbral).0)
}

struct EscenaCast<'a> {
    tabla: &'a Tabla,
    similitudes: &'a [Vec<f64>],
    umbral: f64,
}

fn afinidades_en_texto(escena: &EscenaCast, afinidad: &[f64]) -> Vec<String> {
    afinidad
        .iter()
        .enumerate()
        .map(|(elemento, valor)| {
            format!("{}: {}", escena.tabla.nombres[elemento], numero(*valor, 2))
        })
        .collect()
}

fn formar_grupo_narrado(
    escena: &EscenaCast,
    libre: &[bool],
    comienzo: usize,
    eventos_a_narrar: usize,
    narracion: &mut Narracion,
) -> Vec<usize> {
    let similitudes = escena.similitudes;
    let umbral = escena.umbral;
    let cantidad = similitudes.len();
    let mut dentro = vec![false; cantidad];
    let mut afinidad = vec![0.0; cantidad];
    dentro[comienzo] = true;
    agrupamiento::sumar_similitudes(&mut afinidad, similitudes, comienzo, 1.0);
    let mut miembros = 1;
    let mut cambios_permitidos = 2 * libre.iter().filter(|&&esta| esta).count();
    let mut eventos = 0;
    let integrantes = |dentro: &[bool]| -> Vec<usize> {
        (0..cantidad).filter(|&elemento| dentro[elemento]).collect()
    };
    if eventos_a_narrar > 0 {
        narracion.agregar(
            lista_de(
                Pizarra::nueva()
                    .variable("umbral", numero(umbral, 2))
                    .variable("comienzo", comienzo)
                    .variable("miembros", miembros),
                "Afinidad con el grupo (suma de correlaciones con sus miembros)",
                &afinidades_en_texto(escena, &afinidad),
                &[comienzo],
            )
            .paso(
                "dentro[comienzo] = true;",
                format!(
                    "Un grupo nuevo empieza con {}. La afinidad de cada fila es la suma de su correlación con los miembros del grupo; por ahora, solo con {}.",
                    escena.tabla.nombres[comienzo], escena.tabla.nombres[comienzo]
                ),
            ),
        );
    }
    while cambios_permitidos > 0 {
        cambios_permitidos -= 1;
        if let Some(cercano) = agrupamiento::mas_parecido_afuera(&afinidad, &dentro, libre)
            && afinidad[cercano] / miembros as f64 >= umbral
        {
            let media = afinidad[cercano] / miembros as f64;
            dentro[cercano] = true;
            miembros += 1;
            agrupamiento::sumar_similitudes(&mut afinidad, similitudes, cercano, 1.0);
            if eventos < eventos_a_narrar {
                narracion.agregar(
                    lista_de(
                        Pizarra::nueva()
                            .variable("umbral", numero(umbral, 2))
                            .variable("comienzo", comienzo)
                            .variable("miembros", miembros)
                            .variable("cercano", cercano),
                        "Afinidad con el grupo (suma de correlaciones con sus miembros)",
                        &afinidades_en_texto(escena, &afinidad),
                        &[cercano],
                    )
                    .paso(
                        "dentro[cercano] = true;",
                        format!(
                            "{} es la fila libre más parecida al grupo: su correlación media con él es {}, que llega al umbral {}. Entra, y el grupo es {}.",
                            escena.tabla.nombres[cercano],
                            numero(media, 2),
                            numero(umbral, 2),
                            nombres_de(escena.tabla, &integrantes(&dentro))
                        ),
                    ),
                );
            }
            eventos += 1;
            continue;
        }
        if let Some(lejano) = agrupamiento::menos_parecido_adentro(&afinidad, &dentro)
            && miembros > 1
            && afinidad[lejano] / (miembros as f64) < umbral
        {
            let media = afinidad[lejano] / miembros as f64;
            dentro[lejano] = false;
            miembros -= 1;
            agrupamiento::sumar_similitudes(&mut afinidad, similitudes, lejano, -1.0);
            if eventos < eventos_a_narrar {
                narracion.agregar(
                    lista_de(
                        Pizarra::nueva()
                            .variable("umbral", numero(umbral, 2))
                            .variable("comienzo", comienzo)
                            .variable("miembros", miembros)
                            .variable("lejano", lejano),
                        "Afinidad con el grupo (suma de correlaciones con sus miembros)",
                        &afinidades_en_texto(escena, &afinidad),
                        &[lejano],
                    )
                    .paso(
                        "dentro[lejano] = false;",
                        format!(
                            "{} ya no encaja: su correlación media con el grupo bajó a {}, menos que el umbral {}. Sale del grupo, que ahora es {}.",
                            escena.tabla.nombres[lejano],
                            numero(media, 2),
                            numero(umbral, 2),
                            nombres_de(escena.tabla, &integrantes(&dentro))
                        ),
                    ),
                );
            }
            eventos += 1;
            continue;
        }
        if eventos_a_narrar > 0 {
            let candidato = agrupamiento::mas_parecido_afuera(&afinidad, &dentro, libre);
            let explicacion = match candidato {
                Some(otro) => format!(
                    "Nadie más cumple: la fila libre más parecida, {}, tiene correlación media {} con el grupo y no llega al umbral {}, y ningún miembro está por debajo. El grupo queda cerrado: {}.",
                    escena.tabla.nombres[otro],
                    numero(afinidad[otro] / miembros as f64, 2),
                    numero(umbral, 2),
                    nombres_de(escena.tabla, &integrantes(&dentro))
                ),
                None => format!(
                    "No quedan filas libres afuera. El grupo queda cerrado: {}.",
                    nombres_de(escena.tabla, &integrantes(&dentro))
                ),
            };
            narracion.agregar(
                Pizarra::nueva()
                    .variable("umbral", numero(umbral, 2))
                    .variable("comienzo", comienzo)
                    .variable("miembros", miembros)
                    .paso("break;", explicacion),
            );
        }
        break;
    }
    integrantes(&dentro)
}

fn narrar_cast(tabla: &Tabla, umbral: f64) -> (Vec<Paso>, Vec<Vec<usize>>) {
    let similitudes = agrupamiento::matriz_de_correlaciones(&tabla.filas);
    let escena = EscenaCast {
        tabla,
        similitudes: &similitudes,
        umbral,
    };
    let cantidad = similitudes.len();
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        tabla_de_numeros(
            Pizarra::nueva()
                .variable("cantidad", cantidad)
                .variable("umbral", numero(umbral, 2)),
            "Correlación entre cada par de filas",
            &tabla.nombres,
            &tabla.nombres,
            &similitudes,
            2,
            &[],
        )
        .paso(
            "let mut libre = vec![true; cantidad];",
            format!(
                "Primero se calcula la correlación entre cada par de {} filas: 1 es misma forma, 0 sin relación, −1 formas opuestas. Dos filas son «parecidas» si su correlación llega al umbral θ = {}. Todas las filas empiezan libres (sin grupo).",
                miles(cantidad),
                numero(umbral, 2)
            ),
        ),
    );
    let mut libre = vec![true; cantidad];
    let mut grado: Vec<usize> = (0..cantidad)
        .map(|elemento| agrupamiento::vecinos_parecidos(elemento, &similitudes, umbral, &libre))
        .collect();
    let grados: Vec<String> = grado
        .iter()
        .enumerate()
        .map(|(elemento, vecinos)| {
            format!(
                "{}: {} {}",
                tabla.nombres[elemento],
                vecinos,
                segun_cantidad(*vecinos, "parecida", "parecidas")
            )
        })
        .collect();
    let mas_conectado = agrupamiento::libre_de_mayor_grado(&grado, &libre);
    narracion.agregar(
        lista_de(
            Pizarra::nueva().variable("umbral", numero(umbral, 2)),
            "Filas parecidas a cada fila (su grado)",
            &grados,
            &[mas_conectado],
        )
        .paso(
            ".map(|elemento| vecinos_parecidos(elemento, similitudes, umbral, &libre))",
            format!(
                "Para cada fila contamos cuántas otras le son parecidas (su grado). {} tiene el mayor ({}): es la que tiene más «amigas», y con ella empezará el primer grupo.",
                tabla.nombres[mas_conectado], grado[mas_conectado]
            ),
        ),
    );
    let mut grupos: Vec<Vec<usize>> = Vec::new();
    let mut filas_de_grupos: Vec<Vec<String>> = Vec::new();
    let mut quedan = cantidad;
    let salto = salto_de(cantidad);
    while quedan > 0 {
        let orden = grupos.len();
        let comienzo = agrupamiento::libre_de_mayor_grado(&grado, &libre);
        let eventos = if orden == 0 {
            EVENTOS_DEL_PRIMER_GRUPO
        } else if orden < PASOS_DETALLADOS {
            EVENTOS_DE_LOS_SIGUIENTES
        } else {
            0
        };
        if orden < PASOS_DETALLADOS {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("quedan", quedan)
                    .variable("comienzo", comienzo)
                    .paso(
                        "let comienzo = libre_de_mayor_grado(&grado, &libre);",
                        format!(
                            "Quedan {} filas sin grupo. Para el grupo {} empezamos con la fila libre de mayor grado: {} ({}).",
                            miles(quedan),
                            orden + 1,
                            tabla.nombres[comienzo],
                            grado[comienzo]
                        ),
                    ),
            );
        }
        let grupo = formar_grupo_narrado(&escena, &libre, comienzo, eventos, &mut narracion);
        for &elemento in &grupo {
            libre[elemento] = false;
            quedan -= 1;
        }
        for &elemento in &grupo {
            for otro in 0..cantidad {
                if libre[otro] && similitudes[otro][elemento] >= umbral {
                    grado[otro] -= 1;
                }
            }
        }
        if orden == 0 {
            narracion.agregar(
                lista_de(
                    Pizarra::nueva().variable("quedan", quedan),
                    "Filas parecidas que siguen libres (su grado)",
                    &(0..cantidad)
                        .map(|elemento| {
                            if libre[elemento] {
                                format!("{}: {}", tabla.nombres[elemento], grado[elemento])
                            } else {
                                format!("{}: ya tiene grupo", tabla.nombres[elemento])
                            }
                        })
                        .collect::<Vec<String>>(),
                    &grupo,
                )
                .paso(
                    "grado[otro] -= 1;",
                    "Las filas del grupo ya no están libres. Cada fila libre parecida a alguna de ellas pierde una «amiga» disponible, así que su grado baja. Los grados que ves son los nuevos.",
                ),
            );
        }
        filas_de_grupos.push(vec![
            format!("{}", orden + 1),
            miles(grupo.len()),
            nombres_de(tabla, &grupo),
        ]);
        let miembros_del_grupo = grupo.len();
        grupos.push(grupo);
        if toca_narrar(orden, salto) {
            narracion.agregar(
                tabla_de_textos(
                    Pizarra::nueva()
                        .variable("quedan", quedan)
                        .variable("umbral", numero(umbral, 2)),
                    "Grupos formados",
                    &["grupo", "filas", "miembros"],
                    &filas_de_grupos,
                    Some(orden),
                )
                .paso(
                    "grupos.push(grupo);",
                    format!(
                        "El grupo {} queda con {} {}. Esas filas ya no están libres y las parecidas a ellas pierden una amiga disponible. Quedan {} filas sin grupo.",
                        orden + 1,
                        miembros_del_grupo,
                        segun_cantidad(miembros_del_grupo, "fila", "filas"),
                        miles(quedan)
                    ),
                ),
            );
        }
    }
    let solitarios = grupos.iter().filter(|grupo| grupo.len() == 1).count();
    let ultimo = tabla_de_textos(
        Pizarra::nueva().variable("umbral", numero(umbral, 2)),
        "Grupos formados",
        &["grupo", "filas", "miembros"],
        &filas_de_grupos,
        None,
    )
    .paso(
        "Ok(grupos)",
        format!(
            "Listo: {} {} con umbral {}, sin haber elegido de antemano cuántos habría. {} {} quedó {} solo. Abajo, «Qué significa» explica cómo leerlos.",
            miles(grupos.len()),
            segun_cantidad(grupos.len(), "grupo", "grupos"),
            numero(umbral, 2),
            miles(solitarios),
            segun_cantidad(solitarios, "fila", "filas"),
            segun_cantidad(solitarios, "quedó", "quedaron")
        ),
    );
    (narracion.terminar(ultimo), grupos)
}

fn pasos_componentes(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let tabla = leer_tabla(entrada, "datos", "p")?;
    exigir_filas(&tabla, 2)?;
    if tabla.dimension() < 2 {
        return Err("La tabla tiene una sola columna: no hay nada que proyectar.".to_string());
    }
    Ok(narrar_componentes(&tabla)?.0)
}

fn potencia_narrada(
    matriz: &[Vec<f64>],
    distinto_de: &[f64],
    columnas: &[String],
    narracion: &mut Narracion,
) -> Vec<f64> {
    let dimension = matriz.len();
    let es_primera = distinto_de.is_empty();
    let escala: f64 = (0..dimension).map(|eje| matriz[eje][eje].abs()).sum();
    let mut azar = Azar::nuevo(SEMILLA_DE_POTENCIA);
    let mut vector: Vec<f64> = (0..dimension).map(|_| azar.normal()).collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("escala", numero(escala, 3))
            .variable("vector", vector_corto(&vector))
            .paso(
                "let mut vector: Vec<f64> = (0..dimension).map(|_| azar.normal()).collect();",
                "El método de la potencia empieza con una dirección cualquiera, sacada del azar (con una semilla fija, así siempre sale lo mismo). Luego la «empuja» una y otra vez con la matriz de covarianza: la dirección en que los datos más varían es la que sobrevive.",
            ),
    );
    agrupamiento::quitar_direccion(&mut vector, distinto_de);
    if !es_primera {
        narracion.agregar(
            Pizarra::nueva()
                .variable("vector", vector_corto(&vector))
                .paso(
                    "quitar_direccion(&mut vector, distinto_de);",
                    "Para la segunda dirección le quitamos al vector la parte que apunta hacia la primera. Así la segunda queda perpendicular a la primera y no repite lo mismo.",
                ),
        );
    }
    if agrupamiento::normalizar(&mut vector) <= PARTE_DESPRECIABLE {
        return vec![0.0; dimension];
    }
    let detalle = if es_primera { PASOS_DETALLADOS } else { 1 };
    let mut vueltas = 0;
    let mut estabilizo = false;
    let mut proxima_pausa = 10;
    for vuelta in 0..VUELTAS_DE_POTENCIA {
        let mut siguiente = agrupamiento::multiplicar(matriz, &vector);
        agrupamiento::quitar_direccion(&mut siguiente, distinto_de);
        if agrupamiento::normalizar(&mut siguiente) <= escala * PARTE_DESPRECIABLE {
            break;
        }
        let cambio = agrupamiento::distancia_al_cuadrado(&siguiente, &vector);
        vector = siguiente;
        vueltas = vuelta + 1;
        let toca = vuelta < detalle || (es_primera && vueltas == proxima_pausa);
        if vueltas == proxima_pausa {
            proxima_pausa *= 3;
        }
        if toca {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("escala", numero(escala, 3))
                    .variable("vector", vector_corto(&vector))
                    .variable("cambio", cientifico(cambio))
                    .paso(
                        "let mut siguiente = multiplicar(matriz, &vector);",
                        format!(
                            "Vuelta {}: multiplicamos el vector por la matriz y lo dejamos de largo 1. El vector se movió {} (la distancia al cuadrado entre el de antes y el de ahora): {}.",
                            miles(vueltas),
                            cientifico(cambio),
                            if vuelta == 0 {
                                "si cambia poco, ya casi apunta a la dirección buscada"
                            } else {
                                "cada vuelta lo acerca más a la dirección buscada"
                            }
                        ),
                    ),
            );
        }
        if cambio < CAMBIO_DESPRECIABLE {
            estabilizo = true;
            break;
        }
    }
    orientar_narrado(&mut vector, columnas, vueltas, estabilizo, narracion);
    vector
}

fn orientar_narrado(
    vector: &mut [f64],
    columnas: &[String],
    vueltas: usize,
    estabilizo: bool,
    narracion: &mut Narracion,
) {
    agrupamiento::orientar(vector);
    let pesos: Vec<Vec<f64>> = vector.iter().map(|&peso| vec![peso]).collect();
    narracion.agregar(
        tabla_de_numeros(
            Pizarra::nueva().variable("vector", vector_corto(vector)),
            "Peso de cada columna en esta dirección",
            columnas,
            &["peso".to_string()],
            &pesos,
            3,
            &[],
        )
        .paso(
            "orientar(&mut vector);",
            format!(
                "{} tras {} {}, el vector ya no cambia y esa es la dirección. Solo falta darle un signo fijo para que el resultado siempre salga igual. Cada peso dice cuánto cuenta esa columna: más grande en valor absoluto, más influye.",
                if estabilizo { "Listo:" } else { "Se detuvo:" },
                miles(vueltas),
                segun_cantidad(vueltas, "vuelta", "vueltas")
            ),
        ),
    );
}

fn narrar_componentes(tabla: &Tabla) -> Result<(Vec<Paso>, Componentes), String> {
    let puntos = &tabla.filas;
    let dimension = agrupamiento::revisar_puntos(puntos)?;
    let cantidad = puntos.len();
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        tabla_de_numeros(
            Pizarra::nueva().variable("dimension", dimension),
            "Tus datos",
            &tabla.nombres,
            &tabla.columnas,
            puntos,
            2,
            &[],
        )
        .paso(
            "let dimension = revisar_puntos(puntos)?;",
            format!(
                "Tienes {} filas con {} números cada una: para dibujarlas harían falta {} ejes. La meta es encontrar las 2 direcciones en que los datos más varían y usarlas como ejes de un plano donde sí se pueda dibujar, perdiendo lo menos posible.",
                miles(cantidad),
                miles(dimension),
                miles(dimension)
            ),
        ),
    );
    let media = agrupamiento::centro_de_gravedad(puntos, dimension);
    narracion.agregar(
        lista_de(
            Pizarra::nueva()
                .variable("dimension", dimension)
                .variable("media", vector_corto(&media)),
            "Promedio de cada columna",
            &tabla
                .columnas
                .iter()
                .zip(&media)
                .map(|(columna, promedio)| format!("{columna}: {}", numero(*promedio, 2)))
                .collect::<Vec<String>>(),
            &[],
        )
        .paso(
            "let media = centro_de_gravedad(puntos, dimension);",
            "Primero se calcula el promedio de cada columna: es el centro de la nube de puntos. Las direcciones se miden a partir de ese centro.",
        ),
    );
    let covarianza = agrupamiento::matriz_de_covarianza(puntos, &media);
    narracion.agregar(
        tabla_de_numeros(
            Pizarra::nueva().variable("dimension", dimension),
            "Covarianza entre columnas",
            &tabla.columnas,
            &tabla.columnas,
            &covarianza,
            2,
            &[],
        )
        .paso(
            "let covarianza = matriz_de_covarianza(puntos, &media);",
            "La matriz de covarianza dice, para cada par de columnas, si suben y bajan juntas (número positivo grande), al revés (negativo) o sin relación (cerca de 0). Es la información que resume cómo se dispersan los datos.",
        ),
    );
    let varianza_total: f64 = (0..dimension).map(|eje| covarianza[eje][eje]).sum();
    let diagonal: Vec<(usize, usize)> = (0..dimension).map(|eje| (eje, eje)).collect();
    narracion.agregar(
        tabla_de_numeros(
            Pizarra::nueva()
                .variable("dimension", dimension)
                .variable("varianza_total", numero(varianza_total, 3)),
            "Covarianza entre columnas",
            &tabla.columnas,
            &tabla.columnas,
            &covarianza,
            2,
            &diagonal,
        )
        .paso(
            "let varianza_total = (0..dimension).map(|eje| covarianza[eje][eje]).sum();",
            format!(
                "La diagonal es la varianza de cada columna (cuánto se dispersa por sí sola). Su suma, {}, es toda la variación que hay en los datos: la que queremos conservar lo más posible en el dibujo.",
                numero(varianza_total, 3)
            ),
        ),
    );
    narracion.agregar(
        Pizarra::nueva().variable("varianza_total", numero(varianza_total, 3)).paso(
            "let primera = vector_propio_dominante(&covarianza, &[]);",
            "Ahora buscamos la primera dirección: la que más varianza captura. Para eso se usa el método de la potencia.",
        ),
    );
    let primera = potencia_narrada(&covarianza, &[], &tabla.columnas, &mut narracion);
    narracion.agregar(
        Pizarra::nueva().variable("varianza_total", numero(varianza_total, 3)).paso(
            "let segunda = vector_propio_dominante(&covarianza, &primera);",
            "La segunda dirección se busca igual, pero obligándola a ser perpendicular a la primera.",
        ),
    );
    let segunda = potencia_narrada(&covarianza, &primera, &tabla.columnas, &mut narracion);
    let varianza_primera = agrupamiento::cociente_de_rayleigh(&covarianza, &primera);
    let varianza_segunda = agrupamiento::cociente_de_rayleigh(&covarianza, &segunda);
    let parte = |varianza: f64| {
        if varianza_total > 0.0 {
            porcentaje(varianza / varianza_total, 0)
        } else {
            porcentaje(0.0, 0)
        }
    };
    narracion.agregar(
        Pizarra::nueva()
            .variable("varianza_total", numero(varianza_total, 3))
            .variable("varianza_primera", numero(varianza_primera, 3))
            .paso(
                "let varianza_primera = cociente_de_rayleigh(&covarianza, &primera);",
                format!(
                    "Cuánta variación captura la primera dirección: {} de {}, o sea {} del total.",
                    numero(varianza_primera, 3),
                    numero(varianza_total, 3),
                    parte(varianza_primera)
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("varianza_total", numero(varianza_total, 3))
            .variable("varianza_segunda", numero(varianza_segunda, 3))
            .paso(
                "let varianza_segunda = cociente_de_rayleigh(&covarianza, &segunda);",
                format!(
                    "La segunda captura {}, o sea {} del total. Juntas conservan {} de la variación.",
                    numero(varianza_segunda, 3),
                    parte(varianza_segunda),
                    parte(varianza_primera + varianza_segunda)
                ),
            ),
    );
    let mut componentes = Componentes {
        media,
        primera,
        segunda,
        varianza_primera,
        varianza_segunda,
        varianza_total,
        coordenadas: Vec::with_capacity(cantidad),
    };
    let primer_punto = &puntos[0];
    let coordenada_primera = agrupamiento::proyectar(primer_punto, &componentes);
    let mut productos: Vec<Vec<String>> = (0..dimension.min(VALORES_A_LA_VISTA))
        .map(|columna| {
            vec![
                tabla.columnas[columna].clone(),
                numero(primer_punto[columna] - componentes.media[columna], 2),
                numero(componentes.primera[columna], 3),
                numero(componentes.segunda[columna], 3),
            ]
        })
        .collect();
    productos.push(vec![
        "coordenadas".to_string(),
        String::new(),
        numero(coordenada_primera[0], 3),
        numero(coordenada_primera[1], 3),
    ]);
    let ultima_fila = productos.len() - 1;
    narracion.agregar(
        tabla_de_textos(
            Pizarra::nueva().variable("coordenada", vector_corto(&coordenada_primera)),
            &format!("Proyección de {}", tabla.nombres[0]),
            &[
                "columna",
                "valor − promedio",
                "peso en la dirección 1",
                "peso en la dirección 2",
            ],
            &productos,
            Some(ultima_fila),
        )
        .paso(
            "let coordenada = proyectar(punto, &componentes);",
            format!(
                "Para dibujar {} le restamos el promedio y multiplicamos, columna por columna, por los pesos de cada dirección, y sumamos. Sus coordenadas en el plano son ({}, {}).",
                tabla.nombres[0],
                numero(coordenada_primera[0], 3),
                numero(coordenada_primera[1], 3)
            ),
        ),
    );
    for punto in puntos {
        let coordenada = agrupamiento::proyectar(punto, &componentes);
        componentes.coordenadas.push(coordenada);
    }
    let lugares: Vec<String> = componentes
        .coordenadas
        .iter()
        .enumerate()
        .map(|(indice, lugar)| {
            format!(
                "{}: ({}, {})",
                tabla.nombres[indice],
                numero(lugar[0], 2),
                numero(lugar[1], 2)
            )
        })
        .collect();
    narracion.agregar(
        lista_de(
            Pizarra::nueva().variable("dimension", dimension),
            "Coordenadas en el plano",
            &lugares,
            &[0],
        )
        .paso(
            "for punto in puntos {",
            format!(
                "Lo mismo se hace con las {} filas: cada una pasa de {} números a solo 2 coordenadas. Filas con perfiles parecidos quedan cerca en el plano.",
                miles(cantidad),
                miles(dimension)
            ),
        ),
    );
    let conservado = if varianza_total > 0.0 {
        (varianza_primera + varianza_segunda) / varianza_total
    } else {
        0.0
    };
    let ultimo = Pizarra::nueva()
        .variable("varianza_primera", numero(varianza_primera, 3))
        .variable("varianza_segunda", numero(varianza_segunda, 3))
        .paso(
            "Ok(componentes)",
            format!(
                "Listo: el dibujo en 2D conserva {} de la variación de tus datos. {} «Qué significa» te enseña a leer el dibujo.",
                porcentaje(conservado, 0),
                if conservado >= VARIANZA_FIEL {
                    "Es bastante: lo que veas en el plano es una buena imagen de los datos."
                } else if conservado >= VARIANZA_POBRE {
                    "Es una parte razonable, pero se perdió algo de información."
                } else {
                    "Es poco: el plano es solo una imagen parcial de los datos."
                }
            ),
        );
    Ok((narracion.terminar(ultimo), componentes))
}

#[cfg(test)]
mod pruebas_del_paso_a_paso {
    use super::*;

    fn ejemplo_de(id: &str) -> Entrada {
        Entrada::desde_json(crate::ejemplo(id).unwrap()).unwrap()
    }

    fn tabla_de(id: &str) -> Tabla {
        leer_tabla(&ejemplo_de(id), "datos", "p").unwrap()
    }

    fn pasos_de(id: &str) -> Vec<Paso> {
        crate::pasos(id, crate::ejemplo(id).unwrap()).unwrap()
    }

    #[test]
    fn la_narracion_de_distancias_coincide_con_el_nucleo() {
        let tabla = tabla_de("agrupamiento.distancias");
        let (pasos, primeras) = narrar_distancias(&tabla);
        let euclidea = agrupamiento::matriz_de_distancias(&tabla.filas, Medida::Euclidea);
        let pearson = agrupamiento::matriz_de_distancias(&tabla.filas, Medida::Pearson);
        assert_eq!(primeras[0], euclidea[0][1]);
        assert_eq!(primeras[1], pearson[0][1]);
        assert!(pasos.len() > 10 && pasos.len() <= MAXIMO_PASOS);
    }

    #[test]
    fn la_narracion_del_mas_lejano_coincide_con_el_nucleo() {
        let id = "agrupamiento.mas_lejano_primero";
        let tabla = tabla_de(id);
        let k = leer_k(&ejemplo_de(id), "k", tabla.cantidad()).unwrap();
        let (pasos, centros) = narrar_mas_lejano(&tabla, k);
        assert_eq!(
            centros,
            agrupamiento::mas_lejano_primero(&tabla.filas, k).unwrap()
        );
        assert!(pasos.len() > 5);
    }

    #[test]
    fn la_narracion_de_la_distorsion_coincide_con_el_nucleo() {
        let id = "agrupamiento.distorsion";
        let entrada = ejemplo_de(id);
        let tabla = leer_tabla(&entrada, "datos", "p").unwrap();
        let centros = leer_tabla(&entrada, "centros", "centro ").unwrap();
        let (pasos, valor) = narrar_distorsion(&tabla, &centros);
        assert_eq!(
            valor,
            agrupamiento::distorsion_cuadratica(&tabla.filas, &centros.filas)
        );
        assert!(pasos.len() > 5);
    }

    #[test]
    fn la_narracion_del_codo_coincide_con_el_nucleo() {
        let id = "agrupamiento.codo";
        let entrada = ejemplo_de(id);
        let tabla = leer_tabla(&entrada, "datos", "p").unwrap();
        let k_maximo = entrada.natural("k_maximo").unwrap();
        let semilla = leer_semilla(&entrada).unwrap();
        let (pasos, curva) = narrar_codo(&tabla, k_maximo, semilla).unwrap();
        assert_eq!(
            curva,
            agrupamiento::curva_del_codo(&tabla.filas, k_maximo, semilla).unwrap()
        );
        assert!(pasos.len() > 5);
    }

    #[test]
    fn la_narracion_de_k_medias_suaves_coincide_con_el_nucleo() {
        let id = "agrupamiento.k_medias_suaves";
        let entrada = ejemplo_de(id);
        let tabla = leer_tabla(&entrada, "datos", "p").unwrap();
        let k = leer_k(&entrada, "k", tabla.cantidad()).unwrap();
        let rigidez = entrada.decimal("rigidez").unwrap();
        let iteraciones = leer_iteraciones(&entrada, "iteraciones").unwrap();
        let (pasos, resultado) = narrar_suaves(&tabla, k, rigidez, iteraciones);
        assert_eq!(
            resultado,
            agrupamiento::k_medias_suaves(&tabla.filas, k, rigidez, iteraciones).unwrap()
        );
        assert!(pasos.len() > 5);
    }

    #[test]
    fn la_narracion_jerarquica_coincide_con_el_nucleo() {
        let id = "agrupamiento.jerarquico";
        let tabla = tabla_de(id);
        let distancias = agrupamiento::matriz_de_distancias(&tabla.filas, Medida::Pearson);
        let (pasos, uniones) = narrar_jerarquico(&tabla, &distancias, MEDIDA_PEARSON);
        assert_eq!(
            uniones,
            agrupamiento::agrupamiento_jerarquico(&distancias).unwrap()
        );
        assert!(pasos.len() > 10);
    }

    #[test]
    fn la_narracion_de_cast_coincide_con_el_nucleo() {
        let id = "agrupamiento.cast";
        let entrada = ejemplo_de(id);
        let tabla = leer_tabla(&entrada, "datos", "p").unwrap();
        let umbral = entrada.decimal("umbral").unwrap();
        let (pasos, grupos) = narrar_cast(&tabla, umbral);
        let similitudes = agrupamiento::matriz_de_correlaciones(&tabla.filas);
        assert_eq!(grupos, agrupamiento::cast(&similitudes, umbral).unwrap());
        assert!(pasos.len() > 5);
    }

    #[test]
    fn la_narracion_de_componentes_coincide_con_el_nucleo() {
        let tabla = tabla_de("agrupamiento.componentes_principales");
        let (pasos, componentes) = narrar_componentes(&tabla).unwrap();
        assert_eq!(
            componentes,
            agrupamiento::componentes_principales(&tabla.filas).unwrap()
        );
        assert!(pasos.len() > 5);
    }

    #[test]
    fn todas_las_entradas_del_modulo_narran() {
        for algoritmo in algoritmos() {
            assert!(algoritmo.pasos.is_some(), "{}", algoritmo.id);
            let pasos = pasos_de(algoritmo.id);
            assert!(!pasos.is_empty() && pasos.len() <= MAXIMO_PASOS);
        }
    }
}
