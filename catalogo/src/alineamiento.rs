use crate::formato::{cientifico, decimal, miles, porcentaje};
use crate::pizarra::{Marca as Resalte, Narracion, Pizarra};
use crate::{Algoritmo, Arista, Bloque, Campo, Entrada, FilaAlineada, Marca, Paso, Salida, Tipo};
use nucleo::alineamiento::{self, Alineamiento, Puntuacion};
use nucleo::azar::Azar;
use nucleo::matrices::{blosum62, pam250};
use nucleo::secuencia;
use serde_json::{Value, json};
use std::collections::{HashMap, HashSet, VecDeque};

const CAPITULO: &str = "alineamiento";
const PUNTUACIONES: &[&str] = &["BLOSUM62", "PAM250", "ADN +1/-1"];
const LIMITE_CUADRATICO: usize = 25_000_000;
const LIMITE_LINEAL: usize = 400_000_000;
const LIMITE_TRES: usize = 8_000_000;
const LIMITE_DINERO: i64 = 1_000_000;
const LIMITE_MONEDAS: usize = 100;
const LIMITE_ESQUINAS: usize = 1_000_000;
const LIMITE_NODOS: usize = 100_000;
const LIMITE_ARISTAS: usize = 500_000;
const LIMITE_GRAFO_VISIBLE: usize = 60;
const LIMITE_PASO_A_PASO: usize = 15;
const MAXIMO_PASOS: usize = 400;
const BARAJADOS: usize = 5;
const LIMITE_PESO: f64 = 1e12;
const MAXIMO_LETRAS_VISIBLES: usize = 20_000;
const GENOMA_ECOLI: f64 = 4_641_652.0;

const CONSEJO_CUADRATICO: &str = "Para secuencias así de largas usa «Alineamiento en espacio lineal», que no guarda la tabla entera.";
const CONSEJO_LINEAL: &str = "Recorta las secuencias a la región que te interesa.";

const PUNTUACION_BLOSUM: Campo = Campo {
    id: "puntuacion",
    etiqueta: "Puntuación de cada par de letras",
    tipo: Tipo::Opcion(PUNTUACIONES),
    ejemplo: "BLOSUM62",
};
const PUNTUACION_PAM: Campo = Campo {
    id: "puntuacion",
    etiqueta: "Puntuación de cada par de letras",
    tipo: Tipo::Opcion(PUNTUACIONES),
    ejemplo: "PAM250",
};
const PUNTUACION_ADN: Campo = Campo {
    id: "puntuacion",
    etiqueta: "Puntuación de cada par de letras",
    tipo: Tipo::Opcion(PUNTUACIONES),
    ejemplo: "ADN +1/-1",
};
const INDEL_CINCO: Campo = Campo {
    id: "indel",
    etiqueta: "Castigo por cada hueco (indel)",
    tipo: Tipo::Entero,
    ejemplo: "5",
};

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![
        Algoritmo {
            id: "alineamiento.cambio",
            capitulo: CAPITULO,
            nombre: "Cambio con el mínimo de monedas",
            resumen: "Cuántas monedas hacen falta, como mínimo, para pagar una cantidad exacta.",
            campos: &[
                Campo {
                    id: "dinero",
                    etiqueta: "Cantidad a pagar",
                    tipo: Tipo::Entero,
                    ejemplo: "40",
                },
                Campo {
                    id: "monedas",
                    etiqueta: "Valores de las monedas",
                    tipo: Tipo::Numeros,
                    ejemplo: "1 5 10 20 25",
                },
            ],
            funcion: "alineamiento::cambio_minimo",
            ejecutar: cambio,
            pasos: Some(pasos_cambio),
        },
        Algoritmo {
            id: "alineamiento.manhattan",
            capitulo: CAPITULO,
            nombre: "Turista de Manhattan",
            resumen: "El recorrido más valioso por una cuadrícula de calles que solo bajan o van a la derecha.",
            campos: &[
                Campo {
                    id: "abajo",
                    etiqueta: "Pesos de las calles hacia abajo",
                    tipo: Tipo::Matriz,
                    ejemplo: "4 0 3 2 4\n6 8 0 0 3\n0 8 9 5 4\n4 7 5 7 4",
                },
                Campo {
                    id: "derecha",
                    etiqueta: "Pesos de las calles hacia la derecha",
                    tipo: Tipo::Matriz,
                    ejemplo: "5 1 0 0\n1 1 1 2\n9 7 0 4\n2 4 8 8\n1 4 9 8",
                },
            ],
            funcion: "alineamiento::turista_de_manhattan",
            ejecutar: manhattan,
            pasos: Some(pasos_manhattan),
        },
        Algoritmo {
            id: "alineamiento.orden_topologico",
            capitulo: CAPITULO,
            nombre: "Orden topológico",
            resumen: "Ordenar los nodos de un grafo sin ciclos para que cada flecha apunte hacia adelante.",
            campos: &[Campo {
                id: "aristas",
                etiqueta: "Flechas, una por línea (a -> b, o a -> b, c)",
                tipo: Tipo::Texto,
                ejemplo: "cultivar bacterias -> extraer ADN\nextraer ADN -> medir concentración\nextraer ADN -> PCR\ndiseñar cebadores -> PCR\nmedir concentración -> PCR\nPCR -> electroforesis\nPCR -> purificar\npurificar -> secuenciar\nelectroforesis -> secuenciar",
            }],
            funcion: "alineamiento::orden_topologico",
            ejecutar: orden_topologico,
            pasos: Some(pasos_orden_topologico),
        },
        Algoritmo {
            id: "alineamiento.camino_dag",
            capitulo: CAPITULO,
            nombre: "Camino más largo en un grafo sin ciclos",
            resumen: "El camino de mayor peso entre dos nodos, recorriendo los nodos en orden topológico.",
            campos: &[
                Campo {
                    id: "aristas",
                    etiqueta: "Flechas con peso, una por línea (a -> b: peso)",
                    tipo: Tipo::Texto,
                    ejemplo: "inicio -> a: 3\ninicio -> b: 2\ninicio -> c: 5\na -> d: 4\nb -> d: 6\nb -> e: 1\nc -> e: 2\nc -> f: 1\nd -> f: 3\ne -> f: 4\nd -> fin: 2\ne -> fin: 7\nf -> fin: 5",
                },
                Campo {
                    id: "fuente",
                    etiqueta: "Nodo de salida",
                    tipo: Tipo::Texto,
                    ejemplo: "inicio",
                },
                Campo {
                    id: "sumidero",
                    etiqueta: "Nodo de llegada",
                    tipo: Tipo::Texto,
                    ejemplo: "fin",
                },
            ],
            funcion: "alineamiento::camino_mas_largo",
            ejecutar: camino_dag,
            pasos: Some(pasos_camino_dag),
        },
        Algoritmo {
            id: "alineamiento.subsecuencia_comun",
            capitulo: CAPITULO,
            nombre: "Subsecuencia común más larga",
            resumen: "Las letras que aparecen en el mismo orden en las dos secuencias, no necesariamente seguidas.",
            campos: &[
                Campo {
                    id: "v",
                    etiqueta: "Secuencia v",
                    tipo: Tipo::Secuencia,
                    ejemplo: "@gen:tufA",
                },
                Campo {
                    id: "w",
                    etiqueta: "Secuencia w",
                    tipo: Tipo::Secuencia,
                    ejemplo: "@gen:tufB",
                },
                Campo {
                    id: "semilla",
                    etiqueta: "Semilla para barajar",
                    tipo: Tipo::Entero,
                    ejemplo: "1",
                },
            ],
            funcion: "alineamiento::subsecuencia_comun_mas_larga",
            ejecutar: subsecuencia_comun,
            pasos: Some(pasos_subsecuencia_comun),
        },
        Algoritmo {
            id: "alineamiento.global",
            capitulo: CAPITULO,
            nombre: "Alineamiento global",
            resumen: "Alinear dos secuencias de punta a punta con el mejor puntaje posible.",
            campos: &[
                Campo {
                    id: "v",
                    etiqueta: "Secuencia v",
                    tipo: Tipo::Proteina,
                    ejemplo: ">GyrA 75-88\nGKYHPHGDSAVYDT",
                },
                Campo {
                    id: "w",
                    etiqueta: "Secuencia w",
                    tipo: Tipo::Proteina,
                    ejemplo: ">ParC 72-84\nGKYHPHGDSACYE",
                },
                PUNTUACION_BLOSUM,
                INDEL_CINCO,
            ],
            funcion: "alineamiento::alineamiento_global",
            ejecutar: global,
            pasos: Some(pasos_global),
        },
        Algoritmo {
            id: "alineamiento.local",
            capitulo: CAPITULO,
            nombre: "Alineamiento local",
            resumen: "El tramo de cada secuencia que mejor se parece al de la otra, ignorando el resto.",
            campos: &[
                Campo {
                    id: "v",
                    etiqueta: "Secuencia v",
                    tipo: Tipo::Proteina,
                    ejemplo: "@proteina:gyrA",
                },
                Campo {
                    id: "w",
                    etiqueta: "Secuencia w",
                    tipo: Tipo::Proteina,
                    ejemplo: "@proteina:parC",
                },
                PUNTUACION_PAM,
                INDEL_CINCO,
            ],
            funcion: "alineamiento::alineamiento_local",
            ejecutar: local,
            pasos: Some(pasos_local),
        },
        Algoritmo {
            id: "alineamiento.edicion",
            capitulo: CAPITULO,
            nombre: "Distancia de edición",
            resumen: "El menor número de sustituciones, inserciones y borrados que convierten una secuencia en la otra.",
            campos: &[
                Campo {
                    id: "v",
                    etiqueta: "Secuencia v",
                    tipo: Tipo::Secuencia,
                    ejemplo: "@gen:tufA",
                },
                Campo {
                    id: "w",
                    etiqueta: "Secuencia w",
                    tipo: Tipo::Secuencia,
                    ejemplo: "@gen:tufB",
                },
            ],
            funcion: "alineamiento::distancia_de_edicion",
            ejecutar: edicion,
            pasos: Some(pasos_edicion),
        },
        Algoritmo {
            id: "alineamiento.ajustado",
            capitulo: CAPITULO,
            nombre: "Alineamiento ajustado",
            resumen: "Dónde encaja mejor una secuencia corta, entera, dentro de una larga.",
            campos: &[
                Campo {
                    id: "v",
                    etiqueta: "Secuencia larga v",
                    tipo: Tipo::Secuencia,
                    ejemplo: "@ecoli:3925700-3926000",
                },
                Campo {
                    id: "w",
                    etiqueta: "Secuencia corta w",
                    tipo: Tipo::Secuencia,
                    ejemplo: "TTATCCACA",
                },
                PUNTUACION_ADN,
                Campo {
                    id: "indel",
                    etiqueta: "Castigo por cada hueco (indel)",
                    tipo: Tipo::Entero,
                    ejemplo: "1",
                },
            ],
            funcion: "alineamiento::alineamiento_ajustado",
            ejecutar: ajustado,
            pasos: Some(pasos_ajustado),
        },
        Algoritmo {
            id: "alineamiento.solapado",
            capitulo: CAPITULO,
            nombre: "Alineamiento solapado",
            resumen: "Cuánto se superponen el final de una secuencia y el comienzo de otra, como dos lecturas vecinas.",
            campos: &[
                Campo {
                    id: "v",
                    etiqueta: "Secuencia v (su final)",
                    tipo: Tipo::Secuencia,
                    ejemplo: "@ecoli:3883000-3883080",
                },
                Campo {
                    id: "w",
                    etiqueta: "Secuencia w (su comienzo)",
                    tipo: Tipo::Secuencia,
                    ejemplo: "GGTAGTAGCATTTAAACTCTTCGATCGCGTTGTTTTGCAGGGCTTTAACCATGTCCTGAACAAAGCGCTCGGAGTGCATA",
                },
                PUNTUACION_ADN,
                Campo {
                    id: "indel",
                    etiqueta: "Castigo por cada hueco (indel)",
                    tipo: Tipo::Entero,
                    ejemplo: "2",
                },
            ],
            funcion: "alineamiento::alineamiento_solapado",
            ejecutar: solapado,
            pasos: Some(pasos_solapado),
        },
        Algoritmo {
            id: "alineamiento.afin",
            capitulo: CAPITULO,
            nombre: "Alineamiento con huecos afines",
            resumen: "Alineamiento global donde abrir un hueco cuesta más que alargarlo, con tres tablas a la vez.",
            campos: &[
                Campo {
                    id: "v",
                    etiqueta: "Secuencia v",
                    tipo: Tipo::Proteina,
                    ejemplo: "@proteina:gyrB",
                },
                Campo {
                    id: "w",
                    etiqueta: "Secuencia w",
                    tipo: Tipo::Proteina,
                    ejemplo: "@proteina:parE",
                },
                PUNTUACION_BLOSUM,
                Campo {
                    id: "apertura",
                    etiqueta: "Castigo por abrir un hueco",
                    tipo: Tipo::Entero,
                    ejemplo: "11",
                },
                Campo {
                    id: "extension",
                    etiqueta: "Castigo por alargarlo una letra más",
                    tipo: Tipo::Entero,
                    ejemplo: "1",
                },
            ],
            funcion: "alineamiento::alineamiento_afin",
            ejecutar: afin,
            pasos: Some(pasos_afin),
        },
        Algoritmo {
            id: "alineamiento.arista_central",
            capitulo: CAPITULO,
            nombre: "Arista central",
            resumen: "Por dónde cruza el alineamiento óptimo la columna del medio, calculado sin guardar la tabla.",
            campos: &[
                Campo {
                    id: "v",
                    etiqueta: "Secuencia v",
                    tipo: Tipo::Proteina,
                    ejemplo: "@proteina:gyrA",
                },
                Campo {
                    id: "w",
                    etiqueta: "Secuencia w",
                    tipo: Tipo::Proteina,
                    ejemplo: "@proteina:parC",
                },
                PUNTUACION_BLOSUM,
                INDEL_CINCO,
            ],
            funcion: "alineamiento::arista_central",
            ejecutar: arista_central,
            pasos: Some(pasos_arista_central),
        },
        Algoritmo {
            id: "alineamiento.espacio_lineal",
            capitulo: CAPITULO,
            nombre: "Alineamiento en espacio lineal",
            resumen: "El mismo alineamiento global, partiendo el problema por la arista central para usar poca memoria.",
            campos: &[
                Campo {
                    id: "v",
                    etiqueta: "Secuencia v",
                    tipo: Tipo::Proteina,
                    ejemplo: "@proteina:lacZ",
                },
                Campo {
                    id: "w",
                    etiqueta: "Secuencia w",
                    tipo: Tipo::Proteina,
                    ejemplo: "@proteina:ebgA",
                },
                PUNTUACION_BLOSUM,
                INDEL_CINCO,
            ],
            funcion: "alineamiento::alineamiento_espacio_lineal",
            ejecutar: espacio_lineal,
            pasos: Some(pasos_espacio_lineal),
        },
        Algoritmo {
            id: "alineamiento.tres",
            capitulo: CAPITULO,
            nombre: "Subsecuencia común de tres secuencias",
            resumen: "Alinear tres secuencias a la vez premiando las columnas donde las tres letras coinciden.",
            campos: &[
                Campo {
                    id: "u",
                    etiqueta: "Secuencia u",
                    tipo: Tipo::Proteina,
                    ejemplo: ">ThrA 1-45\nMRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITN",
                },
                Campo {
                    id: "v",
                    etiqueta: "Secuencia v",
                    tipo: Tipo::Proteina,
                    ejemplo: ">MetL 1-45\nMSVIAQAGAKGRQLHKFGGSSLADVKCYLRVAGIMAEYSQPDDMM",
                },
                Campo {
                    id: "w",
                    etiqueta: "Secuencia w",
                    tipo: Tipo::Proteina,
                    ejemplo: ">LysC 1-45\nMSEIVVSKFGGTSVADFDAMNRSADIVLSDANVRLVVLSASAGIT",
                },
            ],
            funcion: "alineamiento::subsecuencia_comun_de_tres",
            ejecutar: tres,
            pasos: Some(pasos_tres),
        },
    ]
}

fn cambio(entrada: &Entrada) -> Result<Salida, String> {
    let (dinero, monedas) = leer_cambio(entrada)?;
    let tabla = alineamiento::tabla_de_cambio(dinero, &monedas);
    let usadas = alineamiento::cambio_minimo(dinero, &monedas);
    let codicioso = cambio_codicioso(dinero, &monedas);
    let mut salida = Salida::default();
    match &usadas {
        Some(lista) => salida.agregar(Bloque::valor("Monedas necesarias", miles(lista.len()))),
        None => salida.agregar(Bloque::valor("Monedas necesarias", "imposible")),
    }
    salida.agregar(Bloque::valor(
        "Con la regla codiciosa (la moneda más grande que cabe)",
        match codicioso {
            Some(cuenta) => miles(cuenta),
            None => "se atasca".to_string(),
        },
    ));
    if let Some(lista) = &usadas
        && !lista.is_empty()
    {
        let mut categorias = Vec::new();
        let mut valores = Vec::new();
        for &moneda in &monedas {
            let veces = lista.iter().filter(|&&usada| usada == moneda).count();
            if veces > 0 {
                categorias.push(miles(moneda));
                valores.push(veces as f64);
            }
        }
        salida.agregar(Bloque::barras(
            "Monedas de la solución",
            "Valor de la moneda",
            "Cuántas",
            categorias,
            valores,
        ));
    }
    let puntos: Vec<[f64; 2]> = tabla
        .iter()
        .enumerate()
        .filter_map(|(cantidad, minimo)| minimo.map(|cuenta| [cantidad as f64, cuenta as f64]))
        .collect();
    salida.agregar(Bloque::linea(
        "Mínimo de monedas para cada cantidad",
        "Cantidad",
        "Monedas",
        &puntos,
        Vec::new(),
    ));
    let mut peores = 0;
    let mut primera_peor = None;
    let mut imposibles = 0;
    let mut primera_imposible = None;
    for (cantidad, minimo) in tabla.iter().enumerate().skip(1) {
        match minimo {
            None => {
                imposibles += 1;
                primera_imposible.get_or_insert(cantidad);
            }
            Some(cuenta) => {
                if cambio_codicioso(cantidad, &monedas) != Some(*cuenta) {
                    peores += 1;
                    primera_peor.get_or_insert(cantidad);
                }
            }
        }
    }
    match &usadas {
        Some(lista) if lista.is_empty() => {
            salida.explicar("Para pagar 0 no hace falta ninguna moneda.");
        }
        Some(lista) => {
            salida.explicar(format!(
                "Para pagar {}, el mínimo es {}: {}.",
                miles(dinero),
                cuantos(lista.len(), "moneda", "monedas"),
                desglose(lista, &monedas)
            ));
            match codicioso {
                Some(cuenta) if cuenta > lista.len() => salida.explicar(format!(
                    "Si eliges siempre la moneda más grande que cabe usarías {cuenta}. La programación dinámica gana porque prueba, para cada cantidad, todas las monedas posibles como última."
                )),
                Some(_) => salida.explicar(
                    "Aquí la regla codiciosa (la moneda más grande que cabe) también llega al mínimo.",
                ),
                None => salida.explicar(
                    "La regla codiciosa se atasca: tomando siempre la moneda más grande queda un resto que ninguna moneda paga, aunque sí había solución.",
                ),
            }
        }
        None => {
            salida.explicar(format!(
                "Ninguna combinación de estas monedas suma exactamente {}.",
                miles(dinero)
            ));
        }
    }
    if let Some(cantidad) = primera_peor {
        salida.explicar(format!(
            "Entre 1 y {}, la regla codiciosa falla en {} (la primera es {}): por eso hace falta una tabla y no basta con una regla simple.",
            miles(dinero),
            cuantos(peores, "cantidad", "cantidades"),
            miles(cantidad)
        ));
    } else if dinero > 0 && imposibles == 0 {
        salida.explicar(format!(
            "Para todas las cantidades de 1 a {} la regla codiciosa da el mínimo: con este juego de monedas (hasta esa cantidad) basta con ella.",
            miles(dinero)
        ));
    }
    if let Some(cantidad) = primera_imposible {
        salida.explicar(format!(
            "Entre 1 y {}, {} no se puede pagar con estas monedas (la primera es {}). Con una moneda de 1 todas serían posibles.",
            miles(dinero),
            cuantos(imposibles, "cantidad", "cantidades"),
            miles(cantidad)
        ));
    }
    salida.explicar(format!(
        "La tabla tiene {} casillas y en cada una se prueba{}: unos {} pasos, en vez de probar todas las combinaciones.",
        miles(dinero + 1),
        if monedas.len() == 1 {
            " una moneda".to_string()
        } else {
            format!("n {} monedas", monedas.len())
        },
        miles((dinero + 1) * monedas.len())
    ));
    Ok(salida)
}

fn cambio_codicioso(dinero: usize, monedas_de_mayor_a_menor: &[usize]) -> Option<usize> {
    let mut resto = dinero;
    let mut cuenta = 0;
    for &moneda in monedas_de_mayor_a_menor {
        cuenta += resto / moneda;
        resto %= moneda;
    }
    if resto == 0 { Some(cuenta) } else { None }
}

fn desglose(usadas: &[usize], monedas: &[usize]) -> String {
    let mut partes = Vec::new();
    for &moneda in monedas {
        let veces = usadas.iter().filter(|&&usada| usada == moneda).count();
        if veces > 0 {
            partes.push(format!("{veces} × {}", miles(moneda)));
        }
    }
    partes.join(" + ")
}

fn manhattan(entrada: &Entrada) -> Result<Salida, String> {
    let (abajo, derecha) = leer_manhattan(entrada)?;
    let recorrido = alineamiento::turista_de_manhattan(&abajo, &derecha)?;
    let filas = recorrido.mejor.len();
    let columnas = recorrido.mejor[0].len();
    let codicioso = turista_codicioso(&abajo, &derecha);
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        "Peso del mejor recorrido",
        numero(recorrido.longitud),
    ));
    salida.agregar(Bloque::valor("Turista codicioso", numero(codicioso)));
    let mut flechas = String::new();
    let mut bajadas = 0;
    for par in recorrido.camino.windows(2) {
        if par[1].0 > par[0].0 {
            flechas.push('↓');
            bajadas += 1;
        } else {
            flechas.push('→');
        }
    }
    salida.agregar(Bloque::texto("Recorrido desde (0, 0)", flechas));
    if filas <= LIMITE_GRAFO_VISIBLE && columnas <= LIMITE_GRAFO_VISIBLE {
        salida.agregar(Bloque::calor(
            "Mejor peso para llegar a cada esquina",
            (0..filas).map(|i| i.to_string()).collect(),
            (0..columnas).map(|j| j.to_string()).collect(),
            recorrido.mejor.clone(),
        ));
    }
    salida.agregar(Bloque::lista(
        "Esquinas del recorrido (fila, columna), contando desde 0",
        recorrido.camino.iter().map(|(i, j)| format!("({i}, {j})")),
    ));
    salida.explicar(format!(
        "El recorrido más valioso suma {}: baja {} veces y va a la derecha {} veces, como todo recorrido de esta cuadrícula.",
        numero(recorrido.longitud),
        bajadas,
        recorrido.camino.len() - 1 - bajadas
    ));
    if codicioso < recorrido.longitud {
        salida.explicar(format!(
            "Un turista codicioso, que en cada esquina toma la calle más valiosa, sumaría {}: pierde {} por no mirar más allá del paso siguiente.",
            numero(codicioso),
            numero(recorrido.longitud - codicioso)
        ));
    } else {
        salida.explicar(
            "Aquí el turista codicioso tiene suerte y llega al mismo valor; con otros pesos no pasa.",
        );
    }
    salida.explicar(
        "Cada esquina guarda el mejor valor para llegar a ella. Como solo se llega desde arriba o desde la izquierda, basta comparar dos números por esquina, fila por fila.",
    );
    Ok(salida)
}

fn revisar_pesos(abajo: &[Vec<f64>], derecha: &[Vec<f64>]) -> Result<(), String> {
    let demasiado_grande = abajo
        .iter()
        .chain(derecha)
        .flatten()
        .any(|peso| peso.abs() > LIMITE_PESO);
    if demasiado_grande {
        return Err(format!(
            "Los pesos deben estar entre -{} y {}: más grandes harían que las sumas se desborden.",
            miles(LIMITE_PESO as usize),
            miles(LIMITE_PESO as usize)
        ));
    }
    Ok(())
}

fn turista_codicioso(abajo: &[Vec<f64>], derecha: &[Vec<f64>]) -> f64 {
    let ultima_fila = derecha.len() - 1;
    let ultima_columna = derecha[0].len();
    let (mut i, mut j) = (0, 0);
    let mut suma = 0.0;
    while i < ultima_fila || j < ultima_columna {
        let puede_bajar = i < ultima_fila;
        let puede_girar = j < ultima_columna;
        if puede_bajar && (!puede_girar || abajo[i][j] >= derecha[i][j]) {
            suma += abajo[i][j];
            i += 1;
        } else {
            suma += derecha[i][j];
            j += 1;
        }
    }
    suma
}

fn numero(valor: f64) -> String {
    if valor.fract() == 0.0 && valor.abs() < 1e15 {
        let signo = if valor < 0.0 { "-" } else { "" };
        format!("{signo}{}", miles(valor.abs() as usize))
    } else {
        decimal(valor, 2)
    }
}

struct Grafo {
    nombres: Vec<String>,
    aristas: Vec<(usize, usize, f64)>,
}

fn leer_grafo(entrada: &Entrada, con_pesos: bool) -> Result<Grafo, String> {
    let texto = entrada.texto("aristas")?;
    let mut nombres: Vec<String> = Vec::new();
    let mut indices: HashMap<String, usize> = HashMap::new();
    let mut aristas = Vec::new();
    for (numero_de_linea, linea) in texto.lines().enumerate() {
        let linea = linea.trim();
        if linea.is_empty() {
            continue;
        }
        let (origen, resto) = linea.split_once("->").ok_or_else(|| {
            format!(
                "En la línea {}, falta la flecha «->»: «{linea}».",
                numero_de_linea + 1
            )
        })?;
        let (destinos, peso) = if con_pesos {
            let (destino, peso) = resto.rsplit_once(':').ok_or_else(|| {
                format!(
                    "En la línea {}, escribe el peso después de dos puntos, como «a -> b: 3».",
                    numero_de_linea + 1
                )
            })?;
            let peso: f64 = peso.trim().parse().map_err(|_| {
                format!(
                    "En la línea {}, «{}» no es un número.",
                    numero_de_linea + 1,
                    peso.trim()
                )
            })?;
            if !peso.is_finite() || peso.abs() > LIMITE_PESO {
                return Err(format!(
                    "En la línea {}, el peso debe estar entre -{} y {}.",
                    numero_de_linea + 1,
                    miles(LIMITE_PESO as usize),
                    miles(LIMITE_PESO as usize)
                ));
            }
            (vec![destino], peso)
        } else {
            (resto.split(',').collect(), 0.0)
        };
        let desde = registrar_nodo(origen, &mut nombres, &mut indices, numero_de_linea)?;
        for destino in destinos {
            let hasta = registrar_nodo(destino, &mut nombres, &mut indices, numero_de_linea)?;
            aristas.push((desde, hasta, peso));
        }
        if nombres.len() > LIMITE_NODOS || aristas.len() > LIMITE_ARISTAS {
            return Err(format!(
                "El grafo es demasiado grande: el límite es {} nodos y {} flechas.",
                miles(LIMITE_NODOS),
                miles(LIMITE_ARISTAS)
            ));
        }
    }
    if aristas.is_empty() {
        return Err("Escribe al menos una flecha, como «a -> b».".to_string());
    }
    Ok(Grafo { nombres, aristas })
}

fn registrar_nodo(
    texto: &str,
    nombres: &mut Vec<String>,
    indices: &mut HashMap<String, usize>,
    numero_de_linea: usize,
) -> Result<usize, String> {
    let nombre = texto.trim();
    if nombre.is_empty() {
        return Err(format!(
            "En la línea {}, falta el nombre de un nodo.",
            numero_de_linea + 1
        ));
    }
    if let Some(&indice) = indices.get(nombre) {
        return Ok(indice);
    }
    nombres.push(nombre.to_string());
    indices.insert(nombre.to_string(), nombres.len() - 1);
    Ok(nombres.len() - 1)
}

fn buscar_nodo(grafo: &Grafo, entrada: &Entrada, campo: &str) -> Result<usize, String> {
    let texto = entrada.texto(campo)?;
    let nombre = texto.trim();
    grafo
        .nombres
        .iter()
        .position(|candidato| candidato == nombre)
        .ok_or_else(|| format!("El nodo «{nombre}» no aparece en ninguna flecha."))
}

fn bloque_de_grafo(grafo: &Grafo, con_pesos: bool) -> Bloque {
    Bloque::grafo(
        "Grafo",
        grafo.nombres.clone(),
        grafo
            .aristas
            .iter()
            .map(|&(desde, hasta, peso)| Arista {
                desde,
                hasta,
                etiqueta: if con_pesos {
                    numero(peso)
                } else {
                    String::new()
                },
            })
            .collect(),
    )
}

fn orden_topologico(entrada: &Entrada) -> Result<Salida, String> {
    let grafo = leer_grafo(entrada, false)?;
    let pares: Vec<(usize, usize)> = grafo
        .aristas
        .iter()
        .map(|&(desde, hasta, _)| (desde, hasta))
        .collect();
    let orden = alineamiento::orden_topologico(grafo.nombres.len(), &pares)?;
    let conocidas: HashSet<(usize, usize)> = pares.iter().copied().collect();
    let unico = orden
        .windows(2)
        .all(|par| conocidas.contains(&(par[0], par[1])));
    let mut tiene_entrada = vec![false; grafo.nombres.len()];
    let mut tiene_salida = vec![false; grafo.nombres.len()];
    for &(desde, hasta) in &pares {
        tiene_salida[desde] = true;
        tiene_entrada[hasta] = true;
    }
    let comienzos: Vec<&str> = (0..grafo.nombres.len())
        .filter(|&nodo| !tiene_entrada[nodo])
        .map(|nodo| grafo.nombres[nodo].as_str())
        .collect();
    let finales: Vec<&str> = (0..grafo.nombres.len())
        .filter(|&nodo| !tiene_salida[nodo])
        .map(|nodo| grafo.nombres[nodo].as_str())
        .collect();
    let mut salida = Salida::default();
    salida.agregar(Bloque::lista(
        "Orden topológico",
        orden.iter().map(|&nodo| grafo.nombres[nodo].clone()),
    ));
    salida.agregar(Bloque::valor("Nodos", miles(grafo.nombres.len())));
    salida.agregar(Bloque::valor("Flechas", miles(pares.len())));
    if grafo.nombres.len() <= LIMITE_GRAFO_VISIBLE {
        salida.agregar(bloque_de_grafo(&grafo, false));
    }
    salida.explicar(format!(
        "En la lista, cada una de las {} flechas va de un nodo a otro que aparece más adelante: puedes hacer las tareas en ese orden sin que falte nada antes.",
        miles(pares.len())
    ));
    if unico {
        salida.explicar(
            "Este orden es el único posible: cada nodo de la lista tiene una flecha hacia el siguiente.",
        );
    } else {
        salida.explicar(
            "Hay otros órdenes válidos: en algún momento había más de un nodo listo y se eligió uno. Lo que no puede cambiar es el sentido de las flechas.",
        );
    }
    salida.explicar(format!(
        "Se puede empezar por {} (nada llega a ellos) y se termina en {} (no salen flechas).",
        enumerar(&comienzos),
        enumerar(&finales)
    ));
    salida.explicar(
        "Recorrer los nodos en este orden es lo que permite calcular el camino más largo, y en un alineamiento la tabla se llena fila por fila justamente porque ese es un orden topológico.",
    );
    Ok(salida)
}

fn enumerar(nombres: &[&str]) -> String {
    const MAXIMO: usize = 5;
    let mut partes: Vec<String> = nombres
        .iter()
        .take(MAXIMO)
        .map(|nombre| format!("«{nombre}»"))
        .collect();
    if nombres.len() > MAXIMO {
        partes.push(format!("y {} más", nombres.len() - MAXIMO));
    }
    partes.join(", ")
}

fn camino_dag(entrada: &Entrada) -> Result<Salida, String> {
    let grafo = leer_grafo(entrada, true)?;
    let fuente = buscar_nodo(&grafo, entrada, "fuente")?;
    let sumidero = buscar_nodo(&grafo, entrada, "sumidero")?;
    let camino =
        alineamiento::camino_mas_largo(grafo.nombres.len(), &grafo.aristas, fuente, sumidero)?;
    let nombres_del_camino: Vec<&str> = camino
        .nodos
        .iter()
        .map(|&nodo| grafo.nombres[nodo].as_str())
        .collect();
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Peso del camino", numero(camino.longitud)));
    salida.agregar(Bloque::texto("Camino", nombres_del_camino.join(" → ")));
    salida.agregar(Bloque::valor("Flechas usadas", camino.nodos.len() - 1));
    if grafo.nombres.len() <= LIMITE_GRAFO_VISIBLE {
        salida.agregar(bloque_de_grafo(&grafo, true));
    }
    salida.explicar(format!(
        "El camino más pesado de «{}» a «{}» suma {} y usa {} de las {} flechas del grafo.",
        grafo.nombres[fuente],
        grafo.nombres[sumidero],
        numero(camino.longitud),
        camino.nodos.len() - 1,
        grafo.aristas.len()
    ));
    if grafo.aristas.iter().any(|&(_, _, peso)| peso < 0.0) {
        salida.explicar(
            "Hay pesos negativos y el método funciona igual, porque sin ciclos ningún camino puede repetir nodos para sumar sin fin.",
        );
    }
    salida.explicar(
        "Cada nodo toma el mejor valor entre sus flechas de entrada, y por eso los nodos se visitan en orden topológico: cuando llega el turno de uno, todos sus anteriores ya están listos.",
    );
    salida.explicar(
        "Un alineamiento es este mismo problema en un grafo con forma de cuadrícula: cada camino de la esquina inicial a la final es un alineamiento posible.",
    );
    Ok(salida)
}

fn leer_secuencia(entrada: &Entrada, campo: &str) -> Result<(Vec<u8>, String), String> {
    let secuencia = entrada.proteina(campo)?;
    Ok((secuencia, nombre_de(entrada, campo)))
}

fn nombre_de(entrada: &Entrada, campo: &str) -> String {
    let texto = entrada.texto(campo).unwrap_or_default();
    let limpio = texto.trim();
    if let Some(nombre) = limpio
        .strip_prefix("@proteina:")
        .or_else(|| limpio.strip_prefix("@gen:"))
    {
        return nombre.trim().to_string();
    }
    if let Some(inicio) = inicio_en_ecoli(entrada, campo) {
        let largo = entrada.proteina(campo).map(|s| s.len()).unwrap_or(1);
        return format!(
            "E. coli {}–{}",
            miles(inicio),
            miles(inicio + largo.saturating_sub(1))
        );
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

fn inicio_en_ecoli(entrada: &Entrada, campo: &str) -> Option<usize> {
    let texto = entrada.texto(campo).ok()?;
    let rango = texto.trim().strip_prefix("@ecoli:")?;
    let (inicio, _) = rango.split_once('-')?;
    inicio.trim().parse().ok()
}

fn leer_puntuacion(entrada: &Entrada) -> Result<(Puntuacion, &'static str), String> {
    let elegida = entrada.opcion("puntuacion", PUNTUACIONES)?;
    let puntuar: Puntuacion = match elegida.as_str() {
        "BLOSUM62" => blosum62,
        "PAM250" => pam250,
        _ => alineamiento::puntaje_adn,
    };
    let nombre = PUNTUACIONES
        .iter()
        .find(|&&opcion| opcion == elegida)
        .copied()
        .unwrap_or("ADN +1/-1");
    Ok((puntuar, nombre))
}

fn leer_castigo(entrada: &Entrada, campo: &str) -> Result<i64, String> {
    let valor = entrada.entero(campo)?;
    if !(0..=1000).contains(&valor) {
        return Err(format!(
            "El castigo «{campo}» debe estar entre 0 y 1000: es lo que se resta del puntaje."
        ));
    }
    Ok(valor)
}

fn revisar_tamano(v: &[u8], w: &[u8], limite: usize, consejo: &str) -> Result<(), String> {
    let celdas = v.len().saturating_mul(w.len());
    if celdas > limite {
        return Err(format!(
            "La tabla tendría {} × {} = {} celdas y el límite de este método es {}. {consejo}",
            miles(v.len()),
            miles(w.len()),
            miles(celdas),
            miles(limite)
        ));
    }
    Ok(())
}

fn parece_adn(secuencia: &[u8]) -> bool {
    !secuencia.is_empty()
        && secuencia
            .iter()
            .all(|letra| matches!(letra, b'A' | b'C' | b'G' | b'T' | b'U' | b'N'))
}

struct Columnas {
    total: usize,
    pares: usize,
    identicas: usize,
    parecidas: usize,
    huecos: usize,
    tramos: usize,
    huecos_en_extremos: usize,
}

fn contar_columnas(fila_v: &[u8], fila_w: &[u8], puntuar: Puntuacion) -> Columnas {
    let mut cuenta = Columnas {
        total: fila_v.len().min(fila_w.len()),
        pares: 0,
        identicas: 0,
        parecidas: 0,
        huecos: 0,
        tramos: 0,
        huecos_en_extremos: 0,
    };
    let mut hueco_anterior_en_v = false;
    let mut hueco_anterior_en_w = false;
    for (&a, &b) in fila_v.iter().zip(fila_w) {
        let hueco_en_v = a == b'-';
        let hueco_en_w = b == b'-';
        if hueco_en_v || hueco_en_w {
            cuenta.huecos += 1;
            if (hueco_en_v && !hueco_anterior_en_v) || (hueco_en_w && !hueco_anterior_en_w) {
                cuenta.tramos += 1;
            }
        } else {
            cuenta.pares += 1;
            if a == b {
                cuenta.identicas += 1;
            }
            if puntuar(a, b) > 0 {
                cuenta.parecidas += 1;
            }
        }
        hueco_anterior_en_v = hueco_en_v;
        hueco_anterior_en_w = hueco_en_w;
    }
    let es_hueco = |(a, b): (&u8, &u8)| *a == b'-' || *b == b'-';
    let al_principio = fila_v
        .iter()
        .zip(fila_w)
        .take_while(|&par| es_hueco(par))
        .count();
    let al_final = fila_v
        .iter()
        .rev()
        .zip(fila_w.iter().rev())
        .take_while(|&par| es_hueco(par))
        .count();
    cuenta.huecos_en_extremos = if al_principio == cuenta.total {
        cuenta.total
    } else {
        al_principio + al_final
    };
    cuenta
}

fn fraccion(parte: usize, total: usize) -> f64 {
    if total == 0 {
        0.0
    } else {
        parte as f64 / total as f64
    }
}

fn entero(valor: i64) -> String {
    let signo = if valor < 0 { "-" } else { "" };
    format!("{signo}{}", miles(valor.unsigned_abs() as usize))
}

fn tramo(inicio: usize, largo: usize) -> String {
    if largo == 0 {
        "ninguno".to_string()
    } else {
        format!("{}–{}", miles(inicio + 1), miles(inicio + largo))
    }
}

fn letras(fila: &[u8]) -> usize {
    fila.iter().filter(|&&letra| letra != b'-').count()
}

fn cuantos(cantidad: usize, singular: &str, plural: &str) -> String {
    if cantidad == 1 {
        format!("1 {singular}")
    } else {
        format!("{} {plural}", miles(cantidad))
    }
}

fn nombre_con_tramo(
    entrada: &Entrada,
    campo: &str,
    nombre: &str,
    inicio: usize,
    largo: usize,
) -> String {
    match inicio_en_ecoli(entrada, campo) {
        Some(comienzo) if largo > 0 => format!(
            "E. coli {}–{}",
            miles(comienzo + inicio),
            miles(comienzo + inicio + largo - 1)
        ),
        _ => format!("{nombre} {}", tramo(inicio, largo)),
    }
}

fn agregar_alineamiento(
    salida: &mut Salida,
    alineado: &Alineamiento,
    nombre_v: &str,
    nombre_w: &str,
) {
    salida.agregar(Bloque::alineamiento(
        "Alineamiento",
        vec![
            fila_visible(nombre_v, &alineado.fila_v),
            fila_visible(nombre_w, &alineado.fila_w),
        ],
    ));
    avisar_si_se_recorta(salida, alineado.fila_v.len());
}

fn fila_visible(nombre: &str, fila: &[u8]) -> FilaAlineada {
    FilaAlineada {
        nombre: nombre.to_string(),
        secuencia: texto_visible(fila),
    }
}

fn texto_visible(letras: &[u8]) -> String {
    let visibles = &letras[..letras.len().min(MAXIMO_LETRAS_VISIBLES)];
    String::from_utf8_lossy(visibles).into_owned()
}

fn avisar_si_se_recorta(salida: &mut Salida, largo: usize) {
    if largo > MAXIMO_LETRAS_VISIBLES {
        salida.explicar(format!(
            "Se muestran las primeras {} de {} columnas; las cuentas usan el resultado completo.",
            miles(MAXIMO_LETRAS_VISIBLES),
            miles(largo)
        ));
    }
}

fn agregar_medidas(salida: &mut Salida, columnas: &Columnas, con_matriz: bool) {
    salida.agregar(Bloque::valor(
        "Identidad",
        format!(
            "{} ({} de {} columnas)",
            porcentaje(fraccion(columnas.identicas, columnas.total), 1),
            miles(columnas.identicas),
            miles(columnas.total)
        ),
    ));
    if con_matriz {
        salida.agregar(Bloque::valor(
            "Similitud (puntaje positivo)",
            format!(
                "{} ({} columnas)",
                porcentaje(fraccion(columnas.parecidas, columnas.total), 1),
                miles(columnas.parecidas)
            ),
        ));
    }
    salida.agregar(Bloque::valor(
        "Huecos",
        format!(
            "{} ({} en {})",
            porcentaje(fraccion(columnas.huecos, columnas.total), 1),
            cuantos(columnas.huecos, "columna", "columnas"),
            cuantos(columnas.tramos, "tramo", "tramos")
        ),
    ));
}

fn explicar_puntuacion(salida: &mut Salida, nombre_puntuacion: &str, v: &[u8], w: &[u8]) {
    let ambas_adn = parece_adn(v) && parece_adn(w);
    if nombre_puntuacion != "ADN +1/-1" && ambas_adn {
        salida.explicar(format!(
            "Las dos secuencias parecen ADN, pero {nombre_puntuacion} es una matriz para proteínas: sus puntajes no tienen sentido para bases. Elige «ADN +1/-1»."
        ));
    }
    if nombre_puntuacion == "ADN +1/-1" && !ambas_adn {
        salida.explicar(
            "Alguna secuencia no parece ADN: con «ADN +1/-1» todo cambio vale lo mismo, mientras que BLOSUM62 o PAM250 premian los cambios entre aminoácidos parecidos.",
        );
    }
}

fn explicar_parecido(salida: &mut Salida, columnas: &Columnas, adn: bool) {
    if columnas.total == 0 {
        salida.explicar(
            "No quedó ninguna columna alineada: con estos puntajes no hay ningún tramo que valga la pena emparejar.",
        );
        return;
    }
    let identidad = fraccion(columnas.identicas, columnas.total);
    let texto_identidad = porcentaje(identidad, 1);
    if adn {
        if identidad >= 0.9 {
            salida.explicar(format!(
                "Identidad de {texto_identidad}: son casi idénticas; difieren en {} de {} columnas.",
                miles(columnas.total - columnas.identicas),
                miles(columnas.total)
            ));
        } else if identidad >= 0.6 {
            salida.explicar(format!(
                "Identidad de {texto_identidad}: por encima de lo que da el azar. Dos ADN sin relación, alineados con huecos, llegan a coincidir en cerca de la mitad de las columnas."
            ));
        } else {
            salida.explicar(format!(
                "Identidad de {texto_identidad}: dentro de lo que da el azar, porque dos ADN sin relación, alineados con huecos, ya coinciden en cerca de la mitad de las columnas. Si son genes, compara mejor sus proteínas."
            ));
        }
    } else if identidad >= 0.4 && columnas.pares >= 70 {
        salida.explicar(format!(
            "Identidad de {texto_identidad} sobre {} residuos emparejados (40 % o más sobre al menos 70): casi seguro son homólogas, es decir, descienden de un ancestro común.",
            miles(columnas.pares)
        ));
    } else if identidad >= 0.4 {
        salida.explicar(format!(
            "Identidad alta ({texto_identidad}) pero solo en {} residuos emparejados: en tramos de menos de 70 residuos un 40 % puede salir por azar. Mira si el parecido sigue en un tramo más largo.",
            miles(columnas.pares)
        ));
    } else if identidad >= 0.35 {
        salida.explicar(format!(
            "Identidad de {texto_identidad}: justo por encima de la zona de penumbra (20–35 %). La homología es probable, sobre todo si el alineamiento es largo."
        ));
    } else if identidad >= 0.2 {
        salida.explicar(format!(
            "Identidad de {texto_identidad}: zona de penumbra (20–35 %), donde el parecido puede ser herencia o azar. Necesitas más evidencia: motivos conservados, la estructura, o comparar el puntaje con el de secuencias barajadas."
        ));
    } else {
        salida.explicar(format!(
            "Identidad de {texto_identidad}, menos de 20 %: este alineamiento no da evidencia de que estén emparentadas."
        ));
    }
    if !adn && columnas.parecidas > columnas.identicas {
        salida.explicar(format!(
            "Además, {} letras distintas pero con puntaje positivo (cambios conservadores, como I por V): contando también esas, la similitud es {}.",
            cuantos(
                columnas.parecidas - columnas.identicas,
                "columna tiene",
                "columnas tienen"
            ),
            porcentaje(fraccion(columnas.parecidas, columnas.total), 1)
        ));
    }
    if columnas.huecos > 0 {
        let mut texto = format!(
            "Los huecos ocupan {} de las columnas, en {}.",
            porcentaje(fraccion(columnas.huecos, columnas.total), 1),
            cuantos(columnas.tramos, "tramo", "tramos")
        );
        let unidad = if adn { "bases" } else { "aminoácidos" };
        if columnas.huecos_en_extremos * 2 >= columnas.huecos {
            if !adn && columnas.huecos_en_extremos >= 40 {
                texto.push_str(&format!(
                    " La mayoría está en los extremos: una secuencia tiene {} {unidad} de más en sus puntas, lo que puede ser un dominio entero (suelen medir de 40 a 300 aminoácidos).",
                    miles(columnas.huecos_en_extremos)
                ));
            } else {
                texto.push_str(
                    " Están sobre todo en los extremos: una secuencia empieza o termina un poco antes que la otra.",
                );
            }
        }
        salida.explicar(texto);
    }
}

fn global(entrada: &Entrada) -> Result<Salida, String> {
    let (v, nombre_v) = leer_secuencia(entrada, "v")?;
    let (w, nombre_w) = leer_secuencia(entrada, "w")?;
    let (puntuar, nombre_puntuacion) = leer_puntuacion(entrada)?;
    let indel = leer_castigo(entrada, "indel")?;
    revisar_tamano(&v, &w, LIMITE_CUADRATICO, CONSEJO_CUADRATICO)?;
    let alineado = alineamiento::alineamiento_global(&v, &w, puntuar, indel);
    Ok(resultado_de_extremo_a_extremo(
        &alineado,
        (&v, &w),
        (&nombre_v, &nombre_w),
        puntuar,
        nombre_puntuacion,
    ))
}

fn resultado_de_extremo_a_extremo(
    alineado: &Alineamiento,
    secuencias: (&[u8], &[u8]),
    nombres: (&str, &str),
    puntuar: Puntuacion,
    nombre_puntuacion: &str,
) -> Salida {
    let (v, w) = secuencias;
    let columnas = contar_columnas(&alineado.fila_v, &alineado.fila_w, puntuar);
    let con_matriz = nombre_puntuacion != "ADN +1/-1";
    let mut salida = Salida::default();
    agregar_alineamiento(&mut salida, alineado, nombres.0, nombres.1);
    salida.agregar(Bloque::valor("Puntaje", entero(alineado.puntaje)));
    agregar_medidas(&mut salida, &columnas, con_matriz);
    salida.agregar(Bloque::valor(
        "Longitudes",
        format!("{} y {} letras", miles(v.len()), miles(w.len())),
    ));
    explicar_puntuacion(&mut salida, nombre_puntuacion, v, w);
    explicar_parecido(&mut salida, &columnas, !con_matriz && parece_adn(v));
    let corta = v.len().min(w.len());
    let larga = v.len().max(w.len());
    if corta * 10 < larga * 7 {
        salida.explicar(
            "Una secuencia es mucho más larga que la otra (más de 30 % de diferencia): un alineamiento global las obliga a emparejarse de punta a punta; si solo comparten una parte, prueba el alineamiento local.",
        );
    }
    salida
}

fn local(entrada: &Entrada) -> Result<Salida, String> {
    let (v, nombre_v) = leer_secuencia(entrada, "v")?;
    let (w, nombre_w) = leer_secuencia(entrada, "w")?;
    let (puntuar, nombre_puntuacion) = leer_puntuacion(entrada)?;
    let indel = leer_castigo(entrada, "indel")?;
    revisar_tamano(&v, &w, LIMITE_CUADRATICO, CONSEJO_CUADRATICO)?;
    let alineado = alineamiento::alineamiento_local(&v, &w, puntuar, indel);
    let columnas = contar_columnas(&alineado.fila_v, &alineado.fila_w, puntuar);
    let con_matriz = nombre_puntuacion != "ADN +1/-1";
    let largo_v = letras(&alineado.fila_v);
    let largo_w = letras(&alineado.fila_w);
    let mut salida = Salida::default();
    agregar_alineamiento(
        &mut salida,
        &alineado,
        &nombre_con_tramo(entrada, "v", &nombre_v, alineado.inicio_v, largo_v),
        &nombre_con_tramo(entrada, "w", &nombre_w, alineado.inicio_w, largo_w),
    );
    salida.agregar(Bloque::valor("Puntaje", entero(alineado.puntaje)));
    salida.agregar(Bloque::valor(
        "Tramo de v (contando desde 1)",
        tramo(alineado.inicio_v, largo_v),
    ));
    salida.agregar(Bloque::valor(
        "Tramo de w (contando desde 1)",
        tramo(alineado.inicio_w, largo_w),
    ));
    agregar_medidas(&mut salida, &columnas, con_matriz);
    explicar_puntuacion(&mut salida, nombre_puntuacion, &v, &w);
    if columnas.total == 0 {
        salida.explicar(
            "No hay ningún tramo con puntaje positivo: con esta puntuación, ninguna parte de una se parece a la otra.",
        );
        return Ok(salida);
    }
    let cubre_v = fraccion(largo_v, v.len());
    let cubre_w = fraccion(largo_w, w.len());
    salida.explicar(format!(
        "El tramo parecido cubre {} de v (letras {}) y {} de w (letras {}), contando desde 1.",
        porcentaje(cubre_v, 0),
        tramo(alineado.inicio_v, largo_v),
        porcentaje(cubre_w, 0),
        tramo(alineado.inicio_w, largo_w)
    ));
    if cubre_v.min(cubre_w) < 0.5 {
        salida.explicar(
            "Menos de la mitad de alguna de las dos entra en el tramo: comparten una parte (por ejemplo un dominio), no la secuencia entera.",
        );
    } else if cubre_v.min(cubre_w) >= 0.9 {
        salida.explicar(
            "El tramo abarca casi todo de las dos: el parecido recorre las secuencias completas y un alineamiento global daría algo muy similar.",
        );
    }
    explicar_parecido(&mut salida, &columnas, !con_matriz && parece_adn(&v));
    Ok(salida)
}

fn edicion(entrada: &Entrada) -> Result<Salida, String> {
    let (v, nombre_v) = leer_secuencia(entrada, "v")?;
    let (w, nombre_w) = leer_secuencia(entrada, "w")?;
    revisar_tamano(&v, &w, LIMITE_LINEAL, CONSEJO_LINEAL)?;
    let distancia = alineamiento::distancia_de_edicion(&v, &w);
    let larga = v.len().max(w.len());
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Distancia de edición", miles(distancia)));
    salida.agregar(Bloque::valor(
        "Cambios por cada 100 letras",
        decimal(100.0 * fraccion(distancia, larga), 2),
    ));
    salida.explicar(format!(
        "La distancia es {}: el menor número de cambios (sustituir, insertar o borrar una letra) que convierten v en w.",
        miles(distancia)
    ));
    if distancia == 0 {
        salida.explicar("Las dos secuencias son idénticas.");
        return Ok(salida);
    }
    salida.explicar(format!(
        "Son {} cambios por cada 100 letras de la secuencia más larga.",
        decimal(100.0 * fraccion(distancia, larga), 2)
    ));
    if v.len().saturating_mul(w.len()) > LIMITE_CUADRATICO {
        salida.explicar(
            "Las secuencias son largas: se calculó solo la distancia, con memoria lineal, sin mostrar el alineamiento.",
        );
        return Ok(salida);
    }
    let alineado = alineamiento::alineamiento_global(&v, &w, alineamiento::puntaje_edicion, 1);
    agregar_alineamiento(&mut salida, &alineado, &nombre_v, &nombre_w);
    let mut sustituciones = 0;
    let mut inserciones = 0;
    let mut borrados = 0;
    for (&a, &b) in alineado.fila_v.iter().zip(&alineado.fila_w) {
        if a == b'-' {
            inserciones += 1;
        } else if b == b'-' {
            borrados += 1;
        } else if a != b {
            sustituciones += 1;
        }
    }
    salida.agregar(Bloque::barras(
        "Tipos de cambio",
        "Cambio",
        "Cuántos",
        vec![
            "Sustituciones".to_string(),
            "Inserciones".to_string(),
            "Borrados".to_string(),
        ],
        vec![sustituciones as f64, inserciones as f64, borrados as f64],
    ));
    salida.explicar(format!(
        "Un camino de cambios mínimo: {}, {} y {} (puede haber otros con el mismo total).",
        cuantos(sustituciones, "sustitución", "sustituciones"),
        cuantos(inserciones, "letra insertada", "letras insertadas"),
        cuantos(borrados, "letra borrada", "letras borradas")
    ));
    let solo_sustituciones = inserciones == 0 && borrados == 0 && sustituciones > 0;
    if parece_adn(&v) && parece_adn(&w) && solo_sustituciones && v.len() % 3 == 0 {
        explicar_codones(&mut salida, &v, &w, sustituciones);
    }
    let tramos = tramos_de_huecos(&alineado.fila_v, &alineado.fila_w);
    if parece_adn(&v) && parece_adn(&w) && !tramos.is_empty() {
        if tramos.iter().all(|largo| largo % 3 == 0) {
            salida.explicar(
                "Todas las inserciones y borrados vienen en tramos de 3, 6, 9… bases: si son genes, el marco de lectura se conserva y la proteína solo gana o pierde aminoácidos.",
            );
        } else {
            salida.explicar(
                "Algún tramo insertado o borrado no mide un múltiplo de 3: si son genes, eso correría el marco de lectura y cambiaría toda la proteína desde ahí.",
            );
        }
    }
    Ok(salida)
}

fn explicar_codones(salida: &mut Salida, v: &[u8], w: &[u8], sustituciones: usize) {
    let mut por_posicion = [0usize; 3];
    for (indice, (a, b)) in v.iter().zip(w).enumerate() {
        if a != b {
            por_posicion[indice % 3] += 1;
        }
    }
    let proteina_v = secuencia::traducir(v);
    let proteina_w = secuencia::traducir(w);
    let es_arranque = |adn: &[u8]| matches!(adn.get(..3), Some(b"ATG" | b"GTG" | b"TTG"));
    let primer_codon = if es_arranque(v) && es_arranque(w) {
        1
    } else {
        0
    };
    let aminoacidos_distintos = proteina_v
        .iter()
        .zip(&proteina_w)
        .skip(primer_codon)
        .filter(|(a, b)| a != b)
        .count();
    salida.agregar(Bloque::barras(
        "Sustituciones según la posición en el codón",
        "Posición en el codón",
        "Sustituciones",
        vec!["1.ª".to_string(), "2.ª".to_string(), "3.ª".to_string()],
        por_posicion.iter().map(|&cuenta| cuenta as f64).collect(),
    ));
    salida.explicar(format!(
        "Solo hay sustituciones y el largo es múltiplo de 3. Si son genes leídos desde su primer codón, {} de las {} caen en la tercera base del codón y las proteínas que codifican difieren en {} (un codón de inicio ATG, GTG o TTG siempre se lee como metionina).",
        miles(por_posicion[2]),
        miles(sustituciones),
        cuantos(aminoacidos_distintos, "aminoácido", "aminoácidos")
    ));
    if aminoacidos_distintos * 2 < sustituciones {
        salida.explicar(
            "La mayoría de los cambios son silenciosos: el código genético es redundante, sobre todo en la tercera base, y por eso dos genes divergen más que las proteínas que codifican.",
        );
    }
}

fn tramos_de_huecos(fila_v: &[u8], fila_w: &[u8]) -> Vec<usize> {
    let mut tramos = Vec::new();
    for fila in [fila_v, fila_w] {
        let mut largo = 0;
        for &letra in fila {
            if letra == b'-' {
                largo += 1;
            } else if largo > 0 {
                tramos.push(largo);
                largo = 0;
            }
        }
        if largo > 0 {
            tramos.push(largo);
        }
    }
    tramos
}

fn ajustado(entrada: &Entrada) -> Result<Salida, String> {
    let (v, nombre_v) = leer_secuencia(entrada, "v")?;
    let (w, nombre_w) = leer_secuencia(entrada, "w")?;
    let (puntuar, nombre_puntuacion) = leer_puntuacion(entrada)?;
    let indel = leer_castigo(entrada, "indel")?;
    revisar_tamano(&v, &w, LIMITE_CUADRATICO, CONSEJO_CUADRATICO)?;
    let alineado = alineamiento::alineamiento_ajustado(&v, &w, puntuar, indel);
    let columnas = contar_columnas(&alineado.fila_v, &alineado.fila_w, puntuar);
    let con_matriz = nombre_puntuacion != "ADN +1/-1";
    let largo_v = letras(&alineado.fila_v);
    let mut salida = Salida::default();
    agregar_alineamiento(
        &mut salida,
        &alineado,
        &nombre_con_tramo(entrada, "v", &nombre_v, alineado.inicio_v, largo_v),
        &nombre_w,
    );
    salida.agregar(Bloque::valor("Puntaje", entero(alineado.puntaje)));
    salida.agregar(Bloque::valor(
        "Dónde encaja w en v (contando desde 1)",
        tramo(alineado.inicio_v, largo_v),
    ));
    agregar_medidas(&mut salida, &columnas, con_matriz);
    explicar_puntuacion(&mut salida, nombre_puntuacion, &v, &w);
    if w.len() > v.len() {
        salida.explicar(
            "w es más larga que v: en un alineamiento ajustado la corta debe caber dentro de la larga. Prueba intercambiarlas.",
        );
    }
    if largo_v > 0 {
        salida.explicar(format!(
            "w encaja en las letras {} de v (contando desde 1).",
            tramo(alineado.inicio_v, largo_v)
        ));
        if let Some(inicio) = inicio_en_ecoli(entrada, "v") {
            salida.explicar(format!(
                "Como v empieza en la posición {} del genoma de E. coli, eso corresponde a las posiciones {}–{} del genoma (la letra 1 de v es la {}).",
                miles(inicio),
                miles(inicio + alineado.inicio_v),
                miles(inicio + alineado.inicio_v + largo_v - 1),
                miles(inicio)
            ));
        }
    }
    if columnas.identicas == w.len() && columnas.total == w.len() && !w.is_empty() {
        salida.explicar("Es una copia exacta de w, sin cambios ni huecos.");
    } else {
        salida.explicar(format!(
            "No es una copia exacta: de {} columnas, {} con letras distintas o con hueco.",
            miles(columnas.total),
            cuantos(columnas.total - columnas.identicas, "tiene", "tienen")
        ));
    }
    if parece_adn(&v) && parece_adn(&w) && !w.is_empty() {
        let esperadas = (v.len().saturating_sub(w.len()) + 1) as f64 * 0.25f64.powi(w.len() as i32);
        salida.explicar(format!(
            "En ADN al azar con las cuatro bases igual de probables, una copia exacta de {} bases aparecería unas {} veces en una secuencia del largo de v.",
            w.len(),
            cientifico(esperadas)
        ));
    }
    Ok(salida)
}

fn solapado(entrada: &Entrada) -> Result<Salida, String> {
    let (v, nombre_v) = leer_secuencia(entrada, "v")?;
    let (w, nombre_w) = leer_secuencia(entrada, "w")?;
    let (puntuar, nombre_puntuacion) = leer_puntuacion(entrada)?;
    let indel = leer_castigo(entrada, "indel")?;
    revisar_tamano(&v, &w, LIMITE_CUADRATICO, CONSEJO_CUADRATICO)?;
    let alineado = alineamiento::alineamiento_solapado(&v, &w, puntuar, indel);
    let columnas = contar_columnas(&alineado.fila_v, &alineado.fila_w, puntuar);
    let con_matriz = nombre_puntuacion != "ADN +1/-1";
    let largo_v = letras(&alineado.fila_v);
    let largo_w = letras(&alineado.fila_w);
    let mut salida = Salida::default();
    agregar_alineamiento(
        &mut salida,
        &alineado,
        &nombre_con_tramo(entrada, "v", &nombre_v, alineado.inicio_v, largo_v),
        &nombre_con_tramo(entrada, "w", &nombre_w, 0, largo_w),
    );
    salida.agregar(Bloque::valor("Puntaje", entero(alineado.puntaje)));
    salida.agregar(Bloque::valor(
        "Largo del solapamiento",
        miles(columnas.total),
    ));
    agregar_medidas(&mut salida, &columnas, con_matriz);
    explicar_puntuacion(&mut salida, nombre_puntuacion, &v, &w);
    if columnas.total == 0 {
        salida.explicar(
            "No hay solapamiento con puntaje positivo: el final de v no se parece al comienzo de w.",
        );
        return Ok(salida);
    }
    let mut unida = v[..alineado.inicio_v + largo_v].to_vec();
    unida.extend_from_slice(&w[largo_w..]);
    salida.agregar(Bloque::texto(
        "Secuencia unida (v y luego lo que sobra de w)",
        texto_visible(&unida),
    ));
    if unida.len() > MAXIMO_LETRAS_VISIBLES {
        salida.explicar(format!(
            "Se muestran las primeras {} de las {} letras de la secuencia unida.",
            miles(MAXIMO_LETRAS_VISIBLES),
            miles(unida.len())
        ));
    }
    let diferencias = columnas.total - columnas.identicas;
    salida.explicar(format!(
        "Las últimas {} letras de v se superponen con las primeras {} de w, con {}. Unidas forman una secuencia de {} letras.",
        miles(largo_v),
        miles(largo_w),
        cuantos(diferencias, "diferencia", "diferencias"),
        miles(unida.len())
    ));
    if parece_adn(&v) && parece_adn(&w) {
        let probabilidad = probabilidad_de_azar(columnas.total, diferencias);
        salida.explicar(format!(
            "Si fuera azar (cuatro bases igual de probables), un tramo de {} bases coincidiría con a lo sumo {} con probabilidad {}.",
            miles(columnas.total),
            cuantos(diferencias, "diferencia", "diferencias"),
            cientifico(probabilidad)
        ));
        if probabilidad * GENOMA_ECOLI < 0.01 {
            salida.explicar(
                "Aun probando contra cada una de las 4.6 millones de posiciones del genoma de E. coli, esperarías menos de 0.01 coincidencias así por azar: lo más probable es que las dos lecturas vengan del mismo lugar y se puedan unir.",
            );
        } else {
            salida.explicar(
                "Es un solapamiento corto o con muchas diferencias: en un genoma del tamaño del de E. coli algo así aparece por azar. Para unir lecturas con confianza hacen falta solapamientos más largos.",
            );
        }
    }
    explicar_parecido(&mut salida, &columnas, !con_matriz && parece_adn(&v));
    Ok(salida)
}

fn probabilidad_de_azar(largo: usize, diferencias: usize) -> f64 {
    let mut total = 0.0;
    let mut logaritmo_combinaciones = 0.0;
    for k in 0..=diferencias.min(largo) {
        if k > 0 {
            logaritmo_combinaciones += ((largo - k + 1) as f64).ln() - (k as f64).ln();
        }
        let logaritmo = logaritmo_combinaciones + k as f64 * 3f64.ln() - largo as f64 * 4f64.ln();
        total += logaritmo.exp();
    }
    total.min(1.0)
}

fn afin(entrada: &Entrada) -> Result<Salida, String> {
    let (v, nombre_v) = leer_secuencia(entrada, "v")?;
    let (w, nombre_w) = leer_secuencia(entrada, "w")?;
    let (puntuar, nombre_puntuacion) = leer_puntuacion(entrada)?;
    let apertura = leer_castigo(entrada, "apertura")?;
    let extension = leer_castigo(entrada, "extension")?;
    if extension > apertura {
        return Err(format!(
            "La extensión ({extension}) no puede ser mayor que la apertura ({apertura}): si alargar un hueco costara más que abrir otro, ya no serían huecos afines."
        ));
    }
    revisar_tamano(&v, &w, LIMITE_CUADRATICO, CONSEJO_CUADRATICO)?;
    let alineado = alineamiento::alineamiento_afin(&v, &w, puntuar, apertura, extension);
    let mut salida = resultado_de_extremo_a_extremo(
        &alineado,
        (&v, &w),
        (&nombre_v, &nombre_w),
        puntuar,
        nombre_puntuacion,
    );
    let tramos = tramos_de_huecos(&alineado.fila_v, &alineado.fila_w);
    if !tramos.is_empty() {
        let total: usize = tramos.iter().sum();
        let mas_largo = tramos.iter().copied().max().unwrap_or(0);
        salida.explicar(format!(
            "De media cada tramo de huecos mide {} letras y el más largo {}. Cada tramo cuesta {apertura} por abrirse y {extension} por cada letra más.",
            decimal(total as f64 / tramos.len() as f64, 1),
            miles(mas_largo)
        ));
        if extension < apertura && total > tramos.len() {
            salida.explicar(
                "Como alargar un hueco cuesta menos que abrir otro, el método prefiere pocos huecos largos a muchos cortos, que es lo que pasa en la evolución: una inserción o un borrado suele afectar a varias letras seguidas.",
            );
        }
    }
    if (apertura, extension) != (11, 1) && nombre_puntuacion == "BLOSUM62" {
        salida.explicar(
            "Lo habitual con BLOSUM62 es apertura 11 y extensión 1; con otros valores el método cambia cuántos huecos prefiere y cuán largos.",
        );
    }
    Ok(salida)
}

fn arista_central(entrada: &Entrada) -> Result<Salida, String> {
    let (v, _) = leer_secuencia(entrada, "v")?;
    let (w, _) = leer_secuencia(entrada, "w")?;
    let (puntuar, nombre_puntuacion) = leer_puntuacion(entrada)?;
    let indel = leer_castigo(entrada, "indel")?;
    revisar_tamano(&v, &w, LIMITE_LINEAL, CONSEJO_LINEAL)?;
    let arista = alineamiento::arista_central(&v, &w, puntuar, indel)
        .ok_or_else(|| "Hacen falta letras en alguna de las dos secuencias.".to_string())?;
    let (i, j) = arista.desde;
    let (siguiente_i, siguiente_j) = arista.hasta;
    let mejor = arista.longitudes[i];
    let empates = arista
        .longitudes
        .iter()
        .filter(|&&longitud| longitud == mejor)
        .count();
    let tipo = if siguiente_i > i && siguiente_j > j {
        format!(
            "diagonal: empareja la letra {} de v («{}») con la letra {} de w («{}»)",
            i + 1,
            v[i] as char,
            j + 1,
            w[j] as char
        )
    } else if siguiente_j > j {
        format!(
            "horizontal: la letra {} de w («{}») va frente a un hueco",
            j + 1,
            w[j] as char
        )
    } else {
        format!(
            "vertical: la letra {} de v («{}») va frente a un hueco",
            i + 1,
            v[i] as char
        )
    };
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        "Arista central (nodos contando desde 0)",
        format!("({i}, {j}) → ({siguiente_i}, {siguiente_j})"),
    ));
    salida.agregar(Bloque::valor(
        "Columna central",
        format!("{j} de 0 a {}", w.len()),
    ));
    salida.agregar(Bloque::valor(
        "Puntaje del alineamiento global",
        entero(mejor),
    ));
    let puntos: Vec<[f64; 2]> = arista
        .longitudes
        .iter()
        .enumerate()
        .map(|(fila, &longitud)| [fila as f64, longitud as f64])
        .collect();
    salida.agregar(Bloque::linea(
        "Mejor camino que pasa por cada nodo de la columna central",
        "Fila del nodo (contando desde 0)",
        "Puntaje del camino",
        &puntos,
        vec![Marca {
            x: i as f64,
            texto: format!("nodo central ({i}, {j})"),
        }],
    ));
    explicar_puntuacion(&mut salida, nombre_puntuacion, &v, &w);
    salida.explicar(format!(
        "La columna central es la {j} (la mitad de las {} letras de w). De sus {} nodos, el de la fila {i} es por donde pasa el mejor camino: su puntaje, {}, es el del alineamiento global completo.",
        miles(w.len()),
        miles(v.len() + 1),
        entero(mejor)
    ));
    salida.explicar(format!("La arista que sale de ahí es {tipo}."));
    if empates > 1 {
        salida.explicar(format!(
            "{} nodos de la columna empatan con ese puntaje: hay más de un alineamiento óptimo.",
            miles(empates)
        ));
    }
    salida.explicar(format!(
        "Para saberlo bastó guardar dos columnas de {} números, no la tabla de {} celdas: este es el paso que repite el alineamiento en espacio lineal en cada mitad.",
        miles(v.len() + 1),
        miles(v.len().saturating_mul(w.len()))
    ));
    Ok(salida)
}

fn espacio_lineal(entrada: &Entrada) -> Result<Salida, String> {
    let (v, nombre_v) = leer_secuencia(entrada, "v")?;
    let (w, nombre_w) = leer_secuencia(entrada, "w")?;
    let (puntuar, nombre_puntuacion) = leer_puntuacion(entrada)?;
    let indel = leer_castigo(entrada, "indel")?;
    revisar_tamano(&v, &w, LIMITE_LINEAL, CONSEJO_LINEAL)?;
    let alineado = alineamiento::alineamiento_espacio_lineal(&v, &w, puntuar, indel);
    let mut salida = resultado_de_extremo_a_extremo(
        &alineado,
        (&v, &w),
        (&nombre_v, &nombre_w),
        puntuar,
        nombre_puntuacion,
    );
    let tabla = (v.len() + 1).saturating_mul(w.len() + 1);
    let columnas_guardadas = 4 * (v.len() + 1) * 8;
    salida.agregar(Bloque::valor("Celdas de una tabla completa", miles(tabla)));
    salida.explicar(format!(
        "Una tabla completa de flechas tendría {} celdas (unos {} MB, a un byte por celda). Este método solo guarda unas pocas columnas de {} números (unos {} kB) y llega al mismo puntaje, a cambio de hacer cerca del doble de cuentas.",
        miles(tabla),
        decimal(tabla as f64 / 1_048_576.0, 1),
        miles(v.len() + 1),
        miles(columnas_guardadas.div_ceil(1024))
    ));
    Ok(salida)
}

fn tres(entrada: &Entrada) -> Result<Salida, String> {
    let (u, nombre_u) = leer_secuencia(entrada, "u")?;
    let (v, nombre_v) = leer_secuencia(entrada, "v")?;
    let (w, nombre_w) = leer_secuencia(entrada, "w")?;
    let celdas = (u.len() + 1)
        .saturating_mul(v.len() + 1)
        .saturating_mul(w.len() + 1);
    if celdas > LIMITE_TRES {
        return Err(format!(
            "El cubo tendría {} celdas y el límite es {}: con tres secuencias el trabajo crece con el producto de los tres largos. Usa tramos de unas 200 letras.",
            miles(celdas),
            miles(LIMITE_TRES)
        ));
    }
    let triple = alineamiento::subsecuencia_comun_de_tres(&u, &v, &w);
    let columnas = triple.fila_u.len();
    let mut mejor_inicio = 0;
    let mut mejor_largo = 0;
    let mut largo = 0;
    for k in 0..columnas {
        let iguales = triple.fila_u[k] != b'-'
            && triple.fila_u[k] == triple.fila_v[k]
            && triple.fila_v[k] == triple.fila_w[k];
        if iguales {
            largo += 1;
            if largo > mejor_largo {
                mejor_largo = largo;
                mejor_inicio = k + 1 - largo;
            }
        } else {
            largo = 0;
        }
    }
    let mut salida = Salida::default();
    salida.agregar(Bloque::alineamiento(
        "Alineamiento de las tres",
        vec![
            fila_visible(&nombre_u, &triple.fila_u),
            fila_visible(&nombre_v, &triple.fila_v),
            fila_visible(&nombre_w, &triple.fila_w),
        ],
    ));
    avisar_si_se_recorta(&mut salida, columnas);
    salida.agregar(Bloque::valor(
        "Columnas con las tres letras iguales",
        miles(triple.puntaje),
    ));
    salida.agregar(Bloque::valor("Columnas del alineamiento", miles(columnas)));
    let corta = u.len().min(v.len()).min(w.len());
    salida.explicar(format!(
        "Las tres comparten {} letras en el mismo orden ({} de la más corta): cada una es una columna donde coinciden las tres.",
        miles(triple.puntaje),
        porcentaje(fraccion(triple.puntaje, corta), 0)
    ));
    let minimo_seguido = if parece_adn(&u) && parece_adn(&v) && parece_adn(&w) {
        8
    } else {
        3
    };
    if mejor_largo >= minimo_seguido {
        let bloque: String =
            String::from_utf8_lossy(&triple.fila_u[mejor_inicio..mejor_inicio + mejor_largo])
                .into_owned();
        salida.explicar(format!(
            "El tramo seguido más largo donde las tres coinciden es «{bloque}» ({mejor_largo} columnas). Un bloque conservado así en tres parientes suele marcar un sitio importante para la función."
        ));
    } else if triple.puntaje > 0 {
        salida.explicar(
            "Las coincidencias están sueltas, sin tramos seguidos lo bastante largos para no ser casualidad: puede ser solo azar.",
        );
    }
    if [&triple.fila_u, &triple.fila_v, &triple.fila_w]
        .iter()
        .any(|fila| fila.contains(&b'-'))
    {
        salida.explicar(
            "Aquí los huecos no cuestan nada, por eso aparecen siempre que ayudan a sumar una columna más.",
        );
    }
    salida.explicar(format!(
        "El cubo tiene {} celdas: con cada secuencia que se agrega el trabajo se multiplica, y por eso para muchas secuencias se usan métodos aproximados.",
        miles(celdas)
    ));
    Ok(salida)
}

fn subsecuencia_comun(entrada: &Entrada) -> Result<Salida, String> {
    let (v, _) = leer_secuencia(entrada, "v")?;
    let (w, _) = leer_secuencia(entrada, "w")?;
    let semilla = entrada.natural("semilla")? as u64;
    revisar_tamano(&v, &w, LIMITE_CUADRATICO, CONSEJO_LINEAL)?;
    let comun = alineamiento::subsecuencia_comun_mas_larga(&v, &w);
    let mut azar = Azar::nuevo(semilla);
    let mut barajadas = Vec::with_capacity(BARAJADOS);
    for _ in 0..BARAJADOS {
        let v_barajada = barajar(&v, &mut azar);
        let w_barajada = barajar(&w, &mut azar);
        barajadas.push(alineamiento::subsecuencia_comun_mas_larga(&v_barajada, &w_barajada).len());
    }
    let promedio = barajadas.iter().sum::<usize>() as f64 / BARAJADOS as f64;
    let maximo = barajadas.iter().copied().max().unwrap_or(0);
    let corta = v.len().min(w.len());
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Longitud", miles(comun.len())));
    salida.agregar(Bloque::valor(
        "Parte de la secuencia más corta",
        porcentaje(fraccion(comun.len(), corta), 1),
    ));
    salida.agregar(Bloque::barras(
        "Comparación con secuencias barajadas",
        "Pareja",
        "Longitud de la subsecuencia común",
        std::iter::once("Real".to_string())
            .chain((1..=BARAJADOS).map(|vuelta| format!("Barajada {vuelta}")))
            .collect(),
        std::iter::once(comun.len() as f64)
            .chain(barajadas.iter().map(|&largo| largo as f64))
            .collect(),
    ));
    salida.agregar(Bloque::texto("Subsecuencia común", texto_visible(&comun)));
    salida.explicar(format!(
        "La subsecuencia común más larga tiene {} letras: {} de la secuencia más corta.",
        miles(comun.len()),
        porcentaje(fraccion(comun.len(), corta), 1)
    ));
    if comun.len() == corta && corta > 0 {
        salida.explicar("La más corta entra entera en la otra, letra por letra y en orden.");
    }
    salida.explicar(format!(
        "Al barajar las letras de las dos ({BARAJADOS} veces, semilla {semilla}), que conserva cuántas hay de cada una pero rompe el orden, sale en promedio {} y como mucho {}.",
        decimal(promedio, 1),
        miles(maximo)
    ));
    if comun.len() > maximo && promedio > 0.0 {
        let exceso = comun.len() as f64 / promedio - 1.0;
        if exceso >= 0.2 {
            salida.explicar(format!(
                "La real es {} más larga que el promedio barajado: las secuencias comparten orden, no solo composición, lo que apunta a un origen común.",
                porcentaje(exceso, 0)
            ));
        } else {
            salida.explicar(format!(
                "La real supera a todas las barajadas, pero solo por {}: la señal es débil. Un alineamiento con matriz de puntajes y castigo por huecos separa mejor el parecido del azar.",
                porcentaje(exceso, 0)
            ));
        }
    } else {
        salida.explicar(
            "La real no supera a las barajadas: una subsecuencia común así de larga sale por azar con letras en cualquier orden. Con un alfabeto de 4 letras, al azar ya se comparte cerca de dos tercios.",
        );
    }
    Ok(salida)
}

fn barajar(letras: &[u8], azar: &mut Azar) -> Vec<u8> {
    let mut copia = letras.to_vec();
    for i in (1..copia.len()).rev() {
        let j = azar.entero(i + 1);
        copia.swap(i, j);
    }
    copia
}

struct Tablero {
    v: Vec<u8>,
    w: Vec<u8>,
    nombre_v: String,
    nombre_w: String,
    puntuacion: &'static str,
    indel: i64,
    matriz: Vec<Vec<Option<i64>>>,
    flechas: Vec<Vec<Option<&'static str>>>,
    camino: Vec<(usize, usize)>,
    fila_v: Vec<u8>,
    fila_w: Vec<u8>,
    puntaje: Option<i64>,
}

impl Tablero {
    fn estado(&self, fase: &str, celda: Option<(usize, usize)>, detalle: Value) -> Value {
        let al_derecho =
            |fila: &[u8]| -> String { fila.iter().rev().map(|&letra| letra as char).collect() };
        json!({
            "fase": fase,
            "v": String::from_utf8_lossy(&self.v),
            "w": String::from_utf8_lossy(&self.w),
            "nombre_v": self.nombre_v,
            "nombre_w": self.nombre_w,
            "puntuacion": self.puntuacion,
            "indel": self.indel,
            "matriz": self.matriz,
            "flechas": self.flechas,
            "celda": celda.map(|(i, j)| vec![i, j]),
            "candidatos": detalle.get("candidatos").cloned().unwrap_or(Value::Null),
            "par": detalle.get("par").cloned().unwrap_or(Value::Null),
            "ganador": detalle.get("ganador").cloned().unwrap_or(Value::Null),
            "camino": self.camino.iter().map(|&(i, j)| vec![i, j]).collect::<Vec<_>>(),
            "fila_v": al_derecho(&self.fila_v),
            "fila_w": al_derecho(&self.fila_w),
            "puntaje": self.puntaje,
        })
    }
}

fn pasos_global(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (v, nombre_v) = leer_secuencia(entrada, "v")?;
    let (w, nombre_w) = leer_secuencia(entrada, "w")?;
    if v.len() > LIMITE_PASO_A_PASO || w.len() > LIMITE_PASO_A_PASO {
        return Err(format!(
            "El paso a paso dibuja la tabla entera, así que solo acepta secuencias de hasta {LIMITE_PASO_A_PASO} letras (estas tienen {} y {}). Usa un tramo corto.",
            miles(v.len()),
            miles(w.len())
        ));
    }
    let (puntuar, puntuacion) = leer_puntuacion(entrada)?;
    let indel = leer_castigo(entrada, "indel")?;
    let (n, m) = (v.len(), w.len());
    let mut tablero = Tablero {
        v: v.clone(),
        w: w.clone(),
        nombre_v,
        nombre_w,
        puntuacion,
        indel,
        matriz: vec![vec![None; m + 1]; n + 1],
        flechas: vec![vec![None; m + 1]; n + 1],
        camino: Vec::new(),
        fila_v: Vec::new(),
        fila_w: Vec::new(),
        puntaje: None,
    };
    let mut pasos = Vec::new();
    tablero.matriz[0][0] = Some(0);
    pasos.push(Paso::nuevo(
        "let mut flechas = Flechas::nuevas(v.len() + 1, w.len() + 1);",
        format!(
            "Una tabla de {} × {}: la fila i y la columna j guardan el mejor puntaje para alinear las primeras i letras de v con las primeras j de w. La esquina (0, 0) vale 0.",
            n + 1,
            m + 1
        ),
        tablero.estado("inicio", Some((0, 0)), Value::Null),
    ));
    for j in 1..=m {
        tablero.matriz[0][j] = Some(-(j as i64) * indel);
        tablero.flechas[0][j] = Some("izquierda");
    }
    pasos.push(Paso::nuevo(
        "anterior[j] = anterior[j - 1] - indel;",
        format!(
            "Fila 0: alinear el comienzo de w solo contra huecos cuesta {indel} por cada letra."
        ),
        tablero.estado("primera_fila", None, Value::Null),
    ));
    for i in 1..=n {
        tablero.matriz[i][0] = Some(-(i as i64) * indel);
        tablero.flechas[i][0] = Some("arriba");
        pasos.push(Paso::nuevo(
            "actual[0] = anterior[0] - indel;",
            format!(
                "Fila {i}, columna 0: las primeras {i} letras de v contra huecos valen {}.",
                -(i as i64) * indel
            ),
            tablero.estado("primera_columna", Some((i, 0)), Value::Null),
        ));
        for j in 1..=m {
            let par = puntuar(v[i - 1], w[j - 1]) as i64;
            let diagonal = tablero.matriz[i - 1][j - 1].unwrap_or(0) + par;
            let arriba = tablero.matriz[i - 1][j].unwrap_or(0) - indel;
            let izquierda = tablero.matriz[i][j - 1].unwrap_or(0) - indel;
            let (puntaje, ganador) = if diagonal >= arriba && diagonal >= izquierda {
                (diagonal, "diagonal")
            } else if arriba >= izquierda {
                (arriba, "arriba")
            } else {
                (izquierda, "izquierda")
            };
            tablero.matriz[i][j] = Some(puntaje);
            tablero.flechas[i][j] = Some(ganador);
            let explicacion = match ganador {
                "diagonal" => format!(
                    "emparejar «{}» con «{}» ({par})",
                    v[i - 1] as char,
                    w[j - 1] as char
                ),
                "arriba" => format!("poner «{}» frente a un hueco", v[i - 1] as char),
                _ => format!("poner un hueco frente a «{}»", w[j - 1] as char),
            };
            pasos.push(Paso::nuevo(
                "let (puntaje, flecha) = elegir_mejor(diagonal, arriba, izquierda);",
                format!(
                    "Celda ({i}, {j}): diagonal {diagonal}, arriba {arriba}, izquierda {izquierda}. Gana {ganador} con {puntaje}: {explicacion}."
                ),
                tablero.estado(
                    "celda",
                    Some((i, j)),
                    json!({
                        "candidatos": {"diagonal": diagonal, "arriba": arriba, "izquierda": izquierda},
                        "par": {"letra_v": (v[i - 1] as char).to_string(), "letra_w": (w[j - 1] as char).to_string(), "puntaje": par},
                        "ganador": ganador,
                    }),
                ),
            ));
        }
    }
    let final_puntaje = tablero.matriz[n][m].unwrap_or(0);
    tablero.puntaje = Some(final_puntaje);
    tablero.camino.push((n, m));
    pasos.push(Paso::nuevo(
        "Modo::Global => (ultima_fila[w.len()], v.len(), w.len()),",
        format!(
            "La tabla está llena. El puntaje del mejor alineamiento global está en la última esquina ({n}, {m}): {final_puntaje}. Ahora se sigue el camino de flechas hacia atrás."
        ),
        tablero.estado("puntaje_final", Some((n, m)), Value::Null),
    ));
    let (mut i, mut j) = (n, m);
    while i > 0 || j > 0 {
        let flecha = tablero.flechas[i][j].unwrap_or("diagonal");
        let (linea, nota) = match flecha {
            "diagonal" => {
                tablero.fila_v.push(v[i - 1]);
                tablero.fila_w.push(w[j - 1]);
                let nota = format!(
                    "Desde ({i}, {j}) la flecha es diagonal: «{}» va con «{}».",
                    v[i - 1] as char,
                    w[j - 1] as char
                );
                i -= 1;
                j -= 1;
                ("DIAGONAL => {", nota)
            }
            "arriba" => {
                tablero.fila_v.push(v[i - 1]);
                tablero.fila_w.push(b'-');
                let nota = format!(
                    "Desde ({i}, {j}) la flecha viene de arriba: «{}» va frente a un hueco.",
                    v[i - 1] as char
                );
                i -= 1;
                ("ARRIBA => {", nota)
            }
            _ => {
                tablero.fila_v.push(b'-');
                tablero.fila_w.push(w[j - 1]);
                let nota = format!(
                    "Desde ({i}, {j}) la flecha viene de la izquierda: un hueco va frente a «{}».",
                    w[j - 1] as char
                );
                j -= 1;
                ("IZQUIERDA => {", nota)
            }
        };
        tablero.camino.push((i, j));
        pasos.push(Paso::nuevo(
            linea,
            nota,
            tablero.estado("camino", Some((i, j)), Value::Null),
        ));
    }
    let al_derecho_v: String = tablero.fila_v.iter().rev().map(|&l| l as char).collect();
    let al_derecho_w: String = tablero.fila_w.iter().rev().map(|&l| l as char).collect();
    pasos.push(Paso::nuevo(
        "fila_v.reverse();",
        format!(
            "Las letras se juntaron de atrás hacia adelante; al darlas vuelta queda el alineamiento {al_derecho_v} / {al_derecho_w}, con puntaje {final_puntaje}."
        ),
        tablero.estado("fin", Some((0, 0)), Value::Null),
    ));
    Ok(reducir_pasos(pasos))
}

fn reducir_pasos(pasos: Vec<Paso>) -> Vec<Paso> {
    let total = pasos.len();
    if total <= MAXIMO_PASOS {
        return pasos;
    }
    let mut elegidos = Vec::with_capacity(MAXIMO_PASOS);
    let mut siguiente = 0;
    for (indice, paso) in pasos.into_iter().enumerate() {
        if indice == siguiente * (total - 1) / (MAXIMO_PASOS - 1) {
            elegidos.push(paso);
            siguiente += 1;
        }
    }
    elegidos
}

const LIMITE_NARRADO: usize = 2_000_000;
const FILAS_VISIBLES: usize = 8;
const COLUMNAS_VISIBLES: usize = 8;

struct Ritmo {
    detalle: usize,
    cada: usize,
}

impl Ritmo {
    fn nuevo(total: usize, detalle: usize, objetivo: usize) -> Ritmo {
        Ritmo {
            detalle,
            cada: (total.saturating_sub(detalle) / objetivo.max(1)).max(1),
        }
    }

    fn toca(&self, indice: usize) -> bool {
        indice < self.detalle || (indice - self.detalle).is_multiple_of(self.cada)
    }
}

fn corto(texto: &[u8]) -> String {
    let visible: String = texto.iter().take(24).map(|&letra| letra as char).collect();
    if texto.len() > 24 {
        format!("{visible}…")
    } else {
        visible
    }
}

fn letra(secuencia: &[u8], posicion: usize) -> char {
    secuencia.get(posicion).map_or('·', |&valor| valor as char)
}

fn revisar_narrado(v: &[u8], w: &[u8]) -> Result<(), String> {
    revisar_tamano(
        v,
        w,
        LIMITE_NARRADO,
        "El paso a paso guarda la tabla entera para enseñártela, así que usa un tramo más corto.",
    )
}

struct Cuadricula<'a> {
    titulo: &'a str,
    tamano: (usize, usize),
    rotulo_fila: &'a dyn Fn(usize) -> String,
    rotulo_columna: &'a dyn Fn(usize) -> String,
    valor: &'a dyn Fn(usize, usize) -> String,
}

impl Cuadricula<'_> {
    fn dibujar(
        &self,
        pizarra: Pizarra,
        centro: (usize, usize),
        resaltadas: &[(usize, usize)],
    ) -> Pizarra {
        let (filas, columnas) = self.tamano;
        let desde_fila = centro
            .0
            .saturating_sub(FILAS_VISIBLES * 2 / 3)
            .min(filas.saturating_sub(FILAS_VISIBLES));
        let hasta_fila = (desde_fila + FILAS_VISIBLES).min(filas);
        let desde_columna = centro
            .1
            .saturating_sub(COLUMNAS_VISIBLES * 2 / 3)
            .min(columnas.saturating_sub(COLUMNAS_VISIBLES));
        let hasta_columna = (desde_columna + COLUMNAS_VISIBLES).min(columnas);
        let mut encabezados = vec![String::new()];
        encabezados.extend((desde_columna..hasta_columna).map(|j| (self.rotulo_columna)(j)));
        let mut cuerpo = Vec::new();
        for i in desde_fila..hasta_fila {
            let mut fila = vec![(self.rotulo_fila)(i)];
            fila.extend((desde_columna..hasta_columna).map(|j| (self.valor)(i, j)));
            cuerpo.push(fila);
        }
        let mut marcadas = Vec::new();
        for &(i, j) in std::iter::once(&centro).chain(resaltadas) {
            if (desde_fila..hasta_fila).contains(&i) && (desde_columna..hasta_columna).contains(&j)
            {
                marcadas.push((i - desde_fila, j - desde_columna + 1));
            }
        }
        let referencias: Vec<&str> = encabezados.iter().map(String::as_str).collect();
        pizarra.tabla(self.titulo, &referencias, cuerpo, &marcadas)
    }
}

fn rotulo_de_letra(secuencia: &[u8]) -> impl Fn(usize) -> String + '_ {
    move |indice| {
        if indice == 0 {
            "0".to_string()
        } else {
            format!("{indice} {}", letra(secuencia, indice - 1))
        }
    }
}

fn rotulo_numerado(indice: usize) -> String {
    indice.to_string()
}

fn texto_de_celda(valor: i64) -> String {
    if valor <= IMPOSIBLE_VISIBLE {
        "−∞".to_string()
    } else {
        valor.to_string()
    }
}

const IMPOSIBLE_VISIBLE: i64 = i64::MIN / 8;

fn leer_cambio(entrada: &Entrada) -> Result<(usize, Vec<usize>), String> {
    let dinero = entrada.entero("dinero")?;
    if !(0..=LIMITE_DINERO).contains(&dinero) {
        return Err(format!(
            "La cantidad debe estar entre 0 y {}.",
            miles(LIMITE_DINERO as usize)
        ));
    }
    let dinero = dinero as usize;
    let mut monedas = Vec::new();
    for valor in entrada.enteros("monedas")? {
        if !(1..=LIMITE_DINERO).contains(&valor) {
            return Err(format!(
                "Cada moneda debe valer entre 1 y {}; «{valor}» no sirve.",
                miles(LIMITE_DINERO as usize)
            ));
        }
        monedas.push(valor as usize);
    }
    monedas.sort_unstable_by(|a, b| b.cmp(a));
    monedas.dedup();
    if monedas.len() > LIMITE_MONEDAS {
        return Err(format!(
            "Escribe como mucho {LIMITE_MONEDAS} valores de moneda distintos."
        ));
    }
    Ok((dinero, monedas))
}

fn tabla_de_monedas(pizarra: Pizarra, menos_monedas: &[Option<usize>], actual: usize) -> Pizarra {
    let desde = actual.saturating_sub(10);
    let filas: Vec<Vec<String>> = (desde..=actual)
        .map(|cantidad| {
            vec![
                cantidad.to_string(),
                menos_monedas[cantidad].map_or("None".to_string(), |valor| valor.to_string()),
            ]
        })
        .collect();
    let ultima = filas.len() - 1;
    pizarra.tabla(
        "menos_monedas",
        &["cantidad", "menos monedas"],
        filas,
        &[(ultima, 1)],
    )
}

fn pasos_cambio(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (dinero, monedas) = leer_cambio(entrada)?;
    Ok(narrar_cambio(dinero, &monedas).0)
}

fn narrar_cambio(dinero: usize, monedas: &[usize]) -> (Vec<Paso>, Option<Vec<usize>>) {
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let lista_monedas = monedas
        .iter()
        .map(usize::to_string)
        .collect::<Vec<_>>()
        .join(" ");
    let base = || {
        Pizarra::nueva()
            .variable("dinero", dinero)
            .variable("monedas", format!("[{lista_monedas}]"))
    };
    let mut menos_monedas: Vec<Option<usize>> = vec![None; dinero + 1];
    narracion.agregar(base().paso(
        "let mut menos_monedas = vec![None; dinero + 1];",
        format!(
            "Quieres pagar {} con monedas de {lista_monedas}. Se prepara una casilla para cada cantidad de 0 a {}, donde irá el mínimo de monedas que la paga. Ahora todas están vacías (None).",
            miles(dinero),
            miles(dinero)
        ),
    ));
    menos_monedas[0] = Some(0);
    narracion.agregar(
        tabla_de_monedas(base(), &menos_monedas, 0).paso(
            "menos_monedas[0] = Some(0);",
            "Pagar 0 no necesita ninguna moneda: es el punto de partida. Todas las demás cantidades se calculan a partir de las anteriores.",
        ),
    );
    let ritmo = Ritmo::nuevo(dinero, 10, 40);
    let mut ya_detallada = false;
    for cantidad in 1..=dinero {
        let utiles = if ya_detallada {
            0
        } else {
            monedas
                .iter()
                .filter(|&&moneda| moneda <= cantidad && menos_monedas[cantidad - moneda].is_some())
                .count()
        };
        let detallar = (!ya_detallada && utiles >= 2) || cantidad == dinero;
        ya_detallada |= detallar;
        let mostrar = detallar || ritmo.toca(cantidad - 1);
        let mut pruebas: Vec<(usize, usize, usize)> = Vec::new();
        let mut pasos_de_moneda = 0;
        for &moneda in monedas {
            if moneda == 0 || moneda > cantidad {
                if detallar && pasos_de_moneda < 3 && moneda > cantidad {
                    pasos_de_moneda += 1;
                    narracion.agregar(
                        tabla_de_monedas(
                            base()
                                .variable("cantidad", cantidad)
                                .variable("moneda", moneda),
                            &menos_monedas,
                            cantidad,
                        )
                        .paso(
                            "if moneda == 0 || moneda > cantidad {",
                            format!(
                                "Para pagar {cantidad}, la moneda de {moneda} es demasiado grande: no cabe y se salta."
                            ),
                        ),
                    );
                }
                continue;
            }
            let Some(antes) = menos_monedas[cantidad - moneda] else {
                if detallar && pasos_de_moneda < 8 {
                    pasos_de_moneda += 1;
                    narracion.agregar(
                        tabla_de_monedas(
                            base()
                                .variable("cantidad", cantidad)
                                .variable("moneda", moneda),
                            &menos_monedas,
                            cantidad,
                        )
                        .paso(
                            "if let Some(antes) = menos_monedas[cantidad - moneda] {",
                            format!(
                                "Si la última moneda fuera de {moneda}, quedarían {} por pagar, y esa casilla está vacía: no se puede. Esta moneda no sirve aquí.",
                                cantidad - moneda
                            ),
                        ),
                    );
                }
                continue;
            };
            let con_esta_moneda = antes + 1;
            pruebas.push((moneda, antes, con_esta_moneda));
            let mejora = menos_monedas[cantidad].is_none_or(|actual| con_esta_moneda < actual);
            let anterior_mejor = menos_monedas[cantidad];
            if mejora {
                menos_monedas[cantidad] = Some(con_esta_moneda);
            }
            if detallar && pasos_de_moneda < 8 {
                pasos_de_moneda += 1;
                let nota = if mejora {
                    format!(
                        "Si la última moneda es de {moneda}, quedan {} por pagar, que necesitan {antes}; con esta moneda son {con_esta_moneda}. {} Se anota en la casilla.",
                        cantidad - moneda,
                        match anterior_mejor {
                            None => "Es la primera forma que encuentro.".to_string(),
                            Some(previo) => format!("Mejora los {previo} de antes."),
                        }
                    )
                } else {
                    format!(
                        "Si la última moneda es de {moneda}, quedan {} por pagar, que necesitan {antes}; con esta moneda son {con_esta_moneda}. No mejora los {} que ya había, así que no se cambia nada.",
                        cantidad - moneda,
                        anterior_mejor.unwrap_or(0)
                    )
                };
                narracion.agregar(
                    tabla_de_monedas(
                        base()
                            .variable("cantidad", cantidad)
                            .variable("moneda", moneda)
                            .variable("antes", antes)
                            .variable("con_esta_moneda", con_esta_moneda),
                        &menos_monedas,
                        cantidad,
                    )
                    .paso(
                        if mejora {
                            "menos_monedas[cantidad] = Some(con_esta_moneda);"
                        } else {
                            "if menos_monedas[cantidad].is_none_or(|actual| con_esta_moneda < actual) {"
                        },
                        nota,
                    ),
                );
            }
        }
        if mostrar && !detallar {
            let pizarra = base().variable("cantidad", cantidad);
            match menos_monedas[cantidad] {
                Some(mejor) => {
                    let (moneda, antes, _) = pruebas
                        .iter()
                        .copied()
                        .find(|&(_, _, total)| total == mejor)
                        .unwrap_or((0, 0, mejor));
                    let intentos = pruebas
                        .iter()
                        .take(4)
                        .map(|&(moneda, _, total)| format!("{moneda} → {total}"))
                        .collect::<Vec<_>>()
                        .join(", ");
                    let resto = if pruebas.len() > 4 { ", …" } else { "" };
                    narracion.agregar(
                        tabla_de_monedas(
                            pizarra
                                .variable("moneda", moneda)
                                .variable("antes", antes)
                                .variable("con_esta_moneda", mejor),
                            &menos_monedas,
                            cantidad,
                        )
                        .paso(
                            "menos_monedas[cantidad] = Some(con_esta_moneda);",
                            format!(
                                "Para pagar {cantidad} pruebo cada moneda que cabe ({intentos}{resto}). La mejor deja {} por pagar, que necesitan {antes}: en total {mejor}.",
                                cantidad - moneda
                            ),
                        ),
                    );
                }
                None => narracion.agregar(
                    tabla_de_monedas(pizarra, &menos_monedas, cantidad).paso(
                        "for cantidad in 1..=dinero {",
                        format!(
                            "Para pagar {cantidad} ninguna moneda sirve: todas dejan un resto imposible. La casilla se queda vacía (None)."
                        ),
                    ),
                ),
            }
        }
    }
    narracion.agregar(
        tabla_de_monedas(base(), &menos_monedas, dinero).paso(
            "let menos_monedas = tabla_de_cambio(dinero, monedas);",
            match menos_monedas[dinero] {
                Some(total) => format!(
                    "La tabla está completa: pagar {} necesita {total} monedas como mínimo. Falta saber cuáles son; se averigua caminando hacia atrás.",
                    miles(dinero)
                ),
                None => format!(
                    "La tabla está completa, pero la casilla de {} quedó vacía: ninguna combinación suma exactamente esa cantidad.",
                    miles(dinero)
                ),
            },
        ),
    );
    let mut usadas = Vec::new();
    let mut resto = dinero;
    let mut elecciones = Vec::new();
    while resto > 0 {
        let Some(necesarias) = menos_monedas[resto] else {
            break;
        };
        let elegida = monedas.iter().copied().find(|&moneda| {
            moneda > 0 && moneda <= resto && menos_monedas[resto - moneda] == Some(necesarias - 1)
        });
        let Some(elegida) = elegida else {
            break;
        };
        elecciones.push((resto, necesarias, elegida));
        usadas.push(elegida);
        resto -= elegida;
    }
    let completo = resto == 0;
    narracion.agregar(base().variable("resto", dinero).paso(
        "let mut resto = dinero;",
        if completo {
            format!(
                "Ahora se reconstruye la respuesta: se parte de resto = {} y se busca, una por una, una moneda que acerque el resto a 0 gastando una moneda menos.",
                miles(dinero)
            )
        } else {
            "No hay respuesta que reconstruir: la casilla de la cantidad está vacía.".to_string()
        },
    ));
    let ritmo_de_vuelta = Ritmo::nuevo(elecciones.len(), 6, 12);
    let mut acumuladas: Vec<String> = Vec::new();
    for (indice, &(antes_del_resto, necesarias, elegida)) in elecciones.iter().enumerate() {
        acumuladas.push(elegida.to_string());
        if ritmo_de_vuelta.toca(indice) {
            let visibles = acumuladas.len().saturating_sub(12);
            narracion.agregar(
                base()
                    .variable("resto", antes_del_resto - elegida)
                    .variable("necesarias", necesarias)
                    .variable("elegida", elegida)
                    .lista(
                        "usadas",
                        acumuladas[visibles..].to_vec(),
                        &[acumuladas.len() - visibles - 1],
                    )
                    .paso(
                        "usadas.push(elegida);",
                        format!(
                            "Quedaban {} por pagar y hacen falta {necesarias} monedas. Busco una moneda tras la cual el resto necesite una menos: la de {elegida} (queda {}). La agrego a la respuesta.",
                            miles(antes_del_resto),
                            miles(antes_del_resto - elegida)
                        ),
                    ),
            );
        }
    }
    if !completo && dinero > 0 {
        narracion.agregar(base().variable("resto", resto).paso(
            "let necesarias = menos_monedas[resto]?;",
            "La casilla está vacía, así que el signo ? termina la función y devuelve None: no hay forma de pagar esa cantidad.",
        ));
    }
    let resultado = if completo { Some(usadas.clone()) } else { None };
    let ultimo = match &resultado {
        Some(lista) if lista.is_empty() => base().paso(
            "Some(usadas)",
            "Pagar 0 no necesita monedas: la respuesta es la lista vacía. Mira «Qué significa» para la lectura completa.",
        ),
        Some(lista) => base()
            .variable("resto", 0)
            .lista(
                "usadas",
                lista.iter().take(24).map(usize::to_string).collect(),
                &[],
            )
            .paso(
                "Some(usadas)",
                format!(
                    "El resto llegó a 0. La respuesta son {}: {}. Mira «Qué significa» para ver cómo se compara con elegir siempre la moneda más grande.",
                    cuantos(lista.len(), "moneda", "monedas"),
                    desglose(lista, monedas)
                ),
            ),
        None => base().paso(
            "Some(usadas)",
            "No hubo respuesta: la función devolvió None. Mira «Qué significa».",
        ),
    };
    (narracion.terminar(ultimo), resultado)
}

type Calles = Vec<Vec<f64>>;

fn leer_manhattan(entrada: &Entrada) -> Result<(Calles, Calles), String> {
    let abajo = entrada.matriz("abajo")?;
    let derecha = entrada.matriz("derecha")?;
    revisar_pesos(&abajo, &derecha)?;
    let esquinas = (abajo.len() + 1).saturating_mul(derecha[0].len() + 1);
    if esquinas > LIMITE_ESQUINAS {
        return Err(format!(
            "La cuadrícula tendría {} esquinas; el límite es {}.",
            miles(esquinas),
            miles(LIMITE_ESQUINAS)
        ));
    }
    Ok((abajo, derecha))
}

fn pasos_manhattan(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (abajo, derecha) = leer_manhattan(entrada)?;
    alineamiento::turista_de_manhattan(&abajo, &derecha)?;
    Ok(narrar_manhattan(&abajo, &derecha).0)
}

fn flechas_del_camino(camino: &[(usize, usize)]) -> String {
    camino
        .windows(2)
        .map(|par| if par[1].0 > par[0].0 { '↓' } else { '→' })
        .take(40)
        .collect()
}

fn narrar_manhattan(
    abajo: &[Vec<f64>],
    derecha: &[Vec<f64>],
) -> (Vec<Paso>, f64, Vec<(usize, usize)>) {
    let filas = derecha.len();
    let columnas = derecha[0].len() + 1;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let mut mejor = vec![vec![0.0; columnas]; filas];
    let dibujar = |pizarra: Pizarra,
                   mejor: &[Vec<f64>],
                   centro: (usize, usize),
                   resaltadas: &[(usize, usize)]| {
        let valor = |i: usize, j: usize| numero(mejor[i][j]);
        Cuadricula {
            titulo: "mejor (peso máximo para llegar a cada esquina)",
            tamano: (filas, columnas),
            rotulo_fila: &rotulo_numerado,
            rotulo_columna: &rotulo_numerado,
            valor: &valor,
        }
        .dibujar(pizarra, centro, resaltadas)
    };
    narracion.agregar(
        dibujar(
            Pizarra::nueva()
                .variable("filas", filas)
                .variable("columnas", columnas),
            &mejor,
            (0, 0),
            &[],
        )
        .paso(
            "let mut mejor = vec![vec![0.0; columnas]; filas];",
            format!(
                "La ciudad tiene {filas} × {columnas} esquinas. El turista sale de la esquina (0, 0) y solo puede bajar o ir a la derecha. mejor[i][j] guardará el mayor peso con que se llega a cada esquina: la de salida vale 0."
            ),
        ),
    );
    let ritmo_borde = Ritmo::nuevo(filas - 1, 3, 3);
    for i in 1..filas {
        mejor[i][0] = mejor[i - 1][0] + abajo[i - 1][0];
        if ritmo_borde.toca(i - 1) || i == filas - 1 {
            narracion.agregar(
                dibujar(
                    Pizarra::nueva()
                        .variable("filas", filas)
                        .variable("columnas", columnas)
                        .variable("i", i),
                    &mejor,
                    (i, 0),
                    &[],
                )
                .paso(
                    "mejor[i][0] = mejor[i - 1][0] + abajo[i - 1][0];",
                    format!(
                        "En la primera columna solo se llega bajando: {} + la calle de peso {} = {}.",
                        numero(mejor[i - 1][0]),
                        numero(abajo[i - 1][0]),
                        numero(mejor[i][0])
                    ),
                ),
            );
        }
    }
    let ritmo_borde = Ritmo::nuevo(columnas - 1, 3, 3);
    for j in 1..columnas {
        mejor[0][j] = mejor[0][j - 1] + derecha[0][j - 1];
        if ritmo_borde.toca(j - 1) || j == columnas - 1 {
            narracion.agregar(
                dibujar(
                    Pizarra::nueva()
                        .variable("filas", filas)
                        .variable("columnas", columnas)
                        .variable("j", j),
                    &mejor,
                    (0, j),
                    &[],
                )
                .paso(
                    "mejor[0][j] = mejor[0][j - 1] + derecha[0][j - 1];",
                    format!(
                        "En la primera fila solo se llega yendo a la derecha: {} + la calle de peso {} = {}.",
                        numero(mejor[0][j - 1]),
                        numero(derecha[0][j - 1]),
                        numero(mejor[0][j])
                    ),
                ),
            );
        }
    }
    let total = (filas - 1) * (columnas - 1);
    let ritmo = Ritmo::nuevo(total, 5, 20);
    for i in 1..filas {
        for j in 1..columnas {
            let bajando = mejor[i - 1][j] + abajo[i - 1][j];
            let girando = mejor[i][j - 1] + derecha[i][j - 1];
            mejor[i][j] = bajando.max(girando);
            let indice = (i - 1) * (columnas - 1) + (j - 1);
            if ritmo.toca(indice) || indice == total - 1 {
                let veredicto = if bajando > girando {
                    "Llegar bajando vale más: gana la calle de arriba."
                } else if girando > bajando {
                    "Llegar girando vale más: gana la calle de la izquierda."
                } else {
                    "Empatan: da igual por dónde llegues."
                };
                narracion.agregar(
                    dibujar(
                        Pizarra::nueva()
                            .variable("filas", filas)
                            .variable("columnas", columnas)
                            .variable("i", i)
                            .variable("j", j)
                            .variable("bajando", numero(bajando))
                            .variable("girando", numero(girando)),
                        &mejor,
                        (i, j),
                        &[(i - 1, j), (i, j - 1)],
                    )
                    .paso(
                        "mejor[i][j] = bajando.max(girando);",
                        format!(
                            "Esquina ({i}, {j}). Desde arriba: {} + {} = {}. Desde la izquierda: {} + {} = {}. {veredicto} Se guarda {}.",
                            numero(mejor[i - 1][j]),
                            numero(abajo[i - 1][j]),
                            numero(bajando),
                            numero(mejor[i][j - 1]),
                            numero(derecha[i][j - 1]),
                            numero(girando),
                            numero(mejor[i][j])
                        ),
                    ),
                );
            }
        }
    }
    let longitud = mejor[filas - 1][columnas - 1];
    narracion.agregar(
        dibujar(
            Pizarra::nueva()
                .variable("filas", filas)
                .variable("columnas", columnas),
            &mejor,
            (filas - 1, columnas - 1),
            &[],
        )
        .paso(
            "let camino = camino_del_turista(&mejor, abajo, derecha);",
            format!(
                "La tabla está llena: el mejor recorrido pesa {}, lo que dice la última esquina. Para saber por dónde pasa, se camina hacia atrás desde ella, volviendo siempre por la calle que dio el mayor valor.",
                numero(longitud)
            ),
        ),
    );
    let mut i = filas - 1;
    let mut j = columnas - 1;
    let mut camino = vec![(i, j)];
    let mut recorridos = Vec::new();
    while i > 0 || j > 0 {
        let llego_bajando = i > 0
            && (j == 0 || mejor[i - 1][j] + abajo[i - 1][j] >= mejor[i][j - 1] + derecha[i][j - 1]);
        recorridos.push((i, j, llego_bajando));
        if llego_bajando {
            i -= 1;
        } else {
            j -= 1;
        }
        camino.push((i, j));
    }
    let ritmo_vuelta = Ritmo::nuevo(recorridos.len(), 5, 12);
    for (indice, &(desde_i, desde_j, bajando)) in recorridos.iter().enumerate() {
        if !ritmo_vuelta.toca(indice) {
            continue;
        }
        let hasta = if bajando {
            (desde_i - 1, desde_j)
        } else {
            (desde_i, desde_j - 1)
        };
        narracion.agregar(
            dibujar(
                Pizarra::nueva()
                    .variable("i", desde_i)
                    .variable("j", desde_j)
                    .variable("llego_bajando", bajando),
                &mejor,
                hasta,
                &camino[..=indice],
            )
            .paso(
                if bajando { "i -= 1;" } else { "j -= 1;" },
                if bajando {
                    format!(
                        "A la esquina ({desde_i}, {desde_j}) se llegó bajando desde ({}, {desde_j}), porque por ahí el valor era mayor o igual. Retrocedo hacia arriba.",
                        desde_i - 1
                    )
                } else {
                    format!(
                        "A la esquina ({desde_i}, {desde_j}) se llegó desde la izquierda, ({desde_i}, {}). Retrocedo hacia la izquierda.",
                        desde_j - 1
                    )
                },
            ),
        );
    }
    camino.reverse();
    let ultimo = Pizarra::nueva()
        .variable("longitud", numero(longitud))
        .lista(
            "camino (fila, columna)",
            camino
                .iter()
                .take(24)
                .map(|(i, j)| format!("({i}, {j})"))
                .collect(),
            &[],
        )
        .paso(
            "longitud: mejor[filas - 1][columnas - 1],",
            format!(
                "Dado vuelta, el recorrido sale de (0, 0): {}. Su peso total es {}. Mira «Qué significa» para compararlo con un turista que solo mira la calle de al lado.",
                flechas_del_camino(&camino),
                numero(longitud)
            ),
        );
    (narracion.terminar(ultimo), longitud, camino)
}

fn nombre_corto(nombre: &str) -> String {
    let visible: String = nombre.chars().take(22).collect();
    if nombre.chars().count() > 22 {
        format!("{visible}…")
    } else {
        visible
    }
}

fn tabla_de_nodos(
    pizarra: Pizarra,
    titulo: &str,
    columnas: &[&str],
    nodos: &[usize],
    nombres: &[String],
    foco: Option<usize>,
    valor: &dyn Fn(usize) -> Vec<String>,
) -> Pizarra {
    let posicion = foco
        .and_then(|nodo| nodos.iter().position(|&candidato| candidato == nodo))
        .unwrap_or(0);
    let desde = posicion
        .saturating_sub(8)
        .min(nodos.len().saturating_sub(12));
    let hasta = (desde + 12).min(nodos.len());
    let filas: Vec<Vec<String>> = nodos[desde..hasta]
        .iter()
        .map(|&nodo| {
            let mut fila = vec![nombre_corto(&nombres[nodo])];
            fila.extend(valor(nodo));
            fila
        })
        .collect();
    let marcadas: Vec<(usize, usize)> = foco
        .and_then(|nodo| {
            nodos[desde..hasta]
                .iter()
                .position(|&candidato| candidato == nodo)
        })
        .map(|fila| vec![(fila, 1)])
        .unwrap_or_default();
    pizarra.tabla(titulo, columnas, filas, &marcadas)
}

fn pasos_orden_topologico(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let grafo = leer_grafo(entrada, false)?;
    let pares: Vec<(usize, usize)> = grafo
        .aristas
        .iter()
        .map(|&(desde, hasta, _)| (desde, hasta))
        .collect();
    alineamiento::orden_topologico(grafo.nombres.len(), &pares)?;
    Ok(narrar_orden_topologico(&grafo.nombres, &pares).0)
}

fn narrar_orden_topologico(
    nombres: &[String],
    aristas: &[(usize, usize)],
) -> (Vec<Paso>, Vec<usize>) {
    let nodos = nombres.len();
    let todos: Vec<usize> = (0..nodos).collect();
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let mut entrantes = vec![0usize; nodos];
    let mut salientes = vec![Vec::new(); nodos];
    for &(desde, hasta) in aristas {
        salientes[desde].push(hasta);
        entrantes[hasta] += 1;
    }
    let vista_entrantes = |entrantes: &[usize], nodo: usize| vec![entrantes[nodo].to_string()];
    narracion.agregar(
        Pizarra::nueva()
            .variable("nodos", nodos)
            .variable("aristas", aristas.len())
            .lista(
                "flechas",
                aristas
                    .iter()
                    .take(14)
                    .map(|&(desde, hasta)| {
                        format!("{} → {}", nombre_corto(&nombres[desde]), nombre_corto(&nombres[hasta]))
                    })
                    .collect(),
                &[],
            )
            .paso(
                "let mut entrantes = vec![0; nodos];",
                format!(
                    "Hay {} nodos y {} flechas. Una tarea solo puede hacerse cuando ya se hicieron todas las que apuntan hacia ella. Por eso primero se cuenta, para cada nodo, cuántas flechas le llegan.",
                    miles(nodos),
                    miles(aristas.len())
                ),
            ),
    );
    let (primera_desde, primera_hasta) = aristas[0];
    narracion.agregar(
        Pizarra::nueva()
            .variable("nodos", nodos)
            .variable("desde", nombre_corto(&nombres[primera_desde]))
            .variable("hasta", nombre_corto(&nombres[primera_hasta]))
            .paso(
                "salientes[desde].push(hasta);",
                format!(
                    "Con la primera flecha, «{}» → «{}»: se anota que de «{}» se puede ir a «{}», y a «{}» le llega una flecha más.",
                    nombre_corto(&nombres[primera_desde]),
                    nombre_corto(&nombres[primera_hasta]),
                    nombre_corto(&nombres[primera_desde]),
                    nombre_corto(&nombres[primera_hasta]),
                    nombre_corto(&nombres[primera_hasta])
                ),
            ),
    );
    narracion.agregar(
        tabla_de_nodos(
            Pizarra::nueva().variable("nodos", nodos),
            "entrantes",
            &["nodo", "flechas que le llegan"],
            &todos,
            nombres,
            None,
            &|nodo| vista_entrantes(&entrantes, nodo),
        )
        .paso(
            "entrantes[hasta] += 1;",
            "Se repite con todas las flechas. Ahora cada nodo sabe cuántas flechas le llegan: los que tienen 0 no dependen de nadie y pueden hacerse primero.",
        ),
    );
    let mut sin_entrantes: VecDeque<usize> =
        (0..nodos).filter(|&nodo| entrantes[nodo] == 0).collect();
    let nombres_de_cola = |cola: &VecDeque<usize>| -> Vec<String> {
        cola.iter()
            .take(14)
            .map(|&nodo| nombre_corto(&nombres[nodo]))
            .collect()
    };
    narracion.agregar(
        Pizarra::nueva()
            .variable("nodos", nodos)
            .lista("sin_entrantes (la cola)", nombres_de_cola(&sin_entrantes), &[])
            .paso(
                "let mut sin_entrantes: VecDeque<usize> =",
                format!(
                    "La cola empieza con los nodos que no tienen flechas entrantes: {}. Ellos pueden ir primero.",
                    sin_entrantes
                        .iter()
                        .take(5)
                        .map(|&nodo| format!("«{}»", nombre_corto(&nombres[nodo])))
                        .collect::<Vec<_>>()
                        .join(", ")
                ),
            ),
    );
    let ritmo = Ritmo::nuevo(nodos, 8, 14);
    let mut orden: Vec<usize> = Vec::with_capacity(nodos);
    while let Some(nodo) = sin_entrantes.pop_front() {
        let indice = orden.len();
        orden.push(nodo);
        let detallar = indice < 8;
        let mostrar = detallar || ritmo.toca(indice);
        if mostrar {
            let ultimos: Vec<String> = orden
                .iter()
                .rev()
                .take(14)
                .rev()
                .map(|&hecho| nombre_corto(&nombres[hecho]))
                .collect();
            let resaltado = ultimos.len() - 1;
            narracion.agregar(
                Pizarra::nueva()
                    .variable("nodo", nombre_corto(&nombres[nodo]))
                    .variable("orden.len()", orden.len())
                    .lista("orden", ultimos, &[resaltado])
                    .paso(
                        "orden.push(nodo);",
                        format!(
                            "Saco «{}» de la cola y lo pongo en el orden (lugar {}). Ahora «quito» las {} flechas que salen de él.",
                            nombre_corto(&nombres[nodo]),
                            orden.len(),
                            salientes[nodo].len()
                        ),
                    ),
            );
        }
        let mut liberados = 0;
        for &siguiente in &salientes[nodo] {
            entrantes[siguiente] -= 1;
            let llega_a_cero = entrantes[siguiente] == 0;
            if llega_a_cero {
                sin_entrantes.push_back(siguiente);
                liberados += 1;
            }
            if detallar && liberados < 4 {
                let nota = if llega_a_cero {
                    format!(
                        "Se quita la flecha «{}» → «{}»: a «{}» ya no le falta nada, entra a la cola.",
                        nombre_corto(&nombres[nodo]),
                        nombre_corto(&nombres[siguiente]),
                        nombre_corto(&nombres[siguiente])
                    )
                } else {
                    format!(
                        "Se quita la flecha «{}» → «{}»: a «{}» todavía le faltan {} flechas, así que espera.",
                        nombre_corto(&nombres[nodo]),
                        nombre_corto(&nombres[siguiente]),
                        nombre_corto(&nombres[siguiente]),
                        entrantes[siguiente]
                    )
                };
                narracion.agregar(
                    tabla_de_nodos(
                        Pizarra::nueva()
                            .variable("nodo", nombre_corto(&nombres[nodo]))
                            .variable("siguiente", nombre_corto(&nombres[siguiente])),
                        "entrantes",
                        &["nodo", "flechas que le llegan"],
                        &todos,
                        nombres,
                        Some(siguiente),
                        &|otro| vista_entrantes(&entrantes, otro),
                    )
                    .lista(
                        "sin_entrantes (la cola)",
                        nombres_de_cola(&sin_entrantes),
                        &[],
                    )
                    .paso(
                        if llega_a_cero {
                            "sin_entrantes.push_back(siguiente);"
                        } else {
                            "entrantes[siguiente] -= 1;"
                        },
                        nota,
                    ),
                );
            }
        }
    }
    narracion.agregar(
        Pizarra::nueva()
            .variable("orden.len()", orden.len())
            .variable("nodos", nodos)
            .paso(
                "if orden.len() < nodos {",
                format!(
                    "Se colocaron los {} nodos: si hubiera un ciclo, algunos nunca llegarían a 0 flechas entrantes y faltarían. No falta ninguno, así que no hay ciclo.",
                    miles(nodos)
                ),
            ),
    );
    let ultimo = Pizarra::nueva()
        .lista(
            "orden",
            orden
                .iter()
                .take(24)
                .map(|&nodo| nombre_corto(&nombres[nodo]))
                .collect(),
            &[],
        )
        .paso(
            "Ok(orden)",
            "Este es el orden topológico: cada flecha apunta hacia adelante en la lista. Mira «Qué significa» para ver si es el único posible.",
        );
    (narracion.terminar(ultimo), orden)
}

fn pasos_camino_dag(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let grafo = leer_grafo(entrada, true)?;
    let fuente = buscar_nodo(&grafo, entrada, "fuente")?;
    let sumidero = buscar_nodo(&grafo, entrada, "sumidero")?;
    alineamiento::camino_mas_largo(grafo.nombres.len(), &grafo.aristas, fuente, sumidero)?;
    narrar_camino_dag(&grafo, fuente, sumidero).map(|resultado| resultado.0)
}

fn narrar_camino_dag(
    grafo: &Grafo,
    fuente: usize,
    sumidero: usize,
) -> Result<(Vec<Paso>, f64, Vec<usize>), String> {
    let nodos = grafo.nombres.len();
    let nombres = &grafo.nombres;
    let pares: Vec<(usize, usize)> = grafo
        .aristas
        .iter()
        .map(|&(desde, hasta, _)| (desde, hasta))
        .collect();
    let orden = alineamiento::orden_topologico(nodos, &pares)?;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let base = || {
        Pizarra::nueva()
            .variable("fuente", nombre_corto(&nombres[fuente]))
            .variable("sumidero", nombre_corto(&nombres[sumidero]))
    };
    narracion.agregar(
        base()
            .variable("nodos", nodos)
            .lista(
                "orden (topológico)",
                orden.iter().take(24).map(|&nodo| nombre_corto(&nombres[nodo])).collect(),
                &[],
            )
            .paso(
                "let orden = orden_topologico(nodos, &pares)?;",
                format!(
                    "Buscamos el camino de mayor peso de «{}» a «{}». Primero se ordenan los nodos para que cada flecha apunte hacia adelante. Así, cuando llegue el turno de un nodo, todos los que le llegan ya están calculados.",
                    nombre_corto(&nombres[fuente]),
                    nombre_corto(&nombres[sumidero])
                ),
            ),
    );
    let mut entrantes: Vec<Vec<(usize, f64)>> = vec![Vec::new(); nodos];
    for &(desde, hasta, peso) in &grafo.aristas {
        entrantes[hasta].push((desde, peso));
    }
    narracion.agregar(
        base()
            .variable("nodos", nodos)
            .variable("aristas", grafo.aristas.len())
            .paso(
                "entrantes[hasta].push((desde, peso));",
                "Se anota, para cada nodo, de dónde le llegan flechas y cuánto pesa cada una. Servirá para mirar «hacia atrás» desde cada nodo.",
            ),
    );
    let mut mejor: Vec<Option<f64>> = vec![None; nodos];
    let mut anterior: Vec<Option<usize>> = vec![None; nodos];
    mejor[fuente] = Some(0.0);
    let vista = |mejor: &[Option<f64>], anterior: &[Option<usize>], nodo: usize| {
        vec![
            mejor[nodo].map_or("None".to_string(), numero),
            anterior[nodo].map_or("None".to_string(), |previo| nombre_corto(&nombres[previo])),
        ]
    };
    narracion.agregar(
        tabla_de_nodos(
            base(),
            "mejor",
            &["nodo", "mejor", "anterior"],
            &orden,
            nombres,
            Some(fuente),
            &|nodo| vista(&mejor, &anterior, nodo),
        )
        .paso(
            "mejor[fuente] = Some(0.0);",
            format!(
                "El nodo de salida «{}» vale 0 porque ahí aún no hay recorrido. Los demás están vacíos (None) hasta que se calculen.",
                nombre_corto(&nombres[fuente])
            ),
        ),
    );
    let ritmo = Ritmo::nuevo(orden.len(), 12, 25);
    for (indice, &nodo) in orden.iter().enumerate() {
        let detallar = indice < 12;
        let mostrar = detallar || ritmo.toca(indice);
        let mut mejoro_con: Option<(usize, f64, f64)> = None;
        let mut intentos = Vec::new();
        for &(desde, peso) in &entrantes[nodo] {
            let Some(llegada) = mejor[desde] else {
                if detallar && mostrar {
                    narracion.agregar(
                        tabla_de_nodos(
                            base().variable("nodo", nombre_corto(&nombres[nodo])),
                            "mejor",
                            &["nodo", "mejor", "anterior"],
                            &orden,
                            nombres,
                            Some(nodo),
                            &|otro| vista(&mejor, &anterior, otro),
                        )
                        .paso(
                            "let Some(llegada) = mejor[desde] else {",
                            format!(
                                "A «{}» le llega una flecha desde «{}», pero a «{}» no se llega desde la salida: no sirve y se salta.",
                                nombre_corto(&nombres[nodo]),
                                nombre_corto(&nombres[desde]),
                                nombre_corto(&nombres[desde])
                            ),
                        ),
                    );
                }
                continue;
            };
            let candidato = llegada + peso;
            let mejora = mejor[nodo].is_none_or(|actual| candidato > actual);
            let valor_previo = mejor[nodo];
            if mejora {
                mejor[nodo] = Some(candidato);
                anterior[nodo] = Some(desde);
                mejoro_con = Some((desde, peso, candidato));
            }
            intentos.push((desde, candidato));
            if detallar && mostrar && intentos.len() <= 6 {
                narracion.agregar(
                    tabla_de_nodos(
                        base()
                            .variable("nodo", nombre_corto(&nombres[nodo]))
                            .variable("llegada", numero(llegada))
                            .variable("peso", numero(peso))
                            .variable("candidato", numero(candidato)),
                        "mejor",
                        &["nodo", "mejor", "anterior"],
                        &orden,
                        nombres,
                        Some(nodo),
                        &|otro| vista(&mejor, &anterior, otro),
                    )
                    .paso(
                        if mejora {
                            "mejor[nodo] = Some(candidato);"
                        } else {
                            "let candidato = llegada + peso;"
                        },
                        format!(
                            "Llegando a «{}» desde «{}»: {} + {} = {}. {}",
                            nombre_corto(&nombres[nodo]),
                            nombre_corto(&nombres[desde]),
                            numero(llegada),
                            numero(peso),
                            numero(candidato),
                            match (mejora, valor_previo) {
                                (true, None) =>
                                    "Es la primera forma de llegar: se anota.".to_string(),
                                (true, Some(previo)) =>
                                    format!("Supera los {} de antes: se anota.", numero(previo)),
                                (false, previo) => format!(
                                    "No supera los {} que ya había.",
                                    numero(previo.unwrap_or(0.0))
                                ),
                            }
                        ),
                    ),
                );
            }
        }
        if mostrar && !detallar {
            let nota = match mejoro_con {
                Some((desde, peso, candidato)) => format!(
                    "En «{}» se compararon {} flechas de entrada. La mejor viene de «{}» con peso {}, y da {}.",
                    nombre_corto(&nombres[nodo]),
                    intentos.len(),
                    nombre_corto(&nombres[desde]),
                    numero(peso),
                    numero(candidato)
                ),
                None => format!(
                    "«{}» no se alcanza desde la salida, así que queda vacío.",
                    nombre_corto(&nombres[nodo])
                ),
            };
            narracion.agregar(
                tabla_de_nodos(
                    base().variable("nodo", nombre_corto(&nombres[nodo])),
                    "mejor",
                    &["nodo", "mejor", "anterior"],
                    &orden,
                    nombres,
                    Some(nodo),
                    &|otro| vista(&mejor, &anterior, otro),
                )
                .paso("for &nodo in &orden {", nota),
            );
        }
    }
    let longitud =
        mejor[sumidero].ok_or("No hay ningún camino desde la fuente hasta el sumidero.")?;
    narracion.agregar(
        tabla_de_nodos(
            base(),
            "mejor",
            &["nodo", "mejor", "anterior"],
            &orden,
            nombres,
            Some(sumidero),
            &|nodo| vista(&mejor, &anterior, nodo),
        )
        .paso(
            "let longitud =",
            format!(
                "Todos los nodos están calculados. El mejor camino hasta «{}» pesa {}. Falta saber por dónde va: cada nodo recuerda de cuál vino (la columna «anterior»).",
                nombre_corto(&nombres[sumidero]),
                numero(longitud)
            ),
        ),
    );
    let mut camino = vec![sumidero];
    let mut nodo = sumidero;
    let mut pasos_atras = Vec::new();
    while let Some(previo) = anterior[nodo] {
        pasos_atras.push((nodo, previo));
        camino.push(previo);
        nodo = previo;
    }
    camino.reverse();
    let ritmo_vuelta = Ritmo::nuevo(pasos_atras.len(), 6, 12);
    for (indice, &(nodo, previo)) in pasos_atras.iter().enumerate() {
        if ritmo_vuelta.toca(indice) {
            let recorrido: Vec<String> = pasos_atras[..=indice]
                .iter()
                .rev()
                .take(12)
                .map(|&(_, antes)| nombre_corto(&nombres[antes]))
                .collect();
            narracion.agregar(
                base()
                    .variable("nodo", nombre_corto(&nombres[nodo]))
                    .variable("previo", nombre_corto(&nombres[previo]))
                    .lista("camino (de atrás hacia adelante)", recorrido.into_iter().rev().collect(), &[])
                    .paso(
                        "camino.push(previo);",
                        format!(
                            "«{}» recuerda que su mejor valor vino de «{}»: se agrega al camino y se sigue desde ahí.",
                            nombre_corto(&nombres[nodo]),
                            nombre_corto(&nombres[previo])
                        ),
                    ),
            );
        }
    }
    let ultimo = base()
        .variable("longitud", numero(longitud))
        .lista(
            "camino",
            camino.iter().take(24).map(|&nodo| nombre_corto(&nombres[nodo])).collect(),
            &[],
        )
        .paso(
            "camino.reverse();",
            format!(
                "Al darlo vuelta, el camino va de «{}» a «{}» pasando por {} nodos y pesa {}. Mira «Qué significa».",
                nombre_corto(&nombres[fuente]),
                nombre_corto(&nombres[sumidero]),
                camino.len(),
                numero(longitud)
            ),
        );
    Ok((narracion.terminar(ultimo), longitud, camino))
}

const FLECHA_DIAGONAL: u8 = 1;
const FLECHA_ARRIBA: u8 = 2;
const FLECHA_IZQUIERDA: u8 = 3;

fn dibujo_de_flecha(flecha: u8) -> &'static str {
    match flecha {
        FLECHA_DIAGONAL => "↖",
        FLECHA_ARRIBA => "↑",
        FLECHA_IZQUIERDA => "←",
        _ => "",
    }
}

fn nombre_de_flecha(flecha: u8) -> &'static str {
    match flecha {
        FLECHA_DIAGONAL => "diagonal",
        FLECHA_ARRIBA => "arriba",
        _ => "izquierda",
    }
}

fn pasos_subsecuencia_comun(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (v, _) = leer_secuencia(entrada, "v")?;
    let (w, _) = leer_secuencia(entrada, "w")?;
    entrada.natural("semilla")?;
    revisar_tamano(&v, &w, LIMITE_CUADRATICO, CONSEJO_LINEAL)?;
    revisar_narrado(&v, &w)?;
    Ok(narrar_subsecuencia(&v, &w).0)
}

fn narrar_subsecuencia(v: &[u8], w: &[u8]) -> (Vec<Paso>, Vec<u8>) {
    let (n, m) = (v.len(), w.len());
    let columnas = m + 1;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let mut largos = vec![0u32; (n + 1) * columnas];
    let mut flechas = vec![0u8; (n + 1) * columnas];
    let dibujar = |pizarra: Pizarra,
                   largos: &[u32],
                   flechas: &[u8],
                   hecho: (usize, usize),
                   centro: (usize, usize),
                   resaltadas: &[(usize, usize)]| {
        let valor = |i: usize, j: usize| {
            if i == 0 || j == 0 || i < hecho.0 || (i == hecho.0 && j <= hecho.1) {
                format!(
                    "{}{}",
                    largos[i * columnas + j],
                    dibujo_de_flecha(flechas[i * columnas + j])
                )
            } else {
                "·".to_string()
            }
        };
        Cuadricula {
            titulo: "largos (la flecha dice de dónde vino cada valor)",
            tamano: (n + 1, m + 1),
            rotulo_fila: &rotulo_de_letra(v),
            rotulo_columna: &rotulo_de_letra(w),
            valor: &valor,
        }
        .dibujar(pizarra, centro, resaltadas)
    };
    let base = || {
        Pizarra::nueva()
            .variable("v", corto(v))
            .variable("w", corto(w))
    };
    narracion.agregar(
        dibujar(base(), &largos, &flechas, (0, 0), (0, 0), &[]).paso(
            "let mut flechas = Flechas::nuevas(v.len() + 1, w.len() + 1);",
            format!(
                "Buscamos las letras que aparecen en el mismo orden en las dos secuencias (v tiene {} letras, w tiene {}). Se arma una tabla de {} × {}: la casilla (i, j) guardará el largo de la mejor subsecuencia común entre las primeras i letras de v y las primeras j de w. La fila 0 y la columna 0 valen 0 porque ahí una de las dos está vacía.",
                miles(n),
                miles(m),
                miles(n + 1),
                miles(m + 1)
            ),
        ),
    );
    let cada_fila = (n / 8).max(1);
    let mut coincidencias_mostradas = 0;
    let mut celdas_mostradas = 0;
    for i in 1..=n {
        for j in 1..=m {
            let coincide = v[i - 1] == w[j - 1];
            let diagonal = largos[(i - 1) * columnas + (j - 1)];
            let arriba = largos[(i - 1) * columnas + j];
            let izquierda = largos[i * columnas + (j - 1)];
            let (valor, flecha) = if coincide {
                (diagonal + 1, FLECHA_DIAGONAL)
            } else if arriba >= izquierda {
                (arriba, FLECHA_ARRIBA)
            } else {
                (izquierda, FLECHA_IZQUIERDA)
            };
            largos[i * columnas + j] = valor;
            flechas[i * columnas + j] = flecha;
            let detalle = celdas_mostradas < 6;
            let evento = coincide && coincidencias_mostradas < 8;
            if detalle || evento {
                celdas_mostradas += 1;
                if coincide {
                    coincidencias_mostradas += 1;
                }
                let (linea, nota) = if coincide {
                    (
                        "actual[j] = anterior[j - 1] + 1;",
                        format!(
                            "Las letras «{}» y «{}» son iguales: se pueden agregar a la subsecuencia. Se toma el valor de la diagonal, {diagonal}, y se suma 1: {valor}.",
                            v[i - 1] as char,
                            w[j - 1] as char
                        ),
                    )
                } else if flecha == FLECHA_ARRIBA {
                    (
                        "actual[j] = anterior[j];",
                        format!(
                            "«{}» y «{}» son distintas. Se hereda el mejor de los vecinos: de arriba viene {arriba} y de la izquierda {izquierda}. Gana arriba (empate o más): {valor}.",
                            v[i - 1] as char,
                            w[j - 1] as char
                        ),
                    )
                } else {
                    (
                        "actual[j] = actual[j - 1];",
                        format!(
                            "«{}» y «{}» son distintas. Se hereda el mejor de los vecinos: de arriba viene {arriba} y de la izquierda {izquierda}. Gana la izquierda: {valor}.",
                            v[i - 1] as char,
                            w[j - 1] as char
                        ),
                    )
                };
                narracion.agregar(
                    dibujar(
                        base()
                            .variable("i", i)
                            .variable("j", j)
                            .variable("v[i - 1]", v[i - 1] as char)
                            .variable("w[j - 1]", w[j - 1] as char)
                            .secuencia(v, &[(i - 1, i, Resalte::Actual)]),
                        &largos,
                        &flechas,
                        (i, j),
                        (i, j),
                        &[(i - 1, j - 1), (i - 1, j), (i, j - 1)],
                    )
                    .paso(linea, nota),
                );
            }
        }
        if i.is_multiple_of(cada_fila) && i < n {
            narracion.agregar(
                dibujar(
                    base().variable("i", i).secuencia(v, &[(i - 1, i, Resalte::Actual)]),
                    &largos,
                    &flechas,
                    (i, m),
                    (i, m),
                    &[],
                )
                .paso(
                    "std::mem::swap(&mut anterior, &mut actual);",
                    format!(
                        "Terminó la fila {i} de {n}: con las primeras {i} letras de v y todas las de w, la mejor subsecuencia común mide {}. La fila siguiente solo necesita esta, por eso se guarda para el próximo turno.",
                        largos[i * columnas + m]
                    ),
                ),
            );
        }
    }
    let total = largos[n * columnas + m];
    narracion.agregar(
        dibujar(base(), &largos, &flechas, (n, m), (n, m), &[]).paso(
            "let mut comun = Vec::with_capacity(anterior[w.len()]);",
            format!(
                "La tabla está llena. La esquina final dice el largo de la respuesta: {}. Para saber cuáles letras son, se sigue el camino de las flechas desde esa esquina hacia atrás.",
                miles(total as usize)
            ),
        ),
    );
    let mut comun = Vec::new();
    let (mut i, mut j) = (n, m);
    let mut movimientos = Vec::new();
    while i > 0 && j > 0 {
        let flecha = flechas[i * columnas + j];
        movimientos.push((i, j, flecha));
        match flecha {
            FLECHA_DIAGONAL => {
                comun.push(v[i - 1]);
                i -= 1;
                j -= 1;
            }
            FLECHA_ARRIBA => i -= 1,
            _ => j -= 1,
        }
    }
    let ritmo = Ritmo::nuevo(movimientos.len(), 6, 14);
    let mut guardadas = 0;
    for (indice, &(i, j, flecha)) in movimientos.iter().enumerate() {
        if flecha == FLECHA_DIAGONAL {
            guardadas += 1;
        }
        if !ritmo.toca(indice) {
            continue;
        }
        let (linea, nota) = match flecha {
            FLECHA_DIAGONAL => (
                "comun.push(v[i - 1]);",
                format!(
                    "En ({i}, {j}) la flecha es diagonal: «{}» está en las dos secuencias y entra a la respuesta. Subo en diagonal.",
                    v[i - 1] as char
                ),
            ),
            FLECHA_ARRIBA => (
                "ARRIBA => i -= 1,",
                format!(
                    "En ({i}, {j}) la flecha apunta arriba: «{}» de v no entra. Subo una fila.",
                    v[i - 1] as char
                ),
            ),
            _ => (
                "_ => j -= 1,",
                format!(
                    "En ({i}, {j}) la flecha apunta a la izquierda: «{}» de w no entra. Me muevo una columna a la izquierda.",
                    w[j - 1] as char
                ),
            ),
        };
        let reunidas: Vec<String> = comun[..guardadas.min(comun.len())]
            .iter()
            .rev()
            .take(20)
            .rev()
            .map(|&letra| (letra as char).to_string())
            .collect();
        let resaltado = reunidas.len().saturating_sub(1);
        narracion.agregar(
            dibujar(
                base()
                    .variable("i", i)
                    .variable("j", j)
                    .lista("comun (de atrás hacia adelante)", reunidas, &[resaltado])
                    .secuencia(v, &[(i - 1, i, Resalte::Actual)]),
                &largos,
                &flechas,
                (n, m),
                (i, j),
                &[],
            )
            .paso(linea, nota),
        );
    }
    comun.reverse();
    let ultimo = base()
        .lista(
            "comun",
            comun.iter().take(24).map(|&letra| (letra as char).to_string()).collect(),
            &[],
        )
        .paso(
            "comun.reverse();",
            format!(
                "Las letras se juntaron de atrás hacia adelante; al darlas vuelta queda la subsecuencia común más larga: «{}», de {} letras. Mira «Qué significa» para compararla con secuencias barajadas.",
                corto(&comun),
                miles(comun.len())
            ),
        );
    (narracion.terminar(ultimo), comun)
}

fn pasos_edicion(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (v, _) = leer_secuencia(entrada, "v")?;
    let (w, _) = leer_secuencia(entrada, "w")?;
    revisar_tamano(&v, &w, LIMITE_LINEAL, CONSEJO_LINEAL)?;
    revisar_narrado(&v, &w)?;
    Ok(narrar_edicion(&v, &w).0)
}

fn narrar_edicion(v: &[u8], w: &[u8]) -> (Vec<Paso>, usize) {
    let (n, m) = (v.len(), w.len());
    let columnas = m + 1;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let mut distancias = vec![0usize; (n + 1) * columnas];
    let dibujar = |pizarra: Pizarra,
                   distancias: &[usize],
                   hecho: (usize, usize),
                   centro: (usize, usize),
                   resaltadas: &[(usize, usize)]| {
        let valor = |i: usize, j: usize| {
            if i < hecho.0 || (i == hecho.0 && j <= hecho.1) {
                distancias[i * columnas + j].to_string()
            } else {
                "·".to_string()
            }
        };
        Cuadricula {
            titulo: "distancias (cambios mínimos entre los comienzos de v y de w)",
            tamano: (n + 1, m + 1),
            rotulo_fila: &rotulo_de_letra(v),
            rotulo_columna: &rotulo_de_letra(w),
            valor: &valor,
        }
        .dibujar(pizarra, centro, resaltadas)
    };
    let base = || {
        Pizarra::nueva()
            .variable("v", corto(v))
            .variable("w", corto(w))
    };
    for (j, celda) in distancias.iter_mut().enumerate().take(columnas) {
        *celda = j;
    }
    narracion.agregar(
        dibujar(base(), &distancias, (0, m), (0, 0), &[]).paso(
            "let mut anterior: Vec<usize> = (0..=w.len()).collect();",
            format!(
                "Queremos el menor número de cambios (sustituir, insertar o borrar una letra) para convertir v en w. La tabla {} × {} guarda, en (i, j), los cambios mínimos para convertir las primeras i letras de v en las primeras j de w. La fila 0 es 0, 1, 2…: convertir «nada» en j letras cuesta j inserciones.",
                miles(n + 1),
                miles(m + 1)
            ),
        ),
    );
    let cada_fila = (n / 8).max(1);
    let mut mostradas = 0;
    for i in 1..=n {
        distancias[i * columnas] = i;
        if i == 2 {
            narracion.agregar(
                dibujar(
                    base().variable("i", i),
                    &distancias,
                    (i, 0),
                    (i, 0),
                    &[],
                )
                .paso(
                    "actual[0] = i;",
                    format!(
                        "La primera casilla de cada fila vale i: convertir las primeras {i} letras de v en «nada» cuesta {i} borrados."
                    ),
                ),
            );
        }
        for j in 1..=m {
            let cambio = if v[i - 1] == w[j - 1] { 0 } else { 1 };
            let sustituir = distancias[(i - 1) * columnas + (j - 1)] + cambio;
            let borrar = distancias[(i - 1) * columnas + j] + 1;
            let insertar = distancias[i * columnas + (j - 1)] + 1;
            let mejor = sustituir.min(borrar).min(insertar);
            distancias[i * columnas + j] = mejor;
            let toca = mostradas < 6 || (i == 2 && j <= 3);
            if toca {
                mostradas += 1;
                let eleccion = if mejor == sustituir {
                    if cambio == 0 {
                        "Gana la diagonal sin costo: las letras coinciden."
                    } else {
                        "Gana la diagonal: se sustituye una letra por otra."
                    }
                } else if mejor == borrar {
                    "Gana borrar una letra de v."
                } else {
                    "Gana insertar una letra de w."
                };
                narracion.agregar(
                    dibujar(
                        base()
                            .variable("i", i)
                            .variable("j", j)
                            .variable("cambio", cambio)
                            .variable("sustituir", sustituir)
                            .variable("borrar", borrar)
                            .variable("insertar", insertar)
                            .variable("actual[j]", mejor)
                            .secuencia(v, &[(i - 1, i, Resalte::Actual)]),
                        &distancias,
                        (i, j),
                        (i, j),
                        &[(i - 1, j - 1), (i - 1, j), (i, j - 1)],
                    )
                    .paso(
                        "actual[j] = sustituir.min(borrar).min(insertar);",
                        format!(
                            "Letras «{}» y «{}» (cambio = {cambio}). Sustituir cuesta {sustituir}, borrar {borrar}, insertar {insertar}. Se queda el menor: {mejor}. {eleccion}",
                            v[i - 1] as char,
                            w[j - 1] as char
                        ),
                    ),
                );
            }
        }
        if i.is_multiple_of(cada_fila) && i < n {
            narracion.agregar(
                dibujar(
                    base().variable("i", i).secuencia(v, &[(i - 1, i, Resalte::Actual)]),
                    &distancias,
                    (i, m),
                    (i, m),
                    &[],
                )
                .paso(
                    "std::mem::swap(&mut anterior, &mut actual);",
                    format!(
                        "Terminó la fila {i} de {n}. Convertir las primeras {i} letras de v en todo w cuesta {} cambios. Esta fila se guarda como «anterior» para calcular la siguiente.",
                        miles(distancias[i * columnas + m])
                    ),
                ),
            );
        }
    }
    let distancia = distancias[n * columnas + m];
    let ultimo = dibujar(base(), &distancias, (n, m), (n, m), &[]).paso(
        "anterior[w.len()]",
        format!(
            "La última casilla responde: se necesitan {} cambios como mínimo para convertir v en w. Mira «Qué significa» para ver de qué tipo son.",
            miles(distancia)
        ),
    );
    (narracion.terminar(ultimo), distancia)
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Variante {
    Local,
    Ajustado,
    Solapado,
}

struct DatosDeTabla {
    v: Vec<u8>,
    w: Vec<u8>,
    puntuar: Puntuacion,
    puntuacion: &'static str,
    indel: i64,
}

fn leer_tabla(entrada: &Entrada) -> Result<DatosDeTabla, String> {
    let (v, _) = leer_secuencia(entrada, "v")?;
    let (w, _) = leer_secuencia(entrada, "w")?;
    let (puntuar, puntuacion) = leer_puntuacion(entrada)?;
    let indel = leer_castigo(entrada, "indel")?;
    revisar_tamano(&v, &w, LIMITE_CUADRATICO, CONSEJO_CUADRATICO)?;
    revisar_narrado(&v, &w)?;
    Ok(DatosDeTabla {
        v,
        w,
        puntuar,
        puntuacion,
        indel,
    })
}

fn pasos_local(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let datos = leer_tabla(entrada)?;
    Ok(narrar_tabla(Variante::Local, &datos).0)
}

fn pasos_ajustado(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let datos = leer_tabla(entrada)?;
    Ok(narrar_tabla(Variante::Ajustado, &datos).0)
}

fn pasos_solapado(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let datos = leer_tabla(entrada)?;
    Ok(narrar_tabla(Variante::Solapado, &datos).0)
}

fn elegir_mejor_narrado(diagonal: i64, arriba: i64, izquierda: i64) -> (i64, u8) {
    if diagonal >= arriba && diagonal >= izquierda {
        (diagonal, FLECHA_DIAGONAL)
    } else if arriba >= izquierda {
        (arriba, FLECHA_ARRIBA)
    } else {
        (izquierda, FLECHA_IZQUIERDA)
    }
}

fn posicion_del_mayor(valores: &[i64]) -> usize {
    let mut mejor = 0;
    for (posicion, &valor) in valores.iter().enumerate() {
        if valor > valores[mejor] {
            mejor = posicion;
        }
    }
    mejor
}

fn lista_de_valores(valores: &[i64], centro: usize) -> (Vec<String>, usize) {
    let desde = centro
        .saturating_sub(6)
        .min(valores.len().saturating_sub(13));
    let hasta = (desde + 13).min(valores.len());
    (
        (desde..hasta)
            .map(|posicion| format!("{posicion}: {}", valores[posicion]))
            .collect(),
        centro - desde,
    )
}

fn narrar_tabla(variante: Variante, datos: &DatosDeTabla) -> (Vec<Paso>, Alineamiento) {
    let (v, w) = (datos.v.as_slice(), datos.w.as_slice());
    let (puntuar, indel) = (datos.puntuar, datos.indel);
    let (n, m) = (v.len(), w.len());
    let columnas = m + 1;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let mut valores = vec![0i64; (n + 1) * columnas];
    let mut flechas = vec![0u8; (n + 1) * columnas];
    let dibujar = |pizarra: Pizarra,
                   valores: &[i64],
                   flechas: &[u8],
                   hecho: (usize, usize),
                   centro: (usize, usize),
                   resaltadas: &[(usize, usize)]| {
        let valor = |i: usize, j: usize| {
            if i == 0 || j == 0 || i < hecho.0 || (i == hecho.0 && j <= hecho.1) {
                format!(
                    "{}{}",
                    valores[i * columnas + j],
                    dibujo_de_flecha(flechas[i * columnas + j])
                )
            } else {
                "·".to_string()
            }
        };
        Cuadricula {
            titulo: "tabla de puntajes (la flecha dice de dónde vino cada valor)",
            tamano: (n + 1, m + 1),
            rotulo_fila: &rotulo_de_letra(v),
            rotulo_columna: &rotulo_de_letra(w),
            valor: &valor,
        }
        .dibujar(pizarra, centro, resaltadas)
    };
    let base = || {
        Pizarra::nueva()
            .variable("v", corto(v))
            .variable("w", corto(w))
            .variable("indel", indel)
    };
    let presentacion = match variante {
        Variante::Local => format!(
            "Buscamos el tramo de v y el tramo de w que mejor se parecen, ignorando el resto. Se arma una tabla de {} × {}: la casilla (i, j) guarda el mejor puntaje de un alineamiento que termina en la letra i de v y la letra j de w. Si un puntaje saliera negativo, se reinicia en 0: un alineamiento local puede empezar donde quiera.",
            miles(n + 1),
            miles(m + 1)
        ),
        Variante::Ajustado => format!(
            "Buscamos dónde encaja w entera dentro de v. Se arma una tabla de {} × {}: la casilla (i, j) guarda el mejor puntaje de alinear las primeras j letras de w con un tramo de v que termina en su letra i. Empezar en cualquier lugar de v no cuesta nada.",
            miles(n + 1),
            miles(m + 1)
        ),
        Variante::Solapado => format!(
            "Buscamos cuánto se superponen el final de v y el comienzo de w. Se arma una tabla de {} × {}: la casilla (i, j) guarda el mejor puntaje de alinear un tramo de v que termina en su letra i con las primeras j letras de w. Empezar en cualquier lugar de v no cuesta nada.",
            miles(n + 1),
            miles(m + 1)
        ),
    };
    narracion.agregar(
        dibujar(
            base()
                .variable("puntuar", datos.puntuacion)
                .secuencia(v, &[]),
            &valores,
            &flechas,
            (0, 0),
            (0, 0),
            &[],
        )
        .paso(
            "let mut flechas = Flechas::nuevas(v.len() + 1, w.len() + 1);",
            presentacion,
        ),
    );
    if variante == Variante::Local {
        narracion.agregar(
            dibujar(base(), &valores, &flechas, (0, m), (0, 0), &[]).paso(
                "if modo != Modo::Local {",
                "En el modo local no se castiga el comienzo: la fila 0 y la columna 0 valen 0, porque el tramo parecido puede empezar en cualquier letra de las dos secuencias.",
            ),
        );
    } else {
        for j in 1..=m {
            valores[j] = valores[j - 1] - indel;
            flechas[j] = FLECHA_IZQUIERDA;
        }
        narracion.agregar(
            dibujar(base(), &valores, &flechas, (0, m), (0, m.min(3)), &[]).paso(
                "anterior[j] = anterior[j - 1] - indel;",
                format!(
                    "Fila 0: si todavía no se usó ninguna letra de v, las primeras j letras de w van contra huecos y cada una resta {indel}. La columna 0 vale 0 porque en v se puede empezar donde sea."
                ),
            ),
        );
    }
    let cada_fila = (n / 8).max(1);
    let mut mejor_local = (0i64, 0usize, 0usize);
    let mut mostradas = 0;
    let mut maximos_mostrados = 0;
    for i in 1..=n {
        for j in 1..=m {
            let par = puntuar(v[i - 1], w[j - 1]) as i64;
            let diagonal = valores[(i - 1) * columnas + (j - 1)] + par;
            let arriba = valores[(i - 1) * columnas + j] - indel;
            let izquierda = valores[i * columnas + (j - 1)] - indel;
            let (puntaje, flecha) = elegir_mejor_narrado(diagonal, arriba, izquierda);
            let reinicia = variante == Variante::Local && puntaje <= 0;
            let nuevo_maximo = variante == Variante::Local && puntaje > mejor_local.0;
            if reinicia {
                valores[i * columnas + j] = 0;
                flechas[i * columnas + j] = 0;
            } else {
                valores[i * columnas + j] = puntaje;
                flechas[i * columnas + j] = flecha;
            }
            if nuevo_maximo {
                mejor_local = (puntaje, i, j);
            }
            let detalle = mostradas < 6;
            let evento = nuevo_maximo && maximos_mostrados < 12;
            if !(detalle || evento) {
                continue;
            }
            mostradas += 1;
            if nuevo_maximo {
                maximos_mostrados += 1;
            }
            let resumen = format!(
                "Casilla ({i}, {j}), letras «{}» de v y «{}» de w. Emparejarlas vale {par}: {diagonal} en total. Un hueco en w (viniendo de arriba) da {arriba}. Un hueco en v (viniendo de la izquierda) da {izquierda}.",
                v[i - 1] as char,
                w[j - 1] as char
            );
            let (linea, cierre) = if reinicia {
                (
                    "if modo == Modo::Local && puntaje <= 0 {",
                    format!(
                        "La mejor opción vale {puntaje}, que no es positiva: el modo local reinicia aquí en 0 y no deja flecha, porque seguir con un tramo malo solo restaría."
                    ),
                )
            } else if nuevo_maximo {
                (
                    "mejor_local = (puntaje, i, j);",
                    format!(
                        "Gana {} con {puntaje}. Es el mayor puntaje visto hasta ahora en toda la tabla: se anota como el mejor local, porque ahí puede terminar el tramo parecido.",
                        nombre_de_flecha(flecha)
                    ),
                )
            } else {
                (
                    "let (puntaje, flecha) = elegir_mejor(diagonal, arriba, izquierda);",
                    format!(
                        "Gana {} con {puntaje}: ese valor y su flecha se guardan en la casilla.",
                        nombre_de_flecha(flecha)
                    ),
                )
            };
            let mut pizarra = base()
                .variable("i", i)
                .variable("j", j)
                .variable("diagonal", diagonal)
                .variable("arriba", arriba)
                .variable("izquierda", izquierda)
                .variable("puntaje", puntaje);
            if variante == Variante::Local {
                pizarra = pizarra.variable(
                    "mejor_local",
                    format!("({}, {}, {})", mejor_local.0, mejor_local.1, mejor_local.2),
                );
            }
            let mut marcas = vec![(i - 1, i, Resalte::Actual)];
            if nuevo_maximo {
                marcas.push((i - 1, i, Resalte::Hallado));
            }
            narracion.agregar(
                dibujar(
                    pizarra.secuencia(v, &marcas),
                    &valores,
                    &flechas,
                    (i, j),
                    (i, j),
                    &[(i - 1, j - 1), (i - 1, j), (i, j - 1)],
                )
                .paso(linea, format!("{resumen} {cierre}")),
            );
        }
        if i.is_multiple_of(cada_fila) && i < n {
            let cierre = match variante {
                Variante::Local => format!(
                    "El mejor puntaje local hasta ahora es {} en ({}, {}).",
                    mejor_local.0, mejor_local.1, mejor_local.2
                ),
                _ => format!(
                    "Usando w completa y terminando en esta letra de v, el puntaje es {}.",
                    valores[i * columnas + m]
                ),
            };
            narracion.agregar(
                dibujar(
                    base().variable("i", i).secuencia(v, &[(i - 1, i, Resalte::Actual)]),
                    &valores,
                    &flechas,
                    (i, m),
                    (i, m),
                    &[],
                )
                .paso(
                    "std::mem::swap(&mut anterior, &mut actual);",
                    format!(
                        "Terminó la fila {i} de {n}. Se guarda como «anterior» para calcular la siguiente. {cierre}"
                    ),
                ),
            );
        }
    }
    let ultima_columna: Vec<i64> = (0..=n).map(|i| valores[i * columnas + m]).collect();
    let ultima_fila: Vec<i64> = valores[n * columnas..].to_vec();
    let (puntaje, fin_v, fin_w) = match variante {
        Variante::Local => mejor_local,
        Variante::Ajustado => {
            let i = posicion_del_mayor(&ultima_columna);
            (ultima_columna[i], i, m)
        }
        Variante::Solapado => {
            let j = posicion_del_mayor(&ultima_fila);
            (ultima_fila[j], n, j)
        }
    };
    let pizarra_final =
        |pizarra: Pizarra| dibujar(pizarra, &valores, &flechas, (n, m), (fin_v, fin_w), &[]);
    match variante {
        Variante::Local => narracion.agregar(
            pizarra_final(
                base()
                    .variable(
                        "mejor_local",
                        format!("({puntaje}, {fin_v}, {fin_w})"),
                    )
                    .secuencia(v, &[(fin_v.saturating_sub(1), fin_v, Resalte::Hallado)]),
            )
            .paso(
                "Modo::Local => mejor_local,",
                if puntaje > 0 {
                    format!(
                        "La tabla está llena. El mejor puntaje de toda ella es {puntaje}, en la casilla ({fin_v}, {fin_w}): ahí termina el mejor alineamiento local. Se recorre hacia atrás hasta llegar a una casilla sin flecha."
                    )
                } else {
                    "La tabla está llena y ninguna casilla fue positiva: con estos puntajes ningún tramo de v se parece a uno de w, así que no hay alineamiento local.".to_string()
                },
            ),
        ),
        Variante::Ajustado => {
            let (elementos, resaltado) = lista_de_valores(&ultima_columna, fin_v);
            narracion.agregar(
                pizarra_final(
                    base()
                        .lista("ultima_columna (fila: puntaje)", elementos, &[resaltado])
                        .secuencia(v, &[(fin_v.saturating_sub(1), fin_v, Resalte::Hallado)]),
                )
                .paso(
                    "let i = posicion_del_maximo(&ultima_columna);",
                    format!(
                        "La última columna dice, para cada fila, cuánto vale encajar w entera terminando en esa letra de v. El mayor valor es {puntaje}, en la fila {fin_v}: ahí termina el mejor encaje."
                    ),
                ),
            );
        }
        Variante::Solapado => {
            let (elementos, resaltado) = lista_de_valores(&ultima_fila, fin_w);
            narracion.agregar(
                pizarra_final(
                    base()
                        .lista("ultima_fila (columna: puntaje)", elementos, &[resaltado])
                        .secuencia(v, &[]),
                )
                .paso(
                    "let j = posicion_del_maximo(&ultima_fila);",
                    format!(
                        "La última fila dice, para cada columna, cuánto vale que el final de v coincida con las primeras j letras de w. El mayor valor es {puntaje}, con j = {fin_w}: el solapamiento usa las primeras {fin_w} letras de w."
                    ),
                ),
            );
        }
    }
    let mut i = fin_v;
    let mut j = fin_w;
    let mut fila_v = Vec::new();
    let mut fila_w = Vec::new();
    let mut movimientos = Vec::new();
    loop {
        let flecha = flechas[i * columnas + j];
        if !(1..=3).contains(&flecha) {
            break;
        }
        movimientos.push((i, j, flecha));
        match flecha {
            FLECHA_DIAGONAL => {
                fila_v.push(v[i - 1]);
                fila_w.push(w[j - 1]);
                i -= 1;
                j -= 1;
            }
            FLECHA_ARRIBA => {
                fila_v.push(v[i - 1]);
                fila_w.push(b'-');
                i -= 1;
            }
            _ => {
                fila_v.push(b'-');
                fila_w.push(w[j - 1]);
                j -= 1;
            }
        }
    }
    let ritmo = Ritmo::nuevo(movimientos.len(), 6, 14);
    for (indice, &(desde_i, desde_j, flecha)) in movimientos.iter().enumerate() {
        if !ritmo.toca(indice) {
            continue;
        }
        let (linea, nota) = match flecha {
            FLECHA_DIAGONAL => (
                "DIAGONAL => {",
                format!(
                    "En ({desde_i}, {desde_j}) la flecha es diagonal: «{}» de v va con «{}» de w. Subo en diagonal.",
                    v[desde_i - 1] as char,
                    w[desde_j - 1] as char
                ),
            ),
            FLECHA_ARRIBA => (
                "ARRIBA => {",
                format!(
                    "En ({desde_i}, {desde_j}) la flecha viene de arriba: «{}» de v va frente a un hueco. Subo una fila.",
                    v[desde_i - 1] as char
                ),
            ),
            _ => (
                "IZQUIERDA => {",
                format!(
                    "En ({desde_i}, {desde_j}) la flecha viene de la izquierda: un hueco va frente a «{}» de w. Me muevo a la izquierda.",
                    w[desde_j - 1] as char
                ),
            ),
        };
        let recorrido = &movimientos[..=indice];
        let hasta_ahora = recorrido.len();
        let marcas = [
            (desde_i.saturating_sub(1), desde_i, Resalte::Actual),
            (desde_i, fin_v, Resalte::Hallado),
        ];
        let camino: Vec<(usize, usize)> = recorrido.iter().map(|&(a, b, _)| (a, b)).collect();
        narracion.agregar(
            dibujar(
                base()
                    .variable("i", desde_i)
                    .variable("j", desde_j)
                    .variable("columnas recorridas", hasta_ahora)
                    .secuencia(v, &marcas),
                &valores,
                &flechas,
                (n, m),
                (desde_i, desde_j),
                &camino,
            )
            .paso(linea, nota),
        );
    }
    narracion.agregar(
        dibujar(
            base()
                .variable("i", i)
                .variable("j", j)
                .secuencia(v, &[(i, fin_v, Resalte::Hallado)]),
            &valores,
            &flechas,
            (n, m),
            (i, j),
            &[],
        )
        .paso(
            "_ => break,",
            match variante {
                Variante::Local => format!(
                    "La casilla ({i}, {j}) no tiene flecha: vale 0 y ahí empieza el tramo parecido, que arranca después de la letra {i} de v y de la letra {j} de w."
                ),
                Variante::Ajustado => format!(
                    "Se llegó a la columna 0, la que no tiene flecha: w empieza a encajar justo después de la letra {i} de v."
                ),
                Variante::Solapado => format!(
                    "Se llegó a la columna 0, la que no tiene flecha: el solapamiento empieza después de la letra {i} de v."
                ),
            },
        ),
    );
    fila_v.reverse();
    fila_w.reverse();
    let ultimo = base()
        .variable("puntaje", puntaje)
        .lista(
            "alineamiento",
            vec![corto(&fila_v), corto(&fila_w)],
            &[],
        )
        .paso(
            "fila_v.reverse();",
            format!(
                "Las letras se juntaron de atrás hacia adelante; al darlas vuelta queda el alineamiento de {} columnas, con puntaje {puntaje}. Mira «Qué significa» para interpretarlo.",
                miles(fila_v.len())
            ),
        );
    let alineado = Alineamiento {
        puntaje,
        fila_v,
        fila_w,
        inicio_v: i,
        inicio_w: j,
    };
    (narracion.terminar(ultimo), alineado)
}

const HUECO_ABIERTO: u8 = 4;
const HUECO_ALARGADO: u8 = 5;
const NADA: i64 = i64::MIN / 4;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Nivel {
    Inferior,
    Medio,
    Superior,
}

fn elegir_hueco_narrado(alargar: i64, abrir: i64) -> (i64, u8) {
    if alargar >= abrir {
        (alargar, HUECO_ALARGADO)
    } else {
        (abrir, HUECO_ABIERTO)
    }
}

struct DatosAfin {
    v: Vec<u8>,
    w: Vec<u8>,
    puntuar: Puntuacion,
    apertura: i64,
    extension: i64,
}

fn leer_afin(entrada: &Entrada) -> Result<DatosAfin, String> {
    let (v, _) = leer_secuencia(entrada, "v")?;
    let (w, _) = leer_secuencia(entrada, "w")?;
    let (puntuar, _) = leer_puntuacion(entrada)?;
    let apertura = leer_castigo(entrada, "apertura")?;
    let extension = leer_castigo(entrada, "extension")?;
    if extension > apertura {
        return Err(format!(
            "La extensión ({extension}) no puede ser mayor que la apertura ({apertura}): si alargar un hueco costara más que abrir otro, ya no serían huecos afines."
        ));
    }
    revisar_tamano(&v, &w, LIMITE_CUADRATICO, CONSEJO_CUADRATICO)?;
    revisar_narrado(&v, &w)?;
    Ok(DatosAfin {
        v,
        w,
        puntuar,
        apertura,
        extension,
    })
}

fn pasos_afin(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let datos = leer_afin(entrada)?;
    Ok(narrar_afin(&datos).0)
}

fn texto_de_nivel(nivel: Nivel) -> &'static str {
    match nivel {
        Nivel::Inferior => "inferior",
        Nivel::Medio => "medio",
        Nivel::Superior => "superior",
    }
}

fn narrar_afin(datos: &DatosAfin) -> (Vec<Paso>, Alineamiento) {
    let (v, w) = (datos.v.as_slice(), datos.w.as_slice());
    let (puntuar, apertura) = (datos.puntuar, datos.apertura);
    let extension = datos.extension.min(apertura);
    let (n, m) = (v.len(), w.len());
    let columnas = m + 1;
    let tamano = (n + 1) * columnas;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let mut inferior = vec![NADA; tamano];
    let mut medio = vec![NADA; tamano];
    let mut superior = vec![NADA; tamano];
    let mut flechas_inferior = vec![0u8; tamano];
    let mut flechas_medio = vec![0u8; tamano];
    let mut flechas_superior = vec![0u8; tamano];
    let dibujar = |pizarra: Pizarra,
                   medio: &[i64],
                   flechas: &[u8],
                   hecho: (usize, usize),
                   centro: (usize, usize),
                   resaltadas: &[(usize, usize)]| {
        let valor = |i: usize, j: usize| {
            if i < hecho.0 || (i == hecho.0 && j <= hecho.1) {
                format!(
                    "{}{}",
                    texto_de_celda(medio[i * columnas + j]),
                    dibujo_de_flecha(flechas[i * columnas + j])
                )
            } else {
                "·".to_string()
            }
        };
        Cuadricula {
            titulo: "tabla del medio (la mejor, la que termina con letra contra letra)",
            tamano: (n + 1, m + 1),
            rotulo_fila: &rotulo_de_letra(v),
            rotulo_columna: &rotulo_de_letra(w),
            valor: &valor,
        }
        .dibujar(pizarra, centro, resaltadas)
    };
    let base = || {
        Pizarra::nueva()
            .variable("v", corto(v))
            .variable("w", corto(w))
            .variable("apertura", apertura)
            .variable("extension", extension)
    };
    narracion.agregar(
        dibujar(base().secuencia(v, &[]), &medio, &flechas_medio, (0, 0), (0, 0), &[]).paso(
            "let extension = extension.min(apertura);",
            format!(
                "Abrir un hueco cuesta {apertura} y cada letra más del mismo hueco cuesta solo {extension}. Para recordar si un hueco ya está abierto se usan tres tablas: la de abajo (termina con una letra de v frente a un hueco), la de arriba (termina con un hueco frente a una letra de w) y la del medio, la mejor de todas, que es la que da el puntaje."
            ),
        ),
    );
    medio[0] = 0;
    for j in 1..columnas {
        let (valor, flecha) =
            elegir_hueco_narrado(superior[j - 1] - extension, medio[j - 1] - apertura);
        superior[j] = valor;
        flechas_superior[j] = flecha;
        medio[j] = valor;
        flechas_medio[j] = FLECHA_IZQUIERDA;
    }
    narracion.agregar(
        dibujar(base(), &medio, &flechas_medio, (0, m), (0, m.min(3)), &[]).paso(
            "let (valor, flecha) = elegir_hueco(superior[j - 1] - extension, medio[j - 1] - apertura);",
            format!(
                "Fila 0: alinear solo huecos contra las primeras j letras de w. La primera letra abre el hueco ({}) y cada letra más lo alarga ({} cada una). La columna 0 se llena igual, hacia abajo.",
                -apertura,
                -extension
            ),
        ),
    );
    let cada_fila = (n / 8).max(1);
    let mut mostradas = 0;
    let mut alargues_inferior = 0;
    let mut alargues_superior = 0;
    for i in 1..=n {
        let (valor, flecha) = elegir_hueco_narrado(
            inferior[(i - 1) * columnas] - extension,
            medio[(i - 1) * columnas] - apertura,
        );
        inferior[i * columnas] = valor;
        flechas_inferior[i * columnas] = flecha;
        medio[i * columnas] = valor;
        flechas_medio[i * columnas] = FLECHA_ARRIBA;
        for j in 1..=m {
            let celda = i * columnas + j;
            let (valor_inferior, flecha_inferior) = elegir_hueco_narrado(
                inferior[celda - columnas] - extension,
                medio[celda - columnas] - apertura,
            );
            inferior[celda] = valor_inferior;
            flechas_inferior[celda] = flecha_inferior;
            let (valor_superior, flecha_superior) =
                elegir_hueco_narrado(superior[celda - 1] - extension, medio[celda - 1] - apertura);
            superior[celda] = valor_superior;
            flechas_superior[celda] = flecha_superior;
            let par = puntuar(v[i - 1], w[j - 1]) as i64;
            let diagonal = medio[celda - columnas - 1] + par;
            let (valor_medio, flecha_medio) =
                elegir_mejor_narrado(diagonal, valor_inferior, valor_superior);
            medio[celda] = valor_medio;
            flechas_medio[celda] = flecha_medio;
            let alarga_abajo = flecha_inferior == HUECO_ALARGADO
                && valor_inferior > IMPOSIBLE_VISIBLE
                && alargues_inferior < 3;
            let alarga_arriba = flecha_superior == HUECO_ALARGADO
                && valor_superior > IMPOSIBLE_VISIBLE
                && alargues_superior < 3;
            let detalle = mostradas < 5;
            if !(detalle || alarga_abajo || alarga_arriba) {
                continue;
            }
            mostradas += 1;
            let (linea, nota) = if alarga_abajo && !detalle {
                alargues_inferior += 1;
                (
                    "let (valor_inferior, flecha_inferior) = elegir_hueco(",
                    format!(
                        "Casilla ({i}, {j}): para terminar con «{}» de v frente a un hueco hay dos caminos. Alargar un hueco que ya estaba abierto arriba cuesta solo {extension} y da {}; abrir uno nuevo desde la tabla del medio cuesta {apertura} y daría {}. Alargar es más barato: por eso los huecos largos salen baratos.",
                        v[i - 1] as char,
                        texto_de_celda(inferior[celda - columnas] - extension),
                        texto_de_celda(medio[celda - columnas] - apertura)
                    ),
                )
            } else if alarga_arriba && !detalle {
                alargues_superior += 1;
                (
                    "let (valor_superior, flecha_superior) =",
                    format!(
                        "Casilla ({i}, {j}): para terminar con un hueco frente a «{}» de w, alargar el hueco de la izquierda cuesta solo {extension} y da {}; abrir uno nuevo cuesta {apertura} y daría {}. Gana alargar.",
                        w[j - 1] as char,
                        texto_de_celda(superior[celda - 1] - extension),
                        texto_de_celda(medio[celda - 1] - apertura)
                    ),
                )
            } else {
                (
                    "let (valor_medio, flecha_medio) =",
                    format!(
                        "Casilla ({i}, {j}), letras «{}» y «{}». Emparejarlas vale {par}: {} en total. Terminar con un hueco en w da {}, y con un hueco en v da {}. Gana {}: {}.",
                        v[i - 1] as char,
                        w[j - 1] as char,
                        texto_de_celda(diagonal),
                        texto_de_celda(valor_inferior),
                        texto_de_celda(valor_superior),
                        nombre_de_flecha(flecha_medio),
                        valor_medio
                    ),
                )
            };
            narracion.agregar(
                dibujar(
                    base()
                        .variable("i", i)
                        .variable("j", j)
                        .variable("valor_inferior", texto_de_celda(valor_inferior))
                        .variable("valor_superior", texto_de_celda(valor_superior))
                        .variable("diagonal", texto_de_celda(diagonal))
                        .variable("valor_medio", valor_medio)
                        .lista(
                            "las tres tablas en esta casilla",
                            vec![
                                format!("inferior: {}", texto_de_celda(valor_inferior)),
                                format!("medio: {valor_medio}"),
                                format!("superior: {}", texto_de_celda(valor_superior)),
                            ],
                            &[1],
                        )
                        .secuencia(v, &[(i - 1, i, Resalte::Actual)]),
                    &medio,
                    &flechas_medio,
                    (i, j),
                    (i, j),
                    &[(i - 1, j - 1), (i - 1, j), (i, j - 1)],
                )
                .paso(linea, nota),
            );
        }
        if i.is_multiple_of(cada_fila) && i < n {
            narracion.agregar(
                dibujar(
                    base().variable("i", i).secuencia(v, &[(i - 1, i, Resalte::Actual)]),
                    &medio,
                    &flechas_medio,
                    (i, m),
                    (i, m),
                    &[],
                )
                .paso(
                    "std::mem::swap(&mut medio, &mut medio_anterior);",
                    format!(
                        "Terminó la fila {i} de {n}. Alinear las primeras {i} letras de v con todo w vale {} hasta aquí.",
                        texto_de_celda(medio[i * columnas + m])
                    ),
                ),
            );
        }
    }
    let puntaje = medio[n * columnas + m];
    narracion.agregar(
        dibujar(base(), &medio, &flechas_medio, (n, m), (n, m), &[]).paso(
            "puntaje: medio[w.len()],",
            format!(
                "Las tres tablas están llenas. El puntaje del mejor alineamiento está en la esquina final de la tabla del medio: {puntaje}. Ahora se reconstruye el alineamiento hacia atrás, pasando de una tabla a otra cuando hay un hueco."
            ),
        ),
    );
    let mut fila_v = Vec::new();
    let mut fila_w = Vec::new();
    let (mut i, mut j) = (n, m);
    let mut nivel = Nivel::Medio;
    let mut eventos: Vec<(usize, usize, Nivel, Nivel, char, char)> = Vec::new();
    while i > 0 || j > 0 {
        let celda = i * columnas + j;
        let desde = nivel;
        match nivel {
            Nivel::Medio => match flechas_medio[celda] {
                FLECHA_DIAGONAL => {
                    eventos.push((i, j, desde, nivel, v[i - 1] as char, w[j - 1] as char));
                    fila_v.push(v[i - 1]);
                    fila_w.push(w[j - 1]);
                    i -= 1;
                    j -= 1;
                }
                FLECHA_ARRIBA => {
                    nivel = Nivel::Inferior;
                    eventos.push((i, j, desde, nivel, ' ', ' '));
                }
                _ => {
                    nivel = Nivel::Superior;
                    eventos.push((i, j, desde, nivel, ' ', ' '));
                }
            },
            Nivel::Inferior => {
                fila_v.push(v[i - 1]);
                fila_w.push(b'-');
                if flechas_inferior[celda] == HUECO_ABIERTO {
                    nivel = Nivel::Medio;
                }
                eventos.push((i, j, desde, nivel, v[i - 1] as char, '-'));
                i -= 1;
            }
            Nivel::Superior => {
                fila_v.push(b'-');
                fila_w.push(w[j - 1]);
                if flechas_superior[celda] == HUECO_ABIERTO {
                    nivel = Nivel::Medio;
                }
                eventos.push((i, j, desde, nivel, '-', w[j - 1] as char));
                j -= 1;
            }
        }
    }
    let ritmo = Ritmo::nuevo(eventos.len(), 6, 14);
    let mut cambios_mostrados = 0;
    for (indice, &(desde_i, desde_j, desde, hasta, letra_v, letra_w)) in eventos.iter().enumerate()
    {
        let cambio = desde != hasta;
        let importante = cambio && cambios_mostrados < 10;
        if !(ritmo.toca(indice) || importante) {
            continue;
        }
        if cambio {
            cambios_mostrados += 1;
        }
        let (linea, nota) = match (desde, hasta) {
            (Nivel::Medio, Nivel::Medio) => (
                "Nivel::Medio => match flechas_medio.leer(i, j) {",
                format!(
                    "En ({desde_i}, {desde_j}) de la tabla del medio la flecha es diagonal: «{letra_v}» va con «{letra_w}»."
                ),
            ),
            (Nivel::Medio, Nivel::Inferior) => (
                "ARRIBA => nivel = Nivel::Inferior,",
                format!(
                    "En ({desde_i}, {desde_j}) la mejor opción termina con un hueco frente a una letra de v: me paso a la tabla inferior para ver cómo se formó ese hueco."
                ),
            ),
            (Nivel::Medio, _) => (
                "_ => nivel = Nivel::Superior,",
                format!(
                    "En ({desde_i}, {desde_j}) la mejor opción termina con un hueco frente a una letra de w: me paso a la tabla superior."
                ),
            ),
            (Nivel::Inferior, destino) => (
                "Nivel::Inferior => {",
                format!(
                    "Tabla inferior en ({desde_i}, {desde_j}): «{letra_v}» de v va frente a un hueco. {}",
                    if destino == Nivel::Medio {
                        "Aquí se abrió el hueco, así que vuelvo a la tabla del medio."
                    } else {
                        "El hueco se alargó desde arriba: sigo en esta tabla."
                    }
                ),
            ),
            (Nivel::Superior, destino) => (
                "Nivel::Superior => {",
                format!(
                    "Tabla superior en ({desde_i}, {desde_j}): un hueco va frente a «{letra_w}» de w. {}",
                    if destino == Nivel::Medio {
                        "Aquí se abrió el hueco, así que vuelvo a la tabla del medio."
                    } else {
                        "El hueco se alargó desde la izquierda: sigo en esta tabla."
                    }
                ),
            ),
        };
        let camino: Vec<(usize, usize)> = eventos[..=indice].iter().map(|e| (e.0, e.1)).collect();
        narracion.agregar(
            dibujar(
                base()
                    .variable("i", desde_i)
                    .variable("j", desde_j)
                    .variable("nivel", texto_de_nivel(desde))
                    .secuencia(v, &[(desde_i.saturating_sub(1), desde_i, Resalte::Actual)]),
                &medio,
                &flechas_medio,
                (n, m),
                (desde_i, desde_j),
                &camino,
            )
            .paso(linea, nota),
        );
    }
    fila_v.reverse();
    fila_w.reverse();
    let ultimo = base()
        .variable("puntaje", puntaje)
        .lista(
            "alineamiento",
            vec![corto(&fila_v), corto(&fila_w)],
            &[],
        )
        .paso(
            "fila_v.reverse();",
            format!(
                "Las letras se juntaron de atrás hacia adelante; al darlas vuelta queda el alineamiento de {} columnas con puntaje {puntaje}. Mira «Qué significa» para interpretarlo.",
                miles(fila_v.len())
            ),
        );
    let alineado = Alineamiento {
        puntaje,
        fila_v,
        fila_w,
        inicio_v: 0,
        inicio_w: 0,
    };
    (narracion.terminar(ultimo), alineado)
}

const CONSEJO_NARRADO: &str = "Para ver el paso a paso usa un tramo más corto.";

fn muestra_de_columna(valores: &[i64], resaltada: Option<usize>) -> (Vec<String>, Vec<usize>) {
    let ultimo = valores.len() - 1;
    let mut filas: Vec<usize> = (0..8).map(|k| k * ultimo / 7).collect();
    if let Some(fila) = resaltada {
        filas.push(fila);
    }
    filas.sort_unstable();
    filas.dedup();
    let resaltos = resaltada
        .and_then(|fila| filas.iter().position(|&candidata| candidata == fila))
        .map(|posicion| vec![posicion])
        .unwrap_or_default();
    (
        filas
            .iter()
            .map(|&fila| format!("fila {fila}: {}", valores[fila]))
            .collect(),
        resaltos,
    )
}

fn pasos_arista_central(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (v, _) = leer_secuencia(entrada, "v")?;
    let (w, _) = leer_secuencia(entrada, "w")?;
    let (puntuar, _) = leer_puntuacion(entrada)?;
    let indel = leer_castigo(entrada, "indel")?;
    revisar_tamano(&v, &w, LIMITE_LINEAL, CONSEJO_LINEAL)?;
    revisar_tamano(&v, &w, LIMITE_CUADRATICO, CONSEJO_NARRADO)?;
    narrar_arista_central(&v, &w, puntuar, indel)
        .map(|resultado| resultado.0)
        .ok_or_else(|| "Hacen falta letras en alguna de las dos secuencias.".to_string())
}

fn narrar_arista_central(
    v: &[u8],
    w: &[u8],
    puntuar: Puntuacion,
    indel: i64,
) -> Option<(Vec<Paso>, alineamiento::AristaCentral)> {
    let n = v.len();
    let m = w.len();
    let centro = m / 2;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let base = || {
        Pizarra::nueva()
            .variable("v", corto(v))
            .variable("w", corto(w))
            .variable("indel", indel)
    };
    narracion.agregar(
        base()
            .variable("centro", centro)
            .secuencia(
                w,
                &[(0, centro, Resalte::Comparado), (centro, m, Resalte::Actual)],
            )
            .paso(
                "let centro = w.len() / 2;",
                format!(
                    "Se parte w por la mitad: la primera mitad son sus {centro} primeras letras y la segunda mitad las {} restantes. Queremos saber por qué fila cruza el mejor alineamiento la columna {centro} de la tabla, sin llenar la tabla entera: bastan dos recorridos de una columna a la vez, uno desde el comienzo y otro desde el final.",
                    miles(m - centro)
                ),
            ),
    );
    let mut columna: Vec<i64> = (0..=n).map(|i| -(i as i64) * indel).collect();
    let mut izquierda = vec![0i64; n + 1];
    let (texto, resaltos) = muestra_de_columna(&columna, None);
    narracion.agregar(
        base()
            .lista("columna (fila: puntaje)", texto, &resaltos)
            .paso(
                "let mut columna: Vec<i64> = (0..=v.len()).map(|i| -(i as i64) * indel).collect();",
                format!(
                    "Desde el comienzo: la columna 0 es v entera contra huecos, así que la fila i vale -i × {indel}. Desde ella se calcula la columna siguiente, y así hasta la columna {centro}."
                ),
            ),
    );
    let ritmo = Ritmo::nuevo(centro, 2, 5);
    for (indice, &letra_w) in w[..centro].iter().enumerate() {
        std::mem::swap(&mut columna, &mut izquierda);
        columna[0] = izquierda[0] - indel;
        for i in 1..=n {
            let diagonal = izquierda[i - 1] + puntuar(v[i - 1], letra_w) as i64;
            let horizontal = izquierda[i] - indel;
            let vertical = columna[i - 1] - indel;
            columna[i] = diagonal.max(horizontal).max(vertical);
            if indice == 0 && i <= 3 {
                narracion.agregar(
                    base()
                        .variable("i", i)
                        .variable("diagonal", diagonal)
                        .variable("horizontal", horizontal)
                        .variable("vertical", vertical)
                        .variable("columna[i]", columna[i])
                        .secuencia(w, &[(indice, indice + 1, Resalte::Actual)])
                        .paso(
                            "let vertical = columna[i - 1] - indel;",
                            format!(
                                "Primera columna, fila {i}: emparejar «{}» con «{}» da {diagonal}; venir de la izquierda (hueco en v) da {horizontal}; venir de arriba (hueco en w) da {vertical}. Se queda el mayor: {}.",
                                v[i - 1] as char,
                                letra_w as char,
                                columna[i]
                            ),
                        ),
                );
            }
        }
        if ritmo.toca(indice) && centro > 0 {
            let (texto, resaltos) = muestra_de_columna(&columna, None);
            narracion.agregar(
                base()
                    .variable("letra_w", letra_w as char)
                    .lista("columna (fila: puntaje)", texto, &resaltos)
                    .secuencia(w, &[(indice, indice + 1, Resalte::Actual)])
                    .paso(
                        "for &letra_w in w {",
                        format!(
                            "Se terminó la columna {} de {centro} (letra «{}» de w). Solo se guardan la columna actual y la anterior; las demás se olvidan.",
                            indice + 1,
                            letra_w as char
                        ),
                    ),
            );
        }
    }
    let desde_la_fuente = columna;
    let (texto, resaltos) = muestra_de_columna(&desde_la_fuente, None);
    narracion.agregar(
        base()
            .lista("desde_la_fuente (fila: puntaje)", texto, &resaltos)
            .paso(
                "let desde_la_fuente = puntajes_desde_la_fuente(v, &w[..centro], puntajes, indel);",
                format!(
                    "Listo el primer recorrido: para cada fila de la columna {centro}, cuánto vale el mejor camino que llega hasta ahí desde el comienzo de la tabla."
                ),
            ),
    );
    let mut columna: Vec<i64> = (0..=n).map(|i| -((n - i) as i64) * indel).collect();
    let mut derecha = vec![0i64; n + 1];
    let (texto, resaltos) = muestra_de_columna(&columna, None);
    narracion.agregar(
        base()
            .lista("columna (fila: puntaje)", texto, &resaltos)
            .paso(
                "let mut columna: Vec<i64> = (0..=n).map(|i| -((n - i) as i64) * indel).collect();",
                "Ahora el segundo recorrido, al revés: desde el final de la tabla hacia la columna central. En la última columna, cada fila solo puede seguir con huecos hasta el final.",
            ),
    );
    let segunda = &w[centro..];
    let ritmo = Ritmo::nuevo(segunda.len(), 2, 5);
    for (indice, &letra_w) in segunda.iter().rev().enumerate() {
        std::mem::swap(&mut columna, &mut derecha);
        columna[n] = derecha[n] - indel;
        for i in (0..n).rev() {
            let diagonal = derecha[i + 1] + puntuar(v[i], letra_w) as i64;
            let horizontal = derecha[i] - indel;
            let vertical = columna[i + 1] - indel;
            columna[i] = diagonal.max(horizontal).max(vertical);
            if indice == 0 && n - i <= 3 {
                narracion.agregar(
                    base()
                        .variable("i", i)
                        .variable("diagonal", diagonal)
                        .variable("horizontal", horizontal)
                        .variable("vertical", vertical)
                        .variable("columna[i]", columna[i])
                        .secuencia(w, &[(m - 1 - indice, m - indice, Resalte::Actual)])
                        .paso(
                            "let vertical = columna[i + 1] - indel;",
                            format!(
                                "Última columna, fila {i} (se va de abajo hacia arriba): emparejar «{}» con «{}» da {diagonal}; seguir por la derecha da {horizontal}; seguir hacia abajo da {vertical}. Se queda el mayor: {}.",
                                v[i] as char,
                                letra_w as char,
                                columna[i]
                            ),
                        ),
                );
            }
        }
        if ritmo.toca(indice) {
            let (texto, resaltos) = muestra_de_columna(&columna, None);
            narracion.agregar(
                base()
                    .variable("letra_w", letra_w as char)
                    .lista("columna (fila: puntaje)", texto, &resaltos)
                    .secuencia(w, &[(m - 1 - indice, m - indice, Resalte::Actual)])
                    .paso(
                        "for &letra_w in w.iter().rev() {",
                        format!(
                            "Se terminó la columna {} de {} contando desde el final (letra «{}» de w).",
                            indice + 1,
                            segunda.len(),
                            letra_w as char
                        ),
                    ),
            );
        }
    }
    let hasta_el_sumidero = columna;
    let hasta_el_sumidero_a_la_derecha = derecha;
    let (texto, resaltos) = muestra_de_columna(&hasta_el_sumidero, None);
    narracion.agregar(
        base()
            .lista("hasta_el_sumidero (fila: puntaje)", texto, &resaltos)
            .paso(
                "puntajes_hasta_el_sumidero(v, &w[centro..], puntajes, indel);",
                format!(
                    "Listo el segundo recorrido: para cada fila de la columna {centro}, cuánto vale el mejor camino que sigue desde ahí hasta el final de la tabla."
                ),
            ),
    );
    let longitudes: Vec<i64> = desde_la_fuente
        .iter()
        .zip(&hasta_el_sumidero)
        .map(|(antes, despues)| antes + despues)
        .collect();
    let fila = posicion_del_mayor(&longitudes);
    let empates = longitudes
        .iter()
        .filter(|&&longitud| longitud == longitudes[fila])
        .count();
    let (texto, resaltos) = muestra_de_columna(&longitudes, Some(fila));
    narracion.agregar(
        base()
            .lista("longitudes (fila: puntaje)", texto, &resaltos)
            .paso(
                "let longitudes: Vec<i64> = desde_la_fuente",
                "Se suman las dos cosas: el mejor camino que pasa por cada fila de la columna central vale lo que cuesta llegar ahí más lo que cuesta salir. Esa suma es el puntaje de un alineamiento completo que cruza la columna por esa fila.",
            ),
    );
    let (texto, resaltos) = muestra_de_columna(&longitudes, Some(fila));
    narracion.agregar(
        base()
            .variable("fila", fila)
            .variable("longitudes[fila]", longitudes[fila])
            .lista("longitudes (fila: puntaje)", texto, &resaltos)
            .paso(
                "let fila = posicion_del_maximo(&longitudes);",
                format!(
                    "El mayor puntaje es {}, en la fila {fila}: por ahí cruza la columna central el mejor alineamiento{}. Ese valor es el puntaje del alineamiento global completo.",
                    longitudes[fila],
                    if empates > 1 {
                        format!(" (otras {} filas empatan, así que hay más de uno)", empates - 1)
                    } else {
                        String::new()
                    }
                ),
            ),
    );
    let resto = hasta_el_sumidero[fila];
    narracion.agregar(
        base()
            .variable("fila", fila)
            .variable("resto", resto)
            .paso(
                "let resto = hasta_el_sumidero[fila];",
                format!(
                    "Desde el nodo ({fila}, {centro}) falta recorrer un tramo que vale {resto}. Se busca cuál de los tres movimientos posibles (diagonal, derecha, abajo) conserva justo ese valor: ese movimiento es la arista central."
                ),
            ),
    );
    let hay_columna_siguiente = centro < m;
    let hay_fila_siguiente = fila < n;
    let diagonal_vale = hay_fila_siguiente
        && hay_columna_siguiente
        && puntuar(v[fila], w[centro]) as i64 + hasta_el_sumidero_a_la_derecha[fila + 1] == resto;
    let derecha_vale =
        hay_columna_siguiente && hasta_el_sumidero_a_la_derecha[fila] - indel == resto;
    let abajo_vale = hay_fila_siguiente && hasta_el_sumidero[fila + 1] - indel == resto;
    let (hasta, linea, nota) = if diagonal_vale {
        (
            (fila + 1, centro + 1),
            "(fila + 1, centro + 1)",
            format!(
                "Emparejar «{}» de v con «{}» de w y seguir desde ({}, {}) conserva el valor {resto}: la arista central es diagonal.",
                v[fila] as char,
                w[centro] as char,
                fila + 1,
                centro + 1
            ),
        )
    } else if derecha_vale {
        (
            (fila, centro + 1),
            "(fila, centro + 1)",
            format!(
                "Poner «{}» de w frente a un hueco y seguir desde ({fila}, {}) conserva el valor {resto}: la arista central es horizontal.",
                w[centro] as char,
                centro + 1
            ),
        )
    } else if abajo_vale {
        (
            (fila + 1, centro),
            "(fila + 1, centro)",
            format!(
                "Poner «{}» de v frente a un hueco y seguir desde ({}, {centro}) conserva el valor {resto}: la arista central es vertical.",
                v[fila] as char,
                fila + 1
            ),
        )
    } else {
        return None;
    };
    narracion.agregar(
        base()
            .variable("fila", fila)
            .variable("hay_columna_siguiente", hay_columna_siguiente)
            .variable("hay_fila_siguiente", hay_fila_siguiente)
            .paso(linea, nota),
    );
    let arista = alineamiento::AristaCentral {
        desde: (fila, centro),
        hasta,
        longitudes,
    };
    let ultimo = base()
        .variable("desde", format!("({}, {})", arista.desde.0, arista.desde.1))
        .variable("hasta", format!("({}, {})", arista.hasta.0, arista.hasta.1))
        .paso(
            "Some(AristaCentral {",
            format!(
                "La arista central va de ({}, {}) a ({}, {}) y es parte de un alineamiento óptimo. En el alineamiento en espacio lineal se repite esto en cada mitad. Mira «Qué significa».",
                arista.desde.0, arista.desde.1, arista.hasta.0, arista.hasta.1
            ),
        );
    Some((narracion.terminar(ultimo), arista))
}

struct Mitades<'a> {
    v: &'a [u8],
    w: &'a [u8],
    puntuar: Puntuacion,
    indel: i64,
    narracion: Narracion,
    fila_v: Vec<u8>,
    fila_w: Vec<u8>,
    llamadas: usize,
    cada: usize,
    emisiones: usize,
}

type Tramo = (usize, usize);

impl Mitades<'_> {
    fn pizarra(&self, tramo_v: Tramo, tramo_w: Tramo) -> Pizarra {
        let columnas: Vec<String> = self
            .fila_v
            .iter()
            .zip(&self.fila_w)
            .rev()
            .take(10)
            .rev()
            .map(|(&a, &b)| format!("{}/{}", a as char, b as char))
            .collect();
        let resaltadas: Vec<usize> = columnas.len().checked_sub(1).into_iter().collect();
        Pizarra::nueva()
            .variable("v.len()", tramo_v.1 - tramo_v.0)
            .variable("w.len()", tramo_w.1 - tramo_w.0)
            .lista("columnas ya escritas (las últimas)", columnas, &resaltadas)
            .secuencia(self.v, &[(tramo_v.0, tramo_v.1, Resalte::Actual)])
    }

    fn mitades(&mut self, tramo_v: Tramo, tramo_w: Tramo, profundidad: usize) {
        let total_v = self.v;
        let total_w = self.w;
        let v = &total_v[tramo_v.0..tramo_v.1];
        let w = &total_w[tramo_w.0..tramo_w.1];
        let indice = self.llamadas;
        self.llamadas += 1;
        let mostrar = indice < 6 || indice.is_multiple_of(self.cada);
        if v.len() <= 1 || w.len() <= 1 {
            let directo = alineamiento::alineamiento_global(v, w, self.puntuar, self.indel);
            if mostrar {
                let paso = self.pizarra(tramo_v, tramo_w).paso(
                    "let directo = alinear(v, w, puntajes, indel, Modo::Global);",
                    format!(
                        "El trozo v[{}..{}] × w[{}..{}] ya es tan pequeño (una de las dos tiene una letra o ninguna) que se alinea de una vez con una tabla completa: aporta {} columnas. Profundidad {profundidad}.",
                        tramo_v.0,
                        tramo_v.1,
                        tramo_w.0,
                        tramo_w.1,
                        directo.fila_v.len()
                    ),
                );
                self.narracion.agregar(paso);
            }
            self.fila_v.extend_from_slice(&directo.fila_v);
            self.fila_w.extend_from_slice(&directo.fila_w);
            return;
        }
        let Some(arista) = alineamiento::arista_central(v, w, self.puntuar, self.indel) else {
            return;
        };
        let (i, j) = arista.desde;
        let (siguiente_i, siguiente_j) = arista.hasta;
        if mostrar {
            let paso = self
                .pizarra(tramo_v, tramo_w)
                .variable("i", i)
                .variable("j", j)
                .variable("siguiente_i", siguiente_i)
                .variable("siguiente_j", siguiente_j)
                .paso(
                    "let Some(arista) = buscar_arista_central(v, w, puntajes, indel) else {",
                    format!(
                        "Trozo v[{}..{}] × w[{}..{}] (profundidad {profundidad}). Su arista central cruza la columna {j} por la fila {i} y sigue hacia ({siguiente_i}, {siguiente_j}). Esa arista es una columna del alineamiento; la parte de arriba-izquierda y la de abajo-derecha se resuelven aparte, y cada una es más chica.",
                        tramo_v.0, tramo_v.1, tramo_w.0, tramo_w.1
                    ),
                );
            self.narracion.agregar(paso);
        }
        self.mitades(
            (tramo_v.0, tramo_v.0 + i),
            (tramo_w.0, tramo_w.0 + j),
            profundidad + 1,
        );
        let letra_v = if siguiente_i > i { v[i] } else { b'-' };
        let letra_w = if siguiente_j > j { w[j] } else { b'-' };
        self.fila_v.push(letra_v);
        self.fila_w.push(letra_w);
        if self.emisiones < 4 {
            self.emisiones += 1;
            let paso = self
                .pizarra(tramo_v, tramo_w)
                .variable("i", i)
                .variable("j", j)
                .paso(
                    "fila_v.push(if siguiente_i > i { v[i] } else { b'-' });",
                    format!(
                        "Ya se resolvió la parte de arriba-izquierda de este trozo. Ahora se escribe la columna de la arista central: «{}» de v frente a «{}» de w.",
                        letra_v as char, letra_w as char
                    ),
                );
            self.narracion.agregar(paso);
        }
        self.mitades(
            (tramo_v.0 + siguiente_i, tramo_v.1),
            (tramo_w.0 + siguiente_j, tramo_w.1),
            profundidad + 1,
        );
    }
}

fn pasos_espacio_lineal(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (v, _) = leer_secuencia(entrada, "v")?;
    let (w, _) = leer_secuencia(entrada, "w")?;
    let (puntuar, _) = leer_puntuacion(entrada)?;
    let indel = leer_castigo(entrada, "indel")?;
    revisar_tamano(&v, &w, LIMITE_LINEAL, CONSEJO_LINEAL)?;
    revisar_tamano(&v, &w, LIMITE_CUADRATICO, CONSEJO_NARRADO)?;
    Ok(narrar_espacio_lineal(&v, &w, puntuar, indel).0)
}

fn narrar_espacio_lineal(
    v: &[u8],
    w: &[u8],
    puntuar: Puntuacion,
    indel: i64,
) -> (Vec<Paso>, Alineamiento) {
    let mut estado = Mitades {
        v,
        w,
        puntuar,
        indel,
        narracion: Narracion::nueva(MAXIMO_PASOS),
        fila_v: Vec::with_capacity(v.len() + w.len()),
        fila_w: Vec::with_capacity(v.len() + w.len()),
        llamadas: 0,
        cada: ((v.len() + w.len()) / 4).max(1),
        emisiones: 0,
    };
    let paso = estado.pizarra((0, v.len()), (0, w.len())).paso(
        "let mut fila_v = Vec::with_capacity(v.len() + w.len());",
        format!(
            "Hay que alinear v ({} letras) con w ({} letras) sin guardar la tabla entera. La idea: se busca la arista central, es decir, por dónde cruza el mejor alineamiento la columna del medio. Esa arista es una columna del resultado. Luego se hace lo mismo con el trozo de arriba-izquierda y con el de abajo-derecha, cada vez más pequeños, hasta que alineen de una vez.",
            miles(v.len()),
            miles(w.len())
        ),
    );
    estado.narracion.agregar(paso);
    estado.mitades((0, v.len()), (0, w.len()), 0);
    let Mitades {
        mut narracion,
        fila_v,
        fila_w,
        llamadas,
        ..
    } = estado;
    narracion.agregar(
        Pizarra::nueva()
            .variable("fila_v.len()", fila_v.len())
            .variable("fila_w.len()", fila_w.len())
            .paso(
                "alinear_por_mitades(v, w, &puntajes, indel, &mut fila_v, &mut fila_w);",
                format!(
                    "Terminaron todas las mitades: se resolvieron {} trozos en total y las filas quedaron con {} columnas, escritas de izquierda a derecha. En ningún momento se guardó más que un par de columnas de la tabla.",
                    miles(llamadas),
                    miles(fila_v.len())
                ),
            ),
    );
    let puntaje = alineamiento::puntuar_alineamiento(&fila_v, &fila_w, puntuar, indel);
    let ultimo = Pizarra::nueva()
        .variable("puntaje", puntaje)
        .lista("alineamiento", vec![corto(&fila_v), corto(&fila_w)], &[])
        .paso(
            "let puntaje = puntuar_alineamiento(&fila_v, &fila_w, puntuar, indel);",
            format!(
                "Se suma el puntaje de las columnas, una por una, y da {puntaje}: el mismo del alineamiento global, con mucha menos memoria. Mira «Qué significa» para interpretarlo."
            ),
        );
    let alineado = Alineamiento {
        puntaje,
        fila_v,
        fila_w,
        inicio_v: 0,
        inicio_w: 0,
    };
    (narracion.terminar(ultimo), alineado)
}

const MOVIMIENTOS_DE_TRES: &[(usize, usize, usize)] = &[
    (1, 1, 1),
    (1, 1, 0),
    (1, 0, 1),
    (0, 1, 1),
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
];

fn pasos_tres(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (u, _) = leer_secuencia(entrada, "u")?;
    let (v, _) = leer_secuencia(entrada, "v")?;
    let (w, _) = leer_secuencia(entrada, "w")?;
    let celdas = (u.len() + 1)
        .saturating_mul(v.len() + 1)
        .saturating_mul(w.len() + 1);
    if celdas > LIMITE_TRES {
        return Err(format!(
            "El cubo tendría {} celdas y el límite es {}: con tres secuencias el trabajo crece con el producto de los tres largos. Usa tramos de unas 200 letras.",
            miles(celdas),
            miles(LIMITE_TRES)
        ));
    }
    if celdas > LIMITE_NARRADO {
        return Err(format!(
            "El cubo tendría {} celdas y el paso a paso lo guarda entero para enseñártelo: el límite es {}. Usa tramos de unas 120 letras.",
            miles(celdas),
            miles(LIMITE_NARRADO)
        ));
    }
    Ok(narrar_tres(&u, &v, &w).0)
}

fn llegada_de_tres(
    cubo: &[u32],
    dimensiones: (usize, usize),
    secuencias: (&[u8], &[u8], &[u8]),
    celda: (usize, usize, usize),
    movimiento: (usize, usize, usize),
) -> Option<u32> {
    let (u, v, w) = secuencias;
    let (i, j, k) = celda;
    let (paso_u, paso_v, paso_w) = movimiento;
    if i < paso_u || j < paso_v || k < paso_w {
        return None;
    }
    let (ancho, fondo) = dimensiones;
    let antes = cubo[((i - paso_u) * ancho + (j - paso_v)) * fondo + (k - paso_w)];
    let tres_iguales = movimiento == (1, 1, 1) && u[i - 1] == v[j - 1] && v[j - 1] == w[k - 1];
    Some(if tres_iguales { antes + 1 } else { antes })
}

fn nombre_de_movimiento(movimiento: (usize, usize, usize)) -> String {
    let partes: Vec<&str> = [
        (movimiento.0, "u"),
        (movimiento.1, "v"),
        (movimiento.2, "w"),
    ]
    .iter()
    .filter(|(avanza, _)| *avanza == 1)
    .map(|(_, nombre)| *nombre)
    .collect();
    format!(
        "({}, {}, {}) avanza {}",
        movimiento.0,
        movimiento.1,
        movimiento.2,
        partes.join(" y ")
    )
}

fn narrar_tres(u: &[u8], v: &[u8], w: &[u8]) -> (Vec<Paso>, alineamiento::AlineamientoTriple) {
    let (nu, nv, nw) = (u.len(), v.len(), w.len());
    let (ancho, fondo) = (nv + 1, nw + 1);
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let mut cubo = vec![0u32; (nu + 1) * ancho * fondo];
    let secuencias = (u, v, w);
    let base = || {
        Pizarra::nueva()
            .variable("u", corto(u))
            .variable("v", corto(v))
            .variable("w", corto(w))
    };
    let dibujar = |pizarra: Pizarra,
                   cubo: &[u32],
                   capa: usize,
                   hecho: (usize, usize),
                   centro: (usize, usize),
                   resaltadas: &[(usize, usize)]| {
        let valor = |j: usize, k: usize| {
            if j < hecho.0 || (j == hecho.0 && k <= hecho.1) {
                cubo[(capa * ancho + j) * fondo + k].to_string()
            } else {
                "·".to_string()
            }
        };
        Cuadricula {
            titulo: "capa i del cubo (filas: letras de v, columnas: letras de w)",
            tamano: (nv + 1, nw + 1),
            rotulo_fila: &rotulo_de_letra(v),
            rotulo_columna: &rotulo_de_letra(w),
            valor: &valor,
        }
        .dibujar(pizarra, centro, resaltadas)
    };
    narracion.agregar(
        dibujar(base().secuencia(u, &[]), &cubo, 0, (0, 0), (0, 0), &[]).paso(
            "let mut cubo = Cubo::de_ceros(u.len() + 1, v.len() + 1, w.len() + 1);",
            format!(
                "Con tres secuencias la tabla pasa a ser un cubo de {} × {} × {} = {} celdas: la celda (i, j, k) guarda cuántas columnas con las tres letras iguales caben usando las primeras i letras de u, j de v y k de w. Abajo ves una capa del cubo, la de la letra i de u.",
                miles(nu + 1),
                miles(nv + 1),
                miles(nw + 1),
                miles((nu + 1) * ancho * fondo)
            ),
        ),
    );
    let cada_capa = (nu / 8).max(1);
    let mut mostradas = 0;
    let mut coincidencias = 0;
    for i in 0..=nu {
        for j in 0..=nv {
            for k in 0..=nw {
                let mut mejor = 0;
                let mut candidatas: Vec<((usize, usize, usize), Option<u32>)> = Vec::new();
                for &movimiento in MOVIMIENTOS_DE_TRES {
                    let llegada =
                        llegada_de_tres(&cubo, (ancho, fondo), secuencias, (i, j, k), movimiento);
                    if let Some(valor) = llegada {
                        mejor = mejor.max(valor);
                    }
                    candidatas.push((movimiento, llegada));
                }
                cubo[(i * ancho + j) * fondo + k] = mejor;
                let tres_iguales =
                    i > 0 && j > 0 && k > 0 && u[i - 1] == v[j - 1] && v[j - 1] == w[k - 1];
                let es_detalle = i > 0 && j > 0 && k > 0 && mostradas < 3;
                let es_evento = tres_iguales && coincidencias < 6;
                if !(es_detalle || es_evento) {
                    continue;
                }
                mostradas += 1;
                if tres_iguales {
                    coincidencias += 1;
                }
                let elementos: Vec<String> = candidatas
                    .iter()
                    .map(|&(movimiento, valor)| {
                        format!(
                            "{}: {}",
                            nombre_de_movimiento(movimiento),
                            valor.map_or("fuera del cubo".to_string(), |v| v.to_string())
                        )
                    })
                    .collect();
                let ganador = candidatas
                    .iter()
                    .position(|&(_, valor)| valor == Some(mejor))
                    .unwrap_or(0);
                let nota = if tres_iguales {
                    format!(
                        "Celda ({i}, {j}, {k}): las tres letras son «{}». Avanzar las tres a la vez suma una columna con las tres iguales: la celda anterior valía {}, ahora {mejor}. Es el único movimiento que suma; los demás solo heredan.",
                        u[i - 1] as char,
                        mejor.saturating_sub(1)
                    )
                } else {
                    format!(
                        "Celda ({i}, {j}, {k}): letras «{}», «{}» y «{}», no son las tres iguales. Se prueban los siete movimientos posibles desde las celdas vecinas y ninguno suma, solo heredan: se queda el mayor, {mejor}.",
                        u[i - 1] as char,
                        v[j - 1] as char,
                        w[k - 1] as char
                    )
                };
                narracion.agregar(
                    dibujar(
                        base()
                            .variable("i", i)
                            .variable("j", j)
                            .variable("k", k)
                            .variable("mejor", mejor)
                            .lista(
                                "los siete movimientos (de dónde llega)",
                                elementos,
                                &[ganador],
                            )
                            .secuencia(u, &[(i - 1, i, Resalte::Actual)]),
                        &cubo,
                        i,
                        (j, k),
                        (j, k),
                        &[(j - 1, k - 1), (j - 1, k), (j, k - 1)],
                    )
                    .paso(
                        if tres_iguales {
                            "Some(if tres_iguales { antes + 1 } else { antes })"
                        } else {
                            "mejor = mejor.max(llegada);"
                        },
                        nota,
                    ),
                );
            }
        }
        if i.is_multiple_of(cada_capa) && i > 0 && i < nu {
            narracion.agregar(
                dibujar(
                    base()
                        .variable("i", i)
                        .secuencia(u, &[(i - 1, i, Resalte::Actual)]),
                    &cubo,
                    i,
                    (nv, nw),
                    (nv, nw),
                    &[],
                )
                .paso(
                    "for i in 0..=u.len() {",
                    format!(
                        "Terminó la capa {i} de {nu}. Usando las primeras {i} letras de u con todo v y todo w, caben {} columnas con las tres letras iguales.",
                        cubo[(i * ancho + nv) * fondo + nw]
                    ),
                ),
            );
        }
    }
    let puntaje = cubo[(nu * ancho + nv) * fondo + nw] as usize;
    narracion.agregar(
        dibujar(
            base().variable("puntaje", puntaje),
            &cubo,
            nu,
            (nv, nw),
            (nv, nw),
            &[],
        )
        .paso(
            "let (mut i, mut j, mut k) = (u.len(), v.len(), w.len());",
            format!(
                "El cubo está lleno: la esquina final dice que se pueden alinear {puntaje} columnas con las tres letras iguales. Para saber cuáles, se vuelve desde esa esquina hacia atrás, buscando en cada celda el movimiento que dio su valor."
            ),
        ),
    );
    let mut fila_u = Vec::new();
    let mut fila_v = Vec::new();
    let mut fila_w = Vec::new();
    let (mut i, mut j, mut k) = (nu, nv, nw);
    let mut movimientos = Vec::new();
    while i > 0 || j > 0 || k > 0 {
        let actual = cubo[(i * ancho + j) * fondo + k];
        let Some(&(paso_u, paso_v, paso_w)) = MOVIMIENTOS_DE_TRES.iter().find(|&&movimiento| {
            llegada_de_tres(&cubo, (ancho, fondo), secuencias, (i, j, k), movimiento)
                == Some(actual)
        }) else {
            break;
        };
        movimientos.push((i, j, k, (paso_u, paso_v, paso_w)));
        fila_u.push(if paso_u == 1 { u[i - 1] } else { b'-' });
        fila_v.push(if paso_v == 1 { v[j - 1] } else { b'-' });
        fila_w.push(if paso_w == 1 { w[k - 1] } else { b'-' });
        i -= paso_u;
        j -= paso_v;
        k -= paso_w;
    }
    let ritmo = Ritmo::nuevo(movimientos.len(), 6, 14);
    for (indice, &(desde_i, desde_j, desde_k, movimiento)) in movimientos.iter().enumerate() {
        if !ritmo.toca(indice) {
            continue;
        }
        let columna = (
            if movimiento.0 == 1 {
                u[desde_i - 1] as char
            } else {
                '-'
            },
            if movimiento.1 == 1 {
                v[desde_j - 1] as char
            } else {
                '-'
            },
            if movimiento.2 == 1 {
                w[desde_k - 1] as char
            } else {
                '-'
            },
        );
        narracion.agregar(
            dibujar(
                base()
                    .variable("i", desde_i)
                    .variable("j", desde_j)
                    .variable("k", desde_k)
                    .lista(
                        "columna que se agrega",
                        vec![columna.0.to_string(), columna.1.to_string(), columna.2.to_string()],
                        &[],
                    )
                    .secuencia(u, &[(desde_i.saturating_sub(1), desde_i, Resalte::Actual)]),
                &cubo,
                desde_i,
                (nv, nw),
                (desde_j, desde_k),
                &[],
            )
            .paso(
                "fila_u.push(if paso_u == 1 { u[i - 1] } else { b'-' });",
                format!(
                    "En la celda ({desde_i}, {desde_j}, {desde_k}), el movimiento {} da el mismo valor que la celda: se agrega la columna «{}» / «{}» / «{}» (un guion es un hueco) y se retrocede.",
                    nombre_de_movimiento(movimiento),
                    columna.0,
                    columna.1,
                    columna.2
                ),
            ),
        );
    }
    fila_u.reverse();
    fila_v.reverse();
    fila_w.reverse();
    let ultimo = base()
        .variable("puntaje", puntaje)
        .lista(
            "alineamiento de las tres",
            vec![corto(&fila_u), corto(&fila_v), corto(&fila_w)],
            &[],
        )
        .paso(
            "fila_u.reverse();",
            format!(
                "Las columnas se juntaron de atrás hacia adelante; al darlas vuelta queda el alineamiento de las tres, con {} columnas y {puntaje} de ellas con las tres letras iguales. Mira «Qué significa».",
                miles(fila_u.len())
            ),
        );
    let triple = alineamiento::AlineamientoTriple {
        puntaje,
        fila_u,
        fila_v,
        fila_w,
    };
    (narracion.terminar(ultimo), triple)
}

#[cfg(test)]
mod pruebas {
    use super::*;

    #[test]
    fn la_probabilidad_de_azar_no_se_desborda_con_solapamientos_largos() {
        let probabilidad = probabilidad_de_azar(1400, 660);
        assert!(probabilidad.is_finite());
        assert!(probabilidad < 1e-10);
        assert!((probabilidad_de_azar(4, 4) - 1.0).abs() < 1e-12);
        assert!((probabilidad_de_azar(1, 0) - 0.25).abs() < 1e-12);
    }

    #[test]
    fn los_pesos_enormes_se_rechazan() {
        let entrada = Entrada::nueva(&[
            ("abajo", "1e308 1e308 1e308\n1e308 1e308 1e308"),
            ("derecha", "1 1\n1 1\n1 1"),
        ]);
        assert!(manhattan(&entrada).is_err());
        let grafo = Entrada::nueva(&[
            ("aristas", "a -> b: 1e308\nb -> c: 1e308"),
            ("fuente", "a"),
            ("sumidero", "c"),
        ]);
        assert!(camino_dag(&grafo).is_err());
    }

    #[test]
    fn la_extension_mayor_que_la_apertura_se_rechaza() {
        let entrada = Entrada::nueva(&[
            ("v", "ACGTACGT"),
            ("w", "ACGACGT"),
            ("puntuacion", "ADN +1/-1"),
            ("apertura", "2"),
            ("extension", "3"),
        ]);
        assert!(afin(&entrada).is_err());
    }

    fn ejemplo_de(id: &str) -> Entrada {
        Entrada::desde_json(crate::ejemplo(id).unwrap()).unwrap()
    }

    #[test]
    fn el_cambio_narrado_da_la_misma_respuesta() {
        let (dinero, monedas) = leer_cambio(&ejemplo_de("alineamiento.cambio")).unwrap();
        let (pasos, usadas) = narrar_cambio(dinero, &monedas);
        assert_eq!(usadas, alineamiento::cambio_minimo(dinero, &monedas));
        assert!(pasos.len() > 10);
        let imposible = narrar_cambio(7, &[5, 10]);
        assert_eq!(imposible.1, None);
        let cero = narrar_cambio(0, &[1, 5]);
        assert_eq!(cero.1, Some(Vec::new()));
    }

    #[test]
    fn manhattan_narrado_da_el_mismo_recorrido() {
        let (abajo, derecha) = leer_manhattan(&ejemplo_de("alineamiento.manhattan")).unwrap();
        let (pasos, longitud, camino) = narrar_manhattan(&abajo, &derecha);
        let real = alineamiento::turista_de_manhattan(&abajo, &derecha).unwrap();
        assert_eq!(longitud, real.longitud);
        assert_eq!(camino, real.camino);
        assert!(pasos.len() > 10);
    }

    #[test]
    fn el_orden_topologico_narrado_da_el_mismo_orden() {
        let entrada = ejemplo_de("alineamiento.orden_topologico");
        let grafo = leer_grafo(&entrada, false).unwrap();
        let pares: Vec<(usize, usize)> = grafo
            .aristas
            .iter()
            .map(|&(desde, hasta, _)| (desde, hasta))
            .collect();
        let (pasos, orden) = narrar_orden_topologico(&grafo.nombres, &pares);
        assert_eq!(
            orden,
            alineamiento::orden_topologico(grafo.nombres.len(), &pares).unwrap()
        );
        assert!(pasos.len() > 10);
    }

    #[test]
    fn el_camino_en_el_grafo_narrado_da_el_mismo_camino() {
        let entrada = ejemplo_de("alineamiento.camino_dag");
        let grafo = leer_grafo(&entrada, true).unwrap();
        let fuente = buscar_nodo(&grafo, &entrada, "fuente").unwrap();
        let sumidero = buscar_nodo(&grafo, &entrada, "sumidero").unwrap();
        let (pasos, longitud, camino) = narrar_camino_dag(&grafo, fuente, sumidero).unwrap();
        let real =
            alineamiento::camino_mas_largo(grafo.nombres.len(), &grafo.aristas, fuente, sumidero)
                .unwrap();
        assert_eq!(longitud, real.longitud);
        assert_eq!(camino, real.nodos);
        assert!(pasos.len() > 10);
    }

    #[test]
    fn la_subsecuencia_narrada_da_la_misma_respuesta() {
        let entrada = ejemplo_de("alineamiento.subsecuencia_comun");
        let (v, _) = leer_secuencia(&entrada, "v").unwrap();
        let (w, _) = leer_secuencia(&entrada, "w").unwrap();
        let (pasos, comun) = narrar_subsecuencia(&v, &w);
        assert_eq!(comun, alineamiento::subsecuencia_comun_mas_larga(&v, &w));
        assert!(pasos.len() > 10);
    }

    #[test]
    fn la_edicion_narrada_da_la_misma_distancia() {
        let entrada = ejemplo_de("alineamiento.edicion");
        let (v, _) = leer_secuencia(&entrada, "v").unwrap();
        let (w, _) = leer_secuencia(&entrada, "w").unwrap();
        let (pasos, distancia) = narrar_edicion(&v, &w);
        assert_eq!(distancia, alineamiento::distancia_de_edicion(&v, &w));
        assert!(pasos.len() > 10);
    }

    fn comprobar_tabla(id: &str, variante: Variante) {
        let entrada = ejemplo_de(id);
        let datos = leer_tabla(&entrada).unwrap();
        let (pasos, alineado) = narrar_tabla(variante, &datos);
        let real = match variante {
            Variante::Local => {
                alineamiento::alineamiento_local(&datos.v, &datos.w, datos.puntuar, datos.indel)
            }
            Variante::Ajustado => {
                alineamiento::alineamiento_ajustado(&datos.v, &datos.w, datos.puntuar, datos.indel)
            }
            Variante::Solapado => {
                alineamiento::alineamiento_solapado(&datos.v, &datos.w, datos.puntuar, datos.indel)
            }
        };
        assert_eq!(alineado, real);
        assert!(pasos.len() > 10);
    }

    #[test]
    fn las_tablas_narradas_dan_el_mismo_alineamiento() {
        comprobar_tabla("alineamiento.local", Variante::Local);
        comprobar_tabla("alineamiento.ajustado", Variante::Ajustado);
        comprobar_tabla("alineamiento.solapado", Variante::Solapado);
    }

    #[test]
    fn el_alineamiento_afin_narrado_da_el_mismo_alineamiento() {
        let datos = leer_afin(&ejemplo_de("alineamiento.afin")).unwrap();
        let (pasos, alineado) = narrar_afin(&datos);
        let real = alineamiento::alineamiento_afin(
            &datos.v,
            &datos.w,
            datos.puntuar,
            datos.apertura,
            datos.extension,
        );
        assert_eq!(alineado, real);
        assert!(pasos.len() > 10);
    }

    #[test]
    fn la_arista_central_narrada_da_la_misma_arista() {
        let entrada = ejemplo_de("alineamiento.arista_central");
        let (v, _) = leer_secuencia(&entrada, "v").unwrap();
        let (w, _) = leer_secuencia(&entrada, "w").unwrap();
        let (puntuar, _) = leer_puntuacion(&entrada).unwrap();
        let indel = leer_castigo(&entrada, "indel").unwrap();
        let (pasos, arista) = narrar_arista_central(&v, &w, puntuar, indel).unwrap();
        assert_eq!(
            Some(arista),
            alineamiento::arista_central(&v, &w, puntuar, indel)
        );
        assert!(pasos.len() > 10);
    }

    #[test]
    fn el_espacio_lineal_narrado_da_el_mismo_alineamiento() {
        let entrada = ejemplo_de("alineamiento.espacio_lineal");
        let (v, _) = leer_secuencia(&entrada, "v").unwrap();
        let (w, _) = leer_secuencia(&entrada, "w").unwrap();
        let (puntuar, _) = leer_puntuacion(&entrada).unwrap();
        let indel = leer_castigo(&entrada, "indel").unwrap();
        let (pasos, alineado) = narrar_espacio_lineal(&v, &w, puntuar, indel);
        assert_eq!(
            alineado,
            alineamiento::alineamiento_espacio_lineal(&v, &w, puntuar, indel)
        );
        assert!(pasos.len() > 10);
    }

    #[test]
    fn el_cubo_narrado_da_el_mismo_alineamiento() {
        let entrada = ejemplo_de("alineamiento.tres");
        let (u, _) = leer_secuencia(&entrada, "u").unwrap();
        let (v, _) = leer_secuencia(&entrada, "v").unwrap();
        let (w, _) = leer_secuencia(&entrada, "w").unwrap();
        let (pasos, triple) = narrar_tres(&u, &v, &w);
        assert_eq!(triple, alineamiento::subsecuencia_comun_de_tres(&u, &v, &w));
        assert!(pasos.len() > 10);
    }

    #[test]
    fn cada_paso_narrado_apunta_a_una_linea_del_codigo_mostrado() {
        for algoritmo in algoritmos() {
            let Some(trazador) = algoritmo.pasos else {
                continue;
            };
            let entrada = ejemplo_de(algoritmo.id);
            let pasos = trazador(&entrada).unwrap();
            let fragmentos = crate::codigo(algoritmo.id).unwrap();
            let texto: Vec<&str> = fragmentos
                .iter()
                .flat_map(|fragmento| fragmento.fuente.lines())
                .collect();
            for paso in &pasos {
                let coincidencias = texto
                    .iter()
                    .filter(|linea| linea.contains(paso.linea.as_str()))
                    .count();
                assert!(
                    coincidencias >= 1,
                    "{}: la línea «{}» no está en el código",
                    algoritmo.id,
                    paso.linea
                );
            }
            assert!(pasos.len() <= MAXIMO_PASOS);
        }
    }
}
