use crate::formato::{miles, porcentaje};
use crate::pizarra::{Marca, Narracion, Pizarra};
use crate::{Algoritmo, Arista, Bloque, Campo, Entrada, Paso, Punto, Salida, Tipo};
use nucleo::reordenamientos::{
    self, BloqueSintenico, Casilla, DosRoturas, Genoma, ParCompartido, VentanaDeKmeros,
};
use serde_json::{Value, json};
use std::collections::{HashMap, HashSet};

const MAXIMO_BLOQUES: usize = 100_000;
const MAXIMO_LISTADO: usize = 2_000;
const MAXIMO_ANIMACION: usize = 150;
const MAXIMO_PASOS: usize = 400;
const MAXIMO_NODOS_DIBUJO: usize = 60;
const MAXIMO_CARACTERES: usize = 400_000;
const MAXIMO_FILAS: usize = 2_000;
const MAXIMO_PUNTOS: usize = 2_000;
const MAXIMO_A: usize = 5_000_000;
const MAXIMO_B: usize = 2_000_000;
const MAXIMO_EN_FRASE: usize = 20;

const CAMPO_PERMUTACION: &[Campo] = &[Campo {
    id: "permutacion",
    etiqueta: "Permutación con signo, por ejemplo (+3 -1 +2)",
    tipo: Tipo::Texto,
    ejemplo: "(+5 -1 +2 +3 -6 +4)",
}];

const CAMPO_RUPTURAS: &[Campo] = &[Campo {
    id: "permutacion",
    etiqueta: "Permutación con signo, por ejemplo (+3 -1 +2)",
    tipo: Tipo::Texto,
    ejemplo: "(+4 +5 -3 -2 +1 +6 +7)",
}];

const CAMPO_CROMOSOMA: &[Campo] = &[Campo {
    id: "cromosoma",
    etiqueta: "Cromosoma circular con bloques del 1 al n",
    tipo: Tipo::Texto,
    ejemplo: "(+2 -4 -1 +3)",
}];

const CAMPO_CICLO: &[Campo] = &[Campo {
    id: "ciclo",
    etiqueta: "Nodos del ciclo, dos por bloque",
    tipo: Tipo::Texto,
    ejemplo: "(4 3 5 6 2 1 8 7)",
}];

const CAMPO_GENOMA: &[Campo] = &[Campo {
    id: "genoma",
    etiqueta: "Genoma: cromosomas circulares entre paréntesis",
    tipo: Tipo::Texto,
    ejemplo: "(+1 -3 +2)(-4 +5)",
}];

const CAMPO_ARISTAS: &[Campo] = &[Campo {
    id: "aristas",
    etiqueta: "Aristas de color, por ejemplo (2, 3), (4, 1)",
    tipo: Tipo::Texto,
    ejemplo: "(2, 4), (3, 6), (5, 1), (8, 10), (9, 7)",
}];

const CAMPOS_DISTANCIA: &[Campo] = &[
    Campo {
        id: "p",
        etiqueta: "Genoma P",
        tipo: Tipo::Texto,
        ejemplo: "(+1 -2 +3 +4 -5 +6 +7 +8)",
    },
    Campo {
        id: "q",
        etiqueta: "Genoma Q",
        tipo: Tipo::Texto,
        ejemplo: "(+1 +2 +3 +4)(-7 -5 +6 +8)",
    },
];

const CAMPOS_DOS_ROTURAS_GRAFO: &[Campo] = &[
    Campo {
        id: "aristas",
        etiqueta: "Aristas de color del genoma",
        tipo: Tipo::Texto,
        ejemplo: "(2, 4), (3, 8), (7, 5), (6, 1)",
    },
    Campo {
        id: "nodos",
        etiqueta: "Los cuatro nodos i1, i2, i3, i4",
        tipo: Tipo::Texto,
        ejemplo: "1, 6, 3, 8",
    },
];

const CAMPOS_DOS_ROTURAS_GENOMA: &[Campo] = &[
    Campo {
        id: "genoma",
        etiqueta: "Genoma",
        tipo: Tipo::Texto,
        ejemplo: "(+1 +2 -4 +3 +5 -6)",
    },
    Campo {
        id: "nodos",
        etiqueta: "Los cuatro nodos i1, i2, i3, i4",
        tipo: Tipo::Texto,
        ejemplo: "4, 8, 6, 9",
    },
];

const CAMPOS_TRANSFORMACION: &[Campo] = &[
    Campo {
        id: "p",
        etiqueta: "Genoma P (el de partida)",
        tipo: Tipo::Texto,
        ejemplo: "(+1 +2 +3 +4 +5 +6 +7 +8)",
    },
    Campo {
        id: "q",
        etiqueta: "Genoma Q (el de llegada)",
        tipo: Tipo::Texto,
        ejemplo: "(+1 -4 -3 +7)(+2 +5 -8)(+6)",
    },
];

const CAMPOS_KMEROS: &[Campo] = &[
    Campo {
        id: "a",
        etiqueta: "Secuencia A",
        tipo: Tipo::Secuencia,
        ejemplo: "@ecoli:4166000-4172500",
    },
    Campo {
        id: "b",
        etiqueta: "Secuencia B",
        tipo: Tipo::Secuencia,
        ejemplo: "@ecoli:4034500-4041000",
    },
    Campo {
        id: "k",
        etiqueta: "Longitud k de los k-meros (1 a 32)",
        tipo: Tipo::Entero,
        ejemplo: "15",
    },
];

const CAMPOS_SINTENIA: &[Campo] = &[
    Campo {
        id: "a",
        etiqueta: "Secuencia A",
        tipo: Tipo::Secuencia,
        ejemplo: "@ecoli:223000-230000",
    },
    Campo {
        id: "b",
        etiqueta: "Secuencia B",
        tipo: Tipo::Secuencia,
        ejemplo: "@ecoli:2725000-2732000",
    },
    Campo {
        id: "k",
        etiqueta: "Longitud k de los k-meros (1 a 32)",
        tipo: Tipo::Entero,
        ejemplo: "15",
    },
    Campo {
        id: "distancia",
        etiqueta: "Distancia para unir k-meros vecinos (pb)",
        tipo: Tipo::Entero,
        ejemplo: "100",
    },
    Campo {
        id: "minimo",
        etiqueta: "Mínimo de k-meros por bloque",
        tipo: Tipo::Entero,
        ejemplo: "20",
    },
];

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![
        Algoritmo {
            id: "reordenamientos.ordenamiento_voraz",
            capitulo: "reordenamientos",
            nombre: "Ordenamiento voraz por inversiones",
            resumen: "Ordena una permutación con signo colocando el 1, luego el 2, y así, con una inversión por vez.",
            campos: CAMPO_PERMUTACION,
            funcion: "reordenamientos::ordenamiento_voraz",
            ejecutar: ordenamiento_voraz,
            pasos: Some(pasos_ordenamiento_voraz),
        },
        Algoritmo {
            id: "reordenamientos.rupturas",
            capitulo: "reordenamientos",
            nombre: "Número de rupturas",
            resumen: "Cuenta los pares de vecinos que no son consecutivos y da una cota inferior de inversiones.",
            campos: CAMPO_RUPTURAS,
            funcion: "reordenamientos::numero_de_rupturas",
            ejecutar: numero_de_rupturas,
            pasos: Some(pasos_numero_de_rupturas),
        },
        Algoritmo {
            id: "reordenamientos.cromosoma_a_ciclo",
            capitulo: "reordenamientos",
            nombre: "De cromosoma a ciclo",
            resumen: "Convierte cada bloque con signo en sus dos extremos, cola y cabeza, para formar un ciclo de nodos.",
            campos: CAMPO_CROMOSOMA,
            funcion: "reordenamientos::cromosoma_a_ciclo",
            ejecutar: cromosoma_a_ciclo,
            pasos: Some(pasos_cromosoma_a_ciclo),
        },
        Algoritmo {
            id: "reordenamientos.ciclo_a_cromosoma",
            capitulo: "reordenamientos",
            nombre: "De ciclo a cromosoma",
            resumen: "Lee un ciclo de nodos de dos en dos y recupera el cromosoma con sus signos.",
            campos: CAMPO_CICLO,
            funcion: "reordenamientos::ciclo_a_cromosoma",
            ejecutar: ciclo_a_cromosoma,
            pasos: Some(pasos_ciclo_a_cromosoma),
        },
        Algoritmo {
            id: "reordenamientos.aristas_de_color",
            capitulo: "reordenamientos",
            nombre: "Aristas de color",
            resumen: "Las uniones entre bloques vecinos de un genoma, escritas como pares de nodos.",
            campos: CAMPO_GENOMA,
            funcion: "reordenamientos::aristas_de_color",
            ejecutar: aristas_de_color,
            pasos: Some(pasos_aristas_de_color),
        },
        Algoritmo {
            id: "reordenamientos.grafo_a_genoma",
            capitulo: "reordenamientos",
            nombre: "De grafo a genoma",
            resumen: "Recorre las aristas de color alternando con los bloques y reconstruye los cromosomas.",
            campos: CAMPO_ARISTAS,
            funcion: "reordenamientos::grafo_a_genoma",
            ejecutar: grafo_a_genoma,
            pasos: Some(pasos_grafo_a_genoma),
        },
        Algoritmo {
            id: "reordenamientos.distancia_dos_roturas",
            capitulo: "reordenamientos",
            nombre: "Distancia de 2-roturas",
            resumen: "El menor número de cortes y empalmes que transforma un genoma en otro: bloques menos ciclos.",
            campos: CAMPOS_DISTANCIA,
            funcion: "reordenamientos::distancia_dos_roturas",
            ejecutar: distancia_dos_roturas,
            pasos: Some(pasos_distancia_dos_roturas),
        },
        Algoritmo {
            id: "reordenamientos.dos_roturas_en_grafo",
            capitulo: "reordenamientos",
            nombre: "2-rotura en el grafo del genoma",
            resumen: "Quita dos aristas de color y une sus extremos de otra manera.",
            campos: CAMPOS_DOS_ROTURAS_GRAFO,
            funcion: "reordenamientos::dos_roturas_en_grafo",
            ejecutar: dos_roturas_en_grafo,
            pasos: Some(pasos_dos_roturas_en_grafo),
        },
        Algoritmo {
            id: "reordenamientos.dos_roturas_en_genoma",
            capitulo: "reordenamientos",
            nombre: "2-rotura en el genoma",
            resumen: "Aplica una 2-rotura a un genoma y muestra si fue inversión, fisión o fusión.",
            campos: CAMPOS_DOS_ROTURAS_GENOMA,
            funcion: "reordenamientos::dos_roturas_en_genoma",
            ejecutar: dos_roturas_en_genoma,
            pasos: Some(pasos_dos_roturas_en_genoma),
        },
        Algoritmo {
            id: "reordenamientos.transformacion_mas_corta",
            capitulo: "reordenamientos",
            nombre: "Transformación más corta con 2-roturas",
            resumen: "Una ruta de 2-roturas de P a Q con la menor cantidad posible de pasos.",
            campos: CAMPOS_TRANSFORMACION,
            funcion: "reordenamientos::transformacion_mas_corta",
            ejecutar: transformacion_mas_corta,
            pasos: Some(pasos_transformacion),
        },
        Algoritmo {
            id: "reordenamientos.kmeros_compartidos",
            capitulo: "reordenamientos",
            nombre: "k-meros compartidos",
            resumen: "Todas las posiciones donde A y B tienen el mismo k-mero, en la misma hebra o en la contraria.",
            campos: CAMPOS_KMEROS,
            funcion: "reordenamientos::kmeros_compartidos",
            ejecutar: kmeros_compartidos,
            pasos: Some(pasos_kmeros_compartidos),
        },
        Algoritmo {
            id: "reordenamientos.bloques_de_sintenia",
            capitulo: "reordenamientos",
            nombre: "Bloques de sintenia",
            resumen: "Agrupa los k-meros compartidos cercanos en bloques y los ordena como una permutación con signo.",
            campos: CAMPOS_SINTENIA,
            funcion: "reordenamientos::bloques_de_sintenia",
            ejecutar: bloques_de_sintenia,
            pasos: Some(pasos_bloques_de_sintenia),
        },
    ]
}

fn separar_piezas(texto: &str) -> Vec<String> {
    texto
        .replace('(', " ( ")
        .replace(')', " ) ")
        .replace(['−', '–'], "-")
        .replace([',', ';'], " ")
        .split_whitespace()
        .map(str::to_string)
        .collect()
}

fn leer_genoma(entrada: &Entrada, campo: &str, limite: usize) -> Result<Vec<Vec<i64>>, String> {
    let texto = entrada.texto(campo)?;
    let piezas = separar_piezas(&texto);
    let con_parentesis = piezas.iter().any(|pieza| pieza == "(" || pieza == ")");
    let mut genoma: Vec<Vec<i64>> = Vec::new();
    let mut abierto: Option<Vec<i64>> = None;
    let mut sueltos = Vec::new();
    let mut cantidad = 0;
    for pieza in &piezas {
        match pieza.as_str() {
            "(" => {
                if abierto.is_some() {
                    return Err(format!(
                        "En «{campo}» hay un «(» dentro de otro: cierra cada cromosoma antes de abrir el siguiente."
                    ));
                }
                abierto = Some(Vec::new());
            }
            ")" => {
                let cromosoma = abierto
                    .take()
                    .ok_or_else(|| format!("En «{campo}» sobra un «)»."))?;
                if cromosoma.is_empty() {
                    return Err(format!("En «{campo}» hay un cromosoma vacío «()»."));
                }
                genoma.push(cromosoma);
            }
            _ => {
                let bloque = pieza.parse::<i64>().map_err(|_| {
                    format!(
                        "En «{campo}», «{pieza}» no es un bloque: escribe enteros con signo, como +3 o -2."
                    )
                })?;
                cantidad += 1;
                if cantidad > limite {
                    return Err(format!(
                        "En «{campo}» hay más de {} bloques; aquí se admiten hasta {}.",
                        miles(limite),
                        miles(limite)
                    ));
                }
                match abierto.as_mut() {
                    Some(cromosoma) => cromosoma.push(bloque),
                    None if con_parentesis => {
                        return Err(format!(
                            "En «{campo}», «{pieza}» quedó fuera de los paréntesis."
                        ));
                    }
                    None => sueltos.push(bloque),
                }
            }
        }
    }
    if abierto.is_some() {
        return Err(format!("En «{campo}» falta cerrar un paréntesis."));
    }
    if !sueltos.is_empty() {
        genoma.push(sueltos);
    }
    reordenamientos::validar_genoma(&genoma)
        .map_err(|error| format!("En «{campo}», {}", minuscula_inicial(&error)))?;
    Ok(genoma)
}

fn minuscula_inicial(texto: &str) -> String {
    let mut letras = texto.chars();
    match letras.next() {
        Some(primera) => primera.to_lowercase().chain(letras).collect(),
        None => String::new(),
    }
}

fn leer_permutacion(entrada: &Entrada, campo: &str, limite: usize) -> Result<Vec<i64>, String> {
    let mut genoma = leer_genoma(entrada, campo, limite)?;
    if genoma.len() > 1 {
        return Err(format!(
            "En «{campo}» va una sola permutación, pero hay {} grupos entre paréntesis.",
            genoma.len()
        ));
    }
    Ok(genoma.remove(0))
}

fn leer_nodos(entrada: &Entrada, campo: &str) -> Result<Vec<usize>, String> {
    let texto = entrada.texto(campo)?;
    let mut nodos = Vec::new();
    for pieza in separar_piezas(&texto) {
        if pieza == "(" || pieza == ")" {
            continue;
        }
        let nodo = pieza.parse::<usize>().map_err(|_| {
            format!("En «{campo}», «{pieza}» no es un nodo: los nodos son enteros positivos.")
        })?;
        nodos.push(nodo);
        if nodos.len() > 2 * MAXIMO_BLOQUES {
            return Err(format!(
                "En «{campo}» hay más de {} nodos; aquí se admiten hasta {}.",
                miles(2 * MAXIMO_BLOQUES),
                miles(2 * MAXIMO_BLOQUES)
            ));
        }
    }
    if nodos.is_empty() {
        return Err(format!("El campo «{campo}» no tiene nodos."));
    }
    Ok(nodos)
}

fn leer_aristas(entrada: &Entrada, campo: &str) -> Result<Vec<(usize, usize)>, String> {
    let nodos = leer_nodos(entrada, campo)?;
    if nodos.len() % 2 == 1 {
        return Err(format!(
            "En «{campo}» las aristas van de a pares de nodos, pero hay {}.",
            cuantos(nodos.len(), "número", "números")
        ));
    }
    Ok(nodos.chunks(2).map(|par| (par[0], par[1])).collect())
}

fn leer_cuatro_nodos(entrada: &Entrada, campo: &str) -> Result<[usize; 4], String> {
    let nodos = leer_nodos(entrada, campo)?;
    if nodos.len() != 4 {
        return Err(format!(
            "En «{campo}» van exactamente cuatro nodos (i1, i2, i3, i4), pero hay {}.",
            nodos.len()
        ));
    }
    Ok([nodos[0], nodos[1], nodos[2], nodos[3]])
}

fn texto_permutacion(bloques: &[i64]) -> String {
    let piezas: Vec<String> = bloques.iter().map(|bloque| format!("{bloque:+}")).collect();
    format!("({})", piezas.join(" "))
}

fn texto_genoma(genoma: &[Vec<i64>]) -> String {
    genoma
        .iter()
        .map(|cromosoma| texto_permutacion(cromosoma))
        .collect()
}

fn texto_nodos(nodos: &[usize]) -> String {
    let piezas: Vec<String> = nodos.iter().map(|nodo| nodo.to_string()).collect();
    format!("({})", piezas.join(" "))
}

fn texto_aristas(aristas: &[(usize, usize)]) -> String {
    let piezas: Vec<String> = aristas.iter().map(|(a, b)| format!("({a}, {b})")).collect();
    piezas.join(", ")
}

fn cuantos(cantidad: usize, singular: &str, plural: &str) -> String {
    if cantidad == 1 {
        format!("1 {singular}")
    } else {
        format!("{} {plural}", miles(cantidad))
    }
}

fn contar_bloques(genoma: &[Vec<i64>]) -> usize {
    genoma.iter().map(|cromosoma| cromosoma.len()).sum()
}

fn lista_acotada(etiqueta: &str, elementos: impl Iterator<Item = String>, total: usize) -> Bloque {
    let mut visibles = Vec::new();
    let mut caracteres = 0;
    for elemento in elementos {
        caracteres += elemento.len();
        if caracteres > MAXIMO_CARACTERES || visibles.len() == MAXIMO_FILAS {
            break;
        }
        visibles.push(elemento);
    }
    Bloque::Lista {
        etiqueta: etiqueta.to_string(),
        elementos: visibles,
        total,
    }
}

fn muestrear<T: Clone>(elementos: &[T], maximo: usize) -> Vec<T> {
    if elementos.len() <= maximo || maximo < 2 {
        return elementos.to_vec();
    }
    let ultimo = elementos.len() - 1;
    (0..maximo)
        .map(|indice| elementos[indice * ultimo / (maximo - 1)].clone())
        .collect()
}

fn aristas_de_bloques(bloques: usize) -> Vec<Arista> {
    (1..=bloques)
        .map(|bloque| Arista {
            desde: 2 * bloque - 2,
            hasta: 2 * bloque - 1,
            etiqueta: format!("bloque {bloque}"),
        })
        .collect()
}

fn aristas_dibujables(aristas: &[(usize, usize)], etiqueta: &str) -> Vec<Arista> {
    aristas
        .iter()
        .map(|&(a, b)| Arista {
            desde: a - 1,
            hasta: b - 1,
            etiqueta: etiqueta.to_string(),
        })
        .collect()
}

fn nombres_de_nodos(bloques: usize) -> Vec<String> {
    (1..=2 * bloques).map(|nodo| nodo.to_string()).collect()
}

fn agregar_grafo_del_genoma(salida: &mut Salida, etiqueta: &str, aristas: &[(usize, usize)]) {
    let bloques = aristas.len();
    if 2 * bloques > MAXIMO_NODOS_DIBUJO {
        return;
    }
    let mut todas = aristas_de_bloques(bloques);
    todas.extend(aristas_dibujables(aristas, "unión"));
    salida.agregar(Bloque::grafo(etiqueta, nombres_de_nodos(bloques), todas));
}

fn efecto_de_la_rotura(antes: usize, despues: usize) -> &'static str {
    if despues > antes {
        "fisión"
    } else if despues < antes {
        "fusión"
    } else {
        "inversión"
    }
}

fn explicar_efecto(salida: &mut Salida, antes: usize, despues: usize) {
    match efecto_de_la_rotura(antes, despues) {
        "fisión" => salida.explicar(format!(
            "Fue una fisión: los dos cortes cayeron en el mismo cromosoma y el genoma pasó de {antes} a {}.",
            cuantos(despues, "cromosoma", "cromosomas")
        )),
        "fusión" => salida.explicar(format!(
            "Fue una fusión: los cortes cayeron en dos cromosomas distintos y se unieron en uno; el genoma pasó de {antes} a {}.",
            cuantos(despues, "cromosoma", "cromosomas")
        )),
        _ => salida.explicar(format!(
            "Fue una inversión: los dos cortes cayeron en el mismo cromosoma y el tramo entre ellos quedó dado vuelta; el genoma sigue con {}.",
            cuantos(antes, "cromosoma", "cromosomas")
        )),
    }
}

fn ordenamiento_voraz(entrada: &Entrada) -> Result<Salida, String> {
    let permutacion = leer_permutacion(entrada, "permutacion", MAXIMO_LISTADO)?;
    let permutaciones = reordenamientos::ordenamiento_voraz(&permutacion)?;
    let rupturas_iniciales = reordenamientos::numero_de_rupturas(&permutacion)?;
    let cota = rupturas_iniciales.div_ceil(2);
    let mut curva = vec![[0.0, rupturas_iniciales as f64]];
    let mut filas = Vec::new();
    let mut anterior = permutacion.as_slice();
    let mut volteos = 0;
    let mut subidas = 0;
    let mut rupturas_antes = rupturas_iniciales;
    for (numero, siguiente) in permutaciones.iter().enumerate() {
        let rupturas = reordenamientos::numero_de_rupturas(siguiente)?;
        curva.push([(numero + 1) as f64, rupturas as f64]);
        if rupturas > rupturas_antes {
            subidas += 1;
        }
        rupturas_antes = rupturas;
        let inicio = (0..anterior.len())
            .find(|&i| anterior[i] != siguiente[i])
            .unwrap_or(0);
        let fin = (0..anterior.len())
            .rfind(|&i| anterior[i] != siguiente[i])
            .unwrap_or(0);
        if inicio == fin {
            volteos += 1;
        }
        filas.push(vec![
            miles(numero + 1),
            format!("{}–{}", inicio + 1, fin + 1),
            miles(fin - inicio + 1),
            miles(rupturas),
        ]);
        anterior = siguiente;
    }
    let total = permutaciones.len();
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Inversiones usadas", miles(total)));
    salida.agregar(Bloque::valor(
        "Rupturas al inicio",
        miles(rupturas_iniciales),
    ));
    salida.agregar(Bloque::valor(
        "Mínimo posible: rupturas ÷ 2, hacia arriba",
        miles(cota),
    ));
    if total > 0 {
        salida.agregar(Bloque::linea(
            "Rupturas que quedan tras cada inversión",
            "Inversión",
            "Rupturas",
            &curva,
            Vec::new(),
        ));
        salida.agregar(Bloque::tabla(
            "Cada inversión (posiciones contando desde 1)",
            &[
                "Inversión",
                "Posiciones",
                "Bloques invertidos",
                "Rupturas después",
            ],
            filas,
        ));
        salida.agregar(lista_acotada(
            "Permutación tras cada inversión",
            permutaciones.iter().map(|una| texto_permutacion(una)),
            total,
        ));
    }
    if total == 0 {
        salida.explicar(
            "La permutación ya es la identidad (+1 +2 … +n): no hace falta ninguna inversión.",
        );
        return Ok(salida);
    }
    salida.explicar(format!(
        "El método voraz usó {}: pone el 1 en su lugar, luego el 2, y así, sin mirar el conjunto.",
        cuantos(total, "inversión", "inversiones")
    ));
    salida.explicar(format!(
        "Una inversión corta en dos puntos, así que elimina como mucho 2 rupturas: con {} hacen falta al menos {}. La distancia de inversiones real está entre {} y {}.",
        cuantos(rupturas_iniciales, "ruptura", "rupturas"),
        miles(cota),
        miles(cota),
        miles(total)
    ));
    if total == cota {
        salida.explicar("Aquí el voraz llegó justo a ese mínimo, así que su resultado es óptimo.");
    } else if total >= 2 * cota {
        salida.explicar(
            "El voraz gastó al menos el doble del mínimo: es fácil de seguir, pero no busca la ruta más corta.",
        );
    }
    if subidas > 0 {
        salida.explicar(format!(
            "En {} las rupturas aumentaron (mira la curva): el voraz no intenta reducirlas, por eso se aleja del mínimo.",
            cuantos(subidas, "paso", "pasos")
        ));
    }
    if volteos == 1 {
        salida.explicar(
            "Una de las inversiones solo cambió el signo de un bloque que ya estaba en su lugar pero al revés.",
        );
    } else if volteos > 1 {
        salida.explicar(format!(
            "{} de las inversiones solo cambiaron el signo de un bloque que ya estaba en su lugar pero al revés.",
            miles(volteos)
        ));
    }
    Ok(salida)
}

fn numero_de_rupturas(entrada: &Entrada) -> Result<Salida, String> {
    let permutacion = leer_permutacion(entrada, "permutacion", MAXIMO_BLOQUES)?;
    let rupturas = reordenamientos::rupturas(&permutacion)?;
    let cantidad = reordenamientos::numero_de_rupturas(&permutacion)?;
    let bloques = permutacion.len();
    let pares = bloques + 1;
    let adyacencias = pares - cantidad;
    let cota = cantidad.div_ceil(2);
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Rupturas", miles(cantidad)));
    salida.agregar(Bloque::valor("Adyacencias", miles(adyacencias)));
    salida.agregar(Bloque::valor(
        "Mínimo de inversiones: rupturas ÷ 2, hacia arriba",
        miles(cota),
    ));
    let filas: Vec<Vec<String>> = rupturas
        .iter()
        .take(MAXIMO_FILAS)
        .map(|&(izquierda, derecha)| {
            let izquierda_escrita = vecino_escrito(izquierda, pares);
            let derecha_escrita = vecino_escrito(derecha, pares);
            vec![
                format!("({izquierda_escrita} {derecha_escrita})"),
                format!(
                    "{derecha_escrita} − ({izquierda_escrita}) = {}",
                    derecha - izquierda
                ),
            ]
        })
        .collect();
    salida.agregar(Bloque::Tabla {
        etiqueta: "Rupturas: vecinos que no suman 1 al avanzar".to_string(),
        columnas: vec!["Par de vecinos".to_string(), "Diferencia".to_string()],
        filas,
        total: cantidad,
    });
    if bloques <= MAXIMO_PUNTOS {
        let puntos = permutacion
            .iter()
            .enumerate()
            .map(|(posicion, &bloque)| Punto {
                x: (posicion + 1) as f64,
                y: bloque.abs() as f64,
                grupo: if bloque > 0 { 0 } else { 1 },
                nombre: format!("{bloque:+}"),
            })
            .collect();
        salida.agregar(Bloque::dispersion(
            "Bloque en cada posición (grupo 0: signo +, grupo 1: signo −)",
            "Posición (contando desde 1)",
            "Número de bloque",
            puntos,
        ));
    }
    if cantidad == 0 {
        salida.explicar("No hay rupturas: es la permutación identidad, ya ordenada.");
        return Ok(salida);
    }
    salida.explicar(format!(
        "Con el 0 al principio y el {} al final hay {} pares de vecinos: {} y {}.",
        pares,
        miles(pares),
        cuantos(
            adyacencias,
            "adyacencia (el de la derecha vale 1 más)",
            "adyacencias (el de la derecha vale 1 más)"
        ),
        cuantos(cantidad, "ruptura", "rupturas")
    ));
    salida.explicar(format!(
        "Una inversión corta en dos puntos, así que elimina como mucho 2 rupturas: necesitas al menos {} para ordenar esta permutación.",
        cuantos(cota, "inversión", "inversiones")
    ));
    if cantidad == pares {
        salida.explicar(
            "Ningún par de vecinos está en orden: no queda ningún tramo intacto de la identidad.",
        );
    }
    let negativas: Vec<&[i64]> = permutacion
        .windows(2)
        .filter(|par| par[0] < 0 && par[1] - par[0] == 1)
        .collect();
    if let Some(ejemplo) = negativas.first() {
        let (a, b) = (ejemplo[0], ejemplo[1]);
        if negativas.len() == 1 {
            salida.explicar(format!(
                "La adyacencia ({a:+} {b:+}) es de bloques negativos: también cuenta, porque leída al revés es ({:+} {:+}).",
                -b, -a
            ));
        } else {
            salida.explicar(format!(
                "{} adyacencias son de bloques negativos, como ({a:+} {b:+}): también cuentan, porque leídas al revés son ({:+} {:+}).",
                miles(negativas.len()),
                -b,
                -a
            ));
        }
    }
    Ok(salida)
}

fn vecino_escrito(valor: i64, extremo: usize) -> String {
    if valor == 0 || valor == extremo as i64 {
        valor.to_string()
    } else {
        format!("{valor:+}")
    }
}

fn cromosoma_a_ciclo(entrada: &Entrada) -> Result<Salida, String> {
    let cromosoma = leer_permutacion(entrada, "cromosoma", MAXIMO_BLOQUES)?;
    let ciclo = reordenamientos::cromosoma_a_ciclo(&cromosoma)?;
    let bloques = cromosoma.len();
    let mut salida = Salida::default();
    salida.agregar(Bloque::texto("Ciclo de nodos", texto_nodos(&ciclo)));
    let filas: Vec<Vec<String>> = cromosoma
        .iter()
        .take(MAXIMO_FILAS)
        .enumerate()
        .map(|(posicion, &bloque)| {
            let numero = bloque.unsigned_abs() as usize;
            vec![
                format!("{bloque:+}"),
                miles(2 * numero - 1),
                miles(2 * numero),
                format!("{} → {}", ciclo[2 * posicion], ciclo[2 * posicion + 1]),
            ]
        })
        .collect();
    salida.agregar(Bloque::Tabla {
        etiqueta: "Cada bloque y sus dos extremos".to_string(),
        columnas: vec![
            "Bloque".to_string(),
            "Cola (2i−1)".to_string(),
            "Cabeza (2i)".to_string(),
            "Orden en el ciclo".to_string(),
        ],
        filas,
        total: bloques,
    });
    if 2 * bloques <= MAXIMO_NODOS_DIBUJO {
        let nodos: Vec<String> = ciclo.iter().map(|nodo| nodo.to_string()).collect();
        let mut aristas = Vec::new();
        for (posicion, bloque) in cromosoma.iter().enumerate() {
            aristas.push(Arista {
                desde: 2 * posicion,
                hasta: 2 * posicion + 1,
                etiqueta: format!("{bloque:+}"),
            });
            aristas.push(Arista {
                desde: 2 * posicion + 1,
                hasta: (2 * posicion + 2) % (2 * bloques),
                etiqueta: "unión".to_string(),
            });
        }
        salida.agregar(Bloque::grafo("El cromosoma como ciclo", nodos, aristas));
    }
    let negativos = cromosoma.iter().filter(|&&bloque| bloque < 0).count();
    salida.explicar(format!(
        "Cada bloque i se convierte en dos nodos: su cola 2i−1 y su cabeza 2i. Con {} salen {} nodos.",
        cuantos(bloques, "bloque", "bloques"),
        miles(2 * bloques)
    ));
    if negativos == 0 {
        salida.explicar(
            "Todos los bloques tienen signo +: en el ciclo cada cola va antes que su cabeza.",
        );
    } else {
        salida.explicar(format!(
            "{} de los {} {} signo −: el ciclo recorre esos bloques al revés, primero la cabeza y luego la cola.",
            miles(negativos),
            miles(bloques),
            if negativos == 1 {
                "bloques tiene"
            } else {
                "bloques tienen"
            }
        ));
    }
    salida.explicar(
        "Escribir los cromosomas así, con extremos de bloque, es lo que permite construir el grafo de rupturas y comparar dos genomas.",
    );
    Ok(salida)
}

fn ciclo_a_cromosoma(entrada: &Entrada) -> Result<Salida, String> {
    let nodos = leer_nodos(entrada, "ciclo")?;
    let cromosoma = reordenamientos::ciclo_a_cromosoma(&nodos)?;
    let mut salida = Salida::default();
    salida.agregar(Bloque::texto("Cromosoma", texto_permutacion(&cromosoma)));
    let filas: Vec<Vec<String>> = cromosoma
        .iter()
        .take(MAXIMO_FILAS)
        .enumerate()
        .map(|(posicion, &bloque)| {
            vec![
                format!("{} → {}", nodos[2 * posicion], nodos[2 * posicion + 1]),
                format!("{bloque:+}"),
                if bloque > 0 {
                    "de cola a cabeza".to_string()
                } else {
                    "de cabeza a cola".to_string()
                },
            ]
        })
        .collect();
    salida.agregar(Bloque::Tabla {
        etiqueta: "Cada par de nodos y su bloque".to_string(),
        columnas: vec![
            "Nodos".to_string(),
            "Bloque".to_string(),
            "Sentido".to_string(),
        ],
        filas,
        total: cromosoma.len(),
    });
    let negativos = cromosoma.iter().filter(|&&bloque| bloque < 0).count();
    salida.explicar(format!(
        "El ciclo tiene {} nodos, así que el cromosoma tiene {}: cada par (2i−1, 2i) es el bloque +i y cada par (2i, 2i−1) es el bloque −i.",
        miles(nodos.len()),
        cuantos(cromosoma.len(), "bloque", "bloques")
    ));
    if negativos == 1 {
        salida.explicar("Un par va de la cabeza a la cola, por eso ese bloque lleva signo −.");
    } else if negativos > 1 {
        salida.explicar(format!(
            "{} pares van de la cabeza a la cola, por eso esos bloques llevan signo −.",
            miles(negativos)
        ));
    } else {
        salida.explicar("Todos los pares van de cola a cabeza: todos los bloques llevan signo +.");
    }
    Ok(salida)
}

fn aristas_de_color(entrada: &Entrada) -> Result<Salida, String> {
    let genoma = leer_genoma(entrada, "genoma", MAXIMO_BLOQUES)?;
    let aristas = reordenamientos::aristas_de_color(&genoma)?;
    let mut salida = Salida::default();
    salida.agregar(Bloque::texto("Aristas de color", texto_aristas(&aristas)));
    salida.agregar(Bloque::valor("Aristas de color", miles(aristas.len())));
    salida.agregar(Bloque::valor("Cromosomas", miles(genoma.len())));
    agregar_grafo_del_genoma(&mut salida, "Grafo del genoma: bloques y uniones", &aristas);
    salida.explicar(format!(
        "Cada cromosoma circular de m bloques tiene m uniones entre vecinos, contando la que cierra el círculo: {} {} {}.",
        cuantos(aristas.len(), "bloque", "bloques"),
        if aristas.len() == 1 { "da" } else { "dan" },
        cuantos(aristas.len(), "arista de color", "aristas de color")
    ));
    salida.explicar(
        "Una arista (x, y) une la salida de un bloque con la entrada del siguiente: en el grafo del genoma se alternan aristas de bloque y aristas de color.",
    );
    if genoma.iter().any(|cromosoma| cromosoma.len() == 1) {
        salida.explicar(
            "Un cromosoma de un solo bloque da una arista de color que une la cabeza y la cola de ese mismo bloque.",
        );
    }
    if genoma.len() > 1 {
        salida.explicar(format!(
            "El genoma tiene {} cromosomas, así que el grafo tiene {} ciclos separados.",
            miles(genoma.len()),
            miles(genoma.len())
        ));
    }
    Ok(salida)
}

fn grafo_a_genoma(entrada: &Entrada) -> Result<Salida, String> {
    let aristas = leer_aristas(entrada, "aristas")?;
    let genoma = reordenamientos::grafo_a_genoma(&aristas)?;
    let mut salida = Salida::default();
    salida.agregar(Bloque::texto("Genoma", texto_genoma(&genoma)));
    salida.agregar(Bloque::valor("Cromosomas", miles(genoma.len())));
    salida.agregar(Bloque::valor("Bloques", miles(aristas.len())));
    agregar_grafo_del_genoma(&mut salida, "Grafo del genoma: bloques y uniones", &aristas);
    salida.explicar(format!(
        "Con {} hay {}. El recorrido alterna: cruza un bloque de un extremo al otro y salta por la arista de color al siguiente.",
        cuantos(aristas.len(), "arista de color", "aristas de color"),
        cuantos(aristas.len(), "bloque", "bloques")
    ));
    salida.explicar(format!(
        "Cada vez que el recorrido vuelve al punto de partida se cierra un cromosoma circular: {}.",
        if genoma.len() == 1 {
            "salió 1".to_string()
        } else {
            format!("salieron {}", miles(genoma.len()))
        }
    ));
    salida.explicar(
        "Cada cromosoma se escribe empezando por su bloque de menor número con signo +: un círculo se puede leer desde cualquier punto y en cualquier sentido, así que es el mismo cromosoma.",
    );
    Ok(salida)
}

fn distancia_dos_roturas(entrada: &Entrada) -> Result<Salida, String> {
    let p = leer_genoma(entrada, "p", MAXIMO_BLOQUES)?;
    let q = leer_genoma(entrada, "q", MAXIMO_BLOQUES)?;
    let ciclos = reordenamientos::ciclos_del_grafo_de_rupturas(&p, &q)?;
    let distancia = reordenamientos::distancia_dos_roturas(&p, &q)?;
    let bloques = contar_bloques(&p);
    let triviales = ciclos.iter().filter(|ciclo| ciclo.len() == 2).count();
    let mayor = ciclos.iter().map(|ciclo| ciclo.len()).max().unwrap_or(0);
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Distancia de 2-roturas", miles(distancia)));
    salida.agregar(Bloque::valor("Bloques", miles(bloques)));
    salida.agregar(Bloque::valor("Ciclos rojo-azul", miles(ciclos.len())));
    salida.agregar(Bloque::valor("Ciclos triviales", miles(triviales)));
    let mut ordenados: Vec<&Vec<usize>> = ciclos.iter().collect();
    ordenados.sort_by_key(|ciclo| std::cmp::Reverse(ciclo.len()));
    let filas: Vec<Vec<String>> = ordenados
        .iter()
        .take(MAXIMO_FILAS)
        .map(|ciclo| {
            vec![
                texto_nodos(ciclo),
                miles(ciclo.len()),
                if ciclo.len() == 2 { "sí" } else { "no" }.to_string(),
            ]
        })
        .collect();
    salida.agregar(Bloque::Tabla {
        etiqueta: "Ciclos del grafo de rupturas (rojo: P, azul: Q)".to_string(),
        columnas: vec![
            "Nodos del ciclo".to_string(),
            "Aristas".to_string(),
            "Trivial".to_string(),
        ],
        filas,
        total: ciclos.len(),
    });
    if 2 * bloques <= MAXIMO_NODOS_DIBUJO {
        let mut aristas = aristas_de_bloques(bloques);
        aristas.extend(aristas_dibujables(
            &reordenamientos::aristas_de_color(&p)?,
            "P",
        ));
        aristas.extend(aristas_dibujables(
            &reordenamientos::aristas_de_color(&q)?,
            "Q",
        ));
        salida.agregar(Bloque::grafo(
            "Grafo de rupturas: aristas de P y de Q",
            nombres_de_nodos(bloques),
            aristas,
        ));
    }
    salida.explicar(format!(
        "d(P, Q) = bloques − ciclos = {} − {} = {}: es el menor número de 2-roturas (cortar dos uniones y volver a pegar los extremos) que lleva de P a Q.",
        miles(bloques),
        miles(ciclos.len()),
        miles(distancia)
    ));
    if distancia == 0 {
        salida.explicar("P y Q son el mismo genoma: todos los ciclos son triviales.");
        return Ok(salida);
    }
    if triviales == 0 {
        salida.explicar(
            "Ningún ciclo es trivial: P y Q no comparten ni una sola unión entre bloques.",
        );
    } else if triviales == 1 {
        salida.explicar(format!(
            "1 de los {} ciclos es trivial (una arista roja igual a una azul): es una unión entre bloques que P y Q ya comparten.",
            miles(ciclos.len())
        ));
    } else {
        salida.explicar(format!(
            "{} de los {} ciclos son triviales (una arista roja igual a una azul): son uniones entre bloques que P y Q ya comparten.",
            miles(triviales),
            miles(ciclos.len())
        ));
    }
    if mayor > 2 {
        salida.explicar(format!(
            "El ciclo más largo tiene {} aristas. Cada 2-rotura bien elegida parte un ciclo en dos, y así la distancia baja de a 1.",
            miles(mayor)
        ));
    }
    if p.len() != q.len() {
        let diferencia = p.len().abs_diff(q.len());
        salida.explicar(format!(
            "P tiene {} y Q tiene {}: como cada 2-rotura cambia esa cuenta en 1 como mucho, {}.",
            cuantos(p.len(), "cromosoma", "cromosomas"),
            cuantos(q.len(), "cromosoma", "cromosomas"),
            if diferencia == 1 {
                "al menos un paso tiene que ser una fusión o una fisión".to_string()
            } else {
                format!(
                    "al menos {} pasos tienen que ser fusiones o fisiones",
                    miles(diferencia)
                )
            }
        ));
    }
    salida.explicar(
        "Para ver una ruta concreta paso a paso, usa «Transformación más corta con 2-roturas».",
    );
    Ok(salida)
}

fn dos_roturas_en_grafo(entrada: &Entrada) -> Result<Salida, String> {
    let aristas = leer_aristas(entrada, "aristas")?;
    let [i1, i2, i3, i4] = leer_cuatro_nodos(entrada, "nodos")?;
    let antes = reordenamientos::grafo_a_genoma(&aristas)?;
    let nuevas = reordenamientos::dos_roturas_en_grafo(&aristas, i1, i2, i3, i4)?;
    let despues = reordenamientos::grafo_a_genoma(&nuevas)?;
    let mut salida = Salida::default();
    salida.agregar(Bloque::texto(
        "Aristas de color nuevas",
        texto_aristas(&nuevas),
    ));
    salida.agregar(Bloque::tabla(
        "Cambio",
        &["Se quitan", "Se agregan"],
        vec![
            vec![format!("({i1}, {i2})"), format!("({i1}, {i3})")],
            vec![format!("({i3}, {i4})"), format!("({i2}, {i4})")],
        ],
    ));
    salida.agregar(Bloque::texto("Genoma antes", texto_genoma(&antes)));
    salida.agregar(Bloque::texto("Genoma después", texto_genoma(&despues)));
    agregar_grafo_del_genoma(
        &mut salida,
        "Grafo del genoma después de la 2-rotura",
        &nuevas,
    );
    salida.explicar(format!(
        "La 2-rotura corta las uniones ({i1}, {i2}) y ({i3}, {i4}) y vuelve a pegar los cuatro extremos de otra forma: ({i1}, {i3}) y ({i2}, {i4})."
    ));
    explicar_efecto(&mut salida, antes.len(), despues.len());
    salida.explicar("Los bloques no cambian: una 2-rotura solo cambia qué extremos están unidos.");
    Ok(salida)
}

fn dos_roturas_en_genoma(entrada: &Entrada) -> Result<Salida, String> {
    let genoma = leer_genoma(entrada, "genoma", MAXIMO_BLOQUES)?;
    let [i1, i2, i3, i4] = leer_cuatro_nodos(entrada, "nodos")?;
    let resultado = reordenamientos::dos_roturas_en_genoma(&genoma, i1, i2, i3, i4)?;
    let nuevas = reordenamientos::aristas_de_color(&resultado)?;
    let mut salida = Salida::default();
    salida.agregar(Bloque::texto("Genoma después", texto_genoma(&resultado)));
    salida.agregar(Bloque::valor("Cromosomas antes", miles(genoma.len())));
    salida.agregar(Bloque::valor("Cromosomas después", miles(resultado.len())));
    salida.agregar(Bloque::valor(
        "Tipo de cambio",
        efecto_de_la_rotura(genoma.len(), resultado.len()),
    ));
    agregar_grafo_del_genoma(&mut salida, "Grafo del genoma después", &nuevas);
    salida.explicar(format!(
        "Primero el genoma se pasa a aristas de color; luego se cortan ({i1}, {i2}) y ({i3}, {i4}), se pegan ({i1}, {i3}) y ({i2}, {i4}), y se vuelve a leer como cromosomas."
    ));
    explicar_efecto(&mut salida, genoma.len(), resultado.len());
    salida.explicar(
        "En genomas circulares no hay otra posibilidad: dos cortes en un mismo cromosoma dan una inversión o una fisión, y dos cortes en cromosomas distintos siempre dan una fusión.",
    );
    Ok(salida)
}

fn transformacion_mas_corta(entrada: &Entrada) -> Result<Salida, String> {
    let p = leer_genoma(entrada, "p", MAXIMO_LISTADO)?;
    let q = leer_genoma(entrada, "q", MAXIMO_LISTADO)?;
    let roturas = reordenamientos::transformacion_mas_corta(&p, &q)?;
    let ciclos = reordenamientos::ciclos_del_grafo_de_rupturas(&p, &q)?.len();
    let bloques = contar_bloques(&p);
    let distancia = roturas.len();
    let mut inversiones = 0;
    let mut fisiones = 0;
    let mut fusiones = 0;
    let mut filas = Vec::new();
    let mut cromosomas_antes = p.len();
    for (numero, rotura) in roturas.iter().enumerate() {
        let [i1, i2, i3, i4] = rotura.nodos;
        let efecto = efecto_de_la_rotura(cromosomas_antes, rotura.genoma.len());
        match efecto {
            "fisión" => fisiones += 1,
            "fusión" => fusiones += 1,
            _ => inversiones += 1,
        }
        cromosomas_antes = rotura.genoma.len();
        if filas.len() < MAXIMO_FILAS {
            filas.push(vec![
                miles(numero + 1),
                format!("({i1}, {i2}) y ({i3}, {i4})"),
                format!("({i1}, {i3}) y ({i2}, {i4})"),
                efecto.to_string(),
                miles(distancia - numero - 1),
            ]);
        }
    }
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("2-roturas", miles(distancia)));
    salida.agregar(Bloque::valor("Bloques", miles(bloques)));
    salida.agregar(Bloque::valor("Ciclos al inicio", miles(ciclos)));
    salida.agregar(lista_acotada(
        "Genomas desde P hasta Q",
        std::iter::once(texto_genoma(&p))
            .chain(roturas.iter().map(|rotura| texto_genoma(&rotura.genoma))),
        distancia + 1,
    ));
    if distancia > 0 {
        salida.agregar(Bloque::Tabla {
            etiqueta: "Cada 2-rotura".to_string(),
            columnas: vec![
                "Paso".to_string(),
                "Se cortan".to_string(),
                "Se pegan".to_string(),
                "Tipo".to_string(),
                "Distancia que queda".to_string(),
            ],
            filas,
            total: distancia,
        });
        salida.agregar(Bloque::barras(
            "Tipos de 2-rotura en esta ruta",
            "Tipo",
            "Cantidad",
            vec![
                "Inversiones".to_string(),
                "Fisiones".to_string(),
                "Fusiones".to_string(),
            ],
            vec![inversiones as f64, fisiones as f64, fusiones as f64],
        ));
    }
    if distancia == 0 {
        salida.explicar("P y Q ya son el mismo genoma: no hace falta ninguna 2-rotura.");
        return Ok(salida);
    }
    salida.explicar(format!(
        "La ruta tiene {}, justo la distancia bloques − ciclos = {} − {}: no existe una ruta más corta.",
        cuantos(distancia, "paso", "pasos"),
        miles(bloques),
        miles(ciclos)
    ));
    salida.explicar(
        "En cada paso se elige una arista azul (una unión de Q) que P todavía no tiene y se crea con una 2-rotura: aparece un ciclo trivial nuevo y la distancia baja en 1.",
    );
    salida.explicar(format!(
        "Entre esos pasos hay {}, {} y {}.",
        cuantos(inversiones, "inversión", "inversiones"),
        cuantos(fisiones, "fisión", "fisiones"),
        cuantos(fusiones, "fusión", "fusiones")
    ));
    salida.explicar(
        "Suele haber muchas rutas igual de cortas; esta elige siempre el nodo de menor número que todavía no coincide con Q.",
    );
    Ok(salida)
}

fn estado_de_transformacion(
    paso: usize,
    fase: &str,
    genoma: &[Vec<i64>],
    q: &[Vec<i64>],
    ciclos: usize,
    rotura: Option<[usize; 4]>,
    efecto: Option<&str>,
) -> Result<Value, String> {
    let bloques = contar_bloques(genoma);
    let rojas = reordenamientos::aristas_de_color(genoma)?;
    let azules = reordenamientos::aristas_de_color(q)?;
    Ok(json!({
        "paso": paso,
        "fase": fase,
        "genoma": genoma,
        "texto": texto_genoma(genoma),
        "objetivo": q,
        "texto_objetivo": texto_genoma(q),
        "bloques": bloques,
        "ciclos": ciclos,
        "distancia": bloques - ciclos,
        "rotura": rotura,
        "efecto": efecto,
        "rojas": rojas,
        "azules": azules,
    }))
}

fn pasos_transformacion(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let p = leer_genoma(entrada, "p", MAXIMO_LISTADO)?;
    let q = leer_genoma(entrada, "q", MAXIMO_LISTADO)?;
    let bloques = contar_bloques(&p);
    if bloques > MAXIMO_ANIMACION {
        return Err(format!(
            "Para el paso a paso usa genomas de hasta {MAXIMO_ANIMACION} bloques; P tiene {}.",
            miles(bloques)
        ));
    }
    let roturas: Vec<DosRoturas> = reordenamientos::transformacion_mas_corta(&p, &q)?;
    let ciclos_iniciales = reordenamientos::ciclos_del_grafo_de_rupturas(&p, &q)?.len();
    let mut pasos = vec![Paso::nuevo(
        "let azules = pareja_de_cada_nodo(&aristas_entre_vecinos(q), bloques);",
        format!(
            "Grafo de rupturas: aristas rojas de P y azules de Q. Hay {ciclos_iniciales} ciclos, así que la distancia es {}.",
            bloques - ciclos_iniciales
        ),
        estado_de_transformacion(0, "inicio", &p, &q, ciclos_iniciales, None, None)?,
    )];
    let mut actual = p.clone();
    for (numero, rotura) in roturas.iter().enumerate() {
        let [i1, i2, i3, i4] = rotura.nodos;
        let ciclos_antes = ciclos_iniciales + numero;
        pasos.push(Paso::nuevo(
            "if rojas[nodo] != azules[nodo] {",
            format!(
                "El nodo {i1} está unido a {i2} en P pero a {i3} en Q: su ciclo no es trivial."
            ),
            estado_de_transformacion(
                numero,
                "elige",
                &actual,
                &q,
                ciclos_antes,
                Some(rotura.nodos),
                None,
            )?,
        ));
        let efecto = efecto_de_la_rotura(actual.len(), rotura.genoma.len());
        pasos.push(Paso::nuevo(
            "reconectar(&mut rojas, i1, i2, i3, i4);",
            format!(
                "2-rotura ({efecto}): se cortan ({i1}, {i2}) y ({i3}, {i4}) y se pegan ({i1}, {i3}) y ({i2}, {i4}). Ahora hay {} ciclos.",
                ciclos_antes + 1
            ),
            estado_de_transformacion(
                numero + 1,
                "aplica",
                &rotura.genoma,
                &q,
                ciclos_antes + 1,
                Some(rotura.nodos),
                Some(efecto),
            )?,
        ));
        actual = rotura.genoma.clone();
    }
    pasos.push(Paso::nuevo(
        "Ok(roturas)",
        format!(
            "Todos los ciclos son triviales: P se convirtió en Q en {}.",
            cuantos(roturas.len(), "paso", "pasos")
        ),
        estado_de_transformacion(
            roturas.len(),
            "fin",
            &actual,
            &q,
            ciclos_iniciales + roturas.len(),
            None,
            None,
        )?,
    ));
    Ok(muestrear(&pasos, MAXIMO_PASOS))
}

fn inicio_en_ecoli(entrada: &Entrada, campo: &str) -> Option<usize> {
    let texto = entrada.texto(campo).ok()?;
    let limpio = texto.trim();
    if limpio == "@ecoli" {
        return Some(1);
    }
    let rango = limpio.strip_prefix("@ecoli:")?;
    let (inicio, _) = rango.split_once('-')?;
    inicio.trim().parse().ok()
}

fn leer_secuencias_y_k(entrada: &Entrada) -> Result<(Vec<u8>, Vec<u8>, usize), String> {
    let a = entrada.adn("a")?;
    let b = entrada.adn("b")?;
    let k = entrada.natural("k")?;
    if k == 0 || k > 32 {
        return Err(format!("k debe estar entre 1 y 32, y es {k}."));
    }
    if a.len() > MAXIMO_A {
        return Err(format!(
            "La secuencia A puede tener hasta {} pb y tiene {}.",
            miles(MAXIMO_A),
            miles(a.len())
        ));
    }
    if b.len() > MAXIMO_B {
        return Err(format!(
            "La secuencia B se guarda en un índice y puede tener hasta {} pb; tiene {}. Si quieres comparar contra todo el genoma, ponlo en A.",
            miles(MAXIMO_B),
            miles(b.len())
        ));
    }
    for (nombre, secuencia) in [("A", &a), ("B", &b)] {
        if k > secuencia.len() {
            return Err(format!(
                "k = {k} es mayor que la longitud de la secuencia {nombre} ({} pb).",
                miles(secuencia.len())
            ));
        }
    }
    Ok((a, b, k))
}

fn posicion_mostrada(posicion: usize, inicio: Option<usize>) -> String {
    match inicio {
        Some(base) => format!("{} (E. coli {})", miles(posicion), miles(base + posicion)),
        None => miles(posicion),
    }
}

fn kmeros_compartidos(entrada: &Entrada) -> Result<Salida, String> {
    let (a, b, k) = leer_secuencias_y_k(entrada)?;
    let pares = reordenamientos::kmeros_compartidos(&a, &b, k)?;
    let inicio_a = inicio_en_ecoli(entrada, "a");
    let inicio_b = inicio_en_ecoli(entrada, "b");
    let invertidos = pares.iter().filter(|par| par.invertido).count();
    let directos = pares.len() - invertidos;
    let mut posiciones_de_a: Vec<usize> = pares.iter().map(|par| par.posicion_a).collect();
    posiciones_de_a.dedup();
    let kmeros_de_a = a.len() - k + 1;
    let kmeros_de_b = b.len() - k + 1;
    let esperadas = 2.0 * kmeros_de_a as f64 * kmeros_de_b as f64 / 4f64.powi(k as i32);
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        "Parejas de k-meros compartidos",
        miles(pares.len()),
    ));
    salida.agregar(Bloque::valor("En la misma hebra", miles(directos)));
    salida.agregar(Bloque::valor("En la hebra contraria", miles(invertidos)));
    let puntos: Vec<Punto> = muestrear(&pares, MAXIMO_PUNTOS)
        .iter()
        .map(|par| Punto {
            x: par.posicion_a as f64,
            y: par.posicion_b as f64,
            grupo: if par.invertido { 1 } else { 0 },
            nombre: String::from_utf8_lossy(&a[par.posicion_a..par.posicion_a + k]).to_string(),
        })
        .collect();
    salida.agregar(Bloque::dispersion(
        "Gráfico de puntos (grupo 0: misma hebra, grupo 1: hebra contraria)",
        "Posición en A (contando desde 0)",
        "Posición en B (contando desde 0)",
        puntos,
    ));
    let filas: Vec<Vec<String>> = pares
        .iter()
        .take(MAXIMO_FILAS)
        .map(|par| {
            vec![
                posicion_mostrada(par.posicion_a, inicio_a),
                posicion_mostrada(par.posicion_b, inicio_b),
                String::from_utf8_lossy(&a[par.posicion_a..par.posicion_a + k]).to_string(),
                if par.invertido { "contraria" } else { "misma" }.to_string(),
            ]
        })
        .collect();
    salida.agregar(Bloque::Tabla {
        etiqueta: "Parejas (posiciones contando desde 0)".to_string(),
        columnas: vec![
            "Posición en A".to_string(),
            "Posición en B".to_string(),
            "k-mero en A".to_string(),
            "Hebra en B".to_string(),
        ],
        filas,
        total: pares.len(),
    });
    salida.explicar(format!(
        "Hay {} parejas de {k}-meros compartidos: {} en la misma hebra y {} en la hebra contraria (B tiene el complemento inverso).",
        miles(pares.len()),
        miles(directos),
        miles(invertidos)
    ));
    let umbral = esperadas + 3.0 * esperadas.sqrt();
    if pares.is_empty() {
        salida.explicar(format!(
            "Si A y B fueran ADN al azar esperarías unas {} parejas por casualidad y no hay ninguna: no comparten ningún {k}-mero. Prueba con un k menor o con regiones que sí se parezcan.",
            formato_esperadas(esperadas)
        ));
    } else if pares.len() as f64 > umbral {
        salida.explicar(format!(
            "Si A y B fueran ADN al azar esperarías unas {} parejas por casualidad; hay muchas más (más de 3 desviaciones), así que comparten tramos de verdad.",
            formato_esperadas(esperadas)
        ));
    } else {
        salida.explicar(format!(
            "Si A y B fueran ADN al azar esperarías unas {} parejas por casualidad: lo observado no se distingue del azar con este k.",
            formato_esperadas(esperadas)
        ));
    }
    if !pares.is_empty() {
        if invertidos == 0 {
            salida.explicar(
                "Todos los puntos son de la misma hebra: forman diagonales que suben hacia la derecha, señal de tramos en el mismo orden y sentido en A y en B.",
            );
        } else if directos == 0 {
            salida.explicar(
                "Todos los puntos son de la hebra contraria: forman antidiagonales que bajan hacia la derecha, la marca de un tramo invertido.",
            );
        } else if directos == invertidos {
            salida.explicar(
                "Hay tantos puntos de la misma hebra como de la hebra contraria: mira el gráfico, porque las diagonales que suben son tramos en el mismo sentido y las antidiagonales que bajan son tramos invertidos.",
            );
        } else if directos > invertidos {
            salida.explicar(format!(
                "La mayoría ({}) son de la misma hebra y dibujan diagonales que suben (mismo sentido); los de la hebra contraria dibujan antidiagonales (tramos invertidos).",
                porcentaje(directos as f64 / pares.len() as f64, 1)
            ));
        } else {
            salida.explicar(format!(
                "La mayoría ({}) son de la hebra contraria y dibujan antidiagonales que bajan (tramos invertidos); los de la misma hebra dibujan diagonales.",
                porcentaje(invertidos as f64 / pares.len() as f64, 1)
            ));
        }
        salida.explicar(format!(
            "El {} de los k-meros de A aparece en B. Un hueco en una línea es un tramo que cambió; para agrupar los puntos en bloques usa «Bloques de sintenia».",
            porcentaje(posiciones_de_a.len() as f64 / kmeros_de_a as f64, 1)
        ));
    }
    Ok(salida)
}

fn formato_esperadas(esperadas: f64) -> String {
    if esperadas >= 10.0 {
        miles(esperadas.round() as usize)
    } else if esperadas >= 0.01 {
        format!("{esperadas:.2}")
    } else {
        format!("{esperadas:.1e}")
    }
}

fn se_solapan(bloques: &[BloqueSintenico]) -> bool {
    let mut por_a: Vec<(usize, usize)> = bloques
        .iter()
        .map(|bloque| (bloque.inicio_a, bloque.fin_a))
        .collect();
    let mut por_b: Vec<(usize, usize)> = bloques
        .iter()
        .map(|bloque| (bloque.inicio_b, bloque.fin_b))
        .collect();
    por_a.sort();
    por_b.sort();
    let solapa = |tramos: &Vec<(usize, usize)>| tramos.windows(2).any(|par| par[1].0 <= par[0].1);
    solapa(&por_a) || solapa(&por_b)
}

fn largo_cubierto(mut tramos: Vec<(usize, usize)>) -> usize {
    tramos.sort();
    let mut total = 0;
    let mut hasta_ahora: Option<(usize, usize)> = None;
    for (inicio, fin) in tramos {
        match hasta_ahora {
            Some((abierto, cerrado)) if inicio <= cerrado + 1 => {
                hasta_ahora = Some((abierto, cerrado.max(fin)));
            }
            Some((abierto, cerrado)) => {
                total += cerrado - abierto + 1;
                hasta_ahora = Some((inicio, fin));
            }
            None => hasta_ahora = Some((inicio, fin)),
        }
    }
    if let Some((abierto, cerrado)) = hasta_ahora {
        total += cerrado - abierto + 1;
    }
    total
}

fn bloques_de_sintenia(entrada: &Entrada) -> Result<Salida, String> {
    let (a, b, k) = leer_secuencias_y_k(entrada)?;
    let distancia = entrada.natural("distancia")?;
    let minimo = entrada.natural("minimo")?;
    if distancia == 0 {
        return Err("La distancia para unir k-meros vecinos debe ser al menos 1 pb.".to_string());
    }
    let sintenia = reordenamientos::bloques_de_sintenia(&a, &b, k, distancia, minimo)?;
    let bloques = &sintenia.bloques;
    let inicio_a = inicio_en_ecoli(entrada, "a");
    let inicio_b = inicio_en_ecoli(entrada, "b");
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Bloques de sintenia", miles(bloques.len())));
    salida.agregar(Bloque::valor(
        "k-meros compartidos",
        miles(sintenia.pares.len()),
    ));
    let permutacion = reordenamientos::permutacion_de_bloques(bloques);
    if !bloques.is_empty() && bloques.len() <= MAXIMO_LISTADO {
        salida.agregar(Bloque::texto(
            "Bloques en el orden de B (numerados según A)",
            texto_permutacion(&permutacion),
        ));
    }
    let filas: Vec<Vec<String>> = bloques
        .iter()
        .take(MAXIMO_FILAS)
        .enumerate()
        .map(|(indice, bloque)| {
            vec![
                miles(indice + 1),
                posicion_mostrada(bloque.inicio_a, inicio_a),
                posicion_mostrada(bloque.fin_a, inicio_a),
                posicion_mostrada(bloque.inicio_b, inicio_b),
                posicion_mostrada(bloque.fin_b, inicio_b),
                format!("{} pb", miles(bloque.fin_a - bloque.inicio_a + 1)),
                if bloque.invertido { "−" } else { "+" }.to_string(),
                miles(bloque.kmeros),
            ]
        })
        .collect();
    salida.agregar(Bloque::Tabla {
        etiqueta: "Bloques (posiciones contando desde 0)".to_string(),
        columnas: vec![
            "Bloque".to_string(),
            "Inicio en A".to_string(),
            "Fin en A".to_string(),
            "Inicio en B".to_string(),
            "Fin en B".to_string(),
            "Largo en A".to_string(),
            "Sentido".to_string(),
            "k-meros".to_string(),
        ],
        filas,
        total: bloques.len(),
    });
    let muestra = muestrear(&sintenia.pares, MAXIMO_PUNTOS);
    let grupos = reordenamientos::bloque_de_cada_par(&muestra, bloques);
    let puntos: Vec<Punto> = muestra
        .iter()
        .zip(grupos)
        .map(|(par, grupo)| Punto {
            x: par.posicion_a as f64,
            y: par.posicion_b as f64,
            grupo,
            nombre: if grupo == 0 {
                "fuera de bloques".to_string()
            } else {
                format!("bloque {grupo}")
            },
        })
        .collect();
    salida.agregar(Bloque::dispersion(
        "Gráfico de puntos por bloque (grupo 0: fuera de todo bloque)",
        "Posición en A (contando desde 0)",
        "Posición en B (contando desde 0)",
        puntos,
    ));
    if bloques.is_empty() {
        salida.explicar(format!(
            "No se formó ningún bloque con al menos {} k-meros. Prueba con un k menor, una distancia mayor o un mínimo más bajo.",
            miles(minimo)
        ));
        return Ok(salida);
    }
    salida.explicar(format!(
        "{}: grupos de al menos {} k-meros compartidos que están a menos de unas {} pb entre sí en el gráfico de puntos.",
        if bloques.len() == 1 {
            "Se formó 1 bloque de sintenia".to_string()
        } else {
            format!("Se formaron {} bloques de sintenia", miles(bloques.len()))
        },
        miles(minimo),
        miles(distancia)
    ));
    let invertidos = bloques.iter().filter(|bloque| bloque.invertido).count();
    if invertidos == bloques.len() {
        salida.explicar(
            "Todos los bloques están en la hebra contraria de B (signo −): en el gráfico son antidiagonales que bajan.",
        );
    } else if invertidos > 0 {
        salida.explicar(format!(
            "{} de los {} bloques están en la hebra contraria de B (signo −): en el gráfico son antidiagonales que bajan.",
            miles(invertidos),
            miles(bloques.len())
        ));
    }
    let rupturas = reordenamientos::numero_de_rupturas(&permutacion)?;
    let orden = if permutacion.len() <= MAXIMO_EN_FRASE {
        texto_permutacion(&permutacion)
    } else if permutacion.len() <= MAXIMO_LISTADO {
        "que ves arriba".to_string()
    } else {
        "de todos los bloques".to_string()
    };
    if rupturas == 0 {
        salida.explicar(format!(
            "Leídos en el orden de B, los bloques forman la permutación {orden}: sin rupturas, A y B conservan el mismo orden y sentido."
        ));
    } else {
        salida.explicar(format!(
            "Leídos en el orden de B, los bloques forman la permutación {orden}, con {}: {} para explicar ese orden.",
            cuantos(rupturas, "ruptura", "rupturas"),
            if rupturas.div_ceil(2) == 1 {
                "hace falta al menos 1 inversión".to_string()
            } else {
                format!(
                    "hacen falta al menos {} inversiones",
                    miles(rupturas.div_ceil(2))
                )
            }
        ));
    }
    if se_solapan(bloques) {
        salida.explicar(
            "Algunos bloques se solapan en A o en B: son copias repetidas de un mismo tramo, y entonces el orden de bloques no es una permutación limpia.",
        );
    }
    let cubierto = largo_cubierto(
        bloques
            .iter()
            .map(|bloque| (bloque.inicio_a, bloque.fin_a))
            .collect(),
    );
    salida.explicar(format!(
        "Los bloques cubren el {} de A; el resto no se parece a B a esta escala, y cada hueco entre bloques marca dónde las dos secuencias dejan de coincidir.",
        porcentaje(cubierto as f64 / a.len() as f64, 1)
    ));
    Ok(salida)
}

const VENTANA_LISTA: usize = 24;
const ELEMENTOS_EN_VARIABLE: usize = 12;

fn con_signo(bloque: &i64) -> String {
    format!("{bloque:+}")
}

fn entero_escrito(valor: &i64) -> String {
    valor.to_string()
}

fn nodo_escrito(nodo: &usize) -> String {
    nodo.to_string()
}

fn arista_escrita(arista: &(usize, usize)) -> String {
    format!("({}, {})", arista.0, arista.1)
}

fn resumen<T>(elementos: &[T], escribir: fn(&T) -> String, separador: &str) -> String {
    let piezas: Vec<String> = elementos
        .iter()
        .take(ELEMENTOS_EN_VARIABLE)
        .map(escribir)
        .collect();
    let cola = if elementos.len() > ELEMENTOS_EN_VARIABLE {
        " …"
    } else {
        ""
    };
    format!("{}{cola}", piezas.join(separador))
}

fn corto(texto: &str) -> String {
    let mut letras = texto.chars();
    let inicio: String = letras.by_ref().take(48).collect();
    if letras.next().is_some() {
        format!("{inicio}…")
    } else {
        inicio
    }
}

fn base_de(variables: &[(&str, String)]) -> Pizarra {
    variables
        .iter()
        .fold(Pizarra::nueva(), |pizarra, (nombre, valor)| {
            pizarra.variable(nombre, valor)
        })
}

fn con_lista<T>(
    pizarra: Pizarra,
    titulo: &str,
    elementos: &[T],
    escribir: fn(&T) -> String,
    centro: usize,
    resaltar: &[usize],
) -> Pizarra {
    let total = elementos.len();
    let largo = total.min(VENTANA_LISTA);
    let desde = centro.saturating_sub(VENTANA_LISTA / 3).min(total - largo);
    let hasta = desde + largo;
    let titulo = if largo < total {
        format!(
            "{titulo} (posiciones {}–{} de {})",
            miles(desde + 1),
            miles(hasta),
            miles(total)
        )
    } else {
        titulo.to_string()
    };
    let visibles: Vec<String> = elementos[desde..hasta].iter().map(escribir).collect();
    let marcas: Vec<usize> = resaltar
        .iter()
        .filter(|&&indice| indice >= desde && indice < hasta)
        .map(|indice| indice - desde)
        .collect();
    pizarra.lista(&titulo, visibles, &marcas)
}

enum Toque {
    Detalle,
    Resumen,
    Nada,
}

struct Muestreo {
    visto: usize,
    detalle: usize,
    salto: usize,
}

impl Muestreo {
    fn nuevo(total: usize, detalle: usize, objetivo: usize) -> Muestreo {
        Muestreo {
            visto: 0,
            detalle,
            salto: (total / objetivo.max(1)).max(1),
        }
    }

    fn toca(&mut self) -> Toque {
        self.visto += 1;
        if self.visto <= self.detalle {
            Toque::Detalle
        } else if self.visto.is_multiple_of(self.salto) {
            Toque::Resumen
        } else {
            Toque::Nada
        }
    }
}

fn parentesis(texto: String) -> String {
    format!("({texto})")
}

fn tablero_voraz(
    actual: &[i64],
    guardadas: usize,
    foco: (Option<usize>, Option<i64>, Option<usize>),
    resaltar: &[usize],
) -> Pizarra {
    let mut variables = vec![
        ("actual", parentesis(resumen(actual, con_signo, " "))),
        ("permutaciones", format!("{} guardadas", miles(guardadas))),
    ];
    if let Some(lugar) = foco.0 {
        variables.push(("lugar", lugar.to_string()));
    }
    if let Some(bloque) = foco.1 {
        variables.push(("bloque", bloque.to_string()));
    }
    if let Some(donde) = foco.2 {
        variables.push(("donde", donde.to_string()));
    }
    let centro = resaltar.first().copied().unwrap_or(0);
    con_lista(
        base_de(&variables),
        "permutación actual",
        actual,
        con_signo,
        centro,
        resaltar,
    )
}

fn narrar_voraz(permutacion: &[i64]) -> Result<(Vec<Paso>, Vec<i64>, usize), String> {
    let esperado = reordenamientos::ordenamiento_voraz(permutacion)?;
    let n = permutacion.len();
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let mut muestreo = Muestreo::nuevo(n, 3, 80);
    let mut actual = permutacion.to_vec();
    let mut guardadas = 0;
    narracion.agregar(
        tablero_voraz(&actual, 0, (None, None, None), &[]).paso(
            "validar_bloques(permutacion)?;",
            format!(
                "Tienes {} en este orden. La meta es la identidad (+1 +2 … +{n}). El método voraz arregla un lugar a la vez, de izquierda a derecha. Antes, el código comprueba que sean los números del 1 al {n}, sin repetir: aquí lo son.",
                cuantos(n, "bloque", "bloques")
            ),
        ),
    );
    narracion.agregar(
        tablero_voraz(&actual, 0, (None, None, None), &[]).paso(
            "let mut actual = permutacion.to_vec();",
            "Trabajamos sobre una copia, actual. Cada vez que la cambiemos, guardaremos una foto en permutaciones.",
        ),
    );
    for lugar in 0..n {
        let bloque = lugar as i64 + 1;
        let toque = muestreo.toca();
        let detalle = matches!(toque, Toque::Detalle);
        let avisar = !matches!(toque, Toque::Nada);
        let estaba = actual[lugar].abs() == bloque;
        if detalle {
            narracion.agregar(
                tablero_voraz(&actual, guardadas, (Some(lugar), Some(bloque), None), &[lugar]).paso(
                    "for lugar in 0..actual.len() {",
                    format!(
                        "Toca el lugar {lugar}, la posición {}: ahí tiene que quedar el bloque +{bloque}.",
                        lugar + 1
                    ),
                ),
            );
            let nota = if estaba {
                format!(
                    "Ahí ya está el bloque {bloque} (con signo {}): no hace falta traerlo.",
                    if actual[lugar] > 0 { "+" } else { "−" }
                )
            } else {
                format!(
                    "Ahí hay un {}, no el {bloque}: hay que traer el {bloque} hasta este lugar.",
                    con_signo(&actual[lugar])
                )
            };
            narracion.agregar(
                tablero_voraz(
                    &actual,
                    guardadas,
                    (Some(lugar), Some(bloque), None),
                    &[lugar],
                )
                .paso("if actual[lugar].abs() != bloque {", nota),
            );
        }
        let mut hubo_cambio = false;
        if !estaba {
            let donde = reordenamientos::posicion_del_bloque(&actual, bloque, lugar);
            if detalle {
                narracion.agregar(
                    tablero_voraz(
                        &actual,
                        guardadas,
                        (Some(lugar), Some(bloque), Some(donde)),
                        &[lugar, donde],
                    )
                    .paso(
                        "let donde = posicion_del_bloque(&actual, bloque, lugar);",
                        format!(
                            "Buscamos hacia la derecha: el bloque {bloque} está en la posición {} (en el código, donde = {donde}).",
                            donde + 1
                        ),
                    ),
                );
            }
            reordenamientos::invertir_tramo(&mut actual, lugar, donde);
            guardadas += 1;
            hubo_cambio = true;
            if avisar {
                let marcados: Vec<usize> = (lugar..=donde).collect();
                let nota = if detalle {
                    format!(
                        "Damos vuelta el tramo de la posición {} a la {}: se invierte su orden y cada bloque cambia de signo. Ahora el {bloque} está en su lugar, aunque puede haber quedado negativo. Se guarda esta permutación.",
                        lugar + 1,
                        donde + 1
                    )
                } else {
                    format!(
                        "Lugar {lugar}: el {bloque} estaba en la posición {}; se da vuelta ese tramo de {} y el {bloque} queda en su lugar.",
                        donde + 1,
                        cuantos(donde - lugar + 1, "bloque", "bloques")
                    )
                };
                narracion.agregar(
                    tablero_voraz(
                        &actual,
                        guardadas,
                        (Some(lugar), Some(bloque), Some(donde)),
                        &marcados,
                    )
                    .paso("invertir_tramo(&mut actual, lugar, donde);", nota),
                );
            }
        }
        if detalle {
            let negativo = actual[lugar] == -bloque;
            narracion.agregar(
                tablero_voraz(&actual, guardadas, (Some(lugar), Some(bloque), None), &[lugar]).paso(
                    "if actual[lugar] == -bloque {",
                    if negativo {
                        format!(
                            "Quedó el −{bloque}: está en su lugar pero al revés, así que falta cambiarle el signo."
                        )
                    } else {
                        format!("El {bloque} quedó con signo +: este lugar está listo.")
                    },
                ),
            );
        }
        if actual[lugar] == -bloque {
            actual[lugar] = bloque;
            guardadas += 1;
            hubo_cambio = true;
            if avisar {
                narracion.agregar(
                    tablero_voraz(&actual, guardadas, (Some(lugar), Some(bloque), None), &[lugar])
                        .paso(
                            "actual[lugar] = bloque;",
                            format!(
                                "Se le cambia el signo al bloque del lugar {lugar}: ahora es +{bloque}. Eso cuenta como una inversión más, de un solo bloque, y se guarda."
                            ),
                        ),
                );
            }
        }
        if avisar && !hubo_cambio && !detalle {
            narracion.agregar(
                tablero_voraz(
                    &actual,
                    guardadas,
                    (Some(lugar), Some(bloque), None),
                    &[lugar],
                )
                .paso(
                    "for lugar in 0..actual.len() {",
                    format!("Lugar {lugar}: el +{bloque} ya estaba en su sitio, no se hace nada."),
                ),
            );
        }
    }
    let ultimo = tablero_voraz(&actual, guardadas, (None, None, None), &[]).paso(
        "Ok(permutaciones)",
        format!(
            "Terminamos: {} bastaron para dejar la permutación en +1 +2 … +{n}. Compáralo con el mínimo teórico en «Qué significa»: el voraz no siempre lo alcanza.",
            cuantos(guardadas, "inversión", "inversiones")
        ),
    );
    debug_assert_eq!(guardadas, esperado.len());
    Ok((narracion.terminar(ultimo), actual, guardadas))
}

fn pasos_ordenamiento_voraz(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let permutacion = leer_permutacion(entrada, "permutacion", MAXIMO_LISTADO)?;
    Ok(narrar_voraz(&permutacion)?.0)
}

fn par_escrito(par: Option<(i64, i64)>) -> String {
    match par {
        Some((izquierda, derecha)) => format!("({izquierda}, {derecha})"),
        None => "—".to_string(),
    }
}

fn tablero_rupturas(
    permutacion: &[i64],
    extendida: &[i64],
    encontradas: usize,
    vecinos: Option<(i64, i64)>,
    resaltar: &[usize],
) -> Pizarra {
    let mut variables = vec![(
        "permutacion",
        parentesis(resumen(permutacion, con_signo, " ")),
    )];
    if !extendida.is_empty() {
        variables.push((
            "extendida",
            parentesis(resumen(extendida, entero_escrito, " ")),
        ));
        variables.push(("encontradas", miles(encontradas)));
    }
    if vecinos.is_some() {
        variables.push(("vecinos", par_escrito(vecinos)));
    }
    let centro = resaltar.first().copied().unwrap_or(0);
    if extendida.is_empty() {
        con_lista(
            base_de(&variables),
            "permutación",
            permutacion,
            con_signo,
            0,
            &[],
        )
    } else {
        con_lista(
            base_de(&variables),
            "permutación con bordes",
            extendida,
            entero_escrito,
            centro,
            resaltar,
        )
    }
}

type Aristas = Vec<(usize, usize)>;
type NarradoDeRupturas = Result<(Vec<Paso>, Vec<(i64, i64)>), String>;
type NarradoDeAristas = Result<(Vec<Paso>, Aristas), String>;
type NarradoDeDistancia = Result<(Vec<Paso>, usize, Vec<Vec<usize>>), String>;

fn narrar_rupturas(permutacion: &[i64]) -> NarradoDeRupturas {
    let esperadas = reordenamientos::rupturas(permutacion)?;
    let n = permutacion.len();
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        tablero_rupturas(permutacion, &[], 0, None, &[]).paso(
            "validar_bloques(permutacion)?;",
            format!(
                "Tienes {}. Una ruptura es un par de vecinos que no son consecutivos (el de la derecha no vale 1 más que el de la izquierda). Primero se comprueba que sean los números del 1 al {n}, sin repetir.",
                cuantos(n, "bloque", "bloques")
            ),
        ),
    );
    let mut extendida: Vec<i64> = vec![0];
    narracion.agregar(
        tablero_rupturas(permutacion, &extendida, 0, None, &[0]).paso(
            "extendida.push(0);",
            "Se pone un 0 al principio: es un bloque imaginario que sirve de borde izquierdo.",
        ),
    );
    extendida.extend_from_slice(permutacion);
    narracion.agregar(
        tablero_rupturas(permutacion, &extendida, 0, None, &[]).paso(
            "extendida.extend_from_slice(permutacion);",
            "Después se copian los bloques de la permutación, con sus signos.",
        ),
    );
    extendida.push(n as i64 + 1);
    narracion.agregar(
        tablero_rupturas(permutacion, &extendida, 0, None, &[n + 1]).paso(
            "extendida.push(permutacion.len() as i64 + 1);",
            format!(
                "Y un {} al final, como borde derecho. En la permutación ya ordenada cada vecino vale 1 más que el de su izquierda: 0 1 2 … {}.",
                n + 1,
                n + 1
            ),
        ),
    );
    let mut encontradas: Vec<(i64, i64)> = Vec::new();
    let mut muestreo = Muestreo::nuevo(esperadas.len(), 0, 40);
    for (indice, vecinos) in extendida.windows(2).enumerate() {
        let par = (vecinos[0], vecinos[1]);
        let diferencia = par.1 - par.0;
        let ruptura = diferencia != 1;
        let marcados = [indice, indice + 1];
        if indice < 5 {
            narracion.agregar(
                tablero_rupturas(
                    permutacion,
                    &extendida,
                    encontradas.len(),
                    Some(par),
                    &marcados,
                )
                .paso(
                    "for vecinos in extendida.windows(2) {",
                    format!(
                        "Miramos el par de vecinos número {}: ({}, {}).",
                        indice + 1,
                        par.0,
                        par.1
                    ),
                ),
            );
            narracion.agregar(
                tablero_rupturas(
                    permutacion,
                    &extendida,
                    encontradas.len(),
                    Some(par),
                    &marcados,
                )
                .paso(
                    "if vecinos[1] - vecinos[0] != 1 {",
                    format!(
                        "{} − ({}) = {diferencia}. {}",
                        par.1,
                        par.0,
                        if ruptura {
                            "No es 1: hay una ruptura."
                        } else {
                            "Es 1: son consecutivos, así que es una adyacencia y no se cuenta."
                        }
                    ),
                ),
            );
        }
        if ruptura {
            encontradas.push(par);
            let toque = if indice < 5 {
                Toque::Detalle
            } else {
                muestreo.toca()
            };
            if !matches!(toque, Toque::Nada) {
                narracion.agregar(
                    tablero_rupturas(
                        permutacion,
                        &extendida,
                        encontradas.len(),
                        Some(par),
                        &marcados,
                    )
                    .paso(
                        "encontradas.push((vecinos[0], vecinos[1]));",
                        format!(
                            "El par ({}, {}) de la posición {} es una ruptura: se guarda. Llevas {}.",
                            par.0,
                            par.1,
                            indice + 1,
                            cuantos(encontradas.len(), "ruptura", "rupturas")
                        ),
                    ),
                );
            }
        }
    }
    let ultimo = tablero_rupturas(permutacion, &extendida, encontradas.len(), None, &[]).paso(
        "Ok(rupturas(permutacion)?.len())",
        format!(
            "Se recorrieron los {} pares de vecinos: {} y {}. numero_de_rupturas solo devuelve cuántas son: {}. Cada inversión elimina como mucho 2 rupturas; mira «Qué significa» para la cota.",
            miles(n + 1),
            cuantos(encontradas.len(), "ruptura", "rupturas"),
            cuantos(n + 1 - encontradas.len(), "adyacencia", "adyacencias"),
            miles(encontradas.len())
        ),
    );
    Ok((narracion.terminar(ultimo), encontradas))
}

fn pasos_numero_de_rupturas(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let permutacion = leer_permutacion(entrada, "permutacion", MAXIMO_BLOQUES)?;
    Ok(narrar_rupturas(&permutacion)?.0)
}

fn tablero_ciclo(
    cromosoma: &[i64],
    nodos: &[usize],
    foco: Option<(i64, Option<usize>, Option<usize>)>,
    resaltar: &[usize],
) -> Pizarra {
    let mut variables = vec![("cromosoma", parentesis(resumen(cromosoma, con_signo, " ")))];
    if !nodos.is_empty() || foco.is_some() {
        variables.push(("nodos", parentesis(resumen(nodos, nodo_escrito, " "))));
    }
    if let Some((bloque, cola, cabeza)) = foco {
        variables.push(("bloque", format!("{bloque:+}")));
        if let Some(cola) = cola {
            variables.push(("cola", cola.to_string()));
        }
        if let Some(cabeza) = cabeza {
            variables.push(("cabeza", cabeza.to_string()));
        }
    }
    let centro = resaltar.first().copied().unwrap_or(0);
    con_lista(
        base_de(&variables),
        "nodos del ciclo",
        nodos,
        nodo_escrito,
        centro,
        resaltar,
    )
}

fn narrar_cromosoma_a_ciclo(cromosoma: &[i64]) -> Result<(Vec<Paso>, Vec<usize>), String> {
    let esperado = reordenamientos::cromosoma_a_ciclo(cromosoma)?;
    let n = cromosoma.len();
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        tablero_ciclo(cromosoma, &[], None, &[]).paso(
            "validar_bloques(cromosoma)?;",
            format!(
                "Tienes un cromosoma de {}. Vamos a convertir cada bloque i en dos nodos: su cola, 2i−1, y su cabeza, 2i. Antes se comprueba que sean los números del 1 al {n}, sin repetir.",
                cuantos(n, "bloque", "bloques")
            ),
        ),
    );
    let mut nodos: Vec<usize> = Vec::with_capacity(2 * n);
    narracion.agregar(
        tablero_ciclo(cromosoma, &nodos, None, &[]).paso(
            "let mut nodos = Vec::with_capacity(2 * cromosoma.len());",
            format!(
                "Preparamos la lista de nodos, vacía por ahora. Al final tendrá {} nodos, dos por bloque.",
                miles(2 * n)
            ),
        ),
    );
    let mut muestreo = Muestreo::nuevo(n, 3, 60);
    for (posicion, &bloque) in cromosoma.iter().enumerate() {
        let cola = 2 * bloque.unsigned_abs() as usize - 1;
        let cabeza = 2 * bloque.unsigned_abs() as usize;
        let toque = muestreo.toca();
        let foco = Some((bloque, Some(cola), Some(cabeza)));
        let anterior = nodos.len();
        if matches!(toque, Toque::Detalle) {
            narracion.agregar(
                tablero_ciclo(cromosoma, &nodos, Some((bloque, None, None)), &[]).paso(
                    "for &bloque in cromosoma {",
                    format!(
                        "Siguiente bloque: el {bloque:+}, en la posición {} del cromosoma.",
                        posicion + 1
                    ),
                ),
            );
            narracion.agregar(
                tablero_ciclo(cromosoma, &nodos, Some((bloque, Some(cola), None)), &[]).paso(
                    "let cola = nodo_cola(bloque);",
                    format!(
                        "Su cola es 2·{} − 1 = {cola}: el signo no importa, solo el número del bloque.",
                        bloque.unsigned_abs()
                    ),
                ),
            );
            narracion.agregar(tablero_ciclo(cromosoma, &nodos, foco, &[]).paso(
                "let cabeza = nodo_cabeza(bloque);",
                format!("Su cabeza es 2·{} = {cabeza}.", bloque.unsigned_abs()),
            ));
            narracion.agregar(
                tablero_ciclo(cromosoma, &nodos, foco, &[]).paso(
                    "if bloque > 0 {",
                    if bloque > 0 {
                        "El signo es +: el bloque se lee de cola a cabeza, así que entra primero la cola."
                    } else {
                        "El signo es −: el bloque se lee al revés, de cabeza a cola, así que entra primero la cabeza."
                    },
                ),
            );
        }
        if bloque > 0 {
            nodos.push(cola);
            nodos.push(cabeza);
        } else {
            nodos.push(cabeza);
            nodos.push(cola);
        }
        let marcados = [anterior, anterior + 1];
        let linea = if bloque > 0 {
            "nodos.push(cola);"
        } else {
            "nodos.push(cabeza);"
        };
        match toque {
            Toque::Detalle => {
                narracion.agregar(tablero_ciclo(cromosoma, &nodos, foco, &marcados).paso(
                    linea,
                    format!(
                        "Entran los dos nodos del {bloque:+}: {} y {}. La lista ya tiene {}.",
                        nodos[anterior],
                        nodos[anterior + 1],
                        cuantos(nodos.len(), "nodo", "nodos")
                    ),
                ))
            }
            Toque::Resumen => {
                narracion.agregar(tablero_ciclo(cromosoma, &nodos, foco, &marcados).paso(
                    linea,
                    format!(
                        "Bloque {bloque:+} (posición {}): entran {} y {}.",
                        posicion + 1,
                        nodos[anterior],
                        nodos[anterior + 1]
                    ),
                ))
            }
            Toque::Nada => {}
        }
    }
    let ultimo = tablero_ciclo(cromosoma, &nodos, None, &[]).paso(
        "Ok(nodos_del_cromosoma(cromosoma))",
        format!(
            "Listo: {} bloques dieron {} nodos, el ciclo ({}). Cada bloque + usa su cola y luego su cabeza; cada bloque − al revés. Mira «Qué significa» para ver cómo se usa este ciclo.",
            miles(n),
            miles(nodos.len()),
            resumen(&nodos, nodo_escrito, " ")
        ),
    );
    debug_assert_eq!(nodos, esperado);
    Ok((narracion.terminar(ultimo), nodos))
}

fn pasos_cromosoma_a_ciclo(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let cromosoma = leer_permutacion(entrada, "cromosoma", MAXIMO_BLOQUES)?;
    Ok(narrar_cromosoma_a_ciclo(&cromosoma)?.0)
}

fn tablero_nodos(
    nodos: &[usize],
    cromosoma: &[i64],
    extra: &[(&str, String)],
    resaltar: &[usize],
) -> Pizarra {
    let mut variables = vec![("nodos", parentesis(resumen(nodos, nodo_escrito, " ")))];
    variables.push((extra[0].0, extra[0].1.clone()));
    variables.push(("cromosoma", parentesis(resumen(cromosoma, con_signo, " "))));
    for (nombre, valor) in &extra[1..] {
        variables.push((nombre, valor.clone()));
    }
    let centro = resaltar.first().copied().unwrap_or(0);
    con_lista(
        base_de(&variables),
        "nodos del ciclo",
        nodos,
        nodo_escrito,
        centro,
        resaltar,
    )
}

fn narrar_ciclo_a_cromosoma(nodos: &[usize]) -> Result<(Vec<Paso>, Vec<i64>), String> {
    let esperado = reordenamientos::ciclo_a_cromosoma(nodos)?;
    let cantidad = nodos.len();
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let cantidad_texto = vec![("cantidad", miles(cantidad))];
    narracion.agregar(
        tablero_nodos(nodos, &[], &cantidad_texto, &[]).paso(
            "let cantidad = nodos.len();",
            format!(
                "El ciclo tiene {}. Se leen de dos en dos: cada par es un bloque, y el orden dentro del par dice si lleva signo + o −.",
                cuantos(cantidad, "nodo", "nodos")
            ),
        ),
    );
    narracion.agregar(
        tablero_nodos(nodos, &[], &cantidad_texto, &[]).paso(
            "if cantidad == 0 || cantidad % 2 == 1 {",
            format!(
                "{} es par y no es cero: dos nodos por bloque, así que sale un número exacto de bloques.",
                miles(cantidad)
            ),
        ),
    );
    narracion.agregar(tablero_nodos(nodos, &[], &cantidad_texto, &[]).paso(
        "if let Some(fuera) = nodos.iter().find(|&&nodo| nodo == 0 || nodo > cantidad) {",
        format!(
            "Se busca algún nodo que no esté entre 1 y {}: no hay ninguno.",
            miles(cantidad)
        ),
    ));
    narracion.agregar(tablero_nodos(nodos, &[], &cantidad_texto, &[]).paso(
        "let mut cromosoma = Vec::with_capacity(cantidad / 2);",
        format!(
            "Preparamos el cromosoma, vacío: tendrá {}.",
            cuantos(cantidad / 2, "bloque", "bloques")
        ),
    ));
    let mut cromosoma: Vec<i64> = Vec::with_capacity(cantidad / 2);
    let mut muestreo = Muestreo::nuevo(cantidad / 2, 3, 60);
    for (posicion, extremos) in nodos.chunks(2).enumerate() {
        let entrada = extremos[0];
        let salida = extremos[1];
        let otro = reordenamientos::otro_extremo(entrada);
        let bloque = reordenamientos::bloque_del_nodo(entrada) as i64;
        let cola = reordenamientos::es_cola(entrada);
        let toque = muestreo.toca();
        let marcados = [2 * posicion, 2 * posicion + 1];
        let extras = |con_bloque: bool| {
            let mut lista = vec![
                ("cantidad", miles(cantidad)),
                ("extremos", format!("({entrada}, {salida})")),
                ("entrada", entrada.to_string()),
                ("salida", salida.to_string()),
            ];
            if con_bloque {
                lista.push(("bloque", bloque.to_string()));
            }
            lista
        };
        if matches!(toque, Toque::Detalle) {
            narracion.agregar(
                tablero_nodos(nodos, &cromosoma, &extras(false)[..2], &marcados).paso(
                    "for extremos in nodos.chunks(2) {",
                    format!(
                        "Tomamos el par número {}: los nodos ({entrada}, {salida}).",
                        posicion + 1
                    ),
                ),
            );
            narracion.agregar(
                tablero_nodos(nodos, &cromosoma, &extras(false), &marcados).paso(
                    "if salida != otro_extremo(entrada) {",
                    format!(
                        "El otro extremo del nodo {entrada} es {otro}, y la salida es {salida}: coinciden. Los dos extremos del bloque i son 2i−1 y 2i, así que este par sí es un bloque."
                    ),
                ),
            );
            narracion.agregar(
                tablero_nodos(nodos, &cromosoma, &extras(true), &marcados).paso(
                    "let bloque = bloque_del_nodo(entrada) as i64;",
                    format!(
                        "El nodo {entrada} pertenece al bloque {bloque}: es el entero hacia arriba de {entrada} ÷ 2."
                    ),
                ),
            );
            narracion.agregar(
                tablero_nodos(nodos, &cromosoma, &extras(true), &marcados).paso(
                    "if es_cola(entrada) {",
                    if cola {
                        format!("{entrada} es impar, o sea una cola: entramos por la cola, así que el bloque lleva signo +.")
                    } else {
                        format!("{entrada} es par, o sea una cabeza: entramos por la cabeza, así que el bloque va al revés y lleva signo −.")
                    },
                ),
            );
        }
        cromosoma.push(if cola { bloque } else { -bloque });
        let linea = if cola {
            "cromosoma.push(bloque);"
        } else {
            "cromosoma.push(-bloque);"
        };
        let nota = match toque {
            Toque::Detalle => Some(format!(
                "Se agrega el {:+} al cromosoma, que ya va en {}.",
                cromosoma[posicion],
                cuantos(cromosoma.len(), "bloque", "bloques")
            )),
            Toque::Resumen => Some(format!(
                "Par {} = ({entrada}, {salida}): bloque {:+}.",
                posicion + 1,
                cromosoma[posicion]
            )),
            Toque::Nada => None,
        };
        if let Some(nota) = nota {
            narracion.agregar(
                tablero_nodos(nodos, &cromosoma, &extras(true), &marcados).paso(linea, nota),
            );
        }
    }
    narracion.agregar(
        tablero_nodos(nodos, &cromosoma, &cantidad_texto, &[]).paso(
            "validar_bloques(&cromosoma)?;",
            format!(
                "Se comprueba que los bloques obtenidos sean los números del 1 al {}, sin repetir: lo son.",
                miles(cromosoma.len())
            ),
        ),
    );
    let ultimo = tablero_nodos(nodos, &cromosoma, &cantidad_texto, &[]).paso(
        "Ok(cromosoma)",
        format!(
            "Resultado: el cromosoma ({}). Los pares que iban de cabeza a cola salieron con signo −. «Qué significa» te cuenta cuántos bloques son negativos.",
            resumen(&cromosoma, con_signo, " ")
        ),
    );
    debug_assert_eq!(cromosoma, esperado);
    Ok((narracion.terminar(ultimo), cromosoma))
}

fn pasos_ciclo_a_cromosoma(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let nodos = leer_nodos(entrada, "ciclo")?;
    Ok(narrar_ciclo_a_cromosoma(&nodos)?.0)
}

fn tablero_aristas(genoma: &str, aristas: &[(usize, usize)], extra: &[(&str, String)]) -> Pizarra {
    let mut variables = vec![
        ("genoma", genoma.to_string()),
        ("aristas", format!("{} aristas", miles(aristas.len()))),
    ];
    for (nombre, valor) in extra {
        variables.push((nombre, valor.clone()));
    }
    base_de(&variables)
}

fn narrar_aristas_de_color(genoma: &[Vec<i64>]) -> NarradoDeAristas {
    let esperadas = reordenamientos::aristas_de_color(genoma)?;
    let total_bloques = contar_bloques(genoma);
    let genoma_texto = corto(&texto_genoma(genoma));
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        tablero_aristas(&genoma_texto, &[], &[]).paso(
            "validar_genoma(genoma)?;",
            format!(
                "Tu genoma es {genoma_texto}: {} con {}. Una arista de color une el final de un bloque con el comienzo del bloque que le sigue, y el último con el primero porque el cromosoma es un círculo.",
                cuantos(genoma.len(), "cromosoma", "cromosomas"),
                cuantos(total_bloques, "bloque", "bloques")
            ),
        ),
    );
    let mut aristas: Vec<(usize, usize)> = Vec::with_capacity(total_bloques);
    let mut por_cromosoma = Muestreo::nuevo(genoma.len(), 2, 15);
    let mut por_arista = Muestreo::nuevo(total_bloques, 4, 60);
    for cromosoma in genoma {
        let nodos = reordenamientos::nodos_del_cromosoma(cromosoma);
        let total = nodos.len();
        let cromosoma_texto = parentesis(resumen(cromosoma, con_signo, " "));
        let nodos_texto = parentesis(resumen(&nodos, nodo_escrito, " "));
        let toque = por_cromosoma.toca();
        let extra_cromosoma = vec![("cromosoma", cromosoma_texto.clone())];
        if !matches!(toque, Toque::Nada) {
            narracion.agregar(
                tablero_aristas(&genoma_texto, &aristas, &extra_cromosoma).paso(
                    "for cromosoma in genoma {",
                    format!(
                        "Tomamos el cromosoma {cromosoma_texto}, de {}.",
                        cuantos(cromosoma.len(), "bloque", "bloques")
                    ),
                ),
            );
            let con_nodos = vec![
                ("cromosoma", cromosoma_texto.clone()),
                ("nodos", nodos_texto.clone()),
            ];
            narracion.agregar(
                con_lista(
                    tablero_aristas(&genoma_texto, &aristas, &con_nodos),
                    "nodos del cromosoma",
                    &nodos,
                    nodo_escrito,
                    0,
                    &[],
                )
                .paso(
                    "let nodos = nodos_del_cromosoma(cromosoma);",
                    "Cada bloque pasa a ser dos nodos, cola y cabeza (un bloque − entra al revés). Con esos nodos ya podemos hablar de extremos.",
                ),
            );
            let con_total = vec![
                ("cromosoma", cromosoma_texto.clone()),
                ("nodos", nodos_texto.clone()),
                ("total", total.to_string()),
            ];
            narracion.agregar(
                tablero_aristas(&genoma_texto, &aristas, &con_total).paso(
                    "let total = nodos.len();",
                    format!(
                        "Hay {total} nodos. Lo usaremos para volver al primero al llegar al último, cerrando el círculo."
                    ),
                ),
            );
        }
        for j in 0..cromosoma.len() {
            let salida = nodos[2 * j + 1];
            let entrada_del_siguiente = nodos[(2 * j + 2) % total];
            let toque = por_arista.toca();
            let base = vec![
                ("cromosoma", cromosoma_texto.clone()),
                ("nodos", nodos_texto.clone()),
                ("total", total.to_string()),
                ("j", j.to_string()),
            ];
            if matches!(toque, Toque::Detalle) {
                let mut con_salida = base.clone();
                con_salida.push(("salida", salida.to_string()));
                narracion.agregar(
                    con_lista(
                        tablero_aristas(&genoma_texto, &aristas, &con_salida),
                        "nodos del cromosoma",
                        &nodos,
                        nodo_escrito,
                        2 * j,
                        &[2 * j + 1],
                    )
                    .paso(
                        "let salida = nodos[2 * j + 1];",
                        format!(
                            "El bloque número {} termina en el nodo {salida}: esa es su salida.",
                            j + 1
                        ),
                    ),
                );
                let mut con_entrada = con_salida.clone();
                con_entrada.push(("entrada_del_siguiente", entrada_del_siguiente.to_string()));
                narracion.agregar(
                    con_lista(
                        tablero_aristas(&genoma_texto, &aristas, &con_entrada),
                        "nodos del cromosoma",
                        &nodos,
                        nodo_escrito,
                        2 * j,
                        &[2 * j + 1, (2 * j + 2) % total],
                    )
                    .paso(
                        "let entrada_del_siguiente = nodos[(2 * j + 2) % total];",
                        if j + 1 == cromosoma.len() {
                            format!(
                                "El siguiente de este bloque es el primero, porque el cromosoma se cierra en círculo: su entrada es el nodo {entrada_del_siguiente}."
                            )
                        } else {
                            format!(
                                "El bloque que sigue empieza en el nodo {entrada_del_siguiente}: esa es su entrada."
                            )
                        },
                    ),
                );
            }
            aristas.push((salida, entrada_del_siguiente));
            if !matches!(toque, Toque::Nada) {
                let mut con_todo = base.clone();
                con_todo.push(("salida", salida.to_string()));
                con_todo.push(("entrada_del_siguiente", entrada_del_siguiente.to_string()));
                narracion.agregar(
                    con_lista(
                        tablero_aristas(&genoma_texto, &aristas, &con_todo),
                        "aristas de color hasta ahora",
                        &aristas,
                        arista_escrita,
                        aristas.len() - 1,
                        &[aristas.len() - 1],
                    )
                    .paso(
                        "aristas.push((salida, entrada_del_siguiente));",
                        format!(
                            "Nace la arista de color ({salida}, {entrada_del_siguiente}): une la salida de un bloque con la entrada del siguiente. Van {}.",
                            cuantos(aristas.len(), "arista", "aristas")
                        ),
                    ),
                );
            }
        }
    }
    let ultimo = con_lista(
        tablero_aristas(&genoma_texto, &aristas, &[]),
        "aristas de color",
        &aristas,
        arista_escrita,
        0,
        &[],
    )
    .paso(
        "Ok(aristas_entre_vecinos(genoma))",
        format!(
            "Resultado: {} de color, una por bloque: {}. Mira «Qué significa» para ver cómo se leen en el grafo del genoma.",
            miles(aristas.len()),
            corto(&texto_aristas(&aristas))
        ),
    );
    debug_assert_eq!(aristas, esperadas);
    Ok((narracion.terminar(ultimo), aristas))
}

fn pasos_aristas_de_color(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let genoma = leer_genoma(entrada, "genoma", MAXIMO_BLOQUES)?;
    Ok(narrar_aristas_de_color(&genoma)?.0)
}

fn tabla_de_parejas(
    pareja: &[usize],
    resaltar: &[usize],
) -> (Vec<Vec<String>>, Vec<(usize, usize)>) {
    let ultimo = (pareja.len() - 1).min(12);
    let filas: Vec<Vec<String>> = (1..=ultimo)
        .map(|nodo| vec![nodo.to_string(), pareja[nodo].to_string()])
        .collect();
    let marcas: Vec<(usize, usize)> = resaltar
        .iter()
        .filter(|&&nodo| nodo >= 1 && nodo <= ultimo)
        .map(|&nodo| (nodo - 1, 1))
        .collect();
    (filas, marcas)
}

struct EstadoRecorrido {
    pareja: Vec<usize>,
    visitado: Vec<bool>,
    cromosomas: usize,
    cromosoma: Vec<i64>,
    primero: usize,
    entrada: usize,
    bloque: usize,
    salida: usize,
}

fn tablero_recorrido(
    aristas_texto: &str,
    bloques: usize,
    estado: &EstadoRecorrido,
    hasta: usize,
    resaltar: &[usize],
) -> Pizarra {
    let mut variables = vec![
        ("aristas", aristas_texto.to_string()),
        ("bloques", bloques.to_string()),
        (
            "pareja",
            parentesis(resumen(&estado.pareja[1..], nodo_escrito, " ")),
        ),
    ];
    if hasta >= 1 {
        let vistos = estado.visitado.iter().filter(|&&visto| visto).count();
        variables.push(("visitado", format!("{vistos} de {bloques} bloques")));
        variables.push(("genoma", format!("{} cromosomas", estado.cromosomas)));
    }
    if hasta >= 2 {
        variables.push(("primero", estado.primero.to_string()));
        variables.push((
            "cromosoma",
            parentesis(resumen(&estado.cromosoma, con_signo, " ")),
        ));
    }
    if hasta >= 3 {
        variables.push(("entrada", estado.entrada.to_string()));
    }
    if hasta >= 4 {
        variables.push(("bloque", estado.bloque.to_string()));
    }
    if hasta >= 5 {
        variables.push(("salida", estado.salida.to_string()));
    }
    let centro = estado.cromosoma.len().saturating_sub(1);
    con_lista(
        base_de(&variables),
        "cromosoma en construcción",
        &estado.cromosoma,
        con_signo,
        centro,
        resaltar,
    )
}

fn narrar_grafo_a_genoma(aristas: &[(usize, usize)]) -> Result<(Vec<Paso>, Genoma), String> {
    let esperado = reordenamientos::grafo_a_genoma(aristas)?;
    let bloques = aristas.len();
    let nodos = 2 * bloques;
    let aristas_texto = corto(&texto_aristas(aristas));
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let pareja = reordenamientos::pareja_de_cada_nodo(aristas, bloques);
    let mut estado = EstadoRecorrido {
        pareja: vec![0; nodos + 1],
        visitado: vec![false; bloques + 1],
        cromosomas: 0,
        cromosoma: Vec::new(),
        primero: 0,
        entrada: 0,
        bloque: 0,
        salida: 0,
    };
    narracion.agregar(
        tablero_recorrido(&aristas_texto, bloques, &estado, 0, &[]).paso(
            "let bloques = aristas.len();",
            format!(
                "Tienes {}: {aristas_texto}. Cada arista de color es una unión entre dos nodos; hay un bloque por arista, así que son {} y {} nodos. Vamos a recorrer el grafo para leer los cromosomas.",
                cuantos(bloques, "arista de color", "aristas de color"),
                cuantos(bloques, "bloque", "bloques"),
                miles(nodos)
            ),
        ),
    );
    narracion.agregar(
        tablero_recorrido(&aristas_texto, bloques, &estado, 0, &[]).paso(
            "for &(a, b) in aristas {",
            format!(
                "Se revisa cada arista: ningún nodo puede ser el mismo en sus dos puntas, repetirse en otra arista ni quedar fuera de 1..{}. Todo en orden.",
                miles(nodos)
            ),
        ),
    );
    if let Some(&(a, b)) = aristas.first() {
        let mut parcial = estado.pareja.clone();
        parcial[a] = b;
        parcial[b] = a;
        let (filas, marcas) = tabla_de_parejas(&parcial, &[a, b]);
        let parcial_estado = EstadoRecorrido {
            pareja: parcial,
            visitado: estado.visitado.clone(),
            cromosomas: 0,
            cromosoma: Vec::new(),
            primero: 0,
            entrada: 0,
            bloque: 0,
            salida: 0,
        };
        narracion.agregar(
            tablero_recorrido(&aristas_texto, bloques, &parcial_estado, 0, &[])
                .tabla("pareja de cada nodo", &["nodo", "pareja"], filas, &marcas)
                .paso(
                    "pareja[a] = b;",
                    format!(
                        "Con la primera arista ({a}, {b}): la pareja del {a} es el {b} y la del {b} es el {a}. Así, desde cualquier nodo sabemos a dónde nos lleva su arista."
                    ),
                ),
        );
    }
    estado.pareja = pareja.clone();
    let (filas, marcas) = tabla_de_parejas(&estado.pareja, &[]);
    narracion.agregar(
        tablero_recorrido(&aristas_texto, bloques, &estado, 0, &[])
            .tabla("pareja de cada nodo", &["nodo", "pareja"], filas, &marcas)
            .paso(
                "let pareja = pareja_de_cada_nodo(aristas, bloques);",
                "Repetimos eso con todas las aristas: la tabla pareja dice, para cada nodo, con cuál está unido.",
            ),
    );
    narracion.agregar(
        tablero_recorrido(&aristas_texto, bloques, &estado, 1, &[]).paso(
            "let mut visitado = vec![false; bloques + 1];",
            "Anotamos qué bloques ya usamos, todos en falso. Cada bloque pertenece a un solo cromosoma.",
        ),
    );
    let mut genoma: Genoma = Vec::new();
    let mut por_cromosoma = Muestreo::nuevo(esperado.len(), 2, 15);
    let mut por_vuelta = Muestreo::nuevo(bloques, 5, 60);
    for primero in 1..=bloques {
        if estado.visitado[primero] {
            if genoma.len() <= 1 && bloques <= 30 {
                estado.primero = primero;
                narracion.agregar(
                    tablero_recorrido(&aristas_texto, bloques, &estado, 2, &[]).paso(
                        "if visitado[primero] {",
                        format!(
                            "El bloque {primero} ya fue usado en un cromosoma anterior: se salta."
                        ),
                    ),
                );
            }
            continue;
        }
        estado.primero = primero;
        estado.cromosoma = Vec::new();
        let toque_cromosoma = por_cromosoma.toca();
        let mostrar_cromosoma = !matches!(toque_cromosoma, Toque::Nada);
        if mostrar_cromosoma {
            narracion.agregar(
                tablero_recorrido(&aristas_texto, bloques, &estado, 2, &[]).paso(
                    "for primero in 1..=bloques {",
                    format!(
                        "El bloque {primero} no se ha usado: empieza un cromosoma nuevo, entrando por su cola, el nodo {}.",
                        2 * primero - 1
                    ),
                ),
            );
        }
        estado.entrada = 2 * primero - 1;
        if mostrar_cromosoma {
            narracion.agregar(
                tablero_recorrido(&aristas_texto, bloques, &estado, 3, &[]).paso(
                    "let mut entrada = 2 * primero - 1;",
                    format!("Entramos por el nodo {}.", estado.entrada),
                ),
            );
        }
        while !estado.visitado[reordenamientos::bloque_del_nodo(estado.entrada)] {
            let toque = por_vuelta.toca();
            let detalle = matches!(toque, Toque::Detalle) && mostrar_cromosoma;
            let avisar = detalle || matches!(toque, Toque::Resumen);
            estado.bloque = reordenamientos::bloque_del_nodo(estado.entrada);
            if detalle {
                narracion.agregar(
                    tablero_recorrido(&aristas_texto, bloques, &estado, 4, &[]).paso(
                        "let bloque = bloque_del_nodo(entrada);",
                        format!(
                            "El nodo {} pertenece al bloque {}.",
                            estado.entrada, estado.bloque
                        ),
                    ),
                );
            }
            estado.visitado[estado.bloque] = true;
            let entra_por_cola = reordenamientos::es_cola(estado.entrada);
            estado.cromosoma.push(if entra_por_cola {
                estado.bloque as i64
            } else {
                -(estado.bloque as i64)
            });
            let ultimo_indice = estado.cromosoma.len() - 1;
            if avisar {
                narracion.agregar(
                    tablero_recorrido(&aristas_texto, bloques, &estado, 4, &[ultimo_indice]).paso(
                        if entra_por_cola {
                            "cromosoma.push(bloque as i64);"
                        } else {
                            "cromosoma.push(-(bloque as i64));"
                        },
                        if entra_por_cola {
                            format!(
                                "Entramos por la cola ({}, impar), así que lo anotamos con signo +: +{}.",
                                estado.entrada, estado.bloque
                            )
                        } else {
                            format!(
                                "Entramos por la cabeza ({}, par), así que lo anotamos al revés: −{}.",
                                estado.entrada, estado.bloque
                            )
                        },
                    ),
                );
            }
            estado.salida = reordenamientos::otro_extremo(estado.entrada);
            let siguiente = estado.pareja[estado.salida];
            if detalle {
                narracion.agregar(
                    tablero_recorrido(&aristas_texto, bloques, &estado, 5, &[ultimo_indice]).paso(
                        "let salida = otro_extremo(entrada);",
                        format!(
                            "Cruzamos el bloque hasta su otro extremo: salimos por el nodo {}.",
                            estado.salida
                        ),
                    ),
                );
            }
            let anterior_entrada = estado.entrada;
            estado.entrada = siguiente;
            if avisar {
                narracion.agregar(
                    tablero_recorrido(&aristas_texto, bloques, &estado, 5, &[ultimo_indice]).paso(
                        "entrada = pareja[salida];",
                        format!(
                            "Desde el nodo {} saltamos por su arista de color hasta el nodo {siguiente}: esa es la nueva entrada (venimos de {anterior_entrada}).",
                            estado.salida
                        ),
                    ),
                );
            }
        }
        if mostrar_cromosoma {
            narracion.agregar(
                tablero_recorrido(&aristas_texto, bloques, &estado, 5, &[]).paso(
                    "while !visitado[bloque_del_nodo(entrada)] {",
                    format!(
                        "El nodo {} es del bloque {}, que ya usamos: el círculo se cerró.",
                        estado.entrada,
                        reordenamientos::bloque_del_nodo(estado.entrada)
                    ),
                ),
            );
        }
        genoma.push(std::mem::take(&mut estado.cromosoma));
        estado.cromosomas = genoma.len();
        if mostrar_cromosoma {
            let ultimo_cromosoma = genoma.last().cloned().unwrap_or_default();
            narracion.agregar(
                tablero_recorrido(&aristas_texto, bloques, &estado, 5, &[]).paso(
                    "genoma.push(cromosoma);",
                    format!(
                        "Cromosoma terminado: {}. Van {}.",
                        parentesis(resumen(&ultimo_cromosoma, con_signo, " ")),
                        cuantos(genoma.len(), "cromosoma", "cromosomas")
                    ),
                ),
            );
        }
    }
    let ultimo = tablero_recorrido(&aristas_texto, bloques, &estado, 1, &[]).paso(
        "Ok(recorrer_genoma(&pareja, bloques))",
        format!(
            "Resultado: {}, {} en total. Cada cromosoma se cerró al volver a un bloque ya visto. En «Qué significa» ves el grafo y por qué un círculo se puede escribir empezando en cualquier punto.",
            corto(&texto_genoma(&genoma)),
            cuantos(genoma.len(), "cromosoma", "cromosomas")
        ),
    );
    debug_assert_eq!(genoma, esperado);
    Ok((narracion.terminar(ultimo), genoma))
}

fn pasos_grafo_a_genoma(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let aristas = leer_aristas(entrada, "aristas")?;
    Ok(narrar_grafo_a_genoma(&aristas)?.0)
}

struct EstadoDistancia {
    rojas: Vec<usize>,
    azules: Vec<usize>,
    visitados: usize,
    ciclos: usize,
    inicio: usize,
    ciclo: Vec<usize>,
    nodo: usize,
    vecino_rojo: usize,
}

fn filas_rojo_azul(rojas: &[usize], azules: Option<&[usize]>) -> Vec<Vec<String>> {
    let ultimo = (rojas.len() - 1).min(12);
    (1..=ultimo)
        .map(|nodo| match azules {
            Some(azules) => vec![
                nodo.to_string(),
                rojas[nodo].to_string(),
                azules[nodo].to_string(),
            ],
            None => vec![nodo.to_string(), rojas[nodo].to_string()],
        })
        .collect()
}

fn tablero_distancia(
    textos: (&str, &str),
    bloques: Option<usize>,
    estado: &EstadoDistancia,
    hasta: usize,
) -> Pizarra {
    let mut variables = vec![("p", textos.0.to_string()), ("q", textos.1.to_string())];
    if let Some(bloques) = bloques {
        variables.push(("bloques", bloques.to_string()));
    }
    if hasta >= 1 {
        variables.push((
            "rojas",
            parentesis(resumen(&estado.rojas[1..], nodo_escrito, " ")),
        ));
    }
    if hasta >= 2 {
        variables.push((
            "azules",
            parentesis(resumen(&estado.azules[1..], nodo_escrito, " ")),
        ));
    }
    if hasta >= 3 {
        variables.push((
            "visitado",
            format!("{} de {} nodos", estado.visitados, estado.rojas.len() - 1),
        ));
        variables.push(("ciclos", format!("{} ciclos", estado.ciclos)));
    }
    if hasta >= 4 {
        variables.push(("inicio", estado.inicio.to_string()));
        variables.push((
            "ciclo",
            parentesis(resumen(&estado.ciclo, nodo_escrito, " ")),
        ));
    }
    if hasta >= 5 {
        variables.push(("nodo", estado.nodo.to_string()));
    }
    if hasta >= 6 {
        variables.push(("vecino_rojo", estado.vecino_rojo.to_string()));
    }
    let centro = estado.ciclo.len().saturating_sub(1);
    let pizarra = base_de(&variables);
    if hasta >= 4 {
        con_lista(
            pizarra,
            "ciclo en construcción",
            &estado.ciclo,
            nodo_escrito,
            centro,
            &[centro],
        )
    } else {
        pizarra
    }
}

fn narrar_distancia(p: &[Vec<i64>], q: &[Vec<i64>]) -> NarradoDeDistancia {
    let esperados = reordenamientos::ciclos_del_grafo_de_rupturas(p, q)?;
    let bloques = contar_bloques(p);
    let p_texto = corto(&texto_genoma(p));
    let q_texto = corto(&texto_genoma(q));
    let textos = (p_texto.as_str(), q_texto.as_str());
    let mut estado = EstadoDistancia {
        rojas: reordenamientos::pareja_de_cada_nodo(
            &reordenamientos::aristas_entre_vecinos(p),
            bloques,
        ),
        azules: reordenamientos::pareja_de_cada_nodo(
            &reordenamientos::aristas_entre_vecinos(q),
            bloques,
        ),
        visitados: 0,
        ciclos: 0,
        inicio: 0,
        ciclo: Vec::new(),
        nodo: 0,
        vecino_rojo: 0,
    };
    let mut visitado = vec![false; estado.rojas.len()];
    let mut ciclos: Vec<Vec<usize>> = Vec::new();
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        tablero_distancia(textos, None, &estado, 0).paso(
            "let ciclos = ciclos_del_grafo_de_rupturas(p, q)?;",
            format!(
                "Queremos la distancia entre P = {p_texto} y Q = {q_texto}. Para eso se arma un grafo con las uniones de P (rojas) y las de Q (azules) y se cuentan sus ciclos: la distancia es bloques menos ciclos."
            ),
        ),
    );
    narracion.agregar(tablero_distancia(textos, Some(bloques), &estado, 0).paso(
        "let bloques = validar_mismos_bloques(p, q)?;",
        format!(
            "Se comprueba que P y Q tengan los mismos bloques, sin repetir: los dos tienen {}.",
            cuantos(bloques, "bloque", "bloques")
        ),
    ));
    let filas_rojas = filas_rojo_azul(&estado.rojas, None);
    let marcas_rojas: Vec<(usize, usize)> = (0..filas_rojas.len()).map(|fila| (fila, 1)).collect();
    narracion.agregar(
        tablero_distancia(textos, Some(bloques), &estado, 1)
            .tabla(
                "pareja de cada nodo en P",
                &["nodo", "rojo (P)"],
                filas_rojas,
                &marcas_rojas,
            )
            .paso(
                "let rojas = pareja_de_cada_nodo(&aristas_entre_vecinos(p), bloques);",
                "Para cada nodo anotamos con cuál está unido en P: es su pareja roja. Las uniones de P salen de poner los bloques de cada cromosoma en círculo.",
            ),
    );
    let filas_azules = filas_rojo_azul(&estado.rojas, Some(&estado.azules));
    let marcas_azules: Vec<(usize, usize)> =
        (0..filas_azules.len()).map(|fila| (fila, 2)).collect();
    narracion.agregar(
        tablero_distancia(textos, Some(bloques), &estado, 2)
            .tabla(
                "pareja de cada nodo",
                &["nodo", "rojo (P)", "azul (Q)"],
                filas_azules,
                &marcas_azules,
            )
            .paso(
                "let azules = pareja_de_cada_nodo(&aristas_entre_vecinos(q), bloques);",
                "Lo mismo con Q: la pareja azul de cada nodo. Ahora cada nodo tiene una arista roja y una azul, así que el grafo se parte en ciclos que alternan los dos colores.",
            ),
    );
    narracion.agregar(
        tablero_distancia(textos, Some(bloques), &estado, 3).paso(
            "let mut visitado = vec![false; rojas.len()];",
            "Anotamos qué nodos ya pertenecen a un ciclo. Vamos a buscar ciclos empezando en el nodo 1, luego en el primero que no se haya usado, y así.",
        ),
    );
    let mut por_ciclo = Muestreo::nuevo(esperados.len(), 3, 25);
    let mut por_vuelta = Muestreo::nuevo(bloques, 6, 50);
    let mut saltos_mostrados = 0;
    for inicio in 1..estado.rojas.len() {
        if visitado[inicio] {
            if saltos_mostrados < 2 && ciclos.len() < 3 {
                saltos_mostrados += 1;
                estado.inicio = inicio;
                narracion.agregar(tablero_distancia(textos, Some(bloques), &estado, 3).paso(
                    "if visitado[inicio] {",
                    format!("El nodo {inicio} ya está en un ciclo anterior: se salta."),
                ));
            }
            continue;
        }
        let toque_ciclo = por_ciclo.toca();
        let mostrar_ciclo = !matches!(toque_ciclo, Toque::Nada);
        estado.inicio = inicio;
        estado.ciclo = Vec::new();
        estado.nodo = inicio;
        if mostrar_ciclo {
            narracion.agregar(tablero_distancia(textos, Some(bloques), &estado, 5).paso(
                "let mut nodo = inicio;",
                format!("El nodo {inicio} no está en ningún ciclo: empezamos uno nuevo desde ahí."),
            ));
        }
        loop {
            let toque = por_vuelta.toca();
            let detalle = matches!(toque, Toque::Detalle) && mostrar_ciclo;
            let avisar = detalle || matches!(toque, Toque::Resumen);
            estado.vecino_rojo = estado.rojas[estado.nodo];
            visitado[estado.nodo] = true;
            visitado[estado.vecino_rojo] = true;
            estado.visitados = visitado.iter().filter(|&&visto| visto).count();
            estado.ciclo.push(estado.nodo);
            estado.ciclo.push(estado.vecino_rojo);
            let desde = estado.nodo;
            if detalle {
                narracion.agregar(
                    tablero_distancia(textos, Some(bloques), &estado, 6).paso(
                        "let vecino_rojo = rojas[nodo];",
                        format!(
                            "Desde el nodo {desde} cruzamos su arista roja (de P) hasta el {}. Los dos nodos entran al ciclo.",
                            estado.vecino_rojo
                        ),
                    ),
                );
            }
            estado.nodo = estado.azules[estado.vecino_rojo];
            if avisar {
                narracion.agregar(tablero_distancia(textos, Some(bloques), &estado, 6).paso(
                    "nodo = azules[vecino_rojo];",
                    format!(
                        "Desde el {} cruzamos su arista azul (de Q) hasta el {}.",
                        estado.vecino_rojo, estado.nodo
                    ),
                ));
            }
            let cierra = estado.nodo == inicio;
            if detalle {
                narracion.agregar(
                    tablero_distancia(textos, Some(bloques), &estado, 6).paso(
                        "if nodo == inicio {",
                        if cierra {
                            format!("Llegamos al nodo {inicio}, donde empezamos: el ciclo se cierra.")
                        } else {
                            format!("Llegamos al {}, que no es el {inicio}: seguimos alternando rojo y azul.", estado.nodo)
                        },
                    ),
                );
            }
            if cierra {
                break;
            }
        }
        ciclos.push(std::mem::take(&mut estado.ciclo));
        estado.ciclos = ciclos.len();
        if mostrar_ciclo {
            let hecho = ciclos.last().cloned().unwrap_or_default();
            narracion.agregar(tablero_distancia(textos, Some(bloques), &estado, 6).paso(
                "ciclos.push(ciclo);",
                format!(
                    "Ciclo {} guardado: {}. {}",
                    ciclos.len(),
                    parentesis(resumen(&hecho, nodo_escrito, " ")),
                    if hecho.len() == 2 {
                        "Tiene solo 2 nodos: es un ciclo trivial, una unión que P y Q ya comparten."
                    } else {
                        "Tiene más de 2 nodos: P y Q difieren en estas uniones."
                    }
                ),
            ));
        }
    }
    let distancia = bloques - ciclos.len();
    narracion.agregar(tablero_distancia(textos, Some(bloques), &estado, 3).paso(
        "Ok(ciclos_rojo_azul(&rojas, &azules))",
        format!(
            "No quedan nodos sin visitar: en total hay {}.",
            cuantos(ciclos.len(), "ciclo", "ciclos")
        ),
    ));
    narracion.agregar(tablero_distancia(textos, Some(bloques), &estado, 3).paso(
        "let bloques: usize = p.iter().map(|cromosoma| cromosoma.len()).sum();",
        format!(
            "Contamos los bloques de P sumando los de cada cromosoma: {}.",
            miles(bloques)
        ),
    ));
    let ultimo = tablero_distancia(textos, Some(bloques), &estado, 3).paso(
        "Ok(bloques - ciclos.len())",
        format!(
            "Distancia = bloques − ciclos = {} − {} = {}: ese es el mínimo de 2-roturas para pasar de P a Q. Cada 2-rotura sube la cuenta de ciclos en 1 como máximo. «Transformación más corta» te muestra una ruta real.",
            miles(bloques),
            miles(ciclos.len()),
            miles(distancia)
        ),
    );
    debug_assert_eq!(ciclos, esperados);
    Ok((narracion.terminar(ultimo), distancia, ciclos))
}

fn pasos_distancia_dos_roturas(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let p = leer_genoma(entrada, "p", MAXIMO_BLOQUES)?;
    let q = leer_genoma(entrada, "q", MAXIMO_BLOQUES)?;
    Ok(narrar_distancia(&p, &q)?.0)
}

fn tablero_rotura(
    aristas: &[(usize, usize)],
    nodos: [usize; 4],
    extra: &[(&str, String)],
    resaltar: &[usize],
) -> Pizarra {
    let mut variables = vec![(
        "aristas",
        parentesis(resumen(aristas, arista_escrita, ", ")),
    )];
    variables.push((
        "nodos",
        format!(
            "i1 = {}, i2 = {}, i3 = {}, i4 = {}",
            nodos[0], nodos[1], nodos[2], nodos[3]
        ),
    ));
    for (nombre, valor) in extra {
        variables.push((nombre, valor.clone()));
    }
    let centro = resaltar.first().copied().unwrap_or(0);
    con_lista(
        base_de(&variables),
        "aristas de color",
        aristas,
        arista_escrita,
        centro,
        resaltar,
    )
}

struct DatosRotura<'a> {
    aristas: &'a [(usize, usize)],
    nodos: [usize; 4],
    primera: usize,
    segunda: usize,
}

fn pasos_de_la_rotura(narracion: &mut Narracion, datos: &DatosRotura) -> Vec<(usize, usize)> {
    let [i1, i2, i3, i4] = datos.nodos;
    let primera_texto = vec![("primera", datos.primera.to_string())];
    let ambas_texto = vec![
        ("primera", datos.primera.to_string()),
        ("segunda", datos.segunda.to_string()),
    ];
    narracion.agregar(
        tablero_rotura(
            datos.aristas,
            datos.nodos,
            &primera_texto,
            &[datos.primera],
        )
        .paso(
            "let primera = buscar_arista(aristas, i1, i2)",
            format!(
                "Buscamos la arista ({i1}, {i2}) en cualquier orden: está en la posición {} de la lista (en el código, primera = {}).",
                datos.primera + 1,
                datos.primera
            ),
        ),
    );
    narracion.agregar(
        tablero_rotura(
            datos.aristas,
            datos.nodos,
            &ambas_texto,
            &[datos.segunda, datos.primera],
        )
        .paso(
            "let segunda = buscar_arista(aristas, i3, i4)",
            format!(
                "Y la arista ({i3}, {i4}): está en la posición {} (segunda = {}).",
                datos.segunda + 1,
                datos.segunda
            ),
        ),
    );
    narracion.agregar(
        tablero_rotura(
            datos.aristas,
            datos.nodos,
            &ambas_texto,
            &[datos.primera, datos.segunda],
        )
        .paso(
            "if primera == segunda {",
            format!(
                "Las posiciones {} y {} son distintas: son dos aristas diferentes, así que se pueden romper.",
                datos.primera, datos.segunda
            ),
        ),
    );
    let mut nuevas = datos.aristas.to_vec();
    narracion.agregar(
        tablero_rotura(&nuevas, datos.nodos, &ambas_texto, &[datos.primera]).paso(
            "let mut nuevas = aristas.to_vec();",
            "Se trabaja sobre una copia de las aristas, para no tocar las originales.",
        ),
    );
    nuevas[datos.primera] = (i1, i3);
    narracion.agregar(
        tablero_rotura(&nuevas, datos.nodos, &ambas_texto, &[datos.primera])
            .tabla(
                "cambio",
                &["Se quitan", "Se agregan"],
                vec![vec![format!("({i1}, {i2})"), format!("({i1}, {i3})")]],
                &[(0, 1)],
            )
            .paso(
                "nuevas[primera] = (i1, i3);",
                format!(
                    "La arista ({i1}, {i2}) se rompe y en su lugar nace ({i1}, {i3}): el nodo {i1} se une ahora con el {i3}."
                ),
            ),
    );
    nuevas[datos.segunda] = (i2, i4);
    narracion.agregar(
        tablero_rotura(&nuevas, datos.nodos, &ambas_texto, &[datos.segunda, datos.primera])
            .tabla(
                "cambio",
                &["Se quitan", "Se agregan"],
                vec![
                    vec![format!("({i1}, {i2})"), format!("({i1}, {i3})")],
                    vec![format!("({i3}, {i4})"), format!("({i2}, {i4})")],
                ],
                &[(1, 1)],
            )
            .paso(
                "nuevas[segunda] = (i2, i4);",
                format!(
                    "La arista ({i3}, {i4}) se rompe y nace ({i2}, {i4}): los dos extremos que quedaron sueltos, {i2} y {i4}, se unen entre sí."
                ),
            ),
    );
    nuevas
}

fn narrar_dos_roturas_en_grafo(aristas: &[(usize, usize)], nodos: [usize; 4]) -> NarradoDeAristas {
    let [i1, i2, i3, i4] = nodos;
    let antes = reordenamientos::grafo_a_genoma(aristas)?;
    let esperadas = reordenamientos::dos_roturas_en_grafo(aristas, i1, i2, i3, i4)?;
    let despues = reordenamientos::grafo_a_genoma(&esperadas)?;
    let primera = reordenamientos::buscar_arista(aristas, i1, i2).unwrap_or(0);
    let segunda = reordenamientos::buscar_arista(aristas, i3, i4).unwrap_or(0);
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        tablero_rotura(aristas, nodos, &[], &[]).paso(
            "pub fn dos_roturas_en_grafo(",
            format!(
                "Tienes {} y quieres romper ({i1}, {i2}) y ({i3}, {i4}). Una 2-rotura quita esas dos aristas y pega sus cuatro extremos de otra forma: ({i1}, {i3}) y ({i2}, {i4}).",
                cuantos(aristas.len(), "arista de color", "aristas de color")
            ),
        ),
    );
    let datos = DatosRotura {
        aristas,
        nodos,
        primera,
        segunda,
    };
    let nuevas = pasos_de_la_rotura(&mut narracion, &datos);
    let ultimo = tablero_rotura(&nuevas, nodos, &[], &[primera, segunda]).paso(
        "Ok(nuevas)",
        format!(
            "Resultado: {}. El genoma pasó de {} a {}: {}. Los bloques no cambian, solo qué extremos están unidos; «Qué significa» explica el tipo de cambio.",
            corto(&texto_aristas(&nuevas)),
            cuantos(antes.len(), "cromosoma", "cromosomas"),
            cuantos(despues.len(), "cromosoma", "cromosomas"),
            efecto_de_la_rotura(antes.len(), despues.len())
        ),
    );
    debug_assert_eq!(nuevas, esperadas);
    Ok((narracion.terminar(ultimo), nuevas))
}

fn pasos_dos_roturas_en_grafo(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let aristas = leer_aristas(entrada, "aristas")?;
    let nodos = leer_cuatro_nodos(entrada, "nodos")?;
    Ok(narrar_dos_roturas_en_grafo(&aristas, nodos)?.0)
}

fn narrar_dos_roturas_en_genoma(
    genoma: &[Vec<i64>],
    nodos: [usize; 4],
) -> Result<(Vec<Paso>, Genoma), String> {
    let [i1, i2, i3, i4] = nodos;
    let esperado = reordenamientos::dos_roturas_en_genoma(genoma, i1, i2, i3, i4)?;
    let aristas = reordenamientos::aristas_entre_vecinos(genoma);
    let primera = reordenamientos::buscar_arista(&aristas, i1, i2).unwrap_or(0);
    let segunda = reordenamientos::buscar_arista(&aristas, i3, i4).unwrap_or(0);
    let genoma_texto = corto(&texto_genoma(genoma));
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let variables_iniciales = vec![
        ("genoma", genoma_texto.clone()),
        (
            "nodos",
            format!("i1 = {i1}, i2 = {i2}, i3 = {i3}, i4 = {i4}"),
        ),
    ];
    narracion.agregar(
        base_de(&variables_iniciales).paso(
            "validar_genoma(genoma)?;",
            format!(
                "Tu genoma es {genoma_texto}: {} con {}. Vas a hacer una 2-rotura en las uniones ({i1}, {i2}) y ({i3}, {i4}). Primero se comprueba que el genoma sea válido: bloques del 1 al {n} sin repetir.",
                cuantos(genoma.len(), "cromosoma", "cromosomas"),
                cuantos(contar_bloques(genoma), "bloque", "bloques"),
                n = contar_bloques(genoma)
            ),
        ),
    );
    let mut con_aristas = variables_iniciales.clone();
    con_aristas.push((
        "aristas",
        parentesis(resumen(&aristas, arista_escrita, ", ")),
    ));
    narracion.agregar(
        con_lista(
            base_de(&con_aristas),
            "aristas de color del genoma",
            &aristas,
            arista_escrita,
            primera,
            &[primera, segunda],
        )
        .paso(
            "let aristas = aristas_entre_vecinos(genoma);",
            format!(
                "El genoma se pasa a aristas de color: {}. Una 2-rotura trabaja sobre estas uniones.",
                corto(&texto_aristas(&aristas))
            ),
        ),
    );
    narracion.agregar(
        con_lista(
            base_de(&con_aristas),
            "aristas de color del genoma",
            &aristas,
            arista_escrita,
            primera,
            &[primera, segunda],
        )
        .paso(
            "let nuevas = dos_roturas_en_grafo(&aristas, i1, i2, i3, i4)?;",
            format!(
                "Ahora se rompen las aristas ({i1}, {i2}), en la posición {}, y ({i3}, {i4}), en la posición {}.",
                primera + 1,
                segunda + 1
            ),
        ),
    );
    let datos = DatosRotura {
        aristas: &aristas,
        nodos,
        primera,
        segunda,
    };
    let nuevas = pasos_de_la_rotura(&mut narracion, &datos);
    let antes = genoma.len();
    let despues = reordenamientos::grafo_a_genoma(&nuevas)?;
    let mut con_nuevas = variables_iniciales.clone();
    con_nuevas.push((
        "aristas",
        parentesis(resumen(&aristas, arista_escrita, ", ")),
    ));
    con_nuevas.push(("nuevas", parentesis(resumen(&nuevas, arista_escrita, ", "))));
    narracion.agregar(
        con_lista(
            base_de(&con_nuevas),
            "aristas de color nuevas",
            &nuevas,
            arista_escrita,
            primera,
            &[primera, segunda],
        )
        .paso(
            "Ok(recorrer_genoma(&pareja, bloques))",
            format!(
                "Con las aristas nuevas se recorre otra vez el grafo, alternando bloque y arista de color, hasta cerrar cada círculo: salen {}.",
                cuantos(despues.len(), "cromosoma", "cromosomas")
            ),
        ),
    );
    let mut con_resultado = con_nuevas.clone();
    con_resultado.push(("resultado", corto(&texto_genoma(&despues))));
    let ultimo = base_de(&con_resultado).paso(
        "grafo_a_genoma(&nuevas)",
        format!(
            "Resultado: {}. Pasó de {} a {}: {}. «Qué significa» explica por qué dos cortes en un cromosoma circular solo pueden dar una inversión o una fisión, y en dos cromosomas distintos, una fusión.",
            corto(&texto_genoma(&despues)),
            cuantos(antes, "cromosoma", "cromosomas"),
            cuantos(despues.len(), "cromosoma", "cromosomas"),
            efecto_de_la_rotura(antes, despues.len())
        ),
    );
    debug_assert_eq!(despues, esperado);
    Ok((narracion.terminar(ultimo), despues))
}

fn pasos_dos_roturas_en_genoma(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let genoma = leer_genoma(entrada, "genoma", MAXIMO_BLOQUES)?;
    let nodos = leer_cuatro_nodos(entrada, "nodos")?;
    Ok(narrar_dos_roturas_en_genoma(&genoma, nodos)?.0)
}

struct KmeroVisto {
    fin: usize,
    posicion: usize,
    codigo: u64,
    directo: bool,
    numero_directo: u64,
    numero_reverso: u64,
}

fn primeros_kmeros(secuencia: &[u8], k: usize, cuantos_quiero: usize) -> Vec<KmeroVisto> {
    let mut ventana = VentanaDeKmeros::nueva(k);
    let mut vistos = Vec::new();
    for (fin, &base) in secuencia.iter().enumerate() {
        if let Some((codigo, directo)) = ventana.agregar(base) {
            vistos.push(KmeroVisto {
                fin,
                posicion: fin + 1 - k,
                codigo,
                directo,
                numero_directo: ventana.directo,
                numero_reverso: ventana.reverso,
            });
            if vistos.len() == cuantos_quiero {
                break;
            }
        }
    }
    vistos
}

fn kmero_en(secuencia: &[u8], posicion: usize, k: usize) -> Option<KmeroVisto> {
    let mut ventana = VentanaDeKmeros::nueva(k);
    let mut resultado = None;
    for (desplazamiento, &base) in secuencia[posicion..posicion + k].iter().enumerate() {
        if let Some((codigo, directo)) = ventana.agregar(base) {
            resultado = Some(KmeroVisto {
                fin: posicion + desplazamiento,
                posicion,
                codigo,
                directo,
                numero_directo: ventana.directo,
                numero_reverso: ventana.reverso,
            });
        }
    }
    resultado
}

fn complemento_inverso(texto: &str) -> String {
    texto
        .chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' | 'U' => 'A',
            otra => otra,
        })
        .collect()
}

fn texto_kmero(secuencia: &[u8], posicion: usize, k: usize) -> String {
    String::from_utf8_lossy(&secuencia[posicion..posicion + k]).to_uppercase()
}

fn contar_kmeros(secuencia: &[u8], k: usize) -> usize {
    let mut ventana = VentanaDeKmeros::nueva(k);
    secuencia
        .iter()
        .filter(|&&base| ventana.agregar(base).is_some())
        .count()
}

fn forma_escrita(directo: bool) -> &'static str {
    if directo { "directa" } else { "reversa" }
}

struct DatosDeKmeros<'a> {
    a: &'a [u8],
    b: &'a [u8],
    k: usize,
    inicio_a: Option<usize>,
    inicio_b: Option<usize>,
}

fn marcas_de(posicion: usize, k: usize, marca: Marca) -> [(usize, usize, Marca); 1] {
    [(posicion, posicion + k, marca)]
}

fn narrar_kmeros(datos: &DatosDeKmeros) -> Result<(Vec<Paso>, Vec<ParCompartido>), String> {
    let (a, b, k) = (datos.a, datos.b, datos.k);
    let pares = reordenamientos::kmeros_compartidos(a, b, k)?;
    let indice = reordenamientos::indexar_kmeros(b, k);
    let kmeros_de_b = contar_kmeros(b, k);
    let kmeros_de_a = contar_kmeros(a, k);
    let invertidos = pares.iter().filter(|par| par.invertido).count();
    let directos = pares.len() - invertidos;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let largos = format!("{} pb en A, {} pb en B", miles(a.len()), miles(b.len()));
    let primer_kmero_a = a.len().min(k);
    narracion.agregar(
        base_de(&[("k", k.to_string())])
            .secuencia(a, &marcas_de(0, primer_kmero_a, Marca::Actual))
            .paso(
                "if k == 0 || k > 32 {",
                format!(
                    "Tienes dos secuencias ({largos}) y k = {k}. Buscamos todas las posiciones donde A y B comparten un k-mero de {k} bases, leído igual o como complemento inverso (la otra hebra). Se pide k entre 1 y 32 porque cada base ocupa 2 bits y el k-mero entero cabe en un número de 64 bits."
                ),
            ),
    );
    narracion.agregar(
        base_de(&[("k", k.to_string())])
            .secuencia(b, &marcas_de(0, k, Marca::Actual))
            .paso(
                "let indice = indexar_kmeros(b, k);",
                format!(
                    "Para no comparar cada k-mero de A con todos los de B, primero hacemos un índice de B: un diccionario que dice en qué posiciones aparece cada k-mero. Así, buscar uno cuesta casi nada. B tiene {} k-meros.",
                    miles(kmeros_de_b)
                ),
            ),
    );
    narrar_indice_de_b(&mut narracion, datos, &indice, kmeros_de_b);
    narracion.agregar(
        base_de(&[
            ("k", k.to_string()),
            (
                "indice",
                format!("{} k-meros distintos", miles(indice.primera_y_ultima.len())),
            ),
        ])
        .secuencia(a, &marcas_de(0, primer_kmero_a, Marca::Actual))
        .paso(
            "let mut ventana = VentanaDeKmeros::nueva(k);",
            "Con el índice listo, empieza la ventana para recorrer A: otra ventana vacía, base por base.",
        ),
    );
    narrar_recorrido_de_a(&mut narracion, datos, &indice, &pares);
    let con_pareja = {
        let mut posiciones: Vec<usize> = pares.iter().map(|par| par.posicion_a).collect();
        posiciones.dedup();
        posiciones.len()
    };
    narracion.agregar(
        base_de(&[
            ("k", k.to_string()),
            ("pares", format!("{} parejas", miles(pares.len()))),
        ])
        .paso(
            "for (fin, &base) in a.iter().enumerate() {",
            format!(
                "Seguimos así hasta el final de A: {} k-meros revisados, {} de ellos con al menos una pareja en B.",
                miles(kmeros_de_a),
                miles(con_pareja)
            ),
        ),
    );
    let ultimo = base_de(&[
        ("k", k.to_string()),
        ("pares", format!("{} parejas", miles(pares.len()))),
    ])
    .paso(
        "Ok(pares)",
        format!(
            "Resultado: {} parejas de {k}-meros compartidos, {} en la misma hebra y {} en la hebra contraria. En «Qué significa» comparamos esa cifra con lo que esperarías por puro azar.",
            miles(pares.len()),
            miles(directos),
            miles(invertidos)
        ),
    );
    Ok((narracion.terminar(ultimo), pares))
}

fn narrar_indice_de_b(
    narracion: &mut Narracion,
    datos: &DatosDeKmeros,
    indice: &reordenamientos::IndiceDeKmeros,
    kmeros_de_b: usize,
) {
    let (b, k) = (datos.b, datos.k);
    let mut ventana = VentanaDeKmeros::nueva(k);
    let primera = ventana.agregar(b[0]);
    if primera.is_none() {
        narracion.agregar(
            base_de(&[
                ("k", k.to_string()),
                ("fin", "0".to_string()),
                ("base", (b[0] as char).to_string()),
                ("ventana.directo", ventana.directo.to_string()),
                ("ventana.reverso", ventana.reverso.to_string()),
                ("ventana.bases_seguidas", ventana.bases_seguidas.to_string()),
            ])
            .secuencia(b, &marcas_de(0, 1, Marca::Actual))
            .paso(
                "let Some((codigo, directo)) = ventana.agregar(base) else {",
                format!(
                    "Entra la primera base de B, {}. La ventana lleva {} de las {k} que necesita, así que todavía no hay k-mero: devuelve nada y el ciclo sigue con la siguiente base.",
                    b[0] as char,
                    ventana.bases_seguidas
                ),
            ),
        );
    }
    let tres = primeros_kmeros(b, k, 3);
    let hasta = tres.last().map(|kmero| kmero.posicion).unwrap_or(0);
    for (orden, kmero) in tres.iter().enumerate() {
        let texto = texto_kmero(b, kmero.posicion, k);
        let reverso = complemento_inverso(&texto);
        let base_comun = vec![
            ("k", k.to_string()),
            ("fin", kmero.fin.to_string()),
            ("base", (b[kmero.fin] as char).to_string()),
            ("ventana.directo", kmero.numero_directo.to_string()),
            ("ventana.reverso", kmero.numero_reverso.to_string()),
        ];
        let marcas = marcas_de(kmero.posicion, k, Marca::Actual);
        narracion.agregar(
            base_de(&base_comun).secuencia(b, &marcas).paso(
                "let Some((codigo, directo)) = ventana.agregar(base) else {",
                if orden == 0 {
                    format!(
                        "Tras {k} bases la ventana ya tiene un k-mero completo: {texto}. En directo = {} está ese k-mero leído como número (cada base ocupa 2 bits: A=0, C=1, G=2, T=3). En reverso = {} está su complemento inverso, {reverso}, leído igual.",
                        kmero.numero_directo, kmero.numero_reverso
                    )
                } else {
                    format!(
                        "La ventana avanza una base: se descarta la más vieja y entra {}. Nuevo k-mero {texto}, directo = {}, reverso = {} ({reverso}).",
                        b[kmero.fin] as char,
                        kmero.numero_directo,
                        kmero.numero_reverso
                    )
                },
            ),
        );
        let mut con_codigo = base_comun.clone();
        con_codigo.push(("codigo", kmero.codigo.to_string()));
        con_codigo.push(("directo", kmero.directo.to_string()));
        narracion.agregar(
            base_de(&con_codigo).secuencia(b, &marcas).paso(
                "Some((self.directo.min(self.reverso), self.directo <= self.reverso))",
                format!(
                    "Nos quedamos con el menor de los dos números: codigo = {}. Así el mismo k-mero recibe el mismo código se lea por una hebra o por la otra. directo = {} nos dice que lo leímos en forma {}.",
                    kmero.codigo,
                    kmero.directo,
                    forma_escrita(kmero.directo)
                ),
            ),
        );
        let mut con_posicion = con_codigo.clone();
        con_posicion.push(("posicion", kmero.posicion.to_string()));
        narracion.agregar(
            base_de(&con_posicion).secuencia(b, &marcas).paso(
                "indice.directo[posicion] = directo;",
                format!(
                    "Guardamos en qué forma se leyó el k-mero de la posición {} (fin + 1 − k = {}), para saber luego si una pareja es de la misma hebra o de la contraria.",
                    kmero.posicion, kmero.posicion
                ),
            ),
        );
        let primera_vez = indice.primera_y_ultima[&kmero.codigo].0 == kmero.posicion;
        narracion.agregar(
            base_de(&con_posicion).secuencia(b, &marcas).paso(
                "match indice.primera_y_ultima.get_mut(&codigo) {",
                if primera_vez {
                    format!("Buscamos el código {} en el diccionario: todavía no está.", kmero.codigo)
                } else {
                    format!(
                        "Buscamos el código {} en el diccionario: ya estaba, este k-mero se repite en B.",
                        kmero.codigo
                    )
                },
            ),
        );
        if primera_vez {
            narracion.agregar(
                base_de(&con_posicion).secuencia(b, &marcas).paso(
                    "indice.primera_y_ultima.insert(codigo, (posicion, posicion));",
                    format!(
                        "Se agrega al diccionario: el código {} aparece por primera y última vez (hasta ahora) en la posición {}.",
                        kmero.codigo, kmero.posicion
                    ),
                ),
            );
        }
    }
    let repetido = indice
        .siguiente
        .iter()
        .enumerate()
        .filter(|&(_, &despues)| despues != reordenamientos::NINGUNA)
        .map(|(antes, &despues)| (despues, antes))
        .min();
    if let Some((despues, antes)) = repetido
        && despues > hasta
        && let Some(kmero) = kmero_en(b, despues, k)
    {
        let texto = texto_kmero(b, despues, k);
        let marcas = [
            (despues, despues + k, Marca::Actual),
            (antes, antes + k, Marca::Comparado),
        ];
        let variables = vec![
            ("k", k.to_string()),
            ("codigo", kmero.codigo.to_string()),
            ("directo", kmero.directo.to_string()),
            ("posicion", despues.to_string()),
        ];
        narracion.agregar(
            base_de(&variables).secuencia(b, &marcas).paso(
                "match indice.primera_y_ultima.get_mut(&codigo) {",
                format!(
                    "Más adelante, en la posición {despues}, aparece {texto} con el código {}, que ya estaba en el diccionario (primera aparición en la posición {antes}, marcada en amarillo).",
                    kmero.codigo
                ),
            ),
        );
        narracion.agregar(
            base_de(&variables).secuencia(b, &marcas).paso(
                "indice.siguiente[*ultima] = posicion;",
                format!(
                    "En vez de guardar otro diccionario, se encadena: la aparición de la posición {antes} apunta a la de la posición {despues}. Así cada k-mero tiene una lista enlazada con todas sus posiciones."
                ),
            ),
        );
    }
    narracion.agregar(
        base_de(&[
            ("k", k.to_string()),
            (
                "indice",
                format!("{} k-meros distintos", miles(indice.primera_y_ultima.len())),
            ),
        ])
        .secuencia(b, &marcas_de(0, k, Marca::Actual))
        .paso(
            "let Some((codigo, directo)) = ventana.agregar(base) else {",
            format!(
                "Así se sigue hasta el final de B: {} k-meros en total, {} de ellos distintos y {} repetidos. El índice está listo.",
                miles(kmeros_de_b),
                miles(indice.primera_y_ultima.len()),
                miles(kmeros_de_b - indice.primera_y_ultima.len())
            ),
        ),
    );
}

fn narrar_recorrido_de_a(
    narracion: &mut Narracion,
    datos: &DatosDeKmeros,
    indice: &reordenamientos::IndiceDeKmeros,
    pares: &[ParCompartido],
) {
    let (a, k) = (datos.a, datos.k);
    let resumen_indice = format!("{} k-meros distintos", miles(indice.primera_y_ultima.len()));
    for kmero in primeros_kmeros(a, k, 3) {
        let texto = texto_kmero(a, kmero.posicion, k);
        let marcas = marcas_de(kmero.posicion, k, Marca::Actual);
        let variables = vec![
            ("k", k.to_string()),
            ("indice", resumen_indice.clone()),
            ("fin", kmero.fin.to_string()),
            ("base", (a[kmero.fin] as char).to_string()),
            ("codigo", kmero.codigo.to_string()),
            ("directo_a", kmero.directo.to_string()),
        ];
        narracion.agregar(
            base_de(&variables).secuencia(a, &marcas).paso(
                "let Some((codigo, directo_a)) = ventana.agregar(base) else {",
                format!(
                    "El k-mero de A en la posición {} es {texto}; su código es {} y se leyó en forma {}.",
                    kmero.posicion,
                    kmero.codigo,
                    forma_escrita(kmero.directo)
                ),
            ),
        );
        let hallado = indice.primera_y_ultima.get(&kmero.codigo);
        let mut con_primera = variables.clone();
        if let Some(&(primera, _)) = hallado {
            con_primera.push(("primera", primera.to_string()));
        }
        narracion.agregar(
            base_de(&con_primera).secuencia(a, &marcas).paso(
                "let Some(&(primera, _)) = indice.primera_y_ultima.get(&codigo) else {",
                match hallado {
                    Some(&(primera, _)) => format!(
                        "Ese código sí está en el índice de B: su primera aparición es la posición {primera}. Hay pareja."
                    ),
                    None => "Ese código no está en el índice de B: el k-mero no existe en B ni en su otra hebra. No hay pareja y se pasa a la siguiente base.".to_string(),
                },
            ),
        );
    }
    let mut elegidos: Vec<usize> = Vec::new();
    if let Some(primero) = pares.first() {
        elegidos.push(primero.posicion_a);
    }
    if let Some(contrario) = pares.iter().find(|par| par.invertido) {
        elegidos.push(contrario.posicion_a);
    }
    let mut anterior: Option<usize> = None;
    for par in pares {
        if anterior == Some(par.posicion_a) {
            elegidos.push(par.posicion_a);
            break;
        }
        anterior = Some(par.posicion_a);
    }
    elegidos.dedup();
    let mut vistos: Vec<usize> = Vec::new();
    for posicion_a in elegidos {
        if vistos.contains(&posicion_a) {
            continue;
        }
        vistos.push(posicion_a);
        narrar_pareja_de_a(narracion, datos, indice, pares, posicion_a);
    }
}

fn narrar_pareja_de_a(
    narracion: &mut Narracion,
    datos: &DatosDeKmeros,
    indice: &reordenamientos::IndiceDeKmeros,
    pares: &[ParCompartido],
    posicion_a: usize,
) {
    let (a, b, k) = (datos.a, datos.b, datos.k);
    let desde = pares.partition_point(|par| par.posicion_a < posicion_a);
    let grupo: Vec<&ParCompartido> = pares[desde..]
        .iter()
        .take_while(|par| par.posicion_a == posicion_a)
        .collect();
    let Some(kmero) = kmero_en(a, posicion_a, k) else {
        return;
    };
    let texto = texto_kmero(a, posicion_a, k);
    let resumen_indice = format!("{} k-meros distintos", miles(indice.primera_y_ultima.len()));
    let primera = grupo[0].posicion_b;
    let mostrado_a = posicion_mostrada(posicion_a, datos.inicio_a);
    let marcas_a = marcas_de(posicion_a, k, Marca::Actual);
    let variables = vec![
        ("k", k.to_string()),
        ("indice", resumen_indice.clone()),
        ("fin", kmero.fin.to_string()),
        ("codigo", kmero.codigo.to_string()),
        ("directo_a", kmero.directo.to_string()),
        ("primera", primera.to_string()),
        ("posicion_a", posicion_a.to_string()),
    ];
    narracion.agregar(
        base_de(&variables).secuencia(a, &marcas_a).paso(
            "let posicion_a = fin + 1 - k;",
            format!(
                "Veamos una pareja real. El k-mero {texto} de A termina en la base {}, así que empieza en la posición {mostrado_a}.",
                kmero.fin
            ),
        ),
    );
    let mut con_b = variables.clone();
    con_b.push(("posicion_b", primera.to_string()));
    narracion.agregar(
        base_de(&con_b)
            .secuencia(b, &marcas_de(primera, k, Marca::Hallado))
            .paso(
                "let mut posicion_b = primera;",
                format!(
                    "En B ese código aparece por primera vez en la posición {}. Empezamos a recorrer la cadena de apariciones desde ahí.",
                    posicion_mostrada(primera, datos.inicio_b)
                ),
            ),
    );
    for (orden, par) in grupo.iter().take(2).enumerate() {
        let posicion_b = par.posicion_b;
        let marcas_b = marcas_de(posicion_b, k, Marca::Hallado);
        let mut variables_b = variables.clone();
        variables_b.push(("posicion_b", posicion_b.to_string()));
        let directo_b = indice.directo[posicion_b];
        let texto_b = texto_kmero(b, posicion_b, k);
        narracion.agregar(base_de(&variables_b).secuencia(b, &marcas_b).paso(
            "invertido: directo_a != indice.directo[posicion_b],",
            format!(
                "En A se leyó en forma {} y en B en forma {}: {}. En B el texto es {texto_b}{}.",
                forma_escrita(kmero.directo),
                forma_escrita(directo_b),
                if par.invertido {
                    "son distintas, así que B lo tiene en la hebra contraria"
                } else {
                    "son iguales, así que están en la misma hebra"
                },
                if par.invertido {
                    format!(" (el complemento inverso de {texto})")
                } else {
                    String::new()
                }
            ),
        ));
        let numero = desde + orden + 1;
        let mut con_pares = variables_b.clone();
        con_pares.push(("pares", format!("{} parejas", miles(numero))));
        narracion.agregar(
            base_de(&con_pares).secuencia(b, &marcas_b).paso(
                "pares.push(ParCompartido {",
                format!(
                    "Se guarda la pareja: A {mostrado_a}, B {}, hebra {}. Es la pareja número {} de la lista.",
                    posicion_mostrada(posicion_b, datos.inicio_b),
                    if par.invertido { "contraria" } else { "misma" },
                    miles(numero)
                ),
            ),
        );
        if let Some(siguiente) = grupo.get(orden + 1) {
            narracion.agregar(
                base_de(&con_pares)
                    .secuencia(b, &marcas_de(siguiente.posicion_b, k, Marca::Hallado))
                    .paso(
                        "posicion_b = indice.siguiente[posicion_b];",
                        format!(
                            "Este k-mero aparece otra vez en B: la cadena nos lleva a la posición {}. Se repite el proceso y sale otra pareja para el mismo k-mero de A.",
                            posicion_mostrada(siguiente.posicion_b, datos.inicio_b)
                        ),
                    ),
            );
        } else {
            narracion.agregar(
                base_de(&con_pares).secuencia(b, &marcas_b).paso(
                    "while posicion_b != NINGUNA {",
                    "La cadena terminó (siguiente = NINGUNA): no hay más apariciones de este k-mero en B, y se sigue con la próxima base de A.",
                ),
            );
        }
    }
}

fn pasos_kmeros_compartidos(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (a, b, k) = leer_secuencias_y_k(entrada)?;
    let datos = DatosDeKmeros {
        a: &a,
        b: &b,
        k,
        inicio_a: inicio_en_ecoli(entrada, "a"),
        inicio_b: inicio_en_ecoli(entrada, "b"),
    };
    Ok(narrar_kmeros(&datos)?.0)
}

fn escribir_bloque_fila(posicion: usize, bloque: &BloqueSintenico) -> Vec<String> {
    vec![
        posicion.to_string(),
        miles(bloque.inicio_a),
        miles(bloque.fin_a),
        miles(bloque.inicio_b),
        miles(bloque.fin_b),
        if bloque.invertido { "−" } else { "+" }.to_string(),
        miles(bloque.kmeros),
    ]
}

fn tabla_de_bloques(pizarra: Pizarra, guardados: &[BloqueSintenico], ultimo: bool) -> Pizarra {
    let desde = guardados.len().saturating_sub(8);
    let filas: Vec<Vec<String>> = guardados[desde..]
        .iter()
        .enumerate()
        .map(|(indice, bloque)| escribir_bloque_fila(desde + indice + 1, bloque))
        .collect();
    let resaltar: Vec<(usize, usize)> = if ultimo && !filas.is_empty() {
        vec![(filas.len() - 1, 6)]
    } else {
        Vec::new()
    };
    pizarra.tabla(
        "bloques que se van guardando",
        &[
            "#", "inicio A", "fin A", "inicio B", "fin B", "sentido", "k-meros",
        ],
        filas,
        &resaltar,
    )
}

struct DatosDeSintenia<'a> {
    a: &'a [u8],
    b: &'a [u8],
    k: usize,
    distancia: usize,
    minimo: usize,
    inicio_a: Option<usize>,
    inicio_b: Option<usize>,
}

fn casilla_escrita(casilla: &Casilla) -> String {
    format!(
        "({}, {}, {})",
        casilla.0,
        casilla.1,
        if casilla.2 { "contraria" } else { "misma" }
    )
}

fn narrar_casillas(
    narracion: &mut Narracion,
    pares: &[ParCompartido],
    datos: &DatosDeSintenia,
    casillas: &HashMap<Casilla, Vec<usize>>,
) {
    let distancia = datos.distancia;
    let variables = vec![
        ("k", datos.k.to_string()),
        ("distancia", distancia.to_string()),
        ("minimo", datos.minimo.to_string()),
        ("pares", format!("{} parejas", miles(pares.len()))),
    ];
    narracion.agregar(
        base_de(&variables).paso(
            "let mut casillas: HashMap<Casilla, Vec<usize>> = HashMap::new();",
            format!(
                "Para agrupar, el plano de A contra B se divide en casillas cuadradas de {} pb de lado. Cada pareja caerá en una casilla según su posición en A, su posición en B y su sentido. Ahora el diccionario de casillas está vacío.",
                miles(distancia)
            ),
        ),
    );
    let mut parciales: Vec<(Casilla, usize)> = Vec::new();
    for (numero, par) in pares.iter().take(3).enumerate() {
        let casilla = reordenamientos::casilla_de(par, distancia);
        match parciales.iter_mut().find(|(otra, _)| *otra == casilla) {
            Some((_, cantidad)) => *cantidad += 1,
            None => parciales.push((casilla, 1)),
        }
        let filas: Vec<Vec<String>> = parciales
            .iter()
            .map(|(casilla, cantidad)| vec![casilla_escrita(casilla), cantidad.to_string()])
            .collect();
        let resaltar: Vec<(usize, usize)> = parciales
            .iter()
            .enumerate()
            .filter(|(_, (otra, _))| *otra == casilla)
            .map(|(fila, _)| (fila, 1))
            .collect();
        let mut con_par = variables.clone();
        con_par.push(("numero", numero.to_string()));
        narracion.agregar(
            base_de(&con_par)
                .tabla(
                    "casillas llenas hasta ahora (fila, columna, sentido)",
                    &["casilla", "parejas"],
                    filas,
                    &resaltar,
                )
                .paso(
                    "par.posicion_a / distancia,",
                    format!(
                        "La pareja {} está en A {}, en B {}, {}. Dividiendo entre {}: casilla fila {}, columna {}.",
                        numero + 1,
                        posicion_mostrada(par.posicion_a, datos.inicio_a),
                        posicion_mostrada(par.posicion_b, datos.inicio_b),
                        if par.invertido {
                            "hebra contraria"
                        } else {
                            "misma hebra"
                        },
                        miles(distancia),
                        casilla.0,
                        casilla.1
                    ),
                ),
        );
    }
    let mut llenas: Vec<(&Casilla, usize)> = casillas
        .iter()
        .map(|(casilla, numeros)| (casilla, numeros.len()))
        .collect();
    llenas.sort_by(|x, y| y.1.cmp(&x.1).then(x.0.cmp(y.0)));
    let mayor = llenas.first().map(|(_, cantidad)| *cantidad).unwrap_or(0);
    let filas: Vec<Vec<String>> = llenas
        .iter()
        .take(8)
        .map(|(casilla, cantidad)| vec![casilla_escrita(casilla), miles(*cantidad)])
        .collect();
    let mut con_casillas = variables.clone();
    con_casillas.push((
        "casillas",
        format!("{} casillas llenas", miles(casillas.len())),
    ));
    narracion.agregar(
        base_de(&con_casillas)
            .tabla(
                "casillas más llenas (fila, columna, sentido)",
                &["casilla", "parejas"],
                filas,
                &[(0, 1)],
            )
            .paso(
                ".push(numero);",
                format!(
                    "Así se reparten las {} parejas: {} casillas con al menos una pareja, y la más llena tiene {}. Las casillas con muchas parejas son zonas donde A y B se parecen.",
                    miles(pares.len()),
                    miles(casillas.len()),
                    cuantos(mayor, "pareja", "parejas")
                ),
            ),
    );
}

fn resumen_de_bloque(bloque: &BloqueSintenico) -> String {
    if bloque.kmeros == 0 {
        "vacío".to_string()
    } else {
        format!(
            "{} k-meros, A {}–{}, B {}–{}",
            miles(bloque.kmeros),
            miles(bloque.inicio_a),
            miles(bloque.fin_a),
            miles(bloque.inicio_b),
            miles(bloque.fin_b)
        )
    }
}

struct ContextoDeBloques<'a> {
    datos: &'a DatosDeSintenia<'a>,
    pares: &'a [ParCompartido],
    casillas: &'a HashMap<Casilla, Vec<usize>>,
}

fn variables_de_bloques(
    contexto: &ContextoDeBloques,
    visitadas: usize,
    guardados: usize,
) -> Vec<(&'static str, String)> {
    vec![
        ("k", contexto.datos.k.to_string()),
        ("distancia", contexto.datos.distancia.to_string()),
        ("minimo", contexto.datos.minimo.to_string()),
        ("pares", format!("{} parejas", miles(contexto.pares.len()))),
        (
            "casillas",
            format!("{} casillas llenas", miles(contexto.casillas.len())),
        ),
        ("visitadas", format!("{} casillas", miles(visitadas))),
        ("bloques", format!("{} guardados", miles(guardados))),
    ]
}

fn narrar_primer_bloque(
    narracion: &mut Narracion,
    contexto: &ContextoDeBloques,
    inicio: Casilla,
) -> BloqueSintenico {
    let (pares, casillas, k) = (contexto.pares, contexto.casillas, contexto.datos.k);
    let mut visitadas: HashSet<Casilla> = HashSet::new();
    visitadas.insert(inicio);
    let mut variables = variables_de_bloques(contexto, 1, 0);
    variables.push(("inicio", casilla_escrita(&inicio)));
    narracion.agregar(
        base_de(&variables).paso(
            "let mut pendientes = vec![inicio];",
            format!(
                "Empieza a crecer un bloque desde la casilla {}. Su lista de pendientes tiene solo esa casilla; iremos sacando casillas y añadiendo sus vecinas con parejas.",
                casilla_escrita(&inicio)
            ),
        ),
    );
    let mut bloque = BloqueSintenico {
        inicio_a: usize::MAX,
        fin_a: 0,
        inicio_b: usize::MAX,
        fin_b: 0,
        invertido: inicio.2,
        kmeros: 0,
    };
    let mut pendientes = vec![inicio];
    let mut sacadas = 0;
    while let Some((fila, columna, invertido)) = pendientes.pop() {
        sacadas += 1;
        let detalle = sacadas <= 3;
        let actual = (fila, columna, invertido);
        let mut variables = variables_de_bloques(contexto, visitadas.len(), 0);
        variables.push(("inicio", casilla_escrita(&inicio)));
        variables.push(("pendientes", format!("{} casillas", pendientes.len() + 1)));
        variables.push(("fila", fila.to_string()));
        variables.push(("columna", columna.to_string()));
        if detalle {
            narracion.agregar(
                base_de(&variables).paso(
                    "while let Some((fila, columna, invertido)) = pendientes.pop() {",
                    format!(
                        "Sacamos la casilla {} de los pendientes. Todas las parejas que tiene dentro se van a sumar al bloque.",
                        casilla_escrita(&actual)
                    ),
                ),
            );
        }
        for &numero in &casillas[&actual] {
            reordenamientos::sumar_al_bloque(&mut bloque, &pares[numero], k);
        }
        if detalle {
            let mut con_bloque = variables.clone();
            con_bloque.push(("bloque", resumen_de_bloque(&bloque)));
            narracion.agregar(
                base_de(&con_bloque).paso(
                    "sumar_al_bloque(&mut bloque, &pares[numero], k);",
                    format!(
                        "Esta casilla trae {}: el bloque ahora cubre {}. Los extremos se ajustan al mínimo y máximo de las posiciones, y cada pareja suma un k-mero.",
                        cuantos(casillas[&actual].len(), "pareja", "parejas"),
                        resumen_de_bloque(&bloque)
                    ),
                ),
            );
        }
        let mut nuevas = Vec::new();
        for otra_fila in fila.saturating_sub(1)..=fila + 1 {
            for otra_columna in columna.saturating_sub(1)..=columna + 1 {
                let vecina = (otra_fila, otra_columna, invertido);
                if casillas.contains_key(&vecina) && visitadas.insert(vecina) {
                    pendientes.push(vecina);
                    nuevas.push(vecina);
                }
            }
        }
        if detalle {
            let mut con_vecinas = variables.clone();
            con_vecinas.push(("pendientes", format!("{} casillas", pendientes.len())));
            narracion.agregar(
                base_de(&con_vecinas).paso(
                    "for otra_fila in fila.saturating_sub(1)..=fila + 1 {",
                    if nuevas.is_empty() {
                        "Se miran las 8 casillas vecinas y la propia (mismo sentido): ninguna tiene parejas sin visitar.".to_string()
                    } else {
                        format!(
                            "Se miran las 8 casillas vecinas: {} con parejas que no habíamos visitado ({}). Pasan a los pendientes.",
                            nuevas.len(),
                            nuevas
                                .iter()
                                .take(3)
                                .map(casilla_escrita)
                                .collect::<Vec<_>>()
                                .join(", ")
                        )
                    },
                ),
            );
        }
    }
    let mut final_variables = variables_de_bloques(contexto, visitadas.len(), 0);
    final_variables.push(("inicio", casilla_escrita(&inicio)));
    final_variables.push(("pendientes", "0 casillas".to_string()));
    final_variables.push(("bloque", resumen_de_bloque(&bloque)));
    narracion.agregar(
        base_de(&final_variables).paso(
            "while let Some((fila, columna, invertido)) = pendientes.pop() {",
            format!(
                "Se repitió con {} y la lista de pendientes se vació: el bloque está completo y cubre {}.",
                cuantos(sacadas, "casilla", "casillas"),
                resumen_de_bloque(&bloque)
            ),
        ),
    );
    bloque
}

fn narrar_sintenia(datos: &DatosDeSintenia) -> Result<(Vec<Paso>, Vec<BloqueSintenico>), String> {
    let sintenia = reordenamientos::bloques_de_sintenia(
        datos.a,
        datos.b,
        datos.k,
        datos.distancia,
        datos.minimo,
    )?;
    let pares = sintenia.pares.as_slice();
    let (distancia, minimo, k) = (datos.distancia, datos.minimo, datos.k);
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        base_de(&[
            ("k", k.to_string()),
            ("distancia", distancia.to_string()),
            ("minimo", minimo.to_string()),
        ])
        .paso(
            "if distancia == 0 {",
            format!(
                "Vamos a agrupar los k-meros compartidos entre A ({} pb) y B ({} pb) en bloques de sintenia. Se agrupan los que estén a menos de unas {} pb entre sí, y se descartan los bloques con menos de {} k-meros. La distancia es {}: mayor que cero, así que se puede dividir.",
                miles(datos.a.len()),
                miles(datos.b.len()),
                miles(distancia),
                miles(minimo),
                miles(distancia)
            ),
        ),
    );
    let con_pares = vec![
        ("k", k.to_string()),
        ("distancia", distancia.to_string()),
        ("minimo", minimo.to_string()),
        ("pares", format!("{} parejas", miles(pares.len()))),
    ];
    let invertidos = pares.iter().filter(|par| par.invertido).count();
    narracion.agregar(
        base_de(&con_pares).paso(
            "let pares = kmeros_compartidos(a, b, k)?;",
            format!(
                "Primero se buscan los k-meros compartidos, como en «k-meros compartidos»: salen {}, {} en la misma hebra y {} en la contraria.",
                cuantos(pares.len(), "pareja", "parejas"),
                miles(pares.len() - invertidos),
                miles(invertidos)
            ),
        ),
    );
    narracion.agregar(
        base_de(&con_pares).paso(
            "let bloques = agrupar_en_bloques(&pares, k, distancia, minimo);",
            if pares.is_empty() {
                "No hay ninguna pareja, así que no hay nada que agrupar: no habrá bloques.".to_string()
            } else {
                "Ahora se agrupan: parejas cercanas en el gráfico de puntos (posición en A contra posición en B) forman un bloque.".to_string()
            },
        ),
    );
    let mut casillas: HashMap<Casilla, Vec<usize>> = HashMap::new();
    for (numero, par) in pares.iter().enumerate() {
        casillas
            .entry(reordenamientos::casilla_de(par, distancia))
            .or_default()
            .push(numero);
    }
    let mut candidatos: Vec<(Casilla, BloqueSintenico, usize)> = Vec::new();
    let mut visitadas: HashSet<Casilla> = HashSet::new();
    for par in pares {
        let inicio = reordenamientos::casilla_de(par, distancia);
        if visitadas.insert(inicio) {
            let bloque = reordenamientos::juntar_casillas_vecinas(
                inicio,
                &casillas,
                &mut visitadas,
                pares,
                k,
            );
            candidatos.push((inicio, bloque, visitadas.len()));
        }
    }
    let contexto = ContextoDeBloques {
        datos,
        pares,
        casillas: &casillas,
    };
    let mut guardados: Vec<BloqueSintenico> = Vec::new();
    if !pares.is_empty() {
        narrar_casillas(&mut narracion, pares, datos, &casillas);
        narracion.agregar(
            base_de(&variables_de_bloques(&contexto, 0, 0)).paso(
                "let mut visitadas: HashSet<Casilla> = HashSet::new();",
                "Ahora se recorren las parejas y se agrupan las casillas vecinas. Para no contar una casilla dos veces se anotan las que ya están en algún bloque; al principio, ninguna.",
            ),
        );
        let primera_casilla = candidatos[0].0;
        let mut con_inicio = variables_de_bloques(&contexto, 0, 0);
        con_inicio.push(("inicio", casilla_escrita(&primera_casilla)));
        narracion.agregar(base_de(&variables_de_bloques(&contexto, 0, 0)).paso(
            "for par in pares {",
            format!(
                "Se toma la primera pareja, la de A {} y B {}.",
                posicion_mostrada(pares[0].posicion_a, datos.inicio_a),
                posicion_mostrada(pares[0].posicion_b, datos.inicio_b)
            ),
        ));
        narracion.agregar(base_de(&con_inicio).paso(
            "let inicio = casilla_de(par, distancia);",
            format!(
                "Su casilla es {}: fila = posición en A ÷ {}, columna = posición en B ÷ {}.",
                casilla_escrita(&primera_casilla),
                miles(distancia),
                miles(distancia)
            ),
        ));
        narracion.agregar(
            base_de(&con_inicio).paso(
                "if visitadas.insert(inicio) {",
                "Esa casilla todavía no está en ningún bloque (insert devuelve verdadero): aquí empieza uno nuevo.",
            ),
        );
        narracion.agregar(
            base_de(&con_inicio).paso(
                "let bloque = juntar_casillas_vecinas(inicio, &casillas, &mut visitadas, pares, k);",
                "Se llama a la función que hace crecer el bloque desde esa casilla, uniendo las casillas vecinas que también tengan parejas.",
            ),
        );
    }
    let mut muestreo = Muestreo::nuevo(candidatos.len(), 0, 40);
    let total_candidatos = candidatos.len();
    for (numero, (inicio, real, visitadas_tras)) in candidatos.iter().enumerate() {
        let (bloque, toque) = if numero == 0 {
            let mirado = narrar_primer_bloque(&mut narracion, &contexto, *inicio);
            (mirado, Toque::Detalle)
        } else {
            (real.clone(), muestreo.toca())
        };
        let queda = bloque.kmeros >= minimo;
        let mut variables = variables_de_bloques(&contexto, *visitadas_tras, guardados.len());
        variables.push(("inicio", casilla_escrita(inicio)));
        variables.push(("bloque", resumen_de_bloque(&bloque)));
        if !matches!(toque, Toque::Nada) {
            narracion.agregar(
                tabla_de_bloques(base_de(&variables), &guardados, false).paso(
                    "if bloque.kmeros >= minimo {",
                    format!(
                        "Bloque candidato {} de {} (desde la casilla {}): {}. {}",
                        miles(numero + 1),
                        miles(total_candidatos),
                        casilla_escrita(inicio),
                        resumen_de_bloque(&bloque),
                        if queda {
                            format!("Tiene al menos {} k-meros: se queda.", miles(minimo))
                        } else {
                            format!(
                                "Tiene menos de {} k-meros: es muy pequeño, probablemente casualidad, y se descarta.",
                                miles(minimo)
                            )
                        }
                    ),
                ),
            );
        }
        if queda {
            guardados.push(bloque.clone());
            if !matches!(toque, Toque::Nada) {
                let mut despues = variables_de_bloques(&contexto, *visitadas_tras, guardados.len());
                despues.push(("inicio", casilla_escrita(inicio)));
                despues.push(("bloque", resumen_de_bloque(&bloque)));
                narracion.agregar(tabla_de_bloques(base_de(&despues), &guardados, true).paso(
                    "bloques.push(bloque);",
                    format!(
                        "Se guarda: ya hay {} en la lista.",
                        cuantos(guardados.len(), "bloque", "bloques")
                    ),
                ));
            }
        }
    }
    guardados.sort_by_key(|bloque| (bloque.inicio_a, bloque.inicio_b));
    let con_orden = variables_de_bloques(&contexto, visitadas.len(), guardados.len());
    narracion.agregar(
        tabla_de_bloques(base_de(&con_orden), &guardados, false).paso(
            "bloques.sort_by_key(|bloque| (bloque.inicio_a, bloque.inicio_b));",
            format!(
                "Se recorrieron {} y se guardaron {}. Ahora se ordenan por su inicio en A (y, si empatan, en B); esa es la tabla que ves en el resultado.",
                cuantos(total_candidatos, "candidato", "candidatos"),
                cuantos(guardados.len(), "bloque", "bloques")
            ),
        ),
    );
    let permutacion = reordenamientos::permutacion_de_bloques(&guardados);
    let nota_final = if guardados.is_empty() {
        format!(
            "No quedó ningún bloque con al menos {} k-meros. Prueba con un k menor, una distancia mayor o un mínimo más bajo.",
            miles(minimo)
        )
    } else {
        format!(
            "Resultado: {} de sintenia. Leídos en el orden de B, forman la permutación {}: eso es lo que luego alimenta los algoritmos de reordenamientos. «Qué significa» explica cuántas inversiones harían falta.",
            cuantos(guardados.len(), "bloque", "bloques"),
            parentesis(resumen(&permutacion, con_signo, " "))
        )
    };
    let ultimo = tabla_de_bloques(
        base_de(&[
            ("k", k.to_string()),
            ("pares", format!("{} parejas", miles(pares.len()))),
            ("bloques", format!("{} bloques", miles(guardados.len()))),
        ]),
        &guardados,
        false,
    )
    .paso("Ok(Sintenia { pares, bloques })", nota_final);
    debug_assert_eq!(guardados, sintenia.bloques);
    Ok((narracion.terminar(ultimo), guardados))
}

fn pasos_bloques_de_sintenia(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (a, b, k) = leer_secuencias_y_k(entrada)?;
    let distancia = entrada.natural("distancia")?;
    let minimo = entrada.natural("minimo")?;
    if distancia == 0 {
        return Err("La distancia para unir k-meros vecinos debe ser al menos 1 pb.".to_string());
    }
    let datos = DatosDeSintenia {
        a: &a,
        b: &b,
        k,
        distancia,
        minimo,
        inicio_a: inicio_en_ecoli(entrada, "a"),
        inicio_b: inicio_en_ecoli(entrada, "b"),
    };
    Ok(narrar_sintenia(&datos)?.0)
}

#[cfg(test)]
mod pruebas {
    use super::*;
    use crate::{codigo, ejemplo, pasos};

    fn sorteo(semilla: &mut u64) -> u64 {
        *semilla = semilla
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *semilla >> 33
    }

    fn permutacion_al_azar(n: usize, semilla: &mut u64) -> Vec<i64> {
        let mut bloques: Vec<i64> = (1..=n as i64).collect();
        for i in (1..n).rev() {
            let j = sorteo(semilla) as usize % (i + 1);
            bloques.swap(i, j);
        }
        for bloque in &mut bloques {
            if sorteo(semilla).is_multiple_of(2) {
                *bloque = -*bloque;
            }
        }
        bloques
    }

    fn genoma_al_azar(n: usize, cromosomas: usize, semilla: &mut u64) -> Vec<Vec<i64>> {
        let bloques = permutacion_al_azar(n, semilla);
        let tamano = n.div_ceil(cromosomas);
        bloques.chunks(tamano).map(|trozo| trozo.to_vec()).collect()
    }

    fn secuencia_al_azar(largo: usize, semilla: &mut u64) -> Vec<u8> {
        (0..largo)
            .map(|_| b"ACGT"[sorteo(semilla) as usize % 4])
            .collect()
    }

    fn entrada_de(id: &str) -> Entrada {
        Entrada::desde_json(ejemplo(id).unwrap()).unwrap()
    }

    #[test]
    fn cada_paso_cita_una_linea_del_codigo_mostrado() {
        for algoritmo in algoritmos() {
            let id = algoritmo.id;
            let lista = pasos(id, ejemplo(id).unwrap()).unwrap();
            let texto: String = codigo(id)
                .unwrap()
                .iter()
                .map(|fragmento| fragmento.fuente.as_str())
                .collect::<Vec<_>>()
                .join("\n");
            assert!(!lista.is_empty() && lista.len() <= MAXIMO_PASOS, "{id}");
            for paso in &lista {
                assert!(
                    texto
                        .lines()
                        .any(|linea| linea.contains(paso.linea.as_str())),
                    "{id}: «{}» no está en el código",
                    paso.linea
                );
                assert!(!paso.nota.is_empty(), "{id}: nota vacía");
            }
        }
    }

    #[test]
    fn el_voraz_narrado_llega_al_mismo_resultado() {
        let permutacion = leer_permutacion(
            &entrada_de("reordenamientos.ordenamiento_voraz"),
            "permutacion",
            MAXIMO_LISTADO,
        )
        .unwrap();
        let (pasos, final_, guardadas) = narrar_voraz(&permutacion).unwrap();
        let esperado = reordenamientos::ordenamiento_voraz(&permutacion).unwrap();
        assert_eq!(guardadas, esperado.len());
        assert_eq!(final_, esperado.last().unwrap().clone());
        assert!(pasos.len() > 10);
        let mut semilla = 7;
        for n in [1, 2, 9, 150, 700] {
            let permutacion = permutacion_al_azar(n, &mut semilla);
            let (pasos, final_, guardadas) = narrar_voraz(&permutacion).unwrap();
            let esperado = reordenamientos::ordenamiento_voraz(&permutacion).unwrap();
            assert_eq!(guardadas, esperado.len());
            assert_eq!(final_, esperado.last().cloned().unwrap_or(permutacion));
            assert!(pasos.len() <= MAXIMO_PASOS);
        }
    }

    #[test]
    fn las_rupturas_narradas_coinciden() {
        let mut semilla = 11;
        for n in [1, 6, 80, 5000] {
            let permutacion = permutacion_al_azar(n, &mut semilla);
            let (pasos, encontradas) = narrar_rupturas(&permutacion).unwrap();
            assert_eq!(
                encontradas,
                reordenamientos::rupturas(&permutacion).unwrap()
            );
            assert!(pasos.len() <= MAXIMO_PASOS);
        }
        let identidad: Vec<i64> = (1..=5).collect();
        assert!(narrar_rupturas(&identidad).unwrap().1.is_empty());
    }

    #[test]
    fn los_cambios_entre_cromosoma_y_ciclo_narrados_coinciden() {
        let mut semilla = 3;
        for n in [1, 4, 90, 3000] {
            let cromosoma = permutacion_al_azar(n, &mut semilla);
            let (pasos, ciclo) = narrar_cromosoma_a_ciclo(&cromosoma).unwrap();
            assert_eq!(
                ciclo,
                reordenamientos::cromosoma_a_ciclo(&cromosoma).unwrap()
            );
            assert!(pasos.len() <= MAXIMO_PASOS);
            let (pasos, vuelta) = narrar_ciclo_a_cromosoma(&ciclo).unwrap();
            assert_eq!(vuelta, cromosoma);
            assert!(pasos.len() <= MAXIMO_PASOS);
        }
    }

    #[test]
    fn las_aristas_y_el_grafo_narrados_coinciden() {
        let mut semilla = 5;
        for (n, cromosomas) in [(1, 1), (6, 2), (60, 3), (900, 4), (400, 400)] {
            let genoma = genoma_al_azar(n, cromosomas, &mut semilla);
            let (pasos, aristas) = narrar_aristas_de_color(&genoma).unwrap();
            assert_eq!(aristas, reordenamientos::aristas_de_color(&genoma).unwrap());
            assert!(pasos.len() <= MAXIMO_PASOS);
            let (pasos, vuelta) = narrar_grafo_a_genoma(&aristas).unwrap();
            assert_eq!(vuelta, reordenamientos::grafo_a_genoma(&aristas).unwrap());
            assert!(pasos.len() <= MAXIMO_PASOS);
        }
    }

    #[test]
    fn la_distancia_narrada_coincide() {
        let mut semilla = 13;
        for (n, cromosomas) in [(1, 1), (8, 2), (70, 3), (800, 5)] {
            let p = genoma_al_azar(n, cromosomas, &mut semilla);
            let mut q = genoma_al_azar(n, cromosomas + 1, &mut semilla);
            if n == 1 {
                q = p.clone();
            }
            let (pasos, distancia, ciclos) = narrar_distancia(&p, &q).unwrap();
            assert_eq!(
                distancia,
                reordenamientos::distancia_dos_roturas(&p, &q).unwrap()
            );
            assert_eq!(
                ciclos,
                reordenamientos::ciclos_del_grafo_de_rupturas(&p, &q).unwrap()
            );
            assert!(pasos.len() <= MAXIMO_PASOS);
        }
        let mismo = vec![vec![1, 2, 3]];
        assert_eq!(narrar_distancia(&mismo, &mismo).unwrap().1, 0);
    }

    #[test]
    fn las_dos_roturas_narradas_coinciden() {
        let aristas = vec![(2, 4), (3, 8), (7, 5), (6, 1)];
        let (pasos, nuevas) = narrar_dos_roturas_en_grafo(&aristas, [1, 6, 3, 8]).unwrap();
        assert_eq!(
            nuevas,
            reordenamientos::dos_roturas_en_grafo(&aristas, 1, 6, 3, 8).unwrap()
        );
        assert!(pasos.len() >= 8);
        let genoma = vec![vec![1, 2, -4, 3, 5, -6]];
        let (pasos, resultado) = narrar_dos_roturas_en_genoma(&genoma, [4, 8, 6, 9]).unwrap();
        assert_eq!(
            resultado,
            reordenamientos::dos_roturas_en_genoma(&genoma, 4, 8, 6, 9).unwrap()
        );
        assert!(pasos.len() >= 8);
        assert!(narrar_dos_roturas_en_grafo(&aristas, [1, 6, 1, 6]).is_err());
    }

    #[test]
    fn los_kmeros_y_la_sintenia_narrados_coinciden() {
        let mut semilla = 17;
        let comun = secuencia_al_azar(300, &mut semilla);
        let mut a = secuencia_al_azar(200, &mut semilla);
        a.extend_from_slice(&comun);
        a.extend(secuencia_al_azar(200, &mut semilla));
        let mut b = secuencia_al_azar(100, &mut semilla);
        b.extend(complemento_inverso(&String::from_utf8_lossy(&comun)).bytes());
        b.extend(secuencia_al_azar(150, &mut semilla));
        b.extend_from_slice(&comun[..120]);
        for k in [1, 5, 11] {
            let datos = DatosDeKmeros {
                a: &a,
                b: &b,
                k,
                inicio_a: None,
                inicio_b: Some(10),
            };
            let (pasos, pares) = narrar_kmeros(&datos).unwrap();
            assert_eq!(
                pares,
                reordenamientos::kmeros_compartidos(&a, &b, k).unwrap()
            );
            assert!(pasos.len() <= MAXIMO_PASOS);
        }
        for (k, distancia, minimo) in [(11, 50, 5), (11, 50, 100_000), (15, 1, 1)] {
            let datos = DatosDeSintenia {
                a: &a,
                b: &b,
                k,
                distancia,
                minimo,
                inicio_a: None,
                inicio_b: None,
            };
            let (pasos, bloques) = narrar_sintenia(&datos).unwrap();
            assert_eq!(
                bloques,
                reordenamientos::bloques_de_sintenia(&a, &b, k, distancia, minimo)
                    .unwrap()
                    .bloques
            );
            assert!(pasos.len() <= MAXIMO_PASOS);
        }
        let ninguno = DatosDeKmeros {
            a: b"AAAAAAAA",
            b: b"CCCCCCCC",
            k: 4,
            inicio_a: None,
            inicio_b: None,
        };
        assert!(narrar_kmeros(&ninguno).unwrap().1.is_empty());
    }
}
