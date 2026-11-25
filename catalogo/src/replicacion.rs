use crate::formato::{cientifico, decimal, miles, porcentaje};
use crate::pizarra::{Marca as Resalte, Narracion, Pizarra};
use crate::{Algoritmo, Arista, Bloque, Campo, Entrada, FilaAlineada, Marca, Paso, Salida, Tipo};
use nucleo::replicacion::{self, MAXIMO_D, MAXIMO_K_ENUMERABLE, MAXIMO_K_NUMERICO};
use nucleo::{ecoli, secuencia};
use serde_json::{Value, json};
use std::collections::HashMap;

const ORIC_INICIO: usize = 3_925_744;
const ORIC_FIN: usize = 3_925_975;
const CAJA_DNAA: &[u8] = b"TTATCCACA";
const MAXIMO_PATRON_EXACTO: usize = 1000;
const MAXIMO_PATRON_APROXIMADO: usize = 100;
const MAXIMO_K_TABLA: usize = 64;
const MAXIMO_FILAS: usize = 2000;
const MAXIMO_TEXTO: usize = 100_000;
const MAXIMO_PASOS: usize = 400;
const MAXIMO_PUNTOS_DEL_PASO: usize = 200;
const MAXIMO_MINIMOS_DEL_PASO: usize = 20;
const MAXIMO_ANIMABLE: usize = 20_000_000;
const PUNTOS_DE_CURVA: usize = 2000;
const MAXIMO_NODOS: usize = 60;
const MEJORES: usize = 20;
const TRAMOS: usize = 20;
const SIGNIFICANCIA: f64 = 0.01;
const GENOMA_MINIMO: usize = 100_000;

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![
        Algoritmo {
            id: "replicacion.contar_patron",
            capitulo: "replicacion",
            nombre: "Contar un patrón",
            resumen: "Cuántas veces aparece un patrón en un texto, incluidas las apariciones que se solapan.",
            campos: &[
                Campo {
                    id: "texto",
                    etiqueta: "Texto de ADN",
                    tipo: Tipo::Secuencia,
                    ejemplo: "@ecoli:3925744-3925975",
                },
                Campo {
                    id: "patron",
                    etiqueta: "Patrón",
                    tipo: Tipo::Secuencia,
                    ejemplo: "GATC",
                },
            ],
            funcion: "replicacion::contar_patron",
            ejecutar: contar_patron,
            pasos: Some(pasos_contar_patron),
        },
        Algoritmo {
            id: "replicacion.palabras_frecuentes",
            capitulo: "replicacion",
            nombre: "Palabras más frecuentes",
            resumen: "Los k-meros que más se repiten en un texto, contados en una sola pasada con una tabla de frecuencias.",
            campos: &[
                Campo {
                    id: "texto",
                    etiqueta: "Texto de ADN",
                    tipo: Tipo::Secuencia,
                    ejemplo: "@ecoli:3925744-3925975",
                },
                Campo {
                    id: "k",
                    etiqueta: "k (largo de cada palabra)",
                    tipo: Tipo::Entero,
                    ejemplo: "4",
                },
            ],
            funcion: "replicacion::palabras_frecuentes",
            ejecutar: palabras_frecuentes,
            pasos: Some(pasos_palabras_frecuentes),
        },
        Algoritmo {
            id: "replicacion.complemento_inverso",
            capitulo: "replicacion",
            nombre: "Complemento inverso",
            resumen: "La otra hebra del ADN leída de 5' a 3': se da vuelta la secuencia y se cambia A por T y C por G.",
            campos: &[Campo {
                id: "adn",
                etiqueta: "Secuencia de ADN",
                tipo: Tipo::Secuencia,
                ejemplo: "TTATCCACA",
            }],
            funcion: "secuencia::complemento_inverso",
            ejecutar: complemento_inverso,
            pasos: Some(pasos_complemento_inverso),
        },
        Algoritmo {
            id: "replicacion.posiciones_del_patron",
            capitulo: "replicacion",
            nombre: "Dónde aparece un patrón",
            resumen: "Todas las posiciones, contando desde 0, donde empieza una copia exacta del patrón.",
            campos: &[
                Campo {
                    id: "texto",
                    etiqueta: "Texto de ADN",
                    tipo: Tipo::Secuencia,
                    ejemplo: "@ecoli",
                },
                Campo {
                    id: "patron",
                    etiqueta: "Patrón",
                    tipo: Tipo::Secuencia,
                    ejemplo: "TTATCCACA",
                },
            ],
            funcion: "replicacion::posiciones_del_patron",
            ejecutar: posiciones_del_patron,
            pasos: Some(pasos_posiciones_del_patron),
        },
        Algoritmo {
            id: "replicacion.grumos",
            capitulo: "replicacion",
            nombre: "Grumos de k-meros",
            resumen: "Los k-meros que aparecen al menos t veces dentro de alguna ventana de L bases (clumps), con una ventana que se desliza.",
            campos: &[
                Campo {
                    id: "genoma",
                    etiqueta: "Genoma",
                    tipo: Tipo::Secuencia,
                    ejemplo: "@ecoli",
                },
                Campo {
                    id: "k",
                    etiqueta: "k (largo de cada palabra)",
                    tipo: Tipo::Entero,
                    ejemplo: "9",
                },
                Campo {
                    id: "ventana",
                    etiqueta: "L (largo de la ventana)",
                    tipo: Tipo::Entero,
                    ejemplo: "500",
                },
                Campo {
                    id: "veces",
                    etiqueta: "t (apariciones mínimas)",
                    tipo: Tipo::Entero,
                    ejemplo: "3",
                },
            ],
            funcion: "replicacion::encontrar_grumos",
            ejecutar: grumos,
            pasos: Some(pasos_grumos),
        },
        Algoritmo {
            id: "replicacion.contenido_gc",
            capitulo: "replicacion",
            nombre: "Contenido GC",
            resumen: "Qué fracción de la secuencia son bases G o C.",
            campos: &[Campo {
                id: "adn",
                etiqueta: "Secuencia de ADN",
                tipo: Tipo::Secuencia,
                ejemplo: "@gen:dnaA",
            }],
            funcion: "secuencia::contenido_gc",
            ejecutar: contenido_gc,
            pasos: Some(pasos_contenido_gc),
        },
        Algoritmo {
            id: "replicacion.sesgo_minimo",
            capitulo: "replicacion",
            nombre: "Sesgo GC y su mínimo",
            resumen: "Recorre el genoma sumando 1 por cada G y restando 1 por cada C; donde el total es más bajo suele estar el origen de replicación.",
            campos: &[Campo {
                id: "genoma",
                etiqueta: "Genoma",
                tipo: Tipo::Secuencia,
                ejemplo: "@ecoli",
            }],
            funcion: "replicacion::posiciones_de_sesgo_minimo",
            ejecutar: sesgo_minimo,
            pasos: Some(pasos_sesgo_minimo),
        },
        Algoritmo {
            id: "replicacion.distancia_hamming",
            capitulo: "replicacion",
            nombre: "Distancia de Hamming",
            resumen: "En cuántas posiciones difieren dos secuencias del mismo largo.",
            campos: &[
                Campo {
                    id: "primera",
                    etiqueta: "Primera secuencia",
                    tipo: Tipo::Secuencia,
                    ejemplo: "TTATCCACA",
                },
                Campo {
                    id: "segunda",
                    etiqueta: "Segunda secuencia",
                    tipo: Tipo::Secuencia,
                    ejemplo: "@ecoli:3925907-3925915",
                },
            ],
            funcion: "replicacion::distancia_hamming",
            ejecutar: distancia_hamming,
            pasos: Some(pasos_distancia_hamming),
        },
        Algoritmo {
            id: "replicacion.posiciones_aproximadas",
            capitulo: "replicacion",
            nombre: "Búsqueda aproximada",
            resumen: "Dónde aparece un patrón con como mucho d letras distintas.",
            campos: &[
                Campo {
                    id: "texto",
                    etiqueta: "Texto de ADN",
                    tipo: Tipo::Secuencia,
                    ejemplo: "@ecoli:3925744-3925975",
                },
                Campo {
                    id: "patron",
                    etiqueta: "Patrón",
                    tipo: Tipo::Secuencia,
                    ejemplo: "TTATCCACA",
                },
                Campo {
                    id: "d",
                    etiqueta: "d (diferencias permitidas)",
                    tipo: Tipo::Entero,
                    ejemplo: "1",
                },
            ],
            funcion: "replicacion::posiciones_aproximadas",
            ejecutar: posiciones_aproximadas,
            pasos: Some(pasos_posiciones_aproximadas),
        },
        Algoritmo {
            id: "replicacion.contar_aproximado",
            capitulo: "replicacion",
            nombre: "Conteo aproximado",
            resumen: "Cuántas ventanas del texto se parecen al patrón con como mucho d diferencias.",
            campos: &[
                Campo {
                    id: "texto",
                    etiqueta: "Texto de ADN",
                    tipo: Tipo::Secuencia,
                    ejemplo: "@ecoli",
                },
                Campo {
                    id: "patron",
                    etiqueta: "Patrón",
                    tipo: Tipo::Secuencia,
                    ejemplo: "TGTGGATAA",
                },
                Campo {
                    id: "d",
                    etiqueta: "d (diferencias permitidas)",
                    tipo: Tipo::Entero,
                    ejemplo: "1",
                },
            ],
            funcion: "replicacion::contar_aproximado",
            ejecutar: contar_aproximado,
            pasos: Some(pasos_contar_aproximado),
        },
        Algoritmo {
            id: "replicacion.frecuentes_con_errores",
            capitulo: "replicacion",
            nombre: "Frecuentes con errores",
            resumen: "Los k-meros que más aparecen si cada aparición puede tener hasta d letras distintas.",
            campos: &[
                Campo {
                    id: "texto",
                    etiqueta: "Texto de ADN",
                    tipo: Tipo::Secuencia,
                    ejemplo: "@ecoli:3925744-3925975",
                },
                Campo {
                    id: "k",
                    etiqueta: "k (largo de cada palabra)",
                    tipo: Tipo::Entero,
                    ejemplo: "9",
                },
                Campo {
                    id: "d",
                    etiqueta: "d (diferencias permitidas)",
                    tipo: Tipo::Entero,
                    ejemplo: "1",
                },
            ],
            funcion: "replicacion::palabras_frecuentes_con_errores",
            ejecutar: frecuentes_con_errores,
            pasos: Some(pasos_frecuentes_con_errores),
        },
        Algoritmo {
            id: "replicacion.frecuentes_con_complementos",
            capitulo: "replicacion",
            nombre: "Frecuentes con errores en las dos hebras",
            resumen: "Como el anterior, pero cada k-mero suma también las apariciones de su complemento inverso.",
            campos: &[
                Campo {
                    id: "texto",
                    etiqueta: "Texto de ADN",
                    tipo: Tipo::Secuencia,
                    ejemplo: "@ecoli:3925744-3925975",
                },
                Campo {
                    id: "k",
                    etiqueta: "k (largo de cada palabra)",
                    tipo: Tipo::Entero,
                    ejemplo: "9",
                },
                Campo {
                    id: "d",
                    etiqueta: "d (diferencias permitidas)",
                    tipo: Tipo::Entero,
                    ejemplo: "1",
                },
            ],
            funcion: "replicacion::palabras_frecuentes_con_errores_y_complementos",
            ejecutar: frecuentes_con_complementos,
            pasos: Some(pasos_frecuentes_con_complementos),
        },
        Algoritmo {
            id: "replicacion.arreglo_de_frecuencias",
            capitulo: "replicacion",
            nombre: "Arreglo de frecuencias",
            resumen: "Una casilla por cada uno de los 4^k k-meros posibles, en orden alfabético, con cuántas veces aparece cada uno.",
            campos: &[
                Campo {
                    id: "texto",
                    etiqueta: "Texto de ADN",
                    tipo: Tipo::Secuencia,
                    ejemplo: "@gen:dnaA",
                },
                Campo {
                    id: "k",
                    etiqueta: "k (largo de cada palabra)",
                    tipo: Tipo::Entero,
                    ejemplo: "2",
                },
            ],
            funcion: "replicacion::arreglo_de_frecuencias",
            ejecutar: arreglo_de_frecuencias,
            pasos: Some(pasos_arreglo_de_frecuencias),
        },
        Algoritmo {
            id: "replicacion.numero_de_patron",
            capitulo: "replicacion",
            nombre: "De patrón a número",
            resumen: "El lugar de un k-mero en el orden alfabético, leyéndolo como un número en base 4 (A=0, C=1, G=2, T=3).",
            campos: &[Campo {
                id: "patron",
                etiqueta: "Patrón (hasta 31 letras)",
                tipo: Tipo::Secuencia,
                ejemplo: "GATTACA",
            }],
            funcion: "replicacion::numero_de_patron",
            ejecutar: numero_de_patron,
            pasos: Some(pasos_numero_de_patron),
        },
        Algoritmo {
            id: "replicacion.patron_de_numero",
            capitulo: "replicacion",
            nombre: "De número a patrón",
            resumen: "El camino inverso: escribe el número en base 4 con k cifras y cambia cada cifra por su letra.",
            campos: &[
                Campo {
                    id: "numero",
                    etiqueta: "Número",
                    tipo: Tipo::Entero,
                    ejemplo: "9000",
                },
                Campo {
                    id: "k",
                    etiqueta: "k (letras del patrón)",
                    tipo: Tipo::Entero,
                    ejemplo: "7",
                },
            ],
            funcion: "replicacion::patron_de_numero",
            ejecutar: patron_de_numero,
            pasos: Some(pasos_patron_de_numero),
        },
        Algoritmo {
            id: "replicacion.vecindario",
            capitulo: "replicacion",
            nombre: "Vecindario de un patrón",
            resumen: "Todos los patrones que difieren del dado en como mucho d letras.",
            campos: &[
                Campo {
                    id: "patron",
                    etiqueta: "Patrón (hasta 31 letras)",
                    tipo: Tipo::Secuencia,
                    ejemplo: "TTATCCACA",
                },
                Campo {
                    id: "d",
                    etiqueta: "d (diferencias permitidas)",
                    tipo: Tipo::Entero,
                    ejemplo: "1",
                },
            ],
            funcion: "replicacion::vecindario",
            ejecutar: vecindario,
            pasos: Some(pasos_vecindario),
        },
        Algoritmo {
            id: "replicacion.buscar_ori",
            capitulo: "replicacion",
            nombre: "Buscar el origen de replicación",
            resumen: "Todo el capítulo en cadena: mínimo del sesgo GC, una ventana de 500 pb y los 9-meros más frecuentes con errores en las dos hebras, candidatos a caja DnaA.",
            campos: &[
                Campo {
                    id: "genoma",
                    etiqueta: "Genoma bacteriano completo",
                    tipo: Tipo::Secuencia,
                    ejemplo: "@ecoli",
                },
                Campo {
                    id: "ancho",
                    etiqueta: "Ancho de la ventana (pb)",
                    tipo: Tipo::Entero,
                    ejemplo: "500",
                },
                Campo {
                    id: "k",
                    etiqueta: "k (largo de la caja)",
                    tipo: Tipo::Entero,
                    ejemplo: "9",
                },
                Campo {
                    id: "d",
                    etiqueta: "d (diferencias permitidas)",
                    tipo: Tipo::Entero,
                    ejemplo: "1",
                },
            ],
            funcion: "replicacion::buscar_ori",
            ejecutar: buscar_ori,
            pasos: Some(pasos_buscar_ori),
        },
    ]
}

fn contar_patron(entrada: &Entrada) -> Result<Salida, String> {
    let texto = entrada.secuencia("texto")?;
    let patron = leer_patron(entrada, "patron", MAXIMO_PATRON_EXACTO)?;
    let veces = replicacion::contar_patron(&texto, &patron);
    let posiciones = replicacion::posiciones_del_patron(&texto, &patron);
    let inverso = secuencia::complemento_inverso(&patron);
    let es_palindromo = inverso == patron;
    let en_la_otra = if es_palindromo {
        0
    } else {
        replicacion::contar_patron(&texto, &inverso)
    };
    let esperadas = esperadas_por_composicion(&texto, &patron);
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Apariciones", miles(veces)));
    salida.agregar(Bloque::valor(
        "Esperadas si las bases estuvieran al azar",
        cantidad_esperada(esperadas),
    ));
    let mut categorias = vec![format!("{} (observadas)", letras(&patron))];
    let mut valores = vec![veces as f64];
    if !es_palindromo {
        salida.agregar(Bloque::valor(
            &format!("Apariciones de {} (la otra hebra)", letras(&inverso)),
            miles(en_la_otra),
        ));
        categorias.push(format!("{} (otra hebra)", letras(&inverso)));
        valores.push(en_la_otra as f64);
    }
    categorias.push("Esperadas por azar".to_string());
    valores.push(esperadas);
    salida.agregar(Bloque::barras(
        "Observado frente a esperado",
        "",
        "Apariciones",
        categorias,
        valores,
    ));
    comparar_con_el_azar(&mut salida, veces, esperadas, MODELO_COMPOSICION);
    if let Some((primera, segunda)) = primer_solapamiento(&posiciones, patron.len()) {
        salida.explicar(format!(
            "Algunas apariciones se solapan (por ejemplo, las que empiezan en {} y {}, contando desde 0): se cuentan todas, como pide el libro.",
            miles(primera),
            miles(segunda)
        ));
    }
    if es_palindromo && patron.len() >= 4 {
        explicar_palindromo(&mut salida, &patron);
    } else if en_la_otra > 0 {
        salida.explicar(format!(
            "En la otra hebra, leída como {}, aparece {} veces más: una proteína que reconoce {} lo encuentra en cualquiera de las dos hebras.",
            letras(&inverso),
            miles(en_la_otra),
            letras(&patron)
        ));
    }
    explicar_senales_conocidas(&mut salida, &patron);
    Ok(salida)
}

fn palabras_frecuentes(entrada: &Entrada) -> Result<Salida, String> {
    let texto = entrada.adn("texto")?;
    let k = leer_k(entrada, "k", texto.len(), MAXIMO_K_TABLA)?;
    let ResumenDeConteo {
        frecuentes,
        maximo,
        distintos,
        mejores,
        histograma: (categorias, valores),
    } = resumir_conteo(&texto, k)?;
    let ventanas = texto.len() - k + 1;
    let proporcion = proporciones(&texto);
    let mut salida = Salida::default();
    salida.agregar(Bloque::lista(
        "Más frecuentes",
        frecuentes.iter().map(|kmero| letras(kmero)),
    ));
    salida.agregar(Bloque::valor("Veces que aparece cada uno", miles(maximo)));
    salida.agregar(Bloque::valor(
        "k-meros distintos",
        format!("{} de {} ventanas", miles(distintos), miles(ventanas)),
    ));
    let filas = mejores
        .iter()
        .map(|(kmero, veces)| {
            vec![
                letras(kmero),
                miles(*veces),
                cantidad_esperada(ventanas as f64 * probabilidad_del_patron(kmero, &proporcion)),
            ]
        })
        .collect();
    salida.agregar(Bloque::tabla(
        "Los más repetidos",
        &["k-mero", "Veces", "Esperadas por azar"],
        filas,
    ));
    salida.agregar(Bloque::barras(
        "Cuántos k-meros distintos aparecen n veces",
        "Veces",
        "k-meros distintos",
        categorias,
        valores,
    ));
    let lista = unir(frecuentes.iter().take(5).map(|kmero| letras(kmero)));
    if maximo <= 1 {
        salida.explicar(format!(
            "Ningún {k}-mero se repite: cada uno de los {} aparece una sola vez. Con k = {k} en un texto de {} letras eso es lo normal; prueba con un k más pequeño o busca permitiendo errores.",
            miles(distintos),
            miles(texto.len())
        ));
    } else {
        let posibles = 4f64.powi(k as i32);
        let por_azar = posibles * poisson_al_menos(ventanas as f64 / posibles, maximo);
        if por_azar < SIGNIFICANCIA {
            salida.explicar(format!(
                "En un texto al azar de este largo (con las cuatro bases igual de probables) esperarías {} k-meros distintos que aparezcan {} veces o más. Que {lista} lo haga es una señal, no casualidad.",
                cantidad_esperada(por_azar),
                miles(maximo)
            ));
        } else if por_azar >= 1.0 {
            salida.explicar(format!(
                "Por azar ya esperarías unos {} k-meros que aparezcan {} veces o más: con k = {k} y {} letras las repeticiones son normales. Prueba con un k más grande.",
                cantidad_esperada(por_azar),
                miles(maximo),
                miles(texto.len())
            ));
        } else {
            salida.explicar(format!(
                "Por azar esperarías {} k-meros que aparezcan {} veces o más: es algo raro, pero no lo bastante como para asegurar que {lista} sea una señal.",
                cantidad_esperada(por_azar),
                miles(maximo)
            ));
        }
    }
    if maximo > 1 && frecuentes.len() >= 5 {
        salida.explicar(format!(
            "Hay {} k-meros empatados en el primer lugar: ninguno destaca sobre los demás.",
            miles(frecuentes.len())
        ));
    }
    if maximo > 1
        && let Some(simple) = frecuentes.iter().find(|kmero| es_repeticion_simple(kmero))
    {
        salida.explicar(format!(
            "{} es una repetición simple (una o dos letras que se repiten): se acumula porque se solapa consigo mismo, no necesariamente porque cumpla una función.",
            letras(simple)
        ));
    }
    if maximo > 1
        && let Some(primero) = frecuentes.first()
    {
        explicar_senales_conocidas(&mut salida, primero);
    }
    Ok(salida)
}

struct ResumenDeConteo {
    frecuentes: Vec<Vec<u8>>,
    maximo: usize,
    distintos: usize,
    mejores: Vec<(Vec<u8>, usize)>,
    histograma: (Vec<String>, Vec<f64>),
}

fn resumir_conteo(texto: &[u8], k: usize) -> Result<ResumenDeConteo, String> {
    if replicacion::conviene_un_arreglo(texto.len(), k) {
        let arreglo = replicacion::arreglo_de_frecuencias(texto, k)?;
        let mut mejores = Vec::new();
        for (numero, veces) in mejores_del_arreglo(&arreglo, MEJORES) {
            mejores.push((
                replicacion::patron_de_numero(numero as u64, k)?,
                veces as usize,
            ));
        }
        let presentes = arreglo.iter().filter(|&&veces| veces > 0);
        return Ok(ResumenDeConteo {
            frecuentes: replicacion::los_mas_frecuentes(&arreglo, k),
            maximo: arreglo.iter().copied().max().unwrap_or(0) as usize,
            distintos: presentes.clone().count(),
            mejores,
            histograma: histograma_de_veces(presentes.map(|&veces| veces as usize)),
        });
    }
    let tabla = replicacion::tabla_de_frecuencias(texto, k);
    let maximo = tabla.values().copied().max().unwrap_or(0);
    let mut frecuentes: Vec<Vec<u8>> = tabla
        .iter()
        .filter(|(_, veces)| **veces == maximo)
        .map(|(kmero, _)| kmero.to_vec())
        .collect();
    frecuentes.sort();
    Ok(ResumenDeConteo {
        frecuentes,
        maximo,
        distintos: tabla.len(),
        mejores: mejores_de_la_tabla(&tabla, MEJORES)
            .into_iter()
            .map(|(kmero, veces)| (kmero.to_vec(), veces))
            .collect(),
        histograma: histograma_de_veces(tabla.values().copied()),
    })
}

fn complemento_inverso(entrada: &Entrada) -> Result<Salida, String> {
    let adn = entrada.secuencia("adn")?;
    let inverso = secuencia::complemento_inverso(&adn);
    let mut salida = Salida::default();
    salida.agregar(Bloque::texto(
        "Complemento inverso (de 5' a 3')",
        recortar(&inverso),
    ));
    salida.agregar(Bloque::valor(
        "Longitud",
        format!("{} pb", miles(adn.len())),
    ));
    if adn.len() <= 200 {
        let complemento: Vec<u8> = adn
            .iter()
            .map(|&base| secuencia::complemento(base))
            .collect();
        salida.agregar(Bloque::alineamiento(
            "Las dos hebras, base con base",
            vec![
                FilaAlineada {
                    nombre: "Tu hebra, 5'→3'".to_string(),
                    secuencia: letras(&adn),
                },
                FilaAlineada {
                    nombre: "Hebra de enfrente, 3'←5'".to_string(),
                    secuencia: letras(&complemento),
                },
            ],
        ));
    }
    if inverso == adn {
        explicar_palindromo(&mut salida, &adn);
    } else if adn.len() <= 50 {
        salida.explicar(format!(
            "Las dos hebras corren en sentidos opuestos: lo que en una hebra se lee {}, en la otra se lee {} (las dos de 5' a 3'). Por eso, al buscar una señal en el ADN, también hay que buscar su complemento inverso.",
            letras(&adn),
            letras(&inverso)
        ));
    } else {
        salida.explicar(
            "En una secuencia larga, el complemento inverso es la hebra de enfrente leída de 5' a 3': así se leen, por ejemplo, los genes que están en la hebra −.",
        );
    }
    let raras = adn
        .iter()
        .filter(|base| !matches!(base, b'A' | b'C' | b'G' | b'T'))
        .count();
    if raras > 0 {
        salida.explicar(format!(
            "{} no {} A, C, G ni T (como N): se dejan como están, porque no tienen una pareja definida.",
            plural(raras, "letra", "letras"),
            if raras == 1 { "es" } else { "son" }
        ));
    }
    if adn == CAJA_DNAA || inverso == CAJA_DNAA {
        salida.explicar(
            "TTATCCACA es la caja DnaA de E. coli y TGTGGATAA es la misma caja vista desde la otra hebra; en oriC hay cajas en las dos orientaciones.",
        );
    }
    if adn.len() > MAXIMO_TEXTO {
        salida.explicar(format!(
            "Se muestran las primeras {} letras de {}.",
            miles(MAXIMO_TEXTO),
            miles(adn.len())
        ));
    }
    Ok(salida)
}

fn posiciones_del_patron(entrada: &Entrada) -> Result<Salida, String> {
    let texto = entrada.secuencia("texto")?;
    let patron = leer_patron(entrada, "patron", MAXIMO_PATRON_EXACTO)?;
    let posiciones = replicacion::posiciones_del_patron(&texto, &patron);
    let inverso = secuencia::complemento_inverso(&patron);
    let en_la_otra = if inverso == patron {
        Vec::new()
    } else {
        replicacion::posiciones_del_patron(&texto, &inverso)
    };
    let inicio = inicio_en_ecoli(entrada, "texto");
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Apariciones", miles(posiciones.len())));
    salida.agregar(Bloque::lista(
        "Posiciones (contando desde 0)",
        posiciones.iter().map(|&posicion| miles(posicion)),
    ));
    if let Some(inicio) = inicio {
        let filas = posiciones
            .iter()
            .take(MAXIMO_FILAS)
            .map(|&posicion| vec![miles(posicion), miles(inicio + posicion)])
            .collect();
        salida.agregar(tabla_con_total(
            "En el genoma de E. coli",
            &["Posición (desde 0)", "Coordenada en E. coli (desde 1)"],
            filas,
            posiciones.len(),
        ));
    }
    if !posiciones.is_empty() {
        let (categorias, valores) = contar_por_tramos(posiciones.iter().copied(), texto.len());
        salida.agregar(Bloque::barras(
            "Apariciones a lo largo del texto",
            "Tramo",
            "Apariciones",
            categorias,
            valores,
        ));
    }
    let esperadas = esperadas_por_composicion(&texto, &patron);
    comparar_con_el_azar(&mut salida, posiciones.len(), esperadas, MODELO_COMPOSICION);
    let (juntas, desde) = mas_juntas(&posiciones, 500usize.saturating_sub(patron.len() - 1));
    if juntas >= 3 {
        salida.explicar(format!(
            "En una ventana de 500 pb que empieza en la posición {} llegan a juntarse {} apariciones: eso es un grumo, lo que busca el algoritmo «Grumos de k-meros».",
            miles(desde),
            juntas
        ));
    }
    if let Some(inicio) = inicio
        && contiene_oric(inicio, texto.len())
    {
        let coordenadas: Vec<usize> = posiciones
            .iter()
            .map(|&posicion| inicio + posicion)
            .collect();
        let dentro = coordenadas
            .iter()
            .filter(|&&coordenada| (ORIC_INICIO..=ORIC_FIN).contains(&coordenada))
            .count();
        if dentro > 0 {
            salida.explicar(format!(
                "{} de ellas {} dentro de oriC (coordenadas {}–{}), el origen de replicación anotado de E. coli.",
                dentro,
                if dentro == 1 { "cae" } else { "caen" },
                miles(ORIC_INICIO),
                miles(ORIC_FIN)
            ));
        } else if let Some(cercana) = coordenadas
            .iter()
            .min_by_key(|&&coordenada| distancia_a_oric(coordenada))
        {
            salida.explicar(format!(
                "Ninguna cae dentro de oriC ({}–{}); la más cercana empieza en la coordenada {}, a {} pb. Las cajas reales del origen no son copias exactas: por eso el libro pasa a buscar con errores.",
                miles(ORIC_INICIO),
                miles(ORIC_FIN),
                miles(*cercana),
                miles(distancia_a_oric(*cercana))
            ));
        }
    }
    if !en_la_otra.is_empty() {
        salida.explicar(format!(
            "Su complemento inverso, {}, aparece {} veces: son copias del mismo patrón en la otra hebra.",
            letras(&inverso),
            miles(en_la_otra.len())
        ));
    }
    if let Some(inicio) = inicio {
        salida.explicar(format!(
            "Las posiciones cuentan desde 0, como en el libro; la anotación de E. coli cuenta desde 1, por eso en la tabla cada coordenada es la posición + {}.",
            miles(inicio)
        ));
    }
    Ok(salida)
}

fn grumos(entrada: &Entrada) -> Result<Salida, String> {
    let genoma = entrada.secuencia("genoma")?;
    let k = leer_k(entrada, "k", genoma.len(), MAXIMO_K_ENUMERABLE)?;
    let ventana = entrada.natural("ventana")?;
    if ventana < k || ventana > genoma.len() {
        return Err(format!(
            "La ventana L debe medir entre k ({k}) y el largo del genoma ({} pb).",
            miles(genoma.len())
        ));
    }
    let caben = ventana - k + 1;
    let veces = entrada.natural("veces")?;
    if veces == 0 || veces > caben {
        return Err(format!(
            "t debe estar entre 1 y {}: en una ventana de {} pb caben {} k-meros de {k} letras.",
            miles(caben),
            miles(ventana),
            miles(caben)
        ));
    }
    let grumos = replicacion::encontrar_grumos(&genoma, k, ventana, veces)?;
    let inicio = inicio_en_ecoli(entrada, "genoma");
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        "k-meros que forman grumos",
        miles(grumos.len()),
    ));
    let mut columnas = vec!["k-mero", "Primera ventana (desde 0)"];
    if inicio.is_some() {
        columnas.push("Coordenada en E. coli (desde 1)");
    }
    let filas = grumos
        .iter()
        .take(MAXIMO_FILAS)
        .map(|grumo| {
            let mut fila = vec![letras(&grumo.patron), miles(grumo.ventana)];
            if let Some(inicio) = inicio {
                fila.push(miles(inicio + grumo.ventana));
            }
            fila
        })
        .collect();
    salida.agregar(tabla_con_total("Grumos", &columnas, filas, grumos.len()));
    let (categorias, valores) =
        contar_por_tramos(grumos.iter().map(|grumo| grumo.ventana), genoma.len());
    let tramos_con_grumos = valores.iter().filter(|&&valor| valor > 0.0).count();
    if !grumos.is_empty() {
        salida.agregar(Bloque::barras(
            "Dónde empieza la primera ventana de cada grumo",
            "Tramo",
            "Grumos",
            categorias.clone(),
            valores,
        ));
    }
    if grumos.is_empty() {
        salida.explicar(format!(
            "Ningún {k}-mero aparece {} dentro de una misma ventana de {} pb. Prueba con un t menor o con una ventana más grande.",
            plural(veces, "vez", "veces"),
            miles(ventana)
        ));
    }
    let simples: Vec<&replicacion::Grumo> = grumos
        .iter()
        .filter(|grumo| es_repeticion_simple(&grumo.patron))
        .collect();
    if let Some(ejemplo) = simples.first() {
        salida.explicar(format!(
            "{} de los {} son repeticiones simples (como {}): se juntan porque se solapan consigo mismos, no necesariamente porque sean una señal.",
            miles(simples.len()),
            miles(grumos.len()),
            letras(&ejemplo.patron)
        ));
    }
    if genoma.len() >= GENOMA_MINIMO && tramos_con_grumos * 4 >= categorias.len() * 3 {
        salida.explicar(format!(
            "Hay grumos repartidos por casi todo el genoma ({} de {} tramos tienen alguno): un grumo solo no alcanza para ubicar el origen de replicación. Por eso el libro pasa a mirar el sesgo GC.",
            tramos_con_grumos,
            categorias.len()
        ));
    } else if !grumos.is_empty() && genoma.len() >= 10_000 && tramos_con_grumos <= 3 {
        salida.explicar(format!(
            "Los grumos se concentran en {} de {} tramos: mira esas zonas primero.",
            tramos_con_grumos,
            categorias.len()
        ));
    }
    if inicio == Some(1) && genoma.len() == ecoli::genoma().len() {
        let desde = ORIC_INICIO.saturating_sub(ventana);
        let hasta = (ORIC_FIN - 1 + ventana).min(genoma.len());
        let alrededor = replicacion::encontrar_grumos(&genoma[desde..hasta], k, ventana, veces)?;
        if alrededor.is_empty() {
            salida.explicar(format!(
                "Ninguna ventana de {} pb que toque oriC tiene un {k}-mero repetido {}: las cajas DnaA del origen no son copias exactas, y por eso hacen falta las búsquedas con errores.",
                miles(ventana),
                plural(veces, "vez", "veces")
            ));
        } else {
            salida.explicar(format!(
                "Alrededor de oriC hay {} k-meros que forman grumos, entre ellos {}.",
                miles(alrededor.len()),
                unir(alrededor.iter().take(3).map(|grumo| letras(&grumo.patron)))
            ));
        }
    }
    Ok(salida)
}

fn contenido_gc(entrada: &Entrada) -> Result<Salida, String> {
    let adn = entrada.secuencia("adn")?;
    let gc = secuencia::contenido_gc(&adn);
    let gc_ecoli = secuencia::contenido_gc(ecoli::genoma());
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Contenido GC", porcentaje(gc, 1)));
    salida.agregar(Bloque::valor(
        "Longitud",
        format!("{} pb", miles(adn.len())),
    ));
    if adn.len() < 100 {
        salida.explicar(
            "La secuencia es muy corta (menos de 100 pb): unas pocas bases cambian mucho el porcentaje.",
        );
    }
    let diferencia = gc - gc_ecoli;
    if diferencia.abs() <= 0.05 {
        salida.explicar(format!(
            "Está cerca del promedio del genoma de E. coli ({}). Nada fuera de lo común.",
            porcentaje(gc_ecoli, 1)
        ));
    } else {
        let sentido = if diferencia > 0.0 {
            "más rica"
        } else {
            "más pobre"
        };
        salida.explicar(format!(
            "Es {sentido} en GC que el promedio de E. coli ({}). Una región con GC muy distinto al resto del genoma puede venir de otro organismo (transferencia horizontal).",
            porcentaje(gc_ecoli, 1)
        ));
    }
    if gc > 0.6 {
        salida.explicar(
            "Más de 60 % de GC: el ADN es más estable (G y C se unen con tres puentes de hidrógeno) y puede costar amplificarlo por PCR.",
        );
    }
    if gc < 0.4 {
        salida.explicar(
            "Menos de 40 % de GC: es típico de regiones ricas en A y T, como promotores u orígenes de replicación, que se abren con facilidad.",
        );
    }
    Ok(salida)
}

fn sesgo_minimo(entrada: &Entrada) -> Result<Salida, String> {
    let genoma = entrada.secuencia("genoma")?;
    let minimos = replicacion::posiciones_de_sesgo_minimo(&genoma);
    let recorrido = recorrer_sesgo(&genoma);
    let inicio = inicio_en_ecoli(entrada, "genoma");
    let primero = minimos[0];
    let ultimo = minimos[minimos.len() - 1];
    let mut marcas = vec![Marca {
        x: primero as f64,
        texto: "mínimo".to_string(),
    }];
    marcas.push(Marca {
        x: recorrido.posicion_del_maximo as f64,
        texto: "máximo".to_string(),
    });
    if let Some(inicio) = inicio
        && contiene_oric(inicio, genoma.len())
    {
        marcas.push(Marca {
            x: (ORIC_INICIO - inicio) as f64,
            texto: "oriC anotado".to_string(),
        });
    }
    let mut salida = Salida::default();
    salida.agregar(Bloque::linea(
        "Sesgo G − C acumulado",
        "Posición (pb, desde 0)",
        "G − C",
        &recorrido.puntos,
        marcas,
    ));
    salida.agregar(Bloque::valor("Sesgo mínimo", con_signo(recorrido.minimo)));
    salida.agregar(Bloque::lista(
        "Posiciones del mínimo (contando desde 0)",
        minimos.iter().map(|&posicion| miles(posicion)),
    ));
    salida.agregar(Bloque::valor(
        "Sesgo máximo",
        format!(
            "{} en la posición {}",
            con_signo(recorrido.maximo),
            miles(recorrido.posicion_del_maximo)
        ),
    ));
    if let Some(inicio) = inicio {
        salida.agregar(Bloque::valor(
            "El mínimo en coordenadas de E. coli (desde 1)",
            format!("justo después de la base {}", miles(inicio - 1 + primero)),
        ));
    }
    let longitud = genoma.len();
    let gc = secuencia::contenido_gc(&genoma);
    let escala_del_azar = (longitud as f64 * gc).sqrt().max(1.0);
    let amplitud = (recorrido.maximo - recorrido.minimo) as f64;
    let veces_el_azar = amplitud / escala_del_azar;
    if longitud < 10_000 {
        salida.explicar(format!(
            "La secuencia es corta ({} pb): el sesgo GC solo dice algo sobre el origen cuando recorres un genoma bacteriano entero, de millones de bases; aquí manda el azar.",
            miles(longitud)
        ));
    }
    if veces_el_azar >= 5.0 && longitud >= GENOMA_MINIMO {
        salida.explicar(format!(
            "La curva sube y baja {} unidades, unas {} veces lo que daría el azar (alrededor de {}): hay un sesgo real entre G y C.",
            miles(amplitud as usize),
            decimal(veces_el_azar, 0),
            miles(escala_del_azar as usize)
        ));
        salida.explicar(
            "Del origen al término, la hebra que lees sirve de molde a la hebra rezagada y pasa más tiempo sin pareja; ahí las C se convierten en T con más facilidad, así que quedan más G que C y el sesgo sube. Del término al origen pasa lo contrario: por eso el mínimo señala el origen (ori) y el máximo, el término (ter).",
        );
        let separacion = primero.abs_diff(recorrido.posicion_del_maximo);
        let vuelta = separacion.min(longitud - separacion) as f64 / longitud as f64;
        if (0.4..=0.5).contains(&vuelta) {
            salida.explicar(format!(
                "El máximo queda casi en el lado opuesto del cromosoma circular ({} de vuelta desde el mínimo): ahí suele estar el término, donde se encuentran las dos horquillas de replicación.",
                porcentaje(vuelta, 0)
            ));
        }
    } else if veces_el_azar < 3.0 {
        salida.explicar(format!(
            "La curva sube y baja {} unidades, y por puro azar esperarías unas {}: el mínimo no es una señal fiable.",
            miles(amplitud as usize),
            miles(escala_del_azar as usize)
        ));
    } else {
        salida.explicar(format!(
            "La curva se aleja del azar ({} veces la variación típica), pero menos que en un genoma bacteriano completo: tómalo como una pista débil.",
            decimal(veces_el_azar, 1)
        ));
    }
    if let Some(inicio) = inicio
        && contiene_oric(inicio, longitud)
    {
        let coordenada = inicio - 1 + primero;
        let distancia = distancia_a_oric(coordenada);
        if distancia == 0 {
            salida.explicar(format!(
                "El mínimo cae dentro de oriC, el origen anotado de E. coli ({}–{}).",
                miles(ORIC_INICIO),
                miles(ORIC_FIN)
            ));
        } else {
            salida.explicar(format!(
                "El mínimo queda justo después de la base {} (contando desde 1), a {} pb de oriC, el origen anotado ({}–{}). El sesgo acierta la zona, no la base exacta: el siguiente paso es buscar cajas DnaA cerca.",
                miles(coordenada),
                miles(distancia),
                miles(ORIC_INICIO),
                miles(ORIC_FIN)
            ));
        }
    }
    if minimos.len() > 1 {
        if minimos.len() == ultimo - primero + 1 {
            salida.explicar(format!(
                "El mínimo se repite en {} posiciones seguidas (de {} a {}): entre ellas no hay ni G ni C, que son las únicas letras que mueven el sesgo.",
                miles(minimos.len()),
                miles(primero),
                miles(ultimo)
            ));
        } else if ultimo - primero > longitud / 100 {
            salida.explicar(format!(
                "El mínimo se repite en {} posiciones separadas hasta por {} pb: el origen podría estar en cualquiera de esos tramos.",
                miles(minimos.len()),
                miles(ultimo - primero)
            ));
        } else {
            salida.explicar(format!(
                "El mínimo se repite en {} posiciones muy juntas (de {} a {}): el sesgo baja a ese valor, sube un poco y vuelve a tocarlo, así que el origen está en ese tramo.",
                miles(minimos.len()),
                miles(primero),
                miles(ultimo)
            ));
        }
    }
    Ok(salida)
}

struct Momento {
    desde: usize,
    hasta: usize,
    sesgo_antes: i64,
    sesgo: i64,
    minimo: i64,
    nuevo_minimo: bool,
    empate: bool,
    primera_del_minimo: usize,
}

fn pasos_sesgo_minimo(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let genoma = entrada.secuencia("genoma")?;
    let longitud = genoma.len();
    if longitud > MAXIMO_ANIMABLE {
        return Err(format!(
            "El paso a paso admite hasta {} bases y esta secuencia tiene {}.",
            miles(MAXIMO_ANIMABLE),
            miles(longitud)
        ));
    }
    let recorrido = recorrer_sesgo(&genoma);
    let cantidad = (longitud + 1).min(MAXIMO_PASOS);
    let por_base = cantidad == longitud + 1;
    let secuencia_completa = if por_base {
        Value::String(letras(&genoma))
    } else {
        Value::Null
    };
    let mut pasos = Vec::with_capacity(cantidad);
    let mut curva: Vec<[i64; 2]> = Vec::with_capacity(cantidad);
    let mut sesgo_actual: i64 = 0;
    let mut minimo: i64 = 0;
    let mut primeras: Vec<usize> = vec![0];
    let mut total_minimos: usize = 1;
    let mut leidas = 0;
    for indice in 0..cantidad {
        let marca = if cantidad == 1 {
            0
        } else {
            indice * longitud / (cantidad - 1)
        };
        let desde = leidas;
        let sesgo_antes = sesgo_actual;
        let minimo_antes = minimo;
        let total_antes = total_minimos;
        while leidas < marca {
            sesgo_actual += replicacion::cambio_de_sesgo(genoma[leidas]);
            leidas += 1;
            if sesgo_actual < minimo {
                minimo = sesgo_actual;
                primeras.clear();
                primeras.push(leidas);
                total_minimos = 1;
            } else if sesgo_actual == minimo {
                if primeras.len() < MAXIMO_MINIMOS_DEL_PASO {
                    primeras.push(leidas);
                }
                total_minimos += 1;
            }
        }
        curva.push([marca as i64, sesgo_actual]);
        let nuevo_minimo = minimo < minimo_antes;
        let momento = Momento {
            desde,
            hasta: marca,
            sesgo_antes,
            sesgo: sesgo_actual,
            minimo,
            nuevo_minimo,
            empate: !nuevo_minimo && total_minimos > total_antes,
            primera_del_minimo: primeras[0],
        };
        let (linea, mut nota) = if indice == 0 {
            (
                "let mut posiciones = vec![0];",
                "Antes de leer ninguna base el sesgo vale 0, y la posición 0 es el primer mínimo."
                    .to_string(),
            )
        } else if por_base {
            describir_base(genoma[marca - 1], &momento)
        } else {
            describir_tramo(&momento)
        };
        if indice + 1 == cantidad {
            let donde = if total_minimos == 1 {
                format!("en la posición {}", miles(primeras[0]))
            } else {
                format!("en {} posiciones", miles(total_minimos))
            };
            nota.push_str(&format!(
                " Fin del recorrido: el mínimo es {} y está {donde}.",
                con_signo(minimo)
            ));
        }
        let base = if por_base && marca > 0 {
            Value::String((genoma[marca - 1] as char).to_string())
        } else {
            Value::Null
        };
        let estado = json!({
            "longitud": longitud,
            "posicion": marca,
            "bases_leidas_en_el_paso": marca - desde,
            "base": base,
            "sesgo": sesgo_actual,
            "minimo": minimo,
            "posiciones_minimo": primeras,
            "total_minimos": total_minimos,
            "eje_y": [recorrido.minimo, recorrido.maximo],
            "curva": muestra_pareja(&curva, MAXIMO_PUNTOS_DEL_PASO),
            "secuencia": secuencia_completa,
        });
        pasos.push(Paso::nuevo(linea, nota, estado));
    }
    Ok(pasos)
}

fn describir_base(base: u8, momento: &Momento) -> (&'static str, String) {
    let letra = base as char;
    let posicion = miles(momento.hasta);
    let sesgo = con_signo(momento.sesgo);
    if momento.nuevo_minimo {
        return (
            "minimo = sesgo_actual;",
            format!(
                "Lees una {letra}: el sesgo baja a {sesgo}, el valor más bajo hasta ahora (posición {posicion})."
            ),
        );
    }
    if momento.empate {
        return (
            "} else if sesgo_actual == minimo {",
            format!(
                "Lees una {letra}: el sesgo queda en {sesgo} e iguala el mínimo; la posición {posicion} se suma a la lista."
            ),
        );
    }
    match replicacion::cambio_de_sesgo(base) {
        1 => (
            "b'G' | b'g' => 1,",
            format!("Lees una G: el sesgo sube a {sesgo}."),
        ),
        -1 => (
            "b'C' | b'c' => -1,",
            format!(
                "Lees una C: el sesgo baja a {sesgo}, pero no pasa del mínimo ({}).",
                con_signo(momento.minimo)
            ),
        ),
        _ => (
            "_ => 0,",
            format!("Lees una {letra}: no es G ni C, así que el sesgo sigue en {sesgo}."),
        ),
    }
}

fn describir_tramo(momento: &Momento) -> (&'static str, String) {
    let desde = miles(momento.desde);
    let hasta = miles(momento.hasta);
    if momento.nuevo_minimo {
        return (
            "minimo = sesgo_actual;",
            format!(
                "Entre las posiciones {desde} y {hasta} el sesgo bajó a un mínimo nuevo: {} (por primera vez en la posición {}).",
                con_signo(momento.minimo),
                miles(momento.primera_del_minimo)
            ),
        );
    }
    if momento.empate {
        return (
            "} else if sesgo_actual == minimo {",
            format!(
                "Entre las posiciones {desde} y {hasta} el sesgo volvió a tocar el mínimo ({}); esas posiciones se suman a la lista.",
                con_signo(momento.minimo)
            ),
        );
    }
    (
        "sesgo_actual += cambio_de_sesgo(base);",
        format!(
            "Entre las posiciones {desde} y {hasta} el sesgo pasó de {} a {}; el mínimo sigue en {}.",
            con_signo(momento.sesgo_antes),
            con_signo(momento.sesgo),
            con_signo(momento.minimo)
        ),
    )
}

fn distancia_hamming(entrada: &Entrada) -> Result<Salida, String> {
    let primera = entrada.secuencia("primera")?;
    let segunda = entrada.secuencia("segunda")?;
    let distancia = replicacion::distancia_hamming(&primera, &segunda)?;
    let longitud = primera.len();
    let mut transiciones = 0;
    let mut transversiones = 0;
    let mut filas = Vec::new();
    for (posicion, (&a, &b)) in primera.iter().zip(&segunda).enumerate() {
        if a == b {
            continue;
        }
        let tipo = tipo_de_cambio(a, b);
        match tipo {
            "transición" => transiciones += 1,
            "transversión" => transversiones += 1,
            _ => {}
        }
        if filas.len() < MAXIMO_FILAS {
            filas.push(vec![
                miles(posicion),
                (a as char).to_string(),
                (b as char).to_string(),
                tipo.to_string(),
            ]);
        }
    }
    let proporcion = distancia as f64 / longitud as f64;
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Distancia de Hamming", miles(distancia)));
    salida.agregar(Bloque::valor(
        "Posiciones distintas",
        format!(
            "{} de {} ({})",
            miles(distancia),
            miles(longitud),
            porcentaje(proporcion, 1)
        ),
    ));
    if longitud <= MAXIMO_FILAS {
        salida.agregar(Bloque::alineamiento(
            "Letra con letra",
            vec![
                FilaAlineada {
                    nombre: "Primera".to_string(),
                    secuencia: letras(&primera),
                },
                FilaAlineada {
                    nombre: "Segunda".to_string(),
                    secuencia: letras(&segunda),
                },
            ],
        ));
    }
    if distancia > 0 {
        salida.agregar(tabla_con_total(
            "Dónde difieren",
            &["Posición (desde 0)", "Primera", "Segunda", "Tipo de cambio"],
            filas,
            distancia,
        ));
    }
    if distancia == 0 {
        salida.explicar("Son idénticas: no hay ninguna posición distinta.");
    } else if longitud < 20 {
        salida.explicar(format!(
            "Con secuencias tan cortas cada letra pesa mucho: {} en {} letras es el {}.",
            plural(distancia, "diferencia", "diferencias"),
            miles(longitud),
            porcentaje(proporcion, 0)
        ));
    } else if proporcion >= 0.6 {
        salida.explicar(format!(
            "Difieren en el {}: casi como dos secuencias sin relación (con cuatro letras al azar se espera que difieran en el 75 %). Si una tiene una letra de más o de menos, la de Hamming exagera: conviene alinearlas (capítulo 5).",
            porcentaje(proporcion, 0)
        ));
    } else if proporcion <= 0.1 {
        salida.explicar(format!(
            "Difieren en solo el {}: es muy probable que compartan origen o función.",
            porcentaje(proporcion, 1)
        ));
    } else {
        salida.explicar(format!(
            "Difieren en el {}: más parecidas que dos secuencias al azar (75 %), pero con bastantes cambios.",
            porcentaje(proporcion, 0)
        ));
    }
    if transiciones + transversiones >= 10 {
        let razon = transiciones as f64 / transversiones.max(1) as f64;
        if razon > 1.0 {
            salida.explicar(format!(
                "Hay más transiciones ({transiciones}, A↔G o C↔T) que transversiones ({transversiones}): es lo típico de mutaciones reales, porque cambiar una purina por otra purina es más fácil. Al azar habría el doble de transversiones."
            ));
        } else if razon <= 0.6 {
            salida.explicar(format!(
                "Transiciones ({transiciones}) y transversiones ({transversiones}) están en la proporción que daría el azar, 1 a 2: no se nota la huella de la mutación natural."
            ));
        }
    }
    let inversa = secuencia::complemento_inverso(CAJA_DNAA);
    if longitud == CAJA_DNAA.len()
        && [&primera, &segunda]
            .iter()
            .any(|cual| cual.as_slice() == CAJA_DNAA || **cual == inversa)
    {
        salida.explicar(
            "TTATCCACA es el consenso de la caja DnaA de E. coli. Las cajas reales del origen suelen diferir en una o dos letras, y DnaA igual se une a ellas: por eso el libro busca patrones permitiendo errores.",
        );
    }
    Ok(salida)
}

fn posiciones_aproximadas(entrada: &Entrada) -> Result<Salida, String> {
    let texto = entrada.secuencia("texto")?;
    let patron = leer_patron(entrada, "patron", MAXIMO_PATRON_APROXIMADO)?;
    let d = leer_d(entrada, "d", patron.len(), patron.len())?;
    revisar_largo(&patron, &texto)?;
    let posiciones = replicacion::posiciones_aproximadas(&texto, &patron, d);
    let inverso = secuencia::complemento_inverso(&patron);
    let es_palindromo = inverso == patron;
    let en_la_otra = if es_palindromo {
        Vec::new()
    } else {
        replicacion::posiciones_aproximadas(&texto, &inverso, d)
    };
    let inicio = inicio_en_ecoli(entrada, "texto");
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        &format!(
            "Apariciones con hasta {}",
            plural(d, "diferencia", "diferencias")
        ),
        miles(posiciones.len()),
    ));
    let mut columnas = vec!["Posición (desde 0)"];
    if inicio.is_some() {
        columnas.push("Coordenada en E. coli (desde 1)");
    }
    columnas.extend(["Ventana", "Diferencias"]);
    let mut por_diferencias = vec![0.0; d + 1];
    let mut filas = Vec::new();
    for &posicion in &posiciones {
        let ventana = &texto[posicion..posicion + patron.len()];
        let diferencias = diferencias_entre(ventana, &patron);
        por_diferencias[diferencias] += 1.0;
        if filas.len() < MAXIMO_FILAS {
            let mut fila = vec![miles(posicion)];
            if let Some(inicio) = inicio {
                fila.push(miles(inicio + posicion));
            }
            fila.push(letras(ventana));
            fila.push(diferencias.to_string());
            filas.push(fila);
        }
    }
    salida.agregar(tabla_con_total("Dónde", &columnas, filas, posiciones.len()));
    salida.agregar(Bloque::barras(
        "Cuántas tienen 0, 1, … diferencias",
        "Diferencias",
        "Apariciones",
        (0..=d).map(|cantidad| cantidad.to_string()).collect(),
        por_diferencias,
    ));
    if !es_palindromo {
        salida.agregar(Bloque::valor(
            &format!("En la otra hebra (como {})", letras(&inverso)),
            miles(en_la_otra.len()),
        ));
    }
    let ventanas = texto.len() - patron.len() + 1;
    let probabilidad = probabilidad_con_errores(patron.len(), d);
    let hebras = if es_palindromo { 1.0 } else { 2.0 };
    comparar_con_el_azar(
        &mut salida,
        posiciones.len() + en_la_otra.len(),
        ventanas as f64 * probabilidad * hebras,
        if es_palindromo {
            MODELO_UNIFORME
        } else {
            MODELO_UNIFORME_DOS_HEBRAS
        },
    );
    if !es_palindromo && !en_la_otra.is_empty() {
        salida.explicar(format!(
            "Buscando su complemento inverso, {}, aparecen {} más: son parecidos al patrón en la otra hebra.",
            letras(&inverso),
            miles(en_la_otra.len())
        ));
    }
    if d * 2 >= patron.len() {
        salida.explicar(format!(
            "Con d = {d} en un patrón de {} letras se puede cambiar la mitad o más: casi cualquier ventana se parece y el resultado pierde sentido.",
            patron.len()
        ));
    }
    explicar_senales_conocidas(&mut salida, &patron);
    Ok(salida)
}

fn contar_aproximado(entrada: &Entrada) -> Result<Salida, String> {
    let texto = entrada.secuencia("texto")?;
    let patron = leer_patron(entrada, "patron", MAXIMO_PATRON_APROXIMADO)?;
    let d = leer_d(entrada, "d", patron.len(), patron.len())?;
    revisar_largo(&patron, &texto)?;
    let veces = replicacion::contar_aproximado(&texto, &patron, d);
    let exactas = replicacion::contar_patron(&texto, &patron);
    let ventanas = texto.len() - patron.len() + 1;
    let esperadas = ventanas as f64 * probabilidad_con_errores(patron.len(), d);
    let mut por_diferencias = vec![0.0; d + 1];
    for posicion in replicacion::posiciones_aproximadas(&texto, &patron, d) {
        let ventana = &texto[posicion..posicion + patron.len()];
        por_diferencias[diferencias_entre(ventana, &patron)] += 1.0;
    }
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        &format!(
            "Ventanas con hasta {}",
            plural(d, "diferencia", "diferencias")
        ),
        miles(veces),
    ));
    salida.agregar(Bloque::valor("De ellas, exactas", miles(exactas)));
    salida.agregar(Bloque::valor(
        "Esperadas si las bases estuvieran al azar",
        cantidad_esperada(esperadas),
    ));
    salida.agregar(Bloque::barras(
        "Cuántas tienen 0, 1, … diferencias",
        "Diferencias",
        "Ventanas",
        (0..=d).map(|cantidad| cantidad.to_string()).collect(),
        por_diferencias,
    ));
    comparar_con_el_azar(&mut salida, veces, esperadas, MODELO_UNIFORME);
    if veces > exactas && d > 0 {
        salida.explicar(format!(
            "Solo {} de las {} son copias exactas: buscando sin errores verías apenas una parte de las señales.",
            miles(exactas),
            miles(veces)
        ));
    }
    if d * 2 >= patron.len() {
        salida.explicar(format!(
            "Con d = {d} en un patrón de {} letras se puede cambiar la mitad o más: casi cualquier ventana se parece y el número pierde sentido.",
            patron.len()
        ));
    }
    explicar_senales_conocidas(&mut salida, &patron);
    Ok(salida)
}

fn frecuentes_con_errores(entrada: &Entrada) -> Result<Salida, String> {
    let texto = entrada.secuencia("texto")?;
    let k = leer_k(entrada, "k", texto.len(), MAXIMO_K_ENUMERABLE)?;
    let d = leer_d(entrada, "d", k, MAXIMO_D)?;
    let cuentas = replicacion::cuentas_con_errores(&texto, k, d)?;
    let frecuentes = replicacion::los_mas_frecuentes(&cuentas, k);
    let exactas = replicacion::arreglo_de_frecuencias(&texto, k)?;
    let mut salida = describir_frecuentes(&texto, k, d, &cuentas, &frecuentes, &exactas, false)?;
    let sin_pareja = frecuentes.iter().all(|patron| {
        frecuentes
            .binary_search(&secuencia::complemento_inverso(patron))
            .is_err()
    });
    if !frecuentes.is_empty() && sin_pareja {
        salida.explicar(
            "Ninguno de los ganadores tiene a su complemento inverso entre los ganadores. Como DnaA se une a sus cajas en las dos hebras, prueba «Frecuentes con errores en las dos hebras».",
        );
    }
    Ok(salida)
}

fn frecuentes_con_complementos(entrada: &Entrada) -> Result<Salida, String> {
    let texto = entrada.secuencia("texto")?;
    let k = leer_k(entrada, "k", texto.len(), MAXIMO_K_ENUMERABLE)?;
    let d = leer_d(entrada, "d", k, MAXIMO_D)?;
    let cuentas = replicacion::cuentas_en_ambas_hebras(&texto, k, d)?;
    let frecuentes = replicacion::los_mas_frecuentes(&cuentas, k);
    let mut exactas = replicacion::arreglo_de_frecuencias(&texto, k)?;
    let otra_hebra =
        replicacion::arreglo_de_frecuencias(&secuencia::complemento_inverso(&texto), k)?;
    for (cuenta, otra) in exactas.iter_mut().zip(otra_hebra) {
        *cuenta += otra;
    }
    let mut salida = describir_frecuentes(&texto, k, d, &cuentas, &frecuentes, &exactas, true)?;
    if frecuentes.len() > 1 {
        salida.explicar(
            "Los ganadores vienen en parejas: un k-mero y su complemento inverso suman lo mismo, porque son la misma señal vista desde cada hebra.",
        );
    }
    Ok(salida)
}

fn describir_frecuentes(
    texto: &[u8],
    k: usize,
    d: usize,
    cuentas: &[u32],
    frecuentes: &[Vec<u8>],
    exactas: &[u32],
    dos_hebras: bool,
) -> Result<Salida, String> {
    let maximo = cuentas.iter().copied().max().unwrap_or(0);
    let mejores = mejores_del_arreglo(cuentas, MEJORES);
    let mut salida = Salida::default();
    salida.agregar(Bloque::lista(
        &format!(
            "Más frecuentes con hasta {}",
            plural(d, "diferencia", "diferencias")
        ),
        frecuentes.iter().map(|patron| letras(patron)),
    ));
    salida.agregar(Bloque::valor(
        "Apariciones de cada uno",
        miles(maximo as usize),
    ));
    let mut columnas = vec!["k-mero"];
    if dos_hebras {
        columnas.push("Complemento inverso");
    }
    columnas.extend([
        "Con hasta d diferencias",
        "Exactas",
        "Diferencias con la caja DnaA",
    ]);
    let mut filas = Vec::new();
    let mut categorias = Vec::new();
    let mut valores = Vec::new();
    for &(numero, veces) in &mejores {
        let patron = replicacion::patron_de_numero(numero as u64, k)?;
        let mut fila = vec![letras(&patron)];
        if dos_hebras {
            fila.push(letras(&secuencia::complemento_inverso(&patron)));
        }
        fila.push(miles(veces as usize));
        fila.push(miles(exactas[numero] as usize));
        fila.push(
            diferencias_con_dnaa(&patron)
                .map(|diferencias| diferencias.to_string())
                .unwrap_or_else(|| "—".to_string()),
        );
        filas.push(fila);
        categorias.push(letras(&patron));
        valores.push(veces as f64);
    }
    salida.agregar(Bloque::tabla(
        "Los 20 con más apariciones",
        &columnas,
        filas,
    ));
    if !categorias.is_empty() {
        salida.agregar(Bloque::barras(
            "Los 20 con más apariciones",
            "k-mero",
            "Apariciones",
            categorias,
            valores,
        ));
    }
    if maximo == 0 {
        salida.explicar("No hay ninguna ventana con solo A, C, G y T: no se pudo contar nada.");
        return Ok(salida);
    }
    let ventanas = texto.len() - k + 1;
    let hebras = if dos_hebras { 2.0 } else { 1.0 };
    let media = ventanas as f64 * probabilidad_con_errores(k, d) * hebras;
    let posibles = 4f64.powi(k as i32);
    let por_azar = posibles * poisson_al_menos(media, maximo as usize);
    let lista = unir(frecuentes.iter().take(4).map(|patron| letras(patron)));
    salida.explicar(format!(
        "Con hasta {}, cada {k}-mero aparecería por azar unas {} veces en un texto de este largo; {lista} {} a {}.",
        plural(d, "diferencia", "diferencias"),
        cantidad_esperada(media),
        if frecuentes.len() == 1 { "llega" } else { "llegan" },
        miles(maximo as usize)
    ));
    if por_azar < SIGNIFICANCIA {
        salida.explicar(format!(
            "Entre los {} k-meros posibles, por azar esperarías {} que lleguen a {}: el ganador es un candidato serio a señal.",
            miles(posibles as usize),
            cantidad_esperada(por_azar),
            miles(maximo as usize)
        ));
    } else if por_azar >= 1.0 {
        salida.explicar(format!(
            "Entre los {} k-meros posibles, por azar ya esperarías unos {} que lleguen a {}: el ganador no destaca de verdad.",
            miles(posibles as usize),
            cantidad_esperada(por_azar),
            miles(maximo as usize)
        ));
    }
    if frecuentes.len() >= 6 {
        salida.explicar(format!(
            "Hay {} k-meros empatados en el primer lugar: el conteo solo no elige uno, y ahí ayuda saber cómo es la señal que buscas.",
            miles(frecuentes.len())
        ));
    }
    let parecidos: Vec<String> = frecuentes
        .iter()
        .filter(|patron| diferencias_con_dnaa(patron).is_some_and(|diferencias| diferencias <= 1))
        .map(|patron| letras(patron))
        .collect();
    if !parecidos.is_empty() {
        salida.explicar(format!(
            "{} {} a 1 letra o menos de la caja DnaA de E. coli (TTATCCACA o TGTGGATAA): buenos candidatos a cajas DnaA.",
            unir(parecidos.iter().take(6).cloned()),
            if parecidos.len() == 1 { "está" } else { "están" }
        ));
    }
    Ok(salida)
}

fn arreglo_de_frecuencias(entrada: &Entrada) -> Result<Salida, String> {
    let texto = entrada.secuencia("texto")?;
    let k = leer_k(entrada, "k", texto.len(), MAXIMO_K_ENUMERABLE)?;
    let arreglo = replicacion::arreglo_de_frecuencias(&texto, k)?;
    let posibles = arreglo.len();
    let presentes = arreglo.iter().filter(|&&veces| veces > 0).count();
    let contadas: usize = arreglo.iter().map(|&veces| veces as usize).sum();
    let ventanas = texto.len() - k + 1;
    let mut salida = Salida::default();
    if posibles <= 4096 {
        let numeros: Vec<String> = arreglo.iter().map(|veces| veces.to_string()).collect();
        salida.agregar(Bloque::texto(
            "Arreglo de frecuencias (de la casilla 0 a la última)",
            numeros.join(" "),
        ));
    }
    if k == 2 {
        let nombres: Vec<String> = ["A", "C", "G", "T"].iter().map(|b| b.to_string()).collect();
        let valores = arreglo
            .chunks(4)
            .map(|fila| fila.iter().map(|&veces| veces as f64).collect())
            .collect();
        salida.agregar(Bloque::calor(
            "Dinucleótidos: primera letra (fila) y segunda (columna)",
            nombres.clone(),
            nombres,
            valores,
        ));
    }
    if k <= 3 {
        let mut categorias = Vec::new();
        for numero in 0..posibles {
            categorias.push(letras(&replicacion::patron_de_numero(numero as u64, k)?));
        }
        salida.agregar(Bloque::barras(
            "Cuántas veces aparece cada k-mero",
            "k-mero",
            "Veces",
            categorias,
            arreglo.iter().map(|&veces| veces as f64).collect(),
        ));
    } else {
        let mut filas = Vec::new();
        for (numero, &veces) in arreglo.iter().enumerate() {
            if veces > 0 && filas.len() < MAXIMO_FILAS {
                filas.push(vec![
                    miles(numero),
                    letras(&replicacion::patron_de_numero(numero as u64, k)?),
                    miles(veces as usize),
                ]);
            }
        }
        salida.agregar(tabla_con_total(
            "Casillas que no están en 0",
            &["Casilla (número)", "k-mero", "Veces"],
            filas,
            presentes,
        ));
    }
    salida.agregar(Bloque::valor("Casillas (4^k)", miles(posibles)));
    salida.agregar(Bloque::valor(
        "Casillas que no están en 0",
        miles(presentes),
    ));
    let primero = letras(&replicacion::patron_de_numero(0, k)?);
    let ultimo = letras(&replicacion::patron_de_numero(posibles as u64 - 1, k)?);
    salida.explicar(format!(
        "La casilla i guarda cuántas veces aparece el k-mero número i del orden alfabético (A=0, C=1, G=2, T=3): la primera es {primero} y la última {ultimo}."
    ));
    if ventanas < posibles {
        salida.explicar(format!(
            "Hay {} k-meros posibles y solo {} ventanas en el texto: la mayoría de las casillas tiene que quedar en 0. Con k grande, una tabla con solo los que aparecen gasta mucha menos memoria.",
            miles(posibles),
            miles(ventanas)
        ));
    } else if presentes == posibles {
        salida.explicar(format!(
            "Aparecen los {} k-meros posibles.",
            miles(posibles)
        ));
    } else {
        let faltan = posibles - presentes;
        let faltarian = posibles as f64 * (-(contadas as f64) / posibles as f64).exp();
        if faltarian < 1.0 && faltan >= 5 {
            salida.explicar(format!(
                "Faltan {} k-meros, y si las bases estuvieran al azar faltarían unos {}: esos k-meros parecen evitados por el genoma.",
                miles(faltan),
                cantidad_esperada(faltarian)
            ));
        } else {
            salida.explicar(format!(
                "Faltan {} de los {} k-meros posibles, algo esperable con {} ventanas.",
                miles(faltan),
                miles(posibles),
                miles(contadas)
            ));
        }
    }
    if k == 2 && texto.len() >= 200 {
        let proporcion = proporciones(&texto);
        let esperado_cg = proporcion[1] * proporcion[2] * contadas as f64;
        if esperado_cg > 0.0 {
            let razon = arreglo[6] as f64 / esperado_cg;
            if razon < 0.6 {
                salida.explicar(format!(
                    "CG aparece {} veces lo esperado según la cantidad de C y de G: está empobrecido, como en los genomas de vertebrados, donde la C de los CG suele estar metilada y muta a T.",
                    decimal(razon, 2)
                ));
            } else {
                salida.explicar(format!(
                    "CG aparece {} veces lo esperado según la cantidad de C y de G: no está empobrecido, como es normal en bacterias (en vertebrados suele quedar por debajo de 0.6).",
                    decimal(razon, 2)
                ));
            }
        }
    }
    if k == 1 && contadas > 0 {
        salida.explicar(format!(
            "Con k = 1 el arreglo es la composición de bases: el contenido GC es {}.",
            porcentaje((arreglo[1] + arreglo[2]) as f64 / contadas as f64, 1)
        ));
    }
    if contadas < ventanas {
        salida.explicar(format!(
            "{} con letras distintas de A, C, G y T no se cuentan.",
            plural(ventanas - contadas, "ventana", "ventanas")
        ));
    }
    Ok(salida)
}

fn numero_de_patron(entrada: &Entrada) -> Result<Salida, String> {
    let patron = entrada.adn("patron")?;
    let numero = replicacion::numero_de_patron(&patron)?;
    let k = patron.len();
    let posibles: u64 = 1 << (2 * k);
    let mut acumulado: u64 = 0;
    let mut filas = Vec::new();
    for &base in &patron {
        let valor = replicacion::valor_de_base(base).unwrap_or(0);
        let siguiente = acumulado * 4 + valor;
        filas.push(vec![
            (base as char).to_string(),
            valor.to_string(),
            format!(
                "4 × {} + {} = {}",
                miles(acumulado as usize),
                valor,
                miles(siguiente as usize)
            ),
        ]);
        acumulado = siguiente;
    }
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Número", miles(numero as usize)));
    salida.agregar(Bloque::tabla(
        "Cómo se calcula, letra por letra",
        &["Letra", "Valor", "Número hasta aquí"],
        filas,
    ));
    salida.agregar(Bloque::valor(
        &format!("Patrones posibles con {k} letras"),
        miles(posibles as usize),
    ));
    salida.explicar(format!(
        "Leer el patrón como un número en base 4 (A=0, C=1, G=2, T=3) da su lugar en la lista alfabética de los {} patrones de {k} letras, empezando en 0.",
        miles(posibles as usize)
    ));
    let cuarto = ["primer", "segundo", "tercer", "último"][(numero / (posibles / 4)) as usize];
    salida.explicar(format!(
        "Como empieza con {}, cae en el {cuarto} cuarto de la lista: los primeros {} números son los que empiezan con A, los siguientes con C, y así.",
        patron[0] as char,
        miles((posibles / 4) as usize)
    ));
    if numero == 0 {
        salida.explicar("Es el primero de la lista: todas sus letras son A.");
    } else if numero == posibles - 1 {
        salida.explicar("Es el último de la lista: todas sus letras son T.");
    }
    if k <= MAXIMO_K_ENUMERABLE {
        salida.explicar(format!(
            "Este número es la casilla del arreglo de frecuencias donde se cuenta {}.",
            letras(&patron)
        ));
    }
    Ok(salida)
}

fn patron_de_numero(entrada: &Entrada) -> Result<Salida, String> {
    let numero = entrada.natural("numero")? as u64;
    let k = entrada.natural("k")?;
    if k == 0 {
        return Err("k debe ser al menos 1.".to_string());
    }
    if k > MAXIMO_K_NUMERICO {
        return Err(format!(
            "k puede valer como mucho {MAXIMO_K_NUMERICO}: un número de 64 bits no guarda patrones más largos."
        ));
    }
    let patron = replicacion::patron_de_numero(numero, k)?;
    let mut filas = Vec::new();
    let mut resto = numero;
    for posicion in (0..k).rev() {
        let cifra = resto % 4;
        filas.push(vec![
            miles(resto as usize),
            miles((resto / 4) as usize),
            cifra.to_string(),
            (patron[posicion] as char).to_string(),
        ]);
        resto /= 4;
    }
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Patrón", letras(&patron)));
    salida.agregar(Bloque::tabla(
        "Divisiones entre 4 (de la última letra a la primera)",
        &["Número", "Cociente", "Resto", "Letra"],
        filas,
    ));
    salida.explicar(format!(
        "Cada resto de dividir entre 4 es una letra (0=A, 1=C, 2=G, 3=T), empezando por la última: así se escribe {} en base 4 con {k} cifras.",
        miles(numero as usize)
    ));
    let cuarta_parte: u64 = 1 << (2 * (k - 1));
    if numero < cuarta_parte {
        salida.explicar(format!(
            "Empieza con A porque {} es menor que {}: los números pequeños se rellenan con A a la izquierda, como los ceros en 0042.",
            miles(numero as usize),
            miles(cuarta_parte as usize)
        ));
    }
    if k < MAXIMO_K_NUMERICO {
        salida.explicar(format!(
            "El mismo número con k = {} da A{}: el k es parte del dato, porque dice cuántas cifras tiene el patrón.",
            k + 1,
            letras(&patron)
        ));
    }
    Ok(salida)
}

fn vecindario(entrada: &Entrada) -> Result<Salida, String> {
    let patron = entrada.adn("patron")?;
    if patron.len() > MAXIMO_K_NUMERICO {
        return Err(format!(
            "El patrón puede tener hasta {MAXIMO_K_NUMERICO} letras, y tiene {}.",
            miles(patron.len())
        ));
    }
    let d = entrada.natural("d")?;
    let vecinos = replicacion::vecindario(&patron, d)?;
    let k = patron.len();
    let mut por_diferencias = vec![0.0; d + 1];
    for vecino in &vecinos {
        por_diferencias[diferencias_entre(vecino, &patron)] += 1.0;
    }
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        "Patrones en el vecindario",
        miles(vecinos.len()),
    ));
    salida.agregar(Bloque::barras(
        "Cuántos difieren en 0, 1, … letras",
        "Diferencias",
        "Patrones",
        (0..=d).map(|cantidad| cantidad.to_string()).collect(),
        por_diferencias,
    ));
    if vecinos.len() <= MAXIMO_NODOS {
        let mut aristas = Vec::new();
        for (uno, primero) in vecinos.iter().enumerate() {
            for (otro, segundo) in vecinos.iter().enumerate().skip(uno + 1) {
                if let Some(posicion) = unica_diferencia(primero, segundo) {
                    aristas.push(Arista {
                        desde: uno,
                        hasta: otro,
                        etiqueta: posicion.to_string(),
                    });
                }
            }
        }
        salida.agregar(Bloque::grafo(
            "Vecinos unidos si difieren en una letra",
            vecinos.iter().map(|vecino| letras(vecino)).collect(),
            aristas,
        ));
    }
    salida.agregar(Bloque::lista(
        "Vecindario (en orden alfabético)",
        vecinos.iter().map(|vecino| letras(vecino)),
    ));
    salida.explicar(format!(
        "Para cada i de 0 a {d} hay que elegir qué i posiciones cambiar y, en cada una, una de las otras 3 letras: sumando todo da {} patrones.",
        miles(vecinos.len())
    ));
    if d < k {
        salida.explicar(format!(
            "Con d = {} serían {}: el vecindario crece muy rápido con d, por eso aquí se limita a un millón de patrones.",
            d + 1,
            miles(replicacion::tamano_del_vecindario(k, d + 1) as usize)
        ));
    }
    if vecinos.len() <= MAXIMO_NODOS && d > 0 {
        salida.explicar(
            "En el grafo, dos patrones están unidos si difieren en una sola letra; el número de cada línea es la posición que cambia, contando desde 0.",
        );
    }
    if patron == CAJA_DNAA && d > 0 {
        salida.explicar(format!(
            "Este es el conjunto de todas las cajas DnaA posibles con hasta {}: buscar «con errores» es buscar todos estos patrones a la vez.",
            plural(d, "error", "errores")
        ));
    }
    Ok(salida)
}

fn buscar_ori(entrada: &Entrada) -> Result<Salida, String> {
    let genoma = entrada.secuencia("genoma")?;
    let ancho = entrada.natural("ancho")?;
    if ancho == 0 || ancho > genoma.len() {
        return Err(format!(
            "El ancho de la ventana debe estar entre 1 y el largo del genoma ({} pb).",
            miles(genoma.len())
        ));
    }
    let k = leer_k(entrada, "k", ancho, MAXIMO_K_ENUMERABLE)?;
    let d = leer_d(entrada, "d", k, MAXIMO_D)?;
    let ori = replicacion::buscar_ori(&genoma, ancho, k, d)?;
    let inicio = inicio_en_ecoli(entrada, "genoma");
    let recorrido = recorrer_sesgo(&genoma);
    let mut marcas = vec![
        Marca {
            x: ori.minimos[0] as f64,
            texto: "mínimo del sesgo".to_string(),
        },
        Marca {
            x: ori.inicio as f64,
            texto: "ventana".to_string(),
        },
    ];
    if let Some(inicio) = inicio
        && contiene_oric(inicio, genoma.len())
    {
        marcas.push(Marca {
            x: (ORIC_INICIO - inicio) as f64,
            texto: "oriC anotado".to_string(),
        });
    }
    let mut salida = Salida::default();
    salida.agregar(Bloque::linea(
        "Paso 1: sesgo G − C",
        "Posición (pb, desde 0)",
        "G − C",
        &recorrido.puntos,
        marcas,
    ));
    salida.agregar(Bloque::valor(
        "Paso 1: mínimo del sesgo (desde 0)",
        miles(ori.minimos[0]),
    ));
    let ventana_texto = match inicio {
        Some(inicio) => format!(
            "posiciones {}–{} (desde 0), coordenadas {}–{} de E. coli",
            miles(ori.inicio),
            miles(ori.fin - 1),
            miles(inicio + ori.inicio),
            miles(inicio + ori.fin - 1)
        ),
        None => format!(
            "posiciones {}–{} (desde 0)",
            miles(ori.inicio),
            miles(ori.fin - 1)
        ),
    };
    salida.agregar(Bloque::valor("Paso 2: ventana analizada", ventana_texto));
    let mut columnas = vec![
        "k-mero",
        "Complemento inverso",
        "Apariciones (dos hebras, con errores)",
    ];
    if k == CAJA_DNAA.len() {
        columnas.push("Diferencias con la caja DnaA");
    }
    let mut filas = Vec::new();
    for candidato in &ori.candidatos {
        let inverso = secuencia::complemento_inverso(candidato);
        if inverso < *candidato && ori.candidatos.binary_search(&inverso).is_ok() {
            continue;
        }
        let mut fila = vec![
            letras(candidato),
            letras(&inverso),
            miles(ori.veces as usize),
        ];
        if let Some(diferencias) = diferencias_con_dnaa(candidato) {
            fila.push(diferencias.to_string());
        }
        filas.push(fila);
    }
    salida.agregar(Bloque::tabla(
        "Paso 3: k-meros más frecuentes en la ventana",
        &columnas,
        filas,
    ));
    let parecidos: Vec<&Vec<u8>> = ori
        .candidatos
        .iter()
        .filter(|candidato| {
            diferencias_con_dnaa(candidato).is_some_and(|diferencias| diferencias <= 1)
        })
        .collect();
    let elegidos: Vec<&Vec<u8>> = if parecidos.is_empty() {
        ori.candidatos.iter().take(1).collect()
    } else {
        parecidos.clone()
    };
    let ventana = &genoma[ori.inicio..ori.fin];
    let mut quien: Vec<u32> = vec![u32::MAX; ventana.len()];
    for (indice, candidato) in elegidos.iter().enumerate() {
        for posicion in replicacion::posiciones_aproximadas(ventana, candidato, d) {
            if quien[posicion] == u32::MAX {
                quien[posicion] = indice as u32;
            }
        }
    }
    let mut filas = Vec::new();
    let mut total_filas = 0;
    for (posicion, dueno) in quien.iter().enumerate() {
        if *dueno == u32::MAX {
            continue;
        }
        total_filas += 1;
        if filas.len() >= MAXIMO_FILAS {
            continue;
        }
        let candidato = elegidos[*dueno as usize];
        let absoluta = ori.inicio + posicion;
        let donde = match inicio {
            Some(inicio) => miles(inicio + absoluta),
            None => miles(absoluta),
        };
        filas.push(vec![
            donde,
            letras(&ventana[posicion..posicion + k]),
            letras(candidato),
            diferencias_entre(&ventana[posicion..posicion + k], candidato).to_string(),
        ]);
    }
    let columna_donde = if inicio.is_some() {
        "Coordenada en E. coli (desde 1)"
    } else {
        "Posición (desde 0)"
    };
    salida.agregar(tabla_con_total(
        "Paso 4: dónde están las cajas candidatas en la ventana",
        &[columna_donde, "Secuencia", "Candidato", "Diferencias"],
        filas,
        total_filas,
    ));
    salida.explicar(format!(
        "La ventana mide {} pb y empieza en el mínimo del sesgo: el mínimo solo marca la zona del origen, y unos cientos de pares de bases alcanzan para un origen bacteriano típico (el de E. coli mide 232 pb).",
        miles(ancho)
    ));
    if ancho > 5000 {
        salida.explicar(format!(
            "Una ventana de {} pb es mucho más grande que un origen típico: además del origen entran miles de bases que no lo son, y el conteo pierde precisión.",
            miles(ancho)
        ));
    }
    if let Some(inicio) = inicio
        && contiene_oric(inicio, genoma.len())
    {
        let coordenada = inicio - 1 + ori.minimos[0];
        let primera = inicio + ori.inicio;
        let ultima = inicio + ori.fin - 1;
        let cubre = primera <= ORIC_INICIO && ultima >= ORIC_FIN;
        salida.explicar(format!(
            "El mínimo del sesgo queda a {} pb de oriC, el origen anotado de E. coli ({}–{}), y la ventana {} oriC entero.",
            miles(distancia_a_oric(coordenada)),
            miles(ORIC_INICIO),
            miles(ORIC_FIN),
            if cubre { "contiene" } else { "no contiene" }
        ));
    }
    let media_por_azar = (ancho - k + 1) as f64 * 2.0 * probabilidad_con_errores(k, d);
    let esperados_por_azar =
        4f64.powi(k as i32) * poisson_al_menos(media_por_azar, ori.veces as usize);
    if esperados_por_azar >= 1.0 {
        salida.explicar(format!(
            "Por azar ya esperarías unos {} {k}-meros que lleguen a {} en una ventana de este largo: el conteo solo es una señal débil. Agranda la ventana o revisa la zona a mano.",
            cantidad_esperada(esperados_por_azar),
            plural(ori.veces as usize, "aparición", "apariciones")
        ));
    }
    if ori.candidatos.len() >= 10 {
        salida.explicar(format!(
            "Hay {} k-meros empatados con {} apariciones: el conteo solo no elige uno, y ahí ayuda saber cómo es la caja DnaA.",
            miles(ori.candidatos.len()),
            ori.veces
        ));
    }
    if !parecidos.is_empty() {
        salida.explicar(format!(
            "{} {} a 1 letra o menos de la caja DnaA de E. coli (TTATCCACA): son los mejores candidatos a cajas DnaA, donde la proteína DnaA se une para abrir el ADN y empezar a copiarlo.",
            unir(parecidos.iter().map(|candidato| letras(candidato))),
            if parecidos.len() == 1 { "está" } else { "están" }
        ));
    } else if k == CAJA_DNAA.len() {
        salida.explicar(
            "Ningún candidato se parece a la caja DnaA de E. coli: en otras bacterias la caja puede ser distinta, o el origen puede estar fuera de la ventana.",
        );
    }
    if genoma.len() < GENOMA_MINIMO {
        salida.explicar(format!(
            "Con menos de {} pb el sesgo no es fiable para ubicar el origen: esta búsqueda tiene sentido con un genoma bacteriano completo.",
            miles(GENOMA_MINIMO)
        ));
    }
    Ok(salida)
}

const MODELO_COMPOSICION: &str =
    "Si las bases estuvieran mezcladas al azar, con la misma proporción de A, C, G y T,";
const MODELO_UNIFORME: &str = "Si las cuatro bases fueran igual de probables,";
const MODELO_UNIFORME_DOS_HEBRAS: &str =
    "Sumando las dos hebras, y si las cuatro bases fueran igual de probables,";

fn comparar_con_el_azar(salida: &mut Salida, observadas: usize, esperadas: f64, modelo: &str) {
    let cantidad = observadas as f64;
    if cantidad > esperadas {
        let probabilidad = poisson_al_menos(esperadas, observadas);
        if probabilidad < SIGNIFICANCIA {
            salida.explicar(format!(
                "{modelo} esperarías unas {} apariciones y hay {}. La probabilidad de ver tantas por casualidad es {}: es muy difícil que sea azar.",
                cantidad_esperada(esperadas),
                miles(observadas),
                una_en(probabilidad)
            ));
            return;
        }
    } else if cantidad < esperadas {
        let probabilidad = poisson_como_mucho(esperadas, observadas);
        if probabilidad < SIGNIFICANCIA {
            salida.explicar(format!(
                "{modelo} esperarías unas {} apariciones y hay solo {}: aparece menos de lo normal (la probabilidad de ver tan pocas es {}). Los genomas evitan algunas secuencias; las bacterias, por ejemplo, las que cortarían sus propias enzimas de restricción.",
                cantidad_esperada(esperadas),
                miles(observadas),
                una_en(probabilidad)
            ));
            return;
        }
    }
    salida.explicar(format!(
        "{modelo} esperarías unas {} apariciones; ver {} es compatible con la casualidad.",
        cantidad_esperada(esperadas),
        miles(observadas)
    ));
}

fn explicar_palindromo(salida: &mut Salida, patron: &[u8]) {
    salida.explicar(format!(
        "{} es un palíndromo inverso: es igual a su complemento inverso, así que se lee igual en las dos hebras. Muchas proteínas que trabajan en pareja, como las enzimas de restricción, reconocen secuencias así.",
        letras(patron)
    ));
}

fn explicar_senales_conocidas(salida: &mut Salida, patron: &[u8]) {
    if patron == b"GATC" {
        salida.explicar(
            "GATC es el sitio que metila la enzima Dam en E. coli, y oriC tiene muchos más de lo normal: mientras el origen recién copiado está metilado en una sola hebra, la proteína SeqA lo bloquea y evita que la replicación vuelva a arrancar demasiado pronto.",
        );
    } else if diferencias_con_dnaa(patron).is_some_and(|diferencias| diferencias <= 1) {
        salida.explicar(
            "Se parece a la caja DnaA de E. coli (TTATCCACA, o TGTGGATAA en la otra hebra): la proteína DnaA se une a estas cajas para abrir el ADN en el origen.",
        );
    }
}

fn leer_patron(entrada: &Entrada, campo: &str, maximo: usize) -> Result<Vec<u8>, String> {
    let patron = entrada.adn(campo)?;
    if patron.len() > maximo {
        return Err(format!(
            "El patrón de «{campo}» tiene {} letras; aquí se admiten hasta {}.",
            miles(patron.len()),
            miles(maximo)
        ));
    }
    Ok(patron)
}

fn leer_k(entrada: &Entrada, campo: &str, longitud: usize, maximo: usize) -> Result<usize, String> {
    let k = entrada.natural(campo)?;
    if k == 0 {
        return Err(format!("«{campo}» debe ser al menos 1."));
    }
    if k > maximo {
        return Err(format!(
            "«{campo}» puede valer como mucho {maximo} en este algoritmo, y vale {}.",
            miles(k)
        ));
    }
    if k > longitud {
        return Err(format!(
            "«{campo}» = {k} es más largo que la secuencia ({} letras).",
            miles(longitud)
        ));
    }
    Ok(k)
}

fn leer_d(entrada: &Entrada, campo: &str, k: usize, maximo: usize) -> Result<usize, String> {
    let d = entrada.natural(campo)?;
    if d > k {
        return Err(format!(
            "«{campo}» = {} no puede ser mayor que el largo del patrón ({k}).",
            miles(d)
        ));
    }
    if d > maximo {
        return Err(format!(
            "«{campo}» puede valer como mucho {maximo} aquí: con más errores el número de vecinos explota."
        ));
    }
    Ok(d)
}

fn revisar_largo(patron: &[u8], texto: &[u8]) -> Result<(), String> {
    if patron.len() > texto.len() {
        return Err(format!(
            "El patrón ({} letras) es más largo que el texto ({}).",
            miles(patron.len()),
            miles(texto.len())
        ));
    }
    Ok(())
}

fn inicio_en_ecoli(entrada: &Entrada, campo: &str) -> Option<usize> {
    let texto = entrada.texto(campo).ok()?;
    let texto = texto.trim();
    if texto == "@ecoli" {
        return Some(1);
    }
    let (inicio, _) = texto.strip_prefix("@ecoli:")?.split_once('-')?;
    inicio.trim().parse().ok()
}

fn contiene_oric(inicio: usize, longitud: usize) -> bool {
    inicio <= ORIC_INICIO && inicio + longitud > ORIC_FIN
}

fn distancia_a_oric(coordenada: usize) -> usize {
    if coordenada < ORIC_INICIO {
        ORIC_INICIO - coordenada
    } else {
        coordenada.saturating_sub(ORIC_FIN)
    }
}

fn diferencias_entre(uno: &[u8], otro: &[u8]) -> usize {
    uno.iter().zip(otro).filter(|(a, b)| a != b).count()
}

fn unica_diferencia(uno: &[u8], otro: &[u8]) -> Option<usize> {
    let mut encontrada = None;
    for (posicion, (a, b)) in uno.iter().zip(otro).enumerate() {
        if a != b {
            if encontrada.is_some() {
                return None;
            }
            encontrada = Some(posicion);
        }
    }
    encontrada
}

fn diferencias_con_dnaa(patron: &[u8]) -> Option<usize> {
    if patron.len() != CAJA_DNAA.len() {
        return None;
    }
    let inversa = secuencia::complemento_inverso(CAJA_DNAA);
    Some(diferencias_entre(patron, CAJA_DNAA).min(diferencias_entre(patron, &inversa)))
}

fn tipo_de_cambio(a: u8, b: u8) -> &'static str {
    let purina = |base: u8| matches!(base, b'A' | b'G');
    let pirimidina = |base: u8| matches!(base, b'C' | b'T');
    if (purina(a) && purina(b)) || (pirimidina(a) && pirimidina(b)) {
        "transición"
    } else if (purina(a) && pirimidina(b)) || (pirimidina(a) && purina(b)) {
        "transversión"
    } else {
        "otra letra"
    }
}

fn es_repeticion_simple(kmero: &[u8]) -> bool {
    kmero.len() >= 4 && (2..kmero.len()).all(|i| kmero[i] == kmero[i - 2])
}

fn primer_solapamiento(posiciones: &[usize], largo: usize) -> Option<(usize, usize)> {
    posiciones
        .windows(2)
        .find(|par| par[1] - par[0] < largo)
        .map(|par| (par[0], par[1]))
}

fn mas_juntas(posiciones: &[usize], alcance: usize) -> (usize, usize) {
    let mut mejor = (0, 0);
    if alcance == 0 {
        return mejor;
    }
    let mut izquierda = 0;
    for derecha in 0..posiciones.len() {
        while posiciones[derecha] - posiciones[izquierda] >= alcance {
            izquierda += 1;
        }
        let juntas = derecha - izquierda + 1;
        if juntas > mejor.0 {
            mejor = (juntas, posiciones[izquierda]);
        }
    }
    mejor
}

fn proporciones(adn: &[u8]) -> [f64; 4] {
    let mut cuentas = [0usize; 4];
    for &base in adn {
        if let Some(valor) = replicacion::valor_de_base(base) {
            cuentas[valor as usize] += 1;
        }
    }
    let total: usize = cuentas.iter().sum();
    if total == 0 {
        return [0.25; 4];
    }
    cuentas.map(|cuenta| cuenta as f64 / total as f64)
}

fn probabilidad_del_patron(patron: &[u8], proporcion: &[f64; 4]) -> f64 {
    patron
        .iter()
        .map(|&base| {
            replicacion::valor_de_base(base)
                .map(|valor| proporcion[valor as usize])
                .unwrap_or(0.0)
        })
        .product()
}

fn esperadas_por_composicion(texto: &[u8], patron: &[u8]) -> f64 {
    if patron.len() > texto.len() {
        return 0.0;
    }
    let ventanas = (texto.len() - patron.len() + 1) as f64;
    ventanas * probabilidad_del_patron(patron, &proporciones(texto))
}

fn probabilidad_con_errores(k: usize, d: usize) -> f64 {
    let mut total = 0.0;
    let mut combinaciones = 1.0;
    for cambios in 0..=d.min(k) {
        if cambios > 0 {
            combinaciones = combinaciones * (k - cambios + 1) as f64 / cambios as f64;
        }
        total += combinaciones * 0.75f64.powi(cambios as i32) * 0.25f64.powi((k - cambios) as i32);
    }
    total.min(1.0)
}

fn logaritmo_factorial(n: usize) -> f64 {
    if n < 1000 {
        return (2..=n).map(|i| (i as f64).ln()).sum();
    }
    let x = n as f64;
    x * x.ln() - x + 0.5 * (2.0 * std::f64::consts::PI * x).ln() + 1.0 / (12.0 * x)
}

fn termino_de_poisson(media: f64, veces: usize) -> f64 {
    (-media + veces as f64 * media.ln() - logaritmo_factorial(veces)).exp()
}

fn sumar_hacia_arriba(media: f64, desde: usize) -> f64 {
    let mut termino = termino_de_poisson(media, desde);
    let mut suma = 0.0;
    let mut veces = desde;
    loop {
        suma += termino;
        veces += 1;
        termino *= media / veces as f64;
        if termino <= suma * 1e-16 || veces > desde + 10_000_000 {
            return suma;
        }
    }
}

fn sumar_hacia_abajo(media: f64, hasta: usize) -> f64 {
    let mut termino = termino_de_poisson(media, hasta);
    let mut suma = 0.0;
    let mut veces = hasta;
    loop {
        suma += termino;
        if veces == 0 || termino <= suma * 1e-16 {
            return suma;
        }
        termino *= veces as f64 / media;
        veces -= 1;
    }
}

fn poisson_al_menos(media: f64, veces: usize) -> f64 {
    if veces == 0 {
        return 1.0;
    }
    if media <= 0.0 {
        return 0.0;
    }
    if veces as f64 > media {
        sumar_hacia_arriba(media, veces).min(1.0)
    } else {
        (1.0 - sumar_hacia_abajo(media, veces - 1)).max(0.0)
    }
}

fn poisson_como_mucho(media: f64, veces: usize) -> f64 {
    if media <= 0.0 {
        return 1.0;
    }
    if (veces as f64) < media {
        sumar_hacia_abajo(media, veces).min(1.0)
    } else {
        (1.0 - sumar_hacia_arriba(media, veces + 1)).max(0.0)
    }
}

fn una_en(probabilidad: f64) -> String {
    if probabilidad < 1e-6 {
        "menor que una en un millón".to_string()
    } else {
        format!(
            "de alrededor de una en {}",
            miles((1.0 / probabilidad).round() as usize)
        )
    }
}

fn cantidad_esperada(valor: f64) -> String {
    if valor >= 10.0 {
        miles(valor.round() as usize)
    } else if valor >= 1.0 {
        decimal(valor, 1)
    } else if valor >= 0.01 {
        decimal(valor, 2)
    } else {
        cientifico(valor)
    }
}

fn con_signo(numero: i64) -> String {
    if numero < 0 {
        format!("−{}", miles(numero.unsigned_abs() as usize))
    } else {
        miles(numero as usize)
    }
}

fn plural(cantidad: usize, singular: &str, plural: &str) -> String {
    if cantidad == 1 {
        format!("1 {singular}")
    } else {
        format!("{} {plural}", miles(cantidad))
    }
}

fn letras(secuencia: &[u8]) -> String {
    String::from_utf8_lossy(secuencia).into_owned()
}

fn recortar(secuencia: &[u8]) -> String {
    if secuencia.len() <= MAXIMO_TEXTO {
        letras(secuencia)
    } else {
        format!("{}…", letras(&secuencia[..MAXIMO_TEXTO]))
    }
}

fn unir(partes: impl Iterator<Item = String>) -> String {
    let partes: Vec<String> = partes.collect();
    match partes.len() {
        0 => String::new(),
        1 => partes[0].clone(),
        cantidad => format!(
            "{} y {}",
            partes[..cantidad - 1].join(", "),
            partes[cantidad - 1]
        ),
    }
}

fn tabla_con_total(
    etiqueta: &str,
    columnas: &[&str],
    filas: Vec<Vec<String>>,
    total: usize,
) -> Bloque {
    Bloque::Tabla {
        etiqueta: etiqueta.to_string(),
        columnas: columnas.iter().map(|columna| columna.to_string()).collect(),
        filas,
        total,
    }
}

fn contar_por_tramos(
    posiciones: impl Iterator<Item = usize>,
    longitud: usize,
) -> (Vec<String>, Vec<f64>) {
    let tramos = TRAMOS.min(longitud.max(1));
    let ancho = longitud.div_ceil(tramos).max(1);
    let mut valores = vec![0.0; tramos];
    for posicion in posiciones {
        valores[(posicion / ancho).min(tramos - 1)] += 1.0;
    }
    let categorias = (0..tramos)
        .map(|tramo| {
            if longitud >= 20_000 {
                format!("{} kb", miles(tramo * ancho / 1000))
            } else {
                miles(tramo * ancho)
            }
        })
        .collect();
    (categorias, valores)
}

fn histograma_de_veces(cuentas: impl Iterator<Item = usize>) -> (Vec<String>, Vec<f64>) {
    const TOPE: usize = 15;
    let mut valores = vec![0.0; TOPE];
    let mut mayor = 0;
    for veces in cuentas {
        let casilla = veces.min(TOPE) - 1;
        valores[casilla] += 1.0;
        mayor = mayor.max(casilla + 1);
    }
    valores.truncate(mayor);
    let categorias = (1..=valores.len())
        .map(|veces| {
            if veces == TOPE {
                format!("{TOPE} o más")
            } else {
                veces.to_string()
            }
        })
        .collect();
    (categorias, valores)
}

fn mejores_de_la_tabla<'a>(
    tabla: &HashMap<&'a [u8], usize>,
    cantidad: usize,
) -> Vec<(&'a [u8], usize)> {
    let mut mejores: Vec<(&[u8], usize)> = Vec::with_capacity(cantidad + 1);
    for (&kmero, &veces) in tabla {
        let va_antes =
            |otro: &(&[u8], usize)| veces > otro.1 || (veces == otro.1 && kmero < otro.0);
        if mejores.len() == cantidad && !va_antes(&mejores[cantidad - 1]) {
            continue;
        }
        let lugar = mejores.iter().position(va_antes).unwrap_or(mejores.len());
        mejores.insert(lugar, (kmero, veces));
        mejores.truncate(cantidad);
    }
    mejores
}

fn mejores_del_arreglo(cuentas: &[u32], cantidad: usize) -> Vec<(usize, u32)> {
    let mut mejores: Vec<(usize, u32)> = Vec::with_capacity(cantidad + 1);
    for (numero, &veces) in cuentas.iter().enumerate() {
        if veces == 0 || (mejores.len() == cantidad && veces <= mejores[cantidad - 1].1) {
            continue;
        }
        let lugar = mejores
            .iter()
            .position(|&(_, otra)| veces > otra)
            .unwrap_or(mejores.len());
        mejores.insert(lugar, (numero, veces));
        mejores.truncate(cantidad);
    }
    mejores
}

struct RecorridoDelSesgo {
    puntos: Vec<[f64; 2]>,
    minimo: i64,
    maximo: i64,
    posicion_del_maximo: usize,
}

fn recorrer_sesgo(genoma: &[u8]) -> RecorridoDelSesgo {
    let por_grupo = (genoma.len() + 1).div_ceil(PUNTOS_DE_CURVA / 2).max(1);
    let mut puntos = Vec::with_capacity(PUNTOS_DE_CURVA + 2);
    let mut actual: i64 = 0;
    let mut minimo: i64 = 0;
    let mut maximo: i64 = 0;
    let mut posicion_del_maximo = 0;
    let mut menor = (0, 0);
    let mut mayor = (0, 0);
    let mut en_el_grupo = 0;
    for posicion in 0..=genoma.len() {
        if posicion > 0 {
            actual += replicacion::cambio_de_sesgo(genoma[posicion - 1]);
        }
        minimo = minimo.min(actual);
        if actual > maximo {
            maximo = actual;
            posicion_del_maximo = posicion;
        }
        if en_el_grupo == 0 || actual < menor.1 {
            menor = (posicion, actual);
        }
        if en_el_grupo == 0 || actual > mayor.1 {
            mayor = (posicion, actual);
        }
        en_el_grupo += 1;
        if en_el_grupo == por_grupo || posicion == genoma.len() {
            let (primero, segundo) = if menor.0 <= mayor.0 {
                (menor, mayor)
            } else {
                (mayor, menor)
            };
            puntos.push([primero.0 as f64, primero.1 as f64]);
            if segundo.0 != primero.0 {
                puntos.push([segundo.0 as f64, segundo.1 as f64]);
            }
            en_el_grupo = 0;
        }
    }
    RecorridoDelSesgo {
        puntos,
        minimo,
        maximo,
        posicion_del_maximo,
    }
}

fn muestra_pareja(puntos: &[[i64; 2]], maximo: usize) -> Vec<[i64; 2]> {
    let salto = puntos.len().div_ceil(maximo.max(1)).max(1);
    let mut elegidos: Vec<[i64; 2]> = puntos.iter().step_by(salto).copied().collect();
    if let Some(&ultimo) = puntos.last()
        && elegidos.last() != Some(&ultimo)
    {
        elegidos.push(ultimo);
    }
    elegidos
}

const DETALLE: usize = 6;
const DETALLE_APROXIMADO: usize = 4;
const MAXIMO_HALLAZGOS: usize = 40;
const MAXIMO_EVENTOS: usize = 14;
const CONTROLES: usize = 5;
const FILAS_DEL_PASO: usize = 16;

fn revisar_animable(longitud: usize) -> Result<(), String> {
    if longitud > MAXIMO_ANIMABLE {
        return Err(format!(
            "El paso a paso admite hasta {} bases y esta secuencia tiene {}.",
            miles(MAXIMO_ANIMABLE),
            miles(longitud)
        ));
    }
    Ok(())
}

fn vista(secuencia: &[u8]) -> String {
    if secuencia.len() <= 20 {
        letras(secuencia)
    } else {
        format!(
            "{}… ({} letras)",
            letras(&secuencia[..12]),
            miles(secuencia.len())
        )
    }
}

fn vista_corta(secuencia: &[u8]) -> String {
    if secuencia.len() <= 14 {
        letras(secuencia)
    } else {
        format!("{}…", letras(&secuencia[..12]))
    }
}

fn lista_corta(valores: &[usize]) -> String {
    let mut partes: Vec<String> = valores.iter().take(6).map(|&valor| miles(valor)).collect();
    if valores.len() > 6 {
        partes.push(format!("… ({} en total)", miles(valores.len())));
    }
    format!("[{}]", partes.join(", "))
}

fn puntos_de_control(total: usize) -> Vec<usize> {
    (1..=CONTROLES)
        .map(|vuelta| total * vuelta / (CONTROLES + 1))
        .collect()
}

fn letra_de(base: u8) -> String {
    (base as char).to_string()
}

struct Barrido<'a> {
    texto: &'a [u8],
    patron: &'a [u8],
    d: Option<usize>,
    acumulado: &'static str,
    con_posicion: bool,
    inicio: &'static str,
    bucle: &'static str,
    prueba: &'static str,
    acierto: &'static str,
}

impl Barrido<'_> {
    fn coincide(&self, ventana: &[u8]) -> bool {
        match self.d {
            Some(d) => replicacion::casi_iguales(ventana, self.patron, d),
            None => ventana == self.patron,
        }
    }

    fn acumulado_texto(&self, hallazgos: &[usize]) -> String {
        if self.con_posicion {
            lista_corta(hallazgos)
        } else {
            miles(hallazgos.len())
        }
    }

    fn pizarra(
        &self,
        indice: Option<usize>,
        diferencias: Option<usize>,
        hallazgos: &[usize],
        marcas: &[(usize, usize, Resalte)],
    ) -> Pizarra {
        let mut pizarra = Pizarra::nueva().variable("patron", vista_corta(self.patron));
        if let Some(d) = self.d {
            pizarra = pizarra.variable("d", d);
        }
        if let Some(indice) = indice {
            if self.con_posicion {
                pizarra = pizarra.variable("posicion", miles(indice));
            }
            let ventana = &self.texto[indice..indice + self.patron.len()];
            pizarra = pizarra.variable("ventana", vista_corta(ventana));
        }
        if let Some(diferencias) = diferencias {
            pizarra = pizarra.variable("diferencias", diferencias);
        }
        pizarra
            .variable(self.acumulado, self.acumulado_texto(hallazgos))
            .secuencia(self.texto, marcas)
    }
}

type LetraPorLetra = (Vec<Vec<String>>, Vec<(usize, usize)>, usize, bool);

fn comparar_letras(b: &Barrido, indice: usize) -> LetraPorLetra {
    let d = b.d.unwrap_or(0);
    let ventana = &b.texto[indice..indice + b.patron.len()];
    let mut filas = Vec::new();
    let mut resaltar = Vec::new();
    let mut diferencias = 0;
    let mut aceptada = true;
    for (lugar, (a, c)) in ventana.iter().zip(b.patron).enumerate() {
        let distinta = a != c;
        if filas.len() < FILAS_DEL_PASO {
            if distinta {
                resaltar.push((filas.len(), 3));
            }
            filas.push(vec![
                (lugar + 1).to_string(),
                letra_de(*c),
                letra_de(*a),
                if distinta { "difiere" } else { "igual" }.to_string(),
            ]);
        }
        if distinta {
            diferencias += 1;
            if diferencias > d {
                aceptada = false;
                break;
            }
        }
    }
    (filas, resaltar, diferencias, aceptada)
}

fn barrer(b: &Barrido) -> (Vec<Paso>, Vec<usize>) {
    let k = b.patron.len();
    let ventanas = if k <= b.texto.len() {
        b.texto.len() - k + 1
    } else {
        0
    };
    let detalle = if b.d.is_some() {
        DETALLE_APROXIMADO
    } else {
        DETALLE
    };
    let controles = puntos_de_control(ventanas);
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let mut hallazgos: Vec<usize> = Vec::new();
    let mut avisado = false;
    let patron_texto = letras(b.patron);
    let objetivo = match (b.con_posicion, b.d.is_some()) {
        (false, false) => format!("cuántas veces aparece «{}»", vista(b.patron)),
        (true, false) => format!("en qué posiciones aparece «{}»", vista(b.patron)),
        (false, true) => format!("cuántas ventanas se parecen a «{}»", vista(b.patron)),
        (true, true) => format!(
            "en qué posiciones hay algo parecido a «{}»",
            vista(b.patron)
        ),
    };
    let tolerancia = match b.d {
        Some(d) => format!(
            " Una ventana cuenta si difiere del patrón en {} como mucho.",
            plural(d, "letra", "letras")
        ),
        None => String::new(),
    };
    let cabida = if ventanas == 0 {
        "Como el patrón es más largo que el texto, no cabe ninguna ventana y no hay nada que hallar."
            .to_string()
    } else {
        format!(
            "Su ventana de {k} letras se puede colocar en {} sitios, desde la posición 0 hasta la {}.",
            miles(ventanas),
            miles(ventanas - 1)
        )
    };
    let apertura = if b.con_posicion {
        "Empezamos con la lista `posiciones` vacía."
    } else {
        "Empezamos con el contador `veces` en 0."
    };
    narracion.agregar(b.pizarra(None, None, &hallazgos, &[]).paso(
        b.inicio,
        format!(
            "Queremos saber {objetivo} dentro de un texto de {} letras.{tolerancia} {cabida} {apertura}",
            miles(b.texto.len())
        ),
    ));
    for indice in 0..ventanas {
        let ventana = &b.texto[indice..indice + k];
        let coincide = b.coincide(ventana);
        let actual = [(indice, indice + k, Resalte::Actual)];
        let hallada = [(indice, indice + k, Resalte::Hallado)];
        if indice == 0 {
            narracion.agregar(b.pizarra(Some(0), None, &hallazgos, &actual).paso(
                b.bucle,
                format!(
                    "El bucle toma una ventana: un trozo del texto de {k} letras. La primera empieza en la posición 0 y es «{}»; en cada vuelta se corre una letra a la derecha.",
                    letras(ventana)
                ),
            ));
        }
        if b.d.is_none() && indice < detalle {
            let (marcas, nota) = if coincide {
                (
                    &hallada,
                    format!(
                        "Ventana {} = «{}». ¿Es igual a «{patron_texto}»? Sí, letra por letra.",
                        miles(indice),
                        letras(ventana)
                    ),
                )
            } else {
                let lugar = ventana
                    .iter()
                    .zip(b.patron)
                    .position(|(a, c)| a != c)
                    .unwrap_or(0);
                (
                    &actual,
                    format!(
                        "Ventana {} = «{}». ¿Es igual a «{patron_texto}»? No: en la letra {} hay {} y el patrón pide {}. Pasamos a la siguiente.",
                        miles(indice),
                        letras(ventana),
                        lugar + 1,
                        letra_de(ventana[lugar]),
                        letra_de(b.patron[lugar])
                    ),
                )
            };
            narracion.agregar(
                b.pizarra(Some(indice), None, &hallazgos, marcas)
                    .paso(b.prueba, nota),
            );
        }
        if b.d.is_some() && indice < detalle {
            let d = b.d.unwrap_or(0);
            let (filas, resaltar, diferencias, aceptada) = comparar_letras(b, indice);
            let nota = if aceptada {
                format!(
                    "Comparamos la ventana «{}» con «{patron_texto}» letra a letra: hay {}.",
                    letras(ventana),
                    plural(diferencias, "diferencia", "diferencias")
                )
            } else {
                format!(
                    "Comparamos «{}» con «{patron_texto}» letra a letra: ya van {diferencias} diferencias y solo se permiten {d}.",
                    letras(ventana)
                )
            };
            narracion.agregar(
                b.pizarra(Some(indice), Some(diferencias), &hallazgos, &actual)
                    .tabla(
                        "Letra por letra",
                        &["letra", "patrón", "ventana", "resultado"],
                        filas,
                        &resaltar,
                    )
                    .paso("for (a, b) in ventana.iter().zip(patron) {", nota),
            );
            if aceptada {
                narracion.agregar(
                    b.pizarra(Some(indice), Some(diferencias), &hallazgos, &hallada)
                        .paso(
                            b.prueba,
                            format!(
                                "Con {} (se permiten hasta {d}) la ventana cuenta como coincidencia.",
                                plural(diferencias, "diferencia", "diferencias")
                            ),
                        ),
                );
            } else {
                narracion.agregar(
                    b.pizarra(Some(indice), Some(diferencias), &hallazgos, &actual)
                        .paso(
                            "return false;",
                            format!(
                                "Ya pasamos de las {d} diferencias permitidas: la comparación se corta aquí, sin mirar las letras que faltan, y la ventana no cuenta."
                            ),
                        ),
                );
            }
        }
        if coincide {
            hallazgos.push(indice);
            if hallazgos.len() <= MAXIMO_HALLAZGOS {
                let parecido = match b.d {
                    Some(_) => format!(
                        " Difiere del patrón en {}.",
                        plural(diferencias_entre(ventana, b.patron), "letra", "letras")
                    ),
                    None => String::new(),
                };
                let nota = if b.con_posicion {
                    format!(
                        "Coincide: guardamos la posición {} en la lista, que ahora tiene {}.{parecido}",
                        miles(indice),
                        plural(hallazgos.len(), "posición", "posiciones")
                    )
                } else {
                    format!(
                        "Coincide: `veces` sube a {}. Esta aparición empieza en la posición {}.{parecido}",
                        miles(hallazgos.len()),
                        miles(indice)
                    )
                };
                narracion.agregar(
                    b.pizarra(Some(indice), None, &hallazgos, &hallada)
                        .paso(b.acierto, nota),
                );
            } else if !avisado {
                avisado = true;
                narracion.agregar(b.pizarra(Some(indice), None, &hallazgos, &hallada).paso(
                    b.bucle,
                    "Ya viste cómo se anota cada coincidencia. Para no cansarte, seguimos en silencio y te aviso de cuánto llevamos en unos cuantos puntos del recorrido.",
                ));
            }
        }
        if indice >= detalle && controles.contains(&indice) {
            narracion.agregar(b.pizarra(Some(indice), None, &hallazgos, &actual).paso(
                b.bucle,
                format!(
                    "Vamos por la ventana {} de {}: hasta aquí hay {}.",
                    miles(indice),
                    miles(ventanas),
                    if b.con_posicion {
                        plural(hallazgos.len(), "posición hallada", "posiciones halladas")
                    } else {
                        plural(hallazgos.len(), "aparición", "apariciones")
                    }
                ),
            ));
        }
    }
    let marcas: Vec<(usize, usize, Resalte)> = hallazgos
        .iter()
        .take(30)
        .map(|&posicion| (posicion, posicion + k, Resalte::Hallado))
        .collect();
    let resultado = if b.con_posicion {
        if hallazgos.is_empty() {
            format!(
                "Terminamos: miramos las {} ventanas y ninguna coincide, así que la lista queda vacía.",
                miles(ventanas)
            )
        } else {
            format!(
                "Terminamos: miramos las {} ventanas y hay {}: {}. Son posiciones contando desde 0.",
                miles(ventanas),
                plural(hallazgos.len(), "coincidencia", "coincidencias"),
                lista_corta(&hallazgos)
            )
        }
    } else {
        format!(
            "Terminamos: miramos las {} ventanas y el resultado es {}. En «Qué significa» verás si es mucho o poco para este texto.",
            miles(ventanas),
            miles(hallazgos.len())
        )
    };
    let ultimo = b
        .pizarra(None, None, &hallazgos, &marcas)
        .paso(b.bucle, resultado);
    (narracion.terminar(ultimo), hallazgos)
}

fn preparar_exacto(entrada: &Entrada) -> Result<(Vec<u8>, Vec<u8>), String> {
    let texto = entrada.secuencia("texto")?;
    let patron = leer_patron(entrada, "patron", MAXIMO_PATRON_EXACTO)?;
    revisar_animable(texto.len())?;
    Ok((texto, patron))
}

fn preparar_aproximado(entrada: &Entrada) -> Result<(Vec<u8>, Vec<u8>, usize), String> {
    let texto = entrada.secuencia("texto")?;
    let patron = leer_patron(entrada, "patron", MAXIMO_PATRON_APROXIMADO)?;
    let d = leer_d(entrada, "d", patron.len(), patron.len())?;
    revisar_largo(&patron, &texto)?;
    revisar_animable(texto.len())?;
    Ok((texto, patron, d))
}

fn barrido_de_contar<'a>(texto: &'a [u8], patron: &'a [u8]) -> Barrido<'a> {
    Barrido {
        texto,
        patron,
        d: None,
        acumulado: "veces",
        con_posicion: false,
        inicio: "let mut veces = 0;",
        bucle: "for ventana in texto.windows(patron.len()) {",
        prueba: "if ventana == patron {",
        acierto: "veces += 1;",
    }
}

fn barrido_de_posiciones<'a>(texto: &'a [u8], patron: &'a [u8]) -> Barrido<'a> {
    Barrido {
        texto,
        patron,
        d: None,
        acumulado: "posiciones",
        con_posicion: true,
        inicio: "let mut posiciones = Vec::new();",
        bucle: "for (posicion, ventana) in texto.windows(patron.len()).enumerate() {",
        prueba: "if ventana == patron {",
        acierto: "posiciones.push(posicion);",
    }
}

fn barrido_de_contar_aproximado<'a>(texto: &'a [u8], patron: &'a [u8], d: usize) -> Barrido<'a> {
    Barrido {
        texto,
        patron,
        d: Some(d),
        acumulado: "veces",
        con_posicion: false,
        inicio: "let mut veces = 0;",
        bucle: "for ventana in texto.windows(patron.len()) {",
        prueba: "if casi_iguales(ventana, patron, d) {",
        acierto: "veces += 1;",
    }
}

fn barrido_de_posiciones_aproximadas<'a>(
    texto: &'a [u8],
    patron: &'a [u8],
    d: usize,
) -> Barrido<'a> {
    Barrido {
        texto,
        patron,
        d: Some(d),
        acumulado: "posiciones",
        con_posicion: true,
        inicio: "let mut posiciones = Vec::new();",
        bucle: "for (posicion, ventana) in texto.windows(patron.len()).enumerate() {",
        prueba: "if casi_iguales(ventana, patron, d) {",
        acierto: "posiciones.push(posicion);",
    }
}

fn pasos_contar_patron(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (texto, patron) = preparar_exacto(entrada)?;
    Ok(barrer(&barrido_de_contar(&texto, &patron)).0)
}

fn pasos_posiciones_del_patron(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (texto, patron) = preparar_exacto(entrada)?;
    Ok(barrer(&barrido_de_posiciones(&texto, &patron)).0)
}

fn pasos_contar_aproximado(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (texto, patron, d) = preparar_aproximado(entrada)?;
    Ok(barrer(&barrido_de_contar_aproximado(&texto, &patron, d)).0)
}

fn pasos_posiciones_aproximadas(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (texto, patron, d) = preparar_aproximado(entrada)?;
    Ok(barrer(&barrido_de_posiciones_aproximadas(&texto, &patron, d)).0)
}

fn trazar_complemento_inverso(entrada: &Entrada) -> Result<(Vec<Paso>, Vec<u8>), String> {
    let adn = entrada.secuencia("adn")?;
    revisar_animable(adn.len())?;
    let longitud = adn.len();
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        Pizarra::nueva()
            .variable("adn", vista(&adn))
            .secuencia(&adn, &[])
            .paso(
                "pub fn complemento_inverso(adn: &[u8]) -> Vec<u8> {",
                format!(
                    "Queremos la hebra de enfrente de «{}», leída de 5' a 3'. Son dos movimientos: leer la secuencia desde el final hacia el principio y cambiar cada letra por su pareja (A con T, C con G).",
                    vista(&adn)
                ),
            ),
    );
    let al_reves: Vec<u8> = adn.iter().rev().copied().collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("adn", vista(&adn))
            .variable("adn.iter().rev()", vista(&al_reves))
            .secuencia(&adn, &[(longitud - 1, longitud, Resalte::Actual)])
            .paso(
                "adn.iter().rev().map(|&base| complemento(base)).collect()",
                format!(
                    "`.rev()` recorre la secuencia al revés: primero la última letra ({}), luego la penúltima, y así hasta la primera. Leída al revés queda «{}».",
                    letra_de(adn[longitud - 1]),
                    vista(&al_reves)
                ),
            ),
    );
    let mut resultado: Vec<u8> = Vec::with_capacity(longitud);
    let controles = puntos_de_control(longitud);
    for (vuelta, &base) in adn.iter().rev().enumerate() {
        let pareja = nucleo::secuencia::complemento(base);
        resultado.push(pareja);
        let lugar = longitud - 1 - vuelta;
        if vuelta < 10 {
            let (linea, explicacion) = match base {
                b'A' => ("b'A' => b'T',", "A se une con T."),
                b'T' => ("b'T' => b'A',", "T se une con A."),
                b'U' => ("b'U' => b'A',", "U (el ARN) se une con A."),
                b'C' => ("b'C' => b'G',", "C se une con G."),
                b'G' => ("b'G' => b'C',", "G se une con C."),
                b'a' => ("b'a' => b't',", "a minúscula se une con t minúscula."),
                b't' => ("b't' => b'a',", "t minúscula se une con a minúscula."),
                b'u' => ("b'u' => b'a',", "u minúscula se une con a minúscula."),
                b'c' => ("b'c' => b'g',", "c minúscula se une con g minúscula."),
                b'g' => ("b'g' => b'c',", "g minúscula se une con c minúscula."),
                _ => (
                    "otra => otra,",
                    "No es una base conocida (como N): se deja igual porque no tiene pareja definida.",
                ),
            };
            narracion.agregar(
                Pizarra::nueva()
                    .variable("base", letra_de(base))
                    .variable("resultado", vista(&resultado))
                    .secuencia(&adn, &[(lugar, lugar + 1, Resalte::Actual)])
                    .paso(
                        linea,
                        format!(
                            "Letra {} desde el final, en la posición {}: «{}». {explicacion} La agregamos al resultado, que ahora es «{}».",
                            vuelta + 1,
                            miles(lugar),
                            letra_de(base),
                            vista(&resultado)
                        ),
                    ),
            );
        } else if controles.contains(&vuelta) {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("base", letra_de(base))
                    .variable("resultado", vista(&resultado))
                    .secuencia(&adn, &[(lugar, lugar + 1, Resalte::Actual)])
                    .paso(
                        "adn.iter().rev().map(|&base| complemento(base)).collect()",
                        format!(
                            "Seguimos igual con cada letra: llevamos {} de {}.",
                            miles(vuelta + 1),
                            miles(longitud)
                        ),
                    ),
            );
        }
    }
    let ultimo = Pizarra::nueva()
        .variable("resultado", vista(&resultado))
        .secuencia(&adn, &[])
        .paso(
            "adn.iter().rev().map(|&base| complemento(base)).collect()",
            format!(
                "Listo: el complemento inverso es «{}», de {} letras. Es la otra hebra leída de 5' a 3'; en «Qué significa» ves cómo se relaciona con la tuya.",
                vista(&resultado),
                miles(resultado.len())
            ),
        );
    Ok((narracion.terminar(ultimo), resultado))
}

fn pasos_complemento_inverso(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    Ok(trazar_complemento_inverso(entrada)?.0)
}

fn trazar_contenido_gc(entrada: &Entrada) -> Result<(Vec<Paso>, f64), String> {
    let adn = entrada.secuencia("adn")?;
    revisar_animable(adn.len())?;
    let longitud = adn.len();
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        Pizarra::nueva()
            .variable("adn", vista(&adn))
            .variable("adn.len()", miles(longitud))
            .secuencia(&adn, &[])
            .paso(
                "if adn.is_empty() {",
                format!(
                    "Queremos saber qué fracción de las {} letras son G o C. La secuencia no está vacía, así que podemos seguir.",
                    miles(longitud)
                ),
            ),
    );
    let mut gc = 0usize;
    let controles = puntos_de_control(longitud);
    for (lugar, &base) in adn.iter().enumerate() {
        let cuenta = matches!(base, b'G' | b'C' | b'g' | b'c');
        if cuenta {
            gc += 1;
        }
        let marca = if cuenta {
            Resalte::Hallado
        } else {
            Resalte::Actual
        };
        if lugar < 10 {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("base", letra_de(base))
                    .variable("gc", gc)
                    .secuencia(&adn, &[(lugar, lugar + 1, marca)])
                    .paso(
                        "matches!(base, b'G' | b'C' | b'g' | b'c')",
                        if cuenta {
                            format!(
                                "La letra {} es «{}»: es G o C, así que `filter` la deja pasar y el conteo sube a {}.",
                                lugar + 1,
                                letra_de(base),
                                gc
                            )
                        } else {
                            format!(
                                "La letra {} es «{}»: no es G ni C, así que `filter` la descarta y el conteo se queda en {}.",
                                lugar + 1,
                                letra_de(base),
                                gc
                            )
                        },
                    ),
            );
        } else if controles.contains(&lugar) {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("base", letra_de(base))
                    .variable("gc", miles(gc))
                    .secuencia(&adn, &[(lugar, lugar + 1, marca)])
                    .paso(
                        "matches!(base, b'G' | b'C' | b'g' | b'c')",
                        format!(
                            "Vamos por la letra {} de {}: han pasado {} letras G o C.",
                            miles(lugar + 1),
                            miles(longitud),
                            miles(gc)
                        ),
                    ),
            );
        }
    }
    let fraccion = gc as f64 / longitud as f64;
    narracion.agregar(
        Pizarra::nueva()
            .variable("gc", miles(gc))
            .variable("adn.len()", miles(longitud))
            .paso(
                "let gc = adn",
                format!(
                    "`.count()` cuenta lo que dejó pasar el filtro: {} letras G o C.",
                    miles(gc)
                ),
            ),
    );
    let ultimo = Pizarra::nueva()
        .variable("gc", miles(gc))
        .variable("adn.len()", miles(longitud))
        .paso(
            "gc as f64 / adn.len() as f64",
            format!(
                "Dividimos {} entre {}: {}. Ese es el contenido GC; en «Qué significa» lo comparamos con el de E. coli.",
                miles(gc),
                miles(longitud),
                porcentaje(fraccion, 1)
            ),
        );
    Ok((narracion.terminar(ultimo), fraccion))
}

fn pasos_contenido_gc(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    Ok(trazar_contenido_gc(entrada)?.0)
}

fn tabla_de_vecinas(
    primera: &[u8],
    segunda: &[u8],
    hasta: usize,
) -> (Vec<Vec<String>>, Vec<(usize, usize)>) {
    let desde = hasta.saturating_sub(5);
    let mut filas = Vec::new();
    let mut resaltar = Vec::new();
    for posicion in desde..=hasta {
        let distintas = primera[posicion] != segunda[posicion];
        if distintas {
            resaltar.push((filas.len(), 1));
            resaltar.push((filas.len(), 2));
        }
        filas.push(vec![
            miles(posicion),
            letra_de(primera[posicion]),
            letra_de(segunda[posicion]),
        ]);
    }
    (filas, resaltar)
}

fn trazar_hamming(entrada: &Entrada) -> Result<(Vec<Paso>, usize), String> {
    let primera = entrada.secuencia("primera")?;
    let segunda = entrada.secuencia("segunda")?;
    replicacion::distancia_hamming(&primera, &segunda)?;
    revisar_animable(primera.len())?;
    let longitud = primera.len();
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        Pizarra::nueva()
            .variable("primera", vista(&primera))
            .variable("segunda", vista(&segunda))
            .variable("diferencias", 0)
            .secuencia(&primera, &[])
            .paso(
                "if primera.len() != segunda.len() {",
                format!(
                    "La distancia de Hamming solo se puede medir entre secuencias del mismo largo. Las dos miden {} letras, así que seguimos: las pondremos una frente a la otra y contaremos en cuántas posiciones no coinciden.",
                    miles(longitud)
                ),
            ),
    );
    let mut diferencias = 0usize;
    let mut avisado = false;
    let controles = puntos_de_control(longitud);
    for posicion in 0..longitud {
        let (a, b) = (primera[posicion], segunda[posicion]);
        let distintas = a != b;
        let (filas, resaltar) = tabla_de_vecinas(&primera, &segunda, posicion);
        let marca = if distintas {
            Resalte::Hallado
        } else {
            Resalte::Actual
        };
        if posicion < DETALLE {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("a", letra_de(a))
                    .variable("b", letra_de(b))
                    .variable("diferencias", diferencias)
                    .secuencia(&primera, &[(posicion, posicion + 1, marca)])
                    .tabla(
                        "Últimas posiciones comparadas",
                        &["posición", "primera", "segunda"],
                        filas.clone(),
                        &resaltar,
                    )
                    .paso(
                        "if a != b {",
                        if distintas {
                            format!(
                                "Posición {posicion}: la primera tiene «{}» y la segunda «{}». Son distintas.",
                                letra_de(a),
                                letra_de(b)
                            )
                        } else {
                            format!(
                                "Posición {posicion}: las dos tienen «{}». Son iguales, así que no se cuenta nada.",
                                letra_de(a)
                            )
                        },
                    ),
            );
        }
        if distintas {
            diferencias += 1;
            if diferencias <= MAXIMO_HALLAZGOS {
                narracion.agregar(
                    Pizarra::nueva()
                        .variable("a", letra_de(a))
                        .variable("b", letra_de(b))
                        .variable("diferencias", diferencias)
                        .secuencia(&primera, &[(posicion, posicion + 1, marca)])
                        .tabla(
                            "Últimas posiciones comparadas",
                            &["posición", "primera", "segunda"],
                            filas,
                            &resaltar,
                        )
                        .paso(
                            "diferencias += 1;",
                            format!(
                                "Posición {}: «{}» frente a «{}». Es una diferencia más: llevamos {diferencias}.",
                                miles(posicion),
                                letra_de(a),
                                letra_de(b)
                            ),
                        ),
                );
            } else if !avisado {
                avisado = true;
                narracion.agregar(
                    Pizarra::nueva()
                        .variable("diferencias", diferencias)
                        .secuencia(&primera, &[(posicion, posicion + 1, marca)])
                        .paso(
                            "if a != b {",
                            "Ya viste cómo se cuenta cada diferencia. Seguimos en silencio y te aviso de cuántas llevamos en unos cuantos puntos del recorrido.",
                        ),
                );
            }
        } else if posicion >= DETALLE && controles.contains(&posicion) {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("diferencias", diferencias)
                    .secuencia(&primera, &[(posicion, posicion + 1, marca)])
                    .paso(
                        "if a != b {",
                        format!(
                            "Vamos por la posición {} de {}: hasta aquí hay {}.",
                            miles(posicion),
                            miles(longitud),
                            plural(diferencias, "diferencia", "diferencias")
                        ),
                    ),
            );
        }
    }
    let ultimo = Pizarra::nueva()
        .variable("diferencias", miles(diferencias))
        .paso(
            "Ok(diferencias)",
            format!(
                "Terminamos: las dos secuencias difieren en {} de {} posiciones. En «Qué significa» ves qué proporción es y de qué tipo son los cambios.",
                miles(diferencias),
                miles(longitud)
            ),
        );
    Ok((narracion.terminar(ultimo), diferencias))
}

fn pasos_distancia_hamming(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    Ok(trazar_hamming(entrada)?.0)
}

fn trazar_numero_de_patron(entrada: &Entrada) -> Result<(Vec<Paso>, u64), String> {
    let patron = entrada.adn("patron")?;
    replicacion::numero_de_patron(&patron)?;
    let k = patron.len();
    let tabla_de_valores = vec![
        vec!["A".to_string(), "0".to_string()],
        vec!["C".to_string(), "1".to_string()],
        vec!["G".to_string(), "2".to_string()],
        vec!["T".to_string(), "3".to_string()],
    ];
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        Pizarra::nueva()
            .variable("patron", letras(&patron))
            .variable("numero", 0)
            .secuencia(&patron, &[])
            .tabla("Valor de cada letra", &["letra", "valor"], tabla_de_valores.clone(), &[])
            .paso(
                "let mut numero = 0;",
                format!(
                    "Queremos escribir «{}» como un número. Cada letra vale un dígito en base 4 (A=0, C=1, G=2, T=3) y el patrón se lee de izquierda a derecha, como un número escrito en el sistema decimal. Empezamos con numero = 0.",
                    letras(&patron)
                ),
            ),
    );
    let mut numero: u64 = 0;
    for (lugar, &base) in patron.iter().enumerate() {
        let valor = replicacion::valor_de_base(base).unwrap_or(0);
        let antes = numero;
        numero = numero * 4 + valor;
        narracion.agregar(
            Pizarra::nueva()
                .variable("patron", letras(&patron))
                .variable("base", letra_de(base))
                .variable("valor", valor)
                .variable("numero", numero)
                .secuencia(&patron, &[(lugar, lugar + 1, Resalte::Actual)])
                .tabla(
                    "Valor de cada letra",
                    &["letra", "valor"],
                    tabla_de_valores.clone(),
                    &[(valor as usize, 1)],
                )
                .paso(
                    "numero = numero * 4 + valor;",
                    format!(
                        "Leemos «{}», que vale {valor}. Lo que llevábamos ({}) se corre un lugar a la izquierda, o sea se multiplica por 4, y se le suma el valor de la letra: 4 × {} + {valor} = {}.",
                        letra_de(base),
                        miles(antes as usize),
                        miles(antes as usize),
                        miles(numero as usize)
                    ),
                ),
        );
    }
    let ultimo = Pizarra::nueva()
        .variable("numero", miles(numero as usize))
        .paso(
            "Ok(numero)",
            format!(
                "Terminamos: «{}» es el número {}. Es su lugar en el orden alfabético entre los 4^{k} patrones de {k} letras, contando desde 0. En «Qué significa» ves cuántos hay en total.",
                letras(&patron),
                miles(numero as usize)
            ),
        );
    Ok((narracion.terminar(ultimo), numero))
}

fn pasos_numero_de_patron(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    Ok(trazar_numero_de_patron(entrada)?.0)
}

fn trazar_patron_de_numero(entrada: &Entrada) -> Result<(Vec<Paso>, Vec<u8>), String> {
    let numero = entrada.natural("numero")? as u64;
    let k = entrada.natural("k")?;
    if k == 0 {
        return Err("k debe ser al menos 1.".to_string());
    }
    if k > MAXIMO_K_NUMERICO {
        return Err(format!(
            "k puede valer como mucho {MAXIMO_K_NUMERICO}: un número de 64 bits no guarda patrones más largos."
        ));
    }
    replicacion::patron_de_numero(numero, k)?;
    let cantidad: u64 = 1 << (2 * k);
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        Pizarra::nueva()
            .variable("numero", miles(numero as usize))
            .variable("k", k)
            .variable("cantidad", miles(cantidad as usize))
            .paso(
                "let cantidad: u64 = 1 << (2 * k);",
                format!(
                    "Con {k} letras hay 4^{k} = {} patrones, numerados del 0 al {}. El {} está dentro de ese rango, así que existe un patrón con ese número.",
                    miles(cantidad as usize),
                    miles(cantidad as usize - 1),
                    miles(numero as usize)
                ),
            ),
    );
    let mut patron = vec![b'A'; k];
    narracion.agregar(
        Pizarra::nueva()
            .variable("numero", miles(numero as usize))
            .variable("k", k)
            .variable("patron", letras(&patron))
            .variable("resto", miles(numero as usize))
            .paso(
                "let mut patron = vec![b'A'; k];",
                format!(
                    "Preparamos un patrón de {k} letras, todas A por ahora. Vamos a escribir {} en base 4, de derecha a izquierda: cada vuelta decide una letra.",
                    miles(numero as usize)
                ),
            ),
    );
    let mut resto = numero;
    let mut historial: Vec<Vec<String>> = Vec::new();
    for (vuelta, posicion) in (0..k).rev().enumerate() {
        let cifra = (resto % 4) as usize;
        let antes = resto;
        patron[posicion] = replicacion::BASES[cifra];
        historial.push(vec![
            posicion.to_string(),
            miles(antes as usize),
            cifra.to_string(),
            letra_de(patron[posicion]),
        ]);
        let desde = historial.len().saturating_sub(FILAS_DEL_PASO);
        let filas: Vec<Vec<String>> = historial[desde..].to_vec();
        let ultima = filas.len() - 1;
        let detalle = vuelta < 6;
        if detalle || vuelta + 1 == k {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("posicion", posicion)
                    .variable("resto", miles(antes as usize))
                    .variable("resto % 4", cifra)
                    .variable("patron", letras(&patron))
                    .tabla(
                        "Cifras escritas",
                        &["posición", "resto", "resto % 4", "letra"],
                        filas,
                        &[(ultima, 2), (ultima, 3)],
                    )
                    .paso(
                        "patron[posicion] = BASES[(resto % 4) as usize];",
                        format!(
                            "Posición {posicion}: {} dividido entre 4 deja de resto {cifra}, y la cifra {cifra} es la letra «{}». La escribimos en el patrón, que ahora es «{}».",
                            miles(antes as usize),
                            letra_de(patron[posicion]),
                            letras(&patron)
                        ),
                    ),
            );
        }
        resto /= 4;
        if vuelta < 3 {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("posicion", posicion)
                    .variable("resto", miles(resto as usize))
                    .variable("patron", letras(&patron))
                    .paso(
                        "resto /= 4;",
                        format!(
                            "Dividimos entre 4 (división entera, sin decimales): queda {}. Con eso pasamos a la cifra de la izquierda.",
                            miles(resto as usize)
                        ),
                    ),
            );
        }
    }
    let ultimo = Pizarra::nueva()
        .variable("patron", letras(&patron))
        .paso(
            "Ok(patron)",
            format!(
                "Terminamos: el número {} escrito con {k} letras es «{}». Es el camino inverso de «De patrón a número».",
                miles(numero as usize),
                letras(&patron)
            ),
        );
    Ok((narracion.terminar(ultimo), patron))
}

fn pasos_patron_de_numero(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    Ok(trazar_patron_de_numero(entrada)?.0)
}

fn tabla_de_casillas(
    tocadas: &[(Vec<u8>, u64)],
    frecuencias: &[u32],
    actual: Option<u64>,
) -> (Vec<Vec<String>>, Vec<(usize, usize)>) {
    let mut filas = Vec::new();
    let mut resaltar = Vec::new();
    for (fila, (patron, numero)) in tocadas.iter().enumerate() {
        if Some(*numero) == actual {
            resaltar.push((fila, 2));
        }
        filas.push(vec![
            letras(patron),
            miles(*numero as usize),
            miles(frecuencias[*numero as usize] as usize),
        ]);
    }
    (filas, resaltar)
}

fn llenar_arreglo(narracion: &mut Narracion, texto: &[u8], k: usize) -> Result<Vec<u32>, String> {
    let mut frecuencias = replicacion::arreglo_de_frecuencias(&[], k)?;
    let ventanas = if k <= texto.len() {
        texto.len() - k + 1
    } else {
        0
    };
    let controles = puntos_de_control(ventanas);
    narracion.agregar(
        Pizarra::nueva()
            .variable("texto", vista(texto))
            .variable("k", k)
            .variable(
                "frecuencias",
                format!("{} casillas, todas en 0", miles(frecuencias.len())),
            )
            .secuencia(texto, &[])
            .paso(
                "let mut frecuencias = vec![0u32; 1 << (2 * k)];",
                format!(
                    "Con k = {k} hay 4^{k} = {} k-meros posibles. Preparamos una casilla para cada uno, todas en 0 y en orden alfabético: la casilla 0 es {} y la última es {}. Para saber a qué casilla va cada k-mero lo convertiremos en un número.",
                    miles(frecuencias.len()),
                    "A".repeat(k),
                    "T".repeat(k)
                ),
            ),
    );
    let mut tocadas: Vec<(Vec<u8>, u64)> = Vec::new();
    let mut distintas = 0usize;
    for (indice, kmero) in texto.windows(k).enumerate() {
        let actual = [(indice, indice + k, Resalte::Actual)];
        match replicacion::numero_de_patron(kmero) {
            Ok(numero) => {
                if indice == 0 {
                    let mut acumulado = 0u64;
                    for &base in kmero {
                        let valor = replicacion::valor_de_base(base).unwrap_or(0);
                        let antes = acumulado;
                        acumulado = acumulado * 4 + valor;
                        narracion.agregar(
                            Pizarra::nueva()
                                .variable("kmero", letras(kmero))
                                .variable("base", letra_de(base))
                                .variable("valor", valor)
                                .variable("numero", acumulado)
                                .secuencia(texto, &actual)
                                .paso(
                                    "numero = numero * 4 + valor;",
                                    format!(
                                        "Convertimos el primer k-mero, «{}», en un número. Letra «{}» = {valor}: 4 × {} + {valor} = {}.",
                                        letras(kmero),
                                        letra_de(base),
                                        miles(antes as usize),
                                        miles(acumulado as usize)
                                    ),
                                ),
                        );
                    }
                }
                frecuencias[numero as usize] += 1;
                let veces = frecuencias[numero as usize];
                if veces == 1 {
                    distintas += 1;
                }
                if tocadas.len() < 8 && !tocadas.iter().any(|(_, otro)| *otro == numero) {
                    tocadas.push((kmero.to_vec(), numero));
                }
                if indice < DETALLE {
                    let (filas, resaltar) = tabla_de_casillas(&tocadas, &frecuencias, Some(numero));
                    let nota = if veces == 1 {
                        format!(
                            "«{}» es el número {}: sumamos 1 a esa casilla, que pasa de 0 a 1. Es la primera vez que vemos este k-mero.",
                            letras(kmero),
                            miles(numero as usize)
                        )
                    } else {
                        format!(
                            "«{}» es el número {}: sumamos 1 a esa casilla, que pasa de {} a {veces}. Ya lo habíamos visto antes.",
                            letras(kmero),
                            miles(numero as usize),
                            veces - 1
                        )
                    };
                    narracion.agregar(
                        Pizarra::nueva()
                            .variable("kmero", letras(kmero))
                            .variable("numero", numero)
                            .variable("frecuencias[numero]", veces)
                            .secuencia(texto, &actual)
                            .tabla(
                                "Casillas distintas de 0 (en el orden en que las vimos)",
                                &["k-mero", "casilla", "veces"],
                                filas,
                                &resaltar,
                            )
                            .paso("frecuencias[numero as usize] += 1;", nota),
                    );
                }
            }
            Err(_) => {
                if indice < DETALLE {
                    narracion.agregar(
                        Pizarra::nueva()
                            .variable("kmero", letras(kmero))
                            .secuencia(texto, &actual)
                            .paso(
                                "if let Ok(numero) = numero_de_patron(kmero) {",
                                format!(
                                    "«{}» tiene una letra que no es A, C, G ni T: no se puede convertir en número y se salta.",
                                    letras(kmero)
                                ),
                            ),
                    );
                }
            }
        }
        if indice >= DETALLE && controles.contains(&indice) {
            let (filas, resaltar) = tabla_de_casillas(&tocadas, &frecuencias, None);
            narracion.agregar(
                Pizarra::nueva()
                    .variable("kmero", letras(kmero))
                    .secuencia(texto, &actual)
                    .tabla(
                        "Casillas distintas de 0 (en el orden en que las vimos)",
                        &["k-mero", "casilla", "veces"],
                        filas,
                        &resaltar,
                    )
                    .paso(
                        "if let Ok(numero) = numero_de_patron(kmero) {",
                        format!(
                            "Vamos por el k-mero {} de {}: ya hay {} casillas distintas de 0. Cada vuelta hace lo mismo: convertir, y sumar 1.",
                            miles(indice + 1),
                            miles(ventanas),
                            miles(distintas)
                        ),
                    ),
            );
        }
    }
    Ok(frecuencias)
}

fn trazar_arreglo(entrada: &Entrada) -> Result<(Vec<Paso>, Vec<u32>), String> {
    let texto = entrada.secuencia("texto")?;
    let k = leer_k(entrada, "k", texto.len(), MAXIMO_K_ENUMERABLE)?;
    revisar_animable(texto.len())?;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let frecuencias = llenar_arreglo(&mut narracion, &texto, k)?;
    let presentes = frecuencias.iter().filter(|&&veces| veces > 0).count();
    let mut elementos = Vec::new();
    for (numero, veces) in mejores_del_arreglo(&frecuencias, 8) {
        elementos.push(format!(
            "{} (casilla {}): {} veces",
            letras(&replicacion::patron_de_numero(numero as u64, k)?),
            miles(numero),
            miles(veces as usize)
        ));
    }
    let resaltar: Vec<usize> = if elementos.is_empty() {
        Vec::new()
    } else {
        vec![0]
    };
    let ultimo = Pizarra::nueva()
        .variable("k", k)
        .variable(
            "frecuencias",
            format!(
                "{} casillas, {} distintas de 0",
                miles(frecuencias.len()),
                miles(presentes)
            ),
        )
        .lista("Las casillas más altas", elementos, &resaltar)
        .paso(
            "Ok(frecuencias)",
            format!(
                "Terminamos de leer el texto: de las {} casillas, {} tienen al menos una aparición. La lista de abajo muestra las más altas; en «Qué significa» ves el arreglo completo.",
                miles(frecuencias.len()),
                miles(presentes)
            ),
        );
    Ok((narracion.terminar(ultimo), frecuencias))
}

fn pasos_arreglo_de_frecuencias(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    Ok(trazar_arreglo(entrada)?.0)
}

fn trazar_palabras_frecuentes(entrada: &Entrada) -> Result<(Vec<Paso>, Vec<Vec<u8>>), String> {
    let texto = entrada.adn("texto")?;
    let k = leer_k(entrada, "k", texto.len(), MAXIMO_K_TABLA)?;
    revisar_animable(texto.len())?;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let casillas = 4f64.powi(k as i32);
    let tope = texto.len().saturating_mul(4);
    let casillas_texto = if casillas < 1e15 {
        miles(casillas as usize)
    } else {
        cientifico(casillas)
    };
    narracion.agregar(
        Pizarra::nueva()
            .variable("texto", vista(&texto))
            .variable("k", k)
            .secuencia(&texto, &[])
            .paso(
                "pub fn palabras_frecuentes(texto: &[u8], k: usize) -> Vec<Vec<u8>> {",
                format!(
                    "Queremos los k-meros de {k} letras que más se repiten en un texto de {} letras. Hay dos maneras de contarlos y el código elige una según el tamaño del problema.",
                    miles(texto.len())
                ),
            ),
    );
    let con_arreglo = replicacion::conviene_un_arreglo(texto.len(), k);
    narracion.agregar(
        Pizarra::nueva()
            .variable("longitud", miles(texto.len()))
            .variable("k", k)
            .variable("4^k", &casillas_texto)
            .variable("longitud * 4", miles(tope))
            .paso(
                "(1 << (2 * k)) <= longitud.saturating_mul(4)",
                if con_arreglo {
                    format!(
                        "Un arreglo necesita 4^{k} = {casillas_texto} casillas, y 4 veces el largo del texto es {}. Cabe de sobra, así que usamos un arreglo con una casilla por cada k-mero posible.",
                        miles(tope)
                    )
                } else {
                    format!(
                        "Un arreglo necesitaría 4^{k} = {casillas_texto} casillas (o k pasa de {MAXIMO_K_ENUMERABLE}), y 4 veces el largo del texto es solo {}. Sería desperdiciar memoria, así que usamos una tabla que guarda únicamente los k-meros que de verdad aparecen.",
                        miles(tope)
                    )
                },
            ),
    );
    if con_arreglo {
        let frecuencias = llenar_arreglo(&mut narracion, &texto, k)?;
        let maximo = frecuencias.iter().copied().max().unwrap_or(0);
        let mejores = mejores_del_arreglo(&frecuencias, 8);
        let mut elementos = Vec::new();
        for (numero, veces) in &mejores {
            elementos.push(format!(
                "{}: {} veces",
                letras(&replicacion::patron_de_numero(*numero as u64, k)?),
                miles(*veces as usize)
            ));
        }
        let resaltar: Vec<usize> = mejores
            .iter()
            .enumerate()
            .filter(|(_, (_, veces))| *veces == maximo)
            .map(|(lugar, _)| lugar)
            .collect();
        narracion.agregar(
            Pizarra::nueva()
                .variable("k", k)
                .variable("maximo", miles(maximo as usize))
                .lista("Las casillas más altas", elementos, &resaltar)
                .paso(
                    "let maximo = cuentas.iter().copied().max().unwrap_or(0);",
                    format!(
                        "Ya está el arreglo lleno. Buscamos la casilla más alta: el máximo es {}. Todas las casillas que lleguen a ese valor son las ganadoras.",
                        miles(maximo as usize)
                    ),
                ),
        );
        let frecuentes = replicacion::los_mas_frecuentes(&frecuencias, k);
        let ganadores: Vec<String> = frecuentes.iter().take(20).map(|p| letras(p)).collect();
        narracion.agregar(
            Pizarra::nueva()
                .variable("maximo", miles(maximo as usize))
                .variable("patrones", ganadores.len())
                .lista(
                    "Ganadores",
                    ganadores.clone(),
                    &(0..ganadores.len()).collect::<Vec<_>>(),
                )
                .paso(
                    "patrones.push(patron);",
                    format!(
                        "Recorremos las casillas en orden y, por cada una que vale {}, convertimos su número otra vez en letras y la guardamos. Hay {}.",
                        miles(maximo as usize),
                        plural(frecuentes.len(), "ganador", "ganadores")
                    ),
                ),
        );
        let ultimo = Pizarra::nueva()
            .variable("maximo", miles(maximo as usize))
            .lista(
                "Los más frecuentes",
                ganadores.clone(),
                &(0..ganadores.len()).collect::<Vec<_>>(),
            )
            .paso(
                "return los_mas_frecuentes(&frecuencias, k);",
                format!(
                    "Terminamos: {} de {k} letras se repite{} {} veces, el máximo. Salen en orden alfabético porque es el orden de las casillas.",
                    if frecuentes.len() == 1 {
                        "un k-mero".to_string()
                    } else {
                        format!("{} k-meros", miles(frecuentes.len()))
                    },
                    if frecuentes.len() == 1 { "" } else { "n" },
                    miles(maximo as usize)
                ),
            );
        return Ok((narracion.terminar(ultimo), frecuentes));
    }
    let ventanas = texto.len() - k + 1;
    let controles = puntos_de_control(ventanas);
    let mut tabla: HashMap<&[u8], usize> = HashMap::new();
    let mut orden: Vec<&[u8]> = Vec::new();
    for (indice, kmero) in texto.windows(k).enumerate() {
        let actual = [(indice, indice + k, Resalte::Actual)];
        let valido = kmero
            .iter()
            .all(|&base| replicacion::valor_de_base(base).is_some());
        if valido {
            let entrada_de_la_tabla = tabla.entry(kmero).or_insert(0);
            *entrada_de_la_tabla += 1;
            let veces = *entrada_de_la_tabla;
            if veces == 1 && orden.len() < 8 {
                orden.push(kmero);
            }
            if indice < DETALLE {
                let filas: Vec<Vec<String>> = orden
                    .iter()
                    .map(|patron| vec![letras(patron), miles(tabla[patron])])
                    .collect();
                let resaltar: Vec<(usize, usize)> = orden
                    .iter()
                    .enumerate()
                    .filter(|(_, patron)| **patron == kmero)
                    .map(|(fila, _)| (fila, 1))
                    .collect();
                narracion.agregar(
                    Pizarra::nueva()
                        .variable("kmero", letras(kmero))
                        .variable("tabla.len()", miles(tabla.len()))
                        .secuencia(&texto, &actual)
                        .tabla(
                            "Tabla (k-mero y veces)",
                            &["k-mero", "veces"],
                            filas,
                            &resaltar,
                        )
                        .paso(
                            "*tabla.entry(kmero).or_insert(0) += 1;",
                            if veces == 1 {
                                format!(
                                    "«{}» todavía no estaba en la tabla: se agrega con 0 y se le suma 1.",
                                    letras(kmero)
                                )
                            } else {
                                format!(
                                    "«{}» ya estaba en la tabla: su cuenta sube a {veces}.",
                                    letras(kmero)
                                )
                            },
                        ),
                );
            }
        } else if indice < DETALLE {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("kmero", letras(kmero))
                    .secuencia(&texto, &actual)
                    .paso(
                        "if kmero.iter().all(|&base| valor_de_base(base).is_some()) {",
                        format!(
                            "«{}» tiene una letra que no es A, C, G ni T: no se cuenta.",
                            letras(kmero)
                        ),
                    ),
            );
        }
        if indice >= DETALLE && controles.contains(&indice) {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("kmero", letras(kmero))
                    .variable("tabla.len()", miles(tabla.len()))
                    .secuencia(&texto, &actual)
                    .paso(
                        "if kmero.iter().all(|&base| valor_de_base(base).is_some()) {",
                        format!(
                            "Vamos por el k-mero {} de {}: la tabla ya tiene {} k-meros distintos.",
                            miles(indice + 1),
                            miles(ventanas),
                            miles(tabla.len())
                        ),
                    ),
            );
        }
    }
    let maximo = tabla.values().copied().max().unwrap_or(0);
    narracion.agregar(
        Pizarra::nueva()
            .variable("tabla.len()", miles(tabla.len()))
            .variable("maximo", miles(maximo))
            .paso(
                "let maximo = tabla.values().copied().max().unwrap_or(0);",
                format!(
                    "Terminamos de contar: la tabla tiene {} k-meros distintos y el que más se repite aparece {} veces.",
                    miles(tabla.len()),
                    miles(maximo)
                ),
            ),
    );
    let mut frecuentes: Vec<Vec<u8>> = Vec::new();
    for (kmero, veces) in &tabla {
        if *veces == maximo {
            frecuentes.push(kmero.to_vec());
        }
    }
    narracion.agregar(
        Pizarra::nueva()
            .variable("maximo", miles(maximo))
            .variable("frecuentes", miles(frecuentes.len()))
            .paso(
                "if veces == maximo {",
                format!(
                    "Revisamos la tabla y guardamos los k-meros cuya cuenta es {}: hay {}.",
                    miles(maximo),
                    plural(frecuentes.len(), "ganador", "ganadores")
                ),
            ),
    );
    frecuentes.sort();
    let ganadores: Vec<String> = frecuentes.iter().take(20).map(|p| letras(p)).collect();
    let ultimo = Pizarra::nueva()
        .variable("maximo", miles(maximo))
        .lista(
            "Los más frecuentes",
            ganadores.clone(),
            &(0..ganadores.len()).collect::<Vec<_>>(),
        )
        .paso(
            "frecuentes.sort();",
            "Como la tabla no guarda ningún orden, ordenamos los ganadores alfabéticamente para que el resultado siempre salga igual. Ese es el resultado: mira «Qué significa».",
        );
    Ok((narracion.terminar(ultimo), frecuentes))
}

fn pasos_palabras_frecuentes(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    Ok(trazar_palabras_frecuentes(entrada)?.0)
}

fn trazar_grumos(entrada: &Entrada) -> Result<(Vec<Paso>, Vec<replicacion::Grumo>), String> {
    let genoma = entrada.secuencia("genoma")?;
    let k = leer_k(entrada, "k", genoma.len(), MAXIMO_K_ENUMERABLE)?;
    let largo = entrada.natural("ventana")?;
    if largo < k || largo > genoma.len() {
        return Err(format!(
            "La ventana L debe medir entre k ({k}) y el largo del genoma ({} pb).",
            miles(genoma.len())
        ));
    }
    let caben = largo - k + 1;
    let veces = entrada.natural("veces")?;
    if veces == 0 || veces > caben {
        return Err(format!(
            "t debe estar entre 1 y {}: en una ventana de {} pb caben {} k-meros de {k} letras.",
            miles(caben),
            miles(largo),
            miles(caben)
        ));
    }
    revisar_animable(genoma.len())?;
    let texto = &genoma[..];
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", k)
            .variable("largo", miles(largo))
            .variable("veces", veces)
            .secuencia(texto, &[(0, largo, Resalte::Comparado)])
            .paso(
                "if veces == 0 || largo < k || texto.len() < largo {",
                format!(
                    "Buscamos k-meros de {k} letras que aparezcan al menos {veces} veces dentro de alguna ventana de {} letras del genoma ({} letras en total). Los datos tienen sentido, así que seguimos: primero contamos todo en la primera ventana y luego la deslizamos una letra a la vez, sin volver a contar desde cero.",
                    miles(largo),
                    miles(texto.len())
                ),
            ),
    );
    let mut en_la_ventana = replicacion::arreglo_de_frecuencias(&texto[..largo], k)?;
    let distintos = en_la_ventana.iter().filter(|&&veces| veces > 0).count();
    let mejor = mejores_del_arreglo(&en_la_ventana, 1);
    let masivo = match mejor.first() {
        Some(&(numero, repeticiones)) => format!(
            "El más repetido es {} con {} apariciones.",
            letras(&replicacion::patron_de_numero(numero as u64, k)?),
            miles(repeticiones as usize)
        ),
        None => String::new(),
    };
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", k)
            .variable("largo", miles(largo))
            .variable(
                "en_la_ventana",
                format!("{} casillas distintas de 0", miles(distintos)),
            )
            .secuencia(texto, &[(0, largo, Resalte::Comparado)])
            .paso(
                "let mut en_la_ventana = arreglo_de_frecuencias(&texto[..largo], k)?;",
                format!(
                    "Contamos los k-meros de las primeras {} letras en un arreglo con una casilla por k-mero: hay {} k-meros distintos. {masivo}",
                    miles(largo),
                    miles(distintos)
                ),
            ),
    );
    let mut primera_ventana: HashMap<u64, usize> = HashMap::new();
    let mut orden: Vec<u64> = Vec::new();
    for kmero in texto[..largo].windows(k) {
        if let Ok(numero) = replicacion::numero_de_patron(kmero)
            && en_la_ventana[numero as usize] as usize >= veces
            && !primera_ventana.contains_key(&numero)
        {
            primera_ventana.insert(numero, 0);
            orden.push(numero);
        }
    }
    let nombres: Vec<String> = orden
        .iter()
        .take(12)
        .map(|&numero| letras(&replicacion::patron_de_numero(numero, k).unwrap_or_default()))
        .collect();
    let resaltar_nombres: Vec<usize> = (0..nombres.len()).collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("veces", veces)
            .variable(
                "primera_ventana",
                format!("{} k-meros", miles(primera_ventana.len())),
            )
            .secuencia(texto, &[(0, largo, Resalte::Comparado)])
            .lista("Grumos desde la primera ventana", nombres, &resaltar_nombres)
            .paso(
                "primera_ventana.entry(numero).or_insert(0);",
                if orden.is_empty() {
                    format!(
                        "Ningún k-mero llega a {veces} apariciones en la primera ventana. Si aparece alguno más adelante, lo anotaremos con la posición de la ventana donde lo logró."
                    )
                } else {
                    format!(
                        "Ya en la primera ventana hay {} que llegan a {veces} apariciones: son grumos desde la posición 0. Los anotamos con 0, que es la ventana donde los vimos por primera vez.",
                        plural(orden.len(), "k-mero", "k-meros")
                    )
                },
            ),
    );
    let deslizamientos = texto.len() - largo;
    let controles = puntos_de_control(deslizamientos);
    let mut eventos = 0usize;
    let mut avisado = false;
    for inicio in 1..=deslizamientos {
        let saliente = &texto[inicio - 1..inicio - 1 + k];
        let ultimo = inicio + largo - k;
        let entrante = &texto[ultimo..ultimo + k];
        let marcas = [
            (inicio, inicio + largo, Resalte::Comparado),
            (inicio - 1, inicio - 1 + k, Resalte::Actual),
            (ultimo, ultimo + k, Resalte::Actual),
        ];
        let numero_saliente = replicacion::numero_de_patron(saliente);
        if let Ok(numero) = numero_saliente {
            en_la_ventana[numero as usize] -= 1;
        }
        let numero_entrante = replicacion::numero_de_patron(entrante);
        let mut nuevo = None;
        if let Ok(numero) = numero_entrante {
            en_la_ventana[numero as usize] += 1;
            if en_la_ventana[numero as usize] as usize >= veces
                && !primera_ventana.contains_key(&numero)
            {
                primera_ventana.insert(numero, inicio);
                orden.push(numero);
                nuevo = Some(numero);
            }
        }
        let variables = |pizarra: Pizarra| {
            pizarra
                .variable("largo", miles(largo))
                .variable("veces", veces)
                .variable("inicio", miles(inicio))
                .variable("saliente", letras(saliente))
                .variable("entrante", letras(entrante))
                .variable(
                    "primera_ventana",
                    format!("{} k-meros", miles(primera_ventana.len())),
                )
                .secuencia(texto, &marcas)
        };
        if inicio <= 3 {
            if let Ok(numero) = numero_saliente {
                narracion.agregar(variables(Pizarra::nueva()).variable("numero", numero).paso(
                    "en_la_ventana[numero as usize] -= 1;",
                    format!(
                        "La ventana se corre una letra: sale «{}», el k-mero que empezaba en la posición {}. Su casilla baja a {}.",
                        letras(saliente),
                        miles(inicio - 1),
                        en_la_ventana[numero as usize]
                    ),
                ));
            }
            if let Ok(numero) = numero_entrante {
                narracion.agregar(variables(Pizarra::nueva()).variable("numero", numero).paso(
                    "en_la_ventana[numero as usize] += 1;",
                    format!(
                        "Y entra «{}», el nuevo k-mero del final. Su casilla sube a {}{}.",
                        letras(entrante),
                        en_la_ventana[numero as usize],
                        if nuevo.is_some() {
                            format!(", que ya llega a las {veces} apariciones que buscamos")
                        } else {
                            String::new()
                        }
                    ),
                ));
            }
        }
        if let Some(numero) = nuevo {
            if eventos < MAXIMO_EVENTOS {
                eventos += 1;
                narracion.agregar(variables(Pizarra::nueva()).variable("numero", numero).paso(
                    "primera_ventana.entry(numero).or_insert(inicio);",
                    format!(
                        "¡Nuevo grumo! «{}» acaba de llegar a {veces} apariciones dentro de la ventana que empieza en la posición {}. Lo anotamos con esa posición, y ya no se vuelve a anotar aunque siga apareciendo.",
                        letras(entrante),
                        miles(inicio)
                    ),
                ));
            } else if !avisado {
                avisado = true;
                narracion.agregar(variables(Pizarra::nueva()).variable("numero", numero).paso(
                    "primera_ventana.entry(numero).or_insert(inicio);",
                    "Seguimos hallando grumos nuevos igual que este. Para no cansarte, de aquí en adelante los anoto en silencio y te aviso del total en unos cuantos puntos del recorrido.",
                ));
            }
        }
        if inicio > 3 && controles.contains(&inicio) {
            narracion.agregar(variables(Pizarra::nueva()).paso(
                "for inicio in 1..=texto.len() - largo {",
                format!(
                    "La ventana va en la posición {} de {}: hasta aquí hay {} que forman grumos.",
                    miles(inicio),
                    miles(deslizamientos),
                    plural(primera_ventana.len(), "k-mero", "k-meros")
                ),
            ));
        }
    }
    let mut grumos = Vec::new();
    for (numero, ventana) in primera_ventana {
        grumos.push(replicacion::Grumo {
            patron: replicacion::patron_de_numero(numero, k)?,
            ventana,
        });
    }
    grumos.sort_by(|uno, otro| uno.patron.cmp(&otro.patron));
    let elementos: Vec<String> = grumos
        .iter()
        .take(20)
        .map(|grumo| {
            format!(
                "{} (ventana {})",
                letras(&grumo.patron),
                miles(grumo.ventana)
            )
        })
        .collect();
    let resaltar_todos: Vec<usize> = (0..elementos.len()).collect();
    let ultimo = Pizarra::nueva()
        .variable("grumos", miles(grumos.len()))
        .lista("k-meros que forman grumos", elementos, &resaltar_todos)
        .paso(
            "grumos.sort_by(|uno, otro| uno.patron.cmp(&otro.patron));",
            format!(
                "Terminamos de deslizar la ventana: {} forman grumos. Los ordenamos alfabéticamente y, junto a cada uno, la primera ventana donde lo vimos. En «Qué significa» verás si es mucho para un genoma de este tamaño.",
                plural(grumos.len(), "k-mero", "k-meros")
            ),
        );
    Ok((narracion.terminar(ultimo), grumos))
}

fn pasos_grumos(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    Ok(trazar_grumos(entrada)?.0)
}

fn tabla_de_mejores(
    cuentas: &[u32],
    exactas: Option<&[u32]>,
    k: usize,
    cantidad: usize,
) -> Result<Vec<Vec<String>>, String> {
    let mut filas = Vec::new();
    for (numero, veces) in mejores_del_arreglo(cuentas, cantidad) {
        let mut fila = vec![letras(&replicacion::patron_de_numero(numero as u64, k)?)];
        if let Some(exactas) = exactas {
            fila.push(miles(exactas[numero] as usize));
        }
        fila.push(miles(veces as usize));
        filas.push(fila);
    }
    Ok(filas)
}

fn trazar_frecuentes(
    entrada: &Entrada,
    dos_hebras: bool,
) -> Result<(Vec<Paso>, Vec<Vec<u8>>), String> {
    let texto = entrada.secuencia("texto")?;
    let k = leer_k(entrada, "k", texto.len(), MAXIMO_K_ENUMERABLE)?;
    let d = leer_d(entrada, "d", k, MAXIMO_D)?;
    revisar_animable(texto.len())?;
    let cuentas = if dos_hebras {
        replicacion::cuentas_en_ambas_hebras(&texto, k, d)?
    } else {
        replicacion::cuentas_con_errores(&texto, k, d)?
    };
    let ventanas = texto.len() - k + 1;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let linea_de_inicio = if dos_hebras {
        "let cuentas = cuentas_en_ambas_hebras(texto, k, d)?;"
    } else {
        "let cuentas = cuentas_con_errores(texto, k, d)?;"
    };
    let hebras = if dos_hebras {
        " Además cada k-mero suma las apariciones de su complemento inverso, porque la señal puede estar en cualquiera de las dos hebras."
    } else {
        ""
    };
    narracion.agregar(
        Pizarra::nueva()
            .variable("texto", vista(&texto))
            .variable("k", k)
            .variable("d", d)
            .secuencia(&texto, &[])
            .paso(
                linea_de_inicio,
                format!(
                    "Queremos los k-meros de {k} letras que más aparecen en un texto de {} letras si cada aparición puede tener hasta {} distinta{}.{hebras} Lo haremos en tres etapas: 1) contar las apariciones exactas, 2) repartir cada cuenta entre los k-meros parecidos, 3) quedarnos con los de mayor cuenta.",
                    miles(texto.len()),
                    plural(d, "letra", "letras"),
                    if d == 1 { "" } else { "s" }
                ),
            ),
    );
    let mut exactas = replicacion::arreglo_de_frecuencias(&texto, k)?;
    let distintos = exactas.iter().filter(|&&veces| veces > 0).count();
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", k)
            .variable("exactas", format!("{} k-meros distintos", miles(distintos)))
            .tabla(
                "Los k-meros más repetidos (exactos)",
                &["k-mero", "veces"],
                tabla_de_mejores(&exactas, None, k, 6)?,
                &[(0, 1)],
            )
            .paso(
                "exactas = arreglo_de_frecuencias(texto, k)?;",
                format!(
                    "Etapa 1: contamos cada k-mero tal cual, sin errores. Las {} ventanas del texto dan {} k-meros distintos. Así se ve la parte de arriba de la cuenta.",
                    miles(ventanas),
                    miles(distintos)
                ),
            ),
    );
    if dos_hebras {
        let inverso = secuencia::complemento_inverso(&texto);
        let otra_hebra = replicacion::arreglo_de_frecuencias(&inverso, k)?;
        let otras = otra_hebra.iter().filter(|&&veces| veces > 0).count();
        narracion.agregar(
            Pizarra::nueva()
                .variable("k", k)
                .variable(
                    "otra_hebra",
                    format!("{} k-meros distintos", miles(otras)),
                )
                .tabla(
                    "Los k-meros más repetidos en la otra hebra",
                    &["k-mero", "veces"],
                    tabla_de_mejores(&otra_hebra, None, k, 6)?,
                    &[(0, 1)],
                )
                .paso(
                    "let otra_hebra = arreglo_de_frecuencias(&complemento_inverso(texto), k)?;",
                    "Hacemos la misma cuenta sobre la otra hebra: el complemento inverso del texto. Un k-mero que aparece allí es una copia de su complemento inverso aquí.",
                ),
        );
        for (cuenta, otra) in exactas.iter_mut().zip(&otra_hebra) {
            *cuenta += otra;
        }
        narracion.agregar(
            Pizarra::nueva()
                .variable("k", k)
                .tabla(
                    "Suma de las dos hebras",
                    &["k-mero", "veces"],
                    tabla_de_mejores(&exactas, None, k, 6)?,
                    &[(0, 1)],
                )
                .paso(
                    "*cuenta += otra;",
                    "Sumamos casilla por casilla las dos cuentas. Ahora cada k-mero y su complemento inverso valen lo mismo, porque son la misma señal vista desde cada hebra.",
                ),
        );
    }
    let tamano = replicacion::tamano_del_vecindario(k, d);
    let visitas = distintos as f64 * tamano;
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", k)
            .variable("d", d)
            .variable("distintos", miles(distintos))
            .variable("visitas", miles(visitas as usize))
            .paso(
                "let visitas = distintos as f64 * tamano_del_vecindario(k, d);",
                format!(
                    "Etapa 2: cada k-mero con alguna aparición reparte su cuenta entre todos sus vecinos, o sea los k-meros que difieren de él en {} como mucho. Cada k-mero de {k} letras tiene {} vecinos (contándose a sí mismo), y hay {} k-meros que repartir: serán {} visitas. El código revisa que no sean demasiadas.",
                    plural(d, "letra", "letras"),
                    miles(tamano as usize),
                    miles(distintos),
                    miles(visitas as usize)
                ),
            ),
    );
    let mut aproximadas: HashMap<Vec<u8>, u32> = HashMap::new();
    let mut mostrados = 0usize;
    for (numero, &veces) in exactas.iter().enumerate() {
        if veces == 0 {
            continue;
        }
        if mostrados == 3 {
            break;
        }
        mostrados += 1;
        let patron = replicacion::patron_de_numero(numero as u64, k)?;
        let vecinos = replicacion::vecindario(&patron, d)?;
        for vecino in &vecinos {
            *aproximadas.entry(vecino.clone()).or_insert(0) += veces;
        }
        let filas: Vec<Vec<String>> = vecinos
            .iter()
            .take(6)
            .map(|vecino| {
                vec![
                    letras(vecino),
                    miles(diferencias_entre(vecino, &patron)),
                    miles(aproximadas[vecino] as usize),
                ]
            })
            .collect();
        narracion.agregar(
            Pizarra::nueva()
                .variable("numero", numero)
                .variable("veces", veces)
                .variable("k-mero", letras(&patron))
                .tabla(
                    "Vecinos del k-mero y su cuenta hasta ahora",
                    &["vecino", "diferencias", "cuenta"],
                    filas,
                    &[(0, 2)],
                )
                .paso(
                    "aproximadas[vecino as usize] += veces;",
                    format!(
                        "Tomamos el k-mero «{}», que apareció {} {}. Sus {} vecinos reciben esa cuenta: a cada uno se le suma {veces}. Por eso un k-mero que casi nunca aparece exacto puede acabar con una cuenta alta, si tiene muchos vecinos frecuentes.",
                        letras(&patron),
                        veces,
                        if veces == 1 { "vez" } else { "veces" },
                        miles(vecinos.len())
                    ),
                ),
        );
    }
    narracion.agregar(
        Pizarra::nueva()
            .variable("distintos", miles(distintos))
            .variable("visitas", miles(visitas as usize))
            .tabla(
                "Los k-meros con más cuenta (con errores)",
                if dos_hebras {
                    &["k-mero", "exactas, dos hebras", "con errores"]
                } else {
                    &["k-mero", "exactas", "con errores"]
                },
                tabla_de_mejores(&cuentas, Some(&exactas[..]), k, 8)?,
                &[(0, 2)],
            )
            .paso(
                "repartir_entre_vecinos(&exactas, k, d)",
                format!(
                    "Repetimos lo mismo con los {} k-meros distintos. Lo que viste con los primeros tres pasa con todos: se hicieron {} visitas. Al final cada casilla tiene su cuenta con errores.",
                    miles(distintos),
                    miles(visitas as usize)
                ),
            ),
    );
    let maximo = cuentas.iter().copied().max().unwrap_or(0);
    narracion.agregar(
        Pizarra::nueva()
            .variable("maximo", miles(maximo as usize))
            .paso(
                "let maximo = cuentas.iter().copied().max().unwrap_or(0);",
                format!(
                    "Etapa 3: buscamos la cuenta más alta de todas las casillas. Es {}.",
                    miles(maximo as usize)
                ),
            ),
    );
    let frecuentes = replicacion::los_mas_frecuentes(&cuentas, k);
    let ganadores: Vec<String> = frecuentes.iter().take(20).map(|p| letras(p)).collect();
    let resaltar_todos: Vec<usize> = (0..ganadores.len()).collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("maximo", miles(maximo as usize))
            .variable("patrones", miles(frecuentes.len()))
            .lista("Ganadores", ganadores.clone(), &resaltar_todos)
            .paso(
                "patrones.push(patron);",
                format!(
                    "Recorremos las casillas y guardamos todos los k-meros que valen {}: hay {}.",
                    miles(maximo as usize),
                    plural(frecuentes.len(), "ganador", "ganadores")
                ),
            ),
    );
    let cierre = if dos_hebras {
        "Terminamos. Fíjate si los ganadores vienen en parejas, un k-mero y su complemento inverso: es lo esperado. En «Qué significa» ves qué tan raro es ese número de apariciones."
    } else {
        "Terminamos. Estos son los k-meros más frecuentes con errores. En «Qué significa» ves qué tan raro es ese número de apariciones."
    };
    let ultimo = Pizarra::nueva()
        .variable("maximo", miles(maximo as usize))
        .lista("Los más frecuentes", ganadores, &resaltar_todos)
        .paso("Ok(los_mas_frecuentes(&cuentas, k))", cierre);
    Ok((narracion.terminar(ultimo), frecuentes))
}

fn pasos_frecuentes_con_errores(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    Ok(trazar_frecuentes(entrada, false)?.0)
}

fn pasos_frecuentes_con_complementos(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    Ok(trazar_frecuentes(entrada, true)?.0)
}

struct RecorridoDeVecinos<'a> {
    narracion: Narracion,
    patron: &'a [u8],
    k: usize,
    d: usize,
    tamano: usize,
    numeros: Vec<u64>,
    detalles: usize,
}

impl RecorridoDeVecinos<'_> {
    fn nombre(&self, numero: u64) -> String {
        letras(&replicacion::patron_de_numero(numero, self.k).unwrap_or_default())
    }

    fn lista_reciente(&self) -> (Vec<String>, Vec<usize>) {
        let desde = self.numeros.len().saturating_sub(12);
        let elementos: Vec<String> = self.numeros[desde..]
            .iter()
            .map(|&numero| self.nombre(numero))
            .collect();
        let ultimo = elementos.len().saturating_sub(1);
        (elementos, vec![ultimo])
    }

    fn visitar(&mut self, numero: u64) {
        self.numeros.push(numero);
        if self.detalles == 0 {
            return;
        }
        self.detalles -= 1;
        let nombre = self.nombre(numero);
        let diferencias = diferencias_entre(nombre.as_bytes(), self.patron);
        let (elementos, resaltar) = self.lista_reciente();
        let nota = if self.numeros.len() == 1 {
            format!(
                "Primero visitamos el propio patrón, «{nombre}»: es el vecino con 0 cambios. Cada visita guarda un número en la lista."
            )
        } else {
            format!(
                "Visitamos «{nombre}», que difiere de «{}» en {}. Ya llevamos {} de {}.",
                letras(self.patron),
                plural(diferencias, "letra", "letras"),
                plural(self.numeros.len(), "vecino", "vecinos"),
                miles(self.tamano)
            )
        };
        self.narracion.agregar(
            Pizarra::nueva()
                .variable("numero", numero)
                .variable("vecinos hallados", miles(self.numeros.len()))
                .lista("Vecinos hallados (los últimos)", elementos, &resaltar)
                .paso("visitar(numero);", nota),
        );
    }

    fn recorrer(&mut self, numero: u64, desde: usize, cambios_restantes: usize) {
        self.visitar(numero);
        if cambios_restantes == 0 {
            return;
        }
        for posicion in desde..self.k {
            let peso = 1u64 << (2 * (self.k - 1 - posicion));
            let original = numero / peso % 4;
            for base in 0..4 {
                if base != original {
                    let vecino = numero - original * peso + base * peso;
                    if self.detalles > 0 && self.numeros.len() < 8 {
                        self.detalles -= 1;
                        let antes = self.nombre(numero);
                        let despues = self.nombre(vecino);
                        self.narracion.agregar(
                            Pizarra::nueva()
                                .variable("numero", numero)
                                .variable("posicion", posicion)
                                .variable("peso", peso)
                                .variable("original", original)
                                .variable("base", base)
                                .variable("vecino", vecino)
                                .paso(
                                    "let vecino = numero - original * peso + base * peso;",
                                    format!(
                                        "Cambiamos la letra de la posición {posicion}: «{}» pasa a ser «{}». En números, a {numero} le restamos lo que valía la letra vieja ({original} × {peso}) y le sumamos lo que vale la nueva ({base} × {peso}): sale {vecino}, que es «{despues}» (antes «{antes}»).",
                                        letra_de(replicacion::BASES[original as usize]),
                                        letra_de(replicacion::BASES[base as usize])
                                    ),
                                ),
                        );
                    }
                    self.recorrer(vecino, posicion + 1, cambios_restantes - 1);
                }
            }
            if cambios_restantes == self.d {
                let (elementos, resaltar) = self.lista_reciente();
                self.narracion.agregar(
                    Pizarra::nueva()
                        .variable("posicion", posicion)
                        .variable("vecinos hallados", miles(self.numeros.len()))
                        .lista("Vecinos hallados (los últimos)", elementos, &resaltar)
                        .paso(
                            "for posicion in desde..k {",
                            format!(
                                "Terminamos de probar todos los cambios en la posición {posicion} (y los que se pueden combinar con ellos más a la derecha). Llevamos {} de {}.",
                                plural(self.numeros.len(), "vecino", "vecinos"),
                                miles(self.tamano)
                            ),
                        ),
                );
            }
        }
    }
}

fn trazar_vecindario(entrada: &Entrada) -> Result<(Vec<Paso>, Vec<Vec<u8>>), String> {
    let patron = entrada.adn("patron")?;
    if patron.len() > MAXIMO_K_NUMERICO {
        return Err(format!(
            "El patrón puede tener hasta {MAXIMO_K_NUMERICO} letras, y tiene {}.",
            miles(patron.len())
        ));
    }
    let d = entrada.natural("d")?;
    replicacion::vecindario(&patron, d)?;
    let k = patron.len();
    let numero = replicacion::numero_de_patron(&patron)?;
    let tamano = replicacion::tamano_del_vecindario(k, d) as usize;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        Pizarra::nueva()
            .variable("patron", letras(&patron))
            .variable("d", d)
            .variable("k", k)
            .variable("numero", numero)
            .secuencia(&patron, &[])
            .paso(
                "let numero = numero_de_patron(patron)?;",
                format!(
                    "Queremos todos los patrones que difieren de «{}» en {} como mucho. Primero convertimos el patrón en un número en base 4 (A=0, C=1, G=2, T=3): «{}» es el {}. Cambiar una letra será tan simple como sumar y restar.",
                    letras(&patron),
                    plural(d, "letra", "letras"),
                    letras(&patron),
                    miles(numero as usize)
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", k)
            .variable("d", d)
            .variable("tamano", miles(tamano))
            .paso(
                "let tamano = tamano_del_vecindario(k, d);",
                format!(
                    "Antes de empezar contamos cuántos vecinos habrá: para cada número de cambios (0, 1, … {d}) se cuentan las formas de elegir qué letras cambiar y las 3 letras distintas que puede tomar cada una. En total son {} patrones.",
                    miles(tamano)
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("numero", numero)
            .variable("k", k)
            .variable("d", d)
            .paso(
                "recorrer_vecinos(numero, k, 0, d, &mut |vecino| numeros.push(vecino));",
                "Ahora recorremos las posibilidades. La función se llama a sí misma: visita el patrón actual, y luego, para cada posición y cada letra distinta, hace un cambio y vuelve a empezar con un cambio menos permitido.",
            ),
    );
    let mut recorrido = RecorridoDeVecinos {
        narracion,
        patron: &patron,
        k,
        d,
        tamano,
        numeros: Vec::new(),
        detalles: 14,
    };
    recorrido.recorrer(numero, 0, d);
    let mut narracion = recorrido.narracion;
    let mut numeros = recorrido.numeros;
    narracion.agregar(
        Pizarra::nueva()
            .variable("numeros", format!("{} números", miles(numeros.len())))
            .paso(
                "numeros.sort();",
                format!(
                    "La búsqueda anotó los {} números en el orden en que los visitó. Los ordenamos de menor a mayor, que es lo mismo que el orden alfabético de los patrones.",
                    miles(numeros.len())
                ),
            ),
    );
    numeros.sort();
    let mut vecinos = Vec::with_capacity(numeros.len());
    for &vecino in &numeros {
        vecinos.push(replicacion::patron_de_numero(vecino, k)?);
    }
    let elementos: Vec<String> = vecinos
        .iter()
        .take(24)
        .map(|vecino| letras(vecino))
        .collect();
    let propio: Vec<usize> = vecinos
        .iter()
        .take(24)
        .enumerate()
        .filter(|(_, vecino)| **vecino == patron)
        .map(|(lugar, _)| lugar)
        .collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("vecinos", miles(vecinos.len()))
            .lista("Vecinos en orden alfabético (los primeros)", elementos.clone(), &propio)
            .paso(
                "vecinos.push(patron_de_numero(vecino, k)?);",
                "Por último, cada número se convierte otra vez en letras (el camino inverso de «De patrón a número»).",
            ),
    );
    let ultimo = Pizarra::nueva()
        .variable("vecinos", miles(vecinos.len()))
        .lista("Vecinos en orden alfabético (los primeros)", elementos, &propio)
        .paso(
            "Ok(vecinos)",
            format!(
                "Terminamos: hay {} patrones a {} como mucho de «{}», incluido él mismo. En «Qué significa» ves cuántos hay con cada número de diferencias.",
                miles(vecinos.len()),
                plural(d, "letra", "letras"),
                letras(&patron)
            ),
        );
    Ok((narracion.terminar(ultimo), vecinos))
}

fn pasos_vecindario(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    Ok(trazar_vecindario(entrada)?.0)
}

fn trazar_buscar_ori(
    entrada: &Entrada,
) -> Result<(Vec<Paso>, replicacion::OrigenProbable), String> {
    let genoma = entrada.secuencia("genoma")?;
    let ancho = entrada.natural("ancho")?;
    if ancho == 0 || ancho > genoma.len() {
        return Err(format!(
            "El ancho de la ventana debe estar entre 1 y el largo del genoma ({} pb).",
            miles(genoma.len())
        ));
    }
    let k = leer_k(entrada, "k", ancho, MAXIMO_K_ENUMERABLE)?;
    let d = leer_d(entrada, "d", k, MAXIMO_D)?;
    revisar_animable(genoma.len())?;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        Pizarra::nueva()
            .variable("genoma", vista(&genoma))
            .variable("ancho", miles(ancho))
            .variable("k", k)
            .variable("d", d)
            .secuencia(&genoma, &[])
            .paso(
                "if ancho < k || ancho > genoma.len() {",
                format!(
                    "Buscamos el origen de replicación de un genoma de {} letras, en tres etapas: 1) el sesgo GC señala una región candidata, 2) tomamos una ventana de {} letras ahí, 3) en esa ventana buscamos los k-meros de {k} letras más frecuentes, con hasta {} y en las dos hebras. Los datos tienen sentido, así que empezamos.",
                    miles(genoma.len()),
                    miles(ancho),
                    plural(d, "diferencia", "diferencias")
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("sesgo_actual", 0)
            .variable("minimo", 0)
            .variable("posiciones", "[0]")
            .secuencia(&genoma, &[])
            .paso(
                "let mut sesgo_actual = 0;",
                "Etapa 1: el sesgo GC. Recorremos el genoma sumando 1 por cada G y restando 1 por cada C; las demás letras no cambian nada. Antes de leer ninguna base el sesgo vale 0, y la posición 0 es el primer mínimo.",
            ),
    );
    let controles = puntos_de_control(genoma.len());
    let mut sesgo_actual: i64 = 0;
    let mut minimo: i64 = 0;
    let mut minimos: Vec<usize> = vec![0];
    let mut eventos = 0usize;
    for (indice, &base) in genoma.iter().enumerate() {
        sesgo_actual += replicacion::cambio_de_sesgo(base);
        let posicion = indice + 1;
        let mut nuevo_minimo = false;
        if sesgo_actual < minimo {
            minimo = sesgo_actual;
            minimos.clear();
            minimos.push(posicion);
            nuevo_minimo = true;
        } else if sesgo_actual == minimo {
            minimos.push(posicion);
        }
        let marca = [(indice, indice + 1, Resalte::Actual)];
        if nuevo_minimo && eventos < 6 {
            eventos += 1;
            narracion.agregar(
                Pizarra::nueva()
                    .variable("base", letra_de(base))
                    .variable("sesgo_actual", con_signo(sesgo_actual))
                    .variable("minimo", con_signo(minimo))
                    .variable("posicion", miles(posicion))
                    .secuencia(&genoma, &marca)
                    .paso(
                        "minimo = sesgo_actual;",
                        format!(
                            "En la posición {} el sesgo baja a {}, más bajo que nunca: es el nuevo mínimo. Borramos los mínimos anteriores y guardamos esta posición.",
                            miles(posicion),
                            con_signo(sesgo_actual)
                        ),
                    ),
            );
        } else if indice < 3 {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("base", letra_de(base))
                    .variable("sesgo_actual", con_signo(sesgo_actual))
                    .variable("minimo", con_signo(minimo))
                    .secuencia(&genoma, &marca)
                    .paso(
                        "sesgo_actual += cambio_de_sesgo(base);",
                        format!(
                            "Leemos «{}»: {}. El sesgo queda en {}.",
                            letra_de(base),
                            match base {
                                b'G' | b'g' => "una G suma 1",
                                b'C' | b'c' => "una C resta 1",
                                _ => "no es ni G ni C, así que no cambia nada",
                            },
                            con_signo(sesgo_actual)
                        ),
                    ),
            );
        } else if controles.contains(&posicion) {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("base", letra_de(base))
                    .variable("sesgo_actual", con_signo(sesgo_actual))
                    .variable("minimo", con_signo(minimo))
                    .secuencia(&genoma, &marca)
                    .paso(
                        "sesgo_actual += cambio_de_sesgo(base);",
                        format!(
                            "Vamos por la posición {} de {}: el sesgo vale {} y el mínimo hasta ahora es {}.",
                            miles(posicion),
                            miles(genoma.len()),
                            con_signo(sesgo_actual),
                            con_signo(minimo)
                        ),
                    ),
            );
        }
    }
    narracion.agregar(
        Pizarra::nueva()
            .variable("minimo", con_signo(minimo))
            .variable("minimos", lista_corta(&minimos))
            .paso(
                "let minimos = posiciones_de_sesgo_minimo(genoma);",
                format!(
                    "Terminamos de recorrer el genoma. El sesgo más bajo es {} y se alcanza en {}. La primera de esas posiciones, {}, es el candidato a origen de replicación.",
                    con_signo(minimo),
                    plural(minimos.len(), "posición", "posiciones"),
                    miles(minimos[0])
                ),
            ),
    );
    let inicio = minimos[0].min(genoma.len() - ancho);
    narracion.agregar(
        Pizarra::nueva()
            .variable("minimos[0]", miles(minimos[0]))
            .variable("ancho", miles(ancho))
            .variable("inicio", miles(inicio))
            .secuencia(&genoma, &[(inicio, inicio + 1, Resalte::Actual)])
            .paso(
                "let inicio = minimos[0].min(genoma.len() - ancho);",
                if inicio == minimos[0] {
                    format!(
                        "Etapa 2: la ventana empieza justo en el primer mínimo, la posición {}, porque desde ahí caben {} letras sin salirse del genoma.",
                        miles(inicio),
                        miles(ancho)
                    )
                } else {
                    format!(
                        "Etapa 2: desde el mínimo ({}) no caben {} letras sin salirse del genoma, así que la ventana se corre hacia atrás y empieza en {}.",
                        miles(minimos[0]),
                        miles(ancho),
                        miles(inicio)
                    )
                },
            ),
    );
    let fin = inicio + ancho;
    narracion.agregar(
        Pizarra::nueva()
            .variable("inicio", miles(inicio))
            .variable("ancho", miles(ancho))
            .variable("fin", miles(fin))
            .secuencia(&genoma, &[(inicio, fin, Resalte::Hallado)])
            .paso(
                "let fin = inicio + ancho;",
                format!(
                    "La ventana va de la posición {} a la {} (sin incluir la última): son {} letras. Solo mirando este trozo buscaremos la señal.",
                    miles(inicio),
                    miles(fin),
                    miles(ancho)
                ),
            ),
    );
    let ventana = &genoma[inicio..fin];
    let cuentas = replicacion::cuentas_en_ambas_hebras(ventana, k, d)?;
    let exactas_aqui = replicacion::arreglo_de_frecuencias(ventana, k)?;
    let otra_hebra =
        replicacion::arreglo_de_frecuencias(&secuencia::complemento_inverso(ventana), k)?;
    let mut ambas = exactas_aqui.clone();
    for (cuenta, otra) in ambas.iter_mut().zip(&otra_hebra) {
        *cuenta += otra;
    }
    let distintos = ambas.iter().filter(|&&veces| veces > 0).count();
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", k)
            .variable("d", d)
            .tabla(
                "Los k-meros más repetidos de la ventana (exactos)",
                &["k-mero", "veces"],
                tabla_de_mejores(&exactas_aqui, None, k, 5)?,
                &[(0, 1)],
            )
            .paso(
                "let cuentas = cuentas_en_ambas_hebras(&genoma[inicio..fin], k, d)?;",
                format!(
                    "Etapa 3: en la ventana de {} letras contamos los k-meros de {k} letras. Primero, tal cual aparecen en esta hebra.",
                    miles(ancho)
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", k)
            .tabla(
                "Los k-meros más repetidos en la otra hebra",
                &["k-mero", "veces"],
                tabla_de_mejores(&otra_hebra, None, k, 5)?,
                &[(0, 1)],
            )
            .paso(
                "let otra_hebra = arreglo_de_frecuencias(&complemento_inverso(texto), k)?;",
                "Luego hacemos la misma cuenta en la otra hebra, que es el complemento inverso de la ventana: así una señal escrita en cualquiera de las dos hebras se suma.",
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", k)
            .variable("exactas", format!("{} k-meros distintos", miles(distintos)))
            .tabla(
                "Suma de las dos hebras",
                &["k-mero", "veces"],
                tabla_de_mejores(&ambas, None, k, 5)?,
                &[(0, 1)],
            )
            .paso(
                "*cuenta += otra;",
                "Sumamos las dos cuentas casilla por casilla: cada k-mero y su complemento inverso quedan con el mismo total.",
            ),
    );
    let visitas = distintos as f64 * replicacion::tamano_del_vecindario(k, d);
    narracion.agregar(
        Pizarra::nueva()
            .variable("distintos", miles(distintos))
            .variable("visitas", miles(visitas as usize))
            .tabla(
                "Los k-meros con más cuenta (con errores)",
                &["k-mero", "exactas, dos hebras", "con errores"],
                tabla_de_mejores(&cuentas, Some(&ambas[..]), k, 6)?,
                &[(0, 2)],
            )
            .paso(
                "let visitas = distintos as f64 * tamano_del_vecindario(k, d);",
                format!(
                    "Ahora cada uno de los {} k-meros distintos reparte su cuenta entre sus vecinos, los que difieren en {} como mucho: {} visitas en total. Así un k-mero cuenta también las veces que aparece con un error.",
                    miles(distintos),
                    plural(d, "letra", "letras"),
                    miles(visitas as usize)
                ),
            ),
    );
    let veces = cuentas.iter().copied().max().unwrap_or(0);
    narracion.agregar(
        Pizarra::nueva()
            .variable("veces", miles(veces as usize))
            .paso(
                "let veces = cuentas.iter().copied().max().unwrap_or(0);",
                format!(
                    "La cuenta más alta con errores y en las dos hebras es {}.",
                    miles(veces as usize)
                ),
            ),
    );
    let candidatos = replicacion::los_mas_frecuentes(&cuentas, k);
    let nombres: Vec<String> = candidatos.iter().take(20).map(|p| letras(p)).collect();
    let resaltar_todos: Vec<usize> = (0..nombres.len()).collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("veces", miles(veces as usize))
            .variable("candidatos", miles(candidatos.len()))
            .lista("Candidatos a caja DnaA", nombres.clone(), &resaltar_todos)
            .paso(
                "let candidatos = los_mas_frecuentes(&cuentas, k);",
                format!(
                    "Nos quedamos con los k-meros que alcanzan esa cuenta: {}.",
                    plural(candidatos.len(), "candidato", "candidatos")
                ),
            ),
    );
    let ultimo = Pizarra::nueva()
        .variable("inicio", miles(inicio))
        .variable("fin", miles(fin))
        .variable("veces", miles(veces as usize))
        .lista("Candidatos a caja DnaA", nombres, &resaltar_todos)
        .secuencia(&genoma, &[(inicio, fin, Resalte::Hallado)])
        .paso(
            "Ok(OrigenProbable {",
            format!(
                "Terminamos: la región candidata es la ventana de la posición {} a la {}, y sus k-meros más frecuentes aparecen {} veces. Esos son los candidatos a caja DnaA, la señal que reconoce la proteína que abre el ADN. En «Qué significa» ves qué tan fuerte es.",
                miles(inicio),
                miles(fin),
                miles(veces as usize)
            ),
        );
    let origen = replicacion::OrigenProbable {
        minimos,
        inicio,
        fin,
        candidatos,
        veces,
    };
    Ok((narracion.terminar(ultimo), origen))
}

fn pasos_buscar_ori(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    Ok(trazar_buscar_ori(entrada)?.0)
}

#[cfg(test)]
mod pruebas {
    use super::*;
    use std::time::Instant;

    fn genoma_de_prueba(longitud: usize, semilla: u64) -> String {
        let mut estado = semilla;
        let mut texto = String::with_capacity(longitud);
        for _ in 0..longitud {
            estado = estado
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            texto.push(['A', 'C', 'G', 'T'][(estado >> 33) as usize % 4]);
        }
        texto
    }

    fn interpretacion(salida: &Salida) -> String {
        salida.interpretacion.join(" ")
    }

    #[test]
    fn buscar_ori_con_ventana_enorme_no_se_cuelga() {
        let genoma = genoma_de_prueba(300_000, 7);
        let largo = genoma.len().to_string();
        let entrada = Entrada::nueva(&[
            ("genoma", genoma.as_str()),
            ("ancho", largo.as_str()),
            ("k", "1"),
            ("d", "1"),
        ]);
        let inicio = Instant::now();
        let salida = buscar_ori(&entrada).unwrap();
        assert!(inicio.elapsed().as_secs() < 20);
        assert!(!salida.interpretacion.is_empty());
        assert!(interpretacion(&salida).contains("mucho más grande que un origen típico"));
    }

    #[test]
    fn minimos_seguidos_y_minimos_juntos_se_explican_distinto() {
        let seguidos = format!("C{}GG", "A".repeat(5));
        let salida = sesgo_minimo(&Entrada::nueva(&[("genoma", seguidos.as_str())])).unwrap();
        assert!(interpretacion(&salida).contains("posiciones seguidas"));

        let juntos = format!("CGCGG{}", "A".repeat(1000));
        let salida = sesgo_minimo(&Entrada::nueva(&[("genoma", juntos.as_str())])).unwrap();
        let texto = interpretacion(&salida);
        assert!(texto.contains("muy juntas"));
        assert!(!texto.contains("no hay ni G ni C"));
        assert!(!texto.contains("solo hay A o T"));
    }

    #[test]
    fn un_doce_mero_repetido_tres_veces_no_es_senal_debil() {
        let mut fondo = genoma_de_prueba(1000, 11).into_bytes();
        let caja = b"ACGGTCAGCTTA";
        for inicio in [100, 400, 800] {
            fondo[inicio..inicio + caja.len()].copy_from_slice(caja);
        }
        let genoma = String::from_utf8(fondo).unwrap();
        let entrada = Entrada::nueva(&[
            ("genoma", genoma.as_str()),
            ("ancho", "1000"),
            ("k", "12"),
            ("d", "0"),
        ]);
        let salida = buscar_ori(&entrada).unwrap();
        assert!(!interpretacion(&salida).contains("señal débil"));
    }

    fn entrada_de_ejemplo(id: &str) -> Entrada {
        Entrada::desde_json(crate::ejemplo(id).unwrap()).unwrap()
    }

    fn con_campo(id: &str, campo: &str, valor: &str) -> Entrada {
        let mut ejemplo = crate::ejemplo(id).unwrap();
        ejemplo[campo] = serde_json::Value::String(valor.to_string());
        Entrada::desde_json(ejemplo).unwrap()
    }

    fn revisar_lineas(id: &str, pasos: &[Paso]) {
        let fragmentos = crate::codigo(id).unwrap();
        let codigo: String = fragmentos
            .iter()
            .map(|fragmento| fragmento.fuente.as_str())
            .collect::<Vec<_>>()
            .join("\n");
        assert!(!pasos.is_empty() && pasos.len() <= 400, "{id}");
        for paso in pasos {
            assert!(
                codigo
                    .lines()
                    .any(|linea| linea.contains(paso.linea.as_str())),
                "{id}: la línea «{}» no está en el código",
                paso.linea
            );
            assert!(!paso.nota.is_empty(), "{id}: nota vacía");
        }
    }

    #[test]
    fn los_barridos_dan_lo_mismo_que_el_nucleo() {
        let texto = entrada_de_ejemplo("replicacion.contar_patron")
            .secuencia("texto")
            .unwrap();
        let (_, hallazgos) = barrer(&barrido_de_contar(&texto, b"GATC"));
        assert_eq!(hallazgos.len(), replicacion::contar_patron(&texto, b"GATC"));
        let (pasos, posiciones) = barrer(&barrido_de_posiciones(&texto, b"TTATCCACA"));
        assert_eq!(
            posiciones,
            replicacion::posiciones_del_patron(&texto, b"TTATCCACA")
        );
        revisar_lineas("replicacion.posiciones_del_patron", &pasos);
        for d in 0..=2 {
            let (_, hallazgos) = barrer(&barrido_de_contar_aproximado(&texto, b"TTATCCACA", d));
            assert_eq!(
                hallazgos.len(),
                replicacion::contar_aproximado(&texto, b"TTATCCACA", d)
            );
            let (pasos, posiciones) =
                barrer(&barrido_de_posiciones_aproximadas(&texto, b"TTATCCACA", d));
            assert_eq!(
                posiciones,
                replicacion::posiciones_aproximadas(&texto, b"TTATCCACA", d)
            );
            revisar_lineas("replicacion.posiciones_aproximadas", &pasos);
        }
        let (_, ninguno) = barrer(&barrido_de_contar(b"ACGT", b"ACGTACGT"));
        assert!(ninguno.is_empty());
    }

    #[test]
    fn el_complemento_y_el_gc_dan_lo_mismo_que_el_nucleo() {
        let entrada = entrada_de_ejemplo("replicacion.complemento_inverso");
        let (pasos, resultado) = trazar_complemento_inverso(&entrada).unwrap();
        assert_eq!(
            resultado,
            secuencia::complemento_inverso(&entrada.secuencia("adn").unwrap())
        );
        revisar_lineas("replicacion.complemento_inverso", &pasos);
        let entrada = con_campo("replicacion.complemento_inverso", "adn", "acgtNUu");
        let (pasos, resultado) = trazar_complemento_inverso(&entrada).unwrap();
        assert_eq!(
            resultado,
            secuencia::complemento_inverso(&entrada.secuencia("adn").unwrap())
        );
        revisar_lineas("replicacion.complemento_inverso", &pasos);
        let entrada = entrada_de_ejemplo("replicacion.contenido_gc");
        let (pasos, fraccion) = trazar_contenido_gc(&entrada).unwrap();
        assert_eq!(
            fraccion,
            secuencia::contenido_gc(&entrada.secuencia("adn").unwrap())
        );
        revisar_lineas("replicacion.contenido_gc", &pasos);
    }

    #[test]
    fn hamming_y_los_numeros_dan_lo_mismo_que_el_nucleo() {
        let entrada = entrada_de_ejemplo("replicacion.distancia_hamming");
        let (pasos, distancia) = trazar_hamming(&entrada).unwrap();
        assert_eq!(
            distancia,
            replicacion::distancia_hamming(
                &entrada.secuencia("primera").unwrap(),
                &entrada.secuencia("segunda").unwrap()
            )
            .unwrap()
        );
        revisar_lineas("replicacion.distancia_hamming", &pasos);
        let entrada = Entrada::nueva(&[("primera", "ACGTACGT"), ("segunda", "ACGAACGA")]);
        assert_eq!(trazar_hamming(&entrada).unwrap().1, 2);
        let entrada = entrada_de_ejemplo("replicacion.numero_de_patron");
        let (pasos, numero) = trazar_numero_de_patron(&entrada).unwrap();
        assert_eq!(numero, replicacion::numero_de_patron(b"GATTACA").unwrap());
        revisar_lineas("replicacion.numero_de_patron", &pasos);
        let entrada = entrada_de_ejemplo("replicacion.patron_de_numero");
        let (pasos, patron) = trazar_patron_de_numero(&entrada).unwrap();
        assert_eq!(patron, replicacion::patron_de_numero(9000, 7).unwrap());
        revisar_lineas("replicacion.patron_de_numero", &pasos);
    }

    #[test]
    fn los_conteos_dan_lo_mismo_que_el_nucleo() {
        let entrada = entrada_de_ejemplo("replicacion.arreglo_de_frecuencias");
        let (pasos, arreglo) = trazar_arreglo(&entrada).unwrap();
        let texto = entrada.secuencia("texto").unwrap();
        assert_eq!(
            arreglo,
            replicacion::arreglo_de_frecuencias(&texto, 2).unwrap()
        );
        revisar_lineas("replicacion.arreglo_de_frecuencias", &pasos);
        let entrada = entrada_de_ejemplo("replicacion.palabras_frecuentes");
        let (pasos, frecuentes) = trazar_palabras_frecuentes(&entrada).unwrap();
        let texto = entrada.adn("texto").unwrap();
        assert_eq!(frecuentes, replicacion::palabras_frecuentes(&texto, 4));
        revisar_lineas("replicacion.palabras_frecuentes", &pasos);
        let entrada = con_campo("replicacion.palabras_frecuentes", "k", "20");
        let (pasos, frecuentes) = trazar_palabras_frecuentes(&entrada).unwrap();
        assert_eq!(frecuentes, replicacion::palabras_frecuentes(&texto, 20));
        revisar_lineas("replicacion.palabras_frecuentes", &pasos);
        let entrada = con_campo("replicacion.palabras_frecuentes", "k", "3");
        let (_, frecuentes) = trazar_palabras_frecuentes(&entrada).unwrap();
        assert_eq!(frecuentes, replicacion::palabras_frecuentes(&texto, 3));
    }

    #[test]
    fn los_grumos_y_los_vecinos_dan_lo_mismo_que_el_nucleo() {
        let genoma = genoma_de_prueba(5000, 3);
        let entrada = Entrada::nueva(&[
            ("genoma", genoma.as_str()),
            ("k", "4"),
            ("ventana", "60"),
            ("veces", "4"),
        ]);
        let (pasos, grumos) = trazar_grumos(&entrada).unwrap();
        assert_eq!(
            grumos,
            replicacion::encontrar_grumos(genoma.as_bytes(), 4, 60, 4).unwrap()
        );
        assert!(!grumos.is_empty());
        revisar_lineas("replicacion.grumos", &pasos);
        let entrada = entrada_de_ejemplo("replicacion.vecindario");
        let (pasos, vecinos) = trazar_vecindario(&entrada).unwrap();
        assert_eq!(vecinos, replicacion::vecindario(b"TTATCCACA", 1).unwrap());
        revisar_lineas("replicacion.vecindario", &pasos);
        let entrada = con_campo("replicacion.vecindario", "d", "3");
        let (_, vecinos) = trazar_vecindario(&entrada).unwrap();
        assert_eq!(vecinos, replicacion::vecindario(b"TTATCCACA", 3).unwrap());
    }

    #[test]
    fn los_frecuentes_y_el_origen_dan_lo_mismo_que_el_nucleo() {
        for dos_hebras in [false, true] {
            let id = if dos_hebras {
                "replicacion.frecuentes_con_complementos"
            } else {
                "replicacion.frecuentes_con_errores"
            };
            let entrada = entrada_de_ejemplo(id);
            let (pasos, frecuentes) = trazar_frecuentes(&entrada, dos_hebras).unwrap();
            let texto = entrada.secuencia("texto").unwrap();
            let esperado = if dos_hebras {
                replicacion::palabras_frecuentes_con_errores_y_complementos(&texto, 9, 1)
            } else {
                replicacion::palabras_frecuentes_con_errores(&texto, 9, 1)
            }
            .unwrap();
            assert_eq!(frecuentes, esperado);
            revisar_lineas(id, &pasos);
        }
        let genoma = genoma_de_prueba(20_000, 5);
        let entrada = Entrada::nueva(&[
            ("genoma", genoma.as_str()),
            ("ancho", "300"),
            ("k", "5"),
            ("d", "1"),
        ]);
        let (pasos, origen) = trazar_buscar_ori(&entrada).unwrap();
        assert_eq!(
            origen,
            replicacion::buscar_ori(genoma.as_bytes(), 300, 5, 1).unwrap()
        );
        revisar_lineas("replicacion.buscar_ori", &pasos);
    }
}
