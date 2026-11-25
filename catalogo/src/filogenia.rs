use crate::formato::{decimal, miles, porcentaje};
use crate::pizarra::{Marca, Narracion, Pizarra};
use crate::{Algoritmo, Bloque, Campo, Entrada, FilaAlineada, Paso, Punto, Salida, Tipo};
use nucleo::filogenia::{self, Arbol, Parsimonia};
use serde_json::{Value, json};
use std::collections::HashSet;

const MAXIMO_DE_ESPECIES: usize = 200;
const MAXIMO_DE_POSICIONES: usize = 2000;
const MAXIMO_PARA_ANIMAR: usize = 30;
const MAXIMO_DE_MOVIMIENTOS: usize = 200;
const MAXIMO_EN_TABLA: usize = 20;
const MAXIMO_EN_ALINEAMIENTO: usize = 120;
const INICIOS: &[&str] = &["unión de vecinos", "al azar"];

const ARBOL_CON_PESOS: &str = "Nubaria->n1:2
Ferolis->n1:7
n1->n2:2
Cardenix->n2:3
n2->n3:2
Tavira->n3:2
n3->n4:3
Ostrelo->n4:4
Pumbra->n4:1";

const MATRIZ_ADITIVA: &str = "Nubaria Ferolis Cardenix Tavira Ostrelo Pumbra
0 5 11 9 14 11.5
5 0 12 10 15 12.5
11 12 0 8 13 10.5
9 10 8 0 9 6.5
14 15 13 9 0 5.5
11.5 12.5 10.5 6.5 5.5 0";

const MATRIZ_CON_RUIDO: &str = "Nubaria Ferolis Cardenix Tavira Ostrelo Pumbra
0 8 5 6 8.5 7
8 0 11 11.5 14.5 13
5 11 0 4.5 7 6
6 11.5 4.5 0 6 4.5
8.5 14.5 7 6 0 3.5
7 13 6 4.5 3.5 0";

const SECUENCIAS_ALINEADAS: &str = "Nubaria=TATGCTCATGCAACTCAAAACCATGTCCGTAATGTAGGCGAAAAAGTA
Ferolis=TATCATCATGCACTTCAAAACCATGCGCGTAATGTACGCGGAAAATTA
Cardenix=CTTCCTGCTGCAATTCAAAA---TGTACGTAATGTAGGCGAAATAATA
Tavira=TTTCCTCCTGCTATTCAAAACCATGTACATAAGGTAGGCGAAAAAGTA
Ostrelo=ATTCCTCCCGCTACTCACAACCATGTTCGTAATGTAGGCGAGAAAGTA
Pumbra=TTTCCTCCCGCCATTCAAAACCATGTCCGTAATGTAGGCGAGAAAGTA";

const ARBOL_CON_RAIZ: &str = "((Nubaria,Ferolis),(Cardenix,(Tavira,Ostrelo)));";

const ARBOL_SIN_RAIZ: &str = "(Nubaria,Ferolis,(Cardenix,(Tavira,Ostrelo)));";

const SECUENCIAS_DE_HOJAS: &str = "Nubaria=GGAGTGGCCGTACTAAATGATGAC
Ferolis=GGAGTTCCCGTACTAAATGATGAT
Cardenix=GAAGTTGCCCTTCTAAATTATTAC
Tavira=GACGTTGCCGTTCTAGATTATGAC
Ostrelo=GACGTTGCCGTTCTTTATTATAAC";

const ARBOL_BINARIO: &str = "Nubaria->n1
Ferolis->n1
n1->n2
Cardenix->n2
n2->n3
Tavira->n3
n3->n4
Ostrelo->n4
Pumbra->n4";

const OCHO_ESPECIES: &str = "Nubaria=AAGTATGTCTCAATACGTGACTAAAGACAGGCAACACGAGGGTCCGATTA
Ferolis=ACGTATTTTACAATACGTGACTAAAGACAGGCAACACGAGGATCCGATTA
Cardenix=AAGTATGTTTCAATAAGTGACTAAAGACAGGCGACACGAGTGCCCGATTA
Tavira=AAGTATGCGTCGAAAGGTGACTAAAGACAGGCATCGTTAGGCTCCGATTA
Ostrelo=GAGTATGCGTCGAAAGGTGACTCAAGACAGGCAACGTGAGGCTCCGATTA
Pumbra=AAGAATGTTTCGATAGGTGAATAAAAACCGGCAACGAGAGGCTCCGATTA
Silvana=AAAAATGTTTCGATAGGTGCCTAAAGACTGGCAACGAGAGGCTCCTATTA
Quirela=AAAAATGTTTCGATAGGTGCCTGAAGACAGGTAACGAGAGGCTCCGATGA";

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![
        Algoritmo {
            id: "filogenia.distancias_entre_hojas",
            capitulo: "filogenia",
            nombre: "Distancias entre hojas",
            resumen: "Suma las longitudes del camino entre cada par de hojas de un árbol con pesos.",
            campos: &[Campo {
                id: "arbol",
                etiqueta: "Árbol con pesos (una arista por línea: a->b:longitud)",
                tipo: Tipo::Texto,
                ejemplo: ARBOL_CON_PESOS,
            }],
            funcion: "filogenia::distancias_entre_hojas",
            ejecutar: distancias_entre_hojas,
            pasos: Some(pasos_distancias_entre_hojas),
        },
        Algoritmo {
            id: "filogenia.longitud_de_rama",
            capitulo: "filogenia",
            nombre: "Longitud de rama",
            resumen: "Cuánto mide la rama propia de una hoja, calculado solo desde la matriz de distancias.",
            campos: &[
                Campo {
                    id: "matriz",
                    etiqueta: "Matriz de distancias (primera línea opcional con nombres)",
                    tipo: Tipo::Matriz,
                    ejemplo: MATRIZ_ADITIVA,
                },
                Campo {
                    id: "hoja",
                    etiqueta: "Hoja (nombre o número contando desde 0)",
                    tipo: Tipo::Texto,
                    ejemplo: "Cardenix",
                },
            ],
            funcion: "filogenia::longitud_de_rama",
            ejecutar: longitud_de_rama,
            pasos: Some(pasos_longitud_de_rama),
        },
        Algoritmo {
            id: "filogenia.filogenia_aditiva",
            capitulo: "filogenia",
            nombre: "Filogenia aditiva",
            resumen: "Reconstruye el único árbol que reproduce exactamente una matriz aditiva, hoja por hoja.",
            campos: &[Campo {
                id: "matriz",
                etiqueta: "Matriz de distancias (primera línea opcional con nombres)",
                tipo: Tipo::Matriz,
                ejemplo: MATRIZ_ADITIVA,
            }],
            funcion: "filogenia::filogenia_aditiva",
            ejecutar: filogenia_aditiva,
            pasos: Some(pasos_filogenia_aditiva),
        },
        Algoritmo {
            id: "filogenia.upgma",
            capitulo: "filogenia",
            nombre: "UPGMA",
            resumen: "Árbol con raíz que une siempre los dos grupos más cercanos, suponiendo un reloj molecular.",
            campos: &[Campo {
                id: "matriz",
                etiqueta: "Matriz de distancias (primera línea opcional con nombres)",
                tipo: Tipo::Matriz,
                ejemplo: MATRIZ_CON_RUIDO,
            }],
            funcion: "filogenia::upgma",
            ejecutar: upgma,
            pasos: Some(pasos_upgma),
        },
        Algoritmo {
            id: "filogenia.union_de_vecinos",
            capitulo: "filogenia",
            nombre: "Unión de vecinos",
            resumen: "Árbol sin raíz que une los pares que minimizan la matriz D*, sin suponer reloj molecular.",
            campos: &[Campo {
                id: "matriz",
                etiqueta: "Matriz de distancias (primera línea opcional con nombres)",
                tipo: Tipo::Matriz,
                ejemplo: MATRIZ_CON_RUIDO,
            }],
            funcion: "filogenia::union_de_vecinos",
            ejecutar: union_de_vecinos,
            pasos: Some(pasos_union_de_vecinos),
        },
        Algoritmo {
            id: "filogenia.distancias_evolutivas",
            capitulo: "filogenia",
            nombre: "Distancias p y Jukes–Cantor",
            resumen: "Matriz de distancias entre secuencias alineadas: diferencias observadas y corregidas por cambios ocultos.",
            campos: &[Campo {
                id: "secuencias",
                etiqueta: "Secuencias alineadas (nombre=SECUENCIA por línea, o FASTA)",
                tipo: Tipo::Texto,
                ejemplo: SECUENCIAS_ALINEADAS,
            }],
            funcion: "filogenia::matriz_jukes_cantor",
            ejecutar: distancias_evolutivas,
            pasos: Some(pasos_distancias_evolutivas),
        },
        Algoritmo {
            id: "filogenia.arbol_desde_secuencias",
            capitulo: "filogenia",
            nombre: "Árbol desde secuencias",
            resumen: "De secuencias alineadas a un árbol: distancias de Jukes–Cantor y luego unión de vecinos.",
            campos: &[Campo {
                id: "secuencias",
                etiqueta: "Secuencias alineadas (nombre=SECUENCIA por línea, o FASTA)",
                tipo: Tipo::Texto,
                ejemplo: SECUENCIAS_ALINEADAS,
            }],
            funcion: "filogenia::arbol_desde_secuencias",
            ejecutar: arbol_desde_secuencias,
            pasos: Some(pasos_arbol_desde_secuencias),
        },
        Algoritmo {
            id: "filogenia.parsimonia_pequena",
            capitulo: "filogenia",
            nombre: "Parsimonia pequeña (árbol con raíz)",
            resumen: "Letras de los antepasados que explican las hojas con el menor número de mutaciones.",
            campos: &[
                Campo {
                    id: "arbol",
                    etiqueta: "Árbol con raíz en formato Newick",
                    tipo: Tipo::Texto,
                    ejemplo: ARBOL_CON_RAIZ,
                },
                Campo {
                    id: "secuencias",
                    etiqueta: "Secuencias de las hojas (nombre=SECUENCIA por línea, o FASTA)",
                    tipo: Tipo::Texto,
                    ejemplo: SECUENCIAS_DE_HOJAS,
                },
            ],
            funcion: "filogenia::parsimonia_pequena",
            ejecutar: parsimonia_pequena,
            pasos: Some(pasos_parsimonia_pequena),
        },
        Algoritmo {
            id: "filogenia.parsimonia_sin_raiz",
            capitulo: "filogenia",
            nombre: "Parsimonia pequeña (árbol sin raíz)",
            resumen: "El mismo problema en un árbol sin raíz: el puntaje no depende de dónde se ponga la raíz.",
            campos: &[
                Campo {
                    id: "arbol",
                    etiqueta: "Árbol sin raíz en formato Newick",
                    tipo: Tipo::Texto,
                    ejemplo: ARBOL_SIN_RAIZ,
                },
                Campo {
                    id: "secuencias",
                    etiqueta: "Secuencias de las hojas (nombre=SECUENCIA por línea, o FASTA)",
                    tipo: Tipo::Texto,
                    ejemplo: SECUENCIAS_DE_HOJAS,
                },
            ],
            funcion: "filogenia::parsimonia_sin_raiz",
            ejecutar: parsimonia_sin_raiz,
            pasos: Some(pasos_parsimonia_sin_raiz),
        },
        Algoritmo {
            id: "filogenia.vecinos_cercanos",
            capitulo: "filogenia",
            nombre: "Vecinos cercanos de un árbol",
            resumen: "Los dos árboles que salen de intercambiar subárboles a ambos lados de una arista interna.",
            campos: &[
                Campo {
                    id: "arbol",
                    etiqueta: "Árbol binario sin raíz (una arista por línea: a->b o a->b:longitud)",
                    tipo: Tipo::Texto,
                    ejemplo: ARBOL_BINARIO,
                },
                Campo {
                    id: "arista",
                    etiqueta: "Arista interna (dos nombres)",
                    tipo: Tipo::Texto,
                    ejemplo: "n2 n3",
                },
            ],
            funcion: "filogenia::vecinos_cercanos",
            ejecutar: vecinos_cercanos,
            pasos: Some(pasos_vecinos_cercanos),
        },
        Algoritmo {
            id: "filogenia.parsimonia_grande",
            capitulo: "filogenia",
            nombre: "Parsimonia grande por intercambio de vecinos",
            resumen: "Busca el árbol más parsimonioso moviéndose al mejor vecino cercano mientras el puntaje baje.",
            campos: &[
                Campo {
                    id: "secuencias",
                    etiqueta: "Secuencias alineadas (nombre=SECUENCIA por línea, o FASTA)",
                    tipo: Tipo::Texto,
                    ejemplo: OCHO_ESPECIES,
                },
                Campo {
                    id: "inicio",
                    etiqueta: "Árbol inicial",
                    tipo: Tipo::Opcion(INICIOS),
                    ejemplo: "al azar",
                },
                Campo {
                    id: "semilla",
                    etiqueta: "Semilla (para el árbol al azar)",
                    tipo: Tipo::Entero,
                    ejemplo: "7",
                },
                Campo {
                    id: "movimientos",
                    etiqueta: "Máximo de movimientos",
                    tipo: Tipo::Entero,
                    ejemplo: "100",
                },
            ],
            funcion: "filogenia::parsimonia_grande",
            ejecutar: parsimonia_grande,
            pasos: Some(pasos_parsimonia_grande),
        },
    ]
}

fn numero(valor: f64) -> String {
    if valor.is_infinite() {
        return "∞".to_string();
    }
    if valor.is_nan() {
        return "—".to_string();
    }
    let texto = decimal(valor, 3);
    let texto = texto.trim_end_matches('0').trim_end_matches('.');
    if texto == "-0" {
        return "0".to_string();
    }
    texto.to_string()
}

fn mutaciones(cantidad: usize) -> String {
    if cantidad == 1 {
        "1 mutación".to_string()
    } else {
        format!("{} mutaciones", miles(cantidad))
    }
}

fn redondear(valor: f64) -> f64 {
    let redondeado = (valor * 1e6).round() / 1e6;
    if redondeado == 0.0 { 0.0 } else { redondeado }
}

fn casi_cero(valor: f64) -> bool {
    valor.abs() < 1e-9
}

fn piezas(linea: &str) -> Vec<&str> {
    linea
        .split(|c: char| c.is_whitespace() || c == ',' || c == ';')
        .filter(|pieza| !pieza.is_empty())
        .collect()
}

fn recortar(texto: &str) -> String {
    let mut corto: String = texto.chars().take(40).collect();
    if texto.chars().count() > 40 {
        corto.push('…');
    }
    corto
}

fn revisar_nombres_unicos(nombres: &[String]) -> Result<(), String> {
    let mut vistos = HashSet::new();
    for nombre in nombres {
        if !filogenia::es_nombre_valido(nombre) {
            return Err(format!(
                "«{}» no es un nombre válido: usa letras, números, «_» o «.».",
                recortar(nombre)
            ));
        }
        if !vistos.insert(nombre) {
            return Err(format!("El nombre «{nombre}» está repetido."));
        }
    }
    Ok(())
}

fn leer_matriz(entrada: &Entrada, campo: &str) -> Result<(Vec<String>, Vec<Vec<f64>>), String> {
    let texto = entrada.texto(campo)?;
    let mut lineas: Vec<&str> = texto
        .lines()
        .map(str::trim)
        .filter(|linea| !linea.is_empty())
        .collect();
    let mut nombres: Option<Vec<String>> = None;
    if let Some(primera) = lineas.first() {
        let encabezado = piezas(primera);
        if encabezado.iter().any(|pieza| pieza.parse::<f64>().is_err()) {
            nombres = Some(encabezado.iter().map(|pieza| pieza.to_string()).collect());
            lineas.remove(0);
        } else if encabezado.len() == 1
            && lineas.len() > 1
            && encabezado[0].parse::<usize>() == Ok(lineas.len() - 1)
        {
            lineas.remove(0);
        }
    }
    if lineas.len() > MAXIMO_DE_ESPECIES {
        return Err(format!(
            "La matriz tiene {} filas; el máximo aquí es {MAXIMO_DE_ESPECIES} especies.",
            miles(lineas.len())
        ));
    }
    let mut filas = Vec::with_capacity(lineas.len());
    for (numero_de_fila, linea) in lineas.iter().enumerate() {
        let mut fila = Vec::new();
        for pieza in piezas(linea) {
            let valor = pieza.parse::<f64>().ok().filter(|valor| valor.is_finite());
            let valor = valor.ok_or_else(|| {
                format!(
                    "En la fila {numero_de_fila} de números (contando desde 0), «{}» no es un número.",
                    recortar(pieza)
                )
            })?;
            fila.push(valor);
        }
        filas.push(fila);
    }
    if filas.is_empty() {
        return Err("La matriz no tiene filas de números.".to_string());
    }
    let nombres = match nombres {
        Some(nombres) => {
            if nombres.len() != filas.len() {
                return Err(format!(
                    "La primera línea tiene {} nombres pero hay {} filas de números.",
                    nombres.len(),
                    filas.len()
                ));
            }
            revisar_nombres_unicos(&nombres)?;
            nombres
        }
        None => (0..filas.len()).map(|indice| indice.to_string()).collect(),
    };
    filogenia::validar_distancias(&filas)?;
    if filas.len() < 2 {
        return Err("Hacen falta al menos dos especies.".to_string());
    }
    Ok((nombres, filas))
}

fn secuencia_de_texto(valor: &str) -> Result<Vec<u8>, String> {
    if valor.starts_with('@') {
        return Entrada::nueva(&[("secuencia", valor)]).secuencia("secuencia");
    }
    Ok(valor
        .bytes()
        .filter(|letra| !letra.is_ascii_whitespace())
        .map(|letra| letra.to_ascii_uppercase())
        .collect())
}

fn leer_secuencias(
    entrada: &Entrada,
    campo: &str,
    admite_huecos: bool,
) -> Result<(Vec<String>, Vec<Vec<u8>>), String> {
    let texto = entrada.texto(campo)?;
    let mut nombres = Vec::new();
    let mut secuencias = Vec::new();
    if texto.trim_start().starts_with('>') {
        for registro in entrada.registros(campo)? {
            let nombre = registro
                .nombre
                .split_whitespace()
                .next()
                .unwrap_or_default();
            nombres.push(nombre.to_string());
            secuencias.push(registro.secuencia);
        }
    } else {
        for linea in texto
            .lines()
            .map(str::trim)
            .filter(|linea| !linea.is_empty())
        {
            let (nombre, valor) = linea.split_once('=').ok_or_else(|| {
                format!(
                    "Escribe cada secuencia en su línea como nombre=SECUENCIA (o usa FASTA): a «{}» le falta el «=».",
                    recortar(linea)
                )
            })?;
            nombres.push(nombre.trim().to_string());
            secuencias.push(secuencia_de_texto(valor.trim())?);
            if nombres.len() > MAXIMO_DE_ESPECIES {
                break;
            }
        }
    }
    if nombres.len() > MAXIMO_DE_ESPECIES {
        return Err(format!(
            "Hay más de {MAXIMO_DE_ESPECIES} secuencias; ese es el máximo aquí."
        ));
    }
    revisar_nombres_unicos(&nombres)?;
    for (nombre, secuencia) in nombres.iter().zip(&secuencias) {
        if secuencia.is_empty() {
            return Err(format!("«{nombre}» no tiene secuencia."));
        }
        if secuencia.len() > MAXIMO_DE_POSICIONES {
            return Err(format!(
                "La secuencia de «{nombre}» mide {} posiciones; el máximo aquí es {}.",
                miles(secuencia.len()),
                miles(MAXIMO_DE_POSICIONES)
            ));
        }
        let valida = |letra: &u8| {
            matches!(letra, b'A' | b'C' | b'G' | b'T') || (admite_huecos && *letra == b'-')
        };
        if let Some(letra) = secuencia.iter().find(|letra| !valida(letra)) {
            let permitidas = if admite_huecos {
                "A, C, G, T y «-» para los huecos"
            } else {
                "A, C, G y T"
            };
            return Err(format!(
                "La secuencia de «{nombre}» tiene «{}»: aquí solo valen {permitidas}.",
                *letra as char
            ));
        }
    }
    Ok((nombres, secuencias))
}

fn revisar_mismo_largo(nombres: &[String], secuencias: &[Vec<u8>]) -> Result<usize, String> {
    let largo = secuencias.first().map(Vec::len).unwrap_or(0);
    for (nombre, secuencia) in nombres.iter().zip(secuencias) {
        if secuencia.len() != largo {
            return Err(format!(
                "Las secuencias deben estar alineadas: «{}» mide {} y «{nombre}» mide {}.",
                nombres[0],
                miles(largo),
                miles(secuencia.len())
            ));
        }
    }
    Ok(largo)
}

fn buscar_especie(texto: &str, nombres: &[String]) -> Result<usize, String> {
    let texto = texto.trim();
    if let Some(posicion) = nombres.iter().position(|nombre| nombre == texto) {
        return Ok(posicion);
    }
    match texto.parse::<usize>() {
        Ok(indice) if indice < nombres.len() => Ok(indice),
        _ => Err(format!(
            "«{}» no es una de las especies de la matriz ni su número (contando desde 0).",
            recortar(texto)
        )),
    }
}

fn nombres_unicos(arbol: &mut Arbol, raiz: Option<usize>) {
    let mut usados: HashSet<String> = HashSet::new();
    for nodo in 0..arbol.nombres.len() {
        if !arbol.nombres[nodo].is_empty() && arbol.grado(nodo) <= 1 {
            usados.insert(arbol.nombres[nodo].clone());
        }
    }
    for nodo in 0..arbol.nombres.len() {
        let es_hoja = !arbol.nombres[nodo].is_empty() && arbol.grado(nodo) <= 1;
        if es_hoja {
            continue;
        }
        let mut nombre = arbol.nombres[nodo].clone();
        if nombre.is_empty() {
            nombre = if Some(nodo) == raiz {
                "raíz".to_string()
            } else {
                filogenia::nombre_interno(nodo)
            };
        }
        while usados.contains(&nombre) {
            nombre.push('_');
        }
        usados.insert(nombre.clone());
        arbol.nombres[nodo] = nombre;
    }
}

fn raiz_para_dibujar(arbol: &Arbol) -> usize {
    (0..arbol.nombres.len())
        .rev()
        .find(|&nodo| arbol.grado(nodo) >= 2)
        .unwrap_or(0)
}

fn newick_para_dibujar(arbol: &Arbol, raiz: usize) -> String {
    let mut sin_negativas = arbol.clone();
    for vecinos in sin_negativas.vecinos.iter_mut() {
        for (_, longitud) in vecinos.iter_mut() {
            *longitud = longitud.max(0.0);
        }
    }
    filogenia::escribir_newick(&sin_negativas, raiz)
}

fn lista_de_adyacencia(arbol: &Arbol) -> Vec<String> {
    let mut lineas = Vec::new();
    for (nodo, vecinos) in arbol.vecinos.iter().enumerate() {
        let mut ordenados = vecinos.clone();
        ordenados.sort_by_key(|&(vecino, _)| vecino);
        for (vecino, longitud) in ordenados {
            lineas.push(format!(
                "{}->{}:{}",
                arbol.nombres[nodo],
                arbol.nombres[vecino],
                numero(longitud)
            ));
        }
    }
    lineas
}

fn matriz_en_texto(nombres: &[String], matriz: &[Vec<f64>]) -> String {
    let mut texto = nombres.join(" ");
    for fila in matriz {
        texto.push('\n');
        let valores: Vec<String> = fila.iter().map(|&valor| numero(valor)).collect();
        texto.push_str(&valores.join(" "));
    }
    texto
}

fn agregar_matriz(salida: &mut Salida, etiqueta: &str, nombres: &[String], matriz: &[Vec<f64>]) {
    if nombres.len() <= MAXIMO_EN_TABLA {
        let mut columnas = vec![""];
        columnas.extend(nombres.iter().map(String::as_str));
        let filas = nombres
            .iter()
            .zip(matriz)
            .map(|(nombre, fila)| {
                let mut celdas = vec![nombre.clone()];
                celdas.extend(fila.iter().map(|&valor| numero(valor)));
                celdas
            })
            .collect();
        salida.agregar(Bloque::tabla(etiqueta, &columnas, filas));
    }
    salida.agregar(Bloque::calor(
        etiqueta,
        nombres.to_vec(),
        nombres.to_vec(),
        matriz_redondeada(matriz),
    ));
}

fn pares_extremos(matriz: &[Vec<f64>]) -> ((usize, usize), (usize, usize)) {
    let mut cercano = (0, 1);
    let mut lejano = (0, 1);
    for i in 0..matriz.len() {
        for j in (i + 1)..matriz.len() {
            if matriz[i][j] < matriz[cercano.0][cercano.1] {
                cercano = (i, j);
            }
            if matriz[i][j] > matriz[lejano.0][lejano.1] {
                lejano = (i, j);
            }
        }
    }
    (cercano, lejano)
}

fn mayor_diferencia(arbol: &Arbol, matriz: &[Vec<f64>]) -> Option<(usize, usize, f64)> {
    let hojas: Vec<usize> = (0..matriz.len()).collect();
    let del_arbol = filogenia::distancias_entre_hojas(arbol, &hojas);
    let escala = matriz
        .iter()
        .flatten()
        .fold(1.0_f64, |mayor, &valor| mayor.max(valor.abs()));
    let mut peor = None;
    let mut mayor = 1e-6 * escala;
    for i in 0..matriz.len() {
        for j in (i + 1)..matriz.len() {
            let diferencia = (del_arbol[i][j] - matriz[i][j]).abs();
            if diferencia > mayor {
                mayor = diferencia;
                peor = Some((i, j, del_arbol[i][j]));
            }
        }
    }
    peor
}

fn explicar_ramas_negativas(salida: &mut Salida, arbol: &Arbol) {
    let negativas: Vec<f64> = arbol
        .aristas()
        .iter()
        .map(|arista| arista.2)
        .filter(|&longitud| longitud < -1e-9)
        .collect();
    if let Some(menor) = negativas.iter().copied().reduce(f64::min) {
        salida.explicar(format!(
            "Hay {} rama(s) con longitud negativa (la menor mide {}): no tienen sentido biológico y aparecen cuando las distancias se alejan de un árbol. En el dibujo valen 0.",
            negativas.len(),
            numero(menor)
        ));
    }
}

fn explicar_rama_mas_larga(salida: &mut Salida, arbol: &Arbol, cantidad_de_hojas: usize) {
    let mas_larga = (0..cantidad_de_hojas)
        .filter_map(|hoja| arbol.vecinos[hoja].first().map(|&(_, largo)| (hoja, largo)))
        .max_by(|a, b| a.1.total_cmp(&b.1));
    if let Some((hoja, largo)) = mas_larga
        && largo > 1e-9
    {
        salida.explicar(format!(
            "La rama terminal más larga es la de «{}» ({}): es la especie que más cambió por su cuenta.",
            arbol.nombres[hoja],
            numero(largo)
        ));
    }
}

fn distancias_entre_hojas(entrada: &Entrada) -> Result<Salida, String> {
    let arbol = filogenia::leer_aristas(&entrada.texto("arbol")?)?;
    let hojas = arbol.hojas();
    if hojas.len() < 2 {
        return Err("El árbol necesita al menos dos hojas.".to_string());
    }
    if hojas.len() > MAXIMO_DE_ESPECIES {
        return Err(format!(
            "El árbol tiene {} hojas; el máximo aquí es {MAXIMO_DE_ESPECIES}.",
            miles(hojas.len())
        ));
    }
    let matriz = filogenia::distancias_entre_hojas(&arbol, &hojas);
    let nombres: Vec<String> = hojas
        .iter()
        .map(|&hoja| arbol.nombres[hoja].clone())
        .collect();
    let internos = arbol.nombres.len() - hojas.len();
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Hojas", miles(hojas.len())));
    salida.agregar(Bloque::valor("Nodos internos", miles(internos)));
    agregar_matriz(&mut salida, "Distancias entre hojas", &nombres, &matriz);
    salida.agregar(Bloque::texto(
        "Matriz para copiar",
        matriz_en_texto(&nombres, &matriz),
    ));
    salida.agregar(Bloque::arbol(
        "Árbol",
        newick_para_dibujar(&arbol, raiz_para_dibujar(&arbol)),
    ));
    let ((a, b), (c, d)) = pares_extremos(&matriz);
    salida.explicar(format!(
        "La distancia entre dos hojas es la suma de las longitudes del único camino que las une. Las más alejadas son «{}» y «{}» ({}); las más cercanas, «{}» y «{}» ({}).",
        nombres[c],
        nombres[d],
        numero(matriz[c][d]),
        nombres[a],
        nombres[b],
        numero(matriz[a][b])
    ));
    if casi_cero(matriz[a][b]) {
        salida.explicar(format!(
            "«{}» y «{}» están a distancia 0: el árbol no las distingue.",
            nombres[a], nombres[b]
        ));
    }
    let de_paso = (0..arbol.nombres.len())
        .filter(|&nodo| arbol.grado(nodo) == 2)
        .count();
    let ramas_nulas = arbol
        .aristas()
        .iter()
        .filter(|arista| casi_cero(arista.2))
        .count();
    if de_paso > 0 {
        salida.explicar(format!(
            "Hay {de_paso} nodo(s) interno(s) con solo dos vecinos: no cambian ninguna distancia, así que «Filogenia aditiva» los omitirá al reconstruir el árbol desde esta matriz."
        ));
    } else if ramas_nulas > 0 {
        salida.explicar(format!(
            "Hay {ramas_nulas} rama(s) de longitud 0: al reconstruir desde la matriz, sus dos extremos se funden en un solo nodo."
        ));
    } else {
        salida.explicar(
            "Toda matriz que sale de un árbol es aditiva: copia la matriz en «Filogenia aditiva» y recuperarás este mismo árbol.",
        );
    }
    Ok(salida)
}

fn longitud_de_rama(entrada: &Entrada) -> Result<Salida, String> {
    let (nombres, matriz) = leer_matriz(entrada, "matriz")?;
    let hoja = buscar_especie(&entrada.texto("hoja")?, &nombres)?;
    let (rama, primera, segunda) = filogenia::longitud_de_rama_con_par(&matriz, hoja)?;
    let mut categorias = Vec::new();
    let mut valores = Vec::new();
    for otra in 0..nombres.len() {
        if otra == hoja || otra == primera {
            continue;
        }
        categorias.push(nombres[otra].clone());
        valores.push(redondear(
            (matriz[primera][hoja] + matriz[hoja][otra] - matriz[primera][otra]) / 2.0,
        ));
    }
    let cercana = (0..nombres.len())
        .filter(|&otra| otra != hoja)
        .min_by(|&a, &b| matriz[hoja][a].total_cmp(&matriz[hoja][b]))
        .unwrap_or(primera);
    let nombre = &nombres[hoja];
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        &format!("Longitud de la rama de «{nombre}»"),
        numero(rama),
    ));
    salida.agregar(Bloque::valor(
        "Pariente más cercano",
        format!(
            "{} (distancia {})",
            nombres[cercana],
            numero(matriz[hoja][cercana])
        ),
    ));
    if !valores.is_empty() {
        salida.agregar(Bloque::barras(
            &format!(
                "Candidatas (D({0},{nombre}) + D({nombre},k) − D({0},k)) / 2",
                nombres[primera]
            ),
            "Otra hoja k",
            "Candidata",
            categorias,
            valores,
        ));
    }
    let significado = if rama < -1e-9 {
        String::new()
    } else {
        format!(": es lo que «{nombre}» cambió por su cuenta desde que se separó del resto")
    };
    if primera == segunda {
        salida.explicar(format!(
            "La rama de «{nombre}» mide {}{significado}. Es el menor valor de (D(i, «{nombre}») + D(«{nombre}», k) − D(i, k)) / 2 entre todas las parejas de otras hojas i y k, y se alcanza con i = k = «{}», que da simplemente la distancia entre ambas.",
            numero(rama),
            nombres[primera]
        ));
    } else {
        salida.explicar(format!(
            "La rama de «{nombre}» mide {}{significado}. Es el menor valor de (D(i, «{nombre}») + D(«{nombre}», k) − D(i, k)) / 2 entre todas las parejas de otras hojas i y k, y se alcanza con i = «{}» y k = «{}».",
            numero(rama),
            nombres[primera],
            nombres[segunda]
        ));
    }
    if !filogenia::es_aditiva(&matriz) {
        salida.explicar(
            "La matriz no es aditiva (ningún árbol reproduce exactamente sus distancias): no todas las parejas de hojas dan el mismo valor, y se toma el menor, que es solo una aproximación.",
        );
    }
    if nombres.len() == 2 {
        salida.explicar(format!(
            "Con solo dos hojas no hay nada con qué comparar: la rama de «{nombre}» es toda la distancia entre ambas."
        ));
    } else if rama < -1e-9 {
        salida.explicar(
            "Salió negativa, algo imposible en un árbol real: confirma que estas distancias no vienen de un árbol.",
        );
    } else if casi_cero(rama) {
        salida.explicar(format!(
            "Mide 0: «{nombre}» cae justo sobre un nodo interno, como si fuera el antepasado de otras especies."
        ));
    } else if rama > matriz[hoja][cercana] / 2.0 {
        salida.explicar(format!(
            "Es más de la mitad de la distancia a su pariente más cercano («{}», {}): «{nombre}» tiene una rama larga, acumuló muchos cambios propios.",
            nombres[cercana],
            numero(matriz[hoja][cercana])
        ));
    }
    Ok(salida)
}

fn filogenia_aditiva(entrada: &Entrada) -> Result<Salida, String> {
    let (nombres, matriz) = leer_matriz(entrada, "matriz")?;
    match filogenia::filogenia_aditiva(&matriz, &nombres) {
        Ok(arbol) => Ok(informe_aditivo(arbol, &nombres)),
        Err(error) => explicar_no_aditiva(&nombres, &matriz).ok_or(error),
    }
}

fn explicar_no_aditiva(nombres: &[String], matriz: &[Vec<f64>]) -> Option<Salida> {
    let mut salida = Salida::default();
    if let Some((i, j, k)) = filogenia::trio_no_metrico(matriz) {
        salida.agregar(Bloque::valor("¿Es aditiva?", "No"));
        let mut lados = [(i, j, k), (i, k, j), (j, k, i)];
        lados.sort_by(|a, b| matriz[a.0][a.1].total_cmp(&matriz[b.0][b.1]));
        let (a, b, c) = lados[2];
        salida.agregar(Bloque::tabla(
            "Trío que falla",
            &["Distancia", "Valor"],
            vec![
                vec![
                    format!("d({}, {})", nombres[a], nombres[b]),
                    numero(matriz[a][b]),
                ],
                vec![
                    format!(
                        "d({}, {}) + d({}, {})",
                        nombres[a], nombres[c], nombres[c], nombres[b]
                    ),
                    numero(matriz[a][c] + matriz[c][b]),
                ],
            ],
        ));
        salida.explicar(format!(
            "No es aditiva: d(«{}», «{}») = {} es mayor que pasar por «{}» ({}). Eso rompe la desigualdad del triángulo, así que ningún árbol puede dar estas distancias.",
            nombres[a],
            nombres[b],
            numero(matriz[a][b]),
            nombres[c],
            numero(matriz[a][c] + matriz[c][b])
        ));
        salida.explicar("Revisa los datos o prueba «Unión de vecinos», que busca el árbol que mejor se aproxima.");
        return Some(salida);
    }
    let (i, j, k, l) = filogenia::cuarteto_no_aditivo(matriz)?;
    let (primera, segunda, tercera) = filogenia::sumas_de_cuarteto(matriz, i, j, k, l);
    let sumas = [primera, segunda, tercera];
    let etiquetas = [
        format!(
            "d({}, {}) + d({}, {})",
            nombres[i], nombres[j], nombres[k], nombres[l]
        ),
        format!(
            "d({}, {}) + d({}, {})",
            nombres[i], nombres[k], nombres[j], nombres[l]
        ),
        format!(
            "d({}, {}) + d({}, {})",
            nombres[i], nombres[l], nombres[j], nombres[k]
        ),
    ];
    salida.agregar(Bloque::valor("¿Es aditiva?", "No"));
    salida.agregar(Bloque::tabla(
        "Cuarteto que falla",
        &["Suma", "Valor"],
        etiquetas
            .iter()
            .zip(sumas)
            .map(|(etiqueta, suma)| vec![etiqueta.clone(), numero(suma)])
            .collect(),
    ));
    salida.explicar(format!(
            "No es aditiva: con «{}», «{}», «{}» y «{}» las tres sumas de pares valen {}, {} y {}, y las dos mayores no coinciden. En un árbol siempre coinciden (condición de los cuatro puntos), así que ningún árbol reproduce exactamente estas distancias.",
            nombres[i],
            nombres[j],
            nombres[k],
            nombres[l],
            numero(sumas[0]),
            numero(sumas[1]),
            numero(sumas[2])
        ));
    salida.explicar(
        "Prueba «Unión de vecinos»: no necesita aditividad y busca el árbol que mejor se aproxima.",
    );
    Some(salida)
}

fn informe_aditivo(mut arbol: Arbol, nombres: &[String]) -> Salida {
    nombres_unicos(&mut arbol, None);
    let cantidad = nombres.len();
    let internos = arbol.nombres.len() - cantidad;
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("¿Es aditiva?", "Sí"));
    salida.agregar(Bloque::valor("Nodos internos", miles(internos)));
    salida.agregar(Bloque::arbol(
        "Árbol aditivo",
        newick_para_dibujar(&arbol, raiz_para_dibujar(&arbol)),
    ));
    salida.agregar(Bloque::lista(
        "Lista de adyacencia",
        lista_de_adyacencia(&arbol),
    ));
    salida.explicar(format!(
        "Cumple la condición de los cuatro puntos: en cada grupo de cuatro especies, de las tres formas de sumar distancias por pares, las dos mayores son iguales. Por eso existe un único árbol simple que reproduce exactamente {}.",
        if cantidad == 2 {
            "la distancia entre las dos especies".to_string()
        } else {
            format!("las {} distancias", miles(cantidad * (cantidad - 1) / 2))
        }
    ));
    if cantidad >= 3 && internos == cantidad - 2 {
        salida.explicar("Es binario: cada nodo interno une exactamente tres ramas.");
    } else if cantidad >= 3 {
        salida.explicar(format!(
            "Tiene {internos} nodo(s) interno(s) en vez de {}: algún nodo une más de tres ramas, o alguna especie quedó dentro del árbol.",
            cantidad - 2
        ));
    }
    let dentro: Vec<&str> = (0..cantidad)
        .filter(|&hoja| arbol.grado(hoja) >= 2)
        .map(|hoja| nombres[hoja].as_str())
        .collect();
    if !dentro.is_empty() {
        salida.explicar(format!(
            "«{}» quedó como nodo interno: su rama mide 0, así que en este árbol es el antepasado de otras especies.",
            dentro.join("», «")
        ));
    }
    salida.explicar(
        "Con datos reales la aditividad exacta es rara; cuando falla, usa «Unión de vecinos».",
    );
    salida
}

fn alturas(arbol: &Arbol, raiz: usize) -> Vec<f64> {
    let desde_la_raiz = filogenia::distancias_desde(arbol, raiz);
    let altura_total = desde_la_raiz
        .iter()
        .zip(&arbol.vecinos)
        .filter(|(_, vecinos)| vecinos.len() <= 1)
        .map(|(distancia, _)| *distancia)
        .fold(0.0, f64::max);
    desde_la_raiz
        .iter()
        .map(|distancia| altura_total - distancia)
        .collect()
}

fn upgma(entrada: &Entrada) -> Result<Salida, String> {
    let (nombres, matriz) = leer_matriz(entrada, "matriz")?;
    let resultado = filogenia::upgma(&matriz, &nombres)?;
    let raiz = resultado.raiz;
    let mut arbol = resultado.arbol;
    nombres_unicos(&mut arbol, None);
    let cantidad = nombres.len();
    let altura = alturas(&arbol, raiz);
    let desde_la_raiz = filogenia::distancias_desde(&arbol, raiz);
    let mut filas = Vec::new();
    for nodo in cantidad..arbol.nombres.len() {
        let hijos: Vec<&str> = arbol.vecinos[nodo]
            .iter()
            .filter(|&&(vecino, _)| desde_la_raiz[vecino] > desde_la_raiz[nodo])
            .map(|&(vecino, _)| arbol.nombres[vecino].as_str())
            .collect();
        filas.push(vec![
            arbol.nombres[nodo].clone(),
            hijos.join(" + "),
            numero(altura[nodo]),
        ]);
    }
    let mut salida = Salida::default();
    salida.agregar(Bloque::arbol(
        "Árbol UPGMA (con raíz)",
        newick_para_dibujar(&arbol, raiz),
    ));
    salida.agregar(Bloque::valor("Altura de la raíz", numero(altura[raiz])));
    salida.agregar(Bloque::tabla(
        "Uniones en orden",
        &["Nodo nuevo", "Une a", "Altura"],
        filas.clone(),
    ));
    salida.agregar(Bloque::lista(
        "Lista de adyacencia",
        lista_de_adyacencia(&arbol),
    ));
    salida.explicar(format!(
        "UPGMA supone un reloj molecular: todas las especies cambian al mismo ritmo, por eso todas las hojas quedan a la misma altura ({}) bajo la raíz.",
        numero(altura[raiz])
    ));
    if let Some(primera) = filas.first() {
        salida.explicar(format!(
            "La primera unión fue {} a altura {}: es el par más parecido de la matriz.",
            primera[1], primera[2]
        ));
    }
    match filogenia::trio_no_ultrametrico(&matriz) {
        None => salida.explicar(
            "La matriz es ultramétrica (en cada trío las dos distancias mayores son iguales): el reloj molecular se cumple y el árbol reproduce exactamente las distancias.",
        ),
        Some((i, j, k)) => salida.explicar(format!(
            "La matriz no es ultramétrica: en el trío «{}», «{}», «{}» las distancias son {}, {} y {}, y las dos mayores no son iguales. Las especies no cambiaron al mismo ritmo y UPGMA puede agrupar mal; compara con «Unión de vecinos», que no supone reloj.",
            nombres[i],
            nombres[j],
            nombres[k],
            numero(matriz[i][j]),
            numero(matriz[i][k]),
            numero(matriz[j][k])
        )),
    }
    if let Some((i, j, en_el_arbol)) = mayor_diferencia(&arbol, &matriz) {
        salida.explicar(format!(
            "La mayor diferencia entre el árbol y la matriz está entre «{}» y «{}»: {} en la matriz contra {} en el árbol.",
            nombres[i],
            nombres[j],
            numero(matriz[i][j]),
            numero(en_el_arbol)
        ));
    }
    Ok(salida)
}

fn union_de_vecinos(entrada: &Entrada) -> Result<Salida, String> {
    let (nombres, matriz) = leer_matriz(entrada, "matriz")?;
    let mut estado = filogenia::empezar_union(&matriz, &nombres)?;
    let mut pares = Vec::new();
    while let Some(par) = filogenia::unir_par_de_vecinos(&mut estado) {
        pares.push(par);
    }
    filogenia::unir_ultimo_par(&mut estado);
    let mut arbol = estado.arbol;
    nombres_unicos(&mut arbol, None);
    let filas = pares
        .iter()
        .map(|par| {
            vec![
                arbol.nombres[par.nuevo].clone(),
                format!(
                    "{} + {}",
                    arbol.nombres[par.primero], arbol.nombres[par.segundo]
                ),
                numero(par.rama_primero),
                numero(par.rama_segundo),
            ]
        })
        .collect();
    let cantidad = nombres.len();
    let mut salida = Salida::default();
    salida.agregar(Bloque::arbol(
        "Árbol de unión de vecinos (sin raíz)",
        newick_para_dibujar(&arbol, raiz_para_dibujar(&arbol)),
    ));
    salida.agregar(Bloque::tabla(
        "Uniones en orden",
        &[
            "Nodo nuevo",
            "Une a",
            "Rama del primero",
            "Rama del segundo",
        ],
        filas,
    ));
    salida.agregar(Bloque::lista(
        "Lista de adyacencia",
        lista_de_adyacencia(&arbol),
    ));
    salida.explicar(
        "Unión de vecinos no supone reloj molecular: cada rama puede medir distinto, así que el árbol no tiene raíz.",
    );
    if filogenia::es_aditiva(&matriz) {
        salida.explicar(
            "La matriz es aditiva: en ese caso unión de vecinos reconstruye exactamente el árbol que la produjo (es un teorema).",
        );
    } else {
        let falla = match filogenia::cuarteto_no_aditivo(&matriz) {
            Some((i, j, k, l)) => format!(
                " (falla la condición de los cuatro puntos con «{}», «{}», «{}» y «{}»)",
                nombres[i], nombres[j], nombres[k], nombres[l]
            ),
            None => String::new(),
        };
        let diferencia = match mayor_diferencia(&arbol, &matriz) {
            Some((i, j, en_el_arbol)) => format!(
                " La mayor diferencia está entre «{}» y «{}»: {} en la matriz contra {} en el árbol.",
                nombres[i],
                nombres[j],
                numero(matriz[i][j]),
                numero(en_el_arbol)
            ),
            None => String::new(),
        };
        salida.explicar(format!(
            "La matriz no es aditiva{falla}, así que el árbol es una aproximación.{diferencia}"
        ));
    }
    explicar_ramas_negativas(&mut salida, &arbol);
    if filogenia::trio_no_ultrametrico(&matriz).is_some() {
        salida.explicar(
            "La matriz no es ultramétrica: UPGMA podría equivocarse con estos datos, y unión de vecinos es la opción adecuada.",
        );
    }
    explicar_rama_mas_larga(&mut salida, &arbol, cantidad);
    Ok(salida)
}

fn matriz_redondeada(matriz: &[Vec<f64>]) -> Vec<Vec<f64>> {
    matriz
        .iter()
        .map(|fila| fila.iter().map(|&valor| redondear(valor)).collect())
        .collect()
}

fn aristas_en_json(arbol: &Arbol, etiquetas: &[String]) -> Vec<Value> {
    let mut aristas = Vec::new();
    for (nodo, vecinos) in arbol.vecinos.iter().enumerate() {
        for &(vecino, longitud) in vecinos {
            if vecino < nodo {
                aristas.push(json!({
                    "desde": etiquetas[nodo],
                    "hasta": etiquetas[vecino],
                    "longitud": redondear(longitud),
                }));
            }
        }
    }
    aristas
}

fn foto_de_union(
    fase: &str,
    ronda: usize,
    estado: &filogenia::EstadoDeUnion,
    etiquetas: &[String],
    extra: Value,
) -> Value {
    let nombres: Vec<&str> = estado
        .nodos
        .iter()
        .map(|&nodo| etiquetas[nodo].as_str())
        .collect();
    let mut foto = json!({
        "fase": fase,
        "ronda": ronda,
        "nombres": nombres,
        "distancias": matriz_redondeada(&estado.distancias),
        "totales": null,
        "matriz_de_union": null,
        "par": null,
        "minimo": null,
        "nuevo": null,
        "delta": null,
        "ramas": null,
        "aristas": aristas_en_json(&estado.arbol, etiquetas),
        "newick": null,
    });
    if let (Value::Object(campos), Value::Object(agregados)) = (&mut foto, extra) {
        for (clave, valor) in agregados {
            campos.insert(clave, valor);
        }
    }
    foto
}

fn pasos_union_de_vecinos(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (nombres, matriz) = leer_matriz(entrada, "matriz")?;
    if nombres.len() > MAXIMO_PARA_ANIMAR {
        return Err(format!(
            "El paso a paso admite hasta {MAXIMO_PARA_ANIMAR} especies (la matriz tiene {}); usa «Ejecutar» para matrices más grandes.",
            nombres.len()
        ));
    }
    let cantidad = nombres.len();
    let mut completo = filogenia::union_de_vecinos(&matriz, &nombres)?;
    nombres_unicos(&mut completo, None);
    let etiquetas = completo.nombres.clone();
    let mut estado = filogenia::empezar_union(&matriz, &nombres)?;
    let mut pasos = vec![Paso::nuevo(
        "let mut estado = empezar_union(distancias, nombres)?;",
        format!("Empezamos con las distancias entre las {cantidad} especies."),
        foto_de_union("inicio", 0, &estado, &etiquetas, json!({})),
    )];
    let mut ronda = 0;
    while estado.distancias.len() > 2 {
        ronda += 1;
        let antes = estado.clone();
        let totales = filogenia::distancias_totales(&antes.distancias);
        let matriz_estrella = filogenia::matriz_de_union(&antes.distancias);
        let Some(par) = filogenia::unir_par_de_vecinos(&mut estado) else {
            break;
        };
        let i = antes
            .nodos
            .iter()
            .position(|&nodo| nodo == par.primero)
            .unwrap_or(0);
        let j = antes
            .nodos
            .iter()
            .position(|&nodo| nodo == par.segundo)
            .unwrap_or(0);
        let primero = &etiquetas[par.primero];
        let segundo = &etiquetas[par.segundo];
        let nuevo = &etiquetas[par.nuevo];
        let n = antes.distancias.len();
        let delta = (totales[i] - totales[j]) / (n - 2) as f64;
        let totales_redondeados: Vec<f64> = totales.iter().map(|&total| redondear(total)).collect();
        pasos.push(Paso::nuevo(
            "let totales = distancias_totales(&estado.distancias);",
            format!("Se suma cada fila: la distancia total de cada uno de los {n} nodos a todos los demás."),
            foto_de_union(
                "totales",
                ronda,
                &antes,
                &etiquetas,
                json!({ "totales": totales_redondeados }),
            ),
        ));
        pasos.push(Paso::nuevo(
            "let (i, j) = par_con_menor_valor_de_union(&estado.distancias, &totales);",
            format!(
                "D*(i,j) = ({n}−2)·D(i,j) − total(i) − total(j) es mínimo ({}) para «{primero}» y «{segundo}»: son vecinos.",
                numero(matriz_estrella[i][j])
            ),
            foto_de_union(
                "eleccion",
                ronda,
                &antes,
                &etiquetas,
                json!({
                    "totales": totales_redondeados,
                    "matriz_de_union": matriz_redondeada(&matriz_estrella),
                    "par": [i, j],
                    "minimo": redondear(matriz_estrella[i][j]),
                }),
            ),
        ));
        let mut con_ramas = antes.clone();
        con_ramas.arbol = estado.arbol.clone();
        pasos.push(Paso::nuevo(
            "let rama_i = (estado.distancias[i][j] + delta) / 2.0;",
            format!(
                "Con Δ = (total(i) − total(j)) / ({n}−2) = {}, el nodo nuevo «{nuevo}» cuelga a «{primero}» con una rama de (D(i,j) + Δ) / 2 = {} y a «{segundo}» con una de (D(i,j) − Δ) / 2 = {}.",
                numero(delta),
                numero(par.rama_primero),
                numero(par.rama_segundo)
            ),
            foto_de_union(
                "ramas",
                ronda,
                &con_ramas,
                &etiquetas,
                json!({
                    "totales": totales_redondeados,
                    "par": [i, j],
                    "nuevo": nuevo,
                    "delta": redondear(delta),
                    "ramas": [redondear(par.rama_primero), redondear(par.rama_segundo)],
                }),
            ),
        ));
        pasos.push(Paso::nuevo(
            "reemplazar_par(&mut estado.distancias, i, j, fila_nueva);",
            format!(
                "«{primero}» y «{segundo}» salen de la matriz y entra «{nuevo}», a (D(k,i) + D(k,j) − D(i,j)) / 2 de cada nodo k."
            ),
            foto_de_union(
                "actualizacion",
                ronda,
                &estado,
                &etiquetas,
                json!({ "nuevo": nuevo }),
            ),
        ));
    }
    let antes_del_final = estado.clone();
    if let Some((a, b, longitud)) = filogenia::unir_ultimo_par(&mut estado) {
        let mut cierre = antes_del_final;
        cierre.arbol = estado.arbol.clone();
        let newick = newick_para_dibujar(&completo, raiz_para_dibujar(&completo));
        pasos.push(Paso::nuevo(
            ".conectar(estado.nodos[0], estado.nodos[1], longitud);",
            format!(
                "Quedan dos nodos, «{}» y «{}»: se unen con una rama de {} y el árbol está completo.",
                etiquetas[a],
                etiquetas[b],
                numero(longitud)
            ),
            foto_de_union(
                "final",
                ronda + 1,
                &cierre,
                &etiquetas,
                json!({ "par": [0, 1], "newick": newick }),
            ),
        ));
    }
    Ok(pasos)
}

fn distancias_evolutivas(entrada: &Entrada) -> Result<Salida, String> {
    let (nombres, secuencias) = leer_secuencias(entrada, "secuencias", true)?;
    if nombres.len() < 2 {
        return Err("Hacen falta al menos dos secuencias.".to_string());
    }
    let largo = revisar_mismo_largo(&nombres, &secuencias)?;
    let distancias_p = filogenia::matriz_de_distancias_p(&secuencias)?;
    let jukes_cantor = filogenia::corregir_con_jukes_cantor(&distancias_p);
    let cantidad = nombres.len();
    let mut filas = Vec::new();
    let mut puntos = Vec::new();
    let mut saturados = Vec::new();
    let mut mas_distinto = (0, 1);
    for i in 0..cantidad {
        for j in (i + 1)..cantidad {
            let (comparadas, diferentes) =
                filogenia::contar_diferencias(&secuencias[i], &secuencias[j]);
            filas.push(vec![
                nombres[i].clone(),
                nombres[j].clone(),
                miles(comparadas),
                miles(diferentes),
                numero(distancias_p[i][j]),
                numero(jukes_cantor[i][j]),
            ]);
            if jukes_cantor[i][j].is_finite() {
                puntos.push(Punto {
                    x: redondear(distancias_p[i][j]),
                    y: redondear(jukes_cantor[i][j]),
                    grupo: 0,
                    nombre: format!("{}–{}", nombres[i], nombres[j]),
                });
                if distancias_p[i][j] > distancias_p[mas_distinto.0][mas_distinto.1]
                    || jukes_cantor[mas_distinto.0][mas_distinto.1].is_infinite()
                {
                    mas_distinto = (i, j);
                }
            } else {
                saturados.push((i, j));
            }
        }
    }
    let con_hueco = (0..largo)
        .filter(|&posicion| {
            secuencias
                .iter()
                .any(|secuencia| secuencia[posicion] == b'-')
        })
        .count();
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Secuencias", miles(cantidad)));
    salida.agregar(Bloque::valor("Posiciones alineadas", miles(largo)));
    agregar_matriz(&mut salida, "Distancia p", &nombres, &distancias_p);
    salida.agregar(Bloque::tabla(
        "Pares",
        &[
            "Una",
            "Otra",
            "Comparadas",
            "Distintas",
            "p",
            "Jukes–Cantor",
        ],
        filas,
    ));
    salida.agregar(Bloque::dispersion(
        "Corrección de Jukes–Cantor",
        "p (fracción de posiciones distintas)",
        "Jukes–Cantor (sustituciones por posición)",
        puntos,
    ));
    if saturados.is_empty() {
        agregar_matriz(
            &mut salida,
            "Distancia de Jukes–Cantor",
            &nombres,
            &jukes_cantor,
        );
        salida.agregar(Bloque::texto(
            "Matriz de Jukes–Cantor para copiar",
            matriz_en_texto(&nombres, &jukes_cantor),
        ));
    }
    if let Some(&(i, j)) = saturados.first() {
        salida.explicar(format!(
            "{} par(es) difieren en el 75 % o más de las posiciones comparadas (por ejemplo «{}» y «{}»): ahí Jukes–Cantor da infinito, porque dos secuencias sin ningún parentesco ya coinciden en 1 de cada 4 posiciones por azar. Están demasiado divergentes para medir su distancia así.",
            saturados.len(),
            nombres[i],
            nombres[j]
        ));
    }
    let (i, j) = mas_distinto;
    let p = distancias_p[i][j];
    let corregida = jukes_cantor[i][j];
    if corregida.is_finite() && p > 0.0 {
        salida.explicar(format!(
            "El par más distinto con distancia finita es «{}»–«{}»: p = {} ({} de posiciones distintas) y Jukes–Cantor = {} sustituciones por posición. La corrección suma {} porque estima los cambios ocultos: dos cambios en el mismo sitio se ven como uno o como ninguno.",
            nombres[i],
            nombres[j],
            numero(p),
            porcentaje(p, 1),
            numero(corregida),
            porcentaje(corregida / p - 1.0, 1)
        ));
    }
    let mayor_p = distancias_p
        .iter()
        .flatten()
        .fold(0.0_f64, |mayor, &valor| mayor.max(valor));
    if mayor_p == 0.0 {
        salida.explicar(
            "Todas las secuencias son idénticas en las posiciones comparadas: todas las distancias valen 0 y no hay señal para ordenarlas en un árbol.",
        );
    } else if mayor_p < 0.05 {
        salida.explicar(
            "Todas las diferencias están por debajo del 5 %: p y Jukes–Cantor casi coinciden (la corrección suma menos del 4 %).",
        );
    }
    if con_hueco > 0 {
        salida.explicar(format!(
            "Hay {} posición(es) con hueco («-») en alguna secuencia: cada par se compara solo donde ninguna de las dos tiene hueco.",
            miles(con_hueco)
        ));
    }
    if saturados.is_empty() && mayor_p > 0.0 {
        salida.explicar(
            "Copia la matriz de Jukes–Cantor en «Unión de vecinos», o usa «Árbol desde secuencias» para hacerlo de una vez.",
        );
    }
    Ok(salida)
}

fn arbol_desde_secuencias(entrada: &Entrada) -> Result<Salida, String> {
    let (nombres, secuencias) = leer_secuencias(entrada, "secuencias", true)?;
    if nombres.len() < 2 {
        return Err("Hacen falta al menos dos secuencias.".to_string());
    }
    revisar_mismo_largo(&nombres, &secuencias)?;
    let jukes_cantor = filogenia::matriz_jukes_cantor(&secuencias)?;
    let cantidad = nombres.len();
    let mut salida = Salida::default();
    if let Some((i, j)) = filogenia::par_saturado(&jukes_cantor) {
        let mut filas = Vec::new();
        for a in 0..cantidad {
            for b in (a + 1)..cantidad {
                if jukes_cantor[a][b].is_infinite() {
                    let p = filogenia::distancia_p(&secuencias[a], &secuencias[b]).unwrap_or(1.0);
                    filas.push(vec![nombres[a].clone(), nombres[b].clone(), numero(p)]);
                }
            }
        }
        salida.agregar(Bloque::valor("Árbol", "No se puede construir"));
        salida.agregar(Bloque::tabla(
            "Pares saturados",
            &["Una", "Otra", "p"],
            filas,
        ));
        salida.explicar(format!(
            "«{}» y «{}» difieren en el 75 % o más de sus posiciones: su distancia de Jukes–Cantor es infinita, porque a ese nivel ya no se distingue el parentesco del puro azar. Las secuencias están demasiado divergentes para este método; quita las más lejanas o usa un gen más conservado.",
            nombres[i], nombres[j]
        ));
        return Ok(salida);
    }
    let mut arbol = filogenia::arbol_desde_secuencias(&secuencias, &nombres)?;
    nombres_unicos(&mut arbol, None);
    salida.agregar(Bloque::arbol(
        "Árbol de unión de vecinos con distancias de Jukes–Cantor",
        newick_para_dibujar(&arbol, raiz_para_dibujar(&arbol)),
    ));
    salida.agregar(Bloque::lista(
        "Lista de adyacencia",
        lista_de_adyacencia(&arbol),
    ));
    agregar_matriz(
        &mut salida,
        "Distancias de Jukes–Cantor",
        &nombres,
        &jukes_cantor,
    );
    if cantidad <= MAXIMO_EN_ALINEAMIENTO {
        salida.agregar(Bloque::alineamiento(
            "Secuencias",
            nombres
                .iter()
                .zip(&secuencias)
                .map(|(nombre, secuencia)| FilaAlineada {
                    nombre: nombre.clone(),
                    secuencia: String::from_utf8_lossy(secuencia).to_string(),
                })
                .collect(),
        ));
    }
    salida.explicar(
        "Se contaron las diferencias de cada par (p), se corrigieron con Jukes–Cantor y con esas distancias se armó el árbol por unión de vecinos: las longitudes son sustituciones por posición.",
    );
    let ((a, b), (c, d)) = pares_extremos(&jukes_cantor);
    if casi_cero(jukes_cantor[c][d]) {
        salida.explicar(
            "Todas las secuencias son idénticas en las posiciones comparadas: todas las ramas miden 0 y el árbol no dice nada sobre el parentesco.",
        );
    } else {
        salida.explicar(format!(
            "Las más parecidas son «{}» y «{}» ({} sustituciones por posición).",
            nombres[a],
            nombres[b],
            numero(jukes_cantor[a][b])
        ));
    }
    if !filogenia::es_aditiva(&jukes_cantor) {
        salida.explicar(
            "Estas distancias no son aditivas, como casi siempre con datos reales: por eso se usa unión de vecinos y no «Filogenia aditiva».",
        );
    }
    explicar_ramas_negativas(&mut salida, &arbol);
    explicar_rama_mas_larga(&mut salida, &arbol, cantidad);
    Ok(salida)
}

fn secuencias_en_el_arbol(
    arbol: &Arbol,
    nombres: &[String],
    secuencias: Vec<Vec<u8>>,
) -> Result<Vec<Vec<u8>>, String> {
    let hojas = arbol.hojas().len();
    if hojas > MAXIMO_DE_ESPECIES {
        return Err(format!(
            "El árbol tiene {} hojas; el máximo aquí es {MAXIMO_DE_ESPECIES}.",
            miles(hojas)
        ));
    }
    if nombres.len() < 2 {
        return Err("Hacen falta al menos dos hojas con secuencia.".to_string());
    }
    let mut por_nodo = vec![Vec::new(); arbol.nombres.len()];
    for (nombre, secuencia) in nombres.iter().zip(secuencias) {
        let nodo = arbol
            .buscar_nodo(nombre)
            .ok_or_else(|| format!("«{nombre}» tiene secuencia pero no aparece en el árbol."))?;
        if arbol.grado(nodo) > 1 {
            return Err(format!(
                "«{nombre}» es un nodo interno del árbol: solo las hojas llevan secuencia."
            ));
        }
        por_nodo[nodo] = secuencia;
    }
    Ok(por_nodo)
}

fn minimo_de_cambios(secuencias_de_hojas: &[&Vec<u8>], largo: usize) -> usize {
    let mut minimo = 0;
    for posicion in 0..largo {
        let mut presentes = [false; 4];
        for secuencia in secuencias_de_hojas {
            if let Some(indice) = nucleo::secuencia::indice_base(secuencia[posicion]) {
                presentes[indice] = true;
            }
        }
        let distintas = presentes.iter().filter(|&&presente| presente).count();
        minimo += distintas.saturating_sub(1);
    }
    minimo
}

fn informe_de_parsimonia(
    salida: &mut Salida,
    resultado: &Parsimonia,
    arbol: &Arbol,
    por_nodo: &[Vec<u8>],
    raiz: usize,
) {
    let largo = resultado.cambios_por_posicion.len();
    let hojas: Vec<&Vec<u8>> = por_nodo
        .iter()
        .filter(|secuencia| !secuencia.is_empty())
        .collect();
    let minimo = minimo_de_cambios(&hojas, largo);
    let puntaje = resultado.puntaje;
    salida.agregar(Bloque::valor("Puntaje de parsimonia", mutaciones(puntaje)));
    salida.agregar(Bloque::arbol(
        "Árbol con el número de cambios en cada rama",
        filogenia::escribir_newick(arbol, raiz),
    ));
    salida.agregar(Bloque::lista(
        "Lista de adyacencia (cambios en cada rama)",
        lista_de_adyacencia(arbol),
    ));
    if arbol.nombres.len() <= MAXIMO_EN_ALINEAMIENTO {
        let (orden, _) = orden_desde(arbol, raiz);
        salida.agregar(Bloque::alineamiento(
            "Secuencias de hojas y antepasados",
            orden
                .iter()
                .map(|&nodo| FilaAlineada {
                    nombre: arbol.nombres[nodo].clone(),
                    secuencia: String::from_utf8_lossy(&resultado.secuencias[nodo]).to_string(),
                })
                .collect(),
        ));
    }
    let cambios: Vec<f64> = resultado
        .cambios_por_posicion
        .iter()
        .map(|&cambio| cambio as f64)
        .collect();
    if largo <= 60 {
        salida.agregar(Bloque::barras(
            "Mutaciones por posición (contando desde 0)",
            "Posición",
            "Mutaciones",
            (0..largo).map(|posicion| posicion.to_string()).collect(),
            cambios,
        ));
    } else {
        let puntos: Vec<[f64; 2]> = cambios
            .iter()
            .enumerate()
            .map(|(posicion, &cambio)| [posicion as f64, cambio])
            .collect();
        salida.agregar(Bloque::linea(
            "Mutaciones por posición (contando desde 0)",
            "Posición",
            "Mutaciones",
            &puntos,
            Vec::new(),
        ));
    }
    let mut filas = Vec::new();
    for (a, b, cambios_en_rama) in arbol.aristas() {
        if cambios_en_rama > 0.0 {
            let posiciones: Vec<String> = (0..largo)
                .filter(|&posicion| {
                    resultado.secuencias[a][posicion] != resultado.secuencias[b][posicion]
                })
                .map(|posicion| posicion.to_string())
                .collect();
            filas.push(vec![
                arbol.nombres[a].clone(),
                arbol.nombres[b].clone(),
                miles(cambios_en_rama as usize),
                posiciones.join(", "),
            ]);
        }
    }
    if !filas.is_empty() {
        salida.agregar(Bloque::tabla(
            "Ramas con mutaciones (posiciones contando desde 0)",
            &["Un extremo", "Otro extremo", "Mutaciones", "Posiciones"],
            filas,
        ));
        salida.explicar(format!(
            "El puntaje es {}: ninguna forma de poner letras en los nodos internos explica estas hojas con menos mutaciones en este árbol.",
            miles(puntaje)
        ));
    }
    if puntaje == 0 {
        salida.explicar(
            "Todas las hojas tienen la misma letra en cada posición: no hace falta ninguna mutación y los antepasados son iguales a ellas.",
        );
    } else if puntaje == minimo {
        salida.explicar(
            "Índice de consistencia = 1: cada letra nueva aparece una sola vez en el árbol, sin cambios repetidos ni reversiones (no hay homoplasia).",
        );
    } else {
        salida.explicar(format!(
            "Índice de consistencia = {} / {} = {}: hay {} mutación(es) de más que obligan a repetir un mismo cambio en dos ramas o a revertirlo (homoplasia). Otro árbol podría explicar mejor los datos: compáralo en «Parsimonia grande».",
            miles(minimo),
            miles(puntaje),
            numero(minimo as f64 / puntaje as f64),
            miles(puntaje - minimo)
        ));
    }
    let conservadas = resultado
        .cambios_por_posicion
        .iter()
        .filter(|&&cambio| cambio == 0)
        .count();
    if let Some((posicion, &mayor)) = resultado
        .cambios_por_posicion
        .iter()
        .enumerate()
        .max_by_key(|&(posicion, &cambio)| (cambio, std::cmp::Reverse(posicion)))
        && mayor > 0
    {
        salida.explicar(format!(
            "{} de {} posiciones no cambian nunca; la que más cambia es la {} (contando desde 0), con {}.",
            miles(conservadas),
            miles(largo),
            posicion,
            mutaciones(mayor)
        ));
    }
}

fn orden_desde(arbol: &Arbol, raiz: usize) -> (Vec<usize>, Vec<Option<usize>>) {
    let mut orden = Vec::new();
    let mut padres = vec![None; arbol.nombres.len()];
    let mut visitados = vec![false; arbol.nombres.len()];
    let mut pila = vec![raiz];
    visitados[raiz] = true;
    while let Some(nodo) = pila.pop() {
        orden.push(nodo);
        for &(vecino, _) in arbol.vecinos[nodo].iter().rev() {
            if !visitados[vecino] {
                visitados[vecino] = true;
                padres[vecino] = Some(nodo);
                pila.push(vecino);
            }
        }
    }
    (orden, padres)
}

fn parsimonia_pequena(entrada: &Entrada) -> Result<Salida, String> {
    let leido = filogenia::leer_newick(&entrada.texto("arbol")?)?;
    let (nombres, secuencias) = leer_secuencias(entrada, "secuencias", false)?;
    let por_nodo = secuencias_en_el_arbol(&leido.arbol, &nombres, secuencias)?;
    let resultado = filogenia::parsimonia_pequena(&leido.arbol, leido.raiz, &por_nodo)?;
    let mut arbol = resultado.arbol.clone();
    nombres_unicos(&mut arbol, Some(leido.raiz));
    let mut salida = Salida::default();
    informe_de_parsimonia(&mut salida, &resultado, &arbol, &por_nodo, leido.raiz);
    salida.explicar(format!(
        "La secuencia de «{}» es una reconstrucción del antepasado común de todas las hojas; donde dos letras empatan, cualquiera daría el mismo puntaje.",
        arbol.nombres[leido.raiz]
    ));
    Ok(salida)
}

fn parsimonia_sin_raiz(entrada: &Entrada) -> Result<Salida, String> {
    let leido = filogenia::leer_newick(&entrada.texto("arbol")?)?;
    let tenia_raiz = leido.arbol.grado(leido.raiz) == 2;
    let arbol = filogenia::quitar_raiz(&leido);
    let (nombres, secuencias) = leer_secuencias(entrada, "secuencias", false)?;
    let por_nodo = secuencias_en_el_arbol(&arbol, &nombres, secuencias)?;
    let resultado = filogenia::parsimonia_sin_raiz(&arbol, &por_nodo)?;
    let mut con_nombres = resultado.arbol.clone();
    nombres_unicos(&mut con_nombres, None);
    let mut salida = Salida::default();
    let raiz = raiz_para_dibujar(&con_nombres);
    informe_de_parsimonia(&mut salida, &resultado, &con_nombres, &por_nodo, raiz);
    salida.explicar(
        "En un árbol sin raíz el puntaje no depende de dónde se ponga la raíz: se calculó enraizando en un nodo interno cualquiera, y cualquier otro daría lo mismo.",
    );
    if tenia_raiz {
        salida.explicar(
            "El Newick tenía raíz (dos ramas arriba): se quitó, uniendo esas dos ramas en una sola.",
        );
    }
    Ok(salida)
}

fn leer_arista(texto: &str, arbol: &Arbol) -> Result<(usize, usize), String> {
    let nombres: Vec<&str> = texto
        .split(|c: char| c.is_whitespace() || c == '-' || c == '>' || c == ',' || c == ';')
        .filter(|pieza| !pieza.is_empty())
        .collect();
    if nombres.len() != 2 {
        return Err("Escribe la arista como dos nombres de nodo, por ejemplo «n2 n3».".to_string());
    }
    let mut nodos = [0; 2];
    for (posicion, nombre) in nombres.iter().enumerate() {
        nodos[posicion] = arbol
            .buscar_nodo(nombre)
            .ok_or_else(|| format!("El árbol no tiene un nodo llamado «{}».", recortar(nombre)))?;
    }
    Ok((nodos[0], nodos[1]))
}

fn lado_de(arbol: &Arbol, nodo: usize, excluido: usize) -> Vec<usize> {
    arbol.vecinos[nodo]
        .iter()
        .map(|&(vecino, _)| vecino)
        .filter(|&vecino| vecino != excluido)
        .collect()
}

fn vecinos_cercanos(entrada: &Entrada) -> Result<Salida, String> {
    let arbol = filogenia::leer_aristas(&entrada.texto("arbol")?)?;
    if arbol.hojas().len() > MAXIMO_DE_ESPECIES {
        return Err(format!(
            "El árbol tiene más de {MAXIMO_DE_ESPECIES} hojas; ese es el máximo aquí."
        ));
    }
    let (a, b) = leer_arista(&entrada.texto("arista")?, &arbol)?;
    let (primero, segundo) = filogenia::vecinos_cercanos(&arbol, a, b)?;
    let lado_a = lado_de(&arbol, a, b);
    let lado_b = lado_de(&arbol, b, a);
    let nombre = |nodo: usize| arbol.nombres[nodo].as_str();
    let (w, x, y, z) = (
        nombre(lado_a[0]),
        nombre(lado_a[1]),
        nombre(lado_b[0]),
        nombre(lado_b[1]),
    );
    let mut salida = Salida::default();
    salida.agregar(Bloque::arbol(
        "Árbol original",
        newick_para_dibujar(&arbol, a),
    ));
    salida.agregar(Bloque::arbol(
        &format!("Vecino 1: «{x}» cambia de lugar con «{y}»"),
        newick_para_dibujar(&primero, a),
    ));
    salida.agregar(Bloque::arbol(
        &format!("Vecino 2: «{x}» cambia de lugar con «{z}»"),
        newick_para_dibujar(&segundo, a),
    ));
    salida.agregar(Bloque::lista(
        "Lista de adyacencia del vecino 1",
        lista_de_adyacencia(&primero),
    ));
    salida.agregar(Bloque::lista(
        "Lista de adyacencia del vecino 2",
        lista_de_adyacencia(&segundo),
    ));
    salida.explicar(format!(
        "La arista «{}»–«{}» separa cuatro subárboles: los que cuelgan de «{w}» y «{x}» de un lado, y de «{y}» y «{z}» del otro. Cada vecino cercano cambia un subárbol de un lado por uno del otro; mover «{w}» en vez de «{x}» daría estos mismos dos árboles.",
        nombre(a),
        nombre(b)
    ));
    let internas = arbol
        .aristas()
        .iter()
        .filter(|&&(u, v, _)| arbol.grado(u) > 1 && arbol.grado(v) > 1)
        .count();
    salida.explicar(format!(
        "Este árbol tiene {internas} arista(s) interna(s), así que tiene {} vecinos cercanos en total: son los que revisa la parsimonia grande en cada paso.",
        2 * internas
    ));
    if arbol.aristas().iter().any(|arista| arista.2 != 1.0) {
        salida.explicar("Cada subárbol se lleva consigo la longitud de su rama.");
    }
    Ok(salida)
}

fn arboles_posibles(hojas: usize) -> String {
    if hojas < 3 {
        return "1".to_string();
    }
    let mut logaritmo = 0.0;
    let mut exacto: u128 = 1;
    for impar in (3..=(2 * hojas - 5)).step_by(2) {
        logaritmo += (impar as f64).log10();
        exacto = exacto.saturating_mul(impar as u128);
    }
    if logaritmo < 15.0 {
        miles(exacto as usize)
    } else {
        format!("más de 10^{}", logaritmo.floor() as u64)
    }
}

fn parsimonia_grande(entrada: &Entrada) -> Result<Salida, String> {
    let (nombres, secuencias) = leer_secuencias(entrada, "secuencias", false)?;
    if nombres.len() < 3 {
        return Err("Hacen falta al menos tres secuencias.".to_string());
    }
    revisar_mismo_largo(&nombres, &secuencias)?;
    let inicio = entrada.opcion("inicio", INICIOS)?;
    let semilla = entrada.natural("semilla")? as u64;
    let maximo = entrada.natural("movimientos")?;
    if maximo > MAXIMO_DE_MOVIMIENTOS {
        return Err(format!(
            "El máximo de movimientos debe estar entre 0 y {MAXIMO_DE_MOVIMIENTOS}."
        ));
    }
    let al_azar = inicio == "al azar";
    let inicial = if al_azar {
        filogenia::arbol_al_azar(&nombres, semilla)?
    } else {
        let distancias_p = filogenia::matriz_de_distancias_p(&secuencias)?;
        filogenia::union_de_vecinos(&distancias_p, &nombres)?
    };
    let mut por_nodo = secuencias;
    por_nodo.resize(inicial.nombres.len(), Vec::new());
    let busqueda = filogenia::parsimonia_grande(&inicial, &por_nodo, maximo)?;
    let cantidad = nombres.len();
    let mut arbol_inicial = inicial.clone();
    nombres_unicos(&mut arbol_inicial, None);
    let mut puntos = vec![[0.0, busqueda.puntaje_inicial as f64]];
    let mut filas = Vec::new();
    for (numero_de_movimiento, movimiento) in busqueda.movimientos.iter().enumerate() {
        puntos.push([(numero_de_movimiento + 1) as f64, movimiento.puntaje as f64]);
        let (a, b) = movimiento.arista;
        let (x, y) = movimiento.intercambio;
        let nombre = |nodo: usize| arbol_inicial.nombres[nodo].as_str();
        filas.push(vec![
            (numero_de_movimiento + 1).to_string(),
            format!("{}–{}", nombre(a), nombre(b)),
            format!("{} ↔ {}", nombre(x), nombre(y)),
            miles(movimiento.puntaje),
        ]);
    }
    let puntaje_final = busqueda.parsimonia.puntaje;
    let mut final_con_cambios = busqueda.parsimonia.arbol.clone();
    nombres_unicos(&mut final_con_cambios, None);
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        "Puntaje inicial",
        mutaciones(busqueda.puntaje_inicial),
    ));
    salida.agregar(Bloque::valor("Puntaje final", mutaciones(puntaje_final)));
    salida.agregar(Bloque::valor(
        "Movimientos que mejoraron",
        miles(busqueda.movimientos.len()),
    ));
    salida.agregar(Bloque::linea(
        "Puntaje después de cada movimiento",
        "Movimiento",
        "Mutaciones",
        &puntos,
        Vec::new(),
    ));
    salida.agregar(Bloque::arbol(
        "Árbol final (cambios en cada rama)",
        filogenia::escribir_newick(&final_con_cambios, raiz_para_dibujar(&final_con_cambios)),
    ));
    salida.agregar(Bloque::lista(
        "Lista de adyacencia del árbol final (cambios en cada rama)",
        lista_de_adyacencia(&final_con_cambios),
    ));
    salida.agregar(Bloque::arbol(
        "Árbol inicial",
        newick_para_dibujar(&arbol_inicial, raiz_para_dibujar(&arbol_inicial)),
    ));
    salida.agregar(Bloque::tabla(
        "Movimientos",
        &["Movimiento", "Arista", "Intercambio", "Puntaje"],
        filas,
    ));
    let origen = if al_azar {
        format!("un árbol al azar (semilla {semilla})")
    } else {
        "unión de vecinos sobre las distancias p".to_string()
    };
    if busqueda.movimientos.is_empty() && maximo > 0 {
        salida.explicar(format!(
            "Ningún vecino cercano del árbol inicial ({origen}) baja el puntaje: ya era un óptimo local, con {}.",
            mutaciones(puntaje_final)
        ));
    } else if !busqueda.movimientos.is_empty() {
        salida.explicar(format!(
            "Desde {origen}, la búsqueda hizo {} intercambio(s) de vecinos que bajaron el puntaje de {} a {}.",
            busqueda.movimientos.len(),
            miles(busqueda.puntaje_inicial),
            mutaciones(puntaje_final)
        ));
    }
    if maximo == 0 {
        salida.explicar("Con 0 movimientos solo se midió el árbol inicial.");
    } else if busqueda.limite_alcanzado {
        salida.explicar(format!(
            "Se detuvo al llegar al límite de {maximo} movimientos: quizás quedaban mejoras; sube el límite."
        ));
    }
    match cantidad {
        3 => salida.explicar("Con 3 especies hay un solo árbol sin raíz posible."),
        4 if maximo > 0 => salida.explicar(
            "Con 4 especies hay solo 3 árboles sin raíz posibles y los vecinos cercanos los cubren todos: este es el óptimo global.",
        ),
        _ => salida.explicar(format!(
            "Con {cantidad} especies hay {} árboles sin raíz posibles. La búsqueda solo mira vecinos cercanos, así que puede quedar atrapada en un óptimo local: prueba «al azar» con otras semillas y compara puntajes.",
            arboles_posibles(cantidad)
        )),
    }
    let hojas: Vec<&Vec<u8>> = por_nodo
        .iter()
        .filter(|secuencia| !secuencia.is_empty())
        .collect();
    let minimo = minimo_de_cambios(&hojas, busqueda.parsimonia.cambios_por_posicion.len());
    if puntaje_final > 0 {
        salida.explicar(format!(
            "Índice de consistencia del árbol final = {} / {} = {}{}",
            miles(minimo),
            miles(puntaje_final),
            numero(minimo as f64 / puntaje_final as f64),
            if minimo == puntaje_final {
                ": no queda homoplasia."
            } else {
                ": aún hay cambios repetidos o revertidos que ningún vecino cercano elimina."
            }
        ));
    }
    if !al_azar {
        salida.explicar("La semilla solo se usa con el inicio «al azar».");
    }
    Ok(salida)
}

const MAXIMO_DE_PASOS: usize = 400;
const LIMITE_DE_TABLA: usize = 10;
const LIMITE_DE_LISTA: usize = 14;
const COSTO_IMPOSIBLE: u32 = 1_000_000_000;

fn abreviar(texto: &str, maximo: usize) -> String {
    let mut corto: String = texto.chars().take(maximo).collect();
    if texto.chars().count() > maximo {
        corto.push('…');
    }
    corto
}

fn lista_en_texto(elementos: &[String], maximo: usize) -> String {
    let mut partes: Vec<String> = elementos
        .iter()
        .take(maximo)
        .map(|elemento| abreviar(elemento, 10))
        .collect();
    if elementos.len() > maximo {
        partes.push("…".to_string());
    }
    format!("[{}]", partes.join(", "))
}

fn conteos(valores: &[usize]) -> String {
    abreviar(&format!("{valores:?}"), 40)
}

fn mostrar(indice: usize, total: usize, detallados: usize, resumen: usize) -> bool {
    if indice < detallados || indice + 1 >= total {
        return true;
    }
    let restantes = total - detallados;
    let salto = restantes.div_ceil(resumen.max(1)).max(1);
    (indice - detallados).is_multiple_of(salto)
}

fn con_matriz(
    pizarra: Pizarra,
    titulo: &str,
    nombres: &[String],
    matriz: &[Vec<f64>],
    resaltar: &[(usize, usize)],
) -> Pizarra {
    let visibles = nombres.len().min(matriz.len()).min(LIMITE_DE_TABLA);
    let mut encabezado = vec![String::new()];
    encabezado.extend(
        nombres
            .iter()
            .take(visibles)
            .map(|nombre| abreviar(nombre, 9)),
    );
    let columnas: Vec<&str> = encabezado.iter().map(String::as_str).collect();
    let filas: Vec<Vec<String>> = (0..visibles)
        .map(|i| {
            let mut fila = vec![abreviar(&nombres[i], 9)];
            fila.extend(matriz[i].iter().take(visibles).map(|&valor| numero(valor)));
            fila
        })
        .collect();
    let marcadas: Vec<(usize, usize)> = resaltar
        .iter()
        .filter(|&&(i, j)| i < visibles && j < visibles)
        .map(|&(i, j)| (i, j + 1))
        .collect();
    let titulo = if visibles < matriz.len() {
        format!("{titulo} (primeras {visibles} de {})", matriz.len())
    } else {
        titulo.to_string()
    };
    pizarra.tabla(&titulo, &columnas, filas, &marcadas)
}

fn lista_de_aristas(
    arbol: &Arbol,
    etiquetas: &[String],
    incluir: &dyn Fn(usize, usize) -> bool,
    resaltar: &dyn Fn(usize, usize) -> bool,
    maximo: usize,
) -> (Vec<String>, Vec<usize>) {
    let mut lineas = Vec::new();
    let mut marcadas = Vec::new();
    for (a, b, longitud) in arbol.aristas() {
        if lineas.len() == maximo {
            lineas.push("…".to_string());
            break;
        }
        if incluir(a, b) {
            if resaltar(a, b) {
                marcadas.push(lineas.len());
            }
            lineas.push(format!(
                "{}–{}: {}",
                abreviar(&etiquetas[a], 10),
                abreviar(&etiquetas[b], 10),
                numero(longitud)
            ));
        }
    }
    (lineas, marcadas)
}

fn letras_en_texto(bits: u8) -> String {
    let letras: Vec<String> = "ACGT"
        .chars()
        .enumerate()
        .filter(|&(indice, _)| bits & (1 << indice) != 0)
        .map(|(_, letra)| letra.to_string())
        .collect();
    format!("{{{}}}", letras.join(","))
}

fn es_letra(letra: u8) -> bool {
    matches!(letra, b'A' | b'C' | b'G' | b'T')
}

fn recorrido_desde(arbol: &Arbol, origen: usize) -> (Vec<f64>, Vec<(usize, usize, f64)>) {
    let mut distancias = vec![f64::INFINITY; arbol.nombres.len()];
    let mut visitados = vec![false; arbol.nombres.len()];
    let mut eventos = Vec::new();
    distancias[origen] = 0.0;
    visitados[origen] = true;
    let mut pila = vec![origen];
    while let Some(nodo) = pila.pop() {
        for &(vecino, longitud) in &arbol.vecinos[nodo] {
            if !visitados[vecino] {
                visitados[vecino] = true;
                distancias[vecino] = distancias[nodo] + longitud;
                eventos.push((nodo, vecino, longitud));
                pila.push(vecino);
            }
        }
    }
    (distancias, eventos)
}

fn narrar_distancias_entre_hojas(entrada: &Entrada) -> Result<(Vec<Paso>, Vec<Vec<f64>>), String> {
    let arbol = filogenia::leer_aristas(&entrada.texto("arbol")?)?;
    let hojas = arbol.hojas();
    if hojas.len() < 2 {
        return Err("El árbol necesita al menos dos hojas.".to_string());
    }
    if hojas.len() > MAXIMO_DE_ESPECIES {
        return Err(format!(
            "El árbol tiene {} hojas; el máximo aquí es {MAXIMO_DE_ESPECIES}.",
            miles(hojas.len())
        ));
    }
    let nombres: Vec<String> = hojas
        .iter()
        .map(|&hoja| arbol.nombres[hoja].clone())
        .collect();
    let total = hojas.len();
    let mut narracion = Narracion::nueva(MAXIMO_DE_PASOS);
    let (aristas, _) = lista_de_aristas(
        &arbol,
        &arbol.nombres,
        &|_, _| true,
        &|_, _| false,
        LIMITE_DE_LISTA * 2,
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable(
                "arbol",
                format!(
                    "{} nodos, {} aristas",
                    arbol.nombres.len(),
                    arbol.aristas().len()
                ),
            )
            .variable("hojas", lista_en_texto(&nombres, 6))
            .lista("Aristas del árbol (nodo–nodo: longitud)", aristas, &[])
            .paso(
                "let mut matriz = Vec::with_capacity(hojas.len());",
                format!(
                    "Tu árbol tiene {total} hojas. Queremos una tabla de {total}×{total} con la distancia entre cada par, y esa distancia es la suma de las longitudes de las ramas del camino que las une. La vamos a llenar fila por fila, una hoja a la vez."
                ),
            ),
    );
    let mut matriz: Vec<Vec<f64>> = Vec::with_capacity(total);
    for (indice, &hoja) in hojas.iter().enumerate() {
        let (distancias, visitas) = recorrido_desde(&arbol, hoja);
        let nombre = abreviar(&arbol.nombres[hoja], 12);
        let hechas = format!("{indice} de {total} filas");
        if indice == 0 {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("matriz", &hechas)
                    .variable("hoja", &nombre)
                    .paso(
                        "let desde_la_hoja = distancias_desde(arbol, hoja);",
                        format!(
                            "Para la fila de «{nombre}» recorremos el árbol entero partiendo de ella, y vamos anotando a qué distancia queda cada nodo."
                        ),
                    ),
            );
            let mut conocidas = vec![format!("{nombre}: 0")];
            narracion.agregar(
                Pizarra::nueva()
                    .variable("origen", &nombre)
                    .lista("Distancia desde la hoja", conocidas.clone(), &[0])
                    .paso(
                        "distancias[origen] = 0.0;",
                        format!(
                            "«{nombre}» está a distancia 0 de sí misma. Es el punto de partida."
                        ),
                    ),
            );
            let mut hasta_ahora = vec![f64::INFINITY; arbol.nombres.len()];
            hasta_ahora[hoja] = 0.0;
            for &(nodo, vecino, longitud) in visitas.iter().take(12) {
                hasta_ahora[vecino] = hasta_ahora[nodo] + longitud;
                conocidas.push(format!(
                    "{}: {}",
                    abreviar(&arbol.nombres[vecino], 10),
                    numero(hasta_ahora[vecino])
                ));
                let ultimo = conocidas.len() - 1;
                let (de_nodo, de_vecino) = (
                    abreviar(&arbol.nombres[nodo], 12),
                    abreviar(&arbol.nombres[vecino], 12),
                );
                narracion.agregar(
                    Pizarra::nueva()
                        .variable("origen", &nombre)
                        .variable("nodo", &de_nodo)
                        .variable("vecino", &de_vecino)
                        .variable("longitud", numero(longitud))
                        .lista("Distancia desde la hoja", conocidas.clone(), &[ultimo])
                        .paso(
                            "distancias[vecino] = distancias[nodo] + longitud;",
                            format!(
                                "«{de_vecino}» queda a {} + {} = {} de «{nombre}»: lo que ya sabíamos hasta «{de_nodo}» más la rama que las une.",
                                numero(hasta_ahora[nodo]),
                                numero(longitud),
                                numero(hasta_ahora[vecino])
                            ),
                        ),
                );
            }
            if visitas.len() > 12 {
                narracion.agregar(
                    Pizarra::nueva()
                        .variable("origen", &nombre)
                        .lista("Distancia desde la hoja", conocidas.clone(), &[])
                        .paso(
                            "distancias[vecino] = distancias[nodo] + longitud;",
                            format!(
                                "El recorrido sigue igual, nodo por nodo, hasta tocar los {} nodos del árbol.",
                                miles(arbol.nombres.len())
                            ),
                        ),
                );
            }
        }
        let fila: Vec<f64> = hojas.iter().map(|&otra| distancias[otra]).collect();
        if mostrar(indice, total, 2, 10) && narracion.hay_lugar() {
            let elementos: Vec<String> = nombres
                .iter()
                .zip(&fila)
                .take(LIMITE_DE_LISTA)
                .map(|(otra, &distancia)| format!("{}: {}", abreviar(otra, 10), numero(distancia)))
                .collect();
            let propia = indice.min(LIMITE_DE_LISTA - 1);
            narracion.agregar(
                Pizarra::nueva()
                    .variable("matriz", &hechas)
                    .variable("hoja", &nombre)
                    .lista(
                        &format!("Fila de «{nombre}»: distancia a cada hoja"),
                        elementos,
                        &[propia],
                    )
                    .paso(
                        "let fila: Vec<f64> = hojas.iter().map(|&otra| desde_la_hoja[otra]).collect();",
                        format!(
                            "De todo lo que anotamos nos quedamos solo con las hojas: esa es la fila de «{nombre}». Su distancia consigo misma es 0, y las demás vienen del recorrido."
                        ),
                    ),
            );
        }
        matriz.push(fila);
    }
    let ultimo = con_matriz(
        Pizarra::nueva().variable("matriz", format!("{total} de {total} filas")),
        "Distancias entre hojas",
        &nombres,
        &matriz,
        &[],
    )
    .paso(
        "matriz.push(fila);",
        format!(
            "Con la última fila la tabla de {total}×{total} está completa. Es simétrica: la distancia de A a B es la de B a A. Mira «Qué significa» para ver qué pares quedaron más cerca y más lejos."
        ),
    );
    Ok((narracion.terminar(ultimo), matriz))
}

fn pasos_distancias_entre_hojas(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    narrar_distancias_entre_hojas(entrada).map(|(pasos, _)| pasos)
}

type RamaConPareja = (f64, usize, usize);

fn narrar_longitud_de_rama(entrada: &Entrada) -> Result<(Vec<Paso>, RamaConPareja), String> {
    let (nombres, matriz) = leer_matriz(entrada, "matriz")?;
    let hoja = buscar_especie(&entrada.texto("hoja")?, &nombres)?;
    filogenia::longitud_de_rama_con_par(&matriz, hoja)?;
    let cantidad = nombres.len();
    let nombre = abreviar(&nombres[hoja], 12);
    let fila_y_columna: Vec<(usize, usize)> = (0..cantidad)
        .flat_map(|otra| [(hoja, otra), (otra, hoja)])
        .collect();
    let mut narracion = Narracion::nueva(MAXIMO_DE_PASOS);
    narracion.agregar(
        con_matriz(
            Pizarra::nueva()
                .variable("hoja", format!("{hoja} («{nombre}»)"))
                .variable("cantidad", cantidad),
            "Distancias (fila y columna de la hoja)",
            &nombres,
            &matriz,
            &fila_y_columna,
        )
        .paso(
            "validar_distancias(distancias)?;",
            format!(
                "Tienes {cantidad} especies y quieres la rama de «{nombre}»: la parte de su camino que es solo suya, lo que cambió por su cuenta. Primero se comprueba que la tabla sea cuadrada, simétrica y con ceros en la diagonal; después se calcula todo con esos números, sin dibujar ningún árbol."
            ),
        ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("hoja", format!("{hoja} («{nombre}»)"))
            .variable("mejor", "(∞, hoja, hoja)")
            .paso(
                "let mut mejor = (f64::INFINITY, hoja, hoja);",
                "Probaremos muchas parejas de otras hojas y guardaremos la candidata más chica. Empieza en ∞ para que la primera que probemos siempre la supere.",
            ),
    );
    let mut mejor = (f64::INFINITY, hoja, hoja);
    let mut contadas = 0usize;
    let mut mejoras_mostradas = 0usize;
    let mut candidatas: Vec<String> = Vec::new();
    let mut posicion_de_la_mejor: Option<usize> = None;
    for i in (0..cantidad).filter(|&i| i != hoja) {
        for k in (i..cantidad).filter(|&k| k != hoja) {
            let candidata = (matriz[i][hoja] + matriz[hoja][k] - matriz[i][k]) / 2.0;
            let anterior = mejor;
            let mejora = candidata < mejor.0;
            if contadas < LIMITE_DE_LISTA {
                candidatas.push(format!(
                    "i={}, k={}: {}",
                    abreviar(&nombres[i], 8),
                    abreviar(&nombres[k], 8),
                    numero(candidata)
                ));
            }
            if mejora {
                mejor = (candidata, i, k);
                if contadas < LIMITE_DE_LISTA {
                    posicion_de_la_mejor = Some(contadas);
                }
            }
            let mut resaltadas = Vec::new();
            if contadas < LIMITE_DE_LISTA {
                resaltadas.push(contadas);
            }
            resaltadas.extend(posicion_de_la_mejor);
            let nombre_i = abreviar(&nombres[i], 12);
            let nombre_k = abreviar(&nombres[k], 12);
            if contadas < 6 {
                let veredicto = if mejora && anterior.0.is_infinite() {
                    "Es la primera, así que por ahora es la mejor.".to_string()
                } else if mejora {
                    format!("Es menor que la mejor anterior ({}).", numero(anterior.0))
                } else {
                    format!("No supera a la mejor ({}).", numero(anterior.0))
                };
                narracion.agregar(
                    con_matriz(
                        Pizarra::nueva()
                            .variable("hoja", format!("{hoja} («{nombre}»)"))
                            .variable("mejor", format!("({}, i, k)", numero(anterior.0)))
                            .variable("i", format!("{i} («{nombre_i}»)"))
                            .variable("k", format!("{k} («{nombre_k}»)"))
                            .variable("candidata", numero(candidata)),
                        "Distancias que intervienen",
                        &nombres,
                        &matriz,
                        &[(i, hoja), (hoja, k), (i, k)],
                    )
                    .lista("Candidatas probadas", candidatas.clone(), &resaltadas)
                    .paso(
                        "let candidata = (distancias[i][hoja] + distancias[hoja][k] - distancias[i][k]) / 2.0;",
                        format!(
                            "Con i = «{nombre_i}» y k = «{nombre_k}»: (D(i,hoja) + D(hoja,k) − D(i,k)) / 2 = ({} + {} − {}) / 2 = {}. {veredicto}",
                            numero(matriz[i][hoja]),
                            numero(matriz[hoja][k]),
                            numero(matriz[i][k]),
                            numero(candidata)
                        ),
                    ),
                );
            }
            if mejora && (contadas < 6 || mejoras_mostradas < 20) {
                mejoras_mostradas += 1;
                narracion.agregar(
                    Pizarra::nueva()
                        .variable("hoja", format!("{hoja} («{nombre}»)"))
                        .variable("mejor", format!("({}, {i}, {k})", numero(candidata)))
                        .variable("i", format!("{i} («{nombre_i}»)"))
                        .variable("k", format!("{k} («{nombre_k}»)"))
                        .variable("candidata", numero(candidata))
                        .lista("Candidatas probadas", candidatas.clone(), &resaltadas)
                        .paso(
                            "mejor = (candidata, i, k);",
                            format!(
                                "{} es la menor candidata hasta ahora, así que la guardamos junto con la pareja que la dio: i = «{nombre_i}», k = «{nombre_k}».",
                                numero(candidata)
                            ),
                        ),
                );
            }
            contadas += 1;
        }
    }
    let nombre_i = abreviar(&nombres[mejor.1], 12);
    let nombre_k = abreviar(&nombres[mejor.2], 12);
    narracion.agregar(
        Pizarra::nueva()
            .variable("hoja", format!("{hoja} («{nombre}»)"))
            .variable("mejor", format!("({}, {}, {})", numero(mejor.0), mejor.1, mejor.2))
            .lista("Candidatas probadas", candidatas.clone(), &posicion_de_la_mejor.into_iter().collect::<Vec<_>>())
            .paso(
                "for k in (i..cantidad).filter(|&k| k != hoja) {",
                format!(
                    "Ya probamos las {} parejas posibles. La menor candidata es {}, con i = «{nombre_i}» y k = «{nombre_k}».",
                    miles(contadas),
                    numero(mejor.0)
                ),
            ),
    );
    let rama = mejor.0;
    let significado = if rama < -1e-9 {
        "Salió negativa, algo imposible en un árbol real: señal de que estas distancias no vienen de un árbol."
    } else if casi_cero(rama) {
        "Es 0: la hoja cae justo sobre un nodo interno."
    } else {
        "Es lo que esa especie cambió por su cuenta desde que se separó del resto."
    };
    let ultimo = Pizarra::nueva()
        .variable("hoja", format!("{hoja} («{nombre}»)"))
        .variable("mejor", format!("({}, {}, {})", numero(mejor.0), mejor.1, mejor.2))
        .paso(
            "Ok(rama_entre_las_primeras(distancias, hoja, cantidad))",
            format!(
                "La rama de «{nombre}» mide {}. {significado} Mira «Qué significa» para compararla con su pariente más cercano.",
                numero(rama)
            ),
        );
    Ok((narracion.terminar(ultimo), mejor))
}

fn pasos_longitud_de_rama(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    narrar_longitud_de_rama(entrada).map(|(pasos, _)| pasos)
}

fn narrar_filogenia_aditiva(entrada: &Entrada) -> Result<(Vec<Paso>, (Arbol, bool)), String> {
    let (nombres, matriz) = leer_matriz(entrada, "matriz")?;
    let cantidad = nombres.len();
    let mut narracion = Narracion::nueva(MAXIMO_DE_PASOS);
    narracion.agregar(
        con_matriz(
            Pizarra::nueva().variable("cantidad", cantidad),
            "Matriz de distancias",
            &nombres,
            &matriz,
            &[],
        )
        .paso(
            "validar_distancias(distancias)?;",
            format!(
                "Tienes la tabla de distancias de {cantidad} especies. Queremos un árbol cuyas ramas sumen exactamente esas distancias. La idea: quitar la última hoja, resolver el problema más chico y volver a colgarla; así hasta quedarnos con dos hojas."
            ),
        ),
    );
    let mut arbol = Arbol::default();
    for nombre in &nombres {
        arbol.agregar_nodo(nombre);
    }
    let mut recortada = matriz.clone();
    let mut pendientes: Vec<(usize, f64, usize, usize, f64)> = Vec::new();
    let vueltas = cantidad.saturating_sub(2);
    for (vuelta, hoja) in (2..cantidad).rev().enumerate() {
        let nombre = abreviar(&nombres[hoja], 12);
        let (rama, _, _) = filogenia::rama_entre_las_primeras(&recortada, hoja, hoja + 1);
        let propias: Vec<(usize, usize)> = (0..=hoja)
            .flat_map(|otra| [(hoja, otra), (otra, hoja)])
            .collect();
        let detallada = vuelta < 2;
        if detallada {
            narracion.agregar(
                con_matriz(
                    Pizarra::nueva()
                        .variable("recortada", format!("{cantidad}×{cantidad}"))
                        .variable("hoja", format!("{hoja} («{nombre}»)"))
                        .variable("rama", numero(rama)),
                    "Distancias recortadas",
                    &nombres,
                    &recortada,
                    &propias,
                )
                .paso(
                    "let (rama, _, _) = rama_entre_las_primeras(&recortada, hoja, hoja + 1);",
                    format!(
                        "{} «{nombre}». Su rama propia es lo menor de (D(i,hoja) + D(hoja,k) − D(i,k)) / 2 entre las hojas anteriores: aquí {}.",
                        if vuelta == 0 {
                            "Empezamos por la última hoja,"
                        } else {
                            "Seguimos con la anterior,"
                        },
                        numero(rama)
                    ),
                ),
            );
        }
        for fila in recortada.iter_mut().take(hoja) {
            fila[hoja] -= rama;
        }
        for distancia in recortada[hoja].iter_mut().take(hoja) {
            *distancia -= rama;
        }
        let (desde, hacia) = filogenia::hojas_que_rodean(&recortada, hoja);
        let distancia = recortada[desde][hoja];
        let nombre_desde = abreviar(&nombres[desde], 12);
        let nombre_hacia = abreviar(&nombres[hacia], 12);
        if detallada {
            narracion.agregar(
                con_matriz(
                    Pizarra::nueva()
                        .variable("recortada", format!("{cantidad}×{cantidad}"))
                        .variable("hoja", format!("{hoja} («{nombre}»)"))
                        .variable("rama", numero(rama)),
                    "Distancias recortadas",
                    &nombres,
                    &recortada,
                    &propias,
                )
                .paso(
                    "fila[hoja] -= rama;",
                    format!(
                        "Le restamos {} a toda la fila y la columna de «{nombre}»: es como cortarle su rama propia. Lo que queda es la distancia hasta el punto donde se engancha al resto.",
                        numero(rama)
                    ),
                ),
            );
            narracion.agregar(
                con_matriz(
                    Pizarra::nueva()
                        .variable("hoja", format!("{hoja} («{nombre}»)"))
                        .variable("rama", numero(rama))
                        .variable("desde", format!("{desde} («{nombre_desde}»)"))
                        .variable("hacia", format!("{hacia} («{nombre_hacia}»)")),
                    "Distancias recortadas",
                    &nombres,
                    &recortada,
                    &[(desde, hoja), (hoja, hacia), (desde, hacia)],
                )
                .paso(
                    "let (desde, hacia) = hojas_que_rodean(&recortada, hoja);",
                    format!(
                        "Buscamos dos hojas tales que D(desde,hoja) + D(hoja,hacia) = D(desde,hacia): entonces «{nombre}» cuelga del camino que las une. Salen «{nombre_desde}» y «{nombre_hacia}», y cuelga a {} de «{nombre_desde}».",
                        numero(distancia)
                    ),
                ),
            );
        } else if mostrar(vuelta, vueltas, 2, 8) {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("hoja", format!("{hoja} («{nombre}»)"))
                    .variable("rama", numero(rama))
                    .variable("desde", format!("{desde} («{nombre_desde}»)"))
                    .variable("hacia", format!("{hacia} («{nombre_hacia}»)"))
                    .paso(
                        "pendientes.push(HojaPendiente {",
                        format!(
                            "Lo mismo con «{nombre}»: rama propia {}, y se engancha entre «{nombre_desde}» y «{nombre_hacia}», a {} de «{nombre_desde}». Queda anotado para colgarla después.",
                            numero(rama),
                            numero(distancia)
                        ),
                    ),
            );
        }
        pendientes.push((hoja, rama, desde, hacia, distancia));
    }
    arbol.conectar(0, 1, matriz[0][1]);
    let (primeras, _) = lista_de_aristas(
        &arbol,
        &arbol.nombres,
        &|_, _| true,
        &|_, _| true,
        LIMITE_DE_LISTA,
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("arbol", format!("{} nodos", arbol.nombres.len()))
            .lista("Aristas del árbol (nodo–nodo: longitud)", primeras, &[0])
            .paso(
                "arbol.conectar(0, 1, distancias[0][1]);",
                format!(
                    "Ya no quedan más hojas por recortar. Con las dos primeras, «{}» y «{}», el árbol es una sola rama, tan larga como su distancia: {}. Ahora volvemos a colgar las hojas en orden inverso.",
                    abreviar(&nombres[0], 12),
                    abreviar(&nombres[1], 12),
                    numero(matriz[0][1])
                ),
            ),
    );
    let total = pendientes.len();
    for (orden, &(hoja, rama, desde, hacia, distancia)) in pendientes.iter().rev().enumerate() {
        let nombre = abreviar(&nombres[hoja], 12);
        let nombre_desde = abreviar(&arbol.nombres[desde], 12);
        let nombre_hacia = abreviar(&arbol.nombres[hacia], 12);
        let antes = arbol.nombres.len();
        let punto = filogenia::punto_en_el_camino(&mut arbol, desde, hacia, distancia);
        let nombre_punto = abreviar(&arbol.nombres[punto], 12);
        let nuevo = arbol.nombres.len() > antes;
        let (aristas, marcadas) = lista_de_aristas(
            &arbol,
            &arbol.nombres,
            &|_, _| true,
            &|a, b| a == punto || b == punto,
            LIMITE_DE_LISTA,
        );
        let visible = orden < 3 || mostrar(orden, total, 3, 8);
        if orden < 3 && narracion.hay_lugar() {
            let aviso = if nuevo {
                format!(
                    "Cae en medio de una rama: la partimos en dos y nace el nodo interno «{nombre_punto}»."
                )
            } else {
                format!("Cae justo sobre el nodo «{nombre_punto}», que ya existía.")
            };
            narracion.agregar(
                Pizarra::nueva()
                    .variable("pendiente.hoja", format!("{hoja} («{nombre}»)"))
                    .variable("pendiente.desde", format!("{desde} («{nombre_desde}»)"))
                    .variable("pendiente.hacia", format!("{hacia} («{nombre_hacia}»)"))
                    .variable("punto", &nombre_punto)
                    .lista("Aristas del árbol (nodo–nodo: longitud)", aristas.clone(), &marcadas)
                    .paso(
                        "let punto = punto_en_el_camino(",
                        format!(
                            "Para colgar «{nombre}» caminamos por el árbol de «{nombre_desde}» hacia «{nombre_hacia}» una distancia de {}. {aviso}",
                            numero(distancia)
                        ),
                    ),
            );
        }
        arbol.conectar(punto, hoja, rama);
        if visible {
            let (aristas, marcadas) = lista_de_aristas(
                &arbol,
                &arbol.nombres,
                &|_, _| true,
                &|a, b| (a == punto && b == hoja) || (a == hoja && b == punto),
                LIMITE_DE_LISTA,
            );
            narracion.agregar(
                Pizarra::nueva()
                    .variable("pendiente.hoja", format!("{hoja} («{nombre}»)"))
                    .variable("pendiente.rama", numero(rama))
                    .variable("punto", &nombre_punto)
                    .variable("arbol", format!("{} nodos", arbol.nombres.len()))
                    .lista("Aristas del árbol (nodo–nodo: longitud)", aristas, &marcadas)
                    .paso(
                        "arbol.conectar(punto, pendiente.hoja, pendiente.rama);",
                        format!(
                            "«{nombre}» queda colgada de «{nombre_punto}» con su rama propia de {}. El árbol ya tiene {} hojas puestas.",
                            numero(rama),
                            orden + 3
                        ),
                    ),
            );
        }
    }
    let reproduce = filogenia::reproduce_las_distancias(&arbol, &matriz);
    let mut etiquetas = arbol.clone();
    nombres_unicos(&mut etiquetas, None);
    let newick = abreviar(
        &newick_para_dibujar(&etiquetas, raiz_para_dibujar(&etiquetas)),
        60,
    );
    let ultimo = if reproduce {
        Pizarra::nueva()
            .variable("arbol", newick)
            .paso(
                "Ok(arbol)",
                format!(
                    "Comprobamos el resultado: las distancias entre hojas del árbol coinciden con tu tabla en las {} parejas. La matriz es aditiva y este es su único árbol. Mira «Qué significa» para ver cuántos nodos internos salieron.",
                    miles(cantidad * (cantidad - 1) / 2)
                ),
            )
    } else {
        let falla = match filogenia::cuarteto_no_aditivo(&matriz) {
            Some((i, j, k, l)) => format!(
                " Un grupo que lo delata es «{}», «{}», «{}» y «{}».",
                nombres[i], nombres[j], nombres[k], nombres[l]
            ),
            None => String::new(),
        };
        Pizarra::nueva()
            .variable("arbol", newick)
            .paso(
                "if !reproduce_las_distancias(&arbol, distancias) {",
                format!(
                    "Comprobamos el resultado y no coincide: alguna distancia del árbol difiere de tu tabla, o una rama salió negativa. La matriz no es aditiva y ningún árbol la reproduce exactamente.{falla} Mira «Qué significa» y prueba «Unión de vecinos»."
                ),
            )
    };
    Ok((narracion.terminar(ultimo), (arbol, reproduce)))
}

fn pasos_filogenia_aditiva(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    narrar_filogenia_aditiva(entrada).map(|(pasos, _)| pasos)
}

fn narrar_upgma(entrada: &Entrada) -> Result<(Vec<Paso>, filogenia::ArbolConRaiz), String> {
    let (nombres, matriz) = leer_matriz(entrada, "matriz")?;
    let mut completo = filogenia::upgma(&matriz, &nombres)?;
    nombres_unicos(&mut completo.arbol, None);
    let etiquetas = completo.arbol.nombres.clone();
    let cantidad = nombres.len();
    let mut narracion = Narracion::nueva(MAXIMO_DE_PASOS);
    narracion.agregar(
        con_matriz(
            Pizarra::nueva().variable("distancias", format!("{cantidad}×{cantidad}")),
            "Distancias",
            &nombres,
            &matriz,
            &[],
        )
        .paso(
            "validar_distancias(distancias)?;",
            format!(
                "Tienes las distancias entre {cantidad} especies. UPGMA arma el árbol de abajo hacia arriba: junta siempre los dos grupos más cercanos en un grupo nuevo, y repite hasta que quede uno solo. Supone que todas las especies cambian al mismo ritmo (reloj molecular)."
            ),
        ),
    );
    let mut arbol = Arbol::default();
    let mut grupos: Vec<usize> = nombres
        .iter()
        .map(|nombre| arbol.agregar_nodo(nombre))
        .collect();
    let mut hojas_por_grupo = vec![1; nombres.len()];
    let mut edades = vec![0.0; nombres.len()];
    let mut actual = matriz.clone();
    narracion.agregar(
        Pizarra::nueva()
            .variable("grupos", lista_en_texto(&etiquetas[..cantidad], 6))
            .variable("hojas_por_grupo", conteos(&hojas_por_grupo))
            .variable("edades", "todas 0")
            .variable("actual", format!("{cantidad}×{cantidad}"))
            .paso(
                "let mut grupos: Vec<usize> = nombres",
                "Al principio cada especie es su propio grupo, con una sola hoja y edad 0: todavía no ha pasado tiempo desde que cada una existe.",
            ),
    );
    let rondas = cantidad - 1;
    let mut ronda = 0;
    while actual.len() > 1 {
        let en_detalle = ronda < 2;
        let visible = en_detalle || mostrar(ronda, rondas, 2, 14);
        let (i, j) = filogenia::par_mas_cercano(&actual);
        let etiquetas_actuales: Vec<String> = grupos
            .iter()
            .map(|&grupo| etiquetas[grupo].clone())
            .collect();
        let (nombre_i, nombre_j) = (
            abreviar(&etiquetas_actuales[i], 12),
            abreviar(&etiquetas_actuales[j], 12),
        );
        let distancia = actual[i][j];
        if visible && narracion.hay_lugar() {
            narracion.agregar(
                con_matriz(
                    Pizarra::nueva()
                        .variable("grupos", lista_en_texto(&etiquetas_actuales, 6))
                        .variable("actual", format!("{0}×{0}", actual.len()))
                        .variable("i", format!("{i} («{nombre_i}»)"))
                        .variable("j", format!("{j} («{nombre_j}»)")),
                    "Distancias entre grupos",
                    &etiquetas_actuales,
                    &actual,
                    &[(i, j), (j, i)],
                )
                .paso(
                    "let (i, j) = par_mas_cercano(&actual);",
                    format!(
                        "Ronda {}: los dos grupos más cercanos son «{nombre_i}» y «{nombre_j}», a {}. Los vamos a juntar.",
                        ronda + 1,
                        numero(distancia)
                    ),
                ),
            );
        }
        let nuevo = arbol.agregar_nodo(&filogenia::nombre_interno(arbol.nombres.len()));
        let nombre_nuevo = abreviar(&etiquetas[nuevo], 12);
        let edad = actual[i][j] / 2.0;
        edades.push(edad);
        if en_detalle && narracion.hay_lugar() {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("nuevo", &nombre_nuevo)
                    .variable("edad", numero(edad))
                    .paso(
                        "let edad = actual[i][j] / 2.0;",
                        format!(
                            "El grupo nuevo «{nombre_nuevo}» nace a la edad {} / 2 = {}: si los dos cambiaron al mismo ritmo desde su antepasado, cada uno acumuló la mitad de la distancia.",
                            numero(distancia),
                            numero(edad)
                        ),
                    ),
            );
        }
        let rama_i = edad - edades[grupos[i]];
        let rama_j = edad - edades[grupos[j]];
        arbol.conectar(nuevo, grupos[i], rama_i);
        arbol.conectar(nuevo, grupos[j], rama_j);
        if visible && narracion.hay_lugar() {
            let elementos = vec![
                format!("{nombre_nuevo}–{nombre_i}: {}", numero(rama_i)),
                format!("{nombre_nuevo}–{nombre_j}: {}", numero(rama_j)),
            ];
            narracion.agregar(
                Pizarra::nueva()
                    .variable("nuevo", &nombre_nuevo)
                    .variable("edad", numero(edad))
                    .variable("edades", format!("{} grupos con edad", edades.len()))
                    .lista("Ramas nuevas", elementos, &[0, 1])
                    .paso(
                        "arbol.conectar(nuevo, grupos[i], edad - edades[grupos[i]]);",
                        format!(
                            "Cada rama mide la edad del grupo nuevo menos la edad del grupo que cuelga de ella: «{nombre_i}» queda a {} − {} = {} y «{nombre_j}» a {} − {} = {}.",
                            numero(edad),
                            numero(edades[grupos[i]]),
                            numero(rama_i),
                            numero(edad),
                            numero(edades[grupos[j]]),
                            numero(rama_j)
                        ),
                    ),
            );
        }
        let fila_nueva = filogenia::promedios_hacia_el_grupo(&actual, &hojas_por_grupo, i, j);
        if en_detalle && narracion.hay_lugar() && !fila_nueva.is_empty() {
            let otro = (0..actual.len()).find(|&k| k != i && k != j).unwrap_or(0);
            let nombre_otro = abreviar(&etiquetas_actuales[otro], 12);
            let elementos: Vec<String> = (0..actual.len())
                .filter(|&k| k != i && k != j)
                .zip(&fila_nueva)
                .take(LIMITE_DE_LISTA)
                .map(|(k, &valor)| {
                    format!(
                        "«{nombre_nuevo}» a «{}»: {}",
                        abreviar(&etiquetas_actuales[k], 10),
                        numero(valor)
                    )
                })
                .collect();
            narracion.agregar(
                Pizarra::nueva()
                    .variable("hojas_por_grupo", conteos(&hojas_por_grupo))
                    .variable("fila_nueva", lista_en_texto(&elementos, 4))
                    .lista("Distancias del grupo nuevo", elementos.clone(), &[0])
                    .paso(
                        "let fila_nueva = promedios_hacia_el_grupo(&actual, &hojas_por_grupo, i, j);",
                        format!(
                            "La distancia del grupo nuevo a cada otro es el promedio de las distancias de sus dos partes, contando cuántas hojas tiene cada una. Por ejemplo, hasta «{nombre_otro}»: ({} × {} + {} × {}) / ({} + {}) = {}.",
                            numero(actual[otro][i]),
                            hojas_por_grupo[i],
                            numero(actual[otro][j]),
                            hojas_por_grupo[j],
                            hojas_por_grupo[i],
                            hojas_por_grupo[j],
                            numero(fila_nueva[0])
                        ),
                    ),
            );
        }
        filogenia::reemplazar_par(&mut actual, i, j, fila_nueva);
        let hojas_del_nuevo = hojas_por_grupo[i] + hojas_por_grupo[j];
        filogenia::quitar_par(&mut hojas_por_grupo, i, j);
        hojas_por_grupo.push(hojas_del_nuevo);
        filogenia::quitar_par(&mut grupos, i, j);
        grupos.push(nuevo);
        if visible && actual.len() > 1 && narracion.hay_lugar() {
            let etiquetas_nuevas: Vec<String> = grupos
                .iter()
                .map(|&grupo| etiquetas[grupo].clone())
                .collect();
            let ultimo = etiquetas_nuevas.len() - 1;
            narracion.agregar(
                con_matriz(
                    Pizarra::nueva()
                        .variable("grupos", lista_en_texto(&etiquetas_nuevas, 6))
                        .variable("hojas_por_grupo", conteos(&hojas_por_grupo))
                        .variable("actual", format!("{0}×{0}", actual.len())),
                    "Distancias entre grupos",
                    &etiquetas_nuevas,
                    &actual,
                    &(0..=ultimo).flat_map(|k| [(ultimo, k), (k, ultimo)]).collect::<Vec<_>>(),
                )
                .paso(
                    "reemplazar_par(&mut actual, i, j, fila_nueva);",
                    format!(
                        "«{nombre_i}» y «{nombre_j}» salen de la tabla y entra «{nombre_nuevo}», con {hojas_del_nuevo} hojas. La tabla queda de {0}×{0}.",
                        actual.len()
                    ),
                ),
            );
        }
        ronda += 1;
    }
    let raiz = grupos[0];
    let altura = edades[raiz];
    let newick = abreviar(&newick_para_dibujar(&completo.arbol, completo.raiz), 60);
    let ultimo = Pizarra::nueva()
        .variable("grupos", format!("[{}]", abreviar(&etiquetas[raiz], 12)))
        .variable("edades", format!("raíz a {}", numero(altura)))
        .variable("arbol", newick)
        .paso(
            "raiz: grupos[0],",
            format!(
                "Solo queda un grupo, «{}»: es la raíz del árbol, a la edad {}. Como cada unión puso a los dos grupos a la misma altura, todas las hojas quedaron a igual distancia de la raíz. Mira «Qué significa» para ver si tus datos cumplen ese reloj.",
                abreviar(&etiquetas[raiz], 12),
                numero(altura)
            ),
        );
    let resultado = filogenia::ArbolConRaiz { arbol, raiz };
    Ok((narracion.terminar(ultimo), resultado))
}

fn pasos_upgma(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    narrar_upgma(entrada).map(|(pasos, _)| pasos)
}

fn marcas_de_diferencias(
    primera: &[u8],
    segunda: &[u8],
    hasta: usize,
) -> Vec<(usize, usize, Marca)> {
    primera
        .iter()
        .zip(segunda)
        .take(hasta)
        .enumerate()
        .filter(|&(_, (&a, &b))| es_letra(a) && es_letra(b) && a != b)
        .map(|(posicion, _)| (posicion, posicion + 1, Marca::Hallado))
        .collect()
}

fn narrar_un_par(
    narracion: &mut Narracion,
    nombres: &[String],
    secuencias: &[Vec<u8>],
    distancia_p: f64,
    detalle: usize,
) {
    let (a, b) = (&secuencias[0], &secuencias[1]);
    let (nombre_a, nombre_b) = (abreviar(&nombres[0], 12), abreviar(&nombres[1], 12));
    let largo = a.len();
    let mut comparadas = 0usize;
    let mut diferentes = 0usize;
    narracion.agregar(
        Pizarra::nueva()
            .variable("primera", format!("«{nombre_a}»"))
            .variable("segunda", format!("«{nombre_b}»"))
            .secuencia(a, &[])
            .paso(
                ".map(|j| distancia_p(&secuencias[i], &secuencias[j]))",
                format!(
                    "Empezamos con el primer par, «{nombre_a}» contra «{nombre_b}»: las ponemos una sobre otra y vamos comparando columna por columna, {} en total.",
                    miles(largo)
                ),
            ),
    );
    for posicion in 0..largo.min(detalle) {
        let (letra_a, letra_b) = (a[posicion], b[posicion]);
        let valida = es_letra(letra_a) && es_letra(letra_b);
        if valida {
            comparadas += 1;
            if letra_a != letra_b {
                diferentes += 1;
            }
        }
        let mut marcas = marcas_de_diferencias(a, b, posicion);
        marcas.push((posicion, posicion + 1, Marca::Actual));
        let (linea, nota) = if !valida {
            (
                "if es_base(a) && es_base(b) {",
                format!(
                    "Posición {posicion}: hay un hueco («-») en una de las dos. Esa columna no se puede comparar, así que no cuenta."
                ),
            )
        } else if letra_a == letra_b {
            (
                "comparadas += 1;",
                format!(
                    "Posición {posicion}: las dos tienen {}. Cuenta como comparada y no como distinta.",
                    letra_a as char
                ),
            )
        } else {
            (
                "diferentes += 1;",
                format!(
                    "Posición {posicion}: «{nombre_a}» tiene {} y «{nombre_b}» tiene {}. Cuenta como comparada y como distinta.",
                    letra_a as char, letra_b as char
                ),
            )
        };
        narracion.agregar(
            Pizarra::nueva()
                .variable("a", letra_a as char)
                .variable("b", letra_b as char)
                .variable("comparadas", comparadas)
                .variable("diferentes", diferentes)
                .secuencia(a, &marcas)
                .lista(
                    "Letras de esta columna",
                    vec![
                        format!("{nombre_a}: {}", letra_a as char),
                        format!("{nombre_b}: {}", letra_b as char),
                    ],
                    &[],
                )
                .paso(linea, nota),
        );
    }
    let (total_comparadas, total_diferentes) = filogenia::contar_diferencias(a, b);
    if largo > detalle {
        narracion.agregar(
            Pizarra::nueva()
                .variable("comparadas", total_comparadas)
                .variable("diferentes", total_diferentes)
                .secuencia(a, &marcas_de_diferencias(a, b, largo))
                .paso(
                    "for (&a, &b) in primera.iter().zip(segunda) {",
                    format!(
                        "Así sigue hasta el final de las {} columnas. Resultado: {} comparadas y {} distintas (las letras marcadas).",
                        miles(largo),
                        miles(total_comparadas),
                        miles(total_diferentes)
                    ),
                ),
        );
    }
    narracion.agregar(
        Pizarra::nueva()
            .variable("comparadas", total_comparadas)
            .variable("diferentes", total_diferentes)
            .variable("p", numero(distancia_p))
            .paso(
                "Some(diferentes as f64 / comparadas as f64)",
                format!(
                    "La distancia p es la fracción de columnas distintas: {} / {} = {}, es decir {}.",
                    miles(total_diferentes),
                    miles(total_comparadas),
                    numero(distancia_p),
                    porcentaje(distancia_p, 1)
                ),
            ),
    );
    let corregida = filogenia::distancia_jukes_cantor(distancia_p);
    let (linea, nota) = if distancia_p >= 0.75 {
        (
            "if p >= 0.75 {",
            "Con p de 75 % o más, dos secuencias sin parentesco ya coinciden por azar en 1 de cada 4 columnas: la fórmula no tiene solución y la distancia es infinita (∞).".to_string(),
        )
    } else if distancia_p <= 0.0 {
        (
            "if p <= 0.0 {",
            "Sin diferencias no hay nada que corregir: la distancia es 0.".to_string(),
        )
    } else {
        (
            "-0.75 * (1.0 - 4.0 * p / 3.0).ln()",
            format!(
                "Jukes–Cantor corrige p por los cambios ocultos (una misma columna pudo cambiar dos veces): d = −0,75 · ln(1 − 4p/3) = {}. Siempre sale un poco mayor que p.",
                numero(corregida)
            ),
        )
    };
    narracion.agregar(
        Pizarra::nueva()
            .variable("p", numero(distancia_p))
            .variable("d", numero(corregida))
            .paso(linea, nota),
    );
}

fn narrar_distancias_evolutivas(entrada: &Entrada) -> Result<(Vec<Paso>, Vec<Vec<f64>>), String> {
    let (nombres, secuencias) = leer_secuencias(entrada, "secuencias", true)?;
    if nombres.len() < 2 {
        return Err("Hacen falta al menos dos secuencias.".to_string());
    }
    let largo = revisar_mismo_largo(&nombres, &secuencias)?;
    let distancias_p = filogenia::matriz_de_distancias_p(&secuencias)?;
    let jukes_cantor = filogenia::corregir_con_jukes_cantor(&distancias_p);
    let cantidad = nombres.len();
    let mut narracion = Narracion::nueva(MAXIMO_DE_PASOS);
    narracion.agregar(
        Pizarra::nueva()
            .variable("secuencias", cantidad)
            .secuencia(&secuencias[0], &[])
            .paso(
                "pub fn matriz_jukes_cantor(",
                format!(
                    "Tienes {cantidad} secuencias ya alineadas, de {} columnas cada una. Queremos una tabla con la distancia evolutiva entre cada par: primero contamos cuánto se diferencian y después corregimos ese número.",
                    miles(largo)
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("longitud", miles(largo))
            .paso(
                "let longitud = primera.len();",
                format!(
                    "Antes de comparar se comprueba que todas midan lo mismo ({} columnas). Si no, las columnas no se corresponderían y la comparación no tendría sentido.",
                    miles(largo)
                ),
            ),
    );
    narrar_un_par(&mut narracion, &nombres, &secuencias, distancias_p[0][1], 8);
    if cantidad > 2 {
        let mut filas = Vec::new();
        'pares: for i in 0..cantidad {
            for j in (i + 1)..cantidad {
                if filas.len() == LIMITE_DE_LISTA {
                    break 'pares;
                }
                let (comparadas, diferentes) =
                    filogenia::contar_diferencias(&secuencias[i], &secuencias[j]);
                filas.push(vec![
                    abreviar(&nombres[i], 9),
                    abreviar(&nombres[j], 9),
                    miles(comparadas),
                    miles(diferentes),
                    numero(distancias_p[i][j]),
                ]);
            }
        }
        let pares = cantidad * (cantidad - 1) / 2;
        narracion.agregar(
            Pizarra::nueva()
                .variable("distancias_p", format!("{cantidad}×{cantidad}"))
                .tabla(
                    "Pares (los primeros)",
                    &["Una", "Otra", "Comparadas", "Distintas", "p"],
                    filas,
                    &[(0, 4)],
                )
                .paso(
                    "matriz[i][j] = p;",
                    format!(
                        "Se hace lo mismo con los otros {} pares, repartidos entre los núcleos de tu computadora. Cada p se guarda en la tabla dos veces, en [i][j] y en [j][i].",
                        miles(pares - 1)
                    ),
                ),
        );
    }
    let saturados = jukes_cantor
        .iter()
        .flatten()
        .filter(|valor| valor.is_infinite())
        .count()
        / 2;
    let aviso = if saturados > 0 {
        format!(
            " {} par(es) salieron con ∞: difieren en 75 % o más de las columnas y no se puede medir su distancia así.",
            miles(saturados)
        )
    } else {
        String::new()
    };
    narracion.agregar(
        con_matriz(
            Pizarra::nueva().variable("distancias_p", format!("{cantidad}×{cantidad}")),
            "Distancia de Jukes–Cantor",
            &nombres,
            &jukes_cantor,
            &[(0, 1), (1, 0)],
        )
        .paso(
            ".map(|fila| fila.iter().map(|&p| distancia_jukes_cantor(p)).collect())",
            format!("Aplicamos la misma corrección a cada p de la tabla y sale la tabla corregida.{aviso}"),
        ),
    );
    let ultimo = Pizarra::nueva()
        .variable("distancias_p", format!("{cantidad}×{cantidad}"))
        .paso(
            "Ok(corregir_con_jukes_cantor(&distancias_p))",
            "Eso es todo: dos tablas, p (lo que se ve) y Jukes–Cantor (lo que se estima que pasó, en sustituciones por columna). Mira «Qué significa» para ver el par más distinto y cuánto añadió la corrección.",
        );
    Ok((narracion.terminar(ultimo), jukes_cantor))
}

fn pasos_distancias_evolutivas(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    narrar_distancias_evolutivas(entrada).map(|(pasos, _)| pasos)
}

fn narrar_union_de_vecinos(
    narracion: &mut Narracion,
    matriz: &[Vec<f64>],
    nombres: &[String],
) -> Result<Arbol, String> {
    let mut completo = filogenia::union_de_vecinos(matriz, nombres)?;
    nombres_unicos(&mut completo, None);
    let etiquetas = completo.nombres.clone();
    let mut estado = filogenia::empezar_union(matriz, nombres)?;
    let cantidad = nombres.len();
    narracion.agregar(
        con_matriz(
            Pizarra::nueva().variable("estado", format!("{cantidad} nodos")),
            "Distancias entre nodos",
            nombres,
            matriz,
            &[],
        )
        .paso(
            "let mut estado = empezar_union(distancias, nombres)?;",
            format!(
                "Ahora armamos el árbol por unión de vecinos. Cada especie es un nodo suelto. En cada ronda buscamos dos nodos que sean vecinos en el árbol, los unimos a un nodo nuevo y los reemplazamos en la tabla por ese nodo. Harán falta {} rondas.",
                cantidad.saturating_sub(2)
            ),
        ),
    );
    let rondas = cantidad.saturating_sub(2);
    let mut ronda = 0;
    while estado.distancias.len() > 2 {
        let en_detalle = ronda < 2;
        let visible = (en_detalle || mostrar(ronda, rondas, 2, 12)) && narracion.hay_lugar();
        let antes = estado.clone();
        let n = antes.distancias.len();
        let Some(par) = filogenia::unir_par_de_vecinos(&mut estado) else {
            break;
        };
        if visible {
            let totales = filogenia::distancias_totales(&antes.distancias);
            let estrella = filogenia::matriz_de_union(&antes.distancias);
            let nombres_antes: Vec<String> = antes
                .nodos
                .iter()
                .map(|&nodo| etiquetas[nodo].clone())
                .collect();
            let i = antes
                .nodos
                .iter()
                .position(|&nodo| nodo == par.primero)
                .unwrap_or(0);
            let j = antes
                .nodos
                .iter()
                .position(|&nodo| nodo == par.segundo)
                .unwrap_or(0);
            let primero = abreviar(&etiquetas[par.primero], 12);
            let segundo = abreviar(&etiquetas[par.segundo], 12);
            let nuevo = abreviar(&etiquetas[par.nuevo], 12);
            let delta = (totales[i] - totales[j]) / (n - 2) as f64;
            if en_detalle {
                let elementos: Vec<String> = nombres_antes
                    .iter()
                    .zip(&totales)
                    .take(LIMITE_DE_LISTA)
                    .map(|(nombre, &total)| format!("{}: {}", abreviar(nombre, 10), numero(total)))
                    .collect();
                narracion.agregar(
                    Pizarra::nueva()
                        .variable("n", n)
                        .variable("totales", lista_en_texto(&elementos, 4))
                        .lista("Distancia total de cada nodo a todos los demás", elementos, &[])
                        .paso(
                            "let totales = distancias_totales(&estado.distancias);",
                            format!(
                                "Ronda {}: sumamos cada fila de la tabla. Es la distancia total de cada uno de los {n} nodos a todos los demás; sirve para no unir por error a dos nodos que solo parecen cercanos porque todo está cerca de ellos.",
                                ronda + 1
                            ),
                        ),
                );
            }
            narracion.agregar(
                con_matriz(
                    Pizarra::nueva()
                        .variable("n", n)
                        .variable("i", format!("{i} («{primero}»)"))
                        .variable("j", format!("{j} («{segundo}»)")),
                    "Matriz D* (el menor valor marca a los vecinos)",
                    &nombres_antes,
                    &estrella,
                    &[(i, j), (j, i)],
                )
                .paso(
                    "let (i, j) = par_con_menor_valor_de_union(&estado.distancias, &totales);",
                    format!(
                        "Ronda {}: para cada par calculamos D*(i,j) = (n−2)·D(i,j) − total(i) − total(j). El menor es {}, entre «{primero}» y «{segundo}»: son vecinos y los unimos.",
                        ronda + 1,
                        numero(estrella[i][j])
                    ),
                ),
            );
            if en_detalle {
                narracion.agregar(
                    Pizarra::nueva()
                        .variable("delta", numero(delta))
                        .variable("rama_i", numero(par.rama_primero))
                        .variable("nuevo", &nuevo)
                        .lista(
                            "Ramas nuevas",
                            vec![
                                format!("{nuevo}–{primero}: {}", numero(par.rama_primero)),
                                format!("{nuevo}–{segundo}: {}", numero(par.rama_segundo)),
                            ],
                            &[0, 1],
                        )
                        .paso(
                            "let rama_i = (estado.distancias[i][j] + delta) / 2.0;",
                            format!(
                                "Δ = (total(i) − total(j)) / (n−2) = {}. La rama de «{primero}» al nodo nuevo «{nuevo}» mide (D(i,j) + Δ) / 2 = {} y la de «{segundo}» (D(i,j) − Δ) / 2 = {}.",
                                numero(delta),
                                numero(par.rama_primero),
                                numero(par.rama_segundo)
                            ),
                        ),
                );
            }
            let nombres_despues: Vec<String> = estado
                .nodos
                .iter()
                .map(|&nodo| etiquetas[nodo].clone())
                .collect();
            let ultimo = nombres_despues.len() - 1;
            let del_nuevo: Vec<(usize, usize)> = (0..=ultimo)
                .flat_map(|k| [(ultimo, k), (k, ultimo)])
                .collect();
            narracion.agregar(
                con_matriz(
                    Pizarra::nueva()
                        .variable("nuevo", &nuevo)
                        .variable("estado", format!("{} nodos", estado.nodos.len())),
                    "Distancias entre nodos",
                    &nombres_despues,
                    &estado.distancias,
                    &del_nuevo,
                )
                .paso(
                    "reemplazar_par(&mut estado.distancias, i, j, fila_nueva);",
                    format!(
                        "«{primero}» y «{segundo}» salen de la tabla y entra «{nuevo}», cuya distancia a cada otro nodo k es (D(k,i) + D(k,j) − D(i,j)) / 2. Ahora la tabla es de {0}×{0}.",
                        estado.distancias.len()
                    ),
                ),
            );
        }
        ronda += 1;
    }
    if let Some((a, b, longitud)) = filogenia::unir_ultimo_par(&mut estado) {
        let (aristas, marcadas) = lista_de_aristas(
            &estado.arbol,
            &etiquetas,
            &|_, _| true,
            &|u, v| (u == a && v == b) || (u == b && v == a),
            LIMITE_DE_LISTA,
        );
        narracion.agregar(
            Pizarra::nueva()
                .variable("longitud", numero(longitud))
                .lista("Aristas del árbol (nodo–nodo: longitud)", aristas, &marcadas)
                .paso(
                    ".conectar(estado.nodos[0], estado.nodos[1], longitud);",
                    format!(
                        "Quedan dos nodos, «{}» y «{}»: se unen con una rama de {} y el árbol está completo.",
                        abreviar(&etiquetas[a], 12),
                        abreviar(&etiquetas[b], 12),
                        numero(longitud)
                    ),
                ),
        );
    }
    Ok(estado.arbol)
}

fn narrar_arbol_desde_secuencias(entrada: &Entrada) -> Result<(Vec<Paso>, Option<Arbol>), String> {
    let (nombres, secuencias) = leer_secuencias(entrada, "secuencias", true)?;
    if nombres.len() < 2 {
        return Err("Hacen falta al menos dos secuencias.".to_string());
    }
    let largo = revisar_mismo_largo(&nombres, &secuencias)?;
    let distancias_p = filogenia::matriz_de_distancias_p(&secuencias)?;
    let jukes_cantor = filogenia::corregir_con_jukes_cantor(&distancias_p);
    let cantidad = nombres.len();
    let mut narracion = Narracion::nueva(MAXIMO_DE_PASOS);
    narracion.agregar(
        Pizarra::nueva()
            .variable("secuencias", cantidad)
            .secuencia(&secuencias[0], &[])
            .paso(
                "let distancias = matriz_jukes_cantor(secuencias)?;",
                format!(
                    "Tienes {cantidad} secuencias alineadas de {} columnas. Hacer el árbol toma dos etapas: primero convertimos las secuencias en una tabla de distancias (Jukes–Cantor) y después, con esa tabla, armamos el árbol por unión de vecinos.",
                    miles(largo)
                ),
            ),
    );
    narrar_un_par(&mut narracion, &nombres, &secuencias, distancias_p[0][1], 6);
    narracion.agregar(
        con_matriz(
            Pizarra::nueva().variable("distancias_p", format!("{cantidad}×{cantidad}")),
            "Distancia p",
            &nombres,
            &distancias_p,
            &[],
        )
        .paso(
            "let distancias_p = matriz_de_distancias_p(secuencias)?;",
            "Lo mismo se hace con todos los pares, y sale esta tabla de fracciones de columnas distintas (p).",
        ),
    );
    narracion.agregar(
        con_matriz(
            Pizarra::nueva().variable("distancias", format!("{cantidad}×{cantidad}")),
            "Distancia de Jukes–Cantor",
            &nombres,
            &jukes_cantor,
            &[],
        )
        .paso(
            "Ok(corregir_con_jukes_cantor(&distancias_p))",
            "Y cada p se corrige con Jukes–Cantor para estimar los cambios ocultos. Esta tabla es la que usa el árbol.",
        ),
    );
    if let Some((i, j)) = filogenia::par_saturado(&jukes_cantor) {
        let ultimo = Pizarra::nueva()
            .variable("distancias", format!("{cantidad}×{cantidad}"))
            .paso(
                "if let Some((i, j)) = par_saturado(&distancias) {",
                format!(
                    "Antes de armar el árbol se busca alguna distancia infinita, y aparece una: «{}» y «{}» difieren en 75 % o más de las columnas. A ese nivel no se distingue el parentesco del azar, así que el árbol no se puede construir con estos datos. Mira «Qué significa»: quita las secuencias más lejanas o usa un gen más conservado.",
                    abreviar(&nombres[i], 12),
                    abreviar(&nombres[j], 12)
                ),
            );
        return Ok((narracion.terminar(ultimo), None));
    }
    narracion.agregar(
        Pizarra::nueva()
            .variable("distancias", format!("{cantidad}×{cantidad}"))
            .paso(
                "if let Some((i, j)) = par_saturado(&distancias) {",
                "Se revisa que ninguna distancia sea infinita: todas son finitas, así que el árbol se puede construir.",
            ),
    );
    let arbol = narrar_union_de_vecinos(&mut narracion, &jukes_cantor, &nombres)?;
    let mut con_nombres = arbol.clone();
    nombres_unicos(&mut con_nombres, None);
    let newick = abreviar(
        &newick_para_dibujar(&con_nombres, raiz_para_dibujar(&con_nombres)),
        60,
    );
    let ultimo = Pizarra::nueva()
        .variable("distancias", format!("{cantidad}×{cantidad}"))
        .variable("arbol", newick)
        .paso(
            "union_de_vecinos(&distancias, nombres)",
            "Listo: ese es el árbol, con ramas medidas en sustituciones por columna. Mira «Qué significa» para ver el par más parecido y la rama más larga.",
        );
    Ok((narracion.terminar(ultimo), Some(arbol)))
}

fn pasos_arbol_desde_secuencias(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    narrar_arbol_desde_secuencias(entrada).map(|(pasos, _)| pasos)
}

fn letra_mas_barata_de(costos: &[u32; 4]) -> usize {
    let mut mejor = 0;
    for letra in 1..4 {
        if costos[letra] < costos[mejor] {
            mejor = letra;
        }
    }
    mejor
}

fn letra_para_el_hijo_de(costos_del_hijo: &[u32; 4], letra_del_padre: usize) -> usize {
    let mas_barata = letra_mas_barata_de(costos_del_hijo);
    if costos_del_hijo[letra_del_padre] <= costos_del_hijo[mas_barata] + 1 {
        letra_del_padre
    } else {
        mas_barata
    }
}

fn costos_de_nodo(
    arbol: &Arbol,
    costos: &[[u32; 4]],
    padres: &[Option<usize>],
    secuencias: &[Vec<u8>],
    nodo: usize,
    posicion: usize,
) -> [u32; 4] {
    let mut suma = match secuencias[nodo]
        .get(posicion)
        .and_then(|&letra| nucleo::secuencia::indice_base(letra))
    {
        Some(propia) => {
            let mut fija = [COSTO_IMPOSIBLE; 4];
            fija[propia] = 0;
            fija
        }
        None => [0; 4],
    };
    for &(hijo, _) in &arbol.vecinos[nodo] {
        if padres[nodo] == Some(hijo) {
            continue;
        }
        let del_hijo = costos[hijo];
        let mas_barato = del_hijo[letra_mas_barata_de(&del_hijo)];
        for (letra, total) in suma.iter_mut().enumerate() {
            *total += del_hijo[letra].min(mas_barato + 1);
        }
    }
    suma
}

fn costo_en_texto(costo: u32) -> String {
    if costo >= COSTO_IMPOSIBLE {
        "∞".to_string()
    } else {
        costo.to_string()
    }
}

fn letra_de(indice: usize) -> char {
    b"ACGT"[indice] as char
}

fn tabla_de_costos(
    pizarra: Pizarra,
    etiquetas: &[String],
    costos: &[[u32; 4]],
    nodos: &[usize],
    resaltado: usize,
) -> Pizarra {
    let filas: Vec<Vec<String>> = nodos
        .iter()
        .map(|&nodo| {
            let mut fila = vec![abreviar(&etiquetas[nodo], 10)];
            fila.extend(costos[nodo].iter().map(|&costo| costo_en_texto(costo)));
            fila
        })
        .collect();
    let mejor = letra_mas_barata_de(&costos[nodos[resaltado]]);
    pizarra.tabla(
        "Costo mínimo de cada letra en cada nodo",
        &["Nodo", "A", "C", "G", "T"],
        filas,
        &[(resaltado, mejor + 1)],
    )
}

fn hijos_de(arbol: &Arbol, padres: &[Option<usize>], nodo: usize) -> Vec<usize> {
    arbol.vecinos[nodo]
        .iter()
        .map(|&(hijo, _)| hijo)
        .filter(|&hijo| padres[nodo] != Some(hijo))
        .collect()
}

type ResultadoDeParsimonia = (usize, Vec<Vec<u8>>, Vec<usize>);

fn narrar_parsimonia(
    narracion: &mut Narracion,
    arbol: &Arbol,
    raiz: usize,
    por_nodo: &[Vec<u8>],
    etiquetas: &[String],
) -> ResultadoDeParsimonia {
    let cantidad = arbol.nombres.len();
    let hojas: Vec<usize> = (0..cantidad)
        .filter(|&nodo| !por_nodo[nodo].is_empty())
        .collect();
    let Some(&primera) = hojas.first() else {
        return (0, Vec::new(), Vec::new());
    };
    let longitud = por_nodo[primera].len();
    let nombre_raiz = abreviar(&etiquetas[raiz], 12);
    let elementos: Vec<String> = hojas
        .iter()
        .take(LIMITE_DE_LISTA)
        .map(|&nodo| {
            format!(
                "{}: {}",
                abreviar(&etiquetas[nodo], 10),
                abreviar(&String::from_utf8_lossy(&por_nodo[nodo]), 24)
            )
        })
        .collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("longitud", miles(longitud))
            .variable("cantidad", cantidad)
            .secuencia(&por_nodo[primera], &[])
            .lista("Secuencias de las hojas", elementos, &[])
            .paso(
                "let longitud = revisar_hojas_con_secuencia(arbol, raiz, secuencias)?;",
                format!(
                    "Tienes un árbol de {cantidad} nodos con {} hojas, cada una con su secuencia de {} letras. Queremos poner letras en los nodos internos (los antepasados) de modo que hagan falta las menos mutaciones posibles. Cada posición se resuelve por separado.",
                    hojas.len(),
                    miles(longitud)
                ),
            ),
    );
    let (orden, padres) = filogenia::recorrer_desde(arbol, raiz);
    let nombres_orden: Vec<String> = orden.iter().map(|&nodo| etiquetas[nodo].clone()).collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("raiz", &nombre_raiz)
            .variable("orden", lista_en_texto(&nombres_orden, 6))
            .lista(
                "Orden de recorrido (de la raíz hacia las hojas)",
                nombres_orden.iter().take(LIMITE_DE_LISTA).cloned().collect(),
                &[0],
            )
            .paso(
                "let (orden, padres) = recorrer_desde(arbol, raiz);",
                format!(
                    "Recorremos el árbol desde la raíz «{nombre_raiz}» y anotamos el orden y el padre de cada nodo. Para calcular iremos al revés: de las hojas hacia la raíz."
                ),
            ),
    );
    let mut costos = vec![[0u32; 4]; cantidad];
    let mut letras = vec![0usize; cantidad];
    let mut inferidas = vec![vec![b'A'; longitud]; cantidad];
    let mut cambios_por_posicion: Vec<usize> = Vec::with_capacity(longitud);
    let mut compactas = 0usize;
    for posicion in 0..longitud {
        let marcas = [(posicion, posicion + 1, Marca::Actual)];
        let letras_de_hojas: Vec<String> = hojas
            .iter()
            .take(LIMITE_DE_LISTA)
            .map(|&nodo| {
                format!(
                    "{}: {}",
                    abreviar(&etiquetas[nodo], 10),
                    por_nodo[nodo][posicion] as char
                )
            })
            .collect();
        let detalle = posicion == 0;
        if detalle {
            for &hoja in hojas.iter().take(LIMITE_DE_TABLA) {
                costos[hoja] = costos_de_nodo(arbol, &costos, &padres, por_nodo, hoja, posicion);
            }
            let visibles: Vec<usize> = hojas.iter().copied().take(6).collect();
            narracion.agregar(
                tabla_de_costos(
                    Pizarra::nueva().variable("posicion", posicion),
                    etiquetas,
                    &costos,
                    &visibles,
                    0,
                )
                .secuencia(&por_nodo[primera], &marcas)
                .paso(
                    "fija[propia] = 0;",
                    "Posición 0. Cada hoja conoce su letra: costo 0 para ella e ∞ (imposible) para las otras tres, porque la hoja no puede ser otra cosa.",
                ),
            );
        }
        let mut internos_mostrados = 0usize;
        for &nodo in orden.iter().rev() {
            costos[nodo] = costos_de_nodo(arbol, &costos, &padres, por_nodo, nodo, posicion);
            let hijos = hijos_de(arbol, &padres, nodo);
            if detalle && !hijos.is_empty() && internos_mostrados < 4 && narracion.hay_lugar() {
                internos_mostrados += 1;
                let mut filas_de_nodos: Vec<usize> = hijos.iter().copied().take(4).collect();
                filas_de_nodos.push(nodo);
                let propio = filas_de_nodos.len() - 1;
                let nombre = abreviar(&etiquetas[nodo], 12);
                let barato = letra_mas_barata_de(&costos[nodo]);
                let nombres_hijos: Vec<String> = hijos
                    .iter()
                    .take(4)
                    .map(|&hijo| format!("«{}»", abreviar(&etiquetas[hijo], 12)))
                    .collect();
                narracion.agregar(
                    tabla_de_costos(
                        Pizarra::nueva()
                            .variable("posicion", posicion)
                            .variable("nodo", &nombre),
                        etiquetas,
                        &costos,
                        &filas_de_nodos,
                        propio,
                    )
                    .paso(
                        "costos[nodo] = costos_del_nodo(arbol, &costos, &padres, secuencias, nodo, posicion);",
                        format!(
                            "«{nombre}» une a {}. Para cada letra suma lo que cuesta cada hijo: lo que ya cuesta con esa letra, o lo mejor que puede hacer el hijo más 1 mutación. La más barata es {} con costo {}.",
                            nombres_hijos.join(" y "),
                            letra_de(barato),
                            costos[nodo][barato]
                        ),
                    ),
                );
            }
        }
        letras[raiz] = letra_mas_barata_de(&costos[raiz]);
        cambios_por_posicion.push(costos[raiz][letras[raiz]] as usize);
        let cambios = cambios_por_posicion[posicion];
        if detalle && narracion.hay_lugar() {
            narracion.agregar(
                tabla_de_costos(
                    Pizarra::nueva()
                        .variable("posicion", posicion)
                        .variable("letras[raiz]", letra_de(letras[raiz])),
                    etiquetas,
                    &costos,
                    &[raiz],
                    0,
                )
                .paso(
                    "letras[raiz] = letra_mas_barata(&costos[raiz]);",
                    format!(
                        "Llegamos a la raíz «{nombre_raiz}». Elegimos su letra más barata: {} (costo {}). Esa es la letra del antepasado común en esta posición.",
                        letra_de(letras[raiz]),
                        cambios
                    ),
                ),
            );
            narracion.agregar(
                Pizarra::nueva()
                    .variable("posicion", posicion)
                    .variable("letras[raiz]", letra_de(letras[raiz]))
                    .variable("cambios_por_posicion", conteos(&cambios_por_posicion))
                    .paso(
                        "cambios_por_posicion.push(costos[raiz].por_letra[letras[raiz]] as usize);",
                        format!(
                            "El costo de la raíz con su letra es el número de mutaciones que hacen falta en esta posición: {cambios}. Lo guardamos."
                        ),
                    ),
            );
        }
        let mut bajados = 0usize;
        for &nodo in &orden {
            if let Some(padre) = padres[nodo] {
                letras[nodo] = letra_para_el_hijo_de(&costos[nodo], letras[padre]);
                if detalle && bajados < 4 && por_nodo[nodo].is_empty() && narracion.hay_lugar() {
                    bajados += 1;
                    let nombre = abreviar(&etiquetas[nodo], 12);
                    let nombre_padre = abreviar(&etiquetas[padre], 12);
                    let elegidas: Vec<String> = orden
                        .iter()
                        .take_while(|&&otro| otro != nodo)
                        .chain(std::iter::once(&nodo))
                        .take(LIMITE_DE_LISTA)
                        .map(|&otro| {
                            format!(
                                "{}: {}",
                                abreviar(&etiquetas[otro], 10),
                                letra_de(letras[otro])
                            )
                        })
                        .collect();
                    let ultimo = elegidas.len() - 1;
                    narracion.agregar(
                        Pizarra::nueva()
                            .variable("posicion", posicion)
                            .variable("nodo", &nombre)
                            .variable("padre", &nombre_padre)
                            .lista("Letras elegidas hasta ahora", elegidas, &[ultimo])
                            .paso(
                                "letras[nodo] = letra_para_el_hijo(&costos[nodo], letras[padre]);",
                                format!(
                                    "Ahora bajamos hacia las hojas. «{nombre}» hereda la letra de su padre «{nombre_padre}» ({}) si eso le cuesta a lo sumo 1 más que su mejor opción; si no, toma su letra más barata. Se queda con {}.",
                                    letra_de(letras[padre]),
                                    letra_de(letras[nodo])
                                ),
                            ),
                    );
                }
            }
        }
        for (secuencia, &letra) in inferidas.iter_mut().zip(&letras) {
            secuencia[posicion] = b"ACGT"[letra];
        }
        if detalle && narracion.hay_lugar() {
            let todas: Vec<String> = orden
                .iter()
                .take(LIMITE_DE_LISTA)
                .map(|&otro| {
                    format!(
                        "{}: {}",
                        abreviar(&etiquetas[otro], 10),
                        letra_de(letras[otro])
                    )
                })
                .collect();
            narracion.agregar(
                Pizarra::nueva()
                    .variable("posicion", posicion)
                    .variable("cambios_por_posicion", conteos(&cambios_por_posicion))
                    .lista("Letras de todos los nodos en esta posición", todas, &[0])
                    .paso(
                        "anotar_letras(&mut inferidas, &letras, posicion);",
                        format!(
                            "La posición 0 está lista: copiamos la letra de cada nodo, antepasados incluidos, a su secuencia reconstruida. Esta posición necesitó {}. Las demás posiciones se resuelven igual.",
                            mutaciones(cambios)
                        ),
                    ),
            );
        } else if !detalle && (cambios > 0 || posicion == 1) && compactas < 24 {
            compactas += 1;
            let hijos = hijos_de(arbol, &padres, raiz);
            let mut filas_de_nodos: Vec<usize> = hijos.iter().copied().take(3).collect();
            filas_de_nodos.push(raiz);
            let propio = filas_de_nodos.len() - 1;
            let aviso = if cambios == 0 {
                "Todas las hojas coinciden o se pueden explicar sin mutar: 0 mutaciones."
                    .to_string()
            } else {
                if cambios == 1 {
                    "Hace falta 1 mutación.".to_string()
                } else {
                    format!("Hacen falta {cambios} mutaciones.")
                }
            };
            narracion.agregar(
                tabla_de_costos(
                    Pizarra::nueva()
                        .variable("posicion", posicion)
                        .variable("letras[raiz]", letra_de(letras[raiz]))
                        .variable("cambios_por_posicion", conteos(&cambios_por_posicion)),
                    etiquetas,
                    &costos,
                    &filas_de_nodos,
                    propio,
                )
                .secuencia(&por_nodo[primera], &marcas)
                .lista("Letras de las hojas en esta posición", letras_de_hojas, &[])
                .paso(
                    "cambios_por_posicion.push(costos[raiz].por_letra[letras[raiz]] as usize);",
                    format!(
                        "Posición {posicion}: mismo cálculo, de las hojas a la raíz y de vuelta. La letra más barata de la raíz es {}. {aviso}",
                        letra_de(letras[raiz])
                    ),
                ),
            );
        }
    }
    let puntaje: usize = cambios_por_posicion.iter().sum();
    let con_cambios = cambios_por_posicion
        .iter()
        .filter(|&&cambio| cambio > 0)
        .count();
    let por_posicion: Vec<String> = cambios_por_posicion
        .iter()
        .enumerate()
        .take(LIMITE_DE_LISTA)
        .map(|(posicion, cambio)| format!("posición {posicion}: {cambio}"))
        .collect();
    let resaltadas: Vec<usize> = (0..por_posicion.len())
        .filter(|&posicion| cambios_por_posicion[posicion] > 0)
        .collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("cambios_por_posicion", conteos(&cambios_por_posicion))
            .variable("puntaje", puntaje)
            .lista("Mutaciones en cada posición (las primeras)", por_posicion, &resaltadas)
            .paso(
                "let puntaje = cambios_por_posicion.iter().sum();",
                format!(
                    "El puntaje de parsimonia es la suma de las mutaciones de todas las posiciones: {} mutaciones. {} de las {} posiciones necesitaron al menos una; las otras no cambian nunca y no las mostramos.",
                    miles(puntaje),
                    miles(con_cambios),
                    miles(longitud)
                ),
            ),
    );
    (puntaje, inferidas, cambios_por_posicion)
}

fn ultimo_paso_de_parsimonia(
    arbol: &Arbol,
    etiquetas: &[String],
    resultado: &ResultadoDeParsimonia,
    aclaracion: &str,
) -> Paso {
    let (puntaje, inferidas, _) = resultado;
    let mut ramas = Vec::new();
    for (a, b, _) in arbol.aristas() {
        let mutaciones_en_la_rama = filogenia::distancia_de_hamming(&inferidas[a], &inferidas[b]);
        if mutaciones_en_la_rama > 0 {
            ramas.push(format!(
                "{}–{}: {}",
                abreviar(&etiquetas[a], 10),
                abreviar(&etiquetas[b], 10),
                mutaciones(mutaciones_en_la_rama)
            ));
        }
    }
    let cantidad = ramas.len();
    let visibles: Vec<String> = ramas.iter().take(LIMITE_DE_LISTA).cloned().collect();
    let todas: Vec<usize> = (0..visibles.len()).collect();
    Pizarra::nueva()
        .variable("puntaje", puntaje)
        .variable("con_cambios", format!("{cantidad} ramas con mutaciones"))
        .lista("Ramas con mutaciones", visibles, &todas)
        .paso(
            "let con_cambios = arbol_con_cambios(arbol, &inferidas);",
            format!(
                "Cada rama recibe el número de posiciones donde sus dos extremos tienen letras distintas: así se ve dónde ocurrieron las {} mutaciones. {aclaracion} Mira «Qué significa» para el índice de consistencia y la posición que más cambia.",
                miles(*puntaje)
            ),
        )
}

fn narrar_parsimonia_pequena(
    entrada: &Entrada,
) -> Result<(Vec<Paso>, ResultadoDeParsimonia), String> {
    let leido = filogenia::leer_newick(&entrada.texto("arbol")?)?;
    let (nombres, secuencias) = leer_secuencias(entrada, "secuencias", false)?;
    let por_nodo = secuencias_en_el_arbol(&leido.arbol, &nombres, secuencias)?;
    let resultado = filogenia::parsimonia_pequena(&leido.arbol, leido.raiz, &por_nodo)?;
    let mut con_nombres = resultado.arbol.clone();
    nombres_unicos(&mut con_nombres, Some(leido.raiz));
    let etiquetas = con_nombres.nombres;
    let mut narracion = Narracion::nueva(MAXIMO_DE_PASOS);
    let calculado = narrar_parsimonia(
        &mut narracion,
        &leido.arbol,
        leido.raiz,
        &por_nodo,
        &etiquetas,
    );
    let ultimo = ultimo_paso_de_parsimonia(
        &leido.arbol,
        &etiquetas,
        &calculado,
        "La secuencia de la raíz es una reconstrucción del antepasado común de todas las hojas.",
    );
    Ok((narracion.terminar(ultimo), calculado))
}

fn pasos_parsimonia_pequena(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    narrar_parsimonia_pequena(entrada).map(|(pasos, _)| pasos)
}

fn narrar_parsimonia_sin_raiz(
    entrada: &Entrada,
) -> Result<(Vec<Paso>, ResultadoDeParsimonia), String> {
    let leido = filogenia::leer_newick(&entrada.texto("arbol")?)?;
    let tenia_raiz = leido.arbol.grado(leido.raiz) == 2;
    let arbol = filogenia::quitar_raiz(&leido);
    let (nombres, secuencias) = leer_secuencias(entrada, "secuencias", false)?;
    let por_nodo = secuencias_en_el_arbol(&arbol, &nombres, secuencias)?;
    let resultado = filogenia::parsimonia_sin_raiz(&arbol, &por_nodo)?;
    let mut con_nombres = resultado.arbol.clone();
    nombres_unicos(&mut con_nombres, None);
    let etiquetas = con_nombres.nombres;
    let raiz = (0..arbol.nombres.len())
        .find(|&nodo| arbol.grado(nodo) >= 2)
        .unwrap_or(0);
    let mut narracion = Narracion::nueva(MAXIMO_DE_PASOS);
    let aviso = if tenia_raiz {
        " Tu Newick traía raíz (dos ramas arriba): se quitó, uniendo esas dos ramas en una sola."
    } else {
        ""
    };
    narracion.agregar(
        Pizarra::nueva()
            .variable("hojas", arbol.hojas().len())
            .paso(
                "for hoja in arbol.hojas() {",
                format!(
                    "Tu árbol no tiene raíz: no se sabe dónde está el antepasado más antiguo.{aviso} Primero se comprueba que cada una de sus {} hojas traiga secuencia.",
                    arbol.hojas().len()
                ),
            ),
    );
    let nombre_raiz = abreviar(&etiquetas[raiz], 12);
    narracion.agregar(
        Pizarra::nueva()
            .variable("raiz", format!("{raiz} («{nombre_raiz}»)"))
            .paso(
                ".find(|&nodo| arbol.grado(nodo) >= 2)",
                format!(
                    "Para poder calcular se elige como raíz el primer nodo interno, «{nombre_raiz}». Es solo un punto de partida: el puntaje sale igual poniendo la raíz en cualquier otro lugar."
                ),
            ),
    );
    let calculado = narrar_parsimonia(&mut narracion, &arbol, raiz, &por_nodo, &etiquetas);
    let ultimo = ultimo_paso_de_parsimonia(
        &arbol,
        &etiquetas,
        &calculado,
        "Como no hay raíz, las letras de los antepasados dependen de la raíz elegida, pero el puntaje no.",
    );
    Ok((narracion.terminar(ultimo), calculado))
}

fn pasos_parsimonia_sin_raiz(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    narrar_parsimonia_sin_raiz(entrada).map(|(pasos, _)| pasos)
}

fn narrar_vecinos_cercanos(entrada: &Entrada) -> Result<(Vec<Paso>, (Arbol, Arbol)), String> {
    let arbol = filogenia::leer_aristas(&entrada.texto("arbol")?)?;
    if arbol.hojas().len() > MAXIMO_DE_ESPECIES {
        return Err(format!(
            "El árbol tiene más de {MAXIMO_DE_ESPECIES} hojas; ese es el máximo aquí."
        ));
    }
    let (a, b) = leer_arista(&entrada.texto("arista")?, &arbol)?;
    filogenia::vecinos_cercanos(&arbol, a, b)?;
    let lado_a = lado_de(&arbol, a, b);
    let lado_b = lado_de(&arbol, b, a);
    let (w, x) = (lado_a[0], lado_a[1]);
    let (y, z) = (lado_b[0], lado_b[1]);
    let etiquetas = &arbol.nombres;
    let nombre = |nodo: usize| abreviar(&etiquetas[nodo], 12);
    let (na, nb, nw, nx, ny, nz) = (
        nombre(a),
        nombre(b),
        nombre(w),
        nombre(x),
        nombre(y),
        nombre(z),
    );
    let alrededor = |a_: usize, b_: usize| [a, b].contains(&a_) || [a, b].contains(&b_);
    let mut narracion = Narracion::nueva(MAXIMO_DE_PASOS);
    let (aristas, marcadas) = lista_de_aristas(
        &arbol,
        etiquetas,
        &alrededor,
        &|u, v| (u == a && v == b) || (u == b && v == a),
        LIMITE_DE_LISTA,
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("a", &na)
            .variable("b", &nb)
            .lista("Aristas alrededor de «a» y «b»", aristas, &marcadas)
            .paso(
                "revisar_arista_interna(arbol, a, b)?;",
                format!(
                    "Elegiste la arista «{na}»–«{nb}». Para mover subárboles tiene que ser interna: cada extremo debe tener tres vecinos. «{na}» tiene a «{nw}», «{nx}» y «{nb}»; «{nb}» tiene a «{ny}», «{nz}» y «{na}». Sí lo es."
                ),
            ),
    );
    let (aristas, marcadas) = lista_de_aristas(
        &arbol,
        etiquetas,
        &alrededor,
        &|u, v| (u == a && (v == w || v == x)) || (v == a && (u == w || u == x)),
        LIMITE_DE_LISTA,
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("a", &na)
            .variable("b", &nb)
            .variable("x", &nx)
            .lista("Aristas alrededor de «a» y «b»", aristas, &marcadas)
            .paso(
                "let (_, x) = otros_dos_vecinos(arbol, a, b);",
                format!(
                    "Los otros dos vecinos de «{na}» (sin contar a «{nb}») son «{nw}» y «{nx}». Cada uno es la raíz de un subárbol. Nos quedamos con «{nx}», que será el que se mueva."
                ),
            ),
    );
    let (aristas, marcadas) = lista_de_aristas(
        &arbol,
        etiquetas,
        &alrededor,
        &|u, v| (u == b && (v == y || v == z)) || (v == b && (u == y || u == z)),
        LIMITE_DE_LISTA,
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("a", &na)
            .variable("b", &nb)
            .variable("x", &nx)
            .variable("y", &ny)
            .variable("z", &nz)
            .lista("Aristas alrededor de «a» y «b»", aristas, &marcadas)
            .paso(
                "let (y, z) = otros_dos_vecinos(arbol, b, a);",
                format!(
                    "Del otro lado, los vecinos de «{nb}» (sin contar a «{na}») son «{ny}» y «{nz}». Habrá dos árboles vecinos: en uno «{nx}» cambia de lugar con «{ny}», y en el otro con «{nz}»."
                ),
            ),
    );
    let mut primero = arbol.clone();
    let (aristas, _) = lista_de_aristas(
        &primero,
        etiquetas,
        &alrededor,
        &|_, _| false,
        LIMITE_DE_LISTA,
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("primero", "copia del árbol")
            .lista("Aristas alrededor de «a» y «b»", aristas, &[])
            .paso(
                "let mut primero = arbol.clone();",
                "Hacemos una copia del árbol para el primer vecino. El original no se toca.",
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("a", &na)
            .variable("x", &nx)
            .variable("b", &nb)
            .variable("y", &ny)
            .paso(
                "intercambiar(&mut primero, a, x, b, y);",
                format!(
                    "Intercambiar significa: cortar la rama de «{na}» a «{nx}», cortar la de «{nb}» a «{ny}», y volver a pegar cada subárbol en el lado contrario. Lo hacemos en cuatro pasos."
                ),
            ),
    );
    let largo_x = primero.desconectar(a, x).unwrap_or(1.0);
    let (aristas, _) = lista_de_aristas(
        &primero,
        etiquetas,
        &alrededor,
        &|_, _| false,
        LIMITE_DE_LISTA,
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("largo_x", numero(largo_x))
            .lista("Aristas alrededor de «a» y «b»", aristas, &[])
            .paso(
                "let largo_x = arbol.desconectar(a, x).unwrap_or(1.0);",
                format!(
                    "Cortamos la rama «{na}»–«{nx}» y guardamos su longitud, {}: el subárbol de «{nx}» se lleva su rama consigo.",
                    numero(largo_x)
                ),
            ),
    );
    let largo_y = primero.desconectar(b, y).unwrap_or(1.0);
    let (aristas, _) = lista_de_aristas(
        &primero,
        etiquetas,
        &alrededor,
        &|_, _| false,
        LIMITE_DE_LISTA,
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("largo_x", numero(largo_x))
            .variable("largo_y", numero(largo_y))
            .lista("Aristas alrededor de «a» y «b»", aristas, &[])
            .paso(
                "let largo_y = arbol.desconectar(b, y).unwrap_or(1.0);",
                format!(
                    "Cortamos también «{nb}»–«{ny}», de longitud {}. Ahora los dos subárboles están sueltos.",
                    numero(largo_y)
                ),
            ),
    );
    primero.conectar(b, x, largo_x);
    let (aristas, marcadas) = lista_de_aristas(
        &primero,
        etiquetas,
        &alrededor,
        &|u, v| (u == b && v == x) || (u == x && v == b),
        LIMITE_DE_LISTA,
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("largo_x", numero(largo_x))
            .lista("Aristas alrededor de «a» y «b»", aristas, &marcadas)
            .paso(
                "arbol.conectar(b, x, largo_x);",
                format!(
                    "Pegamos «{nx}» a «{nb}», con la misma longitud que tenía ({}).",
                    numero(largo_x)
                ),
            ),
    );
    primero.conectar(a, y, largo_y);
    let (aristas, marcadas) = lista_de_aristas(
        &primero,
        etiquetas,
        &alrededor,
        &|u, v| {
            (u == b && v == x) || (u == x && v == b) || (u == a && v == y) || (u == y && v == a)
        },
        LIMITE_DE_LISTA,
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("largo_y", numero(largo_y))
            .lista("Aristas alrededor de «a» y «b»", aristas, &marcadas)
            .paso(
                "arbol.conectar(a, y, largo_y);",
                format!(
                    "Y pegamos «{ny}» a «{na}». Primer vecino listo: «{nx}» y «{ny}» cambiaron de lado."
                ),
            ),
    );
    let mut segundo = arbol.clone();
    filogenia::intercambiar(&mut segundo, a, x, b, z);
    let (aristas, marcadas) = lista_de_aristas(
        &segundo,
        etiquetas,
        &alrededor,
        &|u, v| {
            (u == b && v == x) || (u == x && v == b) || (u == a && v == z) || (u == z && v == a)
        },
        LIMITE_DE_LISTA,
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("segundo", "otra copia del árbol")
            .lista("Aristas alrededor de «a» y «b»", aristas, &marcadas)
            .paso(
                "intercambiar(&mut segundo, a, x, b, z);",
                format!(
                    "Con otra copia hacemos lo mismo, pero ahora «{nx}» cambia de lugar con «{nz}»: se cortan «{na}»–«{nx}» y «{nb}»–«{nz}» y se pegan «{nb}»–«{nx}» y «{na}»–«{nz}»."
                ),
            ),
    );
    let newick = |arbol_dibujado: &Arbol| abreviar(&newick_para_dibujar(arbol_dibujado, a), 120);
    let ultimo = Pizarra::nueva()
        .variable("primero", abreviar(&newick(&primero), 40))
        .variable("segundo", abreviar(&newick(&segundo), 40))
        .lista(
            "Árbol original y sus dos vecinos (Newick)",
            vec![newick(&arbol), newick(&primero), newick(&segundo)],
            &[1, 2],
        )
        .paso(
            "Ok((primero, segundo))",
            "Esos son los dos vecinos cercanos de la arista: mismo conjunto de hojas, con un subárbol de cada lado cambiado de lugar. Mira «Qué significa» para ver cuántos vecinos tiene el árbol en total.",
        );
    Ok((narracion.terminar(ultimo), (primero, segundo)))
}

fn pasos_vecinos_cercanos(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    narrar_vecinos_cercanos(entrada).map(|(pasos, _)| pasos)
}

type FilaDeFitch = (usize, usize, usize, u8, usize, bool);

fn fitch_de_posicion(
    arbol: &Arbol,
    conjuntos: &[Vec<u8>],
    posicion: usize,
) -> Option<(Vec<FilaDeFitch>, usize)> {
    let cantidad = arbol.nombres.len();
    let raiz = (0..cantidad).find(|&nodo| arbol.grado(nodo) == 1)?;
    let (orden, padres) = filogenia::recorrer_desde(arbol, raiz);
    let mut letras = vec![0u8; cantidad];
    let mut cambios = vec![0usize; cantidad];
    let mut filas = Vec::new();
    for &nodo in orden.iter().rev() {
        match padres[nodo] {
            Some(padre) if arbol.grado(nodo) == 3 => {
                let hijos = lado_de(arbol, nodo, padre);
                let (primero, segundo) = (hijos[0], hijos[1]);
                let comunes = letras[primero] & letras[segundo];
                let sin_comunes = comunes == 0;
                letras[nodo] = if sin_comunes {
                    letras[primero] | letras[segundo]
                } else {
                    comunes
                };
                cambios[nodo] = cambios[primero] + cambios[segundo] + usize::from(sin_comunes);
                filas.push((
                    nodo,
                    primero,
                    segundo,
                    letras[nodo],
                    cambios[nodo],
                    sin_comunes,
                ));
            }
            _ => {
                letras[nodo] = conjuntos[nodo][posicion];
                cambios[nodo] = 0;
            }
        }
    }
    let hijo = arbol.vecinos[raiz][0].0;
    let total = cambios[hijo] + usize::from(conjuntos[raiz][posicion] & letras[hijo] == 0);
    Some((filas, total))
}

type CandidatoDeVecino = (usize, usize, String, String, usize);

fn vecinos_con_puntaje(
    arbol: &Arbol,
    conjuntos: &[Vec<u8>],
    etiquetas: &[String],
) -> Vec<CandidatoDeVecino> {
    let mut filas = Vec::new();
    for (u, v, _) in arbol.aristas() {
        if arbol.grado(u) != 3 || arbol.grado(v) != 3 {
            continue;
        }
        let Ok((primero, segundo)) = filogenia::vecinos_cercanos(arbol, u, v) else {
            continue;
        };
        let x = lado_de(arbol, u, v)[1];
        let lado = lado_de(arbol, v, u);
        let arista = format!(
            "{}–{}",
            abreviar(&etiquetas[u], 10),
            abreviar(&etiquetas[v], 10)
        );
        for (otro, vecino) in [(lado[0], primero), (lado[1], segundo)] {
            filas.push((
                u,
                v,
                arista.clone(),
                format!(
                    "{} ↔ {}",
                    abreviar(&etiquetas[x], 10),
                    abreviar(&etiquetas[otro], 10)
                ),
                filogenia::puntaje_de_fitch(&vecino, conjuntos),
            ));
        }
    }
    filas
}

type ResultadoDeBusqueda = (Arbol, Vec<filogenia::Movimiento>, usize);

fn narrar_parsimonia_grande(entrada: &Entrada) -> Result<(Vec<Paso>, ResultadoDeBusqueda), String> {
    let (nombres, secuencias) = leer_secuencias(entrada, "secuencias", false)?;
    if nombres.len() < 3 {
        return Err("Hacen falta al menos tres secuencias.".to_string());
    }
    revisar_mismo_largo(&nombres, &secuencias)?;
    let inicio = entrada.opcion("inicio", INICIOS)?;
    let semilla = entrada.natural("semilla")? as u64;
    let maximo = entrada.natural("movimientos")?;
    if maximo > MAXIMO_DE_MOVIMIENTOS {
        return Err(format!(
            "El máximo de movimientos debe estar entre 0 y {MAXIMO_DE_MOVIMIENTOS}."
        ));
    }
    let al_azar = inicio == "al azar";
    let inicial = if al_azar {
        filogenia::arbol_al_azar(&nombres, semilla)?
    } else {
        let distancias_p = filogenia::matriz_de_distancias_p(&secuencias)?;
        filogenia::union_de_vecinos(&distancias_p, &nombres)?
    };
    let mut por_nodo = secuencias;
    por_nodo.resize(inicial.nombres.len(), Vec::new());
    let longitud = filogenia::revisar_arbol_binario(&inicial, &por_nodo)?;
    let cantidad = nombres.len();
    let mut con_nombres = inicial.clone();
    nombres_unicos(&mut con_nombres, None);
    let etiquetas = con_nombres.nombres.clone();
    let origen = if al_azar {
        format!("un árbol al azar (semilla {semilla})")
    } else {
        "el que arma unión de vecinos con las distancias p".to_string()
    };
    let mut narracion = Narracion::nueva(MAXIMO_DE_PASOS);
    let elementos: Vec<String> = (0..cantidad)
        .take(LIMITE_DE_LISTA)
        .map(|nodo| {
            format!(
                "{}: {}",
                abreviar(&etiquetas[nodo], 10),
                abreviar(&String::from_utf8_lossy(&por_nodo[nodo]), 24)
            )
        })
        .collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("longitud", miles(longitud))
            .variable(
                "inicial",
                abreviar(&newick_para_dibujar(&con_nombres, raiz_para_dibujar(&con_nombres)), 40),
            )
            .variable("maximo_de_movimientos", maximo)
            .lista("Secuencias", elementos, &[])
            .paso(
                "let longitud = revisar_arbol_binario(inicial, secuencias)?;",
                format!(
                    "Tienes {cantidad} especies con secuencias de {} letras. Con tantas especies hay {} árboles sin raíz posibles: demasiados para probarlos todos. Así que partimos de {origen} y buscamos mejoras pequeñas, hasta {maximo} movimientos.",
                    miles(longitud),
                    arboles_posibles(cantidad)
                ),
            ),
    );
    let conjuntos = filogenia::conjuntos_de_letras(&por_nodo, longitud);
    let conjuntos_vistos: Vec<String> = (0..cantidad)
        .take(4)
        .map(|nodo| {
            format!(
                "{}: {}",
                abreviar(&etiquetas[nodo], 10),
                (0..longitud.min(6))
                    .map(|posicion| letras_en_texto(conjuntos[nodo][posicion]))
                    .collect::<Vec<_>>()
                    .join(" ")
            )
        })
        .collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("longitud", miles(longitud))
            .variable("conjuntos", format!("{} nodos × {} posiciones", por_nodo.len(), miles(longitud)))
            .lista("Conjuntos de letras (primeras posiciones)", conjuntos_vistos, &[])
            .paso(
                "let conjuntos = conjuntos_de_letras(secuencias, longitud);",
                "Cada letra se guarda como un conjunto de letras posibles: {A}, {C}, {G} o {T}. En las hojas el conjunto tiene una sola letra. En los nodos internos, que no tienen secuencia, empiezan vacíos y se irán calculando.",
            ),
    );
    let mut arbol = inicial.clone();
    let puntaje_inicial = filogenia::puntaje_de_fitch(&arbol, &conjuntos);
    let mut puntaje = puntaje_inicial;
    if let Some((filas, total)) = fitch_de_posicion(&arbol, &conjuntos, 0) {
        let tabla: Vec<Vec<String>> = filas
            .iter()
            .take(8)
            .map(|&(nodo, primero, segundo, letras, cambios, _)| {
                vec![
                    abreviar(&etiquetas[nodo], 10),
                    format!(
                        "{} + {}",
                        abreviar(&etiquetas[primero], 8),
                        abreviar(&etiquetas[segundo], 8)
                    ),
                    letras_en_texto(letras),
                    cambios.to_string(),
                ]
            })
            .collect();
        let resaltadas: Vec<(usize, usize)> = filas
            .iter()
            .take(8)
            .enumerate()
            .filter(|(_, fila)| fila.5)
            .map(|(indice, _)| (indice, 2))
            .collect();
        narracion.agregar(
            Pizarra::nueva()
                .variable("puntaje_inicial", puntaje_inicial)
                .variable("puntaje", puntaje)
                .tabla(
                    "Posición 0, de las hojas hacia arriba",
                    &["Nodo", "Junta a", "Letras", "Cambios hasta aquí"],
                    tabla,
                    &resaltadas,
                )
                .paso(
                    "let puntaje_inicial = puntaje_de_fitch(&arbol, &conjuntos);",
                    format!(
                        "Para medir un árbol se cuenta, posición por posición, cuántas mutaciones hacen falta. En cada nodo se juntan los conjuntos de sus dos hijos: si comparten alguna letra se queda con las compartidas, y si no, toma la unión y suma 1 mutación (las celdas marcadas). La posición 0 necesita {total}; sumando todas las posiciones, el árbol inicial vale {} mutaciones.",
                        miles(puntaje_inicial)
                    ),
                ),
        );
    }
    let mut movimientos: Vec<filogenia::Movimiento> = Vec::new();
    let mut termino_sin_mejora = false;
    while movimientos.len() < maximo {
        let ronda = movimientos.len();
        let detallada = ronda < 3;
        let siguiente = filogenia::mejor_intercambio(&arbol, &conjuntos, puntaje);
        let candidatos = if detallada && cantidad <= 24 {
            vecinos_con_puntaje(&arbol, &conjuntos, &etiquetas)
        } else {
            Vec::new()
        };
        let Some(movimiento) = siguiente else {
            let tabla: Vec<Vec<String>> = candidatos
                .iter()
                .take(LIMITE_DE_LISTA)
                .map(|(_, _, arista, intercambio, valor)| {
                    vec![arista.clone(), intercambio.clone(), valor.to_string()]
                })
                .collect();
            let pizarra = Pizarra::nueva()
                .variable("puntaje", puntaje)
                .variable("movimientos", movimientos.len());
            let pizarra = if tabla.is_empty() {
                pizarra
            } else {
                pizarra.tabla(
                    "Vecinos cercanos y su puntaje",
                    &["Arista", "Intercambio", "Puntaje"],
                    tabla,
                    &[],
                )
            };
            narracion.agregar(pizarra.paso(
                "break;",
                format!(
                    "Revisamos todos los vecinos cercanos del árbol actual y ninguno baja de {} mutaciones. Estamos en un óptimo local: ningún cambio pequeño mejora, así que la búsqueda se detiene.",
                    miles(puntaje)
                ),
            ));
            termino_sin_mejora = true;
            break;
        };
        let (a, b) = movimiento.arista;
        let (x, y) = movimiento.intercambio;
        let (na, nb, nx, ny) = (
            abreviar(&etiquetas[a], 12),
            abreviar(&etiquetas[b], 12),
            abreviar(&etiquetas[x], 12),
            abreviar(&etiquetas[y], 12),
        );
        if detallada && narracion.hay_lugar() {
            let mejor_fila = candidatos.iter().position(|(u, v, _, _, valor)| {
                ((*u == a && *v == b) || (*u == b && *v == a)) && *valor == movimiento.puntaje
            });
            let tabla: Vec<Vec<String>> = candidatos
                .iter()
                .take(LIMITE_DE_LISTA)
                .map(|(_, _, arista, intercambio, valor)| {
                    vec![arista.clone(), intercambio.clone(), valor.to_string()]
                })
                .collect();
            let resaltadas: Vec<(usize, usize)> = mejor_fila
                .filter(|&fila| fila < LIMITE_DE_LISTA)
                .map(|fila| (fila, 2))
                .into_iter()
                .collect();
            let pizarra = Pizarra::nueva()
                .variable("puntaje", puntaje)
                .variable("movimientos", movimientos.len());
            let pizarra = if tabla.is_empty() {
                pizarra
            } else {
                pizarra.tabla(
                    "Vecinos cercanos y su puntaje",
                    &["Arista", "Intercambio", "Puntaje"],
                    tabla,
                    &resaltadas,
                )
            };
            narracion.agregar(pizarra.paso(
                "let Some(movimiento) = mejor_intercambio(&arbol, &conjuntos, puntaje) else {",
                format!(
                    "Movimiento {}: probamos los dos vecinos cercanos de cada arista interna y medimos su puntaje (el árbol actual vale {}). El mejor está en la arista «{na}»–«{nb}», cambiando «{nx}» por «{ny}», y baja el puntaje a {}.",
                    ronda + 1,
                    miles(puntaje),
                    miles(movimiento.puntaje)
                ),
            ));
        }
        filogenia::intercambiar(&mut arbol, a, x, b, y);
        let viejo = puntaje;
        puntaje = movimiento.puntaje;
        movimientos.push(movimiento);
        let visible = detallada || (ronda + 1).is_multiple_of(10);
        if visible && narracion.hay_lugar() {
            let (aristas, marcadas) = lista_de_aristas(
                &arbol,
                &etiquetas,
                &|u, v| [a, b].contains(&u) || [a, b].contains(&v),
                &|u, v| {
                    (u == b && v == x)
                        || (u == x && v == b)
                        || (u == a && v == y)
                        || (u == y && v == a)
                },
                LIMITE_DE_LISTA,
            );
            narracion.agregar(
                Pizarra::nueva()
                    .variable("a", &na)
                    .variable("b", &nb)
                    .variable("x", &nx)
                    .variable("y", &ny)
                    .lista("Aristas alrededor de «a» y «b»", aristas, &marcadas)
                    .paso(
                        "intercambiar(&mut arbol, a, x, b, y);",
                        format!(
                            "Hacemos el intercambio: «{nx}» pasa del lado de «{na}» al de «{nb}», y «{ny}» al revés, cada uno con su subárbol y su rama."
                        ),
                    ),
            );
        }
        if detallada && narracion.hay_lugar() {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("puntaje", puntaje)
                    .variable("movimientos", movimientos.len())
                    .paso(
                        "puntaje = movimiento.puntaje;",
                        format!(
                            "El árbol nuevo necesita {} mutaciones, {} menos que antes. Lo adoptamos y buscamos otra mejora a partir de él.",
                            miles(puntaje),
                            miles(viejo - puntaje)
                        ),
                    ),
            );
        }
    }
    if !termino_sin_mejora {
        let nota = if maximo == 0 {
            "Pediste 0 movimientos, así que solo se midió el árbol inicial.".to_string()
        } else {
            format!(
                "Llegamos al límite de {maximo} movimientos y nos detenemos. Quizá quedaban mejoras: sube el límite para comprobarlo."
            )
        };
        narracion.agregar(
            Pizarra::nueva()
                .variable("puntaje", puntaje)
                .variable("movimientos", movimientos.len())
                .paso("while movimientos.len() < maximo_de_movimientos {", nota),
        );
    }
    let parsimonia = filogenia::parsimonia_sin_raiz(&arbol, &por_nodo)?;
    let hechos: Vec<String> = movimientos
        .iter()
        .enumerate()
        .take(LIMITE_DE_LISTA)
        .map(|(numero_de_movimiento, movimiento)| {
            format!(
                "{}: {} ↔ {} → {}",
                numero_de_movimiento + 1,
                abreviar(&etiquetas[movimiento.intercambio.0], 8),
                abreviar(&etiquetas[movimiento.intercambio.1], 8),
                movimiento.puntaje
            )
        })
        .collect();
    let ultimo = Pizarra::nueva()
        .variable("puntaje_inicial", puntaje_inicial)
        .variable("puntaje", puntaje)
        .variable("movimientos", movimientos.len())
        .lista("Movimientos hechos (intercambio → puntaje)", hechos, &[])
        .paso(
            "let parsimonia = parsimonia_sin_raiz(&arbol, secuencias)?;",
            format!(
                "Para terminar se vuelve a medir el árbol final con la parsimonia pequeña, que también da las letras de los antepasados: {} mutaciones, contra {} del árbol inicial, tras {} movimientos. Mira «Qué significa» para ver si puede haber un árbol mejor.",
                miles(parsimonia.puntaje),
                miles(puntaje_inicial),
                movimientos.len()
            ),
        );
    Ok((
        narracion.terminar(ultimo),
        (arbol, movimientos, parsimonia.puntaje),
    ))
}

fn pasos_parsimonia_grande(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    narrar_parsimonia_grande(entrada).map(|(pasos, _)| pasos)
}

#[cfg(test)]
mod pruebas {
    use serde_json::json;

    fn interpretacion(id: &str, entrada: serde_json::Value) -> String {
        crate::ejecutar(id, entrada)
            .unwrap()
            .interpretacion
            .join("\n")
    }

    #[test]
    fn la_rama_de_dos_hojas_no_se_llama_larga() {
        let texto = interpretacion(
            "filogenia.longitud_de_rama",
            json!({ "matriz": "0 4\n4 0", "hoja": "0" }),
        );
        assert!(texto.contains("solo dos hojas"));
        assert!(!texto.contains("rama larga"));
    }

    #[test]
    fn la_rama_usa_la_mejor_pareja_de_todas() {
        let texto = interpretacion(
            "filogenia.longitud_de_rama",
            json!({ "matriz": "0 5 5 5\n5 0 10 1\n5 10 0 1\n5 1 1 0", "hoja": "3" }),
        );
        assert!(texto.contains("mide -4"));
        assert!(texto.contains("i = «1» y k = «2»"));
    }

    #[test]
    fn dos_especies_aditivas_hablan_en_singular() {
        let texto = interpretacion(
            "filogenia.filogenia_aditiva",
            json!({ "matriz": "0 1\n1 0" }),
        );
        assert!(!texto.contains("las 1 distancias"));
    }

    #[test]
    fn secuencias_identicas_no_tienen_rama_mas_larga() {
        let texto = interpretacion(
            "filogenia.arbol_desde_secuencias",
            json!({ "secuencias": "a=ACGT\nb=ACGT\nc=ACGT" }),
        );
        assert!(texto.contains("idénticas"));
        assert!(!texto.contains("más larga"));
        assert!(!texto.contains("Las más parecidas"));
    }

    #[test]
    fn cuatro_especies_sin_movimientos_no_son_optimo_global() {
        let texto = interpretacion(
            "filogenia.parsimonia_grande",
            json!({
                "secuencias": "a=ACGT\nb=ACGA\nc=ACCA\nd=TCCA",
                "inicio": "al azar",
                "semilla": "1",
                "movimientos": "0"
            }),
        );
        assert!(!texto.contains("óptimo global"));
    }

    use super::*;

    fn entrada(valores: serde_json::Value) -> Entrada {
        Entrada::desde_json(valores).unwrap()
    }

    fn de_ejemplo(id: &str) -> serde_json::Value {
        crate::ejemplo(id).unwrap()
    }

    fn revisar_pasos(id: &str, valores: serde_json::Value) -> Vec<Paso> {
        let pasos = crate::pasos(id, valores).unwrap();
        let codigo: String = crate::codigo(id)
            .unwrap()
            .iter()
            .map(|fragmento| fragmento.fuente.as_str())
            .collect::<Vec<_>>()
            .join("\n");
        assert!(
            !pasos.is_empty() && pasos.len() <= 400,
            "{id}: {} pasos",
            pasos.len()
        );
        for paso in &pasos {
            assert!(
                codigo
                    .lines()
                    .any(|linea| linea.contains(paso.linea.as_str())),
                "{id}: la línea «{}» no está en el código",
                paso.linea
            );
            assert!(!paso.nota.is_empty(), "{id}: nota vacía");
        }
        assert!(serde_json::to_string(&pasos).unwrap().len() < 2_000_000);
        pasos
    }

    fn numeros_falsos(semilla: u64, cantidad: usize) -> Vec<f64> {
        let mut estado = semilla;
        (0..cantidad)
            .map(|_| {
                estado = estado
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                ((estado >> 33) % 1000) as f64 / 10.0
            })
            .collect()
    }

    fn matriz_en_linea(puntos: &[f64], ruido: &[f64]) -> String {
        let cantidad = puntos.len();
        let mut texto = String::new();
        for i in 0..cantidad {
            let fila: Vec<String> = (0..cantidad)
                .map(|j| {
                    if i == j {
                        "0".to_string()
                    } else {
                        let (a, b) = (i.min(j), i.max(j));
                        format!(
                            "{}",
                            (puntos[a] - puntos[b]).abs() + ruido[a * cantidad + b]
                        )
                    }
                })
                .collect();
            texto.push_str(&fila.join(" "));
            texto.push('\n');
        }
        texto
    }

    fn alineamiento_falso(cantidad: usize, largo: usize, semilla: u64) -> String {
        let letras = numeros_falsos(semilla, largo * (cantidad + 1));
        let base: Vec<u8> = (0..largo)
            .map(|posicion| b"ACGT"[(letras[posicion] as usize) % 4])
            .collect();
        let mut texto = String::new();
        for especie in 0..cantidad {
            let mut secuencia = base.clone();
            for posicion in 0..largo {
                let tirada = letras[(especie + 1) * largo + posicion] as usize;
                if tirada.is_multiple_of(11) {
                    secuencia[posicion] = b"ACGT"[tirada % 4];
                }
            }
            texto.push_str(&format!(
                "e{especie}={}\n",
                String::from_utf8_lossy(&secuencia)
            ));
        }
        texto
    }

    #[test]
    fn los_ejemplos_narran_con_lineas_del_codigo() {
        for id in [
            "filogenia.distancias_entre_hojas",
            "filogenia.longitud_de_rama",
            "filogenia.filogenia_aditiva",
            "filogenia.upgma",
            "filogenia.distancias_evolutivas",
            "filogenia.arbol_desde_secuencias",
            "filogenia.parsimonia_pequena",
            "filogenia.parsimonia_sin_raiz",
            "filogenia.vecinos_cercanos",
            "filogenia.parsimonia_grande",
        ] {
            let pasos = revisar_pasos(id, de_ejemplo(id));
            assert!(pasos.len() >= 5, "{id}: muy pocos pasos");
        }
    }

    #[test]
    fn otros_casos_narran_con_lineas_del_codigo() {
        let puntos = numeros_falsos(3, 40);
        let ruido = numeros_falsos(4, 40 * 40);
        let casos = [
            (
                "filogenia.filogenia_aditiva",
                json!({ "matriz": MATRIZ_CON_RUIDO }),
            ),
            (
                "filogenia.filogenia_aditiva",
                json!({ "matriz": "0 3\n3 0" }),
            ),
            (
                "filogenia.filogenia_aditiva",
                json!({ "matriz": matriz_en_linea(&puntos, &vec![0.0; 1600]) }),
            ),
            (
                "filogenia.filogenia_aditiva",
                json!({ "matriz": "0 1 5\n1 0 1\n5 1 0" }),
            ),
            ("filogenia.upgma", json!({ "matriz": "0 3\n3 0" })),
            (
                "filogenia.upgma",
                json!({ "matriz": matriz_en_linea(&puntos, &ruido) }),
            ),
            (
                "filogenia.union_de_vecinos",
                json!({ "matriz": "0 3\n3 0" }),
            ),
            (
                "filogenia.longitud_de_rama",
                json!({ "matriz": "0 4\n4 0", "hoja": "0" }),
            ),
            (
                "filogenia.longitud_de_rama",
                json!({ "matriz": matriz_en_linea(&puntos, &ruido), "hoja": "17" }),
            ),
            (
                "filogenia.distancias_evolutivas",
                json!({ "secuencias": "a=AC-GT\nb=ACGGT\nc=TCAGA" }),
            ),
            (
                "filogenia.distancias_evolutivas",
                json!({ "secuencias": "a=ACGT\nb=ACGT" }),
            ),
            (
                "filogenia.distancias_evolutivas",
                json!({ "secuencias": "a=AAAA\nb=CCCC\nc=GGGG" }),
            ),
            (
                "filogenia.distancias_evolutivas",
                json!({ "secuencias": alineamiento_falso(30, 400, 5) }),
            ),
            (
                "filogenia.arbol_desde_secuencias",
                json!({ "secuencias": "a=AAAA\nb=CCCC\nc=GGGG" }),
            ),
            (
                "filogenia.arbol_desde_secuencias",
                json!({ "secuencias": "a=ACGT\nb=ACGT\nc=ACGT" }),
            ),
            (
                "filogenia.arbol_desde_secuencias",
                json!({ "secuencias": "a=ACGTACGT\nb=ACGTACGA" }),
            ),
            (
                "filogenia.arbol_desde_secuencias",
                json!({ "secuencias": alineamiento_falso(40, 300, 6) }),
            ),
            (
                "filogenia.parsimonia_pequena",
                json!({ "arbol": "(a,b,c,d);", "secuencias": "a=AC\nb=AG\nc=TC\nd=TG" }),
            ),
            (
                "filogenia.parsimonia_pequena",
                json!({ "arbol": "(a,b);", "secuencias": "a=AC\nb=AG" }),
            ),
            (
                "filogenia.parsimonia_sin_raiz",
                json!({ "arbol": ARBOL_CON_RAIZ, "secuencias": SECUENCIAS_DE_HOJAS }),
            ),
            (
                "filogenia.parsimonia_sin_raiz",
                json!({ "arbol": "(a,b,c);", "secuencias": "a=AC\nb=AG\nc=TC" }),
            ),
            (
                "filogenia.vecinos_cercanos",
                json!({ "arbol": ARBOL_BINARIO, "arista": "n1 n2" }),
            ),
            (
                "filogenia.vecinos_cercanos",
                json!({ "arbol": ARBOL_BINARIO, "arista": "n3 n4" }),
            ),
            (
                "filogenia.parsimonia_grande",
                json!({ "secuencias": OCHO_ESPECIES, "inicio": "unión de vecinos", "semilla": "1", "movimientos": "100" }),
            ),
            (
                "filogenia.parsimonia_grande",
                json!({ "secuencias": OCHO_ESPECIES, "inicio": "al azar", "semilla": "3", "movimientos": "0" }),
            ),
            (
                "filogenia.parsimonia_grande",
                json!({ "secuencias": "a=ACGT\nb=ACGA\nc=ACCA", "inicio": "al azar", "semilla": "1", "movimientos": "5" }),
            ),
            (
                "filogenia.parsimonia_grande",
                json!({ "secuencias": "a=ACGT\nb=ACGA\nc=ACCA\nd=TCCA", "inicio": "al azar", "semilla": "1", "movimientos": "5" }),
            ),
            (
                "filogenia.parsimonia_grande",
                json!({ "secuencias": alineamiento_falso(40, 200, 8), "inicio": "al azar", "semilla": "2", "movimientos": "200" }),
            ),
        ];
        for (id, valores) in casos {
            revisar_pasos(id, valores);
        }
        let arbol_grande: String = (0..60)
            .map(|hoja| format!("h{hoja}->r{}:{}\n", hoja / 2, 1 + hoja % 3))
            .chain((0..30).map(|nodo| format!("r{nodo}->s{}:2\n", nodo / 2)))
            .chain((0..15).map(|nodo| format!("s{nodo}->t{}:1\n", nodo / 2)))
            .chain((0..8).map(|nodo| format!("t{nodo}->z:1\n")))
            .collect();
        revisar_pasos(
            "filogenia.distancias_entre_hojas",
            json!({ "arbol": arbol_grande }),
        );
        let secuencias_largas = alineamiento_falso(12, 1500, 9);
        let mut arbol_de_secuencias = String::from("(");
        for especie in 0..12 {
            if especie > 0 {
                arbol_de_secuencias.push(',');
            }
            arbol_de_secuencias.push_str(&format!("e{especie}"));
        }
        arbol_de_secuencias.push_str(");");
        revisar_pasos(
            "filogenia.parsimonia_pequena",
            json!({ "arbol": arbol_de_secuencias, "secuencias": secuencias_largas }),
        );
    }

    #[test]
    fn la_narracion_de_distancias_entre_hojas_da_la_misma_matriz() {
        for texto in [ARBOL_CON_PESOS, "a->b:1\nb->c:2"] {
            let valores = entrada(json!({ "arbol": texto }));
            let (_, matriz) = narrar_distancias_entre_hojas(&valores).unwrap();
            let arbol = filogenia::leer_aristas(texto).unwrap();
            let hojas = arbol.hojas();
            assert_eq!(matriz, filogenia::distancias_entre_hojas(&arbol, &hojas));
            for hoja in hojas {
                assert_eq!(
                    recorrido_desde(&arbol, hoja).0,
                    filogenia::distancias_desde(&arbol, hoja)
                );
            }
        }
    }

    #[test]
    fn la_narracion_de_la_rama_da_la_misma_rama() {
        let puntos = numeros_falsos(11, 12);
        let ruido = numeros_falsos(12, 144);
        for (matriz, hoja) in [
            (MATRIZ_ADITIVA.to_string(), "Cardenix".to_string()),
            (MATRIZ_CON_RUIDO.to_string(), "2".to_string()),
            ("0 4\n4 0".to_string(), "1".to_string()),
            (matriz_en_linea(&puntos, &ruido), "7".to_string()),
        ] {
            let valores = entrada(json!({ "matriz": matriz, "hoja": hoja }));
            let (_, mejor) = narrar_longitud_de_rama(&valores).unwrap();
            let (nombres, distancias) = leer_matriz(&valores, "matriz").unwrap();
            let indice = buscar_especie(&hoja, &nombres).unwrap();
            assert_eq!(
                mejor,
                filogenia::longitud_de_rama_con_par(&distancias, indice).unwrap()
            );
        }
    }

    #[test]
    fn la_narracion_aditiva_da_el_mismo_arbol() {
        let puntos = numeros_falsos(13, 15);
        for matriz in [
            MATRIZ_ADITIVA.to_string(),
            MATRIZ_CON_RUIDO.to_string(),
            "0 3\n3 0".to_string(),
            matriz_en_linea(&puntos, &vec![0.0; 225]),
            "0 1 5\n1 0 1\n5 1 0".to_string(),
        ] {
            let valores = entrada(json!({ "matriz": matriz }));
            let (_, (arbol, reproduce)) = narrar_filogenia_aditiva(&valores).unwrap();
            let (nombres, distancias) = leer_matriz(&valores, "matriz").unwrap();
            match filogenia::filogenia_aditiva(&distancias, &nombres) {
                Ok(esperado) => {
                    assert!(reproduce);
                    assert_eq!(arbol, esperado);
                }
                Err(_) => assert!(!reproduce),
            }
        }
    }

    #[test]
    fn la_narracion_de_upgma_da_el_mismo_arbol() {
        let puntos = numeros_falsos(14, 20);
        let ruido = numeros_falsos(15, 400);
        for matriz in [
            MATRIZ_CON_RUIDO.to_string(),
            MATRIZ_ADITIVA.to_string(),
            "0 3\n3 0".to_string(),
            matriz_en_linea(&puntos, &ruido),
        ] {
            let valores = entrada(json!({ "matriz": matriz }));
            let (_, resultado) = narrar_upgma(&valores).unwrap();
            let (nombres, distancias) = leer_matriz(&valores, "matriz").unwrap();
            assert_eq!(resultado, filogenia::upgma(&distancias, &nombres).unwrap());
        }
    }

    #[test]
    fn la_narracion_de_secuencias_da_las_mismas_distancias_y_arbol() {
        let largo = alineamiento_falso(15, 300, 16);
        for texto in [
            SECUENCIAS_ALINEADAS.to_string(),
            "a=AC-GT\nb=ACGGT\nc=TCAGA".to_string(),
            "a=ACGT\nb=ACGT\nc=ACGT".to_string(),
            "a=AAAA\nb=CCCC\nc=GGGG".to_string(),
            largo,
        ] {
            let valores = entrada(json!({ "secuencias": texto }));
            let (_, jukes_cantor) = narrar_distancias_evolutivas(&valores).unwrap();
            let (nombres, secuencias) = leer_secuencias(&valores, "secuencias", true).unwrap();
            assert_eq!(
                jukes_cantor,
                filogenia::matriz_jukes_cantor(&secuencias).unwrap()
            );
            let (_, arbol) = narrar_arbol_desde_secuencias(&valores).unwrap();
            assert_eq!(
                arbol,
                filogenia::arbol_desde_secuencias(&secuencias, &nombres).ok()
            );
        }
    }

    #[test]
    fn la_narracion_de_parsimonia_da_el_mismo_puntaje_y_letras() {
        for (arbol, secuencias) in [
            (ARBOL_CON_RAIZ, SECUENCIAS_DE_HOJAS.to_string()),
            ("(a,b,c,d);", "a=AC\nb=AG\nc=TC\nd=TG".to_string()),
            ("(a,b);", "a=AC\nb=AG".to_string()),
            (
                "((a,b),(c,d));",
                "a=AAAA\nb=AAAT\nc=TTTT\nd=TTTA".to_string(),
            ),
        ] {
            let valores = entrada(json!({ "arbol": arbol, "secuencias": secuencias }));
            let (_, (puntaje, inferidas, cambios)) = narrar_parsimonia_pequena(&valores).unwrap();
            let leido = filogenia::leer_newick(arbol).unwrap();
            let (nombres, hojas) = leer_secuencias(&valores, "secuencias", false).unwrap();
            let por_nodo = secuencias_en_el_arbol(&leido.arbol, &nombres, hojas).unwrap();
            let esperado =
                filogenia::parsimonia_pequena(&leido.arbol, leido.raiz, &por_nodo).unwrap();
            assert_eq!(puntaje, esperado.puntaje);
            assert_eq!(inferidas, esperado.secuencias);
            assert_eq!(cambios, esperado.cambios_por_posicion);
        }
        for (arbol, secuencias) in [
            (ARBOL_SIN_RAIZ, SECUENCIAS_DE_HOJAS.to_string()),
            (ARBOL_CON_RAIZ, SECUENCIAS_DE_HOJAS.to_string()),
            ("(a,b,c);", "a=AC\nb=AG\nc=TC".to_string()),
        ] {
            let valores = entrada(json!({ "arbol": arbol, "secuencias": secuencias }));
            let (_, (puntaje, inferidas, cambios)) = narrar_parsimonia_sin_raiz(&valores).unwrap();
            let leido = filogenia::leer_newick(arbol).unwrap();
            let sin_raiz = filogenia::quitar_raiz(&leido);
            let (nombres, hojas) = leer_secuencias(&valores, "secuencias", false).unwrap();
            let por_nodo = secuencias_en_el_arbol(&sin_raiz, &nombres, hojas).unwrap();
            let esperado = filogenia::parsimonia_sin_raiz(&sin_raiz, &por_nodo).unwrap();
            assert_eq!(puntaje, esperado.puntaje);
            assert_eq!(inferidas, esperado.secuencias);
            assert_eq!(cambios, esperado.cambios_por_posicion);
        }
    }

    #[test]
    fn la_narracion_de_vecinos_da_los_mismos_arboles() {
        for arista in ["n2 n3", "n1 n2", "n3 n4"] {
            let valores = entrada(json!({ "arbol": ARBOL_BINARIO, "arista": arista }));
            let (_, (primero, segundo)) = narrar_vecinos_cercanos(&valores).unwrap();
            let arbol = filogenia::leer_aristas(ARBOL_BINARIO).unwrap();
            let (a, b) = leer_arista(arista, &arbol).unwrap();
            let (esperado_primero, esperado_segundo) =
                filogenia::vecinos_cercanos(&arbol, a, b).unwrap();
            assert_eq!(primero, esperado_primero);
            assert_eq!(segundo, esperado_segundo);
        }
    }

    #[test]
    fn la_narracion_de_la_busqueda_da_el_mismo_resultado() {
        let largo = alineamiento_falso(14, 120, 21);
        for (secuencias, inicio, semilla, movimientos) in [
            (OCHO_ESPECIES.to_string(), "al azar", "7", "100"),
            (OCHO_ESPECIES.to_string(), "unión de vecinos", "7", "100"),
            (OCHO_ESPECIES.to_string(), "al azar", "2", "0"),
            (OCHO_ESPECIES.to_string(), "al azar", "2", "1"),
            ("a=ACGT\nb=ACGA\nc=ACCA".to_string(), "al azar", "1", "9"),
            (largo, "al azar", "5", "200"),
        ] {
            let valores = entrada(json!({
                "secuencias": secuencias,
                "inicio": inicio,
                "semilla": semilla,
                "movimientos": movimientos
            }));
            let (_, (arbol, hechos, puntaje)) = narrar_parsimonia_grande(&valores).unwrap();
            let (nombres, hojas) = leer_secuencias(&valores, "secuencias", false).unwrap();
            let inicial = if inicio == "al azar" {
                filogenia::arbol_al_azar(&nombres, semilla.parse().unwrap()).unwrap()
            } else {
                let p = filogenia::matriz_de_distancias_p(&hojas).unwrap();
                filogenia::union_de_vecinos(&p, &nombres).unwrap()
            };
            let mut por_nodo = hojas;
            por_nodo.resize(inicial.nombres.len(), Vec::new());
            let esperado =
                filogenia::parsimonia_grande(&inicial, &por_nodo, movimientos.parse().unwrap())
                    .unwrap();
            assert_eq!(arbol, esperado.arbol);
            assert_eq!(hechos, esperado.movimientos);
            assert_eq!(puntaje, esperado.parsimonia.puntaje);
        }
    }

    #[test]
    fn el_mejor_vecino_de_la_tabla_coincide_con_el_movimiento() {
        let valores = entrada(de_ejemplo("filogenia.parsimonia_grande"));
        let (nombres, hojas) = leer_secuencias(&valores, "secuencias", false).unwrap();
        let inicial = filogenia::arbol_al_azar(&nombres, 7).unwrap();
        let mut por_nodo = hojas;
        por_nodo.resize(inicial.nombres.len(), Vec::new());
        let longitud = filogenia::revisar_arbol_binario(&inicial, &por_nodo).unwrap();
        let conjuntos = filogenia::conjuntos_de_letras(&por_nodo, longitud);
        let mut arbol = inicial.clone();
        let mut puntaje = filogenia::puntaje_de_fitch(&arbol, &conjuntos);
        let mut etiquetas = arbol.clone();
        nombres_unicos(&mut etiquetas, None);
        let mut rondas = 0;
        while let Some(movimiento) = filogenia::mejor_intercambio(&arbol, &conjuntos, puntaje) {
            let candidatos = vecinos_con_puntaje(&arbol, &conjuntos, &etiquetas.nombres);
            let menor = candidatos.iter().map(|fila| fila.4).min().unwrap();
            assert_eq!(menor, movimiento.puntaje);
            let (a, b) = movimiento.arista;
            let (x, y) = movimiento.intercambio;
            filogenia::intercambiar(&mut arbol, a, x, b, y);
            assert_eq!(
                filogenia::puntaje_de_fitch(&arbol, &conjuntos),
                movimiento.puntaje
            );
            puntaje = movimiento.puntaje;
            rondas += 1;
        }
        assert!(rondas > 0);
    }

    #[test]
    fn la_tabla_de_fitch_suma_el_puntaje_de_toda_la_secuencia() {
        let valores = entrada(de_ejemplo("filogenia.parsimonia_grande"));
        let (nombres, hojas) = leer_secuencias(&valores, "secuencias", false).unwrap();
        let inicial = filogenia::arbol_al_azar(&nombres, 7).unwrap();
        let mut por_nodo = hojas;
        por_nodo.resize(inicial.nombres.len(), Vec::new());
        let longitud = filogenia::revisar_arbol_binario(&inicial, &por_nodo).unwrap();
        let conjuntos = filogenia::conjuntos_de_letras(&por_nodo, longitud);
        let suma: usize = (0..longitud)
            .map(|posicion| fitch_de_posicion(&inicial, &conjuntos, posicion).unwrap().1)
            .sum();
        assert_eq!(suma, filogenia::puntaje_de_fitch(&inicial, &conjuntos));
    }
}
