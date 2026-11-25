use crate::formato::{miles, porcentaje};
use crate::pizarra::{Marca as Resalte, Narracion, Pizarra};
use crate::{Algoritmo, Arista, Bloque, Campo, Entrada, FilaAlineada, Marca, Paso, Salida, Tipo};
use nucleo::azar::Azar;
use nucleo::ensamblaje::{self, Grafo, ParDeLecturas};
use nucleo::{ecoli, secuencia};
use serde_json::{Value, json};
use std::collections::{HashMap, HashSet};

const MAXIMO_KMEROS: usize = 2_000_000;
const MAXIMO_LARGO_KMERO: usize = 1000;
const MAXIMO_SOLAPAMIENTO: usize = 100_000;
const MAXIMO_PARES: usize = 20_000;
const MAXIMO_D: usize = 100_000;
const MAXIMO_ARISTAS_DE_GRAFO: usize = 5_000_000;
const MAXIMO_CONTIGS: usize = 1_000_000;
const MAXIMO_NODOS_DIBUJO: usize = 60;
const MAXIMO_ARISTAS_DIBUJO: usize = 600;
const MAXIMO_LINEAS: usize = 1000;
const MAXIMO_DESTINOS_POR_LINEA: usize = 20;
const MAXIMO_ETIQUETA_EN_LISTA: usize = 100;
const MAXIMO_TEXTO_VISIBLE: usize = 10_000;
const MAXIMO_FILAS: usize = 2000;
const MAXIMO_FILAS_ESCALONADAS: usize = 60;
const MAXIMO_LARGO_ESCALONADO: usize = 300;
const MAXIMO_BARRAS: usize = 60;
const LARGO_MINIMO_BUSQUEDA: usize = 20;
const LARGO_MAXIMO_BUSQUEDA: usize = 100_000;
const MAXIMO_ARISTAS_ALTERNATIVAS: usize = 200_000;
const INTENTOS_ALTERNATIVOS: usize = 20;
const MAXIMO_KMEROS_ANIMADOS: usize = 60;
const MAXIMO_K_ANIMADO: usize = 40;
const MAXIMO_PASOS: usize = 400;
const LARGO_MINIMO_GENBANK: usize = 200;
const TOPE_DE_VECES: usize = 10;

const ORDENES: &[&str] = &["alfabético", "de aparición"];
const COPIAS: &[&str] = &["juntar las iguales", "una arista por copia"];

const EJEMPLO_CAMINO: &str = "ATGACCATGA
TGACCATGAT
GACCATGATT
ACCATGATTA
CCATGATTAC
CATGATTACG
ATGATTACGG
TGATTACGGA
GATTACGGAT
ATTACGGATT
TTACGGATTC
TACGGATTCA
ACGGATTCAC
CGGATTCACT
GGATTCACTG
GATTCACTGG
ATTCACTGGC
TTCACTGGCC
TCACTGGCCG
CACTGGCCGT
ACTGGCCGTC";
const EJEMPLO_SOLAPAMIENTO: &str = "AATTTA GTTCAA TGACTA CCACCA CACCAT ATCATC ATGACA ACCATC
GACACG TGACAC TTAAAC CATCAT ACCATC ATCACC ACACGC ATCCTG
TTTAAA TCCTGA CTGACT TAAACA AACACC GACTAG TTCAAT TCATCA
CATCAT CATCCT CCTGAC ACACCA CCATCA AAACAC GCGTTC CACCAT
TCACCA TCATCC CGTTCA ACGCGT CCATCA CAATTT CATCAC CACGCG
ACCACC CGCGTT CACCAC ATCATC TCAATT ATTTAA";
const EJEMPLO_DE_BRUIJN: &str = "GCATT TTTTT TACCG TGAAA TTTAC TTCCC ATGAA CGTTT
CCCCT TCTTC CGCAT AAACA CCGTT TCTTT ATTCT TTTCT
CATAC ACATA TTCGC TTTTT TACCT CACAT TTTTA ACCTT
CTTTT GTTTT TTTTC CCTGA CCCTG TTCTT GAAAC CTTCG
TCGCA TTTTT AACAC CTTCC TTCTT ACACA ACCGT TTACC
CCTTC ATACC TCCCC CATTC";
const EJEMPLO_CICLO: &str = "0 -> 1
1 -> 2,5
2 -> 3
3 -> 1,4
4 -> 0
5 -> 6
6 -> 3,7
7 -> 6";
const EJEMPLO_CAMINO_EULERIANO: &str = "1 -> 2
2 -> 3,6
3 -> 4
4 -> 2
6 -> 7
7 -> 8,9
8 -> 7";
const EJEMPLO_RECONSTRUCCION: &str =
    "CGCATTAGCA TTACCACCAC CCATTACCAC TCACCATTAC GGTAACGGTG ATGAAACGCA
CACCATTACC GTGCGGGCTG CATTAGCACC CCACCATCAC AGGTAACGGT CCACAGGTAA
ACCACAGGTA ACCACCACCA AAACGCATTA TAACGGTGCG CCATTACCAC AGCACCACCA
CCACCACCAT CACCACCATT CACCATCACC TAGCACCACC ACGCATTAGC TTAGCACCAC
AACGGTGCGG ACCATCACCA CCACCATTAC ACAGGTAACG CCATCACCAT ATTACCACCA
ACCACCATTA TGAAACGCAT ATTAGCACCA AACGCATTAG CACCATTACC ACCATTACCA
GAAACGCATT ATCACCATTA TTACCACAGG TGCGGGCTGA GTAACGGTGC GCATTAGCAC
CAGGTAACGG ACGGTGCGGG CATTACCACC ACCACCATCA ACCATTACCA CATTACCACA
GCACCACCAT CACCACCATC CACAGGTAAC ATTACCACAG GGTGCGGGCT TACCACCACC
CGGTGCGGGC TACCACAGGT CATCACCATT";
const EJEMPLO_HUECOS: &str = "ATGAA|AATTT TGAAA|ATTTT GAAAG|TTTTC AAAGC|TTTCG
AAGCA|TTCGT AGCAA|TCGTA GCAAT|CGTAC CAATT|GTACT
AATTT|TACTG ATTTT|ACTGA TTTTC|CTGAA TTTCG|TGAAA
TTCGT|GAAAG TCGTA|AAAGG CGTAC|AAGGT GTACT|AGGTT
TACTG|GGTTG ACTGA|GTTGG CTGAA|TTGGT TGAAA|TGGTG
GAAAG|GGTGG AAAGG|GTGGC AAGGT|TGGCG AGGTT|GGCGC
GGTTG|GCGCA GTTGG|CGCAC TTGGT|GCACT TGGTG|CACTT
GGTGG|ACTTC GTGGC|CTTCC TGGCG|TTCCT GGCGC|TCCTG
GCGCA|CCTGA";
const EJEMPLO_PARES: &str = "TCACCA|CCACAG CGCATT|ACCACC ACCACC|ACCACC ACCATT|ACCACC
TACCAC|TAACGG TTACCA|CCATCA GGTAAC|GCGGGC CATTAG|CACCAT
TAACGG|GGGCTG ACCATC|ATTACC CCACCA|CCACCA CCATTA|CCACCA
GCACCA|TTACCA GAAACG|TAGCAC CACAGG|CGGTGC TTACCA|GTAACG
ACCACA|AACGGT CCATTA|CAGGTA CCACAG|ACGGTG GCATTA|CCACCA
CACCAC|CACCAT TACCAC|CATCAC ATTAGC|ACCATT CACCAT|CACAGG
CACCAT|CATTAC TGAAAC|TTAGCA ACCACC|ACCATT CCACCA|TCACCA
AGCACC|ATTACC TTAGCA|CCATTA AAACGC|AGCACC CATTAC|AGGTAA
ACCACC|ATCACC GTAACG|CGGGCT CAGGTA|GTGCGG CCATCA|TTACCA
AACGGT|GGCTGA CATCAC|TACCAC ATCACC|ACCACA ACCATT|ACAGGT
ATTACC|GGTAAC AACGCA|GCACCA CACCAT|CACCAC ACGCAT|CACCAC
CACCAC|TACCAC AGGTAA|TGCGGG TAGCAC|CATTAC ATTACC|ACCATC
CCACCA|CCATTA ACAGGT|GGTGCG ATGAAA|ATTAGC CATTAC|CACCAT";
const EJEMPLO_NO_RAMIFICADOS: &str = "1 -> 2
2 -> 3
3 -> 4,7
4 -> 5
5 -> 6
7 -> 6
6 -> 8
9 -> 10
10 -> 11
11 -> 9";
const EJEMPLO_CONTIGS: &str =
    "TTAGCGTG GCGTGGTG TCTACGAG AGCCCTTC GGGGCTGC AAAGGCTT ACTTGGAC TAGCCTGG
ACGAGTGA ATTATCCC CGTGATTA TGATTATT CCCTTCTA ATGACAGC AAGGCTTA CGTGGTGG
AGCGTGGT GATTAGCG GGTGGTGA CTTCTACG TACGAGTG TGGTCGTG CGTGCGGG TTATCCCA
CTACGAGT TAGCGTGG TTATTATC ACGAGGAA GACGAGGA AGGCTTAG GAGGAAAG TGCGGGGC
CGAGGAAA CACCGTGC CCTTCTAC GCCCTTCT TTAGCCTG TGCACTTG GTGATTAT CGAGTGAT
TTGGACGA GTGGTGAT CCTGGTCG CCACCGTG CTGGTCGT AGGAAAGG ATCCCACC GGACGAGG
CAGCCCTT TATTATCC ACAGCCCT TCGTGATT GCGGGGCT CTTGGACG ATTAGCGT ATTATTAT
GCTGCACT GGTGATTA GTGCGGGG TGACAGCC GAGTGATT ATTAGCCT TGGACGAG GTGATTAG
GGAAAGGC TTCTACGA GAAAGGCT CTGCACTT CGGGGCTG TGATTAGC GTGGTGGT GATTAGCC
CCCACCGT AGCCTGGT GACAGCCC GCACTTGG TGGTGATT CCGTGCGG CACTTGGA TGATTAGC
GATTATTA TGGTGGTG TCCCACCG AGTGATTA GTGATTAG ACCGTGCG GTCGTGAT GGCTGCAC
GCCTGGTC TATCCCAC GGTCGTGA GGGCTGCA";
const EJEMPLO_N50: &str = "870 226 60 515 640 213 69 1380 362 1526
1580 183 592 17 708 149 14 99 4166 15
1455 614 1839 1282 597 56 37 4478 143 3435
74 259 765 269 14 980 779";

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![
        Algoritmo {
            id: "ensamblaje.composicion",
            capitulo: "ensamblaje",
            nombre: "Composición en k-meros",
            resumen: "Todos los k-meros de un texto, como las lecturas perfectas que entregaría una máquina secuenciadora.",
            campos: &[
                Campo {
                    id: "texto",
                    etiqueta: "Texto",
                    tipo: Tipo::Secuencia,
                    ejemplo: "@gen:leuL",
                },
                Campo {
                    id: "k",
                    etiqueta: "k (largo de cada k-mero)",
                    tipo: Tipo::Entero,
                    ejemplo: "5",
                },
                Campo {
                    id: "orden",
                    etiqueta: "Orden de la lista",
                    tipo: Tipo::Opcion(ORDENES),
                    ejemplo: "alfabético",
                },
            ],
            funcion: "ensamblaje::composicion",
            ejecutar: composicion,
            pasos: Some(pasos_composicion),
        },
        Algoritmo {
            id: "ensamblaje.camino_de_kmeros",
            capitulo: "ensamblaje",
            nombre: "Cadena deletreada por un camino",
            resumen: "Une k-meros consecutivos, donde cada uno empieza una letra después del anterior, en una sola cadena.",
            campos: &[Campo {
                id: "kmeros",
                etiqueta: "k-meros en orden",
                tipo: Tipo::Secuencias,
                ejemplo: EJEMPLO_CAMINO,
            }],
            funcion: "ensamblaje::deletrear_camino",
            ejecutar: camino_de_kmeros,
            pasos: Some(pasos_camino_de_kmeros),
        },
        Algoritmo {
            id: "ensamblaje.solapamiento",
            capitulo: "ensamblaje",
            nombre: "Grafo de solapamiento",
            resumen: "Une cada k-mero con los que empiezan con sus últimas k − 1 letras.",
            campos: &[Campo {
                id: "kmeros",
                etiqueta: "k-meros",
                tipo: Tipo::Secuencias,
                ejemplo: EJEMPLO_SOLAPAMIENTO,
            }],
            funcion: "ensamblaje::grafo_de_solapamiento",
            ejecutar: solapamiento,
            pasos: Some(pasos_solapamiento),
        },
        Algoritmo {
            id: "ensamblaje.de_bruijn_texto",
            capitulo: "ensamblaje",
            nombre: "Grafo de De Bruijn de un texto",
            resumen: "Cada k-mero del texto es una arista entre su prefijo y su sufijo, y los (k − 1)-meros iguales se pegan en un solo nodo.",
            campos: &[
                Campo {
                    id: "texto",
                    etiqueta: "Texto",
                    tipo: Tipo::Secuencia,
                    ejemplo: "@gen:thrL",
                },
                Campo {
                    id: "k",
                    etiqueta: "k (largo de cada k-mero)",
                    tipo: Tipo::Entero,
                    ejemplo: "4",
                },
            ],
            funcion: "ensamblaje::de_bruijn_de_texto",
            ejecutar: de_bruijn_texto,
            pasos: Some(pasos_de_bruijn_texto),
        },
        Algoritmo {
            id: "ensamblaje.de_bruijn_kmeros",
            capitulo: "ensamblaje",
            nombre: "Grafo de De Bruijn de k-meros",
            resumen: "El mismo grafo, armado directamente desde una lista desordenada de k-meros.",
            campos: &[Campo {
                id: "kmeros",
                etiqueta: "k-meros",
                tipo: Tipo::Secuencias,
                ejemplo: EJEMPLO_DE_BRUIJN,
            }],
            funcion: "ensamblaje::de_bruijn_de_kmeros",
            ejecutar: de_bruijn_kmeros,
            pasos: Some(pasos_de_bruijn_kmeros),
        },
        Algoritmo {
            id: "ensamblaje.ciclo_euleriano",
            capitulo: "ensamblaje",
            nombre: "Ciclo euleriano",
            resumen: "Un recorrido que usa cada arista una sola vez y vuelve al punto de partida, con el algoritmo de Hierholzer.",
            campos: &[Campo {
                id: "grafo",
                etiqueta: "Grafo (una línea por nodo: origen -> destino1,destino2)",
                tipo: Tipo::Texto,
                ejemplo: EJEMPLO_CICLO,
            }],
            funcion: "ensamblaje::ciclo_euleriano",
            ejecutar: ciclo,
            pasos: Some(pasos_ciclo),
        },
        Algoritmo {
            id: "ensamblaje.camino_euleriano",
            capitulo: "ensamblaje",
            nombre: "Camino euleriano",
            resumen: "Un recorrido que usa cada arista una sola vez, de un nodo inicial a uno final.",
            campos: &[Campo {
                id: "grafo",
                etiqueta: "Grafo (una línea por nodo: origen -> destino1,destino2)",
                tipo: Tipo::Texto,
                ejemplo: EJEMPLO_CAMINO_EULERIANO,
            }],
            funcion: "ensamblaje::camino_euleriano",
            ejecutar: camino,
            pasos: Some(pasos_camino),
        },
        Algoritmo {
            id: "ensamblaje.reconstruccion",
            capitulo: "ensamblaje",
            nombre: "Reconstruir una cadena desde sus k-meros",
            resumen: "Arma el grafo de De Bruijn de los k-meros y recorre un camino euleriano para deletrear la cadena.",
            campos: &[Campo {
                id: "kmeros",
                etiqueta: "k-meros (en cualquier orden)",
                tipo: Tipo::Secuencias,
                ejemplo: EJEMPLO_RECONSTRUCCION,
            }],
            funcion: "ensamblaje::reconstruir_desde_kmeros",
            ejecutar: reconstruccion,
            pasos: Some(pasos_reconstruccion),
        },
        Algoritmo {
            id: "ensamblaje.universal",
            capitulo: "ensamblaje",
            nombre: "Cadena k-universal circular",
            resumen: "La cadena binaria circular más corta que contiene cada k-mero de ceros y unos exactamente una vez.",
            campos: &[Campo {
                id: "k",
                etiqueta: "k (de 1 a 20)",
                tipo: Tipo::Entero,
                ejemplo: "4",
            }],
            funcion: "ensamblaje::cadena_universal_circular",
            ejecutar: universal,
            pasos: Some(pasos_universal),
        },
        Algoritmo {
            id: "ensamblaje.pares",
            capitulo: "ensamblaje",
            nombre: "Reconstruir desde pares de lecturas",
            resumen: "Grafo de De Bruijn pareado y un camino euleriano cuyas dos mitades encajen a la distancia d.",
            campos: &[
                Campo {
                    id: "pares",
                    etiqueta: "Pares en cualquier orden (primera|segunda)",
                    tipo: Tipo::Texto,
                    ejemplo: EJEMPLO_PARES,
                },
                Campo {
                    id: "d",
                    etiqueta: "d (letras entre las dos lecturas)",
                    tipo: Tipo::Entero,
                    ejemplo: "3",
                },
                Campo {
                    id: "semilla",
                    etiqueta: "Semilla (para barajar si hay que reintentar)",
                    tipo: Tipo::Entero,
                    ejemplo: "1",
                },
            ],
            funcion: "ensamblaje::reconstruir_desde_pares",
            ejecutar: pares,
            pasos: Some(pasos_pares),
        },
        Algoritmo {
            id: "ensamblaje.contigs",
            capitulo: "ensamblaje",
            nombre: "Contigs",
            resumen: "Los trozos de genoma que se pueden armar sin dudas: caminos no ramificados del grafo de De Bruijn.",
            campos: &[
                Campo {
                    id: "kmeros",
                    etiqueta: "k-meros",
                    tipo: Tipo::Secuencias,
                    ejemplo: EJEMPLO_CONTIGS,
                },
                Campo {
                    id: "copias",
                    etiqueta: "k-meros repetidos",
                    tipo: Tipo::Opcion(COPIAS),
                    ejemplo: "juntar las iguales",
                },
            ],
            funcion: "ensamblaje::generar_contigs",
            ejecutar: contigs,
            pasos: Some(pasos_contigs),
        },
        Algoritmo {
            id: "ensamblaje.camino_con_huecos",
            capitulo: "ensamblaje",
            nombre: "Cadena deletreada por un camino con huecos",
            resumen: "Une pares de lecturas consecutivos, separados por d letras desconocidas, en una sola cadena.",
            campos: &[
                Campo {
                    id: "pares",
                    etiqueta: "Pares en orden (primera|segunda)",
                    tipo: Tipo::Texto,
                    ejemplo: EJEMPLO_HUECOS,
                },
                Campo {
                    id: "d",
                    etiqueta: "d (letras entre las dos lecturas)",
                    tipo: Tipo::Entero,
                    ejemplo: "3",
                },
            ],
            funcion: "ensamblaje::deletrear_camino_con_huecos",
            ejecutar: camino_con_huecos,
            pasos: Some(pasos_camino_con_huecos),
        },
        Algoritmo {
            id: "ensamblaje.no_ramificados",
            capitulo: "ensamblaje",
            nombre: "Caminos máximos no ramificados",
            resumen: "Los tramos más largos de un grafo en los que cada nodo interior tiene una sola entrada y una sola salida.",
            campos: &[Campo {
                id: "grafo",
                etiqueta: "Grafo (una línea por nodo: origen -> destino1,destino2)",
                tipo: Tipo::Texto,
                ejemplo: EJEMPLO_NO_RAMIFICADOS,
            }],
            funcion: "ensamblaje::caminos_no_ramificados",
            ejecutar: no_ramificados,
            pasos: Some(pasos_no_ramificados),
        },
        Algoritmo {
            id: "ensamblaje.n50",
            capitulo: "ensamblaje",
            nombre: "N50 de un ensamblaje",
            resumen: "El largo del contig más corto entre los más largos que, juntos, suman la mitad del ensamblaje.",
            campos: &[Campo {
                id: "contigs",
                etiqueta: "Contigs (sus secuencias o solo sus largos)",
                tipo: Tipo::Secuencias,
                ejemplo: EJEMPLO_N50,
            }],
            funcion: "ensamblaje::medir_n50",
            ejecutar: n50,
            pasos: Some(pasos_n50),
        },
    ]
}

fn letras(secuencia: &[u8]) -> String {
    String::from_utf8_lossy(secuencia).into_owned()
}

fn recortar(secuencia: &[u8]) -> String {
    if secuencia.len() <= MAXIMO_TEXTO_VISIBLE {
        return letras(secuencia);
    }
    format!(
        "{}… (se muestran {} de {} letras)",
        letras(&secuencia[..MAXIMO_TEXTO_VISIBLE]),
        miles(MAXIMO_TEXTO_VISIBLE),
        miles(secuencia.len())
    )
}

fn etiqueta_corta(etiqueta: &[u8]) -> String {
    if etiqueta.len() <= 24 {
        return letras(etiqueta);
    }
    format!(
        "{}…{}",
        letras(&etiqueta[..10]),
        letras(&etiqueta[etiqueta.len() - 10..])
    )
}

fn etiqueta_de_lista(etiqueta: &[u8]) -> String {
    if etiqueta.len() <= MAXIMO_ETIQUETA_EN_LISTA {
        return letras(etiqueta);
    }
    format!("{}…", letras(&etiqueta[..MAXIMO_ETIQUETA_EN_LISTA]))
}

fn cantidad(numero: usize, singular: &str, plural: &str) -> String {
    if numero == 1 {
        format!("1 {singular}")
    } else {
        format!("{} {plural}", miles(numero))
    }
}

fn es_adn(secuencia: &[u8]) -> bool {
    secuencia
        .iter()
        .all(|letra| matches!(letra, b'A' | b'C' | b'G' | b'T'))
}

fn leer_kmeros(entrada: &Entrada, campo: &str) -> Result<Vec<Vec<u8>>, String> {
    let kmeros = entrada.secuencias(campo)?;
    if kmeros.len() > MAXIMO_KMEROS {
        return Err(format!(
            "Hay {} k-meros; aquí se admiten hasta {}.",
            miles(kmeros.len()),
            miles(MAXIMO_KMEROS)
        ));
    }
    for (numero, kmero) in kmeros.iter().enumerate() {
        if kmero.len() > MAXIMO_LARGO_KMERO {
            return Err(format!(
                "El k-mero {} mide {} letras; aquí cada k-mero puede tener hasta {}. Para partir una secuencia larga en k-meros usa «Composición en k-meros».",
                numero + 1,
                miles(kmero.len()),
                miles(MAXIMO_LARGO_KMERO)
            ));
        }
        if let Some(letra) = kmero.iter().find(|letra| !letra.is_ascii_alphanumeric()) {
            return Err(format!(
                "El k-mero {} («{}») tiene «{}», que no es una letra ni un dígito.",
                numero + 1,
                etiqueta_corta(kmero),
                *letra as char
            ));
        }
    }
    Ok(kmeros)
}

fn rebanadas(kmeros: &[Vec<u8>]) -> Vec<&[u8]> {
    kmeros.iter().map(|kmero| kmero.as_slice()).collect()
}

fn leer_k(entrada: &Entrada, campo: &str, minimo: usize, largo: usize) -> Result<usize, String> {
    let k = entrada.natural(campo)?;
    if k < minimo || k > largo {
        return Err(format!(
            "k debe estar entre {minimo} y el largo del texto ({}).",
            miles(largo)
        ));
    }
    Ok(k)
}

fn leer_d(entrada: &Entrada) -> Result<usize, String> {
    let d = entrada.natural("d")?;
    if d > MAXIMO_D {
        return Err(format!("d puede ser como mucho {}.", miles(MAXIMO_D)));
    }
    Ok(d)
}

fn leer_grafo(entrada: &Entrada, campo: &str) -> Result<Grafo, String> {
    let texto = entrada.texto(campo)?;
    let grafo = ensamblaje::leer_lista_de_adyacencia(&texto)?;
    if ensamblaje::cantidad_de_aristas(&grafo) > MAXIMO_ARISTAS_DE_GRAFO {
        return Err(format!(
            "El grafo tiene más de {} aristas.",
            miles(MAXIMO_ARISTAS_DE_GRAFO)
        ));
    }
    Ok(grafo)
}

fn leer_pares(entrada: &Entrada, campo: &str) -> Result<Vec<ParDeLecturas>, String> {
    let texto = entrada.texto(campo)?;
    let pares = ensamblaje::leer_pares(&texto)?;
    if pares.len() > MAXIMO_PARES {
        return Err(format!(
            "Hay {} pares; aquí se admiten hasta {}.",
            miles(pares.len()),
            miles(MAXIMO_PARES)
        ));
    }
    if pares
        .iter()
        .any(|par| par.primera.len() > MAXIMO_LARGO_KMERO || par.segunda.len() > MAXIMO_LARGO_KMERO)
    {
        return Err(format!(
            "Cada lectura puede tener hasta {} letras.",
            miles(MAXIMO_LARGO_KMERO)
        ));
    }
    Ok(pares)
}

fn ramificados(grafo: &Grafo, entradas: &[usize]) -> Vec<usize> {
    (0..grafo.salidas.len())
        .filter(|&nodo| entradas[nodo] > 1 || grafo.salidas[nodo].len() > 1)
        .collect()
}

fn se_puede_dibujar(grafo: &Grafo) -> bool {
    grafo.nodos.len() <= MAXIMO_NODOS_DIBUJO
        && ensamblaje::cantidad_de_aristas(grafo) <= MAXIMO_ARISTAS_DIBUJO
}

fn dibujo(etiqueta: &str, grafo: &Grafo, aristas: Vec<(usize, usize, String)>) -> Bloque {
    let mut orden: Vec<(usize, usize)> = Vec::new();
    let mut nombres: HashMap<(usize, usize), Vec<String>> = HashMap::new();
    for (desde, hasta, nombre) in aristas {
        let lista = nombres.entry((desde, hasta)).or_default();
        if lista.is_empty() {
            orden.push((desde, hasta));
        }
        lista.push(nombre);
    }
    let mut dibujadas = Vec::with_capacity(orden.len());
    for par in orden {
        let lista = &nombres[&par];
        let todas_iguales = lista.iter().all(|nombre| *nombre == lista[0]);
        let texto = if lista.len() == 1 {
            lista[0].clone()
        } else if todas_iguales {
            format!("{} ×{}", lista[0], lista.len()).trim().to_string()
        } else {
            lista.join(", ")
        };
        dibujadas.push(Arista {
            desde: par.0,
            hasta: par.1,
            etiqueta: texto,
        });
    }
    let nodos = grafo
        .nodos
        .iter()
        .map(|nodo| etiqueta_corta(nodo))
        .collect();
    Bloque::grafo(etiqueta, nodos, dibujadas)
}

fn nombre_de_arista(desde: &[u8], hasta: &[u8]) -> String {
    let separar = |etiqueta: &[u8]| -> Vec<Vec<u8>> {
        etiqueta
            .split(|&letra| letra == b'|')
            .map(|parte| parte.to_vec())
            .collect()
    };
    let partes_desde = separar(desde);
    let partes_hasta = separar(hasta);
    if partes_desde.len() != partes_hasta.len() {
        return String::new();
    }
    let mut partes = Vec::new();
    for (inicio, fin) in partes_desde.iter().zip(&partes_hasta) {
        let mut parte = inicio.clone();
        if let Some(&ultima) = fin.last() {
            parte.push(ultima);
        }
        partes.push(etiqueta_corta(&parte));
    }
    partes.join("|")
}

fn aristas_de_bruijn(grafo: &Grafo) -> Vec<(usize, usize, String)> {
    let mut aristas = Vec::new();
    for (desde, destinos) in grafo.salidas.iter().enumerate() {
        for &hasta in destinos {
            let nombre = nombre_de_arista(&grafo.nodos[desde], &grafo.nodos[hasta]);
            aristas.push((desde, hasta, nombre));
        }
    }
    aristas
}

fn aristas_sin_nombre(grafo: &Grafo) -> Vec<(usize, usize, String)> {
    let mut aristas = Vec::new();
    for (desde, destinos) in grafo.salidas.iter().enumerate() {
        for &hasta in destinos {
            aristas.push((desde, hasta, String::new()));
        }
    }
    aristas
}

fn aristas_numeradas(camino: &[usize]) -> Vec<(usize, usize, String)> {
    camino
        .windows(2)
        .enumerate()
        .map(|(orden, par)| (par[0], par[1], (orden + 1).to_string()))
        .collect()
}

fn lista_de_adyacencia(etiqueta: &str, grafo: &Grafo) -> Bloque {
    let mut lineas = Vec::new();
    let mut total = 0;
    for (desde, destinos) in grafo.salidas.iter().enumerate() {
        if destinos.is_empty() {
            continue;
        }
        total += 1;
        if lineas.len() >= MAXIMO_LINEAS {
            continue;
        }
        let mut nombres: Vec<String> = destinos
            .iter()
            .take(MAXIMO_DESTINOS_POR_LINEA)
            .map(|&hasta| etiqueta_de_lista(&grafo.nodos[hasta]))
            .collect();
        if destinos.len() > MAXIMO_DESTINOS_POR_LINEA {
            nombres.push(format!(
                "… y {} más",
                miles(destinos.len() - MAXIMO_DESTINOS_POR_LINEA)
            ));
        }
        lineas.push(format!(
            "{} -> {}",
            etiqueta_de_lista(&grafo.nodos[desde]),
            nombres.join(",")
        ));
    }
    Bloque::Lista {
        etiqueta: etiqueta.to_string(),
        elementos: lineas,
        total,
    }
}

fn tabla_de_nodos(etiqueta: &str, grafo: &Grafo, entradas: &[usize], nodos: &[usize]) -> Bloque {
    let filas = nodos
        .iter()
        .take(MAXIMO_FILAS)
        .map(|&nodo| {
            vec![
                etiqueta_de_lista(&grafo.nodos[nodo]),
                miles(entradas[nodo]),
                miles(grafo.salidas[nodo].len()),
            ]
        })
        .collect();
    Bloque::Tabla {
        etiqueta: etiqueta.to_string(),
        columnas: vec![
            "Nodo".to_string(),
            "Entran".to_string(),
            "Salen".to_string(),
        ],
        filas,
        total: nodos.len(),
    }
}

fn camino_como_texto(grafo: &Grafo, camino: &[usize]) -> String {
    let mut partes: Vec<String> = camino
        .iter()
        .take(MAXIMO_FILAS)
        .map(|&nodo| etiqueta_de_lista(&grafo.nodos[nodo]))
        .collect();
    if camino.len() > MAXIMO_FILAS {
        partes.push(format!("… ({} nodos en total)", miles(camino.len())));
    }
    partes.join("->")
}

fn fila_escalonada(nombre: String, desde: usize, contenido: &[u8], largo: usize) -> FilaAlineada {
    let mut fila = "-".repeat(desde);
    fila.push_str(&letras(contenido));
    let resto = largo.saturating_sub(fila.len());
    fila.push_str(&"-".repeat(resto));
    FilaAlineada {
        nombre,
        secuencia: fila,
    }
}

struct Hallazgo {
    inicio: usize,
    fin: usize,
    directa: bool,
    veces: usize,
}

fn buscar_en_ecoli(cadena: &[u8]) -> Option<Hallazgo> {
    if cadena.len() < LARGO_MINIMO_BUSQUEDA
        || cadena.len() > LARGO_MAXIMO_BUSQUEDA
        || !es_adn(cadena)
    {
        return None;
    }
    let inversa = secuencia::complemento_inverso(cadena);
    let mut primero = None;
    let mut veces = 0;
    for (posicion, ventana) in ecoli::genoma().windows(cadena.len()).enumerate() {
        let directa = ventana[0] == cadena[0] && ventana == cadena;
        let opuesta = ventana[0] == inversa[0] && ventana == inversa.as_slice();
        if !directa && !opuesta {
            continue;
        }
        veces += 1;
        if primero.is_none() {
            primero = Some((posicion, directa));
        }
    }
    let (posicion, directa) = primero?;
    Some(Hallazgo {
        inicio: posicion + 1,
        fin: posicion + cadena.len(),
        directa,
        veces,
    })
}

fn explicar_ubicacion(salida: &mut Salida, sujeto: &str, cadena: &[u8], si_no_esta: &str) {
    if cadena.len() < LARGO_MINIMO_BUSQUEDA
        || cadena.len() > LARGO_MAXIMO_BUSQUEDA
        || !es_adn(cadena)
    {
        return;
    }
    let Some(hallazgo) = buscar_en_ecoli(cadena) else {
        if !si_no_esta.is_empty() {
            salida.explicar(si_no_esta);
        }
        return;
    };
    let hebra = if hallazgo.directa {
        "directa"
    } else {
        "complementaria"
    };
    let mut texto = format!(
        "{sujeto} está en el genoma de E. coli, en las posiciones {} a {} de la hebra {hebra} (contando desde 1, como en la anotación de NC_000913.3)",
        miles(hallazgo.inicio),
        miles(hallazgo.fin)
    );
    let dentro = ecoli::genes()
        .iter()
        .find(|anotado| anotado.inicio <= hallazgo.inicio && anotado.fin >= hallazgo.fin);
    let toca = ecoli::genes()
        .iter()
        .find(|anotado| anotado.inicio <= hallazgo.fin && anotado.fin >= hallazgo.inicio);
    if let Some(anotado) = dentro {
        texto.push_str(&format!(
            ", dentro del gen {} ({})",
            anotado.nombre, anotado.producto
        ));
    } else if let Some(anotado) = toca {
        texto.push_str(&format!(
            ", y toca el gen {} ({})",
            anotado.nombre, anotado.producto
        ));
    }
    texto.push('.');
    if hallazgo.veces > 1 {
        texto.push_str(&format!(
            " Aparece {} veces en el genoma: es una repetición.",
            miles(hallazgo.veces)
        ));
    }
    salida.explicar(texto);
}

fn contar_iguales<'a>(ordenados: &[&'a [u8]]) -> Vec<(&'a [u8], usize)> {
    let mut grupos: Vec<(&[u8], usize)> = Vec::new();
    for &pieza in ordenados {
        match grupos.last_mut() {
            Some((anterior, veces)) if *anterior == pieza => *veces += 1,
            _ => grupos.push((pieza, 1)),
        }
    }
    grupos
}

fn cuantos_se_repiten(texto: &[u8], largo: usize) -> usize {
    let mut ventanas = ensamblaje::composicion(texto, largo);
    ventanas.sort_unstable();
    contar_iguales(&ventanas)
        .iter()
        .filter(|(_, veces)| *veces > 1)
        .count()
}

fn composicion(entrada: &Entrada) -> Result<Salida, String> {
    let texto = entrada.secuencia("texto")?;
    let k = leer_k(entrada, "k", 1, texto.len())?;
    let orden = entrada.opcion("orden", ORDENES)?;
    let en_orden = ensamblaje::composicion(&texto, k);
    let mut alfabeticos = en_orden.clone();
    alfabeticos.sort_unstable();
    let grupos = contar_iguales(&alfabeticos);
    let mut mas_repetido: (&[u8], usize) = (&[], 0);
    for &(kmero, veces) in &grupos {
        if veces > mas_repetido.1 {
            mas_repetido = (kmero, veces);
        }
    }
    let repetidos = grupos.iter().filter(|(_, veces)| *veces > 1).count();
    let tope = mas_repetido.1.min(TOPE_DE_VECES);
    let mut por_veces = vec![0.0; tope];
    for &(_, veces) in &grupos {
        por_veces[veces.min(TOPE_DE_VECES) - 1] += 1.0;
    }
    let categorias: Vec<String> = (1..=tope)
        .map(|veces| {
            if veces == TOPE_DE_VECES && mas_repetido.1 > TOPE_DE_VECES {
                format!("{TOPE_DE_VECES} o más")
            } else {
                veces.to_string()
            }
        })
        .collect();
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("k-meros", miles(en_orden.len())));
    salida.agregar(Bloque::valor("k-meros distintos", miles(grupos.len())));
    let lista = if orden == ORDENES[0] {
        &alfabeticos
    } else {
        &en_orden
    };
    salida.agregar(Bloque::lista(
        &format!("Composición (orden {orden})"),
        lista.iter().map(|kmero| String::from_utf8_lossy(kmero)),
    ));
    salida.agregar(Bloque::barras(
        "Cuántas veces aparece cada k-mero",
        "Veces que aparece",
        "k-meros distintos",
        categorias,
        por_veces,
    ));
    salida.explicar(format!(
        "Un texto de {} letras tiene {} − {k} + 1 = {} k-meros: uno por cada posición donde cabe una ventana de {k} letras.",
        miles(texto.len()),
        miles(texto.len()),
        miles(en_orden.len())
    ));
    if repetidos > 0 {
        salida.explicar(format!(
            "{} más de una vez; el más repetido, «{}», aparece {} veces. En la composición no queda dicho de qué parte del texto viene cada copia.",
            cantidad(repetidos, "k-mero aparece", "k-meros aparecen"),
            etiqueta_corta(mas_repetido.0),
            miles(mas_repetido.1)
        ));
    } else {
        salida.explicar("Ningún k-mero se repite.");
    }
    if k >= 2 {
        let repetidos_cortos = cuantos_se_repiten(&texto, k - 1);
        if repetidos_cortos == 0 {
            salida.explicar(
                "Tampoco se repite ningún (k − 1)-mero, así que solo hay una cadena con esta composición: el texto se puede reconstruir sin dudas.",
            );
        } else {
            salida.explicar(format!(
                "{} más de una vez: en el grafo de De Bruijn serán nodos con varias entradas o salidas, y puede haber otras cadenas con esta misma composición.",
                cantidad(
                    repetidos_cortos,
                    "trozo de largo k − 1 aparece",
                    "trozos de largo k − 1 aparecen"
                )
            ));
        }
    }
    if orden == ORDENES[0] {
        salida.explicar(
            "En orden alfabético ya no se sabe dónde estaba cada k-mero: así llegan las lecturas de una máquina secuenciadora, y ensamblar es recuperar ese orden.",
        );
    } else {
        salida.explicar(format!(
            "En orden de aparición, cada k-mero empieza una letra después del anterior y comparte con él {}: el texto es el camino que los une en fila.",
            cantidad(k - 1, "letra", "letras")
        ));
    }
    if es_adn(&texto) && k < 32 && 4u64.pow(k as u32) < en_orden.len() as u64 {
        salida.explicar(format!(
            "Con k = {k} solo existen 4^{k} = {} k-meros de ADN distintos, menos que las {} ventanas: por fuerza se repiten. Usa un k más grande para distinguir posiciones.",
            miles(4usize.pow(k as u32)),
            miles(en_orden.len())
        ));
    }
    Ok(salida)
}

fn camino_de_kmeros(entrada: &Entrada) -> Result<Salida, String> {
    let kmeros = leer_kmeros(entrada, "kmeros")?;
    let lista = rebanadas(&kmeros);
    let cadena = ensamblaje::deletrear_camino(&lista)?;
    let k = lista[0].len();
    let mut salida = Salida::default();
    salida.agregar(Bloque::texto("Cadena", recortar(&cadena)));
    salida.agregar(Bloque::valor(
        "Largo",
        format!("{} letras", miles(cadena.len())),
    ));
    if lista.len() <= MAXIMO_FILAS_ESCALONADAS && cadena.len() <= MAXIMO_LARGO_ESCALONADO {
        let mut filas = vec![fila_escalonada(
            "cadena".to_string(),
            0,
            &cadena,
            cadena.len(),
        )];
        for (posicion, kmero) in lista.iter().enumerate() {
            filas.push(fila_escalonada(
                format!("k-mero {}", posicion + 1),
                posicion,
                kmero,
                cadena.len(),
            ));
        }
        salida.agregar(Bloque::alineamiento("Cada k-mero en su lugar", filas));
    }
    if k == 1 {
        salida.explicar(
            "Con k = 1 no hay solapamiento: cada letra se pone a continuación de la anterior.",
        );
    } else {
        salida.explicar(format!(
            "Cada k-mero agrega solo su última letra, porque {} al final del anterior: {} de largo {k} {} {} + {k} − 1 = {} letras.",
            if k == 2 {
                "la otra letra ya estaba".to_string()
            } else {
                format!("las otras {} ya estaban", k - 1)
            },
            cantidad(lista.len(), "k-mero", "k-meros"),
            if lista.len() == 1 { "deletrea" } else { "deletrean" },
            miles(lista.len()),
            miles(cadena.len())
        ));
    }
    explicar_ubicacion(
        &mut salida,
        "Esta cadena",
        &cadena,
        "No aparece tal cual en el genoma de E. coli, en ninguna de sus dos hebras.",
    );
    Ok(salida)
}

fn solapamiento(entrada: &Entrada) -> Result<Salida, String> {
    let kmeros = leer_kmeros(entrada, "kmeros")?;
    if kmeros.len() > MAXIMO_SOLAPAMIENTO {
        return Err(format!(
            "El grafo de solapamiento admite hasta {} k-meros; aquí hay {}.",
            miles(MAXIMO_SOLAPAMIENTO),
            miles(kmeros.len())
        ));
    }
    let lista = rebanadas(&kmeros);
    let grafo = ensamblaje::grafo_de_solapamiento(&lista)?;
    let k = lista[0].len();
    let nodos = grafo.nodos.len();
    let aristas = ensamblaje::cantidad_de_aristas(&grafo);
    let entradas = ensamblaje::grados_de_entrada(&grafo);
    let sin_entrada: Vec<usize> = (0..nodos).filter(|&nodo| entradas[nodo] == 0).collect();
    let sin_salida = (0..nodos)
        .filter(|&nodo| grafo.salidas[nodo].is_empty())
        .count();
    let varias_salidas = (0..nodos)
        .filter(|&nodo| grafo.salidas[nodo].len() > 1)
        .count();
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Nodos (lecturas)", miles(nodos)));
    salida.agregar(Bloque::valor("Aristas (solapamientos)", miles(aristas)));
    if se_puede_dibujar(&grafo) {
        salida.agregar(dibujo(
            "Grafo de solapamiento",
            &grafo,
            aristas_sin_nombre(&grafo),
        ));
    }
    salida.agregar(lista_de_adyacencia("Lista de adyacencia", &grafo));
    salida.explicar(format!(
        "Hay {} y {} de {}. Un genoma que explique todas las lecturas es un camino que pasa exactamente una vez por cada nodo (un camino hamiltoniano), y no se conoce ningún método rápido para encontrarlo en grafos grandes: por eso los ensambladores usan el grafo de De Bruijn.",
        cantidad(nodos, "lectura", "lecturas"),
        cantidad(aristas, "solapamiento", "solapamientos"),
        cantidad(k - 1, "letra", "letras")
    ));
    if varias_salidas == 0
        && sin_entrada.len() == 1
        && sin_salida == 1
        && es_una_fila(&grafo, sin_entrada[0])
    {
        salida.explicar(
            "Cada lectura tiene a lo sumo una que la sigue: el grafo es una sola fila y el orden de las lecturas es único.",
        );
    } else if varias_salidas > 0 {
        salida.explicar(format!(
            "{} más de una posible sucesora: sus últimas {} aparecen al comienzo de varias lecturas, señal de una repetición. Ahí el camino puede seguir por más de un lado.",
            cantidad(varias_salidas, "lectura tiene", "lecturas tienen"),
            cantidad(k - 1, "letra", "letras")
        ));
    }
    match sin_entrada.len() {
        0 => salida.explicar(
            "Todas las lecturas tienen alguna que las precede: el genoma podría ser circular.",
        ),
        1 => salida.explicar(format!(
            "Solo «{}» no tiene ninguna lectura que la preceda: es la candidata a ser el comienzo del genoma.",
            etiqueta_corta(&grafo.nodos[sin_entrada[0]])
        )),
        otras => salida.explicar(format!(
            "{} lecturas no tienen ninguna que las preceda. En un genoma lineal bien cubierto habría una sola, la primera; así que falta cobertura o hay errores, y el ensamblaje quedará en varios trozos.",
            miles(otras)
        )),
    }
    Ok(salida)
}

fn es_una_fila(grafo: &Grafo, inicio: usize) -> bool {
    let mut visitados = 1;
    let mut actual = inicio;
    while let Some(&siguiente) = grafo.salidas[actual].first() {
        if visitados > grafo.nodos.len() {
            return false;
        }
        visitados += 1;
        actual = siguiente;
    }
    visitados == grafo.nodos.len()
}

struct ResumenDeBruijn {
    entradas: Vec<usize>,
    ramificados: Vec<usize>,
    aristas: usize,
}

fn mostrar_de_bruijn(salida: &mut Salida, grafo: &Grafo) -> ResumenDeBruijn {
    let entradas = ensamblaje::grados_de_entrada(grafo);
    let ramificados = ramificados(grafo, &entradas);
    let aristas = ensamblaje::cantidad_de_aristas(grafo);
    salida.agregar(Bloque::valor(
        "Nodos ((k − 1)-meros distintos)",
        miles(grafo.nodos.len()),
    ));
    salida.agregar(Bloque::valor("Aristas (k-meros)", miles(aristas)));
    if se_puede_dibujar(grafo) {
        salida.agregar(dibujo(
            "Grafo de De Bruijn",
            grafo,
            aristas_de_bruijn(grafo),
        ));
    }
    salida.agregar(lista_de_adyacencia("Lista de adyacencia", grafo));
    if !ramificados.is_empty() {
        salida.agregar(tabla_de_nodos(
            "Nodos que se ramifican",
            grafo,
            &entradas,
            &ramificados,
        ));
    }
    ResumenDeBruijn {
        entradas,
        ramificados,
        aristas,
    }
}

fn explicar_ramificados(
    salida: &mut Salida,
    resumen: &ResumenDeBruijn,
    origen: &str,
    hay_recorrido: bool,
) {
    if resumen.ramificados.is_empty() && hay_recorrido {
        salida.explicar(
            "Ningún nodo se ramifica: el grafo no tiene bifurcaciones y hay una única forma de recorrer todas sus aristas.",
        );
    } else if resumen.ramificados.is_empty() {
        salida.explicar(
            "Ningún nodo se ramifica, pero las aristas no forman un solo recorrido: quedan trozos separados que no se pueden unir sin más k-meros.",
        );
    } else {
        salida.explicar(format!(
            "{} más de una entrada o de una salida: son (k − 1)-meros que aparecen varias veces {origen}. Prueba con un k más grande: cuando k − 1 supera el largo de una repetición, esa ramificación desaparece.",
            cantidad(resumen.ramificados.len(), "nodo tiene", "nodos tienen")
        ));
    }
}

fn de_bruijn_texto(entrada: &Entrada) -> Result<Salida, String> {
    let texto = entrada.secuencia("texto")?;
    let k = leer_k(entrada, "k", 2, texto.len())?;
    let grafo = ensamblaje::de_bruijn_de_texto(&texto, k)?;
    let mut salida = Salida::default();
    let resumen = mostrar_de_bruijn(&mut salida, &grafo);
    salida.explicar(format!(
        "El texto tiene {} y cada uno es una arista que va de su prefijo a su sufijo; los nodos ({}) son los (k − 1)-meros distintos. Recorrer las aristas en el orden del texto lo vuelve a deletrear.",
        cantidad(resumen.aristas, "k-mero", "k-meros"),
        miles(grafo.nodos.len())
    ));
    if grafo.nodos.len() < resumen.aristas + 1 {
        salida.explicar(format!(
            "El camino original tenía {} nodos, uno por posición; al pegar los (k − 1)-meros iguales quedaron {}.",
            miles(resumen.aristas + 1),
            miles(grafo.nodos.len())
        ));
    }
    explicar_ramificados(&mut salida, &resumen, "en el texto", true);
    if texto[..k - 1] == texto[texto.len() - (k - 1)..] {
        salida.explicar(
            "El texto empieza y termina con el mismo (k − 1)-mero, así que el recorrido vuelve a su punto de partida: es un ciclo.",
        );
    }
    Ok(salida)
}

fn de_bruijn_kmeros(entrada: &Entrada) -> Result<Salida, String> {
    let kmeros = leer_kmeros(entrada, "kmeros")?;
    let lista = rebanadas(&kmeros);
    let grafo = ensamblaje::de_bruijn_de_kmeros(&lista)?;
    let mut salida = Salida::default();
    let resumen = mostrar_de_bruijn(&mut salida, &grafo);
    let balance = ensamblaje::balance_de_nodos(&grafo);
    let desbalanceados: Vec<usize> = (0..balance.len())
        .filter(|&nodo| balance[nodo] != 0)
        .collect();
    if !desbalanceados.is_empty() {
        salida.agregar(tabla_de_nodos(
            "Nodos desbalanceados",
            &grafo,
            &resumen.entradas,
            &desbalanceados,
        ));
    }
    let recorrido = ensamblaje::camino_euleriano(&grafo);
    let hay_recorrido = recorrido.is_ok();
    match recorrido {
        Ok(camino) if desbalanceados.is_empty() => salida.explicar(format!(
            "Todos los nodos están balanceados (entran tantas aristas como salen), así que hay un ciclo euleriano: los k-meros podrían venir de un genoma circular, y la cadena que deletrea mide {} letras.",
            miles(camino.len() - 1 + grafo.nodos[0].len())
        )),
        Ok(camino) => salida.explicar(format!(
            "De «{}» sale una arista más de las que entran y a «{}» llega una más de las que salen: hay un camino euleriano que empieza en el primero y termina en el segundo, y deletrea una cadena de {} letras con esta composición.",
            etiqueta_corta(&grafo.nodos[camino[0]]),
            etiqueta_corta(&grafo.nodos[camino[camino.len() - 1]]),
            miles(camino.len() - 1 + grafo.nodos[0].len())
        )),
        Err(motivo) => {
            salida.explicar(motivo);
            salida.explicar(
                "Con lecturas reales esto pasa cuando faltan k-meros (zonas sin cobertura) o hay errores de secuenciación: el genoma no se puede deletrear de una sola vez.",
            );
        }
    }
    explicar_ramificados(&mut salida, &resumen, "en el genoma", hay_recorrido);
    let mut vistas = HashSet::new();
    let copias = lista.iter().filter(|kmero| !vistas.insert(**kmero)).count();
    if copias > 0 {
        salida.explicar(format!(
            "La lista trae {} repetidos: cada copia es otra arista entre los mismos dos nodos.",
            cantidad(copias, "k-mero", "k-meros")
        ));
    }
    Ok(salida)
}

fn explicar_desbalance(salida: &mut Salida, grafo: &Grafo) {
    let entradas = ensamblaje::grados_de_entrada(grafo);
    let balance = ensamblaje::balance_de_nodos(grafo);
    let desbalanceados: Vec<usize> = (0..balance.len())
        .filter(|&nodo| balance[nodo] != 0)
        .collect();
    if !desbalanceados.is_empty() {
        salida.agregar(tabla_de_nodos(
            "Nodos desbalanceados",
            grafo,
            &entradas,
            &desbalanceados,
        ));
    }
    if se_puede_dibujar(grafo) {
        salida.agregar(dibujo("Grafo", grafo, aristas_sin_nombre(grafo)));
    }
}

fn explicar_elecciones(salida: &mut Salida, grafo: &Grafo) {
    let con_eleccion = grafo
        .salidas
        .iter()
        .filter(|destinos| destinos.len() > 1)
        .count();
    if con_eleccion == 0 {
        salida.explicar(
            "De cada nodo sale una sola arista, así que no hubo nada que elegir: el recorrido es el único posible.",
        );
    } else {
        salida.explicar(format!(
            "De {} sale más de una arista. El algoritmo de Hierholzer toma la primera sin usar y, si se atasca antes de terminar, retrocede por la pila y empalma un ciclo nuevo; con otro orden de aristas saldría otro recorrido igual de válido.",
            cantidad(con_eleccion, "nodo", "nodos")
        ));
    }
}

fn ciclo(entrada: &Entrada) -> Result<Salida, String> {
    let grafo = leer_grafo(entrada, "grafo")?;
    let mut salida = Salida::default();
    match ensamblaje::ciclo_euleriano(&grafo) {
        Ok(recorrido) => {
            salida.agregar(Bloque::texto(
                "Ciclo euleriano",
                camino_como_texto(&grafo, &recorrido),
            ));
            salida.agregar(Bloque::valor(
                "Aristas recorridas",
                miles(recorrido.len() - 1),
            ));
            if se_puede_dibujar(&grafo) {
                salida.agregar(dibujo(
                    "Orden en que se recorre cada arista",
                    &grafo,
                    aristas_numeradas(&recorrido),
                ));
            }
            salida.explicar(format!(
                "El ciclo recorre {}, cada una una sola vez, y vuelve a «{}», donde empezó. Si empezaras en otro nodo del ciclo obtendrías el mismo recorrido rotado.",
                cantidad(recorrido.len() - 1, "arista", "aristas"),
                etiqueta_corta(&grafo.nodos[recorrido[0]])
            ));
            explicar_elecciones(&mut salida, &grafo);
        }
        Err(motivo) => {
            salida.agregar(Bloque::valor("¿Tiene ciclo euleriano?", "No"));
            explicar_desbalance(&mut salida, &grafo);
            salida.explicar(motivo);
            salida.explicar(
                "Para que haya un ciclo euleriano, en cada nodo deben entrar tantas aristas como salen y todas las aristas deben estar conectadas entre sí.",
            );
        }
    }
    Ok(salida)
}

fn camino(entrada: &Entrada) -> Result<Salida, String> {
    let grafo = leer_grafo(entrada, "grafo")?;
    let mut salida = Salida::default();
    match ensamblaje::camino_euleriano(&grafo) {
        Ok(recorrido) => {
            let inicio = recorrido[0];
            let fin = recorrido[recorrido.len() - 1];
            salida.agregar(Bloque::texto(
                "Camino euleriano",
                camino_como_texto(&grafo, &recorrido),
            ));
            salida.agregar(Bloque::valor(
                "Aristas recorridas",
                miles(recorrido.len() - 1),
            ));
            if se_puede_dibujar(&grafo) {
                salida.agregar(dibujo(
                    "Orden en que se recorre cada arista",
                    &grafo,
                    aristas_numeradas(&recorrido),
                ));
            }
            if inicio == fin {
                salida.explicar(format!(
                    "Todos los nodos están balanceados, así que el camino es en realidad un ciclo: empieza y termina en «{}».",
                    etiqueta_corta(&grafo.nodos[inicio])
                ));
            } else {
                salida.explicar(format!(
                    "De «{}» sale una arista más de las que entran, así que el camino tiene que empezar ahí; a «{}» llega una más de las que salen, y ahí termina. Todos los demás nodos están balanceados.",
                    etiqueta_corta(&grafo.nodos[inicio]),
                    etiqueta_corta(&grafo.nodos[fin])
                ));
            }
            explicar_elecciones(&mut salida, &grafo);
        }
        Err(motivo) => {
            salida.agregar(Bloque::valor("¿Tiene camino euleriano?", "No"));
            explicar_desbalance(&mut salida, &grafo);
            salida.explicar(motivo);
            salida.explicar(
                "Para que haya un camino euleriano, todos los nodos deben estar balanceados salvo, a lo sumo, uno con una salida de más (el inicio) y otro con una entrada de más (el final), y las aristas deben estar conectadas.",
            );
        }
    }
    Ok(salida)
}

fn deletrear_recorrido(grafo: &Grafo, recorrido: &[usize]) -> Result<Vec<u8>, String> {
    let etiquetas: Vec<&[u8]> = recorrido
        .iter()
        .map(|&nodo| grafo.nodos[nodo].as_slice())
        .collect();
    ensamblaje::deletrear_camino(&etiquetas)
}

fn otra_reconstruccion(grafo: &Grafo, original: &[u8], con_eleccion: &[usize]) -> Option<Vec<u8>> {
    let mut variantes: Vec<&[usize]> = con_eleccion
        .iter()
        .take(INTENTOS_ALTERNATIVOS)
        .map(std::slice::from_ref)
        .collect();
    variantes.push(con_eleccion);
    for nodos in variantes {
        let mut copia = grafo.clone();
        for &nodo in nodos {
            copia.salidas[nodo].reverse();
        }
        let Ok(recorrido) = ensamblaje::camino_euleriano(&copia) else {
            continue;
        };
        if let Ok(otra) = deletrear_recorrido(&copia, &recorrido)
            && otra != original
        {
            return Some(otra);
        }
    }
    None
}

fn reconstruccion(entrada: &Entrada) -> Result<Salida, String> {
    let kmeros = leer_kmeros(entrada, "kmeros")?;
    let lista = rebanadas(&kmeros);
    let cadena = ensamblaje::reconstruir_desde_kmeros(&lista)?;
    let grafo = ensamblaje::de_bruijn_de_kmeros(&lista)?;
    let k = lista[0].len();
    let entradas = ensamblaje::grados_de_entrada(&grafo);
    let nodos_ramificados = ramificados(&grafo, &entradas);
    let mut salida = Salida::default();
    salida.agregar(Bloque::texto("Cadena reconstruida", recortar(&cadena)));
    salida.agregar(Bloque::valor(
        "Largo",
        format!("{} letras", miles(cadena.len())),
    ));
    if se_puede_dibujar(&grafo)
        && let Ok(recorrido) = ensamblaje::camino_euleriano(&grafo)
    {
        salida.agregar(dibujo(
            "Grafo de De Bruijn con el orden del camino euleriano",
            &grafo,
            aristas_numeradas(&recorrido),
        ));
    }
    if !nodos_ramificados.is_empty() {
        salida.agregar(tabla_de_nodos(
            "Nodos que se ramifican",
            &grafo,
            &entradas,
            &nodos_ramificados,
        ));
    }
    let con_eleccion: Vec<usize> = (0..grafo.salidas.len())
        .filter(|&nodo| grafo.salidas[nodo].len() > 1)
        .collect();
    let muy_grande = ensamblaje::cantidad_de_aristas(&grafo) > MAXIMO_ARISTAS_ALTERNATIVAS;
    let alternativa = if con_eleccion.is_empty() || muy_grande {
        None
    } else {
        otra_reconstruccion(&grafo, &cadena, &con_eleccion)
    };
    if let Some(otra) = &alternativa {
        salida.agregar(Bloque::texto(
            "Otra cadena con los mismos k-meros",
            recortar(otra),
        ));
    }
    let aportes = if lista.len() == 1 {
        "el único k-mero agrega una más".to_string()
    } else {
        format!(
            "cada uno de los {} k-meros agrega una más",
            miles(lista.len())
        )
    };
    salida.explicar(format!(
        "La cadena mide {} letras: el primer nodo aporta k − 1 = {} y {aportes}.",
        miles(cadena.len()),
        k - 1
    ));
    if nodos_ramificados.is_empty() {
        salida.explicar(
            "Ningún nodo se ramifica, así que hay un solo camino euleriano: esta es la única cadena posible con estos k-meros.",
        );
    } else if alternativa.is_some() {
        salida.explicar(format!(
            "{} más de una entrada o salida: son repeticiones de al menos k − 1 = {} letras. Recorriéndolos en otro orden sale otra cadena distinta con exactamente los mismos k-meros (la ves arriba), así que los k-meros solos no alcanzan para saber cuál es la verdadera; lecturas más largas o pares de lecturas sí ayudan.",
            cantidad(nodos_ramificados.len(), "nodo tiene", "nodos tienen"),
            k - 1
        ));
    } else if con_eleccion.is_empty() {
        salida.explicar(format!(
            "{} más de una entrada, pero de ningún nodo sale más de una arista: el recorrido nunca tiene que elegir, así que esta es la única cadena posible.",
            cantidad(nodos_ramificados.len(), "nodo tiene", "nodos tienen")
        ));
    } else if muy_grande {
        salida.explicar(format!(
            "{} más de una entrada o salida, así que podría haber otras cadenas con los mismos k-meros; con más de {} aristas no se buscan aquí.",
            cantidad(nodos_ramificados.len(), "nodo tiene", "nodos tienen"),
            miles(MAXIMO_ARISTAS_ALTERNATIVAS)
        ));
    } else {
        salida.explicar(format!(
            "{} más de una entrada o salida, pero al invertir el orden de las salidas en cada nodo que se ramifica salió siempre la misma cadena: las repeticiones no se cruzan entre sí y probablemente esta es la única respuesta.",
            cantidad(nodos_ramificados.len(), "nodo tiene", "nodos tienen")
        ));
    }
    if ensamblaje::balance_de_nodos(&grafo)
        .iter()
        .all(|&diferencia| diferencia == 0)
    {
        salida.explicar(
            "Todos los nodos están balanceados: los k-meros podrían venir de un genoma circular, y empezar en otro nodo daría la misma cadena rotada.",
        );
    }
    explicar_ubicacion(
        &mut salida,
        "Esta cadena",
        &cadena,
        "No aparece tal cual en E. coli. Si tus k-meros venían de E. coli, el recorrido tomó otro camino en alguna repetición.",
    );
    Ok(salida)
}

fn universal(entrada: &Entrada) -> Result<Salida, String> {
    let k = entrada.natural("k")?;
    if k == 0 || k > ensamblaje::MAXIMO_K_UNIVERSAL {
        return Err(format!(
            "k debe estar entre 1 y {}: con k más grande la cadena tendría más de un millón de letras.",
            ensamblaje::MAXIMO_K_UNIVERSAL
        ));
    }
    let cadena = ensamblaje::cadena_universal_circular(k)?;
    let mut vuelta = cadena.clone();
    vuelta.extend_from_slice(&cadena[..k - 1]);
    let en_circulo = ensamblaje::composicion(&vuelta, k);
    let distintos: HashSet<&[u8]> = en_circulo.iter().copied().collect();
    let mut salida = Salida::default();
    salida.agregar(Bloque::texto("Cadena circular", recortar(&cadena)));
    salida.agregar(Bloque::valor(
        "Largo",
        format!("{} letras", miles(cadena.len())),
    ));
    if (2..=6).contains(&k) {
        let grafo = ensamblaje::de_bruijn_de_kmeros(&en_circulo)?;
        let mut numeradas = Vec::new();
        let mut usadas = vec![0; grafo.nodos.len()];
        let mut posicion_de: HashMap<&[u8], usize> = HashMap::new();
        for (indice, nodo) in grafo.nodos.iter().enumerate() {
            posicion_de.insert(nodo.as_slice(), indice);
        }
        for (orden, kmero) in en_circulo.iter().enumerate() {
            let desde = posicion_de[&kmero[..k - 1]];
            let hasta = posicion_de[&kmero[1..]];
            usadas[desde] += 1;
            numeradas.push((desde, hasta, (orden + 1).to_string()));
        }
        salida.agregar(dibujo(
            "Grafo de De Bruijn binario, con el orden del ciclo",
            &grafo,
            numeradas,
        ));
    }
    let filas: Vec<Vec<String>> = en_circulo
        .iter()
        .take(MAXIMO_FILAS)
        .enumerate()
        .map(|(posicion, kmero)| vec![posicion.to_string(), letras(kmero)])
        .collect();
    salida.agregar(Bloque::Tabla {
        etiqueta: "Cada k-mero y dónde empieza (contando desde 0, leyendo en círculo)".to_string(),
        columnas: vec!["Posición".to_string(), "k-mero".to_string()],
        filas,
        total: en_circulo.len(),
    });
    if distintos.len() == cadena.len() {
        salida.explicar(format!(
            "Leyendo en círculo, la cadena de {} letras contiene cada uno de los 2^{k} = {} k-meros binarios exactamente una vez. No puede ser más corta: cada posición es el comienzo de un solo k-mero.",
            miles(cadena.len()),
            miles(cadena.len())
        ));
    }
    salida.explicar(
        "Se construye como un ciclo euleriano del grafo de De Bruijn de todos los k-meros binarios. Ese grafo siempre está balanceado (a cada nodo llegan dos aristas y salen dos), por eso el ciclo existe para cualquier k.",
    );
    let exponente = (1usize << (k - 1)) - k;
    let cuantas = if exponente <= 30 {
        miles(1usize << exponente)
    } else {
        format!("2^{}", miles(exponente))
    };
    if exponente == 0 {
        salida.explicar(format!(
            "Con k = {k} hay 2^(2^(k − 1) − k) = 1 sola cadena así, sin contar sus rotaciones: esta."
        ));
    } else {
        salida.explicar(format!(
            "Hay 2^(2^(k − 1) − k) = {cuantas} cadenas así distintas sin contar sus rotaciones; esta es solo una de ellas, la que sale de recorrer las aristas en orden."
        ));
    }
    Ok(salida)
}

fn camino_con_huecos(entrada: &Entrada) -> Result<Salida, String> {
    let pares = leer_pares(entrada, "pares")?;
    let d = leer_d(entrada)?;
    let genoma = ensamblaje::deletrear_camino_con_huecos(&pares, d)?;
    let k = pares[0].primera.len();
    let salto = k + d;
    let mut salida = Salida::default();
    salida.agregar(Bloque::texto("Cadena", recortar(&genoma)));
    salida.agregar(Bloque::valor(
        "Largo",
        format!("{} letras", miles(genoma.len())),
    ));
    if pares.len() <= MAXIMO_FILAS_ESCALONADAS && genoma.len() <= MAXIMO_LARGO_ESCALONADO {
        let mut filas = vec![fila_escalonada(
            "cadena".to_string(),
            0,
            &genoma,
            genoma.len(),
        )];
        for (posicion, par) in pares.iter().enumerate() {
            let mut contenido = par.primera.clone();
            contenido.extend(std::iter::repeat_n(b'-', d));
            contenido.extend_from_slice(&par.segunda);
            filas.push(fila_escalonada(
                format!("par {}", posicion + 1),
                posicion,
                &contenido,
                genoma.len(),
            ));
        }
        salida.agregar(Bloque::alineamiento("Cada par en su lugar", filas));
    }
    let prefijo = k + pares.len() - 1;
    let superpuestas = prefijo - salto;
    let encaje = if superpuestas == 0 {
        "Como hay justo d + 1 pares, los dos tramos no se superponen: el segundo empieza donde termina el primero".to_string()
    } else {
        format!(
            "Coinciden en {} donde se superponen",
            cantidad(superpuestas, "letra", "letras")
        )
    };
    salida.explicar(format!(
        "Las primeras lecturas deletrean un tramo de {} letras y las segundas otro del mismo largo que empieza k + d = {} letras más adelante. {encaje}, y juntos dan 2k + d + n − 1 = {} letras.",
        miles(prefijo),
        miles(salto),
        miles(genoma.len())
    ));
    explicar_ubicacion(
        &mut salida,
        "Esta cadena",
        &genoma,
        "No aparece tal cual en el genoma de E. coli, en ninguna de sus dos hebras.",
    );
    Ok(salida)
}

fn pares(entrada: &Entrada) -> Result<Salida, String> {
    let lista = leer_pares(entrada, "pares")?;
    let d = leer_d(entrada)?;
    let semilla = entrada.natural("semilla")? as u64;
    let armado = ensamblaje::reconstruir_desde_pares(&lista, d, semilla)?;
    let grafo = ensamblaje::de_bruijn_de_pares(&lista)?;
    let k = lista[0].primera.len();
    let entradas = ensamblaje::grados_de_entrada(&grafo);
    let ramificados_pares = ramificados(&grafo, &entradas).len();
    let mut vistas = HashSet::new();
    let mut sueltas: Vec<&[u8]> = Vec::new();
    for par in &lista {
        for lectura in [par.primera.as_slice(), par.segunda.as_slice()] {
            if vistas.insert(lectura) {
                sueltas.push(lectura);
            }
        }
    }
    let grafo_suelto = ensamblaje::de_bruijn_de_kmeros(&sueltas)?;
    let entradas_sueltas = ensamblaje::grados_de_entrada(&grafo_suelto);
    let ramificados_sueltos = ramificados(&grafo_suelto, &entradas_sueltas).len();
    let mut salida = Salida::default();
    salida.agregar(Bloque::texto(
        "Cadena reconstruida",
        recortar(&armado.genoma),
    ));
    salida.agregar(Bloque::valor(
        "Largo",
        format!("{} letras", miles(armado.genoma.len())),
    ));
    salida.agregar(Bloque::valor("Intentos", miles(armado.intentos)));
    salida.agregar(Bloque::valor(
        "Nodos que se ramifican (lecturas sueltas)",
        miles(ramificados_sueltos),
    ));
    salida.agregar(Bloque::valor(
        "Nodos que se ramifican (pares)",
        miles(ramificados_pares),
    ));
    if se_puede_dibujar(&grafo) {
        salida.agregar(dibujo(
            "Grafo de De Bruijn pareado",
            &grafo,
            aristas_sin_nombre(&grafo),
        ));
    }
    salida.explicar(format!(
        "Cada par cubre un tramo de 2k + d = {} letras con un hueco de {d} en medio. En el grafo pareado cada nodo es un par de (k − 1)-meros, así que dos lugares del genoma solo se confunden si se repiten las dos mitades a la vez.",
        miles(2 * k + d)
    ));
    if ramificados_pares < ramificados_sueltos {
        salida.explicar(format!(
            "Con las lecturas sueltas, {} se ramificarían; con los pares quedan {}. Los pares resuelven las repeticiones más cortas que el tramo de 2k + d letras que cubre cada par.",
            cantidad(ramificados_sueltos, "nodo", "nodos"),
            miles(ramificados_pares)
        ));
    } else if ramificados_sueltos > 0 {
        salida.explicar(format!(
            "Los pares no redujeron las ramificaciones ({} con lecturas sueltas y {} con pares): las repeticiones son más largas que el tramo que cubre cada par.",
            miles(ramificados_sueltos),
            miles(ramificados_pares)
        ));
    }
    if armado.intentos == 1 {
        salida.explicar(
            "El primer camino euleriano ya deletreó una cadena coherente: las primeras y las segundas lecturas encajaron a la distancia d.",
        );
    } else {
        salida.explicar(format!(
            "Hicieron falta {} intentos: no todo camino euleriano del grafo pareado deletrea un genoma coherente, así que se barajó el orden de las aristas con la semilla {semilla} y se volvió a probar.",
            miles(armado.intentos)
        ));
    }
    explicar_ubicacion(
        &mut salida,
        "Esta cadena",
        &armado.genoma,
        "No aparece tal cual en el genoma de E. coli, en ninguna de sus dos hebras.",
    );
    Ok(salida)
}

fn no_ramificados(entrada: &Entrada) -> Result<Salida, String> {
    let grafo = leer_grafo(entrada, "grafo")?;
    let caminos = ensamblaje::caminos_no_ramificados(&grafo);
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Caminos", miles(caminos.len())));
    let total = caminos.len();
    let lineas: Vec<String> = caminos
        .iter()
        .take(MAXIMO_LINEAS)
        .map(|recorrido| camino_como_texto(&grafo, recorrido))
        .collect();
    salida.agregar(Bloque::Lista {
        etiqueta: "Caminos máximos no ramificados".to_string(),
        elementos: lineas,
        total,
    });
    if se_puede_dibujar(&grafo) {
        let mut numeradas = Vec::new();
        for (numero, recorrido) in caminos.iter().enumerate() {
            for par in recorrido.windows(2) {
                numeradas.push((par[0], par[1], format!("c{}", numero + 1)));
            }
        }
        salida.agregar(dibujo(
            "Cada arista con el número de su camino",
            &grafo,
            numeradas,
        ));
    }
    if caminos.is_empty() {
        salida.explicar("El grafo no tiene aristas, así que no hay ningún camino.");
        return Ok(salida);
    }
    if caminos.len() == 1 {
        salida.explicar(
            "Hay un solo camino máximo no ramificado: avanza mientras los nodos tengan exactamente una entrada y una salida, y se corta en cuanto llega a uno que se ramifica.",
        );
    } else {
        salida.explicar(format!(
            "Hay {} caminos máximos no ramificados: cada uno avanza mientras los nodos tengan exactamente una entrada y una salida, y se corta en cuanto llega a uno que se ramifica.",
            miles(caminos.len())
        ));
    }
    let entradas = ensamblaje::grados_de_entrada(&grafo);
    let ciclos = caminos
        .iter()
        .filter(|recorrido| {
            let primero = recorrido[0];
            recorrido.last() == Some(&primero)
                && entradas[primero] == 1
                && grafo.salidas[primero].len() == 1
        })
        .count();
    if ciclos == 1 {
        salida.explicar(
            "Uno de ellos es un ciclo aislado: todos sus nodos tienen una entrada y una salida, y no se conecta con el resto del grafo.",
        );
    } else if ciclos > 1 {
        salida.explicar(format!(
            "{} de ellos son ciclos aislados: todos sus nodos tienen una entrada y una salida, y no se conectan con el resto del grafo.",
            miles(ciclos)
        ));
    }
    let mas_largo = caminos
        .iter()
        .map(|recorrido| recorrido.len() - 1)
        .max()
        .unwrap_or(0);
    if mas_largo == 1 {
        salida.explicar(
            "Todos los caminos tienen una sola arista: en ningún nodo se encadenan dos aristas sin que haya una bifurcación, un inicio o un final.",
        );
    } else {
        salida.explicar(format!(
            "El más largo tiene {} aristas. Si el grafo fuera de De Bruijn, cada camino así deletrearía un contig: un trozo del genoma que se arma sin dudas.",
            miles(mas_largo)
        ));
    }
    Ok(salida)
}

fn contigs(entrada: &Entrada) -> Result<Salida, String> {
    let kmeros = leer_kmeros(entrada, "kmeros")?;
    let modo = entrada.opcion("copias", COPIAS)?;
    let mut lista = rebanadas(&kmeros);
    let leidos = lista.len();
    if modo == COPIAS[0] {
        let mut vistos = HashSet::new();
        lista.retain(|kmero| vistos.insert(*kmero));
    }
    let repetidos = leidos - lista.len();
    let piezas = ensamblaje::generar_contigs(&lista)?;
    let grafo = ensamblaje::de_bruijn_de_kmeros(&lista)?;
    let entradas = ensamblaje::grados_de_entrada(&grafo);
    let nodos_ramificados = ramificados(&grafo, &entradas).len();
    let k = lista[0].len();
    let mut orden: Vec<usize> = (0..piezas.len()).collect();
    orden.sort_by(|&uno, &otro| piezas[otro].len().cmp(&piezas[uno].len()));
    let longitudes: Vec<usize> = orden.iter().map(|&indice| piezas[indice].len()).collect();
    let medida = ensamblaje::medir_n50(&longitudes);
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Contigs", miles(piezas.len())));
    salida.agregar(Bloque::valor(
        "El más largo",
        format!("{} letras", miles(longitudes[0])),
    ));
    salida.agregar(Bloque::valor(
        "N50",
        format!("{} letras", miles(medida.n50)),
    ));
    let mostradas = longitudes.len().min(MAXIMO_BARRAS);
    salida.agregar(Bloque::barras(
        if longitudes.len() > MAXIMO_BARRAS {
            "Largo de los contigs más largos"
        } else {
            "Largo de cada contig"
        },
        "Contig (del más largo al más corto)",
        "Letras",
        (1..=mostradas).map(|numero| numero.to_string()).collect(),
        longitudes[..mostradas]
            .iter()
            .map(|&largo| largo as f64)
            .collect(),
    ));
    let filas: Vec<Vec<String>> = orden
        .iter()
        .take(MAXIMO_FILAS)
        .enumerate()
        .map(|(posicion, &indice)| {
            vec![
                (posicion + 1).to_string(),
                miles(piezas[indice].len()),
                etiqueta_de_lista(&piezas[indice]),
            ]
        })
        .collect();
    salida.agregar(Bloque::Tabla {
        etiqueta: "Contigs del más largo al más corto".to_string(),
        columnas: vec![
            "Contig".to_string(),
            "Largo".to_string(),
            "Secuencia".to_string(),
        ],
        filas,
        total: piezas.len(),
    });
    salida.explicar(format!(
        "Salen {}; el más largo mide {} letras y el N50 es {}: la mitad de todas las letras está en contigs de ese largo o más.",
        cantidad(piezas.len(), "contig", "contigs"),
        miles(longitudes[0]),
        miles(medida.n50)
    ));
    if piezas.len() == 1 {
        salida.explicar(
            "Un solo contig: ningún nodo se ramifica y los k-meros alcanzan para armar todo de una pieza.",
        );
    } else if nodos_ramificados == 0 {
        salida.explicar(
            "Ningún nodo se ramifica, pero los contigs no se tocan: faltan los k-meros que los unirían (zonas sin cobertura).",
        );
    } else {
        let donde = if nodos_ramificados == 1 {
            "en el único nodo que se ramifica".to_string()
        } else {
            format!(
                "en cada uno de los {} nodos que se ramifican",
                miles(nodos_ramificados)
            )
        };
        salida.explicar(format!(
            "El ensamblaje se corta {donde}: ahí hay una repetición de al menos k − 1 = {} letras, o falta cobertura, y no se sabe con certeza cómo seguir.",
            k - 1
        ));
        let cortos = longitudes.iter().filter(|&&largo| largo == k).count();
        if cortos > 0 {
            salida.explicar(format!(
                "{} solo k letras: son una sola arista entre dos puntos de corte, típicas de zonas repetidas.",
                cantidad(cortos, "contig mide", "contigs miden")
            ));
        }
        salida.explicar(
            "Con lecturas más largas (un k mayor), las repeticiones más cortas que k − 1 dejan de ramificar el grafo y los contigs se alargan.",
        );
    }
    if repetidos > 0 {
        salida.explicar(format!(
            "La lista traía {} repetidos y se usó una sola copia de cada uno, como hacen los ensambladores. Con «una arista por copia», cada copia de una repetición formaría su propio contig.",
            cantidad(repetidos, "k-mero", "k-meros")
        ));
    }
    explicar_ubicacion(&mut salida, "El contig más largo", &piezas[orden[0]], "");
    Ok(salida)
}

fn leer_longitudes(entrada: &Entrada, campo: &str) -> Result<Vec<usize>, String> {
    let texto = entrada.texto(campo)?;
    let piezas: Vec<&str> = texto
        .split(|letra: char| letra.is_whitespace() || letra == ',' || letra == ';')
        .filter(|pieza| !pieza.is_empty())
        .collect();
    let numeros: Vec<Option<usize>> = piezas.iter().map(|pieza| pieza.parse().ok()).collect();
    let longitudes: Vec<usize> = if !numeros.is_empty() && numeros.iter().all(Option::is_some) {
        numeros.into_iter().flatten().collect()
    } else {
        let secuencias = entrada.secuencias(campo)?;
        for (numero, contig) in secuencias.iter().enumerate() {
            if let Some(letra) = contig.iter().find(|letra| !letra.is_ascii_alphabetic()) {
                return Err(format!(
                    "El contig {} tiene «{}»: escribe solo secuencias de letras o solo números (sus largos).",
                    numero + 1,
                    *letra as char
                ));
            }
        }
        secuencias.iter().map(|contig| contig.len()).collect()
    };
    if longitudes.len() > MAXIMO_CONTIGS {
        return Err(format!(
            "Hay {} contigs; aquí se admiten hasta {}.",
            miles(longitudes.len()),
            miles(MAXIMO_CONTIGS)
        ));
    }
    let mut total: usize = 0;
    for &largo in &longitudes {
        total = total
            .checked_add(largo)
            .filter(|&suma| suma <= usize::MAX / 2)
            .ok_or_else(|| "La suma de los largos es demasiado grande.".to_string())?;
    }
    if total == 0 {
        return Err("Todos los contigs miden 0: no hay nada que medir.".to_string());
    }
    Ok(longitudes)
}

fn n50(entrada: &Entrada) -> Result<Salida, String> {
    let longitudes = leer_longitudes(entrada, "contigs")?;
    let medida = ensamblaje::medir_n50(&longitudes);
    let mut ordenadas = longitudes.clone();
    ordenadas.sort_unstable_by(|uno, otro| otro.cmp(uno));
    let mut acumulado = 0.0;
    let mut puntos = Vec::with_capacity(ordenadas.len() + 1);
    puntos.push([0.0, 0.0]);
    for (posicion, &largo) in ordenadas.iter().enumerate() {
        acumulado += largo as f64;
        puntos.push([(posicion + 1) as f64, acumulado]);
    }
    let promedio = medida.total as f64 / longitudes.len() as f64;
    let cortos = longitudes
        .iter()
        .filter(|&&largo| largo < LARGO_MINIMO_GENBANK)
        .count();
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        "N50",
        format!("{} letras", miles(medida.n50)),
    ));
    salida.agregar(Bloque::valor(
        "L50 (contigs que hacen falta)",
        miles(medida.l50),
    ));
    salida.agregar(Bloque::valor("Contigs", miles(longitudes.len())));
    salida.agregar(Bloque::valor(
        "Total",
        format!("{} letras", miles(medida.total)),
    ));
    salida.agregar(Bloque::valor(
        "El más largo",
        format!("{} letras", miles(ordenadas[0])),
    ));
    salida.agregar(Bloque::linea(
        "Letras acumuladas, sumando del contig más largo al más corto",
        "Contigs sumados",
        "Letras acumuladas",
        &puntos,
        vec![Marca {
            x: medida.l50 as f64,
            texto: format!("N50 = {}", miles(medida.n50)),
        }],
    ));
    let mostradas = ordenadas.len().min(MAXIMO_BARRAS);
    salida.agregar(Bloque::barras(
        if ordenadas.len() > MAXIMO_BARRAS {
            "Largo de los contigs más largos"
        } else {
            "Largo de cada contig"
        },
        "Contig (del más largo al más corto)",
        "Letras",
        (1..=mostradas).map(|numero| numero.to_string()).collect(),
        ordenadas[..mostradas]
            .iter()
            .map(|&largo| largo as f64)
            .collect(),
    ));
    if medida.l50 == 1 {
        salida.explicar(format!(
            "Ordenando de mayor a menor, el contig más largo ya tiene la mitad de las {} letras, así que su largo, {}, es el N50.",
            miles(medida.total),
            miles(medida.n50)
        ));
    } else {
        salida.explicar(format!(
            "Ordenando de mayor a menor, los {} contigs más largos ya suman la mitad de las {} letras; el más corto de ellos mide {}: ese es el N50.",
            miles(medida.l50),
            miles(medida.total),
            miles(medida.n50)
        ));
    }
    if longitudes.len() > 1 {
        let hacen_falta = if medida.l50 == 1 {
            "basta 1".to_string()
        } else {
            format!("hacen falta {}", miles(medida.l50))
        };
        salida.explicar(format!(
            "Para llegar a la mitad {hacen_falta} de los {} contigs ({}), y el largo promedio es {} letras. Cuantos menos contigs hagan falta, más continuo es el ensamblaje.",
            miles(longitudes.len()),
            porcentaje(medida.l50 as f64 / longitudes.len() as f64, 0),
            miles(promedio.round() as usize)
        ));
    }
    if cortos == 1 {
        salida.explicar(format!(
            "1 contig ({}) mide menos de {LARGO_MINIMO_GENBANK} letras, el mínimo que GenBank acepta para depositar un contig.",
            porcentaje(1.0 / longitudes.len() as f64, 0)
        ));
    } else if cortos > 1 {
        salida.explicar(format!(
            "{} contigs ({}) miden menos de {LARGO_MINIMO_GENBANK} letras, el mínimo que GenBank acepta para depositar un contig.",
            miles(cortos),
            porcentaje(cortos as f64 / longitudes.len() as f64, 0)
        ));
    }
    salida.explicar(
        "El N50 mide continuidad, no exactitud: un error que une dos trozos que no van juntos también lo hace subir.",
    );
    Ok(salida)
}

struct Traza {
    k: usize,
    fase: &'static str,
    nodos: Vec<String>,
    desde: Vec<usize>,
    hasta: Vec<usize>,
    kmeros: Vec<String>,
    usada: Vec<bool>,
    nodos_visibles: usize,
    aristas_visibles: usize,
    balance: Option<Vec<i64>>,
    inicio: Option<usize>,
    pila: Vec<usize>,
    camino: Vec<usize>,
    camino_al_reves: bool,
    arista: Option<usize>,
    nodo: Option<usize>,
    texto: Option<String>,
    pasos: Vec<Paso>,
}

impl Traza {
    fn anotar(&mut self, linea: &str, nota: String) {
        let aristas: Vec<Value> = (0..self.aristas_visibles)
            .map(|arista| {
                json!({
                    "desde": self.desde[arista],
                    "hasta": self.hasta[arista],
                    "kmero": self.kmeros[arista],
                    "usada": self.usada[arista],
                })
            })
            .collect();
        let estado = json!({
            "fase": self.fase,
            "k": self.k,
            "nodos": &self.nodos[..self.nodos_visibles],
            "aristas": aristas,
            "balance": self.balance,
            "inicio": self.inicio,
            "pila": self.pila,
            "camino": self.camino,
            "camino_al_reves": self.camino_al_reves,
            "arista": self.arista,
            "nodo": self.nodo,
            "texto": self.texto,
        });
        self.pasos.push(Paso::nuevo(linea, nota, estado));
    }
}

fn pasos_reconstruccion(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let kmeros = leer_kmeros(entrada, "kmeros")?;
    if kmeros.len() > MAXIMO_KMEROS_ANIMADOS {
        return Err(format!(
            "El paso a paso admite hasta {MAXIMO_KMEROS_ANIMADOS} k-meros, para que el grafo se pueda dibujar; aquí hay {}.",
            miles(kmeros.len())
        ));
    }
    if kmeros[0].len() > MAXIMO_K_ANIMADO {
        return Err(format!(
            "El paso a paso admite k-meros de hasta {MAXIMO_K_ANIMADO} letras."
        ));
    }
    let lista = rebanadas(&kmeros);
    let esperado = ensamblaje::reconstruir_desde_kmeros(&lista)?;
    let grafo = ensamblaje::de_bruijn_de_kmeros(&lista)?;
    let k = lista[0].len();
    let mut traza = Traza {
        k,
        fase: "grafo",
        nodos: grafo.nodos.iter().map(|nodo| letras(nodo)).collect(),
        desde: Vec::new(),
        hasta: Vec::new(),
        kmeros: lista.iter().map(|kmero| letras(kmero)).collect(),
        usada: vec![false; lista.len()],
        nodos_visibles: 0,
        aristas_visibles: 0,
        balance: None,
        inicio: None,
        pila: Vec::new(),
        camino: Vec::new(),
        camino_al_reves: true,
        arista: None,
        nodo: None,
        texto: None,
        pasos: Vec::new(),
    };
    let aristas_de_nodo = construir_en_pasos(&mut traza, &grafo, &lista);
    let inicio = elegir_inicio_en_pasos(&mut traza, &grafo);
    recorrer_en_pasos(&mut traza, &grafo, &aristas_de_nodo, inicio);
    let texto = deletrear_en_pasos(&mut traza);
    if texto.as_bytes() != esperado.as_slice() {
        return Err("El paso a paso no llegó a la misma cadena que el algoritmo.".to_string());
    }
    Ok(muestrear(traza.pasos, MAXIMO_PASOS))
}

fn construir_en_pasos(traza: &mut Traza, grafo: &Grafo, lista: &[&[u8]]) -> Vec<Vec<usize>> {
    let k = traza.k;
    let mut indice_de: HashMap<&[u8], usize> = HashMap::new();
    for (indice, nodo) in grafo.nodos.iter().enumerate() {
        indice_de.insert(nodo.as_slice(), indice);
    }
    let mut aristas_de_nodo = vec![Vec::new(); grafo.nodos.len()];
    for (numero, kmero) in lista.iter().enumerate() {
        let desde = indice_de[&kmero[..k - 1]];
        let hasta = indice_de[&kmero[1..]];
        let ya_habia = traza.nodos_visibles;
        let repetida = aristas_de_nodo[desde]
            .iter()
            .any(|&otra| traza.hasta[otra] == hasta);
        aristas_de_nodo[desde].push(numero);
        traza.desde.push(desde);
        traza.hasta.push(hasta);
        traza.nodos_visibles = ya_habia.max(desde + 1).max(hasta + 1);
        traza.aristas_visibles = numero + 1;
        traza.arista = Some(numero);
        traza.nodo = Some(desde);
        let prefijo = &traza.nodos[desde];
        let sufijo = &traza.nodos[hasta];
        let mut nota = format!(
            "El k-mero «{}» es una arista de su prefijo «{prefijo}» a su sufijo «{sufijo}».",
            traza.kmeros[numero]
        );
        let nuevo_desde = desde >= ya_habia;
        let nuevo_hasta = hasta >= ya_habia && hasta != desde;
        if desde == hasta {
            nota.push_str(&format!(
                " El prefijo y el sufijo son iguales: la arista es un bucle que sale de «{prefijo}» y vuelve a él."
            ));
            if repetida {
                nota.push_str(" Ya había un bucle igual: este k-mero está repetido.");
            }
        } else if nuevo_desde && nuevo_hasta {
            nota.push_str(" Los dos nodos son nuevos.");
        } else if nuevo_desde || nuevo_hasta {
            nota.push_str(
                " Uno de los nodos es nuevo; el otro ya estaba y la arista se pega a él.",
            );
        } else if repetida {
            nota.push_str(
                " Los dos nodos ya estaban y ya había una arista igual: este k-mero está repetido.",
            );
        } else {
            nota.push_str(" Los dos nodos ya estaban: la arista se pega a ellos.");
        }
        traza.anotar("grafo.salidas[desde].push(hasta);", nota);
    }
    aristas_de_nodo
}

fn elegir_inicio_en_pasos(traza: &mut Traza, grafo: &Grafo) -> usize {
    let balance = ensamblaje::balance_de_nodos(grafo);
    let con_salida_de_mas = balance.iter().position(|&diferencia| diferencia == 1);
    traza.fase = "inicio";
    traza.balance = Some(balance);
    traza.arista = None;
    traza.nodo = None;
    traza.anotar(
        "let balance = balance_de_nodos(grafo);",
        "Para cada nodo se resta cuántas aristas entran de cuántas salen.".to_string(),
    );
    match con_salida_de_mas {
        Some(nodo) => {
            traza.inicio = Some(nodo);
            traza.nodo = Some(nodo);
            let nota = format!(
                "De «{}» sale una arista más de las que entran: el camino euleriano tiene que empezar ahí.",
                traza.nodos[nodo]
            );
            traza.anotar("inicio = Some(nodo);", nota);
            nodo
        }
        None => {
            let nodo = grafo
                .salidas
                .iter()
                .position(|destinos| !destinos.is_empty())
                .unwrap_or(0);
            traza.inicio = Some(nodo);
            traza.nodo = Some(nodo);
            let nota = format!(
                "Todos los nodos están balanceados: el recorrido es un ciclo y empieza en «{}», el primer nodo con aristas.",
                traza.nodos[nodo]
            );
            traza.anotar(
                "let inicio = inicio.unwrap_or_else(|| primer_nodo_con_salidas(grafo));",
                nota,
            );
            nodo
        }
    }
}

fn recorrer_en_pasos(
    traza: &mut Traza,
    grafo: &Grafo,
    aristas_de_nodo: &[Vec<usize>],
    inicio: usize,
) {
    traza.fase = "recorrido";
    traza.pila = vec![inicio];
    let nota = format!(
        "La pila empieza con «{}». El camino todavía está vacío.",
        traza.nodos[inicio]
    );
    traza.anotar("let mut pila = vec![inicio];", nota);
    let mut aristas_usadas = vec![0; grafo.salidas.len()];
    while let Some(&actual) = traza.pila.last() {
        let usadas = aristas_usadas[actual];
        if usadas < grafo.salidas[actual].len() {
            let destino = grafo.salidas[actual][usadas];
            let arista = aristas_de_nodo[actual][usadas];
            aristas_usadas[actual] += 1;
            traza.usada[arista] = true;
            traza.pila.push(destino);
            traza.arista = Some(arista);
            traza.nodo = Some(destino);
            let nota = format!(
                "Desde «{}» se toma la arista «{}», todavía sin usar, y se avanza a «{}», que entra a la pila.",
                traza.nodos[actual], traza.kmeros[arista], traza.nodos[destino]
            );
            traza.anotar("pila.push(destino);", nota);
        } else {
            traza.pila.pop();
            traza.camino.push(actual);
            traza.arista = None;
            traza.nodo = Some(actual);
            let mut nota = format!(
                "«{}» ya no tiene aristas sin usar: sale de la pila y se anota en el camino, que se arma de atrás hacia adelante.",
                traza.nodos[actual]
            );
            match traza.pila.last() {
                None => nota.push_str(" La pila quedó vacía: se usaron todas las aristas."),
                Some(&anterior) if aristas_usadas[anterior] < grafo.salidas[anterior].len() => {
                    nota.push_str(&format!(
                        " «{}» todavía tiene aristas sin usar: desde ahí se recorrerá un ciclo que se empalma en este punto.",
                        traza.nodos[anterior]
                    ));
                }
                Some(_) => {}
            }
            traza.anotar("camino.push(actual);", nota);
        }
    }
    traza.camino.reverse();
    traza.camino_al_reves = false;
    traza.nodo = None;
    traza.anotar(
        "camino.reverse();",
        "Se da vuelta el camino: ahora va del nodo inicial al final y usa cada arista exactamente una vez."
            .to_string(),
    );
}

fn deletrear_en_pasos(traza: &mut Traza) -> String {
    traza.fase = "deletreo";
    let camino = traza.camino.clone();
    let mut texto = traza.nodos[camino[0]].clone();
    traza.texto = Some(texto.clone());
    traza.nodo = Some(camino[0]);
    let nota = format!(
        "La cadena empieza con el primer nodo completo, «{}».",
        traza.nodos[camino[0]]
    );
    traza.anotar("texto.extend_from_slice(kmeros[0]);", nota);
    for posicion in 1..camino.len() {
        let nodo = camino[posicion];
        let ultima = traza.nodos[nodo].chars().last().unwrap_or(' ');
        texto.push(ultima);
        traza.texto = Some(texto.clone());
        traza.nodo = Some(nodo);
        let mut nota = format!(
            "Del nodo «{}» se agrega solo su última letra, «{ultima}».",
            traza.nodos[nodo]
        );
        if posicion + 1 == camino.len() {
            nota.push_str(&format!(
                " Listo: la cadena mide {} letras.",
                miles(texto.chars().count())
            ));
        }
        traza.anotar("texto.push(kmero[k - 1]);", nota);
    }
    texto
}

fn muestrear(pasos: Vec<Paso>, maximo: usize) -> Vec<Paso> {
    if pasos.len() <= maximo || maximo < 2 {
        return pasos;
    }
    let ultimo = pasos.len() - 1;
    let elegidos: HashSet<usize> = (0..maximo)
        .map(|numero| numero * ultimo / (maximo - 1))
        .collect();
    pasos
        .into_iter()
        .enumerate()
        .filter(|(indice, _)| elegidos.contains(indice))
        .map(|(_, paso)| paso)
        .collect()
}

const DETALLE_INICIAL: usize = 4;
const FILAS_DE_PIZARRA: usize = 8;
const DISTINTO_DEL_ALGORITMO: &str = "El paso a paso no llegó al mismo resultado que el algoritmo.";

fn se_muestra(indice: usize, total: usize, cupo: usize) -> bool {
    indice < DETALLE_INICIAL
        || indice + 1 >= total
        || indice.is_multiple_of((total / cupo.max(1)).max(1))
}

fn con_variables(valores: &[(&str, String)]) -> Pizarra {
    let mut pizarra = Pizarra::nueva();
    for (nombre, valor) in valores {
        pizarra = pizarra.variable(nombre, valor);
    }
    pizarra
}

fn nombre_de_nodo(grafo: &Grafo, nodo: usize) -> String {
    etiqueta_corta(&grafo.nodos[nodo])
}

fn nombres_de_nodos(grafo: &Grafo, nodos: &[usize]) -> Vec<String> {
    nodos
        .iter()
        .map(|&nodo| nombre_de_nodo(grafo, nodo))
        .collect()
}

fn nodos_con_puntos(grafo: &Grafo, nodos: &[usize], total: usize) -> String {
    let mut texto = nombres_de_nodos(grafo, nodos).join(", ");
    if total > nodos.len() {
        texto = format!("…, {texto}");
    }
    texto
}

fn cola_de_nodos(nodos: &[usize], cuantos: usize) -> &[usize] {
    &nodos[nodos.len().saturating_sub(cuantos)..]
}

fn camino_corto(grafo: &Grafo, camino: &[usize]) -> String {
    let mostrados = camino.len().min(6);
    let mut texto = nombres_de_nodos(grafo, &camino[..mostrados]).join(" -> ");
    if camino.len() > mostrados {
        texto.push_str(&format!(" -> … ({} nodos en total)", miles(camino.len())));
    }
    texto
}

fn marcas_de_solapamiento(
    numero: usize,
    k: usize,
    con_nueva: bool,
) -> Vec<(usize, usize, Resalte)> {
    let mut marcas = vec![(numero - 1, numero, Resalte::Actual)];
    if k > 1 {
        marcas.push((numero, numero + k - 1, Resalte::Comparado));
    }
    if con_nueva {
        marcas.push((numero + k - 1, numero + k, Resalte::Hallado));
    }
    marcas
}

fn elegir_nodos(importantes: &[usize], visibles: usize) -> Vec<usize> {
    let mut elegidos: Vec<usize> = importantes
        .iter()
        .copied()
        .filter(|&nodo| nodo < visibles)
        .collect();
    let tope = importantes
        .iter()
        .max()
        .map_or(FILAS_DE_PIZARRA, |&mayor| mayor + 1)
        .min(visibles);
    let mut siguiente = tope;
    while elegidos.len() < FILAS_DE_PIZARRA && siguiente > 0 {
        siguiente -= 1;
        if !elegidos.contains(&siguiente) {
            elegidos.push(siguiente);
        }
    }
    elegidos.sort_unstable();
    elegidos.dedup();
    elegidos
}

fn lista_de_destinos(grafo: &Grafo, nodo: usize) -> String {
    let destinos = &grafo.salidas[nodo];
    let mut texto = nombres_de_nodos(grafo, &destinos[..destinos.len().min(4)]).join(", ");
    if destinos.len() > 4 {
        texto.push_str(", …");
    }
    texto
}

fn tabla_de_grafo(
    pizarra: Pizarra,
    titulo: &str,
    grafo: &Grafo,
    visibles: usize,
    importantes: &[usize],
) -> Pizarra {
    let elegidos = elegir_nodos(importantes, visibles);
    let filas: Vec<Vec<String>> = elegidos
        .iter()
        .map(|&nodo| {
            vec![
                nodo.to_string(),
                nombre_de_nodo(grafo, nodo),
                lista_de_destinos(grafo, nodo),
            ]
        })
        .collect();
    let resaltar: Vec<(usize, usize)> = importantes
        .iter()
        .filter_map(|nodo| elegidos.iter().position(|elegido| elegido == nodo))
        .map(|fila| (fila, 1))
        .collect();
    pizarra.tabla(titulo, &["n.º", "nodo", "flechas hacia"], filas, &resaltar)
}

fn tabla_de_grados(
    pizarra: Pizarra,
    titulo: &str,
    grafo: &Grafo,
    entradas: &[usize],
    importantes: &[usize],
) -> Pizarra {
    let elegidos = elegir_nodos(importantes, grafo.nodos.len());
    let filas: Vec<Vec<String>> = elegidos
        .iter()
        .map(|&nodo| {
            let salen = grafo.salidas[nodo].len();
            vec![
                nombre_de_nodo(grafo, nodo),
                entradas[nodo].to_string(),
                salen.to_string(),
                (salen as i64 - entradas[nodo] as i64).to_string(),
            ]
        })
        .collect();
    let resaltar: Vec<(usize, usize)> = importantes
        .iter()
        .filter_map(|nodo| elegidos.iter().position(|elegido| elegido == nodo))
        .map(|fila| (fila, 3))
        .collect();
    pizarra.tabla(
        titulo,
        &["nodo", "entran", "salen", "salen − entran"],
        filas,
        &resaltar,
    )
}

fn pasos_composicion(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let texto = entrada.secuencia("texto")?;
    let k = leer_k(entrada, "k", 1, texto.len())?;
    let orden = entrada.opcion("orden", ORDENES)?;
    let kmeros = ensamblaje::composicion(&texto, k);
    let total = kmeros.len();
    if total != texto.len() - k + 1 {
        return Err(DISTINTO_DEL_ALGORITMO.to_string());
    }
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        con_variables(&[("k", k.to_string()), ("texto.len()", miles(texto.len()))])
            .secuencia(&texto, &[])
            .paso(
                "if k == 0 || k > texto.len() {",
                format!(
                    "Partimos de un texto de {} letras y de k = {k}. La meta es sacar todos los k-meros: cada trozo de {k} letras seguidas. Primero se revisa que k no sea 0 ni mayor que el texto; aquí está bien.",
                    miles(texto.len())
                ),
            ),
    );
    for posicion in 0..total {
        if !se_muestra(posicion, total, 60) {
            continue;
        }
        let desde = (posicion + 1).saturating_sub(FILAS_DE_PIZARRA);
        let ultimos: Vec<String> = kmeros[desde..=posicion]
            .iter()
            .map(|kmero| letras(kmero))
            .collect();
        let resaltado = ultimos.len() - 1;
        let ventana = letras(kmeros[posicion]);
        let nota = if posicion == 0 {
            format!(
                "La ventana empieza en la primera letra y abarca {k}: es «{ventana}». Se guarda como primer k-mero."
            )
        } else if posicion < DETALLE_INICIAL {
            format!(
                "La ventana se corre una letra a la derecha: sale «{}» y entra «{}». El nuevo k-mero es «{ventana}».",
                texto[posicion - 1] as char,
                texto[posicion + k - 1] as char
            )
        } else {
            format!(
                "Ventana número {}, en las letras {} a {}: «{ventana}». Las ventanas del medio hacen lo mismo y no se muestran todas.",
                miles(posicion + 1),
                miles(posicion + 1),
                miles(posicion + k)
            )
        };
        narracion.agregar(
            con_variables(&[
                ("k", k.to_string()),
                ("texto.len()", miles(texto.len())),
                ("posicion", miles(posicion)),
                ("ventana", ventana),
            ])
            .secuencia(&texto, &[(posicion, posicion + k, Resalte::Actual)])
            .lista(
                &format!("k-meros sacados hasta aquí (van {})", miles(posicion + 1)),
                ultimos,
                &[resaltado],
            )
            .paso("texto.windows(k).collect()", nota),
        );
    }
    let orden_final = if orden == ORDENES[0] {
        "Para mostrártelos se ordenan alfabéticamente; así llegarían de una máquina secuenciadora, sin saber dónde estaba cada uno."
    } else {
        "Se muestran en el orden en que salieron: cada uno empieza una letra después del anterior."
    };
    Ok(narracion.terminar(
        con_variables(&[
            ("k", k.to_string()),
            ("texto.len()", miles(texto.len())),
            ("resultado.len()", miles(total)),
        ])
        .secuencia(&texto, &[])
        .paso(
            "texto.windows(k).collect()",
            format!(
                "Listo: {} k-meros, que son {} − {k} + 1: uno por cada lugar donde cabe la ventana. {orden_final} En «Qué significa» se cuenta cuáles se repiten.",
                miles(total),
                miles(texto.len())
            ),
        ),
    ))
}

fn narrar_deletreo(
    narracion: &mut Narracion,
    kmeros: &[&[u8]],
    cadena: &[u8],
    cupo: usize,
    unidad: (&str, &str),
) -> Result<(), String> {
    let total = kmeros.len();
    let k = kmeros[0].len();
    if cadena.len() != k + total - 1 || &cadena[..k] != kmeros[0] {
        return Err(DISTINTO_DEL_ALGORITMO.to_string());
    }
    let primeros: Vec<String> = kmeros
        .iter()
        .take(FILAS_DE_PIZARRA)
        .map(|kmero| letras(kmero))
        .collect();
    narracion.agregar(
        con_variables(&[("kmeros.len()", miles(total)), ("k", k.to_string())])
            .lista(&format!("{} en el orden dado (los primeros)", unidad.1), primeros, &[0])
            .paso(
                "let k = largo_comun(kmeros)?;",
                format!(
                    "Hay {} y todos miden k = {k} letras. Cada uno debería empezar una letra después del anterior, así que de uno al siguiente solo cambia una letra. La meta: pegarlos en una sola cadena.",
                    cantidad(total, unidad.0, unidad.1)
                ),
            ),
    );
    narracion.agregar(
        con_variables(&[
            ("kmeros.len()", miles(total)),
            ("k", k.to_string()),
            ("texto", etiqueta_corta(kmeros[0])),
        ])
        .secuencia(&cadena[..k], &[(0, k, Resalte::Hallado)])
        .paso(
            "texto.extend_from_slice(kmeros[0]);",
            format!(
                "El texto empieza con el primero completo, «{}»: {k} letras.",
                etiqueta_corta(kmeros[0])
            ),
        ),
    );
    for numero in 1..total {
        if !se_muestra(numero - 1, total - 1, cupo) {
            continue;
        }
        let anterior = kmeros[numero - 1];
        let kmero = kmeros[numero];
        if anterior[1..] != kmero[..k - 1] || cadena[numero + k - 1] != kmero[k - 1] {
            return Err(DISTINTO_DEL_ALGORITMO.to_string());
        }
        let antes = &cadena[..numero + k - 1];
        let despues = &cadena[..numero + k];
        let letra = kmero[k - 1] as char;
        if numero < DETALLE_INICIAL {
            let nota = if k == 1 {
                format!(
                    "Con k = 1 no hay nada que comparar: «{}» se pega tal cual.",
                    letras(kmero)
                )
            } else {
                format!(
                    "Antes de pegar «{}» se comprueba que encaje: las últimas {} letras de «{}» son «{}» y las primeras {} de «{}» son «{}». Son iguales, así que sí continúa al anterior.",
                    letras(kmero),
                    k - 1,
                    letras(anterior),
                    letras(&anterior[1..]),
                    k - 1,
                    letras(kmero),
                    letras(&kmero[..k - 1])
                )
            };
            narracion.agregar(
                con_variables(&[
                    ("kmeros.len()", miles(total)),
                    ("k", k.to_string()),
                    ("texto", etiqueta_corta(antes)),
                    ("numero", numero.to_string()),
                    ("anterior", etiqueta_corta(anterior)),
                    ("kmero", etiqueta_corta(kmero)),
                ])
                .secuencia(antes, &marcas_de_solapamiento(numero, k, false))
                .paso("if anterior[1..] != kmero[..k - 1] {", nota),
            );
        }
        let nota = if k == 1 {
            format!(
                "Con k = 1 no hay solapamiento: la letra «{letra}» se pega a continuación. El texto mide {} letras.",
                miles(despues.len())
            )
        } else {
            format!(
                "«{}» ya comparte sus primeras {} letras con el final del texto, así que solo aporta la última, «{letra}». El texto mide ahora {} letras.",
                letras(kmero),
                k - 1,
                miles(despues.len())
            )
        };
        let nota = if numero >= DETALLE_INICIAL {
            format!("Vuelta número {}: {nota}", miles(numero + 1))
        } else {
            nota
        };
        narracion.agregar(
            con_variables(&[
                ("kmeros.len()", miles(total)),
                ("k", k.to_string()),
                ("texto", etiqueta_corta(despues)),
                ("numero", numero.to_string()),
                ("anterior", etiqueta_corta(anterior)),
                ("kmero", etiqueta_corta(kmero)),
            ])
            .secuencia(despues, &marcas_de_solapamiento(numero, k, true))
            .paso("texto.push(kmero[k - 1]);", nota),
        );
    }
    Ok(())
}

fn pasos_camino_de_kmeros(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let kmeros = leer_kmeros(entrada, "kmeros")?;
    let lista = rebanadas(&kmeros);
    let cadena = ensamblaje::deletrear_camino(&lista)?;
    let k = lista[0].len();
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narrar_deletreo(&mut narracion, &lista, &cadena, 80, ("k-mero", "k-meros"))?;
    Ok(narracion.terminar(
        con_variables(&[
            ("kmeros.len()", miles(lista.len())),
            ("k", k.to_string()),
            ("texto.len()", miles(cadena.len())),
        ])
        .secuencia(&cadena, &[(0, cadena.len(), Resalte::Hallado)])
        .paso(
            "Ok(texto)",
            format!(
                "Listo: {} de largo {k} deletrean una cadena de {} letras, que son {} + {k} − 1. En «Qué significa» se explica por qué cada k-mero aporta una sola letra.",
                cantidad(lista.len(), "k-mero", "k-meros"),
                miles(cadena.len()),
                miles(lista.len())
            ),
        ),
    ))
}

fn pasos_solapamiento(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let kmeros = leer_kmeros(entrada, "kmeros")?;
    if kmeros.len() > MAXIMO_SOLAPAMIENTO {
        return Err(format!(
            "El grafo de solapamiento admite hasta {} k-meros; aquí hay {}.",
            miles(MAXIMO_SOLAPAMIENTO),
            miles(kmeros.len())
        ));
    }
    let lista = rebanadas(&kmeros);
    let grafo = ensamblaje::grafo_de_solapamiento(&lista)?;
    let total = lista.len();
    let k = lista[0].len();
    let aristas_reales = ensamblaje::cantidad_de_aristas(&grafo);
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let primeros: Vec<String> = lista
        .iter()
        .take(FILAS_DE_PIZARRA)
        .map(|kmero| letras(kmero))
        .collect();
    narracion.agregar(
        con_variables(&[("kmeros.len()", miles(total)), ("k", k.to_string())])
            .lista("k-meros (los primeros)", primeros, &[])
            .paso(
                "let k = largo_comun(kmeros)?;",
                format!(
                    "Hay {} y todos miden k = {k} letras. Dos k-meros se solapan cuando las últimas {} letras de uno son las primeras {} del otro. La meta: dibujar una flecha de cada k-mero a los que pueden ir justo después.",
                    cantidad(total, "k-mero", "k-meros"),
                    k - 1,
                    k - 1
                ),
            ),
    );
    let mut por_prefijo: HashMap<&[u8], Vec<usize>> = HashMap::with_capacity(total);
    let mut orden_de_prefijos: Vec<&[u8]> = Vec::new();
    for (indice, kmero) in lista.iter().enumerate() {
        let prefijo = &kmero[..k - 1];
        let ya_estaba = por_prefijo.contains_key(prefijo);
        por_prefijo.entry(prefijo).or_default().push(indice);
        if !ya_estaba {
            orden_de_prefijos.push(prefijo);
        }
        if !se_muestra(indice, total, 25) {
            continue;
        }
        let mut filas_elegidas: Vec<&[u8]> = cola_de_prefijos(&orden_de_prefijos, 5);
        if !filas_elegidas.contains(&prefijo) {
            filas_elegidas.push(prefijo);
        }
        let filas: Vec<Vec<String>> = filas_elegidas
            .iter()
            .map(|clave| {
                let indices = &por_prefijo[clave];
                let mut texto: Vec<String> = indices
                    .iter()
                    .take(6)
                    .map(|&numero| format!("k-mero {}", numero + 1))
                    .collect();
                if indices.len() > 6 {
                    texto.push("…".to_string());
                }
                vec![letras(clave), texto.join(", ")]
            })
            .collect();
        let resaltada = filas_elegidas
            .iter()
            .position(|clave| *clave == prefijo)
            .unwrap_or(0);
        let nota = if ya_estaba {
            format!(
                "«{}» empieza con «{}», sus primeras {} letras. Esa libreta ya tenía ese prefijo: ahora {} k-meros empiezan igual.",
                letras(kmero),
                letras(prefijo),
                k - 1,
                por_prefijo[prefijo].len()
            )
        } else {
            format!(
                "«{}» empieza con «{}», sus primeras {} letras. Es un prefijo nuevo: se abre una página de la libreta y se anota este k-mero.",
                letras(kmero),
                letras(prefijo),
                k - 1
            )
        };
        narracion.agregar(
            con_variables(&[
                ("kmeros.len()", miles(total)),
                ("k", k.to_string()),
                ("indice", indice.to_string()),
                ("kmero", etiqueta_corta(kmero)),
            ])
            .secuencia(kmero, &[(0, k - 1, Resalte::Actual)])
            .tabla(
                "por_prefijo: qué k-meros empiezan con cada prefijo",
                &["prefijo", "k-meros"],
                filas,
                &[(resaltada, 1)],
            )
            .paso(
                "por_prefijo.entry(&kmero[..k - 1]).or_default().push(indice);",
                nota,
            ),
        );
    }
    let mut aristas = 0;
    for (indice, kmero) in lista.iter().enumerate() {
        let sufijo = &kmero[1..];
        let siguientes = por_prefijo.get(sufijo);
        let cuantas = siguientes.map_or(0, Vec::len);
        aristas += cuantas;
        if !se_muestra(indice, total, 12) {
            continue;
        }
        let mut nombres: Vec<String> = siguientes
            .map(|indices| {
                indices
                    .iter()
                    .take(FILAS_DE_PIZARRA)
                    .map(|&otro| letras(lista[otro]))
                    .collect()
            })
            .unwrap_or_default();
        let resaltados: Vec<usize> = (0..nombres.len()).collect();
        if nombres.is_empty() {
            nombres.push("(ninguno)".to_string());
        }
        let nota = if cuantas == 0 {
            format!(
                "Las últimas {} letras de «{}» son «{}». Ningún k-mero empieza así: de este no sale ninguna flecha. Llevamos {} solapamientos.",
                k - 1,
                letras(kmero),
                letras(sufijo),
                miles(aristas)
            )
        } else {
            format!(
                "Las últimas {} letras de «{}» son «{}». Se busca en la libreta: {} empiezan así y pueden ir justo después. Llevamos {} solapamientos.",
                k - 1,
                letras(kmero),
                letras(sufijo),
                cantidad(cuantas, "k-mero", "k-meros"),
                miles(aristas)
            )
        };
        narracion.agregar(
            con_variables(&[
                ("kmeros.len()", miles(total)),
                ("k", k.to_string()),
                ("indice", indice.to_string()),
                ("kmero", etiqueta_corta(kmero)),
                ("aristas", miles(aristas)),
            ])
            .secuencia(kmero, &[(1, k, Resalte::Actual)])
            .lista("k-meros que pueden seguirle", nombres, &resaltados)
            .paso("aristas += siguientes.len();", nota),
        );
    }
    narracion.agregar(
        con_variables(&[
            ("kmeros.len()", miles(total)),
            ("aristas", miles(aristas)),
        ])
        .paso(
            "if aristas > MAXIMO_ARISTAS_DE_SOLAPAMIENTO {",
            format!(
                "Con {} solapamientos contados el grafo cabe: no se pasa del máximo de {}. Un k-mero puede contarse con sí mismo si se solapa consigo (como AAAA); esa flecha se descarta al armar el grafo.",
                miles(aristas),
                miles(ensamblaje::MAXIMO_ARISTAS_DE_SOLAPAMIENTO)
            ),
        ),
    );
    let mut dibujadas = 0;
    for (indice, kmero) in lista.iter().enumerate() {
        let destinos = &grafo.salidas[indice];
        dibujadas += destinos.len();
        if !se_muestra(indice, total, 25) {
            continue;
        }
        let desde = (indice + 1).saturating_sub(6);
        let lineas: Vec<String> = (desde..=indice)
            .map(|nodo| {
                let hacia: Vec<String> = grafo.salidas[nodo]
                    .iter()
                    .take(4)
                    .map(|&otro| letras(lista[otro]))
                    .collect();
                let mut hacia = hacia.join(", ");
                if grafo.salidas[nodo].len() > 4 {
                    hacia.push_str(", …");
                }
                if hacia.is_empty() {
                    hacia = "(ninguno)".to_string();
                }
                format!("{} -> {hacia}", letras(lista[nodo]))
            })
            .collect();
        let resaltado = lineas.len() - 1;
        let nota = if destinos.is_empty() {
            format!(
                "«{}» no tiene sucesores: es un posible final del genoma, o falta un k-mero.",
                letras(kmero)
            )
        } else {
            format!(
                "Se dibuja una flecha de «{}» a cada uno de sus {}: {}. Llevamos {} flechas.",
                letras(kmero),
                cantidad(destinos.len(), "sucesor", "sucesores"),
                nombres_cortos(destinos, &lista),
                miles(dibujadas)
            )
        };
        narracion.agregar(
            con_variables(&[
                ("kmeros.len()", miles(total)),
                ("indice", indice.to_string()),
                ("kmero", etiqueta_corta(kmero)),
                ("siguiente", nombres_cortos(destinos, &lista)),
            ])
            .lista("salidas de los últimos k-meros", lineas, &[resaltado])
            .paso("grafo.salidas[indice].push(siguiente);", nota),
        );
    }
    if dibujadas != aristas_reales {
        return Err(DISTINTO_DEL_ALGORITMO.to_string());
    }
    Ok(narracion.terminar(
        con_variables(&[
            ("grafo.nodos.len()", miles(total)),
            ("aristas", miles(aristas_reales)),
        ])
        .paso(
            "Ok(grafo)",
            format!(
                "Listo: el grafo tiene {} (uno por k-mero) y {}. En «Qué significa» se dice cuál lectura puede ser el comienzo y dónde hay ramas.",
                cantidad(total, "nodo", "nodos"),
                cantidad(aristas_reales, "flecha", "flechas")
            ),
        ),
    ))
}

fn cola_de_prefijos<'a>(orden: &[&'a [u8]], cuantos: usize) -> Vec<&'a [u8]> {
    orden[orden.len().saturating_sub(cuantos)..].to_vec()
}

fn nombres_cortos(destinos: &[usize], lista: &[&[u8]]) -> String {
    let mut nombres: Vec<String> = destinos
        .iter()
        .take(4)
        .map(|&otro| letras(lista[otro]))
        .collect();
    if destinos.len() > 4 {
        nombres.push("…".to_string());
    }
    if nombres.is_empty() {
        return "(ninguno)".to_string();
    }
    nombres.join(", ")
}

fn narrar_nodos(
    narracion: &mut Narracion,
    kmeros: &[&[u8]],
    texto: Option<&[u8]>,
    cupo: usize,
) -> Result<Grafo, String> {
    let k = kmeros[0].len();
    let total = kmeros.len();
    let mut grafo = Grafo::default();
    let mut indices: HashMap<&[u8], usize> = HashMap::with_capacity(total);
    for (numero, kmero) in kmeros.iter().enumerate() {
        let prefijo = &kmero[..k - 1];
        let sufijo = &kmero[1..];
        let nodos_antes = grafo.nodos.len();
        let desde = nodo_para_en_pasos(prefijo, &mut indices, &mut grafo);
        let hasta = nodo_para_en_pasos(sufijo, &mut indices, &mut grafo);
        let repetida = se_muestra(numero, total, cupo) && grafo.salidas[desde].contains(&hasta);
        grafo.salidas[desde].push(hasta);
        if !se_muestra(numero, total, cupo) {
            continue;
        }
        let ventana = match texto {
            Some(texto) => (texto, numero, numero + k),
            None => (*kmero, 0, k),
        };
        let marcas = [(ventana.1, ventana.2, Resalte::Actual)];
        let nodos_ahora = grafo.nodos.len();
        let nuevo_desde = desde >= nodos_antes;
        let nuevo_hasta = hasta >= nodos_antes && hasta != desde;
        let estado = |nodo: usize, es_nuevo: bool| {
            if es_nuevo {
                format!("«{}» es nuevo", nombre_de_nodo(&grafo, nodo))
            } else {
                format!("«{}» ya estaba", nombre_de_nodo(&grafo, nodo))
            }
        };
        let base = || {
            vec![
                ("k", k.to_string()),
                ("kmero", etiqueta_corta(kmero)),
                ("desde", desde.to_string()),
                ("hasta", hasta.to_string()),
            ]
        };
        if numero < 2 {
            narracion.agregar(
                tabla_de_grafo(
                    con_variables(&base()[..2]).secuencia(ventana.0, &marcas),
                    "nodos del grafo hasta ahora",
                    &grafo,
                    grafo.nodos.len(),
                    &[desde],
                )
                .paso(
                    "let desde = nodo_para(&kmero[..k - 1], &mut indices, &mut grafo);",
                    format!(
                        "El k-mero «{}» sale del nodo que lleva su prefijo, «{}» (sus primeras {} letras). {}; {}.",
                        etiqueta_corta(kmero),
                        etiqueta_corta(prefijo),
                        k - 1,
                        estado(desde, nuevo_desde),
                        if nuevo_desde {
                            "se crea con el número siguiente"
                        } else {
                            "se reutiliza el mismo nodo"
                        }
                    ),
                ),
            );
            narracion.agregar(
                tabla_de_grafo(
                    con_variables(&base()[..3]).secuencia(ventana.0, &marcas),
                    "nodos del grafo hasta ahora",
                    &grafo,
                    grafo.nodos.len(),
                    &[desde, hasta],
                )
                .paso(
                    "let hasta = nodo_para(&kmero[1..], &mut indices, &mut grafo);",
                    format!(
                        "Y llega al nodo de su sufijo, «{}» (sus últimas {} letras). {}.",
                        etiqueta_corta(sufijo),
                        k - 1,
                        estado(hasta, nuevo_hasta)
                    ),
                ),
            );
        }
        let nota = if desde == hasta {
            format!(
                "El prefijo y el sufijo son iguales: la flecha del k-mero «{}» es un bucle de «{}» a sí mismo.",
                etiqueta_corta(kmero),
                nombre_de_nodo(&grafo, desde)
            )
        } else if repetida {
            format!(
                "Se pega la flecha «{}» → «{}». Ya había una igual: este k-mero está repetido, y por eso salen dos flechas paralelas.",
                nombre_de_nodo(&grafo, desde),
                nombre_de_nodo(&grafo, hasta)
            )
        } else {
            format!(
                "Se pega la flecha «{}» → «{}», que representa al k-mero «{}». {}; {}.",
                nombre_de_nodo(&grafo, desde),
                nombre_de_nodo(&grafo, hasta),
                etiqueta_corta(kmero),
                estado(desde, nuevo_desde),
                estado(hasta, nuevo_hasta)
            )
        };
        let mut variables = base();
        variables.push(("grafo.nodos.len()", miles(nodos_ahora)));
        variables.push(("k-meros pegados", miles(numero + 1)));
        narracion.agregar(
            tabla_de_grafo(
                con_variables(&variables).secuencia(ventana.0, &marcas),
                "nodos del grafo hasta ahora",
                &grafo,
                grafo.nodos.len(),
                &[desde, hasta],
            )
            .paso("grafo.salidas[desde].push(hasta);", nota),
        );
    }
    Ok(grafo)
}

fn nodo_para_en_pasos<'a>(
    etiqueta: &'a [u8],
    indices: &mut HashMap<&'a [u8], usize>,
    grafo: &mut Grafo,
) -> usize {
    *indices.entry(etiqueta).or_insert_with(|| {
        grafo.nodos.push(etiqueta.to_vec());
        grafo.salidas.push(Vec::new());
        grafo.nodos.len() - 1
    })
}

fn pasos_de_bruijn_texto(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let texto = entrada.secuencia("texto")?;
    let k = leer_k(entrada, "k", 2, texto.len())?;
    let real = ensamblaje::de_bruijn_de_texto(&texto, k)?;
    let kmeros = ensamblaje::composicion(&texto, k);
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        con_variables(&[("k", k.to_string()), ("texto.len()", miles(texto.len()))])
            .secuencia(&texto, &[])
            .paso(
                "if k < 2 || k > texto.len() {",
                format!(
                    "El texto mide {} letras y k = {k}. Se revisa que k esté entre 2 y el largo del texto: sí. La meta: un grafo donde cada k-mero sea una flecha entre sus primeras {} letras y las últimas {}.",
                    miles(texto.len()),
                    k - 1,
                    k - 1
                ),
            ),
    );
    let primeros: Vec<String> = kmeros
        .iter()
        .take(FILAS_DE_PIZARRA)
        .map(|kmero| letras(kmero))
        .collect();
    narracion.agregar(
        con_variables(&[
            ("k", k.to_string()),
            ("texto.len()", miles(texto.len())),
            ("kmeros.len()", miles(kmeros.len())),
        ])
        .secuencia(&texto, &[(0, k, Resalte::Actual)])
        .lista("k-meros del texto (los primeros)", primeros, &[0])
        .paso(
            "de_bruijn_de_kmeros(&composicion(texto, k))",
            format!(
                "Primero el texto se parte en sus {} (la composición). Luego esa lista se convierte en grafo, un k-mero a la vez.",
                cantidad(kmeros.len(), "k-mero", "k-meros")
            ),
        ),
    );
    let grafo = narrar_nodos(&mut narracion, &kmeros, Some(&texto), 60)?;
    if grafo != real {
        return Err(DISTINTO_DEL_ALGORITMO.to_string());
    }
    let aristas = ensamblaje::cantidad_de_aristas(&grafo);
    Ok(narracion.terminar(
        con_variables(&[
            ("grafo.nodos.len()", miles(grafo.nodos.len())),
            ("aristas", miles(aristas)),
        ])
        .paso(
            "grafo.salidas[desde].push(hasta);",
            format!(
                "Listo: {} se pegaron como aristas entre {}. Hay menos nodos que aristas + 1 cuando un trozo de {} letras se repite en el texto. «Qué significa» te muestra los nodos que se ramifican.",
                cantidad(aristas, "k-mero", "k-meros"),
                cantidad(grafo.nodos.len(), "nodo", "nodos"),
                k - 1
            ),
        ),
    ))
}

fn pasos_de_bruijn_kmeros(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let kmeros = leer_kmeros(entrada, "kmeros")?;
    let lista = rebanadas(&kmeros);
    let real = ensamblaje::de_bruijn_de_kmeros(&lista)?;
    let k = lista[0].len();
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let primeros: Vec<String> = lista
        .iter()
        .take(FILAS_DE_PIZARRA)
        .map(|kmero| letras(kmero))
        .collect();
    narracion.agregar(
        con_variables(&[("kmeros.len()", miles(lista.len())), ("k", k.to_string())])
            .lista("k-meros (los primeros)", primeros, &[])
            .paso(
                "let k = largo_comun(kmeros)?;",
                format!(
                    "Hay {} y todos miden k = {k} letras. Están en desorden. La meta: que cada k-mero sea una flecha entre su prefijo y su sufijo, de {} letras cada uno, pegando los nodos que sean iguales.",
                    cantidad(lista.len(), "k-mero", "k-meros"),
                    k - 1
                ),
            ),
    );
    let grafo = narrar_nodos(&mut narracion, &lista, None, 60)?;
    if grafo != real {
        return Err(DISTINTO_DEL_ALGORITMO.to_string());
    }
    let aristas = ensamblaje::cantidad_de_aristas(&grafo);
    Ok(narracion.terminar(
        con_variables(&[
            ("grafo.nodos.len()", miles(grafo.nodos.len())),
            ("aristas", miles(aristas)),
        ])
        .paso(
            "grafo.salidas[desde].push(hasta);",
            format!(
                "Listo: {} y {}. Los k-meros repetidos dan aristas paralelas. En «Qué significa» se ve si el grafo tiene un camino que use todas las aristas.",
                cantidad(grafo.nodos.len(), "nodo", "nodos"),
                cantidad(aristas, "arista", "aristas")
            ),
        ),
    ))
}

fn pizarra_de_recorrido(
    grafo: &Grafo,
    inicio: usize,
    pila: &[usize],
    camino: &[usize],
    extras: &[(&str, String)],
) -> Pizarra {
    let elementos = nombres_de_nodos(grafo, cola_de_nodos(pila, FILAS_DE_PIZARRA));
    let en_camino = if camino.is_empty() {
        "(vacío)".to_string()
    } else {
        nodos_con_puntos(grafo, cola_de_nodos(camino, 4), camino.len())
    };
    let mut variables = vec![
        ("inicio", nombre_de_nodo(grafo, inicio)),
        ("pila.len()", miles(pila.len())),
        ("camino", en_camino),
    ];
    variables.extend(extras.iter().cloned());
    let resaltado = [elementos.len().saturating_sub(1)];
    let resaltar: &[usize] = if elementos.is_empty() {
        &[]
    } else {
        &resaltado
    };
    con_variables(&variables).lista("pila (el último es el nodo actual)", elementos, resaltar)
}

fn narrar_hierholzer(
    narracion: &mut Narracion,
    grafo: &Grafo,
    inicio: usize,
    cupo: usize,
) -> Vec<usize> {
    let aristas = ensamblaje::cantidad_de_aristas(grafo);
    let eventos = 2 * aristas + 2;
    let mut aristas_usadas = vec![0usize; grafo.salidas.len()];
    let mut pila = vec![inicio];
    let mut camino: Vec<usize> = Vec::with_capacity(aristas + 1);
    narracion.agregar(
        pizarra_de_recorrido(grafo, inicio, &pila, &camino, &[]).paso(
            "let mut pila = vec![inicio];",
            format!(
                "La pila empieza con el nodo inicial «{}» y el camino está vacío. La idea: avanzar por aristas que todavía no usaste hasta quedar atascado; entonces ese nodo se anota en el camino y se vuelve un paso atrás.",
                nombre_de_nodo(grafo, inicio)
            ),
        ),
    );
    let mut evento = 0;
    while let Some(&actual) = pila.last() {
        evento += 1;
        let usadas = aristas_usadas[actual];
        let salen = grafo.salidas[actual].len();
        if usadas < salen {
            let destino = grafo.salidas[actual][usadas];
            aristas_usadas[actual] += 1;
            pila.push(destino);
            if se_muestra(evento, eventos, cupo) {
                narracion.agregar(
                    pizarra_de_recorrido(
                        grafo,
                        inicio,
                        &pila,
                        &camino,
                        &[
                            ("actual", nombre_de_nodo(grafo, actual)),
                            ("usadas", usadas.to_string()),
                            ("destino", nombre_de_nodo(grafo, destino)),
                        ],
                    )
                    .paso(
                        "pila.push(destino);",
                        format!(
                            "«{}» todavía tiene aristas sin usar (lleva {usadas} de {salen}): se toma la siguiente, que va hacia «{}», y ese nodo entra a la pila.",
                            nombre_de_nodo(grafo, actual),
                            nombre_de_nodo(grafo, destino)
                        ),
                    ),
                );
            }
        } else {
            pila.pop();
            camino.push(actual);
            if se_muestra(evento, eventos, cupo) {
                let mut nota = format!(
                    "«{}» ya usó todas sus aristas ({salen}): sale de la pila y se anota en el camino, que se va armando al revés.",
                    nombre_de_nodo(grafo, actual)
                );
                match pila.last() {
                    None => nota.push_str(
                        " La pila quedó vacía: ya se recorrió todo lo que se podía alcanzar.",
                    ),
                    Some(&anterior) if aristas_usadas[anterior] < grafo.salidas[anterior].len() => {
                        nota.push_str(&format!(
                            " «{}» todavía tiene aristas sin usar: desde ahí se recorrerá otro ciclo que se empalma en este punto.",
                            nombre_de_nodo(grafo, anterior)
                        ));
                    }
                    Some(_) => {}
                }
                narracion.agregar(
                    pizarra_de_recorrido(
                        grafo,
                        inicio,
                        &pila,
                        &camino,
                        &[("actual", nombre_de_nodo(grafo, actual))],
                    )
                    .paso("camino.push(actual);", nota),
                );
            }
        }
    }
    camino.reverse();
    narracion.agregar(
        pizarra_de_recorrido(grafo, inicio, &pila, &camino, &[]).paso(
            "camino.reverse();",
            "Se da vuelta el camino: ahora va del nodo inicial al final y usa cada arista una sola vez.",
        ),
    );
    camino
}

fn tabla_de_balance(
    pizarra: Pizarra,
    grafo: &Grafo,
    entradas: &[usize],
    balance: &[i64],
) -> Pizarra {
    let desbalanceados: Vec<usize> = (0..balance.len())
        .filter(|&nodo| balance[nodo] != 0)
        .take(4)
        .collect();
    let titulo = if balance.len() > FILAS_DE_PIZARRA {
        "entran y salen de cada nodo (algunos nodos)"
    } else {
        "entran y salen de cada nodo"
    };
    tabla_de_grados(pizarra, titulo, grafo, entradas, &desbalanceados)
}

fn explicar_balance(balance: &[i64], grafo: &Grafo, para_ciclo: bool) -> String {
    let desbalanceados: Vec<usize> = (0..balance.len())
        .filter(|&nodo| balance[nodo] != 0)
        .collect();
    if desbalanceados.is_empty() {
        return if para_ciclo {
            "Para cada nodo se resta cuántas aristas entran de cuántas salen. En todos da 0: cada vez que el recorrido entra a un nodo puede salir, que es lo que pide un ciclo.".to_string()
        } else {
            "Para cada nodo se resta cuántas aristas entran de cuántas salen. En todos da 0: el recorrido puede volver a donde empezó.".to_string()
        };
    }
    let nombres: Vec<String> = desbalanceados
        .iter()
        .take(4)
        .map(|&nodo| format!("«{}» ({:+})", nombre_de_nodo(grafo, nodo), balance[nodo]))
        .collect();
    format!(
        "Para cada nodo se resta cuántas aristas entran de cuántas salen. No en todos da 0: {}{}.",
        nombres.join(", "),
        if desbalanceados.len() > 4 {
            format!(" y {} más", desbalanceados.len() - 4)
        } else {
            String::new()
        }
    )
}

fn pasos_ciclo(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let grafo = leer_grafo(entrada, "grafo")?;
    let real = ensamblaje::ciclo_euleriano(&grafo);
    let aristas = ensamblaje::cantidad_de_aristas(&grafo);
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        tabla_de_grafo(
            con_variables(&[
                ("grafo.nodos.len()", miles(grafo.nodos.len())),
                ("aristas", miles(aristas)),
            ]),
            "lista de adyacencia (algunos nodos)",
            &grafo,
            grafo.nodos.len(),
            &[],
        )
        .paso(
            "let aristas = cantidad_de_aristas(grafo);",
            format!(
                "El grafo tiene {} y {}. Un ciclo euleriano pasa por cada arista exactamente una vez y vuelve al nodo de partida. Primero se cuentan las aristas.",
                cantidad(grafo.nodos.len(), "nodo", "nodos"),
                cantidad(aristas, "arista", "aristas")
            ),
        ),
    );
    if aristas == 0 {
        let motivo = real.err().unwrap_or_default();
        return Ok(
            narracion.terminar(con_variables(&[("aristas", "0".to_string())]).paso(
                "if aristas == 0 {",
                format!("{motivo} Sin aristas no hay nada que recorrer: el algoritmo se detiene."),
            )),
        );
    }
    let entradas = ensamblaje::grados_de_entrada(&grafo);
    let balance = ensamblaje::balance_de_nodos(&grafo);
    narracion.agregar(
        tabla_de_balance(
            con_variables(&[("aristas", miles(aristas))]),
            &grafo,
            &entradas,
            &balance,
        )
        .paso(
            "let balance = balance_de_nodos(grafo);",
            explicar_balance(&balance, &grafo, true),
        ),
    );
    if balance.iter().any(|&diferencia| diferencia != 0) {
        let motivo = real.err().unwrap_or_default();
        return Ok(narracion.terminar(
            tabla_de_balance(
                con_variables(&[("aristas", miles(aristas))]),
                &grafo,
                &entradas,
                &balance,
            )
            .paso(
                "if balance.iter().any(|&diferencia| diferencia != 0) {",
                format!("{motivo} El algoritmo se detiene aquí. En «Qué significa» se dice qué nodos hay que equilibrar."),
            ),
        ));
    }
    let inicio = ensamblaje::primer_nodo_con_salidas(&grafo);
    narracion.agregar(
        con_variables(&[
            ("aristas", miles(aristas)),
            ("inicio", nombre_de_nodo(&grafo, inicio)),
        ])
        .paso(
            "let inicio = primer_nodo_con_salidas(grafo);",
            format!(
                "Como todos están balanceados, el ciclo puede empezar en cualquier nodo con salidas. Se elige el primero: «{}».",
                nombre_de_nodo(&grafo, inicio)
            ),
        ),
    );
    let ciclo = narrar_hierholzer(&mut narracion, &grafo, inicio, 50);
    terminar_recorrido(narracion, &grafo, aristas, ciclo, real, "ciclo")
}

fn terminar_recorrido(
    mut narracion: Narracion,
    grafo: &Grafo,
    aristas: usize,
    recorrido: Vec<usize>,
    real: Result<Vec<usize>, String>,
    nombre: &str,
) -> Result<Vec<Paso>, String> {
    let linea_de_revision = format!("if {nombre}.len() != aristas + 1 {{");
    let variables = [
        (
            if nombre == "ciclo" {
                "ciclo.len()"
            } else {
                "camino.len()"
            },
            miles(recorrido.len()),
        ),
        ("aristas", miles(aristas)),
    ];
    match real {
        Err(motivo) => {
            if recorrido.len() == aristas + 1 {
                return Err(DISTINTO_DEL_ALGORITMO.to_string());
            }
            Ok(narracion.terminar(con_variables(&variables).paso(
                &linea_de_revision,
                format!(
                    "El recorrido solo alcanzó {} nodos y hacían falta {}. {motivo}",
                    miles(recorrido.len()),
                    miles(aristas + 1)
                ),
            )))
        }
        Ok(esperado) => {
            if esperado != recorrido {
                return Err(DISTINTO_DEL_ALGORITMO.to_string());
            }
            narracion.agregar(con_variables(&variables).paso(
                &linea_de_revision,
                format!(
                    "El recorrido tiene {} nodos, que son las {} aristas + 1: se usaron todas, así que el grafo era una sola pieza.",
                    miles(recorrido.len()),
                    miles(aristas)
                ),
            ));
            let primeros: Vec<String> =
                nombres_de_nodos(grafo, &recorrido[..recorrido.len().min(FILAS_DE_PIZARRA)]);
            Ok(narracion.terminar(
                con_variables(&variables)
                    .lista(
                        &format!("{nombre} (los primeros nodos)"),
                        primeros,
                        &[0],
                    )
                    .paso(
                        &format!("Ok({nombre})"),
                        format!(
                            "Listo: {} nodos en el recorrido, empezando en «{}»: {}. En «Qué significa» se explica qué decisiones quedaron abiertas.",
                            miles(recorrido.len()),
                            nombre_de_nodo(grafo, recorrido[0]),
                            camino_corto(grafo, &recorrido)
                        ),
                    ),
            ))
        }
    }
}

fn pasos_camino(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let grafo = leer_grafo(entrada, "grafo")?;
    let real = ensamblaje::camino_euleriano(&grafo);
    let aristas = ensamblaje::cantidad_de_aristas(&grafo);
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        tabla_de_grafo(
            con_variables(&[
                ("grafo.nodos.len()", miles(grafo.nodos.len())),
                ("aristas", miles(aristas)),
            ]),
            "lista de adyacencia (algunos nodos)",
            &grafo,
            grafo.nodos.len(),
            &[],
        )
        .paso(
            "let aristas = cantidad_de_aristas(grafo);",
            format!(
                "El grafo tiene {} y {}. Un camino euleriano pasa por cada arista exactamente una vez, y puede terminar en un nodo distinto del que empezó. Primero se cuentan las aristas.",
                cantidad(grafo.nodos.len(), "nodo", "nodos"),
                cantidad(aristas, "arista", "aristas")
            ),
        ),
    );
    if aristas == 0 {
        let motivo = real.err().unwrap_or_default();
        return Ok(
            narracion.terminar(con_variables(&[("aristas", "0".to_string())]).paso(
                "if aristas == 0 {",
                format!("{motivo} Sin aristas no hay nada que recorrer: el algoritmo se detiene."),
            )),
        );
    }
    let entradas = ensamblaje::grados_de_entrada(&grafo);
    let balance = ensamblaje::balance_de_nodos(&grafo);
    narracion.agregar(
        tabla_de_balance(
            con_variables(&[("aristas", miles(aristas))]),
            &grafo,
            &entradas,
            &balance,
        )
        .paso(
            "let balance = balance_de_nodos(grafo);",
            format!(
                "{} En un camino puede haber un nodo con una salida de más (el comienzo) y uno con una entrada de más (el final).",
                explicar_balance(&balance, &grafo, false)
            ),
        ),
    );
    let mut inicio = None;
    let mut fin = None;
    let mut otros = 0;
    for (nodo, &diferencia) in balance.iter().enumerate() {
        if diferencia == 1 && inicio.is_none() {
            inicio = Some(nodo);
        } else if diferencia == -1 && fin.is_none() {
            fin = Some(nodo);
        } else if diferencia != 0 {
            otros += 1;
        }
    }
    if let Some(nodo) = inicio {
        narracion.agregar(
            tabla_de_balance(
                con_variables(&[("inicio", nombre_de_nodo(&grafo, nodo))]),
                &grafo,
                &entradas,
                &balance,
            )
            .paso(
                "if diferencia == 1 && inicio.is_none() {",
                format!(
                    "«{}» tiene una salida de más: es el candidato a comienzo, porque el camino sale de ahí una vez más de las que entra.",
                    nombre_de_nodo(&grafo, nodo)
                ),
            ),
        );
    }
    if let Some(nodo) = fin {
        narracion.agregar(
            tabla_de_balance(
                con_variables(&[("fin", nombre_de_nodo(&grafo, nodo))]),
                &grafo,
                &entradas,
                &balance,
            )
            .paso(
                "else if diferencia == -1 && fin.is_none() {",
                format!(
                    "«{}» tiene una entrada de más: es el candidato a final, porque el camino llega ahí una vez más de las que sale.",
                    nombre_de_nodo(&grafo, nodo)
                ),
            ),
        );
    }
    narracion.agregar(
        con_variables(&[("otros_desbalanceados", otros.to_string())]).paso(
            "if otros_desbalanceados > 0 {",
            if otros > 0 {
                "Hay nodos desbalanceados que no son ni el comienzo ni el final: así no se puede."
                    .to_string()
            } else {
                "No hay más nodos desbalanceados que el comienzo y el final: se puede seguir."
                    .to_string()
            },
        ),
    );
    if otros > 0 {
        let motivo = real.err().unwrap_or_default();
        return Ok(narracion.terminar(
            con_variables(&[("otros_desbalanceados", otros.to_string())]).paso(
                "if otros_desbalanceados > 0 {",
                format!("{motivo} El algoritmo se detiene aquí."),
            ),
        ));
    }
    let primero = ensamblaje::primer_nodo_con_salidas(&grafo);
    let comienzo = inicio.unwrap_or(primero);
    narracion.agregar(
        con_variables(&[("inicio", nombre_de_nodo(&grafo, comienzo))]).paso(
            "let inicio = inicio.unwrap_or_else(|| primer_nodo_con_salidas(grafo));",
            if inicio.is_some() {
                format!(
                    "Se empieza en «{}», el nodo con la salida de más.",
                    nombre_de_nodo(&grafo, comienzo)
                )
            } else {
                format!(
                    "Ningún nodo tiene la salida de más: el camino será un ciclo y puede empezar en el primer nodo con salidas, «{}».",
                    nombre_de_nodo(&grafo, comienzo)
                )
            },
        ),
    );
    let camino = narrar_hierholzer(&mut narracion, &grafo, comienzo, 50);
    terminar_recorrido(narracion, &grafo, aristas, camino, real, "camino")
}

fn pasos_universal(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let k = entrada.natural("k")?;
    if k == 0 || k > ensamblaje::MAXIMO_K_UNIVERSAL {
        return Err(format!(
            "k debe estar entre 1 y {}: con k más grande la cadena tendría más de un millón de letras.",
            ensamblaje::MAXIMO_K_UNIVERSAL
        ));
    }
    let real = ensamblaje::cadena_universal_circular(k)?;
    let combinaciones = 1usize << k;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        con_variables(&[("k", k.to_string())]).paso(
            "if k == 0 || k > MAXIMO_K_UNIVERSAL {",
            format!(
                "Con k = {k} existen 2^{k} = {} combinaciones de ceros y unos. La meta: una cadena circular, lo más corta posible, donde cada combinación aparezca exactamente una vez. k = {k} está dentro del límite (de 1 a {}).",
                miles(combinaciones),
                ensamblaje::MAXIMO_K_UNIVERSAL
            ),
        ),
    );
    if k == 1 {
        return Ok(narracion.terminar(
            con_variables(&[("k", "1".to_string()), ("cadena", "01".to_string())]).paso(
                "return Ok(b\"01\".to_vec());",
                "Con k = 1 las combinaciones son «0» y «1», y la cadena «01» ya las tiene a las dos: no hace falta armar ningún grafo.",
            ),
        ));
    }
    let cantidad_de_nodos = 1usize << (k - 1);
    narracion.agregar(
        con_variables(&[
            ("k", k.to_string()),
            ("cantidad", miles(cantidad_de_nodos)),
        ])
        .paso(
            "let cantidad = 1usize << (k - 1);",
            format!(
                "Los nodos del grafo serán las combinaciones de k − 1 = {} letras: 2^{} = {}. Cada arista será una combinación de {k} letras que une dos de esos nodos.",
                k - 1,
                k - 1,
                miles(cantidad_de_nodos)
            ),
        ),
    );
    let grafo = ensamblaje::de_bruijn_binario(k);
    for numero in 0..cantidad_de_nodos {
        let sin_el_primero = (numero * 2) % cantidad_de_nodos;
        if grafo.salidas[numero] != [sin_el_primero, sin_el_primero + 1] {
            return Err(DISTINTO_DEL_ALGORITMO.to_string());
        }
        if !se_muestra(numero, cantidad_de_nodos, 10) {
            continue;
        }
        narracion.agregar(
            tabla_de_grafo(
                con_variables(&[
                    ("k", k.to_string()),
                    ("cantidad", miles(cantidad_de_nodos)),
                    ("numero", numero.to_string()),
                    ("sin_el_primero", sin_el_primero.to_string()),
                ]),
                "nodos del grafo hasta ahora",
                &grafo,
                numero + 1,
                &[numero],
            )
            .paso(
                "let sin_el_primero = (numero * 2) % cantidad;",
                format!(
                    "El nodo {numero} se escribe «{}». Si le quitas la primera letra y le agregas un 0 o un 1 al final llegas a «{}» (nodo {sin_el_primero}) y a «{}» (nodo {}): son sus dos aristas.",
                    nombre_de_nodo(&grafo, numero),
                    nombre_de_nodo(&grafo, sin_el_primero),
                    nombre_de_nodo(&grafo, sin_el_primero + 1),
                    sin_el_primero + 1
                ),
            ),
        );
    }
    let entradas = ensamblaje::grados_de_entrada(&grafo);
    let balance = ensamblaje::balance_de_nodos(&grafo);
    narracion.agregar(
        tabla_de_balance(
            con_variables(&[("k", k.to_string())]),
            &grafo,
            &entradas,
            &balance,
        )
        .paso(
            "let ciclo = ciclo_euleriano(&grafo)?;",
            "Ahora se busca un ciclo euleriano: un recorrido que use cada arista una sola vez. Es posible porque en todos los nodos entran 2 aristas y salen 2. Cada arista usada es una combinación de k letras que aparece una vez.",
        ),
    );
    let inicio = ensamblaje::primer_nodo_con_salidas(&grafo);
    let ciclo = narrar_hierholzer(&mut narracion, &grafo, inicio, 40);
    if ciclo != ensamblaje::ciclo_euleriano(&grafo)? {
        return Err(DISTINTO_DEL_ALGORITMO.to_string());
    }
    let etiquetas = ensamblaje::etiquetas_del_camino(&grafo, &ciclo);
    let mut cadena = ensamblaje::deletrear_camino(&etiquetas)?;
    narracion.agregar(
        con_variables(&[
            ("k", k.to_string()),
            ("ciclo.len()", miles(ciclo.len())),
        ])
        .paso(
            "let mut cadena = deletrear_camino(&etiquetas_del_camino(&grafo, &ciclo))?;",
            format!(
                "El ciclo tiene {} nodos. Para convertirlo en letras se escribe el primer nodo y, de cada nodo siguiente, solo su última letra.",
                miles(ciclo.len())
            ),
        ),
    );
    narrar_deletreo(&mut narracion, &etiquetas, &cadena, 25, ("nodo", "nodos"))?;
    let antes = cadena.len();
    cadena.truncate(cadena.len() - (k - 1));
    if cadena != real {
        return Err(DISTINTO_DEL_ALGORITMO.to_string());
    }
    narracion.agregar(
        con_variables(&[
            ("k", k.to_string()),
            ("cadena.len()", miles(cadena.len())),
        ])
        .secuencia(&real, &[(0, real.len(), Resalte::Hallado)])
        .paso(
            "cadena.truncate(cadena.len() - (k - 1));",
            format!(
                "La cadena es circular: sus últimas {} letras ya son el comienzo de la vuelta siguiente, así que se cortan. Pasa de {} a {} letras.",
                k - 1,
                miles(antes),
                miles(real.len())
            ),
        ),
    );
    Ok(narracion.terminar(
        con_variables(&[
            ("k", k.to_string()),
            ("cadena.len()", miles(real.len())),
            ("cadena", etiqueta_corta(&real)),
        ])
        .secuencia(&real, &[(0, real.len(), Resalte::Hallado)])
        .paso(
            "Ok(cadena)",
            format!(
                "Listo: una cadena circular de {} letras en la que cada una de las {} combinaciones de {k} letras aparece una vez. «Qué significa» lo comprueba dando la vuelta.",
                miles(real.len()),
                miles(combinaciones)
            ),
        ),
    ))
}

fn narrar_union(
    narracion: &mut Narracion,
    primeras: &[&[u8]],
    segundas: &[&[u8]],
    d: usize,
    cupo: usize,
) -> Result<Vec<u8>, String> {
    let genoma = ensamblaje::unir_prefijo_y_sufijo(primeras, segundas, d)?;
    let prefijo = ensamblaje::deletrear_camino(primeras)?;
    let sufijo = ensamblaje::deletrear_camino(segundas)?;
    let total = primeras.len();
    let k = primeras[0].len();
    narrar_deletreo(narracion, primeras, &prefijo, cupo, ("lectura", "lecturas"))?;
    narracion.agregar(
        con_variables(&[
            ("primeras.len()", miles(total)),
            ("prefijo.len()", miles(prefijo.len())),
            ("prefijo", etiqueta_corta(&prefijo)),
        ])
        .secuencia(&prefijo, &[(0, prefijo.len(), Resalte::Hallado)])
        .paso(
            "let prefijo = deletrear_camino(primeras)?;",
            format!(
                "Pegadas las primeras lecturas queda la cadena «{}» de {} letras.",
                etiqueta_corta(&prefijo),
                miles(prefijo.len())
            ),
        ),
    );
    narracion.agregar(
        con_variables(&[
            ("primeras.len()", miles(total)),
            ("prefijo.len()", miles(prefijo.len())),
            ("sufijo.len()", miles(sufijo.len())),
            ("sufijo", etiqueta_corta(&sufijo)),
        ])
        .secuencia(&sufijo, &[(0, sufijo.len(), Resalte::Hallado)])
        .paso(
            "let sufijo = deletrear_camino(segundas)?;",
            format!(
                "Con las segundas lecturas se hace lo mismo: queda «{}» de {} letras. Esta cadena va más a la derecha que la anterior.",
                etiqueta_corta(&sufijo),
                miles(sufijo.len())
            ),
        ),
    );
    let salto = k + d;
    narracion.agregar(
        con_variables(&[
            ("k", k.to_string()),
            ("d", d.to_string()),
            ("salto", salto.to_string()),
        ])
        .paso(
            "let salto = k.saturating_add(d);",
            format!(
                "Cada lectura mide {k} y entre las dos hay {d} letras sin leer: la cadena de las segundas empieza {k} + {d} = {salto} letras después de la de las primeras. Desde ahí las dos se superponen y deben decir lo mismo."
            ),
        ),
    );
    let comparables = prefijo.len() - salto;
    for (numero, &letra) in sufijo.iter().enumerate().take(comparables) {
        if !se_muestra(numero, comparables, 20) {
            continue;
        }
        let posicion = salto + numero;
        if prefijo[posicion] != letra {
            return Err(DISTINTO_DEL_ALGORITMO.to_string());
        }
        narracion.agregar(
            con_variables(&[
                ("salto", salto.to_string()),
                ("posicion", posicion.to_string()),
                ("prefijo[posicion]", (prefijo[posicion] as char).to_string()),
                (
                    "sufijo[posicion - salto]",
                    (letra as char).to_string(),
                ),
            ])
            .secuencia(&prefijo, &[(posicion, posicion + 1, Resalte::Actual)])
            .paso(
                "if prefijo[posicion] != sufijo[posicion - salto] {",
                format!(
                    "En la posición {posicion} la cadena de las primeras lecturas dice «{}» y la de las segundas, en su posición {numero}, también dice «{}»: coinciden.",
                    prefijo[posicion] as char,
                    letra as char
                ),
            ),
        );
    }
    if comparables == 0 {
        narracion.agregar(con_variables(&[("salto", salto.to_string())]).paso(
            "for posicion in salto..prefijo.len() {",
            "Las dos cadenas no se superponen: no hay nada que comparar.",
        ));
    }
    let ya_cubierto = prefijo.len() - salto;
    if ya_cubierto > sufijo.len() || prefijo.len() + sufijo.len() - ya_cubierto != genoma.len() {
        return Err(DISTINTO_DEL_ALGORITMO.to_string());
    }
    narracion.agregar(
        con_variables(&[
            ("genoma.len()", miles(genoma.len())),
            ("ya_cubierto", miles(ya_cubierto)),
        ])
        .secuencia(
            &genoma,
            &[
                (0, prefijo.len(), Resalte::Comparado),
                (prefijo.len(), genoma.len(), Resalte::Hallado),
            ],
        )
        .paso(
            "genoma.extend_from_slice(&sufijo[ya_cubierto..]);",
            format!(
                "A las {} letras de las primeras lecturas se les agregan las que solo aportan las segundas, a partir de la letra {} de su cadena: {} más.",
                miles(prefijo.len()),
                miles(ya_cubierto + 1),
                miles(genoma.len() - prefijo.len())
            ),
        ),
    );
    Ok(genoma)
}

fn pasos_camino_con_huecos(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let pares = leer_pares(entrada, "pares")?;
    let d = leer_d(entrada)?;
    let real = ensamblaje::deletrear_camino_con_huecos(&pares, d)?;
    let k = pares[0].primera.len();
    let primeras: Vec<&[u8]> = pares.iter().map(|par| par.primera.as_slice()).collect();
    let segundas: Vec<&[u8]> = pares.iter().map(|par| par.segunda.as_slice()).collect();
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let mostrados: Vec<String> = pares
        .iter()
        .take(FILAS_DE_PIZARRA)
        .map(|par| format!("{}|{}", letras(&par.primera), letras(&par.segunda)))
        .collect();
    narracion.agregar(
        con_variables(&[
            ("pares.len()", miles(pares.len())),
            ("k", k.to_string()),
            ("d", d.to_string()),
        ])
        .lista("pares (los primeros)", mostrados, &[0])
        .paso(
            "let primeras: Vec<&[u8]> = pares.iter().map(|par| par.primera.as_slice()).collect();",
            format!(
                "Hay {} de lecturas de {k} letras; en cada par, la segunda lectura queda a {d} letras de la primera. Se separan las primeras de las segundas: las primeras, en orden, forman una cadena; las segundas, otra.",
                cantidad(pares.len(), "par", "pares")
            ),
        ),
    );
    let genoma = narrar_union(&mut narracion, &primeras, &segundas, d, 14)?;
    if genoma != real {
        return Err(DISTINTO_DEL_ALGORITMO.to_string());
    }
    Ok(narracion.terminar(
        con_variables(&[
            ("genoma.len()", miles(real.len())),
            ("genoma", etiqueta_corta(&real)),
        ])
        .secuencia(&real, &[(0, real.len(), Resalte::Hallado)])
        .paso(
            "Ok(genoma)",
            format!(
                "Listo: la cadena completa tiene {} letras y todas las lecturas encajaron. En «Qué significa» se cuenta qué parte quedó cubierta por las dos mitades.",
                miles(real.len())
            ),
        ),
    ))
}

fn pasos_n50(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let longitudes = leer_longitudes(entrada, "contigs")?;
    let real = ensamblaje::medir_n50(&longitudes);
    let mut ordenadas = longitudes.clone();
    ordenadas.sort_unstable_by(|una, otra| otra.cmp(una));
    let total = ordenadas.len();
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let primeras = |desde: usize| -> Vec<String> {
        ordenadas[desde..(desde + FILAS_DE_PIZARRA).min(total)]
            .iter()
            .map(|&largo| miles(largo))
            .collect()
    };
    narracion.agregar(
        con_variables(&[("longitudes.len()", miles(total))]).paso(
            "pub fn medir_n50(longitudes: &[usize]) -> MedidaN50 {",
            format!(
                "Hay {} con sus largos. El N50 responde: ¿cuál es el largo del contig más corto entre los más largos que, juntos, suman la mitad de todo el ensamblaje?",
                cantidad(total, "contig", "contigs")
            ),
        ),
    );
    narracion.agregar(
        con_variables(&[("ordenadas.len()", miles(total))])
            .lista("largos de mayor a menor (los primeros)", primeras(0), &[0])
            .paso(
                "ordenadas.sort_unstable_by(|una, otra| otra.cmp(una));",
                "Se ordenan los largos de mayor a menor: los contigs más largos van primero.",
            ),
    );
    let mut suma: usize = 0;
    for (posicion, &largo) in ordenadas.iter().enumerate() {
        suma = suma.saturating_add(largo);
        if se_muestra(posicion, total, 6) {
            narracion.agregar(
                con_variables(&[
                    ("total", miles(suma)),
                    ("longitud", miles(largo)),
                    ("posicion", posicion.to_string()),
                ])
                .lista(
                    &format!("largos de mayor a menor (desde el n.º {})", posicion + 1),
                    primeras(posicion),
                    &[0],
                )
                .paso(
                    "total = total.saturating_add(longitud);",
                    format!(
                        "Se suma el largo {} al total: lleva {}.",
                        miles(largo),
                        miles(suma)
                    ),
                ),
            );
        }
    }
    if suma != real.total {
        return Err(DISTINTO_DEL_ALGORITMO.to_string());
    }
    narracion.agregar(
        con_variables(&[("total", miles(real.total))]).paso(
            "let mut medida = MedidaN50 {",
            format!(
                "Todos los contigs suman {} letras. La mitad es {}: hay que ver cuántos contigs largos hacen falta para llegar a ella.",
                miles(real.total),
                miles(real.total / 2)
            ),
        ),
    );
    let mut acumulado: usize = 0;
    let mut encontrado = None;
    for (posicion, &largo) in ordenadas.iter().enumerate() {
        acumulado = acumulado.saturating_add(largo);
        let llega = real.total > 0 && acumulado >= real.total - acumulado;
        if llega {
            encontrado = Some((posicion, largo));
        }
        if llega || se_muestra(posicion, real.l50.max(1), 25) {
            let desde = posicion.saturating_sub(3);
            narracion.agregar(
                con_variables(&[
                    ("total", miles(real.total)),
                    ("acumulado", miles(acumulado)),
                    ("longitud", miles(largo)),
                    ("posicion", posicion.to_string()),
                ])
                .lista(
                    &format!("largos de mayor a menor (desde el n.º {})", desde + 1),
                    primeras(desde),
                    &[posicion - desde],
                )
                .paso(
                    "if total > 0 && acumulado >= total - acumulado {",
                    if llega {
                        format!(
                            "Los {} contigs más largos suman {}, y lo que queda es {}. Ya llegamos a la mitad (acumulado ≥ lo que falta): el último que sumamos, de {} letras, es el N50.",
                            miles(posicion + 1),
                            miles(acumulado),
                            miles(real.total - acumulado),
                            miles(largo)
                        )
                    } else {
                        format!(
                            "Con los {} más largos llevamos {} y faltan {}: todavía no alcanza la mitad, se sigue con el siguiente.",
                            miles(posicion + 1),
                            miles(acumulado),
                            miles(real.total - acumulado)
                        )
                    },
                ),
            );
        }
        if llega {
            break;
        }
    }
    let (posicion, largo) = encontrado.unwrap_or((0, 0));
    if largo != real.n50 || posicion + 1 != real.l50 {
        return Err(DISTINTO_DEL_ALGORITMO.to_string());
    }
    Ok(narracion.terminar(
        con_variables(&[
            ("n50", miles(real.n50)),
            ("l50", miles(real.l50)),
            ("total", miles(real.total)),
        ])
        .paso(
            "medida.n50 = longitud;",
            format!(
                "Listo: N50 = {} letras y L50 = {}, es decir, bastan {} para cubrir la mitad del ensamblaje. Un N50 grande significa pocos trozos largos. En «Qué significa» se compara con el ensamblaje ideal.",
                miles(real.n50),
                miles(real.l50),
                cantidad(real.l50, "contig", "contigs")
            ),
        ),
    ))
}

fn pasos_pares(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let lista = leer_pares(entrada, "pares")?;
    let d = leer_d(entrada)?;
    let semilla = entrada.natural("semilla")? as u64;
    let real = ensamblaje::reconstruir_desde_pares(&lista, d, semilla)?;
    let k = lista[0].primera.len();
    let total = lista.len();
    let mut grafo = ensamblaje::de_bruijn_de_pares(&lista)?;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let mostrados: Vec<String> = lista
        .iter()
        .take(FILAS_DE_PIZARRA)
        .map(|par| format!("{}|{}", letras(&par.primera), letras(&par.segunda)))
        .collect();
    narracion.agregar(
        con_variables(&[
            ("pares.len()", miles(total)),
            ("d", d.to_string()),
            ("semilla", semilla.to_string()),
        ])
        .lista("pares (los primeros)", mostrados, &[0])
        .paso(
            "let mut grafo = de_bruijn_de_pares(pares)?;",
            format!(
                "Hay {} de lecturas de {k} letras, y la segunda queda a {d} letras de la primera. La meta: reconstruir la cadena de la que salieron. Primero se arma un grafo de De Bruijn cuyos nodos son parejas de trozos de {} letras, uno de cada lectura.",
                cantidad(total, "par", "pares"),
                k - 1
            ),
        ),
    );
    let mut parcial = Grafo::default();
    let mut indices: HashMap<Vec<u8>, usize> = HashMap::new();
    for (numero, par) in lista.iter().enumerate() {
        let prefijo = ensamblaje::unir_par(&par.primera[..k - 1], &par.segunda[..k - 1]);
        let sufijo = ensamblaje::unir_par(&par.primera[1..], &par.segunda[1..]);
        let mostrar = se_muestra(numero, total, 20);
        let nodos_antes = parcial.nodos.len();
        let desde = ubicar_pareja(&prefijo, &mut indices, &mut parcial);
        let hasta = ubicar_pareja(&sufijo, &mut indices, &mut parcial);
        parcial.salidas[desde].push(hasta);
        if !mostrar {
            continue;
        }
        let par_texto = format!("{}|{}", letras(&par.primera), letras(&par.segunda));
        let marcas = [
            (0, k, Resalte::Actual),
            (k + 1, 2 * k + 1, Resalte::Comparado),
        ];
        let base = vec![
            ("k", k.to_string()),
            ("par", par_texto.clone()),
            ("prefijo", etiqueta_corta(&prefijo)),
            ("sufijo", etiqueta_corta(&sufijo)),
        ];
        if numero < 2 {
            narracion.agregar(
                con_variables(&base[..3])
                    .secuencia(par_texto.as_bytes(), &marcas)
                    .paso(
                        "prefijos.push(unir_par(&par.primera[..k - 1], &par.segunda[..k - 1]));",
                        format!(
                            "Del par «{par_texto}» se toman las primeras {} letras de cada lectura y se juntan: «{}». Ese será el nodo de partida de la arista de este par.",
                            k - 1,
                            etiqueta_corta(&prefijo)
                        ),
                    ),
            );
            narracion.agregar(
                con_variables(&base)
                    .secuencia(par_texto.as_bytes(), &marcas)
                    .paso(
                        "sufijos.push(unir_par(&par.primera[1..], &par.segunda[1..]));",
                        format!(
                            "Y las últimas {} letras de cada lectura dan «{}», el nodo de llegada.",
                            k - 1,
                            etiqueta_corta(&sufijo)
                        ),
                    ),
            );
        }
        let mut variables = base.clone();
        variables.push(("desde", desde.to_string()));
        variables.push(("hasta", hasta.to_string()));
        let nuevos = parcial.nodos.len() - nodos_antes;
        narracion.agregar(
            tabla_de_grafo(
                con_variables(&variables).secuencia(par_texto.as_bytes(), &marcas),
                "nodos del grafo hasta ahora (cada nodo es una pareja)",
                &parcial,
                parcial.nodos.len(),
                &[desde, hasta],
            )
            .paso(
                "grafo.salidas[desde].push(hasta);",
                format!(
                    "Se pega la arista «{}» → «{}». {}",
                    nombre_de_nodo(&parcial, desde),
                    nombre_de_nodo(&parcial, hasta),
                    match nuevos {
                        0 => "Los dos nodos ya estaban.".to_string(),
                        1 => "Uno de los nodos es nuevo.".to_string(),
                        _ => "Los dos nodos son nuevos.".to_string(),
                    }
                ),
            ),
        );
    }
    if parcial != grafo {
        return Err(DISTINTO_DEL_ALGORITMO.to_string());
    }
    let aristas = ensamblaje::cantidad_de_aristas(&grafo);
    narracion.agregar(
        con_variables(&[
            ("pares.len()", miles(total)),
            ("d", d.to_string()),
        ])
        .paso(
            "if pares.len() < d.saturating_add(1) {",
            format!(
                "Con {} y d = {d} el hueco entre las dos lecturas queda cubierto: hacen falta al menos d + 1 = {} pares y hay de sobra.",
                cantidad(total, "par", "pares"),
                d + 1
            ),
        ),
    );
    let entradas = ensamblaje::grados_de_entrada(&grafo);
    let balance = ensamblaje::balance_de_nodos(&grafo);
    narracion.agregar(
        tabla_de_balance(
            con_variables(&[
                ("grafo.nodos.len()", miles(grafo.nodos.len())),
                ("aristas", miles(aristas)),
            ]),
            &grafo,
            &entradas,
            &balance,
        )
        .paso(
            "camino_euleriano(&grafo)?;",
            format!(
                "Se comprueba que el grafo de {} y {} admite un camino que use cada arista una vez. Es lo que debe pasar si las lecturas vienen de una sola cadena.",
                cantidad(grafo.nodos.len(), "nodo", "nodos"),
                cantidad(aristas, "arista", "aristas")
            ),
        ),
    );
    let inicio_fijo = balance.iter().position(|&diferencia| diferencia == 1);
    let candidatos: Vec<usize> = (0..grafo.salidas.len())
        .filter(|&nodo| !grafo.salidas[nodo].is_empty())
        .collect();
    narracion.agregar(
        con_variables(&[
            ("inicio_fijo", match inicio_fijo {
                Some(nodo) => nombre_de_nodo(&grafo, nodo),
                None => "None".to_string(),
            }),
            ("candidatos.len()", miles(candidatos.len())),
        ])
        .paso(
            "let inicio_fijo = balance.iter().position(|&diferencia| diferencia == 1);",
            match inicio_fijo {
                Some(nodo) => format!(
                    "Hay un nodo con una salida de más, «{}»: el camino tiene que empezar ahí.",
                    nombre_de_nodo(&grafo, nodo)
                ),
                None => format!(
                    "Ningún nodo tiene una salida de más: el recorrido es un ciclo y puede empezar en cualquiera de los {}. Se irán probando.",
                    cantidad(candidatos.len(), "nodo con aristas", "nodos con aristas")
                ),
            },
        ),
    );
    let mut azar = Azar::nuevo(semilla);
    let desplazamiento = azar.entero(candidatos.len());
    narracion.agregar(
        con_variables(&[
            ("semilla", semilla.to_string()),
            ("desplazamiento", desplazamiento.to_string()),
        ])
        .paso(
            "let mut azar = Azar::nuevo(semilla);",
            format!(
                "Si el primer camino no deletrea una cadena coherente (porque el grafo tiene ramas y hay varios caminos posibles), se baraja el orden de las salidas y se prueba otro, hasta {} intentos. La semilla {semilla} hace que el barajado sea siempre el mismo: el mismo resultado cada vez.",
                ensamblaje::INTENTOS_CON_PARES
            ),
        ),
    );
    let mut intentos = 0;
    let mut motivo = String::new();
    let mut armado: Option<Vec<u8>> = None;
    while intentos < ensamblaje::INTENTOS_CON_PARES {
        let inicio = match inicio_fijo {
            Some(nodo) => nodo,
            None => candidatos[(intentos + desplazamiento) % candidatos.len()],
        };
        intentos += 1;
        let camino = ensamblaje::recorrido_de_hierholzer(&grafo, inicio);
        let completo = camino.len() == aristas + 1;
        let mut resultado: Option<Result<Vec<u8>, String>> = None;
        if completo {
            let mut primeras = Vec::with_capacity(camino.len());
            let mut segundas = Vec::with_capacity(camino.len());
            for &nodo in &camino {
                let (primera, segunda) = ensamblaje::separar_par(&grafo.nodos[nodo]);
                primeras.push(primera);
                segundas.push(segunda);
            }
            resultado = Some(ensamblaje::unir_prefijo_y_sufijo(
                &primeras,
                &segundas,
                d + 1,
            ));
        }
        let exito = matches!(resultado, Some(Ok(_)));
        if intentos == 1 || exito {
            let rehecho = narrar_hierholzer(&mut narracion, &grafo, inicio, 30);
            if rehecho != camino {
                return Err(DISTINTO_DEL_ALGORITMO.to_string());
            }
        }
        let numero_de_intento = intentos.to_string();
        let variables = [
            ("intentos", numero_de_intento),
            ("inicio", nombre_de_nodo(&grafo, inicio)),
            ("camino.len()", miles(camino.len())),
            ("aristas + 1", miles(aristas + 1)),
        ];
        if let Some(Ok(genoma)) = resultado {
            if intentos != real.intentos || genoma != real.genoma {
                return Err(DISTINTO_DEL_ALGORITMO.to_string());
            }
            if intentos == 1 {
                narracion.agregar(
                    con_variables(&variables).paso(
                        "if camino.len() == aristas + 1 {",
                        format!(
                            "El camino usó las {} aristas ({} nodos): está completo. Cada nodo del camino es una pareja de trozos; con las primeras partes y las segundas partes se intenta deletrear la cadena.",
                            miles(aristas),
                            miles(camino.len())
                        ),
                    ),
                );
            } else {
                narracion.agregar(
                    con_variables(&variables).paso(
                        "if camino.len() == aristas + 1 {",
                        format!(
                            "Intento {intentos}: el camino está completo ({} nodos) y esta vez las dos mitades encajan.",
                            miles(camino.len())
                        ),
                    ),
                );
            }
            let (primeras_textos, segundas_textos): (Vec<String>, Vec<String>) = camino
                .iter()
                .take(FILAS_DE_PIZARRA)
                .map(|&nodo| {
                    let (primera, segunda) = ensamblaje::separar_par(&grafo.nodos[nodo]);
                    (letras(primera), letras(segunda))
                })
                .unzip();
            let lineas: Vec<String> = primeras_textos
                .iter()
                .zip(&segundas_textos)
                .map(|(primera, segunda)| format!("{primera} | {segunda}"))
                .collect();
            let primer_nodo = camino[0];
            let (primera, segunda) = ensamblaje::separar_par(&grafo.nodos[primer_nodo]);
            narracion.agregar(
                con_variables(&[
                    ("camino.len()", miles(camino.len())),
                    ("primera", letras(primera)),
                    ("segunda", letras(segunda)),
                ])
                .lista("nodos del camino, partidos en (primera | segunda)", lineas, &[0])
                .paso(
                    "let (primera, segunda) = separar_par(&grafo.nodos[nodo]);",
                    format!(
                        "Cada nodo del camino es una pareja; se parte en la mitad que viene de la primera lectura y la que viene de la segunda. Por ejemplo, «{}» se parte en «{}» y «{}».",
                        nombre_de_nodo(&grafo, primer_nodo),
                        letras(primera),
                        letras(segunda)
                    ),
                ),
            );
            let mut primeras = Vec::with_capacity(camino.len());
            let mut segundas = Vec::with_capacity(camino.len());
            for &nodo in &camino {
                let (primera, segunda) = ensamblaje::separar_par(&grafo.nodos[nodo]);
                primeras.push(primera);
                segundas.push(segunda);
            }
            narracion.agregar(
                con_variables(&[
                    ("d", d.to_string()),
                    ("d + 1", (d + 1).to_string()),
                    ("primeras.len()", miles(primeras.len())),
                ])
                .paso(
                    "match unir_prefijo_y_sufijo(&primeras, &segundas, d + 1) {",
                    format!(
                        "Los trozos de las parejas miden {} letras, una menos que las lecturas, así que el hueco entre ellos es d + 1 = {}. Se pegan las primeras partes, se pegan las segundas y se comprueba que se superpongan bien.",
                        k - 1,
                        d + 1
                    ),
                ),
            );
            let genoma = narrar_union(&mut narracion, &primeras, &segundas, d + 1, 12)?;
            if genoma != real.genoma {
                return Err(DISTINTO_DEL_ALGORITMO.to_string());
            }
            armado = Some(genoma);
            break;
        }
        if let Some(Err(error)) = resultado {
            motivo = error;
        }
        if se_muestra(intentos - 1, real.intentos, 10) {
            narracion.agregar(
                con_variables(&variables).paso(
                    "if camino.len() == aristas + 1 {",
                    if completo {
                        format!(
                            "Intento {intentos}: el camino está completo, pero al deletrear la cadena las dos mitades no encajan ({motivo}). Hay que probar otro camino."
                        )
                    } else {
                        format!(
                            "Intento {intentos}: el camino solo alcanzó {} de los {} nodos. Hay que probar otro."
                            , miles(camino.len()), miles(aristas + 1)
                        )
                    },
                ),
            );
        }
        if !ensamblaje::hay_bifurcaciones(&grafo)
            && (inicio_fijo.is_some() || candidatos.len() == 1)
        {
            break;
        }
        ensamblaje::barajar_salidas(&mut grafo, &mut azar);
        if se_muestra(intentos - 1, real.intentos, 10) {
            narracion.agregar(
                con_variables(&[("intentos", intentos.to_string())]).paso(
                    "barajar_salidas(&mut grafo, &mut azar);",
                    "En cada nodo se baraja el orden de sus salidas: el siguiente intento recorrerá las ramas en otro orden.",
                ),
            );
        }
    }
    let Some(genoma) = armado else {
        return Err(DISTINTO_DEL_ALGORITMO.to_string());
    };
    Ok(narracion.terminar(
        con_variables(&[
            ("intentos", intentos.to_string()),
            ("genoma.len()", miles(genoma.len())),
            ("genoma", etiqueta_corta(&genoma)),
        ])
        .secuencia(&genoma, &[(0, genoma.len(), Resalte::Hallado)])
        .paso(
            "Ok(ArmadoConPares { genoma, intentos })",
            format!(
                "Listo: la cadena tiene {} letras y salió en el intento {intentos}. En «Qué significa» se compara con lo que se obtiene usando las lecturas sueltas, sin parejas.",
                miles(genoma.len())
            ),
        ),
    ))
}

fn ubicar_pareja(
    etiqueta: &[u8],
    indices: &mut HashMap<Vec<u8>, usize>,
    grafo: &mut Grafo,
) -> usize {
    *indices.entry(etiqueta.to_vec()).or_insert_with(|| {
        grafo.nodos.push(etiqueta.to_vec());
        grafo.salidas.push(Vec::new());
        grafo.nodos.len() - 1
    })
}

fn narrar_caminos(narracion: &mut Narracion, grafo: &Grafo, caminos: &[Vec<usize>], cupo: usize) {
    let entradas = ensamblaje::grados_de_entrada(grafo);
    let cantidad_de_nodos = grafo.nodos.len();
    let de_paso: Vec<bool> = (0..cantidad_de_nodos)
        .map(|nodo| entradas[nodo] == 1 && grafo.salidas[nodo].len() == 1)
        .collect();
    let cuantos_de_paso = de_paso.iter().filter(|&&es| es).count();
    let puntos_de_corte: Vec<usize> = (0..cantidad_de_nodos)
        .filter(|&nodo| !de_paso[nodo])
        .collect();
    let primeros_cortes: Vec<usize> = puntos_de_corte.iter().copied().take(4).collect();
    narracion.agregar(
        tabla_de_grados(
            con_variables(&[("grafo.salidas.len()", miles(cantidad_de_nodos))]),
            "entran y salen de cada nodo (algunos nodos)",
            grafo,
            &entradas,
            &primeros_cortes,
        )
        .paso(
            "let entradas = grados_de_entrada(grafo);",
            "Se cuenta cuántas aristas entran a cada nodo; cuántas salen ya se sabe por la lista de flechas. La tabla muestra algunos nodos con sus dos números.",
        ),
    );
    narracion.agregar(
        tabla_de_grados(
            con_variables(&[
                ("cantidad", miles(cantidad_de_nodos)),
                ("nodos de paso", miles(cuantos_de_paso)),
            ]),
            "entran y salen de cada nodo (en las filas marcadas hay que decidir)",
            grafo,
            &entradas,
            &primeros_cortes,
        )
        .paso(
            "uno_a_uno[nodo] = entradas[nodo] == 1 && grafo.salidas[nodo].len() == 1;",
            format!(
                "Un nodo es «de paso» si le entra una sola arista y le sale una sola: ahí no hay nada que decidir. De {} nodos, {} son de paso y {} son puntos de corte (ramificados o extremos): de ahí empiezan y ahí terminan los caminos.",
                miles(cantidad_de_nodos),
                miles(cuantos_de_paso),
                miles(puntos_de_corte.len())
            ),
        ),
    );
    let total = caminos.len();
    for (numero, camino) in caminos.iter().enumerate() {
        if !se_muestra(numero, total, cupo) {
            continue;
        }
        let desde = (numero + 1).saturating_sub(6);
        let lineas: Vec<String> = caminos[desde..=numero]
            .iter()
            .map(|otro| camino_corto(grafo, otro))
            .collect();
        let resaltado = lineas.len() - 1;
        let inicio = camino[0];
        let fin = camino[camino.len() - 1];
        if de_paso[inicio] {
            narracion.agregar(
                con_variables(&[
                    ("nodo", nombre_de_nodo(grafo, inicio)),
                    ("actual", nombre_de_nodo(grafo, fin)),
                    ("caminos.len()", miles(numero + 1)),
                ])
                .lista("caminos encontrados (los últimos)", lineas, &[resaltado])
                .paso(
                    "caminos.push(ciclo);",
                    format!(
                        "Todos los nodos de este grupo son de paso: al seguir la flecha se da la vuelta y se regresa a «{}». Es un ciclo cerrado de {} que nadie conecta con el resto del grafo.",
                        nombre_de_nodo(grafo, inicio),
                        cantidad(camino.len() - 1, "nodo", "nodos")
                    ),
                ),
            );
            continue;
        }
        let intermedios = camino.len() - 2;
        let entran = entradas[inicio];
        let salen = grafo.salidas[inicio].len();
        let por_que = if entran == 0 {
            "no le entra ninguna arista: es un comienzo".to_string()
        } else if salen == 0 {
            "no le sale ninguna".to_string()
        } else {
            format!("le entran {entran} y le salen {salen}: es un cruce")
        };
        narracion.agregar(
            con_variables(&[
                ("nodo", nombre_de_nodo(grafo, inicio)),
                ("siguiente", nombre_de_nodo(grafo, camino[1])),
                ("actual", nombre_de_nodo(grafo, fin)),
                ("caminos.len()", miles(numero + 1)),
            ])
            .lista("caminos encontrados (los últimos)", lineas, &[resaltado])
            .paso(
                "caminos.push(camino);",
                format!(
                    "«{}» {por_que}. Desde ahí se avanza por {} de paso hasta «{}», donde el camino se acaba o se ramifica: camino número {}.",
                    nombre_de_nodo(grafo, inicio),
                    cantidad(intermedios, "nodo", "nodos"),
                    nombre_de_nodo(grafo, fin),
                    miles(numero + 1)
                ),
            ),
        );
    }
}

fn pasos_no_ramificados(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let grafo = leer_grafo(entrada, "grafo")?;
    let caminos = ensamblaje::caminos_no_ramificados(&grafo);
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        tabla_de_grafo(
            con_variables(&[
                ("grafo.salidas.len()", miles(grafo.nodos.len())),
                ("aristas", miles(ensamblaje::cantidad_de_aristas(&grafo))),
            ]),
            "lista de adyacencia (algunos nodos)",
            &grafo,
            grafo.nodos.len(),
            &[],
        )
        .paso(
            "pub fn caminos_no_ramificados(grafo: &Grafo) -> Vec<Vec<usize>> {",
            format!(
                "El grafo tiene {} y {}. La meta: encontrar los tramos más largos que se pueden recorrer sin tener que elegir nunca entre dos salidas.",
                cantidad(grafo.nodos.len(), "nodo", "nodos"),
                cantidad(ensamblaje::cantidad_de_aristas(&grafo), "arista", "aristas")
            ),
        ),
    );
    narrar_caminos(&mut narracion, &grafo, &caminos, 40);
    let aristas_cubiertas: usize = caminos.iter().map(|camino| camino.len() - 1).sum();
    if aristas_cubiertas != ensamblaje::cantidad_de_aristas(&grafo) {
        return Err(DISTINTO_DEL_ALGORITMO.to_string());
    }
    let hay_ciclos = caminos
        .iter()
        .any(|camino| camino[0] == camino[camino.len() - 1]);
    Ok(narracion.terminar(
        con_variables(&[
            ("caminos.len()", miles(caminos.len())),
            ("aristas", miles(aristas_cubiertas)),
        ])
        .paso(
            if hay_ciclos {
                "caminos.push(ciclo);"
            } else {
                "caminos.push(camino);"
            },
            format!(
                "Listo: {} que, juntos, cubren las {} del grafo, cada una una sola vez. En «Qué significa» se explica por qué cada uno es un trozo del que no hay duda.",
                cantidad(caminos.len(), "camino máximo", "caminos máximos"),
                cantidad(aristas_cubiertas, "arista", "aristas")
            ),
        ),
    ))
}

fn pasos_contigs(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let kmeros = leer_kmeros(entrada, "kmeros")?;
    let modo = entrada.opcion("copias", COPIAS)?;
    let mut lista = rebanadas(&kmeros);
    let leidos = lista.len();
    if modo == COPIAS[0] {
        let mut vistos = HashSet::new();
        lista.retain(|kmero| vistos.insert(*kmero));
    }
    let repetidos = leidos - lista.len();
    let real = ensamblaje::generar_contigs(&lista)?;
    let k = lista[0].len();
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let primeros: Vec<String> = lista
        .iter()
        .take(FILAS_DE_PIZARRA)
        .map(|kmero| letras(kmero))
        .collect();
    narracion.agregar(
        con_variables(&[
            ("kmeros.len()", miles(lista.len())),
            ("k", k.to_string()),
        ])
        .lista("k-meros (los primeros)", primeros, &[])
        .paso(
            "pub fn generar_contigs(kmeros: &[&[u8]]) -> Result<Vec<Vec<u8>>, String> {",
            format!(
                "Se leyeron {}{}. Un contig es un trozo del genoma que se puede armar sin dudas. La receta: armar el grafo de De Bruijn, cortarlo en los tramos sin ramificaciones y deletrear cada tramo.",
                cantidad(leidos, "k-mero", "k-meros"),
                if repetidos > 0 {
                    format!(
                        "; {} estaban repetidos y se juntaron, así que quedan {}",
                        miles(repetidos),
                        miles(lista.len())
                    )
                } else {
                    String::new()
                }
            ),
        ),
    );
    narracion.agregar(
        con_variables(&[("kmeros.len()", miles(lista.len()))]).paso(
            "let grafo = de_bruijn_de_kmeros(kmeros)?;",
            "Primero se arma el grafo de De Bruijn: cada k-mero es una flecha entre su prefijo y su sufijo.",
        ),
    );
    let grafo = narrar_nodos(&mut narracion, &lista, None, 28)?;
    if grafo != ensamblaje::de_bruijn_de_kmeros(&lista)? {
        return Err(DISTINTO_DEL_ALGORITMO.to_string());
    }
    let caminos = ensamblaje::caminos_no_ramificados(&grafo);
    narrar_caminos(&mut narracion, &grafo, &caminos, 24);
    if caminos.len() != real.len() {
        return Err(DISTINTO_DEL_ALGORITMO.to_string());
    }
    let total = caminos.len();
    for (numero, camino) in caminos.iter().enumerate() {
        if !se_muestra(numero, total, 20) {
            continue;
        }
        if real[numero].len() != camino.len() + k - 2 {
            return Err(DISTINTO_DEL_ALGORITMO.to_string());
        }
        let desde = (numero + 1).saturating_sub(FILAS_DE_PIZARRA);
        let lineas: Vec<String> = real[desde..=numero]
            .iter()
            .map(|contig| etiqueta_corta(contig))
            .collect();
        let resaltado = lineas.len() - 1;
        narracion.agregar(
            con_variables(&[
                ("camino.len()", miles(camino.len())),
                ("contigs.len()", miles(numero + 1)),
            ])
            .secuencia(&real[numero], &[(0, real[numero].len(), Resalte::Hallado)])
            .lista("contigs hasta ahora (los últimos)", lineas, &[resaltado])
            .paso(
                "contigs.push(deletrear_camino(&etiquetas_del_camino(&grafo, &camino))?);",
                format!(
                    "El camino de {} se deletrea como una cadena: el primer nodo completo y, de cada nodo siguiente, solo su última letra. Sale el contig «{}», de {} letras.",
                    cantidad(camino.len(), "nodo", "nodos"),
                    etiqueta_corta(&real[numero]),
                    miles(real[numero].len())
                ),
            ),
        );
    }
    let mas_largo = real.iter().map(Vec::len).max().unwrap_or(0);
    Ok(narracion.terminar(
        con_variables(&[
            ("contigs.len()", miles(real.len())),
            ("el más largo", format!("{} letras", miles(mas_largo))),
        ])
        .paso(
            "Ok(contigs)",
            format!(
                "Listo: {}, el más largo de {} letras. Cada contig es un tramo sin ramificaciones. En «Qué significa» se mide su calidad con el N50.",
                cantidad(real.len(), "contig", "contigs"),
                miles(mas_largo)
            ),
        ),
    ))
}

#[cfg(test)]
mod pruebas {
    use super::*;

    fn con_ejemplo(id: &str, cambios: &[(&str, &str)]) -> Entrada {
        let algoritmo = algoritmos()
            .into_iter()
            .find(|algoritmo| algoritmo.id == id)
            .unwrap();
        let mut pares: Vec<(&str, &str)> = algoritmo
            .campos
            .iter()
            .map(|campo| (campo.id, campo.ejemplo))
            .collect();
        for (campo, valor) in cambios {
            for par in pares.iter_mut() {
                if par.0 == *campo {
                    par.1 = valor;
                }
            }
        }
        Entrada::nueva(&pares)
    }

    fn texto_de(salida: &Salida, etiqueta_buscada: &str) -> String {
        for bloque in &salida.bloques {
            if let Bloque::Texto { etiqueta, texto } = bloque
                && etiqueta == etiqueta_buscada
            {
                return texto.clone();
            }
        }
        panic!("no hay bloque «{etiqueta_buscada}»");
    }

    fn del_gen(nombre: &str) -> String {
        letras(&ecoli::secuencia_gen(ecoli::buscar_gen(nombre).unwrap()))
    }

    #[test]
    fn el_ejemplo_de_reconstruccion_deletrea_thrl() {
        let salida = reconstruccion(&con_ejemplo("ensamblaje.reconstruccion", &[])).unwrap();
        assert_eq!(texto_de(&salida, "Cadena reconstruida"), del_gen("thrL"));
        assert!(salida
            .bloques
            .iter()
            .any(|bloque| matches!(bloque, Bloque::Texto { etiqueta, .. } if etiqueta == "Otra cadena con los mismos k-meros")));
        assert!(
            salida
                .interpretacion
                .iter()
                .any(|frase| frase.contains("thrL"))
        );
    }

    #[test]
    fn el_paso_a_paso_termina_igual_que_nucleo() {
        let entrada = con_ejemplo("ensamblaje.reconstruccion", &[]);
        let pasos = pasos_reconstruccion(&entrada).unwrap();
        assert!(pasos.len() <= MAXIMO_PASOS);
        let ultimo = pasos.last().unwrap();
        assert_eq!(ultimo.estado["texto"], Value::String(del_gen("thrL")));
        assert_eq!(ultimo.estado["fase"], "deletreo");
        let aristas = ultimo.estado["aristas"].as_array().unwrap();
        assert!(aristas.iter().all(|arista| arista["usada"] == true));
        let primero = &pasos[0];
        assert_eq!(primero.estado["aristas"].as_array().unwrap().len(), 1);
    }

    #[test]
    fn el_paso_a_paso_de_un_ciclo_y_sus_limites() {
        let entrada = Entrada::nueva(&[("kmeros", "CGT TAC ACG GTA")]);
        let pasos = pasos_reconstruccion(&entrada).unwrap();
        assert!(
            pasos
                .iter()
                .any(|paso| paso.linea.contains("unwrap_or_else"))
        );
        let muchos = vec!["ACGT"; MAXIMO_KMEROS_ANIMADOS + 1].join(" ");
        assert!(pasos_reconstruccion(&Entrada::nueva(&[("kmeros", &muchos)])).is_err());
    }

    #[test]
    fn muestrear_conserva_el_primero_y_el_ultimo() {
        let pasos: Vec<Paso> = (0..1000)
            .map(|numero| Paso::nuevo("x", numero.to_string(), Value::Null))
            .collect();
        let pocos = muestrear(pasos, 400);
        assert_eq!(pocos.len(), 400);
        assert_eq!(pocos[0].nota, "0");
        assert_eq!(pocos[399].nota, "999");
    }

    #[test]
    fn los_ejemplos_de_genes_se_ubican_en_ecoli() {
        let camino = camino_de_kmeros(&con_ejemplo("ensamblaje.camino_de_kmeros", &[])).unwrap();
        assert!(
            camino
                .interpretacion
                .iter()
                .any(|frase| frase.contains("lacZ"))
        );
        let huecos = camino_con_huecos(&con_ejemplo("ensamblaje.camino_con_huecos", &[])).unwrap();
        assert_eq!(texto_de(&huecos, "Cadena"), del_gen("trpL"));
        let pares = pares(&con_ejemplo("ensamblaje.pares", &[])).unwrap();
        assert_eq!(texto_de(&pares, "Cadena reconstruida"), del_gen("thrL"));
    }

    #[test]
    fn ciclo_sin_balance_se_explica_sin_error() {
        let entrada = Entrada::nueva(&[("grafo", "a -> b\na -> c")]);
        let salida = ciclo(&entrada).unwrap();
        assert!(
            salida
                .interpretacion
                .iter()
                .any(|frase| frase.contains("«a»"))
        );
        let ejemplo = ciclo(&con_ejemplo("ensamblaje.ciclo_euleriano", &[])).unwrap();
        assert_eq!(
            texto_de(&ejemplo, "Ciclo euleriano"),
            "0->1->2->3->1->5->6->7->6->3->4->0"
        );
        let camino = camino(&con_ejemplo("ensamblaje.camino_euleriano", &[])).unwrap();
        assert_eq!(
            texto_de(&camino, "Camino euleriano"),
            "1->2->3->4->2->6->7->8->7->9"
        );
    }

    #[test]
    fn n50_acepta_largos_o_secuencias() {
        let largos = n50(&Entrada::nueva(&[("contigs", "2 3 4 5 6 7 8 9 10")])).unwrap();
        assert!(matches!(&largos.bloques[0], Bloque::Valor { valor, .. } if valor == "8 letras"));
        let secuencias = n50(&Entrada::nueva(&[("contigs", "ACGT\nAC\nACGTACGT")])).unwrap();
        assert!(
            matches!(&secuencias.bloques[0], Bloque::Valor { valor, .. } if valor == "8 letras")
        );
        assert!(n50(&Entrada::nueva(&[("contigs", "0 0")])).is_err());
        assert!(n50(&Entrada::nueva(&[("contigs", "18446744073709551615 5")])).is_err());
    }

    fn dice(salida: &Salida, trozo: &str) -> bool {
        salida
            .interpretacion
            .iter()
            .any(|frase| frase.contains(trozo))
    }

    #[test]
    fn los_singulares_se_escriben_bien() {
        let composicion = composicion(&Entrada::nueva(&[
            ("texto", "@gen:thrL"),
            ("k", "12"),
            ("orden", "de aparición"),
        ]))
        .unwrap();
        assert!(dice(&composicion, "1 k-mero aparece más de una vez"));
        assert!(!dice(&composicion, "aparecen más de una vez; el"));
        let reconstruida = reconstruccion(&Entrada::nueva(&[("kmeros", "@gen:thrL")])).unwrap();
        assert!(!dice(&reconstruida, "los 1 k-meros"));
        let medida = n50(&Entrada::nueva(&[("contigs", "100")])).unwrap();
        assert!(!dice(&medida, "los un contig"));
        assert!(!dice(&medida, "miden (") && !dice(&medida, "mide ("));
        let uno = universal(&Entrada::nueva(&[("k", "2")])).unwrap();
        assert!(!dice(&uno, "= 1 cadenas"));
        let aislado = no_ramificados(&Entrada::nueva(&[(
            "grafo",
            "a -> b\nb -> c\nx -> y\ny -> x",
        )]))
        .unwrap();
        assert!(dice(&aislado, "Uno de ellos es un ciclo aislado"));
    }

    #[test]
    fn no_dice_una_sola_fila_si_hay_trozos_separados() {
        let separados =
            de_bruijn_kmeros(&Entrada::nueva(&[("kmeros", "ACGT CGTA TTTT GGGG")])).unwrap();
        assert!(!dice(&separados, "única forma de recorrer"));
        assert!(dice(&separados, "trozos separados"));
        let con_ciclo =
            solapamiento(&Entrada::nueva(&[("kmeros", "AAC ACG CTT TTC TCT")])).unwrap();
        assert!(!dice(&con_ciclo, "una sola fila"));
        let fila = solapamiento(&Entrada::nueva(&[("kmeros", "AAC ACG CGT")])).unwrap();
        assert!(dice(&fila, "una sola fila"));
    }

    #[test]
    fn contigs_sin_ramificaciones_hablan_de_cobertura() {
        let salida = contigs(&Entrada::nueva(&[
            ("kmeros", "ACGT CGTA TTTT GGGA"),
            ("copias", "juntar las iguales"),
        ]))
        .unwrap();
        assert!(dice(&salida, "cobertura"));
        assert!(!dice(&salida, "se ramifican"));
        assert!(!dice(&salida, "un k mayor"));
    }

    #[test]
    fn acepta_kmeros_binarios() {
        let salida = reconstruccion(&Entrada::nueva(&[(
            "kmeros",
            "0001 0011 0110 1100 1000 0000",
        )]))
        .unwrap();
        assert_eq!(texto_de(&salida, "Cadena reconstruida").len(), 9);
        assert!(reconstruccion(&Entrada::nueva(&[("kmeros", "AC-T CTGA")])).is_err());
    }

    #[test]
    fn camino_con_huecos_sin_superposicion() {
        let salida = camino_con_huecos(&Entrada::nueva(&[
            ("pares", "ACG|GCA CGT|CAA GTT|AAG"),
            ("d", "2"),
        ]))
        .unwrap();
        assert_eq!(texto_de(&salida, "Cadena"), "ACGTTGCAAG");
        assert!(dice(&salida, "no se superponen"));
    }

    #[test]
    fn el_paso_a_paso_explica_los_bucles() {
        let pasos = pasos_reconstruccion(&Entrada::nueva(&[("kmeros", "AAAA AAAC")])).unwrap();
        assert!(pasos[0].nota.contains("bucle"), "{}", pasos[0].nota);
    }

    #[test]
    fn universal_verifica_cada_kmero() {
        let salida = universal(&Entrada::nueva(&[("k", "3")])).unwrap();
        assert_eq!(texto_de(&salida, "Cadena circular"), "00010111");
        assert!(universal(&Entrada::nueva(&[("k", "21")])).is_err());
        assert!(universal(&Entrada::nueva(&[("k", "1")])).is_ok());
    }

    fn trazar(id: &str, cambios: &[(&str, &str)]) -> Vec<Paso> {
        let algoritmo = algoritmos()
            .into_iter()
            .find(|algoritmo| algoritmo.id == id)
            .unwrap();
        let trazador = algoritmo.pasos.unwrap();
        let pasos = trazador(&con_ejemplo(id, cambios)).unwrap();
        assert!(!pasos.is_empty() && pasos.len() <= MAXIMO_PASOS);
        let texto: String = crate::codigo(id)
            .unwrap()
            .iter()
            .map(|fragmento| fragmento.fuente.as_str())
            .collect::<Vec<_>>()
            .join("\n");
        for paso in &pasos {
            assert!(
                texto
                    .lines()
                    .any(|linea| linea.contains(paso.linea.as_str())),
                "{id}: la línea «{}» no está en el código",
                paso.linea
            );
            assert!(!paso.nota.is_empty());
        }
        let tamano = serde_json::to_string(&pasos).unwrap().len();
        assert!(tamano < 2_000_000, "{id}: {tamano} bytes");
        pasos
    }

    fn valor_final(pasos: &[Paso], nombre: &str) -> String {
        let ultimo = pasos.last().unwrap();
        ultimo.estado["variables"]
            .as_array()
            .unwrap()
            .iter()
            .find(|par| par[0] == nombre)
            .unwrap_or_else(|| panic!("no hay la variable «{nombre}»"))[1]
            .as_str()
            .unwrap()
            .to_string()
    }

    fn kmeros_de(entrada: &Entrada) -> Vec<Vec<u8>> {
        leer_kmeros(entrada, "kmeros").unwrap()
    }

    #[test]
    fn la_narracion_de_composicion_da_el_mismo_resultado() {
        let entrada = con_ejemplo("ensamblaje.composicion", &[]);
        let texto = entrada.secuencia("texto").unwrap();
        let pasos = trazar("ensamblaje.composicion", &[]);
        let esperado = ensamblaje::composicion(&texto, 5);
        assert_eq!(
            valor_final(&pasos, "resultado.len()"),
            miles(esperado.len())
        );
        assert!(pasos.len() > 10);
    }

    #[test]
    fn la_narracion_de_composicion_resume_un_genoma_entero() {
        let pasos = trazar(
            "ensamblaje.composicion",
            &[("texto", "@ecoli"), ("k", "12")],
        );
        assert!(pasos.len() > 50);
        assert_eq!(valor_final(&pasos, "resultado.len()"), miles(4_641_641));
    }

    #[test]
    fn la_narracion_del_camino_de_kmeros_da_el_mismo_texto() {
        let entrada = con_ejemplo("ensamblaje.camino_de_kmeros", &[]);
        let kmeros = kmeros_de(&entrada);
        let esperado = ensamblaje::deletrear_camino(&rebanadas(&kmeros)).unwrap();
        let pasos = trazar("ensamblaje.camino_de_kmeros", &[]);
        assert_eq!(valor_final(&pasos, "texto.len()"), miles(esperado.len()));
    }

    #[test]
    fn la_narracion_de_un_camino_roto_da_el_mismo_error() {
        let entrada = Entrada::nueva(&[("kmeros", "ACG GTA")]);
        assert!(camino_de_kmeros(&entrada).is_err());
        assert!(pasos_camino_de_kmeros(&entrada).is_err());
    }

    #[test]
    fn la_narracion_del_solapamiento_da_las_mismas_aristas() {
        let entrada = con_ejemplo("ensamblaje.solapamiento", &[]);
        let kmeros = kmeros_de(&entrada);
        let grafo = ensamblaje::grafo_de_solapamiento(&rebanadas(&kmeros)).unwrap();
        let pasos = trazar("ensamblaje.solapamiento", &[]);
        assert_eq!(
            valor_final(&pasos, "aristas"),
            miles(ensamblaje::cantidad_de_aristas(&grafo))
        );
    }

    #[test]
    fn la_narracion_de_de_bruijn_da_el_mismo_grafo() {
        let entrada = con_ejemplo("ensamblaje.de_bruijn_texto", &[]);
        let texto = entrada.secuencia("texto").unwrap();
        let grafo = ensamblaje::de_bruijn_de_texto(&texto, 4).unwrap();
        let pasos = trazar("ensamblaje.de_bruijn_texto", &[]);
        assert_eq!(
            valor_final(&pasos, "grafo.nodos.len()"),
            miles(grafo.nodos.len())
        );
        assert_eq!(
            valor_final(&pasos, "aristas"),
            miles(ensamblaje::cantidad_de_aristas(&grafo))
        );
        let entrada = con_ejemplo("ensamblaje.de_bruijn_kmeros", &[]);
        let kmeros = kmeros_de(&entrada);
        let grafo = ensamblaje::de_bruijn_de_kmeros(&rebanadas(&kmeros)).unwrap();
        let pasos = trazar("ensamblaje.de_bruijn_kmeros", &[]);
        assert_eq!(
            valor_final(&pasos, "grafo.nodos.len()"),
            miles(grafo.nodos.len())
        );
        assert_eq!(
            valor_final(&pasos, "aristas"),
            miles(ensamblaje::cantidad_de_aristas(&grafo))
        );
    }

    #[test]
    fn la_narracion_del_ciclo_da_el_mismo_recorrido() {
        let entrada = con_ejemplo("ensamblaje.ciclo_euleriano", &[]);
        let grafo = leer_grafo(&entrada, "grafo").unwrap();
        let ciclo = ensamblaje::ciclo_euleriano(&grafo).unwrap();
        let pasos = trazar("ensamblaje.ciclo_euleriano", &[]);
        assert_eq!(valor_final(&pasos, "ciclo.len()"), miles(ciclo.len()));
    }

    #[test]
    fn la_narracion_del_ciclo_explica_por_que_no_hay() {
        let pasos = trazar(
            "ensamblaje.ciclo_euleriano",
            &[("grafo", "1 -> 2\n2 -> 3\n3 -> 2")],
        );
        assert!(
            pasos
                .last()
                .unwrap()
                .nota
                .contains("No hay ciclo euleriano")
        );
        let pasos = trazar(
            "ensamblaje.ciclo_euleriano",
            &[("grafo", "1 -> 2\n2 -> 1\n3 -> 4\n4 -> 3")],
        );
        assert!(pasos.last().unwrap().nota.contains("grupos separados"));
    }

    #[test]
    fn la_narracion_del_camino_euleriano_da_el_mismo_recorrido() {
        let entrada = con_ejemplo("ensamblaje.camino_euleriano", &[]);
        let grafo = leer_grafo(&entrada, "grafo").unwrap();
        let camino = ensamblaje::camino_euleriano(&grafo).unwrap();
        let pasos = trazar("ensamblaje.camino_euleriano", &[]);
        assert_eq!(valor_final(&pasos, "camino.len()"), miles(camino.len()));
        let pasos = trazar(
            "ensamblaje.camino_euleriano",
            &[("grafo", "1 -> 2,3\n2 -> 4\n3 -> 4")],
        );
        assert!(
            pasos
                .last()
                .unwrap()
                .nota
                .contains("No hay camino euleriano")
        );
        let pasos = trazar(
            "ensamblaje.camino_euleriano",
            &[("grafo", "1 -> 2\n2 -> 1")],
        );
        assert!(pasos.last().unwrap().nota.contains("Listo"));
    }

    #[test]
    fn la_narracion_de_la_cadena_universal_da_la_misma_cadena() {
        let pasos = trazar("ensamblaje.universal", &[]);
        let esperado = ensamblaje::cadena_universal_circular(4).unwrap();
        assert_eq!(valor_final(&pasos, "cadena"), letras(&esperado));
        for k in ["1", "2", "3", "10", "20"] {
            let pasos = trazar("ensamblaje.universal", &[("k", k)]);
            let esperado = ensamblaje::cadena_universal_circular(k.parse().unwrap()).unwrap();
            if k == "1" {
                assert_eq!(valor_final(&pasos, "cadena"), "01");
            } else {
                assert_eq!(valor_final(&pasos, "cadena.len()"), miles(esperado.len()));
            }
        }
    }

    #[test]
    fn la_narracion_de_pares_da_el_mismo_genoma() {
        let entrada = con_ejemplo("ensamblaje.pares", &[]);
        let lista = leer_pares(&entrada, "pares").unwrap();
        let esperado = ensamblaje::reconstruir_desde_pares(&lista, 3, 1).unwrap();
        let pasos = trazar("ensamblaje.pares", &[]);
        assert_eq!(
            valor_final(&pasos, "genoma.len()"),
            miles(esperado.genoma.len())
        );
        assert_eq!(
            valor_final(&pasos, "intentos"),
            esperado.intentos.to_string()
        );
        for semilla in ["0", "2", "7"] {
            let pasos = trazar("ensamblaje.pares", &[("semilla", semilla)]);
            let esperado =
                ensamblaje::reconstruir_desde_pares(&lista, 3, semilla.parse().unwrap()).unwrap();
            assert_eq!(
                valor_final(&pasos, "intentos"),
                esperado.intentos.to_string()
            );
        }
    }

    #[test]
    fn la_narracion_de_contigs_da_los_mismos_contigs() {
        let entrada = con_ejemplo("ensamblaje.contigs", &[]);
        let kmeros = kmeros_de(&entrada);
        let mut lista = rebanadas(&kmeros);
        let mut vistos = HashSet::new();
        lista.retain(|kmero| vistos.insert(*kmero));
        let esperado = ensamblaje::generar_contigs(&lista).unwrap();
        let pasos = trazar("ensamblaje.contigs", &[]);
        assert_eq!(valor_final(&pasos, "contigs.len()"), miles(esperado.len()));
        let pasos = trazar("ensamblaje.contigs", &[("copias", COPIAS[1])]);
        assert!(pasos.len() > 10);
    }

    #[test]
    fn la_narracion_con_huecos_da_la_misma_cadena() {
        let entrada = con_ejemplo("ensamblaje.camino_con_huecos", &[]);
        let pares = leer_pares(&entrada, "pares").unwrap();
        let esperado = ensamblaje::deletrear_camino_con_huecos(&pares, 3).unwrap();
        let pasos = trazar("ensamblaje.camino_con_huecos", &[]);
        assert_eq!(valor_final(&pasos, "genoma.len()"), miles(esperado.len()));
        let entrada = Entrada::nueva(&[("pares", "AAA|CCC CCC|GGG"), ("d", "3")]);
        assert!(camino_con_huecos(&entrada).is_err());
        assert!(pasos_camino_con_huecos(&entrada).is_err());
    }

    #[test]
    fn la_narracion_de_caminos_no_ramificados_da_los_mismos_caminos() {
        let entrada = con_ejemplo("ensamblaje.no_ramificados", &[]);
        let grafo = leer_grafo(&entrada, "grafo").unwrap();
        let esperado = ensamblaje::caminos_no_ramificados(&grafo);
        let pasos = trazar("ensamblaje.no_ramificados", &[]);
        assert_eq!(valor_final(&pasos, "caminos.len()"), miles(esperado.len()));
    }

    #[test]
    fn la_narracion_de_n50_da_la_misma_medida() {
        let entrada = con_ejemplo("ensamblaje.n50", &[]);
        let longitudes = leer_longitudes(&entrada, "contigs").unwrap();
        let esperado = ensamblaje::medir_n50(&longitudes);
        let pasos = trazar("ensamblaje.n50", &[]);
        assert_eq!(valor_final(&pasos, "n50"), miles(esperado.n50));
        assert_eq!(valor_final(&pasos, "l50"), miles(esperado.l50));
        assert_eq!(valor_final(&pasos, "total"), miles(esperado.total));
        let muchas: Vec<String> = (1..=5000).map(|largo| largo.to_string()).collect();
        let pasos = trazar("ensamblaje.n50", &[("contigs", &muchas.join(" "))]);
        assert!(pasos.len() <= MAXIMO_PASOS);
    }
}
