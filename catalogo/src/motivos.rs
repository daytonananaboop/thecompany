use crate::formato::{cientifico, decimal, miles, porcentaje};
use crate::pizarra::{Marca as MarcaDeSecuencia, Narracion, Pizarra};
use crate::{Algoritmo, Bloque, Campo, Entrada, FilaAlineada, Marca, Paso, Salida, Tipo};
use nucleo::azar::Azar;
use nucleo::motivos::{self, AnalisisMotivos, Mediana, MotivosEncontrados};
use nucleo::{ecoli, secuencia};
use serde_json::{Value, json};
use std::collections::BTreeMap;

const REGIONES_SOS: &str = "@ecoli:4257015-4257114
@ecoli:2823770-2823869
@ecoli:1020920-1021019
@ecoli:1230667-1230766
@ecoli:2751695-2751794
@ecoli:3997883-3997982
@ecoli:3817660-3817759
@ecoli:2081263-2081362";

const CAJAS_SOS: &str = "@ecoli:4257070-4257085
@ecoli:2823831-2823846
@ecoli:1020944-1020959
@ecoli:1230730-1230745
@ecoli:2751729-2751744
@ecoli:3997909-3997924
@ecoli:3817699-3817714
@ecoli:2081279-2081294";

const REGIONES_CON_CONTROLES: &str = "@ecoli:4257015-4257114
@ecoli:2823770-2823869
@ecoli:1020920-1021019
@ecoli:1230667-1230766
@ecoli:2751695-2751794
@ecoli:813426-813525
@ecoli:250798-250897
@ecoli:3212947-3213046
@ecoli:1292923-1293022";

const PERFIL_SOS: &str = "0.083 0.083 0.083 0.083 0.667 0.083 0.500 0.083 0.667 0.333 0.583 0.250 0.583 0.083 0.750 0.083
0.750 0.083 0.083 0.083 0.083 0.083 0.083 0.167 0.083 0.250 0.083 0.250 0.250 0.750 0.083 0.083
0.083 0.083 0.750 0.083 0.083 0.083 0.250 0.167 0.083 0.083 0.167 0.083 0.083 0.083 0.083 0.750
0.083 0.750 0.083 0.750 0.167 0.750 0.167 0.583 0.167 0.333 0.167 0.417 0.083 0.083 0.083 0.083";

const CAJA_LEXA: &str = "CTGTATATATATACAG";
const CAJA_LEXA_EXTENDIDA: &[u8] = b"TACTGTATATATATACAGTA";
const LETRAS: [&str; 4] = ["A", "C", "G", "T"];

const K_MAXIMA_EXHAUSTIVA: usize = 10;
const D_MAXIMA: usize = 3;
const K_MAXIMA_PERFIL: usize = 50;
const LARGO_MAXIMO_DE_MOTIVO: usize = 100;
const MAXIMO_DE_SECUENCIAS: usize = 1000;
const PRESUPUESTO: f64 = 1.5e9;
const PRESUPUESTO_ENUMERACION: f64 = 3e8;
const PRESUPUESTO_LINEAL: f64 = 2e9;
const ITERACIONES_TIPICAS: f64 = 10.0;
const MOTIVOS_EN_TABLA: usize = 20;
const DISTANCIA_A_UN_GEN: usize = 300;
const AZAR_TOLERADO: f64 = 0.01;
const MAXIMO_PASOS: usize = 400;
const SECUENCIAS_ANIMABLES: usize = 20;
const BASES_ANIMABLES: usize = 3000;

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![
        Algoritmo {
            id: "motivos.enumeracion",
            capitulo: "motivos",
            nombre: "Enumeración de motivos (k, d)",
            resumen: "Todos los k-mers que aparecen en cada secuencia con como mucho d diferencias.",
            campos: &[
                Campo {
                    id: "adn",
                    etiqueta: "Secuencias de ADN",
                    tipo: Tipo::Secuencias,
                    ejemplo: REGIONES_SOS,
                },
                Campo {
                    id: "k",
                    etiqueta: "Longitud del motivo (k)",
                    tipo: Tipo::Entero,
                    ejemplo: "8",
                },
                Campo {
                    id: "d",
                    etiqueta: "Diferencias permitidas (d)",
                    tipo: Tipo::Entero,
                    ejemplo: "1",
                },
            ],
            funcion: "motivos::enumerar_motivos",
            ejecutar: enumeracion,
            pasos: Some(pasos_enumeracion),
        },
        Algoritmo {
            id: "motivos.perfil",
            capitulo: "motivos",
            nombre: "Perfil, consenso, puntuación y entropía",
            resumen: "Cuenta las bases de cada columna de una matriz de motivos y mide cuán conservada está.",
            campos: &[Campo {
                id: "motivos",
                etiqueta: "Motivos (todos de la misma longitud)",
                tipo: Tipo::Secuencias,
                ejemplo: CAJAS_SOS,
            }],
            funcion: "motivos::analizar_motivos",
            ejecutar: perfil,
            pasos: Some(pasos_perfil),
        },
        Algoritmo {
            id: "motivos.distancia",
            capitulo: "motivos",
            nombre: "Distancia entre un patrón y varias secuencias",
            resumen: "Suma, secuencia por secuencia, las diferencias con el k-mer más parecido al patrón.",
            campos: &[
                Campo {
                    id: "patron",
                    etiqueta: "Patrón",
                    tipo: Tipo::Secuencia,
                    ejemplo: CAJA_LEXA,
                },
                Campo {
                    id: "adn",
                    etiqueta: "Secuencias de ADN",
                    tipo: Tipo::Secuencias,
                    ejemplo: REGIONES_CON_CONTROLES,
                },
            ],
            funcion: "motivos::distancia_patron_textos",
            ejecutar: distancia,
            pasos: Some(pasos_distancia),
        },
        Algoritmo {
            id: "motivos.mediana",
            capitulo: "motivos",
            nombre: "Cadena mediana",
            resumen: "El k-mer que minimiza la suma de distancias a todas las secuencias, probando los 4^k posibles.",
            campos: &[
                Campo {
                    id: "adn",
                    etiqueta: "Secuencias de ADN",
                    tipo: Tipo::Secuencias,
                    ejemplo: REGIONES_SOS,
                },
                Campo {
                    id: "k",
                    etiqueta: "Longitud del patrón (k)",
                    tipo: Tipo::Entero,
                    ejemplo: "8",
                },
            ],
            funcion: "motivos::cadena_mediana",
            ejecutar: mediana,
            pasos: Some(pasos_mediana),
        },
        Algoritmo {
            id: "motivos.mas_probable",
            capitulo: "motivos",
            nombre: "k-mer más probable según un perfil",
            resumen: "Recorre el texto y elige la ventana con mayor probabilidad según la matriz de perfil.",
            campos: &[
                Campo {
                    id: "texto",
                    etiqueta: "Texto de ADN",
                    tipo: Tipo::Secuencia,
                    ejemplo: "@ecoli:813426-813525",
                },
                Campo {
                    id: "perfil",
                    etiqueta: "Perfil (4 filas A, C, G, T; una columna por posición)",
                    tipo: Tipo::Matriz,
                    ejemplo: PERFIL_SOS,
                },
            ],
            funcion: "motivos::posicion_mas_probable",
            ejecutar: mas_probable,
            pasos: Some(pasos_mas_probable),
        },
        Algoritmo {
            id: "motivos.codiciosa",
            capitulo: "motivos",
            nombre: "Búsqueda codiciosa de motivos",
            resumen: "Prueba cada k-mer de la primera secuencia y completa los motivos con el k-mer más probable de cada fila.",
            campos: &[
                Campo {
                    id: "adn",
                    etiqueta: "Secuencias de ADN",
                    tipo: Tipo::Secuencias,
                    ejemplo: REGIONES_SOS,
                },
                Campo {
                    id: "k",
                    etiqueta: "Longitud del motivo (k)",
                    tipo: Tipo::Entero,
                    ejemplo: "16",
                },
            ],
            funcion: "motivos::busqueda_codiciosa",
            ejecutar: codiciosa,
            pasos: Some(pasos_codiciosa),
        },
        Algoritmo {
            id: "motivos.codiciosa_laplace",
            capitulo: "motivos",
            nombre: "Búsqueda codiciosa con pseudoconteos",
            resumen: "La misma búsqueda codiciosa, sumando 1 a cada conteo (regla de Laplace) para que ninguna probabilidad sea cero.",
            campos: &[
                Campo {
                    id: "adn",
                    etiqueta: "Secuencias de ADN",
                    tipo: Tipo::Secuencias,
                    ejemplo: REGIONES_SOS,
                },
                Campo {
                    id: "k",
                    etiqueta: "Longitud del motivo (k)",
                    tipo: Tipo::Entero,
                    ejemplo: "16",
                },
            ],
            funcion: "motivos::busqueda_codiciosa",
            ejecutar: codiciosa_laplace,
            pasos: Some(pasos_codiciosa_laplace),
        },
        Algoritmo {
            id: "motivos.aleatoria",
            capitulo: "motivos",
            nombre: "Búsqueda aleatoria de motivos",
            resumen: "Empieza en k-mers al azar y mejora con el perfil hasta que deja de mejorar; repite y se queda con lo mejor.",
            campos: &[
                Campo {
                    id: "adn",
                    etiqueta: "Secuencias de ADN",
                    tipo: Tipo::Secuencias,
                    ejemplo: REGIONES_SOS,
                },
                Campo {
                    id: "k",
                    etiqueta: "Longitud del motivo (k)",
                    tipo: Tipo::Entero,
                    ejemplo: "16",
                },
                Campo {
                    id: "repeticiones",
                    etiqueta: "Repeticiones",
                    tipo: Tipo::Entero,
                    ejemplo: "200",
                },
                Campo {
                    id: "semilla",
                    etiqueta: "Semilla del azar",
                    tipo: Tipo::Entero,
                    ejemplo: "7",
                },
            ],
            funcion: "motivos::busqueda_aleatoria",
            ejecutar: aleatoria,
            pasos: Some(pasos_aleatoria),
        },
        Algoritmo {
            id: "motivos.gibbs",
            capitulo: "motivos",
            nombre: "Muestreo de Gibbs",
            resumen: "Cambia un motivo por vez, sorteando su reemplazo según el perfil de los demás.",
            campos: &[
                Campo {
                    id: "adn",
                    etiqueta: "Secuencias de ADN",
                    tipo: Tipo::Secuencias,
                    ejemplo: REGIONES_SOS,
                },
                Campo {
                    id: "k",
                    etiqueta: "Longitud del motivo (k)",
                    tipo: Tipo::Entero,
                    ejemplo: "16",
                },
                Campo {
                    id: "iteraciones",
                    etiqueta: "Iteraciones por reinicio (N)",
                    tipo: Tipo::Entero,
                    ejemplo: "500",
                },
                Campo {
                    id: "reinicios",
                    etiqueta: "Reinicios",
                    tipo: Tipo::Entero,
                    ejemplo: "20",
                },
                Campo {
                    id: "semilla",
                    etiqueta: "Semilla del azar",
                    tipo: Tipo::Entero,
                    ejemplo: "7",
                },
            ],
            funcion: "motivos::muestreo_gibbs",
            ejecutar: gibbs,
            pasos: Some(pasos_gibbs),
        },
    ]
}

struct Conjunto {
    textos: Vec<Vec<u8>>,
    nombres: Vec<String>,
    inicios: Vec<Option<usize>>,
}

fn texto(secuencia: &[u8]) -> String {
    String::from_utf8_lossy(secuencia).into_owned()
}

fn revisar_adn(nombre: &str, secuencia: &[u8]) -> Result<(), String> {
    if let Some(letra) = secuencia
        .iter()
        .find(|letra| !matches!(letra, b'A' | b'C' | b'G' | b'T'))
    {
        return Err(format!(
            "La secuencia «{nombre}» solo puede tener A, C, G y T, pero contiene «{}».",
            *letra as char
        ));
    }
    Ok(())
}

fn gen_cercano(inicio: usize, fin: usize) -> Option<&'static str> {
    let mut elegido: Option<(usize, &'static str)> = None;
    for candidato in ecoli::genes() {
        let distancia = if candidato.hebra == b'+' {
            if candidato.inicio < inicio {
                continue;
            }
            candidato.inicio.saturating_sub(fin)
        } else {
            if candidato.fin > fin {
                continue;
            }
            inicio.saturating_sub(candidato.fin)
        };
        if distancia <= DISTANCIA_A_UN_GEN && elegido.is_none_or(|(mejor, _)| distancia < mejor) {
            elegido = Some((distancia, candidato.nombre.as_str()));
        }
    }
    elegido.map(|(_, nombre)| nombre)
}

fn describir_origen(pieza: &str, por_defecto: String) -> (String, Option<usize>) {
    if pieza == "@ecoli" {
        return ("genoma de E. coli".to_string(), Some(1));
    }
    if let Some(rango) = pieza.strip_prefix("@ecoli:")
        && let Some((inicio, fin)) = rango.split_once('-')
        && let (Ok(inicio), Ok(fin)) = (inicio.trim().parse::<usize>(), fin.trim().parse::<usize>())
    {
        let nombre = match gen_cercano(inicio, fin) {
            Some(gen_vecino) => format!("{gen_vecino} ({inicio}-{fin})"),
            None => format!("E. coli {inicio}-{fin}"),
        };
        return (nombre, Some(inicio));
    }
    if let Some(nombre) = pieza.strip_prefix("@gen:") {
        return (nombre.trim().to_string(), None);
    }
    (por_defecto, None)
}

fn leer_conjunto(entrada: &Entrada, campo: &str, minimo: usize) -> Result<Conjunto, String> {
    let registros = entrada.registros(campo)?;
    if registros.len() < minimo {
        return Err(format!(
            "El campo «{campo}» necesita al menos {minimo} secuencias y tiene {}.",
            registros.len()
        ));
    }
    if registros.len() > MAXIMO_DE_SECUENCIAS {
        return Err(format!(
            "El campo «{campo}» tiene {} secuencias; el máximo es {MAXIMO_DE_SECUENCIAS}.",
            miles(registros.len())
        ));
    }
    for registro in &registros {
        revisar_adn(&registro.nombre, &registro.secuencia)?;
    }
    let crudo = entrada.texto(campo)?;
    let es_fasta = crudo.trim_start().starts_with('>');
    let piezas: Vec<&str> = crudo
        .split(|c: char| c.is_whitespace() || c == ',' || c == ';')
        .map(str::trim)
        .filter(|pieza| !pieza.is_empty())
        .collect();
    let mut conjunto = Conjunto {
        textos: Vec::new(),
        nombres: Vec::new(),
        inicios: Vec::new(),
    };
    for (indice, registro) in registros.into_iter().enumerate() {
        let pieza = if es_fasta {
            ""
        } else {
            piezas.get(indice).copied().unwrap_or("")
        };
        let por_defecto = if es_fasta {
            registro.nombre.clone()
        } else {
            format!("s{}", indice + 1)
        };
        let (nombre, inicio) = describir_origen(pieza, por_defecto);
        conjunto.textos.push(registro.secuencia);
        conjunto.nombres.push(nombre);
        conjunto.inicios.push(inicio);
    }
    Ok(conjunto)
}

fn leer_texto(entrada: &Entrada, campo: &str) -> Result<(Vec<u8>, Option<usize>), String> {
    let adn = entrada.adn(campo)?;
    let crudo = entrada.texto(campo)?;
    let (_, inicio) = describir_origen(crudo.trim(), String::new());
    Ok((adn, inicio))
}

fn leer_k(entrada: &Entrada, textos: &[Vec<u8>]) -> Result<usize, String> {
    let k = entrada.natural("k")?;
    if k == 0 {
        return Err("k debe ser al menos 1.".to_string());
    }
    let mas_corta = textos.iter().map(Vec::len).min().unwrap_or(0);
    if k > mas_corta {
        return Err(format!(
            "k = {k} es mayor que la secuencia más corta ({} pb).",
            miles(mas_corta)
        ));
    }
    Ok(k)
}

fn leer_k_de_perfil(entrada: &Entrada, textos: &[Vec<u8>]) -> Result<usize, String> {
    let k = leer_k(entrada, textos)?;
    if k > K_MAXIMA_PERFIL {
        return Err(format!(
            "k = {k} es demasiado largo: los motivos regulatorios suelen medir entre 6 y 30 pb y aquí el máximo es {K_MAXIMA_PERFIL}."
        ));
    }
    Ok(k)
}

fn leer_semilla(entrada: &Entrada) -> Result<u64, String> {
    let texto = entrada.texto("semilla")?;
    texto.trim().parse::<u64>().map_err(|_| {
        "La semilla debe ser un número entero entre 0 y 18 446 744 073 709 551 615.".to_string()
    })
}

fn revisar_costo(costo: f64, limite: f64) -> Result<(), String> {
    if costo > limite {
        return Err(format!(
            "Esta entrada pediría unas {} operaciones y el límite es {} para que no tarde más de unos segundos. Usa menos secuencias, secuencias más cortas o un k menor.",
            cientifico(costo),
            cientifico(limite)
        ));
    }
    Ok(())
}

fn ventanas(longitud: usize, k: usize) -> f64 {
    (longitud + 1).saturating_sub(k) as f64
}

fn total_de_ventanas(textos: &[Vec<u8>], k: usize) -> f64 {
    textos.iter().map(|texto| ventanas(texto.len(), k)).sum()
}

fn combinaciones(n: usize, r: usize) -> f64 {
    let mut resultado = 1.0;
    for i in 0..r.min(n) {
        resultado = resultado * (n - i) as f64 / (i + 1) as f64;
    }
    resultado
}

fn binomial_hasta(n: usize, d: usize, p: f64) -> f64 {
    let mut suma = 0.0;
    for i in 0..=d.min(n) {
        suma += combinaciones(n, i) * p.powi(i as i32) * (1.0 - p).powi((n - i) as i32);
    }
    suma.min(1.0)
}

fn probabilidad_en_ventana(k: usize, d: usize) -> f64 {
    binomial_hasta(k, d, 0.75)
}

fn probabilidad_en_texto(k: usize, d: usize, longitud: usize) -> f64 {
    let p = probabilidad_en_ventana(k, d);
    if p >= 1.0 {
        return 1.0;
    }
    1.0 - (ventanas(longitud, k) * (-p).ln_1p()).exp()
}

fn cola_de_poisson(media: f64, minimo: usize) -> f64 {
    if minimo == 0 || minimo as f64 <= media {
        return 1.0;
    }
    if minimo as f64 > media + 12.0 * media.sqrt() + 50.0 {
        return 0.0;
    }
    let mut logaritmo_del_termino = -media;
    let mut acumulado = 0.0;
    for i in 0..minimo {
        acumulado += logaritmo_del_termino.exp();
        logaritmo_del_termino += media.ln() - ((i + 1) as f64).ln();
    }
    (1.0 - acumulado).max(0.0)
}

fn tamano_del_vecindario(k: usize, d: usize) -> f64 {
    let mut total = 0.0;
    for i in 0..=d.min(k) {
        total += combinaciones(k, i) * 3f64.powi(i as i32);
    }
    total
}

fn peor_puntuacion(t: usize, k: usize) -> usize {
    k * (t - t.div_ceil(4))
}

fn informacion(entropia: &[f64]) -> Vec<f64> {
    entropia.iter().map(|bits| (2.0 - bits).max(0.0)).collect()
}

fn diferencias_toleradas(k: usize, comparaciones: usize) -> usize {
    let mut toleradas = 0;
    for d in 0..=k / 4 {
        if comparaciones as f64 * probabilidad_en_ventana(k, d) <= AZAR_TOLERADO {
            toleradas = d;
        }
    }
    toleradas
}

fn parecido_a_lexa(consenso: &[u8]) -> Option<usize> {
    let k = consenso.len();
    if k < 8 {
        return None;
    }
    let (corta, larga) = if k <= CAJA_LEXA_EXTENDIDA.len() {
        (consenso, CAJA_LEXA_EXTENDIDA)
    } else {
        (CAJA_LEXA_EXTENDIDA, consenso)
    };
    let mut menor = corta.len();
    for ventana in larga.windows(corta.len()) {
        menor = menor.min(motivos::distancia_hamming(corta, ventana));
    }
    let comparaciones = larga.len() - corta.len() + 1;
    if menor <= diferencias_toleradas(corta.len(), comparaciones) {
        Some(menor)
    } else {
        None
    }
}

fn pares_no_complementarios(consenso: &[u8]) -> usize {
    let inverso = secuencia::complemento_inverso(consenso);
    motivos::distancia_hamming(consenso, &inverso) / 2
}

fn es_casi_palindromo(consenso: &[u8]) -> bool {
    let pares = consenso.len() / 2;
    if pares < 4 {
        return false;
    }
    let mut tolerados = 0;
    for m in 0..=pares {
        if binomial_hasta(pares, m, 0.75) <= AZAR_TOLERADO {
            tolerados = m;
        }
    }
    binomial_hasta(pares, 0, 0.75) <= AZAR_TOLERADO
        && pares_no_complementarios(consenso) <= tolerados
}

fn aproximado(valor: f64) -> String {
    if valor == 0.0 {
        return "0".to_string();
    }
    if !(0.01..10_000.0).contains(&valor) {
        return cientifico(valor);
    }
    if valor >= 10.0 {
        return miles(valor.round() as usize);
    }
    let cifras = if valor >= 0.1 { 1 } else { 2 };
    decimal(valor, cifras)
}

fn contar(cantidad: usize, singular: &str, plural: &str) -> String {
    let palabra = if cantidad == 1 { singular } else { plural };
    format!("{} {palabra}", miles(cantidad))
}

fn explicar_lexa_y_palindromo(salida: &mut Salida, consenso: &[u8], quien: &str) {
    let k = consenso.len();
    let palindromo = es_casi_palindromo(consenso);
    match parecido_a_lexa(consenso) {
        Some(diferencias) => {
            let caja = if k < CAJA_LEXA.len() {
                format!(
                    "un trozo de la caja SOS de E. coli ({}, cuyo núcleo de 16 pb es {CAJA_LEXA})",
                    texto(CAJA_LEXA_EXTENDIDA)
                )
            } else {
                format!("la caja SOS de E. coli (16 pb, consenso {CAJA_LEXA})")
            };
            let comparadas = k.min(CAJA_LEXA_EXTENDIDA.len());
            let mut frase = match (diferencias, k > CAJA_LEXA_EXTENDIDA.len()) {
                (0, false) => format!(
                    "{quien} {} coincide con {caja}, el sitio donde se une el represor LexA.",
                    texto(consenso)
                ),
                (0, true) => format!(
                    "{quien} {} contiene {caja}, el sitio donde se une el represor LexA.",
                    texto(consenso)
                ),
                (_, false) => format!(
                    "{quien} {} se parece a {caja}, el sitio donde se une el represor LexA: difiere en {diferencias} de {comparadas} letras, algo que al azar pasaría menos de 1 vez en 100.",
                    texto(consenso)
                ),
                (_, true) => format!(
                    "{quien} {} contiene algo muy parecido a {caja}, el sitio donde se une el represor LexA: difiere en {diferencias} de {comparadas} letras, algo que al azar pasaría menos de 1 vez en 100.",
                    texto(consenso)
                ),
            };
            if palindromo {
                frase.push_str(" Como ella, es casi un palíndromo: LexA se une como dímero, una mitad sobre cada brazo CTGT…ACAG.");
            }
            salida.explicar(frase);
        }
        None if palindromo => {
            let parecido = match pares_no_complementarios(consenso) {
                0 => "es igual a su complemento inverso, un palíndromo perfecto".to_string(),
                1 => "es casi igual a su complemento inverso (solo 1 par no encaja)".to_string(),
                pares => {
                    format!(
                        "es casi igual a su complemento inverso (solo {pares} pares no encajan)"
                    )
                }
            };
            salida.explicar(format!(
                "{quien} {} {parecido}: los palíndromos así suelen ser sitios donde se une un dímero de proteína, como LexA en la caja SOS.",
                texto(consenso)
            ));
        }
        None => {}
    }
}

fn explicar_puntuacion(salida: &mut Salida, analisis: &AnalisisMotivos, t: usize) {
    let k = analisis.consenso.len();
    let peor = peor_puntuacion(t, k);
    let fraccion = if peor == 0 {
        0.0
    } else {
        analisis.puntos as f64 / peor as f64
    };
    salida.explicar(format!(
        "La puntuación es {}: son las letras que no coinciden con el consenso {}. Con {t} motivos de {}, el peor caso posible (bases lo más repartidas que se puede) sería {}, así que estás en el {} de ese peor caso.",
        miles(analisis.puntos),
        texto(&analisis.consenso),
        contar(k, "letra", "letras"),
        miles(peor),
        porcentaje(fraccion, 0)
    ));
}

fn explicar_conservacion(salida: &mut Salida, analisis: &AnalisisMotivos, t: usize) {
    let k = analisis.consenso.len();
    let total: f64 = analisis.entropia.iter().sum();
    let bits = informacion(&analisis.entropia);
    let conservadas: Vec<String> = bits
        .iter()
        .enumerate()
        .filter(|(_, valor)| **valor >= 1.0)
        .map(|(columna, _)| (columna + 1).to_string())
        .collect();
    let mut frase = format!(
        "La entropía total es {} bits de un máximo de {} (2 por columna): cuanto más baja, más conservado está el motivo.",
        decimal(total, 2),
        2 * k
    );
    if conservadas.is_empty() {
        frase.push_str(" Ninguna columna llega a 1 bit de información (la mitad del máximo): el motivo está poco conservado o los motivos no están alineados.");
    } else if conservadas.len() == k {
        frase.push_str(" Todas las columnas pasan de 1 bit de información (la mitad del máximo): el motivo está muy conservado.");
    } else {
        frase.push_str(&format!(
            " Las columnas {} (contando desde 1) pasan de 1 bit de información, la mitad del máximo: ahí está el núcleo del motivo; las demás varían más.",
            conservadas.join(", ")
        ));
    }
    salida.explicar(frase);
    let ruido = 3.0 / (2.0 * std::f64::consts::LN_2 * t as f64);
    if ruido >= 0.1 {
        salida.explicar(format!(
            "Con solo {t} motivos, una columna sin ninguna preferencia ya mostraría unos {} bits por puro azar (corrección de muestra pequeña de Schneider): no te fíes de columnas que apenas superan esa cifra.",
            decimal(ruido, 2)
        ));
    }
}

fn agregar_logo_y_perfil(salida: &mut Salida, analisis: &AnalisisMotivos) {
    let categorias: Vec<String> = analisis
        .consenso
        .iter()
        .enumerate()
        .map(|(columna, &letra)| format!("{} {}", columna + 1, letra as char))
        .collect();
    salida.agregar(Bloque::barras(
        "Logo: información de cada columna",
        "Columna (desde 1) y letra del consenso",
        "Bits (2 = siempre la misma letra, 0 = al azar)",
        categorias,
        informacion(&analisis.entropia),
    ));
    let columnas: Vec<String> = (1..=analisis.consenso.len())
        .map(|columna| columna.to_string())
        .collect();
    salida.agregar(Bloque::calor(
        "Perfil: probabilidad de cada base por columna",
        LETRAS.iter().map(|letra| letra.to_string()).collect(),
        columnas,
        analisis.perfil.clone(),
    ));
}

fn coordenada(inicio: Option<usize>, posicion: usize) -> String {
    match inicio {
        Some(inicio) => miles(inicio + posicion),
        None => "—".to_string(),
    }
}

fn mostrar_motivos(
    salida: &mut Salida,
    conjunto: &Conjunto,
    posiciones: &[usize],
    k: usize,
) -> Result<AnalisisMotivos, String> {
    let encontrados = motivos::motivos_en(&conjunto.textos, posiciones, k);
    let analisis = motivos::analizar_motivos(&encontrados)?;
    salida.agregar(Bloque::valor("Consenso", texto(&analisis.consenso)));
    salida.agregar(Bloque::valor(
        "Puntuación",
        format!(
            "{} (peor caso posible: {})",
            miles(analisis.puntos),
            miles(peor_puntuacion(encontrados.len(), k))
        ),
    ));
    let mut filas_alineadas = Vec::new();
    let mut filas = Vec::new();
    for (indice, motivo) in encontrados.iter().enumerate() {
        filas_alineadas.push(FilaAlineada {
            nombre: conjunto.nombres[indice].clone(),
            secuencia: texto(motivo),
        });
        filas.push(vec![
            conjunto.nombres[indice].clone(),
            miles(posiciones[indice]),
            coordenada(conjunto.inicios[indice], posiciones[indice]),
            texto(motivo),
            motivos::distancia_hamming(motivo, &analisis.consenso).to_string(),
        ]);
    }
    salida.agregar(Bloque::alineamiento("Motivos encontrados", filas_alineadas));
    salida.agregar(Bloque::tabla(
        "Dónde está cada motivo",
        &[
            "Secuencia",
            "Posición (contando desde 0)",
            "En E. coli (contando desde 1)",
            "Motivo",
            "Diferencias con el consenso",
        ],
        filas,
    ));
    agregar_logo_y_perfil(salida, &analisis);
    Ok(analisis)
}

fn leer_enumeracion(entrada: &Entrada) -> Result<(Conjunto, usize, usize), String> {
    let conjunto = leer_conjunto(entrada, "adn", 2)?;
    let k = leer_k(entrada, &conjunto.textos)?;
    if k > K_MAXIMA_EXHAUSTIVA {
        return Err(format!(
            "k = {k} es demasiado grande: la enumeración marca cada uno de los 4^k patrones posibles y el máximo es k = {K_MAXIMA_EXHAUSTIVA} (4^10 = 1 048 576)."
        ));
    }
    let d = entrada.natural("d")?;
    if d > D_MAXIMA {
        return Err(format!(
            "d = {d} es demasiado grande: el vecindario de cada k-mer crece como C(k, d) · 3^d y el máximo es d = {D_MAXIMA}."
        ));
    }
    if d > k {
        return Err(format!("d = {d} no puede ser mayor que k = {k}."));
    }
    let vecindario = tamano_del_vecindario(k, d);
    let posibles = 4f64.powi(k as i32);
    let mas_larga = conjunto.textos.iter().map(Vec::len).max().unwrap_or(0) as f64;
    revisar_costo(
        total_de_ventanas(&conjunto.textos, k) * vecindario
            + conjunto.textos.len() as f64 * posibles
            + MOTIVOS_EN_TABLA as f64 * conjunto.textos.len() as f64 * mas_larga * k as f64,
        PRESUPUESTO_ENUMERACION,
    )?;
    Ok((conjunto, k, d))
}

fn enumeracion(entrada: &Entrada) -> Result<Salida, String> {
    let (conjunto, k, d) = leer_enumeracion(entrada)?;
    let vecindario = tamano_del_vecindario(k, d);
    let posibles = 4f64.powi(k as i32);
    let encontrados = motivos::enumerar_motivos(&conjunto.textos, k, d)?;
    let mut esperados = posibles;
    for texto_actual in &conjunto.textos {
        esperados *= probabilidad_en_texto(k, d, texto_actual.len());
    }
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        "Motivos (k, d) encontrados",
        miles(encontrados.len()),
    ));
    salida.agregar(Bloque::valor(
        "Esperados por puro azar",
        aproximado(esperados),
    ));
    salida.agregar(Bloque::lista(
        "Motivos (k, d)",
        encontrados.iter().map(|motivo| texto(motivo)),
    ));
    let mut filas = Vec::new();
    for motivo in encontrados.iter().take(MOTIVOS_EN_TABLA) {
        let mut exactas = 0;
        let mut total = 0;
        for texto_actual in &conjunto.textos {
            let diferencias = motivos::distancia_patron_texto(motivo, texto_actual);
            total += diferencias;
            if diferencias == 0 {
                exactas += 1;
            }
        }
        filas.push(vec![
            texto(motivo),
            format!("{exactas} de {}", conjunto.textos.len()),
            total.to_string(),
        ]);
    }
    if !filas.is_empty() {
        salida.agregar(Bloque::tabla(
            "Cómo aparece cada motivo",
            &[
                "Motivo",
                "Secuencias con copia exacta",
                "Diferencias sumadas en todas",
            ],
            filas,
        ));
    }
    if encontrados.is_empty() {
        salida.explicar(format!(
            "Ningún {k}-mer aparece en todas las secuencias con como mucho {}. Prueba a subir d o a bajar k.",
            contar(d, "diferencia", "diferencias")
        ));
    } else if cola_de_poisson(esperados, encontrados.len()) < AZAR_TOLERADO {
        salida.explicar(format!(
            "Encontraste {} y por puro azar esperarías unos {} (bases al azar, cada una con probabilidad 1/4): encontrar tantos o más saldría por azar menos de 1 vez en 100 (p < 0.01), así que estas secuencias probablemente comparten una señal real.",
            contar(encontrados.len(), "motivo", "motivos"),
            aproximado(esperados)
        ));
    } else {
        salida.explicar(format!(
            "Encontraste {} y por puro azar esperarías unos {}: es compatible con el azar. Con estos k y d casi cualquier patrón aparece en todas; baja d o sube k para que el resultado signifique algo.",
            contar(encontrados.len(), "motivo", "motivos"),
            aproximado(esperados)
        ));
    }
    salida.explicar(format!(
        "Alrededor de cada {k}-mer hay {} a como mucho {} (él incluido), de {} posibles (4^{k}); por eso subir d dispara el número de motivos que pasan el filtro.",
        contar(vecindario as usize, "patrón", "patrones"),
        contar(d, "diferencia", "diferencias"),
        miles(posibles as usize)
    ));
    let mas_parecido_a_lexa = encontrados
        .iter()
        .filter_map(|motivo| parecido_a_lexa(motivo).map(|diferencias| (diferencias, motivo)))
        .min_by_key(|(diferencias, _)| *diferencias);
    if let Some((_, motivo)) = mas_parecido_a_lexa {
        explicar_lexa_y_palindromo(&mut salida, motivo, "El motivo");
    }
    if encontrados.len() > MOTIVOS_EN_TABLA {
        salida.explicar(format!(
            "La tabla muestra solo los primeros {MOTIVOS_EN_TABLA} en orden alfabético; la lista los tiene todos (hasta 2000)."
        ));
    }
    Ok(salida)
}

fn leer_matriz_de_motivos(entrada: &Entrada) -> Result<Conjunto, String> {
    let conjunto = leer_conjunto(entrada, "motivos", 2)?;
    let k = conjunto.textos[0].len();
    if conjunto.textos.iter().any(|motivo| motivo.len() != k) {
        return Err("Todos los motivos deben medir lo mismo.".to_string());
    }
    if k > LARGO_MAXIMO_DE_MOTIVO {
        return Err(format!(
            "Cada motivo mide {} pb; un motivo regulatorio rara vez pasa de unas decenas de bases y aquí el máximo es {LARGO_MAXIMO_DE_MOTIVO}.",
            miles(k)
        ));
    }
    Ok(conjunto)
}

fn perfil(entrada: &Entrada) -> Result<Salida, String> {
    let conjunto = leer_matriz_de_motivos(entrada)?;
    let k = conjunto.textos[0].len();
    let vistas: Vec<&[u8]> = conjunto.textos.iter().map(Vec::as_slice).collect();
    let analisis = motivos::analizar_motivos(&vistas)?;
    let t = vistas.len();
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Consenso", texto(&analisis.consenso)));
    salida.agregar(Bloque::valor(
        "Puntuación",
        format!(
            "{} (peor caso posible: {})",
            miles(analisis.puntos),
            miles(peor_puntuacion(t, k))
        ),
    ));
    salida.agregar(Bloque::valor(
        "Entropía total",
        format!(
            "{} bits (máximo {})",
            decimal(analisis.entropia.iter().sum(), 2),
            2 * k
        ),
    ));
    salida.agregar(Bloque::alineamiento(
        "Matriz de motivos",
        conjunto
            .nombres
            .iter()
            .zip(&conjunto.textos)
            .map(|(nombre, motivo)| FilaAlineada {
                nombre: nombre.clone(),
                secuencia: texto(motivo),
            })
            .collect(),
    ));
    let mut encabezados = vec!["Base".to_string()];
    encabezados.extend((1..=k).map(|columna| columna.to_string()));
    let columnas: Vec<&str> = encabezados.iter().map(String::as_str).collect();
    let filas: Vec<Vec<String>> = LETRAS
        .iter()
        .zip(&analisis.conteos)
        .map(|(letra, conteos_de_la_base)| {
            let mut celdas = vec![letra.to_string()];
            celdas.extend(conteos_de_la_base.iter().map(|conteo| conteo.to_string()));
            celdas
        })
        .collect();
    salida.agregar(Bloque::tabla(
        "Conteos por columna (desde 1)",
        &columnas,
        filas,
    ));
    agregar_logo_y_perfil(&mut salida, &analisis);
    explicar_puntuacion(&mut salida, &analisis, t);
    explicar_conservacion(&mut salida, &analisis, t);
    explicar_lexa_y_palindromo(&mut salida, &analisis.consenso, "El consenso");
    Ok(salida)
}

fn leer_distancia(entrada: &Entrada) -> Result<(Vec<u8>, Conjunto), String> {
    let patron = entrada.adn("patron")?;
    let conjunto = leer_conjunto(entrada, "adn", 1)?;
    let k = patron.len();
    if let Some(indice) = conjunto.textos.iter().position(|texto| texto.len() < k) {
        return Err(format!(
            "El patrón mide {} pb y la secuencia «{}» solo {} pb.",
            miles(k),
            conjunto.nombres[indice],
            miles(conjunto.textos[indice].len())
        ));
    }
    revisar_costo(
        total_de_ventanas(&conjunto.textos, k) * k as f64 * 2.0,
        PRESUPUESTO_LINEAL,
    )?;
    Ok((patron, conjunto))
}

fn distancia(entrada: &Entrada) -> Result<Salida, String> {
    let (patron, conjunto) = leer_distancia(entrada)?;
    let k = patron.len();
    let total = motivos::distancia_patron_textos(&patron, &conjunto.textos)?;
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        "Distancia d(patrón, secuencias)",
        miles(total),
    ));
    let mut filas = Vec::new();
    let mut alineadas = vec![FilaAlineada {
        nombre: "patrón".to_string(),
        secuencia: texto(&patron),
    }];
    let mut categorias = Vec::new();
    let mut valores = Vec::new();
    let mut raras = Vec::new();
    let mut exactas = 0;
    for (indice, texto_actual) in conjunto.textos.iter().enumerate() {
        let (posicion, diferencias) = motivos::coincidencia_mas_cercana(&patron, texto_actual);
        let cercano = &texto_actual[posicion..posicion + k];
        let azar = probabilidad_en_texto(k, diferencias, texto_actual.len());
        if diferencias == 0 {
            exactas += 1;
        }
        if azar < AZAR_TOLERADO {
            raras.push(conjunto.nombres[indice].clone());
        }
        filas.push(vec![
            conjunto.nombres[indice].clone(),
            miles(posicion),
            coordenada(conjunto.inicios[indice], posicion),
            texto(cercano),
            diferencias.to_string(),
            cientifico(azar),
        ]);
        alineadas.push(FilaAlineada {
            nombre: conjunto.nombres[indice].clone(),
            secuencia: texto(cercano),
        });
        categorias.push(conjunto.nombres[indice].clone());
        valores.push(diferencias as f64);
    }
    salida.agregar(Bloque::tabla(
        "El k-mer más parecido en cada secuencia",
        &[
            "Secuencia",
            "Posición (contando desde 0)",
            "En E. coli (contando desde 1)",
            "k-mer más parecido",
            "Diferencias",
            "Probabilidad de algo así por azar",
        ],
        filas,
    ));
    salida.agregar(Bloque::alineamiento("Patrón y coincidencias", alineadas));
    salida.agregar(Bloque::barras(
        "Diferencias en cada secuencia",
        "Secuencia",
        "Diferencias con el k-mer más parecido",
        categorias,
        valores,
    ));
    let t = conjunto.textos.len();
    salida.explicar(format!(
        "La distancia es {}: la suma de las diferencias con el k-mer más parecido al patrón en {} (en promedio {} por secuencia, sobre {}).",
        miles(total),
        contar(t, "secuencia", "secuencias"),
        decimal(total as f64 / t as f64, 1),
        contar(k, "letra", "letras")
    ));
    if exactas > 0 {
        salida.explicar(format!(
            "El patrón aparece sin ningún error en {exactas} de {}.",
            contar(t, "secuencia", "secuencias")
        ));
    }
    if raras.is_empty() {
        salida.explicar("En ninguna secuencia la coincidencia es mejor de lo que daría el azar con p < 0.01 (bases al azar con probabilidad 1/4): no hay pruebas de que el patrón esté realmente ahí.");
    } else if raras.len() == t {
        salida.explicar("En cada secuencia, una coincidencia tan buena saldría por azar menos de 1 vez en 100 (p < 0.01, bases al azar con probabilidad 1/4): probablemente todas llevan una copia real del patrón.");
    } else {
        let verbo = if raras.len() == 1 { "lleva" } else { "llevan" };
        salida.explicar(format!(
            "En {} de {t} secuencias ({}) una coincidencia tan buena saldría por azar menos de 1 vez en 100 (p < 0.01): probablemente {verbo} una copia real. En las demás, la coincidencia no se distingue del azar.",
            raras.len(),
            raras.join(", ")
        ));
    }
    explicar_lexa_y_palindromo(&mut salida, &patron, "El patrón");
    Ok(salida)
}

fn leer_mediana(entrada: &Entrada) -> Result<(Conjunto, usize), String> {
    let conjunto = leer_conjunto(entrada, "adn", 2)?;
    let k = leer_k(entrada, &conjunto.textos)?;
    if k > K_MAXIMA_EXHAUSTIVA {
        return Err(format!(
            "k = {k} es demasiado grande: la cadena mediana revisa los 4^k patrones posibles y cada letra más multiplica el trabajo por 4. Con k = 10 ya son 1 048 576; el máximo es k = {K_MAXIMA_EXHAUSTIVA}."
        ));
    }
    let posibles = 4f64.powi(k as i32);
    revisar_costo(
        conjunto.textos.len() as f64 * posibles * 3.0 * k as f64
            + total_de_ventanas(&conjunto.textos, k) * k as f64,
        PRESUPUESTO,
    )?;
    Ok((conjunto, k))
}

fn mediana(entrada: &Entrada) -> Result<Salida, String> {
    let (conjunto, k) = leer_mediana(entrada)?;
    let posibles = 4f64.powi(k as i32);
    let resultado = motivos::cadena_mediana(&conjunto.textos, k)?;
    let elegida = resultado.patrones[0].clone();
    let t = conjunto.textos.len();
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Cadena mediana", texto(&elegida)));
    salida.agregar(Bloque::valor("Distancia total", miles(resultado.distancia)));
    salida.agregar(Bloque::valor(
        "Patrones revisados",
        miles(posibles as usize),
    ));
    salida.agregar(Bloque::lista(
        "Medianas empatadas (orden alfabético)",
        resultado.patrones.iter().map(|patron| texto(patron)),
    ));
    let mut filas = Vec::new();
    let mut alineadas = vec![FilaAlineada {
        nombre: "mediana".to_string(),
        secuencia: texto(&elegida),
    }];
    for (indice, texto_actual) in conjunto.textos.iter().enumerate() {
        let (posicion, diferencias) = motivos::coincidencia_mas_cercana(&elegida, texto_actual);
        let cercano = &texto_actual[posicion..posicion + k];
        filas.push(vec![
            conjunto.nombres[indice].clone(),
            miles(posicion),
            coordenada(conjunto.inicios[indice], posicion),
            texto(cercano),
            diferencias.to_string(),
        ]);
        alineadas.push(FilaAlineada {
            nombre: conjunto.nombres[indice].clone(),
            secuencia: texto(cercano),
        });
    }
    salida.agregar(Bloque::tabla(
        "Dónde está la mediana en cada secuencia",
        &[
            "Secuencia",
            "Posición (contando desde 0)",
            "En E. coli (contando desde 1)",
            "k-mer más parecido",
            "Diferencias",
        ],
        filas,
    ));
    salida.agregar(Bloque::alineamiento("Mediana y sus apariciones", alineadas));
    salida.explicar(format!(
        "Se revisaron los {} patrones posibles de {} (4^{k}); cada letra más multiplica ese número por 4, por eso este método exacto solo llega hasta k = {K_MAXIMA_EXHAUSTIVA}.",
        miles(posibles as usize),
        contar(k, "letra", "letras")
    ));
    let promedio = resultado.distancia as f64 / t as f64;
    if resultado.distancia == 0 {
        salida.explicar(format!(
            "La mediana {} aparece exacta en las {t} secuencias.",
            texto(&elegida)
        ));
    } else {
        salida.explicar(format!(
            "La mediana {} está a {} en total: en promedio {} por secuencia, sobre {}.",
            texto(&elegida),
            contar(resultado.distancia, "diferencia", "diferencias"),
            decimal(promedio, 1),
            contar(k, "letra", "letras")
        ));
    }
    if resultado.empatados >= 10 {
        salida.explicar(format!(
            "Hay {} patrones empatados con la misma distancia: con tantos empates no hay un único motivo claro; prueba con otro k. La lista muestra los primeros en orden alfabético.",
            miles(resultado.empatados)
        ));
    } else if resultado.empatados > 1 {
        salida.explicar(format!(
            "Hay {} patrones empatados con la misma distancia; el libro acepta cualquiera y aquí se elige el primero en orden alfabético.",
            resultado.empatados
        ));
    }
    explicar_lexa_y_palindromo(&mut salida, &elegida, "La mediana");
    Ok(salida)
}

fn leer_perfil(entrada: &Entrada) -> Result<Vec<Vec<f64>>, String> {
    let perfil = entrada.matriz("perfil")?;
    if perfil.len() != 4 {
        return Err(format!(
            "El perfil debe tener 4 filas (A, C, G y T) y tiene {}.",
            perfil.len()
        ));
    }
    let k = perfil[0].len();
    if k > K_MAXIMA_PERFIL {
        return Err(format!(
            "El perfil tiene {k} columnas; el máximo es {K_MAXIMA_PERFIL}."
        ));
    }
    for columna in 0..k {
        if perfil
            .iter()
            .any(|fila| !(0.0..=1.0).contains(&fila[columna]))
        {
            return Err(format!(
                "En la columna {} del perfil hay un valor fuera de 0 a 1.",
                columna + 1
            ));
        }
        let suma: f64 = perfil.iter().map(|fila| fila[columna]).sum();
        if (suma - 1.0).abs() > 0.02 {
            return Err(format!(
                "La columna {} del perfil suma {} y debería sumar 1: cada columna es una distribución de probabilidad.",
                columna + 1,
                decimal(suma, 3)
            ));
        }
    }
    Ok(perfil)
}

type LecturaMasProbable = (Vec<u8>, Option<usize>, Vec<Vec<f64>>);

fn leer_mas_probable(entrada: &Entrada) -> Result<LecturaMasProbable, String> {
    let perfil = leer_perfil(entrada)?;
    let (adn, inicio) = leer_texto(entrada, "texto")?;
    let k = perfil[0].len();
    if adn.len() < k {
        return Err(format!(
            "El texto mide {} pb y el perfil tiene {k} columnas.",
            miles(adn.len())
        ));
    }
    revisar_costo(ventanas(adn.len(), k) * k as f64 * 2.0, PRESUPUESTO_LINEAL)?;
    Ok((adn, inicio, perfil))
}

fn mas_probable(entrada: &Entrada) -> Result<Salida, String> {
    let (adn, inicio, perfil) = leer_mas_probable(entrada)?;
    let k = perfil[0].len();
    let posicion = motivos::posicion_mas_probable(&adn, &perfil);
    let elegido = &adn[posicion..posicion + k];
    let mayor = motivos::probabilidad_kmer(elegido, &perfil);
    let cantidad_de_ventanas = adn.len() - k + 1;
    let grupo = cantidad_de_ventanas.div_ceil(2000).max(1);
    let mut puntos: Vec<[f64; 2]> = Vec::new();
    let mut mejores: Vec<(f64, usize)> = Vec::new();
    let mut empatadas = 0;
    let mut segunda: f64 = 0.0;
    let mut maximo_del_grupo: Option<[f64; 2]> = None;
    for (indice, ventana) in adn.windows(k).enumerate() {
        let probabilidad = motivos::probabilidad_kmer(ventana, &perfil);
        if probabilidad == mayor {
            empatadas += 1;
        } else {
            segunda = segunda.max(probabilidad);
        }
        if probabilidad > 0.0 {
            let punto = [indice as f64, probabilidad.log10()];
            if maximo_del_grupo.is_none_or(|actual| punto[1] > actual[1]) {
                maximo_del_grupo = Some(punto);
            }
        }
        let cierra_grupo = (indice + 1) % grupo == 0 || indice + 1 == cantidad_de_ventanas;
        if cierra_grupo && let Some(punto) = maximo_del_grupo.take() {
            puntos.push(punto);
        }
        if mejores.len() < 10 || probabilidad > mejores[mejores.len() - 1].0 {
            mejores.push((probabilidad, indice));
            mejores.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
            mejores.truncate(10);
        }
    }
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("k-mer más probable", texto(elegido)));
    salida.agregar(Bloque::valor(
        "Posición (contando desde 0)",
        miles(posicion),
    ));
    if let Some(inicio) = inicio {
        salida.agregar(Bloque::valor(
            "Posición en E. coli (contando desde 1)",
            miles(inicio + posicion),
        ));
    }
    salida.agregar(Bloque::valor("Probabilidad", cientifico(mayor)));
    salida.agregar(Bloque::linea(
        "Probabilidad de cada ventana según el perfil",
        "Posición de la ventana (contando desde 0)",
        "log10 de la probabilidad",
        &puntos,
        vec![Marca {
            x: posicion as f64,
            texto: texto(elegido),
        }],
    ));
    salida.agregar(Bloque::tabla(
        "Las ventanas más probables",
        &[
            "Posición (contando desde 0)",
            "k-mer",
            "Probabilidad",
            "Veces la del azar (0.25 por letra)",
        ],
        mejores
            .iter()
            .map(|(probabilidad, indice)| {
                vec![
                    miles(*indice),
                    texto(&adn[*indice..*indice + k]),
                    cientifico(*probabilidad),
                    aproximado(probabilidad / 0.25f64.powi(k as i32)),
                ]
            })
            .collect(),
    ));
    if mayor <= 0.0 {
        salida.explicar("Todas las ventanas tienen probabilidad 0: el perfil tiene ceros y cada ventana contiene alguna letra imposible en alguna columna. Se devuelve la primera, como manda el libro; para evitarlo, construye el perfil con pseudoconteos.");
        return Ok(salida);
    }
    salida.explicar(format!(
        "El k-mer más probable es {} en la posición {} (contando desde 0), con probabilidad {}: {} veces más que un k-mer cualquiera si cada letra tuviera probabilidad 0.25.",
        texto(elegido),
        miles(posicion),
        cientifico(mayor),
        aproximado(mayor / 0.25f64.powi(k as i32))
    ));
    if empatadas > 1 {
        salida.explicar(format!(
            "Hay {empatadas} ventanas empatadas con la probabilidad máxima; se elige la primera, como en el libro."
        ));
    } else if segunda <= 0.0 || mayor >= 10.0 * segunda {
        salida.explicar(format!(
            "Destaca con claridad: es al menos 10 veces más probable que la segunda mejor ventana ({}). Es un buen candidato a ser una copia real del motivo.",
            cientifico(segunda)
        ));
    } else {
        salida.explicar(format!(
            "No destaca mucho: la segunda mejor ventana tiene probabilidad {}, menos de 10 veces por debajo. Mira la tabla de las ventanas más probables antes de quedarte con una.",
            cientifico(segunda)
        ));
    }
    if let Some(inicio) = inicio {
        salida.explicar(format!(
            "En el genoma de E. coli, que se numera desde 1, esa ventana empieza en la posición {} (inicio de la región {} + posición {}).",
            miles(inicio + posicion),
            miles(inicio),
            miles(posicion)
        ));
    }
    explicar_lexa_y_palindromo(&mut salida, elegido, "El k-mer");
    Ok(salida)
}

fn costo_codicioso(textos: &[Vec<u8>], k: usize) -> f64 {
    let resto: f64 = total_de_ventanas(&textos[1..], k);
    ventanas(textos[0].len(), k) * (resto * k as f64 + (textos.len() * textos.len() * k) as f64)
}

fn salida_codiciosa(
    conjunto: &Conjunto,
    k: usize,
    hallado: &MotivosEncontrados,
) -> Result<(Salida, AnalisisMotivos), String> {
    let mut salida = Salida::default();
    let analisis = mostrar_motivos(&mut salida, conjunto, &hallado.posiciones, k)?;
    let puntos: Vec<[f64; 2]> = hallado
        .puntos_por_intento
        .iter()
        .enumerate()
        .map(|(inicio, &puntos)| [inicio as f64, puntos as f64])
        .collect();
    salida.agregar(Bloque::linea(
        "Puntuación según el k-mer inicial de la primera secuencia",
        "Posición del k-mer inicial (contando desde 0)",
        "Puntuación de los motivos",
        &puntos,
        vec![Marca {
            x: hallado.posiciones[0] as f64,
            texto: "mejor".to_string(),
        }],
    ));
    Ok((salida, analisis))
}

fn filas_sin_probabilidad(textos: &[Vec<u8>], posiciones: &[usize], k: usize) -> usize {
    let mut cuantas = 0;
    for i in 1..textos.len() {
        let elegidos = motivos::motivos_en(&textos[..i], &posiciones[..i], k);
        let perfil_previo = motivos::perfil_de_motivos(&elegidos, 0.0);
        let motivo = &textos[i][posiciones[i]..posiciones[i] + k];
        if motivos::probabilidad_kmer(motivo, &perfil_previo) <= 0.0 {
            cuantas += 1;
        }
    }
    cuantas
}

fn codiciosa(entrada: &Entrada) -> Result<Salida, String> {
    let conjunto = leer_conjunto(entrada, "adn", 2)?;
    let k = leer_k_de_perfil(entrada, &conjunto.textos)?;
    revisar_costo(costo_codicioso(&conjunto.textos, k), PRESUPUESTO)?;
    let hallado = motivos::busqueda_codiciosa(&conjunto.textos, k, 0.0)?;
    let (mut salida, analisis) = salida_codiciosa(&conjunto, k, &hallado)?;
    let t = conjunto.textos.len();
    let sin_probabilidad = filas_sin_probabilidad(&conjunto.textos, &hallado.posiciones, k);
    if sin_probabilidad > 0 {
        salida.explicar(format!(
            "Sin pseudoconteos, en {sin_probabilidad} de {} (todas menos la primera), cada ventana tenía probabilidad 0 porque llevaba alguna letra nunca vista en esa columna, y se tomó la primera ventana por defecto. Compara con la versión con pseudoconteos.",
            contar(t - 1, "secuencia completada", "secuencias completadas")
        ));
    } else {
        salida.explicar("En los motivos elegidos, ninguna secuencia se quedó con todas sus ventanas en probabilidad 0, así que los ceros del perfil no forzaron ninguna elección.");
    }
    explicar_puntuacion(&mut salida, &analisis, t);
    explicar_conservacion(&mut salida, &analisis, t);
    explicar_lexa_y_palindromo(&mut salida, &analisis.consenso, "El consenso");
    Ok(salida)
}

fn codiciosa_laplace(entrada: &Entrada) -> Result<Salida, String> {
    let conjunto = leer_conjunto(entrada, "adn", 2)?;
    let k = leer_k_de_perfil(entrada, &conjunto.textos)?;
    revisar_costo(2.0 * costo_codicioso(&conjunto.textos, k), PRESUPUESTO)?;
    let hallado = motivos::busqueda_codiciosa(&conjunto.textos, k, 1.0)?;
    let sin_pseudoconteos = motivos::busqueda_codiciosa(&conjunto.textos, k, 0.0)?;
    let (mut salida, analisis) = salida_codiciosa(&conjunto, k, &hallado)?;
    let t = conjunto.textos.len();
    if sin_pseudoconteos.puntos > hallado.puntos {
        salida.explicar(format!(
            "Sin pseudoconteos, la misma búsqueda solo llega a una puntuación de {} (aquí {}): sumar 1 a cada conteo evita que un único cero en el perfil descarte un k-mer bueno.",
            miles(sin_pseudoconteos.puntos),
            miles(hallado.puntos)
        ));
    } else if sin_pseudoconteos.puntos == hallado.puntos {
        salida.explicar(format!(
            "Sin pseudoconteos se llega a la misma puntuación ({}): en estas secuencias los ceros del perfil no estorbaron.",
            miles(hallado.puntos)
        ));
    } else {
        salida.explicar(format!(
            "Sin pseudoconteos se llega incluso a {} (aquí {}): la búsqueda codiciosa es rápida pero no garantiza el mejor resultado posible.",
            miles(sin_pseudoconteos.puntos),
            miles(hallado.puntos)
        ));
    }
    explicar_puntuacion(&mut salida, &analisis, t);
    explicar_conservacion(&mut salida, &analisis, t);
    explicar_lexa_y_palindromo(&mut salida, &analisis.consenso, "El consenso");
    Ok(salida)
}

fn histograma_de_intentos(salida: &mut Salida, puntos_por_intento: &[usize], eje_x: &str) {
    let mut cuentas: BTreeMap<usize, usize> = BTreeMap::new();
    for &puntos in puntos_por_intento {
        *cuentas.entry(puntos).or_insert(0) += 1;
    }
    salida.agregar(Bloque::barras(
        "Puntuación final de cada intento",
        eje_x,
        "Intentos",
        cuentas.keys().map(|puntos| puntos.to_string()).collect(),
        cuentas.values().map(|&veces| veces as f64).collect(),
    ));
}

fn explicar_azar(salida: &mut Salida, hallado: &MotivosEncontrados, semilla: u64, intentos: &str) {
    let total = hallado.puntos_por_intento.len();
    let veces = hallado
        .puntos_por_intento
        .iter()
        .filter(|&&puntos| puntos == hallado.puntos)
        .count();
    if total == 1 {
        salida.explicar(format!(
            "Este método depende del azar y aquí hubo un solo intento (semilla {semilla}): con otra semilla el resultado puede cambiar. Haz varios intentos y quédate con el mejor."
        ));
    } else if veces == 1 {
        salida.explicar(format!(
            "Este método depende del azar: con la semilla {semilla}, la mejor puntuación ({}) salió en solo 1 de {total} {intentos}. Con otra semilla podrías no verla: haz más intentos o prueba varias semillas.",
            miles(hallado.puntos)
        ));
    } else {
        salida.explicar(format!(
            "Este método depende del azar: con la semilla {semilla}, la mejor puntuación ({}) salió en {veces} de {total} {intentos}. Que varios arranques al azar lleguen al mismo sitio es buena señal, pero otra semilla puede dar otro resultado: repite con varias.",
            miles(hallado.puntos)
        ));
    }
}

fn leer_aleatoria(entrada: &Entrada) -> Result<(Conjunto, usize, usize, u64), String> {
    let conjunto = leer_conjunto(entrada, "adn", 2)?;
    let k = leer_k_de_perfil(entrada, &conjunto.textos)?;
    let repeticiones = entrada.natural("repeticiones")?;
    if repeticiones == 0 {
        return Err("Hace falta al menos una repetición.".to_string());
    }
    let semilla = leer_semilla(entrada)?;
    revisar_costo(
        repeticiones as f64
            * ITERACIONES_TIPICAS
            * total_de_ventanas(&conjunto.textos, k)
            * k as f64,
        PRESUPUESTO,
    )?;
    Ok((conjunto, k, repeticiones, semilla))
}

fn aleatoria(entrada: &Entrada) -> Result<Salida, String> {
    let (conjunto, k, repeticiones, semilla) = leer_aleatoria(entrada)?;
    let hallado = motivos::busqueda_aleatoria(&conjunto.textos, k, repeticiones, semilla)?;
    let mut salida = Salida::default();
    let analisis = mostrar_motivos(&mut salida, &conjunto, &hallado.posiciones, k)?;
    histograma_de_intentos(&mut salida, &hallado.puntos_por_intento, "Puntuación");
    let t = conjunto.textos.len();
    explicar_azar(&mut salida, &hallado, semilla, "repeticiones");
    explicar_puntuacion(&mut salida, &analisis, t);
    explicar_conservacion(&mut salida, &analisis, t);
    explicar_lexa_y_palindromo(&mut salida, &analisis.consenso, "El consenso");
    Ok(salida)
}

type LecturaGibbs = (Conjunto, usize, usize, usize, u64);

fn leer_gibbs(entrada: &Entrada) -> Result<LecturaGibbs, String> {
    let conjunto = leer_conjunto(entrada, "adn", 2)?;
    let k = leer_k_de_perfil(entrada, &conjunto.textos)?;
    let iteraciones = entrada.natural("iteraciones")?;
    let reinicios = entrada.natural("reinicios")?;
    if reinicios == 0 {
        return Err("Hace falta al menos un reinicio.".to_string());
    }
    let semilla = leer_semilla(entrada)?;
    let t = conjunto.textos.len();
    let mas_larga = conjunto.textos.iter().map(Vec::len).max().unwrap_or(0);
    revisar_costo(
        reinicios as f64
            * (iteraciones as f64 + 1.0)
            * (ventanas(mas_larga, k) + 3.0 * t as f64)
            * k as f64,
        PRESUPUESTO,
    )?;
    Ok((conjunto, k, iteraciones, reinicios, semilla))
}

fn gibbs(entrada: &Entrada) -> Result<Salida, String> {
    let (conjunto, k, iteraciones, reinicios, semilla) = leer_gibbs(entrada)?;
    let t = conjunto.textos.len();
    let hallado = motivos::muestreo_gibbs(&conjunto.textos, k, iteraciones, reinicios, semilla)?;
    let mut salida = Salida::default();
    let analisis = mostrar_motivos(&mut salida, &conjunto, &hallado.posiciones, k)?;
    histograma_de_intentos(&mut salida, &hallado.puntos_por_intento, "Puntuación");
    explicar_azar(&mut salida, &hallado, semilla, "reinicios");
    explicar_puntuacion(&mut salida, &analisis, t);
    explicar_conservacion(&mut salida, &analisis, t);
    explicar_lexa_y_palindromo(&mut salida, &analisis.consenso, "El consenso");
    Ok(salida)
}

fn redondear(valor: f64) -> f64 {
    format!("{valor:.2e}").parse().unwrap_or(valor)
}

fn perfil_redondeado(perfil_actual: &[Vec<f64>]) -> Vec<Vec<f64>> {
    perfil_actual
        .iter()
        .map(|fila| fila.iter().map(|&valor| redondear(valor)).collect())
        .collect()
}

fn kmers_en(textos: &[Vec<u8>], posiciones: &[usize], k: usize) -> Vec<String> {
    motivos::motivos_en(textos, posiciones, k)
        .iter()
        .map(|motivo| texto(motivo))
        .collect()
}

fn rondas_a_mostrar(
    puntos_por_intento: &[usize],
    puntos_iniciales: usize,
    maximo: usize,
) -> Vec<bool> {
    let rondas = puntos_por_intento.len();
    if rondas <= maximo {
        return vec![true; rondas];
    }
    let mut mostrar = vec![false; rondas];
    let mut elegidas = 0;
    let mut mejor = puntos_iniciales;
    for (ronda, &puntos) in puntos_por_intento.iter().enumerate() {
        if puntos < mejor {
            mejor = puntos;
            if elegidas < maximo {
                mostrar[ronda] = true;
                elegidas += 1;
            }
        }
    }
    let faltan = maximo - elegidas;
    for j in 0..faltan {
        mostrar[j * rondas / faltan] = true;
    }
    mostrar
}

struct Traza<'a> {
    textos: &'a [Vec<u8>],
    nombres: &'a [String],
    k: usize,
    rondas: usize,
    mejores: Vec<usize>,
    mejor_puntuacion: usize,
}

impl Traza<'_> {
    fn estado(
        &self,
        fase: &str,
        ronda: Option<usize>,
        posiciones: &[usize],
        perfil_actual: Option<&[Vec<f64>]>,
        extra: Value,
    ) -> Value {
        let mut estado = json!({
            "fase": fase,
            "k": self.k,
            "pseudoconteo": 1,
            "secuencias": self.textos.iter().map(|secuencia| texto(secuencia)).collect::<Vec<_>>(),
            "nombres": self.nombres,
            "ronda": ronda,
            "rondas": self.rondas,
            "fila": Value::Null,
            "posiciones": posiciones,
            "kmers": kmers_en(self.textos, posiciones, self.k),
            "perfil": perfil_actual.map(perfil_redondeado),
            "probabilidades": Value::Null,
            "elegido": Value::Null,
            "puntos": Value::Null,
            "mejoro": Value::Null,
            "mejores": {
                "posiciones": self.mejores,
                "kmers": kmers_en(self.textos, &self.mejores, self.k),
                "puntos": self.mejor_puntuacion,
            },
        });
        if let (Value::Object(campos), Value::Object(nuevos)) = (&mut estado, extra) {
            for (clave, valor) in nuevos {
                campos.insert(clave, valor);
            }
        }
        estado
    }
}

fn pasos_codiciosa_laplace(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let conjunto = leer_conjunto(entrada, "adn", 2)?;
    let k = leer_k_de_perfil(entrada, &conjunto.textos)?;
    let t = conjunto.textos.len();
    if t > SECUENCIAS_ANIMABLES {
        return Err(format!(
            "Para el paso a paso usa como mucho {SECUENCIAS_ANIMABLES} secuencias; aquí hay {t}."
        ));
    }
    let bases: usize = conjunto.textos.iter().map(Vec::len).sum();
    if bases > BASES_ANIMABLES {
        return Err(format!(
            "Para el paso a paso usa como mucho {} pb en total; aquí hay {}.",
            miles(BASES_ANIMABLES),
            miles(bases)
        ));
    }
    let textos = &conjunto.textos;
    let hallado = motivos::busqueda_codiciosa(textos, k, 1.0)?;
    let rondas = textos[0].len() - k + 1;
    let mut traza = Traza {
        textos,
        nombres: &conjunto.nombres,
        k,
        rondas,
        mejores: vec![0; t],
        mejor_puntuacion: 0,
    };
    traza.mejor_puntuacion = motivos::puntuacion(&motivos::motivos_en(textos, &traza.mejores, k));
    let maximo_de_rondas = (MAXIMO_PASOS - 2) / (t + 1);
    let mostrar = rondas_a_mostrar(
        &hallado.puntos_por_intento,
        traza.mejor_puntuacion,
        maximo_de_rondas,
    );
    let mut pasos = vec![Paso::nuevo(
        "let mut mejores = vec![0; adn.len()];",
        format!(
            "Los mejores motivos empiezan siendo el primer k-mer de cada secuencia, con puntuación {}.",
            traza.mejor_puntuacion
        ),
        traza.estado(
            "inicio",
            None,
            &traza.mejores.clone(),
            None,
            json!({ "puntos": traza.mejor_puntuacion }),
        ),
    )];
    let mut posiciones = vec![0; t];
    for inicio in 0..rondas {
        posiciones[0] = inicio;
        let visible = mostrar[inicio];
        if visible {
            pasos.push(Paso::nuevo(
                "posiciones[0] = inicio;",
                format!(
                    "Ronda {inicio}: se fija el k-mer de la posición {inicio} (contando desde 0) de la primera secuencia."
                ),
                traza.estado(
                    "primer_kmer",
                    Some(inicio),
                    &posiciones[..1],
                    None,
                    json!({
                        "fila": 0,
                        "elegido": { "fila": 0, "posicion": inicio, "kmer": texto(&textos[0][inicio..inicio + k]), "probabilidad": Value::Null },
                    }),
                ),
            ));
        }
        for i in 1..t {
            let elegidos = motivos::motivos_en(&textos[..i], &posiciones[..i], k);
            let perfil_actual = motivos::perfil_de_motivos(&elegidos, 1.0);
            posiciones[i] = motivos::posicion_mas_probable(&textos[i], &perfil_actual);
            if visible {
                let probabilidades: Vec<f64> = textos[i]
                    .windows(k)
                    .map(|ventana| redondear(motivos::probabilidad_kmer(ventana, &perfil_actual)))
                    .collect();
                let elegido = &textos[i][posiciones[i]..posiciones[i] + k];
                let base_del_perfil = if i == 1 {
                    "Con el perfil de la primera secuencia".to_string()
                } else {
                    format!("Con el perfil de las {i} secuencias anteriores")
                };
                pasos.push(Paso::nuevo(
                    "posiciones[i] = posicion_mas_probable(&adn[i], &perfil);",
                    format!(
                        "{base_del_perfil}, el k-mer más probable de la secuencia {} es {} (posición {}, contando desde 0).",
                        i + 1,
                        texto(elegido),
                        posiciones[i]
                    ),
                    traza.estado(
                        "elegir",
                        Some(inicio),
                        &posiciones[..=i],
                        Some(&perfil_actual),
                        json!({
                            "fila": i,
                            "probabilidades": probabilidades,
                            "elegido": {
                                "fila": i,
                                "posicion": posiciones[i],
                                "kmer": texto(elegido),
                                "probabilidad": redondear(motivos::probabilidad_kmer(elegido, &perfil_actual)),
                            },
                        }),
                    ),
                ));
            }
        }
        let actuales = motivos::motivos_en(textos, &posiciones, k);
        let puntos = motivos::puntuacion(&actuales);
        let mejoro = puntos < traza.mejor_puntuacion;
        let anterior = traza.mejor_puntuacion;
        if mejoro {
            traza.mejor_puntuacion = puntos;
            traza.mejores.copy_from_slice(&posiciones);
        }
        if visible {
            let perfil_completo = motivos::perfil_de_motivos(&actuales, 1.0);
            let (linea, nota) = if mejoro {
                (
                    "mejores.copy_from_slice(&posiciones);",
                    format!(
                        "Puntuación {puntos}: mejora la mejor hasta ahora ({anterior}) y la reemplaza."
                    ),
                )
            } else {
                (
                    "if puntos < mejor_puntuacion {",
                    format!("Puntuación {puntos}: no mejora la mejor hasta ahora ({anterior})."),
                )
            };
            pasos.push(Paso::nuevo(
                linea,
                nota,
                traza.estado(
                    "comparar",
                    Some(inicio),
                    &posiciones,
                    Some(&perfil_completo),
                    json!({ "puntos": puntos, "mejoro": mejoro }),
                ),
            ));
        }
    }
    let finales = motivos::motivos_en(textos, &traza.mejores, k);
    let perfil_final = motivos::perfil_de_motivos(&finales, 1.0);
    pasos.push(Paso::nuevo(
        "Ok(MotivosEncontrados {",
        format!(
            "Fin: los mejores motivos tienen puntuación {} y consenso {}.",
            traza.mejor_puntuacion,
            texto(&motivos::consenso(&finales))
        ),
        traza.estado(
            "final",
            None,
            &traza.mejores.clone(),
            Some(&perfil_final),
            json!({ "puntos": traza.mejor_puntuacion }),
        ),
    ));
    Ok(pasos)
}

const COLUMNAS_VISIBLES: usize = 24;
const ELEMENTOS_VISIBLES: usize = 20;

fn recortar(secuencia: &[u8], maximo: usize) -> String {
    if secuencia.len() <= maximo {
        texto(secuencia)
    } else {
        format!("{}…", texto(&secuencia[..maximo]))
    }
}

fn fila_de_letra(base: u8) -> usize {
    LETRAS
        .iter()
        .position(|letra| letra.as_bytes()[0] == base)
        .unwrap_or(0)
}

fn con_tabla_numerada(
    pizarra: Pizarra,
    titulo: &str,
    filas: Vec<Vec<String>>,
    resaltar: &[(usize, usize)],
) -> Pizarra {
    let columnas = filas.first().map_or(0, Vec::len);
    let mut encabezados = vec![String::new()];
    encabezados.extend((1..columnas).map(|columna| columna.to_string()));
    let vistos: Vec<&str> = encabezados.iter().map(String::as_str).collect();
    let visibles: Vec<(usize, usize)> = resaltar
        .iter()
        .copied()
        .filter(|&(fila, columna)| fila < filas.len() && columna < columnas)
        .collect();
    pizarra.tabla(titulo, &vistos, filas, &visibles)
}

fn filas_de_perfil(perfil: &[Vec<f64>]) -> Vec<Vec<String>> {
    LETRAS
        .iter()
        .zip(perfil)
        .map(|(letra, fila)| {
            let mut celdas = vec![letra.to_string()];
            celdas.extend(fila.iter().take(COLUMNAS_VISIBLES).map(|&v| decimal(v, 2)));
            celdas
        })
        .collect()
}

fn filas_de_conteos(conteos: &[Vec<usize>]) -> Vec<Vec<String>> {
    LETRAS
        .iter()
        .zip(conteos)
        .map(|(letra, fila)| {
            let mut celdas = vec![letra.to_string()];
            celdas.extend(
                fila.iter()
                    .take(COLUMNAS_VISIBLES)
                    .map(|conteo| conteo.to_string()),
            );
            celdas
        })
        .collect()
}

fn lista_de_kmers(
    pizarra: Pizarra,
    titulo: &str,
    textos: &[Vec<u8>],
    posiciones: &[usize],
    k: usize,
    resaltar: &[usize],
) -> Pizarra {
    let visibles = textos.len().min(ELEMENTOS_VISIBLES);
    let elementos = kmers_en(&textos[..visibles], posiciones, k);
    let resaltados: Vec<usize> = resaltar
        .iter()
        .copied()
        .filter(|&indice| indice < elementos.len())
        .collect();
    pizarra.lista(titulo, elementos, &resaltados)
}

fn tabla_de_intentos(
    pizarra: Pizarra,
    titulo: &str,
    columna: &str,
    puntos_por_intento: &[usize],
) -> Pizarra {
    let mut cuentas: BTreeMap<usize, usize> = BTreeMap::new();
    for &puntos in puntos_por_intento {
        *cuentas.entry(puntos).or_insert(0) += 1;
    }
    let filas: Vec<Vec<String>> = cuentas
        .iter()
        .take(8)
        .map(|(puntos, veces)| vec![puntos.to_string(), miles(*veces)])
        .collect();
    pizarra.tabla(titulo, &["Puntuación", columna], filas, &[(0, 0), (0, 1)])
}

fn vecinos_de(numero: usize, k: usize, d: usize, aparte: &mut [bool]) -> Vec<usize> {
    motivos::marcar_vecinos(numero, k, 0, d, aparte);
    let mut encontrados = Vec::new();
    for (indice, marca) in aparte.iter_mut().enumerate() {
        if *marca {
            encontrados.push(indice);
            *marca = false;
        }
    }
    encontrados
}

fn contar_verdaderos(marcas: &[bool]) -> usize {
    marcas.iter().filter(|&&marca| marca).count()
}

fn narrar_enumeracion(conjunto: &Conjunto, k: usize, d: usize) -> (Vec<Paso>, Vec<Vec<u8>>) {
    let adn = &conjunto.textos;
    let t = adn.len();
    let cantidad_de_patrones = 1usize << (2 * k);
    let resumen_adn = format!("{} secuencias", miles(t));
    let toleradas = contar(d, "diferencia", "diferencias");
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        Pizarra::nueva()
            .variable("adn", &resumen_adn)
            .variable("k", k)
            .variable("d", d)
            .secuencia(&adn[0], &[(0, k, MarcaDeSecuencia::Actual)])
            .lista(
                "Secuencias",
                conjunto
                    .nombres
                    .iter()
                    .take(ELEMENTOS_VISIBLES)
                    .cloned()
                    .collect(),
                &[],
            )
            .paso(
                "pub fn enumerar_motivos",
                format!(
                    "Buscamos los patrones de {k} letras que aparecen en TODAS las {} secuencias, aunque sea con {toleradas} como mucho. Una secuencia «contiene» un patrón si alguna de sus ventanas de {k} letras se le parece así.",
                    miles(t)
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("adn", &resumen_adn)
            .variable("k", k)
            .variable("d", d)
            .variable("cantidad_de_patrones", miles(cantidad_de_patrones))
            .paso(
                "let cantidad_de_patrones = 1 << (2 * k);",
                format!(
                    "Cada letra tiene 4 opciones, así que hay 4^{k} = {} patrones de {k} letras. Cada patrón tendrá su propia casilla en las listas de verdadero/falso que siguen.",
                    miles(cantidad_de_patrones)
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", k)
            .variable("d", d)
            .variable("cantidad_de_patrones", miles(cantidad_de_patrones))
            .variable(
                "en_todas",
                format!(
                    "{} true de {}",
                    miles(cantidad_de_patrones),
                    miles(cantidad_de_patrones)
                ),
            )
            .paso(
                "let mut en_todas = vec![true; cantidad_de_patrones];",
                "Al principio todos los patrones son candidatos: toda casilla de en_todas empieza en true. Iremos tachando (poniendo false) los que falten en alguna secuencia.",
            ),
    );
    let mut en_todas = vec![true; cantidad_de_patrones];
    let mut en_esta = vec![false; cantidad_de_patrones];
    let mut aparte = vec![false; cantidad_de_patrones];
    let mut historial: Vec<Vec<String>> = Vec::new();
    let salto = (t / 12).max(1);
    for (indice, texto_actual) in adn.iter().enumerate() {
        let mostrar = indice < 6 || (indice + 1) % salto == 0 || indice + 1 == t;
        en_esta.fill(false);
        if indice < 2 {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("k", k)
                    .variable("d", d)
                    .variable("en_esta", format!("0 true de {}", miles(cantidad_de_patrones)))
                    .secuencia(texto_actual, &[(0, k, MarcaDeSecuencia::Actual)])
                    .paso(
                        "en_esta.fill(false);",
                        format!(
                            "Secuencia {} ({}): empezamos con un marcador en blanco, en_esta, donde todo vale false. Aquí anotaremos qué patrones aparecen en esta secuencia.",
                            indice + 1,
                            conjunto.nombres[indice]
                        ),
                    ),
            );
        }
        for (posicion, ventana) in texto_actual.windows(k).enumerate() {
            if let Some(numero) = motivos::patron_a_numero(ventana) {
                if indice == 0 && posicion < 2 && narracion.hay_lugar() {
                    narracion.agregar(
                        Pizarra::nueva()
                            .variable("k", k)
                            .variable("d", d)
                            .variable("numero", miles(numero))
                            .secuencia(
                                texto_actual,
                                &[(posicion, posicion + k, MarcaDeSecuencia::Actual)],
                            )
                            .paso(
                                "if let Some(numero) = patron_a_numero(ventana) {",
                                format!(
                                    "La ventana {} (posición {posicion}) se convierte en el número {}: A=0, C=1, G=2, T=3, leído como un número en base 4. Así cada patrón sabe cuál es su casilla.",
                                    texto(ventana),
                                    miles(numero)
                                ),
                            ),
                    );
                    let vecinos = vecinos_de(numero, k, d, &mut aparte);
                    let propio = vecinos
                        .iter()
                        .take(ELEMENTOS_VISIBLES)
                        .position(|&vecino| vecino == numero);
                    let elementos: Vec<String> = vecinos
                        .iter()
                        .take(ELEMENTOS_VISIBLES)
                        .map(|&vecino| texto(&motivos::numero_a_patron(vecino, k)))
                        .collect();
                    let resaltados: Vec<usize> = propio.into_iter().collect();
                    narracion.agregar(
                        Pizarra::nueva()
                            .variable("k", k)
                            .variable("d", d)
                            .variable("numero", miles(numero))
                            .variable("en_esta", format!("{} true", miles(vecinos.len())))
                            .secuencia(
                                texto_actual,
                                &[(posicion, posicion + k, MarcaDeSecuencia::Actual)],
                            )
                            .lista("Casillas marcadas (primeras)", elementos, &resaltados)
                            .paso(
                                "marcar_vecinos(numero, k, 0, d, &mut en_esta);",
                                format!(
                                    "Ponemos true en la casilla de {} y en la de sus {} vecinos: los patrones que se le parecen con {toleradas} como máximo. En total {} casillas; la lista muestra las primeras, en orden alfabético.",
                                    texto(ventana),
                                    miles(vecinos.len() - 1),
                                    miles(vecinos.len())
                                ),
                            ),
                    );
                }
                motivos::marcar_vecinos(numero, k, 0, d, &mut en_esta);
            }
        }
        if mostrar && narracion.hay_lugar() {
            let marcados = contar_verdaderos(&en_esta);
            let antes = contar_verdaderos(&en_todas);
            for (todas, esta) in en_todas.iter_mut().zip(&en_esta) {
                *todas = *todas && *esta;
            }
            let vivos = contar_verdaderos(&en_todas);
            historial.push(vec![
                format!("Secuencia {}", indice + 1),
                miles(marcados),
                miles(vivos),
            ]);
            let ultima = historial.len() - 1;
            narracion.agregar(
                Pizarra::nueva()
                    .variable("k", k)
                    .variable("d", d)
                    .variable(
                        "en_esta",
                        format!("{} true de {}", miles(marcados), miles(cantidad_de_patrones)),
                    )
                    .variable(
                        "en_todas",
                        format!("{} true de {}", miles(vivos), miles(cantidad_de_patrones)),
                    )
                    .tabla(
                        "Cómo se encoge en_todas",
                        &["Secuencia", "Patrones en en_esta", "Quedan en en_todas"],
                        historial.clone(),
                        &[(ultima, 2)],
                    )
                    .paso(
                        "*todas = *todas && *esta;",
                        format!(
                            "Secuencia {}: {} patrones aparecen en ella (con {toleradas} como máximo). Un patrón sigue vivo solo si estaba en true en las dos listas: de {} candidatos quedan {}.",
                            indice + 1,
                            miles(marcados),
                            miles(antes),
                            miles(vivos)
                        ),
                    ),
            );
        } else {
            for (todas, esta) in en_todas.iter_mut().zip(&en_esta) {
                *todas = *todas && *esta;
            }
        }
    }
    let supervivientes = contar_verdaderos(&en_todas);
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", k)
            .variable("d", d)
            .variable(
                "en_todas",
                format!(
                    "{} true de {}",
                    miles(supervivientes),
                    miles(cantidad_de_patrones)
                ),
            )
            .paso(
                "for (numero, &presente) in en_todas.iter().enumerate() {",
                format!(
                    "Ya pasaron las {} secuencias. Recorremos en_todas: cada casilla que quedó en true es un motivo ({k}, {d}). Quedaron {}.",
                    miles(t),
                    miles(supervivientes)
                ),
            ),
    );
    let mut motivos_hallados = Vec::new();
    for (numero, &presente) in en_todas.iter().enumerate() {
        if presente {
            motivos_hallados.push(motivos::numero_a_patron(numero, k));
        }
    }
    if let Some(primero) = motivos_hallados.first() {
        let elementos: Vec<String> = motivos_hallados
            .iter()
            .take(ELEMENTOS_VISIBLES)
            .map(|patron| texto(patron))
            .collect();
        narracion.agregar(
            Pizarra::nueva()
                .variable("k", k)
                .variable("d", d)
                .lista("Motivos (primeros)", elementos, &[0])
                .paso(
                    "motivos.push(numero_a_patron(numero, k));",
                    format!(
                        "Cada casilla en true se vuelve a escribir como letras: el número {} es el patrón {}. Los vamos juntando en una lista, en orden alfabético.",
                        miles(motivos::patron_a_numero(primero).unwrap_or(0)),
                        texto(primero)
                    ),
                ),
        );
    }
    let nota_final = if motivos_hallados.is_empty() {
        format!(
            "Ningún patrón de {k} letras aparece en todas las secuencias con {toleradas} como máximo. Mira «Qué significa» para ver qué probar."
        )
    } else {
        format!(
            "Resultado: {} ({k}, {d}). Mira «Qué significa» para saber si es más de lo que daría el puro azar.",
            contar(motivos_hallados.len(), "motivo", "motivos")
        )
    };
    let ultimo = Pizarra::nueva()
        .variable("k", k)
        .variable("d", d)
        .variable("motivos", miles(motivos_hallados.len()))
        .paso("Ok(motivos)", nota_final);
    (narracion.terminar(ultimo), motivos_hallados)
}

fn pasos_enumeracion(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (conjunto, k, d) = leer_enumeracion(entrada)?;
    Ok(narrar_enumeracion(&conjunto, k, d).0)
}

fn narrar_perfil(vistas: &[&[u8]]) -> (Vec<Paso>, AnalisisMotivos) {
    let t = vistas.len();
    let k = vistas[0].len();
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let resumen = format!("{} motivos", miles(t));
    let todos_los_motivos: Vec<String> = vistas
        .iter()
        .take(ELEMENTOS_VISIBLES)
        .map(|motivo| recortar(motivo, 40))
        .collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("motivos", &resumen)
            .variable("k", k)
            .lista("Motivos", todos_los_motivos.clone(), &[])
            .paso(
                "pub fn analizar_motivos",
                format!(
                    "Tenemos {} motivos de {k} letras, uno debajo del otro. Vamos a leer la matriz columna por columna para ver qué letras se repiten, que es lo que delata un motivo de verdad.",
                    miles(t)
                ),
            ),
    );
    narracion.agregar(con_tabla_numerada(
        Pizarra::nueva()
            .variable("motivos", &resumen)
            .variable("k", k),
        "conteos (todo en cero)",
        filas_de_conteos(&vec![vec![0; k]; 4]),
        &[],
    )
    .paso(
        "let mut conteos = vec![vec![0; k]; 4];",
        "Preparamos una tabla de 4 filas (A, C, G, T) y una columna por posición del motivo. Cada celda cuenta cuántos motivos tienen esa letra en esa columna; empieza en 0.",
    ));
    let mostrados = t.min(5);
    for n in 1..=mostrados {
        let parcial = motivos::contar_bases(&vistas[..n]);
        let motivo = vistas[n - 1];
        let resaltar: Vec<(usize, usize)> = motivo
            .iter()
            .take(COLUMNAS_VISIBLES)
            .enumerate()
            .map(|(columna, &base)| (fila_de_letra(base), columna + 1))
            .collect();
        let nota = if n == 1 {
            format!(
                "Leemos el primer motivo, {}: por cada letra sumamos 1 en su fila y en su columna. Las celdas resaltadas son las que acaban de subir.",
                recortar(motivo, 40)
            )
        } else {
            format!(
                "Leemos el motivo {n} ({}) y sumamos 1 en la celda de cada una de sus letras. Fíjate cómo crecen los números donde los motivos se parecen.",
                recortar(motivo, 40)
            )
        };
        narracion.agregar(
            con_tabla_numerada(
                Pizarra::nueva()
                    .variable("motivos", &resumen)
                    .variable("k", k)
                    .lista("Motivos", todos_los_motivos.clone(), &[n - 1]),
                "conteos",
                filas_de_conteos(&parcial),
                &resaltar,
            )
            .paso("conteos[fila][columna] += 1;", nota),
        );
    }
    let conteos = motivos::contar_bases(vistas);
    if t > mostrados {
        narracion.agregar(con_tabla_numerada(
            Pizarra::nueva()
                .variable("motivos", &resumen)
                .variable("k", k),
            "conteos (todos los motivos)",
            filas_de_conteos(&conteos),
            &[],
        )
        .paso(
            "for motivo in motivos {",
            format!(
                "Los otros {} motivos se suman igual, uno tras otro. Al final cada columna suma {t}: una letra por motivo.",
                miles(t - mostrados)
            ),
        ));
    }
    let perfil = motivos::perfil_de_motivos(vistas, 0.0);
    let fila_mayor_inicial = motivos::fila_mayor(&conteos, 0);
    narracion.agregar(con_tabla_numerada(
        Pizarra::nueva()
            .variable("motivos", &resumen)
            .variable("k", k)
            .variable("pseudoconteo", 0),
        "perfil (probabilidades)",
        filas_de_perfil(&perfil),
        &[(fila_mayor_inicial, 1)],
    )
    .paso(
        ".map(|&conteo| (conteo as f64 + pseudoconteo) / total)",
        format!(
            "Dividimos cada conteo entre {t} (el número de motivos) para obtener probabilidades: en la columna 1, {} tiene {} de {t}, o sea {}. Cada columna suma 1.",
            LETRAS[fila_mayor_inicial],
            conteos[fila_mayor_inicial][0],
            decimal(perfil[fila_mayor_inicial][0], 2)
        ),
    ));
    let letras_del_consenso = motivos::consenso(vistas);
    let mayores: Vec<usize> = (0..k)
        .map(|columna| motivos::fila_mayor(&conteos, columna))
        .collect();
    let mut detalladas: Vec<usize> = vec![0];
    for columna in 1..k {
        if detalladas.len() < 3 && conteos[mayores[columna]][columna] < t {
            detalladas.push(columna);
        }
    }
    for columna in 1..k {
        if detalladas.len() < 3 && !detalladas.contains(&columna) {
            detalladas.push(columna);
        }
    }
    detalladas.sort_unstable();
    for &columna in &detalladas {
        let mayor = mayores[columna];
        let votos = (0..4)
            .map(|fila| format!("{}={}", LETRAS[fila], conteos[fila][columna]))
            .collect::<Vec<_>>()
            .join(", ");
        narracion.agregar(con_tabla_numerada(
            Pizarra::nueva()
                .variable("k", k)
                .variable("columna", columna)
                .variable(
                    "letras",
                    texto(&letras_del_consenso[..=columna]),
                ),
            "conteos",
            filas_de_conteos(&conteos),
            &[(mayor, columna + 1)],
        )
        .paso(
            "letras.push(BASES[fila_mayor(&conteos, columna)]);",
            format!(
                "Columna {} (contando desde 1): {votos}. Gana {} con {} votos y esa es la letra del consenso en esta posición.",
                columna + 1,
                LETRAS[mayor],
                conteos[mayor][columna]
            ),
        ));
    }
    let resaltar_mayores: Vec<(usize, usize)> = mayores
        .iter()
        .enumerate()
        .map(|(columna, &fila)| (fila, columna + 1))
        .collect();
    narracion.agregar(con_tabla_numerada(
        Pizarra::nueva()
            .variable("k", k)
            .variable("letras", recortar(&letras_del_consenso, 40)),
        "conteos (celda ganadora de cada columna)",
        filas_de_conteos(&conteos),
        &resaltar_mayores,
    )
    .paso(
        "letras.push(BASES[fila_mayor(&conteos, columna)]);",
        format!(
            "Repetimos con todas las columnas. El consenso, la letra ganadora de cada una, es {}.",
            recortar(&letras_del_consenso, 40)
        ),
    ));
    let puntos = motivos::puntuacion(vistas);
    for &columna in &detalladas {
        let coinciden = conteos[mayores[columna]][columna];
        let acumulado: usize = (0..=columna)
            .map(|anterior| t - conteos[mayores[anterior]][anterior])
            .sum();
        narracion.agregar(con_tabla_numerada(
            Pizarra::nueva()
                .variable("k", k)
                .variable("columna", columna)
                .variable("diferencias", acumulado),
            "conteos",
            filas_de_conteos(&conteos),
            &[(mayores[columna], columna + 1)],
        )
        .paso(
            "diferencias += motivos.len() - conteos[mas_comun][columna];",
            format!(
                "Columna {}: de {t} motivos, {coinciden} tienen la letra del consenso y {} no. {} La puntuación, contando desde la columna 1 hasta esta, va en {acumulado}.",
                columna + 1,
                t - coinciden,
                if coinciden == t {
                    "No suma nada.".to_string()
                } else {
                    format!("Esas {} letras suman a la puntuación.", t - coinciden)
                }
            ),
        ));
    }
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", k)
            .variable("diferencias", puntos)
            .paso(
                "diferencias += motivos.len() - conteos[mas_comun][columna];",
                format!(
                    "Sumando las {k} columnas, la puntuación es {}. Cuanto más baja, más se parecen los motivos entre sí; 0 significaría que son idénticos.",
                    miles(puntos)
                ),
            ),
    );
    let entropia = motivos::entropia_por_columna(vistas);
    for &columna in detalladas.iter().take(2) {
        let mut partes = Vec::new();
        for (fila, letra) in LETRAS.iter().enumerate() {
            let conteo = conteos[fila][columna];
            if conteo > 0 {
                let p = conteo as f64 / t as f64;
                partes.push(format!(
                    "{letra}: p={} aporta {}",
                    decimal(p, 2),
                    decimal((-p * p.log2()).abs(), 3)
                ));
            }
        }
        narracion.agregar(
            Pizarra::nueva()
                .variable("k", k)
                .variable("columna", columna)
                .variable("entropia", decimal(entropia[columna], 3))
                .paso(
                    "entropia -= p * p.log2();",
                    format!(
                        "Columna {}: {}. Sumadas dan {} bits. 0 bits es una columna donde todos coinciden; 2 bits, una donde las 4 letras salen por igual.",
                        columna + 1,
                        partes.join("; "),
                        decimal(entropia[columna], 3)
                    ),
                ),
        );
    }
    let entropia_total: f64 = entropia.iter().sum();
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", k)
            .variable("entropia", decimal(entropia_total, 2))
            .paso(
                "entropias.push(entropia);",
                format!(
                    "Hacemos lo mismo con cada columna. La entropía total es {} bits de un máximo de {}: cuanto más baja, más conservado está el motivo.",
                    decimal(entropia_total, 2),
                    2 * k
                ),
            ),
    );
    let analisis = AnalisisMotivos {
        conteos,
        perfil,
        consenso: letras_del_consenso,
        puntos,
        entropia,
    };
    let ultimo = Pizarra::nueva()
        .variable("consenso", recortar(&analisis.consenso, 40))
        .variable("puntos", analisis.puntos)
        .variable("entropia", decimal(entropia_total, 2))
        .paso(
            "Ok(AnalisisMotivos {",
            format!(
                "Resultado: consenso {}, puntuación {} y entropía total {} bits. Mira «Qué significa» para saber si es un motivo de verdad.",
                recortar(&analisis.consenso, 40),
                miles(analisis.puntos),
                decimal(entropia_total, 2)
            ),
        );
    (narracion.terminar(ultimo), analisis)
}

fn pasos_perfil(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let conjunto = leer_matriz_de_motivos(entrada)?;
    let vistas: Vec<&[u8]> = conjunto.textos.iter().map(Vec::as_slice).collect();
    Ok(narrar_perfil(&vistas).0)
}

fn tabla_de_comparacion(
    pizarra: Pizarra,
    titulo: &str,
    arriba: (&str, &[u8]),
    abajo: (&str, &[u8]),
) -> Pizarra {
    let columnas = arriba.1.len().min(COLUMNAS_VISIBLES);
    let mut primera = vec![arriba.0.to_string()];
    let mut segunda = vec![abajo.0.to_string()];
    let mut distintas = Vec::new();
    for columna in 0..columnas {
        primera.push((arriba.1[columna] as char).to_string());
        segunda.push((abajo.1[columna] as char).to_string());
        if arriba.1[columna] != abajo.1[columna] {
            distintas.push((1, columna + 1));
        }
    }
    con_tabla_numerada(pizarra, titulo, vec![primera, segunda], &distintas)
}

fn narrar_distancia(patron: &[u8], conjunto: &Conjunto) -> (Vec<Paso>, usize) {
    let adn = &conjunto.textos;
    let k = patron.len();
    let t = adn.len();
    let resumen = format!("{} secuencias", miles(t));
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        Pizarra::nueva()
            .variable("patron", recortar(patron, 40))
            .variable("adn", &resumen)
            .secuencia(&adn[0], &[(0, k, MarcaDeSecuencia::Actual)])
            .paso(
                "pub fn distancia_patron_textos",
                format!(
                    "Medimos qué tan lejos está el patrón {} de {}. En cada secuencia buscamos la ventana de {k} letras que más se le parece y contamos sus diferencias; luego sumamos todas.",
                    recortar(patron, 40),
                    resumen
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("patron", recortar(patron, 40))
            .variable("total", 0)
            .paso(
                "let mut total = 0;",
                "El total empieza en 0 y crecerá con la distancia de cada secuencia.",
            ),
    );
    let mut total = 0;
    let mut historial: Vec<Vec<String>> = Vec::new();
    let salto = (t / 12).max(1);
    for (indice, texto_actual) in adn.iter().enumerate() {
        let detallada = indice == 0;
        let mut mejor_posicion = 0;
        let mut menor_distancia = k + 1;
        let mut mejoras_mostradas = 0;
        for (posicion, ventana) in texto_actual.windows(k).enumerate() {
            let distancia = motivos::distancia_hamming(patron, ventana);
            if detallada && posicion < 3 && narracion.hay_lugar() {
                let mut marcas = vec![(posicion, posicion + k, MarcaDeSecuencia::Actual)];
                if menor_distancia <= k {
                    marcas.push((
                        mejor_posicion,
                        mejor_posicion + k,
                        MarcaDeSecuencia::Hallado,
                    ));
                }
                narracion.agregar(tabla_de_comparacion(
                    Pizarra::nueva()
                        .variable("patron", recortar(patron, 40))
                        .variable("total", total)
                        .variable("distancia", distancia)
                        .secuencia(texto_actual, &marcas),
                    "Letras que difieren (resaltadas)",
                    ("patrón", patron),
                    ("ventana", ventana),
                )
                .paso(
                    "let distancia = distancia_hamming(patron, ventana);",
                    format!(
                        "Secuencia {}, ventana {posicion}: {}. Comparamos letra por letra con el patrón y hay {distancia} {} (las celdas resaltadas).",
                        indice + 1,
                        recortar(ventana, 40),
                        if distancia == 1 { "diferencia" } else { "diferencias" }
                    ),
                ));
            }
            if distancia < menor_distancia {
                menor_distancia = distancia;
                mejor_posicion = posicion;
                if detallada && posicion >= 3 && mejoras_mostradas < 4 && narracion.hay_lugar() {
                    mejoras_mostradas += 1;
                    narracion.agregar(tabla_de_comparacion(
                        Pizarra::nueva()
                            .variable("patron", recortar(patron, 40))
                            .variable("total", total)
                            .variable("menor_distancia", menor_distancia)
                            .variable("mejor_posicion", mejor_posicion)
                            .secuencia(
                                texto_actual,
                                &[(posicion, posicion + k, MarcaDeSecuencia::Hallado)],
                            ),
                        "Letras que difieren (resaltadas)",
                        ("patrón", patron),
                        ("ventana", ventana),
                    )
                    .paso(
                        "menor_distancia = distancia;",
                        format!(
                            "La ventana de la posición {posicion} ({}) solo difiere en {distancia}: es la más parecida hasta ahora, así que la guardamos.",
                            recortar(ventana, 40)
                        ),
                    ));
                }
                if distancia == 0 {
                    break;
                }
            }
        }
        total += menor_distancia;
        historial.push(vec![
            format!("Secuencia {}", indice + 1),
            miles(mejor_posicion),
            menor_distancia.to_string(),
        ]);
        if (indice < 8 || (indice + 1) % salto == 0) && narracion.hay_lugar() {
            let ultima = historial.len() - 1;
            let visibles = historial.len().saturating_sub(8);
            narracion.agregar(
                Pizarra::nueva()
                    .variable("patron", recortar(patron, 40))
                    .variable("total", total)
                    .secuencia(
                        texto_actual,
                        &[(
                            mejor_posicion,
                            mejor_posicion + k,
                            MarcaDeSecuencia::Hallado,
                        )],
                    )
                    .tabla(
                        "Distancia por secuencia",
                        &["Secuencia", "Mejor posición", "Diferencias"],
                        historial[visibles..].to_vec(),
                        &[(ultima - visibles, 2)],
                    )
                    .paso(
                        "total += distancia_patron_texto(patron, texto);",
                        format!(
                            "Secuencia {}: la ventana más parecida está en la posición {} y difiere en {menor_distancia}. {} Total acumulado: {}.",
                            indice + 1,
                            miles(mejor_posicion),
                            if menor_distancia == 0 {
                                "Es una copia exacta."
                            } else {
                                "Se suma al total."
                            },
                            miles(total)
                        ),
                    ),
            );
        }
    }
    let exactas = historial.iter().filter(|fila| fila[2] == "0").count();
    let ultimo = Pizarra::nueva()
        .variable("patron", recortar(patron, 40))
        .variable("total", total)
        .paso(
            "Ok(total)",
            format!(
                "Resultado: la distancia entre el patrón y las {} secuencias es {}; {} de ellas tienen una copia exacta. Mira «Qué significa» para saber si es poco o mucho.",
                miles(t),
                miles(total),
                miles(exactas)
            ),
        );
    (narracion.terminar(ultimo), total)
}

fn pasos_distancia(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (patron, conjunto) = leer_distancia(entrada)?;
    Ok(narrar_distancia(&patron, &conjunto).0)
}

fn niveles_de_distancia(distancias: &[u8], hasta: usize) -> Vec<usize> {
    let mut cuentas = vec![0usize; hasta + 1];
    for &distancia in distancias {
        if let Some(cuenta) = cuentas.get_mut(distancia as usize) {
            *cuenta += 1;
        }
    }
    cuentas
}

fn lista_de_niveles(cuentas: &[usize], hasta: usize) -> Vec<String> {
    cuentas
        .iter()
        .take(hasta + 1)
        .enumerate()
        .map(|(nivel, cuenta)| format!("distancia {nivel}: {} patrones", miles(*cuenta)))
        .collect()
}

fn narrar_mediana(conjunto: &Conjunto, k: usize) -> Result<(Vec<Paso>, Mediana), String> {
    let adn = &conjunto.textos;
    let t = adn.len();
    let cantidad_de_patrones = 1usize << (2 * k);
    let resumen = format!("{} secuencias", miles(t));
    let patrones_texto = miles(cantidad_de_patrones);
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        Pizarra::nueva()
            .variable("adn", &resumen)
            .variable("k", k)
            .secuencia(&adn[0], &[(0, k, MarcaDeSecuencia::Actual)])
            .paso(
                "pub fn cadena_mediana",
                format!(
                    "Buscamos el patrón de {k} letras que, sumando la distancia a cada una de las {resumen}, queda lo más cerca posible de todas. No adivinamos: probamos los {patrones_texto} patrones posibles."
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", k)
            .variable("cantidad_de_patrones", &patrones_texto)
            .paso(
                "let cantidad_de_patrones = 1 << (2 * k);",
                format!(
                    "Con {k} letras y 4 opciones por letra hay 4^{k} = {patrones_texto} patrones. Cada uno tendrá su casilla numerada."
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", k)
            .variable("cantidad_de_patrones", &patrones_texto)
            .variable("suma_de_distancias", format!("{patrones_texto} ceros"))
            .paso(
                "let mut suma_de_distancias = vec![0usize; cantidad_de_patrones];",
                "En suma_de_distancias, cada patrón empieza con 0. Por cada secuencia le sumaremos qué tan lejos está de ella.",
            ),
    );
    let mut suma_de_distancias = vec![0usize; cantidad_de_patrones];
    let mut distancias = vec![u8::MAX; cantidad_de_patrones];
    let salto = (t / 12).max(1);
    for (indice, texto_actual) in adn.iter().enumerate() {
        if !texto_actual
            .windows(k)
            .any(|ventana| motivos::patron_a_numero(ventana).is_some())
        {
            return Err(format!(
                "La secuencia {} no tiene ningún k-mer formado solo por A, C, G y T.",
                indice + 1
            ));
        }
        motivos::distancias_a_un_texto(texto_actual, k, &mut distancias);
        let mostrar = indice < 3 || (indice + 1) % salto == 0 || indice + 1 == t;
        if indice == 0 && narracion.hay_lugar() {
            let mayor = distancias
                .iter()
                .copied()
                .filter(|&distancia| distancia != u8::MAX)
                .max()
                .unwrap_or(0) as usize;
            let cuentas = niveles_de_distancia(&distancias, mayor);
            narracion.agregar(
                Pizarra::nueva()
                    .variable("k", k)
                    .variable("distancias", format!("{patrones_texto} sin visitar"))
                    .secuencia(texto_actual, &[(0, k, MarcaDeSecuencia::Actual)])
                    .paso(
                        "distancias.fill(SIN_VISITAR);",
                        format!(
                            "Secuencia 1: en distancias anotamos, para cada patrón, qué tan lejos está de esta secuencia. Al principio no sabemos nada de ninguno ({patrones_texto} sin visitar)."
                        ),
                    ),
            );
            narracion.agregar(
                Pizarra::nueva()
                    .variable("k", k)
                    .variable("distancia", 0)
                    .lista("Patrones por distancia", lista_de_niveles(&cuentas, 0), &[0])
                    .secuencia(texto_actual, &[(0, k, MarcaDeSecuencia::Hallado)])
                    .paso(
                        "distancias[numero] = 0;",
                        format!(
                            "Los patrones que aparecen tal cual en la secuencia están a distancia 0: son {} distintos (cada ventana de {k} letras cuenta una vez).",
                            miles(cuentas[0])
                        ),
                    ),
            );
            for nivel in 1..=mayor {
                if !narracion.hay_lugar() {
                    break;
                }
                narracion.agregar(
                    Pizarra::nueva()
                        .variable("k", k)
                        .variable("distancia", nivel)
                        .lista(
                            "Patrones por distancia",
                            lista_de_niveles(&cuentas, nivel),
                            &[nivel],
                        )
                        .paso(
                            "distancia += 1;",
                            format!(
                                "Todo patrón que se obtiene cambiando una sola letra de uno ya visitado queda a distancia {nivel}: {} patrones nuevos. Así la distancia se va propagando como una mancha de aceite.",
                                miles(cuentas[nivel])
                            ),
                        ),
                );
            }
        }
        for (suma, &distancia) in suma_de_distancias.iter_mut().zip(&distancias) {
            *suma += distancia as usize;
        }
        if mostrar && narracion.hay_lugar() {
            let mejor_hasta_ahora = suma_de_distancias.iter().copied().min().unwrap_or(0);
            narracion.agregar(
                Pizarra::nueva()
                    .variable("k", k)
                    .variable(
                        "suma_de_distancias",
                        format!("mínimo {}", miles(mejor_hasta_ahora)),
                    )
                    .secuencia(texto_actual, &[(0, k, MarcaDeSecuencia::Actual)])
                    .paso(
                        "*suma += distancia as usize;",
                        format!(
                            "Secuencia {} lista: sumamos sus distancias a la cuenta de cada patrón. Con {} vistas, el mejor patrón suma {} en total.",
                            indice + 1,
                            contar(indice + 1, "secuencia", "secuencias"),
                            miles(mejor_hasta_ahora)
                        ),
                    ),
            );
        }
    }
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", k)
            .variable("suma_de_distancias", patrones_texto.clone())
            .paso(
                "for (numero, &distancia) in suma_de_distancias.iter().enumerate() {",
                format!(
                    "Ya tenemos, para cada uno de los {patrones_texto} patrones, su distancia total a las {} secuencias. Ahora los recorremos buscando la suma más pequeña.",
                    miles(t)
                ),
            ),
    );
    let mut mediana = Mediana {
        distancia: usize::MAX,
        patrones: Vec::new(),
        empatados: 0,
    };
    let mut minimos: Vec<(usize, usize)> = Vec::new();
    for (numero, &distancia) in suma_de_distancias.iter().enumerate() {
        if distancia < mediana.distancia {
            mediana.distancia = distancia;
            mediana.patrones.clear();
            mediana.empatados = 0;
            minimos.push((numero, distancia));
        }
        if distancia == mediana.distancia {
            mediana.empatados += 1;
            if mediana.patrones.len() < motivos::EMPATES_GUARDADOS {
                mediana.patrones.push(motivos::numero_a_patron(numero, k));
            }
        }
    }
    let cantidad_de_minimos = minimos.len();
    for (orden, &(numero, distancia)) in minimos.iter().enumerate() {
        let al_inicio = orden < 2;
        let al_final = orden + 2 >= cantidad_de_minimos;
        if !(al_inicio || al_final) || !narracion.hay_lugar() {
            continue;
        }
        let aviso = if orden == 2 || (orden + 2 == cantidad_de_minimos && !al_inicio) {
            format!(
                " (hubo {} nuevos mínimos en total; mostramos los primeros y los últimos)",
                miles(cantidad_de_minimos)
            )
        } else {
            String::new()
        };
        narracion.agregar(
            Pizarra::nueva()
                .variable("k", k)
                .variable("numero", miles(numero))
                .variable("distancia", miles(distancia))
                .variable("mediana", miles(distancia))
                .paso(
                    "mediana.distancia = distancia;",
                    format!(
                        "El patrón {} (casilla {}) suma {}: es el menor visto hasta ahora, así que pasa a ser el mejor{aviso}.",
                        texto(&motivos::numero_a_patron(numero, k)),
                        miles(numero),
                        miles(distancia)
                    ),
                ),
        );
    }
    let elegidos: Vec<String> = mediana
        .patrones
        .iter()
        .take(ELEMENTOS_VISIBLES)
        .map(|patron| texto(patron))
        .collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", k)
            .variable("mediana", format!("{} patrones empatados", miles(mediana.empatados)))
            .lista("Medianas (alfabético)", elegidos, &[0])
            .paso(
                "mediana.patrones.push(numero_a_patron(numero, k));",
                format!(
                    "Recorridas todas las casillas, el menor total es {}. {} patrón(es) lo alcanzan; guardamos los primeros en orden alfabético.",
                    miles(mediana.distancia),
                    miles(mediana.empatados)
                ),
            ),
    );
    let nota = format!(
        "Resultado: la cadena mediana es {} con distancia total {}.{} Mira «Qué significa» para ver qué tan cerca queda de cada secuencia.",
        texto(&mediana.patrones[0]),
        miles(mediana.distancia),
        if mediana.empatados > 1 {
            format!(
                " Hay {} patrones empatados y mostramos el primero.",
                miles(mediana.empatados)
            )
        } else {
            String::new()
        }
    );
    let ultimo = Pizarra::nueva()
        .variable("k", k)
        .variable("mediana", texto(&mediana.patrones[0]))
        .paso("Ok(mediana)", nota);
    Ok((narracion.terminar(ultimo), mediana))
}

fn pasos_mediana(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (conjunto, k) = leer_mediana(entrada)?;
    Ok(narrar_mediana(&conjunto, k)?.0)
}

fn narrar_mas_probable(adn: &[u8], perfil: &[Vec<f64>]) -> (Vec<Paso>, usize) {
    let k = perfil[0].len();
    let ventanas_totales = adn.len() - k + 1;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(con_tabla_numerada(
        Pizarra::nueva()
            .variable("texto", format!("{} pb", miles(adn.len())))
            .secuencia(adn, &[(0, k, MarcaDeSecuencia::Actual)]),
        "perfil (una columna por posición del motivo)",
        filas_de_perfil(perfil),
        &[],
    )
    .paso(
        "pub fn posicion_mas_probable",
        format!(
            "Buscamos en el texto la ventana de {k} letras que el perfil considera más probable. El perfil dice, para cada columna, qué tan probable es cada letra en esa posición del motivo."
        ),
    ));
    narracion.agregar(
        Pizarra::nueva()
            .variable("texto", format!("{} pb", miles(adn.len())))
            .variable("k", k)
            .paso(
                "let Some(k) = columnas_del_perfil(perfil) else {",
                format!(
                    "El perfil tiene 4 filas (A, C, G, T) y {k} columnas, así que buscamos ventanas de {k} letras. En un texto de {} letras hay {} ventanas.",
                    miles(adn.len()),
                    miles(ventanas_totales)
                ),
            ),
    );
    let mut mejor_posicion = 0;
    let mut mayor_probabilidad = -1.0;
    let mut ceros = 0usize;
    let mut mejoras = 0usize;
    for (posicion, ventana) in adn.windows(k).enumerate() {
        let probabilidad = motivos::probabilidad_kmer(ventana, perfil);
        if probabilidad == 0.0 {
            ceros += 1;
        }
        if posicion < 3 && narracion.hay_lugar() {
            let factores: Vec<f64> = ventana
                .iter()
                .enumerate()
                .map(|(columna, &base)| perfil[fila_de_letra(base)][columna])
                .collect();
            let mut celdas = vec!["probabilidad".to_string()];
            celdas.extend(
                factores
                    .iter()
                    .take(COLUMNAS_VISIBLES)
                    .map(|&factor| decimal(factor, 2)),
            );
            let mut letras = vec!["ventana".to_string()];
            letras.extend(
                ventana
                    .iter()
                    .take(COLUMNAS_VISIBLES)
                    .map(|&base| (base as char).to_string()),
            );
            let ceros_visibles: Vec<(usize, usize)> = factores
                .iter()
                .take(COLUMNAS_VISIBLES)
                .enumerate()
                .filter(|(_, factor)| **factor == 0.0)
                .map(|(columna, _)| (1, columna + 1))
                .collect();
            let nota = if probabilidad == 0.0 {
                format!(
                    "Ventana {posicion} ({}): multiplicamos la probabilidad de cada letra en su columna. Hay una letra con probabilidad 0 en el perfil, y un 0 anula todo el producto: la ventana vale {}.",
                    recortar(ventana, 40),
                    cientifico(probabilidad)
                )
            } else {
                format!(
                    "Ventana {posicion} ({}): miramos la probabilidad de cada letra en su columna y las multiplicamos. Sale {}.",
                    recortar(ventana, 40),
                    cientifico(probabilidad)
                )
            };
            narracion.agregar(
                con_tabla_numerada(
                    Pizarra::nueva()
                        .variable("k", k)
                        .variable("mejor_posicion", mejor_posicion)
                        .variable("mayor_probabilidad", cientifico(mayor_probabilidad))
                        .variable("probabilidad", cientifico(probabilidad))
                        .secuencia(adn, &[(posicion, posicion + k, MarcaDeSecuencia::Actual)]),
                    "Probabilidad de cada letra de la ventana",
                    vec![letras, celdas],
                    &ceros_visibles,
                )
                .paso(
                    "let probabilidad = probabilidad_kmer(ventana, perfil);",
                    nota,
                ),
            );
        }
        if probabilidad > mayor_probabilidad {
            let anterior = mayor_probabilidad;
            mayor_probabilidad = probabilidad;
            mejor_posicion = posicion;
            mejoras += 1;
            if (posicion < 3 || mejoras <= 8) && narracion.hay_lugar() {
                let nota = if anterior < 0.0 {
                    format!(
                        "Es la primera ventana, así que por ahora es la mejor: {}.",
                        cientifico(probabilidad)
                    )
                } else {
                    format!(
                        "La ventana de la posición {posicion} ({}) vale {}, más que la mejor anterior ({}). Pasa a ser la mejor.",
                        recortar(ventana, 40),
                        cientifico(probabilidad),
                        cientifico(anterior)
                    )
                };
                narracion.agregar(
                    Pizarra::nueva()
                        .variable("k", k)
                        .variable("mejor_posicion", mejor_posicion)
                        .variable("mayor_probabilidad", cientifico(mayor_probabilidad))
                        .variable("probabilidad", cientifico(probabilidad))
                        .secuencia(adn, &[(posicion, posicion + k, MarcaDeSecuencia::Hallado)])
                        .paso("mejor_posicion = posicion;", nota),
                );
            }
        } else if posicion < 3 && narracion.hay_lugar() {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("k", k)
                    .variable("mejor_posicion", mejor_posicion)
                    .variable("mayor_probabilidad", cientifico(mayor_probabilidad))
                    .variable("probabilidad", cientifico(probabilidad))
                    .secuencia(
                        adn,
                        &[
                            (posicion, posicion + k, MarcaDeSecuencia::Actual),
                            (
                                mejor_posicion,
                                mejor_posicion + k,
                                MarcaDeSecuencia::Hallado,
                            ),
                        ],
                    )
                    .paso(
                        "if probabilidad > mayor_probabilidad {",
                        format!(
                            "{} no supera a la mejor ({}), así que la mejor se queda donde estaba.",
                            cientifico(probabilidad),
                            cientifico(mayor_probabilidad)
                        ),
                    ),
            );
        }
    }
    let elegido = &adn[mejor_posicion..mejor_posicion + k];
    let nota_resumen = format!(
        "Revisamos las {} ventanas. {} tenían probabilidad 0 (alguna letra imposible según el perfil) y la mejor cambió {} veces.",
        miles(ventanas_totales),
        miles(ceros),
        miles(mejoras)
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", k)
            .variable("mejor_posicion", mejor_posicion)
            .variable("mayor_probabilidad", cientifico(mayor_probabilidad))
            .secuencia(
                adn,
                &[(
                    mejor_posicion,
                    mejor_posicion + k,
                    MarcaDeSecuencia::Hallado,
                )],
            )
            .paso(
                "for (posicion, ventana) in texto.windows(k).enumerate() {",
                nota_resumen,
            ),
    );
    let ultimo = Pizarra::nueva()
        .variable("mejor_posicion", mejor_posicion)
        .variable("mayor_probabilidad", cientifico(mayor_probabilidad))
        .secuencia(
            adn,
            &[(mejor_posicion, mejor_posicion + k, MarcaDeSecuencia::Hallado)],
        )
        .paso(
            "mejor_posicion = posicion;",
            format!(
                "Resultado: la ventana más probable es {} en la posición {} (contando desde 0), con probabilidad {}. Mira «Qué significa» para ver qué tan buena es.",
                recortar(elegido, 40),
                miles(mejor_posicion),
                cientifico(mayor_probabilidad)
            ),
        );
    (narracion.terminar(ultimo), mejor_posicion)
}

fn pasos_mas_probable(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (adn, _, perfil) = leer_mas_probable(entrada)?;
    Ok(narrar_mas_probable(&adn, &perfil).0)
}

fn nota_de_pseudoconteo(pseudoconteo: f64) -> &'static str {
    if pseudoconteo == 0.0 {
        "Sin pseudoconteos, una letra que no aparece en una columna tiene probabilidad 0."
    } else {
        "Con pseudoconteos, ninguna letra tiene probabilidad 0."
    }
}

fn narrar_codiciosa(
    textos: &[Vec<u8>],
    k: usize,
    pseudoconteo: f64,
) -> (Vec<Paso>, MotivosEncontrados) {
    let t = textos.len();
    let rondas = textos[0].len() - k + 1;
    let resumen = format!("{} secuencias", miles(t));
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let mut mejores = vec![0; t];
    let mut mejor_puntuacion = motivos::puntuacion(&motivos::motivos_en(textos, &mejores, k));
    narracion.agregar(
        Pizarra::nueva()
            .variable("adn", &resumen)
            .variable("k", k)
            .variable("pseudoconteo", pseudoconteo)
            .secuencia(&textos[0], &[(0, k, MarcaDeSecuencia::Actual)])
            .paso(
                "pub fn busqueda_codiciosa(",
                format!(
                    "Buscamos un k-mer de {k} letras en cada una de las {resumen} para que, puestos uno debajo del otro, se parezcan lo más posible. La idea codiciosa: fijar un k-mer de la primera secuencia y completar cada fila con el k-mer que mejor encaje con las filas anteriores. {}",
                    nota_de_pseudoconteo(pseudoconteo)
                ),
            ),
    );
    narracion.agregar(
        lista_de_kmers(
            Pizarra::nueva()
                .variable("k", k)
                .variable("mejores", format!("{} ceros", miles(t))),
            "k-mers de «mejores»",
            textos,
            &mejores,
            k,
            &[],
        )
        .paso(
            "let mut mejores = vec![0; adn.len()];",
            "«mejores» guarda la posición elegida en cada secuencia. Empezamos suponiendo que lo mejor es el primer k-mer (posición 0) de cada una.",
        ),
    );
    narracion.agregar(
        lista_de_kmers(
            Pizarra::nueva()
                .variable("k", k)
                .variable("mejor_puntuacion", mejor_puntuacion),
            "k-mers de «mejores»",
            textos,
            &mejores,
            k,
            &[],
        )
        .paso(
            "let mut mejor_puntuacion = puntuacion(&motivos_en(adn, &mejores, k));",
            format!(
                "La puntuación cuenta, columna por columna, cuántas letras no coinciden con la más repetida de esa columna: 0 sería motivos idénticos. Con esta primera elección es {}; cualquier cosa mejor la reemplazará.",
                miles(mejor_puntuacion)
            ),
        ),
    );
    let mut posiciones = vec![0; t];
    let mut puntos_por_intento = Vec::new();
    let mut mejoras = 0usize;
    for inicio in 0..rondas {
        posiciones[0] = inicio;
        let detallada = inicio == 0;
        if detallada {
            narracion.agregar(
                lista_de_kmers(
                    Pizarra::nueva()
                        .variable("k", k)
                        .variable("mejor_puntuacion", mejor_puntuacion)
                        .variable("inicio", inicio)
                        .secuencia(
                            &textos[0],
                            &[(inicio, inicio + k, MarcaDeSecuencia::Actual)],
                        ),
                    "k-mers elegidos en esta ronda",
                    textos,
                    &posiciones[..1],
                    k,
                    &[0],
                )
                .paso(
                    "posiciones[0] = inicio;",
                    format!(
                        "Ronda {inicio} de {}: fijamos el k-mer de la posición {inicio} (contando desde 0) de la primera secuencia. Con él vamos a completar las demás filas.",
                        miles(rondas)
                    ),
                ),
            );
        }
        for i in 1..t {
            let elegidos = motivos::motivos_en(&textos[..i], &posiciones[..i], k);
            let perfil = motivos::perfil_de_motivos(&elegidos, pseudoconteo);
            posiciones[i] = motivos::posicion_mas_probable(&textos[i], &perfil);
            if detallada && narracion.hay_lugar() {
                if i <= 2 {
                    narracion.agregar(con_tabla_numerada(
                        Pizarra::nueva()
                            .variable("k", k)
                            .variable("pseudoconteo", pseudoconteo)
                            .variable("inicio", inicio)
                            .variable("i", i),
                        "perfil de las filas ya elegidas",
                        filas_de_perfil(&perfil),
                        &[],
                    )
                    .paso(
                        "let perfil = perfil_de_motivos(&elegidos, pseudoconteo);",
                        format!(
                            "Con {} armamos el perfil: en cada columna, qué fracción de ellos tiene cada letra. {}",
                            contar(i, "k-mer ya elegido", "k-mers ya elegidos"),
                            nota_de_pseudoconteo(pseudoconteo)
                        ),
                    ));
                }
                if i <= 5 || i + 1 == t {
                    let elegido = &textos[i][posiciones[i]..posiciones[i] + k];
                    let probabilidad = motivos::probabilidad_kmer(elegido, &perfil);
                    let sin_probabilidad = textos[i]
                        .windows(k)
                        .all(|ventana| motivos::probabilidad_kmer(ventana, &perfil) == 0.0);
                    let salto = if i > 6 && i + 1 == t {
                        format!("Las secuencias 7 a {i} se completan igual. ")
                    } else {
                        String::new()
                    };
                    let nota = if sin_probabilidad {
                        format!(
                            "{salto}En la secuencia {} ninguna ventana tiene probabilidad mayor que 0 con este perfil, así que el código se queda con la primera (posición 0): {}. Es lo que arreglan los pseudoconteos.",
                            i + 1,
                            recortar(elegido, 40)
                        )
                    } else {
                        format!(
                            "{salto}En la secuencia {}, la ventana más probable según ese perfil es {} (posición {}, probabilidad {}). Esa fila queda elegida.",
                            i + 1,
                            recortar(elegido, 40),
                            posiciones[i],
                            cientifico(probabilidad)
                        )
                    };
                    narracion.agregar(
                        lista_de_kmers(
                            Pizarra::nueva()
                                .variable("k", k)
                                .variable("inicio", inicio)
                                .variable("i", i)
                                .secuencia(
                                    &textos[i],
                                    &[(
                                        posiciones[i],
                                        posiciones[i] + k,
                                        MarcaDeSecuencia::Hallado,
                                    )],
                                ),
                            "k-mers elegidos en esta ronda",
                            textos,
                            &posiciones[..=i],
                            k,
                            &[i],
                        )
                        .paso(
                            "posiciones[i] = posicion_mas_probable(&adn[i], &perfil);",
                            nota,
                        ),
                    );
                }
            }
        }
        let actuales = motivos::motivos_en(textos, &posiciones, k);
        let puntos = motivos::puntuacion(&actuales);
        puntos_por_intento.push(puntos);
        let anterior = mejor_puntuacion;
        let mejoro = puntos < mejor_puntuacion;
        if mejoro {
            mejor_puntuacion = puntos;
            mejores.copy_from_slice(&posiciones);
            mejoras += 1;
        }
        let decision = if mejoro {
            format!(
                "Puntuación {puntos}: mejora a la mejor hasta ahora ({}) y la reemplaza.",
                miles(anterior)
            )
        } else {
            format!(
                "Puntuación {puntos}: no mejora a la mejor hasta ahora ({}).",
                miles(anterior)
            )
        };
        if detallada && narracion.hay_lugar() {
            narracion.agregar(
                lista_de_kmers(
                    Pizarra::nueva()
                        .variable("k", k)
                        .variable("mejor_puntuacion", anterior)
                        .variable("inicio", inicio)
                        .variable("puntos", puntos),
                    "k-mers de esta ronda",
                    textos,
                    &posiciones,
                    k,
                    &[],
                )
                .paso(
                    "let puntos = puntuacion(&motivos_en(adn, &posiciones, k));",
                    format!(
                        "Con los {t} k-mers elegidos medimos la puntuación de esta ronda: {}. Menos es mejor.",
                        miles(puntos)
                    ),
                ),
            );
            if narracion.hay_lugar() {
                let linea = if mejoro {
                    "mejores.copy_from_slice(&posiciones);"
                } else {
                    "if puntos < mejor_puntuacion {"
                };
                narracion.agregar(
                    lista_de_kmers(
                        Pizarra::nueva()
                            .variable("k", k)
                            .variable("mejor_puntuacion", mejor_puntuacion)
                            .variable("inicio", inicio)
                            .variable("puntos", puntos),
                        "k-mers de «mejores»",
                        textos,
                        &mejores,
                        k,
                        &[],
                    )
                    .paso(linea, decision),
                );
            }
        } else if inicio < 3 && narracion.hay_lugar() {
            let linea = if mejoro {
                "mejores.copy_from_slice(&posiciones);"
            } else {
                "if puntos < mejor_puntuacion {"
            };
            narracion.agregar(
                lista_de_kmers(
                    Pizarra::nueva()
                        .variable("k", k)
                        .variable("mejor_puntuacion", mejor_puntuacion)
                        .variable("inicio", inicio)
                        .variable("puntos", puntos)
                        .secuencia(
                            &textos[0],
                            &[(inicio, inicio + k, MarcaDeSecuencia::Actual)],
                        ),
                    "k-mers de esta ronda",
                    textos,
                    &posiciones,
                    k,
                    &[],
                )
                .paso(
                    linea,
                    format!("Ronda {inicio}: fijamos el k-mer {} y completamos las demás filas. {decision}", recortar(&textos[0][inicio..inicio + k], 40)),
                ),
            );
        } else if mejoro && mejoras <= 25 && narracion.hay_lugar() {
            narracion.agregar(
                lista_de_kmers(
                    Pizarra::nueva()
                        .variable("k", k)
                        .variable("mejor_puntuacion", mejor_puntuacion)
                        .variable("inicio", inicio)
                        .variable("puntos", puntos)
                        .secuencia(
                            &textos[0],
                            &[(inicio, inicio + k, MarcaDeSecuencia::Hallado)],
                        ),
                    "k-mers de «mejores»",
                    textos,
                    &mejores,
                    k,
                    &[],
                )
                .paso(
                    "mejores.copy_from_slice(&posiciones);",
                    format!("Ronda {inicio}: fijamos el k-mer {} y completamos las demás filas. {decision}", recortar(&textos[0][inicio..inicio + k], 40)),
                ),
            );
        }
    }
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", k)
            .variable("mejor_puntuacion", mejor_puntuacion)
            .paso(
                "for inicio in 0..=adn[0].len() - k {",
                format!(
                    "Probamos las {} posiciones posibles de la primera secuencia. {} rondas mejoraron la mejor puntuación; las demás no la superaron y no las mostramos una por una.",
                    miles(rondas),
                    miles(mejoras)
                ),
            ),
    );
    let consenso = motivos::consenso(&motivos::motivos_en(textos, &mejores, k));
    let ultimo = lista_de_kmers(
        Pizarra::nueva()
            .variable("k", k)
            .variable("mejor_puntuacion", mejor_puntuacion)
            .variable("consenso", recortar(&consenso, 40)),
        "Motivos encontrados",
        textos,
        &mejores,
        k,
        &[],
    )
    .paso(
        "Ok(MotivosEncontrados {",
        format!(
            "Resultado: los mejores motivos tienen puntuación {} y consenso {}. Mira «Qué significa» para saber si es una señal real.",
            miles(mejor_puntuacion),
            recortar(&consenso, 40)
        ),
    );
    (
        narracion.terminar(ultimo),
        MotivosEncontrados {
            posiciones: mejores,
            puntos: mejor_puntuacion,
            puntos_por_intento,
        },
    )
}

fn pasos_codiciosa(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let conjunto = leer_conjunto(entrada, "adn", 2)?;
    let k = leer_k_de_perfil(entrada, &conjunto.textos)?;
    revisar_costo(costo_codicioso(&conjunto.textos, k), PRESUPUESTO)?;
    Ok(narrar_codiciosa(&conjunto.textos, k, 0.0).0)
}

fn narrar_aleatoria(
    textos: &[Vec<u8>],
    k: usize,
    repeticiones: usize,
    semilla: u64,
) -> (Vec<Paso>, MotivosEncontrados) {
    let t = textos.len();
    let resumen = format!("{} secuencias", miles(t));
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        Pizarra::nueva()
            .variable("adn", &resumen)
            .variable("k", k)
            .variable("repeticiones", miles(repeticiones))
            .variable("semilla", semilla)
            .paso(
                "pub fn busqueda_aleatoria(",
                format!(
                    "Buscamos un k-mer de {k} letras en cada una de las {resumen}. Receta: empezar con posiciones al azar, mejorarlas con el perfil hasta que no haya mejora, y repetir {} veces quedándonos con el mejor resultado. La semilla {semilla} hace que el «azar» sea siempre el mismo.",
                    miles(repeticiones)
                ),
            ),
    );
    let mut azar = Azar::nuevo(semilla);
    let mut mejores = Vec::new();
    let mut mejor_puntuacion = usize::MAX;
    let mut puntos_por_intento = Vec::with_capacity(repeticiones);
    let mut records = 0usize;
    for repeticion in 0..repeticiones {
        let detallada = repeticion == 0;
        let mut mejores_aqui = motivos::posiciones_al_azar(textos, k, &mut azar);
        let mut mejor_aqui = motivos::puntuacion(&motivos::motivos_en(textos, &mejores_aqui, k));
        let inicial = mejor_aqui;
        if detallada {
            narracion.agregar(
                lista_de_kmers(
                    Pizarra::nueva()
                        .variable("k", k)
                        .variable("semilla", semilla),
                    "k-mers sorteados",
                    textos,
                    &mejores_aqui,
                    k,
                    &[],
                )
                .paso(
                    "let mut mejores = posiciones_al_azar(adn, k, azar);",
                    "Repetición 1: sorteamos una posición al azar en cada secuencia. Es un punto de partida sin ninguna inteligencia.",
                ),
            );
            narracion.agregar(
                lista_de_kmers(
                    Pizarra::nueva()
                        .variable("k", k)
                        .variable("mejor_puntuacion", mejor_aqui),
                    "k-mers sorteados",
                    textos,
                    &mejores_aqui,
                    k,
                    &[],
                )
                .paso(
                    "let mut mejor_puntuacion = puntuacion(&motivos_en(adn, &mejores, k));",
                    format!(
                        "La puntuación de este arranque al azar es {}. Menos es mejor: 0 serían motivos idénticos.",
                        miles(mejor_aqui)
                    ),
                ),
            );
        }
        let mut rondas = 0usize;
        loop {
            let perfil = motivos::perfil_de_motivos(
                &motivos::motivos_en(textos, &mejores_aqui, k),
                motivos::PSEUDOCONTEO_DE_LAPLACE,
            );
            let posiciones = motivos::posiciones_mas_probables(textos, &perfil);
            let puntos = motivos::puntuacion(&motivos::motivos_en(textos, &posiciones, k));
            rondas += 1;
            let mostrar = detallada && rondas <= 6 && narracion.hay_lugar();
            if mostrar && rondas <= 2 {
                narracion.agregar(con_tabla_numerada(
                    Pizarra::nueva()
                        .variable("k", k)
                        .variable("mejor_puntuacion", mejor_aqui),
                    "perfil de los k-mers actuales",
                    filas_de_perfil(&perfil),
                    &[],
                )
                .paso(
                    "let perfil = perfil_de_motivos(&motivos_en(adn, &mejores, k), PSEUDOCONTEO_DE_LAPLACE);",
                    format!(
                        "Ronda de mejora {rondas}: con los k-mers actuales armamos un perfil (sumando 1 a cada conteo para que ninguna letra valga 0). Dice qué letras dominan en cada columna."
                    ),
                ));
                narracion.agregar(
                    lista_de_kmers(
                        Pizarra::nueva()
                            .variable("k", k)
                            .variable("mejor_puntuacion", mejor_aqui),
                        "k-mers más probables",
                        textos,
                        &posiciones,
                        k,
                        &[],
                    )
                    .paso(
                        "let posiciones = posiciones_mas_probables(adn, &perfil);",
                        "En cada secuencia elegimos ahora la ventana más probable según ese perfil: los k-mers se parecen más entre sí.",
                    ),
                );
            }
            if puntos >= mejor_aqui {
                if mostrar && narracion.hay_lugar() {
                    narracion.agregar(
                        Pizarra::nueva()
                            .variable("k", k)
                            .variable("mejor_puntuacion", mejor_aqui)
                            .variable("puntos", puntos)
                            .paso(
                                "return (mejores, mejor_puntuacion);",
                                format!(
                                    "La nueva puntuación ({}) no es mejor que {}: dejamos de mejorar y nos quedamos con los k-mers anteriores. Esta repetición terminó.",
                                    miles(puntos),
                                    miles(mejor_aqui)
                                ),
                            ),
                    );
                }
                break;
            }
            if mostrar && narracion.hay_lugar() {
                narracion.agregar(
                    Pizarra::nueva()
                        .variable("k", k)
                        .variable("mejor_puntuacion", mejor_aqui)
                        .variable("puntos", puntos)
                        .paso(
                            "let puntos = puntuacion(&motivos_en(adn, &posiciones, k));",
                            format!(
                                "La puntuación baja de {} a {}: mejoró, así que adoptamos estos k-mers y volvemos a armar el perfil.",
                                miles(mejor_aqui),
                                miles(puntos)
                            ),
                        ),
                );
            }
            mejores_aqui = posiciones;
            mejor_aqui = puntos;
        }
        puntos_por_intento.push(mejor_aqui);
        let nuevo_record = mejor_aqui < mejor_puntuacion;
        if nuevo_record {
            mejor_puntuacion = mejor_aqui;
            mejores = mejores_aqui;
            records += 1;
        }
        if nuevo_record && repeticion > 0 && records <= 12 && narracion.hay_lugar() {
            narracion.agregar(
                lista_de_kmers(
                    Pizarra::nueva()
                        .variable("k", k)
                        .variable("mejor_puntuacion", mejor_puntuacion)
                        .variable("puntos", mejor_aqui),
                    "mejores k-mers hasta ahora",
                    textos,
                    &mejores,
                    k,
                    &[],
                )
                .paso(
                    "mejor_puntuacion = puntos;",
                    format!(
                        "Repetición {}: partió de {} y mejoró en {} rondas hasta {}. Es el mejor resultado de todas las repeticiones hasta ahora.",
                        miles(repeticion + 1),
                        miles(inicial),
                        rondas - 1,
                        miles(mejor_aqui)
                    ),
                ),
            );
        } else if repeticion > 0 && repeticion < 3 && narracion.hay_lugar() {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("k", k)
                    .variable("mejor_puntuacion", mejor_puntuacion)
                    .variable("puntos", mejor_aqui)
                    .paso(
                        "let (posiciones, puntos) = mejorar_desde_el_azar(adn, k, &mut azar);",
                        format!(
                            "Repetición {}: partió de {} y mejoró en {} rondas hasta {}. No supera a la mejor ({}).",
                            repeticion + 1,
                            miles(inicial),
                            rondas - 1,
                            miles(mejor_aqui),
                            miles(mejor_puntuacion)
                        ),
                    ),
            );
        }
    }
    let veces = puntos_por_intento
        .iter()
        .filter(|&&puntos| puntos == mejor_puntuacion)
        .count();
    narracion.agregar(tabla_de_intentos(
        Pizarra::nueva()
            .variable("k", k)
            .variable("mejor_puntuacion", mejor_puntuacion),
        "Puntuación final de cada repetición (las mejores)",
        "Repeticiones",
        &puntos_por_intento,
    )
    .paso(
        "for _ in 0..repeticiones {",
        format!(
            "Hicimos {} repeticiones. La mejor puntuación, {}, salió en {} de ellas. Que varias lleguen al mismo sitio es buena señal; con otra semilla podrían salir otras.",
            miles(repeticiones),
            miles(mejor_puntuacion),
            miles(veces)
        ),
    ));
    let consenso = motivos::consenso(&motivos::motivos_en(textos, &mejores, k));
    let ultimo = lista_de_kmers(
        Pizarra::nueva()
            .variable("k", k)
            .variable("mejor_puntuacion", mejor_puntuacion)
            .variable("consenso", recortar(&consenso, 40)),
        "Motivos encontrados",
        textos,
        &mejores,
        k,
        &[],
    )
    .paso(
        "Ok(MotivosEncontrados {",
        format!(
            "Resultado: los mejores motivos tienen puntuación {} y consenso {}. Mira «Qué significa» para saber si es una señal real.",
            miles(mejor_puntuacion),
            recortar(&consenso, 40)
        ),
    );
    (
        narracion.terminar(ultimo),
        MotivosEncontrados {
            posiciones: mejores,
            puntos: mejor_puntuacion,
            puntos_por_intento,
        },
    )
}

fn pasos_aleatoria(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (conjunto, k, repeticiones, semilla) = leer_aleatoria(entrada)?;
    Ok(narrar_aleatoria(&conjunto.textos, k, repeticiones, semilla).0)
}

fn narrar_gibbs(
    textos: &[Vec<u8>],
    k: usize,
    iteraciones: usize,
    reinicios: usize,
    semilla: u64,
) -> (Vec<Paso>, MotivosEncontrados) {
    let t = textos.len();
    let resumen = format!("{} secuencias", miles(t));
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        Pizarra::nueva()
            .variable("adn", &resumen)
            .variable("k", k)
            .variable("iteraciones", miles(iteraciones))
            .variable("reinicios", miles(reinicios))
            .variable("semilla", semilla)
            .paso(
                "pub fn muestreo_gibbs(",
                format!(
                    "Buscamos un k-mer de {k} letras en cada una de las {resumen}. Receta: partir de posiciones al azar y, {} veces, quitar un k-mer, armar el perfil con los demás y sortear su reemplazo según ese perfil. Se repite desde cero {} veces (reinicios). La semilla {semilla} fija el azar.",
                    miles(iteraciones),
                    miles(reinicios)
                ),
            ),
    );
    let mut azar = Azar::nuevo(semilla);
    let mut mejores = Vec::new();
    let mut mejor_puntuacion = usize::MAX;
    let mut puntos_por_intento = Vec::with_capacity(reinicios);
    let mut records = 0usize;
    for reinicio in 0..reinicios {
        let detallado = reinicio == 0;
        let mut posiciones = motivos::posiciones_al_azar(textos, k, &mut azar);
        let mut mejores_aqui = posiciones.clone();
        let mut mejor_aqui = motivos::puntuacion(&motivos::motivos_en(textos, &mejores_aqui, k));
        let inicial = mejor_aqui;
        let mut pesos = Vec::new();
        let mut mejoras_aqui = 0usize;
        if detallado {
            narracion.agregar(
                lista_de_kmers(
                    Pizarra::nueva()
                        .variable("k", k)
                        .variable("semilla", semilla),
                    "k-mers sorteados",
                    textos,
                    &posiciones,
                    k,
                    &[],
                )
                .paso(
                    "let mut posiciones = posiciones_al_azar(adn, k, azar);",
                    "Reinicio 1: sorteamos una posición al azar en cada secuencia. Es el punto de partida.",
                ),
            );
            narracion.agregar(
                Pizarra::nueva()
                    .variable("k", k)
                    .variable("mejor_puntuacion", mejor_aqui)
                    .paso(
                        "let mut mejor_puntuacion = puntuacion(&motivos_en(adn, &mejores, k));",
                        format!(
                            "La puntuación del arranque es {}. Menos es mejor: 0 serían motivos idénticos.",
                            miles(mejor_aqui)
                        ),
                    ),
            );
        }
        for iteracion in 0..iteraciones {
            let i = azar.entero(t);
            let mut otros = motivos::motivos_en(textos, &posiciones, k);
            otros.remove(i);
            let perfil = motivos::perfil_de_motivos(&otros, motivos::PSEUDOCONTEO_DE_LAPLACE);
            pesos.clear();
            for ventana in textos[i].windows(k) {
                pesos.push(motivos::probabilidad_kmer(ventana, &perfil));
            }
            let mostrar = detallado && iteracion < 3 && narracion.hay_lugar();
            let completo = mostrar && iteracion == 0;
            let anterior = posiciones[i];
            let previas: Vec<usize> = if completo {
                posiciones.clone()
            } else {
                Vec::new()
            };
            let kmers_otros: Vec<String> = if completo {
                otros
                    .iter()
                    .take(ELEMENTOS_VISIBLES)
                    .map(|m| texto(m))
                    .collect()
            } else {
                Vec::new()
            };
            posiciones[i] = azar.elegir_con_pesos(&pesos);
            let puntos = motivos::puntuacion(&motivos::motivos_en(textos, &posiciones, k));
            if mostrar {
                let numero_de_iteracion = iteracion + 1;
                if completo {
                    narracion.agregar(
                        lista_de_kmers(
                            Pizarra::nueva()
                                .variable("k", k)
                                .variable("i", i)
                                .secuencia(
                                    &textos[i],
                                    &[(anterior, anterior + k, MarcaDeSecuencia::Actual)],
                                ),
                            "k-mers actuales",
                            textos,
                            &previas,
                            k,
                            &[i],
                        )
                        .paso(
                            "let i = azar.entero(adn.len());",
                            format!(
                                "Iteración {numero_de_iteracion}: sorteamos la secuencia {} y la vamos a reconsiderar. Su k-mer actual es {}.",
                                i + 1,
                                recortar(&textos[i][anterior..anterior + k], 40)
                            ),
                        ),
                    );
                    narracion.agregar(
                        Pizarra::nueva()
                            .variable("k", k)
                            .variable("i", i)
                            .lista("k-mers de las demás (otros)", kmers_otros, &[])
                            .paso(
                                "otros.remove(i);",
                                format!(
                                    "Sacamos el k-mer de la secuencia {} y nos quedamos con los de las otras {}. Así el perfil no puede «copiar» lo que ya tenía.",
                                    i + 1,
                                    t - 1
                                ),
                            ),
                    );
                    narracion.agregar(con_tabla_numerada(
                        Pizarra::nueva().variable("k", k).variable("i", i),
                        "perfil de las demás",
                        filas_de_perfil(&perfil),
                        &[],
                    )
                    .paso(
                        "let perfil = perfil_de_motivos(&otros, PSEUDOCONTEO_DE_LAPLACE);",
                        "Con los k-mers de las demás armamos el perfil, sumando 1 a cada conteo para que ninguna letra valga 0.",
                    ));
                    let mut orden: Vec<usize> = (0..pesos.len()).collect();
                    orden.sort_by(|&a, &b| pesos[b].total_cmp(&pesos[a]));
                    let filas: Vec<Vec<String>> = orden
                        .iter()
                        .take(4)
                        .map(|&posicion| {
                            vec![
                                recortar(&textos[i][posicion..posicion + k], 40),
                                posicion.to_string(),
                                cientifico(pesos[posicion]),
                            ]
                        })
                        .collect();
                    narracion.agregar(
                        Pizarra::nueva()
                            .variable("k", k)
                            .variable("i", i)
                            .variable("pesos", format!("{} pesos", miles(pesos.len())))
                            .tabla(
                                "Las ventanas con más peso",
                                &["Ventana", "Posición", "Peso"],
                                filas,
                                &[(0, 2)],
                            )
                            .paso(
                                "pesos.push(probabilidad_kmer(ventana, &perfil));",
                                format!(
                                    "Cada una de las {} ventanas de la secuencia {} recibe un peso: su probabilidad según el perfil. La tabla muestra las que más pesan.",
                                    miles(pesos.len()),
                                    i + 1
                                ),
                            ),
                    );
                }
                let total: f64 = pesos.iter().sum();
                let probabilidad_de_salir = if total > 0.0 {
                    pesos[posiciones[i]] / total
                } else {
                    1.0 / pesos.len() as f64
                };
                narracion.agregar(
                    lista_de_kmers(
                        Pizarra::nueva()
                            .variable("k", k)
                            .variable("i", i)
                            .secuencia(
                                &textos[i],
                                &[
                                    (anterior, anterior + k, MarcaDeSecuencia::Comparado),
                                    (posiciones[i], posiciones[i] + k, MarcaDeSecuencia::Hallado),
                                ],
                            ),
                        "k-mers actuales",
                        textos,
                        &posiciones,
                        k,
                        &[i],
                    )
                    .paso(
                        "posiciones[i] = azar.elegir_con_pesos(&pesos);",
                        format!(
                            "Iteración {} (secuencia {}): sorteamos con esos pesos y salió {} (posición {}), que tenía {} de probabilidad de salir. Las ventanas con más peso salen más, pero cualquiera con peso mayor que 0 puede salir.",
                            iteracion + 1,
                            i + 1,
                            recortar(&textos[i][posiciones[i]..posiciones[i] + k], 40),
                            posiciones[i],
                            porcentaje(probabilidad_de_salir, 1)
                        ),
                    ),
                );
                if narracion.hay_lugar() {
                    narracion.agregar(
                        Pizarra::nueva()
                            .variable("k", k)
                            .variable("i", i)
                            .variable("puntos", puntos)
                            .variable("mejor_puntuacion", mejor_aqui)
                            .paso(
                                "let puntos = puntuacion(&motivos_en(adn, &posiciones, k));",
                                if puntos < mejor_aqui {
                                    format!(
                                        "Con el cambio la puntuación baja a {}: mejora a la mejor de este reinicio ({}) y la reemplaza.",
                                        miles(puntos),
                                        miles(mejor_aqui)
                                    )
                                } else {
                                    format!(
                                        "Con el cambio la puntuación es {}: no mejora a la mejor de este reinicio ({}). Gibbs acepta el cambio de todos modos y sigue.",
                                        miles(puntos),
                                        miles(mejor_aqui)
                                    )
                                },
                            ),
                    );
                }
            }
            if puntos < mejor_aqui {
                mejor_aqui = puntos;
                mejores_aqui.copy_from_slice(&posiciones);
                mejoras_aqui += 1;
                if detallado && iteracion >= 3 && mejoras_aqui <= 8 && narracion.hay_lugar() {
                    narracion.agregar(
                        lista_de_kmers(
                            Pizarra::nueva()
                                .variable("k", k)
                                .variable("i", i)
                                .variable("mejor_puntuacion", mejor_aqui)
                                .variable("puntos", puntos),
                            "mejores k-mers de este reinicio",
                            textos,
                            &mejores_aqui,
                            k,
                            &[],
                        )
                        .paso(
                            "mejores.copy_from_slice(&posiciones);",
                            format!(
                                "Iteración {}: la puntuación baja a {} y guardamos estas posiciones como las mejores de este reinicio.",
                                miles(iteracion + 1),
                                miles(puntos)
                            ),
                        ),
                    );
                }
            }
        }
        puntos_por_intento.push(mejor_aqui);
        let nuevo_record = mejor_aqui < mejor_puntuacion;
        if nuevo_record {
            mejor_puntuacion = mejor_aqui;
            mejores = mejores_aqui;
            records += 1;
        }
        if nuevo_record && reinicio > 0 && records <= 12 && narracion.hay_lugar() {
            narracion.agregar(
                lista_de_kmers(
                    Pizarra::nueva()
                        .variable("k", k)
                        .variable("mejor_puntuacion", mejor_puntuacion)
                        .variable("puntos", mejor_aqui),
                    "mejores k-mers hasta ahora",
                    textos,
                    &mejores,
                    k,
                    &[],
                )
                .paso(
                    "mejores = posiciones;",
                    format!(
                        "Reinicio {}: partió de {} y tras {} iteraciones llegó a {}. Es el mejor de todos los reinicios hasta ahora.",
                        miles(reinicio + 1),
                        miles(inicial),
                        miles(iteraciones),
                        miles(mejor_aqui)
                    ),
                ),
            );
        } else if (reinicio < 3 || nuevo_record) && narracion.hay_lugar() {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("k", k)
                    .variable("mejor_puntuacion", mejor_puntuacion)
                    .variable("puntos", mejor_aqui)
                    .paso(
                        "let (posiciones, puntos) = cadena_de_gibbs(adn, k, iteraciones, &mut azar);",
                        format!(
                            "Reinicio {}: partió de {} y tras {} iteraciones su mejor puntuación fue {}. {}",
                            miles(reinicio + 1),
                            miles(inicial),
                            miles(iteraciones),
                            miles(mejor_aqui),
                            if nuevo_record {
                                "Es la mejor hasta ahora.".to_string()
                            } else {
                                format!("La mejor hasta ahora es {}.", miles(mejor_puntuacion))
                            }
                        ),
                    ),
            );
        }
    }
    let veces = puntos_por_intento
        .iter()
        .filter(|&&puntos| puntos == mejor_puntuacion)
        .count();
    narracion.agregar(tabla_de_intentos(
        Pizarra::nueva()
            .variable("k", k)
            .variable("mejor_puntuacion", mejor_puntuacion),
        "Puntuación final de cada reinicio (las mejores)",
        "Reinicios",
        &puntos_por_intento,
    )
    .paso(
        "for _ in 0..reinicios {",
        format!(
            "Hicimos {} reinicios. La mejor puntuación, {}, salió en {} de ellos. Que varios reinicios lleguen al mismo sitio es buena señal; con otra semilla podrían salir otros.",
            miles(reinicios),
            miles(mejor_puntuacion),
            miles(veces)
        ),
    ));
    let consenso = motivos::consenso(&motivos::motivos_en(textos, &mejores, k));
    let ultimo = lista_de_kmers(
        Pizarra::nueva()
            .variable("k", k)
            .variable("mejor_puntuacion", mejor_puntuacion)
            .variable("consenso", recortar(&consenso, 40)),
        "Motivos encontrados",
        textos,
        &mejores,
        k,
        &[],
    )
    .paso(
        "Ok(MotivosEncontrados {",
        format!(
            "Resultado: los mejores motivos tienen puntuación {} y consenso {}. Mira «Qué significa» para saber si es una señal real.",
            miles(mejor_puntuacion),
            recortar(&consenso, 40)
        ),
    );
    (
        narracion.terminar(ultimo),
        MotivosEncontrados {
            posiciones: mejores,
            puntos: mejor_puntuacion,
            puntos_por_intento,
        },
    )
}

fn pasos_gibbs(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (conjunto, k, iteraciones, reinicios, semilla) = leer_gibbs(entrada)?;
    Ok(narrar_gibbs(&conjunto.textos, k, iteraciones, reinicios, semilla).0)
}

#[cfg(test)]
mod pruebas {
    use super::*;

    #[test]
    fn el_paso_a_paso_termina_igual_que_nucleo() {
        let entrada = Entrada::nueva(&[("adn", REGIONES_SOS), ("k", "16")]);
        let pasos = pasos_codiciosa_laplace(&entrada).unwrap();
        let conjunto = leer_conjunto(&entrada, "adn", 2).unwrap();
        let hallado = motivos::busqueda_codiciosa(&conjunto.textos, 16, 1.0).unwrap();
        let ultimo = &pasos.last().unwrap().estado;
        assert_eq!(ultimo["mejores"]["puntos"], json!(hallado.puntos));
        assert_eq!(ultimo["mejores"]["posiciones"], json!(hallado.posiciones));
        assert!(pasos.len() <= MAXIMO_PASOS);
        assert!(
            pasos
                .iter()
                .any(|paso| paso.estado["mejoro"] == json!(true))
        );
    }

    #[test]
    fn el_ejemplo_encuentra_las_cajas_sos() {
        let salida = crate::ejecutar(
            "motivos.codiciosa_laplace",
            json!({ "adn": REGIONES_SOS, "k": "16" }),
        )
        .unwrap();
        let texto_completo = salida.interpretacion.join(" ");
        assert!(texto_completo.contains("caja SOS"));
    }

    #[test]
    fn rechaza_k_grande_en_la_mediana() {
        let error = crate::ejecutar("motivos.mediana", json!({ "adn": REGIONES_SOS, "k": "11" }))
            .unwrap_err();
        assert!(error.contains("4^k"));
    }

    #[test]
    fn la_cola_de_poisson_distingue_lo_raro_de_lo_normal() {
        assert!(cola_de_poisson(0.05, 1) > AZAR_TOLERADO);
        assert!((cola_de_poisson(0.05, 1) - 0.0488).abs() < 0.001);
        assert!(cola_de_poisson(1.4e-7, 1) < AZAR_TOLERADO);
        assert!(cola_de_poisson(3.0, 10) < AZAR_TOLERADO);
        assert!(cola_de_poisson(3.0, 5) > AZAR_TOLERADO);
        assert_eq!(cola_de_poisson(2.0, 0), 1.0);
        assert_eq!(cola_de_poisson(1_000_000.0, 1_000_000), 1.0);
        assert!(cola_de_poisson(1_000_000.0, 1_002_000) < 0.05);
        assert_eq!(cola_de_poisson(0.0, 1), 0.0);
    }

    #[test]
    fn parecido_a_lexa_con_umbral_del_azar() {
        assert_eq!(parecido_a_lexa(b"CTGTATATAAATACAG"), Some(1));
        assert_eq!(parecido_a_lexa(b"ACTGTATA"), Some(0));
        assert_eq!(parecido_a_lexa(b"GGGCCCGGGCCCGGGC"), None);
        assert!(es_casi_palindromo(b"CTGTATATAAATACAG"));
        assert!(!es_casi_palindromo(b"AAAAAAAAAAAAAAAA"));
    }
    fn entrada_de(id: &str) -> Entrada {
        Entrada::desde_json(crate::ejemplo(id).unwrap()).unwrap()
    }

    fn revisar_tamano(pasos: &[Paso]) {
        assert!(!pasos.is_empty() && pasos.len() <= MAXIMO_PASOS);
        assert!(pasos.iter().all(|paso| !paso.nota.is_empty()));
    }

    #[test]
    fn la_enumeracion_narrada_da_lo_mismo_que_el_nucleo() {
        let entrada = entrada_de("motivos.enumeracion");
        let (conjunto, k, d) = leer_enumeracion(&entrada).unwrap();
        let (pasos, hallados) = narrar_enumeracion(&conjunto, k, d);
        revisar_tamano(&pasos);
        assert_eq!(
            hallados,
            motivos::enumerar_motivos(&conjunto.textos, k, d).unwrap()
        );
        assert!(!hallados.is_empty());
        let otra = Entrada::nueva(&[
            ("adn", "ACGTACGTACGGA\nTTACGTTCGACCA\nGGACGTACCTTAA"),
            ("k", "4"),
            ("d", "2"),
        ]);
        let (conjunto, k, d) = leer_enumeracion(&otra).unwrap();
        let (pasos, hallados) = narrar_enumeracion(&conjunto, k, d);
        revisar_tamano(&pasos);
        assert_eq!(
            hallados,
            motivos::enumerar_motivos(&conjunto.textos, k, d).unwrap()
        );
    }

    #[test]
    fn el_perfil_narrado_da_lo_mismo_que_el_nucleo() {
        let entrada = entrada_de("motivos.perfil");
        let conjunto = leer_matriz_de_motivos(&entrada).unwrap();
        let vistas: Vec<&[u8]> = conjunto.textos.iter().map(Vec::as_slice).collect();
        let (pasos, analisis) = narrar_perfil(&vistas);
        revisar_tamano(&pasos);
        assert_eq!(analisis, motivos::analizar_motivos(&vistas).unwrap());
    }

    #[test]
    fn la_distancia_narrada_da_lo_mismo_que_el_nucleo() {
        let entrada = entrada_de("motivos.distancia");
        let (patron, conjunto) = leer_distancia(&entrada).unwrap();
        let (pasos, total) = narrar_distancia(&patron, &conjunto);
        revisar_tamano(&pasos);
        assert_eq!(
            total,
            motivos::distancia_patron_textos(&patron, &conjunto.textos).unwrap()
        );
    }

    #[test]
    fn la_mediana_narrada_da_lo_mismo_que_el_nucleo() {
        let entrada = entrada_de("motivos.mediana");
        let (conjunto, k) = leer_mediana(&entrada).unwrap();
        let (pasos, mediana) = narrar_mediana(&conjunto, k).unwrap();
        revisar_tamano(&pasos);
        assert_eq!(
            mediana,
            motivos::cadena_mediana(&conjunto.textos, k).unwrap()
        );
    }

    #[test]
    fn el_mas_probable_narrado_da_lo_mismo_que_el_nucleo() {
        let entrada = entrada_de("motivos.mas_probable");
        let (adn, _, perfil) = leer_mas_probable(&entrada).unwrap();
        let (pasos, posicion) = narrar_mas_probable(&adn, &perfil);
        revisar_tamano(&pasos);
        assert_eq!(posicion, motivos::posicion_mas_probable(&adn, &perfil));
    }

    #[test]
    fn la_codiciosa_narrada_da_lo_mismo_que_el_nucleo() {
        let entrada = entrada_de("motivos.codiciosa");
        let conjunto = leer_conjunto(&entrada, "adn", 2).unwrap();
        let (pasos, hallado) = narrar_codiciosa(&conjunto.textos, 16, 0.0);
        revisar_tamano(&pasos);
        assert_eq!(
            hallado,
            motivos::busqueda_codiciosa(&conjunto.textos, 16, 0.0).unwrap()
        );
        let (pasos, hallado) = narrar_codiciosa(&conjunto.textos, 6, 0.0);
        revisar_tamano(&pasos);
        assert_eq!(
            hallado,
            motivos::busqueda_codiciosa(&conjunto.textos, 6, 0.0).unwrap()
        );
    }

    #[test]
    fn la_aleatoria_narrada_da_lo_mismo_que_el_nucleo() {
        let entrada = entrada_de("motivos.aleatoria");
        let (conjunto, k, repeticiones, semilla) = leer_aleatoria(&entrada).unwrap();
        let (pasos, hallado) = narrar_aleatoria(&conjunto.textos, k, repeticiones, semilla);
        revisar_tamano(&pasos);
        assert_eq!(
            hallado,
            motivos::busqueda_aleatoria(&conjunto.textos, k, repeticiones, semilla).unwrap()
        );
        let (pasos, hallado) = narrar_aleatoria(&conjunto.textos, 8, 1, 99);
        revisar_tamano(&pasos);
        assert_eq!(
            hallado,
            motivos::busqueda_aleatoria(&conjunto.textos, 8, 1, 99).unwrap()
        );
    }

    #[test]
    fn gibbs_narrado_da_lo_mismo_que_el_nucleo() {
        let entrada = entrada_de("motivos.gibbs");
        let (conjunto, k, iteraciones, reinicios, semilla) = leer_gibbs(&entrada).unwrap();
        let (pasos, hallado) = narrar_gibbs(&conjunto.textos, k, iteraciones, reinicios, semilla);
        revisar_tamano(&pasos);
        assert_eq!(
            hallado,
            motivos::muestreo_gibbs(&conjunto.textos, k, iteraciones, reinicios, semilla).unwrap()
        );
        let (pasos, hallado) = narrar_gibbs(&conjunto.textos, 8, 0, 2, 3);
        revisar_tamano(&pasos);
        assert_eq!(
            hallado,
            motivos::muestreo_gibbs(&conjunto.textos, 8, 0, 2, 3).unwrap()
        );
    }

    #[test]
    fn los_pasos_cuentan_con_pocas_secuencias_y_muchas_repeticiones() {
        let pasos = crate::pasos(
            "motivos.gibbs",
            json!({ "adn": REGIONES_SOS, "k": "8", "iteraciones": "2000", "reinicios": "300", "semilla": "5" }),
        )
        .unwrap();
        revisar_tamano(&pasos);
        let pasos = crate::pasos(
            "motivos.aleatoria",
            json!({ "adn": REGIONES_SOS, "k": "8", "repeticiones": "5000", "semilla": "5" }),
        )
        .unwrap();
        revisar_tamano(&pasos);
    }
}
