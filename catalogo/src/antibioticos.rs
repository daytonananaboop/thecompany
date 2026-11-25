use crate::formato::{cientifico, decimal, miles, porcentaje};
use crate::pizarra::{Marca as Resalte, Narracion, Pizarra};
use crate::{Algoritmo, Bloque, Campo, Entrada, Marca, Paso, Punto, Salida, Tipo};
use nucleo::antibioticos::{self, Lideres, MASAS_DISTINTAS, MASAS_ENTERAS};
use nucleo::ecoli;
use nucleo::secuencia::{aminoacido, complemento_inverso};
use serde_json::{Value, json};
use std::collections::{BTreeMap, HashSet};
use std::sync::OnceLock;

const ESPECTRO_PALMERA: &str = "0 71 71 97 113 129 131 156 168 168 184 227 239 244 260 281 285 315 324 352 356 373 395 412 416 444 453 483 487 508 524 529 541 584 600 600 612 637 639 655 671 697 697 768";
const ESPECTRO_PALMERA_CON_RUIDO: &str = "0 71 71 97 113 129 131 156 168 168 174 184 239 244 281 285 315 324 352 356 373 395 412 415 416 444 449 453 483 508 524 529 541 584 600 600 612 637 639 655 671 697 697 768";
const ESPECTRO_CON_RUIDO: &str = "0 71 71 71 113 114 115 131 185 186 202 202 227 227 269 273 298 300 300 316 316 340 342 359 371 383 387 387 429 431 454 454 455 456 488 498 500 501 502 502 525 527 569 573 585 616 640 640 656 658 683 684 687 696 727 729 729 754 754 770 771 771 800 825 841 842 842 843 855 855 885 885 885 956";
const ESPECTRO_NO_ESTANDAR: &str = "0 57 69 69 71 83 87 103 113 128 140 147 156 158 160 182 196 197 216 218 227 231 233 243 253 265 269 274 287 300 305 322 344 352 352 356 357 371 374 374 382 409 415 423 425 437 440 443 447 480 484 487 502 504 512 518 540 556 570 571 575 583 587 605 622 627 639 640 684 690 700 709 711 730 731 745 767 769 771 780 787 799 814 824 840 844 856 858 858 867 870 927";
const DISTANCIAS_DE_CORTES: &str =
    "1827 2501 4461 5122 6089 6949 6962 7916 9583 11410 12084 13038 13911 17499 20000";

const FORMAS: &[&str] = &["lineal", "cíclico"];
const MAXIMO_PEPTIDO: usize = 10_000;
const MAXIMO_PEPTIDO_ESPECTRO: usize = 1000;
const MAXIMO_SUBPEPTIDOS_EN_TABLA: usize = 30;
const MAXIMO_TEXTO_DE_PROTEINA: usize = 100_000;
const MAXIMO_FILAS: usize = 2000;
const MAXIMO_PASOS: usize = 400;
const MAXIMO_MOSTRADOS: usize = 30;
const MAXIMO_N: usize = 1000;
const MAXIMO_M: usize = 144;
const MAXIMO_MASAS_EXPERIMENTALES: usize = 100_000;
const MAXIMO_MASAS_LIDERES: usize = 5000;
const MASA_MADRE_MAXIMA_LIDERES: usize = 4000;
const MAXIMO_PEPTIDOS_A_RECORTAR: usize = 10_000;
const MAXIMO_TRABAJO_DE_RECORTE: usize = 20_000_000;
const MAXIMO_BARRAS: usize = 60;
const MAXIMO_DISTANCIAS: usize = 100_000;
const ABUNDANCIA_CARBONO_13: f64 = 0.0107;
const ABUNDANCIA_NITROGENO_15: f64 = 0.00364;
const PEPTIDOS_POR_SEGUNDO: f64 = 1e6;
const SEGUNDOS_POR_DIA: f64 = 86_400.0;
const DIAS_POR_VUELTA_AL_SOL: f64 = 365.25;
const MASA_DEL_AGUA: usize = 18;

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![
        Algoritmo {
            id: "antibioticos.traducir",
            capitulo: "antibioticos",
            nombre: "Traducir ARN a proteína",
            resumen: "Lee el ARN de tres en tres bases con el código genético y se detiene en el primer codón de parada.",
            campos: &[Campo {
                id: "arn",
                etiqueta: "ARN (o ADN de la hebra codificante)",
                tipo: Tipo::Secuencia,
                ejemplo: "@gen:thrL",
            }],
            funcion: "antibioticos::traducir_arn",
            ejecutar: traducir,
            pasos: Some(pasos_traducir),
        },
        Algoritmo {
            id: "antibioticos.codificar",
            capitulo: "antibioticos",
            nombre: "Qué ADN codifica un péptido",
            resumen: "Todos los tramos de un texto de ADN, en cualquiera de las dos hebras, que se traducen en el péptido.",
            campos: &[
                Campo {
                    id: "texto",
                    etiqueta: "Texto de ADN",
                    tipo: Tipo::Secuencia,
                    ejemplo: "@gen:hisL",
                },
                Campo {
                    id: "peptido",
                    etiqueta: "Péptido",
                    tipo: Tipo::Proteina,
                    ejemplo: "HHH",
                },
            ],
            funcion: "antibioticos::codificaciones_del_peptido",
            ejecutar: codificar,
            pasos: Some(pasos_codificar),
        },
        Algoritmo {
            id: "antibioticos.masa",
            capitulo: "antibioticos",
            nombre: "Masa de un péptido",
            resumen: "Fórmula química, masa entera y masa monoisotópica de un péptido lineal o cíclico.",
            campos: &[
                Campo {
                    id: "peptido",
                    etiqueta: "Péptido",
                    tipo: Tipo::Proteina,
                    ejemplo: "@proteina:trpL",
                },
                Campo {
                    id: "forma",
                    etiqueta: "Forma",
                    tipo: Tipo::Opcion(FORMAS),
                    ejemplo: "lineal",
                },
            ],
            funcion: "antibioticos::masa_monoisotopica",
            ejecutar: masa,
            pasos: Some(pasos_masa),
        },
        Algoritmo {
            id: "antibioticos.espectro_ciclico",
            capitulo: "antibioticos",
            nombre: "Espectro teórico de un péptido cíclico",
            resumen: "Las masas de todos los pedazos contiguos de un anillo de aminoácidos, como las vería un espectrómetro ideal.",
            campos: &[Campo {
                id: "peptido",
                etiqueta: "Péptido cíclico",
                tipo: Tipo::Proteina,
                ejemplo: "PALMERA",
            }],
            funcion: "antibioticos::espectro_ciclico",
            ejecutar: espectro_ciclico,
            pasos: Some(pasos_espectro_ciclico),
        },
        Algoritmo {
            id: "antibioticos.espectro_lineal",
            capitulo: "antibioticos",
            nombre: "Espectro de un péptido lineal",
            resumen: "Las masas de todos los pedazos contiguos de un péptido abierto, calculadas como restas de masas de prefijos.",
            campos: &[Campo {
                id: "peptido",
                etiqueta: "Péptido lineal",
                tipo: Tipo::Proteina,
                ejemplo: "@proteina:pheL",
            }],
            funcion: "antibioticos::espectro_lineal",
            ejecutar: espectro_lineal,
            pasos: Some(pasos_espectro_lineal),
        },
        Algoritmo {
            id: "antibioticos.contar",
            capitulo: "antibioticos",
            nombre: "Cuántos péptidos tienen una masa",
            resumen: "Cuenta con programación dinámica las secuencias de masas de aminoácido que suman exactamente una masa.",
            campos: &[Campo {
                id: "masa",
                etiqueta: "Masa (Da)",
                tipo: Tipo::Entero,
                ejemplo: "1000",
            }],
            funcion: "antibioticos::peptidos_por_masa",
            ejecutar: contar,
            pasos: Some(pasos_contar),
        },
        Algoritmo {
            id: "antibioticos.ciclopeptido",
            capitulo: "antibioticos",
            nombre: "Secuenciar un péptido cíclico (ramificación y poda)",
            resumen: "Todos los péptidos cuyo espectro cíclico es idéntico al dado: se alargan candidatos y se podan los inconsistentes.",
            campos: &[Campo {
                id: "espectro",
                etiqueta: "Espectro (masas enteras)",
                tipo: Tipo::Numeros,
                ejemplo: ESPECTRO_PALMERA,
            }],
            funcion: "antibioticos::secuenciar_ciclopeptido",
            ejecutar: ciclopeptido,
            pasos: Some(pasos_ciclopeptido),
        },
        Algoritmo {
            id: "antibioticos.puntaje_ciclico",
            capitulo: "antibioticos",
            nombre: "Puntaje de un péptido cíclico",
            resumen: "Cuántas masas comparten el espectro cíclico de un péptido y un espectro experimental, contando las repetidas.",
            campos: &[
                Campo {
                    id: "peptido",
                    etiqueta: "Péptido cíclico",
                    tipo: Tipo::Proteina,
                    ejemplo: "PALMERA",
                },
                Campo {
                    id: "espectro",
                    etiqueta: "Espectro experimental",
                    tipo: Tipo::Numeros,
                    ejemplo: ESPECTRO_PALMERA_CON_RUIDO,
                },
            ],
            funcion: "antibioticos::puntaje_ciclico",
            ejecutar: puntaje_ciclico,
            pasos: Some(pasos_puntaje_ciclico),
        },
        Algoritmo {
            id: "antibioticos.puntaje_lineal",
            capitulo: "antibioticos",
            nombre: "Puntaje de un péptido lineal",
            resumen: "Como el puntaje cíclico, pero con el espectro lineal: sin los pedazos que dan la vuelta al anillo.",
            campos: &[
                Campo {
                    id: "peptido",
                    etiqueta: "Péptido",
                    tipo: Tipo::Proteina,
                    ejemplo: "PALMERA",
                },
                Campo {
                    id: "espectro",
                    etiqueta: "Espectro experimental",
                    tipo: Tipo::Numeros,
                    ejemplo: ESPECTRO_PALMERA_CON_RUIDO,
                },
            ],
            funcion: "antibioticos::puntaje_lineal",
            ejecutar: puntaje_lineal,
            pasos: Some(pasos_puntaje_lineal),
        },
        Algoritmo {
            id: "antibioticos.recortar",
            capitulo: "antibioticos",
            nombre: "Recortar la tabla de líderes",
            resumen: "Se queda con los N péptidos de mejor puntaje lineal y con todos los que empatan con el último.",
            campos: &[
                Campo {
                    id: "peptidos",
                    etiqueta: "Péptidos (letras, o masas unidas por guiones)",
                    tipo: Tipo::Texto,
                    ejemplo: "PALM PAIM MERA AMLA PAMR WAG GGGG",
                },
                Campo {
                    id: "espectro",
                    etiqueta: "Espectro experimental",
                    tipo: Tipo::Numeros,
                    ejemplo: ESPECTRO_PALMERA_CON_RUIDO,
                },
                Campo {
                    id: "n",
                    etiqueta: "N (cuántos se quedan)",
                    tipo: Tipo::Entero,
                    ejemplo: "3",
                },
            ],
            funcion: "antibioticos::recortar",
            ejecutar: recortar,
            pasos: Some(pasos_recortar),
        },
        Algoritmo {
            id: "antibioticos.lideres",
            capitulo: "antibioticos",
            nombre: "Secuenciar con tabla de líderes",
            resumen: "Secuencia un péptido cíclico desde un espectro con ruido guardando en cada ronda solo los N mejores candidatos.",
            campos: &[
                Campo {
                    id: "espectro",
                    etiqueta: "Espectro experimental",
                    tipo: Tipo::Numeros,
                    ejemplo: ESPECTRO_CON_RUIDO,
                },
                Campo {
                    id: "n",
                    etiqueta: "N (tamaño de la tabla)",
                    tipo: Tipo::Entero,
                    ejemplo: "10",
                },
            ],
            funcion: "antibioticos::secuenciar_con_lideres",
            ejecutar: lideres,
            pasos: Some(pasos_lideres),
        },
        Algoritmo {
            id: "antibioticos.convolucion",
            capitulo: "antibioticos",
            nombre: "Convolución espectral",
            resumen: "Todas las diferencias positivas entre masas del espectro: las que más se repiten suelen ser masas de aminoácidos.",
            campos: &[Campo {
                id: "espectro",
                etiqueta: "Espectro experimental",
                tipo: Tipo::Numeros,
                ejemplo: ESPECTRO_NO_ESTANDAR,
            }],
            funcion: "antibioticos::convolucion_espectral",
            ejecutar: convolucion,
            pasos: Some(pasos_convolucion),
        },
        Algoritmo {
            id: "antibioticos.convolucion_lideres",
            capitulo: "antibioticos",
            nombre: "Secuenciar con convolución y tabla de líderes",
            resumen: "Toma de la convolución las M masas más frecuentes entre 57 y 200 Da y secuencia con ellas, aunque no sean aminoácidos estándar.",
            campos: &[
                Campo {
                    id: "espectro",
                    etiqueta: "Espectro experimental",
                    tipo: Tipo::Numeros,
                    ejemplo: ESPECTRO_NO_ESTANDAR,
                },
                Campo {
                    id: "m",
                    etiqueta: "M (masas que se toman de la convolución)",
                    tipo: Tipo::Entero,
                    ejemplo: "18",
                },
                Campo {
                    id: "n",
                    etiqueta: "N (tamaño de la tabla)",
                    tipo: Tipo::Entero,
                    ejemplo: "40",
                },
            ],
            funcion: "antibioticos::secuenciar_por_convolucion",
            ejecutar: convolucion_lideres,
            pasos: Some(pasos_convolucion_lideres),
        },
        Algoritmo {
            id: "antibioticos.autopista",
            capitulo: "antibioticos",
            nombre: "Puntos desde sus distancias (autopista de peaje)",
            resumen: "Reconstruye puntos sobre una recta a partir de todas sus distancias de a pares, con búsqueda con retroceso.",
            campos: &[Campo {
                id: "distancias",
                etiqueta: "Distancias entre todos los pares de puntos",
                tipo: Tipo::Numeros,
                ejemplo: DISTANCIAS_DE_CORTES,
            }],
            funcion: "antibioticos::puntos_desde_distancias",
            ejecutar: autopista,
            pasos: Some(pasos_autopista),
        },
    ]
}

fn nombre_aminoacido(letra: u8) -> &'static str {
    match letra {
        b'G' => "glicina",
        b'A' => "alanina",
        b'S' => "serina",
        b'P' => "prolina",
        b'V' => "valina",
        b'T' => "treonina",
        b'C' => "cisteína",
        b'I' => "isoleucina",
        b'L' => "leucina",
        b'N' => "asparagina",
        b'D' => "ácido aspártico",
        b'K' => "lisina",
        b'Q' => "glutamina",
        b'E' => "ácido glutámico",
        b'M' => "metionina",
        b'H' => "histidina",
        b'F' => "fenilalanina",
        b'R' => "arginina",
        b'Y' => "tirosina",
        b'W' => "triptófano",
        b'*' => "parada",
        _ => "desconocido",
    }
}

fn letras_de_masa(masa: usize) -> String {
    let letras: Vec<String> = MASAS_ENTERAS
        .iter()
        .filter(|(_, otra)| *otra == masa)
        .map(|(letra, _)| (*letra as char).to_string())
        .collect();
    match letras.len() {
        0 => format!("[{masa}]"),
        1 => letras[0].clone(),
        _ => format!("({})", letras.join("/")),
    }
}

fn letras_de_peptido(masas: &[usize]) -> String {
    masas.iter().map(|&masa| letras_de_masa(masa)).collect()
}

fn texto_de_masas(masas: &[usize]) -> String {
    masas
        .iter()
        .map(|masa| masa.to_string())
        .collect::<Vec<_>>()
        .join("-")
}

fn lista_de_masas(masas: &[usize]) -> String {
    masas
        .iter()
        .map(|masa| masa.to_string())
        .collect::<Vec<_>>()
        .join(", ")
}

fn como_texto(letras: &[u8]) -> String {
    String::from_utf8_lossy(letras).to_string()
}

fn como_arn(codon: &[u8]) -> String {
    codon
        .iter()
        .map(|&base| if base == b'T' { 'U' } else { base as char })
        .collect()
}

fn en_codones(adn: &[u8]) -> String {
    adn.chunks(3).map(como_texto).collect::<Vec<_>>().join(" ")
}

fn miles_grande(numero: u128) -> String {
    let digitos = numero.to_string();
    let mut resultado = String::new();
    for (posicion, digito) in digitos.chars().enumerate() {
        if posicion > 0 && (digitos.len() - posicion).is_multiple_of(3) {
            resultado.push('\u{a0}');
        }
        resultado.push(digito);
    }
    resultado
}

fn numero_grande(numero: f64) -> String {
    if !numero.is_finite() {
        return format!("más de {}", cientifico(f64::MAX));
    }
    if numero < 1e15 {
        miles(numero as usize)
    } else {
        cientifico(numero)
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
        filas: filas.into_iter().take(MAXIMO_FILAS).collect(),
        total,
    }
}

fn leer_peptido(entrada: &Entrada, campo: &str, maximo: usize) -> Result<Vec<u8>, String> {
    let peptido = entrada.proteina(campo)?;
    if peptido.len() > maximo {
        return Err(format!(
            "El péptido tiene {} aminoácidos y aquí el límite es {}.",
            miles(peptido.len()),
            miles(maximo)
        ));
    }
    antibioticos::masas_del_peptido(&peptido)?;
    Ok(peptido)
}

fn leer_espectro(entrada: &Entrada, campo: &str) -> Result<Vec<usize>, String> {
    let numeros = entrada.enteros(campo)?;
    if numeros.len() > MAXIMO_MASAS_EXPERIMENTALES {
        return Err(format!(
            "El campo «{campo}» tiene {} masas y el límite es {}.",
            miles(numeros.len()),
            miles(MAXIMO_MASAS_EXPERIMENTALES)
        ));
    }
    let mut masas = Vec::with_capacity(numeros.len());
    for numero in numeros {
        if numero < 0 {
            return Err(format!(
                "En el campo «{campo}», {numero} es negativo y una masa no puede serlo."
            ));
        }
        masas.push(numero as usize);
    }
    Ok(masas)
}

fn espectro_con_cero(mut espectro: Vec<usize>) -> (Vec<usize>, bool) {
    if espectro.contains(&0) {
        (espectro, false)
    } else {
        espectro.insert(0, 0);
        (espectro, true)
    }
}

fn leer_n(entrada: &Entrada, campo: &str, maximo: usize) -> Result<usize, String> {
    let n = entrada.natural(campo)?;
    if n == 0 || n > maximo {
        return Err(format!(
            "El campo «{campo}» debe estar entre 1 y {}.",
            miles(maximo)
        ));
    }
    Ok(n)
}

fn revisar_para_lideres(espectro: &[usize]) -> Result<(), String> {
    if espectro.len() > MAXIMO_MASAS_LIDERES {
        return Err(format!(
            "El espectro tiene {} masas; aquí el límite es {} para que la búsqueda termine en segundos.",
            miles(espectro.len()),
            miles(MAXIMO_MASAS_LIDERES)
        ));
    }
    let mayor = espectro.iter().copied().max().unwrap_or(0);
    if mayor > MASA_MADRE_MAXIMA_LIDERES {
        return Err(format!(
            "La masa mayor del espectro es {} Da; aquí el límite es {} Da (unos 35 aminoácidos, más que la nisina, un antibiótico peptídico de unos 3 350 Da) para que la búsqueda termine en segundos.",
            miles(mayor),
            miles(MASA_MADRE_MAXIMA_LIDERES)
        ));
    }
    Ok(())
}

fn forma_canonica(peptido: &[usize]) -> Vec<usize> {
    let mut mejor = peptido.to_vec();
    let mut invertido = peptido.to_vec();
    invertido.reverse();
    for mut lectura in [peptido.to_vec(), invertido] {
        for _ in 0..lectura.len() {
            lectura.rotate_left(1);
            if lectura < mejor {
                mejor = lectura.clone();
            }
        }
    }
    mejor
}

fn anillos_distintos(peptidos: &[Vec<usize>]) -> usize {
    peptidos
        .iter()
        .map(|peptido| forma_canonica(peptido))
        .collect::<HashSet<_>>()
        .len()
}

fn lecturas_del_anillo(peptido: &[usize]) -> usize {
    let mut lecturas = HashSet::new();
    let mut invertido = peptido.to_vec();
    invertido.reverse();
    for mut lectura in [peptido.to_vec(), invertido] {
        for _ in 0..lectura.len() {
            lectura.rotate_left(1);
            lecturas.insert(lectura.clone());
        }
    }
    lecturas.len()
}

fn explicar_ambiguedad(salida: &mut Salida, peptidos: &[&[usize]]) {
    let tiene_113 = peptidos.iter().any(|peptido| peptido.contains(&113));
    let tiene_128 = peptidos.iter().any(|peptido| peptido.contains(&128));
    if tiene_113 || tiene_128 {
        salida.explicar(
            "Con masas enteras, 113 puede ser I o L (son isómeros: la misma fórmula, C6H11NO) y 128 puede ser K o Q. El espectro no distingue esos pares, así que cada resultado representa varias secuencias de letras.",
        );
    }
}

fn comparar_espectros(
    teorico: &[usize],
    experimental: &[usize],
) -> BTreeMap<usize, (usize, usize)> {
    let mut tabla: BTreeMap<usize, (usize, usize)> = BTreeMap::new();
    for &masa in teorico {
        tabla.entry(masa).or_insert((0, 0)).0 += 1;
    }
    for &masa in experimental {
        tabla.entry(masa).or_insert((0, 0)).1 += 1;
    }
    tabla
}

fn faltantes_y_sobrantes(
    comparacion: &BTreeMap<usize, (usize, usize)>,
) -> (Vec<String>, Vec<String>) {
    let mut faltan = Vec::new();
    let mut sobran = Vec::new();
    for (&masa, &(teoricas, experimentales)) in comparacion {
        for _ in experimentales..teoricas {
            faltan.push(masa.to_string());
        }
        for _ in teoricas..experimentales {
            sobran.push(masa.to_string());
        }
    }
    (faltan, sobran)
}

fn grafica_de_espectro(etiqueta: &str, espectro: &[usize], marcas: Vec<Marca>) -> Bloque {
    let mut puntos = Vec::new();
    for grupo in espectro.chunk_by(|a, b| a == b) {
        let masa = grupo[0] as f64;
        puntos.push([masa - 0.3, 0.0]);
        puntos.push([masa, grupo.len() as f64]);
        puntos.push([masa + 0.3, 0.0]);
    }
    Bloque::linea(etiqueta, "Masa (Da)", "Veces que aparece", &puntos, marcas)
}

fn fraccion_en_e_coli(letra: u8) -> f64 {
    static FRECUENCIAS: OnceLock<[f64; 26]> = OnceLock::new();
    let frecuencias = FRECUENCIAS.get_or_init(|| {
        let mut cuentas = [0usize; 26];
        let mut total = 0usize;
        for elegido in ecoli::genes() {
            for &residuo in &elegido.proteina {
                if residuo.is_ascii_uppercase() {
                    cuentas[(residuo - b'A') as usize] += 1;
                    total += 1;
                }
            }
        }
        let mut frecuencias = [0.0; 26];
        for (frecuencia, cuenta) in frecuencias.iter_mut().zip(cuentas) {
            *frecuencia = cuenta as f64 / total.max(1) as f64;
        }
        frecuencias
    });
    if letra.is_ascii_uppercase() {
        frecuencias[(letra - b'A') as usize]
    } else {
        0.0
    }
}

fn composicion(peptido: &[u8]) -> Vec<(u8, usize)> {
    let mut cuentas: BTreeMap<u8, usize> = BTreeMap::new();
    for &letra in peptido {
        *cuentas.entry(letra).or_insert(0) += 1;
    }
    let mut lista: Vec<(u8, usize)> = cuentas.into_iter().collect();
    lista.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    lista
}

fn traducir(entrada: &Entrada) -> Result<Salida, String> {
    let leida = entrada.secuencia("arn")?;
    let traduccion = antibioticos::traducir_arn(&leida)?;
    let peptido = &traduccion.peptido;
    let codones_leidos = match traduccion.parada {
        Some(posicion) => posicion / 3 + 1,
        None => leida.len() / 3,
    };
    let mut salida = Salida::default();
    let texto = if peptido.is_empty() {
        "(ningún aminoácido)".to_string()
    } else if peptido.len() > MAXIMO_TEXTO_DE_PROTEINA {
        format!("{}…", como_texto(&peptido[..MAXIMO_TEXTO_DE_PROTEINA]))
    } else {
        como_texto(peptido)
    };
    salida.agregar(Bloque::texto("Proteína", texto));
    salida.agregar(Bloque::valor("Aminoácidos", miles(peptido.len())));
    let parada = match traduccion.parada {
        Some(posicion) => format!(
            "{} en la posición {} (contando desde 0)",
            como_arn(&leida[posicion..posicion + 3]),
            miles(posicion)
        ),
        None => "ninguno".to_string(),
    };
    salida.agregar(Bloque::valor("Codón de parada", parada));
    let mut filas = Vec::new();
    for numero in 0..codones_leidos.min(MAXIMO_FILAS) {
        let codon = &leida[numero * 3..numero * 3 + 3];
        let letra = aminoacido(codon);
        let simbolo = if letra == b'*' {
            "parada".to_string()
        } else {
            (letra as char).to_string()
        };
        filas.push(vec![
            miles(numero * 3),
            como_arn(codon),
            simbolo,
            nombre_aminoacido(letra).to_string(),
        ]);
    }
    salida.agregar(tabla_con_total(
        "Codones leídos",
        &[
            "Posición (contando desde 0)",
            "Codón",
            "Aminoácido",
            "Nombre",
        ],
        filas,
        codones_leidos,
    ));
    if !peptido.is_empty() {
        let lista = composicion(peptido);
        salida.agregar(Bloque::barras(
            "Composición de la proteína",
            "Aminoácido",
            "Veces",
            lista
                .iter()
                .map(|(letra, _)| (*letra as char).to_string())
                .collect(),
            lista.iter().map(|(_, veces)| *veces as f64).collect(),
        ));
    }
    explicar_traduccion(&mut salida, &leida, &traduccion);
    Ok(salida)
}

fn explicar_traduccion(salida: &mut Salida, leida: &[u8], traduccion: &antibioticos::Traduccion) {
    let peptido = &traduccion.peptido;
    if leida.contains(&b'T') {
        salida.explicar(
            "Tu secuencia trae T, así que es ADN: la leí como el ARN mensajero que se copia de ella, con U en lugar de T.",
        );
    }
    match traduccion.parada {
        Some(posicion) if posicion + 3 == leida.len() => salida.explicar(format!(
            "La última tripleta, {}, es un codón de parada: la traducción terminó justo al final, como en un gen completo que da {} aminoácidos.",
            como_arn(&leida[posicion..]),
            miles(peptido.len())
        )),
        Some(posicion) => salida.explicar(format!(
            "La traducción se detuvo en el codón de parada {} (posición {}, contando desde 0) y {} bases quedaron sin leer: el ribosoma se suelta en la primera parada, así que lo que sigue no forma parte de esta proteína.",
            como_arn(&leida[posicion..posicion + 3]),
            miles(posicion),
            miles(leida.len() - posicion - 3)
        )),
        None => {
            let mut texto = format!(
                "No apareció ningún codón de parada en {} codones: puede ser un trozo de un gen, o un marco de lectura distinto al del gen.",
                miles(leida.len() / 3)
            );
            let sobrantes = leida.len() % 3;
            if sobrantes > 0 {
                texto.push_str(&format!(
                    " Al final sobran {sobrantes} bases que no alcanzan para un codón."
                ));
            }
            salida.explicar(texto);
        }
    }
    if leida.len() >= 3 && traduccion.parada != Some(0) {
        let inicio = como_arn(&leida[..3]);
        match inicio.as_str() {
            "AUG" => {}
            "GUG" | "UUG" => salida.explicar(format!(
                "Empieza con {inicio} y no con AUG: en bacterias GUG y UUG también sirven de codón de inicio, y en ese caso el ribosoma pone metionina (M) aunque la tabla diga {}.",
                (aminoacido(&leida[..3]) as char)
            )),
            _ => salida.explicar(format!(
                "No empieza con AUG (empieza con {inicio}): si esperabas una proteína completa, quizá el marco de lectura está corrido; prueba empezando una o dos bases más adelante."
            )),
        }
    }
    if peptido.len() >= 8
        && let Some(anotado) = ecoli::genes().iter().find(|candidato| {
            candidato.proteina.len() == peptido.len() && candidato.proteina[1..] == peptido[1..]
        })
    {
        salida.explicar(format!(
            "Es idéntica a la proteína anotada del gen {} de E. coli ({}).",
            anotado.nombre, anotado.producto
        ));
    }
    if peptido.len() >= 10 {
        let lista = composicion(peptido);
        let (letra, veces) = lista[0];
        let fraccion = veces as f64 / peptido.len() as f64;
        let esperada = fraccion_en_e_coli(letra);
        if veces >= 3 && fraccion > 2.0 * esperada {
            salida.explicar(format!(
                "El {} de la proteína es {} ({}), contra {} en el conjunto de proteínas de E. coli: está muy enriquecida en ese aminoácido (más del doble de lo habitual).",
                porcentaje(fraccion, 1),
                nombre_aminoacido(letra),
                letra as char,
                porcentaje(esperada, 1)
            ));
        }
    }
}

fn inicio_en_e_coli(entrada: &Entrada, campo: &str) -> Option<usize> {
    let texto = entrada.texto(campo).ok()?;
    let limpio = texto.trim();
    if limpio == "@ecoli" {
        return Some(1);
    }
    let rango = limpio.strip_prefix("@ecoli:")?;
    let (inicio, _) = rango.split_once('-')?;
    inicio.trim().parse().ok()
}

fn gen_del_texto(entrada: &Entrada, campo: &str) -> Option<String> {
    let texto = entrada.texto(campo).ok()?;
    let nombre = texto.trim().strip_prefix("@gen:")?;
    Some(nombre.trim().to_string())
}

fn hay_solape(encontradas: &[antibioticos::Codificacion], largo: usize) -> bool {
    b"+-".iter().any(|&hebra| {
        let inicios: Vec<usize> = encontradas
            .iter()
            .filter(|codificacion| codificacion.hebra == hebra)
            .map(|codificacion| codificacion.inicio)
            .collect();
        inicios.windows(2).any(|par| par[1] - par[0] < largo)
    })
}

fn codificar(entrada: &Entrada) -> Result<Salida, String> {
    let texto = entrada.adn("texto")?;
    let peptido = leer_peptido(entrada, "peptido", MAXIMO_PEPTIDO)?;
    let encontradas = antibioticos::codificaciones_del_peptido(&texto, &peptido)?;
    let largo = 3 * peptido.len();
    let en_e_coli = inicio_en_e_coli(entrada, "texto");
    let mut columnas = vec![
        "Inicio (contando desde 0)",
        "Hebra",
        "ADN tal como está en el texto",
        "Codones que lee el ribosoma",
    ];
    if en_e_coli.is_some() {
        columnas.push("En E. coli (contando desde 1)");
    }
    let mut filas = Vec::new();
    let mut codificantes = HashSet::new();
    for codificacion in &encontradas {
        let codificante = if codificacion.hebra == b'+' {
            codificacion.fragmento.clone()
        } else {
            complemento_inverso(&codificacion.fragmento)
        };
        if filas.len() < MAXIMO_FILAS {
            let mut fila = vec![
                miles(codificacion.inicio),
                if codificacion.hebra == b'+' {
                    "+".to_string()
                } else {
                    "−".to_string()
                },
                como_texto(&codificacion.fragmento),
                en_codones(&codificante),
            ];
            if let Some(base) = en_e_coli {
                fila.push(format!(
                    "{}–{}",
                    miles(base + codificacion.inicio),
                    miles(base + codificacion.inicio + largo - 1)
                ));
            }
            filas.push(fila);
        }
        codificantes.insert(codificante);
    }
    let posibles = antibioticos::secuencias_que_codifican(&peptido);
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        "Tramos que codifican el péptido",
        miles(encontradas.len()),
    ));
    salida.agregar(tabla_con_total(
        "Codificaciones",
        &columnas,
        filas,
        encontradas.len(),
    ));
    salida.agregar(Bloque::valor(
        "Secuencias de ADN posibles para este péptido",
        numero_grande(posibles),
    ));
    let nombre = como_texto(&peptido);
    if encontradas.is_empty() {
        salida.explicar(format!(
            "Ningún tramo de {} bases, en ninguna de las dos hebras, se traduce en {nombre}. Si el péptido existe en la célula, puede que no lo fabrique el ribosoma: muchos antibióticos peptídicos, como los de este capítulo, los arman enzimas (síntesis no ribosómica) sin un gen que los deletree.",
            miles(largo)
        ));
    } else {
        let directas = encontradas
            .iter()
            .filter(|codificacion| codificacion.hebra == b'+')
            .count();
        let inversas = encontradas.len() - directas;
        salida.explicar(format!(
            "Hay {} tramos de {} bases que codifican {nombre}: {} en la hebra + (se leen tal cual en tu texto) y {} en la hebra − (se leen en el complemento inverso, de derecha a izquierda).",
            miles(encontradas.len()),
            miles(largo),
            miles(directas),
            miles(inversas)
        ));
        if codificantes.len() > 1 {
            salida.explicar(format!(
                "Entre ellos hay {} secuencias de ADN distintas para el mismo péptido: el código genético es redundante y varios codones sinónimos dan el mismo aminoácido.",
                miles(codificantes.len())
            ));
        }
        let se_solapan = hay_solape(&encontradas, largo);
        if se_solapan {
            salida.explicar(format!(
                "Algunos tramos se solapan (empiezan a menos de {} bases uno del otro): pasa cuando el péptido está hecho de repeticiones, como varias histidinas seguidas.",
                miles(largo)
            ));
        }
        if let Some(nombre_gen) = gen_del_texto(entrada, "texto") {
            let fuera = encontradas
                .iter()
                .filter(|codificacion| codificacion.hebra == b'-' || codificacion.inicio % 3 != 0)
                .count();
            if fuera == 0 {
                salida.explicar(format!(
                    "Todas están en el marco de lectura del gen {nombre_gen} (empiezan en posiciones múltiplos de 3, en la hebra +): son las que lee el ribosoma al fabricar su proteína."
                ));
            } else {
                salida.explicar(format!(
                    "{} están fuera del marco de lectura del gen {nombre_gen} (no empiezan en un múltiplo de 3 o están en la otra hebra): son coincidencias que el ribosoma no lee al fabricar esa proteína.",
                    miles(fuera)
                ));
            }
        }
        if let Some(base) = en_e_coli {
            salida.explicar(format!(
                "Las posiciones de la tabla cuentan desde 0 dentro de tu texto; en coordenadas de E. coli (que cuentan desde 1), la primera ocupa {}–{}.",
                miles(base + encontradas[0].inicio),
                miles(base + encontradas[0].inicio + largo - 1)
            ));
        }
    }
    salida.explicar(format!(
        "En teoría hay {} secuencias de ADN de {} bases que codifican este péptido (se multiplica cuántos codones tiene cada aminoácido): buscarlo en un genoma es buscar cualquiera de ellas.",
        numero_grande(posibles),
        miles(largo)
    ));
    Ok(salida)
}

fn elemento_en_formula(simbolo: &str, cuantos: u32) -> String {
    match cuantos {
        0 => String::new(),
        1 => simbolo.to_string(),
        _ => format!("{simbolo}{cuantos}"),
    }
}

fn texto_de_formula(formula: &antibioticos::Formula) -> String {
    [
        elemento_en_formula("C", formula.carbono),
        elemento_en_formula("H", formula.hidrogeno),
        elemento_en_formula("N", formula.nitrogeno),
        elemento_en_formula("O", formula.oxigeno),
        elemento_en_formula("S", formula.azufre),
    ]
    .concat()
}

fn masa(entrada: &Entrada) -> Result<Salida, String> {
    let peptido = leer_peptido(entrada, "peptido", MAXIMO_PEPTIDO)?;
    let forma = entrada.opcion("forma", FORMAS)?;
    let ciclico = forma == "cíclico";
    let formula = antibioticos::formula_del_peptido(&peptido, ciclico)?;
    let monoisotopica = antibioticos::masa_monoisotopica(&peptido, ciclico)?;
    let masas = antibioticos::masas_del_peptido(&peptido)?;
    let entera = antibioticos::masa_total(&masas) + if ciclico { 0 } else { MASA_DEL_AGUA };
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Fórmula", texto_de_formula(&formula)));
    salida.agregar(Bloque::valor(
        "Masa monoisotópica",
        format!("{} Da", decimal(monoisotopica, 4)),
    ));
    salida.agregar(Bloque::valor(
        "Masa entera",
        format!("{} Da", miles(entera)),
    ));
    salida.agregar(Bloque::valor("Aminoácidos", miles(peptido.len())));
    let mut filas = Vec::new();
    for (letra, veces) in composicion(&peptido) {
        let residuo = antibioticos::masa_monoisotopica(&[letra], true)?;
        filas.push(vec![
            (letra as char).to_string(),
            nombre_aminoacido(letra).to_string(),
            miles(veces),
            antibioticos::masa_entera(letra).unwrap_or(0).to_string(),
            decimal(residuo, 4),
        ]);
    }
    salida.agregar(Bloque::tabla(
        "Aminoácidos del péptido",
        &[
            "Letra",
            "Nombre",
            "Veces",
            "Masa entera del residuo",
            "Masa monoisotópica del residuo",
        ],
        filas,
    ));
    salida.explicar(format!(
        "La masa entera ({} Da) cuenta cada átomo con un número redondo (C 12, H 1, N 14, O 16, S 32). La monoisotópica ({} Da) usa la masa real del isótopo más común de cada elemento: la diferencia, {} Da, sale de que el hidrógeno pesa algo más de 1 (1.0078) y el oxígeno algo menos de 16 (15.9949).",
        miles(entera),
        decimal(monoisotopica, 4),
        decimal(monoisotopica - entera as f64, 4)
    ));
    if ciclico {
        salida.explicar(
            "Como es cíclico no se suma agua: cada enlace peptídico libera una molécula de agua y al cerrar el anillo se pierde una más que en el péptido lineal.",
        );
    } else {
        salida.explicar(
            "Como es lineal, a la suma de los residuos se le agrega una molécula de agua (18.0106 Da): un H en el extremo amino y un OH en el extremo carboxilo.",
        );
    }
    if peptido.contains(&b'K') || peptido.contains(&b'Q') {
        salida.explicar(
            "Tiene K o Q: con masa entera los dos pesan 128, pero su masa exacta difiere en 0.036 Da (lisina 128.095, glutamina 128.059). Un espectrómetro de alta resolución sí los separa.",
        );
    }
    if peptido.contains(&b'I') || peptido.contains(&b'L') {
        salida.explicar(
            "Tiene I o L: son isómeros (los dos C6H11NO), así que pesan exactamente lo mismo y ninguna medida de masa los distingue.",
        );
    }
    let pesados = ABUNDANCIA_CARBONO_13 * formula.carbono as f64
        + ABUNDANCIA_NITROGENO_15 * formula.nitrogeno as f64;
    let sin_pesados = (-pesados).exp();
    if pesados >= 1.0 {
        salida.explicar(format!(
            "Con {} carbonos y {} nitrógenos, cada molécula trae en promedio {} átomos de isótopos pesados (sobre todo carbono 13) y solo el {} no trae ninguno: en el espectro, el pico de un dalton más será al menos tan alto como el monoisotópico, que deja de ser el más visible.",
            miles(formula.carbono as usize),
            miles(formula.nitrogeno as usize),
            decimal(pesados, 2),
            porcentaje(sin_pesados, 0)
        ));
    } else {
        salida.explicar(format!(
            "Con {} carbonos y {} nitrógenos, cada molécula trae en promedio {} átomos de isótopos pesados y el {} no trae ninguno: el pico monoisotópico es el más alto de su grupo.",
            miles(formula.carbono as usize),
            miles(formula.nitrogeno as usize),
            decimal(pesados, 2),
            porcentaje(sin_pesados, 0)
        ));
    }
    if peptido.len() >= 20
        && peptido
            .iter()
            .all(|letra| matches!(letra, b'A' | b'C' | b'G' | b'T'))
    {
        salida.explicar(
            "Tu péptido solo tiene las letras A, C, G y T: ¿no será ADN? Como proteína, esas letras son alanina, cisteína, glicina y treonina.",
        );
    }
    Ok(salida)
}

fn subpeptidos_con_letras(peptido: &[u8], ciclico: bool) -> Vec<(String, usize)> {
    let n = peptido.len();
    let mut lista = vec![("(vacío)".to_string(), 0)];
    let largo_maximo = if ciclico { n.saturating_sub(1) } else { n };
    for inicio in 0..n {
        for largo in 1..=largo_maximo {
            if !ciclico && inicio + largo > n {
                break;
            }
            let mut letras = Vec::with_capacity(largo);
            let mut suma = 0;
            for desplazamiento in 0..largo {
                let letra = peptido[(inicio + desplazamiento) % n];
                letras.push(letra);
                suma += antibioticos::masa_entera(letra).unwrap_or(0);
            }
            lista.push((como_texto(&letras), suma));
        }
    }
    if ciclico && n > 0 {
        let total = peptido
            .iter()
            .map(|&letra| antibioticos::masa_entera(letra).unwrap_or(0))
            .sum();
        lista.push((como_texto(peptido), total));
    }
    lista.sort_by_key(|(_, masa)| *masa);
    lista
}

fn espectro_ciclico(entrada: &Entrada) -> Result<Salida, String> {
    espectro_teorico(entrada, true)
}

fn espectro_lineal(entrada: &Entrada) -> Result<Salida, String> {
    espectro_teorico(entrada, false)
}

fn espectro_teorico(entrada: &Entrada, ciclico: bool) -> Result<Salida, String> {
    let peptido = leer_peptido(entrada, "peptido", MAXIMO_PEPTIDO_ESPECTRO)?;
    let masas = antibioticos::masas_del_peptido(&peptido)?;
    let espectro = if ciclico {
        antibioticos::espectro_ciclico(&masas)
    } else {
        antibioticos::espectro_lineal(&masas)
    };
    let total = antibioticos::masa_total(&masas);
    let n = peptido.len();
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Masas en el espectro", miles(espectro.len())));
    salida.agregar(Bloque::valor(
        "Masa del péptido entero",
        format!("{} Da", miles(total)),
    ));
    salida.agregar(grafica_de_espectro(
        "Espectro",
        &espectro,
        vec![Marca {
            x: total as f64,
            texto: "péptido entero".to_string(),
        }],
    ));
    salida.agregar(Bloque::lista(
        "Espectro (masas de menor a mayor)",
        espectro.iter().map(|masa| masa.to_string()),
    ));
    if n <= MAXIMO_SUBPEPTIDOS_EN_TABLA {
        let filas = subpeptidos_con_letras(&peptido, ciclico)
            .into_iter()
            .map(|(letras, masa)| vec![letras, masa.to_string()])
            .collect();
        salida.agregar(Bloque::tabla(
            "Subpéptidos y sus masas",
            &["Subpéptido", "Masa (Da)"],
            filas,
        ));
    }
    if n == 0 {
        salida.explicar("El péptido no tiene aminoácidos: su espectro es solo la masa 0.");
        return Ok(salida);
    }
    if ciclico {
        if n == 1 {
            salida.explicar(
                "Con un solo aminoácido el espectro tiene dos masas: 0 (el pedazo vacío) y la del aminoácido.",
            );
        } else {
            salida.explicar(format!(
                "Un anillo de {n} aminoácidos tiene {n}·{} = {} pedazos que no son ni vacíos ni el anillo entero (cada uno de los {n} puntos de partida, con largos de 1 a {}). Sumando el vacío (masa 0) y el anillo completo, el espectro tiene {} masas.",
                n - 1,
                miles(n * (n - 1)),
                n - 1,
                miles(espectro.len())
            ));
        }
    } else {
        let prefijos: Vec<String> = (0..=n.min(5))
            .map(|largo| antibioticos::masa_total(&masas[..largo]).to_string())
            .collect();
        salida.explicar(format!(
            "Un péptido lineal de {n} aminoácidos tiene {n}·{}/2 = {} pedazos contiguos; con el vacío son {} masas. Si fuera cíclico serían {}, porque también contarían los pedazos que dan la vuelta.",
            n + 1,
            miles(n * (n + 1) / 2),
            miles(espectro.len()),
            miles(n * (n - 1) + 2)
        ));
        salida.explicar(format!(
            "Basta con las masas de los prefijos ({}{}) para obtener todo el espectro: la masa de cada pedazo es la resta de dos prefijos.",
            prefijos.join(", "),
            if n > 5 { ", …" } else { "" }
        ));
    }
    let distintas = espectro.chunk_by(|a, b| a == b).count();
    if distintas < espectro.len() {
        let mas_repetida = espectro
            .chunk_by(|a, b| a == b)
            .max_by_key(|grupo| grupo.len())
            .unwrap_or(&[0]);
        salida.explicar(format!(
            "Solo {} masas son distintas: subpéptidos diferentes pueden pesar lo mismo, por eso el espectro es una lista con repeticiones y no un conjunto. La más repetida, {}, aparece {} veces.",
            miles(distintas),
            mas_repetida[0],
            miles(mas_repetida.len())
        ));
    }
    explicar_ambiguedad(&mut salida, &[&masas]);
    if ciclico && n >= 2 {
        let lecturas = lecturas_del_anillo(&masas);
        let mut texto = format!(
            "Hay {} maneras de leer este anillo (empezando en cada aminoácido y en los dos sentidos) y todas dan exactamente este espectro: el espectro no dice dónde empieza el anillo ni hacia qué lado se lee.",
            miles(lecturas)
        );
        if lecturas < 2 * n {
            texto.push_str(&format!(
                " Son menos de 2·{n} = {} porque el anillo tiene repeticiones o simetría.",
                2 * n
            ));
        }
        salida.explicar(texto);
        if espectro.len() > antibioticos::MAXIMO_MASAS_EXACTAS {
            salida.explicar(format!(
                "Este espectro tiene más de {} masas: la ramificación y poda no lo acepta; para reconstruirlo usa la tabla de líderes.",
                antibioticos::MAXIMO_MASAS_EXACTAS
            ));
        }
    }
    Ok(salida)
}

fn contar(entrada: &Entrada) -> Result<Salida, String> {
    let masa = entrada.natural("masa")?;
    let cuentas = antibioticos::peptidos_por_masa(masa)?;
    let total = cuentas[masa];
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        &format!("Péptidos de {} Da", miles(masa)),
        miles_grande(total),
    ));
    let puntos: Vec<[f64; 2]> = cuentas
        .iter()
        .enumerate()
        .filter(|(_, cuenta)| **cuenta > 0)
        .map(|(otra, cuenta)| [otra as f64, (*cuenta as f64).log10()])
        .collect();
    salida.agregar(Bloque::linea(
        "Cuántos péptidos hay de cada masa",
        "Masa (Da)",
        "log₁₀ del número de péptidos",
        &puntos,
        vec![Marca {
            x: masa as f64,
            texto: format!("{masa} Da"),
        }],
    ));
    let mut filas = Vec::new();
    for otra in (0..=masa).step_by(100) {
        filas.push(vec![miles(otra), miles_grande(cuentas[otra])]);
    }
    if !masa.is_multiple_of(100) {
        filas.push(vec![miles(masa), miles_grande(total)]);
    }
    salida.agregar(Bloque::tabla(
        "Cada 100 Da",
        &["Masa (Da)", "Péptidos"],
        filas,
    ));
    if masa == 0 {
        salida.explicar(
            "Con masa 0 solo existe el péptido vacío: por eso la cuenta de cada masa empieza en 1.",
        );
    } else if total == 0 {
        if masa < antibioticos::MASA_MINIMA_DE_AMINOACIDO {
            salida.explicar(format!(
                "Ningún péptido pesa {masa} Da: el aminoácido más liviano, la glicina, ya pesa 57."
            ));
        } else {
            salida.explicar(format!(
                "Ninguna secuencia de las 18 masas suma exactamente {} Da: no toda masa entera se puede armar con aminoácidos.",
                miles(masa)
            ));
        }
    }
    if masa >= 400 && total > 0 && cuentas[masa - 100] > 0 {
        let razon = (total as f64 / cuentas[masa - 100] as f64).powf(0.01);
        salida.explicar(format!(
            "Cada dalton de más multiplica la cantidad de péptidos por unas {} veces; cada 100 Da, por unas {}. Es un crecimiento exponencial: en escala logarítmica la curva es casi una recta.",
            decimal(razon, 3),
            decimal(razon.powi(100), 0)
        ));
    }
    let dias = total as f64 / PEPTIDOS_POR_SEGUNDO / SEGUNDOS_POR_DIA;
    if dias >= DIAS_POR_VUELTA_AL_SOL {
        salida.explicar(format!(
            "Aunque revisaras un millón de péptidos por segundo, tardarías unos {} años en verlos todos: por eso el capítulo pasa a la ramificación y poda, que descarta ramas enteras sin mirarlas.",
            numero_grande(dias / DIAS_POR_VUELTA_AL_SOL)
        ));
    }
    salida.explicar(
        "Se cuentan 18 masas y no 20 letras porque I y L pesan 113 y K y Q pesan 128: cada péptido contado es una secuencia de masas. El orden importa: GA y AG cuentan como dos.",
    );
    Ok(salida)
}

struct Ronda {
    numero: usize,
    expandidos: usize,
    inconsistentes: usize,
    llegaron: usize,
    nuevos: usize,
    conservados: usize,
    encontrados: usize,
    muestra_expandidos: Vec<Vec<usize>>,
    muestra_conservados: Vec<Vec<usize>>,
}

struct Busqueda {
    masas: Vec<usize>,
    masa_madre: usize,
    bloques: Vec<usize>,
    rondas: Vec<Ronda>,
    encontrados: Vec<Vec<usize>>,
}

fn buscar_por_ramificacion(espectro: &[usize]) -> Result<Busqueda, String> {
    let mut experimental = antibioticos::espectro_para_ramificar(espectro)?;
    let bloques = antibioticos::aminoacidos_presentes(&experimental);
    let masa_madre = experimental.masa_madre;
    let mut candidatos: Vec<Vec<usize>> = vec![Vec::new()];
    let mut encontrados: Vec<Vec<usize>> = Vec::new();
    let mut rondas = Vec::new();
    while !candidatos.is_empty() {
        let expandidos = antibioticos::expandir(&candidatos, &bloques)?;
        let llegaron = expandidos
            .iter()
            .filter(|peptido| antibioticos::masa_total(peptido) == masa_madre)
            .count();
        let cantidad = expandidos.len();
        let antes = encontrados.len();
        let muestra_expandidos = expandidos.iter().take(MAXIMO_MOSTRADOS).cloned().collect();
        candidatos =
            antibioticos::conservar_consistentes(expandidos, &mut experimental, &mut encontrados);
        rondas.push(Ronda {
            numero: rondas.len() + 1,
            expandidos: cantidad,
            inconsistentes: cantidad - llegaron - candidatos.len(),
            llegaron,
            nuevos: encontrados.len() - antes,
            conservados: candidatos.len(),
            encontrados: encontrados.len(),
            muestra_expandidos,
            muestra_conservados: candidatos.iter().take(MAXIMO_MOSTRADOS).cloned().collect(),
        });
    }
    Ok(Busqueda {
        masas: experimental.masas.clone(),
        masa_madre,
        bloques,
        rondas,
        encontrados,
    })
}

fn ciclopeptido(entrada: &Entrada) -> Result<Salida, String> {
    let (espectro, agregado) = espectro_con_cero(leer_espectro(entrada, "espectro")?);
    let busqueda = buscar_por_ramificacion(&espectro)?;
    let encontrados = &busqueda.encontrados;
    let anillos = anillos_distintos(encontrados);
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        "Péptidos encontrados",
        miles(encontrados.len()),
    ));
    salida.agregar(Bloque::valor("Anillos distintos", miles(anillos)));
    let filas = encontrados
        .iter()
        .map(|peptido| {
            vec![
                texto_de_masas(peptido),
                letras_de_peptido(peptido),
                miles(antibioticos::masa_total(peptido)),
            ]
        })
        .collect();
    salida.agregar(Bloque::tabla(
        "Péptidos que explican el espectro",
        &["Masas", "Letras posibles", "Masa total (Da)"],
        filas,
    ));
    let puntos: Vec<[f64; 2]> = busqueda
        .rondas
        .iter()
        .map(|ronda| [ronda.numero as f64, ronda.expandidos as f64])
        .collect();
    let marcas = busqueda
        .rondas
        .iter()
        .filter(|ronda| ronda.llegaron > 0)
        .map(|ronda| Marca {
            x: ronda.numero as f64,
            texto: "llegan a la masa madre".to_string(),
        })
        .collect();
    salida.agregar(Bloque::linea(
        "Candidatos revisados en cada ronda",
        "Ronda (largo del candidato)",
        "Candidatos",
        &puntos,
        marcas,
    ));
    let filas = busqueda
        .rondas
        .iter()
        .map(|ronda| {
            vec![
                miles(ronda.numero),
                miles(ronda.expandidos),
                miles(ronda.inconsistentes),
                miles(ronda.llegaron),
                miles(ronda.nuevos),
                miles(ronda.conservados),
            ]
        })
        .collect();
    salida.agregar(Bloque::tabla(
        "Rondas",
        &[
            "Ronda",
            "Candidatos",
            "Podados por inconsistentes",
            "Llegan a la masa madre",
            "Encontrados",
            "Siguen",
        ],
        filas,
    ));
    if agregado {
        salida.explicar(
            "Agregué la masa 0 (el subpéptido vacío): todo espectro teórico la incluye y sin ella ningún candidato sería consistente.",
        );
    }
    if encontrados.is_empty() {
        let mut texto = "Ningún péptido tiene un espectro cíclico idéntico a este. Un espectro real casi siempre trae masas de más (ruido) o de menos, y basta una para que la búsqueda exacta falle: prueba la secuenciación con tabla de líderes, que tolera errores.".to_string();
        if busqueda.bloques.is_empty() {
            texto.push_str(
                " Además, ninguna masa del espectro es la de un aminoácido: no hay con qué empezar.",
            );
        }
        salida.explicar(texto);
    } else {
        let largo = encontrados[0].len();
        if anillos == 1 {
            salida.explicar(format!(
                "Los {} péptidos son el mismo anillo leído desde distintos puntos de partida y en los dos sentidos (con {largo} aminoácidos hay como mucho {} lecturas). El espectro determina el péptido cíclico, salvo por dónde empieza.",
                miles(encontrados.len()),
                2 * largo
            ));
        } else {
            salida.explicar(format!(
                "Hay {} anillos distintos con exactamente este espectro: el espectro solo no alcanza para decidir cuál es el verdadero.",
                miles(anillos)
            ));
        }
        let vistas: Vec<&[usize]> = encontrados
            .iter()
            .map(|peptido| peptido.as_slice())
            .collect();
        explicar_ambiguedad(&mut salida, &vistas);
    }
    let revisados: usize = busqueda.rondas.iter().map(|ronda| ronda.expandidos).sum();
    if busqueda.masa_madre <= antibioticos::MASA_MAXIMA_PARA_CONTAR {
        let fuerza_bruta =
            antibioticos::peptidos_por_masa(busqueda.masa_madre)?[busqueda.masa_madre];
        if fuerza_bruta > revisados as u128 {
            salida.explicar(format!(
                "La búsqueda revisó {} candidatos en {} rondas. Por fuerza bruta habría que revisar los {} péptidos que pesan {} Da: la poda descarta casi todos sin mirarlos.",
                miles(revisados),
                miles(busqueda.rondas.len()),
                miles_grande(fuerza_bruta),
                miles(busqueda.masa_madre)
            ));
        }
    }
    Ok(salida)
}

struct Momento<'a> {
    ronda: usize,
    fase: &'a str,
    candidatos: &'a [Vec<usize>],
    total_candidatos: usize,
    expandidos: usize,
    inconsistentes: usize,
    llegaron: usize,
    nuevos: usize,
    encontrados: usize,
}

fn momento_de_poda<'a>(ronda: &'a Ronda, fase: &'a str) -> Momento<'a> {
    Momento {
        ronda: ronda.numero,
        fase,
        candidatos: &ronda.muestra_conservados,
        total_candidatos: ronda.conservados,
        expandidos: ronda.expandidos,
        inconsistentes: ronda.inconsistentes,
        llegaron: ronda.llegaron,
        nuevos: ronda.nuevos,
        encontrados: ronda.encontrados,
    }
}

fn peptido_en_json(peptido: &[usize]) -> Value {
    json!({
        "masas": peptido,
        "masa": antibioticos::masa_total(peptido),
        "letras": letras_de_peptido(peptido),
    })
}

fn estado_de_busqueda(busqueda: &Busqueda, momento: &Momento) -> Value {
    let encontrados: Vec<Value> = busqueda.encontrados[..momento.encontrados]
        .iter()
        .take(MAXIMO_MOSTRADOS)
        .map(|peptido| peptido_en_json(peptido))
        .collect();
    json!({
        "ronda": momento.ronda,
        "fase": momento.fase,
        "masa_madre": busqueda.masa_madre,
        "espectro": busqueda.masas,
        "bloques": busqueda
            .bloques
            .iter()
            .map(|&masa| json!({"masa": masa, "letras": letras_de_masa(masa)}))
            .collect::<Vec<_>>(),
        "candidatos": momento
            .candidatos
            .iter()
            .take(MAXIMO_MOSTRADOS)
            .map(|peptido| peptido_en_json(peptido))
            .collect::<Vec<_>>(),
        "total_candidatos": momento.total_candidatos,
        "expandidos": momento.expandidos,
        "inconsistentes": momento.inconsistentes,
        "llegaron_a_masa_madre": momento.llegaron,
        "descartados": momento.inconsistentes + momento.llegaron - momento.nuevos,
        "nuevos_encontrados": momento.nuevos,
        "encontrados": encontrados,
        "total_encontrados": momento.encontrados,
    })
}

fn contados(numero: usize, singular: &str, plural: &str) -> String {
    if numero == 1 {
        format!("1 {singular}")
    } else {
        format!("{} {plural}", miles(numero))
    }
}

fn pasos_ciclopeptido(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (espectro, agregado) = espectro_con_cero(leer_espectro(entrada, "espectro")?);
    let busqueda = buscar_por_ramificacion(&espectro)?;
    let mut pasos = Vec::new();
    let vacio = vec![Vec::new()];
    let lista_de_bloques = if busqueda.bloques.is_empty() {
        "ninguna".to_string()
    } else {
        busqueda
            .bloques
            .iter()
            .map(|&masa| masa_con_letras(masa))
            .collect::<Vec<_>>()
            .join(", ")
    };
    let mut nota = format!(
        "El espectro tiene {} masas; la mayor, {} Da, es la del péptido entero (la masa madre). De las 18 masas de aminoácido aparecen {}: {lista_de_bloques}. Solo esas pueden formar el péptido, así que son los bloques para alargar a los candidatos.",
        miles(busqueda.masas.len()),
        miles(busqueda.masa_madre),
        busqueda.bloques.len()
    );
    if agregado {
        nota.push_str(" Agregué la masa 0, que todo espectro teórico incluye.");
    }
    pasos.push(Paso::nuevo(
        "let bloques = aminoacidos_presentes(&experimental);",
        nota,
        estado_de_busqueda(
            &busqueda,
            &Momento {
                ronda: 0,
                fase: "inicio",
                candidatos: &vacio,
                total_candidatos: 1,
                expandidos: 0,
                inconsistentes: 0,
                llegaron: 0,
                nuevos: 0,
                encontrados: 0,
            },
        ),
    ));
    let mut anteriores_candidatos = 1;
    let mut anteriores_encontrados = 0;
    for ronda in &busqueda.rondas {
        let origen = if ronda.numero == 1 {
            "el péptido vacío se alarga".to_string()
        } else if anteriores_candidatos == 1 {
            "el único candidato que quedó se alarga".to_string()
        } else {
            format!(
                "cada uno de los {} candidatos se alarga",
                miles(anteriores_candidatos)
            )
        };
        pasos.push(Paso::nuevo(
            "candidatos = expandir(&candidatos, &bloques)?;",
            format!(
                "Ronda {}: {origen} con cada uno de los {} bloques y salen {} de {}.",
                ronda.numero,
                busqueda.bloques.len(),
                contados(ronda.expandidos, "péptido", "péptidos"),
                contados(ronda.numero, "aminoácido", "aminoácidos")
            ),
            estado_de_busqueda(
                &busqueda,
                &Momento {
                    ronda: ronda.numero,
                    fase: "expandir",
                    candidatos: &ronda.muestra_expandidos,
                    total_candidatos: ronda.expandidos,
                    expandidos: ronda.expandidos,
                    inconsistentes: 0,
                    llegaron: 0,
                    nuevos: 0,
                    encontrados: anteriores_encontrados,
                },
            ),
        ));
        let mut nota = if ronda.inconsistentes == 0 {
            "No se poda ningún candidato: el espectro lineal de cada uno cabe en el experimental."
                .to_string()
        } else {
            format!(
                "Se podan {}: su espectro lineal pide masas que el espectro no tiene (o más copias de las que tiene), y alargarlos no lo arreglaría.",
                contados(ronda.inconsistentes, "candidato", "candidatos")
            )
        };
        if ronda.llegaron > 0 {
            nota.push_str(&format!(
                " {} a la masa madre y salen de la lista: {} el espectro exactamente.",
                contados(ronda.llegaron, "llegó", "llegaron"),
                if ronda.nuevos == 1 {
                    "1 explica".to_string()
                } else {
                    format!("{} explican", miles(ronda.nuevos))
                }
            ));
        }
        nota.push_str(&format!(
            " {} para la ronda siguiente.",
            match ronda.conservados {
                0 => "No queda ninguno".to_string(),
                1 => "Queda 1".to_string(),
                otros => format!("Quedan {}", miles(otros)),
            }
        ));
        pasos.push(Paso::nuevo(
            "} else if experimental.es_consistente(&peptido) {",
            nota,
            estado_de_busqueda(&busqueda, &momento_de_poda(ronda, "podar")),
        ));
        if ronda.nuevos > 0 {
            pasos.push(Paso::nuevo(
                "encontrados.push(peptido);",
                format!(
                    "{} de {} Da con un espectro cíclico idéntico al experimental: se guardan como solución.",
                    contados(ronda.nuevos, "péptido", "péptidos"),
                    miles(busqueda.masa_madre)
                ),
                estado_de_busqueda(&busqueda, &momento_de_poda(ronda, "encontrados")),
            ));
        }
        anteriores_candidatos = ronda.conservados;
        anteriores_encontrados = ronda.encontrados;
    }
    let rondas = busqueda.rondas.len();
    pasos.push(Paso::nuevo(
        "Ok(encontrados)",
        format!(
            "No quedan candidatos: la búsqueda termina tras {} con {}.",
            contados(rondas, "ronda", "rondas"),
            contados(
                busqueda.encontrados.len(),
                "péptido encontrado",
                "péptidos encontrados"
            )
        ),
        estado_de_busqueda(
            &busqueda,
            &Momento {
                ronda: rondas,
                fase: "fin",
                candidatos: &[],
                total_candidatos: 0,
                expandidos: 0,
                inconsistentes: 0,
                llegaron: 0,
                nuevos: 0,
                encontrados: busqueda.encontrados.len(),
            },
        ),
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

fn puntaje_ciclico(entrada: &Entrada) -> Result<Salida, String> {
    puntaje(entrada, true)
}

fn puntaje_lineal(entrada: &Entrada) -> Result<Salida, String> {
    puntaje(entrada, false)
}

fn puntaje(entrada: &Entrada, ciclico: bool) -> Result<Salida, String> {
    let peptido = leer_peptido(entrada, "peptido", MAXIMO_PEPTIDO_ESPECTRO)?;
    let masas = antibioticos::masas_del_peptido(&peptido)?;
    let espectro = leer_espectro(entrada, "espectro")?;
    let calcular = |candidato: &[usize]| {
        if ciclico {
            antibioticos::puntaje_ciclico(candidato, &espectro)
        } else {
            antibioticos::puntaje_lineal(candidato, &espectro)
        }
    };
    let obtenido = calcular(&masas)?;
    let teorico = if ciclico {
        antibioticos::espectro_ciclico(&masas)
    } else {
        antibioticos::espectro_lineal(&masas)
    };
    let comparacion = comparar_espectros(&teorico, &espectro);
    let (faltan, sobran) = faltantes_y_sobrantes(&comparacion);
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        "Puntaje",
        format!(
            "{} de {} masas del espectro",
            miles(obtenido),
            miles(espectro.len())
        ),
    ));
    salida.agregar(Bloque::valor(
        if ciclico {
            "Masas del espectro cíclico teórico"
        } else {
            "Masas del espectro lineal teórico"
        },
        miles(teorico.len()),
    ));
    let filas = comparacion
        .iter()
        .map(|(&masa, &(teoricas, experimentales))| {
            vec![
                masa.to_string(),
                miles(teoricas),
                miles(experimentales),
                miles(teoricas.min(experimentales)),
            ]
        })
        .collect();
    salida.agregar(Bloque::tabla(
        "Comparación masa por masa",
        &[
            "Masa (Da)",
            "En el teórico",
            "En el experimental",
            "Compartidas",
        ],
        filas,
    ));
    salida.agregar(Bloque::lista(
        "Masas teóricas que faltan en el experimento",
        faltan.iter(),
    ));
    salida.agregar(Bloque::lista(
        "Masas del experimento que el péptido no explica",
        sobran.iter(),
    ));
    if faltan.is_empty() && sobran.is_empty() {
        salida.explicar(
            "Puntaje perfecto: el espectro teórico y el experimental son idénticos, así que este péptido explica todas las masas sin que sobre ni falte ninguna.",
        );
    } else {
        salida.explicar(format!(
            "El péptido explica {} de las {} masas del espectro experimental ({}). Sobran {} masas sin explicación (ruido, o señal de que el péptido no es el correcto) y faltan {} masas teóricas que el espectrómetro no registró.",
            miles(obtenido),
            miles(espectro.len()),
            porcentaje(obtenido as f64 / espectro.len() as f64, 1),
            miles(sobran.len()),
            miles(faltan.len())
        ));
    }
    if masas.len() >= 3 {
        let mut ordenado = masas.clone();
        ordenado.sort_unstable();
        if forma_canonica(&ordenado) != forma_canonica(&masas) {
            let otro = calcular(&ordenado)?;
            if otro < obtenido {
                salida.explicar(format!(
                    "Los mismos aminoácidos ordenados de menor a mayor masa ({}) obtienen {}: el orden de tu péptido explica {} masas más, así que el espectro sí apoya ese orden.",
                    letras_de_peptido(&ordenado),
                    miles(otro),
                    miles(obtenido - otro)
                ));
            } else {
                salida.explicar(format!(
                    "Los mismos aminoácidos ordenados de menor a mayor masa ({}) obtienen {}, igual o más: el espectro no apoya en particular el orden de tu péptido.",
                    letras_de_peptido(&ordenado),
                    miles(otro)
                ));
            }
        }
    }
    if !espectro.contains(&0) {
        salida.explicar(
            "Tu espectro no trae la masa 0. El espectro teórico siempre la incluye (el subpéptido vacío), así que esa masa no puede coincidir.",
        );
    }
    if !ciclico {
        let en_anillo = antibioticos::puntaje_ciclico(&masas, &espectro)?;
        salida.explicar(format!(
            "El puntaje lineal no cuenta los pedazos que dan la vuelta al anillo: sirve para comparar péptidos todavía abiertos, como los candidatos de la tabla de líderes. Si cerraras este péptido en un anillo, su puntaje cíclico sería {}.",
            miles(en_anillo)
        ));
    }
    Ok(salida)
}

fn leer_lista_de_peptidos(texto: &str) -> Result<Vec<(String, Vec<usize>)>, String> {
    let mut peptidos = Vec::new();
    let mut trabajo = 0usize;
    for pieza in texto.split(|c: char| c.is_whitespace() || c == ',' || c == ';') {
        if pieza.is_empty() {
            continue;
        }
        if peptidos.len() == MAXIMO_PEPTIDOS_A_RECORTAR {
            return Err(format!(
                "Hay más de {} péptidos; ese es el límite.",
                miles(MAXIMO_PEPTIDOS_A_RECORTAR)
            ));
        }
        let masas = if pieza.chars().all(|c| c.is_ascii_digit() || c == '-') {
            let mut masas = Vec::new();
            for numero in pieza.split('-') {
                let masa: usize = numero
                    .parse()
                    .map_err(|_| format!("«{pieza}» no es una lista de masas como 113-128-186."))?;
                if masa == 0 || masa > antibioticos::MASA_MAXIMA {
                    return Err(format!(
                        "En «{pieza}», cada masa debe estar entre 1 y {} Da.",
                        miles(antibioticos::MASA_MAXIMA)
                    ));
                }
                masas.push(masa);
            }
            masas
        } else {
            let letras: Vec<u8> = pieza
                .bytes()
                .map(|letra| letra.to_ascii_uppercase())
                .collect();
            antibioticos::masas_del_peptido(&letras)
                .map_err(|error| format!("En «{pieza}»: {error}"))?
        };
        if masas.len() > MAXIMO_PEPTIDO_ESPECTRO {
            return Err(format!(
                "«{}…» tiene más de {} aminoácidos, el límite.",
                pieza.chars().take(12).collect::<String>(),
                miles(MAXIMO_PEPTIDO_ESPECTRO)
            ));
        }
        trabajo += masas.len() * (masas.len() + 1) / 2;
        if trabajo > MAXIMO_TRABAJO_DE_RECORTE {
            return Err(format!(
                "Los péptidos de la lista son demasiado largos en conjunto: sumando los pedazos de cada uno hay más de {} y ese es el límite. Usa menos péptidos o más cortos.",
                miles(MAXIMO_TRABAJO_DE_RECORTE)
            ));
        }
        peptidos.push((pieza.to_string(), masas));
    }
    if peptidos.is_empty() {
        return Err("No hay ningún péptido en la lista.".to_string());
    }
    Ok(peptidos)
}

fn recortar(entrada: &Entrada) -> Result<Salida, String> {
    let texto = entrada.texto("peptidos")?;
    let peptidos = leer_lista_de_peptidos(&texto)?;
    let espectro = leer_espectro(entrada, "espectro")?;
    let n = entrada.natural("n")?;
    if n == 0 {
        return Err("N debe ser al menos 1.".to_string());
    }
    let tabla: Vec<Vec<usize>> = peptidos.iter().map(|(_, masas)| masas.clone()).collect();
    let conservados = antibioticos::recortar(tabla.clone(), &espectro, n)?;
    let mut puntajes = Vec::with_capacity(tabla.len());
    for masas in &tabla {
        puntajes.push(antibioticos::puntaje_lineal(masas, &espectro)?);
    }
    let mut se_queda = vec![false; tabla.len()];
    let mut siguiente = 0;
    for (indice, masas) in tabla.iter().enumerate() {
        if siguiente < conservados.len() && conservados[siguiente] == *masas {
            se_queda[indice] = true;
            siguiente += 1;
        }
    }
    let mut orden: Vec<usize> = (0..tabla.len()).collect();
    orden.sort_by(|a, b| puntajes[*b].cmp(&puntajes[*a]));
    let corte = orden
        .iter()
        .filter(|&&indice| se_queda[indice])
        .map(|&indice| puntajes[indice])
        .min();
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor(
        "Se quedan",
        format!("{} de {}", miles(conservados.len()), miles(tabla.len())),
    ));
    if let Some(corte) = corte {
        salida.agregar(Bloque::valor("Puntaje de corte", miles(corte)));
    }
    let filas = orden
        .iter()
        .map(|&indice| {
            vec![
                peptidos[indice].0.clone(),
                texto_de_masas(&tabla[indice]),
                miles(puntajes[indice]),
                if se_queda[indice] { "sí" } else { "no" }.to_string(),
            ]
        })
        .collect();
    salida.agregar(Bloque::tabla(
        "Péptidos de mayor a menor puntaje lineal",
        &["Péptido", "Masas", "Puntaje lineal", "¿Se queda?"],
        filas,
    ));
    if tabla.len() <= MAXIMO_BARRAS {
        salida.agregar(Bloque::barras(
            "Puntaje lineal de cada péptido",
            "Péptido",
            "Puntaje",
            orden
                .iter()
                .map(|&indice| peptidos[indice].0.clone())
                .collect(),
            orden
                .iter()
                .map(|&indice| puntajes[indice] as f64)
                .collect(),
        ));
    }
    if n >= tabla.len() {
        salida.explicar(format!(
            "N = {} es mayor o igual que la cantidad de péptidos ({}): no se descarta ninguno.",
            miles(n),
            miles(tabla.len())
        ));
    } else {
        salida.explicar(format!(
            "Se ordenan por puntaje lineal y el corte queda en el puntaje del {}.º, {}: se quedan todos los que tienen al menos eso.",
            miles(n),
            miles(corte.unwrap_or(0))
        ));
        if conservados.len() > n {
            salida.explicar(format!(
                "Pediste {}, pero se quedan {}: los que empatan con el último también entran, para no descartar al azar candidatos igual de buenos.",
                miles(n),
                miles(conservados.len())
            ));
        }
    }
    let mut gemelos = None;
    for i in 0..peptidos.len() {
        for j in i + 1..peptidos.len() {
            if peptidos[i].1 == peptidos[j].1 && !peptidos[i].0.eq_ignore_ascii_case(&peptidos[j].0)
            {
                gemelos = Some((i, j));
                break;
            }
        }
        if gemelos.is_some() {
            break;
        }
    }
    if let Some((i, j)) = gemelos {
        salida.explicar(format!(
            "«{}» y «{}» tienen exactamente las mismas masas (I y L pesan 113; K y Q, 128): para el espectro son el mismo péptido y siempre empatan.",
            peptidos[i].0, peptidos[j].0
        ));
    }
    Ok(salida)
}

fn mostrar_lideres(salida: &mut Salida, resultado: &Lideres, espectro: &[usize]) {
    let Some(lider) = resultado.peptidos.first() else {
        salida.agregar(Bloque::valor("Péptido líder", "ninguno"));
        return;
    };
    salida.agregar(Bloque::valor("Péptido líder", texto_de_masas(lider)));
    salida.agregar(Bloque::valor("Letras posibles", letras_de_peptido(lider)));
    salida.agregar(Bloque::valor(
        "Puntaje",
        format!(
            "{} de {} masas",
            miles(resultado.puntaje),
            miles(espectro.len())
        ),
    ));
    let filas = resultado
        .peptidos
        .iter()
        .take(MAXIMO_FILAS)
        .map(|peptido| vec![texto_de_masas(peptido), letras_de_peptido(peptido)])
        .collect();
    salida.agregar(tabla_con_total(
        "Empatados con el mejor puntaje",
        &["Masas", "Letras posibles"],
        filas,
        resultado.peptidos.len(),
    ));
    let teorico = antibioticos::espectro_ciclico(lider);
    let mut ordenado = espectro.to_vec();
    ordenado.sort_unstable();
    salida.agregar(grafica_de_espectro(
        "Espectro experimental",
        &ordenado,
        Vec::new(),
    ));
    let (faltan, sobran) = faltantes_y_sobrantes(&comparar_espectros(&teorico, espectro));
    salida.agregar(Bloque::lista(
        "Masas del espectro que el líder no explica",
        sobran.iter(),
    ));
    salida.agregar(Bloque::lista(
        "Masas del líder que faltan en el espectro",
        faltan.iter(),
    ));
}

fn explicar_lideres(salida: &mut Salida, resultado: &Lideres, espectro: &[usize]) {
    let Some(lider) = resultado.peptidos.first() else {
        salida.explicar(format!(
            "Ningún candidato llegó a pesar exactamente {} Da, la masa mayor del espectro. Revisa que esa masa sea la del péptido completo: la tabla de líderes la usa como meta.",
            miles(espectro.iter().copied().max().unwrap_or(0))
        ));
        return;
    };
    let teorico = antibioticos::espectro_ciclico(lider);
    let faltan = teorico.len() - resultado.puntaje;
    let sobran = espectro.len() - resultado.puntaje;
    if faltan == 0 && sobran == 0 {
        salida.explicar(
            "Puntaje perfecto: el espectro cíclico del líder es idéntico al experimental, así que lo explica completo.",
        );
    } else {
        salida.explicar(format!(
            "El líder explica {} de las {} masas ({}). {} masas del espectro quedan sin explicar (ruido) y faltan {} que el líder predice pero el espectrómetro no registró: con un espectro real eso es lo esperable, y por eso aquí no se exige un puntaje perfecto.",
            miles(resultado.puntaje),
            miles(espectro.len()),
            porcentaje(resultado.puntaje as f64 / espectro.len() as f64, 1),
            miles(sobran),
            miles(faltan)
        ));
    }
    if resultado.peptidos.len() > 1 {
        let anillos = anillos_distintos(&resultado.peptidos);
        if anillos == 1 {
            salida.explicar(format!(
                "Los {} empatados son el mismo anillo leído desde otro punto de partida o al revés: el resultado no es ambiguo.",
                miles(resultado.peptidos.len())
            ));
        } else {
            explicar_empate(salida, &resultado.peptidos, anillos);
        }
    }
    explicar_ambiguedad(salida, &[lider]);
}

fn masas_que_sobran(uno: &[usize], otro: &[usize]) -> Vec<usize> {
    let mut restantes = otro.to_vec();
    let mut sobran = Vec::new();
    for &masa in uno {
        match restantes.iter().position(|&otra| otra == masa) {
            Some(posicion) => {
                restantes.swap_remove(posicion);
            }
            None => sobran.push(masa),
        }
    }
    sobran.sort_unstable();
    sobran
}

fn como_suma(masas: &[usize]) -> String {
    masas
        .iter()
        .map(|masa| masa.to_string())
        .collect::<Vec<_>>()
        .join(" + ")
}

fn explicar_empate(salida: &mut Salida, empatados: &[Vec<usize>], anillos: usize) {
    let corto = empatados.iter().min_by_key(|peptido| peptido.len());
    let largo = empatados.iter().max_by_key(|peptido| peptido.len());
    let (Some(corto), Some(largo)) = (corto, largo) else {
        return;
    };
    if corto.len() == largo.len() {
        salida.explicar(format!(
            "Empatan {} anillos distintos del mismo largo: con este ruido el espectro no alcanza para elegir uno. Un N más grande deja más candidatos en carrera en cada ronda y a veces encuentra uno que explica más masas.",
            miles(anillos)
        ));
        return;
    }
    salida.explicar(format!(
        "Empatan {} anillos distintos, y no todos tienen el mismo largo ({} y {} aminoácidos): {} en uno pesa lo mismo que {} en otro. Un N más grande no cambia el puntaje de esos anillos, así que seguirán empatados salvo que aparezca uno mejor; hace falta un espectro con menos ruido para decidir.",
        miles(anillos),
        corto.len(),
        largo.len(),
        como_suma(&masas_que_sobran(corto, largo)),
        como_suma(&masas_que_sobran(largo, corto))
    ));
}

fn lideres(entrada: &Entrada) -> Result<Salida, String> {
    let espectro = leer_espectro(entrada, "espectro")?;
    revisar_para_lideres(&espectro)?;
    let n = leer_n(entrada, "n", MAXIMO_N)?;
    let resultado = antibioticos::secuenciar_con_lideres(&espectro, n).map_err(|error| {
        format!("{error} Prueba con un N más chico: los empates con el último de la tabla también se quedan y la hacen crecer.")
    })?;
    let mut salida = Salida::default();
    mostrar_lideres(&mut salida, &resultado, &espectro);
    explicar_lideres(&mut salida, &resultado, &espectro);
    Ok(salida)
}

fn nombre_de_diferencia(masa: usize) -> String {
    let letras: Vec<String> = MASAS_ENTERAS
        .iter()
        .filter(|(_, otra)| *otra == masa)
        .map(|(letra, _)| (*letra as char).to_string())
        .collect();
    if letras.is_empty() {
        "no estándar".to_string()
    } else {
        letras.join(" o ")
    }
}

fn masa_con_letras(masa: usize) -> String {
    format!("{masa} ({})", nombre_de_diferencia(masa))
}

fn convolucion(entrada: &Entrada) -> Result<Salida, String> {
    let espectro = leer_espectro(entrada, "espectro")?;
    if espectro.len() > MAXIMO_MASAS_LIDERES {
        return Err(format!(
            "El espectro tiene {} masas; aquí el límite es {}.",
            miles(espectro.len()),
            miles(MAXIMO_MASAS_LIDERES)
        ));
    }
    let resultado = antibioticos::convolucion_espectral(&espectro)?;
    let total: usize = resultado.iter().map(|&(_, veces)| veces).sum();
    let mut en_rango: Vec<(usize, usize)> = resultado
        .iter()
        .copied()
        .filter(|(diferencia, _)| {
            (antibioticos::MASA_MINIMA_DE_AMINOACIDO..=antibioticos::MASA_MAXIMA_DE_AMINOACIDO)
                .contains(diferencia)
        })
        .collect();
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Diferencias positivas", miles(total)));
    salida.agregar(Bloque::valor(
        "Diferencias distintas entre 57 y 200 Da",
        miles(en_rango.len()),
    ));
    let mut por_masa = en_rango.clone();
    por_masa.sort_unstable();
    salida.agregar(Bloque::barras(
        "Diferencias entre 57 y 200 Da",
        "Diferencia (Da)",
        "Veces",
        por_masa
            .iter()
            .map(|(diferencia, _)| {
                if MASAS_DISTINTAS.contains(diferencia) {
                    format!("{diferencia} {}", letras_de_masa(*diferencia))
                } else {
                    diferencia.to_string()
                }
            })
            .collect(),
        por_masa.iter().map(|(_, veces)| *veces as f64).collect(),
    ));
    let filas = resultado
        .iter()
        .take(MAXIMO_FILAS)
        .map(|&(diferencia, veces)| {
            vec![
                diferencia.to_string(),
                miles(veces),
                if (antibioticos::MASA_MINIMA_DE_AMINOACIDO
                    ..=antibioticos::MASA_MAXIMA_DE_AMINOACIDO)
                    .contains(&diferencia)
                {
                    nombre_de_diferencia(diferencia)
                } else {
                    "fuera de rango".to_string()
                },
            ]
        })
        .collect();
    salida.agregar(tabla_con_total(
        "Convolución, de la más a la menos frecuente",
        &["Diferencia (Da)", "Veces", "Aminoácido"],
        filas,
        resultado.len(),
    ));
    if resultado.is_empty() {
        salida.explicar(
            "Todas las masas son iguales (o hay una sola): no hay ninguna diferencia positiva.",
        );
        return Ok(salida);
    }
    if en_rango.is_empty() {
        salida.explicar(
            "Ninguna diferencia cae entre 57 y 200 Da, el rango de las masas de aminoácidos: este espectro no sugiere ningún aminoácido.",
        );
        return Ok(salida);
    }
    en_rango.truncate(10);
    let (primera, veces_primera) = en_rango[0];
    salida.explicar(format!(
        "Si dos masas del espectro corresponden a un pedazo del péptido y al mismo pedazo con un aminoácido más, su diferencia es la masa de ese aminoácido. Por eso las masas verdaderas se repiten mucho: la diferencia más frecuente entre 57 y 200 Da es {primera} ({}), que aparece {} veces.",
        nombre_de_diferencia(primera),
        miles(veces_primera)
    ));
    let mejores: Vec<String> = en_rango
        .iter()
        .take(5)
        .map(|&(diferencia, veces)| {
            format!(
                "{diferencia} ×{} ({})",
                miles(veces),
                nombre_de_diferencia(diferencia)
            )
        })
        .collect();
    salida.explicar(format!(
        "Las más frecuentes en ese rango: {}.",
        mejores.join(", ")
    ));
    let raras: Vec<String> = en_rango
        .iter()
        .filter(|(diferencia, _)| !MASAS_DISTINTAS.contains(diferencia))
        .map(|(diferencia, _)| diferencia.to_string())
        .collect();
    if !raras.is_empty() {
        salida.explicar(format!(
            "{} de las {} más frecuentes ({}) no son masas de los 20 aminoácidos estándar: pueden ser aminoácidos modificados, comunes en péptidos no ribosómicos, o coincidencias del ruido. La secuenciación con convolución las prueba igual.",
            raras.len(),
            en_rango.len(),
            raras.join(", ")
        ));
    }
    Ok(salida)
}

fn convolucion_lideres(entrada: &Entrada) -> Result<Salida, String> {
    let espectro = leer_espectro(entrada, "espectro")?;
    revisar_para_lideres(&espectro)?;
    let m = leer_n(entrada, "m", MAXIMO_M)?;
    let n = leer_n(entrada, "n", MAXIMO_N)?;
    let resultado = antibioticos::secuenciar_por_convolucion(&espectro, m, n).map_err(|error| {
        format!("{error} Prueba con un M o un N más chico: cada masa del alfabeto multiplica los candidatos de cada ronda.")
    })?;
    let estandar = antibioticos::secuenciar_con_lideres(&espectro, n)?;
    let mut salida = Salida::default();
    mostrar_lideres(&mut salida, &resultado, &espectro);
    salida.agregar(Bloque::lista(
        "Alfabeto elegido por la convolución",
        resultado.alfabeto.iter().map(|&masa| masa_con_letras(masa)),
    ));
    salida.agregar(Bloque::valor(
        "Puntaje con las 18 masas estándar",
        if estandar.peptidos.is_empty() {
            "sin péptido".to_string()
        } else {
            miles(estandar.puntaje)
        },
    ));
    explicar_lideres(&mut salida, &resultado, &espectro);
    if resultado.alfabeto.len() > m {
        salida.explicar(format!(
            "La convolución aportó {} masas aunque pediste M = {m}: las que empatan con la {m}.ª en frecuencia también entran.",
            resultado.alfabeto.len()
        ));
    }
    if let Some(lider) = resultado.peptidos.first() {
        let mut raras: Vec<usize> = lider
            .iter()
            .copied()
            .filter(|masa| !MASAS_DISTINTAS.contains(masa))
            .collect();
        raras.sort_unstable();
        raras.dedup();
        if !raras.is_empty() {
            salida.explicar(format!(
                "El líder usa masas que no son de ninguno de los 20 aminoácidos estándar ({}): la tabla de líderes con las 18 masas de siempre nunca podría encontrarlo; la convolución las descubrió en los datos.",
                lista_de_masas(&raras)
            ));
            let mut conocidas = Vec::new();
            if raras.contains(&69) {
                conocidas.push("69 Da es la deshidroalanina");
            }
            if raras.contains(&83) {
                conocidas.push("83 Da es la deshidrobutirina");
            }
            if !conocidas.is_empty() {
                salida.explicar(format!(
                    "{}: aminoácidos modificados que aparecen en antibióticos peptídicos como la nisina.",
                    conocidas.join(" y ")
                ));
            }
        }
        if !estandar.peptidos.is_empty() {
            if resultado.puntaje > estandar.puntaje {
                salida.explicar(format!(
                    "Con las 18 masas estándar el mejor puntaje es {}; con el alfabeto de la convolución sube a {}.",
                    miles(estandar.puntaje),
                    miles(resultado.puntaje)
                ));
            } else {
                salida.explicar(format!(
                    "Con las 18 masas estándar el mejor puntaje es {}, no menor que el de la convolución ({}): aquí el alfabeto estándar alcanza.",
                    miles(estandar.puntaje),
                    miles(resultado.puntaje)
                ));
            }
        }
    }
    Ok(salida)
}

fn autopista(entrada: &Entrada) -> Result<Salida, String> {
    let distancias = entrada.enteros("distancias")?;
    if distancias.len() > MAXIMO_DISTANCIAS {
        return Err(format!(
            "Hay {} distancias y el límite es {}.",
            miles(distancias.len()),
            miles(MAXIMO_DISTANCIAS)
        ));
    }
    let ignoradas = distancias
        .iter()
        .filter(|&&distancia| distancia <= 0)
        .count();
    let resultado = antibioticos::puntos_desde_distancias(&distancias)?;
    let mut salida = Salida::default();
    if ignoradas > 0 {
        salida.explicar(format!(
            "Ignoré {} valores menores o iguales a cero: si das las diferencias de todos con todos (incluidos ceros y negativos), alcanza con las positivas.",
            miles(ignoradas)
        ));
    }
    let Some(puntos) = resultado else {
        salida.agregar(Bloque::valor("Puntos", "no hay solución"));
        salida.explicar(
            "No existe ningún conjunto de puntos con exactamente estas distancias. En un experimento real, como una digestión parcial con enzimas de restricción, basta una distancia perdida o de más para que no haya solución exacta.",
        );
        return Ok(salida);
    };
    let ancho = puntos[puntos.len() - 1];
    salida.agregar(Bloque::valor("Puntos", miles(puntos.len())));
    salida.agregar(Bloque::valor(
        "Distancia entre el primero y el último",
        ancho,
    ));
    salida.agregar(Bloque::lista(
        "Posiciones",
        puntos.iter().map(|punto| punto.to_string()),
    ));
    salida.agregar(Bloque::dispersion(
        "Puntos sobre la recta",
        "Posición",
        "",
        puntos
            .iter()
            .map(|&punto| Punto {
                x: punto as f64,
                y: 0.0,
                grupo: 0,
                nombre: punto.to_string(),
            })
            .collect(),
    ));
    let filas = puntos
        .windows(2)
        .map(|par| {
            vec![
                par[0].to_string(),
                par[1].to_string(),
                (par[1] - par[0]).to_string(),
            ]
        })
        .collect();
    salida.agregar(Bloque::tabla(
        "Tramos entre puntos vecinos",
        &["Desde", "Hasta", "Largo"],
        filas,
    ));
    salida.explicar(format!(
        "Se reconstruyeron {} puntos entre 0 y {ancho}. El primero es 0 por convención: las distancias dicen cuánto se separan los puntos, pero no dónde está la recta.",
        miles(puntos.len())
    ));
    let mut espejo: Vec<i64> = puntos.iter().map(|&punto| ancho - punto).collect();
    espejo.sort_unstable();
    if espejo == puntos {
        salida.explicar("El conjunto es simétrico: su espejo (ancho − x) es el mismo conjunto.");
    } else {
        let texto: Vec<String> = espejo.iter().map(|punto| punto.to_string()).collect();
        salida.explicar(format!(
            "Su espejo ({}) tiene exactamente las mismas distancias: el problema siempre tiene al menos dos soluciones, una al derecho y otra al revés, y la búsqueda con retroceso devuelve la primera que encuentra.",
            texto.join(", ")
        ));
    }
    salida.explicar(
        "Así se arma un mapa de restricción por digestión parcial: los puntos son los sitios de corte (y los extremos) y las distancias, los largos de todos los fragmentos posibles.",
    );
    Ok(salida)
}

fn opcional(valor: Option<usize>) -> String {
    valor.map_or_else(|| "–".to_string(), |numero| numero.to_string())
}

fn recortado(texto: &str, maximo: usize) -> String {
    if texto.chars().count() <= maximo {
        texto.to_string()
    } else {
        format!("{}…", texto.chars().take(maximo).collect::<String>())
    }
}

fn final_de_letras(letras: &[u8], maximo: usize) -> String {
    if letras.len() <= maximo {
        como_texto(letras)
    } else {
        format!("…{}", como_texto(&letras[letras.len() - maximo..]))
    }
}

fn masas_cortas(masas: &[usize], maximo: usize) -> String {
    let mut texto = masas
        .iter()
        .take(maximo)
        .map(|masa| masa.to_string())
        .collect::<Vec<_>>()
        .join(" ");
    if masas.len() > maximo {
        texto.push_str(" …");
    }
    texto
}

fn cada_cuanto(total: usize, deseados: usize) -> usize {
    (total / deseados.max(1)).max(1)
}

fn toca(indice: usize, detalle: usize, salto: usize) -> bool {
    indice < detalle || indice.is_multiple_of(salto)
}

fn ventana<T: Clone>(elementos: &[T], actual: usize, ancho: usize, puntos: T) -> (Vec<T>, usize) {
    if elementos.is_empty() {
        return (Vec::new(), 0);
    }
    let actual = actual.min(elementos.len() - 1);
    let inicio = actual
        .saturating_sub(ancho / 2)
        .min(elementos.len().saturating_sub(ancho));
    let fin = (inicio + ancho).min(elementos.len());
    let mut visibles = Vec::with_capacity(fin - inicio + 2);
    let mut relativo = actual - inicio;
    if inicio > 0 {
        visibles.push(puntos.clone());
        relativo += 1;
    }
    visibles.extend_from_slice(&elementos[inicio..fin]);
    if fin < elementos.len() {
        visibles.push(puntos);
    }
    (visibles, relativo)
}

fn las_ultimas(masas: &[usize], cuantas: usize) -> (Vec<String>, Vec<usize>) {
    let desde = masas.len().saturating_sub(cuantas);
    let elementos: Vec<String> = masas[desde..].iter().map(|masa| masa.to_string()).collect();
    let ultimo = elementos.len().saturating_sub(1);
    (elementos, vec![ultimo])
}

struct EstadoTraduccion {
    numero: Option<usize>,
    codon: String,
    letra: String,
    peptido: String,
}

fn pizarra_de_traduccion(
    arn: &[u8],
    estado: &EstadoTraduccion,
    marcas: &[(usize, usize, Resalte)],
) -> Pizarra {
    Pizarra::nueva()
        .variable("arn", format!("{} bases", miles(arn.len())))
        .variable("numero", opcional(estado.numero))
        .variable("codon", &estado.codon)
        .variable("letra", &estado.letra)
        .variable("peptido", &estado.peptido)
        .secuencia(arn, marcas)
}

fn pasos_traducir(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let leida = entrada.secuencia("arn")?;
    antibioticos::traducir_arn(&leida)?;
    Ok(narrar_traduccion(&leida)?.0)
}

fn narrar_traduccion(arn: &[u8]) -> Result<(Vec<Paso>, antibioticos::Traduccion), String> {
    let mostrado: Vec<u8> = arn
        .iter()
        .map(|&base| if base == b'T' { b'U' } else { base })
        .collect();
    let codones = arn.len() / 3;
    let salto = cada_cuanto(codones, 40);
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let mut peptido: Vec<u8> = Vec::new();
    let mut estado = EstadoTraduccion {
        numero: None,
        codon: "–".to_string(),
        letra: "–".to_string(),
        peptido: "(vacío)".to_string(),
    };
    let mut nota = format!(
        "Tu texto tiene {} bases: {} de 3 bases cada uno. La proteína empieza vacía y el código va a leer de tres en tres, de izquierda a derecha, hasta encontrar un codón de parada.",
        miles(arn.len()),
        contados(codones, "codón completo", "codones completos")
    );
    if !arn.len().is_multiple_of(3) {
        nota.push_str(&format!(
            " Al final sobran {} bases que no alcanzan para un codón.",
            arn.len() % 3
        ));
    }
    if arn.contains(&b'T') {
        nota.push_str(" Trae T, así que es ADN: la T cuenta como U.");
    }
    let marcas_iniciales = if codones > 0 {
        vec![(0, 3, Resalte::Actual)]
    } else {
        Vec::new()
    };
    narracion.agregar(
        pizarra_de_traduccion(&mostrado, &estado, &marcas_iniciales)
            .paso("let mut peptido = Vec::with_capacity(arn.len() / 3);", nota),
    );
    let mut parada = None;
    for (numero, codon) in arn.chunks_exact(3).enumerate() {
        let letra = aminoacido(codon);
        if letra == b'X' {
            return Err(format!(
                "El codón «{}» (posición {}, contando desde 0) tiene letras que no son A, C, G ni U.",
                como_texto(codon),
                numero * 3
            ));
        }
        let inicio = numero * 3;
        let detallado = numero < 3;
        let registrar = detallado || letra == b'*' || numero.is_multiple_of(salto);
        if letra != b'*' {
            peptido.push(letra);
        }
        if !registrar {
            continue;
        }
        let texto_codon = como_arn(codon);
        estado.numero = Some(numero);
        estado.codon = texto_codon.clone();
        let marcas = [(inicio, inicio + 3, Resalte::Actual)];
        if detallado {
            narracion.agregar(pizarra_de_traduccion(&mostrado, &estado, &marcas).paso(
                "for (numero, codon) in arn.chunks_exact(3).enumerate() {",
                format!(
                    "Toma el codón número {numero}: las bases {inicio} a {}, que dicen {texto_codon}.",
                    inicio + 2
                ),
            ));
        }
        if detallado || letra == b'*' {
            estado.letra = if letra == b'*' {
                "* (parada)".to_string()
            } else {
                (letra as char).to_string()
            };
            let nota = if letra == b'*' {
                format!(
                    "En la tabla del código genético, {texto_codon} no lleva aminoácido: es un codón de parada."
                )
            } else {
                format!(
                    "En la tabla del código genético, {texto_codon} codifica {} ({}).",
                    nombre_aminoacido(letra),
                    letra as char
                )
            };
            narracion.agregar(
                pizarra_de_traduccion(&mostrado, &estado, &marcas)
                    .paso("let letra = aminoacido(codon);", nota),
            );
        }
        if letra == b'*' {
            parada = Some(inicio);
            let marcas = [(inicio, inicio + 3, Resalte::Hallado)];
            narracion.agregar(pizarra_de_traduccion(&mostrado, &estado, &marcas).paso(
                "if letra == b'*' {",
                "Es un codón de parada (UAA, UAG o UGA): le dice al ribosoma que la proteína terminó. Aquí se deja de leer.",
            ));
            break;
        }
        estado.letra = (letra as char).to_string();
        estado.peptido = final_de_letras(&peptido, 30);
        let nota = if detallado {
            format!(
                "Agrego {} al final de la proteína, que ya lleva {}.",
                letra as char,
                contados(peptido.len(), "aminoácido", "aminoácidos")
            )
        } else {
            format!(
                "Codón {numero}: {texto_codon} da {}. La proteína ya lleva {}{}.",
                letra as char,
                contados(peptido.len(), "aminoácido", "aminoácidos"),
                if salto > 1 {
                    " (los codones que no te enseño se leen igual)"
                } else {
                    ""
                }
            )
        };
        narracion.agregar(
            pizarra_de_traduccion(&mostrado, &estado, &marcas).paso("peptido.push(letra);", nota),
        );
    }
    estado.peptido = final_de_letras(&peptido, 60);
    let traduccion = antibioticos::Traduccion { peptido, parada };
    let ultimo = match parada {
        Some(posicion) => pizarra_de_traduccion(
            &mostrado,
            &estado,
            &[(posicion, posicion + 3, Resalte::Hallado)],
        )
        .paso(
            "return Ok(Traduccion {",
            format!(
                "Resultado: {} (proteína de {}). La lectura se detuvo en la posición {}{}. En «Qué significa» te digo si parece una proteína completa.",
                estado.peptido,
                contados(
                    traduccion.peptido.len(),
                    "aminoácido",
                    "aminoácidos"
                ),
                miles(posicion),
                if arn.len() == posicion + 3 {
                    ", justo en el último codón".to_string()
                } else {
                    format!(
                        "; {} que siguen no se leen",
                        contados(arn.len() - posicion - 3, "base", "bases")
                    )
                }
            ),
        ),
        None => pizarra_de_traduccion(&mostrado, &estado, &[]).paso(
            "Ok(Traduccion {",
            format!(
                "Se acabaron los codones sin encontrar ninguna parada. Resultado: {} (proteína de {}). En «Qué significa» te cuento qué puede querer decir eso.",
                estado.peptido,
                contados(
                    traduccion.peptido.len(),
                    "aminoácido",
                    "aminoácidos"
                )
            ),
        ),
    };
    Ok((narracion.terminar(ultimo), traduccion))
}

struct EstadoCodificacion {
    ultimo_inicio: Option<usize>,
    inicio: Option<usize>,
    inicio_en_la_inversa: Option<usize>,
    encontradas: usize,
}

fn pizarra_de_codificacion(
    texto: &[u8],
    peptido: &[u8],
    largo: usize,
    estado: &EstadoCodificacion,
    marcas: &[(usize, usize, Resalte)],
) -> Pizarra {
    Pizarra::nueva()
        .variable("texto", format!("{} bases", miles(texto.len())))
        .variable("peptido", recortado(&como_texto(peptido), 24))
        .variable("largo", largo)
        .variable("ultimo_inicio", opcional(estado.ultimo_inicio))
        .variable("inicio", opcional(estado.inicio))
        .variable(
            "inicio_en_la_inversa",
            opcional(estado.inicio_en_la_inversa),
        )
        .variable("encontradas", estado.encontradas)
        .secuencia(texto, marcas)
}

fn aminoacidos_en_lista(aminoacidos: &[u8], cuantos: usize) -> Vec<String> {
    aminoacidos
        .iter()
        .take(cuantos)
        .enumerate()
        .map(|(posicion, &letra)| format!("{posicion}: {}", letra as char))
        .collect()
}

fn primer_fallo(aminoacidos: &[u8], inicio: usize, peptido: &[u8]) -> Option<usize> {
    for (indice, &letra) in peptido.iter().enumerate() {
        if aminoacidos[inicio + 3 * indice] != letra {
            return Some(indice);
        }
    }
    None
}

fn filas_de_comparacion(
    aminoacidos: &[u8],
    inicio: usize,
    peptido: &[u8],
    fallo: Option<usize>,
) -> Vec<Vec<String>> {
    let fin = fallo.map_or(peptido.len(), |indice| indice + 1);
    (fin.saturating_sub(6)..fin)
        .map(|indice| {
            let vista = aminoacidos[inicio + 3 * indice];
            vec![
                (peptido[indice] as char).to_string(),
                (vista as char).to_string(),
                if vista == peptido[indice] {
                    "sí"
                } else {
                    "no"
                }
                .to_string(),
            ]
        })
        .collect()
}

fn marcas_de_comparacion(
    inicio: usize,
    largo: usize,
    fallo: Option<usize>,
    cuantos: usize,
) -> Vec<(usize, usize, Resalte)> {
    let fin = fallo.map_or(cuantos, |indice| indice + 1);
    let mut marcas = vec![(inicio, inicio + largo, Resalte::Actual)];
    for indice in fin.saturating_sub(6)..fin {
        marcas.push((
            inicio + 3 * indice,
            inicio + 3 * indice + 3,
            Resalte::Comparado,
        ));
    }
    marcas
}

fn tripleta(texto: &[u8], desde: usize) -> String {
    let hasta = (desde + 3).min(texto.len());
    como_texto(&texto[desde.min(hasta)..hasta])
}

fn pasos_codificar(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let texto = entrada.adn("texto")?;
    let peptido = leer_peptido(entrada, "peptido", MAXIMO_PEPTIDO)?;
    antibioticos::codificaciones_del_peptido(&texto, &peptido)?;
    Ok(narrar_codificacion(&texto, &peptido).0)
}

fn narrar_codificacion(
    texto: &[u8],
    peptido: &[u8],
) -> (Vec<Paso>, Vec<antibioticos::Codificacion>) {
    let largo = 3 * peptido.len();
    let letras_del_peptido = como_texto(peptido);
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let mut estado = EstadoCodificacion {
        ultimo_inicio: None,
        inicio: None,
        inicio_en_la_inversa: None,
        encontradas: 0,
    };
    let mut encontradas: Vec<antibioticos::Codificacion> = Vec::new();
    narracion.agregar(
        pizarra_de_codificacion(texto, peptido, largo, &estado, &[]).paso(
            "let largo = 3 * peptido.len();",
            format!(
                "Buscas el péptido {} en un texto de {} bases. Cada aminoácido se escribe con un codón de 3 bases, así que el péptido ocupa {} bases en el ADN. Hay que mirar las dos hebras, porque el gen puede estar en cualquiera.",
                recortado(&letras_del_peptido, 24),
                miles(texto.len()),
                miles(largo)
            ),
        ),
    );
    if largo > texto.len() {
        let ultimo = pizarra_de_codificacion(texto, peptido, largo, &estado, &[]).paso(
            "return Ok(encontradas);",
            format!(
                "El péptido necesita {} bases y el texto solo tiene {}: no cabe en ninguna parte, así que no hay coincidencias.",
                miles(largo),
                miles(texto.len())
            ),
        );
        return (narracion.terminar(ultimo), encontradas);
    }
    let directa = texto.windows(3).map(aminoacido).collect::<Vec<u8>>();
    let marcas_de_codones = [
        (0, 3, Resalte::Actual),
        (1, 4, Resalte::Comparado),
        (2, 5, Resalte::Comparado),
    ];
    narracion.agregar(
        pizarra_de_codificacion(texto, peptido, largo, &estado, &marcas_de_codones)
            .lista(
                "Aminoácido que empieza en cada posición",
                aminoacidos_en_lista(&directa, 18),
                &[],
            )
            .paso(
                "let directa = aminoacido_en_cada_posicion(texto);",
                format!(
                    "Traduzco el texto desde CADA posición, no solo de 3 en 3: desde la 0 leo {}, desde la 1 leo {}, desde la 2 leo {}… Así los tres marcos de lectura quedan listos a la vez. Obtengo {} aminoácidos.",
                    tripleta(texto, 0),
                    tripleta(texto, 1),
                    tripleta(texto, 2),
                    miles(directa.len())
                ),
            ),
    );
    let inversa_texto = complemento_inverso(texto);
    let inversa = inversa_texto
        .windows(3)
        .map(aminoacido)
        .collect::<Vec<u8>>();
    narracion.agregar(
        pizarra_de_codificacion(&inversa_texto, peptido, largo, &estado, &marcas_de_codones)
            .lista(
                "Aminoácido que empieza en cada posición (hebra inversa)",
                aminoacidos_en_lista(&inversa, 18),
                &[],
            )
            .paso(
                "let inversa = aminoacido_en_cada_posicion(&complemento_inverso(texto));",
                "La otra hebra es el complemento inverso: se lee al revés y cambiando A con T y C con G. En la pizarra ves ese texto. Lo traduzco igual, posición por posición.",
            ),
    );
    let ultimo_inicio = texto.len() - largo;
    estado.ultimo_inicio = Some(ultimo_inicio);
    narracion.agregar(
        pizarra_de_codificacion(texto, peptido, largo, &estado, &[]).paso(
            "let ultimo_inicio = texto.len() - largo;",
            format!(
                "El tramo de {} bases puede empezar en las posiciones 0 a {} (después ya no cabe).",
                miles(largo),
                miles(ultimo_inicio)
            ),
        ),
    );
    let mut aciertos = 0;
    for inicio in 0..=ultimo_inicio {
        let fallo = primer_fallo(&directa, inicio, peptido);
        let acierto = fallo.is_none();
        let detallado = inicio < 3;
        if !detallado && !acierto {
            continue;
        }
        estado.inicio = Some(inicio);
        if acierto {
            aciertos += 1;
            encontradas.push(antibioticos::Codificacion {
                inicio,
                hebra: b'+',
                fragmento: texto[inicio..inicio + largo].to_vec(),
            });
            estado.encontradas = encontradas.len();
        }
        if detallado {
            let filas = filas_de_comparacion(&directa, inicio, peptido, fallo);
            let nota = match fallo {
                Some(indice) => format!(
                    "Desde la posición {inicio}: el aminoácido {} del péptido es {} pero el texto da {}. No coincide, así que ni sigo mirando.",
                    indice + 1,
                    peptido[indice] as char,
                    directa[inicio + 3 * indice] as char
                ),
                None => format!(
                    "Desde la posición {inicio}: todos los aminoácidos coinciden con {}.",
                    recortado(&letras_del_peptido, 24)
                ),
            };
            let resaltar: Vec<(usize, usize)> = filas
                .iter()
                .enumerate()
                .map(|(fila, _)| (fila, 2))
                .collect();
            narracion.agregar(
                pizarra_de_codificacion(
                    texto,
                    peptido,
                    largo,
                    &estado,
                    &marcas_de_comparacion(inicio, largo, fallo, peptido.len()),
                )
                .tabla(
                    "Péptido contra texto (hebra directa)",
                    &["Pide el péptido", "Da el texto", "¿Igual?"],
                    filas,
                    &resaltar,
                )
                .paso("if codifica_desde(&directa, inicio, peptido) {", nota),
            );
            if inicio == 2 && ultimo_inicio > 2 {
                narracion.agregar(
                    pizarra_de_codificacion(texto, peptido, largo, &estado, &[]).paso(
                        "for inicio in 0..=ultimo_inicio {",
                        format!(
                            "Con las otras {} posiciones se hace exactamente lo mismo. Desde aquí solo te enseño las que coinciden.",
                            miles(ultimo_inicio - 2)
                        ),
                    ),
                );
            }
        }
        if acierto && aciertos <= 6 {
            narracion.agregar(
                pizarra_de_codificacion(
                    texto,
                    peptido,
                    largo,
                    &estado,
                    &[(inicio, inicio + largo, Resalte::Hallado)],
                )
                .paso(
                    "encontradas.push(Codificacion {",
                    format!(
                        "¡Coincide desde la posición {inicio}! Los codones {} dan {}. Se guarda como hallazgo en la hebra +.",
                        en_codones(&texto[inicio..inicio + largo]),
                        recortado(&letras_del_peptido, 24)
                    ),
                ),
            );
        }
    }
    let en_la_directa = encontradas.len();
    estado.inicio = None;
    narracion.agregar(
        pizarra_de_codificacion(texto, peptido, largo, &estado, &[]).paso(
            "for inicio_en_la_inversa in 0..=ultimo_inicio {",
            format!(
                "Terminó la hebra directa: {} en {} posiciones. Ahora la inversa, con el mismo método.",
                contados(en_la_directa, "coincidencia", "coincidencias"),
                miles(ultimo_inicio + 1)
            ),
        ),
    );
    let mut aciertos = 0;
    for inicio_en_la_inversa in 0..=ultimo_inicio {
        let fallo = primer_fallo(&inversa, inicio_en_la_inversa, peptido);
        let acierto = fallo.is_none();
        let detallado = inicio_en_la_inversa < 3;
        if !detallado && !acierto {
            continue;
        }
        estado.inicio_en_la_inversa = Some(inicio_en_la_inversa);
        if detallado {
            let filas = filas_de_comparacion(&inversa, inicio_en_la_inversa, peptido, fallo);
            let nota = match fallo {
                Some(indice) => format!(
                    "Hebra inversa, posición {inicio_en_la_inversa}: el aminoácido {} del péptido es {} pero la hebra inversa da {}. No coincide.",
                    indice + 1,
                    peptido[indice] as char,
                    inversa[inicio_en_la_inversa + 3 * indice] as char
                ),
                None => format!(
                    "Hebra inversa, posición {inicio_en_la_inversa}: todos los aminoácidos coinciden con {}.",
                    recortado(&letras_del_peptido, 24)
                ),
            };
            let resaltar: Vec<(usize, usize)> = filas
                .iter()
                .enumerate()
                .map(|(fila, _)| (fila, 2))
                .collect();
            narracion.agregar(
                pizarra_de_codificacion(
                    &inversa_texto,
                    peptido,
                    largo,
                    &estado,
                    &marcas_de_comparacion(inicio_en_la_inversa, largo, fallo, peptido.len()),
                )
                .tabla(
                    "Péptido contra texto (hebra inversa)",
                    &["Pide el péptido", "Da el texto", "¿Igual?"],
                    filas,
                    &resaltar,
                )
                .paso(
                    "if codifica_desde(&inversa, inicio_en_la_inversa, peptido) {",
                    nota,
                ),
            );
        }
        if acierto {
            aciertos += 1;
            let inicio = ultimo_inicio - inicio_en_la_inversa;
            encontradas.push(antibioticos::Codificacion {
                inicio,
                hebra: b'-',
                fragmento: texto[inicio..inicio + largo].to_vec(),
            });
            estado.inicio = Some(inicio);
            estado.encontradas = encontradas.len();
            if aciertos <= 6 {
                let tramo = &texto[inicio..inicio + largo];
                narracion.agregar(
                    pizarra_de_codificacion(
                        texto,
                        peptido,
                        largo,
                        &estado,
                        &[(inicio, inicio + largo, Resalte::Hallado)],
                    )
                    .paso(
                        "let inicio = ultimo_inicio - inicio_en_la_inversa;",
                        format!(
                            "¡Coincide en la hebra inversa! Esa posición ({inicio_en_la_inversa}) de la hebra inversa es el tramo que empieza en {inicio} del texto original: {}. Su complemento inverso, {}, da {}.",
                            como_texto(tramo),
                            en_codones(&complemento_inverso(tramo)),
                            recortado(&letras_del_peptido, 24)
                        ),
                    ),
                );
            }
        }
        if inicio_en_la_inversa == 2 && ultimo_inicio > 2 {
            narracion.agregar(
                pizarra_de_codificacion(texto, peptido, largo, &estado, &[]).paso(
                    "for inicio_en_la_inversa in 0..=ultimo_inicio {",
                    "Las demás posiciones de la hebra inversa se revisan igual. Desde aquí solo te enseño las que coinciden.",
                ),
            );
        }
    }
    encontradas.sort_by_key(|codificacion| codificacion.inicio);
    estado.inicio = None;
    estado.inicio_en_la_inversa = None;
    estado.encontradas = encontradas.len();
    narracion.agregar(
        pizarra_de_codificacion(texto, peptido, largo, &estado, &[]).paso(
            "encontradas.sort_by_key(|codificacion| codificacion.inicio);",
            format!(
                "Junto lo de las dos hebras ({} en la directa y {} en la inversa) y lo ordeno por posición en el texto.",
                en_la_directa,
                encontradas.len() - en_la_directa
            ),
        ),
    );
    let marcas: Vec<(usize, usize, Resalte)> = encontradas
        .iter()
        .map(|c| (c.inicio, c.inicio + largo, Resalte::Hallado))
        .collect();
    let posiciones: Vec<String> = encontradas
        .iter()
        .take(12)
        .map(|c| format!("{} ({})", c.inicio, c.hebra as char))
        .collect();
    let ultimo = pizarra_de_codificacion(texto, peptido, largo, &estado, &marcas)
        .lista("Dónde está (posición y hebra)", posiciones, &[])
        .paso(
            "Ok(encontradas)",
            format!(
                "Resultado: {} de {}. Mira «Qué significa» para saber si es una casualidad o un gen de verdad.",
                contados(encontradas.len(), "coincidencia", "coincidencias"),
                recortado(&letras_del_peptido, 24)
            ),
        );
    (narracion.terminar(ultimo), encontradas)
}

fn atomos_en_fila(nombre: &str, formula: &antibioticos::Formula) -> Vec<String> {
    vec![
        nombre.to_string(),
        formula.carbono.to_string(),
        formula.hidrogeno.to_string(),
        formula.nitrogeno.to_string(),
        formula.oxigeno.to_string(),
        formula.azufre.to_string(),
    ]
}

fn sumar_formulas(total: &mut antibioticos::Formula, residuo: &antibioticos::Formula) {
    total.carbono += residuo.carbono;
    total.hidrogeno += residuo.hidrogeno;
    total.nitrogeno += residuo.nitrogeno;
    total.oxigeno += residuo.oxigeno;
    total.azufre += residuo.azufre;
}

fn pizarra_de_masa(
    ciclico: bool,
    total: &antibioticos::Formula,
    letra: Option<u8>,
    residuo: Option<&antibioticos::Formula>,
) -> Pizarra {
    let mut filas = Vec::new();
    if let Some(residuo) = residuo {
        filas.push(atomos_en_fila("residuo", residuo));
    }
    filas.push(atomos_en_fila("total", total));
    let ultima = filas.len() - 1;
    let resaltar: Vec<(usize, usize)> = match residuo {
        Some(residuo) => [
            residuo.carbono,
            residuo.hidrogeno,
            residuo.nitrogeno,
            residuo.oxigeno,
            residuo.azufre,
        ]
        .iter()
        .enumerate()
        .filter(|(_, cuantos)| **cuantos > 0)
        .map(|(columna, _)| (ultima, columna + 1))
        .collect(),
        None => Vec::new(),
    };
    Pizarra::nueva()
        .variable("ciclico", ciclico)
        .variable("total", texto_de_formula(total))
        .variable(
            "letra",
            letra.map_or_else(|| "–".to_string(), |letra| (letra as char).to_string()),
        )
        .variable(
            "residuo",
            residuo.map_or_else(|| "–".to_string(), texto_de_formula),
        )
        .tabla("Átomos", &["", "C", "H", "N", "O", "S"], filas, &resaltar)
}

fn pasos_masa(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let peptido = leer_peptido(entrada, "peptido", MAXIMO_PEPTIDO)?;
    let forma = entrada.opcion("forma", FORMAS)?;
    let ciclico = forma == "cíclico";
    antibioticos::masa_monoisotopica(&peptido, ciclico)?;
    Ok(narrar_masa(&peptido, ciclico)?.0)
}

fn narrar_masa(
    peptido: &[u8],
    ciclico: bool,
) -> Result<(Vec<Paso>, antibioticos::Formula, f64), String> {
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let vacia = antibioticos::Formula {
        carbono: 0,
        hidrogeno: 0,
        nitrogeno: 0,
        oxigeno: 0,
        azufre: 0,
    };
    let mut total = if ciclico { vacia } else { antibioticos::AGUA };
    let nota = if ciclico {
        "Un péptido cíclico es un anillo: no tiene extremos libres, así que no sobra ninguna molécula de agua. El total empieza sin ningún átomo."
    } else {
        "Un péptido lineal tiene dos extremos libres (uno con H y otro con OH), que juntos son una molécula de agua, H2O. El total empieza con esos átomos."
    };
    narracion.agregar(pizarra_de_masa(ciclico, &total, None, None).paso(
        "let mut total = if ciclico {",
        format!(
            "Péptido {} de {}. {nota}",
            if ciclico { "cíclico" } else { "lineal" },
            contados(peptido.len(), "aminoácido", "aminoácidos")
        ),
    ));
    let cada = cada_cuanto(peptido.len(), 14);
    for (indice, &letra) in peptido.iter().enumerate() {
        let residuo = antibioticos::formula_del_peptido(&[letra], true)?;
        let detallado = indice < 3;
        if detallado {
            narracion.agregar(
                pizarra_de_masa(ciclico, &total, Some(letra), Some(&residuo)).paso(
                    "let residuo = formula_de_residuo(letra).ok_or_else(|| no_es_aminoacido(letra))?;",
                    format!(
                        "Aminoácido {} ({}, {}): en la tabla, al unirse a una cadena aporta {} (el residuo: el aminoácido sin su agua).",
                        indice + 1,
                        letra as char,
                        nombre_aminoacido(letra),
                        texto_de_formula(&residuo)
                    ),
                ),
            );
        }
        let antes = texto_de_formula(&total);
        sumar_formulas(&mut total, &residuo);
        if detallado || indice == peptido.len() - 1 || indice.is_multiple_of(cada) {
            let nota = if detallado {
                format!(
                    "Sumo, átomo por átomo (C, H, N, O, S), el residuo al total: {antes} se vuelve {}.",
                    texto_de_formula(&total)
                )
            } else {
                format!(
                    "Aminoácido {} ({}): aporta {}. Con los {} que llevo, el total es {}. Los que no te enseño se suman igual.",
                    miles(indice + 1),
                    letra as char,
                    texto_de_formula(&residuo),
                    miles(indice + 1),
                    texto_de_formula(&total)
                )
            };
            narracion.agregar(
                pizarra_de_masa(ciclico, &total, Some(letra), Some(&residuo))
                    .paso("total.carbono += residuo.carbono;", nota),
            );
        }
    }
    let formula = total;
    narracion.agregar(
        pizarra_de_masa(ciclico, &formula, None, None).paso(
            "Ok(total)",
            format!(
                "Ya sumé los {}: la fórmula del péptido es {}. Ahora hay que convertir cada átomo en masa.",
                contados(peptido.len(), "aminoácido", "aminoácidos"),
                texto_de_formula(&formula)
            ),
        ),
    );
    let terminos = [
        (
            "carbono (C)",
            formula.carbono,
            antibioticos::MASA_CARBONO_12,
            "Ok(atomos.carbono as f64 * MASA_CARBONO_12",
        ),
        (
            "hidrógeno (H)",
            formula.hidrogeno,
            antibioticos::MASA_HIDROGENO_1,
            "+ atomos.hidrogeno as f64 * MASA_HIDROGENO_1",
        ),
        (
            "nitrógeno (N)",
            formula.nitrogeno,
            antibioticos::MASA_NITROGENO_14,
            "+ atomos.nitrogeno as f64 * MASA_NITROGENO_14",
        ),
        (
            "oxígeno (O)",
            formula.oxigeno,
            antibioticos::MASA_OXIGENO_16,
            "+ atomos.oxigeno as f64 * MASA_OXIGENO_16",
        ),
        (
            "azufre (S)",
            formula.azufre,
            antibioticos::MASA_AZUFRE_32,
            "+ atomos.azufre as f64 * MASA_AZUFRE_32)",
        ),
    ];
    let mut suma = 0.0;
    let mut filas: Vec<Vec<String>> = terminos
        .iter()
        .map(|(nombre, cuantos, unitaria, _)| {
            vec![
                nombre.to_string(),
                cuantos.to_string(),
                decimal(*unitaria, 5),
                "·".to_string(),
            ]
        })
        .collect();
    let mut ultimo = None;
    for (posicion, (nombre, cuantos, unitaria, linea)) in terminos.iter().enumerate() {
        let aporte = *cuantos as f64 * unitaria;
        suma += aporte;
        filas[posicion][3] = decimal(aporte, 4);
        let nota = if posicion + 1 == terminos.len() {
            format!(
                "{nombre}: {cuantos} × {} = {} Da. Suma final: {} Da. Es la masa monoisotópica: la del péptido hecho solo con los isótopos más ligeros y abundantes (C-12, H-1, N-14, O-16, S-32), la que verías en un espectrómetro de alta resolución. En «Qué significa» está lo que quiere decir.",
                decimal(*unitaria, 5),
                decimal(aporte, 4),
                decimal(suma, 4)
            )
        } else {
            format!(
                "{nombre}: {cuantos} átomos × {} Da = {} Da. Llevo sumado {} Da.",
                decimal(*unitaria, 5),
                decimal(aporte, 4),
                decimal(suma, 4)
            )
        };
        let paso = Pizarra::nueva()
            .variable("ciclico", ciclico)
            .variable("atomos", texto_de_formula(&formula))
            .tabla(
                "Masa de cada tipo de átomo",
                &["Átomo", "Cantidad", "Masa de uno (Da)", "Aporte (Da)"],
                filas.clone(),
                &[(posicion, 3)],
            )
            .paso(linea, nota);
        if posicion + 1 == terminos.len() {
            ultimo = Some(paso);
        } else {
            narracion.agregar(paso);
        }
    }
    let ultimo = ultimo.ok_or_else(|| "Sin pasos finales.".to_string())?;
    Ok((narracion.terminar(ultimo), formula, suma))
}

struct EstadoEspectro {
    total: Option<usize>,
    n: Option<usize>,
    inicio: Option<usize>,
    suma: Option<usize>,
    fin: Option<usize>,
    masa: Option<usize>,
    cuantas: usize,
}

fn pizarra_de_espectro(
    letras: &[u8],
    masas: &[usize],
    ciclico: bool,
    estado: &EstadoEspectro,
    marcas: &[(usize, usize, Resalte)],
) -> Pizarra {
    Pizarra::nueva()
        .variable("peptido", format!("[{}]", masas_cortas(masas, 8)))
        .variable("ciclico", ciclico)
        .variable("masas", format!("{} masas", miles(estado.cuantas)))
        .variable("total", opcional(estado.total))
        .variable("n", opcional(estado.n))
        .variable("inicio", opcional(estado.inicio))
        .variable("suma", opcional(estado.suma))
        .variable("fin", opcional(estado.fin))
        .variable("masa", opcional(estado.masa))
        .secuencia(letras, marcas)
}

fn pasos_espectro_ciclico(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    pasos_de_espectro(entrada, true)
}

fn pasos_espectro_lineal(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    pasos_de_espectro(entrada, false)
}

fn pasos_de_espectro(entrada: &Entrada, ciclico: bool) -> Result<Vec<Paso>, String> {
    let peptido = leer_peptido(entrada, "peptido", MAXIMO_PEPTIDO_ESPECTRO)?;
    let masas = antibioticos::masas_del_peptido(&peptido)?;
    Ok(narrar_espectro(&peptido, &masas, ciclico).0)
}

fn narrar_espectro(letras: &[u8], masas: &[usize], ciclico: bool) -> (Vec<Paso>, Vec<usize>) {
    let n = masas.len();
    let total = antibioticos::masa_total(masas);
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let mut espectro: Vec<usize> = Vec::new();
    let mut estado = EstadoEspectro {
        total: None,
        n: None,
        inicio: None,
        suma: None,
        fin: None,
        masa: None,
        cuantas: 0,
    };
    let desglose: Vec<String> = letras
        .iter()
        .zip(masas)
        .take(8)
        .map(|(&letra, masa)| format!("{} = {masa}", letra as char))
        .collect();
    let linea_inicial = if ciclico {
        "masas_de_subpeptidos(peptido, true, &mut masas);"
    } else {
        "masas_de_subpeptidos(peptido, false, &mut masas);"
    };
    narracion.agregar(
        pizarra_de_espectro(letras, masas, ciclico, &estado, &[]).paso(
            linea_inicial,
            format!(
                "Cada aminoácido pesa una masa entera ({}{}). El espectro se arma sumando las masas de todos los pedazos seguidos del péptido{}.",
                desglose.join(", "),
                if n > 8 { ", …" } else { "" },
                if ciclico {
                    ", y también de los que dan la vuelta al anillo"
                } else {
                    ""
                }
            ),
        ),
    );
    espectro.push(0);
    estado.cuantas = 1;
    let (elementos, resaltar) = las_ultimas(&espectro, 10);
    narracion.agregar(
        pizarra_de_espectro(letras, masas, ciclico, &estado, &[])
            .lista("Masas anotadas (las últimas)", elementos, &resaltar)
            .paso(
                "masas.push(0);",
                "Primero anoto 0: es el «pedazo vacío», que todo espectro incluye.",
            ),
    );
    estado.total = Some(total);
    narracion.agregar(
        pizarra_de_espectro(letras, masas, ciclico, &estado, &[]).paso(
            "let total = masa_total(peptido);",
            format!(
                "La masa del péptido entero es la suma de todas sus masas: {} = {}.",
                masas_cortas(masas, 8).replace(' ', " + "),
                miles(total)
            ),
        ),
    );
    estado.n = Some(n);
    narracion.agregar(
        pizarra_de_espectro(letras, masas, ciclico, &estado, &[]).paso(
            "let n = peptido.len();",
            format!(
                "El péptido tiene {}.",
                contados(n, "aminoácido", "aminoácidos")
            ),
        ),
    );
    let cada = cada_cuanto(n, 25);
    for inicio in 0..n {
        let detallado = inicio < 2;
        estado.inicio = Some(inicio);
        estado.fin = None;
        estado.masa = None;
        let mut suma = 0;
        estado.suma = Some(0);
        if detallado {
            let (elementos, resaltar) = las_ultimas(&espectro, 10);
            narracion.agregar(
                pizarra_de_espectro(
                    letras,
                    masas,
                    ciclico,
                    &estado,
                    &[(inicio, inicio + 1, Resalte::Actual)],
                )
                .lista("Masas anotadas (las últimas)", elementos, &resaltar)
                .paso(
                    "for inicio in 0..n {",
                    format!(
                        "Empiezo en el aminoácido {} ({}) con la suma en 0 y voy agregando los que siguen, de uno en uno.",
                        inicio + 1,
                        letras[inicio] as char
                    ),
                ),
            );
        }
        let antes = espectro.len();
        for (fin, &masa) in masas.iter().enumerate().skip(inicio) {
            suma += masa;
            espectro.push(suma);
            let complemento = ciclico && inicio > 0 && fin + 1 < n;
            if detallado && fin - inicio < 5 {
                estado.fin = Some(fin);
                estado.masa = Some(masa);
                estado.suma = Some(suma);
                estado.cuantas = espectro.len();
                let (elementos, resaltar) = las_ultimas(&espectro, 10);
                narracion.agregar(
                    pizarra_de_espectro(
                        letras,
                        masas,
                        ciclico,
                        &estado,
                        &[(inicio, fin + 1, Resalte::Actual)],
                    )
                    .lista("Masas anotadas (las últimas)", elementos, &resaltar)
                    .paso(
                        "masas.push(suma);",
                        format!(
                            "Sumo el aminoácido {} (masa {masa}): el pedazo {} pesa {suma}. Lo anoto.",
                            fin + 1,
                            recortado(&como_texto(&letras[inicio..=fin]), 14)
                        ),
                    ),
                );
            }
            if complemento {
                espectro.push(total - suma);
                if detallado && fin - inicio < 5 {
                    estado.cuantas = espectro.len();
                    let (elementos, resaltar) = las_ultimas(&espectro, 10);
                    narracion.agregar(
                        pizarra_de_espectro(
                            letras,
                            masas,
                            ciclico,
                            &estado,
                            &[
                                (inicio, fin + 1, Resalte::Actual),
                                (fin + 1, n, Resalte::Comparado),
                                (0, inicio, Resalte::Comparado),
                            ],
                        )
                        .lista("Masas anotadas (las últimas)", elementos, &resaltar)
                        .paso(
                            "masas.push(total - suma);",
                            format!(
                                "En un anillo, lo que NO está en ese pedazo también es un pedazo (el resto que da la vuelta): pesa {} − {suma} = {}. Lo anoto también.",
                                miles(total),
                                total - suma
                            ),
                        ),
                    );
                }
            }
        }
        let anotadas = espectro.len() - antes;
        let faltaron = detallado && n - inicio > 5;
        if faltaron || (!detallado && (inicio == n - 1 || inicio.is_multiple_of(cada))) {
            estado.fin = Some(n - 1);
            estado.masa = Some(masas[n - 1]);
            estado.suma = Some(suma);
            estado.cuantas = espectro.len();
            let (elementos, resaltar) = las_ultimas(&espectro, 10);
            narracion.agregar(
                pizarra_de_espectro(
                    letras,
                    masas,
                    ciclico,
                    &estado,
                    &[(inicio, n, Resalte::Actual)],
                )
                .lista("Masas anotadas (las últimas)", elementos, &resaltar)
                .paso(
                    "masas.push(suma);",
                    format!(
                        "Empezando en el aminoácido {} ({}) sumo hasta el final: anoté {}. La última, {suma}, es la del pedazo que llega hasta el último aminoácido. Los pedazos que no te enseño se calculan igual.",
                        inicio + 1,
                        letras[inicio] as char,
                        contados(anotadas, "masa", "masas")
                    ),
                ),
            );
        }
    }
    estado.cuantas = espectro.len();
    espectro.sort_unstable();
    let esperadas = if ciclico {
        n * (n - 1) + 2
    } else {
        n * (n + 1) / 2 + 1
    };
    let mut mostradas: Vec<String> = espectro
        .iter()
        .take(40)
        .map(|masa| masa.to_string())
        .collect();
    if espectro.len() > 40 {
        mostradas.push("…".to_string());
    }
    let ultimo = pizarra_de_espectro(letras, masas, ciclico, &estado, &[])
        .lista("Espectro ordenado", mostradas, &[])
        .paso(
            "masas.sort_unstable();",
            format!(
                "Ordeno las masas de menor a mayor y listo: el espectro {} tiene {} masas ({} = {}), de 0 a {}. En «Qué significa» está cómo se usa.",
                if ciclico { "cíclico" } else { "lineal" },
                miles(espectro.len()),
                if ciclico { "n·(n−1) + 2" } else { "n·(n+1)/2 + 1" },
                miles(esperadas),
                miles(total)
            ),
        );
    (narracion.terminar(ultimo), espectro)
}

struct EstadoConteo {
    total: Option<usize>,
    valor: Option<u128>,
}

fn pizarra_de_conteo(masa: usize, estado: &EstadoConteo, no_nulos: &[String]) -> Pizarra {
    let desde = no_nulos.len().saturating_sub(8);
    let visibles = no_nulos[desde..].to_vec();
    let ultimo = visibles.len().saturating_sub(1);
    Pizarra::nueva()
        .variable("masa", miles(masa))
        .variable("total", opcional(estado.total))
        .variable(
            "cuantos[total]",
            estado.valor.map_or_else(|| "–".to_string(), miles_grande),
        )
        .lista(
            "Masas que ya tienen péptidos (las últimas): masa → cuántos",
            visibles,
            &[ultimo],
        )
}

fn pasos_contar(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let masa = entrada.natural("masa")?;
    antibioticos::peptidos_por_masa(masa)?;
    Ok(narrar_conteo(masa).0)
}

fn narrar_conteo(masa: usize) -> (Vec<Paso>, Vec<u128>) {
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let mut estado = EstadoConteo {
        total: None,
        valor: None,
    };
    let mut no_nulos: Vec<String> = Vec::new();
    let mut cuantos = vec![0u128; masa + 1];
    narracion.agregar(
        Pizarra::nueva()
            .variable("masa", miles(masa))
            .variable("total", "–")
            .variable("cuantos[total]", "–")
            .lista(
                "Las 18 masas de aminoácido distintas (Da)",
                MASAS_DISTINTAS.iter().map(|masa| masa.to_string()).collect(),
                &[],
            )
            .paso(
                "let mut cuantos = vec![0u128; masa + 1];",
                format!(
                    "Quieres saber cuántos péptidos (cadenas ordenadas de aminoácidos) pesan exactamente {} Da. Armo una tabla: cuantos[t] será el número de péptidos que pesan t. De momento todo vale 0. Solo hay 18 masas distintas de aminoácido, porque I y L pesan igual (113), y también K y Q (128).",
                    miles(masa)
                ),
            ),
    );
    cuantos[0] = 1;
    no_nulos.push("0 → 1".to_string());
    estado.total = Some(0);
    estado.valor = Some(1);
    narracion.agregar(
        pizarra_de_conteo(masa, &estado, &no_nulos).paso(
            "cuantos[0] = 1;",
            "Hay exactamente una manera de pesar 0 Da: no poner ningún aminoácido (el péptido vacío). Todo el conteo se apoya en ese 1.",
        ),
    );
    if masa >= 1 {
        narracion.agregar(
            pizarra_de_conteo(masa, &estado, &no_nulos).paso(
                "for total in 1..=masa {",
                format!(
                    "Lleno la tabla de 1 hasta {}. Para cada total pruebo cada masa de aminoácido como ÚLTIMA pieza: si cabe, todos los péptidos que pesan «total − última» pueden terminar con ella, y se suman. Hasta 56 Da ninguna pieza cabe (la más ligera pesa 57), así que esos totales se quedan en 0.",
                    miles(masa)
                ),
            ),
        );
    }
    let cada = cada_cuanto(masa, 22);
    let mut mostrados = 0;
    for total in 1..=masa {
        let mut valor = 0u128;
        for ultima in MASAS_DISTINTAS {
            if ultima <= total {
                valor += cuantos[total - ultima];
            }
        }
        cuantos[total] = valor;
        if valor == 0 {
            continue;
        }
        no_nulos.push(format!("{} → {}", miles(total), miles_grande(valor)));
        let registrar = mostrados < 8 || total.is_multiple_of(cada) || total == masa;
        if !registrar {
            continue;
        }
        mostrados += 1;
        estado.total = Some(total);
        estado.valor = Some(valor);
        let mut filas = Vec::new();
        let mut resaltar = Vec::new();
        let mut aportes = Vec::new();
        for ultima in MASAS_DISTINTAS {
            if ultima <= total {
                let resto = total - ultima;
                filas.push(vec![
                    masa_con_letras(ultima),
                    miles(resto),
                    miles_grande(cuantos[resto]),
                ]);
                if cuantos[resto] > 0 {
                    resaltar.push((filas.len() - 1, 2));
                    aportes.push(format!("{ultima} aporta {}", miles_grande(cuantos[resto])));
                }
            }
        }
        narracion.agregar(
            pizarra_de_conteo(masa, &estado, &no_nulos)
                .tabla(
                    &format!("Para total = {}: última pieza posible", miles(total)),
                    &["Última pieza", "Queda (total − última)", "cuantos[queda]"],
                    filas,
                    &resaltar,
                )
                .paso(
                    "cuantos[total] += cuantos[total - ultima];",
                    format!(
                        "Para {} Da: {}. Sumo todo y cuantos[{}] = {}.",
                        miles(total),
                        recortado(&aportes.join(", "), 120),
                        miles(total),
                        miles_grande(valor)
                    ),
                ),
        );
    }
    estado.total = Some(masa);
    estado.valor = Some(cuantos[masa]);
    let nota = if cuantos[masa] == 0 {
        format!(
            "Ninguna combinación de aminoácidos pesa exactamente {} Da, así que la respuesta es 0 (cuantos[{}] = 0).",
            miles(masa),
            miles(masa)
        )
    } else {
        format!(
            "La tabla está llena. La respuesta es cuantos[{}] = {}: ese es el número de péptidos que pesan exactamente {} Da. En «Qué significa» lo comparo con cosas que puedas imaginar.",
            miles(masa),
            miles_grande(cuantos[masa]),
            miles(masa)
        )
    };
    let ultimo = pizarra_de_conteo(masa, &estado, &no_nulos).paso("Ok(cuantos)", nota);
    (narracion.terminar(ultimo), cuantos)
}

fn pasos_puntaje_ciclico(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    pasos_de_puntaje(entrada, true)
}

fn pasos_puntaje_lineal(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    pasos_de_puntaje(entrada, false)
}

fn pasos_de_puntaje(entrada: &Entrada, ciclico: bool) -> Result<Vec<Paso>, String> {
    let peptido = leer_peptido(entrada, "peptido", MAXIMO_PEPTIDO_ESPECTRO)?;
    let masas = antibioticos::masas_del_peptido(&peptido)?;
    let espectro = leer_espectro(entrada, "espectro")?;
    if ciclico {
        antibioticos::puntaje_ciclico(&masas, &espectro)?;
    } else {
        antibioticos::puntaje_lineal(&masas, &espectro)?;
    }
    Ok(narrar_puntaje(&peptido, &masas, &espectro, ciclico)?.0)
}

struct EstadoPuntaje {
    masa: Option<usize>,
    compartidas: usize,
    conteo: Option<usize>,
    usados: Option<u32>,
}

fn pizarra_de_puntaje(
    letras: &[u8],
    ciclico: bool,
    experimentales: usize,
    teoricas: usize,
    estado: &EstadoPuntaje,
) -> Pizarra {
    Pizarra::nueva()
        .variable("peptido", recortado(&como_texto(letras), 24))
        .variable("ciclico", ciclico)
        .variable("espectro", format!("{} masas", miles(experimentales)))
        .variable("teorico", format!("{} masas", miles(teoricas)))
        .variable("masa", opcional(estado.masa))
        .variable("conteo[masa]", opcional(estado.conteo))
        .variable(
            "usados[masa]",
            estado
                .usados
                .map_or_else(|| "–".to_string(), |usados| usados.to_string()),
        )
        .variable("compartidas", estado.compartidas)
}

fn narrar_puntaje(
    letras: &[u8],
    masas: &[usize],
    espectro: &[usize],
    ciclico: bool,
) -> Result<(Vec<Paso>, usize), String> {
    let experimental = antibioticos::Espectro::nuevo(espectro)?;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let mut estado = EstadoPuntaje {
        masa: None,
        compartidas: 0,
        conteo: None,
        usados: None,
    };
    let mut distintas_del_experimental: Vec<(usize, usize)> = Vec::new();
    for grupo in experimental.masas.chunk_by(|a, b| a == b) {
        distintas_del_experimental.push((grupo[0], grupo.len()));
    }
    let filas_iniciales: Vec<Vec<String>> = distintas_del_experimental
        .iter()
        .take(14)
        .map(|(masa, veces)| vec![masa.to_string(), veces.to_string()])
        .collect();
    let mut nota = format!(
        "Tienes el péptido {} y un espectro medido de {}, de 0 a {} Da. Voy a contar cuántas masas comparten el espectro teórico del péptido y el medido. Primero ordeno el espectro medido y anoto cuántas veces aparece cada masa (conteo).",
        recortado(&como_texto(letras), 24),
        contados(experimental.masas.len(), "masa", "masas"),
        miles(experimental.masa_madre)
    );
    if distintas_del_experimental.len() < experimental.masas.len() {
        nota.push_str(" Hay masas repetidas: cada copia solo puede emparejarse una vez.");
    }
    narracion.agregar(
        pizarra_de_puntaje(letras, ciclico, experimental.masas.len(), 0, &estado)
            .tabla(
                "Espectro medido: cuántas veces aparece cada masa",
                &["Masa (Da)", "Veces"],
                filas_iniciales,
                &[],
            )
            .paso("let mut experimental = Espectro::nuevo(espectro)?;", nota),
    );
    let mut teorico = Vec::new();
    antibioticos::masas_de_subpeptidos(masas, ciclico, &mut teorico);
    let mut distintas: Vec<usize> = teorico
        .iter()
        .copied()
        .chain(experimental.masas.iter().copied())
        .collect();
    distintas.sort_unstable();
    distintas.dedup();
    let mut usados: Vec<u32> = vec![0; experimental.masa_madre + 1];
    narracion.agregar(
        pizarra_de_puntaje(
            letras,
            ciclico,
            experimental.masas.len(),
            teorico.len(),
            &estado,
        )
        .lista(
            "Espectro teórico (en el orden en que lo calcula el código)",
            teorico.iter().take(30).map(|masa| masa.to_string()).collect(),
            &[],
        )
        .paso(
            "masas_de_subpeptidos(peptido, ciclico, &mut self.teorico);",
            format!(
                "Calculo el espectro teórico del péptido: las masas de todos sus pedazos seguidos{}. Salen {}. El orden no importa para contar.",
                if ciclico {
                    ", incluidos los que dan la vuelta al anillo"
                } else {
                    " (sin dar la vuelta)"
                },
                contados(teorico.len(), "masa", "masas")
            ),
        ),
    );
    let cada = cada_cuanto(teorico.len(), 20);
    let mut compartidas = 0;
    for (indice, &masa) in teorico.iter().enumerate() {
        let veces = experimental.veces(masa);
        let usadas = if masa < usados.len() { usados[masa] } else { 0 };
        let coincide = masa < usados.len() && (usadas as usize) < veces;
        if coincide {
            usados[masa] += 1;
            compartidas += 1;
        }
        if !toca(indice, 8, cada) {
            continue;
        }
        estado.masa = Some(masa);
        estado.compartidas = compartidas;
        estado.conteo = Some(veces);
        estado.usados = Some(if coincide { usadas + 1 } else { usadas });
        let posicion = distintas.binary_search(&masa).unwrap_or(0);
        let desde = posicion.saturating_sub(5);
        let hasta = (desde + 11).min(distintas.len());
        let filas: Vec<Vec<String>> = distintas[desde..hasta]
            .iter()
            .map(|otra| {
                let tiene = experimental.veces(*otra);
                let gastadas = if *otra < usados.len() {
                    usados[*otra]
                } else {
                    0
                };
                vec![otra.to_string(), tiene.to_string(), gastadas.to_string()]
            })
            .collect();
        let fila_actual = posicion - desde;
        let (nota, linea) = if coincide {
            (
                format!(
                    "Masa {masa} (la {} del teórico): en el medido aparece {} y aún no he gastado ninguna de esas copias, así que coincide. Gasto una copia (usados[{masa}] pasa a {}) y el puntaje sube a {compartidas}.",
                    indice + 1,
                    contados(veces, "vez", "veces"),
                    usadas + 1
                ),
                "compartidas += 1;",
            )
        } else if veces == 0 {
            (
                format!(
                    "Masa {masa} (la {} del teórico): el espectro medido no la tiene, así que no suma. El puntaje se queda en {compartidas}.",
                    indice + 1
                ),
                "if masa < self.conteo.len() && self.usados[masa] < self.conteo[masa] {",
            )
        } else {
            (
                format!(
                    "Masa {masa} (la {} del teórico): el medido la tiene {} pero ya gasté todas con masas teóricas anteriores iguales, así que no queda con quién emparejarla. El puntaje se queda en {compartidas}.",
                    indice + 1,
                    contados(veces, "vez", "veces")
                ),
                "if masa < self.conteo.len() && self.usados[masa] < self.conteo[masa] {",
            )
        };
        narracion.agregar(
            pizarra_de_puntaje(
                letras,
                ciclico,
                experimental.masas.len(),
                teorico.len(),
                &estado,
            )
            .tabla(
                "Masas cercanas: medido y copias ya usadas",
                &["Masa", "En el medido", "Usadas"],
                filas,
                &[(fila_actual, 1), (fila_actual, 2)],
            )
            .paso(linea, nota),
        );
    }
    estado.masa = None;
    estado.conteo = None;
    estado.usados = None;
    estado.compartidas = compartidas;
    narracion.agregar(
        pizarra_de_puntaje(
            letras,
            ciclico,
            experimental.masas.len(),
            teorico.len(),
            &estado,
        )
        .paso(
            "self.usados[masa] = 0;",
            "Con el conteo hecho, pongo en cero el registro de copias usadas: es solo limpieza para que la próxima comparación empiece de cero.",
        ),
    );
    let ultimo = pizarra_de_puntaje(
        letras,
        ciclico,
        experimental.masas.len(),
        teorico.len(),
        &estado,
    )
    .paso(
        if ciclico {
            "Ok(experimental.compartidas(peptido, true))"
        } else {
            "Ok(experimental.compartidas(peptido, false))"
        },
        format!(
            "Puntaje final: {} de las {} masas del espectro teórico (y de las {} del medido) coinciden. En «Qué significa» te digo cómo leer ese número.",
            miles(compartidas),
            miles(teorico.len()),
            miles(experimental.masas.len())
        ),
    );
    Ok((narracion.terminar(ultimo), compartidas))
}

fn pasos_recortar(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let texto = entrada.texto("peptidos")?;
    let peptidos = leer_lista_de_peptidos(&texto)?;
    let espectro = leer_espectro(entrada, "espectro")?;
    let n = entrada.natural("n")?;
    if n == 0 {
        return Err("N debe ser al menos 1.".to_string());
    }
    let tabla: Vec<Vec<usize>> = peptidos.iter().map(|(_, masas)| masas.clone()).collect();
    antibioticos::recortar(tabla, &espectro, n)?;
    Ok(narrar_recorte(&peptidos, &espectro, n)?.0)
}

fn pizarra_de_recorte(
    n: usize,
    peptidos: usize,
    peptido: Option<&str>,
    puntajes: usize,
    corte: Option<usize>,
    conservados: usize,
) -> Pizarra {
    Pizarra::nueva()
        .variable("n", n)
        .variable("tabla", contados(peptidos, "péptido", "péptidos"))
        .variable("peptido", peptido.unwrap_or("–"))
        .variable("puntajes", format!("{} puntajes", miles(puntajes)))
        .variable("corte", opcional(corte))
        .variable("conservados", conservados)
}

fn narrar_recorte(
    peptidos: &[(String, Vec<usize>)],
    espectro: &[usize],
    n: usize,
) -> Result<(Vec<Paso>, Vec<Vec<usize>>), String> {
    let mut experimental = antibioticos::Espectro::nuevo(espectro)?;
    let cuantos = peptidos.len();
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let mut filas: Vec<Vec<String>> = peptidos
        .iter()
        .map(|(nombre, masas)| {
            vec![
                recortado(nombre, 16),
                recortado(&texto_de_masas(masas), 22),
                "?".to_string(),
                "·".to_string(),
            ]
        })
        .collect();
    let columnas = ["Péptido", "Masas", "Puntaje lineal", "¿Se queda?"];
    let (visibles, _) = ventana(&filas, 0, 10, vec!["…".to_string(); 4]);
    narracion.agregar(
        pizarra_de_recorte(n, cuantos, None, 0, None, 0)
            .tabla("La tabla de péptidos", &columnas, visibles, &[])
            .paso(
                "let mut experimental = Espectro::nuevo(espectro)?;",
                format!(
                    "La tabla tiene {} y solo caben {} (más los que empaten con el último). «Mejor» quiere decir más masas en común entre el espectro lineal del péptido y el espectro medido, que tiene {}. Primero ordeno el espectro medido y anoto cuántas veces sale cada masa.",
                    contados(cuantos, "péptido", "péptidos"),
                    miles(n),
                    contados(experimental.masas.len(), "masa", "masas")
                ),
            ),
    );
    if cuantos <= n {
        narracion.agregar(pizarra_de_recorte(n, cuantos, None, 0, None, 0).paso(
            "if tabla.len() <= n {",
            format!(
                "La tabla tiene {} y caben {}: no hace falta descartar a nadie.",
                miles(cuantos),
                miles(n)
            ),
        ));
        let ultimo = pizarra_de_recorte(n, cuantos, None, 0, None, cuantos).paso(
            "return tabla;",
            "Se devuelve la tabla entera, igual que como llegó. Prueba con un N más chico para ver el recorte.",
        );
        let tabla: Vec<Vec<usize>> = peptidos.iter().map(|(_, masas)| masas.clone()).collect();
        return Ok((narracion.terminar(ultimo), tabla));
    }
    narracion.agregar(pizarra_de_recorte(n, cuantos, None, 0, None, 0).paso(
        "if tabla.len() <= n {",
        format!(
            "{} no caben en {}: hay que recortar. Empiezo puntuando a cada péptido.",
            contados(cuantos, "péptido", "péptidos"),
            miles(n)
        ),
    ));
    let cada = cada_cuanto(cuantos, 8);
    let mut puntajes: Vec<usize> = Vec::with_capacity(cuantos);
    for (indice, (nombre, masas)) in peptidos.iter().enumerate() {
        let puntaje = experimental.compartidas(masas, false);
        puntajes.push(puntaje);
        filas[indice][2] = puntaje.to_string();
        if !toca(indice, 5, cada) {
            continue;
        }
        let (visibles, relativo) = ventana(&filas, indice, 10, vec!["…".to_string(); 4]);
        let teoricas = masas.len() * (masas.len() + 1) / 2 + 1;
        narracion.agregar(
            pizarra_de_recorte(n, cuantos, Some(&recortado(nombre, 16)), puntajes.len(), None, 0)
                .tabla("La tabla de péptidos", &columnas, visibles, &[(relativo, 2)])
                .paso(
                    "puntajes.push(experimental.compartidas(peptido, false));",
                    format!(
                        "{}: su espectro lineal tiene {} masas y {} están también en el espectro medido. Puntaje: {puntaje}.",
                        recortado(nombre, 16),
                        miles(teoricas),
                        miles(puntaje)
                    ),
                ),
        );
    }
    let mut de_mayor_a_menor = puntajes.clone();
    de_mayor_a_menor.sort_unstable_by(|a, b| b.cmp(a));
    let textos: Vec<String> = de_mayor_a_menor.iter().map(|p| p.to_string()).collect();
    let (visibles, relativo) = ventana(&textos, n - 1, 14, "…".to_string());
    narracion.agregar(
        pizarra_de_recorte(n, cuantos, None, puntajes.len(), None, 0)
            .lista("Puntajes de mayor a menor", visibles, &[relativo])
            .paso(
                "de_mayor_a_menor.sort_unstable_by(|a, b| b.cmp(a));",
                format!(
                    "Ordeno los puntajes de mayor a menor. El puesto {} (el último que cabe) está resaltado.",
                    miles(n)
                ),
            ),
    );
    let corte = de_mayor_a_menor[n - 1];
    narracion.agregar(
        pizarra_de_recorte(n, cuantos, None, puntajes.len(), Some(corte), 0).paso(
            "let corte = de_mayor_a_menor[n - 1];",
            format!(
                "El puntaje del puesto {} es {corte}: ese es el corte. Se queda quien tenga {corte} o más, y por eso también los que empatan con el último.",
                miles(n)
            ),
        ),
    );
    let mut conservados: Vec<Vec<usize>> = Vec::new();
    for (indice, ((nombre, masas), &puntaje)) in peptidos.iter().zip(&puntajes).enumerate() {
        let se_queda = puntaje >= corte;
        filas[indice][3] = if se_queda { "sí" } else { "no" }.to_string();
        if se_queda {
            conservados.push(masas.clone());
        }
        if !toca(indice, 5, cada) {
            continue;
        }
        let (visibles, relativo) = ventana(&filas, indice, 10, vec!["…".to_string(); 4]);
        narracion.agregar(
            pizarra_de_recorte(
                n,
                cuantos,
                Some(&recortado(nombre, 16)),
                puntajes.len(),
                Some(corte),
                conservados.len(),
            )
            .tabla(
                "La tabla de péptidos",
                &columnas,
                visibles,
                &[(relativo, 3)],
            )
            .paso(
                if se_queda {
                    "conservados.push(peptido);"
                } else {
                    "if puntaje >= corte {"
                },
                if se_queda {
                    format!(
                        "{}: puntaje {puntaje} es mayor o igual que el corte ({corte}). Se queda.",
                        recortado(nombre, 16)
                    )
                } else {
                    format!(
                        "{}: puntaje {puntaje} es menor que el corte ({corte}). Se descarta.",
                        recortado(nombre, 16)
                    )
                },
            ),
        );
    }
    let nombres_de_los_que_quedan: Vec<String> = peptidos
        .iter()
        .zip(&puntajes)
        .filter(|(_, puntaje)| **puntaje >= corte)
        .map(|((nombre, _), _)| recortado(nombre, 16))
        .collect();
    let (visibles, relativo) = ventana(&filas, 0, 10, vec!["…".to_string(); 4]);
    let _ = relativo;
    let ultimo = pizarra_de_recorte(
        n,
        cuantos,
        None,
        puntajes.len(),
        Some(corte),
        conservados.len(),
    )
    .tabla("La tabla de péptidos", &columnas, visibles, &[])
    .lista(
        "Se quedan",
        nombres_de_los_que_quedan.into_iter().take(20).collect(),
        &[],
    )
    .paso(
        "Ok(recortar_tabla(tabla, &mut experimental, n))",
        format!(
            "Quedan {} de {} (pedías {}). Mira «Qué significa» para ver qué pasó con los empates.",
            miles(conservados.len()),
            miles(cuantos),
            miles(n)
        ),
    );
    Ok((narracion.terminar(ultimo), conservados))
}

struct EstadoLideres {
    tabla: usize,
    revisados: usize,
    masa: Option<usize>,
    puntaje: Option<usize>,
}

fn pizarra_de_lideres(
    n: usize,
    alfabeto: usize,
    madre: usize,
    estado: &EstadoLideres,
    lideres: &Lideres,
) -> Pizarra {
    Pizarra::nueva()
        .variable("n", n)
        .variable("alfabeto", format!("{alfabeto} masas"))
        .variable("masa_madre", miles(madre))
        .variable("tabla", contados(estado.tabla, "péptido", "péptidos"))
        .variable("revisados", miles(estado.revisados))
        .variable("masa", opcional(estado.masa))
        .variable("puntaje", opcional(estado.puntaje))
        .variable("lideres.puntaje", lideres.puntaje)
        .variable("lideres.peptidos", lideres.peptidos.len())
}

fn candidatos_en_lista(candidatos: &[Vec<usize>], cuantos: usize) -> Vec<String> {
    candidatos
        .iter()
        .take(cuantos)
        .map(|peptido| recortado(&texto_de_masas(peptido), 30))
        .collect()
}

fn narrar_lideres(
    narracion: &mut Narracion,
    espectro: &[usize],
    n: usize,
    alfabeto: &[usize],
    linea_inicial: &str,
    presentacion: &str,
) -> Result<(Lideres, Paso), String> {
    if alfabeto.is_empty() || alfabeto.contains(&0) {
        return Err("Hace falta al menos una masa de aminoácido mayor que cero.".to_string());
    }
    let mut experimental = antibioticos::Espectro::nuevo(espectro)?;
    let madre = experimental.masa_madre;
    let mut lideres = Lideres {
        alfabeto: alfabeto.to_vec(),
        puntaje: 0,
        peptidos: Vec::new(),
    };
    let mut estado = EstadoLideres {
        tabla: 1,
        revisados: 0,
        masa: None,
        puntaje: None,
    };
    narracion.agregar(
        pizarra_de_lideres(n, alfabeto.len(), madre, &estado, &lideres)
            .lista(
                "Alfabeto: masas de aminoácido que se pueden usar",
                alfabeto.iter().map(|masa| masa_con_letras(*masa)).collect(),
                &[],
            )
            .paso(
                linea_inicial,
                format!(
                    "{presentacion} El espectro tiene {} y la mayor, {} Da, es la masa del péptido entero (la masa madre). La idea: hacer crecer péptidos un aminoácido por ronda y, en cada ronda, quedarme solo con los {} mejores.",
                    contados(experimental.masas.len(), "masa", "masas"),
                    miles(madre),
                    miles(n)
                ),
            ),
    );
    let mut tabla: Vec<Vec<usize>> = vec![Vec::new()];
    narracion.agregar(
        pizarra_de_lideres(n, alfabeto.len(), madre, &estado, &lideres).paso(
            "let mut tabla: Vec<Vec<usize>> = vec![Vec::new()];",
            "La tabla empieza con un solo péptido: el vacío, sin ningún aminoácido.",
        ),
    );
    let mut ronda: usize = 0;
    let mut evaluaciones = 0;
    let mut nuevos_lideres = 0;
    let mut empates = 0;
    while !tabla.is_empty() {
        ronda += 1;
        let anteriores = tabla.len();
        tabla = antibioticos::expandir(&tabla, alfabeto)?;
        estado.revisados += tabla.len();
        estado.tabla = tabla.len();
        if estado.revisados > antibioticos::MAXIMO_CANDIDATOS_EN_TOTAL {
            return Err(
                "La tabla de líderes revisó más de dos millones de candidatos en total, el límite: hay muchos empates y la tabla no deja de crecer."
                    .to_string(),
            );
        }
        if ronda <= 3 || ronda.is_multiple_of(5) {
            narracion.agregar(
                pizarra_de_lideres(n, alfabeto.len(), madre, &estado, &lideres)
                    .lista(
                        "Candidatos nuevos (los primeros)",
                        candidatos_en_lista(&tabla, 10),
                        &[],
                    )
                    .paso(
                        "tabla = expandir(&tabla, alfabeto)?;",
                        format!(
                            "Ronda {ronda}: {} de la tabla se alarga con cada una de las {} masas del alfabeto: salen {}.",
                            if anteriores == 1 {
                                "el único péptido".to_string()
                            } else {
                                format!("cada uno de los {} péptidos", miles(anteriores))
                            },
                            alfabeto.len(),
                            contados(tabla.len(), "candidato", "candidatos")
                        ),
                    ),
            );
        }
        let llegaron = tabla
            .iter()
            .filter(|peptido| antibioticos::masa_total(peptido) == madre)
            .count();
        let pasados = tabla
            .iter()
            .filter(|peptido| antibioticos::masa_total(peptido) > madre)
            .count();
        if ronda <= 2 || llegaron > 0 {
            narracion.agregar(
                pizarra_de_lideres(n, alfabeto.len(), madre, &estado, &lideres).paso(
                    "} else if masa < experimental.masa_madre {",
                    format!(
                        "Reviso el peso de cada candidato: {} pesan menos que la masa madre ({} Da) y siguen en carrera; {} pesan justo eso y ya son anillos completos, que puntúo contra el espectro entero; {} se pasaron de peso y se descartan.",
                        miles(tabla.len() - llegaron - pasados),
                        miles(madre),
                        miles(llegaron),
                        miles(pasados)
                    ),
                ),
            );
        }
        let mut en_carrera = Vec::with_capacity(tabla.len());
        for peptido in tabla {
            let masa = antibioticos::masa_total(&peptido);
            if masa == madre {
                let puntaje = experimental.compartidas(&peptido, true);
                estado.masa = Some(masa);
                estado.puntaje = Some(puntaje);
                if evaluaciones < 3 {
                    narracion.agregar(
                        pizarra_de_lideres(n, alfabeto.len(), madre, &estado, &lideres).paso(
                            "let puntaje = experimental.compartidas(&peptido, true);",
                            format!(
                                "El candidato {} pesa justo {} Da. Comparo su espectro cíclico con el medido: comparten {} masas.",
                                recortado(&texto_de_masas(&peptido), 40),
                                miles(madre),
                                miles(puntaje)
                            ),
                        ),
                    );
                }
                evaluaciones += 1;
                if lideres.peptidos.is_empty() || puntaje > lideres.puntaje {
                    let anterior = lideres.puntaje;
                    let habia = !lideres.peptidos.is_empty();
                    lideres.puntaje = puntaje;
                    lideres.peptidos.clear();
                    lideres.peptidos.push(peptido.to_vec());
                    nuevos_lideres += 1;
                    if nuevos_lideres <= 5 {
                        narracion.agregar(
                            pizarra_de_lideres(n, alfabeto.len(), madre, &estado, &lideres).paso(
                                "lideres.puntaje = puntaje;",
                                if habia {
                                    format!(
                                        "¡Nuevo líder! Su puntaje, {puntaje}, supera al mejor anterior ({anterior}): la lista de líderes se vacía y empieza con este péptido."
                                    )
                                } else {
                                    format!(
                                        "Es el primer anillo completo que aparece, así que es el líder por ahora, con puntaje {puntaje}."
                                    )
                                },
                            ),
                        );
                    }
                } else if puntaje == lideres.puntaje {
                    lideres.peptidos.push(peptido.to_vec());
                    empates += 1;
                    if empates <= 3 {
                        narracion.agregar(
                            pizarra_de_lideres(n, alfabeto.len(), madre, &estado, &lideres).paso(
                                "lideres.peptidos.push(peptido.to_vec());",
                                format!(
                                    "Empata con el líder (puntaje {puntaje}): se agrega a la lista. Ya hay {}.",
                                    contados(lideres.peptidos.len(), "líder", "líderes")
                                ),
                            ),
                        );
                    }
                }
                en_carrera.push(peptido);
            } else if masa < madre {
                en_carrera.push(peptido);
            }
        }
        estado.masa = None;
        estado.puntaje = None;
        let antes = en_carrera.len();
        let registrar_recorte = ronda <= 6 || ronda.is_multiple_of(3) || llegaron > 0;
        let mut mejores: Vec<Vec<String>> = Vec::new();
        let mut corte = None;
        if registrar_recorte {
            let mut con_puntaje: Vec<(usize, usize)> = en_carrera
                .iter()
                .enumerate()
                .map(|(indice, peptido)| (indice, experimental.compartidas(peptido, false)))
                .collect();
            if antes > n {
                let mut ordenados: Vec<usize> =
                    con_puntaje.iter().map(|(_, puntaje)| *puntaje).collect();
                ordenados.sort_unstable_by(|a, b| b.cmp(a));
                corte = Some(ordenados[n - 1]);
            }
            con_puntaje.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
            for (indice, puntaje) in con_puntaje.iter().take(6) {
                mejores.push(vec![
                    recortado(&letras_de_peptido(&en_carrera[*indice]), 22),
                    antibioticos::masa_total(&en_carrera[*indice]).to_string(),
                    puntaje.to_string(),
                ]);
            }
        }
        tabla = antibioticos::recortar_tabla(en_carrera, &mut experimental, n);
        estado.tabla = tabla.len();
        if registrar_recorte {
            let nota = match corte {
                Some(corte) => format!(
                    "Ronda {ronda}: siguen en carrera {} y la tabla solo guarda {}. Puntúo a cada uno con su espectro lineal; el puesto {} tiene {}, ese es el corte, y los empates también se quedan. Quedan {}.",
                    miles(antes),
                    miles(n),
                    miles(n),
                    contados(corte, "punto", "puntos"),
                    miles(tabla.len())
                ),
                None => format!(
                    "Ronda {ronda}: siguen en carrera {} y caben todos en la tabla de {}, así que no se descarta a nadie.",
                    miles(antes),
                    miles(n)
                ),
            };
            narracion.agregar(
                pizarra_de_lideres(n, alfabeto.len(), madre, &estado, &lideres)
                    .tabla(
                        "Los mejores de la ronda (puntaje lineal)",
                        &["Péptido", "Masa (Da)", "Puntaje"],
                        mejores,
                        &[],
                    )
                    .paso(
                        "tabla = recortar_tabla(en_carrera, &mut experimental, n);",
                        nota,
                    ),
            );
        }
    }
    let anillos = anillos_distintos(&lideres.peptidos);
    let resumen: Vec<String> = lideres
        .peptidos
        .iter()
        .take(10)
        .map(|peptido| recortado(&texto_de_masas(peptido), 30))
        .collect();
    let nota = if lideres.peptidos.is_empty() {
        format!(
            "La tabla se vació sin que ningún candidato pesara justo {} Da: no hay péptido. En «Qué significa» te digo qué probar.",
            miles(madre)
        )
    } else {
        format!(
            "No quedan candidatos: la búsqueda terminó en {} rondas. El mejor puntaje cíclico es {} y lo alcanzan {} ({}, porque un mismo anillo se puede leer desde varios puntos). Mira «Qué significa».",
            ronda,
            miles(lideres.puntaje),
            contados(lideres.peptidos.len(), "péptido", "péptidos"),
            contados(anillos, "anillo distinto", "anillos distintos")
        )
    };
    estado.tabla = 0;
    let ultimo = pizarra_de_lideres(n, alfabeto.len(), madre, &estado, &lideres)
        .lista("Péptidos líderes", resumen, &[])
        .paso("Ok(lideres)", nota);
    Ok((lideres, ultimo))
}

fn pasos_lideres(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let espectro = leer_espectro(entrada, "espectro")?;
    revisar_para_lideres(&espectro)?;
    let n = leer_n(entrada, "n", MAXIMO_N)?;
    antibioticos::secuenciar_con_lideres(&espectro, n)?;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let (_, ultimo) = narrar_lideres(
        &mut narracion,
        &espectro,
        n,
        &MASAS_DISTINTAS,
        "secuenciar_con_alfabeto(espectro, n, &MASAS_DISTINTAS)",
        "Voy a secuenciar un péptido cíclico desde un espectro con ruido, usando como alfabeto las 18 masas de aminoácido estándar.",
    )?;
    Ok(narracion.terminar(ultimo))
}

fn pasos_convolucion(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let espectro = leer_espectro(entrada, "espectro")?;
    if espectro.len() > MAXIMO_MASAS_LIDERES {
        return Err(format!(
            "El espectro tiene {} masas; aquí el límite es {}.",
            miles(espectro.len()),
            miles(MAXIMO_MASAS_LIDERES)
        ));
    }
    antibioticos::convolucion_espectral(&espectro)?;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let (_, ultimo) = narrar_convolucion(&mut narracion, &espectro)?;
    Ok(narracion.terminar(ultimo))
}

fn tabla_de_convolucion(convolucion: &[(usize, usize)], cuantas: usize) -> Vec<Vec<String>> {
    convolucion
        .iter()
        .take(cuantas)
        .map(|(diferencia, veces)| {
            vec![
                diferencia.to_string(),
                miles(*veces),
                nombre_de_diferencia(*diferencia),
            ]
        })
        .collect()
}

fn narrar_convolucion(
    narracion: &mut Narracion,
    espectro: &[usize],
) -> Result<(Vec<(usize, usize)>, Paso), String> {
    let experimental = antibioticos::Espectro::nuevo(espectro)?;
    let masas = &experimental.masas;
    let madre = experimental.masa_madre;
    let (visibles, _) = ventana(masas, 0, 24, 0);
    narracion.agregar(
        Pizarra::nueva()
            .variable("masa_madre", miles(madre))
            .variable("mayor", "–")
            .variable("menor", "–")
            .variable("diferencia", "–")
            .variable("repeticiones[diferencia]", "–")
            .lista(
                "Espectro ordenado (las primeras masas)",
                visibles.iter().map(|masa| masa.to_string()).collect(),
                &[],
            )
            .paso(
                "let experimental = Espectro::nuevo(espectro)?;",
                format!(
                    "El espectro tiene {}, de {} a {} Da, y ya está ordenado. La convolución resta TODAS las masas de a pares: si dos masas difieren justo en lo que pesa un aminoácido, es porque un pedazo del péptido tiene ese aminoácido de más.",
                    contados(masas.len(), "masa", "masas"),
                    miles(masas[0]),
                    miles(madre)
                ),
            ),
    );
    let mut repeticiones = vec![0usize; madre + 1];
    narracion.agregar(
        Pizarra::nueva()
            .variable("masa_madre", miles(madre))
            .variable("mayor", "–")
            .variable("menor", "–")
            .variable("diferencia", "–")
            .variable("repeticiones[diferencia]", "–")
            .paso(
                "let mut repeticiones = vec![0usize; experimental.masa_madre + 1];",
                format!(
                    "Preparo una casilla por cada diferencia posible (de 0 a {}), todas en 0. Ahí iré contando cuántas veces sale cada diferencia.",
                    miles(madre)
                ),
            ),
    );
    let cada = cada_cuanto(masas.len(), 18);
    let mut pares = 0usize;
    let mut pasos_de_pares = 0;
    for (posicion, &mayor) in masas.iter().enumerate() {
        for &menor in &masas[..posicion] {
            let diferencia = mayor - menor;
            if diferencia > 0 {
                repeticiones[diferencia] += 1;
            }
            pares += 1;
            if pasos_de_pares < 6 && posicion <= 3 {
                pasos_de_pares += 1;
                let (nota, linea) = if diferencia > 0 {
                    (
                        format!(
                            "{mayor} − {menor} = {diferencia}. Sumo 1 a la casilla {diferencia}, que ya lleva {}.",
                            repeticiones[diferencia]
                        ),
                        "repeticiones[diferencia] += 1;",
                    )
                } else {
                    (
                        format!(
                            "{mayor} − {menor} = 0: son masas iguales, y una diferencia de 0 no dice nada, así que no se cuenta."
                        ),
                        "if diferencia > 0 {",
                    )
                };
                let (elementos, relativo) = ventana(
                    &masas
                        .iter()
                        .map(|masa| masa.to_string())
                        .collect::<Vec<_>>(),
                    posicion,
                    20,
                    "…".to_string(),
                );
                narracion.agregar(
                    Pizarra::nueva()
                        .variable("masa_madre", miles(madre))
                        .variable("mayor", mayor)
                        .variable("menor", menor)
                        .variable("diferencia", diferencia)
                        .variable("repeticiones[diferencia]", repeticiones[diferencia])
                        .lista("Espectro ordenado", elementos, &[relativo])
                        .paso(linea, nota),
                );
            }
        }
        if posicion > 3 && (posicion.is_multiple_of(cada) || posicion == masas.len() - 1) {
            let mut mejores: Vec<(usize, usize)> = repeticiones
                .iter()
                .enumerate()
                .filter(|(_, veces)| **veces > 0)
                .map(|(diferencia, veces)| (diferencia, *veces))
                .collect();
            mejores.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
            let (elementos, relativo) = ventana(
                &masas
                    .iter()
                    .map(|masa| masa.to_string())
                    .collect::<Vec<_>>(),
                posicion,
                20,
                "…".to_string(),
            );
            narracion.agregar(
                Pizarra::nueva()
                    .variable("masa_madre", miles(madre))
                    .variable("mayor", mayor)
                    .variable("menor", masas[posicion - 1])
                    .variable("diferencia", mayor - masas[posicion - 1])
                    .variable(
                        "repeticiones[diferencia]",
                        repeticiones[mayor - masas[posicion - 1]],
                    )
                    .lista("Espectro ordenado", elementos, &[relativo])
                    .tabla(
                        "Diferencias más repetidas hasta ahora",
                        &["Diferencia", "Veces", "Aminoácido"],
                        tabla_de_convolucion(&mejores, 5),
                        &[],
                    )
                    .paso(
                        "for &menor in &masas[..posicion] {",
                        format!(
                            "Con mayor = {mayor} resto cada una de las {} masas anteriores. Van {} restas hechas. Las demás masas se procesan igual; lo que importa es cuáles diferencias se repiten más.",
                            miles(posicion),
                            miles(pares)
                        ),
                    ),
            );
        }
    }
    let mut convolucion = Vec::new();
    for (diferencia, &cuantas) in repeticiones.iter().enumerate() {
        if cuantas > 0 {
            convolucion.push((diferencia, cuantas));
        }
    }
    let total: usize = convolucion.iter().map(|&(_, veces)| veces).sum();
    narracion.agregar(
        Pizarra::nueva()
            .variable("masa_madre", miles(madre))
            .variable("convolucion", format!("{} diferencias", miles(convolucion.len())))
            .paso(
                "convolucion.push((diferencia, cuantas));",
                format!(
                    "Terminé de restar: {} restas positivas. Paso la tabla a una lista con solo las diferencias que salieron al menos una vez: {}.",
                    miles(total),
                    contados(convolucion.len(), "diferencia distinta", "diferencias distintas")
                ),
            ),
    );
    convolucion.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    let filas = tabla_de_convolucion(&convolucion, 10);
    narracion.agregar(
        Pizarra::nueva()
            .variable("masa_madre", miles(madre))
            .variable("convolucion", format!("{} diferencias", miles(convolucion.len())))
            .tabla(
                "Las diferencias más repetidas",
                &["Diferencia", "Veces", "Aminoácido"],
                filas.clone(),
                &[],
            )
            .paso(
                "convolucion.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));",
                "Ordeno de la que más veces salió a la que menos (a igualdad, la menor primero). Si una masa de aminoácido sale muy repetida, es buena señal de que el péptido la contiene.",
            ),
    );
    let nota = format!(
        "Resultado: {} diferencias distintas; la más repetida es {} ({} veces). En «Qué significa» te cuento cómo usarlas para adivinar los aminoácidos.",
        miles(convolucion.len()),
        convolucion.first().map_or(0, |(diferencia, _)| *diferencia),
        miles(convolucion.first().map_or(0, |(_, veces)| *veces))
    );
    let ultimo = Pizarra::nueva()
        .variable("masa_madre", miles(madre))
        .variable(
            "convolucion",
            format!("{} diferencias", miles(convolucion.len())),
        )
        .tabla(
            "Las diferencias más repetidas",
            &["Diferencia", "Veces", "Aminoácido"],
            filas,
            &[],
        )
        .paso("Ok(convolucion)", nota);
    Ok((convolucion, ultimo))
}

fn pasos_convolucion_lideres(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let espectro = leer_espectro(entrada, "espectro")?;
    revisar_para_lideres(&espectro)?;
    let m = leer_n(entrada, "m", MAXIMO_M)?;
    let n = leer_n(entrada, "n", MAXIMO_N)?;
    antibioticos::secuenciar_por_convolucion(&espectro, m, n)?;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let (_, ultimo) = narrar_por_convolucion(&mut narracion, &espectro, m, n)?;
    Ok(narracion.terminar(ultimo))
}

fn pizarra_de_alfabeto(
    m: usize,
    convolucion: usize,
    en_rango: Option<usize>,
    corte: Option<usize>,
    alfabeto: Option<usize>,
) -> Pizarra {
    Pizarra::nueva()
        .variable("m", m)
        .variable("convolucion", format!("{} diferencias", miles(convolucion)))
        .variable("en_rango", opcional(en_rango))
        .variable("corte", opcional(corte))
        .variable("alfabeto", opcional(alfabeto))
}

fn narrar_por_convolucion(
    narracion: &mut Narracion,
    espectro: &[usize],
    m: usize,
    n: usize,
) -> Result<(Lideres, Paso), String> {
    narracion.agregar(
        pizarra_de_alfabeto(m, 0, None, None, None).paso(
            "let alfabeto = alfabeto_por_convolucion(espectro, m)?;",
            format!(
                "Aquí no se sabe de qué aminoácidos está hecho el péptido: pueden ser raros, no estándar. Así que primero le pregunto al espectro, con la convolución, qué masas parecen piezas. Pido las M = {} más repetidas y después buscaré con una tabla de {} péptidos (N).",
                miles(m),
                miles(n)
            ),
        ),
    );
    let (convolucion, ultimo_de_la_convolucion) = narrar_convolucion(narracion, espectro)?;
    narracion.agregar(ultimo_de_la_convolucion);
    let mut en_rango: Vec<(usize, usize)> = Vec::new();
    for &(diferencia, cuantas) in &convolucion {
        if (antibioticos::MASA_MINIMA_DE_AMINOACIDO..=antibioticos::MASA_MAXIMA_DE_AMINOACIDO)
            .contains(&diferencia)
        {
            en_rango.push((diferencia, cuantas));
        }
    }
    if en_rango.is_empty() {
        return Err(format!(
            "La convolución no tiene ninguna diferencia entre {} y {} Da, así que no sugiere ningún aminoácido.",
            antibioticos::MASA_MINIMA_DE_AMINOACIDO,
            antibioticos::MASA_MAXIMA_DE_AMINOACIDO
        ));
    }
    narracion.agregar(
        pizarra_de_alfabeto(m, convolucion.len(), Some(en_rango.len()), None, None)
            .tabla(
                "Diferencias entre 57 y 200 Da (las más repetidas)",
                &["Diferencia", "Veces", "Aminoácido"],
                tabla_de_convolucion(&en_rango, 12),
                &[],
            )
            .paso(
                "en_rango.push((diferencia, cuantas));",
                format!(
                    "De todas las diferencias solo me sirven las que caben en un aminoácido: entre {} y {} Da. Quedan {}, ya ordenadas de más a menos repetida.",
                    antibioticos::MASA_MINIMA_DE_AMINOACIDO,
                    antibioticos::MASA_MAXIMA_DE_AMINOACIDO,
                    contados(en_rango.len(), "candidata", "candidatas")
                ),
            ),
    );
    let puesto = m.min(en_rango.len());
    let corte = en_rango[puesto - 1].1;
    narracion.agregar(
        pizarra_de_alfabeto(m, convolucion.len(), Some(en_rango.len()), Some(corte), None)
            .tabla(
                "Diferencias entre 57 y 200 Da (las más repetidas)",
                &["Diferencia", "Veces", "Aminoácido"],
                tabla_de_convolucion(&en_rango, 12),
                &[(puesto - 1, 1)],
            )
            .paso(
                "let corte = en_rango[m.min(en_rango.len()) - 1].1;",
                format!(
                    "La diferencia del puesto {} es {} y sale {} veces: ese número es el corte. Entran todas las que salgan al menos esas veces, así que los empates con ella también entran.",
                    miles(puesto),
                    en_rango[puesto - 1].0,
                    miles(corte)
                ),
            ),
    );
    let mut alfabeto = Vec::new();
    for &(diferencia, cuantas) in &en_rango {
        if cuantas >= corte {
            alfabeto.push(diferencia);
        }
    }
    alfabeto.sort_unstable();
    narracion.agregar(
        pizarra_de_alfabeto(
            m,
            convolucion.len(),
            Some(en_rango.len()),
            Some(corte),
            Some(alfabeto.len()),
        )
        .lista(
            "Alfabeto elegido por la convolución",
            alfabeto.iter().map(|masa| masa_con_letras(*masa)).collect(),
            &[],
        )
        .paso(
            "alfabeto.sort_unstable();",
            format!(
                "Me quedo con {} y las ordeno de menor a mayor: ese es el alfabeto de este péptido.{}",
                contados(alfabeto.len(), "masa", "masas"),
                if alfabeto.len() > m {
                    format!(
                        " Son más que las {m} que pediste porque hubo empates con la última."
                    )
                } else {
                    String::new()
                }
            ),
        ),
    );
    narrar_lideres(
        narracion,
        espectro,
        n,
        &alfabeto,
        "secuenciar_con_alfabeto(espectro, n, &alfabeto)",
        "Con ese alfabeto, ahora hago la búsqueda con tabla de líderes.",
    )
}

struct Retroceso<'a> {
    narracion: &'a mut Narracion,
    pendientes: BTreeMap<i64, usize>,
    puntos: Vec<i64>,
    ancho: i64,
    intentos: usize,
    eventos: usize,
    cantidad_de_puntos: usize,
}

fn quitar_pendiente(pendientes: &mut BTreeMap<i64, usize>, distancia: i64) -> bool {
    match pendientes.get_mut(&distancia) {
        Some(cuantas) => {
            *cuantas -= 1;
            if *cuantas == 0 {
                pendientes.remove(&distancia);
            }
            true
        }
        None => false,
    }
}

impl Retroceso<'_> {
    fn pizarra(&self, mayor: Option<i64>, candidato: Option<i64>) -> Pizarra {
        let lista: Vec<String> = self
            .pendientes
            .iter()
            .rev()
            .take(10)
            .map(|(distancia, veces)| format!("{distancia} (×{veces})"))
            .collect();
        let resaltar: Vec<usize> = if mayor.is_some() && !lista.is_empty() {
            vec![0]
        } else {
            Vec::new()
        };
        Pizarra::nueva()
            .variable("cantidad_de_puntos", self.cantidad_de_puntos)
            .variable("ancho", self.ancho)
            .variable("puntos", recortado(&format!("{:?}", self.puntos), 60))
            .variable("intentos", self.intentos)
            .variable(
                "mayor",
                mayor.map_or_else(|| "–".to_string(), |v| v.to_string()),
            )
            .variable(
                "candidato",
                candidato.map_or_else(|| "–".to_string(), |v| v.to_string()),
            )
            .lista(
                "Distancias que faltan por explicar (de mayor a menor)",
                lista,
                &resaltar,
            )
    }

    fn registrar(&mut self) -> bool {
        self.eventos += 1;
        let evento = self.eventos;
        evento <= 24 || evento.is_multiple_of((evento.next_power_of_two() / 8).max(1))
    }

    fn probar(&mut self, candidato: i64, grabar: bool) -> (bool, Vec<Vec<String>>, Option<usize>) {
        let mut quitadas = Vec::with_capacity(self.puntos.len());
        let mut filas = Vec::new();
        for (indice, &punto) in self.puntos.iter().enumerate() {
            let distancia = (candidato - punto).abs();
            let hay = quitar_pendiente(&mut self.pendientes, distancia);
            if grabar && filas.len() < 12 {
                filas.push(vec![
                    punto.to_string(),
                    distancia.to_string(),
                    if hay {
                        "sí, la tacho"
                    } else {
                        "ya no queda ninguna"
                    }
                    .to_string(),
                ]);
            }
            if hay {
                quitadas.push(distancia);
            } else {
                for devuelta in quitadas {
                    *self.pendientes.entry(devuelta).or_insert(0) += 1;
                }
                return (false, filas, Some(indice));
            }
        }
        (true, filas, None)
    }

    fn colocar(&mut self) -> Result<bool, String> {
        let Some(&mayor) = self.pendientes.keys().next_back() else {
            let ancho = self.ancho;
            let paso = self.pizarra(None, None).paso(
                "return Ok(true);",
                format!(
                    "No queda ninguna distancia por explicar: todas las distancias de tu lista salen de los puntos colocados. ¡Hay solución! ({} puntos entre 0 y {ancho}).",
                    self.puntos.len()
                ),
            );
            self.narracion.agregar(paso);
            return Ok(true);
        };
        self.intentos += 1;
        if self.intentos > antibioticos::MAXIMO_INTENTOS {
            return Err(
                "La búsqueda con retroceso probó un millón de posiciones (el límite) sin terminar: estas distancias obligan a retroceder demasiadas veces. Revisa si sobra o falta alguna."
                    .to_string(),
            );
        }
        let mut opciones = vec![mayor];
        if self.ancho - mayor != mayor {
            opciones.push(self.ancho - mayor);
        }
        if self.registrar() {
            let nota = if opciones.len() == 2 {
                format!(
                    "La mayor distancia que falta es {mayor}. Tiene que ser la distancia de un punto nuevo a uno de los extremos: si es al 0, el punto está en {mayor}; si es al {}, está en {}. Pruebo las dos opciones.",
                    self.ancho,
                    self.ancho - mayor
                )
            } else {
                format!(
                    "La mayor distancia que falta es {mayor}, justo la mitad del ancho: el punto nuevo solo puede estar en {mayor}. Hay una única opción."
                )
            };
            let paso = self
                .pizarra(Some(mayor), None)
                .paso("let mut opciones = vec![mayor];", nota);
            self.narracion.agregar(paso);
        }
        for candidato in opciones {
            let grabar = self.registrar();
            let (cabe, filas, fallo) = self.probar(candidato, grabar);
            if cabe {
                self.puntos.push(candidato);
                if grabar {
                    let resaltar: Vec<(usize, usize)> =
                        (0..filas.len()).map(|fila| (fila, 2)).collect();
                    let paso = self
                        .pizarra(Some(mayor), Some(candidato))
                        .tabla(
                            "Distancias del candidato a los puntos ya puestos",
                            &["Punto", "Distancia", "¿Estaba pendiente?"],
                            filas,
                            &resaltar,
                        )
                        .paso(
                            "puntos.push(candidato);",
                            format!(
                                "El candidato {candidato}: todas sus distancias a los puntos ya puestos estaban pendientes. Las tacho de la tabla y lo coloco. Sigo con la mayor distancia que falte."
                            ),
                        );
                    self.narracion.agregar(paso);
                }
                if self.colocar()? {
                    return Ok(true);
                }
                self.puntos.pop();
                for &punto in &self.puntos {
                    *self
                        .pendientes
                        .entry((candidato - punto).abs())
                        .or_insert(0) += 1;
                }
                if self.registrar() {
                    let paso = self.pizarra(Some(mayor), Some(candidato)).paso(
                        "puntos.pop();",
                        format!(
                            "Por este camino no hubo solución: quito el punto {candidato} y devuelvo sus distancias a la tabla. Esto es el retroceso: deshacer lo último y probar otra cosa."
                        ),
                    );
                    self.narracion.agregar(paso);
                }
            } else if grabar {
                let punto = fallo.map_or(0, |indice| self.puntos[indice]);
                let resaltar: Vec<(usize, usize)> =
                    fallo.map(|fila| (fila, 2)).into_iter().collect();
                let paso = self
                    .pizarra(Some(mayor), Some(candidato))
                    .tabla(
                        "Distancias del candidato a los puntos ya puestos",
                        &["Punto", "Distancia", "¿Estaba pendiente?"],
                        filas,
                        &resaltar,
                    )
                    .paso(
                        "if quitar_distancias_a(candidato, puntos, pendientes) {",
                        format!(
                            "El candidato {candidato} quedaría a distancia {} del punto {punto}, pero esa distancia ya no está en la lista. No puede ser: lo descarto y devuelvo lo que había tachado.",
                            (candidato - punto).abs()
                        ),
                    );
                self.narracion.agregar(paso);
            }
        }
        if self.registrar() {
            let paso = self.pizarra(Some(mayor), None).paso(
                "Ok(false)",
                "Ninguna de las opciones sirve con los puntos que hay: esta rama no tiene solución y se vuelve un paso atrás.",
            );
            self.narracion.agregar(paso);
        }
        Ok(false)
    }
}

fn pasos_autopista(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let distancias = entrada.enteros("distancias")?;
    if distancias.len() > MAXIMO_DISTANCIAS {
        return Err(format!(
            "Hay {} distancias y el límite es {}.",
            miles(distancias.len()),
            miles(MAXIMO_DISTANCIAS)
        ));
    }
    antibioticos::puntos_desde_distancias(&distancias)?;
    Ok(narrar_autopista(&distancias)?.0)
}

fn pendientes_en_lista(pendientes: &BTreeMap<i64, usize>) -> Vec<String> {
    pendientes
        .iter()
        .rev()
        .take(10)
        .map(|(distancia, veces)| format!("{distancia} (×{veces})"))
        .collect()
}

fn narrar_autopista(distancias: &[i64]) -> Result<(Vec<Paso>, Option<Vec<i64>>), String> {
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let mut pendientes: BTreeMap<i64, usize> = BTreeMap::new();
    let mut cantidad = 0;
    for &distancia in distancias {
        if distancia > 0 {
            *pendientes.entry(distancia).or_insert(0) += 1;
            cantidad += 1;
        }
    }
    if cantidad == 0 {
        return Err("Hace falta al menos una distancia mayor que cero.".to_string());
    }
    let ignoradas = distancias.len() - cantidad;
    let mut nota = format!(
        "Tienes {} distancias positivas entre puntos de una recta. Las cuento en una tabla de pendientes: cada distancia con las veces que todavía no he explicado.",
        miles(cantidad)
    );
    if ignoradas > 0 {
        nota.push_str(&format!(
            " Ignoré {} valores menores o iguales a cero.",
            miles(ignoradas)
        ));
    }
    narracion.agregar(
        Pizarra::nueva()
            .variable("cantidad", cantidad)
            .lista(
                "Distancias pendientes (de mayor a menor)",
                pendientes_en_lista(&pendientes),
                &[],
            )
            .paso("*pendientes.entry(distancia).or_insert(0) += 1;", nota),
    );
    let Some(cantidad_de_puntos) = antibioticos::puntos_para(cantidad) else {
        return Err(format!(
            "Con n puntos hay n(n−1)/2 distancias (1, 3, 6, 10, 15…), pero aquí hay {cantidad} distancias positivas."
        ));
    };
    if cantidad_de_puntos > antibioticos::MAXIMO_PUNTOS {
        return Err(format!(
            "Esas distancias corresponden a {cantidad_de_puntos} puntos y el límite es {}: la búsqueda con retroceso podría tardar demasiado.",
            antibioticos::MAXIMO_PUNTOS
        ));
    }
    narracion.agregar(
        Pizarra::nueva()
            .variable("cantidad", cantidad)
            .variable("cantidad_de_puntos", cantidad_de_puntos)
            .paso(
                "let Some(cantidad_de_puntos) = puntos_para(cantidad) else {",
                format!(
                    "Con n puntos hay n·(n−1)/2 distancias (todos contra todos). {} = {}·{}/2, así que buscas {} puntos.",
                    miles(cantidad),
                    cantidad_de_puntos,
                    cantidad_de_puntos - 1,
                    cantidad_de_puntos
                ),
            ),
    );
    let ancho = *pendientes.keys().next_back().unwrap_or(&0);
    narracion.agregar(
        Pizarra::nueva()
            .variable("cantidad", cantidad)
            .variable("cantidad_de_puntos", cantidad_de_puntos)
            .variable("ancho", ancho)
            .lista(
                "Distancias pendientes (de mayor a menor)",
                pendientes_en_lista(&pendientes),
                &[0],
            )
            .paso(
                "let ancho = *pendientes.keys().next_back().unwrap_or(&0);",
                format!(
                    "La mayor distancia, {ancho}, solo puede ser la que hay entre el primer y el último punto: el ancho de toda la recta."
                ),
            ),
    );
    quitar_pendiente(&mut pendientes, ancho);
    let mut retroceso = Retroceso {
        narracion: &mut narracion,
        pendientes,
        puntos: vec![0, ancho],
        ancho,
        intentos: 0,
        eventos: 0,
        cantidad_de_puntos,
    };
    let paso = retroceso.pizarra(None, None).paso(
        "let mut puntos = vec![0, ancho];",
        format!(
            "Pongo los dos extremos: 0 y {ancho}, y tacho esa distancia de la tabla. Faltan {} puntos por colocar.",
            cantidad_de_puntos - 2
        ),
    );
    retroceso.narracion.agregar(paso);
    let paso = retroceso.pizarra(None, None).paso(
        "if colocar(&mut pendientes, &mut puntos, ancho, &mut intentos)? {",
        "Empieza la búsqueda con retroceso: coloco un punto a la vez, siempre explicando la mayor distancia que falte; si me equivoco, deshago y pruebo otra cosa.",
    );
    retroceso.narracion.agregar(paso);
    let encontrado = retroceso.colocar()?;
    let intentos = retroceso.intentos;
    let mut puntos = retroceso.puntos;
    let resultado = if encontrado {
        puntos.sort_unstable();
        Some(puntos)
    } else {
        None
    };
    let ultimo = match &resultado {
        Some(puntos) => Pizarra::nueva()
            .variable("cantidad_de_puntos", cantidad_de_puntos)
            .variable("ancho", ancho)
            .variable("puntos", recortado(&format!("{puntos:?}"), 80))
            .variable("intentos", intentos)
            .paso(
                "Ok(Some(puntos))",
                format!(
                    "Solución: {}. Tras {} el retroceso encontró un conjunto de puntos cuyas distancias son exactamente las tuyas. En «Qué significa» te cuento si es la única.",
                    recortado(&format!("{puntos:?}"), 120),
                    contados(intentos, "intento", "intentos")
                ),
            ),
        None => Pizarra::nueva()
            .variable("cantidad_de_puntos", cantidad_de_puntos)
            .variable("ancho", ancho)
            .variable("intentos", intentos)
            .paso(
                "Ok(None)",
                format!(
                    "Se probaron todas las posibilidades ({}) y ninguna explica todas las distancias: no hay solución. Mira «Qué significa» para ver qué suele fallar.",
                    contados(intentos, "intento", "intentos")
                ),
            ),
    };
    Ok((narracion.terminar(ultimo), resultado))
}

#[cfg(test)]
mod pruebas {
    use super::*;

    #[test]
    fn la_traza_encuentra_lo_mismo_que_el_nucleo() {
        let espectro: Vec<usize> = ESPECTRO_PALMERA
            .split_whitespace()
            .map(|masa| masa.parse().unwrap())
            .collect();
        let busqueda = buscar_por_ramificacion(&espectro).unwrap();
        let directo = antibioticos::secuenciar_ciclopeptido(&espectro).unwrap();
        assert_eq!(busqueda.encontrados, directo);
        assert_eq!(busqueda.encontrados.len(), 14);
        assert_eq!(anillos_distintos(&busqueda.encontrados), 1);
    }

    #[test]
    fn letras_ambiguas_y_no_estandar() {
        assert_eq!(letras_de_peptido(&[97, 113, 128, 69]), "P(I/L)(K/Q)[69]");
        assert_eq!(lecturas_del_anillo(&[103, 71, 131, 71]), 4);
        assert_eq!(lecturas_del_anillo(&[57, 71, 113]), 6);
    }

    #[test]
    fn formula_sin_subindices_de_uno() {
        let glicina = antibioticos::formula_del_peptido(b"G", false).unwrap();
        assert_eq!(texto_de_formula(&glicina), "C2H5NO2");
        let metionina = antibioticos::formula_del_peptido(b"M", false).unwrap();
        assert_eq!(texto_de_formula(&metionina), "C5H11NO2S");
    }

    #[test]
    fn el_estado_de_cada_poda_cuadra() {
        let entrada = Entrada::nueva(&[("espectro", ESPECTRO_PALMERA)]);
        let pasos = pasos_ciclopeptido(&entrada).unwrap();
        assert_eq!(pasos[0].estado["fase"], "inicio");
        let mut podas = 0;
        for paso in &pasos {
            let estado = &paso.estado;
            assert!(estado["candidatos"].as_array().unwrap().len() <= MAXIMO_MOSTRADOS);
            if estado["fase"] != "podar" {
                continue;
            }
            podas += 1;
            let numero = |campo: &str| estado[campo].as_u64().unwrap();
            assert_eq!(
                numero("inconsistentes")
                    + numero("llegaron_a_masa_madre")
                    + numero("total_candidatos"),
                numero("expandidos")
            );
            assert_eq!(
                numero("descartados") + numero("nuevos_encontrados") + numero("total_candidatos"),
                numero("expandidos")
            );
        }
        assert_eq!(podas, 7);
        let ultimo = pasos.last().unwrap();
        assert_eq!(ultimo.estado["fase"], "fin");
        assert_eq!(ultimo.estado["total_encontrados"], 14);
        assert_eq!(ultimo.estado["encontrados"].as_array().unwrap().len(), 14);
    }

    #[test]
    fn empate_entre_anillos_de_distinto_largo() {
        let resultado = Lideres {
            alfabeto: vec![57, 71, 128],
            puntaje: 3,
            peptidos: vec![vec![128, 87, 99], vec![57, 71, 87, 99]],
        };
        let mut salida = Salida::default();
        explicar_lideres(&mut salida, &resultado, &[0, 87, 99, 314]);
        let texto = salida.interpretacion.join(" ");
        assert!(
            texto.contains("128 en uno pesa lo mismo que 57 + 71"),
            "{texto}"
        );
        assert!(!texto.contains("Prueba con un N más grande"), "{texto}");
    }

    #[test]
    fn isotopos_segun_poisson() {
        let corto = masa(&Entrada::nueva(&[("peptido", "GAG"), ("forma", "lineal")])).unwrap();
        let texto = corto.interpretacion.join(" ");
        assert!(texto.contains("no trae ninguno"), "{texto}");
        assert!(texto.contains("el más alto de su grupo"), "{texto}");
        let largo = masa(&Entrada::nueva(&[
            ("peptido", "@proteina:dnaA"),
            ("forma", "lineal"),
        ]))
        .unwrap();
        let texto = largo.interpretacion.join(" ");
        assert!(texto.contains("al menos tan alto"), "{texto}");
    }

    #[test]
    fn agrega_el_cero_si_falta() {
        let entrada = Entrada::nueva(&[("espectro", "57 71 113 128 170 184 241")]);
        let salida = ciclopeptido(&entrada).unwrap();
        assert!(salida.interpretacion[0].contains("masa 0"));
        let pasos = pasos_ciclopeptido(&entrada).unwrap();
        assert!(pasos.len() >= 3);
    }

    #[test]
    fn la_lista_de_recorte_tiene_un_limite_de_trabajo() {
        let larga = "A".repeat(1000);
        let lista = vec![larga; 100].join(" ");
        let error = leer_lista_de_peptidos(&lista).unwrap_err();
        assert!(error.contains("demasiado largos en conjunto"), "{error}");
        assert!(leer_lista_de_peptidos("PALM PAIM MERA").is_ok());
    }

    #[test]
    fn el_solape_se_detecta_aunque_otra_hebra_este_en_medio() {
        let hallazgo = |inicio: usize, hebra: u8| antibioticos::Codificacion {
            inicio,
            hebra,
            fragmento: Vec::new(),
        };
        let revueltas = [hallazgo(0, b'+'), hallazgo(1, b'-'), hallazgo(2, b'+')];
        assert!(hay_solape(&revueltas, 3));
        let separadas = [hallazgo(0, b'+'), hallazgo(1, b'-'), hallazgo(3, b'+')];
        assert!(!hay_solape(&separadas, 3));
    }

    fn del_ejemplo(id: &str) -> Entrada {
        Entrada::desde_json(crate::ejemplo(id).unwrap()).unwrap()
    }

    fn masas_de(texto: &str) -> Vec<usize> {
        texto
            .split_whitespace()
            .map(|masa| masa.parse().unwrap())
            .collect()
    }

    #[test]
    fn narrar_la_traduccion_da_lo_mismo_que_el_nucleo() {
        let mut textos = vec![
            del_ejemplo("antibioticos.traducir")
                .secuencia("arn")
                .unwrap(),
        ];
        for arn in [
            "AUGGCUUAAGGG",
            "AUGGCUGCU",
            "UAA",
            "AUGGCUGCUGA",
            "ATGAAATAG",
        ] {
            textos.push(arn.as_bytes().to_vec());
        }
        for arn in textos {
            let esperada = antibioticos::traducir_arn(&arn).unwrap();
            let (pasos, obtenida) = narrar_traduccion(&arn).unwrap();
            assert_eq!(obtenida, esperada);
            assert!(!pasos.is_empty() && pasos.len() <= MAXIMO_PASOS);
        }
    }

    #[test]
    fn narrar_la_codificacion_da_lo_mismo_que_el_nucleo() {
        let entrada = del_ejemplo("antibioticos.codificar");
        let texto = entrada.adn("texto").unwrap();
        let peptido = entrada.proteina("peptido").unwrap();
        let mut casos = vec![(texto, peptido)];
        for (texto, peptido) in [
            ("CATCATCAT", "HHH"),
            ("ATGATGATGCATCATCAT", "MMM"),
            ("ATGATGATGCATCATCAT", "HHH"),
            ("ATGGCCATGGCGCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA", "MA"),
            ("ACGT", "MMMM"),
            ("AAA", "K"),
        ] {
            casos.push((texto.as_bytes().to_vec(), peptido.as_bytes().to_vec()));
        }
        for (texto, peptido) in casos {
            let esperadas = antibioticos::codificaciones_del_peptido(&texto, &peptido).unwrap();
            let (pasos, obtenidas) = narrar_codificacion(&texto, &peptido);
            assert_eq!(obtenidas, esperadas);
            assert!(!pasos.is_empty() && pasos.len() <= MAXIMO_PASOS);
        }
    }

    #[test]
    fn narrar_la_masa_da_lo_mismo_que_el_nucleo() {
        let ejemplo = del_ejemplo("antibioticos.masa")
            .proteina("peptido")
            .unwrap();
        for peptido in [ejemplo, b"PALMERA".to_vec(), b"G".to_vec(), b"MCW".to_vec()] {
            for ciclico in [false, true] {
                let formula = antibioticos::formula_del_peptido(&peptido, ciclico).unwrap();
                let masa = antibioticos::masa_monoisotopica(&peptido, ciclico).unwrap();
                let (_, obtenida, suma) = narrar_masa(&peptido, ciclico).unwrap();
                assert_eq!(obtenida, formula);
                assert_eq!(suma, masa);
            }
        }
    }

    #[test]
    fn narrar_el_espectro_da_lo_mismo_que_el_nucleo() {
        let lineal = del_ejemplo("antibioticos.espectro_lineal")
            .proteina("peptido")
            .unwrap();
        let ciclico = del_ejemplo("antibioticos.espectro_ciclico")
            .proteina("peptido")
            .unwrap();
        for peptido in [
            lineal,
            ciclico,
            b"G".to_vec(),
            b"GA".to_vec(),
            b"PALMERAPALMERA".to_vec(),
        ] {
            let masas = antibioticos::masas_del_peptido(&peptido).unwrap();
            let (pasos, obtenido) = narrar_espectro(&peptido, &masas, true);
            assert_eq!(obtenido, antibioticos::espectro_ciclico(&masas));
            assert!(!pasos.is_empty() && pasos.len() <= MAXIMO_PASOS);
            let (_, obtenido) = narrar_espectro(&peptido, &masas, false);
            assert_eq!(obtenido, antibioticos::espectro_lineal(&masas));
        }
    }

    #[test]
    fn narrar_el_conteo_da_lo_mismo_que_el_nucleo() {
        let ejemplo = del_ejemplo("antibioticos.contar").natural("masa").unwrap();
        for masa in [ejemplo, 0, 1, 56, 57, 58, 114, 500, 3000] {
            let (pasos, obtenido) = narrar_conteo(masa);
            assert_eq!(obtenido, antibioticos::peptidos_por_masa(masa).unwrap());
            assert!(!pasos.is_empty() && pasos.len() <= MAXIMO_PASOS);
        }
    }

    #[test]
    fn narrar_el_puntaje_da_lo_mismo_que_el_nucleo() {
        let espectro = masas_de(ESPECTRO_PALMERA_CON_RUIDO);
        for peptido in ["PALMERA", "A", "GG", "WWW", "PAL"] {
            let letras = peptido.as_bytes();
            let masas = antibioticos::masas_del_peptido(letras).unwrap();
            let (pasos, ciclico) = narrar_puntaje(letras, &masas, &espectro, true).unwrap();
            assert_eq!(
                ciclico,
                antibioticos::puntaje_ciclico(&masas, &espectro).unwrap()
            );
            assert!(!pasos.is_empty() && pasos.len() <= MAXIMO_PASOS);
            let (_, lineal) = narrar_puntaje(letras, &masas, &espectro, false).unwrap();
            assert_eq!(
                lineal,
                antibioticos::puntaje_lineal(&masas, &espectro).unwrap()
            );
        }
        let con_repetidas = masas_de("0 71 71 71 71 142 142");
        let masas = antibioticos::masas_del_peptido(b"AAA").unwrap();
        let (_, puntaje) = narrar_puntaje(b"AAA", &masas, &con_repetidas, true).unwrap();
        assert_eq!(
            puntaje,
            antibioticos::puntaje_ciclico(&masas, &con_repetidas).unwrap()
        );
    }

    #[test]
    fn narrar_el_recorte_da_lo_mismo_que_el_nucleo() {
        let entrada = del_ejemplo("antibioticos.recortar");
        let peptidos = leer_lista_de_peptidos(&entrada.texto("peptidos").unwrap()).unwrap();
        let espectro = leer_espectro(&entrada, "espectro").unwrap();
        let tabla: Vec<Vec<usize>> = peptidos.iter().map(|(_, masas)| masas.clone()).collect();
        for n in [1, 2, 3, 5, 7, 20] {
            let esperada = antibioticos::recortar(tabla.clone(), &espectro, n).unwrap();
            let (pasos, obtenida) = narrar_recorte(&peptidos, &espectro, n).unwrap();
            assert_eq!(obtenida, esperada);
            assert!(!pasos.is_empty() && pasos.len() <= MAXIMO_PASOS);
        }
    }

    #[test]
    fn narrar_los_lideres_da_lo_mismo_que_el_nucleo() {
        let espectro = masas_de(ESPECTRO_CON_RUIDO);
        for n in [1, 3, 10] {
            let mut narracion = Narracion::nueva(MAXIMO_PASOS);
            let (obtenido, ultimo) = narrar_lideres(
                &mut narracion,
                &espectro,
                n,
                &MASAS_DISTINTAS,
                "secuenciar_con_alfabeto(espectro, n, &MASAS_DISTINTAS)",
                "Prueba.",
            )
            .unwrap();
            assert_eq!(
                obtenido,
                antibioticos::secuenciar_con_lideres(&espectro, n).unwrap()
            );
            let pasos = narracion.terminar(ultimo);
            assert!(!pasos.is_empty() && pasos.len() <= MAXIMO_PASOS);
        }
    }

    #[test]
    fn narrar_la_convolucion_da_lo_mismo_que_el_nucleo() {
        for texto in [ESPECTRO_NO_ESTANDAR, ESPECTRO_PALMERA, "0 57", "0 0 0 71"] {
            let espectro = masas_de(texto);
            let mut narracion = Narracion::nueva(MAXIMO_PASOS);
            let (obtenida, ultimo) = narrar_convolucion(&mut narracion, &espectro).unwrap();
            assert_eq!(
                obtenida,
                antibioticos::convolucion_espectral(&espectro).unwrap()
            );
            let pasos = narracion.terminar(ultimo);
            assert!(!pasos.is_empty() && pasos.len() <= MAXIMO_PASOS);
        }
    }

    #[test]
    fn narrar_la_convolucion_con_lideres_da_lo_mismo_que_el_nucleo() {
        let espectro = masas_de(ESPECTRO_NO_ESTANDAR);
        for (m, n) in [(18, 40), (5, 10), (1, 5)] {
            let mut narracion = Narracion::nueva(MAXIMO_PASOS);
            let (obtenido, ultimo) =
                narrar_por_convolucion(&mut narracion, &espectro, m, n).unwrap();
            assert_eq!(
                obtenido,
                antibioticos::secuenciar_por_convolucion(&espectro, m, n).unwrap()
            );
            let pasos = narracion.terminar(ultimo);
            assert!(!pasos.is_empty() && pasos.len() <= MAXIMO_PASOS);
        }
    }

    #[test]
    fn narrar_la_autopista_da_lo_mismo_que_el_nucleo() {
        let ejemplo = del_ejemplo("antibioticos.autopista")
            .enteros("distancias")
            .unwrap();
        let sin_solucion = vec![1, 2, 4];
        let con_ceros = vec![0, -3, 2, 2, 4];
        let pequena = vec![5];
        for distancias in [ejemplo, sin_solucion, con_ceros, pequena] {
            let esperado = antibioticos::puntos_desde_distancias(&distancias).unwrap();
            let (pasos, obtenido) = narrar_autopista(&distancias).unwrap();
            assert_eq!(obtenido, esperado);
            assert!(!pasos.is_empty() && pasos.len() <= MAXIMO_PASOS);
        }
    }

    #[test]
    fn cada_entrada_narrada_tiene_pasos_con_pizarra() {
        let ids = [
            "traducir",
            "codificar",
            "masa",
            "espectro_ciclico",
            "espectro_lineal",
            "contar",
            "puntaje_ciclico",
            "puntaje_lineal",
            "recortar",
            "lideres",
            "convolucion",
            "convolucion_lideres",
            "autopista",
        ];
        for id in ids {
            let completo = format!("antibioticos.{id}");
            let pasos = crate::pasos(&completo, crate::ejemplo(&completo).unwrap()).unwrap();
            assert!(pasos.len() >= 3, "{completo}: {} pasos", pasos.len());
            assert!(pasos.len() <= MAXIMO_PASOS);
            let tamano: usize = pasos.iter().map(|paso| paso.estado.to_string().len()).sum();
            assert!(tamano < 2_000_000, "{completo}: {tamano} bytes");
            for paso in &pasos {
                assert_eq!(paso.estado["tipo"], "pizarra");
            }
        }
    }
}
