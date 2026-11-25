use crate::formato::{decimal, miles};
use crate::pizarra::{Marca, Narracion, Pizarra};
use crate::{Algoritmo, Bloque, Campo, Entrada, Paso, Salida, Tipo};
use nucleo::filogenia::{self, Arbol};
use nucleo::verosimilitud::{
    self, Calculo, Datos, Eleccion, Modelo, Referencia, Resultado, Sustitucion,
};
use std::collections::HashSet;

const MODELOS: &[&str] = &[
    "el mejor por AIC",
    "JC69",
    "JC69+Γ",
    "K80",
    "K80+Γ",
    "HKY85",
    "HKY85+Γ",
];
const REFERENCIAS: &[&str] = &["máxima verosimilitud", "unión de vecinos"];
const SEMILLA: u64 = 1;
const UMBRAL_DE_APOYO: f64 = 70.0;
const NODOS_EN_DETALLE: usize = 12;
const PATRONES_EN_TABLA: usize = 12;
const COLUMNAS_EN_DETALLE: usize = 20;
const NOMBRES_EN_GRUPO: usize = 6;

const ESPECIES: &str = "Nubaria=TCTTTCAACCCAGTCTTGTTGGGCTCGTGGACAGTGAGACATAGAATTGGTGGAAAACACGGTATTGCTGTCTTACCTAGGTAAGTGGTATAGTACATACGACAAACAAAGGTAGGGCAAGAGACTTCGTGAGTCTATGGCACCTTACGGTAGAGTAATAGCATACCTTTAGGAGCTGAAGCGTTATTGTTATTACATGCATAATGGGACCCCTTGTTAAAAAAAAAGCCTGGGAGTCCC
Ferolis=TCTTTCGACCCAATCTTGTTGAGCTCGTGTACAGTGAGACATGGAATTGATGGAAAACACGGATTTGCTATCTTGTCTAGGTAACTGGTATAGTACATACGACAAACAAAGGCAGGGCAAGAAATTTCGTGAATCTGTGGCACTCTGCGTTTTAGTAATAGCATACCTTTAGAAGCTGTAGCGTTATCGTTATTATCTGCATAATGGGACCCCTTGTTAAATAAAAACCCTGGGAGTCCC
Cardenix=TTTTTCGACCCAATCTTATTCAGCTCACGTACAGTGAGACATGGCATTGATGGAGAACACGCATTTACTGTCTAACCTAGGTAATTGGTATAGTACATACGACAAGCAAAGGTAGGGCAAAAAAATTTGTAACTCTATGGCACACTGCGATATAATAATCGCATACCTTTAAAAGTTGTAAGGTTGTCGTTATAATCCGCATGATGGTACTTCTAGTTAAACACAAATCCTGGGAGTCCC
Tavira=TCTTTCAACCCAATCTTGTTGAGCTCGCGCACAGTAAGACATGGAATTGGTGGAAAACTCGCGTTTACTGTCTGTCCTAGGTAATTGGTAAAGTACATACGACAAACAATGGCAGGGCAAAAAACTTCGTAACTCTGTGGCACACTACGATATAGTAATAGCATACCTTTAGAAGCTGTAGCGTTATCGTTATAACCTACATAATAGGACCGCTTGTTAAATACAAACCCTAGGAGTCCC
Ostrelo=TCTTTCCACCCAATCTTATTGAACTCGCGGGCAGAAAGACATGGAGTTGGTGGAAAACACGCGTTTACTGTCTGGCCTAGGTAATTGGTAAAGTACATACGACAGACAAAGGTAGGGCAAAAAACTTCGTAACTCTGTGGCACACTGCGATATAGTAATAGCATACCTTTAGAAGCTGTAGCGTTAACGTTATAACCTGCATAATGGGACCGCTTGTTAAAAACAAACCCTGGGAGTCTC
Pumbra=TCTTTCAACCCACCCCTGTTAAGCTCACGTACAGTGAGACATGGAATTGGTAGATAACGTGTGTTTACTGTCTTACCTAGCCAATTGATATAGTTTATATGACCTCCAAAAGCAGGGCAAAAAACTTCGCAACTCTATGGCACAATACGATATAGTAATAGCATACCTTTAGTAGTTGTAATGTTGACGTTATAACTTGTGGAATAGGACGTCTAGTTAAAAAAAAACGCTGGGAGTCCC
Silvana=TCTTTCTACCCACCCTTATTAAGCTCACGCACAGTGAGACATGGAATTGGTAGATAACGTGTGCTTACTTTCTTACCTAGGCAGGTGATATAGTTCATATGACCTCCAAAAGCAGGGCAACAAACTTCGTAACTCTATGGTATACTGCGATATAGTAATAGCATACCTTTAGAAGTTGTAGAGTTTACGTTATAACCTGTGGAACAGGACATCTAGTTAAAAAAAAACGCTGGGAGTCCC
Quirela=TCTTTCAACCCAACCTTGTTCAGCTCATGCATAGTAAGACATGGAATTGTTCGAAAACGGGCGGTTACTGTCTTACCTAGGTAATTATTATAGTTCATACGACAAACAACGGTAGGGCAAAAGACTTCGCTACTCTGTGGCACACTGCGATATAGTAATAGCATACCTTTAGGAGTTGTAGCGTTATCGTTATAATCTATGTAATGGGACCCATGGTTAAATAAAAGCGCTGGGAGTCCC";

const CAMPO_SECUENCIAS: Campo = Campo {
    id: "secuencias",
    etiqueta: "Secuencias alineadas (nombre=SECUENCIA por línea, o FASTA)",
    tipo: Tipo::Texto,
    ejemplo: ESPECIES,
};

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![
        Algoritmo {
            id: "filogenia.maxima_verosimilitud",
            capitulo: "filogenia",
            nombre: "Árbol de máxima verosimilitud",
            resumen: "El árbol y el modelo de sustitución que hacen más probable el alineamiento (poda de Felsenstein, JC69/K80/HKY85 con o sin Γ, elegido por AIC, y búsqueda por intercambio de vecinos como PhyML); a diferencia de unión de vecinos, usa cada columna y no solo distancias resumidas.",
            campos: &[
                CAMPO_SECUENCIAS,
                Campo {
                    id: "modelo",
                    etiqueta: "Modelo de sustitución",
                    tipo: Tipo::Opcion(MODELOS),
                    ejemplo: "el mejor por AIC",
                },
            ],
            funcion: "verosimilitud::maxima_verosimilitud",
            ejecutar: maxima_verosimilitud,
            pasos: Some(pasos_maxima_verosimilitud),
        },
        Algoritmo {
            id: "filogenia.bootstrap",
            capitulo: "filogenia",
            nombre: "Bootstrap de un árbol",
            resumen: "Cuánto apoyan los datos cada rama: se remuestrean las columnas con reemplazo, se rehace el árbol (máxima verosimilitud o unión de vecinos) en paralelo y se cuenta cuántas veces reaparece cada grupo; los árboles clásicos no dicen qué ramas son firmes.",
            campos: &[
                CAMPO_SECUENCIAS,
                Campo {
                    id: "replicas",
                    etiqueta: "Réplicas (1 a 200)",
                    tipo: Tipo::Entero,
                    ejemplo: "100",
                },
                Campo {
                    id: "arbol",
                    etiqueta: "Árbol que se rehace en cada réplica",
                    tipo: Tipo::Opcion(REFERENCIAS),
                    ejemplo: "máxima verosimilitud",
                },
            ],
            funcion: "verosimilitud::bootstrap",
            ejecutar: bootstrap,
            pasos: Some(pasos_bootstrap),
        },
    ]
}

fn recortar(texto: &str) -> String {
    let mut corto: String = texto.chars().take(40).collect();
    if texto.chars().count() > 40 {
        corto.push('…');
    }
    corto
}

fn leer_fasta_alineado(texto: &str) -> Vec<(String, Vec<u8>)> {
    let mut registros: Vec<(String, Vec<u8>)> = Vec::new();
    for linea in texto
        .lines()
        .map(str::trim)
        .filter(|linea| !linea.is_empty())
    {
        if let Some(encabezado) = linea.strip_prefix('>') {
            let nombre = encabezado.split_whitespace().next().unwrap_or_default();
            registros.push((nombre.to_string(), Vec::new()));
        } else if let Some((_, secuencia)) = registros.last_mut() {
            secuencia.extend(
                linea
                    .bytes()
                    .filter(|letra| !letra.is_ascii_whitespace())
                    .map(|letra| letra.to_ascii_uppercase()),
            );
        }
    }
    registros.retain(|(_, secuencia)| !secuencia.is_empty());
    registros
}

fn leer_alineamiento(entrada: &Entrada) -> Result<(Vec<String>, Vec<Vec<u8>>), String> {
    let texto = entrada.texto("secuencias")?;
    let mut nombres = Vec::new();
    let mut secuencias = Vec::new();
    if texto.trim_start().starts_with('>') {
        for (nombre, secuencia) in leer_fasta_alineado(&texto) {
            nombres.push(nombre);
            secuencias.push(secuencia);
            if nombres.len() > verosimilitud::MAXIMO_DE_SECUENCIAS {
                break;
            }
        }
        if nombres.is_empty() {
            return Err("El campo «secuencias» no tiene secuencias.".to_string());
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
            secuencias.push(
                valor
                    .bytes()
                    .filter(|letra| !letra.is_ascii_whitespace())
                    .map(|letra| letra.to_ascii_uppercase())
                    .collect(),
            );
            if nombres.len() > verosimilitud::MAXIMO_DE_SECUENCIAS {
                break;
            }
        }
    }
    if nombres.len() > verosimilitud::MAXIMO_DE_SECUENCIAS {
        return Err(format!(
            "Hay más de {} secuencias; ese es el máximo para la búsqueda de máxima verosimilitud.",
            verosimilitud::MAXIMO_DE_SECUENCIAS
        ));
    }
    for nombre in &nombres {
        if !filogenia::es_nombre_valido(nombre) {
            return Err(format!(
                "«{}» no es un nombre válido: usa letras, números, «_» o «.».",
                recortar(nombre)
            ));
        }
    }
    verosimilitud::revisar_alineamiento(&secuencias, &nombres)?;
    Ok((nombres, secuencias))
}

fn leer_eleccion(entrada: &Entrada) -> Result<Eleccion, String> {
    let opcion = entrada.opcion("modelo", MODELOS)?;
    let sustitucion = match opcion.trim_end_matches("+Γ") {
        "JC69" => Sustitucion::Jc69,
        "K80" => Sustitucion::K80,
        "HKY85" => Sustitucion::Hky85,
        _ => return Ok(Eleccion::MejorPorAic),
    };
    Ok(Eleccion::Fija(sustitucion, opcion.ends_with("+Γ")))
}

fn leer_referencia(entrada: &Entrada) -> Result<Referencia, String> {
    let opcion = entrada.opcion("arbol", REFERENCIAS)?;
    if opcion == REFERENCIAS[0] {
        Ok(Referencia::MaximaVerosimilitud)
    } else {
        Ok(Referencia::UnionDeVecinos)
    }
}

fn leer_replicas(entrada: &Entrada) -> Result<usize, String> {
    let replicas = entrada.natural("replicas")?;
    if replicas == 0 || replicas > verosimilitud::MAXIMO_DE_REPLICAS {
        return Err(format!(
            "Las réplicas deben ir de 1 a {}.",
            verosimilitud::MAXIMO_DE_REPLICAS
        ));
    }
    Ok(replicas)
}

fn cantidad(numero: usize, singular: &str, plural: &str) -> String {
    if numero == 1 {
        format!("1 {singular}")
    } else {
        format!("{numero} {plural}")
    }
}

fn ln(valor: f64) -> String {
    decimal(valor, 2)
}

fn corto(valor: f64) -> String {
    if valor == 0.0 {
        return "0".to_string();
    }
    if valor.abs() >= 0.001 && valor.abs() < 10_000.0 {
        let texto = decimal(valor, 4);
        return texto
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string();
    }
    format!("{valor:.3e}")
}

fn rama(valor: f64) -> String {
    filogenia::longitud_en_texto(valor)
}

fn newick(arbol: &Arbol, hojas: usize) -> String {
    verosimilitud::newick_con_apoyos(arbol, hojas, &[])
}

fn grupo(mascara: u64, nombres: &[String]) -> String {
    let dentro: Vec<&str> = nombres
        .iter()
        .enumerate()
        .filter(|(indice, _)| mascara & (1u64 << indice) != 0)
        .map(|(_, nombre)| nombre.as_str())
        .collect();
    let mut texto = dentro
        .iter()
        .take(NOMBRES_EN_GRUPO)
        .copied()
        .collect::<Vec<_>>()
        .join(", ");
    if dentro.len() > NOMBRES_EN_GRUPO {
        texto.push_str(&format!(", … (+{})", dentro.len() - NOMBRES_EN_GRUPO));
    }
    format!("{{{texto}}}")
}

fn hojas_bajo(arbol: &Arbol, nodo: usize, padre: usize, hojas: usize) -> u64 {
    if nodo < hojas {
        return 1u64 << nodo;
    }
    verosimilitud::hijos(arbol, nodo, padre)
        .iter()
        .fold(0, |mascara, &hijo| {
            mascara | hojas_bajo(arbol, hijo, nodo, hojas)
        })
}

fn descripcion_del_modelo(modelo: &Modelo) -> String {
    let mut partes = vec![modelo.nombre()];
    if modelo.sustitucion != Sustitucion::Jc69 {
        partes.push(format!("κ = {}", decimal(modelo.kappa, 2)));
    }
    if let Some(alfa) = modelo.alfa {
        partes.push(format!("α = {}", decimal(alfa, 2)));
    }
    if modelo.sustitucion == Sustitucion::Hky85 {
        let pi: Vec<String> = modelo.frecuencias.iter().map(|f| decimal(*f, 3)).collect();
        partes.push(format!("π(A, C, G, T) = {}", pi.join(", ")));
    }
    partes.join(" · ")
}

fn ordenados_por_aic(resultado: &Resultado) -> Vec<usize> {
    let mut orden: Vec<usize> = (0..resultado.comparacion.len()).collect();
    orden.sort_by(|&a, &b| {
        resultado.comparacion[a]
            .aic
            .total_cmp(&resultado.comparacion[b].aic)
    });
    orden
}

fn filas_de_modelos(resultado: &Resultado) -> Vec<Vec<String>> {
    let mejor = ordenados_por_aic(resultado)[0];
    let aic_minimo = resultado.comparacion[mejor].aic;
    resultado
        .comparacion
        .iter()
        .enumerate()
        .map(|(indice, opcion)| {
            let mut nombre = opcion.modelo.nombre();
            if indice == resultado.elegido {
                nombre.push_str(" ←");
            }
            vec![
                nombre,
                opcion.parametros.to_string(),
                ln(opcion.log_verosimilitud),
                ln(opcion.aic),
                ln(opcion.aic - aic_minimo),
            ]
        })
        .collect()
}

fn ramas_del_arbol(arbol: &Arbol, nombres: &[String]) -> Vec<Vec<String>> {
    let hojas = nombres.len();
    let mut filas = Vec::new();
    for (nodo, padre) in verosimilitud::orden_posterior(arbol) {
        let largo = verosimilitud::longitud(arbol, padre, nodo);
        let etiqueta = if padre == 0 {
            nombres[0].clone()
        } else if nodo < hojas {
            nombres[nodo].clone()
        } else {
            grupo(hojas_bajo(arbol, nodo, padre, hojas), nombres)
        };
        filas.push(vec![etiqueta, rama(largo)]);
    }
    filas
}

fn maxima_verosimilitud(entrada: &Entrada) -> Result<Salida, String> {
    let (nombres, secuencias) = leer_alineamiento(entrada)?;
    let eleccion = leer_eleccion(entrada)?;
    let resultado = verosimilitud::maxima_verosimilitud(&secuencias, &nombres, eleccion)?;
    let hojas = nombres.len();
    let mut salida = Salida::default();
    let texto_newick = newick(&resultado.arbol, hojas);
    salida.agregar(Bloque::arbol(
        &format!(
            "Árbol de máxima verosimilitud ({})",
            resultado.modelo.nombre()
        ),
        texto_newick.clone(),
    ));
    salida.agregar(Bloque::valor(
        "ln L (log-verosimilitud del árbol)",
        ln(resultado.log_verosimilitud),
    ));
    salida.agregar(Bloque::valor(
        "Modelo",
        descripcion_del_modelo(&resultado.modelo),
    ));
    salida.agregar(Bloque::tabla(
        "Comparación de modelos en el árbol inicial (← el usado)",
        &["Modelo", "Parámetros libres", "ln L", "AIC", "ΔAIC"],
        filas_de_modelos(&resultado),
    ));
    let mejor_aic = resultado.comparacion[ordenados_por_aic(&resultado)[0]].aic;
    salida.agregar(Bloque::barras(
        "ΔAIC de cada modelo (0 = el mejor)",
        "modelo",
        "ΔAIC",
        resultado
            .comparacion
            .iter()
            .map(|opcion| opcion.modelo.nombre())
            .collect(),
        resultado
            .comparacion
            .iter()
            .map(|opcion| opcion.aic - mejor_aic)
            .collect(),
    ));
    let puntos: Vec<[f64; 2]> = resultado
        .rondas
        .iter()
        .map(|ronda| [ronda.numero as f64, ronda.log_verosimilitud])
        .collect();
    salida.agregar(Bloque::linea(
        "ln L en cada ronda de intercambio de vecinos",
        "ronda",
        "ln L",
        &puntos,
        Vec::new(),
    ));
    salida.agregar(Bloque::tabla(
        "Rondas de búsqueda",
        &["Ronda", "ln L", "Intercambios que mejoraban", "Aplicados"],
        resultado
            .rondas
            .iter()
            .map(|ronda| {
                vec![
                    ronda.numero.to_string(),
                    ln(ronda.log_verosimilitud),
                    ronda.candidatos.to_string(),
                    ronda.aplicados.to_string(),
                ]
            })
            .collect(),
    ));
    salida.agregar(Bloque::tabla(
        "Ramas (sustituciones por posición)",
        &["Rama hacia", "Longitud"],
        ramas_del_arbol(&resultado.arbol, &nombres),
    ));
    salida.agregar(Bloque::arbol(
        "Árbol inicial: unión de vecinos, con ramas reajustadas por el modelo",
        newick(&resultado.arbol_de_union, hojas),
    ));
    salida.agregar(Bloque::valor(
        "Distancia de Robinson–Foulds con unión de vecinos",
        format!(
            "{} de {} posibles",
            resultado.robinson_foulds,
            2 * (hojas - 3)
        ),
    ));
    salida.agregar(Bloque::valor(
        "Posiciones y patrones distintos",
        format!(
            "{} posiciones, {} patrones",
            miles(resultado.sitios),
            miles(resultado.patrones)
        ),
    ));
    salida.agregar(Bloque::texto("Newick", texto_newick));
    explicar_maxima_verosimilitud(&mut salida, &resultado, eleccion, hojas);
    Ok(salida)
}

fn explicar_maxima_verosimilitud(
    salida: &mut Salida,
    resultado: &Resultado,
    eleccion: Eleccion,
    hojas: usize,
) {
    let orden = ordenados_por_aic(resultado);
    let mejor = &resultado.comparacion[orden[0]];
    let segundo = &resultado.comparacion[orden[1]];
    let diferencia = segundo.aic - mejor.aic;
    if diferencia < 2.0 {
        salida.explicar(format!(
            "El menor AIC es de {}, pero {} queda a solo {} unidades: con ΔAIC < 2 los dos explican los datos casi igual de bien (criterio de Burnham y Anderson).",
            mejor.modelo.nombre(),
            segundo.modelo.nombre(),
            decimal(diferencia, 2)
        ));
    } else {
        let fuerza = if diferencia > 10.0 {
            "contundente (más de 10)"
        } else {
            "clara (entre 2 y 10)"
        };
        salida.explicar(format!(
            "El AIC elige {}: el siguiente, {}, queda {} unidades atrás, una ventaja {fuerza} según Burnham y Anderson.",
            mejor.modelo.nombre(),
            segundo.modelo.nombre(),
            decimal(diferencia, 1)
        ));
    }
    if eleccion != Eleccion::MejorPorAic && resultado.elegido != orden[0] {
        salida.explicar(format!(
            "Pediste {}, que tiene un AIC {} unidades peor que el mejor: el árbol puede cambiar con un modelo mejor.",
            resultado.comparacion[resultado.elegido].modelo.nombre(),
            decimal(
                resultado.comparacion[resultado.elegido].aic - mejor.aic,
                1
            )
        ));
    }
    let modelo = &resultado.modelo;
    let mut detalles = Vec::new();
    if modelo.sustitucion != Sustitucion::Jc69 && modelo.kappa > 2.0 {
        detalles.push(format!(
            "las transiciones (A↔G, C↔T) ocurren {} veces más que cada transversión",
            decimal(modelo.kappa, 1)
        ));
    }
    if let Some(alfa) = modelo.alfa {
        if alfa < 1.0 {
            detalles.push(format!(
                "α = {} indica tasas muy desiguales: muchas posiciones casi no cambian y unas pocas cambian mucho",
                decimal(alfa, 2)
            ));
        } else if alfa > 5.0 {
            detalles.push(format!(
                "α = {} indica que todas las posiciones cambian a ritmo parecido",
                decimal(alfa, 1)
            ));
        }
    }
    if !detalles.is_empty() {
        salida.explicar(format!("En estos datos {}.", detalles.join("; y ")));
    }
    let ganancia = resultado.log_verosimilitud - resultado.log_verosimilitud_de_union;
    if resultado.robinson_foulds == 0 {
        salida.explicar(
            "La topología es la misma que da «Árbol desde secuencias» (unión de vecinos con Jukes–Cantor): ambos métodos coinciden; lo que cambia es que aquí las ramas se miden con el modelo elegido y no con distancias de Jukes–Cantor.",
        );
    } else {
        salida.explicar(format!(
            "El árbol difiere del de «Árbol desde secuencias» (unión de vecinos) en {} (Robinson–Foulds = {} de {}); la máxima verosimilitud mejoró ln L en {} sobre esa topología. Las ramas que cambian suelen ser cortas: mide su apoyo con «Bootstrap de un árbol».",
            cantidad(resultado.robinson_foulds / 2, "rama interna", "ramas internas"),
            resultado.robinson_foulds,
            2 * (hojas.max(3) - 3),
            decimal(ganancia, 2)
        ));
    }
    let rondas = resultado.rondas.len().saturating_sub(1);
    salida.explicar(format!(
        "La búsqueda hizo {} de intercambio de vecinos y paró cuando ningún intercambio subía ln L: es un óptimo local, como en PhyML; RAxML-NG (movimientos SPR) e IQ-TREE (perturbaciones al azar y varios reinicios) exploran más para escapar de ellos.",
        cantidad(rondas, "ronda", "rondas")
    ));
    let mas_larga = resultado
        .arbol
        .aristas()
        .iter()
        .map(|&(_, _, largo)| largo)
        .fold(0.0, f64::max);
    if resultado.sitios < 200 {
        salida.explicar(format!(
            "Con solo {} posiciones las longitudes y la topología tienen mucha incertidumbre; usa más secuencia si puedes.",
            resultado.sitios
        ));
    } else if mas_larga > 1.0 {
        salida.explicar(format!(
            "Hay una rama de {} sustituciones por posición: a esa distancia las secuencias están cerca de saturarse y la posición de esa rama es poco fiable.",
            decimal(mas_larga, 2)
        ));
    }
}

fn bootstrap(entrada: &Entrada) -> Result<Salida, String> {
    let (nombres, secuencias) = leer_alineamiento(entrada)?;
    let replicas = leer_replicas(entrada)?;
    let referencia = leer_referencia(entrada)?;
    let resultado = verosimilitud::bootstrap(&secuencias, &nombres, replicas, referencia, SEMILLA)?;
    let hojas = nombres.len();
    let mut salida = Salida::default();
    let texto_newick = verosimilitud::newick_con_apoyos(&resultado.arbol, hojas, &resultado.apoyos);
    let titulo = match referencia {
        Referencia::MaximaVerosimilitud => "Árbol de máxima verosimilitud con apoyo bootstrap (%)",
        Referencia::UnionDeVecinos => "Árbol de unión de vecinos con apoyo bootstrap (%)",
    };
    salida.agregar(Bloque::arbol(titulo, texto_newick.clone()));
    let mut apoyos = resultado.apoyos.clone();
    apoyos.sort_by(|a, b| b.porcentaje.total_cmp(&a.porcentaje));
    salida.agregar(Bloque::tabla(
        "Apoyo de cada rama interna",
        &["Grupo", "Réplicas que lo tienen", "Apoyo"],
        apoyos
            .iter()
            .map(|apoyo| {
                vec![
                    grupo(apoyo.particion.hojas, &nombres),
                    format!("{} de {}", apoyo.veces, replicas),
                    format!("{:.0} %", apoyo.porcentaje),
                ]
            })
            .collect(),
    ));
    salida.agregar(Bloque::barras(
        "Apoyo bootstrap por rama (%)",
        "grupo",
        "%",
        apoyos
            .iter()
            .map(|apoyo| grupo(apoyo.particion.hojas, &nombres))
            .collect(),
        apoyos.iter().map(|apoyo| apoyo.porcentaje).collect(),
    ));
    if let Some(modelo) = &resultado.modelo {
        salida.agregar(Bloque::valor("Modelo", descripcion_del_modelo(modelo)));
    }
    if let Some(valor) = resultado.log_verosimilitud {
        salida.agregar(Bloque::valor("ln L del árbol", ln(valor)));
    }
    let firmes = apoyos
        .iter()
        .filter(|apoyo| apoyo.porcentaje >= UMBRAL_DE_APOYO)
        .count();
    salida.agregar(Bloque::valor(
        "Ramas con apoyo ≥ 70 %",
        format!("{firmes} de {}", apoyos.len()),
    ));
    salida.agregar(Bloque::texto("Newick con apoyos", texto_newick));
    explicar_bootstrap(&mut salida, &apoyos, &nombres, replicas, referencia);
    Ok(salida)
}

fn explicar_bootstrap(
    salida: &mut Salida,
    apoyos: &[verosimilitud::Apoyo],
    nombres: &[String],
    replicas: usize,
    referencia: Referencia,
) {
    let firmes = apoyos
        .iter()
        .filter(|apoyo| apoyo.porcentaje >= UMBRAL_DE_APOYO)
        .count();
    salida.explicar(format!(
        "{firmes} de las {} ramas internas reaparecen en al menos el 70 % de las réplicas, el umbral usual para llamar a una rama bien apoyada (Hillis y Bull, 1993).",
        apoyos.len()
    ));
    if let Some(debil) = apoyos.last() {
        if debil.porcentaje == 0.0 {
            salida.explicar(format!(
                "La rama más débil agrupa {} y no aparece en ninguna réplica: ninguna columna la distingue (sus ramas miden casi cero), así que los datos no la resuelven.",
                grupo(debil.particion.hojas, nombres)
            ));
        } else if debil.porcentaje < UMBRAL_DE_APOYO {
            salida.explicar(format!(
                "La rama más débil agrupa {} y solo aparece en el {:.0} % de las réplicas: no la tomes como segura; cambia según qué columnas toquen en suerte.",
                grupo(debil.particion.hojas, nombres),
                debil.porcentaje
            ));
        } else if debil.porcentaje >= 95.0 {
            salida.explicar(
                "Todas las ramas superan el 95 %: con estos datos el árbol es muy estable.",
            );
        }
    }
    match referencia {
        Referencia::MaximaVerosimilitud => salida.explicar(
            "Cada réplica repitió la búsqueda de máxima verosimilitud desde su propio árbol de unión de vecinos, con los parámetros del modelo (κ, α y frecuencias) fijos en los del alineamiento original para ir más rápido.",
        ),
        Referencia::UnionDeVecinos => salida.explicar(
            "Cada réplica rehízo el árbol de unión de vecinos del capítulo 7: es mucho más rápido que máxima verosimilitud, pero usa solo distancias, así que un apoyo alto aquí no garantiza el mismo apoyo con un modelo.",
        ),
    }
    if replicas < 100 {
        salida.explicar(format!(
            "Con {replicas} réplicas cada porcentaje tiene un error de varios puntos; usa 100 o más para publicar."
        ));
    }
    salida.explicar(
        "El bootstrap mide cuánto depende una rama de las columnas que se usaron, no la probabilidad de que la rama sea verdadera.",
    );
}

fn primer_sitio_variable(datos: &Datos) -> usize {
    for (sitio, &patron) in datos.patron_de_sitio.iter().enumerate() {
        let mut vistas = HashSet::new();
        for estados in &datos.estados {
            if estados[patron] < verosimilitud::FALTANTE {
                vistas.insert(estados[patron]);
            }
        }
        if vistas.len() > 1 {
            return sitio;
        }
    }
    0
}

fn columna_en_texto(datos: &Datos, patron: usize) -> String {
    datos
        .estados
        .iter()
        .map(|estados| ['A', 'C', 'G', 'T', '-'][estados[patron].min(4) as usize])
        .collect()
}

fn categoria_central(tasas: &[f64]) -> usize {
    let mut elegida = 0;
    for (indice, tasa) in tasas.iter().enumerate() {
        if (tasa - 1.0).abs() < (tasas[elegida] - 1.0).abs() {
            elegida = indice;
        }
    }
    elegida
}

fn fila_de_vector(nombre: &str, vector: &[f64; 4]) -> Vec<String> {
    let mut fila = vec![nombre.to_string()];
    fila.extend(vector.iter().map(|&valor| corto(valor)));
    fila
}

fn vector_visible(calculo: &Calculo, nodo: usize, categoria: usize) -> [f64; 4] {
    let parcial = &calculo.abajo[nodo];
    let escala = parcial.escalas.first().copied().unwrap_or(0);
    let factor = (-(escala as f64) * 256.0 * std::f64::consts::LN_2).exp();
    parcial.valores[categoria].map(|valor| valor * factor)
}

fn pasos_maxima_verosimilitud(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (nombres, secuencias) = leer_alineamiento(entrada)?;
    let eleccion = leer_eleccion(entrada)?;
    let resultado = verosimilitud::maxima_verosimilitud(&secuencias, &nombres, eleccion)?;
    let hojas = nombres.len();
    let datos = verosimilitud::comprimir_en_patrones(&secuencias, &nombres);
    let mut narracion = Narracion::nueva(400);
    narracion.agregar(
        Pizarra::nueva()
            .variable("sitios", resultado.sitios)
            .variable("especies", hojas)
            .secuencia(&secuencias[0], &[(0, resultado.sitios.min(60), Marca::Actual)])
            .lista("Especies", nombres.clone(), &[])
            .paso(
                "let sitios = revisar_alineamiento(secuencias, nombres)?;",
                format!(
                    "Tienes {hojas} secuencias alineadas de {} posiciones. La meta: el árbol, con sus ramas y su modelo, que hace más probable ver exactamente este alineamiento.",
                    resultado.sitios
                ),
            ),
    );
    let filas_de_patrones: Vec<Vec<String>> = (0..datos.pesos.len().min(PATRONES_EN_TABLA))
        .map(|patron| {
            vec![
                (patron + 1).to_string(),
                columna_en_texto(&datos, patron),
                decimal(datos.pesos[patron], 0),
            ]
        })
        .collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("sitios", resultado.sitios)
            .variable("patrones", datos.pesos.len())
            .tabla(
                "Columnas distintas (de arriba abajo, una letra por especie)",
                &["patrón", "columna", "peso"],
                filas_de_patrones,
                &[],
            )
            .paso(
                "let datos = comprimir_en_patrones(secuencias, nombres);",
                format!(
                    "Columnas iguales dan la misma probabilidad, así que las {} posiciones se juntan en {} patrones distintos; cada patrón se calcula una vez y cuenta tantas veces como su peso.",
                    resultado.sitios,
                    datos.pesos.len()
                ),
            ),
    );
    let inicial = verosimilitud::arbol_inicial(&datos)?;
    narracion.agregar(
        Pizarra::nueva()
            .variable("patrones", datos.pesos.len())
            .variable("inicial", "unión de vecinos")
            .lista("Árbol inicial (Newick)", vec![newick(&inicial, hojas)], &[0])
            .paso(
                "let inicial = arbol_inicial(&datos)?;",
                "El punto de partida es el árbol de unión de vecinos con distancias de Jukes–Cantor, el mismo de «Árbol desde secuencias»; las ramas negativas o nulas se suben a 0.001 para poder optimizarlas.",
            ),
    );
    narrar_poda(&mut narracion, &resultado, &datos, hojas)?;
    narrar_newton(&mut narracion, &resultado, &datos, &inicial);
    narrar_modelos(&mut narracion, &resultado);
    narrar_busqueda(&mut narracion, &resultado, &datos, &nombres);
    let texto_newick = newick(&resultado.arbol, hojas);
    narracion.agregar(
        Pizarra::nueva()
            .variable("ln L", ln(resultado.log_verosimilitud))
            .variable("modelo", descripcion_del_modelo(&resultado.modelo))
            .paso(
                "let ajuste = afinar_modelo(&datos, &busqueda.arbol, inicio.modelo);",
                format!(
                    "En el árbol final se reajustan juntos las ramas, κ y α hasta que ln L deja de subir: queda en {}.",
                    ln(resultado.log_verosimilitud)
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("ln L", ln(resultado.log_verosimilitud))
            .variable("robinson_foulds", resultado.robinson_foulds)
            .lista(
                "Árbol inicial y árbol final",
                vec![newick(&resultado.arbol_de_union, hojas), texto_newick.clone()],
                &[1],
            )
            .paso(
                "let robinson_foulds = distancia_robinson_foulds(",
                if resultado.robinson_foulds == 0 {
                    "Se comparan los grupos que separa cada rama interna en ambos árboles: son los mismos (Robinson–Foulds = 0), así que unión de vecinos ya tenía la topología correcta.".to_string()
                } else {
                    format!(
                        "Se comparan los grupos que separa cada rama interna: {} grupos están en uno solo de los dos árboles (Robinson–Foulds = {}), así que la búsqueda cambió {} respecto a unión de vecinos.",
                        resultado.robinson_foulds,
                        resultado.robinson_foulds,
                        cantidad(resultado.robinson_foulds / 2, "rama", "ramas")
                    )
                },
            ),
    );
    let ultimo = Pizarra::nueva()
        .variable("ln L", ln(resultado.log_verosimilitud))
        .variable("robinson_foulds", resultado.robinson_foulds)
        .lista("Árbol de máxima verosimilitud (Newick)", vec![texto_newick], &[0])
        .paso(
            "Ok(Resultado {",
            format!(
                "Listo: árbol de máxima verosimilitud con {} y ln L = {}. En «Qué significa» se explica qué modelo ganó, cuánto se separa de unión de vecinos y cuándo confiar en él.",
                resultado.modelo.nombre(),
                ln(resultado.log_verosimilitud)
            ),
        );
    Ok(narracion.terminar(ultimo))
}

fn narrar_poda(
    narracion: &mut Narracion,
    resultado: &Resultado,
    datos: &Datos,
    hojas: usize,
) -> Result<(), String> {
    let inicio = &resultado.comparacion[resultado.elegido];
    let arbol = &inicio.arbol;
    let sitio = primer_sitio_variable(datos);
    let del_sitio = verosimilitud::datos_de_un_sitio(datos, sitio);
    let mut calculo = Calculo::nuevo(&del_sitio, inicio.modelo, arbol.nombres.len());
    calculo.recalcular_abajo(arbol);
    let categoria = categoria_central(&calculo.tasas);
    let tasa = calculo.tasas[categoria];
    let columna = columna_en_texto(&del_sitio, 0);
    let mut filas: Vec<Vec<String>> = (0..hojas)
        .map(|hoja| {
            fila_de_vector(
                &arbol.nombres[hoja],
                &vector_visible(&calculo, hoja, categoria),
            )
        })
        .collect();
    let aviso_gamma = if calculo.tasas.len() > 1 {
        format!(
            " Con Γ todo se repite para las {} tasas; aquí ves la de tasa {}.",
            calculo.tasas.len(),
            corto(tasa)
        )
    } else {
        String::new()
    };
    narracion.agregar(
        Pizarra::nueva()
            .variable("sitio", sitio + 1)
            .variable("columna", &columna)
            .variable("modelo", inicio.modelo.nombre())
            .tabla("Vectores L(x) de cada nodo", &["nodo", "A", "C", "G", "T"], filas.clone(), &[])
            .paso(
                "vector[estado as usize] = 1.0;",
                format!(
                    "Poda de Felsenstein en la posición {} (columna {columna}). Cada hoja empieza con un vector: 1 en su base y 0 en las otras (un «-» o una N valen 1 en las cuatro).{aviso_gamma}",
                    sitio + 1
                ),
            ),
    );
    let orden: Vec<(usize, usize)> = verosimilitud::orden_posterior(arbol)
        .into_iter()
        .filter(|&(nodo, _)| nodo >= hojas)
        .collect();
    if let Some(&(nodo, padre)) = orden.first() {
        let hijos = verosimilitud::hijos(arbol, nodo, padre);
        let largo = verosimilitud::longitud(arbol, nodo, hijos[0]);
        let matriz = calculo.matrices(largo)[categoria];
        let bases = ["A", "C", "G", "T"];
        let filas_de_matriz: Vec<Vec<String>> = matriz
            .iter()
            .zip(bases)
            .map(|(fila, base)| fila_de_vector(&format!("desde {base}"), fila))
            .collect();
        narracion.agregar(
            Pizarra::nueva()
                .variable("rama", rama(largo))
                .variable("tasa", corto(tasa))
                .tabla("P(t): probabilidad de terminar en cada base", &["", "A", "C", "G", "T"], filas_de_matriz, &[(0, 1), (1, 2), (2, 3), (3, 4)])
                .paso(
                    "let peso = (valor * tiempo).exp();",
                    format!(
                        "Para bajar por una rama de largo {} se usa P(t) = Σ e^(λ·t)·M, la forma cerrada de {}: la diagonal es la probabilidad de no cambiar y fuera de ella, la de cada cambio.",
                        rama(largo),
                        inicio.modelo.sustitucion.nombre()
                    ),
                ),
        );
    }
    let mut detallados = 0;
    for &(nodo, padre) in &orden {
        let hijos = verosimilitud::hijos(arbol, nodo, padre);
        if hijos.len() != 2 {
            continue;
        }
        filas.push(fila_de_vector(
            &arbol.nombres[nodo],
            &vector_visible(&calculo, nodo, categoria),
        ));
        if detallados >= NODOS_EN_DETALLE {
            continue;
        }
        detallados += 1;
        let fila = filas.len() - 1;
        narracion.agregar(
            Pizarra::nueva()
                .variable("nodo", &arbol.nombres[nodo])
                .variable("rama_primero", rama(verosimilitud::longitud(arbol, nodo, hijos[0])))
                .variable("rama_segundo", rama(verosimilitud::longitud(arbol, nodo, hijos[1])))
                .tabla(
                    "Vectores L(x) de cada nodo",
                    &["nodo", "A", "C", "G", "T"],
                    filas.clone(),
                    &[(fila, 1), (fila, 2), (fila, 3), (fila, 4)],
                )
                .paso(
                    "desde_a[0] * desde_b[0],",
                    format!(
                        "El nodo {} junta a {} y {}: para cada base x que pudo tener, multiplica la probabilidad de producir lo de abajo por un hijo (Σ P(x→y)·L(y)) por la del otro.",
                        arbol.nombres[nodo], arbol.nombres[hijos[0]], arbol.nombres[hijos[1]]
                    ),
                ),
        );
    }
    if orden.len() > detallados {
        narracion.agregar(Pizarra::nueva().variable("nodo", "…").paso(
            "for (nodo, padre) in orden_posterior(arbol) {",
            format!(
                "Los otros {} nodos internos se calculan igual, siempre de las hojas hacia arriba.",
                orden.len() - detallados
            ),
        ));
    }
    let ln_sitio = calculo.log_verosimilitud_en_raiz(arbol);
    narracion.agregar(
        Pizarra::nueva()
            .variable("sitio", sitio + 1)
            .variable("L(sitio)", corto(ln_sitio.exp()))
            .variable("ln L(sitio)", decimal(ln_sitio, 4))
            .paso(
                "let coeficientes = self.coeficientes(&self.abajo[0], &self.abajo[hijo]);",
                format!(
                    "Arriba de todo está la hoja {}: se suman sus cuatro bases pesadas por π y por la última rama. La probabilidad de esta columna es {} (ln = {}).",
                    arbol.nombres[0],
                    corto(ln_sitio.exp()),
                    decimal(ln_sitio, 4)
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("patrones", datos.pesos.len())
            .variable("ln L", ln(inicio.log_verosimilitud))
            .paso(
                "derivadas.valor += peso * sitio.ln();",
                format!(
                    "Lo mismo con los {} patrones; sus logaritmos, cada uno por su peso, se suman: ln L = {} en el árbol inicial con ese modelo.",
                    datos.pesos.len(),
                    ln(inicio.log_verosimilitud)
                ),
            ),
    );
    Ok(())
}

fn narrar_newton(narracion: &mut Narracion, resultado: &Resultado, datos: &Datos, inicial: &Arbol) {
    let modelo = resultado.comparacion[resultado.elegido].modelo;
    let mut calculo = Calculo::nuevo(datos, modelo, inicial.nombres.len());
    calculo.recalcular_abajo(inicial);
    let hijo = inicial.vecinos[0][0].0;
    let coeficientes = calculo.coeficientes(&calculo.abajo[0], &calculo.abajo[hijo]);
    let correccion = calculo.correccion_de_escala(&calculo.abajo[0], &calculo.abajo[hijo]);
    let inicio = verosimilitud::longitud(inicial, 0, hijo);
    let (optima, _) = calculo.optimizar_longitud(&coeficientes, inicio);
    let mut actual_rama = inicio;
    for iteracion in 0..8 {
        let actual = calculo.evaluar_rama(&coeficientes, actual_rama);
        narracion.agregar(
            Pizarra::nueva()
                .variable("rama", rama(actual_rama))
                .variable("actual.valor", ln(actual.valor + correccion))
                .variable("actual.primera", corto(actual.primera))
                .variable("actual.segunda", corto(actual.segunda))
                .paso(
                    "rama - actual.primera / actual.segunda",
                    if iteracion == 0 {
                        format!(
                            "Las ramas se optimizan una por una con Newton–Raphson. Primera rama ({}): empieza en {}, con pendiente {} y curvatura {}; el siguiente largo es rama − pendiente / curvatura.",
                            inicial.nombres[0],
                            rama(actual_rama),
                            corto(actual.primera),
                            corto(actual.segunda)
                        )
                    } else {
                        format!(
                            "Ahora la rama mide {} y ln L = {}: la pendiente ({}) se acerca a 0, señal de que estamos cerca del máximo.",
                            rama(actual_rama),
                            ln(actual.valor + correccion),
                            corto(actual.primera)
                        )
                    },
                ),
        );
        if (actual_rama - optima).abs() < 1e-6 * (1.0 + optima) || actual.segunda >= 0.0 {
            break;
        }
        actual_rama = (actual_rama - actual.primera / actual.segunda).clamp(1e-8, 10.0);
    }
    narracion.agregar(
        Pizarra::nueva()
            .variable("rama", rama(optima))
            .paso(
                "let nuevo = self.una_pasada_de_ramas(arbol);",
                format!(
                    "Esa rama queda en {}. Una pasada recorre todas las ramas así y se repite hasta que ln L sube menos de 0.01.",
                    rama(optima)
                ),
            ),
    );
}

fn narrar_modelos(narracion: &mut Narracion, resultado: &Resultado) {
    let inicio = &resultado.comparacion[resultado.elegido];
    if inicio.modelo.sustitucion != Sustitucion::Jc69 {
        narracion.agregar(
            Pizarra::nueva()
                .variable("kappa", decimal(inicio.modelo.kappa, 3))
                .paso(
                    "calculo.optimizar_kappa(&arbol);",
                    format!(
                        "κ (cuánto más frecuentes son las transiciones) se busca con el método de Brent: prueba valores, ajusta una parábola y se queda con el que da mayor ln L: κ = {}.",
                        decimal(inicio.modelo.kappa, 3)
                    ),
                ),
        );
    }
    if let Some(alfa) = inicio.modelo.alfa {
        let tasas = inicio.modelo.tasas();
        narracion.agregar(
            Pizarra::nueva()
                .variable("alfa", decimal(alfa, 3))
                .lista(
                    "Tasas de las 4 categorías Γ",
                    tasas.iter().map(|&tasa| corto(tasa)).collect(),
                    &[],
                )
                .paso(
                    "tasas.push((acumulado - acumulado_anterior) * categorias as f64);",
                    format!(
                        "α = {} da la forma de la gamma; se corta en 4 partes igual de probables y cada categoría usa la tasa media de su parte (Yang, 1994).",
                        decimal(alfa, 3)
                    ),
                ),
        );
    }
    let filas = filas_de_modelos(resultado);
    let mejor = ordenados_por_aic(resultado)[0];
    narracion.agregar(
        Pizarra::nueva()
            .variable("elegido", resultado.comparacion[resultado.elegido].modelo.nombre())
            .tabla(
                "Seis modelos ajustados en el árbol inicial",
                &["Modelo", "k", "ln L", "AIC", "ΔAIC"],
                filas,
                &[(mejor, 3), (mejor, 4)],
            )
            .paso(
                "2.0 * parametros as f64 - 2.0 * log_verosimilitud",
                format!(
                    "Cada modelo se ajusta igual y se califica con AIC = 2k − 2 ln L (k = ramas + parámetros del modelo): más parámetros solo ganan si suben ln L lo suficiente. El menor AIC es de {}.",
                    resultado.comparacion[mejor].modelo.nombre()
                ),
            ),
    );
}

fn narrar_busqueda(
    narracion: &mut Narracion,
    resultado: &Resultado,
    datos: &Datos,
    nombres: &[String],
) {
    let hojas = nombres.len();
    let inicio = &resultado.comparacion[resultado.elegido];
    let mut arbol = inicio.arbol.clone();
    let mut calculo = Calculo::nuevo(datos, inicio.modelo, arbol.nombres.len());
    calculo.ajustar_ramas(&mut arbol, verosimilitud::MAXIMO_DE_PASADAS);
    let candidatos = calculo.intercambios_que_mejoran(&arbol);
    let mut filas = vec![vec![
        "0".to_string(),
        ln(resultado.rondas[0].log_verosimilitud),
        "—".to_string(),
        "—".to_string(),
    ]];
    narracion.agregar(
        Pizarra::nueva()
            .variable("ronda", 0)
            .variable("actual", ln(resultado.rondas[0].log_verosimilitud))
            .tabla("Rondas", &["ronda", "ln L", "mejoran", "aplicados"], filas.clone(), &[(0, 1)])
            .paso(
                "let mut actual = calculo.ajustar_ramas(&mut arbol, MAXIMO_DE_PASADAS);",
                format!(
                    "Empieza la búsqueda de topología con el modelo elegido: en el árbol inicial ln L = {}.",
                    ln(resultado.rondas[0].log_verosimilitud)
                ),
            ),
    );
    let filas_de_candidatos: Vec<Vec<String>> = candidatos
        .iter()
        .take(6)
        .map(|candidato| {
            let (nodo, hijo) = candidato.extremos;
            let (hermano, nieto) = candidato.subarboles;
            vec![
                grupo(hojas_bajo(&arbol, hermano, nodo, hojas), nombres),
                grupo(hojas_bajo(&arbol, nieto, hijo, hojas), nombres),
                format!("+{}", decimal(candidato.ganancia, 2)),
            ]
        })
        .collect();
    let nota = if candidatos.is_empty() {
        "En cada rama interna hay dos maneras de intercambiar los subárboles de sus extremos; ninguna sube ln L, así que el árbol inicial ya es un óptimo local.".to_string()
    } else {
        format!(
            "En cada rama interna se prueban los dos intercambios de vecinos posibles, reoptimizando solo la rama central; {} ln L. La mejor ganancia es +{}.",
            if candidatos.len() == 1 {
                "uno sube".to_string()
            } else {
                format!("{} suben", candidatos.len())
            },
            decimal(candidatos[0].ganancia, 2)
        )
    };
    narracion.agregar(
        Pizarra::nueva()
            .variable("ronda", 1)
            .variable("candidatos", candidatos.len())
            .tabla(
                "Intercambios que mejoran (ronda 1)",
                &["sale de un lado", "sale del otro", "ganancia en ln L"],
                filas_de_candidatos,
                &[(0, 2)],
            )
            .paso(
                "let candidatos = calculo.intercambios_que_mejoran(&arbol);",
                nota,
            ),
    );
    for (anterior, ronda) in resultado.rondas.iter().zip(resultado.rondas.iter().skip(1)) {
        filas.push(vec![
            ronda.numero.to_string(),
            ln(ronda.log_verosimilitud),
            ronda.candidatos.to_string(),
            ronda.aplicados.to_string(),
        ]);
        let fila = filas.len() - 1;
        narracion.agregar(
            Pizarra::nueva()
                .variable("ronda", ronda.numero)
                .variable("actual", ln(ronda.log_verosimilitud))
                .tabla("Rondas", &["ronda", "ln L", "mejoran", "aplicados"], filas.clone(), &[(fila, 1)])
                .paso(
                    "aplicar_intercambios(&mut arbol, &elegidos);",
                    format!(
                        "Ronda {}: de {} que mejoraban se aplican {} que no tocan los mismos nodos, se reoptimizan las ramas y ln L sube de {} a {}.",
                        ronda.numero,
                        cantidad(ronda.candidatos, "intercambio", "intercambios"),
                        ronda.aplicados,
                        ln(anterior.log_verosimilitud),
                        ln(ronda.log_verosimilitud)
                    ),
                ),
        );
    }
    let numero = resultado.rondas.len();
    narracion.agregar(
        Pizarra::nueva()
            .variable("ronda", numero)
            .variable("candidatos", 0)
            .tabla("Rondas", &["ronda", "ln L", "mejoran", "aplicados"], filas, &[])
            .paso(
                "if candidatos.is_empty() {",
                format!(
                    "En la ronda {numero} ningún intercambio sube ln L en más de 0.001: la búsqueda se detiene en este árbol."
                ),
            ),
    );
}

fn pasos_bootstrap(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (nombres, secuencias) = leer_alineamiento(entrada)?;
    let replicas = leer_replicas(entrada)?;
    let referencia = leer_referencia(entrada)?;
    let resultado = verosimilitud::bootstrap(&secuencias, &nombres, replicas, referencia, SEMILLA)?;
    let hojas = nombres.len();
    let datos = verosimilitud::comprimir_en_patrones(&secuencias, &nombres);
    let sitios = datos.patron_de_sitio.len();
    let mut narracion = Narracion::nueva(400);
    narracion.agregar(
        Pizarra::nueva()
            .variable("replicas", replicas)
            .variable("especies", hojas)
            .variable("sitios", sitios)
            .paso(
                "revisar_alineamiento(secuencias, nombres)?;",
                format!(
                    "{hojas} especies y {sitios} posiciones. Idea del bootstrap (Felsenstein, 1985): fingir {replicas} alineamientos nuevos sorteando columnas, rehacer el árbol con cada uno y ver qué ramas se repiten."
                ),
            ),
    );
    let grupos: Vec<String> = resultado
        .apoyos
        .iter()
        .map(|apoyo| grupo(apoyo.particion.hojas, &nombres))
        .collect();
    let (linea, nota) = match referencia {
        Referencia::MaximaVerosimilitud => (
            "let resultado = maxima_verosimilitud(secuencias, nombres, Eleccion::MejorPorAic)?;",
            format!(
                "Primero el árbol de referencia con los datos originales: máxima verosimilitud con {}. Sus {} ramas internas separan estos grupos.",
                resultado
                    .modelo
                    .map(|modelo| modelo.nombre())
                    .unwrap_or_default(),
                grupos.len()
            ),
        ),
        Referencia::UnionDeVecinos => (
            "Referencia::UnionDeVecinos => (arbol_de_union(&datos)?, None, None),",
            format!(
                "Primero el árbol de referencia con los datos originales: unión de vecinos con distancias de Jukes–Cantor. Sus {} ramas internas separan estos grupos.",
                grupos.len()
            ),
        ),
    };
    narracion.agregar(
        Pizarra::nueva()
            .variable("replicas", replicas)
            .lista("Grupos del árbol de referencia", grupos.clone(), &[])
            .paso(linea, nota),
    );
    let columnas = verosimilitud::columnas_de_replica(sitios, 0, SEMILLA);
    let mostradas: Vec<String> = columnas
        .iter()
        .take(COLUMNAS_EN_DETALLE)
        .map(|columna| (columna + 1).to_string())
        .collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("replica", 1)
            .variable("sitios", sitios)
            .lista("Columnas sorteadas para la réplica 1", mostradas, &[0])
            .paso(
                "(0..sitios).map(|_| azar.entero(sitios)).collect()",
                format!(
                    "La réplica 1 sortea {sitios} columnas con reemplazo: su primera columna es la {} del original, la segunda la {}, y así.",
                    columnas[0] + 1,
                    columnas.get(1).map_or(columnas[0], |columna| *columna) + 1
                ),
            ),
    );
    let mut veces_por_columna = vec![0usize; sitios];
    for &columna in &columnas {
        veces_por_columna[columna] += 1;
    }
    let ausentes = veces_por_columna
        .iter()
        .filter(|&&veces| veces == 0)
        .count();
    let filas_de_columnas: Vec<Vec<String>> = veces_por_columna
        .iter()
        .take(COLUMNAS_EN_DETALLE)
        .enumerate()
        .map(|(columna, veces)| vec![(columna + 1).to_string(), veces.to_string()])
        .collect();
    let muestra = verosimilitud::remuestrear(&datos, &columnas);
    narracion.agregar(
        Pizarra::nueva()
            .variable("replica", 1)
            .variable("ausentes", ausentes)
            .variable("patrones", muestra.pesos.len())
            .tabla("Veces que salió cada columna original", &["columna", "veces"], filas_de_columnas, &[])
            .paso(
                "cuentas[datos.patron_de_sitio[columna]] += 1.0;",
                format!(
                    "Algunas columnas salen dos o tres veces y otras ninguna: {ausentes} de {sitios} quedaron fuera (lo esperado es cerca del 37 %). Así cambian los pesos de los patrones."
                ),
            ),
    );
    let modelo = resultado.modelo;
    let arbol_de_la_replica = verosimilitud::arbol_de_replica(&datos, 0, SEMILLA, modelo)?;
    let linea_del_arbol = match referencia {
        Referencia::MaximaVerosimilitud => {
            "Ok(buscar_con_intercambios(&muestra, &arbol, modelo).arbol)"
        }
        Referencia::UnionDeVecinos => {
            "let mut arbol = filogenia::union_de_vecinos(&distancias, &muestra.nombres)?;"
        }
    };
    narracion.agregar(
        Pizarra::nueva()
            .variable("replica", 1)
            .lista(
                "Árbol de la réplica 1 (Newick)",
                vec![newick(&arbol_de_la_replica, hojas)],
                &[0],
            )
            .paso(
                linea_del_arbol,
                "Con esas columnas se rehace el árbol entero, con el mismo método que el de referencia.",
            ),
    );
    let de_la_primera: HashSet<u64> =
        verosimilitud::particiones_resueltas(&arbol_de_la_replica, hojas)
            .iter()
            .map(|particion| particion.hojas)
            .collect();
    let filas_de_presencia: Vec<Vec<String>> = resultado
        .apoyos
        .iter()
        .zip(&grupos)
        .map(|(apoyo, nombre)| {
            let esta = de_la_primera.contains(&apoyo.particion.hojas);
            vec![nombre.clone(), if esta { "sí" } else { "no" }.to_string()]
        })
        .collect();
    let presentes = filas_de_presencia
        .iter()
        .filter(|fila| fila[1] == "sí")
        .count();
    narracion.agregar(
        Pizarra::nueva()
            .variable("replica", 1)
            .variable("presentes", presentes)
            .tabla("¿Está el grupo en la réplica 1?", &["grupo", "¿está?"], filas_de_presencia, &[])
            .paso(
                ".filter(|de_la_replica| de_la_replica.contains(&particion.hojas))",
                format!(
                    "Cada rama del árbol de referencia se busca en el de la réplica como un grupo de hojas (una rama de longitud casi cero no cuenta, porque no está resuelta): {presentes} de {} reaparecen.",
                    grupos.len()
                ),
            ),
    );
    let mut veces = vec![0usize; resultado.apoyos.len()];
    for (indice, particiones) in resultado.particiones_por_replica.iter().enumerate() {
        for (cuenta, apoyo) in veces.iter_mut().zip(&resultado.apoyos) {
            if particiones.contains(&apoyo.particion.hojas) {
                *cuenta += 1;
            }
        }
        let numero = indice + 1;
        let contar = (2..=5).contains(&numero) || numero % 10 == 0;
        if !contar || numero == replicas {
            continue;
        }
        let filas: Vec<Vec<String>> = grupos
            .iter()
            .zip(&veces)
            .map(|(nombre, cuenta)| {
                vec![
                    nombre.clone(),
                    cuenta.to_string(),
                    format!("{:.0} %", 100.0 * *cuenta as f64 / numero as f64),
                ]
            })
            .collect();
        narracion.agregar(
            Pizarra::nueva()
                .variable("replica", numero)
                .tabla("Cuentas hasta ahora", &["grupo", "veces", "%"], filas, &[])
                .paso(
                    "let arbol = arbol_de_replica(&datos, replica, semilla, modelo)?;",
                    if numero == 2 {
                        "Las demás réplicas se hacen igual, en paralelo en todos los núcleos y cada una con su propia semilla. Tras 2 réplicas así van las cuentas.".to_string()
                    } else {
                        format!("Tras {numero} réplicas así van las cuentas: los porcentajes se van estabilizando.")
                    },
                ),
        );
    }
    let filas_finales: Vec<Vec<String>> = resultado
        .apoyos
        .iter()
        .zip(&grupos)
        .map(|(apoyo, nombre)| {
            vec![
                nombre.clone(),
                apoyo.veces.to_string(),
                format!("{:.0} %", apoyo.porcentaje),
            ]
        })
        .collect();
    let resaltar: Vec<(usize, usize)> = resultado
        .apoyos
        .iter()
        .enumerate()
        .filter(|(_, apoyo)| apoyo.porcentaje >= UMBRAL_DE_APOYO)
        .map(|(fila, _)| (fila, 2))
        .collect();
    let firmes = resaltar.len();
    narracion.agregar(
        Pizarra::nueva()
            .variable("replica", replicas)
            .variable("firmes", firmes)
            .tabla("Apoyo final", &["grupo", "veces", "%"], filas_finales, &resaltar)
            .paso(
                "porcentaje: 100.0 * veces as f64 / replicas as f64,",
                format!(
                    "Tras las {replicas} réplicas, el apoyo de cada rama es el porcentaje de réplicas que la tienen; {firmes} de {} llegan al 70 % (resaltadas).",
                    grupos.len()
                ),
            ),
    );
    let ultimo = Pizarra::nueva()
        .variable("firmes", firmes)
        .lista(
            "Newick con apoyos",
            vec![verosimilitud::newick_con_apoyos(&resultado.arbol, hojas, &resultado.apoyos)],
            &[0],
        )
        .paso(
            "let apoyos = contar_apoyos(&particiones(&arbol, hojas), &particiones_por_replica);",
            "Los porcentajes se escriben como etiquetas de los nodos internos del Newick. En «Qué significa» se explica qué ramas son firmes y cuáles no.",
        );
    Ok(narracion.terminar(ultimo))
}

#[cfg(test)]
mod pruebas {
    use super::*;

    fn entrada_de_ejemplo(extra: &[(&str, &str)]) -> Entrada {
        let mut pares = vec![("secuencias", ESPECIES)];
        pares.extend_from_slice(extra);
        Entrada::nueva(&pares)
    }

    #[test]
    fn la_narracion_termina_en_el_mismo_arbol() {
        let entrada = entrada_de_ejemplo(&[("modelo", "el mejor por AIC")]);
        let salida = maxima_verosimilitud(&entrada).unwrap();
        let Bloque::Arbol {
            newick: esperado, ..
        } = &salida.bloques[0]
        else {
            panic!("el primer bloque debe ser el árbol");
        };
        let pasos = pasos_maxima_verosimilitud(&entrada).unwrap();
        let ultimo = pasos.last().unwrap();
        assert_eq!(ultimo.estado["lista"]["elementos"][0], esperado.as_str());
        assert!(pasos.len() <= 400 && pasos.len() > 10);
    }

    #[test]
    fn la_poda_narrada_da_la_verosimilitud_del_arbol() {
        let (nombres, secuencias) = leer_alineamiento(&entrada_de_ejemplo(&[])).unwrap();
        let resultado =
            verosimilitud::maxima_verosimilitud(&secuencias, &nombres, Eleccion::MejorPorAic)
                .unwrap();
        let datos = verosimilitud::comprimir_en_patrones(&secuencias, &nombres);
        let inicio = &resultado.comparacion[resultado.elegido];
        let mut total = 0.0;
        for sitio in 0..datos.patron_de_sitio.len() {
            let del_sitio = verosimilitud::datos_de_un_sitio(&datos, sitio);
            let mut calculo = Calculo::nuevo(&del_sitio, inicio.modelo, inicio.arbol.nombres.len());
            total += calculo.log_verosimilitud(&inicio.arbol);
        }
        assert!((total - inicio.log_verosimilitud).abs() < 1e-6);
    }

    #[test]
    fn el_newton_narrado_llega_a_la_rama_optima() {
        let entrada = entrada_de_ejemplo(&[("modelo", "HKY85+Γ")]);
        let pasos = pasos_maxima_verosimilitud(&entrada).unwrap();
        let newton: Vec<&Paso> = pasos
            .iter()
            .filter(|paso| paso.linea == "rama - actual.primera / actual.segunda")
            .collect();
        let final_ = pasos
            .iter()
            .find(|paso| paso.linea == "let nuevo = self.una_pasada_de_ramas(arbol);")
            .unwrap();
        let ultima = &newton.last().unwrap().estado["variables"][0][1];
        assert_eq!(ultima, &final_.estado["variables"][0][1]);
    }

    #[test]
    fn el_bootstrap_narrado_cuenta_lo_mismo() {
        for referencia in REFERENCIAS {
            let entrada = entrada_de_ejemplo(&[("replicas", "30"), ("arbol", referencia)]);
            let salida = bootstrap(&entrada).unwrap();
            let pasos = pasos_bootstrap(&entrada).unwrap();
            let ultimo = pasos.last().unwrap();
            let Bloque::Arbol { newick, .. } = &salida.bloques[0] else {
                panic!("el primer bloque debe ser el árbol");
            };
            assert_eq!(ultimo.estado["lista"]["elementos"][0], newick.as_str());
            let (nombres, secuencias) = leer_alineamiento(&entrada).unwrap();
            let referencia = leer_referencia(&entrada).unwrap();
            let resultado =
                verosimilitud::bootstrap(&secuencias, &nombres, 30, referencia, SEMILLA).unwrap();
            let datos = verosimilitud::comprimir_en_patrones(&secuencias, &nombres);
            let primera =
                verosimilitud::arbol_de_replica(&datos, 0, SEMILLA, resultado.modelo).unwrap();
            let mut esperadas: Vec<u64> =
                verosimilitud::particiones_resueltas(&primera, nombres.len())
                    .iter()
                    .map(|particion| particion.hojas)
                    .collect();
            let mut obtenidas = resultado.particiones_por_replica[0].clone();
            esperadas.sort();
            obtenidas.sort();
            assert_eq!(esperadas, obtenidas);
        }
    }

    #[test]
    fn el_fasta_conserva_los_huecos_y_los_datos_faltantes() {
        let fasta = ">a uno\nAC?T\nGG-A\n>b\nACGTGGTA\n>c\nAC-TNGCA\n>d\nTCGAGGCA\n";
        let entrada = Entrada::nueva(&[("secuencias", fasta)]);
        let (nombres, secuencias) = leer_alineamiento(&entrada).unwrap();
        assert_eq!(nombres, ["a", "b", "c", "d"]);
        assert_eq!(secuencias[0], b"AC?TGG-A");
        assert!(secuencias.iter().all(|secuencia| secuencia.len() == 8));
        assert!(leer_alineamiento(&Entrada::nueva(&[("secuencias", ">solo\n")])).is_err());
    }

    #[test]
    fn sin_senal_el_bootstrap_dice_que_los_datos_no_resuelven_la_rama() {
        let iguales = "a=ACGTACGT\nb=ACGTACGT\nc=ACGTACGT\nd=ACGTACGT\ne=ACGTACGT";
        for referencia in REFERENCIAS {
            let entrada = Entrada::nueva(&[
                ("secuencias", iguales),
                ("replicas", "20"),
                ("arbol", referencia),
            ]);
            let texto = bootstrap(&entrada).unwrap().interpretacion.join(" ");
            assert!(texto.contains("no aparece en ninguna réplica"), "{texto}");
            assert!(texto.starts_with("0 de las 2 ramas"), "{texto}");
        }
    }

    #[test]
    fn el_ejemplo_muestra_lo_que_promete() {
        let salida =
            maxima_verosimilitud(&entrada_de_ejemplo(&[("modelo", "el mejor por AIC")])).unwrap();
        let texto = salida.interpretacion.join(" ");
        assert!(
            texto.contains("difiere del de «Árbol desde secuencias»"),
            "{texto}"
        );
        let apoyo = bootstrap(&entrada_de_ejemplo(&[
            ("replicas", "100"),
            ("arbol", "máxima verosimilitud"),
        ]))
        .unwrap();
        let texto = apoyo.interpretacion.join(" ");
        assert!(texto.contains("La rama más débil"), "{texto}");
    }
}
