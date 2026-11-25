use crate::formato::{cientifico, decimal, miles, porcentaje};
use crate::pizarra::{Marca, Narracion, Pizarra};
use crate::{Algoritmo, Bloque, Campo, Entrada, Paso, Salida, Tipo};
use nucleo::motivos_em::{self, Modelo, MotivoEm, Opciones, Parametros, Sitio};
use nucleo::{ecoli, motivos, secuencia};
use rayon::prelude::*;
use std::time::Instant;

const REGIONES_SOS_Y_CONTROLES: &str = "@ecoli:4257015-4257114
@ecoli:2823770-2823869
@ecoli:1020920-1021019
@ecoli:1230667-1230766
@ecoli:2751695-2751794
@ecoli:3997883-3997982
@ecoli:3817660-3817759
@ecoli:2081263-2081362
@ecoli:3212947-3213046
@ecoli:1292923-1293022";

const MODELOS: &[&str] = &["OOPS", "ZOOPS"];
const HEBRAS: &[&str] = &["una", "ambas"];
const LETRAS: [&str; 4] = ["A", "C", "G", "T"];
const CAJA_LEXA: &str = "CTGTATATATATACAG";
const CAJA_LEXA_EXTENDIDA: &[u8] = b"TACTGTATATATATACAGTA";

const K_MAXIMA: usize = 50;
const MAXIMO_DE_SECUENCIAS: usize = 1000;
const PSEUDOCONTEO: f64 = 1.0;
const ARRANQUES_MAXIMOS: usize = 1000;
const ARRANQUES_MINIMOS: usize = 50;
const PRESUPUESTO: f64 = 1e9;
const PRESUPUESTO_CLASICOS: f64 = 1.5e9;
const GIBBS_ITERACIONES: usize = 500;
const GIBBS_REINICIOS: usize = 20;
const GIBBS_SEMILLA: u64 = 7;
const UMBRAL_POSTERIOR: f64 = 0.5;
const AZAR_TOLERADO: f64 = 0.001;
const DISTANCIA_A_UN_GEN: usize = 300;
const MAXIMO_PASOS: usize = 400;
const SECUENCIAS_ANIMABLES: usize = 20;
const BASES_ANIMABLES: usize = 3000;
const SEMILLAS_EN_DETALLE: usize = 3;
const MEJORAS_MOSTRADAS: usize = 12;
const ITERACIONES_EN_DETALLE: usize = 10;
const CADA_CUANTAS_ITERACIONES: usize = 10;

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![Algoritmo {
        id: "motivos.meme",
        capitulo: "motivos",
        nombre: "MEME: motivos por maximización de la esperanza",
        resumen: "En vez de elegir un k-mer por secuencia como Gibbs o la búsqueda codiciosa, reparte la probabilidad del sitio entre todas las ventanas (modelos OOPS y ZOOPS, en una o ambas hebras) y ajusta la matriz de pesos por EM desde los k-mers más prometedores.",
        campos: &[
            Campo {
                id: "adn",
                etiqueta: "Secuencias de ADN",
                tipo: Tipo::Secuencias,
                ejemplo: REGIONES_SOS_Y_CONTROLES,
            },
            Campo {
                id: "k",
                etiqueta: "Longitud del motivo (k)",
                tipo: Tipo::Entero,
                ejemplo: "16",
            },
            Campo {
                id: "modelo",
                etiqueta: "Modelo (OOPS: un sitio por secuencia; ZOOPS: cero o uno)",
                tipo: Tipo::Opcion(MODELOS),
                ejemplo: "ZOOPS",
            },
            Campo {
                id: "hebras",
                etiqueta: "Hebras donde buscar",
                tipo: Tipo::Opcion(HEBRAS),
                ejemplo: "ambas",
            },
        ],
        funcion: "motivos_em::meme",
        ejecutar: meme,
        pasos: Some(pasos_meme),
    }]
}

struct Conjunto {
    textos: Vec<Vec<u8>>,
    nombres: Vec<String>,
    inicios: Vec<Option<usize>>,
}

struct Pedido {
    conjunto: Conjunto,
    opciones: Opciones,
    recortados: bool,
}

fn texto(secuencia: &[u8]) -> String {
    String::from_utf8_lossy(secuencia).into_owned()
}

fn recortar(secuencia: &[u8], maximo: usize) -> String {
    if secuencia.len() <= maximo {
        texto(secuencia)
    } else {
        format!("{}…", texto(&secuencia[..maximo]))
    }
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

fn leer_conjunto(entrada: &Entrada) -> Result<Conjunto, String> {
    let registros = entrada.registros("adn")?;
    if registros.len() < 2 {
        return Err(format!(
            "El campo «adn» necesita al menos 2 secuencias y tiene {}.",
            registros.len()
        ));
    }
    if registros.len() > MAXIMO_DE_SECUENCIAS {
        return Err(format!(
            "El campo «adn» tiene {} secuencias; el máximo es {}.",
            miles(registros.len()),
            miles(MAXIMO_DE_SECUENCIAS)
        ));
    }
    for registro in &registros {
        if let Some(letra) = registro
            .secuencia
            .iter()
            .find(|letra| !matches!(letra, b'A' | b'C' | b'G' | b'T'))
        {
            return Err(format!(
                "La secuencia «{}» solo puede tener A, C, G y T, pero contiene «{}».",
                registro.nombre, *letra as char
            ));
        }
    }
    let crudo = entrada.texto("adn")?;
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
        let (pieza, por_defecto) = if es_fasta {
            ("", registro.nombre.clone())
        } else {
            (
                piezas.get(indice).copied().unwrap_or(""),
                format!("s{}", indice + 1),
            )
        };
        let (nombre, inicio) = describir_origen(pieza, por_defecto);
        conjunto.textos.push(registro.secuencia);
        conjunto.nombres.push(nombre);
        conjunto.inicios.push(inicio);
    }
    Ok(conjunto)
}

fn leer_k(entrada: &Entrada, textos: &[Vec<u8>]) -> Result<usize, String> {
    let k = entrada.natural("k")?;
    if k < 2 {
        return Err("k debe ser al menos 2.".to_string());
    }
    if k > K_MAXIMA {
        return Err(format!(
            "k = {k} es demasiado largo: los motivos regulatorios suelen medir entre 6 y 30 pb y aquí el máximo es {K_MAXIMA}."
        ));
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

fn costo_por_iteracion(textos: &[Vec<u8>], k: usize, ambas_hebras: bool) -> f64 {
    let ventanas: usize = textos.iter().map(|texto| texto.len() + 1 - k).sum();
    let hebras = if ambas_hebras { 2.0 } else { 1.0 };
    2.0 * ventanas as f64 * hebras * k as f64
}

fn leer_pedido(entrada: &Entrada) -> Result<Pedido, String> {
    let conjunto = leer_conjunto(entrada)?;
    let k = leer_k(entrada, &conjunto.textos)?;
    let modelo = match entrada.opcion("modelo", MODELOS)?.as_str() {
        "ZOOPS" => Modelo::Zoops,
        _ => Modelo::Oops,
    };
    let ambas_hebras = entrada.opcion("hebras", HEBRAS)? == "ambas";
    let costo = costo_por_iteracion(&conjunto.textos, k, ambas_hebras);
    let iteraciones_por_arranque = (motivos_em::ITERACIONES_DE_PRUEBA + 1) as f64;
    let posibles = PRESUPUESTO / (costo * iteraciones_por_arranque);
    if posibles < ARRANQUES_MINIMOS as f64 {
        return Err(format!(
            "Esta entrada pediría unas {} operaciones por iteración y no alcanzaría para probar ni {ARRANQUES_MINIMOS} puntos de arranque en unos segundos. Usa menos secuencias, secuencias más cortas o un k menor.",
            cientifico(costo)
        ));
    }
    let arranques = (posibles as usize).min(ARRANQUES_MAXIMOS);
    let distintos: usize = conjunto
        .textos
        .iter()
        .map(|texto| texto.len() + 1 - k)
        .sum();
    Ok(Pedido {
        conjunto,
        opciones: Opciones {
            k,
            modelo,
            ambas_hebras,
            pseudoconteo: PSEUDOCONTEO,
            arranques,
        },
        recortados: arranques < ARRANQUES_MAXIMOS && arranques < distintos,
    })
}

fn nombre_del_modelo(modelo: Modelo) -> &'static str {
    match modelo {
        Modelo::Oops => "OOPS",
        Modelo::Zoops => "ZOOPS",
    }
}

fn sitio_orientado(textos: &[Vec<u8>], indice: usize, sitio: &Sitio, k: usize) -> Vec<u8> {
    let ventana = &textos[indice][sitio.posicion..sitio.posicion + k];
    if sitio.hebra == b'-' {
        secuencia::complemento_inverso(ventana)
    } else {
        ventana.to_vec()
    }
}

fn coordenada(inicio: Option<usize>, posicion: usize) -> String {
    match inicio {
        Some(inicio) => miles(inicio + posicion),
        None => "—".to_string(),
    }
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

fn menor_distancia(corta: &[u8], larga: &[u8]) -> usize {
    let mut menor = corta.len();
    for ventana in larga.windows(corta.len()) {
        menor = menor.min(motivos::distancia_hamming(corta, ventana));
    }
    menor
}

fn parecido_a_lexa(consenso: &[u8]) -> Option<(usize, usize)> {
    let k = consenso.len();
    if k < 8 {
        return None;
    }
    let inverso = secuencia::complemento_inverso(consenso);
    let (menor, comparadas, comparaciones) = if k <= CAJA_LEXA_EXTENDIDA.len() {
        let menor = menor_distancia(consenso, CAJA_LEXA_EXTENDIDA)
            .min(menor_distancia(&inverso, CAJA_LEXA_EXTENDIDA));
        (menor, k, 2 * (CAJA_LEXA_EXTENDIDA.len() - k + 1))
    } else {
        let menor = menor_distancia(CAJA_LEXA_EXTENDIDA, consenso)
            .min(menor_distancia(CAJA_LEXA_EXTENDIDA, &inverso));
        let largo = CAJA_LEXA_EXTENDIDA.len();
        (menor, largo, 2 * (k - largo + 1))
    };
    if comparaciones as f64 * binomial_hasta(comparadas, menor, 0.75) <= AZAR_TOLERADO {
        Some((menor, comparadas))
    } else {
        None
    }
}

fn pares_no_complementarios(consenso: &[u8]) -> usize {
    let inverso = secuencia::complemento_inverso(consenso);
    motivos::distancia_hamming(consenso, &inverso) / 2
}

struct Comparado {
    metodo: &'static str,
    consenso: Vec<u8>,
    puntos: usize,
    bits: f64,
    milisegundos: f64,
}

fn medir_motivos(motivos_elegidos: &[&[u8]]) -> (Vec<u8>, usize, f64) {
    let entropia = motivos::entropia_por_columna(motivos_elegidos);
    let bits: f64 = entropia.iter().map(|valor| (2.0 - valor).max(0.0)).sum();
    (
        motivos::consenso(motivos_elegidos),
        motivos::puntuacion(motivos_elegidos),
        bits,
    )
}

fn costo_de_gibbs(textos: &[Vec<u8>], k: usize) -> f64 {
    let t = textos.len() as f64;
    let mas_larga = textos.iter().map(Vec::len).max().unwrap_or(0);
    GIBBS_REINICIOS as f64
        * (GIBBS_ITERACIONES as f64 + 1.0)
        * ((mas_larga + 1 - k) as f64 + 3.0 * t)
        * k as f64
}

fn costo_codicioso(textos: &[Vec<u8>], k: usize) -> f64 {
    let resto: usize = textos[1..].iter().map(|texto| texto.len() + 1 - k).sum();
    let t = textos.len();
    (textos[0].len() + 1 - k) as f64 * (resto as f64 * k as f64 + (t * t * k) as f64)
}

fn comparar_con_clasicos(textos: &[Vec<u8>], k: usize) -> Result<Vec<Comparado>, String> {
    let mut filas = Vec::new();
    if costo_de_gibbs(textos, k) <= PRESUPUESTO_CLASICOS {
        let reloj = Instant::now();
        let hallado =
            motivos::muestreo_gibbs(textos, k, GIBBS_ITERACIONES, GIBBS_REINICIOS, GIBBS_SEMILLA)?;
        let milisegundos = reloj.elapsed().as_secs_f64() * 1000.0;
        let elegidos = motivos::motivos_en(textos, &hallado.posiciones, k);
        let (consenso, puntos, bits) = medir_motivos(&elegidos);
        filas.push(Comparado {
            metodo: "Muestreo de Gibbs",
            consenso,
            puntos,
            bits,
            milisegundos,
        });
    }
    if costo_codicioso(textos, k) <= PRESUPUESTO_CLASICOS {
        let reloj = Instant::now();
        let hallado = motivos::busqueda_codiciosa(textos, k, 1.0)?;
        let milisegundos = reloj.elapsed().as_secs_f64() * 1000.0;
        let elegidos = motivos::motivos_en(textos, &hallado.posiciones, k);
        let (consenso, puntos, bits) = medir_motivos(&elegidos);
        filas.push(Comparado {
            metodo: "Búsqueda codiciosa con pseudoconteos",
            consenso,
            puntos,
            bits,
            milisegundos,
        });
    }
    Ok(filas)
}

fn es_mejor(a: &Comparado, b: &Comparado) -> bool {
    a.puntos < b.puntos || (a.puntos == b.puntos && a.bits > b.bits + 1e-9)
}

fn perfil_en_filas(perfil: &[[f64; 4]]) -> Vec<Vec<f64>> {
    let mut filas: Vec<Vec<f64>> = (0..4).map(|_| Vec::with_capacity(perfil.len())).collect();
    for columna in perfil {
        for (base, fila) in filas.iter_mut().enumerate() {
            fila.push((columna[base] * 1000.0).round() / 1000.0);
        }
    }
    filas
}

fn agregar_bloques(salida: &mut Salida, pedido: &Pedido, motivo: &MotivoEm, milisegundos: f64) {
    let conjunto = &pedido.conjunto;
    let opciones = &pedido.opciones;
    let k = opciones.k;
    let consenso = if opciones.ambas_hebras {
        format!(
            "{} (en la otra hebra: {})",
            texto(&motivo.consenso),
            texto(&secuencia::complemento_inverso(&motivo.consenso))
        )
    } else {
        texto(&motivo.consenso)
    };
    salida.agregar(Bloque::valor("Consenso", consenso));
    let t = conjunto.textos.len();
    salida.agregar(Bloque::valor(
        "Log-verosimilitud frente al fondo",
        format!(
            "{} bits ({} por secuencia)",
            decimal(motivo.verosimilitud, 2),
            decimal(motivo.verosimilitud / t as f64, 2)
        ),
    ));
    let categorias: Vec<String> = motivo
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
        motivo.informacion.clone(),
    ));
    let mut filas = Vec::with_capacity(t);
    for (indice, sitio) in motivo.sitios.iter().enumerate() {
        let mut fila = vec![
            conjunto.nombres[indice].clone(),
            miles(sitio.posicion),
            coordenada(conjunto.inicios[indice], sitio.posicion),
            (sitio.hebra as char).to_string(),
            texto(&sitio_orientado(&conjunto.textos, indice, sitio, k)),
            decimal(sitio.probabilidad, 3),
        ];
        if opciones.modelo == Modelo::Zoops {
            fila.push(decimal(sitio.probabilidad_de_sitio, 3));
        }
        filas.push(fila);
    }
    let mut columnas = vec![
        "Secuencia",
        "Posición (contando desde 0)",
        "En E. coli (contando desde 1)",
        "Hebra",
        "Sitio (leído en su hebra)",
        "Probabilidad posterior del sitio",
    ];
    if opciones.modelo == Modelo::Zoops {
        columnas.push("Probabilidad de tener algún sitio");
    }
    salida.agregar(Bloque::tabla(
        "Mejor sitio de cada secuencia",
        &columnas,
        filas,
    ));
    salida.agregar(Bloque::calor(
        "Matriz de pesos (probabilidad de cada base por columna)",
        LETRAS.iter().map(|letra| letra.to_string()).collect(),
        (1..=k).map(|columna| columna.to_string()).collect(),
        perfil_en_filas(&motivo.perfil),
    ));
    let puntos: Vec<[f64; 2]> = motivo
        .historial
        .iter()
        .enumerate()
        .map(|(iteracion, &bits)| [iteracion as f64, bits])
        .collect();
    salida.agregar(Bloque::linea(
        "Log-verosimilitud en cada iteración de EM",
        "Iteración (0 = perfil de arranque)",
        "Bits frente al fondo",
        &puntos,
        Vec::new(),
    ));
    let mut detalles = format!(
        "{} k-mers probados como arranque con {} iteraciones cada uno; el mejor fue {}. Desde él, EM hizo {} iteraciones",
        miles(motivo.arranques_probados),
        motivos_em::ITERACIONES_DE_PRUEBA,
        texto(&motivo.arranque),
        miles(motivo.iteraciones)
    );
    if opciones.modelo == Modelo::Zoops {
        detalles.push_str(&format!(
            "; γ (fracción de secuencias con sitio) = {}",
            decimal(motivo.gamma, 3)
        ));
    }
    detalles.push_str(&format!(
        ". Fondo: A {}, C {}, G {}, T {}. Tiempo: {} ms.",
        decimal(motivo.fondo[0], 3),
        decimal(motivo.fondo[1], 3),
        decimal(motivo.fondo[2], 3),
        decimal(motivo.fondo[3], 3),
        decimal(milisegundos, 1)
    ));
    salida.agregar(Bloque::texto("Cómo llegó", detalles));
}

fn explicar_sitios(salida: &mut Salida, pedido: &Pedido, motivo: &MotivoEm) {
    let t = motivo.sitios.len();
    let seguros = motivo
        .sitios
        .iter()
        .filter(|sitio| sitio.probabilidad > UMBRAL_POSTERIOR)
        .count();
    let mut frase = format!(
        "En {seguros} de {t} secuencias el mejor sitio tiene probabilidad posterior mayor que {}: ahí el modelo concentra su confianza en una sola ventana, aunque eso no prueba que el motivo sea real: mira también los bits de información.",
        decimal(UMBRAL_POSTERIOR, 1)
    );
    if pedido.opciones.ambas_hebras {
        frase.push_str(" La probabilidad suma las dos hebras de la misma posición, porque en un sitio palindrómico el modelo la reparte entre ambas.");
    }
    if pedido.opciones.modelo == Modelo::Zoops {
        let sin_sitio: Vec<String> = motivo
            .sitios
            .iter()
            .enumerate()
            .filter(|(_, sitio)| sitio.probabilidad_de_sitio < UMBRAL_POSTERIOR)
            .map(|(indice, _)| pedido.conjunto.nombres[indice].clone())
            .collect();
        if sin_sitio.is_empty() {
            frase.push_str(&format!(
                " ZOOPS estima que el {} de las secuencias lleva el motivo y no descarta ninguna.",
                porcentaje(motivo.gamma, 0)
            ));
        } else {
            frase.push_str(&format!(
                " ZOOPS estima que el {} de las secuencias lleva el motivo y deja sin sitio (probabilidad de tenerlo menor que {}) a: {}.",
                porcentaje(motivo.gamma, 0),
                decimal(UMBRAL_POSTERIOR, 1),
                sin_sitio.join(", ")
            ));
        }
    } else if seguros < t {
        frase.push_str(" OOPS obliga a poner un sitio en cada secuencia; si sospechas que alguna no tiene el motivo, prueba el modelo ZOOPS.");
    }
    salida.explicar(frase);
}

fn explicar_informacion(salida: &mut Salida, motivo: &MotivoEm) {
    let k = motivo.informacion.len();
    let total: f64 = motivo.informacion.iter().sum();
    let fuertes = motivo
        .informacion
        .iter()
        .filter(|&&bits| bits >= 1.0)
        .count();
    let mut frase = format!(
        "El motivo reúne {} bits de información de un máximo de {} (2 por columna); {fuertes} de sus {k} columnas pasan de 1 bit.",
        decimal(total, 1),
        2 * k
    );
    if total / (k as f64) < 1.0 {
        frase.push_str(" Con menos de 1 bit por columna en promedio el motivo es dudoso: en secuencias al azar EM también encuentra patrones de entre 0.5 y 1 bit por columna, así que puede ser ruido o k puede ser demasiado largo.");
    } else if fuertes < k {
        let debiles: Vec<usize> = motivo
            .informacion
            .iter()
            .enumerate()
            .filter(|(_, bits)| **bits < 1.0)
            .map(|(columna, _)| columna + 1)
            .collect();
        let mut en_los_bordes = true;
        for (orden, &columna) in debiles.iter().enumerate() {
            let desde_el_inicio = columna == orden + 1;
            let hasta_el_final = columna + (debiles.len() - orden) == k + 1;
            if !desde_el_inicio && !hasta_el_final {
                en_los_bordes = false;
            }
        }
        let lista: Vec<String> = debiles.iter().map(usize::to_string).collect();
        if en_los_bordes {
            frase.push_str(&format!(
                " Las columnas débiles ({}) están en los bordes: el motivo real puede ser algo más corto que k.",
                lista.join(", ")
            ));
        } else {
            frase.push_str(&format!(
                " Las columnas débiles ({}) están dentro del motivo: son posiciones que varían entre sitios, como el espaciador central de los sitios palindrómicos.",
                lista.join(", ")
            ));
        }
    }
    salida.explicar(frase);
}

fn explicar_verosimilitud(salida: &mut Salida, motivo: &MotivoEm, t: usize) {
    let inicial = motivo.historial.first().copied().unwrap_or(0.0);
    let mut frase = format!(
        "La log-verosimilitud frente al fondo subió de {} a {} bits en {} iteraciones de EM ({} bits por secuencia): cuanto más alta, mejor explica el motivo tus secuencias que el puro fondo.",
        decimal(inicial, 1),
        decimal(motivo.verosimilitud, 1),
        miles(motivo.iteraciones),
        decimal(motivo.verosimilitud / t as f64, 1)
    );
    if motivo.iteraciones >= motivos_em::ITERACIONES_MAXIMAS {
        frase.push_str(&format!(
            " EM llegó al tope de {} iteraciones sin estabilizarse del todo: tómalo como aproximado.",
            motivos_em::ITERACIONES_MAXIMAS
        ));
    } else {
        frase.push_str(
            " EM se detuvo porque la matriz ya no cambiaba (menos de 0.000001 por casilla).",
        );
    }
    salida.explicar(frase);
}

fn explicar_lexa(salida: &mut Salida, consenso: &[u8]) {
    let Some((diferencias, comparadas)) = parecido_a_lexa(consenso) else {
        return;
    };
    let mut frase = if diferencias == 0 {
        format!(
            "El consenso {} coincide con la caja SOS de E. coli (núcleo de 16 pb {CAJA_LEXA}), donde se une el represor LexA.",
            texto(consenso)
        )
    } else {
        format!(
            "El consenso {} se parece a la caja SOS de E. coli (núcleo de 16 pb {CAJA_LEXA}), donde se une el represor LexA: difiere en {diferencias} de {comparadas} letras, algo que al azar pasaría menos de 1 vez en 1 000.",
            texto(consenso)
        )
    };
    if consenso.len() >= 8 && pares_no_complementarios(consenso) <= consenso.len() / 8 {
        frase.push_str(" Como ella, es casi un palíndromo: LexA se une como dímero, una mitad sobre cada brazo CTGT…ACAG.");
    }
    salida.explicar(frase);
}

fn explicar_comparacion(
    salida: &mut Salida,
    propio: &Comparado,
    clasicos: &[Comparado],
    ambas: bool,
) {
    if clasicos.is_empty() {
        salida.explicar("La entrada es demasiado grande para correr también el muestreo de Gibbs y la búsqueda codiciosa en unos segundos, así que esta vez no se comparan.");
        return;
    }
    let mut mejor = propio;
    for otro in clasicos {
        if es_mejor(otro, mejor) {
            mejor = otro;
        }
    }
    let resumen: Vec<String> = clasicos
        .iter()
        .map(|otro| {
            format!(
                "«{}» {} (puntuación {})",
                otro.metodo,
                texto(&otro.consenso),
                miles(otro.puntos)
            )
        })
        .collect();
    let mut frase = format!(
        "Con la misma entrada, MEME da puntuación {} (letras que no coinciden con su consenso; menos es mejor) frente a {}.",
        miles(propio.puntos),
        resumen.join(" y ")
    );
    if std::ptr::eq(mejor, propio) {
        frase.push_str(" MEME encontró el motivo más conservado.");
    } else {
        frase.push_str(&format!(
            " Aquí «{}» encontró un motivo más conservado por esa medida; recuerda que EM maximiza la verosimilitud y no la puntuación.",
            mejor.metodo
        ));
    }
    if ambas {
        frase.push_str(" Los métodos del capítulo 2 solo miran la hebra directa; MEME pudo leer sitios en la otra.");
    }
    salida.explicar(frase);
}

fn meme(entrada: &Entrada) -> Result<Salida, String> {
    let pedido = leer_pedido(entrada)?;
    let k = pedido.opciones.k;
    let textos = &pedido.conjunto.textos;
    let reloj = Instant::now();
    let motivo = motivos_em::meme(textos, &pedido.opciones)?;
    let milisegundos = reloj.elapsed().as_secs_f64() * 1000.0;
    let mut salida = Salida::default();
    agregar_bloques(&mut salida, &pedido, &motivo, milisegundos);
    let orientados: Vec<Vec<u8>> = motivo
        .sitios
        .iter()
        .enumerate()
        .map(|(indice, sitio)| sitio_orientado(textos, indice, sitio, k))
        .collect();
    let vistas: Vec<&[u8]> = orientados.iter().map(Vec::as_slice).collect();
    let (consenso_de_sitios, puntos, bits) = medir_motivos(&vistas);
    let propio = Comparado {
        metodo: "MEME (EM)",
        consenso: consenso_de_sitios,
        puntos,
        bits,
        milisegundos,
    };
    let clasicos = comparar_con_clasicos(textos, k)?;
    let mut filas = Vec::new();
    for fila in std::iter::once(&propio).chain(clasicos.iter()) {
        filas.push(vec![
            fila.metodo.to_string(),
            texto(&fila.consenso),
            miles(fila.puntos),
            decimal(fila.bits, 1),
            decimal(fila.milisegundos, 1),
        ]);
    }
    salida.agregar(Bloque::tabla(
        "Comparación con los métodos del capítulo 2 (mismo k, mismas secuencias)",
        &[
            "Método",
            "Consenso de sus sitios",
            "Puntuación (menos es mejor)",
            "Información de sus sitios (bits)",
            "Tiempo (ms)",
        ],
        filas,
    ));
    let t = textos.len();
    explicar_sitios(&mut salida, &pedido, &motivo);
    explicar_informacion(&mut salida, &motivo);
    explicar_verosimilitud(&mut salida, &motivo, t);
    explicar_lexa(&mut salida, &motivo.consenso);
    explicar_comparacion(
        &mut salida,
        &propio,
        &clasicos,
        pedido.opciones.ambas_hebras,
    );
    if pedido.recortados {
        salida.explicar(format!(
            "Para no pasar de unos segundos se probaron solo los primeros {} k-mers distintos como arranque; con secuencias más cortas o menos secuencias se prueban más.",
            miles(pedido.opciones.arranques)
        ));
    }
    Ok(salida)
}

fn tabla_de_perfil(perfil: &[[f64; 4]], cifras: usize) -> (Vec<Vec<String>>, Vec<(usize, usize)>) {
    let mut filas = Vec::with_capacity(4);
    for (base, letra) in LETRAS.iter().enumerate() {
        let mut fila = vec![letra.to_string()];
        for columna in perfil {
            fila.push(decimal(columna[base], cifras));
        }
        filas.push(fila);
    }
    let mut resaltar = Vec::with_capacity(perfil.len());
    for (columna, probabilidades) in perfil.iter().enumerate() {
        let mut mayor = 0;
        for base in 1..4 {
            if probabilidades[base] > probabilidades[mayor] {
                mayor = base;
            }
        }
        resaltar.push((mayor, columna + 1));
    }
    (filas, resaltar)
}

fn columnas_de_perfil(k: usize) -> Vec<String> {
    let mut columnas = vec!["base".to_string()];
    for columna in 1..=k {
        columnas.push(columna.to_string());
    }
    columnas
}

fn con_perfil(pizarra: Pizarra, titulo: &str, perfil: &[[f64; 4]], cifras: usize) -> Pizarra {
    let (filas, resaltar) = tabla_de_perfil(perfil, cifras);
    let columnas = columnas_de_perfil(perfil.len());
    let referencias: Vec<&str> = columnas.iter().map(String::as_str).collect();
    pizarra.tabla(titulo, &referencias, filas, &resaltar)
}

fn nombre_de_candidato(candidato: usize, ambas: bool) -> (usize, char) {
    if ambas {
        (
            candidato / 2,
            if candidato.is_multiple_of(2) {
                '+'
            } else {
                '-'
            },
        )
    } else {
        (candidato, '+')
    }
}

fn lista_de_ventanas(z: &[f64], ambas: bool) -> (Vec<String>, Vec<usize>) {
    let mut elementos = Vec::with_capacity(z.len());
    let mut mayor = 0;
    for (candidato, &probabilidad) in z.iter().enumerate() {
        let (posicion, hebra) = nombre_de_candidato(candidato, ambas);
        elementos.push(format!(
            "{posicion} ({hebra}): {}",
            decimal(probabilidad, 4)
        ));
        if probabilidad > z[mayor] {
            mayor = candidato;
        }
    }
    (elementos, vec![mayor])
}

fn fondo_en_texto(fondo: &[f64; 4]) -> String {
    format!(
        "A {} · C {} · G {} · T {}",
        decimal(fondo[0], 3),
        decimal(fondo[1], 3),
        decimal(fondo[2], 3),
        decimal(fondo[3], 3)
    )
}

fn iteraciones_visibles(iteracion: usize) -> bool {
    iteracion <= ITERACIONES_EN_DETALLE || iteracion.is_multiple_of(CADA_CUANTAS_ITERACIONES)
}

struct Escena<'a> {
    textos: &'a [Vec<u8>],
    nombres: &'a [String],
    opciones: Opciones,
    fondo: [f64; 4],
}

impl Escena<'_> {
    fn pizarra(&self, iteraciones: usize, verosimilitud: f64, gamma: f64) -> Pizarra {
        let mut pizarra = Pizarra::nueva()
            .variable("k", self.opciones.k)
            .variable("modelo", nombre_del_modelo(self.opciones.modelo))
            .variable("ambas_hebras", self.opciones.ambas_hebras)
            .variable("fondo", fondo_en_texto(&self.fondo))
            .variable("iteraciones", iteraciones)
            .variable(
                "verosimilitud",
                format!("{} bits", decimal(verosimilitud, 2)),
            );
        if self.opciones.modelo == Modelo::Zoops {
            pizarra = pizarra.variable("gamma", decimal(gamma, 3));
        }
        pizarra
    }

    fn paso_e_de_una_secuencia(
        &self,
        narracion: &mut Narracion,
        indice: usize,
        z: &[f64],
        iteraciones: usize,
        verosimilitud: f64,
        parametros: &Parametros,
    ) {
        let k = self.opciones.k;
        let ambas = self.opciones.ambas_hebras;
        let sitio = motivos_em::mejor_sitio(z, ambas);
        let (elementos, resaltar) = lista_de_ventanas(z, ambas);
        let ventana = &self.textos[indice][sitio.posicion..sitio.posicion + k];
        let leida = if sitio.hebra == b'-' {
            secuencia::complemento_inverso(ventana)
        } else {
            ventana.to_vec()
        };
        let mut nota = format!(
            "Paso E en {}: cada ventana recibe su probabilidad de ser el sitio, proporcional a cuánto más probable es según la matriz que según el fondo. La mejor es la posición {} (hebra {}), {}, con probabilidad {}.",
            self.nombres[indice],
            sitio.posicion,
            sitio.hebra as char,
            texto(&leida),
            decimal(sitio.probabilidad, 3)
        );
        if self.opciones.modelo == Modelo::Zoops {
            nota.push_str(&format!(
                " Entre todas suman {}: es la probabilidad de que esta secuencia tenga sitio.",
                decimal(sitio.probabilidad_de_sitio, 3)
            ));
        } else {
            nota.push_str(" En OOPS todas suman 1: la secuencia tiene exactamente un sitio.");
        }
        narracion.agregar(
            self.pizarra(iteraciones, verosimilitud, parametros.gamma)
                .variable("i", indice)
                .secuencia(
                    &self.textos[indice],
                    &[(sitio.posicion, sitio.posicion + k, Marca::Actual)],
                )
                .lista(
                    "Probabilidad de cada ventana: posición (hebra)",
                    elementos,
                    &resaltar,
                )
                .paso("*valor *= escala;", nota),
        );
    }
}

fn pasos_meme(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let pedido = leer_pedido(entrada)?;
    let t = pedido.conjunto.textos.len();
    if t > SECUENCIAS_ANIMABLES {
        return Err(format!(
            "Para el paso a paso usa como mucho {SECUENCIAS_ANIMABLES} secuencias; aquí hay {t}."
        ));
    }
    let bases: usize = pedido.conjunto.textos.iter().map(Vec::len).sum();
    if bases > BASES_ANIMABLES {
        return Err(format!(
            "Para el paso a paso usa como mucho {} pb en total; aquí hay {}.",
            miles(BASES_ANIMABLES),
            miles(bases)
        ));
    }
    let opciones = pedido.opciones;
    let k = opciones.k;
    let ambas = opciones.ambas_hebras;
    let codigos = motivos_em::codificar(&pedido.conjunto.textos)?;
    let fondo = motivos_em::frecuencias_de_fondo(&codigos, ambas);
    let escena = Escena {
        textos: &pedido.conjunto.textos,
        nombres: &pedido.conjunto.nombres,
        opciones,
        fondo,
    };
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let hebras_en_texto = if ambas {
        "en las dos hebras"
    } else {
        "solo en la hebra directa"
    };
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", k)
            .variable("modelo", nombre_del_modelo(opciones.modelo))
            .variable("ambas_hebras", ambas)
            .secuencia(&escena.textos[0], &[])
            .paso(
                "let codigos = codificar(adn)?;",
                format!(
                    "Tienes {t} secuencias ({} pb en total) y buscas un motivo de {k} letras con el modelo {}, {hebras_en_texto}. Primero cada letra se vuelve un número: A = 0, C = 1, G = 2, T = 3.",
                    miles(bases),
                    nombre_del_modelo(opciones.modelo)
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("k", k)
            .variable("modelo", nombre_del_modelo(opciones.modelo))
            .variable("ambas_hebras", ambas)
            .variable("fondo", fondo_en_texto(&fondo))
            .paso(
                "conteos[codigo as usize] += 1.0;",
                if ambas {
                    "Se cuentan las bases de todas las secuencias (más 1 de cada una) para el modelo de fondo: lo que esperas sin motivo. Como buscas en las dos hebras, A y T comparten frecuencia, y C y G también."
                } else {
                    "Se cuentan las bases de todas las secuencias (más 1 de cada una) para el modelo de fondo: lo que esperas sin motivo."
                },
            ),
    );
    let semillas = motivos_em::arranques(&codigos, k, opciones.arranques);
    let primeras: Vec<String> = semillas
        .iter()
        .take(12)
        .map(|semilla| texto(&motivos_em::decodificar(semilla)))
        .collect();
    narracion.agregar(
        escena
            .pizarra(0, 0.0, 0.0)
            .variable("semillas", miles(semillas.len()))
            .lista("Primeros k-mers de arranque", primeras, &[])
            .paso(
                "if vistos.insert(ventana) {",
                format!(
                    "Como MEME, cada k-mer distinto de tus secuencias es un posible punto de arranque: hay {}. Cada uno se convertirá en una matriz inicial.",
                    miles(semillas.len())
                ),
            ),
    );
    let puntajes: Vec<f64> = semillas
        .par_iter()
        .map(|semilla| motivos_em::probar_arranque(&codigos, semilla, &fondo, &opciones))
        .collect();
    for (indice, semilla) in semillas.iter().take(SEMILLAS_EN_DETALLE).enumerate() {
        narracion.agregar(
            escena
                .pizarra(0, puntajes[indice], 0.0)
                .variable("semilla", texto(&motivos_em::decodificar(semilla)))
                .variable("puntaje", format!("{} bits", decimal(puntajes[indice], 2)))
                .paso(
                    "let recorrido = iterar(codigos, inicial, fondo, opciones, ITERACIONES_DE_PRUEBA);",
                    format!(
                        "Arranque {}: desde {} se hacen solo {} iteraciones de EM y se anota la log-verosimilitud alcanzada, {} bits. Es una prueba rápida para no gastar tiempo en arranques malos.",
                        indice + 1,
                        texto(&motivos_em::decodificar(semilla)),
                        motivos_em::ITERACIONES_DE_PRUEBA,
                        decimal(puntajes[indice], 2)
                    ),
                ),
        );
    }
    let mut mejor = 0;
    let mut mejoras = 0;
    for (indice, &puntaje) in puntajes.iter().enumerate() {
        if puntaje > puntajes[mejor] {
            mejor = indice;
            mejoras += 1;
            if mejoras <= MEJORAS_MOSTRADAS {
                narracion.agregar(
                    escena
                        .pizarra(0, puntaje, 0.0)
                        .variable("mejor", indice)
                        .variable("puntaje", format!("{} bits", decimal(puntaje, 2)))
                        .paso(
                            "if puntaje > puntajes[mejor] {",
                            format!(
                                "El arranque {} ({}) llega a {} bits: supera al mejor anterior y pasa a ser el elegido.",
                                indice + 1,
                                texto(&motivos_em::decodificar(semillas[indice])),
                                decimal(puntaje, 2)
                            ),
                        ),
                );
            }
        }
    }
    let mut orden: Vec<usize> = (0..puntajes.len()).collect();
    orden.sort_by(|a, b| puntajes[*b].total_cmp(&puntajes[*a]).then(a.cmp(b)));
    let podio: Vec<String> = orden
        .iter()
        .take(10)
        .map(|&indice| {
            format!(
                "{}: {} bits",
                texto(&motivos_em::decodificar(semillas[indice])),
                decimal(puntajes[indice], 2)
            )
        })
        .collect();
    narracion.agregar(
        escena
            .pizarra(0, puntajes[mejor], 0.0)
            .variable("mejor", mejor)
            .variable("puntaje", format!("{} bits", decimal(puntajes[mejor], 2)))
            .lista("Los 10 mejores arranques", podio, &[0])
            .paso(
                "let mejor = mejor_arranque(&puntajes);",
                format!(
                    "Tras probar los {} arranques ({} mejoras del récord), gana {} con {} bits. Solo este se llevará hasta el final.",
                    miles(semillas.len()),
                    miles(mejoras),
                    texto(&motivos_em::decodificar(semillas[mejor])),
                    decimal(puntajes[mejor], 2)
                ),
            ),
    );
    let mut parametros = motivos_em::parametros_de_arranque(semillas[mejor], opciones.modelo);
    narracion.agregar(con_perfil(
        escena.pizarra(0, 0.0, parametros.gamma),
        "Matriz de arranque",
        &parametros.perfil,
        2,
    )
    .paso(
        "columna[codigo as usize] = PROBABILIDAD_DE_ARRANQUE;",
        format!(
            "La matriz inicial pone probabilidad 0.5 a la letra del arranque en cada columna y reparte el resto (0.17) entre las otras tres. Es una apuesta suave: EM puede corregirla.{}",
            if opciones.modelo == Modelo::Zoops {
                " En ZOOPS además γ, la fracción de secuencias con sitio, arranca en 0.5."
            } else {
                ""
            }
        ),
    ));
    let mut posteriores = motivos_em::posteriores_vacios(&codigos, &opciones);
    let mut historial: Vec<f64> = Vec::new();
    let mut iteraciones = 0;
    while iteraciones < motivos_em::ITERACIONES_MAXIMAS {
        let verosimilitud =
            motivos_em::paso_e(&codigos, &parametros, &fondo, &opciones, &mut posteriores);
        historial.push(verosimilitud / std::f64::consts::LN_2);
        let bits = verosimilitud / std::f64::consts::LN_2;
        let nuevos = motivos_em::paso_m(&codigos, &posteriores, &fondo, &opciones);
        let cambio = motivos_em::cambio_entre(&parametros, &nuevos);
        let primera = iteraciones == 0;
        if primera {
            let pesos = motivos_em::pesos_logaritmicos(&parametros.perfil, &fondo);
            let primera_ventana = &codigos[0][..k];
            let mut puntaje = 0.0;
            for (columna, &codigo) in pesos.iter().zip(primera_ventana) {
                puntaje += columna[codigo as usize];
            }
            narracion.agregar(
                escena
                    .pizarra(iteraciones, 0.0, parametros.gamma)
                    .variable("i", 0)
                    .secuencia(&escena.textos[0], &[(0, k, Marca::Actual)])
                    .paso(
                        "puntajes.push(puntaje_directo(pesos, ventana));",
                        format!(
                            "Paso E, primera iteración. La primera ventana de {}, {}, se puntúa sumando, columna por columna, el logaritmo de (probabilidad en la matriz / probabilidad en el fondo): da {}. Positivo = se parece más al motivo que al fondo.{}",
                            escena.nombres[0],
                            texto(&escena.textos[0][..k]),
                            decimal(puntaje, 2),
                            if ambas {
                                " Lo mismo se hace con su complemento inverso, para la otra hebra."
                            } else {
                                ""
                            }
                        ),
                    ),
            );
        }
        if primera || iteraciones == 1 {
            for (indice, z) in posteriores.iter().enumerate() {
                if indice == 0 || primera {
                    escena.paso_e_de_una_secuencia(
                        &mut narracion,
                        indice,
                        z,
                        iteraciones,
                        bits,
                        &parametros,
                    );
                }
            }
        }
        if primera {
            let (conteos, sitios_esperados) =
                motivos_em::conteos_esperados(&codigos, &posteriores, &opciones);
            narracion.agregar(con_perfil(
                escena.pizarra(iteraciones, bits, parametros.gamma),
                "Conteos esperados (cada ventana pesa su probabilidad)",
                &conteos,
                2,
            )
            .paso(
                "sumar_directo(&mut conteos, ventana, peso);",
                format!(
                    "Paso M: en vez de quedarse con un solo k-mer por secuencia, cada ventana suma sus letras pesadas por su probabilidad. En total entran {} sitios esperados.",
                    decimal(sitios_esperados, 2)
                ),
            ));
            narracion.agregar(con_perfil(
                escena.pizarra(iteraciones, bits, nuevos.gamma),
                "Matriz nueva",
                &nuevos.perfil,
                2,
            )
            .paso(
                "probabilidades[base] = (columna[base] + beta * fondo[base]) / total;",
                format!(
                    "Cada columna se normaliza tras sumar un pseudoconteo de {} repartido como el fondo, para que ninguna probabilidad sea 0. El consenso pasa a ser {}.",
                    decimal(opciones.pseudoconteo, 1),
                    texto(&motivos_em::consenso_del_perfil(&nuevos.perfil))
                ),
            ));
            if opciones.modelo == Modelo::Zoops {
                narracion.agregar(
                    escena.pizarra(iteraciones, bits, nuevos.gamma).paso(
                        "Modelo::Zoops => (sitios_esperados / codigos.len() as f64).clamp(GAMMA_MINIMO, 1.0),",
                        format!(
                            "γ se actualiza con los sitios esperados entre el número de secuencias: {} / {t} = {}.",
                            decimal(sitios_esperados, 2),
                            decimal(nuevos.gamma, 3)
                        ),
                    ),
                );
            }
        }
        let iteracion = iteraciones + 1;
        if iteraciones_visibles(iteracion) {
            let historial_en_texto: Vec<String> = historial
                .iter()
                .enumerate()
                .map(|(numero, valor)| {
                    format!("iteración {}: {} bits", numero + 1, decimal(*valor, 3))
                })
                .collect();
            let ultimo = historial_en_texto.len() - 1;
            let subida = if iteraciones == 0 {
                "Es el punto de partida.".to_string()
            } else {
                let diferencia = historial[iteraciones] - historial[iteraciones - 1];
                if diferencia >= 0.0 {
                    format!(
                        "Subió {} bits respecto a la iteración anterior: cada vuelta de EM mejora el ajuste.",
                        cientifico(diferencia)
                    )
                } else {
                    format!(
                        "Bajó apenas {} bits: el pseudoconteo tira un poco de la matriz hacia el fondo.",
                        cientifico(-diferencia)
                    )
                }
            };
            narracion.agregar(
                escena
                    .pizarra(iteracion, bits, nuevos.gamma)
                    .variable("cambio", cientifico(cambio))
                    .variable("consenso", texto(&motivos_em::consenso_del_perfil(&nuevos.perfil)))
                    .lista("Log-verosimilitud por iteración", historial_en_texto, &[ultimo])
                    .paso(
                        "historial.push(verosimilitud / std::f64::consts::LN_2);",
                        format!(
                            "Iteración {iteracion}: con la matriz que entró, tus secuencias tienen {} bits de log-verosimilitud por encima del fondo. {subida} La matriz cambió como mucho {} en una casilla.",
                            decimal(bits, 3),
                            cientifico(cambio)
                        ),
                    ),
            );
        }
        parametros = nuevos;
        iteraciones += 1;
        if cambio < motivos_em::CAMBIO_MINIMO {
            narracion.agregar(
                escena
                    .pizarra(iteraciones, bits, parametros.gamma)
                    .variable("cambio", cientifico(cambio))
                    .paso(
                        "if cambio < CAMBIO_MINIMO {",
                        format!(
                            "El cambio ({}) ya es menor que 0.000001: EM convergió tras {iteraciones} iteraciones.",
                            cientifico(cambio)
                        ),
                    ),
            );
            break;
        }
    }
    let verosimilitud =
        motivos_em::paso_e(&codigos, &parametros, &fondo, &opciones, &mut posteriores);
    let bits = verosimilitud / std::f64::consts::LN_2;
    historial.push(bits);
    let mut filas = Vec::with_capacity(t);
    let mut resaltar = Vec::new();
    let mut sitios = Vec::with_capacity(t);
    for (indice, z) in posteriores.iter().enumerate() {
        let sitio = motivos_em::mejor_sitio(z, ambas);
        let ventana = &escena.textos[indice][sitio.posicion..sitio.posicion + k];
        let leida = if sitio.hebra == b'-' {
            secuencia::complemento_inverso(ventana)
        } else {
            ventana.to_vec()
        };
        if sitio.probabilidad > UMBRAL_POSTERIOR {
            resaltar.push((indice, 4));
        }
        filas.push(vec![
            escena.nombres[indice].clone(),
            sitio.posicion.to_string(),
            (sitio.hebra as char).to_string(),
            recortar(&leida, 30),
            decimal(sitio.probabilidad, 3),
            decimal(sitio.probabilidad_de_sitio, 3),
        ]);
        sitios.push(sitio);
    }
    let marcas: Vec<(usize, usize, Marca)> =
        vec![(sitios[0].posicion, sitios[0].posicion + k, Marca::Hallado)];
    narracion.agregar(
        escena
            .pizarra(iteraciones, bits, parametros.gamma)
            .secuencia(&escena.textos[0], &marcas)
            .tabla(
                "Mejor sitio de cada secuencia",
                &["secuencia", "posición", "hebra", "sitio", "probabilidad", "tiene sitio"],
                filas,
                &resaltar,
            )
            .paso(
                "sitios.push(mejor_sitio(z, opciones.ambas_hebras));",
                format!(
                    "Con la matriz final, un último paso E da la probabilidad de cada ventana; de cada secuencia se guarda la posición más probable (sumando sus dos hebras). {} de {t} superan 0.5 (resaltadas).",
                    resaltar.len()
                ),
            ),
    );
    let consenso = motivos_em::consenso_del_perfil(&parametros.perfil);
    let informacion: f64 = motivos_em::informacion_por_columna(&parametros.perfil)
        .iter()
        .sum();
    let ultimo = con_perfil(
        escena
            .pizarra(iteraciones, bits, parametros.gamma)
            .variable("consenso", texto(&consenso)),
        "Matriz final",
        &parametros.perfil,
        2,
    )
    .paso(
        "Ok(MotivoEm {",
        format!(
            "Fin: el motivo tiene consenso {}, {} bits de información en total y una log-verosimilitud de {} bits frente al fondo. En «Qué significa» se interpretan estos números y se comparan con Gibbs y la búsqueda codiciosa.",
            texto(&consenso),
            decimal(informacion, 1),
            decimal(bits, 2)
        ),
    );
    Ok(narracion.terminar(ultimo))
}

#[cfg(test)]
mod pruebas {
    use super::*;
    use serde_json::json;

    fn ejemplo() -> Entrada {
        Entrada::nueva(&[
            ("adn", REGIONES_SOS_Y_CONTROLES),
            ("k", "16"),
            ("modelo", "ZOOPS"),
            ("hebras", "ambas"),
        ])
    }

    #[test]
    fn el_paso_a_paso_termina_igual_que_nucleo() {
        let entrada = ejemplo();
        let pedido = leer_pedido(&entrada).unwrap();
        let motivo = motivos_em::meme(&pedido.conjunto.textos, &pedido.opciones).unwrap();
        let pasos = pasos_meme(&entrada).unwrap();
        assert!(pasos.len() <= MAXIMO_PASOS);
        let ultimo = pasos.last().unwrap();
        let consenso = texto(&motivo.consenso);
        assert!(ultimo.nota.contains(&consenso));
        assert!(ultimo.nota.contains(&decimal(motivo.verosimilitud, 2)));
        let variables = ultimo.estado["variables"].as_array().unwrap();
        assert!(variables.contains(&json!(["iteraciones", motivo.iteraciones.to_string()])));
        assert!(variables.contains(&json!(["consenso", consenso])));
    }

    #[test]
    fn el_ejemplo_encuentra_la_caja_sos_y_los_controles_sin_sitio() {
        let salida = crate::ejecutar(
            "motivos.meme",
            json!({ "adn": REGIONES_SOS_Y_CONTROLES, "k": "16", "modelo": "ZOOPS", "hebras": "ambas" }),
        )
        .unwrap();
        let todo = salida.interpretacion.join(" ");
        assert!(todo.contains("caja SOS"));
        assert!(todo.contains("deja sin sitio"));
        assert!(todo.contains("Muestreo de Gibbs"));
    }

    #[test]
    fn parecido_a_lexa_en_ambas_orientaciones() {
        assert_eq!(parecido_a_lexa(b"CTGTATATATATACAG"), Some((0, 16)));
        assert_eq!(parecido_a_lexa(b"CTGTATATAAATACAG"), Some((1, 16)));
        assert_eq!(parecido_a_lexa(b"GGGCCCGGGCCCGGGC"), None);
        assert_eq!(parecido_a_lexa(b"ACGT"), None);
        assert_eq!(parecido_a_lexa(b"TACTGCAT"), None);
        assert_eq!(parecido_a_lexa(b"TACTGTAT"), Some((0, 8)));
    }

    fn entrada_pequena(modelo: &str, hebras: &str, k: &str) -> Entrada {
        Entrada::nueva(&[
            (
                "adn",
                "TTGACATGCAGTCAGGATCCATTAGCAT\nGGATCCATAGCGCGATATTCGATCAG\nCGCGATTGACATAGCAGGATCCATTA\nAGCTTGGATCCATGCGCAATCGCTAG",
            ),
            ("k", k),
            ("modelo", modelo),
            ("hebras", hebras),
        ])
    }

    #[test]
    fn el_paso_a_paso_termina_igual_que_nucleo_en_todos_los_modos() {
        for modelo in ["OOPS", "ZOOPS"] {
            for hebras in ["una", "ambas"] {
                for k in ["3", "6", "9"] {
                    let entrada = entrada_pequena(modelo, hebras, k);
                    let pedido = leer_pedido(&entrada).unwrap();
                    let motivo =
                        motivos_em::meme(&pedido.conjunto.textos, &pedido.opciones).unwrap();
                    let pasos = pasos_meme(&entrada).unwrap();
                    assert!(pasos.len() <= MAXIMO_PASOS);
                    let ultimo = pasos.last().unwrap();
                    let consenso = texto(&motivo.consenso);
                    assert!(ultimo.nota.contains(&consenso));
                    assert!(ultimo.nota.contains(&decimal(motivo.verosimilitud, 2)));
                    let variables = ultimo.estado["variables"].as_array().unwrap();
                    assert!(
                        variables.contains(&json!(["iteraciones", motivo.iteraciones.to_string()]))
                    );
                    assert!(variables.contains(&json!(["consenso", consenso])));
                }
            }
        }
    }

    #[test]
    fn el_azar_no_se_presenta_como_motivo_firme() {
        let salida = crate::ejecutar(
            "motivos.meme",
            json!({
                "adn": "ACGTTGCAAGCTTCGATCGGATTACAGCTAGGCTAATCGGCTA\nTTAGCCGATATCGCGGCATTACGGATCTAGCTTAGCGAT\nGGCATCGATTACGCGATAGCTTACGATCGGCTAATCGAT\nCGATTAGCGGCTATACGCTAGCATCGGATTCGAGCTAGC",
                "k": "8",
                "modelo": "OOPS",
                "hebras": "una"
            }),
        )
        .unwrap();
        let todo = salida.interpretacion.join(" ");
        assert!(todo.contains("dudoso"));
        assert!(!todo.contains("está seguro"));
    }

    #[test]
    fn rechaza_entradas_sin_sentido() {
        let error = crate::ejecutar(
            "motivos.meme",
            json!({ "adn": "ACGTACGT ACGTAC", "k": "1", "modelo": "OOPS", "hebras": "una" }),
        )
        .unwrap_err();
        assert!(error.contains("al menos 2"));
        let error = crate::ejecutar(
            "motivos.meme",
            json!({ "adn": "ACGTACGT", "k": "4", "modelo": "OOPS", "hebras": "una" }),
        )
        .unwrap_err();
        assert!(error.contains("2 secuencias"));
    }
}
