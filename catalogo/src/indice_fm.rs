use crate::formato::{cientifico, decimal, miles};
use crate::pizarra::{Marca, Narracion, Pizarra};
use crate::{Algoritmo, Bloque, Campo, Entrada, Paso, Salida, Tipo};
use nucleo::ecoli;
use nucleo::indice_fm::{
    self, Acierto, Bloque as BloqueDelIndice, CENTINELA, IndiceDoble, IndiceFm, Intervalo,
    LETRAS_DE_UN_BYTE, LETRAS_DEL_INDICE, PASO_DE_CONTROL, PASO_DE_MUESTRA, Rama, SIMBOLOS, VACIO,
};
use nucleo::replicacion;
use nucleo::secuencia::complemento_inverso;
use std::borrow::Cow;
use std::sync::OnceLock;
use std::time::{Duration, Instant};

const DIFERENCIAS: &[&str] = &["0", "1", "2"];
const MAXIMO_PATRONES: usize = 1000;
const MAXIMO_LARGO_PATRON: usize = 1000;
const MAXIMO_TEXTO_INDICE: usize = 12_000_000;
const MAXIMO_TEXTO_BWT: usize = 5_000_000;
const MAXIMO_ROTACIONES: usize = 40;
const MAXIMO_COMPARADOS: usize = 10;
const MAXIMO_FILAS: usize = 2000;
const MAXIMO_MOSTRADO: usize = 2000;
const TRAMOS: usize = 20;
const SIGNIFICANCIA: f64 = 0.01;
const TIEMPO_DE_MEDIDA: Duration = Duration::from_millis(2);
const TIEMPO_TOTAL_DE_MEDIDA: Duration = Duration::from_millis(1000);
const TIEMPO_DEL_ESCANEO: Duration = Duration::from_millis(2500);
const MAXIMO_REPETICIONES: usize = 200;
const MAXIMO_PASOS: usize = 400;
const LETRAS_DETALLADAS: usize = 45;
const RAMAS_DETALLADAS: usize = 25;
const PATRONES_DETALLADOS: usize = 30;
const FILAS_DETALLADAS: usize = 20;
const LARGO_DEL_SUFIJO: usize = 16;
const PRIMERAS_FILAS: usize = 8;
const CHI: &[u8] = b"GCTGGTGG";
const LARGO_DE_REPETICION: usize = 16;
const FACTOR_DE_AZAR: f64 = 1.5;

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![
        Algoritmo {
            id: "replicacion.indice_fm",
            capitulo: "replicacion",
            nombre: "Índice FM: buscar sin recorrer el texto",
            resumen: "Ordena el texto una sola vez (arreglo de sufijos con SA-IS, transformada de Burrows-Wheeler y puntos de control) y después encuentra cada patrón, en ambas hebras y con hasta 2 diferencias, en microsegundos: lo nuevo frente al escaneo clásico es que el tiempo de búsqueda ya no depende del largo del texto.",
            campos: &[
                Campo {
                    id: "patrones",
                    etiqueta: "Patrones (separados por espacios)",
                    tipo: Tipo::Secuencias,
                    ejemplo: "GTGCCAGCAGCCGCGGTAAT GCTGGTGG TTATCCACA",
                },
                Campo {
                    id: "texto",
                    etiqueta: "Texto de ADN",
                    tipo: Tipo::Secuencia,
                    ejemplo: "@ecoli",
                },
                Campo {
                    id: "d",
                    etiqueta: "d (diferencias permitidas)",
                    tipo: Tipo::Opcion(DIFERENCIAS),
                    ejemplo: "1",
                },
            ],
            funcion: "indice_fm::buscar_patrones",
            ejecutar: buscar_con_indice,
            pasos: Some(pasos_indice_fm),
        },
        Algoritmo {
            id: "replicacion.bwt",
            capitulo: "replicacion",
            nombre: "Transformada de Burrows-Wheeler",
            resumen: "La última columna de las rotaciones ordenadas del texto, calculada en tiempo lineal con un arreglo de sufijos (SA-IS) y deshecha letra por letra con el mapeo LF: lo nuevo es que reordena el texto sin perder nada y junta letras iguales, la base de los compresores por bloques y de los índices de los mapeadores de lecturas.",
            campos: &[Campo {
                id: "texto",
                etiqueta: "Texto",
                tipo: Tipo::Secuencia,
                ejemplo: "@ecoli:1198457-1198496",
            }],
            funcion: "indice_fm::transformada_e_inversa",
            ejecutar: transformada,
            pasos: Some(pasos_bwt),
        },
    ]
}

pub struct IndiceListo {
    indice: IndiceDoble,
    construccion: Duration,
}

pub fn preparar() -> &'static IndiceListo {
    static INDICE: OnceLock<IndiceListo> = OnceLock::new();
    INDICE.get_or_init(|| {
        let inicio = Instant::now();
        let indice =
            IndiceDoble::construir(ecoli::genoma()).expect("el genoma incluido no está vacío");
        IndiceListo {
            indice,
            construccion: inicio.elapsed(),
        }
    })
}

struct Pedido {
    patrones: Vec<Vec<u8>>,
    texto: Cow<'static, [u8]>,
    d: usize,
    completo: bool,
    inicio_ecoli: Option<usize>,
}

enum IndiceUsado {
    DeEcoli(&'static IndiceListo),
    Propio(Box<IndiceDoble>, Duration),
}

impl IndiceUsado {
    fn doble(&self) -> &IndiceDoble {
        match self {
            IndiceUsado::DeEcoli(listo) => &listo.indice,
            IndiceUsado::Propio(indice, _) => indice,
        }
    }

    fn directo(&self) -> &IndiceFm {
        &self.doble().directo
    }

    fn construccion(&self) -> Duration {
        match self {
            IndiceUsado::DeEcoli(listo) => listo.construccion,
            IndiceUsado::Propio(_, tiempo) => *tiempo,
        }
    }

    fn ya_estaba(&self) -> bool {
        matches!(self, IndiceUsado::DeEcoli(_))
    }
}

fn leer_pedido(entrada: &Entrada) -> Result<Pedido, String> {
    let d = entrada
        .opcion("d", DIFERENCIAS)?
        .parse::<usize>()
        .map_err(|_| "«d» debe ser 0, 1 o 2.".to_string())?;
    let patrones = entrada.secuencias("patrones")?;
    if patrones.len() > MAXIMO_PATRONES {
        return Err(format!(
            "Escribiste {} patrones; aquí se admiten hasta {}.",
            miles(patrones.len()),
            miles(MAXIMO_PATRONES)
        ));
    }
    if let Some(largo) = patrones
        .iter()
        .map(Vec::len)
        .find(|&largo| largo > MAXIMO_LARGO_PATRON)
    {
        return Err(format!(
            "Un patrón tiene {} letras; aquí se admiten hasta {}.",
            miles(largo),
            miles(MAXIMO_LARGO_PATRON)
        ));
    }
    indice_fm::revisar_patrones(&patrones, d)?;
    let completo = entrada
        .texto("texto")
        .map(|texto| texto.trim() == "@ecoli")
        .unwrap_or(false);
    let texto: Cow<'static, [u8]> = if completo {
        Cow::Borrowed(ecoli::genoma())
    } else {
        Cow::Owned(entrada.secuencia("texto")?)
    };
    if texto.len() > MAXIMO_TEXTO_INDICE {
        return Err(format!(
            "El texto tiene {} letras; aquí se indexan hasta {}.",
            miles(texto.len()),
            miles(MAXIMO_TEXTO_INDICE)
        ));
    }
    Ok(Pedido {
        patrones,
        texto,
        d,
        completo,
        inicio_ecoli: inicio_en_ecoli(entrada, "texto"),
    })
}

fn preparar_indice(pedido: &Pedido) -> Result<IndiceUsado, String> {
    if pedido.completo {
        return Ok(IndiceUsado::DeEcoli(preparar()));
    }
    let inicio = Instant::now();
    let indice = IndiceDoble::construir(&pedido.texto)?;
    Ok(IndiceUsado::Propio(Box::new(indice), inicio.elapsed()))
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

struct Tiempos {
    contar: Duration,
    ubicar: Duration,
}

struct Comparacion {
    patrones: usize,
    indice: Duration,
    clasico: Duration,
    clasico_de_cada: Vec<Duration>,
    coinciden: bool,
}

fn tiempo_de_medida(patrones: usize) -> Duration {
    (TIEMPO_TOTAL_DE_MEDIDA / (2 * patrones.max(1) as u32)).min(TIEMPO_DE_MEDIDA)
}

fn medir(
    limite: Duration,
    mut trabajo: impl FnMut() -> Result<(), String>,
) -> Result<Duration, String> {
    let inicio = Instant::now();
    let mut veces: u32 = 0;
    while (veces as usize) < MAXIMO_REPETICIONES && (veces == 0 || inicio.elapsed() < limite) {
        trabajo()?;
        veces += 1;
    }
    Ok(inicio.elapsed() / veces)
}

fn medir_indice(
    indice: &IndiceDoble,
    patron: &[u8],
    d: usize,
    limite: Duration,
) -> Result<Tiempos, String> {
    let solo_este = [patron.to_vec()];
    let contar = medir(limite, || {
        std::hint::black_box(indice_fm::buscar_en_ambas_hebras(indice, patron, d));
        Ok(())
    })?;
    let ubicar = medir(limite, || {
        std::hint::black_box(indice_fm::buscar_patrones(indice, &solo_este, d)?);
        Ok(())
    })?;
    Ok(Tiempos { contar, ubicar })
}

fn escaneo_clasico(texto: &[u8], patron: &[u8], d: usize) -> (Vec<usize>, Vec<usize>) {
    let buscar = |buscado: &[u8]| {
        if d == 0 {
            replicacion::posiciones_del_patron(texto, buscado)
        } else {
            replicacion::posiciones_aproximadas(texto, buscado, d)
        }
    };
    let inverso = complemento_inverso(patron);
    let directas = buscar(patron);
    let inversas = if inverso == patron {
        Vec::new()
    } else {
        buscar(&inverso)
    };
    (directas, inversas)
}

fn posiciones_de_hebra(aciertos: &[Acierto], hebra: u8) -> Vec<usize> {
    aciertos
        .iter()
        .filter(|acierto| acierto.hebra == hebra)
        .map(|acierto| acierto.posicion)
        .collect()
}

fn comparar_con_el_escaneo(
    texto: &[u8],
    patrones: &[Vec<u8>],
    resultados: &[Vec<Acierto>],
    tiempos: &[Tiempos],
    d: usize,
    presupuesto: Duration,
) -> Comparacion {
    let mut comparacion = Comparacion {
        patrones: 0,
        indice: Duration::ZERO,
        clasico: Duration::ZERO,
        clasico_de_cada: Vec::new(),
        coinciden: true,
    };
    for ((patron, aciertos), tiempo) in patrones
        .iter()
        .zip(resultados)
        .zip(tiempos)
        .take(MAXIMO_COMPARADOS)
    {
        if comparacion.patrones > 0 && comparacion.clasico >= presupuesto {
            break;
        }
        let inicio = Instant::now();
        let (directas, inversas) = escaneo_clasico(texto, patron, d);
        let clasico = inicio.elapsed();
        comparacion.clasico += clasico;
        comparacion.clasico_de_cada.push(clasico);
        comparacion.indice += tiempo.ubicar;
        comparacion.patrones += 1;
        if posiciones_de_hebra(aciertos, b'+') != directas
            || posiciones_de_hebra(aciertos, b'-') != inversas
        {
            comparacion.coinciden = false;
        }
    }
    comparacion
}

fn patrones_comparados(cantidad: usize) -> String {
    if cantidad == 1 {
        "en el único patrón comparado".to_string()
    } else {
        format!("en los {} patrones comparados", miles(cantidad))
    }
}

fn veces_con_decimales(veces: f64) -> String {
    if veces >= 10.0 {
        format!("{} veces", miles(veces.round() as usize))
    } else {
        format!("{} veces", decimal(veces, 1))
    }
}

fn cantidad_esperada(valor: f64) -> String {
    if valor >= 100.0 {
        miles(valor.round() as usize)
    } else if valor >= 0.01 {
        decimal(valor, 2)
    } else {
        cientifico(valor)
    }
}

fn promedio(tiempos: impl Iterator<Item = Duration>) -> Duration {
    let mut suma = Duration::ZERO;
    let mut cantidad: u32 = 0;
    for tiempo in tiempos {
        suma += tiempo;
        cantidad += 1;
    }
    suma / cantidad.max(1)
}

fn duracion(tiempo: Duration) -> String {
    let segundos = tiempo.as_secs_f64();
    if segundos < 1e-3 {
        format!("{}\u{a0}µs", decimal(segundos * 1e6, 1))
    } else if segundos < 1.0 {
        format!("{}\u{a0}ms", decimal(segundos * 1e3, 1))
    } else {
        format!("{}\u{a0}s", decimal(segundos, 2))
    }
}

fn veces_mas_rapido(lento: Duration, rapido: Duration) -> f64 {
    lento.as_secs_f64() / rapido.as_secs_f64().max(1e-9)
}

fn proporciones(texto: &[u8]) -> [f64; 4] {
    let mut cuentas = [0usize; 4];
    for &letra in texto {
        match letra {
            b'A' => cuentas[0] += 1,
            b'C' => cuentas[1] += 1,
            b'G' => cuentas[2] += 1,
            b'T' => cuentas[3] += 1,
            _ => {}
        }
    }
    let total = texto.len().max(1) as f64;
    [
        cuentas[0] as f64 / total,
        cuentas[1] as f64 / total,
        cuentas[2] as f64 / total,
        cuentas[3] as f64 / total,
    ]
}

fn probabilidad_de_letra(letra: u8, proporcion: &[f64; 4]) -> f64 {
    match letra {
        b'A' => proporcion[0],
        b'C' => proporcion[1],
        b'G' => proporcion[2],
        _ => proporcion[3],
    }
}

fn probabilidad_con_diferencias(patron: &[u8], d: usize, proporcion: &[f64; 4]) -> f64 {
    let mut con_errores = vec![0.0; d + 1];
    con_errores[0] = 1.0;
    for &letra in patron {
        let acierta = probabilidad_de_letra(letra, proporcion);
        for errores in (0..=d).rev() {
            let desde_un_error = if errores > 0 {
                con_errores[errores - 1] * (1.0 - acierta)
            } else {
                0.0
            };
            con_errores[errores] = con_errores[errores] * acierta + desde_un_error;
        }
    }
    con_errores.iter().sum()
}

fn esperadas_al_azar(texto: &[u8], patron: &[u8], d: usize, proporcion: &[f64; 4]) -> f64 {
    if patron.len() > texto.len() {
        return 0.0;
    }
    let ventanas = (texto.len() - patron.len() + 1) as f64;
    let hebras = if complemento_inverso(patron) == patron {
        1.0
    } else {
        2.0
    };
    ventanas * hebras * probabilidad_con_diferencias(patron, d, proporcion)
}

fn logaritmo_factorial(n: usize) -> f64 {
    if n < 30 {
        return (2..=n).map(|k| (k as f64).ln()).sum();
    }
    let x = n as f64;
    x * x.ln() - x + 0.5 * (2.0 * std::f64::consts::PI * x).ln() + 1.0 / (12.0 * x)
}

fn poisson_al_menos(media: f64, veces: usize) -> f64 {
    if veces == 0 {
        return 1.0;
    }
    if media <= 0.0 {
        return 0.0;
    }
    let mut termino = (-media + veces as f64 * media.ln() - logaritmo_factorial(veces)).exp();
    let mut suma = 0.0;
    let mut k = veces;
    while termino > 0.0 && termino >= suma * 1e-15 && k < veces + 1_000_000 {
        suma += termino;
        k += 1;
        termino *= media / k as f64;
    }
    suma.min(1.0)
}

fn letra_del_simbolo(simbolo: usize) -> char {
    LETRAS_DEL_INDICE[simbolo] as char
}

fn hebra_texto(hebra: u8) -> &'static str {
    if hebra == b'+' { "+" } else { "−" }
}

fn texto_hallado(texto: &[u8], patron: &[u8], acierto: &Acierto) -> String {
    let fin = (acierto.posicion + patron.len()).min(texto.len());
    let tramo = &texto[acierto.posicion..fin];
    let leido = if acierto.hebra == b'-' {
        complemento_inverso(tramo)
    } else {
        tramo.to_vec()
    };
    leido
        .iter()
        .zip(patron)
        .map(|(&letra, &esperada)| {
            if letra == esperada {
                letra as char
            } else {
                letra.to_ascii_lowercase() as char
            }
        })
        .collect()
}

fn recortar(secuencia: &[u8], maximo: usize) -> String {
    if secuencia.len() <= maximo {
        String::from_utf8_lossy(secuencia).into_owned()
    } else {
        format!("{}…", String::from_utf8_lossy(&secuencia[..maximo]))
    }
}

fn memoria_de(indice: &IndiceFm) -> usize {
    indice.bloques.len() * std::mem::size_of::<BloqueDelIndice>() + indice.muestras.len() * 4
}

fn megabytes(bytes: usize) -> String {
    if bytes < 1_048_576 {
        format!("{}\u{a0}kB", decimal(bytes as f64 / 1024.0, 1))
    } else {
        format!("{}\u{a0}MB", decimal(bytes as f64 / 1_048_576.0, 1))
    }
}

fn con_diferencias(d: usize) -> String {
    match d {
        0 => "sin diferencias".to_string(),
        1 => "con hasta 1 diferencia".to_string(),
        _ => format!("con hasta {d} diferencias"),
    }
}

struct Sobrerrepresentado {
    numero: usize,
    razon: f64,
    esperadas: f64,
}

struct Medidas {
    total: usize,
    contar: Duration,
    ubicar: Duration,
    por_aparicion: Option<Duration>,
    aceleracion: f64,
}

fn buscar_con_indice(entrada: &Entrada) -> Result<Salida, String> {
    let pedido = leer_pedido(entrada)?;
    let usado = preparar_indice(&pedido)?;
    let indice = usado.doble();
    let resultados = indice_fm::buscar_patrones(indice, &pedido.patrones, pedido.d)?;
    let mut tiempos = Vec::with_capacity(pedido.patrones.len());
    let limite = tiempo_de_medida(pedido.patrones.len());
    for patron in &pedido.patrones {
        tiempos.push(medir_indice(indice, patron, pedido.d, limite)?);
    }
    let comparacion = comparar_con_el_escaneo(
        &pedido.texto,
        &pedido.patrones,
        &resultados,
        &tiempos,
        pedido.d,
        TIEMPO_DEL_ESCANEO,
    );
    let total: usize = resultados.iter().map(Vec::len).sum();
    let suma_contar: Duration = tiempos.iter().map(|tiempo| tiempo.contar).sum();
    let suma_ubicar: Duration = tiempos.iter().map(|tiempo| tiempo.ubicar).sum();
    let medidas = Medidas {
        total,
        contar: promedio(tiempos.iter().map(|tiempo| tiempo.contar)),
        ubicar: promedio(tiempos.iter().map(|tiempo| tiempo.ubicar)),
        por_aparicion: if total > 0 {
            Some(suma_ubicar.saturating_sub(suma_contar) / total as u32)
        } else {
            None
        },
        aceleracion: veces_mas_rapido(comparacion.clasico, comparacion.indice),
    };
    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("Apariciones (ambas hebras)", miles(total)));
    salida.agregar(Bloque::valor(
        "Contar con el índice",
        format!(
            "{} por patrón en promedio (el más rápido: {})",
            duracion(medidas.contar),
            duracion(
                tiempos
                    .iter()
                    .map(|tiempo| tiempo.contar)
                    .min()
                    .unwrap_or_default()
            )
        ),
    ));
    salida.agregar(Bloque::valor(
        "Contar y ubicar con el índice",
        format!("{} por patrón en promedio", duracion(medidas.ubicar)),
    ));
    if let Some(por_aparicion) = medidas.por_aparicion {
        salida.agregar(Bloque::valor(
            "Ubicar cada aparición (LF hasta una muestra)",
            duracion(por_aparicion),
        ));
    }
    let nombre_clasico = nombre_del_clasico(pedido.d);
    salida.agregar(Bloque::valor(
        &format!("Escaneo clásico («{nombre_clasico}», ambas hebras)"),
        format!(
            "{} por patrón (medido en {} {})",
            duracion(promedio(comparacion.clasico_de_cada.iter().copied())),
            miles(comparacion.patrones),
            if comparacion.patrones == 1 {
                "patrón"
            } else {
                "patrones"
            }
        ),
    ));
    salida.agregar(Bloque::valor(
        "Aceleración del índice (contar y ubicar contra el escaneo)",
        veces_con_decimales(medidas.aceleracion),
    ));
    salida.agregar(Bloque::valor(
        "Construcción del índice (texto y texto al revés)",
        if usado.ya_estaba() {
            format!(
                "{} una sola vez, al preparar la app",
                duracion(usado.construccion())
            )
        } else {
            duracion(usado.construccion())
        },
    ));
    salida.agregar(Bloque::valor(
        "¿Mismas posiciones que el escaneo clásico?",
        if comparacion.coinciden {
            format!("Sí, {}", patrones_comparados(comparacion.patrones))
        } else {
            "No".to_string()
        },
    ));
    let proporcion = proporciones(&pedido.texto);
    let mut filas_de_resumen = Vec::new();
    let mut sobrerrepresentados = Vec::new();
    for (numero, (patron, aciertos)) in pedido.patrones.iter().zip(&resultados).enumerate() {
        let esperadas = esperadas_al_azar(&pedido.texto, patron, pedido.d, &proporcion);
        let observadas = aciertos.len();
        let razon = observadas as f64 / esperadas.max(1e-12);
        if observadas as f64 > esperadas
            && razon >= 2.0
            && poisson_al_menos(esperadas, observadas) < SIGNIFICANCIA
        {
            sobrerrepresentados.push(Sobrerrepresentado {
                numero,
                razon,
                esperadas,
            });
        }
        let clasico = comparacion
            .clasico_de_cada
            .get(numero)
            .map(|tiempo| duracion(*tiempo))
            .unwrap_or_else(|| "—".to_string());
        filas_de_resumen.push(vec![
            recortar(patron, 40),
            miles(patron.len()),
            miles(posiciones_de_hebra(aciertos, b'+').len()),
            miles(posiciones_de_hebra(aciertos, b'-').len()),
            cantidad_esperada(esperadas),
            duracion(tiempos[numero].contar),
            duracion(tiempos[numero].ubicar),
            clasico,
        ]);
    }
    sobrerrepresentados.sort_by(|a, b| b.razon.total_cmp(&a.razon));
    salida.agregar(Bloque::tabla(
        "Por patrón",
        &[
            "Patrón",
            "Largo",
            "Hebra +",
            "Hebra −",
            "Esperadas al azar",
            "Contar (índice)",
            "Contar y ubicar (índice)",
            "Escaneo clásico",
        ],
        filas_de_resumen,
    ));
    agregar_apariciones(&mut salida, &pedido, &resultados, total);
    agregar_por_dentro(&mut salida, &usado, pedido.texto.len());
    if total > 0 {
        agregar_tramos(&mut salida, &resultados, pedido.texto.len());
    }
    explicar_busqueda(
        &mut salida,
        &pedido,
        &usado,
        &medidas,
        &comparacion,
        &sobrerrepresentados,
        &resultados,
    );
    Ok(salida)
}

fn nombre_del_clasico(d: usize) -> &'static str {
    if d == 0 {
        "Dónde aparece un patrón"
    } else {
        "Búsqueda aproximada"
    }
}

fn agregar_apariciones(
    salida: &mut Salida,
    pedido: &Pedido,
    resultados: &[Vec<Acierto>],
    total: usize,
) {
    let mut filas = Vec::new();
    'patrones: for (patron, aciertos) in pedido.patrones.iter().zip(resultados) {
        for acierto in aciertos {
            if filas.len() == MAXIMO_FILAS {
                break 'patrones;
            }
            let mut fila = vec![
                recortar(patron, 40),
                hebra_texto(acierto.hebra).to_string(),
                miles(acierto.posicion),
            ];
            if let Some(inicio) = pedido.inicio_ecoli {
                fila.push(miles(inicio + acierto.posicion));
            }
            fila.push(acierto.diferencias.to_string());
            fila.push(recortar(
                texto_hallado(&pedido.texto, patron, acierto).as_bytes(),
                60,
            ));
            filas.push(fila);
        }
    }
    let columnas: &[&str] = if pedido.inicio_ecoli.is_some() {
        &[
            "Patrón",
            "Hebra",
            "Posición (desde 0)",
            "Coordenada en E. coli (desde 1)",
            "Diferencias",
            "Lo que dice el texto (minúscula = diferencia)",
        ]
    } else {
        &[
            "Patrón",
            "Hebra",
            "Posición (desde 0)",
            "Diferencias",
            "Lo que dice el texto (minúscula = diferencia)",
        ]
    };
    salida.agregar(Bloque::tabla_con_total(
        "Apariciones",
        columnas,
        filas,
        total,
    ));
}

fn agregar_por_dentro(salida: &mut Salida, usado: &IndiceUsado, longitud: usize) {
    let directo = usado.directo();
    let memoria = memoria_de(directo) + memoria_de(&usado.doble().del_reves);
    salida.agregar(Bloque::tabla(
        "El índice por dentro",
        &["Parte", "Valor"],
        vec![
            vec![
                "Filas (sufijos de texto$)".to_string(),
                miles(directo.filas),
            ],
            vec![
                format!("Bloques de {PASO_DE_CONTROL} filas con puntos de control"),
                miles(directo.bloques.len()),
            ],
            vec![
                format!("Muestras del arreglo de sufijos (1 de cada {PASO_DE_MUESTRA} posiciones)"),
                miles(directo.muestras.len()),
            ],
            vec![
                "Memoria del índice (texto y texto al revés)".to_string(),
                megabytes(memoria),
            ],
            vec![
                "Un arreglo de sufijos completo ocuparía".to_string(),
                megabytes(4 * (longitud + 1)),
            ],
        ],
    ));
    let mut filas = Vec::new();
    for simbolo in 0..SIMBOLOS {
        filas.push(vec![
            letra_del_simbolo(simbolo).to_string(),
            miles(directo.primera_fila[simbolo]),
            miles(directo.primera_fila[simbolo + 1] - directo.primera_fila[simbolo]),
        ]);
    }
    salida.agregar(Bloque::tabla(
        "Columna F: dónde empieza cada letra (primera_fila)",
        &["Letra", "Primera fila", "Filas"],
        filas,
    ));
}

fn agregar_tramos(salida: &mut Salida, resultados: &[Vec<Acierto>], longitud: usize) {
    let tramos = TRAMOS.min(longitud.max(1));
    let ancho = longitud.div_ceil(tramos).max(1);
    let mut cuentas = vec![0.0; tramos];
    for acierto in resultados.iter().flatten() {
        cuentas[(acierto.posicion / ancho).min(tramos - 1)] += 1.0;
    }
    let categorias = (0..tramos)
        .map(|tramo| {
            format!(
                "{}–{}",
                miles(tramo * ancho),
                miles(((tramo + 1) * ancho).min(longitud))
            )
        })
        .collect();
    salida.agregar(Bloque::barras(
        "Apariciones a lo largo del texto",
        "Tramo (posición desde 0)",
        "Apariciones",
        categorias,
        cuentas,
    ));
}

fn explicar_busqueda(
    salida: &mut Salida,
    pedido: &Pedido,
    usado: &IndiceUsado,
    medidas: &Medidas,
    comparacion: &Comparacion,
    sobrerrepresentados: &[Sobrerrepresentado],
    resultados: &[Vec<Acierto>],
) {
    let nombre_clasico = nombre_del_clasico(pedido.d);
    let clasico = promedio(comparacion.clasico_de_cada.iter().copied());
    salida.explicar(format!(
        "El índice encontró {} apariciones de {} {} en ambas hebras, {}, sin recorrer el texto: cada letra del patrón solo achica el intervalo [arriba, abajo) de filas de las rotaciones ordenadas.",
        miles(medidas.total),
        miles(pedido.patrones.len()),
        if pedido.patrones.len() == 1 { "patrón" } else { "patrones" },
        con_diferencias(pedido.d)
    ));
    let por_aparicion = match medidas.por_aparicion {
        Some(tiempo) => format!(", y ubicar cada aparición unos {}", duracion(tiempo)),
        None => String::new(),
    };
    if medidas.aceleracion >= 10.0 {
        salida.explicar(format!(
            "Contar con el índice tomó {} por patrón{por_aparicion}; el escaneo de «{nombre_clasico}» tardó {} por patrón. Contando y ubicando, el índice fue {} veces más rápido: su costo depende del largo del patrón y de cuántas veces aparece, nunca del largo del texto.",
            duracion(medidas.contar),
            duracion(clasico),
            miles(medidas.aceleracion.round() as usize)
        ));
    } else if medidas.aceleracion >= 1.0 {
        salida.explicar(format!(
            "Contar con el índice tomó {} por patrón{por_aparicion}, y en total fue solo {} veces más rápido que «{nombre_clasico}» ({} por patrón): con textos cortos o patrones muy frecuentes pesa más ubicar cada aparición (caminar con LF hasta una fila con muestra) que recorrer el texto.",
            duracion(medidas.contar),
            decimal(medidas.aceleracion, 1),
            duracion(clasico)
        ));
    } else {
        salida.explicar(format!(
            "Aquí el escaneo de «{nombre_clasico}» fue más rápido ({} contra {} por patrón): el texto es tan corto, o las apariciones tantas, que ubicarlas una por una cuesta más que recorrerlo. Contar sin ubicar tomó {}.",
            duracion(clasico),
            duracion(medidas.ubicar),
            duracion(medidas.contar)
        ));
    }
    if usado.ya_estaba() {
        salida.explicar(format!(
            "El índice de E. coli (texto y texto al revés) se construyó una sola vez al preparar la app, en {}, y se reutiliza en cada búsqueda: así trabajan BWA y Bowtie con millones de lecturas.",
            duracion(usado.construccion())
        ));
    } else {
        let ahorro = clasico.as_secs_f64() - medidas.ubicar.as_secs_f64();
        if ahorro > 0.0 {
            let necesarios = (usado.construccion().as_secs_f64() / ahorro)
                .ceil()
                .max(1.0);
            salida.explicar(format!(
                "Construir el índice de tu texto tomó {}; con lo que ahorra cada búsqueda, compensa a partir de unos {} patrones. Por eso los mapeadores lo construyen una vez y lo reutilizan.",
                duracion(usado.construccion()),
                miles(necesarios as usize)
            ));
        }
    }
    if comparacion.coinciden {
        salida.explicar(format!(
            "Las posiciones coinciden exactamente con el escaneo clásico {}: el índice es exacto, no una aproximación.",
            patrones_comparados(comparacion.patrones)
        ));
    } else {
        salida.explicar(
            "Las posiciones NO coinciden con el escaneo clásico: eso no debería pasar; revisa el texto y los patrones.",
        );
    }
    if let Some(primero) = sobrerrepresentados.first() {
        let patron = &pedido.patrones[primero.numero];
        let mut frase = format!(
            "{} aparece {} veces contra {} esperadas si las bases se mezclaran al azar con la misma composición (Poisson, p < 0.01): está sobrerrepresentado.",
            recortar(patron, 40),
            miles(resultados[primero.numero].len()),
            cantidad_esperada(primero.esperadas)
        );
        if patron == CHI && pedido.inicio_ecoli.is_some() {
            frase.push_str(" Es el sitio Chi, la señal que la enzima RecBCD reconoce para reparar el ADN por recombinación.");
        } else if patron.len() >= LARGO_DE_REPETICION {
            frase.push_str(" Un patrón tan largo repetido no es azar: son copias de un mismo elemento del genoma (en E. coli, por ejemplo, los siete operones de ARN ribosomal).");
        }
        salida.explicar(frase);
    }
    let chi = sobrerrepresentados
        .iter()
        .skip(1)
        .find(|elegido| pedido.patrones[elegido.numero] == CHI);
    if let Some(elegido) = chi
        && pedido.inicio_ecoli.is_some()
        && salida.interpretacion.len() < 6
    {
        salida.explicar(format!(
            "También {} está sobrerrepresentado ({} veces lo esperado): es el sitio Chi, la señal que la enzima RecBCD reconoce para reparar el ADN por recombinación.",
            String::from_utf8_lossy(CHI),
            decimal(elegido.razon, 1)
        ));
    }
    if pedido.d > 0 && salida.interpretacion.len() < 6 {
        salida.explicar(format!(
            "Con d = {} el retroceso prueba cambiar letras del patrón, y el índice del texto al revés da una cota mínima de diferencias para cada prefijo (el arreglo D de BWA) que poda ramas sin futuro; aun así las ramas crecen como (3·m)^d, por eso los mapeadores pasan a semillas exactas cuando hay más errores.",
            pedido.d
        ));
    }
}

fn transformada(entrada: &Entrada) -> Result<Salida, String> {
    let texto = leer_texto_bwt(entrada)?;
    let inicio = Instant::now();
    let (bwt, recuperado) = indice_fm::transformada_e_inversa(&texto)?;
    let tiempo = inicio.elapsed();
    let mut salida = Salida::default();
    salida.agregar(Bloque::texto(
        "BWT (última columna de las rotaciones ordenadas)",
        recortar(&bwt, MAXIMO_MOSTRADO),
    ));
    salida.agregar(Bloque::valor(
        "¿La inversa devolvió el texto?",
        if recuperado == texto {
            "Sí, letra por letra"
        } else {
            "No"
        },
    ));
    if texto.len() <= MAXIMO_ROTACIONES {
        let sufijos = indice_fm::arreglo_de_sufijos(&texto)?;
        let filas = sufijos
            .iter()
            .enumerate()
            .map(|(fila, &inicio)| {
                let girada = rotacion(&texto, inicio as usize);
                vec![
                    fila.to_string(),
                    inicio.to_string(),
                    girada.clone(),
                    girada.chars().next().unwrap_or('$').to_string(),
                    girada.chars().last().unwrap_or('$').to_string(),
                ]
            })
            .collect();
        salida.agregar(Bloque::tabla(
            "Rotaciones de texto$ ordenadas",
            &["Fila", "Empieza en", "Rotación", "F", "L (BWT)"],
            filas,
        ));
    }
    let estadisticas = Rachas::medir(&texto, &bwt);
    salida.agregar(Bloque::tabla(
        "Rachas de letras iguales",
        &["Medida", "Texto", "BWT"],
        vec![
            vec!["Letras".to_string(), miles(texto.len()), miles(bwt.len())],
            vec![
                "Rachas".to_string(),
                miles(estadisticas.en_el_texto),
                miles(estadisticas.en_la_bwt),
            ],
            vec![
                "Largo medio de racha".to_string(),
                decimal(estadisticas.medio_del_texto, 2),
                decimal(estadisticas.medio_de_la_bwt, 2),
            ],
            vec![
                "Racha más larga".to_string(),
                miles(racha_mas_larga(&texto)),
                miles(racha_mas_larga(&bwt)),
            ],
            vec![
                "Largo medio si las letras fueran al azar".to_string(),
                decimal(estadisticas.al_azar, 2),
                decimal(estadisticas.al_azar, 2),
            ],
        ],
    ));
    let (categorias, valores) = histograma_de_rachas(&bwt);
    salida.agregar(Bloque::barras(
        "Rachas de la BWT según su largo",
        "Largo de la racha",
        "Rachas",
        categorias,
        valores,
    ));
    salida.agregar(Bloque::valor(
        "Tiempo (transformada e inversa)",
        duracion(tiempo),
    ));
    explicar_bwt(&mut salida, &texto, &recuperado, &estadisticas, tiempo);
    Ok(salida)
}

fn leer_texto_bwt(entrada: &Entrada) -> Result<Vec<u8>, String> {
    let texto = entrada.secuencia("texto")?;
    if texto.len() > MAXIMO_TEXTO_BWT {
        return Err(format!(
            "El texto tiene {} letras; aquí se admiten hasta {}.",
            miles(texto.len()),
            miles(MAXIMO_TEXTO_BWT)
        ));
    }
    if texto.contains(&CENTINELA) {
        return Err("El texto no puede llevar «$»: ese símbolo marca el final.".to_string());
    }
    Ok(texto)
}

fn rotacion(texto: &[u8], inicio: usize) -> String {
    let mut rotacion = String::with_capacity(texto.len() + 1);
    rotacion.push_str(&String::from_utf8_lossy(&texto[inicio.min(texto.len())..]));
    rotacion.push('$');
    rotacion.push_str(&String::from_utf8_lossy(&texto[..inicio.min(texto.len())]));
    rotacion
}

struct Rachas {
    en_el_texto: usize,
    en_la_bwt: usize,
    medio_del_texto: f64,
    medio_de_la_bwt: f64,
    al_azar: f64,
}

impl Rachas {
    fn medir(texto: &[u8], bwt: &[u8]) -> Rachas {
        let en_el_texto = indice_fm::rachas(texto);
        let en_la_bwt = indice_fm::rachas(bwt);
        let mut cuentas = vec![0usize; LETRAS_DE_UN_BYTE];
        for &letra in texto {
            cuentas[letra as usize] += 1;
        }
        let total = texto.len().max(1) as f64;
        let repetir: f64 = cuentas
            .iter()
            .map(|&cuenta| (cuenta as f64 / total).powi(2))
            .sum();
        Rachas {
            en_el_texto,
            en_la_bwt,
            medio_del_texto: texto.len() as f64 / en_el_texto.max(1) as f64,
            medio_de_la_bwt: bwt.len() as f64 / en_la_bwt.max(1) as f64,
            al_azar: (1.0 / (1.0 - repetir)).min(texto.len().max(1) as f64),
        }
    }
}

fn racha_mas_larga(secuencia: &[u8]) -> usize {
    let mut mejor = 0;
    let mut actual = 0;
    for (posicion, &letra) in secuencia.iter().enumerate() {
        if posicion > 0 && secuencia[posicion - 1] == letra {
            actual += 1;
        } else {
            actual = 1;
        }
        mejor = mejor.max(actual);
    }
    mejor
}

fn histograma_de_rachas(secuencia: &[u8]) -> (Vec<String>, Vec<f64>) {
    let categorias = ["1", "2", "3", "4", "5–8", "9–16", "17 o más"];
    let mut valores = vec![0.0; categorias.len()];
    let mut actual = 0;
    for (posicion, &letra) in secuencia.iter().enumerate() {
        actual += 1;
        let termina = posicion + 1 == secuencia.len() || secuencia[posicion + 1] != letra;
        if termina {
            let casilla = match actual {
                1..=4 => actual - 1,
                5..=8 => 4,
                9..=16 => 5,
                _ => 6,
            };
            valores[casilla] += 1.0;
            actual = 0;
        }
    }
    (
        categorias
            .iter()
            .map(|categoria| categoria.to_string())
            .collect(),
        valores,
    )
}

fn explicar_bwt(
    salida: &mut Salida,
    texto: &[u8],
    recuperado: &[u8],
    rachas: &Rachas,
    tiempo: Duration,
) {
    if recuperado == texto {
        salida.explicar(
            "La BWT es la última columna de las rotaciones de texto$ ordenadas; la inversa la deshizo con el mapeo LF (la k-ésima G de la columna L es la k-ésima G de la columna F) y devolvió exactamente el texto: reordenar así no pierde información.",
        );
    } else {
        salida.explicar("La inversa NO devolvió el texto: eso no debería pasar.");
    }
    salida.explicar(format!(
        "La BWT tiene {} rachas de letras iguales contra {} del texto: {} letras por racha en promedio contra {}.",
        miles(rachas.en_la_bwt),
        miles(rachas.en_el_texto),
        decimal(rachas.medio_de_la_bwt, 2),
        decimal(rachas.medio_del_texto, 2)
    ));
    if rachas.medio_de_la_bwt >= FACTOR_DE_AZAR * rachas.al_azar {
        salida.explicar(format!(
            "Eso supera con creces lo que darían letras al azar con la misma composición ({} por racha): los contextos repetidos del texto ponen letras iguales juntas, y un compresor por rachas (como bzip2 después de la BWT) guardaría solo {} pares letra-largo.",
            decimal(rachas.al_azar, 2),
            miles(rachas.en_la_bwt)
        ));
    } else {
        salida.explicar(format!(
            "Es casi lo que darían letras al azar con la misma composición ({} por racha): este texto tiene pocas repeticiones y la BWT apenas ayuda a comprimirlo. Pasa con un solo genoma; en colecciones de genomas parecidos las rachas crecen mucho más lento que el texto, y en eso se basa el r-index (Gagie, Navarro y Prezza, 2018).",
            decimal(rachas.al_azar, 2)
        ));
    }
    salida.explicar(format!(
        "Ordenar los {} sufijos con SA-IS y deshacer la transformada tomó {}: el tiempo crece en línea recta con el largo, sin comparar rotaciones completas. Con puntos de control y una muestra del arreglo de sufijos, esta columna es el «Índice FM: buscar sin recorrer el texto», que encuentra patrones sin el recorrido de «Dónde aparece un patrón».",
        miles(texto.len() + 1),
        duracion(tiempo)
    ));
}

struct Candidato {
    simbolo: usize,
    arriba: usize,
    abajo: usize,
    diferencias: usize,
    apilado: bool,
}

struct Visita {
    rama: Rama,
    camino: Vec<u8>,
    esperado: usize,
    cota: usize,
    candidatos: Vec<Candidato>,
    pendientes: usize,
}

fn recorrer_ramas(
    indice: &IndiceFm,
    patron: &[u8],
    d: usize,
    cotas: &[usize],
    mut al_visitar: impl FnMut(&Visita),
) -> Vec<(usize, usize, usize)> {
    let mut hallados = Vec::new();
    let mut pendientes = vec![Rama {
        faltan: patron.len(),
        arriba: 0,
        abajo: indice.filas,
        diferencias: 0,
    }];
    let mut caminos: Vec<Vec<u8>> = vec![Vec::new()];
    while let (Some(rama), Some(camino)) = (pendientes.pop(), caminos.pop()) {
        if rama.faltan == 0 {
            hallados.push((rama.arriba, rama.abajo, rama.diferencias));
            al_visitar(&Visita {
                rama,
                camino,
                esperado: 0,
                cota: 0,
                candidatos: Vec::new(),
                pendientes: pendientes.len(),
            });
            continue;
        }
        let esperado = indice_fm::simbolo_de(patron[rama.faltan - 1]) as usize;
        let cota = if rama.faltan >= 2 {
            cotas[rama.faltan - 2]
        } else {
            0
        };
        let antes_de_arriba = indice.ocurrencias_de_todas(rama.arriba);
        let antes_de_abajo = indice.ocurrencias_de_todas(rama.abajo);
        let mut candidatos = Vec::new();
        let otros = (1..SIMBOLOS).rev().filter(|&simbolo| simbolo != esperado);
        for simbolo in otros.chain([esperado]) {
            let diferencias = rama.diferencias + usize::from(simbolo != esperado);
            let arriba = indice.primera_fila[simbolo] + antes_de_arriba[simbolo];
            let abajo = indice.primera_fila[simbolo] + antes_de_abajo[simbolo];
            let apilado = diferencias + cota <= d && arriba < abajo;
            if apilado {
                pendientes.push(Rama {
                    faltan: rama.faltan - 1,
                    arriba,
                    abajo,
                    diferencias,
                });
                let mut nuevo = vec![LETRAS_DEL_INDICE[simbolo]];
                nuevo.extend_from_slice(&camino);
                caminos.push(nuevo);
            }
            candidatos.push(Candidato {
                simbolo,
                arriba,
                abajo,
                diferencias,
                apilado,
            });
        }
        al_visitar(&Visita {
            rama,
            camino,
            esperado,
            cota,
            candidatos,
            pendientes: pendientes.len(),
        });
    }
    hallados
}

fn camino_marcado(camino: &[u8], patron: &[u8]) -> String {
    let desde = patron.len() - camino.len();
    let marcado: String = camino
        .iter()
        .zip(&patron[desde..])
        .map(|(&letra, &esperada)| {
            if letra == esperada {
                letra as char
            } else {
                letra.to_ascii_lowercase() as char
            }
        })
        .collect();
    if marcado.is_empty() {
        "(nada)".to_string()
    } else {
        marcado
    }
}

fn plural(cantidad: usize, singular: &str, varias: &str) -> String {
    if cantidad == 1 {
        format!("1 {singular}")
    } else {
        format!("{} {varias}", miles(cantidad))
    }
}

fn sufijo_corto(texto: &[u8], inicio: usize) -> String {
    if inicio >= texto.len() {
        return "$".to_string();
    }
    let fin = (inicio + LARGO_DEL_SUFIJO).min(texto.len());
    let mut sufijo = String::from_utf8_lossy(&texto[inicio..fin]).into_owned();
    if fin == texto.len() {
        sufijo.push('$');
    } else {
        sufijo.push('…');
    }
    sufijo
}

fn filas_a_mostrar(filas: usize, arriba: usize, abajo: usize) -> Vec<Option<usize>> {
    if filas <= MAXIMO_ROTACIONES + 1 {
        return (0..filas).map(Some).collect();
    }
    if abajo <= arriba {
        let centro = arriba.min(filas - 1);
        return (centro.saturating_sub(2)..(centro + 3).min(filas))
            .map(Some)
            .collect();
    }
    if abajo - arriba <= 8 {
        return (arriba.saturating_sub(2)..(abajo + 2).min(filas))
            .map(Some)
            .collect();
    }
    let mut lista: Vec<Option<usize>> = Vec::new();
    if arriba > 0 {
        lista.push(Some(arriba - 1));
    }
    lista.extend((arriba..arriba + 3).map(Some));
    lista.push(None);
    lista.extend((abajo - 3..abajo).map(Some));
    if abajo < filas {
        lista.push(Some(abajo));
    }
    lista
}

fn tabla_de_filas(
    indice: &IndiceFm,
    texto: &[u8],
    arriba: usize,
    abajo: usize,
    siguiente: Option<usize>,
) -> (Vec<Vec<String>>, Vec<(usize, usize)>) {
    let corto = indice.filas <= MAXIMO_ROTACIONES + 1;
    let mut filas = Vec::new();
    let mut resaltar = Vec::new();
    for (lugar, fila) in filas_a_mostrar(indice.filas, arriba, abajo)
        .into_iter()
        .enumerate()
    {
        let Some(fila) = fila else {
            filas.push(vec![
                "…".to_string(),
                "…".to_string(),
                format!("{} filas más", miles(abajo - arriba - 6)),
                "…".to_string(),
            ]);
            continue;
        };
        let inicio = indice.posicion_de_fila(fila);
        let sufijo = if corto {
            rotacion(texto, inicio)
        } else {
            sufijo_corto(texto, inicio)
        };
        let letra_l = indice.simbolo_de_l(fila) as usize;
        filas.push(vec![
            miles(fila),
            letra_del_simbolo(indice.simbolo_de_f(fila) as usize).to_string(),
            sufijo,
            letra_del_simbolo(letra_l).to_string(),
        ]);
        if (arriba..abajo).contains(&fila) {
            resaltar.push((lugar, 1));
            resaltar.push((lugar, 2));
            if siguiente == Some(letra_l) {
                resaltar.push((lugar, 3));
            }
        }
    }
    (filas, resaltar)
}

fn columnas_de_filas(indice: &IndiceFm) -> [&'static str; 4] {
    if indice.filas <= MAXIMO_ROTACIONES + 1 {
        ["fila", "F", "rotación", "L"]
    } else {
        ["fila", "F", "sufijo (primeras letras)", "L"]
    }
}

fn marcas_del_patron(largo: usize, faltan: usize) -> Vec<(usize, usize, Marca)> {
    let mut marcas = Vec::new();
    if faltan < largo {
        marcas.push((faltan, largo, Marca::Hallado));
    }
    if faltan > 0 {
        marcas.push((faltan - 1, faltan, Marca::Actual));
    }
    marcas
}

fn pasos_indice_fm(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let pedido = leer_pedido(entrada)?;
    let usado = preparar_indice(&pedido)?;
    let indice = usado.doble();
    let directo = &indice.directo;
    let texto: &[u8] = &pedido.texto;
    let patron = &pedido.patrones[0];
    let d = pedido.d;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    narracion.agregar(
        Pizarra::nueva()
            .variable("patrones", miles(pedido.patrones.len()))
            .variable("patron", recortar(patron, 40))
            .variable("n", miles(texto.len()))
            .variable("d", d)
            .secuencia(patron, &[])
            .paso(
                "revisar_patrones(patrones, d)?;",
                format!(
                    "Vas a buscar {} {} en un texto de {} letras, en ambas hebras y {}. En vez de recorrer el texto, el código lo ordena una vez (el índice FM) y después cada letra del patrón solo achica un intervalo de filas. Sigues el primero: {}.",
                    miles(pedido.patrones.len()),
                    if pedido.patrones.len() == 1 { "patrón" } else { "patrones" },
                    miles(texto.len()),
                    con_diferencias(d),
                    recortar(patron, 40)
                ),
            ),
    );
    narrar_construccion(&mut narracion, &usado, texto);
    let directos = if d == 0 {
        narrar_busqueda_exacta(&mut narracion, directo, texto, patron)
    } else {
        let cotas = indice.cotas_de_diferencias(patron);
        narrar_cotas(&mut narracion, patron, &cotas, d);
        narrar_retroceso(&mut narracion, directo, texto, patron, d, &cotas)
    };
    let mut intervalos: Vec<Intervalo> = directos
        .iter()
        .map(|&(arriba, abajo, diferencias)| Intervalo {
            hebra: b'+',
            arriba,
            abajo,
            diferencias,
        })
        .collect();
    let inverso = complemento_inverso(patron);
    if inverso == *patron {
        narracion.agregar(
            Pizarra::nueva()
                .variable("inverso", recortar(&inverso, 40))
                .paso(
                    "if inverso != patron {",
                    format!(
                        "{} es igual a su complemento inverso (un palíndromo): cada aparición ya está en las dos hebras, así que no se busca dos veces.",
                        recortar(patron, 40)
                    ),
                ),
        );
    } else {
        let inversos = indice_fm::buscar_en_una_hebra(indice, &inverso, d, b'-');
        let filas: usize = inversos.iter().map(|i| i.abajo - i.arriba).sum();
        narracion.agregar(
            Pizarra::nueva()
                .variable("inverso", recortar(&inverso, 40))
                .variable("intervalos", inversos.len())
                .variable("filas", miles(filas))
                .secuencia(&inverso, &[])
                .paso(
                    "let inverso = complemento_inverso(patron);",
                    format!(
                        "Para la hebra − haces lo mismo con el complemento inverso, {}: deja {} {} con {} {} en total.",
                        recortar(&inverso, 40),
                        miles(inversos.len()),
                        if inversos.len() == 1 { "intervalo" } else { "intervalos" },
                        miles(filas),
                        if filas == 1 { "fila" } else { "filas" }
                    ),
                ),
        );
        intervalos.extend(inversos);
    }
    if let Some(primero) = intervalos.first() {
        narrar_ubicacion(&mut narracion, directo, &pedido, patron, primero.arriba);
    }
    narrar_otros_patrones(&mut narracion, indice, &pedido);
    let resultados = indice_fm::buscar_patrones(indice, &pedido.patrones, d)?;
    let total: usize = resultados.iter().map(Vec::len).sum();
    let ultimo = Pizarra::nueva()
        .variable("patrones", miles(pedido.patrones.len()))
        .variable("apariciones", miles(total))
        .variable("del_primero", miles(resultados[0].len()))
        .paso(
            "resultado.push(ubicar(&indice.directo, intervalos));",
            format!(
                "Listo: cada fila de cada intervalo se ubicó con LF y las muestras, y quedaron {} apariciones en total ({} del primer patrón), ordenadas por posición. En «Qué significa» ves cuánto más rápido fue que el escaneo clásico.",
                miles(total),
                miles(resultados[0].len())
            ),
        );
    Ok(narracion.terminar(ultimo))
}

fn narrar_construccion(narracion: &mut Narracion, usado: &IndiceUsado, texto: &[u8]) {
    let directo = usado.directo();
    let corto = directo.filas <= MAXIMO_ROTACIONES + 1;
    let (filas, _) = tabla_de_filas(
        directo,
        texto,
        0,
        if corto { directo.filas } else { PRIMERAS_FILAS },
        None,
    );
    let como = if usado.ya_estaba() {
        format!(
            "El índice de E. coli ya estaba en memoria: se construyó una sola vez al preparar la app, en {}.",
            duracion(usado.construccion())
        )
    } else {
        format!(
            "Construir el índice de tu texto (y el de tu texto al revés, en paralelo) tomó {}.",
            duracion(usado.construccion())
        )
    };
    let filas_mostradas = if corto {
        "Aquí ves todas las rotaciones".to_string()
    } else {
        "Aquí ves las primeras filas".to_string()
    };
    narracion.agregar(
        Pizarra::nueva()
            .variable("filas", miles(directo.filas))
            .tabla(
                "Rotaciones de texto$ ordenadas",
                &columnas_de_filas(directo),
                filas,
                &[],
            )
            .paso(
                "let sufijos = ordenar_sufijos(&simbolos, SIMBOLOS);",
                format!(
                    "{como} SA-IS ordenó los {} sufijos de texto$ en tiempo lineal; ordenar sufijos es lo mismo que ordenar las rotaciones, porque $ es único y el menor. {filas_mostradas}: F es la primera letra y L la última (la BWT).",
                    miles(directo.filas)
                ),
            ),
    );
    let mut filas = Vec::new();
    for simbolo in 0..SIMBOLOS {
        filas.push(vec![
            letra_del_simbolo(simbolo).to_string(),
            miles(directo.primera_fila[simbolo]),
            miles(directo.primera_fila[simbolo + 1] - directo.primera_fila[simbolo]),
        ]);
    }
    narracion.agregar(
        Pizarra::nueva()
            .variable("filas", miles(directo.filas))
            .tabla("primera_fila", &["letra", "primera fila", "filas"], filas, &[])
            .paso(
                "primera_fila: primeras_filas(&simbolos),",
                format!(
                    "Como F está ordenada, basta guardar en qué fila empieza cada letra: las A ocupan desde la fila {} hasta antes de la {}, y así con las demás.",
                    miles(directo.primera_fila[1]),
                    miles(directo.primera_fila[2])
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("bloques", miles(directo.bloques.len()))
            .variable("muestras", miles(directo.muestras.len()))
            .variable("memoria", megabytes(memoria_de(directo)))
            .paso(
                "bloques: armar_bloques(&simbolos, &sufijos),",
                format!(
                    "La columna L se guarda en bloques de {PASO_DE_CONTROL} filas: cada bloque recuerda cuántas letras de cada tipo había antes (puntos de control) y cuáles de sus filas guardan su posición en el texto (1 de cada {PASO_DE_MUESTRA}). Contar letras hasta una fila cuesta como mucho {} comparaciones.",
                    PASO_DE_CONTROL - 1
                ),
            ),
    );
}

fn narrar_busqueda_exacta(
    narracion: &mut Narracion,
    indice: &IndiceFm,
    texto: &[u8],
    patron: &[u8],
) -> Vec<(usize, usize, usize)> {
    let largo = patron.len();
    let mut arriba = 0;
    let mut abajo = indice.filas;
    let siguiente = Some(indice_fm::simbolo_de(patron[largo - 1]) as usize);
    let (filas, resaltar) = tabla_de_filas(indice, texto, arriba, abajo, siguiente);
    narracion.agregar(
        Pizarra::nueva()
            .variable("letra", "—")
            .variable("arriba", miles(arriba))
            .variable("abajo", miles(abajo))
            .variable("coincidencias", miles(abajo - arriba))
            .secuencia(patron, &marcas_del_patron(largo, largo))
            .tabla("filas", &columnas_de_filas(indice), filas, &resaltar)
            .paso(
                "let mut abajo = self.filas;",
                "Empiezas con todas las filas: el intervalo [arriba, abajo) contiene los sufijos que empiezan con lo ya leído del patrón, y al principio no has leído nada. El patrón se lee de derecha a izquierda.",
            ),
    );
    let mut resumidas = 0;
    for (leidas, &letra) in patron.iter().rev().enumerate() {
        let simbolo = indice_fm::simbolo_de(letra);
        let antes = (arriba, abajo);
        (arriba, abajo) = indice.paso_atras(simbolo, arriba, abajo);
        let faltan = largo - leidas - 1;
        let vacio = arriba >= abajo;
        if leidas >= LETRAS_DETALLADAS && !vacio && faltan > 0 {
            resumidas += 1;
            continue;
        }
        let primera = indice.primera_fila[simbolo as usize];
        let siguiente = if faltan > 0 {
            Some(indice_fm::simbolo_de(patron[faltan - 1]) as usize)
        } else {
            None
        };
        let (filas, resaltar) = tabla_de_filas(indice, texto, arriba, abajo, siguiente);
        let mut nota = format!(
            "{} Sus filas empiezan en primera_fila = {}; en L hay {} {} antes de la fila {} y {} antes de la {}, así que el nuevo intervalo es [{}, {}): {} {}.",
            if leidas == 0 {
                format!("Lees la última letra del patrón, {}.", letra as char)
            } else {
                format!(
                    "Lees la letra anterior, {}: ahora buscas «{}».",
                    letra as char,
                    String::from_utf8_lossy(&patron[faltan..])
                )
            },
            miles(primera),
            miles(arriba.saturating_sub(primera)),
            letra as char,
            miles(antes.0),
            miles(abajo.saturating_sub(primera)),
            miles(antes.1),
            miles(arriba),
            miles(abajo),
            miles(abajo.saturating_sub(arriba)),
            if abajo.saturating_sub(arriba) == 1 {
                "fila"
            } else {
                "filas"
            }
        );
        if resumidas > 0 {
            nota = format!(
                "(Se resumieron {} letras que siguieron achicando el intervalo.) {nota}",
                miles(resumidas)
            );
            resumidas = 0;
        }
        if vacio {
            nota.push_str(" El intervalo quedó vacío: el patrón no está en esta hebra.");
        } else if faltan == 0 {
            nota.push_str(" Ya leíste todo el patrón: cada fila del intervalo es una aparición.");
        }
        narracion.agregar(
            Pizarra::nueva()
                .variable("letra", letra as char)
                .variable("arriba", miles(arriba))
                .variable("abajo", miles(abajo))
                .variable("coincidencias", miles(abajo.saturating_sub(arriba)))
                .secuencia(patron, &marcas_del_patron(largo, faltan))
                .tabla("filas", &columnas_de_filas(indice), filas, &resaltar)
                .paso(
                    "(arriba, abajo) = self.paso_atras(simbolo_de(letra), arriba, abajo);",
                    nota,
                ),
        );
        if vacio {
            return Vec::new();
        }
    }
    vec![(arriba, abajo, 0)]
}

fn narrar_cotas(narracion: &mut Narracion, patron: &[u8], cotas: &[usize], d: usize) {
    let elementos: Vec<String> = cotas
        .iter()
        .enumerate()
        .map(|(posicion, cota)| {
            format!(
                "hasta la letra {} ({}): {}",
                posicion + 1,
                patron[posicion] as char,
                cota
            )
        })
        .take(MAXIMO_ROTACIONES * 2)
        .collect();
    let cambios: Vec<usize> = cotas
        .iter()
        .enumerate()
        .filter(|(posicion, cota)| **cota > 0 && (*posicion == 0 || cotas[posicion - 1] < **cota))
        .map(|(posicion, _)| posicion)
        .filter(|&posicion| posicion < MAXIMO_ROTACIONES * 2)
        .collect();
    let cota_final = cotas.last().copied().unwrap_or(0);
    narracion.agregar(
        Pizarra::nueva()
            .variable("cota", cota_final)
            .variable("d", d)
            .lista("cotas", elementos, &cambios)
            .paso(
                "cotas.push(cota);",
                format!(
                    "Antes de buscar con diferencias, el código lee el patrón de izquierda a derecha en el índice del texto al revés: cada vez que un tramo ya no aparece en el texto, ese tramo necesita al menos una diferencia y la cota sube (el arreglo D de BWA). {}",
                    if cota_final == 0 {
                        "Aquí la cota queda en 0 porque el patrón entero aparece tal cual: no podará nada, pero con patrones que no están en el texto ahorra muchas ramas.".to_string()
                    } else {
                        format!("Al final la cota es {cota_final}: una rama que ya gastó diferencias y no le alcanzan para lo que falta se poda sin explorarla.")
                    }
                ),
            ),
    );
}

fn nota_de_rama(visita: &Visita, patron: &[u8], d: usize) -> String {
    let rama = visita.rama;
    let filas = rama.abajo - rama.arriba;
    let marcado = camino_marcado(&visita.camino, patron);
    if rama.faltan == 0 {
        return format!(
            "Esta rama leyó todo el patrón: «{marcado}» (minúscula = letra cambiada) aparece en {} [{}, {}), con {}. Sacas de la pila la rama siguiente.",
            plural(filas, "fila", "filas"),
            miles(rama.arriba),
            miles(rama.abajo),
            plural(rama.diferencias, "diferencia", "diferencias")
        );
    }
    let letra = patron[rama.faltan - 1] as char;
    let hijo = visita
        .candidatos
        .iter()
        .find(|candidato| candidato.simbolo == visita.esperado);
    let filas_del_hijo = hijo.map_or(0, |candidato| {
        candidato.abajo.saturating_sub(candidato.arriba)
    });
    let variantes = visita
        .candidatos
        .iter()
        .filter(|candidato| candidato.apilado && candidato.simbolo != visita.esperado)
        .count();
    let sigue = if filas_del_hijo == 0 {
        format!("con la {letra} no queda ninguna fila")
    } else if hijo.is_some_and(|candidato| !candidato.apilado) {
        format!(
            "con la {letra} quedarían {}, pero la cota dice que lo que falta necesita más diferencias de las que quedan",
            plural(filas_del_hijo, "fila", "filas")
        )
    } else {
        format!(
            "con la {letra} quedan {}",
            plural(filas_del_hijo, "fila", "filas")
        )
    };
    let otras = if variantes > 0 {
        format!(
            "; {} con otra letra (una diferencia más) van a la pila",
            plural(variantes, "variante", "variantes")
        )
    } else if rama.diferencias == d {
        "; otra letra ya no cabe: la rama gastó sus d diferencias".to_string()
    } else {
        "; ninguna otra letra deja filas que valga la pena seguir".to_string()
    };
    if rama.faltan == patron.len() {
        return format!(
            "Empiezas con todas las filas, [0, {}). El patrón se lee de derecha a izquierda y su última letra es {letra}. Con los puntos de control cuentas a la vez cuántas A, C, G y T hay en L dentro del intervalo: {sigue}{otras}.",
            miles(rama.abajo)
        );
    }
    if rama.diferencias == 0 {
        return format!(
            "Las {} de [{}, {}) son los sufijos que empiezan con «{marcado}». Toca la {letra} (letra {} del patrón): {sigue}{otras}.",
            plural(filas, "fila", "filas"),
            miles(rama.arriba),
            miles(rama.abajo),
            rama.faltan
        );
    }
    format!(
        "Sacas de la pila la rama «{marcado}» ({}; minúscula = letra cambiada): {} de [{}, {}). Toca la {letra}: {sigue}{otras}.",
        plural(rama.diferencias, "diferencia", "diferencias"),
        plural(filas, "fila", "filas"),
        miles(rama.arriba),
        miles(rama.abajo)
    )
}

fn narrar_retroceso(
    narracion: &mut Narracion,
    indice: &IndiceFm,
    texto: &[u8],
    patron: &[u8],
    d: usize,
    cotas: &[usize],
) -> Vec<(usize, usize, usize)> {
    let largo = patron.len();
    let mut exactas = 0;
    let mut con_cambios = 0;
    let mut resumidas = 0;
    let mut hallazgos = 0;
    let mut visitadas = 0;
    let hallados = recorrer_ramas(indice, patron, d, cotas, |visita| {
        visitadas += 1;
        let rama = visita.rama;
        let detallar = if rama.diferencias == 0 {
            exactas += 1;
            exactas <= LETRAS_DETALLADAS || rama.faltan == 0
        } else {
            con_cambios += 1;
            con_cambios <= RAMAS_DETALLADAS
        };
        if rama.faltan == 0 {
            hallazgos += 1;
        }
        if !detallar || !narracion.hay_lugar() {
            resumidas += 1;
            return;
        }
        let mut nota = nota_de_rama(visita, patron, d);
        if resumidas > 0 {
            nota = format!(
                "(Se resumieron {} parecidas.) {nota}",
                plural(resumidas, "rama", "ramas")
            );
            resumidas = 0;
        }
        let pizarra = Pizarra::nueva()
            .variable("faltan", rama.faltan)
            .variable("diferencias", rama.diferencias)
            .variable("arriba", miles(rama.arriba))
            .variable("abajo", miles(rama.abajo))
            .variable("pendientes", visita.pendientes)
            .secuencia(patron, &marcas_del_patron(largo, rama.faltan));
        if rama.faltan == 0 {
            narracion.agregar(pizarra.paso("hallados.push(Intervalo {", nota));
            return;
        }
        let mut elementos = Vec::new();
        let mut resaltar = Vec::new();
        let mut candidatos: Vec<&Candidato> = visita
            .candidatos
            .iter()
            .filter(|candidato| {
                indice.primera_fila[candidato.simbolo] < indice.primera_fila[candidato.simbolo + 1]
            })
            .collect();
        candidatos.sort_by_key(|candidato| candidato.simbolo);
        for candidato in candidatos {
            let filas = candidato.abajo.saturating_sub(candidato.arriba);
            let decision = if candidato.apilado {
                "a la pila"
            } else if filas == 0 {
                "vacío"
            } else if candidato.diferencias > d {
                "se pasa de d"
            } else {
                "la cota la poda"
            };
            if candidato.apilado {
                resaltar.push(elementos.len());
            }
            elementos.push(format!(
                "{}{}: {}, {} dif. → {}",
                letra_del_simbolo(candidato.simbolo),
                if candidato.simbolo == visita.esperado {
                    " (la del patrón)"
                } else {
                    ""
                },
                plural(filas, "fila", "filas"),
                candidato.diferencias,
                decision
            ));
        }
        let (filas, resaltar_filas) = tabla_de_filas(
            indice,
            texto,
            rama.arriba,
            rama.abajo,
            Some(visita.esperado),
        );
        narracion.agregar(
            pizarra
                .variable("cota", visita.cota)
                .tabla("filas", &columnas_de_filas(indice), filas, &resaltar_filas)
                .lista("letras probadas", elementos, &resaltar)
                .paso("for simbolo in otros.chain([esperado]) {", nota),
        );
    });
    let filas: usize = hallados
        .iter()
        .map(|(arriba, abajo, _)| abajo - arriba)
        .sum();
    narracion.agregar(
        Pizarra::nueva()
            .variable("ramas", miles(visitadas))
            .variable("intervalos", hallados.len())
            .variable("filas", miles(filas))
            .paso(
                "while let Some(rama) = pendientes.pop() {",
                format!(
                    "La pila se vació: en la hebra + el retroceso visitó {} y {} el patrón, con {} en total.{}",
                    plural(visitadas, "rama", "ramas"),
                    if hallazgos == 1 {
                        "1 terminó".to_string()
                    } else {
                        format!("{} terminaron", miles(hallazgos))
                    },
                    plural(filas, "fila", "filas"),
                    if resumidas > 0 {
                        format!(" (Las últimas {} se resumieron.)", plural(resumidas, "rama", "ramas"))
                    } else {
                        String::new()
                    }
                ),
            ),
    );
    hallados
}

fn narrar_ubicacion(
    narracion: &mut Narracion,
    indice: &IndiceFm,
    pedido: &Pedido,
    patron: &[u8],
    fila: usize,
) {
    let texto: &[u8] = &pedido.texto;
    let posicion = indice.posicion_de_fila(fila);
    let fin = (posicion + patron.len()).min(texto.len());
    let mut actual = fila;
    let mut pasos = 0;
    while !indice.esta_marcada(actual) {
        let letra = letra_del_simbolo(indice.simbolo_de_l(actual) as usize);
        let siguiente = indice.lf(actual);
        let donde = posicion - pasos - 1;
        let nota = if pasos == 0 {
            format!(
                "Ya sabes cuántas veces aparece, falta saber dónde. La fila {} no guarda su posición (solo 1 de cada {PASO_DE_MUESTRA} la guarda). Su letra en L es {letra}, la que va justo antes en el texto: con LF saltas a la fila {} del sufijo que empieza una letra antes, y cuentas un paso.",
                miles(actual),
                miles(siguiente)
            )
        } else {
            format!(
                "Tampoco la fila {} tiene muestra: L dice {letra}, saltas a la fila {} y ya van {} hacia la izquierda.",
                miles(actual),
                miles(siguiente),
                plural(pasos + 1, "paso", "pasos")
            )
        };
        narracion.agregar(
            Pizarra::nueva()
                .variable("fila", miles(fila))
                .variable("actual", miles(siguiente))
                .variable("pasos", pasos + 1)
                .secuencia(
                    texto,
                    &[
                        (donde, donde + 1, Marca::Actual),
                        (posicion, fin, Marca::Hallado),
                    ],
                )
                .paso("actual = self.lf(actual);", nota),
        );
        actual = siguiente;
        pasos += 1;
    }
    let muestra = indice.muestras[indice.marcadas_hasta(actual)] as usize;
    let coordenada = match pedido.inicio_ecoli {
        Some(inicio) => format!(" (coordenada {} de E. coli)", miles(inicio + posicion)),
        None => String::new(),
    };
    narracion.agregar(
        Pizarra::nueva()
            .variable("fila", miles(fila))
            .variable("actual", miles(actual))
            .variable("pasos", pasos)
            .variable("posicion", miles(posicion))
            .secuencia(
                texto,
                &[
                    (muestra, muestra + 1, Marca::Comparado),
                    (posicion, fin, Marca::Hallado),
                ],
            )
            .paso(
                "self.muestras[self.marcadas_hasta(actual)] as usize + pasos",
                format!(
                    "La fila {} sí tiene muestra: su sufijo empieza en {}. Le sumas {} y sabes que la fila {} es la posición {}{coordenada}: ahí está la primera aparición.",
                    miles(actual),
                    miles(muestra),
                    plural(pasos, "paso", "pasos"),
                    miles(fila),
                    miles(posicion)
                ),
            ),
    );
}

fn narrar_otros_patrones(narracion: &mut Narracion, indice: &IndiceDoble, pedido: &Pedido) {
    let mut resumidos = 0;
    let mut filas_resumidas = 0;
    for (numero, patron) in pedido.patrones.iter().enumerate().skip(1) {
        let intervalos = indice_fm::buscar_en_ambas_hebras(indice, patron, pedido.d);
        let mas: usize = intervalos
            .iter()
            .filter(|i| i.hebra == b'+')
            .map(|i| i.abajo - i.arriba)
            .sum();
        let menos: usize = intervalos
            .iter()
            .filter(|i| i.hebra == b'-')
            .map(|i| i.abajo - i.arriba)
            .sum();
        if numero > PATRONES_DETALLADOS {
            resumidos += 1;
            filas_resumidas += mas + menos;
            continue;
        }
        narracion.agregar(
            Pizarra::nueva()
                .variable("patron", recortar(patron, 40))
                .variable("hebra_mas", miles(mas))
                .variable("hebra_menos", miles(menos))
                .secuencia(patron, &[])
                .paso(
                    "let intervalos = buscar_en_ambas_hebras(indice, patron, d);",
                    format!(
                        "El patrón {} ({}) se busca igual, letra por letra desde el final: {} {} en la hebra + y {} en la hebra −.",
                        numero + 1,
                        recortar(patron, 40),
                        miles(mas),
                        if mas == 1 { "aparición" } else { "apariciones" },
                        miles(menos)
                    ),
                ),
        );
    }
    if resumidos > 0 {
        narracion.agregar(
            Pizarra::nueva()
                .variable("patrones", miles(resumidos))
                .variable("filas", miles(filas_resumidas))
                .paso(
                    "for patron in patrones {",
                    format!(
                        "Los {} patrones restantes se buscan igual y suman {} apariciones.",
                        miles(resumidos),
                        miles(filas_resumidas)
                    ),
                ),
        );
    }
}

fn pasos_bwt(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let texto = leer_texto_bwt(entrada)?;
    let n = texto.len();
    let corto = n <= MAXIMO_ROTACIONES;
    let mut narracion = Narracion::nueva(MAXIMO_PASOS);
    let mut simbolos: Vec<u32> = texto.iter().map(|&letra| letra as u32 + 1).collect();
    simbolos.push(0);
    let mut inicio = Pizarra::nueva()
        .variable("texto", recortar(&texto, 40))
        .variable("n", miles(n + 1));
    if corto {
        let filas = (0..=n)
            .map(|posicion| vec![posicion.to_string(), rotacion(&texto, posicion)])
            .collect();
        inicio = inicio.tabla(
            "rotaciones sin ordenar",
            &["empieza en", "rotación"],
            filas,
            &[],
        );
    } else {
        inicio = inicio.secuencia(&texto, &[]);
    }
    narracion.agregar(inicio.paso(
        "simbolos.push(0);",
        format!(
            "Agregas $ al final del texto: un símbolo único y menor que todas las letras (en el código vale 0 y cada letra vale su byte más 1). Quedan {} rotaciones; la BWT es la última letra de cada una cuando están ordenadas.",
            miles(n + 1)
        ),
    ));
    let sufijos = narrar_sa_is(&mut narracion, &simbolos, &texto, corto);
    let bwt = indice_fm::ultima_columna(&texto, &sufijos);
    narrar_lectura(&mut narracion, &texto, &sufijos, &bwt, corto);
    let recuperado = narrar_inversa(&mut narracion, &bwt, corto);
    let rachas = Rachas::medir(&texto, &bwt);
    let ultimo = Pizarra::nueva()
        .variable("bwt", recortar(&bwt, 60))
        .variable("recuperado", recortar(&recuperado, 40))
        .variable("rachas_bwt", miles(rachas.en_la_bwt))
        .variable("rachas_texto", miles(rachas.en_el_texto))
        .paso(
            "Ok((bwt, recuperado))",
            format!(
                "La BWT es {} y la inversa devolvió {} el texto. La BWT tiene {} rachas de letras iguales contra {} del texto: en «Qué significa» ves si eso ayuda a comprimir.",
                recortar(&bwt, 60),
                if recuperado == texto { "exactamente" } else { "algo distinto de" },
                miles(rachas.en_la_bwt),
                miles(rachas.en_el_texto)
            ),
        );
    Ok(narracion.terminar(ultimo))
}

fn letra_de_simbolo_de_byte(simbolo: u32) -> String {
    if simbolo == 0 {
        "$".to_string()
    } else {
        ((simbolo - 1) as u8 as char).to_string()
    }
}

fn lista_de_sufijos(sufijos: &[u32]) -> Vec<String> {
    sufijos
        .iter()
        .take(MAXIMO_ROTACIONES * 2)
        .map(|&posicion| {
            if posicion == VACIO {
                "—".to_string()
            } else {
                posicion.to_string()
            }
        })
        .collect()
}

fn narrar_sa_is(
    narracion: &mut Narracion,
    simbolos: &[u32],
    texto: &[u8],
    corto: bool,
) -> Vec<u32> {
    let es_s = indice_fm::tipos_de_sufijo(simbolos);
    let lms = indice_fm::posiciones_lms(&es_s);
    let cuantos_s = es_s.iter().filter(|&&tipo| tipo).count();
    let mut pizarra = Pizarra::nueva()
        .variable("tipo_s", miles(cuantos_s))
        .variable("tipo_l", miles(es_s.len() - cuantos_s));
    if corto {
        let filas = (0..simbolos.len())
            .map(|posicion| {
                vec![
                    posicion.to_string(),
                    letra_de_simbolo_de_byte(simbolos[posicion]),
                    if es_s[posicion] { "S" } else { "L" }.to_string(),
                    if indice_fm::es_lms(&es_s, posicion) {
                        "LMS"
                    } else {
                        ""
                    }
                    .to_string(),
                ]
            })
            .collect();
        let resaltar: Vec<(usize, usize)> = lms.iter().map(|&p| (p as usize, 3)).collect();
        pizarra = pizarra.tabla(
            "tipos",
            &["posición", "letra", "tipo", "¿LMS?"],
            filas,
            &resaltar,
        );
    } else {
        let tipos: Vec<u8> = es_s
            .iter()
            .take(240)
            .map(|&tipo| if tipo { b'S' } else { b'L' })
            .collect();
        pizarra = pizarra.secuencia(&tipos, &[]);
    }
    narracion.agregar(pizarra.paso(
        "es_s[i] = texto[i] < texto[i + 1] || (texto[i] == texto[i + 1] && es_s[i + 1]);",
        "Ordenar rotaciones de texto$ es ordenar sus sufijos, y el código lo hace con SA-IS (Nong, Zhang y Chan, 2009). Primero, de derecha a izquierda, marca cada sufijo: S si es menor que el que empieza una letra después, L si es mayor; con letras iguales hereda el tipo del siguiente.",
    ));
    narracion.agregar(
        Pizarra::nueva()
            .variable("lms", miles(lms.len()))
            .lista(
                "posiciones LMS",
                lms.iter().take(MAXIMO_ROTACIONES * 2).map(|p| p.to_string()).collect(),
                &[],
            )
            .paso(
                "lms.push(posicion as u32);",
                format!(
                    "Un LMS es un sufijo S justo después de uno L: hay {} y son las semillas. Si sabes ordenar los LMS, el orden de todos los demás sale solo (se «induce»).",
                    miles(lms.len())
                ),
            ),
    );
    let tamanos = indice_fm::tamanos_de_cubeta(simbolos, LETRAS_DE_UN_BYTE);
    let inicios = indice_fm::inicios_de_cubeta(&tamanos);
    let finales = indice_fm::finales_de_cubeta(&tamanos);
    let filas: Vec<Vec<String>> = (0..LETRAS_DE_UN_BYTE)
        .filter(|&letra| tamanos[letra] > 0)
        .map(|letra| {
            vec![
                letra_de_simbolo_de_byte(letra as u32),
                tamanos[letra].to_string(),
                inicios[letra].to_string(),
                finales[letra].to_string(),
            ]
        })
        .collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("letras_distintas", filas.len())
            .tabla("cubetas", &["letra", "tamaño", "inicio", "fin"], filas, &[])
            .paso(
                "tamanos[letra as usize] += 1;",
                "Cuentas cuántas veces está cada letra: en el arreglo final, los sufijos que empiezan con la misma letra ocupan un tramo seguido (su cubeta), en orden alfabético.",
            ),
    );
    let mut sufijos = vec![VACIO; simbolos.len()];
    indice_fm::colocar_lms(simbolos, &tamanos, &lms, &mut sufijos);
    narracion.agregar(
        Pizarra::nueva()
            .lista("sufijos", lista_de_sufijos(&sufijos), &[])
            .paso(
                "sufijos[finales[letra] as usize] = posicion;",
                "Pones cada LMS al final de la cubeta de su primera letra (— es un lugar vacío). Su orden dentro de la cubeta todavía es provisional.",
            ),
    );
    indice_fm::inducir_tipo_l(simbolos, &es_s, &tamanos, &mut sufijos);
    narracion.agregar(
        Pizarra::nueva()
            .lista("sufijos", lista_de_sufijos(&sufijos), &[])
            .paso(
                "sufijos[inicios[letra] as usize] = anterior as u32;",
                "Recorres el arreglo de izquierda a derecha: si el sufijo que empieza una letra antes de uno ya colocado es L, va al primer lugar libre de su cubeta. Así se colocan todos los L.",
            ),
    );
    indice_fm::inducir_tipo_s(simbolos, &es_s, &tamanos, &mut sufijos);
    narracion.agregar(
        Pizarra::nueva()
            .lista("sufijos", lista_de_sufijos(&sufijos), &[])
            .paso(
                "sufijos[finales[letra] as usize] = anterior as u32;",
                "Ahora de derecha a izquierda, igual con los S, llenando cada cubeta desde el final. Con esto las subcadenas LMS (de un LMS al siguiente) quedan ordenadas.",
            ),
    );
    let ordenados = indice_fm::lms_ordenados(&es_s, &sufijos);
    let (reducido, nombres) = indice_fm::nombrar_lms(simbolos, &es_s, &lms, &ordenados);
    let mut pizarra = Pizarra::nueva()
        .variable("lms", miles(lms.len()))
        .variable("nombres", miles(nombres));
    if corto {
        let filas = lms
            .iter()
            .zip(&reducido)
            .map(|(&posicion, nombre)| {
                vec![
                    posicion.to_string(),
                    subcadena_lms(simbolos, &es_s, posicion as usize),
                    nombre.to_string(),
                ]
            })
            .collect();
        pizarra = pizarra.tabla("nombres", &["LMS en", "subcadena", "nombre"], filas, &[]);
    }
    narracion.agregar(pizarra.paso(
        "nombre += 1;",
        format!(
            "Recorres los LMS en el orden que salió y comparas cada subcadena con la anterior: si cambia, el nombre sube. Hay {} nombres para {} LMS.",
            miles(nombres),
            miles(lms.len())
        ),
    ));
    let lms_en_orden: Vec<u32> = if nombres < lms.len() {
        let orden_reducido = indice_fm::ordenar_sufijos(&reducido, nombres);
        narracion.agregar(
            Pizarra::nueva()
                .variable("reducido", miles(reducido.len()))
                .lista(
                    "texto reducido (nombres)",
                    reducido.iter().take(MAXIMO_ROTACIONES * 2).map(|n| n.to_string()).collect(),
                    &[],
                )
                .paso(
                    "let orden_reducido = ordenar_sufijos(&reducido, nombres);",
                    format!(
                        "Hay nombres repetidos, así que el orden aún no es seguro: el código escribe los nombres como un texto nuevo de {} letras (la mitad o menos) y se llama a sí mismo para ordenarlo. Ese orden es el de los LMS.",
                        miles(reducido.len())
                    ),
                ),
        );
        orden_reducido.iter().map(|&i| lms[i as usize]).collect()
    } else {
        narracion.agregar(
            Pizarra::nueva()
                .variable("nombres", miles(nombres))
                .paso(
                    "let lms_en_orden: Vec<u32> = if nombres < lms.len() {",
                    "Todos los nombres son distintos: el orden de los LMS ya es el definitivo y no hace falta recursión.",
                ),
        );
        ordenados
    };
    indice_fm::colocar_lms(simbolos, &tamanos, &lms_en_orden, &mut sufijos);
    indice_fm::inducir_tipo_l(simbolos, &es_s, &tamanos, &mut sufijos);
    indice_fm::inducir_tipo_s(simbolos, &es_s, &tamanos, &mut sufijos);
    let mostrado = if texto.len() > MAXIMO_ROTACIONES * 2 {
        format!(
            " Ves las primeras {} posiciones de {}.",
            MAXIMO_ROTACIONES * 2,
            miles(sufijos.len())
        )
    } else {
        String::new()
    };
    narracion.agregar(
        Pizarra::nueva()
            .lista("sufijos", lista_de_sufijos(&sufijos), &[])
            .paso(
                "colocar_lms(texto, &tamanos, &lms_en_orden, &mut sufijos);",
                format!(
                    "Con los LMS en su orden final, los colocas otra vez e induces L y S: sale el arreglo de sufijos completo, sin haber comparado nunca dos rotaciones enteras.{mostrado}"
                ),
            ),
    );
    sufijos
}

fn subcadena_lms(simbolos: &[u32], es_s: &[bool], inicio: usize) -> String {
    let mut subcadena = letra_de_simbolo_de_byte(simbolos[inicio]);
    let mut posicion = inicio + 1;
    while posicion < simbolos.len() {
        subcadena.push_str(&letra_de_simbolo_de_byte(simbolos[posicion]));
        if indice_fm::es_lms(es_s, posicion) {
            break;
        }
        posicion += 1;
    }
    subcadena
}

fn narrar_lectura(
    narracion: &mut Narracion,
    texto: &[u8],
    sufijos: &[u32],
    bwt: &[u8],
    corto: bool,
) {
    let mostradas = if corto {
        sufijos.len()
    } else {
        PRIMERAS_FILAS.min(sufijos.len())
    };
    let filas: Vec<Vec<String>> = sufijos
        .iter()
        .take(mostradas)
        .enumerate()
        .map(|(fila, &inicio)| {
            let completa = if corto {
                rotacion(texto, inicio as usize)
            } else {
                sufijo_corto(texto, inicio as usize)
            };
            vec![
                fila.to_string(),
                inicio.to_string(),
                completa,
                (bwt[fila] as char).to_string(),
            ]
        })
        .collect();
    let resaltar: Vec<(usize, usize)> = (0..mostradas).map(|fila| (fila, 3)).collect();
    narracion.agregar(
        Pizarra::nueva()
            .variable("filas", miles(sufijos.len()))
            .tabla(
                "rotaciones ordenadas",
                &["fila", "empieza en", if corto { "rotación" } else { "sufijo" }, "L"],
                filas.clone(),
                &resaltar,
            )
            .paso(
                "for &inicio in sufijos {",
                "Estas son las rotaciones ya ordenadas. La última letra de cada una (columna L) es la letra del texto que va justo antes de donde empieza: eso es lo que el código lee.",
            ),
    );
    let detalladas = if corto {
        sufijos.len()
    } else {
        FILAS_DETALLADAS.min(sufijos.len())
    };
    let mut racha = 0;
    let mut explicada = false;
    for (fila, &inicio) in sufijos.iter().enumerate().take(detalladas) {
        if fila > 0 && bwt[fila] == bwt[fila - 1] {
            racha += 1;
        } else {
            racha = 1;
        }
        let (linea, mut nota) = if inicio == 0 {
            (
                "columna.push(CENTINELA);",
                format!(
                    "Fila {fila}: la rotación empieza en 0, el texto entero; antes de ella solo está el $, así que escribes $."
                ),
            )
        } else {
            (
                "columna.push(texto[inicio as usize - 1]);",
                format!(
                    "Fila {fila}: el sufijo empieza en {inicio}; la letra anterior, texto[{}], es {}.",
                    inicio - 1,
                    texto[inicio as usize - 1] as char
                ),
            )
        };
        if racha >= 2 && !explicada {
            nota.push_str(&format!(
                " Igual que en la fila anterior: las rotaciones vecinas empiezan parecido, así que suelen tener la misma letra antes, y la racha de {} ya mide {racha}. Por eso la BWT junta letras iguales.",
                bwt[fila] as char
            ));
            explicada = true;
        } else if racha >= 2 {
            nota.push_str(&format!(
                " La racha de {} ya mide {racha}.",
                bwt[fila] as char
            ));
        }
        let marcas = if inicio == 0 {
            Vec::new()
        } else {
            vec![(inicio as usize - 1, inicio as usize, Marca::Actual)]
        };
        let mut pizarra = Pizarra::nueva()
            .variable("fila", fila)
            .variable("inicio", inicio)
            .variable("letra", bwt[fila] as char)
            .variable("bwt", recortar(&bwt[..=fila], 60))
            .secuencia(texto, &marcas);
        if corto {
            pizarra = pizarra.tabla(
                "rotaciones ordenadas",
                &["fila", "empieza en", "rotación", "L"],
                filas.clone(),
                &[(fila, 3)],
            );
        }
        narracion.agregar(pizarra.paso(linea, nota));
    }
    if detalladas < sufijos.len() {
        narracion.agregar(
            Pizarra::nueva()
                .variable("bwt", recortar(bwt, 60))
                .paso(
                    "columna.push(texto[inicio as usize - 1]);",
                    format!(
                        "Las {} filas restantes se leen igual. La BWT completa mide {} letras, las mismas del texto más el $.",
                        miles(sufijos.len() - detalladas),
                        miles(bwt.len())
                    ),
                ),
        );
    }
}

fn narrar_inversa(narracion: &mut Narracion, bwt: &[u8], corto: bool) -> Vec<u8> {
    let primera_fila = indice_fm::primeras_filas_de_letras(bwt);
    let rango = indice_fm::rangos(bwt);
    let mut filas = Vec::new();
    for orden in 0..LETRAS_DE_UN_BYTE {
        let cuantas = if orden + 1 < LETRAS_DE_UN_BYTE {
            primera_fila[orden + 1] - primera_fila[orden]
        } else {
            bwt.len() - primera_fila[orden]
        };
        if cuantas > 0 {
            filas.push(vec![
                letra_de_simbolo_de_byte(orden as u32),
                primera_fila[orden].to_string(),
                cuantas.to_string(),
            ]);
        }
    }
    narracion.agregar(
        Pizarra::nueva()
            .variable("bwt", recortar(bwt, 60))
            .tabla("primera_fila", &["letra", "primera fila", "cuántas"], filas, &[])
            .paso(
                "let primera_fila = primeras_filas_de_letras(bwt);",
                "Ahora la deshaces usando solo la BWT. Si la ordenas obtienes la columna F, y para eso basta contar: primera_fila dice en qué fila de F empieza cada letra.",
            ),
    );
    let elementos: Vec<String> = bwt
        .iter()
        .zip(&rango)
        .take(MAXIMO_ROTACIONES * 2)
        .enumerate()
        .map(|(fila, (&letra, &r))| format!("fila {fila}: {} n.º {}", letra as char, r + 1))
        .collect();
    narracion.agregar(
        Pizarra::nueva()
            .lista("rango en L", elementos, &[])
            .paso(
                "rango.push(vistas[orden]);",
                "Para cada fila anotas cuántas letras iguales hubo antes en L. La k-ésima G de L y la k-ésima G de F son la misma G del texto: ese es el mapeo LF.",
            ),
    );
    let mut texto = vec![0u8; bwt.len() - 1];
    let mut fila = 0;
    let detalladas = if corto {
        texto.len()
    } else {
        FILAS_DETALLADAS.min(texto.len())
    };
    for (paso, posicion) in (0..texto.len()).rev().enumerate() {
        let letra = bwt[fila];
        texto[posicion] = letra;
        let siguiente = primera_fila[indice_fm::orden_de_letra(letra)] + rango[fila] as usize;
        if paso < detalladas {
            narracion.agregar(
                Pizarra::nueva()
                    .variable("fila", fila)
                    .variable("letra", letra as char)
                    .variable("posicion", posicion)
                    .variable("recuperado", recortar(&texto[posicion..], 40))
                    .paso(
                        "texto[posicion] = letra;",
                        format!(
                            "La fila {fila} termina en {}: esa letra va justo antes de lo ya reconstruido, en la posición {posicion}. LF te lleva a la fila {siguiente}, donde esa misma {} está en F.",
                            letra as char,
                            letra as char
                        ),
                    ),
            );
        }
        fila = siguiente;
    }
    if detalladas < texto.len() {
        narracion.agregar(
            Pizarra::nueva()
                .variable("recuperado", recortar(&texto, 40))
                .paso(
                    "fila = primera_fila[orden_de_letra(letra)] + rango[fila] as usize;",
                    format!(
                        "Las {} letras restantes salen igual, una por salto LF, de atrás hacia adelante.",
                        miles(texto.len() - detalladas)
                    ),
                ),
        );
    }
    texto
}

#[cfg(test)]
mod pruebas {
    use super::*;
    use nucleo::azar::Azar;

    fn texto_al_azar(azar: &mut Azar, largo: usize) -> Vec<u8> {
        (0..largo).map(|_| b"ACGT"[azar.entero(4)]).collect()
    }

    fn lineas_validas(id: &str, pasos: &[Paso]) {
        let codigo: String = crate::codigo(id)
            .unwrap()
            .iter()
            .map(|fragmento| fragmento.fuente.clone())
            .collect::<Vec<_>>()
            .join("\n");
        for paso in pasos {
            assert!(
                codigo
                    .lines()
                    .any(|linea| linea.contains(paso.linea.as_str())),
                "{id}: falta la línea «{}»",
                paso.linea
            );
        }
        assert!(!pasos.is_empty() && pasos.len() <= MAXIMO_PASOS);
    }

    #[test]
    fn el_recorrido_narrado_es_el_del_nucleo() {
        let mut azar = Azar::nuevo(3);
        for _ in 0..150 {
            let largo = 1 + azar.entero(300);
            let texto = texto_al_azar(&mut azar, largo);
            let indice = IndiceDoble::construir(&texto).unwrap();
            let d = 1 + azar.entero(2);
            let largo_del_patron = d + 1 + azar.entero(10);
            let patron = texto_al_azar(&mut azar, largo_del_patron);
            let cotas = indice.cotas_de_diferencias(&patron);
            let narrados = recorrer_ramas(&indice.directo, &patron, d, &cotas, |_| {});
            let del_nucleo: Vec<(usize, usize, usize)> =
                indice_fm::buscar_en_una_hebra(&indice, &patron, d, b'+')
                    .iter()
                    .map(|i| (i.arriba, i.abajo, i.diferencias))
                    .collect();
            assert_eq!(narrados, del_nucleo);
        }
    }

    #[test]
    fn la_busqueda_exacta_narrada_es_la_del_nucleo() {
        let mut azar = Azar::nuevo(9);
        for _ in 0..100 {
            let largo = 5 + azar.entero(200);
            let texto = texto_al_azar(&mut azar, largo);
            let indice = IndiceDoble::construir(&texto).unwrap();
            let inicio = azar.entero(largo - 4);
            let patron = texto[inicio..inicio + 1 + azar.entero(4)].to_vec();
            let mut narracion = Narracion::nueva(MAXIMO_PASOS);
            let narrados = narrar_busqueda_exacta(&mut narracion, &indice.directo, &texto, &patron);
            let (arriba, abajo) = indice.directo.intervalo(&patron);
            assert_eq!(narrados, vec![(arriba, abajo, 0)]);
        }
    }

    #[test]
    fn sa_is_narrado_es_el_del_nucleo() {
        let mut azar = Azar::nuevo(13);
        for _ in 0..100 {
            let largo = azar.entero(150);
            let texto = texto_al_azar(&mut azar, largo);
            let simbolos: Vec<u32> = texto
                .iter()
                .map(|&l| l as u32 + 1)
                .chain(std::iter::once(0))
                .collect();
            let mut narracion = Narracion::nueva(MAXIMO_PASOS);
            let narrado = narrar_sa_is(&mut narracion, &simbolos, &texto, largo <= 40);
            assert_eq!(narrado, indice_fm::arreglo_de_sufijos(&texto).unwrap());
            let bwt = indice_fm::ultima_columna(&texto, &narrado);
            let recuperado = narrar_inversa(&mut narracion, &bwt, largo <= 40);
            assert_eq!(recuperado, texto);
        }
    }

    #[test]
    fn pasos_del_ejemplo_y_de_otras_entradas() {
        let casos = [
            ("GTGCCAGCAGCCGCGGTAAT GCTGGTGG TTATCCACA", "@ecoli", "1"),
            ("GTGCCAGCAGCCGCGGTAAT", "@ecoli", "0"),
            ("TTATCCACA GATC", "@ecoli", "1"),
            ("TTATCCACA", "@ecoli:3900000-3950000", "2"),
            ("ACAT", "GATTACATACATTT", "0"),
            ("GGGG", "ACATACATAC", "0"),
            ("ATGT ACAT", "ACATACATAC", "1"),
            ("CCGG", "AAACCGTT", "2"),
            (
                "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
                "@ecoli:1-200000",
                "0",
            ),
        ];
        for (patrones, texto, d) in casos {
            let entrada = Entrada::nueva(&[("patrones", patrones), ("texto", texto), ("d", d)]);
            let pasos = pasos_indice_fm(&entrada).unwrap();
            lineas_validas("replicacion.indice_fm", &pasos);
            let salida = buscar_con_indice(&entrada).unwrap();
            assert!(!salida.interpretacion.is_empty());
        }
        let demasiadas = Entrada::nueva(&[("patrones", "GATC"), ("texto", "@ecoli"), ("d", "2")]);
        assert!(pasos_indice_fm(&demasiadas).is_err());
        assert!(buscar_con_indice(&demasiadas).is_err());
    }

    #[test]
    fn pasos_de_la_bwt() {
        for texto in [
            "@ecoli:1198457-1198496",
            "ACATG",
            "A",
            "GATTACAGATTACAGATTACA",
            "@ecoli:1-5000",
        ] {
            let entrada = Entrada::nueva(&[("texto", texto)]);
            let pasos = pasos_bwt(&entrada).unwrap();
            lineas_validas("replicacion.bwt", &pasos);
            let salida = transformada(&entrada).unwrap();
            assert!(!salida.interpretacion.is_empty());
        }
    }

    #[test]
    fn el_ejemplo_encuentra_los_siete_operones() {
        let entrada = Entrada::nueva(&[
            ("patrones", "GTGCCAGCAGCCGCGGTAAT"),
            ("texto", "@ecoli"),
            ("d", "0"),
        ]);
        let pedido = leer_pedido(&entrada).unwrap();
        let usado = preparar_indice(&pedido).unwrap();
        let resultados =
            indice_fm::buscar_patrones(usado.doble(), &pedido.patrones, pedido.d).unwrap();
        assert_eq!(posiciones_de_hebra(&resultados[0], b'+').len(), 5);
        assert_eq!(posiciones_de_hebra(&resultados[0], b'-').len(), 2);
    }

    #[test]
    fn la_poisson_y_el_azar() {
        assert!((poisson_al_menos(2.0, 0) - 1.0).abs() < 1e-12);
        let a_mano = 1.0 - (-2.0f64).exp() * (1.0 + 2.0);
        assert!((poisson_al_menos(2.0, 2) - a_mano).abs() < 1e-9);
        let uniforme = [0.25; 4];
        let exacta = probabilidad_con_diferencias(b"ACGT", 0, &uniforme);
        assert!((exacta - 1.0 / 256.0).abs() < 1e-12);
        let con_una = probabilidad_con_diferencias(b"ACGT", 1, &uniforme);
        assert!((con_una - 13.0 / 256.0).abs() < 1e-12);
    }

    #[test]
    fn el_tiempo_de_medida_se_reparte_entre_muchos_patrones() {
        assert_eq!(tiempo_de_medida(1), TIEMPO_DE_MEDIDA);
        assert_eq!(tiempo_de_medida(250), TIEMPO_DE_MEDIDA);
        assert_eq!(tiempo_de_medida(1000), Duration::from_micros(500));
        assert!(tiempo_de_medida(1000) * 2000 <= TIEMPO_TOTAL_DE_MEDIDA * 2);
        assert_eq!(tiempo_de_medida(0), TIEMPO_DE_MEDIDA);
    }

    #[test]
    fn el_escaneo_clasico_se_corta_al_pasar_su_presupuesto() {
        let mut azar = Azar::nuevo(11);
        let texto = texto_al_azar(&mut azar, 5000);
        let indice = IndiceDoble::construir(&texto).unwrap();
        let patrones: Vec<Vec<u8>> = (0..8).map(|_| texto_al_azar(&mut azar, 12)).collect();
        let resultados = indice_fm::buscar_patrones(&indice, &patrones, 1).unwrap();
        let tiempos: Vec<Tiempos> = patrones
            .iter()
            .map(|_| Tiempos {
                contar: Duration::ZERO,
                ubicar: Duration::ZERO,
            })
            .collect();
        let corto =
            comparar_con_el_escaneo(&texto, &patrones, &resultados, &tiempos, 1, Duration::ZERO);
        assert_eq!(corto.patrones, 1);
        assert!(corto.coinciden);
        let largo = comparar_con_el_escaneo(
            &texto,
            &patrones,
            &resultados,
            &tiempos,
            1,
            Duration::from_secs(600),
        );
        assert_eq!(largo.patrones, 8);
        assert!(largo.coinciden);
    }

    #[test]
    fn los_textos_de_la_comparacion_concuerdan_en_numero() {
        assert_eq!(patrones_comparados(1), "en el único patrón comparado");
        assert_eq!(patrones_comparados(3), "en los 3 patrones comparados");
        assert_eq!(veces_con_decimales(0.3), "0.3 veces");
        assert_eq!(veces_con_decimales(23.4), "23 veces");
    }
}
