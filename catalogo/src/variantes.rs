use crate::formato::{cientifico, decimal, miles, porcentaje};
use crate::pizarra::{Marca, Narracion, Pizarra};
use crate::{Algoritmo, Bloque, Campo, Entrada, Paso, Salida, Tipo};
use nucleo::ecoli::{self, Gen};
use nucleo::formatos::{self, RegistroSam};
use nucleo::secuencia::{aminoacido, indice_base};
use nucleo::variantes::{
    self, Anotacion, Apilamiento, Cambio, Comparacion, Efecto, IndiceDeGenes, Informe, Llamada,
    Umbrales, Vecino,
};
use std::fs::{self, File};
use std::io::{BufRead, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

const PROFUNDIDAD_MAXIMA: usize = 100_000;
const CALIDAD_MAXIMA: usize = 60;
const DISTANCIA_DE_PROMOTOR: usize = 150;
const CUBIERTA_SUFICIENTE: f64 = 0.9;
const PROFUNDIDAD_RECOMENDADA: f64 = 10.0;
const SNPS_PARA_TS_TV: usize = 10;
const QUAL_BAJA: f64 = 20.0;
static CONTADOR: AtomicU64 = AtomicU64::new(0);

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![Algoritmo {
        id: "puerto.llamar_variantes",
        capitulo: "puerto",
        nombre: "Llamado de variantes",
        resumen: "Encuentra SNPs e indels en lecturas alineadas a un genoma haploide y dice qué le hacen a cada gen.",
        campos: &[
            Campo {
                id: "alineamientos",
                etiqueta: "Lecturas alineadas (SAM o SAM.gz, puede estar desordenado)",
                tipo: Tipo::Archivo,
                ejemplo: "@ejemplo:alineamientos",
            },
            Campo {
                id: "referencia",
                etiqueta: "Genoma de referencia (@ecoli, una secuencia o la ruta de un FASTA)",
                tipo: Tipo::Base,
                ejemplo: "@ecoli",
            },
            Campo {
                id: "profundidad",
                etiqueta: "Profundidad mínima (lecturas)",
                tipo: Tipo::Entero,
                ejemplo: "5",
            },
            Campo {
                id: "fraccion",
                etiqueta: "Fracción mínima de lecturas con el alelo alternativo",
                tipo: Tipo::Decimal,
                ejemplo: "0.8",
            },
            Campo {
                id: "calidad",
                etiqueta: "Calidad mínima de base y de mapeo (Phred)",
                tipo: Tipo::Entero,
                ejemplo: "20",
            },
            Campo {
                id: "verdad",
                etiqueta: "Mutaciones verdaderas para comparar (VCF, opcional)",
                tipo: Tipo::Archivo,
                ejemplo: "@ejemplo:mutaciones",
            },
        ],
        funcion: "variantes::llamar_variantes",
        ejecutar: llamar_variantes,
        pasos: Some(pasos_llamado),
    }]
}

struct Conteo {
    transiciones: usize,
    transversiones: usize,
    inserciones: usize,
    deleciones: usize,
}

impl Conteo {
    fn de(llamadas: &[Llamada]) -> Conteo {
        let mut conteo = Conteo {
            transiciones: 0,
            transversiones: 0,
            inserciones: 0,
            deleciones: 0,
        };
        for llamada in llamadas {
            let cambio = &llamada.cambio;
            if cambio.es_snp() {
                if variantes::es_transicion(cambio.referencia[0], cambio.alternativa[0]) {
                    conteo.transiciones += 1;
                } else {
                    conteo.transversiones += 1;
                }
            } else if cambio.es_insercion() {
                conteo.inserciones += 1;
            } else if cambio.es_delecion() {
                conteo.deleciones += 1;
            }
        }
        conteo
    }

    fn snps(&self) -> usize {
        self.transiciones + self.transversiones
    }

    fn indels(&self) -> usize {
        self.inserciones + self.deleciones
    }
}

fn llamar_variantes(entrada: &Entrada) -> Result<Salida, String> {
    let ruta_sam = entrada.ruta("alineamientos")?;
    let referencia = leer_referencia(entrada)?;
    let umbrales = leer_umbrales(entrada)?;
    let es_ecoli = referencia.as_slice() == ecoli::genoma();
    let genes: &[Gen] = if es_ecoli { ecoli::genes() } else { &[] };
    let informe = variantes::llamar_variantes(&ruta_sam, &referencia, genes, &umbrales)?;
    let ruta_vcf = guardar_vcf(&ruta_sam, &informe, referencia.len())?;
    let comparacion = if entrada.tiene("verdad") {
        Some(comparar_con_verdad(entrada, &informe, &referencia)?)
    } else {
        None
    };
    let conteo = Conteo::de(&informe.llamadas);

    let mut salida = Salida::default();
    salida.agregar(Bloque::valor("SNPs", miles(conteo.snps())));
    salida.agregar(Bloque::valor("Inserciones", miles(conteo.inserciones)));
    salida.agregar(Bloque::valor("Deleciones", miles(conteo.deleciones)));
    if conteo.transversiones > 0 {
        salida.agregar(Bloque::valor(
            "Transiciones / transversiones (Ts/Tv)",
            decimal(cociente_ts_tv(&conteo), 2),
        ));
    }
    salida.agregar(Bloque::valor(
        "Lecturas usadas",
        format!(
            "{} de {}",
            miles(informe.lecturas_usadas),
            miles(informe.lecturas_usadas + informe.lecturas_descartadas)
        ),
    ));
    salida.agregar(Bloque::valor(
        &format!(
            "Referencia con al menos {} lecturas",
            umbrales.profundidad_minima
        ),
        porcentaje(fraccion_cubierta(&informe, referencia.len()), 1),
    ));
    salida.agregar(Bloque::valor(
        "Profundidad media donde hay lecturas",
        format!("{}×", decimal(informe.profundidad.media, 1)),
    ));
    salida.agregar(Bloque::tabla(
        "Variantes",
        &[
            "Posición (desde 1)",
            "REF",
            "ALT",
            "Profundidad",
            "Fracción",
            "QUAL",
            "Gen",
            "Efecto",
            "Cambio en la proteína",
            "Nota",
        ],
        informe.llamadas.iter().map(fila_de_variante).collect(),
    ));
    salida.agregar(Bloque::barras(
        "Tipos de variante",
        "Tipo",
        "Cantidad",
        vec![
            "Transiciones".to_string(),
            "Transversiones".to_string(),
            "Inserciones".to_string(),
            "Deleciones".to_string(),
        ],
        vec![
            conteo.transiciones as f64,
            conteo.transversiones as f64,
            conteo.inserciones as f64,
            conteo.deleciones as f64,
        ],
    ));
    if es_ecoli && !informe.llamadas.is_empty() {
        let (categorias, valores) = contar_efectos(&informe.llamadas);
        salida.agregar(Bloque::barras(
            "Efecto sobre los genes",
            "Efecto",
            "Variantes",
            categorias,
            valores,
        ));
    }
    salida.agregar(Bloque::archivo("Variantes (VCF)", &ruta_vcf));
    if let Some(comparacion) = &comparacion {
        agregar_comparacion(&mut salida, comparacion);
    }

    explicar_resultado(&mut salida, &informe, &umbrales, &conteo, es_ecoli);
    explicar_cobertura(&mut salida, &informe, &umbrales, referencia.len());
    explicar_resistencias(&mut salida, &informe.llamadas);
    explicar_inactivaciones(&mut salida, &informe.llamadas);
    if let Some(comparacion) = &comparacion {
        explicar_comparacion(&mut salida, comparacion);
    }
    Ok(salida)
}

fn leer_referencia(entrada: &Entrada) -> Result<Vec<u8>, String> {
    let texto = entrada.texto("referencia")?;
    let limpio = texto.trim();
    let ruta = Path::new(limpio);
    if !limpio.starts_with('@') && !limpio.contains('\n') && ruta.is_file() {
        return formatos::leer_fasta_archivo(ruta)?
            .into_iter()
            .map(|registro| registro.secuencia)
            .find(|secuencia| !secuencia.is_empty())
            .ok_or_else(|| format!("{} no tiene ninguna secuencia.", ruta.display()));
    }
    entrada.secuencia("referencia")
}

fn leer_umbrales(entrada: &Entrada) -> Result<Umbrales, String> {
    let profundidad = entrada.natural("profundidad")?;
    if !(1..=PROFUNDIDAD_MAXIMA).contains(&profundidad) {
        return Err(format!(
            "La profundidad mínima debe estar entre 1 y {} lecturas.",
            miles(PROFUNDIDAD_MAXIMA)
        ));
    }
    let fraccion = entrada.decimal("fraccion")?;
    if !(0.5..=1.0).contains(&fraccion) {
        return Err(
            "En un genoma haploide la fracción mínima debe estar entre 0.5 y 1: el alelo alternativo tiene que estar en la mayoría de las lecturas."
                .to_string(),
        );
    }
    let calidad = entrada.natural("calidad")?;
    if calidad > CALIDAD_MAXIMA {
        return Err(format!(
            "La calidad mínima debe estar entre 0 y {CALIDAD_MAXIMA} (escala Phred)."
        ));
    }
    Ok(Umbrales {
        calidad_mapeo_minima: calidad as u8,
        calidad_base_minima: calidad as u8,
        profundidad_minima: profundidad as u32,
        fraccion_minima: fraccion,
        ..Umbrales::default()
    })
}

fn guardar_vcf(ruta_sam: &Path, informe: &Informe, longitud: usize) -> Result<PathBuf, String> {
    let junto_al_sam = crate::ruta_de_salida(ruta_sam, "variantes", "vcf");
    if escribir_vcf(&junto_al_sam, informe, longitud).is_ok() {
        return Ok(junto_al_sam);
    }
    let nombre = junto_al_sam
        .file_name()
        .map(|nombre| nombre.to_os_string())
        .unwrap_or_else(|| "variantes.vcf".into());
    let en_la_carpeta_de_trabajo = crate::carpeta_de_ejemplos()?.join(nombre);
    escribir_vcf(&en_la_carpeta_de_trabajo, informe, longitud)?;
    Ok(en_la_carpeta_de_trabajo)
}

fn escribir_vcf(ruta: &Path, informe: &Informe, longitud: usize) -> Result<(), String> {
    let numero = CONTADOR.fetch_add(1, Ordering::Relaxed);
    let temporal = ruta.with_extension(format!("tmp{}-{numero}", std::process::id()));
    let archivo = File::create(&temporal)
        .map_err(|error| format!("No se pudo crear {}: {error}", temporal.display()))?;
    let mut destino = BufWriter::new(archivo);
    let escrito = variantes::escribir_llamadas_vcf(
        &mut destino,
        &informe.nombre_referencia,
        longitud,
        &informe.llamadas,
    )
    .and_then(|_| destino.flush());
    drop(destino);
    if let Err(error) = escrito {
        let _ = fs::remove_file(&temporal);
        return Err(format!("No se pudo escribir {}: {error}", ruta.display()));
    }
    fs::rename(&temporal, ruta)
        .map_err(|error| format!("No se pudo guardar {}: {error}", ruta.display()))
}

fn comparar_con_verdad(
    entrada: &Entrada,
    informe: &Informe,
    referencia: &[u8],
) -> Result<Comparacion, String> {
    let ruta = entrada.ruta("verdad")?;
    let verdad = variantes::leer_cambios_vcf(formatos::abrir(&ruta)?)?;
    let llamadas: Vec<Cambio> = informe
        .llamadas
        .iter()
        .map(|llamada| llamada.cambio.clone())
        .collect();
    Ok(variantes::comparar_con_verdad(
        &llamadas, &verdad, referencia,
    ))
}

fn cociente_ts_tv(conteo: &Conteo) -> f64 {
    conteo.transiciones as f64 / conteo.transversiones.max(1) as f64
}

fn fraccion_cubierta(informe: &Informe, longitud: usize) -> f64 {
    informe.profundidad.suficientes as f64 / longitud.max(1) as f64
}

fn alelo(bases: &[u8]) -> String {
    if bases.len() <= 15 {
        return String::from_utf8_lossy(bases).to_string();
    }
    format!(
        "{}… ({} pb)",
        String::from_utf8_lossy(&bases[..12]),
        bases.len()
    )
}

fn nombre_del_efecto(efecto: &Efecto) -> &'static str {
    match efecto {
        Efecto::Intergenica { .. } => "intergénica",
        Efecto::NoCodificante => "en un gen sin proteína",
        Efecto::Sinonima => "sinónima",
        Efecto::CambioDeAminoacido => "cambio de aminoácido",
        Efecto::ParadaPrematura => "parada prematura",
        Efecto::PerdidaDeParada => "pérdida del codón de parada",
        Efecto::PerdidaDeInicio => "pérdida del codón de inicio",
        Efecto::CorrimientoDeMarco => "corrimiento del marco de lectura",
        Efecto::InsercionEnMarco => "inserción sin corrimiento",
        Efecto::DelecionEnMarco => "deleción sin corrimiento",
    }
}

fn corta_la_proteina(efecto: &Efecto) -> bool {
    matches!(
        efecto,
        Efecto::ParadaPrematura | Efecto::CorrimientoDeMarco | Efecto::PerdidaDeInicio
    )
}

fn nombre_del_lugar(anotacion: &Anotacion) -> String {
    let Efecto::Intergenica { izquierda, derecha } = &anotacion.efecto else {
        return anotacion.nombre_gen.clone();
    };
    match (izquierda, derecha) {
        (Some(antes), Some(despues)) => format!("entre {} y {}", antes.nombre, despues.nombre),
        (Some(antes), None) => format!("después de {}", antes.nombre),
        (None, Some(despues)) => format!("antes de {}", despues.nombre),
        (None, None) => "—".to_string(),
    }
}

fn promotor_cercano(anotacion: &Anotacion) -> Option<&Vecino> {
    let Efecto::Intergenica { izquierda, derecha } = &anotacion.efecto else {
        return None;
    };
    [izquierda, derecha]
        .into_iter()
        .flatten()
        .filter(|vecino| vecino.corriente_arriba && vecino.distancia <= DISTANCIA_DE_PROMOTOR)
        .min_by_key(|vecino| vecino.distancia)
}

fn nota_de_anotacion(anotacion: &Anotacion) -> Option<String> {
    if let Some(antibiotico) = anotacion.resistencia {
        return Some(format!("resistencia a {antibiotico}"));
    }
    if corta_la_proteina(&anotacion.efecto) {
        return Some(format!("probablemente inactiva {}", anotacion.nombre_gen));
    }
    promotor_cercano(anotacion).map(|vecino| {
        format!(
            "{} pb antes del inicio de {}: posible promotor",
            vecino.distancia, vecino.nombre
        )
    })
}

fn unir(textos: Vec<String>) -> String {
    let utiles: Vec<String> = textos
        .into_iter()
        .filter(|texto| !texto.is_empty())
        .collect();
    if utiles.is_empty() {
        return "—".to_string();
    }
    utiles.join(" / ")
}

fn fila_de_variante(llamada: &Llamada) -> Vec<String> {
    let anotaciones = &llamada.anotaciones;
    let mut notas: Vec<String> = anotaciones.iter().filter_map(nota_de_anotacion).collect();
    if llamada.calidad < QUAL_BAJA {
        notas.push("QUAL baja: podría ser un error".to_string());
    }
    vec![
        miles(llamada.cambio.posicion + 1),
        alelo(&llamada.cambio.referencia),
        alelo(&llamada.cambio.alternativa),
        miles(llamada.profundidad as usize),
        porcentaje(llamada.fraccion, 0),
        decimal(llamada.calidad, 0),
        unir(anotaciones.iter().map(nombre_del_lugar).collect()),
        unir(
            anotaciones
                .iter()
                .map(|anotacion| nombre_del_efecto(&anotacion.efecto).to_string())
                .collect(),
        ),
        unir(
            anotaciones
                .iter()
                .map(|anotacion| anotacion.proteina.clone())
                .collect(),
        ),
        if notas.is_empty() {
            "—".to_string()
        } else {
            notas.join("; ")
        },
    ]
}

fn contar_efectos(llamadas: &[Llamada]) -> (Vec<String>, Vec<f64>) {
    let mut categorias: Vec<String> = Vec::new();
    let mut valores: Vec<f64> = Vec::new();
    for anotacion in llamadas.iter().flat_map(|llamada| &llamada.anotaciones) {
        let nombre = nombre_del_efecto(&anotacion.efecto);
        match categorias.iter().position(|categoria| categoria == nombre) {
            Some(indice) => valores[indice] += 1.0,
            None => {
                categorias.push(nombre.to_string());
                valores.push(1.0);
            }
        }
    }
    (categorias, valores)
}

fn agregar_comparacion(salida: &mut Salida, comparacion: &Comparacion) {
    if let Some(sensibilidad) = comparacion.sensibilidad() {
        salida.agregar(Bloque::valor("Sensibilidad", porcentaje(sensibilidad, 1)));
    }
    if let Some(precision) = comparacion.precision() {
        salida.agregar(Bloque::valor("Precisión", porcentaje(precision, 1)));
    }
    let mut filas: Vec<(Cambio, &str)> = Vec::new();
    for cambio in &comparacion.verdaderos_positivos {
        filas.push((cambio.clone(), "verdadero positivo"));
    }
    for cambio in &comparacion.falsos_positivos {
        filas.push((cambio.clone(), "falso positivo"));
    }
    for cambio in &comparacion.falsos_negativos {
        filas.push((cambio.clone(), "falso negativo (no encontrada)"));
    }
    filas.sort();
    salida.agregar(Bloque::tabla(
        "Comparación con las mutaciones verdaderas",
        &["Posición (desde 1)", "REF", "ALT", "Resultado"],
        filas
            .into_iter()
            .map(|(cambio, resultado)| {
                vec![
                    miles(cambio.posicion + 1),
                    alelo(&cambio.referencia),
                    alelo(&cambio.alternativa),
                    resultado.to_string(),
                ]
            })
            .collect(),
    ));
}

fn plural(cantidad: usize, singular: &str, varios: &str) -> String {
    if cantidad == 1 {
        format!("1 {singular}")
    } else {
        format!("{} {varios}", miles(cantidad))
    }
}

fn explicar_resultado(
    salida: &mut Salida,
    informe: &Informe,
    umbrales: &Umbrales,
    conteo: &Conteo,
    es_ecoli: bool,
) {
    if informe.lecturas_usadas == 0 {
        salida.explicar(format!(
            "Ninguna lectura pasó los filtros: se usan solo lecturas mapeadas a esta referencia, primarias, no duplicadas y con calidad de mapeo de al menos Q{}. Revisa que el SAM venga de mapear contra esta misma referencia; sin lecturas no se puede decir nada sobre variantes.",
            umbrales.calidad_mapeo_minima
        ));
        return;
    }
    let mut texto = if informe.llamadas.is_empty() {
        format!(
            "No apareció ninguna variante con al menos {} lecturas: donde hay cobertura, la muestra es igual a la referencia.",
            umbrales.profundidad_minima
        )
    } else {
        format!(
            "Encontraste {} y {} ({} y {}). Un SNP necesita al menos {} lecturas y {} de ellas con la otra base; a un indel le basta la mitad, porque las lecturas con indels se alinean peor.",
            plural(conteo.snps(), "SNP", "SNPs"),
            plural(conteo.indels(), "indel", "indels"),
            plural(conteo.inserciones, "inserción", "inserciones"),
            plural(conteo.deleciones, "deleción", "deleciones"),
            umbrales.profundidad_minima,
            porcentaje(umbrales.fraccion_minima, 0)
        )
    };
    if !es_ecoli {
        texto.push_str(
            " No se anotaron genes porque la referencia no es el genoma de E. coli incluido.",
        );
    }
    salida.explicar(texto);
    if conteo.snps() >= SNPS_PARA_TS_TV {
        let cociente = cociente_ts_tv(conteo);
        let valor = if conteo.transversiones == 0 {
            "infinito (no hay transversiones)".to_string()
        } else {
            decimal(cociente, 2)
        };
        let lectura = if conteo.transversiones > 0 && cociente < 0.8 {
            "Si los cambios cayeran al azar saldría cerca de 0.5, porque cada base tiene dos transversiones posibles y una sola transición. Un valor así sugiere mutaciones al azar (como las de una simulación) o errores, no la huella típica de la evolución bacteriana."
        } else if conteo.transversiones == 0 || cociente > 1.5 {
            "Hay más transiciones (A↔G, C↔T) de las que daría el azar (≈ 0.5): es la huella típica de mutaciones reales en bacterias, porque las transiciones ocurren con más facilidad."
        } else {
            "Está entre el 0.5 del azar y los valores mayores que 1.5 de las mutaciones reales en bacterias: puede ser una mezcla de mutaciones reales y errores."
        };
        salida.explicar(format!("Ts/Tv = {valor}. {lectura}"));
    }
}

fn explicar_cobertura(
    salida: &mut Salida,
    informe: &Informe,
    umbrales: &Umbrales,
    longitud: usize,
) {
    if informe.lecturas_usadas == 0 {
        return;
    }
    let mut partes: Vec<String> = Vec::new();
    let cubierta = fraccion_cubierta(informe, longitud);
    if cubierta < CUBIERTA_SUFICIENTE {
        partes.push(format!(
            "Solo el {} de la referencia tiene al menos {} lecturas: en el resto no puedes afirmar que no haya variantes, solo que no las viste.",
            porcentaje(cubierta, 1),
            umbrales.profundidad_minima
        ));
    }
    if informe.profundidad.media < PROFUNDIDAD_RECOMENDADA {
        partes.push(format!(
            "La profundidad media es {}×: con menos de 10× muchas variantes reales no juntan lecturas suficientes para confirmarlas.",
            decimal(informe.profundidad.media, 1)
        ));
    }
    if informe.profundidad.sin_confirmar > 0 {
        partes.push(format!(
            "Hay {} donde casi todas las lecturas muestran otra base, pero son menos de {}: son candidatas que habría que confirmar con más cobertura.",
            plural(
                informe.profundidad.sin_confirmar,
                "posición",
                "posiciones"
            ),
            umbrales.profundidad_minima
        ));
    }
    if !partes.is_empty() {
        salida.explicar(partes.join(" "));
    }
}

fn explicar_resistencias(salida: &mut Salida, llamadas: &[Llamada]) {
    let mut frases: Vec<String> = Vec::new();
    for anotacion in llamadas.iter().flat_map(|llamada| &llamada.anotaciones) {
        let Some(antibiotico) = anotacion.resistencia else {
            continue;
        };
        let sitio = match anotacion.nombre_gen.as_str() {
            "rpoB" => {
                "la región de rpoB que determina la resistencia (codones 507–533), justo donde la rifampicina se une a la ARN polimerasa"
            }
            "gyrA" => {
                "un aminoácido clave de la girasa (codones 83 y 87 de gyrA), donde se unen las quinolonas"
            }
            _ => {
                "un aminoácido clave de la topoisomerasa IV (codones 80 y 84 de parC), donde se unen las quinolonas"
            }
        };
        frases.push(format!(
            "{} {} cae en {sitio}: esta cepa probablemente es resistente a {}.",
            anotacion.nombre_gen,
            anotacion.proteina,
            if antibiotico == "rifampicina" {
                "la rifampicina"
            } else {
                "las quinolonas (como la ciprofloxacina)"
            }
        ));
    }
    if frases.is_empty() {
        return;
    }
    frases.push(
        "Es una predicción a partir del genoma: un antibiograma en el laboratorio lo confirma."
            .to_string(),
    );
    salida.explicar(frases.join(" "));
}

fn explicar_inactivaciones(salida: &mut Salida, llamadas: &[Llamada]) {
    let mut genes: Vec<String> = Vec::new();
    for anotacion in llamadas.iter().flat_map(|llamada| &llamada.anotaciones) {
        if corta_la_proteina(&anotacion.efecto) && !genes.contains(&anotacion.nombre_gen) {
            genes.push(anotacion.nombre_gen.clone());
        }
    }
    if genes.is_empty() {
        return;
    }
    let cantidad = genes.len();
    let mut nombres = genes.iter().take(5).cloned().collect::<Vec<_>>().join(", ");
    if cantidad > 5 {
        nombres.push_str(&format!(" y {} más", cantidad - 5));
    }
    salida.explicar(format!(
        "En {nombres} hay un codón de parada prematuro, un corrimiento del marco de lectura o se perdió el codón de inicio: la proteína sale cortada o no se fabrica, así que {} probablemente quedó inactivado. Cuanto más cerca del inicio del gen esté el cambio, más seguro es.",
        if cantidad == 1 {
            "ese gen"
        } else {
            "cada uno de esos genes"
        }
    ));
}

fn explicar_comparacion(salida: &mut Salida, comparacion: &Comparacion) {
    let verdaderos = comparacion.verdaderos_positivos.len();
    let falsos_positivos = comparacion.falsos_positivos.len();
    let falsos_negativos = comparacion.falsos_negativos.len();
    let (Some(sensibilidad), Some(precision)) =
        (comparacion.sensibilidad(), comparacion.precision())
    else {
        salida.explicar(
            "No se pudo calcular sensibilidad y precisión: el VCF de verdad o el resultado no tienen variantes.",
        );
        return;
    };
    let mut texto = format!(
        "Contra las mutaciones verdaderas: sensibilidad {} (encontraste {} de {}) y precisión {} ({} de tus {} llamadas son reales).",
        porcentaje(sensibilidad, 1),
        verdaderos,
        verdaderos + falsos_negativos,
        porcentaje(precision, 1),
        verdaderos,
        verdaderos + falsos_positivos
    );
    if falsos_negativos == 0 && falsos_positivos == 0 {
        texto.push_str(" No se perdió ninguna mutación y ninguna llamada es falsa.");
    }
    if falsos_negativos > 0 {
        let indels = comparacion
            .falsos_negativos
            .iter()
            .filter(|cambio| !cambio.es_snp())
            .count();
        let cuales = match indels {
            0 => "ninguna es un indel".to_string(),
            1 => "1 es un indel".to_string(),
            varios => format!("{varios} son indels"),
        };
        texto.push_str(&format!(
            " {} {}: {cuales}. Los indels cuestan más porque las lecturas que los cruzan se alinean peor, y cualquier mutación se pierde si cae donde hay pocas lecturas.",
            if falsos_negativos == 1 {
                "Se perdió"
            } else {
                "Se perdieron"
            },
            plural(falsos_negativos, "mutación", "mutaciones")
        ));
    }
    if falsos_positivos > 0 {
        texto.push_str(
            " Los falsos positivos suelen venir de errores de secuenciación repetidos o de lecturas mal mapeadas en zonas repetidas del genoma.",
        );
    }
    salida.explicar(texto);
}

const LINEAS_PARA_BUSCAR_EJEMPLOS: usize = 200_000;
const EJEMPLOS_DE_ANOTACION: usize = 3;
const BASES_EN_LA_TABLA: usize = 6;

struct EjemplosDeLecturas {
    usada: Option<RegistroSam>,
    descartada: Option<(RegistroSam, String)>,
    con_indel: Option<RegistroSam>,
}

struct Escaneo {
    con_lecturas: usize,
    pocas: usize,
    no_acgt: usize,
    coinciden: usize,
    ruido: usize,
    ejemplo_pocas: Option<usize>,
    ejemplo_sin_lecturas: Option<usize>,
    ejemplo_coinciden: Option<usize>,
    ejemplo_ruido: Option<usize>,
}

struct CandidatoIndel {
    cambio: Cambio,
    a_favor: u32,
    cruzan: u32,
    profundidad: u32,
    fraccion: f64,
    pasa: bool,
}

fn corto(texto: &str, maximo: usize) -> String {
    if texto.chars().count() <= maximo {
        return texto.to_string();
    }
    let inicio: String = texto.chars().take(maximo).collect();
    format!("{inicio}…")
}

fn cambio_en_texto(cambio: &Cambio) -> String {
    format!(
        "{} → {}",
        alelo(&cambio.referencia),
        alelo(&cambio.alternativa)
    )
}

fn lugar(posicion: usize) -> String {
    format!("posición {} (contando desde 1)", miles(posicion + 1))
}

fn motivo_de_descarte(
    lectura: &RegistroSam,
    nombre_referencia: &str,
    longitud: usize,
    umbrales: &Umbrales,
) -> String {
    if !lectura.mapeado() {
        "no se pudo mapear en ningún lugar del genoma".to_string()
    } else if lectura.posicion > longitud {
        "dice que empieza fuera de la referencia".to_string()
    } else if lectura.bandera & 0xF00 != 0 {
        "está marcada como secundaria, suplementaria, duplicada o que no pasó el control de calidad"
            .to_string()
    } else if lectura.calidad_mapeo < umbrales.calidad_mapeo_minima {
        format!(
            "su calidad de mapeo es {} y el mínimo es {}: no hay seguridad de que esté en el lugar correcto",
            lectura.calidad_mapeo, umbrales.calidad_mapeo_minima
        )
    } else if lectura.referencia != nombre_referencia {
        "viene de otra referencia distinta a la que elegiste".to_string()
    } else {
        "no trae secuencia o CIGAR".to_string()
    }
}

fn tiene_indel(lectura: &RegistroSam) -> bool {
    lectura
        .cigar
        .iter()
        .any(|&(_, operacion)| operacion == b'I' || operacion == b'D')
}

fn buscar_ejemplos(
    ruta_sam: &Path,
    nombre_referencia: &str,
    longitud: usize,
    umbrales: &Umbrales,
) -> Result<EjemplosDeLecturas, String> {
    let mut ejemplos = EjemplosDeLecturas {
        usada: None,
        descartada: None,
        con_indel: None,
    };
    let fuente = formatos::abrir(ruta_sam)?;
    for (numero, linea) in fuente.split(b'\n').enumerate() {
        if numero >= LINEAS_PARA_BUSCAR_EJEMPLOS {
            break;
        }
        let linea = linea.map_err(|error| format!("No se pudo leer el SAM: {error}"))?;
        let Ok(Some(lectura)) = formatos::leer_sam(&linea) else {
            continue;
        };
        if variantes::sirve(&lectura, nombre_referencia, longitud, umbrales) {
            if ejemplos.con_indel.is_none() && tiene_indel(&lectura) {
                ejemplos.con_indel = Some(lectura.clone());
            }
            if ejemplos.usada.is_none() {
                ejemplos.usada = Some(lectura);
            }
        } else if ejemplos.descartada.is_none() {
            let motivo = motivo_de_descarte(&lectura, nombre_referencia, longitud, umbrales);
            ejemplos.descartada = Some((lectura, motivo));
        }
        if ejemplos.usada.is_some() && ejemplos.descartada.is_some() && ejemplos.con_indel.is_some()
        {
            break;
        }
    }
    Ok(ejemplos)
}

fn tramo_en_referencia(lectura: &RegistroSam) -> usize {
    lectura
        .cigar
        .iter()
        .filter(|&&(_, operacion)| matches!(operacion, b'M' | b'=' | b'X' | b'D' | b'N'))
        .map(|&(largo, _)| largo as usize)
        .sum()
}

fn cigar_en_lista(lectura: &RegistroSam) -> Vec<String> {
    lectura
        .cigar
        .iter()
        .map(|&(largo, operacion)| format!("{largo}{}", operacion as char))
        .collect()
}

fn primer_tramo_alineado(lectura: &RegistroSam) -> Option<(usize, usize, usize)> {
    let mut en_referencia = lectura.posicion.checked_sub(1)?;
    let mut en_lectura = 0;
    for &(largo, operacion) in &lectura.cigar {
        let largo = largo as usize;
        match operacion {
            b'M' | b'=' | b'X' => return Some((en_referencia, en_lectura, largo)),
            b'I' | b'S' => en_lectura += largo,
            b'D' | b'N' => en_referencia += largo,
            _ => {}
        }
    }
    None
}

fn bases_que_cuentan(lectura: &RegistroSam, umbrales: &Umbrales) -> usize {
    let mut cuentan = 0;
    let mut en_lectura = 0;
    for &(largo, operacion) in &lectura.cigar {
        let largo = largo as usize;
        match operacion {
            b'M' | b'=' | b'X' => {
                for i in en_lectura..en_lectura + largo {
                    let calidad = calidad_de(lectura, i);
                    let base = lectura.secuencia.get(i).copied().unwrap_or(b'N');
                    if calidad >= umbrales.calidad_base_minima && indice_base(base).is_some() {
                        cuentan += 1;
                    }
                }
                en_lectura += largo;
            }
            b'I' | b'S' => en_lectura += largo,
            _ => {}
        }
    }
    cuentan
}

fn calidad_de(lectura: &RegistroSam, indice: usize) -> u8 {
    match lectura.calidad.get(indice) {
        Some(&simbolo) => simbolo.saturating_sub(33),
        None => 30,
    }
}

struct PrimerIndel {
    operacion: u8,
    largo: usize,
    en_referencia: usize,
    en_lectura: usize,
}

fn primer_indel(lectura: &RegistroSam) -> Option<PrimerIndel> {
    let inicio = lectura.posicion.checked_sub(1)?;
    let mut en_referencia = inicio;
    let mut en_lectura = 0;
    for &(largo, operacion) in &lectura.cigar {
        let largo = largo as usize;
        match operacion {
            b'M' | b'=' | b'X' => {
                en_referencia += largo;
                en_lectura += largo;
            }
            b'I' => {
                if en_referencia > inicio {
                    return Some(PrimerIndel {
                        operacion,
                        largo,
                        en_referencia,
                        en_lectura,
                    });
                }
                en_lectura += largo;
            }
            b'D' => {
                if en_referencia > inicio {
                    return Some(PrimerIndel {
                        operacion,
                        largo,
                        en_referencia,
                        en_lectura,
                    });
                }
                en_referencia += largo;
            }
            b'N' => en_referencia += largo,
            b'S' => en_lectura += largo,
            _ => {}
        }
    }
    None
}

fn cambio_del_indel(
    lectura: &RegistroSam,
    indel: &PrimerIndel,
    referencia: &[u8],
) -> Option<Cambio> {
    let ancla = indel.en_referencia - 1;
    let base_ancla = *referencia.get(ancla)?;
    if indel.operacion == b'I' {
        let insertadas = lectura
            .secuencia
            .get(indel.en_lectura..indel.en_lectura + indel.largo)?;
        let mut alternativa = vec![base_ancla];
        alternativa.extend_from_slice(insertadas);
        return Some(Cambio {
            posicion: ancla,
            referencia: vec![base_ancla],
            alternativa,
        });
    }
    let fin = ancla + indel.largo + 1;
    Some(Cambio {
        posicion: ancla,
        referencia: referencia.get(ancla..fin)?.to_vec(),
        alternativa: vec![base_ancla],
    })
}

fn escanear(apilamiento: &Apilamiento, referencia: &[u8], umbrales: &Umbrales) -> Escaneo {
    let mut escaneo = Escaneo {
        con_lecturas: 0,
        pocas: 0,
        no_acgt: 0,
        coinciden: 0,
        ruido: 0,
        ejemplo_pocas: None,
        ejemplo_sin_lecturas: None,
        ejemplo_coinciden: None,
        ejemplo_ruido: None,
    };
    for (posicion, &base_referencia) in referencia.iter().enumerate().take(apilamiento.bases.len())
    {
        let profundidad = apilamiento.profundidad(posicion);
        if profundidad > 0 {
            escaneo.con_lecturas += 1;
        }
        if profundidad < umbrales.profundidad_minima {
            escaneo.pocas += 1;
            if profundidad > 0 && escaneo.ejemplo_pocas.is_none() {
                escaneo.ejemplo_pocas = Some(posicion);
            }
            if escaneo.ejemplo_sin_lecturas.is_none() {
                escaneo.ejemplo_sin_lecturas = Some(posicion);
            }
            continue;
        }
        let Some(indice_referencia) = indice_base(base_referencia) else {
            escaneo.no_acgt += 1;
            continue;
        };
        let (_, a_favor) =
            variantes::alternativa_mas_frecuente(&apilamiento.bases[posicion], indice_referencia);
        if a_favor == 0 {
            escaneo.coinciden += 1;
            if escaneo.ejemplo_coinciden.is_none() {
                escaneo.ejemplo_coinciden = Some(posicion);
            }
            continue;
        }
        if (a_favor as f64 / profundidad as f64) < umbrales.fraccion_minima {
            escaneo.ruido += 1;
            if escaneo.ejemplo_ruido.is_none() {
                escaneo.ejemplo_ruido = Some(posicion);
            }
        }
    }
    escaneo
}

fn candidatos_de_indel(apilamiento: &Apilamiento, umbrales: &Umbrales) -> Vec<CandidatoIndel> {
    let mut candidatos: Vec<CandidatoIndel> = apilamiento
        .indels
        .iter()
        .filter_map(|(cambio, &a_favor)| {
            let &cruzan = apilamiento.cruzan.get(cambio.posicion)?;
            let profundidad = cruzan.max(a_favor);
            let fraccion = a_favor as f64 / profundidad as f64;
            let pasa = profundidad >= umbrales.profundidad_minima
                && fraccion >= umbrales.fraccion_minima_indel;
            Some(CandidatoIndel {
                cambio: cambio.clone(),
                a_favor,
                cruzan,
                profundidad,
                fraccion,
                pasa,
            })
        })
        .collect();
    candidatos.sort_by(|una, otra| una.cambio.cmp(&otra.cambio));
    candidatos
}

fn tabla_de_pila(
    pizarra: Pizarra,
    apilamiento: &Apilamiento,
    posicion: usize,
    resaltar: &[(usize, usize)],
) -> Pizarra {
    let conteos = apilamiento.bases[posicion];
    pizarra.tabla(
        "Lecturas apiladas en esta posición",
        &["A", "C", "G", "T"],
        vec![conteos.iter().map(|n| miles(*n as usize)).collect()],
        resaltar,
    )
}

fn gen_que_contiene<'a>(genes: &'a [Gen], nombre: &str, posicion: usize) -> Option<&'a Gen> {
    genes.iter().find(|gen_buscado| {
        gen_buscado.nombre == nombre
            && gen_buscado.inicio - 1 <= posicion
            && posicion < gen_buscado.fin
    })
}

fn posiciones_del_codon(elegido: &Gen, numero: usize) -> (usize, usize) {
    if elegido.hebra == b'-' {
        let ultima = elegido.fin - 1 - (numero - 1) * 3;
        (ultima - 2, ultima + 1)
    } else {
        let primera = elegido.inicio - 1 + (numero - 1) * 3;
        (primera, primera + 3)
    }
}

fn pasos_llamado(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let ruta_sam = entrada.ruta("alineamientos")?;
    let referencia = leer_referencia(entrada)?;
    let umbrales = leer_umbrales(entrada)?;
    let es_ecoli = referencia.as_slice() == ecoli::genoma();
    let genes: &[Gen] = if es_ecoli { ecoli::genes() } else { &[] };
    let verdad = if entrada.tiene("verdad") {
        let ruta = entrada.ruta("verdad")?;
        Some(variantes::leer_cambios_vcf(formatos::abrir(&ruta)?)?)
    } else {
        None
    };
    let (pasos, _) = narrar(&ruta_sam, &referencia, genes, &umbrales, verdad.as_deref())?;
    Ok(pasos)
}

fn narrar(
    ruta_sam: &Path,
    referencia: &[u8],
    genes: &[Gen],
    umbrales: &Umbrales,
    verdad: Option<&[Cambio]>,
) -> Result<(Vec<Paso>, Vec<Llamada>), String> {
    if referencia.is_empty() {
        return Err("La referencia está vacía.".to_string());
    }
    let mayusculas = referencia.to_ascii_uppercase();
    let referencia: &[u8] = &mayusculas;
    let apilamiento = variantes::apilar(formatos::abrir(ruta_sam)?, referencia, umbrales)?;
    let ejemplos = buscar_ejemplos(
        ruta_sam,
        &apilamiento.nombre_referencia,
        referencia.len(),
        umbrales,
    )?;
    let escaneo = escanear(&apilamiento, referencia, umbrales);
    let snps = variantes::llamar_snps(&apilamiento, referencia, umbrales);
    let candidatos = candidatos_de_indel(&apilamiento, umbrales);
    let mut llamadas = snps.clone();
    llamadas.extend(variantes::llamar_indels(&apilamiento, umbrales));
    let total_antes_de_ordenar = llamadas.len();
    llamadas.sort_by(|una, otra| una.cambio.cmp(&otra.cambio));
    let indice = IndiceDeGenes::nuevo(genes, referencia.len());
    for llamada in &mut llamadas {
        llamada.anotaciones = variantes::anotar(&llamada.cambio, referencia, &indice);
    }
    let conteo = Conteo::de(&llamadas);
    let largo = format!("{} pb", miles(referencia.len()));
    let profundidad_minima = umbrales.profundidad_minima;
    let mut narracion = Narracion::nueva(400);

    narracion.agregar(
        Pizarra::nueva()
            .variable("referencia", &largo)
            .variable("umbrales.profundidad_minima", profundidad_minima)
            .variable("umbrales.fraccion_minima", decimal(umbrales.fraccion_minima, 2))
            .variable(
                "umbrales.fraccion_minima_indel",
                decimal(umbrales.fraccion_minima_indel, 2),
            )
            .variable("umbrales.calidad_mapeo_minima", umbrales.calidad_mapeo_minima)
            .variable("umbrales.calidad_base_minima", umbrales.calidad_base_minima)
            .secuencia(referencia, &[(0, referencia.len().min(30), Marca::Actual)])
            .paso(
                "revisar_umbrales(referencia, umbrales)?;",
                format!(
                    "Tienes una referencia de {largo} y un archivo de lecturas ya alineadas. El objetivo es encontrar los lugares donde las lecturas dicen otra cosa que la referencia. Lo primero que hace el código es revisar que los umbrales (cuántas lecturas, qué fracción, qué calidad) tengan sentido."
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("referencia", &largo)
            .secuencia(referencia, &[(0, referencia.len().min(30), Marca::Actual)])
            .paso(
                "let referencia = &referencia.to_ascii_uppercase();",
                "Pasa toda la referencia a mayúsculas para que «a» y «A» cuenten como la misma base al comparar.",
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("referencia", &largo)
            .variable("umbrales.calidad_mapeo_minima", umbrales.calidad_mapeo_minima)
            .variable("umbrales.calidad_base_minima", umbrales.calidad_base_minima)
            .paso(
                "let apilamiento = apilar(formatos::abrir(ruta_sam)?, referencia, umbrales)?;",
                format!(
                    "Primera etapa: apilar. Por cada posición del genoma el código lleva cuatro contadores (A, C, G, T) y va leyendo las lecturas una por una, sumando 1 en el contador de lo que cada lectura vio en cada posición. Al terminar tendrás una pila de {} lecturas utilizables y {} descartadas.",
                    miles(apilamiento.lecturas_usadas),
                    miles(apilamiento.lecturas_descartadas)
                ),
            ),
    );
    let mut calidades: Vec<usize> = vec![10, 20, 30, 40, umbrales.calidad_base_minima as usize];
    calidades.sort_unstable();
    calidades.dedup();
    let resaltada = calidades
        .iter()
        .position(|&calidad| calidad == umbrales.calidad_base_minima as usize);
    narracion.agregar(
        Pizarra::nueva()
            .variable("umbrales.calidad_base_minima", umbrales.calidad_base_minima)
            .lista(
                "Q (calidad Phred) → probabilidad de que la base esté mal",
                calidades
                    .iter()
                    .map(|&calidad| {
                        format!(
                            "Q{calidad} → {}",
                            porcentaje(variantes::probabilidad_de_error(calidad as u8), 2)
                        )
                    })
                    .collect(),
                resaltada.as_slice(),
            )
            .paso(
                "let errores = tabla_de_errores();",
                format!(
                    "Antes de leer, el código calcula una vez una tablita: cada calidad Q de una base se convierte en una probabilidad de error (Q20 es 1 error cada 100, Q30 uno cada 1 000). Tu mínimo es Q{}: las bases con menos que eso no se cuentan.",
                    umbrales.calidad_base_minima
                ),
            ),
    );

    if let Some(lectura) = &ejemplos.usada {
        let inicio = lectura.posicion - 1;
        let marca = [(
            inicio,
            (inicio + tramo_en_referencia(lectura)).min(referencia.len()),
            Marca::Actual,
        )];
        narracion.agregar(
            Pizarra::nueva()
                .variable("lectura.posicion", lectura.posicion)
                .variable("lectura.bandera", lectura.bandera)
                .variable("lectura.calidad_mapeo", lectura.calidad_mapeo)
                .variable("umbrales.calidad_mapeo_minima", umbrales.calidad_mapeo_minima)
                .secuencia(referencia, &marca)
                .paso(
                    "if !sirve(",
                    format!(
                        "Veamos la primera lectura útil ({}). Está mapeada en la {}, ocupa {} bases del genoma y su calidad de mapeo es {}: pasa todos los filtros, así que se usa.",
                        corto(&String::from_utf8_lossy(&lectura.nombre), 24),
                        lugar(inicio),
                        tramo_en_referencia(lectura),
                        lectura.calidad_mapeo
                    ),
                ),
        );
        narracion.agregar(
            Pizarra::nueva()
                .variable("lectura.posicion", lectura.posicion)
                .lista(
                    "CIGAR de la lectura (M = base alineada, I = insertada, D = borrada, S = recortada)",
                    cigar_en_lista(lectura),
                    &[],
                )
                .secuencia(referencia, &marca)
                .paso(
                    "for &(largo, operacion) in &lectura.cigar {",
                    "El CIGAR es la receta del alineamiento: dice, tramo por tramo, cómo se corresponde la lectura con la referencia. El código lo recorre en orden y hace algo distinto según la letra.",
                ),
        );
        if let Some((en_referencia, en_lectura, largo_tramo)) = primer_tramo_alineado(lectura) {
            let filas: Vec<Vec<String>> = (0..largo_tramo.min(BASES_EN_LA_TABLA))
                .map(|i| {
                    let base_lectura = lectura
                        .secuencia
                        .get(en_lectura + i)
                        .copied()
                        .unwrap_or(b'N');
                    let base_referencia =
                        referencia.get(en_referencia + i).copied().unwrap_or(b'N');
                    let calidad = calidad_de(lectura, en_lectura + i);
                    let cuenta = calidad >= umbrales.calidad_base_minima
                        && indice_base(base_lectura).is_some();
                    vec![
                        miles(en_referencia + i + 1),
                        (base_referencia as char).to_string(),
                        (base_lectura as char).to_string(),
                        format!("Q{calidad}"),
                        if cuenta { "sí" } else { "no" }.to_string(),
                    ]
                })
                .collect();
            let diferentes: Vec<(usize, usize)> = filas
                .iter()
                .enumerate()
                .filter(|(_, fila)| fila[1] != fila[2])
                .map(|(i, _)| (i, 2))
                .collect();
            narracion.agregar(
                Pizarra::nueva()
                    .variable("lectura.posicion", lectura.posicion)
                    .variable("umbrales.calidad_base_minima", umbrales.calidad_base_minima)
                    .tabla(
                        "Primeras bases de la lectura contra la referencia",
                        &["posición", "referencia", "lectura", "calidad", "¿cuenta?"],
                        filas,
                        &diferentes,
                    )
                    .secuencia(referencia, &marca)
                    .paso(
                        "contar_alineadas(",
                        format!(
                            "En el tramo alineado (M) el código mira base por base. Si la calidad llega al mínimo, suma 1 en el contador de esa letra para esa posición del genoma. De las {} bases de esta lectura cuentan {}.",
                            lectura.secuencia.len(),
                            bases_que_cuentan(lectura, umbrales)
                        ),
                    ),
            );
        }
    }
    if let Some((lectura, motivo)) = &ejemplos.descartada {
        narracion.agregar(
            Pizarra::nueva()
                .variable("lectura.posicion", lectura.posicion)
                .variable("lectura.bandera", lectura.bandera)
                .variable("lectura.calidad_mapeo", lectura.calidad_mapeo)
                .variable("umbrales.calidad_mapeo_minima", umbrales.calidad_mapeo_minima)
                .paso(
                    "if !sirve(",
                    format!(
                        "Esta lectura ({}) no pasa el filtro: {motivo}. Se cuenta como descartada y no aporta nada a la pila.",
                        corto(&String::from_utf8_lossy(&lectura.nombre), 24)
                    ),
                ),
        );
    }
    if let Some(lectura) = &ejemplos.con_indel
        && let Some(indel) = primer_indel(lectura)
        && let Some(cambio) = cambio_del_indel(lectura, &indel, referencia)
    {
        let es_insercion = indel.operacion == b'I';
        let posicion_en_cigar = lectura
            .cigar
            .iter()
            .position(|&(_, operacion)| operacion == indel.operacion)
            .map(|i| [i])
            .unwrap_or([0]);
        let marca = [(
            cambio.posicion,
            (cambio.posicion + 1).min(referencia.len()),
            Marca::Actual,
        )];
        narracion.agregar(
            Pizarra::nueva()
                .variable("lectura.posicion", lectura.posicion)
                .lista("CIGAR de esta lectura", cigar_en_lista(lectura), &posicion_en_cigar)
                .secuencia(referencia, &marca)
                .paso(
                    if es_insercion { "b'I' => {" } else { "b'D' => {" },
                    if es_insercion {
                        format!(
                            "Otra lectura trae una I en su CIGAR: {} bases que están en la lectura pero no en la referencia (una inserción). No suman a los contadores de A, C, G, T; el código las anota aparte, como un posible indel anclado en la {}.",
                            indel.largo,
                            lugar(cambio.posicion)
                        )
                    } else {
                        format!(
                            "Otra lectura trae una D en su CIGAR: {} bases de la referencia que la lectura se salta (una deleción). El código las anota como un posible indel anclado en la {} y además suma 1 a los «borrados» de cada posición saltada.",
                            indel.largo,
                            lugar(cambio.posicion)
                        )
                    },
                ),
        );
        let normalizado = variantes::normalizar(&cambio, referencia);
        let apoyo = apilamiento.indels.get(&normalizado).copied().unwrap_or(0);
        narracion.agregar(
            Pizarra::nueva()
                .variable("ancla", cambio.posicion)
                .variable("cambio", cambio_en_texto(&cambio))
                .secuencia(referencia, &marca)
                .paso(
                    if es_insercion {
                        "registrar_insercion(apilamiento, referencia, en_referencia - 1, insertadas);"
                    } else {
                        "registrar_delecion(apilamiento, referencia, en_referencia - 1, largo);"
                    },
                    format!(
                        "El indel se escribe como lo hace un VCF: la base de la referencia justo antes (el ancla, la {}) y lo que cambia: {}.",
                        lugar(cambio.posicion),
                        cambio_en_texto(&cambio)
                    ),
                ),
        );
        narracion.agregar(
            Pizarra::nueva()
                .variable("cambio", cambio_en_texto(&cambio))
                .variable("normalizado", format!("{} en la {}", cambio_en_texto(&normalizado), lugar(normalizado.posicion)))
                .variable("apoyo", apoyo)
                .secuencia(
                    referencia,
                    &[(
                        normalizado.posicion,
                        (normalizado.posicion + 1).min(referencia.len()),
                        Marca::Hallado,
                    )],
                )
                .paso(
                    ".entry(normalizar(cambio, referencia))",
                    format!(
                        "Un mismo indel en una zona repetida se puede escribir en varios lugares. Normalizar lo mueve siempre al mismo sitio para que todas las lecturas que lo ven sumen en la misma cuenta. Al final, {} lecturas apoyan este indel.",
                        miles(apoyo as usize)
                    ),
                ),
        );
    }
    if let Some(posicion) = escaneo.ejemplo_coinciden {
        let profundidad = apilamiento.profundidad(posicion);
        let indice_referencia = indice_base(referencia[posicion]).unwrap_or(0);
        narracion.agregar(
            tabla_de_pila(
                Pizarra::nueva()
                    .variable("posicion", posicion)
                    .variable("profundidad", profundidad),
                &apilamiento,
                posicion,
                &[(0, indice_referencia)],
            )
            .secuencia(referencia, &[(posicion, posicion + 1, Marca::Actual)])
            .paso(
                "columna[indice as usize] += 1;",
                format!(
                    "Así queda la pila en la {} cuando ya pasaron todas las lecturas: {} lecturas cubren este lugar y todas dicen la misma base que la referencia ({}). Esa tablita de cuatro contadores, repetida para cada posición del genoma, es todo lo que el resto del algoritmo necesita.",
                    lugar(posicion),
                    miles(profundidad as usize),
                    referencia[posicion] as char
                ),
            ),
        );
    }
    narracion.agregar(
        Pizarra::nueva()
            .variable("apilamiento.lecturas_usadas", miles(apilamiento.lecturas_usadas))
            .variable(
                "apilamiento.lecturas_descartadas",
                miles(apilamiento.lecturas_descartadas),
            )
            .variable("apilamiento.nombre_referencia", &apilamiento.nombre_referencia)
            .paso(
                "Ok(apilamiento)",
                format!(
                    "La pila está lista: se usaron {} lecturas y se descartaron {}. {} de las {} posiciones tienen al menos una lectura.",
                    miles(apilamiento.lecturas_usadas),
                    miles(apilamiento.lecturas_descartadas),
                    miles(escaneo.con_lecturas),
                    miles(referencia.len())
                ),
            ),
    );

    narracion.agregar(
        Pizarra::nueva()
            .variable("umbrales.profundidad_minima", profundidad_minima)
            .variable("umbrales.fraccion_minima", decimal(umbrales.fraccion_minima, 2))
            .secuencia(referencia, &[(0, referencia.len().min(30), Marca::Actual)])
            .paso(
                "for (posicion, &base_referencia) in referencia.iter().enumerate().take(apilamiento.bases.len())",
                "Segunda etapa: buscar SNPs (cambios de una sola base). El código recorre las posiciones del genoma de la primera a la última y en cada una se hace tres preguntas, en orden: ¿hay suficientes lecturas? ¿alguna base distinta a la referencia? ¿esa base distinta está en la mayoría de las lecturas?",
            ),
    );
    if let Some(posicion) = escaneo.ejemplo_pocas.or(escaneo.ejemplo_sin_lecturas) {
        let profundidad = apilamiento.profundidad(posicion);
        narracion.agregar(
            Pizarra::nueva()
                .variable("posicion", posicion)
                .variable("profundidad", profundidad)
                .variable("umbrales.profundidad_minima", profundidad_minima)
                .secuencia(referencia, &[(posicion, posicion + 1, Marca::Actual)])
                .paso(
                    "if profundidad < umbrales.profundidad_minima {",
                    format!(
                        "Pregunta 1 en la {}: solo hay {} y pediste al menos {}. Con tan pocas no se puede distinguir un cambio real de un error, así que se salta. Esto pasa en {} posiciones.",
                        lugar(posicion),
                        plural(profundidad as usize, "lectura", "lecturas"),
                        profundidad_minima,
                        miles(escaneo.pocas)
                    ),
                ),
        );
    }
    if let Some(posicion) = escaneo.ejemplo_coinciden {
        let profundidad = apilamiento.profundidad(posicion);
        let indice_referencia = indice_base(referencia[posicion]).unwrap_or(0);
        narracion.agregar(
            tabla_de_pila(
                Pizarra::nueva()
                    .variable("posicion", posicion)
                    .variable("profundidad", profundidad)
                    .variable("a_favor", 0),
                &apilamiento,
                posicion,
                &[(0, indice_referencia)],
            )
            .secuencia(referencia, &[(posicion, posicion + 1, Marca::Actual)])
            .paso(
                "if a_favor == 0 {",
                format!(
                    "Pregunta 2 en la {}: hay lecturas de sobra, pero ninguna base distinta a la referencia ({}). No hay nada que llamar. Esto es lo normal: {} posiciones son así.",
                    lugar(posicion),
                    referencia[posicion] as char,
                    miles(escaneo.coinciden)
                ),
            ),
        );
    }
    if let Some(posicion) = escaneo.ejemplo_ruido {
        let profundidad = apilamiento.profundidad(posicion);
        let indice_referencia = indice_base(referencia[posicion]).unwrap_or(0);
        let (indice_alternativo, a_favor) =
            variantes::alternativa_mas_frecuente(&apilamiento.bases[posicion], indice_referencia);
        let fraccion = a_favor as f64 / profundidad as f64;
        narracion.agregar(
            tabla_de_pila(
                Pizarra::nueva()
                    .variable("posicion", posicion)
                    .variable("profundidad", profundidad)
                    .variable("a_favor", a_favor)
                    .variable("fraccion", decimal(fraccion, 2)),
                &apilamiento,
                posicion,
                &[(0, indice_referencia), (0, indice_alternativo)],
            )
            .secuencia(referencia, &[(posicion, posicion + 1, Marca::Actual)])
            .paso(
                "if fraccion < umbrales.fraccion_minima {",
                format!(
                    "Pregunta 3 en la {}: {} de {} lecturas dicen otra base ({}), pero tu mínimo es {}. Una minoría así suele ser un error de secuenciación, no una mutación, así que se descarta. Esto pasa en {} posiciones.",
                    lugar(posicion),
                    a_favor,
                    profundidad,
                    porcentaje(fraccion, 0),
                    porcentaje(umbrales.fraccion_minima, 0),
                    miles(escaneo.ruido)
                ),
            ),
        );
    }
    for (orden, llamada) in snps.iter().take(3).enumerate() {
        let posicion = llamada.cambio.posicion;
        let indice_referencia = indice_base(llamada.cambio.referencia[0]).unwrap_or(0);
        let indice_alternativo = indice_base(llamada.cambio.alternativa[0]).unwrap_or(0);
        narracion.agregar(
            tabla_de_pila(
                Pizarra::nueva()
                    .variable("posicion", posicion)
                    .variable("base_referencia", llamada.cambio.referencia[0] as char)
                    .variable("profundidad", llamada.profundidad)
                    .variable("a_favor", llamada.a_favor)
                    .variable("fraccion", decimal(llamada.fraccion, 2)),
                &apilamiento,
                posicion,
                &[(0, indice_referencia), (0, indice_alternativo)],
            )
            .secuencia(referencia, &[(posicion, posicion + 1, Marca::Actual)])
            .paso(
                "llamadas.push(Llamada {",
                format!(
                    "{} En la {}, {} de {} lecturas ({}) dicen {} donde la referencia tiene {}. Pasa las tres preguntas: es un SNP y se guarda.",
                    if orden == 0 { "¡Encontramos uno!" } else { "Otro SNP." },
                    lugar(posicion),
                    llamada.a_favor,
                    llamada.profundidad,
                    porcentaje(llamada.fraccion, 0),
                    llamada.cambio.alternativa[0] as char,
                    llamada.cambio.referencia[0] as char
                ),
            ),
        );
        if orden == 0 {
            let error_hacia_esa_base =
                variantes::error_medio(&apilamiento, posicion, umbrales) / 3.0;
            narracion.agregar(
                Pizarra::nueva()
                    .variable("posicion", posicion)
                    .variable("profundidad", llamada.profundidad)
                    .variable("a_favor", llamada.a_favor)
                    .variable("error_hacia_esa_base", cientifico(error_hacia_esa_base))
                    .variable("calidad", decimal(llamada.calidad, 1))
                    .secuencia(referencia, &[(posicion, posicion + 1, Marca::Actual)])
                    .paso(
                        "calidad: calidad_phred(a_favor, profundidad, error_hacia_esa_base),",
                        format!(
                            "QUAL mide qué tan seguro es el SNP: es la probabilidad de que {} de {} lecturas se equivoquen todas hacia la misma base por puro error de secuenciación, escrita en escala Phred. Aquí sale {}: cuanto más alto, menos probable que sea un error.",
                            llamada.a_favor,
                            llamada.profundidad,
                            decimal(llamada.calidad, 1)
                        ),
                    ),
            );
        }
    }
    narracion.agregar(
        Pizarra::nueva()
            .variable("llamadas", format!("{} SNPs", miles(snps.len())))
            .paso(
                "let mut llamadas = llamar_snps(&apilamiento, referencia, umbrales);",
                format!(
                    "Terminó el recorrido. De {} posiciones: {} con pocas lecturas, {} donde todos coinciden con la referencia, {} con una minoría sospechosa y {} SNPs confirmados.",
                    miles(referencia.len()),
                    miles(escaneo.pocas + escaneo.no_acgt),
                    miles(escaneo.coinciden),
                    miles(escaneo.ruido),
                    miles(snps.len())
                ),
            ),
    );

    narracion.agregar(
        Pizarra::nueva()
            .variable("llamadas", format!("{} SNPs", miles(snps.len())))
            .variable("apilamiento.indels", miles(apilamiento.indels.len()))
            .paso(
                "for (cambio, &a_favor) in &apilamiento.indels {",
                format!(
                    "Tercera etapa: indels (inserciones y deleciones). Durante el apilado el código fue anotando cada indel que vio en el CIGAR de las lecturas, y cuántas lecturas lo apoyan. Hay {} candidatos distintos y se revisa cada uno.",
                    miles(apilamiento.indels.len())
                ),
            ),
    );
    if let Some(candidato) = candidatos.iter().find(|candidato| candidato.pasa) {
        narracion.agregar(
            Pizarra::nueva()
                .variable("cambio", cambio_en_texto(&candidato.cambio))
                .variable("a_favor", candidato.a_favor)
                .variable("cruzan", candidato.cruzan)
                .variable("profundidad", candidato.profundidad)
                .variable("fraccion", decimal(candidato.fraccion, 2))
                .secuencia(
                    referencia,
                    &[(candidato.cambio.posicion, candidato.cambio.posicion + 1, Marca::Actual)],
                )
                .paso(
                    "let fraccion = a_favor as f64 / profundidad as f64;",
                    format!(
                        "Candidato en la {}: {}. {} lecturas lo traen y {} lecturas cruzan ese punto, así que la fracción es {}. A un indel le basta {} porque las lecturas con indels se alinean peor. Pasa.",
                        lugar(candidato.cambio.posicion),
                        cambio_en_texto(&candidato.cambio),
                        candidato.a_favor,
                        candidato.cruzan,
                        porcentaje(candidato.fraccion, 0),
                        porcentaje(umbrales.fraccion_minima_indel, 0)
                    ),
                ),
        );
    }
    if let Some(candidato) = candidatos.iter().find(|candidato| !candidato.pasa) {
        narracion.agregar(
            Pizarra::nueva()
                .variable("cambio", cambio_en_texto(&candidato.cambio))
                .variable("a_favor", candidato.a_favor)
                .variable("cruzan", candidato.cruzan)
                .variable("profundidad", candidato.profundidad)
                .variable("fraccion", decimal(candidato.fraccion, 2))
                .secuencia(
                    referencia,
                    &[(candidato.cambio.posicion, candidato.cambio.posicion + 1, Marca::Actual)],
                )
                .paso(
                    "if profundidad < umbrales.profundidad_minima || fraccion < umbrales.fraccion_minima_indel {",
                    format!(
                        "Otro candidato, en la {}: {}. Solo {} de {} lecturas lo traen ({}), no alcanza los umbrales. Se descarta.",
                        lugar(candidato.cambio.posicion),
                        cambio_en_texto(&candidato.cambio),
                        candidato.a_favor,
                        candidato.profundidad,
                        porcentaje(candidato.fraccion, 0)
                    ),
                ),
        );
    }
    narracion.agregar(
        Pizarra::nueva()
            .variable("llamadas", format!("{} variantes", miles(total_antes_de_ordenar)))
            .paso(
                "llamadas.dedup_by_key(|llamada| llamada.cambio.posicion);",
                format!(
                    "Si dos indels distintos caen en la misma posición, solo se queda el que tiene más lecturas. De {} candidatos que pasaron los umbrales quedan {} indels.",
                    miles(candidatos.iter().filter(|candidato| candidato.pasa).count()),
                    miles(conteo.indels())
                ),
            ),
    );
    narracion.agregar(
        Pizarra::nueva()
            .variable("llamadas", format!("{} variantes", miles(llamadas.len())))
            .paso(
                "llamadas.sort_by(|una, otra| una.cambio.cmp(&otra.cambio));",
                format!(
                    "Se juntan los {} SNPs con los {} indels y se ordenan por posición en el genoma.",
                    miles(conteo.snps()),
                    miles(conteo.indels())
                ),
            ),
    );

    narracion.agregar(
        Pizarra::nueva()
            .variable("genes", miles(genes.len()))
            .variable("llamadas", format!("{} variantes", miles(llamadas.len())))
            .paso(
                "let indice = IndiceDeGenes::nuevo(genes, referencia.len());",
                if genes.is_empty() {
                    "Cuarta etapa: anotar, decir qué le hace cada variante a los genes. Como la referencia no es el genoma de E. coli incluido, no hay genes que consultar: las variantes quedan sin anotar.".to_string()
                } else {
                    format!(
                        "Cuarta etapa: anotar, decir qué le hace cada variante a los genes. Primero se ordenan los {} genes de E. coli por posición para encontrar rápido cuál toca cada variante.",
                        miles(genes.len())
                    )
                },
            ),
    );
    let mut vistos: Vec<&'static str> = Vec::new();
    let mut elegidas: Vec<&Llamada> = Vec::new();
    for llamada in &llamadas {
        let Some(anotacion) = llamada.anotaciones.first() else {
            continue;
        };
        let nombre = nombre_del_efecto(&anotacion.efecto);
        if vistos.contains(&nombre) || elegidas.len() >= EJEMPLOS_DE_ANOTACION {
            continue;
        }
        vistos.push(nombre);
        elegidas.push(llamada);
    }
    if let Some(resistente) = llamadas.iter().find(|llamada| {
        llamada
            .anotaciones
            .iter()
            .any(|anotacion| anotacion.resistencia.is_some())
    }) && !elegidas.contains(&resistente)
    {
        elegidas.push(resistente);
    }
    let mut codones_mostrados = false;
    for llamada in elegidas {
        let anotacion = &llamada.anotaciones[0];
        let posicion = llamada.cambio.posicion;
        if !codones_mostrados
            && llamada.cambio.es_snp()
            && matches!(
                anotacion.efecto,
                Efecto::CambioDeAminoacido | Efecto::Sinonima | Efecto::ParadaPrematura
            )
            && let Some(elegido) = gen_que_contiene(genes, &anotacion.nombre_gen, posicion)
        {
            codones_mostrados = true;
            let dentro = if elegido.hebra == b'-' {
                elegido.fin - 1 - posicion
            } else {
                posicion - (elegido.inicio - 1)
            };
            let numero = dentro / 3 + 1;
            let (desde, hasta) = posiciones_del_codon(elegido, numero);
            let antes = variantes::codon(referencia, elegido, numero, None);
            let despues = variantes::codon(referencia, elegido, numero, Some(&llamada.cambio));
            if let (Some(antes), Some(despues)) = (antes, despues) {
                let hebra = if elegido.hebra == b'-' {
                    " El gen está en la hebra reversa, por eso se lee el complemento."
                } else {
                    ""
                };
                narracion.agregar(
                    Pizarra::nueva()
                        .variable("elegido.nombre", &elegido.nombre)
                        .variable("numero", numero)
                        .variable("antes", String::from_utf8_lossy(&antes))
                        .secuencia(
                            referencia,
                            &[(desde, hasta, Marca::Comparado), (posicion, posicion + 1, Marca::Actual)],
                        )
                        .paso(
                            "codon(referencia, elegido, numero, None),",
                            format!(
                                "La {} cae dentro del gen {}, en el codón {numero}. Con la referencia, ese codón es {} ({}).{hebra}",
                                lugar(posicion),
                                elegido.nombre,
                                String::from_utf8_lossy(&antes),
                                aminoacido(&antes) as char
                            ),
                        ),
                );
                narracion.agregar(
                    Pizarra::nueva()
                        .variable("elegido.nombre", &elegido.nombre)
                        .variable("numero", numero)
                        .variable("antes", String::from_utf8_lossy(&antes))
                        .variable("despues", String::from_utf8_lossy(&despues))
                        .secuencia(
                            referencia,
                            &[(desde, hasta, Marca::Comparado), (posicion, posicion + 1, Marca::Actual)],
                        )
                        .paso(
                            "codon(referencia, elegido, numero, Some(cambio)),",
                            format!(
                                "Se arma el mismo codón pero con la base cambiada: ahora es {} ({}). Comparar los dos aminoácidos dice si la proteína cambia.",
                                String::from_utf8_lossy(&despues),
                                aminoacido(&despues) as char
                            ),
                        ),
                );
            }
        }
        narracion.agregar(
            Pizarra::nueva()
                .variable("posicion", posicion)
                .variable("cambio", cambio_en_texto(&llamada.cambio))
                .variable("anotacion.nombre_gen", nombre_del_lugar(anotacion))
                .variable("anotacion.efecto", nombre_del_efecto(&anotacion.efecto))
                .variable(
                    "anotacion.proteina",
                    if anotacion.proteina.is_empty() {
                        "—"
                    } else {
                        &anotacion.proteina
                    },
                )
                .variable("anotacion.resistencia", anotacion.resistencia.unwrap_or("—"))
                .secuencia(referencia, &[(posicion, posicion + 1, Marca::Actual)])
                .paso(
                    "llamada.anotaciones = anotar(&llamada.cambio, referencia, &indice);",
                    format!(
                        "Variante en la {} ({}). {}{}{}",
                        lugar(posicion),
                        cambio_en_texto(&llamada.cambio),
                        match &anotacion.efecto {
                            Efecto::Intergenica { .. } => format!(
                                "Cae {}: no toca ninguna proteína, pero podría estar en un promotor.",
                                nombre_del_lugar(anotacion)
                            ),
                            efecto => format!(
                                "Cae en el gen {} con efecto «{}»{}.",
                                anotacion.nombre_gen,
                                nombre_del_efecto(efecto),
                                if anotacion.proteina.is_empty() {
                                    String::new()
                                } else {
                                    format!(" (cambio {})", anotacion.proteina)
                                }
                            ),
                        },
                        if anotacion.resistencia.is_some() {
                            " Está en un sitio conocido de resistencia a antibióticos."
                        } else {
                            ""
                        },
                        if corta_la_proteina(&anotacion.efecto) {
                            " Corta la proteína, así que probablemente inactiva el gen."
                        } else {
                            ""
                        }
                    ),
                ),
        );
    }

    let mut ultimo = Pizarra::nueva()
        .variable("llamadas", format!("{} variantes", miles(llamadas.len())))
        .variable("lecturas_usadas", miles(apilamiento.lecturas_usadas))
        .variable(
            "lecturas_descartadas",
            miles(apilamiento.lecturas_descartadas),
        );
    let mut nota = format!(
        "Resultado: {} SNPs y {} indels ({} y {}), a partir de {} lecturas. Mira «Qué significa» para saber cómo leer cada número y la tabla de variantes.",
        miles(conteo.snps()),
        miles(conteo.indels()),
        plural(conteo.inserciones, "inserción", "inserciones"),
        plural(conteo.deleciones, "deleción", "deleciones"),
        miles(apilamiento.lecturas_usadas)
    );
    if let Some(verdad) = verdad {
        let cambios: Vec<Cambio> = llamadas
            .iter()
            .map(|llamada| llamada.cambio.clone())
            .collect();
        let comparacion = variantes::comparar_con_verdad(&cambios, verdad, referencia);
        if let (Some(sensibilidad), Some(precision)) =
            (comparacion.sensibilidad(), comparacion.precision())
        {
            ultimo = ultimo
                .variable("sensibilidad", porcentaje(sensibilidad, 1))
                .variable("precision", porcentaje(precision, 1));
            nota.push_str(&format!(
                " Contra las mutaciones verdaderas: sensibilidad {} y precisión {}.",
                porcentaje(sensibilidad, 1),
                porcentaje(precision, 1)
            ));
        }
    }
    let informe_final = ultimo.paso("Ok(Informe {", nota);
    Ok((narracion.terminar(informe_final), llamadas))
}

#[cfg(test)]
mod pruebas {
    use super::*;

    #[test]
    fn la_narracion_llega_al_mismo_resultado_que_el_nucleo() {
        let entrada =
            Entrada::desde_json(crate::ejemplo("puerto.llamar_variantes").unwrap()).unwrap();
        let ruta_sam = entrada.ruta("alineamientos").unwrap();
        let referencia = leer_referencia(&entrada).unwrap();
        let umbrales = leer_umbrales(&entrada).unwrap();
        let genes = ecoli::genes();
        let verdad =
            variantes::leer_cambios_vcf(formatos::abrir(&entrada.ruta("verdad").unwrap()).unwrap())
                .unwrap();
        let (pasos, llamadas) =
            narrar(&ruta_sam, &referencia, genes, &umbrales, Some(&verdad)).unwrap();
        let informe =
            variantes::llamar_variantes(&ruta_sam, &referencia, genes, &umbrales).unwrap();
        assert_eq!(llamadas.len(), informe.llamadas.len());
        for (propia, del_nucleo) in llamadas.iter().zip(&informe.llamadas) {
            assert_eq!(propia, del_nucleo);
        }
        assert!(!pasos.is_empty() && pasos.len() <= 400);
    }
}
