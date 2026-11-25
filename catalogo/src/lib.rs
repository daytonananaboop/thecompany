mod ejemplos;
mod entrada;
pub mod formato;
mod fuente;
pub mod pizarra;
mod salida;

mod agrupamiento;
mod alineamiento;
mod antibioticos;
mod busqueda;
mod calidad;
mod cebadores;
mod ensamblaje;
mod filogenia;
mod genes;
mod mapeo;
mod motivos;
mod reordenamientos;
mod replicacion;
mod simulacion;
mod variantes;

mod comunidades;
mod dcj;
mod ensamblador;
mod frente_de_onda;
mod indice_fm;
mod minhash;
mod motivos_em;
mod proteomica;
mod verosimilitud;

pub use entrada::Entrada;
pub use fuente::{Fragmento, fragmentos};
pub use salida::{Arista, Bloque, FilaAlineada, Marca, Paso, Punto, Salida};

use serde::Serialize;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

#[derive(Serialize, Debug, Clone, Copy)]
pub struct Capitulo {
    pub id: &'static str,
    pub numero: u8,
    pub titulo: &'static str,
    pub resumen: &'static str,
}

#[derive(Serialize, Debug, Clone, Copy, PartialEq)]
#[serde(tag = "clase", content = "opciones", rename_all = "snake_case")]
pub enum Tipo {
    Secuencia,
    Secuencias,
    Proteina,
    Texto,
    Entero,
    Decimal,
    Numeros,
    Matriz,
    Archivo,
    Base,
    Opcion(&'static [&'static str]),
}

#[derive(Serialize, Debug, Clone, Copy)]
pub struct Campo {
    pub id: &'static str,
    pub etiqueta: &'static str,
    pub tipo: Tipo,
    pub ejemplo: &'static str,
}

pub type Ejecutor = fn(&Entrada) -> Result<Salida, String>;
pub type Trazador = fn(&Entrada) -> Result<Vec<Paso>, String>;

#[derive(Clone, Copy)]
pub struct Algoritmo {
    pub id: &'static str,
    pub capitulo: &'static str,
    pub nombre: &'static str,
    pub resumen: &'static str,
    pub campos: &'static [Campo],
    pub funcion: &'static str,
    pub ejecutar: Ejecutor,
    pub pasos: Option<Trazador>,
}

#[derive(Serialize, Debug, Clone)]
pub struct Ficha {
    pub id: &'static str,
    pub capitulo: &'static str,
    pub nombre: &'static str,
    pub resumen: &'static str,
    pub campos: &'static [Campo],
    pub funcion: &'static str,
    pub paso_a_paso: bool,
    pub vanguardia: Option<&'static str>,
}

const VANGUARDIA: &[(&str, &str)] = &[
    (
        "replicacion.indice_fm",
        "Ferragina y Manzini (2000), FOCS; arreglo de sufijos SA-IS: Nong, Zhang y Chan (2009), DCC.",
    ),
    (
        "replicacion.bwt",
        "Burrows y Wheeler (1994), SRC Research Report 124, Digital Equipment Corporation.",
    ),
    (
        "motivos.meme",
        "Bailey y Elkan (1994), ISMB 2:28–36; modelos OOPS y ZOOPS: Bailey y Elkan (1995), Machine Learning 21:51–80.",
    ),
    (
        "ensamblaje.ensamblador",
        "Zerbino y Birney (2008), Genome Research 18:821–829 (Velvet); Bankevich et al. (2012), Journal of Computational Biology 19:455–477 (SPAdes).",
    ),
    (
        "antibioticos.digestion",
        "Olsen, Ong y Mann (2004), Molecular & Cellular Proteomics 3:608–614; masas: Eng, McCormack y Yates (1994), J. Am. Soc. Mass Spectrom. 5:976–989.",
    ),
    (
        "antibioticos.busqueda_de_espectros",
        "Eng, McCormack y Yates (1994), J. Am. Soc. Mass Spectrom. 5:976–989; Craig y Beavis (2004), Bioinformatics 20:1466–1467; Elias y Gygi (2007), Nature Methods 4:207–214.",
    ),
    (
        "alineamiento.wfa",
        "Marco-Sola, Moure, Moreto y Espinosa (2021), Bioinformatics 37:456–463.",
    ),
    (
        "alineamiento.myers",
        "Myers (1999), Journal of the ACM 46:395–415; banda que se duplica: Ukkonen (1985), Information and Control 64:100–118.",
    ),
    (
        "reordenamientos.dcj",
        "Yancopoulos, Attie y Friedberg (2005), Bioinformatics 21:3340–3346; Bergeron, Mixtacki y Stoye (2006), WABI, LNCS 4175:163–173.",
    ),
    (
        "filogenia.maxima_verosimilitud",
        "Felsenstein (1981), Journal of Molecular Evolution 17:368–376; Hasegawa, Kishino y Yano (1985); Yang (1994); Guindon y Gascuel (2003), Systematic Biology 52:696–704.",
    ),
    (
        "filogenia.bootstrap",
        "Felsenstein (1985), Evolution 39:783–791.",
    ),
    (
        "agrupamiento.louvain",
        "Blondel, Guillaume, Lambiotte y Lefebvre (2008), Journal of Statistical Mechanics P10008; grafo de vecinos compartidos: Levine et al. (2015), Cell 162:184–197; resolución: Reichardt y Bornholdt (2006), Physical Review E 74:016110.",
    ),
    (
        "agrupamiento.tsne",
        "van der Maaten y Hinton (2008), Journal of Machine Learning Research 9:2579–2605; tasa de aprendizaje: Belkina et al. (2019), Nature Communications 10:5415.",
    ),
    (
        "puerto.minhash",
        "Broder (1997), Compression and Complexity of Sequences; Ondov et al. (2016), Genome Biology 17:132 (Mash).",
    ),
];

pub fn referencia_de_vanguardia(id: &str) -> Option<&'static str> {
    VANGUARDIA
        .iter()
        .find(|(clave, _)| *clave == id)
        .map(|(_, referencia)| *referencia)
}

pub fn nucleos() -> usize {
    rayon::current_num_threads()
}

const CAPITULOS: &[Capitulo] = &[
    Capitulo {
        id: "replicacion",
        numero: 1,
        titulo: "Origen de replicación",
        resumen: "Dónde empieza a copiarse un genoma bacteriano: patrones frecuentes, sesgo GC y cajas DnaA.",
    },
    Capitulo {
        id: "motivos",
        numero: 2,
        titulo: "Motivos regulatorios",
        resumen: "Señales cortas y parecidas escondidas en varias secuencias, buscadas con perfiles y azar.",
    },
    Capitulo {
        id: "ensamblaje",
        numero: 3,
        titulo: "Ensamblaje de genomas",
        resumen: "Reconstruir un genoma a partir de fragmentos con grafos de De Bruijn y caminos eulerianos.",
    },
    Capitulo {
        id: "antibioticos",
        numero: 4,
        titulo: "Secuenciación de péptidos",
        resumen: "Qué hebra de ADN codifica un péptido y cómo reconstruirlo desde su espectro de masas.",
    },
    Capitulo {
        id: "alineamiento",
        numero: 5,
        titulo: "Alineamiento de secuencias",
        resumen: "Comparar secuencias con programación dinámica: global, local, con huecos afines y en espacio lineal.",
    },
    Capitulo {
        id: "reordenamientos",
        numero: 6,
        titulo: "Reordenamientos genómicos",
        resumen: "Cuántas inversiones y roturas separan dos genomas, y dónde están los bloques compartidos.",
    },
    Capitulo {
        id: "filogenia",
        numero: 7,
        titulo: "Filogenia",
        resumen: "Árboles evolutivos construidos desde distancias o desde caracteres.",
    },
    Capitulo {
        id: "agrupamiento",
        numero: 8,
        titulo: "Expresión génica y agrupamiento",
        resumen: "Agrupar genes que cambian juntos su expresión.",
    },
    Capitulo {
        id: "puerto",
        numero: 9,
        titulo: "Herramientas de laboratorio",
        resumen: "Para datos reales: calidad de lecturas, búsqueda de parecidos, mapeo, variantes, genes y cebadores.",
    },
];

pub fn capitulos() -> &'static [Capitulo] {
    CAPITULOS
}

pub fn todos() -> &'static [Algoritmo] {
    static TODOS: OnceLock<Vec<Algoritmo>> = OnceLock::new();
    TODOS.get_or_init(|| {
        let mut lista = Vec::new();
        lista.extend(replicacion::algoritmos());
        lista.extend(indice_fm::algoritmos());
        lista.extend(motivos::algoritmos());
        lista.extend(motivos_em::algoritmos());
        lista.extend(ensamblaje::algoritmos());
        lista.extend(ensamblador::algoritmos());
        lista.extend(antibioticos::algoritmos());
        lista.extend(proteomica::algoritmos());
        lista.extend(alineamiento::algoritmos());
        lista.extend(frente_de_onda::algoritmos());
        lista.extend(reordenamientos::algoritmos());
        lista.extend(dcj::algoritmos());
        lista.extend(filogenia::algoritmos());
        lista.extend(verosimilitud::algoritmos());
        lista.extend(agrupamiento::algoritmos());
        lista.extend(comunidades::algoritmos());
        lista.extend(simulacion::algoritmos());
        lista.extend(calidad::algoritmos());
        lista.extend(busqueda::algoritmos());
        lista.extend(mapeo::algoritmos());
        lista.extend(variantes::algoritmos());
        lista.extend(genes::algoritmos());
        lista.extend(cebadores::algoritmos());
        lista.extend(minhash::algoritmos());
        lista
    })
}

pub fn preparar() {
    std::thread::scope(|hilos| {
        hilos.spawn(|| {
            nucleo::ecoli::genoma();
            nucleo::ecoli::genes();
        });
        hilos.spawn(busqueda::preparar);
        hilos.spawn(mapeo::preparar);
        hilos.spawn(genes::preparar);
        hilos.spawn(|| {
            indice_fm::preparar();
        });
    });
}

pub fn buscar(id: &str) -> Option<&'static Algoritmo> {
    todos().iter().find(|algoritmo| algoritmo.id == id)
}

pub fn fichas() -> Vec<Ficha> {
    todos()
        .iter()
        .map(|algoritmo| Ficha {
            id: algoritmo.id,
            capitulo: algoritmo.capitulo,
            nombre: algoritmo.nombre,
            resumen: algoritmo.resumen,
            campos: algoritmo.campos,
            funcion: algoritmo.funcion,
            paso_a_paso: algoritmo.pasos.is_some(),
            vanguardia: referencia_de_vanguardia(algoritmo.id),
        })
        .collect()
}

pub fn ejemplo(id: &str) -> Result<serde_json::Value, String> {
    let algoritmo = encontrar(id)?;
    let mut valores = serde_json::Map::new();
    for campo in algoritmo.campos {
        valores.insert(
            campo.id.to_string(),
            serde_json::Value::String(campo.ejemplo.to_string()),
        );
    }
    Ok(serde_json::Value::Object(valores))
}

fn encontrar(id: &str) -> Result<&'static Algoritmo, String> {
    buscar(id).ok_or_else(|| format!("No existe el algoritmo «{id}»."))
}

fn proteger<T>(trabajo: impl FnOnce() -> Result<T, String>) -> Result<T, String> {
    match catch_unwind(AssertUnwindSafe(trabajo)) {
        Ok(resultado) => resultado,
        Err(error) => {
            let mensaje = error
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| error.downcast_ref::<&str>().map(|texto| texto.to_string()))
                .unwrap_or_default();
            Err(format!("Error interno del algoritmo. {mensaje}"))
        }
    }
}

pub fn ejecutar(id: &str, entrada: serde_json::Value) -> Result<Salida, String> {
    let algoritmo = encontrar(id)?;
    let entrada = Entrada::desde_json(entrada)?;
    proteger(|| (algoritmo.ejecutar)(&entrada))
}

pub fn pasos(id: &str, entrada: serde_json::Value) -> Result<Vec<Paso>, String> {
    let algoritmo = encontrar(id)?;
    let trazador = algoritmo
        .pasos
        .ok_or_else(|| format!("«{}» no tiene paso a paso.", algoritmo.nombre))?;
    let entrada = Entrada::desde_json(entrada)?;
    proteger(|| trazador(&entrada))
}

pub fn codigo(id: &str) -> Result<Vec<Fragmento>, String> {
    let algoritmo = encontrar(id)?;
    fuente::fragmentos(algoritmo.funcion)
}

pub fn carpeta_de_ejemplos() -> Result<PathBuf, String> {
    let carpeta = std::env::temp_dir().join("la-compania");
    std::fs::create_dir_all(&carpeta)
        .map_err(|error| format!("No se pudo crear {}: {error}", carpeta.display()))?;
    Ok(carpeta)
}

pub fn ruta_de_salida(entrada: &Path, sufijo: &str, extension: &str) -> PathBuf {
    let carpeta = entrada.parent().map(Path::to_path_buf).unwrap_or_default();
    let mut base = entrada
        .file_name()
        .map(|nombre| nombre.to_string_lossy().to_string())
        .unwrap_or_else(|| "resultado".to_string());
    for terminacion in [
        ".gz", ".fastq", ".fq", ".fasta", ".fa", ".fna", ".sam", ".vcf",
    ] {
        if let Some(recortado) = base.strip_suffix(terminacion) {
            base = recortado.to_string();
        }
    }
    carpeta.join(format!("{base}.{sufijo}.{extension}"))
}
