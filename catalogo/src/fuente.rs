use serde::Serialize;
use std::collections::{HashMap, HashSet, VecDeque};

const FUENTES: &[(&str, &str)] = &[
    (
        "agrupamiento",
        include_str!("../../nucleo/src/agrupamiento.rs"),
    ),
    (
        "alineamiento",
        include_str!("../../nucleo/src/alineamiento.rs"),
    ),
    (
        "antibioticos",
        include_str!("../../nucleo/src/antibioticos.rs"),
    ),
    ("azar", include_str!("../../nucleo/src/azar.rs")),
    ("busqueda", include_str!("../../nucleo/src/busqueda.rs")),
    ("calidad", include_str!("../../nucleo/src/calidad.rs")),
    ("cebadores", include_str!("../../nucleo/src/cebadores.rs")),
    (
        "comunidades",
        include_str!("../../nucleo/src/comunidades.rs"),
    ),
    ("dcj", include_str!("../../nucleo/src/dcj.rs")),
    ("ecoli", include_str!("../../nucleo/src/ecoli.rs")),
    (
        "ensamblador",
        include_str!("../../nucleo/src/ensamblador.rs"),
    ),
    ("ensamblaje", include_str!("../../nucleo/src/ensamblaje.rs")),
    ("filogenia", include_str!("../../nucleo/src/filogenia.rs")),
    ("formatos", include_str!("../../nucleo/src/formatos.rs")),
    (
        "frente_de_onda",
        include_str!("../../nucleo/src/frente_de_onda.rs"),
    ),
    ("genes", include_str!("../../nucleo/src/genes.rs")),
    ("indice_fm", include_str!("../../nucleo/src/indice_fm.rs")),
    ("mapeo", include_str!("../../nucleo/src/mapeo.rs")),
    ("matrices", include_str!("../../nucleo/src/matrices.rs")),
    ("minhash", include_str!("../../nucleo/src/minhash.rs")),
    ("motivos", include_str!("../../nucleo/src/motivos.rs")),
    ("motivos_em", include_str!("../../nucleo/src/motivos_em.rs")),
    ("proteomica", include_str!("../../nucleo/src/proteomica.rs")),
    (
        "reordenamientos",
        include_str!("../../nucleo/src/reordenamientos.rs"),
    ),
    (
        "replicacion",
        include_str!("../../nucleo/src/replicacion.rs"),
    ),
    ("secuencia", include_str!("../../nucleo/src/secuencia.rs")),
    ("simulacion", include_str!("../../nucleo/src/simulacion.rs")),
    ("variantes", include_str!("../../nucleo/src/variantes.rs")),
    (
        "verosimilitud",
        include_str!("../../nucleo/src/verosimilitud.rs"),
    ),
];

#[derive(Serialize, Debug, Clone, PartialEq)]
pub struct Fragmento {
    pub nombre: String,
    pub fuente: String,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Clase {
    Funcion,
    Tipo,
    Implementacion,
    Constante,
    Importacion,
    Otra,
}

#[derive(Debug, Clone)]
struct Definicion {
    nombre: String,
    clase: Clase,
    texto: String,
    metodos: Vec<String>,
}

struct Modulo {
    definiciones: Vec<Definicion>,
    importados: HashMap<String, String>,
}

pub fn fragmentos(funcion: &str) -> Result<Vec<Fragmento>, String> {
    let (modulo_inicial, nombre_inicial) = funcion
        .split_once("::")
        .ok_or_else(|| format!("La función «{funcion}» debe escribirse como modulo::nombre."))?;
    let mut modulos: HashMap<&str, Modulo> = HashMap::new();
    for (nombre, codigo) in FUENTES {
        modulos.insert(nombre, analizar(codigo));
    }
    let principal = modulos
        .get(modulo_inicial)
        .and_then(|modulo| buscar(modulo, nombre_inicial))
        .ok_or_else(|| format!("No encuentro el código de «{funcion}»."))?;
    let mut resultado = Vec::new();
    let mut vistos: HashSet<(String, usize)> = HashSet::new();
    let mut pendientes: VecDeque<(String, usize)> = VecDeque::new();
    pendientes.push_back((modulo_inicial.to_string(), principal));
    while let Some((nombre_modulo, indice)) = pendientes.pop_front() {
        if !vistos.insert((nombre_modulo.clone(), indice)) {
            continue;
        }
        let modulo = &modulos[nombre_modulo.as_str()];
        let definicion = &modulo.definiciones[indice];
        resultado.push(Fragmento {
            nombre: format!("{nombre_modulo}::{}", definicion.nombre),
            fuente: definicion.texto.clone(),
        });
        for (anterior, identificador) in identificadores(&definicion.texto) {
            let destino = match anterior {
                Some(prefijo) if modulos.contains_key(prefijo.as_str()) => Some(prefijo),
                _ => match modulo.importados.get(&identificador) {
                    Some(origen) => Some(origen.clone()),
                    None => Some(nombre_modulo.clone()),
                },
            };
            let Some(destino) = destino else { continue };
            let Some(otro) = modulos.get(destino.as_str()) else {
                continue;
            };
            for encontrado in referencias(otro, &identificador) {
                if !vistos.contains(&(destino.clone(), encontrado)) {
                    pendientes.push_back((destino.clone(), encontrado));
                }
            }
        }
    }
    Ok(resultado)
}

fn buscar(modulo: &Modulo, nombre: &str) -> Option<usize> {
    modulo
        .definiciones
        .iter()
        .position(|definicion| definicion.clase == Clase::Funcion && definicion.nombre == nombre)
}

fn referencias(modulo: &Modulo, identificador: &str) -> Vec<usize> {
    let mut encontrados = Vec::new();
    for (indice, definicion) in modulo.definiciones.iter().enumerate() {
        let coincide = match definicion.clase {
            Clase::Funcion | Clase::Tipo | Clase::Constante => definicion.nombre == identificador,
            Clase::Implementacion => {
                definicion.nombre == identificador
                    || definicion
                        .metodos
                        .iter()
                        .any(|metodo| metodo == identificador)
            }
            Clase::Importacion | Clase::Otra => false,
        };
        if coincide {
            encontrados.push(indice);
        }
    }
    encontrados
}

fn analizar(codigo: &str) -> Modulo {
    let mut definiciones = Vec::new();
    let mut importados = HashMap::new();
    for texto in elementos(codigo) {
        let definicion = describir(&texto);
        if definicion.clase == Clase::Importacion {
            registrar_importacion(&texto, &mut importados);
        }
        definiciones.push(definicion);
    }
    Modulo {
        definiciones,
        importados,
    }
}

fn registrar_importacion(texto: &str, importados: &mut HashMap<String, String>) {
    let Some(resto) = texto
        .trim()
        .trim_start_matches("pub ")
        .strip_prefix("use crate::")
    else {
        return;
    };
    let resto = resto.trim_end_matches(';').trim();
    let Some((modulo, nombres)) = resto.split_once("::") else {
        return;
    };
    let nombres = nombres.trim_start_matches('{').trim_end_matches('}');
    for nombre in nombres.split(',') {
        let nombre = nombre.trim();
        if !nombre.is_empty() && nombre != "*" && nombre != "self" {
            importados.insert(nombre.to_string(), modulo.to_string());
        }
    }
}

fn elementos(codigo: &str) -> Vec<String> {
    let letras: Vec<char> = codigo.chars().collect();
    let mut elementos = Vec::new();
    let mut profundidad = 0usize;
    let mut agrupado = 0usize;
    let mut inicio: Option<usize> = None;
    let mut i = 0;
    while i < letras.len() {
        let letra = letras[i];
        if inicio.is_none() && !letra.is_whitespace() {
            inicio = Some(i);
        }
        if let Some(fin) = saltar_literal(&letras, i) {
            i = fin;
            continue;
        }
        match letra {
            '(' | '[' => agrupado += 1,
            ')' | ']' => agrupado = agrupado.saturating_sub(1),
            '{' => profundidad += 1,
            '}' => {
                profundidad = profundidad.saturating_sub(1);
                if profundidad == 0 {
                    let mut fin = i + 1;
                    if letras.get(fin) == Some(&';') {
                        fin += 1;
                    }
                    if let Some(comienzo) = inicio.take() {
                        elementos.push(letras[comienzo..fin].iter().collect());
                    }
                    i = fin;
                    continue;
                }
            }
            ';' if profundidad == 0 && agrupado == 0 => {
                if let Some(comienzo) = inicio.take() {
                    elementos.push(letras[comienzo..=i].iter().collect());
                }
            }
            _ => {}
        }
        i += 1;
    }
    elementos
}

fn saltar_literal(letras: &[char], i: usize) -> Option<usize> {
    let letra = letras[i];
    if letra == 'r' && (i == 0 || !es_letra_de_identificador(letras[i - 1])) {
        let mut j = i + 1;
        let mut almohadillas = 0;
        while letras.get(j) == Some(&'#') {
            almohadillas += 1;
            j += 1;
        }
        if letras.get(j) == Some(&'"') {
            j += 1;
            while j < letras.len() {
                if letras[j] == '"' && (1..=almohadillas).all(|k| letras.get(j + k) == Some(&'#')) {
                    return Some(j + 1 + almohadillas);
                }
                j += 1;
            }
            return Some(letras.len());
        }
        return None;
    }
    if letra == '"' {
        let mut j = i + 1;
        while j < letras.len() {
            match letras[j] {
                '\\' => j += 2,
                '"' => return Some(j + 1),
                _ => j += 1,
            }
        }
        return Some(letras.len());
    }
    if letra == '\'' {
        if letras.get(i + 1) == Some(&'\\') {
            let mut j = i + 2;
            while j < letras.len() && letras[j] != '\'' {
                j += 1;
            }
            return Some(j + 1);
        }
        if letras.get(i + 2) == Some(&'\'') {
            return Some(i + 3);
        }
    }
    None
}

fn es_letra_de_identificador(letra: char) -> bool {
    letra.is_alphanumeric() || letra == '_'
}

fn describir(texto: &str) -> Definicion {
    let palabras = identificadores(texto);
    let mut clase = Clase::Otra;
    let mut nombre = String::new();
    for (posicion, (_, palabra)) in palabras.iter().enumerate() {
        let siguiente = palabras
            .get(posicion + 1)
            .map(|(_, p)| p.clone())
            .unwrap_or_default();
        match palabra.as_str() {
            "pub" | "crate" | "derive" | "allow" | "cfg" | "test" => continue,
            "fn" => {
                clase = Clase::Funcion;
                nombre = siguiente;
            }
            "struct" | "enum" | "type" | "trait" => {
                clase = Clase::Tipo;
                nombre = siguiente;
            }
            "const" | "static" if siguiente != "fn" && siguiente != "unsafe" => {
                clase = Clase::Constante;
                nombre = siguiente;
            }
            "const" | "static" => continue,
            "use" => {
                clase = Clase::Importacion;
                nombre = siguiente;
            }
            "impl" => {
                clase = Clase::Implementacion;
                nombre = nombre_implementado(texto);
            }
            _ => continue,
        }
        break;
    }
    let metodos = if clase == Clase::Implementacion {
        let mut metodos = Vec::new();
        let lista = identificadores(texto);
        for (posicion, (_, palabra)) in lista.iter().enumerate() {
            if palabra == "fn"
                && let Some((_, metodo)) = lista.get(posicion + 1)
            {
                metodos.push(metodo.clone());
            }
        }
        metodos
    } else {
        Vec::new()
    };
    Definicion {
        nombre,
        clase,
        texto: texto.to_string(),
        metodos,
    }
}

fn nombre_implementado(texto: &str) -> String {
    let cabecera = texto.split('{').next().unwrap_or_default();
    let cabecera = match cabecera.split_once(" for ") {
        Some((_, tipo)) => tipo,
        None => cabecera.trim_start().trim_start_matches("impl"),
    };
    identificadores(cabecera)
        .into_iter()
        .map(|(_, palabra)| palabra)
        .find(|palabra| palabra.chars().next().is_some_and(|c| c.is_uppercase()))
        .unwrap_or_default()
}

fn identificadores(texto: &str) -> Vec<(Option<String>, String)> {
    let letras: Vec<char> = texto.chars().collect();
    let mut resultado: Vec<(Option<String>, String)> = Vec::new();
    let mut i = 0;
    while i < letras.len() {
        if let Some(fin) = saltar_literal(&letras, i) {
            i = fin;
            continue;
        }
        let letra = letras[i];
        if letra.is_alphabetic() || letra == '_' {
            let inicio = i;
            while i < letras.len() && es_letra_de_identificador(letras[i]) {
                i += 1;
            }
            let palabra: String = letras[inicio..i].iter().collect();
            let prefijo = if inicio >= 2 && letras[inicio - 1] == ':' && letras[inicio - 2] == ':' {
                resultado.last().map(|(_, anterior)| anterior.clone())
            } else {
                None
            };
            resultado.push((prefijo, palabra));
            continue;
        }
        i += 1;
    }
    resultado
}
