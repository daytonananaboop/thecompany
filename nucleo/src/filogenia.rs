use crate::azar::Azar;
use crate::secuencia::indice_base;
use rayon::prelude::*;

const MAXIMO_DE_NODOS: usize = 2000;
const INFINITO: u32 = 1_000_000_000;
const MAXIMA_DISTANCIA: f64 = 1e12;
const BASES: &[u8] = b"ACGT";

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Arbol {
    pub nombres: Vec<String>,
    pub vecinos: Vec<Vec<(usize, f64)>>,
}

impl Arbol {
    pub fn agregar_nodo(&mut self, nombre: &str) -> usize {
        self.nombres.push(nombre.to_string());
        self.vecinos.push(Vec::new());
        self.nombres.len() - 1
    }

    pub fn conectar(&mut self, a: usize, b: usize, longitud: f64) {
        self.vecinos[a].push((b, longitud));
        self.vecinos[b].push((a, longitud));
    }

    pub fn desconectar(&mut self, a: usize, b: usize) -> Option<f64> {
        let posicion = self.vecinos[a]
            .iter()
            .position(|&(vecino, _)| vecino == b)?;
        let (_, longitud) = self.vecinos[a].remove(posicion);
        self.vecinos[b].retain(|&(vecino, _)| vecino != a);
        Some(longitud)
    }

    pub fn longitud_entre(&self, a: usize, b: usize) -> Option<f64> {
        self.vecinos[a]
            .iter()
            .find(|&&(vecino, _)| vecino == b)
            .map(|&(_, longitud)| longitud)
    }

    pub fn grado(&self, nodo: usize) -> usize {
        self.vecinos[nodo].len()
    }

    pub fn hojas(&self) -> Vec<usize> {
        (0..self.nombres.len())
            .filter(|&nodo| self.grado(nodo) <= 1)
            .collect()
    }

    pub fn aristas(&self) -> Vec<(usize, usize, f64)> {
        let mut aristas = Vec::new();
        for (nodo, vecinos) in self.vecinos.iter().enumerate() {
            for &(vecino, longitud) in vecinos {
                if nodo < vecino {
                    aristas.push((nodo, vecino, longitud));
                }
            }
        }
        aristas
    }

    pub fn buscar_nodo(&self, nombre: &str) -> Option<usize> {
        self.nombres.iter().position(|otro| otro == nombre)
    }

    pub fn quitar_nodo(&mut self, nodo: usize) {
        while let Some(&(vecino, _)) = self.vecinos[nodo].first() {
            self.desconectar(nodo, vecino);
        }
        let ultimo = self.nombres.len() - 1;
        self.nombres.swap_remove(nodo);
        self.vecinos.swap_remove(nodo);
        if nodo == ultimo {
            return;
        }
        for indice in 0..self.vecinos[nodo].len() {
            let vecino = self.vecinos[nodo][indice].0;
            for arista in self.vecinos[vecino].iter_mut() {
                if arista.0 == ultimo {
                    arista.0 = nodo;
                }
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ArbolConRaiz {
    pub arbol: Arbol,
    pub raiz: usize,
}

pub fn nombre_interno(nodo: usize) -> String {
    format!("n{nodo}")
}

pub fn es_nombre_valido(nombre: &str) -> bool {
    !nombre.is_empty()
        && nombre
            .chars()
            .all(|letra| letra.is_alphanumeric() || letra == '_' || letra == '.')
}

fn casi_igual(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * (1.0 + a.abs().max(b.abs()))
}

pub fn recorrer_desde(arbol: &Arbol, raiz: usize) -> (Vec<usize>, Vec<Option<usize>>) {
    let cantidad = arbol.nombres.len();
    let mut orden = Vec::with_capacity(cantidad);
    let mut padres = vec![None; cantidad];
    let mut visitados = vec![false; cantidad];
    let mut pila = vec![raiz];
    visitados[raiz] = true;
    while let Some(nodo) = pila.pop() {
        orden.push(nodo);
        for &(vecino, _) in &arbol.vecinos[nodo] {
            if !visitados[vecino] {
                visitados[vecino] = true;
                padres[vecino] = Some(nodo);
                pila.push(vecino);
            }
        }
    }
    (orden, padres)
}

fn es_un_arbol(arbol: &Arbol) -> bool {
    let cantidad = arbol.nombres.len();
    cantidad > 0
        && arbol.aristas().len() + 1 == cantidad
        && recorrer_desde(arbol, 0).0.len() == cantidad
}

pub fn leer_aristas(texto: &str) -> Result<Arbol, String> {
    let mut arbol = Arbol::default();
    for linea in texto.lines() {
        let linea = linea.trim();
        if linea.is_empty() {
            continue;
        }
        let (desde, resto) = linea
            .split_once("->")
            .ok_or_else(|| format!("La línea «{linea}» debe tener la forma a->b:longitud."))?;
        let (hasta, longitud) = match resto.split_once(':') {
            Some((hasta, longitud)) => (hasta, leer_longitud(longitud)?),
            None => (resto, 1.0),
        };
        let a = nodo_con_nombre(&mut arbol, desde.trim())?;
        let b = nodo_con_nombre(&mut arbol, hasta.trim())?;
        agregar_arista_sin_repetir(&mut arbol, a, b, longitud)?;
    }
    if arbol.nombres.is_empty() {
        return Err("No hay ninguna arista.".to_string());
    }
    if !es_un_arbol(&arbol) {
        return Err("Las aristas no forman un árbol: hay un ciclo o partes sueltas.".to_string());
    }
    Ok(arbol)
}

fn leer_longitud(texto: &str) -> Result<f64, String> {
    let texto = texto.trim();
    match texto.parse::<f64>() {
        Ok(longitud) if longitud.is_finite() && (0.0..=MAXIMA_DISTANCIA).contains(&longitud) => {
            Ok(longitud)
        }
        _ => Err(format!(
            "«{texto}» no es una longitud válida: debe ser un número entre 0 y 1 000 000 000 000."
        )),
    }
}

fn nodo_con_nombre(arbol: &mut Arbol, nombre: &str) -> Result<usize, String> {
    if !es_nombre_valido(nombre) {
        return Err(format!(
            "«{nombre}» no es un nombre válido: usa letras, números, «_» o «.»."
        ));
    }
    if let Some(nodo) = arbol.buscar_nodo(nombre) {
        return Ok(nodo);
    }
    if arbol.nombres.len() == MAXIMO_DE_NODOS {
        return Err(format!("El árbol tiene más de {MAXIMO_DE_NODOS} nodos."));
    }
    Ok(arbol.agregar_nodo(nombre))
}

fn agregar_arista_sin_repetir(
    arbol: &mut Arbol,
    a: usize,
    b: usize,
    longitud: f64,
) -> Result<(), String> {
    if a == b {
        return Err(format!(
            "La arista de «{}» une el nodo consigo mismo.",
            arbol.nombres[a]
        ));
    }
    match arbol.longitud_entre(a, b) {
        None => {
            arbol.conectar(a, b, longitud);
            Ok(())
        }
        Some(anterior) if casi_igual(anterior, longitud) => Ok(()),
        Some(anterior) => Err(format!(
            "La arista entre «{}» y «{}» aparece con dos longitudes distintas: {anterior} y {longitud}.",
            arbol.nombres[a], arbol.nombres[b]
        )),
    }
}

pub fn leer_newick(texto: &str) -> Result<ArbolConRaiz, String> {
    let texto = texto.trim();
    let texto = texto.strip_suffix(';').unwrap_or(texto);
    if texto.contains(';') {
        return Err("Escribe un solo árbol, terminado en «;».".to_string());
    }
    let letras: Vec<char> = texto.chars().collect();
    let mut arbol = Arbol::default();
    let mut padres: Vec<Option<usize>> = Vec::new();
    let mut longitudes: Vec<f64> = Vec::new();
    let mut abiertos: Vec<usize> = Vec::new();
    let mut ultimo: Option<usize> = None;
    let mut recien_cerrado = false;
    let mut i = 0;
    while i < letras.len() {
        let letra = letras[i];
        if letra.is_whitespace() {
            i += 1;
            continue;
        }
        match letra {
            '(' => {
                if ultimo.is_some() {
                    return Err("Falta una coma antes de «(».".to_string());
                }
                let padre = abiertos.last().copied();
                let nodo = nodo_de_newick(&mut arbol, &mut padres, &mut longitudes, padre, "")?;
                abiertos.push(nodo);
                i += 1;
            }
            ',' | ')' => {
                if ultimo.is_none() {
                    return Err("Hay una hoja sin nombre: revisa comas y paréntesis.".to_string());
                }
                ultimo = None;
                recien_cerrado = false;
                if letra == ')' {
                    let nodo = abiertos.pop().ok_or("Sobra un paréntesis de cierre.")?;
                    ultimo = Some(nodo);
                    recien_cerrado = true;
                }
                i += 1;
            }
            ':' => {
                let nodo = ultimo.ok_or("Falta un nombre antes de «:».")?;
                if !longitudes[nodo].is_nan() {
                    return Err(
                        "Cada nodo puede llevar una sola longitud después de «:».".to_string()
                    );
                }
                let fin = fin_de_palabra(&letras, i + 1);
                let palabra: String = letras[i + 1..fin].iter().collect();
                longitudes[nodo] = leer_longitud(&palabra)?;
                recien_cerrado = false;
                i = fin;
            }
            _ => {
                let fin = fin_de_palabra(&letras, i);
                let nombre: String = letras[i..fin].iter().collect();
                if !es_nombre_valido(&nombre) {
                    return Err(format!(
                        "«{nombre}» no es un nombre válido: usa letras, números, «_» o «.»."
                    ));
                }
                match ultimo {
                    Some(nodo) if recien_cerrado => arbol.nombres[nodo] = nombre,
                    Some(_) => return Err(format!("Falta una coma antes de «{nombre}».")),
                    None => {
                        let padre = abiertos.last().copied();
                        let nodo = nodo_de_newick(
                            &mut arbol,
                            &mut padres,
                            &mut longitudes,
                            padre,
                            &nombre,
                        )?;
                        ultimo = Some(nodo);
                    }
                }
                recien_cerrado = false;
                i = fin;
            }
        }
    }
    if !abiertos.is_empty() {
        return Err("Falta cerrar un paréntesis.".to_string());
    }
    if ultimo.is_none() {
        return Err("El árbol está vacío o termina en una coma.".to_string());
    }
    let raices: Vec<usize> = (0..padres.len())
        .filter(|&nodo| padres[nodo].is_none())
        .collect();
    if raices.len() != 1 {
        return Err(
            "Hay varios árboles sueltos: encierra todo en un par de paréntesis.".to_string(),
        );
    }
    revisar_nombres_repetidos(&arbol)?;
    for (nodo, padre) in padres.iter().enumerate() {
        if let Some(padre) = *padre {
            let longitud = if longitudes[nodo].is_nan() {
                1.0
            } else {
                longitudes[nodo]
            };
            arbol.conectar(padre, nodo, longitud);
        }
    }
    Ok(ArbolConRaiz {
        arbol,
        raiz: raices[0],
    })
}

fn nodo_de_newick(
    arbol: &mut Arbol,
    padres: &mut Vec<Option<usize>>,
    longitudes: &mut Vec<f64>,
    padre: Option<usize>,
    nombre: &str,
) -> Result<usize, String> {
    if arbol.nombres.len() == MAXIMO_DE_NODOS {
        return Err(format!("El árbol tiene más de {MAXIMO_DE_NODOS} nodos."));
    }
    padres.push(padre);
    longitudes.push(f64::NAN);
    Ok(arbol.agregar_nodo(nombre))
}

fn fin_de_palabra(letras: &[char], desde: usize) -> usize {
    let mut fin = desde;
    while fin < letras.len() && !letras[fin].is_whitespace() && !"(),:;".contains(letras[fin]) {
        fin += 1;
    }
    fin
}

fn revisar_nombres_repetidos(arbol: &Arbol) -> Result<(), String> {
    let mut vistos = std::collections::HashSet::new();
    for nombre in &arbol.nombres {
        if !nombre.is_empty() && !vistos.insert(nombre) {
            return Err(format!(
                "El nombre «{nombre}» aparece dos veces en el árbol."
            ));
        }
    }
    Ok(())
}

pub fn escribir_newick(arbol: &Arbol, raiz: usize) -> String {
    let mut texto = String::new();
    let mut visitados = vec![false; arbol.nombres.len()];
    if raiz < arbol.nombres.len() {
        escribir_subarbol(arbol, raiz, &mut visitados, &mut texto);
    }
    texto.push(';');
    texto
}

fn escribir_subarbol(arbol: &Arbol, nodo: usize, visitados: &mut [bool], texto: &mut String) {
    visitados[nodo] = true;
    let hijos: Vec<(usize, f64)> = arbol.vecinos[nodo]
        .iter()
        .filter(|&&(vecino, _)| !visitados[vecino])
        .copied()
        .collect();
    if !hijos.is_empty() {
        texto.push('(');
        for (posicion, &(hijo, longitud)) in hijos.iter().enumerate() {
            if posicion > 0 {
                texto.push(',');
            }
            escribir_subarbol(arbol, hijo, visitados, texto);
            texto.push(':');
            texto.push_str(&longitud_en_texto(longitud));
        }
        texto.push(')');
    }
    texto.push_str(&arbol.nombres[nodo]);
}

pub fn longitud_en_texto(longitud: f64) -> String {
    let texto = format!("{longitud:.6}");
    let texto = texto.trim_end_matches('0').trim_end_matches('.');
    if texto == "-0" {
        return "0".to_string();
    }
    texto.to_string()
}

pub fn distancias_desde(arbol: &Arbol, origen: usize) -> Vec<f64> {
    let mut distancias = vec![f64::INFINITY; arbol.nombres.len()];
    let mut visitados = vec![false; arbol.nombres.len()];
    distancias[origen] = 0.0;
    visitados[origen] = true;
    let mut pila = vec![origen];
    while let Some(nodo) = pila.pop() {
        for &(vecino, longitud) in &arbol.vecinos[nodo] {
            if !visitados[vecino] {
                visitados[vecino] = true;
                distancias[vecino] = distancias[nodo] + longitud;
                pila.push(vecino);
            }
        }
    }
    distancias
}

pub fn distancias_entre_hojas(arbol: &Arbol, hojas: &[usize]) -> Vec<Vec<f64>> {
    let mut matriz = Vec::with_capacity(hojas.len());
    for &hoja in hojas {
        let desde_la_hoja = distancias_desde(arbol, hoja);
        let fila: Vec<f64> = hojas.iter().map(|&otra| desde_la_hoja[otra]).collect();
        matriz.push(fila);
    }
    matriz
}

pub fn validar_distancias(distancias: &[Vec<f64>]) -> Result<(), String> {
    let cantidad = distancias.len();
    if cantidad == 0 {
        return Err("La matriz de distancias está vacía.".to_string());
    }
    for (fila, valores) in distancias.iter().enumerate() {
        if valores.len() != cantidad {
            return Err(format!(
                "La matriz debe ser cuadrada: la fila {fila} (contando desde 0) tiene {} números y debería tener {cantidad}.",
                valores.len()
            ));
        }
    }
    for (fila, valores) in distancias.iter().enumerate() {
        for (columna, &valor) in valores.iter().enumerate() {
            if !valor.is_finite() || !(0.0..=MAXIMA_DISTANCIA).contains(&valor) {
                return Err(format!(
                    "La distancia de la fila {fila}, columna {columna} (contando desde 0) debe ser un número entre 0 y 1 000 000 000 000."
                ));
            }
            if fila == columna && valor != 0.0 {
                return Err(format!(
                    "La diagonal debe tener ceros, pero en la fila {fila} (contando desde 0) hay {valor}."
                ));
            }
            if !casi_igual(valor, distancias[columna][fila]) {
                return Err(format!(
                    "La matriz debe ser simétrica: la fila {fila}, columna {columna} vale {valor} y la fila {columna}, columna {fila} vale {} (contando desde 0).",
                    distancias[columna][fila]
                ));
            }
        }
    }
    Ok(())
}

fn revisar_nombres(cantidad: usize, nombres: &[String]) -> Result<(), String> {
    if nombres.len() != cantidad {
        return Err(format!(
            "Hay {} nombres para {cantidad} filas de la matriz.",
            nombres.len()
        ));
    }
    Ok(())
}

pub fn longitud_de_rama(distancias: &[Vec<f64>], hoja: usize) -> Result<f64, String> {
    longitud_de_rama_con_par(distancias, hoja).map(|(longitud, _, _)| longitud)
}

pub fn longitud_de_rama_con_par(
    distancias: &[Vec<f64>],
    hoja: usize,
) -> Result<(f64, usize, usize), String> {
    validar_distancias(distancias)?;
    let cantidad = distancias.len();
    if cantidad < 2 {
        return Err("Hacen falta al menos dos hojas para hablar de una rama.".to_string());
    }
    if hoja >= cantidad {
        return Err(format!(
            "La hoja {hoja} no existe: la matriz tiene {cantidad} hojas, numeradas desde 0."
        ));
    }
    Ok(rama_entre_las_primeras(distancias, hoja, cantidad))
}

pub fn rama_entre_las_primeras(
    distancias: &[Vec<f64>],
    hoja: usize,
    cantidad: usize,
) -> (f64, usize, usize) {
    let mut mejor = (f64::INFINITY, hoja, hoja);
    for i in (0..cantidad).filter(|&i| i != hoja) {
        for k in (i..cantidad).filter(|&k| k != hoja) {
            let candidata = (distancias[i][hoja] + distancias[hoja][k] - distancias[i][k]) / 2.0;
            if candidata < mejor.0 {
                mejor = (candidata, i, k);
            }
        }
    }
    mejor
}

pub fn es_aditiva(distancias: &[Vec<f64>]) -> bool {
    let sin_nombres = vec![String::new(); distancias.len()];
    validar_distancias(distancias).is_ok()
        && reproduce_las_distancias(&armar_arbol_aditivo(distancias, &sin_nombres), distancias)
}

pub fn trio_no_metrico(distancias: &[Vec<f64>]) -> Option<(usize, usize, usize)> {
    let cantidad = distancias.len();
    for i in 0..cantidad {
        for j in (i + 1)..cantidad {
            for k in (j + 1)..cantidad {
                let (a, b, c) = (distancias[i][j], distancias[i][k], distancias[j][k]);
                let mayor = a.max(b).max(c);
                let otros_dos = a + b + c - mayor;
                if mayor > otros_dos && !casi_igual(mayor, otros_dos) {
                    return Some((i, j, k));
                }
            }
        }
    }
    None
}

pub fn cuarteto_no_aditivo(distancias: &[Vec<f64>]) -> Option<(usize, usize, usize, usize)> {
    let cantidad = distancias.len();
    for i in 0..cantidad {
        for j in (i + 1)..cantidad {
            for k in (j + 1)..cantidad {
                for l in (k + 1)..cantidad {
                    let (primera, segunda, tercera) = sumas_de_cuarteto(distancias, i, j, k, l);
                    if !dos_mayores_iguales(primera, segunda, tercera) {
                        return Some((i, j, k, l));
                    }
                }
            }
        }
    }
    None
}

pub fn sumas_de_cuarteto(
    distancias: &[Vec<f64>],
    i: usize,
    j: usize,
    k: usize,
    l: usize,
) -> (f64, f64, f64) {
    (
        distancias[i][j] + distancias[k][l],
        distancias[i][k] + distancias[j][l],
        distancias[i][l] + distancias[j][k],
    )
}

fn dos_mayores_iguales(a: f64, b: f64, c: f64) -> bool {
    let mayor = a.max(b).max(c);
    let menor = a.min(b).min(c);
    let del_medio = a + b + c - mayor - menor;
    casi_igual(del_medio, mayor)
}

pub fn trio_no_ultrametrico(distancias: &[Vec<f64>]) -> Option<(usize, usize, usize)> {
    let cantidad = distancias.len();
    for i in 0..cantidad {
        for j in (i + 1)..cantidad {
            for k in (j + 1)..cantidad {
                if !dos_mayores_iguales(distancias[i][j], distancias[i][k], distancias[j][k]) {
                    return Some((i, j, k));
                }
            }
        }
    }
    None
}

struct HojaPendiente {
    hoja: usize,
    rama: f64,
    desde: usize,
    hacia: usize,
    distancia: f64,
}

pub fn filogenia_aditiva(distancias: &[Vec<f64>], nombres: &[String]) -> Result<Arbol, String> {
    validar_distancias(distancias)?;
    revisar_nombres(distancias.len(), nombres)?;
    let arbol = armar_arbol_aditivo(distancias, nombres);
    if !reproduce_las_distancias(&arbol, distancias) {
        return Err(
            "La matriz no es aditiva: ningún árbol reproduce exactamente todas sus distancias."
                .to_string(),
        );
    }
    Ok(arbol)
}

fn armar_arbol_aditivo(distancias: &[Vec<f64>], nombres: &[String]) -> Arbol {
    let cantidad = distancias.len();
    let mut arbol = Arbol::default();
    for nombre in nombres {
        arbol.agregar_nodo(nombre);
    }
    if cantidad < 2 {
        return arbol;
    }
    let mut recortada = distancias.to_vec();
    let mut pendientes = Vec::new();
    for hoja in (2..cantidad).rev() {
        let (rama, _, _) = rama_entre_las_primeras(&recortada, hoja, hoja + 1);
        for fila in recortada.iter_mut().take(hoja) {
            fila[hoja] -= rama;
        }
        for distancia in recortada[hoja].iter_mut().take(hoja) {
            *distancia -= rama;
        }
        let (desde, hacia) = hojas_que_rodean(&recortada, hoja);
        pendientes.push(HojaPendiente {
            hoja,
            rama,
            desde,
            hacia,
            distancia: recortada[desde][hoja],
        });
    }
    arbol.conectar(0, 1, distancias[0][1]);
    for pendiente in pendientes.iter().rev() {
        let punto = punto_en_el_camino(
            &mut arbol,
            pendiente.desde,
            pendiente.hacia,
            pendiente.distancia,
        );
        arbol.conectar(punto, pendiente.hoja, pendiente.rama);
    }
    arbol
}

pub fn reproduce_las_distancias(arbol: &Arbol, distancias: &[Vec<f64>]) -> bool {
    let sin_ramas_negativas = arbol
        .aristas()
        .iter()
        .all(|&(_, _, longitud)| longitud >= 0.0 || casi_igual(longitud, 0.0));
    let hojas: Vec<usize> = (0..distancias.len()).collect();
    let del_arbol = distancias_entre_hojas(arbol, &hojas);
    let iguales = del_arbol
        .iter()
        .flatten()
        .zip(distancias.iter().flatten())
        .all(|(&en_el_arbol, &en_la_matriz)| casi_igual(en_el_arbol, en_la_matriz));
    sin_ramas_negativas && iguales
}

pub fn hojas_que_rodean(distancias: &[Vec<f64>], hoja: usize) -> (usize, usize) {
    let mut mejor = (0, 1);
    let mut menor_diferencia = f64::INFINITY;
    for i in 0..hoja {
        for k in (i + 1)..hoja {
            let diferencia = (distancias[i][k] - distancias[i][hoja] - distancias[hoja][k]).abs();
            if diferencia < menor_diferencia {
                menor_diferencia = diferencia;
                mejor = (i, k);
            }
        }
    }
    mejor
}

pub fn punto_en_el_camino(arbol: &mut Arbol, desde: usize, hacia: usize, distancia: f64) -> usize {
    let camino = camino_entre(arbol, desde, hacia);
    let mut recorrido = 0.0;
    for par in camino.windows(2) {
        let (u, v) = (par[0], par[1]);
        if casi_igual(recorrido, distancia) {
            return u;
        }
        let largo = arbol.longitud_entre(u, v).unwrap_or(0.0);
        if recorrido + largo > distancia && !casi_igual(recorrido + largo, distancia) {
            let nuevo = arbol.agregar_nodo(&nombre_interno(arbol.nombres.len()));
            arbol.desconectar(u, v);
            arbol.conectar(u, nuevo, distancia - recorrido);
            arbol.conectar(nuevo, v, recorrido + largo - distancia);
            return nuevo;
        }
        recorrido += largo;
    }
    hacia
}

fn camino_entre(arbol: &Arbol, desde: usize, hacia: usize) -> Vec<usize> {
    let (_, padres) = recorrer_desde(arbol, desde);
    let mut camino = vec![hacia];
    let mut actual = hacia;
    while let Some(padre) = padres[actual] {
        camino.push(padre);
        actual = padre;
    }
    camino.reverse();
    camino
}

pub fn reemplazar_par(distancias: &mut Vec<Vec<f64>>, i: usize, j: usize, fila_nueva: Vec<f64>) {
    distancias.remove(j);
    distancias.remove(i);
    for (fila, &valor) in distancias.iter_mut().zip(&fila_nueva) {
        fila.remove(j);
        fila.remove(i);
        fila.push(valor);
    }
    let mut ultima = fila_nueva;
    ultima.push(0.0);
    distancias.push(ultima);
}

pub fn quitar_par(lista: &mut Vec<usize>, i: usize, j: usize) {
    lista.remove(j);
    lista.remove(i);
}

pub fn par_mas_cercano(distancias: &[Vec<f64>]) -> (usize, usize) {
    let mut mejor = (0, 1);
    let mut menor = f64::INFINITY;
    for (i, fila) in distancias.iter().enumerate() {
        for (j, &distancia) in fila.iter().enumerate().skip(i + 1) {
            if distancia < menor {
                menor = distancia;
                mejor = (i, j);
            }
        }
    }
    mejor
}

pub fn upgma(distancias: &[Vec<f64>], nombres: &[String]) -> Result<ArbolConRaiz, String> {
    validar_distancias(distancias)?;
    revisar_nombres(distancias.len(), nombres)?;
    let mut arbol = Arbol::default();
    let mut grupos: Vec<usize> = nombres
        .iter()
        .map(|nombre| arbol.agregar_nodo(nombre))
        .collect();
    let mut hojas_por_grupo = vec![1; nombres.len()];
    let mut edades = vec![0.0; nombres.len()];
    let mut actual = distancias.to_vec();
    while actual.len() > 1 {
        let (i, j) = par_mas_cercano(&actual);
        let nuevo = arbol.agregar_nodo(&nombre_interno(arbol.nombres.len()));
        let edad = actual[i][j] / 2.0;
        edades.push(edad);
        arbol.conectar(nuevo, grupos[i], edad - edades[grupos[i]]);
        arbol.conectar(nuevo, grupos[j], edad - edades[grupos[j]]);
        let fila_nueva = promedios_hacia_el_grupo(&actual, &hojas_por_grupo, i, j);
        reemplazar_par(&mut actual, i, j, fila_nueva);
        let hojas_del_nuevo = hojas_por_grupo[i] + hojas_por_grupo[j];
        quitar_par(&mut hojas_por_grupo, i, j);
        hojas_por_grupo.push(hojas_del_nuevo);
        quitar_par(&mut grupos, i, j);
        grupos.push(nuevo);
    }
    Ok(ArbolConRaiz {
        arbol,
        raiz: grupos[0],
    })
}

pub fn promedios_hacia_el_grupo(
    distancias: &[Vec<f64>],
    hojas_por_grupo: &[usize],
    i: usize,
    j: usize,
) -> Vec<f64> {
    let peso_i = hojas_por_grupo[i] as f64;
    let peso_j = hojas_por_grupo[j] as f64;
    let mut fila = Vec::with_capacity(distancias.len());
    for (k, valores) in distancias.iter().enumerate() {
        if k != i && k != j {
            fila.push((valores[i] * peso_i + valores[j] * peso_j) / (peso_i + peso_j));
        }
    }
    fila
}

#[derive(Debug, Clone, PartialEq)]
pub struct EstadoDeUnion {
    pub distancias: Vec<Vec<f64>>,
    pub nodos: Vec<usize>,
    pub arbol: Arbol,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ParUnido {
    pub primero: usize,
    pub segundo: usize,
    pub nuevo: usize,
    pub rama_primero: f64,
    pub rama_segundo: f64,
}

pub fn union_de_vecinos(distancias: &[Vec<f64>], nombres: &[String]) -> Result<Arbol, String> {
    let mut estado = empezar_union(distancias, nombres)?;
    while estado.distancias.len() > 2 {
        unir_par_de_vecinos(&mut estado);
    }
    unir_ultimo_par(&mut estado);
    Ok(estado.arbol)
}

pub fn empezar_union(distancias: &[Vec<f64>], nombres: &[String]) -> Result<EstadoDeUnion, String> {
    validar_distancias(distancias)?;
    revisar_nombres(distancias.len(), nombres)?;
    let mut arbol = Arbol::default();
    let nodos = nombres
        .iter()
        .map(|nombre| arbol.agregar_nodo(nombre))
        .collect();
    Ok(EstadoDeUnion {
        distancias: distancias.to_vec(),
        nodos,
        arbol,
    })
}

pub fn unir_par_de_vecinos(estado: &mut EstadoDeUnion) -> Option<ParUnido> {
    let n = estado.distancias.len();
    if n < 3 {
        return None;
    }
    let totales = distancias_totales(&estado.distancias);
    let (i, j) = par_con_menor_valor_de_union(&estado.distancias, &totales);
    let delta = (totales[i] - totales[j]) / (n - 2) as f64;
    let rama_i = (estado.distancias[i][j] + delta) / 2.0;
    let rama_j = (estado.distancias[i][j] - delta) / 2.0;
    let nuevo = estado
        .arbol
        .agregar_nodo(&nombre_interno(estado.arbol.nombres.len()));
    estado.arbol.conectar(nuevo, estado.nodos[i], rama_i);
    estado.arbol.conectar(nuevo, estado.nodos[j], rama_j);
    let fila_nueva = distancias_al_nodo_nuevo(&estado.distancias, i, j);
    reemplazar_par(&mut estado.distancias, i, j, fila_nueva);
    let par = ParUnido {
        primero: estado.nodos[i],
        segundo: estado.nodos[j],
        nuevo,
        rama_primero: rama_i,
        rama_segundo: rama_j,
    };
    quitar_par(&mut estado.nodos, i, j);
    estado.nodos.push(nuevo);
    Some(par)
}

pub fn unir_ultimo_par(estado: &mut EstadoDeUnion) -> Option<(usize, usize, f64)> {
    if estado.distancias.len() != 2 {
        return None;
    }
    let longitud = estado.distancias[0][1];
    estado
        .arbol
        .conectar(estado.nodos[0], estado.nodos[1], longitud);
    Some((estado.nodos[0], estado.nodos[1], longitud))
}

pub fn distancias_totales(distancias: &[Vec<f64>]) -> Vec<f64> {
    distancias.iter().map(|fila| fila.iter().sum()).collect()
}

pub fn valor_de_union(distancias: &[Vec<f64>], totales: &[f64], i: usize, j: usize) -> f64 {
    let n = distancias.len() as f64;
    (n - 2.0) * distancias[i][j] - totales[i] - totales[j]
}

pub fn matriz_de_union(distancias: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let totales = distancias_totales(distancias);
    let n = distancias.len();
    let mut matriz = vec![vec![0.0; n]; n];
    for (i, fila) in matriz.iter_mut().enumerate() {
        for (j, valor) in fila.iter_mut().enumerate() {
            if i != j {
                *valor = valor_de_union(distancias, &totales, i, j);
            }
        }
    }
    matriz
}

fn par_con_menor_valor_de_union(distancias: &[Vec<f64>], totales: &[f64]) -> (usize, usize) {
    let n = distancias.len();
    let mut mejor = (0, 1);
    let mut menor = f64::INFINITY;
    for i in 0..n {
        for j in (i + 1)..n {
            let valor = valor_de_union(distancias, totales, i, j);
            if valor < menor {
                menor = valor;
                mejor = (i, j);
            }
        }
    }
    mejor
}

fn distancias_al_nodo_nuevo(distancias: &[Vec<f64>], i: usize, j: usize) -> Vec<f64> {
    let mut fila = Vec::with_capacity(distancias.len());
    for (k, valores) in distancias.iter().enumerate() {
        if k != i && k != j {
            fila.push((valores[i] + valores[j] - distancias[i][j]) / 2.0);
        }
    }
    fila
}

fn es_base(letra: u8) -> bool {
    matches!(letra, b'A' | b'C' | b'G' | b'T')
}

pub fn contar_diferencias(primera: &[u8], segunda: &[u8]) -> (usize, usize) {
    let mut comparadas = 0;
    let mut diferentes = 0;
    for (&a, &b) in primera.iter().zip(segunda) {
        if es_base(a) && es_base(b) {
            comparadas += 1;
            if a != b {
                diferentes += 1;
            }
        }
    }
    (comparadas, diferentes)
}

pub fn distancia_p(primera: &[u8], segunda: &[u8]) -> Option<f64> {
    let (comparadas, diferentes) = contar_diferencias(primera, segunda);
    if comparadas == 0 {
        return None;
    }
    Some(diferentes as f64 / comparadas as f64)
}

pub fn distancia_jukes_cantor(p: f64) -> f64 {
    if p >= 0.75 {
        return f64::INFINITY;
    }
    if p <= 0.0 {
        return 0.0;
    }
    -0.75 * (1.0 - 4.0 * p / 3.0).ln()
}

fn revisar_alineamiento(secuencias: &[Vec<u8>]) -> Result<usize, String> {
    let Some(primera) = secuencias.first() else {
        return Err("No hay secuencias.".to_string());
    };
    let longitud = primera.len();
    if longitud == 0 {
        return Err("Las secuencias están vacías.".to_string());
    }
    for (indice, secuencia) in secuencias.iter().enumerate() {
        if secuencia.len() != longitud {
            return Err(format!(
                "Las secuencias deben estar alineadas y medir lo mismo: la {indice} (contando desde 0) mide {} y la primera mide {longitud}.",
                secuencia.len()
            ));
        }
    }
    Ok(longitud)
}

pub fn matriz_de_distancias_p(secuencias: &[Vec<u8>]) -> Result<Vec<Vec<f64>>, String> {
    revisar_alineamiento(secuencias)?;
    let cantidad = secuencias.len();
    let filas_de_arriba: Vec<Vec<Option<f64>>> = (0..cantidad)
        .into_par_iter()
        .map(|i| {
            ((i + 1)..cantidad)
                .map(|j| distancia_p(&secuencias[i], &secuencias[j]))
                .collect()
        })
        .collect();
    let mut matriz = vec![vec![0.0; cantidad]; cantidad];
    for (i, fila) in filas_de_arriba.iter().enumerate() {
        for (j, &p) in ((i + 1)..cantidad).zip(fila) {
            let p = p.ok_or_else(|| {
                format!(
                    "Las secuencias {i} y {j} (contando desde 0) no comparten ninguna posición sin hueco."
                )
            })?;
            matriz[i][j] = p;
            matriz[j][i] = p;
        }
    }
    Ok(matriz)
}

pub fn corregir_con_jukes_cantor(distancias_p: &[Vec<f64>]) -> Vec<Vec<f64>> {
    distancias_p
        .iter()
        .map(|fila| fila.iter().map(|&p| distancia_jukes_cantor(p)).collect())
        .collect()
}

pub fn matriz_jukes_cantor(secuencias: &[Vec<u8>]) -> Result<Vec<Vec<f64>>, String> {
    let distancias_p = matriz_de_distancias_p(secuencias)?;
    Ok(corregir_con_jukes_cantor(&distancias_p))
}

pub fn par_saturado(distancias: &[Vec<f64>]) -> Option<(usize, usize)> {
    for (i, fila) in distancias.iter().enumerate() {
        for (j, distancia) in fila.iter().enumerate() {
            if distancia.is_infinite() {
                return Some((i, j));
            }
        }
    }
    None
}

pub fn arbol_desde_secuencias(secuencias: &[Vec<u8>], nombres: &[String]) -> Result<Arbol, String> {
    revisar_nombres(secuencias.len(), nombres)?;
    let distancias = matriz_jukes_cantor(secuencias)?;
    if let Some((i, j)) = par_saturado(&distancias) {
        return Err(format!(
            "«{}» y «{}» difieren en el 75 % o más de sus posiciones: la distancia de Jukes–Cantor es infinita y no se puede armar el árbol.",
            nombres[i], nombres[j]
        ));
    }
    union_de_vecinos(&distancias, nombres)
}

#[derive(Debug, Clone, PartialEq)]
pub struct Parsimonia {
    pub puntaje: usize,
    pub secuencias: Vec<Vec<u8>>,
    pub cambios_por_posicion: Vec<usize>,
    pub arbol: Arbol,
}

pub fn parsimonia_pequena(
    arbol: &Arbol,
    raiz: usize,
    secuencias: &[Vec<u8>],
) -> Result<Parsimonia, String> {
    let longitud = revisar_hojas_con_secuencia(arbol, raiz, secuencias)?;
    let cantidad = arbol.nombres.len();
    let (orden, padres) = recorrer_desde(arbol, raiz);
    let mut costos = vec![Costos { por_letra: [0; 4] }; cantidad];
    let mut letras = vec![0usize; cantidad];
    let mut inferidas = vec![vec![b'A'; longitud]; cantidad];
    let mut cambios_por_posicion = Vec::with_capacity(longitud);
    for posicion in 0..longitud {
        for &nodo in orden.iter().rev() {
            costos[nodo] = costos_del_nodo(arbol, &costos, &padres, secuencias, nodo, posicion);
        }
        letras[raiz] = letra_mas_barata(&costos[raiz]);
        cambios_por_posicion.push(costos[raiz].por_letra[letras[raiz]] as usize);
        for &nodo in &orden {
            if let Some(padre) = padres[nodo] {
                letras[nodo] = letra_para_el_hijo(&costos[nodo], letras[padre]);
            }
        }
        anotar_letras(&mut inferidas, &letras, posicion);
    }
    let puntaje = cambios_por_posicion.iter().sum();
    let con_cambios = arbol_con_cambios(arbol, &inferidas);
    Ok(Parsimonia {
        puntaje,
        secuencias: inferidas,
        cambios_por_posicion,
        arbol: con_cambios,
    })
}

fn anotar_letras(inferidas: &mut [Vec<u8>], letras: &[usize], posicion: usize) {
    for (secuencia, &letra) in inferidas.iter_mut().zip(letras) {
        secuencia[posicion] = BASES[letra];
    }
}

#[derive(Debug, Clone, Copy)]
struct Costos {
    por_letra: [u32; 4],
}

fn costos_del_nodo(
    arbol: &Arbol,
    costos: &[Costos],
    padres: &[Option<usize>],
    secuencias: &[Vec<u8>],
    nodo: usize,
    posicion: usize,
) -> Costos {
    let mut suma = match secuencias[nodo]
        .get(posicion)
        .and_then(|&letra| indice_base(letra))
    {
        Some(propia) => {
            let mut fija = [INFINITO; 4];
            fija[propia] = 0;
            fija
        }
        None => [0; 4],
    };
    for &(hijo, _) in &arbol.vecinos[nodo] {
        if padres[nodo] == Some(hijo) {
            continue;
        }
        let del_hijo = costos[hijo].por_letra;
        let mas_barato = del_hijo[letra_mas_barata(&costos[hijo])];
        for (letra, total) in suma.iter_mut().enumerate() {
            *total += del_hijo[letra].min(mas_barato + 1);
        }
    }
    Costos { por_letra: suma }
}

fn letra_mas_barata(costos: &Costos) -> usize {
    let mut mejor = 0;
    for letra in 1..4 {
        if costos.por_letra[letra] < costos.por_letra[mejor] {
            mejor = letra;
        }
    }
    mejor
}

fn letra_para_el_hijo(costos_del_hijo: &Costos, letra_del_padre: usize) -> usize {
    let por_letra = costos_del_hijo.por_letra;
    let mas_barata = letra_mas_barata(costos_del_hijo);
    if por_letra[letra_del_padre] <= por_letra[mas_barata] + 1 {
        letra_del_padre
    } else {
        mas_barata
    }
}

pub fn distancia_de_hamming(primera: &[u8], segunda: &[u8]) -> usize {
    primera.iter().zip(segunda).filter(|(a, b)| a != b).count()
}

fn arbol_con_cambios(arbol: &Arbol, secuencias: &[Vec<u8>]) -> Arbol {
    let mut con_cambios = arbol.clone();
    for (nodo, vecinos) in con_cambios.vecinos.iter_mut().enumerate() {
        for (vecino, longitud) in vecinos.iter_mut() {
            *longitud = distancia_de_hamming(&secuencias[nodo], &secuencias[*vecino]) as f64;
        }
    }
    con_cambios
}

fn revisar_hojas_con_secuencia(
    arbol: &Arbol,
    raiz: usize,
    secuencias: &[Vec<u8>],
) -> Result<usize, String> {
    let cantidad = arbol.nombres.len();
    if raiz >= cantidad {
        return Err("El árbol está vacío o la raíz no existe.".to_string());
    }
    if secuencias.len() != cantidad {
        return Err(format!(
            "Hay {} secuencias para {cantidad} nodos.",
            secuencias.len()
        ));
    }
    if !es_un_arbol(arbol) {
        return Err("Las aristas no forman un árbol: hay un ciclo o partes sueltas.".to_string());
    }
    let mut longitud = None;
    for (nodo, secuencia) in secuencias.iter().enumerate() {
        let nombre = &arbol.nombres[nodo];
        if secuencia.is_empty() {
            let vecinos_que_no_son_padre = if nodo == raiz { 0 } else { 1 };
            if arbol.grado(nodo) <= vecinos_que_no_son_padre {
                return Err(format!("La hoja «{nombre}» no tiene secuencia."));
            }
            continue;
        }
        if let Some(letra) = secuencia.iter().find(|&&letra| !es_base(letra)) {
            return Err(format!(
                "La secuencia de «{nombre}» tiene «{}»: aquí solo valen A, C, G y T.",
                *letra as char
            ));
        }
        match longitud {
            None => longitud = Some(secuencia.len()),
            Some(anterior) if anterior != secuencia.len() => {
                return Err(format!(
                    "Todas las secuencias deben medir lo mismo: la de «{nombre}» mide {} y otra mide {anterior}.",
                    secuencia.len()
                ));
            }
            Some(_) => {}
        }
    }
    longitud.ok_or_else(|| "Ninguna hoja tiene secuencia.".to_string())
}

pub fn parsimonia_sin_raiz(arbol: &Arbol, secuencias: &[Vec<u8>]) -> Result<Parsimonia, String> {
    for hoja in arbol.hojas() {
        if secuencias
            .get(hoja)
            .is_some_and(|secuencia| secuencia.is_empty())
        {
            return Err(format!(
                "La hoja «{}» no tiene secuencia.",
                arbol.nombres[hoja]
            ));
        }
    }
    let raiz = (0..arbol.nombres.len())
        .find(|&nodo| arbol.grado(nodo) >= 2)
        .unwrap_or(0);
    parsimonia_pequena(arbol, raiz, secuencias)
}

pub fn quitar_raiz(enraizado: &ArbolConRaiz) -> Arbol {
    let mut arbol = enraizado.arbol.clone();
    let raiz = enraizado.raiz;
    if raiz >= arbol.nombres.len() || arbol.grado(raiz) != 2 {
        return arbol;
    }
    let (a, largo_a) = arbol.vecinos[raiz][0];
    let (b, largo_b) = arbol.vecinos[raiz][1];
    arbol.conectar(a, b, largo_a + largo_b);
    arbol.quitar_nodo(raiz);
    arbol
}

fn otros_dos_vecinos(arbol: &Arbol, nodo: usize, excluido: usize) -> (usize, usize) {
    let mut otros = [nodo; 2];
    let mut encontrados = 0;
    for &(vecino, _) in &arbol.vecinos[nodo] {
        if vecino != excluido && encontrados < 2 {
            otros[encontrados] = vecino;
            encontrados += 1;
        }
    }
    (otros[0], otros[1])
}

pub fn intercambiar(arbol: &mut Arbol, a: usize, x: usize, b: usize, y: usize) {
    let largo_x = arbol.desconectar(a, x).unwrap_or(1.0);
    let largo_y = arbol.desconectar(b, y).unwrap_or(1.0);
    arbol.conectar(b, x, largo_x);
    arbol.conectar(a, y, largo_y);
}

fn revisar_arista_interna(arbol: &Arbol, a: usize, b: usize) -> Result<(), String> {
    let cantidad = arbol.nombres.len();
    if a >= cantidad || b >= cantidad || a == b || arbol.longitud_entre(a, b).is_none() {
        return Err("Esos dos nodos no están unidos por una arista.".to_string());
    }
    if arbol.grado(a) != 3 || arbol.grado(b) != 3 {
        return Err(format!(
            "La arista «{}»–«{}» no es interna de un árbol binario: cada extremo debe tener exactamente tres vecinos.",
            arbol.nombres[a], arbol.nombres[b]
        ));
    }
    Ok(())
}

pub fn vecinos_cercanos(arbol: &Arbol, a: usize, b: usize) -> Result<(Arbol, Arbol), String> {
    revisar_arista_interna(arbol, a, b)?;
    let (_, x) = otros_dos_vecinos(arbol, a, b);
    let (y, z) = otros_dos_vecinos(arbol, b, a);
    let mut primero = arbol.clone();
    intercambiar(&mut primero, a, x, b, y);
    let mut segundo = arbol.clone();
    intercambiar(&mut segundo, a, x, b, z);
    Ok((primero, segundo))
}

#[derive(Debug, Clone, PartialEq)]
pub struct Movimiento {
    pub arista: (usize, usize),
    pub intercambio: (usize, usize),
    pub puntaje: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BusquedaDeParsimonia {
    pub arbol: Arbol,
    pub parsimonia: Parsimonia,
    pub puntaje_inicial: usize,
    pub movimientos: Vec<Movimiento>,
    pub limite_alcanzado: bool,
}

pub fn parsimonia_grande(
    inicial: &Arbol,
    secuencias: &[Vec<u8>],
    maximo_de_movimientos: usize,
) -> Result<BusquedaDeParsimonia, String> {
    let longitud = revisar_arbol_binario(inicial, secuencias)?;
    let conjuntos = conjuntos_de_letras(secuencias, longitud);
    let mut arbol = inicial.clone();
    let puntaje_inicial = puntaje_de_fitch(&arbol, &conjuntos);
    let mut puntaje = puntaje_inicial;
    let mut movimientos = Vec::new();
    while movimientos.len() < maximo_de_movimientos {
        let Some(movimiento) = mejor_intercambio(&arbol, &conjuntos, puntaje) else {
            break;
        };
        let (a, b) = movimiento.arista;
        let (x, y) = movimiento.intercambio;
        intercambiar(&mut arbol, a, x, b, y);
        puntaje = movimiento.puntaje;
        movimientos.push(movimiento);
    }
    let parsimonia = parsimonia_sin_raiz(&arbol, secuencias)?;
    Ok(BusquedaDeParsimonia {
        limite_alcanzado: movimientos.len() == maximo_de_movimientos,
        arbol,
        parsimonia,
        puntaje_inicial,
        movimientos,
    })
}

#[derive(Debug, Clone, Default)]
struct Subarbol {
    letras: Vec<u8>,
    cambios: usize,
}

fn juntar(primero: &Subarbol, segundo: &Subarbol) -> Subarbol {
    let mut letras = vec![0; primero.letras.len()];
    let nuevos = unir_conjuntos(&primero.letras, &segundo.letras, &mut letras);
    Subarbol {
        letras,
        cambios: primero.cambios + segundo.cambios + nuevos,
    }
}

fn primera_hoja(arbol: &Arbol) -> Option<usize> {
    (0..arbol.nombres.len()).find(|&nodo| arbol.grado(nodo) == 1)
}

fn subarboles_hacia_abajo(
    arbol: &Arbol,
    orden: &[usize],
    padres: &[Option<usize>],
    conjuntos: &[Vec<u8>],
) -> Vec<Subarbol> {
    let mut abajo = vec![Subarbol::default(); arbol.nombres.len()];
    for &nodo in orden.iter().rev() {
        abajo[nodo] = match padres[nodo] {
            Some(padre) if arbol.grado(nodo) == 3 => {
                let (primero, segundo) = otros_dos_vecinos(arbol, nodo, padre);
                juntar(&abajo[primero], &abajo[segundo])
            }
            _ => Subarbol {
                letras: conjuntos[nodo].clone(),
                cambios: 0,
            },
        };
    }
    abajo
}

fn subarboles_hacia_arriba(
    arbol: &Arbol,
    orden: &[usize],
    padres: &[Option<usize>],
    abajo: &[Subarbol],
) -> Vec<Subarbol> {
    let mut arriba = vec![Subarbol::default(); arbol.nombres.len()];
    for &nodo in orden {
        let Some(padre) = padres[nodo] else {
            continue;
        };
        arriba[nodo] = match padres[padre] {
            None => abajo[padre].clone(),
            Some(abuelo) => {
                let (primero, segundo) = otros_dos_vecinos(arbol, padre, abuelo);
                let hermano = if primero == nodo { segundo } else { primero };
                juntar(&arriba[padre], &abajo[hermano])
            }
        };
    }
    arriba
}

pub fn puntaje_de_fitch(arbol: &Arbol, conjuntos: &[Vec<u8>]) -> usize {
    let Some(raiz) = primera_hoja(arbol) else {
        return 0;
    };
    let (orden, padres) = recorrer_desde(arbol, raiz);
    let abajo = subarboles_hacia_abajo(arbol, &orden, &padres, conjuntos);
    let hijo = arbol.vecinos[raiz][0].0;
    abajo[hijo].cambios + contar_sin_letra_comun(&conjuntos[raiz], &abajo[hijo].letras)
}

pub fn mejor_intercambio(
    arbol: &Arbol,
    conjuntos: &[Vec<u8>],
    puntaje_actual: usize,
) -> Option<Movimiento> {
    let raiz = primera_hoja(arbol)?;
    let (orden, padres) = recorrer_desde(arbol, raiz);
    let abajo = subarboles_hacia_abajo(arbol, &orden, &padres, conjuntos);
    let arriba = subarboles_hacia_arriba(arbol, &orden, &padres, &abajo);
    let mut izquierda = vec![0; conjuntos[raiz].len()];
    let mut derecha = izquierda.clone();
    let mut mejor: Option<Movimiento> = None;
    for &b in &orden {
        let Some(a) = padres[b] else {
            continue;
        };
        let Some(padre_de_a) = padres[a] else {
            continue;
        };
        if arbol.grado(a) != 3 || arbol.grado(b) != 3 {
            continue;
        }
        let (primero, segundo) = otros_dos_vecinos(arbol, a, padre_de_a);
        let x = if primero == b { segundo } else { primero };
        let (y, z) = otros_dos_vecinos(arbol, b, a);
        for (sube, queda) in [(y, z), (z, y)] {
            let puntaje = puntaje_tras_el_intercambio(
                (&arriba[a], &abajo[sube]),
                (&abajo[x], &abajo[queda]),
                &mut izquierda,
                &mut derecha,
            );
            let mejor_hasta_ahora = mejor.as_ref().map_or(puntaje_actual, |otro| otro.puntaje);
            if puntaje < mejor_hasta_ahora {
                mejor = Some(Movimiento {
                    arista: (a, b),
                    intercambio: (x, sube),
                    puntaje,
                });
            }
        }
    }
    mejor
}

fn puntaje_tras_el_intercambio(
    lado_de_a: (&Subarbol, &Subarbol),
    lado_de_b: (&Subarbol, &Subarbol),
    izquierda: &mut [u8],
    derecha: &mut [u8],
) -> usize {
    let (primero_de_a, segundo_de_a) = lado_de_a;
    let (primero_de_b, segundo_de_b) = lado_de_b;
    let mut puntaje =
        primero_de_a.cambios + segundo_de_a.cambios + primero_de_b.cambios + segundo_de_b.cambios;
    puntaje += unir_conjuntos(&primero_de_a.letras, &segundo_de_a.letras, izquierda);
    puntaje += unir_conjuntos(&primero_de_b.letras, &segundo_de_b.letras, derecha);
    puntaje + contar_sin_letra_comun(izquierda, derecha)
}

pub fn conjuntos_de_letras(secuencias: &[Vec<u8>], longitud: usize) -> Vec<Vec<u8>> {
    let mut conjuntos = Vec::with_capacity(secuencias.len());
    for secuencia in secuencias {
        if secuencia.is_empty() {
            conjuntos.push(vec![0; longitud]);
        } else {
            let letras = secuencia
                .iter()
                .map(|&letra| 1u8 << indice_base(letra).unwrap_or(0))
                .collect();
            conjuntos.push(letras);
        }
    }
    conjuntos
}

fn unir_conjuntos(primero: &[u8], segundo: &[u8], resultado: &mut [u8]) -> usize {
    let mut cambios: u32 = 0;
    for ((&a, &b), salida) in primero.iter().zip(segundo).zip(resultado.iter_mut()) {
        let comunes = a & b;
        let sin_letra_comun = comunes == 0;
        *salida = if sin_letra_comun { a | b } else { comunes };
        cambios += sin_letra_comun as u32;
    }
    cambios as usize
}

fn contar_sin_letra_comun(primero: &[u8], segundo: &[u8]) -> usize {
    primero
        .iter()
        .zip(segundo)
        .filter(|&(&a, &b)| a & b == 0)
        .count()
}

pub fn revisar_arbol_binario(arbol: &Arbol, secuencias: &[Vec<u8>]) -> Result<usize, String> {
    let hojas = arbol.hojas();
    if hojas.len() < 3 {
        return Err("Hacen falta al menos tres secuencias.".to_string());
    }
    let longitud = revisar_hojas_con_secuencia(arbol, hojas[0], secuencias)?;
    for (nodo, secuencia) in secuencias.iter().enumerate() {
        let grado = arbol.grado(nodo);
        if grado != 1 && grado != 3 {
            return Err(format!(
                "El nodo «{}» tiene {grado} vecinos: en un árbol binario sin raíz cada nodo interno tiene tres.",
                arbol.nombres[nodo]
            ));
        }
        if grado == 1 && secuencia.is_empty() {
            return Err(format!(
                "La hoja «{}» no tiene secuencia.",
                arbol.nombres[nodo]
            ));
        }
        if grado == 3 && !secuencia.is_empty() {
            return Err(format!(
                "El nodo interno «{}» no debe traer secuencia: se deduce.",
                arbol.nombres[nodo]
            ));
        }
    }
    Ok(longitud)
}

pub fn arbol_al_azar(nombres: &[String], semilla: u64) -> Result<Arbol, String> {
    let cantidad = nombres.len();
    if cantidad < 3 {
        return Err("Hacen falta al menos tres hojas.".to_string());
    }
    let mut azar = Azar::nuevo(semilla);
    let mut arbol = Arbol::default();
    for nombre in nombres {
        arbol.agregar_nodo(nombre);
    }
    let centro = arbol.agregar_nodo(&nombre_interno(cantidad));
    for hoja in 0..3 {
        arbol.conectar(centro, hoja, 1.0);
    }
    for hoja in 3..cantidad {
        let aristas = arbol.aristas();
        let (a, b, _) = aristas[azar.entero(aristas.len())];
        let nuevo = arbol.agregar_nodo(&nombre_interno(arbol.nombres.len()));
        arbol.desconectar(a, b);
        arbol.conectar(a, nuevo, 1.0);
        arbol.conectar(nuevo, b, 1.0);
        arbol.conectar(nuevo, hoja, 1.0);
    }
    Ok(arbol)
}
