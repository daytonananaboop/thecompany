use crate::azar::Azar;
use std::collections::HashMap;

pub const MAXIMO_K_UNIVERSAL: usize = 20;
pub const MAXIMO_ARISTAS_DE_SOLAPAMIENTO: usize = 5_000_000;
pub const INTENTOS_CON_PARES: usize = 100;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Grafo {
    pub nodos: Vec<Vec<u8>>,
    pub salidas: Vec<Vec<usize>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParDeLecturas {
    pub primera: Vec<u8>,
    pub segunda: Vec<u8>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ArmadoConPares {
    pub genoma: Vec<u8>,
    pub intentos: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MedidaN50 {
    pub n50: usize,
    pub l50: usize,
    pub total: usize,
}

pub fn composicion(texto: &[u8], k: usize) -> Vec<&[u8]> {
    if k == 0 || k > texto.len() {
        return Vec::new();
    }
    texto.windows(k).collect()
}

pub fn deletrear_camino(kmeros: &[&[u8]]) -> Result<Vec<u8>, String> {
    let k = largo_comun(kmeros)?;
    let mut texto = Vec::with_capacity(k + kmeros.len() - 1);
    texto.extend_from_slice(kmeros[0]);
    for numero in 1..kmeros.len() {
        let anterior = kmeros[numero - 1];
        let kmero = kmeros[numero];
        if anterior[1..] != kmero[..k - 1] {
            return Err(format!(
                "El k-mero {} («{}») no continúa al anterior («{}»): las últimas {} letras de uno deben ser las primeras del otro.",
                numero + 1,
                mostrar(kmero),
                mostrar(anterior),
                k - 1
            ));
        }
        texto.push(kmero[k - 1]);
    }
    Ok(texto)
}

pub fn largo_comun(kmeros: &[&[u8]]) -> Result<usize, String> {
    let Some(primero) = kmeros.first() else {
        return Err("No hay ningún k-mero.".to_string());
    };
    let k = primero.len();
    if k == 0 {
        return Err("Los k-meros no pueden estar vacíos.".to_string());
    }
    for (numero, kmero) in kmeros.iter().enumerate() {
        if kmero.len() != k {
            return Err(format!(
                "Todos los k-meros deben medir lo mismo: el 1 mide {k} y el {} mide {}.",
                numero + 1,
                kmero.len()
            ));
        }
    }
    Ok(k)
}

fn mostrar(letras: &[u8]) -> String {
    String::from_utf8_lossy(letras).into_owned()
}

pub fn grafo_de_solapamiento(kmeros: &[&[u8]]) -> Result<Grafo, String> {
    let k = largo_comun(kmeros)?;
    if k < 2 {
        return Err("Los k-meros deben tener al menos 2 letras para solaparse.".to_string());
    }
    let mut por_prefijo: HashMap<&[u8], Vec<usize>> = HashMap::with_capacity(kmeros.len());
    for (indice, kmero) in kmeros.iter().enumerate() {
        por_prefijo.entry(&kmero[..k - 1]).or_default().push(indice);
    }
    let mut aristas = 0;
    for kmero in kmeros {
        if let Some(siguientes) = por_prefijo.get(&kmero[1..]) {
            aristas += siguientes.len();
        }
    }
    if aristas > MAXIMO_ARISTAS_DE_SOLAPAMIENTO {
        return Err(format!(
            "El grafo de solapamiento tendría más de {MAXIMO_ARISTAS_DE_SOLAPAMIENTO} aristas: hay demasiados k-meros que se solapan entre sí."
        ));
    }
    let mut grafo = Grafo {
        nodos: kmeros.iter().map(|kmero| kmero.to_vec()).collect(),
        salidas: vec![Vec::new(); kmeros.len()],
    };
    for (indice, kmero) in kmeros.iter().enumerate() {
        if let Some(siguientes) = por_prefijo.get(&kmero[1..]) {
            for &siguiente in siguientes {
                if siguiente != indice {
                    grafo.salidas[indice].push(siguiente);
                }
            }
        }
    }
    Ok(grafo)
}

pub fn de_bruijn_de_texto(texto: &[u8], k: usize) -> Result<Grafo, String> {
    if k < 2 || k > texto.len() {
        return Err(format!(
            "k debe estar entre 2 y el largo del texto ({}).",
            texto.len()
        ));
    }
    de_bruijn_de_kmeros(&composicion(texto, k))
}

pub fn de_bruijn_de_kmeros(kmeros: &[&[u8]]) -> Result<Grafo, String> {
    let k = largo_comun(kmeros)?;
    if k < 2 {
        return Err("Los k-meros deben tener al menos 2 letras.".to_string());
    }
    let mut grafo = Grafo::default();
    let mut indices = HashMap::with_capacity(kmeros.len());
    for kmero in kmeros {
        let desde = nodo_para(&kmero[..k - 1], &mut indices, &mut grafo);
        let hasta = nodo_para(&kmero[1..], &mut indices, &mut grafo);
        grafo.salidas[desde].push(hasta);
    }
    Ok(grafo)
}

fn nodo_para<'a>(
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

pub fn cantidad_de_aristas(grafo: &Grafo) -> usize {
    grafo.salidas.iter().map(|destinos| destinos.len()).sum()
}

pub fn grados_de_entrada(grafo: &Grafo) -> Vec<usize> {
    let mut entradas = vec![0; grafo.salidas.len()];
    for destinos in &grafo.salidas {
        for &destino in destinos {
            entradas[destino] += 1;
        }
    }
    entradas
}

pub fn balance_de_nodos(grafo: &Grafo) -> Vec<i64> {
    let entradas = grados_de_entrada(grafo);
    let mut balance = Vec::with_capacity(entradas.len());
    for (nodo, destinos) in grafo.salidas.iter().enumerate() {
        balance.push(destinos.len() as i64 - entradas[nodo] as i64);
    }
    balance
}

pub fn ciclo_euleriano(grafo: &Grafo) -> Result<Vec<usize>, String> {
    let aristas = cantidad_de_aristas(grafo);
    if aristas == 0 {
        return Err("El grafo no tiene aristas.".to_string());
    }
    let balance = balance_de_nodos(grafo);
    if balance.iter().any(|&diferencia| diferencia != 0) {
        return Err(format!(
            "No hay ciclo euleriano: en cada nodo deben entrar tantas aristas como salen, y no pasa en {}.",
            describir_desbalance(grafo, &balance)
        ));
    }
    let inicio = primer_nodo_con_salidas(grafo);
    let ciclo = recorrido_de_hierholzer(grafo, inicio);
    if ciclo.len() != aristas + 1 {
        return Err(
            "No hay ciclo euleriano: las aristas forman grupos separados que no se tocan."
                .to_string(),
        );
    }
    Ok(ciclo)
}

pub fn camino_euleriano(grafo: &Grafo) -> Result<Vec<usize>, String> {
    let aristas = cantidad_de_aristas(grafo);
    if aristas == 0 {
        return Err("El grafo no tiene aristas.".to_string());
    }
    let balance = balance_de_nodos(grafo);
    let mut inicio = None;
    let mut fin = None;
    let mut otros_desbalanceados = 0;
    for (nodo, &diferencia) in balance.iter().enumerate() {
        if diferencia == 1 && inicio.is_none() {
            inicio = Some(nodo);
        } else if diferencia == -1 && fin.is_none() {
            fin = Some(nodo);
        } else if diferencia != 0 {
            otros_desbalanceados += 1;
        }
    }
    if otros_desbalanceados > 0 {
        return Err(format!(
            "No hay camino euleriano: solo el nodo inicial puede tener una salida de más y solo el final una entrada de más, pero hay desbalance en {}.",
            describir_desbalance(grafo, &balance)
        ));
    }
    let inicio = inicio.unwrap_or_else(|| primer_nodo_con_salidas(grafo));
    let camino = recorrido_de_hierholzer(grafo, inicio);
    if camino.len() != aristas + 1 {
        return Err(
            "No hay camino euleriano: las aristas forman grupos separados que no se tocan."
                .to_string(),
        );
    }
    Ok(camino)
}

pub fn primer_nodo_con_salidas(grafo: &Grafo) -> usize {
    grafo
        .salidas
        .iter()
        .position(|destinos| !destinos.is_empty())
        .unwrap_or(0)
}

fn describir_desbalance(grafo: &Grafo, balance: &[i64]) -> String {
    let mut partes = Vec::new();
    let mut total = 0;
    for (nodo, &diferencia) in balance.iter().enumerate() {
        if diferencia == 0 {
            continue;
        }
        total += 1;
        if partes.len() < 4 {
            let salen = grafo.salidas[nodo].len() as i64;
            partes.push(format!(
                "«{}» (entran {}, salen {salen})",
                mostrar(&grafo.nodos[nodo]),
                salen - diferencia
            ));
        }
    }
    if total > partes.len() {
        partes.push(format!("y {} nodos más", total - partes.len()));
    }
    partes.join(", ")
}

pub fn recorrido_de_hierholzer(grafo: &Grafo, inicio: usize) -> Vec<usize> {
    if inicio >= grafo.salidas.len() {
        return Vec::new();
    }
    let mut aristas_usadas = vec![0; grafo.salidas.len()];
    let mut pila = vec![inicio];
    let mut camino = Vec::with_capacity(cantidad_de_aristas(grafo) + 1);
    while let Some(&actual) = pila.last() {
        let usadas = aristas_usadas[actual];
        if usadas < grafo.salidas[actual].len() {
            let destino = grafo.salidas[actual][usadas];
            aristas_usadas[actual] += 1;
            pila.push(destino);
        } else {
            pila.pop();
            camino.push(actual);
        }
    }
    camino.reverse();
    camino
}

pub fn reconstruir_desde_kmeros(kmeros: &[&[u8]]) -> Result<Vec<u8>, String> {
    let grafo = de_bruijn_de_kmeros(kmeros)?;
    let camino = camino_euleriano(&grafo)?;
    deletrear_camino(&etiquetas_del_camino(&grafo, &camino))
}

pub fn etiquetas_del_camino<'a>(grafo: &'a Grafo, camino: &[usize]) -> Vec<&'a [u8]> {
    camino
        .iter()
        .map(|&nodo| grafo.nodos[nodo].as_slice())
        .collect()
}

pub fn cadena_universal_circular(k: usize) -> Result<Vec<u8>, String> {
    if k == 0 || k > MAXIMO_K_UNIVERSAL {
        return Err(format!("k debe estar entre 1 y {MAXIMO_K_UNIVERSAL}."));
    }
    if k == 1 {
        return Ok(b"01".to_vec());
    }
    let grafo = de_bruijn_binario(k);
    let ciclo = ciclo_euleriano(&grafo)?;
    let mut cadena = deletrear_camino(&etiquetas_del_camino(&grafo, &ciclo))?;
    cadena.truncate(cadena.len() - (k - 1));
    Ok(cadena)
}

pub fn de_bruijn_binario(k: usize) -> Grafo {
    let cantidad = 1usize << (k - 1);
    let mut grafo = Grafo {
        nodos: Vec::with_capacity(cantidad),
        salidas: Vec::with_capacity(cantidad),
    };
    for numero in 0..cantidad {
        grafo.nodos.push(en_binario(numero, k - 1));
        let sin_el_primero = (numero * 2) % cantidad;
        grafo.salidas.push(vec![sin_el_primero, sin_el_primero + 1]);
    }
    grafo
}

fn en_binario(numero: usize, largo: usize) -> Vec<u8> {
    let mut letras = Vec::with_capacity(largo);
    for posicion in (0..largo).rev() {
        if (numero >> posicion) & 1 == 1 {
            letras.push(b'1');
        } else {
            letras.push(b'0');
        }
    }
    letras
}

pub fn deletrear_camino_con_huecos(pares: &[ParDeLecturas], d: usize) -> Result<Vec<u8>, String> {
    let primeras: Vec<&[u8]> = pares.iter().map(|par| par.primera.as_slice()).collect();
    let segundas: Vec<&[u8]> = pares.iter().map(|par| par.segunda.as_slice()).collect();
    unir_prefijo_y_sufijo(&primeras, &segundas, d)
}

pub fn unir_prefijo_y_sufijo(
    primeras: &[&[u8]],
    segundas: &[&[u8]],
    d: usize,
) -> Result<Vec<u8>, String> {
    let prefijo = deletrear_camino(primeras)?;
    let sufijo = deletrear_camino(segundas)?;
    let k = primeras[0].len();
    if segundas[0].len() != k {
        return Err("Las dos lecturas de cada par deben medir lo mismo.".to_string());
    }
    let salto = k.saturating_add(d);
    if prefijo.len() < salto {
        return Err(format!(
            "Con {} pares el hueco entre las lecturas queda sin cubrir: hacen falta al menos d + 1 pares.",
            primeras.len()
        ));
    }
    for posicion in salto..prefijo.len() {
        if prefijo[posicion] != sufijo[posicion - salto] {
            return Err(format!(
                "Las lecturas no encajan: en la posición {posicion} (contando desde 0) las primeras dicen «{}» y las segundas «{}».",
                prefijo[posicion] as char,
                sufijo[posicion - salto] as char
            ));
        }
    }
    let mut genoma = prefijo;
    let ya_cubierto = genoma.len() - salto;
    genoma.extend_from_slice(&sufijo[ya_cubierto..]);
    Ok(genoma)
}

pub fn de_bruijn_de_pares(pares: &[ParDeLecturas]) -> Result<Grafo, String> {
    let primeras: Vec<&[u8]> = pares.iter().map(|par| par.primera.as_slice()).collect();
    let k = largo_comun(&primeras)?;
    if pares.iter().any(|par| par.segunda.len() != k) {
        return Err("Las dos lecturas de cada par deben medir lo mismo.".to_string());
    }
    if k < 2 {
        return Err("Las lecturas deben tener al menos 2 letras.".to_string());
    }
    let mut prefijos = Vec::with_capacity(pares.len());
    let mut sufijos = Vec::with_capacity(pares.len());
    for par in pares {
        prefijos.push(unir_par(&par.primera[..k - 1], &par.segunda[..k - 1]));
        sufijos.push(unir_par(&par.primera[1..], &par.segunda[1..]));
    }
    let mut grafo = Grafo::default();
    let mut indices = HashMap::with_capacity(pares.len());
    for (prefijo, sufijo) in prefijos.iter().zip(&sufijos) {
        let desde = nodo_para(prefijo, &mut indices, &mut grafo);
        let hasta = nodo_para(sufijo, &mut indices, &mut grafo);
        grafo.salidas[desde].push(hasta);
    }
    Ok(grafo)
}

pub fn unir_par(primera: &[u8], segunda: &[u8]) -> Vec<u8> {
    let mut etiqueta = Vec::with_capacity(primera.len() + 1 + segunda.len());
    etiqueta.extend_from_slice(primera);
    etiqueta.push(b'|');
    etiqueta.extend_from_slice(segunda);
    etiqueta
}

pub fn separar_par(etiqueta: &[u8]) -> (&[u8], &[u8]) {
    let mitad = etiqueta.len() / 2;
    (&etiqueta[..mitad], &etiqueta[mitad + 1..])
}

pub fn reconstruir_desde_pares(
    pares: &[ParDeLecturas],
    d: usize,
    semilla: u64,
) -> Result<ArmadoConPares, String> {
    let mut grafo = de_bruijn_de_pares(pares)?;
    if pares.len() < d.saturating_add(1) {
        return Err(format!(
            "Con {} pares y d = {d} el hueco entre las lecturas queda sin cubrir: hacen falta al menos d + 1 pares.",
            pares.len()
        ));
    }
    camino_euleriano(&grafo)?;
    let aristas = cantidad_de_aristas(&grafo);
    let balance = balance_de_nodos(&grafo);
    let inicio_fijo = balance.iter().position(|&diferencia| diferencia == 1);
    let candidatos: Vec<usize> = (0..grafo.salidas.len())
        .filter(|&nodo| !grafo.salidas[nodo].is_empty())
        .collect();
    let mut azar = Azar::nuevo(semilla);
    let desplazamiento = azar.entero(candidatos.len());
    let mut motivo = String::new();
    let mut intentos = 0;
    while intentos < INTENTOS_CON_PARES {
        let inicio = match inicio_fijo {
            Some(nodo) => nodo,
            None => candidatos[(intentos + desplazamiento) % candidatos.len()],
        };
        intentos += 1;
        let camino = recorrido_de_hierholzer(&grafo, inicio);
        if camino.len() == aristas + 1 {
            let mut primeras = Vec::with_capacity(camino.len());
            let mut segundas = Vec::with_capacity(camino.len());
            for &nodo in &camino {
                let (primera, segunda) = separar_par(&grafo.nodos[nodo]);
                primeras.push(primera);
                segundas.push(segunda);
            }
            match unir_prefijo_y_sufijo(&primeras, &segundas, d + 1) {
                Ok(genoma) => return Ok(ArmadoConPares { genoma, intentos }),
                Err(error) => motivo = error,
            }
        }
        if !hay_bifurcaciones(&grafo) && (inicio_fijo.is_some() || candidatos.len() == 1) {
            break;
        }
        barajar_salidas(&mut grafo, &mut azar);
    }
    Err(format!(
        "Ningún camino euleriano de los {intentos} intentos deletrea un genoma coherente. El último falló así: {motivo}"
    ))
}

pub fn hay_bifurcaciones(grafo: &Grafo) -> bool {
    grafo.salidas.iter().any(|destinos| destinos.len() > 1)
}

pub fn barajar_salidas(grafo: &mut Grafo, azar: &mut Azar) {
    for destinos in &mut grafo.salidas {
        for posicion in (1..destinos.len()).rev() {
            let otra = azar.entero(posicion + 1);
            destinos.swap(posicion, otra);
        }
    }
}

pub fn caminos_no_ramificados(grafo: &Grafo) -> Vec<Vec<usize>> {
    let entradas = grados_de_entrada(grafo);
    let cantidad = grafo.salidas.len();
    let mut uno_a_uno = vec![false; cantidad];
    for nodo in 0..cantidad {
        uno_a_uno[nodo] = entradas[nodo] == 1 && grafo.salidas[nodo].len() == 1;
    }
    let mut visitado = vec![false; cantidad];
    let mut caminos = Vec::new();
    for nodo in 0..cantidad {
        if uno_a_uno[nodo] {
            continue;
        }
        for &siguiente in &grafo.salidas[nodo] {
            let mut camino = vec![nodo, siguiente];
            let mut actual = siguiente;
            while uno_a_uno[actual] {
                visitado[actual] = true;
                actual = grafo.salidas[actual][0];
                camino.push(actual);
            }
            caminos.push(camino);
        }
    }
    for nodo in 0..cantidad {
        if !uno_a_uno[nodo] || visitado[nodo] {
            continue;
        }
        let mut ciclo = vec![nodo];
        let mut actual = nodo;
        loop {
            visitado[actual] = true;
            actual = grafo.salidas[actual][0];
            ciclo.push(actual);
            if actual == nodo {
                break;
            }
        }
        caminos.push(ciclo);
    }
    caminos
}

pub fn generar_contigs(kmeros: &[&[u8]]) -> Result<Vec<Vec<u8>>, String> {
    let grafo = de_bruijn_de_kmeros(kmeros)?;
    let mut contigs = Vec::new();
    for camino in caminos_no_ramificados(&grafo) {
        contigs.push(deletrear_camino(&etiquetas_del_camino(&grafo, &camino))?);
    }
    Ok(contigs)
}

pub fn medir_n50(longitudes: &[usize]) -> MedidaN50 {
    let mut ordenadas = longitudes.to_vec();
    ordenadas.sort_unstable_by(|una, otra| otra.cmp(una));
    let mut total: usize = 0;
    for &longitud in &ordenadas {
        total = total.saturating_add(longitud);
    }
    let mut medida = MedidaN50 {
        n50: 0,
        l50: 0,
        total,
    };
    let mut acumulado: usize = 0;
    for (posicion, &longitud) in ordenadas.iter().enumerate() {
        acumulado = acumulado.saturating_add(longitud);
        if total > 0 && acumulado >= total - acumulado {
            medida.n50 = longitud;
            medida.l50 = posicion + 1;
            break;
        }
    }
    medida
}

pub fn leer_lista_de_adyacencia(texto: &str) -> Result<Grafo, String> {
    let mut grafo = Grafo::default();
    let mut indices = HashMap::new();
    for (numero, linea) in texto.lines().enumerate() {
        let linea = linea.trim();
        if linea.is_empty() {
            continue;
        }
        let Some((origen, destinos)) = linea.split_once("->").or_else(|| linea.split_once(':'))
        else {
            return Err(format!(
                "La línea {} no tiene «->»: escribe «origen -> destino1,destino2».",
                numero + 1
            ));
        };
        let origen = origen.trim();
        if origen.is_empty() || origen.contains(|letra: char| letra.is_whitespace() || letra == ',')
        {
            return Err(format!(
                "En la línea {} falta el nodo de origen, o tiene espacios o comas.",
                numero + 1
            ));
        }
        if destinos.contains("->") {
            return Err(format!(
                "La línea {} tiene más de una flecha «->»: escribe un solo origen por línea.",
                numero + 1
            ));
        }
        let desde = nodo_para(origen.as_bytes(), &mut indices, &mut grafo);
        for destino in destinos.split(|letra: char| letra == ',' || letra.is_whitespace()) {
            if destino.is_empty() {
                continue;
            }
            let hasta = nodo_para(destino.as_bytes(), &mut indices, &mut grafo);
            grafo.salidas[desde].push(hasta);
        }
    }
    if grafo.nodos.is_empty() {
        return Err("No hay ninguna línea de la forma «origen -> destino».".to_string());
    }
    Ok(grafo)
}

pub fn leer_pares(texto: &str) -> Result<Vec<ParDeLecturas>, String> {
    let mut limpio = texto.replace('\t', " ");
    while limpio.contains(" |") || limpio.contains("| ") {
        limpio = limpio.replace(" |", "|").replace("| ", "|");
    }
    let mut pares = Vec::new();
    for pieza in limpio.split_whitespace() {
        let Some((primera, segunda)) = pieza.split_once('|') else {
            return Err(format!(
                "«{pieza}» no es un par: escribe las dos lecturas separadas por «|»."
            ));
        };
        let es_lectura = |lectura: &str| {
            !lectura.is_empty() && lectura.bytes().all(|letra| letra.is_ascii_alphabetic())
        };
        if !es_lectura(primera) || !es_lectura(segunda) {
            return Err(format!(
                "«{pieza}» no es un par válido: cada lectura debe tener solo letras."
            ));
        }
        pares.push(ParDeLecturas {
            primera: primera.to_ascii_uppercase().into_bytes(),
            segunda: segunda.to_ascii_uppercase().into_bytes(),
        });
    }
    if pares.is_empty() {
        return Err("No hay ningún par de lecturas.".to_string());
    }
    Ok(pares)
}
