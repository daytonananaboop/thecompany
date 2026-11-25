use crate::secuencia::complemento_inverso;
use std::collections::HashMap;

pub const BASES: &[u8] = b"ACGT";
const NO_ES_BASE: u8 = 4;
pub const MAXIMO_K_ENUMERABLE: usize = 12;
pub const MAXIMO_K_NUMERICO: usize = 31;
pub const MAXIMO_D: usize = 3;
pub const MAXIMO_VECINOS: usize = 1_000_000;
pub const MAXIMO_VISITAS: usize = 200_000_000;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Grumo {
    pub patron: Vec<u8>,
    pub ventana: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OrigenProbable {
    pub minimos: Vec<usize>,
    pub inicio: usize,
    pub fin: usize,
    pub candidatos: Vec<Vec<u8>>,
    pub veces: u32,
}

pub fn contar_patron(texto: &[u8], patron: &[u8]) -> usize {
    if patron.is_empty() {
        return 0;
    }
    let mut veces = 0;
    for ventana in texto.windows(patron.len()) {
        if ventana == patron {
            veces += 1;
        }
    }
    veces
}

pub fn tabla_de_frecuencias(texto: &[u8], k: usize) -> HashMap<&[u8], usize> {
    let mut tabla = HashMap::new();
    if k == 0 {
        return tabla;
    }
    for kmero in texto.windows(k) {
        if kmero.iter().all(|&base| valor_de_base(base).is_some()) {
            *tabla.entry(kmero).or_insert(0) += 1;
        }
    }
    tabla
}

pub fn conviene_un_arreglo(longitud: usize, k: usize) -> bool {
    (1..=MAXIMO_K_ENUMERABLE).contains(&k) && (1 << (2 * k)) <= longitud.saturating_mul(4)
}

pub fn palabras_frecuentes(texto: &[u8], k: usize) -> Vec<Vec<u8>> {
    if conviene_un_arreglo(texto.len(), k)
        && let Ok(frecuencias) = arreglo_de_frecuencias(texto, k)
    {
        return los_mas_frecuentes(&frecuencias, k);
    }
    let tabla = tabla_de_frecuencias(texto, k);
    let maximo = tabla.values().copied().max().unwrap_or(0);
    let mut frecuentes = Vec::new();
    for (kmero, veces) in tabla {
        if veces == maximo {
            frecuentes.push(kmero.to_vec());
        }
    }
    frecuentes.sort();
    frecuentes
}

pub fn posiciones_del_patron(texto: &[u8], patron: &[u8]) -> Vec<usize> {
    let mut posiciones = Vec::new();
    if patron.is_empty() {
        return posiciones;
    }
    for (posicion, ventana) in texto.windows(patron.len()).enumerate() {
        if ventana == patron {
            posiciones.push(posicion);
        }
    }
    posiciones
}

pub fn encontrar_grumos(
    texto: &[u8],
    k: usize,
    largo: usize,
    veces: usize,
) -> Result<Vec<Grumo>, String> {
    if veces == 0 || largo < k || texto.len() < largo {
        return Ok(Vec::new());
    }
    let mut en_la_ventana = arreglo_de_frecuencias(&texto[..largo], k)?;
    let mut primera_ventana: HashMap<u64, usize> = HashMap::new();
    for kmero in texto[..largo].windows(k) {
        if let Ok(numero) = numero_de_patron(kmero)
            && en_la_ventana[numero as usize] as usize >= veces
        {
            primera_ventana.entry(numero).or_insert(0);
        }
    }
    for inicio in 1..=texto.len() - largo {
        let saliente = &texto[inicio - 1..inicio - 1 + k];
        if let Ok(numero) = numero_de_patron(saliente) {
            en_la_ventana[numero as usize] -= 1;
        }
        let ultimo = inicio + largo - k;
        let entrante = &texto[ultimo..ultimo + k];
        if let Ok(numero) = numero_de_patron(entrante) {
            en_la_ventana[numero as usize] += 1;
            if en_la_ventana[numero as usize] as usize >= veces {
                primera_ventana.entry(numero).or_insert(inicio);
            }
        }
    }
    let mut grumos = Vec::new();
    for (numero, ventana) in primera_ventana {
        grumos.push(Grumo {
            patron: patron_de_numero(numero, k)?,
            ventana,
        });
    }
    grumos.sort_by(|uno, otro| uno.patron.cmp(&otro.patron));
    Ok(grumos)
}

pub fn cambio_de_sesgo(base: u8) -> i64 {
    match base {
        b'G' | b'g' => 1,
        b'C' | b'c' => -1,
        _ => 0,
    }
}

pub fn sesgo(genoma: &[u8]) -> Vec<i64> {
    let mut valores = Vec::with_capacity(genoma.len() + 1);
    let mut actual = 0;
    valores.push(actual);
    for &base in genoma {
        actual += cambio_de_sesgo(base);
        valores.push(actual);
    }
    valores
}

pub fn posiciones_de_sesgo_minimo(genoma: &[u8]) -> Vec<usize> {
    let mut sesgo_actual = 0;
    let mut minimo = 0;
    let mut posiciones = vec![0];
    for (indice, &base) in genoma.iter().enumerate() {
        sesgo_actual += cambio_de_sesgo(base);
        let posicion = indice + 1;
        if sesgo_actual < minimo {
            minimo = sesgo_actual;
            posiciones.clear();
            posiciones.push(posicion);
        } else if sesgo_actual == minimo {
            posiciones.push(posicion);
        }
    }
    posiciones
}

pub fn distancia_hamming(primera: &[u8], segunda: &[u8]) -> Result<usize, String> {
    if primera.len() != segunda.len() {
        return Err(format!(
            "Las dos secuencias deben medir lo mismo, pero una tiene {} letras y la otra {}.",
            primera.len(),
            segunda.len()
        ));
    }
    let mut diferencias = 0;
    for (a, b) in primera.iter().zip(segunda) {
        if a != b {
            diferencias += 1;
        }
    }
    Ok(diferencias)
}

pub fn casi_iguales(ventana: &[u8], patron: &[u8], d: usize) -> bool {
    let mut diferencias = 0;
    for (a, b) in ventana.iter().zip(patron) {
        if a != b {
            diferencias += 1;
            if diferencias > d {
                return false;
            }
        }
    }
    true
}

pub fn posiciones_aproximadas(texto: &[u8], patron: &[u8], d: usize) -> Vec<usize> {
    let mut posiciones = Vec::new();
    if patron.is_empty() {
        return posiciones;
    }
    for (posicion, ventana) in texto.windows(patron.len()).enumerate() {
        if casi_iguales(ventana, patron, d) {
            posiciones.push(posicion);
        }
    }
    posiciones
}

pub fn contar_aproximado(texto: &[u8], patron: &[u8], d: usize) -> usize {
    if patron.is_empty() {
        return 0;
    }
    let mut veces = 0;
    for ventana in texto.windows(patron.len()) {
        if casi_iguales(ventana, patron, d) {
            veces += 1;
        }
    }
    veces
}

pub fn valor_de_base(base: u8) -> Option<u64> {
    const VALORES: [u8; 256] = {
        let mut valores = [NO_ES_BASE; 256];
        valores[b'A' as usize] = 0;
        valores[b'C' as usize] = 1;
        valores[b'G' as usize] = 2;
        valores[b'T' as usize] = 3;
        valores[b'a' as usize] = 0;
        valores[b'c' as usize] = 1;
        valores[b'g' as usize] = 2;
        valores[b't' as usize] = 3;
        valores
    };
    let valor = VALORES[base as usize];
    if valor == NO_ES_BASE {
        None
    } else {
        Some(valor as u64)
    }
}

pub fn numero_de_patron(patron: &[u8]) -> Result<u64, String> {
    if patron.len() > MAXIMO_K_NUMERICO {
        return Err(format!(
            "Un número de 64 bits guarda patrones de hasta {MAXIMO_K_NUMERICO} letras, y este tiene {}.",
            patron.len()
        ));
    }
    let mut numero = 0;
    for &base in patron {
        let Some(valor) = valor_de_base(base) else {
            return Err(format!("«{}» no es A, C, G ni T.", base as char));
        };
        numero = numero * 4 + valor;
    }
    Ok(numero)
}

pub fn patron_de_numero(numero: u64, k: usize) -> Result<Vec<u8>, String> {
    if k > MAXIMO_K_NUMERICO {
        return Err(format!(
            "Un número de 64 bits guarda patrones de hasta {MAXIMO_K_NUMERICO} letras, y pediste {k}."
        ));
    }
    let cantidad: u64 = 1 << (2 * k);
    if numero >= cantidad {
        return Err(format!(
            "Con {k} letras hay {cantidad} patrones, numerados del 0 al {}; {numero} se pasa.",
            cantidad - 1
        ));
    }
    let mut patron = vec![b'A'; k];
    let mut resto = numero;
    for posicion in (0..k).rev() {
        patron[posicion] = BASES[(resto % 4) as usize];
        resto /= 4;
    }
    Ok(patron)
}

pub fn arreglo_de_frecuencias(texto: &[u8], k: usize) -> Result<Vec<u32>, String> {
    if k == 0 || k > MAXIMO_K_ENUMERABLE {
        return Err(format!(
            "k debe estar entre 1 y {MAXIMO_K_ENUMERABLE}: el arreglo tiene 4^k casillas y con k = {k} no cabe."
        ));
    }
    let mut frecuencias = vec![0u32; 1 << (2 * k)];
    for kmero in texto.windows(k) {
        if let Ok(numero) = numero_de_patron(kmero) {
            frecuencias[numero as usize] += 1;
        }
    }
    Ok(frecuencias)
}

pub fn tamano_del_vecindario(k: usize, d: usize) -> f64 {
    let mut total = 0.0;
    let mut combinaciones = 1.0;
    for cambios in 0..=d.min(k) {
        if cambios > 0 {
            combinaciones = combinaciones * (k - cambios + 1) as f64 / cambios as f64;
        }
        total += combinaciones * 3f64.powi(cambios as i32);
    }
    total
}

fn recorrer_vecinos(
    numero: u64,
    k: usize,
    desde: usize,
    cambios_restantes: usize,
    visitar: &mut impl FnMut(u64),
) {
    visitar(numero);
    if cambios_restantes == 0 {
        return;
    }
    for posicion in desde..k {
        let peso = 1u64 << (2 * (k - 1 - posicion));
        let original = numero / peso % 4;
        for base in 0..4 {
            if base != original {
                let vecino = numero - original * peso + base * peso;
                recorrer_vecinos(vecino, k, posicion + 1, cambios_restantes - 1, visitar);
            }
        }
    }
}

pub fn vecindario(patron: &[u8], d: usize) -> Result<Vec<Vec<u8>>, String> {
    let numero = numero_de_patron(patron)?;
    let k = patron.len();
    if d > k {
        return Err(format!("d = {d} es mayor que el largo del patrón ({k})."));
    }
    let tamano = tamano_del_vecindario(k, d);
    if tamano > MAXIMO_VECINOS as f64 {
        return Err(format!(
            "El vecindario tendría {tamano:.0} patrones; el límite es {MAXIMO_VECINOS}."
        ));
    }
    let mut numeros = Vec::with_capacity(tamano as usize);
    recorrer_vecinos(numero, k, 0, d, &mut |vecino| numeros.push(vecino));
    numeros.sort();
    let mut vecinos = Vec::with_capacity(numeros.len());
    for vecino in numeros {
        vecinos.push(patron_de_numero(vecino, k)?);
    }
    Ok(vecinos)
}

fn repartir_entre_vecinos(exactas: &[u32], k: usize, d: usize) -> Result<Vec<u32>, String> {
    if d > k || d > MAXIMO_D {
        return Err(format!(
            "d debe ser como mucho {} (no más que k ni que {MAXIMO_D}), y es {d}.",
            k.min(MAXIMO_D)
        ));
    }
    let distintos = exactas.iter().filter(|&&veces| veces > 0).count();
    let visitas = distintos as f64 * tamano_del_vecindario(k, d);
    if visitas > MAXIMO_VISITAS as f64 {
        return Err(format!(
            "Habría que visitar {visitas:.0} vecinos (límite {MAXIMO_VISITAS}): usa un texto más corto, un k menor o un d menor."
        ));
    }
    let mut aproximadas = vec![0u32; exactas.len()];
    for (numero, &veces) in exactas.iter().enumerate() {
        if veces > 0 {
            recorrer_vecinos(numero as u64, k, 0, d, &mut |vecino| {
                aproximadas[vecino as usize] += veces;
            });
        }
    }
    Ok(aproximadas)
}

pub fn cuentas_con_errores(texto: &[u8], k: usize, d: usize) -> Result<Vec<u32>, String> {
    let exactas = arreglo_de_frecuencias(texto, k)?;
    repartir_entre_vecinos(&exactas, k, d)
}

pub fn cuentas_en_ambas_hebras(texto: &[u8], k: usize, d: usize) -> Result<Vec<u32>, String> {
    let mut exactas = arreglo_de_frecuencias(texto, k)?;
    let otra_hebra = arreglo_de_frecuencias(&complemento_inverso(texto), k)?;
    for (cuenta, otra) in exactas.iter_mut().zip(otra_hebra) {
        *cuenta += otra;
    }
    repartir_entre_vecinos(&exactas, k, d)
}

pub fn los_mas_frecuentes(cuentas: &[u32], k: usize) -> Vec<Vec<u8>> {
    let maximo = cuentas.iter().copied().max().unwrap_or(0);
    let mut patrones = Vec::new();
    if maximo == 0 {
        return patrones;
    }
    for (numero, &veces) in cuentas.iter().enumerate() {
        if veces == maximo
            && let Ok(patron) = patron_de_numero(numero as u64, k)
        {
            patrones.push(patron);
        }
    }
    patrones
}

pub fn palabras_frecuentes_con_errores(
    texto: &[u8],
    k: usize,
    d: usize,
) -> Result<Vec<Vec<u8>>, String> {
    let cuentas = cuentas_con_errores(texto, k, d)?;
    Ok(los_mas_frecuentes(&cuentas, k))
}

pub fn palabras_frecuentes_con_errores_y_complementos(
    texto: &[u8],
    k: usize,
    d: usize,
) -> Result<Vec<Vec<u8>>, String> {
    let cuentas = cuentas_en_ambas_hebras(texto, k, d)?;
    Ok(los_mas_frecuentes(&cuentas, k))
}

pub fn buscar_ori(
    genoma: &[u8],
    ancho: usize,
    k: usize,
    d: usize,
) -> Result<OrigenProbable, String> {
    if ancho < k || ancho > genoma.len() {
        return Err(format!(
            "La ventana debe medir entre k ({k}) y el largo del genoma ({}), y mide {ancho}.",
            genoma.len()
        ));
    }
    let minimos = posiciones_de_sesgo_minimo(genoma);
    let inicio = minimos[0].min(genoma.len() - ancho);
    let fin = inicio + ancho;
    let cuentas = cuentas_en_ambas_hebras(&genoma[inicio..fin], k, d)?;
    let veces = cuentas.iter().copied().max().unwrap_or(0);
    let candidatos = los_mas_frecuentes(&cuentas, k);
    Ok(OrigenProbable {
        minimos,
        inicio,
        fin,
        candidatos,
        veces,
    })
}
