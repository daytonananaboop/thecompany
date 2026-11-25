use crate::azar::Azar;

pub const BASES: &[u8] = b"ACGT";

const K_MAXIMA_EXHAUSTIVA: usize = 10;
const D_MAXIMA: usize = 3;
pub const EMPATES_GUARDADOS: usize = 100;
pub const PSEUDOCONTEO_DE_LAPLACE: f64 = 1.0;
const SIN_VISITAR: u8 = u8::MAX;
const SIN_FILA: u8 = 4;

fn fila_de(base: u8) -> Option<usize> {
    const FILA_DE_CADA_LETRA: [u8; 256] = {
        let mut filas = [SIN_FILA; 256];
        filas[b'A' as usize] = 0;
        filas[b'C' as usize] = 1;
        filas[b'G' as usize] = 2;
        filas[b'T' as usize] = 3;
        filas[b'a' as usize] = 0;
        filas[b'c' as usize] = 1;
        filas[b'g' as usize] = 2;
        filas[b't' as usize] = 3;
        filas
    };
    let fila = FILA_DE_CADA_LETRA[base as usize];
    if fila == SIN_FILA {
        None
    } else {
        Some(fila as usize)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct MotivosEncontrados {
    pub posiciones: Vec<usize>,
    pub puntos: usize,
    pub puntos_por_intento: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Mediana {
    pub distancia: usize,
    pub patrones: Vec<Vec<u8>>,
    pub empatados: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AnalisisMotivos {
    pub conteos: Vec<Vec<usize>>,
    pub perfil: Vec<Vec<f64>>,
    pub consenso: Vec<u8>,
    pub puntos: usize,
    pub entropia: Vec<f64>,
}

fn validar(adn: &[Vec<u8>], k: usize) -> Result<(), String> {
    if adn.is_empty() {
        return Err("Hace falta al menos una secuencia.".to_string());
    }
    if k == 0 {
        return Err("k debe ser al menos 1.".to_string());
    }
    for (numero, texto) in adn.iter().enumerate() {
        if texto.len() < k {
            return Err(format!(
                "La secuencia {} mide {} pb, menos que k = {k}.",
                numero + 1,
                texto.len()
            ));
        }
    }
    Ok(())
}

pub fn distancia_hamming(a: &[u8], b: &[u8]) -> usize {
    let mut diferencias = 0;
    for (x, y) in a.iter().zip(b) {
        if x != y {
            diferencias += 1;
        }
    }
    diferencias
}

pub fn coincidencia_mas_cercana(patron: &[u8], texto: &[u8]) -> (usize, usize) {
    let k = patron.len();
    if k == 0 || texto.len() < k {
        return (0, k);
    }
    let mut mejor_posicion = 0;
    let mut menor_distancia = k + 1;
    for (posicion, ventana) in texto.windows(k).enumerate() {
        let distancia = distancia_hamming(patron, ventana);
        if distancia < menor_distancia {
            menor_distancia = distancia;
            mejor_posicion = posicion;
            if distancia == 0 {
                break;
            }
        }
    }
    (mejor_posicion, menor_distancia)
}

pub fn distancia_patron_texto(patron: &[u8], texto: &[u8]) -> usize {
    let (_, distancia) = coincidencia_mas_cercana(patron, texto);
    distancia
}

pub fn distancia_patron_textos(patron: &[u8], adn: &[Vec<u8>]) -> Result<usize, String> {
    validar(adn, patron.len())?;
    let mut total = 0;
    for texto in adn {
        total += distancia_patron_texto(patron, texto);
    }
    Ok(total)
}

pub fn patron_a_numero(patron: &[u8]) -> Option<usize> {
    let mut numero = 0;
    for &base in patron {
        numero = numero * 4 + fila_de(base)?;
    }
    Some(numero)
}

pub fn numero_a_patron(numero: usize, k: usize) -> Vec<u8> {
    let mut patron = vec![b'A'; k];
    let mut resto = numero;
    for base in patron.iter_mut().rev() {
        *base = BASES[resto % 4];
        resto /= 4;
    }
    patron
}

pub fn marcar_vecinos(
    numero: usize,
    k: usize,
    desde: usize,
    errores_permitidos: usize,
    marcados: &mut [bool],
) {
    marcados[numero] = true;
    if errores_permitidos == 0 {
        return;
    }
    for posicion in desde..k {
        let desplazamiento = 2 * posicion;
        let original = (numero >> desplazamiento) & 0b11;
        for base in 0..4 {
            if base != original {
                let vecino = (numero & !(0b11 << desplazamiento)) | (base << desplazamiento);
                marcar_vecinos(vecino, k, posicion + 1, errores_permitidos - 1, marcados);
            }
        }
    }
}

fn revisar_busqueda_exhaustiva(k: usize) -> Result<(), String> {
    if k > K_MAXIMA_EXHAUSTIVA {
        return Err(format!(
            "k = {k} es demasiado grande para revisar todos los patrones: hay 4^k de ellos y el máximo es k = {K_MAXIMA_EXHAUSTIVA} (4^10 = 1 048 576)."
        ));
    }
    Ok(())
}

pub fn enumerar_motivos(adn: &[Vec<u8>], k: usize, d: usize) -> Result<Vec<Vec<u8>>, String> {
    validar(adn, k)?;
    revisar_busqueda_exhaustiva(k)?;
    if d > D_MAXIMA {
        return Err(format!(
            "d = {d} es demasiado grande: el vecindario de cada k-mer crece muy rápido y el máximo es d = {D_MAXIMA}."
        ));
    }
    if d > k {
        return Err(format!("d = {d} no puede ser mayor que k = {k}."));
    }
    let cantidad_de_patrones = 1 << (2 * k);
    let mut en_todas = vec![true; cantidad_de_patrones];
    let mut en_esta = vec![false; cantidad_de_patrones];
    for texto in adn {
        en_esta.fill(false);
        for ventana in texto.windows(k) {
            if let Some(numero) = patron_a_numero(ventana) {
                marcar_vecinos(numero, k, 0, d, &mut en_esta);
            }
        }
        for (todas, esta) in en_todas.iter_mut().zip(&en_esta) {
            *todas = *todas && *esta;
        }
    }
    let mut motivos = Vec::new();
    for (numero, &presente) in en_todas.iter().enumerate() {
        if presente {
            motivos.push(numero_a_patron(numero, k));
        }
    }
    Ok(motivos)
}

pub fn distancias_a_un_texto(texto: &[u8], k: usize, distancias: &mut [u8]) {
    distancias.fill(SIN_VISITAR);
    let mut frente = Vec::new();
    for ventana in texto.windows(k) {
        if let Some(numero) = patron_a_numero(ventana)
            && distancias[numero] == SIN_VISITAR
        {
            distancias[numero] = 0;
            frente.push(numero);
        }
    }
    let mut distancia = 0;
    while !frente.is_empty() {
        distancia += 1;
        let mut siguiente_frente = Vec::new();
        for &numero in &frente {
            for posicion in 0..k {
                let desplazamiento = 2 * posicion;
                for base in 0..4 {
                    let vecino = (numero & !(0b11 << desplazamiento)) | (base << desplazamiento);
                    if distancias[vecino] == SIN_VISITAR {
                        distancias[vecino] = distancia;
                        siguiente_frente.push(vecino);
                    }
                }
            }
        }
        frente = siguiente_frente;
    }
}

pub fn cadena_mediana(adn: &[Vec<u8>], k: usize) -> Result<Mediana, String> {
    validar(adn, k)?;
    revisar_busqueda_exhaustiva(k)?;
    let cantidad_de_patrones = 1 << (2 * k);
    let mut suma_de_distancias = vec![0usize; cantidad_de_patrones];
    let mut distancias = vec![SIN_VISITAR; cantidad_de_patrones];
    for (numero, texto) in adn.iter().enumerate() {
        if !texto
            .windows(k)
            .any(|ventana| patron_a_numero(ventana).is_some())
        {
            return Err(format!(
                "La secuencia {} no tiene ningún k-mer formado solo por A, C, G y T.",
                numero + 1
            ));
        }
        distancias_a_un_texto(texto, k, &mut distancias);
        for (suma, &distancia) in suma_de_distancias.iter_mut().zip(&distancias) {
            *suma += distancia as usize;
        }
    }
    let mut mediana = Mediana {
        distancia: usize::MAX,
        patrones: Vec::new(),
        empatados: 0,
    };
    for (numero, &distancia) in suma_de_distancias.iter().enumerate() {
        if distancia < mediana.distancia {
            mediana.distancia = distancia;
            mediana.patrones.clear();
            mediana.empatados = 0;
        }
        if distancia == mediana.distancia {
            mediana.empatados += 1;
            if mediana.patrones.len() < EMPATES_GUARDADOS {
                mediana.patrones.push(numero_a_patron(numero, k));
            }
        }
    }
    Ok(mediana)
}

pub fn contar_bases(motivos: &[&[u8]]) -> Vec<Vec<usize>> {
    let k = motivos.first().map_or(0, |motivo| motivo.len());
    let mut conteos = vec![vec![0; k]; 4];
    for motivo in motivos {
        for (columna, &base) in motivo.iter().take(k).enumerate() {
            if let Some(fila) = fila_de(base) {
                conteos[fila][columna] += 1;
            }
        }
    }
    conteos
}

pub fn perfil_de_motivos(motivos: &[&[u8]], pseudoconteo: f64) -> Vec<Vec<f64>> {
    let total = motivos.len() as f64 + 4.0 * pseudoconteo;
    let mut perfil = Vec::with_capacity(4);
    for conteos_de_la_base in contar_bases(motivos) {
        let probabilidades = conteos_de_la_base
            .iter()
            .map(|&conteo| (conteo as f64 + pseudoconteo) / total)
            .collect();
        perfil.push(probabilidades);
    }
    perfil
}

pub fn fila_mayor(conteos: &[Vec<usize>], columna: usize) -> usize {
    let mut mayor = 0;
    for (fila, conteos_de_la_base) in conteos.iter().enumerate() {
        if conteos_de_la_base[columna] > conteos[mayor][columna] {
            mayor = fila;
        }
    }
    mayor
}

pub fn consenso(motivos: &[&[u8]]) -> Vec<u8> {
    let conteos = contar_bases(motivos);
    let k = conteos[0].len();
    let mut letras = Vec::with_capacity(k);
    for columna in 0..k {
        letras.push(BASES[fila_mayor(&conteos, columna)]);
    }
    letras
}

pub fn puntuacion(motivos: &[&[u8]]) -> usize {
    let conteos = contar_bases(motivos);
    let k = conteos[0].len();
    let mut diferencias = 0;
    for columna in 0..k {
        let mas_comun = fila_mayor(&conteos, columna);
        diferencias += motivos.len() - conteos[mas_comun][columna];
    }
    diferencias
}

pub fn entropia_por_columna(motivos: &[&[u8]]) -> Vec<f64> {
    let conteos = contar_bases(motivos);
    let k = conteos[0].len();
    let total = motivos.len() as f64;
    let mut entropias = Vec::with_capacity(k);
    for columna in 0..k {
        let mut entropia = 0.0;
        for conteos_de_la_base in &conteos {
            let conteo = conteos_de_la_base[columna];
            if conteo > 0 {
                let p = conteo as f64 / total;
                entropia -= p * p.log2();
            }
        }
        entropias.push(entropia);
    }
    entropias
}

pub fn analizar_motivos(motivos: &[&[u8]]) -> Result<AnalisisMotivos, String> {
    let Some(primero) = motivos.first() else {
        return Err("Hace falta al menos un motivo.".to_string());
    };
    if primero.is_empty() {
        return Err("Los motivos no pueden estar vacíos.".to_string());
    }
    if motivos.iter().any(|motivo| motivo.len() != primero.len()) {
        return Err("Todos los motivos deben medir lo mismo.".to_string());
    }
    Ok(AnalisisMotivos {
        conteos: contar_bases(motivos),
        perfil: perfil_de_motivos(motivos, 0.0),
        consenso: consenso(motivos),
        puntos: puntuacion(motivos),
        entropia: entropia_por_columna(motivos),
    })
}

fn columnas_del_perfil(perfil: &[Vec<f64>]) -> Option<usize> {
    if perfil.len() != 4 {
        return None;
    }
    let k = perfil[0].len();
    if perfil.iter().all(|fila| fila.len() == k) {
        Some(k)
    } else {
        None
    }
}

pub fn probabilidad_kmer(kmer: &[u8], perfil: &[Vec<f64>]) -> f64 {
    if perfil.len() != 4 || perfil[0].len() != kmer.len() {
        return 0.0;
    }
    let mut probabilidad = 1.0;
    for (columna, &base) in kmer.iter().enumerate() {
        probabilidad *= match fila_de(base) {
            Some(fila) => perfil[fila].get(columna).copied().unwrap_or(0.0),
            None => 0.0,
        };
    }
    probabilidad
}

pub fn posicion_mas_probable(texto: &[u8], perfil: &[Vec<f64>]) -> usize {
    let Some(k) = columnas_del_perfil(perfil) else {
        return 0;
    };
    if k == 0 || texto.len() < k {
        return 0;
    }
    let mut mejor_posicion = 0;
    let mut mayor_probabilidad = -1.0;
    for (posicion, ventana) in texto.windows(k).enumerate() {
        let probabilidad = probabilidad_kmer(ventana, perfil);
        if probabilidad > mayor_probabilidad {
            mayor_probabilidad = probabilidad;
            mejor_posicion = posicion;
        }
    }
    mejor_posicion
}

pub fn motivos_en<'a>(adn: &'a [Vec<u8>], posiciones: &[usize], k: usize) -> Vec<&'a [u8]> {
    let mut motivos = Vec::with_capacity(posiciones.len());
    for (texto, &inicio) in adn.iter().zip(posiciones) {
        let motivo = texto.get(inicio..inicio.saturating_add(k));
        motivos.push(motivo.unwrap_or_default());
    }
    motivos
}

pub fn busqueda_codiciosa(
    adn: &[Vec<u8>],
    k: usize,
    pseudoconteo: f64,
) -> Result<MotivosEncontrados, String> {
    validar(adn, k)?;
    if !pseudoconteo.is_finite() || pseudoconteo < 0.0 {
        return Err("El pseudoconteo debe ser un número mayor o igual a cero.".to_string());
    }
    let mut mejores = vec![0; adn.len()];
    let mut mejor_puntuacion = puntuacion(&motivos_en(adn, &mejores, k));
    let mut posiciones = vec![0; adn.len()];
    let mut puntos_por_intento = Vec::new();
    for inicio in 0..=adn[0].len() - k {
        posiciones[0] = inicio;
        for i in 1..adn.len() {
            let elegidos = motivos_en(&adn[..i], &posiciones[..i], k);
            let perfil = perfil_de_motivos(&elegidos, pseudoconteo);
            posiciones[i] = posicion_mas_probable(&adn[i], &perfil);
        }
        let puntos = puntuacion(&motivos_en(adn, &posiciones, k));
        puntos_por_intento.push(puntos);
        if puntos < mejor_puntuacion {
            mejor_puntuacion = puntos;
            mejores.copy_from_slice(&posiciones);
        }
    }
    Ok(MotivosEncontrados {
        posiciones: mejores,
        puntos: mejor_puntuacion,
        puntos_por_intento,
    })
}

pub fn posiciones_al_azar(adn: &[Vec<u8>], k: usize, azar: &mut Azar) -> Vec<usize> {
    let mut posiciones = Vec::with_capacity(adn.len());
    for texto in adn {
        posiciones.push(azar.entero(texto.len() - k + 1));
    }
    posiciones
}

pub fn posiciones_mas_probables(adn: &[Vec<u8>], perfil: &[Vec<f64>]) -> Vec<usize> {
    let mut posiciones = Vec::with_capacity(adn.len());
    for texto in adn {
        posiciones.push(posicion_mas_probable(texto, perfil));
    }
    posiciones
}

fn mejorar_desde_el_azar(adn: &[Vec<u8>], k: usize, azar: &mut Azar) -> (Vec<usize>, usize) {
    let mut mejores = posiciones_al_azar(adn, k, azar);
    let mut mejor_puntuacion = puntuacion(&motivos_en(adn, &mejores, k));
    loop {
        let perfil = perfil_de_motivos(&motivos_en(adn, &mejores, k), PSEUDOCONTEO_DE_LAPLACE);
        let posiciones = posiciones_mas_probables(adn, &perfil);
        let puntos = puntuacion(&motivos_en(adn, &posiciones, k));
        if puntos >= mejor_puntuacion {
            return (mejores, mejor_puntuacion);
        }
        mejores = posiciones;
        mejor_puntuacion = puntos;
    }
}

pub fn busqueda_aleatoria(
    adn: &[Vec<u8>],
    k: usize,
    repeticiones: usize,
    semilla: u64,
) -> Result<MotivosEncontrados, String> {
    validar(adn, k)?;
    if repeticiones == 0 {
        return Err("Hace falta al menos una repetición.".to_string());
    }
    let mut azar = Azar::nuevo(semilla);
    let mut mejores = Vec::new();
    let mut mejor_puntuacion = usize::MAX;
    let mut puntos_por_intento = Vec::with_capacity(repeticiones);
    for _ in 0..repeticiones {
        let (posiciones, puntos) = mejorar_desde_el_azar(adn, k, &mut azar);
        puntos_por_intento.push(puntos);
        if puntos < mejor_puntuacion {
            mejor_puntuacion = puntos;
            mejores = posiciones;
        }
    }
    Ok(MotivosEncontrados {
        posiciones: mejores,
        puntos: mejor_puntuacion,
        puntos_por_intento,
    })
}

fn cadena_de_gibbs(
    adn: &[Vec<u8>],
    k: usize,
    iteraciones: usize,
    azar: &mut Azar,
) -> (Vec<usize>, usize) {
    let mut posiciones = posiciones_al_azar(adn, k, azar);
    let mut mejores = posiciones.clone();
    let mut mejor_puntuacion = puntuacion(&motivos_en(adn, &mejores, k));
    let mut pesos = Vec::new();
    for _ in 0..iteraciones {
        let i = azar.entero(adn.len());
        let mut otros = motivos_en(adn, &posiciones, k);
        otros.remove(i);
        let perfil = perfil_de_motivos(&otros, PSEUDOCONTEO_DE_LAPLACE);
        pesos.clear();
        for ventana in adn[i].windows(k) {
            pesos.push(probabilidad_kmer(ventana, &perfil));
        }
        posiciones[i] = azar.elegir_con_pesos(&pesos);
        let puntos = puntuacion(&motivos_en(adn, &posiciones, k));
        if puntos < mejor_puntuacion {
            mejor_puntuacion = puntos;
            mejores.copy_from_slice(&posiciones);
        }
    }
    (mejores, mejor_puntuacion)
}

pub fn muestreo_gibbs(
    adn: &[Vec<u8>],
    k: usize,
    iteraciones: usize,
    reinicios: usize,
    semilla: u64,
) -> Result<MotivosEncontrados, String> {
    validar(adn, k)?;
    if reinicios == 0 {
        return Err("Hace falta al menos un reinicio.".to_string());
    }
    let mut azar = Azar::nuevo(semilla);
    let mut mejores = Vec::new();
    let mut mejor_puntuacion = usize::MAX;
    let mut puntos_por_intento = Vec::with_capacity(reinicios);
    for _ in 0..reinicios {
        let (posiciones, puntos) = cadena_de_gibbs(adn, k, iteraciones, &mut azar);
        puntos_por_intento.push(puntos);
        if puntos < mejor_puntuacion {
            mejor_puntuacion = puntos;
            mejores = posiciones;
        }
    }
    Ok(MotivosEncontrados {
        posiciones: mejores,
        puntos: mejor_puntuacion,
        puntos_por_intento,
    })
}
