use rayon::prelude::*;
use std::collections::HashSet;

pub const BASES: &[u8; 4] = b"ACGT";

const PROBABILIDAD_DE_ARRANQUE: f64 = 0.5;
const GAMMA_INICIAL: f64 = 0.5;
const GAMMA_MINIMO: f64 = 1e-6;
pub const ITERACIONES_DE_PRUEBA: usize = 3;
pub const ITERACIONES_MAXIMAS: usize = 200;
pub const CAMBIO_MINIMO: f64 = 1e-6;
pub const PSEUDOCONTEO_MINIMO: f64 = 1e-6;
pub const PSEUDOCONTEO_MAXIMO: f64 = 1e6;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Modelo {
    Oops,
    Zoops,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Opciones {
    pub k: usize,
    pub modelo: Modelo,
    pub ambas_hebras: bool,
    pub pseudoconteo: f64,
    pub arranques: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Parametros {
    pub perfil: Vec<[f64; 4]>,
    pub gamma: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Recorrido {
    pub parametros: Parametros,
    pub historial: Vec<f64>,
    pub posteriores: Vec<Vec<f64>>,
    pub iteraciones: usize,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sitio {
    pub posicion: usize,
    pub hebra: u8,
    pub probabilidad: f64,
    pub probabilidad_de_sitio: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MotivoEm {
    pub perfil: Vec<[f64; 4]>,
    pub gamma: f64,
    pub fondo: [f64; 4],
    pub consenso: Vec<u8>,
    pub informacion: Vec<f64>,
    pub verosimilitud: f64,
    pub historial: Vec<f64>,
    pub sitios: Vec<Sitio>,
    pub arranque: Vec<u8>,
    pub arranques_probados: usize,
    pub iteraciones: usize,
}

fn codigo_de(base: u8) -> Option<u8> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

pub fn codificar(adn: &[Vec<u8>]) -> Result<Vec<Vec<u8>>, String> {
    let mut codigos = Vec::with_capacity(adn.len());
    for (numero, texto) in adn.iter().enumerate() {
        let mut fila = Vec::with_capacity(texto.len());
        for &base in texto {
            let codigo = codigo_de(base).ok_or_else(|| {
                format!(
                    "La secuencia {} contiene «{}»; solo se admiten A, C, G y T.",
                    numero + 1,
                    base as char
                )
            })?;
            fila.push(codigo);
        }
        codigos.push(fila);
    }
    Ok(codigos)
}

fn validar(adn: &[Vec<u8>], opciones: &Opciones) -> Result<(), String> {
    if adn.is_empty() {
        return Err("Hace falta al menos una secuencia.".to_string());
    }
    if opciones.k == 0 {
        return Err("k debe ser al menos 1.".to_string());
    }
    for (numero, texto) in adn.iter().enumerate() {
        if texto.len() < opciones.k {
            return Err(format!(
                "La secuencia {} mide {} pb, menos que k = {}.",
                numero + 1,
                texto.len(),
                opciones.k
            ));
        }
    }
    if !(PSEUDOCONTEO_MINIMO..=PSEUDOCONTEO_MAXIMO).contains(&opciones.pseudoconteo) {
        return Err(format!(
            "El pseudoconteo debe estar entre {PSEUDOCONTEO_MINIMO} y {PSEUDOCONTEO_MAXIMO}."
        ));
    }
    if opciones.arranques == 0 {
        return Err("Hace falta al menos un punto de arranque.".to_string());
    }
    Ok(())
}

fn hebras(opciones: &Opciones) -> usize {
    if opciones.ambas_hebras { 2 } else { 1 }
}

pub fn frecuencias_de_fondo(codigos: &[Vec<u8>], ambas_hebras: bool) -> [f64; 4] {
    let mut conteos = [1.0; 4];
    for fila in codigos {
        for &codigo in fila {
            conteos[codigo as usize] += 1.0;
        }
    }
    if ambas_hebras {
        let a_o_t = (conteos[0] + conteos[3]) / 2.0;
        let c_o_g = (conteos[1] + conteos[2]) / 2.0;
        conteos = [a_o_t, c_o_g, c_o_g, a_o_t];
    }
    let total: f64 = conteos.iter().sum();
    conteos.map(|conteo| conteo / total)
}

pub fn arranques(codigos: &[Vec<u8>], k: usize, maximo: usize) -> Vec<&[u8]> {
    let mut vistos: HashSet<&[u8]> = HashSet::new();
    let mut semillas = Vec::new();
    for fila in codigos {
        for ventana in fila.windows(k) {
            if semillas.len() == maximo {
                return semillas;
            }
            if vistos.insert(ventana) {
                semillas.push(ventana);
            }
        }
    }
    semillas
}

pub fn parametros_de_arranque(semilla: &[u8], modelo: Modelo) -> Parametros {
    let resto = (1.0 - PROBABILIDAD_DE_ARRANQUE) / 3.0;
    let mut perfil = Vec::with_capacity(semilla.len());
    for &codigo in semilla {
        let mut columna = [resto; 4];
        columna[codigo as usize] = PROBABILIDAD_DE_ARRANQUE;
        perfil.push(columna);
    }
    let gamma = match modelo {
        Modelo::Oops => 1.0,
        Modelo::Zoops => GAMMA_INICIAL,
    };
    Parametros { perfil, gamma }
}

pub fn pesos_logaritmicos(perfil: &[[f64; 4]], fondo: &[f64; 4]) -> Vec<[f64; 4]> {
    let mut pesos = Vec::with_capacity(perfil.len());
    for columna in perfil {
        let mut peso = [0.0; 4];
        for base in 0..4 {
            peso[base] = (columna[base] / fondo[base]).ln();
        }
        pesos.push(peso);
    }
    pesos
}

fn puntaje_directo(pesos: &[[f64; 4]], ventana: &[u8]) -> f64 {
    let mut puntaje = 0.0;
    for (columna, &codigo) in pesos.iter().zip(ventana) {
        puntaje += columna[(codigo & 3) as usize];
    }
    puntaje
}

fn puntaje_inverso(pesos: &[[f64; 4]], ventana: &[u8]) -> f64 {
    let mut puntaje = 0.0;
    for (columna, &codigo) in pesos.iter().zip(ventana.iter().rev()) {
        puntaje += columna[((3 - codigo) & 3) as usize];
    }
    puntaje
}

pub fn puntajes_de_ventanas(
    pesos: &[[f64; 4]],
    fila: &[u8],
    dos_hebras: bool,
    puntajes: &mut Vec<f64>,
) {
    puntajes.clear();
    for ventana in fila.windows(pesos.len()) {
        puntajes.push(puntaje_directo(pesos, ventana));
        if dos_hebras {
            puntajes.push(puntaje_inverso(pesos, ventana));
        }
    }
}

fn sumar_logaritmos(a: f64, b: f64) -> f64 {
    if b == f64::NEG_INFINITY {
        return a;
    }
    let mayor = a.max(b);
    mayor + ((a - mayor).exp() + (b - mayor).exp()).ln()
}

pub fn posteriores_de_secuencia(z: &mut [f64], gamma: f64, modelo: Modelo) -> f64 {
    let maximo = z.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut suma = 0.0;
    for valor in z.iter_mut() {
        *valor = (*valor - maximo).exp();
        suma += *valor;
    }
    let candidatos = z.len() as f64;
    let ln_con_sitio = (gamma / candidatos).ln() + maximo + suma.ln();
    let ln_total = match modelo {
        Modelo::Oops => ln_con_sitio,
        Modelo::Zoops => sumar_logaritmos(ln_con_sitio, (1.0 - gamma).ln()),
    };
    let escala = (ln_con_sitio - ln_total).exp() / suma;
    for valor in z.iter_mut() {
        *valor *= escala;
    }
    ln_total
}

pub fn paso_e(
    codigos: &[Vec<u8>],
    parametros: &Parametros,
    fondo: &[f64; 4],
    opciones: &Opciones,
    posteriores: &mut [Vec<f64>],
) -> f64 {
    let pesos = pesos_logaritmicos(&parametros.perfil, fondo);
    let por_secuencia: Vec<f64> = codigos
        .par_iter()
        .zip(posteriores.par_iter_mut())
        .map(|(fila, z)| {
            puntajes_de_ventanas(&pesos, fila, opciones.ambas_hebras, z);
            posteriores_de_secuencia(z, parametros.gamma, opciones.modelo)
        })
        .collect();
    let mut verosimilitud = 0.0;
    for valor in por_secuencia {
        verosimilitud += valor;
    }
    verosimilitud
}

fn sumar_directo(conteos: &mut [[f64; 4]], ventana: &[u8], peso: f64) {
    for (columna, &codigo) in conteos.iter_mut().zip(ventana) {
        columna[(codigo & 3) as usize] += peso;
    }
}

fn sumar_inverso(conteos: &mut [[f64; 4]], ventana: &[u8], peso: f64) {
    for (columna, &codigo) in conteos.iter_mut().zip(ventana.iter().rev()) {
        columna[((3 - codigo) & 3) as usize] += peso;
    }
}

pub fn conteos_esperados(
    codigos: &[Vec<u8>],
    posteriores: &[Vec<f64>],
    opciones: &Opciones,
) -> (Vec<[f64; 4]>, f64) {
    let k = opciones.k;
    let cuantas_hebras = hebras(opciones);
    let mut conteos = vec![[0.0; 4]; k];
    let mut sitios_esperados = 0.0;
    for (fila, z) in codigos.iter().zip(posteriores) {
        for (candidato, &peso) in z.iter().enumerate() {
            sitios_esperados += peso;
            let inicio = candidato / cuantas_hebras;
            let ventana = &fila[inicio..inicio + k];
            if candidato % cuantas_hebras == 0 {
                sumar_directo(&mut conteos, ventana, peso);
            } else {
                sumar_inverso(&mut conteos, ventana, peso);
            }
        }
    }
    (conteos, sitios_esperados)
}

pub fn paso_m(
    codigos: &[Vec<u8>],
    posteriores: &[Vec<f64>],
    fondo: &[f64; 4],
    opciones: &Opciones,
) -> Parametros {
    let (conteos, sitios_esperados) = conteos_esperados(codigos, posteriores, opciones);
    let beta = opciones.pseudoconteo;
    let mut perfil = Vec::with_capacity(conteos.len());
    for columna in &conteos {
        let total: f64 = columna.iter().sum::<f64>() + beta;
        let mut probabilidades = [0.0; 4];
        for base in 0..4 {
            probabilidades[base] = (columna[base] + beta * fondo[base]) / total;
        }
        perfil.push(probabilidades);
    }
    let gamma = match opciones.modelo {
        Modelo::Oops => 1.0,
        Modelo::Zoops => (sitios_esperados / codigos.len() as f64).clamp(GAMMA_MINIMO, 1.0),
    };
    Parametros { perfil, gamma }
}

pub fn cambio_entre(antes: &Parametros, despues: &Parametros) -> f64 {
    let mut mayor = (antes.gamma - despues.gamma).abs();
    for (columna_antes, columna_despues) in antes.perfil.iter().zip(&despues.perfil) {
        for base in 0..4 {
            mayor = mayor.max((columna_antes[base] - columna_despues[base]).abs());
        }
    }
    mayor
}

pub fn posteriores_vacios(codigos: &[Vec<u8>], opciones: &Opciones) -> Vec<Vec<f64>> {
    codigos
        .iter()
        .map(|fila| Vec::with_capacity((fila.len() - opciones.k + 1) * hebras(opciones)))
        .collect()
}

pub fn iterar(
    codigos: &[Vec<u8>],
    inicial: Parametros,
    fondo: &[f64; 4],
    opciones: &Opciones,
    maximo: usize,
) -> Recorrido {
    let mut posteriores = posteriores_vacios(codigos, opciones);
    let mut parametros = inicial;
    let mut historial = Vec::new();
    let mut iteraciones = 0;
    while iteraciones < maximo {
        let verosimilitud = paso_e(codigos, &parametros, fondo, opciones, &mut posteriores);
        historial.push(verosimilitud / std::f64::consts::LN_2);
        let nuevos = paso_m(codigos, &posteriores, fondo, opciones);
        let cambio = cambio_entre(&parametros, &nuevos);
        parametros = nuevos;
        iteraciones += 1;
        if cambio < CAMBIO_MINIMO {
            break;
        }
    }
    let verosimilitud = paso_e(codigos, &parametros, fondo, opciones, &mut posteriores);
    historial.push(verosimilitud / std::f64::consts::LN_2);
    Recorrido {
        parametros,
        historial,
        posteriores,
        iteraciones,
    }
}

pub fn probar_arranque(
    codigos: &[Vec<u8>],
    semilla: &[u8],
    fondo: &[f64; 4],
    opciones: &Opciones,
) -> f64 {
    let inicial = parametros_de_arranque(semilla, opciones.modelo);
    let recorrido = iterar(codigos, inicial, fondo, opciones, ITERACIONES_DE_PRUEBA);
    recorrido
        .historial
        .last()
        .copied()
        .unwrap_or(f64::NEG_INFINITY)
}

pub fn mejor_arranque(puntajes: &[f64]) -> usize {
    let mut mejor = 0;
    for (indice, &puntaje) in puntajes.iter().enumerate() {
        if puntaje > puntajes[mejor] {
            mejor = indice;
        }
    }
    mejor
}

pub fn mejor_sitio(z: &[f64], dos_hebras: bool) -> Sitio {
    let cuantas_hebras = if dos_hebras { 2 } else { 1 };
    let mut sitio = Sitio {
        posicion: 0,
        hebra: b'+',
        probabilidad: -1.0,
        probabilidad_de_sitio: z.iter().sum(),
    };
    for (posicion, pareja) in z.chunks(cuantas_hebras).enumerate() {
        let probabilidad: f64 = pareja.iter().sum();
        if probabilidad > sitio.probabilidad {
            sitio.posicion = posicion;
            sitio.probabilidad = probabilidad;
            sitio.hebra = if pareja.len() == 2 && pareja[1] > pareja[0] {
                b'-'
            } else {
                b'+'
            };
        }
    }
    sitio
}

pub fn consenso_del_perfil(perfil: &[[f64; 4]]) -> Vec<u8> {
    let mut consenso = Vec::with_capacity(perfil.len());
    for columna in perfil {
        let mut mayor = 0;
        for base in 1..4 {
            if columna[base] > columna[mayor] {
                mayor = base;
            }
        }
        consenso.push(BASES[mayor]);
    }
    consenso
}

pub fn informacion_por_columna(perfil: &[[f64; 4]]) -> Vec<f64> {
    let mut bits = Vec::with_capacity(perfil.len());
    for columna in perfil {
        let mut entropia = 0.0;
        for &p in columna {
            if p > 0.0 {
                entropia -= p * p.log2();
            }
        }
        bits.push(2.0 - entropia);
    }
    bits
}

pub fn decodificar(codigos: &[u8]) -> Vec<u8> {
    codigos
        .iter()
        .map(|&codigo| BASES[codigo as usize])
        .collect()
}

pub fn meme(adn: &[Vec<u8>], opciones: &Opciones) -> Result<MotivoEm, String> {
    validar(adn, opciones)?;
    let codigos = codificar(adn)?;
    let fondo = frecuencias_de_fondo(&codigos, opciones.ambas_hebras);
    let semillas = arranques(&codigos, opciones.k, opciones.arranques);
    let puntajes: Vec<f64> = semillas
        .par_iter()
        .map(|semilla| probar_arranque(&codigos, semilla, &fondo, opciones))
        .collect();
    let mejor = mejor_arranque(&puntajes);
    let inicial = parametros_de_arranque(semillas[mejor], opciones.modelo);
    let recorrido = iterar(&codigos, inicial, &fondo, opciones, ITERACIONES_MAXIMAS);
    let mut sitios = Vec::with_capacity(codigos.len());
    for z in &recorrido.posteriores {
        sitios.push(mejor_sitio(z, opciones.ambas_hebras));
    }
    let perfil = recorrido.parametros.perfil;
    Ok(MotivoEm {
        consenso: consenso_del_perfil(&perfil),
        informacion: informacion_por_columna(&perfil),
        perfil,
        gamma: recorrido.parametros.gamma,
        fondo,
        verosimilitud: recorrido.historial.last().copied().unwrap_or(0.0),
        historial: recorrido.historial,
        sitios,
        arranque: decodificar(semillas[mejor]),
        arranques_probados: semillas.len(),
        iteraciones: recorrido.iteraciones,
    })
}
