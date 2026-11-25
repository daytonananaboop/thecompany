use std::collections::VecDeque;

pub type Puntuacion = fn(u8, u8) -> i32;

const INICIO: u8 = 0;
const DIAGONAL: u8 = 1;
const ARRIBA: u8 = 2;
const IZQUIERDA: u8 = 3;
const ABRE_HUECO: u8 = 4;
const ALARGA_HUECO: u8 = 5;
const IMPOSIBLE: i64 = i64::MIN / 4;

pub fn tabla_de_cambio(dinero: usize, monedas: &[usize]) -> Vec<Option<usize>> {
    let mut menos_monedas = vec![None; dinero + 1];
    menos_monedas[0] = Some(0);
    for cantidad in 1..=dinero {
        for &moneda in monedas {
            if moneda == 0 || moneda > cantidad {
                continue;
            }
            if let Some(antes) = menos_monedas[cantidad - moneda] {
                let con_esta_moneda = antes + 1;
                if menos_monedas[cantidad].is_none_or(|actual| con_esta_moneda < actual) {
                    menos_monedas[cantidad] = Some(con_esta_moneda);
                }
            }
        }
    }
    menos_monedas
}

pub fn cambio_minimo(dinero: usize, monedas: &[usize]) -> Option<Vec<usize>> {
    let menos_monedas = tabla_de_cambio(dinero, monedas);
    let mut usadas = Vec::new();
    let mut resto = dinero;
    while resto > 0 {
        let necesarias = menos_monedas[resto]?;
        let elegida = monedas.iter().copied().find(|&moneda| {
            moneda > 0 && moneda <= resto && menos_monedas[resto - moneda] == Some(necesarias - 1)
        })?;
        usadas.push(elegida);
        resto -= elegida;
    }
    Some(usadas)
}

#[derive(Debug, Clone, PartialEq)]
pub struct Recorrido {
    pub longitud: f64,
    pub mejor: Vec<Vec<f64>>,
    pub camino: Vec<(usize, usize)>,
}

pub fn turista_de_manhattan(abajo: &[Vec<f64>], derecha: &[Vec<f64>]) -> Result<Recorrido, String> {
    let filas = derecha.len();
    if filas == 0 {
        return Err("La matriz de calles hacia la derecha no tiene filas.".to_string());
    }
    let columnas = derecha[0].len() + 1;
    let encajan = abajo.len() == filas - 1
        && abajo.iter().all(|fila| fila.len() == columnas)
        && derecha.iter().all(|fila| fila.len() == columnas - 1);
    if !encajan {
        return Err(format!(
            "Las matrices no encajan: con {filas} filas de esquinas y {columnas} columnas, las calles hacia abajo forman {} filas de {columnas} números y las calles hacia la derecha {filas} filas de {} números.",
            filas - 1,
            columnas - 1
        ));
    }
    let mut mejor = vec![vec![0.0; columnas]; filas];
    for i in 1..filas {
        mejor[i][0] = mejor[i - 1][0] + abajo[i - 1][0];
    }
    for j in 1..columnas {
        mejor[0][j] = mejor[0][j - 1] + derecha[0][j - 1];
    }
    for i in 1..filas {
        for j in 1..columnas {
            let bajando = mejor[i - 1][j] + abajo[i - 1][j];
            let girando = mejor[i][j - 1] + derecha[i][j - 1];
            mejor[i][j] = bajando.max(girando);
        }
    }
    let camino = camino_del_turista(&mejor, abajo, derecha);
    Ok(Recorrido {
        longitud: mejor[filas - 1][columnas - 1],
        mejor,
        camino,
    })
}

fn camino_del_turista(
    mejor: &[Vec<f64>],
    abajo: &[Vec<f64>],
    derecha: &[Vec<f64>],
) -> Vec<(usize, usize)> {
    let mut i = mejor.len() - 1;
    let mut j = mejor[0].len() - 1;
    let mut camino = vec![(i, j)];
    while i > 0 || j > 0 {
        let llego_bajando = i > 0
            && (j == 0 || mejor[i - 1][j] + abajo[i - 1][j] >= mejor[i][j - 1] + derecha[i][j - 1]);
        if llego_bajando {
            i -= 1;
        } else {
            j -= 1;
        }
        camino.push((i, j));
    }
    camino.reverse();
    camino
}

pub fn subsecuencia_comun_mas_larga(v: &[u8], w: &[u8]) -> Vec<u8> {
    let mut flechas = Flechas::nuevas(v.len() + 1, w.len() + 1);
    let mut anterior = vec![0usize; w.len() + 1];
    let mut actual = vec![0usize; w.len() + 1];
    for i in 1..=v.len() {
        for j in 1..=w.len() {
            if v[i - 1] == w[j - 1] {
                actual[j] = anterior[j - 1] + 1;
                flechas.poner(i, j, DIAGONAL);
            } else if anterior[j] >= actual[j - 1] {
                actual[j] = anterior[j];
                flechas.poner(i, j, ARRIBA);
            } else {
                actual[j] = actual[j - 1];
                flechas.poner(i, j, IZQUIERDA);
            }
        }
        std::mem::swap(&mut anterior, &mut actual);
    }
    let mut comun = Vec::with_capacity(anterior[w.len()]);
    let (mut i, mut j) = (v.len(), w.len());
    while i > 0 && j > 0 {
        match flechas.leer(i, j) {
            DIAGONAL => {
                comun.push(v[i - 1]);
                i -= 1;
                j -= 1;
            }
            ARRIBA => i -= 1,
            _ => j -= 1,
        }
    }
    comun.reverse();
    comun
}

pub fn orden_topologico(nodos: usize, aristas: &[(usize, usize)]) -> Result<Vec<usize>, String> {
    let mut entrantes = vec![0; nodos];
    let mut salientes = vec![Vec::new(); nodos];
    for &(desde, hasta) in aristas {
        if desde >= nodos || hasta >= nodos {
            return Err(format!(
                "La arista {desde} → {hasta} usa un nodo que no existe: el grafo tiene {nodos}."
            ));
        }
        salientes[desde].push(hasta);
        entrantes[hasta] += 1;
    }
    let mut sin_entrantes: VecDeque<usize> =
        (0..nodos).filter(|&nodo| entrantes[nodo] == 0).collect();
    let mut orden = Vec::with_capacity(nodos);
    while let Some(nodo) = sin_entrantes.pop_front() {
        orden.push(nodo);
        for &siguiente in &salientes[nodo] {
            entrantes[siguiente] -= 1;
            if entrantes[siguiente] == 0 {
                sin_entrantes.push_back(siguiente);
            }
        }
    }
    if orden.len() < nodos {
        return Err("El grafo tiene un ciclo, así que no existe un orden topológico.".to_string());
    }
    Ok(orden)
}

#[derive(Debug, Clone, PartialEq)]
pub struct CaminoLargo {
    pub longitud: f64,
    pub nodos: Vec<usize>,
}

pub fn camino_mas_largo(
    nodos: usize,
    aristas: &[(usize, usize, f64)],
    fuente: usize,
    sumidero: usize,
) -> Result<CaminoLargo, String> {
    if fuente >= nodos || sumidero >= nodos {
        return Err("La fuente y el sumidero deben ser nodos del grafo.".to_string());
    }
    let pares: Vec<(usize, usize)> = aristas
        .iter()
        .map(|&(desde, hasta, _)| (desde, hasta))
        .collect();
    let orden = orden_topologico(nodos, &pares)?;
    let mut entrantes = vec![Vec::new(); nodos];
    for &(desde, hasta, peso) in aristas {
        entrantes[hasta].push((desde, peso));
    }
    let mut mejor: Vec<Option<f64>> = vec![None; nodos];
    let mut anterior: Vec<Option<usize>> = vec![None; nodos];
    mejor[fuente] = Some(0.0);
    for &nodo in &orden {
        for &(desde, peso) in &entrantes[nodo] {
            let Some(llegada) = mejor[desde] else {
                continue;
            };
            let candidato = llegada + peso;
            if mejor[nodo].is_none_or(|actual| candidato > actual) {
                mejor[nodo] = Some(candidato);
                anterior[nodo] = Some(desde);
            }
        }
    }
    let longitud =
        mejor[sumidero].ok_or("No hay ningún camino desde la fuente hasta el sumidero.")?;
    let mut camino = vec![sumidero];
    let mut nodo = sumidero;
    while let Some(previo) = anterior[nodo] {
        camino.push(previo);
        nodo = previo;
    }
    camino.reverse();
    Ok(CaminoLargo {
        longitud,
        nodos: camino,
    })
}

pub fn puntaje_adn(a: u8, b: u8) -> i32 {
    if a == b { 1 } else { -1 }
}

pub fn puntaje_edicion(a: u8, b: u8) -> i32 {
    if a == b { 0 } else { -1 }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alineamiento {
    pub puntaje: i64,
    pub fila_v: Vec<u8>,
    pub fila_w: Vec<u8>,
    pub inicio_v: usize,
    pub inicio_w: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Modo {
    Global,
    Local,
    Ajustado,
    Solapado,
}

pub fn alineamiento_global(v: &[u8], w: &[u8], puntuar: Puntuacion, indel: i64) -> Alineamiento {
    let puntajes = TablaDePuntajes::para(v, w, puntuar);
    alinear(v, w, &puntajes, indel, Modo::Global)
}

pub fn alineamiento_local(v: &[u8], w: &[u8], puntuar: Puntuacion, indel: i64) -> Alineamiento {
    let puntajes = TablaDePuntajes::para(v, w, puntuar);
    alinear(v, w, &puntajes, indel, Modo::Local)
}

pub fn distancia_de_edicion(v: &[u8], w: &[u8]) -> usize {
    let mut anterior: Vec<usize> = (0..=w.len()).collect();
    let mut actual = vec![0; w.len() + 1];
    for i in 1..=v.len() {
        actual[0] = i;
        for j in 1..=w.len() {
            let cambio = if v[i - 1] == w[j - 1] { 0 } else { 1 };
            let sustituir = anterior[j - 1] + cambio;
            let borrar = anterior[j] + 1;
            let insertar = actual[j - 1] + 1;
            actual[j] = sustituir.min(borrar).min(insertar);
        }
        std::mem::swap(&mut anterior, &mut actual);
    }
    anterior[w.len()]
}

pub fn alineamiento_ajustado(v: &[u8], w: &[u8], puntuar: Puntuacion, indel: i64) -> Alineamiento {
    let puntajes = TablaDePuntajes::para(v, w, puntuar);
    alinear(v, w, &puntajes, indel, Modo::Ajustado)
}

pub fn alineamiento_solapado(v: &[u8], w: &[u8], puntuar: Puntuacion, indel: i64) -> Alineamiento {
    let puntajes = TablaDePuntajes::para(v, w, puntuar);
    alinear(v, w, &puntajes, indel, Modo::Solapado)
}

fn alinear(v: &[u8], w: &[u8], puntajes: &TablaDePuntajes, indel: i64, modo: Modo) -> Alineamiento {
    let mut flechas = Flechas::nuevas(v.len() + 1, w.len() + 1);
    let mut anterior = vec![0; w.len() + 1];
    let mut actual = vec![0; w.len() + 1];
    let mut ultima_columna = vec![0; v.len() + 1];
    let mut mejor_local = (0, 0, 0);
    if modo != Modo::Local {
        for j in 1..=w.len() {
            anterior[j] = anterior[j - 1] - indel;
            flechas.poner(0, j, IZQUIERDA);
        }
    }
    ultima_columna[0] = anterior[w.len()];
    for i in 1..=v.len() {
        let letra_v = v[i - 1];
        if modo == Modo::Global {
            actual[0] = anterior[0] - indel;
            flechas.poner(i, 0, ARRIBA);
        }
        for j in 1..=w.len() {
            let diagonal = anterior[j - 1] + puntajes.de(letra_v, w[j - 1]);
            let arriba = anterior[j] - indel;
            let izquierda = actual[j - 1] - indel;
            let (puntaje, flecha) = elegir_mejor(diagonal, arriba, izquierda);
            if modo == Modo::Local && puntaje <= 0 {
                actual[j] = 0;
                flechas.poner(i, j, INICIO);
                continue;
            }
            actual[j] = puntaje;
            flechas.poner(i, j, flecha);
            if modo == Modo::Local && puntaje > mejor_local.0 {
                mejor_local = (puntaje, i, j);
            }
        }
        ultima_columna[i] = actual[w.len()];
        std::mem::swap(&mut anterior, &mut actual);
    }
    let ultima_fila = anterior;
    let (puntaje, fin_v, fin_w) = match modo {
        Modo::Global => (ultima_fila[w.len()], v.len(), w.len()),
        Modo::Local => mejor_local,
        Modo::Ajustado => {
            let i = posicion_del_maximo(&ultima_columna);
            (ultima_columna[i], i, w.len())
        }
        Modo::Solapado => {
            let j = posicion_del_maximo(&ultima_fila);
            (ultima_fila[j], v.len(), j)
        }
    };
    reconstruir(v, w, &flechas, fin_v, fin_w, puntaje)
}

fn elegir_mejor(diagonal: i64, arriba: i64, izquierda: i64) -> (i64, u8) {
    if diagonal >= arriba && diagonal >= izquierda {
        (diagonal, DIAGONAL)
    } else if arriba >= izquierda {
        (arriba, ARRIBA)
    } else {
        (izquierda, IZQUIERDA)
    }
}

fn posicion_del_maximo(valores: &[i64]) -> usize {
    let mut mejor = 0;
    for (posicion, &valor) in valores.iter().enumerate() {
        if valor > valores[mejor] {
            mejor = posicion;
        }
    }
    mejor
}

fn reconstruir(
    v: &[u8],
    w: &[u8],
    flechas: &Flechas,
    fin_v: usize,
    fin_w: usize,
    puntaje: i64,
) -> Alineamiento {
    let mut fila_v = Vec::new();
    let mut fila_w = Vec::new();
    let (mut i, mut j) = (fin_v, fin_w);
    loop {
        match flechas.leer(i, j) {
            DIAGONAL => {
                fila_v.push(v[i - 1]);
                fila_w.push(w[j - 1]);
                i -= 1;
                j -= 1;
            }
            ARRIBA => {
                fila_v.push(v[i - 1]);
                fila_w.push(b'-');
                i -= 1;
            }
            IZQUIERDA => {
                fila_v.push(b'-');
                fila_w.push(w[j - 1]);
                j -= 1;
            }
            _ => break,
        }
    }
    fila_v.reverse();
    fila_w.reverse();
    Alineamiento {
        puntaje,
        fila_v,
        fila_w,
        inicio_v: i,
        inicio_w: j,
    }
}

struct Flechas {
    columnas: usize,
    celdas: Vec<u8>,
}

impl Flechas {
    fn nuevas(filas: usize, columnas: usize) -> Flechas {
        Flechas {
            columnas,
            celdas: vec![INICIO; filas * columnas],
        }
    }

    fn poner(&mut self, i: usize, j: usize, flecha: u8) {
        self.celdas[i * self.columnas + j] = flecha;
    }

    fn leer(&self, i: usize, j: usize) -> u8 {
        self.celdas[i * self.columnas + j]
    }
}

struct TablaDePuntajes {
    valores: Vec<i64>,
}

impl TablaDePuntajes {
    fn para(v: &[u8], w: &[u8], puntuar: Puntuacion) -> TablaDePuntajes {
        let mut valores = vec![0; 256 * 256];
        let letras_de_w = letras_distintas(w);
        for a in letras_distintas(v) {
            for &b in &letras_de_w {
                valores[a as usize * 256 + b as usize] = puntuar(a, b) as i64;
            }
        }
        TablaDePuntajes { valores }
    }

    fn de(&self, a: u8, b: u8) -> i64 {
        self.valores[a as usize * 256 + b as usize]
    }
}

fn letras_distintas(secuencia: &[u8]) -> Vec<u8> {
    let mut aparece = vec![false; 256];
    for &letra in secuencia {
        aparece[letra as usize] = true;
    }
    (0..=255u8)
        .filter(|&letra| aparece[letra as usize])
        .collect()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Nivel {
    Inferior,
    Medio,
    Superior,
}

pub fn alineamiento_afin(
    v: &[u8],
    w: &[u8],
    puntuar: Puntuacion,
    apertura: i64,
    extension: i64,
) -> Alineamiento {
    let extension = extension.min(apertura);
    let puntajes = TablaDePuntajes::para(v, w, puntuar);
    let columnas = w.len() + 1;
    let mut flechas_inferior = Flechas::nuevas(v.len() + 1, columnas);
    let mut flechas_medio = Flechas::nuevas(v.len() + 1, columnas);
    let mut flechas_superior = Flechas::nuevas(v.len() + 1, columnas);
    let mut inferior_anterior = vec![IMPOSIBLE; columnas];
    let mut medio_anterior = vec![IMPOSIBLE; columnas];
    let mut inferior = vec![IMPOSIBLE; columnas];
    let mut medio = vec![IMPOSIBLE; columnas];
    let mut superior = vec![IMPOSIBLE; columnas];
    medio[0] = 0;
    for j in 1..columnas {
        let (valor, flecha) = elegir_hueco(superior[j - 1] - extension, medio[j - 1] - apertura);
        superior[j] = valor;
        flechas_superior.poner(0, j, flecha);
        medio[j] = valor;
        flechas_medio.poner(0, j, IZQUIERDA);
    }
    for i in 1..=v.len() {
        std::mem::swap(&mut inferior, &mut inferior_anterior);
        std::mem::swap(&mut medio, &mut medio_anterior);
        let (valor, flecha) = elegir_hueco(
            inferior_anterior[0] - extension,
            medio_anterior[0] - apertura,
        );
        inferior[0] = valor;
        flechas_inferior.poner(i, 0, flecha);
        medio[0] = valor;
        flechas_medio.poner(i, 0, ARRIBA);
        superior[0] = IMPOSIBLE;
        for j in 1..columnas {
            let (valor_inferior, flecha_inferior) = elegir_hueco(
                inferior_anterior[j] - extension,
                medio_anterior[j] - apertura,
            );
            inferior[j] = valor_inferior;
            flechas_inferior.poner(i, j, flecha_inferior);
            let (valor_superior, flecha_superior) =
                elegir_hueco(superior[j - 1] - extension, medio[j - 1] - apertura);
            superior[j] = valor_superior;
            flechas_superior.poner(i, j, flecha_superior);
            let diagonal = medio_anterior[j - 1] + puntajes.de(v[i - 1], w[j - 1]);
            let (valor_medio, flecha_medio) =
                elegir_mejor(diagonal, valor_inferior, valor_superior);
            medio[j] = valor_medio;
            flechas_medio.poner(i, j, flecha_medio);
        }
    }
    let mut fila_v = Vec::new();
    let mut fila_w = Vec::new();
    let (mut i, mut j) = (v.len(), w.len());
    let mut nivel = Nivel::Medio;
    while i > 0 || j > 0 {
        match nivel {
            Nivel::Medio => match flechas_medio.leer(i, j) {
                DIAGONAL => {
                    fila_v.push(v[i - 1]);
                    fila_w.push(w[j - 1]);
                    i -= 1;
                    j -= 1;
                }
                ARRIBA => nivel = Nivel::Inferior,
                _ => nivel = Nivel::Superior,
            },
            Nivel::Inferior => {
                fila_v.push(v[i - 1]);
                fila_w.push(b'-');
                if flechas_inferior.leer(i, j) == ABRE_HUECO {
                    nivel = Nivel::Medio;
                }
                i -= 1;
            }
            Nivel::Superior => {
                fila_v.push(b'-');
                fila_w.push(w[j - 1]);
                if flechas_superior.leer(i, j) == ABRE_HUECO {
                    nivel = Nivel::Medio;
                }
                j -= 1;
            }
        }
    }
    fila_v.reverse();
    fila_w.reverse();
    Alineamiento {
        puntaje: medio[w.len()],
        fila_v,
        fila_w,
        inicio_v: 0,
        inicio_w: 0,
    }
}

fn elegir_hueco(alargar: i64, abrir: i64) -> (i64, u8) {
    if alargar >= abrir {
        (alargar, ALARGA_HUECO)
    } else {
        (abrir, ABRE_HUECO)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AristaCentral {
    pub desde: (usize, usize),
    pub hasta: (usize, usize),
    pub longitudes: Vec<i64>,
}

pub fn arista_central(
    v: &[u8],
    w: &[u8],
    puntuar: Puntuacion,
    indel: i64,
) -> Option<AristaCentral> {
    let puntajes = TablaDePuntajes::para(v, w, puntuar);
    buscar_arista_central(v, w, &puntajes, indel)
}

fn buscar_arista_central(
    v: &[u8],
    w: &[u8],
    puntajes: &TablaDePuntajes,
    indel: i64,
) -> Option<AristaCentral> {
    let centro = w.len() / 2;
    let desde_la_fuente = puntajes_desde_la_fuente(v, &w[..centro], puntajes, indel);
    let (hasta_el_sumidero, hasta_el_sumidero_a_la_derecha) =
        puntajes_hasta_el_sumidero(v, &w[centro..], puntajes, indel);
    let longitudes: Vec<i64> = desde_la_fuente
        .iter()
        .zip(&hasta_el_sumidero)
        .map(|(antes, despues)| antes + despues)
        .collect();
    let fila = posicion_del_maximo(&longitudes);
    let resto = hasta_el_sumidero[fila];
    let hay_columna_siguiente = centro < w.len();
    let hay_fila_siguiente = fila < v.len();
    let hasta = if hay_fila_siguiente
        && hay_columna_siguiente
        && puntajes.de(v[fila], w[centro]) + hasta_el_sumidero_a_la_derecha[fila + 1] == resto
    {
        (fila + 1, centro + 1)
    } else if hay_columna_siguiente && hasta_el_sumidero_a_la_derecha[fila] - indel == resto {
        (fila, centro + 1)
    } else if hay_fila_siguiente && hasta_el_sumidero[fila + 1] - indel == resto {
        (fila + 1, centro)
    } else {
        return None;
    };
    Some(AristaCentral {
        desde: (fila, centro),
        hasta,
        longitudes,
    })
}

fn puntajes_desde_la_fuente(
    v: &[u8],
    w: &[u8],
    puntajes: &TablaDePuntajes,
    indel: i64,
) -> Vec<i64> {
    let mut columna: Vec<i64> = (0..=v.len()).map(|i| -(i as i64) * indel).collect();
    let mut izquierda = vec![0; v.len() + 1];
    for &letra_w in w {
        std::mem::swap(&mut columna, &mut izquierda);
        columna[0] = izquierda[0] - indel;
        for i in 1..=v.len() {
            let diagonal = izquierda[i - 1] + puntajes.de(v[i - 1], letra_w);
            let horizontal = izquierda[i] - indel;
            let vertical = columna[i - 1] - indel;
            columna[i] = diagonal.max(horizontal).max(vertical);
        }
    }
    columna
}

fn puntajes_hasta_el_sumidero(
    v: &[u8],
    w: &[u8],
    puntajes: &TablaDePuntajes,
    indel: i64,
) -> (Vec<i64>, Vec<i64>) {
    let n = v.len();
    let mut columna: Vec<i64> = (0..=n).map(|i| -((n - i) as i64) * indel).collect();
    let mut derecha = vec![0; n + 1];
    for &letra_w in w.iter().rev() {
        std::mem::swap(&mut columna, &mut derecha);
        columna[n] = derecha[n] - indel;
        for i in (0..n).rev() {
            let diagonal = derecha[i + 1] + puntajes.de(v[i], letra_w);
            let horizontal = derecha[i] - indel;
            let vertical = columna[i + 1] - indel;
            columna[i] = diagonal.max(horizontal).max(vertical);
        }
    }
    (columna, derecha)
}

pub fn alineamiento_espacio_lineal(
    v: &[u8],
    w: &[u8],
    puntuar: Puntuacion,
    indel: i64,
) -> Alineamiento {
    let puntajes = TablaDePuntajes::para(v, w, puntuar);
    let mut fila_v = Vec::with_capacity(v.len() + w.len());
    let mut fila_w = Vec::with_capacity(v.len() + w.len());
    alinear_por_mitades(v, w, &puntajes, indel, &mut fila_v, &mut fila_w);
    let puntaje = puntuar_alineamiento(&fila_v, &fila_w, puntuar, indel);
    Alineamiento {
        puntaje,
        fila_v,
        fila_w,
        inicio_v: 0,
        inicio_w: 0,
    }
}

fn alinear_por_mitades(
    v: &[u8],
    w: &[u8],
    puntajes: &TablaDePuntajes,
    indel: i64,
    fila_v: &mut Vec<u8>,
    fila_w: &mut Vec<u8>,
) {
    if v.len() <= 1 || w.len() <= 1 {
        let directo = alinear(v, w, puntajes, indel, Modo::Global);
        fila_v.extend_from_slice(&directo.fila_v);
        fila_w.extend_from_slice(&directo.fila_w);
        return;
    }
    let Some(arista) = buscar_arista_central(v, w, puntajes, indel) else {
        return;
    };
    let (i, j) = arista.desde;
    let (siguiente_i, siguiente_j) = arista.hasta;
    alinear_por_mitades(&v[..i], &w[..j], puntajes, indel, fila_v, fila_w);
    fila_v.push(if siguiente_i > i { v[i] } else { b'-' });
    fila_w.push(if siguiente_j > j { w[j] } else { b'-' });
    alinear_por_mitades(
        &v[siguiente_i..],
        &w[siguiente_j..],
        puntajes,
        indel,
        fila_v,
        fila_w,
    );
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlineamientoTriple {
    pub puntaje: usize,
    pub fila_u: Vec<u8>,
    pub fila_v: Vec<u8>,
    pub fila_w: Vec<u8>,
}

const MOVIMIENTOS: &[(usize, usize, usize)] = &[
    (1, 1, 1),
    (1, 1, 0),
    (1, 0, 1),
    (0, 1, 1),
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
];

pub fn subsecuencia_comun_de_tres(u: &[u8], v: &[u8], w: &[u8]) -> AlineamientoTriple {
    let mut cubo = Cubo::de_ceros(u.len() + 1, v.len() + 1, w.len() + 1);
    for i in 0..=u.len() {
        for j in 0..=v.len() {
            for k in 0..=w.len() {
                let mut mejor = 0;
                for &movimiento in MOVIMIENTOS {
                    if let Some(llegada) = llegar_por(&cubo, (u, v, w), (i, j, k), movimiento) {
                        mejor = mejor.max(llegada);
                    }
                }
                cubo.anotar(i, j, k, mejor);
            }
        }
    }
    let mut fila_u = Vec::new();
    let mut fila_v = Vec::new();
    let mut fila_w = Vec::new();
    let (mut i, mut j, mut k) = (u.len(), v.len(), w.len());
    while i > 0 || j > 0 || k > 0 {
        let actual = cubo.mirar(i, j, k);
        let Some(&(paso_u, paso_v, paso_w)) = MOVIMIENTOS.iter().find(|&&movimiento| {
            llegar_por(&cubo, (u, v, w), (i, j, k), movimiento) == Some(actual)
        }) else {
            break;
        };
        fila_u.push(if paso_u == 1 { u[i - 1] } else { b'-' });
        fila_v.push(if paso_v == 1 { v[j - 1] } else { b'-' });
        fila_w.push(if paso_w == 1 { w[k - 1] } else { b'-' });
        i -= paso_u;
        j -= paso_v;
        k -= paso_w;
    }
    fila_u.reverse();
    fila_v.reverse();
    fila_w.reverse();
    AlineamientoTriple {
        puntaje: cubo.mirar(u.len(), v.len(), w.len()) as usize,
        fila_u,
        fila_v,
        fila_w,
    }
}

fn llegar_por(
    cubo: &Cubo,
    secuencias: (&[u8], &[u8], &[u8]),
    celda: (usize, usize, usize),
    movimiento: (usize, usize, usize),
) -> Option<u32> {
    let (u, v, w) = secuencias;
    let (i, j, k) = celda;
    let (paso_u, paso_v, paso_w) = movimiento;
    if i < paso_u || j < paso_v || k < paso_w {
        return None;
    }
    let antes = cubo.mirar(i - paso_u, j - paso_v, k - paso_w);
    let tres_iguales = movimiento == (1, 1, 1) && u[i - 1] == v[j - 1] && v[j - 1] == w[k - 1];
    Some(if tres_iguales { antes + 1 } else { antes })
}

struct Cubo {
    ancho: usize,
    fondo: usize,
    celdas: Vec<u32>,
}

impl Cubo {
    fn de_ceros(alto: usize, ancho: usize, fondo: usize) -> Cubo {
        Cubo {
            ancho,
            fondo,
            celdas: vec![0; alto * ancho * fondo],
        }
    }

    fn mirar(&self, i: usize, j: usize, k: usize) -> u32 {
        self.celdas[(i * self.ancho + j) * self.fondo + k]
    }

    fn anotar(&mut self, i: usize, j: usize, k: usize, valor: u32) {
        self.celdas[(i * self.ancho + j) * self.fondo + k] = valor;
    }
}

pub fn puntuar_alineamiento(fila_v: &[u8], fila_w: &[u8], puntuar: Puntuacion, indel: i64) -> i64 {
    let mut puntaje = 0;
    for (&a, &b) in fila_v.iter().zip(fila_w) {
        if a == b'-' || b == b'-' {
            puntaje -= indel;
        } else {
            puntaje += puntuar(a, b) as i64;
        }
    }
    puntaje
}

pub fn puntuar_alineamiento_afin(
    fila_v: &[u8],
    fila_w: &[u8],
    puntuar: Puntuacion,
    apertura: i64,
    extension: i64,
) -> i64 {
    let extension = extension.min(apertura);
    let mut puntaje = 0;
    let mut hueco_abierto_en_v = false;
    let mut hueco_abierto_en_w = false;
    for (&a, &b) in fila_v.iter().zip(fila_w) {
        if a == b'-' {
            puntaje -= if hueco_abierto_en_v {
                extension
            } else {
                apertura
            };
            hueco_abierto_en_v = true;
            hueco_abierto_en_w = false;
        } else if b == b'-' {
            puntaje -= if hueco_abierto_en_w {
                extension
            } else {
                apertura
            };
            hueco_abierto_en_w = true;
            hueco_abierto_en_v = false;
        } else {
            puntaje += puntuar(a, b) as i64;
            hueco_abierto_en_v = false;
            hueco_abierto_en_w = false;
        }
    }
    puntaje
}
