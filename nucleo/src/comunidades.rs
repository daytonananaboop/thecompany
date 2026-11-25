use crate::agrupamiento::{Medida, distancia_al_cuadrado, revisar_puntos};
use crate::azar::Azar;
use rayon::prelude::*;
use std::collections::BTreeMap;

pub const PESO_MINIMO_SNN: f64 = 1.0 / 15.0;
const MEJORA_MINIMA: f64 = 1e-12;
pub const MAXIMO_PASADAS: usize = 1000;
pub const TOLERANCIA_DE_ENTROPIA: f64 = 1e-5;
pub const INTENTOS_DE_BETA: usize = 100;
const PROBABILIDAD_MINIMA: f64 = 1e-12;
const ESCALA_INICIAL: f64 = 1e-4;
pub const MOMENTO_INICIAL: f64 = 0.5;
pub const MOMENTO_FINAL: f64 = 0.8;
pub const ITERACIONES_EXAGERADAS: usize = 250;
const GANANCIA_MINIMA: f64 = 0.01;
pub const CADA_CUANTO_MEDIR: usize = 25;
const PARTE_DESPRECIABLE: f64 = 1e-12;

#[derive(Debug, Clone, PartialEq)]
pub struct Grafo {
    pub vecinos: Vec<Vec<(usize, f64)>>,
    pub lazos: Vec<f64>,
    pub grados: Vec<f64>,
    pub peso_total: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Particion {
    pub comunidad: Vec<usize>,
    pub grado_de_comunidad: Vec<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Movimiento {
    pub nodo: usize,
    pub desde: usize,
    pub hacia: usize,
    pub ganancia: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Nivel {
    pub nodos: usize,
    pub movimientos: usize,
    pub comunidades: usize,
    pub modularidad: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Louvain {
    pub comunidad: Vec<usize>,
    pub cantidad: usize,
    pub modularidad: f64,
    pub niveles: Vec<Nivel>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Comunidades {
    pub vecinos: Vec<Vec<usize>>,
    pub grafo: Grafo,
    pub louvain: Louvain,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OpcionesTsne {
    pub perplejidad: f64,
    pub exageracion: f64,
    pub iteraciones: usize,
    pub semilla: u64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Calibracion {
    pub beta: f64,
    pub entropia: f64,
    pub intentos: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Tsne {
    pub coordenadas: Vec<[f64; 2]>,
    pub divergencia: f64,
    pub historia: Vec<(usize, f64)>,
    pub sigmas: Vec<f64>,
    pub tasa: f64,
}

pub fn filas_para_medir(puntos: &[Vec<f64>], medida: Medida) -> Vec<Vec<f64>> {
    match medida {
        Medida::Euclidea => puntos.to_vec(),
        Medida::Pearson => {
            let mut filas: Vec<Vec<f64>> = puntos
                .iter()
                .map(|punto| centrar_y_normalizar(punto))
                .collect();
            let es_plana = |fila: &Vec<f64>| fila.iter().all(|&valor| valor == 0.0);
            if filas.iter().any(es_plana) {
                for fila in &mut filas {
                    fila.push(if es_plana(fila) { 1.0 } else { 0.0 });
                }
            }
            filas
        }
    }
}

fn centrar_y_normalizar(valores: &[f64]) -> Vec<f64> {
    let promedio = valores.iter().sum::<f64>() / valores.len() as f64;
    let mut centrados: Vec<f64> = valores.iter().map(|valor| valor - promedio).collect();
    let magnitud = valores
        .iter()
        .fold(0.0f64, |mayor, valor| mayor.max(valor.abs()));
    let largo = centrados
        .iter()
        .map(|valor| valor * valor)
        .sum::<f64>()
        .sqrt();
    if largo <= magnitud * PARTE_DESPRECIABLE || largo == 0.0 {
        return vec![0.0; valores.len()];
    }
    for valor in &mut centrados {
        *valor /= largo;
    }
    centrados
}

pub fn distancia_medida(cuadrado: f64, medida: Medida) -> f64 {
    match medida {
        Medida::Euclidea => cuadrado.sqrt(),
        Medida::Pearson => cuadrado / 2.0,
    }
}

pub fn distancias_desde(filas: &[Vec<f64>], punto: usize) -> Vec<(f64, usize)> {
    let mut distancias = Vec::with_capacity(filas.len());
    llenar_distancias(filas, punto, &mut distancias);
    distancias
}

fn llenar_distancias(filas: &[Vec<f64>], punto: usize, distancias: &mut Vec<(f64, usize)>) {
    distancias.clear();
    for (otro, fila) in filas.iter().enumerate() {
        if otro != punto {
            distancias.push((distancia_al_cuadrado(&filas[punto], fila), otro));
        }
    }
}

fn mas_cerca(a: &(f64, usize), b: &(f64, usize)) -> std::cmp::Ordering {
    a.0.total_cmp(&b.0).then(a.1.cmp(&b.1))
}

pub fn quedarse_con_los_k(distancias: &mut Vec<(f64, usize)>, k: usize) -> Vec<usize> {
    if k < distancias.len() {
        distancias.select_nth_unstable_by(k - 1, mas_cerca);
        distancias.truncate(k);
    }
    distancias.sort_unstable_by(mas_cerca);
    distancias.iter().map(|&(_, otro)| otro).collect()
}

pub fn vecinos_mas_cercanos(filas: &[Vec<f64>], k: usize) -> Vec<Vec<usize>> {
    (0..filas.len())
        .into_par_iter()
        .map_init(
            || Vec::with_capacity(filas.len()),
            |distancias, punto| {
                llenar_distancias(filas, punto, distancias);
                quedarse_con_los_k(distancias, k)
            },
        )
        .collect()
}

pub fn vecindario_con_el_mismo(vecinos: &[usize], punto: usize) -> Vec<usize> {
    let mut vecindario = vecinos.to_vec();
    vecindario.push(punto);
    vecindario.sort_unstable();
    vecindario
}

pub fn compartidos(a: &[usize], b: &[usize]) -> usize {
    let (mut i, mut j, mut cuenta) = (0, 0, 0);
    while i < a.len() && j < b.len() {
        if a[i] == b[j] {
            cuenta += 1;
            i += 1;
            j += 1;
        } else if a[i] < b[j] {
            i += 1;
        } else {
            j += 1;
        }
    }
    cuenta
}

pub fn peso_jaccard(a: &[usize], b: &[usize]) -> f64 {
    let comunes = compartidos(a, b);
    comunes as f64 / (a.len() + b.len() - comunes) as f64
}

pub fn grafo_snn(vecinos: &[Vec<usize>]) -> Grafo {
    let vecindarios: Vec<Vec<usize>> = vecinos
        .iter()
        .enumerate()
        .map(|(punto, lista)| vecindario_con_el_mismo(lista, punto))
        .collect();
    let aristas: Vec<Vec<(usize, usize, f64)>> = (0..vecinos.len())
        .into_par_iter()
        .map(|punto| {
            let mut propias = Vec::new();
            for &otro in &vecinos[punto] {
                let ya_la_tiene_el_otro =
                    otro < punto && vecindarios[otro].binary_search(&punto).is_ok();
                if ya_la_tiene_el_otro {
                    continue;
                }
                let peso = peso_jaccard(&vecindarios[punto], &vecindarios[otro]);
                if peso >= PESO_MINIMO_SNN {
                    propias.push((punto, otro, peso));
                }
            }
            propias
        })
        .collect();
    let mut grafo = grafo_vacio(vecinos.len());
    for (a, b, peso) in aristas.into_iter().flatten() {
        agregar_arista(&mut grafo, a, b, peso);
    }
    for lista in &mut grafo.vecinos {
        lista.sort_unstable_by_key(|&(otro, _)| otro);
    }
    grafo
}

pub fn grafo_vacio(nodos: usize) -> Grafo {
    Grafo {
        vecinos: vec![Vec::new(); nodos],
        lazos: vec![0.0; nodos],
        grados: vec![0.0; nodos],
        peso_total: 0.0,
    }
}

pub fn agregar_arista(grafo: &mut Grafo, a: usize, b: usize, peso: f64) {
    if a == b {
        grafo.lazos[a] += peso;
        grafo.grados[a] += 2.0 * peso;
    } else {
        grafo.vecinos[a].push((b, peso));
        grafo.vecinos[b].push((a, peso));
        grafo.grados[a] += peso;
        grafo.grados[b] += peso;
    }
    grafo.peso_total += peso;
}

pub fn modularidad(grafo: &Grafo, comunidad: &[usize], resolucion: f64) -> f64 {
    let m = grafo.peso_total;
    if m <= 0.0 {
        return 0.0;
    }
    let cuantas = comunidad.iter().max().map_or(0, |mayor| mayor + 1);
    let mut adentro = vec![0.0; cuantas];
    let mut grado = vec![0.0; cuantas];
    for nodo in 0..grafo.vecinos.len() {
        let propia = comunidad[nodo];
        grado[propia] += grafo.grados[nodo];
        adentro[propia] += grafo.lazos[nodo];
        for &(otro, peso) in &grafo.vecinos[nodo] {
            if otro > nodo && comunidad[otro] == propia {
                adentro[propia] += peso;
            }
        }
    }
    let mut total = 0.0;
    for c in 0..cuantas {
        let fraccion = grado[c] / (2.0 * m);
        total += adentro[c] / m - resolucion * fraccion * fraccion;
    }
    total
}

pub fn particion_inicial(grafo: &Grafo) -> Particion {
    Particion {
        comunidad: (0..grafo.vecinos.len()).collect(),
        grado_de_comunidad: grafo.grados.clone(),
    }
}

pub fn mover_un_nodo(
    grafo: &Grafo,
    particion: &mut Particion,
    nodo: usize,
    resolucion: f64,
    peso_hacia: &mut [f64],
    tocadas: &mut Vec<usize>,
) -> Movimiento {
    let m = grafo.peso_total;
    let grado = grafo.grados[nodo];
    let desde = particion.comunidad[nodo];
    tocadas.clear();
    for &(otro, peso) in &grafo.vecinos[nodo] {
        let suya = particion.comunidad[otro];
        if peso_hacia[suya] == 0.0 {
            tocadas.push(suya);
        }
        peso_hacia[suya] += peso;
    }
    particion.grado_de_comunidad[desde] -= grado;
    let valor = |comunidad: usize, grados: &[f64], pesos: &[f64]| {
        pesos[comunidad] - resolucion * grados[comunidad] * grado / (2.0 * m)
    };
    let quedarse = valor(desde, &particion.grado_de_comunidad, peso_hacia);
    let mut hacia = desde;
    let mut mejor = quedarse;
    for &candidata in tocadas.iter() {
        let posible = valor(candidata, &particion.grado_de_comunidad, peso_hacia);
        if posible > mejor + MEJORA_MINIMA {
            mejor = posible;
            hacia = candidata;
        }
    }
    particion.grado_de_comunidad[hacia] += grado;
    particion.comunidad[nodo] = hacia;
    for &comunidad in tocadas.iter() {
        peso_hacia[comunidad] = 0.0;
    }
    Movimiento {
        nodo,
        desde,
        hacia,
        ganancia: (mejor - quedarse) / m,
    }
}

pub fn pasada_local(
    grafo: &Grafo,
    particion: &mut Particion,
    orden: &[usize],
    resolucion: f64,
) -> usize {
    let mut peso_hacia = vec![0.0; grafo.vecinos.len()];
    let mut tocadas = Vec::new();
    let mut movidos = 0;
    for &nodo in orden {
        let movimiento = mover_un_nodo(
            grafo,
            particion,
            nodo,
            resolucion,
            &mut peso_hacia,
            &mut tocadas,
        );
        if movimiento.hacia != movimiento.desde {
            movidos += 1;
        }
    }
    movidos
}

pub fn barajar(cantidad: usize, azar: &mut Azar) -> Vec<usize> {
    let mut orden: Vec<usize> = (0..cantidad).collect();
    for posicion in (1..cantidad).rev() {
        let otra = azar.entero(posicion + 1);
        orden.swap(posicion, otra);
    }
    orden
}

pub fn renumerar(comunidad: &[usize]) -> (Vec<usize>, usize) {
    let etiquetas = comunidad.iter().max().map_or(0, |mayor| mayor + 1);
    let mut nuevo_numero = vec![usize::MAX; etiquetas];
    let mut cantidad = 0;
    let mut renumerada = Vec::with_capacity(comunidad.len());
    for &vieja in comunidad {
        if nuevo_numero[vieja] == usize::MAX {
            nuevo_numero[vieja] = cantidad;
            cantidad += 1;
        }
        renumerada.push(nuevo_numero[vieja]);
    }
    (renumerada, cantidad)
}

pub fn agregar_comunidades(grafo: &Grafo, comunidad: &[usize], cantidad: usize) -> Grafo {
    let mut nuevo = grafo_vacio(cantidad);
    let mut acumulado: Vec<BTreeMap<usize, f64>> = vec![BTreeMap::new(); cantidad];
    for nodo in 0..grafo.vecinos.len() {
        let propia = comunidad[nodo];
        nuevo.lazos[propia] += grafo.lazos[nodo];
        for &(otro, peso) in &grafo.vecinos[nodo] {
            if otro < nodo {
                continue;
            }
            let suya = comunidad[otro];
            if suya == propia {
                nuevo.lazos[propia] += peso;
            } else {
                *acumulado[propia.min(suya)]
                    .entry(propia.max(suya))
                    .or_insert(0.0) += peso;
            }
        }
    }
    let lazos = std::mem::take(&mut nuevo.lazos);
    nuevo.lazos = vec![0.0; cantidad];
    for (c, lazo) in lazos.into_iter().enumerate() {
        if lazo > 0.0 {
            agregar_arista(&mut nuevo, c, c, lazo);
        }
    }
    for (a, destinos) in acumulado.into_iter().enumerate() {
        for (b, peso) in destinos {
            agregar_arista(&mut nuevo, a, b, peso);
        }
    }
    for lista in &mut nuevo.vecinos {
        lista.sort_unstable_by_key(|&(otro, _)| otro);
    }
    nuevo
}

pub fn louvain(grafo: &Grafo, resolucion: f64, semilla: u64) -> Result<Louvain, String> {
    if !resolucion.is_finite() || resolucion <= 0.0 {
        return Err("La resolución debe ser un número mayor que cero.".to_string());
    }
    let mut azar = Azar::nuevo(semilla);
    let mut asignacion: Vec<usize> = (0..grafo.vecinos.len()).collect();
    let mut niveles = Vec::new();
    let mut agregado: Option<Grafo> = None;
    if grafo.peso_total > 0.0 {
        loop {
            let actual = agregado.as_ref().unwrap_or(grafo);
            let mut particion = particion_inicial(actual);
            let orden = barajar(actual.vecinos.len(), &mut azar);
            let mut movimientos = 0;
            for _ in 0..MAXIMO_PASADAS {
                let movidos = pasada_local(actual, &mut particion, &orden, resolucion);
                movimientos += movidos;
                if movidos == 0 {
                    break;
                }
            }
            if movimientos == 0 {
                break;
            }
            let (comunidad, cantidad) = renumerar(&particion.comunidad);
            for nodo in asignacion.iter_mut() {
                *nodo = comunidad[*nodo];
            }
            niveles.push(Nivel {
                nodos: actual.vecinos.len(),
                movimientos,
                comunidades: cantidad,
                modularidad: modularidad(grafo, &asignacion, resolucion),
            });
            let siguiente = agregar_comunidades(actual, &comunidad, cantidad);
            agregado = Some(siguiente);
        }
    }
    let (comunidad, cantidad) = ordenar_por_tamano(&asignacion);
    Ok(Louvain {
        modularidad: modularidad(grafo, &comunidad, resolucion),
        comunidad,
        cantidad,
        niveles,
    })
}

pub fn ordenar_por_tamano(comunidad: &[usize]) -> (Vec<usize>, usize) {
    let (renumerada, cantidad) = renumerar(comunidad);
    let mut tamanos = vec![0usize; cantidad];
    for &c in &renumerada {
        tamanos[c] += 1;
    }
    let mut orden: Vec<usize> = (0..cantidad).collect();
    orden.sort_by_key(|&c| (std::cmp::Reverse(tamanos[c]), c));
    let mut lugar = vec![0; cantidad];
    for (posicion, &c) in orden.iter().enumerate() {
        lugar[c] = posicion;
    }
    (renumerada.iter().map(|&c| lugar[c]).collect(), cantidad)
}

fn revisar_vecinos(puntos: &[Vec<f64>], k: usize) -> Result<(), String> {
    revisar_puntos(puntos)?;
    if puntos.len() < 3 {
        return Err("Hacen falta al menos 3 filas para formar un grafo de vecinos.".to_string());
    }
    if k == 0 || k >= puntos.len() {
        return Err(format!(
            "La cantidad de vecinos k debe estar entre 1 y {} (las filas menos una).",
            puntos.len() - 1
        ));
    }
    Ok(())
}

pub fn agrupar_por_comunidades(
    puntos: &[Vec<f64>],
    medida: Medida,
    k: usize,
    resolucion: f64,
    semilla: u64,
) -> Result<Comunidades, String> {
    revisar_vecinos(puntos, k)?;
    let filas = filas_para_medir(puntos, medida);
    let vecinos = vecinos_mas_cercanos(&filas, k);
    let grafo = grafo_snn(&vecinos);
    let louvain = louvain(&grafo, resolucion, semilla)?;
    Ok(Comunidades {
        vecinos,
        grafo,
        louvain,
    })
}

pub fn indice_rand_ajustado(a: &[usize], b: &[usize]) -> f64 {
    let n = a.len().min(b.len());
    if n < 2 {
        return 1.0;
    }
    let filas = a[..n].iter().max().map_or(0, |mayor| mayor + 1);
    let columnas = b[..n].iter().max().map_or(0, |mayor| mayor + 1);
    let mut tabla = vec![0usize; filas * columnas];
    let mut por_fila = vec![0usize; filas];
    let mut por_columna = vec![0usize; columnas];
    for (&x, &y) in a.iter().zip(b).take(n) {
        tabla[x * columnas + y] += 1;
        por_fila[x] += 1;
        por_columna[y] += 1;
    }
    let parejas = |cuenta: usize| (cuenta * cuenta.saturating_sub(1) / 2) as f64;
    let indice: f64 = tabla.iter().map(|&cuenta| parejas(cuenta)).sum();
    let suma_filas: f64 = por_fila.iter().map(|&cuenta| parejas(cuenta)).sum();
    let suma_columnas: f64 = por_columna.iter().map(|&cuenta| parejas(cuenta)).sum();
    let esperado = suma_filas * suma_columnas / parejas(n);
    let maximo = (suma_filas + suma_columnas) / 2.0;
    if (maximo - esperado).abs() < PARTE_DESPRECIABLE {
        return 1.0;
    }
    (indice - esperado) / (maximo - esperado)
}

pub fn entropia_con_beta(
    distancias: &[f64],
    propio: usize,
    beta: f64,
    probabilidades: &mut [f64],
) -> f64 {
    let (minima, maxima) = distancias
        .iter()
        .enumerate()
        .filter(|&(otro, _)| otro != propio)
        .fold(
            (f64::INFINITY, 0.0f64),
            |(menor, mayor), (_, &distancia)| (menor.min(distancia), mayor.max(distancia)),
        );
    let ruido = maxima * PARTE_DESPRECIABLE;
    let sobre_la_minima = |distancia: f64| {
        let diferencia = distancia - minima;
        if diferencia <= ruido { 0.0 } else { diferencia }
    };
    let mut suma = 0.0;
    for (otro, &distancia) in distancias.iter().enumerate() {
        probabilidades[otro] = if otro == propio {
            0.0
        } else {
            (-sobre_la_minima(distancia) * beta).exp()
        };
        suma += probabilidades[otro];
    }
    let mut ponderada = 0.0;
    for (probabilidad, &distancia) in probabilidades.iter_mut().zip(distancias) {
        *probabilidad /= suma;
        ponderada += *probabilidad * sobre_la_minima(distancia);
    }
    suma.ln() + beta * ponderada
}

pub fn beta_inicial(distancias: &[f64], propio: usize) -> f64 {
    let otros = distancias.len().saturating_sub(1);
    let suma: f64 = distancias
        .iter()
        .enumerate()
        .filter(|&(otro, _)| otro != propio)
        .map(|(_, &distancia)| distancia)
        .sum();
    let media = suma / otros as f64;
    if media.is_finite() && media > 0.0 {
        1.0 / media
    } else {
        1.0
    }
}

pub fn calibrar_un_punto(
    distancias: &[f64],
    propio: usize,
    perplejidad: f64,
    probabilidades: &mut [f64],
) -> Calibracion {
    let objetivo = perplejidad.ln();
    let mut beta = beta_inicial(distancias, propio);
    let mut beta_minima = 0.0;
    let mut beta_maxima = f64::INFINITY;
    let mut entropia = entropia_con_beta(distancias, propio, beta, probabilidades);
    let mut intentos = 1;
    while (entropia - objetivo).abs() > TOLERANCIA_DE_ENTROPIA && intentos < INTENTOS_DE_BETA {
        if entropia > objetivo {
            beta_minima = beta;
            beta = if beta_maxima.is_finite() {
                (beta + beta_maxima) / 2.0
            } else {
                beta * 2.0
            };
        } else {
            beta_maxima = beta;
            beta = (beta + beta_minima) / 2.0;
        }
        entropia = entropia_con_beta(distancias, propio, beta, probabilidades);
        intentos += 1;
    }
    Calibracion {
        beta,
        entropia,
        intentos,
    }
}

pub fn afinidades(filas: &[Vec<f64>], perplejidad: f64) -> (Vec<f64>, Vec<f64>) {
    let n = filas.len();
    let mut condicionales = vec![0.0; n * n];
    let sigmas: Vec<f64> = condicionales
        .par_chunks_mut(n)
        .enumerate()
        .map(|(punto, probabilidades)| {
            let distancias: Vec<f64> = filas
                .iter()
                .map(|otra| distancia_al_cuadrado(&filas[punto], otra))
                .collect();
            let calibracion = calibrar_un_punto(&distancias, punto, perplejidad, probabilidades);
            (1.0 / (2.0 * calibracion.beta)).sqrt()
        })
        .collect();
    (simetrizar(condicionales, n), sigmas)
}

pub fn simetrizar(mut condicionales: Vec<f64>, n: usize) -> Vec<f64> {
    let total = 2.0 * n as f64;
    for i in 0..n {
        condicionales[i * n + i] = 0.0;
        for j in i + 1..n {
            let conjunta = (condicionales[i * n + j] + condicionales[j * n + i]) / total;
            let conjunta = conjunta.max(PROBABILIDAD_MINIMA);
            condicionales[i * n + j] = conjunta;
            condicionales[j * n + i] = conjunta;
        }
    }
    condicionales
}

pub fn coordenadas_iniciales(n: usize, semilla: u64) -> Vec<[f64; 2]> {
    let mut azar = Azar::nuevo(semilla);
    (0..n)
        .map(|_| {
            [
                azar.normal() * ESCALA_INICIAL,
                azar.normal() * ESCALA_INICIAL,
            ]
        })
        .collect()
}

pub fn gradiente(conjuntas: &[f64], coordenadas: &[[f64; 2]], exageracion: f64) -> Vec<[f64; 2]> {
    let n = coordenadas.len();
    let parciales: Vec<([f64; 2], [f64; 2], f64)> = (0..n)
        .into_par_iter()
        .map(|i| {
            let fila = &conjuntas[i * n..(i + 1) * n];
            let [xi, yi] = coordenadas[i];
            let mut atraccion = [0.0, 0.0];
            let mut repulsion = [0.0, 0.0];
            let mut suma = 0.0;
            for (j, &[xj, yj]) in coordenadas.iter().enumerate() {
                if j == i {
                    continue;
                }
                let (dx, dy) = (xi - xj, yi - yj);
                let cercania = 1.0 / (1.0 + dx * dx + dy * dy);
                suma += cercania;
                let tira = fila[j] * cercania;
                atraccion[0] += tira * dx;
                atraccion[1] += tira * dy;
                let empuja = cercania * cercania;
                repulsion[0] += empuja * dx;
                repulsion[1] += empuja * dy;
            }
            (atraccion, repulsion, suma)
        })
        .collect();
    let normalizador: f64 = parciales.iter().map(|&(_, _, suma)| suma).sum();
    parciales
        .into_iter()
        .map(|(atraccion, repulsion, _)| {
            [
                4.0 * (exageracion * atraccion[0] - repulsion[0] / normalizador),
                4.0 * (exageracion * atraccion[1] - repulsion[1] / normalizador),
            ]
        })
        .collect()
}

pub fn divergencia_kl(conjuntas: &[f64], coordenadas: &[[f64; 2]]) -> f64 {
    let n = coordenadas.len();
    let parciales: Vec<(f64, f64, f64)> = (0..n)
        .into_par_iter()
        .map(|i| {
            let fila = &conjuntas[i * n..(i + 1) * n];
            let [xi, yi] = coordenadas[i];
            let mut suma = 0.0;
            let mut masa = 0.0;
            let mut contribucion = 0.0;
            for (j, &[xj, yj]) in coordenadas.iter().enumerate() {
                if j == i {
                    continue;
                }
                let (dx, dy) = (xi - xj, yi - yj);
                let cercania = 1.0 / (1.0 + dx * dx + dy * dy);
                suma += cercania;
                masa += fila[j];
                contribucion += fila[j] * (fila[j] / cercania).ln();
            }
            (contribucion, masa, suma)
        })
        .collect();
    let normalizador: f64 = parciales.iter().map(|&(_, _, suma)| suma).sum();
    let masa: f64 = parciales.iter().map(|&(_, masa, _)| masa).sum();
    let contribucion: f64 = parciales.iter().map(|&(parte, _, _)| parte).sum();
    contribucion + masa * normalizador.ln()
}

pub fn tasa_automatica(n: usize, exageracion: f64) -> f64 {
    (n as f64 / exageracion / 4.0).max(50.0)
}

pub fn descender(
    coordenadas: &mut [[f64; 2]],
    velocidad: &mut [[f64; 2]],
    ganancias: &mut [[f64; 2]],
    gradiente: &[[f64; 2]],
    momento: f64,
    tasa: f64,
) {
    for i in 0..coordenadas.len() {
        for eje in 0..2 {
            let sigue_bajando = velocidad[i][eje] * gradiente[i][eje] < 0.0;
            ganancias[i][eje] = if sigue_bajando {
                ganancias[i][eje] + 0.2
            } else {
                ganancias[i][eje] * 0.8
            };
            ganancias[i][eje] = ganancias[i][eje].max(GANANCIA_MINIMA);
            velocidad[i][eje] =
                momento * velocidad[i][eje] - tasa * ganancias[i][eje] * gradiente[i][eje];
            coordenadas[i][eje] += velocidad[i][eje];
        }
    }
    centrar(coordenadas);
}

fn centrar(coordenadas: &mut [[f64; 2]]) {
    let n = coordenadas.len() as f64;
    let mut media = [0.0, 0.0];
    for punto in coordenadas.iter() {
        media[0] += punto[0] / n;
        media[1] += punto[1] / n;
    }
    for punto in coordenadas.iter_mut() {
        punto[0] -= media[0];
        punto[1] -= media[1];
    }
}

fn revisar_tsne(puntos: &[Vec<f64>], opciones: &OpcionesTsne) -> Result<(), String> {
    revisar_puntos(puntos)?;
    if puntos.len() < 4 {
        return Err("t-SNE necesita al menos 4 filas.".to_string());
    }
    let tope = (puntos.len() - 1) as f64 / 3.0;
    if !opciones.perplejidad.is_finite()
        || opciones.perplejidad < 1.0
        || opciones.perplejidad > tope
    {
        return Err(format!(
            "La perplejidad debe estar entre 1 y {tope:.1}: con {} filas, tres veces la perplejidad no puede pasar de las filas menos una.",
            puntos.len()
        ));
    }
    if !opciones.exageracion.is_finite() || opciones.exageracion < 1.0 {
        return Err("La exageración temprana debe ser un número mayor o igual a 1.".to_string());
    }
    if opciones.iteraciones == 0 {
        return Err("Hace falta al menos una iteración.".to_string());
    }
    Ok(())
}

pub fn tsne(puntos: &[Vec<f64>], medida: Medida, opciones: OpcionesTsne) -> Result<Tsne, String> {
    revisar_tsne(puntos, &opciones)?;
    let n = puntos.len();
    let filas = filas_para_medir(puntos, medida);
    let (conjuntas, sigmas) = afinidades(&filas, opciones.perplejidad);
    let tasa = tasa_automatica(n, opciones.exageracion);
    let mut coordenadas = coordenadas_iniciales(n, opciones.semilla);
    let mut velocidad = vec![[0.0, 0.0]; n];
    let mut ganancias = vec![[1.0, 1.0]; n];
    let mut historia = Vec::new();
    for iteracion in 1..=opciones.iteraciones {
        let exagerada = iteracion <= ITERACIONES_EXAGERADAS;
        let exageracion = if exagerada { opciones.exageracion } else { 1.0 };
        let momento = if exagerada {
            MOMENTO_INICIAL
        } else {
            MOMENTO_FINAL
        };
        if iteracion == ITERACIONES_EXAGERADAS + 1 {
            velocidad.fill([0.0, 0.0]);
            ganancias.fill([1.0, 1.0]);
        }
        let empuje = gradiente(&conjuntas, &coordenadas, exageracion);
        descender(
            &mut coordenadas,
            &mut velocidad,
            &mut ganancias,
            &empuje,
            momento,
            tasa,
        );
        if iteracion % CADA_CUANTO_MEDIR == 0 || iteracion == opciones.iteraciones {
            historia.push((iteracion, divergencia_kl(&conjuntas, &coordenadas)));
        }
    }
    let divergencia = historia.last().map_or(0.0, |&(_, valor)| valor);
    Ok(Tsne {
        coordenadas,
        divergencia,
        historia,
        sigmas,
        tasa,
    })
}

pub fn conservacion_de_vecinos(originales: &[Vec<f64>], coordenadas: &[[f64; 2]], k: usize) -> f64 {
    let n = originales.len().min(coordenadas.len());
    if n < 2 || k == 0 {
        return 0.0;
    }
    let k = k.min(n - 1);
    let plano: Vec<Vec<f64>> = coordenadas[..n]
        .iter()
        .map(|punto| punto.to_vec())
        .collect();
    let antes = vecinos_mas_cercanos(&originales[..n], k);
    let despues = vecinos_mas_cercanos(&plano, k);
    let mut conservados = 0;
    for (lista_antes, lista_despues) in antes.iter().zip(&despues) {
        let mut ordenada = lista_despues.clone();
        ordenada.sort_unstable();
        conservados += lista_antes
            .iter()
            .filter(|vecino| ordenada.binary_search(vecino).is_ok())
            .count();
    }
    conservados as f64 / (n * k) as f64
}
