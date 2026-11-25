use crate::azar::Azar;
use crate::filogenia::{self, Arbol};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

pub const MAXIMO_DE_SECUENCIAS: usize = 60;
pub const MAXIMO_DE_POSICIONES: usize = 5000;
pub const MAXIMO_DE_REPLICAS: usize = 200;
pub const CATEGORIAS_GAMMA: usize = 4;
pub const FALTANTE: u8 = 4;
const RAMA_MINIMA: f64 = 1e-8;
const RAMA_MAXIMA: f64 = 10.0;
const RAMA_INICIAL_MINIMA: f64 = 1e-3;
const ESCALA: f64 = 1.157_920_892_373_162e77;
const UMBRAL_DE_ESCALA: f64 = 1.0 / ESCALA;
const TOLERANCIA: f64 = 1e-2;
const TOLERANCIA_FINAL: f64 = 1e-4;
const MEJORA_MINIMA: f64 = 1e-3;
const MAXIMO_DE_RONDAS: usize = 200;
pub const MAXIMO_DE_PASADAS: usize = 20;
const PASADAS_POR_RONDA: usize = 2;
const MAXIMO_DE_CICLOS: usize = 20;
const MAXIMO_DE_PASOS_DE_NEWTON: usize = 30;
const MAXIMO_DE_ITERACIONES_BRENT: usize = 60;
const KAPPA_MINIMO: f64 = 0.05;
const KAPPA_MAXIMO: f64 = 100.0;
const ALFA_MINIMO: f64 = 0.02;
const ALFA_MAXIMO: f64 = 100.0;
const FRECUENCIA_MINIMA: f64 = 1e-3;
const DISTANCIA_SATURADA: f64 = 5.0;
const RAMA_SIN_APOYO: f64 = 1e-6;
const PATRONES_POR_TAREA: usize = 256;
const PROPORCION_DORADA: f64 = 0.381_966_011_250_105;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Sustitucion {
    Jc69,
    K80,
    Hky85,
}

impl Sustitucion {
    pub fn nombre(self) -> &'static str {
        match self {
            Sustitucion::Jc69 => "JC69",
            Sustitucion::K80 => "K80",
            Sustitucion::Hky85 => "HKY85",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Modelo {
    pub sustitucion: Sustitucion,
    pub kappa: f64,
    pub frecuencias: [f64; 4],
    pub alfa: Option<f64>,
}

impl Modelo {
    pub fn nuevo(sustitucion: Sustitucion, con_gamma: bool, observadas: [f64; 4]) -> Modelo {
        let frecuencias = match sustitucion {
            Sustitucion::Hky85 => observadas,
            _ => [0.25; 4],
        };
        let kappa = match sustitucion {
            Sustitucion::Jc69 => 1.0,
            _ => 2.0,
        };
        Modelo {
            sustitucion,
            kappa,
            frecuencias,
            alfa: con_gamma.then_some(0.5),
        }
    }

    pub fn nombre(&self) -> String {
        match self.alfa {
            Some(_) => format!("{}+Γ", self.sustitucion.nombre()),
            None => self.sustitucion.nombre().to_string(),
        }
    }

    pub fn parametros_libres(&self) -> usize {
        let de_la_sustitucion = match self.sustitucion {
            Sustitucion::Jc69 => 0,
            Sustitucion::K80 => 1,
            Sustitucion::Hky85 => 4,
        };
        de_la_sustitucion + usize::from(self.alfa.is_some())
    }

    pub fn tasas(&self) -> Vec<f64> {
        match self.alfa {
            Some(alfa) => tasas_gamma(alfa, CATEGORIAS_GAMMA),
            None => vec![1.0],
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Espectro {
    pub valores: [f64; 4],
    pub matrices: [[[f64; 4]; 4]; 4],
}

fn es_purina(base: usize) -> bool {
    base == 0 || base == 2
}

pub fn espectro(modelo: &Modelo) -> Espectro {
    let pi = modelo.frecuencias;
    let kappa = modelo.kappa;
    let purinas = pi[0] + pi[2];
    let pirimidinas = pi[1] + pi[3];
    let beta = 1.0 / (2.0 * (purinas * pirimidinas + kappa * (pi[0] * pi[2] + pi[1] * pi[3])));
    let valores = [
        0.0,
        -beta,
        -beta * (1.0 + purinas * (kappa - 1.0)),
        -beta * (1.0 + pirimidinas * (kappa - 1.0)),
    ];
    let parte = |termino: usize, i: usize, j: usize| {
        let grupo_de_j = if es_purina(j) { purinas } else { pirimidinas };
        let mismo_grupo = es_purina(i) == es_purina(j);
        let identidad = if i == j { 1.0 } else { 0.0 };
        match termino {
            0 => pi[j],
            1 if mismo_grupo => pi[j] * (1.0 / grupo_de_j - 1.0),
            1 => -pi[j],
            2 if mismo_grupo && es_purina(i) => identidad - pi[j] / purinas,
            3 if mismo_grupo && !es_purina(i) => identidad - pi[j] / pirimidinas,
            _ => 0.0,
        }
    };
    let matrices = std::array::from_fn(|termino| {
        std::array::from_fn(|i| std::array::from_fn(|j| parte(termino, i, j)))
    });
    Espectro { valores, matrices }
}

pub fn matriz_de_transicion(espectro: &Espectro, tiempo: f64) -> [[f64; 4]; 4] {
    let mut matriz = [[0.0; 4]; 4];
    for (valor, termino) in espectro.valores.iter().zip(&espectro.matrices) {
        let peso = (valor * tiempo).exp();
        for (fila, fila_del_termino) in matriz.iter_mut().zip(termino) {
            for (celda, parte) in fila.iter_mut().zip(fila_del_termino) {
                *celda += peso * parte;
            }
        }
    }
    matriz
}

fn ln_gamma(x: f64) -> f64 {
    let mut desplazado = x;
    let mut correccion = 0.0;
    while desplazado < 15.0 {
        correccion -= desplazado.ln();
        desplazado += 1.0;
    }
    let inverso = 1.0 / desplazado;
    let cuadrado = inverso * inverso;
    let serie = inverso
        * (1.0 / 12.0 - cuadrado * (1.0 / 360.0 - cuadrado * (1.0 / 1260.0 - cuadrado / 1680.0)));
    correccion + (desplazado - 0.5) * desplazado.ln() - desplazado
        + 0.5 * (2.0 * std::f64::consts::PI).ln()
        + serie
}

pub fn gamma_incompleta(a: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x.is_infinite() {
        return 1.0;
    }
    let frente = (a * x.ln() - x - ln_gamma(a)).exp();
    if x < a + 1.0 {
        let mut termino = 1.0 / a;
        let mut suma = termino;
        let mut denominador = a;
        for _ in 0..10_000 {
            denominador += 1.0;
            termino *= x / denominador;
            suma += termino;
            if termino < suma * 1e-16 {
                break;
            }
        }
        return (frente * suma).min(1.0);
    }
    let diminuto = 1e-300;
    let mut b = x + 1.0 - a;
    let mut c = 1.0 / diminuto;
    let mut d = 1.0 / b;
    let mut fraccion = d;
    for i in 1..10_000 {
        let numerador = -(i as f64) * (i as f64 - a);
        b += 2.0;
        d = numerador * d + b;
        if d.abs() < diminuto {
            d = diminuto;
        }
        c = b + numerador / c;
        if c.abs() < diminuto {
            c = diminuto;
        }
        d = 1.0 / d;
        let cambio = d * c;
        fraccion *= cambio;
        if (cambio - 1.0).abs() < 1e-16 {
            break;
        }
    }
    (1.0 - frente * fraccion).max(0.0)
}

pub fn cuantil_gamma(a: f64, probabilidad: f64) -> f64 {
    let mut alto = a.max(1.0);
    while gamma_incompleta(a, alto) < probabilidad {
        alto *= 2.0;
    }
    let mut log_bajo = -690.0;
    let mut log_alto = alto.ln();
    for _ in 0..200 {
        let log_medio = 0.5 * (log_bajo + log_alto);
        if gamma_incompleta(a, log_medio.exp()) < probabilidad {
            log_bajo = log_medio;
        } else {
            log_alto = log_medio;
        }
    }
    (0.5 * (log_bajo + log_alto)).exp()
}

pub fn tasas_gamma(alfa: f64, categorias: usize) -> Vec<f64> {
    let mut tasas = Vec::with_capacity(categorias);
    let mut acumulado_anterior = 0.0;
    for k in 1..=categorias {
        let acumulado = if k == categorias {
            1.0
        } else {
            let limite = cuantil_gamma(alfa, k as f64 / categorias as f64);
            gamma_incompleta(alfa + 1.0, limite)
        };
        tasas.push((acumulado - acumulado_anterior) * categorias as f64);
        acumulado_anterior = acumulado;
    }
    tasas
}

#[derive(Debug, Clone, PartialEq)]
pub struct Datos {
    pub nombres: Vec<String>,
    pub estados: Vec<Vec<u8>>,
    pub pesos: Vec<f64>,
    pub patron_de_sitio: Vec<usize>,
}

pub fn estado_de_letra(letra: u8) -> u8 {
    match letra.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => FALTANTE,
    }
}

pub fn comprimir_en_patrones(secuencias: &[Vec<u8>], nombres: &[String]) -> Datos {
    let hojas = secuencias.len();
    let largo = secuencias.first().map_or(0, Vec::len);
    let mut indice_de_columna: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut estados = vec![Vec::new(); hojas];
    let mut pesos: Vec<f64> = Vec::new();
    let mut patron_de_sitio = Vec::with_capacity(largo);
    let mut columna = vec![FALTANTE; hojas];
    for sitio in 0..largo {
        for (estado, secuencia) in columna.iter_mut().zip(secuencias) {
            *estado = estado_de_letra(secuencia[sitio]);
        }
        let siguiente = pesos.len();
        let patron = *indice_de_columna
            .entry(columna.clone())
            .or_insert(siguiente);
        if patron == siguiente {
            pesos.push(0.0);
            for (estados_de_la_hoja, &estado) in estados.iter_mut().zip(&columna) {
                estados_de_la_hoja.push(estado);
            }
        }
        pesos[patron] += 1.0;
        patron_de_sitio.push(patron);
    }
    Datos {
        nombres: nombres.to_vec(),
        estados,
        pesos,
        patron_de_sitio,
    }
}

pub fn frecuencias_observadas(datos: &Datos) -> [f64; 4] {
    let mut cuentas = [0.0; 4];
    for estados in &datos.estados {
        for (&estado, peso) in estados.iter().zip(&datos.pesos) {
            if estado < FALTANTE {
                cuentas[estado as usize] += peso;
            }
        }
    }
    let total: f64 = cuentas.iter().sum();
    if total <= 0.0 {
        return [0.25; 4];
    }
    let acotadas = cuentas.map(|cuenta| (cuenta / total).max(FRECUENCIA_MINIMA));
    let suma: f64 = acotadas.iter().sum();
    acotadas.map(|frecuencia| frecuencia / suma)
}

pub fn datos_de_un_sitio(datos: &Datos, sitio: usize) -> Datos {
    let patron = datos.patron_de_sitio[sitio];
    Datos {
        nombres: datos.nombres.clone(),
        estados: datos
            .estados
            .iter()
            .map(|estados| vec![estados[patron]])
            .collect(),
        pesos: vec![1.0],
        patron_de_sitio: vec![0],
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Parcial {
    pub valores: Vec<[f64; 4]>,
    pub escalas: Vec<u32>,
}

pub fn vector_de_estado(estado: u8) -> [f64; 4] {
    if estado >= FALTANTE {
        return [1.0; 4];
    }
    let mut vector = [0.0; 4];
    vector[estado as usize] = 1.0;
    vector
}

fn parcial_de_hoja(estados: &[u8], categorias: usize) -> Parcial {
    let mut valores = Vec::with_capacity(estados.len() * categorias);
    for &estado in estados {
        let vector = vector_de_estado(estado);
        for _ in 0..categorias {
            valores.push(vector);
        }
    }
    Parcial {
        valores,
        escalas: vec![0; estados.len()],
    }
}

fn aplicar(matriz: &[[f64; 4]; 4], vector: &[f64; 4]) -> [f64; 4] {
    let mut resultado = [0.0; 4];
    for (salida, fila) in resultado.iter_mut().zip(matriz) {
        *salida =
            fila[0] * vector[0] + fila[1] * vector[1] + fila[2] * vector[2] + fila[3] * vector[3];
    }
    resultado
}

pub fn multiplicar(primero: &Parcial, segundo: &Parcial) -> Parcial {
    let categorias = primero.valores.len() / primero.escalas.len().max(1);
    let mut producto = Parcial {
        valores: vec![[0.0; 4]; primero.valores.len()],
        escalas: vec![0; primero.escalas.len()],
    };
    producto
        .valores
        .par_chunks_mut(categorias.max(1))
        .zip(producto.escalas.par_iter_mut())
        .enumerate()
        .with_min_len(PATRONES_POR_TAREA)
        .for_each(|(patron, (bloque, escala))| {
            let desde = patron * categorias;
            let de_primero = &primero.valores[desde..desde + categorias];
            let de_segundo = &segundo.valores[desde..desde + categorias];
            for ((vector, a), b) in bloque.iter_mut().zip(de_primero).zip(de_segundo) {
                *vector = [a[0] * b[0], a[1] * b[1], a[2] * b[2], a[3] * b[3]];
            }
            *escala = primero.escalas[patron] + segundo.escalas[patron];
            reescalar(bloque, escala);
        });
    producto
}

fn reescalar(bloque: &mut [[f64; 4]], escala: &mut u32) {
    let mayor = bloque
        .iter()
        .flat_map(|vector| vector.iter())
        .fold(0.0f64, |mayor, &valor| mayor.max(valor));
    if mayor > 0.0 && mayor < UMBRAL_DE_ESCALA {
        for vector in bloque.iter_mut() {
            for valor in vector.iter_mut() {
                *valor *= ESCALA;
            }
        }
        *escala += 1;
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Derivadas {
    pub valor: f64,
    pub primera: f64,
    pub segunda: f64,
}

pub fn hijos(arbol: &Arbol, nodo: usize, padre: usize) -> Vec<usize> {
    arbol.vecinos[nodo]
        .iter()
        .map(|&(vecino, _)| vecino)
        .filter(|&vecino| vecino != padre)
        .collect()
}

pub fn longitud(arbol: &Arbol, a: usize, b: usize) -> f64 {
    arbol.longitud_entre(a, b).unwrap_or(RAMA_MINIMA)
}

pub fn fijar_longitud(arbol: &mut Arbol, a: usize, b: usize, nueva: f64) {
    for (vecino, largo) in arbol.vecinos[a].iter_mut() {
        if *vecino == b {
            *largo = nueva;
        }
    }
    for (vecino, largo) in arbol.vecinos[b].iter_mut() {
        if *vecino == a {
            *largo = nueva;
        }
    }
}

pub fn orden_posterior(arbol: &Arbol) -> Vec<(usize, usize)> {
    let mut orden = Vec::with_capacity(arbol.nombres.len());
    let Some(&(primero, _)) = arbol.vecinos.first().and_then(|vecinos| vecinos.first()) else {
        return orden;
    };
    let mut pila = vec![(primero, 0)];
    while let Some((nodo, padre)) = pila.pop() {
        orden.push((nodo, padre));
        for hijo in hijos(arbol, nodo, padre) {
            pila.push((hijo, nodo));
        }
    }
    orden.reverse();
    orden
}

pub struct Calculo<'a> {
    datos: &'a Datos,
    pub modelo: Modelo,
    pub tasas: Vec<f64>,
    pub espectro: Espectro,
    pub abajo: Vec<Parcial>,
    log_escala: f64,
    pub tolerancia: f64,
}

impl<'a> Calculo<'a> {
    pub fn nuevo(datos: &'a Datos, modelo: Modelo, nodos: usize) -> Calculo<'a> {
        let tasas = modelo.tasas();
        let mut abajo = vec![Parcial::default(); nodos];
        for (hoja, estados) in datos.estados.iter().enumerate() {
            abajo[hoja] = parcial_de_hoja(estados, tasas.len());
        }
        Calculo {
            datos,
            modelo,
            tasas,
            espectro: espectro(&modelo),
            abajo,
            log_escala: ESCALA.ln(),
            tolerancia: TOLERANCIA,
        }
    }

    pub fn cambiar_modelo(&mut self, modelo: Modelo) {
        self.modelo = modelo;
        self.espectro = espectro(&modelo);
        self.tasas = modelo.tasas();
    }

    pub fn matrices(&self, rama: f64) -> Vec<[[f64; 4]; 4]> {
        self.tasas
            .iter()
            .map(|tasa| matriz_de_transicion(&self.espectro, tasa * rama))
            .collect()
    }

    pub fn combinar_en(
        &self,
        destino: &mut Parcial,
        primero: &Parcial,
        rama_primero: f64,
        segundo: &Parcial,
        rama_segundo: f64,
    ) {
        let matrices_primero = self.matrices(rama_primero);
        let matrices_segundo = self.matrices(rama_segundo);
        let categorias = self.tasas.len();
        let patrones = self.datos.pesos.len();
        destino.valores.resize(patrones * categorias, [0.0; 4]);
        destino.escalas.resize(patrones, 0);
        destino
            .valores
            .par_chunks_mut(categorias)
            .zip(destino.escalas.par_iter_mut())
            .enumerate()
            .with_min_len(PATRONES_POR_TAREA)
            .for_each(|(patron, (bloque, escala))| {
                let desde = patron * categorias;
                let de_primero = &primero.valores[desde..desde + categorias];
                let de_segundo = &segundo.valores[desde..desde + categorias];
                let matrices = matrices_primero.iter().zip(&matrices_segundo);
                for (((vector, a), b), (matriz_a, matriz_b)) in bloque
                    .iter_mut()
                    .zip(de_primero)
                    .zip(de_segundo)
                    .zip(matrices)
                {
                    let desde_a = aplicar(matriz_a, a);
                    let desde_b = aplicar(matriz_b, b);
                    *vector = [
                        desde_a[0] * desde_b[0],
                        desde_a[1] * desde_b[1],
                        desde_a[2] * desde_b[2],
                        desde_a[3] * desde_b[3],
                    ];
                }
                *escala = primero.escalas[patron] + segundo.escalas[patron];
                reescalar(bloque, escala);
            });
    }

    pub fn combinar(
        &self,
        primero: &Parcial,
        rama_primero: f64,
        segundo: &Parcial,
        rama_segundo: f64,
    ) -> Parcial {
        let mut destino = Parcial::default();
        self.combinar_en(&mut destino, primero, rama_primero, segundo, rama_segundo);
        destino
    }

    pub fn propagar(&self, parcial: &Parcial, rama: f64) -> Parcial {
        let matrices = self.matrices(rama);
        let categorias = self.tasas.len();
        let valores = parcial
            .valores
            .iter()
            .enumerate()
            .map(|(indice, vector)| aplicar(&matrices[indice % categorias], vector))
            .collect();
        Parcial {
            valores,
            escalas: parcial.escalas.clone(),
        }
    }

    pub fn recalcular_nodo(&mut self, arbol: &Arbol, nodo: usize, padre: usize) {
        let hijos = hijos(arbol, nodo, padre);
        if hijos.len() != 2 {
            return;
        }
        let mut destino = std::mem::take(&mut self.abajo[nodo]);
        self.combinar_en(
            &mut destino,
            &self.abajo[hijos[0]],
            longitud(arbol, nodo, hijos[0]),
            &self.abajo[hijos[1]],
            longitud(arbol, nodo, hijos[1]),
        );
        self.abajo[nodo] = destino;
    }

    pub fn recalcular_abajo(&mut self, arbol: &Arbol) {
        for (nodo, padre) in orden_posterior(arbol) {
            self.recalcular_nodo(arbol, nodo, padre);
        }
    }

    pub fn coeficientes(&self, lado_a: &Parcial, lado_b: &Parcial) -> Vec<[f64; 4]> {
        let pi = self.modelo.frecuencias;
        let purinas = pi[0] + pi[2];
        let pirimidinas = pi[1] + pi[3];
        lado_a
            .valores
            .par_iter()
            .zip(&lado_b.valores)
            .with_min_len(PATRONES_POR_TAREA)
            .map(|(a, b)| {
                let purinas_a = pi[0] * a[0] + pi[2] * a[2];
                let pirimidinas_a = pi[1] * a[1] + pi[3] * a[3];
                let purinas_b = pi[0] * b[0] + pi[2] * b[2];
                let pirimidinas_b = pi[1] * b[1] + pi[3] * b[3];
                let iguales_purinas = pi[0] * a[0] * b[0] + pi[2] * a[2] * b[2];
                let iguales_pirimidinas = pi[1] * a[1] * b[1] + pi[3] * a[3] * b[3];
                [
                    (purinas_a + pirimidinas_a) * (purinas_b + pirimidinas_b),
                    purinas_a * ((1.0 / purinas - 1.0) * purinas_b - pirimidinas_b)
                        + pirimidinas_a * ((1.0 / pirimidinas - 1.0) * pirimidinas_b - purinas_b),
                    iguales_purinas - purinas_a * purinas_b / purinas,
                    iguales_pirimidinas - pirimidinas_a * pirimidinas_b / pirimidinas,
                ]
            })
            .collect()
    }

    pub fn correccion_de_escala(&self, lado_a: &Parcial, lado_b: &Parcial) -> f64 {
        let mut correccion = 0.0;
        for ((peso, escala_a), escala_b) in self
            .datos
            .pesos
            .iter()
            .zip(&lado_a.escalas)
            .zip(&lado_b.escalas)
        {
            correccion += peso * (escala_a + escala_b) as f64;
        }
        -correccion * self.log_escala
    }

    pub fn evaluar_rama(&self, coeficientes: &[[f64; 4]], rama: f64) -> Derivadas {
        let categorias = self.tasas.len();
        let peso_de_categoria = 1.0 / categorias as f64;
        let mut factores = vec![[0.0; 4]; categorias];
        let mut exponenciales = vec![[0.0; 4]; categorias];
        for (categoria, tasa) in self.tasas.iter().enumerate() {
            for termino in 0..4 {
                factores[categoria][termino] = self.espectro.valores[termino] * tasa;
                exponenciales[categoria][termino] = (factores[categoria][termino] * rama).exp();
            }
        }
        let mut derivadas = Derivadas {
            valor: 0.0,
            primera: 0.0,
            segunda: 0.0,
        };
        for (de_patron, peso) in coeficientes.chunks_exact(categorias).zip(&self.datos.pesos) {
            let mut sitio = 0.0;
            let mut primera = 0.0;
            let mut segunda = 0.0;
            for ((terminos, exponencial), factor) in
                de_patron.iter().zip(&exponenciales).zip(&factores)
            {
                for termino in 0..4 {
                    let parte = terminos[termino] * exponencial[termino];
                    sitio += parte;
                    primera += parte * factor[termino];
                    segunda += parte * factor[termino] * factor[termino];
                }
            }
            let sitio = (sitio * peso_de_categoria).max(f64::MIN_POSITIVE);
            let primera = primera * peso_de_categoria / sitio;
            let segunda = segunda * peso_de_categoria / sitio;
            derivadas.valor += peso * sitio.ln();
            derivadas.primera += peso * primera;
            derivadas.segunda += peso * (segunda - primera * primera);
        }
        derivadas
    }

    pub fn optimizar_longitud(&self, coeficientes: &[[f64; 4]], inicial: f64) -> (f64, f64) {
        let mut rama = inicial.clamp(RAMA_MINIMA, RAMA_MAXIMA);
        let mut actual = self.evaluar_rama(coeficientes, rama);
        for _ in 0..MAXIMO_DE_PASOS_DE_NEWTON {
            let mut propuesta = if actual.segunda < 0.0 {
                rama - actual.primera / actual.segunda
            } else if actual.primera > 0.0 {
                rama * 2.0
            } else {
                rama / 2.0
            };
            propuesta = propuesta.clamp(RAMA_MINIMA, RAMA_MAXIMA);
            let mut nueva = self.evaluar_rama(coeficientes, propuesta);
            let mut intentos = 0;
            while nueva.valor < actual.valor && intentos < 30 {
                propuesta = 0.5 * (rama + propuesta);
                nueva = self.evaluar_rama(coeficientes, propuesta);
                intentos += 1;
            }
            if nueva.valor < actual.valor {
                break;
            }
            let cambio = (propuesta - rama).abs();
            rama = propuesta;
            actual = nueva;
            if cambio < 1e-6 * (1.0 + rama) {
                break;
            }
        }
        (rama, actual.valor)
    }

    pub fn evaluar_configuracion(
        &self,
        lado_a: &Parcial,
        lado_b: &Parcial,
        inicial: f64,
    ) -> (f64, f64) {
        let coeficientes = self.coeficientes(lado_a, lado_b);
        let (rama, valor) = self.optimizar_longitud(&coeficientes, inicial);
        (rama, valor + self.correccion_de_escala(lado_a, lado_b))
    }

    pub fn log_verosimilitud_en_raiz(&self, arbol: &Arbol) -> f64 {
        let hijo = arbol.vecinos[0][0].0;
        let coeficientes = self.coeficientes(&self.abajo[0], &self.abajo[hijo]);
        let derivadas = self.evaluar_rama(&coeficientes, longitud(arbol, 0, hijo));
        derivadas.valor + self.correccion_de_escala(&self.abajo[0], &self.abajo[hijo])
    }

    pub fn log_verosimilitud(&mut self, arbol: &Arbol) -> f64 {
        self.recalcular_abajo(arbol);
        self.log_verosimilitud_en_raiz(arbol)
    }

    fn optimizar_desde(
        &mut self,
        arbol: &mut Arbol,
        padre: usize,
        nodo: usize,
        de_arriba: &Parcial,
    ) {
        let coeficientes = self.coeficientes(de_arriba, &self.abajo[nodo]);
        let (rama, _) = self.optimizar_longitud(&coeficientes, longitud(arbol, padre, nodo));
        fijar_longitud(arbol, padre, nodo, rama);
        let hijos = hijos(arbol, nodo, padre);
        if hijos.len() != 2 {
            return;
        }
        for (posicion, &hijo) in hijos.iter().enumerate() {
            let hermano = hijos[1 - posicion];
            let arriba_del_hijo = self.combinar(
                de_arriba,
                rama,
                &self.abajo[hermano],
                longitud(arbol, nodo, hermano),
            );
            self.optimizar_desde(arbol, nodo, hijo, &arriba_del_hijo);
        }
        self.recalcular_nodo(arbol, nodo, padre);
    }

    pub fn una_pasada_de_ramas(&mut self, arbol: &mut Arbol) -> f64 {
        let hijo = arbol.vecinos[0][0].0;
        let hoja_de_la_raiz = self.abajo[0].clone();
        self.optimizar_desde(arbol, 0, hijo, &hoja_de_la_raiz);
        self.log_verosimilitud_en_raiz(arbol)
    }

    pub fn ajustar_ramas(&mut self, arbol: &mut Arbol, pasadas: usize) -> f64 {
        let mut actual = self.log_verosimilitud(arbol);
        for _ in 0..pasadas {
            let nuevo = self.una_pasada_de_ramas(arbol);
            let mejora = nuevo - actual;
            actual = nuevo;
            if mejora < self.tolerancia {
                break;
            }
        }
        actual
    }

    pub fn log_verosimilitud_con(&mut self, arbol: &Arbol, modelo: Modelo) -> f64 {
        self.cambiar_modelo(modelo);
        self.log_verosimilitud(arbol)
    }

    fn optimizar_parametro(
        &mut self,
        arbol: &Arbol,
        limites: (f64, f64),
        leer: fn(&Modelo) -> f64,
        escribir: fn(Modelo, f64) -> Modelo,
    ) {
        let inicial = self.modelo;
        let actual = self.log_verosimilitud(arbol);
        let tolerancia = self.tolerancia;
        let (mejor, valor) = minimizar_brent(
            |logaritmo| -self.log_verosimilitud_con(arbol, escribir(inicial, logaritmo.exp())),
            (limites.0.ln(), limites.1.ln()),
            leer(&inicial).clamp(limites.0, limites.1).ln(),
            tolerancia,
        );
        if -valor > actual {
            self.cambiar_modelo(escribir(inicial, mejor.exp()));
        } else {
            self.cambiar_modelo(inicial);
        }
    }

    pub fn optimizar_kappa(&mut self, arbol: &Arbol) {
        self.optimizar_parametro(
            arbol,
            (KAPPA_MINIMO, KAPPA_MAXIMO),
            |modelo| modelo.kappa,
            |modelo, kappa| Modelo { kappa, ..modelo },
        );
    }

    pub fn optimizar_alfa(&mut self, arbol: &Arbol) {
        self.optimizar_parametro(
            arbol,
            (ALFA_MINIMO, ALFA_MAXIMO),
            |modelo| modelo.alfa.unwrap_or(1.0),
            |modelo, alfa| Modelo {
                alfa: Some(alfa),
                ..modelo
            },
        );
    }

    pub fn intercambios_que_mejoran(&mut self, arbol: &Arbol) -> Vec<Intercambio> {
        self.recalcular_abajo(arbol);
        let mut candidatos = Vec::new();
        let hijo = arbol.vecinos[0][0].0;
        let hoja_de_la_raiz = self.abajo[0].clone();
        self.revisar_desde(arbol, 0, hijo, &hoja_de_la_raiz, &mut candidatos);
        candidatos.sort_by(|a, b| b.ganancia.total_cmp(&a.ganancia));
        candidatos
    }

    fn revisar_desde(
        &self,
        arbol: &Arbol,
        padre: usize,
        nodo: usize,
        de_arriba: &Parcial,
        candidatos: &mut Vec<Intercambio>,
    ) {
        let hijos = hijos(arbol, nodo, padre);
        if hijos.len() != 2 {
            return;
        }
        let desde_arriba = self.propagar(de_arriba, longitud(arbol, padre, nodo));
        for (posicion, &hijo) in hijos.iter().enumerate() {
            let hermano = hijos[1 - posicion];
            let desde_hermano = self.propagar(&self.abajo[hermano], longitud(arbol, nodo, hermano));
            if let Some(mejor) = self.mejor_intercambio_en(
                arbol,
                (nodo, hijo),
                hermano,
                &desde_arriba,
                &desde_hermano,
            ) {
                candidatos.push(mejor);
            }
            let arriba_del_hijo = multiplicar(&desde_arriba, &desde_hermano);
            self.revisar_desde(arbol, nodo, hijo, &arriba_del_hijo, candidatos);
        }
    }

    fn mejor_intercambio_en(
        &self,
        arbol: &Arbol,
        (nodo, hijo): (usize, usize),
        hermano: usize,
        desde_arriba: &Parcial,
        desde_hermano: &Parcial,
    ) -> Option<Intercambio> {
        let nietos = hijos(arbol, hijo, nodo);
        if nietos.len() != 2 {
            return None;
        }
        let desde_nietos: Vec<Parcial> = nietos
            .iter()
            .map(|&nieto| self.propagar(&self.abajo[nieto], longitud(arbol, hijo, nieto)))
            .collect();
        let central = longitud(arbol, nodo, hijo);
        let (_, actual) = self.evaluar_configuracion(
            &multiplicar(desde_arriba, desde_hermano),
            &multiplicar(&desde_nietos[0], &desde_nietos[1]),
            central,
        );
        let mut mejor: Option<Intercambio> = None;
        for (posicion, &nieto) in nietos.iter().enumerate() {
            let lado_del_nodo = multiplicar(desde_arriba, &desde_nietos[posicion]);
            let lado_del_hijo = multiplicar(desde_hermano, &desde_nietos[1 - posicion]);
            let (rama, valor) = self.evaluar_configuracion(&lado_del_nodo, &lado_del_hijo, central);
            let ganancia = valor - actual;
            let supera = mejor.is_none_or(|otro| ganancia > otro.ganancia);
            if ganancia > MEJORA_MINIMA && supera {
                mejor = Some(Intercambio {
                    extremos: (nodo, hijo),
                    subarboles: (hermano, nieto),
                    ganancia,
                    rama_central: rama,
                });
            }
        }
        mejor
    }
}

fn minimizar_brent(
    mut funcion: impl FnMut(f64) -> f64,
    (desde, hasta): (f64, f64),
    inicial: f64,
    precision: f64,
) -> (f64, f64) {
    let (mut a, mut b) = (desde, hasta);
    let mut x = inicial.clamp(a, b);
    let (mut w, mut v) = (x, x);
    let mut fx = funcion(x);
    let (mut fw, mut fv) = (fx, fx);
    let mut paso: f64 = 0.0;
    let mut paso_anterior: f64 = 0.0;
    for _ in 0..MAXIMO_DE_ITERACIONES_BRENT {
        let medio = 0.5 * (a + b);
        let tolerancia = precision * (1.0 + x.abs());
        if (x - medio).abs() <= 2.0 * tolerancia - 0.5 * (b - a) {
            break;
        }
        let mut usar_dorado = true;
        if paso_anterior.abs() > tolerancia {
            let r = (x - w) * (fx - fv);
            let mut q = (x - v) * (fx - fw);
            let mut p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if q > 0.0 {
                p = -p;
            } else {
                q = -q;
            }
            if p.abs() < (0.5 * q * paso_anterior).abs() && p > q * (a - x) && p < q * (b - x) {
                paso_anterior = paso;
                paso = p / q;
                let u = x + paso;
                if u - a < 2.0 * tolerancia || b - u < 2.0 * tolerancia {
                    paso = if medio > x { tolerancia } else { -tolerancia };
                }
                usar_dorado = false;
            }
        }
        if usar_dorado {
            paso_anterior = if x >= medio { a - x } else { b - x };
            paso = PROPORCION_DORADA * paso_anterior;
        }
        let u = if paso.abs() >= tolerancia {
            x + paso
        } else if paso > 0.0 {
            x + tolerancia
        } else {
            x - tolerancia
        };
        let fu = funcion(u);
        if fu <= fx {
            if u >= x {
                a = x;
            } else {
                b = x;
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        } else {
            if u < x {
                a = u;
            } else {
                b = u;
            }
            if fu <= fw || w == x {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            } else if fu <= fv || v == x || v == w {
                v = u;
                fv = fu;
            }
        }
    }
    (x, fx)
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Intercambio {
    pub extremos: (usize, usize),
    pub subarboles: (usize, usize),
    pub ganancia: f64,
    pub rama_central: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ronda {
    pub numero: usize,
    pub log_verosimilitud: f64,
    pub candidatos: usize,
    pub aplicados: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Busqueda {
    pub arbol: Arbol,
    pub log_verosimilitud: f64,
    pub rondas: Vec<Ronda>,
}

pub fn sin_conflictos(candidatos: &[Intercambio]) -> Vec<Intercambio> {
    let mut usados = HashSet::new();
    let mut elegidos = Vec::new();
    for candidato in candidatos {
        let nodos = [
            candidato.extremos.0,
            candidato.extremos.1,
            candidato.subarboles.0,
            candidato.subarboles.1,
        ];
        if nodos.iter().any(|nodo| usados.contains(nodo)) {
            continue;
        }
        usados.extend(nodos);
        elegidos.push(*candidato);
    }
    elegidos
}

pub fn aplicar_intercambios(arbol: &mut Arbol, elegidos: &[Intercambio]) {
    for intercambio in elegidos {
        let (nodo, hijo) = intercambio.extremos;
        let (hermano, nieto) = intercambio.subarboles;
        filogenia::intercambiar(arbol, nodo, hermano, hijo, nieto);
        fijar_longitud(arbol, nodo, hijo, intercambio.rama_central);
    }
}

pub fn buscar_con_intercambios(datos: &Datos, inicial: &Arbol, modelo: Modelo) -> Busqueda {
    let mut calculo = Calculo::nuevo(datos, modelo, inicial.nombres.len());
    let mut arbol = inicial.clone();
    let mut actual = calculo.ajustar_ramas(&mut arbol, MAXIMO_DE_PASADAS);
    let mut rondas = vec![Ronda {
        numero: 0,
        log_verosimilitud: actual,
        candidatos: 0,
        aplicados: 0,
    }];
    for numero in 1..=MAXIMO_DE_RONDAS {
        let candidatos = calculo.intercambios_que_mejoran(&arbol);
        if candidatos.is_empty() {
            break;
        }
        let elegidos = sin_conflictos(&candidatos);
        let anterior = arbol.clone();
        aplicar_intercambios(&mut arbol, &elegidos);
        let mut nuevo = calculo.ajustar_ramas(&mut arbol, PASADAS_POR_RONDA);
        let mut aplicados = elegidos.len();
        if nuevo < actual + MEJORA_MINIMA && aplicados > 1 {
            arbol = anterior.clone();
            aplicar_intercambios(&mut arbol, &elegidos[..1]);
            nuevo = calculo.ajustar_ramas(&mut arbol, PASADAS_POR_RONDA);
            aplicados = 1;
        }
        if nuevo < actual + MEJORA_MINIMA {
            arbol = anterior;
            break;
        }
        actual = nuevo;
        rondas.push(Ronda {
            numero,
            log_verosimilitud: actual,
            candidatos: candidatos.len(),
            aplicados,
        });
    }
    Busqueda {
        arbol,
        log_verosimilitud: actual,
        rondas,
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Ajuste {
    pub arbol: Arbol,
    pub modelo: Modelo,
    pub log_verosimilitud: f64,
}

pub fn ajustar_modelo(datos: &Datos, inicial: &Arbol, modelo: Modelo) -> Ajuste {
    ajustar_con_tolerancia(datos, inicial, modelo, TOLERANCIA)
}

pub fn afinar_modelo(datos: &Datos, inicial: &Arbol, modelo: Modelo) -> Ajuste {
    ajustar_con_tolerancia(datos, inicial, modelo, TOLERANCIA_FINAL)
}

fn ajustar_con_tolerancia(
    datos: &Datos,
    inicial: &Arbol,
    modelo: Modelo,
    tolerancia: f64,
) -> Ajuste {
    let mut arbol = inicial.clone();
    let mut calculo = Calculo::nuevo(datos, modelo, arbol.nombres.len());
    calculo.tolerancia = tolerancia;
    let mut actual = calculo.ajustar_ramas(&mut arbol, MAXIMO_DE_PASADAS);
    for _ in 0..MAXIMO_DE_CICLOS {
        if modelo.parametros_libres() == 0 {
            break;
        }
        if modelo.sustitucion != Sustitucion::Jc69 {
            calculo.optimizar_kappa(&arbol);
        }
        if modelo.alfa.is_some() {
            calculo.optimizar_alfa(&arbol);
        }
        let nuevo = calculo.ajustar_ramas(&mut arbol, PASADAS_POR_RONDA);
        let mejora = nuevo - actual;
        actual = nuevo;
        if mejora < tolerancia {
            break;
        }
    }
    Ajuste {
        arbol,
        modelo: calculo.modelo,
        log_verosimilitud: actual,
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Comparacion {
    pub modelo: Modelo,
    pub log_verosimilitud: f64,
    pub parametros: usize,
    pub aic: f64,
    pub arbol: Arbol,
}

pub fn criterio_de_akaike(log_verosimilitud: f64, parametros: usize) -> f64 {
    2.0 * parametros as f64 - 2.0 * log_verosimilitud
}

pub const MODELOS: [(Sustitucion, bool); 6] = [
    (Sustitucion::Jc69, false),
    (Sustitucion::Jc69, true),
    (Sustitucion::K80, false),
    (Sustitucion::K80, true),
    (Sustitucion::Hky85, false),
    (Sustitucion::Hky85, true),
];

pub fn comparar_modelos(datos: &Datos, arbol: &Arbol) -> Vec<Comparacion> {
    let observadas = frecuencias_observadas(datos);
    let ramas = arbol.aristas().len();
    let base = ajustar_modelo(
        datos,
        arbol,
        Modelo::nuevo(Sustitucion::Jc69, false, observadas),
    );
    MODELOS
        .par_iter()
        .map(|&(sustitucion, con_gamma)| {
            let modelo = Modelo::nuevo(sustitucion, con_gamma, observadas);
            let ajuste = if modelo.parametros_libres() == 0 {
                base.clone()
            } else {
                ajustar_modelo(datos, &base.arbol, modelo)
            };
            let parametros = ajuste.modelo.parametros_libres() + ramas;
            Comparacion {
                modelo: ajuste.modelo,
                log_verosimilitud: ajuste.log_verosimilitud,
                parametros,
                aic: criterio_de_akaike(ajuste.log_verosimilitud, parametros),
                arbol: ajuste.arbol,
            }
        })
        .collect()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Eleccion {
    MejorPorAic,
    Fija(Sustitucion, bool),
}

pub fn elegir_modelo(comparacion: &[Comparacion], eleccion: Eleccion) -> usize {
    let mut elegido = 0;
    for (indice, opcion) in comparacion.iter().enumerate() {
        let gana = match eleccion {
            Eleccion::MejorPorAic => opcion.aic < comparacion[elegido].aic,
            Eleccion::Fija(sustitucion, con_gamma) => {
                opcion.modelo.sustitucion == sustitucion
                    && opcion.modelo.alfa.is_some() == con_gamma
            }
        };
        if gana {
            elegido = indice;
        }
    }
    elegido
}

pub fn revisar_alineamiento(secuencias: &[Vec<u8>], nombres: &[String]) -> Result<usize, String> {
    if secuencias.len() != nombres.len() {
        return Err("Debe haber un nombre por secuencia.".to_string());
    }
    if secuencias.len() < 3 {
        return Err("Hacen falta al menos tres secuencias para un árbol sin raíz.".to_string());
    }
    if secuencias.len() > MAXIMO_DE_SECUENCIAS {
        return Err(format!(
            "Hay {} secuencias; la búsqueda de máxima verosimilitud admite hasta {MAXIMO_DE_SECUENCIAS}.",
            secuencias.len()
        ));
    }
    let largo = secuencias[0].len();
    if largo == 0 {
        return Err("Las secuencias están vacías.".to_string());
    }
    if largo > MAXIMO_DE_POSICIONES {
        return Err(format!(
            "Las secuencias miden {largo} posiciones; la búsqueda de máxima verosimilitud admite hasta {MAXIMO_DE_POSICIONES}."
        ));
    }
    let mut vistos = HashSet::new();
    for (nombre, secuencia) in nombres.iter().zip(secuencias) {
        if !vistos.insert(nombre) {
            return Err(format!("El nombre «{nombre}» está repetido."));
        }
        if secuencia.len() != largo {
            return Err(format!(
                "Las secuencias deben estar alineadas: «{}» mide {largo} y «{nombre}» mide {}.",
                nombres[0],
                secuencia.len()
            ));
        }
        let valida = |letra: &u8| {
            matches!(
                letra.to_ascii_uppercase(),
                b'A' | b'C' | b'G' | b'T' | b'N' | b'-' | b'?'
            )
        };
        if let Some(letra) = secuencia.iter().find(|letra| !valida(letra)) {
            return Err(format!(
                "La secuencia de «{nombre}» tiene «{}»: aquí solo valen A, C, G, T y «-», «N» o «?» para lo que falta.",
                *letra as char
            ));
        }
    }
    Ok(largo)
}

pub fn acotar_ramas(arbol: &mut Arbol) {
    for vecinos in arbol.vecinos.iter_mut() {
        for (_, largo) in vecinos.iter_mut() {
            *largo = largo.clamp(RAMA_INICIAL_MINIMA, RAMA_MAXIMA);
        }
    }
}

pub fn arbol_de_union(datos: &Datos) -> Result<Arbol, String> {
    let distancias = distancias_jukes_cantor(datos);
    filogenia::union_de_vecinos(&distancias, &datos.nombres)
}

pub fn arbol_inicial(datos: &Datos) -> Result<Arbol, String> {
    let mut arbol = arbol_de_union(datos)?;
    acotar_ramas(&mut arbol);
    Ok(arbol)
}

#[derive(Debug, Clone, PartialEq)]
pub struct Resultado {
    pub arbol: Arbol,
    pub modelo: Modelo,
    pub log_verosimilitud: f64,
    pub comparacion: Vec<Comparacion>,
    pub elegido: usize,
    pub rondas: Vec<Ronda>,
    pub arbol_de_union: Arbol,
    pub log_verosimilitud_de_union: f64,
    pub robinson_foulds: usize,
    pub sitios: usize,
    pub patrones: usize,
}

pub fn maxima_verosimilitud(
    secuencias: &[Vec<u8>],
    nombres: &[String],
    eleccion: Eleccion,
) -> Result<Resultado, String> {
    let sitios = revisar_alineamiento(secuencias, nombres)?;
    let datos = comprimir_en_patrones(secuencias, nombres);
    let inicial = arbol_inicial(&datos)?;
    let comparacion = comparar_modelos(&datos, &inicial);
    let elegido = elegir_modelo(&comparacion, eleccion);
    let inicio = &comparacion[elegido];
    let busqueda = buscar_con_intercambios(&datos, &inicio.arbol, inicio.modelo);
    let ajuste = afinar_modelo(&datos, &busqueda.arbol, inicio.modelo);
    let robinson_foulds = distancia_robinson_foulds(&ajuste.arbol, &inicio.arbol, nombres.len());
    Ok(Resultado {
        arbol: ajuste.arbol,
        modelo: ajuste.modelo,
        log_verosimilitud: ajuste.log_verosimilitud,
        rondas: busqueda.rondas,
        arbol_de_union: inicio.arbol.clone(),
        log_verosimilitud_de_union: inicio.log_verosimilitud,
        comparacion,
        elegido,
        robinson_foulds,
        sitios,
        patrones: datos.pesos.len(),
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Particion {
    pub extremos: (usize, usize),
    pub hojas: u64,
}

pub fn particiones(arbol: &Arbol, hojas: usize) -> Vec<Particion> {
    let mut mascaras = vec![0u64; arbol.nombres.len()];
    let mut resultado = Vec::new();
    for (nodo, padre) in orden_posterior(arbol) {
        if nodo < hojas {
            mascaras[nodo] = 1u64 << nodo;
            continue;
        }
        mascaras[nodo] = hijos(arbol, nodo, padre)
            .iter()
            .fold(0, |mascara, &hijo| mascara | mascaras[hijo]);
        if padre >= hojas {
            resultado.push(Particion {
                extremos: (padre, nodo),
                hojas: mascaras[nodo],
            });
        }
    }
    resultado
}

pub fn particiones_resueltas(arbol: &Arbol, hojas: usize) -> Vec<Particion> {
    particiones(arbol, hojas)
        .into_iter()
        .filter(|particion| {
            longitud(arbol, particion.extremos.0, particion.extremos.1).abs() > RAMA_SIN_APOYO
        })
        .collect()
}

pub fn distancia_robinson_foulds(primero: &Arbol, segundo: &Arbol, hojas: usize) -> usize {
    let de_primero: HashSet<u64> = particiones(primero, hojas)
        .iter()
        .map(|p| p.hojas)
        .collect();
    let de_segundo: HashSet<u64> = particiones(segundo, hojas)
        .iter()
        .map(|p| p.hojas)
        .collect();
    de_primero.symmetric_difference(&de_segundo).count()
}

pub fn columnas_de_replica(sitios: usize, replica: usize, semilla: u64) -> Vec<usize> {
    let mezcla = (replica as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    let mut azar = Azar::nuevo(semilla ^ mezcla);
    (0..sitios).map(|_| azar.entero(sitios)).collect()
}

pub fn remuestrear(datos: &Datos, columnas: &[usize]) -> Datos {
    let mut cuentas = vec![0.0; datos.pesos.len()];
    for &columna in columnas {
        cuentas[datos.patron_de_sitio[columna]] += 1.0;
    }
    let mut nuevo_indice = vec![usize::MAX; cuentas.len()];
    let mut pesos = Vec::new();
    for (patron, &cuenta) in cuentas.iter().enumerate() {
        if cuenta > 0.0 {
            nuevo_indice[patron] = pesos.len();
            pesos.push(cuenta);
        }
    }
    let estados = datos
        .estados
        .iter()
        .map(|estados| {
            estados
                .iter()
                .zip(&cuentas)
                .filter(|(_, cuenta)| **cuenta > 0.0)
                .map(|(&estado, _)| estado)
                .collect()
        })
        .collect();
    Datos {
        nombres: datos.nombres.clone(),
        estados,
        pesos,
        patron_de_sitio: columnas
            .iter()
            .map(|&columna| nuevo_indice[datos.patron_de_sitio[columna]])
            .collect(),
    }
}

fn distancia_entre_hojas(datos: &Datos, primera: usize, segunda: usize) -> f64 {
    let mut comparadas = 0.0;
    let mut diferentes = 0.0;
    for ((&a, &b), peso) in datos.estados[primera]
        .iter()
        .zip(&datos.estados[segunda])
        .zip(&datos.pesos)
    {
        if a < FALTANTE && b < FALTANTE {
            comparadas += peso;
            if a != b {
                diferentes += peso;
            }
        }
    }
    if comparadas == 0.0 {
        return DISTANCIA_SATURADA;
    }
    filogenia::distancia_jukes_cantor(diferentes / comparadas).min(DISTANCIA_SATURADA)
}

pub fn distancias_jukes_cantor(datos: &Datos) -> Vec<Vec<f64>> {
    let hojas = datos.estados.len();
    let mut distancias = vec![vec![0.0; hojas]; hojas];
    let pares: Vec<(usize, usize)> = (0..hojas)
        .flat_map(|i| ((i + 1)..hojas).map(move |j| (i, j)))
        .collect();
    for (i, j) in pares {
        let distancia = distancia_entre_hojas(datos, i, j);
        distancias[i][j] = distancia;
        distancias[j][i] = distancia;
    }
    distancias
}

pub fn arbol_de_replica(
    datos: &Datos,
    replica: usize,
    semilla: u64,
    modelo: Option<Modelo>,
) -> Result<Arbol, String> {
    let columnas = columnas_de_replica(datos.patron_de_sitio.len(), replica, semilla);
    let muestra = remuestrear(datos, &columnas);
    let distancias = distancias_jukes_cantor(&muestra);
    let mut arbol = filogenia::union_de_vecinos(&distancias, &muestra.nombres)?;
    match modelo {
        Some(modelo) => {
            acotar_ramas(&mut arbol);
            Ok(buscar_con_intercambios(&muestra, &arbol, modelo).arbol)
        }
        None => Ok(arbol),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Referencia {
    MaximaVerosimilitud,
    UnionDeVecinos,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Apoyo {
    pub particion: Particion,
    pub veces: usize,
    pub porcentaje: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Bootstrap {
    pub arbol: Arbol,
    pub apoyos: Vec<Apoyo>,
    pub replicas: usize,
    pub particiones_por_replica: Vec<Vec<u64>>,
    pub modelo: Option<Modelo>,
    pub log_verosimilitud: Option<f64>,
}

pub fn contar_apoyos(referencia: &[Particion], particiones_por_replica: &[Vec<u64>]) -> Vec<Apoyo> {
    let replicas = particiones_por_replica.len().max(1);
    referencia
        .iter()
        .map(|&particion| {
            let veces = particiones_por_replica
                .iter()
                .filter(|de_la_replica| de_la_replica.contains(&particion.hojas))
                .count();
            Apoyo {
                particion,
                veces,
                porcentaje: 100.0 * veces as f64 / replicas as f64,
            }
        })
        .collect()
}

pub fn bootstrap(
    secuencias: &[Vec<u8>],
    nombres: &[String],
    replicas: usize,
    referencia: Referencia,
    semilla: u64,
) -> Result<Bootstrap, String> {
    revisar_alineamiento(secuencias, nombres)?;
    if nombres.len() < 4 {
        return Err("El bootstrap necesita al menos cuatro secuencias: con tres no hay ramas internas que apoyar.".to_string());
    }
    if replicas == 0 || replicas > MAXIMO_DE_REPLICAS {
        return Err(format!(
            "Las réplicas deben ir de 1 a {MAXIMO_DE_REPLICAS}."
        ));
    }
    let datos = comprimir_en_patrones(secuencias, nombres);
    let (arbol, modelo, log_verosimilitud) = match referencia {
        Referencia::MaximaVerosimilitud => {
            let resultado = maxima_verosimilitud(secuencias, nombres, Eleccion::MejorPorAic)?;
            (
                resultado.arbol,
                Some(resultado.modelo),
                Some(resultado.log_verosimilitud),
            )
        }
        Referencia::UnionDeVecinos => (arbol_de_union(&datos)?, None, None),
    };
    let hojas = nombres.len();
    let particiones_por_replica = (0..replicas)
        .into_par_iter()
        .map(|replica| {
            let arbol = arbol_de_replica(&datos, replica, semilla, modelo)?;
            Ok(particiones_resueltas(&arbol, hojas)
                .iter()
                .map(|p| p.hojas)
                .collect())
        })
        .collect::<Result<Vec<Vec<u64>>, String>>()?;
    let apoyos = contar_apoyos(&particiones(&arbol, hojas), &particiones_por_replica);
    Ok(Bootstrap {
        arbol,
        apoyos,
        replicas,
        particiones_por_replica,
        modelo,
        log_verosimilitud,
    })
}

pub fn newick_con_apoyos(arbol: &Arbol, hojas: usize, apoyos: &[Apoyo]) -> String {
    let mut etiquetado = arbol.clone();
    for nombre in etiquetado.nombres.iter_mut().skip(hojas) {
        nombre.clear();
    }
    for apoyo in apoyos {
        etiquetado.nombres[apoyo.particion.extremos.1] = format!("{:.0}", apoyo.porcentaje);
    }
    let raiz = arbol
        .vecinos
        .first()
        .and_then(|vecinos| vecinos.first())
        .map_or(0, |&(vecino, _)| vecino);
    filogenia::escribir_newick(&etiquetado, raiz)
}
