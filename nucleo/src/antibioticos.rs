use crate::secuencia::{aminoacido, complemento_inverso};
use std::collections::BTreeMap;

pub const MASAS_ENTERAS: [(u8, usize); 20] = [
    (b'G', 57),
    (b'A', 71),
    (b'S', 87),
    (b'P', 97),
    (b'V', 99),
    (b'T', 101),
    (b'C', 103),
    (b'I', 113),
    (b'L', 113),
    (b'N', 114),
    (b'D', 115),
    (b'K', 128),
    (b'Q', 128),
    (b'E', 129),
    (b'M', 131),
    (b'H', 137),
    (b'F', 147),
    (b'R', 156),
    (b'Y', 163),
    (b'W', 186),
];

pub const MASAS_DISTINTAS: [usize; 18] = [
    57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186,
];

pub const MASA_CARBONO_12: f64 = 12.0;
pub const MASA_HIDROGENO_1: f64 = 1.00782503207;
pub const MASA_NITROGENO_14: f64 = 14.0030740048;
pub const MASA_OXIGENO_16: f64 = 15.99491461956;
pub const MASA_AZUFRE_32: f64 = 31.97207100;

pub const MASA_MINIMA_DE_AMINOACIDO: usize = 57;
pub const MASA_MAXIMA_DE_AMINOACIDO: usize = 200;
pub const MASA_MAXIMA: usize = 100_000;
pub const MASA_MAXIMA_PARA_CONTAR: usize = 3000;
pub const MAXIMO_MASAS_EXACTAS: usize = 200;
pub const MAXIMO_CANDIDATOS: usize = 1_000_000;
pub const MAXIMO_CANDIDATOS_EN_TOTAL: usize = 2_000_000;
pub const MAXIMO_PUNTOS: usize = 40;
pub const MAXIMO_INTENTOS: usize = 1_000_000;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Formula {
    pub carbono: u32,
    pub hidrogeno: u32,
    pub nitrogeno: u32,
    pub oxigeno: u32,
    pub azufre: u32,
}

pub const AGUA: Formula = Formula {
    carbono: 0,
    hidrogeno: 2,
    nitrogeno: 0,
    oxigeno: 1,
    azufre: 0,
};

#[derive(Debug, Clone, PartialEq)]
pub struct Traduccion {
    pub peptido: Vec<u8>,
    pub parada: Option<usize>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Codificacion {
    pub inicio: usize,
    pub hebra: u8,
    pub fragmento: Vec<u8>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Lideres {
    pub alfabeto: Vec<usize>,
    pub puntaje: usize,
    pub peptidos: Vec<Vec<usize>>,
}

pub fn traducir_arn(arn: &[u8]) -> Result<Traduccion, String> {
    let mut peptido = Vec::with_capacity(arn.len() / 3);
    for (numero, codon) in arn.as_chunks::<3>().0.iter().enumerate() {
        let letra = aminoacido(codon);
        if letra == b'X' {
            return Err(format!(
                "El codón «{}» (posición {}, contando desde 0) tiene letras que no son A, C, G ni U.",
                String::from_utf8_lossy(codon),
                numero * 3
            ));
        }
        if letra == b'*' {
            return Ok(Traduccion {
                peptido,
                parada: Some(numero * 3),
            });
        }
        peptido.push(letra);
    }
    Ok(Traduccion {
        peptido,
        parada: None,
    })
}

pub fn masa_entera(letra_buscada: u8) -> Option<usize> {
    for (letra, masa) in MASAS_ENTERAS {
        if letra == letra_buscada {
            return Some(masa);
        }
    }
    None
}

pub fn masas_del_peptido(peptido: &[u8]) -> Result<Vec<usize>, String> {
    let mut masas = Vec::with_capacity(peptido.len());
    for &letra in peptido {
        match masa_entera(letra) {
            Some(masa) => masas.push(masa),
            None => return Err(no_es_aminoacido(letra)),
        }
    }
    Ok(masas)
}

fn no_es_aminoacido(letra: u8) -> String {
    format!(
        "«{}» no es uno de los 20 aminoácidos estándar.",
        letra as char
    )
}

fn revisar_peptido(peptido: &[u8]) -> Result<(), String> {
    if peptido.is_empty() {
        return Err("El péptido no tiene ningún aminoácido.".to_string());
    }
    masas_del_peptido(peptido)?;
    Ok(())
}

pub fn codificaciones_del_peptido(
    texto: &[u8],
    peptido: &[u8],
) -> Result<Vec<Codificacion>, String> {
    revisar_peptido(peptido)?;
    let largo = 3 * peptido.len();
    let mut encontradas = Vec::new();
    if largo > texto.len() {
        return Ok(encontradas);
    }
    let directa = aminoacido_en_cada_posicion(texto);
    let inversa = aminoacido_en_cada_posicion(&complemento_inverso(texto));
    let ultimo_inicio = texto.len() - largo;
    for inicio in 0..=ultimo_inicio {
        if codifica_desde(&directa, inicio, peptido) {
            encontradas.push(Codificacion {
                inicio,
                hebra: b'+',
                fragmento: texto[inicio..inicio + largo].to_vec(),
            });
        }
    }
    for inicio_en_la_inversa in 0..=ultimo_inicio {
        if codifica_desde(&inversa, inicio_en_la_inversa, peptido) {
            let inicio = ultimo_inicio - inicio_en_la_inversa;
            encontradas.push(Codificacion {
                inicio,
                hebra: b'-',
                fragmento: texto[inicio..inicio + largo].to_vec(),
            });
        }
    }
    encontradas.sort_by_key(|codificacion| codificacion.inicio);
    Ok(encontradas)
}

fn aminoacido_en_cada_posicion(adn: &[u8]) -> Vec<u8> {
    adn.windows(3).map(aminoacido).collect()
}

fn codifica_desde(aminoacidos: &[u8], inicio: usize, peptido: &[u8]) -> bool {
    for (indice, &letra) in peptido.iter().enumerate() {
        if aminoacidos[inicio + 3 * indice] != letra {
            return false;
        }
    }
    true
}

pub fn codones_de(letra_buscada: u8) -> usize {
    let bases = b"ACGT";
    let mut cuantos = 0;
    for &primera in bases {
        for &segunda in bases {
            for &tercera in bases {
                if aminoacido(&[primera, segunda, tercera]) == letra_buscada {
                    cuantos += 1;
                }
            }
        }
    }
    cuantos
}

pub fn secuencias_que_codifican(peptido: &[u8]) -> f64 {
    let mut total = 1.0;
    for &letra in peptido {
        total *= codones_de(letra) as f64;
    }
    total
}

pub fn formula_del_peptido(peptido: &[u8], ciclico: bool) -> Result<Formula, String> {
    let mut total = if ciclico {
        Formula {
            carbono: 0,
            hidrogeno: 0,
            nitrogeno: 0,
            oxigeno: 0,
            azufre: 0,
        }
    } else {
        AGUA
    };
    for &letra in peptido {
        let residuo = formula_de_residuo(letra).ok_or_else(|| no_es_aminoacido(letra))?;
        total.carbono += residuo.carbono;
        total.hidrogeno += residuo.hidrogeno;
        total.nitrogeno += residuo.nitrogeno;
        total.oxigeno += residuo.oxigeno;
        total.azufre += residuo.azufre;
    }
    Ok(total)
}

fn formula_de_residuo(letra: u8) -> Option<Formula> {
    let (carbono, hidrogeno, nitrogeno, oxigeno, azufre) = match letra {
        b'G' => (2, 3, 1, 1, 0),
        b'A' => (3, 5, 1, 1, 0),
        b'S' => (3, 5, 1, 2, 0),
        b'P' => (5, 7, 1, 1, 0),
        b'V' => (5, 9, 1, 1, 0),
        b'T' => (4, 7, 1, 2, 0),
        b'C' => (3, 5, 1, 1, 1),
        b'I' | b'L' => (6, 11, 1, 1, 0),
        b'N' => (4, 6, 2, 2, 0),
        b'D' => (4, 5, 1, 3, 0),
        b'K' => (6, 12, 2, 1, 0),
        b'Q' => (5, 8, 2, 2, 0),
        b'E' => (5, 7, 1, 3, 0),
        b'M' => (5, 9, 1, 1, 1),
        b'H' => (6, 7, 3, 1, 0),
        b'F' => (9, 9, 1, 1, 0),
        b'R' => (6, 12, 4, 1, 0),
        b'Y' => (9, 9, 1, 2, 0),
        b'W' => (11, 10, 2, 1, 0),
        _ => return None,
    };
    Some(Formula {
        carbono,
        hidrogeno,
        nitrogeno,
        oxigeno,
        azufre,
    })
}

pub fn masa_monoisotopica(peptido: &[u8], ciclico: bool) -> Result<f64, String> {
    let atomos = formula_del_peptido(peptido, ciclico)?;
    Ok(atomos.carbono as f64 * MASA_CARBONO_12
        + atomos.hidrogeno as f64 * MASA_HIDROGENO_1
        + atomos.nitrogeno as f64 * MASA_NITROGENO_14
        + atomos.oxigeno as f64 * MASA_OXIGENO_16
        + atomos.azufre as f64 * MASA_AZUFRE_32)
}

pub fn masa_total(peptido: &[usize]) -> usize {
    peptido.iter().sum()
}

pub fn masas_de_subpeptidos(peptido: &[usize], ciclico: bool, masas: &mut Vec<usize>) {
    masas.clear();
    masas.push(0);
    let total = masa_total(peptido);
    let n = peptido.len();
    for inicio in 0..n {
        let mut suma = 0;
        for (fin, &masa) in peptido.iter().enumerate().skip(inicio) {
            suma += masa;
            masas.push(suma);
            if ciclico && inicio > 0 && fin + 1 < n {
                masas.push(total - suma);
            }
        }
    }
}

pub fn espectro_lineal(peptido: &[usize]) -> Vec<usize> {
    let mut masas = Vec::new();
    masas_de_subpeptidos(peptido, false, &mut masas);
    masas.sort_unstable();
    masas
}

pub fn espectro_ciclico(peptido: &[usize]) -> Vec<usize> {
    let mut masas = Vec::new();
    masas_de_subpeptidos(peptido, true, &mut masas);
    masas.sort_unstable();
    masas
}

pub fn peptidos_por_masa(masa: usize) -> Result<Vec<u128>, String> {
    if masa > MASA_MAXIMA_PARA_CONTAR {
        return Err(format!(
            "La masa debe ser como mucho {MASA_MAXIMA_PARA_CONTAR} Da: un poco más arriba el número de péptidos ya no cabe en un entero de 128 bits."
        ));
    }
    let mut cuantos = vec![0u128; masa + 1];
    cuantos[0] = 1;
    for total in 1..=masa {
        for ultima in MASAS_DISTINTAS {
            if ultima <= total {
                cuantos[total] += cuantos[total - ultima];
            }
        }
    }
    Ok(cuantos)
}

pub struct Espectro {
    pub masas: Vec<usize>,
    pub masa_madre: usize,
    conteo: Vec<u32>,
    usados: Vec<u32>,
    teorico: Vec<usize>,
}

impl Espectro {
    pub fn nuevo(masas: &[usize]) -> Result<Espectro, String> {
        if masas.is_empty() {
            return Err("El espectro no tiene ninguna masa.".to_string());
        }
        let mut ordenadas = masas.to_vec();
        ordenadas.sort_unstable();
        let masa_madre = ordenadas[ordenadas.len() - 1];
        if masa_madre > MASA_MAXIMA {
            return Err(format!(
                "La masa {masa_madre} es demasiado grande: el límite es {MASA_MAXIMA} Da."
            ));
        }
        let mut conteo = vec![0; masa_madre + 1];
        for &masa in &ordenadas {
            conteo[masa] += 1;
        }
        Ok(Espectro {
            masas: ordenadas,
            masa_madre,
            conteo,
            usados: vec![0; masa_madre + 1],
            teorico: Vec::new(),
        })
    }

    pub fn veces(&self, masa: usize) -> usize {
        if masa < self.conteo.len() {
            self.conteo[masa] as usize
        } else {
            0
        }
    }

    pub fn compartidas(&mut self, peptido: &[usize], ciclico: bool) -> usize {
        masas_de_subpeptidos(peptido, ciclico, &mut self.teorico);
        let mut compartidas = 0;
        for &masa in &self.teorico {
            if masa < self.conteo.len() && self.usados[masa] < self.conteo[masa] {
                self.usados[masa] += 1;
                compartidas += 1;
            }
        }
        for &masa in &self.teorico {
            if masa < self.usados.len() {
                self.usados[masa] = 0;
            }
        }
        compartidas
    }

    pub fn es_consistente(&mut self, peptido: &[usize]) -> bool {
        self.compartidas(peptido, false) == self.teorico.len()
    }

    pub fn explica_exactamente(&mut self, peptido: &[usize]) -> bool {
        let compartidas = self.compartidas(peptido, true);
        compartidas == self.teorico.len() && compartidas == self.masas.len()
    }
}

pub fn secuenciar_ciclopeptido(espectro: &[usize]) -> Result<Vec<Vec<usize>>, String> {
    let mut experimental = espectro_para_ramificar(espectro)?;
    let bloques = aminoacidos_presentes(&experimental);
    let mut candidatos: Vec<Vec<usize>> = vec![Vec::new()];
    let mut encontrados: Vec<Vec<usize>> = Vec::new();
    while !candidatos.is_empty() {
        candidatos = expandir(&candidatos, &bloques)?;
        candidatos = conservar_consistentes(candidatos, &mut experimental, &mut encontrados);
    }
    Ok(encontrados)
}

pub fn espectro_para_ramificar(espectro: &[usize]) -> Result<Espectro, String> {
    if espectro.len() > MAXIMO_MASAS_EXACTAS {
        return Err(format!(
            "La ramificación y poda acepta espectros de hasta {MAXIMO_MASAS_EXACTAS} masas (un péptido cíclico de 14 aminoácidos). Para espectros más grandes o con ruido usa la tabla de líderes."
        ));
    }
    Espectro::nuevo(espectro)
}

pub fn aminoacidos_presentes(experimental: &Espectro) -> Vec<usize> {
    let mut presentes = Vec::new();
    for masa in MASAS_DISTINTAS {
        if experimental.veces(masa) > 0 {
            presentes.push(masa);
        }
    }
    presentes
}

pub fn expandir(peptidos: &[Vec<usize>], bloques: &[usize]) -> Result<Vec<Vec<usize>>, String> {
    let cantidad = peptidos.len().saturating_mul(bloques.len());
    if cantidad > MAXIMO_CANDIDATOS {
        return Err(
            "La búsqueda pasaría de un millón de candidatos, el límite: el espectro es demasiado grande o demasiado ambiguo."
                .to_string(),
        );
    }
    let mut expandidos = Vec::with_capacity(cantidad);
    for peptido in peptidos {
        for &masa in bloques {
            let mut alargado = Vec::with_capacity(peptido.len() + 1);
            alargado.extend_from_slice(peptido);
            alargado.push(masa);
            expandidos.push(alargado);
        }
    }
    Ok(expandidos)
}

pub fn conservar_consistentes(
    candidatos: Vec<Vec<usize>>,
    experimental: &mut Espectro,
    encontrados: &mut Vec<Vec<usize>>,
) -> Vec<Vec<usize>> {
    let mut conservados = Vec::new();
    for peptido in candidatos {
        if masa_total(&peptido) == experimental.masa_madre {
            if experimental.explica_exactamente(&peptido) {
                encontrados.push(peptido);
            }
        } else if experimental.es_consistente(&peptido) {
            conservados.push(peptido);
        }
    }
    conservados
}

pub fn puntaje_ciclico(peptido: &[usize], espectro: &[usize]) -> Result<usize, String> {
    let mut experimental = Espectro::nuevo(espectro)?;
    Ok(experimental.compartidas(peptido, true))
}

pub fn puntaje_lineal(peptido: &[usize], espectro: &[usize]) -> Result<usize, String> {
    let mut experimental = Espectro::nuevo(espectro)?;
    Ok(experimental.compartidas(peptido, false))
}

pub fn recortar(
    tabla: Vec<Vec<usize>>,
    espectro: &[usize],
    n: usize,
) -> Result<Vec<Vec<usize>>, String> {
    let mut experimental = Espectro::nuevo(espectro)?;
    Ok(recortar_tabla(tabla, &mut experimental, n))
}

pub fn recortar_tabla(
    tabla: Vec<Vec<usize>>,
    experimental: &mut Espectro,
    n: usize,
) -> Vec<Vec<usize>> {
    if tabla.len() <= n {
        return tabla;
    }
    if n == 0 {
        return Vec::new();
    }
    let mut puntajes = Vec::with_capacity(tabla.len());
    for peptido in &tabla {
        puntajes.push(experimental.compartidas(peptido, false));
    }
    let mut de_mayor_a_menor = puntajes.clone();
    de_mayor_a_menor.sort_unstable_by(|a, b| b.cmp(a));
    let corte = de_mayor_a_menor[n - 1];
    let mut conservados = Vec::new();
    for (peptido, puntaje) in tabla.into_iter().zip(puntajes) {
        if puntaje >= corte {
            conservados.push(peptido);
        }
    }
    conservados
}

pub fn secuenciar_con_lideres(espectro: &[usize], n: usize) -> Result<Lideres, String> {
    secuenciar_con_alfabeto(espectro, n, &MASAS_DISTINTAS)
}

pub fn secuenciar_con_alfabeto(
    espectro: &[usize],
    n: usize,
    alfabeto: &[usize],
) -> Result<Lideres, String> {
    if alfabeto.is_empty() || alfabeto.contains(&0) {
        return Err("Hace falta al menos una masa de aminoácido mayor que cero.".to_string());
    }
    let mut experimental = Espectro::nuevo(espectro)?;
    let mut lideres = Lideres {
        alfabeto: alfabeto.to_vec(),
        puntaje: 0,
        peptidos: Vec::new(),
    };
    let mut tabla: Vec<Vec<usize>> = vec![Vec::new()];
    let mut revisados = 0;
    while !tabla.is_empty() {
        tabla = expandir(&tabla, alfabeto)?;
        revisados += tabla.len();
        if revisados > MAXIMO_CANDIDATOS_EN_TOTAL {
            return Err(
                "La tabla de líderes revisó más de dos millones de candidatos en total, el límite: hay muchos empates y la tabla no deja de crecer."
                    .to_string(),
            );
        }
        let mut en_carrera = Vec::with_capacity(tabla.len());
        for peptido in tabla {
            let masa = masa_total(&peptido);
            if masa == experimental.masa_madre {
                let puntaje = experimental.compartidas(&peptido, true);
                anotar_si_es_lider(&mut lideres, &peptido, puntaje);
                en_carrera.push(peptido);
            } else if masa < experimental.masa_madre {
                en_carrera.push(peptido);
            }
        }
        tabla = recortar_tabla(en_carrera, &mut experimental, n);
    }
    Ok(lideres)
}

fn anotar_si_es_lider(lideres: &mut Lideres, peptido: &[usize], puntaje: usize) {
    if lideres.peptidos.is_empty() || puntaje > lideres.puntaje {
        lideres.puntaje = puntaje;
        lideres.peptidos.clear();
        lideres.peptidos.push(peptido.to_vec());
    } else if puntaje == lideres.puntaje {
        lideres.peptidos.push(peptido.to_vec());
    }
}

pub fn convolucion_espectral(espectro: &[usize]) -> Result<Vec<(usize, usize)>, String> {
    let experimental = Espectro::nuevo(espectro)?;
    let masas = &experimental.masas;
    let mut repeticiones = vec![0usize; experimental.masa_madre + 1];
    for (posicion, &mayor) in masas.iter().enumerate() {
        for &menor in &masas[..posicion] {
            let diferencia = mayor - menor;
            if diferencia > 0 {
                repeticiones[diferencia] += 1;
            }
        }
    }
    let mut convolucion = Vec::new();
    for (diferencia, &cuantas) in repeticiones.iter().enumerate() {
        if cuantas > 0 {
            convolucion.push((diferencia, cuantas));
        }
    }
    convolucion.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    Ok(convolucion)
}

pub fn alfabeto_por_convolucion(espectro: &[usize], m: usize) -> Result<Vec<usize>, String> {
    let convolucion = convolucion_espectral(espectro)?;
    let mut en_rango = Vec::new();
    for (diferencia, cuantas) in convolucion {
        if (MASA_MINIMA_DE_AMINOACIDO..=MASA_MAXIMA_DE_AMINOACIDO).contains(&diferencia) {
            en_rango.push((diferencia, cuantas));
        }
    }
    if m == 0 || en_rango.is_empty() {
        return Ok(Vec::new());
    }
    let corte = en_rango[m.min(en_rango.len()) - 1].1;
    let mut alfabeto = Vec::new();
    for (diferencia, cuantas) in en_rango {
        if cuantas >= corte {
            alfabeto.push(diferencia);
        }
    }
    alfabeto.sort_unstable();
    Ok(alfabeto)
}

pub fn secuenciar_por_convolucion(
    espectro: &[usize],
    m: usize,
    n: usize,
) -> Result<Lideres, String> {
    if m == 0 {
        return Err("M debe ser al menos 1.".to_string());
    }
    let alfabeto = alfabeto_por_convolucion(espectro, m)?;
    if alfabeto.is_empty() {
        return Err(format!(
            "La convolución no tiene ninguna diferencia entre {MASA_MINIMA_DE_AMINOACIDO} y {MASA_MAXIMA_DE_AMINOACIDO} Da, así que no sugiere ningún aminoácido."
        ));
    }
    secuenciar_con_alfabeto(espectro, n, &alfabeto)
}

pub fn puntos_desde_distancias(distancias: &[i64]) -> Result<Option<Vec<i64>>, String> {
    let mut pendientes: BTreeMap<i64, usize> = BTreeMap::new();
    let mut cantidad = 0;
    for &distancia in distancias {
        if distancia > 0 {
            *pendientes.entry(distancia).or_insert(0) += 1;
            cantidad += 1;
        }
    }
    if cantidad == 0 {
        return Err("Hace falta al menos una distancia mayor que cero.".to_string());
    }
    let Some(cantidad_de_puntos) = puntos_para(cantidad) else {
        return Err(format!(
            "Con n puntos hay n(n−1)/2 distancias (1, 3, 6, 10, 15…), pero aquí hay {cantidad} distancias positivas."
        ));
    };
    if cantidad_de_puntos > MAXIMO_PUNTOS {
        return Err(format!(
            "Esas distancias corresponden a {cantidad_de_puntos} puntos y el límite es {MAXIMO_PUNTOS}: la búsqueda con retroceso podría tardar demasiado."
        ));
    }
    let ancho = *pendientes.keys().next_back().unwrap_or(&0);
    quitar(&mut pendientes, ancho);
    let mut puntos = vec![0, ancho];
    let mut intentos = 0;
    if colocar(&mut pendientes, &mut puntos, ancho, &mut intentos)? {
        puntos.sort_unstable();
        Ok(Some(puntos))
    } else {
        Ok(None)
    }
}

pub fn puntos_para(distancias: usize) -> Option<usize> {
    let mut n = 2;
    while n * (n - 1) / 2 < distancias {
        n += 1;
    }
    if n * (n - 1) / 2 == distancias {
        Some(n)
    } else {
        None
    }
}

fn colocar(
    pendientes: &mut BTreeMap<i64, usize>,
    puntos: &mut Vec<i64>,
    ancho: i64,
    intentos: &mut usize,
) -> Result<bool, String> {
    let Some(&mayor) = pendientes.keys().next_back() else {
        return Ok(true);
    };
    *intentos += 1;
    if *intentos > MAXIMO_INTENTOS {
        return Err(
            "La búsqueda con retroceso probó un millón de posiciones (el límite) sin terminar: estas distancias obligan a retroceder demasiadas veces. Revisa si sobra o falta alguna."
                .to_string(),
        );
    }
    let mut opciones = vec![mayor];
    if ancho - mayor != mayor {
        opciones.push(ancho - mayor);
    }
    for candidato in opciones {
        if quitar_distancias_a(candidato, puntos, pendientes) {
            puntos.push(candidato);
            if colocar(pendientes, puntos, ancho, intentos)? {
                return Ok(true);
            }
            puntos.pop();
            devolver_distancias_a(candidato, puntos, pendientes);
        }
    }
    Ok(false)
}

fn quitar_distancias_a(
    candidato: i64,
    puntos: &[i64],
    pendientes: &mut BTreeMap<i64, usize>,
) -> bool {
    let mut quitadas = Vec::with_capacity(puntos.len());
    for &punto in puntos {
        let distancia = (candidato - punto).abs();
        if quitar(pendientes, distancia) {
            quitadas.push(distancia);
        } else {
            for devuelta in quitadas {
                *pendientes.entry(devuelta).or_insert(0) += 1;
            }
            return false;
        }
    }
    true
}

fn devolver_distancias_a(candidato: i64, puntos: &[i64], pendientes: &mut BTreeMap<i64, usize>) {
    for &punto in puntos {
        *pendientes.entry((candidato - punto).abs()).or_insert(0) += 1;
    }
}

fn quitar(pendientes: &mut BTreeMap<i64, usize>, distancia: i64) -> bool {
    match pendientes.get_mut(&distancia) {
        Some(cuantas) => {
            *cuantas -= 1;
            if *cuantas == 0 {
                pendientes.remove(&distancia);
            }
            true
        }
        None => false,
    }
}
