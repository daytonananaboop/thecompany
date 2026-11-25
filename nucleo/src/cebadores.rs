use crate::secuencia::{complemento_inverso, indice_base};
use rayon::prelude::*;
use std::cmp::Ordering;

pub const CONSTANTE_DE_LOS_GASES: f64 = 1.987;
pub const CERO_ABSOLUTO: f64 = 273.15;
pub const TEMPERATURA_DE_ESTRUCTURAS: f64 = 37.0;
pub const LONGITUD_MAXIMA_DEL_MODELO: usize = 60;
pub const BASES_DEL_EXTREMO_3: usize = 15;
pub const UMBRAL_DIMERO: f64 = -9.0;
pub const UMBRAL_DIMERO_EXTREMO_3: f64 = -5.0;
pub const UMBRAL_HORQUILLA: f64 = -3.0;
pub const PLANTILLA_MAXIMA: usize = 20_000;
pub const PRODUCTO_MAXIMO: usize = 5_000;
pub const PARES_EVALUADOS_MAXIMOS: usize = 300_000;
pub const BASES_DE_LA_PINZA: usize = 5;
pub const PINZA_MINIMA: usize = 1;
pub const PINZA_MAXIMA: usize = 3;

const ENTROPIA_DE_SIMETRIA: f64 = -1.4;
const ENTROPIA_POR_FOSFATO: f64 = 0.368;
const BITS_DE_LA_MARCA: u32 = 24;
const TALLO_MINIMO: usize = 3;
const LAZO_MINIMO: usize = 3;
const PESO_TEMPERATURA: f64 = 1.0;
const PESO_LONGITUD: f64 = 1.0;
const PESO_DIFERENCIA_TM: f64 = 1.0;
const PESO_ESTRUCTURA: f64 = 1.0;
const PESO_SITIO_EXTRA: f64 = 5.0;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CorreccionSal {
    SantaLucia,
    Owczarzy,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Condiciones {
    pub monovalentes_mm: f64,
    pub magnesio_mm: f64,
    pub dntp_mm: f64,
    pub oligo_nm: f64,
    pub correccion: CorreccionSal,
}

impl Default for Condiciones {
    fn default() -> Condiciones {
        Condiciones {
            monovalentes_mm: 50.0,
            magnesio_mm: 1.5,
            dntp_mm: 0.6,
            oligo_nm: 50.0,
            correccion: CorreccionSal::SantaLucia,
        }
    }
}

pub fn apilamiento(primera: u8, segunda: u8) -> (f64, f64) {
    match (primera, segunda) {
        (b'A', b'A') | (b'T', b'T') => (-7.9, -22.2),
        (b'A', b'T') => (-7.2, -20.4),
        (b'T', b'A') => (-7.2, -21.3),
        (b'C', b'A') | (b'T', b'G') => (-8.5, -22.7),
        (b'G', b'T') | (b'A', b'C') => (-8.4, -22.4),
        (b'C', b'T') | (b'A', b'G') => (-7.8, -21.0),
        (b'G', b'A') | (b'T', b'C') => (-8.2, -22.2),
        (b'C', b'G') => (-10.6, -27.2),
        (b'G', b'C') => (-9.8, -24.4),
        (b'G', b'G') | (b'C', b'C') => (-8.0, -19.9),
        _ => (0.0, 0.0),
    }
}

pub fn iniciacion(base_terminal: u8) -> (f64, f64) {
    match base_terminal {
        b'G' | b'C' => (0.1, -2.8),
        _ => (2.3, 4.1),
    }
}

pub fn se_aparean(una: u8, otra: u8) -> bool {
    matches!(
        (una, otra),
        (b'A', b'T') | (b'T', b'A') | (b'C', b'G') | (b'G', b'C')
    )
}

pub fn es_autocomplementario(oligo: &[u8]) -> bool {
    let longitud = oligo.len();
    if longitud == 0 || !longitud.is_multiple_of(2) {
        return false;
    }
    (0..longitud / 2).all(|i| se_aparean(oligo[i], oligo[longitud - 1 - i]))
}

fn contar_gc(oligo: &[u8]) -> usize {
    oligo
        .iter()
        .filter(|&&base| base == b'G' || base == b'C')
        .count()
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Duplex {
    pub entalpia: f64,
    pub entropia: f64,
    pub longitud: usize,
    pub bases_gc: usize,
    pub autocomplementario: bool,
}

pub fn duplex(oligo: &[u8]) -> Duplex {
    let mut entalpia = 0.0;
    let mut entropia = 0.0;
    for par in oligo.windows(2) {
        let (entalpia_del_par, entropia_del_par) = apilamiento(par[0], par[1]);
        entalpia += entalpia_del_par;
        entropia += entropia_del_par;
    }
    if let (Some(&primera), Some(&ultima)) = (oligo.first(), oligo.last()) {
        let (entalpia_inicio, entropia_inicio) = iniciacion(primera);
        let (entalpia_final, entropia_final) = iniciacion(ultima);
        entalpia += entalpia_inicio + entalpia_final;
        entropia += entropia_inicio + entropia_final;
    }
    let autocomplementario = es_autocomplementario(oligo);
    if autocomplementario {
        entropia += ENTROPIA_DE_SIMETRIA;
    }
    Duplex {
        entalpia,
        entropia,
        longitud: oligo.len(),
        bases_gc: contar_gc(oligo),
        autocomplementario,
    }
}

pub fn magnesio_libre(condiciones: &Condiciones) -> f64 {
    (condiciones.magnesio_mm - condiciones.dntp_mm).max(0.0)
}

pub fn sodio_equivalente(condiciones: &Condiciones) -> f64 {
    condiciones.monovalentes_mm + 120.0 * magnesio_libre(condiciones).sqrt()
}

fn concentracion_efectiva(duplex: &Duplex, condiciones: &Condiciones) -> f64 {
    let oligo = condiciones.oligo_nm * 1e-9;
    if duplex.autocomplementario {
        oligo
    } else {
        oligo / 4.0
    }
}

fn tm_en_un_molar(duplex: &Duplex, condiciones: &Condiciones) -> f64 {
    let concentracion = concentracion_efectiva(duplex, condiciones);
    1000.0 * duplex.entalpia / (duplex.entropia + CONSTANTE_DE_LOS_GASES * concentracion.ln())
}

fn tm_santalucia(duplex: &Duplex, condiciones: &Condiciones) -> f64 {
    let sodio = sodio_equivalente(condiciones) / 1000.0;
    let fosfatos = duplex.longitud.saturating_sub(1) as f64;
    let entropia = duplex.entropia + ENTROPIA_POR_FOSFATO * fosfatos * sodio.ln();
    let concentracion = concentracion_efectiva(duplex, condiciones);
    1000.0 * duplex.entalpia / (entropia + CONSTANTE_DE_LOS_GASES * concentracion.ln())
        - CERO_ABSOLUTO
}

fn correccion_monovalente(fraccion_gc: f64, monovalentes: f64) -> f64 {
    let logaritmo = monovalentes.ln();
    (4.29 * fraccion_gc - 3.95) * 1e-5 * logaritmo + 9.40e-6 * logaritmo * logaritmo
}

fn correccion_owczarzy(duplex: &Duplex, condiciones: &Condiciones) -> f64 {
    let monovalentes = condiciones.monovalentes_mm / 1000.0;
    let magnesio = magnesio_libre(condiciones) / 1000.0;
    let fraccion_gc = duplex.bases_gc as f64 / duplex.longitud.max(1) as f64;
    let proporcion = if monovalentes > 0.0 {
        magnesio.sqrt() / monovalentes
    } else {
        f64::INFINITY
    };
    if proporcion < 0.22 {
        return correccion_monovalente(fraccion_gc, monovalentes);
    }
    let (mut a, b, c, mut d, e, f, mut g) = (
        3.92e-5, -9.11e-6, 6.26e-5, 1.42e-5, -4.82e-4, 5.25e-4, 8.31e-5,
    );
    if proporcion < 6.0 {
        let ln_mon = monovalentes.ln();
        a *= 0.843 - 0.352 * monovalentes.sqrt() * ln_mon;
        d *= 1.279 - 4.03e-3 * ln_mon - 8.03e-3 * ln_mon * ln_mon;
        g *= 0.486 - 0.258 * ln_mon + 5.25e-3 * ln_mon * ln_mon * ln_mon;
    }
    let ln_mg = magnesio.ln();
    let pares = duplex.longitud as f64;
    a + b * ln_mg
        + fraccion_gc * (c + d * ln_mg)
        + (e + f * ln_mg + g * ln_mg * ln_mg) / (2.0 * (pares - 1.0))
}

fn tm_owczarzy(duplex: &Duplex, condiciones: &Condiciones) -> f64 {
    let inversa =
        1.0 / tm_en_un_molar(duplex, condiciones) + correccion_owczarzy(duplex, condiciones);
    1.0 / inversa - CERO_ABSOLUTO
}

pub fn tm_del_duplex(duplex: &Duplex, condiciones: &Condiciones) -> f64 {
    match condiciones.correccion {
        CorreccionSal::SantaLucia => tm_santalucia(duplex, condiciones),
        CorreccionSal::Owczarzy => tm_owczarzy(duplex, condiciones),
    }
}

fn validar_condiciones(condiciones: &Condiciones) -> Result<(), String> {
    let concentraciones = [
        condiciones.monovalentes_mm,
        condiciones.magnesio_mm,
        condiciones.dntp_mm,
        condiciones.oligo_nm,
    ];
    if concentraciones
        .iter()
        .any(|valor| !valor.is_finite() || *valor < 0.0)
    {
        return Err("Las concentraciones deben ser números positivos o cero.".to_string());
    }
    if condiciones.oligo_nm <= 0.0 {
        return Err("La concentración del cebador debe ser mayor que cero.".to_string());
    }
    if sodio_equivalente(condiciones) <= 0.0 {
        return Err(
            "Sin sodio o potasio, y sin magnesio libre (más magnesio que dNTP), el ADN no forma una doble hélice estable."
                .to_string(),
        );
    }
    Ok(())
}

fn validar_oligo(oligo: &[u8]) -> Result<(), String> {
    if oligo.len() < 2 {
        return Err("El oligo necesita al menos 2 bases para tener vecinos.".to_string());
    }
    if oligo.len() > LONGITUD_MAXIMA_DEL_MODELO {
        return Err(format!(
            "El modelo de vecinos cercanos vale para oligos de hasta {LONGITUD_MAXIMA_DEL_MODELO} bases y este mide {}.",
            oligo.len()
        ));
    }
    if let Some(letra) = oligo
        .iter()
        .find(|&&base| !matches!(base, b'A' | b'C' | b'G' | b'T'))
    {
        return Err(format!(
            "El oligo solo puede tener A, C, G y T, pero contiene «{}».",
            *letra as char
        ));
    }
    Ok(())
}

pub fn temperatura_fusion(oligo: &[u8], condiciones: &Condiciones) -> Result<f64, String> {
    validar_oligo(oligo)?;
    validar_condiciones(condiciones)?;
    let tm = tm_del_duplex(&duplex(oligo), condiciones);
    if !tm.is_finite() || tm <= -CERO_ABSOLUTO {
        return Err(
            "Con estas concentraciones el modelo no da una temperatura válida.".to_string(),
        );
    }
    Ok(tm)
}

fn fraccion_en_equilibrio(a: f64) -> f64 {
    if !a.is_finite() {
        return 1.0;
    }
    2.0 * a / (2.0 * a + 1.0 + (4.0 * a + 1.0).sqrt())
}

pub fn curva_de_fusion(oligo: &[u8], condiciones: &Condiciones) -> Result<Vec<(f64, f64)>, String> {
    let tm = temperatura_fusion(oligo, condiciones)?;
    let duplex = duplex(oligo);
    let oligo_molar = condiciones.oligo_nm * 1e-9;
    let entalpia = 1000.0 * duplex.entalpia;
    let entropia = entalpia / (tm + CERO_ABSOLUTO)
        - CONSTANTE_DE_LOS_GASES * concentracion_efectiva(&duplex, condiciones).ln();
    let mut curva = Vec::with_capacity(101);
    for paso in 0..=100 {
        let temperatura = tm - 25.0 + paso as f64 * 0.5;
        let kelvin = temperatura + CERO_ABSOLUTO;
        let energia_libre = entalpia - kelvin * entropia;
        let constante = (-energia_libre / (CONSTANTE_DE_LOS_GASES * kelvin)).exp();
        let a = if duplex.autocomplementario {
            2.0 * constante * oligo_molar
        } else {
            constante * oligo_molar / 2.0
        };
        curva.push((temperatura, fraccion_en_equilibrio(a)));
    }
    Ok(curva)
}

pub fn correccion_de_entropia(condiciones: &Condiciones) -> f64 {
    ENTROPIA_POR_FOSFATO * (sodio_equivalente(condiciones) / 1000.0).ln()
}

fn energia_libre(entalpia: f64, entropia: f64) -> f64 {
    entalpia - (TEMPERATURA_DE_ESTRUCTURAS + CERO_ABSOLUTO) * entropia / 1000.0
}

fn energia_de_apilamiento(primera: u8, segunda: u8, sal: f64) -> f64 {
    let (entalpia, entropia) = apilamiento(primera, segunda);
    energia_libre(entalpia, entropia + sal)
}

fn energia_de_iniciacion(base_terminal: u8) -> f64 {
    let (entalpia, entropia) = iniciacion(base_terminal);
    energia_libre(entalpia, entropia)
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Apareamiento {
    pub energia: f64,
    pub en_a: usize,
    pub en_b: usize,
    pub pares: usize,
}

impl Apareamiento {
    pub fn ninguno() -> Apareamiento {
        Apareamiento {
            energia: 0.0,
            en_a: 0,
            en_b: 0,
            pares: 0,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Dimero {
    pub mas_estable: Apareamiento,
    pub extremo_3: Apareamiento,
}

fn es_base(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T')
}

fn anotar_tramo(
    dimero: &mut Dimero,
    a: &[u8],
    b: &[u8],
    en_a: usize,
    en_espejo: usize,
    pares: usize,
    energia: f64,
) {
    let tramo = Apareamiento {
        energia,
        en_a,
        en_b: b.len() - 1 - en_espejo,
        pares,
    };
    if tramo.energia < dimero.mas_estable.energia {
        dimero.mas_estable = tramo;
    }
    let toca_un_extremo_3 = en_a + pares == a.len() || en_espejo == 0;
    if toca_un_extremo_3 && tramo.energia < dimero.extremo_3.energia {
        dimero.extremo_3 = tramo;
    }
}

fn anotar_bloque(
    dimero: &mut Dimero,
    a: &[u8],
    b: &[u8],
    en_a: usize,
    en_espejo: usize,
    pares: usize,
    sal: f64,
) {
    let mut apilado = 0.0;
    let mut mejor_inicio = (energia_de_iniciacion(a[en_a]), 0);
    for ultimo in 1..pares {
        apilado += energia_de_apilamiento(a[en_a + ultimo - 1], a[en_a + ultimo], sal);
        let final_del_tramo = energia_de_iniciacion(a[en_a + ultimo]);
        let (valor_inicio, inicio) = mejor_inicio;
        let largo = ultimo - inicio + 1;
        let energia = valor_inicio + final_del_tramo + apilado;
        anotar_tramo(
            dimero,
            a,
            b,
            en_a + inicio,
            en_espejo + inicio,
            largo,
            energia,
        );
        if en_espejo == 0 && inicio != 0 {
            let desde_el_borde = energia_de_iniciacion(a[en_a]) + final_del_tramo + apilado;
            anotar_tramo(dimero, a, b, en_a, 0, ultimo + 1, desde_el_borde);
        }
        let candidato = final_del_tramo - apilado;
        if candidato < mejor_inicio.0 {
            mejor_inicio = (candidato, ultimo);
        }
    }
}

fn buscar_dimero(a: &[u8], b: &[u8], sal: f64) -> Dimero {
    let mut dimero = Dimero {
        mas_estable: Apareamiento::ninguno(),
        extremo_3: Apareamiento::ninguno(),
    };
    if a.is_empty() || b.is_empty() {
        return dimero;
    }
    let espejo = complemento_inverso(b);
    for desplazamiento in 0..a.len() + b.len() - 1 {
        let inicio_a = desplazamiento.saturating_sub(b.len() - 1);
        let inicio_espejo = (b.len() - 1).saturating_sub(desplazamiento);
        let enfrentadas = (a.len() - inicio_a).min(b.len() - inicio_espejo);
        let mut racha = 0;
        for i in 0..enfrentadas {
            let base = a[inicio_a + i];
            let iguales = (base == espejo[inicio_espejo + i]) & es_base(base);
            if !iguales & (racha >= 2) {
                let desde = i - racha;
                anotar_bloque(
                    &mut dimero,
                    a,
                    b,
                    inicio_a + desde,
                    inicio_espejo + desde,
                    racha,
                    sal,
                );
            }
            racha = if iguales { racha + 1 } else { 0 };
        }
        if racha >= 2 {
            let desde = enfrentadas - racha;
            anotar_bloque(
                &mut dimero,
                a,
                b,
                inicio_a + desde,
                inicio_espejo + desde,
                racha,
                sal,
            );
        }
    }
    dimero
}

pub fn dimero(a: &[u8], b: &[u8], condiciones: &Condiciones) -> Result<Dimero, String> {
    validar_condiciones(condiciones)?;
    Ok(buscar_dimero(a, b, correccion_de_entropia(condiciones)))
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Horquilla {
    pub energia: f64,
    pub inicio: usize,
    pub fin: usize,
    pub pares: usize,
    pub lazo: usize,
}

impl Horquilla {
    pub fn ninguna() -> Horquilla {
        Horquilla {
            energia: 0.0,
            inicio: 0,
            fin: 0,
            pares: 0,
            lazo: 0,
        }
    }
}

pub fn energia_de_lazo(bases: usize) -> f64 {
    match bases {
        0..=4 => 3.5,
        5 => 3.3,
        6 => 4.0,
        7 => 4.2,
        8 => 4.3,
        9 => 4.5,
        10 => 4.6,
        _ => {
            let kelvin = TEMPERATURA_DE_ESTRUCTURAS + CERO_ABSOLUTO;
            4.6 + 2.44 * CONSTANTE_DE_LOS_GASES * kelvin / 1000.0 * (bases as f64 / 10.0).ln()
        }
    }
}

fn buscar_horquilla(oligo: &[u8], sal: f64) -> Horquilla {
    let mut mejor = Horquilla::ninguna();
    for interior_5 in 0..oligo.len() {
        for interior_3 in interior_5 + LAZO_MINIMO + 1..oligo.len() {
            if !se_aparean(oligo[interior_5], oligo[interior_3]) {
                continue;
            }
            let lazo = interior_3 - interior_5 - 1;
            let mut energia = energia_de_lazo(lazo);
            let mut pares = 1;
            while pares <= interior_5
                && interior_3 + pares < oligo.len()
                && se_aparean(oligo[interior_5 - pares], oligo[interior_3 + pares])
            {
                let base = interior_5 - pares;
                energia += energia_de_apilamiento(oligo[base], oligo[base + 1], sal);
                pares += 1;
                if pares >= TALLO_MINIMO && energia < mejor.energia {
                    mejor = Horquilla {
                        energia,
                        inicio: interior_5 + 1 - pares,
                        fin: interior_3 + pares - 1,
                        pares,
                        lazo,
                    };
                }
            }
        }
    }
    mejor
}

pub fn horquilla(oligo: &[u8], condiciones: &Condiciones) -> Result<Horquilla, String> {
    validar_condiciones(condiciones)?;
    Ok(buscar_horquilla(oligo, correccion_de_entropia(condiciones)))
}

fn codificar(kmer: &[u8]) -> Option<u64> {
    let mut codigo = 0u64;
    for &base in kmer {
        codigo = (codigo << 2) | indice_base(base)? as u64;
    }
    Some(codigo)
}

fn lugar_de_la_marca(codigo: u64) -> (usize, u64) {
    let huella = (codigo & ((1 << BITS_DE_LA_MARCA) - 1)) as usize;
    (huella / 64, 1 << (huella % 64))
}

fn marcar(marcas: &mut [u64], codigo: u64) {
    let (palabra, bit) = lugar_de_la_marca(codigo);
    marcas[palabra] |= bit;
}

fn esta_marcado(marcas: &[u64], codigo: u64) -> bool {
    let (palabra, bit) = lugar_de_la_marca(codigo);
    marcas[palabra] & bit != 0
}

fn contar_en_tramo(tramo: &[u8], k: usize, buscados: &[u64], marcas: &[u64]) -> Vec<usize> {
    let mascara = (1u64 << (2 * k)) - 1;
    let mut conteos = vec![0; buscados.len()];
    let mut codigo = 0u64;
    let mut seguidas = 0;
    for &base in tramo {
        let Some(valor) = indice_base(base) else {
            seguidas = 0;
            continue;
        };
        codigo = ((codigo << 2) | valor as u64) & mascara;
        seguidas += 1;
        if seguidas >= k
            && esta_marcado(marcas, codigo)
            && let Ok(ranura) = buscados.binary_search(&codigo)
        {
            conteos[ranura] += 1;
        }
    }
    conteos
}

fn sumar_conteos(mut uno: Vec<usize>, otro: Vec<usize>) -> Vec<usize> {
    for (total, parcial) in uno.iter_mut().zip(otro) {
        *total += parcial;
    }
    uno
}

pub fn contar_sitios(genoma: &[u8], extremos: &[&[u8]]) -> Vec<usize> {
    let Some(primero) = extremos.first() else {
        return Vec::new();
    };
    let k = primero.len();
    if k == 0 || k > 31 || extremos.iter().any(|extremo| extremo.len() != k) {
        return vec![0; extremos.len()];
    }
    let codigos: Vec<(Option<u64>, Option<u64>)> = extremos
        .iter()
        .map(|extremo| (codificar(extremo), codificar(&complemento_inverso(extremo))))
        .collect();
    let mut buscados: Vec<u64> = Vec::with_capacity(2 * codigos.len());
    for (directo, inverso) in &codigos {
        buscados.extend(directo);
        buscados.extend(inverso);
    }
    buscados.sort_unstable();
    buscados.dedup();
    let mut marcas = vec![0u64; (1 << BITS_DE_LA_MARCA) / 64];
    for &codigo in &buscados {
        marcar(&mut marcas, codigo);
    }
    let partes = 4 * rayon::current_num_threads();
    let largo_de_parte = genoma.len().div_ceil(partes).max(1);
    let conteos = (0..partes)
        .into_par_iter()
        .map(|parte| {
            let desde = (parte * largo_de_parte).min(genoma.len());
            let hasta = (desde + largo_de_parte + k - 1).min(genoma.len());
            contar_en_tramo(&genoma[desde..hasta], k, &buscados, &marcas)
        })
        .reduce(|| vec![0; buscados.len()], sumar_conteos);
    let veces = |codigo: &u64| {
        buscados
            .binary_search(codigo)
            .map_or(0, |ranura| conteos[ranura])
    };
    codigos
        .iter()
        .map(|par| match par {
            (Some(directo), Some(inverso)) if directo == inverso => veces(directo),
            (Some(directo), Some(inverso)) => veces(directo) + veces(inverso),
            _ => 0,
        })
        .collect()
}

#[derive(Debug, Clone, PartialEq)]
pub struct OpcionesDiseno {
    pub longitud_minima: usize,
    pub longitud_optima: usize,
    pub longitud_maxima: usize,
    pub tm_minima: f64,
    pub tm_optima: f64,
    pub tm_maxima: f64,
    pub gc_minimo: f64,
    pub gc_maximo: f64,
    pub repeticion_maxima: usize,
    pub producto_minimo: usize,
    pub producto_maximo: usize,
    pub objetivo: Option<(usize, usize)>,
    pub diferencia_tm_maxima: f64,
    pub pares: usize,
    pub condiciones: Condiciones,
}

impl OpcionesDiseno {
    pub fn recomendadas() -> OpcionesDiseno {
        OpcionesDiseno {
            longitud_minima: 18,
            longitud_optima: 20,
            longitud_maxima: 25,
            tm_minima: 57.0,
            tm_optima: 60.0,
            tm_maxima: 63.0,
            gc_minimo: 0.40,
            gc_maximo: 0.60,
            repeticion_maxima: 4,
            producto_minimo: 100,
            producto_maximo: 300,
            objetivo: None,
            diferencia_tm_maxima: 3.0,
            pares: 5,
            condiciones: Condiciones::default(),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Cebador {
    pub inicio: usize,
    pub longitud: usize,
    pub secuencia: Vec<u8>,
    pub tm: f64,
    pub gc: f64,
    pub autodimero: f64,
    pub autodimero_extremo_3: f64,
    pub horquilla: f64,
    pub sitios: Option<usize>,
    pub penalizacion: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ParDeCebadores {
    pub izquierdo: Cebador,
    pub derecho: Cebador,
    pub producto: usize,
    pub diferencia_tm: f64,
    pub dimero: f64,
    pub dimero_extremo_3: f64,
    pub penalizacion: f64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Descartes {
    pub bases_ambiguas: usize,
    pub repeticiones: usize,
    pub contenido_gc: usize,
    pub pinza_gc: usize,
    pub temperatura: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Diseno {
    pub pares: Vec<ParDeCebadores>,
    pub izquierdos: usize,
    pub derechos: usize,
    pub descartes: Descartes,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Hebra {
    Directa,
    Inversa,
}

enum Motivo {
    BasesAmbiguas,
    Repeticiones,
    ContenidoGc,
    PinzaGc,
    Temperatura,
}

impl Descartes {
    fn anotar(&mut self, motivo: Motivo) {
        match motivo {
            Motivo::BasesAmbiguas => self.bases_ambiguas += 1,
            Motivo::Repeticiones => self.repeticiones += 1,
            Motivo::ContenidoGc => self.contenido_gc += 1,
            Motivo::PinzaGc => self.pinza_gc += 1,
            Motivo::Temperatura => self.temperatura += 1,
        }
    }
}

struct SumasAcumuladas {
    entalpia: Vec<f64>,
    entropia: Vec<f64>,
    gc: Vec<usize>,
    ambiguas: Vec<usize>,
    rachas_largas: Vec<usize>,
}

fn sumas_acumuladas(plantilla: &[u8], repeticion_maxima: usize) -> SumasAcumuladas {
    let n = plantilla.len();
    let mut sumas = SumasAcumuladas {
        entalpia: vec![0.0; n + 1],
        entropia: vec![0.0; n + 1],
        gc: vec![0; n + 1],
        ambiguas: vec![0; n + 1],
        rachas_largas: vec![0; n + 1],
    };
    let mut racha = 0;
    for (posicion, &base) in plantilla.iter().enumerate() {
        if posicion > 0 && plantilla[posicion - 1] == base {
            racha += 1;
        } else {
            racha = 1;
        }
        let es_gc = base == b'G' || base == b'C';
        let es_ambigua = !matches!(base, b'A' | b'C' | b'G' | b'T');
        sumas.gc[posicion + 1] = sumas.gc[posicion] + usize::from(es_gc);
        sumas.ambiguas[posicion + 1] = sumas.ambiguas[posicion] + usize::from(es_ambigua);
        sumas.rachas_largas[posicion + 1] =
            sumas.rachas_largas[posicion] + usize::from(racha > repeticion_maxima);
        let (entalpia, entropia) = match plantilla.get(posicion + 1) {
            Some(&siguiente) => apilamiento(base, siguiente),
            None => (0.0, 0.0),
        };
        sumas.entalpia[posicion + 1] = sumas.entalpia[posicion] + entalpia;
        sumas.entropia[posicion + 1] = sumas.entropia[posicion] + entropia;
    }
    sumas
}

fn duplex_de_ventana(
    plantilla: &[u8],
    sumas: &SumasAcumuladas,
    inicio: usize,
    longitud: usize,
) -> Duplex {
    let fin = inicio + longitud;
    let (entalpia_inicio, entropia_inicio) = iniciacion(plantilla[inicio]);
    let (entalpia_final, entropia_final) = iniciacion(plantilla[fin - 1]);
    let autocomplementario = es_autocomplementario(&plantilla[inicio..fin]);
    let mut entropia =
        sumas.entropia[fin - 1] - sumas.entropia[inicio] + entropia_inicio + entropia_final;
    if autocomplementario {
        entropia += ENTROPIA_DE_SIMETRIA;
    }
    Duplex {
        entalpia: sumas.entalpia[fin - 1] - sumas.entalpia[inicio]
            + entalpia_inicio
            + entalpia_final,
        entropia,
        longitud,
        bases_gc: sumas.gc[fin] - sumas.gc[inicio],
        autocomplementario,
    }
}

struct Ventana {
    inicio: usize,
    longitud: usize,
    tm: f64,
}

fn revisar_ventana(
    plantilla: &[u8],
    sumas: &SumasAcumuladas,
    inicio: usize,
    longitud: usize,
    hebra: Hebra,
    opciones: &OpcionesDiseno,
) -> Result<f64, Motivo> {
    let fin = inicio + longitud;
    if sumas.ambiguas[fin] > sumas.ambiguas[inicio] {
        return Err(Motivo::BasesAmbiguas);
    }
    if sumas.rachas_largas[fin] > sumas.rachas_largas[inicio + opciones.repeticion_maxima] {
        return Err(Motivo::Repeticiones);
    }
    let gc = (sumas.gc[fin] - sumas.gc[inicio]) as f64 / longitud as f64;
    if gc < opciones.gc_minimo || gc > opciones.gc_maximo {
        return Err(Motivo::ContenidoGc);
    }
    let (desde, hasta) = match hebra {
        Hebra::Directa => (fin - BASES_DE_LA_PINZA, fin),
        Hebra::Inversa => (inicio, inicio + BASES_DE_LA_PINZA),
    };
    let pinza = sumas.gc[hasta] - sumas.gc[desde];
    if !(PINZA_MINIMA..=PINZA_MAXIMA).contains(&pinza) {
        return Err(Motivo::PinzaGc);
    }
    let duplex = duplex_de_ventana(plantilla, sumas, inicio, longitud);
    let tm = tm_del_duplex(&duplex, &opciones.condiciones);
    if !(opciones.tm_minima..=opciones.tm_maxima).contains(&tm) {
        return Err(Motivo::Temperatura);
    }
    Ok(tm)
}

fn buscar_ventanas(
    plantilla: &[u8],
    sumas: &SumasAcumuladas,
    hebra: Hebra,
    opciones: &OpcionesDiseno,
    descartes: &mut Descartes,
) -> Vec<Ventana> {
    let mut aceptadas = Vec::new();
    for inicio in 0..plantilla.len() {
        for longitud in opciones.longitud_minima..=opciones.longitud_maxima {
            if inicio + longitud > plantilla.len() {
                break;
            }
            match revisar_ventana(plantilla, sumas, inicio, longitud, hebra, opciones) {
                Ok(tm) => aceptadas.push(Ventana {
                    inicio,
                    longitud,
                    tm,
                }),
                Err(motivo) => descartes.anotar(motivo),
            }
        }
    }
    aceptadas
}

fn exceso(energia: f64, umbral: f64) -> f64 {
    (umbral - energia).max(0.0)
}

fn penalizar_cebador(cebador: &Cebador, opciones: &OpcionesDiseno) -> f64 {
    let estructura = exceso(cebador.autodimero, UMBRAL_DIMERO)
        + exceso(cebador.autodimero_extremo_3, UMBRAL_DIMERO_EXTREMO_3)
        + exceso(cebador.horquilla, UMBRAL_HORQUILLA);
    let sitios_de_mas = cebador.sitios.unwrap_or(1).saturating_sub(1);
    PESO_TEMPERATURA * (cebador.tm - opciones.tm_optima).abs()
        + PESO_LONGITUD * cebador.longitud.abs_diff(opciones.longitud_optima) as f64
        + PESO_ESTRUCTURA * estructura
        + PESO_SITIO_EXTRA * sitios_de_mas as f64
}

fn construir_cebador(
    plantilla: &[u8],
    ventana: &Ventana,
    hebra: Hebra,
    opciones: &OpcionesDiseno,
    sal: f64,
) -> Cebador {
    let tramo = &plantilla[ventana.inicio..ventana.inicio + ventana.longitud];
    let secuencia = match hebra {
        Hebra::Directa => tramo.to_vec(),
        Hebra::Inversa => complemento_inverso(tramo),
    };
    let propio = buscar_dimero(&secuencia, &secuencia, sal);
    let pliegue = buscar_horquilla(&secuencia, sal);
    let mut cebador = Cebador {
        inicio: ventana.inicio,
        longitud: ventana.longitud,
        gc: contar_gc(&secuencia) as f64 / ventana.longitud as f64,
        secuencia,
        tm: ventana.tm,
        autodimero: propio.mas_estable.energia,
        autodimero_extremo_3: propio.extremo_3.energia,
        horquilla: pliegue.energia,
        sitios: None,
        penalizacion: 0.0,
    };
    cebador.penalizacion = penalizar_cebador(&cebador, opciones);
    cebador
}

fn construir_cebadores(
    plantilla: &[u8],
    ventanas: &[Ventana],
    hebra: Hebra,
    opciones: &OpcionesDiseno,
) -> Vec<Cebador> {
    let sal = correccion_de_entropia(&opciones.condiciones);
    ventanas
        .par_iter()
        .map(|ventana| construir_cebador(plantilla, ventana, hebra, opciones, sal))
        .collect()
}

fn extremo_en_la_plantilla(cebador: &Cebador, hebra: Hebra) -> (usize, usize) {
    match hebra {
        Hebra::Directa => {
            let fin = cebador.inicio + cebador.longitud;
            (fin - BASES_DEL_EXTREMO_3, fin)
        }
        Hebra::Inversa => (cebador.inicio, cebador.inicio + BASES_DEL_EXTREMO_3),
    }
}

fn anotar_sitios(
    plantilla: &[u8],
    izquierdos: &mut [Cebador],
    derechos: &mut [Cebador],
    genoma: &[u8],
    opciones: &OpcionesDiseno,
) {
    let mut extremos: Vec<&[u8]> = Vec::with_capacity(izquierdos.len() + derechos.len());
    for cebador in izquierdos.iter() {
        let (desde, hasta) = extremo_en_la_plantilla(cebador, Hebra::Directa);
        extremos.push(&plantilla[desde..hasta]);
    }
    for cebador in derechos.iter() {
        let (desde, hasta) = extremo_en_la_plantilla(cebador, Hebra::Inversa);
        extremos.push(&plantilla[desde..hasta]);
    }
    let sitios = contar_sitios(genoma, &extremos);
    for (cebador, cuenta) in izquierdos.iter_mut().chain(derechos.iter_mut()).zip(sitios) {
        cebador.sitios = Some(cuenta);
        cebador.penalizacion = penalizar_cebador(cebador, opciones);
    }
}

fn producto_valido(
    izquierdo: &Cebador,
    derecho: &Cebador,
    opciones: &OpcionesDiseno,
) -> Option<usize> {
    let fin_izquierdo = izquierdo.inicio + izquierdo.longitud;
    if derecho.inicio < fin_izquierdo {
        return None;
    }
    let producto = derecho.inicio + derecho.longitud - izquierdo.inicio;
    if producto < opciones.producto_minimo || producto > opciones.producto_maximo {
        return None;
    }
    if let Some((inicio, longitud)) = opciones.objetivo
        && (fin_izquierdo > inicio || derecho.inicio < inicio + longitud)
    {
        return None;
    }
    if (izquierdo.tm - derecho.tm).abs() > opciones.diferencia_tm_maxima {
        return None;
    }
    Some(producto)
}

fn redondear(penalizacion: f64) -> f64 {
    (penalizacion * 10_000.0).round() / 10_000.0
}

fn penalizacion_sin_dimero(izquierdo: &Cebador, derecho: &Cebador) -> f64 {
    izquierdo.penalizacion
        + derecho.penalizacion
        + PESO_DIFERENCIA_TM * (izquierdo.tm - derecho.tm).abs()
}

fn evaluar_par(
    izquierdo: &Cebador,
    derecho: &Cebador,
    producto: usize,
    sal: f64,
) -> ParDeCebadores {
    let cruzado = buscar_dimero(&izquierdo.secuencia, &derecho.secuencia, sal);
    let diferencia_tm = (izquierdo.tm - derecho.tm).abs();
    let estructura = exceso(cruzado.mas_estable.energia, UMBRAL_DIMERO)
        + exceso(cruzado.extremo_3.energia, UMBRAL_DIMERO_EXTREMO_3);
    ParDeCebadores {
        izquierdo: izquierdo.clone(),
        derecho: derecho.clone(),
        producto,
        diferencia_tm,
        dimero: cruzado.mas_estable.energia,
        dimero_extremo_3: cruzado.extremo_3.energia,
        penalizacion: redondear(
            penalizacion_sin_dimero(izquierdo, derecho) + PESO_ESTRUCTURA * estructura,
        ),
    }
}

fn comparar_pares(uno: &ParDeCebadores, otro: &ParDeCebadores) -> Ordering {
    uno.penalizacion
        .total_cmp(&otro.penalizacion)
        .then(uno.izquierdo.inicio.cmp(&otro.izquierdo.inicio))
        .then(uno.izquierdo.longitud.cmp(&otro.izquierdo.longitud))
        .then(uno.derecho.inicio.cmp(&otro.derecho.inicio))
        .then(uno.derecho.longitud.cmp(&otro.derecho.longitud))
}

fn guardar_si_es_de_los_mejores(
    mejores: &mut Vec<ParDeCebadores>,
    par: ParDeCebadores,
    cupo: usize,
) {
    let lugar = mejores.partition_point(|otro| comparar_pares(otro, &par) == Ordering::Less);
    if lugar < cupo {
        mejores.insert(lugar, par);
        mejores.truncate(cupo);
    }
}

fn peor_penalizacion(mejores: &[ParDeCebadores], cupo: usize) -> f64 {
    if mejores.len() < cupo {
        return f64::INFINITY;
    }
    mejores.last().map_or(f64::INFINITY, |par| par.penalizacion)
}

fn no_puede_entrar(
    mejores: &[ParDeCebadores],
    cupo: usize,
    penalizacion: f64,
    izquierdo: &Cebador,
    derecho: Option<&Cebador>,
) -> bool {
    if mejores.len() < cupo {
        return false;
    }
    let Some(peor) = mejores.last() else {
        return false;
    };
    let mut orden = penalizacion
        .total_cmp(&peor.penalizacion)
        .then(izquierdo.inicio.cmp(&peor.izquierdo.inicio))
        .then(izquierdo.longitud.cmp(&peor.izquierdo.longitud));
    if let Some(derecho) = derecho {
        orden = orden
            .then(derecho.inicio.cmp(&peor.derecho.inicio))
            .then(derecho.longitud.cmp(&peor.derecho.longitud));
    }
    orden == Ordering::Greater
}

fn derechos_al_alcance<'a>(
    derechos: &'a [Cebador],
    izquierdo: &Cebador,
    opciones: &OpcionesDiseno,
) -> &'a [Cebador] {
    let primer_inicio =
        (izquierdo.inicio + opciones.producto_minimo).saturating_sub(opciones.longitud_maxima);
    let ultimo_inicio = izquierdo.inicio + opciones.producto_maximo;
    let desde = derechos.partition_point(|derecho| derecho.inicio < primer_inicio);
    let hasta = derechos.partition_point(|derecho| derecho.inicio < ultimo_inicio);
    &derechos[desde..hasta.max(desde)]
}

fn emparejar(
    izquierdos: &[Cebador],
    derechos: &[Cebador],
    opciones: &OpcionesDiseno,
) -> Result<Vec<ParDeCebadores>, String> {
    let sal = correccion_de_entropia(&opciones.condiciones);
    let mut evaluados = 0;
    let mut por_penalizacion: Vec<&Cebador> = izquierdos.iter().collect();
    por_penalizacion.sort_by(|uno, otro| uno.penalizacion.total_cmp(&otro.penalizacion));
    let menor_derecho = derechos
        .iter()
        .map(|derecho| derecho.penalizacion)
        .fold(f64::INFINITY, f64::min);
    let mut mejores: Vec<ParDeCebadores> = Vec::new();
    for izquierdo in por_penalizacion {
        let mejor_posible = redondear(izquierdo.penalizacion + menor_derecho);
        if mejor_posible > peor_penalizacion(&mejores, opciones.pares) {
            break;
        }
        if no_puede_entrar(&mejores, opciones.pares, mejor_posible, izquierdo, None) {
            continue;
        }
        for derecho in derechos_al_alcance(derechos, izquierdo, opciones) {
            let Some(producto) = producto_valido(izquierdo, derecho, opciones) else {
                continue;
            };
            let minimo_posible = redondear(penalizacion_sin_dimero(izquierdo, derecho));
            if no_puede_entrar(
                &mejores,
                opciones.pares,
                minimo_posible,
                izquierdo,
                Some(derecho),
            ) {
                continue;
            }
            evaluados += 1;
            if evaluados > PARES_EVALUADOS_MAXIMOS {
                return Err(format!(
                    "Hay demasiados pares de cebadores casi iguales (más de {PARES_EVALUADOS_MAXIMOS} por revisar), algo típico de una zona repetitiva. Recorta la plantilla para dejar fuera la repetición."
                ));
            }
            let par = evaluar_par(izquierdo, derecho, producto, sal);
            guardar_si_es_de_los_mejores(&mut mejores, par, opciones.pares);
        }
    }
    Ok(mejores)
}

fn validar_diseno(plantilla: &[u8], opciones: &OpcionesDiseno) -> Result<(), String> {
    validar_condiciones(&opciones.condiciones)?;
    if plantilla.is_empty() {
        return Err("La plantilla está vacía.".to_string());
    }
    if plantilla.len() > PLANTILLA_MAXIMA {
        return Err(format!(
            "La plantilla mide {} pb y el límite es {PLANTILLA_MAXIMA} pb: recorta la región que quieres amplificar dejando unos cientos de bases a cada lado.",
            plantilla.len()
        ));
    }
    let longitudes_en_orden = BASES_DEL_EXTREMO_3 <= opciones.longitud_minima
        && opciones.longitud_minima <= opciones.longitud_optima
        && opciones.longitud_optima <= opciones.longitud_maxima
        && opciones.longitud_maxima <= LONGITUD_MAXIMA_DEL_MODELO;
    if !longitudes_en_orden {
        return Err(format!(
            "Las longitudes deben cumplir {BASES_DEL_EXTREMO_3} ≤ mínima ≤ óptima ≤ máxima ≤ {LONGITUD_MAXIMA_DEL_MODELO}."
        ));
    }
    let temperaturas = [opciones.tm_minima, opciones.tm_optima, opciones.tm_maxima];
    if temperaturas.iter().any(|valor| !valor.is_finite())
        || opciones.tm_minima > opciones.tm_optima
        || opciones.tm_optima > opciones.tm_maxima
    {
        return Err("Las temperaturas deben cumplir mínima ≤ óptima ≤ máxima.".to_string());
    }
    if !(0.0 <= opciones.gc_minimo
        && opciones.gc_minimo <= opciones.gc_maximo
        && opciones.gc_maximo <= 1.0)
    {
        return Err("El contenido GC debe cumplir 0 ≤ mínimo ≤ máximo ≤ 1.".to_string());
    }
    if opciones.repeticion_maxima == 0 || opciones.repeticion_maxima >= opciones.longitud_minima {
        return Err(
            "La racha máxima de una misma base debe ser al menos 1 y menor que el cebador."
                .to_string(),
        );
    }
    if opciones.producto_minimo == 0 || opciones.producto_minimo > opciones.producto_maximo {
        return Err("El tamaño del producto debe cumplir 1 ≤ mínimo ≤ máximo.".to_string());
    }
    if opciones.producto_maximo > PRODUCTO_MAXIMO {
        return Err(format!(
            "El producto puede medir hasta {PRODUCTO_MAXIMO} pb: la Taq común no copia bien fragmentos más largos."
        ));
    }
    if plantilla.len() < opciones.producto_minimo {
        return Err(format!(
            "La plantilla mide {} pb, menos que el producto mínimo de {} pb.",
            plantilla.len(),
            opciones.producto_minimo
        ));
    }
    if let Some((inicio, longitud)) = opciones.objetivo {
        if longitud == 0 || inicio.saturating_add(longitud) > plantilla.len() {
            return Err(
                "La región objetivo debe tener al menos 1 base y quedar dentro de la plantilla."
                    .to_string(),
            );
        }
        if longitud + 2 * opciones.longitud_minima > opciones.producto_maximo {
            return Err(format!(
                "La región objetivo ({longitud} pb) más dos cebadores no cabe en un producto de hasta {} pb.",
                opciones.producto_maximo
            ));
        }
    }
    if opciones.pares == 0 || opciones.pares > 100 {
        return Err("Pide entre 1 y 100 pares.".to_string());
    }
    if !opciones.diferencia_tm_maxima.is_finite() || opciones.diferencia_tm_maxima < 0.0 {
        return Err("La diferencia de Tm permitida debe ser cero o positiva.".to_string());
    }
    Ok(())
}

pub fn disenar_cebadores(
    plantilla: &[u8],
    opciones: &OpcionesDiseno,
    genoma: Option<&[u8]>,
) -> Result<Diseno, String> {
    validar_diseno(plantilla, opciones)?;
    let sumas = sumas_acumuladas(plantilla, opciones.repeticion_maxima);
    let mut descartes = Descartes::default();
    let ventanas_izquierdas =
        buscar_ventanas(plantilla, &sumas, Hebra::Directa, opciones, &mut descartes);
    let ventanas_derechas =
        buscar_ventanas(plantilla, &sumas, Hebra::Inversa, opciones, &mut descartes);
    let mut izquierdos =
        construir_cebadores(plantilla, &ventanas_izquierdas, Hebra::Directa, opciones);
    let mut derechos = construir_cebadores(plantilla, &ventanas_derechas, Hebra::Inversa, opciones);
    if let Some(genoma) = genoma {
        anotar_sitios(plantilla, &mut izquierdos, &mut derechos, genoma, opciones);
    }
    let pares = emparejar(&izquierdos, &derechos, opciones)?;
    Ok(Diseno {
        pares,
        izquierdos: izquierdos.len(),
        derechos: derechos.len(),
        descartes,
    })
}
