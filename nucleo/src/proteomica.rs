use rayon::prelude::*;

pub const MASA_DEL_AGUA: f64 = 18.010_564_684;
pub const MASA_DEL_PROTON: f64 = 1.007_276_467;
pub const CARBAMIDOMETILO: f64 = 57.021_464;
pub const CARGAS_SI_NO_SE_SABE: [u8; 3] = [1, 2, 3];
pub const CARGA_MAXIMA: u8 = 8;
pub const MAXIMO_CORTES_PERDIDOS: usize = 3;
pub const LARGO_MAXIMO_PERMITIDO: usize = 100;
pub const FDR_ESTANDAR: f64 = 0.01;
pub const PREFIJO_DE_SENUELO: &str = "señuelo_";

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Digestion {
    pub cortes_perdidos: usize,
    pub largo_minimo: usize,
    pub largo_maximo: usize,
    pub carbamidometilo: bool,
}

impl Default for Digestion {
    fn default() -> Digestion {
        Digestion {
            cortes_perdidos: 1,
            largo_minimo: 6,
            largo_maximo: 40,
            carbamidometilo: true,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Busqueda {
    pub tolerancia_ppm: f64,
    pub tolerancia_fragmento: f64,
}

impl Default for Busqueda {
    fn default() -> Busqueda {
        Busqueda {
            tolerancia_ppm: 10.0,
            tolerancia_fragmento: 0.02,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PeptidoDigerido {
    pub inicio: usize,
    pub fin: usize,
    pub cortes_perdidos: usize,
    pub masa: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Proteina {
    pub nombre: String,
    pub secuencia: Vec<u8>,
    pub senuelo: bool,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PeptidoIndexado {
    pub masa: f64,
    pub proteina: u32,
    pub inicio: u32,
    pub fin: u32,
    pub cortes_perdidos: u8,
}

#[derive(Debug, Clone)]
pub struct Indice {
    pub proteinas: Vec<Proteina>,
    pub peptidos: Vec<PeptidoIndexado>,
    pub digestion: Digestion,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Espectro {
    pub titulo: String,
    pub mz_precursor: f64,
    pub cargas: Vec<u8>,
    pub picos: Vec<(f64, f64)>,
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Puntaje {
    pub hiperpuntaje: f64,
    pub iones_b: usize,
    pub iones_y: usize,
    pub intensidad: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Psm {
    pub espectro: usize,
    pub peptido: PeptidoIndexado,
    pub carga: u8,
    pub error_ppm: f64,
    pub puntaje: Puntaje,
    pub senuelo: bool,
    pub candidatos: usize,
    pub q: f64,
}

pub fn masa_de_residuo(letra: u8, carbamidometilo: bool) -> Option<f64> {
    let masa = match letra {
        b'G' => 57.021_464,
        b'A' => 71.037_114,
        b'S' => 87.032_028,
        b'P' => 97.052_764,
        b'V' => 99.068_414,
        b'T' => 101.047_679,
        b'C' => 103.009_185,
        b'I' | b'L' => 113.084_064,
        b'N' => 114.042_927,
        b'D' => 115.026_943,
        b'Q' => 128.058_578,
        b'K' => 128.094_963,
        b'E' => 129.042_593,
        b'M' => 131.040_485,
        b'H' => 137.058_912,
        b'F' => 147.068_414,
        b'R' => 156.101_111,
        b'Y' => 163.063_329,
        b'W' => 186.079_313,
        _ => return None,
    };
    if letra == b'C' && carbamidometilo {
        Some(masa + CARBAMIDOMETILO)
    } else {
        Some(masa)
    }
}

pub fn masa_del_peptido(peptido: &[u8], carbamidometilo: bool) -> Option<f64> {
    let mut total = MASA_DEL_AGUA;
    for &letra in peptido {
        total += masa_de_residuo(letra, carbamidometilo)?;
    }
    Some(total)
}

pub fn mz_con_carga(masa: f64, carga: u8) -> f64 {
    (masa + carga as f64 * MASA_DEL_PROTON) / carga as f64
}

pub fn masa_neutra(mz: f64, carga: u8) -> f64 {
    (mz - MASA_DEL_PROTON) * carga as f64
}

pub fn revisar_digestion(digestion: &Digestion) -> Result<(), String> {
    if digestion.cortes_perdidos > MAXIMO_CORTES_PERDIDOS {
        return Err(format!(
            "Se admiten como mucho {MAXIMO_CORTES_PERDIDOS} cortes perdidos."
        ));
    }
    if digestion.largo_minimo == 0
        || digestion.largo_minimo > digestion.largo_maximo
        || digestion.largo_maximo > LARGO_MAXIMO_PERMITIDO
    {
        return Err(format!(
            "El largo de los péptidos debe ir de 1 a {LARGO_MAXIMO_PERMITIDO} y el mínimo no puede pasar del máximo."
        ));
    }
    Ok(())
}

pub fn corta_despues(proteina: &[u8], posicion: usize) -> bool {
    let letra = proteina[posicion];
    let siguiente = proteina.get(posicion + 1).copied();
    (letra == b'K' || letra == b'R') && siguiente != Some(b'P')
}

pub fn fronteras(proteina: &[u8]) -> Vec<usize> {
    let mut fronteras = vec![0];
    for posicion in 0..proteina.len() {
        if corta_despues(proteina, posicion) && posicion + 1 < proteina.len() {
            fronteras.push(posicion + 1);
        }
    }
    if !proteina.is_empty() {
        fronteras.push(proteina.len());
    }
    fronteras
}

pub fn digerir(proteina: &[u8], digestion: &Digestion) -> Vec<PeptidoDigerido> {
    let fronteras = fronteras(proteina);
    let mut peptidos = Vec::new();
    for primera in 0..fronteras.len().saturating_sub(1) {
        for cortes_perdidos in 0..=digestion.cortes_perdidos {
            let ultima = primera + cortes_perdidos + 1;
            if ultima >= fronteras.len() {
                break;
            }
            let inicio = fronteras[primera];
            let fin = fronteras[ultima];
            let largo = fin - inicio;
            if largo > digestion.largo_maximo {
                break;
            }
            if largo < digestion.largo_minimo {
                continue;
            }
            let peptido = &proteina[inicio..fin];
            if let Some(masa) = masa_del_peptido(peptido, digestion.carbamidometilo) {
                peptidos.push(PeptidoDigerido {
                    inicio,
                    fin,
                    cortes_perdidos,
                    masa,
                });
            }
        }
    }
    peptidos
}

pub fn senuelo(proteina: &[u8]) -> Vec<u8> {
    proteina.iter().rev().copied().collect()
}

fn con_senuelos(objetivos: &[(String, Vec<u8>)]) -> Vec<Proteina> {
    let mut proteinas = Vec::with_capacity(objetivos.len() * 2);
    for (nombre, secuencia) in objetivos {
        proteinas.push(Proteina {
            nombre: nombre.clone(),
            secuencia: secuencia.clone(),
            senuelo: false,
        });
    }
    for (nombre, secuencia) in objetivos {
        proteinas.push(Proteina {
            nombre: format!("{PREFIJO_DE_SENUELO}{nombre}"),
            secuencia: senuelo(secuencia),
            senuelo: true,
        });
    }
    proteinas
}

fn indexar_proteina(
    numero: usize,
    proteina: &Proteina,
    digestion: &Digestion,
) -> Vec<PeptidoIndexado> {
    digerir(&proteina.secuencia, digestion)
        .into_iter()
        .map(|peptido| PeptidoIndexado {
            masa: peptido.masa,
            proteina: numero as u32,
            inicio: peptido.inicio as u32,
            fin: peptido.fin as u32,
            cortes_perdidos: peptido.cortes_perdidos as u8,
        })
        .collect()
}

impl Indice {
    pub fn construir(
        objetivos: &[(String, Vec<u8>)],
        digestion: Digestion,
    ) -> Result<Indice, String> {
        revisar_digestion(&digestion)?;
        if objetivos.is_empty() {
            return Err("No hay proteínas para digerir.".to_string());
        }
        let proteinas = con_senuelos(objetivos);
        let mut peptidos: Vec<PeptidoIndexado> = proteinas
            .par_iter()
            .enumerate()
            .flat_map_iter(|(numero, proteina)| indexar_proteina(numero, proteina, &digestion))
            .collect();
        peptidos.par_sort_unstable_by(|a, b| {
            a.masa
                .total_cmp(&b.masa)
                .then(a.proteina.cmp(&b.proteina))
                .then(a.inicio.cmp(&b.inicio))
                .then(a.fin.cmp(&b.fin))
        });
        Ok(Indice {
            proteinas,
            peptidos,
            digestion,
        })
    }

    pub fn secuencia(&self, peptido: &PeptidoIndexado) -> &[u8] {
        let proteina = &self.proteinas[peptido.proteina as usize];
        &proteina.secuencia[peptido.inicio as usize..peptido.fin as usize]
    }

    pub fn es_senuelo(&self, peptido: &PeptidoIndexado) -> bool {
        self.proteinas[peptido.proteina as usize].senuelo
    }

    pub fn en_ventana(&self, minima: f64, maxima: f64) -> &[PeptidoIndexado] {
        let desde = self
            .peptidos
            .partition_point(|peptido| peptido.masa < minima);
        let hasta = self
            .peptidos
            .partition_point(|peptido| peptido.masa <= maxima);
        &self.peptidos[desde..hasta.max(desde)]
    }
}

struct Borrador {
    titulo: String,
    mz_precursor: Option<f64>,
    cargas: Vec<u8>,
    picos: Vec<(f64, f64)>,
}

pub fn leer_mgf(texto: &str) -> Result<Vec<Espectro>, String> {
    let mut espectros = Vec::new();
    let mut abierto: Option<Borrador> = None;
    for (numero, linea) in texto.lines().enumerate() {
        let linea = linea.trim();
        if linea.is_empty() || linea.starts_with(['#', ';', '!', '/']) {
            continue;
        }
        if linea.eq_ignore_ascii_case("BEGIN IONS") {
            if abierto.is_some() {
                return Err(format!(
                    "Línea {}: empieza un espectro sin cerrar el anterior con END IONS.",
                    numero + 1
                ));
            }
            abierto = Some(Borrador {
                titulo: format!("espectro {}", espectros.len() + 1),
                mz_precursor: None,
                cargas: Vec::new(),
                picos: Vec::new(),
            });
            continue;
        }
        let Some(borrador) = abierto.as_mut() else {
            continue;
        };
        if linea.eq_ignore_ascii_case("END IONS") {
            if let Some(borrador) = abierto.take() {
                espectros.push(cerrar_espectro(borrador)?);
            }
            continue;
        }
        leer_linea_de_espectro(borrador, linea, numero + 1)?;
    }
    if abierto.is_some() {
        return Err("El último espectro no termina con END IONS.".to_string());
    }
    if espectros.is_empty() {
        return Err("No encontré ningún espectro (BEGIN IONS … END IONS).".to_string());
    }
    Ok(espectros)
}

fn leer_linea_de_espectro(
    borrador: &mut Borrador,
    linea: &str,
    numero: usize,
) -> Result<(), String> {
    if let Some((clave, valor)) = linea.split_once('=')
        && clave.starts_with(|letra: char| letra.is_ascii_alphabetic())
    {
        match clave.trim().to_ascii_uppercase().as_str() {
            "TITLE" => borrador.titulo = valor.trim().to_string(),
            "PEPMASS" => borrador.mz_precursor = Some(leer_pepmass(valor, numero)?),
            "CHARGE" => borrador.cargas = leer_cargas(valor, numero)?,
            _ => {}
        }
        return Ok(());
    }
    let mut piezas = linea.split_whitespace();
    let mz = piezas.next().and_then(|pieza| pieza.parse::<f64>().ok());
    let intensidad = match piezas.next() {
        Some(pieza) => pieza.parse::<f64>().ok(),
        None => Some(1.0),
    };
    match (mz, intensidad) {
        (Some(mz), Some(intensidad)) if mz.is_finite() && intensidad.is_finite() => {
            if mz > 0.0 && intensidad > 0.0 {
                borrador.picos.push((mz, intensidad));
            }
            Ok(())
        }
        _ => Err(format!(
            "Línea {numero}: «{linea}» no es un pico (m/z e intensidad)."
        )),
    }
}

fn leer_pepmass(valor: &str, numero: usize) -> Result<f64, String> {
    let mz = valor
        .split_whitespace()
        .next()
        .and_then(|pieza| pieza.parse::<f64>().ok())
        .filter(|mz| mz.is_finite() && *mz > MASA_DEL_PROTON);
    mz.ok_or_else(|| format!("Línea {numero}: PEPMASS debe ser un m/z positivo."))
}

fn leer_cargas(valor: &str, numero: usize) -> Result<Vec<u8>, String> {
    let mut cargas = Vec::new();
    for pieza in valor.split(|letra: char| letra.is_whitespace() || letra == ',') {
        let pieza = pieza.trim();
        if pieza.is_empty() || pieza.eq_ignore_ascii_case("and") {
            continue;
        }
        if pieza.ends_with('-') {
            return Err(format!(
                "Línea {numero}: solo se buscan iones positivos (CHARGE={valor})."
            ));
        }
        let carga = pieza
            .trim_end_matches('+')
            .parse::<u8>()
            .ok()
            .filter(|carga| (1..=CARGA_MAXIMA).contains(carga))
            .ok_or_else(|| format!("Línea {numero}: la carga debe ir de 1+ a {CARGA_MAXIMA}+."))?;
        if !cargas.contains(&carga) {
            cargas.push(carga);
        }
    }
    Ok(cargas)
}

fn cerrar_espectro(borrador: Borrador) -> Result<Espectro, String> {
    let mz_precursor = borrador
        .mz_precursor
        .ok_or_else(|| format!("El espectro «{}» no tiene PEPMASS.", borrador.titulo))?;
    let cargas = if borrador.cargas.is_empty() {
        CARGAS_SI_NO_SE_SABE.to_vec()
    } else {
        borrador.cargas
    };
    let mut picos = borrador.picos;
    picos.sort_by(|a, b| a.0.total_cmp(&b.0));
    Ok(Espectro {
        titulo: borrador.titulo,
        mz_precursor,
        cargas,
        picos,
    })
}

pub fn normalizar(picos: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let maxima = picos.iter().map(|pico| pico.1).fold(0.0, f64::max);
    if maxima <= 0.0 {
        return picos.to_vec();
    }
    picos
        .iter()
        .map(|&(mz, intensidad)| (mz, intensidad * 100.0 / maxima))
        .collect()
}

pub fn llenar_iones(peptido: &[u8], carbamidometilo: bool, b: &mut Vec<f64>, y: &mut Vec<f64>) {
    b.clear();
    y.clear();
    let mut prefijo = MASA_DEL_PROTON;
    for &letra in &peptido[..peptido.len().saturating_sub(1)] {
        prefijo += masa_de_residuo(letra, carbamidometilo).unwrap_or(0.0);
        b.push(prefijo);
    }
    let mut sufijo = MASA_DEL_AGUA + MASA_DEL_PROTON;
    for &letra in peptido[1.min(peptido.len())..].iter().rev() {
        sufijo += masa_de_residuo(letra, carbamidometilo).unwrap_or(0.0);
        y.push(sufijo);
    }
}

pub fn pico_cercano(picos: &[(f64, f64)], mz: f64, tolerancia: f64) -> Option<usize> {
    let derecha = picos.partition_point(|pico| pico.0 < mz);
    let mut mejor: Option<(usize, f64)> = None;
    for indice in [derecha.wrapping_sub(1), derecha] {
        let Some(pico) = picos.get(indice) else {
            continue;
        };
        let distancia = (pico.0 - mz).abs();
        if distancia <= tolerancia && mejor.is_none_or(|(_, otra)| distancia < otra) {
            mejor = Some((indice, distancia));
        }
    }
    mejor.map(|(indice, _)| indice)
}

fn emparejar(iones: &[f64], picos: &[(f64, f64)], tolerancia: f64) -> (usize, f64) {
    let mut emparejados = 0;
    let mut intensidad = 0.0;
    for &ion in iones {
        if let Some(indice) = pico_cercano(picos, ion, tolerancia) {
            emparejados += 1;
            intensidad += picos[indice].1;
        }
    }
    (emparejados, intensidad)
}

pub fn logaritmo_factorial(n: usize) -> f64 {
    (2..=n).map(|k| (k as f64).ln()).sum()
}

pub fn hiperpuntaje(intensidad: f64, iones_b: usize, iones_y: usize) -> f64 {
    if iones_b + iones_y == 0 {
        return 0.0;
    }
    (1.0 + intensidad).ln() + logaritmo_factorial(iones_b) + logaritmo_factorial(iones_y)
}

pub fn puntuar(b: &[f64], y: &[f64], picos: &[(f64, f64)], tolerancia: f64) -> Puntaje {
    let (iones_b, intensidad_b) = emparejar(b, picos, tolerancia);
    let (iones_y, intensidad_y) = emparejar(y, picos, tolerancia);
    let intensidad = intensidad_b + intensidad_y;
    Puntaje {
        hiperpuntaje: hiperpuntaje(intensidad, iones_b, iones_y),
        iones_b,
        iones_y,
        intensidad,
    }
}

pub fn ventana(masa: f64, tolerancia_ppm: f64) -> (f64, f64) {
    let margen = masa * tolerancia_ppm / 1e6;
    (masa - margen, masa + margen)
}

pub fn error_ppm(observada: f64, teorica: f64) -> f64 {
    (observada - teorica) / teorica * 1e6
}

pub fn mejor_psm(
    numero: usize,
    espectro: &Espectro,
    indice: &Indice,
    busqueda: &Busqueda,
) -> Option<Psm> {
    let picos = normalizar(&espectro.picos);
    let carbamidometilo = indice.digestion.carbamidometilo;
    let mut b = Vec::new();
    let mut y = Vec::new();
    let mut mejor: Option<Psm> = None;
    let mut candidatos = 0;
    for &carga in &espectro.cargas {
        let masa = masa_neutra(espectro.mz_precursor, carga);
        let (minima, maxima) = ventana(masa, busqueda.tolerancia_ppm);
        for peptido in indice.en_ventana(minima, maxima) {
            candidatos += 1;
            llenar_iones(indice.secuencia(peptido), carbamidometilo, &mut b, &mut y);
            let puntaje = puntuar(&b, &y, &picos, busqueda.tolerancia_fragmento);
            let record = mejor.map_or(0.0, |psm| psm.puntaje.hiperpuntaje);
            if puntaje.hiperpuntaje > record {
                mejor = Some(Psm {
                    espectro: numero,
                    peptido: *peptido,
                    carga,
                    error_ppm: error_ppm(masa, peptido.masa),
                    puntaje,
                    senuelo: indice.es_senuelo(peptido),
                    candidatos: 0,
                    q: 1.0,
                });
            }
        }
    }
    mejor.map(|psm| Psm { candidatos, ..psm })
}

pub fn tasa_de_falsos(senuelos: usize, objetivos: usize) -> f64 {
    if objetivos == 0 {
        return 1.0;
    }
    (senuelos as f64 / objetivos as f64).min(1.0)
}

pub fn ordenar_por_puntaje(psms: &mut [Psm]) {
    psms.sort_by(|a, b| {
        b.puntaje
            .hiperpuntaje
            .total_cmp(&a.puntaje.hiperpuntaje)
            .then(a.espectro.cmp(&b.espectro))
    });
}

pub fn tasas_por_umbral(psms: &[Psm]) -> Vec<f64> {
    let mut tasas = vec![1.0; psms.len()];
    let mut objetivos = 0;
    let mut senuelos = 0;
    let mut inicio = 0;
    while inicio < psms.len() {
        let umbral = psms[inicio].puntaje.hiperpuntaje;
        let mut fin = inicio;
        while fin < psms.len() && psms[fin].puntaje.hiperpuntaje == umbral {
            if psms[fin].senuelo {
                senuelos += 1;
            } else {
                objetivos += 1;
            }
            fin += 1;
        }
        let tasa = tasa_de_falsos(senuelos, objetivos);
        for valor in &mut tasas[inicio..fin] {
            *valor = tasa;
        }
        inicio = fin;
    }
    tasas
}

pub fn asignar_valores_q(psms: &mut [Psm]) {
    ordenar_por_puntaje(psms);
    let tasas = tasas_por_umbral(psms);
    let mut minimo = 1.0_f64;
    for (psm, tasa) in psms.iter_mut().zip(&tasas).rev() {
        minimo = minimo.min(*tasa);
        psm.q = minimo;
    }
}

pub fn revisar_busqueda(busqueda: &Busqueda) -> Result<(), String> {
    if !(busqueda.tolerancia_ppm > 0.0 && busqueda.tolerancia_ppm <= 1000.0) {
        return Err("La tolerancia del precursor debe ir de más de 0 a 1000 ppm.".to_string());
    }
    if !(busqueda.tolerancia_fragmento > 0.0 && busqueda.tolerancia_fragmento <= 1.0) {
        return Err("La tolerancia de los fragmentos debe ir de más de 0 a 1 Da.".to_string());
    }
    Ok(())
}

pub fn buscar_espectros(
    espectros: &[Espectro],
    indice: &Indice,
    busqueda: &Busqueda,
) -> Result<Vec<Psm>, String> {
    revisar_busqueda(busqueda)?;
    let mut psms: Vec<Psm> = espectros
        .par_iter()
        .enumerate()
        .filter_map(|(numero, espectro)| mejor_psm(numero, espectro, indice, busqueda))
        .collect();
    asignar_valores_q(&mut psms);
    Ok(psms)
}

pub fn aceptados(psms: &[Psm], umbral: f64) -> Vec<Psm> {
    psms.iter()
        .filter(|psm| !psm.senuelo && psm.q <= umbral)
        .copied()
        .collect()
}
