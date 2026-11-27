use crate::azar::Azar;
use crate::ecoli::{self, Gen};
use crate::formatos::{Lectura, RegistroSam};
use crate::secuencia::{aminoacido, complemento, complemento_inverso};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Mutacion {
    pub posicion: usize,
    pub referencia: Vec<u8>,
    pub alternativa: Vec<u8>,
    pub nombre: String,
}

#[derive(Debug, Clone)]
pub struct OpcionesSimulacion {
    pub longitud_lectura: usize,
    pub cobertura: f64,
    pub calidad_inicial: f64,
    pub calidad_final: f64,
    pub fraccion_con_adaptador: f64,
    pub adaptador: Vec<u8>,
    pub regiones: Vec<(usize, usize)>,
    pub semilla: u64,
}

impl Default for OpcionesSimulacion {
    fn default() -> OpcionesSimulacion {
        OpcionesSimulacion {
            longitud_lectura: 150,
            cobertura: 15.0,
            calidad_inicial: 38.0,
            calidad_final: 24.0,
            fraccion_con_adaptador: 0.08,
            adaptador: b"ACGGTAGCTACGATCGGTACTTGCA".to_vec(),
            regiones: Vec::new(),
            semilla: 1,
        }
    }
}

pub struct Simulacion {
    pub lecturas: Vec<Lectura>,
    pub verdad: Vec<RegistroSam>,
    pub mutaciones: Vec<Mutacion>,
}

pub struct GenomaMutado {
    pub bases: Vec<u8>,
    pub origen: Vec<Option<usize>>,
}

pub fn mutacion_en_codon(
    elegido: &Gen,
    numero_codon: usize,
    aminoacido_nuevo: u8,
) -> Option<Mutacion> {
    if numero_codon == 0 {
        return None;
    }
    let inicio_en_gen = (numero_codon - 1) * 3;
    if inicio_en_gen + 3 > elegido.fin + 1 - elegido.inicio {
        return None;
    }
    let posicion_genoma = |indice_en_gen: usize| -> usize {
        if elegido.hebra == b'-' {
            elegido.fin - 1 - indice_en_gen
        } else {
            elegido.inicio - 1 + indice_en_gen
        }
    };
    let genoma = ecoli::genoma();
    let mut codon = [0u8; 3];
    for (indice, base) in codon.iter_mut().enumerate() {
        let original = genoma[posicion_genoma(inicio_en_gen + indice)];
        *base = if elegido.hebra == b'-' {
            complemento(original)
        } else {
            original
        };
    }
    let original = aminoacido(&codon);
    for indice in 0..3 {
        for nueva in *b"ACGT" {
            if nueva == codon[indice] {
                continue;
            }
            let mut cambiado = codon;
            cambiado[indice] = nueva;
            if aminoacido(&cambiado) == aminoacido_nuevo {
                let posicion = posicion_genoma(inicio_en_gen + indice);
                let alternativa = if elegido.hebra == b'-' {
                    complemento(nueva)
                } else {
                    nueva
                };
                return Some(Mutacion {
                    posicion,
                    referencia: vec![genoma[posicion]],
                    alternativa: vec![alternativa],
                    nombre: format!(
                        "{} {}{}{}",
                        elegido.nombre, original as char, numero_codon, aminoacido_nuevo as char
                    ),
                });
            }
        }
    }
    None
}

pub fn mutaciones_al_azar(
    referencia: &[u8],
    regiones: &[(usize, usize)],
    cantidad_snps: usize,
    cantidad_indels: usize,
    ocupadas: &[usize],
    azar: &mut Azar,
) -> Vec<Mutacion> {
    let regiones: Vec<(usize, usize)> = if regiones.is_empty() {
        vec![(0, referencia.len())]
    } else {
        regiones.to_vec()
    };
    let total: usize = regiones.iter().map(|(inicio, fin)| fin - inicio).sum();
    let mut usadas: Vec<usize> = ocupadas.to_vec();
    let mut mutaciones = Vec::new();
    let mut intentos = 0;
    while mutaciones.len() < cantidad_snps + cantidad_indels && intentos < 100_000 {
        intentos += 1;
        let mut desplazamiento = azar.entero(total);
        let mut posicion = 0;
        for (inicio, fin) in &regiones {
            if desplazamiento < fin - inicio {
                posicion = inicio + desplazamiento;
                break;
            }
            desplazamiento -= fin - inicio;
        }
        if posicion < 5 || posicion + 10 >= referencia.len() {
            continue;
        }
        if usadas.iter().any(|&otra| otra.abs_diff(posicion) < 300) {
            continue;
        }
        let es_indel = mutaciones.len() >= cantidad_snps;
        let mutacion = if es_indel {
            let longitud = 1 + azar.entero(3);
            if azar.entero(2) == 0 {
                Mutacion {
                    posicion: posicion - 1,
                    referencia: referencia[posicion - 1..posicion + longitud].to_vec(),
                    alternativa: vec![referencia[posicion - 1]],
                    nombre: format!("deleción de {longitud} pb"),
                }
            } else {
                let mut alternativa = vec![referencia[posicion - 1]];
                for _ in 0..longitud {
                    alternativa.push(b"ACGT"[azar.entero(4)]);
                }
                Mutacion {
                    posicion: posicion - 1,
                    referencia: vec![referencia[posicion - 1]],
                    alternativa,
                    nombre: format!("inserción de {longitud} pb"),
                }
            }
        } else {
            let original = referencia[posicion];
            let mut nueva = original;
            while nueva == original {
                nueva = b"ACGT"[azar.entero(4)];
            }
            Mutacion {
                posicion,
                referencia: vec![original],
                alternativa: vec![nueva],
                nombre: "SNP".to_string(),
            }
        };
        usadas.push(posicion);
        mutaciones.push(mutacion);
    }
    mutaciones.sort_by_key(|mutacion| mutacion.posicion);
    mutaciones
}

pub fn aplicar_mutaciones(referencia: &[u8], mutaciones: &[Mutacion]) -> GenomaMutado {
    let mut bases = Vec::with_capacity(referencia.len() + 1000);
    let mut origen = Vec::with_capacity(referencia.len() + 1000);
    let mut posicion = 0;
    for mutacion in mutaciones {
        while posicion < mutacion.posicion {
            bases.push(referencia[posicion]);
            origen.push(Some(posicion));
            posicion += 1;
        }
        let comunes = mutacion.referencia.len().min(mutacion.alternativa.len());
        for (indice, &alternativa) in mutacion.alternativa[..comunes].iter().enumerate() {
            bases.push(alternativa);
            origen.push(Some(posicion + indice));
        }
        for &insertada in &mutacion.alternativa[comunes..] {
            bases.push(insertada);
            origen.push(None);
        }
        posicion += mutacion.referencia.len();
    }
    while posicion < referencia.len() {
        bases.push(referencia[posicion]);
        origen.push(Some(posicion));
        posicion += 1;
    }
    GenomaMutado { bases, origen }
}

pub fn calidades(opciones: &OpcionesSimulacion, longitud: usize, azar: &mut Azar) -> Vec<u8> {
    let mut calidad = Vec::with_capacity(longitud);
    for indice in 0..longitud {
        let avance = indice as f64 / longitud.max(2).saturating_sub(1) as f64;
        let media =
            opciones.calidad_inicial + (opciones.calidad_final - opciones.calidad_inicial) * avance;
        let valor = (media + azar.normal() * 2.5).round().clamp(2.0, 41.0);
        calidad.push(valor as u8);
    }
    calidad
}

pub fn agregar_errores(bases: &mut [u8], calidad: &[u8], azar: &mut Azar) {
    for (base, &q) in bases.iter_mut().zip(calidad) {
        let probabilidad_error = 10f64.powf(-(q as f64) / 10.0);
        if azar.decimal() < probabilidad_error {
            let original = *base;
            while *base == original {
                *base = b"ACGT"[azar.entero(4)];
            }
        }
    }
}

pub fn cigar_desde_origen(origen: &[Option<usize>]) -> (usize, Vec<(u32, u8)>) {
    let mut operaciones: Vec<u8> = Vec::new();
    let mut primera = None;
    let mut anterior: Option<usize> = None;
    for posicion in origen {
        match posicion {
            Some(actual) => {
                if primera.is_none() {
                    primera = Some(*actual);
                }
                if let Some(previa) = anterior {
                    let borradas = actual.saturating_sub(previa + 1);
                    operaciones.extend(std::iter::repeat_n(b'D', borradas));
                }
                operaciones.push(b'M');
                anterior = Some(*actual);
            }
            None => {
                if primera.is_some() {
                    operaciones.push(b'I');
                } else {
                    operaciones.push(b'S');
                }
            }
        }
    }
    for operacion in operaciones.iter_mut().rev() {
        if *operacion != b'I' {
            break;
        }
        *operacion = b'S';
    }
    let mut cigar: Vec<(u32, u8)> = Vec::new();
    for operacion in operaciones {
        match cigar.last_mut() {
            Some((longitud, ultima)) if *ultima == operacion => *longitud += 1,
            _ => cigar.push((1, operacion)),
        }
    }
    (primera.unwrap_or(0), cigar)
}

pub fn simular_lecturas(
    referencia: &[u8],
    nombre_referencia: &str,
    mutaciones: &[Mutacion],
    opciones: &OpcionesSimulacion,
) -> Result<Simulacion, String> {
    let longitud = opciones.longitud_lectura;
    if !(30..=1000).contains(&longitud) {
        return Err("La longitud de lectura debe estar entre 30 y 1000 bases.".to_string());
    }
    if !(0.1..=200.0).contains(&opciones.cobertura) {
        return Err("La cobertura debe estar entre 0.1 y 200.".to_string());
    }
    if referencia.len() < longitud * 2 {
        return Err(format!(
            "La referencia debe medir al menos {} bases para lecturas de {longitud}.",
            longitud * 2
        ));
    }
    let mut ordenadas = mutaciones.to_vec();
    ordenadas.sort_by_key(|mutacion| mutacion.posicion);
    for par in ordenadas.windows(2) {
        if par[0].posicion + par[0].referencia.len() > par[1].posicion {
            return Err("Hay mutaciones que se superponen.".to_string());
        }
    }
    let genoma = aplicar_mutaciones(referencia, &ordenadas);
    let regiones = regiones_en_mutado(&genoma, &opciones.regiones, referencia.len());
    let largo_total: usize = regiones.iter().map(|(inicio, fin)| fin - inicio).sum();
    if largo_total < longitud {
        return Err("Las regiones son más cortas que una lectura.".to_string());
    }
    let cantidad = ((largo_total as f64 * opciones.cobertura) / longitud as f64).round() as usize;
    if cantidad * longitud > 150_000_000 {
        return Err(
            "La simulación pasaría de 150 millones de bases; baja la cobertura.".to_string(),
        );
    }
    let mut azar = Azar::nuevo(opciones.semilla);
    let mut lecturas = Vec::with_capacity(cantidad);
    let mut verdad = Vec::with_capacity(cantidad);
    for numero in 0..cantidad {
        let con_adaptador = azar.decimal() < opciones.fraccion_con_adaptador;
        let largo_fragmento = if con_adaptador {
            longitud / 3 + azar.entero(longitud - longitud / 3 - 10)
        } else {
            longitud
        };
        let (inicio_region, fin_region) = elegir_region(&regiones, &mut azar);
        if fin_region - inicio_region < largo_fragmento {
            continue;
        }
        let inicio = inicio_region + azar.entero(fin_region - inicio_region - largo_fragmento + 1);
        let fragmento = &genoma.bases[inicio..inicio + largo_fragmento];
        let origen = &genoma.origen[inicio..inicio + largo_fragmento];
        let reversa = azar.entero(2) == 1;
        let mut bases = if reversa {
            complemento_inverso(fragmento)
        } else {
            fragmento.to_vec()
        };
        let mut origen_lectura: Vec<Option<usize>> = if reversa {
            origen.iter().rev().copied().collect()
        } else {
            origen.to_vec()
        };
        let mut relleno = 0;
        while bases.len() < longitud {
            let siguiente = if relleno < opciones.adaptador.len() {
                opciones.adaptador[relleno]
            } else {
                b"ACGT"[azar.entero(4)]
            };
            bases.push(siguiente);
            origen_lectura.push(None);
            relleno += 1;
        }
        let calidad_numerica = calidades(opciones, longitud, &mut azar);
        agregar_errores(&mut bases, &calidad_numerica, &mut azar);
        let calidad: Vec<u8> = calidad_numerica.iter().map(|q| q + 33).collect();
        let nombre = format!(
            "lectura{}_{}_{}",
            numero + 1,
            origen.iter().flatten().next().map(|p| p + 1).unwrap_or(0),
            if reversa { "menos" } else { "mas" }
        )
        .into_bytes();
        let (secuencia_sam, calidad_sam, origen_sam) = if reversa {
            (
                complemento_inverso(&bases),
                calidad.iter().rev().copied().collect(),
                origen_lectura.iter().rev().copied().collect::<Vec<_>>(),
            )
        } else {
            (bases.clone(), calidad.clone(), origen_lectura.clone())
        };
        let (posicion, cigar) = cigar_desde_origen(&origen_sam);
        verdad.push(RegistroSam {
            nombre: nombre.clone(),
            bandera: if reversa { 16 } else { 0 },
            referencia: nombre_referencia.to_string(),
            posicion: posicion + 1,
            calidad_mapeo: 60,
            cigar,
            secuencia: secuencia_sam,
            calidad: calidad_sam,
        });
        lecturas.push(Lectura {
            nombre,
            secuencia: bases,
            calidad,
        });
    }
    Ok(Simulacion {
        lecturas,
        verdad,
        mutaciones: ordenadas,
    })
}

pub fn regiones_en_mutado(
    genoma: &GenomaMutado,
    regiones: &[(usize, usize)],
    longitud_referencia: usize,
) -> Vec<(usize, usize)> {
    if regiones.is_empty() {
        return vec![(0, genoma.bases.len())];
    }
    let mut convertidas = Vec::new();
    for &(inicio, fin) in regiones {
        let fin = fin.min(longitud_referencia);
        if inicio >= fin {
            continue;
        }
        let primera = genoma
            .origen
            .iter()
            .position(|origen| origen.is_some_and(|p| p >= inicio));
        let ultima = genoma
            .origen
            .iter()
            .rposition(|origen| origen.is_some_and(|p| p < fin));
        if let (Some(primera), Some(ultima)) = (primera, ultima)
            && primera <= ultima
        {
            convertidas.push((primera, ultima + 1));
        }
    }
    convertidas
}

pub fn elegir_region(regiones: &[(usize, usize)], azar: &mut Azar) -> (usize, usize) {
    let pesos: Vec<f64> = regiones
        .iter()
        .map(|(inicio, fin)| (fin - inicio) as f64)
        .collect();
    regiones[azar.elegir_con_pesos(&pesos)]
}

pub const REGIONES_DE_EJEMPLO: [(usize, usize); 2] =
    [(2_300_000, 2_400_000), (3_900_000, 4_200_000)];

pub fn mutaciones_de_resistencia() -> Result<Vec<Mutacion>, String> {
    let mut mutaciones = Vec::new();
    for (nombre, codon, nuevo) in [("rpoB", 531, b'F'), ("gyrA", 83, b'L')] {
        let elegido = ecoli::buscar_gen(nombre)
            .ok_or_else(|| format!("No encuentro el gen {nombre} en E. coli."))?;
        let mutacion = mutacion_en_codon(elegido, codon, nuevo)
            .ok_or_else(|| format!("No se puede crear la mutación en {nombre}."))?;
        mutaciones.push(mutacion);
    }
    Ok(mutaciones)
}

pub fn muestra_de_ejemplo() -> Result<Simulacion, String> {
    let conocidas = mutaciones_de_resistencia()?;
    let ocupadas: Vec<usize> = conocidas.iter().map(|mutacion| mutacion.posicion).collect();
    let mut azar = Azar::nuevo(2026);
    let mut todas = mutaciones_al_azar(
        ecoli::genoma(),
        &REGIONES_DE_EJEMPLO,
        20,
        4,
        &ocupadas,
        &mut azar,
    );
    todas.extend(conocidas);
    let opciones = OpcionesSimulacion {
        regiones: REGIONES_DE_EJEMPLO.to_vec(),
        semilla: 2026,
        ..OpcionesSimulacion::default()
    };
    simular_lecturas(ecoli::genoma(), ecoli::IDENTIFICADOR, &todas, &opciones)
}
