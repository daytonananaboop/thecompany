use crate::secuencia::{complemento_inverso, contenido_gc, indice_base};
use std::collections::HashMap;
use std::io::{self, Write};

pub const LONGITUD_MINIMA: usize = 90;
pub const LONGITUD_PARA_ENTRENAR: usize = 300;
pub const SECUENCIA_MINIMA_PARA_ENTRENAR: usize = 20_000;
pub const SOLAPAMIENTO_MISMA_HEBRA: usize = 60;
pub const SOLAPAMIENTO_ENTRE_COLAS: usize = 200;
pub const SOLAPAMIENTO_ENTRE_CABEZAS: usize = 0;
pub const CODONES_DE_INICIO: &[&str] = &["ATG", "GTG", "TTG"];
pub const MOTIVOS_RBS: &[&str] = &[
    "AGGAGG", "AGGAG", "GGAGG", "AGGA", "GGAG", "GAGG", "AGG", "GGA", "GAG",
];
pub const ESPACIADORES_RBS: &[(usize, usize)] = &[(5, 10), (11, 15)];
pub const CLASES_RBS: usize = 9;
pub const RONDAS_DE_ENTRENAMIENTO: usize = 3;
pub const PESO_DEL_INICIO: f64 = 3.0;
const HEXAMEROS: usize = 4096;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Rbs {
    pub motivo: &'static str,
    pub espaciador: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Inicio {
    pub posicion: usize,
    pub codon: &'static str,
    pub rbs: Option<Rbs>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Marco {
    pub hebra: u8,
    pub inicios: Vec<Inicio>,
    pub fin: usize,
    pub tiene_parada: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct GenPredicho {
    pub desde: usize,
    pub hasta: usize,
    pub hebra: u8,
    pub codon_de_inicio: &'static str,
    pub rbs: Option<Rbs>,
    pub puntaje_codificacion: f64,
    pub puntaje_inicio: f64,
    pub puntaje: f64,
    pub parcial_al_inicio: bool,
    pub parcial_al_final: bool,
}

impl GenPredicho {
    pub fn longitud(&self) -> usize {
        self.hasta.saturating_sub(self.desde)
    }

    pub fn extremo_de_parada(&self) -> usize {
        if self.hebra == b'-' {
            self.desde
        } else {
            self.hasta
        }
    }

    pub fn completo(&self) -> bool {
        !self.parcial_al_inicio && !self.parcial_al_final
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Modelo {
    pub codificacion: Vec<f64>,
    pub rbs: [f64; CLASES_RBS],
    pub codones_de_inicio: [f64; 3],
    pub inicio_tipico: f64,
    pub gc: f64,
    pub genes_de_entrenamiento: usize,
}

impl Modelo {
    pub fn neutro(gc: f64) -> Modelo {
        Modelo {
            codificacion: vec![0.0; HEXAMEROS],
            rbs: [0.0; CLASES_RBS],
            codones_de_inicio: [0.0; 3],
            inicio_tipico: 0.0,
            gc,
            genes_de_entrenamiento: 0,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Prediccion {
    pub genes: Vec<GenPredicho>,
    pub modelo: Modelo,
    pub entrenado_con_la_entrada: bool,
}

pub struct Hebra {
    pub signo: u8,
    pub bases: Vec<u8>,
    pub hexameros: Vec<Option<u16>>,
}

pub struct Hebras {
    pub directa: Hebra,
    pub inversa: Hebra,
}

impl Hebras {
    pub fn de(&self, signo: u8) -> &Hebra {
        if signo == b'-' {
            &self.inversa
        } else {
            &self.directa
        }
    }

    pub fn longitud(&self) -> usize {
        self.directa.bases.len()
    }
}

pub fn es_codon_de_parada(codon: &[u8]) -> bool {
    matches!(codon, b"TAA" | b"TAG" | b"TGA")
}

pub fn codon_de_inicio(codon: &[u8]) -> Option<&'static str> {
    CODONES_DE_INICIO
        .iter()
        .copied()
        .find(|nombre| nombre.as_bytes() == codon)
}

pub fn buscar_rbs(bases: &[u8], inicio: usize) -> Option<Rbs> {
    for longitud in (3..=6).rev() {
        for &(menor, mayor) in ESPACIADORES_RBS {
            for espaciador in menor..=mayor {
                let Some(desde) = inicio.checked_sub(espaciador + longitud) else {
                    continue;
                };
                let Some(tramo) = bases.get(desde..desde + longitud) else {
                    continue;
                };
                for &motivo in MOTIVOS_RBS {
                    if motivo.as_bytes() == tramo {
                        return Some(Rbs { motivo, espaciador });
                    }
                }
            }
        }
    }
    None
}

pub fn clase_rbs(rbs: Option<Rbs>) -> usize {
    match rbs {
        None => CLASES_RBS - 1,
        Some(encontrado) => {
            let lejano = encontrado.espaciador > ESPACIADORES_RBS[0].1;
            let largo = encontrado.motivo.len().clamp(3, 6);
            (largo - 3) * 2 + usize::from(lejano)
        }
    }
}

pub fn hexameros_por_posicion(bases: &[u8]) -> Vec<Option<u16>> {
    let mut hexameros = vec![None; bases.len()];
    let mut numero = 0;
    let mut seguidas = 0;
    for (posicion, &base) in bases.iter().enumerate() {
        match indice_base(base) {
            Some(indice) => {
                numero = (numero * 4 + indice) % HEXAMEROS;
                seguidas += 1;
            }
            None => seguidas = 0,
        }
        if seguidas >= 6 {
            hexameros[posicion - 5] = Some(numero as u16);
        }
    }
    hexameros
}

pub fn preparar_hebras(adn: &[u8]) -> Hebras {
    let directa = adn.to_ascii_uppercase();
    let inversa = complemento_inverso(&directa);
    Hebras {
        directa: Hebra {
            signo: b'+',
            hexameros: hexameros_por_posicion(&directa),
            bases: directa,
        },
        inversa: Hebra {
            signo: b'-',
            hexameros: hexameros_por_posicion(&inversa),
            bases: inversa,
        },
    }
}

pub fn marcos_abiertos(hebra: &Hebra) -> Vec<Marco> {
    let bases = &hebra.bases;
    let mut marcos = Vec::new();
    let mut inicios = Vec::new();
    for desfase in 0..3 {
        inicios.clear();
        let mut posicion = desfase;
        while posicion + 3 <= bases.len() {
            let codon = &bases[posicion..posicion + 3];
            if es_codon_de_parada(codon) {
                guardar_marco(&mut marcos, hebra.signo, &inicios, posicion + 3, true);
                inicios.clear();
            } else if let Some(nombre) = codon_de_inicio(codon) {
                inicios.push(Inicio {
                    posicion,
                    codon: nombre,
                    rbs: buscar_rbs(bases, posicion),
                });
            } else if posicion == desfase {
                inicios.push(Inicio {
                    posicion,
                    codon: "",
                    rbs: None,
                });
            }
            posicion += 3;
        }
        guardar_marco(&mut marcos, hebra.signo, &inicios, posicion, false);
    }
    marcos
}

fn guardar_marco(
    marcos: &mut Vec<Marco>,
    signo: u8,
    inicios: &[Inicio],
    fin: usize,
    tiene_parada: bool,
) {
    let suficientemente_largos: Vec<Inicio> = inicios
        .iter()
        .copied()
        .take_while(|inicio| fin - inicio.posicion >= LONGITUD_MINIMA)
        .collect();
    if !suficientemente_largos.is_empty() {
        marcos.push(Marco {
            hebra: signo,
            inicios: suficientemente_largos,
            fin,
            tiene_parada,
        });
    }
}

fn reflejar_si_es_inversa(
    signo: u8,
    desde: usize,
    hasta: usize,
    longitud: usize,
) -> (usize, usize) {
    if signo == b'-' {
        (longitud - hasta, longitud - desde)
    } else {
        (desde, hasta)
    }
}

pub fn hexameros_de_fondo(hebras: &Hebras) -> Vec<f64> {
    let mut conteo = vec![0.0; HEXAMEROS];
    for hebra in [&hebras.directa, &hebras.inversa] {
        for numero in hebra.hexameros.iter().flatten() {
            conteo[*numero as usize] += 1.0;
        }
    }
    conteo
}

pub fn hexameros_en_genes(hebras: &Hebras, genes: &[GenPredicho]) -> Vec<f64> {
    let longitud = hebras.longitud();
    let mut conteo = vec![0.0; HEXAMEROS];
    for predicho in genes {
        let hebra = hebras.de(predicho.hebra);
        let (inicio, fin) =
            reflejar_si_es_inversa(predicho.hebra, predicho.desde, predicho.hasta, longitud);
        let mut posicion = inicio;
        while posicion + 6 <= fin {
            if let Some(numero) = hebra.hexameros[posicion] {
                conteo[numero as usize] += 1.0;
            }
            posicion += 3;
        }
    }
    conteo
}

pub fn razon_logaritmica(observados: &[f64], esperados: &[f64]) -> Vec<f64> {
    let total_observados: f64 = observados.iter().sum::<f64>() + observados.len() as f64;
    let total_esperados: f64 = esperados.iter().sum::<f64>() + esperados.len() as f64;
    observados
        .iter()
        .zip(esperados)
        .map(|(observado, esperado)| {
            let frecuencia_observada = (observado + 1.0) / total_observados;
            let frecuencia_esperada = (esperado + 1.0) / total_esperados;
            (frecuencia_observada / frecuencia_esperada).ln()
        })
        .collect()
}

fn indice_de_codon(codon: &str) -> Option<usize> {
    CODONES_DE_INICIO.iter().position(|nombre| *nombre == codon)
}

pub fn entrenar_inicios(modelo: &mut Modelo, marcos: &[Marco], genes: &[GenPredicho]) {
    let mut rbs_de_todos = [0.0; CLASES_RBS];
    let mut codones_de_todos = [0.0; 3];
    for marco in marcos {
        for inicio in &marco.inicios {
            if let Some(indice) = indice_de_codon(inicio.codon) {
                codones_de_todos[indice] += 1.0;
                rbs_de_todos[clase_rbs(inicio.rbs)] += 1.0;
            }
        }
    }
    let mut rbs_de_genes = [0.0; CLASES_RBS];
    let mut codones_de_genes = [0.0; 3];
    let mut con_codon = Vec::new();
    for predicho in genes {
        if let Some(indice) = indice_de_codon(predicho.codon_de_inicio) {
            codones_de_genes[indice] += 1.0;
            rbs_de_genes[clase_rbs(predicho.rbs)] += 1.0;
            con_codon.push(predicho);
        }
    }
    let rbs = razon_logaritmica(&rbs_de_genes, &rbs_de_todos);
    let codones = razon_logaritmica(&codones_de_genes, &codones_de_todos);
    modelo.rbs.copy_from_slice(&rbs);
    modelo.codones_de_inicio.copy_from_slice(&codones);
    let mut suma = 0.0;
    for predicho in &con_codon {
        suma += puntaje_del_inicio(modelo, predicho.codon_de_inicio, predicho.rbs);
    }
    modelo.inicio_tipico = suma / con_codon.len().max(1) as f64;
}

pub fn puntaje_del_inicio(modelo: &Modelo, codon: &str, rbs: Option<Rbs>) -> f64 {
    match indice_de_codon(codon) {
        Some(indice) => modelo.codones_de_inicio[indice] + modelo.rbs[clase_rbs(rbs)],
        None => modelo.inicio_tipico,
    }
}

pub fn candidatos_del_marco(
    marco: &Marco,
    hebra: &Hebra,
    modelo: &Modelo,
    candidatos: &mut Vec<GenPredicho>,
) {
    let longitud = hebra.bases.len();
    let mut posicion = if marco.tiene_parada {
        marco.fin - 3
    } else {
        marco.fin
    };
    let mut codificacion = 0.0;
    for inicio in marco.inicios.iter().rev() {
        while posicion > inicio.posicion {
            posicion -= 3;
            if let Some(numero) = hebra.hexameros[posicion] {
                codificacion += modelo.codificacion[numero as usize];
            }
        }
        let puntaje_inicio = puntaje_del_inicio(modelo, inicio.codon, inicio.rbs);
        let puntaje = codificacion + PESO_DEL_INICIO * puntaje_inicio;
        if codificacion <= 0.0 || puntaje <= 0.0 {
            continue;
        }
        let (desde, hasta) =
            reflejar_si_es_inversa(marco.hebra, inicio.posicion, marco.fin, longitud);
        candidatos.push(GenPredicho {
            desde,
            hasta,
            hebra: marco.hebra,
            codon_de_inicio: inicio.codon,
            rbs: inicio.rbs,
            puntaje_codificacion: codificacion,
            puntaje_inicio,
            puntaje,
            parcial_al_inicio: inicio.codon.is_empty(),
            parcial_al_final: !marco.tiene_parada,
        });
    }
}

pub fn solapamiento_permitido(antes: &GenPredicho, despues: &GenPredicho) -> usize {
    if antes.hebra == despues.hebra {
        SOLAPAMIENTO_MISMA_HEBRA
    } else if antes.hebra == b'+' {
        SOLAPAMIENTO_ENTRE_COLAS
    } else {
        SOLAPAMIENTO_ENTRE_CABEZAS
    }
}

pub fn pueden_ir_seguidos(antes: &GenPredicho, despues: &GenPredicho) -> bool {
    if antes.desde >= despues.desde || antes.hasta >= despues.hasta {
        return false;
    }
    let solapamiento = antes.hasta.saturating_sub(despues.desde);
    solapamiento <= solapamiento_permitido(antes, despues)
}

fn mismo_final(uno: &GenPredicho, otro: &GenPredicho) -> bool {
    uno.hasta == otro.hasta && uno.hebra == otro.hebra
}

fn fin_de_cada_grupo(ordenados: &[&GenPredicho]) -> Vec<usize> {
    let total = ordenados.len();
    let mut fin_del_grupo = vec![total; total];
    for posicion in (0..total.saturating_sub(1)).rev() {
        if mismo_final(ordenados[posicion], ordenados[posicion + 1]) {
            fin_del_grupo[posicion] = fin_del_grupo[posicion + 1];
        } else {
            fin_del_grupo[posicion] = posicion + 1;
        }
    }
    fin_del_grupo
}

pub fn elegir_genes(candidatos: &[GenPredicho]) -> Vec<usize> {
    let mut orden: Vec<usize> = (0..candidatos.len()).collect();
    orden.sort_by_key(|&indice| (candidatos[indice].hasta, candidatos[indice].desde));
    let ordenados: Vec<&GenPredicho> = orden.iter().map(|&indice| &candidatos[indice]).collect();
    let fines: Vec<usize> = ordenados.iter().map(|candidato| candidato.hasta).collect();
    let fin_del_grupo = fin_de_cada_grupo(&ordenados);
    let mayor_solapamiento = [
        SOLAPAMIENTO_MISMA_HEBRA,
        SOLAPAMIENTO_ENTRE_COLAS,
        SOLAPAMIENTO_ENTRE_CABEZAS,
    ]
    .into_iter()
    .max()
    .unwrap_or(0);
    let mut cadena = vec![0.0; orden.len()];
    let mut anterior: Vec<Option<usize>> = vec![None; orden.len()];
    let mut mejor_hasta_aqui: Vec<Option<usize>> = vec![None; orden.len()];
    let mut mejor_del_grupo: Vec<usize> = (0..orden.len()).collect();
    for actual in 0..orden.len() {
        let candidato = ordenados[actual];
        let sin_tocarse = fines[..actual].partition_point(|&fin| fin <= candidato.desde);
        let rozando = fines[..actual]
            .partition_point(|&fin| fin <= candidato.desde.saturating_add(mayor_solapamiento));
        let mut previo = if sin_tocarse > 0 {
            mejor_hasta_aqui[sin_tocarse - 1]
        } else {
            None
        };
        let mut suma = previo.map_or(0.0, |indice| cadena[indice]);
        let mut grupo = sin_tocarse;
        while grupo < rozando {
            let fin = fin_del_grupo[grupo].min(rozando);
            let empiezan_antes =
                ordenados[grupo..fin].partition_point(|otro| otro.desde < candidato.desde);
            if empiezan_antes > 0 {
                let otro = mejor_del_grupo[grupo + empiezan_antes - 1];
                if cadena[otro] > suma && pueden_ir_seguidos(ordenados[otro], candidato) {
                    suma = cadena[otro];
                    previo = Some(otro);
                }
            }
            grupo = fin;
        }
        cadena[actual] = suma + candidato.puntaje;
        anterior[actual] = previo;
        if actual > 0
            && mismo_final(ordenados[actual - 1], candidato)
            && cadena[mejor_del_grupo[actual - 1]] >= cadena[actual]
        {
            mejor_del_grupo[actual] = mejor_del_grupo[actual - 1];
        }
        let mut mejor = if actual > 0 {
            mejor_hasta_aqui[actual - 1]
        } else {
            None
        };
        let supera = match mejor {
            Some(indice) => cadena[actual] > cadena[indice],
            None => cadena[actual] > 0.0,
        };
        if supera {
            mejor = Some(actual);
        }
        mejor_hasta_aqui[actual] = mejor;
    }
    let mut elegidos = Vec::new();
    let mut paso = mejor_hasta_aqui.last().copied().flatten();
    while let Some(indice) = paso {
        elegidos.push(orden[indice]);
        paso = anterior[indice];
    }
    elegidos.reverse();
    elegidos
}

pub fn marcos_de_ambas_hebras(hebras: &Hebras) -> Vec<Marco> {
    let mut marcos = marcos_abiertos(&hebras.directa);
    marcos.extend(marcos_abiertos(&hebras.inversa));
    marcos
}

pub fn genes_con_modelo(hebras: &Hebras, marcos: &[Marco], modelo: &Modelo) -> Vec<GenPredicho> {
    let mut candidatos = Vec::new();
    for marco in marcos {
        candidatos_del_marco(marco, hebras.de(marco.hebra), modelo, &mut candidatos);
    }
    let mut genes: Vec<GenPredicho> = elegir_genes(&candidatos)
        .into_iter()
        .map(|indice| candidatos[indice])
        .collect();
    genes.sort_by_key(|predicho| (predicho.desde, predicho.hasta));
    genes
}

pub fn marcos_mas_largos(hebras: &Hebras, marcos: &[Marco]) -> Vec<GenPredicho> {
    let longitud = hebras.longitud();
    let mut candidatos = Vec::new();
    for marco in marcos {
        let Some(inicio) = marco.inicios.iter().find(|inicio| !inicio.codon.is_empty()) else {
            continue;
        };
        let largo = marco.fin - inicio.posicion;
        if largo < LONGITUD_PARA_ENTRENAR || !marco.tiene_parada {
            continue;
        }
        let (desde, hasta) =
            reflejar_si_es_inversa(marco.hebra, inicio.posicion, marco.fin, longitud);
        candidatos.push(GenPredicho {
            desde,
            hasta,
            hebra: marco.hebra,
            codon_de_inicio: inicio.codon,
            rbs: inicio.rbs,
            puntaje: largo as f64,
            ..GenPredicho::default()
        });
    }
    elegir_genes(&candidatos)
        .into_iter()
        .map(|indice| candidatos[indice])
        .collect()
}

fn entrenar(hebras: &Hebras, marcos: &[Marco]) -> Modelo {
    let fondo = hexameros_de_fondo(hebras);
    let mut modelo = Modelo::neutro(contenido_gc(&hebras.directa.bases));
    let mut entrenamiento = marcos_mas_largos(hebras, marcos);
    for ronda in 1..=RONDAS_DE_ENTRENAMIENTO {
        if ronda > 1 {
            entrenamiento = genes_con_modelo(hebras, marcos, &modelo);
            entrenamiento.retain(|predicho| {
                predicho.completo() && predicho.longitud() >= LONGITUD_PARA_ENTRENAR
            });
        }
        let en_genes = hexameros_en_genes(hebras, &entrenamiento);
        modelo.codificacion = razon_logaritmica(&en_genes, &fondo);
        entrenar_inicios(&mut modelo, marcos, &entrenamiento);
        modelo.genes_de_entrenamiento = entrenamiento.len();
    }
    modelo
}

pub fn entrenar_modelo(adn: &[u8]) -> Modelo {
    let hebras = preparar_hebras(adn);
    let marcos = marcos_de_ambas_hebras(&hebras);
    entrenar(&hebras, &marcos)
}

pub fn predecir_genes(adn: &[u8], respaldo: &Modelo) -> Prediccion {
    let hebras = preparar_hebras(adn);
    let marcos = marcos_de_ambas_hebras(&hebras);
    let entrenado_con_la_entrada = adn.len() >= SECUENCIA_MINIMA_PARA_ENTRENAR;
    let modelo = if entrenado_con_la_entrada {
        entrenar(&hebras, &marcos)
    } else {
        respaldo.clone()
    };
    let genes = genes_con_modelo(&hebras, &marcos, &modelo);
    Prediccion {
        genes,
        modelo,
        entrenado_con_la_entrada,
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Anotado {
    pub desde: usize,
    pub hasta: usize,
    pub hebra: u8,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Comparacion {
    pub anotados: usize,
    pub encontrados: usize,
    pub con_inicio_exacto: usize,
    pub anotados_largos: usize,
    pub largos_encontrados: usize,
    pub largos_con_inicio_exacto: usize,
    pub predichos_completos: usize,
    pub predichos_sin_anotacion: usize,
    pub pareja: Vec<Option<usize>>,
    pub hallado: Vec<bool>,
}

fn parada_anotada(anotado: &Anotado) -> usize {
    if anotado.hebra == b'-' {
        anotado.desde
    } else {
        anotado.hasta
    }
}

pub fn comparar_con_anotacion(predichos: &[GenPredicho], anotados: &[Anotado]) -> Comparacion {
    let mut por_parada: HashMap<(u8, usize), Vec<usize>> = HashMap::new();
    for (indice, anotado) in anotados.iter().enumerate() {
        por_parada
            .entry((anotado.hebra, parada_anotada(anotado)))
            .or_default()
            .push(indice);
    }
    let mut hallado = vec![false; anotados.len()];
    let mut exacto = vec![false; anotados.len()];
    let mut comparacion = Comparacion {
        anotados: anotados.len(),
        ..Comparacion::default()
    };
    for predicho in predichos {
        let con_su_parada: &[usize] =
            match por_parada.get(&(predicho.hebra, predicho.extremo_de_parada())) {
                Some(indices) if !predicho.parcial_al_final => indices,
                _ => &[],
            };
        let mismo_inicio = con_su_parada.iter().copied().find(|&indice| {
            anotados[indice].desde == predicho.desde && anotados[indice].hasta == predicho.hasta
        });
        for &indice in con_su_parada {
            hallado[indice] = true;
        }
        if let Some(indice) = mismo_inicio {
            exacto[indice] = true;
        }
        comparacion
            .pareja
            .push(mismo_inicio.or(con_su_parada.first().copied()));
        if predicho.completo() {
            comparacion.predichos_completos += 1;
            if con_su_parada.is_empty() {
                comparacion.predichos_sin_anotacion += 1;
            }
        }
    }
    for (indice, anotado) in anotados.iter().enumerate() {
        let largo = anotado.hasta.saturating_sub(anotado.desde) >= LONGITUD_PARA_ENTRENAR;
        comparacion.anotados_largos += usize::from(largo);
        comparacion.encontrados += usize::from(hallado[indice]);
        comparacion.con_inicio_exacto += usize::from(exacto[indice]);
        comparacion.largos_encontrados += usize::from(largo && hallado[indice]);
        comparacion.largos_con_inicio_exacto += usize::from(largo && exacto[indice]);
    }
    comparacion.hallado = hallado;
    comparacion
}

pub fn escribir_gff(
    destino: &mut impl Write,
    nombre: &str,
    longitud: usize,
    genes: &[GenPredicho],
) -> io::Result<()> {
    writeln!(destino, "##gff-version 3")?;
    writeln!(destino, "##sequence-region {nombre} 1 {longitud}")?;
    for (numero, predicho) in genes.iter().enumerate() {
        let (motivo, espaciador) = match predicho.rbs {
            Some(rbs) => (rbs.motivo, rbs.espaciador.to_string()),
            None => ("ninguno", "-".to_string()),
        };
        let codon = if predicho.codon_de_inicio.is_empty() {
            "borde"
        } else {
            predicho.codon_de_inicio
        };
        let (izquierda, derecha) = if predicho.hebra == b'-' {
            (predicho.parcial_al_final, predicho.parcial_al_inicio)
        } else {
            (predicho.parcial_al_inicio, predicho.parcial_al_final)
        };
        writeln!(
            destino,
            "{nombre}\tLaCompania\tCDS\t{}\t{}\t{:.1}\t{}\t0\tID=gen_{};partial={}{};start_type={codon};rbs_motif={motivo};rbs_spacer={espaciador}",
            predicho.desde + 1,
            predicho.hasta,
            predicho.puntaje,
            predicho.hebra as char,
            numero + 1,
            u8::from(izquierda),
            u8::from(derecha),
        )?;
    }
    Ok(())
}
