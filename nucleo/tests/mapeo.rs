use nucleo::azar::Azar;
use nucleo::ecoli;
use nucleo::formatos::{Lectura, RegistroSam, abrir, escribir_fastq, leer_sam};
use nucleo::mapeo::*;
use nucleo::secuencia::complemento_inverso;
use nucleo::simulacion;
use std::collections::HashSet;
use std::fs::{self, File};
use std::io::{BufRead, BufWriter, Write};
use std::path::{Path, PathBuf};

fn al_azar(longitud: usize, semilla: u64) -> Vec<u8> {
    let mut azar = Azar::nuevo(semilla);
    (0..longitud).map(|_| b"ACGT"[azar.entero(4)]).collect()
}

fn carpeta(nombre: &str) -> PathBuf {
    let carpeta =
        std::env::temp_dir().join(format!("mapeo-prueba-{}-{nombre}", std::process::id()));
    fs::create_dir_all(&carpeta).unwrap();
    carpeta
}

fn codigo_canonico(kmer: &[u8]) -> Option<u64> {
    let mut directo = 0u64;
    let mut inverso = 0u64;
    for (indice, &base) in kmer.iter().enumerate() {
        let valor = match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        directo = directo * 4 + valor;
        inverso += (3 - valor) << (2 * indice);
    }
    if directo == inverso {
        return None;
    }
    Some(directo.min(inverso))
}

fn minimizadores_a_mano(secuencia: &[u8], k: usize, w: usize) -> HashSet<usize> {
    let mut elegidos = HashSet::new();
    if secuencia.len() < k {
        return elegidos;
    }
    let cantidad = secuencia.len() - k + 1;
    let valores: Vec<Option<u64>> = (0..cantidad)
        .map(|inicio| codigo_canonico(&secuencia[inicio..inicio + k]).map(barajar))
        .collect();
    for fin in w - 1..cantidad {
        let mut mejor: Option<(u64, usize)> = None;
        let primera = fin + 1 - w;
        for (desplazamiento, candidato) in valores[primera..=fin].iter().enumerate() {
            if let Some(valor) = *candidato
                && mejor.is_none_or(|(menor, _)| valor < menor)
            {
                mejor = Some((valor, primera + desplazamiento));
            }
        }
        if let Some((_, posicion)) = mejor {
            elegidos.insert(posicion);
        }
    }
    elegidos
}

#[test]
fn barajar_no_repite_valores() {
    let mut vistos = HashSet::new();
    for codigo in 0..200_000u64 {
        assert!(vistos.insert(barajar(codigo)));
    }
    assert_ne!(barajar(0), 0);
}

#[test]
fn minimizador_canonico_a_mano() {
    let directo = minimizadores(b"AAC", 3, 1);
    assert_eq!(
        directo,
        vec![Minimizador {
            valor: barajar(1),
            posicion: 0,
            reverso: false
        }]
    );
    let opuesto = minimizadores(b"GTT", 3, 1);
    assert_eq!(
        opuesto,
        vec![Minimizador {
            valor: barajar(1),
            posicion: 0,
            reverso: true
        }]
    );
    assert!(minimizadores(b"ACGT", 4, 1).is_empty());
}

#[test]
fn minimizadores_coinciden_con_la_definicion() {
    let mut secuencia = al_azar(3000, 11);
    secuencia[700] = b'N';
    secuencia[1500] = b'N';
    secuencia[1510] = b'R';
    for (k, w) in [(15, 10), (5, 3), (11, 1), (21, 25)] {
        let calculados = minimizadores(&secuencia, k, w);
        let posiciones: HashSet<usize> = calculados.iter().map(|m| m.posicion).collect();
        assert_eq!(posiciones.len(), calculados.len());
        assert_eq!(
            posiciones,
            minimizadores_a_mano(&secuencia, k, w),
            "k={k} w={w}"
        );
        for minimizador in &calculados {
            let kmer = &secuencia[minimizador.posicion..minimizador.posicion + k];
            assert_eq!(Some(minimizador.valor), codigo_canonico(kmer).map(barajar));
        }
        assert!(
            calculados
                .windows(2)
                .all(|par| par[0].posicion < par[1].posicion)
        );
    }
}

#[test]
fn las_dos_hebras_dan_los_mismos_minimizadores() {
    let secuencia = al_azar(5000, 12);
    let opuesta = complemento_inverso(&secuencia);
    let k = 15;
    let directos: HashSet<(u64, usize, bool)> = minimizadores(&secuencia, k, 10)
        .iter()
        .map(|m| (m.valor, m.posicion, m.reverso))
        .collect();
    let reflejados: HashSet<(u64, usize, bool)> = minimizadores(&opuesta, k, 10)
        .iter()
        .map(|m| (m.valor, secuencia.len() - k - m.posicion, !m.reverso))
        .collect();
    assert_eq!(directos, reflejados);
}

#[test]
fn parametros_malos_no_dan_minimizadores() {
    let secuencia = al_azar(100, 13);
    assert!(minimizadores(&secuencia, 0, 10).is_empty());
    assert!(minimizadores(&secuencia, 15, 0).is_empty());
    assert!(minimizadores(&secuencia, 32, 10).is_empty());
    assert!(minimizadores(&secuencia[..10], 15, 1).is_empty());
    assert!(minimizadores(&secuencia[..20], 15, 10).is_empty());
    assert!(minimizadores(b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", 5, 2).is_empty());
    assert!(minimizadores(b"", 15, 10).is_empty());
}

#[test]
fn indice_rechaza_entradas_malas() {
    let referencia = al_azar(1000, 14);
    assert!(construir_indice(b"", 15, 10, 50).is_err());
    assert!(construir_indice(&referencia, 0, 10, 50).is_err());
    assert!(construir_indice(&referencia, 32, 10, 50).is_err());
    assert!(construir_indice(&referencia, 15, 0, 50).is_err());
    assert!(construir_indice(&referencia, 15, 10, 0).is_err());
    assert!(construir_indice(&referencia[..10], 15, 10, 50).is_err());
    assert!(construir_indice(b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", 15, 10, 50).is_err());
    let indice = construir_indice(&referencia, 15, 10, 50).unwrap();
    assert_eq!(indice.longitud_referencia, 1000);
    assert_eq!(indice.descartados, 0);
}

#[test]
fn indice_encuentra_todas_las_copias_y_descarta_las_frecuentes() {
    let repetido = al_azar(40, 15);
    let mut referencia = Vec::new();
    let mut inicios = Vec::new();
    for copia in 0..5 {
        referencia.extend(al_azar(200, 100 + copia));
        inicios.push(referencia.len());
        referencia.extend(&repetido);
    }
    referencia.extend(al_azar(200, 99));
    let propios = minimizadores(&repetido, 15, 10);
    assert!(!propios.is_empty());
    let completo = construir_indice(&referencia, 15, 10, 5).unwrap();
    let elegido = propios[propios.len() / 2];
    let encontrados = completo.buscar(elegido.valor);
    assert_eq!(encontrados.len(), 5);
    let posiciones: Vec<usize> = encontrados.iter().map(|s| s.posicion).collect();
    let esperadas: Vec<usize> = inicios.iter().map(|i| i + elegido.posicion).collect();
    assert_eq!(posiciones, esperadas);
    assert_eq!(completo.descartados, 0);
    let recortado = construir_indice(&referencia, 15, 10, 4).unwrap();
    assert!(recortado.buscar(elegido.valor).is_empty());
    assert!(recortado.descartados >= 1);
    assert!(recortado.semillas.len() < completo.semillas.len());
    assert!(completo.buscar(barajar(123_456_789)).is_empty());
}

#[test]
fn anclas_de_una_lectura_caen_en_su_diagonal() {
    let referencia = al_azar(20_000, 16);
    let indice = construir_indice(&referencia, 15, 10, 50).unwrap();
    let lectura = &referencia[7_000..7_150];
    let anclas = buscar_anclas(lectura, &indice);
    assert!(anclas.len() >= 10);
    for ancla in &anclas {
        assert!(!ancla.reverso);
        assert_eq!(ancla.en_referencia - ancla.en_lectura, 7_000);
    }
    let inversa = complemento_inverso(lectura);
    let anclas = buscar_anclas(&inversa, &indice);
    assert!(anclas.len() >= 10);
    for ancla in &anclas {
        assert!(ancla.reverso);
        assert_eq!(ancla.en_referencia - ancla.en_lectura, 7_000);
    }
    assert!(buscar_anclas(b"ACGT", &indice).is_empty());
}

fn ancla(reverso: bool, en_referencia: usize, en_lectura: usize) -> Ancla {
    Ancla {
        reverso,
        en_referencia,
        en_lectura,
    }
}

#[test]
fn cadena_colineal_a_mano() {
    let anclas = vec![
        ancla(false, 5000, 3),
        ancla(false, 121, 20),
        ancla(false, 100, 0),
        ancla(false, 110, 10),
    ];
    let cadenas = encadenar(&anclas, 15);
    assert_eq!(cadenas.len(), 2);
    assert_eq!(
        cadenas[0].anclas,
        vec![
            ancla(false, 100, 0),
            ancla(false, 110, 10),
            ancla(false, 121, 20)
        ]
    );
    assert!((cadenas[0].puntaje - 34.85).abs() < 1e-9);
    assert_eq!(cadenas[1].anclas, vec![ancla(false, 5000, 3)]);
    assert_eq!(cadenas[1].puntaje, 15.0);
    assert_eq!(cadenas[0].diagonales(), (100, 101));
}

#[test]
fn cadenas_de_hebras_distintas_no_se_mezclan() {
    let anclas = vec![ancla(false, 100, 0), ancla(true, 110, 10)];
    let cadenas = encadenar(&anclas, 15);
    assert_eq!(cadenas.len(), 2);
    assert!(
        cadenas
            .iter()
            .all(|c| c.anclas.len() == 1 && c.puntaje == 15.0)
    );
}

#[test]
fn anclas_fuera_de_orden_o_lejanas_no_se_encadenan() {
    let cruzadas = vec![ancla(false, 100, 20), ancla(false, 120, 5)];
    assert_eq!(encadenar(&cruzadas, 15).len(), 2);
    let desviadas = vec![ancla(false, 100, 0), ancla(false, 400, 10)];
    assert_eq!(encadenar(&desviadas, 15).len(), 2);
    assert!(encadenar(&[], 15).is_empty());
}

#[test]
fn dos_copias_dan_cadenas_empatadas() {
    let anclas = vec![
        ancla(false, 2000, 0),
        ancla(false, 2010, 10),
        ancla(false, 100, 0),
        ancla(false, 110, 10),
    ];
    let cadenas = encadenar(&anclas, 15);
    assert_eq!(cadenas.len(), 2);
    assert_eq!(cadenas[0].puntaje, 25.0);
    assert_eq!(cadenas[1].puntaje, 25.0);
}

const X: &[u8] = b"GATTACAGTGGTCATGAACT";

fn referencia_con_x() -> Vec<u8> {
    let mut referencia = b"CCCCCCCCCC".to_vec();
    referencia.extend_from_slice(X);
    referencia.extend_from_slice(b"CCCCCCCCCC");
    referencia
}

#[test]
fn alineamiento_exacto() {
    let alineado = alinear_con_recorte(X, &referencia_con_x(), 10, 10).unwrap();
    assert_eq!(alineado.posicion, 10);
    assert_eq!(alineado.cigar, vec![(20, b'M')]);
    assert_eq!(alineado.puntaje, 20);
    assert_eq!(alineado.coincidencias, 20);
    assert_eq!(alineado.diferencias, 0);
}

#[test]
fn alineamiento_con_un_cambio() {
    let mut lectura = X.to_vec();
    lectura[7] = b'C';
    let alineado = alinear_con_recorte(&lectura, &referencia_con_x(), 10, 10).unwrap();
    assert_eq!(alineado.posicion, 10);
    assert_eq!(alineado.cigar, vec![(20, b'M')]);
    assert_eq!(alineado.puntaje, 15);
    assert_eq!(alineado.diferencias, 1);
}

#[test]
fn alineamiento_con_deleccion() {
    let mut lectura = X[..11].to_vec();
    lectura.extend_from_slice(&X[12..]);
    let alineado = alinear_con_recorte(&lectura, &referencia_con_x(), 10, 10).unwrap();
    assert_eq!(alineado.posicion, 10);
    assert_eq!(alineado.cigar, vec![(11, b'M'), (1, b'D'), (8, b'M')]);
    assert_eq!(alineado.puntaje, 12);
    assert_eq!(alineado.borradas, 1);
}

#[test]
fn alineamiento_con_insercion() {
    let mut lectura = X[..11].to_vec();
    lectura.push(b'C');
    lectura.extend_from_slice(&X[11..]);
    let alineado = alinear_con_recorte(&lectura, &referencia_con_x(), 10, 10).unwrap();
    assert_eq!(alineado.posicion, 10);
    assert_eq!(alineado.cigar, vec![(11, b'M'), (1, b'I'), (9, b'M')]);
    assert_eq!(alineado.puntaje, 13);
    assert_eq!(alineado.insertadas, 1);
}

#[test]
fn adaptador_al_final_se_recorta() {
    let mut lectura = X.to_vec();
    lectura.extend_from_slice(b"GGGGGGGGGG");
    let alineado = alinear_con_recorte(&lectura, &referencia_con_x(), 10, 10).unwrap();
    assert_eq!(alineado.posicion, 10);
    assert_eq!(alineado.cigar, vec![(20, b'M'), (10, b'S')]);
    assert_eq!(alineado.puntaje, 20);
}

#[test]
fn adaptador_al_inicio_se_recorta() {
    let mut lectura = b"TTTTTTTTTT".to_vec();
    lectura.extend_from_slice(X);
    let alineado = alinear_con_recorte(&lectura, &referencia_con_x(), 0, 0).unwrap();
    assert_eq!(alineado.posicion, 10);
    assert_eq!(alineado.cigar, vec![(10, b'S'), (20, b'M')]);
    assert_eq!(alineado.puntaje, 20);
}

#[test]
fn un_cambio_en_la_ultima_base_no_se_recorta() {
    let mut lectura = X.to_vec();
    lectura[19] = b'A';
    let alineado = alinear_con_recorte(&lectura, &referencia_con_x(), 10, 10).unwrap();
    assert_eq!(alineado.cigar, vec![(20, b'M')]);
    assert_eq!(alineado.puntaje, 15);
}

#[test]
fn dos_cambios_al_final_se_recortan() {
    let mut lectura = X.to_vec();
    lectura[18] = b'A';
    lectura[19] = b'A';
    let alineado = alinear_con_recorte(&lectura, &referencia_con_x(), 10, 10).unwrap();
    assert_eq!(alineado.cigar, vec![(18, b'M'), (2, b'S')]);
    assert_eq!(alineado.puntaje, 18);
}

#[test]
fn sin_huecos_da_lo_mismo_que_la_tabla_completa() {
    let referencia = al_azar(20_000, 29);
    let mut azar = Azar::nuevo(30);
    for _ in 0..300 {
        let inicio = 50 + azar.entero(19_000);
        let mut lectura = referencia[inicio..inicio + 120].to_vec();
        if azar.entero(2) == 1 {
            let lugar = azar.entero(120);
            lectura[lugar] = if lectura[lugar] == b'G' { b'T' } else { b'G' };
        }
        let diagonal = inicio as i64;
        let directo = alinear_con_recorte(&lectura, &referencia, diagonal, diagonal).unwrap();
        let con_tabla =
            alinear_con_recorte(&lectura, &referencia, diagonal - 1, diagonal + 1).unwrap();
        assert_eq!(directo, con_tabla);
        assert_eq!(directo.posicion, inicio);
    }
}

#[test]
fn alineamiento_imposible_devuelve_nada() {
    let referencia = referencia_con_x();
    assert!(alinear_con_recorte(b"", &referencia, 10, 10).is_none());
    assert!(alinear_con_recorte(X, b"", 0, 0).is_none());
    assert!(alinear_con_recorte(X, &referencia, 1000, 1000).is_none());
    assert!(alinear_con_recorte(X, &referencia, 10, 5).is_none());
    assert!(alinear_con_recorte(b"GGGGGGGG", b"CCCCCCCCCCCC", 0, 0).is_none());
}

#[test]
fn calidad_de_mapeo_a_mano() {
    assert_eq!(calidad_de_mapeo(100, 100, 20), 0);
    assert_eq!(calidad_de_mapeo(100, 0, 20), 60);
    assert_eq!(calidad_de_mapeo(100, 90, 20), 18);
    assert_eq!(calidad_de_mapeo(100, 90, 5), 9);
    assert_eq!(calidad_de_mapeo(0, 0, 10), 0);
    assert_eq!(calidad_de_mapeo(-5, 0, 10), 0);
    assert_eq!(calidad_de_mapeo(50, 80, 10), 0);
    assert_eq!(calidad_de_mapeo(1, 0, 10), 0);
}

fn indice_de(referencia: &[u8]) -> Indice {
    construir_indice(referencia, K, W, FRECUENCIA_MAXIMA).unwrap()
}

#[test]
fn mapea_una_lectura_en_cada_hebra() {
    let referencia = al_azar(50_000, 17);
    let indice = indice_de(&referencia);
    let lectura = referencia[12_345..12_445].to_vec();
    let mapeo = mapear_lectura(&lectura, &indice, &referencia, 0).unwrap();
    assert!(!mapeo.reverso);
    assert_eq!(mapeo.alineamiento.posicion, 12_345);
    assert_eq!(mapeo.alineamiento.cigar, vec![(100, b'M')]);
    assert_eq!(mapeo.calidad_mapeo, 60);
    let inversa = complemento_inverso(&lectura);
    let mapeo = mapear_lectura(&inversa, &indice, &referencia, 0).unwrap();
    assert!(mapeo.reverso);
    assert_eq!(mapeo.alineamiento.posicion, 12_345);
    assert_eq!(mapeo.alineamiento.cigar, vec![(100, b'M')]);
    assert_eq!(mapeo.calidad_mapeo, 60);
}

#[test]
fn mapea_con_cambios_huecos_y_adaptador() {
    let referencia = al_azar(50_000, 18);
    let indice = indice_de(&referencia);
    let mut lectura = referencia[30_000..30_060].to_vec();
    lectura.extend_from_slice(&referencia[30_061..30_120]);
    lectura[20] = if lectura[20] == b'A' { b'C' } else { b'A' };
    lectura.extend_from_slice(b"ACGGTAGCTACGATCGGTACTTGCA");
    let mapeo = mapear_lectura(&lectura, &indice, &referencia, 0).unwrap();
    assert_eq!(mapeo.alineamiento.posicion, 30_000);
    let ultima = *mapeo.alineamiento.cigar.last().unwrap();
    assert_eq!(ultima.1, b'S');
    assert!(ultima.0 >= 20);
    assert_eq!(mapeo.alineamiento.diferencias, 1);
    assert_eq!(mapeo.alineamiento.borradas, 1);
    assert_eq!(mapeo.alineamiento.insertadas, 0);
}

#[test]
fn lectura_repetida_tiene_calidad_cero() {
    let copia = al_azar(500, 19);
    let mut referencia = al_azar(20_000, 20);
    referencia.extend_from_slice(&copia);
    referencia.extend(al_azar(20_000, 21));
    referencia.extend_from_slice(&copia);
    referencia.extend(al_azar(10_000, 22));
    let indice = indice_de(&referencia);
    let lectura = &copia[100..250];
    let mut posiciones = HashSet::new();
    for desempate in 0..20 {
        let mapeo = mapear_lectura(lectura, &indice, &referencia, desempate).unwrap();
        assert_eq!(mapeo.calidad_mapeo, 0);
        posiciones.insert(mapeo.alineamiento.posicion);
    }
    assert_eq!(posiciones, HashSet::from([20_100, 40_600]));
}

#[test]
fn repeticion_en_tandem_corta_tiene_calidad_cero() {
    let unidad = al_azar(40, 31);
    let mut referencia = al_azar(3_000, 32);
    for _ in 0..5 {
        referencia.extend_from_slice(&unidad);
    }
    referencia.extend(al_azar(3_000, 33));
    let indice = indice_de(&referencia);
    let lectura = referencia[3_050..3_150].to_vec();
    let mut posiciones = HashSet::new();
    for desempate in 0..40 {
        let mapeo = mapear_lectura(&lectura, &indice, &referencia, desempate).unwrap();
        assert_eq!(mapeo.calidad_mapeo, 0);
        assert_eq!(mapeo.alineamiento.cigar, vec![(100, b'M')]);
        posiciones.insert(mapeo.alineamiento.posicion);
    }
    assert!(posiciones.len() >= 2);
    assert!(
        posiciones
            .iter()
            .all(|posicion| (posicion - 3_010) % 40 == 0)
    );
}

#[test]
fn lectura_ajena_o_corta_no_se_mapea() {
    let referencia = al_azar(50_000, 23);
    let indice = indice_de(&referencia);
    assert!(mapear_lectura(&al_azar(150, 24), &indice, &referencia, 0).is_none());
    assert!(mapear_lectura(&referencia[100..115], &indice, &referencia, 0).is_none());
    assert!(mapear_lectura(b"", &indice, &referencia, 0).is_none());
    assert!(mapear_lectura(&[b'N'; 150], &indice, &referencia, 0).is_none());
    let larga = &referencia[1000..1000 + LONGITUD_MAXIMA_DE_LECTURA + 1];
    assert!(mapear_lectura(larga, &indice, &referencia, 0).is_none());
}

fn escribir_lecturas(ruta: &Path, lecturas: &[Lectura]) {
    let mut destino = BufWriter::new(File::create(ruta).unwrap());
    for lectura in lecturas {
        escribir_fastq(&mut destino, lectura).unwrap();
    }
    destino.flush().unwrap();
}

fn leer_lineas(ruta: &Path) -> Vec<String> {
    abrir(ruta)
        .unwrap()
        .lines()
        .map(|linea| linea.unwrap())
        .collect()
}

#[test]
fn mapear_fastq_escribe_un_sam_correcto() {
    let referencia = al_azar(30_000, 25);
    let indice = indice_de(&referencia);
    let calidad: Vec<u8> = (0..100).map(|i| b'!' + (i % 40) as u8).collect();
    let directa = referencia[1_000..1_100].to_vec();
    let inversa = complemento_inverso(&referencia[5_000..5_100]);
    let lecturas = vec![
        Lectura {
            nombre: b"uno".to_vec(),
            secuencia: directa.clone(),
            calidad: calidad.clone(),
        },
        Lectura {
            nombre: b"dos con comentario".to_vec(),
            secuencia: inversa.clone(),
            calidad: calidad.clone(),
        },
        Lectura {
            nombre: b"tres".to_vec(),
            secuencia: al_azar(100, 26),
            calidad: calidad.clone(),
        },
        Lectura {
            nombre: b"cuatro".to_vec(),
            secuencia: referencia[1_050..1_150].to_vec(),
            calidad: calidad.clone(),
        },
        Lectura {
            nombre: b"cinco".to_vec(),
            secuencia: referencia[29_900..30_000].to_vec(),
            calidad: calidad.clone(),
        },
    ];
    let carpeta = carpeta("sam");
    let entrada = carpeta.join("lecturas.fastq");
    let salida = carpeta.join("lecturas.sam");
    escribir_lecturas(&entrada, &lecturas);
    let resumen = mapear_fastq(&entrada, &indice, &referencia, "prueba", &salida).unwrap();
    let lineas = leer_lineas(&salida);
    assert_eq!(lineas[0], "@HD\tVN:1.6\tSO:unsorted");
    assert_eq!(lineas[1], "@SQ\tSN:prueba\tLN:30000");
    let registros: Vec<_> = lineas
        .iter()
        .filter_map(|linea| leer_sam(linea.as_bytes()).unwrap())
        .collect();
    assert_eq!(registros.len(), 5);
    assert_eq!(registros[0].nombre, b"uno");
    assert_eq!(registros[0].bandera, 0);
    assert_eq!(registros[0].referencia, "prueba");
    assert_eq!(registros[0].posicion, 1_001);
    assert_eq!(registros[0].calidad_mapeo, 60);
    assert_eq!(registros[0].cigar, vec![(100, b'M')]);
    assert_eq!(registros[0].secuencia, directa);
    assert_eq!(registros[0].calidad, calidad);
    assert_eq!(registros[1].nombre, b"dos");
    assert_eq!(registros[1].bandera, 16);
    assert_eq!(registros[1].posicion, 5_001);
    assert_eq!(registros[1].secuencia, referencia[5_000..5_100].to_vec());
    let al_reves: Vec<u8> = calidad.iter().rev().copied().collect();
    assert_eq!(registros[1].calidad, al_reves);
    assert_eq!(registros[2].bandera, 4);
    assert_eq!(registros[2].posicion, 0);
    assert_eq!(registros[2].referencia, "*");
    assert!(registros[2].cigar.is_empty());
    assert_eq!(registros[2].calidad_mapeo, 0);
    assert_eq!(registros[3].posicion, 1_051);
    assert_eq!(registros[4].posicion, 29_901);
    assert_eq!(registros[4].cigar, vec![(100, b'M')]);
    assert_eq!(resumen.lecturas, 5);
    assert_eq!(resumen.mapeadas, 4);
    assert_eq!(resumen.en_hebra_reversa, 1);
    assert_eq!(resumen.calidades[60], 4);
    assert_eq!(resumen.calidades.iter().sum::<usize>(), 4);
    assert_eq!(resumen.coincidencias, 400);
    assert_eq!(resumen.diferencias, 0);
    assert_eq!(resumen.bases_cubiertas, 350);
    assert_eq!(
        resumen.tramos,
        vec![
            Tramo {
                inicio: 1_000,
                fin: 1_150,
                profundidad_media: 200.0 / 150.0
            },
            Tramo {
                inicio: 5_000,
                fin: 5_100,
                profundidad_media: 1.0
            },
            Tramo {
                inicio: 29_900,
                fin: 30_000,
                profundidad_media: 1.0
            }
        ]
    );
    assert!((resumen.profundidad_media - 400.0 / 350.0).abs() < 1e-12);
    assert_eq!(resumen.largo_ventana, 300);
    assert_eq!(resumen.profundidad_por_ventana.len(), 100);
    assert!((resumen.profundidad_por_ventana[3] - 200.0 / 300.0).abs() < 1e-12);
    assert_eq!(resumen.profundidad_por_ventana[4], 0.0);
    assert!((resumen.profundidad_por_ventana[16] - 100.0 / 300.0).abs() < 1e-12);
    assert!((resumen.profundidad_por_ventana[99] - 100.0 / 300.0).abs() < 1e-12);
    let suma: f64 = resumen.profundidad_por_ventana.iter().sum();
    assert!((suma * 300.0 - 400.0).abs() < 1e-9);
    assert_eq!(resumen.minimizadores_en_indice, indice.semillas.len());
    assert_eq!(resumen.minimizadores_descartados, indice.descartados);
    fs::remove_dir_all(&carpeta).unwrap();
}

#[test]
fn largo_de_ventana_segun_la_referencia() {
    assert_eq!(largo_de_ventana(4_641_652), LARGO_VENTANA);
    assert_eq!(largo_de_ventana(1_000_000), LARGO_VENTANA);
    assert_eq!(largo_de_ventana(999_999), 9_999);
    assert_eq!(largo_de_ventana(30_000), 300);
    assert_eq!(largo_de_ventana(50), 1);
    assert_eq!(largo_de_ventana(0), 1);
}

#[test]
fn mapear_lecturas_construye_el_indice_y_da_lo_mismo() {
    let referencia = al_azar(40_000, 31);
    let mut azar = Azar::nuevo(32);
    let mut lecturas = Vec::new();
    for numero in 0..200 {
        let inicio = azar.entero(39_800);
        let tramo = &referencia[inicio..inicio + 150];
        let secuencia = if numero % 2 == 0 {
            tramo.to_vec()
        } else {
            complemento_inverso(tramo)
        };
        lecturas.push(Lectura {
            nombre: format!("l{numero}").into_bytes(),
            secuencia,
            calidad: vec![b'I'; 150],
        });
    }
    let carpeta = carpeta("envoltura");
    let entrada = carpeta.join("lecturas.fastq");
    let con_indice = carpeta.join("con_indice.sam");
    let directo = carpeta.join("directo.sam");
    escribir_lecturas(&entrada, &lecturas);
    let indice = indice_de(&referencia);
    let esperado = mapear_fastq(&entrada, &indice, &referencia, "r", &con_indice).unwrap();
    let obtenido = mapear_lecturas(&entrada, &referencia, "r", &directo).unwrap();
    assert_eq!(obtenido, esperado);
    assert_eq!(obtenido.mapeadas, 200);
    assert_eq!(leer_lineas(&con_indice), leer_lineas(&directo));
    assert!(mapear_lecturas(&entrada, b"ACGT", "r", &directo).is_err());
    fs::remove_dir_all(&carpeta).unwrap();
}

#[test]
fn lectura_con_n_se_mapea_con_una_diferencia() {
    let referencia = al_azar(50_000, 33);
    let indice = indice_de(&referencia);
    let mut lectura = referencia[8_000..8_150].to_vec();
    lectura[75] = b'N';
    let mapeo = mapear_lectura(&lectura, &indice, &referencia, 0).unwrap();
    assert_eq!(mapeo.alineamiento.posicion, 8_000);
    assert_eq!(mapeo.alineamiento.cigar, vec![(150, b'M')]);
    assert_eq!(mapeo.alineamiento.diferencias, 1);
    assert_eq!(mapeo.calidad_mapeo, 60);
}

#[test]
fn mapear_fastq_rechaza_entradas_malas() {
    let referencia = al_azar(30_000, 27);
    let indice = indice_de(&referencia);
    let carpeta = carpeta("malas");
    let salida = carpeta.join("salida.sam");
    let inexistente = carpeta.join("no_existe.fastq");
    assert!(mapear_fastq(&inexistente, &indice, &referencia, "r", &salida).is_err());
    let roto = carpeta.join("roto.fastq");
    fs::write(&roto, b"@a\nACGT\n+\nII\n").unwrap();
    assert!(mapear_fastq(&roto, &indice, &referencia, "r", &salida).is_err());
    let bueno = carpeta.join("bueno.fastq");
    fs::write(&bueno, b"@a\nACGT\n+\nIIII\n").unwrap();
    let otra = al_azar(10_000, 28);
    assert!(mapear_fastq(&bueno, &indice, &otra, "r", &salida).is_err());
    assert!(mapear_fastq(&bueno, &indice, &referencia, "r", &bueno).is_err());
    assert_eq!(fs::read(&bueno).unwrap(), b"@a\nACGT\n+\nIIII\n");
    let mismo_por_otro_camino = carpeta.join(".").join("bueno.fastq");
    assert!(mapear_fastq(&bueno, &indice, &referencia, "r", &mismo_por_otro_camino).is_err());
    assert_eq!(fs::read(&bueno).unwrap(), b"@a\nACGT\n+\nIIII\n");
    let vacio = carpeta.join("vacio.fastq");
    fs::write(&vacio, b"").unwrap();
    let resumen = mapear_fastq(&vacio, &indice, &referencia, "r", &salida).unwrap();
    assert_eq!(resumen.lecturas, 0);
    assert_eq!(resumen.mapeadas, 0);
    assert_eq!(resumen.profundidad_media, 0.0);
    assert!(resumen.tramos.is_empty());
    fs::remove_dir_all(&carpeta).unwrap();
}

#[test]
fn exactitud_en_la_muestra_de_ejemplo() {
    let muestra = simulacion::muestra_de_ejemplo().unwrap();
    let genoma = ecoli::genoma();
    let indice = indice_de(genoma);
    let carpeta = carpeta("ejemplo");
    let entrada = carpeta.join("ejemplo.fastq");
    let salida = carpeta.join("ejemplo.sam");
    escribir_lecturas(&entrada, &muestra.lecturas);
    let resumen = mapear_fastq(&entrada, &indice, genoma, ecoli::IDENTIFICADOR, &salida).unwrap();
    let registros: Vec<_> = leer_lineas(&salida)
        .iter()
        .filter_map(|linea| leer_sam(linea.as_bytes()).unwrap())
        .collect();
    assert_eq!(registros.len(), muestra.verdad.len());
    let mut correctas = 0;
    let mut en_copia_igual_de_buena = 0;
    let mut erradas = 0;
    let mut con_calidad = 0;
    let mut correctas_con_calidad = 0;
    for (obtenido, esperado) in registros.iter().zip(&muestra.verdad) {
        assert_eq!(obtenido.nombre, esperado.nombre);
        let acierta = obtenido.mapeado()
            && obtenido.hebra_reversa() == esperado.hebra_reversa()
            && obtenido.posicion.abs_diff(esperado.posicion) <= 5;
        if obtenido.mapeado() && obtenido.calidad_mapeo > 0 {
            con_calidad += 1;
            if acierta {
                correctas_con_calidad += 1;
            }
        }
        if acierta {
            correctas += 1;
        } else if obtenido.mapeado()
            && puntaje_en_su_lugar(esperado, genoma) <= puntaje_en_su_lugar(obtenido, genoma)
        {
            en_copia_igual_de_buena += 1;
        } else {
            erradas += 1;
        }
    }
    let total = registros.len() as f64;
    eprintln!(
        "lecturas {} mapeadas {} en su lugar exacto {} ({:.2} %) en una copia igual de buena {} ({:.2} %) erradas {} con MAPQ > 0: {} de {} en su lugar ({:.3} %) MAPQ 0: {}",
        registros.len(),
        resumen.mapeadas,
        correctas,
        100.0 * correctas as f64 / total,
        en_copia_igual_de_buena,
        100.0 * en_copia_igual_de_buena as f64 / total,
        erradas,
        correctas_con_calidad,
        con_calidad,
        100.0 * correctas_con_calidad as f64 / con_calidad as f64,
        resumen.calidades[0]
    );
    assert!(resumen.mapeadas as f64 >= 0.99 * total);
    assert!(correctas as f64 >= 0.96 * total);
    assert!((correctas + en_copia_igual_de_buena) as f64 >= 0.99 * total);
    assert!(erradas as f64 <= 0.001 * total);
    assert!(correctas_con_calidad as f64 >= 0.999 * con_calidad as f64);
    fs::remove_dir_all(&carpeta).unwrap();
}

fn puntaje_en_su_lugar(registro: &RegistroSam, genoma: &[u8]) -> i32 {
    let recorte_inicial = match registro.cigar.first() {
        Some(&(largo, b'S')) => largo as i64,
        _ => 0,
    };
    let diagonal = registro.posicion as i64 - 1 - recorte_inicial;
    alinear_con_recorte(&registro.secuencia, genoma, diagonal - 5, diagonal + 5)
        .map(|alineado| alineado.puntaje)
        .unwrap_or(0)
}
