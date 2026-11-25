use nucleo::azar::Azar;
use nucleo::calidad::{
    Descarte, ModoAdaptador, OpcionesRecorte, Recorte, analizar_fastq, buscar_adaptador,
    contar_secuencias, corte_por_calidad, detectar_adaptador, percentil_de_histograma,
    recortar_fastq, recortar_lectura, secuencias_sobrerrepresentadas, tasa_de_duplicacion,
};
use nucleo::formatos::{Lectura, abrir_fastq, escribir_fastq};
use std::fs;
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

const ADAPTADOR_INVENTADO: &[u8] = b"GTCAGTTCAGGACTTAGCCGTATGAC";
const ADAPTADOR_CORTO: &[u8] = b"GATCCTTAGCAAGTCC";
const INSERTO: &[u8] = b"AAAACCCCGGGGTTTT";

static CONTADOR: AtomicUsize = AtomicUsize::new(0);

fn ruta_temporal(nombre: &str) -> PathBuf {
    let numero = CONTADOR.fetch_add(1, Ordering::Relaxed);
    let carpeta =
        std::env::temp_dir().join(format!("calidad-pruebas-{}-{numero}", std::process::id()));
    fs::create_dir_all(&carpeta).unwrap();
    carpeta.join(nombre)
}

fn escribir_texto(nombre: &str, texto: &str) -> PathBuf {
    let ruta = ruta_temporal(nombre);
    fs::write(&ruta, texto).unwrap();
    ruta
}

fn escribir_lecturas(nombre: &str, lecturas: &[Lectura]) -> PathBuf {
    let ruta = ruta_temporal(nombre);
    let mut destino = BufWriter::new(fs::File::create(&ruta).unwrap());
    for lectura in lecturas {
        escribir_fastq(&mut destino, lectura).unwrap();
    }
    destino.flush().unwrap();
    ruta
}

fn leer_lecturas(ruta: &Path) -> Vec<Lectura> {
    abrir_fastq(ruta).unwrap().map(|l| l.unwrap()).collect()
}

fn simbolos(calidades: &[u8]) -> Vec<u8> {
    calidades.iter().map(|q| q + 33).collect()
}

fn unir(partes: &[&[u8]]) -> Vec<u8> {
    partes.concat()
}

fn lectura(secuencia: &[u8], calidades: &[u8]) -> Lectura {
    Lectura {
        nombre: b"prueba".to_vec(),
        secuencia: secuencia.to_vec(),
        calidad: simbolos(calidades),
    }
}

fn lecturas_sinteticas(
    cantidad: usize,
    longitud: usize,
    fraccion_con_adaptador: f64,
    semilla: u64,
) -> Vec<Lectura> {
    let mut azar = Azar::nuevo(semilla);
    let genoma: Vec<u8> = (0..20_000).map(|_| b"ACGT"[azar.entero(4)]).collect();
    let calidad: Vec<u8> = (0..longitud)
        .map(|i| 38 - (13 * i / (longitud - 1)) as u8)
        .collect();
    let mut lecturas = Vec::with_capacity(cantidad);
    for numero in 0..cantidad {
        let con_adaptador = azar.decimal() < fraccion_con_adaptador;
        let inserto = if con_adaptador {
            30 + azar.entero(61)
        } else {
            longitud
        };
        let inicio = azar.entero(genoma.len() - inserto);
        let mut secuencia = genoma[inicio..inicio + inserto].to_vec();
        secuencia.extend_from_slice(ADAPTADOR_INVENTADO);
        while secuencia.len() < longitud {
            secuencia.push(b"ACGT"[azar.entero(4)]);
        }
        secuencia.truncate(longitud);
        lecturas.push(Lectura {
            nombre: format!("r{numero}_{inserto}").into_bytes(),
            secuencia,
            calidad: simbolos(&calidad),
        });
    }
    lecturas
}

fn numero_e_inserto(nombre: &[u8]) -> (usize, usize) {
    let texto = String::from_utf8_lossy(nombre);
    let (numero, inserto) = texto.trim_start_matches('r').split_once('_').unwrap();
    (numero.parse().unwrap(), inserto.parse().unwrap())
}

#[test]
fn ventana_deslizante_desde_el_extremo_3() {
    assert_eq!(corte_por_calidad(&simbolos(&[40; 8]), 4, 20), 8);
    let cola_mala = simbolos(&[40, 40, 40, 40, 40, 40, 10, 10, 10, 10]);
    assert_eq!(corte_por_calidad(&cola_mala, 4, 20), 8);
    let una_mala_al_final = simbolos(&[40, 40, 40, 40, 40, 2]);
    assert_eq!(corte_por_calidad(&una_mala_al_final, 4, 20), 6);
    let toda_mala = simbolos(&[10; 6]);
    assert_eq!(corte_por_calidad(&toda_mala, 4, 20), 0);
    assert_eq!(corte_por_calidad(&toda_mala, 4, 0), 6);
    assert_eq!(corte_por_calidad(&simbolos(&[20, 20, 20, 20]), 4, 20), 4);
    assert_eq!(corte_por_calidad(b"", 4, 20), 0);
    assert_eq!(corte_por_calidad(&simbolos(&[30, 10]), 4, 20), 2);
    let de_a_una = simbolos(&[40, 40, 5, 40, 5, 5]);
    assert_eq!(corte_por_calidad(&de_a_una, 1, 20), 4);
    assert_eq!(corte_por_calidad(&de_a_una, 0, 20), 4);
    assert_eq!(corte_por_calidad(b"  I", 1, 20), 3);
    assert_eq!(corte_por_calidad(b"I  ", 1, 20), 1);
}

#[test]
fn encuentra_el_adaptador_completo_o_parcial() {
    let completo = unir(&[INSERTO, ADAPTADOR_CORTO]);
    assert_eq!(buscar_adaptador(&completo, ADAPTADOR_CORTO, 3), Some(16));
    let dos_errores = unir(&[INSERTO, b"GATACTTAGCTAGTCC"]);
    assert_eq!(buscar_adaptador(&dos_errores, ADAPTADOR_CORTO, 3), Some(16));
    let tres_errores = unir(&[INSERTO, b"GATACTTAGCTAGACC"]);
    assert_eq!(buscar_adaptador(&tres_errores, ADAPTADOR_CORTO, 3), None);
    let tres_al_final = unir(&[INSERTO, b"GAT"]);
    assert_eq!(
        buscar_adaptador(&tres_al_final, ADAPTADOR_CORTO, 3),
        Some(16)
    );
    assert_eq!(buscar_adaptador(&tres_al_final, ADAPTADOR_CORTO, 10), None);
    let dos_al_final = unir(&[INSERTO, b"GA"]);
    assert_eq!(buscar_adaptador(&dos_al_final, ADAPTADOR_CORTO, 3), None);
    let siete_con_error = unir(&[INSERTO, b"GATCCTA"]);
    assert_eq!(buscar_adaptador(&siete_con_error, ADAPTADOR_CORTO, 3), None);
    let ocho_con_error = unir(&[INSERTO, b"GATCCTTT"]);
    assert_eq!(
        buscar_adaptador(&ocho_con_error, ADAPTADOR_CORTO, 3),
        Some(16)
    );
    let dos_veces = unir(&[ADAPTADOR_CORTO, ADAPTADOR_CORTO]);
    assert_eq!(buscar_adaptador(&dos_veces, ADAPTADOR_CORTO, 3), Some(0));
    assert_eq!(buscar_adaptador(&completo, b"", 3), None);
    assert_eq!(buscar_adaptador(b"", ADAPTADOR_CORTO, 3), None);
    assert_eq!(buscar_adaptador(b"GA", ADAPTADOR_CORTO, 0), Some(0));
}

#[test]
fn recorta_una_lectura_paso_por_paso() {
    let opciones = OpcionesRecorte {
        adaptador: ModoAdaptador::Conocido(ADAPTADOR_CORTO.to_vec()),
        ventana: 4,
        calidad_minima: 20,
        maximo_n: 1,
        longitud_minima: 10,
    };
    let con_adaptador = lectura(&unir(&[INSERTO, ADAPTADOR_CORTO]), &[40; 32]);
    assert_eq!(
        recortar_lectura(&con_adaptador, ADAPTADOR_CORTO, &opciones),
        Recorte {
            inicio_adaptador: Some(16),
            fin: 16,
            descarte: None
        }
    );
    let mut calidades = vec![40; 12];
    calidades.extend([10, 10, 10, 10]);
    let cola_mala = lectura(INSERTO, &calidades);
    assert_eq!(
        recortar_lectura(&cola_mala, ADAPTADOR_CORTO, &opciones),
        Recorte {
            inicio_adaptador: None,
            fin: 14,
            descarte: None
        }
    );
    let con_n = lectura(b"AANNCCCCGGGGTTTT", &[40; 16]);
    assert_eq!(
        recortar_lectura(&con_n, ADAPTADOR_CORTO, &opciones).descarte,
        Some(Descarte::DemasiadasN)
    );
    let n_dentro_del_adaptador = lectura(&unir(&[INSERTO, b"GATCCTTAGCAAGTNN"]), &[40; 32]);
    assert_eq!(
        recortar_lectura(&n_dentro_del_adaptador, ADAPTADOR_CORTO, &opciones),
        Recorte {
            inicio_adaptador: Some(16),
            fin: 16,
            descarte: None
        }
    );
    let exigente = OpcionesRecorte {
        longitud_minima: 20,
        ..opciones.clone()
    };
    assert_eq!(
        recortar_lectura(&con_adaptador, ADAPTADOR_CORTO, &exigente).descarte,
        Some(Descarte::MuyCorta)
    );
    let sin_adaptador = recortar_lectura(&con_adaptador, b"", &opciones);
    assert_eq!(sin_adaptador.inicio_adaptador, None);
    assert_eq!(sin_adaptador.fin, 32);
    let mal_formada = Lectura {
        nombre: Vec::new(),
        secuencia: INSERTO.to_vec(),
        calidad: simbolos(&[40; 3]),
    };
    let recorte = recortar_lectura(&mal_formada, ADAPTADOR_CORTO, &opciones);
    assert!(recorte.fin <= 3);
}

#[test]
fn cuartiles_desde_un_histograma() {
    let histograma = [0, 2, 0, 1, 1];
    assert_eq!(percentil_de_histograma(&histograma, 0.0), 1);
    assert_eq!(percentil_de_histograma(&histograma, 0.25), 1);
    assert_eq!(percentil_de_histograma(&histograma, 0.5), 1);
    assert_eq!(percentil_de_histograma(&histograma, 0.75), 3);
    assert_eq!(percentil_de_histograma(&histograma, 1.0), 4);
    assert_eq!(percentil_de_histograma(&[0, 0, 0], 0.5), 0);
    assert_eq!(percentil_de_histograma(&[], 0.5), 0);
}

fn revisar_informe_pequeno(ruta: &Path) {
    let informe = analizar_fastq(ruta).unwrap();
    assert_eq!(informe.lecturas, 4);
    assert_eq!(informe.bases, 14);
    assert_eq!(informe.longitud_minima, 2);
    assert_eq!(informe.longitud_maxima, 4);
    assert_eq!(informe.lecturas_por_longitud, vec![0, 0, 1, 0, 3]);
    assert_eq!(informe.calidad_por_posicion.len(), 4);
    let primera = &informe.calidad_por_posicion[0];
    assert_eq!(primera.lecturas, 4);
    assert!((primera.media - 27.5).abs() < 1e-9);
    assert_eq!(
        (
            primera.cuartil_inferior,
            primera.mediana,
            primera.cuartil_superior
        ),
        (10, 20, 40)
    );
    assert_eq!(informe.calidad_por_posicion[1], *primera);
    let tercera = &informe.calidad_por_posicion[2];
    assert_eq!(tercera.lecturas, 3);
    assert!((tercera.media - 70.0 / 3.0).abs() < 1e-9);
    assert_eq!(
        (
            tercera.cuartil_inferior,
            tercera.mediana,
            tercera.cuartil_superior
        ),
        (10, 20, 40)
    );
    assert_eq!(
        informe.composicion_por_posicion,
        vec![
            [3, 0, 1, 0, 0],
            [0, 2, 1, 1, 0],
            [0, 1, 2, 0, 0],
            [0, 0, 0, 2, 1]
        ]
    );
    assert_eq!(informe.lecturas_por_gc.len(), 101);
    assert_eq!(informe.lecturas_por_gc[0], 1);
    assert_eq!(informe.lecturas_por_gc[50], 2);
    assert_eq!(informe.lecturas_por_gc[100], 1);
    assert_eq!(informe.lecturas_por_gc.iter().sum::<usize>(), 4);
    assert_eq!(informe.bases_q20, 10);
    assert_eq!(informe.bases_q30, 6);
    assert_eq!(informe.bases_n, 1);
    assert_eq!(informe.bases_gc, 7);
    assert!((informe.fraccion_q30() - 6.0 / 14.0).abs() < 1e-9);
    assert!((informe.fraccion_q20() - 10.0 / 14.0).abs() < 1e-9);
    assert!((informe.fraccion_n() - 1.0 / 14.0).abs() < 1e-9);
    assert!((informe.gc_medio() - 7.0 / 13.0).abs() < 1e-9);
    assert_eq!(informe.lecturas_muestreadas, 4);
    assert!((informe.duplicacion - 0.25).abs() < 1e-9);
    assert_eq!(informe.sobrerrepresentadas.len(), 1);
    assert_eq!(informe.sobrerrepresentadas[0].secuencia, b"ACGT");
    assert_eq!(informe.sobrerrepresentadas[0].veces, 2);
    assert!(informe.adaptador.is_none());
}

#[test]
fn informe_de_un_fastq_pequeno_calculado_a_mano() {
    let texto = "@a\nACGT\n+\nIIII\n@b\nACGT\n+\n5555\n@c\nggcn\n+\n++++\n@d\nAT\n+\nII\n";
    revisar_informe_pequeno(&escribir_texto("pequeno.fastq", texto));
    let windows = texto.replace('\n', "\r\n");
    revisar_informe_pequeno(&escribir_texto("pequeno_crlf.fastq", &windows));
}

#[test]
fn informe_de_un_fastq_comprimido() {
    let texto = "@a\nACGT\n+\nIIII\n@b\nACGT\n+\n5555\n@c\nGGCN\n+\n++++\n@d\nAT\n+\nII\n";
    let ruta = ruta_temporal("pequeno.fastq.gz");
    let mut compresor = flate2::write::GzEncoder::new(
        fs::File::create(&ruta).unwrap(),
        flate2::Compression::default(),
    );
    compresor.write_all(texto.as_bytes()).unwrap();
    compresor.finish().unwrap();
    revisar_informe_pequeno(&ruta);
}

#[test]
fn duplicacion_y_secuencias_sobrerrepresentadas() {
    let muestra: Vec<Vec<u8>> = ["AAA", "AAA", "AAA", "CCC", "CCC", "GGG"]
        .iter()
        .map(|s| s.as_bytes().to_vec())
        .collect();
    let veces = contar_secuencias(&muestra);
    assert_eq!(veces.len(), 3);
    assert!((tasa_de_duplicacion(&veces, muestra.len()) - 0.5).abs() < 1e-9);
    let sobre = secuencias_sobrerrepresentadas(&veces, muestra.len());
    assert_eq!(sobre.len(), 2);
    assert_eq!(sobre[0].secuencia, b"AAA");
    assert_eq!(sobre[0].veces, 3);
    assert!((sobre[0].fraccion - 0.5).abs() < 1e-9);
    assert_eq!(sobre[1].secuencia, b"CCC");

    let empate: Vec<Vec<u8>> = ["TT", "AA", "TT", "AA"]
        .iter()
        .map(|s| s.as_bytes().to_vec())
        .collect();
    let veces = contar_secuencias(&empate);
    let sobre = secuencias_sobrerrepresentadas(&veces, empate.len());
    assert_eq!(sobre[0].secuencia, b"AA");
    assert_eq!(sobre[1].secuencia, b"TT");

    let vacia: Vec<Vec<u8>> = Vec::new();
    let veces = contar_secuencias(&vacia);
    assert_eq!(tasa_de_duplicacion(&veces, 0), 0.0);
    assert!(secuencias_sobrerrepresentadas(&veces, 0).is_empty());

    let mut grande: Vec<Vec<u8>> = (0..2993usize)
        .map(|numero| {
            (0..8)
                .map(|cifra| b"ACGT"[(numero >> (2 * cifra)) & 3])
                .collect()
        })
        .collect();
    for _ in 0..4 {
        grande.push(b"TTTTTTTTTT".to_vec());
    }
    for _ in 0..3 {
        grande.push(b"GGGGGGGGGG".to_vec());
    }
    let veces = contar_secuencias(&grande);
    let sobre = secuencias_sobrerrepresentadas(&veces, grande.len());
    assert_eq!(sobre.len(), 1);
    assert_eq!(sobre[0].secuencia, b"TTTTTTTTTT");
    assert!((tasa_de_duplicacion(&veces, grande.len()) - 5.0 / 3000.0).abs() < 1e-9);
}

#[test]
fn detecta_un_adaptador_inventado_desde_los_datos() {
    let lecturas = lecturas_sinteticas(6000, 100, 0.1, 7);
    let muestra: Vec<Vec<u8>> = lecturas.iter().map(|l| l.secuencia.clone()).collect();
    let detectado = detectar_adaptador(&muestra).expect("debe encontrar el adaptador");
    assert_eq!(detectado.secuencia, ADAPTADOR_INVENTADO);
    let esperadas = lecturas
        .iter()
        .filter(|l| numero_e_inserto(&l.nombre).1 < 100)
        .count();
    assert_eq!(detectado.lecturas, esperadas);
    assert_eq!(detectado.revisadas, 6000);
    assert!((detectado.fraccion - esperadas as f64 / 6000.0).abs() < 1e-9);
}

#[test]
fn no_inventa_adaptadores_en_lecturas_limpias() {
    let lecturas = lecturas_sinteticas(6000, 100, 0.0, 11);
    let muestra: Vec<Vec<u8>> = lecturas.iter().map(|l| l.secuencia.clone()).collect();
    assert!(detectar_adaptador(&muestra).is_none());
    let repetidas = vec![b"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG".to_vec(); 500];
    assert!(detectar_adaptador(&repetidas).is_none());
    let cortas = vec![
        b"ACG".to_vec(),
        Vec::new(),
        b"NNNNNNNNNNNNNNNNNNNNNN".to_vec(),
    ];
    assert!(detectar_adaptador(&cortas).is_none());
}

#[test]
fn detecta_el_adaptador_del_ejemplo_de_ecoli() {
    let muestra = nucleo::simulacion::muestra_de_ejemplo().unwrap();
    let secuencias: Vec<Vec<u8>> = muestra
        .lecturas
        .iter()
        .map(|l| l.secuencia.clone())
        .collect();
    let detectado = detectar_adaptador(&secuencias).expect("debe encontrar el adaptador");
    assert_eq!(detectado.secuencia, b"ACGGTAGCTACGATCGGTACTTGCA");
    assert_eq!(detectado.revisadas, secuencias.len());
    assert!(
        detectado.fraccion > 0.06 && detectado.fraccion < 0.10,
        "{}",
        detectado.fraccion
    );
}

#[test]
fn recorta_un_archivo_y_conserva_el_orden() {
    let lecturas = lecturas_sinteticas(6000, 100, 0.1, 7);
    let entrada = escribir_lecturas("sinteticas.fastq", &lecturas);
    let salida = entrada.with_file_name("recortadas.fastq");
    let opciones = OpcionesRecorte {
        longitud_minima: 20,
        ..OpcionesRecorte::default()
    };
    let resumen = recortar_fastq(&entrada, &salida, &opciones).unwrap();
    assert_eq!(resumen.adaptador, ADAPTADOR_INVENTADO);
    assert_eq!(resumen.lecturas_entrada, 6000);
    assert_eq!(resumen.bases_entrada, 600_000);
    assert_eq!(
        resumen.lecturas_salida + resumen.descartadas_por_n + resumen.descartadas_por_longitud,
        6000
    );
    let recortadas = leer_lecturas(&salida);
    assert_eq!(recortadas.len(), resumen.lecturas_salida);
    assert_eq!(
        recortadas.iter().map(|l| l.secuencia.len()).sum::<usize>(),
        resumen.bases_salida
    );
    let con_adaptador = lecturas
        .iter()
        .filter(|l| numero_e_inserto(&l.nombre).1 < 100)
        .count();
    assert!(resumen.lecturas_con_adaptador >= con_adaptador);
    let mut anterior = None;
    for recortada in &recortadas {
        let (numero, inserto) = numero_e_inserto(&recortada.nombre);
        assert!(anterior.is_none_or(|previo| previo < numero));
        anterior = Some(numero);
        let original = &lecturas[numero];
        assert!(original.secuencia.starts_with(&recortada.secuencia));
        assert!(original.calidad.starts_with(&recortada.calidad));
        if inserto < 100 {
            assert_eq!(recortada.secuencia.len(), inserto);
        } else {
            assert!(recortada.secuencia.len() >= 90);
        }
        assert!(
            !recortada
                .secuencia
                .windows(12)
                .any(|ventana| ventana == &ADAPTADOR_INVENTADO[..12])
        );
    }

    let salida_gz = entrada.with_file_name("recortadas.fastq.gz");
    let resumen_gz = recortar_fastq(&entrada, &salida_gz, &opciones).unwrap();
    assert_eq!(resumen_gz, resumen);
    let firma = fs::read(&salida_gz).unwrap();
    assert_eq!(&firma[..2], &[0x1f, 0x8b]);
    assert_eq!(leer_lecturas(&salida_gz), recortadas);

    let informe = analizar_fastq(&salida).unwrap();
    assert_eq!(informe.lecturas, recortadas.len());
    assert!(informe.adaptador.is_none());
}

#[test]
fn recorte_con_adaptador_escrito_o_sin_adaptador() {
    let lecturas = lecturas_sinteticas(2000, 100, 0.1, 3);
    let entrada = escribir_lecturas("sinteticas.fastq", &lecturas);
    let salida = entrada.with_file_name("recortadas.fastq");
    let escrito = OpcionesRecorte {
        adaptador: ModoAdaptador::Conocido(ADAPTADOR_INVENTADO[..20].to_vec()),
        longitud_minima: 20,
        ..OpcionesRecorte::default()
    };
    let resumen = recortar_fastq(&entrada, &salida, &escrito).unwrap();
    assert_eq!(resumen.adaptador, &ADAPTADOR_INVENTADO[..20]);
    assert!(resumen.bases_de_adaptador > 0);
    let ninguno = OpcionesRecorte {
        adaptador: ModoAdaptador::Ninguno,
        longitud_minima: 20,
        ..OpcionesRecorte::default()
    };
    let resumen = recortar_fastq(&entrada, &salida, &ninguno).unwrap();
    assert!(resumen.adaptador.is_empty());
    assert_eq!(resumen.lecturas_con_adaptador, 0);
    assert_eq!(resumen.bases_de_adaptador, 0);
    assert_eq!(resumen.lecturas_salida, 2000);
}

#[test]
fn descarta_por_n_y_por_longitud() {
    let texto = "@larga\nACGTACGTACGT\n+\nIIIIIIIIIIII\n@con_n\nANNNNNNCGTAC\n+\nIIIIIIIIIIII\n@corta\nACG\n+\nIII\n@mala\nACGTACGTACGT\n+\n++++++++++++\n";
    let entrada = escribir_texto("varias.fastq", texto);
    let salida = entrada.with_file_name("varias.recortadas.fastq");
    let opciones = OpcionesRecorte {
        adaptador: ModoAdaptador::Ninguno,
        ventana: 4,
        calidad_minima: 20,
        maximo_n: 5,
        longitud_minima: 10,
    };
    let resumen = recortar_fastq(&entrada, &salida, &opciones).unwrap();
    assert_eq!(resumen.lecturas_entrada, 4);
    assert_eq!(resumen.lecturas_salida, 1);
    assert_eq!(resumen.descartadas_por_n, 1);
    assert_eq!(resumen.descartadas_por_longitud, 2);
    assert_eq!(resumen.bases_por_calidad, 12);
    assert_eq!(resumen.lecturas_recortadas_por_calidad, 1);
    assert_eq!(resumen.bases_q30_entrada, 27);
    assert_eq!(resumen.bases_q30_salida, 12);
    let recortadas = leer_lecturas(&salida);
    assert_eq!(recortadas.len(), 1);
    assert_eq!(recortadas[0].nombre, b"larga");
}

#[test]
fn separa_los_cortes_cortos_de_adaptador() {
    let corta = format!("{}GAT", String::from_utf8_lossy(INSERTO).repeat(4));
    let larga = String::from_utf8_lossy(&unir(&[INSERTO, ADAPTADOR_CORTO])).to_string();
    let texto = format!(
        "@corta\n{corta}\n+\n{}\n@larga\n{larga}\n+\n{}\n",
        "I".repeat(67),
        "I".repeat(32)
    );
    let entrada = escribir_texto("cortes.fastq", &texto);
    let salida = entrada.with_file_name("cortes.recortadas.fastq");
    let opciones = OpcionesRecorte {
        adaptador: ModoAdaptador::Conocido(ADAPTADOR_CORTO.to_vec()),
        longitud_minima: 10,
        ..OpcionesRecorte::default()
    };
    let resumen = recortar_fastq(&entrada, &salida, &opciones).unwrap();
    assert_eq!(resumen.lecturas_con_adaptador, 2);
    assert_eq!(resumen.lecturas_con_adaptador_corto, 1);
    assert_eq!(resumen.bases_de_adaptador, 19);
    assert_eq!(resumen.bases_salida, 64 + 16);
}

#[test]
fn errores_claros_en_vez_de_panico() {
    let no_existe = ruta_temporal("no_existe.fastq");
    assert!(analizar_fastq(&no_existe).is_err());
    let roto = escribir_texto("roto.fastq", "@r1\nACGT\n+\nII\n");
    assert!(analizar_fastq(&roto).is_err());
    let salida = roto.with_file_name("roto.recortadas.fastq");
    assert!(recortar_fastq(&roto, &salida, &OpcionesRecorte::default()).is_err());
    let sin_arroba = escribir_texto("sin_arroba.fastq", ">r1\nACGT\n");
    assert!(analizar_fastq(&sin_arroba).is_err());

    let bueno = escribir_texto("bueno.fastq", "@r1\nACGTACGT\n+\nIIIIIIII\n");
    let salida = bueno.with_file_name("bueno.recortadas.fastq");
    let ventana_cero = OpcionesRecorte {
        ventana: 0,
        ..OpcionesRecorte::default()
    };
    assert!(recortar_fastq(&bueno, &salida, &ventana_cero).is_err());
    assert!(recortar_fastq(&bueno, &bueno, &OpcionesRecorte::default()).is_err());
    assert_eq!(analizar_fastq(&bueno).unwrap().lecturas, 1);
}

#[test]
fn archivo_vacio_y_lecturas_raras() {
    let vacio = escribir_texto("vacio.fastq", "");
    let informe = analizar_fastq(&vacio).unwrap();
    assert_eq!(informe.lecturas, 0);
    assert_eq!(informe.bases, 0);
    assert!(informe.calidad_por_posicion.is_empty());
    assert_eq!(informe.duplicacion, 0.0);
    assert!(informe.adaptador.is_none());
    assert_eq!(informe.fraccion_q30(), 0.0);
    assert_eq!(informe.gc_medio(), 0.0);
    let salida = vacio.with_file_name("vacio.recortadas.fastq");
    let resumen = recortar_fastq(&vacio, &salida, &OpcionesRecorte::default()).unwrap();
    assert_eq!(resumen.lecturas_entrada, 0);
    assert_eq!(fs::read(&salida).unwrap().len(), 0);

    let raras = "@vacia\n\n+\n\n@n\nNNNN\n+\n!!!!\n@iupac\nRYKM\n+\n~~~~\n@espacios\nAC\n+\n  \n";
    let ruta = escribir_texto("raras.fastq", raras);
    let informe = analizar_fastq(&ruta).unwrap();
    assert_eq!(informe.lecturas, 4);
    assert_eq!(informe.bases, 10);
    assert_eq!(informe.bases_n, 8);
    assert_eq!(informe.lecturas_por_gc[50], 1);
    assert_eq!(informe.calidad_por_posicion[0].cuartil_superior, 41);
    let salida = ruta.with_file_name("raras.recortadas.fastq");
    let resumen = recortar_fastq(&ruta, &salida, &OpcionesRecorte::default()).unwrap();
    assert_eq!(resumen.lecturas_salida, 0);
}

#[test]
fn lecturas_largas_limitan_las_posiciones() {
    let secuencia = "ACGT".repeat(400);
    let calidad = "I".repeat(1600);
    let texto = format!("@larga\n{secuencia}\n+\n{calidad}\n");
    let ruta = escribir_texto("larga.fastq", &texto);
    let informe = analizar_fastq(&ruta).unwrap();
    assert_eq!(informe.bases, 1600);
    assert_eq!(informe.bases_q30, 1600);
    assert_eq!(informe.longitud_maxima, 1600);
    assert_eq!(informe.calidad_por_posicion.len(), 1000);
    assert_eq!(informe.composicion_por_posicion.len(), 1000);
}

#[test]
fn umbral_enorme_no_desborda() {
    assert_eq!(corte_por_calidad(&simbolos(&[40; 6]), 4, usize::MAX), 0);
    assert_eq!(corte_por_calidad(&simbolos(&[40; 6]), 6, usize::MAX / 3), 0);
}

#[test]
fn salida_comprimida_sin_lecturas_sigue_siendo_gzip_valido() {
    let ruta = escribir_texto(
        "todas_malas.fastq",
        "@a\nACGT\n+\n!!!!\n@b\nTTTT\n+\n!!!!\n",
    );
    let salida = ruta.with_file_name("todas_malas.recortadas.fastq.gz");
    let resumen = recortar_fastq(&ruta, &salida, &OpcionesRecorte::default()).unwrap();
    assert_eq!(resumen.lecturas_salida, 0);
    let bytes = fs::read(&salida).unwrap();
    assert_eq!(&bytes[..2], &[0x1f, 0x8b]);
    let mut texto = Vec::new();
    flate2::read::MultiGzDecoder::new(&bytes[..])
        .read_to_end(&mut texto)
        .unwrap();
    assert!(texto.is_empty());
    assert_eq!(analizar_fastq(&salida).unwrap().lecturas, 0);
}
