use flate2::Compression;
use flate2::write::GzEncoder;
use nucleo::azar::Azar;
use nucleo::ecoli::{self, Gen};
use nucleo::formatos::{RegistroSam, escribir_cabecera_sam, escribir_sam};
use nucleo::simulacion::{
    Mutacion, OpcionesSimulacion, muestra_de_ejemplo, mutaciones_de_resistencia, simular_lecturas,
};
use nucleo::variantes::*;
use std::io::Write;
use std::path::PathBuf;

const CABECERA_CORTA: &str = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr\tLN:10\n";

fn linea_sam(
    nombre: &str,
    bandera: u16,
    posicion: usize,
    mapq: u8,
    cigar: &str,
    bases: &str,
) -> String {
    let calidad = "I".repeat(bases.len());
    format!("{nombre}\t{bandera}\tchr\t{posicion}\t{mapq}\t{cigar}\t*\t0\t0\t{bases}\t{calidad}\n")
}

fn sam_de(longitud: usize, lecturas: &[(usize, &str, &str)]) -> Vec<u8> {
    let mut texto = format!("@SQ\tSN:chr\tLN:{longitud}\n");
    for (numero, (posicion, cigar, bases)) in lecturas.iter().enumerate() {
        texto.push_str(&linea_sam(
            &format!("l{numero}"),
            0,
            *posicion,
            60,
            cigar,
            bases,
        ));
    }
    texto.into_bytes()
}

fn cambio(posicion: usize, referencia: &str, alternativa: &str) -> Cambio {
    Cambio {
        posicion,
        referencia: referencia.as_bytes().to_vec(),
        alternativa: alternativa.as_bytes().to_vec(),
    }
}

fn archivo_temporal(nombre: &str) -> PathBuf {
    std::env::temp_dir().join(format!("prueba_variantes_{}_{nombre}", std::process::id()))
}

fn escribir_verdad_sam(ruta: &PathBuf, nombre: &str, longitud: usize, registros: &[RegistroSam]) {
    let mut destino = std::io::BufWriter::new(std::fs::File::create(ruta).unwrap());
    escribir_cabecera_sam(&mut destino, nombre, longitud).unwrap();
    for registro in registros {
        escribir_sam(&mut destino, registro).unwrap();
    }
    destino.flush().unwrap();
}

#[test]
fn apila_bases_borrados_inserciones_y_filtra_lecturas() {
    let referencia = b"ACGTACGTAC";
    let mut texto = CABECERA_CORTA.to_string();
    texto.push_str(&linea_sam("r1", 0, 1, 60, "5M", "ACGTA"));
    texto.push_str("r2\t4\t*\t0\t0\t*\t*\t0\t0\tACGTA\tIIIII\n");
    texto.push_str(&linea_sam("r3", 256, 1, 60, "5M", "ACGTA"));
    texto.push_str(&linea_sam("r4", 0, 1, 10, "5M", "ACGTA"));
    texto.push_str("r5\t16\tchr\t2\t60\t4M\t*\t0\t0\tCTTA\tII#I\n");
    texto.push_str(&linea_sam("r6", 0, 3, 60, "2M1D2M", "GTCG"));
    texto.push_str(&linea_sam("r7", 0, 4, 60, "1M2I2M", "TGGAC"));
    texto.push_str(&linea_sam("r8", 0, 6, 60, "2S3M", "NNCGT"));
    texto.push_str(&linea_sam("r9", 1024, 1, 60, "5M", "ACGTA"));
    texto.push_str(&linea_sam("r10", 2048, 1, 60, "5M", "ACGTA"));
    texto.push_str("r11\t0\totro\t1\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n");
    let apilamiento = apilar(texto.as_bytes(), referencia, &Umbrales::default()).unwrap();
    assert_eq!(apilamiento.nombre_referencia, "chr");
    assert_eq!(apilamiento.lecturas_usadas, 5);
    assert_eq!(apilamiento.lecturas_descartadas, 6);
    let esperadas: [[u32; 4]; 10] = [
        [1, 0, 0, 0],
        [0, 2, 0, 0],
        [0, 0, 2, 1],
        [0, 0, 0, 3],
        [3, 0, 0, 0],
        [0, 3, 0, 0],
        [0, 0, 2, 0],
        [0, 0, 0, 1],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
    ];
    assert_eq!(apilamiento.bases, esperadas.to_vec());
    assert_eq!(apilamiento.borrados, vec![0, 0, 0, 0, 1, 0, 0, 0, 0, 0]);
    assert_eq!(apilamiento.cruzan, vec![1, 2, 3, 4, 2, 2, 1, 0, 0, 0]);
    assert_eq!(apilamiento.indels.len(), 2);
    assert_eq!(apilamiento.indels[&cambio(3, "TA", "T")], 1);
    assert_eq!(apilamiento.indels[&cambio(3, "T", "TGG")], 1);
    assert!((apilamiento.suma_error[1] as f64 - 2e-4).abs() < 1e-7);
    assert!((apilamiento.suma_error[3] as f64 - 3e-4).abs() < 1e-7);
}

#[test]
fn sin_cabecera_usa_la_primera_referencia_mapeada() {
    let texto = "a\t4\t*\t0\t0\t*\t*\t0\t0\tAC\tII\nb\t0\tuno\t1\t60\t2M\t*\t0\t0\tAC\tII\nc\t0\tdos\t1\t60\t2M\t*\t0\t0\tAC\tII\n";
    let apilamiento = apilar(texto.as_bytes(), b"ACGT", &Umbrales::default()).unwrap();
    assert_eq!(apilamiento.nombre_referencia, "uno");
    assert_eq!(apilamiento.lecturas_usadas, 1);
    assert_eq!(apilamiento.lecturas_descartadas, 2);
}

#[test]
fn elige_la_secuencia_de_cabecera_con_la_misma_longitud() {
    let texto = "@SQ\tSN:plasmido\tLN:3\n@SQ\tSN:cromosoma\tLN:4\nb\t0\tcromosoma\t1\t60\t2M\t*\t0\t0\tAC\tII\nc\t0\tplasmido\t1\t60\t2M\t*\t0\t0\tAC\tII\n";
    let apilamiento = apilar(texto.as_bytes(), b"ACGT", &Umbrales::default()).unwrap();
    assert_eq!(apilamiento.nombre_referencia, "cromosoma");
    assert_eq!(apilamiento.lecturas_usadas, 1);
}

#[test]
fn rechaza_un_sam_de_otra_referencia_o_mal_formado() {
    let umbrales = Umbrales::default();
    let otra = "@SQ\tSN:chr\tLN:99\nb\t0\tchr\t1\t60\t2M\t*\t0\t0\tAC\tII\n";
    let error = apilar(otra.as_bytes(), b"ACGT", &umbrales).unwrap_err();
    assert!(error.contains("99"), "{error}");
    let cigar_corto = "b\t0\tchr\t1\t60\t5M\t*\t0\t0\tACGT\tIIII\n";
    assert!(apilar(cigar_corto.as_bytes(), b"ACGTACGT", &umbrales).is_err());
    let calidad_corta = "b\t0\tchr\t1\t60\t4M\t*\t0\t0\tACGT\tII\n";
    assert!(apilar(calidad_corta.as_bytes(), b"ACGTACGT", &umbrales).is_err());
    let pocas_columnas = "b\t0\tchr\t1\t60\n";
    let error = apilar(pocas_columnas.as_bytes(), b"ACGT", &umbrales).unwrap_err();
    assert!(error.contains("Línea 1"), "{error}");
    assert!(apilar(&b"@SQ\tSN:chr\tLN:4\n"[..], b"ACGT", &umbrales).is_ok());
}

#[test]
fn lecturas_raras_no_rompen_nada() {
    let umbrales = Umbrales::default();
    let mut texto = String::new();
    texto.push_str("a\t0\tchr\t3\t60\t6M\t*\t0\t0\tGTACGT\t*\n");
    texto.push_str("b\t0\tchr\t4\t60\t1M3D2M\t*\t0\t0\tTAC\tIII\n");
    texto.push_str("c\t0\tchr\t1\t60\t2I2M\t*\t0\t0\tGGAC\tIIII\n");
    texto.push_str("d\t0\tchr\t1\t60\t1D2M\t*\t0\t0\tCT\tII\n");
    texto.push_str("e\t0\tchr\t2\t60\t1M5N1M2H\t*\t0\t0\tCC\tII\n");
    texto.push_str("f\t0\tchr\t1\t60\t2M\t*\t0\t0\t*\t*\n");
    texto.push_str("g\t0\tchr\t1\t60\t2M\t*\t0\t0\tAN\tII\n");
    texto.push_str("h\t0\tchr\t18446744073709551615\t60\t2M\t*\t0\t0\tAC\tII\n");
    texto.push_str("i\t0\tchr\t1\t60\t1M4294967295D1M\t*\t0\t0\tAC\tII\n");
    texto.push_str("j\t0\tchr\t5\t60\t1M0I0D1M\t*\t0\t0\tAA\tII\n");
    let referencia = b"ACGTA";
    let apilamiento = apilar(texto.as_bytes(), referencia, &umbrales).unwrap();
    assert_eq!(apilamiento.bases.len(), 5);
    assert_eq!(apilamiento.bases[2][2], 1);
    assert_eq!(apilamiento.bases[4][0], 2);
    assert_eq!(apilamiento.lecturas_descartadas, 2);
    assert!(apilamiento.indels.is_empty());
    assert!(apilar(texto.as_bytes(), b"", &umbrales).is_ok());
    let llamadas = llamar_snps(&apilamiento, referencia, &umbrales);
    assert!(llamadas.is_empty());
    let llamadas = llamar_indels(&apilamiento, &umbrales);
    assert!(llamadas.is_empty());
}

#[test]
fn calidad_phred_es_la_cola_binomial() {
    assert!((calidad_phred(1, 1, 0.1) - 10.0).abs() < 1e-9);
    assert!((calidad_phred(2, 2, 0.1) - 20.0).abs() < 1e-9);
    assert!((calidad_phred(1, 2, 0.5) - 1.249_387_366).abs() < 1e-6);
    assert!((calidad_phred(3, 10, 0.01) - 39.436_703_3).abs() < 1e-5);
    assert_eq!(calidad_phred(0, 5, 0.1), 0.0);
    assert_eq!(calidad_phred(3, 3, 1.0), 0.0);
    assert!(calidad_phred(5000, 5000, 1e-4).is_finite());
    assert!((calidad_phred(98, 100, 0.01) - 1_923.140_357_9).abs() < 1e-6);
    assert!((calidad_phred(60, 100, 0.3) - 92.898_868_834).abs() < 1e-6);
}

const DIEZ: &str = "AAAAACCCCCGGGGGTTTTT";

#[test]
fn llama_un_snp_justo_en_los_umbrales() {
    let referencia = DIEZ.as_bytes();
    let mut lecturas = vec![(6, "10M", "CCCCCAGGGG"); 4];
    lecturas.push((6, "10M", "CCCCCGGGGG"));
    let sam = sam_de(20, &lecturas);
    let umbrales = Umbrales::default();
    let apilamiento = apilar(&sam[..], referencia, &umbrales).unwrap();
    let llamadas = llamar_snps(&apilamiento, referencia, &umbrales);
    assert_eq!(llamadas.len(), 1);
    let llamada = &llamadas[0];
    assert_eq!(llamada.cambio, cambio(10, "G", "A"));
    assert_eq!(llamada.profundidad, 5);
    assert_eq!(llamada.a_favor, 4);
    assert_eq!(llamada.fraccion, 0.8);
    assert!(
        (llamada.calidad - 172.095_265_96).abs() < 0.01,
        "{}",
        llamada.calidad
    );
}

#[test]
fn no_llama_snps_con_poca_fraccion_o_poca_profundidad() {
    let referencia = DIEZ.as_bytes();
    let umbrales = Umbrales::default();
    let mut lecturas = vec![(6, "10M", "CCCCCAGGGG"); 3];
    lecturas.extend(vec![(6, "10M", "CCCCCGGGGG"); 2]);
    let apilamiento = apilar(&sam_de(20, &lecturas)[..], referencia, &umbrales).unwrap();
    assert!(llamar_snps(&apilamiento, referencia, &umbrales).is_empty());
    let pocas = vec![(6, "10M", "CCCCCAGGGG"); 4];
    let apilamiento = apilar(&sam_de(20, &pocas)[..], referencia, &umbrales).unwrap();
    assert!(llamar_snps(&apilamiento, referencia, &umbrales).is_empty());
}

#[test]
fn un_borrado_cuenta_en_la_profundidad_del_snp() {
    let referencia = DIEZ.as_bytes();
    let umbrales = Umbrales::default();
    let mut lecturas = vec![(6, "10M", "CCCCCAGGGG"); 4];
    lecturas.push((6, "5M1D4M", "CCCCCGGGG"));
    let apilamiento = apilar(&sam_de(20, &lecturas)[..], referencia, &umbrales).unwrap();
    let llamadas = llamar_snps(&apilamiento, referencia, &umbrales);
    assert_eq!(llamadas.len(), 1);
    assert_eq!(llamadas[0].profundidad, 5);
}

#[test]
fn empate_entre_dos_alternativas_elige_la_primera_en_acgt() {
    let referencia = DIEZ.as_bytes();
    let umbrales = Umbrales {
        profundidad_minima: 4,
        fraccion_minima: 0.5,
        ..Umbrales::default()
    };
    let mut lecturas = vec![(6, "10M", "CCCCCTGGGG"); 2];
    lecturas.extend(vec![(6, "10M", "CCCCCAGGGG"); 2]);
    let apilamiento = apilar(&sam_de(20, &lecturas)[..], referencia, &umbrales).unwrap();
    let llamadas = llamar_snps(&apilamiento, referencia, &umbrales);
    assert_eq!(llamadas.len(), 1);
    assert_eq!(llamadas[0].cambio, cambio(10, "G", "A"));
    assert_eq!(llamadas[0].fraccion, 0.5);
}

#[test]
fn la_calidad_de_base_y_de_mapeo_filtran() {
    let referencia = DIEZ.as_bytes();
    let mut texto = String::new();
    for numero in 0..5 {
        texto.push_str(&format!(
            "b{numero}\t0\tchr\t6\t60\t10M\t*\t0\t0\tCCCCCAGGGG\tIIIII#IIII\n"
        ));
    }
    let umbrales = Umbrales::default();
    let apilamiento = apilar(texto.as_bytes(), referencia, &umbrales).unwrap();
    assert!(llamar_snps(&apilamiento, referencia, &umbrales).is_empty());
    let permisivos = Umbrales {
        calidad_base_minima: 0,
        ..Umbrales::default()
    };
    let apilamiento = apilar(texto.as_bytes(), referencia, &permisivos).unwrap();
    assert_eq!(llamar_snps(&apilamiento, referencia, &permisivos).len(), 1);
    let exigentes = Umbrales {
        calidad_mapeo_minima: 61,
        calidad_base_minima: 0,
        ..Umbrales::default()
    };
    let apilamiento = apilar(texto.as_bytes(), referencia, &exigentes).unwrap();
    assert_eq!(apilamiento.lecturas_usadas, 0);
}

const QUINCE: &str = "AAAAACCCCCGGGGG";

#[test]
fn llama_una_delecion_con_la_mitad_de_las_lecturas() {
    let referencia = DIEZ.as_bytes();
    let umbrales = Umbrales::default();
    let mut lecturas = vec![(1, "12M1D2M", "AAAAACCCCCGGGG"); 2];
    lecturas.extend(vec![(1, "10M1D4M", "AAAAACCCCCGGGG"); 2]);
    lecturas.extend(vec![(1, "15M", QUINCE); 4]);
    let apilamiento = apilar(&sam_de(20, &lecturas)[..], referencia, &umbrales).unwrap();
    assert_eq!(apilamiento.indels[&cambio(9, "CG", "C")], 4);
    let llamadas = llamar_indels(&apilamiento, &umbrales);
    assert_eq!(llamadas.len(), 1);
    assert_eq!(llamadas[0].cambio, cambio(9, "CG", "C"));
    assert_eq!(llamadas[0].a_favor, 4);
    assert_eq!(llamadas[0].profundidad, 8);
    assert_eq!(llamadas[0].fraccion, 0.5);
    assert!(llamadas[0].calidad > 0.0);
    assert!(llamar_snps(&apilamiento, referencia, &umbrales).is_empty());
}

#[test]
fn no_llama_una_delecion_minoritaria() {
    let referencia = DIEZ.as_bytes();
    let umbrales = Umbrales::default();
    let mut lecturas = vec![(1, "12M1D2M", "AAAAACCCCCGGGG"); 3];
    lecturas.extend(vec![(1, "15M", QUINCE); 5]);
    let apilamiento = apilar(&sam_de(20, &lecturas)[..], referencia, &umbrales).unwrap();
    assert!(llamar_indels(&apilamiento, &umbrales).is_empty());
}

#[test]
fn llama_una_insercion_alineada_a_la_izquierda() {
    let referencia = DIEZ.as_bytes();
    let umbrales = Umbrales::default();
    let mut lecturas = vec![(1, "8M1I7M", "AAAAACCCCCCGGGGG"); 3];
    lecturas.extend(vec![(1, "6M1I9M", "AAAAACCCCCCGGGGG"); 2]);
    lecturas.push((1, "15M", QUINCE));
    let apilamiento = apilar(&sam_de(20, &lecturas)[..], referencia, &umbrales).unwrap();
    let llamadas = llamar_indels(&apilamiento, &umbrales);
    assert_eq!(llamadas.len(), 1);
    assert_eq!(llamadas[0].cambio, cambio(4, "A", "AC"));
    assert_eq!(llamadas[0].a_favor, 5);
    assert_eq!(llamadas[0].profundidad, 6);
}

#[test]
fn dos_indels_empatados_en_la_misma_posicion_dan_una_sola_llamada() {
    let referencia = DIEZ.as_bytes();
    let umbrales = Umbrales::default();
    let mut lecturas = vec![(1, "5M1I10M", "AAAAATCCCCCGGGGG"); 3];
    lecturas.extend(vec![(1, "5M1I10M", "AAAAAGCCCCCGGGGG"); 3]);
    let apilamiento = apilar(&sam_de(20, &lecturas)[..], referencia, &umbrales).unwrap();
    let llamadas = llamar_indels(&apilamiento, &umbrales);
    assert_eq!(llamadas.len(), 1);
    assert_eq!(llamadas[0].cambio, cambio(4, "A", "AG"));
}

#[test]
fn normaliza_indels_hacia_la_izquierda() {
    let referencia = b"TTGAAAACTT";
    assert_eq!(
        normalizar(&cambio(5, "AA", "A"), referencia),
        cambio(2, "GA", "G")
    );
    assert_eq!(
        normalizar(&cambio(4, "A", "AA"), referencia),
        cambio(2, "G", "GA")
    );
    assert_eq!(
        normalizar(&cambio(2, "GAA", "GTA"), referencia),
        cambio(3, "A", "T")
    );
    assert_eq!(
        normalizar(&cambio(7, "C", "T"), referencia),
        cambio(7, "C", "T")
    );
    assert_eq!(
        normalizar(&cambio(7, "G", "T"), referencia),
        cambio(7, "G", "T")
    );
    assert_eq!(
        normalizar(&cambio(9, "TA", "T"), referencia),
        cambio(9, "TA", "T")
    );
    assert_eq!(
        normalizar(&cambio(0, "T", "T"), referencia),
        cambio(0, "T", "T")
    );
    let repetido = b"GCACACAT";
    assert_eq!(
        normalizar(&cambio(4, "ACA", "A"), repetido),
        cambio(0, "GCA", "G")
    );
    assert_eq!(
        normalizar(&cambio(0, "TT", "T"), referencia),
        cambio(0, "TT", "T")
    );
}

fn gen_de_prueba(nombre: &str, inicio: usize, fin: usize, hebra: u8, proteina: &str) -> Gen {
    Gen {
        locus: format!("p_{nombre}"),
        nombre: nombre.to_string(),
        inicio,
        fin,
        hebra,
        producto: String::new(),
        proteina: proteina.as_bytes().to_vec(),
    }
}

const GENOMA_DE_PRUEBA: &str = "TTTTATGAAATGGGGCTAACCCCCTTACTGCATGGGG";

fn genes_de_prueba() -> Vec<Gen> {
    vec![
        gen_de_prueba("uno", 5, 19, b'+', "MKWG"),
        gen_de_prueba("dos", 25, 33, b'-', "MQ"),
    ]
}

fn anotacion_unica(cambio_pedido: Cambio) -> Anotacion {
    let genes = genes_de_prueba();
    let referencia = GENOMA_DE_PRUEBA.as_bytes();
    let indice = IndiceDeGenes::nuevo(&genes, referencia.len());
    let mut anotaciones = anotar(&cambio_pedido, referencia, &indice);
    assert_eq!(anotaciones.len(), 1, "{cambio_pedido:?}");
    anotaciones.remove(0)
}

fn revisar(cambio_pedido: Cambio, nombre_gen: &str, efecto: Efecto, proteina: &str) {
    let anotacion = anotacion_unica(cambio_pedido.clone());
    assert_eq!(anotacion.nombre_gen, nombre_gen, "{cambio_pedido:?}");
    assert_eq!(anotacion.efecto, efecto, "{cambio_pedido:?}");
    assert_eq!(anotacion.proteina, proteina, "{cambio_pedido:?}");
    assert_eq!(anotacion.resistencia, None);
}

#[test]
fn anota_snps_en_un_gen_de_la_hebra_directa() {
    revisar(cambio(7, "A", "T"), "uno", Efecto::ParadaPrematura, "K2*");
    revisar(cambio(9, "A", "G"), "uno", Efecto::Sinonima, "K2K");
    revisar(
        cambio(10, "T", "C"),
        "uno",
        Efecto::CambioDeAminoacido,
        "W3R",
    );
    revisar(cambio(4, "A", "C"), "uno", Efecto::PerdidaDeInicio, "M1L");
    revisar(cambio(4, "A", "G"), "uno", Efecto::Sinonima, "M1M");
    revisar(cambio(16, "T", "C"), "uno", Efecto::PerdidaDeParada, "*5Q");
    revisar(cambio(17, "A", "G"), "uno", Efecto::Sinonima, "*5*");
}

#[test]
fn anota_indels_en_un_gen_de_la_hebra_directa() {
    revisar(
        cambio(9, "AT", "A"),
        "uno",
        Efecto::CorrimientoDeMarco,
        "W3fs",
    );
    revisar(
        cambio(9, "ATGG", "A"),
        "uno",
        Efecto::DelecionEnMarco,
        "W3del",
    );
    revisar(
        cambio(10, "TGGG", "T"),
        "uno",
        Efecto::DelecionEnMarco,
        "W3_G4del",
    );
    revisar(
        cambio(12, "G", "GAAA"),
        "uno",
        Efecto::InsercionEnMarco,
        "W3_G4ins",
    );
    revisar(
        cambio(12, "G", "GA"),
        "uno",
        Efecto::CorrimientoDeMarco,
        "G4fs",
    );
    revisar(cambio(3, "TA", "T"), "uno", Efecto::PerdidaDeInicio, "M1?");
    revisar(cambio(4, "A", "AC"), "uno", Efecto::PerdidaDeInicio, "M1?");
}

#[test]
fn anota_cambios_en_un_gen_de_la_hebra_complementaria() {
    revisar(cambio(29, "G", "A"), "dos", Efecto::ParadaPrematura, "Q2*");
    revisar(cambio(27, "C", "T"), "dos", Efecto::Sinonima, "Q2Q");
    revisar(cambio(32, "T", "C"), "dos", Efecto::Sinonima, "M1M");
    revisar(cambio(31, "A", "G"), "dos", Efecto::PerdidaDeInicio, "M1T");
    revisar(
        cambio(27, "CT", "C"),
        "dos",
        Efecto::CorrimientoDeMarco,
        "Q2fs",
    );
    revisar(
        cambio(27, "C", "CA"),
        "dos",
        Efecto::CorrimientoDeMarco,
        "Q2fs",
    );
    revisar(
        cambio(26, "ACTG", "A"),
        "dos",
        Efecto::DelecionEnMarco,
        "Q2del",
    );
}

fn vecino(nombre: &str, distancia: usize, corriente_arriba: bool) -> Option<Vecino> {
    Some(Vecino {
        nombre: nombre.to_string(),
        distancia,
        corriente_arriba,
    })
}

#[test]
fn anota_variantes_intergenicas_con_sus_vecinos() {
    let entre = anotacion_unica(cambio(21, "C", "A"));
    assert_eq!(entre.nombre_gen, "");
    assert_eq!(
        entre.efecto,
        Efecto::Intergenica {
            izquierda: vecino("uno", 3, false),
            derecha: vecino("dos", 3, false),
        }
    );
    assert_eq!(
        anotacion_unica(cambio(1, "T", "G")).efecto,
        Efecto::Intergenica {
            izquierda: None,
            derecha: vecino("uno", 3, true),
        }
    );
    assert_eq!(
        anotacion_unica(cambio(35, "G", "A")).efecto,
        Efecto::Intergenica {
            izquierda: vecino("dos", 3, true),
            derecha: None,
        }
    );
    assert_eq!(
        anotacion_unica(cambio(3, "T", "TA")).efecto,
        Efecto::Intergenica {
            izquierda: None,
            derecha: vecino("uno", 1, true),
        }
    );
}

#[test]
fn genes_superpuestos_no_codificantes_y_sin_genes() {
    let referencia = GENOMA_DE_PRUEBA.as_bytes();
    let mut genes = genes_de_prueba();
    genes.push(gen_de_prueba("arn", 6, 12, b'+', ""));
    genes.push(gen_de_prueba("fuera", 30, 80, b'+', "M"));
    let indice = IndiceDeGenes::nuevo(&genes, referencia.len());
    let anotaciones = anotar(&cambio(10, "T", "C"), referencia, &indice);
    assert_eq!(anotaciones.len(), 2);
    assert_eq!(anotaciones[0].nombre_gen, "uno");
    assert_eq!(anotaciones[1].nombre_gen, "arn");
    assert_eq!(anotaciones[1].efecto, Efecto::NoCodificante);
    let vacio = IndiceDeGenes::nuevo(&[], referencia.len());
    assert!(anotar(&cambio(10, "T", "C"), referencia, &vacio).is_empty());
}

#[test]
fn reconoce_mutaciones_de_resistencia_en_e_coli() {
    let genoma = ecoli::genoma();
    let indice = IndiceDeGenes::nuevo(ecoli::genes(), genoma.len());
    let mutaciones = mutaciones_de_resistencia().unwrap();
    let esperado = [
        ("rpoB", "S531F", "rifampicina"),
        ("gyrA", "S83L", "quinolonas"),
    ];
    for (mutacion, (nombre, proteina, antibiotico)) in mutaciones.iter().zip(esperado) {
        let pedido = Cambio {
            posicion: mutacion.posicion,
            referencia: mutacion.referencia.clone(),
            alternativa: mutacion.alternativa.clone(),
        };
        let anotaciones = anotar(&pedido, genoma, &indice);
        assert_eq!(anotaciones.len(), 1);
        assert_eq!(anotaciones[0].nombre_gen, nombre);
        assert_eq!(anotaciones[0].efecto, Efecto::CambioDeAminoacido);
        assert_eq!(anotaciones[0].proteina, proteina);
        assert_eq!(anotaciones[0].resistencia, Some(antibiotico));
    }
    let rpob = ecoli::buscar_gen("rpoB").unwrap();
    let primera_del_codon_600 = rpob.inicio - 1 + 599 * 3;
    let base = genoma[primera_del_codon_600];
    let otra = if base == b'A' { b'C' } else { b'A' };
    let lejos = Cambio {
        posicion: primera_del_codon_600,
        referencia: vec![base],
        alternativa: vec![otra],
    };
    let anotaciones = anotar(&lejos, genoma, &indice);
    assert_eq!(anotaciones[0].nombre_gen, "rpoB");
    assert_eq!(anotaciones[0].resistencia, None);
}

#[test]
fn escribe_y_lee_vcf() {
    let llamada = Llamada {
        cambio: cambio(3, "A", "AC"),
        profundidad: 12,
        a_favor: 11,
        fraccion: 11.0 / 12.0,
        calidad: 87.25,
        anotaciones: Vec::new(),
    };
    let mut destino = Vec::new();
    escribir_llamadas_vcf(&mut destino, "chr", 100, &[llamada]).unwrap();
    let texto = String::from_utf8(destino).unwrap();
    assert!(texto.starts_with("##fileformat=VCFv4.2\n"));
    assert!(texto.contains("##contig=<ID=chr,length=100>"));
    assert!(
        texto.ends_with("chr\t4\t.\tA\tAC\t87.2\tPASS\tDP=12;AF=0.917\n")
            || texto.ends_with("chr\t4\t.\tA\tAC\t87.3\tPASS\tDP=12;AF=0.917\n"),
        "{texto}"
    );
    let leidos = leer_cambios_vcf(texto.as_bytes()).unwrap();
    assert_eq!(leidos, vec![cambio(3, "A", "AC")]);
    let varios = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr\t10\t.\tA\tC,T\t50\tPASS\t.\nchr\t20\t.\tG\t<DEL>\t50\tPASS\t.\nchr\t30\t.\tG\t.\t50\tPASS\t.\n";
    assert_eq!(
        leer_cambios_vcf(varios.as_bytes()).unwrap(),
        vec![cambio(9, "A", "C"), cambio(9, "A", "T")]
    );
    assert!(leer_cambios_vcf(&b"chr\tx\t.\tA\tC\t1\tPASS\t.\n"[..]).is_err());
    assert!(leer_cambios_vcf(&b"chr\t0\t.\tA\tC\t1\tPASS\t.\n"[..]).is_err());
}

#[test]
fn compara_con_la_verdad_despues_de_normalizar() {
    let referencia = b"TTGAAAACTT";
    let llamadas = vec![cambio(2, "GA", "G"), cambio(8, "T", "C")];
    let verdad = vec![cambio(5, "AA", "A"), cambio(0, "T", "G")];
    let comparacion = comparar_con_verdad(&llamadas, &verdad, referencia);
    assert_eq!(comparacion.verdaderos_positivos, vec![cambio(2, "GA", "G")]);
    assert_eq!(comparacion.falsos_positivos, vec![cambio(8, "T", "C")]);
    assert_eq!(comparacion.falsos_negativos, vec![cambio(0, "T", "G")]);
    assert_eq!(comparacion.sensibilidad(), Some(0.5));
    assert_eq!(comparacion.precision(), Some(0.5));
    let vacia = comparar_con_verdad(&[], &[], referencia);
    assert_eq!(vacia.sensibilidad(), None);
    assert_eq!(vacia.precision(), None);
}

#[test]
fn clasifica_transiciones_y_tipos() {
    assert!(es_transicion(b'A', b'G'));
    assert!(es_transicion(b'T', b'C'));
    assert!(!es_transicion(b'A', b'C'));
    assert!(!es_transicion(b'G', b'T'));
    assert!(cambio(1, "A", "G").es_snp());
    assert!(cambio(1, "A", "AG").es_insercion());
    assert!(cambio(1, "AG", "A").es_delecion());
    assert!(!cambio(1, "AG", "A").es_snp());
}

fn referencia_al_azar(longitud: usize, semilla: u64) -> Vec<u8> {
    let mut azar = Azar::nuevo(semilla);
    (0..longitud).map(|_| b"ACGT"[azar.entero(4)]).collect()
}

fn como_cambio(mutacion: &Mutacion) -> Cambio {
    Cambio {
        posicion: mutacion.posicion,
        referencia: mutacion.referencia.clone(),
        alternativa: mutacion.alternativa.clone(),
    }
}

#[test]
fn encuentra_mutaciones_simuladas_desde_un_sam_comprimido() {
    let referencia = referencia_al_azar(30_000, 11);
    let mut mutaciones = Vec::new();
    for (numero, posicion) in [3_000, 9_000, 15_000, 21_000].into_iter().enumerate() {
        let original = referencia[posicion];
        let nueva =
            b"ACGT"[(nucleo::secuencia::indice_base(original).unwrap() + 1 + numero % 3) % 4];
        mutaciones.push(Mutacion {
            posicion,
            referencia: vec![original],
            alternativa: vec![nueva],
            nombre: "SNP".to_string(),
        });
    }
    mutaciones.push(Mutacion {
        posicion: 12_000,
        referencia: referencia[12_000..12_003].to_vec(),
        alternativa: vec![referencia[12_000]],
        nombre: "deleción".to_string(),
    });
    mutaciones.push(Mutacion {
        posicion: 25_000,
        referencia: vec![referencia[25_000]],
        alternativa: vec![referencia[25_000], b'T', b'T'],
        nombre: "inserción".to_string(),
    });
    let opciones = OpcionesSimulacion {
        cobertura: 20.0,
        fraccion_con_adaptador: 0.05,
        semilla: 31,
        ..OpcionesSimulacion::default()
    };
    let muestra = simular_lecturas(&referencia, "azar", &mutaciones, &opciones).unwrap();
    let ruta = archivo_temporal("comprimido.sam.gz");
    let mut texto = Vec::new();
    escribir_cabecera_sam(&mut texto, "azar", referencia.len()).unwrap();
    for registro in muestra.verdad.iter().rev() {
        escribir_sam(&mut texto, registro).unwrap();
    }
    let mut comprimido = GzEncoder::new(Vec::new(), Compression::fast());
    comprimido.write_all(&texto).unwrap();
    std::fs::write(&ruta, comprimido.finish().unwrap()).unwrap();
    let informe = llamar_variantes(&ruta, &referencia, &[], &Umbrales::default()).unwrap();
    std::fs::remove_file(&ruta).unwrap();
    assert_eq!(informe.nombre_referencia, "azar");
    assert_eq!(informe.lecturas_usadas, muestra.verdad.len());
    let llamadas: Vec<Cambio> = informe.llamadas.iter().map(|l| l.cambio.clone()).collect();
    let verdad: Vec<Cambio> = muestra.mutaciones.iter().map(como_cambio).collect();
    let comparacion = comparar_con_verdad(&llamadas, &verdad, &referencia);
    assert!(
        comparacion.falsos_positivos.is_empty(),
        "{:?}",
        comparacion.falsos_positivos
    );
    assert!(
        comparacion.falsos_negativos.is_empty(),
        "{:?}",
        comparacion.falsos_negativos
    );
    assert!(informe.llamadas.iter().all(|l| l.anotaciones.is_empty()));
    assert!(informe.llamadas.iter().all(|l| l.calidad > 50.0));
    assert!(informe.profundidad.suficientes > 29_000);
    assert!(informe.profundidad.media > 15.0 && informe.profundidad.media < 25.0);
}

#[test]
fn rechaza_umbrales_y_referencias_imposibles() {
    let ruta = archivo_temporal("vacio.sam");
    std::fs::write(&ruta, "@SQ\tSN:chr\tLN:4\n").unwrap();
    let malos = [
        Umbrales {
            profundidad_minima: 0,
            ..Umbrales::default()
        },
        Umbrales {
            fraccion_minima: 0.0,
            ..Umbrales::default()
        },
        Umbrales {
            fraccion_minima: 1.5,
            ..Umbrales::default()
        },
        Umbrales {
            fraccion_minima_indel: f64::NAN,
            ..Umbrales::default()
        },
    ];
    for umbrales in malos {
        assert!(llamar_variantes(&ruta, b"ACGT", &[], &umbrales).is_err());
    }
    assert!(llamar_variantes(&ruta, b"", &[], &Umbrales::default()).is_err());
    let informe = llamar_variantes(&ruta, b"ACGT", &[], &Umbrales::default()).unwrap();
    assert!(informe.llamadas.is_empty());
    assert_eq!(informe.profundidad.con_lecturas, 0);
    std::fs::remove_file(&ruta).unwrap();
    assert!(llamar_variantes(&ruta, b"ACGT", &[], &Umbrales::default()).is_err());
}

#[test]
fn muestra_de_ejemplo_completa() {
    let muestra = muestra_de_ejemplo().unwrap();
    let genoma = ecoli::genoma();
    let ruta = archivo_temporal("ejemplo.sam");
    escribir_verdad_sam(&ruta, ecoli::IDENTIFICADOR, genoma.len(), &muestra.verdad);
    let informe = llamar_variantes(&ruta, genoma, ecoli::genes(), &Umbrales::default()).unwrap();
    std::fs::remove_file(&ruta).unwrap();
    let llamadas: Vec<Cambio> = informe.llamadas.iter().map(|l| l.cambio.clone()).collect();
    let verdad: Vec<Cambio> = muestra.mutaciones.iter().map(como_cambio).collect();
    let comparacion = comparar_con_verdad(&llamadas, &verdad, genoma);
    let snps_verdaderos = comparacion
        .verdaderos_positivos
        .iter()
        .filter(|c| c.es_snp())
        .count();
    let indels_verdaderos = comparacion.verdaderos_positivos.len() - snps_verdaderos;
    assert_eq!(snps_verdaderos, 22);
    assert!(
        comparacion.falsos_positivos.is_empty(),
        "{:?}",
        comparacion.falsos_positivos
    );
    assert!(indels_verdaderos >= 3);
    let proteinas: Vec<(String, String, Option<&str>)> = informe
        .llamadas
        .iter()
        .flat_map(|l| l.anotaciones.iter())
        .map(|a| (a.nombre_gen.clone(), a.proteina.clone(), a.resistencia))
        .collect();
    assert!(proteinas.contains(&("rpoB".to_string(), "S531F".to_string(), Some("rifampicina"))));
    assert!(proteinas.contains(&("gyrA".to_string(), "S83L".to_string(), Some("quinolonas"))));
    let mut vcf = Vec::new();
    escribir_llamadas_vcf(
        &mut vcf,
        &informe.nombre_referencia,
        genoma.len(),
        &informe.llamadas,
    )
    .unwrap();
    let releidos = leer_cambios_vcf(&vcf[..]).unwrap();
    assert_eq!(releidos, llamadas);
}

#[test]
fn un_indel_que_deja_el_codon_de_inicio_intacto_no_pierde_el_inicio() {
    revisar(
        cambio(5, "T", "TG"),
        "uno",
        Efecto::CorrimientoDeMarco,
        "K2fs",
    );
    let genes = vec![gen_de_prueba("repetido", 2, 13, b'+', "MMK")];
    let referencia = b"CATGATGAAATAAC";
    let indice = IndiceDeGenes::nuevo(&genes, referencia.len());
    let anotaciones = anotar(&cambio(0, "CATG", "C"), referencia, &indice);
    assert_eq!(anotaciones.len(), 1);
    assert_eq!(anotaciones[0].efecto, Efecto::DelecionEnMarco);
    assert_eq!(anotaciones[0].proteina, "M1del");
    let anotaciones = anotar(&cambio(0, "CATGA", "C"), referencia, &indice);
    assert_eq!(anotaciones[0].efecto, Efecto::PerdidaDeInicio);
}

#[test]
fn normaliza_indels_en_el_borde_izquierdo_con_ancla_a_la_derecha() {
    let referencia = b"AAAGTT";
    let borrado = cambio(0, "AA", "A");
    assert_eq!(normalizar(&cambio(1, "AA", "A"), referencia), borrado);
    assert_eq!(
        normalizar(&cambio(2, "AG", "A"), referencia),
        cambio(2, "AG", "A")
    );
    assert_eq!(normalizar(&cambio(0, "AAA", "AA"), referencia), borrado);
    let insertado = cambio(0, "A", "AA");
    assert_eq!(normalizar(&cambio(2, "A", "AA"), referencia), insertado);
    let comparacion = comparar_con_verdad(
        &[cambio(1, "AA", "A")],
        &[cambio(0, "AAA", "AA")],
        referencia,
    );
    assert_eq!(comparacion.verdaderos_positivos, vec![borrado]);
    assert!(comparacion.falsos_positivos.is_empty());
    assert!(comparacion.falsos_negativos.is_empty());
}

#[test]
fn las_llamadas_en_un_borde_homopolimerico_se_suman_en_un_solo_indel() {
    let referencia = b"AAAAGCCCCTTTTGGGGACGT";
    let umbrales = Umbrales::default();
    let mut lecturas = vec![(1, "2M1D18M", "AAAGCCCCTTTTGGGGACGT"); 2];
    lecturas.extend(vec![(1, "3M1D17M", "AAAGCCCCTTTTGGGGACGT"); 2]);
    lecturas.push((1, "1M1D19M", "AAAGCCCCTTTTGGGGACGT"));
    let apilamiento = apilar(&sam_de(21, &lecturas)[..], referencia, &umbrales).unwrap();
    let llamadas = llamar_indels(&apilamiento, &umbrales);
    assert_eq!(llamadas.len(), 1);
    assert_eq!(llamadas[0].cambio, cambio(0, "AA", "A"));
    assert_eq!(llamadas[0].a_favor, 5);
}

#[test]
fn rechaza_una_referencia_mas_grande_que_el_limite() {
    let ruta = archivo_temporal("grande.sam");
    std::fs::write(&ruta, "@SQ\tSN:chr\tLN:4\n").unwrap();
    let enorme = vec![b'A'; 16_000_001];
    let error = llamar_variantes(&ruta, &enorme, &[], &Umbrales::default()).unwrap_err();
    assert!(error.contains("16 millones"), "{error}");
    std::fs::remove_file(&ruta).unwrap();
}

#[test]
fn una_referencia_en_minusculas_da_las_mismas_llamadas() {
    let ruta = archivo_temporal("minusculas.sam");
    let referencia = b"TTTACGGGGTTT";
    let mut lecturas = vec![(1, "4M1I8M", "TTTACGGGGTTT"); 0];
    lecturas.extend(vec![(1, "4M1I8M", "TTTAACGGGGTTT"); 6]);
    let texto = sam_de(12, &lecturas);
    std::fs::write(&ruta, &texto).unwrap();
    let en_mayusculas = llamar_variantes(&ruta, referencia, &[], &Umbrales::default()).unwrap();
    let en_minusculas = llamar_variantes(
        &ruta,
        &referencia.to_ascii_lowercase(),
        &[],
        &Umbrales::default(),
    )
    .unwrap();
    std::fs::remove_file(&ruta).unwrap();
    assert_eq!(en_mayusculas.llamadas.len(), 1);
    assert_eq!(en_mayusculas.llamadas[0].cambio, cambio(2, "T", "TA"));
    assert_eq!(en_minusculas.llamadas, en_mayusculas.llamadas);
}

#[test]
fn el_vcf_de_verdad_ignora_alelos_iguales_a_la_referencia_o_con_letras_raras() {
    let texto = "chr\t5\t.\tA\tA\t1\tPASS\t.\nchr\t6\t.\tA\tN\t1\tPASS\t.\nchr\t7\t.\tA\tR,C\t1\tPASS\t.\nchr\t8\t.\ta\tc\t1\tPASS\t.\n";
    assert_eq!(
        leer_cambios_vcf(texto.as_bytes()).unwrap(),
        vec![cambio(6, "A", "C"), cambio(7, "A", "C")]
    );
}
