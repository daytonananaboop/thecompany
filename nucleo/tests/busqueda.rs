use nucleo::azar::Azar;
use nucleo::busqueda::*;
use nucleo::ecoli;
use nucleo::formatos::Registro;
use nucleo::matrices::{AMINOACIDOS, blosum62};
use nucleo::secuencia::complemento_inverso;
use std::sync::OnceLock;

fn proteinas() -> &'static BaseDeProteinas {
    static BASE: OnceLock<BaseDeProteinas> = OnceLock::new();
    BASE.get_or_init(proteinas_de_ecoli)
}

fn genoma() -> &'static BaseDeAdn {
    static BASE: OnceLock<BaseDeAdn> = OnceLock::new();
    BASE.get_or_init(genoma_de_ecoli)
}

fn indice_en_base(nombre: &str) -> usize {
    proteinas()
        .registros
        .iter()
        .position(|registro| registro.nombre.split_whitespace().next() == Some(nombre))
        .unwrap()
}

fn posiciones_al_azar(longitud: usize, cuantas: usize, azar: &mut Azar) -> Vec<usize> {
    let mut posiciones: Vec<usize> = (0..longitud).collect();
    for i in 0..cuantas {
        let j = i + azar.entero(longitud - i);
        posiciones.swap(i, j);
    }
    posiciones.truncate(cuantas);
    posiciones
}

fn mutar(secuencia: &[u8], alfabeto: &[u8], fraccion: f64, semilla: u64) -> Vec<u8> {
    let mut azar = Azar::nuevo(semilla);
    let mut mutada = secuencia.to_vec();
    let cuantas = (secuencia.len() as f64 * fraccion).round() as usize;
    for posicion in posiciones_al_azar(secuencia.len(), cuantas, &mut azar) {
        loop {
            let nueva = alfabeto[azar.entero(alfabeto.len())];
            if nueva != mutada[posicion] {
                mutada[posicion] = nueva;
                break;
            }
        }
    }
    mutada
}

fn al_azar(alfabeto: &[u8], longitud: usize, semilla: u64) -> Vec<u8> {
    let mut azar = Azar::nuevo(semilla);
    (0..longitud)
        .map(|_| alfabeto[azar.entero(alfabeto.len())])
        .collect()
}

fn palabra(codigo: u16) -> Vec<u8> {
    let codigo = codigo as usize;
    vec![
        AMINOACIDOS[codigo / 400],
        AMINOACIDOS[codigo / 20 % 20],
        AMINOACIDOS[codigo % 20],
    ]
}

#[test]
fn puntajes_de_adn_y_de_proteina() {
    let adn = Puntuacion::adn(2, -3);
    assert_eq!(adn.puntaje(b'A', b'A'), 2);
    assert_eq!(adn.puntaje(b'G', b'G'), 2);
    assert_eq!(adn.puntaje(b'A', b'C'), -3);
    assert_eq!(adn.puntaje(b'N', b'N'), -3);
    let proteina = Puntuacion::blosum62();
    assert_eq!(proteina.puntaje(b'W', b'W'), 11);
    assert_eq!(proteina.puntaje(b'A', b'R'), -1);
    assert_eq!(proteina.puntaje(b'I', b'V'), 3);
    assert_eq!(proteina.puntaje(b'*', b'A'), -4);
}

#[test]
fn vecinas_de_aaa_y_de_ggg() {
    let de_aaa = palabras_vecinas(b"AAA", 11);
    assert_eq!(de_aaa.len(), 1);
    assert_eq!(palabra(de_aaa[0]), b"AAA");
    let de_ggg: Vec<Vec<u8>> = palabras_vecinas(b"GGG", 11)
        .into_iter()
        .map(palabra)
        .collect();
    assert_eq!(de_ggg.len(), 13);
    assert!(de_ggg.contains(&b"GGG".to_vec()));
    assert!(de_ggg.contains(&b"DGG".to_vec()));
    assert!(de_ggg.contains(&b"GGS".to_vec()));
    assert!(!de_ggg.contains(&b"GGR".to_vec()));
    assert!(palabras_vecinas(b"AA", 11).is_empty());
    assert!(palabras_vecinas(b"AXA", 11).is_empty());
}

#[test]
fn vecinas_coinciden_con_la_fuerza_bruta() {
    for consulta in [&b"WWW"[..], b"MKV", b"CHY", b"PEG", b"LLL"] {
        for umbral in [11, 13, 20] {
            let mut esperadas = Vec::new();
            for &x in AMINOACIDOS {
                for &y in AMINOACIDOS {
                    for &z in AMINOACIDOS {
                        let puntaje = blosum62(consulta[0], x)
                            + blosum62(consulta[1], y)
                            + blosum62(consulta[2], z);
                        if puntaje >= umbral {
                            esperadas.push(vec![x, y, z]);
                        }
                    }
                }
            }
            let mut obtenidas: Vec<Vec<u8>> = palabras_vecinas(consulta, umbral)
                .into_iter()
                .map(palabra)
                .collect();
            obtenidas.sort();
            esperadas.sort();
            assert_eq!(obtenidas, esperadas);
        }
    }
}

#[test]
fn codigo_de_palabra_es_la_base_veinte() {
    assert_eq!(codigo_de_palabra(b"AAA"), Some(0));
    assert_eq!(codigo_de_palabra(b"ARN"), Some(22));
    assert_eq!(codigo_de_palabra(b"VVV"), Some(7999));
    assert_eq!(codigo_de_palabra(b"AXA"), None);
}

#[test]
fn extension_sin_huecos_se_detiene_con_la_caida() {
    let tramo = extender_sin_huecos(
        b"TTACGTACGTTT",
        b"CCACGTACGTCC",
        2,
        2,
        8,
        &Puntuacion::adn(2, -3),
        5,
    );
    assert_eq!(
        tramo,
        Tramo {
            inicio_consulta: 2,
            inicio_sujeto: 2,
            longitud: 8,
            puntaje: 16
        }
    );
}

#[test]
fn extension_sin_huecos_cruza_un_desapareamiento() {
    let tramo = extender_sin_huecos(
        b"ACGTACGTTAAAA",
        b"ACGTACGTCAAAA",
        0,
        0,
        8,
        &Puntuacion::adn(2, -3),
        5,
    );
    assert_eq!(tramo.longitud, 13);
    assert_eq!(tramo.puntaje, 21);
    let corta = extender_sin_huecos(
        b"ACGTACGTTAAAA",
        b"ACGTACGTCAAAA",
        0,
        0,
        8,
        &Puntuacion::adn(2, -3),
        2,
    );
    assert_eq!(corta.longitud, 8);
    assert_eq!(corta.puntaje, 16);
}

#[test]
fn extension_sin_huecos_no_alarga_con_un_empate() {
    let tramo = extender_sin_huecos(b"WWWAS", b"WWWRT", 0, 0, 3, &Puntuacion::blosum62(), 16);
    assert_eq!(tramo.longitud, 3);
    assert_eq!(tramo.puntaje, 33);
    let hacia_atras = extender_sin_huecos(b"SAWWW", b"TRWWW", 2, 2, 3, &Puntuacion::blosum62(), 16);
    assert_eq!(hacia_atras.inicio_consulta, 2);
    assert_eq!(hacia_atras.longitud, 3);
}

#[test]
fn posiciones_fuera_de_rango_no_rompen_nada() {
    let adn = Puntuacion::adn(2, -3);
    let vacio = extender_sin_huecos(b"ACGT", b"ACGT", 3, 0, 3, &adn, 5);
    assert_eq!((vacio.longitud, vacio.puntaje), (0, 0));
    let afuera = alinear_con_huecos(b"ACGT", b"ACGT", 50, 50, &parametros_blastn(), 110);
    assert_eq!((afuera.inicio_consulta, afuera.fin_consulta), (0, 4));
    assert_eq!(afuera.puntaje, 8);
    let vacias = alinear_con_huecos(b"", b"", 0, 0, &parametros_blastp(), 64);
    assert_eq!(vacias.puntaje, 0);
    assert!(genoma().registro(7).is_empty());
}

#[test]
fn alineamiento_con_un_hueco() {
    let consulta = b"AAAACCCCGGGGTTTTAAAACCCC";
    let sujeto = b"AAAACCCCGGGTTTTAAAACCCC";
    let alineamiento = alinear_con_huecos(consulta, sujeto, 2, 2, &parametros_blastn(), 110);
    assert_eq!(alineamiento.puntaje, 23 * 2 - (5 + 2));
    assert_eq!(alineamiento.inicio_consulta, 0);
    assert_eq!(alineamiento.fin_consulta, 24);
    assert_eq!(alineamiento.inicio_sujeto, 0);
    assert_eq!(alineamiento.fin_sujeto, 23);
    let parejas = alineamiento
        .operaciones
        .iter()
        .filter(|&&operacion| operacion == PAREJA)
        .count();
    let huecos = alineamiento
        .operaciones
        .iter()
        .filter(|&&operacion| operacion == HUECO_EN_SUJETO)
        .count();
    assert_eq!((parejas, huecos), (23, 1));
}

#[test]
fn alineamiento_deja_fuera_los_bordes_distintos() {
    let consulta = b"TTTTTTAAAACCCCGGGGTTTTAAAACCCCTTTTTT";
    let sujeto = b"GGGGGGAAAACCCCGGGGTTTTAAAACCCCGGGGGG";
    let alineamiento = alinear_con_huecos(consulta, sujeto, 15, 15, &parametros_blastn(), 110);
    assert_eq!(alineamiento.puntaje, 48);
    assert_eq!(
        (alineamiento.inicio_consulta, alineamiento.fin_consulta),
        (6, 30)
    );
    assert_eq!(
        (alineamiento.inicio_sujeto, alineamiento.fin_sujeto),
        (6, 30)
    );
    assert!(
        alineamiento
            .operaciones
            .iter()
            .all(|&operacion| operacion == PAREJA)
    );
}

#[test]
fn alineamiento_con_hueco_largo_en_la_consulta() {
    let consulta = b"MKVLHEWTRAQFDNSGIPLCYMKVLHEWTRAQ";
    let sujeto = b"MKVLHEWTRAQFDNSGGGGGIPLCYMKVLHEWTRAQ";
    let alineamiento = alinear_con_huecos(consulta, sujeto, 3, 3, &parametros_blastp(), 64);
    let huecos = alineamiento
        .operaciones
        .iter()
        .filter(|&&operacion| operacion == HUECO_EN_CONSULTA)
        .count();
    assert_eq!(huecos, 4);
    assert_eq!(alineamiento.fin_consulta, consulta.len());
    assert_eq!(alineamiento.fin_sujeto, sujeto.len());
}

#[test]
fn bits_y_evalor_siguen_karlin_altschul() {
    let estadistica = parametros_blastp().estadistica;
    let calculados = bits(100, &estadistica);
    assert!((calculados - 43.128).abs() < 0.01, "{calculados}");
    let esperado = 1000.0 * 1e6 * 2f64.powf(-40.0);
    assert!((evalor(40.0, 1000.0 * 1e6) - esperado).abs() < 1e-12);
    assert!((evalor(40.0, 1e9) - 9.095e-4).abs() < 1e-6);
}

#[test]
fn puntaje_necesario_es_el_minimo_que_cumple_el_evalor() {
    let estadistica = parametros_blastp().estadistica;
    let minimo = puntaje_necesario(1e-3, 1e9, &estadistica);
    assert_eq!(minimo, 92);
    assert!(evalor(bits(minimo, &estadistica), 1e9) <= 1e-3);
    assert!(evalor(bits(minimo - 1, &estadistica), 1e9) > 1e-3);
    assert_eq!(puntaje_necesario(0.0, 1e9, &estadistica), i32::MAX);
}

#[test]
fn ajuste_de_longitud_y_sus_topes() {
    let estadistica = parametros_blastp().estadistica;
    let ajuste = ajuste_de_longitud(300.0, 1e6, 1000.0, &estadistica);
    assert!((ajuste - 116.6).abs() < 0.2, "{ajuste}");
    assert_eq!(ajuste_de_longitud(20.0, 1e6, 1000.0, &estadistica), 0.0);
    let corta = ajuste_de_longitud(300.0, 1000.0, 100.0, &estadistica);
    assert!(corta >= 0.0 && 1000.0 - 100.0 * corta >= 1.0 / estadistica.k - 1e-9);
    let espacio = espacio_de_busqueda(300.0, 1e6, 1000.0, &estadistica);
    assert!((espacio - (300.0 - ajuste) * (1e6 - 1000.0 * ajuste)).abs() < 1.0);
}

#[test]
fn indices_de_ecoli() {
    assert_eq!(proteinas().registros.len(), ecoli::genes().len());
    assert_eq!(proteinas().total_letras, 1_330_036);
    let parc = &proteinas().registros[indice_en_base("parC")];
    assert!(parc.nombre.contains("b3019"));
    assert!(parc.nombre.contains("topoisomerase IV"));
    assert_eq!(genoma().cantidad(), 1);
    assert_eq!(genoma().total_letras, ecoli::genoma().len());
    assert_eq!(genoma().nombres[0], ecoli::IDENTIFICADOR);
    assert_eq!(genoma().registro(0), ecoli::genoma());
}

#[test]
fn blastp_encuentra_parc_y_su_paralogo_gyra() {
    let parc = &ecoli::buscar_gen("parC").unwrap().proteina;
    let aciertos = blastp(parc, proteinas(), 1e-3, 20).unwrap();
    let primero = &aciertos[0];
    assert_eq!(primero.sujeto, indice_en_base("parC"));
    assert_eq!(primero.identidades, parc.len());
    assert_eq!(primero.longitud, parc.len());
    assert!(primero.evalor < 1e-100);
    assert_eq!(primero.consulta_alineada, primero.sujeto_alineado);
    let gyra = aciertos
        .iter()
        .find(|acierto| acierto.sujeto == indice_en_base("gyrA"))
        .expect("gyrA debe aparecer");
    let identidad = gyra.identidades as f64 / gyra.longitud as f64;
    assert!((0.25..0.45).contains(&identidad), "{identidad}");
    assert!(gyra.evalor < 1e-30);
    assert!(gyra.longitud > 300);
    for par in aciertos.windows(2) {
        assert!(par[0].menos_log10_evalor >= par[1].menos_log10_evalor);
    }
}

#[test]
fn blastp_fragmento_exacto_con_coordenadas() {
    let reca = &ecoli::buscar_gen("recA").unwrap().proteina;
    let aciertos = blastp(&reca[50..250], proteinas(), 1e-3, 5).unwrap();
    let primero = &aciertos[0];
    assert_eq!(primero.sujeto, indice_en_base("recA"));
    assert_eq!(primero.marco, 0);
    assert_eq!((primero.inicio_consulta, primero.fin_consulta), (0, 200));
    assert_eq!((primero.inicio_sujeto, primero.fin_sujeto), (50, 250));
    assert_eq!(primero.identidades, 200);
    assert_eq!(primero.desapareamientos, 0);
    assert_eq!(primero.aperturas, 0);
    assert_eq!(primero.consulta_alineada, reca[50..250].to_vec());
    assert_eq!(primero.linea_media, reca[50..250].to_vec());
}

#[test]
fn blastp_con_insercion_abre_un_hueco() {
    let reca = &ecoli::buscar_gen("recA").unwrap().proteina;
    let mut consulta = reca[50..150].to_vec();
    consulta.extend_from_slice(b"GPGPG");
    consulta.extend_from_slice(&reca[150..250]);
    let aciertos = blastp(&consulta, proteinas(), 1e-3, 5).unwrap();
    let primero = &aciertos[0];
    assert_eq!(primero.sujeto, indice_en_base("recA"));
    assert_eq!(primero.identidades, 200);
    assert_eq!(primero.aperturas, 1);
    assert_eq!(primero.huecos, 5);
    assert_eq!(primero.longitud, 205);
    assert_eq!((primero.inicio_sujeto, primero.fin_sujeto), (50, 250));
    assert_eq!(
        primero
            .sujeto_alineado
            .iter()
            .filter(|&&letra| letra == b'-')
            .count(),
        5
    );
}

#[test]
fn linea_media_marca_identicos_y_positivos() {
    let sujeto = al_azar(AMINOACIDOS, 120, 7);
    let mut consulta = sujeto[10..110].to_vec();
    let lugar = (20..60).find(|&i| consulta[i] != b'C').unwrap();
    let original_positivo = consulta[lugar];
    let positivo = *AMINOACIDOS
        .iter()
        .find(|&&otra| otra != original_positivo && blosum62(original_positivo, otra) > 0)
        .unwrap();
    consulta[lugar] = positivo;
    let original_negativo = consulta[60];
    let negativo = *AMINOACIDOS
        .iter()
        .find(|&&otra| blosum62(original_negativo, otra) < 0)
        .unwrap();
    consulta[60] = negativo;
    let base = indexar_proteinas(vec![Registro {
        nombre: "inventada".to_string(),
        secuencia: sujeto.clone(),
    }])
    .unwrap();
    let aciertos = blastp(&consulta, &base, 1e-3, 5).unwrap();
    assert_eq!(aciertos.len(), 1);
    let unico = &aciertos[0];
    assert_eq!((unico.inicio_sujeto, unico.fin_sujeto), (10, 110));
    assert_eq!(unico.longitud, 100);
    assert_eq!(unico.identidades, 98);
    assert_eq!(unico.positivos, 99);
    assert_eq!(unico.desapareamientos, 2);
    assert_eq!(unico.linea_media[lugar], b'+');
    assert_eq!(unico.linea_media[60], b' ');
    assert_eq!(unico.linea_media[0], consulta[0]);
}

#[test]
fn blastp_sensible_con_treinta_por_ciento_de_mutaciones() {
    for (numero, nombre) in ["recA", "parC", "lacZ", "dnaA", "rpoB"].iter().enumerate() {
        let original = &ecoli::buscar_gen(nombre).unwrap().proteina;
        let mutada = mutar(original, AMINOACIDOS, 0.30, numero as u64 + 11);
        let aciertos = blastp(&mutada, proteinas(), 1e-3, 10).unwrap();
        let propio = aciertos
            .iter()
            .find(|acierto| acierto.sujeto == indice_en_base(nombre))
            .unwrap_or_else(|| panic!("{nombre} mutado no se encontró"));
        assert!(propio.evalor < 1e-10, "{nombre}: E = {}", propio.evalor);
        assert_eq!(aciertos[0].sujeto, propio.sujeto, "{nombre}");
    }
}

#[test]
fn blastp_al_azar_no_encuentra_nada() {
    for semilla in 1..=3 {
        let consulta = al_azar(AMINOACIDOS, 300, semilla);
        let aciertos = blastp(&consulta, proteinas(), 1e-3, 10).unwrap();
        assert!(aciertos.is_empty(), "semilla {semilla}: {}", aciertos.len());
    }
}

#[test]
fn blastp_rechaza_entradas_malas() {
    assert!(blastp(b"", proteinas(), 1e-3, 10).is_err());
    assert!(blastp(b"MK", proteinas(), 1e-3, 10).is_err());
    let larguisima = vec![b'A'; LONGITUD_MAXIMA_PROTEINA + 1];
    assert!(blastp(&larguisima, proteinas(), 1e-3, 10).is_err());
    assert!(indexar_proteinas(Vec::new()).is_err());
    let sin_palabras = blastp(b"XXXXX", proteinas(), 1e-3, 10).unwrap();
    assert!(sin_palabras.is_empty());
    let con_minusculas = blastp(b"mkvlhewtra", proteinas(), 10.0, 3);
    assert!(con_minusculas.is_ok());
    assert!(
        blastp(b"WWWWWWWWWW", proteinas(), 1e-3, 0)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn blastn_region_exacta_en_ambas_hebras() {
    let region = ecoli::genoma()[1000..2000].to_vec();
    let directa = blastn(&region, genoma(), 1e-3, 5).unwrap();
    assert_eq!(directa[0].marco, 1);
    assert_eq!(
        (directa[0].inicio_sujeto, directa[0].fin_sujeto),
        (1000, 2000)
    );
    assert_eq!(
        (directa[0].inicio_consulta, directa[0].fin_consulta),
        (0, 1000)
    );
    assert_eq!(directa[0].identidades, 1000);
    assert_eq!(directa[0].linea_media, vec![b'|'; 1000]);
    let inversa = complemento_inverso(&region);
    let contraria = blastn(&inversa, genoma(), 1e-3, 5).unwrap();
    let primero = &contraria[0];
    assert_eq!(primero.marco, -1);
    assert_eq!((primero.inicio_sujeto, primero.fin_sujeto), (1000, 2000));
    assert_eq!((primero.inicio_consulta, primero.fin_consulta), (0, 1000));
    assert_eq!(primero.consulta_alineada, inversa);
    assert_eq!(primero.sujeto_alineado, inversa);
    assert_eq!(primero.puntaje, directa[0].puntaje);
}

#[test]
fn blastn_fragmento_mutado_de_tufa_encuentra_tufa_y_tufb() {
    let tufa = ecoli::buscar_gen("tufA").unwrap();
    let tufb = ecoli::buscar_gen("tufB").unwrap();
    let codificante = ecoli::secuencia_gen(tufa);
    let fragmento = mutar(&codificante[200..700], b"ACGT", 0.05, 3);
    let aciertos = blastn(&fragmento, genoma(), 1e-3, 10).unwrap();
    assert!(aciertos.len() >= 2);
    let dentro = |acierto: &Acierto, elegido: &ecoli::Gen| {
        acierto.inicio_sujeto + 1 >= elegido.inicio && acierto.fin_sujeto <= elegido.fin
    };
    let en_tufa = aciertos.iter().find(|a| dentro(a, tufa)).unwrap();
    let en_tufb = aciertos.iter().find(|a| dentro(a, tufb)).unwrap();
    assert_eq!(en_tufa.marco, -1);
    assert_eq!(en_tufb.marco, 1);
    for acierto in [en_tufa, en_tufb] {
        assert!(acierto.evalor < 1e-100);
        let identidad = acierto.identidades as f64 / acierto.longitud as f64;
        assert!(identidad > 0.9 && identidad < 0.99, "{identidad}");
        assert!(acierto.longitud >= 480);
    }
}

#[test]
fn blastn_base_propia_con_dos_registros() {
    let uno = al_azar(b"ACGT", 3000, 1);
    let dos = al_azar(b"ACGT", 3000, 2);
    let base = indexar_adn(vec![
        Registro {
            nombre: "uno".to_string(),
            secuencia: uno.clone(),
        },
        Registro {
            nombre: "dos".to_string(),
            secuencia: dos.clone(),
        },
    ])
    .unwrap();
    assert_eq!(base.cantidad(), 2);
    assert_eq!(base.total_letras, 6000);
    let aciertos = blastn(&dos[1000..1400], &base, 1e-3, 5).unwrap();
    assert_eq!(aciertos[0].sujeto, 1);
    assert_eq!(
        (aciertos[0].inicio_sujeto, aciertos[0].fin_sujeto),
        (1000, 1400)
    );
    let mut puente = uno[2800..].to_vec();
    puente.extend_from_slice(&dos[..200]);
    let mut partes = blastn(&puente, &base, 1e-3, 5).unwrap();
    partes.sort_by_key(|acierto| acierto.sujeto);
    assert_eq!(partes.len(), 2);
    assert_eq!(
        (
            partes[0].sujeto,
            partes[0].inicio_sujeto,
            partes[0].fin_sujeto
        ),
        (0, 2800, 3000)
    );
    assert_eq!(
        (partes[0].inicio_consulta, partes[0].fin_consulta),
        (0, 200)
    );
    assert_eq!(
        (
            partes[1].sujeto,
            partes[1].inicio_sujeto,
            partes[1].fin_sujeto
        ),
        (1, 0, 200)
    );
    assert_eq!(
        (partes[1].inicio_consulta, partes[1].fin_consulta),
        (200, 400)
    );
}

#[test]
fn blastn_al_azar_no_encuentra_nada() {
    for semilla in 1..=3 {
        let consulta = al_azar(b"ACGT", 500, semilla + 100);
        let aciertos = blastn(&consulta, genoma(), 1e-3, 10).unwrap();
        assert!(aciertos.is_empty(), "semilla {semilla}");
    }
}

#[test]
fn blastn_rechaza_entradas_malas() {
    assert!(blastn(b"", genoma(), 1e-3, 10).is_err());
    assert!(blastn(b"ACGTACGTAC", genoma(), 1e-3, 10).is_err());
    let larguisima = vec![b'A'; LONGITUD_MAXIMA_ADN + 1];
    assert!(blastn(&larguisima, genoma(), 1e-3, 10).is_err());
    assert!(indexar_adn(Vec::new()).is_err());
    let repetida: Vec<u8> = b"GCGCCGCATCCGGCA"
        .iter()
        .copied()
        .cycle()
        .take(LONGITUD_MAXIMA_ADN)
        .collect();
    let error = blastn(&repetida, genoma(), 1e-3, 10).unwrap_err();
    assert!(error.contains("coincidencias"), "{error}");
    let solo_n = blastn(b"NNNNNNNNNNNNNNNNNNNN", genoma(), 1e-3, 10).unwrap();
    assert!(solo_n.is_empty());
    let con_n = ecoli::genoma()[5000..5400]
        .iter()
        .enumerate()
        .map(|(i, &base)| if i % 50 == 25 { b'N' } else { base })
        .collect::<Vec<u8>>();
    let aciertos = blastn(&con_n, genoma(), 1e-3, 10).unwrap();
    assert_eq!(aciertos[0].inicio_sujeto, 5000);
}

#[test]
fn blastx_marcos_positivos_y_negativos() {
    let reca = ecoli::buscar_gen("recA").unwrap();
    let codificante = ecoli::secuencia_gen(reca);
    let sujeto = indice_en_base("recA");
    let marco_uno = blastx(&codificante[150..750], proteinas(), 1e-3, 5).unwrap();
    assert_eq!(marco_uno[0].sujeto, sujeto);
    assert_eq!(marco_uno[0].marco, 1);
    assert_eq!(
        (marco_uno[0].inicio_consulta, marco_uno[0].fin_consulta),
        (0, 600)
    );
    assert_eq!(
        (marco_uno[0].inicio_sujeto, marco_uno[0].fin_sujeto),
        (50, 250)
    );
    assert_eq!(marco_uno[0].identidades, 200);
    let marco_tres = blastx(&codificante[151..751], proteinas(), 1e-3, 5).unwrap();
    assert_eq!(marco_tres[0].marco, 3);
    assert_eq!(
        (marco_tres[0].inicio_consulta, marco_tres[0].fin_consulta),
        (2, 599)
    );
    assert_eq!(
        (marco_tres[0].inicio_sujeto, marco_tres[0].fin_sujeto),
        (51, 250)
    );
    let inversa = complemento_inverso(&codificante[150..750]);
    let negativo = blastx(&inversa, proteinas(), 1e-3, 5).unwrap();
    assert_eq!(negativo[0].marco, -1);
    assert_eq!(
        (negativo[0].inicio_consulta, negativo[0].fin_consulta),
        (0, 600)
    );
    assert_eq!(
        (negativo[0].inicio_sujeto, negativo[0].fin_sujeto),
        (50, 250)
    );
    let corrida = complemento_inverso(&codificante[149..750]);
    let negativo_dos = blastx(&corrida, proteinas(), 1e-3, 5).unwrap();
    assert_eq!(negativo_dos[0].marco, -2);
    assert_eq!(
        (
            negativo_dos[0].inicio_consulta,
            negativo_dos[0].fin_consulta
        ),
        (0, 600)
    );
}

#[test]
fn blastx_cambio_de_marco_da_dos_aciertos() {
    let codificante = ecoli::secuencia_gen(ecoli::buscar_gen("recA").unwrap());
    let mut consulta = codificante[150..450].to_vec();
    consulta.extend_from_slice(&codificante[451..750]);
    let aciertos = blastx(&consulta, proteinas(), 1e-3, 10).unwrap();
    let sujeto = indice_en_base("recA");
    let mut marcos: Vec<i8> = aciertos
        .iter()
        .filter(|acierto| acierto.sujeto == sujeto && acierto.evalor < 1e-20)
        .map(|acierto| acierto.marco)
        .collect();
    marcos.sort();
    marcos.dedup();
    assert_eq!(marcos, vec![1, 3]);
}

#[test]
fn blastx_rechaza_entradas_malas() {
    assert!(blastx(b"", proteinas(), 1e-3, 10).is_err());
    assert!(blastx(b"ACGTACGT", proteinas(), 1e-3, 10).is_err());
    let larguisima = vec![b'A'; LONGITUD_MAXIMA_BLASTX + 1];
    assert!(blastx(&larguisima, proteinas(), 1e-3, 10).is_err());
    assert!(
        blastx(b"NNNNNNNNNNNN", proteinas(), 1e-3, 10)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn consulta_de_una_sola_letra_repetida_termina_pronto() {
    let inicio = std::time::Instant::now();
    for letra in *b"AL" {
        let consulta = vec![letra; LONGITUD_MAXIMA_PROTEINA];
        let aciertos = blastp(&consulta, proteinas(), 1e-3, 20).unwrap();
        assert!(aciertos.len() <= 20);
    }
    assert!(inicio.elapsed().as_secs() < 8);
}

#[test]
fn proteina_de_la_base_demasiado_larga_se_rechaza() {
    let larga = Registro {
        nombre: "genoma".to_string(),
        secuencia: vec![b'A'; LONGITUD_MAXIMA_SUJETO_PROTEINA + 1],
    };
    let error = indexar_proteinas(vec![larga]).err().unwrap();
    assert!(error.contains("genoma"));
    let justa = Registro {
        nombre: "justa".to_string(),
        secuencia: vec![b'A'; 1000],
    };
    assert!(indexar_proteinas(vec![justa]).is_ok());
}
