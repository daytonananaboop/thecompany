use nucleo::azar::Azar;
use nucleo::motivos::*;

fn adn(textos: &[&str]) -> Vec<Vec<u8>> {
    textos
        .iter()
        .map(|texto| texto.as_bytes().to_vec())
        .collect()
}

fn vistas<'a>(motivos: &[&'a str]) -> Vec<&'a [u8]> {
    motivos.iter().map(|motivo| motivo.as_bytes()).collect()
}

fn cerca(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

#[test]
fn hamming_cuenta_diferencias() {
    assert_eq!(distancia_hamming(b"ACGT", b"ACCA"), 2);
    assert_eq!(distancia_hamming(b"", b""), 0);
    assert_eq!(distancia_hamming(b"AAA", b"AAA"), 0);
}

#[test]
fn coincidencia_mas_cercana_elige_la_primera_ventana_empatada() {
    assert_eq!(coincidencia_mas_cercana(b"AAA", b"TTACCTTAAC"), (6, 1));
    assert_eq!(coincidencia_mas_cercana(b"GAT", b"CGATT"), (1, 0));
    assert_eq!(distancia_patron_texto(b"GAT", b"AAGGT"), 1);
    assert_eq!(distancia_patron_texto(b"GAT", b"TTTTT"), 2);
}

#[test]
fn coincidencia_con_entradas_raras_no_rompe() {
    assert_eq!(coincidencia_mas_cercana(b"", b"ACGT"), (0, 0));
    assert_eq!(coincidencia_mas_cercana(b"ACGTA", b"AC"), (0, 5));
}

#[test]
fn distancia_entre_patron_y_varias_secuencias() {
    let textos = adn(&["CGATT", "AAGGT", "TTTTT"]);
    assert_eq!(distancia_patron_textos(b"GAT", &textos), Ok(3));
    assert!(distancia_patron_textos(b"", &textos).is_err());
    assert!(distancia_patron_textos(b"GATTACA", &textos).is_err());
    assert!(distancia_patron_textos(b"GAT", &[]).is_err());
}

#[test]
fn enumeracion_exacta_es_la_interseccion() {
    let textos = adn(&["ACGT", "CGTA"]);
    assert_eq!(enumerar_motivos(&textos, 2, 0).unwrap(), adn(&["CG", "GT"]));
}

#[test]
fn enumeracion_con_un_error() {
    assert_eq!(
        enumerar_motivos(&adn(&["AA", "CC"]), 2, 1).unwrap(),
        adn(&["AC", "CA"])
    );
    assert_eq!(
        enumerar_motivos(&adn(&["AAAT", "TTTT"]), 3, 1).unwrap(),
        adn(&["ATT", "TAT"])
    );
}

#[test]
fn enumeracion_sin_resultados_y_con_todo() {
    assert!(
        enumerar_motivos(&adn(&["AAAA", "CCCC"]), 3, 0)
            .unwrap()
            .is_empty()
    );
    assert_eq!(enumerar_motivos(&adn(&["A", "C"]), 1, 1).unwrap().len(), 4);
}

#[test]
fn enumeracion_ignora_ventanas_con_letras_raras() {
    assert_eq!(
        enumerar_motivos(&adn(&["ANCG", "CGCG"]), 2, 0).unwrap(),
        adn(&["CG"])
    );
}

#[test]
fn enumeracion_rechaza_parametros_imposibles() {
    let textos = adn(&["ACGTACGTACGTACGT", "ACGTACGTACGTACGT"]);
    assert!(enumerar_motivos(&textos, 0, 0).is_err());
    assert!(enumerar_motivos(&textos, 11, 1).is_err());
    assert!(enumerar_motivos(&textos, 8, 4).is_err());
    assert!(enumerar_motivos(&textos, 2, 3).is_err());
    assert!(enumerar_motivos(&textos, 17, 1).is_err());
    assert!(enumerar_motivos(&[], 3, 1).is_err());
}

#[test]
fn mediana_unica() {
    let mediana = cadena_mediana(&adn(&["TACGTT", "GGACGA", "ACGCCC"]), 3).unwrap();
    assert_eq!(mediana.distancia, 0);
    assert_eq!(mediana.patrones, adn(&["ACG"]));
    assert_eq!(mediana.empatados, 1);
}

#[test]
fn mediana_con_empates_en_orden_alfabetico() {
    let mediana = cadena_mediana(&adn(&["AA", "CC"]), 2).unwrap();
    assert_eq!(mediana.distancia, 2);
    assert_eq!(mediana.patrones, adn(&["AA", "AC", "CA", "CC"]));
    assert_eq!(mediana.empatados, 4);
}

#[test]
fn mediana_guarda_pocos_empates_pero_los_cuenta_todos() {
    let mediana = cadena_mediana(&adn(&["AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT"]), 6).unwrap();
    assert_eq!(mediana.distancia, 18);
    assert_eq!(mediana.empatados, 4096);
    assert!(mediana.patrones.len() <= 100);
    assert_eq!(mediana.patrones[0], b"AAAAAA".to_vec());
}

#[test]
fn mediana_rechaza_k_grande_o_imposible() {
    let textos = adn(&["ACGTACGTACGTACGT", "ACGTACGTACGTACGT"]);
    assert!(cadena_mediana(&textos, 11).is_err());
    assert!(cadena_mediana(&textos, 0).is_err());
    assert!(cadena_mediana(&adn(&["ACG"]), 4).is_err());
}

#[test]
fn mediana_coincide_con_la_distancia() {
    let textos = adn(&["CTGATAT", "GGCTGAAT", "TTCTGAAT"]);
    let mediana = cadena_mediana(&textos, 4).unwrap();
    for patron in &mediana.patrones {
        assert_eq!(
            distancia_patron_textos(patron, &textos).unwrap(),
            mediana.distancia
        );
    }
    assert_eq!(mediana.distancia, 0);
    assert_eq!(mediana.patrones, adn(&["CTGA"]));
}

#[test]
fn perfil_sin_y_con_pseudoconteos() {
    let motivos = vistas(&["ACG", "ACT", "CCT"]);
    let perfil = perfil_de_motivos(&motivos, 0.0);
    assert_eq!(perfil.len(), 4);
    assert!(cerca(perfil[0][0], 2.0 / 3.0));
    assert!(cerca(perfil[1][0], 1.0 / 3.0));
    assert!(cerca(perfil[1][1], 1.0));
    assert!(cerca(perfil[2][2], 1.0 / 3.0));
    assert!(cerca(perfil[3][2], 2.0 / 3.0));
    let suavizado = perfil_de_motivos(&motivos, 1.0);
    assert!(cerca(suavizado[0][0], 3.0 / 7.0));
    assert!(cerca(suavizado[1][0], 2.0 / 7.0));
    assert!(cerca(suavizado[2][0], 1.0 / 7.0));
    assert!(cerca(suavizado[3][0], 1.0 / 7.0));
    for columna in 0..3 {
        let suma: f64 = suavizado.iter().map(|fila| fila[columna]).sum();
        assert!(cerca(suma, 1.0));
    }
}

#[test]
fn perfil_vacio_no_rompe() {
    let perfil = perfil_de_motivos(&[], 1.0);
    assert_eq!(perfil.len(), 4);
    assert!(perfil.iter().all(|fila| fila.is_empty()));
    assert!(contar_bases(&[]).iter().all(|fila| fila.is_empty()));
    assert!(consenso(&[]).is_empty());
    assert_eq!(puntuacion(&[]), 0);
    assert!(entropia_por_columna(&[]).is_empty());
}

#[test]
fn consenso_y_puntuacion() {
    let motivos = vistas(&["ACG", "ACT", "CCT"]);
    assert_eq!(consenso(&motivos), b"ACT".to_vec());
    assert_eq!(puntuacion(&motivos), 2);
    let empatados = vistas(&["AC", "CA"]);
    assert_eq!(consenso(&empatados), b"AA".to_vec());
    assert_eq!(puntuacion(&empatados), 2);
    assert_eq!(puntuacion(&vistas(&["GATTACA", "GATTACA"])), 0);
}

#[test]
fn letras_raras_cuentan_como_diferencia() {
    let motivos = vistas(&["AN", "AN"]);
    assert_eq!(
        contar_bases(&motivos),
        vec![vec![2, 0], vec![0, 0], vec![0, 0], vec![0, 0]]
    );
    assert_eq!(puntuacion(&motivos), 2);
}

#[test]
fn entropia_por_columna_en_bits() {
    let entropia = entropia_por_columna(&vistas(&["ACG", "ACT", "CCT"]));
    let esperada = -(2.0 / 3.0 * (2.0f64 / 3.0).log2() + 1.0 / 3.0 * (1.0f64 / 3.0).log2());
    assert!(cerca(entropia[0], esperada));
    assert!(cerca(entropia[1], 0.0));
    assert!(cerca(entropia[2], esperada));
    let uniforme = entropia_por_columna(&vistas(&["A", "C", "G", "T"]));
    assert!(cerca(uniforme[0], 2.0));
}

#[test]
fn analisis_completo_de_una_matriz_de_motivos() {
    let analisis = analizar_motivos(&vistas(&["ACG", "ACT", "CCT"])).unwrap();
    assert_eq!(analisis.conteos[1], vec![1, 3, 0]);
    assert_eq!(analisis.consenso, b"ACT".to_vec());
    assert_eq!(analisis.puntos, 2);
    assert_eq!(analisis.perfil.len(), 4);
    assert_eq!(analisis.perfil[0].len(), 3);
    assert_eq!(analisis.entropia.len(), 3);
    assert!(analizar_motivos(&[]).is_err());
    assert!(analizar_motivos(&vistas(&["ACG", "AC"])).is_err());
    assert!(analizar_motivos(&vistas(&["", ""])).is_err());
}

fn perfil_de_ejemplo() -> Vec<Vec<f64>> {
    vec![
        vec![0.5, 0.1],
        vec![0.5, 0.2],
        vec![0.0, 0.3],
        vec![0.0, 0.4],
    ]
}

#[test]
fn probabilidad_de_un_kmer() {
    let perfil = perfil_de_ejemplo();
    assert!(cerca(probabilidad_kmer(b"AC", &perfil), 0.1));
    assert!(cerca(probabilidad_kmer(b"ct", &perfil), 0.2));
    assert!(cerca(probabilidad_kmer(b"GA", &perfil), 0.0));
    assert!(cerca(probabilidad_kmer(b"NA", &perfil), 0.0));
    assert!(cerca(probabilidad_kmer(b"ACG", &perfil), 0.0));
}

#[test]
fn mas_probable_desempata_con_la_primera() {
    assert_eq!(posicion_mas_probable(b"GACTAT", &perfil_de_ejemplo()), 2);
}

#[test]
fn mas_probable_con_todo_en_cero_elige_la_primera() {
    let perfil = vec![vec![0.0], vec![0.0], vec![1.0], vec![0.0]];
    assert_eq!(posicion_mas_probable(b"AAAA", &perfil), 0);
}

#[test]
fn mas_probable_con_entradas_raras_no_rompe() {
    assert_eq!(posicion_mas_probable(b"ACGT", &[]), 0);
    assert_eq!(posicion_mas_probable(b"A", &vec![vec![0.25; 2]; 4]), 0);
    assert_eq!(posicion_mas_probable(b"ACGT", &vec![vec![0.25; 2]; 3]), 0);
    let desparejo = vec![vec![0.25; 2], vec![0.25; 2], vec![0.25; 2], vec![0.25]];
    assert_eq!(posicion_mas_probable(b"ACGT", &desparejo), 0);
    assert!(cerca(probabilidad_kmer(b"AT", &desparejo), 0.0));
    assert!(cerca(probabilidad_kmer(b"ACG", &desparejo), 0.0));
}

#[test]
fn motivos_en_recorta_cada_secuencia() {
    let textos = adn(&["GGACT", "CCACT"]);
    assert_eq!(
        motivos_en(&textos, &[2, 1], 3),
        vec![&b"ACT"[..], &b"CAC"[..]]
    );
    assert_eq!(motivos_en(&textos, &[9, 0], 3)[0], &b""[..]);
}

#[test]
fn codiciosa_sin_pseudoconteos_a_mano() {
    let textos = adn(&["GGACT", "CCACT", "TTACA"]);
    let hallado = busqueda_codiciosa(&textos, 3, 0.0).unwrap();
    assert_eq!(hallado.posiciones, vec![2, 2, 0]);
    assert_eq!(hallado.puntos, 3);
    assert_eq!(hallado.puntos_por_intento, vec![4, 5, 3]);
}

#[test]
fn codiciosa_con_pseudoconteos_a_mano() {
    let textos = adn(&["GGACT", "CCACT", "TTACA"]);
    let hallado = busqueda_codiciosa(&textos, 3, 1.0).unwrap();
    assert_eq!(hallado.posiciones, vec![2, 2, 2]);
    assert_eq!(hallado.puntos, 1);
    assert_eq!(hallado.puntos_por_intento, vec![3, 2, 1]);
}

#[test]
fn codiciosa_rechaza_entradas_malas() {
    let textos = adn(&["GGACT", "CCACT"]);
    assert!(busqueda_codiciosa(&textos, 0, 1.0).is_err());
    assert!(busqueda_codiciosa(&textos, 6, 1.0).is_err());
    assert!(busqueda_codiciosa(&textos, 3, -1.0).is_err());
    assert!(busqueda_codiciosa(&textos, 3, f64::NAN).is_err());
    assert!(busqueda_codiciosa(&[], 3, 1.0).is_err());
}

#[test]
fn codiciosa_con_una_sola_secuencia_no_rompe() {
    let hallado = busqueda_codiciosa(&adn(&["ACGTT"]), 2, 1.0).unwrap();
    assert_eq!(hallado.puntos, 0);
    assert_eq!(hallado.puntos_por_intento.len(), 4);
}

fn con_motivo_implantado(
    motivo: &[u8],
    cantidad: usize,
    longitud: usize,
    semilla: u64,
) -> (Vec<Vec<u8>>, Vec<usize>) {
    let mut azar = Azar::nuevo(semilla);
    let mut textos = Vec::new();
    let mut lugares = Vec::new();
    for _ in 0..cantidad {
        let mut texto: Vec<u8> = (0..longitud).map(|_| BASES[azar.entero(4)]).collect();
        let lugar = azar.entero(longitud - motivo.len() + 1);
        texto[lugar..lugar + motivo.len()].copy_from_slice(motivo);
        textos.push(texto);
        lugares.push(lugar);
    }
    (textos, lugares)
}

#[test]
fn codiciosa_con_pseudoconteos_encuentra_un_motivo_implantado() {
    let (textos, lugares) = con_motivo_implantado(b"CTGTATATACAG", 6, 60, 11);
    let hallado = busqueda_codiciosa(&textos, 12, 1.0).unwrap();
    assert_eq!(hallado.puntos, 0);
    assert_eq!(hallado.posiciones, lugares);
}

#[test]
fn aleatoria_encuentra_un_motivo_implantado() {
    let (textos, lugares) = con_motivo_implantado(b"TTGACAGCTA", 6, 40, 5);
    let hallado = busqueda_aleatoria(&textos, 10, 300, 1).unwrap();
    assert_eq!(hallado.puntos, 0);
    assert_eq!(hallado.posiciones, lugares);
    assert_eq!(hallado.puntos_por_intento.len(), 300);
    assert_eq!(
        hallado.puntos,
        *hallado.puntos_por_intento.iter().min().unwrap()
    );
}

#[test]
fn aleatoria_es_reproducible_con_la_misma_semilla() {
    let (textos, _) = con_motivo_implantado(b"GATTACA", 5, 50, 2);
    let primera = busqueda_aleatoria(&textos, 7, 20, 42).unwrap();
    let segunda = busqueda_aleatoria(&textos, 7, 20, 42).unwrap();
    assert_eq!(primera, segunda);
}

#[test]
fn aleatoria_rechaza_entradas_malas() {
    let textos = adn(&["GGACT", "CCACT"]);
    assert!(busqueda_aleatoria(&textos, 3, 0, 1).is_err());
    assert!(busqueda_aleatoria(&textos, 9, 5, 1).is_err());
    assert!(busqueda_aleatoria(&[], 3, 5, 1).is_err());
}

#[test]
fn gibbs_encuentra_un_motivo_implantado() {
    let (textos, lugares) = con_motivo_implantado(b"TTGACAGCTA", 6, 40, 5);
    let hallado = muestreo_gibbs(&textos, 10, 300, 20, 3).unwrap();
    assert_eq!(hallado.puntos, 0);
    assert_eq!(hallado.posiciones, lugares);
    assert_eq!(hallado.puntos_por_intento.len(), 20);
}

#[test]
fn gibbs_es_reproducible_y_acepta_cero_iteraciones() {
    let (textos, _) = con_motivo_implantado(b"GATTACA", 5, 50, 2);
    assert_eq!(
        muestreo_gibbs(&textos, 7, 50, 3, 9).unwrap(),
        muestreo_gibbs(&textos, 7, 50, 3, 9).unwrap()
    );
    let quieto = muestreo_gibbs(&textos, 7, 0, 1, 9).unwrap();
    assert_eq!(quieto.posiciones.len(), 5);
}

#[test]
fn gibbs_rechaza_entradas_malas() {
    let textos = adn(&["GGACT", "CCACT"]);
    assert!(muestreo_gibbs(&textos, 3, 10, 0, 1).is_err());
    assert!(muestreo_gibbs(&textos, 0, 10, 1, 1).is_err());
    assert!(muestreo_gibbs(&[], 3, 10, 1, 1).is_err());
}

#[test]
fn gibbs_con_una_sola_secuencia_no_rompe() {
    let hallado = muestreo_gibbs(&adn(&["ACGTTGCA"]), 3, 20, 2, 4).unwrap();
    assert_eq!(hallado.puntos, 0);
}

#[test]
fn las_busquedas_nunca_empeoran_el_punto_de_partida() {
    let (textos, _) = con_motivo_implantado(b"ACGTTGCAAC", 5, 80, 8);
    for pseudoconteo in [0.0, 1.0] {
        let hallado = busqueda_codiciosa(&textos, 10, pseudoconteo).unwrap();
        let iniciales = motivos_en(&textos, &[0; 5], 10);
        assert!(hallado.puntos <= puntuacion(&iniciales));
        assert_eq!(
            hallado.puntos,
            puntuacion(&motivos_en(&textos, &hallado.posiciones, 10))
        );
    }
    let gibbs = muestreo_gibbs(&textos, 10, 200, 5, 1).unwrap();
    assert_eq!(
        gibbs.puntos,
        puntuacion(&motivos_en(&textos, &gibbs.posiciones, 10))
    );
}

#[test]
fn la_mediana_rechaza_un_texto_sin_kmers_validos() {
    let error = cadena_mediana(&adn(&["NNNN", "ACGT"]), 2).unwrap_err();
    assert!(error.contains("secuencia 1"));
    let correcta = cadena_mediana(&adn(&["ACGT", "ACGT"]), 2).unwrap();
    assert_eq!(correcta.distancia, 0);
}
