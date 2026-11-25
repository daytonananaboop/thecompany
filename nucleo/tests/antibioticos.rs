use nucleo::antibioticos::*;
use nucleo::ecoli;

fn ordenados(mut peptidos: Vec<Vec<usize>>) -> Vec<Vec<usize>> {
    peptidos.sort();
    peptidos
}

fn mismas_masas(peptido: &[usize], esperado: &[usize]) -> bool {
    let mut a = peptido.to_vec();
    let mut b = esperado.to_vec();
    a.sort();
    b.sort();
    a == b
}

#[test]
fn traduce_hasta_el_primer_codon_de_parada() {
    let resultado = traducir_arn(b"AUGGCCUAAGGG").unwrap();
    assert_eq!(resultado.peptido, b"MA");
    assert_eq!(resultado.parada, Some(6));
    let sin_parada = traducir_arn(b"AUGUUUCC").unwrap();
    assert_eq!(sin_parada.peptido, b"MF");
    assert_eq!(sin_parada.parada, None);
    let desde_adn = traducir_arn(b"ATGTGGTGA").unwrap();
    assert_eq!(desde_adn.peptido, b"MW");
    assert_eq!(desde_adn.parada, Some(6));
    let vacio = traducir_arn(b"").unwrap();
    assert!(vacio.peptido.is_empty());
    assert_eq!(vacio.parada, None);
    let parada_inicial = traducir_arn(b"UAGAUG").unwrap();
    assert!(parada_inicial.peptido.is_empty());
    assert_eq!(parada_inicial.parada, Some(0));
}

#[test]
fn traducir_rechaza_codones_con_letras_raras() {
    assert!(traducir_arn(b"AUGNNN").is_err());
    assert!(traducir_arn(b"AUGUAANNN").is_ok());
}

#[test]
fn traduce_un_gen_real_de_e_coli() {
    let thr = ecoli::buscar_gen("thrL").unwrap();
    let adn = ecoli::secuencia_gen(thr);
    let resultado = traducir_arn(&adn).unwrap();
    assert_eq!(resultado.peptido, thr.proteina);
    assert_eq!(resultado.parada, Some(adn.len() - 3));
}

#[test]
fn encuentra_codificaciones_en_las_dos_hebras() {
    let encontradas = codificaciones_del_peptido(b"CAAGTTCGAAATTTG", b"KF").unwrap();
    let resumen: Vec<(usize, u8)> = encontradas.iter().map(|c| (c.inicio, c.hebra)).collect();
    assert_eq!(resumen, vec![(1, b'+'), (8, b'+'), (8, b'-')]);
    assert_eq!(encontradas[0].fragmento, b"AAGTTC");
    assert_eq!(encontradas[1].fragmento, b"AAATTT");
    assert_eq!(encontradas[2].fragmento, b"AAATTT");
}

#[test]
fn una_codificacion_en_la_hebra_menos_se_lee_al_reves() {
    let encontradas = codificaciones_del_peptido(b"TTGAACCATTT", b"MV").unwrap();
    assert_eq!(encontradas.len(), 1);
    assert_eq!(encontradas[0].inicio, 3);
    assert_eq!(encontradas[0].hebra, b'-');
    assert_eq!(encontradas[0].fragmento, b"AACCAT");
}

#[test]
fn codificaciones_que_se_solapan_con_codones_distintos() {
    let encontradas = codificaciones_del_peptido(b"CACCATCAC", b"HH").unwrap();
    let inicios: Vec<usize> = encontradas.iter().map(|c| c.inicio).collect();
    assert_eq!(inicios, vec![0, 3]);
    assert_eq!(encontradas[0].fragmento, b"CACCAT");
    assert_eq!(encontradas[1].fragmento, b"CATCAC");
    assert!(encontradas.iter().all(|c| c.hebra == b'+'));
}

#[test]
fn codificaciones_casos_limite() {
    assert!(codificaciones_del_peptido(b"ACGT", b"").is_err());
    assert!(codificaciones_del_peptido(b"ACGTACGT", b"MB").is_err());
    assert!(
        codificaciones_del_peptido(b"ATG", b"MM")
            .unwrap()
            .is_empty()
    );
    assert!(codificaciones_del_peptido(b"", b"M").unwrap().is_empty());
    let con_n = codificaciones_del_peptido(b"ATGNNNATG", b"M").unwrap();
    assert_eq!(con_n.len(), 2);
}

#[test]
fn codificaciones_en_genes_reales() {
    let his = ecoli::buscar_gen("hisL").unwrap();
    let adn = ecoli::secuencia_gen(his);
    let encontradas = codificaciones_del_peptido(&adn, &his.proteina).unwrap();
    assert_eq!(encontradas.len(), 1);
    assert_eq!(encontradas[0].inicio, 0);
    assert_eq!(encontradas[0].hebra, b'+');
    let trp = ecoli::buscar_gen("trpL").unwrap();
    assert_eq!(trp.hebra, b'-');
    let region = ecoli::region(trp.inicio, trp.fin).unwrap();
    let encontradas = codificaciones_del_peptido(region, &trp.proteina).unwrap();
    assert_eq!(encontradas.len(), 1);
    assert_eq!(encontradas[0].inicio, 3);
    assert_eq!(encontradas[0].hebra, b'-');
}

#[test]
fn cuenta_las_secuencias_que_codifican() {
    assert_eq!(codones_de(b'M'), 1);
    assert_eq!(codones_de(b'L'), 6);
    assert_eq!(codones_de(b'H'), 2);
    assert_eq!(codones_de(b'*'), 3);
    assert_eq!(secuencias_que_codifican(b"MW"), 1.0);
    assert_eq!(secuencias_que_codifican(b"HHH"), 8.0);
    assert_eq!(secuencias_que_codifican(b"LSR"), 216.0);
}

#[test]
fn tabla_de_masas_enteras() {
    assert_eq!(masa_entera(b'G'), Some(57));
    assert_eq!(masa_entera(b'W'), Some(186));
    assert_eq!(masa_entera(b'I'), masa_entera(b'L'));
    assert_eq!(masa_entera(b'K'), masa_entera(b'Q'));
    assert_eq!(masa_entera(b'B'), None);
    assert_eq!(MASAS_DISTINTAS.len(), 18);
    assert!(MASAS_DISTINTAS.windows(2).all(|par| par[0] < par[1]));
    assert_eq!(masas_del_peptido(b"GAL").unwrap(), vec![57, 71, 113]);
    assert!(masas_del_peptido(b"GAX").is_err());
}

#[test]
fn formulas_coinciden_con_masas_enteras() {
    for &(letra, masa) in MASAS_ENTERAS.iter() {
        let formula = formula_del_peptido(&[letra], true).unwrap();
        let entera = 12 * formula.carbono
            + formula.hidrogeno
            + 14 * formula.nitrogeno
            + 16 * formula.oxigeno
            + 32 * formula.azufre;
        assert_eq!(entera as usize, masa, "{}", letra as char);
    }
}

#[test]
fn masa_monoisotopica_de_glicina_y_diferencia_k_q() {
    let formula = formula_del_peptido(b"G", false).unwrap();
    assert_eq!(
        (
            formula.carbono,
            formula.hidrogeno,
            formula.nitrogeno,
            formula.oxigeno,
            formula.azufre
        ),
        (2, 5, 1, 2, 0)
    );
    let glicina = masa_monoisotopica(b"G", false).unwrap();
    assert!((glicina - 75.03203).abs() < 1e-4);
    let residuo = masa_monoisotopica(b"G", true).unwrap();
    assert!((residuo - 57.02146).abs() < 1e-4);
    let diferencia =
        masa_monoisotopica(b"K", true).unwrap() - masa_monoisotopica(b"Q", true).unwrap();
    assert!(diferencia > 0.036 && diferencia < 0.037);
    assert_eq!(
        masa_monoisotopica(b"I", true).unwrap(),
        masa_monoisotopica(b"L", true).unwrap()
    );
    assert!(masa_monoisotopica(b"GZ", false).is_err());
}

#[test]
fn espectro_ciclico_y_lineal_pequenos() {
    assert_eq!(
        espectro_ciclico(&[57, 71, 113]),
        vec![0, 57, 71, 113, 128, 170, 184, 241]
    );
    assert_eq!(
        espectro_lineal(&[57, 71, 113]),
        vec![0, 57, 71, 113, 128, 184, 241]
    );
    assert_eq!(
        espectro_ciclico(&[103, 71, 131, 71]),
        vec![
            0, 71, 71, 103, 131, 174, 174, 202, 202, 245, 273, 305, 305, 376
        ]
    );
    assert_eq!(espectro_ciclico(&[]), vec![0]);
    assert_eq!(espectro_lineal(&[]), vec![0]);
    assert_eq!(espectro_ciclico(&[57]), vec![0, 57]);
    assert_eq!(espectro_lineal(&[57]), vec![0, 57]);
}

#[test]
fn longitud_de_los_espectros() {
    for n in 1..12 {
        let peptido: Vec<usize> = (0..n).map(|i| MASAS_DISTINTAS[i % 18]).collect();
        assert_eq!(espectro_ciclico(&peptido).len(), n * (n - 1) + 2);
        assert_eq!(espectro_lineal(&peptido).len(), n * (n + 1) / 2 + 1);
    }
}

#[test]
fn espectro_ciclico_no_cambia_al_rotar_ni_invertir() {
    let peptido = vec![97, 71, 113, 131, 129, 156, 71];
    let original = espectro_ciclico(&peptido);
    for giro in 0..peptido.len() {
        let mut rotado = peptido.clone();
        rotado.rotate_left(giro);
        assert_eq!(espectro_ciclico(&rotado), original);
        rotado.reverse();
        assert_eq!(espectro_ciclico(&rotado), original);
    }
}

#[test]
fn cuenta_peptidos_por_masa() {
    let cuentas = peptidos_por_masa(171).unwrap();
    assert_eq!(cuentas[0], 1);
    assert_eq!(cuentas[56], 0);
    assert_eq!(cuentas[57], 1);
    assert_eq!(cuentas[114], 2);
    assert_eq!(cuentas[128], 3);
    assert_eq!(cuentas[171], 3);
    assert_eq!(cuentas.len(), 172);
    assert!(peptidos_por_masa(3001).is_err());
    let grande = peptidos_por_masa(3000).unwrap();
    assert!(grande[3000] > grande[2000]);
    assert!(grande[2000] > grande[1000]);
}

#[test]
fn secuencia_un_tripeptido_ciclico() {
    let encontrados = secuenciar_ciclopeptido(&[0, 57, 71, 113, 128, 170, 184, 241]).unwrap();
    assert_eq!(
        ordenados(encontrados),
        vec![
            vec![57, 71, 113],
            vec![57, 113, 71],
            vec![71, 57, 113],
            vec![71, 113, 57],
            vec![113, 57, 71],
            vec![113, 71, 57],
        ]
    );
}

#[test]
fn secuencia_un_peptido_simetrico_y_uno_repetido() {
    let encontrados = secuenciar_ciclopeptido(&espectro_ciclico(&[103, 71, 131, 71])).unwrap();
    assert_eq!(
        ordenados(encontrados),
        vec![
            vec![71, 103, 71, 131],
            vec![71, 131, 71, 103],
            vec![103, 71, 131, 71],
            vec![131, 71, 103, 71],
        ]
    );
    let repetido = secuenciar_ciclopeptido(&[0, 57, 57, 114]).unwrap();
    assert_eq!(repetido, vec![vec![57, 57]]);
}

#[test]
fn secuenciar_acepta_el_espectro_desordenado() {
    let encontrados = secuenciar_ciclopeptido(&[241, 113, 0, 184, 57, 170, 128, 71]).unwrap();
    assert_eq!(encontrados.len(), 6);
}

#[test]
fn espectro_para_ramificar_revisa_los_limites() {
    let preparado = espectro_para_ramificar(&[241, 0, 57, 57]).unwrap();
    assert_eq!(preparado.masa_madre, 241);
    assert_eq!(preparado.masas, vec![0, 57, 57, 241]);
    assert_eq!(preparado.veces(57), 2);
    assert_eq!(preparado.veces(71), 0);
    assert_eq!(preparado.veces(5000), 0);
    let muchas: Vec<usize> = (0..=MAXIMO_MASAS_EXACTAS).collect();
    assert!(espectro_para_ramificar(&muchas).is_err());
    let justas: Vec<usize> = (0..MAXIMO_MASAS_EXACTAS).collect();
    assert!(espectro_para_ramificar(&justas).is_ok());
    assert!(espectro_para_ramificar(&[]).is_err());
}

#[test]
fn conservar_consistentes_separa_soluciones_y_poda() {
    let mut experimental = Espectro::nuevo(&[0, 57, 71, 113, 128, 170, 184, 241]).unwrap();
    assert_eq!(aminoacidos_presentes(&experimental), vec![57, 71, 113, 128]);
    let candidatos = vec![
        vec![57, 71],
        vec![71, 99],
        vec![57, 71, 113],
        vec![113, 57, 71],
        vec![57, 128],
        vec![71, 71],
        vec![113, 128],
    ];
    let mut encontrados = Vec::new();
    let conservados = conservar_consistentes(candidatos, &mut experimental, &mut encontrados);
    assert_eq!(conservados, vec![vec![57, 71]]);
    assert_eq!(encontrados, vec![vec![57, 71, 113], vec![113, 57, 71]]);
}

#[test]
fn secuenciar_casos_limite() {
    assert!(secuenciar_ciclopeptido(&[]).is_err());
    assert!(secuenciar_ciclopeptido(&[0, MASA_MAXIMA + 1]).is_err());
    let muchas: Vec<usize> = (0..=MAXIMO_MASAS_EXACTAS).collect();
    assert!(secuenciar_ciclopeptido(&muchas).is_err());
    assert!(secuenciar_ciclopeptido(&[57, 71, 128]).unwrap().is_empty());
    assert!(secuenciar_ciclopeptido(&[0]).unwrap().is_empty());
    let con_ruido = secuenciar_ciclopeptido(&[0, 57, 71, 99, 113, 128, 184, 241]).unwrap();
    assert!(con_ruido.is_empty());
}

#[test]
fn expandir_alarga_y_respeta_el_limite() {
    let expandidos = expandir(&[vec![57], vec![71]], &[57, 113]).unwrap();
    assert_eq!(
        expandidos,
        vec![vec![57, 57], vec![57, 113], vec![71, 57], vec![71, 113]]
    );
    let muchos = vec![Vec::new(); MAXIMO_CANDIDATOS / 18 + 1];
    assert!(expandir(&muchos, &MASAS_DISTINTAS).is_err());
}

#[test]
fn puntajes_cuentan_masas_compartidas_con_repeticion() {
    let exacto = [0, 57, 71, 113, 128, 170, 184, 241];
    assert_eq!(puntaje_ciclico(&[57, 71, 113], &exacto).unwrap(), 8);
    assert_eq!(puntaje_lineal(&[57, 71, 113], &exacto).unwrap(), 7);
    let ruidoso = [0, 57, 57, 71, 99, 128, 241];
    assert_eq!(puntaje_ciclico(&[57, 71, 113], &ruidoso).unwrap(), 5);
    assert_eq!(puntaje_lineal(&[57, 71, 113], &ruidoso).unwrap(), 5);
    assert_eq!(
        puntaje_ciclico(&[103, 71, 131, 71], &[0, 71, 103]).unwrap(),
        3
    );
    assert_eq!(
        puntaje_ciclico(&[103, 71, 131, 71], &[0, 71, 71, 71, 103]).unwrap(),
        4
    );
    assert_eq!(puntaje_ciclico(&[57], &[500, 900]).unwrap(), 0);
    assert!(puntaje_ciclico(&[57], &[]).is_err());
}

#[test]
fn recortar_conserva_los_empates() {
    let espectro = [0, 57, 71, 128];
    let tabla = vec![
        vec![57, 71],
        vec![71, 57],
        vec![57, 57],
        vec![113],
        vec![71],
    ];
    assert_eq!(
        recortar(tabla.clone(), &espectro, 2).unwrap(),
        vec![vec![57, 71], vec![71, 57]]
    );
    assert_eq!(
        recortar(tabla.clone(), &espectro, 3).unwrap(),
        vec![vec![57, 71], vec![71, 57], vec![57, 57], vec![71]]
    );
    assert_eq!(recortar(tabla.clone(), &espectro, 10).unwrap(), tabla);
    assert!(recortar(tabla.clone(), &espectro, 0).unwrap().is_empty());
}

#[test]
fn tabla_de_lideres_tolera_ruido() {
    let ruidoso = [0, 57, 71, 99, 113, 128, 184, 241];
    let lideres = secuenciar_con_lideres(&ruidoso, 10).unwrap();
    assert_eq!(lideres.puntaje, 7);
    assert_eq!(lideres.peptidos.len(), 6);
    assert!(
        lideres
            .peptidos
            .iter()
            .all(|peptido| mismas_masas(peptido, &[57, 71, 113]))
    );
    assert_eq!(lideres.alfabeto, MASAS_DISTINTAS.to_vec());
}

#[test]
fn tabla_de_lideres_casos_limite() {
    assert!(secuenciar_con_lideres(&[], 5).is_err());
    let sin_espacio = secuenciar_con_lideres(&[0, 57, 71, 128], 0).unwrap();
    assert_eq!(sin_espacio.peptidos, vec![vec![128]]);
    assert_eq!(sin_espacio.puntaje, 2);
    let inalcanzable = secuenciar_con_lideres(&[1, 2, 3], 5).unwrap();
    assert!(inalcanzable.peptidos.is_empty());
}

#[test]
fn convolucion_de_un_espectro_pequeno() {
    assert_eq!(
        convolucion_espectral(&[0, 57, 71, 128]).unwrap(),
        vec![(57, 2), (71, 2), (14, 1), (128, 1)]
    );
    assert_eq!(convolucion_espectral(&[57, 0, 57]).unwrap(), vec![(57, 2)]);
    assert!(convolucion_espectral(&[5]).unwrap().is_empty());
    let espectro = espectro_ciclico(&[57, 69, 113]);
    let convolucion = convolucion_espectral(&espectro).unwrap();
    assert_eq!(&convolucion[..3], &[(57, 4), (69, 4), (113, 4)]);
    let total: usize = convolucion.iter().map(|&(_, veces)| veces).sum();
    assert_eq!(total, 28);
}

#[test]
fn alfabeto_por_convolucion_con_empates() {
    let espectro = espectro_ciclico(&[57, 69, 113]);
    let mut tres = alfabeto_por_convolucion(&espectro, 3).unwrap();
    tres.sort();
    assert_eq!(tres, vec![57, 69, 113]);
    let mut cuatro = alfabeto_por_convolucion(&espectro, 4).unwrap();
    cuatro.sort();
    assert_eq!(cuatro, vec![57, 69, 113, 126, 170, 182]);
    let mut uno = alfabeto_por_convolucion(&[0, 57, 71, 128], 1).unwrap();
    uno.sort();
    assert_eq!(uno, vec![57, 71]);
    assert!(alfabeto_por_convolucion(&espectro, 0).unwrap().is_empty());
}

#[test]
fn convolucion_descubre_masas_no_estandar() {
    let espectro = espectro_ciclico(&[57, 69, 113]);
    let lideres = secuenciar_por_convolucion(&espectro, 3, 10).unwrap();
    assert_eq!(lideres.puntaje, 8);
    assert_eq!(lideres.peptidos.len(), 6);
    assert!(
        lideres
            .peptidos
            .iter()
            .all(|peptido| mismas_masas(peptido, &[57, 69, 113]))
    );
    let estandar = secuenciar_con_lideres(&espectro, 10).unwrap();
    assert!(estandar.puntaje < 8);
    assert!(secuenciar_por_convolucion(&[0, 5, 10], 3, 10).is_err());
}

#[test]
fn reconstruye_puntos_desde_distancias() {
    let distancias = [9, 2, 17, 4, 11, 5, 13, 7, 8, 6];
    assert_eq!(
        puntos_desde_distancias(&distancias).unwrap(),
        Some(vec![0, 6, 8, 13, 17])
    );
    let con_ceros_y_negativos = [
        0, -17, 2, 0, -13, 4, -11, 5, 0, -9, 6, -8, 7, 0, -7, 8, -6, 9, -5, 0, -4, 11, -2, 13, 17,
    ];
    assert_eq!(
        puntos_desde_distancias(&con_ceros_y_negativos).unwrap(),
        Some(vec![0, 6, 8, 13, 17])
    );
    assert_eq!(puntos_desde_distancias(&[7]).unwrap(), Some(vec![0, 7]));
}

#[test]
fn autopista_sin_solucion_y_errores() {
    assert_eq!(puntos_desde_distancias(&[1, 1, 1]).unwrap(), None);
    assert!(puntos_desde_distancias(&[]).is_err());
    assert!(puntos_desde_distancias(&[0, -5]).is_err());
    assert!(puntos_desde_distancias(&[1, 2]).is_err());
    let demasiados: Vec<i64> = (1..=(41 * 40 / 2)).collect();
    assert!(puntos_desde_distancias(&demasiados).is_err());
}

fn distancias_de_una_regla(marcas: i64) -> Vec<i64> {
    let puntos: Vec<i64> = (0..marcas).map(|i| i * 10).collect();
    let mut distancias = Vec::new();
    for i in 0..puntos.len() {
        for j in i + 1..puntos.len() {
            distancias.push(puntos[j] - puntos[i]);
        }
    }
    distancias
}

#[test]
fn autopista_con_una_regla_y_con_una_distancia_corrida() {
    let regla: Vec<i64> = (0..40).map(|i| i * 10).collect();
    assert_eq!(
        puntos_desde_distancias(&distancias_de_una_regla(40)).unwrap(),
        Some(regla)
    );
    let mut corrida = distancias_de_una_regla(12);
    corrida[0] += 1;
    assert_eq!(puntos_desde_distancias(&corrida).unwrap(), None);
}

#[test]
fn autopista_se_rinde_con_un_mensaje_si_prueba_demasiado() {
    let mut corrida = distancias_de_una_regla(40);
    corrida[0] += 1;
    let mensaje = puntos_desde_distancias(&corrida).unwrap_err();
    assert!(mensaje.contains("un millón"), "{mensaje}");
}

#[test]
fn autopista_con_distancias_de_puntos_conocidos() {
    let puntos: Vec<i64> = vec![0, 13, 21, 50, 51, 77, 90, 140, 141, 199];
    let mut distancias = Vec::new();
    for i in 0..puntos.len() {
        for j in i + 1..puntos.len() {
            distancias.push(puntos[j] - puntos[i]);
        }
    }
    distancias.reverse();
    let reconstruidos = puntos_desde_distancias(&distancias).unwrap().unwrap();
    assert_eq!(reconstruidos.len(), puntos.len());
    assert_eq!(reconstruidos[0], 0);
    let mut recalculadas = Vec::new();
    for i in 0..reconstruidos.len() {
        for j in i + 1..reconstruidos.len() {
            recalculadas.push(reconstruidos[j] - reconstruidos[i]);
        }
    }
    recalculadas.sort();
    distancias.sort();
    assert_eq!(recalculadas, distancias);
}

#[test]
fn tabla_de_lideres_se_rinde_si_los_empates_la_hacen_crecer() {
    let espectro = espectro_ciclico(&[71; 56]);
    let mensaje = secuenciar_con_lideres(&espectro, 1000).unwrap_err();
    assert!(mensaje.contains("dos millones"), "{mensaje}");
}
