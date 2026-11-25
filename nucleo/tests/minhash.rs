use nucleo::azar::Azar;
use nucleo::minhash::{
    Menores, Ventana, bosquejar, comparar, contar_compartidos, contencion_exacta, distancia_mash,
    estimar_kmeros_distintos, identidad_de_contencion, jaccard_al_azar, jaccard_exacto,
    log10_valor_p, matriz_de_distancias, mezclar, todos_los_kmeros,
};
use nucleo::secuencia::complemento_inverso;
use std::collections::HashSet;

fn al_azar(azar: &mut Azar, longitud: usize) -> Vec<u8> {
    (0..longitud).map(|_| b"ACGT"[azar.entero(4)]).collect()
}

fn mutar(azar: &mut Azar, secuencia: &[u8], fraccion: f64) -> Vec<u8> {
    let mut copia = secuencia.to_vec();
    let cambios = (fraccion * secuencia.len() as f64).round() as usize;
    let mut elegidas = HashSet::new();
    while elegidas.len() < cambios {
        elegidas.insert(azar.entero(secuencia.len()));
    }
    for posicion in elegidas {
        let original = copia[posicion];
        loop {
            let nueva = b"ACGT"[azar.entero(4)];
            if nueva != original {
                copia[posicion] = nueva;
                break;
            }
        }
    }
    copia
}

fn codigo_ingenuo(kmero: &[u8]) -> u64 {
    let mut codigo = 0u64;
    for &base in kmero {
        let valor = match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            _ => 3,
        };
        codigo = codigo * 4 + valor;
    }
    codigo
}

fn kmeros_ingenuos(secuencia: &[u8], k: usize) -> HashSet<u64> {
    let mayusculas: Vec<u8> = secuencia
        .iter()
        .map(|base| match base.to_ascii_uppercase() {
            b'U' => b'T',
            otra => otra,
        })
        .collect();
    let mut conjunto = HashSet::new();
    if mayusculas.len() < k {
        return conjunto;
    }
    for kmero in mayusculas.windows(k) {
        if kmero.iter().all(|base| b"ACGT".contains(base)) {
            let inverso = complemento_inverso(kmero);
            let canonico = codigo_ingenuo(kmero).min(codigo_ingenuo(&inverso));
            conjunto.insert(mezclar(canonico));
        }
    }
    conjunto
}

fn bosquejo_ingenuo(conjunto: &HashSet<u64>, cupo: usize) -> Vec<u64> {
    let mut ordenados: Vec<u64> = conjunto.iter().copied().collect();
    ordenados.sort();
    ordenados.truncate(cupo);
    ordenados
}

fn compartidos_ingenuos(a: &HashSet<u64>, b: &HashSet<u64>, cupo: usize) -> (usize, usize) {
    let mut union: Vec<u64> = a.union(b).copied().collect();
    union.sort();
    union.truncate(cupo);
    let compartidos = union
        .iter()
        .filter(|hash| a.contains(hash) && b.contains(hash))
        .count();
    (compartidos, union.len())
}

fn cola_binomial(compartidos: usize, vistos: usize, azar: f64) -> f64 {
    let mut cola = 0.0;
    for exitos in compartidos..=vistos {
        let mut combinaciones = 1.0;
        for i in 1..=exitos {
            combinaciones *= (vistos - exitos + i) as f64 / i as f64;
        }
        cola +=
            combinaciones * azar.powi(exitos as i32) * (1.0 - azar).powi((vistos - exitos) as i32);
    }
    cola
}

fn una_pieza(secuencia: &[u8]) -> Vec<Vec<u8>> {
    vec![secuencia.to_vec()]
}

#[test]
fn el_hash_coincide_con_el_calculado_a_mano() {
    assert_eq!(mezclar(0), 0x9ca0_66f1_a4ab_2eea);
    assert_eq!(mezclar(6), 0xe8b3_4dc9_dc73_356c);
    assert_eq!(mezclar(44), 0xef02_c688_5383_b23d);
}

#[test]
fn la_ventana_da_el_kmero_canonico() {
    let mut ventana = Ventana::nueva(3);
    let salidas: Vec<Option<u64>> = b"ACGTA".iter().map(|&b| ventana.agregar(b)).collect();
    assert_eq!(
        salidas,
        vec![
            None,
            None,
            Some(mezclar(6)),
            Some(mezclar(6)),
            Some(mezclar(44))
        ]
    );
    assert_eq!(ventana.adelante, 44);
    assert_eq!(ventana.atras, 49);
}

#[test]
fn una_letra_rara_corta_la_ventana() {
    let mut ventana = Ventana::nueva(4);
    let hallados: Vec<u64> = b"ACGTNACGT"
        .iter()
        .filter_map(|&b| ventana.agregar(b))
        .collect();
    assert_eq!(hallados.len(), 2);
    assert_eq!(hallados[0], hallados[1]);
    let bosquejo = bosquejar(&una_pieza(b"ACGTNACGT"), 4, 10);
    assert_eq!(bosquejo.kmeros, 2);
    assert_eq!(bosquejo.hashes.len(), 1);
    assert_eq!(bosquejo.longitud, 9);
}

#[test]
fn minusculas_y_u_cuentan_como_bases() {
    let mayusculas = bosquejar(&una_pieza(b"ACGTTGCAAGT"), 5, 100);
    let minusculas = bosquejar(&una_pieza(b"acguugcaagu"), 5, 100);
    assert_eq!(mayusculas.hashes, minusculas.hashes);
}

#[test]
fn el_bosquejo_de_acgta_a_mano() {
    let bosquejo = bosquejar(&una_pieza(b"ACGTA"), 3, 10);
    assert_eq!(bosquejo.hashes, vec![mezclar(6), mezclar(44)]);
    assert_eq!(bosquejo.kmeros, 3);
    assert_eq!(bosquejo.longitud, 5);
    let recortado = bosquejar(&una_pieza(b"ACGTA"), 3, 1);
    assert_eq!(recortado.hashes, vec![mezclar(6)]);
}

#[test]
fn los_menores_guardan_los_mas_chicos_sin_repetir() {
    let mut menores = Menores::nuevo(3);
    for hash in [50, 10, 40, 10, 30, 60, 20, 30] {
        menores.ofrecer(hash);
    }
    let guardados: Vec<u64> = menores.conjunto.iter().copied().collect();
    assert_eq!(guardados, vec![10, 20, 30]);
    assert_eq!(menores.mayor_guardado, 30);
    assert!(!menores.ofrecer(30));
    assert!(!menores.ofrecer(99));
    assert!(menores.ofrecer(5));
    assert_eq!(menores.mayor_guardado, 20);
}

#[test]
fn contar_compartidos_a_mano() {
    assert_eq!(contar_compartidos(&[1, 3, 5, 7], &[3, 4, 5, 9], 3), (1, 3));
    assert_eq!(contar_compartidos(&[1, 3, 5, 7], &[3, 4, 5, 9], 10), (2, 6));
    assert_eq!(contar_compartidos(&[], &[3, 4], 10), (0, 2));
    assert_eq!(contar_compartidos(&[2, 4], &[2, 4], 10), (2, 2));
}

#[test]
fn distancia_de_mash_a_mano() {
    assert_eq!(distancia_mash(1.0, 21), 0.0);
    assert_eq!(distancia_mash(0.0, 21), 1.0);
    let esperada = (1.5f64).ln() / 21.0;
    assert!((distancia_mash(0.5, 21) - esperada).abs() < 1e-12);
    assert!((distancia_mash(0.5, 21) - 0.019_307_8).abs() < 1e-6);
}

#[test]
fn valor_p_a_mano() {
    assert!((log10_valor_p(1, 3, 0.5) - 0.875f64.log10()).abs() < 1e-12);
    assert!((log10_valor_p(3, 3, 0.5) - 0.125f64.log10()).abs() < 1e-12);
    assert_eq!(log10_valor_p(0, 1000, 1e-6), 0.0);
    assert_eq!(log10_valor_p(5, 10, 1.0), 0.0);
    let diminuto = log10_valor_p(500, 1000, 1e-6);
    assert!(diminuto < -2000.0 && diminuto.is_finite());
}

#[test]
fn valor_p_contra_la_suma_binomial_directa() {
    let mut azar = Azar::nuevo(7);
    for _ in 0..300 {
        let vistos = 1 + azar.entero(60);
        let compartidos = 1 + azar.entero(vistos);
        let probabilidad = 0.001 + 0.9 * azar.decimal();
        let directa = cola_binomial(compartidos, vistos, probabilidad);
        let calculada = 10f64.powf(log10_valor_p(compartidos, vistos, probabilidad));
        assert!(
            (directa - calculada).abs() <= 1e-9 * directa.max(1e-300) + 1e-13,
            "{compartidos}/{vistos} con {probabilidad}: {directa} contra {calculada}"
        );
    }
}

#[test]
fn el_jaccard_al_azar_crece_con_el_genoma() {
    let chico = jaccard_al_azar(5_000, 5_000, 21);
    let grande = jaccard_al_azar(5_000_000, 5_000_000, 21);
    assert!(chico < grande && grande < 1e-5);
    let corto = jaccard_al_azar(5_000_000, 5_000_000, 8);
    assert!(corto > 0.9);
}

#[test]
fn bosquejar_igual_que_la_version_ingenua() {
    let mut azar = Azar::nuevo(11);
    for vuelta in 0..120 {
        let longitud = azar.entero(3000);
        let mut secuencia = al_azar(&mut azar, longitud);
        for base in secuencia.iter_mut() {
            if azar.entero(60) == 0 {
                *base = b'N';
            }
        }
        if vuelta % 7 == 0 {
            secuencia = secuencia.repeat(3);
        }
        let k = 1 + azar.entero(32);
        let cupo = 1 + azar.entero(400);
        let bosquejo = bosquejar(&una_pieza(&secuencia), k, cupo);
        let conjunto = kmeros_ingenuos(&secuencia, k);
        assert_eq!(bosquejo.hashes, bosquejo_ingenuo(&conjunto, cupo), "k={k}");
        let mut exacto: Vec<u64> = conjunto.into_iter().collect();
        exacto.sort();
        assert_eq!(todos_los_kmeros(&una_pieza(&secuencia), k), exacto);
    }
}

#[test]
fn compartidos_igual_que_la_union_ingenua() {
    let mut azar = Azar::nuevo(23);
    for _ in 0..120 {
        let longitud = 500 + azar.entero(4000);
        let base = al_azar(&mut azar, longitud);
        let fraccion = azar.decimal() * 0.2;
        let otra = mutar(&mut azar, &base, fraccion);
        let k = 5 + azar.entero(25);
        let cupo = 1 + azar.entero(1500);
        let a = bosquejar(&una_pieza(&base), k, cupo);
        let b = bosquejar(&una_pieza(&otra), k, cupo);
        let esperado =
            compartidos_ingenuos(&kmeros_ingenuos(&base, k), &kmeros_ingenuos(&otra, k), cupo);
        assert_eq!(contar_compartidos(&a.hashes, &b.hashes, cupo), esperado);
    }
}

#[test]
fn con_pocos_kmeros_el_jaccard_es_exacto() {
    let mut azar = Azar::nuevo(5);
    let a = al_azar(&mut azar, 400);
    let b = mutar(&mut azar, &a, 0.03);
    let resultado = comparar(&[una_pieza(&a), una_pieza(&b)], 21, 1000).unwrap();
    let exacto = jaccard_exacto(
        &todos_los_kmeros(&una_pieza(&a), 21),
        &todos_los_kmeros(&una_pieza(&b), 21),
    );
    assert_eq!(resultado.pares[0].medida.jaccard, exacto);
}

#[test]
fn el_ani_estimado_sigue_a_la_mutacion() {
    let mut azar = Azar::nuevo(2026);
    let original = al_azar(&mut azar, 200_000);
    let uno = mutar(&mut azar, &original, 0.01);
    let cinco = mutar(&mut azar, &original, 0.05);
    let genomas = vec![una_pieza(&original), una_pieza(&uno), una_pieza(&cinco)];
    let resultado = comparar(&genomas, 21, 1000).unwrap();
    let completos: Vec<Vec<u64>> = genomas.iter().map(|g| todos_los_kmeros(g, 21)).collect();
    for par in &resultado.pares {
        let exacto = jaccard_exacto(&completos[par.primero], &completos[par.segundo]);
        let error = (exacto * (1.0 - exacto) / 1000.0).sqrt();
        assert!(
            (par.medida.jaccard - exacto).abs() <= 4.0 * error + 1e-9,
            "{} contra {}",
            par.medida.jaccard,
            exacto
        );
        assert!(par.medida.log10_valor_p < -100.0);
    }
    let ani_uno = resultado.pares[0].medida.identidad;
    let ani_cinco = resultado.pares[1].medida.identidad;
    assert!((ani_uno - 0.99).abs() < 0.003, "{ani_uno}");
    assert!((ani_cinco - 0.95).abs() < 0.01, "{ani_cinco}");
    assert!(ani_cinco < ani_uno);
}

#[test]
fn una_secuencia_y_su_complemento_inverso_son_iguales() {
    let mut azar = Azar::nuevo(3);
    let a = al_azar(&mut azar, 50_000);
    let b = complemento_inverso(&a);
    let resultado = comparar(&[una_pieza(&a), una_pieza(&b)], 21, 1000).unwrap();
    let medida = resultado.pares[0].medida;
    assert_eq!(medida.compartidos, medida.vistos);
    assert_eq!(medida.distancia, 0.0);
    assert_eq!(medida.identidad, 1.0);
}

#[test]
fn secuencias_sin_relacion_tienen_distancia_uno() {
    let mut azar = Azar::nuevo(9);
    let a = al_azar(&mut azar, 50_000);
    let b = al_azar(&mut azar, 50_000);
    let resultado = comparar(&[una_pieza(&a), una_pieza(&b)], 21, 1000).unwrap();
    let medida = resultado.pares[0].medida;
    assert_eq!(medida.compartidos, 0);
    assert_eq!(medida.distancia, 1.0);
    assert_eq!(medida.log10_valor_p, 0.0);
    assert_eq!(
        matriz_de_distancias(&resultado),
        vec![vec![0.0, 1.0], vec![1.0, 0.0]]
    );
}

#[test]
fn la_contencion_ve_un_trozo_dentro_de_un_genoma() {
    let mut azar = Azar::nuevo(77);
    let grande = al_azar(&mut azar, 300_000);
    let chico = grande[100_000..130_000].to_vec();
    let resultado = comparar(&[una_pieza(&chico), una_pieza(&grande)], 21, 1000).unwrap();
    assert_eq!(resultado.contencion[0][1], 1.0);
    assert_eq!(resultado.contencion[0][0], 1.0);
    let exacta = contencion_exacta(
        &todos_los_kmeros(&una_pieza(&grande), 21),
        &todos_los_kmeros(&una_pieza(&chico), 21),
    );
    assert!((resultado.contencion[1][0] - exacta).abs() < 0.05);
    assert!(resultado.pares[0].medida.jaccard < 0.2);
    assert_eq!(identidad_de_contencion(1.0, 21), 1.0);
    let mutado = mutar(&mut azar, &chico, 0.02);
    let otra = comparar(&[una_pieza(&mutado), una_pieza(&grande)], 21, 1000).unwrap();
    let identidad = identidad_de_contencion(otra.contencion[0][1], 21);
    assert!((identidad - 0.98).abs() < 0.004, "{identidad}");
}

#[test]
fn la_contencion_es_la_fraccion_del_bosquejo_presente() {
    let mut azar = Azar::nuevo(31);
    for _ in 0..40 {
        let largo_a = 200 + azar.entero(3000);
        let a = al_azar(&mut azar, largo_a);
        let fraccion = azar.decimal() * 0.1;
        let b = mutar(&mut azar, &a, fraccion);
        let largo_c = 100 + azar.entero(500);
        let c = al_azar(&mut azar, largo_c);
        let k = 8 + azar.entero(20);
        let cupo = 1 + azar.entero(800);
        let genomas = vec![una_pieza(&a), una_pieza(&b), una_pieza(&c)];
        let resultado = comparar(&genomas, k, cupo).unwrap();
        let completos: Vec<HashSet<u64>> = [&a, &b, &c]
            .iter()
            .map(|secuencia| kmeros_ingenuos(secuencia, k))
            .collect();
        for (i, bosquejo) in resultado.bosquejos.iter().enumerate() {
            for (j, completo) in completos.iter().enumerate() {
                let hallados = bosquejo
                    .hashes
                    .iter()
                    .filter(|hash| completo.contains(hash))
                    .count();
                let esperada = hallados as f64 / bosquejo.hashes.len() as f64;
                assert_eq!(resultado.contencion[i][j], esperada);
            }
        }
    }
}

#[test]
fn las_piezas_no_se_pegan() {
    let piezas = vec![b"ACGTACGTAC".to_vec(), b"GGCCTTAAGG".to_vec()];
    let bosquejo = bosquejar(&piezas, 5, 100);
    assert_eq!(bosquejo.kmeros, 12);
    assert_eq!(bosquejo.longitud, 20);
    let pegado = bosquejar(&una_pieza(b"ACGTACGTACGGCCTTAAGG"), 5, 100);
    assert_eq!(pegado.kmeros, 16);
}

#[test]
fn estima_cuantos_kmeros_distintos_hay() {
    let mut azar = Azar::nuevo(41);
    let secuencia = al_azar(&mut azar, 100_000);
    let bosquejo = bosquejar(&una_pieza(&secuencia), 21, 1000);
    let estimado = estimar_kmeros_distintos(&bosquejo);
    let exacto = todos_los_kmeros(&una_pieza(&secuencia), 21).len() as f64;
    assert!(
        (estimado / exacto - 1.0).abs() < 0.12,
        "{estimado} contra {exacto}"
    );
    let chico = bosquejar(&una_pieza(b"ACGTACGTTT"), 3, 1000);
    assert_eq!(estimar_kmeros_distintos(&chico), chico.hashes.len() as f64);
}

#[test]
fn errores_en_espanol() {
    let a = una_pieza(b"ACGTACGTACGTACGTACGTACGTAAAC");
    let b = una_pieza(b"ACGTACGTACGTACGTACGTACGTAAAG");
    assert!(comparar(&[a.clone(), b.clone()], 0, 1000).is_err());
    assert!(comparar(&[a.clone(), b.clone()], 33, 1000).is_err());
    assert!(comparar(&[a.clone(), b.clone()], 21, 0).is_err());
    assert!(comparar(&[a.clone(), b.clone()], 21, 100_001).is_err());
    assert!(comparar(std::slice::from_ref(&a), 21, 1000).is_err());
    assert!(comparar(&[], 21, 1000).is_err());
    let corta = una_pieza(b"ACGTNNACGT");
    let error = comparar(&[a.clone(), corta], 21, 1000).unwrap_err();
    assert!(error.contains("número 2"));
    let muchas = vec![a.clone(); 101];
    assert!(comparar(&muchas, 21, 1000).is_err());
    assert!(comparar(&[a, b], 21, 1000).is_ok());
}

#[test]
fn bosquejar_el_genoma_de_ecoli_es_rapido() {
    let genoma = nucleo::ecoli::genoma().to_vec();
    let inicio = std::time::Instant::now();
    let bosquejo = bosquejar(&[genoma], 21, 1000);
    assert!(inicio.elapsed().as_secs_f64() < 5.0);
    assert_eq!(bosquejo.hashes.len(), 1000);
    let estimado = estimar_kmeros_distintos(&bosquejo);
    assert!(estimado > 4.0e6 && estimado < 5.4e6, "{estimado}");
}

#[test]
fn validar_limita_los_hashes_en_total() {
    assert!(nucleo::minhash::validar(100, 21, 10_000).is_ok());
    assert!(nucleo::minhash::validar(100, 21, 10_001).is_err());
    assert!(nucleo::minhash::validar(10, 21, 100_000).is_ok());
    assert!(nucleo::minhash::validar(11, 21, 100_000).is_err());
}
