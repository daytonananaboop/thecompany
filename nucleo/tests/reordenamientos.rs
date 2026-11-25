use nucleo::azar::Azar;
use nucleo::ecoli;
use nucleo::reordenamientos::*;
use nucleo::secuencia::complemento_inverso;

fn permutacion_al_azar(azar: &mut Azar, longitud: usize) -> Vec<i64> {
    let mut bloques: Vec<i64> = (1..=longitud as i64).collect();
    for indice in (1..longitud).rev() {
        let otro = azar.entero(indice + 1);
        bloques.swap(indice, otro);
    }
    for bloque in bloques.iter_mut() {
        if azar.entero(2) == 1 {
            *bloque = -*bloque;
        }
    }
    bloques
}

fn genoma_al_azar(azar: &mut Azar, longitud: usize) -> Vec<Vec<i64>> {
    let bloques = permutacion_al_azar(azar, longitud);
    let mut genoma = vec![Vec::new()];
    for bloque in bloques {
        let ultimo = genoma.len() - 1;
        if !genoma[ultimo].is_empty() && azar.entero(4) == 0 {
            genoma.push(Vec::new());
        }
        let ultimo = genoma.len() - 1;
        genoma[ultimo].push(bloque);
    }
    genoma
}

fn aristas_normalizadas(aristas: &[(usize, usize)]) -> Vec<(usize, usize)> {
    let mut normalizadas: Vec<(usize, usize)> =
        aristas.iter().map(|&(a, b)| (a.min(b), a.max(b))).collect();
    normalizadas.sort();
    normalizadas
}

fn mismas_aristas(p: &[Vec<i64>], q: &[Vec<i64>]) -> bool {
    aristas_normalizadas(&aristas_de_color(p).unwrap())
        == aristas_normalizadas(&aristas_de_color(q).unwrap())
}

fn adn_al_azar(azar: &mut Azar, longitud: usize) -> Vec<u8> {
    (0..longitud).map(|_| b"ACGT"[azar.entero(4)]).collect()
}

#[test]
fn voraz_ordena_un_caso_corto() {
    let pasos = ordenamiento_voraz(&[2, -1, 3]).unwrap();
    assert_eq!(pasos, vec![vec![1, -2, 3], vec![1, 2, 3]]);
}

#[test]
fn voraz_voltea_signos_por_separado() {
    let pasos = ordenamiento_voraz(&[-3, 1, 2]).unwrap();
    assert_eq!(
        pasos,
        vec![
            vec![-1, 3, 2],
            vec![1, 3, 2],
            vec![1, -2, -3],
            vec![1, 2, -3],
            vec![1, 2, 3],
        ]
    );
}

#[test]
fn voraz_con_identidad_y_un_bloque() {
    assert!(ordenamiento_voraz(&[1, 2, 3]).unwrap().is_empty());
    assert_eq!(ordenamiento_voraz(&[-1]).unwrap(), vec![vec![1]]);
}

#[test]
fn voraz_rechaza_permutaciones_invalidas() {
    assert!(ordenamiento_voraz(&[]).is_err());
    assert!(ordenamiento_voraz(&[0]).is_err());
    assert!(ordenamiento_voraz(&[1, 1]).is_err());
    assert!(ordenamiento_voraz(&[1, -1]).is_err());
    assert!(ordenamiento_voraz(&[1, 3]).is_err());
    assert!(ordenamiento_voraz(&[i64::MIN]).is_err());
    assert!(ordenamiento_voraz(&[i64::MAX, 1]).is_err());
}

#[test]
fn voraz_siempre_termina_en_la_identidad_con_inversiones() {
    let mut azar = Azar::nuevo(6);
    for longitud in 1..40 {
        let inicial = permutacion_al_azar(&mut azar, longitud);
        let pasos = ordenamiento_voraz(&inicial).unwrap();
        assert!(pasos.len() <= 2 * longitud);
        let identidad: Vec<i64> = (1..=longitud as i64).collect();
        if pasos.is_empty() {
            assert_eq!(inicial, identidad);
            continue;
        }
        assert_eq!(pasos.last().unwrap(), &identidad);
        let mut anterior = inicial.clone();
        for siguiente in &pasos {
            let inicio = (0..longitud)
                .find(|&i| anterior[i] != siguiente[i])
                .unwrap();
            let fin = (0..longitud)
                .rfind(|&i| anterior[i] != siguiente[i])
                .unwrap();
            for desplazamiento in 0..=fin - inicio {
                assert_eq!(
                    siguiente[inicio + desplazamiento],
                    -anterior[fin - desplazamiento]
                );
            }
            anterior = siguiente.clone();
        }
    }
}

#[test]
fn cuenta_rupturas() {
    assert_eq!(numero_de_rupturas(&[2, 3, -1, 4]).unwrap(), 3);
    assert_eq!(numero_de_rupturas(&[-2, -1]).unwrap(), 2);
    assert_eq!(numero_de_rupturas(&[1, 2, 3, 4, 5]).unwrap(), 0);
    assert_eq!(numero_de_rupturas(&[-1]).unwrap(), 2);
    assert_eq!(numero_de_rupturas(&[3, 2, 1]).unwrap(), 4);
    assert!(numero_de_rupturas(&[]).is_err());
    assert!(numero_de_rupturas(&[2]).is_err());
}

#[test]
fn lista_las_rupturas_con_sus_vecinos() {
    assert_eq!(rupturas(&[-2, -1]).unwrap(), vec![(0, -2), (-1, 3)]);
    assert_eq!(
        rupturas(&[2, 3, -1, 4]).unwrap(),
        vec![(0, 2), (3, -1), (-1, 4)]
    );
}

#[test]
fn las_rupturas_acotan_al_voraz() {
    let mut azar = Azar::nuevo(11);
    for longitud in 1..30 {
        let permutacion = permutacion_al_azar(&mut azar, longitud);
        let rupturas = numero_de_rupturas(&permutacion).unwrap();
        assert!(rupturas <= longitud + 1);
        let pasos = ordenamiento_voraz(&permutacion).unwrap();
        assert!(pasos.len() >= rupturas.div_ceil(2));
        if let Some(ultima) = pasos.last() {
            assert_eq!(numero_de_rupturas(ultima).unwrap(), 0);
        }
    }
}

#[test]
fn cromosoma_a_ciclo_y_vuelta() {
    assert_eq!(
        cromosoma_a_ciclo(&[2, -1, -3]).unwrap(),
        vec![3, 4, 2, 1, 6, 5]
    );
    assert_eq!(
        ciclo_a_cromosoma(&[3, 4, 2, 1, 6, 5]).unwrap(),
        vec![2, -1, -3]
    );
    assert_eq!(ciclo_a_cromosoma(&[2, 1]).unwrap(), vec![-1]);
    let mut azar = Azar::nuevo(3);
    for longitud in 1..25 {
        let cromosoma = permutacion_al_azar(&mut azar, longitud);
        let ciclo = cromosoma_a_ciclo(&cromosoma).unwrap();
        assert_eq!(ciclo.len(), 2 * longitud);
        assert_eq!(ciclo_a_cromosoma(&ciclo).unwrap(), cromosoma);
    }
}

#[test]
fn ciclo_a_cromosoma_rechaza_ciclos_malos() {
    assert!(ciclo_a_cromosoma(&[]).is_err());
    assert!(ciclo_a_cromosoma(&[1]).is_err());
    assert!(ciclo_a_cromosoma(&[1, 3]).is_err());
    assert!(ciclo_a_cromosoma(&[2, 3]).is_err());
    assert!(ciclo_a_cromosoma(&[1, 2, 1, 2]).is_err());
    assert!(ciclo_a_cromosoma(&[3, 4]).is_err());
    assert!(ciclo_a_cromosoma(&[0, 1]).is_err());
    assert!(ciclo_a_cromosoma(&[usize::MAX, usize::MAX - 1]).is_err());
    assert!(cromosoma_a_ciclo(&[0]).is_err());
    assert!(cromosoma_a_ciclo(&[]).is_err());
}

#[test]
fn aristas_de_color_de_dos_cromosomas() {
    let genoma = vec![vec![2, -1, -3], vec![4]];
    assert_eq!(
        aristas_de_color(&genoma).unwrap(),
        vec![(4, 2), (1, 6), (5, 3), (8, 7)]
    );
    assert!(aristas_de_color(&[]).is_err());
    assert!(aristas_de_color(&[vec![1], vec![]]).is_err());
    assert!(aristas_de_color(&[vec![1, 2], vec![2]]).is_err());
    assert!(aristas_de_color(&[vec![1, 3]]).is_err());
}

#[test]
fn grafo_a_genoma_empieza_cada_cromosoma_por_su_menor_bloque() {
    let genoma = grafo_a_genoma(&[(4, 2), (1, 6), (5, 3), (8, 7)]).unwrap();
    assert_eq!(genoma, vec![vec![1, -2, 3], vec![4]]);
    assert_eq!(grafo_a_genoma(&[(2, 1)]).unwrap(), vec![vec![1]]);
}

#[test]
fn grafo_a_genoma_rechaza_grafos_malos() {
    assert!(grafo_a_genoma(&[]).is_err());
    assert!(grafo_a_genoma(&[(1, 1)]).is_err());
    assert!(grafo_a_genoma(&[(1, 3)]).is_err());
    assert!(grafo_a_genoma(&[(0, 1)]).is_err());
    assert!(grafo_a_genoma(&[(1, 2), (1, 3)]).is_err());
    assert!(grafo_a_genoma(&[(usize::MAX, 1)]).is_err());
}

#[test]
fn ida_y_vuelta_entre_genoma_y_grafo() {
    let mut azar = Azar::nuevo(21);
    for longitud in 1..40 {
        let genoma = genoma_al_azar(&mut azar, longitud);
        let aristas = aristas_de_color(&genoma).unwrap();
        assert_eq!(aristas.len(), longitud);
        let reconstruido = grafo_a_genoma(&aristas).unwrap();
        assert_eq!(reconstruido.len(), genoma.len());
        assert!(mismas_aristas(&genoma, &reconstruido));
        for cromosoma in &reconstruido {
            let menor = cromosoma.iter().map(|b| b.abs()).min().unwrap();
            assert_eq!(cromosoma[0], menor);
        }
    }
}

#[test]
fn dos_roturas_en_el_grafo() {
    let aristas = [(2, 4), (3, 8), (7, 5), (6, 1)];
    assert_eq!(
        dos_roturas_en_grafo(&aristas, 1, 6, 3, 8).unwrap(),
        vec![(2, 4), (6, 8), (7, 5), (1, 3)]
    );
    assert!(dos_roturas_en_grafo(&aristas, 1, 2, 3, 8).is_err());
    assert!(dos_roturas_en_grafo(&aristas, 1, 6, 1, 6).is_err());
    assert!(dos_roturas_en_grafo(&aristas, 1, 6, 99, 8).is_err());
}

#[test]
fn dos_roturas_en_el_genoma_parte_un_cromosoma() {
    let genoma = vec![vec![1, -2, -4, 3]];
    assert_eq!(
        aristas_de_color(&genoma).unwrap(),
        vec![(2, 4), (3, 8), (7, 5), (6, 1)]
    );
    assert_eq!(
        dos_roturas_en_genoma(&genoma, 1, 6, 3, 8).unwrap(),
        vec![vec![1, -2], vec![3, -4]]
    );
    assert_eq!(
        dos_roturas_en_genoma(&[vec![1, -2], vec![3, -4]], 1, 3, 6, 8).unwrap(),
        vec![vec![1, -2, -4, 3]]
    );
    assert!(dos_roturas_en_genoma(&genoma, 1, 2, 3, 4).is_err());
    assert!(dos_roturas_en_genoma(&[vec![0]], 1, 2, 3, 4).is_err());
}

#[test]
fn distancia_de_dos_roturas() {
    let p = vec![vec![1, 2, 3, 4]];
    assert_eq!(
        distancia_dos_roturas(&p, &[vec![1, -3], vec![2, 4]]).unwrap(),
        3
    );
    assert_eq!(
        distancia_dos_roturas(&p, &[vec![1, 2], vec![3, 4]]).unwrap(),
        1
    );
    assert_eq!(
        distancia_dos_roturas(&p, &[vec![-2, -1, -4, -3]]).unwrap(),
        0
    );
    assert_eq!(distancia_dos_roturas(&p, &p).unwrap(), 0);
    assert!(distancia_dos_roturas(&[vec![1, 2]], &[vec![1, 2, 3]]).is_err());
    assert!(distancia_dos_roturas(&[vec![1, 2]], &[vec![1, 1]]).is_err());
    let ciclos = ciclos_del_grafo_de_rupturas(&p, &[vec![1, 2], vec![3, 4]]).unwrap();
    assert_eq!(ciclos, vec![vec![1, 8, 5, 4], vec![2, 3], vec![6, 7]]);
}

#[test]
fn transformacion_mas_corta_paso_a_paso() {
    let p = vec![vec![1, 2, 3, 4]];
    let q = vec![vec![1, -3], vec![2, 4]];
    let pasos = transformacion_mas_corta(&p, &q).unwrap();
    let nodos: Vec<[usize; 4]> = pasos.iter().map(|paso| paso.nodos).collect();
    assert_eq!(nodos, vec![[1, 8, 5, 4], [2, 3, 6, 7], [3, 7, 8, 4]]);
    let genomas: Vec<Vec<Vec<i64>>> = pasos.iter().map(|paso| paso.genoma.clone()).collect();
    assert_eq!(
        genomas,
        vec![
            vec![vec![1, 2, -4, -3]],
            vec![vec![1, -3], vec![2, -4]],
            vec![vec![1, -3], vec![2, 4]],
        ]
    );
    assert!(transformacion_mas_corta(&p, &p).unwrap().is_empty());
    assert!(transformacion_mas_corta(&p, &[vec![1, 2, 3]]).is_err());
}

#[test]
fn transformacion_mas_corta_usa_exactamente_la_distancia() {
    let mut azar = Azar::nuevo(44);
    for longitud in 1..40 {
        let p = genoma_al_azar(&mut azar, longitud);
        let q = genoma_al_azar(&mut azar, longitud);
        let distancia = distancia_dos_roturas(&p, &q).unwrap();
        assert_eq!(distancia, distancia_dos_roturas(&q, &p).unwrap());
        let pasos = transformacion_mas_corta(&p, &q).unwrap();
        assert_eq!(pasos.len(), distancia);
        let mut anterior = p.clone();
        for (numero, paso) in pasos.iter().enumerate() {
            let [i1, i2, i3, i4] = paso.nodos;
            assert_eq!(
                dos_roturas_en_genoma(&anterior, i1, i2, i3, i4).unwrap(),
                paso.genoma
            );
            assert_eq!(
                distancia_dos_roturas(&paso.genoma, &q).unwrap(),
                distancia - numero - 1
            );
            anterior = paso.genoma.clone();
        }
        assert!(mismas_aristas(&anterior, &q));
    }
}

#[test]
fn kmeros_compartidos_en_ambas_hebras() {
    let pares = kmeros_compartidos(b"CATGCAT", b"ATGCCAT", 3).unwrap();
    let vistos: Vec<(usize, usize, bool)> = pares
        .iter()
        .map(|par| (par.posicion_a, par.posicion_b, par.invertido))
        .collect();
    assert_eq!(
        vistos,
        vec![
            (0, 0, true),
            (0, 4, false),
            (1, 0, false),
            (1, 4, true),
            (2, 1, false),
            (3, 1, true),
            (4, 0, true),
            (4, 4, false),
        ]
    );
}

#[test]
fn kmeros_compartidos_con_la_hebra_contraria() {
    let pares = kmeros_compartidos(b"GATTACA", b"TGTAATC", 3).unwrap();
    let vistos: Vec<(usize, usize, bool)> = pares
        .iter()
        .map(|par| (par.posicion_a, par.posicion_b, par.invertido))
        .collect();
    assert_eq!(
        vistos,
        vec![
            (0, 4, true),
            (1, 3, true),
            (2, 2, true),
            (3, 1, true),
            (4, 0, true)
        ]
    );
}

#[test]
fn kmeros_compartidos_saltan_letras_raras_y_casos_borde() {
    let pares = kmeros_compartidos(b"ACGNACG", b"ACG", 3).unwrap();
    let vistos: Vec<(usize, usize, bool)> = pares
        .iter()
        .map(|par| (par.posicion_a, par.posicion_b, par.invertido))
        .collect();
    assert_eq!(vistos, vec![(0, 0, false), (4, 0, false)]);
    assert!(kmeros_compartidos(b"ACGT", b"ACGT", 0).is_err());
    assert!(kmeros_compartidos(b"ACGT", b"ACGT", 33).is_err());
    assert!(kmeros_compartidos(b"ACG", b"ACGT", 5).unwrap().is_empty());
    assert!(kmeros_compartidos(b"", b"", 3).unwrap().is_empty());
    let muchos = vec![b'A'; 3000];
    assert!(kmeros_compartidos(&muchos, &muchos, 1).is_err());
}

#[test]
fn kmeros_de_32_contra_el_complemento_inverso() {
    let a = ecoli::region(1, 40).unwrap();
    let b = complemento_inverso(a);
    let pares = kmeros_compartidos(a, &b, 32).unwrap();
    assert_eq!(pares.len(), 9);
    for (x, par) in pares.iter().enumerate() {
        assert_eq!(
            (par.posicion_a, par.posicion_b, par.invertido),
            (x, 8 - x, true)
        );
    }
}

#[test]
fn kmeros_compartidos_coinciden_con_la_comparacion_directa() {
    let mut azar = Azar::nuevo(9);
    let a = adn_al_azar(&mut azar, 300);
    let b = adn_al_azar(&mut azar, 250);
    let k = 4;
    let mut esperados = Vec::new();
    for x in 0..=a.len() - k {
        for y in 0..=b.len() - k {
            let kmero_a = &a[x..x + k];
            let kmero_b = &b[y..y + k];
            if kmero_a == kmero_b {
                esperados.push((x, y, false));
            } else if kmero_a == complemento_inverso(kmero_b).as_slice() {
                esperados.push((x, y, true));
            }
        }
    }
    let obtenidos: Vec<(usize, usize, bool)> = kmeros_compartidos(&a, &b, k)
        .unwrap()
        .iter()
        .map(|par| (par.posicion_a, par.posicion_b, par.invertido))
        .collect();
    assert_eq!(obtenidos, esperados);
}

#[test]
fn bloques_de_sintenia_encuentran_una_inversion() {
    let mut azar = Azar::nuevo(17);
    let x = adn_al_azar(&mut azar, 1000);
    let y = adn_al_azar(&mut azar, 1000);
    let z = adn_al_azar(&mut azar, 1000);
    let a = [x.clone(), y.clone(), z.clone()].concat();
    let b = [z, complemento_inverso(&y), x].concat();
    let k = 16;
    let sintenia = bloques_de_sintenia(&a, &b, k, 50, 20).unwrap();
    let esperados = [
        (0, 999, 2000, 2999, false),
        (1000, 1999, 1000, 1999, true),
        (2000, 2999, 0, 999, false),
    ];
    assert_eq!(sintenia.bloques.len(), esperados.len());
    let cerca = |obtenido: usize, esperado: usize| obtenido.abs_diff(esperado) < k;
    for (bloque, esperado) in sintenia.bloques.iter().zip(esperados) {
        assert!(cerca(bloque.inicio_a, esperado.0), "{bloque:?}");
        assert!(cerca(bloque.fin_a, esperado.1), "{bloque:?}");
        assert!(cerca(bloque.inicio_b, esperado.2), "{bloque:?}");
        assert!(cerca(bloque.fin_b, esperado.3), "{bloque:?}");
        assert_eq!(bloque.invertido, esperado.4);
        assert!(
            bloque.kmeros >= 985 && bloque.kmeros < 985 + k,
            "{bloque:?}"
        );
    }
    assert_eq!(permutacion_de_bloques(&sintenia.bloques), vec![3, -2, 1]);
    assert!(sintenia.pares.len() >= 3 * 985);
}

#[test]
fn bloques_de_sintenia_casos_borde() {
    assert!(bloques_de_sintenia(b"ACGT", b"ACGT", 2, 0, 1).is_err());
    assert!(bloques_de_sintenia(b"ACGT", b"ACGT", 0, 10, 1).is_err());
    let vacia = bloques_de_sintenia(b"ACGT", b"TTTT", 3, 10, 1).unwrap();
    assert!(vacia.bloques.is_empty());
    let todo = bloques_de_sintenia(b"ACGTACGT", b"ACGTACGT", 3, usize::MAX, 0).unwrap();
    assert!(!todo.bloques.is_empty());
    assert!(permutacion_de_bloques(&[]).is_empty());
}

fn bloque_de_un_par_directo(par: &ParCompartido, bloques: &[BloqueSintenico]) -> usize {
    for (indice, bloque) in bloques.iter().enumerate() {
        let dentro = par.invertido == bloque.invertido
            && (bloque.inicio_a..=bloque.fin_a).contains(&par.posicion_a)
            && (bloque.inicio_b..=bloque.fin_b).contains(&par.posicion_b);
        if dentro {
            return indice + 1;
        }
    }
    0
}

#[test]
fn bloque_de_cada_par_coincide_con_la_busqueda_directa() {
    let mut azar = Azar::nuevo(23);
    for ronda in 0..30 {
        let x = adn_al_azar(&mut azar, 200);
        let y = adn_al_azar(&mut azar, 200);
        let a = [x.clone(), y.clone(), adn_al_azar(&mut azar, 50)].concat();
        let b = [complemento_inverso(&y), adn_al_azar(&mut azar, 30), x].concat();
        let k = 6 + ronda % 5;
        let sintenia = bloques_de_sintenia(&a, &b, k, 5 + ronda, ronda % 4).unwrap();
        let obtenidos = bloque_de_cada_par(&sintenia.pares, &sintenia.bloques);
        assert_eq!(obtenidos.len(), sintenia.pares.len());
        for (par, obtenido) in sintenia.pares.iter().zip(obtenidos) {
            assert_eq!(obtenido, bloque_de_un_par_directo(par, &sintenia.bloques));
        }
    }
    assert!(bloque_de_cada_par(&[], &[]).is_empty());
}

#[test]
fn bloque_de_cada_par_es_rapido_con_muchos_bloques() {
    let mut azar = Azar::nuevo(5);
    let a = adn_al_azar(&mut azar, 300_000);
    let b = adn_al_azar(&mut azar, 300_000);
    let sintenia = bloques_de_sintenia(&a, &b, 9, 1, 0).unwrap();
    assert!(sintenia.bloques.len() > 100_000);
    let muestra: Vec<ParCompartido> = sintenia.pares.iter().step_by(100).copied().collect();
    let inicio = std::time::Instant::now();
    let obtenidos = bloque_de_cada_par(&muestra, &sintenia.bloques);
    assert!(inicio.elapsed().as_secs() < 2);
    assert!(obtenidos.iter().all(|&bloque| bloque > 0));
}

#[test]
fn el_error_de_demasiadas_parejas_no_pide_un_k_mayor_cuando_k_es_32() {
    let repetida = vec![b'A'; 2000];
    let con_32 = kmeros_compartidos(&repetida, &repetida, 32).unwrap_err();
    assert!(!con_32.contains("k mayor"), "{con_32}");
    let con_10 = kmeros_compartidos(&repetida, &repetida, 10).unwrap_err();
    assert!(con_10.contains("k mayor"), "{con_10}");
}
