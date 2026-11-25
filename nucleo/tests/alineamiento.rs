use nucleo::alineamiento::*;
use nucleo::azar::Azar;
use nucleo::matrices::{blosum62, pam250};

fn sin_huecos(fila: &[u8]) -> Vec<u8> {
    fila.iter()
        .copied()
        .filter(|&letra| letra != b'-')
        .collect()
}

fn al_azar(azar: &mut Azar, longitud_maxima: usize, letras: &[u8]) -> Vec<u8> {
    let longitud = azar.entero(longitud_maxima + 1);
    (0..longitud)
        .map(|_| letras[azar.entero(letras.len())])
        .collect()
}

fn global_simple(v: &[u8], w: &[u8], puntuar: Puntuacion, indel: i64) -> i64 {
    let mut tabla = vec![vec![0i64; w.len() + 1]; v.len() + 1];
    for (i, letra_v) in std::iter::once(0u8).chain(v.iter().copied()).enumerate() {
        for (j, letra_w) in std::iter::once(0u8).chain(w.iter().copied()).enumerate() {
            if i == 0 && j == 0 {
                continue;
            }
            let mut mejor = i64::MIN;
            if i > 0 && j > 0 {
                mejor = mejor.max(tabla[i - 1][j - 1] + puntuar(letra_v, letra_w) as i64);
            }
            if i > 0 {
                mejor = mejor.max(tabla[i - 1][j] - indel);
            }
            if j > 0 {
                mejor = mejor.max(tabla[i][j - 1] - indel);
            }
            tabla[i][j] = mejor;
        }
    }
    tabla[v.len()][w.len()]
}

fn afin_simple(v: &[u8], w: &[u8], puntuar: Puntuacion, apertura: i64, extension: i64) -> i64 {
    let nada = i64::MIN / 4;
    let mut parejas = vec![vec![nada; w.len() + 1]; v.len() + 1];
    let mut huecos_w = vec![vec![nada; w.len() + 1]; v.len() + 1];
    let mut huecos_v = vec![vec![nada; w.len() + 1]; v.len() + 1];
    parejas[0][0] = 0;
    for i in 0..=v.len() {
        for j in 0..=w.len() {
            if i > 0 {
                huecos_w[i][j] = (parejas[i - 1][j] - apertura)
                    .max(huecos_w[i - 1][j] - extension)
                    .max(huecos_v[i - 1][j] - apertura);
            }
            if j > 0 {
                huecos_v[i][j] = (parejas[i][j - 1] - apertura)
                    .max(huecos_v[i][j - 1] - extension)
                    .max(huecos_w[i][j - 1] - apertura);
            }
            if i > 0 && j > 0 {
                let antes = parejas[i - 1][j - 1]
                    .max(huecos_w[i - 1][j - 1])
                    .max(huecos_v[i - 1][j - 1]);
                parejas[i][j] = antes + puntuar(v[i - 1], w[j - 1]) as i64;
            }
        }
    }
    let (i, j) = (v.len(), w.len());
    parejas[i][j].max(huecos_w[i][j]).max(huecos_v[i][j])
}

fn comun_simple(v: &[u8], w: &[u8]) -> usize {
    let mut tabla = vec![vec![0usize; w.len() + 1]; v.len() + 1];
    for i in 1..=v.len() {
        for j in 1..=w.len() {
            tabla[i][j] = if v[i - 1] == w[j - 1] {
                tabla[i - 1][j - 1] + 1
            } else {
                tabla[i - 1][j].max(tabla[i][j - 1])
            };
        }
    }
    tabla[v.len()][w.len()]
}

fn comun_de_tres_simple(u: &[u8], v: &[u8], w: &[u8]) -> usize {
    let mut tabla = vec![vec![vec![0usize; w.len() + 1]; v.len() + 1]; u.len() + 1];
    for i in 0..=u.len() {
        for j in 0..=v.len() {
            for k in 0..=w.len() {
                let mut mejor = 0;
                if i > 0 {
                    mejor = mejor.max(tabla[i - 1][j][k]);
                }
                if j > 0 {
                    mejor = mejor.max(tabla[i][j - 1][k]);
                }
                if k > 0 {
                    mejor = mejor.max(tabla[i][j][k - 1]);
                }
                if i > 0 && j > 0 && k > 0 && u[i - 1] == v[j - 1] && v[j - 1] == w[k - 1] {
                    mejor = mejor.max(tabla[i - 1][j - 1][k - 1] + 1);
                }
                tabla[i][j][k] = mejor;
            }
        }
    }
    tabla[u.len()][v.len()][w.len()]
}

fn es_subsecuencia(corta: &[u8], larga: &[u8]) -> bool {
    let mut restantes = larga.iter();
    corta
        .iter()
        .all(|letra| restantes.any(|candidata| candidata == letra))
}

fn revisar(
    alineamiento: &Alineamiento,
    v: &[u8],
    w: &[u8],
    puntuar: Puntuacion,
    indel: i64,
) -> (usize, usize) {
    assert_eq!(alineamiento.fila_v.len(), alineamiento.fila_w.len());
    let parte_v = sin_huecos(&alineamiento.fila_v);
    let parte_w = sin_huecos(&alineamiento.fila_w);
    let fin_v = alineamiento.inicio_v + parte_v.len();
    let fin_w = alineamiento.inicio_w + parte_w.len();
    assert_eq!(&v[alineamiento.inicio_v..fin_v], parte_v.as_slice());
    assert_eq!(&w[alineamiento.inicio_w..fin_w], parte_w.as_slice());
    assert_eq!(
        puntuar_alineamiento(&alineamiento.fila_v, &alineamiento.fila_w, puntuar, indel),
        alineamiento.puntaje
    );
    (fin_v, fin_w)
}

#[test]
fn cambio_minimo_supera_a_la_regla_codiciosa() {
    let monedas = [1, 5, 10, 20, 25];
    let mut usadas = cambio_minimo(40, &monedas).unwrap();
    usadas.sort();
    assert_eq!(usadas, vec![20, 20]);
    let usadas = cambio_minimo(11, &[1, 5, 7]).unwrap();
    assert_eq!(usadas.len(), 3);
    assert_eq!(usadas.iter().sum::<usize>(), 11);
    assert_eq!(
        tabla_de_cambio(6, &[1, 3, 4]),
        vec![
            Some(0),
            Some(1),
            Some(2),
            Some(1),
            Some(1),
            Some(2),
            Some(2)
        ]
    );
}

#[test]
fn cambio_imposible_o_trivial() {
    assert_eq!(cambio_minimo(7, &[2, 4]), None);
    assert_eq!(cambio_minimo(0, &[3]), Some(vec![]));
    assert_eq!(cambio_minimo(5, &[0, 5]), Some(vec![5]));
    assert_eq!(cambio_minimo(5, &[]), None);
    let tabla = tabla_de_cambio(5, &[2]);
    assert_eq!(tabla, vec![Some(0), None, Some(1), None, Some(2), None]);
}

#[test]
fn turista_de_manhattan_toma_el_camino_mas_pesado() {
    let abajo = vec![vec![1.0, 0.0, 2.0], vec![4.0, 3.0, 1.0]];
    let derecha = vec![vec![3.0, 2.0], vec![0.0, 4.0], vec![1.0, 5.0]];
    let recorrido = turista_de_manhattan(&abajo, &derecha).unwrap();
    assert_eq!(recorrido.longitud, 11.0);
    assert_eq!(recorrido.mejor[1][2], 7.0);
    assert_eq!(recorrido.mejor[2][1], 6.0);
    assert_eq!(recorrido.camino.first(), Some(&(0, 0)));
    assert_eq!(recorrido.camino.last(), Some(&(2, 2)));
    let mut suma = 0.0;
    for par in recorrido.camino.windows(2) {
        let ((i, j), (siguiente_i, siguiente_j)) = (par[0], par[1]);
        if siguiente_i == i + 1 && siguiente_j == j {
            suma += abajo[i][j];
        } else {
            assert_eq!((siguiente_i, siguiente_j), (i, j + 1));
            suma += derecha[i][j];
        }
    }
    assert_eq!(suma, 11.0);
}

#[test]
fn turista_con_cuadricula_minima_y_matrices_que_no_encajan() {
    let recorrido = turista_de_manhattan(&[], &[vec![]]).unwrap();
    assert_eq!(recorrido.longitud, 0.0);
    assert_eq!(recorrido.camino, vec![(0, 0)]);
    let abajo = vec![vec![1.0, 2.0]];
    assert!(turista_de_manhattan(&abajo, &[vec![1.0]]).is_err());
    assert!(turista_de_manhattan(&abajo, &[]).is_err());
    assert!(turista_de_manhattan(&abajo, &[vec![1.0], vec![2.0, 3.0]]).is_err());
}

#[test]
fn subsecuencia_comun_mas_larga_a_mano() {
    let comun = subsecuencia_comun_mas_larga(b"TTAGCCGA", b"TAGGCA");
    assert_eq!(comun.len(), 5);
    assert!(es_subsecuencia(&comun, b"TTAGCCGA"));
    assert!(es_subsecuencia(&comun, b"TAGGCA"));
    assert_eq!(subsecuencia_comun_mas_larga(b"ACGT", b"TGCA").len(), 1);
    assert_eq!(subsecuencia_comun_mas_larga(b"", b"ACGT"), b"");
    assert_eq!(subsecuencia_comun_mas_larga(b"AAAA", b"CCCC"), b"");
    assert_eq!(
        subsecuencia_comun_mas_larga(b"GATTACA", b"GATTACA"),
        b"GATTACA"
    );
}

#[test]
fn orden_topologico_respeta_cada_arista() {
    let aristas = [(0, 1), (1, 2), (0, 2), (3, 2), (4, 0)];
    let orden = orden_topologico(6, &aristas).unwrap();
    assert_eq!(orden.len(), 6);
    let mut lugar = [0; 6];
    for (posicion, &nodo) in orden.iter().enumerate() {
        lugar[nodo] = posicion;
    }
    for &(desde, hasta) in &aristas {
        assert!(lugar[desde] < lugar[hasta]);
    }
    assert_eq!(orden_topologico(0, &[]).unwrap(), Vec::<usize>::new());
}

#[test]
fn orden_topologico_detecta_ciclos_y_nodos_raros() {
    assert!(orden_topologico(3, &[(0, 1), (1, 2), (2, 0)]).is_err());
    assert!(orden_topologico(2, &[(1, 1)]).is_err());
    assert!(orden_topologico(2, &[(0, 5)]).is_err());
}

#[test]
fn camino_mas_largo_en_un_grafo_aciclico() {
    let aristas = [
        (0, 1, 2.0),
        (0, 2, 1.0),
        (1, 3, 4.0),
        (2, 3, 7.0),
        (3, 4, 1.0),
        (1, 4, 3.0),
    ];
    let camino = camino_mas_largo(5, &aristas, 0, 4).unwrap();
    assert_eq!(camino.longitud, 9.0);
    assert_eq!(camino.nodos, vec![0, 2, 3, 4]);
    let quieto = camino_mas_largo(5, &aristas, 3, 3).unwrap();
    assert_eq!(quieto.longitud, 0.0);
    assert_eq!(quieto.nodos, vec![3]);
    let negativo = camino_mas_largo(3, &[(0, 1, -2.0), (1, 2, -3.0), (0, 2, -6.0)], 0, 2).unwrap();
    assert_eq!(negativo.longitud, -5.0);
}

#[test]
fn camino_mas_largo_avisa_de_errores() {
    let aristas = [(0, 1, 1.0), (2, 1, 1.0)];
    assert!(camino_mas_largo(3, &aristas, 0, 2).is_err());
    assert!(camino_mas_largo(3, &aristas, 0, 7).is_err());
    assert!(camino_mas_largo(2, &[(0, 1, 1.0), (1, 0, 1.0)], 0, 1).is_err());
}

#[test]
fn alineamiento_global_a_mano() {
    let adn = alineamiento_global(b"ACGT", b"AGT", puntaje_adn, 1);
    assert_eq!(adn.puntaje, 2);
    assert_eq!(adn.fila_v, b"ACGT");
    assert_eq!(adn.fila_w, b"A-GT");
    let proteina = alineamiento_global(b"WHAT", b"WHT", blosum62, 5);
    assert_eq!(proteina.puntaje, 19);
    assert_eq!(proteina.fila_w, b"WH-T");
    let vacia = alineamiento_global(b"", b"ACG", puntaje_adn, 2);
    assert_eq!(vacia.puntaje, -6);
    assert_eq!(vacia.fila_v, b"---");
    assert_eq!(vacia.fila_w, b"ACG");
    let nada = alineamiento_global(b"", b"", blosum62, 5);
    assert_eq!(nada.puntaje, 0);
    assert!(nada.fila_v.is_empty());
}

#[test]
fn alineamiento_local_a_mano() {
    let local = alineamiento_local(b"TTTACGTAAA", b"GGACGTCC", puntaje_adn, 2);
    assert_eq!(local.puntaje, 4);
    assert_eq!(local.fila_v, b"ACGT");
    assert_eq!(local.fila_w, b"ACGT");
    assert_eq!((local.inicio_v, local.inicio_w), (3, 2));
    let nada = alineamiento_local(b"AAA", b"CCC", puntaje_adn, 2);
    assert_eq!(nada.puntaje, 0);
    assert!(nada.fila_v.is_empty());
    let con_pam = alineamiento_local(b"DDDWWCWDDD", b"LLWWCWLL", pam250, 5);
    assert_eq!(con_pam.fila_v, b"WWCW");
    assert_eq!(con_pam.puntaje, 17 + 17 + 12 + 17);
}

#[test]
fn distancia_de_edicion_a_mano() {
    assert_eq!(distancia_de_edicion(b"ACGT", b"ACGT"), 0);
    assert_eq!(distancia_de_edicion(b"ACGT", b"AGT"), 1);
    assert_eq!(distancia_de_edicion(b"AAAA", b"TTTT"), 4);
    assert_eq!(distancia_de_edicion(b"ACGTACGT", b"CGTACGTA"), 2);
    assert_eq!(distancia_de_edicion(b"", b"ACG"), 3);
    assert_eq!(distancia_de_edicion(b"ACG", b""), 3);
}

#[test]
fn ajustado_y_solapado_a_mano() {
    let exacto = alineamiento_ajustado(b"TTTTACGTTTTT", b"ACG", puntaje_adn, 1);
    assert_eq!(exacto.puntaje, 3);
    assert_eq!(exacto.inicio_v, 4);
    assert_eq!(exacto.fila_v, b"ACG");
    let con_hueco = alineamiento_ajustado(b"TTTTACGTTTT", b"ACCG", puntaje_adn, 1);
    assert_eq!(con_hueco.puntaje, 2);
    assert_eq!(sin_huecos(&con_hueco.fila_w), b"ACCG");
    let solape = alineamiento_solapado(b"GGGGACGT", b"ACGTCCCC", puntaje_adn, 2);
    assert_eq!(solape.puntaje, 4);
    assert_eq!(solape.fila_v, b"ACGT");
    assert_eq!(solape.fila_w, b"ACGT");
    assert_eq!((solape.inicio_v, solape.inicio_w), (4, 0));
    let sin_solape = alineamiento_solapado(b"AAAA", b"CCCC", puntaje_adn, 2);
    assert_eq!(sin_solape.puntaje, 0);
}

#[test]
fn huecos_afines_prefieren_un_solo_hueco_largo() {
    let adn = alineamiento_afin(b"AAAACCCC", b"AAAATTTTCCCC", puntaje_adn, 3, 1);
    assert_eq!(adn.puntaje, 2);
    assert_eq!(adn.fila_v, b"AAAA----CCCC");
    let proteina = alineamiento_afin(b"MKVLAAGIV", b"MKVLGIV", blosum62, 11, 1);
    assert_eq!(proteina.puntaje, 20);
    assert_eq!(proteina.fila_w, b"MKVL--GIV");
    assert_eq!(
        puntuar_alineamiento_afin(&proteina.fila_v, &proteina.fila_w, blosum62, 11, 1),
        20
    );
    let vacio = alineamiento_afin(b"", b"ACG", puntaje_adn, 3, 1);
    assert_eq!(vacio.puntaje, -5);
    assert_eq!(vacio.fila_v, b"---");
    let nada = alineamiento_afin(b"", b"", puntaje_adn, 3, 1);
    assert_eq!(nada.puntaje, 0);
}

#[test]
fn arista_central_a_mano() {
    let arista = arista_central(b"ACGT", b"ACGT", puntaje_adn, 1).unwrap();
    assert_eq!(arista.desde, (2, 2));
    assert_eq!(arista.hasta, (3, 3));
    assert_eq!(arista.longitudes.len(), 5);
    assert_eq!(arista.longitudes[2], 4);
    assert!(arista_central(b"", b"", puntaje_adn, 1).is_none());
    let solo_v = arista_central(b"AC", b"", puntaje_adn, 1).unwrap();
    assert_eq!(solo_v.desde.1, 0);
    assert_eq!(solo_v.hasta.1, 0);
    assert_eq!(solo_v.hasta.0, solo_v.desde.0 + 1);
    let solo_w = arista_central(b"", b"ACG", puntaje_adn, 1).unwrap();
    assert_eq!(solo_w.desde, (0, 1));
    assert_eq!(solo_w.hasta, (0, 2));
}

#[test]
fn espacio_lineal_da_lo_mismo_que_la_tabla_completa() {
    let v = b"MSDLAREITPVNIEEELKSSYLDYAMSVIVGRALPDVRDGLKPVHRRVLYAM";
    let w = b"MSDMAERLALHEFTENAYLNYSMYVIMDRALPFIGDGLKPVQRRIVYAM";
    let lineal = alineamiento_espacio_lineal(v, w, blosum62, 5);
    let completo = alineamiento_global(v, w, blosum62, 5);
    assert_eq!(lineal.puntaje, completo.puntaje);
    revisar(&lineal, v, w, blosum62, 5);
    let vacio = alineamiento_espacio_lineal(b"", b"AC", puntaje_adn, 1);
    assert_eq!(vacio.puntaje, -2);
    assert_eq!(vacio.fila_v, b"--");
}

#[test]
fn tres_secuencias_a_mano() {
    let iguales = subsecuencia_comun_de_tres(b"ACGT", b"ACGT", b"ACGT");
    assert_eq!(iguales.puntaje, 4);
    assert_eq!(iguales.fila_u, b"ACGT");
    let distintas = subsecuencia_comun_de_tres(b"AAA", b"CCC", b"GGG");
    assert_eq!(distintas.puntaje, 0);
    let mezcla = subsecuencia_comun_de_tres(b"ACG", b"AG", b"CAG");
    assert_eq!(mezcla.puntaje, 2);
    assert_eq!(sin_huecos(&mezcla.fila_u), b"ACG");
    assert_eq!(sin_huecos(&mezcla.fila_v), b"AG");
    assert_eq!(sin_huecos(&mezcla.fila_w), b"CAG");
    let vacia = subsecuencia_comun_de_tres(b"", b"AC", b"A");
    assert_eq!(vacia.puntaje, 0);
    assert_eq!(vacia.fila_u, b"--");
}

#[test]
fn puntajes_de_columnas() {
    assert_eq!(puntaje_adn(b'A', b'A'), 1);
    assert_eq!(puntaje_adn(b'A', b'C'), -1);
    assert_eq!(puntaje_edicion(b'G', b'G'), 0);
    assert_eq!(puntaje_edicion(b'G', b'T'), -1);
    assert_eq!(
        puntuar_alineamiento(b"AC-T", b"A-GT", puntaje_adn, 2),
        1 - 2 - 2 + 1
    );
    assert_eq!(
        puntuar_alineamiento_afin(b"A---CT", b"AGGGC-", puntaje_adn, 3, 1),
        1 - 5 + 1 - 3
    );
}

#[test]
fn todas_las_variantes_coinciden_con_versiones_simples() {
    let mut azar = Azar::nuevo(2024);
    for vuelta in 0..300 {
        let letras: &[u8] = if vuelta % 2 == 0 {
            b"ACGT"
        } else {
            b"ARNDCQEGHILKMFPSTWYV"
        };
        let puntuar: Puntuacion = match vuelta % 3 {
            0 => puntaje_adn,
            1 => blosum62,
            _ => pam250,
        };
        let indel = 1 + azar.entero(6) as i64;
        let v = al_azar(&mut azar, 7, letras);
        let w = al_azar(&mut azar, 7, letras);

        let global = alineamiento_global(&v, &w, puntuar, indel);
        assert_eq!(global.puntaje, global_simple(&v, &w, puntuar, indel));
        assert_eq!(revisar(&global, &v, &w, puntuar, indel), (v.len(), w.len()));
        assert_eq!((global.inicio_v, global.inicio_w), (0, 0));

        let mut mejor_local = 0;
        let mut mejor_ajustado = i64::MIN;
        let mut mejor_solapado = i64::MIN;
        for a in 0..=v.len() {
            for b in a..=v.len() {
                mejor_ajustado = mejor_ajustado.max(global_simple(&v[a..b], &w, puntuar, indel));
                for c in 0..=w.len() {
                    for d in c..=w.len() {
                        mejor_local =
                            mejor_local.max(global_simple(&v[a..b], &w[c..d], puntuar, indel));
                    }
                }
            }
            for d in 0..=w.len() {
                mejor_solapado =
                    mejor_solapado.max(global_simple(&v[a..], &w[..d], puntuar, indel));
            }
        }
        let local = alineamiento_local(&v, &w, puntuar, indel);
        assert_eq!(local.puntaje, mejor_local);
        revisar(&local, &v, &w, puntuar, indel);

        let ajustado = alineamiento_ajustado(&v, &w, puntuar, indel);
        assert_eq!(ajustado.puntaje, mejor_ajustado);
        assert_eq!(ajustado.inicio_w, 0);
        assert_eq!(revisar(&ajustado, &v, &w, puntuar, indel).1, w.len());

        let solapado = alineamiento_solapado(&v, &w, puntuar, indel);
        assert_eq!(solapado.puntaje, mejor_solapado);
        assert_eq!(solapado.inicio_w, 0);
        assert_eq!(revisar(&solapado, &v, &w, puntuar, indel).0, v.len());

        let lineal = alineamiento_espacio_lineal(&v, &w, puntuar, indel);
        assert_eq!(lineal.puntaje, global.puntaje);
        revisar(&lineal, &v, &w, puntuar, indel);

        if let Some(arista) = arista_central(&v, &w, puntuar, indel) {
            let (i, j) = arista.desde;
            let (siguiente_i, siguiente_j) = arista.hasta;
            assert_eq!(j, w.len() / 2);
            let paso = if siguiente_i > i && siguiente_j > j {
                puntuar(v[i], w[j]) as i64
            } else {
                assert!(siguiente_i + siguiente_j == i + j + 1);
                -indel
            };
            let antes = global_simple(&v[..i], &w[..j], puntuar, indel);
            let despues = global_simple(&v[siguiente_i..], &w[siguiente_j..], puntuar, indel);
            assert_eq!(antes + paso + despues, global.puntaje);
            assert_eq!(arista.longitudes.iter().max(), Some(&global.puntaje));
        } else {
            assert!(v.is_empty() && w.is_empty());
        }

        let apertura = indel + azar.entero(5) as i64;
        let extension = 1 + azar.entero(indel as usize) as i64;
        let afin = alineamiento_afin(&v, &w, puntuar, apertura, extension);
        assert_eq!(
            afin.puntaje,
            afin_simple(&v, &w, puntuar, apertura, extension)
        );
        assert_eq!(sin_huecos(&afin.fila_v), v);
        assert_eq!(sin_huecos(&afin.fila_w), w);
        assert_eq!(
            puntuar_alineamiento_afin(&afin.fila_v, &afin.fila_w, puntuar, apertura, extension),
            afin.puntaje
        );

        let distancia = distancia_de_edicion(&v, &w);
        assert_eq!(distancia as i64, -global_simple(&v, &w, puntaje_edicion, 1));

        let comun = subsecuencia_comun_mas_larga(&v, &w);
        assert_eq!(comun.len(), comun_simple(&v, &w));
        assert!(es_subsecuencia(&comun, &v) && es_subsecuencia(&comun, &w));

        let u = al_azar(&mut azar, 6, b"ACG");
        let triple = subsecuencia_comun_de_tres(&u, &v, &w);
        assert_eq!(triple.puntaje, comun_de_tres_simple(&u, &v, &w));
        assert_eq!(sin_huecos(&triple.fila_u), u);
        assert_eq!(sin_huecos(&triple.fila_v), v);
        assert_eq!(sin_huecos(&triple.fila_w), w);
        let columnas = triple.fila_u.len();
        assert!(triple.fila_v.len() == columnas && triple.fila_w.len() == columnas);
        let iguales = (0..columnas)
            .filter(|&k| {
                triple.fila_u[k] != b'-'
                    && triple.fila_u[k] == triple.fila_v[k]
                    && triple.fila_v[k] == triple.fila_w[k]
            })
            .count();
        assert_eq!(iguales, triple.puntaje);
    }
}

#[test]
fn cambio_coincide_con_busqueda_exhaustiva() {
    let mut azar = Azar::nuevo(11);
    for _ in 0..100 {
        let mut monedas: Vec<usize> = (0..1 + azar.entero(4))
            .map(|_| 1 + azar.entero(9))
            .collect();
        monedas.dedup();
        let dinero = azar.entero(30);
        let mut menos = vec![usize::MAX; dinero + 1];
        menos[0] = 0;
        for cantidad in 1..=dinero {
            for &moneda in &monedas {
                if moneda <= cantidad && menos[cantidad - moneda] != usize::MAX {
                    menos[cantidad] = menos[cantidad].min(menos[cantidad - moneda] + 1);
                }
            }
        }
        match cambio_minimo(dinero, &monedas) {
            Some(usadas) => {
                assert_eq!(usadas.len(), menos[dinero]);
                assert_eq!(usadas.iter().sum::<usize>(), dinero);
                assert!(usadas.iter().all(|moneda| monedas.contains(moneda)));
            }
            None => assert_eq!(menos[dinero], usize::MAX),
        }
    }
}

#[test]
fn afin_con_extension_mayor_que_apertura_da_un_puntaje_que_coincide_con_sus_filas() {
    let mut azar = Azar::nuevo(77);
    for _ in 0..300 {
        let v = al_azar(&mut azar, 6, b"ACGT");
        let w = al_azar(&mut azar, 6, b"ACGT");
        let apertura = azar.entero(5) as i64;
        let extension = apertura + 1 + azar.entero(4) as i64;
        let resultado = alineamiento_afin(&v, &w, puntaje_adn, apertura, extension);
        let reevaluado = puntuar_alineamiento_afin(
            &resultado.fila_v,
            &resultado.fila_w,
            puntaje_adn,
            apertura,
            extension,
        );
        assert_eq!(
            resultado.puntaje, reevaluado,
            "{v:?} {w:?} {apertura} {extension}"
        );
        assert_eq!(sin_huecos(&resultado.fila_v), v);
        assert_eq!(sin_huecos(&resultado.fila_w), w);
    }
}
