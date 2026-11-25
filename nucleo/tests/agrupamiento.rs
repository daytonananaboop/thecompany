use nucleo::agrupamiento::*;
use nucleo::azar::Azar;

fn cerca(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

fn puntos(filas: &[&[f64]]) -> Vec<Vec<f64>> {
    filas.iter().map(|fila| fila.to_vec()).collect()
}

#[test]
fn distancia_euclidea_y_al_cuadrado() {
    assert_eq!(distancia_euclidea(&[0.0, 0.0], &[3.0, 4.0]), 5.0);
    assert_eq!(distancia_al_cuadrado(&[0.0, 0.0], &[3.0, 4.0]), 25.0);
    assert_eq!(
        distancia_al_cuadrado(&[1.0, 2.0, 3.0], &[1.0, 2.0, 3.0]),
        0.0
    );
    assert_eq!(distancia_euclidea(&[], &[]), 0.0);
}

#[test]
fn correlacion_de_pearson_calculada_a_mano() {
    assert!(cerca(
        correlacion_pearson(&[1.0, 2.0, 3.0, 4.0], &[1.0, 3.0, 2.0, 4.0]),
        0.8
    ));
    assert!(cerca(
        correlacion_pearson(&[1.0, 2.0, 3.0], &[2.0, 4.0, 6.0]),
        1.0
    ));
    assert!(cerca(
        correlacion_pearson(&[1.0, 2.0, 3.0], &[30.0, 20.0, 10.0]),
        -1.0
    ));
    assert_eq!(correlacion_pearson(&[5.0, 5.0, 5.0], &[1.0, 2.0, 3.0]), 0.0);
    assert_eq!(correlacion_pearson(&[], &[]), 0.0);
    assert!(cerca(
        distancia_pearson(&[1.0, 2.0, 3.0], &[3.0, 2.0, 1.0]),
        2.0
    ));
    assert!(cerca(
        distancia_pearson(&[1.0, 2.0, 3.0], &[10.0, 20.0, 30.0]),
        0.0
    ));
}

#[test]
fn pearson_mira_la_forma_y_euclidea_la_magnitud() {
    let pequeno = [1.0, 2.0, 3.0];
    let grande = [10.0, 20.0, 30.0];
    let plano = [1.5, 1.5, 1.5];
    assert!(distancia_euclidea(&pequeno, &plano) < distancia_euclidea(&pequeno, &grande));
    assert!(distancia_pearson(&pequeno, &grande) < distancia_pearson(&pequeno, &plano));
}

#[test]
fn matrices_de_distancias_y_correlaciones() {
    let datos = puntos(&[&[0.0, 0.0], &[3.0, 4.0], &[6.0, 8.0]]);
    assert_eq!(
        matriz_de_distancias(&datos, Medida::Euclidea),
        vec![
            vec![0.0, 5.0, 10.0],
            vec![5.0, 0.0, 5.0],
            vec![10.0, 5.0, 0.0]
        ]
    );
    let perfiles = puntos(&[
        &[1.0, 2.0, 3.0],
        &[10.0, 20.0, 30.0],
        &[3.0, 2.0, 1.0],
        &[2.0, 2.0, 2.0],
    ]);
    let pearson = matriz_de_distancias(&perfiles, Medida::Pearson);
    assert!(cerca(pearson[0][1], 0.0));
    assert!(cerca(pearson[0][2], 2.0));
    assert!(cerca(pearson[0][3], 1.0));
    for (indice, fila) in pearson.iter().enumerate() {
        assert_eq!(fila[indice], 0.0);
    }
    let similitudes = matriz_de_correlaciones(&perfiles);
    assert!(cerca(similitudes[1][2], -1.0));
    assert!(cerca(similitudes[2][1], -1.0));
    assert!(cerca(similitudes[0][0], 1.0));
}

#[test]
fn las_matrices_son_simetricas_y_coinciden_con_cada_par() {
    let mut azar = Azar::nuevo(5);
    let datos: Vec<Vec<f64>> = (0..23)
        .map(|_| (0..6).map(|_| azar.normal()).collect())
        .collect();
    let euclidea = matriz_de_distancias(&datos, Medida::Euclidea);
    let pearson = matriz_de_distancias(&datos, Medida::Pearson);
    let correlaciones = matriz_de_correlaciones(&datos);
    for fila in 0..datos.len() {
        assert_eq!(euclidea[fila].len(), datos.len());
        assert_eq!(euclidea[fila][fila], 0.0);
        assert_eq!(pearson[fila][fila], 0.0);
        for columna in 0..datos.len() {
            assert_eq!(euclidea[fila][columna], euclidea[columna][fila]);
            assert_eq!(pearson[fila][columna], pearson[columna][fila]);
            assert_eq!(correlaciones[fila][columna], correlaciones[columna][fila]);
            if fila != columna {
                assert!(cerca(
                    euclidea[fila][columna],
                    distancia_euclidea(&datos[fila], &datos[columna])
                ));
                assert!(cerca(
                    pearson[fila][columna],
                    distancia_pearson(&datos[fila], &datos[columna])
                ));
            }
        }
    }
    assert!(matriz_de_distancias(&[], Medida::Euclidea).is_empty());
    assert!(matriz_de_correlaciones(&[]).is_empty());
}

#[test]
fn revisar_puntos_rechaza_datos_rotos() {
    assert!(revisar_puntos(&[]).is_err());
    assert!(revisar_puntos(&[vec![]]).is_err());
    assert!(revisar_puntos(&puntos(&[&[1.0, 2.0], &[3.0]])).is_err());
    assert!(revisar_puntos(&puntos(&[&[1.0, f64::NAN]])).is_err());
    assert!(revisar_puntos(&puntos(&[&[1.0, f64::INFINITY]])).is_err());
    assert_eq!(revisar_puntos(&puntos(&[&[1.0, 2.0], &[3.0, 4.0]])), Ok(2));
}

#[test]
fn el_mas_lejano_primero_elige_centros_alejados() {
    let datos = puntos(&[
        &[0.0, 0.0],
        &[1.0, 0.0],
        &[10.0, 0.0],
        &[10.0, 1.0],
        &[5.0, 5.0],
    ]);
    let centros = mas_lejano_primero(&datos, 3).unwrap();
    assert_eq!(centros, vec![0, 3, 4]);
    let elegidos: Vec<Vec<f64>> = centros.iter().map(|&i| datos[i].clone()).collect();
    assert!(cerca(radio_de_los_centros(&datos, &elegidos), 1.0));
    assert_eq!(mas_lejano_primero(&datos, 1).unwrap(), vec![0]);
}

#[test]
fn el_mas_lejano_primero_desempata_y_no_repite() {
    let datos = puntos(&[&[0.0], &[2.0], &[-2.0]]);
    assert_eq!(mas_lejano_primero(&datos, 2).unwrap(), vec![0, 1]);
    let iguales = puntos(&[&[1.0, 1.0], &[1.0, 1.0], &[1.0, 1.0]]);
    assert_eq!(mas_lejano_primero(&iguales, 3).unwrap(), vec![0, 1, 2]);
    assert!(mas_lejano_primero(&datos, 0).is_err());
    assert!(mas_lejano_primero(&datos, 4).is_err());
    assert!(mas_lejano_primero(&[], 1).is_err());
}

#[test]
fn distorsion_es_el_promedio_de_distancias_al_cuadrado() {
    let datos = puntos(&[&[0.0, 0.0], &[2.0, 0.0], &[10.0, 0.0]]);
    let centros = puntos(&[&[1.0, 0.0], &[10.0, 0.0]]);
    assert!(cerca(distorsion_cuadratica(&datos, &centros), 2.0 / 3.0));
    assert_eq!(distorsion_cuadratica(&[], &centros), 0.0);
    assert_eq!(centro_mas_cercano(&[9.0, 0.0], &centros), 1);
    assert_eq!(centro_mas_cercano(&[5.5, 0.0], &centros), 0);
}

#[test]
fn lloyd_con_los_primeros_puntos_como_en_el_libro() {
    let datos = puntos(&[&[1.0], &[2.0], &[3.0], &[10.0], &[11.0], &[12.0]]);
    let resultado = k_medias(&datos, 2, Inicio::PrimerosPuntos, 100).unwrap();
    assert_eq!(resultado.centros, vec![vec![2.0], vec![11.0]]);
    assert_eq!(resultado.grupos, vec![0, 0, 0, 1, 1, 1]);
    assert_eq!(resultado.iteraciones, 2);
    assert!(resultado.convergio);
    assert!(cerca(resultado.distorsion, 4.0 / 6.0));
}

#[test]
fn lloyd_se_detiene_en_el_maximo_de_iteraciones() {
    let datos = puntos(&[&[1.0], &[2.0], &[3.0], &[10.0], &[11.0], &[12.0]]);
    let resultado = k_medias(&datos, 2, Inicio::PrimerosPuntos, 1).unwrap();
    assert!(!resultado.convergio);
    assert_eq!(resultado.iteraciones, 1);
    assert!(cerca(resultado.centros[1][0], 7.6));
    assert_eq!(resultado.grupos, vec![0, 0, 0, 1, 1, 1]);
    assert!(cerca(resultado.distorsion, 41.68 / 6.0));
    let quieto = k_medias(&datos, 2, Inicio::PrimerosPuntos, 0).unwrap();
    assert_eq!(quieto.centros, vec![vec![1.0], vec![2.0]]);
    assert_eq!(quieto.grupos, vec![0, 1, 1, 1, 1, 1]);
}

#[test]
fn lloyd_sale_de_un_mal_comienzo() {
    let datos = puntos(&[&[0.0], &[0.1], &[10.0], &[10.1]]);
    let resultado = k_medias(&datos, 2, Inicio::PrimerosPuntos, 100).unwrap();
    assert_eq!(resultado.grupos, vec![0, 0, 1, 1]);
    assert!(cerca(resultado.centros[0][0], 0.05));
    assert!(cerca(resultado.centros[1][0], 10.05));
}

#[test]
fn lloyd_paso_a_paso_coincide_con_k_medias() {
    let mut azar = Azar::nuevo(7);
    let datos: Vec<Vec<f64>> = (0..60)
        .map(|_| vec![azar.normal() * 3.0, azar.normal() * 3.0, azar.normal()])
        .collect();
    let resultado = k_medias(&datos, 4, Inicio::PrimerosPuntos, 100).unwrap();
    let mut centros = centros_iniciales(&datos, 4, Inicio::PrimerosPuntos);
    let mut grupos = de_centros_a_grupos(&datos, &centros);
    let mut vueltas = 0;
    loop {
        vueltas += 1;
        centros = de_grupos_a_centros(&datos, &grupos, &centros);
        let nuevos = de_centros_a_grupos(&datos, &centros);
        if nuevos == grupos {
            break;
        }
        grupos = nuevos;
    }
    assert_eq!(resultado.centros, centros);
    assert_eq!(resultado.grupos, grupos);
    assert_eq!(resultado.iteraciones, vueltas);
}

#[test]
fn un_grupo_vacio_conserva_su_centro() {
    let datos = puntos(&[&[0.0], &[1.0]]);
    let centros = puntos(&[&[0.5], &[100.0]]);
    let grupos = de_centros_a_grupos(&datos, &centros);
    assert_eq!(grupos, vec![0, 0]);
    let movidos = de_grupos_a_centros(&datos, &grupos, &centros);
    assert_eq!(movidos, vec![vec![0.5], vec![100.0]]);
}

#[test]
fn k_medias_mas_mas_reparte_los_centros_entre_nubes_lejanas() {
    let mut datos = Vec::new();
    for nube in [0.0, 100.0, 200.0] {
        for paso in 0..5 {
            datos.push(vec![nube + paso as f64 * 0.1]);
        }
    }
    for semilla in 0..20 {
        let elegidos = k_medias_mas_mas(&datos, 3, semilla);
        let mut nubes: Vec<usize> = elegidos.iter().map(|&i| i / 5).collect();
        nubes.sort();
        assert_eq!(nubes, vec![0, 1, 2], "semilla {semilla}");
        let resultado = k_medias(&datos, 3, Inicio::KMediasMasMas(semilla), 100).unwrap();
        for nube in 0..3 {
            let grupo = resultado.grupos[nube * 5];
            assert!(
                resultado.grupos[nube * 5..nube * 5 + 5]
                    .iter()
                    .all(|&g| g == grupo)
            );
        }
    }
    assert_eq!(
        k_medias_mas_mas(&datos, 3, 9),
        k_medias_mas_mas(&datos, 3, 9)
    );
    assert!(k_medias_mas_mas(&[], 3, 1).is_empty());
}

#[test]
fn k_medias_rechaza_entradas_imposibles() {
    let datos = puntos(&[&[1.0], &[2.0]]);
    assert!(k_medias(&datos, 0, Inicio::PrimerosPuntos, 10).is_err());
    assert!(k_medias(&datos, 3, Inicio::PrimerosPuntos, 10).is_err());
    assert!(k_medias(&[], 1, Inicio::PrimerosPuntos, 10).is_err());
    let rotos = puntos(&[&[1.0, 2.0], &[1.0]]);
    assert!(k_medias(&rotos, 1, Inicio::PrimerosPuntos, 10).is_err());
}

#[test]
fn k_medias_suaves_una_vuelta_calculada_a_mano() {
    let datos = puntos(&[&[0.0], &[2.0]]);
    let resultado = k_medias_suaves(&datos, 2, 1.0, 1).unwrap();
    let cerca_del_propio = 1.0 / (1.0 + (-2.0f64).exp());
    let esperado = 2.0 * (1.0 - cerca_del_propio);
    assert!(cerca(resultado.centros[0][0], esperado));
    assert!(cerca(resultado.centros[1][0], 2.0 - esperado));
    for fila in &resultado.responsabilidades {
        assert!(cerca(fila.iter().sum::<f64>(), 1.0));
    }
    let separacion = 2.0 - 2.0 * esperado;
    let propia = 1.0 / (1.0 + (-separacion).exp());
    assert!(cerca(resultado.responsabilidades[0][0], propia));
    assert!(cerca(resultado.responsabilidades[1][1], propia));
}

#[test]
fn k_medias_suaves_con_rigidez_extrema_no_da_nan() {
    let datos = puntos(&[&[0.0], &[1000.0], &[2000.0], &[1.0]]);
    let dura = k_medias_suaves(&datos, 2, 1e6, 50).unwrap();
    for fila in &dura.responsabilidades {
        assert!(fila.iter().all(|valor| valor.is_finite()));
        assert!(fila.iter().any(|&valor| valor > 0.999));
    }
    let blanda = k_medias_suaves(&datos, 2, 1e-12, 5).unwrap();
    for fila in &blanda.responsabilidades {
        assert!(fila.iter().all(|&valor| (valor - 0.5).abs() < 1e-6));
    }
    assert!(k_medias_suaves(&datos, 2, -1.0, 5).is_err());
    assert!(k_medias_suaves(&datos, 2, f64::NAN, 5).is_err());
    assert!(k_medias_suaves(&datos, 5, 1.0, 5).is_err());
    assert!(k_medias_suaves(&datos, 0, 1.0, 5).is_err());
}

#[test]
fn k_medias_suaves_con_cero_vueltas_deja_los_primeros_puntos() {
    let datos = puntos(&[&[0.0, 1.0], &[4.0, 4.0], &[9.0, 9.0]]);
    let resultado = k_medias_suaves(&datos, 2, 1.0, 0).unwrap();
    assert_eq!(resultado.centros, vec![vec![0.0, 1.0], vec![4.0, 4.0]]);
    assert_eq!(resultado.responsabilidades.len(), 3);
}

fn matriz(filas: &[&[f64]]) -> Vec<Vec<f64>> {
    puntos(filas)
}

#[test]
fn jerarquico_con_distancia_promedio_calculado_a_mano() {
    let distancias = matriz(&[
        &[0.0, 2.0, 6.0, 10.0],
        &[2.0, 0.0, 5.0, 9.0],
        &[6.0, 5.0, 0.0, 4.0],
        &[10.0, 9.0, 4.0, 0.0],
    ]);
    let uniones = agrupamiento_jerarquico(&distancias).unwrap();
    assert_eq!(uniones.len(), 3);
    assert_eq!((uniones[0].primero, uniones[0].segundo), (0, 1));
    assert_eq!(uniones[0].miembros, vec![0, 1]);
    assert!(cerca(uniones[0].distancia, 2.0));
    assert_eq!((uniones[1].primero, uniones[1].segundo), (2, 3));
    assert_eq!(uniones[1].miembros, vec![2, 3]);
    assert!(cerca(uniones[1].distancia, 4.0));
    assert_eq!((uniones[2].primero, uniones[2].segundo), (4, 5));
    assert_eq!(uniones[2].miembros, vec![0, 1, 2, 3]);
    assert!(cerca(uniones[2].distancia, 7.5));
    let nombres: Vec<String> = ["a", "b", "c", "d"].iter().map(|n| n.to_string()).collect();
    assert_eq!(
        arbol_newick(&uniones, &nombres),
        "((a:1,b:1):2.75,(c:2,d:2):1.75);"
    );
}

#[test]
fn jerarquico_bordes() {
    let uno = agrupamiento_jerarquico(&matriz(&[&[0.0]])).unwrap();
    assert!(uno.is_empty());
    assert_eq!(arbol_newick(&uno, &["solo".to_string()]), "solo;");
    assert!(agrupamiento_jerarquico(&[]).is_err());
    assert!(agrupamiento_jerarquico(&matriz(&[&[0.0, 1.0], &[1.0]])).is_err());
    assert!(agrupamiento_jerarquico(&matriz(&[&[0.0, f64::NAN], &[1.0, 0.0]])).is_err());
    let empates = matriz(&[&[0.0, 1.0, 1.0], &[1.0, 0.0, 1.0], &[1.0, 1.0, 0.0]]);
    let uniones = agrupamiento_jerarquico(&empates).unwrap();
    assert_eq!(uniones[0].miembros, vec![0, 1]);
    assert_eq!(uniones[1].miembros, vec![0, 1, 2]);
    let raros = vec!["a b".to_string(), "c(d)".to_string(), "e:f".to_string()];
    let arbol = arbol_newick(&uniones, &raros);
    assert_eq!(arbol, "((a_b:0.5,c_d_:0.5):0,e_f:0.5);");
}

fn jerarquico_ingenuo(distancias: &[Vec<f64>]) -> Vec<(Vec<usize>, f64)> {
    let mut grupos: Vec<Vec<usize>> = (0..distancias.len()).map(|i| vec![i]).collect();
    let mut resultado = Vec::new();
    while grupos.len() > 1 {
        let mut mejor = (0, 1, f64::INFINITY);
        for a in 0..grupos.len() {
            for b in a + 1..grupos.len() {
                let mut suma = 0.0;
                for &x in &grupos[a] {
                    for &y in &grupos[b] {
                        suma += distancias[x][y];
                    }
                }
                let promedio = suma / (grupos[a].len() * grupos[b].len()) as f64;
                if promedio < mejor.2 {
                    mejor = (a, b, promedio);
                }
            }
        }
        let segundo = grupos.remove(mejor.1);
        let primero = grupos.remove(mejor.0);
        let mut unido = primero;
        unido.extend(segundo);
        let mut ordenado = unido.clone();
        ordenado.sort();
        resultado.push((ordenado, mejor.2));
        grupos.push(unido);
    }
    resultado
}

#[test]
fn jerarquico_coincide_con_la_version_lenta() {
    let mut azar = Azar::nuevo(42);
    for tamano in [2, 3, 5, 9, 17, 30] {
        let valores: Vec<f64> = (0..tamano * tamano)
            .map(|_| azar.decimal() * 10.0)
            .collect();
        let distancias: Vec<Vec<f64>> = (0..tamano)
            .map(|i| {
                (0..tamano)
                    .map(|j| match i.cmp(&j) {
                        std::cmp::Ordering::Less => valores[i * tamano + j],
                        std::cmp::Ordering::Greater => valores[j * tamano + i],
                        std::cmp::Ordering::Equal => 0.0,
                    })
                    .collect()
            })
            .collect();
        let rapido = agrupamiento_jerarquico(&distancias).unwrap();
        let lento = jerarquico_ingenuo(&distancias);
        assert_eq!(rapido.len(), lento.len());
        for (union, (miembros, distancia)) in rapido.iter().zip(&lento) {
            let mut ordenados = union.miembros.clone();
            ordenados.sort();
            assert_eq!(&ordenados, miembros);
            assert!((union.distancia - distancia).abs() < 1e-9);
        }
    }
}

#[test]
fn cast_con_bloques_calculado_a_mano() {
    let similitudes = matriz(&[
        &[1.0, 0.9, 0.8, 0.1, 0.0],
        &[0.9, 1.0, 0.85, 0.2, 0.1],
        &[0.8, 0.85, 1.0, 0.0, 0.2],
        &[0.1, 0.2, 0.0, 1.0, 0.9],
        &[0.0, 0.1, 0.2, 0.9, 1.0],
    ]);
    assert_eq!(
        cast(&similitudes, 0.5).unwrap(),
        vec![vec![0, 1, 2], vec![3, 4]]
    );
    assert_eq!(
        cast(&similitudes, 0.95).unwrap(),
        vec![vec![0], vec![1], vec![2], vec![3], vec![4]]
    );
    assert_eq!(cast(&similitudes, 0.0).unwrap().len(), 1);
}

#[test]
fn cast_empieza_por_el_de_mas_vecinos() {
    let similitudes = matriz(&[
        &[1.0, 0.1, 0.1, 0.1],
        &[0.1, 1.0, 0.9, 0.9],
        &[0.1, 0.9, 1.0, 0.9],
        &[0.1, 0.9, 0.9, 1.0],
    ]);
    assert_eq!(
        cast(&similitudes, 0.6).unwrap(),
        vec![vec![1, 2, 3], vec![0]]
    );
}

#[test]
fn cast_quita_al_que_queda_lejos_del_grupo() {
    let similitudes = matriz(&[
        &[1.0, 0.95, 0.9, 0.9, 0.9],
        &[0.95, 1.0, 0.33, 0.33, 0.33],
        &[0.9, 0.33, 1.0, 0.9, 0.9],
        &[0.9, 0.33, 0.9, 1.0, 0.9],
        &[0.9, 0.33, 0.9, 0.9, 1.0],
    ]);
    assert_eq!(
        cast(&similitudes, 0.6).unwrap(),
        vec![vec![0, 2, 3, 4], vec![1]]
    );
    let ciclo = matriz(&[
        &[1.0, 0.7, 0.7, -1.0],
        &[0.7, 1.0, 0.7, 0.6],
        &[0.7, 0.7, 1.0, 0.6],
        &[-1.0, 0.6, 0.6, 1.0],
    ]);
    let grupos = cast(&ciclo, 0.6).unwrap();
    let mut todos: Vec<usize> = grupos.iter().flatten().copied().collect();
    todos.sort();
    assert_eq!(todos, vec![0, 1, 2, 3]);
}

#[test]
fn cast_rechaza_matrices_rotas() {
    assert!(cast(&[], 0.5).is_err());
    assert!(cast(&matriz(&[&[1.0, 0.5], &[0.5]]), 0.5).is_err());
    assert!(cast(&matriz(&[&[1.0]]), f64::NAN).is_err());
    let raro = cast(&matriz(&[&[0.0, 0.0], &[0.0, 0.0]]), 2.0).unwrap();
    assert_eq!(raro, vec![vec![0], vec![1]]);
}

#[test]
fn cast_reparte_todos_los_genes_una_sola_vez() {
    let mut azar = Azar::nuevo(3);
    let mut perfiles = Vec::new();
    for patron in 0..3 {
        for _ in 0..8 {
            let perfil: Vec<f64> = (0..6)
                .map(|t| {
                    let base = match patron {
                        0 => t as f64,
                        1 => -(t as f64),
                        _ => -((t as f64 - 2.5).powi(2)),
                    };
                    base + azar.normal() * 0.2
                })
                .collect();
            perfiles.push(perfil);
        }
    }
    let similitudes = matriz_de_correlaciones(&perfiles);
    let grupos = cast(&similitudes, 0.7).unwrap();
    let mut todos: Vec<usize> = grupos.iter().flatten().copied().collect();
    todos.sort();
    assert_eq!(todos, (0..24).collect::<Vec<_>>());
    assert_eq!(grupos.len(), 3);
    for grupo in &grupos {
        let patron = grupo[0] / 8;
        assert!(grupo.iter().all(|&elemento| elemento / 8 == patron));
    }
}

#[test]
fn componentes_principales_de_una_recta() {
    let datos = puntos(&[&[0.0, 0.0], &[1.0, 2.0], &[2.0, 4.0], &[3.0, 6.0]]);
    let componentes = componentes_principales(&datos).unwrap();
    let raiz = 5.0f64.sqrt();
    assert!(cerca(componentes.primera[0], 1.0 / raiz));
    assert!(cerca(componentes.primera[1], 2.0 / raiz));
    assert!(cerca(componentes.varianza_primera, 25.0 / 3.0));
    assert!(componentes.varianza_segunda.abs() < 1e-9);
    assert!(cerca(componentes.varianza_total, 25.0 / 3.0));
    assert!(cerca(componentes.coordenadas[0][0], -1.5 * raiz));
    assert!(cerca(componentes.coordenadas[3][0], 1.5 * raiz));
    assert!(
        componentes
            .coordenadas
            .iter()
            .all(|punto| punto[1].abs() < 1e-6)
    );
    let centro = proyectar(&[1.5, 3.0], &componentes);
    assert!(cerca(centro[0], 0.0) && cerca(centro[1], 0.0));
}

#[test]
fn componentes_principales_de_una_cruz() {
    let datos = puntos(&[&[-1.0, 0.0], &[1.0, 0.0], &[0.0, -0.5], &[0.0, 0.5]]);
    let componentes = componentes_principales(&datos).unwrap();
    assert!(cerca(componentes.primera[0].abs(), 1.0));
    assert!(cerca(componentes.segunda[1].abs(), 1.0));
    assert!(cerca(componentes.varianza_primera, 2.0 / 3.0));
    assert!(cerca(componentes.varianza_segunda, 0.5 / 3.0));
    assert!(componentes.primera[0] > 0.0 && componentes.segunda[1] > 0.0);
}

#[test]
fn componentes_principales_bordes() {
    let uno = componentes_principales(&puntos(&[&[3.0, 4.0, 5.0]])).unwrap();
    assert_eq!(uno.coordenadas, vec![[0.0, 0.0]]);
    let recta = componentes_principales(&puntos(&[&[1.0], &[3.0], &[5.0]])).unwrap();
    assert!(cerca(recta.coordenadas[0][0], -2.0));
    assert_eq!(recta.coordenadas[2][1], 0.0);
    assert!(componentes_principales(&[]).is_err());
    let iguales = componentes_principales(&puntos(&[&[1.0, 1.0], &[1.0, 1.0]])).unwrap();
    assert!(
        iguales
            .coordenadas
            .iter()
            .all(|punto| punto[0] == 0.0 && punto[1] == 0.0)
    );
}

#[test]
fn codo_de_tres_nubes() {
    let datos = puntos(&[
        &[0.0],
        &[0.1],
        &[0.2],
        &[10.0],
        &[10.1],
        &[10.2],
        &[20.0],
        &[20.1],
        &[20.2],
    ]);
    let curva = curva_del_codo(&datos, 5, 1).unwrap();
    assert_eq!(curva.len(), 5);
    assert!(curva.windows(2).all(|par| par[1] <= par[0] + 1e-12));
    assert_eq!(k_del_codo(&curva), Some(3));
    assert_eq!(curva_del_codo(&datos, 50, 1).unwrap().len(), 9);
    assert!(curva_del_codo(&datos, 0, 1).is_err());
}

#[test]
fn k_del_codo_a_mano() {
    assert_eq!(k_del_codo(&[100.0, 40.0, 5.0, 4.0, 3.5]), Some(3));
    assert_eq!(k_del_codo(&[100.0, 20.0, 5.0, 4.0]), Some(3));
    assert_eq!(k_del_codo(&[100.0, 10.0, 5.0, 2.5]), Some(2));
    assert_eq!(k_del_codo(&[9.0, 4.0, 1.0, 0.9, 0.8, 0.0]), Some(3));
    assert_eq!(k_del_codo(&[10.0, 2.0, 3.0, 1.0]), Some(2));
    assert_eq!(k_del_codo(&[10.0, 9.0, 8.0, 7.0]), None);
    assert_eq!(k_del_codo(&[10.0, 5.0, 0.0]), None);
    assert_eq!(k_del_codo(&[10.0, 9.9, 0.0]), None);
    assert_eq!(k_del_codo(&[10.0, 9.0]), None);
    assert_eq!(k_del_codo(&[10.0]), None);
    assert_eq!(k_del_codo(&[]), None);
    assert_eq!(k_del_codo(&[0.0, 0.0, 0.0]), None);
    assert_eq!(k_del_codo(&[1.0, 2.0, 3.0]), None);
}

#[test]
fn el_codo_no_se_deja_llevar_por_un_grupo_muy_lejano() {
    let curva = [15.0, 3.8, 1.5, 1.08, 0.72, 0.52, 0.48, 0.41, 0.38, 0.33];
    assert_eq!(k_del_codo(&curva), Some(3));
}

#[test]
fn componentes_principales_en_cualquier_direccion() {
    let datos = puntos(&[&[0.0, 0.0], &[3.0, -2.0], &[6.0, -4.0]]);
    let componentes = componentes_principales(&datos).unwrap();
    let largo = 13.0f64.sqrt();
    assert!(cerca(componentes.primera[0], 3.0 / largo));
    assert!(cerca(componentes.primera[1], -2.0 / largo));
    let mut azar = Azar::nuevo(11);
    let nube: Vec<Vec<f64>> = (0..200)
        .map(|_| {
            let a = azar.normal();
            let b = azar.normal();
            vec![
                a + 0.3 * b,
                a - 0.3 * b,
                1.2 * azar.normal(),
                0.1 * azar.normal(),
            ]
        })
        .collect();
    let componentes = componentes_principales(&nube).unwrap();
    let producto: f64 = componentes
        .primera
        .iter()
        .zip(&componentes.segunda)
        .map(|(x, y)| x * y)
        .sum();
    assert!(producto.abs() < 1e-9);
    assert!(componentes.varianza_primera >= componentes.varianza_segunda);
    assert!(
        componentes.varianza_primera + componentes.varianza_segunda
            <= componentes.varianza_total + 1e-9
    );
    assert!(componentes.primera[0] > 0.6 && componentes.primera[1] > 0.6);
    assert!(componentes.segunda[2].abs() > 0.9);
}

#[test]
fn perfil_constante_con_decimales_no_se_correlaciona_con_nada() {
    for constante in [0.1, 3.3, 1e5 + 0.1, 1e-3, 7.7] {
        let plano = vec![constante; 6];
        let otro = [1.0, 2.0, 3.5, 0.0, 5.0, 7.0];
        assert_eq!(correlacion_pearson(&plano, &otro), 0.0);
        assert_eq!(correlacion_pearson(&plano, &plano), 0.0);
        assert_eq!(correlacion_pearson(&plano, &[0.3; 6]), 0.0);
        assert_eq!(distancia_pearson(&plano, &otro), 1.0);
    }
    let perfiles = puntos(&[&[0.1, 0.1, 0.1], &[3.3, 3.3, 3.3], &[1.0, 2.0, 3.0]]);
    let correlaciones = matriz_de_correlaciones(&perfiles);
    assert_eq!(correlaciones[0][1], 0.0);
    assert_eq!(correlaciones[0][2], 0.0);
    assert!(cerca(correlaciones[2][2], 1.0));
}

#[test]
fn correlacion_de_valores_diminutos_sigue_siendo_valida() {
    let a = [1e-200, 2e-200, 3e-200];
    let b = [3.0, 2.0, 1.0];
    assert!(cerca(correlacion_pearson(&a, &a), 1.0));
    assert!(cerca(correlacion_pearson(&a, &b), -1.0));
}
