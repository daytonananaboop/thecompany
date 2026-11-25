use nucleo::agrupamiento::Medida;
use nucleo::azar::Azar;
use nucleo::comunidades::*;

fn cerca(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

fn puntos(filas: &[&[f64]]) -> Vec<Vec<f64>> {
    filas.iter().map(|fila| fila.to_vec()).collect()
}

fn grafo_de(nodos: usize, aristas: &[(usize, usize, f64)]) -> Grafo {
    let mut grafo = grafo_vacio(nodos);
    for &(a, b, peso) in aristas {
        agregar_arista(&mut grafo, a, b, peso);
    }
    for lista in &mut grafo.vecinos {
        lista.sort_unstable_by_key(|&(otro, _)| otro);
    }
    grafo
}

fn dos_triangulos() -> Grafo {
    grafo_de(
        6,
        &[
            (0, 1, 1.0),
            (1, 2, 1.0),
            (0, 2, 1.0),
            (3, 4, 1.0),
            (4, 5, 1.0),
            (3, 5, 1.0),
            (2, 3, 1.0),
        ],
    )
}

fn grafo_al_azar(azar: &mut Azar, nodos: usize) -> Grafo {
    let mut aristas = Vec::new();
    for a in 0..nodos {
        for b in a + 1..nodos {
            if azar.decimal() < 0.4 {
                aristas.push((a, b, 0.1 + azar.decimal()));
            }
        }
    }
    if aristas.is_empty() {
        aristas.push((0, 1, 1.0));
    }
    grafo_de(nodos, &aristas)
}

fn nubes(azar: &mut Azar, por_nube: usize, centros: &[[f64; 3]]) -> (Vec<Vec<f64>>, Vec<usize>) {
    let mut datos = Vec::new();
    let mut verdad = Vec::new();
    for (grupo, centro) in centros.iter().enumerate() {
        for _ in 0..por_nube {
            datos.push(
                centro
                    .iter()
                    .map(|valor| valor + 0.3 * azar.normal())
                    .collect(),
            );
            verdad.push(grupo);
        }
    }
    (datos, verdad)
}

fn todas_las_particiones(nodos: usize) -> Vec<Vec<usize>> {
    let mut resultado = Vec::new();
    let mut actual = vec![0usize; nodos];
    fn crecer(
        posicion: usize,
        mayor: usize,
        actual: &mut Vec<usize>,
        resultado: &mut Vec<Vec<usize>>,
    ) {
        if posicion == actual.len() {
            resultado.push(actual.clone());
            return;
        }
        for etiqueta in 0..=mayor + 1 {
            actual[posicion] = etiqueta;
            crecer(posicion + 1, mayor.max(etiqueta), actual, resultado);
        }
    }
    if nodos > 0 {
        crecer(1, 0, &mut actual, &mut resultado);
    }
    resultado
}

#[test]
fn vecinos_calculados_a_mano() {
    let datos = puntos(&[&[0.0], &[1.0], &[3.0], &[7.0]]);
    let vecinos = vecinos_mas_cercanos(&datos, 2);
    assert_eq!(
        vecinos,
        vec![vec![1, 2], vec![0, 2], vec![1, 0], vec![2, 1]]
    );
}

#[test]
fn los_empates_se_resuelven_por_el_indice_menor() {
    let datos = puntos(&[&[0.0], &[1.0], &[-1.0], &[2.0], &[-2.0]]);
    let vecinos = vecinos_mas_cercanos(&datos, 3);
    assert_eq!(vecinos[0], vec![1, 2, 3]);
}

#[test]
fn pearson_junta_perfiles_de_la_misma_forma() {
    let datos = puntos(&[
        &[1.0, 2.0, 3.0, 4.0],
        &[10.0, 20.0, 30.0, 40.0],
        &[1.2, 1.2, 1.2, 1.3],
        &[4.0, 3.0, 2.0, 1.0],
    ]);
    let filas = filas_para_medir(&datos, Medida::Pearson);
    assert_eq!(vecinos_mas_cercanos(&filas, 1)[0], vec![1]);
    let euclideas = filas_para_medir(&datos, Medida::Euclidea);
    assert_eq!(vecinos_mas_cercanos(&euclideas, 1)[0], vec![2]);
    let cuadrado = nucleo::agrupamiento::distancia_al_cuadrado(&filas[0], &filas[3]);
    assert!(cerca(distancia_medida(cuadrado, Medida::Pearson), 2.0));
    let plano = filas_para_medir(&puntos(&[&[5.0, 5.0, 5.0]]), Medida::Pearson);
    assert_eq!(plano[0], vec![0.0, 0.0, 0.0, 1.0]);
}

#[test]
fn jaccard_de_vecindarios() {
    assert_eq!(compartidos(&[1, 3, 5, 7], &[2, 3, 4, 7]), 2);
    assert_eq!(compartidos(&[], &[1]), 0);
    assert!(cerca(peso_jaccard(&[1, 3, 5, 7], &[2, 3, 4, 7]), 2.0 / 6.0));
    assert!(cerca(peso_jaccard(&[1, 2], &[1, 2]), 1.0));
    assert_eq!(vecindario_con_el_mismo(&[5, 1, 3], 2), vec![1, 2, 3, 5]);
}

#[test]
fn el_grafo_snn_une_solo_a_los_que_comparten_vecinos() {
    let datos = puntos(&[
        &[0.0, 0.0],
        &[0.1, 0.0],
        &[0.0, 0.1],
        &[10.0, 10.0],
        &[10.1, 10.0],
        &[10.0, 10.1],
    ]);
    let vecinos = vecinos_mas_cercanos(&datos, 2);
    let grafo = grafo_snn(&vecinos);
    for nodo in 0..3 {
        assert!(grafo.vecinos[nodo].iter().all(|&(otro, _)| otro < 3));
    }
    for nodo in 3..6 {
        assert!(grafo.vecinos[nodo].iter().all(|&(otro, _)| otro >= 3));
    }
    assert!(grafo.vecinos[0].iter().all(|&(_, peso)| cerca(peso, 1.0)));
    assert!(cerca(grafo.peso_total, 6.0));
    let grados: f64 = grafo.grados.iter().sum();
    assert!(cerca(grados, 2.0 * grafo.peso_total));
}

#[test]
fn modularidad_de_dos_triangulos_a_mano() {
    let grafo = dos_triangulos();
    let q = modularidad(&grafo, &[0, 0, 0, 1, 1, 1], 1.0);
    assert!(cerca(q, 6.0 / 7.0 - 0.5));
    assert!(cerca(modularidad(&grafo, &[0; 6], 1.0), 0.0));
    let q2 = modularidad(&grafo, &[0, 0, 0, 1, 1, 1], 2.0);
    assert!(cerca(q2, 6.0 / 7.0 - 1.0));
    assert_eq!(modularidad(&grafo_vacio(3), &[0, 1, 2], 1.0), 0.0);
}

#[test]
fn louvain_separa_los_dos_triangulos() {
    let resultado = louvain(&dos_triangulos(), 1.0, 1).unwrap();
    assert_eq!(resultado.cantidad, 2);
    assert_eq!(resultado.comunidad[0], resultado.comunidad[1]);
    assert_eq!(resultado.comunidad[1], resultado.comunidad[2]);
    assert_eq!(resultado.comunidad[3], resultado.comunidad[5]);
    assert_ne!(resultado.comunidad[0], resultado.comunidad[3]);
    assert!(cerca(resultado.modularidad, 6.0 / 7.0 - 0.5));
    assert!(!resultado.niveles.is_empty());
}

#[test]
fn louvain_sin_aristas_deja_cada_nodo_solo() {
    let resultado = louvain(&grafo_vacio(4), 1.0, 3).unwrap();
    assert_eq!(resultado.cantidad, 4);
    assert_eq!(resultado.modularidad, 0.0);
    assert!(resultado.niveles.is_empty());
}

#[test]
fn la_resolucion_debe_ser_positiva() {
    assert!(louvain(&dos_triangulos(), 0.0, 1).is_err());
    assert!(louvain(&dos_triangulos(), f64::NAN, 1).is_err());
}

#[test]
fn la_ganancia_coincide_con_la_modularidad_recalculada() {
    let mut azar = Azar::nuevo(21);
    for _ in 0..200 {
        let nodos = 3 + azar.entero(10);
        let grafo = grafo_al_azar(&mut azar, nodos);
        let mut particion = particion_inicial(&grafo);
        let mut peso_hacia = vec![0.0; nodos];
        let mut tocadas = Vec::new();
        for _ in 0..nodos * 2 {
            let nodo = azar.entero(nodos);
            let antes = modularidad(&grafo, &particion.comunidad, 1.3);
            let movimiento = mover_un_nodo(
                &grafo,
                &mut particion,
                nodo,
                1.3,
                &mut peso_hacia,
                &mut tocadas,
            );
            let despues = modularidad(&grafo, &particion.comunidad, 1.3);
            assert!((despues - antes - movimiento.ganancia).abs() < 1e-9);
            assert!(movimiento.ganancia >= 0.0);
            assert!(peso_hacia.iter().all(|&peso| peso == 0.0));
            let mut grados = vec![0.0; nodos];
            for (otro, &c) in particion.comunidad.iter().enumerate() {
                grados[c] += grafo.grados[otro];
            }
            for (calculado, guardado) in grados.iter().zip(&particion.grado_de_comunidad) {
                assert!((calculado - guardado).abs() < 1e-9);
            }
        }
    }
}

#[test]
fn agregar_comunidades_conserva_la_modularidad() {
    let mut azar = Azar::nuevo(8);
    for _ in 0..200 {
        let nodos = 2 + azar.entero(12);
        let grafo = grafo_al_azar(&mut azar, nodos);
        let etiquetas: Vec<usize> = (0..nodos).map(|_| azar.entero(4)).collect();
        let (comunidad, cantidad) = renumerar(&etiquetas);
        let agregado = agregar_comunidades(&grafo, &comunidad, cantidad);
        let solos: Vec<usize> = (0..cantidad).collect();
        assert!((agregado.peso_total - grafo.peso_total).abs() < 1e-9);
        assert!(
            (modularidad(&grafo, &comunidad, 0.8) - modularidad(&agregado, &solos, 0.8)).abs()
                < 1e-9
        );
    }
}

#[test]
fn louvain_contra_todas_las_particiones_posibles() {
    let mut azar = Azar::nuevo(99);
    for _ in 0..60 {
        let nodos = 3 + azar.entero(6);
        let grafo = grafo_al_azar(&mut azar, nodos);
        let resolucion = [0.5, 1.0, 1.5][azar.entero(3)];
        let resultado = louvain(&grafo, resolucion, azar.siguiente()).unwrap();
        let mejor = todas_las_particiones(nodos)
            .iter()
            .map(|particion| modularidad(&grafo, particion, resolucion))
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(resultado.modularidad <= mejor + 1e-9);
        assert!(resultado.modularidad >= mejor - 0.1);
        assert!(cerca(
            resultado.modularidad,
            modularidad(&grafo, &resultado.comunidad, resolucion)
        ));
        let mut tamanos = vec![0; resultado.cantidad];
        for &c in &resultado.comunidad {
            tamanos[c] += 1;
        }
        assert!(tamanos.windows(2).all(|par| par[0] >= par[1]));
        for primera in 0..resultado.cantidad {
            for segunda in primera + 1..resultado.cantidad {
                let unidas: Vec<usize> = resultado
                    .comunidad
                    .iter()
                    .map(|&c| if c == segunda { primera } else { c })
                    .collect();
                assert!(modularidad(&grafo, &unidas, resolucion) <= resultado.modularidad + 1e-9);
            }
        }
    }
}

#[test]
fn las_comunidades_recuperan_nubes_separadas() {
    let mut azar = Azar::nuevo(4);
    let (datos, verdad) = nubes(
        &mut azar,
        30,
        &[[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [0.0, 5.0, 5.0]],
    );
    let resultado = agrupar_por_comunidades(&datos, Medida::Euclidea, 10, 1.0, 1).unwrap();
    assert!(resultado.louvain.cantidad >= 3);
    for comunidad in 0..resultado.louvain.cantidad {
        let nubes_dentro: std::collections::HashSet<usize> = verdad
            .iter()
            .zip(&resultado.louvain.comunidad)
            .filter(|&(_, &c)| c == comunidad)
            .map(|(&nube, _)| nube)
            .collect();
        assert_eq!(nubes_dentro.len(), 1);
    }
    assert!(resultado.louvain.modularidad > 0.5);
    let otra = agrupar_por_comunidades(&datos, Medida::Euclidea, 10, 1.0, 1).unwrap();
    assert_eq!(otra, resultado);
}

#[test]
fn mas_resolucion_da_mas_comunidades() {
    let mut azar = Azar::nuevo(6);
    let (datos, _) = nubes(&mut azar, 40, &[[0.0, 0.0, 0.0], [6.0, 0.0, 0.0]]);
    let pocas = agrupar_por_comunidades(&datos, Medida::Euclidea, 10, 0.2, 1).unwrap();
    let muchas = agrupar_por_comunidades(&datos, Medida::Euclidea, 10, 4.0, 1).unwrap();
    assert!(muchas.louvain.cantidad > pocas.louvain.cantidad);
}

#[test]
fn errores_de_entrada_de_las_comunidades() {
    let datos = puntos(&[&[0.0], &[1.0], &[2.0], &[3.0]]);
    assert!(agrupar_por_comunidades(&datos, Medida::Euclidea, 0, 1.0, 1).is_err());
    assert!(agrupar_por_comunidades(&datos, Medida::Euclidea, 4, 1.0, 1).is_err());
    assert!(agrupar_por_comunidades(&datos[..2], Medida::Euclidea, 1, 1.0, 1).is_err());
    assert!(agrupar_por_comunidades(&datos, Medida::Euclidea, 2, -1.0, 1).is_err());
    assert!(agrupar_por_comunidades(&[], Medida::Euclidea, 2, 1.0, 1).is_err());
    let malos = puntos(&[&[0.0], &[f64::NAN], &[2.0], &[3.0]]);
    assert!(agrupar_por_comunidades(&malos, Medida::Euclidea, 2, 1.0, 1).is_err());
    let iguales = puntos(&[&[1.0], &[1.0], &[1.0], &[1.0]]);
    assert!(agrupar_por_comunidades(&iguales, Medida::Euclidea, 2, 1.0, 1).is_ok());
}

#[test]
fn rand_ajustado_a_mano() {
    assert!(cerca(
        indice_rand_ajustado(&[0, 0, 1, 1], &[1, 1, 0, 0]),
        1.0
    ));
    assert!(cerca(
        indice_rand_ajustado(&[0, 0, 1, 1], &[0, 1, 0, 1]),
        -0.5
    ));
    assert!(cerca(
        indice_rand_ajustado(&[0, 0, 0, 1, 1, 1], &[0, 0, 1, 1, 2, 2]),
        0.24242424242424243
    ));
    assert!(cerca(indice_rand_ajustado(&[0, 0, 0], &[0, 0, 0]), 1.0));
}

#[test]
fn la_calibracion_alcanza_la_perplejidad() {
    let mut azar = Azar::nuevo(12);
    for _ in 0..50 {
        let n = 10 + azar.entero(40);
        let distancias: Vec<f64> = (0..n).map(|_| azar.decimal() * 20.0).collect();
        let perplejidad = 2.0 + azar.decimal() * (n as f64 / 3.0 - 2.0);
        let mut probabilidades = vec![0.0; n];
        let calibracion = calibrar_un_punto(&distancias, 0, perplejidad, &mut probabilidades);
        assert!((calibracion.entropia.exp() - perplejidad).abs() < 1e-3 * perplejidad);
        assert_eq!(probabilidades[0], 0.0);
        assert!(cerca(probabilidades.iter().sum::<f64>(), 1.0));
        let mut otra = vec![0.0; n];
        let h = entropia_con_beta(&distancias, 0, calibracion.beta, &mut otra);
        let directa: f64 = otra
            .iter()
            .filter(|&&p| p > 0.0)
            .map(|&p| -p * p.ln())
            .sum();
        assert!((h - directa).abs() < 1e-9);
    }
}

#[test]
fn las_afinidades_son_simetricas_y_suman_uno() {
    let mut azar = Azar::nuevo(3);
    let (datos, _) = nubes(&mut azar, 10, &[[0.0, 0.0, 0.0], [4.0, 4.0, 0.0]]);
    let n = datos.len();
    let (conjuntas, sigmas) = afinidades(&datos, 5.0);
    assert_eq!(sigmas.len(), n);
    assert!(sigmas.iter().all(|&sigma| sigma > 0.0 && sigma.is_finite()));
    let suma: f64 = conjuntas.iter().sum();
    assert!((suma - 1.0).abs() < 1e-6);
    for i in 0..n {
        assert_eq!(conjuntas[i * n + i], 0.0);
        for j in 0..n {
            assert_eq!(conjuntas[i * n + j], conjuntas[j * n + i]);
        }
    }
}

#[test]
fn el_gradiente_coincide_con_derivar_la_divergencia() {
    let mut azar = Azar::nuevo(17);
    let (datos, _) = nubes(&mut azar, 6, &[[0.0, 0.0, 0.0], [3.0, 0.0, 1.0]]);
    let (conjuntas, _) = afinidades(&datos, 3.0);
    let coordenadas: Vec<[f64; 2]> = (0..datos.len())
        .map(|_| [azar.normal(), azar.normal()])
        .collect();
    let analitico = gradiente(&conjuntas, &coordenadas, 1.0);
    let h = 1e-6;
    for i in 0..coordenadas.len() {
        for eje in 0..2 {
            let mut mas = coordenadas.clone();
            mas[i][eje] += h;
            let mut menos = coordenadas.clone();
            menos[i][eje] -= h;
            let numerico =
                (divergencia_kl(&conjuntas, &mas) - divergencia_kl(&conjuntas, &menos)) / (2.0 * h);
            assert!((numerico - analitico[i][eje]).abs() < 1e-5);
        }
    }
}

#[test]
fn tsne_separa_nubes_y_baja_la_divergencia() {
    let mut azar = Azar::nuevo(5);
    let (datos, verdad) = nubes(
        &mut azar,
        20,
        &[[0.0, 0.0, 0.0], [6.0, 0.0, 0.0], [0.0, 6.0, 6.0]],
    );
    let opciones = OpcionesTsne {
        perplejidad: 10.0,
        exageracion: 12.0,
        iteraciones: 600,
        semilla: 2,
    };
    let resultado = tsne(&datos, Medida::Euclidea, opciones).unwrap();
    assert_eq!(resultado.coordenadas.len(), datos.len());
    assert!(resultado.historia.first().unwrap().1 > resultado.divergencia);
    assert_eq!(resultado.historia.last().unwrap().0, 600);
    assert!(resultado.divergencia > 0.0 && resultado.divergencia < 1.0);
    let plano: Vec<Vec<f64>> = resultado.coordenadas.iter().map(|p| p.to_vec()).collect();
    let mas_cercano = vecinos_mas_cercanos(&plano, 1);
    for (punto, lista) in mas_cercano.iter().enumerate() {
        assert_eq!(verdad[lista[0]], verdad[punto]);
    }
    assert!(conservacion_de_vecinos(&datos, &resultado.coordenadas, 5) > 0.5);
    let otra = tsne(&datos, Medida::Euclidea, opciones).unwrap();
    assert_eq!(otra, resultado);
    let media: f64 = resultado.coordenadas.iter().map(|p| p[0]).sum::<f64>();
    assert!(media.abs() < 1e-6);
}

#[test]
fn errores_de_entrada_de_tsne() {
    let datos: Vec<Vec<f64>> = (0..10).map(|i| vec![i as f64, (i * i) as f64]).collect();
    let bien = OpcionesTsne {
        perplejidad: 3.0,
        exageracion: 12.0,
        iteraciones: 10,
        semilla: 1,
    };
    assert!(tsne(&datos, Medida::Euclidea, bien).is_ok());
    let alta = OpcionesTsne {
        perplejidad: 3.5,
        ..bien
    };
    assert!(tsne(&datos, Medida::Euclidea, alta).is_err());
    let baja = OpcionesTsne {
        perplejidad: 0.5,
        ..bien
    };
    assert!(tsne(&datos, Medida::Euclidea, baja).is_err());
    let sin_iteraciones = OpcionesTsne {
        iteraciones: 0,
        ..bien
    };
    assert!(tsne(&datos, Medida::Euclidea, sin_iteraciones).is_err());
    let sin_exagerar = OpcionesTsne {
        exageracion: 0.5,
        ..bien
    };
    assert!(tsne(&datos, Medida::Euclidea, sin_exagerar).is_err());
    assert!(tsne(&datos[..3], Medida::Euclidea, bien).is_err());
    let iguales = vec![vec![1.0, 1.0]; 10];
    let resultado = tsne(&iguales, Medida::Pearson, bien).unwrap();
    assert!(
        resultado
            .coordenadas
            .iter()
            .all(|p| p[0].is_finite() && p[1].is_finite())
    );
}

#[test]
fn conservar_vecinos_de_una_copia_exacta_es_uno() {
    let datos: Vec<Vec<f64>> = (0..12)
        .map(|i| vec![(i * 7 % 12) as f64, (i * i % 5) as f64])
        .collect();
    let mismo: Vec<[f64; 2]> = datos.iter().map(|fila| [fila[0], fila[1]]).collect();
    assert!(cerca(conservacion_de_vecinos(&datos, &mismo, 3), 1.0));
}

#[test]
fn con_velocidad_cero_las_ganancias_no_dependen_del_signo_del_gradiente() {
    let mut coordenadas = vec![[0.0, 0.0], [1.0, 1.0]];
    let mut velocidad = vec![[0.0, 0.0]; 2];
    let mut ganancias = vec![[1.0, 1.0]; 2];
    let gradiente = vec![[1.0, -1.0], [-2.0, 3.0]];
    descender(
        &mut coordenadas,
        &mut velocidad,
        &mut ganancias,
        &gradiente,
        0.5,
        10.0,
    );
    for fila in &ganancias {
        for &ganancia in fila {
            assert!(cerca(ganancia, 0.8));
        }
    }
    descender(
        &mut coordenadas,
        &mut velocidad,
        &mut ganancias,
        &gradiente,
        0.5,
        10.0,
    );
    for fila in &ganancias {
        for &ganancia in fila {
            assert!(cerca(ganancia, 1.0));
        }
    }
}

#[test]
fn copias_a_otra_escala_comparten_la_afinidad_por_igual_con_pearson() {
    let base = [1.0, 2.0, 4.0, 3.0, 5.0, 0.5];
    let mut azar = Azar::nuevo(11);
    let mut datos: Vec<Vec<f64>> = (1..=20)
        .map(|c| {
            base.iter()
                .map(|valor| valor * c as f64 * 0.37 + 0.1 * c as f64)
                .collect()
        })
        .collect();
    for _ in 0..15 {
        datos.push((0..6).map(|_| azar.normal()).collect());
    }
    let filas = filas_para_medir(&datos, Medida::Pearson);
    let (conjuntas, _) = afinidades(&filas, 5.0);
    let n = datos.len();
    let mut menor = f64::INFINITY;
    let mut mayor = 0.0f64;
    for i in 0..20 {
        for j in 0..20 {
            if i != j {
                menor = menor.min(conjuntas[i * n + j]);
                mayor = mayor.max(conjuntas[i * n + j]);
            }
        }
    }
    assert!(mayor / menor < 1.0 + 1e-6, "{menor} {mayor}");
}

#[test]
fn un_gen_plano_con_pearson_queda_a_correlacion_cero_de_todos() {
    let datos = puntos(&[
        &[1.0, 2.0, 3.0],
        &[3.0, 2.0, 1.0],
        &[5.0, 5.0, 5.0],
        &[7.0, 7.0, 7.0],
        &[1.0, 0.0, 0.0],
    ]);
    let filas = filas_para_medir(&datos, Medida::Pearson);
    let lejos = distancias_desde(&filas, 2);
    for (cuadrado, otro) in lejos {
        if otro == 3 {
            assert!(cerca(cuadrado, 0.0));
        } else {
            assert!(cerca(cuadrado, 2.0));
            assert!(cerca(distancia_medida(cuadrado, Medida::Pearson), 1.0));
        }
    }
    let normales = distancias_desde(&filas, 0);
    let con_pearson = |otro: usize| normales.iter().find(|&&(_, o)| o == otro).unwrap().0;
    assert!(cerca(con_pearson(1), 4.0));
}

#[test]
fn las_afinidades_no_dependen_de_la_escala_de_los_datos() {
    let mut azar = Azar::nuevo(21);
    let (datos, _) = nubes(&mut azar, 15, &[[0.0, 0.0, 0.0], [5.0, 5.0, 0.0]]);
    let (base, _) = afinidades(&datos, 8.0);
    for escala in [1e-40, 1e-20, 1e-8, 1e6, 1e11] {
        let cambiados: Vec<Vec<f64>> = datos
            .iter()
            .map(|fila| fila.iter().map(|valor| valor * escala).collect())
            .collect();
        let (otras, sigmas) = afinidades(&cambiados, 8.0);
        let diferencia = base
            .iter()
            .zip(&otras)
            .fold(0.0f64, |mayor, (a, b)| mayor.max((a - b).abs()));
        assert!(diferencia < 1e-6, "escala {escala}: {diferencia}");
        assert!(sigmas.iter().all(|sigma| sigma.is_finite() && *sigma > 0.0));
    }
}
