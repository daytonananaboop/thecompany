use nucleo::azar::Azar;
use nucleo::dcj::*;
use nucleo::reordenamientos::distancia_dos_roturas;
use std::collections::{HashMap, VecDeque};

fn genoma(texto: &str) -> Vec<Cromosoma> {
    leer_genoma(texto, 1_000_000).unwrap()
}

fn distancia(a: &str, b: &str) -> usize {
    distancia_dcj(&genoma(a), &genoma(b)).unwrap().distancia
}

fn tipos(a: &str, b: &str) -> Vec<TipoDcj> {
    ordenar(&genoma(a), &genoma(b))
        .unwrap()
        .operaciones
        .iter()
        .map(|operacion| operacion.tipo)
        .collect()
}

fn vertices(pareja: &[usize]) -> Vec<Vertice> {
    let mut lista = Vec::new();
    for (p, &q) in pareja.iter().enumerate().skip(1) {
        if q == TELOMERO {
            lista.push((p, TELOMERO));
        } else if p < q {
            lista.push((p, q));
        }
    }
    lista
}

fn con_cambio(pareja: &[usize], cortadas: [Vertice; 2], nuevas: [Vertice; 2]) -> Vec<usize> {
    let mut copia = pareja.to_vec();
    aplicar(&mut copia, &cortadas, &nuevas);
    copia
}

fn vecinos(pareja: &[usize]) -> Vec<Vec<usize>> {
    let lista = vertices(pareja);
    let mut salida = Vec::new();
    for (i, &u) in lista.iter().enumerate() {
        if u.1 != TELOMERO {
            salida.push(con_cambio(
                pareja,
                [u, NADA],
                [(u.0, TELOMERO), (u.1, TELOMERO)],
            ));
        }
        for &v in &lista[i + 1..] {
            match (u.1 == TELOMERO, v.1 == TELOMERO) {
                (true, true) => salida.push(con_cambio(pareja, [u, v], [(u.0, v.0), NADA])),
                (false, true) => {
                    salida.push(con_cambio(pareja, [u, v], [(u.0, v.0), (u.1, TELOMERO)]));
                    salida.push(con_cambio(pareja, [u, v], [(u.1, v.0), (u.0, TELOMERO)]));
                }
                (true, false) => {
                    salida.push(con_cambio(pareja, [u, v], [(v.0, u.0), (v.1, TELOMERO)]));
                    salida.push(con_cambio(pareja, [u, v], [(v.1, u.0), (v.0, TELOMERO)]));
                }
                (false, false) => {
                    salida.push(con_cambio(pareja, [u, v], [(u.0, v.0), (u.1, v.1)]));
                    salida.push(con_cambio(pareja, [u, v], [(u.0, v.1), (u.1, v.0)]));
                }
            }
        }
    }
    salida
}

fn distancias_por_bfs(inicio: &[usize]) -> HashMap<Vec<usize>, usize> {
    let mut distancias = HashMap::new();
    distancias.insert(inicio.to_vec(), 0);
    let mut cola = VecDeque::from([inicio.to_vec()]);
    while let Some(actual) = cola.pop_front() {
        let siguiente = distancias[&actual] + 1;
        for vecino in vecinos(&actual) {
            if !distancias.contains_key(&vecino) {
                distancias.insert(vecino.clone(), siguiente);
                cola.push_back(vecino);
            }
        }
    }
    distancias
}

fn genoma_al_azar(azar: &mut Azar, bloques: usize) -> Vec<Cromosoma> {
    let mut orden: Vec<i64> = (1..=bloques as i64).collect();
    for i in (1..bloques).rev() {
        let j = azar.entero(i + 1);
        orden.swap(i, j);
    }
    let mut genoma: Vec<Cromosoma> = Vec::new();
    for bloque in orden {
        let con_signo = if azar.entero(2) == 0 { bloque } else { -bloque };
        let nuevo = genoma.is_empty() || azar.entero(4) == 0;
        if nuevo {
            genoma.push(Cromosoma {
                lineal: azar.entero(2) == 0,
                bloques: vec![con_signo],
            });
        } else {
            let ultimo = genoma.len() - 1;
            genoma[ultimo].bloques.push(con_signo);
        }
    }
    genoma
}

fn etiquetas_seguidas(bloques: usize) -> Vec<u64> {
    (1..=bloques as u64).collect()
}

fn resumen_de_cromosomas(genoma: &[Cromosoma]) -> (usize, usize) {
    (
        genoma.len(),
        genoma.iter().filter(|cromosoma| cromosoma.lineal).count(),
    )
}

fn revisar_escenario(a: &[Cromosoma], b: &[Cromosoma]) {
    let ordenamiento = ordenar(a, b).unwrap();
    let etiquetas = &ordenamiento.etiquetas;
    assert_eq!(
        ordenamiento.operaciones.len(),
        ordenamiento.conteo.distancia
    );
    let mut pareja = ordenamiento.pareja_a.clone();
    let mut restante = ordenamiento.conteo.distancia;
    for operacion in &ordenamiento.operaciones {
        let genoma_antes = genoma_desde_adyacencias(&pareja, etiquetas);
        aplicar(&mut pareja, &operacion.cortadas, &operacion.nuevas);
        let genoma_despues = genoma_desde_adyacencias(&pareja, etiquetas);
        let nueva = distancia_dcj(&genoma_despues, b).unwrap().distancia;
        assert_eq!(nueva + 1, restante);
        restante = nueva;
        let (total_antes, lineales_antes) = resumen_de_cromosomas(&genoma_antes);
        let (total_despues, lineales_despues) = resumen_de_cromosomas(&genoma_despues);
        match operacion.tipo {
            TipoDcj::Fision => assert_eq!(total_despues, total_antes + 1),
            TipoDcj::Fusion => assert_eq!(total_despues + 1, total_antes),
            TipoDcj::Linealizacion => {
                assert_eq!(total_despues, total_antes);
                assert_eq!(lineales_despues, lineales_antes + 1);
            }
            TipoDcj::Circularizacion => {
                assert_eq!(total_despues, total_antes);
                assert_eq!(lineales_despues + 1, lineales_antes);
            }
            TipoDcj::Inversion | TipoDcj::Translocacion => {
                assert_eq!(total_despues, total_antes);
                assert_eq!(lineales_despues, lineales_antes);
            }
        }
    }
    assert_eq!(pareja, ordenamiento.pareja_b);
}

#[test]
fn lee_lineales_y_circulares() {
    let leido = genoma("[+1 -2 3] (+4, −5)");
    assert_eq!(
        leido,
        vec![
            Cromosoma {
                lineal: true,
                bloques: vec![1, -2, 3]
            },
            Cromosoma {
                lineal: false,
                bloques: vec![4, -5]
            },
        ]
    );
}

#[test]
fn rechaza_textos_mal_escritos() {
    for malo in [
        "",
        "   ",
        "+1 +2",
        "[+1 +2",
        "[+1 +2)",
        "(+1 [+2])",
        "[]",
        "[+1 x]",
        "]",
        "[+1] +2",
    ] {
        assert!(leer_genoma(malo, 100).is_err(), "{malo}");
    }
    assert!(leer_genoma("[+1 +2 +3]", 2).is_err());
}

#[test]
fn rechaza_bloques_invalidos() {
    let error = distancia_dcj(&genoma("[+1 0]"), &genoma("[+1 +2]")).unwrap_err();
    assert!(error.contains("bloque 0"), "{error}");
    let error = distancia_dcj(&genoma("[+1 -1]"), &genoma("[+1 +2]")).unwrap_err();
    assert!(error.contains("más de una vez"), "{error}");
    let error = distancia_dcj(&genoma("[+1 +2]"), &genoma("[+1 +3]")).unwrap_err();
    assert!(error.contains("bloque 2"), "{error}");
    let error = distancia_dcj(&genoma("[+1 +2]"), &genoma("[+1 +2 +3]")).unwrap_err();
    assert!(error.contains("bloque 3"), "{error}");
    assert!(ordenar(&[], &genoma("[+1]")).is_err());
}

#[test]
fn genomas_iguales_estan_a_distancia_cero() {
    let a = "[+1 +2 -3] (+4 +5) [+6]";
    assert_eq!(distancia(a, a), 0);
    assert_eq!(distancia(a, "[+3 -2 -1] (+5 +4) [-6]"), 0);
    assert!(tipos(a, a).is_empty());
}

#[test]
fn una_inversion() {
    let conteo = distancia_dcj(&genoma("[+1 +2 +3]"), &genoma("[+1 -2 +3]")).unwrap();
    assert_eq!(
        conteo,
        Conteo {
            bloques: 3,
            ciclos: 1,
            caminos_impares: 2,
            caminos_pares: 0,
            distancia: 1
        }
    );
    assert_eq!(tipos("[+1 +2 +3]", "[+1 -2 +3]"), vec![TipoDcj::Inversion]);
}

#[test]
fn fision_y_fusion() {
    let conteo = distancia_dcj(&genoma("[+1 +2]"), &genoma("[+1] [+2]")).unwrap();
    assert_eq!(conteo.ciclos, 0);
    assert_eq!(conteo.caminos_impares, 2);
    assert_eq!(conteo.caminos_pares, 1);
    assert_eq!(conteo.distancia, 1);
    assert_eq!(tipos("[+1 +2]", "[+1] [+2]"), vec![TipoDcj::Fision]);
    assert_eq!(tipos("[+1] [+2]", "[+1 +2]"), vec![TipoDcj::Fusion]);
}

#[test]
fn circularizacion_y_linealizacion() {
    let conteo = distancia_dcj(&genoma("[+1 +2]"), &genoma("(+1 +2)")).unwrap();
    assert_eq!(conteo.ciclos, 1);
    assert_eq!(conteo.caminos_impares, 0);
    assert_eq!(conteo.caminos_pares, 1);
    assert_eq!(conteo.distancia, 1);
    assert_eq!(tipos("[+1 +2]", "(+1 +2)"), vec![TipoDcj::Circularizacion]);
    assert_eq!(tipos("(+1 +2)", "[+1 +2]"), vec![TipoDcj::Linealizacion]);
    assert_eq!(tipos("(+1 +2)", "[+2 +1]"), vec![TipoDcj::Linealizacion]);
}

#[test]
fn una_translocacion() {
    assert_eq!(distancia("[+1 +2] [+3 +4]", "[+1 +4] [+3 +2]"), 1);
    assert_eq!(
        tipos("[+1 +2] [+3 +4]", "[+1 +4] [+3 +2]"),
        vec![TipoDcj::Translocacion]
    );
}

#[test]
fn escision_de_un_circulo() {
    assert_eq!(distancia("[+1 +2 +3]", "[+1 +3] (+2)"), 1);
    assert_eq!(tipos("[+1 +2 +3]", "[+1 +3] (+2)"), vec![TipoDcj::Fision]);
    assert_eq!(tipos("[+1 +3] (+2)", "[+1 +2 +3]"), vec![TipoDcj::Fusion]);
}

#[test]
fn etiquetas_salteadas() {
    assert_eq!(distancia("[+10 +20 +30]", "[+10 -20 +30]"), 1);
    let ordenamiento = ordenar(&genoma("[+10 +20 +30]"), &genoma("(+30 +20 +10)")).unwrap();
    assert_eq!(ordenamiento.etiquetas, vec![10, 20, 30]);
    assert_eq!(nombre_del_extremo(1, &ordenamiento.etiquetas), "10t");
    assert_eq!(nombre_del_extremo(6, &ordenamiento.etiquetas), "30h");
    revisar_escenario(&genoma("[+10 +20 +30]"), &genoma("(+30 +20 +10)"));
}

#[test]
fn un_caso_mezclado_hecho_a_mano() {
    let a = genoma("[+1 +2 +3 +4] (+5 +6)");
    let b = genoma("[+1 -3 -2] (+4) [+5 +6]");
    let ordenamiento = ordenar(&a, &b).unwrap();
    assert_eq!(ordenamiento.conteo.bloques, 6);
    assert_eq!(ordenamiento.conteo.distancia, 3);
    revisar_escenario(&a, &b);
}

#[test]
fn el_grafo_reparte_todos_los_extremos() {
    let mut azar = Azar::nuevo(5);
    for _ in 0..200 {
        let bloques = 1 + azar.entero(30);
        let a = genoma_al_azar(&mut azar, bloques);
        let b = genoma_al_azar(&mut azar, bloques);
        let etiquetas = etiquetas_seguidas(bloques);
        let pareja_a = adyacencias(&a, &etiquetas).unwrap();
        let pareja_b = adyacencias(&b, &etiquetas).unwrap();
        let componentes = grafo_de_adyacencias(&pareja_a, &pareja_b);
        let mut vistos = vec![0; 2 * bloques + 1];
        for componente in &componentes {
            for &extremo in &componente.extremos {
                vistos[extremo] += 1;
            }
            let impar = componente.extremos.len() % 2 == 1;
            assert_eq!(impar, componente.clase == Clase::CaminoAB);
        }
        assert!(vistos[1..].iter().all(|&veces| veces == 1));
        let conteo = contar(bloques, &componentes);
        assert_eq!(conteo.caminos_impares % 2, 0);
    }
}

#[test]
fn ida_y_vuelta_entre_genoma_y_adyacencias() {
    let mut azar = Azar::nuevo(11);
    for _ in 0..300 {
        let bloques = 1 + azar.entero(40);
        let a = genoma_al_azar(&mut azar, bloques);
        let etiquetas = etiquetas_seguidas(bloques);
        let pareja = adyacencias(&a, &etiquetas).unwrap();
        let reconstruido = genoma_desde_adyacencias(&pareja, &etiquetas);
        assert_eq!(adyacencias(&reconstruido, &etiquetas).unwrap(), pareja);
        assert_eq!(distancia_dcj(&a, &reconstruido).unwrap().distancia, 0);
    }
}

#[test]
fn coincide_con_la_busqueda_exhaustiva_hasta_cinco_bloques() {
    let mut azar = Azar::nuevo(23);
    for bloques in 1..=5 {
        let etiquetas = etiquetas_seguidas(bloques);
        for _ in 0..3 {
            let a = genoma_al_azar(&mut azar, bloques);
            let pareja_a = adyacencias(&a, &etiquetas).unwrap();
            let todas = distancias_por_bfs(&pareja_a);
            for (pareja_b, esperada) in &todas {
                let b = genoma_desde_adyacencias(pareja_b, &etiquetas);
                assert_eq!(
                    distancia_dcj(&a, &b).unwrap().distancia,
                    *esperada,
                    "{a:?} {b:?}"
                );
            }
        }
    }
}

#[test]
fn coincide_con_la_busqueda_exhaustiva_con_seis_bloques() {
    let mut azar = Azar::nuevo(29);
    let etiquetas = etiquetas_seguidas(6);
    let a = genoma_al_azar(&mut azar, 6);
    let pareja_a = adyacencias(&a, &etiquetas).unwrap();
    let todas = distancias_por_bfs(&pareja_a);
    assert_eq!(todas.len(), 140_152);
    let mut revisados = 0;
    for (pareja_b, esperada) in &todas {
        let b = genoma_desde_adyacencias(pareja_b, &etiquetas);
        assert_eq!(distancia_dcj(&a, &b).unwrap().distancia, *esperada);
        if revisados < 2_000 {
            revisar_escenario(&a, &b);
            revisados += 1;
        }
    }
}

#[test]
fn cada_paso_del_escenario_baja_la_distancia_en_uno() {
    let mut azar = Azar::nuevo(31);
    for _ in 0..300 {
        let bloques = 1 + azar.entero(25);
        let a = genoma_al_azar(&mut azar, bloques);
        let b = genoma_al_azar(&mut azar, bloques);
        revisar_escenario(&a, &b);
    }
}

#[test]
fn con_todo_circular_coincide_con_las_dos_roturas() {
    let mut azar = Azar::nuevo(37);
    for _ in 0..300 {
        let bloques = 1 + azar.entero(40);
        let mut a = genoma_al_azar(&mut azar, bloques);
        let mut b = genoma_al_azar(&mut azar, bloques);
        for cromosoma in a.iter_mut().chain(b.iter_mut()) {
            cromosoma.lineal = false;
        }
        let p: Vec<Vec<i64>> = a.iter().map(|c| c.bloques.clone()).collect();
        let q: Vec<Vec<i64>> = b.iter().map(|c| c.bloques.clone()).collect();
        let conteo = distancia_dcj(&a, &b).unwrap();
        assert_eq!(conteo.caminos_impares + conteo.caminos_pares, 0);
        assert_eq!(conteo.distancia, distancia_dos_roturas(&p, &q).unwrap());
    }
}

#[test]
fn genomas_grandes_sin_problemas() {
    let mut azar = Azar::nuevo(41);
    let a = genoma_al_azar(&mut azar, 200_000);
    let b = genoma_al_azar(&mut azar, 200_000);
    let conteo = distancia_dcj(&a, &b).unwrap();
    assert!(conteo.distancia <= 200_000);
    let a = genoma_al_azar(&mut azar, 3_000);
    let b = genoma_al_azar(&mut azar, 3_000);
    let ordenamiento = ordenar(&a, &b).unwrap();
    assert_eq!(
        ordenamiento.operaciones.len(),
        ordenamiento.conteo.distancia
    );
}

#[test]
fn el_bloque_mas_pequeno_de_i64_se_rechaza_y_el_vecino_funciona() {
    let a = genoma("[-9223372036854775808 +1]");
    let b = genoma("[+1 -9223372036854775808]");
    let error = distancia_dcj(&a, &b).unwrap_err();
    assert!(error.contains("demasiado grande"), "{error}");
    let a = genoma("[-9223372036854775807 +1]");
    let b = genoma("[+1 -9223372036854775807]");
    let ordenamiento = ordenar(&a, &b).unwrap();
    assert_eq!(
        ordenamiento.operaciones.len(),
        ordenamiento.conteo.distancia
    );
    let final_ = genoma_desde_adyacencias(&ordenamiento.pareja_b, &ordenamiento.etiquetas);
    assert_eq!(final_, b);
}
