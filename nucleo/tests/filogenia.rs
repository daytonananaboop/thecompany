use nucleo::azar::Azar;
use nucleo::filogenia::{self, Arbol};

fn nombres(lista: &[&str]) -> Vec<String> {
    lista.iter().map(|nombre| nombre.to_string()).collect()
}

fn numerados(cantidad: usize) -> Vec<String> {
    (0..cantidad).map(|numero| numero.to_string()).collect()
}

fn matriz_abcd() -> Vec<Vec<f64>> {
    vec![
        vec![0.0, 5.0, 7.0, 11.0],
        vec![5.0, 0.0, 8.0, 12.0],
        vec![7.0, 8.0, 0.0, 6.0],
        vec![11.0, 12.0, 6.0, 0.0],
    ]
}

fn hermanas(arbol: &Arbol, primera: &str, segunda: &str) -> bool {
    let a = arbol.buscar_nodo(primera).unwrap();
    let b = arbol.buscar_nodo(segunda).unwrap();
    arbol.grado(a) == 1 && arbol.grado(b) == 1 && arbol.vecinos[a][0].0 == arbol.vecinos[b][0].0
}

fn parecidas(primera: &[Vec<f64>], segunda: &[Vec<f64>]) -> bool {
    primera.len() == segunda.len()
        && primera.iter().zip(segunda).all(|(fila_a, fila_b)| {
            fila_a.len() == fila_b.len()
                && fila_a.iter().zip(fila_b).all(|(a, b)| (a - b).abs() < 1e-9)
        })
}

fn arbol_con_pesos(cantidad: usize, semilla: u64) -> Arbol {
    let base = filogenia::arbol_al_azar(&numerados(cantidad), semilla).unwrap();
    let mut azar = Azar::nuevo(semilla + 1000);
    let mut arbol = Arbol::default();
    for nombre in &base.nombres {
        arbol.agregar_nodo(nombre);
    }
    for (a, b, _) in base.aristas() {
        arbol.conectar(a, b, 1.0 + azar.entero(9) as f64);
    }
    arbol
}

fn secuencias_al_azar(cantidad: usize, longitud: usize, semilla: u64) -> Vec<Vec<u8>> {
    let mut azar = Azar::nuevo(semilla);
    (0..cantidad)
        .map(|_| (0..longitud).map(|_| b"ACGT"[azar.entero(4)]).collect())
        .collect()
}

#[test]
fn lee_aristas_y_calcula_distancias_entre_hojas() {
    let arbol = filogenia::leer_aristas("a->x:2\nb->x:3\nx->y:4\nc->y:1\nd->y:5").unwrap();
    let hojas = arbol.hojas();
    let nombres_de_hojas: Vec<&str> = hojas
        .iter()
        .map(|&hoja| arbol.nombres[hoja].as_str())
        .collect();
    assert_eq!(nombres_de_hojas, ["a", "b", "c", "d"]);
    assert_eq!(
        filogenia::distancias_entre_hojas(&arbol, &hojas),
        matriz_abcd()
    );
}

#[test]
fn acepta_cada_arista_escrita_en_los_dos_sentidos() {
    let arbol =
        filogenia::leer_aristas("a->x:2\r\nx->a:2\nb->x:3\nx->b:3\n\n  x->c:1.5 \nc->x:1.5")
            .unwrap();
    assert_eq!(arbol.aristas().len(), 3);
    assert_eq!(arbol.longitud_entre(0, 1), Some(2.0));
    let sin_longitud = filogenia::leer_aristas("a->b\nb->c").unwrap();
    assert_eq!(sin_longitud.longitud_entre(0, 1), Some(1.0));
}

#[test]
fn rechaza_aristas_que_no_forman_un_arbol() {
    assert!(filogenia::leer_aristas("a->b:1\nb->c:1\nc->a:1").is_err());
    assert!(filogenia::leer_aristas("a->b:1\nc->d:1").is_err());
    assert!(filogenia::leer_aristas("a->b:1\nb->a:2").is_err());
    assert!(filogenia::leer_aristas("a->a:1").is_err());
    assert!(filogenia::leer_aristas("a->b:-1").is_err());
    assert!(filogenia::leer_aristas("a->b:x").is_err());
    assert!(filogenia::leer_aristas("a-b").is_err());
    assert!(filogenia::leer_aristas("a b->c").is_err());
    assert!(filogenia::leer_aristas("").is_err());
    assert!(filogenia::leer_aristas("  \n ").is_err());
}

#[test]
fn longitud_de_rama_de_cada_hoja() {
    let distancias = matriz_abcd();
    let ramas: Vec<f64> = (0..4)
        .map(|hoja| filogenia::longitud_de_rama(&distancias, hoja).unwrap())
        .collect();
    assert_eq!(ramas, [2.0, 3.0, 1.0, 5.0]);
    assert!(filogenia::longitud_de_rama(&distancias, 4).is_err());
    let dos = vec![vec![0.0, 7.0], vec![7.0, 0.0]];
    assert_eq!(filogenia::longitud_de_rama(&dos, 1).unwrap(), 7.0);
    assert!(filogenia::longitud_de_rama(&[vec![0.0]], 0).is_err());
    assert!(filogenia::longitud_de_rama(&[], 0).is_err());
}

#[test]
fn valida_la_matriz_de_distancias() {
    assert!(filogenia::validar_distancias(&matriz_abcd()).is_ok());
    assert!(filogenia::validar_distancias(&[]).is_err());
    assert!(filogenia::validar_distancias(&[vec![0.0, 1.0]]).is_err());
    assert!(filogenia::validar_distancias(&[vec![0.0, 1.0], vec![2.0, 0.0]]).is_err());
    assert!(filogenia::validar_distancias(&[vec![1.0, 1.0], vec![1.0, 0.0]]).is_err());
    assert!(filogenia::validar_distancias(&[vec![0.0, -1.0], vec![-1.0, 0.0]]).is_err());
    assert!(filogenia::validar_distancias(&[vec![0.0, f64::NAN], vec![f64::NAN, 0.0]]).is_err());
}

#[test]
fn condicion_de_los_cuatro_puntos() {
    let distancias = matriz_abcd();
    assert!(filogenia::es_aditiva(&distancias));
    assert_eq!(
        filogenia::sumas_de_cuarteto(&distancias, 0, 1, 2, 3),
        (11.0, 19.0, 19.0)
    );
    let mut alterada = distancias.clone();
    alterada[0][3] = 12.0;
    alterada[3][0] = 12.0;
    assert_eq!(
        filogenia::cuarteto_no_aditivo(&alterada),
        Some((0, 1, 2, 3))
    );
    assert!(!filogenia::es_aditiva(&alterada));
    let sin_triangulo = vec![
        vec![0.0, 1.0, 5.0],
        vec![1.0, 0.0, 1.0],
        vec![5.0, 1.0, 0.0],
    ];
    assert_eq!(filogenia::trio_no_metrico(&sin_triangulo), Some((0, 1, 2)));
    assert!(!filogenia::es_aditiva(&sin_triangulo));
}

#[test]
fn filogenia_aditiva_reconstruye_el_arbol() {
    let distancias = matriz_abcd();
    let arbol = filogenia::filogenia_aditiva(&distancias, &nombres(&["a", "b", "c", "d"])).unwrap();
    assert_eq!(arbol.nombres.len(), 6);
    assert_eq!(arbol.aristas().len(), 5);
    assert_eq!(arbol.nombres[4], "n4");
    assert_eq!(arbol.longitud_entre(0, 4), Some(2.0));
    assert_eq!(arbol.longitud_entre(4, 5), Some(4.0));
    assert!(hermanas(&arbol, "a", "b"));
    assert!(hermanas(&arbol, "c", "d"));
    assert_eq!(
        filogenia::distancias_entre_hojas(&arbol, &[0, 1, 2, 3]),
        distancias
    );
}

#[test]
fn filogenia_aditiva_con_pocas_hojas_y_entradas_malas() {
    let una = filogenia::filogenia_aditiva(&[vec![0.0]], &nombres(&["a"])).unwrap();
    assert_eq!(una.nombres.len(), 1);
    let dos =
        filogenia::filogenia_aditiva(&[vec![0.0, 3.0], vec![3.0, 0.0]], &nombres(&["a", "b"]))
            .unwrap();
    assert_eq!(dos.longitud_entre(0, 1), Some(3.0));
    let mut alterada = matriz_abcd();
    alterada[0][3] = 12.0;
    alterada[3][0] = 12.0;
    let cuatro = nombres(&["a", "b", "c", "d"]);
    assert!(filogenia::filogenia_aditiva(&alterada, &cuatro).is_err());
    assert!(filogenia::filogenia_aditiva(&matriz_abcd(), &nombres(&["a"])).is_err());
    assert!(filogenia::filogenia_aditiva(&[], &[]).is_err());
}

#[test]
fn filogenia_aditiva_con_ramas_de_longitud_cero() {
    let distancias = vec![
        vec![0.0, 2.0, 3.0],
        vec![2.0, 0.0, 5.0],
        vec![3.0, 5.0, 0.0],
    ];
    let arbol = filogenia::filogenia_aditiva(&distancias, &nombres(&["a", "b", "c"])).unwrap();
    assert!(parecidas(
        &filogenia::distancias_entre_hojas(&arbol, &[0, 1, 2]),
        &distancias
    ));
}

#[test]
fn arboles_al_azar_dan_matrices_que_se_reconstruyen() {
    for semilla in 0..30 {
        let cantidad = 3 + (semilla as usize % 12);
        let arbol = arbol_con_pesos(cantidad, semilla);
        let hojas: Vec<usize> = (0..cantidad).collect();
        let distancias = filogenia::distancias_entre_hojas(&arbol, &hojas);
        assert!(filogenia::es_aditiva(&distancias));
        let nombres_de_hojas = numerados(cantidad);
        let aditivo = filogenia::filogenia_aditiva(&distancias, &nombres_de_hojas).unwrap();
        assert_eq!(aditivo.nombres.len(), 2 * cantidad - 2);
        assert!(parecidas(
            &filogenia::distancias_entre_hojas(&aditivo, &hojas),
            &distancias
        ));
        let vecinos = filogenia::union_de_vecinos(&distancias, &nombres_de_hojas).unwrap();
        assert_eq!(vecinos.nombres.len(), 2 * cantidad - 2);
        assert!(parecidas(
            &filogenia::distancias_entre_hojas(&vecinos, &hojas),
            &distancias
        ));
    }
}

#[test]
fn upgma_une_por_promedios_y_mide_alturas() {
    let distancias = vec![
        vec![0.0, 2.0, 4.0, 12.0],
        vec![2.0, 0.0, 4.0, 12.0],
        vec![4.0, 4.0, 0.0, 6.0],
        vec![12.0, 12.0, 6.0, 0.0],
    ];
    let resultado = filogenia::upgma(&distancias, &nombres(&["a", "b", "c", "d"])).unwrap();
    let arbol = &resultado.arbol;
    assert_eq!(resultado.raiz, 6);
    assert_eq!(arbol.longitud_entre(0, 4), Some(1.0));
    assert_eq!(arbol.longitud_entre(1, 4), Some(1.0));
    assert_eq!(arbol.longitud_entre(2, 5), Some(2.0));
    assert_eq!(arbol.longitud_entre(4, 5), Some(1.0));
    assert_eq!(arbol.longitud_entre(3, 6), Some(5.0));
    assert_eq!(arbol.longitud_entre(5, 6), Some(3.0));
    assert_eq!(
        filogenia::trio_no_ultrametrico(&distancias),
        Some((0, 2, 3))
    );
}

#[test]
fn upgma_reproduce_una_matriz_ultrametrica() {
    let distancias = vec![
        vec![0.0, 2.0, 6.0, 6.0],
        vec![2.0, 0.0, 6.0, 6.0],
        vec![6.0, 6.0, 0.0, 4.0],
        vec![6.0, 6.0, 4.0, 0.0],
    ];
    assert_eq!(filogenia::trio_no_ultrametrico(&distancias), None);
    let resultado = filogenia::upgma(&distancias, &nombres(&["a", "b", "c", "d"])).unwrap();
    assert_eq!(
        filogenia::distancias_entre_hojas(&resultado.arbol, &[0, 1, 2, 3]),
        distancias
    );
    let una = filogenia::upgma(&[vec![0.0]], &nombres(&["a"])).unwrap();
    assert_eq!(una.raiz, 0);
    assert!(filogenia::upgma(&[], &[]).is_err());
}

#[test]
fn union_de_vecinos_recupera_un_arbol_aditivo() {
    let arbol =
        filogenia::union_de_vecinos(&matriz_abcd(), &nombres(&["a", "b", "c", "d"])).unwrap();
    assert_eq!(arbol.longitud_entre(0, 4), Some(2.0));
    assert_eq!(arbol.longitud_entre(1, 4), Some(3.0));
    assert_eq!(arbol.longitud_entre(2, 5), Some(1.0));
    assert_eq!(arbol.longitud_entre(3, 5), Some(5.0));
    assert_eq!(arbol.longitud_entre(4, 5), Some(4.0));
}

#[test]
fn union_de_vecinos_paso_a_paso() {
    let mut estado =
        filogenia::empezar_union(&matriz_abcd(), &nombres(&["a", "b", "c", "d"])).unwrap();
    assert_eq!(
        filogenia::distancias_totales(&estado.distancias),
        [23.0, 25.0, 21.0, 29.0]
    );
    let union = filogenia::matriz_de_union(&estado.distancias);
    assert_eq!(union[0][1], -38.0);
    assert_eq!(union[2][3], -38.0);
    assert_eq!(union[0][2], -30.0);
    assert_eq!(union[1][1], 0.0);
    let primero = filogenia::unir_par_de_vecinos(&mut estado).unwrap();
    assert_eq!((primero.primero, primero.segundo, primero.nuevo), (0, 1, 4));
    assert_eq!((primero.rama_primero, primero.rama_segundo), (2.0, 3.0));
    assert_eq!(estado.nodos, [2, 3, 4]);
    assert_eq!(
        estado.distancias,
        vec![
            vec![0.0, 6.0, 5.0],
            vec![6.0, 0.0, 9.0],
            vec![5.0, 9.0, 0.0]
        ]
    );
    let segundo = filogenia::unir_par_de_vecinos(&mut estado).unwrap();
    assert_eq!((segundo.primero, segundo.segundo, segundo.nuevo), (2, 3, 5));
    assert_eq!((segundo.rama_primero, segundo.rama_segundo), (1.0, 5.0));
    assert!(filogenia::unir_par_de_vecinos(&mut estado).is_none());
    assert_eq!(estado.nodos, [4, 5]);
    assert_eq!(filogenia::unir_ultimo_par(&mut estado), Some((4, 5, 4.0)));
    assert_eq!(estado.arbol.longitud_entre(4, 5), Some(4.0));
}

#[test]
fn union_de_vecinos_con_pocas_hojas_y_entradas_malas() {
    let una = filogenia::union_de_vecinos(&[vec![0.0]], &nombres(&["a"])).unwrap();
    assert_eq!(una.nombres.len(), 1);
    let dos = filogenia::union_de_vecinos(&[vec![0.0, 4.0], vec![4.0, 0.0]], &nombres(&["a", "b"]))
        .unwrap();
    assert_eq!(dos.longitud_entre(0, 1), Some(4.0));
    assert!(filogenia::union_de_vecinos(&[], &[]).is_err());
    assert!(filogenia::union_de_vecinos(&[vec![0.0, 1.0]], &nombres(&["a"])).is_err());
    assert!(
        filogenia::union_de_vecinos(&[vec![0.0, 1.0], vec![2.0, 0.0]], &nombres(&["a", "b"]))
            .is_err()
    );
}

#[test]
fn distancia_p_y_jukes_cantor() {
    assert_eq!(
        filogenia::contar_diferencias(b"ACGTACGTAC", b"ACGTTCGTAA"),
        (10, 2)
    );
    assert_eq!(
        filogenia::distancia_p(b"ACGTACGTAC", b"ACGTTCGTAA"),
        Some(0.2)
    );
    assert_eq!(filogenia::distancia_p(b"AC-T", b"ACGA"), Some(1.0 / 3.0));
    assert_eq!(filogenia::distancia_p(b"--", b"AC"), None);
    assert!((filogenia::distancia_jukes_cantor(0.2) - 0.232616).abs() < 1e-6);
    assert_eq!(filogenia::distancia_jukes_cantor(0.0), 0.0);
    assert!(filogenia::distancia_jukes_cantor(0.75).is_infinite());
    assert!(filogenia::distancia_jukes_cantor(1.0).is_infinite());
}

#[test]
fn matrices_de_distancias_desde_secuencias() {
    let secuencias = vec![
        b"ACGTACGTAC".to_vec(),
        b"ACGTTCGTAA".to_vec(),
        b"TGCATGCATG".to_vec(),
    ];
    let p = filogenia::matriz_de_distancias_p(&secuencias).unwrap();
    assert_eq!(p[0][1], 0.2);
    assert_eq!(p[1][0], 0.2);
    assert_eq!(p[0][0], 0.0);
    assert_eq!(p[0][2], 1.0);
    let jukes_cantor = filogenia::matriz_jukes_cantor(&secuencias).unwrap();
    assert!((jukes_cantor[0][1] - 0.232616).abs() < 1e-6);
    assert!(jukes_cantor[0][2].is_infinite());
    assert!(filogenia::matriz_de_distancias_p(&[]).is_err());
    assert!(filogenia::matriz_de_distancias_p(&[b"ACG".to_vec(), b"AC".to_vec()]).is_err());
    assert!(filogenia::matriz_de_distancias_p(&[b"AC".to_vec(), b"--".to_vec()]).is_err());
    assert!(filogenia::matriz_de_distancias_p(&[Vec::new(), Vec::new()]).is_err());
}

#[test]
fn arbol_desde_secuencias_agrupa_a_las_parecidas() {
    let secuencias = vec![
        b"AAAAAAAAAACCCCCCCCCC".to_vec(),
        b"AAAAAAAAATCCCCCCCCCC".to_vec(),
        b"AAAAAGGGGGCCCCCTTTTT".to_vec(),
        b"AAAAAGGGGGCCCCATTTTT".to_vec(),
    ];
    let arbol =
        filogenia::arbol_desde_secuencias(&secuencias, &nombres(&["a", "b", "c", "d"])).unwrap();
    assert!(hermanas(&arbol, "a", "b"));
    assert!(hermanas(&arbol, "c", "d"));
    let saturadas = vec![b"ACGT".to_vec(), b"CATG".to_vec(), b"ACGA".to_vec()];
    assert!(filogenia::arbol_desde_secuencias(&saturadas, &nombres(&["a", "b", "c"])).is_err());
}

#[test]
fn newick_ida_y_vuelta() {
    let leido = filogenia::leer_newick("((a:1,b:2.5):3,c:4);").unwrap();
    assert_eq!(leido.raiz, 0);
    assert_eq!(leido.arbol.nombres.len(), 5);
    assert_eq!(
        filogenia::escribir_newick(&leido.arbol, leido.raiz),
        "((a:1,b:2.5):3,c:4);"
    );
    let con_nombres = filogenia::leer_newick(" ( (a , b)ab ,\n c ) raiz ; ").unwrap();
    assert_eq!(con_nombres.arbol.nombres[con_nombres.raiz], "raiz");
    assert!(con_nombres.arbol.buscar_nodo("ab").is_some());
    assert_eq!(
        filogenia::escribir_newick(&con_nombres.arbol, con_nombres.raiz),
        "((a:1,b:1)ab:1,c:1)raiz;"
    );
    let sola = filogenia::leer_newick("a;").unwrap();
    assert_eq!(sola.arbol.nombres, ["a"]);
    assert_eq!(filogenia::escribir_newick(&sola.arbol, 0), "a;");
}

#[test]
fn newick_mal_escrito() {
    for malo in [
        "",
        ";",
        "(a,b",
        "a,b",
        "(a,,b);",
        "(a,b));",
        "(a,a);",
        "(a:x,b);",
        "(a:-1,b);",
        "((a,b)c d);",
        "(a,b);c",
        "a(b,c);",
        "(a,b)(c);",
        "();",
        "(a,);",
        "(a b,c);",
        "(a,b):;",
        "(@ecoli,b);",
    ] {
        assert!(filogenia::leer_newick(malo).is_err(), "aceptó «{malo}»");
    }
}

#[test]
fn quitar_la_raiz_une_sus_dos_ramas() {
    let leido = filogenia::leer_newick("((a:1,b:1):2,(c:1,d:1):3);").unwrap();
    let arbol = filogenia::quitar_raiz(&leido);
    assert_eq!(arbol.nombres.len(), 6);
    assert_eq!(arbol.aristas().len(), 5);
    assert!(hermanas(&arbol, "a", "b"));
    let a = arbol.buscar_nodo("a").unwrap();
    let c = arbol.buscar_nodo("c").unwrap();
    assert_eq!(filogenia::distancias_desde(&arbol, a)[c], 7.0);
    let ya_sin_raiz = filogenia::leer_newick("(a,b,(c,d));").unwrap();
    assert_eq!(filogenia::quitar_raiz(&ya_sin_raiz), ya_sin_raiz.arbol);
}

fn secuencias_por_nodo(arbol: &Arbol, hojas: &[(&str, &str)]) -> Vec<Vec<u8>> {
    let mut secuencias = vec![Vec::new(); arbol.nombres.len()];
    for (nombre, secuencia) in hojas {
        secuencias[arbol.buscar_nodo(nombre).unwrap()] = secuencia.as_bytes().to_vec();
    }
    secuencias
}

#[test]
fn parsimonia_pequena_en_un_arbol_con_raiz() {
    let leido = filogenia::leer_newick("((a,b),(c,d));").unwrap();
    let secuencias = secuencias_por_nodo(
        &leido.arbol,
        &[("a", "ACG"), ("b", "ATG"), ("c", "TCA"), ("d", "TCG")],
    );
    let resultado = filogenia::parsimonia_pequena(&leido.arbol, leido.raiz, &secuencias).unwrap();
    assert_eq!(resultado.puntaje, 3);
    assert_eq!(resultado.cambios_por_posicion, [1, 1, 1]);
    assert_eq!(resultado.secuencias[0], b"ACG");
    assert_eq!(resultado.secuencias[1], b"ACG");
    assert_eq!(resultado.secuencias[4], b"TCG");
    assert_eq!(resultado.secuencias[2], b"ACG");
    assert_eq!(resultado.arbol.longitud_entre(0, 4), Some(1.0));
    assert_eq!(resultado.arbol.longitud_entre(1, 3), Some(1.0));
    assert_eq!(resultado.arbol.longitud_entre(0, 1), Some(0.0));
    let total: f64 = resultado
        .arbol
        .aristas()
        .iter()
        .map(|arista| arista.2)
        .sum();
    assert_eq!(total, 3.0);
}

#[test]
fn parsimonia_pequena_rechaza_datos_incompletos() {
    let leido = filogenia::leer_newick("((a,b),(c,d));").unwrap();
    let arbol = &leido.arbol;
    let sin_d = secuencias_por_nodo(arbol, &[("a", "ACG"), ("b", "ATG"), ("c", "TCA")]);
    assert!(filogenia::parsimonia_pequena(arbol, 0, &sin_d).is_err());
    let distinto_largo = secuencias_por_nodo(
        arbol,
        &[("a", "ACG"), ("b", "AT"), ("c", "TCA"), ("d", "TCG")],
    );
    assert!(filogenia::parsimonia_pequena(arbol, 0, &distinto_largo).is_err());
    let con_hueco = secuencias_por_nodo(
        arbol,
        &[("a", "AC-"), ("b", "ATG"), ("c", "TCA"), ("d", "TCG")],
    );
    assert!(filogenia::parsimonia_pequena(arbol, 0, &con_hueco).is_err());
    assert!(filogenia::parsimonia_pequena(arbol, 9, &sin_d).is_err());
    assert!(filogenia::parsimonia_pequena(arbol, 0, &[]).is_err());
}

#[test]
fn parsimonia_sin_raiz_no_depende_de_donde_se_empiece() {
    let leido = filogenia::leer_newick("(a,b,(c,d));").unwrap();
    let secuencias = secuencias_por_nodo(
        &leido.arbol,
        &[("a", "ACG"), ("b", "ATG"), ("c", "TCA"), ("d", "TCG")],
    );
    let resultado = filogenia::parsimonia_sin_raiz(&leido.arbol, &secuencias).unwrap();
    assert_eq!(resultado.puntaje, 3);
    for raiz in 0..leido.arbol.nombres.len() {
        let otra = filogenia::parsimonia_pequena(&leido.arbol, raiz, &secuencias).unwrap();
        assert_eq!(otra.puntaje, 3);
    }
    let dos = filogenia::leer_newick("(a,b);").unwrap();
    let pareja = secuencias_por_nodo(&dos.arbol, &[("a", "AAC"), ("b", "ATG")]);
    let sin_raiz = filogenia::quitar_raiz(&dos);
    let pareja_sin_raiz = secuencias_por_nodo(&sin_raiz, &[("a", "AAC"), ("b", "ATG")]);
    assert_eq!(
        filogenia::parsimonia_pequena(&dos.arbol, dos.raiz, &pareja)
            .unwrap()
            .puntaje,
        2
    );
    assert_eq!(
        filogenia::parsimonia_sin_raiz(&sin_raiz, &pareja_sin_raiz)
            .unwrap()
            .puntaje,
        2
    );
}

#[test]
fn parsimonia_sin_raiz_exige_secuencia_en_cada_hoja() {
    let leido = filogenia::leer_newick("((a,b),(c,d));").unwrap();
    let arbol = filogenia::quitar_raiz(&leido);
    let sin_d = secuencias_por_nodo(&arbol, &[("a", "ACG"), ("b", "ATG"), ("c", "TCA")]);
    let error = filogenia::parsimonia_sin_raiz(&arbol, &sin_d).unwrap_err();
    assert!(error.contains("«d»"), "{error}");
    let completas = secuencias_por_nodo(
        &arbol,
        &[("a", "ACG"), ("b", "ATG"), ("c", "TCA"), ("d", "TCG")],
    );
    let resultado = filogenia::parsimonia_sin_raiz(&arbol, &completas).unwrap();
    assert_eq!(resultado.puntaje, 3);
    let d = arbol.buscar_nodo("d").unwrap();
    assert_eq!(resultado.secuencias[d], b"TCG");
    let pareja = filogenia::leer_aristas("a->b").unwrap();
    let solo_a = secuencias_por_nodo(&pareja, &[("a", "ACG")]);
    assert!(filogenia::parsimonia_sin_raiz(&pareja, &solo_a).is_err());
}

#[test]
fn vecinos_cercanos_intercambian_subarboles() {
    let arbol = filogenia::leer_aristas("a->u\nb->u\nu->v\nc->v\nd->v").unwrap();
    let u = arbol.buscar_nodo("u").unwrap();
    let v = arbol.buscar_nodo("v").unwrap();
    let (primero, segundo) = filogenia::vecinos_cercanos(&arbol, u, v).unwrap();
    assert!(hermanas(&primero, "a", "c"));
    assert!(hermanas(&primero, "b", "d"));
    assert!(hermanas(&segundo, "a", "d"));
    assert!(hermanas(&segundo, "b", "c"));
    assert_eq!(primero.aristas().len(), 5);
    let a = arbol.buscar_nodo("a").unwrap();
    assert!(filogenia::vecinos_cercanos(&arbol, a, u).is_err());
    assert!(filogenia::vecinos_cercanos(&arbol, u, u).is_err());
    assert!(filogenia::vecinos_cercanos(&arbol, u, 99).is_err());
}

#[test]
fn vecinos_cercanos_conservan_longitudes() {
    let arbol = filogenia::leer_aristas("a->u:1\nb->u:2\nu->v:3\nc->v:4\nd->v:5").unwrap();
    let u = arbol.buscar_nodo("u").unwrap();
    let v = arbol.buscar_nodo("v").unwrap();
    let (primero, _) = filogenia::vecinos_cercanos(&arbol, u, v).unwrap();
    let b = primero.buscar_nodo("b").unwrap();
    let c = primero.buscar_nodo("c").unwrap();
    assert_eq!(primero.longitud_entre(b, v), Some(2.0));
    assert_eq!(primero.longitud_entre(c, u), Some(4.0));
}

#[test]
fn parsimonia_grande_mejora_con_intercambios() {
    let arbol = filogenia::leer_aristas("a->u\nc->u\nu->v\nb->v\nd->v").unwrap();
    let secuencias = secuencias_por_nodo(
        &arbol,
        &[("a", "AAAA"), ("b", "AAAA"), ("c", "CCCC"), ("d", "CCCC")],
    );
    let busqueda = filogenia::parsimonia_grande(&arbol, &secuencias, 10).unwrap();
    assert_eq!(busqueda.puntaje_inicial, 8);
    assert_eq!(busqueda.movimientos.len(), 1);
    assert_eq!(busqueda.movimientos[0].puntaje, 4);
    assert_eq!(busqueda.parsimonia.puntaje, 4);
    assert!(!busqueda.limite_alcanzado);
    assert!(hermanas(&busqueda.arbol, "a", "b"));
    let sin_movimientos = filogenia::parsimonia_grande(&arbol, &secuencias, 0).unwrap();
    assert_eq!(sin_movimientos.parsimonia.puntaje, 8);
    assert!(sin_movimientos.limite_alcanzado);
}

#[test]
fn parsimonia_grande_rechaza_arboles_no_binarios() {
    let estrella = filogenia::leer_aristas("a->u\nb->u\nc->u\nd->u").unwrap();
    let secuencias =
        secuencias_por_nodo(&estrella, &[("a", "A"), ("b", "A"), ("c", "C"), ("d", "C")]);
    assert!(filogenia::parsimonia_grande(&estrella, &secuencias, 5).is_err());
    let pareja = filogenia::leer_aristas("a->b").unwrap();
    let dos = secuencias_por_nodo(&pareja, &[("a", "A"), ("b", "C")]);
    assert!(filogenia::parsimonia_grande(&pareja, &dos, 5).is_err());
}

#[test]
fn fitch_y_sankoff_dan_el_mismo_puntaje() {
    for semilla in 0..20 {
        let cantidad = 3 + (semilla as usize % 10);
        let arbol = filogenia::arbol_al_azar(&numerados(cantidad), semilla).unwrap();
        let mut secuencias = secuencias_al_azar(cantidad, 25, semilla + 7);
        secuencias.resize(arbol.nombres.len(), Vec::new());
        let sankoff = filogenia::parsimonia_sin_raiz(&arbol, &secuencias).unwrap();
        let busqueda = filogenia::parsimonia_grande(&arbol, &secuencias, 0).unwrap();
        assert_eq!(busqueda.puntaje_inicial, sankoff.puntaje);
        let mejorada = filogenia::parsimonia_grande(&arbol, &secuencias, 50).unwrap();
        assert!(mejorada.parsimonia.puntaje <= sankoff.puntaje);
        let ultimo = mejorada
            .movimientos
            .last()
            .map(|movimiento| movimiento.puntaje)
            .unwrap_or(mejorada.puntaje_inicial);
        assert_eq!(ultimo, mejorada.parsimonia.puntaje);
        let suma: f64 = mejorada
            .parsimonia
            .arbol
            .aristas()
            .iter()
            .map(|arista| arista.2)
            .sum();
        assert_eq!(suma as usize, mejorada.parsimonia.puntaje);
    }
}

fn todos_los_vecinos_cercanos(arbol: &Arbol) -> Vec<Arbol> {
    let mut vecinos = Vec::new();
    for (a, b, _) in arbol.aristas() {
        if let Ok((primero, segundo)) = filogenia::vecinos_cercanos(arbol, a, b) {
            vecinos.push(primero);
            vecinos.push(segundo);
        }
    }
    vecinos
}

fn puntaje_de(arbol: &Arbol, secuencias: &[Vec<u8>]) -> usize {
    filogenia::parsimonia_sin_raiz(arbol, secuencias)
        .unwrap()
        .puntaje
}

#[test]
fn cada_movimiento_es_el_mejor_vecino_y_su_puntaje_es_exacto() {
    for semilla in 0..25 {
        let cantidad = 5 + (semilla as usize % 9);
        let inicial = filogenia::arbol_al_azar(&numerados(cantidad), semilla + 300).unwrap();
        let mut secuencias = secuencias_al_azar(cantidad, 30, semilla + 900);
        secuencias.resize(inicial.nombres.len(), Vec::new());
        let busqueda = filogenia::parsimonia_grande(&inicial, &secuencias, 30).unwrap();
        assert_eq!(busqueda.puntaje_inicial, puntaje_de(&inicial, &secuencias));
        let mut arbol = inicial.clone();
        let mut puntaje = busqueda.puntaje_inicial;
        for movimiento in &busqueda.movimientos {
            let mejor_vecino = todos_los_vecinos_cercanos(&arbol)
                .iter()
                .map(|vecino| puntaje_de(vecino, &secuencias))
                .min()
                .unwrap();
            assert_eq!(movimiento.puntaje, mejor_vecino);
            assert!(movimiento.puntaje < puntaje);
            let (a, b) = movimiento.arista;
            let (x, y) = movimiento.intercambio;
            filogenia::intercambiar(&mut arbol, a, x, b, y);
            assert_eq!(puntaje_de(&arbol, &secuencias), movimiento.puntaje);
            puntaje = movimiento.puntaje;
        }
        if !busqueda.limite_alcanzado {
            assert!(
                todos_los_vecinos_cercanos(&arbol)
                    .iter()
                    .all(|vecino| puntaje_de(vecino, &secuencias) >= puntaje)
            );
        }
        assert_eq!(busqueda.arbol, arbol);
        assert_eq!(busqueda.parsimonia.puntaje, puntaje);
    }
}

#[test]
fn arbol_al_azar_es_binario_y_reproducible() {
    let primero = filogenia::arbol_al_azar(&numerados(9), 42).unwrap();
    let segundo = filogenia::arbol_al_azar(&numerados(9), 42).unwrap();
    assert_eq!(primero, segundo);
    assert_eq!(primero.nombres.len(), 16);
    assert_eq!(primero.hojas(), (0..9).collect::<Vec<usize>>());
    assert!((9..16).all(|nodo| primero.grado(nodo) == 3));
    assert!(filogenia::arbol_al_azar(&numerados(2), 1).is_err());
}

#[test]
fn aditividad_por_reconstruccion_coincide_con_los_cuatro_puntos() {
    let mut azar = Azar::nuevo(77);
    let mut aditivas = 0;
    for semilla in 0..60 {
        let cantidad = 3 + (semilla as usize % 8);
        let arbol = arbol_con_pesos(cantidad, semilla + 500);
        let hojas: Vec<usize> = (0..cantidad).collect();
        let mut distancias = filogenia::distancias_entre_hojas(&arbol, &hojas);
        if semilla % 2 == 1 {
            let i = azar.entero(cantidad);
            let j = (i + 1 + azar.entero(cantidad - 1)) % cantidad;
            let cambio = 1.0 + azar.entero(3) as f64;
            distancias[i][j] += cambio;
            distancias[j][i] += cambio;
        }
        let por_definicion = filogenia::trio_no_metrico(&distancias).is_none()
            && filogenia::cuarteto_no_aditivo(&distancias).is_none();
        assert_eq!(filogenia::es_aditiva(&distancias), por_definicion);
        if por_definicion {
            aditivas += 1;
        }
    }
    assert!(aditivas > 30 && aditivas < 60);
}

#[test]
fn distancias_desde_termina_aunque_la_suma_desborde() {
    let mut arbol = Arbol::default();
    let nodos: Vec<usize> = (0..4)
        .map(|n| arbol.agregar_nodo(&format!("h{n}")))
        .collect();
    for par in nodos.windows(2) {
        arbol.conectar(par[0], par[1], 1e308);
    }
    let distancias = filogenia::distancias_desde(&arbol, 0);
    assert_eq!(distancias.len(), 4);
    assert!(distancias[3].is_infinite());
}

#[test]
fn distancias_y_longitudes_gigantes_se_rechazan() {
    let enorme = vec![
        vec![0.0, 1e308, 1e308],
        vec![1e308, 0.0, 1e308],
        vec![1e308, 1e308, 0.0],
    ];
    let hojas = nombres(&["a", "b", "c"]);
    assert!(filogenia::upgma(&enorme, &hojas).is_err());
    assert!(filogenia::union_de_vecinos(&enorme, &hojas).is_err());
    assert!(filogenia::filogenia_aditiva(&enorme, &hojas).is_err());
    assert!(filogenia::leer_aristas("a->b:1e308").is_err());
    assert!(filogenia::leer_newick("(a:1e13,b:1);").is_err());
    assert!(filogenia::leer_aristas("a->b:1000000000000").is_ok());
}

#[test]
fn la_longitud_de_rama_es_el_minimo_de_todas_las_parejas() {
    let distancias = vec![
        vec![0.0, 5.0, 5.0, 5.0],
        vec![5.0, 0.0, 10.0, 1.0],
        vec![5.0, 10.0, 0.0, 1.0],
        vec![5.0, 1.0, 1.0, 0.0],
    ];
    let (rama, primera, segunda) = filogenia::longitud_de_rama_con_par(&distancias, 3).unwrap();
    assert_eq!(rama, -4.0);
    assert_eq!((primera, segunda), (1, 2));
    assert_eq!(filogenia::longitud_de_rama(&distancias, 3).unwrap(), -4.0);
    let dos = vec![vec![0.0, 7.0], vec![7.0, 0.0]];
    assert_eq!(
        filogenia::longitud_de_rama_con_par(&dos, 0).unwrap(),
        (7.0, 1, 1)
    );
}

#[test]
fn jukes_cantor_de_cero_es_cero_positivo() {
    let cero = filogenia::distancia_jukes_cantor(0.0);
    assert_eq!(cero, 0.0);
    assert!(cero.is_sign_positive());
    let matriz = filogenia::matriz_jukes_cantor(&[b"ACGT".to_vec(), b"ACGT".to_vec()]).unwrap();
    assert!(
        matriz
            .iter()
            .flatten()
            .all(|valor| valor.is_sign_positive())
    );
}

#[test]
fn newick_rechaza_dos_longitudes_en_un_nodo() {
    assert!(filogenia::leer_newick("(a:1:2,b:1);").is_err());
    assert!(filogenia::leer_newick("(a:1,b:1)r:3:4;").is_err());
    let leido = filogenia::leer_newick("(a,b:2);").unwrap();
    assert_eq!(leido.arbol.longitud_entre(0, 1), Some(1.0));
}
