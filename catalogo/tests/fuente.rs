use catalogo::{codigo, formato};

#[test]
fn muestra_la_funcion_y_lo_que_usa() {
    let fragmentos = codigo("replicacion.contenido_gc").unwrap();
    assert_eq!(fragmentos[0].nombre, "secuencia::contenido_gc");
    assert!(fragmentos[0].fuente.starts_with("pub fn contenido_gc"));
    assert!(fragmentos[0].fuente.trim_end().ends_with('}'));
}

#[test]
fn formato_de_numeros() {
    assert_eq!(formato::miles(4641652), "4\u{a0}641\u{a0}652");
    assert_eq!(formato::miles(999), "999");
    assert_eq!(formato::porcentaje(0.5079, 1), "50.8\u{a0}%");
    assert_eq!(formato::cientifico(3.2e-45), "3.2e-45");
    let puntos: Vec<[f64; 2]> = (0..10_000).map(|i| [i as f64, (i % 7) as f64]).collect();
    let reducidos = formato::reducir_puntos(&puntos, 100);
    assert!(reducidos.len() <= 100);
    assert!(reducidos.windows(2).all(|par| par[0][0] < par[1][0]));
}

#[test]
fn los_tipos_de_arreglo_no_cortan_el_codigo() {
    let fragmentos = catalogo::fragmentos("secuencia::traducir").unwrap();
    let nombres: Vec<&str> = fragmentos.iter().map(|f| f.nombre.as_str()).collect();
    assert!(nombres.contains(&"secuencia::aminoacido"));
    let tabla = fragmentos
        .iter()
        .find(|f| f.nombre == "secuencia::CODIGO_GENETICO")
        .unwrap();
    assert!(tabla.fuente.contains("KNKNTTTT"));
    assert!(tabla.fuente.trim_end().ends_with(';'));
    let matrices = catalogo::fragmentos("matrices::blosum62").unwrap();
    let nombres: Vec<&str> = matrices.iter().map(|f| f.nombre.as_str()).collect();
    assert!(nombres.contains(&"matrices::BLOSUM62"));
    assert!(nombres.contains(&"matrices::indice_aminoacido"));
    assert!(nombres.contains(&"matrices::construir_indices"));
    let blosum = matrices
        .iter()
        .find(|f| f.nombre == "matrices::BLOSUM62")
        .unwrap();
    assert!(blosum.fuente.trim_end().ends_with("];"));
}

#[test]
fn todos_los_modulos_del_nucleo_tienen_su_fuente() {
    let carpeta = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../nucleo/src");
    for entrada in std::fs::read_dir(carpeta).unwrap() {
        let nombre = entrada.unwrap().file_name().to_string_lossy().to_string();
        let Some(modulo) = nombre.strip_suffix(".rs") else {
            continue;
        };
        if modulo == "lib" {
            continue;
        }
        let error = catalogo::fragmentos(&format!("{modulo}::funcion_que_no_existe")).unwrap_err();
        assert!(
            error.contains("No encuentro el código"),
            "{modulo}: {error}"
        );
    }
}
