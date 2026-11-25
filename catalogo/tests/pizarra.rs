use catalogo::pizarra::{Marca, Narracion, Pizarra};

#[test]
fn la_pizarra_arma_el_estado() {
    let genoma: Vec<u8> = (0..1000).map(|i| b"ACGT"[i % 4]).collect();
    let paso = Pizarra::nueva()
        .variable("k", 3)
        .variable("contador", 2)
        .secuencia(
            &genoma,
            &[(500, 503, Marca::Actual), (10, 13, Marca::Hallado)],
        )
        .tabla(
            "cuentas",
            &["patrón", "veces"],
            vec![vec!["ACG".into(), "2".into()]],
            &[(0, 1)],
        )
        .lista("hallados", vec!["10".into()], &[0])
        .paso("let contador = 0;", "Empieza en cero.");
    assert_eq!(paso.estado["tipo"], "pizarra");
    assert_eq!(paso.estado["variables"][1][0], "contador");
    let secuencia = &paso.estado["secuencia"];
    let desplazamiento = secuencia["desplazamiento"].as_u64().unwrap() as usize;
    let texto = secuencia["texto"].as_str().unwrap();
    assert!(texto.len() <= 240);
    assert!(desplazamiento <= 500 && 503 <= desplazamiento + texto.len());
    assert_eq!(secuencia["largo"], 1000);
    assert_eq!(secuencia["marcas"].as_array().unwrap().len(), 1);
}

#[test]
fn la_narracion_respeta_el_maximo() {
    let mut narracion = Narracion::nueva(5);
    for i in 0..10 {
        narracion.agregar(Pizarra::nueva().variable("i", i).paso("x", "paso"));
    }
    let pasos = narracion.terminar(Pizarra::nueva().paso("x", "fin"));
    assert_eq!(pasos.len(), 5);
    assert!(pasos[4].nota.starts_with("Se resumieron 6 pasos"));
}
