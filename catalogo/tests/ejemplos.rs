use catalogo::Entrada;
use nucleo::formatos::{abrir_fastq, leer_sam, leer_vcf};
use std::io::BufRead;

#[test]
fn los_ejemplos_de_archivos_existen_y_se_leen() {
    let entrada = Entrada::nueva(&[
        ("a", "@ejemplo:lecturas"),
        ("b", "@ejemplo:alineamientos"),
        ("c", "@ejemplo:mutaciones"),
        ("d", "@ejemplo:nada"),
    ]);
    let lecturas = entrada.ruta("a").unwrap();
    let cantidad = abrir_fastq(&lecturas)
        .unwrap()
        .filter(|l| l.is_ok())
        .count();
    assert!(cantidad > 35_000);
    let alineamientos = entrada.ruta("b").unwrap();
    let mut registros = 0;
    for linea in nucleo::formatos::abrir(&alineamientos)
        .unwrap()
        .split(b'\n')
    {
        if leer_sam(&linea.unwrap()).unwrap().is_some() {
            registros += 1;
        }
    }
    assert_eq!(registros, cantidad);
    let mutaciones = entrada.ruta("c").unwrap();
    let mut variantes = 0;
    for linea in nucleo::formatos::abrir(&mutaciones).unwrap().split(b'\n') {
        if leer_vcf(&linea.unwrap()).unwrap().is_some() {
            variantes += 1;
        }
    }
    assert_eq!(variantes, 26);
    assert!(entrada.ruta("d").is_err());
}
