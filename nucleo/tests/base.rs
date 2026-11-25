use nucleo::ecoli;
use nucleo::formatos::{self, LectorFastq, Lectura};
use nucleo::matrices::{blosum62, pam250};
use nucleo::secuencia::*;
use std::io::{BufReader, Cursor};

#[test]
fn complemento_inverso_de_adn() {
    assert_eq!(complemento_inverso(b"AAAACCCGGT"), b"ACCGGGTTTT");
    assert_eq!(complemento_inverso(b""), b"");
}

#[test]
fn contenido_gc_cuenta_g_y_c() {
    assert_eq!(contenido_gc(b"GGCCAATT"), 0.5);
    assert_eq!(contenido_gc(b""), 0.0);
}

#[test]
fn traduce_codones_y_paradas() {
    assert_eq!(traducir(b"ATGGCCTAA"), b"MA*");
    assert_eq!(traducir(b"AUGUGGUGA"), b"MW*");
    assert_eq!(aminoacido(b"ANG"), b'X');
    assert_eq!(traducir(b"ATGGC"), b"M");
}

#[test]
fn limpiar_quita_espacios_y_numeros() {
    assert_eq!(limpiar(" acg t\r\n12g"), b"ACGTG");
}

#[test]
fn lee_fasta_con_y_sin_cabecera() {
    let registros = formatos::leer_fasta(b">uno primero\r\nACG\ntt\n>dos\nGGG\n");
    assert_eq!(registros.len(), 2);
    assert_eq!(registros[0].nombre, "uno primero");
    assert_eq!(registros[0].secuencia, b"ACGTT");
    assert_eq!(registros[1].secuencia, b"GGG");
    let sueltos = formatos::leer_fasta(b"acgt\nacgt");
    assert_eq!(sueltos.len(), 1);
    assert_eq!(sueltos[0].secuencia, b"ACGTACGT");
    assert!(formatos::leer_fasta(b"  \n").is_empty());
}

#[test]
fn lee_fastq_y_detecta_errores() {
    let texto = b"@r1\r\nACGT\r\n+\r\nIIII\r\n\n@r2\nac\n+r2\n#I\n".to_vec();
    let lector = LectorFastq::nuevo(Box::new(BufReader::new(Cursor::new(texto))));
    let lecturas: Vec<Lectura> = lector.map(|l| l.unwrap()).collect();
    assert_eq!(lecturas.len(), 2);
    assert_eq!(lecturas[0].nombre, b"r1");
    assert_eq!(lecturas[1].secuencia, b"AC");
    let roto = b"@r1\nACGT\n+\nII\n".to_vec();
    let mut lector = LectorFastq::nuevo(Box::new(BufReader::new(Cursor::new(roto))));
    assert!(lector.next().unwrap().is_err());
    let cortado = b"@r1\nACGT\n".to_vec();
    let mut lector = LectorFastq::nuevo(Box::new(BufReader::new(Cursor::new(cortado))));
    assert!(lector.next().unwrap().is_err());
}

#[test]
fn matrices_conocidas() {
    assert_eq!(blosum62(b'W', b'W'), 11);
    assert_eq!(blosum62(b'A', b'R'), -1);
    assert_eq!(pam250(b'W', b'W'), 17);
    assert_eq!(pam250(b'C', b'C'), 12);
    assert_eq!(blosum62(b'A', b'*'), -4);
}

#[test]
fn genoma_de_ecoli_incluido() {
    let genoma = ecoli::genoma();
    assert_eq!(genoma.len(), 4_641_652);
    assert_eq!(&genoma[..12], b"AGCTTTTCATTC");
    assert!((contenido_gc(genoma) - 0.5079).abs() < 0.001);
}

#[test]
fn genes_de_ecoli_incluidos() {
    let genes = ecoli::genes();
    assert_eq!(genes.len(), 4300);
    let dnaa = ecoli::buscar_gen("dnaA").unwrap();
    assert_eq!(dnaa.hebra, b'-');
    let adn = ecoli::secuencia_gen(dnaa);
    assert_eq!(&adn[..3], b"GTG");
    let proteina = traducir(&adn);
    assert_eq!(&proteina[1..proteina.len() - 1], &dnaa.proteina[1..]);
    let rpob = ecoli::buscar_gen("b3987").unwrap();
    let proteina = traducir(&ecoli::secuencia_gen(rpob));
    assert_eq!(&proteina[..proteina.len() - 1], rpob.proteina.as_slice());
    assert_eq!(*proteina.last().unwrap(), b'*');
    assert_eq!(ecoli::region(1, 4), Some(&b"AGCT"[..]));
    assert_eq!(ecoli::region(0, 4), None);
}

#[test]
fn sam_ida_y_vuelta() {
    use nucleo::formatos::{RegistroSam, escribir_sam, leer_sam, texto_a_cigar};
    let registro = RegistroSam {
        nombre: b"r1".to_vec(),
        bandera: 16,
        referencia: "NC_000913.3".to_string(),
        posicion: 100,
        calidad_mapeo: 60,
        cigar: texto_a_cigar(b"5M1I4M").unwrap(),
        secuencia: b"ACGTACGTAC".to_vec(),
        calidad: b"IIIIIIIIII".to_vec(),
    };
    let mut texto = Vec::new();
    escribir_sam(&mut texto, &registro).unwrap();
    let leido = leer_sam(&texto).unwrap().unwrap();
    assert_eq!(leido, registro);
    assert!(leido.mapeado() && leido.hebra_reversa());
    assert!(leer_sam(b"@SQ\tSN:x").unwrap().is_none());
    assert!(texto_a_cigar(b"5Q").is_err());
    assert!(texto_a_cigar(b"M").is_err());
}

#[test]
fn vcf_ida_y_vuelta() {
    use nucleo::formatos::{Variante, escribir_vcf, leer_vcf};
    let variante = Variante {
        cromosoma: "NC_000913.3".to_string(),
        posicion: 4182813,
        referencia: b"C".to_vec(),
        alternativa: b"T".to_vec(),
        calidad: 99.0,
        profundidad: 20,
        frecuencia: 1.0,
    };
    let mut texto = Vec::new();
    escribir_vcf(&mut texto, &variante).unwrap();
    assert_eq!(leer_vcf(&texto).unwrap().unwrap(), variante);
    assert!(leer_vcf(b"#CHROM").unwrap().is_none());
}

#[test]
fn azar_reproducible() {
    use nucleo::azar::Azar;
    let mut a = Azar::nuevo(7);
    let mut b = Azar::nuevo(7);
    let primeros: Vec<u64> = (0..5).map(|_| a.siguiente()).collect();
    let segundos: Vec<u64> = (0..5).map(|_| b.siguiente()).collect();
    assert_eq!(primeros, segundos);
    let mut c = Azar::nuevo(1);
    for _ in 0..1000 {
        let x = c.decimal();
        assert!((0.0..1.0).contains(&x));
        assert!(c.entero(10) < 10);
    }
    assert_eq!(c.elegir_con_pesos(&[0.0, 1.0, 0.0]), 1);
}
