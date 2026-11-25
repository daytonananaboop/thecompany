use nucleo::azar::Azar;
use nucleo::ecoli;
use nucleo::indice_fm::*;
use nucleo::replicacion::{posiciones_aproximadas, posiciones_del_patron};
use nucleo::secuencia::complemento_inverso;

fn texto_al_azar(azar: &mut Azar, longitud: usize, letras: &[u8]) -> Vec<u8> {
    (0..longitud)
        .map(|_| letras[azar.entero(letras.len())])
        .collect()
}

fn sufijos_a_mano(texto: &[u8]) -> Vec<u32> {
    let mut con_fin: Vec<u16> = texto.iter().map(|&letra| letra as u16 + 1).collect();
    con_fin.push(0);
    let mut orden: Vec<u32> = (0..con_fin.len() as u32).collect();
    orden.sort_by(|&a, &b| con_fin[a as usize..].cmp(&con_fin[b as usize..]));
    orden
}

fn bwt_a_mano(texto: &[u8]) -> Vec<u8> {
    let mut con_fin = texto.to_vec();
    con_fin.push(b'$');
    let n = con_fin.len();
    let mut rotaciones: Vec<Vec<u8>> = (0..n)
        .map(|i| {
            let mut rotacion = con_fin[i..].to_vec();
            rotacion.extend_from_slice(&con_fin[..i]);
            rotacion
        })
        .collect();
    rotaciones.sort_by(|a, b| {
        let clave = |r: &Vec<u8>| -> Vec<u16> {
            r.iter()
                .map(|&l| if l == b'$' { 0 } else { l as u16 + 1 })
                .collect()
        };
        clave(a).cmp(&clave(b))
    });
    rotaciones.iter().map(|r| r[n - 1]).collect()
}

fn hamming(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).filter(|(x, y)| x != y).count()
}

fn aciertos_a_mano(texto: &[u8], patron: &[u8], d: usize) -> Vec<Acierto> {
    let mut aciertos = Vec::new();
    let inverso = complemento_inverso(patron);
    let mut hebras = vec![(b'+', patron.to_vec())];
    if inverso != patron {
        hebras.push((b'-', inverso));
    }
    for (hebra, buscado) in hebras {
        if buscado.len() > texto.len() {
            continue;
        }
        for (posicion, ventana) in texto.windows(buscado.len()).enumerate() {
            let diferencias = hamming(ventana, &buscado);
            if diferencias <= d {
                aciertos.push(Acierto {
                    posicion,
                    hebra,
                    diferencias,
                });
            }
        }
    }
    aciertos.sort();
    aciertos
}

#[test]
fn arreglo_de_sufijos_a_mano() {
    assert_eq!(
        arreglo_de_sufijos(b"ACATG").unwrap(),
        vec![5, 0, 2, 1, 4, 3]
    );
    assert_eq!(arreglo_de_sufijos(b"").unwrap(), vec![0]);
    assert_eq!(arreglo_de_sufijos(b"T").unwrap(), vec![1, 0]);
    assert_eq!(arreglo_de_sufijos(b"GGGG").unwrap(), vec![4, 3, 2, 1, 0]);
}

#[test]
fn sa_is_igual_que_ordenar_a_mano() {
    let mut azar = Azar::nuevo(7);
    let alfabetos: [&[u8]; 5] = [
        b"AB",
        b"ACGT",
        b"ACGTN",
        b"abcdefghijklmnopqrstuvwxyz",
        b"A",
    ];
    for prueba in 0..600 {
        let letras = alfabetos[prueba % alfabetos.len()];
        let longitud = 1 + azar.entero(300);
        let texto = texto_al_azar(&mut azar, longitud, letras);
        assert_eq!(
            arreglo_de_sufijos(&texto).unwrap(),
            sufijos_a_mano(&texto),
            "texto {}",
            String::from_utf8_lossy(&texto)
        );
    }
}

#[test]
fn sa_is_con_textos_repetitivos() {
    let mut textos: Vec<Vec<u8>> = Vec::new();
    textos.push(b"ACG".repeat(200));
    textos.push(b"AAAAAAAAAC".repeat(50));
    textos.push(b"ABAABAAABAAAAB".repeat(20));
    let mut fibonacci = vec![b"A".to_vec(), b"AB".to_vec()];
    for _ in 0..10 {
        let mut siguiente = fibonacci[fibonacci.len() - 1].clone();
        siguiente.extend_from_slice(&fibonacci[fibonacci.len() - 2]);
        fibonacci.push(siguiente);
    }
    textos.push(fibonacci.pop().unwrap());
    textos.push(vec![b'Z'; 1000]);
    for texto in textos {
        assert_eq!(arreglo_de_sufijos(&texto).unwrap(), sufijos_a_mano(&texto));
    }
}

#[test]
fn bwt_a_mano_y_su_inversa() {
    assert_eq!(transformada(b"ACATG").unwrap(), b"G$CATA".to_vec());
    assert_eq!(inversa(b"G$CATA").unwrap(), b"ACATG".to_vec());
    assert_eq!(transformada(b"").unwrap(), b"$".to_vec());
    assert_eq!(inversa(b"$").unwrap(), Vec::<u8>::new());
    let (bwt, recuperado) = transformada_e_inversa(b"TACATACA").unwrap();
    assert_eq!(bwt, bwt_a_mano(b"TACATACA"));
    assert_eq!(recuperado, b"TACATACA".to_vec());
}

#[test]
fn bwt_con_letras_menores_que_el_centinela() {
    let texto = b"a b!c #a b";
    assert_eq!(transformada(texto).unwrap(), bwt_a_mano(texto));
    assert_eq!(
        inversa(&transformada(texto).unwrap()).unwrap(),
        texto.to_vec()
    );
}

#[test]
fn bwt_igual_a_rotaciones_ordenadas_y_vuelve() {
    let mut azar = Azar::nuevo(11);
    for prueba in 0..300 {
        let letras: &[u8] = if prueba % 2 == 0 { b"ACGT" } else { b"ab" };
        let longitud = azar.entero(120);
        let texto = texto_al_azar(&mut azar, longitud, letras);
        let bwt = transformada(&texto).unwrap();
        assert_eq!(bwt, bwt_a_mano(&texto));
        assert_eq!(inversa(&bwt).unwrap(), texto);
    }
}

#[test]
fn errores_de_la_bwt() {
    assert!(transformada(b"AC$GT").is_err());
    assert!(inversa(b"ACGT").is_err());
    assert!(inversa(b"A$C$").is_err());
    assert!(inversa(b"").is_err());
}

#[test]
fn cuenta_rachas() {
    assert_eq!(rachas(b""), 0);
    assert_eq!(rachas(b"A"), 1);
    assert_eq!(rachas(b"AAACCA"), 3);
    assert_eq!(rachas(b"ACGT"), 4);
}

#[test]
fn intervalo_a_mano() {
    let indice = IndiceFm::construir(b"ACATG").unwrap();
    assert_eq!(indice.filas, 6);
    let columna_l: Vec<u8> = (0..6).map(|fila| indice.simbolo_de_l(fila)).collect();
    assert_eq!(columna_l, vec![3, 0, 2, 1, 4, 1]);
    assert_eq!(indice.primera_fila, [0, 1, 3, 4, 5, 6, 6]);
    assert_eq!(indice.paso_atras(4, 0, 6), (5, 6));
    assert_eq!(indice.paso_atras(1, 5, 6), (2, 3));
    assert_eq!(indice.intervalo(b"AT"), (2, 3));
    assert_eq!(indice.contar(b"A"), 2);
    assert_eq!(indice.contar(b"CAT"), 1);
    assert_eq!(indice.contar(b"GA"), 0);
    assert_eq!(indice.posicion_de_fila(2), 2);
    assert_eq!(indice.simbolo_de_f(0), 0);
    assert_eq!(indice.simbolo_de_f(2), 1);
    assert_eq!(indice.simbolo_de_f(5), 4);
    let todas: Vec<usize> = (0..6).map(|fila| indice.posicion_de_fila(fila)).collect();
    assert_eq!(todas, vec![5, 0, 2, 1, 4, 3]);
}

#[test]
fn cotas_de_diferencias_a_mano() {
    let indice = IndiceDoble::construir(b"ACGTACGT").unwrap();
    assert_eq!(indice.cotas_de_diferencias(b"ACGA"), vec![0, 0, 0, 1]);
    assert_eq!(indice.cotas_de_diferencias(b"TTTT"), vec![0, 1, 1, 2]);
    assert_eq!(indice.cotas_de_diferencias(b"GTAC"), vec![0, 0, 0, 0]);
}

#[test]
fn las_cotas_nunca_pasan_de_las_diferencias_reales() {
    let mut azar = Azar::nuevo(31);
    for _ in 0..200 {
        let largo_del_texto = 1 + azar.entero(200);
        let texto = texto_al_azar(&mut azar, largo_del_texto, b"ACGT");
        let indice = IndiceDoble::construir(&texto).unwrap();
        let largo_del_patron = 1 + azar.entero(12);
        let patron = texto_al_azar(&mut azar, largo_del_patron, b"ACGT");
        let cotas = indice.cotas_de_diferencias(&patron);
        for largo in 1..=patron.len() {
            let mejor = texto
                .windows(largo)
                .map(|ventana| hamming(ventana, &patron[..largo]))
                .min()
                .unwrap_or(largo);
            assert!(cotas[largo - 1] <= mejor);
        }
    }
}

#[test]
fn ambas_hebras_a_mano() {
    let indice = IndiceDoble::construir(b"AAACCGTT").unwrap();
    let hallados = buscar_patrones(&indice, &[b"ACGG".to_vec()], 0).unwrap();
    assert_eq!(
        hallados[0],
        vec![Acierto {
            posicion: 3,
            hebra: b'-',
            diferencias: 0
        }]
    );
    let palindromo = buscar_patrones(&indice, &[b"CCGG".to_vec()], 1).unwrap();
    assert!(palindromo[0].iter().all(|acierto| acierto.hebra == b'+'));
    assert_eq!(
        palindromo[0],
        vec![Acierto {
            posicion: 3,
            hebra: b'+',
            diferencias: 1
        }]
    );
}

#[test]
fn fm_igual_que_escaneo_a_mano() {
    let mut azar = Azar::nuevo(23);
    for prueba in 0..400 {
        let letras: &[u8] = if prueba % 5 == 0 { b"ACGTN" } else { b"ACGT" };
        let longitud = 1 + azar.entero(400);
        let texto = texto_al_azar(&mut azar, longitud, letras);
        let indice = IndiceDoble::construir(&texto).unwrap();
        for _ in 0..5 {
            let d = azar.entero(3);
            let largo = d + 1 + azar.entero(8);
            let patron = if azar.entero(2) == 0 && texto.len() > largo {
                let inicio = azar.entero(texto.len() - largo);
                texto[inicio..inicio + largo]
                    .iter()
                    .map(|&l| if l == b'N' { b'A' } else { l })
                    .collect()
            } else {
                texto_al_azar(&mut azar, largo, b"ACGT")
            };
            let hallados = buscar_patrones(&indice, std::slice::from_ref(&patron), d).unwrap();
            assert_eq!(
                hallados[0],
                aciertos_a_mano(&texto, &patron, d),
                "texto {} patron {} d {d}",
                String::from_utf8_lossy(&texto),
                String::from_utf8_lossy(&patron)
            );
            assert_eq!(
                indice.directo.contar(&patron),
                posiciones_del_patron(&texto, &patron).len()
            );
        }
    }
}

#[test]
fn todas_las_filas_se_ubican() {
    let mut azar = Azar::nuevo(5);
    let texto = texto_al_azar(&mut azar, 5000, b"ACGT");
    let indice = IndiceFm::construir(&texto).unwrap();
    let sufijos = arreglo_de_sufijos(&texto).unwrap();
    let simbolos: Vec<u32> = texto
        .iter()
        .map(|&l| simbolo_de(l) as u32)
        .chain(std::iter::once(0))
        .collect();
    assert_eq!(ordenar_sufijos(&simbolos, SIMBOLOS), sufijos);
    for (fila, &inicio) in sufijos.iter().enumerate() {
        assert_eq!(indice.posicion_de_fila(fila), inicio as usize);
    }
    let columna_l: Vec<u8> = (0..indice.filas)
        .map(|fila| indice.simbolo_de_l(fila))
        .collect();
    for fila in 0..=indice.filas {
        for simbolo in 0..SIMBOLOS as u8 {
            let a_mano = columna_l[..fila].iter().filter(|&&s| s == simbolo).count();
            assert_eq!(indice.ocurrencias(simbolo, fila), a_mano);
        }
    }
}

#[test]
fn letras_minusculas_y_desconocidas() {
    let indice = IndiceDoble::construir(b"acgtNNacgt").unwrap();
    assert_eq!(indice.directo.contar(b"ACGT"), 2);
    let hallados = buscar_patrones(&indice, &[b"TAAC".to_vec()], 2).unwrap();
    assert_eq!(hallados[0], aciertos_a_mano(b"ACGTNNACGT", b"TAAC", 2));
}

#[test]
fn errores_de_la_busqueda() {
    let indice = IndiceDoble::construir(b"ACGTACGTAC").unwrap();
    assert!(IndiceFm::construir(b"").is_err());
    assert!(IndiceDoble::construir(b"").is_err());
    assert!(buscar_patrones(&indice, &[], 0).is_err());
    assert!(buscar_patrones(&indice, &[b"ACG".to_vec()], 3).is_err());
    assert!(buscar_patrones(&indice, &[b"ACN".to_vec()], 0).is_err());
    assert!(buscar_patrones(&indice, &[b"AC".to_vec()], 2).is_err());
    assert!(buscar_patrones(&indice, &[Vec::new()], 0).is_err());
    let largo = buscar_patrones(&indice, &[b"ACGTACGTACGTACGT".to_vec()], 2).unwrap();
    assert!(largo[0].is_empty());
}

#[test]
fn demasiadas_apariciones_dan_error() {
    let texto = vec![b'A'; 600_000];
    let indice = IndiceDoble::construir(&texto).unwrap();
    let error = buscar_patrones(&indice, &[b"A".to_vec()], 0).unwrap_err();
    assert!(error.contains("medio millón"));
    assert_eq!(indice.directo.contar(b"AAAA"), 599_997);
}

#[test]
fn ecoli_igual_que_el_escaneo_clasico() {
    let genoma = ecoli::genoma();
    let indice = IndiceDoble::construir(genoma).unwrap();
    assert_eq!(indice.directo.filas, genoma.len() + 1);
    let patrones: Vec<Vec<u8>> = vec![
        b"TTATCCACA".to_vec(),
        b"GCTGGTGG".to_vec(),
        genoma[1_000_000..1_000_020].to_vec(),
    ];
    let hallados = buscar_patrones(&indice, &patrones, 0).unwrap();
    for (patron, aciertos) in patrones.iter().zip(&hallados) {
        let directos = posiciones_del_patron(genoma, patron);
        let inversos = posiciones_del_patron(genoma, &complemento_inverso(patron));
        let mas: Vec<usize> = aciertos
            .iter()
            .filter(|a| a.hebra == b'+')
            .map(|a| a.posicion)
            .collect();
        let menos: Vec<usize> = aciertos
            .iter()
            .filter(|a| a.hebra == b'-')
            .map(|a| a.posicion)
            .collect();
        assert_eq!(mas, directos);
        assert_eq!(menos, inversos);
    }
    assert!(hallados[2].iter().any(|a| a.posicion == 1_000_000));
    let aproximados = buscar_patrones(&indice, &patrones[..1], 1).unwrap();
    let mas: Vec<usize> = aproximados[0]
        .iter()
        .filter(|a| a.hebra == b'+')
        .map(|a| a.posicion)
        .collect();
    assert_eq!(mas, posiciones_aproximadas(genoma, &patrones[0], 1));
}
