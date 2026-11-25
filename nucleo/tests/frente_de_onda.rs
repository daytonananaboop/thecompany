use nucleo::alineamiento::{alineamiento_afin, distancia_de_edicion};
use nucleo::azar::Azar;
use nucleo::frente_de_onda::*;

fn gotoh(v: &[u8], w: &[u8], castigos: Castigos) -> usize {
    let infinito = usize::MAX / 4;
    let abrir = castigos.apertura + castigos.extension;
    let alargar = castigos.extension;
    let columnas = w.len() + 1;
    let mut parejas = vec![vec![infinito; columnas]; v.len() + 1];
    let mut inserciones = vec![vec![infinito; columnas]; v.len() + 1];
    let mut borrados = vec![vec![infinito; columnas]; v.len() + 1];
    parejas[0][0] = 0;
    for i in 0..=v.len() {
        for j in 0..columnas {
            if j > 0 {
                inserciones[i][j] =
                    (parejas[i][j - 1] + abrir).min(inserciones[i][j - 1] + alargar);
            }
            if i > 0 {
                borrados[i][j] = (parejas[i - 1][j] + abrir).min(borrados[i - 1][j] + alargar);
            }
            if i == 0 && j == 0 {
                continue;
            }
            let mut mejor = inserciones[i][j].min(borrados[i][j]);
            if i > 0 && j > 0 {
                let cambio = if v[i - 1] == w[j - 1] {
                    0
                } else {
                    castigos.diferencia
                };
                mejor = mejor.min(parejas[i - 1][j - 1] + cambio);
            }
            parejas[i][j] = mejor;
        }
    }
    parejas[v.len()][w.len()]
}

fn levenshtein(v: &[u8], w: &[u8]) -> usize {
    let mut tabla = vec![vec![0; w.len() + 1]; v.len() + 1];
    for (i, fila) in tabla.iter_mut().enumerate() {
        fila[0] = i;
    }
    for (j, celda) in tabla[0].iter_mut().enumerate() {
        *celda = j;
    }
    for i in 1..=v.len() {
        for j in 1..=w.len() {
            let cambio = usize::from(v[i - 1] != w[j - 1]);
            tabla[i][j] = (tabla[i - 1][j - 1] + cambio)
                .min(tabla[i - 1][j] + 1)
                .min(tabla[i][j - 1] + 1);
        }
    }
    tabla[v.len()][w.len()]
}

fn al_azar(azar: &mut Azar, minimo: usize, extra: usize, letras: &[u8]) -> Vec<u8> {
    let largo = minimo + azar.entero(extra + 1);
    (0..largo)
        .map(|_| letras[azar.entero(letras.len())])
        .collect()
}

fn mutar(azar: &mut Azar, original: &[u8], tasa: f64, letras: &[u8]) -> Vec<u8> {
    let mut copia = Vec::with_capacity(original.len() + 8);
    for &letra in original {
        if azar.decimal() >= tasa {
            copia.push(letra);
            continue;
        }
        match azar.entero(4) {
            0 => {}
            1 => {
                copia.push(letra);
                for _ in 0..=azar.entero(4) {
                    copia.push(letras[azar.entero(letras.len())]);
                }
            }
            _ => copia.push(letras[azar.entero(letras.len())]),
        }
    }
    copia
}

fn sin_huecos(fila: &[u8]) -> Vec<u8> {
    fila.iter()
        .copied()
        .filter(|&letra| letra != b'-')
        .collect()
}

fn comprobar(v: &[u8], w: &[u8], castigos: Castigos) {
    let resultado = alinear_wfa(v, w, castigos).unwrap();
    let esperado = gotoh(v, w, castigos);
    assert_eq!(
        resultado.penalizacion,
        esperado,
        "v={} w={} castigos={castigos:?}",
        String::from_utf8_lossy(v),
        String::from_utf8_lossy(w)
    );
    let operaciones = resultado.operaciones.unwrap();
    assert_eq!(penalizacion_de(&operaciones, castigos), esperado);
    let (fila_v, fila_w) = filas_alineadas(v, w, &operaciones);
    assert_eq!(sin_huecos(&fila_v), v);
    assert_eq!(sin_huecos(&fila_w), w);
    for ((&a, &b), &operacion) in fila_v.iter().zip(&fila_w).zip(&operaciones) {
        match operacion {
            b'=' => assert_eq!(a, b),
            b'X' => assert!(a != b && a != b'-' && b != b'-'),
            b'I' => assert_eq!(a, b'-'),
            b'D' => assert_eq!(b, b'-'),
            _ => panic!("operación desconocida"),
        }
    }
}

fn diferencia_cuatro(a: u8, b: u8) -> i32 {
    if a == b { 0 } else { -4 }
}

#[test]
fn secuencias_identicas_no_cuestan_nada() {
    let resultado = alinear_wfa(b"GATTACA", b"GATTACA", CASTIGOS_DEL_ARTICULO).unwrap();
    assert_eq!(resultado.penalizacion, 0);
    assert_eq!(resultado.frentes, 1);
    assert_eq!(cigar(&resultado.operaciones.unwrap()), "7=");
}

#[test]
fn una_letra_distinta_cuesta_la_diferencia() {
    let resultado = alinear_wfa(b"ACGTACGT", b"ACGAACGT", CASTIGOS_DEL_ARTICULO).unwrap();
    assert_eq!(resultado.penalizacion, 4);
    assert_eq!(cigar(&resultado.operaciones.unwrap()), "3=1X4=");
}

#[test]
fn un_hueco_cuesta_abrir_mas_alargar() {
    let resultado = alinear_wfa(b"CCCCAAAA", b"CCCCGAAAA", CASTIGOS_DEL_ARTICULO).unwrap();
    assert_eq!(resultado.penalizacion, 8);
    assert_eq!(cigar(&resultado.operaciones.unwrap()), "4=1I4=");
    let resultado = alinear_wfa(b"CCCCGGGAAAA", b"CCCCAAAA", CASTIGOS_DEL_ARTICULO).unwrap();
    assert_eq!(resultado.penalizacion, 6 + 3 * 2);
    assert_eq!(cigar(&resultado.operaciones.unwrap()), "4=3D4=");
}

#[test]
fn un_hueco_largo_sale_mas_barato_que_dos_cortos() {
    let resultado = alinear_wfa(b"TTTTTTGGCCCCCC", b"TTTTTTCCCCCC", CASTIGOS_DEL_ARTICULO).unwrap();
    assert_eq!(resultado.penalizacion, 10);
    assert_eq!(cigar(&resultado.operaciones.unwrap()), "6=2D6=");
}

#[test]
fn con_una_secuencia_vacia_todo_es_hueco() {
    let resultado = alinear_wfa(b"", b"", CASTIGOS_DEL_ARTICULO).unwrap();
    assert_eq!(resultado.penalizacion, 0);
    assert_eq!(resultado.operaciones.unwrap(), Vec::<u8>::new());
    let resultado = alinear_wfa(b"", b"ACG", CASTIGOS_DEL_ARTICULO).unwrap();
    assert_eq!(resultado.penalizacion, 12);
    assert_eq!(cigar(&resultado.operaciones.unwrap()), "3I");
    let resultado = alinear_wfa(b"ACGT", b"", CASTIGOS_DEL_ARTICULO).unwrap();
    assert_eq!(resultado.penalizacion, 14);
    assert_eq!(cigar(&resultado.operaciones.unwrap()), "4D");
}

#[test]
fn sin_coincidencias_elige_lo_mas_barato() {
    let resultado = alinear_wfa(b"AAAA", b"CCCC", CASTIGOS_DEL_ARTICULO).unwrap();
    assert_eq!(resultado.penalizacion, 16);
    let caro = Castigos {
        diferencia: 20,
        apertura: 1,
        extension: 1,
    };
    let resultado = alinear_wfa(b"AAAA", b"CCCC", caro).unwrap();
    assert_eq!(resultado.penalizacion, 10);
    let forma = cigar(&resultado.operaciones.unwrap());
    assert!(forma == "4I4D" || forma == "4D4I", "{forma}");
}

#[test]
fn castigos_y_largos_imposibles_se_rechazan() {
    let sin_diferencia = Castigos {
        diferencia: 0,
        apertura: 6,
        extension: 2,
    };
    assert!(alinear_wfa(b"AC", b"AG", sin_diferencia).is_err());
    let sin_extension = Castigos {
        diferencia: 4,
        apertura: 6,
        extension: 0,
    };
    assert!(alinear_wfa(b"AC", b"AG", sin_extension).is_err());
    let enorme = Castigos {
        diferencia: 4,
        apertura: 5000,
        extension: 2,
    };
    assert!(alinear_wfa(b"AC", b"AG", enorme).is_err());
    let largo = vec![b'A'; 60_000];
    assert!(alinear_wfa(b"A", &largo, CASTIGOS_DEL_ARTICULO).is_err());
}

#[test]
fn las_cuentas_de_celdas_crecen_con_el_puntaje() {
    let v = b"ACGTTGCAACGTTGCAACGT";
    let w = b"ACGTAGCAACGTTGCATCGT";
    let resultado = alinear_wfa(v, w, CASTIGOS_DEL_ARTICULO).unwrap();
    assert_eq!(resultado.penalizacion, 8);
    assert!(resultado.frentes >= 3);
    assert!(resultado.celdas < v.len() * w.len());
}

#[test]
fn wfa_igual_a_gotoh_en_casos_al_azar() {
    let mut azar = Azar::nuevo(2021);
    for caso in 0..1500 {
        let letras: &[u8] = match caso % 3 {
            0 => b"AC",
            1 => b"ACGT",
            _ => b"ACDEFGHIKLMNPQRSTVWY",
        };
        let v = al_azar(&mut azar, 0, 39, letras);
        let w = if caso % 2 == 0 {
            mutar(&mut azar, &v, 0.15, letras)
        } else {
            al_azar(&mut azar, 0, 39, letras)
        };
        let castigos = Castigos {
            diferencia: 1 + azar.entero(9),
            apertura: azar.entero(12),
            extension: 1 + azar.entero(5),
        };
        comprobar(&v, &w, castigos);
    }
}

#[test]
fn wfa_igual_a_gotoh_en_secuencias_largas_y_parecidas() {
    let mut azar = Azar::nuevo(7);
    for _ in 0..12 {
        let v = al_azar(&mut azar, 300, 399, b"ACGT");
        let w = mutar(&mut azar, &v, 0.04, b"ACGT");
        comprobar(&v, &w, CASTIGOS_DEL_ARTICULO);
    }
}

#[test]
fn wfa_coincide_con_el_afin_clasico_del_capitulo() {
    let mut azar = Azar::nuevo(99);
    for _ in 0..60 {
        let v = al_azar(&mut azar, 0, 79, b"ACGT");
        let w = mutar(&mut azar, &v, 0.1, b"ACGT");
        let resultado = alinear_wfa(&v, &w, CASTIGOS_DEL_ARTICULO).unwrap();
        let clasico = alineamiento_afin(&v, &w, diferencia_cuatro, 8, 2);
        assert_eq!(resultado.penalizacion as i64, -clasico.puntaje);
    }
}

#[test]
fn sin_memoria_queda_el_puntaje_exacto() {
    let mut azar = Azar::nuevo(31);
    let v = al_azar(&mut azar, 12_000, 0, b"ACGT");
    let w = al_azar(&mut azar, 12_000, 0, b"ACGT");
    let castigos = Castigos {
        diferencia: 1,
        apertura: 0,
        extension: 1,
    };
    let resultado = alinear_wfa(&v, &w, castigos).unwrap();
    assert!(resultado.operaciones.is_none());
    assert_eq!(resultado.penalizacion, distancia_de_edicion(&v, &w));
}

#[test]
fn letras_iguales_cuenta_de_ocho_en_ocho() {
    assert_eq!(letras_iguales(b"ACGTACGTACGTA", b"ACGTACGTACGTT"), 12);
    assert_eq!(letras_iguales(b"ACGTACGTAC", b"ACGTACGTAC"), 10);
    assert_eq!(letras_iguales(b"ACGTTCGTAC", b"ACGTACGTAC"), 4);
    assert_eq!(letras_iguales(b"", b"ACG"), 0);
}

#[test]
fn myers_en_casos_a_mano() {
    assert_eq!(distancia_myers(b"", b"").distancia, 0);
    assert_eq!(distancia_myers(b"", b"ACG").distancia, 3);
    assert_eq!(distancia_myers(b"ACGT", b"").distancia, 4);
    assert_eq!(distancia_myers(b"CASA", b"CAPA").distancia, 1);
    assert_eq!(distancia_myers(b"GATO", b"PATOS").distancia, 2);
    assert_eq!(distancia_myers(b"ACGT", b"TGCA").distancia, 4);
    assert_eq!(distancia_myers(b"AAAA", b"AAAA").distancia, 0);
}

#[test]
fn myers_igual_a_levenshtein_en_los_bordes_de_palabra() {
    let mut azar = Azar::nuevo(1999);
    for largo in [1, 2, 63, 64, 65, 127, 128, 129, 191, 192, 193, 300] {
        for _ in 0..6 {
            let v = al_azar(&mut azar, largo, 0, b"ACGT");
            let w = mutar(&mut azar, &v, 0.2, b"ACGT");
            assert_eq!(
                distancia_myers(&v, &w).distancia,
                levenshtein(&v, &w),
                "largo {largo}"
            );
            let otra = al_azar(&mut azar, 0, 2 * largo, b"ACGT");
            assert_eq!(distancia_myers(&v, &otra).distancia, levenshtein(&v, &otra));
            assert_eq!(distancia_myers(&otra, &v).distancia, levenshtein(&otra, &v));
        }
    }
}

#[test]
fn myers_igual_a_levenshtein_al_azar() {
    let mut azar = Azar::nuevo(64);
    for caso in 0..800 {
        let letras: &[u8] = if caso % 2 == 0 {
            b"ACGT"
        } else {
            b"ACDEFGHIKLMNPQRSTVWY"
        };
        let v = al_azar(&mut azar, 0, 149, letras);
        let w = if caso % 3 == 0 {
            al_azar(&mut azar, 0, 149, letras)
        } else {
            mutar(&mut azar, &v, 0.1, letras)
        };
        assert_eq!(distancia_myers(&v, &w).distancia, levenshtein(&v, &w));
    }
}

#[test]
fn myers_igual_a_la_distancia_clasica_del_capitulo() {
    let mut azar = Azar::nuevo(5);
    for _ in 0..10 {
        let v = al_azar(&mut azar, 1000, 999, b"ACGT");
        let w = mutar(&mut azar, &v, 0.05, b"ACGT");
        assert_eq!(
            distancia_myers(&v, &w).distancia,
            distancia_de_edicion(&v, &w)
        );
    }
}

#[test]
fn tabla_de_iguales_marca_cada_posicion() {
    let tabla = TablaDeIguales::para(b"ACCA");
    assert_eq!(tabla.palabras, 1);
    assert_eq!(tabla.de(b'A'), &[0b1001]);
    assert_eq!(tabla.de(b'C'), &[0b0110]);
    assert_eq!(tabla.de(b'G'), &[0]);
}

#[test]
fn avanzar_palabra_reproduce_una_columna_a_mano() {
    let tabla = TablaDeIguales::para(b"AC");
    let mut palabra = Palabra::inicial();
    let salida = avanzar_palabra(&mut palabra, tabla.de(b'A')[0], 1, 0b10);
    assert_eq!(salida, -1);
    assert_eq!(palabra.vertical_mas & 0b11, 0b10);
    assert_eq!(palabra.vertical_menos & 0b11, 0b01);
}

#[test]
fn cigar_y_penalizacion_de_operaciones() {
    let operaciones = b"===XII==D";
    assert_eq!(cigar(operaciones), "3=1X2I2=1D");
    assert_eq!(
        penalizacion_de(operaciones, CASTIGOS_DEL_ARTICULO),
        4 + 10 + 8
    );
    assert_eq!(cigar(b""), "");
}

#[test]
fn castigos_enormes_siguen_igual_a_gotoh() {
    let mut azar = Azar::nuevo(77);
    let grandes = [
        Castigos {
            diferencia: 1000,
            apertura: 1000,
            extension: 1000,
        },
        Castigos {
            diferencia: 1000,
            apertura: 0,
            extension: 1,
        },
        Castigos {
            diferencia: 1,
            apertura: 1000,
            extension: 1000,
        },
    ];
    for castigos in grandes {
        for _ in 0..10 {
            let v = al_azar(&mut azar, 20, 60, b"ACGT");
            let w = mutar(&mut azar, &v, 0.2, b"ACGT");
            comprobar(&v, &w, castigos);
        }
    }
}

#[test]
fn muchos_puntajes_con_pocas_celdas_no_agotan_la_memoria() {
    let v = vec![b'A'; 3000];
    let w = vec![b'C'; 3000];
    let castigos = Castigos {
        diferencia: 1000,
        apertura: 1000,
        extension: 1000,
    };
    let resultado = alinear_wfa(&v, &w, castigos).unwrap();
    assert_eq!(resultado.penalizacion, 3_000_000);
}
