use nucleo::azar::Azar;
use nucleo::cebadores::*;
use nucleo::secuencia::complemento_inverso;
use std::cmp::Ordering;

fn cerca(medido: f64, esperado: f64, tolerancia: f64) -> bool {
    (medido - esperado).abs() <= tolerancia
}

fn un_molar() -> Condiciones {
    Condiciones {
        monovalentes_mm: 1000.0,
        magnesio_mm: 0.0,
        dntp_mm: 0.0,
        ..Condiciones::default()
    }
}

fn con_owczarzy(condiciones: Condiciones) -> Condiciones {
    Condiciones {
        correccion: CorreccionSal::Owczarzy,
        ..condiciones
    }
}

fn adn_al_azar(longitud: usize, semilla: u64) -> Vec<u8> {
    let mut azar = Azar::nuevo(semilla);
    (0..longitud).map(|_| b"ACGT"[azar.entero(4)]).collect()
}

#[test]
fn condiciones_por_defecto_son_las_de_primer3() {
    let condiciones = Condiciones::default();
    assert_eq!(condiciones.monovalentes_mm, 50.0);
    assert_eq!(condiciones.magnesio_mm, 1.5);
    assert_eq!(condiciones.dntp_mm, 0.6);
    assert_eq!(condiciones.oligo_nm, 50.0);
    assert_eq!(condiciones.correccion, CorreccionSal::SantaLucia);
}

#[test]
fn duplex_suma_vecinos_iniciacion_y_simetria() {
    let acgt = duplex(b"ACGT");
    assert!(cerca(acgt.entalpia, -22.8, 1e-9));
    assert!(cerca(acgt.entropia, -65.2, 1e-9));
    assert_eq!(acgt.longitud, 4);
    assert_eq!(acgt.bases_gc, 2);
    assert!(acgt.autocomplementario);
    let largo = duplex(b"GATTCGACCTAGGCATCAGT");
    assert!(cerca(largo.entalpia, -153.5, 1e-9));
    assert!(cerca(largo.entropia, -416.4, 1e-9));
    assert_eq!(largo.bases_gc, 10);
    assert!(!largo.autocomplementario);
    let vacio = duplex(b"");
    assert_eq!(vacio.entalpia, 0.0);
    assert_eq!(vacio.longitud, 0);
}

#[test]
fn autocomplementario_solo_con_pares_de_watson_y_crick() {
    assert!(es_autocomplementario(b"GAATTC"));
    assert!(es_autocomplementario(b"CCAGTTGAATTCAACTGG"));
    assert!(!es_autocomplementario(b"GAATTCA"));
    assert!(!es_autocomplementario(b"NN"));
    assert!(!es_autocomplementario(b""));
}

#[test]
fn sodio_equivalente_y_magnesio_libre() {
    let condiciones = Condiciones::default();
    assert!(cerca(
        sodio_equivalente(&condiciones),
        50.0 + 120.0 * 0.9f64.sqrt(),
        1e-9
    ));
    let sin_magnesio_util = Condiciones {
        magnesio_mm: 0.5,
        ..condiciones
    };
    assert!(cerca(sodio_equivalente(&sin_magnesio_util), 50.0, 1e-12));
    assert!(cerca(magnesio_libre(&condiciones), 0.9, 1e-12));
    let sin_dntp = Condiciones {
        dntp_mm: 0.0,
        ..condiciones
    };
    assert!(cerca(magnesio_libre(&sin_dntp), 1.5, 1e-12));
    assert_eq!(magnesio_libre(&sin_magnesio_util), 0.0);
    let iguales = Condiciones {
        magnesio_mm: 0.6,
        ..condiciones
    };
    assert_eq!(magnesio_libre(&iguales), 0.0);
}

#[test]
fn tm_santalucia_como_primer3() {
    let condiciones = Condiciones::default();
    let tm = temperatura_fusion(b"GATTCGACCTAGGCATCAGT", &condiciones).unwrap();
    assert!(cerca(tm, 56.8113, 0.001), "{tm}");
    let palindromo = temperatura_fusion(b"CCAGTTGAATTCAACTGG", &condiciones).unwrap();
    assert!(cerca(palindromo, 51.8878, 0.001), "{palindromo}");
    let largo = temperatura_fusion(b"TTTGCAGGCTTACGATAAGCCAGTC", &condiciones).unwrap();
    assert!(cerca(largo, 63.5876, 0.001), "{largo}");
    let en_un_molar = temperatura_fusion(b"GATTCGACCTAGGCATCAGT", &un_molar()).unwrap();
    assert!(cerca(en_un_molar, 66.0327, 0.001), "{en_un_molar}");
}

#[test]
fn el_magnesio_sube_la_tm_y_los_dntp_lo_atrapan() {
    let oligo = b"GATTCGACCTAGGCATCAGT";
    let con = |magnesio_mm: f64, dntp_mm: f64| {
        temperatura_fusion(
            oligo,
            &Condiciones {
                magnesio_mm,
                dntp_mm,
                ..Condiciones::default()
            },
        )
        .unwrap()
    };
    assert!(cerca(con(0.0, 0.6), 51.0284, 0.001));
    assert!(cerca(con(3.0, 0.6), 58.6291, 0.001));
    assert_eq!(con(0.5, 0.6), con(0.0, 0.0));
    assert!(con(3.0, 0.6) > con(1.5, 0.6));
    assert!(con(1.5, 0.6) > con(0.0, 0.6));
}

#[test]
fn tm_owczarzy_2008() {
    let oligo = b"GATTCGACCTAGGCATCAGT";
    let normal = temperatura_fusion(oligo, &con_owczarzy(Condiciones::default())).unwrap();
    assert!(cerca(normal, 56.2644, 0.001), "{normal}");
    let solo_sodio = temperatura_fusion(
        oligo,
        &con_owczarzy(Condiciones {
            magnesio_mm: 0.0,
            dntp_mm: 0.0,
            ..Condiciones::default()
        }),
    )
    .unwrap();
    assert!(cerca(solo_sodio, 50.8210, 0.001), "{solo_sodio}");
    let solo_magnesio = temperatura_fusion(
        oligo,
        &con_owczarzy(Condiciones {
            monovalentes_mm: 0.0,
            magnesio_mm: 10.0,
            dntp_mm: 0.0,
            ..Condiciones::default()
        }),
    )
    .unwrap();
    assert!(cerca(solo_magnesio, 60.3957, 0.001), "{solo_magnesio}");
    let mucho_magnesio = temperatura_fusion(
        oligo,
        &con_owczarzy(Condiciones {
            magnesio_mm: 100.0,
            dntp_mm: 0.0,
            ..Condiciones::default()
        }),
    )
    .unwrap();
    assert!(cerca(mucho_magnesio, 61.2433, 0.001), "{mucho_magnesio}");
    let con = |magnesio_mm: f64, dntp_mm: f64| {
        temperatura_fusion(
            oligo,
            &con_owczarzy(Condiciones {
                monovalentes_mm: 10.0,
                magnesio_mm,
                dntp_mm,
                ..Condiciones::default()
            }),
        )
        .unwrap()
    };
    assert!(cerca(con(0.0, 0.0), 36.3767, 0.001), "{}", con(0.0, 0.0));
    assert_eq!(con(0.5, 0.6), con(0.0, 0.0));
    assert_eq!(con(0.6, 0.6), con(0.0, 0.0));
    assert!(cerca(con(0.61, 0.6), 39.6824, 0.001), "{}", con(0.61, 0.6));
}

#[test]
fn la_tm_es_la_misma_para_el_oligo_y_su_complemento() {
    for semilla in 0..50 {
        let oligo = adn_al_azar(18 + (semilla as usize % 13), semilla);
        let inverso = complemento_inverso(&oligo);
        for condiciones in [Condiciones::default(), con_owczarzy(Condiciones::default())] {
            let directa = temperatura_fusion(&oligo, &condiciones).unwrap();
            let opuesta = temperatura_fusion(&inverso, &condiciones).unwrap();
            assert!(cerca(directa, opuesta, 1e-9));
        }
    }
}

#[test]
fn tm_rechaza_entradas_malas() {
    let normal = Condiciones::default();
    assert!(temperatura_fusion(b"", &normal).is_err());
    assert!(temperatura_fusion(b"A", &normal).is_err());
    assert!(temperatura_fusion(b"ACGTN", &normal).is_err());
    assert!(temperatura_fusion(b"acgt", &normal).is_err());
    assert!(temperatura_fusion(&[b'A'; 61], &normal).is_err());
    assert!(temperatura_fusion(&[b'A'; 60], &normal).is_ok());
    let malas = [
        Condiciones {
            oligo_nm: 0.0,
            ..normal
        },
        Condiciones {
            magnesio_mm: -1.0,
            ..normal
        },
        Condiciones {
            monovalentes_mm: f64::NAN,
            ..normal
        },
        Condiciones {
            dntp_mm: f64::INFINITY,
            ..normal
        },
        Condiciones {
            monovalentes_mm: 0.0,
            magnesio_mm: 0.0,
            ..normal
        },
        Condiciones {
            monovalentes_mm: 0.0,
            magnesio_mm: 0.6,
            dntp_mm: 0.6,
            correccion: CorreccionSal::Owczarzy,
            ..normal
        },
    ];
    for condiciones in malas {
        assert!(
            temperatura_fusion(b"GATTCGACCTAGGCATCAGT", &condiciones).is_err(),
            "{condiciones:?}"
        );
    }
}

#[test]
fn curva_de_fusion_cruza_la_mitad_en_la_tm() {
    let oligo = b"GATTCGACCTAGGCATCAGT";
    for condiciones in [Condiciones::default(), con_owczarzy(Condiciones::default())] {
        let tm = temperatura_fusion(oligo, &condiciones).unwrap();
        let curva = curva_de_fusion(oligo, &condiciones).unwrap();
        assert_eq!(curva.len(), 101);
        assert!(cerca(curva[0].0, tm - 25.0, 1e-9));
        assert!(cerca(curva[50].0, tm, 1e-9));
        assert!(cerca(curva[50].1, 0.5, 1e-9));
        assert!(curva[0].1 > 0.99);
        assert!(curva[100].1 < 0.01);
        assert!(curva.windows(2).all(|par| par[0].1 >= par[1].1));
    }
    let palindromo = b"CCAGTTGAATTCAACTGG";
    let tm = temperatura_fusion(palindromo, &Condiciones::default()).unwrap();
    let curva = curva_de_fusion(palindromo, &Condiciones::default()).unwrap();
    assert!(cerca(curva[50].0, tm, 1e-9));
    assert!(cerca(curva[50].1, 0.5, 1e-9));
    assert!(curva_de_fusion(b"A", &Condiciones::default()).is_err());
}

#[test]
fn autodimero_en_el_extremo_3() {
    let oligo = b"AAAAAAGCGC";
    let resultado = dimero(oligo, oligo, &un_molar()).unwrap();
    let esperado =
        2.0 * (-9.8 + 0.31015 * 24.4) + (-10.6 + 0.31015 * 27.2) + 2.0 * (0.1 + 0.31015 * 2.8);
    assert!(cerca(resultado.mas_estable.energia, esperado, 1e-6));
    assert!(cerca(esperado, -4.6918, 0.001));
    assert_eq!(resultado.mas_estable.pares, 4);
    assert_eq!(resultado.mas_estable.en_a, 6);
    assert_eq!(resultado.mas_estable.en_b, 9);
    assert_eq!(resultado.extremo_3, resultado.mas_estable);
}

#[test]
fn autodimero_interno_no_toca_el_extremo_3() {
    let oligo = b"AAAAGCGCAAAA";
    let resultado = dimero(oligo, oligo, &un_molar()).unwrap();
    assert!(cerca(resultado.mas_estable.energia, -4.6918, 0.001));
    assert_eq!(resultado.mas_estable.pares, 4);
    assert_eq!(resultado.mas_estable.en_a, 4);
    assert_eq!(resultado.mas_estable.en_b, 7);
    assert_eq!(resultado.extremo_3.pares, 0);
    assert_eq!(resultado.extremo_3.energia, 0.0);
}

#[test]
fn dimero_cruzado_completo() {
    let resultado = dimero(b"TCCCCCCT", b"AGGGGGGA", &un_molar()).unwrap();
    let esperado = (-8.2 + 0.31015 * 22.2)
        + 5.0 * (-8.0 + 0.31015 * 19.9)
        + (-7.8 + 0.31015 * 21.0)
        + 2.0 * (2.3 - 0.31015 * 4.1);
    assert!(cerca(resultado.mas_estable.energia, esperado, 1e-6));
    assert!(cerca(esperado, -9.6848, 0.001));
    assert_eq!(resultado.mas_estable.pares, 8);
    assert_eq!(resultado.mas_estable.en_a, 0);
    assert_eq!(resultado.mas_estable.en_b, 7);
    let al_reves = dimero(b"AGGGGGGA", b"TCCCCCCT", &un_molar()).unwrap();
    assert!(cerca(
        al_reves.mas_estable.energia,
        resultado.mas_estable.energia,
        1e-9
    ));
}

#[test]
fn la_sal_debilita_poco_el_dimero() {
    let oligo = b"AAAAAAGCGC";
    let resultado = dimero(oligo, oligo, &Condiciones::default()).unwrap();
    assert!(cerca(resultado.mas_estable.energia, -4.0723, 0.001));
}

#[test]
fn con_poca_sal_un_tramo_mas_corto_puede_ser_el_mas_estable() {
    let poca_sal = Condiciones {
        monovalentes_mm: 5.0,
        magnesio_mm: 0.0,
        dntp_mm: 0.6,
        ..Condiciones::default()
    };
    let resultado = dimero(b"TACAGCTGC", b"GCAGCTGTA", &poca_sal).unwrap();
    assert!(cerca(resultado.mas_estable.energia, -5.1803, 0.001));
    assert!(cerca(resultado.extremo_3.energia, -5.1803, 0.001));
    assert!(resultado.mas_estable.pares < 9);
    let otro = dimero(b"CCCTGAGGTAG", b"CGTTGTACCATGCAGTACAAAACTAG", &poca_sal).unwrap();
    assert!(cerca(otro.mas_estable.energia, -0.0744, 0.001));
}

#[test]
fn sin_complementariedad_no_hay_dimero() {
    let resultado = dimero(b"AAAAAAAA", b"AAAAAAAA", &un_molar()).unwrap();
    assert_eq!(resultado.mas_estable.pares, 0);
    assert_eq!(resultado.mas_estable.energia, 0.0);
    let vacio = dimero(b"", b"ACGT", &un_molar()).unwrap();
    assert_eq!(vacio.mas_estable.pares, 0);
    let con_n = dimero(b"NNNN", b"NNNN", &un_molar()).unwrap();
    assert_eq!(con_n.mas_estable.pares, 0);
    let malas = Condiciones {
        oligo_nm: -1.0,
        ..Condiciones::default()
    };
    assert!(dimero(b"ACGT", b"ACGT", &malas).is_err());
}

#[test]
fn horquilla_con_tallo_de_cuatro_y_lazo_de_cuatro() {
    let resultado = horquilla(b"GCGCAAAAGCGC", &un_molar()).unwrap();
    let tallo = 2.0 * (-9.8 + 0.31015 * 24.4) + (-10.6 + 0.31015 * 27.2);
    assert!(cerca(resultado.energia, tallo + 3.5, 1e-6));
    assert!(cerca(resultado.energia, -3.1286, 0.001));
    assert_eq!(resultado.pares, 4);
    assert_eq!(resultado.lazo, 4);
    assert_eq!(resultado.inicio, 0);
    assert_eq!(resultado.fin, 11);
}

#[test]
fn horquilla_necesita_tallo_y_lazo_de_tres() {
    let corto = horquilla(b"GGGAACCC", &un_molar()).unwrap();
    assert_eq!(corto.pares, 0);
    assert_eq!(corto.energia, 0.0);
    let nada = horquilla(b"AAAAAAAAAAAA", &un_molar()).unwrap();
    assert_eq!(nada.pares, 0);
    let vacio = horquilla(b"", &un_molar()).unwrap();
    assert_eq!(vacio.pares, 0);
    let tres = horquilla(b"GCGCAAGCGC", &un_molar()).unwrap();
    assert_eq!(tres.pares, 3);
    assert_eq!(tres.lazo, 4);
    assert_eq!((tres.inicio, tres.fin), (0, 9));
    assert!(cerca(tres.energia, -2.2323 - 2.1639 + 3.5, 0.001));
}

#[test]
fn lazos_largos_cuestan_mas() {
    let mut anterior = horquilla(b"GGGGGAAAAAACCCCC", &un_molar()).unwrap();
    for lazo in 7..=30 {
        let mut oligo = b"GGGGG".to_vec();
        oligo.extend(std::iter::repeat_n(b'A', lazo));
        oligo.extend(b"CCCCC");
        let actual = horquilla(&oligo, &un_molar()).unwrap();
        assert_eq!(actual.lazo, lazo);
        assert!(actual.energia >= anterior.energia - 0.25, "lazo {lazo}");
        anterior = actual;
    }
}

#[test]
fn cuenta_sitios_en_ambas_hebras() {
    let extremo: &[u8] = b"ACGTTGCAAGGCTTA";
    let mut genoma = b"AAAAAAAAAA".to_vec();
    genoma.extend(extremo);
    genoma.extend(b"AAAAAAAAAA");
    genoma.extend(complemento_inverso(extremo));
    genoma.extend(b"AAAAAAAAAA");
    let ausente: &[u8] = b"CCCCCCCCCCCCCCC";
    let con_n: &[u8] = b"ACGTTGCAAGGCTTN";
    let cuentas = contar_sitios(&genoma, &[extremo, ausente, con_n, extremo]);
    assert_eq!(cuentas, vec![2, 0, 0, 2]);
    assert!(contar_sitios(&genoma, &[]).is_empty());
    assert_eq!(contar_sitios(b"", &[extremo]), vec![0]);
    assert_eq!(
        contar_sitios(&genoma, &[b"ACG".as_slice(), b"ACGT".as_slice()]),
        vec![0, 0]
    );
}

#[test]
fn cuenta_sitios_igual_que_buscando_uno_por_uno() {
    let genoma = adn_al_azar(200_003, 11);
    let mut azar = Azar::nuevo(5);
    let mut extremos: Vec<Vec<u8>> = Vec::new();
    for _ in 0..40 {
        let inicio = azar.entero(genoma.len() - 15);
        extremos.push(genoma[inicio..inicio + 15].to_vec());
    }
    extremos.push(genoma[..15].to_vec());
    extremos.push(genoma[genoma.len() - 15..].to_vec());
    let repetido = genoma[1000..1015].to_vec();
    let mut con_copias = genoma.clone();
    for lugar in [50_000, 100_000, 150_000] {
        con_copias[lugar..lugar + 15].copy_from_slice(&complemento_inverso(&repetido));
    }
    extremos.push(repetido);
    let referencias: Vec<&[u8]> = extremos.iter().map(|e| e.as_slice()).collect();
    let rapido = contar_sitios(&con_copias, &referencias);
    for (extremo, cuenta) in extremos.iter().zip(&rapido) {
        let inverso = complemento_inverso(extremo);
        let lento = con_copias
            .windows(15)
            .filter(|ventana| *ventana == extremo.as_slice() || *ventana == inverso.as_slice())
            .count();
        assert_eq!(*cuenta, lento);
    }
    assert!(*rapido.last().unwrap() >= 4);
}

fn racha_maxima(secuencia: &[u8]) -> usize {
    let mut maxima = 0;
    let mut actual = 0;
    for (posicion, base) in secuencia.iter().enumerate() {
        if posicion > 0 && secuencia[posicion - 1] == *base {
            actual += 1;
        } else {
            actual = 1;
        }
        maxima = maxima.max(actual);
    }
    maxima
}

fn gc_de(secuencia: &[u8]) -> usize {
    secuencia
        .iter()
        .filter(|base| matches!(base, b'G' | b'C'))
        .count()
}

fn exceso(energia: f64, umbral: f64) -> f64 {
    (umbral - energia).max(0.0)
}

struct Candidato {
    inicio: usize,
    longitud: usize,
    secuencia: Vec<u8>,
    tm: f64,
    penalizacion: f64,
}

fn candidatos_a_mano(
    plantilla: &[u8],
    opciones: &OpcionesDiseno,
    directo: bool,
    genoma: Option<&[u8]>,
) -> Vec<Candidato> {
    let mut lista = Vec::new();
    for inicio in 0..plantilla.len() {
        for longitud in opciones.longitud_minima..=opciones.longitud_maxima {
            if inicio + longitud > plantilla.len() {
                continue;
            }
            let tramo = &plantilla[inicio..inicio + longitud];
            if tramo.iter().any(|b| !b"ACGT".contains(b)) {
                continue;
            }
            let secuencia = if directo {
                tramo.to_vec()
            } else {
                complemento_inverso(tramo)
            };
            let gc = gc_de(&secuencia) as f64 / longitud as f64;
            let pinza = gc_de(&secuencia[longitud - 5..]);
            if racha_maxima(&secuencia) > 4
                || !(0.4..=0.6).contains(&gc)
                || !(1..=3).contains(&pinza)
            {
                continue;
            }
            let tm = temperatura_fusion(&secuencia, &opciones.condiciones).unwrap();
            if tm < opciones.tm_minima || tm > opciones.tm_maxima {
                continue;
            }
            let propio = dimero(&secuencia, &secuencia, &opciones.condiciones).unwrap();
            let pliegue = horquilla(&secuencia, &opciones.condiciones).unwrap();
            let mut penalizacion = (tm - opciones.tm_optima).abs()
                + longitud.abs_diff(opciones.longitud_optima) as f64
                + exceso(propio.mas_estable.energia, -9.0)
                + exceso(propio.extremo_3.energia, -5.0)
                + exceso(pliegue.energia, -3.0);
            if let Some(genoma) = genoma {
                let extremo = &secuencia[longitud - 15..];
                let inverso = complemento_inverso(extremo);
                let sitios = genoma
                    .windows(15)
                    .filter(|v| *v == extremo || *v == inverso.as_slice())
                    .count();
                penalizacion += 5.0 * sitios.saturating_sub(1) as f64;
            }
            lista.push(Candidato {
                inicio,
                longitud,
                secuencia,
                tm,
                penalizacion,
            });
        }
    }
    lista
}

fn mejores_a_mano(
    plantilla: &[u8],
    opciones: &OpcionesDiseno,
    genoma: Option<&[u8]>,
) -> Vec<(f64, usize, usize, usize, usize)> {
    let izquierdos = candidatos_a_mano(plantilla, opciones, true, genoma);
    let derechos = candidatos_a_mano(plantilla, opciones, false, genoma);
    let mut pares = Vec::new();
    for izquierdo in &izquierdos {
        for derecho in &derechos {
            let fin_izquierdo = izquierdo.inicio + izquierdo.longitud;
            if derecho.inicio < fin_izquierdo {
                continue;
            }
            let producto = derecho.inicio + derecho.longitud - izquierdo.inicio;
            if producto < opciones.producto_minimo || producto > opciones.producto_maximo {
                continue;
            }
            if let Some((inicio, longitud)) = opciones.objetivo
                && (fin_izquierdo > inicio || derecho.inicio < inicio + longitud)
            {
                continue;
            }
            let diferencia = (izquierdo.tm - derecho.tm).abs();
            if diferencia > opciones.diferencia_tm_maxima {
                continue;
            }
            let cruzado = dimero(
                &izquierdo.secuencia,
                &derecho.secuencia,
                &opciones.condiciones,
            )
            .unwrap();
            let suma = izquierdo.penalizacion
                + derecho.penalizacion
                + diferencia
                + exceso(cruzado.mas_estable.energia, -9.0)
                + exceso(cruzado.extremo_3.energia, -5.0);
            let penalizacion = (suma * 10_000.0).round() / 10_000.0;
            pares.push((
                penalizacion,
                izquierdo.inicio,
                izquierdo.longitud,
                derecho.inicio,
                derecho.longitud,
            ));
        }
    }
    pares.sort_by(|a, b| {
        a.0.total_cmp(&b.0)
            .then(a.1.cmp(&b.1))
            .then(a.2.cmp(&b.2))
            .then(a.3.cmp(&b.3))
            .then(a.4.cmp(&b.4))
    });
    pares.truncate(opciones.pares);
    pares
}

fn comparar_con_busqueda_completa(
    plantilla: &[u8],
    opciones: &OpcionesDiseno,
    genoma: Option<&[u8]>,
) {
    let diseno = disenar_cebadores(plantilla, opciones, genoma).unwrap();
    let esperados = mejores_a_mano(plantilla, opciones, genoma);
    assert!(!esperados.is_empty());
    assert_eq!(diseno.pares.len(), esperados.len());
    for (par, esperado) in diseno.pares.iter().zip(&esperados) {
        assert!(
            cerca(par.penalizacion, esperado.0, 1e-6),
            "{} {}",
            par.penalizacion,
            esperado.0
        );
        assert_eq!(
            (
                par.izquierdo.inicio,
                par.izquierdo.longitud,
                par.derecho.inicio,
                par.derecho.longitud
            ),
            (esperado.1, esperado.2, esperado.3, esperado.4)
        );
    }
}

#[test]
fn el_diseno_coincide_con_la_busqueda_completa() {
    let plantilla = adn_al_azar(420, 3);
    let opciones = OpcionesDiseno::recomendadas();
    comparar_con_busqueda_completa(&plantilla, &opciones, None);
    let estrecho = OpcionesDiseno {
        producto_minimo: 150,
        producto_maximo: 200,
        pares: 8,
        ..OpcionesDiseno::recomendadas()
    };
    comparar_con_busqueda_completa(&plantilla, &estrecho, None);
}

#[test]
fn el_diseno_con_objetivo_coincide_con_la_busqueda_completa() {
    let plantilla = adn_al_azar(500, 21);
    let opciones = OpcionesDiseno {
        objetivo: Some((200, 40)),
        ..OpcionesDiseno::recomendadas()
    };
    comparar_con_busqueda_completa(&plantilla, &opciones, None);
    let diseno = disenar_cebadores(&plantilla, &opciones, None).unwrap();
    for par in &diseno.pares {
        assert!(par.izquierdo.inicio + par.izquierdo.longitud <= 200);
        assert!(par.derecho.inicio >= 240);
    }
}

#[test]
fn el_diseno_penaliza_cebadores_que_se_repiten_en_el_genoma() {
    let plantilla = adn_al_azar(400, 8);
    let mut genoma = adn_al_azar(5_000, 9);
    genoma.extend(&plantilla);
    genoma.extend(adn_al_azar(3_000, 10));
    let opciones = OpcionesDiseno::recomendadas();
    comparar_con_busqueda_completa(&plantilla, &opciones, Some(&genoma));
    let unico = disenar_cebadores(&plantilla, &opciones, Some(&genoma)).unwrap();
    for par in &unico.pares {
        assert_eq!(par.izquierdo.sitios, Some(1));
        assert_eq!(par.derecho.sitios, Some(1));
    }
    let mut doble = genoma.clone();
    doble.extend(complemento_inverso(&plantilla[..200]));
    comparar_con_busqueda_completa(&plantilla, &opciones, Some(&doble));
    let repetido = disenar_cebadores(&plantilla, &opciones, Some(&doble)).unwrap();
    let sin_genoma = disenar_cebadores(&plantilla, &opciones, None).unwrap();
    assert!(
        sin_genoma
            .pares
            .iter()
            .all(|par| par.izquierdo.sitios.is_none())
    );
    assert!(
        repetido
            .pares
            .iter()
            .all(|par| par.izquierdo.inicio + par.izquierdo.longitud > 200
                || par.izquierdo.sitios == Some(2))
    );
}

#[test]
fn cada_par_cumple_todas_las_reglas() {
    let plantilla = adn_al_azar(3_000, 4);
    let opciones = OpcionesDiseno::recomendadas();
    let diseno = disenar_cebadores(&plantilla, &opciones, None).unwrap();
    assert_eq!(diseno.pares.len(), 5);
    assert!(diseno.izquierdos > 0 && diseno.derechos > 0);
    for par in &diseno.pares {
        for (cebador, directo) in [(&par.izquierdo, true), (&par.derecho, false)] {
            let tramo = &plantilla[cebador.inicio..cebador.inicio + cebador.longitud];
            let esperado = if directo {
                tramo.to_vec()
            } else {
                complemento_inverso(tramo)
            };
            assert_eq!(cebador.secuencia, esperado);
            assert!((18..=25).contains(&cebador.longitud));
            assert!(cebador.tm >= 57.0 && cebador.tm <= 63.0);
            assert!(cebador.gc >= 0.4 && cebador.gc <= 0.6);
            assert!(racha_maxima(&cebador.secuencia) <= 4);
            let pinza = gc_de(&cebador.secuencia[cebador.longitud - 5..]);
            assert!((1..=3).contains(&pinza));
            let tm = temperatura_fusion(&cebador.secuencia, &opciones.condiciones).unwrap();
            assert!(cerca(cebador.tm, tm, 1e-6));
        }
        assert_eq!(
            par.producto,
            par.derecho.inicio + par.derecho.longitud - par.izquierdo.inicio
        );
        assert!((100..=300).contains(&par.producto));
        assert!(par.diferencia_tm <= 3.0);
    }
    assert!(
        diseno
            .pares
            .windows(2)
            .all(|dos| dos[0].penalizacion <= dos[1].penalizacion)
    );
    let otra_vez = disenar_cebadores(&plantilla, &opciones, None).unwrap();
    assert_eq!(diseno, otra_vez);
}

#[test]
fn ningun_cebador_pisa_bases_ambiguas() {
    let mut plantilla = adn_al_azar(1_200, 6);
    for posicion in (0..plantilla.len()).step_by(97) {
        plantilla[posicion] = b'N';
    }
    let diseno = disenar_cebadores(&plantilla, &OpcionesDiseno::recomendadas(), None).unwrap();
    assert!(diseno.descartes.bases_ambiguas > 0);
    for par in &diseno.pares {
        for cebador in [&par.izquierdo, &par.derecho] {
            let tramo = &plantilla[cebador.inicio..cebador.inicio + cebador.longitud];
            assert!(!tramo.contains(&b'N'));
        }
    }
}

#[test]
fn sin_candidatos_no_hay_pares_pero_si_motivos() {
    let plantilla = vec![b'A'; 500];
    let diseno = disenar_cebadores(&plantilla, &OpcionesDiseno::recomendadas(), None).unwrap();
    assert!(diseno.pares.is_empty());
    assert_eq!(diseno.izquierdos, 0);
    assert!(diseno.descartes.repeticiones > 0);
    let rica_en_at = b"ATTACGATTAAT".repeat(50);
    let pobre = disenar_cebadores(&rica_en_at, &OpcionesDiseno::recomendadas(), None).unwrap();
    assert!(pobre.pares.is_empty());
    assert_eq!(pobre.izquierdos, 0);
    assert!(pobre.descartes.contenido_gc > 0);
}

#[test]
fn el_diseno_rechaza_opciones_imposibles() {
    let plantilla = adn_al_azar(1_000, 1);
    let normal = OpcionesDiseno::recomendadas();
    assert!(disenar_cebadores(b"", &normal, None).is_err());
    assert!(disenar_cebadores(&plantilla[..80], &normal, None).is_err());
    assert!(disenar_cebadores(&adn_al_azar(20_001, 1), &normal, None).is_err());
    let malas = [
        OpcionesDiseno {
            longitud_minima: 26,
            ..normal.clone()
        },
        OpcionesDiseno {
            longitud_minima: 10,
            ..normal.clone()
        },
        OpcionesDiseno {
            longitud_maxima: 61,
            ..normal.clone()
        },
        OpcionesDiseno {
            tm_minima: 64.0,
            ..normal.clone()
        },
        OpcionesDiseno {
            tm_optima: f64::NAN,
            ..normal.clone()
        },
        OpcionesDiseno {
            gc_maximo: 1.5,
            ..normal.clone()
        },
        OpcionesDiseno {
            repeticion_maxima: 0,
            ..normal.clone()
        },
        OpcionesDiseno {
            producto_minimo: 400,
            ..normal.clone()
        },
        OpcionesDiseno {
            producto_maximo: 5_001,
            ..normal.clone()
        },
        OpcionesDiseno {
            objetivo: Some((990, 20)),
            ..normal.clone()
        },
        OpcionesDiseno {
            objetivo: Some((100, 0)),
            ..normal.clone()
        },
        OpcionesDiseno {
            objetivo: Some((300, 280)),
            ..normal.clone()
        },
        OpcionesDiseno {
            pares: 0,
            ..normal.clone()
        },
        OpcionesDiseno {
            diferencia_tm_maxima: -1.0,
            ..normal.clone()
        },
        OpcionesDiseno {
            condiciones: Condiciones {
                oligo_nm: 0.0,
                ..Condiciones::default()
            },
            ..normal.clone()
        },
    ];
    for opciones in malas {
        assert!(
            disenar_cebadores(&plantilla, &opciones, None).is_err(),
            "{opciones:?}"
        );
    }
}

#[test]
fn los_pares_vienen_ordenados_y_sin_repetir() {
    let plantilla = adn_al_azar(2_000, 12);
    let opciones = OpcionesDiseno {
        pares: 20,
        ..OpcionesDiseno::recomendadas()
    };
    let diseno = disenar_cebadores(&plantilla, &opciones, None).unwrap();
    assert_eq!(diseno.pares.len(), 20);
    for dos in diseno.pares.windows(2) {
        let orden = dos[0].penalizacion.total_cmp(&dos[1].penalizacion);
        assert_ne!(orden, Ordering::Greater);
        let iguales = dos[0].izquierdo.inicio == dos[1].izquierdo.inicio
            && dos[0].izquierdo.longitud == dos[1].izquierdo.longitud
            && dos[0].derecho.inicio == dos[1].derecho.inicio
            && dos[0].derecho.longitud == dos[1].derecho.longitud;
        assert!(!iguales);
    }
}

#[test]
fn una_zona_repetitiva_se_rechaza_sin_colgarse() {
    let plantilla: Vec<u8> = b"ACGTAC".repeat(300);
    let inicio = std::time::Instant::now();
    let resultado = disenar_cebadores(&plantilla, &OpcionesDiseno::recomendadas(), None);
    assert!(resultado.unwrap_err().contains("repetitiva"));
    assert!(inicio.elapsed().as_secs() < 20);
}

#[test]
fn los_empates_se_deciden_por_la_posicion() {
    let mut opciones = OpcionesDiseno::recomendadas();
    opciones.pares = 20;
    let mut plantilla = adn_al_azar(300, 41);
    plantilla.extend_from_slice(&plantilla.clone());
    let diseno = disenar_cebadores(&plantilla, &opciones, None).unwrap();
    for par in diseno.pares.windows(2) {
        assert_ne!(
            comparar_pares_para_prueba(&par[0], &par[1]),
            Ordering::Greater
        );
    }
}

fn comparar_pares_para_prueba(uno: &ParDeCebadores, otro: &ParDeCebadores) -> Ordering {
    uno.penalizacion
        .total_cmp(&otro.penalizacion)
        .then(uno.izquierdo.inicio.cmp(&otro.izquierdo.inicio))
        .then(uno.izquierdo.longitud.cmp(&otro.izquierdo.longitud))
        .then(uno.derecho.inicio.cmp(&otro.derecho.inicio))
        .then(uno.derecho.longitud.cmp(&otro.derecho.longitud))
}
