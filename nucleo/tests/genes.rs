use nucleo::azar::Azar;
use nucleo::ecoli;
use nucleo::genes::*;
use nucleo::secuencia::complemento_inverso;
use std::sync::OnceLock;

fn modelo_ecoli() -> &'static Modelo {
    static MODELO: OnceLock<Modelo> = OnceLock::new();
    MODELO.get_or_init(|| entrenar_modelo(ecoli::genoma()))
}

fn candidato(desde: usize, hasta: usize, hebra: u8, puntaje: f64) -> GenPredicho {
    GenPredicho {
        desde,
        hasta,
        hebra,
        puntaje,
        ..GenPredicho::default()
    }
}

fn borde(posicion: usize) -> Inicio {
    Inicio {
        posicion,
        codon: "",
        rbs: None,
    }
}

#[test]
fn reconoce_codones_de_parada_y_de_inicio() {
    assert!(es_codon_de_parada(b"TAA"));
    assert!(es_codon_de_parada(b"TAG"));
    assert!(es_codon_de_parada(b"TGA"));
    assert!(!es_codon_de_parada(b"TGG"));
    assert!(!es_codon_de_parada(b"TA"));
    assert_eq!(codon_de_inicio(b"ATG"), Some("ATG"));
    assert_eq!(codon_de_inicio(b"GTG"), Some("GTG"));
    assert_eq!(codon_de_inicio(b"TTG"), Some("TTG"));
    assert_eq!(codon_de_inicio(b"CTG"), None);
    assert_eq!(codon_de_inicio(b"AT"), None);
}

#[test]
fn encuentra_el_motivo_de_shine_dalgarno_mas_largo() {
    assert_eq!(
        buscar_rbs(b"AGGAGGAAAAAAAATG", 13),
        Some(Rbs {
            motivo: "AGGAGG",
            espaciador: 7
        })
    );
    assert_eq!(
        buscar_rbs(b"AGGAGGAAAAATG", 10),
        Some(Rbs {
            motivo: "AGGAG",
            espaciador: 5
        })
    );
    assert_eq!(buscar_rbs(b"CCCCCCCCCCCCCCCCCCCCATG", 20), None);
    assert_eq!(buscar_rbs(b"GGATG", 2), None);
}

#[test]
fn prefiere_el_espaciador_cercano_al_lejano() {
    let mut bases = vec![b'C'; 30];
    bases[27..30].copy_from_slice(b"ATG");
    bases[10..14].copy_from_slice(b"AGGA");
    let mut con_los_dos = bases.clone();
    con_los_dos[15..19].copy_from_slice(b"GGAG");
    assert_eq!(
        buscar_rbs(&con_los_dos, 27),
        Some(Rbs {
            motivo: "GGAG",
            espaciador: 8
        })
    );
    assert_eq!(
        buscar_rbs(&bases, 27),
        Some(Rbs {
            motivo: "AGGA",
            espaciador: 13
        })
    );
}

#[test]
fn clasifica_el_rbs_por_largo_y_distancia() {
    let rbs = |motivo, espaciador| Some(Rbs { motivo, espaciador });
    assert_eq!(clase_rbs(None), CLASES_RBS - 1);
    assert_eq!(clase_rbs(rbs("GGA", 5)), 0);
    assert_eq!(clase_rbs(rbs("GAG", 15)), 1);
    assert_eq!(clase_rbs(rbs("GGAG", 10)), 2);
    assert_eq!(clase_rbs(rbs("GGAG", 11)), 3);
    assert_eq!(clase_rbs(rbs("AGGAGG", 7)), 6);
    assert_eq!(clase_rbs(rbs("AGGAGG", 12)), 7);
}

#[test]
fn numera_los_hexameros_de_cada_posicion() {
    assert_eq!(
        hexameros_por_posicion(b"ACGTAC"),
        vec![Some(433), None, None, None, None, None]
    );
    assert_eq!(
        hexameros_por_posicion(b"acgtacG"),
        vec![Some(433), Some(1734), None, None, None, None, None]
    );
    let con_hueco = hexameros_por_posicion(b"AAAAAANAAAAAA");
    let mut esperado = vec![None; 13];
    esperado[0] = Some(0);
    esperado[7] = Some(0);
    assert_eq!(con_hueco, esperado);
    assert!(hexameros_por_posicion(b"").is_empty());
}

fn ejemplo_de_marcos() -> Vec<u8> {
    let mut adn = b"CCCATG".to_vec();
    for _ in 0..30 {
        adn.extend_from_slice(b"GCT");
    }
    adn.extend_from_slice(b"TAACCC");
    adn
}

#[test]
fn arma_los_marcos_abiertos_de_la_hebra_directa() {
    let hebras = preparar_hebras(&ejemplo_de_marcos());
    assert_eq!(hebras.directa.bases.len(), 102);
    let marcos = marcos_abiertos(&hebras.directa);
    assert_eq!(
        marcos,
        vec![
            Marco {
                hebra: b'+',
                inicios: vec![
                    borde(0),
                    Inicio {
                        posicion: 3,
                        codon: "ATG",
                        rbs: None
                    }
                ],
                fin: 99,
                tiene_parada: true,
            },
            Marco {
                hebra: b'+',
                inicios: vec![borde(1)],
                fin: 100,
                tiene_parada: false,
            },
            Marco {
                hebra: b'+',
                inicios: vec![borde(2)],
                fin: 101,
                tiene_parada: false,
            },
        ]
    );
}

#[test]
fn arma_los_marcos_abiertos_de_la_hebra_inversa() {
    let adn = ejemplo_de_marcos();
    let hebras = preparar_hebras(&adn);
    assert_eq!(hebras.inversa.bases, complemento_inverso(&adn));
    assert_eq!(
        marcos_abiertos(&hebras.inversa),
        vec![
            Marco {
                hebra: b'-',
                inicios: vec![borde(0)],
                fin: 102,
                tiene_parada: false,
            },
            Marco {
                hebra: b'-',
                inicios: vec![borde(2)],
                fin: 101,
                tiene_parada: false,
            },
        ]
    );
}

#[test]
fn sin_marcos_en_secuencias_vacias_o_cortas() {
    for adn in [&b""[..], b"ATG", b"ACGTN", &[b'A'; 89]] {
        let hebras = preparar_hebras(adn);
        assert!(marcos_abiertos(&hebras.directa).is_empty());
        assert!(marcos_abiertos(&hebras.inversa).is_empty());
    }
}

#[test]
fn elige_genes_respetando_los_solapamientos() {
    assert!(elegir_genes(&[]).is_empty());
    let separados = [candidato(0, 300, b'+', 5.0), candidato(400, 700, b'+', 7.0)];
    assert_eq!(elegir_genes(&separados), vec![0, 1]);
    let roce_permitido = [
        candidato(0, 300, b'+', 10.0),
        candidato(240, 600, b'+', 8.0),
    ];
    assert_eq!(elegir_genes(&roce_permitido), vec![0, 1]);
    let roce_excesivo = [
        candidato(0, 300, b'+', 10.0),
        candidato(239, 600, b'+', 8.0),
    ];
    assert_eq!(elegir_genes(&roce_excesivo), vec![0]);
    let colas = [
        candidato(0, 300, b'+', 10.0),
        candidato(100, 500, b'-', 8.0),
    ];
    assert_eq!(elegir_genes(&colas), vec![0, 1]);
    let colas_de_mas = [candidato(0, 300, b'+', 10.0), candidato(99, 500, b'-', 8.0)];
    assert_eq!(elegir_genes(&colas_de_mas), vec![0]);
    let cabezas = [
        candidato(0, 300, b'-', 10.0),
        candidato(299, 600, b'+', 8.0),
    ];
    assert_eq!(elegir_genes(&cabezas), vec![0]);
    let cabezas_justas = [
        candidato(0, 300, b'-', 10.0),
        candidato(300, 600, b'+', 8.0),
    ];
    assert_eq!(elegir_genes(&cabezas_justas), vec![0, 1]);
}

#[test]
fn elige_la_mejor_combinacion_y_desempata_por_el_primero() {
    let contenido = [
        candidato(0, 1000, b'+', 10.0),
        candidato(100, 400, b'-', 3.0),
    ];
    assert_eq!(elegir_genes(&contenido), vec![0]);
    let pareja_contra_uno = [
        candidato(0, 300, b'+', 5.0),
        candidato(250, 600, b'+', 5.0),
        candidato(100, 500, b'+', 8.0),
    ];
    assert_eq!(elegir_genes(&pareja_contra_uno), vec![0, 1]);
    let empate = [candidato(0, 300, b'+', 5.0), candidato(30, 300, b'+', 5.0)];
    assert_eq!(elegir_genes(&empate), vec![0]);
    let negativo = [candidato(0, 300, b'+', -1.0)];
    assert!(elegir_genes(&negativo).is_empty());
}

fn mejor_total_probando_todo(candidatos: &[GenPredicho]) -> f64 {
    let mut mejor: f64 = 0.0;
    for mascara in 0u32..(1 << candidatos.len()) {
        let mut elegidos: Vec<&GenPredicho> = (0..candidatos.len())
            .filter(|indice| mascara & (1 << indice) != 0)
            .map(|indice| &candidatos[indice])
            .collect();
        elegidos.sort_by_key(|c| (c.hasta, c.desde));
        let valido = elegidos
            .windows(2)
            .all(|par| pueden_ir_seguidos(par[0], par[1]));
        if valido {
            mejor = mejor.max(elegidos.iter().map(|c| c.puntaje).sum());
        }
    }
    mejor
}

#[test]
fn la_programacion_dinamica_iguala_a_probar_todas_las_combinaciones() {
    let mut azar = Azar::nuevo(7);
    for _ in 0..300 {
        let cantidad = 1 + azar.entero(10);
        let candidatos: Vec<GenPredicho> = (0..cantidad)
            .map(|_| {
                let desde = azar.entero(1500);
                let largo = 90 + azar.entero(600);
                let hebra = if azar.entero(2) == 0 { b'+' } else { b'-' };
                let puntaje = azar.entero(40) as f64 - 5.0;
                candidato(desde, desde + largo, hebra, puntaje)
            })
            .collect();
        let elegidos = elegir_genes(&candidatos);
        let total: f64 = elegidos.iter().map(|&i| candidatos[i].puntaje).sum();
        for par in elegidos.windows(2) {
            assert!(pueden_ir_seguidos(&candidatos[par[0]], &candidatos[par[1]]));
        }
        assert!((total - mejor_total_probando_todo(&candidatos)).abs() < 1e-9);
    }
}

#[test]
fn encuentra_los_genes_del_operon_thr_de_ecoli() {
    let region = ecoli::region(1, 6000).unwrap();
    let prediccion = predecir_genes(region, modelo_ecoli());
    assert!(!prediccion.entrenado_con_la_entrada);
    let completos: Vec<(usize, usize, u8, &str)> = prediccion
        .genes
        .iter()
        .filter(|g| g.completo())
        .map(|g| (g.desde, g.hasta, g.hebra, g.codon_de_inicio))
        .collect();
    for esperado in [
        (336, 2799, b'+', "ATG"),
        (2800, 3733, b'+', "ATG"),
        (3733, 5020, b'+', "ATG"),
    ] {
        assert!(completos.contains(&esperado), "falta {esperado:?}");
    }
}

#[test]
fn encuentra_genes_en_la_hebra_inversa() {
    let region = ecoli::region(5001, 9000).unwrap();
    let prediccion = predecir_genes(region, modelo_ecoli());
    let inversos: Vec<(usize, usize)> = prediccion
        .genes
        .iter()
        .filter(|g| g.hebra == b'-' && g.completo())
        .map(|g| (g.desde, g.hasta))
        .collect();
    assert!(inversos.contains(&(5683 - 5001, 6459 - 5000)));
    assert!(inversos.contains(&(6529 - 5001, 7959 - 5000)));
}

#[test]
fn marca_como_parciales_los_genes_cortados_por_los_bordes() {
    let region = ecoli::region(1001, 4000).unwrap();
    let genes = predecir_genes(region, modelo_ecoli()).genes;
    let primero = genes[0];
    assert_eq!((primero.desde, primero.hasta), (2, 2799 - 1000));
    assert!(primero.parcial_al_inicio && !primero.parcial_al_final);
    assert_eq!(primero.codon_de_inicio, "");
    let ultimo = genes[genes.len() - 1];
    assert_eq!((ultimo.desde, ultimo.hasta), (3734 - 1001, 3000));
    assert!(ultimo.parcial_al_final && !ultimo.parcial_al_inicio);
}

fn revisar_genes(adn: &[u8], genes: &[GenPredicho]) {
    for g in genes {
        assert!(g.hasta <= adn.len());
        assert!(g.longitud() >= LONGITUD_MINIMA);
        assert_eq!(g.longitud() % 3, 0);
        let tramo = &adn[g.desde..g.hasta];
        let en_su_hebra = if g.hebra == b'-' {
            complemento_inverso(tramo)
        } else {
            tramo.to_vec()
        };
        let codones: Vec<&[u8]> = en_su_hebra.chunks(3).collect();
        if !g.parcial_al_inicio {
            assert_eq!(codones[0], g.codon_de_inicio.as_bytes());
        }
        let ultimo = codones.len() - 1;
        assert_eq!(es_codon_de_parada(codones[ultimo]), !g.parcial_al_final);
        for codon in &codones[..ultimo] {
            assert!(!es_codon_de_parada(codon));
        }
        if let Some(rbs) = g.rbs {
            assert!((5..=15).contains(&rbs.espaciador));
        }
    }
}

#[test]
fn cada_gen_va_de_un_inicio_a_una_parada_sin_paradas_internas() {
    let region = ecoli::region(200_001, 260_000).unwrap();
    let prediccion = predecir_genes(region, modelo_ecoli());
    assert!(prediccion.entrenado_con_la_entrada);
    assert!(prediccion.genes.len() > 40);
    revisar_genes(region, &prediccion.genes);
    let mut azar = Azar::nuevo(3);
    let al_azar: Vec<u8> = (0..30_000).map(|_| b"ACGT"[azar.entero(4)]).collect();
    let prediccion = predecir_genes(&al_azar, modelo_ecoli());
    revisar_genes(&al_azar, &prediccion.genes);
}

#[test]
fn no_se_rompe_con_entradas_raras() {
    let neutro = Modelo::neutro(0.5);
    assert!(predecir_genes(b"", &neutro).genes.is_empty());
    assert!(predecir_genes(b"ACGTN", modelo_ecoli()).genes.is_empty());
    assert!(
        predecir_genes(&[b'N'; 500], modelo_ecoli())
            .genes
            .is_empty()
    );
    let solo_a = vec![b'A'; 30_000];
    let prediccion = predecir_genes(&solo_a, modelo_ecoli());
    assert!(prediccion.entrenado_con_la_entrada);
    assert!(prediccion.genes.is_empty());
    let vacio = entrenar_modelo(b"");
    assert!(vacio.codificacion.iter().all(|valor| valor.is_finite()));
    let region = ecoli::region(1, 6000).unwrap();
    let minusculas = region.to_ascii_lowercase();
    assert_eq!(
        predecir_genes(&minusculas, modelo_ecoli()).genes,
        predecir_genes(region, modelo_ecoli()).genes
    );
}

#[test]
fn compara_con_la_anotacion_por_codon_de_parada() {
    let anotados = [
        Anotado {
            desde: 0,
            hasta: 300,
            hebra: b'+',
        },
        Anotado {
            desde: 400,
            hasta: 700,
            hebra: b'-',
        },
        Anotado {
            desde: 800,
            hasta: 900,
            hebra: b'+',
        },
    ];
    let mut cortado = candidato(1500, 1800, b'+', 1.0);
    cortado.parcial_al_final = true;
    let predichos = [
        candidato(0, 300, b'+', 1.0),
        candidato(400, 650, b'-', 1.0),
        candidato(1000, 1300, b'+', 1.0),
        candidato(790, 900, b'+', 1.0),
        cortado,
    ];
    let comparacion = comparar_con_anotacion(&predichos, &anotados);
    assert_eq!(
        comparacion,
        Comparacion {
            anotados: 3,
            encontrados: 3,
            con_inicio_exacto: 1,
            anotados_largos: 2,
            largos_encontrados: 2,
            largos_con_inicio_exacto: 1,
            predichos_completos: 4,
            predichos_sin_anotacion: 1,
            pareja: vec![Some(0), Some(1), None, Some(2), None],
            hallado: vec![true, true, true],
        }
    );
}

#[test]
fn escribe_gff3() {
    let mut primero = candidato(0, 300, b'+', 12.345);
    primero.codon_de_inicio = "ATG";
    primero.rbs = Some(Rbs {
        motivo: "AGGAG",
        espaciador: 7,
    });
    let mut segundo = candidato(400, 700, b'-', 3.0);
    segundo.parcial_al_inicio = true;
    let mut texto = Vec::new();
    escribir_gff(&mut texto, "prueba", 1000, &[primero, segundo]).unwrap();
    assert_eq!(
        String::from_utf8(texto).unwrap(),
        "##gff-version 3\n\
         ##sequence-region prueba 1 1000\n\
         prueba\tLaCompania\tCDS\t1\t300\t12.3\t+\t0\tID=gen_1;partial=00;start_type=ATG;rbs_motif=AGGAG;rbs_spacer=7\n\
         prueba\tLaCompania\tCDS\t401\t700\t3.0\t-\t0\tID=gen_2;partial=01;start_type=borde;rbs_motif=ninguno;rbs_spacer=-\n"
    );
}

#[test]
fn alcanza_la_precision_esperada_en_todo_el_genoma_de_ecoli() {
    let genoma = ecoli::genoma();
    let prediccion = predecir_genes(genoma, modelo_ecoli());
    assert!(prediccion.entrenado_con_la_entrada);
    let anotados: Vec<Anotado> = ecoli::genes()
        .iter()
        .map(|g| Anotado {
            desde: g.inicio - 1,
            hasta: g.fin,
            hebra: g.hebra,
        })
        .collect();
    let comparacion = comparar_con_anotacion(&prediccion.genes, &anotados);
    let largos = comparacion.anotados_largos as f64;
    assert!(comparacion.largos_encontrados as f64 / largos >= 0.90);
    assert!(comparacion.largos_con_inicio_exacto as f64 / largos >= 0.70);
    let extras = comparacion.predichos_sin_anotacion as f64;
    assert!(extras / comparacion.predichos_completos as f64 <= 0.10);
    revisar_genes(genoma, &prediccion.genes);
}

#[test]
fn buscar_rbs_no_se_sale_de_la_secuencia() {
    assert_eq!(buscar_rbs(b"AGGAGGAAAAATG", 100), None);
    assert_eq!(buscar_rbs(b"", 0), None);
    assert_eq!(buscar_rbs(b"AGG", 3), None);
}

#[test]
fn clase_rbs_acepta_motivos_de_cualquier_largo() {
    let corto = Rbs {
        motivo: "AG",
        espaciador: 5,
    };
    let largo = Rbs {
        motivo: "AGGAGGAGG",
        espaciador: 12,
    };
    assert_eq!(clase_rbs(Some(corto)), 0);
    assert_eq!(clase_rbs(Some(largo)), 7);
}

#[test]
fn la_comparacion_cuenta_todos_los_anotados_con_la_misma_parada() {
    let anotados = [
        Anotado {
            desde: 0,
            hasta: 600,
            hebra: b'+',
        },
        Anotado {
            desde: 150,
            hasta: 600,
            hebra: b'+',
        },
        Anotado {
            desde: 1000,
            hasta: 1450,
            hebra: b'-',
        },
        Anotado {
            desde: 1000,
            hasta: 1300,
            hebra: b'-',
        },
    ];
    let predichos = [
        candidato(150, 600, b'+', 1.0),
        candidato(1000, 1300, b'-', 1.0),
    ];
    let comparacion = comparar_con_anotacion(&predichos, &anotados);
    assert_eq!(comparacion.anotados, 4);
    assert_eq!(comparacion.encontrados, 4);
    assert_eq!(comparacion.con_inicio_exacto, 2);
    assert_eq!(comparacion.hallado, vec![true, true, true, true]);
    assert_eq!(comparacion.pareja, vec![Some(1), Some(3)]);
    assert_eq!(comparacion.predichos_sin_anotacion, 0);
}

#[test]
fn elegir_genes_con_miles_de_inicios_en_el_mismo_marco_no_se_cuelga() {
    let mut candidatos = Vec::new();
    for numero in 0..60_000usize {
        let desde = numero * 3;
        candidatos.push(candidato(desde, 200_000, b'+', numero as f64 * 0.001 + 1.0));
    }
    for numero in 0..2_000usize {
        let desde = 200_000 + numero * 3;
        candidatos.push(candidato(desde, 260_000, b'-', 2.0));
    }
    let elegidos = elegir_genes(&candidatos);
    assert_eq!(elegidos.len(), 2);
    assert_eq!(elegidos[0], 59_999);
}
