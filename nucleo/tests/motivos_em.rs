use nucleo::azar::Azar;
use nucleo::motivos_em::*;
use nucleo::secuencia::complemento_inverso;

fn opciones(k: usize, modelo: Modelo, ambas_hebras: bool) -> Opciones {
    Opciones {
        k,
        modelo,
        ambas_hebras,
        pseudoconteo: 1.0,
        arranques: 1000,
    }
}

fn cerca(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9 * (1.0 + a.abs().max(b.abs()))
}

fn al_azar(azar: &mut Azar, largo: usize) -> Vec<u8> {
    (0..largo).map(|_| b"ACGT"[azar.entero(4)]).collect()
}

fn con_motivo(
    azar: &mut Azar,
    motivo: &[u8],
    cuantas: usize,
    largo: usize,
) -> (Vec<Vec<u8>>, Vec<usize>) {
    let mut textos = Vec::new();
    let mut posiciones = Vec::new();
    for _ in 0..cuantas {
        let mut texto = al_azar(azar, largo);
        let posicion = azar.entero(largo - motivo.len() + 1);
        texto[posicion..posicion + motivo.len()].copy_from_slice(motivo);
        textos.push(texto);
        posiciones.push(posicion);
    }
    (textos, posiciones)
}

fn ingenuo(
    codigos: &[Vec<u8>],
    parametros: &Parametros,
    fondo: &[f64; 4],
    opciones: &Opciones,
) -> (f64, Vec<Vec<f64>>) {
    let k = opciones.k;
    let mut verosimilitud = 0.0;
    let mut todos = Vec::new();
    for fila in codigos {
        let solo_fondo: f64 = fila.iter().map(|&c| fondo[c as usize]).product();
        let mut con_sitio = Vec::new();
        for inicio in 0..=fila.len() - k {
            let mut hebras = vec![fila[inicio..inicio + k].to_vec()];
            if opciones.ambas_hebras {
                hebras.push(
                    fila[inicio..inicio + k]
                        .iter()
                        .rev()
                        .map(|&c| 3 - c)
                        .collect(),
                );
            }
            for leida in hebras {
                let mut probabilidad = 1.0;
                for (posicion, &c) in fila.iter().enumerate() {
                    if posicion >= inicio && posicion < inicio + k {
                        probabilidad *=
                            parametros.perfil[posicion - inicio][leida[posicion - inicio] as usize];
                    } else {
                        probabilidad *= fondo[c as usize];
                    }
                }
                con_sitio.push(probabilidad);
            }
        }
        let candidatos = con_sitio.len() as f64;
        let lambda = parametros.gamma / candidatos;
        let mut total: f64 = con_sitio.iter().map(|p| lambda * p).sum();
        if opciones.modelo == Modelo::Zoops {
            total += (1.0 - parametros.gamma) * solo_fondo;
        }
        verosimilitud += (total / solo_fondo).ln();
        todos.push(con_sitio.iter().map(|p| lambda * p / total).collect());
    }
    (verosimilitud, todos)
}

fn perfil_al_azar(azar: &mut Azar, k: usize) -> Vec<[f64; 4]> {
    (0..k)
        .map(|_| {
            let mut columna = [0.0; 4];
            for valor in columna.iter_mut() {
                *valor = 0.05 + azar.decimal();
            }
            let suma: f64 = columna.iter().sum();
            columna.map(|valor| valor / suma)
        })
        .collect()
}

#[test]
fn fondo_cuenta_bases_con_un_pseudoconteo() {
    let codigos = codificar(&[b"AAAC".to_vec(), b"GT".to_vec()]).unwrap();
    let fondo = frecuencias_de_fondo(&codigos, false);
    assert!(cerca(fondo[0], 4.0 / 10.0));
    assert!(cerca(fondo[1], 2.0 / 10.0));
    assert!(cerca(fondo[2], 2.0 / 10.0));
    assert!(cerca(fondo[3], 2.0 / 10.0));
    let simetrico = frecuencias_de_fondo(&codigos, true);
    assert!(cerca(simetrico[0], 3.0 / 10.0));
    assert!(cerca(simetrico[3], 3.0 / 10.0));
    assert!(cerca(simetrico[1], 2.0 / 10.0));
}

#[test]
fn posteriores_oops_y_zoops_a_mano() {
    let mut z = vec![2f64.ln(), 0.0];
    let ln_total = posteriores_de_secuencia(&mut z, 1.0, Modelo::Oops);
    assert!(cerca(z[0], 2.0 / 3.0) && cerca(z[1], 1.0 / 3.0));
    assert!(cerca(ln_total, 1.5f64.ln()));
    let mut z = vec![2f64.ln(), 0.0];
    let ln_total = posteriores_de_secuencia(&mut z, 0.5, Modelo::Zoops);
    assert!(cerca(z[0], 0.4) && cerca(z[1], 0.2));
    assert!(cerca(ln_total, 1.25f64.ln()));
    let mut z = vec![0.0, 0.0];
    posteriores_de_secuencia(&mut z, 1.0, Modelo::Zoops);
    assert!(cerca(z[0], 0.5) && cerca(z[1], 0.5));
}

#[test]
fn las_dos_hebras_puntuan_igual_un_palindromo() {
    let codigos = codificar(&[b"ACGT".to_vec()]).unwrap();
    let perfil = vec![
        [0.7, 0.1, 0.1, 0.1],
        [0.1, 0.7, 0.1, 0.1],
        [0.1, 0.1, 0.7, 0.1],
        [0.1, 0.1, 0.1, 0.7],
    ];
    let pesos = pesos_logaritmicos(&perfil, &[0.25; 4]);
    let mut puntajes = Vec::new();
    puntajes_de_ventanas(&pesos, &codigos[0], true, &mut puntajes);
    assert_eq!(puntajes.len(), 2);
    assert!(cerca(puntajes[0], puntajes[1]));
    assert!(cerca(puntajes[0], 4.0 * 2.8f64.ln()));
    let sitio = mejor_sitio(&[0.25, 0.25], true);
    assert_eq!((sitio.posicion, sitio.hebra), (0, b'+'));
    assert!(cerca(sitio.probabilidad, 0.5));
}

#[test]
fn paso_m_con_posteriores_seguros_cuenta_las_letras() {
    let codigos = codificar(&[b"AAC".to_vec(), b"GAC".to_vec()]).unwrap();
    let opciones = opciones(2, Modelo::Oops, false);
    let posteriores = vec![vec![0.0, 1.0], vec![0.0, 1.0]];
    let fondo = [0.25; 4];
    let nuevos = paso_m(&codigos, &posteriores, &fondo, &opciones);
    assert!(cerca(nuevos.perfil[0][0], 2.25 / 3.0));
    assert!(cerca(nuevos.perfil[0][1], 0.25 / 3.0));
    assert!(cerca(nuevos.perfil[1][1], 2.25 / 3.0));
    assert_eq!(nuevos.gamma, 1.0);
    let ambas = Opciones {
        ambas_hebras: true,
        ..opciones
    };
    let posteriores = vec![vec![0.0, 0.0, 0.0, 1.0], vec![0.0, 0.0, 0.0, 1.0]];
    let reverso = paso_m(&codigos, &posteriores, &fondo, &ambas);
    assert!(cerca(reverso.perfil[0][2], 2.25 / 3.0));
    assert!(cerca(reverso.perfil[1][3], 2.25 / 3.0));
}

#[test]
fn informacion_y_consenso_de_columnas_extremas() {
    let perfil = vec![[1.0, 0.0, 0.0, 0.0], [0.25; 4], [0.0, 0.0, 0.5, 0.5]];
    let bits = informacion_por_columna(&perfil);
    assert!(cerca(bits[0], 2.0) && cerca(bits[1], 0.0) && cerca(bits[2], 1.0));
    assert_eq!(consenso_del_perfil(&perfil), b"AAG".to_vec());
}

#[test]
fn arranques_son_kmers_distintos_en_orden() {
    let codigos = codificar(&[b"AAAA".to_vec(), b"AACA".to_vec()]).unwrap();
    let semillas = arranques(&codigos, 2, 10);
    let texto: Vec<Vec<u8>> = semillas.iter().map(|s| decodificar(s)).collect();
    assert_eq!(texto, vec![b"AA".to_vec(), b"AC".to_vec(), b"CA".to_vec()]);
    assert_eq!(arranques(&codigos, 2, 2).len(), 2);
}

#[test]
fn rechaza_entradas_malas() {
    let bien = vec![b"ACGTACGT".to_vec(), b"TTTTACGT".to_vec()];
    assert!(meme(&[], &opciones(3, Modelo::Oops, false)).is_err());
    assert!(meme(&bien, &opciones(0, Modelo::Oops, false)).is_err());
    assert!(meme(&bien, &opciones(9, Modelo::Oops, false)).is_err());
    assert!(meme(&[b"ACGN".to_vec()], &opciones(2, Modelo::Oops, false)).is_err());
    let sin_pseudoconteo = Opciones {
        pseudoconteo: 0.0,
        ..opciones(3, Modelo::Oops, false)
    };
    assert!(meme(&bien, &sin_pseudoconteo).is_err());
    let sin_arranques = Opciones {
        arranques: 0,
        ..opciones(3, Modelo::Oops, false)
    };
    assert!(meme(&bien, &sin_arranques).is_err());
    let corto = meme(
        &[b"ACG".to_vec(), b"ACG".to_vec()],
        &opciones(3, Modelo::Zoops, true),
    )
    .unwrap();
    assert_eq!(corto.sitios.len(), 2);
}

#[test]
fn oops_encuentra_un_motivo_plantado() {
    let mut azar = Azar::nuevo(11);
    let motivo = b"TTGACAGCTAGC";
    let (textos, posiciones) = con_motivo(&mut azar, motivo, 8, 120);
    let hallado = meme(&textos, &opciones(motivo.len(), Modelo::Oops, false)).unwrap();
    assert_eq!(hallado.consenso, motivo.to_vec());
    for (sitio, &posicion) in hallado.sitios.iter().zip(&posiciones) {
        assert_eq!(sitio.posicion, posicion);
        assert!(sitio.probabilidad > 0.5);
        assert!(cerca(sitio.probabilidad_de_sitio, 1.0));
    }
    assert!(hallado.informacion.iter().all(|&bits| bits > 1.0));
    assert_eq!(hallado.historial.len(), hallado.iteraciones + 1);
    assert!(hallado.historial.last().unwrap() > hallado.historial.first().unwrap());
}

#[test]
fn zoops_deja_sin_sitio_a_las_secuencias_sin_motivo() {
    let mut azar = Azar::nuevo(5);
    let motivo = b"GATCCGGTTACA";
    let (mut textos, _) = con_motivo(&mut azar, motivo, 15, 100);
    for _ in 0..5 {
        textos.push(al_azar(&mut azar, 100));
    }
    let hallado = meme(&textos, &opciones(motivo.len(), Modelo::Zoops, false)).unwrap();
    assert_eq!(hallado.consenso, motivo.to_vec());
    for sitio in &hallado.sitios[..15] {
        assert!(sitio.probabilidad_de_sitio > 0.5);
    }
    let sin_sitio = hallado.sitios[15..]
        .iter()
        .filter(|sitio| sitio.probabilidad_de_sitio < 0.5)
        .count();
    assert!(sin_sitio >= 4);
    assert!((hallado.gamma - 0.75).abs() < 0.06);
}

#[test]
fn ambas_hebras_encuentra_sitios_al_reves() {
    let mut azar = Azar::nuevo(3);
    let motivo = b"AACGGTTGCAAT";
    let (mut textos, posiciones) = con_motivo(&mut azar, motivo, 8, 90);
    for texto in textos.iter_mut().skip(4) {
        *texto = complemento_inverso(texto);
    }
    let hallado = meme(&textos, &opciones(motivo.len(), Modelo::Oops, true)).unwrap();
    let directo = hallado.consenso == motivo.to_vec();
    assert!(directo || hallado.consenso == complemento_inverso(motivo));
    for (indice, sitio) in hallado.sitios.iter().enumerate() {
        let esperada = if indice < 4 {
            posiciones[indice]
        } else {
            90 - motivo.len() - posiciones[indice]
        };
        assert_eq!(sitio.posicion, esperada);
        let al_reves = (indice >= 4) == directo;
        assert_eq!(sitio.hebra, if al_reves { b'-' } else { b'+' });
    }
    let una_hebra = meme(&textos, &opciones(motivo.len(), Modelo::Oops, false)).unwrap();
    assert!(hallado.verosimilitud > una_hebra.verosimilitud);
}

#[test]
fn paso_e_coincide_con_la_probabilidad_completa() {
    let mut azar = Azar::nuevo(42);
    for caso in 0..200 {
        let k = 1 + azar.entero(5);
        let modelo = if caso % 2 == 0 {
            Modelo::Oops
        } else {
            Modelo::Zoops
        };
        let ambas = caso % 3 == 0;
        let opciones = opciones(k, modelo, ambas);
        let textos: Vec<Vec<u8>> = (0..1 + azar.entero(4))
            .map(|_| {
                let largo = k + azar.entero(12);
                al_azar(&mut azar, largo)
            })
            .collect();
        let codigos = codificar(&textos).unwrap();
        let fondo = frecuencias_de_fondo(&codigos, ambas);
        let parametros = Parametros {
            perfil: perfil_al_azar(&mut azar, k),
            gamma: if modelo == Modelo::Oops {
                1.0
            } else {
                0.05 + 0.9 * azar.decimal()
            },
        };
        let mut posteriores = posteriores_vacios(&codigos, &opciones);
        let rapido = paso_e(&codigos, &parametros, &fondo, &opciones, &mut posteriores);
        let (lento, esperados) = ingenuo(&codigos, &parametros, &fondo, &opciones);
        assert!(
            cerca(rapido, lento),
            "caso {caso}: {rapido} frente a {lento}"
        );
        for (fila, esperada) in posteriores.iter().zip(&esperados) {
            assert_eq!(fila.len(), esperada.len());
            for (a, b) in fila.iter().zip(esperada) {
                assert!((a - b).abs() < 1e-9, "caso {caso}");
            }
        }
    }
}

#[test]
fn em_nunca_baja_el_objetivo_con_su_prior() {
    let mut azar = Azar::nuevo(9);
    for caso in 0..30 {
        let modelo = if caso % 2 == 0 {
            Modelo::Oops
        } else {
            Modelo::Zoops
        };
        let opciones = opciones(6, modelo, caso % 3 == 0);
        let (textos, _) = con_motivo(&mut azar, b"ACGGTA", 5, 40);
        let codigos = codificar(&textos).unwrap();
        let fondo = frecuencias_de_fondo(&codigos, opciones.ambas_hebras);
        let semilla = &codigos[0][azar.entero(35)..][..6];
        let mut parametros = parametros_de_arranque(semilla, modelo);
        let mut posteriores = posteriores_vacios(&codigos, &opciones);
        let mut anterior = f64::NEG_INFINITY;
        for _ in 0..25 {
            let verosimilitud = paso_e(&codigos, &parametros, &fondo, &opciones, &mut posteriores);
            let mut prior = 0.0;
            for columna in &parametros.perfil {
                for base in 0..4 {
                    prior += fondo[base] * columna[base].ln();
                }
            }
            let objetivo = verosimilitud + prior;
            assert!(
                objetivo >= anterior - 1e-9,
                "caso {caso}: {objetivo} < {anterior}"
            );
            anterior = objetivo;
            parametros = paso_m(&codigos, &posteriores, &fondo, &opciones);
        }
    }
}

#[test]
fn meme_elige_el_arranque_de_mejor_prueba() {
    let mut azar = Azar::nuevo(21);
    let (textos, _) = con_motivo(&mut azar, b"CCATGGAT", 5, 50);
    let opciones = opciones(8, Modelo::Oops, false);
    let hallado = meme(&textos, &opciones).unwrap();
    let codigos = codificar(&textos).unwrap();
    let fondo = frecuencias_de_fondo(&codigos, false);
    let semillas = arranques(&codigos, 8, 1000);
    let puntajes: Vec<f64> = semillas
        .iter()
        .map(|semilla| probar_arranque(&codigos, semilla, &fondo, &opciones))
        .collect();
    let mejor = mejor_arranque(&puntajes);
    assert_eq!(decodificar(semillas[mejor]), hallado.arranque);
    assert_eq!(hallado.arranques_probados, semillas.len());
    let recorrido = iterar(
        &codigos,
        parametros_de_arranque(semillas[mejor], Modelo::Oops),
        &fondo,
        &opciones,
        ITERACIONES_MAXIMAS,
    );
    assert_eq!(recorrido.historial, hallado.historial);
}

#[test]
fn un_pseudoconteo_extremo_se_rechaza_en_vez_de_dar_nan() {
    let textos = vec![b"AAAAAAAA".to_vec(), b"CCCCCCCC".to_vec()];
    for extremo in [5e-324, 1e-310, 1e-7, 1e7, 1e308, f64::INFINITY, f64::NAN] {
        let pedido = Opciones {
            pseudoconteo: extremo,
            ..opciones(3, Modelo::Zoops, false)
        };
        assert!(meme(&textos, &pedido).is_err());
    }
    for limite in [PSEUDOCONTEO_MINIMO, PSEUDOCONTEO_MAXIMO] {
        let pedido = Opciones {
            pseudoconteo: limite,
            ..opciones(3, Modelo::Zoops, true)
        };
        let hallado = meme(&textos, &pedido).unwrap();
        assert!(hallado.verosimilitud.is_finite());
        assert!(hallado.gamma.is_finite());
        assert!(
            hallado
                .perfil
                .iter()
                .flatten()
                .all(|p| p.is_finite() && *p > 0.0)
        );
        assert!(hallado.sitios.iter().all(|s| s.probabilidad.is_finite()));
    }
}
