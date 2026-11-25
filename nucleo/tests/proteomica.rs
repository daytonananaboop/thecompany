use nucleo::azar::Azar;
use nucleo::proteomica::{
    Busqueda, Digestion, Espectro, Indice, MASA_DEL_AGUA, MASA_DEL_PROTON, PeptidoIndexado, Psm,
    Puntaje, aceptados, asignar_valores_q, buscar_espectros, corta_despues, digerir, fronteras,
    hiperpuntaje, leer_mgf, llenar_iones, masa_de_residuo, masa_del_peptido, masa_neutra,
    mejor_psm, mz_con_carga, normalizar, pico_cercano, puntuar, senuelo, ventana,
};

fn cerca(a: f64, b: f64, margen: f64) -> bool {
    (a - b).abs() <= margen
}

fn sin_cisteina() -> Digestion {
    Digestion {
        carbamidometilo: false,
        ..Digestion::default()
    }
}

#[test]
fn masa_de_peptide_conocida() {
    let masa = masa_del_peptido(b"PEPTIDE", false).unwrap();
    assert!(cerca(masa, 799.359_96, 1e-4), "{masa}");
    assert!(cerca(mz_con_carga(masa, 2), 400.687_26, 1e-4));
    assert!(cerca(masa_neutra(mz_con_carga(masa, 3), 3), masa, 1e-9));
    assert_eq!(masa_del_peptido(b"PEPXIDE", false), None);
    assert_eq!(masa_del_peptido(b"", false), Some(MASA_DEL_AGUA));
}

#[test]
fn la_cisteina_lleva_carbamidometilo_si_se_pide() {
    let libre = masa_de_residuo(b'C', false).unwrap();
    let fija = masa_de_residuo(b'C', true).unwrap();
    assert!(cerca(fija - libre, 57.021_464, 1e-9));
    assert_eq!(masa_de_residuo(b'I', true), masa_de_residuo(b'L', true));
    assert_eq!(masa_de_residuo(b'U', true), None);
}

#[test]
fn tripsina_corta_despues_de_k_y_r_salvo_antes_de_p() {
    let proteina = b"AKPRKGG";
    assert!(!corta_despues(proteina, 1));
    assert!(corta_despues(proteina, 3));
    assert!(corta_despues(proteina, 4));
    assert!(!corta_despues(proteina, 0));
    assert_eq!(fronteras(proteina), vec![0, 4, 5, 7]);
    assert_eq!(fronteras(b"AAK"), vec![0, 3]);
    assert_eq!(fronteras(b""), vec![0]);
}

#[test]
fn digestion_a_mano() {
    let proteina = b"MKAAAAAARGGGGGGKPLLLLLLK";
    let reglas = Digestion {
        cortes_perdidos: 1,
        largo_minimo: 6,
        largo_maximo: 40,
        carbamidometilo: true,
    };
    let peptidos: Vec<(String, usize)> = digerir(proteina, &reglas)
        .iter()
        .map(|p| {
            (
                String::from_utf8(proteina[p.inicio..p.fin].to_vec()).unwrap(),
                p.cortes_perdidos,
            )
        })
        .collect();
    assert_eq!(
        peptidos,
        vec![
            ("MKAAAAAAR".to_string(), 1),
            ("AAAAAAR".to_string(), 0),
            ("AAAAAARGGGGGGKPLLLLLLK".to_string(), 1),
            ("GGGGGGKPLLLLLLK".to_string(), 0),
        ]
    );
    let sin_perdidos = Digestion {
        cortes_perdidos: 0,
        ..reglas
    };
    assert_eq!(digerir(proteina, &sin_perdidos).len(), 2);
    let cortos = Digestion {
        largo_maximo: 10,
        ..reglas
    };
    assert_eq!(digerir(proteina, &cortos).len(), 2);
}

#[test]
fn los_peptidos_con_letras_raras_se_descartan() {
    let peptidos = digerir(b"AAAAUAAKGGGGGGR", &Digestion::default());
    assert_eq!(peptidos.len(), 1);
    assert_eq!(peptidos[0].inicio, 8);
}

#[test]
fn iones_b_y_de_peptide() {
    let mut b = Vec::new();
    let mut y = Vec::new();
    llenar_iones(b"PEPTIDE", false, &mut b, &mut y);
    assert_eq!(b.len(), 6);
    assert_eq!(y.len(), 6);
    assert!(cerca(b[0], 98.060_04, 1e-4));
    assert!(cerca(b[1], 227.102_63, 1e-4));
    assert!(cerca(y[0], 148.060_43, 1e-4));
    assert!(cerca(y[1], 263.087_37, 1e-4));
    let masa = masa_del_peptido(b"PEPTIDE", false).unwrap();
    for i in 0..6 {
        assert!(cerca(b[i] + y[5 - i], masa + 2.0 * MASA_DEL_PROTON, 1e-9));
    }
    llenar_iones(b"K", false, &mut b, &mut y);
    assert!(b.is_empty() && y.is_empty());
    llenar_iones(b"", false, &mut b, &mut y);
    assert!(b.is_empty() && y.is_empty());
}

#[test]
fn el_pico_mas_cercano_dentro_de_la_tolerancia() {
    let picos = [(100.0, 1.0), (100.03, 5.0), (200.0, 2.0)];
    assert_eq!(pico_cercano(&picos, 100.01, 0.02), Some(0));
    assert_eq!(pico_cercano(&picos, 100.025, 0.02), Some(1));
    assert_eq!(pico_cercano(&picos, 100.015, 0.02), Some(0));
    assert_eq!(pico_cercano(&picos, 150.0, 0.02), None);
    assert_eq!(pico_cercano(&picos, 50.0, 0.02), None);
    assert_eq!(pico_cercano(&picos, 200.015, 0.02), Some(2));
    assert_eq!(pico_cercano(&[], 10.0, 1.0), None);
}

#[test]
fn hiperpuntaje_a_mano() {
    assert_eq!(hiperpuntaje(50.0, 0, 0), 0.0);
    let esperado = 51.0_f64.ln() + 6.0_f64.ln() + 2.0_f64.ln();
    assert!(cerca(hiperpuntaje(50.0, 3, 2), esperado, 1e-12));
    let picos = [(98.06, 100.0), (148.06, 50.0), (300.0, 10.0)];
    let puntaje = puntuar(&[98.06, 227.1], &[148.06], &picos, 0.02);
    assert_eq!(puntaje.iones_b, 1);
    assert_eq!(puntaje.iones_y, 1);
    assert!(cerca(puntaje.intensidad, 150.0, 1e-12));
    assert!(cerca(puntaje.hiperpuntaje, 151.0_f64.ln(), 1e-12));
}

#[test]
fn normalizar_lleva_el_maximo_a_cien() {
    let picos = normalizar(&[(1.0, 5.0), (2.0, 20.0)]);
    assert_eq!(picos, vec![(1.0, 25.0), (2.0, 100.0)]);
    assert!(normalizar(&[]).is_empty());
}

#[test]
fn ventana_en_ppm() {
    let (minima, maxima) = ventana(1000.0, 10.0);
    assert!(cerca(minima, 999.99, 1e-9));
    assert!(cerca(maxima, 1000.01, 1e-9));
}

#[test]
fn leer_mgf_completo() {
    let texto = "MASS=Monoisotopic\r\nBEGIN IONS\r\nTITLE=uno\r\nPEPMASS=500.25 1000\r\nCHARGE=2+ and 3+\r\n300.1 10\r\n200.5 40\r\nEND IONS\r\n\r\n# nota\r\nBEGIN IONS\nPEPMASS=400.2\n150.0\nEND IONS\n";
    let espectros = leer_mgf(texto).unwrap();
    assert_eq!(espectros.len(), 2);
    assert_eq!(espectros[0].titulo, "uno");
    assert_eq!(espectros[0].mz_precursor, 500.25);
    assert_eq!(espectros[0].cargas, vec![2, 3]);
    assert_eq!(espectros[0].picos, vec![(200.5, 40.0), (300.1, 10.0)]);
    assert_eq!(espectros[1].titulo, "espectro 2");
    assert_eq!(espectros[1].cargas, vec![1, 2, 3]);
    assert_eq!(espectros[1].picos, vec![(150.0, 1.0)]);
}

#[test]
fn leer_mgf_con_errores_claros() {
    assert!(leer_mgf("").unwrap_err().contains("ningún espectro"));
    assert!(leer_mgf("hola").unwrap_err().contains("ningún espectro"));
    assert!(
        leer_mgf("BEGIN IONS\nPEPMASS=500\n100 1\n")
            .unwrap_err()
            .contains("END IONS")
    );
    assert!(
        leer_mgf("BEGIN IONS\n100 1\nEND IONS")
            .unwrap_err()
            .contains("PEPMASS")
    );
    assert!(
        leer_mgf("BEGIN IONS\nPEPMASS=abc\nEND IONS")
            .unwrap_err()
            .contains("PEPMASS")
    );
    assert!(
        leer_mgf("BEGIN IONS\nPEPMASS=500\nCHARGE=2-\nEND IONS")
            .unwrap_err()
            .contains("positivos")
    );
    assert!(
        leer_mgf("BEGIN IONS\nPEPMASS=500\nCHARGE=40+\nEND IONS")
            .unwrap_err()
            .contains("carga")
    );
    assert!(
        leer_mgf("BEGIN IONS\nPEPMASS=500\n100 x\nEND IONS")
            .unwrap_err()
            .contains("Línea 3")
    );
    assert!(
        leer_mgf("BEGIN IONS\nPEPMASS=500\nBEGIN IONS\nEND IONS")
            .unwrap_err()
            .contains("sin cerrar")
    );
    let ceros = leer_mgf("BEGIN IONS\nPEPMASS=500\n100 0\n-3 5\nEND IONS").unwrap();
    assert!(ceros[0].picos.is_empty());
}

fn proteinas_de_prueba() -> Vec<(String, Vec<u8>)> {
    vec![
        (
            "uno".to_string(),
            b"MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRR".to_vec(),
        ),
        (
            "dos".to_string(),
            b"MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT".to_vec(),
        ),
        (
            "tres".to_string(),
            b"MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVK".to_vec(),
        ),
    ]
}

fn espectro_perfecto(peptido: &[u8], carga: u8, carbamidometilo: bool) -> Espectro {
    let mut b = Vec::new();
    let mut y = Vec::new();
    llenar_iones(peptido, carbamidometilo, &mut b, &mut y);
    let mut picos: Vec<(f64, f64)> = b.iter().chain(&y).map(|&mz| (mz, 50.0)).collect();
    picos.sort_by(|a, b| a.0.total_cmp(&b.0));
    Espectro {
        titulo: "perfecto".to_string(),
        mz_precursor: mz_con_carga(masa_del_peptido(peptido, carbamidometilo).unwrap(), carga),
        cargas: vec![carga],
        picos,
    }
}

#[test]
fn el_indice_tiene_objetivos_y_senuelos_ordenados_por_masa() {
    let indice = Indice::construir(&proteinas_de_prueba(), Digestion::default()).unwrap();
    assert_eq!(indice.proteinas.len(), 6);
    assert!(indice.proteinas[3].senuelo);
    assert_eq!(
        indice.proteinas[3].secuencia,
        senuelo(&indice.proteinas[0].secuencia)
    );
    assert!(indice.proteinas[3].nombre.ends_with("uno"));
    assert!(
        indice
            .peptidos
            .windows(2)
            .all(|par| par[0].masa <= par[1].masa)
    );
    assert!(indice.peptidos.iter().any(|p| indice.es_senuelo(p)));
    for peptido in &indice.peptidos {
        let masa = masa_del_peptido(indice.secuencia(peptido), true).unwrap();
        assert!(cerca(masa, peptido.masa, 1e-9));
    }
    assert!(Indice::construir(&[], Digestion::default()).is_err());
    let malas = Digestion {
        cortes_perdidos: 9,
        ..Digestion::default()
    };
    assert!(Indice::construir(&proteinas_de_prueba(), malas).is_err());
}

#[test]
fn un_espectro_perfecto_encuentra_su_peptido() {
    let indice = Indice::construir(&proteinas_de_prueba(), Digestion::default()).unwrap();
    let espectro = espectro_perfecto(b"GFFYTPK", 2, true);
    let psm = mejor_psm(0, &espectro, &indice, &Busqueda::default()).unwrap();
    assert_eq!(indice.secuencia(&psm.peptido), b"GFFYTPK");
    assert!(!psm.senuelo);
    assert_eq!(psm.puntaje.iones_b, 6);
    assert_eq!(psm.puntaje.iones_y, 6);
    assert_eq!(psm.carga, 2);
    assert!(psm.error_ppm.abs() < 1e-6);
    assert!(psm.candidatos >= 1);
}

#[test]
fn la_cisteina_cambia_la_masa_buscada() {
    let reglas = sin_cisteina();
    let indice = Indice::construir(&proteinas_de_prueba(), reglas).unwrap();
    let espectro = espectro_perfecto(b"LLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGER", 3, false);
    let psm = mejor_psm(0, &espectro, &indice, &Busqueda::default()).unwrap();
    assert_eq!(
        indice.secuencia(&psm.peptido),
        b"LLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGER"
    );
    let con_cisteina = Indice::construir(&proteinas_de_prueba(), Digestion::default()).unwrap();
    let otro = mejor_psm(0, &espectro, &con_cisteina, &Busqueda::default());
    assert!(
        otro.is_none_or(|psm| con_cisteina.secuencia(&psm.peptido)
            != b"LLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGER")
    );
}

fn psm_falso(espectro: usize, puntaje: f64, senuelo: bool) -> Psm {
    Psm {
        espectro,
        peptido: PeptidoIndexado {
            masa: 1000.0,
            proteina: 0,
            inicio: 0,
            fin: 1,
            cortes_perdidos: 0,
        },
        carga: 2,
        error_ppm: 0.0,
        puntaje: Puntaje {
            hiperpuntaje: puntaje,
            ..Puntaje::default()
        },
        senuelo,
        candidatos: 1,
        q: 1.0,
    }
}

#[test]
fn valores_q_a_mano() {
    let mut psms = vec![
        psm_falso(0, 10.0, false),
        psm_falso(1, 9.0, false),
        psm_falso(2, 8.0, true),
        psm_falso(3, 7.0, false),
        psm_falso(4, 6.0, false),
        psm_falso(5, 5.0, false),
        psm_falso(6, 4.0, true),
        psm_falso(7, 3.0, true),
    ];
    psms.reverse();
    asignar_valores_q(&mut psms);
    let q: Vec<f64> = psms.iter().map(|psm| psm.q).collect();
    let esperado = [0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.4, 0.6];
    for (obtenido, esperado) in q.iter().zip(esperado) {
        assert!(cerca(*obtenido, esperado, 1e-12), "{q:?}");
    }
    assert_eq!(aceptados(&psms, 0.01).len(), 2);
    assert_eq!(aceptados(&psms, 0.2).len(), 5);
}

#[test]
fn empates_de_puntaje_se_cuentan_juntos() {
    let mut psms = vec![
        psm_falso(0, 5.0, false),
        psm_falso(1, 5.0, true),
        psm_falso(2, 9.0, false),
    ];
    asignar_valores_q(&mut psms);
    assert_eq!(psms[0].q, 0.0);
    assert!(cerca(psms[1].q, 0.5, 1e-12));
    assert!(cerca(psms[2].q, 0.5, 1e-12));
    let mut solo_senuelos = vec![psm_falso(0, 5.0, true)];
    asignar_valores_q(&mut solo_senuelos);
    assert_eq!(solo_senuelos[0].q, 1.0);
}

#[test]
fn busqueda_rechaza_tolerancias_absurdas() {
    let indice = Indice::construir(&proteinas_de_prueba(), Digestion::default()).unwrap();
    let espectro = espectro_perfecto(b"GFFYTPK", 2, true);
    for (ppm, da) in [
        (0.0, 0.02),
        (-1.0, 0.02),
        (5000.0, 0.02),
        (10.0, 0.0),
        (10.0, 3.0),
        (f64::NAN, 0.02),
    ] {
        let busqueda = Busqueda {
            tolerancia_ppm: ppm,
            tolerancia_fragmento: da,
        };
        assert!(buscar_espectros(std::slice::from_ref(&espectro), &indice, &busqueda).is_err());
    }
}

fn azar_de_proteina(azar: &mut Azar, largo: usize) -> Vec<u8> {
    const LETRAS: &[u8] = b"ACDEFGHIKLMNPQRSTVWYKRKRP";
    (0..largo)
        .map(|_| {
            if azar.entero(300) == 0 {
                b'U'
            } else {
                LETRAS[azar.entero(LETRAS.len())]
            }
        })
        .collect()
}

fn digerir_a_la_fuerza(proteina: &[u8], reglas: &Digestion) -> Vec<(usize, usize, usize)> {
    let mut sitios = vec![0];
    for posicion in 0..proteina.len().saturating_sub(1) {
        if matches!(proteina[posicion], b'K' | b'R') && proteina[posicion + 1] != b'P' {
            sitios.push(posicion + 1);
        }
    }
    sitios.push(proteina.len());
    let mut resultado = Vec::new();
    for &inicio in &sitios {
        for &fin in &sitios {
            if fin <= inicio || inicio == proteina.len() {
                continue;
            }
            let perdidos = sitios.iter().filter(|&&s| s > inicio && s < fin).count();
            let largo = fin - inicio;
            let valido = proteina[inicio..fin]
                .iter()
                .all(|&letra| masa_de_residuo(letra, true).is_some());
            if perdidos <= reglas.cortes_perdidos
                && largo >= reglas.largo_minimo
                && largo <= reglas.largo_maximo
                && valido
            {
                resultado.push((inicio, fin, perdidos));
            }
        }
    }
    resultado.sort();
    resultado
}

#[test]
fn digerir_coincide_con_la_fuerza_bruta() {
    let mut azar = Azar::nuevo(41);
    for _ in 0..300 {
        let largo = azar.entero(200);
        let proteina = azar_de_proteina(&mut azar, largo);
        let minimo = 1 + azar.entero(8);
        let reglas = Digestion {
            cortes_perdidos: azar.entero(4),
            largo_minimo: minimo,
            largo_maximo: minimo + azar.entero(40),
            carbamidometilo: true,
        };
        let mut obtenidos: Vec<(usize, usize, usize)> = digerir(&proteina, &reglas)
            .iter()
            .map(|p| (p.inicio, p.fin, p.cortes_perdidos))
            .collect();
        obtenidos.sort();
        assert_eq!(obtenidos, digerir_a_la_fuerza(&proteina, &reglas));
    }
}

fn mejor_a_la_fuerza(
    espectro: &Espectro,
    indice: &Indice,
    busqueda: &Busqueda,
) -> Option<(f64, PeptidoIndexado)> {
    let picos = normalizar(&espectro.picos);
    let mut mejor: Option<(f64, PeptidoIndexado)> = None;
    for &carga in &espectro.cargas {
        let masa = masa_neutra(espectro.mz_precursor, carga);
        for peptido in &indice.peptidos {
            if (masa - peptido.masa).abs() > masa * busqueda.tolerancia_ppm / 1e6 {
                continue;
            }
            let mut b = Vec::new();
            let mut y = Vec::new();
            llenar_iones(indice.secuencia(peptido), true, &mut b, &mut y);
            let mut nb = 0;
            let mut ny = 0;
            let mut intensidad = 0.0;
            for (iones, contador) in [(&b, &mut nb), (&y, &mut ny)] {
                for &ion in iones.iter() {
                    let mut elegido: Option<(f64, f64)> = None;
                    for &(mz, altura) in &picos {
                        let distancia = (mz - ion).abs();
                        if distancia <= busqueda.tolerancia_fragmento
                            && elegido.is_none_or(|(otra, _)| distancia < otra)
                        {
                            elegido = Some((distancia, altura));
                        }
                    }
                    if let Some((_, altura)) = elegido {
                        *contador += 1;
                        intensidad += altura;
                    }
                }
            }
            let puntaje = hiperpuntaje(intensidad, nb, ny);
            if puntaje > mejor.map_or(0.0, |(record, _)| record) {
                mejor = Some((puntaje, *peptido));
            }
        }
    }
    mejor
}

#[test]
fn la_busqueda_indexada_coincide_con_puntuar_todo() {
    let mut azar = Azar::nuevo(7);
    for _ in 0..40 {
        let objetivos: Vec<(String, Vec<u8>)> = (0..2 + azar.entero(6))
            .map(|numero| {
                (format!("p{numero}"), {
                    let largo = 50 + azar.entero(300);
                    azar_de_proteina(&mut azar, largo)
                })
            })
            .collect();
        let indice = Indice::construir(&objetivos, Digestion::default()).unwrap();
        if indice.peptidos.is_empty() {
            continue;
        }
        let busqueda = Busqueda {
            tolerancia_ppm: 5.0 + azar.decimal() * 200.0,
            tolerancia_fragmento: 0.01 + azar.decimal() * 0.3,
        };
        let mut espectros = Vec::new();
        for _ in 0..10 {
            let elegido = indice.peptidos[azar.entero(indice.peptidos.len())];
            let mut espectro =
                espectro_perfecto(indice.secuencia(&elegido), 1 + azar.entero(3) as u8, true);
            espectro.picos.retain(|_| azar.entero(3) > 0);
            for _ in 0..azar.entero(30) {
                espectro
                    .picos
                    .push((100.0 + azar.decimal() * 2000.0, 1.0 + azar.decimal() * 80.0));
            }
            espectro.picos.sort_by(|a, b| a.0.total_cmp(&b.0));
            espectro.mz_precursor += (azar.decimal() - 0.5) * 0.01;
            if azar.entero(4) == 0 {
                espectro.cargas = vec![1, 2, 3];
            }
            espectros.push(espectro);
        }
        let psms = buscar_espectros(&espectros, &indice, &busqueda).unwrap();
        for (numero, espectro) in espectros.iter().enumerate() {
            let esperado = mejor_a_la_fuerza(espectro, &indice, &busqueda);
            let obtenido = psms.iter().find(|psm| psm.espectro == numero);
            match (esperado, obtenido) {
                (None, None) => {}
                (Some((puntaje, peptido)), Some(psm)) => {
                    assert!(cerca(puntaje, psm.puntaje.hiperpuntaje, 1e-9));
                    assert_eq!(indice.secuencia(&peptido), indice.secuencia(&psm.peptido));
                }
                otro => panic!("no coincide: {otro:?}"),
            }
        }
        assert!(
            psms.windows(2)
                .all(|par| par[0].puntaje.hiperpuntaje >= par[1].puntaje.hiperpuntaje)
        );
        assert!(psms.windows(2).all(|par| par[0].q <= par[1].q));
    }
}

#[test]
fn la_ventana_del_indice_coincide_con_recorrer_todo() {
    let mut azar = Azar::nuevo(3);
    let objetivos: Vec<(String, Vec<u8>)> = (0..20)
        .map(|numero| (format!("p{numero}"), azar_de_proteina(&mut azar, 400)))
        .collect();
    let indice = Indice::construir(&objetivos, Digestion::default()).unwrap();
    for _ in 0..500 {
        let centro = 500.0 + azar.decimal() * 3000.0;
        let (minima, maxima) = ventana(centro, azar.decimal() * 500.0);
        let rapido = indice.en_ventana(minima, maxima).len();
        let lento = indice
            .peptidos
            .iter()
            .filter(|p| p.masa >= minima && p.masa <= maxima)
            .count();
        assert_eq!(rapido, lento);
    }
}
