use nucleo::azar::Azar;
use nucleo::replicacion::*;
use nucleo::secuencia::complemento_inverso;
use std::collections::HashMap;

fn adn_al_azar(azar: &mut Azar, longitud: usize) -> Vec<u8> {
    (0..longitud).map(|_| b"ACGT"[azar.entero(4)]).collect()
}

fn hamming_simple(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).filter(|(x, y)| x != y).count()
}

fn todos_los_kmeros(k: usize) -> Vec<Vec<u8>> {
    let mut lista = vec![Vec::new()];
    for _ in 0..k {
        let mut siguiente = Vec::new();
        for prefijo in &lista {
            for &base in b"ACGT" {
                let mut nuevo = prefijo.clone();
                nuevo.push(base);
                siguiente.push(nuevo);
            }
        }
        lista = siguiente;
    }
    lista
}

#[test]
fn cuenta_apariciones_solapadas() {
    assert_eq!(contar_patron(b"CGCGCGAT", b"CGCG"), 2);
    assert_eq!(contar_patron(b"AAAA", b"AA"), 3);
    assert_eq!(contar_patron(b"ACGT", b"TT"), 0);
}

#[test]
fn contar_patron_en_bordes() {
    assert_eq!(contar_patron(b"ACGT", b""), 0);
    assert_eq!(contar_patron(b"", b"A"), 0);
    assert_eq!(contar_patron(b"AC", b"ACG"), 0);
    assert_eq!(contar_patron(b"ACG", b"ACG"), 1);
}

#[test]
fn tabla_de_frecuencias_cuenta_cada_kmero() {
    let tabla = tabla_de_frecuencias(b"ACGTTACGTAACG", 3);
    assert_eq!(tabla[&b"ACG"[..]], 3);
    assert_eq!(tabla[&b"CGT"[..]], 2);
    assert_eq!(tabla[&b"AAC"[..]], 1);
    assert_eq!(tabla.values().sum::<usize>(), 11);
    assert!(tabla_de_frecuencias(b"ACGT", 0).is_empty());
    assert!(tabla_de_frecuencias(b"ACGT", 5).is_empty());
}

#[test]
fn palabras_frecuentes_con_ganador_unico() {
    assert_eq!(
        palabras_frecuentes(b"ACGTTACGTAACG", 3),
        vec![b"ACG".to_vec()]
    );
}

#[test]
fn palabras_frecuentes_con_empate_ordenado() {
    assert_eq!(
        palabras_frecuentes(b"AACCAACC", 2),
        vec![b"AA".to_vec(), b"AC".to_vec(), b"CC".to_vec()]
    );
}

#[test]
fn palabras_frecuentes_en_bordes() {
    assert!(palabras_frecuentes(b"ACGT", 0).is_empty());
    assert!(palabras_frecuentes(b"ACGT", 9).is_empty());
    assert!(palabras_frecuentes(b"", 3).is_empty());
    assert_eq!(palabras_frecuentes(b"ACG", 3), vec![b"ACG".to_vec()]);
}

#[test]
fn tabla_de_frecuencias_ignora_letras_raras() {
    let tabla = tabla_de_frecuencias(b"ACNAC", 2);
    assert_eq!(tabla.len(), 1);
    assert_eq!(tabla[&b"AC"[..]], 2);
}

#[test]
fn conviene_un_arreglo_segun_el_largo() {
    assert!(conviene_un_arreglo(4, 2));
    assert!(!conviene_un_arreglo(3, 2));
    assert!(!conviene_un_arreglo(1000, 0));
    assert!(!conviene_un_arreglo(usize::MAX, 13));
    assert!(conviene_un_arreglo(usize::MAX, 12));
}

fn frecuentes_con_tabla_simple(texto: &[u8], k: usize) -> Vec<Vec<u8>> {
    let mut cuenta: HashMap<&[u8], usize> = HashMap::new();
    for kmero in texto.windows(k) {
        if kmero.iter().all(|b| b"ACGT".contains(b)) {
            *cuenta.entry(kmero).or_insert(0) += 1;
        }
    }
    let maximo = cuenta.values().copied().max().unwrap_or(0);
    let mut lista: Vec<Vec<u8>> = cuenta
        .into_iter()
        .filter(|(_, n)| *n == maximo)
        .map(|(kmero, _)| kmero.to_vec())
        .collect();
    lista.sort();
    lista
}

#[test]
fn palabras_frecuentes_igual_con_arreglo_o_tabla() {
    let mut azar = Azar::nuevo(17);
    for _ in 0..40 {
        let longitud = 1 + azar.entero(400);
        let mut texto = adn_al_azar(&mut azar, longitud);
        if azar.entero(3) == 0 {
            let donde = azar.entero(longitud);
            texto[donde] = b'N';
        }
        let k = 1 + azar.entero(5);
        assert_eq!(
            palabras_frecuentes(&texto, k),
            frecuentes_con_tabla_simple(&texto, k)
        );
    }
}

#[test]
fn posiciones_exactas_desde_cero() {
    assert_eq!(
        posiciones_del_patron(b"TTGATTGATTGA", b"TTGA"),
        vec![0, 4, 8]
    );
    assert_eq!(posiciones_del_patron(b"AAAAA", b"AAA"), vec![0, 1, 2]);
    assert!(posiciones_del_patron(b"ACGT", b"").is_empty());
    assert!(posiciones_del_patron(b"AC", b"ACGT").is_empty());
}

fn grumo(patron: &[u8], ventana: usize) -> Grumo {
    Grumo {
        patron: patron.to_vec(),
        ventana,
    }
}

#[test]
fn grumos_con_dos_apariciones() {
    let texto = b"CATGCATGTTTTTTTCAT";
    assert_eq!(
        encontrar_grumos(texto, 3, 8, 2),
        Ok(vec![grumo(b"ATG", 0), grumo(b"CAT", 0), grumo(b"TTT", 4)])
    );
}

#[test]
fn grumos_con_tres_apariciones() {
    assert_eq!(
        encontrar_grumos(b"CATGCATGTTTTTTTCAT", 3, 8, 3),
        Ok(vec![grumo(b"TTT", 5)])
    );
}

#[test]
fn grumos_con_una_aparicion_marcan_la_primera_ventana() {
    let grumos = encontrar_grumos(b"CATGCATGTTTTTTTCAT", 3, 8, 1).unwrap();
    let ttc = grumos.iter().find(|g| g.patron == b"TTC").unwrap();
    assert_eq!(ttc.ventana, 8);
    assert_eq!(grumos.len(), 9);
}

#[test]
fn grumos_en_bordes() {
    assert!(encontrar_grumos(b"ACGTACGT", 0, 4, 1).is_err());
    assert!(encontrar_grumos(b"ACGTACGTACGTACGT", 13, 16, 1).is_err());
    assert_eq!(encontrar_grumos(b"ACGTACGT", 3, 2, 1), Ok(Vec::new()));
    assert_eq!(encontrar_grumos(b"ACGTACGT", 3, 20, 1), Ok(Vec::new()));
    assert_eq!(encontrar_grumos(b"ACGTACGT", 3, 6, 0), Ok(Vec::new()));
    assert_eq!(encontrar_grumos(b"ACGTACGT", 3, 6, 9), Ok(Vec::new()));
    assert_eq!(
        encontrar_grumos(b"ACGTACGT", 4, 8, 2),
        Ok(vec![grumo(b"ACGT", 0)])
    );
    assert_eq!(
        encontrar_grumos(b"ACGNACGTACG", 3, 11, 2),
        Ok(vec![grumo(b"ACG", 0)])
    );
}

fn grumos_por_fuerza_bruta(texto: &[u8], k: usize, largo: usize, veces: usize) -> Vec<Grumo> {
    let mut primera: HashMap<Vec<u8>, usize> = HashMap::new();
    for inicio in 0..=texto.len() - largo {
        let mut cuenta: HashMap<&[u8], usize> = HashMap::new();
        for kmero in texto[inicio..inicio + largo].windows(k) {
            *cuenta.entry(kmero).or_insert(0) += 1;
        }
        for (kmero, n) in cuenta {
            if n >= veces {
                primera.entry(kmero.to_vec()).or_insert(inicio);
            }
        }
    }
    let mut lista: Vec<Grumo> = primera
        .into_iter()
        .map(|(patron, ventana)| Grumo { patron, ventana })
        .collect();
    lista.sort_by(|a, b| a.patron.cmp(&b.patron));
    lista
}

#[test]
fn grumos_coinciden_con_fuerza_bruta() {
    let mut azar = Azar::nuevo(11);
    for _ in 0..40 {
        let longitud = 30 + azar.entero(120);
        let texto = adn_al_azar(&mut azar, longitud);
        let k = 1 + azar.entero(4);
        let largo = k + azar.entero(25);
        let veces = 1 + azar.entero(4);
        if largo > texto.len() {
            continue;
        }
        assert_eq!(
            encontrar_grumos(&texto, k, largo, veces).unwrap(),
            grumos_por_fuerza_bruta(&texto, k, largo, veces)
        );
    }
}

#[test]
fn sesgo_paso_a_paso() {
    assert_eq!(sesgo(b"GGCAC"), vec![0, 1, 2, 1, 1, 0]);
    assert_eq!(sesgo(b""), vec![0]);
    assert_eq!(sesgo(b"cgN"), vec![0, -1, 0, 0]);
}

#[test]
fn cambio_de_sesgo_por_base() {
    assert_eq!(cambio_de_sesgo(b'G'), 1);
    assert_eq!(cambio_de_sesgo(b'C'), -1);
    assert_eq!(cambio_de_sesgo(b'A'), 0);
    assert_eq!(cambio_de_sesgo(b'N'), 0);
}

#[test]
fn minimo_del_sesgo_unico() {
    assert_eq!(posiciones_de_sesgo_minimo(b"CCGGAC"), vec![2]);
    assert_eq!(posiciones_de_sesgo_minimo(b"CGCAC"), vec![5]);
}

#[test]
fn minimo_del_sesgo_con_empates() {
    assert_eq!(posiciones_de_sesgo_minimo(b"CGC"), vec![1, 3]);
    assert_eq!(posiciones_de_sesgo_minimo(b"ATAT"), vec![0, 1, 2, 3, 4]);
    assert_eq!(posiciones_de_sesgo_minimo(b"CAAG"), vec![1, 2, 3]);
}

#[test]
fn minimo_del_sesgo_en_bordes() {
    assert_eq!(posiciones_de_sesgo_minimo(b""), vec![0]);
    assert_eq!(posiciones_de_sesgo_minimo(b"GGG"), vec![0]);
    assert_eq!(posiciones_de_sesgo_minimo(b"cg"), vec![1]);
}

#[test]
fn minimo_del_sesgo_coincide_con_el_arreglo() {
    let mut azar = Azar::nuevo(5);
    for _ in 0..30 {
        let longitud = azar.entero(300);
        let genoma = adn_al_azar(&mut azar, longitud);
        let valores = sesgo(&genoma);
        let minimo = *valores.iter().min().unwrap();
        let esperado: Vec<usize> = (0..valores.len())
            .filter(|&i| valores[i] == minimo)
            .collect();
        assert_eq!(posiciones_de_sesgo_minimo(&genoma), esperado);
    }
}

#[test]
fn distancia_de_hamming() {
    assert_eq!(distancia_hamming(b"GATTACA", b"GACTATA"), Ok(2));
    assert_eq!(distancia_hamming(b"", b""), Ok(0));
    assert_eq!(distancia_hamming(b"ACGT", b"TGCA"), Ok(4));
    assert!(distancia_hamming(b"ACG", b"AC").is_err());
}

#[test]
fn casi_iguales_respeta_el_limite() {
    assert!(casi_iguales(b"ACGT", b"ACGA", 1));
    assert!(!casi_iguales(b"ACGT", b"TCGA", 1));
    assert!(casi_iguales(b"ACGT", b"TCGA", 2));
    assert!(casi_iguales(b"ACGT", b"ACGT", 0));
}

#[test]
fn posiciones_aproximadas_segun_d() {
    let texto = b"TTACGTTACCTTAGG";
    assert_eq!(posiciones_aproximadas(texto, b"TTAC", 0), vec![0, 5]);
    assert_eq!(posiciones_aproximadas(texto, b"TTAC", 1), vec![0, 5, 10]);
    assert_eq!(posiciones_aproximadas(texto, b"TTAC", 2), vec![0, 5, 6, 10]);
    assert_eq!(posiciones_aproximadas(texto, b"TTAC", 4).len(), 12);
    assert!(posiciones_aproximadas(texto, b"", 1).is_empty());
    assert!(posiciones_aproximadas(b"TTA", b"TTAC", 1).is_empty());
}

#[test]
fn conteo_aproximado_segun_d() {
    let texto = b"TTACGTTACCTTAGG";
    assert_eq!(contar_aproximado(texto, b"TTAC", 0), 2);
    assert_eq!(contar_aproximado(texto, b"TTAC", 1), 3);
    assert_eq!(contar_aproximado(texto, b"TTAC", 2), 4);
    assert_eq!(contar_aproximado(texto, b"", 2), 0);
}

#[test]
fn aproximadas_coinciden_con_hamming() {
    let mut azar = Azar::nuevo(21);
    for _ in 0..30 {
        let longitud = 20 + azar.entero(200);
        let texto = adn_al_azar(&mut azar, longitud);
        let k = 1 + azar.entero(8);
        let patron = adn_al_azar(&mut azar, k);
        let d = azar.entero(k + 1);
        let esperado: Vec<usize> = (0..=texto.len() - k)
            .filter(|&i| hamming_simple(&texto[i..i + k], &patron) <= d)
            .collect();
        assert_eq!(posiciones_aproximadas(&texto, &patron, d), esperado);
        assert_eq!(contar_aproximado(&texto, &patron, d), esperado.len());
    }
}

#[test]
fn valor_de_cada_base() {
    assert_eq!(valor_de_base(b'A'), Some(0));
    assert_eq!(valor_de_base(b'c'), Some(1));
    assert_eq!(valor_de_base(b'G'), Some(2));
    assert_eq!(valor_de_base(b't'), Some(3));
    assert_eq!(valor_de_base(b'N'), None);
    assert_eq!(valor_de_base(b'U'), None);
    assert_eq!(valor_de_base(0xff), None);
}

#[test]
fn numero_de_patron_en_base_cuatro() {
    assert_eq!(numero_de_patron(b"A"), Ok(0));
    assert_eq!(numero_de_patron(b"T"), Ok(3));
    assert_eq!(numero_de_patron(b"GT"), Ok(11));
    assert_eq!(numero_de_patron(b"CAT"), Ok(19));
    assert_eq!(numero_de_patron(b"TTTT"), Ok(255));
    assert_eq!(numero_de_patron(b""), Ok(0));
    assert_eq!(numero_de_patron(&[b'T'; 31]), Ok((1u64 << 62) - 1));
}

#[test]
fn numero_de_patron_rechaza_lo_invalido() {
    assert!(numero_de_patron(b"ACGN").is_err());
    assert!(numero_de_patron(&[b'A'; 32]).is_err());
}

#[test]
fn patron_de_numero_inverso() {
    assert_eq!(patron_de_numero(19, 3), Ok(b"CAT".to_vec()));
    assert_eq!(patron_de_numero(0, 4), Ok(b"AAAA".to_vec()));
    assert_eq!(patron_de_numero(255, 4), Ok(b"TTTT".to_vec()));
    assert_eq!(patron_de_numero(11, 4), Ok(b"AAGT".to_vec()));
    assert_eq!(patron_de_numero(0, 0), Ok(Vec::new()));
    assert!(patron_de_numero(256, 4).is_err());
    assert!(patron_de_numero(0, 32).is_err());
    assert!(patron_de_numero(u64::MAX, 31).is_err());
}

#[test]
fn ida_y_vuelta_entre_patron_y_numero() {
    let mut azar = Azar::nuevo(3);
    for _ in 0..200 {
        let k = 1 + azar.entero(31);
        let patron = adn_al_azar(&mut azar, k);
        let numero = numero_de_patron(&patron).unwrap();
        assert_eq!(patron_de_numero(numero, k).unwrap(), patron);
    }
}

#[test]
fn arreglo_de_frecuencias_por_numero() {
    assert_eq!(arreglo_de_frecuencias(b"ACGTA", 1), Ok(vec![2, 1, 1, 1]));
    let arreglo = arreglo_de_frecuencias(b"ACGTA", 2).unwrap();
    assert_eq!(arreglo.len(), 16);
    let mut esperado = vec![0u32; 16];
    esperado[1] = 1;
    esperado[6] = 1;
    esperado[11] = 1;
    esperado[12] = 1;
    assert_eq!(arreglo, esperado);
}

#[test]
fn arreglo_de_frecuencias_salta_letras_raras() {
    let arreglo = arreglo_de_frecuencias(b"ACNGT", 2).unwrap();
    assert_eq!(arreglo.iter().sum::<u32>(), 2);
    assert_eq!(arreglo[1], 1);
    assert_eq!(arreglo[11], 1);
}

#[test]
fn arreglo_de_frecuencias_limites() {
    assert!(arreglo_de_frecuencias(b"ACGT", 0).is_err());
    assert!(arreglo_de_frecuencias(b"ACGT", 13).is_err());
    assert_eq!(
        arreglo_de_frecuencias(b"AC", 3)
            .unwrap()
            .iter()
            .sum::<u32>(),
        0
    );
}

#[test]
fn tamano_del_vecindario_por_formula() {
    assert_eq!(tamano_del_vecindario(2, 1), 7.0);
    assert_eq!(tamano_del_vecindario(9, 1), 28.0);
    assert_eq!(tamano_del_vecindario(9, 3), 2620.0);
    assert_eq!(tamano_del_vecindario(12, 3), 6571.0);
    assert_eq!(tamano_del_vecindario(3, 0), 1.0);
    assert_eq!(tamano_del_vecindario(2, 2), 16.0);
    assert!(tamano_del_vecindario(4_000_000, 4_000_000).is_infinite());
}

#[test]
fn vecindario_pequeno() {
    let esperado: Vec<Vec<u8>> = ["AA", "AC", "AG", "AT", "CC", "GC", "TC"]
        .iter()
        .map(|p| p.as_bytes().to_vec())
        .collect();
    assert_eq!(vecindario(b"AC", 1), Ok(esperado));
    assert_eq!(vecindario(b"ACG", 0), Ok(vec![b"ACG".to_vec()]));
    assert_eq!(vecindario(b"AC", 2).unwrap(), todos_los_kmeros(2));
}

#[test]
fn vecindario_cumple_la_formula_y_la_distancia() {
    for (patron, d) in [(&b"TTATCCACA"[..], 1), (b"ACGTAC", 2), (b"GATTACA", 3)] {
        let vecinos = vecindario(patron, d).unwrap();
        assert_eq!(vecinos.len() as f64, tamano_del_vecindario(patron.len(), d));
        assert!(vecinos.windows(2).all(|par| par[0] < par[1]));
        assert!(vecinos.iter().all(|v| hamming_simple(v, patron) <= d));
    }
}

#[test]
fn vecindario_rechaza_lo_imposible() {
    assert!(vecindario(b"ACG", 4).is_err());
    assert!(vecindario(b"ACN", 1).is_err());
    assert!(vecindario(&[b'A'; 32], 1).is_err());
    assert!(vecindario(&[b'A'; 31], 5).is_err());
    assert_eq!(vecindario(&[b'A'; 31], 1).unwrap().len(), 94);
    assert_eq!(vecindario(b"acg", 0), Ok(vec![b"ACG".to_vec()]));
}

#[test]
fn frecuentes_con_errores_a_mano() {
    let esperado: Vec<Vec<u8>> = ["AG", "AT", "CC", "CT", "GC", "GG"]
        .iter()
        .map(|p| p.as_bytes().to_vec())
        .collect();
    assert_eq!(palabras_frecuentes_con_errores(b"ACGT", 2, 1), Ok(esperado));
}

#[test]
fn frecuentes_con_complementos_a_mano() {
    assert_eq!(
        palabras_frecuentes_con_errores(b"AAAA", 2, 0),
        Ok(vec![b"AA".to_vec()])
    );
    assert_eq!(
        palabras_frecuentes_con_errores_y_complementos(b"AAAA", 2, 0),
        Ok(vec![b"AA".to_vec(), b"TT".to_vec()])
    );
    assert_eq!(
        palabras_frecuentes_con_errores_y_complementos(b"AAAA", 2, 1),
        Ok(vec![b"AT".to_vec(), b"TA".to_vec()])
    );
}

#[test]
fn cuentas_con_errores_por_numero() {
    let cuentas = cuentas_con_errores(b"AAAA", 2, 1).unwrap();
    assert_eq!(cuentas[0], 3);
    assert_eq!(cuentas[1], 3);
    assert_eq!(cuentas[5], 0);
    let ambas = cuentas_en_ambas_hebras(b"AAAA", 2, 1).unwrap();
    assert_eq!(ambas[3], 6);
    assert_eq!(ambas[12], 6);
    assert_eq!(ambas[15], 3);
}

fn frecuentes_por_fuerza_bruta(texto: &[u8], k: usize, d: usize, ambas: bool) -> Vec<Vec<u8>> {
    let mut mejores = Vec::new();
    let mut maximo = 0;
    for patron in todos_los_kmeros(k) {
        let mut veces = 0;
        for ventana in texto.windows(k) {
            if hamming_simple(ventana, &patron) <= d {
                veces += 1;
            }
            if ambas && hamming_simple(ventana, &complemento_inverso(&patron)) <= d {
                veces += 1;
            }
        }
        if veces > maximo {
            maximo = veces;
            mejores.clear();
        }
        if veces == maximo && veces > 0 {
            mejores.push(patron);
        }
    }
    mejores
}

#[test]
fn frecuentes_con_errores_coinciden_con_fuerza_bruta() {
    let mut azar = Azar::nuevo(8);
    for _ in 0..25 {
        let longitud = 10 + azar.entero(60);
        let texto = adn_al_azar(&mut azar, longitud);
        let k = 1 + azar.entero(5);
        let d = azar.entero(k.min(3) + 1);
        assert_eq!(
            palabras_frecuentes_con_errores(&texto, k, d).unwrap(),
            frecuentes_por_fuerza_bruta(&texto, k, d, false)
        );
        assert_eq!(
            palabras_frecuentes_con_errores_y_complementos(&texto, k, d).unwrap(),
            frecuentes_por_fuerza_bruta(&texto, k, d, true)
        );
    }
}

#[test]
fn frecuentes_con_errores_limites() {
    assert!(palabras_frecuentes_con_errores(b"ACGT", 0, 0).is_err());
    assert!(palabras_frecuentes_con_errores(b"ACGT", 13, 1).is_err());
    assert!(palabras_frecuentes_con_errores(b"ACGT", 2, 3).is_err());
    assert!(palabras_frecuentes_con_errores(b"ACGTACGTAC", 5, 4).is_err());
    assert_eq!(palabras_frecuentes_con_errores(b"AC", 3, 1), Ok(Vec::new()));
    assert_eq!(
        palabras_frecuentes_con_errores_y_complementos(b"NNNN", 2, 1),
        Ok(Vec::new())
    );
}

#[test]
fn los_mas_frecuentes_salen_en_orden() {
    let mut cuentas = vec![0u32; 16];
    cuentas[9] = 4;
    cuentas[2] = 4;
    cuentas[3] = 1;
    assert_eq!(
        los_mas_frecuentes(&cuentas, 2),
        vec![b"AG".to_vec(), b"GC".to_vec()]
    );
    assert!(los_mas_frecuentes(&[0u32; 16], 2).is_empty());
}

#[test]
fn ori_en_un_genoma_inventado() {
    let genoma = b"CCCCCCCCCCATTAATTAATTAATGGGGGGGGGG";
    let ori = buscar_ori(genoma, 14, 4, 0).unwrap();
    assert_eq!(ori.minimos, (10..=24).collect::<Vec<usize>>());
    assert_eq!(ori.inicio, 10);
    assert_eq!(ori.fin, 24);
    assert_eq!(ori.veces, 6);
    assert_eq!(
        ori.candidatos,
        vec![b"ATTA".to_vec(), b"TAAT".to_vec(), b"TTAA".to_vec()]
    );
}

#[test]
fn ori_con_minimo_cerca_del_final() {
    let genoma = b"GGGGATATCCCC";
    let ori = buscar_ori(genoma, 6, 2, 0).unwrap();
    assert_eq!(ori.minimos, vec![0, 12]);
    assert_eq!(ori.inicio, 0);
    assert_eq!(ori.fin, 6);
    let genoma = b"GGGGGGATATCC";
    let ori = buscar_ori(genoma, 6, 2, 0).unwrap();
    assert_eq!(ori.minimos, vec![0]);
    let genoma = b"GGAAAAAAACCCCC";
    let ori = buscar_ori(genoma, 6, 2, 0).unwrap();
    assert_eq!(ori.minimos, vec![14]);
    assert_eq!(ori.inicio, 8);
    assert_eq!(ori.fin, 14);
}

#[test]
fn ori_rechaza_ventanas_imposibles() {
    assert!(buscar_ori(b"ACGTACGT", 20, 4, 0).is_err());
    assert!(buscar_ori(b"ACGTACGT", 3, 4, 0).is_err());
    assert!(buscar_ori(b"ACGTACGT", 8, 0, 0).is_err());
    assert!(buscar_ori(b"ACGTACGT", 8, 4, 5).is_err());
}

#[test]
fn ecoli_coincide_con_una_verificacion_independiente() {
    let genoma = nucleo::ecoli::genoma();
    assert_eq!(
        posiciones_de_sesgo_minimo(genoma),
        vec![3925597, 3925598, 3925599, 3925600]
    );
    assert_eq!(posiciones_del_patron(genoma, b"TTATCCACA").len(), 22);
    assert_eq!(posiciones_del_patron(genoma, b"TGTGGATAA").len(), 10);
    assert_eq!(encontrar_grumos(genoma, 9, 500, 3).unwrap().len(), 1904);
    let origen = buscar_ori(genoma, 500, 9, 1).unwrap();
    assert_eq!(origen.inicio, 3925597);
    assert_eq!(origen.veces, 4);
    assert!(origen.candidatos.contains(&b"TTATCCACA".to_vec()));
}
