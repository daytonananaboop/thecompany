use nucleo::azar::Azar;
use nucleo::ensamblaje::*;
use std::collections::{HashMap, HashSet};

fn piezas(texto: &str) -> Vec<&[u8]> {
    texto
        .split_whitespace()
        .map(|pieza| pieza.as_bytes())
        .collect()
}

fn etiquetas(grafo: &Grafo, camino: &[usize]) -> Vec<String> {
    camino
        .iter()
        .map(|&nodo| String::from_utf8(grafo.nodos[nodo].clone()).unwrap())
        .collect()
}

fn usa_cada_arista_una_vez(grafo: &Grafo, camino: &[usize]) -> bool {
    let mut restantes: HashMap<(usize, usize), usize> = HashMap::new();
    for (desde, destinos) in grafo.salidas.iter().enumerate() {
        for &hasta in destinos {
            *restantes.entry((desde, hasta)).or_insert(0) += 1;
        }
    }
    for par in camino.windows(2) {
        match restantes.get_mut(&(par[0], par[1])) {
            Some(quedan) if *quedan > 0 => *quedan -= 1,
            _ => return false,
        }
    }
    restantes.values().all(|&quedan| quedan == 0)
}

fn composicion_ordenada(texto: &[u8], k: usize) -> Vec<Vec<u8>> {
    let mut kmeros: Vec<Vec<u8>> = texto.windows(k).map(|kmero| kmero.to_vec()).collect();
    kmeros.sort();
    kmeros
}

fn barajar<T>(lista: &mut [T], semilla: u64) {
    let mut azar = Azar::nuevo(semilla);
    for posicion in (1..lista.len()).rev() {
        let otra = azar.entero(posicion + 1);
        lista.swap(posicion, otra);
    }
}

fn genoma_al_azar(longitud: usize, semilla: u64) -> Vec<u8> {
    let mut azar = Azar::nuevo(semilla);
    (0..longitud).map(|_| b"ACGT"[azar.entero(4)]).collect()
}

fn pares_de(genoma: &[u8], k: usize, d: usize) -> Vec<ParDeLecturas> {
    let mut pares = Vec::new();
    let largo = 2 * k + d;
    for inicio in 0..=genoma.len() - largo {
        pares.push(ParDeLecturas {
            primera: genoma[inicio..inicio + k].to_vec(),
            segunda: genoma[inicio + k + d..inicio + largo].to_vec(),
        });
    }
    pares
}

fn pares_ordenados(pares: &[ParDeLecturas]) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut lista: Vec<(Vec<u8>, Vec<u8>)> = pares
        .iter()
        .map(|par| (par.primera.clone(), par.segunda.clone()))
        .collect();
    lista.sort();
    lista
}

const GENOMA_CON_HUECOS: &[u8] = b"ACGTTGCAAGCT";

fn pares_del_genoma_con_huecos() -> Vec<ParDeLecturas> {
    let lineas = ["ACG|GCA", "CGT|CAA", "GTT|AAG", "TTG|AGC", "TGC|GCT"];
    leer_pares(&lineas.join("\n")).unwrap()
}

#[test]
fn composicion_en_orden_de_aparicion() {
    let kmeros = composicion(b"GATTACAGATT", 4);
    assert_eq!(kmeros, piezas("GATT ATTA TTAC TACA ACAG CAGA AGAT GATT"));
    assert_eq!(composicion(b"GATTACA", 7), piezas("GATTACA"));
}

#[test]
fn composicion_con_k_fuera_de_rango() {
    assert!(composicion(b"GATTACA", 8).is_empty());
    assert!(composicion(b"GATTACA", 0).is_empty());
    assert!(composicion(b"", 3).is_empty());
    assert_eq!(composicion(b"GAT", 1), piezas("G A T"));
}

#[test]
fn deletrea_un_camino_de_kmeros() {
    assert_eq!(
        deletrear_camino(&piezas("GATT ATTA TTAC TACA")).unwrap(),
        b"GATTACA"
    );
    assert_eq!(deletrear_camino(&piezas("CCA")).unwrap(), b"CCA");
    assert_eq!(deletrear_camino(&piezas("A C G")).unwrap(), b"ACG");
}

#[test]
fn deletrear_rechaza_caminos_rotos() {
    assert!(deletrear_camino(&[]).is_err());
    assert!(deletrear_camino(&piezas("GATT ATT")).is_err());
    assert!(deletrear_camino(&[b"".as_slice()]).is_err());
    let error = deletrear_camino(&piezas("GATT ATTA TACA")).unwrap_err();
    assert!(error.contains('3'), "{error}");
    assert!(error.contains("TACA"), "{error}");
}

#[test]
fn grafo_de_solapamiento_une_sufijo_con_prefijo() {
    let grafo = grafo_de_solapamiento(&piezas("ATG TGC GCA TGA GCT")).unwrap();
    assert_eq!(
        grafo.salidas,
        vec![vec![1, 3], vec![2, 4], vec![], vec![], vec![]]
    );
    assert_eq!(grafo.nodos[3], b"TGA");
    assert_eq!(cantidad_de_aristas(&grafo), 4);
}

#[test]
fn solapamiento_con_kmeros_repetidos_no_hace_bucles() {
    let grafo = grafo_de_solapamiento(&piezas("AAA AAA")).unwrap();
    assert_eq!(grafo.salidas, vec![vec![1], vec![0]]);
    let solo = grafo_de_solapamiento(&piezas("AAA")).unwrap();
    assert_eq!(solo.salidas, vec![Vec::<usize>::new()]);
}

#[test]
fn solapamiento_rechaza_entradas_invalidas() {
    assert!(grafo_de_solapamiento(&[]).is_err());
    assert!(grafo_de_solapamiento(&piezas("A C")).is_err());
    assert!(grafo_de_solapamiento(&piezas("ACG AC")).is_err());
}

#[test]
fn solapamiento_rechaza_grafos_gigantes() {
    let mismo = vec![b"AAAA".as_slice(); 3000];
    let error = grafo_de_solapamiento(&mismo).unwrap_err();
    assert!(error.contains("aristas"), "{error}");
}

#[test]
fn de_bruijn_de_un_texto_pega_nodos_iguales() {
    let grafo = de_bruijn_de_texto(b"GATTACAGATT", 4).unwrap();
    let nombres: Vec<&[u8]> = grafo.nodos.iter().map(|nodo| nodo.as_slice()).collect();
    assert_eq!(nombres, piezas("GAT ATT TTA TAC ACA CAG AGA"));
    assert_eq!(
        grafo.salidas,
        vec![
            vec![1, 1],
            vec![2],
            vec![3],
            vec![4],
            vec![5],
            vec![6],
            vec![0]
        ]
    );
    assert_eq!(cantidad_de_aristas(&grafo), 8);
    assert_eq!(grados_de_entrada(&grafo), vec![1, 2, 1, 1, 1, 1, 1]);
    assert_eq!(balance_de_nodos(&grafo), vec![1, -1, 0, 0, 0, 0, 0]);
}

#[test]
fn de_bruijn_de_texto_valida_k() {
    assert!(de_bruijn_de_texto(b"GATTACA", 1).is_err());
    assert!(de_bruijn_de_texto(b"GATTACA", 0).is_err());
    assert!(de_bruijn_de_texto(b"GATTACA", 8).is_err());
    assert!(de_bruijn_de_texto(b"", 2).is_err());
    let entero = de_bruijn_de_texto(b"GATTACA", 7).unwrap();
    assert_eq!(entero.salidas, vec![vec![1], vec![]]);
}

#[test]
fn de_bruijn_de_kmeros_igual_que_del_texto() {
    let texto = b"GATTACAGATT";
    let desde_texto = de_bruijn_de_texto(texto, 4).unwrap();
    let desde_kmeros = de_bruijn_de_kmeros(&composicion(texto, 4)).unwrap();
    assert_eq!(desde_texto, desde_kmeros);
}

#[test]
fn de_bruijn_de_kmeros_desordenados() {
    let grafo = de_bruijn_de_kmeros(&piezas("TTAC GATT AGAT ATTA")).unwrap();
    let nombres: Vec<&[u8]> = grafo.nodos.iter().map(|nodo| nodo.as_slice()).collect();
    assert_eq!(nombres, piezas("TTA TAC GAT ATT AGA"));
    assert_eq!(
        grafo.salidas,
        vec![vec![1], vec![], vec![3], vec![0], vec![2]]
    );
}

#[test]
fn de_bruijn_de_kmeros_rechaza_entradas_invalidas() {
    assert!(de_bruijn_de_kmeros(&[]).is_err());
    assert!(de_bruijn_de_kmeros(&piezas("A C G")).is_err());
    assert!(de_bruijn_de_kmeros(&piezas("ACG ACGT")).is_err());
}

#[test]
fn lee_listas_de_adyacencia_con_y_sin_espacios() {
    let grafo = leer_lista_de_adyacencia("0 -> 1,2\n1->2\r\n2: 0\n\n3 ->\n").unwrap();
    assert_eq!(etiquetas(&grafo, &[0, 1, 2, 3]), ["0", "1", "2", "3"]);
    assert_eq!(grafo.salidas, vec![vec![1, 2], vec![2], vec![0], vec![]]);
    let con_espacios = leer_lista_de_adyacencia("a -> b c\na -> d").unwrap();
    assert_eq!(con_espacios.salidas[0], vec![1, 2, 3]);
}

#[test]
fn rechaza_listas_de_adyacencia_mal_escritas() {
    assert!(leer_lista_de_adyacencia("").is_err());
    assert!(leer_lista_de_adyacencia("   \n  ").is_err());
    assert!(leer_lista_de_adyacencia("a b").is_err());
    assert!(leer_lista_de_adyacencia(" -> b").is_err());
    assert!(leer_lista_de_adyacencia("a c -> b").is_err());
    let error = leer_lista_de_adyacencia("a -> b\nc d").unwrap_err();
    assert!(error.contains('2'), "{error}");
    let doble = leer_lista_de_adyacencia("a -> b\nb -> c -> d").unwrap_err();
    assert!(doble.contains('2'), "{doble}");
}

#[test]
fn ciclo_euleriano_de_un_grafo_balanceado() {
    let grafo = leer_lista_de_adyacencia("a -> b\nb -> c,d\nc -> a\nd -> b").unwrap();
    let ciclo = ciclo_euleriano(&grafo).unwrap();
    assert_eq!(etiquetas(&grafo, &ciclo), ["a", "b", "d", "b", "c", "a"]);
    assert!(usa_cada_arista_una_vez(&grafo, &ciclo));
    let triangulo = leer_lista_de_adyacencia("0 -> 1\n1 -> 2\n2 -> 0").unwrap();
    let vuelta = ciclo_euleriano(&triangulo).unwrap();
    assert_eq!(etiquetas(&triangulo, &vuelta), ["0", "1", "2", "0"]);
}

#[test]
fn ciclo_euleriano_con_bucles_y_aristas_repetidas() {
    let grafo = leer_lista_de_adyacencia("x -> x,y,y\ny -> x,x").unwrap();
    let ciclo = ciclo_euleriano(&grafo).unwrap();
    assert_eq!(ciclo.len(), 6);
    assert_eq!(ciclo.first(), ciclo.last());
    assert!(usa_cada_arista_una_vez(&grafo, &ciclo));
}

#[test]
fn sin_ciclo_euleriano_si_hay_nodos_desbalanceados() {
    let grafo = leer_lista_de_adyacencia("a -> b\na -> c").unwrap();
    let error = ciclo_euleriano(&grafo).unwrap_err();
    assert!(error.contains("«a»"), "{error}");
    assert!(camino_euleriano(&grafo).is_err());
}

#[test]
fn sin_ciclo_euleriano_si_las_aristas_no_se_tocan() {
    let grafo = leer_lista_de_adyacencia("a -> b\nb -> a\nc -> d\nd -> c").unwrap();
    assert!(ciclo_euleriano(&grafo).is_err());
    assert!(camino_euleriano(&grafo).is_err());
}

#[test]
fn sin_aristas_no_hay_recorrido() {
    let grafo = leer_lista_de_adyacencia("a ->").unwrap();
    assert!(ciclo_euleriano(&grafo).is_err());
    assert!(camino_euleriano(&grafo).is_err());
    assert!(recorrido_de_hierholzer(&grafo, 99).is_empty());
}

#[test]
fn camino_euleriano_empieza_donde_sobra_una_salida() {
    let grafo = leer_lista_de_adyacencia("1 -> 2\n2 -> 3\n3 -> 1,4").unwrap();
    let camino = camino_euleriano(&grafo).unwrap();
    assert_eq!(etiquetas(&grafo, &camino), ["3", "1", "2", "3", "4"]);
    let ida_y_vuelta = leer_lista_de_adyacencia("0 -> 1\n1 -> 0").unwrap();
    let cerrado = camino_euleriano(&ida_y_vuelta).unwrap();
    assert_eq!(etiquetas(&ida_y_vuelta, &cerrado), ["0", "1", "0"]);
}

#[test]
fn camino_euleriano_rechaza_dos_inicios() {
    let grafo = leer_lista_de_adyacencia("a -> b\nc -> d").unwrap();
    let error = camino_euleriano(&grafo).unwrap_err();
    assert!(error.contains("«c»") || error.contains("«a»"), "{error}");
}

#[test]
fn ciclos_eulerianos_de_grafos_al_azar() {
    for semilla in 0..40 {
        let mut azar = Azar::nuevo(semilla);
        let nodos = 2 + azar.entero(8);
        let pasos = 1 + azar.entero(40);
        let mut salidas = vec![Vec::new(); nodos];
        let mut actual = 0;
        for _ in 0..pasos {
            let siguiente = azar.entero(nodos);
            salidas[actual].push(siguiente);
            actual = siguiente;
        }
        salidas[actual].push(0);
        for (indice, destinos) in salidas.iter_mut().enumerate() {
            barajar(destinos, semilla * 100 + indice as u64);
        }
        let grafo = Grafo {
            nodos: (0..nodos)
                .map(|nodo| nodo.to_string().into_bytes())
                .collect(),
            salidas,
        };
        let ciclo = ciclo_euleriano(&grafo).unwrap();
        assert_eq!(ciclo.first(), ciclo.last());
        assert!(usa_cada_arista_una_vez(&grafo, &ciclo), "semilla {semilla}");
        let camino = camino_euleriano(&grafo).unwrap();
        assert!(
            usa_cada_arista_una_vez(&grafo, &camino),
            "semilla {semilla}"
        );
    }
}

#[test]
fn reconstruye_desde_kmeros_desordenados() {
    let kmeros = piezas("TACA GATT CAGA ATTA AGAT ACAG GATT TTAC");
    assert_eq!(reconstruir_desde_kmeros(&kmeros).unwrap(), b"GATTACAGATT");
}

#[test]
fn reconstruccion_rechaza_kmeros_que_no_forman_camino() {
    assert!(reconstruir_desde_kmeros(&piezas("AAC GGT")).is_err());
    assert!(reconstruir_desde_kmeros(&[]).is_err());
    assert_eq!(reconstruir_desde_kmeros(&piezas("AC")).unwrap(), b"AC");
}

#[test]
fn reconstruye_genomas_al_azar_sin_repeticiones() {
    for semilla in 0..10 {
        let genoma = genoma_al_azar(300, semilla);
        let mut kmeros = composicion(&genoma, 20);
        barajar(&mut kmeros, semilla + 1000);
        assert_eq!(reconstruir_desde_kmeros(&kmeros).unwrap(), genoma);
    }
}

#[test]
fn con_repeticiones_conserva_la_composicion() {
    for semilla in 0..10 {
        let genoma = genoma_al_azar(200, semilla);
        let mut kmeros = composicion(&genoma, 4);
        barajar(&mut kmeros, semilla + 2000);
        let armado = reconstruir_desde_kmeros(&kmeros).unwrap();
        assert_eq!(armado.len(), genoma.len());
        assert_eq!(
            composicion_ordenada(&armado, 4),
            composicion_ordenada(&genoma, 4)
        );
    }
}

#[test]
fn cadena_universal_de_ejemplos_chicos() {
    assert_eq!(cadena_universal_circular(1).unwrap(), b"01");
    assert_eq!(cadena_universal_circular(2).unwrap(), b"0011");
    assert_eq!(cadena_universal_circular(3).unwrap(), b"00010111");
}

#[test]
fn cadena_universal_contiene_cada_kmero_una_vez() {
    for k in 1..=12 {
        let cadena = cadena_universal_circular(k).unwrap();
        assert_eq!(cadena.len(), 1 << k);
        assert!(cadena.iter().all(|&letra| letra == b'0' || letra == b'1'));
        let mut vuelta = cadena.clone();
        vuelta.extend_from_slice(&cadena[..k - 1]);
        let mut vistos = HashSet::new();
        for kmero in vuelta.windows(k) {
            assert!(vistos.insert(kmero.to_vec()), "k = {k}");
        }
        assert_eq!(vistos.len(), 1 << k);
    }
}

#[test]
fn cadena_universal_del_k_maximo() {
    let k = MAXIMO_K_UNIVERSAL;
    let cadena = cadena_universal_circular(k).unwrap();
    assert_eq!(cadena.len(), 1 << k);
    let mascara = (1usize << k) - 1;
    let mut visto = vec![false; 1 << k];
    let mut numero = 0;
    for posicion in 0..cadena.len() + k - 1 {
        let bit = (cadena[posicion % cadena.len()] - b'0') as usize;
        numero = ((numero << 1) | bit) & mascara;
        if posicion + 1 >= k {
            assert!(!visto[numero], "repetido en {posicion}");
            visto[numero] = true;
        }
    }
    assert!(visto.iter().all(|&marcado| marcado));
}

#[test]
fn cadena_universal_rechaza_k_fuera_de_rango() {
    assert!(cadena_universal_circular(0).is_err());
    assert!(cadena_universal_circular(MAXIMO_K_UNIVERSAL + 1).is_err());
    assert!(cadena_universal_circular(usize::MAX).is_err());
}

#[test]
fn lee_pares_de_lecturas() {
    let pares = leer_pares("acg|GCA CGT | CAA\nGTT|AAG\r\n").unwrap();
    assert_eq!(pares.len(), 3);
    assert_eq!(pares[0].primera, b"ACG");
    assert_eq!(pares[0].segunda, b"GCA");
    assert_eq!(pares[1].segunda, b"CAA");
    assert_eq!(pares[2].primera, b"GTT");
}

#[test]
fn rechaza_pares_mal_escritos() {
    assert!(leer_pares("").is_err());
    assert!(leer_pares("ACG").is_err());
    assert!(leer_pares("ACG|").is_err());
    assert!(leer_pares("A|C|G").is_err());
    assert!(leer_pares("AC1|GCA").is_err());
}

#[test]
fn deletrea_un_camino_con_huecos() {
    let pares = pares_del_genoma_con_huecos();
    assert_eq!(
        deletrear_camino_con_huecos(&pares, 2).unwrap(),
        GENOMA_CON_HUECOS
    );
    assert_eq!(pares_de(GENOMA_CON_HUECOS, 3, 2), pares);
}

#[test]
fn camino_con_huecos_justo_sin_superposicion() {
    let pares = pares_del_genoma_con_huecos();
    assert_eq!(
        deletrear_camino_con_huecos(&pares[..3], 2).unwrap(),
        &GENOMA_CON_HUECOS[..10]
    );
}

#[test]
fn camino_con_huecos_rechaza_lo_que_no_encaja() {
    let pares = pares_del_genoma_con_huecos();
    assert!(deletrear_camino_con_huecos(&pares[..2], 2).is_err());
    assert!(deletrear_camino_con_huecos(&pares, 1).is_err());
    assert!(deletrear_camino_con_huecos(&pares, usize::MAX).is_err());
    assert!(deletrear_camino_con_huecos(&[], 2).is_err());
    let mut cambiado = pares.clone();
    cambiado[0].segunda = b"TCA".to_vec();
    assert!(deletrear_camino_con_huecos(&cambiado, 2).is_err());
    let mut desparejo = pares.clone();
    desparejo[0].segunda = b"GC".to_vec();
    assert!(deletrear_camino_con_huecos(&desparejo, 2).is_err());
}

#[test]
fn de_bruijn_pareado_une_prefijos_y_sufijos() {
    let grafo = de_bruijn_de_pares(&pares_del_genoma_con_huecos()).unwrap();
    assert_eq!(grafo.nodos[0], b"AC|GC");
    assert_eq!(grafo.nodos[1], b"CG|CA");
    assert_eq!(cantidad_de_aristas(&grafo), 5);
    assert_eq!(grafo.nodos.len(), 6);
}

#[test]
fn reconstruye_desde_pares_desordenados() {
    let mut pares = pares_del_genoma_con_huecos();
    barajar(&mut pares, 3);
    let armado = reconstruir_desde_pares(&pares, 2, 7).unwrap();
    assert_eq!(armado.genoma, GENOMA_CON_HUECOS);
    assert_eq!(armado.intentos, 1);
}

#[test]
fn reconstruye_desde_pares_cuando_el_grafo_es_un_ciclo() {
    let genoma = b"GACAGGGTGG";
    let pares = pares_de(genoma, 2, 2);
    let armado = reconstruir_desde_pares(&pares, 2, 136).unwrap();
    let mut esperados = pares_de(genoma, 2, 2);
    let mut obtenidos = pares_de(&armado.genoma, 2, 2);
    esperados
        .sort_by(|una, otra| (&una.primera, &una.segunda).cmp(&(&otra.primera, &otra.segunda)));
    obtenidos
        .sort_by(|una, otra| (&una.primera, &una.segunda).cmp(&(&otra.primera, &otra.segunda)));
    assert_eq!(armado.genoma.len(), genoma.len());
    assert_eq!(obtenidos, esperados);
}

#[test]
fn reconstruccion_con_pares_rechaza_entradas_invalidas() {
    let pares = pares_del_genoma_con_huecos();
    assert!(reconstruir_desde_pares(&[], 2, 1).is_err());
    assert!(reconstruir_desde_pares(&pares[..2], 2, 1).is_err());
    assert!(reconstruir_desde_pares(&pares, usize::MAX, 1).is_err());
    assert!(reconstruir_desde_pares(&pares, 5, 1).is_err());
}

#[test]
fn pares_al_azar_dan_la_misma_composicion() {
    let mut reintentos = 0;
    for semilla in 0..30 {
        let mut genoma = genoma_al_azar(120, semilla);
        let repeticion = genoma[10..16].to_vec();
        genoma.splice(60..66, repeticion);
        let mut pares = pares_de(&genoma, 4, 3);
        barajar(&mut pares, semilla + 500);
        match reconstruir_desde_pares(&pares, 3, semilla) {
            Ok(armado) => {
                assert_eq!(armado.genoma.len(), genoma.len());
                assert_eq!(
                    pares_ordenados(&pares_de(&armado.genoma, 4, 3)),
                    pares_ordenados(&pares)
                );
                if armado.intentos > 1 {
                    reintentos += 1;
                }
            }
            Err(error) => assert!(error.contains("intent"), "{error}"),
        }
    }
    assert!(reintentos > 0);
}

#[test]
fn caminos_no_ramificados_y_ciclos_aislados() {
    let grafo = leer_lista_de_adyacencia(
        "p -> q\nq -> r\nr -> s,t\ns -> u\nt -> u\nu -> v\nw -> x\nx -> w",
    )
    .unwrap();
    let caminos: Vec<Vec<String>> = caminos_no_ramificados(&grafo)
        .iter()
        .map(|camino| etiquetas(&grafo, camino))
        .collect();
    assert_eq!(
        caminos,
        vec![
            vec!["p", "q", "r"],
            vec!["r", "s", "u"],
            vec!["r", "t", "u"],
            vec!["u", "v"],
            vec!["w", "x", "w"]
        ]
    );
}

#[test]
fn caminos_no_ramificados_con_bucle_propio() {
    let grafo = leer_lista_de_adyacencia("a -> a").unwrap();
    assert_eq!(caminos_no_ramificados(&grafo), vec![vec![0, 0]]);
    let vacio = leer_lista_de_adyacencia("a ->").unwrap();
    assert!(caminos_no_ramificados(&vacio).is_empty());
}

#[test]
fn contigs_se_cortan_en_las_repeticiones() {
    let contigs = generar_contigs(&composicion(b"GATTACAGATT", 4)).unwrap();
    let textos: Vec<&[u8]> = contigs.iter().map(|contig| contig.as_slice()).collect();
    assert_eq!(textos, piezas("GATT GATT ATTACAGAT"));
}

#[test]
fn contigs_de_un_genoma_sin_repeticiones_y_de_uno_circular() {
    let genoma = genoma_al_azar(150, 77);
    let mut kmeros = composicion(&genoma, 25);
    barajar(&mut kmeros, 78);
    assert_eq!(generar_contigs(&kmeros).unwrap(), vec![genoma]);
    let circular = generar_contigs(&piezas("CGT TAC ACG GTA")).unwrap();
    assert_eq!(circular, vec![b"CGTACG".to_vec()]);
    assert!(generar_contigs(&[]).is_err());
}

#[test]
fn n50_de_contigs() {
    let medida = medir_n50(&[2, 3, 4, 5, 6, 7, 8, 9, 10]);
    assert_eq!(
        medida,
        MedidaN50 {
            n50: 8,
            l50: 3,
            total: 54
        }
    );
    let empatados = medir_n50(&[5, 5, 5, 5]);
    assert_eq!((empatados.n50, empatados.l50), (5, 2));
    assert_eq!(medir_n50(&[100]).n50, 100);
    assert_eq!(medir_n50(&[1, 1, 1, 1, 50]).l50, 1);
}

#[test]
fn n50_sin_contigs() {
    let vacio = MedidaN50 {
        n50: 0,
        l50: 0,
        total: 0,
    };
    assert_eq!(medir_n50(&[]), vacio);
    assert_eq!(medir_n50(&[0, 0]), vacio);
}

#[test]
fn n50_con_largos_enormes_no_se_desborda() {
    let medida = medir_n50(&[usize::MAX, usize::MAX, 3]);
    assert_eq!((medida.n50, medida.l50), (usize::MAX, 1));
    assert_eq!(medida.total, usize::MAX);
}
