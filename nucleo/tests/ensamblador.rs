use nucleo::azar::Azar;
use nucleo::ensamblador::*;
use nucleo::formatos::{Lectura, escribir_fastq};
use nucleo::secuencia::complemento_inverso;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

static CONTADOR: AtomicUsize = AtomicUsize::new(0);

fn ruta_temporal(extension: &str) -> PathBuf {
    let numero = CONTADOR.fetch_add(1, Ordering::Relaxed);
    std::env::temp_dir().join(format!(
        "ensamblador-prueba-{}-{numero}.{extension}",
        std::process::id()
    ))
}

fn escribir_lecturas(secuencias: &[Vec<u8>]) -> PathBuf {
    let ruta = ruta_temporal("fastq");
    let mut destino = BufWriter::new(File::create(&ruta).unwrap());
    for (numero, secuencia) in secuencias.iter().enumerate() {
        let lectura = Lectura {
            nombre: format!("l{numero}").into_bytes(),
            secuencia: secuencia.clone(),
            calidad: vec![b'I'; secuencia.len()],
        };
        escribir_fastq(&mut destino, &lectura).unwrap();
    }
    destino.flush().unwrap();
    ruta
}

fn genoma_al_azar(largo: usize, azar: &mut Azar) -> Vec<u8> {
    (0..largo).map(|_| b"ACGT"[azar.entero(4)]).collect()
}

fn simular(
    genoma: &[u8],
    largo: usize,
    cobertura: usize,
    error: f64,
    azar: &mut Azar,
) -> Vec<Vec<u8>> {
    let cantidad = genoma.len() * cobertura / largo;
    let mut lecturas = Vec::with_capacity(cantidad);
    for _ in 0..cantidad {
        let inicio = azar.entero(genoma.len() - largo + 1);
        let mut lectura = genoma[inicio..inicio + largo].to_vec();
        if azar.entero(2) == 1 {
            lectura = complemento_inverso(&lectura);
        }
        for base in lectura.iter_mut() {
            if azar.decimal() < error {
                let original = *base;
                while *base == original {
                    *base = b"ACGT"[azar.entero(4)];
                }
            }
        }
        lecturas.push(lectura);
    }
    lecturas
}

fn texto(letras: &[u8]) -> String {
    String::from_utf8(letras.to_vec()).unwrap()
}

fn canonico_texto(kmero: &str) -> String {
    let inverso = texto(&complemento_inverso(kmero.as_bytes()));
    if inverso.as_str() < kmero {
        inverso
    } else {
        kmero.to_string()
    }
}

fn conteo_ingenuo(lecturas: &[Vec<u8>], k: usize) -> HashMap<String, u32> {
    let mut cuentas = HashMap::new();
    for lectura in lecturas {
        let lectura = texto(lectura).to_ascii_uppercase();
        for trozo in lectura.split(|letra: char| !"ACGT".contains(letra)) {
            if trozo.len() < k {
                continue;
            }
            for inicio in 0..=trozo.len() - k {
                *cuentas
                    .entry(canonico_texto(&trozo[inicio..inicio + k]))
                    .or_insert(0) += 1;
            }
        }
    }
    cuentas
}

fn tabla_desde(conteos: &HashMap<String, u32>, k: usize) -> Kmeros {
    let mut kmeros = Kmeros::nuevo(k);
    for (kmero, &veces) in conteos {
        let mut canonico = None;
        recorrer_kmeros(kmero.as_bytes(), k, |_, leido| {
            canonico = Some(leido.canonico())
        });
        kmeros.agregar(canonico.unwrap(), veces);
    }
    kmeros
}

fn kmeros_de_texto(secuencia: &[u8], k: usize) -> Vec<String> {
    let mut lista = Vec::new();
    recorrer_kmeros(secuencia, k, |_, kmero| {
        lista.push(texto(&letras_de(kmero.canonico(), k)))
    });
    lista
}

fn sucesores_ingenuos(conjunto: &HashSet<String>, kmero: &str) -> Vec<String> {
    let mut lista = Vec::new();
    for base in ["A", "C", "G", "T"] {
        let siguiente = format!("{}{base}", &kmero[1..]);
        if conjunto.contains(&canonico_texto(&siguiente)) {
            lista.push(siguiente);
        }
    }
    lista
}

fn predecesores_ingenuos(conjunto: &HashSet<String>, kmero: &str) -> Vec<String> {
    let inverso = texto(&complemento_inverso(kmero.as_bytes()));
    sucesores_ingenuos(conjunto, &inverso)
        .into_iter()
        .map(|otro| texto(&complemento_inverso(otro.as_bytes())))
        .collect()
}

fn particion_ingenua(conjunto: &HashSet<String>) -> Vec<Vec<String>> {
    let mut ordenados: Vec<&String> = conjunto.iter().collect();
    ordenados.sort();
    let mut vistos: HashSet<String> = HashSet::new();
    let mut grupos = Vec::new();
    for semilla in ordenados {
        if !vistos.insert(semilla.clone()) {
            continue;
        }
        let mut grupo = vec![semilla.clone()];
        for inicio in [
            semilla.clone(),
            texto(&complemento_inverso(semilla.as_bytes())),
        ] {
            let mut actual = inicio;
            loop {
                let siguientes = sucesores_ingenuos(conjunto, &actual);
                if siguientes.len() != 1 {
                    break;
                }
                let siguiente = siguientes[0].clone();
                if predecesores_ingenuos(conjunto, &siguiente).len() != 1 {
                    break;
                }
                if !vistos.insert(canonico_texto(&siguiente)) {
                    break;
                }
                grupo.push(canonico_texto(&siguiente));
                actual = siguiente;
            }
        }
        grupo.sort();
        grupos.push(grupo);
    }
    grupos.sort();
    grupos
}

fn particion_de(grafo: &GrafoDeUnitigs) -> Vec<Vec<String>> {
    let mut grupos: Vec<Vec<String>> = grafo
        .unitigs
        .iter()
        .map(|unitig| {
            let mut grupo = kmeros_de_texto(&unitig.secuencia, grafo.k);
            grupo.sort();
            grupo
        })
        .collect();
    grupos.sort();
    grupos
}

fn kmeros_con_errores(azar: &mut Azar, k: usize) -> HashMap<String, u32> {
    let genoma = genoma_al_azar(60 + azar.entero(200), azar);
    let mut lecturas = vec![genoma.clone()];
    for _ in 0..azar.entero(6) {
        let inicio = azar.entero(genoma.len() - k);
        let fin = (inicio + k + azar.entero(2 * k)).min(genoma.len());
        let mut copia = genoma[inicio..fin].to_vec();
        let posicion = azar.entero(copia.len());
        copia[posicion] = b"ACGT"[azar.entero(4)];
        lecturas.push(copia);
    }
    if azar.entero(3) == 0 {
        let inicio = azar.entero(genoma.len() / 2);
        let repeticion = genoma[inicio..inicio + k + azar.entero(k)].to_vec();
        lecturas.push(
            [
                genoma_al_azar(30, azar),
                repeticion,
                genoma_al_azar(30, azar),
            ]
            .concat(),
        );
    }
    conteo_ingenuo(&lecturas, k)
}

#[test]
fn codifica_kmeros_canonicos_a_mano() {
    let mut vistos = Vec::new();
    recorrer_kmeros(b"ACGTT", 5, |posicion, kmero| {
        vistos.push((posicion, kmero))
    });
    assert_eq!(vistos.len(), 1);
    let (posicion, kmero) = vistos[0];
    assert_eq!(posicion, 0);
    assert_eq!(kmero.directo, 0b00_01_10_11_11);
    assert_eq!(kmero.reverso, 0b00_00_01_10_11);
    assert_eq!(kmero.canonico(), 27);
    assert_eq!(letras_de(kmero.canonico(), 5), b"AACGT".to_vec());
    assert_eq!(desde_canonico(27, 5).reverso, kmero.directo);
    assert_eq!(kmero.girado().canonico(), 27);
}

#[test]
fn las_n_cortan_los_kmeros_y_las_minusculas_cuentan() {
    let mut posiciones = Vec::new();
    recorrer_kmeros(b"acgtaNACGTAGG", 5, |posicion, _| posiciones.push(posicion));
    assert_eq!(posiciones, vec![0, 6, 7, 8]);
    let mut nada = 0;
    recorrer_kmeros(b"ACGTACGT", 0, |_, _| nada += 1);
    recorrer_kmeros(b"ACGTACGT", 40, |_, _| nada += 1);
    assert_eq!(nada, 0);
}

#[test]
fn avanzar_mantiene_las_dos_hebras() {
    let mut azar = Azar::nuevo(7);
    for _ in 0..200 {
        let k = 5 + 2 * azar.entero(14);
        let secuencia = genoma_al_azar(k + 20, &mut azar);
        recorrer_kmeros(&secuencia, k, |posicion, kmero| {
            let directo = &secuencia[posicion..posicion + k];
            assert_eq!(letras_de(kmero.directo, k), directo.to_vec());
            assert_eq!(letras_de(kmero.reverso, k), complemento_inverso(directo));
            assert_eq!(
                desde_canonico(kmero.canonico(), k).canonico(),
                kmero.canonico()
            );
        });
    }
}

#[test]
fn el_valle_es_el_primer_minimo() {
    let mut barras = vec![0; ULTIMA_BARRA + 1];
    for (cobertura, veces) in [100, 10, 5, 8, 20, 30, 12, 3].into_iter().enumerate() {
        barras[cobertura + 1] = veces;
    }
    assert_eq!(primer_valle(&barras), Some(3));
    let umbral = elegir_umbral(&barras, 0);
    assert_eq!(umbral.minimo, 3);
    assert!(umbral.automatico);
    let pedido = elegir_umbral(&barras, 5);
    assert_eq!(pedido.minimo, 5);
    assert_eq!(pedido.valle, Some(3));
    assert!(!pedido.automatico);
    let bajando: Vec<usize> = (0..=ULTIMA_BARRA).map(|c| 1000 - c).collect();
    assert_eq!(primer_valle(&bajando), None);
    assert_eq!(elegir_umbral(&bajando, 0).minimo, MINIMO_SIN_VALLE);
    assert_eq!(primer_valle(&[]), None);
    let mut sin_errores = vec![0; ULTIMA_BARRA + 1];
    sin_errores[1] = 2;
    sin_errores[2] = 7;
    assert_eq!(primer_valle(&sin_errores), Some(1));
}

#[test]
fn el_conteo_paralelo_coincide_con_el_ingenuo() {
    let mut azar = Azar::nuevo(11);
    for _ in 0..25 {
        let k = 5 + 2 * azar.entero(8);
        let mut lecturas = Vec::new();
        for _ in 0..azar.entero(3000) + 1 {
            let mut lectura = genoma_al_azar(azar.entero(60) + 1, &mut azar);
            if azar.entero(5) == 0 && !lectura.is_empty() {
                let posicion = azar.entero(lectura.len());
                lectura[posicion] = b'N';
            }
            lecturas.push(lectura);
        }
        let ruta = escribir_lecturas(&lecturas);
        let esperado = conteo_ingenuo(&lecturas, k);
        match contar_kmeros(&ruta, k) {
            Ok(conteo) => {
                assert_eq!(conteo.lecturas, lecturas.len());
                assert_eq!(conteo.kmeros.cantidad(), esperado.len());
                assert_eq!(
                    conteo.kmeros_leidos,
                    esperado.values().map(|&v| v as usize).sum::<usize>()
                );
                for (kmero, &veces) in &esperado {
                    let mut canonico = 0;
                    recorrer_kmeros(kmero.as_bytes(), k, |_, leido| canonico = leido.canonico());
                    assert_eq!(conteo.kmeros.cuenta(canonico), veces);
                }
            }
            Err(_) => assert!(esperado.is_empty()),
        }
        std::fs::remove_file(ruta).ok();
    }
}

#[test]
fn la_compactacion_coincide_con_la_ingenua() {
    let mut azar = Azar::nuevo(23);
    for _ in 0..150 {
        let k = 5 + 2 * azar.entero(4);
        let conteos = kmeros_con_errores(&mut azar, k);
        let kmeros = tabla_desde(&conteos, k);
        let grafo = compactar(&kmeros);
        let conjunto: HashSet<String> = conteos.keys().cloned().collect();
        assert_eq!(particion_de(&grafo), particion_ingenua(&conjunto));
        for (numero, unitig) in grafo.unitigs.iter().enumerate() {
            assert!(unitig.secuencia <= complemento_inverso(&unitig.secuencia));
            let cuentas: Vec<u32> = kmeros_de_texto(&unitig.secuencia, k)
                .iter()
                .map(|kmero| conteos[kmero])
                .collect();
            let media =
                cuentas.iter().map(|&c| c as u64).sum::<u64>() as f64 / cuentas.len() as f64;
            assert_eq!(unitig.cobertura, media);
            for lado in [false, true] {
                let propio = Lado {
                    unitig: numero,
                    al_reves: lado,
                };
                for &vecino in grafo.salidas(propio) {
                    assert!(grafo.salidas(vecino.opuesto()).contains(&propio.opuesto()));
                }
            }
        }
        for par in grafo.unitigs.windows(2) {
            assert!(par[0].secuencia.len() >= par[1].secuencia.len());
        }
    }
}

#[test]
fn las_aristas_unen_extremos_que_se_solapan() {
    let mut azar = Azar::nuevo(5);
    for _ in 0..60 {
        let k = 5 + 2 * azar.entero(4);
        let kmeros = tabla_desde(&kmeros_con_errores(&mut azar, k), k);
        let grafo = compactar(&kmeros);
        let orientada = |lado: Lado| {
            let secuencia = &grafo.unitigs[lado.unitig].secuencia;
            if lado.al_reves {
                complemento_inverso(secuencia)
            } else {
                secuencia.clone()
            }
        };
        for numero in 0..grafo.unitigs.len() {
            for al_reves in [false, true] {
                let lado = Lado {
                    unitig: numero,
                    al_reves,
                };
                let propia = orientada(lado);
                let cola = &propia[propia.len() - (k - 1)..];
                let esperados: usize = (0..4)
                    .filter(|&base| {
                        let mut siguiente = cola.to_vec();
                        siguiente.push(b"ACGT"[base]);
                        let mut canonico = 0;
                        recorrer_kmeros(&siguiente, k, |_, kmero| canonico = kmero.canonico());
                        kmeros.contiene(canonico)
                    })
                    .count();
                assert_eq!(grafo.salidas(lado).len(), esperados);
                for &vecino in grafo.salidas(lado) {
                    let otra = orientada(vecino);
                    assert_eq!(&otra[..k - 1], cola);
                }
            }
        }
    }
}

fn genoma_sin_repeticiones(largo: usize, k: usize, azar: &mut Azar) -> Vec<u8> {
    loop {
        let genoma = genoma_al_azar(largo, azar);
        let mut distintos = HashSet::new();
        let sin_ramas = genoma.windows(k - 1).all(|pieza| {
            let inverso = complemento_inverso(pieza);
            inverso != pieza && distintos.insert(pieza.to_vec().min(inverso))
        });
        if sin_ramas {
            return genoma;
        }
    }
}

fn agregar_texto(kmeros: &mut Kmeros, secuencia: &[u8], veces: u32) {
    let k = kmeros.k;
    recorrer_kmeros(secuencia, k, |_, kmero| {
        kmeros.agregar(kmero.canonico(), veces)
    });
}

fn solo_unitig(limpieza: &Limpieza) -> Vec<u8> {
    assert_eq!(limpieza.grafo.unitigs.len(), 1);
    limpieza.grafo.unitigs[0].secuencia.clone()
}

fn orientado(secuencia: &[u8]) -> Vec<u8> {
    secuencia.to_vec().min(complemento_inverso(secuencia))
}

#[test]
fn quita_una_punta_de_baja_cobertura() {
    let k = 11;
    let mut azar = Azar::nuevo(3);
    let genoma = genoma_sin_repeticiones(300, k, &mut azar);
    let mut kmeros = Kmeros::nuevo(k);
    agregar_texto(&mut kmeros, &genoma, 10);
    let error = 150;
    let mut punta = genoma[error + 1 - k..error + 4].to_vec();
    punta[k - 1] = if genoma[error] == b'A' { b'C' } else { b'A' };
    agregar_texto(&mut kmeros, &punta, 2);
    let antes = compactar(&kmeros);
    assert_eq!(antes.unitigs.len(), 3);
    let puntas = buscar_puntas(&antes);
    assert_eq!(puntas.len(), 1);
    assert_eq!(puntas[0].secuencia, orientado(&punta));
    assert_eq!(puntas[0].cobertura, 2.0);
    assert_eq!(puntas[0].cobertura_rival, 10.0);
    let limpieza = limpiar(&mut kmeros);
    assert_eq!(limpieza.puntas.len(), 1);
    assert_eq!(limpieza.burbujas.len(), 0);
    assert_eq!(limpieza.rondas, 1);
    assert_eq!(solo_unitig(&limpieza), orientado(&genoma));
}

#[test]
fn no_quita_puntas_largas_ni_mas_cubiertas() {
    let k = 11;
    let mut azar = Azar::nuevo(4);
    let genoma = genoma_sin_repeticiones(300, k, &mut azar);
    let error = 150;
    let cambiada = if genoma[error] == b'A' { b'C' } else { b'A' };
    let mut corta = genoma[error + 1 - k..error + 4].to_vec();
    corta[k - 1] = cambiada;
    let mut kmeros = Kmeros::nuevo(k);
    agregar_texto(&mut kmeros, &genoma, 10);
    agregar_texto(&mut kmeros, &corta, 12);
    assert_eq!(compactar(&kmeros).unitigs.len(), 3);
    assert!(buscar_puntas(&compactar(&kmeros)).is_empty());
    let mut larga = genoma[error + 1 - k..error].to_vec();
    larga.push(cambiada);
    larga.extend(genoma_al_azar(2 * k, &mut azar));
    let mut otra = Kmeros::nuevo(k);
    agregar_texto(&mut otra, &genoma, 10);
    agregar_texto(&mut otra, &larga, 2);
    let grafo = compactar(&otra);
    assert!(
        grafo
            .unitigs
            .iter()
            .any(|u| u.secuencia == orientado(&larga))
    );
    assert!(buscar_puntas(&grafo).is_empty());
}

#[test]
fn revienta_una_burbuja_y_conserva_la_rama_fuerte() {
    let k = 11;
    let mut azar = Azar::nuevo(9);
    let genoma = genoma_sin_repeticiones(400, k, &mut azar);
    let mut kmeros = Kmeros::nuevo(k);
    agregar_texto(&mut kmeros, &genoma, 10);
    let snp = 200;
    let mut rama = genoma[snp + 1 - k..snp + k].to_vec();
    rama[k - 1] = if genoma[snp] == b'G' { b'T' } else { b'G' };
    agregar_texto(&mut kmeros, &rama, 3);
    let antes = compactar(&kmeros);
    assert_eq!(antes.unitigs.len(), 4);
    let burbujas = buscar_burbujas(&antes);
    assert_eq!(burbujas.len(), 1);
    assert_eq!(burbujas[0].quitada, orientado(&rama));
    assert_eq!(burbujas[0].cobertura_quitada, 3.0);
    assert_eq!(burbujas[0].cobertura_conservada, 10.0);
    assert_eq!(burbujas[0].unitigs_conservados, 1);
    let principal = genoma[snp + 1 - k..snp + k].to_vec();
    let esperada = if orientado(&rama) == rama {
        principal
    } else {
        complemento_inverso(&principal)
    };
    assert_eq!(burbujas[0].conservada, esperada);
    let limpieza = limpiar(&mut kmeros);
    assert_eq!(limpieza.burbujas.len(), 1);
    assert_eq!(limpieza.puntas.len(), 0);
    assert_eq!(solo_unitig(&limpieza), orientado(&genoma));
}

#[test]
fn una_burbuja_con_ramas_muy_distintas_no_se_revienta() {
    let k = 11;
    let mut azar = Azar::nuevo(10);
    let genoma = genoma_sin_repeticiones(400, k, &mut azar);
    let snp = 200;
    let mut rama = genoma[snp + 1 - k..snp + k].to_vec();
    rama[k - 1] = if genoma[snp] == b'G' { b'T' } else { b'G' };
    let relleno = genoma_al_azar(12, &mut azar);
    rama.splice(k..k, relleno);
    let mut kmeros = Kmeros::nuevo(k);
    agregar_texto(&mut kmeros, &genoma, 10);
    agregar_texto(&mut kmeros, &rama, 3);
    assert!(!largos_parecidos(2 * k - 1, rama.len()));
    assert!(buscar_burbujas(&compactar(&kmeros)).is_empty());
}

#[test]
fn ensambla_lecturas_con_errores_sin_copiarlos() {
    let k = 21;
    let mut azar = Azar::nuevo(2026);
    let genoma = genoma_sin_repeticiones(6000, k, &mut azar);
    let lecturas = simular(&genoma, 100, 40, 0.005, &mut azar);
    let ruta = escribir_lecturas(&lecturas);
    let resultado = ensamblar(&ruta, k, 0).unwrap();
    let sin_filtro = ensamblar(&ruta, k, 1).unwrap();
    let filtro_bajo = ensamblar(&ruta, k, 2).unwrap();
    std::fs::remove_file(&ruta).ok();
    assert_eq!(resultado.lecturas, lecturas.len());
    assert!(resultado.umbral.automatico);
    let valle = resultado.umbral.valle.unwrap();
    assert!((2..=8).contains(&valle), "valle {valle}");
    assert_eq!(sin_filtro.umbral.minimo, 1);
    assert!(
        sin_filtro.puntas.len() > 200,
        "puntas {}",
        sin_filtro.puntas.len()
    );
    assert!(
        sin_filtro.burbujas.len() > 200,
        "burbujas {}",
        sin_filtro.burbujas.len()
    );
    assert!(sin_filtro.contigs.len() * 50 < sin_filtro.unitigs_iniciales);
    assert!(filtro_bajo.puntas.len() + filtro_bajo.burbujas.len() > 5);
    assert_eq!(filtro_bajo.contigs.len(), 1);
    assert!(filtro_bajo.contigs[0].secuencia.len() >= genoma.len() - 100);
    let verdaderos: HashSet<String> = kmeros_de_texto(&genoma, k).into_iter().collect();
    for contig in &resultado.contigs {
        for kmero in kmeros_de_texto(&contig.secuencia, k) {
            assert!(verdaderos.contains(&kmero));
        }
    }
    let mayor = &resultado.contigs[0].secuencia;
    assert!(mayor.len() >= genoma.len() - 200, "mayor {}", mayor.len());
    let en_el_genoma = |pieza: &[u8]| genoma.windows(pieza.len()).any(|ventana| ventana == pieza);
    assert!(en_el_genoma(mayor) || en_el_genoma(&complemento_inverso(mayor)));
    assert_eq!(resultado.medida.n50, mayor.len());
    assert_eq!(
        resultado.medida.total,
        resultado
            .contigs
            .iter()
            .map(|c| c.secuencia.len())
            .sum::<usize>()
    );
    assert!(resultado.kmeros_solidos < resultado.kmeros_distintos);
}

#[test]
fn sin_limpiar_hay_muchos_mas_pedazos() {
    let k = 21;
    let mut azar = Azar::nuevo(77);
    let genoma = genoma_sin_repeticiones(3000, k, &mut azar);
    let lecturas = simular(&genoma, 100, 30, 0.01, &mut azar);
    let ruta = escribir_lecturas(&lecturas);
    let crudos = largos_sin_limpiar(&ruta, k).unwrap();
    let limpio = ensamblar(&ruta, k, 0).unwrap();
    std::fs::remove_file(&ruta).ok();
    assert!(crudos.len() > 10 * limpio.contigs.len().max(1));
}

#[test]
fn lee_gzip_igual_que_texto() {
    let k = 15;
    let mut azar = Azar::nuevo(31);
    let genoma = genoma_sin_repeticiones(2000, k, &mut azar);
    let lecturas = simular(&genoma, 80, 25, 0.004, &mut azar);
    let plano = escribir_lecturas(&lecturas);
    let comprimido = ruta_temporal("fastq.gz");
    let mut codificador = flate2::write::GzEncoder::new(
        File::create(&comprimido).unwrap(),
        flate2::Compression::fast(),
    );
    codificador
        .write_all(&std::fs::read(&plano).unwrap())
        .unwrap();
    codificador.finish().unwrap();
    let uno = ensamblar(&plano, k, 0).unwrap();
    let otro = ensamblar(&comprimido, k, 0).unwrap();
    std::fs::remove_file(&plano).ok();
    std::fs::remove_file(&comprimido).ok();
    assert_eq!(uno.contigs, otro.contigs);
    assert_eq!(uno.histograma, otro.histograma);
}

#[test]
fn rechaza_entradas_malas_en_espanol() {
    let mut azar = Azar::nuevo(1);
    let lecturas = vec![genoma_al_azar(50, &mut azar); 4];
    let ruta = escribir_lecturas(&lecturas);
    for k in [0, 3, 4, 30, 33, 99] {
        let error = ensamblar(&ruta, k, 0).err().unwrap();
        assert!(error.contains("impar"), "{error}");
    }
    assert!(ensamblar(&ruta, 31, COBERTURA_MAXIMA_PEDIDA + 1).is_err());
    let vacio = ensamblar(&ruta, 31, COBERTURA_MAXIMA_PEDIDA).unwrap();
    assert!(vacio.contigs.is_empty());
    assert_eq!(vacio.medida.n50, 0);
    let cortas = escribir_lecturas(&[
        b"ACGTACGT".to_vec(),
        b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN".to_vec(),
    ]);
    assert!(ensamblar(&cortas, 21, 0).err().unwrap().contains("k-meros"));
    let vacio = ruta_temporal("fastq");
    File::create(&vacio).unwrap();
    assert!(
        ensamblar(&vacio, 21, 0)
            .err()
            .unwrap()
            .contains("ninguna lectura")
    );
    let fasta = ruta_temporal("fasta");
    std::fs::write(&fasta, ">s\nACGTACGTACGT\n").unwrap();
    assert!(ensamblar(&fasta, 5, 0).is_err());
    let cortado = ruta_temporal("fastq");
    std::fs::write(&cortado, "@l\nACGTACGTAC\n+\n").unwrap();
    assert!(ensamblar(&cortado, 5, 0).is_err());
    assert!(ensamblar(&ruta_temporal("no-existe"), 21, 0).is_err());
    for archivo in [ruta, cortas, vacio, fasta, cortado] {
        std::fs::remove_file(archivo).ok();
    }
}

#[test]
fn escribe_contigs_en_fasta() {
    let contigs = vec![
        Contig {
            secuencia: b"ACGT".repeat(40),
            cobertura: 12.34,
        },
        Contig {
            secuencia: b"GGA".to_vec(),
            cobertura: 3.0,
        },
    ];
    let mut destino = Vec::new();
    escribir_contigs(&mut destino, &contigs).unwrap();
    let texto = String::from_utf8(destino).unwrap();
    let lineas: Vec<&str> = texto.lines().collect();
    assert_eq!(lineas[0], ">contig_1 largo=160 cobertura_kmeros=12.3");
    assert_eq!(lineas[1].len(), 70);
    assert_eq!(lineas[4], ">contig_2 largo=3 cobertura_kmeros=3.0");
    assert_eq!(lineas[5], "GGA");
}

#[test]
fn compara_contigs_con_la_referencia() {
    let k = 11;
    let mut azar = Azar::nuevo(8);
    let referencia = genoma_sin_repeticiones(1000, k, &mut azar);
    let mut mutado = referencia[100..400].to_vec();
    mutado[150] = if mutado[150] == b'A' { b'C' } else { b'A' };
    let contigs = vec![
        Contig {
            secuencia: complemento_inverso(&mutado),
            cobertura: 10.0,
        },
        Contig {
            secuencia: referencia[600..640].to_vec(),
            cobertura: 10.0,
        },
        Contig {
            secuencia: genoma_al_azar(30, &mut azar),
            cobertura: 1.0,
        },
    ];
    let indice = indexar_referencia(&referencia, &[(50, 450), (590, 700), (990, 2000)], k);
    assert_eq!(indice.regiones, vec![(50, 450), (590, 700), (990, 1000)]);
    assert_eq!(indice.posiciones(indice.kmeros[0].0).len(), 1);
    assert_eq!(indice.ajenos(&referencia[600..640]), 0);
    let comparacion = comparar_con_referencia(&contigs, &indice, 100);
    assert_eq!(comparacion.largo_referencia, 520);
    assert_eq!(comparacion.cubiertas, 299 + 40);
    assert_eq!(comparacion.cubiertas_por_largos, 299);
    assert_eq!(
        comparacion.huecos,
        vec![
            (50, 100),
            (250, 251),
            (400, 450),
            (590, 600),
            (640, 700),
            (990, 1000)
        ]
    );
    assert_eq!(comparacion.ubicaciones[0], Some((100, 400)));
    assert_eq!(comparacion.ubicaciones[1], Some((600, 640)));
    assert_eq!(comparacion.ubicaciones[2], None);
    assert_eq!(comparacion.kmeros_de_contigs, 290 + 30 + 20);
    assert!(comparacion.kmeros_ajenos >= k + 20);
    let invalido = indexar_referencia(&referencia, &[(0, 1000)], 4);
    assert!(invalido.kmeros.is_empty());
    assert_eq!(
        comparar_con_referencia(&contigs, &invalido, 100),
        Comparacion::default()
    );
}

#[test]
fn revienta_burbujas_que_se_cruzan() {
    let k = 11;
    let mut azar = Azar::nuevo(12);
    let genoma = genoma_sin_repeticiones(400, k, &mut azar);
    let mut kmeros = Kmeros::nuevo(k);
    agregar_texto(&mut kmeros, &genoma, 10);
    let mut ramas = Vec::new();
    for snp in [200, 205, 240] {
        let mut rama = genoma[snp + 1 - k..snp + k].to_vec();
        rama[k - 1] = if genoma[snp] == b'C' { b'G' } else { b'C' };
        agregar_texto(&mut kmeros, &rama, 2);
        ramas.push(rama);
    }
    let antes = compactar(&kmeros);
    assert!(antes.unitigs.len() > 6);
    let limpieza = limpiar(&mut kmeros);
    assert_eq!(limpieza.burbujas.len(), 3);
    assert!(limpieza.burbujas.iter().any(|b| b.unitigs_conservados > 1));
    assert!(
        limpieza
            .burbujas
            .iter()
            .all(|b| b.cobertura_conservada == 10.0)
    );
    assert_eq!(solo_unitig(&limpieza), orientado(&genoma));
}

#[test]
fn recompactar_da_lo_mismo_que_compactar_de_nuevo() {
    let mut azar = Azar::nuevo(41);
    let mut revisados = 0;
    for _ in 0..200 {
        let k = 5 + 2 * azar.entero(4);
        let conteos = kmeros_con_errores(&mut azar, k);
        let mut kmeros = tabla_desde(&conteos, k);
        let grafo = compactar(&kmeros);
        let mut quitados: Vec<usize> = buscar_puntas(&grafo)
            .iter()
            .map(|punta| punta.unitig)
            .chain(
                buscar_burbujas(&grafo)
                    .iter()
                    .map(|burbuja| burbuja.unitig_quitado),
            )
            .collect();
        if quitados.is_empty() && grafo.unitigs.len() > 1 {
            quitados.push(azar.entero(grafo.unitigs.len()));
        }
        quitados.sort_unstable();
        quitados.dedup();
        for &numero in &quitados {
            quitar_unitig(&mut kmeros, &grafo.unitigs[numero]);
        }
        let nuevo = recompactar(grafo, &kmeros, &quitados);
        let desde_cero = compactar(&kmeros);
        assert_eq!(particion_de(&nuevo), particion_de(&desde_cero));
        let resumen = |grafo: &GrafoDeUnitigs| {
            let mut lista: Vec<(Vec<u8>, u64, usize, usize)> = grafo
                .unitigs
                .iter()
                .map(|unitig| {
                    (
                        unitig.secuencia.clone(),
                        unitig.cobertura.to_bits(),
                        unitig.al_final.len(),
                        unitig.al_inicio.len(),
                    )
                })
                .collect();
            lista.sort();
            lista
        };
        assert_eq!(resumen(&nuevo), resumen(&desde_cero));
        revisados += quitados.len();
    }
    assert!(revisados > 200);
}

fn cruce_armado(contexto: usize, k: usize) -> (Kmeros, Vec<Vec<u8>>) {
    let mut azar = Azar::nuevo(77);
    loop {
        let izquierdas = [genoma_al_azar(40, &mut azar), genoma_al_azar(40, &mut azar)];
        let derechas = [genoma_al_azar(40, &mut azar), genoma_al_azar(40, &mut azar)];
        let repeticion = genoma_al_azar(40, &mut azar);
        if izquierdas[0][39] == izquierdas[1][39] || derechas[0][0] == derechas[1][0] {
            continue;
        }
        let mut variante = repeticion.clone();
        variante[10] = if repeticion[10] == b'A' { b'C' } else { b'A' };
        let lecturas = vec![
            [&izquierdas[0][..], &repeticion, &derechas[0]].concat(),
            [&izquierdas[1][..], &repeticion, &derechas[1]].concat(),
            [&izquierdas[contexto][..], &variante, &derechas[contexto]].concat(),
        ];
        let mut kmeros = Kmeros::nuevo(k);
        agregar_texto(&mut kmeros, &lecturas[0], 5);
        agregar_texto(&mut kmeros, &lecturas[1], 5);
        agregar_texto(&mut kmeros, &lecturas[2], 20);
        let grafo = compactar(&kmeros);
        let rama = orientado(&repeticion[..21]);
        let con_rama = grafo.unitigs.iter().find(|unitig| unitig.secuencia == rama);
        if con_rama.is_some_and(|unitig| unitig.al_inicio.len() + unitig.al_final.len() == 3) {
            return (kmeros, lecturas);
        }
    }
}

#[test]
fn no_revienta_una_rama_que_tambien_une_otro_contexto() {
    let k = 11;
    for contexto in 0..2 {
        let (mut kmeros, lecturas) = cruce_armado(contexto, k);
        let limpieza = limpiar(&mut kmeros);
        assert_eq!(limpieza.burbujas.len(), 0);
        for lectura in &lecturas[..2] {
            recorrer_kmeros(lectura, k, |_, kmero| {
                assert!(kmeros.contiene(kmero.canonico()));
            });
        }
    }
}

#[test]
fn rechaza_demasiados_kmeros_distintos() {
    assert!(revisar_kmeros_distintos(MAXIMO_DE_KMEROS_DISTINTOS).is_ok());
    let error = revisar_kmeros_distintos(MAXIMO_DE_KMEROS_DISTINTOS + 1).unwrap_err();
    assert!(error.contains("k-meros distintos"));
}

#[test]
fn rechaza_una_lectura_larguisima() {
    let ruta = ruta_temporal("fastq");
    let largo = LARGO_MAXIMO_DE_LECTURA + 1;
    let mut texto = Vec::with_capacity(2 * largo + 20);
    texto.extend_from_slice(b"@larga\n");
    texto.extend(std::iter::repeat_n(b'A', largo));
    texto.extend_from_slice(b"\n+\n");
    texto.extend(std::iter::repeat_n(b'I', largo));
    texto.push(b'\n');
    std::fs::write(&ruta, texto).unwrap();
    let resultado = ensamblar(&ruta, 31, 0);
    std::fs::remove_file(&ruta).ok();
    assert!(resultado.err().unwrap().contains("lectura de más de"));
}

#[test]
fn las_lecturas_largas_se_reparten_en_bloques_chicos() {
    let k = 21;
    let mut azar = Azar::nuevo(5);
    let genoma = genoma_sin_repeticiones(3_000, k, &mut azar);
    let lecturas = vec![genoma.clone(), genoma.clone(), genoma];
    let ruta = escribir_lecturas(&lecturas);
    let conteo = contar_kmeros(&ruta, k).unwrap();
    std::fs::remove_file(&ruta).ok();
    assert_eq!(conteo.lecturas, 3);
    assert_eq!(conteo.kmeros_leidos, 3 * (3_000 - k + 1));
    assert_eq!(conteo.kmeros.cantidad(), 3_000 - k + 1);
}
