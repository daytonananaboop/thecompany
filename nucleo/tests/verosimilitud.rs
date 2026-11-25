use nucleo::azar::Azar;
use nucleo::filogenia::{self, Arbol};
use nucleo::verosimilitud::{
    self, Calculo, Eleccion, Modelo, Referencia, Sustitucion, comprimir_en_patrones,
};

fn indices(cantidad: usize) -> Vec<usize> {
    (0..cantidad).collect()
}

fn nombres(cantidad: usize) -> Vec<String> {
    (0..cantidad).map(|i| format!("e{i}")).collect()
}

fn arbol_con_aristas(hojas: usize, aristas: &[(usize, usize, f64)]) -> Arbol {
    let mut arbol = Arbol::default();
    let internos = aristas
        .iter()
        .map(|&(a, b, _)| a.max(b) + 1)
        .max()
        .unwrap_or(hojas);
    for nodo in 0..internos {
        arbol.agregar_nodo(&format!("n{nodo}"));
    }
    for &(a, b, largo) in aristas {
        arbol.conectar(a, b, largo);
    }
    arbol
}

fn arbol_al_azar(hojas: usize, azar: &mut Azar) -> Arbol {
    let mut arbol = filogenia::arbol_al_azar(&nombres(hojas), azar.siguiente()).unwrap();
    for (a, b, _) in arbol.aristas() {
        let largo = 0.02 + 0.4 * azar.decimal();
        verosimilitud::fijar_longitud(&mut arbol, a, b, largo);
    }
    arbol
}

fn modelo_al_azar(azar: &mut Azar) -> Modelo {
    let sustitucion = [Sustitucion::Jc69, Sustitucion::K80, Sustitucion::Hky85][azar.entero(3)];
    let mut frecuencias = [0.0; 4];
    for frecuencia in frecuencias.iter_mut() {
        *frecuencia = 0.1 + azar.decimal();
    }
    let suma: f64 = frecuencias.iter().sum();
    let mut modelo = Modelo::nuevo(
        sustitucion,
        azar.decimal() < 0.5,
        frecuencias.map(|f| f / suma),
    );
    if sustitucion != Sustitucion::Jc69 {
        modelo.kappa = 0.5 + 6.0 * azar.decimal();
    }
    if modelo.alfa.is_some() {
        modelo.alfa = Some(0.2 + 3.0 * azar.decimal());
    }
    modelo
}

fn elegir(azar: &mut Azar, pesos: &[f64; 4]) -> usize {
    let mut objetivo = azar.decimal() * pesos.iter().sum::<f64>();
    for (indice, peso) in pesos.iter().enumerate() {
        if objetivo < *peso {
            return indice;
        }
        objetivo -= peso;
    }
    3
}

fn simular(
    arbol: &Arbol,
    hojas: usize,
    modelo: &Modelo,
    largo: usize,
    azar: &mut Azar,
) -> Vec<Vec<u8>> {
    let espectro = verosimilitud::espectro(modelo);
    let tasas = modelo.tasas();
    let mut secuencias = vec![vec![b'A'; largo]; hojas];
    let raiz = hojas;
    for sitio in indices(largo) {
        let tasa = tasas[azar.entero(tasas.len())];
        let mut estados = vec![0usize; arbol.nombres.len()];
        estados[raiz] = elegir(azar, &modelo.frecuencias);
        let mut pila = vec![(raiz, usize::MAX)];
        while let Some((nodo, padre)) = pila.pop() {
            for &(hijo, largo_rama) in &arbol.vecinos[nodo] {
                if hijo == padre {
                    continue;
                }
                let matriz = verosimilitud::matriz_de_transicion(&espectro, tasa * largo_rama);
                estados[hijo] = elegir(azar, &matriz[estados[nodo]]);
                pila.push((hijo, nodo));
            }
        }
        for (hoja, secuencia) in secuencias.iter_mut().enumerate() {
            secuencia[sitio] = b"ACGT"[estados[hoja]];
        }
    }
    secuencias
}

fn log_verosimilitud(secuencias: &[Vec<u8>], arbol: &Arbol, modelo: Modelo) -> f64 {
    let datos = comprimir_en_patrones(secuencias, &nombres(secuencias.len()));
    let mut calculo = Calculo::nuevo(&datos, modelo, arbol.nombres.len());
    calculo.log_verosimilitud(arbol)
}

fn fuerza_bruta(secuencias: &[Vec<u8>], arbol: &Arbol, modelo: &Modelo) -> f64 {
    let hojas = secuencias.len();
    let nodos = arbol.nombres.len();
    let internos: Vec<usize> = (hojas..nodos).collect();
    let espectro = verosimilitud::espectro(modelo);
    let tasas = modelo.tasas();
    let aristas = arbol.aristas();
    let mut total = 0.0;
    for sitio in indices(secuencias[0].len()) {
        let mut del_sitio = 0.0;
        for &tasa in &tasas {
            let matrices: Vec<[[f64; 4]; 4]> = aristas
                .iter()
                .map(|&(_, _, largo)| verosimilitud::matriz_de_transicion(&espectro, tasa * largo))
                .collect();
            for asignacion in 0..4usize.pow(internos.len() as u32) {
                let mut estados = vec![0usize; nodos];
                let mut resto = asignacion;
                for &interno in &internos {
                    estados[interno] = resto % 4;
                    resto /= 4;
                }
                let raiz = internos[0];
                let mut producto = modelo.frecuencias[estados[raiz]];
                let mut suma_hojas = 1.0;
                for (indice, &(a, b, _)) in aristas.iter().enumerate() {
                    let (interno, otro) = if a >= hojas { (a, b) } else { (b, a) };
                    if otro >= hojas {
                        let (padre, hijo) = orientar(arbol, raiz, a, b);
                        producto *= matrices[indice][estados[padre]][estados[hijo]];
                        continue;
                    }
                    let letra = secuencias[otro][sitio];
                    let estado = verosimilitud::estado_de_letra(letra);
                    let vector = verosimilitud::vector_de_estado(estado);
                    let mut suma = 0.0;
                    for (destino, valor) in vector.iter().enumerate() {
                        suma += matrices[indice][estados[interno]][destino] * valor;
                    }
                    suma_hojas *= suma;
                }
                del_sitio += producto * suma_hojas / tasas.len() as f64;
            }
        }
        total += del_sitio.ln();
    }
    total
}

fn orientar(arbol: &Arbol, raiz: usize, a: usize, b: usize) -> (usize, usize) {
    let distancias = filogenia::distancias_desde(&sin_largos(arbol), raiz);
    if distancias[a] < distancias[b] {
        (a, b)
    } else {
        (b, a)
    }
}

fn sin_largos(arbol: &Arbol) -> Arbol {
    let mut copia = arbol.clone();
    for vecinos in copia.vecinos.iter_mut() {
        for (_, largo) in vecinos.iter_mut() {
            *largo = 1.0;
        }
    }
    copia
}

fn todas_las_topologias(hojas: usize) -> Vec<Arbol> {
    let mut base = Arbol::default();
    for nombre in nombres(hojas) {
        base.agregar_nodo(&nombre);
    }
    let centro = base.agregar_nodo("c");
    for hoja in indices(3) {
        base.conectar(centro, hoja, 0.1);
    }
    let mut arboles = vec![base];
    for hoja in 3..hojas {
        let mut siguientes = Vec::new();
        for arbol in &arboles {
            for (a, b, _) in arbol.aristas() {
                let mut nuevo = arbol.clone();
                let medio = nuevo.agregar_nodo("m");
                nuevo.desconectar(a, b);
                nuevo.conectar(a, medio, 0.1);
                nuevo.conectar(medio, b, 0.1);
                nuevo.conectar(medio, hoja, 0.1);
                siguientes.push(nuevo);
            }
        }
        arboles = siguientes;
    }
    arboles
}

#[test]
fn la_matriz_de_transicion_es_una_probabilidad() {
    let mut azar = Azar::nuevo(3);
    for _ in 0..50 {
        let modelo = modelo_al_azar(&mut azar);
        let espectro = verosimilitud::espectro(&modelo);
        let identidad = verosimilitud::matriz_de_transicion(&espectro, 0.0);
        for i in indices(4) {
            for j in indices(4) {
                let esperado = if i == j { 1.0 } else { 0.0 };
                assert!((identidad[i][j] - esperado).abs() < 1e-12);
            }
        }
        let t = 0.01 + azar.decimal();
        let matriz = verosimilitud::matriz_de_transicion(&espectro, t);
        for i in indices(4) {
            let suma: f64 = matriz[i].iter().sum();
            assert!((suma - 1.0).abs() < 1e-12);
            for j in indices(4) {
                assert!(matriz[i][j] >= 0.0);
                let ida = modelo.frecuencias[i] * matriz[i][j];
                let vuelta = modelo.frecuencias[j] * matriz[j][i];
                assert!((ida - vuelta).abs() < 1e-12);
            }
        }
        let lejos = verosimilitud::matriz_de_transicion(&espectro, 200.0);
        for fila in lejos {
            for j in indices(4) {
                assert!((fila[j] - modelo.frecuencias[j]).abs() < 1e-9);
            }
        }
    }
}

fn exponencial_por_serie(tasa: &[[f64; 4]; 4], tiempo: f64) -> [[f64; 4]; 4] {
    let mut escalada = [[0.0; 4]; 4];
    let cuadrados = 10;
    let factor = tiempo / 2f64.powi(cuadrados);
    for i in indices(4) {
        for j in indices(4) {
            escalada[i][j] = tasa[i][j] * factor;
        }
    }
    let mut resultado = [[0.0; 4]; 4];
    let mut termino = [[0.0; 4]; 4];
    for i in indices(4) {
        resultado[i][i] = 1.0;
        termino[i][i] = 1.0;
    }
    for k in 1..20 {
        termino = producto(&termino, &escalada);
        for fila in termino.iter_mut() {
            for valor in fila.iter_mut() {
                *valor /= k as f64;
            }
        }
        for i in indices(4) {
            for j in indices(4) {
                resultado[i][j] += termino[i][j];
            }
        }
    }
    for _ in 0..cuadrados {
        resultado = producto(&resultado, &resultado);
    }
    resultado
}

fn producto(a: &[[f64; 4]; 4], b: &[[f64; 4]; 4]) -> [[f64; 4]; 4] {
    let mut c = [[0.0; 4]; 4];
    for i in indices(4) {
        for j in indices(4) {
            for k in indices(4) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

#[test]
fn la_forma_cerrada_de_hky_coincide_con_la_exponencial_de_q() {
    let mut azar = Azar::nuevo(11);
    for _ in 0..40 {
        let modelo = modelo_al_azar(&mut azar);
        let pi = modelo.frecuencias;
        let mut tasa = [[0.0; 4]; 4];
        for i in indices(4) {
            for j in indices(4) {
                if i != j {
                    let transicion = (i + j) % 2 == 0;
                    tasa[i][j] = pi[j] * if transicion { modelo.kappa } else { 1.0 };
                }
            }
        }
        let mut media = 0.0;
        for i in indices(4) {
            let salida: f64 = tasa[i].iter().sum();
            tasa[i][i] = -salida;
            media += pi[i] * salida;
        }
        for fila in tasa.iter_mut() {
            for valor in fila.iter_mut() {
                *valor /= media;
            }
        }
        let tiempo = 0.05 + 2.0 * azar.decimal();
        let esperada = exponencial_por_serie(&tasa, tiempo);
        let obtenida =
            verosimilitud::matriz_de_transicion(&verosimilitud::espectro(&modelo), tiempo);
        for i in indices(4) {
            for j in indices(4) {
                assert!(
                    (esperada[i][j] - obtenida[i][j]).abs() < 1e-10,
                    "{i},{j}: {} contra {}",
                    esperada[i][j],
                    obtenida[i][j]
                );
            }
        }
    }
}

#[test]
fn las_tasas_gamma_de_yang_tienen_media_uno() {
    let tasas = verosimilitud::tasas_gamma(0.5, 4);
    let esperadas = [0.033_388, 0.251_916, 0.820_268, 2.894_428];
    for (tasa, esperada) in tasas.iter().zip(esperadas) {
        assert!((tasa - esperada).abs() < 1e-5, "{tasa} contra {esperada}");
    }
    for alfa in [0.02, 0.1, 0.7, 1.0, 3.0, 20.0, 100.0] {
        let tasas = verosimilitud::tasas_gamma(alfa, 4);
        let media: f64 = tasas.iter().sum::<f64>() / 4.0;
        assert!((media - 1.0).abs() < 1e-9, "α = {alfa}: media {media}");
        assert!(tasas.windows(2).all(|par| par[0] <= par[1]));
    }
    let casi_iguales = verosimilitud::tasas_gamma(100.0, 4);
    assert!(casi_iguales.iter().all(|tasa| (tasa - 1.0).abs() < 0.15));
}

#[test]
fn la_poda_coincide_con_sumar_todas_las_asignaciones() {
    let mut azar = Azar::nuevo(21);
    for caso in 0..60 {
        let hojas = 3 + caso % 4;
        let arbol = arbol_al_azar(hojas, &mut azar);
        let modelo = modelo_al_azar(&mut azar);
        let mut secuencias = simular(&arbol, hojas, &modelo, 12, &mut azar);
        secuencias[0][3] = b'-';
        secuencias[hojas - 1][5] = b'N';
        let rapida = log_verosimilitud(&secuencias, &arbol, modelo);
        let lenta = fuerza_bruta(&secuencias, &arbol, &modelo);
        assert!(
            (rapida - lenta).abs() < 1e-8 * lenta.abs().max(1.0),
            "caso {caso}: {rapida} contra {lenta}"
        );
    }
}

#[test]
fn agrupar_columnas_iguales_no_cambia_la_verosimilitud() {
    let mut azar = Azar::nuevo(5);
    let arbol = arbol_al_azar(7, &mut azar);
    let modelo = modelo_al_azar(&mut azar);
    let secuencias = simular(&arbol, 7, &modelo, 400, &mut azar);
    let datos = comprimir_en_patrones(&secuencias, &nombres(7));
    assert!(datos.pesos.len() < 400);
    assert_eq!(datos.pesos.iter().sum::<f64>(), 400.0);
    let junta = log_verosimilitud(&secuencias, &arbol, modelo);
    let mut por_columna = 0.0;
    for sitio in indices(400) {
        let columna: Vec<Vec<u8>> = secuencias.iter().map(|s| vec![s[sitio]]).collect();
        por_columna += log_verosimilitud(&columna, &arbol, modelo);
    }
    assert!((junta - por_columna).abs() < 1e-8);
}

fn poda_en_logaritmos(secuencias: &[Vec<u8>], arbol: &Arbol, rama_minima: f64) -> f64 {
    let hojas = secuencias.len();
    let modelo = Modelo::nuevo(Sustitucion::Jc69, false, [0.25; 4]);
    let espectro = verosimilitud::espectro(&modelo);
    let mut total = 0.0;
    for sitio in indices(secuencias[0].len()) {
        let raiz = hojas;
        let vector = log_del_nodo(
            secuencias,
            arbol,
            &espectro,
            sitio,
            raiz,
            usize::MAX,
            rama_minima,
        );
        let partes: Vec<f64> = vector.iter().map(|v| 0.25f64.ln() + v).collect();
        total += log_suma(&partes);
    }
    total
}

fn log_suma(partes: &[f64]) -> f64 {
    let mayor = partes.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    mayor + partes.iter().map(|p| (p - mayor).exp()).sum::<f64>().ln()
}

fn log_del_nodo(
    secuencias: &[Vec<u8>],
    arbol: &Arbol,
    espectro: &verosimilitud::Espectro,
    sitio: usize,
    nodo: usize,
    padre: usize,
    rama_minima: f64,
) -> [f64; 4] {
    if nodo < secuencias.len() {
        let estado = verosimilitud::estado_de_letra(secuencias[nodo][sitio]);
        return verosimilitud::vector_de_estado(estado).map(|v| v.ln());
    }
    let mut resultado = [0.0; 4];
    for &(hijo, largo) in &arbol.vecinos[nodo] {
        if hijo == padre {
            continue;
        }
        let debajo = log_del_nodo(secuencias, arbol, espectro, sitio, hijo, nodo, rama_minima);
        let matriz = verosimilitud::matriz_de_transicion(espectro, largo.max(rama_minima));
        for x in indices(4) {
            let partes: Vec<f64> = (0..4).map(|y| matriz[x][y].ln() + debajo[y]).collect();
            resultado[x] += log_suma(&partes);
        }
    }
    resultado
}

#[test]
fn muchas_hojas_no_se_desbordan_gracias_al_reescalado() {
    let mut azar = Azar::nuevo(8);
    let mut arbol = arbol_al_azar(60, &mut azar);
    for (a, b, _) in arbol.aristas() {
        verosimilitud::fijar_longitud(&mut arbol, a, b, 1e-4);
    }
    let secuencias: Vec<Vec<u8>> = (0..60)
        .map(|_| (0..20).map(|_| b"ACGT"[azar.entero(4)]).collect())
        .collect();
    let modelo = Modelo::nuevo(Sustitucion::Jc69, false, [0.25; 4]);
    let valor = log_verosimilitud(&secuencias, &arbol, modelo);
    let esperado = poda_en_logaritmos(&secuencias, &arbol, 0.0);
    assert!(esperado < -20.0 * 177.0, "{esperado}");
    assert!(
        (valor - esperado).abs() < 1e-6 * esperado.abs(),
        "{valor} contra {esperado}"
    );
}

#[test]
fn con_tres_especies_la_suma_de_ramas_es_la_distancia_de_jukes_cantor() {
    let primera = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let mut segunda = primera.clone();
    for posicion in [0, 5, 10, 15, 20, 25, 30, 35] {
        segunda[posicion] = if primera[posicion] == b'A' {
            b'C'
        } else {
            b'A'
        };
    }
    let secuencias = vec![primera.clone(), segunda, primera];
    let datos = comprimir_en_patrones(&secuencias, &nombres(3));
    let arbol = arbol_con_aristas(3, &[(3, 0, 0.1), (3, 1, 0.1), (3, 2, 0.1)]);
    let ajuste = verosimilitud::ajustar_modelo(
        &datos,
        &arbol,
        Modelo::nuevo(Sustitucion::Jc69, false, [0.25; 4]),
    );
    let p: f64 = 8.0 / 40.0;
    let esperada = -0.75 * (1.0 - 4.0 * p / 3.0).ln();
    let a = ajuste.arbol.longitud_entre(3, 0).unwrap();
    let b = ajuste.arbol.longitud_entre(3, 1).unwrap();
    let c = ajuste.arbol.longitud_entre(3, 2).unwrap();
    assert!(
        (a + b - esperada).abs() < 1e-4,
        "{a} + {b} contra {esperada}"
    );
    assert!(a < 1e-4 && c < 1e-4);
    let maxima = 32.0 * (0.25 * (1.0 - p)).ln() + 8.0 * (0.25 * p / 3.0).ln();
    assert!((ajuste.log_verosimilitud - maxima).abs() < 1e-4);
}

#[test]
fn secuencias_iguales_dan_ramas_casi_nulas() {
    let secuencia = b"ACGTTGCAACGTAGCT".to_vec();
    let secuencias = vec![secuencia.clone(); 5];
    let resultado = verosimilitud::maxima_verosimilitud(
        &secuencias,
        &nombres(5),
        Eleccion::Fija(Sustitucion::Jc69, false),
    )
    .unwrap();
    for (_, _, largo) in resultado.arbol.aristas() {
        assert!(largo < 1e-4);
    }
    let esperado = 16.0 * 0.25f64.ln();
    assert!((resultado.log_verosimilitud - esperado).abs() < 1e-3);
}

#[test]
fn cada_rama_queda_en_su_maximo() {
    let mut azar = Azar::nuevo(31);
    for _ in 0..10 {
        let arbol = arbol_al_azar(6, &mut azar);
        let modelo = modelo_al_azar(&mut azar);
        let secuencias = simular(&arbol, 6, &modelo, 300, &mut azar);
        let datos = comprimir_en_patrones(&secuencias, &nombres(6));
        let ajuste = verosimilitud::ajustar_modelo(&datos, &arbol, modelo);
        let base = log_verosimilitud(&secuencias, &ajuste.arbol, ajuste.modelo);
        assert!((base - ajuste.log_verosimilitud).abs() < 1e-6);
        for (a, b, largo) in ajuste.arbol.aristas() {
            for factor in [0.9, 1.1] {
                let mut movido = ajuste.arbol.clone();
                verosimilitud::fijar_longitud(&mut movido, a, b, (largo * factor).max(1e-8));
                let valor = log_verosimilitud(&secuencias, &movido, ajuste.modelo);
                assert!(valor <= base + 2e-3, "{valor} > {base}");
            }
        }
        if let Some(alfa) = ajuste.modelo.alfa {
            for factor in [0.8, 1.25] {
                let otro = Modelo {
                    alfa: Some(alfa * factor),
                    ..ajuste.modelo
                };
                assert!(log_verosimilitud(&secuencias, &ajuste.arbol, otro) <= base + 2e-3);
            }
        }
    }
}

#[test]
fn los_intercambios_encuentran_el_mejor_arbol_de_cinco_especies() {
    let mut azar = Azar::nuevo(41);
    let topologias = todas_las_topologias(5);
    assert_eq!(topologias.len(), 15);
    let mut aciertos = 0;
    let casos = 20;
    for _ in 0..casos {
        let arbol = arbol_al_azar(5, &mut azar);
        let modelo = Modelo::nuevo(Sustitucion::Jc69, false, [0.25; 4]);
        let secuencias = simular(&arbol, 5, &modelo, 150, &mut azar);
        let datos = comprimir_en_patrones(&secuencias, &nombres(5));
        let ajustes: Vec<f64> = topologias
            .iter()
            .map(|topologia| {
                verosimilitud::ajustar_modelo(&datos, topologia, modelo).log_verosimilitud
            })
            .collect();
        let mejor = (0..ajustes.len())
            .max_by(|&a, &b| ajustes[a].total_cmp(&ajustes[b]))
            .unwrap();
        let inicio = arbol_al_azar(5, &mut azar);
        let busqueda = verosimilitud::buscar_con_intercambios(&datos, &inicio, modelo);
        let final_ = verosimilitud::ajustar_modelo(&datos, &busqueda.arbol, modelo);
        assert!(final_.log_verosimilitud <= ajustes[mejor] + 5e-2);
        assert!(
            busqueda
                .rondas
                .windows(2)
                .all(|par| par[1].log_verosimilitud > par[0].log_verosimilitud)
        );
        if verosimilitud::distancia_robinson_foulds(&busqueda.arbol, &topologias[mejor], 5) == 0 {
            aciertos += 1;
        }
    }
    assert!(aciertos >= casos - 2, "{aciertos} de {casos}");
}

#[test]
fn recupera_el_arbol_verdadero_con_senal_clara() {
    let mut azar = Azar::nuevo(51);
    let verdadero = arbol_al_azar(8, &mut azar);
    let mut modelo = Modelo::nuevo(Sustitucion::Hky85, true, [0.3, 0.2, 0.2, 0.3]);
    modelo.kappa = 4.0;
    modelo.alfa = Some(0.8);
    let secuencias = simular(&verdadero, 8, &modelo, 3000, &mut azar);
    let resultado =
        verosimilitud::maxima_verosimilitud(&secuencias, &nombres(8), Eleccion::MejorPorAic)
            .unwrap();
    assert_eq!(
        verosimilitud::distancia_robinson_foulds(&resultado.arbol, &verdadero, 8),
        0
    );
    assert_eq!(resultado.modelo.sustitucion, Sustitucion::Hky85);
    assert!(resultado.modelo.alfa.is_some());
    assert!((resultado.modelo.kappa - 4.0).abs() < 1.0);
    assert!(resultado.log_verosimilitud >= resultado.log_verosimilitud_de_union - 1e-6);
    let mejor = &resultado.comparacion[resultado.elegido];
    assert!(
        resultado
            .comparacion
            .iter()
            .all(|otro| otro.aic >= mejor.aic)
    );
}

#[test]
fn particiones_y_robinson_foulds() {
    let primero = arbol_con_aristas(
        5,
        &[
            (5, 0, 0.1),
            (5, 1, 0.1),
            (5, 6, 0.1),
            (6, 2, 0.1),
            (6, 7, 0.1),
            (7, 3, 0.1),
            (7, 4, 0.1),
        ],
    );
    let segundo = arbol_con_aristas(
        5,
        &[
            (5, 0, 0.1),
            (5, 2, 0.1),
            (5, 6, 0.1),
            (6, 1, 0.1),
            (6, 7, 0.1),
            (7, 3, 0.1),
            (7, 4, 0.1),
        ],
    );
    let mut lados: Vec<u64> = verosimilitud::particiones(&primero, 5)
        .iter()
        .map(|p| p.hojas)
        .collect();
    lados.sort();
    assert_eq!(lados, vec![0b11000, 0b11100]);
    assert_eq!(
        verosimilitud::distancia_robinson_foulds(&primero, &primero, 5),
        0
    );
    assert_eq!(
        verosimilitud::distancia_robinson_foulds(&primero, &segundo, 5),
        2
    );
}

#[test]
fn las_distancias_de_los_patrones_son_las_del_capitulo_siete() {
    let mut azar = Azar::nuevo(61);
    let arbol = arbol_al_azar(6, &mut azar);
    let modelo = Modelo::nuevo(Sustitucion::K80, false, [0.25; 4]);
    let mut secuencias = simular(&arbol, 6, &modelo, 200, &mut azar);
    secuencias[2][7] = b'-';
    let datos = comprimir_en_patrones(&secuencias, &nombres(6));
    let rapidas = verosimilitud::distancias_jukes_cantor(&datos);
    let clasicas = filogenia::matriz_jukes_cantor(&secuencias).unwrap();
    for i in indices(6) {
        for j in indices(6) {
            assert!((rapidas[i][j] - clasicas[i][j]).abs() < 1e-12);
        }
    }
}

#[test]
fn remuestrear_conserva_el_total_de_sitios() {
    let mut azar = Azar::nuevo(71);
    let arbol = arbol_al_azar(5, &mut azar);
    let modelo = Modelo::nuevo(Sustitucion::Jc69, false, [0.25; 4]);
    let secuencias = simular(&arbol, 5, &modelo, 120, &mut azar);
    let datos = comprimir_en_patrones(&secuencias, &nombres(5));
    let columnas = verosimilitud::columnas_de_replica(120, 3, 9);
    assert_eq!(columnas, verosimilitud::columnas_de_replica(120, 3, 9));
    assert_ne!(columnas, verosimilitud::columnas_de_replica(120, 4, 9));
    let muestra = verosimilitud::remuestrear(&datos, &columnas);
    assert_eq!(muestra.pesos.iter().sum::<f64>(), 120.0);
    assert!(muestra.pesos.iter().all(|&peso| peso > 0.0));
    for (sitio, &columna) in columnas.iter().enumerate() {
        let patron = muestra.patron_de_sitio[sitio];
        for hoja in indices(5) {
            assert_eq!(
                muestra.estados[hoja][patron],
                verosimilitud::estado_de_letra(secuencias[hoja][columna])
            );
        }
    }
}

#[test]
fn el_bootstrap_apoya_del_todo_los_grupos_claros() {
    let mut azar = Azar::nuevo(81);
    let verdadero = arbol_con_aristas(
        6,
        &[
            (6, 0, 0.05),
            (6, 1, 0.05),
            (6, 7, 0.3),
            (7, 2, 0.05),
            (7, 8, 0.3),
            (8, 9, 0.3),
            (9, 3, 0.05),
            (9, 4, 0.05),
            (8, 5, 0.05),
        ],
    );
    let modelo = Modelo::nuevo(Sustitucion::Jc69, false, [0.25; 4]);
    let secuencias = simular(&verdadero, 6, &modelo, 600, &mut azar);
    for referencia in [Referencia::MaximaVerosimilitud, Referencia::UnionDeVecinos] {
        let resultado =
            verosimilitud::bootstrap(&secuencias, &nombres(6), 50, referencia, 1).unwrap();
        assert_eq!(resultado.particiones_por_replica.len(), 50);
        assert_eq!(resultado.apoyos.len(), 3);
        for apoyo in &resultado.apoyos {
            assert!(apoyo.porcentaje >= 95.0, "{apoyo:?}");
        }
        let newick = verosimilitud::newick_con_apoyos(&resultado.arbol, 6, &resultado.apoyos);
        assert!(newick.contains(")100:") || newick.contains(")98:") || newick.contains(")96:"));
        let otra = verosimilitud::bootstrap(&secuencias, &nombres(6), 50, referencia, 1).unwrap();
        assert_eq!(otra, resultado);
    }
}

#[test]
fn el_bootstrap_sin_senal_no_apoya_nada() {
    let mut azar = Azar::nuevo(91);
    let secuencias: Vec<Vec<u8>> = (0..7)
        .map(|_| (0..300).map(|_| b"ACGT"[azar.entero(4)]).collect())
        .collect();
    let resultado =
        verosimilitud::bootstrap(&secuencias, &nombres(7), 60, Referencia::UnionDeVecinos, 2)
            .unwrap();
    assert!(resultado.apoyos.iter().all(|apoyo| apoyo.porcentaje < 70.0));
}

#[test]
fn el_bootstrap_de_secuencias_identicas_no_apoya_ramas_inventadas() {
    let secuencias = vec![b"ACGTACGTACGT".to_vec(); 6];
    for referencia in [Referencia::UnionDeVecinos, Referencia::MaximaVerosimilitud] {
        let resultado =
            verosimilitud::bootstrap(&secuencias, &nombres(6), 20, referencia, 3).unwrap();
        assert_eq!(resultado.apoyos.len(), 3);
        assert!(resultado.apoyos.iter().all(|apoyo| apoyo.veces == 0));
    }
}

#[test]
fn las_ramas_de_longitud_cero_no_son_particiones_resueltas() {
    let mut arbol = arbol_con_aristas(
        4,
        &[
            (0, 4, 0.1),
            (1, 4, 0.1),
            (4, 5, 0.0),
            (2, 5, 0.1),
            (3, 5, 0.1),
        ],
    );
    assert_eq!(verosimilitud::particiones(&arbol, 4).len(), 1);
    assert!(verosimilitud::particiones_resueltas(&arbol, 4).is_empty());
    verosimilitud::fijar_longitud(&mut arbol, 4, 5, 0.2);
    assert_eq!(verosimilitud::particiones_resueltas(&arbol, 4).len(), 1);
}

#[test]
fn rechaza_entradas_que_no_sirven() {
    let ok = vec![b"ACGT".to_vec(), b"ACGA".to_vec(), b"ACTT".to_vec()];
    let dos = vec![b"ACGT".to_vec(), b"ACGA".to_vec()];
    let disparejas = vec![b"ACGT".to_vec(), b"ACG".to_vec(), b"ACTT".to_vec()];
    let raras = vec![b"ACGT".to_vec(), b"ACGX".to_vec(), b"ACTT".to_vec()];
    let vacias = vec![Vec::new(), Vec::new(), Vec::new()];
    let elegir = Eleccion::MejorPorAic;
    assert!(verosimilitud::maxima_verosimilitud(&dos, &nombres(2), elegir).is_err());
    assert!(verosimilitud::maxima_verosimilitud(&disparejas, &nombres(3), elegir).is_err());
    assert!(verosimilitud::maxima_verosimilitud(&raras, &nombres(3), elegir).is_err());
    assert!(verosimilitud::maxima_verosimilitud(&vacias, &nombres(3), elegir).is_err());
    assert!(verosimilitud::maxima_verosimilitud(&ok, &nombres(2), elegir).is_err());
    let repetidos = vec!["a".to_string(), "a".to_string(), "b".to_string()];
    assert!(verosimilitud::maxima_verosimilitud(&ok, &repetidos, elegir).is_err());
    let largas = vec![vec![b'A'; 5001]; 3];
    assert!(verosimilitud::maxima_verosimilitud(&largas, &nombres(3), elegir).is_err());
    let muchas = vec![b"ACGT".to_vec(); 61];
    assert!(verosimilitud::maxima_verosimilitud(&muchas, &nombres(61), elegir).is_err());
    assert!(verosimilitud::maxima_verosimilitud(&ok, &nombres(3), elegir).is_ok());
    let cuatro = vec![
        b"ACGT".to_vec(),
        b"ACGA".to_vec(),
        b"ACTT".to_vec(),
        b"TCTT".to_vec(),
    ];
    let referencia = Referencia::UnionDeVecinos;
    assert!(verosimilitud::bootstrap(&ok, &nombres(3), 10, referencia, 1).is_err());
    assert!(verosimilitud::bootstrap(&cuatro, &nombres(4), 0, referencia, 1).is_err());
    assert!(verosimilitud::bootstrap(&cuatro, &nombres(4), 201, referencia, 1).is_err());
    assert!(verosimilitud::bootstrap(&cuatro, &nombres(4), 10, referencia, 1).is_ok());
    let opuestas = vec![b"AAAA".to_vec(), b"CCCC".to_vec(), b"GGGG".to_vec()];
    let saturadas = verosimilitud::maxima_verosimilitud(&opuestas, &nombres(3), elegir).unwrap();
    assert!(saturadas.log_verosimilitud.is_finite());
    let sin_comun = vec![b"AC--".to_vec(), b"--GT".to_vec(), b"ACGT".to_vec()];
    assert!(verosimilitud::maxima_verosimilitud(&sin_comun, &nombres(3), elegir).is_ok());
}
