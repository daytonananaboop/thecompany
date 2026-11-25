use crate::azar::Azar;
use rayon::prelude::*;

const YA_ES_CENTRO: f64 = -1.0;
const INTENTOS_POR_K: u64 = 3;
const ITERACIONES_DEL_CODO: usize = 100;
const VUELTAS_DE_POTENCIA: usize = 1000;
const CAMBIO_DESPRECIABLE: f64 = 1e-24;
const PARTE_DESPRECIABLE: f64 = 1e-12;
const SEMILLA_DE_POTENCIA: u64 = 1;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Medida {
    Euclidea,
    Pearson,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Inicio {
    PrimerosPuntos,
    KMediasMasMas(u64),
}

#[derive(Debug, Clone, PartialEq)]
pub struct KMedias {
    pub centros: Vec<Vec<f64>>,
    pub grupos: Vec<usize>,
    pub iteraciones: usize,
    pub convergio: bool,
    pub distorsion: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct KMediasSuaves {
    pub centros: Vec<Vec<f64>>,
    pub responsabilidades: Vec<Vec<f64>>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Union {
    pub primero: usize,
    pub segundo: usize,
    pub distancia: f64,
    pub miembros: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Componentes {
    pub media: Vec<f64>,
    pub primera: Vec<f64>,
    pub segunda: Vec<f64>,
    pub varianza_primera: f64,
    pub varianza_segunda: f64,
    pub varianza_total: f64,
    pub coordenadas: Vec<[f64; 2]>,
}

pub fn distancia_al_cuadrado(a: &[f64], b: &[f64]) -> f64 {
    let mut suma = 0.0;
    for (x, y) in a.iter().zip(b) {
        let diferencia = x - y;
        suma += diferencia * diferencia;
    }
    suma
}

pub fn distancia_euclidea(a: &[f64], b: &[f64]) -> f64 {
    distancia_al_cuadrado(a, b).sqrt()
}

pub fn producto_punto(a: &[f64], b: &[f64]) -> f64 {
    let mut suma = 0.0;
    for (x, y) in a.iter().zip(b) {
        suma += x * y;
    }
    suma
}

pub fn estandarizar(valores: &[f64]) -> Vec<f64> {
    if valores.is_empty() {
        return Vec::new();
    }
    let promedio = valores.iter().sum::<f64>() / valores.len() as f64;
    let mut centrados: Vec<f64> = valores.iter().map(|valor| valor - promedio).collect();
    let magnitud = valores
        .iter()
        .fold(0.0f64, |mayor, valor| mayor.max(valor.abs()));
    let desvio = centrados
        .iter()
        .fold(0.0f64, |mayor, valor| mayor.max(valor.abs()));
    if desvio <= magnitud * PARTE_DESPRECIABLE {
        return vec![0.0; valores.len()];
    }
    for valor in &mut centrados {
        *valor /= desvio;
    }
    let largo = producto_punto(&centrados, &centrados).sqrt();
    for valor in &mut centrados {
        *valor /= largo;
    }
    centrados
}

pub fn correlacion_pearson(a: &[f64], b: &[f64]) -> f64 {
    let largo = a.len().min(b.len());
    let primero = estandarizar(&a[..largo]);
    let segundo = estandarizar(&b[..largo]);
    producto_punto(&primero, &segundo).clamp(-1.0, 1.0)
}

pub fn distancia_pearson(a: &[f64], b: &[f64]) -> f64 {
    1.0 - correlacion_pearson(a, b)
}

pub fn matriz_de_correlaciones(puntos: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let estandarizados: Vec<Vec<f64>> = puntos.iter().map(|punto| estandarizar(punto)).collect();
    let cantidad = estandarizados.len();
    let mut correlaciones: Vec<Vec<f64>> = (0..cantidad)
        .into_par_iter()
        .map(|fila| {
            let mut valores = vec![0.0; cantidad];
            for columna in fila..cantidad {
                let correlacion = producto_punto(&estandarizados[fila], &estandarizados[columna]);
                valores[columna] = correlacion.clamp(-1.0, 1.0);
            }
            valores
        })
        .collect();
    copiar_al_otro_lado(&mut correlaciones);
    correlaciones
}

pub fn matriz_de_distancias(puntos: &[Vec<f64>], medida: Medida) -> Vec<Vec<f64>> {
    let cantidad = puntos.len();
    let mut distancias: Vec<Vec<f64>> = match medida {
        Medida::Euclidea => (0..cantidad)
            .into_par_iter()
            .map(|fila| {
                let mut valores = vec![0.0; cantidad];
                for columna in fila + 1..cantidad {
                    valores[columna] = distancia_euclidea(&puntos[fila], &puntos[columna]);
                }
                valores
            })
            .collect(),
        Medida::Pearson => matriz_de_correlaciones(puntos)
            .into_iter()
            .map(|fila| {
                fila.into_iter()
                    .map(|correlacion| 1.0 - correlacion)
                    .collect()
            })
            .collect(),
    };
    for (diagonal, fila) in distancias.iter_mut().enumerate() {
        fila[diagonal] = 0.0;
    }
    copiar_al_otro_lado(&mut distancias);
    distancias
}

fn copiar_al_otro_lado(matriz: &mut [Vec<f64>]) {
    for fila in 1..matriz.len() {
        let (filas_de_arriba, desde_esta) = matriz.split_at_mut(fila);
        let esta_fila = &mut desde_esta[0];
        for (columna, de_arriba) in filas_de_arriba.iter().enumerate() {
            esta_fila[columna] = de_arriba[fila];
        }
    }
}

pub fn revisar_puntos(puntos: &[Vec<f64>]) -> Result<usize, String> {
    let Some(primero) = puntos.first() else {
        return Err("No hay puntos.".to_string());
    };
    let dimension = primero.len();
    if dimension == 0 {
        return Err("Los puntos no tienen coordenadas.".to_string());
    }
    for (indice, punto) in puntos.iter().enumerate() {
        if punto.len() != dimension {
            return Err(format!(
                "El punto {indice} (contando desde 0) tiene {} coordenadas y el primero tiene {dimension}.",
                punto.len()
            ));
        }
        if punto.iter().any(|valor| !valor.is_finite()) {
            return Err(format!(
                "El punto {indice} (contando desde 0) tiene un valor que no es un número finito."
            ));
        }
    }
    Ok(dimension)
}

fn revisar_k(k: usize, cantidad: usize) -> Result<(), String> {
    if k == 0 {
        return Err("k debe ser al menos 1.".to_string());
    }
    if k > cantidad {
        return Err(format!(
            "k = {k} es mayor que la cantidad de puntos ({cantidad})."
        ));
    }
    Ok(())
}

fn revisar_matriz_cuadrada(matriz: &[Vec<f64>]) -> Result<(), String> {
    if matriz.is_empty() {
        return Err("La matriz está vacía.".to_string());
    }
    for fila in matriz {
        if fila.len() != matriz.len() {
            return Err(format!(
                "La matriz debe ser cuadrada: tiene {} filas y una fila con {} valores.",
                matriz.len(),
                fila.len()
            ));
        }
        if fila.iter().any(|valor| !valor.is_finite()) {
            return Err("La matriz tiene un valor que no es un número finito.".to_string());
        }
    }
    Ok(())
}

pub fn mas_lejano_primero(puntos: &[Vec<f64>], k: usize) -> Result<Vec<usize>, String> {
    revisar_puntos(puntos)?;
    revisar_k(k, puntos.len())?;
    let mut centros = vec![0];
    let mut distancia_al_centro = vec![f64::INFINITY; puntos.len()];
    acercar_al_nuevo_centro(&mut distancia_al_centro, puntos, &puntos[0]);
    distancia_al_centro[0] = YA_ES_CENTRO;
    while centros.len() < k {
        let mas_lejano = posicion_del_maximo(&distancia_al_centro);
        centros.push(mas_lejano);
        acercar_al_nuevo_centro(&mut distancia_al_centro, puntos, &puntos[mas_lejano]);
        distancia_al_centro[mas_lejano] = YA_ES_CENTRO;
    }
    Ok(centros)
}

pub fn acercar_al_nuevo_centro(
    distancia_al_centro: &mut [f64],
    puntos: &[Vec<f64>],
    centro: &[f64],
) {
    for (actual, punto) in distancia_al_centro.iter_mut().zip(puntos) {
        let nueva = distancia_euclidea(punto, centro);
        if nueva < *actual {
            *actual = nueva;
        }
    }
}

pub fn posicion_del_maximo(valores: &[f64]) -> usize {
    let mut mejor = 0;
    for (indice, &valor) in valores.iter().enumerate() {
        if valor > valores[mejor] {
            mejor = indice;
        }
    }
    mejor
}

pub fn radio_de_los_centros(puntos: &[Vec<f64>], centros: &[Vec<f64>]) -> f64 {
    let mut mayor = 0.0f64;
    for punto in puntos {
        mayor = mayor.max(menor_distancia_al_cuadrado(punto, centros));
    }
    mayor.sqrt()
}

pub fn centro_mas_cercano(punto: &[f64], centros: &[Vec<f64>]) -> usize {
    let mut mejor = 0;
    let mut menor = f64::INFINITY;
    for (indice, centro) in centros.iter().enumerate() {
        let cercania = distancia_al_cuadrado(punto, centro);
        if cercania < menor {
            menor = cercania;
            mejor = indice;
        }
    }
    mejor
}

fn menor_distancia_al_cuadrado(punto: &[f64], centros: &[Vec<f64>]) -> f64 {
    let mut menor = f64::INFINITY;
    for centro in centros {
        menor = menor.min(distancia_al_cuadrado(punto, centro));
    }
    menor
}

pub fn distorsion_cuadratica(puntos: &[Vec<f64>], centros: &[Vec<f64>]) -> f64 {
    if puntos.is_empty() {
        return 0.0;
    }
    let mut suma = 0.0;
    for punto in puntos {
        suma += menor_distancia_al_cuadrado(punto, centros);
    }
    suma / puntos.len() as f64
}

pub fn k_medias(
    puntos: &[Vec<f64>],
    k: usize,
    inicio: Inicio,
    maximo_iteraciones: usize,
) -> Result<KMedias, String> {
    revisar_puntos(puntos)?;
    revisar_k(k, puntos.len())?;
    let centros = centros_iniciales(puntos, k, inicio);
    Ok(lloyd(puntos, centros, maximo_iteraciones))
}

pub fn centros_iniciales(puntos: &[Vec<f64>], k: usize, inicio: Inicio) -> Vec<Vec<f64>> {
    let elegidos = match inicio {
        Inicio::PrimerosPuntos => (0..k.min(puntos.len())).collect(),
        Inicio::KMediasMasMas(semilla) => k_medias_mas_mas(puntos, k, semilla),
    };
    elegidos
        .into_iter()
        .map(|indice| puntos[indice].clone())
        .collect()
}

pub fn k_medias_mas_mas(puntos: &[Vec<f64>], k: usize, semilla: u64) -> Vec<usize> {
    if puntos.is_empty() || k == 0 {
        return Vec::new();
    }
    let mut azar = Azar::nuevo(semilla);
    let primero = azar.entero(puntos.len());
    let mut elegidos = vec![primero];
    let mut distancia_al_centro = vec![f64::INFINITY; puntos.len()];
    acercar_al_nuevo_centro(&mut distancia_al_centro, puntos, &puntos[primero]);
    while elegidos.len() < k.min(puntos.len()) {
        let pesos: Vec<f64> = distancia_al_centro
            .iter()
            .map(|lejania| lejania * lejania)
            .collect();
        let siguiente = azar.elegir_con_pesos(&pesos);
        elegidos.push(siguiente);
        acercar_al_nuevo_centro(&mut distancia_al_centro, puntos, &puntos[siguiente]);
    }
    elegidos
}

fn lloyd(puntos: &[Vec<f64>], mut centros: Vec<Vec<f64>>, maximo_iteraciones: usize) -> KMedias {
    let mut grupos = de_centros_a_grupos(puntos, &centros);
    let mut iteraciones = 0;
    let mut convergio = false;
    while iteraciones < maximo_iteraciones {
        iteraciones += 1;
        centros = de_grupos_a_centros(puntos, &grupos, &centros);
        let grupos_nuevos = de_centros_a_grupos(puntos, &centros);
        if grupos_nuevos == grupos {
            convergio = true;
            break;
        }
        grupos = grupos_nuevos;
    }
    let distorsion = distorsion_cuadratica(puntos, &centros);
    KMedias {
        centros,
        grupos,
        iteraciones,
        convergio,
        distorsion,
    }
}

pub fn de_centros_a_grupos(puntos: &[Vec<f64>], centros: &[Vec<f64>]) -> Vec<usize> {
    puntos
        .par_iter()
        .map(|punto| centro_mas_cercano(punto, centros))
        .collect()
}

pub fn de_grupos_a_centros(
    puntos: &[Vec<f64>],
    grupos: &[usize],
    centros: &[Vec<f64>],
) -> Vec<Vec<f64>> {
    let mut sumas: Vec<Vec<f64>> = centros
        .iter()
        .map(|centro| vec![0.0; centro.len()])
        .collect();
    let mut cantidades = vec![0usize; centros.len()];
    for (punto, &grupo) in puntos.iter().zip(grupos) {
        cantidades[grupo] += 1;
        for (suma, valor) in sumas[grupo].iter_mut().zip(punto) {
            *suma += valor;
        }
    }
    let mut nuevos = Vec::with_capacity(centros.len());
    for (grupo, suma) in sumas.into_iter().enumerate() {
        if cantidades[grupo] == 0 {
            nuevos.push(centros[grupo].clone());
        } else {
            let cantidad = cantidades[grupo] as f64;
            nuevos.push(suma.into_iter().map(|valor| valor / cantidad).collect());
        }
    }
    nuevos
}

pub fn curva_del_codo(
    puntos: &[Vec<f64>],
    k_maximo: usize,
    semilla: u64,
) -> Result<Vec<f64>, String> {
    revisar_puntos(puntos)?;
    if k_maximo == 0 {
        return Err("El k máximo debe ser al menos 1.".to_string());
    }
    let mut curva = Vec::new();
    for k in 1..=k_maximo.min(puntos.len()) {
        let mut mejor = f64::INFINITY;
        for intento in 0..INTENTOS_POR_K {
            let inicio = Inicio::KMediasMasMas(semilla.wrapping_add(intento));
            let resultado = k_medias(puntos, k, inicio, ITERACIONES_DEL_CODO)?;
            mejor = mejor.min(resultado.distorsion);
        }
        curva.push(mejor);
    }
    Ok(curva)
}

pub fn k_del_codo(curva: &[f64]) -> Option<usize> {
    let mut elegido = None;
    let mut mayor_frenada = 0.0;
    for k in 2..curva.len() {
        let caida_al_llegar = caida_relativa(curva[k - 2], curva[k - 1]);
        let caida_al_seguir = caida_relativa(curva[k - 1], curva[k]);
        let frenada = caida_al_llegar - caida_al_seguir;
        if frenada > mayor_frenada {
            mayor_frenada = frenada;
            elegido = Some(k);
        }
    }
    elegido
}

pub fn caida_relativa(antes: f64, despues: f64) -> f64 {
    if antes <= 0.0 {
        return 0.0;
    }
    (antes - despues) / antes
}

pub fn k_medias_suaves(
    puntos: &[Vec<f64>],
    k: usize,
    rigidez: f64,
    iteraciones: usize,
) -> Result<KMediasSuaves, String> {
    revisar_puntos(puntos)?;
    revisar_k(k, puntos.len())?;
    if !rigidez.is_finite() || rigidez < 0.0 {
        return Err("La rigidez β debe ser un número positivo.".to_string());
    }
    let mut centros: Vec<Vec<f64>> = puntos[..k].to_vec();
    for _ in 0..iteraciones {
        let responsabilidades = de_centros_a_grupos_suaves(puntos, &centros, rigidez);
        centros = de_grupos_suaves_a_centros(puntos, &responsabilidades, &centros);
    }
    let responsabilidades = de_centros_a_grupos_suaves(puntos, &centros, rigidez);
    Ok(KMediasSuaves {
        centros,
        responsabilidades,
    })
}

pub fn de_centros_a_grupos_suaves(
    puntos: &[Vec<f64>],
    centros: &[Vec<f64>],
    rigidez: f64,
) -> Vec<Vec<f64>> {
    puntos
        .par_iter()
        .map(|punto| responsabilidades_de_un_punto(punto, centros, rigidez))
        .collect()
}

fn responsabilidades_de_un_punto(punto: &[f64], centros: &[Vec<f64>], rigidez: f64) -> Vec<f64> {
    let distancias: Vec<f64> = centros
        .iter()
        .map(|centro| distancia_euclidea(punto, centro))
        .collect();
    let mas_cerca = distancias.iter().copied().fold(f64::INFINITY, f64::min);
    let mut atracciones: Vec<f64> = distancias
        .iter()
        .map(|lejania| (-rigidez * (lejania - mas_cerca)).exp())
        .collect();
    let total: f64 = atracciones.iter().sum();
    for atraccion in &mut atracciones {
        *atraccion /= total;
    }
    atracciones
}

pub fn de_grupos_suaves_a_centros(
    puntos: &[Vec<f64>],
    responsabilidades: &[Vec<f64>],
    centros: &[Vec<f64>],
) -> Vec<Vec<f64>> {
    centros
        .par_iter()
        .enumerate()
        .map(|(grupo, centro)| {
            let mut suma = vec![0.0; centro.len()];
            let mut peso_total = 0.0;
            for (punto, pesos) in puntos.iter().zip(responsabilidades) {
                let peso = pesos[grupo];
                peso_total += peso;
                for (acumulado, valor) in suma.iter_mut().zip(punto) {
                    *acumulado += peso * valor;
                }
            }
            if peso_total > 0.0 {
                suma.into_iter().map(|valor| valor / peso_total).collect()
            } else {
                centro.clone()
            }
        })
        .collect()
}

pub fn agrupamiento_jerarquico(distancias: &[Vec<f64>]) -> Result<Vec<Union>, String> {
    revisar_matriz_cuadrada(distancias)?;
    let cantidad = distancias.len();
    let mut entre_grupos = distancias.to_vec();
    let mut miembros: Vec<Vec<usize>> = (0..cantidad).map(|punto| vec![punto]).collect();
    let mut nodo: Vec<usize> = (0..cantidad).collect();
    let mut activo = vec![true; cantidad];
    let mut vecino: Vec<usize> = (0..cantidad)
        .map(|grupo| vecino_mas_cercano(grupo, &entre_grupos, &activo))
        .collect();
    let mut uniones = Vec::with_capacity(cantidad - 1);
    for paso in 0..cantidad - 1 {
        let (a, b) = par_mas_cercano(&entre_grupos, &vecino, &activo);
        let distancia = entre_grupos[a][b];
        let peso_a = miembros[a].len() as f64;
        let peso_b = miembros[b].len() as f64;
        for otro in 0..cantidad {
            if activo[otro] && otro != a && otro != b {
                let promedio = (peso_a * entre_grupos[a][otro] + peso_b * entre_grupos[b][otro])
                    / (peso_a + peso_b);
                entre_grupos[a][otro] = promedio;
                entre_grupos[otro][a] = promedio;
            }
        }
        activo[b] = false;
        let absorbidos = std::mem::take(&mut miembros[b]);
        miembros[a].extend(absorbidos);
        uniones.push(Union {
            primero: nodo[a],
            segundo: nodo[b],
            distancia,
            miembros: miembros[a].clone(),
        });
        nodo[a] = cantidad + paso;
        for otro in 0..cantidad {
            if activo[otro] && (otro == a || vecino[otro] == a || vecino[otro] == b) {
                vecino[otro] = vecino_mas_cercano(otro, &entre_grupos, &activo);
            }
        }
    }
    Ok(uniones)
}

pub fn vecino_mas_cercano(grupo: usize, entre_grupos: &[Vec<f64>], activo: &[bool]) -> usize {
    let mut mejor = grupo;
    let mut menor = f64::INFINITY;
    for (otro, &lejania) in entre_grupos[grupo].iter().enumerate() {
        if activo[otro] && otro != grupo && lejania < menor {
            menor = lejania;
            mejor = otro;
        }
    }
    mejor
}

pub fn par_mas_cercano(
    entre_grupos: &[Vec<f64>],
    vecino: &[usize],
    activo: &[bool],
) -> (usize, usize) {
    let mut mejor = (0, 0);
    let mut menor = f64::INFINITY;
    for (grupo, &cercano) in vecino.iter().enumerate() {
        if activo[grupo] && cercano != grupo && entre_grupos[grupo][cercano] < menor {
            menor = entre_grupos[grupo][cercano];
            mejor = (grupo, cercano);
        }
    }
    (mejor.0.min(mejor.1), mejor.0.max(mejor.1))
}

pub fn arbol_newick(uniones: &[Union], nombres: &[String]) -> String {
    let mut textos: Vec<String> = nombres.iter().map(|nombre| nombre_newick(nombre)).collect();
    let mut alturas = vec![0.0; nombres.len()];
    for union in uniones {
        let altura = union.distancia / 2.0;
        let primero = std::mem::take(&mut textos[union.primero]);
        let segundo = std::mem::take(&mut textos[union.segundo]);
        let rama_primero = numero_newick(altura - alturas[union.primero]);
        let rama_segundo = numero_newick(altura - alturas[union.segundo]);
        textos.push(format!(
            "({primero}:{rama_primero},{segundo}:{rama_segundo})"
        ));
        alturas.push(altura);
    }
    format!("{};", textos.last().cloned().unwrap_or_default())
}

fn nombre_newick(nombre: &str) -> String {
    nombre
        .chars()
        .map(|letra| {
            if "()[],:;' \t".contains(letra) {
                '_'
            } else {
                letra
            }
        })
        .collect()
}

fn numero_newick(longitud: f64) -> String {
    let positiva = if longitud > 0.0 { longitud } else { 0.0 };
    let texto = format!("{positiva:.6}");
    texto
        .trim_end_matches('0')
        .trim_end_matches('.')
        .to_string()
}

pub fn cast(similitudes: &[Vec<f64>], umbral: f64) -> Result<Vec<Vec<usize>>, String> {
    revisar_matriz_cuadrada(similitudes)?;
    if !umbral.is_finite() {
        return Err("El umbral θ debe ser un número.".to_string());
    }
    let cantidad = similitudes.len();
    let mut libre = vec![true; cantidad];
    let mut grado: Vec<usize> = (0..cantidad)
        .map(|elemento| vecinos_parecidos(elemento, similitudes, umbral, &libre))
        .collect();
    let mut grupos = Vec::new();
    let mut quedan = cantidad;
    while quedan > 0 {
        let comienzo = libre_de_mayor_grado(&grado, &libre);
        let grupo = formar_grupo_cast(comienzo, similitudes, umbral, &libre);
        for &elemento in &grupo {
            libre[elemento] = false;
            quedan -= 1;
        }
        for &elemento in &grupo {
            for otro in 0..cantidad {
                if libre[otro] && similitudes[otro][elemento] >= umbral {
                    grado[otro] -= 1;
                }
            }
        }
        grupos.push(grupo);
    }
    Ok(grupos)
}

pub fn vecinos_parecidos(
    elemento: usize,
    similitudes: &[Vec<f64>],
    umbral: f64,
    libre: &[bool],
) -> usize {
    let mut cuenta = 0;
    for (otro, &parecido) in similitudes[elemento].iter().enumerate() {
        if libre[otro] && otro != elemento && parecido >= umbral {
            cuenta += 1;
        }
    }
    cuenta
}

pub fn libre_de_mayor_grado(grado: &[usize], libre: &[bool]) -> usize {
    let mut mejor = None;
    for (elemento, &vecinos) in grado.iter().enumerate() {
        if !libre[elemento] {
            continue;
        }
        match mejor {
            Some(actual) if grado[actual] >= vecinos => {}
            _ => mejor = Some(elemento),
        }
    }
    mejor.unwrap_or(0)
}

fn formar_grupo_cast(
    comienzo: usize,
    similitudes: &[Vec<f64>],
    umbral: f64,
    libre: &[bool],
) -> Vec<usize> {
    let cantidad = similitudes.len();
    let mut dentro = vec![false; cantidad];
    let mut afinidad = vec![0.0; cantidad];
    dentro[comienzo] = true;
    sumar_similitudes(&mut afinidad, similitudes, comienzo, 1.0);
    let mut miembros = 1;
    let mut cambios_permitidos = 2 * libre.iter().filter(|&&esta| esta).count();
    while cambios_permitidos > 0 {
        cambios_permitidos -= 1;
        if let Some(cercano) = mas_parecido_afuera(&afinidad, &dentro, libre)
            && afinidad[cercano] / miembros as f64 >= umbral
        {
            dentro[cercano] = true;
            miembros += 1;
            sumar_similitudes(&mut afinidad, similitudes, cercano, 1.0);
            continue;
        }
        if let Some(lejano) = menos_parecido_adentro(&afinidad, &dentro)
            && miembros > 1
            && afinidad[lejano] / (miembros as f64) < umbral
        {
            dentro[lejano] = false;
            miembros -= 1;
            sumar_similitudes(&mut afinidad, similitudes, lejano, -1.0);
            continue;
        }
        break;
    }
    (0..cantidad).filter(|&elemento| dentro[elemento]).collect()
}

pub fn sumar_similitudes(
    afinidad: &mut [f64],
    similitudes: &[Vec<f64>],
    elemento: usize,
    signo: f64,
) {
    for (otro, valor) in afinidad.iter_mut().enumerate() {
        *valor += signo * similitudes[otro][elemento];
    }
}

pub fn mas_parecido_afuera(afinidad: &[f64], dentro: &[bool], libre: &[bool]) -> Option<usize> {
    let mut mejor: Option<usize> = None;
    for (elemento, &valor) in afinidad.iter().enumerate() {
        if dentro[elemento] || !libre[elemento] {
            continue;
        }
        match mejor {
            Some(actual) if afinidad[actual] >= valor => {}
            _ => mejor = Some(elemento),
        }
    }
    mejor
}

pub fn menos_parecido_adentro(afinidad: &[f64], dentro: &[bool]) -> Option<usize> {
    let mut peor: Option<usize> = None;
    for (elemento, &valor) in afinidad.iter().enumerate() {
        if !dentro[elemento] {
            continue;
        }
        match peor {
            Some(actual) if afinidad[actual] <= valor => {}
            _ => peor = Some(elemento),
        }
    }
    peor
}

pub fn componentes_principales(puntos: &[Vec<f64>]) -> Result<Componentes, String> {
    let dimension = revisar_puntos(puntos)?;
    let media = centro_de_gravedad(puntos, dimension);
    let covarianza = matriz_de_covarianza(puntos, &media);
    let varianza_total = (0..dimension).map(|eje| covarianza[eje][eje]).sum();
    let primera = vector_propio_dominante(&covarianza, &[]);
    let segunda = vector_propio_dominante(&covarianza, &primera);
    let varianza_primera = cociente_de_rayleigh(&covarianza, &primera);
    let varianza_segunda = cociente_de_rayleigh(&covarianza, &segunda);
    let mut componentes = Componentes {
        media,
        primera,
        segunda,
        varianza_primera,
        varianza_segunda,
        varianza_total,
        coordenadas: Vec::with_capacity(puntos.len()),
    };
    for punto in puntos {
        let coordenada = proyectar(punto, &componentes);
        componentes.coordenadas.push(coordenada);
    }
    Ok(componentes)
}

pub fn centro_de_gravedad(puntos: &[Vec<f64>], dimension: usize) -> Vec<f64> {
    let mut suma = vec![0.0; dimension];
    for punto in puntos {
        for (acumulado, valor) in suma.iter_mut().zip(punto) {
            *acumulado += valor;
        }
    }
    let cantidad = puntos.len().max(1) as f64;
    suma.into_iter().map(|valor| valor / cantidad).collect()
}

pub fn matriz_de_covarianza(puntos: &[Vec<f64>], media: &[f64]) -> Vec<Vec<f64>> {
    let dimension = media.len();
    let mut covarianza = vec![vec![0.0; dimension]; dimension];
    let mut centrado = vec![0.0; dimension];
    for punto in puntos {
        for ((destino, valor), promedio) in centrado.iter_mut().zip(punto).zip(media) {
            *destino = valor - promedio;
        }
        for (fila, valores) in covarianza.iter_mut().enumerate() {
            let lado = centrado[fila];
            for (columna, acumulado) in valores.iter_mut().enumerate().skip(fila) {
                *acumulado += lado * centrado[columna];
            }
        }
    }
    let divisor = (puntos.len().max(2) - 1) as f64;
    (0..dimension)
        .map(|fila| {
            (0..dimension)
                .map(|columna| covarianza[fila.min(columna)][fila.max(columna)] / divisor)
                .collect()
        })
        .collect()
}

fn vector_propio_dominante(matriz: &[Vec<f64>], distinto_de: &[f64]) -> Vec<f64> {
    let dimension = matriz.len();
    let escala: f64 = (0..dimension).map(|eje| matriz[eje][eje].abs()).sum();
    let mut azar = Azar::nuevo(SEMILLA_DE_POTENCIA);
    let mut vector: Vec<f64> = (0..dimension).map(|_| azar.normal()).collect();
    quitar_direccion(&mut vector, distinto_de);
    if normalizar(&mut vector) <= PARTE_DESPRECIABLE {
        return vec![0.0; dimension];
    }
    for _ in 0..VUELTAS_DE_POTENCIA {
        let mut siguiente = multiplicar(matriz, &vector);
        quitar_direccion(&mut siguiente, distinto_de);
        if normalizar(&mut siguiente) <= escala * PARTE_DESPRECIABLE {
            break;
        }
        let cambio = distancia_al_cuadrado(&siguiente, &vector);
        vector = siguiente;
        if cambio < CAMBIO_DESPRECIABLE {
            break;
        }
    }
    orientar(&mut vector);
    vector
}

pub fn multiplicar(matriz: &[Vec<f64>], vector: &[f64]) -> Vec<f64> {
    matriz
        .iter()
        .map(|fila| producto_punto(fila, vector))
        .collect()
}

pub fn quitar_direccion(vector: &mut [f64], direccion: &[f64]) {
    let parte = producto_punto(vector, direccion);
    for (valor, eje) in vector.iter_mut().zip(direccion) {
        *valor -= parte * eje;
    }
}

pub fn normalizar(vector: &mut [f64]) -> f64 {
    let largo = producto_punto(vector, vector).sqrt();
    if largo > 0.0 && largo.is_finite() {
        for valor in vector.iter_mut() {
            *valor /= largo;
        }
    }
    largo
}

pub fn orientar(vector: &mut [f64]) {
    let mut mayor = 0.0f64;
    for &valor in vector.iter() {
        if valor.abs() > mayor.abs() {
            mayor = valor;
        }
    }
    if mayor < 0.0 {
        for valor in vector.iter_mut() {
            *valor = -*valor;
        }
    }
}

pub fn cociente_de_rayleigh(matriz: &[Vec<f64>], vector: &[f64]) -> f64 {
    producto_punto(vector, &multiplicar(matriz, vector))
}

pub fn proyectar(punto: &[f64], componentes: &Componentes) -> [f64; 2] {
    let centrado: Vec<f64> = punto
        .iter()
        .zip(&componentes.media)
        .map(|(valor, promedio)| valor - promedio)
        .collect();
    [
        producto_punto(&centrado, &componentes.primera),
        producto_punto(&centrado, &componentes.segunda),
    ]
}
