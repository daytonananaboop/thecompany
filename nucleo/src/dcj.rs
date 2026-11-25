pub const TELOMERO: usize = 0;

pub type Vertice = (usize, usize);

pub const NADA: Vertice = (0, 0);

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Cromosoma {
    pub lineal: bool,
    pub bloques: Vec<i64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Clase {
    Ciclo,
    CaminoAB,
    CaminoAA,
    CaminoBB,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Componente {
    pub clase: Clase,
    pub extremos: Vec<usize>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Conteo {
    pub bloques: usize,
    pub ciclos: usize,
    pub caminos_impares: usize,
    pub caminos_pares: usize,
    pub distancia: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TipoDcj {
    Inversion,
    Translocacion,
    Fusion,
    Fision,
    Circularizacion,
    Linealizacion,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Operacion {
    pub cortadas: [Vertice; 2],
    pub nuevas: [Vertice; 2],
    pub tipo: TipoDcj,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ordenamiento {
    pub etiquetas: Vec<u64>,
    pub pareja_a: Vec<usize>,
    pub pareja_b: Vec<usize>,
    pub componentes: Vec<Componente>,
    pub conteo: Conteo,
    pub operaciones: Vec<Operacion>,
}

pub fn leer_genoma(texto: &str, limite: usize) -> Result<Vec<Cromosoma>, String> {
    let mut genoma = Vec::new();
    let mut abierto: Option<Cromosoma> = None;
    let mut cantidad = 0;
    for pieza in piezas(texto) {
        match pieza.as_str() {
            "[" | "(" => {
                if abierto.is_some() {
                    return Err(
                        "Hay un cromosoma dentro de otro: cierra cada uno antes de abrir el siguiente."
                            .to_string(),
                    );
                }
                abierto = Some(Cromosoma {
                    lineal: pieza == "[",
                    bloques: Vec::new(),
                });
            }
            "]" | ")" => {
                let cromosoma = abierto
                    .take()
                    .ok_or_else(|| format!("Sobra un «{pieza}»: no hay cromosoma abierto."))?;
                if cromosoma.lineal != (pieza == "]") {
                    return Err(
                        "Un cromosoma lineal se abre con «[» y se cierra con «]»; uno circular, con «(» y «)»."
                            .to_string(),
                    );
                }
                if cromosoma.bloques.is_empty() {
                    return Err("Hay un cromosoma vacío.".to_string());
                }
                genoma.push(cromosoma);
            }
            _ => {
                let bloque = pieza.parse::<i64>().map_err(|_| {
                    format!("«{pieza}» no es un bloque: escribe enteros con signo, como +3 o -2.")
                })?;
                cantidad += 1;
                if cantidad > limite {
                    return Err(format!(
                        "Hay más de {limite} bloques; aquí se admiten hasta {limite}."
                    ));
                }
                let cromosoma = abierto.as_mut().ok_or_else(|| {
                    format!(
                        "«{pieza}» quedó fuera de un cromosoma: los lineales van entre [ ] y los circulares entre ( )."
                    )
                })?;
                cromosoma.bloques.push(bloque);
            }
        }
    }
    if abierto.is_some() {
        return Err("Falta cerrar un cromosoma.".to_string());
    }
    if genoma.is_empty() {
        return Err(
            "No hay ningún cromosoma: escribe los lineales entre corchetes, como [+1 -2], y los circulares entre paréntesis, como (+3 +4)."
                .to_string(),
        );
    }
    Ok(genoma)
}

fn piezas(texto: &str) -> Vec<String> {
    texto
        .replace('[', " [ ")
        .replace(']', " ] ")
        .replace('(', " ( ")
        .replace(')', " ) ")
        .replace(['−', '–'], "-")
        .replace([',', ';'], " ")
        .split_whitespace()
        .map(str::to_string)
        .collect()
}

pub fn etiquetas_comunes(a: &[Cromosoma], b: &[Cromosoma]) -> Result<Vec<u64>, String> {
    let etiquetas_a = etiquetas_del_genoma(a).map_err(|error| format!("Genoma A. {error}"))?;
    let etiquetas_b = etiquetas_del_genoma(b).map_err(|error| format!("Genoma B. {error}"))?;
    if let Some(error) = primera_diferencia(&etiquetas_a, &etiquetas_b) {
        return Err(error);
    }
    Ok(etiquetas_a)
}

fn etiquetas_del_genoma(genoma: &[Cromosoma]) -> Result<Vec<u64>, String> {
    if genoma.is_empty() {
        return Err("No tiene cromosomas.".to_string());
    }
    let mut etiquetas = Vec::new();
    for cromosoma in genoma {
        if cromosoma.bloques.is_empty() {
            return Err("Hay un cromosoma vacío.".to_string());
        }
        for &bloque in &cromosoma.bloques {
            if bloque == 0 {
                return Err("El bloque 0 no existe: numera los bloques desde 1.".to_string());
            }
            if bloque == i64::MIN {
                return Err(format!(
                    "El bloque {bloque} es demasiado grande: usa números entre -{} y {}.",
                    i64::MAX,
                    i64::MAX
                ));
            }
            etiquetas.push(bloque.unsigned_abs());
        }
    }
    etiquetas.sort_unstable();
    for par in etiquetas.windows(2) {
        if par[0] == par[1] {
            return Err(format!("El bloque {} aparece más de una vez.", par[0]));
        }
    }
    Ok(etiquetas)
}

fn primera_diferencia(a: &[u64], b: &[u64]) -> Option<String> {
    let mut i = 0;
    let mut j = 0;
    while i < a.len() || j < b.len() {
        let de_a = a.get(i).copied().unwrap_or(u64::MAX);
        let de_b = b.get(j).copied().unwrap_or(u64::MAX);
        if de_a < de_b {
            return Some(format!(
                "El bloque {de_a} está en el genoma A pero no en el B: los dos deben tener los mismos bloques."
            ));
        }
        if de_b < de_a {
            return Some(format!(
                "El bloque {de_b} está en el genoma B pero no en el A: los dos deben tener los mismos bloques."
            ));
        }
        i += 1;
        j += 1;
    }
    None
}

pub fn es_cola(extremo: usize) -> bool {
    extremo % 2 == 1
}

pub fn otro_extremo(extremo: usize) -> usize {
    if es_cola(extremo) {
        extremo + 1
    } else {
        extremo - 1
    }
}

pub fn bloque_del_extremo(extremo: usize) -> usize {
    extremo.div_ceil(2)
}

fn indice(bloque: i64, etiquetas: &[u64]) -> Result<usize, String> {
    let numero = bloque.unsigned_abs();
    let seguidas = etiquetas.last() == Some(&(etiquetas.len() as u64));
    if seguidas && numero >= 1 && numero <= etiquetas.len() as u64 {
        return Ok(numero as usize);
    }
    etiquetas
        .binary_search(&numero)
        .map(|posicion| posicion + 1)
        .map_err(|_| format!("El bloque {numero} no está en la lista de bloques."))
}

fn extremo_de_entrada(bloque: i64, etiquetas: &[u64]) -> Result<usize, String> {
    let i = indice(bloque, etiquetas)?;
    Ok(if bloque > 0 { 2 * i - 1 } else { 2 * i })
}

fn extremo_de_salida(bloque: i64, etiquetas: &[u64]) -> Result<usize, String> {
    Ok(otro_extremo(extremo_de_entrada(bloque, etiquetas)?))
}

pub fn adyacencias(genoma: &[Cromosoma], etiquetas: &[u64]) -> Result<Vec<usize>, String> {
    let mut pareja = vec![TELOMERO; 2 * etiquetas.len() + 1];
    for cromosoma in genoma {
        for vecinos in cromosoma.bloques.windows(2) {
            let salida = extremo_de_salida(vecinos[0], etiquetas)?;
            let entrada = extremo_de_entrada(vecinos[1], etiquetas)?;
            pareja[salida] = entrada;
            pareja[entrada] = salida;
        }
        if !cromosoma.lineal {
            let ultimo = cromosoma.bloques[cromosoma.bloques.len() - 1];
            let salida = extremo_de_salida(ultimo, etiquetas)?;
            let entrada = extremo_de_entrada(cromosoma.bloques[0], etiquetas)?;
            pareja[salida] = entrada;
            pareja[entrada] = salida;
        }
    }
    Ok(pareja)
}

pub fn grafo_de_adyacencias(pareja_a: &[usize], pareja_b: &[usize]) -> Vec<Componente> {
    let mut visitado = vec![false; pareja_a.len()];
    let mut componentes = Vec::new();
    for extremo in 1..pareja_a.len() {
        if pareja_a[extremo] == TELOMERO && !visitado[extremo] {
            componentes.push(recorrer_camino(
                extremo,
                true,
                pareja_a,
                pareja_b,
                &mut visitado,
            ));
        }
    }
    for extremo in 1..pareja_b.len() {
        if pareja_b[extremo] == TELOMERO && !visitado[extremo] {
            componentes.push(recorrer_camino(
                extremo,
                false,
                pareja_a,
                pareja_b,
                &mut visitado,
            ));
        }
    }
    for extremo in 1..pareja_a.len() {
        if !visitado[extremo] {
            componentes.push(recorrer_ciclo(extremo, pareja_a, pareja_b, &mut visitado));
        }
    }
    componentes
}

fn recorrer_camino(
    inicio: usize,
    empieza_en_a: bool,
    pareja_a: &[usize],
    pareja_b: &[usize],
    visitado: &mut [bool],
) -> Componente {
    let mut extremos = vec![inicio];
    visitado[inicio] = true;
    let mut actual = inicio;
    let mut buscar_en_a = !empieza_en_a;
    loop {
        let siguiente = if buscar_en_a {
            pareja_a[actual]
        } else {
            pareja_b[actual]
        };
        if siguiente == TELOMERO {
            break;
        }
        visitado[siguiente] = true;
        extremos.push(siguiente);
        actual = siguiente;
        buscar_en_a = !buscar_en_a;
    }
    let clase = match (empieza_en_a, buscar_en_a) {
        (true, true) => Clase::CaminoAA,
        (false, false) => Clase::CaminoBB,
        _ => Clase::CaminoAB,
    };
    Componente { clase, extremos }
}

fn recorrer_ciclo(
    inicio: usize,
    pareja_a: &[usize],
    pareja_b: &[usize],
    visitado: &mut [bool],
) -> Componente {
    let mut extremos = Vec::new();
    let mut actual = inicio;
    let mut buscar_en_a = true;
    loop {
        visitado[actual] = true;
        extremos.push(actual);
        actual = if buscar_en_a {
            pareja_a[actual]
        } else {
            pareja_b[actual]
        };
        buscar_en_a = !buscar_en_a;
        if actual == inicio {
            break;
        }
    }
    Componente {
        clase: Clase::Ciclo,
        extremos,
    }
}

pub fn contar(bloques: usize, componentes: &[Componente]) -> Conteo {
    let mut conteo = Conteo {
        bloques,
        ..Conteo::default()
    };
    for componente in componentes {
        match componente.clase {
            Clase::Ciclo => conteo.ciclos += 1,
            Clase::CaminoAB => conteo.caminos_impares += 1,
            Clase::CaminoAA | Clase::CaminoBB => conteo.caminos_pares += 1,
        }
    }
    conteo.distancia = bloques - conteo.ciclos - conteo.caminos_impares / 2;
    conteo
}

pub fn distancia_dcj(a: &[Cromosoma], b: &[Cromosoma]) -> Result<Conteo, String> {
    let etiquetas = etiquetas_comunes(a, b)?;
    let pareja_a = adyacencias(a, &etiquetas)?;
    let pareja_b = adyacencias(b, &etiquetas)?;
    let componentes = grafo_de_adyacencias(&pareja_a, &pareja_b);
    Ok(contar(etiquetas.len(), &componentes))
}

pub fn ordenar(a: &[Cromosoma], b: &[Cromosoma]) -> Result<Ordenamiento, String> {
    let etiquetas = etiquetas_comunes(a, b)?;
    let pareja_a = adyacencias(a, &etiquetas)?;
    let pareja_b = adyacencias(b, &etiquetas)?;
    let componentes = grafo_de_adyacencias(&pareja_a, &pareja_b);
    let conteo = contar(etiquetas.len(), &componentes);
    let operaciones = escenario(pareja_a.clone(), &pareja_b);
    Ok(Ordenamiento {
        etiquetas,
        pareja_a,
        pareja_b,
        componentes,
        conteo,
        operaciones,
    })
}

pub fn escenario(mut pareja: Vec<usize>, objetivo: &[usize]) -> Vec<Operacion> {
    let mut operaciones = Vec::new();
    for p in 1..objetivo.len() {
        let q = objetivo[p];
        if q != TELOMERO && p < q && pareja[p] != q {
            operaciones.push(unir(&mut pareja, p, q));
        }
    }
    for p in 1..objetivo.len() {
        if objetivo[p] == TELOMERO && pareja[p] != TELOMERO {
            operaciones.push(cortar(&mut pareja, p));
        }
    }
    operaciones
}

fn vertice(p: usize, q: usize) -> Vertice {
    if p == TELOMERO {
        (q, TELOMERO)
    } else if q == TELOMERO {
        (p, TELOMERO)
    } else {
        (p.min(q), p.max(q))
    }
}

fn unir(pareja: &mut [usize], p: usize, q: usize) -> Operacion {
    let x = pareja[p];
    let y = pareja[q];
    let antes = cromosomas_tocados(pareja, &[p, q]);
    let cortadas = [vertice(p, x), vertice(q, y)];
    let nuevas = [vertice(p, q), vertice(x, y)];
    aplicar(pareja, &cortadas, &nuevas);
    let despues = cromosomas_tocados(pareja, &[p, q, x, y]);
    Operacion {
        cortadas,
        nuevas,
        tipo: clasificar(&antes, &despues),
    }
}

fn cortar(pareja: &mut [usize], p: usize) -> Operacion {
    let q = pareja[p];
    let antes = cromosomas_tocados(pareja, &[p]);
    let cortadas = [vertice(p, q), NADA];
    let nuevas = [(p, TELOMERO), (q, TELOMERO)];
    aplicar(pareja, &cortadas, &nuevas);
    let despues = cromosomas_tocados(pareja, &[p, q]);
    Operacion {
        cortadas,
        nuevas,
        tipo: clasificar(&antes, &despues),
    }
}

pub fn aplicar(pareja: &mut [usize], cortadas: &[Vertice; 2], nuevas: &[Vertice; 2]) {
    for &(u, v) in cortadas {
        pareja[u] = TELOMERO;
        pareja[v] = TELOMERO;
    }
    for &(u, v) in nuevas {
        if v != TELOMERO {
            pareja[u] = v;
            pareja[v] = u;
        }
    }
}

fn cromosomas_tocados(pareja: &[usize], extremos: &[usize]) -> Vec<(usize, bool)> {
    let mut tocados = Vec::with_capacity(4);
    for &extremo in extremos {
        if extremo == TELOMERO {
            continue;
        }
        let cromosoma = cromosoma_de(pareja, extremo);
        if !tocados.contains(&cromosoma) {
            tocados.push(cromosoma);
        }
    }
    tocados
}

fn cromosoma_de(pareja: &[usize], inicio: usize) -> (usize, bool) {
    let mut menor = inicio;
    let mut actual = inicio;
    loop {
        let salida = otro_extremo(actual);
        menor = menor.min(salida);
        let siguiente = pareja[salida];
        if siguiente == TELOMERO {
            break;
        }
        if siguiente == inicio {
            return (menor, false);
        }
        menor = menor.min(siguiente);
        actual = siguiente;
    }
    let mut actual = inicio;
    while pareja[actual] != TELOMERO {
        let vecino = pareja[actual];
        actual = otro_extremo(vecino);
        menor = menor.min(vecino).min(actual);
    }
    (menor, true)
}

fn clasificar(antes: &[(usize, bool)], despues: &[(usize, bool)]) -> TipoDcj {
    let lineales_antes = antes.iter().filter(|(_, lineal)| *lineal).count();
    let lineales_despues = despues.iter().filter(|(_, lineal)| *lineal).count();
    if despues.len() > antes.len() {
        TipoDcj::Fision
    } else if despues.len() < antes.len() {
        TipoDcj::Fusion
    } else if lineales_despues > lineales_antes {
        TipoDcj::Linealizacion
    } else if lineales_despues < lineales_antes {
        TipoDcj::Circularizacion
    } else if antes.len() == 1 {
        TipoDcj::Inversion
    } else {
        TipoDcj::Translocacion
    }
}

pub fn genoma_desde_adyacencias(pareja: &[usize], etiquetas: &[u64]) -> Vec<Cromosoma> {
    let mut visitado = vec![false; etiquetas.len() + 1];
    let mut genoma = Vec::new();
    for extremo in 1..pareja.len() {
        if pareja[extremo] == TELOMERO && !visitado[bloque_del_extremo(extremo)] {
            genoma.push(Cromosoma {
                lineal: true,
                bloques: leer_cromosoma(extremo, pareja, etiquetas, &mut visitado),
            });
        }
    }
    for bloque in 1..=etiquetas.len() {
        if !visitado[bloque] {
            genoma.push(Cromosoma {
                lineal: false,
                bloques: leer_cromosoma(2 * bloque - 1, pareja, etiquetas, &mut visitado),
            });
        }
    }
    genoma
}

fn leer_cromosoma(
    inicio: usize,
    pareja: &[usize],
    etiquetas: &[u64],
    visitado: &mut [bool],
) -> Vec<i64> {
    let mut bloques = Vec::new();
    let mut entrada = inicio;
    while entrada != TELOMERO && !visitado[bloque_del_extremo(entrada)] {
        let bloque = bloque_del_extremo(entrada);
        visitado[bloque] = true;
        let etiqueta = etiquetas[bloque - 1] as i64;
        bloques.push(if es_cola(entrada) {
            etiqueta
        } else {
            -etiqueta
        });
        entrada = pareja[otro_extremo(entrada)];
    }
    bloques
}

pub fn nombre_del_extremo(extremo: usize, etiquetas: &[u64]) -> String {
    let etiqueta = etiquetas
        .get(bloque_del_extremo(extremo).wrapping_sub(1))
        .copied()
        .unwrap_or(0);
    if es_cola(extremo) {
        format!("{etiqueta}t")
    } else {
        format!("{etiqueta}h")
    }
}
