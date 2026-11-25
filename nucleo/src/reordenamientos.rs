use crate::secuencia::indice_base;
use std::collections::{HashMap, HashSet};

pub type Genoma = Vec<Vec<i64>>;

pub const MAXIMO_PARES: usize = 2_000_000;

pub const NINGUNA: usize = usize::MAX;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DosRoturas {
    pub nodos: [usize; 4],
    pub genoma: Genoma,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ParCompartido {
    pub posicion_a: usize,
    pub posicion_b: usize,
    pub invertido: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BloqueSintenico {
    pub inicio_a: usize,
    pub fin_a: usize,
    pub inicio_b: usize,
    pub fin_b: usize,
    pub invertido: bool,
    pub kmeros: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sintenia {
    pub pares: Vec<ParCompartido>,
    pub bloques: Vec<BloqueSintenico>,
}

pub fn validar_bloques(bloques: &[i64]) -> Result<(), String> {
    let cantidad = bloques.len();
    if cantidad == 0 {
        return Err("No hay ningún bloque.".to_string());
    }
    let mut visto = vec![false; cantidad + 1];
    for &bloque in bloques {
        let numero = bloque.unsigned_abs();
        if numero == 0 {
            return Err("El bloque 0 no existe: los bloques se numeran desde 1.".to_string());
        }
        if numero > cantidad as u64 {
            return Err(format!(
                "Los bloques deben ser los números del 1 al {cantidad}, tantos como bloques hay, pero aparece el {numero}. ¿Falta algún bloque?"
            ));
        }
        if visto[numero as usize] {
            return Err(format!("El bloque {numero} aparece más de una vez."));
        }
        visto[numero as usize] = true;
    }
    Ok(())
}

pub fn ordenamiento_voraz(permutacion: &[i64]) -> Result<Vec<Vec<i64>>, String> {
    validar_bloques(permutacion)?;
    let mut actual = permutacion.to_vec();
    let mut permutaciones = Vec::new();
    for lugar in 0..actual.len() {
        let bloque = lugar as i64 + 1;
        if actual[lugar].abs() != bloque {
            let donde = posicion_del_bloque(&actual, bloque, lugar);
            invertir_tramo(&mut actual, lugar, donde);
            permutaciones.push(actual.clone());
        }
        if actual[lugar] == -bloque {
            actual[lugar] = bloque;
            permutaciones.push(actual.clone());
        }
    }
    Ok(permutaciones)
}

pub fn posicion_del_bloque(permutacion: &[i64], bloque: i64, desde: usize) -> usize {
    let mut posicion = desde;
    while permutacion[posicion].abs() != bloque {
        posicion += 1;
    }
    posicion
}

pub fn invertir_tramo(permutacion: &mut [i64], inicio: usize, fin: usize) {
    permutacion[inicio..=fin].reverse();
    for valor in &mut permutacion[inicio..=fin] {
        *valor = -*valor;
    }
}

pub fn numero_de_rupturas(permutacion: &[i64]) -> Result<usize, String> {
    Ok(rupturas(permutacion)?.len())
}

pub fn rupturas(permutacion: &[i64]) -> Result<Vec<(i64, i64)>, String> {
    validar_bloques(permutacion)?;
    let mut extendida = Vec::with_capacity(permutacion.len() + 2);
    extendida.push(0);
    extendida.extend_from_slice(permutacion);
    extendida.push(permutacion.len() as i64 + 1);
    let mut encontradas = Vec::new();
    for vecinos in extendida.windows(2) {
        if vecinos[1] - vecinos[0] != 1 {
            encontradas.push((vecinos[0], vecinos[1]));
        }
    }
    Ok(encontradas)
}

pub fn validar_genoma(genoma: &[Vec<i64>]) -> Result<usize, String> {
    if genoma.is_empty() {
        return Err("El genoma no tiene cromosomas.".to_string());
    }
    if genoma.iter().any(|cromosoma| cromosoma.is_empty()) {
        return Err("Hay un cromosoma sin bloques.".to_string());
    }
    let todos = genoma.concat();
    validar_bloques(&todos)?;
    Ok(todos.len())
}

fn validar_mismos_bloques(p: &[Vec<i64>], q: &[Vec<i64>]) -> Result<usize, String> {
    let bloques_p = validar_genoma(p).map_err(|error| format!("Genoma P. {error}"))?;
    let bloques_q = validar_genoma(q).map_err(|error| format!("Genoma Q. {error}"))?;
    if bloques_p != bloques_q {
        return Err(format!(
            "P tiene {bloques_p} bloques y Q tiene {bloques_q}: los dos genomas deben tener los mismos bloques."
        ));
    }
    Ok(bloques_p)
}

fn nodo_cola(bloque: i64) -> usize {
    2 * bloque.unsigned_abs() as usize - 1
}

fn nodo_cabeza(bloque: i64) -> usize {
    2 * bloque.unsigned_abs() as usize
}

pub fn es_cola(nodo: usize) -> bool {
    nodo % 2 == 1
}

pub fn otro_extremo(nodo: usize) -> usize {
    if es_cola(nodo) { nodo + 1 } else { nodo - 1 }
}

pub fn bloque_del_nodo(nodo: usize) -> usize {
    nodo.div_ceil(2)
}

pub fn cromosoma_a_ciclo(cromosoma: &[i64]) -> Result<Vec<usize>, String> {
    validar_bloques(cromosoma)?;
    Ok(nodos_del_cromosoma(cromosoma))
}

pub fn nodos_del_cromosoma(cromosoma: &[i64]) -> Vec<usize> {
    let mut nodos = Vec::with_capacity(2 * cromosoma.len());
    for &bloque in cromosoma {
        let cola = nodo_cola(bloque);
        let cabeza = nodo_cabeza(bloque);
        if bloque > 0 {
            nodos.push(cola);
            nodos.push(cabeza);
        } else {
            nodos.push(cabeza);
            nodos.push(cola);
        }
    }
    nodos
}

pub fn ciclo_a_cromosoma(nodos: &[usize]) -> Result<Vec<i64>, String> {
    let cantidad = nodos.len();
    if cantidad == 0 || cantidad % 2 == 1 {
        return Err("El ciclo debe tener una cantidad par de nodos: dos por bloque.".to_string());
    }
    if let Some(fuera) = nodos.iter().find(|&&nodo| nodo == 0 || nodo > cantidad) {
        return Err(format!(
            "Los nodos deben ser los números del 1 al {cantidad}, tantos como nodos hay, pero aparece el {fuera}."
        ));
    }
    let mut cromosoma = Vec::with_capacity(cantidad / 2);
    for extremos in nodos.chunks(2) {
        let entrada = extremos[0];
        let salida = extremos[1];
        if salida != otro_extremo(entrada) {
            return Err(format!(
                "Los nodos {entrada} y {salida} no son los dos extremos de un mismo bloque (2i−1 y 2i)."
            ));
        }
        let bloque = bloque_del_nodo(entrada) as i64;
        if es_cola(entrada) {
            cromosoma.push(bloque);
        } else {
            cromosoma.push(-bloque);
        }
    }
    validar_bloques(&cromosoma)?;
    Ok(cromosoma)
}

pub fn aristas_de_color(genoma: &[Vec<i64>]) -> Result<Vec<(usize, usize)>, String> {
    validar_genoma(genoma)?;
    Ok(aristas_entre_vecinos(genoma))
}

pub fn aristas_entre_vecinos(genoma: &[Vec<i64>]) -> Vec<(usize, usize)> {
    let mut aristas = Vec::new();
    for cromosoma in genoma {
        let nodos = nodos_del_cromosoma(cromosoma);
        let total = nodos.len();
        for j in 0..cromosoma.len() {
            let salida = nodos[2 * j + 1];
            let entrada_del_siguiente = nodos[(2 * j + 2) % total];
            aristas.push((salida, entrada_del_siguiente));
        }
    }
    aristas
}

pub fn pareja_de_cada_nodo(aristas: &[(usize, usize)], bloques: usize) -> Vec<usize> {
    let mut pareja = vec![0; 2 * bloques + 1];
    for &(a, b) in aristas {
        pareja[a] = b;
        pareja[b] = a;
    }
    pareja
}

pub fn grafo_a_genoma(aristas: &[(usize, usize)]) -> Result<Genoma, String> {
    let bloques = aristas.len();
    if bloques == 0 {
        return Err("No hay aristas de color.".to_string());
    }
    let nodos = 2 * bloques;
    let mut usado = vec![false; nodos + 1];
    for &(a, b) in aristas {
        if a == b {
            return Err(format!("La arista ({a}, {b}) une un nodo consigo mismo."));
        }
        for nodo in [a, b] {
            if nodo == 0 || nodo > nodos {
                return Err(format!(
                    "Los nodos deben ser los números del 1 al {nodos}, dos por arista, pero aparece el {nodo}."
                ));
            }
            if usado[nodo] {
                return Err(format!(
                    "El nodo {nodo} está en más de una arista de color."
                ));
            }
            usado[nodo] = true;
        }
    }
    let pareja = pareja_de_cada_nodo(aristas, bloques);
    Ok(recorrer_genoma(&pareja, bloques))
}

fn recorrer_genoma(pareja: &[usize], bloques: usize) -> Genoma {
    let mut visitado = vec![false; bloques + 1];
    let mut genoma = Vec::new();
    for primero in 1..=bloques {
        if visitado[primero] {
            continue;
        }
        let mut cromosoma = Vec::new();
        let mut entrada = 2 * primero - 1;
        while !visitado[bloque_del_nodo(entrada)] {
            let bloque = bloque_del_nodo(entrada);
            visitado[bloque] = true;
            if es_cola(entrada) {
                cromosoma.push(bloque as i64);
            } else {
                cromosoma.push(-(bloque as i64));
            }
            let salida = otro_extremo(entrada);
            entrada = pareja[salida];
        }
        genoma.push(cromosoma);
    }
    genoma
}

pub fn dos_roturas_en_grafo(
    aristas: &[(usize, usize)],
    i1: usize,
    i2: usize,
    i3: usize,
    i4: usize,
) -> Result<Vec<(usize, usize)>, String> {
    let primera = buscar_arista(aristas, i1, i2)
        .ok_or_else(|| format!("No hay ninguna arista de color ({i1}, {i2})."))?;
    let segunda = buscar_arista(aristas, i3, i4)
        .ok_or_else(|| format!("No hay ninguna arista de color ({i3}, {i4})."))?;
    if primera == segunda {
        return Err("Las dos aristas que se rompen deben ser distintas.".to_string());
    }
    let mut nuevas = aristas.to_vec();
    nuevas[primera] = (i1, i3);
    nuevas[segunda] = (i2, i4);
    Ok(nuevas)
}

pub fn buscar_arista(aristas: &[(usize, usize)], a: usize, b: usize) -> Option<usize> {
    aristas
        .iter()
        .position(|&arista| arista == (a, b) || arista == (b, a))
}

pub fn dos_roturas_en_genoma(
    genoma: &[Vec<i64>],
    i1: usize,
    i2: usize,
    i3: usize,
    i4: usize,
) -> Result<Genoma, String> {
    validar_genoma(genoma)?;
    let aristas = aristas_entre_vecinos(genoma);
    let nuevas = dos_roturas_en_grafo(&aristas, i1, i2, i3, i4)?;
    grafo_a_genoma(&nuevas)
}

pub fn distancia_dos_roturas(p: &[Vec<i64>], q: &[Vec<i64>]) -> Result<usize, String> {
    let ciclos = ciclos_del_grafo_de_rupturas(p, q)?;
    let bloques: usize = p.iter().map(|cromosoma| cromosoma.len()).sum();
    Ok(bloques - ciclos.len())
}

pub fn ciclos_del_grafo_de_rupturas(
    p: &[Vec<i64>],
    q: &[Vec<i64>],
) -> Result<Vec<Vec<usize>>, String> {
    let bloques = validar_mismos_bloques(p, q)?;
    let rojas = pareja_de_cada_nodo(&aristas_entre_vecinos(p), bloques);
    let azules = pareja_de_cada_nodo(&aristas_entre_vecinos(q), bloques);
    Ok(ciclos_rojo_azul(&rojas, &azules))
}

fn ciclos_rojo_azul(rojas: &[usize], azules: &[usize]) -> Vec<Vec<usize>> {
    let mut visitado = vec![false; rojas.len()];
    let mut ciclos = Vec::new();
    for inicio in 1..rojas.len() {
        if visitado[inicio] {
            continue;
        }
        let mut ciclo = Vec::new();
        let mut nodo = inicio;
        loop {
            let vecino_rojo = rojas[nodo];
            visitado[nodo] = true;
            visitado[vecino_rojo] = true;
            ciclo.push(nodo);
            ciclo.push(vecino_rojo);
            nodo = azules[vecino_rojo];
            if nodo == inicio {
                break;
            }
        }
        ciclos.push(ciclo);
    }
    ciclos
}

pub fn transformacion_mas_corta(p: &[Vec<i64>], q: &[Vec<i64>]) -> Result<Vec<DosRoturas>, String> {
    let bloques = validar_mismos_bloques(p, q)?;
    let mut rojas = pareja_de_cada_nodo(&aristas_entre_vecinos(p), bloques);
    let azules = pareja_de_cada_nodo(&aristas_entre_vecinos(q), bloques);
    let mut roturas = Vec::new();
    for nodo in 1..=2 * bloques {
        if rojas[nodo] != azules[nodo] {
            let i1 = nodo;
            let i2 = rojas[i1];
            let i3 = azules[i1];
            let i4 = rojas[i3];
            reconectar(&mut rojas, i1, i2, i3, i4);
            roturas.push(DosRoturas {
                nodos: [i1, i2, i3, i4],
                genoma: recorrer_genoma(&rojas, bloques),
            });
        }
    }
    Ok(roturas)
}

fn reconectar(rojas: &mut [usize], i1: usize, i2: usize, i3: usize, i4: usize) {
    rojas[i1] = i3;
    rojas[i3] = i1;
    rojas[i2] = i4;
    rojas[i4] = i2;
}

pub struct VentanaDeKmeros {
    pub k: usize,
    pub mascara: u64,
    pub directo: u64,
    pub reverso: u64,
    pub bases_seguidas: usize,
}

impl VentanaDeKmeros {
    pub fn nueva(k: usize) -> VentanaDeKmeros {
        let mascara = if k >= 32 {
            u64::MAX
        } else {
            (1 << (2 * k)) - 1
        };
        VentanaDeKmeros {
            k,
            mascara,
            directo: 0,
            reverso: 0,
            bases_seguidas: 0,
        }
    }

    pub fn agregar(&mut self, base: u8) -> Option<(u64, bool)> {
        let Some(valor) = indice_base(base) else {
            self.bases_seguidas = 0;
            return None;
        };
        let valor = valor as u64;
        self.directo = ((self.directo << 2) | valor) & self.mascara;
        self.reverso = (self.reverso >> 2) | ((3 - valor) << (2 * (self.k - 1)));
        self.bases_seguidas += 1;
        if self.bases_seguidas < self.k {
            return None;
        }
        Some((self.directo.min(self.reverso), self.directo <= self.reverso))
    }
}

pub struct IndiceDeKmeros {
    pub primera_y_ultima: HashMap<u64, (usize, usize)>,
    pub siguiente: Vec<usize>,
    pub directo: Vec<bool>,
}

pub fn indexar_kmeros(secuencia: &[u8], k: usize) -> IndiceDeKmeros {
    let mut indice = IndiceDeKmeros {
        primera_y_ultima: HashMap::with_capacity(secuencia.len()),
        siguiente: vec![NINGUNA; secuencia.len()],
        directo: vec![true; secuencia.len()],
    };
    let mut ventana = VentanaDeKmeros::nueva(k);
    for (fin, &base) in secuencia.iter().enumerate() {
        let Some((codigo, directo)) = ventana.agregar(base) else {
            continue;
        };
        let posicion = fin + 1 - k;
        indice.directo[posicion] = directo;
        match indice.primera_y_ultima.get_mut(&codigo) {
            Some((_, ultima)) => {
                indice.siguiente[*ultima] = posicion;
                *ultima = posicion;
            }
            None => {
                indice.primera_y_ultima.insert(codigo, (posicion, posicion));
            }
        }
    }
    indice
}

pub fn kmeros_compartidos(a: &[u8], b: &[u8], k: usize) -> Result<Vec<ParCompartido>, String> {
    if k == 0 || k > 32 {
        return Err(format!("k debe estar entre 1 y 32, y es {k}."));
    }
    let indice = indexar_kmeros(b, k);
    let mut pares = Vec::new();
    let mut ventana = VentanaDeKmeros::nueva(k);
    for (fin, &base) in a.iter().enumerate() {
        let Some((codigo, directo_a)) = ventana.agregar(base) else {
            continue;
        };
        let Some(&(primera, _)) = indice.primera_y_ultima.get(&codigo) else {
            continue;
        };
        let posicion_a = fin + 1 - k;
        let mut posicion_b = primera;
        while posicion_b != NINGUNA {
            if pares.len() == MAXIMO_PARES {
                let consejo = if k < 32 {
                    "Usa un k mayor o regiones más cortas."
                } else {
                    "Usa regiones más cortas."
                };
                return Err(format!(
                    "Hay más de {} millones de parejas de k-meros compartidos: las secuencias son demasiado repetitivas para k = {k}. {consejo}",
                    MAXIMO_PARES / 1_000_000
                ));
            }
            pares.push(ParCompartido {
                posicion_a,
                posicion_b,
                invertido: directo_a != indice.directo[posicion_b],
            });
            posicion_b = indice.siguiente[posicion_b];
        }
    }
    Ok(pares)
}

pub fn bloques_de_sintenia(
    a: &[u8],
    b: &[u8],
    k: usize,
    distancia: usize,
    minimo: usize,
) -> Result<Sintenia, String> {
    if distancia == 0 {
        return Err("La distancia para unir k-meros vecinos debe ser al menos 1.".to_string());
    }
    let pares = kmeros_compartidos(a, b, k)?;
    let bloques = agrupar_en_bloques(&pares, k, distancia, minimo);
    Ok(Sintenia { pares, bloques })
}

pub type Casilla = (usize, usize, bool);

pub fn casilla_de(par: &ParCompartido, distancia: usize) -> Casilla {
    (
        par.posicion_a / distancia,
        par.posicion_b / distancia,
        par.invertido,
    )
}

fn agrupar_en_bloques(
    pares: &[ParCompartido],
    k: usize,
    distancia: usize,
    minimo: usize,
) -> Vec<BloqueSintenico> {
    let mut casillas: HashMap<Casilla, Vec<usize>> = HashMap::new();
    for (numero, par) in pares.iter().enumerate() {
        casillas
            .entry(casilla_de(par, distancia))
            .or_default()
            .push(numero);
    }
    let mut visitadas: HashSet<Casilla> = HashSet::new();
    let mut bloques = Vec::new();
    for par in pares {
        let inicio = casilla_de(par, distancia);
        if visitadas.insert(inicio) {
            let bloque = juntar_casillas_vecinas(inicio, &casillas, &mut visitadas, pares, k);
            if bloque.kmeros >= minimo {
                bloques.push(bloque);
            }
        }
    }
    bloques.sort_by_key(|bloque| (bloque.inicio_a, bloque.inicio_b));
    bloques
}

pub fn juntar_casillas_vecinas(
    inicio: Casilla,
    casillas: &HashMap<Casilla, Vec<usize>>,
    visitadas: &mut HashSet<Casilla>,
    pares: &[ParCompartido],
    k: usize,
) -> BloqueSintenico {
    let mut bloque = BloqueSintenico {
        inicio_a: usize::MAX,
        fin_a: 0,
        inicio_b: usize::MAX,
        fin_b: 0,
        invertido: inicio.2,
        kmeros: 0,
    };
    let mut pendientes = vec![inicio];
    while let Some((fila, columna, invertido)) = pendientes.pop() {
        for &numero in &casillas[&(fila, columna, invertido)] {
            sumar_al_bloque(&mut bloque, &pares[numero], k);
        }
        for otra_fila in fila.saturating_sub(1)..=fila + 1 {
            for otra_columna in columna.saturating_sub(1)..=columna + 1 {
                let vecina = (otra_fila, otra_columna, invertido);
                if casillas.contains_key(&vecina) && visitadas.insert(vecina) {
                    pendientes.push(vecina);
                }
            }
        }
    }
    bloque
}

pub fn sumar_al_bloque(bloque: &mut BloqueSintenico, par: &ParCompartido, k: usize) {
    bloque.inicio_a = bloque.inicio_a.min(par.posicion_a);
    bloque.fin_a = bloque.fin_a.max(par.posicion_a + k - 1);
    bloque.inicio_b = bloque.inicio_b.min(par.posicion_b);
    bloque.fin_b = bloque.fin_b.max(par.posicion_b + k - 1);
    bloque.kmeros += 1;
}

pub fn permutacion_de_bloques(bloques: &[BloqueSintenico]) -> Vec<i64> {
    let mut orden_en_b: Vec<usize> = (0..bloques.len()).collect();
    orden_en_b.sort_by_key(|&indice| (bloques[indice].inicio_b, indice));
    let mut permutacion = Vec::with_capacity(bloques.len());
    for indice in orden_en_b {
        let numero = indice as i64 + 1;
        if bloques[indice].invertido {
            permutacion.push(-numero);
        } else {
            permutacion.push(numero);
        }
    }
    permutacion
}

pub fn bloque_de_cada_par(pares: &[ParCompartido], bloques: &[BloqueSintenico]) -> Vec<usize> {
    let mut orden: Vec<usize> = (0..pares.len()).collect();
    orden.sort_by_key(|&indice| pares[indice].posicion_a);
    let mut resultado = vec![0; pares.len()];
    for (indice_bloque, bloque) in bloques.iter().enumerate() {
        let primero = orden.partition_point(|&indice| pares[indice].posicion_a < bloque.inicio_a);
        for &indice in &orden[primero..] {
            let par = &pares[indice];
            if par.posicion_a > bloque.fin_a {
                break;
            }
            let dentro = resultado[indice] == 0
                && par.invertido == bloque.invertido
                && (bloque.inicio_b..=bloque.fin_b).contains(&par.posicion_b);
            if dentro {
                resultado[indice] = indice_bloque + 1;
            }
        }
    }
    resultado
}
