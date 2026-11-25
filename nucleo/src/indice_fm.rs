use crate::secuencia::complemento_inverso;
use rayon::prelude::*;

pub const CENTINELA: u8 = b'$';
pub const LETRAS_DEL_INDICE: &[u8; 6] = b"$ACGTN";
pub const SIMBOLOS: usize = 6;
pub const LETRAS_DE_UN_BYTE: usize = 257;
pub const PASO_DE_CONTROL: usize = 64;
pub const PASO_DE_MUESTRA: usize = 16;
pub const MINIMO_PARA_REPARTIR: usize = 2048;
pub const MAXIMO_DIFERENCIAS: usize = 2;
pub const MAXIMO_APARICIONES: usize = 500_000;
pub const MAXIMO_LARGO: usize = 1 << 31;
pub const VACIO: u32 = u32::MAX;
pub const SIN_LETRA: u8 = u8::MAX;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Acierto {
    pub posicion: usize,
    pub hebra: u8,
    pub diferencias: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Intervalo {
    pub hebra: u8,
    pub arriba: usize,
    pub abajo: usize,
    pub diferencias: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Rama {
    pub faltan: usize,
    pub arriba: usize,
    pub abajo: usize,
    pub diferencias: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Bloque {
    pub antes: [u32; SIMBOLOS],
    pub marcadas_antes: u32,
    pub marcadas: u64,
    pub letras: [u8; PASO_DE_CONTROL],
}

pub struct IndiceFm {
    pub filas: usize,
    pub primera_fila: [usize; SIMBOLOS + 1],
    pub bloques: Vec<Bloque>,
    pub muestras: Vec<u32>,
}

pub struct IndiceDoble {
    pub directo: IndiceFm,
    pub del_reves: IndiceFm,
}

pub fn simbolo_de(letra: u8) -> u8 {
    match letra.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 3,
        b'T' => 4,
        _ => 5,
    }
}

pub fn tipos_de_sufijo(texto: &[u32]) -> Vec<bool> {
    let n = texto.len();
    let mut es_s = vec![false; n];
    es_s[n - 1] = true;
    for i in (0..n - 1).rev() {
        es_s[i] = texto[i] < texto[i + 1] || (texto[i] == texto[i + 1] && es_s[i + 1]);
    }
    es_s
}

pub fn es_lms(es_s: &[bool], posicion: usize) -> bool {
    posicion > 0 && es_s[posicion] && !es_s[posicion - 1]
}

pub fn posiciones_lms(es_s: &[bool]) -> Vec<u32> {
    let mut lms = Vec::new();
    for posicion in 1..es_s.len() {
        if es_lms(es_s, posicion) {
            lms.push(posicion as u32);
        }
    }
    lms
}

pub fn tamanos_de_cubeta(texto: &[u32], alfabeto: usize) -> Vec<u32> {
    let mut tamanos = vec![0u32; alfabeto];
    for &letra in texto {
        tamanos[letra as usize] += 1;
    }
    tamanos
}

pub fn inicios_de_cubeta(tamanos: &[u32]) -> Vec<u32> {
    let mut inicios = Vec::with_capacity(tamanos.len());
    let mut suma = 0;
    for &tamano in tamanos {
        inicios.push(suma);
        suma += tamano;
    }
    inicios
}

pub fn finales_de_cubeta(tamanos: &[u32]) -> Vec<u32> {
    let mut finales = Vec::with_capacity(tamanos.len());
    let mut suma = 0;
    for &tamano in tamanos {
        suma += tamano;
        finales.push(suma);
    }
    finales
}

pub fn colocar_lms(texto: &[u32], tamanos: &[u32], lms: &[u32], sufijos: &mut [u32]) {
    sufijos.fill(VACIO);
    let mut finales = finales_de_cubeta(tamanos);
    for &posicion in lms.iter().rev() {
        let letra = texto[posicion as usize] as usize;
        finales[letra] -= 1;
        sufijos[finales[letra] as usize] = posicion;
    }
}

pub fn inducir_tipo_l(texto: &[u32], es_s: &[bool], tamanos: &[u32], sufijos: &mut [u32]) {
    let mut inicios = inicios_de_cubeta(tamanos);
    for i in 0..sufijos.len() {
        let posicion = sufijos[i];
        if posicion == VACIO || posicion == 0 {
            continue;
        }
        let anterior = posicion as usize - 1;
        if !es_s[anterior] {
            let letra = texto[anterior] as usize;
            sufijos[inicios[letra] as usize] = anterior as u32;
            inicios[letra] += 1;
        }
    }
}

pub fn inducir_tipo_s(texto: &[u32], es_s: &[bool], tamanos: &[u32], sufijos: &mut [u32]) {
    let mut finales = finales_de_cubeta(tamanos);
    for i in (0..sufijos.len()).rev() {
        let posicion = sufijos[i];
        if posicion == VACIO || posicion == 0 {
            continue;
        }
        let anterior = posicion as usize - 1;
        if es_s[anterior] {
            let letra = texto[anterior] as usize;
            finales[letra] -= 1;
            sufijos[finales[letra] as usize] = anterior as u32;
        }
    }
}

pub fn lms_ordenados(es_s: &[bool], sufijos: &[u32]) -> Vec<u32> {
    let mut ordenados = Vec::new();
    for &posicion in sufijos {
        if posicion != VACIO && es_lms(es_s, posicion as usize) {
            ordenados.push(posicion);
        }
    }
    ordenados
}

pub fn mismas_subcadenas_lms(texto: &[u32], es_s: &[bool], una: usize, otra: usize) -> bool {
    let mut d = 0;
    loop {
        if texto[una + d] != texto[otra + d] || es_s[una + d] != es_s[otra + d] {
            return false;
        }
        if d > 0 {
            let termina_una = es_lms(es_s, una + d);
            let termina_otra = es_lms(es_s, otra + d);
            if termina_una || termina_otra {
                return termina_una && termina_otra;
            }
        }
        d += 1;
    }
}

pub fn nombrar_lms(
    texto: &[u32],
    es_s: &[bool],
    lms: &[u32],
    ordenados: &[u32],
) -> (Vec<u32>, usize) {
    let mut nombre_en = vec![VACIO; texto.len()];
    let mut nombre = 0;
    for (lugar, &posicion) in ordenados.iter().enumerate() {
        if lugar > 0
            && !mismas_subcadenas_lms(
                texto,
                es_s,
                ordenados[lugar - 1] as usize,
                posicion as usize,
            )
        {
            nombre += 1;
        }
        nombre_en[posicion as usize] = nombre;
    }
    let mut reducido = Vec::with_capacity(lms.len());
    for &posicion in lms {
        reducido.push(nombre_en[posicion as usize]);
    }
    (reducido, nombre as usize + 1)
}

pub fn ordenar_sufijos(texto: &[u32], alfabeto: usize) -> Vec<u32> {
    let n = texto.len();
    if n == 1 {
        return vec![0];
    }
    let es_s = tipos_de_sufijo(texto);
    let tamanos = tamanos_de_cubeta(texto, alfabeto);
    let lms = posiciones_lms(&es_s);
    let mut sufijos = vec![VACIO; n];
    colocar_lms(texto, &tamanos, &lms, &mut sufijos);
    inducir_tipo_l(texto, &es_s, &tamanos, &mut sufijos);
    inducir_tipo_s(texto, &es_s, &tamanos, &mut sufijos);
    let ordenados = lms_ordenados(&es_s, &sufijos);
    let (reducido, nombres) = nombrar_lms(texto, &es_s, &lms, &ordenados);
    let lms_en_orden: Vec<u32> = if nombres < lms.len() {
        let orden_reducido = ordenar_sufijos(&reducido, nombres);
        orden_reducido.iter().map(|&i| lms[i as usize]).collect()
    } else {
        ordenados
    };
    colocar_lms(texto, &tamanos, &lms_en_orden, &mut sufijos);
    inducir_tipo_l(texto, &es_s, &tamanos, &mut sufijos);
    inducir_tipo_s(texto, &es_s, &tamanos, &mut sufijos);
    sufijos
}

pub fn revisar_largo(largo: usize) -> Result<(), String> {
    if largo >= MAXIMO_LARGO {
        return Err(
            "El texto es demasiado largo: aquí se admiten hasta dos mil millones de letras."
                .to_string(),
        );
    }
    Ok(())
}

pub fn arreglo_de_sufijos(texto: &[u8]) -> Result<Vec<u32>, String> {
    revisar_largo(texto.len())?;
    let mut simbolos: Vec<u32> = Vec::with_capacity(texto.len() + 1);
    for &letra in texto {
        simbolos.push(letra as u32 + 1);
    }
    simbolos.push(0);
    Ok(ordenar_sufijos(&simbolos, LETRAS_DE_UN_BYTE))
}

pub fn ultima_columna(texto: &[u8], sufijos: &[u32]) -> Vec<u8> {
    let mut columna = Vec::with_capacity(sufijos.len());
    for &inicio in sufijos {
        if inicio == 0 {
            columna.push(CENTINELA);
        } else {
            columna.push(texto[inicio as usize - 1]);
        }
    }
    columna
}

pub fn transformada(texto: &[u8]) -> Result<Vec<u8>, String> {
    if texto.contains(&CENTINELA) {
        return Err("El texto no puede llevar «$»: ese símbolo marca el final.".to_string());
    }
    let sufijos = arreglo_de_sufijos(texto)?;
    Ok(ultima_columna(texto, &sufijos))
}

pub fn orden_de_letra(letra: u8) -> usize {
    if letra == CENTINELA {
        0
    } else {
        letra as usize + 1
    }
}

pub fn primeras_filas_de_letras(bwt: &[u8]) -> Vec<usize> {
    let mut cuentas = vec![0usize; LETRAS_DE_UN_BYTE];
    for &letra in bwt {
        cuentas[orden_de_letra(letra)] += 1;
    }
    let mut primera_fila = Vec::with_capacity(LETRAS_DE_UN_BYTE);
    let mut suma = 0;
    for cuenta in cuentas {
        primera_fila.push(suma);
        suma += cuenta;
    }
    primera_fila
}

pub fn rangos(bwt: &[u8]) -> Vec<u32> {
    let mut vistas = vec![0u32; LETRAS_DE_UN_BYTE];
    let mut rango = Vec::with_capacity(bwt.len());
    for &letra in bwt {
        let orden = orden_de_letra(letra);
        rango.push(vistas[orden]);
        vistas[orden] += 1;
    }
    rango
}

pub fn inversa(bwt: &[u8]) -> Result<Vec<u8>, String> {
    let centinelas = bwt.iter().filter(|&&letra| letra == CENTINELA).count();
    if centinelas != 1 {
        return Err(format!(
            "La BWT debe tener exactamente un «$», y esta tiene {centinelas}."
        ));
    }
    let primera_fila = primeras_filas_de_letras(bwt);
    let rango = rangos(bwt);
    let mut texto = vec![0u8; bwt.len() - 1];
    let mut fila = 0;
    for posicion in (0..texto.len()).rev() {
        let letra = bwt[fila];
        texto[posicion] = letra;
        fila = primera_fila[orden_de_letra(letra)] + rango[fila] as usize;
    }
    Ok(texto)
}

pub fn transformada_e_inversa(texto: &[u8]) -> Result<(Vec<u8>, Vec<u8>), String> {
    let bwt = transformada(texto)?;
    let recuperado = inversa(&bwt)?;
    Ok((bwt, recuperado))
}

pub fn rachas(secuencia: &[u8]) -> usize {
    if secuencia.is_empty() {
        return 0;
    }
    let mut cuenta = 1;
    for par in secuencia.windows(2) {
        if par[0] != par[1] {
            cuenta += 1;
        }
    }
    cuenta
}

fn primeras_filas(simbolos: &[u32]) -> [usize; SIMBOLOS + 1] {
    let mut cuentas = [0usize; SIMBOLOS];
    for &simbolo in simbolos {
        cuentas[simbolo as usize] += 1;
    }
    let mut primera_fila = [0usize; SIMBOLOS + 1];
    for simbolo in 0..SIMBOLOS {
        primera_fila[simbolo + 1] = primera_fila[simbolo] + cuentas[simbolo];
    }
    primera_fila
}

fn letra_l(simbolos: &[u32], inicio: u32) -> u8 {
    if inicio == 0 {
        0
    } else {
        simbolos[inicio as usize - 1] as u8
    }
}

fn bloque_nuevo(antes: [u32; SIMBOLOS], marcadas_antes: u32) -> Bloque {
    Bloque {
        antes,
        marcadas_antes,
        marcadas: 0,
        letras: [SIN_LETRA; PASO_DE_CONTROL],
    }
}

fn armar_bloques(simbolos: &[u32], sufijos: &[u32]) -> Vec<Bloque> {
    let mut bloques = Vec::with_capacity(sufijos.len() / PASO_DE_CONTROL + 1);
    let mut antes = [0u32; SIMBOLOS];
    let mut marcadas_antes = 0;
    for (fila, &inicio) in sufijos.iter().enumerate() {
        if fila.is_multiple_of(PASO_DE_CONTROL) {
            bloques.push(bloque_nuevo(antes, marcadas_antes));
        }
        let simbolo = letra_l(simbolos, inicio);
        let lugar = fila % PASO_DE_CONTROL;
        let bloque = &mut bloques[fila / PASO_DE_CONTROL];
        bloque.letras[lugar] = simbolo;
        antes[simbolo as usize] += 1;
        if (inicio as usize).is_multiple_of(PASO_DE_MUESTRA) {
            bloque.marcadas |= 1 << lugar;
            marcadas_antes += 1;
        }
    }
    if sufijos.len().is_multiple_of(PASO_DE_CONTROL) {
        bloques.push(bloque_nuevo(antes, marcadas_antes));
    }
    bloques
}

fn tomar_muestras(sufijos: &[u32]) -> Vec<u32> {
    let mut muestras = Vec::with_capacity(sufijos.len() / PASO_DE_MUESTRA + 1);
    for &inicio in sufijos {
        if (inicio as usize).is_multiple_of(PASO_DE_MUESTRA) {
            muestras.push(inicio);
        }
    }
    muestras
}

impl IndiceFm {
    pub fn construir(texto: &[u8]) -> Result<IndiceFm, String> {
        if texto.is_empty() {
            return Err("El texto está vacío: no hay nada que indexar.".to_string());
        }
        revisar_largo(texto.len())?;
        let mut simbolos: Vec<u32> = Vec::with_capacity(texto.len() + 1);
        for &letra in texto {
            simbolos.push(simbolo_de(letra) as u32);
        }
        simbolos.push(0);
        let sufijos = ordenar_sufijos(&simbolos, SIMBOLOS);
        Ok(IndiceFm {
            filas: sufijos.len(),
            primera_fila: primeras_filas(&simbolos),
            bloques: armar_bloques(&simbolos, &sufijos),
            muestras: tomar_muestras(&sufijos),
        })
    }

    pub fn simbolo_de_l(&self, fila: usize) -> u8 {
        self.bloques[fila / PASO_DE_CONTROL].letras[fila % PASO_DE_CONTROL]
    }

    pub fn ocurrencias(&self, simbolo: u8, fila: usize) -> usize {
        let bloque = &self.bloques[fila / PASO_DE_CONTROL];
        let dentro_del_bloque = bloque.letras[..fila % PASO_DE_CONTROL]
            .iter()
            .filter(|&&letra| letra == simbolo)
            .count();
        bloque.antes[simbolo as usize] as usize + dentro_del_bloque
    }

    pub fn ocurrencias_de_todas(&self, fila: usize) -> [usize; SIMBOLOS] {
        let bloque = &self.bloques[fila / PASO_DE_CONTROL];
        let mut cuentas = bloque.antes.map(|cuenta| cuenta as usize);
        for &letra in &bloque.letras[..fila % PASO_DE_CONTROL] {
            cuentas[letra as usize] += 1;
        }
        cuentas
    }

    pub fn paso_atras(&self, simbolo: u8, arriba: usize, abajo: usize) -> (usize, usize) {
        let primera = self.primera_fila[simbolo as usize];
        (
            primera + self.ocurrencias(simbolo, arriba),
            primera + self.ocurrencias(simbolo, abajo),
        )
    }

    pub fn intervalo(&self, patron: &[u8]) -> (usize, usize) {
        let mut arriba = 0;
        let mut abajo = self.filas;
        for &letra in patron.iter().rev() {
            (arriba, abajo) = self.paso_atras(simbolo_de(letra), arriba, abajo);
            if arriba >= abajo {
                return (0, 0);
            }
        }
        (arriba, abajo)
    }

    pub fn contar(&self, patron: &[u8]) -> usize {
        let (arriba, abajo) = self.intervalo(patron);
        abajo - arriba
    }

    pub fn simbolo_de_f(&self, fila: usize) -> u8 {
        let mut simbolo = 0;
        while simbolo + 1 < SIMBOLOS && self.primera_fila[simbolo + 1] <= fila {
            simbolo += 1;
        }
        simbolo as u8
    }

    pub fn lf(&self, fila: usize) -> usize {
        let simbolo = self.simbolo_de_l(fila);
        self.primera_fila[simbolo as usize] + self.ocurrencias(simbolo, fila)
    }

    pub fn esta_marcada(&self, fila: usize) -> bool {
        let bloque = &self.bloques[fila / PASO_DE_CONTROL];
        (bloque.marcadas >> (fila % PASO_DE_CONTROL)) & 1 == 1
    }

    pub fn marcadas_hasta(&self, fila: usize) -> usize {
        let bloque = &self.bloques[fila / PASO_DE_CONTROL];
        let anteriores = bloque.marcadas & ((1u64 << (fila % PASO_DE_CONTROL)) - 1);
        bloque.marcadas_antes as usize + anteriores.count_ones() as usize
    }

    pub fn posicion_de_fila(&self, fila: usize) -> usize {
        let mut actual = fila;
        let mut pasos = 0;
        while !self.esta_marcada(actual) {
            actual = self.lf(actual);
            pasos += 1;
        }
        self.muestras[self.marcadas_hasta(actual)] as usize + pasos
    }

    pub fn buscar_con_diferencias(
        &self,
        patron: &[u8],
        d: usize,
        hebra: u8,
        cotas: &[usize],
    ) -> Vec<Intervalo> {
        let mut hallados = Vec::new();
        let mut pendientes = vec![Rama {
            faltan: patron.len(),
            arriba: 0,
            abajo: self.filas,
            diferencias: 0,
        }];
        while let Some(rama) = pendientes.pop() {
            if rama.faltan == 0 {
                hallados.push(Intervalo {
                    hebra,
                    arriba: rama.arriba,
                    abajo: rama.abajo,
                    diferencias: rama.diferencias,
                });
                continue;
            }
            let esperado = simbolo_de(patron[rama.faltan - 1]) as usize;
            let cota_de_lo_que_falta = if rama.faltan >= 2 {
                cotas[rama.faltan - 2]
            } else {
                0
            };
            let antes_de_arriba = self.ocurrencias_de_todas(rama.arriba);
            let antes_de_abajo = self.ocurrencias_de_todas(rama.abajo);
            let otros = (1..SIMBOLOS).rev().filter(|&simbolo| simbolo != esperado);
            for simbolo in otros.chain([esperado]) {
                let diferencias = rama.diferencias + usize::from(simbolo != esperado);
                let arriba = self.primera_fila[simbolo] + antes_de_arriba[simbolo];
                let abajo = self.primera_fila[simbolo] + antes_de_abajo[simbolo];
                if diferencias + cota_de_lo_que_falta <= d && arriba < abajo {
                    pendientes.push(Rama {
                        faltan: rama.faltan - 1,
                        arriba,
                        abajo,
                        diferencias,
                    });
                }
            }
        }
        hallados
    }
}

impl IndiceDoble {
    pub fn construir(texto: &[u8]) -> Result<IndiceDoble, String> {
        let al_reves: Vec<u8> = texto.iter().rev().copied().collect();
        let (directo, del_reves) = rayon::join(
            || IndiceFm::construir(texto),
            || IndiceFm::construir(&al_reves),
        );
        Ok(IndiceDoble {
            directo: directo?,
            del_reves: del_reves?,
        })
    }

    pub fn cotas_de_diferencias(&self, patron: &[u8]) -> Vec<usize> {
        let mut cotas = Vec::with_capacity(patron.len());
        let mut cota = 0;
        let mut arriba = 0;
        let mut abajo = self.del_reves.filas;
        for &letra in patron {
            (arriba, abajo) = self.del_reves.paso_atras(simbolo_de(letra), arriba, abajo);
            if arriba >= abajo {
                cota += 1;
                arriba = 0;
                abajo = self.del_reves.filas;
            }
            cotas.push(cota);
        }
        cotas
    }
}

pub fn buscar_en_una_hebra(
    indice: &IndiceDoble,
    patron: &[u8],
    d: usize,
    hebra: u8,
) -> Vec<Intervalo> {
    if d > 0 {
        let cotas = indice.cotas_de_diferencias(patron);
        return indice
            .directo
            .buscar_con_diferencias(patron, d, hebra, &cotas);
    }
    let (arriba, abajo) = indice.directo.intervalo(patron);
    if arriba < abajo {
        vec![Intervalo {
            hebra,
            arriba,
            abajo,
            diferencias: 0,
        }]
    } else {
        Vec::new()
    }
}

pub fn buscar_en_ambas_hebras(indice: &IndiceDoble, patron: &[u8], d: usize) -> Vec<Intervalo> {
    let mut intervalos = buscar_en_una_hebra(indice, patron, d, b'+');
    let inverso = complemento_inverso(patron);
    if inverso != patron {
        intervalos.extend(buscar_en_una_hebra(indice, &inverso, d, b'-'));
    }
    intervalos
}

pub fn ubicar(indice: &IndiceFm, intervalos: &[Intervalo]) -> Vec<Acierto> {
    let mut filas = Vec::new();
    for intervalo in intervalos {
        for fila in intervalo.arriba..intervalo.abajo {
            filas.push((fila, intervalo));
        }
    }
    let ubicar_fila = |&(fila, intervalo): &(usize, &Intervalo)| Acierto {
        posicion: indice.posicion_de_fila(fila),
        hebra: intervalo.hebra,
        diferencias: intervalo.diferencias,
    };
    let mut aciertos: Vec<Acierto> = if filas.len() < MINIMO_PARA_REPARTIR {
        filas.iter().map(ubicar_fila).collect()
    } else {
        filas.par_iter().map(ubicar_fila).collect()
    };
    aciertos.sort_unstable();
    aciertos
}

pub fn revisar_patrones(patrones: &[Vec<u8>], d: usize) -> Result<(), String> {
    if d > MAXIMO_DIFERENCIAS {
        return Err(format!(
            "Se admiten hasta {MAXIMO_DIFERENCIAS} diferencias, y pediste {d}."
        ));
    }
    if patrones.is_empty() {
        return Err("No hay patrones que buscar.".to_string());
    }
    for patron in patrones {
        if let Some(letra) = patron
            .iter()
            .find(|letra| !matches!(letra, b'A' | b'C' | b'G' | b'T'))
        {
            return Err(format!(
                "Los patrones solo admiten A, C, G y T, pero uno contiene «{}».",
                *letra as char
            ));
        }
        if patron.len() <= d {
            return Err(format!(
                "Con {d} diferencias un patrón de {} letras coincide en cualquier lugar: usa patrones de más de {d} letras.",
                patron.len()
            ));
        }
    }
    Ok(())
}

pub fn buscar_patrones(
    indice: &IndiceDoble,
    patrones: &[Vec<u8>],
    d: usize,
) -> Result<Vec<Vec<Acierto>>, String> {
    revisar_patrones(patrones, d)?;
    let mut intervalos_de_cada = Vec::with_capacity(patrones.len());
    let mut total = 0;
    for patron in patrones {
        let intervalos = buscar_en_ambas_hebras(indice, patron, d);
        for intervalo in &intervalos {
            total += intervalo.abajo - intervalo.arriba;
        }
        if total > MAXIMO_APARICIONES {
            return Err(
                "Los patrones aparecen más de medio millón de veces en total; usa patrones más largos o menos diferencias."
                    .to_string(),
            );
        }
        intervalos_de_cada.push(intervalos);
    }
    let mut resultado = Vec::with_capacity(patrones.len());
    for intervalos in &intervalos_de_cada {
        resultado.push(ubicar(&indice.directo, intervalos));
    }
    Ok(resultado)
}
