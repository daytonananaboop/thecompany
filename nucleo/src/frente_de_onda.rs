use rayon::prelude::*;
use std::ops::{Index, IndexMut};

pub const BITS_POR_PALABRA: usize = 64;
const BIT_ALTO: u64 = 1 << 63;
const NULO: i32 = i32::MIN / 2;
const LARGO_MAXIMO: usize = 20_000_000;
const CASTIGO_MAXIMO: usize = 1000;
const LIMITE_DE_CELDAS: usize = 800_000_000;
const LIMITE_DE_MEMORIA: usize = 60_000_000;
const ANCHO_EN_PARALELO: usize = 16_384;
const TROZO_EN_PARALELO: usize = 4096;
const CELDAS_POR_PUNTAJE: usize = 24;
const PUNTAJES_ANTES_DE_SOLTAR: usize = 4096;

static NINGUNO: Option<Frente> = None;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Castigos {
    pub diferencia: usize,
    pub apertura: usize,
    pub extension: usize,
}

pub const CASTIGOS_DEL_ARTICULO: Castigos = Castigos {
    diferencia: 4,
    apertura: 6,
    extension: 2,
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlineamientoWfa {
    pub penalizacion: usize,
    pub operaciones: Option<Vec<u8>>,
    pub frentes: usize,
    pub celdas: usize,
}

pub fn alinear_wfa(v: &[u8], w: &[u8], castigos: Castigos) -> Result<AlineamientoWfa, String> {
    revisar_entrada(v, w, castigos)?;
    let mut ondas = FrentesDeOnda::nuevos(v, w, castigos);
    ondas.extender();
    while !ondas.llego_al_final() {
        ondas.siguiente()?;
        ondas.extender();
    }
    Ok(AlineamientoWfa {
        penalizacion: ondas.puntaje,
        operaciones: ondas.rastrear(),
        frentes: ondas.frentes,
        celdas: ondas.celdas,
    })
}

pub fn revisar_entrada(v: &[u8], w: &[u8], castigos: Castigos) -> Result<(), String> {
    if castigos.diferencia == 0 || castigos.extension == 0 {
        return Err(
            "El castigo por letra distinta y el de alargar un hueco deben ser al menos 1: WFA avanza de puntaje en puntaje y un castigo 0 no lo dejaría avanzar.".to_string(),
        );
    }
    if castigos.diferencia > CASTIGO_MAXIMO
        || castigos.apertura > CASTIGO_MAXIMO
        || castigos.extension > CASTIGO_MAXIMO
    {
        return Err(format!(
            "Cada castigo debe estar entre 0 y {CASTIGO_MAXIMO}."
        ));
    }
    if v.len() > LARGO_MAXIMO || w.len() > LARGO_MAXIMO {
        return Err("WFA acepta secuencias de hasta 20 millones de letras.".to_string());
    }
    let resta = v.len().abs_diff(w.len());
    if resta.saturating_mul(resta) / 2 > LIMITE_DE_CELDAS {
        return Err(
            "Las secuencias difieren demasiado en largo: solo el hueco que las iguala ya obligaría a WFA a recorrer cientos de millones de celdas. Compara tramos de largo parecido o usa la distancia de edición con vectores de bits (Myers).".to_string(),
        );
    }
    Ok(())
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Frente {
    pub bajo: i32,
    pub alto: i32,
    pub desplazamientos: Vec<i32>,
}

impl Frente {
    pub fn en(&self, diagonal: i32) -> i32 {
        Vista::de(Some(self)).en(diagonal)
    }

    pub fn ancho(&self) -> usize {
        self.desplazamientos.len()
    }

    fn nulo(bajo: i32, alto: i32) -> Frente {
        Frente {
            bajo,
            alto,
            desplazamientos: vec![NULO; (alto - bajo + 1) as usize],
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Historial {
    primero: usize,
    frentes: Vec<Option<Frente>>,
}

impl Historial {
    fn con(inicial: Option<Frente>) -> Historial {
        Historial {
            primero: 0,
            frentes: vec![inicial],
        }
    }

    fn agregar(&mut self, frente: Option<Frente>) {
        self.frentes.push(frente);
    }

    fn soltar_hasta(&mut self, puntaje: usize) {
        if puntaje > self.primero {
            self.frentes.drain(..puntaje - self.primero);
            self.primero = puntaje;
        }
    }
}

impl Index<usize> for Historial {
    type Output = Option<Frente>;

    fn index(&self, puntaje: usize) -> &Option<Frente> {
        match puntaje.checked_sub(self.primero) {
            Some(posicion) => &self.frentes[posicion],
            None => &NINGUNO,
        }
    }
}

impl IndexMut<usize> for Historial {
    fn index_mut(&mut self, puntaje: usize) -> &mut Option<Frente> {
        &mut self.frentes[puntaje - self.primero]
    }
}

#[derive(Clone, Copy)]
struct Vista<'f> {
    bajo: i32,
    desplazamientos: &'f [i32],
}

impl Vista<'_> {
    fn de(frente: Option<&Frente>) -> Vista<'_> {
        match frente {
            Some(frente) => Vista {
                bajo: frente.bajo,
                desplazamientos: &frente.desplazamientos,
            },
            None => Vista {
                bajo: 0,
                desplazamientos: &[],
            },
        }
    }

    fn en(&self, diagonal: i32) -> i32 {
        let posicion = (diagonal - self.bajo) as usize;
        self.desplazamientos.get(posicion).copied().unwrap_or(NULO)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Origen {
    Inicio,
    Diferencia,
    Insercion,
    Borrado,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Componente {
    Parejas,
    Inserciones,
    Borrados,
}

pub struct FrentesDeOnda<'a> {
    v: &'a [u8],
    w: &'a [u8],
    largo_v: i32,
    largo_w: i32,
    castigos: Castigos,
    pub puntaje: usize,
    pub parejas: Historial,
    pub inserciones: Historial,
    pub borrados: Historial,
    pub frentes: usize,
    pub celdas: usize,
    pub solo_puntaje: bool,
    guardadas: usize,
    olvidados: usize,
}

impl<'a> FrentesDeOnda<'a> {
    pub fn nuevos(v: &'a [u8], w: &'a [u8], castigos: Castigos) -> FrentesDeOnda<'a> {
        let origen = Frente {
            bajo: 0,
            alto: 0,
            desplazamientos: vec![0],
        };
        FrentesDeOnda {
            v,
            w,
            largo_v: v.len() as i32,
            largo_w: w.len() as i32,
            castigos,
            puntaje: 0,
            parejas: Historial::con(Some(origen)),
            inserciones: Historial::con(None),
            borrados: Historial::con(None),
            frentes: 1,
            celdas: 1,
            solo_puntaje: false,
            guardadas: 1,
            olvidados: 0,
        }
    }

    pub fn diagonal_final(&self) -> i32 {
        self.largo_w - self.largo_v
    }

    pub fn extender(&mut self) {
        let (v, w) = (self.v, self.w);
        let puntaje = self.puntaje;
        let Some(frente) = self.parejas[puntaje].as_mut() else {
            return;
        };
        let bajo = frente.bajo;
        if frente.ancho() < ANCHO_EN_PARALELO {
            extender_trozo(v, w, bajo, &mut frente.desplazamientos);
            return;
        }
        frente
            .desplazamientos
            .par_chunks_mut(TROZO_EN_PARALELO)
            .enumerate()
            .for_each(|(numero, trozo)| {
                let primera = bajo + (numero * TROZO_EN_PARALELO) as i32;
                extender_trozo(v, w, primera, trozo);
            });
    }

    pub fn llego_al_final(&self) -> bool {
        match &self.parejas[self.puntaje] {
            Some(frente) => frente.en(self.diagonal_final()) == self.largo_w,
            None => false,
        }
    }

    pub fn siguiente(&mut self) -> Result<(), String> {
        self.puntaje += 1;
        let Castigos {
            diferencia,
            apertura,
            extension,
        } = self.castigos;
        let por_diferencia = self.anterior(Componente::Parejas, diferencia);
        let por_apertura = self.anterior(Componente::Parejas, apertura + extension);
        let insercion_previa = self.anterior(Componente::Inserciones, extension);
        let borrado_previo = self.anterior(Componente::Borrados, extension);
        let inserciones = self.calcular_inserciones(por_apertura, insercion_previa);
        let borrados = self.calcular_borrados(por_apertura, borrado_previo);
        let parejas =
            self.calcular_parejas(por_diferencia, inserciones.as_ref(), borrados.as_ref());
        self.guardar(parejas, inserciones, borrados)
    }

    fn anterior(&self, componente: Componente, atras: usize) -> Option<&Frente> {
        self.en(componente, self.puntaje, atras)
    }

    fn en(&self, componente: Componente, puntaje: usize, atras: usize) -> Option<&Frente> {
        if atras > puntaje {
            return None;
        }
        let lista = match componente {
            Componente::Parejas => &self.parejas,
            Componente::Inserciones => &self.inserciones,
            Componente::Borrados => &self.borrados,
        };
        lista[puntaje - atras].as_ref()
    }

    fn calcular_inserciones(
        &self,
        por_apertura: Option<&Frente>,
        previa: Option<&Frente>,
    ) -> Option<Frente> {
        let (bajo, alto) = rango(&[por_apertura, previa])?;
        let (bajo, alto) = self.recortar(bajo + 1, alto + 1)?;
        let mut inserciones = Frente::nulo(bajo, alto);
        self.tomar_el_mayor(&mut inserciones, por_apertura, 1, 1);
        self.tomar_el_mayor(&mut inserciones, previa, 1, 1);
        Some(inserciones)
    }

    fn calcular_borrados(
        &self,
        por_apertura: Option<&Frente>,
        previo: Option<&Frente>,
    ) -> Option<Frente> {
        let (bajo, alto) = rango(&[por_apertura, previo])?;
        let (bajo, alto) = self.recortar(bajo - 1, alto - 1)?;
        let mut borrados = Frente::nulo(bajo, alto);
        self.tomar_el_mayor(&mut borrados, por_apertura, -1, 0);
        self.tomar_el_mayor(&mut borrados, previo, -1, 0);
        Some(borrados)
    }

    fn calcular_parejas(
        &self,
        por_diferencia: Option<&Frente>,
        inserciones: Option<&Frente>,
        borrados: Option<&Frente>,
    ) -> Option<Frente> {
        let (bajo, alto) = rango(&[por_diferencia, inserciones, borrados])?;
        let mut parejas = Frente::nulo(bajo, alto);
        self.tomar_el_mayor(&mut parejas, por_diferencia, 0, 1);
        self.tomar_el_mayor(&mut parejas, inserciones, 0, 0);
        self.tomar_el_mayor(&mut parejas, borrados, 0, 0);
        Some(parejas)
    }

    fn tomar_el_mayor(
        &self,
        destino: &mut Frente,
        fuente: Option<&Frente>,
        salto_de_diagonal: i32,
        letras_de_w: i32,
    ) {
        let Some(fuente) = fuente else {
            return;
        };
        let desde = destino.bajo.max(fuente.bajo + salto_de_diagonal);
        let hasta = destino.alto.min(fuente.alto + salto_de_diagonal);
        if desde > hasta {
            return;
        }
        let nuevos = &mut destino.desplazamientos
            [(desde - destino.bajo) as usize..=(hasta - destino.bajo) as usize];
        let primera_vieja = desde - salto_de_diagonal - fuente.bajo;
        let viejos = &fuente.desplazamientos[primera_vieja as usize..][..nuevos.len()];
        for (posicion, (nuevo, &viejo)) in nuevos.iter_mut().zip(viejos).enumerate() {
            let diagonal = desde + posicion as i32;
            *nuevo = (*nuevo).max(self.validar(diagonal, viejo + letras_de_w));
        }
    }

    fn recortar(&self, bajo: i32, alto: i32) -> Option<(i32, i32)> {
        let bajo = bajo.max(-self.largo_v);
        let alto = alto.min(self.largo_w);
        if bajo > alto {
            None
        } else {
            Some((bajo, alto))
        }
    }

    fn validar(&self, diagonal: i32, columna: i32) -> i32 {
        let fila = columna - diagonal;
        let dentro = columna >= 0 && fila >= 0 && columna <= self.largo_w && fila <= self.largo_v;
        if dentro { columna } else { NULO }
    }

    fn guardar(
        &mut self,
        parejas: Option<Frente>,
        inserciones: Option<Frente>,
        borrados: Option<Frente>,
    ) -> Result<(), String> {
        let nuevas = ancho(&parejas) + ancho(&inserciones) + ancho(&borrados);
        self.celdas += nuevas;
        if self.celdas > LIMITE_DE_CELDAS {
            return Err(
                "WFA llegó a su límite de cientos de millones de celdas sin terminar: las secuencias son demasiado distintas para este método. Usa la distancia de edición con vectores de bits (Myers).".to_string(),
            );
        }
        if parejas.is_some() {
            self.frentes += 1;
        }
        self.parejas.agregar(parejas);
        self.inserciones.agregar(inserciones);
        self.borrados.agregar(borrados);
        self.guardadas += nuevas + CELDAS_POR_PUNTAJE;
        if self.guardadas > LIMITE_DE_MEMORIA {
            self.solo_puntaje = true;
        }
        if self.solo_puntaje {
            self.olvidar_los_viejos();
        }
        Ok(())
    }

    fn olvidar_los_viejos(&mut self) {
        let Castigos {
            diferencia,
            apertura,
            extension,
        } = self.castigos;
        let alcance = diferencia.max(apertura + extension);
        if self.olvidados + alcance + PUNTAJES_ANTES_DE_SOLTAR < self.puntaje {
            self.olvidados = self.puntaje - alcance;
            self.parejas.soltar_hasta(self.olvidados);
            self.inserciones.soltar_hasta(self.olvidados);
            self.borrados.soltar_hasta(self.olvidados);
        }
    }

    pub fn origen_de_pareja(&self, puntaje: usize, diagonal: i32) -> (i32, Origen) {
        if puntaje == 0 {
            return (0, Origen::Inicio);
        }
        let cambio = Vista::de(self.en(Componente::Parejas, puntaje, self.castigos.diferencia));
        let cambiada = self.validar(diagonal, cambio.en(diagonal) + 1);
        let insercion = Vista::de(self.en(Componente::Inserciones, puntaje, 0)).en(diagonal);
        let borrado = Vista::de(self.en(Componente::Borrados, puntaje, 0)).en(diagonal);
        let mejor = cambiada.max(insercion).max(borrado);
        if mejor == cambiada {
            (mejor, Origen::Diferencia)
        } else if mejor == insercion {
            (mejor, Origen::Insercion)
        } else {
            (mejor, Origen::Borrado)
        }
    }

    pub fn rastrear(&self) -> Option<Vec<u8>> {
        if self.solo_puntaje {
            return None;
        }
        let abrir = self.castigos.apertura + self.castigos.extension;
        let alargar = self.castigos.extension;
        let mut operaciones = Vec::new();
        let mut puntaje = self.puntaje;
        let mut diagonal = self.diagonal_final();
        let mut columna = self.largo_w;
        let mut componente = Componente::Parejas;
        loop {
            match componente {
                Componente::Parejas => {
                    let (antes, origen) = self.origen_de_pareja(puntaje, diagonal);
                    let iguales = (columna - antes) as usize;
                    operaciones.extend(std::iter::repeat_n(b'=', iguales));
                    columna = antes;
                    match origen {
                        Origen::Inicio => break,
                        Origen::Diferencia => {
                            operaciones.push(b'X');
                            puntaje -= self.castigos.diferencia;
                            columna -= 1;
                        }
                        Origen::Insercion => componente = Componente::Inserciones,
                        Origen::Borrado => componente = Componente::Borrados,
                    }
                }
                Componente::Inserciones => {
                    operaciones.push(b'I');
                    let desde_pareja = Vista::de(self.en(Componente::Parejas, puntaje, abrir));
                    if self.validar(diagonal, desde_pareja.en(diagonal - 1) + 1) == columna {
                        puntaje -= abrir;
                        componente = Componente::Parejas;
                    } else {
                        puntaje -= alargar;
                    }
                    diagonal -= 1;
                    columna -= 1;
                }
                Componente::Borrados => {
                    operaciones.push(b'D');
                    let desde_pareja = Vista::de(self.en(Componente::Parejas, puntaje, abrir));
                    if self.validar(diagonal, desde_pareja.en(diagonal + 1)) == columna {
                        puntaje -= abrir;
                        componente = Componente::Parejas;
                    } else {
                        puntaje -= alargar;
                    }
                    diagonal += 1;
                }
            }
        }
        operaciones.reverse();
        Some(operaciones)
    }
}

fn extender_trozo(v: &[u8], w: &[u8], primera: i32, desplazamientos: &mut [i32]) {
    for (posicion, columna) in desplazamientos.iter_mut().enumerate() {
        if *columna < 0 {
            continue;
        }
        let diagonal = primera + posicion as i32;
        let j = *columna as usize;
        let i = (*columna - diagonal) as usize;
        *columna += letras_iguales(&v[i..], &w[j..]) as i32;
    }
}

pub fn letras_iguales(a: &[u8], b: &[u8]) -> usize {
    let largo = a.len().min(b.len());
    let mut iguales = 0;
    while iguales + 8 <= largo {
        let mut ocho_de_a = [0u8; 8];
        let mut ocho_de_b = [0u8; 8];
        ocho_de_a.copy_from_slice(&a[iguales..iguales + 8]);
        ocho_de_b.copy_from_slice(&b[iguales..iguales + 8]);
        let distintos = u64::from_le_bytes(ocho_de_a) ^ u64::from_le_bytes(ocho_de_b);
        if distintos != 0 {
            return iguales + (distintos.trailing_zeros() / 8) as usize;
        }
        iguales += 8;
    }
    while iguales < largo && a[iguales] == b[iguales] {
        iguales += 1;
    }
    iguales
}

fn rango(fuentes: &[Option<&Frente>]) -> Option<(i32, i32)> {
    let presentes = fuentes.iter().flatten();
    let bajo = presentes.clone().map(|frente| frente.bajo).min()?;
    let alto = presentes.map(|frente| frente.alto).max()?;
    Some((bajo, alto))
}

fn ancho(frente: &Option<Frente>) -> usize {
    frente.as_ref().map(Frente::ancho).unwrap_or(0)
}

pub fn cigar(operaciones: &[u8]) -> String {
    let mut texto = String::new();
    let mut inicio = 0;
    while inicio < operaciones.len() {
        let letra = operaciones[inicio];
        let mut fin = inicio;
        while fin < operaciones.len() && operaciones[fin] == letra {
            fin += 1;
        }
        texto.push_str(&(fin - inicio).to_string());
        texto.push(letra as char);
        inicio = fin;
    }
    texto
}

pub fn filas_alineadas(v: &[u8], w: &[u8], operaciones: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let mut fila_v = Vec::with_capacity(operaciones.len());
    let mut fila_w = Vec::with_capacity(operaciones.len());
    let (mut i, mut j) = (0, 0);
    for &operacion in operaciones {
        match operacion {
            b'I' => {
                fila_v.push(b'-');
                fila_w.push(w[j]);
                j += 1;
            }
            b'D' => {
                fila_v.push(v[i]);
                fila_w.push(b'-');
                i += 1;
            }
            _ => {
                fila_v.push(v[i]);
                fila_w.push(w[j]);
                i += 1;
                j += 1;
            }
        }
    }
    (fila_v, fila_w)
}

pub fn penalizacion_de(operaciones: &[u8], castigos: Castigos) -> usize {
    let mut total = 0;
    let mut anterior = b'=';
    for &operacion in operaciones {
        match operacion {
            b'X' => total += castigos.diferencia,
            b'I' | b'D' if operacion != anterior => total += castigos.apertura + castigos.extension,
            b'I' | b'D' => total += castigos.extension,
            _ => {}
        }
        anterior = operacion;
    }
    total
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DistanciaMyers {
    pub distancia: usize,
    pub patron_es_v: bool,
    pub palabras_por_columna: usize,
    pub pasadas: Vec<Pasada>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Pasada {
    pub limite: usize,
    pub distancia: Option<usize>,
    pub palabras_calculadas: usize,
}

pub fn distancia_myers(v: &[u8], w: &[u8]) -> DistanciaMyers {
    let patron_es_v = v.len() <= w.len();
    let (patron, texto) = if patron_es_v { (v, w) } else { (w, v) };
    let tabla = TablaDeIguales::para(patron);
    let mut pasadas = Vec::new();
    let mut limite = (texto.len() - patron.len()).max(BITS_POR_PALABRA);
    loop {
        if limite * 4 >= patron.len() {
            limite = patron.len() + texto.len();
        }
        let pasada = pasada_en_banda(patron, texto, &tabla, limite);
        pasadas.push(pasada);
        if let Some(distancia) = pasada.distancia {
            return DistanciaMyers {
                distancia,
                patron_es_v,
                palabras_por_columna: tabla.palabras,
                pasadas,
            };
        }
        limite *= 2;
    }
}

pub fn pasada_en_banda(
    patron: &[u8],
    texto: &[u8],
    tabla: &TablaDeIguales,
    limite: usize,
) -> Pasada {
    if patron.is_empty() {
        return Pasada {
            limite,
            distancia: Some(texto.len()),
            palabras_calculadas: 0,
        };
    }
    let resta = texto.len() - patron.len();
    let filas_arriba = (limite + resta) / 2;
    let filas_abajo = (limite - resta) / 2;
    let ultima_palabra = tabla.palabras - 1;
    let bit_de_la_ultima_fila = 1 << ((patron.len() - 1) % BITS_POR_PALABRA);
    let mut palabras = vec![Palabra::inicial(); tabla.palabras];
    let mut primera = 0;
    let mut ultima = 0;
    let mut fondo = filas_de_la_palabra(patron.len(), 0);
    let mut palabras_calculadas = 0;
    for (columna, &letra) in texto.iter().enumerate() {
        let j = columna + 1;
        let fila_mas_baja = (j + filas_abajo).min(patron.len());
        while ultima < (fila_mas_baja - 1) / BITS_POR_PALABRA {
            ultima += 1;
            fondo += filas_de_la_palabra(patron.len(), ultima);
        }
        let fila_mas_alta = j.saturating_sub(filas_arriba).max(1);
        primera = primera.max((fila_mas_alta - 1) / BITS_POR_PALABRA);
        let iguales = &tabla.de(letra)[primera..=ultima];
        let mut horizontal = 1;
        for (numero, (palabra, &mascara)) in palabras[primera..=ultima]
            .iter_mut()
            .zip(iguales)
            .enumerate()
        {
            let bit_alto = if primera + numero == ultima_palabra {
                bit_de_la_ultima_fila
            } else {
                BIT_ALTO
            };
            horizontal = avanzar_palabra(palabra, mascara, horizontal, bit_alto);
        }
        fondo = fondo.saturating_add_signed(horizontal);
        palabras_calculadas += ultima - primera + 1;
    }
    let distancia = (fondo <= limite).then_some(fondo);
    Pasada {
        limite,
        distancia,
        palabras_calculadas,
    }
}

fn filas_de_la_palabra(filas: usize, numero: usize) -> usize {
    (filas - numero * BITS_POR_PALABRA).min(BITS_POR_PALABRA)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Palabra {
    pub vertical_mas: u64,
    pub vertical_menos: u64,
}

impl Palabra {
    pub fn inicial() -> Palabra {
        Palabra {
            vertical_mas: u64::MAX,
            vertical_menos: 0,
        }
    }
}

pub struct TablaDeIguales {
    pub palabras: usize,
    indice: Vec<usize>,
    mascaras: Vec<u64>,
}

impl TablaDeIguales {
    pub fn para(patron: &[u8]) -> TablaDeIguales {
        let palabras = patron.len().div_ceil(BITS_POR_PALABRA);
        let mut indice = vec![0; 256];
        let mut distintas = 0;
        for &letra in patron {
            if indice[letra as usize] == 0 {
                distintas += 1;
                indice[letra as usize] = distintas;
            }
        }
        let mut mascaras = vec![0; (distintas + 1) * palabras];
        for (posicion, &letra) in patron.iter().enumerate() {
            let fila = indice[letra as usize];
            let numero = posicion / BITS_POR_PALABRA;
            mascaras[fila * palabras + numero] |= 1 << (posicion % BITS_POR_PALABRA);
        }
        TablaDeIguales {
            palabras,
            indice,
            mascaras,
        }
    }

    pub fn de(&self, letra: u8) -> &[u64] {
        let fila = self.indice[letra as usize];
        &self.mascaras[fila * self.palabras..(fila + 1) * self.palabras]
    }
}

pub fn avanzar_palabra(
    palabra: &mut Palabra,
    iguales: u64,
    entrada: isize,
    bit_alto: u64,
) -> isize {
    let Palabra {
        vertical_mas,
        vertical_menos,
    } = *palabra;
    let entra_mas = (entrada > 0) as u64;
    let entra_menos = (entrada < 0) as u64;
    let x_vertical = iguales | vertical_menos;
    let iguales = iguales | entra_menos;
    let x_horizontal =
        ((iguales & vertical_mas).wrapping_add(vertical_mas) ^ vertical_mas) | iguales;
    let horizontal_mas = vertical_menos | !(x_horizontal | vertical_mas);
    let horizontal_menos = vertical_mas & x_horizontal;
    let salida =
        (horizontal_mas & bit_alto != 0) as isize - (horizontal_menos & bit_alto != 0) as isize;
    let horizontal_mas = (horizontal_mas << 1) | entra_mas;
    let horizontal_menos = (horizontal_menos << 1) | entra_menos;
    palabra.vertical_mas = horizontal_menos | !(x_vertical | horizontal_mas);
    palabra.vertical_menos = horizontal_mas & x_vertical;
    salida
}
