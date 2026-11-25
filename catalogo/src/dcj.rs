use crate::formato::miles;
use crate::pizarra::{Narracion, Pizarra};
use crate::{Algoritmo, Arista, Bloque, Campo, Entrada, Paso, Salida, Tipo};
use nucleo::dcj::{
    self, Clase, Componente, Conteo, Cromosoma, NADA, Ordenamiento, TELOMERO, TipoDcj, Vertice,
};
use nucleo::reordenamientos;

const MAXIMO_BLOQUES: usize = 1_000_000;
const MAXIMO_ESCENARIO: usize = 5_000;
const MAXIMO_NARRACION: usize = 200;
const MAXIMO_FILAS: usize = 2_000;
const MAXIMO_CARACTERES: usize = 400_000;
const MAXIMO_NODOS_DIBUJO: usize = 60;
const EXTREMOS_EN_RECORRIDO: usize = 14;
const CROMOSOMAS_DETALLADOS: usize = 6;
const COMPONENTES_DETALLADOS: usize = 12;
const OPERACIONES_DETALLADAS: usize = 6;
const LUGAR_PARA_OPERACIONES: usize = 220;
const SATURACION: f64 = 0.5;
const CERCANIA: f64 = 0.1;

const CAMPOS: &[Campo] = &[
    Campo {
        id: "a",
        etiqueta: "Genoma A (el de partida): cromosomas lineales entre corchetes [+1 -2 +3], circulares entre paréntesis (+4 +5)",
        tipo: Tipo::Texto,
        ejemplo: "[+1 +2 +3 +4 +5 +6] [+7 +8 +9] (+10 +11 +12) [+13 +14]",
    },
    Campo {
        id: "b",
        etiqueta: "Genoma B (el de llegada), con los mismos bloques: [ ] lineal, ( ) circular",
        tipo: Tipo::Texto,
        ejemplo: "[+1 -4 -3 -2 +5] (+6 +7) [+8 +9 +13 +14] [+11 +12 +10]",
    },
];

const L_ORDENAR: &str =
    "pub fn ordenar(a: &[Cromosoma], b: &[Cromosoma]) -> Result<Ordenamiento, String> {";
const L_ETIQUETAS: &str = "let etiquetas = etiquetas_comunes(a, b)?;";
const L_UNION: &str = "pareja[salida] = entrada;";
const L_CIERRE: &str = "if !cromosoma.lineal {";
const L_PAREJA_A: &str = "let pareja_a = adyacencias(a, &etiquetas)?;";
const L_PAREJA_B: &str = "let pareja_b = adyacencias(b, &etiquetas)?;";
const L_GRAFO: &str = "let componentes = grafo_de_adyacencias(&pareja_a, &pareja_b);";
const L_CAMINO_AB: &str = "_ => Clase::CaminoAB,";
const L_CAMINO_AA: &str = "(true, true) => Clase::CaminoAA,";
const L_CAMINO_BB: &str = "(false, false) => Clase::CaminoBB,";
const L_CICLO: &str = "clase: Clase::Ciclo,";
const L_FORMULA: &str = "conteo.distancia = bloques - conteo.ciclos - conteo.caminos_impares / 2;";
const L_ESCENARIO: &str = "let operaciones = escenario(pareja_a.clone(), &pareja_b);";
const L_BUSCA_UNION: &str = "if q != TELOMERO && p < q && pareja[p] != q {";
const L_UNIR: &str = "operaciones.push(unir(&mut pareja, p, q));";
const L_APLICAR: &str = "aplicar(pareja, &cortadas, &nuevas);";
const L_CORTAR: &str = "operaciones.push(cortar(&mut pareja, p));";
const L_FIN: &str = "Ok(Ordenamiento {";

pub fn algoritmos() -> Vec<Algoritmo> {
    vec![Algoritmo {
        id: "reordenamientos.dcj",
        capitulo: "reordenamientos",
        nombre: "Distancia DCJ con cromosomas lineales y circulares",
        resumen: "Un solo tipo de operación, el doble corte y empalme (DCJ), cubre inversiones, translocaciones, fusiones, fisiones y círculos; a diferencia de las 2-roturas del capítulo, acepta cromosomas lineales y da la distancia exacta N − (C + I/2) en tiempo lineal, con un escenario óptimo paso a paso.",
        campos: CAMPOS,
        funcion: "dcj::ordenar",
        ejecutar: ejecutar_dcj,
        pasos: Some(pasos_dcj),
    }]
}

fn leer(entrada: &Entrada, campo: &str, nombre: &str) -> Result<Vec<Cromosoma>, String> {
    let texto = entrada.texto(campo)?;
    dcj::leer_genoma(&texto, MAXIMO_BLOQUES).map_err(|error| format!("Genoma {nombre}. {error}"))
}

fn leer_pareja(entrada: &Entrada) -> Result<(Vec<Cromosoma>, Vec<Cromosoma>), String> {
    let a = leer(entrada, "a", "A")?;
    let b = leer(entrada, "b", "B")?;
    Ok((a, b))
}

fn texto_cromosoma(cromosoma: &Cromosoma) -> String {
    let piezas: Vec<String> = cromosoma
        .bloques
        .iter()
        .map(|bloque| format!("{bloque:+}"))
        .collect();
    if cromosoma.lineal {
        format!("[{}]", piezas.join(" "))
    } else {
        format!("({})", piezas.join(" "))
    }
}

fn texto_genoma(genoma: &[Cromosoma]) -> String {
    let piezas: Vec<String> = genoma.iter().map(texto_cromosoma).collect();
    piezas.join(" ")
}

fn recortar(texto: String, maximo: usize) -> String {
    if texto.chars().count() <= maximo {
        return texto;
    }
    let corto: String = texto.chars().take(maximo).collect();
    format!("{corto}…")
}

fn nombre(extremo: usize, etiquetas: &[u64]) -> String {
    dcj::nombre_del_extremo(extremo, etiquetas)
}

fn texto_par(par: &[Vertice; 2], etiquetas: &[u64]) -> String {
    match par {
        [primero, NADA] | [NADA, primero] => texto_vertice(*primero, etiquetas),
        [primero, segundo] => format!(
            "{} y {}",
            texto_vertice(*primero, etiquetas),
            texto_vertice(*segundo, etiquetas)
        ),
    }
}

fn texto_vertice(vertice: Vertice, etiquetas: &[u64]) -> String {
    match vertice {
        NADA => "nada".to_string(),
        (p, TELOMERO) => format!("{{{}}}", nombre(p, etiquetas)),
        (p, q) => format!("{{{}, {}}}", nombre(p, etiquetas), nombre(q, etiquetas)),
    }
}

fn vertices_del_genoma(pareja: &[usize]) -> Vec<Vertice> {
    let mut vertices = Vec::new();
    for (extremo, &vecino) in pareja.iter().enumerate().skip(1) {
        if vecino == TELOMERO {
            vertices.push((extremo, TELOMERO));
        } else if extremo < vecino {
            vertices.push((extremo, vecino));
        }
    }
    vertices
}

fn textos_de_vertices(pareja: &[usize], etiquetas: &[u64]) -> Vec<String> {
    vertices_del_genoma(pareja)
        .into_iter()
        .map(|vertice| texto_vertice(vertice, etiquetas))
        .collect()
}

fn nombre_del_tipo(tipo: TipoDcj) -> &'static str {
    match tipo {
        TipoDcj::Inversion => "inversión",
        TipoDcj::Translocacion => "translocación",
        TipoDcj::Fusion => "fusión",
        TipoDcj::Fision => "fisión",
        TipoDcj::Circularizacion => "circularización",
        TipoDcj::Linealizacion => "linealización",
    }
}

const TIPOS: [TipoDcj; 6] = [
    TipoDcj::Inversion,
    TipoDcj::Translocacion,
    TipoDcj::Fusion,
    TipoDcj::Fision,
    TipoDcj::Circularizacion,
    TipoDcj::Linealizacion,
];

fn significado_biologico(tipo: TipoDcj) -> &'static str {
    match tipo {
        TipoDcj::Inversion => {
            "Un tramo de un cromosoma se corta en sus dos bordes y se vuelve a pegar dado vuelta: cambia el orden y la hebra de los genes, no cuántos cromosomas hay."
        }
        TipoDcj::Translocacion => {
            "Dos cromosomas lineales intercambian sus extremos: genes que estaban juntos pasan a otro cromosoma, como en muchas translocaciones recíprocas de tumores."
        }
        TipoDcj::Fusion => {
            "Dos cromosomas se unen en uno (por ejemplo, el cromosoma 2 humano viene de fusionar dos cromosomas de un simio ancestral) o un círculo se integra en otro cromosoma, como un plásmido que se inserta."
        }
        TipoDcj::Fision => {
            "Un cromosoma se parte en dos, o de él se escinde un círculo, como un elemento móvil que sale del cromosoma bacteriano."
        }
        TipoDcj::Circularizacion => {
            "Los dos telómeros de un cromosoma lineal se unen y queda circular, como el genoma del fago lambda al entrar en la célula."
        }
        TipoDcj::Linealizacion => {
            "Un cromosoma circular se abre en un punto y queda lineal, con dos telómeros nuevos."
        }
    }
}

fn nombre_de_clase(clase: Clase) -> &'static str {
    match clase {
        Clase::Ciclo => "ciclo",
        Clase::CaminoAB => "camino A–B (impar)",
        Clase::CaminoAA => "camino A–A (par)",
        Clase::CaminoBB => "camino B–B (par)",
    }
}

fn aporte(clase: Clase) -> &'static str {
    match clase {
        Clase::Ciclo => "1",
        Clase::CaminoAB => "½",
        Clase::CaminoAA | Clase::CaminoBB => "0",
    }
}

fn recorrido(componente: &Componente, etiquetas: &[u64]) -> String {
    let mut piezas: Vec<String> = componente
        .extremos
        .iter()
        .take(EXTREMOS_EN_RECORRIDO)
        .map(|&extremo| nombre(extremo, etiquetas))
        .collect();
    if componente.extremos.len() > EXTREMOS_EN_RECORRIDO {
        piezas.push("…".to_string());
    }
    piezas.join(" → ")
}

fn fila_de_componente(numero: usize, componente: &Componente, etiquetas: &[u64]) -> Vec<String> {
    vec![
        miles(numero),
        nombre_de_clase(componente.clase).to_string(),
        miles(componente.extremos.len()),
        recorrido(componente, etiquetas),
        aporte(componente.clase).to_string(),
    ]
}

fn cuantos(cantidad: usize, singular: &str, plural: &str) -> String {
    if cantidad == 1 {
        format!("1 {singular}")
    } else {
        format!("{} {plural}", miles(cantidad))
    }
}

fn todo_circular(genoma: &[Cromosoma]) -> bool {
    genoma.iter().all(|cromosoma| !cromosoma.lineal)
}

fn todo_lineal(genoma: &[Cromosoma]) -> bool {
    genoma.iter().all(|cromosoma| cromosoma.lineal)
}

fn como_permutaciones(genoma: &[Cromosoma], etiquetas: &[u64]) -> Vec<Vec<i64>> {
    genoma
        .iter()
        .map(|cromosoma| {
            cromosoma
                .bloques
                .iter()
                .map(|&bloque| {
                    let posicion =
                        etiquetas.binary_search(&bloque.unsigned_abs()).unwrap_or(0) as i64 + 1;
                    if bloque > 0 { posicion } else { -posicion }
                })
                .collect()
        })
        .collect()
}

fn texto_formula(conteo: &Conteo) -> String {
    format!(
        "d = N − (C + I/2) = {} − ({} + {}/2) = {}",
        miles(conteo.bloques),
        miles(conteo.ciclos),
        miles(conteo.caminos_impares),
        miles(conteo.distancia)
    )
}

fn grafo_de_adyacencias(ordenamiento: &Ordenamiento) -> Option<Bloque> {
    let etiquetas = &ordenamiento.etiquetas;
    let vertices_a = vertices_del_genoma(&ordenamiento.pareja_a);
    let vertices_b = vertices_del_genoma(&ordenamiento.pareja_b);
    if vertices_a.len() + vertices_b.len() > MAXIMO_NODOS_DIBUJO {
        return None;
    }
    let extremos = ordenamiento.pareja_a.len();
    let mut nodo_en_a = vec![0; extremos];
    let mut nodo_en_b = vec![0; extremos];
    let mut nodos = Vec::new();
    for &vertice in &vertices_a {
        nodo_en_a[vertice.0] = nodos.len();
        nodo_en_a[vertice.1] = nodos.len();
        nodos.push(format!("A {}", texto_vertice(vertice, etiquetas)));
    }
    for &vertice in &vertices_b {
        nodo_en_b[vertice.0] = nodos.len();
        nodo_en_b[vertice.1] = nodos.len();
        nodos.push(format!("B {}", texto_vertice(vertice, etiquetas)));
    }
    let aristas = (1..extremos)
        .map(|extremo| Arista {
            desde: nodo_en_a[extremo],
            hasta: nodo_en_b[extremo],
            etiqueta: nombre(extremo, etiquetas),
        })
        .collect();
    Some(Bloque::grafo(
        "Grafo de adyacencias: cada arista une el vértice de A y el de B que comparten un extremo",
        nodos,
        aristas,
    ))
}

fn lista_acotada(etiqueta: &str, elementos: Vec<String>, total: usize) -> Bloque {
    let mut visibles = Vec::new();
    let mut caracteres = 0;
    for elemento in elementos {
        caracteres += elemento.len();
        if caracteres > MAXIMO_CARACTERES || visibles.len() == MAXIMO_FILAS {
            break;
        }
        visibles.push(elemento);
    }
    Bloque::Lista {
        etiqueta: etiqueta.to_string(),
        elementos: visibles,
        total,
    }
}

struct Recorrido {
    genomas: Vec<Vec<Cromosoma>>,
    cuenta: [usize; 6],
    pasa_por_circulos: bool,
}

fn recorrer_escenario(ordenamiento: &Ordenamiento) -> Recorrido {
    let mut pareja = ordenamiento.pareja_a.clone();
    let mut genomas = Vec::new();
    let mut cuenta = [0; 6];
    let mut pasa_por_circulos = false;
    for operacion in &ordenamiento.operaciones {
        dcj::aplicar(&mut pareja, &operacion.cortadas, &operacion.nuevas);
        let genoma = dcj::genoma_desde_adyacencias(&pareja, &ordenamiento.etiquetas);
        pasa_por_circulos |= !todo_lineal(&genoma);
        if let Some(posicion) = TIPOS.iter().position(|tipo| *tipo == operacion.tipo) {
            cuenta[posicion] += 1;
        }
        if genomas.len() < MAXIMO_FILAS {
            genomas.push(genoma);
        }
    }
    Recorrido {
        genomas,
        cuenta,
        pasa_por_circulos,
    }
}

fn agregar_escenario(salida: &mut Salida, ordenamiento: &Ordenamiento, recorrido: &Recorrido) {
    let etiquetas = &ordenamiento.etiquetas;
    let distancia = ordenamiento.conteo.distancia;
    let filas: Vec<Vec<String>> = ordenamiento
        .operaciones
        .iter()
        .zip(&recorrido.genomas)
        .enumerate()
        .map(|(numero, (operacion, genoma))| {
            vec![
                miles(numero + 1),
                texto_par(&operacion.cortadas, etiquetas),
                texto_par(&operacion.nuevas, etiquetas),
                nombre_del_tipo(operacion.tipo).to_string(),
                miles(distancia - numero - 1),
                recortar(texto_genoma(genoma), 200),
            ]
        })
        .collect();
    salida.agregar(Bloque::tabla_con_total(
        "Escenario óptimo: una DCJ por paso",
        &[
            "Paso",
            "Se cortan",
            "Se forman",
            "Tipo",
            "Distancia que queda",
            "Genoma después",
        ],
        filas,
        distancia,
    ));
    let mut genomas = vec![texto_genoma(&dcj::genoma_desde_adyacencias(
        &ordenamiento.pareja_a,
        etiquetas,
    ))];
    genomas.extend(recorrido.genomas.iter().map(|genoma| texto_genoma(genoma)));
    salida.agregar(lista_acotada(
        "Genomas desde A hasta B",
        genomas,
        distancia + 1,
    ));
    let categorias = TIPOS
        .iter()
        .map(|tipo| nombre_del_tipo(*tipo).to_string())
        .collect();
    let valores = recorrido.cuenta.iter().map(|&n| n as f64).collect();
    salida.agregar(Bloque::barras(
        "Tipos de DCJ en este escenario",
        "Tipo",
        "Cantidad",
        categorias,
        valores,
    ));
    let significados: Vec<Vec<String>> = TIPOS
        .iter()
        .zip(recorrido.cuenta)
        .filter(|(_, veces)| *veces > 0)
        .map(|(tipo, veces)| {
            vec![
                nombre_del_tipo(*tipo).to_string(),
                miles(veces),
                significado_biologico(*tipo).to_string(),
            ]
        })
        .collect();
    salida.agregar(Bloque::tabla(
        "Qué es cada tipo de DCJ en biología",
        &["Tipo", "Veces", "Qué le pasa al genoma"],
        significados,
    ));
}

fn ejecutar_dcj(entrada: &Entrada) -> Result<Salida, String> {
    let (a, b) = leer_pareja(entrada)?;
    let etiquetas = dcj::etiquetas_comunes(&a, &b)?;
    let mut salida = Salida::default();
    let con_escenario = etiquetas.len() <= MAXIMO_ESCENARIO;
    let (conteo, ordenamiento) = if con_escenario {
        let ordenamiento = dcj::ordenar(&a, &b)?;
        (ordenamiento.conteo, Some(ordenamiento))
    } else {
        (dcj::distancia_dcj(&a, &b)?, None)
    };
    salida.agregar(Bloque::valor("Distancia DCJ", miles(conteo.distancia)));
    salida.agregar(Bloque::valor("Bloques (N)", miles(conteo.bloques)));
    salida.agregar(Bloque::valor("Ciclos (C)", miles(conteo.ciclos)));
    salida.agregar(Bloque::valor(
        "Caminos impares (I)",
        miles(conteo.caminos_impares),
    ));
    salida.agregar(Bloque::valor("Caminos pares", miles(conteo.caminos_pares)));
    salida.agregar(Bloque::valor("Fórmula", texto_formula(&conteo)));
    let comparable = todo_circular(&a) && todo_circular(&b);
    let dos_roturas = if comparable {
        Some(reordenamientos::distancia_dos_roturas(
            &como_permutaciones(&a, &etiquetas),
            &como_permutaciones(&b, &etiquetas),
        )?)
    } else {
        None
    };
    if let Some(clasica) = dos_roturas {
        salida.agregar(Bloque::valor(
            "Distancia de 2-roturas (capítulo 6)",
            miles(clasica),
        ));
    }
    let mut recorrido = None;
    if let Some(ordenamiento) = &ordenamiento {
        let filas: Vec<Vec<String>> = ordenamiento
            .componentes
            .iter()
            .take(MAXIMO_FILAS)
            .enumerate()
            .map(|(numero, componente)| fila_de_componente(numero + 1, componente, &etiquetas))
            .collect();
        salida.agregar(Bloque::tabla_con_total(
            "Componentes del grafo de adyacencias (t = cola, h = cabeza de cada bloque)",
            &["N.º", "Tipo", "Aristas", "Recorrido", "Aporta a C + I/2"],
            filas,
            ordenamiento.componentes.len(),
        ));
        if let Some(grafo) = grafo_de_adyacencias(ordenamiento) {
            salida.agregar(grafo);
        }
        if conteo.distancia > 0 {
            let hecho = recorrer_escenario(ordenamiento);
            agregar_escenario(&mut salida, ordenamiento, &hecho);
            recorrido = Some(hecho);
        }
    }
    interpretar(
        &mut salida,
        &conteo,
        dos_roturas,
        &a,
        &b,
        recorrido.as_ref(),
    );
    Ok(salida)
}

fn interpretar(
    salida: &mut Salida,
    conteo: &Conteo,
    dos_roturas: Option<usize>,
    a: &[Cromosoma],
    b: &[Cromosoma],
    recorrido: Option<&Recorrido>,
) {
    salida.explicar(format!(
        "{}: es el menor número de operaciones DCJ (cortar dos uniones y volver a pegar los cuatro extremos de otra forma) que convierten A en B.",
        texto_formula(conteo)
    ));
    if conteo.distancia == 0 {
        salida.explicar("A y B tienen exactamente las mismas adyacencias y telómeros: son el mismo genoma escrito quizá de otra forma, y no hace falta ninguna operación.");
        return;
    }
    let fraccion = conteo.distancia as f64 / conteo.bloques as f64;
    if fraccion >= SATURACION {
        salida.explicar(format!(
            "La distancia es {} de N, y la DCJ no puede pasar de casi N: los genomas están muy revueltos. Con tantas operaciones la parsimonia suele quedarse corta (saturación), así que tómala como un mínimo, no como el número real de eventos.",
            crate::formato::porcentaje(fraccion, 0)
        ));
    } else if fraccion <= CERCANIA {
        salida.explicar(format!(
            "La distancia es solo {} de N: los genomas son muy parecidos y un escenario parsimonioso como este es una explicación creíble de su historia.",
            crate::formato::porcentaje(fraccion, 0)
        ));
    }
    if conteo.caminos_impares > 0 || conteo.caminos_pares > 0 {
        salida.explicar(format!(
            "Los cromosomas lineales tienen telómeros, y en el grafo aparecen caminos además de ciclos: hay {} y {}. Un camino impar va de un telómero de A a uno de B y vale medio ciclo: los impares se aprovechan de dos en dos, por eso entra I/2. Un camino par (A–A o B–B) no suma nada, y una sola DCJ puede partirlo en dos caminos impares o cerrarlo en un ciclo.",
            cuantos(conteo.caminos_impares, "camino impar", "caminos impares"),
            cuantos(conteo.caminos_pares, "camino par", "caminos pares")
        ));
    }
    match dos_roturas {
        Some(clasica) if clasica == conteo.distancia => salida.explicar(format!(
            "Como A y B son todos circulares, «Distancia de 2-roturas» del capítulo 6 se puede usar y da lo mismo ({}): sin telómeros no hay caminos, la fórmula queda N − C y una DCJ es exactamente una 2-rotura.",
            miles(clasica)
        )),
        Some(clasica) => salida.explicar(format!(
            "«Distancia de 2-roturas» da {} y la DCJ da {}: no deberían diferir con todo circular; revisa los genomas.",
            miles(clasica),
            miles(conteo.distancia)
        )),
        None => salida.explicar(
            "«Distancia de 2-roturas» del capítulo 6 no sirve aquí porque solo admite cromosomas circulares; la DCJ trata los lineales gracias a los telómeros.",
        ),
    }
    match recorrido {
        Some(hecho) => {
            let partes: Vec<String> = TIPOS
                .iter()
                .zip(hecho.cuenta)
                .filter(|(_, veces)| *veces > 0)
                .map(|(tipo, veces)| format!("{} {}", miles(veces), nombre_del_tipo(*tipo)))
                .collect();
            let mut frase = format!(
                "Este escenario usa {}; hay muchos escenarios igual de cortos y este sigue el orden de Bergeron, Mixtacki y Stoye (primero crea cada adyacencia de B, luego cada telómero). La tabla de tipos dice qué significa cada uno en biología.",
                partes.join(", ")
            );
            if todo_lineal(a) && todo_lineal(b) && hecho.pasa_por_circulos {
                frase.push_str(" Ojo: A y B son lineales pero el escenario pasa por un círculo; si buscas solo inversiones, translocaciones, fusiones y fisiones, la distancia real puede ser algo mayor (la DCJ es una cota inferior).");
            }
            salida.explicar(frase);
        }
        None if conteo.bloques > MAXIMO_ESCENARIO => salida.explicar(format!(
            "Con más de {} bloques solo se calcula la distancia (es lineal en N); el escenario paso a paso se omite.",
            miles(MAXIMO_ESCENARIO)
        )),
        None => {}
    }
}

struct Trazo {
    narracion: Narracion,
    etiquetas: Vec<u64>,
    agregados: usize,
}

impl Trazo {
    fn agregar(&mut self, paso: Paso) {
        self.agregados += 1;
        self.narracion.agregar(paso);
    }
}

fn pasos_dcj(entrada: &Entrada) -> Result<Vec<Paso>, String> {
    let (a, b) = leer_pareja(entrada)?;
    let etiquetas = dcj::etiquetas_comunes(&a, &b)?;
    if etiquetas.len() > MAXIMO_NARRACION {
        return Err(format!(
            "Para el paso a paso usa genomas de hasta {} bloques; estos tienen {}.",
            miles(MAXIMO_NARRACION),
            miles(etiquetas.len())
        ));
    }
    let mut trazo = Trazo {
        narracion: Narracion::nueva(400),
        etiquetas,
        agregados: 0,
    };
    narrar_inicio(&mut trazo, &a, &b);
    let pareja_a = narrar_adyacencias(&mut trazo, &a, "A", L_PAREJA_A)?;
    let pareja_b = narrar_adyacencias(&mut trazo, &b, "B", L_PAREJA_B)?;
    let componentes = narrar_grafo(&mut trazo, &pareja_a, &pareja_b);
    let conteo = dcj::contar(trazo.etiquetas.len(), &componentes);
    trazo.agregar(
        Pizarra::nueva()
            .variable("bloques", miles(conteo.bloques))
            .variable("ciclos", miles(conteo.ciclos))
            .variable("caminos_impares", miles(conteo.caminos_impares))
            .variable("caminos_pares", miles(conteo.caminos_pares))
            .variable("distancia", miles(conteo.distancia))
            .paso(
                L_FORMULA,
                format!(
                    "{}. Cada ciclo vale 1, cada camino impar ½ y cada camino par 0. Esa es la distancia: ninguna ruta de DCJ más corta existe.",
                    texto_formula(&conteo)
                ),
            ),
    );
    let operaciones = dcj::escenario(pareja_a.clone(), &pareja_b);
    narrar_escenario(&mut trazo, &pareja_a, &operaciones, conteo.distancia);
    let final_ = dcj::genoma_desde_adyacencias(&pareja_b, &trazo.etiquetas);
    let ultimo = Pizarra::nueva()
        .variable("distancia", miles(conteo.distancia))
        .variable("operaciones", miles(operaciones.len()))
        .lista(
            "Genoma final (igual a B)",
            final_.iter().map(texto_cromosoma).collect(),
            &[],
        )
        .paso(
            L_FIN,
            format!(
                "Listo: A se convirtió en B con {}, justo la distancia. Mira «Qué significa» para leer los números y los tipos de operación.",
                cuantos(operaciones.len(), "DCJ", "DCJ")
            ),
        );
    Ok(trazo.narracion.terminar(ultimo))
}

fn narrar_inicio(trazo: &mut Trazo, a: &[Cromosoma], b: &[Cromosoma]) {
    trazo.agregar(
        Pizarra::nueva()
            .variable("a", recortar(texto_genoma(a), 120))
            .variable("b", recortar(texto_genoma(b), 120))
            .paso(
                L_ORDENAR,
                "Queremos el menor número de DCJ que convierten A en B. Los corchetes son cromosomas lineales y los paréntesis, circulares.",
            ),
    );
    let lista: Vec<String> = trazo
        .etiquetas
        .iter()
        .take(40)
        .enumerate()
        .map(|(posicion, etiqueta)| {
            format!(
                "bloque {etiqueta}: cola {etiqueta}t = extremo {}, cabeza {etiqueta}h = extremo {}",
                2 * posicion + 1,
                2 * posicion + 2
            )
        })
        .collect();
    trazo.agregar(
        Pizarra::nueva()
            .variable("a", recortar(texto_genoma(a), 120))
            .variable("b", recortar(texto_genoma(b), 120))
            .variable("etiquetas", miles(trazo.etiquetas.len()))
            .lista("Extremos de cada bloque", lista, &[])
            .paso(
                L_ETIQUETAS,
                format!(
                    "A y B tienen los mismos {} bloques, sin repetidos ni ceros. Cada bloque tiene dos extremos: la cola (t) y la cabeza (h); el + significa que se lee de t a h.",
                    miles(trazo.etiquetas.len())
                ),
            ),
    );
}

fn narrar_adyacencias(
    trazo: &mut Trazo,
    genoma: &[Cromosoma],
    letra: &str,
    linea_final: &str,
) -> Result<Vec<usize>, String> {
    let mut parcial: Vec<Cromosoma> = Vec::new();
    for (numero, cromosoma) in genoma.iter().enumerate() {
        parcial.push(cromosoma.clone());
        if numero >= CROMOSOMAS_DETALLADOS {
            continue;
        }
        let pareja = dcj::adyacencias(&parcial, &trazo.etiquetas)?;
        let anteriores = dcj::adyacencias(&parcial[..numero], &trazo.etiquetas)?;
        let textos = textos_de_parcial(&pareja, &parcial, &trazo.etiquetas);
        let viejos = textos_de_parcial(&anteriores, &parcial[..numero], &trazo.etiquetas);
        let resaltar: Vec<usize> = textos
            .iter()
            .enumerate()
            .filter(|(_, texto)| !viejos.contains(texto))
            .map(|(posicion, _)| posicion)
            .collect();
        let uniones = cromosoma.bloques.len() - usize::from(cromosoma.lineal);
        let nota = if cromosoma.lineal {
            format!(
                "Cromosoma lineal {}: cada par de bloques vecinos da una adyacencia (la salida de uno con la entrada del siguiente), {} en total, y sus dos puntas quedan como telómeros.",
                texto_cromosoma(cromosoma),
                miles(uniones)
            )
        } else {
            format!(
                "Cromosoma circular {}: {} adyacencias, contando la que cierra el círculo uniendo el último bloque con el primero. No tiene telómeros.",
                texto_cromosoma(cromosoma),
                miles(uniones)
            )
        };
        trazo.agregar(
            Pizarra::nueva()
                .variable("cromosoma", recortar(texto_cromosoma(cromosoma), 80))
                .variable("lineal", cromosoma.lineal)
                .lista(
                    &format!("Adyacencias y telómeros de {letra}"),
                    textos,
                    &resaltar,
                )
                .paso(if cromosoma.lineal { L_UNION } else { L_CIERRE }, nota),
        );
    }
    let pareja = dcj::adyacencias(genoma, &trazo.etiquetas)?;
    let textos = textos_de_vertices(&pareja, &trazo.etiquetas);
    let telomeros = pareja[1..].iter().filter(|&&v| v == TELOMERO).count();
    trazo.agregar(
        Pizarra::nueva()
            .variable("cromosomas", miles(genoma.len()))
            .variable("telomeros", miles(telomeros))
            .lista(&format!("Adyacencias y telómeros de {letra}"), textos, &[])
            .paso(
                linea_final,
                format!(
                    "El genoma {letra} queda descrito solo por sus adyacencias {{x, y}} y sus telómeros {{x}}: {} telómeros, dos por cromosoma lineal. Esta lista es todo lo que la DCJ necesita saber.",
                    miles(telomeros)
                ),
            ),
    );
    Ok(pareja)
}

fn textos_de_parcial(pareja: &[usize], parcial: &[Cromosoma], etiquetas: &[u64]) -> Vec<String> {
    let mut usados = vec![false; pareja.len()];
    for cromosoma in parcial {
        for &bloque in &cromosoma.bloques {
            let posicion = etiquetas.binary_search(&bloque.unsigned_abs()).unwrap_or(0);
            usados[2 * posicion + 1] = true;
            usados[2 * posicion + 2] = true;
        }
    }
    vertices_del_genoma(pareja)
        .into_iter()
        .filter(|vertice| usados[vertice.0])
        .map(|vertice| texto_vertice(vertice, etiquetas))
        .collect()
}

fn narrar_grafo(trazo: &mut Trazo, pareja_a: &[usize], pareja_b: &[usize]) -> Vec<Componente> {
    let componentes = dcj::grafo_de_adyacencias(pareja_a, pareja_b);
    trazo.agregar(
        Pizarra::nueva()
            .variable("componentes", miles(componentes.len()))
            .paso(
                L_GRAFO,
                "Grafo de adyacencias: un vértice por cada adyacencia o telómero de A y de B, y una arista por cada extremo, que une el vértice de A con el de B donde aparece. Cada vértice toca a lo sumo dos aristas, así que el grafo se parte en ciclos y caminos.",
            ),
    );
    let total = componentes.len();
    let salto = (total / 40).max(1);
    let mut filas: Vec<Vec<String>> = Vec::new();
    let mut ciclos = 0;
    let mut caminos_impares = 0;
    for (numero, componente) in componentes.iter().enumerate() {
        match componente.clase {
            Clase::Ciclo => ciclos += 1,
            Clase::CaminoAB => caminos_impares += 1,
            Clase::CaminoAA | Clase::CaminoBB => {}
        }
        filas.push(fila_de_componente(numero + 1, componente, &trazo.etiquetas));
        let detallado = numero < COMPONENTES_DETALLADOS || numero % salto == 0;
        if !detallado || !trazo.narracion.hay_lugar() {
            continue;
        }
        let visibles: Vec<Vec<String>> = filas.iter().rev().take(10).rev().cloned().collect();
        let fila_actual = visibles.len() - 1;
        trazo.agregar(
            Pizarra::nueva()
                .variable("clase", nombre_de_clase(componente.clase))
                .variable("aristas", miles(componente.extremos.len()))
                .variable("ciclos", miles(ciclos))
                .variable("caminos_impares", miles(caminos_impares))
                .tabla(
                    "Componentes encontrados",
                    &["N.º", "Tipo", "Aristas", "Recorrido", "Aporta"],
                    visibles,
                    &[(fila_actual, 1), (fila_actual, 3)],
                )
                .paso(
                    linea_de_clase(componente.clase),
                    nota_de_componente(componente, &trazo.etiquetas),
                ),
        );
    }
    componentes
}

fn linea_de_clase(clase: Clase) -> &'static str {
    match clase {
        Clase::Ciclo => L_CICLO,
        Clase::CaminoAB => L_CAMINO_AB,
        Clase::CaminoAA => L_CAMINO_AA,
        Clase::CaminoBB => L_CAMINO_BB,
    }
}

fn nota_de_componente(componente: &Componente, etiquetas: &[u64]) -> String {
    let inicio = nombre(componente.extremos[0], etiquetas);
    let aristas = componente.extremos.len();
    match componente.clase {
        Clase::Ciclo if aristas == 2 => format!(
            "Desde {inicio} volví al mismo lugar en 2 aristas: A y B comparten esta adyacencia, es un ciclo trivial y ya está en su sitio."
        ),
        Clase::Ciclo => format!(
            "Desde {inicio} fui saltando de A a B y de B a A hasta volver: un ciclo de {aristas} aristas. Romperlo en ciclos triviales cuesta {} DCJ.",
            miles(aristas / 2 - 1)
        ),
        Clase::CaminoAB => format!(
            "Empecé en el telómero {{{inicio}}} de A y, tras {}, terminé en un telómero de B: camino impar (A–B). Vale ½.",
            cuantos(aristas, "arista", "aristas")
        ),
        Clase::CaminoAA => format!(
            "Empecé en el telómero {{{inicio}}} de A y terminé en otro telómero de A tras {aristas} aristas: camino par (A–A). No suma nada."
        ),
        Clase::CaminoBB => format!(
            "Este camino va de un telómero de B a otro telómero de B ({aristas} aristas), empezando en {{{inicio}}}: camino par (B–B). No suma nada."
        ),
    }
}

fn narrar_escenario(
    trazo: &mut Trazo,
    pareja_a: &[usize],
    operaciones: &[dcj::Operacion],
    distancia: usize,
) {
    trazo.agregar(
        Pizarra::nueva()
            .variable("distancia", miles(distancia))
            .lista(
                "Genoma actual (A)",
                dcj::genoma_desde_adyacencias(pareja_a, &trazo.etiquetas)
                    .iter()
                    .map(texto_cromosoma)
                    .collect(),
                &[],
            )
            .paso(
                L_ESCENARIO,
                "Ahora el escenario: recorro las adyacencias de B y creo cada una que A no tenga; después corto lo que en B es telómero. Cada DCJ sube C + I/2 en 1, así que baja la distancia en exactamente 1.",
            ),
    );
    let mut pareja = pareja_a.to_vec();
    let mut anterior: Vec<String> = dcj::genoma_desde_adyacencias(&pareja, &trazo.etiquetas)
        .iter()
        .map(texto_cromosoma)
        .collect();
    let lugar = LUGAR_PARA_OPERACIONES.min(400usize.saturating_sub(trazo.agregados + 20));
    let salto = operaciones.len().div_ceil(lugar.max(1)).max(1);
    for (numero, operacion) in operaciones.iter().enumerate() {
        dcj::aplicar(&mut pareja, &operacion.cortadas, &operacion.nuevas);
        let restante = distancia - numero - 1;
        let es_ultima = numero + 1 == operaciones.len();
        if numero % salto != 0 && !es_ultima {
            continue;
        }
        let actual: Vec<String> = dcj::genoma_desde_adyacencias(&pareja, &trazo.etiquetas)
            .iter()
            .map(texto_cromosoma)
            .collect();
        let nuevos: Vec<usize> = actual
            .iter()
            .enumerate()
            .filter(|(_, texto)| !anterior.contains(texto))
            .map(|(posicion, _)| posicion)
            .collect();
        let es_union = operacion.cortadas[1] != NADA;
        let cortadas = texto_par(&operacion.cortadas, &trazo.etiquetas);
        let nuevas = texto_par(&operacion.nuevas, &trazo.etiquetas);
        let tipo = nombre_del_tipo(operacion.tipo);
        if numero < OPERACIONES_DETALLADAS {
            let (linea, nota) = if es_union {
                (
                    L_BUSCA_UNION,
                    format!(
                        "B tiene la adyacencia {} y A todavía no: hay que crearla.",
                        texto_vertice(operacion.nuevas[0], &trazo.etiquetas)
                    ),
                )
            } else {
                (
                    L_CORTAR,
                    format!(
                        "En B, {} es un telómero, pero en A sigue pegado en {}: hay que cortar esa adyacencia.",
                        texto_vertice((operacion.nuevas[0].0, TELOMERO), &trazo.etiquetas),
                        texto_vertice(operacion.cortadas[0], &trazo.etiquetas)
                    ),
                )
            };
            trazo.agregar(
                Pizarra::nueva()
                    .variable("paso", miles(numero + 1))
                    .variable("cortadas", &cortadas)
                    .variable("nuevas", &nuevas)
                    .variable("tipo", "?")
                    .variable("distancia", miles(restante + 1))
                    .lista("Genoma actual", anterior.clone(), &[])
                    .paso(linea, nota),
            );
        }
        let linea = if numero < OPERACIONES_DETALLADAS {
            L_APLICAR
        } else if es_union {
            L_UNIR
        } else {
            L_CORTAR
        };
        trazo.agregar(
            Pizarra::nueva()
                .variable("paso", miles(numero + 1))
                .variable("cortadas", &cortadas)
                .variable("nuevas", &nuevas)
                .variable("tipo", tipo)
                .variable("distancia", miles(restante))
                .lista("Genoma actual", actual.clone(), &nuevos)
                .paso(
                    linea,
                    format!(
                        "DCJ {}: corto {cortadas} y formo {nuevas}. Fue una {tipo}. Queda{} {}.",
                        miles(numero + 1),
                        if restante == 1 { "" } else { "n" },
                        cuantos(restante, "operación", "operaciones")
                    ),
                ),
        );
        anterior = actual;
    }
}

#[cfg(test)]
mod pruebas {
    use super::*;

    fn ejemplo() -> Entrada {
        Entrada::nueva(&[("a", CAMPOS[0].ejemplo), ("b", CAMPOS[1].ejemplo)])
    }

    #[test]
    fn el_ejemplo_tiene_todos_los_tipos() {
        let (a, b) = leer_pareja(&ejemplo()).unwrap();
        let ordenamiento = dcj::ordenar(&a, &b).unwrap();
        let hecho = recorrer_escenario(&ordenamiento);
        assert!(
            hecho.cuenta.iter().all(|&veces| veces > 0),
            "{:?}",
            hecho.cuenta
        );
    }

    #[test]
    fn el_paso_a_paso_da_lo_mismo_que_el_nucleo() {
        let entrada = ejemplo();
        let (a, b) = leer_pareja(&entrada).unwrap();
        let ordenamiento = dcj::ordenar(&a, &b).unwrap();
        let pasos = pasos_dcj(&entrada).unwrap();
        assert!(pasos.len() >= 15 && pasos.len() <= 400);
        let ultimo = pasos.last().unwrap();
        let variables = ultimo.estado["variables"].as_array().unwrap();
        assert_eq!(
            variables[0][1].as_str().unwrap(),
            miles(ordenamiento.conteo.distancia)
        );
        assert_eq!(
            variables[1][1].as_str().unwrap(),
            miles(ordenamiento.operaciones.len())
        );
        let codigo: String = crate::codigo("reordenamientos.dcj")
            .unwrap()
            .iter()
            .map(|fragmento| fragmento.fuente.clone())
            .collect::<Vec<_>>()
            .join("\n");
        for paso in &pasos {
            assert!(
                codigo
                    .lines()
                    .any(|linea| linea.contains(paso.linea.as_str())),
                "{}",
                paso.linea
            );
        }
    }

    #[test]
    fn genomas_grandes_al_azar_no_pasan_de_400_pasos() {
        let mut azar = nucleo::azar::Azar::nuevo(3);
        for bloques in [1, 2, 7, 50, 200] {
            let mut orden: Vec<i64> = (1..=bloques as i64).collect();
            for i in (1..bloques).rev() {
                let j = azar.entero(i + 1);
                orden.swap(i, j);
            }
            let piezas: Vec<String> = orden
                .iter()
                .map(|b| {
                    if azar.entero(2) == 0 {
                        format!("+{b}")
                    } else {
                        format!("-{b}")
                    }
                })
                .collect();
            let a: Vec<String> = (1..=bloques).map(|b| format!("+{b}")).collect();
            let entrada = Entrada::nueva(&[
                ("a", &format!("[{}]", a.join(" "))),
                ("b", &format!("({})", piezas.join(" "))),
            ]);
            let pasos = pasos_dcj(&entrada).unwrap();
            assert!(pasos.len() <= 400);
            let salida = ejecutar_dcj(&entrada).unwrap();
            assert!(!salida.interpretacion.is_empty());
        }
    }

    #[test]
    fn todo_circular_coincide_con_dos_roturas() {
        let entrada =
            Entrada::nueva(&[("a", "(+1 +2 +3 +4)(+5 +6)"), ("b", "(+1 -3 -2)(+4 +6 +5)")]);
        let salida = ejecutar_dcj(&entrada).unwrap();
        let mut dcj = String::new();
        let mut clasica = String::new();
        for bloque in &salida.bloques {
            if let Bloque::Valor { etiqueta, valor } = bloque {
                if etiqueta == "Distancia DCJ" {
                    dcj = valor.clone();
                }
                if etiqueta == "Distancia de 2-roturas (capítulo 6)" {
                    clasica = valor.clone();
                }
            }
        }
        assert!(!dcj.is_empty());
        assert_eq!(dcj, clasica);
    }

    #[test]
    fn errores_en_castellano() {
        for (a, b) in [
            ("[+1 +2]", "[+1 +3]"),
            ("[+1 0]", "[+1 +2]"),
            ("[+1 +1]", "[+1 +2]"),
            ("+1 +2", "[+1 +2]"),
            ("[+1 +2", "[+1 +2]"),
        ] {
            let entrada = Entrada::nueva(&[("a", a), ("b", b)]);
            assert!(ejecutar_dcj(&entrada).is_err());
            assert!(pasos_dcj(&entrada).is_err());
        }
    }

    #[test]
    fn el_bloque_mas_grande_no_hace_saltar_nada() {
        let entrada = Entrada::nueva(&[
            ("a", "[-9223372036854775808 +1] (+5)"),
            ("b", "[+1 -9223372036854775808] (-5)"),
        ]);
        assert!(ejecutar_dcj(&entrada).is_err());
        assert!(pasos_dcj(&entrada).is_err());
        let entrada = Entrada::nueva(&[
            ("a", "[-9223372036854775807 +1] (+5)"),
            ("b", "[+1 -9223372036854775807] (-5)"),
        ]);
        assert!(ejecutar_dcj(&entrada).is_ok());
        assert!(pasos_dcj(&entrada).is_ok());
    }
}
