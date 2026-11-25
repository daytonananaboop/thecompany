use catalogo::{Tipo, capitulos, codigo, ejecutar, ejemplo, fichas, pasos, todos};
use serde_json::Value;
use std::collections::HashSet;
use std::sync::mpsc;
use std::thread;
use std::time::Duration;

fn con_limite<T: Send + 'static>(
    segundos: u64,
    trabajo: impl FnOnce() -> T + Send + 'static,
) -> Option<T> {
    let (envio, recepcion) = mpsc::channel();
    thread::spawn(move || {
        let _ = envio.send(trabajo());
    });
    recepcion.recv_timeout(Duration::from_secs(segundos)).ok()
}

#[test]
fn identificadores_unicos_y_capitulos_validos() {
    let capitulos: HashSet<&str> = capitulos().iter().map(|c| c.id).collect();
    let mut vistos = HashSet::new();
    for algoritmo in todos() {
        assert!(vistos.insert(algoritmo.id), "id repetido: {}", algoritmo.id);
        assert!(
            capitulos.contains(algoritmo.capitulo),
            "{}: capítulo desconocido",
            algoritmo.id
        );
        assert!(
            algoritmo
                .id
                .starts_with(&format!("{}.", algoritmo.capitulo))
                || algoritmo.capitulo == "puerto",
            "{}: el id debe empezar con el capítulo",
            algoritmo.id
        );
        let mut campos = HashSet::new();
        for campo in algoritmo.campos {
            assert!(campos.insert(campo.id), "{}: campo repetido", algoritmo.id);
            assert!(!campo.etiqueta.is_empty());
        }
        assert!(!algoritmo.nombre.is_empty() && !algoritmo.resumen.is_empty());
    }
    assert_eq!(fichas().len(), todos().len());
}

#[test]
fn cada_algoritmo_funciona_con_su_ejemplo() {
    for algoritmo in todos() {
        let id = algoritmo.id;
        let entrada = ejemplo(id).unwrap();
        let copia = entrada.clone();
        let salida = con_limite(60, move || ejecutar(id, copia))
            .unwrap_or_else(|| panic!("{id}: tardó más de 60 s con su ejemplo"))
            .unwrap_or_else(|error| panic!("{id}: falló con su ejemplo: {error}"));
        assert!(!salida.bloques.is_empty(), "{id}: sin resultados");
        assert!(
            !salida.interpretacion.is_empty(),
            "{id}: sin interpretación"
        );
        let fragmentos = codigo(id).unwrap_or_else(|error| panic!("{id}: {error}"));
        assert!(!fragmentos.is_empty(), "{id}: sin código");
        assert!(
            fragmentos[0].fuente.contains("fn "),
            "{id}: el primer fragmento no es una función"
        );
        if algoritmo.pasos.is_some() {
            let copia = entrada.clone();
            let lista = con_limite(60, move || pasos(id, copia))
                .unwrap_or_else(|| panic!("{id}: el paso a paso tardó más de 60 s"))
                .unwrap_or_else(|error| panic!("{id}: paso a paso falló: {error}"));
            assert!(
                !lista.is_empty() && lista.len() <= 400,
                "{id}: {} pasos",
                lista.len()
            );
            let texto: String = fragmentos
                .iter()
                .map(|f| f.fuente.as_str())
                .collect::<Vec<_>>()
                .join("\n");
            for paso in &lista {
                assert!(!paso.linea.trim().is_empty(), "{id}: paso sin línea");
                assert!(
                    texto
                        .lines()
                        .any(|linea| linea.contains(paso.linea.as_str())),
                    "{id}: la línea «{}» no está en el código mostrado",
                    paso.linea
                );
                assert!(!paso.nota.is_empty(), "{id}: paso sin nota");
            }
        }
    }
}

#[test]
fn cada_capitulo_con_algoritmos_tiene_paso_a_paso() {
    for capitulo in capitulos() {
        if capitulo.id == "puerto" {
            continue;
        }
        let del_capitulo: Vec<_> = todos()
            .iter()
            .filter(|a| a.capitulo == capitulo.id)
            .collect();
        if del_capitulo.len() > 1 {
            assert!(
                del_capitulo.iter().any(|a| a.pasos.is_some()),
                "el capítulo {} no tiene paso a paso",
                capitulo.id
            );
        }
    }
}

const MALOS: &[&str] = &[
    "", "   ", "0", "-5", "99999999", "xyz", "ACGTN", ">s\nAC", "A", "1 2\n3", "@ecoli", "@nada",
];

#[test]
fn ninguna_entrada_mala_rompe_ni_cuelga() {
    for algoritmo in todos() {
        let id = algoritmo.id;
        for campo in algoritmo.campos {
            if campo.tipo == Tipo::Archivo {
                continue;
            }
            for malo in MALOS {
                let mut entrada = ejemplo(id).unwrap();
                if let Value::Object(valores) = &mut entrada {
                    valores.insert(campo.id.to_string(), Value::String(malo.to_string()));
                }
                let copia = entrada.clone();
                let resultado = con_limite(30, move || ejecutar(id, copia));
                let resultado = resultado.unwrap_or_else(|| {
                    panic!("{id}: se colgó con «{malo}» en el campo «{}»", campo.id)
                });
                if let Err(error) = &resultado {
                    assert!(
                        !error.starts_with("Error interno"),
                        "{id}: entró en pánico con «{malo}» en «{}»: {error}",
                        campo.id
                    );
                }
                if algoritmo.pasos.is_some() {
                    let copia = entrada.clone();
                    let resultado = con_limite(30, move || pasos(id, copia)).unwrap_or_else(|| {
                        panic!(
                            "{id}: el paso a paso se colgó con «{malo}» en «{}»",
                            campo.id
                        )
                    });
                    if let Err(error) = &resultado {
                        assert!(
                            !error.starts_with("Error interno"),
                            "{id}: el paso a paso entró en pánico con «{malo}» en «{}»: {error}",
                            campo.id
                        );
                    }
                    if let Ok(lista) = resultado {
                        assert!(lista.len() <= 400, "{id}: demasiados pasos");
                    }
                }
            }
        }
    }
}

#[test]
fn todos_los_algoritmos_tienen_paso_a_paso() {
    for algoritmo in todos() {
        assert!(
            algoritmo.pasos.is_some(),
            "{} no tiene paso a paso",
            algoritmo.id
        );
    }
}

#[test]
fn cada_referencia_de_vanguardia_es_un_algoritmo() {
    let con_referencia = fichas().iter().filter(|f| f.vanguardia.is_some()).count();
    assert_eq!(con_referencia, 14);
}
