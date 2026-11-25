use crate::Paso;
use serde_json::{Value, json};

const VENTANA_DE_SECUENCIA: usize = 240;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Marca {
    Actual,
    Hallado,
    Comparado,
}

impl Marca {
    fn nombre(self) -> &'static str {
        match self {
            Marca::Actual => "actual",
            Marca::Hallado => "hallado",
            Marca::Comparado => "comparado",
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct Pizarra {
    variables: Vec<(String, String)>,
    secuencia: Option<Value>,
    tabla: Option<Value>,
    lista: Option<Value>,
}

impl Pizarra {
    pub fn nueva() -> Pizarra {
        Pizarra::default()
    }

    pub fn variable(mut self, nombre: &str, valor: impl ToString) -> Pizarra {
        self.variables.push((nombre.to_string(), valor.to_string()));
        self
    }

    pub fn secuencia(mut self, texto: &[u8], marcas: &[(usize, usize, Marca)]) -> Pizarra {
        let centro = marcas
            .iter()
            .find(|(_, _, tipo)| *tipo == Marca::Actual)
            .or_else(|| marcas.first())
            .map(|(desde, _, _)| *desde)
            .unwrap_or(0);
        let desde = centro
            .saturating_sub(VENTANA_DE_SECUENCIA / 3)
            .min(texto.len().saturating_sub(VENTANA_DE_SECUENCIA));
        let hasta = (desde + VENTANA_DE_SECUENCIA).min(texto.len());
        let visibles: Vec<Value> = marcas
            .iter()
            .filter(|(inicio, fin, _)| *fin > desde && *inicio < hasta)
            .map(|(inicio, fin, tipo)| {
                json!({ "desde": inicio, "hasta": fin, "tipo": tipo.nombre() })
            })
            .collect();
        self.secuencia = Some(json!({
            "texto": String::from_utf8_lossy(&texto[desde..hasta]),
            "desplazamiento": desde,
            "largo": texto.len(),
            "marcas": visibles,
        }));
        self
    }

    pub fn tabla(
        mut self,
        titulo: &str,
        columnas: &[&str],
        filas: Vec<Vec<String>>,
        resaltar: &[(usize, usize)],
    ) -> Pizarra {
        self.tabla = Some(json!({
            "titulo": titulo,
            "columnas": columnas,
            "filas": filas,
            "resaltar": resaltar.iter().map(|(fila, columna)| json!([fila, columna])).collect::<Vec<_>>(),
        }));
        self
    }

    pub fn lista(mut self, titulo: &str, elementos: Vec<String>, resaltar: &[usize]) -> Pizarra {
        self.lista = Some(json!({
            "titulo": titulo,
            "elementos": elementos,
            "resaltar": resaltar,
        }));
        self
    }

    pub fn paso(self, linea: &str, nota: impl Into<String>) -> Paso {
        let mut estado = json!({
            "tipo": "pizarra",
            "variables": self.variables.iter().map(|(nombre, valor)| json!([nombre, valor])).collect::<Vec<_>>(),
        });
        if let Some(secuencia) = self.secuencia {
            estado["secuencia"] = secuencia;
        }
        if let Some(tabla) = self.tabla {
            estado["tabla"] = tabla;
        }
        if let Some(lista) = self.lista {
            estado["lista"] = lista;
        }
        Paso::nuevo(linea, nota, estado)
    }
}

pub struct Narracion {
    pasos: Vec<Paso>,
    maximo: usize,
    omitidos: usize,
}

impl Narracion {
    pub fn nueva(maximo: usize) -> Narracion {
        Narracion {
            pasos: Vec::new(),
            maximo: maximo.clamp(2, 400),
            omitidos: 0,
        }
    }

    pub fn hay_lugar(&self) -> bool {
        self.pasos.len() + 1 < self.maximo
    }

    pub fn agregar(&mut self, paso: Paso) {
        if self.hay_lugar() {
            self.pasos.push(paso);
        } else {
            self.omitidos += 1;
        }
    }

    pub fn terminar(mut self, mut ultimo: Paso) -> Vec<Paso> {
        if self.omitidos > 0 {
            ultimo.nota = format!(
                "Se resumieron {} pasos que repetían lo mismo. {}",
                crate::formato::miles(self.omitidos),
                ultimo.nota
            );
        }
        self.pasos.push(ultimo);
        self.pasos
    }
}
