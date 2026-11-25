use crate::formato::reducir_puntos;
use serde::Serialize;

const MAXIMO_ELEMENTOS: usize = 2000;
const MAXIMO_PUNTOS: usize = 2000;

#[derive(Serialize, Debug, Clone, Default)]
pub struct Salida {
    pub bloques: Vec<Bloque>,
    pub interpretacion: Vec<String>,
}

impl Salida {
    pub fn agregar(&mut self, bloque: Bloque) {
        self.bloques.push(bloque);
    }

    pub fn explicar(&mut self, texto: impl Into<String>) {
        self.interpretacion.push(texto.into());
    }
}

#[derive(Serialize, Debug, Clone, PartialEq)]
pub struct Marca {
    pub x: f64,
    pub texto: String,
}

#[derive(Serialize, Debug, Clone, PartialEq)]
pub struct Punto {
    pub x: f64,
    pub y: f64,
    pub grupo: usize,
    pub nombre: String,
}

#[derive(Serialize, Debug, Clone, PartialEq)]
pub struct Arista {
    pub desde: usize,
    pub hasta: usize,
    pub etiqueta: String,
}

#[derive(Serialize, Debug, Clone, PartialEq)]
pub struct FilaAlineada {
    pub nombre: String,
    pub secuencia: String,
}

#[derive(Serialize, Debug, Clone, PartialEq)]
#[serde(tag = "tipo", rename_all = "snake_case")]
pub enum Bloque {
    Valor {
        etiqueta: String,
        valor: String,
    },
    Texto {
        etiqueta: String,
        texto: String,
    },
    Lista {
        etiqueta: String,
        elementos: Vec<String>,
        total: usize,
    },
    Tabla {
        etiqueta: String,
        columnas: Vec<String>,
        filas: Vec<Vec<String>>,
        total: usize,
    },
    Linea {
        etiqueta: String,
        eje_x: String,
        eje_y: String,
        puntos: Vec<[f64; 2]>,
        marcas: Vec<Marca>,
    },
    Barras {
        etiqueta: String,
        eje_x: String,
        eje_y: String,
        categorias: Vec<String>,
        valores: Vec<f64>,
    },
    Dispersion {
        etiqueta: String,
        eje_x: String,
        eje_y: String,
        puntos: Vec<Punto>,
    },
    Calor {
        etiqueta: String,
        filas: Vec<String>,
        columnas: Vec<String>,
        valores: Vec<Vec<f64>>,
    },
    Arbol {
        etiqueta: String,
        newick: String,
    },
    Grafo {
        etiqueta: String,
        nodos: Vec<String>,
        aristas: Vec<Arista>,
    },
    Alineamiento {
        etiqueta: String,
        filas: Vec<FilaAlineada>,
    },
    Archivo {
        etiqueta: String,
        ruta: String,
    },
}

impl Bloque {
    pub fn valor(etiqueta: &str, valor: impl ToString) -> Bloque {
        Bloque::Valor {
            etiqueta: etiqueta.to_string(),
            valor: valor.to_string(),
        }
    }

    pub fn texto(etiqueta: &str, texto: impl Into<String>) -> Bloque {
        Bloque::Texto {
            etiqueta: etiqueta.to_string(),
            texto: texto.into(),
        }
    }

    pub fn lista<T: ToString>(etiqueta: &str, elementos: impl IntoIterator<Item = T>) -> Bloque {
        let mut total = 0;
        let mut visibles = Vec::new();
        for elemento in elementos {
            total += 1;
            if visibles.len() < MAXIMO_ELEMENTOS {
                visibles.push(elemento.to_string());
            }
        }
        Bloque::Lista {
            etiqueta: etiqueta.to_string(),
            elementos: visibles,
            total,
        }
    }

    pub fn tabla(etiqueta: &str, columnas: &[&str], filas: Vec<Vec<String>>) -> Bloque {
        let total = filas.len();
        Bloque::tabla_con_total(etiqueta, columnas, filas, total)
    }

    pub fn tabla_con_total(
        etiqueta: &str,
        columnas: &[&str],
        filas: Vec<Vec<String>>,
        total: usize,
    ) -> Bloque {
        Bloque::Tabla {
            etiqueta: etiqueta.to_string(),
            columnas: columnas.iter().map(|columna| columna.to_string()).collect(),
            filas: filas.into_iter().take(MAXIMO_ELEMENTOS).collect(),
            total,
        }
    }

    pub fn linea(
        etiqueta: &str,
        eje_x: &str,
        eje_y: &str,
        puntos: &[[f64; 2]],
        marcas: Vec<Marca>,
    ) -> Bloque {
        Bloque::Linea {
            etiqueta: etiqueta.to_string(),
            eje_x: eje_x.to_string(),
            eje_y: eje_y.to_string(),
            puntos: reducir_puntos(puntos, MAXIMO_PUNTOS),
            marcas,
        }
    }

    pub fn barras(
        etiqueta: &str,
        eje_x: &str,
        eje_y: &str,
        categorias: Vec<String>,
        valores: Vec<f64>,
    ) -> Bloque {
        Bloque::Barras {
            etiqueta: etiqueta.to_string(),
            eje_x: eje_x.to_string(),
            eje_y: eje_y.to_string(),
            categorias,
            valores,
        }
    }

    pub fn dispersion(etiqueta: &str, eje_x: &str, eje_y: &str, puntos: Vec<Punto>) -> Bloque {
        Bloque::Dispersion {
            etiqueta: etiqueta.to_string(),
            eje_x: eje_x.to_string(),
            eje_y: eje_y.to_string(),
            puntos: muestra_pareja(puntos, MAXIMO_PUNTOS),
        }
    }

    pub fn calor(
        etiqueta: &str,
        filas: Vec<String>,
        columnas: Vec<String>,
        valores: Vec<Vec<f64>>,
    ) -> Bloque {
        Bloque::Calor {
            etiqueta: etiqueta.to_string(),
            filas,
            columnas,
            valores,
        }
    }

    pub fn arbol(etiqueta: &str, newick: impl Into<String>) -> Bloque {
        Bloque::Arbol {
            etiqueta: etiqueta.to_string(),
            newick: newick.into(),
        }
    }

    pub fn grafo(etiqueta: &str, nodos: Vec<String>, aristas: Vec<Arista>) -> Bloque {
        Bloque::Grafo {
            etiqueta: etiqueta.to_string(),
            nodos,
            aristas,
        }
    }

    pub fn alineamiento(etiqueta: &str, filas: Vec<FilaAlineada>) -> Bloque {
        Bloque::Alineamiento {
            etiqueta: etiqueta.to_string(),
            filas,
        }
    }

    pub fn archivo(etiqueta: &str, ruta: &std::path::Path) -> Bloque {
        Bloque::Archivo {
            etiqueta: etiqueta.to_string(),
            ruta: ruta.display().to_string(),
        }
    }
}

fn muestra_pareja<T>(elementos: Vec<T>, maximo: usize) -> Vec<T> {
    if elementos.len() <= maximo || maximo == 0 {
        return elementos;
    }
    let total = elementos.len();
    let mut elegidos = Vec::with_capacity(maximo);
    let mut siguiente = 0;
    for (indice, elemento) in elementos.into_iter().enumerate() {
        if indice == siguiente * total / maximo {
            elegidos.push(elemento);
            siguiente += 1;
        }
    }
    elegidos
}

#[derive(Serialize, Debug, Clone, PartialEq)]
pub struct Paso {
    pub linea: String,
    pub nota: String,
    pub estado: serde_json::Value,
}

impl Paso {
    pub fn nuevo(linea: &str, nota: impl Into<String>, estado: serde_json::Value) -> Paso {
        Paso {
            linea: linea.to_string(),
            nota: nota.into(),
            estado,
        }
    }
}
