use nucleo::ecoli;
use nucleo::formatos::{Registro, leer_fasta};
use serde_json::{Map, Value};
use std::path::PathBuf;

pub struct Entrada {
    valores: Map<String, Value>,
}

impl Entrada {
    pub fn desde_json(valor: Value) -> Result<Entrada, String> {
        match valor {
            Value::Object(valores) => Ok(Entrada { valores }),
            Value::Null => Ok(Entrada {
                valores: Map::new(),
            }),
            _ => Err("La entrada debe ser un objeto con un valor por campo.".to_string()),
        }
    }

    pub fn nueva(pares: &[(&str, &str)]) -> Entrada {
        let mut valores = Map::new();
        for (campo, valor) in pares {
            valores.insert(campo.to_string(), Value::String(valor.to_string()));
        }
        Entrada { valores }
    }

    pub fn tiene(&self, campo: &str) -> bool {
        self.crudo(campo)
            .map(|texto| !texto.trim().is_empty())
            .unwrap_or(false)
    }

    fn crudo(&self, campo: &str) -> Result<String, String> {
        match self.valores.get(campo) {
            Some(Value::String(texto)) => Ok(texto.clone()),
            Some(Value::Number(numero)) => Ok(numero.to_string()),
            Some(Value::Bool(valor)) => Ok(valor.to_string()),
            Some(Value::Null) | None => Err(format!("Falta el campo «{campo}».")),
            Some(otro) => Ok(otro.to_string()),
        }
    }

    pub fn texto(&self, campo: &str) -> Result<String, String> {
        let texto = self.crudo(campo)?;
        if texto.trim().is_empty() {
            return Err(format!("El campo «{campo}» está vacío."));
        }
        Ok(texto)
    }

    pub fn secuencia(&self, campo: &str) -> Result<Vec<u8>, String> {
        let texto = self.texto(campo)?;
        if let Some(especial) = resolver_especial(texto.trim())? {
            return Ok(especial);
        }
        let registros = leer_fasta(texto.as_bytes());
        let secuencia: Vec<u8> = registros
            .into_iter()
            .next()
            .map(|registro| registro.secuencia)
            .unwrap_or_default()
            .into_iter()
            .filter(|letra| letra.is_ascii_alphabetic())
            .collect();
        if secuencia.is_empty() {
            return Err(format!(
                "El campo «{campo}» no tiene ninguna letra de secuencia."
            ));
        }
        Ok(secuencia)
    }

    pub fn adn(&self, campo: &str) -> Result<Vec<u8>, String> {
        let secuencia = self.secuencia(campo)?;
        if let Some(letra) = secuencia
            .iter()
            .find(|letra| !matches!(letra, b'A' | b'C' | b'G' | b'T'))
        {
            return Err(format!(
                "El campo «{campo}» solo admite A, C, G y T, pero contiene «{}».",
                *letra as char
            ));
        }
        Ok(secuencia)
    }

    pub fn proteina(&self, campo: &str) -> Result<Vec<u8>, String> {
        let texto = self.texto(campo)?;
        let limpio = texto.trim();
        if let Some(nombre) = limpio.strip_prefix("@proteina:") {
            let elegido = ecoli::buscar_gen(nombre.trim())
                .ok_or_else(|| format!("E. coli no tiene un gen llamado «{}».", nombre.trim()))?;
            return Ok(elegido.proteina.clone());
        }
        self.secuencia(campo)
    }

    pub fn registros(&self, campo: &str) -> Result<Vec<Registro>, String> {
        let texto = self.texto(campo)?;
        if texto.trim_start().starts_with('>') {
            let registros: Vec<Registro> = leer_fasta(texto.as_bytes())
                .into_iter()
                .filter(|registro| !registro.secuencia.is_empty())
                .collect();
            if registros.is_empty() {
                return Err(format!("El campo «{campo}» no tiene secuencias."));
            }
            return Ok(registros);
        }
        let mut registros = Vec::new();
        for pieza in texto.split(|c: char| c.is_whitespace() || c == ',' || c == ';') {
            let pieza = pieza.trim();
            if pieza.is_empty() {
                continue;
            }
            let secuencia = match resolver_especial(pieza)? {
                Some(especial) => especial,
                None => pieza
                    .bytes()
                    .map(|letra| letra.to_ascii_uppercase())
                    .collect(),
            };
            registros.push(Registro {
                nombre: format!("s{}", registros.len() + 1),
                secuencia,
            });
        }
        if registros.is_empty() {
            return Err(format!("El campo «{campo}» no tiene secuencias."));
        }
        Ok(registros)
    }

    pub fn secuencias(&self, campo: &str) -> Result<Vec<Vec<u8>>, String> {
        Ok(self
            .registros(campo)?
            .into_iter()
            .map(|registro| registro.secuencia)
            .collect())
    }

    pub fn natural(&self, campo: &str) -> Result<usize, String> {
        let texto = self.texto(campo)?;
        texto.trim().parse::<usize>().map_err(|_| {
            format!("El campo «{campo}» debe ser un número entero mayor o igual a cero.")
        })
    }

    pub fn entero(&self, campo: &str) -> Result<i64, String> {
        let texto = self.texto(campo)?;
        texto
            .trim()
            .parse::<i64>()
            .map_err(|_| format!("El campo «{campo}» debe ser un número entero."))
    }

    pub fn decimal(&self, campo: &str) -> Result<f64, String> {
        let texto = self.texto(campo)?;
        let numero = texto
            .trim()
            .replace(',', ".")
            .parse::<f64>()
            .map_err(|_| format!("El campo «{campo}» debe ser un número."))?;
        if !numero.is_finite() {
            return Err(format!("El campo «{campo}» debe ser un número finito."));
        }
        Ok(numero)
    }

    pub fn numeros(&self, campo: &str) -> Result<Vec<f64>, String> {
        let texto = self.texto(campo)?;
        let mut numeros = Vec::new();
        for pieza in texto.split(|c: char| c.is_whitespace() || c == ',' || c == ';') {
            if pieza.is_empty() {
                continue;
            }
            let numero = pieza
                .parse::<f64>()
                .map_err(|_| format!("En el campo «{campo}», «{pieza}» no es un número."))?;
            if !numero.is_finite() {
                return Err(format!(
                    "En el campo «{campo}», «{pieza}» no es un número finito."
                ));
            }
            numeros.push(numero);
        }
        if numeros.is_empty() {
            return Err(format!("El campo «{campo}» no tiene números."));
        }
        Ok(numeros)
    }

    pub fn enteros(&self, campo: &str) -> Result<Vec<i64>, String> {
        let texto = self.texto(campo)?;
        let mut numeros = Vec::new();
        for pieza in
            texto.split(|c: char| c.is_whitespace() || c == ',' || c == ';' || c == '(' || c == ')')
        {
            if pieza.is_empty() {
                continue;
            }
            let numero = pieza
                .parse::<i64>()
                .map_err(|_| format!("En el campo «{campo}», «{pieza}» no es un número entero."))?;
            numeros.push(numero);
        }
        if numeros.is_empty() {
            return Err(format!("El campo «{campo}» no tiene números."));
        }
        Ok(numeros)
    }

    pub fn matriz(&self, campo: &str) -> Result<Vec<Vec<f64>>, String> {
        let texto = self.texto(campo)?;
        let mut filas = Vec::new();
        for linea in texto.lines() {
            let mut fila = Vec::new();
            for pieza in linea.split(|c: char| c.is_whitespace() || c == ',' || c == ';') {
                if pieza.is_empty() {
                    continue;
                }
                let numero = pieza
                    .parse::<f64>()
                    .map_err(|_| format!("En el campo «{campo}», «{pieza}» no es un número."))?;
                if !numero.is_finite() {
                    return Err(format!(
                        "En el campo «{campo}», «{pieza}» no es un número finito."
                    ));
                }
                fila.push(numero);
            }
            if !fila.is_empty() {
                filas.push(fila);
            }
        }
        if filas.is_empty() {
            return Err(format!("El campo «{campo}» no tiene números."));
        }
        let columnas = filas[0].len();
        if filas.iter().any(|fila| fila.len() != columnas) {
            return Err(format!(
                "En el campo «{campo}», todas las filas deben tener la misma cantidad de números."
            ));
        }
        Ok(filas)
    }

    pub fn ruta(&self, campo: &str) -> Result<PathBuf, String> {
        let texto = self.texto(campo)?;
        let limpio = texto.trim();
        if let Some(nombre) = limpio.strip_prefix("@ejemplo:") {
            return crate::ejemplos::generar(nombre.trim());
        }
        let ruta = PathBuf::from(limpio);
        if !ruta.is_file() {
            return Err(format!("No existe el archivo {}.", ruta.display()));
        }
        Ok(ruta)
    }

    pub fn opcion(&self, campo: &str, opciones: &[&str]) -> Result<String, String> {
        let texto = self.texto(campo)?;
        let limpio = texto.trim();
        if opciones.contains(&limpio) {
            Ok(limpio.to_string())
        } else {
            Err(format!(
                "El campo «{campo}» debe ser una de estas opciones: {}.",
                opciones.join(", ")
            ))
        }
    }
}

fn resolver_especial(texto: &str) -> Result<Option<Vec<u8>>, String> {
    if texto == "@ecoli" {
        return Ok(Some(ecoli::genoma().to_vec()));
    }
    if let Some(rango) = texto.strip_prefix("@ecoli:") {
        let (inicio, fin) = rango
            .split_once('-')
            .ok_or_else(|| "Escribe la región como @ecoli:INICIO-FIN.".to_string())?;
        let inicio: usize = inicio
            .trim()
            .parse()
            .map_err(|_| "El inicio de la región debe ser un número.".to_string())?;
        let fin: usize = fin
            .trim()
            .parse()
            .map_err(|_| "El fin de la región debe ser un número.".to_string())?;
        let region = ecoli::region(inicio, fin).ok_or_else(|| {
            format!(
                "La región debe estar entre 1 y {} y el inicio no puede pasar del fin.",
                ecoli::genoma().len()
            )
        })?;
        return Ok(Some(region.to_vec()));
    }
    if let Some(nombre) = texto.strip_prefix("@gen:") {
        let elegido = ecoli::buscar_gen(nombre.trim())
            .ok_or_else(|| format!("E. coli no tiene un gen llamado «{}».", nombre.trim()))?;
        return Ok(Some(ecoli::secuencia_gen(elegido)));
    }
    if texto.starts_with('@') && !texto.contains('\n') {
        return Err(format!("No conozco la referencia «{texto}»."));
    }
    Ok(None)
}
