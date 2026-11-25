use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Registro {
    pub nombre: String,
    pub secuencia: Vec<u8>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Lectura {
    pub nombre: Vec<u8>,
    pub secuencia: Vec<u8>,
    pub calidad: Vec<u8>,
}

pub fn leer_fasta(texto: &[u8]) -> Vec<Registro> {
    let mut registros: Vec<Registro> = Vec::new();
    for linea in texto.split(|&byte| byte == b'\n') {
        let linea = recortar(linea);
        if linea.is_empty() {
            continue;
        }
        if linea[0] == b'>' {
            let nombre = String::from_utf8_lossy(&linea[1..]).trim().to_string();
            registros.push(Registro {
                nombre,
                secuencia: Vec::new(),
            });
            continue;
        }
        if registros.is_empty() {
            registros.push(Registro {
                nombre: "secuencia".to_string(),
                secuencia: Vec::new(),
            });
        }
        let actual = registros.last_mut().expect("hay un registro");
        actual.secuencia.extend(
            linea
                .iter()
                .filter(|byte| byte.is_ascii_alphabetic() || **byte == b'*' || **byte == b'-')
                .map(|byte| byte.to_ascii_uppercase()),
        );
    }
    registros
}

fn recortar(linea: &[u8]) -> &[u8] {
    let mut inicio = 0;
    let mut fin = linea.len();
    while inicio < fin && linea[inicio].is_ascii_whitespace() {
        inicio += 1;
    }
    while fin > inicio && linea[fin - 1].is_ascii_whitespace() {
        fin -= 1;
    }
    &linea[inicio..fin]
}

pub fn abrir(ruta: &Path) -> Result<Box<dyn BufRead + Send>, String> {
    let mut archivo = File::open(ruta)
        .map_err(|error| format!("No se pudo abrir {}: {error}", ruta.display()))?;
    let mut firma = [0u8; 2];
    let leidos = archivo
        .read(&mut firma)
        .map_err(|error| format!("No se pudo leer {}: {error}", ruta.display()))?;
    let archivo = File::open(ruta)
        .map_err(|error| format!("No se pudo abrir {}: {error}", ruta.display()))?;
    if leidos == 2 && firma == [0x1f, 0x8b] {
        Ok(Box::new(BufReader::with_capacity(
            1 << 20,
            MultiGzDecoder::new(archivo),
        )))
    } else {
        Ok(Box::new(BufReader::with_capacity(1 << 20, archivo)))
    }
}

pub fn leer_fasta_archivo(ruta: &Path) -> Result<Vec<Registro>, String> {
    let mut contenido = Vec::new();
    abrir(ruta)?
        .read_to_end(&mut contenido)
        .map_err(|error| format!("No se pudo leer {}: {error}", ruta.display()))?;
    Ok(leer_fasta(&contenido))
}

pub struct LectorFastq {
    fuente: Box<dyn BufRead + Send>,
    linea: Vec<u8>,
    numero: usize,
}

pub fn abrir_fastq(ruta: &Path) -> Result<LectorFastq, String> {
    Ok(LectorFastq::nuevo(abrir(ruta)?))
}

impl LectorFastq {
    pub fn nuevo(fuente: Box<dyn BufRead + Send>) -> LectorFastq {
        LectorFastq {
            fuente,
            linea: Vec::new(),
            numero: 0,
        }
    }

    fn siguiente_linea(&mut self) -> Result<bool, String> {
        self.linea.clear();
        let leidos = self
            .fuente
            .read_until(b'\n', &mut self.linea)
            .map_err(|error| format!("Error al leer el FASTQ: {error}"))?;
        if leidos == 0 {
            return Ok(false);
        }
        self.numero += 1;
        while matches!(self.linea.last(), Some(b'\n') | Some(b'\r')) {
            self.linea.pop();
        }
        Ok(true)
    }

    fn leer_registro(&mut self) -> Result<Option<Lectura>, String> {
        loop {
            if !self.siguiente_linea()? {
                return Ok(None);
            }
            if !self.linea.is_empty() {
                break;
            }
        }
        if self.linea[0] != b'@' {
            return Err(format!(
                "Línea {} del FASTQ: se esperaba '@' al inicio de una lectura.",
                self.numero
            ));
        }
        let nombre = self.linea[1..].to_vec();
        if !self.siguiente_linea()? {
            return Err("El FASTQ termina a mitad de una lectura.".to_string());
        }
        let secuencia: Vec<u8> = self.linea.iter().map(|b| b.to_ascii_uppercase()).collect();
        if !self.siguiente_linea()? || self.linea.first() != Some(&b'+') {
            return Err(format!(
                "Línea {} del FASTQ: se esperaba '+' después de la secuencia.",
                self.numero
            ));
        }
        if !self.siguiente_linea()? {
            return Err("El FASTQ termina a mitad de una lectura.".to_string());
        }
        let calidad = self.linea.clone();
        if calidad.len() != secuencia.len() {
            return Err(format!(
                "Línea {} del FASTQ: la calidad y la secuencia no tienen la misma longitud.",
                self.numero
            ));
        }
        Ok(Some(Lectura {
            nombre,
            secuencia,
            calidad,
        }))
    }
}

impl Iterator for LectorFastq {
    type Item = Result<Lectura, String>;

    fn next(&mut self) -> Option<Self::Item> {
        self.leer_registro().transpose()
    }
}

pub fn escribir_fastq(destino: &mut impl Write, lectura: &Lectura) -> io::Result<()> {
    destino.write_all(b"@")?;
    destino.write_all(&lectura.nombre)?;
    destino.write_all(b"\n")?;
    destino.write_all(&lectura.secuencia)?;
    destino.write_all(b"\n+\n")?;
    destino.write_all(&lectura.calidad)?;
    destino.write_all(b"\n")
}

pub fn escribir_fasta(destino: &mut impl Write, nombre: &str, secuencia: &[u8]) -> io::Result<()> {
    writeln!(destino, ">{nombre}")?;
    for trozo in secuencia.chunks(70) {
        destino.write_all(trozo)?;
        destino.write_all(b"\n")?;
    }
    Ok(())
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct RegistroSam {
    pub nombre: Vec<u8>,
    pub bandera: u16,
    pub referencia: String,
    pub posicion: usize,
    pub calidad_mapeo: u8,
    pub cigar: Vec<(u32, u8)>,
    pub secuencia: Vec<u8>,
    pub calidad: Vec<u8>,
}

impl RegistroSam {
    pub fn mapeado(&self) -> bool {
        self.bandera & 4 == 0 && self.posicion > 0
    }

    pub fn hebra_reversa(&self) -> bool {
        self.bandera & 16 != 0
    }
}

pub fn cigar_a_texto(cigar: &[(u32, u8)]) -> String {
    if cigar.is_empty() {
        return "*".to_string();
    }
    cigar
        .iter()
        .map(|(longitud, operacion)| format!("{longitud}{}", *operacion as char))
        .collect()
}

pub fn texto_a_cigar(texto: &[u8]) -> Result<Vec<(u32, u8)>, String> {
    if texto == b"*" {
        return Ok(Vec::new());
    }
    let mut cigar = Vec::new();
    let mut numero: u32 = 0;
    let mut hay_numero = false;
    for &letra in texto {
        if letra.is_ascii_digit() {
            numero = numero
                .checked_mul(10)
                .and_then(|n| n.checked_add((letra - b'0') as u32))
                .ok_or_else(|| "CIGAR con un número demasiado grande.".to_string())?;
            hay_numero = true;
        } else if b"MIDNSHP=X".contains(&letra) && hay_numero {
            cigar.push((numero, letra));
            numero = 0;
            hay_numero = false;
        } else {
            return Err(format!(
                "CIGAR inválido: {}",
                String::from_utf8_lossy(texto)
            ));
        }
    }
    if hay_numero {
        return Err(format!(
            "CIGAR inválido: {}",
            String::from_utf8_lossy(texto)
        ));
    }
    Ok(cigar)
}

pub fn escribir_cabecera_sam(
    destino: &mut impl Write,
    referencia: &str,
    longitud: usize,
) -> io::Result<()> {
    writeln!(destino, "@HD\tVN:1.6\tSO:unsorted")?;
    writeln!(destino, "@SQ\tSN:{referencia}\tLN:{longitud}")
}

pub fn escribir_sam(destino: &mut impl Write, registro: &RegistroSam) -> io::Result<()> {
    let referencia = if registro.referencia.is_empty() {
        "*"
    } else {
        registro.referencia.as_str()
    };
    destino.write_all(&registro.nombre)?;
    write!(
        destino,
        "\t{}\t{referencia}\t{}\t{}\t{}\t*\t0\t0\t",
        registro.bandera,
        registro.posicion,
        registro.calidad_mapeo,
        cigar_a_texto(&registro.cigar)
    )?;
    if registro.secuencia.is_empty() {
        destino.write_all(b"*")?;
    } else {
        destino.write_all(&registro.secuencia)?;
    }
    destino.write_all(b"\t")?;
    if registro.calidad.is_empty() {
        destino.write_all(b"*")?;
    } else {
        destino.write_all(&registro.calidad)?;
    }
    destino.write_all(b"\n")
}

pub fn leer_sam(linea: &[u8]) -> Result<Option<RegistroSam>, String> {
    let linea = recortar(linea);
    if linea.is_empty() || linea[0] == b'@' {
        return Ok(None);
    }
    let columnas: Vec<&[u8]> = linea.split(|&byte| byte == b'\t').collect();
    if columnas.len() < 11 {
        return Err("Una línea del SAM tiene menos de 11 columnas.".to_string());
    }
    let numero = |texto: &[u8], que: &str| -> Result<usize, String> {
        std::str::from_utf8(texto)
            .ok()
            .and_then(|t| t.parse::<usize>().ok())
            .ok_or_else(|| format!("En el SAM, {que} no es un número."))
    };
    let bandera = numero(columnas[1], "la bandera")?;
    let calidad_mapeo = numero(columnas[4], "la calidad de mapeo")?;
    let secuencia = if columnas[9] == b"*" {
        Vec::new()
    } else {
        columnas[9].to_ascii_uppercase()
    };
    let calidad = if columnas[10] == b"*" {
        Vec::new()
    } else {
        columnas[10].to_vec()
    };
    Ok(Some(RegistroSam {
        nombre: columnas[0].to_vec(),
        bandera: bandera.min(u16::MAX as usize) as u16,
        referencia: String::from_utf8_lossy(columnas[2]).to_string(),
        posicion: numero(columnas[3], "la posición")?,
        calidad_mapeo: calidad_mapeo.min(255) as u8,
        cigar: texto_a_cigar(columnas[5])?,
        secuencia,
        calidad,
    }))
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Variante {
    pub cromosoma: String,
    pub posicion: usize,
    pub referencia: Vec<u8>,
    pub alternativa: Vec<u8>,
    pub calidad: f64,
    pub profundidad: u32,
    pub frecuencia: f64,
}

pub fn escribir_cabecera_vcf(
    destino: &mut impl Write,
    referencia: &str,
    longitud: usize,
) -> io::Result<()> {
    writeln!(destino, "##fileformat=VCFv4.2")?;
    writeln!(destino, "##contig=<ID={referencia},length={longitud}>")?;
    writeln!(
        destino,
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Profundidad\">"
    )?;
    writeln!(
        destino,
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Frecuencia del alelo alternativo\">"
    )?;
    writeln!(destino, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
}

pub fn escribir_vcf(destino: &mut impl Write, variante: &Variante) -> io::Result<()> {
    writeln!(
        destino,
        "{}\t{}\t.\t{}\t{}\t{:.1}\tPASS\tDP={};AF={:.3}",
        variante.cromosoma,
        variante.posicion,
        String::from_utf8_lossy(&variante.referencia),
        String::from_utf8_lossy(&variante.alternativa),
        variante.calidad,
        variante.profundidad,
        variante.frecuencia
    )
}

pub fn leer_vcf(linea: &[u8]) -> Result<Option<Variante>, String> {
    let linea = recortar(linea);
    if linea.is_empty() || linea[0] == b'#' {
        return Ok(None);
    }
    let columnas: Vec<&[u8]> = linea.split(|&byte| byte == b'\t').collect();
    if columnas.len() < 8 {
        return Err("Una línea del VCF tiene menos de 8 columnas.".to_string());
    }
    let texto = |bytes: &[u8]| String::from_utf8_lossy(bytes).to_string();
    let mut variante = Variante {
        cromosoma: texto(columnas[0]),
        posicion: texto(columnas[1])
            .parse()
            .map_err(|_| "En el VCF, la posición no es un número.".to_string())?,
        referencia: columnas[3].to_ascii_uppercase(),
        alternativa: columnas[4].to_ascii_uppercase(),
        calidad: texto(columnas[5]).parse().unwrap_or(0.0),
        profundidad: 0,
        frecuencia: 0.0,
    };
    for dato in texto(columnas[7]).split(';') {
        if let Some(valor) = dato.strip_prefix("DP=") {
            variante.profundidad = valor.parse().unwrap_or(0);
        }
        if let Some(valor) = dato.strip_prefix("AF=") {
            variante.frecuencia = valor.parse().unwrap_or(0.0);
        }
    }
    Ok(Some(variante))
}
