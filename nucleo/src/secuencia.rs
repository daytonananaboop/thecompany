pub fn complemento(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'u' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        otra => otra,
    }
}

pub fn complemento_inverso(adn: &[u8]) -> Vec<u8> {
    adn.iter().rev().map(|&base| complemento(base)).collect()
}

pub fn contenido_gc(adn: &[u8]) -> f64 {
    if adn.is_empty() {
        return 0.0;
    }
    let gc = adn
        .iter()
        .filter(|&&base| matches!(base, b'G' | b'C' | b'g' | b'c'))
        .count();
    gc as f64 / adn.len() as f64
}

pub fn transcribir(adn: &[u8]) -> Vec<u8> {
    adn.iter()
        .map(|&base| match base {
            b'T' => b'U',
            b't' => b'u',
            otra => otra,
        })
        .collect()
}

pub fn indice_base(base: u8) -> Option<usize> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

const CODIGO_GENETICO: &[u8; 64] =
    b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

pub fn aminoacido(codon: &[u8]) -> u8 {
    if codon.len() != 3 {
        return b'X';
    }
    let mut numero = 0;
    for &base in codon {
        match indice_base(base) {
            Some(indice) => numero = numero * 4 + indice,
            None => return b'X',
        }
    }
    CODIGO_GENETICO[numero]
}

pub fn traducir(secuencia: &[u8]) -> Vec<u8> {
    secuencia
        .as_chunks::<3>()
        .0
        .iter()
        .map(|codon| aminoacido(codon))
        .collect()
}

pub fn limpiar(texto: &str) -> Vec<u8> {
    texto
        .bytes()
        .filter(|byte| byte.is_ascii_alphabetic())
        .map(|byte| byte.to_ascii_uppercase())
        .collect()
}
