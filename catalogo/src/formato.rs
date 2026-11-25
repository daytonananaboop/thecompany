pub fn miles(numero: usize) -> String {
    let digitos = numero.to_string();
    let mut resultado = String::new();
    for (posicion, digito) in digitos.chars().enumerate() {
        if posicion > 0 && (digitos.len() - posicion).is_multiple_of(3) {
            resultado.push('\u{a0}');
        }
        resultado.push(digito);
    }
    resultado
}

pub fn decimal(numero: f64, cifras: usize) -> String {
    format!("{numero:.cifras$}")
}

pub fn porcentaje(fraccion: f64, cifras: usize) -> String {
    format!("{}\u{a0}%", decimal(fraccion * 100.0, cifras))
}

pub fn cientifico(numero: f64) -> String {
    if numero == 0.0 {
        return "0".to_string();
    }
    if numero.abs() >= 0.001 && numero.abs() < 10_000.0 {
        return decimal(numero, 3);
    }
    format!("{numero:.1e}")
}

pub fn reducir_puntos(puntos: &[[f64; 2]], maximo: usize) -> Vec<[f64; 2]> {
    if puntos.len() <= maximo || maximo < 4 {
        return puntos.to_vec();
    }
    let grupos = maximo / 2;
    let por_grupo = puntos.len().div_ceil(grupos);
    let mut reducidos = Vec::with_capacity(maximo);
    for grupo in puntos.chunks(por_grupo) {
        let mut menor = 0;
        let mut mayor = 0;
        for (indice, punto) in grupo.iter().enumerate() {
            if punto[1] < grupo[menor][1] {
                menor = indice;
            }
            if punto[1] > grupo[mayor][1] {
                mayor = indice;
            }
        }
        let primero = menor.min(mayor);
        let segundo = menor.max(mayor);
        reducidos.push(grupo[primero]);
        if segundo != primero {
            reducidos.push(grupo[segundo]);
        }
    }
    reducidos
}
