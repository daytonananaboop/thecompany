use nucleo::azar::Azar;
use nucleo::ecoli;
use nucleo::formatos::RegistroSam;
use nucleo::secuencia::{complemento_inverso, traducir};
use nucleo::simulacion::*;

fn referencia_al_azar(longitud: usize, semilla: u64) -> Vec<u8> {
    let mut azar = Azar::nuevo(semilla);
    (0..longitud).map(|_| b"ACGT"[azar.entero(4)]).collect()
}

fn diferencias_con_referencia(registro: &RegistroSam, referencia: &[u8]) -> (usize, usize) {
    let mut en_referencia = registro.posicion - 1;
    let mut en_lectura = 0;
    let mut diferencias = 0;
    let mut comparadas = 0;
    for &(longitud, operacion) in &registro.cigar {
        let longitud = longitud as usize;
        match operacion {
            b'M' => {
                for k in 0..longitud {
                    comparadas += 1;
                    if registro.secuencia[en_lectura + k] != referencia[en_referencia + k] {
                        diferencias += 1;
                    }
                }
                en_lectura += longitud;
                en_referencia += longitud;
            }
            b'I' | b'S' => en_lectura += longitud,
            b'D' => en_referencia += longitud,
            _ => {}
        }
    }
    assert_eq!(en_lectura, registro.secuencia.len());
    (diferencias, comparadas)
}

#[test]
fn mutacion_de_resistencia_en_rpob_y_gyra() {
    let mutaciones = mutaciones_de_resistencia().unwrap();
    assert_eq!(mutaciones.len(), 2);
    assert_eq!(mutaciones[0].nombre, "rpoB S531F");
    assert_eq!(mutaciones[1].nombre, "gyrA S83L");
    for (mutacion, nombre, codon, nuevo) in [
        (&mutaciones[0], "rpoB", 531, b'F'),
        (&mutaciones[1], "gyrA", 83, b'L'),
    ] {
        let elegido = ecoli::buscar_gen(nombre).unwrap();
        let mut genoma = ecoli::genoma().to_vec();
        assert_eq!(genoma[mutacion.posicion], mutacion.referencia[0]);
        genoma[mutacion.posicion] = mutacion.alternativa[0];
        let tramo = &genoma[elegido.inicio - 1..elegido.fin];
        let adn = if elegido.hebra == b'-' {
            complemento_inverso(tramo)
        } else {
            tramo.to_vec()
        };
        let proteina = traducir(&adn);
        assert_eq!(proteina[codon - 1], nuevo);
        assert_eq!(ecoli::buscar_gen(nombre).unwrap().proteina[codon - 1], b'S');
    }
}

#[test]
fn lecturas_sin_mutaciones_coinciden_con_la_referencia() {
    let referencia = referencia_al_azar(20_000, 3);
    let opciones = OpcionesSimulacion {
        cobertura: 5.0,
        calidad_inicial: 40.0,
        calidad_final: 40.0,
        fraccion_con_adaptador: 0.2,
        semilla: 9,
        ..OpcionesSimulacion::default()
    };
    let simulacion = simular_lecturas(&referencia, "prueba", &[], &opciones).unwrap();
    let esperadas = (20_000.0 * 5.0 / 150.0) as usize;
    assert!(simulacion.lecturas.len().abs_diff(esperadas) <= 2);
    let mut diferencias = 0;
    let mut comparadas = 0;
    let mut recortadas = 0;
    for (lectura, registro) in simulacion.lecturas.iter().zip(&simulacion.verdad) {
        assert_eq!(lectura.secuencia.len(), 150);
        assert_eq!(lectura.calidad.len(), 150);
        if registro
            .cigar
            .iter()
            .any(|(_, operacion)| *operacion == b'S')
        {
            recortadas += 1;
        }
        let (d, c) = diferencias_con_referencia(registro, &referencia);
        diferencias += d;
        comparadas += c;
    }
    assert!((diferencias as f64) / (comparadas as f64) < 0.001);
    let fraccion = recortadas as f64 / simulacion.lecturas.len() as f64;
    assert!(
        (0.15..0.25).contains(&fraccion),
        "fracción con adaptador {fraccion}"
    );
}

#[test]
fn las_mutaciones_aparecen_en_las_lecturas_y_en_el_cigar() {
    let referencia = referencia_al_azar(10_000, 4);
    let mutaciones = vec![
        Mutacion {
            posicion: 2000,
            referencia: vec![referencia[2000]],
            alternativa: vec![if referencia[2000] == b'A' { b'C' } else { b'A' }],
            nombre: "SNP".to_string(),
        },
        Mutacion {
            posicion: 5000,
            referencia: referencia[5000..5003].to_vec(),
            alternativa: vec![referencia[5000]],
            nombre: "deleción".to_string(),
        },
        Mutacion {
            posicion: 8000,
            referencia: vec![referencia[8000]],
            alternativa: vec![referencia[8000], b'G', b'G'],
            nombre: "inserción".to_string(),
        },
    ];
    let opciones = OpcionesSimulacion {
        cobertura: 30.0,
        calidad_inicial: 40.0,
        calidad_final: 40.0,
        fraccion_con_adaptador: 0.0,
        semilla: 5,
        ..OpcionesSimulacion::default()
    };
    let simulacion = simular_lecturas(&referencia, "prueba", &mutaciones, &opciones).unwrap();
    let mut con_alternativa = 0;
    let mut sobre_snp = 0;
    let mut con_deleciones = 0;
    let mut con_inserciones = 0;
    for registro in &simulacion.verdad {
        let inicio = registro.posicion - 1;
        let largo: usize = registro
            .cigar
            .iter()
            .filter(|(_, o)| matches!(o, b'M' | b'D'))
            .map(|(l, _)| *l as usize)
            .sum();
        if inicio <= 2000 && 2000 < inicio + largo {
            sobre_snp += 1;
            let indice = 2000 - inicio;
            if registro.secuencia[indice] == mutaciones[0].alternativa[0] {
                con_alternativa += 1;
            }
        }
        if registro.cigar.iter().any(|(_, o)| *o == b'D') {
            con_deleciones += 1;
        }
        if registro.cigar.iter().any(|(_, o)| *o == b'I') {
            con_inserciones += 1;
        }
        let (d, c) = diferencias_con_referencia(registro, &referencia);
        assert!(
            d * 20 < c,
            "demasiadas diferencias en {:?}",
            String::from_utf8_lossy(&registro.nombre)
        );
    }
    assert!(sobre_snp > 10);
    assert_eq!(con_alternativa, sobre_snp);
    assert!(con_deleciones > 10);
    assert!(con_inserciones > 10);
}

#[test]
fn la_misma_semilla_da_lo_mismo() {
    let referencia = referencia_al_azar(5_000, 1);
    let opciones = OpcionesSimulacion::default();
    let a = simular_lecturas(&referencia, "x", &[], &opciones).unwrap();
    let b = simular_lecturas(&referencia, "x", &[], &opciones).unwrap();
    assert_eq!(a.lecturas, b.lecturas);
}

#[test]
fn rechaza_opciones_imposibles() {
    let referencia = referencia_al_azar(1_000, 1);
    let mut opciones = OpcionesSimulacion {
        longitud_lectura: 10,
        ..OpcionesSimulacion::default()
    };
    assert!(simular_lecturas(&referencia, "x", &[], &opciones).is_err());
    opciones.longitud_lectura = 600;
    assert!(simular_lecturas(&referencia, "x", &[], &opciones).is_err());
    opciones.longitud_lectura = 100;
    opciones.cobertura = 1000.0;
    assert!(simular_lecturas(&referencia, "x", &[], &opciones).is_err());
}

#[test]
fn muestra_de_ejemplo_cubre_las_regiones() {
    let muestra = muestra_de_ejemplo().unwrap();
    assert!(muestra.lecturas.len() > 35_000);
    assert_eq!(muestra.mutaciones.len(), 26);
    for registro in muestra.verdad.iter().take(2000) {
        let dentro = REGIONES_DE_EJEMPLO
            .iter()
            .any(|(inicio, fin)| registro.posicion > *inicio && registro.posicion <= *fin);
        assert!(dentro);
    }
}
