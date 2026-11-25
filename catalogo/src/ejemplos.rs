use nucleo::formatos::{
    Variante, escribir_cabecera_sam, escribir_cabecera_vcf, escribir_fastq, escribir_sam,
    escribir_vcf,
};
use nucleo::simulacion::{self, Simulacion};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicU64, Ordering};

static CANDADO: Mutex<()> = Mutex::new(());
static CONTADOR: AtomicU64 = AtomicU64::new(0);

pub fn generar(nombre: &str) -> Result<PathBuf, String> {
    let carpeta = crate::carpeta_de_ejemplos()?;
    let lecturas = carpeta.join("ejemplo_v1_lecturas.fastq");
    let alineamientos = carpeta.join("ejemplo_v1_lecturas.verdad.sam");
    let mutaciones = carpeta.join("ejemplo_v1_mutaciones.vcf");
    let pedido = match nombre {
        "lecturas" => lecturas.clone(),
        "alineamientos" => alineamientos.clone(),
        "mutaciones" => mutaciones.clone(),
        _ => return Err(format!("No conozco el ejemplo «{nombre}».")),
    };
    let _turno = CANDADO
        .lock()
        .unwrap_or_else(|envenenado| envenenado.into_inner());
    if lecturas.is_file() && alineamientos.is_file() && mutaciones.is_file() {
        return Ok(pedido);
    }
    let muestra = simulacion::muestra_de_ejemplo()?;
    escribir_simulacion(
        &muestra,
        nucleo::ecoli::IDENTIFICADOR,
        nucleo::ecoli::genoma().len(),
        &lecturas,
        &alineamientos,
        &mutaciones,
    )?;
    Ok(pedido)
}

pub fn escribir_simulacion(
    muestra: &Simulacion,
    referencia: &str,
    longitud: usize,
    ruta_lecturas: &Path,
    ruta_alineamientos: &Path,
    ruta_mutaciones: &Path,
) -> Result<(), String> {
    escribir_seguro(ruta_lecturas, |destino| {
        for lectura in &muestra.lecturas {
            escribir_fastq(destino, lectura)?;
        }
        Ok(())
    })?;
    escribir_seguro(ruta_alineamientos, |destino| {
        escribir_cabecera_sam(destino, referencia, longitud)?;
        for registro in &muestra.verdad {
            escribir_sam(destino, registro)?;
        }
        Ok(())
    })?;
    escribir_seguro(ruta_mutaciones, |destino| {
        escribir_cabecera_vcf(destino, referencia, longitud)?;
        for mutacion in &muestra.mutaciones {
            let variante = Variante {
                cromosoma: referencia.to_string(),
                posicion: mutacion.posicion + 1,
                referencia: mutacion.referencia.clone(),
                alternativa: mutacion.alternativa.clone(),
                calidad: 60.0,
                profundidad: 0,
                frecuencia: 1.0,
            };
            escribir_vcf(destino, &variante)?;
        }
        Ok(())
    })
}

fn escribir_seguro(
    ruta: &Path,
    contenido: impl FnOnce(&mut BufWriter<File>) -> std::io::Result<()>,
) -> Result<(), String> {
    let numero = CONTADOR.fetch_add(1, Ordering::Relaxed);
    let temporal = ruta.with_extension(format!("tmp{}-{numero}", std::process::id()));
    let archivo = File::create(&temporal)
        .map_err(|error| format!("No se pudo crear {}: {error}", temporal.display()))?;
    let mut destino = BufWriter::new(archivo);
    contenido(&mut destino)
        .and_then(|_| destino.flush())
        .map_err(|error| format!("No se pudo escribir {}: {error}", temporal.display()))?;
    drop(destino);
    fs::rename(&temporal, ruta)
        .map_err(|error| format!("No se pudo guardar {}: {error}", ruta.display()))
}
