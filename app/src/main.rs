#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use serde_json::{Value, json};

#[tauri::command]
fn capitulos() -> Value {
    json!(catalogo::capitulos())
}

#[tauri::command]
fn fichas() -> Value {
    json!(catalogo::fichas())
}

#[tauri::command]
fn genes_ecoli() -> Value {
    let lista: Vec<Value> = nucleo::ecoli::genes()
        .iter()
        .map(|elegido| json!([elegido.nombre, elegido.producto]))
        .collect();
    Value::Array(lista)
}

#[tauri::command]
fn nucleos() -> usize {
    catalogo::nucleos()
}

#[tauri::command]
fn guardar_texto(ruta: String, contenido: String) -> Result<(), String> {
    std::fs::write(&ruta, contenido).map_err(|error| format!("No se pudo guardar {ruta}: {error}"))
}

#[tauri::command]
fn guardar_bytes(ruta: String, datos: Vec<u8>) -> Result<(), String> {
    std::fs::write(&ruta, datos).map_err(|error| format!("No se pudo guardar {ruta}: {error}"))
}

#[tauri::command]
fn ejemplo(id: String) -> Result<Value, String> {
    catalogo::ejemplo(&id)
}

#[tauri::command]
fn codigo(id: String) -> Result<Value, String> {
    catalogo::codigo(&id).map(|fragmentos| json!(fragmentos))
}

#[tauri::command]
async fn ejecutar(id: String, entrada: Value) -> Result<Value, String> {
    tauri::async_runtime::spawn_blocking(move || {
        catalogo::ejecutar(&id, entrada).map(|salida| json!(salida))
    })
    .await
    .map_err(|error| format!("No se pudo ejecutar: {error}"))?
}

#[tauri::command]
async fn pasos(id: String, entrada: Value) -> Result<Value, String> {
    tauri::async_runtime::spawn_blocking(move || {
        catalogo::pasos(&id, entrada).map(|lista| json!(lista))
    })
    .await
    .map_err(|error| format!("No se pudo preparar el paso a paso: {error}"))?
}

fn main() {
    std::thread::spawn(catalogo::preparar);
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_opener::init())
        .invoke_handler(tauri::generate_handler![
            capitulos,
            fichas,
            genes_ecoli,
            nucleos,
            guardar_texto,
            guardar_bytes,
            ejemplo,
            codigo,
            ejecutar,
            pasos
        ])
        .run(tauri::generate_context!())
        .expect("no se pudo abrir la ventana");
}
