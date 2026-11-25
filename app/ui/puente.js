export function invocar(comando, argumentos = {}) {
  const tauri = window.__TAURI__;
  if (!tauri || !tauri.core) {
    return Promise.reject("La Compañía debe abrirse como aplicación de escritorio.");
  }
  return tauri.core.invoke(comando, argumentos);
}

export async function elegirArchivo() {
  const tauri = window.__TAURI__;
  if (!tauri || !tauri.dialog) {
    return null;
  }
  const elegido = await tauri.dialog.open({ multiple: false, directory: false });
  return typeof elegido === "string" ? elegido : null;
}

export async function guardarComo(nombre, contenido) {
  const tauri = window.__TAURI__;
  if (!tauri || !tauri.dialog) {
    return null;
  }
  const ruta = await tauri.dialog.save({ defaultPath: nombre });
  if (typeof ruta !== "string") {
    return null;
  }
  if (typeof contenido === "string") {
    await tauri.core.invoke("guardar_texto", { ruta, contenido });
  } else {
    await tauri.core.invoke("guardar_bytes", { ruta, datos: Array.from(contenido) });
  }
  return ruta;
}

export async function mostrarEnCarpeta(ruta) {
  const tauri = window.__TAURI__;
  if (tauri && tauri.opener) {
    await tauri.opener.revealItemInDir(ruta);
  }
}

export function recordar(clave, valor) {
  try {
    localStorage.setItem("la-compania." + clave, valor);
  } catch (error) {
    return;
  }
}

export function recuperar(clave) {
  try {
    return localStorage.getItem("la-compania." + clave);
  } catch (error) {
    return null;
  }
}
