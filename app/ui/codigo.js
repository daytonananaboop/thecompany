import { el } from "./bloques.js";

const CLAVES = new Set([
  "as", "break", "const", "continue", "crate", "else", "enum", "false", "fn", "for", "if", "impl", "in",
  "let", "loop", "match", "mod", "move", "mut", "pub", "ref", "return", "self", "Self", "static", "struct",
  "super", "trait", "true", "type", "use", "where", "while",
]);

function pintarLinea(texto) {
  const linea = el("span", { class: "l" });
  let i = 0;
  let suelto = "";
  const soltar = () => {
    if (suelto) {
      linea.append(document.createTextNode(suelto));
      suelto = "";
    }
  };
  while (i < texto.length) {
    const letra = texto[i];
    if (letra === '"') {
      let fin = i + 1;
      while (fin < texto.length && texto[fin] !== '"') {
        fin += texto[fin] === "\\" ? 2 : 1;
      }
      soltar();
      linea.append(el("span", { class: "rs-texto", texto: texto.slice(i, fin + 1) }));
      i = fin + 1;
      continue;
    }
    if (letra === "'" && i + 2 < texto.length && (texto[i + 2] === "'" || (texto[i + 1] === "\\" && texto[i + 3] === "'"))) {
      const fin = texto[i + 1] === "\\" ? i + 4 : i + 3;
      soltar();
      linea.append(el("span", { class: "rs-texto", texto: texto.slice(i, fin) }));
      i = fin;
      continue;
    }
    if (/[0-9]/.test(letra) && (i === 0 || !/[A-Za-z_]/.test(texto[i - 1]))) {
      let fin = i;
      while (fin < texto.length && /[0-9_.a-fA-Fxu]/.test(texto[fin]) && !(texto[fin] === "." && texto[fin + 1] === ".")) {
        fin += 1;
      }
      soltar();
      linea.append(el("span", { class: "rs-numero", texto: texto.slice(i, fin) }));
      i = fin;
      continue;
    }
    if (/[A-Za-z_]/.test(letra)) {
      let fin = i;
      while (fin < texto.length && /[A-Za-z0-9_]/.test(texto[fin])) {
        fin += 1;
      }
      const palabra = texto.slice(i, fin);
      if (CLAVES.has(palabra)) {
        soltar();
        linea.append(el("span", { class: "rs-clave", texto: palabra }));
      } else if (/^[A-Z]/.test(palabra)) {
        soltar();
        linea.append(el("span", { class: "rs-tipo", texto: palabra }));
      } else {
        suelto += palabra;
      }
      i = fin;
      continue;
    }
    suelto += letra;
    i += 1;
  }
  soltar();
  if (!linea.childNodes.length) {
    linea.append(document.createTextNode(" "));
  }
  return linea;
}

export function pintarCodigo(fragmentos) {
  const lineas = [];
  const caja = el("div", { class: "fragmentos" });
  for (const fragmento of fragmentos) {
    const pre = el("pre", { class: "rust", tabindex: 0 });
    for (const texto of fragmento.fuente.split("\n")) {
      const linea = pintarLinea(texto);
      lineas.push({ elemento: linea, texto, fragmento: fragmento.nombre });
      pre.append(linea);
    }
    caja.append(el("div", { class: "fragmento" }, [el("h4", { texto: fragmento.nombre }), pre]));
  }
  return { caja, lineas };
}
