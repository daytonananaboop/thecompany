import { invocar } from "./puente.js";
import { el, formatear } from "./bloques.js";
import { pintarCodigo } from "./codigo.js";
import { DIBUJOS, PIZARRA } from "./dibujos.js";

const VELOCIDADES = [
  ["Lento", 1400],
  ["Normal", 700],
  ["Rápido", 250],
];

function dibujoGenerico(lienzo, estado) {
  const filas = [];
  for (const [clave, valor] of Object.entries(estado || {})) {
    let contenido;
    if (Array.isArray(valor) && valor.length && Array.isArray(valor[0])) {
      contenido = el("table", { class: "mini-tabla" }, valor.slice(0, 30).map((fila) => el("tr", {}, fila.slice(0, 30).map((c) => el("td", { texto: c === null ? "" : typeof c === "number" ? formatear(c) : String(c) })))));
    } else if (typeof valor === "object" && valor !== null) {
      contenido = el("code", { texto: JSON.stringify(valor).slice(0, 400) });
    } else {
      contenido = el("code", { texto: typeof valor === "number" ? formatear(valor) : String(valor) });
    }
    filas.push(el("div", { class: "valor" }, [el("span", { class: "v-etiqueta", texto: clave.replaceAll("_", " ") }), contenido]));
  }
  lienzo.replaceChildren(el("div", { class: "valores" }, filas));
}

export async function abrirPasoAPaso(contenedor, ficha, entrada, fragmentos) {
  contenedor.hidden = false;
  contenedor.replaceChildren(el("p", { class: "nota-paso", texto: "Preparando el paso a paso…" }));
  let lista;
  try {
    lista = await invocar("pasos", { id: ficha.id, entrada });
  } catch (error) {
    contenedor.replaceChildren(el("p", { class: "problema", texto: String(error) }));
    return;
  }
  if (!lista.length) {
    contenedor.replaceChildren(el("p", { class: "problema", texto: "El algoritmo no produjo pasos con esta entrada." }));
    return;
  }
  const dibujar = DIBUJOS[ficha.id] || dibujoGenerico;
  const lienzo = el("div", { class: "lienzo" });
  const { caja, lineas } = pintarCodigo(fragmentos);
  const panelCodigo = el("div", { class: "codigo-paso" }, [caja]);
  const nota = el("p", { class: "nota-paso", "aria-live": "polite" });
  const contador = el("span", { class: "contador" });
  const regla = el("input", { type: "range", min: 0, max: lista.length - 1, value: 0, "aria-label": "Paso" });
  const reproducir = el("button", { type: "button", class: "boton", texto: "Reproducir" });
  const velocidad = el("select", { "aria-label": "Velocidad" }, VELOCIDADES.map(([nombre, ms], i) => el("option", { value: ms, selected: i === 1, texto: nombre })));
  const boton = (texto, accion, etiqueta) => {
    const b = el("button", { type: "button", class: "boton", "aria-label": etiqueta, texto });
    b.addEventListener("click", accion);
    return b;
  };
  let actual = 0;
  let reloj = null;
  const detener = () => {
    clearInterval(reloj);
    reloj = null;
    reproducir.textContent = "Reproducir";
  };
  const ir = (indice) => {
    actual = Math.max(0, Math.min(lista.length - 1, indice));
    const paso = lista[actual];
    try {
      const propio = paso.estado && paso.estado.tipo === "pizarra" ? PIZARRA : dibujar;
      propio(lienzo, paso.estado, actual, lista);
    } catch (error) {
      dibujoGenerico(lienzo, paso.estado);
    }
    nota.textContent = paso.nota;
    contador.textContent = `Paso ${actual + 1} de ${lista.length}`;
    regla.value = actual;
    for (const linea of lineas) {
      linea.elemento.classList.remove("activa");
    }
    const encontrada = lineas.find((linea) => linea.texto.includes(paso.linea));
    if (encontrada) {
      encontrada.elemento.classList.add("activa");
      const arriba = encontrada.elemento.offsetTop - panelCodigo.offsetTop;
      if (arriba < panelCodigo.scrollTop || arriba > panelCodigo.scrollTop + panelCodigo.clientHeight - 40) {
        panelCodigo.scrollTop = Math.max(0, arriba - panelCodigo.clientHeight / 3);
      }
    }
  };
  reproducir.addEventListener("click", () => {
    if (reloj) {
      detener();
      return;
    }
    if (actual >= lista.length - 1) {
      ir(0);
    }
    reproducir.textContent = "Pausa";
    reloj = setInterval(() => {
      if (actual >= lista.length - 1) {
        detener();
        return;
      }
      ir(actual + 1);
    }, Number(velocidad.value));
  });
  velocidad.addEventListener("change", () => {
    if (reloj) {
      detener();
      reproducir.click();
    }
  });
  regla.addEventListener("input", () => {
    detener();
    ir(Number(regla.value));
  });
  const controles = el("div", { class: "controles" }, [
    boton("Inicio", () => { detener(); ir(0); }, "Primer paso"),
    boton("Atrás", () => { detener(); ir(actual - 1); }, "Paso anterior"),
    reproducir,
    boton("Adelante", () => { detener(); ir(actual + 1); }, "Paso siguiente"),
    boton("Final", () => { detener(); ir(lista.length - 1); }, "Último paso"),
    regla,
    velocidad,
    contador,
  ]);
  const cerrar = boton("Cerrar", () => { detener(); contenedor.hidden = true; contenedor.replaceChildren(); }, "Cerrar el paso a paso");
  const marco = el("div", { class: "seccion" }, [
    el("div", { class: "acciones" }, [el("h3", { class: "oculto-visual", texto: "Paso a paso" }), cerrar]),
    el("div", { class: "lienzo-y-codigo" }, [lienzo, panelCodigo]),
    nota,
    controles,
  ]);
  contenedor.replaceChildren(marco);
  marco.addEventListener("keydown", (evento) => {
    if (evento.target.tagName === "SELECT" || evento.target.tagName === "INPUT") {
      return;
    }
    if (evento.key === "ArrowRight") {
      detener();
      ir(actual + 1);
    } else if (evento.key === "ArrowLeft") {
      detener();
      ir(actual - 1);
    }
  });
  ir(0);
  contenedor.scrollIntoView({ behavior: "smooth", block: "start" });
}
