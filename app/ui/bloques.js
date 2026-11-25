import { mostrarEnCarpeta, guardarComo } from "./puente.js";

const SVG = "http://www.w3.org/2000/svg";
const FORMAS = ["circulo", "cuadro", "triangulo", "rombo", "circulo", "cuadro", "triangulo", "rombo"];

export function el(etiqueta, atributos = {}, hijos = []) {
  const elemento = document.createElement(etiqueta);
  for (const [clave, valor] of Object.entries(atributos)) {
    if (clave === "texto") {
      elemento.textContent = valor;
    } else if (clave.startsWith("al")) {
      elemento.addEventListener(clave.slice(2).toLowerCase(), valor);
    } else if (valor !== undefined && valor !== null && valor !== false) {
      elemento.setAttribute(clave, valor === true ? "" : valor);
    }
  }
  for (const hijo of hijos) {
    if (hijo !== null && hijo !== undefined) {
      elemento.append(hijo);
    }
  }
  return elemento;
}

export function sv(etiqueta, atributos = {}, hijos = []) {
  const elemento = document.createElementNS(SVG, etiqueta);
  for (const [clave, valor] of Object.entries(atributos)) {
    if (clave === "texto") {
      elemento.textContent = valor;
    } else if (valor !== undefined && valor !== null) {
      elemento.setAttribute(clave, valor);
    }
  }
  for (const hijo of hijos) {
    if (hijo) {
      elemento.append(hijo);
    }
  }
  return elemento;
}

export function color(variable) {
  return getComputedStyle(document.documentElement).getPropertyValue(variable).trim();
}

export function colorSerie(indice) {
  if (indice < 0) {
    return color("--eje");
  }
  return color(`--serie-${(indice % 8) + 1}`);
}

export function formatear(numero) {
  if (typeof numero !== "number" || !Number.isFinite(numero)) {
    return String(numero);
  }
  const absoluto = Math.abs(numero);
  if (absoluto !== 0 && (absoluto < 0.001 || absoluto >= 1e9)) {
    return numero.toExponential(2);
  }
  const cifras = absoluto >= 100 ? 0 : absoluto >= 10 ? 1 : absoluto >= 1 ? 2 : 3;
  const [entera, decimal] = numero.toFixed(cifras).split(".");
  const signo = entera.startsWith("-") ? "-" : "";
  const digitos = signo ? entera.slice(1) : entera;
  let agrupado = "";
  for (let i = 0; i < digitos.length; i += 1) {
    if (i > 0 && (digitos.length - i) % 3 === 0) {
      agrupado += " ";
    }
    agrupado += digitos[i];
  }
  const limpio = decimal ? decimal.replace(/0+$/, "") : "";
  return signo + agrupado + (limpio ? "." + limpio : "");
}

const punta = () => document.getElementById("punta");

export function mostrarPunta(evento, contenido) {
  const caja = punta();
  caja.textContent = contenido;
  caja.hidden = false;
  const margen = 14;
  let x = evento.clientX + margen;
  let y = evento.clientY + margen;
  const ancho = caja.offsetWidth;
  const alto = caja.offsetHeight;
  if (x + ancho > window.innerWidth - 8) {
    x = evento.clientX - ancho - margen;
  }
  if (y + alto > window.innerHeight - 8) {
    y = evento.clientY - alto - margen;
  }
  caja.style.left = x + "px";
  caja.style.top = y + "px";
}

export function ocultarPunta() {
  punta().hidden = true;
}

async function copiar(texto, boton) {
  try {
    await navigator.clipboard.writeText(texto);
  } catch (error) {
    const area = el("textarea", {}, []);
    area.value = texto;
    document.body.append(area);
    area.select();
    document.execCommand("copy");
    area.remove();
  }
  const original = boton.textContent;
  boton.textContent = "Copiado";
  setTimeout(() => {
    boton.textContent = original;
  }, 1400);
}

function botonCopiar(obtener) {
  const boton = el("button", { type: "button", class: "enlace", texto: "Copiar" });
  boton.addEventListener("click", () => copiar(obtener(), boton));
  return boton;
}

function figura(etiqueta, nota, cuerpo, extra = []) {
  return el("figure", { class: "bloque" }, [
    el("figcaption", {}, [el("span", { texto: etiqueta }), nota ? el("span", { class: "nota-bloque", texto: nota }) : null, ...extra]),
    ...cuerpo,
  ]);
}

function pasos(minimo, maximo, cantidad) {
  if (!Number.isFinite(minimo) || !Number.isFinite(maximo)) {
    return [0];
  }
  if (minimo === maximo) {
    return [minimo];
  }
  const crudo = (maximo - minimo) / cantidad;
  const magnitud = Math.pow(10, Math.floor(Math.log10(crudo)));
  const normal = crudo / magnitud;
  const paso = (normal < 1.5 ? 1 : normal < 3 ? 2 : normal < 7 ? 5 : 10) * magnitud;
  const lista = [];
  for (let v = Math.ceil(minimo / paso) * paso; v <= maximo + paso * 1e-9; v += paso) {
    lista.push(Math.abs(v) < paso * 1e-9 ? 0 : v);
  }
  return lista;
}

function escala(dominio, rango) {
  const [d0, d1] = dominio;
  const [r0, r1] = rango;
  const ancho = d1 - d0 || 1;
  return (valor) => r0 + ((valor - d0) / ancho) * (r1 - r0);
}

function ejes(grupo, x, y, dominioX, dominioY, area, etiquetas) {
  const rejilla = sv("g", { class: "rejilla" });
  for (const valor of pasos(dominioY[0], dominioY[1], 5)) {
    const yy = y(valor);
    rejilla.append(sv("line", { x1: area.izq, x2: area.der, y1: yy, y2: yy }));
    grupo.append(sv("text", { x: area.izq - 8, y: yy + 4, "text-anchor": "end", texto: formatear(valor) }));
  }
  grupo.prepend(rejilla);
  if (dominioX) {
    for (const valor of pasos(dominioX[0], dominioX[1], 6)) {
      const xx = x(valor);
      grupo.append(sv("line", { class: "eje", x1: xx, x2: xx, y1: area.abajo, y2: area.abajo + 5 }));
      grupo.append(sv("text", { x: xx, y: area.abajo + 18, "text-anchor": "middle", texto: formatear(valor) }));
    }
  }
  grupo.append(sv("line", { class: "eje", x1: area.izq, x2: area.der, y1: area.abajo, y2: area.abajo }));
  grupo.append(sv("line", { class: "eje", x1: area.izq, x2: area.izq, y1: area.arriba, y2: area.abajo }));
  if (etiquetas.x) {
    grupo.append(sv("text", { class: "titulo-eje", x: (area.izq + area.der) / 2, y: area.abajo + 36, "text-anchor": "middle", texto: etiquetas.x }));
  }
  if (etiquetas.y) {
    const cy = (area.arriba + area.abajo) / 2;
    grupo.append(sv("text", { class: "titulo-eje", x: 14, y: cy, "text-anchor": "middle", transform: `rotate(-90 14 ${cy})`, texto: etiquetas.y }));
  }
}

function rango(valores) {
  let minimo = Infinity;
  let maximo = -Infinity;
  for (const v of valores) {
    if (Number.isFinite(v)) {
      minimo = Math.min(minimo, v);
      maximo = Math.max(maximo, v);
    }
  }
  if (!Number.isFinite(minimo)) {
    return [0, 1];
  }
  if (minimo === maximo) {
    return [minimo - 1, maximo + 1];
  }
  return [minimo, maximo];
}

export function graficaLinea(bloque, opciones = {}) {
  const contenedor = el("div", { class: "grafica" });
  const completo = rango((bloque.puntos || []).map((p) => p[0]));
  const dibujar = (dominioX) => {
    const ancho = 760;
    const alto = opciones.alto || 280;
    const area = { izq: 66, der: ancho - 16, arriba: 18, abajo: alto - 46 };
    const puntos = (bloque.puntos || []).filter((p) => p[0] >= dominioX[0] && p[0] <= dominioX[1]);
    const dominioY = rango(puntos.map((p) => p[1]));
    const margenY = (dominioY[1] - dominioY[0]) * 0.06;
    dominioY[0] -= margenY;
    dominioY[1] += margenY;
    const x = escala(dominioX, [area.izq, area.der]);
    const y = escala(dominioY, [area.abajo, area.arriba]);
    const lienzo = sv("svg", { viewBox: `0 0 ${ancho} ${alto}`, role: "img", "aria-label": bloque.etiqueta });
    const grupo = sv("g");
    lienzo.append(grupo);
    ejes(grupo, x, y, dominioX, dominioY, area, { x: bloque.eje_x, y: bloque.eje_y });
    for (const marca of bloque.marcas || []) {
      if (marca.x < dominioX[0] || marca.x > dominioX[1]) {
        continue;
      }
      const xx = x(marca.x);
      grupo.append(sv("line", { x1: xx, x2: xx, y1: area.arriba, y2: area.abajo, stroke: color("--rojo"), "stroke-dasharray": "4 4", "stroke-width": 1.5 }));
      grupo.append(sv("text", { x: Math.min(xx + 4, area.der - 4), y: area.arriba + 10, "text-anchor": xx > area.der - 120 ? "end" : "start", fill: color("--rojo"), texto: marca.texto }));
    }
    if (puntos.length) {
      const d = puntos.map((p, i) => (i ? "L" : "M") + x(p[0]).toFixed(1) + " " + y(p[1]).toFixed(1)).join(" ");
      grupo.append(sv("path", { d, fill: "none", stroke: color("--serie-1"), "stroke-width": 2, "stroke-linejoin": "round", "stroke-linecap": "round" }));
    }
    const cruz = sv("line", { y1: area.arriba, y2: area.abajo, stroke: color("--eje"), "stroke-width": 1, visibility: "hidden" });
    const marcador = sv("circle", { r: 4.5, fill: color("--serie-1"), stroke: color("--hoja"), "stroke-width": 2, visibility: "hidden" });
    const seleccion = sv("rect", { y: area.arriba, height: area.abajo - area.arriba, fill: color("--serie-1"), opacity: 0.15, visibility: "hidden" });
    const captura = sv("rect", { x: area.izq, y: area.arriba, width: area.der - area.izq, height: area.abajo - area.arriba, fill: "transparent", style: "cursor: crosshair" });
    grupo.append(seleccion, cruz, marcador, captura);
    const enDatos = (evento) => {
      const caja = lienzo.getBoundingClientRect();
      const xx = ((evento.clientX - caja.left) / caja.width) * ancho;
      return { xx, valor: dominioX[0] + ((xx - area.izq) / (area.der - area.izq)) * (dominioX[1] - dominioX[0]) };
    };
    let arrastre = null;
    captura.addEventListener("pointerdown", (evento) => {
      arrastre = enDatos(evento);
      captura.setPointerCapture(evento.pointerId);
    });
    captura.addEventListener("pointerup", (evento) => {
      if (!arrastre) {
        return;
      }
      const final = enDatos(evento);
      const desde = Math.min(arrastre.valor, final.valor);
      const hasta = Math.max(arrastre.valor, final.valor);
      const recorrido = Math.abs(final.xx - arrastre.xx);
      arrastre = null;
      seleccion.setAttribute("visibility", "hidden");
      if (recorrido > 8 && hasta > desde) {
        dibujar([desde, hasta]);
      }
    });
    captura.addEventListener("dblclick", () => dibujar(completo));
    captura.addEventListener("pointermove", (evento) => {
      const { xx, valor } = enDatos(evento);
      if (arrastre) {
        seleccion.setAttribute("x", Math.min(xx, arrastre.xx));
        seleccion.setAttribute("width", Math.abs(xx - arrastre.xx));
        seleccion.setAttribute("visibility", "visible");
      }
      if (!puntos.length) {
        return;
      }
      let bajo = 0;
      let alto2 = puntos.length - 1;
      while (bajo < alto2) {
        const medio = (bajo + alto2) >> 1;
        if (puntos[medio][0] < valor) {
          bajo = medio + 1;
        } else {
          alto2 = medio;
        }
      }
      let elegido = bajo;
      if (bajo > 0 && Math.abs(puntos[bajo - 1][0] - valor) < Math.abs(puntos[bajo][0] - valor)) {
        elegido = bajo - 1;
      }
      const p = puntos[elegido];
      cruz.setAttribute("x1", x(p[0]));
      cruz.setAttribute("x2", x(p[0]));
      marcador.setAttribute("cx", x(p[0]));
      marcador.setAttribute("cy", y(p[1]));
      cruz.setAttribute("visibility", "visible");
      marcador.setAttribute("visibility", "visible");
      mostrarPunta(evento, `${bloque.eje_x || "x"}: ${formatear(p[0])}\n${bloque.eje_y || "y"}: ${formatear(p[1])}`);
    });
    captura.addEventListener("pointerleave", () => {
      cruz.setAttribute("visibility", "hidden");
      marcador.setAttribute("visibility", "hidden");
      ocultarPunta();
    });
    const acercado = dominioX[0] > completo[0] || dominioX[1] < completo[1];
    contenedor.replaceChildren(lienzo, el("p", { class: "ayuda-grafica", texto: acercado ? "Doble clic para ver todo." : "Arrastra sobre la gráfica para acercar." }));
  };
  dibujar(completo);
  return contenedor;
}

const PROPIEDADES_VISIBLES = ["fill", "stroke", "stroke-width", "stroke-dasharray", "opacity", "font-family", "font-size", "font-weight", "text-anchor"];

function copiarEstilos(original, copia) {
  const estilo = getComputedStyle(original);
  const reglas = PROPIEDADES_VISIBLES.map((propiedad) => `${propiedad}:${estilo.getPropertyValue(propiedad)}`).join(";");
  copia.setAttribute("style", reglas);
  for (let i = 0; i < original.children.length; i += 1) {
    copiarEstilos(original.children[i], copia.children[i]);
  }
}

function svgIndependiente(svg) {
  const copia = svg.cloneNode(true);
  copiarEstilos(svg, copia);
  copia.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  const caja = svg.viewBox.baseVal;
  const fondo = sv("rect", { x: caja.x, y: caja.y, width: caja.width, height: caja.height, fill: color("--hoja") });
  copia.insertBefore(fondo, copia.firstChild);
  return { texto: new XMLSerializer().serializeToString(copia), ancho: caja.width, alto: caja.height };
}

function nombreDeArchivo(etiqueta, extension) {
  const base = etiqueta.normalize("NFD").replace(/[\u0300-\u036f]/g, "").toLowerCase().replace(/[^a-z0-9]+/g, "-").replace(/^-|-$/g, "") || "resultado";
  return `${base}.${extension}`;
}

async function svgComoPng(svg) {
  const { texto, ancho, alto } = svgIndependiente(svg);
  const imagen = new Image();
  const url = URL.createObjectURL(new Blob([texto], { type: "image/svg+xml" }));
  await new Promise((resolver, rechazar) => {
    imagen.onload = resolver;
    imagen.onerror = rechazar;
    imagen.src = url;
  });
  const escalaPng = 2;
  const lienzo = document.createElement("canvas");
  lienzo.width = ancho * escalaPng;
  lienzo.height = alto * escalaPng;
  const contexto = lienzo.getContext("2d");
  contexto.scale(escalaPng, escalaPng);
  contexto.drawImage(imagen, 0, 0, ancho, alto);
  URL.revokeObjectURL(url);
  const blob = await new Promise((resolver) => lienzo.toBlob(resolver, "image/png"));
  return new Uint8Array(await blob.arrayBuffer());
}

function botonesDeExportar(etiqueta, contenedor) {
  const svg = () => contenedor.querySelector("svg");
  const boton = (texto, accion) => {
    const b = el("button", { type: "button", class: "enlace", texto });
    b.addEventListener("click", async () => {
      try {
        const ruta = await accion();
        if (ruta) {
          const original = b.textContent;
          b.textContent = "Guardado";
          setTimeout(() => {
            b.textContent = original;
          }, 1400);
        }
      } catch (error) {
        b.textContent = "No se pudo guardar";
      }
    });
    return b;
  };
  return [
    boton("Guardar SVG", () => guardarComo(nombreDeArchivo(etiqueta, "svg"), svgIndependiente(svg()).texto)),
    boton("Guardar PNG", async () => guardarComo(nombreDeArchivo(etiqueta, "png"), await svgComoPng(svg()))),
  ];
}

function figuraGrafica(etiqueta, nota, grafica) {
  return figura(etiqueta, nota, [grafica], botonesDeExportar(etiqueta, grafica));
}

function csv(columnas, filas) {
  const celda = (valor) => {
    const texto = String(valor ?? "");
    return /[",\n]/.test(texto) ? `"${texto.replaceAll('"', '""')}"` : texto;
  };
  return [columnas, ...filas].map((fila) => fila.map(celda).join(",")).join("\n") + "\n";
}

function barra(x, y0, y1, ancho) {
  const arriba = Math.min(y0, y1);
  const abajo = Math.max(y0, y1);
  const alto = abajo - arriba;
  const radio = Math.min(4, ancho / 2, alto);
  if (y1 <= y0) {
    return `M${x} ${abajo} L${x} ${arriba + radio} Q${x} ${arriba} ${x + radio} ${arriba} L${x + ancho - radio} ${arriba} Q${x + ancho} ${arriba} ${x + ancho} ${arriba + radio} L${x + ancho} ${abajo} Z`;
  }
  return `M${x} ${arriba} L${x + ancho} ${arriba} L${x + ancho} ${abajo - radio} Q${x + ancho} ${abajo} ${x + ancho - radio} ${abajo} L${x + radio} ${abajo} Q${x} ${abajo} ${x} ${abajo - radio} Z`;
}

export function graficaBarras(bloque) {
  const ancho = 760;
  const alto = 280;
  const categorias = bloque.categorias || [];
  const valores = bloque.valores || [];
  const largaEtiqueta = Math.max(0, ...categorias.map((c) => String(c).length));
  const giradas = categorias.length > 12 || largaEtiqueta > 8;
  const area = { izq: 66, der: ancho - 16, arriba: 18, abajo: alto - (giradas ? 70 : 46) };
  const dominioY = rango([0, ...valores]);
  const x0 = area.izq;
  const paso = (area.der - area.izq) / Math.max(1, categorias.length);
  const anchoBarra = Math.max(1, paso - 2);
  const y = escala(dominioY, [area.abajo, area.arriba]);
  const lienzo = sv("svg", { viewBox: `0 0 ${ancho} ${alto}`, role: "img", "aria-label": bloque.etiqueta });
  const grupo = sv("g");
  lienzo.append(grupo);
  ejes(grupo, null, y, null, dominioY, area, { x: bloque.eje_x, y: bloque.eje_y });
  const cada = Math.max(1, Math.ceil(categorias.length / 40));
  categorias.forEach((categoria, i) => {
    const valor = valores[i] ?? 0;
    const forma = sv("path", { d: barra(x0 + i * paso + 1, y(0), y(valor), anchoBarra), fill: color("--serie-1") });
    const zona = sv("rect", { x: x0 + i * paso, y: area.arriba, width: paso, height: area.abajo - area.arriba, fill: "transparent" });
    zona.addEventListener("pointermove", (evento) => mostrarPunta(evento, `${categoria}\n${bloque.eje_y || "valor"}: ${formatear(valor)}`));
    zona.addEventListener("pointerleave", ocultarPunta);
    grupo.append(forma, zona);
    if (i % cada === 0) {
      const cx = x0 + i * paso + paso / 2;
      const etiqueta = String(categoria).length > 14 ? String(categoria).slice(0, 13) + "…" : String(categoria);
      grupo.append(sv("text", giradas
        ? { x: cx, y: area.abajo + 10, "text-anchor": "end", transform: `rotate(-50 ${cx} ${area.abajo + 10})`, texto: etiqueta }
        : { x: cx, y: area.abajo + 18, "text-anchor": "middle", texto: etiqueta }));
    }
  });
  if (bloque.eje_x) {
    grupo.append(sv("text", { class: "titulo-eje", x: (area.izq + area.der) / 2, y: alto - 4, "text-anchor": "middle", texto: bloque.eje_x }));
  }
  return el("div", { class: "grafica" }, [lienzo]);
}

function marcadorForma(forma, cx, cy, r, relleno, borde) {
  if (forma === "cuadro") {
    return sv("rect", { x: cx - r, y: cy - r, width: r * 2, height: r * 2, fill: relleno, stroke: borde, "stroke-width": 1.5 });
  }
  if (forma === "triangulo") {
    return sv("path", { d: `M${cx} ${cy - r * 1.2} L${cx + r * 1.1} ${cy + r * 0.8} L${cx - r * 1.1} ${cy + r * 0.8} Z`, fill: relleno, stroke: borde, "stroke-width": 1.5 });
  }
  if (forma === "rombo") {
    return sv("path", { d: `M${cx} ${cy - r * 1.25} L${cx + r * 1.25} ${cy} L${cx} ${cy + r * 1.25} L${cx - r * 1.25} ${cy} Z`, fill: relleno, stroke: borde, "stroke-width": 1.5 });
  }
  return sv("circle", { cx, cy, r, fill: relleno, stroke: borde, "stroke-width": 1.5 });
}

export function leyendaGrupos(grupos) {
  return el("div", { class: "leyenda" }, grupos.map((g) => {
    const muestra = sv("svg", { viewBox: "0 0 12 12", "aria-hidden": "true" }, [marcadorForma(FORMAS[g % 8], 6, 6, 4.2, colorSerie(g), "none")]);
    return el("span", {}, [muestra, document.createTextNode(`Grupo ${g + 1}`)]);
  }));
}

export function graficaDispersion(bloque, opciones = {}) {
  const ancho = opciones.ancho || 760;
  const alto = opciones.alto || 340;
  const area = { izq: 66, der: ancho - 16, arriba: 16, abajo: alto - 46 };
  const puntos = bloque.puntos || [];
  const dominioX = rango(puntos.map((p) => p.x));
  const dominioY = rango(puntos.map((p) => p.y));
  for (const d of [dominioX, dominioY]) {
    const margen = (d[1] - d[0]) * 0.06;
    d[0] -= margen;
    d[1] += margen;
  }
  const x = escala(dominioX, [area.izq, area.der]);
  const y = escala(dominioY, [area.abajo, area.arriba]);
  const lienzo = sv("svg", { viewBox: `0 0 ${ancho} ${alto}`, role: "img", "aria-label": bloque.etiqueta });
  const grupo = sv("g");
  lienzo.append(grupo);
  ejes(grupo, x, y, dominioX, dominioY, area, { x: bloque.eje_x, y: bloque.eje_y });
  const grupos = [...new Set(puntos.map((p) => p.grupo))].filter((g) => g >= 0).sort((a, b) => a - b);
  for (const p of puntos) {
    grupo.append(marcadorForma(FORMAS[p.grupo % 8], x(p.x), y(p.y), 4.5, colorSerie(p.grupo), color("--hoja")));
  }
  for (const [anterior, centro] of opciones.movimientos || []) {
    grupo.append(sv("line", { x1: x(anterior.x), y1: y(anterior.y), x2: x(centro.x), y2: y(centro.y), stroke: color("--tinta"), "stroke-width": 1.5, "stroke-dasharray": "4 3" }));
    grupo.append(sv("circle", { cx: x(anterior.x), cy: y(anterior.y), r: 6, fill: "none", stroke: color("--eje"), "stroke-width": 1.5 }));
  }
  for (const centro of opciones.centros || []) {
    grupo.append(sv("path", { d: `M${x(centro.x) - 8} ${y(centro.y)} L${x(centro.x) + 8} ${y(centro.y)} M${x(centro.x)} ${y(centro.y) - 8} L${x(centro.x)} ${y(centro.y) + 8}`, stroke: colorSerie(centro.grupo), "stroke-width": 3 }));
    grupo.append(sv("circle", { cx: x(centro.x), cy: y(centro.y), r: 9, fill: "none", stroke: color("--tinta"), "stroke-width": 1.5 }));
  }
  const captura = sv("rect", { x: area.izq, y: area.arriba, width: area.der - area.izq, height: area.abajo - area.arriba, fill: "transparent" });
  grupo.append(captura);
  captura.addEventListener("pointermove", (evento) => {
    const caja = lienzo.getBoundingClientRect();
    const px = ((evento.clientX - caja.left) / caja.width) * ancho;
    const py = ((evento.clientY - caja.top) / caja.height) * alto;
    let mejor = null;
    let distancia = 14 * 14;
    for (const p of puntos) {
      const d = (x(p.x) - px) ** 2 + (y(p.y) - py) ** 2;
      if (d < distancia) {
        distancia = d;
        mejor = p;
      }
    }
    if (mejor) {
      mostrarPunta(evento, `${mejor.nombre || "punto"}${mejor.grupo >= 0 ? `\nGrupo ${mejor.grupo + 1}` : ""}\n(${formatear(mejor.x)}, ${formatear(mejor.y)})`);
    } else {
      ocultarPunta();
    }
  });
  captura.addEventListener("pointerleave", ocultarPunta);
  const partes = [lienzo];
  if (grupos.length > 1) {
    partes.unshift(leyendaGrupos(grupos));
  }
  return el("div", { class: "grafica" }, partes);
}

function mezclar(a, b, t) {
  const leer = (hex) => [1, 3, 5].map((i) => parseInt(hex.slice(i, i + 2), 16));
  const [ra, ga, ba] = leer(a);
  const [rb, gb, bb] = leer(b);
  const c = (u, v) => Math.round(u + (v - u) * t).toString(16).padStart(2, "0");
  return `#${c(ra, rb)}${c(ga, gb)}${c(ba, bb)}`;
}

export function colorCalor(valor, minimo, maximo) {
  const neutro = color("--neutro");
  if (minimo < 0 && maximo > 0) {
    const tope = Math.max(-minimo, maximo);
    const t = Math.min(1, Math.abs(valor) / tope);
    return mezclar(neutro, valor < 0 ? color("--serie-1") : color("--serie-2"), t);
  }
  const t = maximo === minimo ? 1 : (valor - minimo) / (maximo - minimo);
  return mezclar(neutro, color("--serie-1"), Math.max(0, Math.min(1, t)));
}

export function graficaCalor(bloque) {
  const filas = bloque.filas || [];
  const columnas = bloque.columnas || [];
  const valores = bloque.valores || [];
  const todos = valores.flat().filter(Number.isFinite);
  const [minimo, maximo] = todos.length ? [Math.min(...todos), Math.max(...todos)] : [0, 1];
  const izquierda = Math.min(150, 12 + Math.max(4, ...filas.map((f) => String(f).length)) * 7);
  const arriba = Math.min(90, 14 + Math.max(2, ...columnas.map((c) => String(c).length)) * 6);
  const ancho = 760;
  const celda = Math.max(4, Math.min(32, (ancho - izquierda - 10) / Math.max(1, columnas.length)));
  const altoCelda = Math.max(4, Math.min(22, 520 / Math.max(1, filas.length)));
  const alto = arriba + filas.length * altoCelda + 46;
  const lienzo = sv("svg", { viewBox: `0 0 ${ancho} ${alto}`, role: "img", "aria-label": bloque.etiqueta });
  columnas.forEach((columna, j) => {
    const cx = izquierda + j * celda + celda / 2;
    lienzo.append(sv("text", { x: cx, y: arriba - 6, "text-anchor": "start", transform: `rotate(-50 ${cx} ${arriba - 6})`, texto: String(columna) }));
  });
  filas.forEach((fila, i) => {
    if (altoCelda >= 9) {
      lienzo.append(sv("text", { x: izquierda - 6, y: arriba + i * altoCelda + altoCelda / 2 + 4, "text-anchor": "end", texto: String(fila) }));
    }
    (valores[i] || []).forEach((valor, j) => {
      const rect = sv("rect", { x: izquierda + j * celda + 1, y: arriba + i * altoCelda + 1, width: Math.max(1, celda - 2), height: Math.max(1, altoCelda - 2), fill: colorCalor(valor, minimo, maximo) });
      rect.addEventListener("pointermove", (evento) => mostrarPunta(evento, `${fila} × ${columnas[j]}\n${formatear(valor)}`));
      rect.addEventListener("pointerleave", ocultarPunta);
      lienzo.append(rect);
    });
  });
  const yLeyenda = arriba + filas.length * altoCelda + 18;
  const gradiente = sv("linearGradient", { id: `g${Math.random().toString(36).slice(2)}` });
  for (let k = 0; k <= 10; k += 1) {
    gradiente.append(sv("stop", { offset: k / 10, "stop-color": colorCalor(minimo + ((maximo - minimo) * k) / 10, minimo, maximo) }));
  }
  lienzo.append(sv("defs", {}, [gradiente]));
  lienzo.append(sv("rect", { x: izquierda, y: yLeyenda, width: 180, height: 10, fill: `url(#${gradiente.id})` }));
  lienzo.append(sv("text", { x: izquierda, y: yLeyenda + 24, texto: formatear(minimo) }));
  lienzo.append(sv("text", { x: izquierda + 180, y: yLeyenda + 24, "text-anchor": "end", texto: formatear(maximo) }));
  return el("div", { class: "grafica" }, [lienzo]);
}

export function leerNewick(texto) {
  let posicion = 0;
  const limpio = texto.trim();
  const leerNombre = () => {
    let nombre = "";
    if (limpio[posicion] === "'") {
      posicion += 1;
      while (posicion < limpio.length && limpio[posicion] !== "'") {
        nombre += limpio[posicion];
        posicion += 1;
      }
      posicion += 1;
      return nombre;
    }
    while (posicion < limpio.length && !"(),:;".includes(limpio[posicion])) {
      nombre += limpio[posicion];
      posicion += 1;
    }
    return nombre.trim();
  };
  const leerNodo = () => {
    const nodo = { nombre: "", largo: null, hijos: [] };
    if (limpio[posicion] === "(") {
      posicion += 1;
      nodo.hijos.push(leerNodo());
      while (limpio[posicion] === ",") {
        posicion += 1;
        nodo.hijos.push(leerNodo());
      }
      if (limpio[posicion] !== ")") {
        throw new Error("Newick sin cerrar");
      }
      posicion += 1;
    }
    nodo.nombre = leerNombre();
    if (limpio[posicion] === ":") {
      posicion += 1;
      const inicio = posicion;
      while (posicion < limpio.length && !"(),;".includes(limpio[posicion])) {
        posicion += 1;
      }
      const largo = parseFloat(limpio.slice(inicio, posicion));
      nodo.largo = Number.isFinite(largo) ? largo : null;
    }
    return nodo;
  };
  return leerNodo();
}

export function graficaArbol(newick) {
  let raiz;
  try {
    raiz = leerNewick(newick);
  } catch (error) {
    return el("pre", { class: "tipeado", texto: newick });
  }
  const hojas = [];
  const conLargos = { valor: false };
  const recorrer = (nodo, profundidad, distancia) => {
    nodo.profundidad = profundidad;
    if (nodo.largo !== null) {
      conLargos.valor = true;
    }
    nodo.distancia = distancia + (nodo.largo ?? 1);
    for (const hijo of nodo.hijos) {
      recorrer(hijo, profundidad + 1, nodo.distancia);
    }
    if (!nodo.hijos.length) {
      hojas.push(nodo);
    }
  };
  recorrer(raiz, 0, 0);
  raiz.distancia = 0;
  const reajustar = (nodo, distancia) => {
    nodo.x = conLargos.valor ? distancia : nodo.profundidad;
    for (const hijo of nodo.hijos) {
      reajustar(hijo, (conLargos.valor ? distancia + (hijo.largo ?? 0) : 0));
    }
  };
  reajustar(raiz, 0);
  if (!conLargos.valor) {
    const maxima = Math.max(...hojas.map((h) => h.profundidad));
    for (const hoja of hojas) {
      hoja.x = maxima;
    }
  }
  hojas.forEach((hoja, i) => {
    hoja.y = i;
  });
  const ubicar = (nodo) => {
    if (nodo.hijos.length) {
      nodo.hijos.forEach(ubicar);
      nodo.y = nodo.hijos.reduce((suma, h) => suma + h.y, 0) / nodo.hijos.length;
    }
  };
  ubicar(raiz);
  const etiquetaMaxima = Math.max(4, ...hojas.map((h) => h.nombre.length));
  const ancho = 760;
  const derecha = ancho - Math.min(220, etiquetaMaxima * 7.5 + 16);
  const paso = Math.max(16, Math.min(26, 520 / Math.max(1, hojas.length)));
  const alto = hojas.length * paso + 50;
  const todosX = [];
  const juntar = (nodo) => {
    todosX.push(nodo.x);
    nodo.hijos.forEach(juntar);
  };
  juntar(raiz);
  const maximoX = Math.max(...todosX) || 1;
  const x = (v) => 20 + (v / maximoX) * (derecha - 20);
  const y = (v) => 16 + v * paso;
  const lienzo = sv("svg", { viewBox: `0 0 ${ancho} ${alto}`, role: "img", "aria-label": "Árbol" });
  const tinta = color("--tinta");
  const dibujar = (nodo) => {
    if (nodo.hijos.length) {
      const ys = nodo.hijos.map((h) => y(h.y));
      lienzo.append(sv("line", { x1: x(nodo.x), x2: x(nodo.x), y1: Math.min(...ys), y2: Math.max(...ys), stroke: tinta, "stroke-width": 1.6 }));
      for (const hijo of nodo.hijos) {
        const linea = sv("line", { x1: x(nodo.x), x2: x(hijo.x), y1: y(hijo.y), y2: y(hijo.y), stroke: tinta, "stroke-width": 1.6 });
        if (hijo.largo !== null) {
          linea.addEventListener("pointermove", (evento) => mostrarPunta(evento, `largo de la rama: ${formatear(hijo.largo)}`));
          linea.addEventListener("pointerleave", ocultarPunta);
          linea.setAttribute("stroke-width", 2.4);
        }
        lienzo.append(linea);
        dibujar(hijo);
      }
      if (nodo.nombre && !/^[0-9.]+$/.test(nodo.nombre)) {
        lienzo.append(sv("text", { x: x(nodo.x) + 3, y: y(nodo.y) - 4, texto: nodo.nombre }));
      }
    } else {
      lienzo.append(sv("circle", { cx: x(nodo.x), cy: y(nodo.y), r: 2.5, fill: tinta }));
      lienzo.append(sv("text", { x: x(nodo.x) + 6, y: y(nodo.y) + 4, style: "font-family: var(--cuerpo); font-size: 13px; fill: var(--tinta)", texto: nodo.nombre }));
    }
  };
  dibujar(raiz);
  if (conLargos.valor) {
    const escalaValor = pasos(0, maximoX, 4)[1] || maximoX;
    const largoBarra = x(escalaValor) - x(0);
    lienzo.append(sv("line", { x1: 20, x2: 20 + largoBarra, y1: alto - 18, y2: alto - 18, stroke: tinta, "stroke-width": 1.6 }));
    lienzo.append(sv("text", { x: 20 + largoBarra + 6, y: alto - 14, texto: formatear(escalaValor) }));
  }
  return el("div", { class: "grafica" }, [lienzo]);
}

export function ubicarGrafo(cantidad, aristas, ancho, alto) {
  const posiciones = [];
  for (let i = 0; i < cantidad; i += 1) {
    const angulo = (2 * Math.PI * i) / Math.max(1, cantidad);
    posiciones.push({ x: ancho / 2 + Math.cos(angulo) * ancho * 0.35, y: alto / 2 + Math.sin(angulo) * alto * 0.35 });
  }
  if (cantidad < 3) {
    return posiciones;
  }
  const ideal = Math.sqrt((ancho * alto) / cantidad) * 0.75;
  let temperatura = ancho / 8;
  for (let vuelta = 0; vuelta < 300; vuelta += 1) {
    const fuerzas = posiciones.map(() => ({ x: 0, y: 0 }));
    for (let i = 0; i < cantidad; i += 1) {
      for (let j = i + 1; j < cantidad; j += 1) {
        const dx = posiciones[i].x - posiciones[j].x;
        const dy = posiciones[i].y - posiciones[j].y;
        const d = Math.max(0.01, Math.hypot(dx, dy));
        const empuje = (ideal * ideal) / d;
        fuerzas[i].x += (dx / d) * empuje;
        fuerzas[i].y += (dy / d) * empuje;
        fuerzas[j].x -= (dx / d) * empuje;
        fuerzas[j].y -= (dy / d) * empuje;
      }
    }
    for (const arista of aristas) {
      if (arista.desde === arista.hasta) {
        continue;
      }
      const a = posiciones[arista.desde];
      const b = posiciones[arista.hasta];
      if (!a || !b) {
        continue;
      }
      const dx = a.x - b.x;
      const dy = a.y - b.y;
      const d = Math.max(0.01, Math.hypot(dx, dy));
      const tiro = (d * d) / ideal;
      fuerzas[arista.desde].x -= (dx / d) * tiro;
      fuerzas[arista.desde].y -= (dy / d) * tiro;
      fuerzas[arista.hasta].x += (dx / d) * tiro;
      fuerzas[arista.hasta].y += (dy / d) * tiro;
    }
    posiciones.forEach((p, i) => {
      const f = fuerzas[i];
      const magnitud = Math.max(0.01, Math.hypot(f.x, f.y));
      p.x += (f.x / magnitud) * Math.min(magnitud, temperatura);
      p.y += (f.y / magnitud) * Math.min(magnitud, temperatura);
      p.x = Math.max(40, Math.min(ancho - 40, p.x));
      p.y = Math.max(24, Math.min(alto - 24, p.y));
    });
    temperatura *= 0.985;
  }
  return posiciones;
}

export function graficaGrafo(nodos, aristas, opciones = {}) {
  const ancho = opciones.ancho || 760;
  const alto = opciones.alto || Math.max(260, Math.min(520, 120 + nodos.length * 14));
  const posiciones = opciones.posiciones || ubicarGrafo(nodos.length, aristas, ancho, alto);
  const tinta = color("--tinta");
  const rojo = color("--rojo");
  const id = `flecha${Math.random().toString(36).slice(2)}`;
  const lienzo = sv("svg", { viewBox: `0 0 ${ancho} ${alto}`, role: "img", "aria-label": opciones.etiqueta || "Grafo" });
  lienzo.append(sv("defs", {}, [
    sv("marker", { id, viewBox: "0 0 10 10", refX: 10, refY: 5, markerWidth: 7, markerHeight: 7, orient: "auto-start-reverse" }, [sv("path", { d: "M0 0 L10 5 L0 10 Z", fill: tinta })]),
    sv("marker", { id: id + "r", viewBox: "0 0 10 10", refX: 10, refY: 5, markerWidth: 7, markerHeight: 7, orient: "auto-start-reverse" }, [sv("path", { d: "M0 0 L10 5 L0 10 Z", fill: rojo })]),
  ]));
  const anchoNodo = (texto) => Math.max(18, String(texto).length * 7 + 10);
  const vistas = new Map();
  aristas.forEach((arista, indice) => {
    const a = posiciones[arista.desde];
    const b = posiciones[arista.hasta];
    if (!a || !b) {
      return;
    }
    const resaltada = opciones.resaltadas && opciones.resaltadas.has(indice);
    const usada = opciones.usadas && opciones.usadas.has(indice);
    const trazo = resaltada ? rojo : tinta;
    const marca = resaltada ? `url(#${id}r)` : `url(#${id})`;
    const grosor = resaltada ? 3 : usada ? 2.2 : 1.3;
    const opacidad = opciones.usadas && !usada && !resaltada ? 0.45 : 1;
    if (arista.desde === arista.hasta) {
      const d = `M${a.x - 6} ${a.y - 10} C ${a.x - 30} ${a.y - 50}, ${a.x + 30} ${a.y - 50}, ${a.x + 6} ${a.y - 10}`;
      lienzo.append(sv("path", { d, fill: "none", stroke: trazo, "stroke-width": grosor, opacity: opacidad, "marker-end": marca }));
      return;
    }
    const clave = [Math.min(arista.desde, arista.hasta), Math.max(arista.desde, arista.hasta)].join("-");
    const repeticion = vistas.get(clave) || 0;
    vistas.set(clave, repeticion + 1);
    const dx = b.x - a.x;
    const dy = b.y - a.y;
    const d = Math.max(1, Math.hypot(dx, dy));
    const recorte = (anchoNodo(nodos[arista.hasta]) / 2) * Math.abs(dx / d) + 11 * Math.abs(dy / d) + 3;
    const salida = (anchoNodo(nodos[arista.desde]) / 2) * Math.abs(dx / d) + 11 * Math.abs(dy / d);
    const curva = (repeticion + (arista.desde > arista.hasta ? 1 : 0)) * 22;
    const mx = (a.x + b.x) / 2 - (dy / d) * curva;
    const my = (a.y + b.y) / 2 + (dx / d) * curva;
    const inicioX = a.x + (dx / d) * salida;
    const inicioY = a.y + (dy / d) * salida;
    const finX = b.x - (dx / d) * recorte;
    const finY = b.y - (dy / d) * recorte;
    lienzo.append(sv("path", { d: `M${inicioX} ${inicioY} Q${mx} ${my} ${finX} ${finY}`, fill: "none", stroke: trazo, "stroke-width": grosor, opacity: opacidad, "marker-end": marca }));
    if (arista.etiqueta && aristas.length <= 40) {
      lienzo.append(sv("text", { x: (inicioX + 2 * mx + finX) / 4, y: (inicioY + 2 * my + finY) / 4 - 4, "text-anchor": "middle", style: "font-size: 10px", texto: arista.etiqueta }));
    }
  });
  nodos.forEach((nombre, i) => {
    const p = posiciones[i];
    const ancho2 = anchoNodo(nombre);
    const activo = opciones.activos && opciones.activos.has(i);
    const grupo = sv("g", {}, [
      sv("rect", { x: p.x - ancho2 / 2, y: p.y - 11, width: ancho2, height: 22, rx: 3, fill: activo ? rojo : color("--hoja"), stroke: tinta, "stroke-width": 1.2 }),
      sv("text", { x: p.x, y: p.y + 4, "text-anchor": "middle", style: `fill: ${activo ? color("--hoja") : tinta}; font-size: 11px`, texto: String(nombre) }),
    ]);
    const entrada = aristas.filter((a) => a.hasta === i).length;
    const salidaN = aristas.filter((a) => a.desde === i).length;
    grupo.addEventListener("pointermove", (evento) => mostrarPunta(evento, `${nombre}\nentran ${entrada}, salen ${salidaN}`));
    grupo.addEventListener("pointerleave", ocultarPunta);
    lienzo.append(grupo);
  });
  return el("div", { class: "grafica" }, [lienzo]);
}

function esLineaMedia(secuencia) {
  return /^[|:+.* ]*$/.test(secuencia) && /[|:+.]/.test(secuencia);
}

export function vistaAlineamiento(filas) {
  const nombres = filas.map((f) => f.nombre);
  const ancho = Math.min(18, Math.max(4, ...nombres.map((n) => n.length)));
  const lineas = [];
  let media = null;
  let reales = filas;
  if (filas.length === 3 && esLineaMedia(filas[1].secuencia)) {
    media = filas[1].secuencia;
    reales = [filas[0], filas[2]];
  } else if (filas.length === 2) {
    const [a, b] = filas.map((f) => f.secuencia);
    media = "";
    for (let i = 0; i < Math.max(a.length, b.length); i += 1) {
      media += a[i] === "-" || b[i] === "-" ? " " : a[i] === b[i] ? "|" : ".";
    }
  }
  const largo = Math.max(...filas.map((f) => f.secuencia.length));
  const bloque = 60;
  const pre = el("pre", { class: "alineamiento" });
  for (let inicio = 0; inicio < largo; inicio += bloque) {
    const fin = Math.min(largo, inicio + bloque);
    const pintarFila = (fila, otra) => {
      pre.append(el("span", { class: "nombre-fila", texto: fila.nombre.slice(0, ancho).padEnd(ancho + 2) }));
      for (let i = inicio; i < fin; i += 1) {
        const letra = fila.secuencia[i] ?? " ";
        if (letra === "-") {
          pre.append(el("span", { class: "hueco", texto: letra }));
        } else if (otra && otra.secuencia[i] && otra.secuencia[i] !== "-" && otra.secuencia[i] !== letra) {
          pre.append(el("span", { class: "distinto", texto: letra }));
        } else {
          pre.append(document.createTextNode(letra));
        }
      }
      pre.append(document.createTextNode(`  ${fin}\n`));
    };
    if (reales.length === 2) {
      pintarFila(reales[0], reales[1]);
      if (media !== null) {
        pre.append(document.createTextNode(" ".repeat(ancho + 2) + media.slice(inicio, fin) + "\n"));
      }
      pintarFila(reales[1], reales[0]);
    } else {
      for (const fila of reales) {
        pintarFila(fila, null);
      }
    }
    pre.append(document.createTextNode("\n"));
  }
  return pre;
}

function vistaTabla(bloque) {
  const cuerpo = el("tbody", {}, bloque.filas.map((fila) => el("tr", {}, fila.map((celda) => el("td", { texto: celda })))));
  const tabla = el("table", { class: "tabla" }, [el("thead", {}, [el("tr", {}, bloque.columnas.map((c) => el("th", { scope: "col", texto: c })))]), cuerpo]);
  return el("div", { class: "tabla-caja", tabindex: 0 }, [tabla]);
}

function notaTotal(mostrados, total, palabra) {
  return total > mostrados ? `mostrando ${formatear(mostrados)} de ${formatear(total)} ${palabra}` : `${formatear(total)} ${palabra}`;
}

export function pintarBloques(bloques) {
  const partes = [];
  let valores = null;
  for (const bloque of bloques) {
    if (bloque.tipo === "valor") {
      if (!valores) {
        valores = el("div", { class: "valores" });
        partes.push(valores);
      }
      valores.append(el("div", { class: "valor" }, [el("span", { class: "v-etiqueta", texto: bloque.etiqueta }), el("span", { class: "v-dato", texto: bloque.valor })]));
      continue;
    }
    valores = null;
    partes.push(pintarBloque(bloque));
  }
  return partes;
}

export function pintarBloque(bloque) {
  switch (bloque.tipo) {
    case "texto":
      return figura(bloque.etiqueta, `${formatear(bloque.texto.length)} caracteres`, [el("pre", { class: "tipeado", texto: bloque.texto })], [botonCopiar(() => bloque.texto)]);
    case "lista":
      return figura(bloque.etiqueta, notaTotal(bloque.elementos.length, bloque.total, "elementos"), [el("ul", { class: "lista-elementos" }, bloque.elementos.map((e) => el("li", { texto: e })))], [botonCopiar(() => bloque.elementos.join("\n"))]);
    case "tabla": {
      const guardar = el("button", { type: "button", class: "enlace", texto: "Guardar CSV" });
      guardar.addEventListener("click", () => guardarComo(nombreDeArchivo(bloque.etiqueta, "csv"), csv(bloque.columnas, bloque.filas)));
      return figura(bloque.etiqueta, notaTotal(bloque.filas.length, bloque.total, "filas"), [vistaTabla(bloque)], [botonCopiar(() => [bloque.columnas, ...bloque.filas].map((f) => f.join("\t")).join("\n")), guardar]);
    }
    case "linea":
      return figuraGrafica(bloque.etiqueta, null, graficaLinea(bloque));
    case "barras":
      return figuraGrafica(bloque.etiqueta, null, graficaBarras(bloque));
    case "dispersion":
      return figuraGrafica(bloque.etiqueta, `${formatear(bloque.puntos.length)} puntos`, graficaDispersion(bloque));
    case "calor":
      return figuraGrafica(bloque.etiqueta, `${bloque.filas.length} × ${bloque.columnas.length}`, graficaCalor(bloque));
    case "arbol":
    {
      const arbol = graficaArbol(bloque.newick);
      return figura(bloque.etiqueta, null, [arbol, el("pre", { class: "tipeado", texto: bloque.newick })], [botonCopiar(() => bloque.newick), ...botonesDeExportar(bloque.etiqueta, arbol)]);
    }
    case "grafo":
      return figuraGrafica(bloque.etiqueta, `${bloque.nodos.length} nodos, ${bloque.aristas.length} aristas`, graficaGrafo(bloque.nodos, bloque.aristas, { etiqueta: bloque.etiqueta }));
    case "alineamiento":
      return figura(bloque.etiqueta, null, [vistaAlineamiento(bloque.filas)], [botonCopiar(() => bloque.filas.map((f) => `>${f.nombre}\n${f.secuencia}`).join("\n"))]);
    case "archivo": {
      const mostrar = el("button", { type: "button", class: "enlace", texto: "Mostrar en la carpeta" });
      mostrar.addEventListener("click", () => mostrarEnCarpeta(bloque.ruta));
      return figura(bloque.etiqueta, null, [el("div", { class: "archivo" }, [el("code", { texto: bloque.ruta }), mostrar, botonCopiar(() => bloque.ruta)])]);
    }
    default:
      return figura(bloque.etiqueta || "Resultado", null, [el("pre", { class: "tipeado", texto: JSON.stringify(bloque, null, 2) })]);
  }
}
