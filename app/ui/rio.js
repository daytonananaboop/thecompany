const SVG = "http://www.w3.org/2000/svg";

const NOMBRES_CORTOS = {
  replicacion: "Replicación",
  motivos: "Motivos",
  ensamblaje: "Ensamblaje",
  antibioticos: "Péptidos",
  alineamiento: "Alineamiento",
  reordenamientos: "Reordenamientos",
  filogenia: "Filogenia",
  agrupamiento: "Agrupamiento",
  puerto: "El Puerto",
};

const ROMANOS = ["", "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX"];

function nodo(nombre, atributos = {}, hijos = []) {
  const elemento = document.createElementNS(SVG, nombre);
  for (const [clave, valor] of Object.entries(atributos)) {
    elemento.setAttribute(clave, valor);
  }
  for (const hijo of hijos) {
    elemento.append(hijo);
  }
  return elemento;
}

function texto(contenido, atributos) {
  const elemento = nodo("text", atributos);
  elemento.textContent = contenido;
  return elemento;
}

function plural(cantidad, uno, varios) {
  return `${cantidad} ${cantidad === 1 ? uno : varios}`;
}

function alturaDelRio(x) {
  return 176 + 58 * Math.sin(((x - 170) / 1010) * Math.PI * 2 * 1.75);
}

function trazoDelRio() {
  let d = "";
  for (let x = 150; x <= 1180; x += 10) {
    d += (x === 150 ? "M" : " L") + x + " " + alturaDelRio(x).toFixed(1);
  }
  return d;
}

function azar(semilla) {
  let estado = semilla;
  return () => {
    estado = (estado * 1103515245 + 12345) % 2147483648;
    return estado / 2147483648;
  };
}

function selva(ocupado) {
  const grupo = nodo("g", { class: "selva", "aria-hidden": "true" });
  const siguiente = azar(7);
  for (let i = 0; i < 230; i += 1) {
    const x = 170 + siguiente() * 1010;
    const y = 18 + siguiente() * 310;
    if (Math.abs(y - alturaDelRio(x)) < 46) {
      continue;
    }
    if (ocupado.some(([cx, cy, ancho, alto]) => Math.abs(x - cx) < ancho && Math.abs(y - cy) < alto)) {
      continue;
    }
    grupo.append(nodo("path", { d: `M${x - 4} ${y + 3} L${x} ${y - 4} L${x + 4} ${y + 3} M${x} ${y - 4} L${x} ${y + 6}` }));
  }
  return grupo;
}

function choza(x, y, arriba) {
  const base = arriba ? y - 30 : y + 36;
  return nodo("g", {}, [
    nodo("line", { x1: x - 4, y1: arriba ? y - 17 : y + 17, x2: x - 4, y2: arriba ? base + 12 : base + 2, class: "s-tenue", "stroke-dasharray": "2 3" }),
    nodo("path", { class: "choza", d: `M${x - 9} ${base + 12} L${x - 9} ${base + 3} L${x} ${base - 5} L${x + 9} ${base + 3} L${x + 9} ${base + 12} Z` }),
    nodo("line", { x1: x + 6, y1: base - 1, x2: x + 6, y2: base - 16, class: "s-tinta", "stroke-width": 1.4 }),
    nodo("path", { class: "bandera", d: `M${x + 6.5} ${base - 16} L${x + 17} ${base - 12.5} L${x + 6.5} ${base - 9} Z` }),
  ]);
}

function muelle(x, y) {
  return nodo("g", {}, [
    nodo("rect", { class: "choza", x: x - 14, y: y - 8, width: 28, height: 14 }),
    nodo("path", { class: "choza", d: `M${x - 17} ${y - 8} L${x} ${y - 20} L${x + 17} ${y - 8} Z` }),
    nodo("line", { x1: x + 10, y1: y - 20, x2: x + 10, y2: y - 38, class: "s-tinta", "stroke-width": 1.4 }),
    nodo("path", { class: "bandera", d: `M${x + 10.5} ${y - 38} L${x + 23} ${y - 34} L${x + 10.5} ${y - 30} Z` }),
    nodo("line", { x1: x - 22, y1: y + 10, x2: x + 22, y2: y + 10, class: "s-tinta", "stroke-width": 3 }),
    nodo("line", { x1: x - 18, y1: y + 10, x2: x - 18, y2: y + 20, class: "s-tinta", "stroke-width": 2 }),
    nodo("line", { x1: x + 18, y1: y + 10, x2: x + 18, y2: y + 20, class: "s-tinta", "stroke-width": 2 }),
  ]);
}

export function dibujarCarta(contenedor, capitulos, conteos, elegido, alElegir) {
  const lienzo = nodo("svg", { viewBox: "0 0 1200 340", role: "group", "aria-label": "Carta de los capítulos" });
  const defs = nodo("defs", {}, [
    nodo("radialGradient", { id: "niebla" }, [nodo("stop", { offset: "0", class: "niebla-0" }), nodo("stop", { offset: "1", class: "niebla-1" })]),
    nodo("pattern", { id: "olas", width: 24, height: 12, patternUnits: "userSpaceOnUse" }, [
      nodo("path", { d: "M0 6 Q6 0 12 6 T24 6", fill: "none", class: "s-tenue", "stroke-width": 1.2 }),
    ]),
  ]);
  lienzo.append(defs);
  lienzo.append(nodo("rect", { x: 0, y: 0, width: 120, height: 340, fill: "url(#olas)" }));
  lienzo.append(nodo("path", { d: "M120 0 C 128 60, 112 120, 124 170 S 118 280, 122 340", fill: "none", class: "s-tinta", "stroke-width": 1.5 }));
  lienzo.append(texto("El mar", { class: "rotulo-mar", x: 58, y: 130, transform: "rotate(-90 58 130)", "text-anchor": "middle" }));

  const estaciones = capitulos.filter((c) => c.id !== "puerto");
  const posiciones = estaciones.map((capitulo, indice) => {
    const x = 285 + indice * 122;
    return { capitulo, x, y: alturaDelRio(x), arriba: indice % 2 === 0 };
  });
  const ocupado = posiciones.map((p) => [p.x, p.arriba ? p.y - 66 : p.y + 70, 82, 40]);
  ocupado.push([250, 40, 120, 30], [175, 230, 70, 60]);
  lienzo.append(selva(ocupado));

  lienzo.append(nodo("path", { class: "rio-cauce", d: trazoDelRio() }));
  lienzo.append(nodo("path", { class: "rio-agua", d: trazoDelRio() }));
  lienzo.append(nodo("path", { class: "rio-brillo", d: trazoDelRio() }));
  lienzo.append(nodo("ellipse", { class: "niebla-banco", cx: 1140, cy: 150, rx: 140, ry: 120, fill: "url(#niebla)" }));

  const total = Object.values(conteos).reduce((a, b) => a + b, 0);
  lienzo.append(nodo("rect", { x: 150, y: 14, width: 230, height: 46, class: "f-hoja s-tinta", "stroke-width": 1 }));
  lienzo.append(nodo("rect", { x: 154, y: 18, width: 222, height: 38, fill: "none", class: "s-tinta", "stroke-width": 0.6 }));
  lienzo.append(texto("Carta de los capítulos", { class: "cartela", x: 265, y: 37, "text-anchor": "middle" }));
  lienzo.append(texto(`${total} algoritmos en ${capitulos.length} estaciones`, { class: "cartela-llano", x: 265, y: 51, "text-anchor": "middle" }));

  const estacionesSvg = [];
  const crearEstacion = (capitulo, dibujo, etiquetas) => {
    const grupo = nodo("g", {
      class: "estacion" + (capitulo.id === elegido ? " elegida" : ""),
      tabindex: 0,
      role: "button",
      "aria-pressed": capitulo.id === elegido ? "true" : "false",
      "aria-label": `${capitulo.titulo}: ${plural(conteos[capitulo.id] || 0, "algoritmo", "algoritmos")}`,
    });
    grupo.append(dibujo, ...etiquetas);
    grupo.addEventListener("click", () => alElegir(capitulo.id));
    grupo.addEventListener("keydown", (evento) => {
      if (evento.key === "Enter" || evento.key === " ") {
        evento.preventDefault();
        alElegir(capitulo.id);
      }
    });
    estacionesSvg.push(grupo);
    lienzo.append(grupo);
  };

  for (const { capitulo, x, y, arriba } of posiciones) {
    const yNombre = arriba ? y - 70 : y + 70;
    const etiquetas = [
      nodo("rect", { class: "zona", x: x - 70, y: arriba ? y - 108 : y + 12, width: 140, height: 96 }),
      texto(ROMANOS[capitulo.numero] || String(capitulo.numero), { class: "numero-estacion", x, y: arriba ? yNombre - 17 : yNombre + 30, "text-anchor": "middle" }),
      texto(NOMBRES_CORTOS[capitulo.id] || capitulo.titulo, { class: "nombre-estacion", x, y: yNombre, "text-anchor": "middle" }),
      texto(plural(conteos[capitulo.id] || 0, "algoritmo", "algoritmos"), { class: "cuenta-estacion", x, y: yNombre + 14, "text-anchor": "middle" }),
    ];
    crearEstacion(capitulo, choza(x, y, arriba), etiquetas);
  }

  const puerto = capitulos.find((c) => c.id === "puerto");
  if (puerto) {
    const x = 178;
    const y = alturaDelRio(150) + 52;
    crearEstacion(puerto, muelle(x, y), [
      nodo("rect", { class: "zona", x: x - 60, y: y - 44, width: 130, height: 110 }),
      texto("El Puerto", { class: "nombre-estacion", x, y: y + 40, "text-anchor": "middle" }),
      texto(plural(conteos.puerto || 0, "herramienta", "herramientas"), { class: "cuenta-estacion", x, y: y + 54, "text-anchor": "middle" }),
    ]);
  }

  contenedor.replaceChildren(lienzo);
  return estacionesSvg;
}
