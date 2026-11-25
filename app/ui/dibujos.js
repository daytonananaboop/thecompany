import { el, sv, color, colorSerie, formatear, colorCalor, graficaGrafo, ubicarGrafo, graficaDispersion, graficaArbol } from "./bloques.js";

function poner(lienzo, ...partes) {
  lienzo.replaceChildren(...partes.filter(Boolean));
}

function valores(pares) {
  return el("div", { class: "valores" }, pares.filter(Boolean).map(([etiqueta, valor]) => el("div", { class: "valor" }, [el("span", { class: "v-etiqueta", texto: etiqueta }), el("span", { class: "v-dato", texto: valor })])));
}

function escala(d0, d1, r0, r1) {
  const ancho = d1 - d0 || 1;
  return (v) => r0 + ((v - d0) / ancho) * (r1 - r0);
}

function sesgo(lienzo, estado) {
  const ancho = 520;
  const alto = 240;
  const area = { izq: 60, der: ancho - 14, arriba: 14, abajo: alto - 34 };
  const [ymin, ymax] = estado.eje_y;
  const x = escala(0, Math.max(1, estado.longitud), area.izq, area.der);
  const y = escala(ymin - 0.5, ymax + 0.5, area.abajo, area.arriba);
  const tinta = color("--tinta");
  const svg = sv("svg", { viewBox: `0 0 ${ancho} ${alto}`, role: "img", "aria-label": "Sesgo G − C acumulado" });
  svg.append(sv("line", { x1: area.izq, x2: area.der, y1: y(0), y2: y(0), stroke: color("--eje"), "stroke-dasharray": "3 3" }));
  svg.append(sv("line", { x1: area.izq, x2: area.der, y1: area.abajo, y2: area.abajo, stroke: color("--eje") }));
  svg.append(sv("line", { x1: area.izq, x2: area.izq, y1: area.arriba, y2: area.abajo, stroke: color("--eje") }));
  for (const valor of [ymin, 0, ymax]) {
    svg.append(sv("text", { x: area.izq - 6, y: y(valor) + 4, "text-anchor": "end", texto: formatear(valor) }));
  }
  svg.append(sv("text", { x: area.der, y: alto - 8, "text-anchor": "end", texto: `${formatear(estado.longitud)} pb` }));
  svg.append(sv("text", { x: area.izq, y: alto - 8, texto: "0" }));
  svg.append(sv("line", { x1: area.izq, x2: area.der, y1: y(estado.minimo), y2: y(estado.minimo), stroke: color("--rojo"), "stroke-width": 1, "stroke-dasharray": "2 4" }));
  if (estado.curva.length) {
    const d = estado.curva.map((p, i) => (i ? "L" : "M") + x(p[0]).toFixed(1) + " " + y(p[1]).toFixed(1)).join(" ");
    svg.append(sv("path", { d, fill: "none", stroke: color("--serie-1"), "stroke-width": 2, "stroke-linejoin": "round" }));
  }
  for (const posicion of estado.posiciones_minimo) {
    svg.append(sv("path", { d: `M${x(posicion)} ${y(estado.minimo) + 4} l-5 9 l10 0 z`, fill: color("--rojo") }));
  }
  svg.append(sv("circle", { cx: x(estado.posicion), cy: y(estado.sesgo), r: 5, fill: color("--rojo"), stroke: color("--hoja"), "stroke-width": 2 }));
  const partes = [
    valores([
      ["Posición", formatear(estado.posicion)],
      ["Sesgo", formatear(estado.sesgo)],
      ["Mínimo hasta aquí", formatear(estado.minimo)],
      ["Veces en el mínimo", formatear(estado.total_minimos)],
      estado.base ? ["Base leída", estado.base] : null,
    ]),
    el("div", { class: "grafica" }, [svg]),
  ];
  if (estado.secuencia) {
    const pre = el("pre", { class: "alineamiento" });
    [...estado.secuencia].forEach((letra, i) => {
      const leida = i < estado.posicion;
      const actual = i === estado.posicion - 1;
      const clase = actual ? "distinto" : leida ? "" : "hueco";
      pre.append(clase ? el("span", { class: clase, texto: letra }) : document.createTextNode(letra));
      if ((i + 1) % 60 === 0) {
        pre.append(document.createTextNode("\n"));
      }
    });
    partes.push(pre);
  }
  poner(lienzo, ...partes);
}

function tablaPerfil(perfil) {
  if (!perfil) {
    return null;
  }
  const bases = ["A", "C", "G", "T"];
  const k = perfil[0].length;
  return el("table", { class: "mini-tabla", "aria-label": "Perfil" }, [
    el("tr", {}, [el("th", { texto: "" }), ...Array.from({ length: k }, (_, j) => el("th", { texto: String(j + 1) }))]),
    ...perfil.map((fila, i) => el("tr", {}, [
      el("th", { texto: bases[i] }),
      ...fila.map((p) => el("td", { texto: p.toFixed(2), style: `background:${colorCalor(p, 0, 1)}` })),
    ])),
  ]);
}

function motivos(lienzo, estado) {
  const pre = el("pre", { class: "alineamiento" });
  const largoNombre = Math.min(16, Math.max(2, ...estado.nombres.map((n) => n.length)));
  estado.secuencias.forEach((secuencia, fila) => {
    const posicion = fila < estado.posiciones.length ? estado.posiciones[fila] : null;
    const actual = estado.fila === fila;
    pre.append(el("span", { class: "nombre-fila", texto: (actual ? "▸ " : "  ") + estado.nombres[fila].slice(0, largoNombre).padEnd(largoNombre + 1) }));
    if (posicion === null) {
      pre.append(el("span", { class: "hueco", texto: secuencia }));
    } else {
      pre.append(document.createTextNode(secuencia.slice(0, posicion)));
      pre.append(el("span", { class: "distinto", texto: secuencia.slice(posicion, posicion + estado.k) }));
      pre.append(document.createTextNode(secuencia.slice(posicion + estado.k)));
    }
    pre.append(document.createTextNode("\n"));
  });
  const partes = [
    valores([
      ["Ronda", estado.ronda === null ? "—" : `${estado.ronda + 1} de ${estado.rondas}`],
      ["k", String(estado.k)],
      estado.puntos !== null ? ["Puntuación", formatear(estado.puntos)] : null,
      ["Mejor puntuación", formatear(estado.mejores.puntos)],
      estado.elegido ? ["Elegido", estado.elegido.kmer] : null,
    ]),
    pre,
  ];
  if (estado.probabilidades) {
    const ancho = 520;
    const alto = 70;
    const maximo = Math.max(...estado.probabilidades, 1e-12);
    const paso = ancho / estado.probabilidades.length;
    const svg = sv("svg", { viewBox: `0 0 ${ancho} ${alto}`, role: "img", "aria-label": "Probabilidad de cada ventana" });
    estado.probabilidades.forEach((p, i) => {
      const h = (p / maximo) * (alto - 16);
      const elegida = estado.elegido && estado.elegido.posicion === i;
      svg.append(sv("rect", { x: i * paso + 0.5, y: alto - 12 - h, width: Math.max(1, paso - 1), height: Math.max(0.5, h), fill: elegida ? color("--rojo") : color("--serie-1") }));
    });
    svg.append(sv("text", { x: 0, y: alto - 1, texto: "probabilidad de cada ventana de la fila actual" }));
    partes.push(el("div", { class: "grafica" }, [svg]));
  }
  const perfil = tablaPerfil(estado.perfil);
  if (perfil) {
    partes.push(el("div", { style: "overflow-x:auto" }, [perfil]));
  }
  poner(lienzo, ...partes);
}

const FLECHAS = { diagonal: "↖", arriba: "↑", izquierda: "←" };

function alineamientoGlobal(lienzo, estado) {
  const camino = new Set(estado.camino.map(([i, j]) => `${i},${j}`));
  const filas = [];
  filas.push(el("tr", {}, [el("th", { texto: "" }), el("th", { texto: "−" }), ...[...estado.w].map((letra) => el("th", { texto: letra }))]));
  estado.matriz.forEach((fila, i) => {
    const celdas = [el("th", { texto: i === 0 ? "−" : estado.v[i - 1] })];
    fila.forEach((valor, j) => {
      const flecha = estado.flechas[i][j];
      const clases = [];
      if (valor === null) {
        clases.push("vacia");
      }
      if (camino.has(`${i},${j}`)) {
        clases.push("camino");
      }
      if (estado.celda && estado.celda[0] === i && estado.celda[1] === j) {
        clases.push("actual");
      }
      celdas.push(el("td", { class: clases.join(" "), texto: valor === null ? "·" : `${flecha ? FLECHAS[flecha] : ""}${valor}` }));
    });
    filas.push(el("tr", {}, celdas));
  });
  const partes = [
    valores([
      ["Puntuación", estado.puntuacion],
      ["Hueco", `−${estado.indel}`],
      estado.par ? ["Par", `${estado.par.letra_v} / ${estado.par.letra_w}: ${estado.par.puntaje >= 0 ? "+" : ""}${estado.par.puntaje}`] : null,
      estado.puntaje !== null ? ["Puntaje final", formatear(estado.puntaje)] : null,
    ]),
    el("div", { style: "overflow-x:auto" }, [el("table", { class: "mini-tabla", "aria-label": "Tabla de programación dinámica" }, filas)]),
  ];
  if (estado.candidatos) {
    partes.push(el("div", { class: "valores" }, ["diagonal", "arriba", "izquierda"].map((origen) => el("div", { class: "valor" }, [
      el("span", { class: "v-etiqueta", texto: `${FLECHAS[origen]} ${origen}` }),
      el("span", { class: "v-dato", style: origen === estado.ganador ? `color:${color("--rojo")}` : "", texto: formatear(estado.candidatos[origen]) }),
    ]))));
  }
  if (estado.fila_v || estado.fila_w) {
    partes.push(el("pre", { class: "alineamiento", texto: `${estado.nombre_v.slice(0, 12).padEnd(13)}${estado.fila_v}\n${estado.nombre_w.slice(0, 12).padEnd(13)}${estado.fila_w}` }));
  }
  poner(lienzo, ...partes);
}

function cromosomasCirculares(genoma, marcados) {
  const ancho = 520;
  const columnas = Math.min(3, Math.max(1, genoma.length));
  const radio = Math.min(90, ancho / columnas / 2 - 34);
  const filas = Math.ceil(genoma.length / columnas);
  const alto = filas * (radio * 2 + 60) + 10;
  const svg = sv("svg", { viewBox: `0 0 ${ancho} ${alto}`, role: "img", "aria-label": "Cromosomas circulares" });
  const tinta = color("--tinta");
  genoma.forEach((cromosoma, c) => {
    const cx = (c % columnas) * (ancho / columnas) + ancho / columnas / 2;
    const cy = Math.floor(c / columnas) * (radio * 2 + 60) + radio + 30;
    svg.append(sv("circle", { cx, cy, r: radio, fill: "none", stroke: color("--tinta-tenue"), "stroke-width": 1 }));
    const n = cromosoma.length;
    cromosoma.forEach((bloque, b) => {
      const a0 = (2 * Math.PI * b) / n - Math.PI / 2 + 0.12;
      const a1 = (2 * Math.PI * (b + 1)) / n - Math.PI / 2 - 0.12;
      const [desde, hasta] = bloque > 0 ? [a0, a1] : [a1, a0];
      const x0 = cx + radio * Math.cos(desde);
      const y0 = cy + radio * Math.sin(desde);
      const x1 = cx + radio * Math.cos(hasta);
      const y1 = cy + radio * Math.sin(hasta);
      const barrido = bloque > 0 ? 1 : 0;
      const marcado = marcados.has(Math.abs(bloque));
      const trazo = marcado ? color("--rojo") : colorSerie(bloque > 0 ? 0 : 1);
      const flecha = `fl${c}-${b}`;
      svg.append(sv("defs", {}, [sv("marker", { id: flecha, viewBox: "0 0 10 10", refX: 6, refY: 5, markerWidth: 5, markerHeight: 5, orient: "auto" }, [sv("path", { d: "M0 0 L10 5 L0 10 Z", fill: trazo })])]));
      svg.append(sv("path", { d: `M${x0} ${y0} A${radio} ${radio} 0 ${Math.abs(a1 - a0) > Math.PI ? 1 : 0} ${barrido} ${x1} ${y1}`, fill: "none", stroke: trazo, "stroke-width": marcado ? 6 : 4, "marker-end": `url(#${flecha})` }));
      const medio = (a0 + a1) / 2;
      svg.append(sv("text", { x: cx + (radio + 16) * Math.cos(medio), y: cy + (radio + 16) * Math.sin(medio) + 4, "text-anchor": "middle", style: `fill:${tinta}; font-size:12px`, texto: (bloque > 0 ? "+" : "−") + Math.abs(bloque) }));
    });
  });
  return el("div", { class: "grafica" }, [svg]);
}

function dosRupturas(lienzo, estado) {
  const marcados = new Set((estado.rotura || []).map((nodo) => Math.ceil(nodo / 2)));
  const descripcion = estado.rotura
    ? `Se cortan (${estado.rotura[0]}, ${estado.rotura[1]}) y (${estado.rotura[2]}, ${estado.rotura[3]}); se unen (${estado.rotura[0]}, ${estado.rotura[2]}) y (${estado.rotura[1]}, ${estado.rotura[3]}).`
    : "";
  poner(lienzo, 
    valores([
      ["2-rupturas hechas", formatear(estado.paso)],
      ["Ciclos", formatear(estado.ciclos)],
      ["Distancia que falta", formatear(estado.distancia)],
      estado.efecto ? ["Efecto", estado.efecto] : null,
    ]),
    el("pre", { class: "alineamiento", texto: `P ahora   ${estado.texto}\nobjetivo  ${estado.texto_objetivo}` }),
    cromosomasCirculares(estado.genoma, marcados),
    descripcion ? el("p", { class: "resumen", texto: descripcion }) : null,
  );
}

const DISPOSICIONES = new WeakMap();

function disposicionFija(lista, ancho) {
  if (!DISPOSICIONES.has(lista)) {
    const completo = lista.reduce((mejor, paso) => (paso.estado.nodos.length > mejor.nodos.length ? paso.estado : mejor), lista[0].estado);
    const final = [...lista].reverse().find((paso) => paso.estado.camino.length && !paso.estado.camino_al_reves);
    const orden = [];
    const visto = new Set();
    for (const nodo of final ? final.estado.camino : []) {
      if (!visto.has(nodo)) {
        visto.add(nodo);
        orden.push(nodo);
      }
    }
    completo.nodos.forEach((_, i) => {
      if (!visto.has(i)) {
        orden.push(i);
      }
    });
    const anchoNodo = Math.max(18, Math.max(...completo.nodos.map((n) => n.length)) * 7 + 10);
    const porFila = Math.max(2, Math.floor((ancho - 20) / (anchoNodo + 36)));
    const paso = (ancho - 20) / porFila;
    const posiciones = [];
    orden.forEach((nodo, i) => {
      const fila = Math.floor(i / porFila);
      const columna = fila % 2 === 0 ? i % porFila : porFila - 1 - (i % porFila);
      posiciones[nodo] = { x: 10 + paso * columna + paso / 2, y: 40 + fila * 74 };
    });
    DISPOSICIONES.set(lista, { posiciones, alto: 40 + Math.ceil(orden.length / porFila) * 74 });
  }
  return DISPOSICIONES.get(lista);
}

function deBruijn(lienzo, estado, indice, lista) {
  const { posiciones, alto } = disposicionFija(lista, 520);
  const aristas = estado.aristas.map((a) => ({ desde: a.desde, hasta: a.hasta, etiqueta: a.kmero }));
  const usadas = new Set(estado.aristas.flatMap((a, i) => (a.usada ? [i] : [])));
  const grafo = graficaGrafo(estado.nodos, aristas, {
    posiciones,
    alto,
    ancho: 520,
    usadas: estado.fase === "grafo" ? null : usadas,
    resaltadas: new Set(estado.arista === null ? [] : [estado.arista]),
    activos: new Set(estado.nodo === null ? [] : [estado.nodo]),
    etiqueta: "Grafo de De Bruijn",
  });
  const nombre = (i) => estado.nodos[i] ?? "?";
  const lineas = [];
  if (estado.balance) {
    const inicios = estado.balance.flatMap((b, i) => (b > 0 ? [`${nombre(i)} (+${b})`] : []));
    const finales = estado.balance.flatMap((b, i) => (b < 0 ? [`${nombre(i)} (${b})`] : []));
    lineas.push(`salidas − entradas: ${inicios.length ? "empieza en " + inicios.join(", ") : "todo balanceado (ciclo)"}${finales.length ? "; termina en " + finales.join(", ") : ""}`);
  }
  if (estado.pila.length) {
    lineas.push(`pila: ${estado.pila.map(nombre).join(" → ")}`);
  }
  if (estado.camino.length) {
    lineas.push(`camino${estado.camino_al_reves ? " (todavía al revés)" : ""}: ${estado.camino.map(nombre).join(" → ")}`);
  }
  if (estado.texto) {
    lineas.push(`texto: ${estado.texto}`);
  }
  poner(lienzo, 
    valores([
      ["Fase", { grafo: "armar el grafo", inicio: "elegir dónde empezar", recorrido: "recorrer (Hierholzer)", deletreo: "deletrear" }[estado.fase] || estado.fase],
      ["k", String(estado.k)],
      ["Nodos", formatear(estado.nodos.length)],
      ["Aristas usadas", `${usadas.size} de ${estado.aristas.length}`],
    ]),
    grafo,
    lineas.length ? el("pre", { class: "alineamiento", texto: lineas.join("\n") }) : null,
  );
}

function tiraDeEspectro(estado) {
  const ancho = 520;
  const alto = 64;
  const x = escala(0, Math.max(1, estado.masa_madre), 10, ancho - 10);
  const svg = sv("svg", { viewBox: `0 0 ${ancho} ${alto}`, role: "img", "aria-label": "Espectro experimental" });
  svg.append(sv("line", { x1: 10, x2: ancho - 10, y1: 34, y2: 34, stroke: color("--eje") }));
  const bloques = new Set(estado.bloques.map((b) => b.masa));
  for (const masa of estado.espectro) {
    const esBloque = bloques.has(masa);
    svg.append(sv("line", { x1: x(masa), x2: x(masa), y1: esBloque ? 14 : 22, y2: 34, stroke: esBloque ? color("--rojo") : color("--serie-1"), "stroke-width": esBloque ? 2.5 : 1.5 }));
  }
  svg.append(sv("text", { x: 10, y: 54, texto: "0" }));
  svg.append(sv("text", { x: ancho - 10, y: 54, "text-anchor": "end", texto: `${formatear(estado.masa_madre)} Da` }));
  svg.append(sv("text", { x: ancho / 2, y: 54, "text-anchor": "middle", texto: "rojo: un solo aminoácido" }));
  return el("div", { class: "grafica" }, [svg]);
}

function ciclopeptido(lienzo, estado) {
  const fila = (p) => el("tr", {}, [el("td", { texto: p.letras || "(vacío)" }), el("td", { texto: p.masas.join("-") || "—" }), el("td", { texto: formatear(p.masa) })]);
  const tabla = (titulo, lista, total) => el("figure", { class: "bloque" }, [
    el("figcaption", {}, [el("span", { texto: titulo }), el("span", { class: "nota-bloque", texto: total > lista.length ? `mostrando ${lista.length} de ${formatear(total)}` : `${formatear(total)}` })]),
    lista.length
      ? el("div", { class: "tabla-caja" }, [el("table", { class: "tabla" }, [el("thead", {}, [el("tr", {}, ["Péptido", "Masas", "Total"].map((t) => el("th", { texto: t })))]), el("tbody", {}, lista.map(fila))])])
      : el("p", { class: "resumen", texto: "Ninguno todavía." }),
  ]);
  poner(lienzo, 
    valores([
      ["Ronda", formatear(estado.ronda)],
      ["Expandidos", formatear(estado.expandidos)],
      ["Inconsistentes", formatear(estado.inconsistentes)],
      ["Descartados", formatear(estado.descartados)],
      ["Siguen", formatear(estado.total_candidatos)],
      ["Encontrados", formatear(estado.total_encontrados)],
    ]),
    tiraDeEspectro(estado),
    el("p", { class: "resumen", texto: `Bloques: ${estado.bloques.map((b) => `${b.letras} ${b.masa}`).join(", ")}` }),
    tabla(estado.fase === "expandir" ? "Candidatos recién expandidos" : "Candidatos que siguen", estado.candidatos, estado.total_candidatos),
    estado.total_encontrados ? tabla("Péptidos encontrados", estado.encontrados, estado.total_encontrados) : null,
  );
}

function matrizConNombres(nombres, matriz, resaltar, minimo) {
  return el("div", { style: "overflow-x:auto" }, [el("table", { class: "mini-tabla" }, [
    el("tr", {}, [el("th", { texto: "" }), ...nombres.map((n) => el("th", { texto: n }))]),
    ...matriz.map((fila, i) => el("tr", {}, [
      el("th", { texto: nombres[i] }),
      ...fila.map((valor, j) => {
        const elegido = resaltar && ((i === resaltar[0] && j === resaltar[1]) || (i === resaltar[1] && j === resaltar[0]));
        const enFila = resaltar && !elegido && (resaltar.includes(i) || resaltar.includes(j));
        return el("td", { class: elegido ? "actual" : enFila ? "camino" : "", texto: i === j ? "—" : formatear(valor) });
      }),
    ])),
  ])]);
}

function unionDeVecinos(lienzo, estado) {
  const enEleccion = estado.fase === "eleccion" && estado.matriz_de_union;
  const par = estado.par;
  const titulo = enEleccion ? "Matriz D* (se une el par con el menor valor)" : "Matriz de distancias D";
  const aristas = estado.aristas.map((a) => `${a.desde} — ${a.hasta}: ${formatear(a.longitud)}`).join("\n");
  poner(
    lienzo,
    valores([
      ["Ronda", formatear(estado.ronda)],
      estado.nuevo ? ["Nodo nuevo", estado.nuevo] : null,
      par ? ["Par", `${estado.nombres[par[0]]} + ${estado.nombres[par[1]]}`] : null,
      estado.minimo !== null ? ["D* mínimo", formatear(estado.minimo)] : null,
      estado.delta !== null ? ["Δ", formatear(estado.delta)] : null,
      estado.ramas ? ["Ramas", `${formatear(estado.ramas[0])} y ${formatear(estado.ramas[1])}`] : null,
    ]),
    el("p", { class: "resumen", texto: titulo }),
    matrizConNombres(estado.nombres, enEleccion ? estado.matriz_de_union : estado.distancias, par, estado.minimo),
    estado.totales ? el("p", { class: "resumen", texto: `Totales por fila: ${estado.nombres.map((n, i) => `${n} ${formatear(estado.totales[i])}`).join(", ")}` }) : null,
    estado.newick ? graficaArbol(estado.newick) : aristas ? el("pre", { class: "alineamiento", texto: `Ramas del árbol hasta ahora:\n${aristas}` }) : null,
  );
}

function lloyd(lienzo, estado) {
  const puntos = estado.puntos.map(([x, y], i) => ({ x, y, grupo: estado.grupos ? estado.grupos[i] : -1, nombre: estado.nombres[i] }));
  const centros = estado.centros.map(([x, y], i) => ({ x, y, grupo: i }));
  const movimientos = estado.centros_anteriores ? estado.centros_anteriores.map(([x, y], i) => [{ x, y }, centros[i]]) : [];
  const fase = { inicio: "elegir centros", asignar: "asignar cada punto a su centro", mover: "mover cada centro al promedio", fin: "fin" }[estado.fase] || estado.fase;
  poner(
    lienzo,
    valores([
      ["Iteración", formatear(estado.iteracion)],
      ["Paso", fase],
      ["Distorsión", formatear(estado.distorsion)],
      estado.cambiaron !== null ? ["Cambiaron de grupo", formatear(estado.cambiaron)] : null,
      estado.cantidades ? ["Tamaños", estado.cantidades.join(" / ")] : null,
    ]),
    graficaDispersion({ etiqueta: "k-medias", eje_x: estado.ejes[0], eje_y: estado.ejes[1], puntos }, { centros, movimientos, alto: 360, ancho: 560 }),
  );
}

const PRIORIDAD_DE_MARCA = { actual: 3, hallado: 2, comparado: 1 };

function vistaDeSecuencia(secuencia) {
  const pre = el("pre", { class: "alineamiento pizarra-secuencia" });
  const texto = secuencia.texto;
  const inicio = secuencia.desplazamiento;
  const porLinea = 60;
  for (let fila = 0; fila < texto.length; fila += porLinea) {
    pre.append(el("span", { class: "nombre-fila", texto: String(inicio + fila).padStart(9) + "  " }));
    for (let i = fila; i < Math.min(texto.length, fila + porLinea); i += 1) {
      const posicion = inicio + i;
      let tipo = null;
      for (const marca of secuencia.marcas) {
        if (posicion >= marca.desde && posicion < marca.hasta && (!tipo || PRIORIDAD_DE_MARCA[marca.tipo] > PRIORIDAD_DE_MARCA[tipo])) {
          tipo = marca.tipo;
        }
      }
      pre.append(tipo ? el("span", { class: `marca-${tipo}`, texto: texto[i] }) : document.createTextNode(texto[i]));
    }
    pre.append(document.createTextNode("\n"));
  }
  const fin = inicio + texto.length;
  return el("figure", { class: "bloque" }, [
    el("figcaption", {}, [el("span", { texto: "Secuencia" }), el("span", { class: "nota-bloque", texto: secuencia.largo > texto.length ? `posiciones ${formatear(inicio)} a ${formatear(fin - 1)} de ${formatear(secuencia.largo)}, contando desde 0` : "contando desde 0" })]),
    pre,
  ]);
}

function pizarra(lienzo, estado, indice, lista) {
  const anterior = indice > 0 && lista[indice - 1].estado && lista[indice - 1].estado.tipo === "pizarra" ? lista[indice - 1].estado : null;
  const previos = new Map((anterior ? anterior.variables : []).map(([nombre, valor]) => [nombre, valor]));
  const variables = (estado.variables || []).map(([nombre, valor]) => {
    const cambio = anterior !== null && previos.get(nombre) !== valor;
    return el("div", { class: "variable" + (cambio ? " cambio" : "") }, [
      el("span", { class: "nombre-variable", texto: nombre }),
      el("span", { class: "valor-variable", texto: valor }),
    ]);
  });
  const partes = [];
  if (variables.length) {
    partes.push(el("figure", { class: "bloque" }, [
      el("figcaption", {}, [el("span", { texto: "Variables" }), el("span", { class: "nota-bloque", texto: "las que acaban de cambiar están marcadas" })]),
      el("div", { class: "pizarra-variables" }, variables),
    ]));
  }
  if (estado.secuencia) {
    partes.push(vistaDeSecuencia(estado.secuencia));
  }
  if (estado.tabla) {
    const resaltadas = new Set((estado.tabla.resaltar || []).map(([f, c]) => `${f},${c}`));
    partes.push(el("figure", { class: "bloque" }, [
      el("figcaption", {}, [el("span", { texto: estado.tabla.titulo || "Tabla" })]),
      el("div", { style: "overflow-x:auto" }, [el("table", { class: "mini-tabla" }, [
        el("tr", {}, estado.tabla.columnas.map((c) => el("th", { texto: c }))),
        ...estado.tabla.filas.map((fila, i) => el("tr", {}, fila.map((celda, j) => el("td", { class: resaltadas.has(`${i},${j}`) ? "actual" : "", texto: celda })))),
      ])]),
    ]));
  }
  if (estado.lista) {
    const resaltados = new Set(estado.lista.resaltar || []);
    partes.push(el("figure", { class: "bloque" }, [
      el("figcaption", {}, [el("span", { texto: estado.lista.titulo || "Lista" })]),
      el("ul", { class: "lista-elementos" }, estado.lista.elementos.map((e, i) => el("li", { class: resaltados.has(i) ? "marca-actual" : "", texto: e }))),
    ]));
  }
  poner(lienzo, ...partes);
}

export const PIZARRA = pizarra;

export const DIBUJOS = {
  "replicacion.sesgo_minimo": sesgo,
  "motivos.codiciosa_laplace": motivos,
  "alineamiento.global": alineamientoGlobal,
  "reordenamientos.transformacion_mas_corta": dosRupturas,
  "ensamblaje.reconstruccion": deBruijn,
  "antibioticos.ciclopeptido": ciclopeptido,
  "filogenia.union_de_vecinos": unionDeVecinos,
  "agrupamiento.lloyd": lloyd,
};

export const AYUDANTES = { valores, graficaGrafo, ubicarGrafo, graficaDispersion, graficaArbol };
