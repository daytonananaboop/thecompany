import { invocar, elegirArchivo, recordar, recuperar } from "./puente.js";
import { el, pintarBloques } from "./bloques.js";
import { dibujarCarta } from "./rio.js";
import { pintarCodigo } from "./codigo.js";
import { abrirPasoAPaso } from "./pasos.js";

const ROMANOS = ["", "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX"];

const estado = {
  capitulos: [],
  fichas: [],
  capitulo: null,
  algoritmo: null,
  busqueda: "",
  valores: new Map(),
  salida: null,
  vanguardia: false,
  nucleos: null,
};

function visibles() {
  return estado.fichas.filter((ficha) => estado.vanguardia || !ficha.vanguardia);
}

function aplicarVanguardia() {
  $("vanguardia").setAttribute("aria-pressed", estado.vanguardia ? "true" : "false");
  $("vanguardia-llano").textContent = estado.vanguardia ? "Algoritmos recientes a la vista" : "Ver algoritmos recientes";
}

const $ = (id) => document.getElementById(id);

function normalizar(texto) {
  return texto.normalize("NFD").replace(/[̀-ͯ]/g, "").toLowerCase();
}

function temaActual() {
  const elegido = document.documentElement.dataset.tema;
  if (elegido) {
    return elegido;
  }
  return window.matchMedia("(prefers-color-scheme: dark)").matches ? "noche" : "dia";
}

function aplicarTema(tema) {
  if (tema === "dia" || tema === "noche") {
    document.documentElement.dataset.tema = tema;
  }
  const noche = temaActual() === "noche";
  $("tema").setAttribute("aria-pressed", noche ? "true" : "false");
  $("tema-rio").textContent = noche ? "El río de día" : "El río de noche";
  $("tema-llano").textContent = noche ? "Modo claro" : "Modo oscuro";
}

function conteos() {
  const cuenta = {};
  for (const ficha of visibles()) {
    cuenta[ficha.capitulo] = (cuenta[ficha.capitulo] || 0) + 1;
  }
  return cuenta;
}

function capituloDe(id) {
  return estado.capitulos.find((c) => c.id === id);
}

function fichaDe(id) {
  return estado.fichas.find((f) => f.id === id);
}

function rotuloCapitulo(capitulo) {
  if (!capitulo) {
    return "";
  }
  if (capitulo.id === "puerto") {
    return "El Puerto: herramientas de laboratorio";
  }
  return `Capítulo ${ROMANOS[capitulo.numero] || capitulo.numero}: ${capitulo.titulo}`;
}

function elegirCapitulo(id) {
  estado.capitulo = id;
  estado.busqueda = "";
  $("buscar").value = "";
  const primera = visibles().find((f) => f.capitulo === id);
  if (primera && (!estado.algoritmo || fichaDe(estado.algoritmo).capitulo !== id)) {
    elegirAlgoritmo(primera.id, false);
  }
  pintarCarta();
  pintarIndice();
}

function elegirAlgoritmo(id, enfocar = true) {
  estado.algoritmo = id;
  estado.salida = null;
  const ficha = fichaDe(id);
  if (ficha && ficha.capitulo !== estado.capitulo) {
    estado.capitulo = ficha.capitulo;
    pintarCarta();
  }
  recordar("algoritmo", id);
  pintarIndice();
  pintarPliego();
  if (enfocar) {
    $("pliego").querySelector("h2")?.focus({ preventScroll: true });
  }
}

function pintarCarta() {
  dibujarCarta($("carta"), estado.capitulos, conteos(), estado.capitulo, elegirCapitulo);
}

function botonAlgoritmo(ficha, conCapitulo) {
  const boton = el("button", { type: "button", "aria-current": ficha.id === estado.algoritmo ? "true" : "false" }, [
    el("span", { class: "nombre-algoritmo" }, [document.createTextNode(ficha.nombre), ficha.vanguardia ? el("span", { class: "insignia", texto: "vanguardia" }) : null]),
    conCapitulo ? el("span", { class: "capitulo-algoritmo", texto: rotuloCapitulo(capituloDe(ficha.capitulo)) }) : null,

  ]);
  boton.addEventListener("click", () => elegirAlgoritmo(ficha.id));
  return el("li", {}, [boton]);
}

function pintarIndice() {
  const indice = $("indice");
  if (estado.busqueda) {
    const consulta = normalizar(estado.busqueda);
    const encontradas = visibles().filter((f) => normalizar(`${f.nombre} ${f.resumen} ${f.id} ${f.vanguardia || ""}`).includes(consulta));
    indice.replaceChildren(
      el("h2", { texto: "Búsqueda" }),
      el("p", { class: "resumen-capitulo", texto: encontradas.length ? `${encontradas.length} algoritmos coinciden con «${estado.busqueda}».` : `Nada coincide con «${estado.busqueda}». Prueba con otra palabra, como «alineamiento» o «GC».` }),
      el("ol", {}, encontradas.map((f) => botonAlgoritmo(f, true))),
    );
    return;
  }
  const capitulo = capituloDe(estado.capitulo);
  const lista = visibles().filter((f) => f.capitulo === estado.capitulo);
  indice.replaceChildren(
    el("h2", { texto: rotuloCapitulo(capitulo) }),
    capitulo ? el("p", { class: "resumen-capitulo", texto: capitulo.resumen }) : null,
    el("ol", {}, lista.map((f) => botonAlgoritmo(f, false))),
  );
}

function valoresDe(ficha) {
  if (!estado.valores.has(ficha.id)) {
    const iniciales = {};
    for (const campo of ficha.campos) {
      iniciales[campo.id] = campo.ejemplo;
    }
    estado.valores.set(ficha.id, iniciales);
  }
  return estado.valores.get(ficha.id);
}

function atajo(texto, accion) {
  const boton = el("button", { type: "button", class: "enlace", texto });
  boton.addEventListener("click", accion);
  return boton;
}

function campoGen(prefijo, etiqueta, poner) {
  const entrada = el("input", { type: "text", list: "genes-ecoli", placeholder: "dnaA", "aria-label": etiqueta, spellcheck: "false" });
  entrada.addEventListener("change", () => {
    const nombre = entrada.value.trim();
    if (nombre) {
      poner(`${prefijo}${nombre}`);
    }
  });
  return el("label", {}, [document.createTextNode(etiqueta + " "), entrada]);
}

function controlDeCampo(ficha, campo, valores) {
  const id = `campo-${campo.id}`;
  const clase = campo.tipo.clase;
  const corto = clase === "entero" || clase === "decimal" || clase === "opcion";
  const enLinea = corto || clase === "numeros" || clase === "archivo" || clase === "base";
  let control;
  if (clase === "opcion") {
    control = el("select", { id }, campo.tipo.opciones.map((o) => el("option", { value: o, selected: o === valores[campo.id], texto: o })));
  } else if (enLinea) {
    control = el("input", { id, type: "text", spellcheck: "false", inputmode: corto && clase !== "opcion" ? "decimal" : null });
    control.value = valores[campo.id] ?? "";
  } else {
    const lineas = String(valores[campo.id] ?? "").split("\n").length;
    control = el("textarea", { id, spellcheck: "false", rows: Math.min(12, Math.max(clase === "matriz" ? 5 : 3, lineas + 1)) });
    control.value = valores[campo.id] ?? "";
  }
  const guardar = () => {
    valores[campo.id] = control.value;
  };
  control.addEventListener("input", guardar);
  control.addEventListener("change", guardar);
  const poner = (valor) => {
    control.value = valor;
    guardar();
  };
  const atajos = [];
  if (clase === "secuencia" || clase === "secuencias") {
    atajos.push(atajo("Genoma completo de E. coli", () => poner("@ecoli")));
    atajos.push(campoGen("@gen:", "Gen de E. coli:", poner));
  }
  if (clase === "proteina") {
    atajos.push(campoGen("@proteina:", "Proteína de E. coli:", poner));
  }
  if (clase === "archivo" || clase === "base") {
    atajos.push(atajo("Elegir archivo…", async () => {
      const ruta = await elegirArchivo();
      if (ruta) {
        poner(ruta);
        if (clase === "base") {
          recordarBase(ruta);
        }
      }
    }));
  }
  if (clase === "base") {
    atajos.unshift(atajo("E. coli incluida", () => poner("@ecoli")));
    for (const ruta of basesRecientes()) {
      atajos.push(atajo(ruta.split(/[\\/]/).pop(), () => poner(ruta)));
    }
  }
  if (campo.ejemplo && !corto && !(clase === "base" && campo.ejemplo === "@ecoli")) {
    atajos.push(atajo("Usar el ejemplo", () => poner(campo.ejemplo)));
  }
  return {
    corto,
    nodo: el("div", { class: "campo" }, [
      el("label", { for: id, texto: campo.etiqueta }),
      control,
      atajos.length ? el("div", { class: "atajos" }, atajos) : null,
    ]),
  };
}

function pintarPliego() {
  const pliego = $("pliego");
  const ficha = fichaDe(estado.algoritmo);
  if (!ficha) {
    pliego.replaceChildren(el("div", { class: "vacio" }, [el("b", { texto: "Elige una estación en la carta" }), el("span", { texto: "Cada estación del río es un capítulo. El Puerto guarda las herramientas para datos reales." })]));
    return;
  }
  const valores = valoresDe(ficha);
  const largos = [];
  const cortos = [];
  for (const campo of ficha.campos) {
    const { corto, nodo } = controlDeCampo(ficha, campo, valores);
    (corto ? cortos : largos).push(nodo);
  }
  const ejecutar = el("button", { type: "submit", class: "boton principal", texto: "Ejecutar" });
  const tiempo = el("span", { class: "tiempo", "aria-live": "polite" });
  const ejemplo = atajo("Volver al ejemplo", () => {
    estado.valores.delete(ficha.id);
    pintarPliego();
  });
  const verPasos = ficha.paso_a_paso ? el("button", { type: "button", class: "boton", texto: "Explícame paso a paso" }) : null;
  const formulario = el("form", { class: "formulario" }, [
    ...largos,
    cortos.length ? el("div", { class: "fila-campos" }, cortos) : null,
    el("div", { class: "acciones" }, [ejecutar, verPasos, ejemplo, tiempo]),
  ]);
  const resultado = el("div", { class: "resultados" });
  const pasoAPaso = el("section", { class: "paso-a-paso", hidden: true, "aria-label": "Paso a paso" });
  const codigo = el("details", { class: "codigo-rust" }, [el("summary", {}, [document.createTextNode("Cómo funciona"), el("span", { class: "llano", texto: ficha.funcion })])]);
  let fragmentos = [];
  invocar("codigo", { id: ficha.id }).then((lista) => {
    fragmentos = lista;
    const total = lista.reduce((suma, f) => suma + f.fuente.split("\n").length, 0);
    codigo.querySelector(".llano").textContent = `${ficha.funcion}, ${total} líneas de Rust`;
    codigo.append(pintarCodigo(lista).caja);
  }).catch((error) => codigo.append(el("p", { class: "problema", texto: String(error) })));

  const entradaActual = () => ({ ...valores });
  formulario.addEventListener("submit", async (evento) => {
    evento.preventDefault();
    ejecutar.disabled = true;
    tiempo.textContent = "Calculando…";
    const inicio = performance.now();
    try {
      const salida = await invocar("ejecutar", { id: ficha.id, entrada: entradaActual() });
      const ms = performance.now() - inicio;
      estado.salida = salida;
      const duracion = ms < 1000 ? `${Math.max(1, Math.round(ms))} ms` : `${(ms / 1000).toFixed(1)} s`;
      tiempo.textContent = estado.nucleos ? `Listo en ${duracion}, calculado en tu equipo con ${estado.nucleos} núcleos` : `Listo en ${duracion}`;
      pintarSalida(resultado, salida);
    } catch (error) {
      tiempo.textContent = "";
      resultado.replaceChildren(el("p", { class: "problema", role: "alert", texto: String(error) }));
    } finally {
      ejecutar.disabled = false;
    }
  });
  formulario.addEventListener("keydown", (evento) => {
    if (evento.key === "Enter" && (evento.ctrlKey || evento.metaKey)) {
      evento.preventDefault();
      formulario.requestSubmit();
    }
  });
  if (verPasos) {
    verPasos.addEventListener("click", () => abrirPasoAPaso(pasoAPaso, ficha, entradaActual(), fragmentos));
  }
  const capitulo = capituloDe(ficha.capitulo);
  pliego.replaceChildren(
    el("header", { class: "pliego-cabeza" }, [
      el("p", { class: "capitulo", texto: rotuloCapitulo(capitulo) }),
      el("h2", { tabindex: -1, texto: ficha.nombre }),
      el("p", { class: "resumen", texto: ficha.resumen }),
      ficha.vanguardia ? el("p", { class: "referencia", texto: `Vanguardia. ${ficha.vanguardia}` }) : null,
    ]),
    formulario,
    resultado,
    pasoAPaso,
    codigo,
  );
}

function pintarSalida(contenedor, salida) {
  contenedor.replaceChildren(
    el("section", { class: "seccion" }, [el("h3", { texto: "Resultado" }), ...pintarBloques(salida.bloques)]),
    el("section", { class: "seccion" }, [el("h3", { texto: "Qué significa" }), el("ul", { class: "interpretacion" }, salida.interpretacion.map((t) => el("li", { texto: t })))]),
  );
}

function basesRecientes() {
  try {
    const lista = JSON.parse(recuperar("bases") || "[]");
    return Array.isArray(lista) ? lista.filter((ruta) => typeof ruta === "string").slice(0, 5) : [];
  } catch (error) {
    return [];
  }
}

function recordarBase(ruta) {
  const lista = [ruta, ...basesRecientes().filter((otra) => otra !== ruta)].slice(0, 5);
  recordar("bases", JSON.stringify(lista));
}

async function iniciar() {
  aplicarTema(recuperar("tema"));
  estado.vanguardia = recuperar("vanguardia") === "si";
  aplicarVanguardia();
  $("vanguardia").addEventListener("click", () => {
    estado.vanguardia = !estado.vanguardia;
    recordar("vanguardia", estado.vanguardia ? "si" : "no");
    aplicarVanguardia();
    const actual = fichaDe(estado.algoritmo);
    if (actual && actual.vanguardia && !estado.vanguardia) {
      const primera = visibles().find((f) => f.capitulo === estado.capitulo) || visibles()[0];
      if (primera) {
        elegirAlgoritmo(primera.id, false);
      }
    }
    pintarCarta();
    pintarIndice();
  });
  invocar("nucleos").then((cantidad) => {
    estado.nucleos = cantidad;
  }).catch(() => {});
  fetch("licencias.html").then((respuesta) => {
    if (respuesta.ok) {
      $("enlace-licencias").hidden = false;
    }
  }).catch(() => {});
  $("tema").addEventListener("click", () => {
    const nuevo = temaActual() === "noche" ? "dia" : "noche";
    recordar("tema", nuevo);
    aplicarTema(nuevo);
    pintarCarta();
    const resultado = document.querySelector(".resultados");
    if (estado.salida && resultado) {
      pintarSalida(resultado, estado.salida);
    }
  });
  $("buscar").addEventListener("input", (evento) => {
    estado.busqueda = evento.target.value.trim();
    pintarIndice();
  });
  try {
    const [capitulos, fichas] = await Promise.all([invocar("capitulos"), invocar("fichas")]);
    estado.capitulos = capitulos;
    estado.fichas = fichas;
  } catch (error) {
    $("pliego").replaceChildren(el("p", { class: "problema", texto: String(error) }));
    return;
  }
  const guardado = recuperar("algoritmo");
  const recordada = fichaDe(guardado);
  const inicial = recordada && (estado.vanguardia || !recordada.vanguardia) ? recordada : visibles()[0];
  if (inicial) {
    estado.capitulo = inicial.capitulo;
    estado.algoritmo = inicial.id;
  } else if (estado.capitulos.length) {
    estado.capitulo = estado.capitulos[0].id;
  }
  pintarCarta();
  pintarIndice();
  pintarPliego();
  invocar("genes_ecoli").then((genes) => {
    $("genes-ecoli").replaceChildren(...genes.map(([nombre, producto]) => el("option", { value: nombre, texto: producto })));
  }).catch(() => {});
}

iniciar();
