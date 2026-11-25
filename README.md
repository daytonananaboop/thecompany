# La Compañía

Algoritmos de bioinformática escritos en Rust, en una aplicación de escritorio en español. Cada algoritmo muestra su resultado, una explicación de qué significa y el código que lo calculó. No usa inteligencia artificial ni internet: el genoma, los genes y las proteínas de *Escherichia coli* K-12 MG1655 vienen incluidos.

- Capítulos I–VIII: replicación, motivos, ensamblaje, péptidos, alineamiento, reordenamientos, filogenia y agrupamiento.
- El Puerto: simulación de lecturas, calidad y recorte, búsqueda tipo BLAST, mapeo, llamado de variantes, predicción de genes y diseño de cebadores.

## Instalar

Descarga el instalador de tu sistema en [Releases](../../releases):

- Windows: `La.Compania_x.y.z_x64-setup.exe`. Si aparece «Windows protegió tu PC», elige «Más información» y luego «Ejecutar de todas formas».
- macOS: `La.Compania_x.y.z_universal.dmg`. La primera vez, abre la aplicación con clic derecho y «Abrir».
- Linux: `.deb` o `.AppImage`.

## Compilar

```bash
cargo test
```

```bash
cargo install tauri-cli --version "^2"
```

```bash
cd app && cargo tauri build
```

En Linux hacen falta `libwebkit2gtk-4.1-dev`, `libsoup-3.0-dev` y `librsvg2-dev`.

## Licencia

Código bajo CC0 1.0. Fuentes tipográficas bajo SIL Open Font License (`app/ui/fuentes`). Las licencias de las bibliotecas de terceros vienen dentro de la aplicación. Datos de NCBI RefSeq NC_000913.3.

Uso educativo. No sirve para decisiones clínicas ni diagnósticas.
