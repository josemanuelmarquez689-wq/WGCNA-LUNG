# directed_attack_network

Robustez global de una red de co-expresión WGCNA (24 módulos) ante **ataques dirigidos a genes hub**.

Pregunta: *cuando se quitan los hubs de un módulo, ¿la red completa sigue unida y comunicada? ¿El daño es mayor que quitar genes cualesquiera?*

## Estructura

```
directed_attack_network/
├── data/
│   └── datExpr_colores_24Modulos.RData   -> liga simbólica a ../WGCNA/datExpr_colores_24Modulos.RData
├── scripts/
│   ├── 01_robustez_global_integridad.R   # robustez GLOBAL: ataques dirigidos + control al azar de toda la red
│   ├── 02_control_mismo_modulo.R         # robustez GLOBAL: control estricto, azar dentro del mismo módulo
│   ├── 03_conservacion_modulo_jaccard.R  # robustez DEL MÓDULO: re-clustering + Jaccard máximo
│   └── 04_figuras_red.R                  # dendrograma, red de genes y red de módulos
├── run_all.sh                            # corre 04 -> 03 -> 01 -> 02 y junta las figuras
├── figures/                              # figuras finales (Fig1-Fig9)
├── results/                              # tablas y salidas completas de cada script
└── logs/                                 # logs de run_all.sh
```

## Datos

`datExpr_colores_24Modulos.RData` contiene:

- `datExpr_5000`: matriz de expresión, 276 muestras × 5000 genes.
- `colores_base`: módulo asignado a cada gen. Son 24 módulos más `grey` (2840 genes sin módulo). El más grande es turquoise (265) y el más chico darkgrey (30).

Parámetros de la red: `power = 12`, red `signed`, `TOMType = "signed"`.

> **Ojo con el power:** `../WGCNA/Ataque_masivo_2.R` dice `POWER_OPTIMO <- 9`, pero eso corresponde a otra red. Los resultados de esta red de 24 módulos se reproducen exactamente con **power 12** (BLACK 1H = 0.244318) y no con 9 (0.4645).

## Método

### Ataques (24 módulos × 19 ataques = 456)
Los genes de cada módulo se ordenan por conectividad intramodular (`kWithin`) y se quitan:
- `1 H` … `15 H`: los 1 a 15 hubs principales.
- `5%`, `10%`, `20%`, `90%`: ese porcentaje del módulo, empezando por los hubs.

### Robustez del módulo (script 03): Jaccard máximo
1. Se quitan los hubs y se reconstruye la red desde cero: adyacencia, TOM, `hclust` (average) y `cutreeDynamic` (deepSplit 2, minClusterSize 30).
2. Se calcula una matriz de Jaccard entre los módulos originales y los de la red alterada, usando solo los genes sobrevivientes y sin `grey`.
3. **Conservación del módulo atacado = Jaccard máximo de su fila**, es decir, el módulo nuevo con el que mejor coincide. Se usa el máximo porque al re-clusterizar los colores no corresponden a los originales.
4. Criticidad = mediana de la conservación en hubs 3–7: < 0.3 extremadamente crítico, < 0.5 crítico, < 0.7 moderado, ≥ 0.7 robusto.

Este script reproduce exactamente `Matriz_Conservacion_Real.csv` del backup validado (`../WGCNA/Heatmaps_Robustez_NuevaRed_BACKUP_24Modulos_20260923_095150/`).

### Robustez global (scripts 01 y 02): red y métricas
1. Grafo: dos genes están conectados si su TOM es ≥ **0.0226**, el umbral del top 1% de pares de la red original. Con ese umbral la red original tiene eficiencia global 0.164, un componente gigante de 3723 genes y 49 componentes.
2. En cada ataque se recalcula la TOM sin los genes removidos y se aplica **el mismo umbral**.
3. Métricas, relativas a la red original y normalizadas por el N original:
   - **Eficiencia global relativa**: qué tan cortos siguen siendo los caminos entre genes (1 = igual que la original).
   - **Componente gigante relativo**: tamaño del bloque principal conectado después del ataque / 3723. Compara tamaños; no es un Jaccard.

### Controles (30 repeticiones cada uno)
- **Script 01, azar de toda la red**: quitar el mismo número de genes al azar entre los 5000. Es un control permisivo, porque los hubs siempre están más conectados que un gen promedio.
- **Script 02, azar del mismo módulo** (control estricto): quitar el mismo número de genes al azar **del módulo atacado**. Responde si los hubs de X dañan más que cualquier gen de X.
- Se reportan como z = (media_azar − ataque) / sd_azar y como razón de daño = (1 − Ef_hubs) / (1 − Ef_azar). Se usa el valor-p empírico; `*` en los heatmaps = peor que las 30 repeticiones.

## Cómo correr

Todo de una vez (~1.5 h; los logs quedan en `logs/`):

```bash
cd /STORAGE/genut/Manuel/directed_attack_network
bash run_all.sh
```

O cada script por separado, desde la raíz del proyecto. El 02 necesita las salidas del 01.

```bash
/opt/R/current/bin/Rscript scripts/04_figuras_red.R                                       # ~1 min
NUCLEOS=20 /opt/R/current/bin/Rscript scripts/03_conservacion_modulo_jaccard.R            # ~5 min (456 re-clusterings)
NUCLEOS=20 N_REP=30 /opt/R/current/bin/Rscript scripts/01_robustez_global_integridad.R   # ~12 min (1806 redes)
NUCLEOS=24 N_REP=30 /opt/R/current/bin/Rscript scripts/02_control_mismo_modulo.R         # ~75 min (11910 redes)
```

Prueba rápida, para verificar que todo corre (pocos módulos y pocas repeticiones; escribe en `results/prueba_*`):

```bash
PRUEBA=1 NUCLEOS=6 /opt/R/current/bin/Rscript scripts/03_conservacion_modulo_jaccard.R
PRUEBA=1 N_REP=4 NUCLEOS=8 /opt/R/current/bin/Rscript scripts/01_robustez_global_integridad.R
PRUEBA=1 N_REP=4 NUCLEOS=8 /opt/R/current/bin/Rscript scripts/02_control_mismo_modulo.R
```

Las pruebas escriben en `results/prueba_*`. Requisitos: paquetes de R `WGCNA`, `igraph`, `pheatmap` y `parallel` (librería en `~/R/x86_64-pc-linux-gnu-library/4.6`).

## Salidas

### Figuras (`figures/`)

| Figura | Contenido |
|---|---|
| `Fig1_Dendrograma_Modulos.pdf` | Dendrograma de genes (1 − TOM) con la barra de colores de módulo |
| `Fig2_Red_Genes.pdf` / `.png` | Grafo de los genes en módulos (aristas TOM top 1%), coloreado por módulo, con el hub principal de cada módulo resaltado y anotado con su símbolo |
| `Fig3_Red_Modulos.pdf` | Red de módulos: tamaño = n genes, grosor = densidad de aristas fuertes entre módulos |
| `Fig4_Conservacion_Modulo_Jaccard.pdf` | Robustez del módulo: Jaccard máximo por ataque, con criticidad |
| `Fig5_Eficiencia_Global_Relativa.pdf` | Robustez global: eficiencia relativa |
| `Fig6_Componente_Gigante_Relativo.pdf` | Robustez global: tamaño del componente gigante relativo (no es Jaccard) |
| `Fig7_Dano_vs_Azar_Red.pdf` | z contra azar de toda la red |
| `Fig8_Dano_vs_Azar_MismoModulo.pdf` | z contra azar del mismo módulo (control estricto) |
| `Fig9_Razon_Dano_Hubs_vs_MismoModulo.pdf` | Veces más daño de los hubs que de genes al azar del módulo |

En la Fig1 los colores no forman bloques perfectamente continuos bajo el árbol. Es esperable: el árbol se calcula con `hclust` sobre la TOM completa, mientras que la asignación de módulos viene de la reconstrucción original de la red.

### Tablas

`results/red/`: `Hubs_Top5_Por_Modulo.csv` (top 5 hubs por módulo con símbolo, kWithin y kTotal) y las densidades de aristas internas y entre módulos.

`results/conservacion/`: `Matriz_Conservacion_Real.csv`, `Conservacion_Largo.csv` (incluye el módulo de la red alterada con el que coincide mejor), `Jaccard_Maximo_Por_Modulo.csv`, `Criticidad_Modulos.csv` y `Matrices_Jaccard_Por_Ataque.rds` (las 456 matrices completas).

`results/integridad/`:

| Archivo | Contenido |
|---|---|
| `Heatmap_Eficiencia_Global_Relativa.pdf` | Eficiencia global después del ataque / original |
| `Heatmap_Componente_Gigante_Relativo.pdf` | Tamaño del componente gigante después del ataque / original |
| `Heatmap_Dano_vs_Azar_Eficiencia.pdf` | z contra azar de toda la red (color saturado en ±10) |
| `Heatmap_Dano_vs_Azar_MismoModulo.pdf` | z contra azar del mismo módulo (control estricto) |
| `Heatmap_Razon_Dano_Hubs_vs_MismoModulo.pdf` | Veces más daño de los hubs que de genes al azar del módulo |
| `Ranking_Robustez_Global.csv`, `Ranking_Control_MismoModulo.csv` | Resumen por módulo |
| `Resultados_*_Largo.csv`, `Controles_*.csv` | Datos por ataque y por repetición |

## Resultados principales (corrida del 12-oct-2026)

- **BROWN** es el módulo que más daña la red: su hub principal (1H) hace 27 veces más daño que un gen al azar de BROWN, con z = 31 contra azar de toda la red y z = 13 contra azar del mismo módulo. Es significativo en las 19 columnas.
- **PURPLE** tiene un "super-hub": el 1H hace 26 veces más daño que un gen al azar de PURPLE; a partir del segundo hub la razón baja a 2–4.
- **RED, BLUE, YELLOW y PINK**: sus hubs hacen ~2.5–3 veces más daño que genes al azar de su módulo.
- **TURQUOISE** importa como bloque, no por sus hubs (razón 1.7): sus aristas fuertes son redundantes dentro del módulo.
- **Robustos**: ROYALBLUE, TAN, DARKRED, LIGHTGREEN, CYAN y GREENYELLOW (razón ≈ 1–1.4).

### Cuidados de interpretación
- **Columna 90%**: el azar casi no varía, así que el z se infla aunque la razón sea ~1.1. Ahí conviene leer la razón, no el z.
- **Columnas 3H, 7H y 11H** del control de toda la red: el z se infla porque en esos tamaños el azar varió poco. Conviene leer la fila completa o el ranking.
- **Métricas descartadas**: Jaccard_Global, ARI, correlación de la TOM completa y Jaccard de aristas. Las tres primeras se diluyen al promediar sobre toda la red. El Jaccard de aristas mide sobre todo lo que se quita, no la respuesta de la red. Esos análisis siguen en `../WGCNA/`.
