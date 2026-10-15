# Guía de Ejecución y Resultados de Benchmarks (DLK-PETSc)

Este directorio contiene los resultados experimentales, figuras para publicación y tablas comparativas generadas por la suite de solvers **DLK-GMRES** y **DLK-LGMRES** implementados sobre **PETSc**.

---

## 1. Organización del Directorio

```text
results/
├── README.md                      <- Este documento de referencia y guía de ejecución
├── figures/                       <- Gráficas de publicación en alta resolución (PNG 300 DPI y EPS vectorial)
│   ├── benchmark_bars_iterations.png / .eps  <- Gráfico de barras comparativo de iteraciones
│   ├── benchmark_bars_time.png / .eps        <- Gráfico de barras comparativo de tiempo de CPU
│   ├── Comparisonofmethodssherman5.png / .eps<- Curvas de convergencia de métodos (estilo SIAM)
│   ├── PerformanceofDMDcLQRGMRESmsherman5_with_m.png / .eps <- Rendimiento en 3 paneles (Residuo, Error DMDc, m)
│   └── ...
└── tables/                        <- Datos tabulares en formato CSV y Markdown
    ├── benchmark_comparison.csv   <- Tabla de rendimiento entre solvers en toda la suite
    ├── benchmark_comparison.md    <- Resumen formateado en Markdown
    ├── cycles_sherman5_dlk.csv    <- Trayectoria de ciclos, residuo, error DMDc y m
    └── ...
```

---

## 2. Cómo Ejecutar los Solvers con Cualquier Matriz `.bin`

El driver ejecutable [`test_dlk`](file:///home/carlos_benitez/DLK-rest/test_dlk) permite cargar cualquier sistema lineal de `data/bin/` y resolverlo con cualquier solver nativo o adaptativo de PETSc.

### Sintaxis Básica
```bash
./test_dlk -matrix data/bin/<nombre_matriz>.bin [opciones_ksp] [opciones_dlk] [opciones_pc]
```

### Opciones Principales de la Línea de Comandos

| Parámetro | Tipo | Valor por Defecto | Descripción |
| :--- | :---: | :---: | :--- |
| `-matrix <path>` | String | Requerido | Ruta al archivo binario de la matriz (`data/bin/<nombre>.bin`) |
| `-ksp_type <tipo>` | String | `dlk` | Tipo de solver: `dlk`, `dlk_lgmres`, `gmres`, `lgmres`, `bcgs` |
| `-tol <val>` / `-ksp_rtol <val>` | Real | `1e-8` | Tolerancia de convergencia para el residuo relativo |
| `-ksp_max_it <n>` | Entero | `10000` | Número máximo de iteraciones permitidas |
| `-pc_type <tipo>` | String | `ilu` | Precondicionador: `ilu` (ILU(0)), `none` (sin precondicionador), `jacobi` |
| `-cycle_file <path.csv>` | String | Opcional | Exporta el historial de ciclos exteriores ($m$, residuo relativo, error DMDc) |
| `-history_file <path.csv>` | String | Opcional | Exporta el historial paso a paso de cada producto matriz-vector |

### Opciones Específicas del Solver DLK (`-ksp_type dlk` o `dlk_lgmres`)

| Opción DLK | Tipo | Defecto | Descripción |
| :--- | :---: | :---: | :--- |
| `-ksp_dlk_controller` | Enum | `lqr` | Tipo de controlador adaptativo: `lqr`, `pd`, `none` |
| `-ksp_dlk_m0` | Entero | `30` | Dimensión inicial de reinicio ($m_0$) |
| `-ksp_dlk_m_min` | Entero | `2` | Cota mínima permitida para $m$ |
| `-ksp_dlk_m_max` | Entero | `60` | Cota máxima permitida para $m$ |
| `-ksp_dlk_alpha` | Real | `1e-4` | Parámetro de penalización $\alpha$ en la ecuación CARE de LQR |
| `-ksp_dlk_aug` | Entero | `2` | Número de vectores de enriquecimiento $l$ (en `dlk_lgmres`) |
| `-ksp_dlk_p` | Entero | `5` | Ventana de snapshots para DMDc |

---

## 3. Ejemplos Prácticos de Ejecución

### Ejemplo 1: DLK-GMRES Adaptativo (LQR) en `sherman5`
```bash
./test_dlk -matrix data/bin/sherman5.bin \
           -ksp_type dlk \
           -ksp_dlk_controller lqr \
           -ksp_dlk_m0 30 \
           -ksp_dlk_m_min 30 \
           -ksp_dlk_m_max 331 \
           -ksp_dlk_alpha 0.001 \
           -pc_type none \
           -tol 1e-12 \
           -cycle_file results/tables/cycles_sherman5_dlk.csv
```

### Ejemplo 2: DLK-LGMRES con Enriquecimiento ($l=2$) en `cavity05`
```bash
./test_dlk -matrix data/bin/cavity05.bin \
           -ksp_type dlk_lgmres \
           -ksp_dlk_aug 2 \
           -ksp_dlk_m0 20 \
           -ksp_dlk_m_max 118 \
           -pc_type none \
           -tol 1e-8
```

### Ejemplo 3: Comparación Rápida con GMRES Estándar
```bash
./test_dlk -matrix data/bin/sherman1.bin -ksp_type gmres -ksp_gmres_restart 30 -pc_type ilu -tol 1e-8
```

---

## 4. Generación de Tablas y Gráficos Comparativos

### Ejecutar el Benchmark Completo
Para ejecutar automáticamente la batería de solvers en el conjunto de matrices de prueba y actualizar las tablas:
```bash
python3 scripts/run_benchmark_comparison.py
```
Esto genera:
* `results/tables/benchmark_comparison.csv`
* `results/tables/benchmark_comparison.md`

### Generar los Gráficos de Barras Comparativos
Para construir los gráficos de barras a partir del archivo CSV:
```bash
python3 scripts/plot_benchmark_bars.py
```
Salidas generadas:
* `results/figures/benchmark_bars_iterations.png` (y `.eps` vectorial para LaTeX)
* `results/figures/benchmark_bars_time.png` (y `.eps` vectorial para LaTeX)

### Generar las Figuras Estilo SIAM de `sherman5`
```bash
# Curvas comparativas de convergencia:
python3 scripts/plot_comparison_of_methods.py

# Gráfica de 3 paneles con la dinámica adaptativa de m:
python3 scripts/plot_performance_with_m.py \
    --cycle_file results/tables/cycles_sherman5_dlk.csv \
    --output_prefix results/figures/PerformanceofDMDcLQRGMRESmsherman5_with_m \
    --m_max 331 --m_min 30
```
