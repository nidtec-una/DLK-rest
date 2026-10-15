# DLK-rest: Solvers Krylov con Reinicio Adaptativo (DMDc-LQR) en PETSc

[![Language](https://img.shields.io/badge/Language-C11-blue.svg)](https://en.wikipedia.org/wiki/C11_(C_standard_revision))
[![Framework](https://img.shields.io/badge/Framework-PETSc%203.24+-red.svg)](https://petsc.org/)
[![MPI](https://img.shields.io/badge/Parallel-OpenMPI%20%2F%20MPICH-informational.svg)](https://www.open-mpi.org/)
[![Linear Algebra](https://img.shields.io/badge/BLAS%2FLAPACK-Optimized-success.svg)](http://www.netlib.org/lapack/)
[![License](https://img.shields.io/badge/License-MIT%20%2F%20Academic-lightgrey.svg)](LICENSE)

Implementación nativa, modular y de alto rendimiento en **PETSc (C / MPI)** de la familia de solucionadores de subespacios de Krylov con reinicio adaptativo basados en control óptimo y modelado reducido impulsado por datos, correspondiente a la investigación:

> **"Optimal Control for Restarting Minimal Residual Krylov Solvers via Data-Driven Residual Dynamics"**  
> *Carlos Benitez, Juan C. Cabral, Christian E. Schaerer.*

---

## 1. Motivación y Fundamento Teórico

El método de residuo mínimo generalizado con reinicio, **$\text{GMRES}(m)$**, enfrenta un dilema fundamental en sistemas lineales de gran escala $A x = b$:
* **Valores pequeños de $m$:** Bajo coste computacional y de memoria por ciclo, pero alto riesgo de **estancamiento severo** o convergencia lenta.
* **Valores grandes de $m$:** Mayor probabilidad de convergencia, pero a costa de un **crecimiento cuadrático en ortogonalizaciones** y almacenamiento en memoria ($m$ vectores de dimensión $N$).

### La Solución DLK (Dynamic Mode Decomposition with Control)
En lugar de fijar un parámetro estático $m$ o recurrir a heurísticas ad-hoc, el método **DLK**:
1. Recolecta una ventana deslizante de $p$ snapshots normalizados del residuo $r_k = b - A x_k$.
2. Aplica **DMDc (Dynamic Mode Decomposition with Control)** para identificar un modelo lineal en espacio de estados de orden reducido:
   $$\hat{x}_{k+1} = \bar{H} \hat{x}_k + \bar{B} u_k$$
   donde la variable de control $u_k = m_k$ es la dimensión del subespacio Krylov para el siguiente ciclo.
3. Resuelve una **Ecuación Algebraica de Riccati en Tiempo Continuo (CARE)** mediante el formalismo LQR sobre el subespacio controlable:
   $$\bar{H}^T P + P \bar{H} - P \bar{B} R^{-1} \bar{B}^T P + Q = 0, \quad K = R^{-1} \bar{B}^T P$$
4. Selecciona de forma óptima y adaptativa el parámetro de reinicio $m_{k+1}$ acotado en $[m_{\min}, m_{\max}]$, rompiendo mesetas de convergencia con el menor coste computacional.

El paquete implementa tanto la variante estándar **`DLK-GMRES`** como la variante enriquecida con aproximaciones de error acumuladas **`DLK-LGMRES`**.

---

## 2. Componentes del Repositorio

* **`src/` e `include/`**: Implementación nativa del solver KSP (`kspdlk.c`), rutinas de identificación DMDc y solución de Riccati (CARE) en C (`dmdc_lqr.c`), y encabezados públicos de PETSc (`petsckspdlk.h`).
* **`data/`**: Colección de 22 matrices SuiteSparse y sus binarios optimizados en `data/bin/`, junto con el script de conversión `convert_mat_to_bin.py`.
* **`scripts/`**: Automatización del benchmark comparativo (`run_benchmark_comparison.py`) y generación de figuras para publicación (`plot_benchmark_bars.py`, `plot_comparison_of_methods.py`, `plot_performance_with_m.py`).
* **`tests/`**: Driver de validación y benchmarking `test_dlk.c`.
* **`results/`**: Tablas experimentales en `results/tables/` (`matrix_info.csv`, `cycles_experiment.csv`, `avg_relativa_time_experiment.csv`) y figuras vectoriales en `results/figures/`.

---

## 3. Requisitos y Compilación

### Requisitos del Sistema
* **Compilador C:** GCC 9+ / Clang 11+ compatible con C11.
* **PETSc:** Versión 3.20 o superior configurada con soporte para números reales (`--with-scalar-type=real`).
* **Librerías BLAS/LAPACK:** OpenBLAS, Intel MKL o Netlib LAPACK.
* **MPI:** OpenMPI o MPICH.
* **Python 3:** (Para scripts de datos y figuras) `scipy`, `numpy`, `matplotlib`.

### Compilación
```bash
# Compilar el ejecutable test_dlk:
make clean && make
```
El `Makefile` detecta automáticamente la instalación de PETSc mediante `pkg-config` o variables de entorno del sistema (`PETSC_DIR` / `PETSC_ARCH`).

---

## 4. Guía de Uso del Solver

Los solvers se integran de forma nativa en el sistema de opciones de PETSc mediante el prefijo `-ksp_dlk_*`:

### Solvers Registrados (`-ksp_type`)
* **`dlk`**: GMRES con adaptación dinámica de $m$ vía DMDc-LQR o PD.
* **`dlk_lgmres`**: LGMRES adaptativo que combina subespacio Krylov de tamaño variable con $\ell$ vectores de enriquecimiento.

### Opciones de la Línea de Comandos

| Opción CLI | Tipo | Defecto | Descripción |
| :--- | :---: | :---: | :--- |
| `-matrix <path.bin>` | String | Requerido | Archivo binario PETSc con matriz $A$ y vector $b$. |
| `-ksp_type <tipo>` | String | `dlk` | Solver: `dlk`, `dlk_lgmres`, `gmres`, `lgmres`, `bcgs`. |
| `-ksp_dlk_controller` | Enum | `lqr` | Algoritmo de control: `lqr`, `pd`, `none`. |
| `-ksp_dlk_m0` | Entero | `30` | Dimensión de reinicio inicial ($m_0$). |
| `-ksp_dlk_m_min` | Entero | `2` | Cota mínima permitida para $m$. |
| `-ksp_dlk_m_max` | Entero | `60` | Cota máxima permitida para $m$. |
| `-ksp_dlk_alpha` | Real | `1e-4` | Penalización $\alpha$ en CARE ($R = \alpha I$). |
| `-ksp_dlk_aug` | Entero | `2` | Número de vectores de enriquecimiento $\ell$ (en `dlk_lgmres`). |
| `-pc_type <tipo>` | String | `ilu` | Precondicionador: `ilu`, `none`, `jacobi`, `bjacobi`. |
| `-tol <val>` | Real | `1e-8` | Tolerancia de convergencia relativa $\|r_k\|/\|r_0\|$. |
| `-cycle_file <path.csv>` | String | Opcional | Exporta historial de ciclos ($m$, residuo relativo, error DMDc). |

---

## 5. Ejemplos de Ejecución Rápida

### A. Ejecutar Sherman5 con DLK-GMRES Adaptativo (Sin Precondicionador)
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

### B. Ejecutar Cavity05 con DLK-LGMRES ($\ell = 2$)
```bash
./test_dlk -matrix data/bin/cavity05.bin \
           -ksp_type dlk_lgmres \
           -ksp_dlk_aug 2 \
           -ksp_dlk_m0 20 \
           -ksp_dlk_m_max 118 \
           -pc_type none \
           -tol 1e-8
```

### C. Ejecutar en Paralelo con MPI (4 procesos)
```bash
mpirun -np 4 ./test_dlk -m 80 -n 80 -ksp_type dlk_lgmres -ksp_dlk_controller lqr -pc_type bjacobi
```

---

## 6. Resultados Experimentales y Benchmarks

Evaluación comparativa integral de los cuatro solvers Krylov sobre el conjunto completo de **22 matrices de la colección SuiteSparse** bajo las mismas condiciones del artículo (`tol = 1e-9`, sin precondicionador, ejecución secuencial `OMP_NUM_THREADS=1`). Se denota como `NC` los métodos que no alcanzaron la tolerancia en el límite de ciclos.

### Tabla 1: Información y Propiedades de las Matrices de Prueba (`results/tables/matrix_info.csv`)

| Problem | Matrix | Dimensión ($n$) | No ceros ($nnz$) | Condicionamiento $\text{cond}(A)$ |
| :---: | :--- | :---: | :---: | :---: |
| 1 | `add20` | 2,395 | 13,151 | $1.7637 \times 10^4$ |
| 2 | `cavity05` | 1,182 | 32,632 | $9.1829 \times 10^5$ |
| 3 | `cavity10` | 2,597 | 76,171 | $4.4641 \times 10^6$ |
| 4 | `cdde1` | 961 | 4,681 | $4.0763 \times 10^3$ |
| 5 | `circuit2` | 4,510 | 21,199 | $7.1181 \times 10^6$ |
| 6 | `ex40` | 7,740 | 456,188 | $1.2019 \times 10^7$ |
| 7 | `fpgatrans01` | 1,220 | 7,382 | $2.8784 \times 10^4$ |
| 8 | `memplus` | 17,758 | 99,147 | $2.6661 \times 10^5$ |
| 9 | `orsirr1` | 1,030 | 6,858 | $1.6720 \times 10^5$ |
| 10 | `orsreg1` | 2,205 | 14,133 | $1.5394 \times 10^4$ |
| 11 | `pde2961` | 2,961 | 14,585 | $9.4909 \times 10^2$ |
| 12 | `raefsky1` | 3,242 | 293,409 | $3.1593 \times 10^4$ |
| 13 | `raefsky2` | 3,242 | 293,551 | $1.0823 \times 10^4$ |
| 14 | `rdb2048` | 2,048 | 12,032 | $2.0918 \times 10^3$ |
| 15 | `sherman1` | 1,000 | 3,750 | $2.2575 \times 10^4$ |
| 16 | `sherman4` | 1,104 | 3,786 | $7.1579 \times 10^3$ |
| 17 | `sherman5` | 3,312 | 20,793 | $3.9023 \times 10^5$ |
| 18 | `steam2` | 600 | 5,660 | $3.5481 \times 10^6$ |
| 19 | `wang2` | 2,903 | 19,093 | $3.1844 \times 10^4$ |
| 20 | `wang4` | 26,068 | 177,196 | $4.9125 \times 10^4$ |
| 21 | `watt1` | 1,856 | 11,360 | $5.3847 \times 10^9$ |
| 22 | `young3c` | 841 | 3,988 | $1.1532 \times 10^4$ |

### Tabla 2: Número de Ciclos de Reinicio hasta Convergencia (`results/tables/cycles_experiment.csv`)

| Problem | Matrix | GMRES(30) | LGMRES(28, 2) | DMDc–LQR–GMRES($m$) | DMDc–LQR–LGMRES($m, 2$) |
| :---: | :--- | :---: | :---: | :---: | :---: |
| 1 | `add20` | 35 | 21 | **8** | **9** |
| 2 | `cavity05` | NC | 80 | **9** | **11** |
| 3 | `cavity10` | NC | 127 | **28** | **12** |
| 4 | `cdde1` | 33 | 28 | 88 | **8** |
| 5 | `circuit2` | 25 | 20 | **7** | 23 |
| 6 | `ex40` | NC | NC | **83** | **58** |
| 7 | `fpgatrans01` | 104 | 35 | **14** | 28 |
| 8 | `memplus` | 93 | 49 | **8** | **7** |
| 9 | `orsirr1` | 158 | 78 | 52 | **19** |
| 10 | `orsreg1` | 23 | 19 | **18** | 43 |
| 11 | `pde2961` | 13 | 22 | 32 | 22 |
| 12 | `raefsky1` | 172 | 27 | 33 | **13** |
| 13 | `raefsky2` | 149 | 83 | 78 | **17** |
| 14 | `rdb2048` | 16 | **13** | 19 | 47 |
| 15 | `sherman1` | 123 | 27 | **18** | **19** |
| 16 | `sherman4` | 26 | 16 | 12 | **7** |
| 17 | `sherman5` | NC | NC | **52** | **20** |
| 18 | `steam2` | 14 | 23 | **13** | 17 |
| 19 | `wang2` | 93 | 29 | 24 | **15** |
| 20 | `wang4` | 94 | NC | **36** | **21** |
| 21 | `watt1` | 137 | 26 | 25 | **13** |
| 22 | `young3c` | 342 | 104 | 92 | **48** |

### Tabla 3: Tiempo Computacional Promedio en Segundos (`results/tables/avg_relativa_time_experiment.csv`)

| Problem | Matrix | GMRES(30) | LGMRES(28, 2) | DMDc–LQR–GMRES($m$) | DMDc–LQR–LGMRES($m, 2$) |
| :---: | :--- | :---: | :---: | :---: | :---: |
| 1 | `add20` | 0.0631 | 0.0459 | 0.0321 | 0.0192 |
| 2 | `cavity05` | NC | 0.1235 | 0.0162 | 0.0185 |
| 3 | `cavity10` | NC | 0.4472 | 0.3883 | 0.1917 |
| 4 | `cdde1` | 0.0239 | 0.0185 | 0.0261 | 0.0128 |
| 5 | `circuit2` | 0.0826 | 0.0595 | 0.0717 | 0.0964 |
| 6 | `ex40` | NC | NC | 83.7608 | 64.1836 |
| 7 | `fpgatrans01` | 0.1022 | 0.0323 | 0.0248 | 0.0263 |
| 8 | `memplus` | 1.4365 | 0.7388 | 0.8232 | 0.3724 |
| 9 | `orsirr1` | 0.1358 | 0.0611 | 0.0895 | 0.0725 |
| 10 | `orsreg1` | 0.0450 | 0.0330 | 0.0390 | 0.0365 |
| 11 | `pde2961` | 0.0293 | 0.0444 | 0.0323 | 0.0426 |
| 12 | `raefsky1` | 2.4398 | 0.3245 | 0.5384 | 0.3547 |
| 13 | `raefsky2` | 2.2124 | 1.0479 | 0.7282 | 0.6007 |
| 14 | `rdb2048` | 0.0255 | 0.0195 | 0.0243 | 0.0275 |
| 15 | `sherman1` | 0.0817 | 0.0185 | 0.0378 | 0.0261 |
| 16 | `sherman4` | 0.0216 | 0.0131 | 0.0138 | 0.0114 |
| 17 | `sherman5` | NC | NC | 3.0936 | 0.0268 |
| 18 | `steam2` | 0.0074 | 0.0114 | 0.0121 | 0.0085 |
| 19 | `wang2` | 0.2134 | 0.0558 | 0.0801 | 0.0905 |
| 20 | `wang4` | 2.1696 | NC | 3.8352 | 2.7772 |
| 21 | `watt1` | 0.1986 | 0.0335 | 0.0371 | 0.0709 |
| 22 | `young3c` | 0.2516 | 0.0635 | 0.0917 | 0.0943 |

---

## 7. Figuras Comparativas para Publicación

Todas las gráficas se generan automáticamente en formato vectorial EPS (para inclusión directa en LaTeX) y PNG a 300 DPI, siendo análogas a las del artículo científico:

### A. Gráficos de Barras Comparativos (Ciclos y Tiempo Computacional Relativo)
Generados automáticamente por [`scripts/plot_benchmark_bars.py`](file:///home/carlos_benitez/DLK-rest/scripts/plot_benchmark_bars.py) evaluando las razones logarítmicas de convergencia $\log_{10}(\text{GMRES}(30)/\text{Krylov Solver})$ con respecto a GMRES(30) como base (LGMRES, DMDc-LQR-GMRES y DMDc-LQR-LGMRES, excluyendo CG):

#### Reducción en Número de Ciclos
![Number of Cycles](results/figures/bar_plot_Number_of_cycles.png)

#### Reducción en Tiempo Computacional Relativo
![Relative Computation Time](results/figures/bar_plot_Average_computation_relative_time.png)

### B. Curvas de Convergencia Estilo SIAM (`sherman5`)
Generada con `python3 scripts/plot_comparison_of_methods.py`:

![Comparison of Methods sherman5](results/figures/Comparisonofmethodssherman5.png)

### C. Rendimiento Multietapa y Dinámica de $m$ (`sherman5`)
Generada con `python3 scripts/plot_performance_with_m.py`:

![Performance DLK-GMRES sherman5](results/figures/PerformanceofDMDcLQRGMRESmsherman5_with_m.png)

---

## 8. Reproducibilidad Automatizada

Para reproducir todos los experimentos y actualizar tablas y gráficas desde cero:

```bash
# 1. Compilar el binario C:
make clean && make

# 2. Convertir y verificar el dataset de matrices:
python3 data/convert_mat_to_bin.py

# 3. Ejecutar el benchmark comparativo completo:
python3 scripts/run_benchmark_comparison.py

# 4. Generar los gráficos de barras comparativos:
python3 scripts/plot_benchmark_bars.py

# 5. Generar las curvas de convergencia SIAM para sherman5:
python3 scripts/plot_comparison_of_methods.py
python3 scripts/plot_performance_with_m.py \
    --cycle_file results/tables/cycles_sherman5_dlk.csv \
    --output_prefix results/figures/PerformanceofDMDcLQRGMRESmsherman5_with_m \
    --m_max 331 --m_min 30
```

---

## 9. Cita y Referencia Académica

Si utilizas este software o los algoritmos DLK-GMRES / DLK-LGMRES en tus investigaciones, por favor cita:

```bibtex
@article{benitez2026dlk,
  title   = {Optimal Control for Restarting Minimal Residual Krylov Solvers via Data-Driven Residual Dynamics},
  author  = {Benitez, Carlos and Cabral, Juan C. and Schaerer, Christian E.},
  journal = {arXiv preprint},
  year    = {2026}
}
```
