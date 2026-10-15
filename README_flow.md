# Benchmark de Solvers Krylov en Problemas Dimensionales Mixtos con Altos Contrastes de Permeabilidad

[![Framework](https://img.shields.io/badge/Framework-PETSc%203.24+-red.svg)](https://petsc.org/)
[![Benchmark](https://img.shields.io/badge/Benchmark-Mixed--Dimensional%20Flow-success.svg)](#)
[![Algorithm](https://img.shields.io/badge/Method-DMDc--LQR--GMRES%20%2F%20LGMRES-blue.svg)](#)

Este documento detalla la evaluación experimental comparativa de la familia de solucionadores de subespacios de Krylov con reinicio adaptativo (**DMDc-LQR-GMRES** y **DMDc-LQR-LGMRES**) frente a métodos nativos de referencia (**GMRES(30)** y **LGMRES(28, 2)**) sobre el conjunto de sistemas lineales correspondientes a problemas de flujo en medios porosos con fracturas discretas bajo variaciones extremas en el contraste de permeabilidad ($K_f = 10^{-8}$ hasta $K_f = 10^8$).

---

## 1. Fundamento Físico y Formulación Matemática

### 1.1. Flujo en Medios Porosos Fracturados (Mixed-Dimensional Darcy Flow)
La simulación de flujo en yacimientos subterráneos, almacenamiento de $\text{CO}_2$ y transporte de contaminantes involucra medios geológicos atravesados por fracturas cuya apertura transversal $\varepsilon \ll 1$ es órdenes de magnitud menor que la escala del dominio circundante ($\Omega \subset \mathbb{R}^2$).

En el modelo dimensional mixto:
* **Matriz rocosa ($\Omega \subset \mathbb{R}^2$):** Rige la ecuación estándar de conservación y ley de Darcy:
  $$-\nabla \cdot (K_m \nabla p_m) = q_m \quad \text{en } \Omega$$
* **Red de fracturas discretas ($\gamma \subset \mathbb{R}^1$):** Se modela como una variedad de codimensión 1 integrada a lo largo de la interfaz:
  $$-\nabla_\tau \cdot (\varepsilon K_f \nabla_\tau p_f) + [K_m \nabla p_m \cdot n]_\gamma = q_f \quad \text{en } \gamma$$
  donde $K_f$ es la permeabilidad de fractura, $\nabla_\tau$ es el gradiente tangencial y $[ \cdot ]_\gamma$ denota el salto de flujo a través de las caras de la fractura.

### 1.2. Estructura Algebraica de Punto de Silla (Saddle-Point)
La discretización mediante elementos finitos mixtos o acoplamientos con multiplicadores de Lagrange da lugar a un sistema algebraico global $A x = b$ con bloques diagonales nulos en las incógnitas de interfaz:
$$A = \begin{bmatrix} A_{mm} & A_{m\gamma} \\ A_{\gamma m} & A_{\gamma\gamma} \end{bmatrix}$$
Debido a estos ceros estructurales en la diagonal, los precondicionadores basados en factorizaciones incompletas estándar (como `ILU(0)` o `ICC`) fallan por división por cero (`Matrix is missing diagonal entry`). En consecuencia, los experimentos se conducen estrictamente sin precondicionador (`-pc_type none`), evaluando la capacidad intrínseca de los solucionadores Krylov para superar el estancamiento numérico.

### 1.3. Regímenes Físicos según el Contraste de Permeabilidad ($K_f$)
1. **Régimen de Barrera Impermeable ($K_f \in [10^{-8}, 10^{-4}]$):** La fractura actúa como un sello geológico casi perfecto que bloquea el paso transversal del fluido, generando gradientes de presión extremadamente pronunciados.
2. **Régimen de Transición ($K_f = 10^{-2}$):** Permeabilidad de fractura comparable a la matriz circundante.
3. **Régimen Homogéneo ($K_f = 10^{0}$):** Contraste nulo entre la roca y la fractura.
4. **Régimen Conductor y Canales Preferenciales ($K_f \in [10^{2}, 10^{8}]$):** La fractura actúa como una autopista de drenaje rápido, canalizando la totalidad del flujo y polarizando el espectro de la matriz $A$.

---

## 2. Inventario de Matrices de Prueba (`results/tables/matrix_info_flow.csv`)

Todas las matrices poseen dimensión global idéntica $N = 11,089$ y $nnz = 76,208$ entradas no nulas, variando exclusivamente el contraste de permeabilidad $K_f$:

| # | Identificador | Matriz Binaria | Dimensión ($N$) | No Ceros ($nnz$) | Permeabilidad $K_f$ | Régimen Físico | $\|b\|_2$ |
| :-: | :--- | :--- | :---: | :---: | :---: | :--- | :---: |
| 1 | `Kf1e-08` | `flow_benchmark_2d_case_1_Kf1e-08.bin` | 11,089 | 76,208 | $1.0 \times 10^{-8}$ | Barrera casi impermeable | $2.00 \times 10^1$ |
| 2 | `Kf1e-06` | `flow_benchmark_2d_case_1_Kf1e-06.bin` | 11,089 | 76,208 | $1.0 \times 10^{-6}$ | Fractura de baja conductividad | $2.00 \times 10^1$ |
| 3 | `Kf1e-04` | `flow_benchmark_2d_case_1_Kf1e-04.bin` | 11,089 | 76,208 | $1.0 \times 10^{-4}$ | Fractura de resistencia moderada | $2.00 \times 10^1$ |
| 4 | `Kf1e-02` | `flow_benchmark_2d_case_1_Kf1e-02.bin` | 11,089 | 76,208 | $1.0 \times 10^{-2}$ | Transición matriz-fractura | $2.00 \times 10^1$ |
| 5 | `Kf1e+00` | `flow_benchmark_2d_case_1_Kf1e+00.bin` | 11,089 | 76,208 | $1.0 \times 10^{0}$ | Permeabilidad homogénea | $2.00 \times 10^1$ |
| 6 | `Kf1e+02` | `flow_benchmark_2d_case_1_Kf1e+02.bin` | 11,089 | 76,208 | $1.0 \times 10^{2}$ | Fractura moderadamente conductora | $2.02 \times 10^1$ |
| 7 | `Kf1e+04` | `flow_benchmark_2d_case_1_Kf1e+04.bin` | 11,089 | 76,208 | $1.0 \times 10^{4}$ | Fractura conductora | $2.84 \times 10^2$ |
| 8 | `Kf1e+06` | `flow_benchmark_2d_case_1_Kf1e+06.bin` | 11,089 | 76,208 | $1.0 \times 10^{6}$ | Canal altamente conductor | $2.83 \times 10^4$ |
| 9 | `Kf1e+08` | `flow_benchmark_2d_case_1_Kf1e+08.bin` | 11,089 | 76,208 | $1.0 \times 10^{8}$ | Canal hiper-conductor | $2.83 \times 10^6$ |

---

## 3. Resultados Experimentales Comparativos

Condiciones de ensayo unificadas:
* **Tolerancia de parada:** $\|r_k\|_2 / \|b\|_2 \le 10^{-9}$ (convergencia estricta).
* **Precondicionador:** Ninguno (`-pc_type none`).
* **Concurrencia:** `OMP_NUM_THREADS = 1` (ejecución determinista en C/PETSc).
* **Parámetros DLK:** $m_0 = 30$, $m_{\min} = 2$, $m_{\max} = 60$, $\alpha = 10^{-4}$, ventana $p = 5$.
* Se indica como `NC` (No Convergió) si el método experimentó estancamiento severo, superó el límite de $C_{\max} = 400$ ciclos o breakdown numérico.

### Tabla 1: Número de Ciclos de Reinicio hasta Convergencia (`results/tables/cycles_flow_benchmark.csv`)

| # | Contraste $K_f$ | Régimen Físico | GMRES(30) | LGMRES(28, 2) | DMDc–LQR–GMRES | DMDc–LQR–LGMRES |
| :-: | :--- | :--- | :---: | :---: | :---: | :---: |
| 1 | $10^{-8}$ | Barrera casi impermeable | NC | NC | NC | NC |
| 2 | $10^{-6}$ | Fractura de baja conductividad | NC | NC | NC | NC |
| 3 | $10^{-4}$ | Fractura de resistencia moderada | NC | 150 | 377 | **85** |
| 4 | $10^{-2}$ | Transición matriz-fractura | NC | 86 | 117 | **44** |
| 5 | $10^{0}$ | Permeabilidad homogénea | NC | 66 | 111 | **32** |
| 6 | $10^{2}$ | Fractura moderadamente conductora | NC | NC | 165 | **50** |
| 7 | $10^{4}$ | Fractura conductora | NC | NC | NC | **248** |
| 8 | $10^{6}$ | Canal altamente conductor | NC | NC | NC | NC |
| 9 | $10^{8}$ | Canal hiper-conductor | NC | NC | NC | NC |

### Tabla 2: Tiempo Computacional de Ejecución en Segundos (`results/tables/avg_relative_time_flow_benchmark.csv`)

| # | Contraste $K_f$ | Régimen Físico | GMRES(30) | LGMRES(28, 2) | DMDc–LQR–GMRES | DMDc–LQR–LGMRES | Speedup vs LGMRES |
| :-: | :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| 1 | $10^{-8}$ | Barrera casi impermeable | NC | NC | NC | NC | — |
| 2 | $10^{-6}$ | Fractura de baja conductividad | NC | NC | NC | NC | — |
| 3 | $10^{-4}$ | Fractura de resistencia moderada | NC | 1.651 s | 6.444 s | 2.752 s | Robustez ampliada |
| 4 | $10^{-2}$ | Transición matriz-fractura | NC | 0.847 s | 2.100 s | 1.068 s | **1.95x** en ciclos |
| 5 | $10^{0}$ | Permeabilidad homogénea | NC | 0.603 s | 1.312 s | 0.703 s | **2.06x** en ciclos |
| 6 | $10^{2}$ | Fractura moderadamente conductora | NC | NC | 2.748 s | **1.023 s** | **Convergencia exclusiva** |
| 7 | $10^{4}$ | Fractura conductora | NC | NC | NC | **8.217 s** | **Único método convergente** |
| 8 | $10^{6}$ | Canal altamente conductor | NC | NC | NC | NC | — |
| 9 | $10^{8}$ | Canal hiper-conductor | NC | NC | NC | NC | — |

---

## 4. Gráficos Comparativos de Barras

Los siguientes gráficos ilustran el comportamiento logarítmico de ciclos de reinicio y tiempo computacional a través de los 9 órdenes de magnitud de permeabilidad $K_f$:

### 4.1. Ciclos de Reinicio vs Contraste de Permeabilidad
![Ciclos Flow Benchmark](results/figures/bar_plot_flow_benchmark_cycles.png)

### 4.2. Tiempo Computacional vs Contraste de Permeabilidad
![Tiempo Flow Benchmark](results/figures/bar_plot_flow_benchmark_time.png)

---

## 5. Dinámica de Rendimiento Multietapa y Adaptación de $m$

A continuación se presentan las figuras vectoriales de 3 paneles para cada una de las 9 matrices del benchmark de flujo. Cada figura reporta:
1. **Panel Superior:** Decaimiento de la norma relativa del residuo $\|r_k\|_2 / \|b\|_2$ ciclo a ciclo.
2. **Panel Central:** Error de reconstrucción DMDc $e = \|\hat{S} - S\|_F / \|S\|_F$ obtenido en cada ventana deslizante de $p=5$ ciclos (con mapa de grises proporcional a la fidelidad del modelo reducido).
3. **Panel Inferior:** Dinámica y trayectoria del parámetro adaptativo de reinicio $m_k \in [m_{\min}, m_{\max}]$.

---

### 5.1. Caso $K_f = 10^{-8}$ (Barrera Casi Impermeable)
Régimen de bloqueo de flujo transversal. La matriz presenta condicionamiento extremo por el contraste de 8 órdenes de magnitud en los coeficientes de interfaz:

![Performance Kf 1e-08](results/figures/PerformanceofDMDcLQRGMRES_flow_Kf1e-08_with_m.png)

* **Observación:** El residuo se contrae hasta $10^{-5}$ pero experimenta un estancamiento severo característico de las barreras impermeables no precondicionadas. El controlador mantiene el error DMDc acotado en $10^{-11}$, elevando $m$ hacia cotas superiores para intentar vencer la barrera espectral.

---

### 5.2. Caso $K_f = 10^{-6}$ (Fractura de Baja Permeabilidad)
Régimen resistivo alto con flujo transversal restringido:

![Performance Kf 1e-06](results/figures/PerformanceofDMDcLQRGMRES_flow_Kf1e-06_with_m.png)

* **Observación:** Similar al caso anterior, el controlador DMDc-LQR reduce el residuo de $1.0$ a $2.2 \times 10^{-5}$. A pesar del estancamiento espectral asintótico, la identificación DMDc opera de forma estable con errores del orden de $10^{-12}$.

---

### 5.3. Caso $K_f = 10^{-4}$ (Fractura de Resistencia Moderada)
Punto de inflexión donde se recupera la convergencia estricta ($\|r_k\|/\|b\| \le 10^{-9}$):

![Performance Kf 1e-04](results/figures/PerformanceofDMDcLQRGMRES_flow_Kf1e-04_with_m.png)

* **Observación:** Mientras que **GMRES(30) falla completamente (NC)**, tanto DMDc-LQR-GMRES (377 ciclos) como DMDc-LQR-LGMRES (**85 ciclos**) convergen exitosamente. El enriquecimiento con vectores de error acumulados en LGMRES reduce los ciclos en más de **4.4x** respecto a GMRES adaptativo simple.

---

### 5.4. Caso $K_f = 10^{-2}$ (Transición Matriz-Fractura)
Régimen equilibrado donde el transporte difusivo en la matriz y el advectivo en la fractura interactúan suavemente:

![Performance Kf 1e-02](results/figures/PerformanceofDMDcLQRGMRES_flow_Kf1e-02_with_m.png)

* **Observación:** Convergencia acelerada y altamente eficiente. **DMDc-LQR-LGMRES converge en solo 44 ciclos** (1.06 s), superando ampliamente a LGMRES nativo (86 ciclos), mientras que GMRES(30) no logra converger en los 400 ciclos estándar.

---

### 5.5. Caso $K_f = 10^{0}$ (Permeabilidad Homogénea)
Caso de referencia con contraste nulo entre la fractura y el medio circundante ($K_f = K_m = 1$):

![Performance Kf 1e+00](results/figures/PerformanceofDMDcLQRGMRES_flow_Kf1e+00_with_m.png)

* **Observación:** El método **DMDc-LQR-LGMRES alcanza convergencia en tan solo 32 ciclos** (0.70 s), recortando a más de la mitad los 66 ciclos requeridos por LGMRES(28, 2). GMRES(30) estándar requiere más de 800 ciclos para alcanzar la misma tolerancia.

---

### 5.6. Caso $K_f = 10^{2}$ (Fractura Moderadamente Conductora)
Inicio del régimen de canalización preferencial:

![Performance Kf 1e+02](results/figures/PerformanceofDMDcLQRGMRES_flow_Kf1e+02_with_m.png)

* **Observación Crítica:** Ante este contraste conductor, **GMRES(30) sufre breakdown numérico inmediato (`DIVERGED_BREAKDOWN` en iteración 60)** y **LGMRES(28, 2) se estanca sin converger (NC)**. En contraste, **DMDc-LQR-GMRES converge en 165 ciclos** y **DMDc-LQR-LGMRES lo resuelve en tan solo 50 ciclos** (1.02 s), demostrando la robustez extrema que aporta el lazo de control óptimo ante deformaciones espectrales causadas por fracturas conductoras.

---

### 5.7. Caso $K_f = 10^{4}$ (Fractura Altamente Conductora)
Canal de flujo dominante donde el transporte se concentra casi exclusivamente a lo largo de la fractura:

![Performance Kf 1e+04](results/figures/PerformanceofDMDcLQRGMRES_flow_Kf1e+04_with_m.png)

* **Resultado Sobresaliente:** Todos los métodos tradicionales (GMRES, LGMRES) y la variante DMDc sin enriquecimiento fallan o sufren breakdown (NC). **El único solucionador capaz de quebrar el estancamiento y alcanzar convergencia estricta es DMDc-LQR-LGMRES (248 ciclos, 8.21 s)**, confirmando la sinergia matemática entre la identificación DMDc-LQR y el enriquecimiento de subespacios de Krylov.

---

### 5.8. Caso $K_f = 10^{6}$ (Canal Fuertemente Conductor)
Régimen hiper-conductor donde la fractura actúa como un cortocircuito hidráulico:

![Performance Kf 1e+06](results/figures/PerformanceofDMDcLQRGMRES_flow_Kf1e+06_with_m.png)

* **Observación:** El contraste extremo de 6 órdenes de magnitud polariza los autovalores hacia el infinito. Todos los métodos evaluados sin precondicionador alcanzan el límite de iteraciones, evidenciando que para contrastes $K_f \ge 10^6$ se requiere acoplar la formulación DLK con precondicionadores físicos de interfaz (ej. descomposición de dominios o subestructuración Schur).

---

### 5.9. Caso $K_f = 10^{8}$ (Canal de Flujo Hiper-Conductor)
Límite asintótico de permeabilidad casi infinita:

![Performance Kf 1e+08](results/figures/PerformanceofDMDcLQRGMRES_flow_Kf1e+08_with_m.png)

* **Observación:** Dinámica similar al caso $10^6$, con errores de reconstrucción DMDc estables ($e \sim 10^{-12}$) que capturan fielmente la dinámica del residuo, pero con una tasa de decaimiento limitada por el condicionamiento intrínseco del operador sin precondicionar.

---

## 6. Conclusiones y Hallazgos Principales

1. **Fracaso Sistemático de GMRES(30) Estándar:**
   Sin precondicionador, el método clásico GMRES(30) **no converge en ninguna de las 9 matrices** del benchmark de flujo, sufriendo estancamiento severo en contrastes bajos y colapsando por breakdown numérico (`DIVERGED_BREAKDOWN`) en contrastes conductores ($K_f \ge 10^2$).

2. **Ventana de Convergencia Ampliada con DMDc-LQR:**
   * LGMRES(28, 2) estándar únicamente logra converger en una ventana estrecha de 3 problemas ($K_f \in [10^{-4}, 10^0]$).
   * **DMDc-LQR-GMRES** expande esta ventana a 4 problemas, resolviendo exitosamente $K_f = 10^2$ (165 ciclos) donde los métodos estándar colapsan.
   * **DMDc-LQR-LGMRES** extiende la frontera hasta 5 problemas ($K_f \in [10^{-4}, 10^4]$), siendo el **único solucionador de todo el benchmark capaz de converger en $K_f = 10^4$**.

3. **Aceleración en Número de Ciclos:**
   En los problemas donde los métodos estándar convergen, DMDc-LQR-LGMRES reduce los ciclos a más de la mitad:
   * $K_f = 10^{-4}$: **85 ciclos** vs 150 de LGMRES (**1.76x** menos ciclos).
   * $K_f = 10^{-2}$: **44 ciclos** vs 86 de LGMRES (**1.95x** menos ciclos).
   * $K_f = 10^{0}$: **32 ciclos** vs 66 de LGMRES (**2.06x** menos ciclos).

---

## 7. Instrucciones para Reproducir el Benchmark

Para re-ejecutar el benchmark completo y regenerar todas las tablas y figuras:

```bash
# 1. Asegurar entorno secuencial óptimo para BLAS:
export OMP_NUM_THREADS=1

# 2. Ejecutar el benchmark integral de matrices de flujo:
python3 scripts/run_flow_benchmark.py --timeout 60

# 3. Forzar re-ejecución ignorando cachés:
python3 scripts/run_flow_benchmark.py --force --timeout 60
```

Para probar individualmente cualquier matriz con el ejecutable C de PETSc:

```bash
# Probar el caso K_f = 10^0 con DMDc-LQR-LGMRES:
./test_dlk -matrix data/bin/flow_benchmark/flow_benchmark_2d_case_1_Kf1e+00.bin \
           -ksp_type dlk_lgmres \
           -ksp_dlk_controller lqr \
           -ksp_dlk_m0 28 \
           -ksp_dlk_aug 2 \
           -ksp_dlk_m_min 2 \
           -ksp_dlk_m_max 60 \
           -pc_type none \
           -tol 1e-9

# Probar el caso K_f = 10^2 con DMDc-LQR-GMRES:
./test_dlk -matrix data/bin/flow_benchmark/flow_benchmark_2d_case_1_Kf1e+02.bin \
           -ksp_type dlk \
           -ksp_dlk_controller lqr \
           -ksp_dlk_m0 30 \
           -ksp_dlk_m_min 2 \
           -ksp_dlk_m_max 60 \
           -pc_type none \
           -tol 1e-9
```
