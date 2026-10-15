# DLK-rest: Solvers Krylov con Control Adaptativo (DMDc-LQR / PD) en PETSc

Implementación nativa y de alto rendimiento en **PETSc (C / MPI)** de la familia de solvers de Krylov con reinicio adaptativo basados en dinámica de residuos impulsada por datos, correspondiente a la investigación:
> **"Optimal Control for Restarting Minimal Residual Krylov Solvers via Data-Driven Residual Dynamics"**  
> *Carlos Benitez, Juan C. Cabral, Christian E. Schaerer.*

El paquete implementa tanto la variante **DLK-GMRES** como la variante enriquecida **DLK-LGMRES** como tipos de solvers `KSP` nativos de primera clase en PETSc.

---

## 1. Arquitectura y Solvers Disponibles

Los solvers están formalmente registrados en el catálogo de PETSc y pueden invocarse con `-ksp_type`:

| Tipo KSP (`-ksp_type`) | Algoritmo Base | Subespacio de Enriquecimiento ($\ell$) | Descripción |
| :--- | :--- | :--- | :--- |
| `dlk` | **GMRES** | $\ell = 0$ (Estándar) | GMRES adaptativo con control DMDc-LQR o PD sobre el parámetro de reinicio $m_k$. |
| `dlk_lgmres` | **LGMRES** | $\ell = 2$ (Configurable con `-ksp_dlk_aug`) | Loose GMRES que enriquece el subespacio con aproximaciones de error de ciclos previos más control adaptativo. |

---

## 2. Parámetros de Configuración por Línea de Comandos

Todas las opciones son configurables en tiempo de ejecución sin necesidad de recompilar:

| Opción | Tipo | Valor por Defecto | Descripción |
| :--- | :--- | :--- | :--- |
| `-ksp_type` | string | `dlk` | Seleccionar `dlk` o `dlk_lgmres`. |
| `-ksp_dlk_controller` | enum | `lqr` | Tipo de controlador adaptativo: `lqr`, `pd` o `none`. |
| `-ksp_dlk_m0` | int | `30` | Parámetro de reinicio inicial $m_0$. |
| `-ksp_dlk_m_min` | int | `2` | Cota mínima de reinicio $m$. |
| `-ksp_dlk_m_max` | int | $\lfloor 0.1 \times N \rfloor$ | Cota máxima de reinicio $m$ (escala automáticamente según $N$). |
| `-ksp_dlk_aug` | int | `2` (`dlk_lgmres`) / `0` (`dlk`) | Número de vectores de enriquecimiento $\ell$ para LGMRES. |
| `-ksp_dlk_alpha` | real | `1e-4` | Factor de penalización $\alpha$ en la función de costo de LQR. |
| `-ksp_dlk_p` | int | `5` | Tamaño de ventana de snapshots $p$ para DMDc. |
| `-ksp_dlk_q` | int | `5` | Rango de truncamiento $q$ para entradas en DMDc. |
| `-ksp_dlk_t` | int | `6` | Rango de truncamiento $t$ para estados reducidos en DMDc. |
| `-cycle_file` | string | `""` | Exporta la evolución ciclo a ciclo ($m_k$, residuo, error DMDc) a un archivo CSV. |
| `-history_file` | string | `""` | Exporta la trayectoria del residuo relativo en cada iteración interna a un archivo CSV. |

---

## 3. Compilación y Construcción

Para compilar el proyecto (requiere PETSc con BLAS/LAPACK y MPI instalados):

```bash
make clean && make
```
Genera el binario ejecutable `test_dlk`.

---

## 4. Resultados Numéricos: Matriz Sherman5 (SuiteSparse)

La matriz **`sherman5`** ($N = 3312 \times 3312$, $\text{nnz} = 20\,793$, simulación de yacimientos petrolíferos) representa un caso de prueba desafiante donde los métodos de reinicio fijos sufren estancamiento severo sin precondicionamiento.

### A. Rendimiento sin Precondicionador (`-pc_type none`, $\text{tol} = 10^{-9}$)

| Solver (KSP) | Configuración | Ciclos ($k$) | Iteraciones Totales | Residuo Final $\|b - Ax\|_2$ | Residuo Relativo | Tiempo (s) | Estado |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **GMRES(30)** | Fijo $m = 30$ | >333 | 10 000 | $5.0321 \times 10^{1}$ | $8.11 \times 10^{-1}$ | 1.01 s | **Estancado** |
| **GMRES(110)** | Fijo $m = 110$ | >90 | 10 000 | $8.4280 \times 10^{-5}$ | $1.36 \times 10^{-6}$ | 2.14 s | **No converge** |
| **LGMRES(107, 3)** | Fijo $m = 107, \ell = 3$ | 50 | 5 299 | $6.1691 \times 10^{-8}$ | $9.94 \times 10^{-10}$ | 1.36 s | Convergió |
| **DLK-LQR-GMRES** | Adaptativo ($m_0 = 30, m_{\max} = 331$) | **28** | 7 602 | $6.2051 \times 10^{-8}$ | **$9.995 \times 10^{-10}$** | 4.65 s | **Convergió** |
| **DLK-LQR-LGMRES** | Adaptativo ($\ell = 3, m_0 = 30, m_{\max} = 331$) | **47** | **4 925** | $6.2069 \times 10^{-8}$ | **$9.998 \times 10^{-10}$** | **1.71 s** | **Convergió (Óptimo)** |

### B. Rendimiento con Precondicionador ILU(0)

| Solver (KSP) | Algoritmo / Control | Iteraciones Totales | Residuo Final $\|b - Ax\|_2$ | Tiempo (s) |
| :--- | :--- | :---: | :---: | :---: |
| **DLK-LQR-LGMRES** | LGMRES Adaptativo ($\ell = 2, m_0 = 30$) | **50** | $1.2623 \times 10^{-6}$ | 0.0111 s |
| **DLK-LQR-GMRES** | GMRES Adaptativo ($m_0 = 30$) | **50** | $1.2623 \times 10^{-6}$ | 0.0097 s |
| **LGMRES Nativo** | Fijo ($m = 30, \ell = 2$) | 51 | $1.5987 \times 10^{-6}$ | 0.0088 s |
| **GMRES Nativo** | Fijo ($m = 30$) | 50 | $1.2623 \times 10^{-6}$ | 0.0102 s |
| **BiCGStab** | Baseline Krylov | 28 | $1.2425 \times 10^{-7}$ | 0.0060 s |

---

## 5. Comparación de Métodos de Krylov para Publicación (Estilo SIAM)

El script `plot_comparison_of_methods.py` reproduce la figura estándar de publicación científica que contrasta las trayectorias de convergencia en escala de grises con marcadores diferenciados:

![Comparison of Methods sherman5](results/Comparisonofmethodssherman5.png)

*(Archivos vectoriales y de alta resolución: `results/Comparisonofmethodssherman5.png` y `results/Comparisonofmethodssherman5.eps`)*

---

## 6. Curvas de Rendimiento con Dinámica de $m$ (Estilo SIAM)

El script `plot_performance_with_m.py` genera las figuras en formato SIAM de 3 paneles (curva de residuo, barras de error de aproximación DMDc con gradiente y trayectoria adaptativa de $m$):

### Dinámica de DLK-LQR-LGMRES en Sherman5:
![Performance DLK-LGMRES sherman5](results/PerformanceofDMDcLQRLGMRESmellsherman5_with_m.png)

### Dinámica de DLK-LQR-GMRES en Sherman5:
![Performance DLK-GMRES sherman5](results/PerformanceofDMDcLQRGMRESmsherman5_with_m.png)

---

## 7. Ejemplos de Ejecución

### A. Resolver matriz de SuiteSparse con DLK-LQR-GMRES:
```bash
./test_dlk -matrix data/sherman5.bin -ksp_type dlk -ksp_dlk_m0 30 -ksp_dlk_controller lqr -pc_type none -tol 1e-9 -ksp_max_it 50000 -cycle_file results/cycles_sherman5_dlk.csv
```

### B. Resolver con DLK-LQR-LGMRES ($\ell = 3$):
```bash
./test_dlk -matrix data/sherman5.bin -ksp_type dlk_lgmres -ksp_dlk_m0 30 -ksp_dlk_aug 3 -ksp_dlk_controller lqr -pc_type none -tol 1e-9 -ksp_max_it 50000 -cycle_file results/cycles_sherman5_dlklgmres.csv
```

### C. Generar gráfica de comparación de métodos (PNG y EPS):
```bash
python3 plot_comparison_of_methods.py --output_prefix results/Comparisonofmethodssherman5 --title sherman5 --tol 1e-9
```

### D. Generar gráficas de dinámica adaptativa SIAM en 3 paneles (PNG y EPS):
```bash
python3 plot_performance_with_m.py --cycle_file results/cycles_sherman5_dlklgmres.csv --output_prefix results/PerformanceofDMDcLQRLGMRESmellsherman5_with_m --title "DMDc-LQR-LGMRES - sherman5" --m_max 331
```

### E. Ejecución en Paralelo con MPI (4 procesos):
```bash
mpirun -np 4 ./test_dlk -m 80 -n 80 -ksp_type dlk_lgmres -ksp_dlk_controller lqr -pc_type bjacobi
```

---

## 8. Estructura del Repositorio

* `include/petsckspdlk.h`: Prototipos de la API pública de PETSc para `KSPDLK` y `KSPDLKLGMRES`.
* `include/kspdlkimpl.h`: Estructura interna de datos (`KSP_DLK`), base de Arnoldi, vectores de enriquecimiento y buffer de snapshots.
* `include/dmdc_lqr.h`: Declaración de los algoritmos de DMDc y Riccati CARE continuo con reducción por controlabilidad.
* `src/dmdc_lqr.c`: Identificación DMDc basada en Gram SVD, resolución CARE de Laub (1979) y cálculo del error de reconstrucción de Frobenius.
* `src/kspdlk.c`: Implementación del ciclo de vida KSP de PETSc (`Create`, `SetUp`, `Solve`, `Destroy`, `SetFromOptions`, `View`).
* `tests/test_dlk.c`: Driver CLI con soporte de carga de matrices binarias y exportación de historiales.
* `plot_comparison_of_methods.py`: Generador de figuras de comparación de métodos estilo SIAM (PNG y EPS vectorial).
* `plot_performance_with_m.py`: Generador de figuras de rendimiento en 3 paneles (PNG y EPS vectorial).
* `data/sherman5.bin`: Matriz real SuiteSparse convertida a binario de PETSc.
* `results/`: Tablas de resultados, archivos CSV de ciclos y figuras vectoriales/ráster.
