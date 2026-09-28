# KSPDLK: Solver Nativo PETSc de GMRES Reiniciado con Control Adaptativo (DMDc-LQR / PD)

Implementación pura y nativa en **PETSc (C / MPI)** del algoritmo de reinicio adaptativo presentado en la investigación *"Optimal Control for Restarting Minimal Residual Krylov Solvers via Data-Driven Residual Dynamics"*.

El solver está registrado formalmente como un tipo `KSP` de PETSc (`#define KSPDLK "dlk"`), operando de manera idéntica a solvers estándar como `KSPGMRES` o `KSPLGMRES`.

---

## 1. Características Principales

* **Tipo KSP Nativo:** Registrado mediante `KSPRegister("dlk", KSPCreate_DLK)`.
* **Configurable en tiempo de ejecución:** Todas las opciones se controlan con flags estándar de PETSc (`-ksp_type dlk`, `-ksp_dlk_controller lqr`, etc.).
* **Precondicionadores Transparentes:** Compatible con cualquier precondicionador (`-pc_type ilu`, `-pc_type jacobi`, `-pc_type bjacobi`, `-pc_type gamg`).
* **Soporte Paralelo MPI:** Implementado sobre vectores y matrices distribuidas (`mpiaij`, `mpivec`), resolviendo los problemas reducidos de DMDc/DARE de forma ultrarrápida en memoria local.
* **Monitores y Perfilado:** Compatible con `-ksp_monitor`, `-ksp_converged_reason` y `-log_view`.

---

## 2. Opciones de Línea de Comandos (`-help`)

| Opción | Tipo | Valor por Defecto | Descripción |
| :--- | :--- | :--- | :--- |
| `-ksp_type` | string | `gmres` | Seleccionar `dlk` para usar este solver. |
| `-ksp_dlk_m0` | int | `30` | Parámetro de reinicio inicial $m_0$. |
| `-ksp_dlk_m_min` | int | `2` | Cota inferior permitida para el reinicio $m$. |
| `-ksp_dlk_m_max` | int | `60` | Cota superior permitida para el reinicio $m$. |
| `-ksp_dlk_controller` | enum | `lqr` | Tipo de controlador: `lqr`, `pd`, o `none`. |
| `-ksp_dlk_alpha` | real | `1e-4` | Factor de penalización $\alpha$ en la función de costo de LQR. |
| `-ksp_dlk_p` | int | `5` | Número de ciclos por ventana de snapshots $p$. |
| `-ksp_dlk_q` | int | `5` | Rango de truncamiento $q$ en espacio de entradas. |
| `-ksp_dlk_t` | int | `6` | Rango de truncamiento $t$ en espacio de estados reducidos. |

---

## 3. Compilación

```bash
cd petsc_solver
make clean && make
```
Genera el binario ejecutable `test_dlk`.

---

## 4. Ejemplos de Uso

### A. Ejecución básica con monitor iterativo
```bash
./test_dlk -m 50 -n 50 -ksp_type dlk -ksp_dlk_controller lqr -ksp_monitor
```

### B. Uso de controlador PD
```bash
./test_dlk -m 50 -n 50 -ksp_type dlk -ksp_dlk_controller pd -ksp_dlk_m0 20 -ksp_monitor
```

### C. Ejecución paralela con MPI (4 procesos)
```bash
mpirun -np 4 ./test_dlk -m 60 -n 60 -ksp_type dlk -ksp_dlk_controller lqr -pc_type bjacobi
```

### D. Batería de experimentos y curvas de convergencia
```bash
./tests/run_experiments.sh 60 60
```
Genera la tabla comparativa y la imagen en `results/convergence_dlk_vs_gmres.png`.
Visualizable desde Windows en:
`\\wsl.localhost\Ubuntu\home\carlos_benitez\DLK-rest\petsc_solver\results\convergence_dlk_vs_gmres.png`

---

## 5. Estructura del Código

* `include/petsckspdlk.h`: Cabecera pública para incluir en aplicaciones de usuario.
* `include/kspdlkimpl.h`: Definición de la estructura interna `KSP_DLK`.
* `include/dmdc_lqr.h`: Declaración de los motores matemáticos de DMDc, DARE y PD.
* `src/dmdc_lqr.c`: Descomposición SVD basada en la matriz de Gram de snapshots y resolución de Riccati por iteración de punto fijo.
* `src/kspdlk.c`: Métodos virtuales de PETSc (`KSPCreate`, `KSPSetUp`, `KSPSolve`, `KSPDestroy`, `KSPSetFromOptions`, `KSPView`).
* `tests/test_dlk.c`: Programa de prueba con problemas convección-difusión no simétricos.
