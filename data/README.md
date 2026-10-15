# Colección de Matrices Dispersas para DLK-PETSc

Este directorio contiene los sistemas lineales dispersos utilizados para evaluar y comparar los algoritmos de reinicio adaptativo **DLK-GMRES** y **DLK-LGMRES** frente a métodos Krylov estándar en **PETSc**.

Las matrices provienen de la [SuiteSparse Matrix Collection](https://sparse.tamu.edu/) (anteriormente *University of Florida Sparse Matrix Collection*), representando problemas reales de dinámica de fluidos computacional (CFD), análisis estructural, simulación de circuitos integrados y ecuaciones diferenciales parciales (EDP).

---

## 1. Organización del Directorio

```text
data/
├── README.md                     <- Este archivo de documentación
├── convert_mat_to_bin.py         <- Script para convertir archivos .mat a .bin de PETSc
├── bin/                          <- Matrices serializadas en formato binario nativo PETSc
│   ├── add20.bin
│   ├── cavity05.bin
│   ├── ...
│   └── young3c.bin
├── mat_collection_test/          <- Conjunto inicial de pruebas (.mat)
└── mat_collection_test_full/     <- Colección completa de matrices de evaluación (.mat)
```

---

## 2. Formato Binario de PETSc (`.bin`)

Para lograr máxima velocidad de I/O en entornos paralelos y secuenciales de PETSc, las matrices se convierten desde MATLAB/SciPy a formato binario nativo mediante `PetscBinaryIO`.

### Estructura del Archivo Binario
Cada archivo `.bin` empaqueta de forma continua en precisión doble (64 bits IEEE 754):
1. **Encabezado y Matriz $A$ (`MatSeqAIJ`):** Dimensiones $M \times N$, punteros de fila (CSR `indptr`), índices de columna (`indices`) y valores no nulos (`data`).
2. **Encabezado y Vector $b$ (`Vec`):** Dimensión $N$ y valores del término independiente.

### Convención del Vector Lado Derecho ($b$)
* **Matrices con $b$ propio:** Se preserva el vector original provisto en el archivo `.mat` (ej. `sherman5`, `raefsky1`).
* **Matrices sin $b$:** Aquellas matrices que no definían un lado derecho fueron completadas automáticamente con un vector unitario:
  $$b = \begin{bmatrix} 1 \\ 1 \\ \vdots \\ 1 \end{bmatrix} \in \mathbb{R}^N$$
  siguiendo la convención estándar utilizada en los experimentos de referencia en MATLAB.

---

## 3. Inventario de Matrices Disponibles (`data/bin/`)

| # | Matriz | Dimensión ($N$) | No Ceros ($nnz$) | Dominio de Aplicación | Origen de $b$ | Archivo Binario |
| :-: | :--- | :---: | :---: | :--- | :---: | :--- |
| 1 | `add20` | 2,395 | 13,151 | Simulación de Circuitos | Original | `data/bin/add20.bin` |
| 2 | `cavity05` | 1,182 | 32,632 | Dinámica de Fluidos (CFD) | Original | `data/bin/cavity05.bin` |
| 3 | `cavity10` | 2,597 | 76,171 | Dinámica de Fluidos (CFD) | Original | `data/bin/cavity10.bin` |
| 4 | `cdde1` | 961 | 4,681 | Convección-Difusión 2D | Generado ($b = \mathbf{1}$) | `data/bin/cdde1.bin` |
| 5 | `circuit_2` | 4,510 | 21,199 | Simulación de Circuitos | Original | `data/bin/circuit_2.bin` |
| 6 | `ex40` | 7,740 | 456,188 | Mecánica Estructural 3D | Generado ($b = \mathbf{1}$) | `data/bin/ex40.bin` |
| 7 | `fpga_trans_01` | 1,220 | 7,382 | Circuitos FPGA | Original | `data/bin/fpga_trans_01.bin` |
| 8 | `memplus` | 17,758 | 99,147 | Memoria / VLSI | Original | `data/bin/memplus.bin` |
| 9 | `orsirr_1` | 1,030 | 6,858 | Modelado de Reservorios | Generado ($b = \mathbf{1}$) | `data/bin/orsirr_1.bin` |
| 10 | `orsreg_1` | 2,205 | 14,133 | Modelado de Reservorios | Generado ($b = \mathbf{1}$) | `data/bin/orsreg_1.bin` |
| 11 | `pde2961` | 2,961 | 14,585 | EDPs Elípticas 2D | Generado ($b = \mathbf{1}$) | `data/bin/pde2961.bin` |
| 12 | `raefsky1` | 3,242 | 293,409 | Flujo Compresible | Original | `data/bin/raefsky1.bin` |
| 13 | `raefsky2` | 3,242 | 293,551 | Flujo Incompresible | Original | `data/bin/raefsky2.bin` |
| 14 | `rdb2048` | 2,048 | 12,032 | Reacción-Difusión | Generado ($b = \mathbf{1}$) | `data/bin/rdb2048.bin` |
| 15 | `sherman1` | 1,000 | 3,750 | Ingeniería de Petróleo | Original | `data/bin/sherman1.bin` |
| 16 | `sherman4` | 1,104 | 3,786 | Ingeniería de Petróleo | Original | `data/bin/sherman4.bin` |
| 17 | `sherman5` | 3,312 | 20,793 | Ingeniería de Petróleo | Original | `data/bin/sherman5.bin` |
| 18 | `steam2` | 600 | 5,660 | Simulación Geotérmica | Generado ($b = \mathbf{1}$) | `data/bin/steam2.bin` |
| 19 | `wang2` | 2,903 | 19,093 | Dispositivo Semiconductor | Generado ($b = \mathbf{1}$) | `data/bin/wang2.bin` |
| 20 | `wang4` | 26,068 | 177,196 | Dispositivo Semiconductor | Generado ($b = \mathbf{1}$) | `data/bin/wang4.bin` |
| 21 | `watt_1` | 1,856 | 11,360 | Flujo en Medios Porosos | Generado ($b = \mathbf{1}$) | `data/bin/watt_1.bin` |
| 22 | `young3c` | 841 | 3,988 | Aerodinámica / Acústica | Generado ($b = \mathbf{1}$) | `data/bin/young3c.bin` |

---

## 4. Cómo Regenerar o Agregar Nuevas Matrices

Para volver a generar o actualizar todos los archivos `.bin` desde los `.mat`:

```bash
# Ejecutar desde la raíz del repositorio:
python3 data/convert_mat_to_bin.py
```

El script buscará automáticamente cualquier archivo `.mat` dentro de `data/mat_collection_test_full/` y `data/mat_collection_test/`, detectará si contiene vector $b$ o lo generará si hace falta, y exportará el archivo compilado a `data/bin/<nombre>.bin`.

---

## 5. Uso Rápido en C / PETSc

Para ejecutar cualquier matriz con el ejecutable `test_dlk`:

```bash
# Resolver sherman5 con DLK-GMRES (LQR adaptativo):
./test_dlk -matrix data/bin/sherman5.bin -ksp_type dlk -ksp_dlk_controller lqr -pc_type none -tol 1e-9

# Resolver raefsky1 con DLK-LGMRES (LQR adaptativo con 2 vectores de enriquecimiento):
./test_dlk -matrix data/bin/raefsky1.bin -ksp_type dlk_lgmres -ksp_dlk_aug 2 -ksp_dlk_controller lqr -pc_type ilu -tol 1e-9
```
