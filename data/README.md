# Colección de Matrices Dispersas para DLK-PETSc

Este directorio contiene los sistemas lineales dispersos utilizados para evaluar y comparar los algoritmos de reinicio adaptativo **DLK-GMRES** y **DLK-LGMRES** frente a métodos Krylov estándar en **PETSc**.

Las matrices provienen de la [SuiteSparse Matrix Collection](https://sparse.tamu.edu/) (anteriormente *University of Florida Sparse Matrix Collection*), representando problemas reales de dinámica de fluidos computacional (CFD), análisis estructural, simulación de circuitos integrados y ecuaciones diferenciales parciales (EDP).

---

## 1. Estructura y Contenido del Directorio

Los archivos binarios se organizan en subdirectorios temáticos dentro de `data/bin/`, manteniendo enlaces simbólicos en la raíz para garantizar compatibilidad total con scripts anteriores:

* **`bin/suitesparse/`**: 22 matrices de la colección SuiteSparse utilizadas en los experimentos y benchmarks comparativos del artículo (con alias para compatibilidad como `circuit2.bin`, `orsirr1.bin`, etc.).
* **`bin/flow_benchmark/`**: 9 matrices de problemas dimensionales mixtos de flujo en medios fracturados con contrastes de permeabilidad que varían entre $K_f = 10^{-8}$ y $K_f = 10^8$.
* **`bin/*.bin`**: Enlaces simbólicos relativos hacia los archivos en los subdirectorios anteriores, permitiendo ejecutar `./test_dlk -matrix data/bin/<matriz>.bin` de forma transparente.
* **`convert_mat_to_bin.py`**: Script en Python para serializar la colección SuiteSparse `.mat` a formato `.bin` en `data/bin/suitesparse/`.
* **`../scripts/convert_mixed_dimensional.py`**: Script para serializar las matrices de flujo fracturado a `data/bin/flow_benchmark/`.
* **`mat_collection_test/` y `mat_collection_test_full/`**: Fuentes originales de matrices en formato `.mat` de MATLAB / SuiteSparse.

---

## 2. Formato Binario de PETSc (`.bin`)

Para lograr máxima velocidad de I/O en entornos paralelos y secuenciales de PETSc, las matrices se convierten desde MATLAB/SciPy a formato binario nativo mediante `PetscBinaryIO`.

### Estructura del Archivo Binario
Cada archivo `.bin` empaqueta de forma continua en precisión doble (64 bits IEEE 754):
1. **Encabezado y Matriz $A$ (`MatSeqAIJ`):** Dimensiones $M \times N$, punteros de fila (CSR `indptr`), índices de columna (`indices`) y valores no nulos (`data`).
2. **Encabezado y Vector $b$ (`Vec`):** Dimensión $N$ y valores del término independiente.

### Convención del Vector Lado Derecho ($b$)
* **Matrices con $b$ propio:** Se preserva el vector original provisto en el archivo `.mat` (ej. `sherman5`, `raefsky1`, matrices de flujo fracturado).
* **Matrices sin $b$:** Aquellas matrices que no definían un lado derecho fueron completadas automáticamente con un vector unitario:
  $$b = \begin{bmatrix} 1 \\ 1 \\ \vdots \\ 1 \end{bmatrix} \in \mathbb{R}^N$$
  siguiendo la convención estándar utilizada en los experimentos de referencia en MATLAB.

---

## 3. Inventario de Matrices Disponibles

### 3.1. Colección SuiteSparse (`data/bin/suitesparse/`)

| # | Matriz | Dimensión ($N$) | No Ceros ($nnz$) | Dominio de Aplicación | Origen de $b$ | Ruta Canónica |
| :-: | :--- | :---: | :---: | :--- | :---: | :--- |
| 1 | `add20` | 2,395 | 13,151 | Simulación de Circuitos | Original | `data/bin/suitesparse/add20.bin` |
| 2 | `cavity05` | 1,182 | 32,632 | Dinámica de Fluidos (CFD) | Original | `data/bin/suitesparse/cavity05.bin` |
| 3 | `cavity10` | 2,597 | 76,171 | Dinámica de Fluidos (CFD) | Original | `data/bin/suitesparse/cavity10.bin` |
| 4 | `cdde1` | 961 | 4,681 | Convección-Difusión 2D | Generado ($b = \mathbf{1}$) | `data/bin/suitesparse/cdde1.bin` |
| 5 | `circuit_2` | 4,510 | 21,199 | Simulación de Circuitos | Original | `data/bin/suitesparse/circuit_2.bin` |
| 6 | `ex40` | 7,740 | 456,188 | Mecánica Estructural 3D | Generado ($b = \mathbf{1}$) | `data/bin/suitesparse/ex40.bin` |
| 7 | `fpga_trans_01` | 1,220 | 7,382 | Circuitos FPGA | Original | `data/bin/suitesparse/fpga_trans_01.bin` |
| 8 | `memplus` | 17,758 | 99,147 | Memoria / VLSI | Original | `data/bin/suitesparse/memplus.bin` |
| 9 | `orsirr_1` | 1,030 | 6,858 | Modelado de Reservorios | Generado ($b = \mathbf{1}$) | `data/bin/suitesparse/orsirr_1.bin` |
| 10 | `orsreg_1` | 2,205 | 14,133 | Modelado de Reservorios | Generado ($b = \mathbf{1}$) | `data/bin/suitesparse/orsreg_1.bin` |
| 11 | `pde2961` | 2,961 | 14,585 | EDPs Elípticas 2D | Generado ($b = \mathbf{1}$) | `data/bin/suitesparse/pde2961.bin` |
| 12 | `raefsky1` | 3,242 | 293,409 | Flujo Compresible | Original | `data/bin/suitesparse/raefsky1.bin` |
| 13 | `raefsky2` | 3,242 | 293,551 | Flujo Incompresible | Original | `data/bin/suitesparse/raefsky2.bin` |
| 14 | `rdb2048` | 2,048 | 12,032 | Reacción-Difusión | Generado ($b = \mathbf{1}$) | `data/bin/suitesparse/rdb2048.bin` |
| 15 | `sherman1` | 1,000 | 3,750 | Ingeniería de Petróleo | Original | `data/bin/suitesparse/sherman1.bin` |
| 16 | `sherman4` | 1,104 | 3,786 | Ingeniería de Petróleo | Original | `data/bin/suitesparse/sherman4.bin` |
| 17 | `sherman5` | 3,312 | 20,793 | Ingeniería de Petróleo | Original | `data/bin/suitesparse/sherman5.bin` |
| 18 | `steam2` | 600 | 5,660 | Simulación Geotérmica | Generado ($b = \mathbf{1}$) | `data/bin/suitesparse/steam2.bin` |
| 19 | `wang2` | 2,903 | 19,093 | Dispositivo Semiconductor | Generado ($b = \mathbf{1}$) | `data/bin/suitesparse/wang2.bin` |
| 20 | `wang4` | 26,068 | 177,196 | Dispositivo Semiconductor | Generado ($b = \mathbf{1}$) | `data/bin/suitesparse/wang4.bin` |
| 21 | `watt_1` | 1,856 | 11,360 | Flujo en Medios Porosos | Generado ($b = \mathbf{1}$) | `data/bin/suitesparse/watt_1.bin` |
| 22 | `young3c` | 841 | 3,988 | Aerodinámica / Acústica | Generado ($b = \mathbf{1}$) | `data/bin/suitesparse/young3c.bin` |

*Nota:* También se encuentran disponibles los alias directos `circuit2.bin`, `fpgatrans01.bin`, `orsirr1.bin`, `orsreg1.bin` y `watt1.bin`.

### 3.2. Problemas Dimensionales Mixtos en Medios Fracturados (`data/bin/flow_benchmark/`)

Conjunto de sistemas lineales correspondientes a problemas de flujo en medios porosos 2D con fracturas discretas (Mixed-Dimensional Flow) sujetos a variaciones extremas en la permeabilidad de fractura ($K_f$ desde $10^{-8}$ hasta $10^8$). Para un análisis exhaustivo de los resultados experimentales, tablas comparativas y gráficos de rendimiento multietapa, consultar **[`../README_flow.md`](../README_flow.md)**:

| # | Archivo Binario | Dimensión ($N$) | No Ceros ($nnz$) | Permeabilidad $K_f$ | Régimen Físico | Origen de $b$ |
| :-: | :--- | :---: | :---: | :---: | :--- | :---: |
| 1 | `flow_benchmark_2d_case_1_Kf1e-08.bin` | 11,089 | 76,208 | $10^{-8}$ | Fractura altamente impermeable (barrera) | Original |
| 2 | `flow_benchmark_2d_case_1_Kf1e-06.bin` | 11,089 | 76,208 | $10^{-6}$ | Fractura de baja permeabilidad | Original |
| 3 | `flow_benchmark_2d_case_1_Kf1e-04.bin` | 11,089 | 76,208 | $10^{-4}$ | Fractura moderadamente resistente | Original |
| 4 | `flow_benchmark_2d_case_1_Kf1e-02.bin` | 11,089 | 76,208 | $10^{-2}$ | Transición matriz-fractura | Original |
| 5 | `flow_benchmark_2d_case_1_Kf1e+00.bin` | 11,089 | 76,208 | $10^{0}$ | Permeabilidad homogénea | Original |
| 6 | `flow_benchmark_2d_case_1_Kf1e+02.bin` | 11,089 | 76,208 | $10^{2}$ | Fractura moderadamente conductora | Original |
| 7 | `flow_benchmark_2d_case_1_Kf1e+04.bin` | 11,089 | 76,208 | $10^{4}$ | Fractura conductora | Original |
| 8 | `flow_benchmark_2d_case_1_Kf1e+06.bin` | 11,089 | 76,208 | $10^{6}$ | Fractura altamente conductora | Original |
| 9 | `flow_benchmark_2d_case_1_Kf1e+08.bin` | 11,089 | 76,208 | $10^{8}$ | Canal de flujo casi ideal | Original |

---

## 4. Cómo Regenerar o Agregar Nuevas Matrices

Para volver a generar o actualizar todos los archivos `.bin` desde los `.mat`:

```bash
# 1. Regenerar matrices de SuiteSparse:
python3 data/convert_mat_to_bin.py

# 2. Regenerar matrices de flujo fracturado:
python3 scripts/convert_mixed_dimensional.py
```

Ambos scripts empaquetan en sus respectivos subdirectorios (`suitesparse/` y `flow_benchmark/`) y refrescan automáticamente los enlaces simbólicos en `data/bin/` para total compatibilidad.

---

## 5. Uso Rápido en C / PETSc

Para ejecutar cualquier matriz con el ejecutable `test_dlk`:

```bash
# Resolver sherman5 con DLK-GMRES (LQR adaptativo) usando la ruta agrupada:
./test_dlk -matrix data/bin/suitesparse/sherman5.bin -ksp_type dlk -ksp_dlk_controller lqr -pc_type none -tol 1e-9

# O bien mediante el enlace de retrocompatibilidad:
./test_dlk -matrix data/bin/sherman5.bin -ksp_type dlk -ksp_dlk_controller lqr -pc_type none -tol 1e-9

# Resolver una matriz de flujo fracturado (con -pc_type none dado que la matriz posee bloques interfaz tipo punto de silla):
./test_dlk -matrix data/bin/flow_benchmark/flow_benchmark_2d_case_1_Kf1e+00.bin -ksp_type dlk -ksp_dlk_controller lqr -pc_type none -tol 1e-9
```

