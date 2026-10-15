# Tabla Comparativa de Rendimiento (PETSc Benchmark)

Resultados obtenidos ejecutando solvers Krylov en sistemas lineales SuiteSparse:

| Matriz | Dimensión ($N$) | Solver | Iteraciones | Tiempo (s) | Residuo $\|b-Ax\|$ | Estado |
| :--- | :---: | :--- | :---: | :---: | :---: | :--- |
| `sherman5` | 3312 | **GMRES(m)** | 41 | 0.0085 | 1.18e-06 | No convergió |
| `sherman5` | 3312 | **LGMRES(m,2)** | 47 | 0.0072 | 1.01e-06 | No convergió |
| `sherman5` | 3312 | **DLK-GMRES** | 41 | 0.0344 | 1.18e-06 | No convergió |
| `sherman5` | 3312 | **DLK-LGMRES** | 41 | 0.0158 | 1.18e-06 | No convergió |
| `sherman5` | 3312 | **BiCGStab** | 24 | 0.0049 | 3.20e-07 | Convergió |
| `sherman1` | 1000 | **GMRES(m)** | 61 | 0.0030 | 2.82e-08 | Convergió |
| `sherman1` | 1000 | **LGMRES(m,2)** | 60 | 0.0029 | 3.13e-08 | Convergió |
| `sherman1` | 1000 | **DLK-GMRES** | 61 | 0.0032 | 2.82e-08 | Convergió |
| `sherman1` | 1000 | **DLK-LGMRES** | 63 | 0.0039 | 2.67e-08 | Convergió |
| `sherman1` | 1000 | **BiCGStab** | 34 | 0.0019 | 1.41e-08 | Convergió |
| `sherman4` | 1104 | **GMRES(m)** | 61 | 0.0031 | 9.13e-07 | Convergió |
| `sherman4` | 1104 | **LGMRES(m,2)** | 48 | 0.0026 | 7.71e-07 | Convergió |
| `sherman4` | 1104 | **DLK-GMRES** | 61 | 0.0035 | 9.13e-07 | Convergió |
| `sherman4` | 1104 | **DLK-LGMRES** | 61 | 0.0034 | 9.06e-07 | Convergió |
| `sherman4` | 1104 | **BiCGStab** | 28 | 0.0015 | 1.27e-07 | Convergió |
| `watt_1` | 1856 | **GMRES(m)** | 68 | 0.0061 | 9.40e-07 | Convergió |
| `watt_1` | 1856 | **LGMRES(m,2)** | 60 | 0.0066 | 8.48e-07 | Convergió |
| `watt_1` | 1856 | **DLK-GMRES** | 68 | 0.0072 | 9.40e-07 | Convergió |
| `watt_1` | 1856 | **DLK-LGMRES** | 71 | 0.0083 | 9.07e-07 | Convergió |
| `watt_1` | 1856 | **BiCGStab** | 31 | 0.0050 | 4.63e-07 | Convergió |
| `orsreg_1` | 2205 | **GMRES(m)** | 54 | 0.0058 | 8.83e-06 | No convergió |
| `orsreg_1` | 2205 | **LGMRES(m,2)** | 54 | 0.0066 | 2.76e-06 | No convergió |
| `orsreg_1` | 2205 | **DLK-GMRES** | 54 | 0.0087 | 8.83e-06 | No convergió |
| `orsreg_1` | 2205 | **DLK-LGMRES** | 55 | 0.0091 | 8.83e-06 | No convergió |
| `orsreg_1` | 2205 | **BiCGStab** | 35 | 0.0050 | 2.25e-06 | No convergió |
| `cdde1` | 961 | **GMRES(m)** | 75 | 0.0036 | 5.03e-07 | Convergió |
| `cdde1` | 961 | **LGMRES(m,2)** | 60 | 0.0031 | 4.17e-07 | Convergió |
| `cdde1` | 961 | **DLK-GMRES** | 75 | 0.0043 | 5.03e-07 | Convergió |
| `cdde1` | 961 | **DLK-LGMRES** | 80 | 0.0037 | 4.53e-07 | Convergió |
| `cdde1` | 961 | **BiCGStab** | 29 | 0.0019 | 1.24e-07 | Convergió |
| `steam2` | 600 | **GMRES(m)** | 2 | 0.0027 | 7.52e-12 | Convergió |
| `steam2` | 600 | **LGMRES(m,2)** | 2 | 0.0007 | 7.52e-12 | Convergió |
| `steam2` | 600 | **DLK-GMRES** | 2 | 0.0012 | 7.28e-12 | Convergió |
| `steam2` | 600 | **DLK-LGMRES** | 2 | 0.0011 | 7.28e-12 | Convergió |
| `steam2` | 600 | **BiCGStab** | 1 | 0.0005 | 6.13e-12 | Convergió |
| `cavity05` | 1182 | **GMRES(m)** | 10000 | 0.6803 | 1.50e-16 | No convergió |
| `cavity05` | 1182 | **LGMRES(m,2)** | 1796 | 0.0984 | 1.04e-19 | Convergió |
| `cavity05` | 1182 | **DLK-GMRES** | 297 | 0.0202 | 4.51e-15 | Convergió |
| `cavity05` | 1182 | **DLK-LGMRES** | 264 | 0.0203 | 5.62e-15 | Convergió |
| `cavity05` | 1182 | **BiCGStab** | 695 | 0.0562 | 3.31e-20 | Convergió |
