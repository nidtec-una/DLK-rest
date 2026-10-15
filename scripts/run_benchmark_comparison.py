#!/usr/bin/env python3
"""
scripts/run_benchmark_comparison.py

Ejecuta el conjunto de pruebas comparativo de solvers Krylov en PETSc sobre una selección
de matrices representativas de data/bin/:
- GMRES(30)
- LGMRES(30, 2)
- DLK-GMRES (Adaptativo DMDc-LQR)
- DLK-LGMRES (Adaptativo DMDc-LQR + Augmentation l=2)
- BiCGStab

Genera:
1. results/tables/benchmark_comparison.csv (para posterior generación de gráficos de barras)
2. results/tables/benchmark_comparison.md  (tabla formateada en Markdown)
"""

import os
import sys
import subprocess
import csv
import re

BENCHMARK_MATRICES = [
    {"name": "sherman5",      "m0": 30, "m_max": 331, "tol": 1e-8, "pc_type": "ilu"},
    {"name": "sherman1",      "m0": 20, "m_max": 100, "tol": 1e-8, "pc_type": "ilu"},
    {"name": "sherman4",      "m0": 20, "m_max": 110, "tol": 1e-8, "pc_type": "ilu"},
    {"name": "watt_1",        "m0": 20, "m_max": 185, "tol": 1e-8, "pc_type": "ilu"},
    {"name": "orsreg_1",      "m0": 20, "m_max": 220, "tol": 1e-8, "pc_type": "ilu"},
    {"name": "cdde1",         "m0": 20, "m_max": 96,  "tol": 1e-8, "pc_type": "ilu"},
    {"name": "steam2",        "m0": 15, "m_max": 60,  "tol": 1e-8, "pc_type": "ilu"},
    {"name": "cavity05",      "m0": 20, "m_max": 118, "tol": 1e-8, "pc_type": "none"},
]

SOLVER_CONFIGS = [
    {
        "id": "GMRES",
        "label": "GMRES(m)",
        "cmd_args": lambda m0, m_max: ["-ksp_type", "gmres", "-ksp_gmres_restart", str(m0)]
    },
    {
        "id": "LGMRES",
        "label": "LGMRES(m,2)",
        "cmd_args": lambda m0, m_max: ["-ksp_type", "lgmres", "-ksp_lgmres_restart", str(m0), "-ksp_lgmres_augment", "2"]
    },
    {
        "id": "DLK-GMRES",
        "label": "DLK-GMRES",
        "cmd_args": lambda m0, m_max: [
            "-ksp_type", "dlk",
            "-ksp_dlk_controller", "lqr",
            "-ksp_dlk_m0", str(m0),
            "-ksp_dlk_m_min", str(m0 // 2 if m0 > 10 else 5),
            "-ksp_dlk_m_max", str(m_max),
            "-ksp_dlk_alpha", "0.001"
        ]
    },
    {
        "id": "DLK-LGMRES",
        "label": "DLK-LGMRES",
        "cmd_args": lambda m0, m_max: [
            "-ksp_type", "dlk_lgmres",
            "-ksp_dlk_controller", "lqr",
            "-ksp_dlk_m0", str(m0),
            "-ksp_dlk_m_min", str(m0 // 2 if m0 > 10 else 5),
            "-ksp_dlk_m_max", str(m_max),
            "-ksp_dlk_aug", "2",
            "-ksp_dlk_alpha", "0.001"
        ]
    },
    {
        "id": "BiCGStab",
        "label": "BiCGStab",
        "cmd_args": lambda m0, m_max: ["-ksp_type", "bcgs"]
    }
]

def run_test_dlk(bin_file, extra_args, tol=1e-8, max_it=10000, pc_type="ilu"):
    cmd = ["./test_dlk", "-matrix", bin_file, "-tol", str(tol), "-ksp_max_it", str(max_it), "-pc_type", pc_type] + extra_args
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=120)
        output = res.stdout
    except subprocess.TimeoutExpired:
        return None, None, None, None, "Timeout (>120s)"
    except Exception as e:
        return None, None, None, None, f"Error: {e}"

    dim_m = re.search(r"Dimensión del Sistema\s*:\s*N\s*=\s*(\d+)", output)
    its_m = re.search(r"Iteraciones Totales\s*:\s*(\d+)", output)
    res_m = re.search(r"Residuo \|\|b - Ax\|\|\s*:\s*([0-9.eE+-]+)", output)
    time_m = re.search(r"Tiempo de Ejecución\s*:\s*([0-9.]+)", output)

    dim = int(dim_m.group(1)) if dim_m else 0
    its = int(its_m.group(1)) if its_m else -1
    residual = float(res_m.group(1)) if res_m else 1.0
    elapsed = float(time_m.group(1)) if time_m else 0.0

    converged = (its < max_it) and (residual <= tol * 100)
    status = "Convergió" if converged else "No convergió"

    return dim, its, residual, elapsed, status

def main():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    os.chdir(repo_root)

    out_csv = "results/tables/benchmark_comparison.csv"
    out_md = "results/tables/benchmark_comparison.md"
    os.makedirs("results/tables", exist_ok=True)

    print("=" * 100)
    print(" BENCHMARK COMPARATIVO DE SOLVERS KRYLOV (DLK vs NATIVOS PETSc)")
    print("=" * 100)
    print(f"{'Matriz':<12} | {'Solver':<12} | {'N':<7} | {'Iter':<7} | {'Tiempo (s)':<11} | {'Residuo':<14} | {'Estado'}")
    print("-" * 100)

    rows = []

    for mat_info in BENCHMARK_MATRICES:
        mat_name = mat_info["name"]
        bin_path = f"data/bin/{mat_name}.bin"
        if not os.path.exists(bin_path):
            print(f"[Aviso] Matriz {bin_path} no encontrada, omitiendo.")
            continue

        m0 = mat_info["m0"]
        m_max = mat_info["m_max"]
        tol = mat_info["tol"]
        pc_type = mat_info.get("pc_type", "ilu")

        for sc in SOLVER_CONFIGS:
            args = sc["cmd_args"](m0, m_max)
            dim, its, residual, elapsed, status = run_test_dlk(bin_path, args, tol=tol, pc_type=pc_type)

            row = {
                "Matrix": mat_name,
                "Dimension_N": dim if dim else 0,
                "Solver_ID": sc["id"],
                "Solver_Label": sc["label"],
                "Iterations": its if its is not None else -1,
                "Time_sec": elapsed if elapsed is not None else -1.0,
                "Residual": residual if residual is not None else -1.0,
                "Tolerance": tol,
                "Status": status
            }
            rows.append(row)

            its_str = str(its) if its >= 0 else "N/A"
            time_str = f"{elapsed:.4f}" if elapsed >= 0 else "N/A"
            res_str = f"{residual:.2e}" if residual >= 0 else "N/A"

            print(f"{mat_name:<12} | {sc['id']:<12} | {dim:<7} | {its_str:<7} | {time_str:<11} | {res_str:<14} | {status}")

    # Guardar CSV
    fieldnames = ["Matrix", "Dimension_N", "Solver_ID", "Solver_Label", "Iterations", "Time_sec", "Residual", "Tolerance", "Status"]
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    # Guardar Markdown
    with open(out_md, "w") as f:
        f.write("# Tabla Comparativa de Rendimiento (PETSc Benchmark)\n\n")
        f.write("Resultados obtenidos ejecutando solvers Krylov en sistemas lineales SuiteSparse:\n\n")
        f.write(r"| Matriz | Dimensión ($N$) | Solver | Iteraciones | Tiempo (s) | Residuo $\|b-Ax\|$ | Estado |" + "\n")
        f.write("| :--- | :---: | :--- | :---: | :---: | :---: | :--- |\n")
        for r in rows:
            its_str = str(r['Iterations']) if r['Iterations'] >= 0 else "N/A"
            time_str = f"{r['Time_sec']:.4f}" if r['Time_sec'] >= 0 else "N/A"
            res_str = f"{r['Residual']:.2e}" if r['Residual'] >= 0 else "N/A"
            f.write(f"| `{r['Matrix']}` | {r['Dimension_N']} | **{r['Solver_Label']}** | {its_str} | {time_str} | {res_str} | {r['Status']} |\n")

    print("-" * 100)
    print(f"[✓] Benchmark finalizado exitosamente.")
    print(f"    CSV: {out_csv}")
    print(f"    Markdown: {out_md}")
    print("=" * 100)

if __name__ == "__main__":
    main()
