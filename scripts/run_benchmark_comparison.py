#!/usr/bin/env python3
"""
scripts/run_benchmark_comparison.py

Ejecuta el benchmark de solvers Krylov en PETSc sobre las matrices de la colección SuiteSparse
(data/bin/) y exporta exactamente los archivos CSV que utiliza el artículo en LaTeX:
1. results/tables/matrix_info.csv
2. results/tables/cycles_experiment.csv
3. results/tables/avg_relativa_time_experiment.csv

Los 4 métodos evaluados son:
- GMRES (Cycles_GMRES, avg_time_GMRES)
- LGMRES (Cycles_LGMRES, avg_time_LGMRES)
- GMRES_DMDc_LQR (Cycles_GMRES_DMDc_LQR, avg_time_GMRES_DMDc_LQR)
- LGMRES_DMDc_LQR (Cycles_LGMRES_DMDc_LQR, avg_time_LGMRES_DMDc_LQR)
"""

import os
import sys
import subprocess
import csv
import re
import argparse

def run_solver(bin_file, solver_type, m0, m_max, tol=1e-8, max_it=10000, pc_type="none", aug_l=2):
    cmd = ["./test_dlk", "-matrix", bin_file, "-tol", str(tol), "-ksp_max_it", str(max_it), "-pc_type", pc_type]

    if solver_type == "GMRES":
        cmd += ["-ksp_type", "gmres", "-ksp_gmres_restart", str(m0)]
    elif solver_type == "LGMRES":
        cmd += ["-ksp_type", "lgmres", "-ksp_lgmres_restart", str(m0), "-ksp_lgmres_augment", str(aug_l)]
    elif solver_type == "GMRES_DMDc_LQR":
        cmd += [
            "-ksp_type", "dlk",
            "-ksp_dlk_controller", "lqr",
            "-ksp_dlk_m0", str(m0),
            "-ksp_dlk_m_min", str(max(2, m0 // 2)),
            "-ksp_dlk_m_max", str(m_max),
            "-ksp_dlk_alpha", "0.001"
        ]
    elif solver_type == "LGMRES_DMDc_LQR":
        cmd += [
            "-ksp_type", "dlk_lgmres",
            "-ksp_dlk_controller", "lqr",
            "-ksp_dlk_m0", str(m0),
            "-ksp_dlk_m_min", str(max(2, m0 // 2)),
            "-ksp_dlk_m_max", str(m_max),
            "-ksp_dlk_aug", str(aug_l),
            "-ksp_dlk_alpha", "0.001"
        ]

    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=120)
        output = res.stdout
    except subprocess.TimeoutExpired:
        return "NC", "NC"
    except Exception:
        return "NC", "NC"

    its_m = re.search(r"Iteraciones Totales\s*:\s*(\d+)", output)
    cyc_m = re.search(r"Ciclos Totales\s*:\s*(\d+)", output)
    res_m = re.search(r"Residuo \|\|b - Ax\|\|\s*:\s*([0-9.eE+-]+)", output)
    time_m = re.search(r"Tiempo de Ejecución\s*:\s*([0-9.]+)", output)

    its = int(its_m.group(1)) if its_m else -1
    cycles = int(cyc_m.group(1)) if cyc_m else -1
    residual = float(res_m.group(1)) if res_m else 1.0
    elapsed = float(time_m.group(1)) if time_m else 0.0

    # Criterio de convergencia: residuo alcanzado y no superó límite de iteraciones
    converged = (its < max_it) and (residual <= tol * 50)
    if not converged or cycles <= 0:
        return "NC", "NC"

    return cycles, round(elapsed, 5)

def main():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    os.chdir(repo_root)

    matrix_info_path = "results/tables/matrix_info.csv"
    cycles_csv_path = "results/tables/cycles_experiment.csv"
    time_csv_path = "results/tables/avg_relativa_time_experiment.csv"

    parser = argparse.ArgumentParser(description="Ejecutar benchmark comparativo para artículo LaTeX")
    parser.add_argument('--problems', type=str, default="", help="Lista de números de problema separados por comas (ej. 1,4,15,16,17,18,22)")
    parser.add_argument('--all', action='store_true', help="Ejecutar todos los 22 problemas")
    args = parser.parse_args()

    # Cargar información de matrices
    problems = []
    with open(matrix_info_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            problems.append({
                "Problem": int(row["Problem"]),
                "Matrix": row["Matrix"],
                "n": int(row["n"]),
                "nnz": int(row["nnz"]),
                "cond_A": row["cond_A"]
            })

    # Cargar resultados previos si existen
    cycles_data = {}
    time_data = {}
    if os.path.exists(cycles_csv_path):
        with open(cycles_csv_path, 'r') as f:
            for row in csv.DictReader(f):
                cycles_data[int(row["Problem"])] = row
    if os.path.exists(time_csv_path):
        with open(time_csv_path, 'r') as f:
            for row in csv.DictReader(f):
                time_data[int(row["Problem"])] = row

    # Determinar qué problemas ejecutar
    selected_probs = set()
    if args.all:
        selected_probs = set(p["Problem"] for p in problems)
    elif args.problems:
        selected_probs = set(int(x.strip()) for x in args.problems.split(","))
    else:
        # Por defecto ejecutar conjunto representativo si no se especifica
        selected_probs = {4, 15, 16, 17, 18, 22}

    print("=" * 110)
    print(" BENCHMARK SUITESPARSE PETSc (Formato de Tablas LaTeX)")
    print("=" * 110)
    print(f"{'#':<3} | {'Matriz':<14} | {'N':<6} | {'Cycles GMRES':<13} | {'Cycles LGMRES':<14} | {'Cycles DLK':<11} | {'Cycles DLK-L':<13} | {'Tiempo DLK (s)'}")
    print("-" * 110)

    for p in problems:
        pid = p["Problem"]
        mat_name = p["Matrix"]
        n = p["n"]

        if pid not in selected_probs and pid in cycles_data and pid in time_data:
            # Preservar datos existentes
            c_row = cycles_data[pid]
            t_row = time_data[pid]
            print(f"{pid:<3} | {mat_name:<14} | {n:<6} | {c_row['Cycles_GMRES']:<13} | {c_row['Cycles_LGMRES']:<14} | {c_row['Cycles_GMRES_DMDc_LQR']:<11} | {c_row['Cycles_LGMRES_DMDc_LQR']:<13} | {t_row['avg_time_GMRES_DMDc_LQR']} (previo)")
            continue

        bin_path = f"data/bin/{mat_name}.bin"
        if not os.path.exists(bin_path):
            print(f"[Aviso] No existe {bin_path}, omitiendo problema {pid}.")
            continue

        m0 = min(30, max(10, int(round(0.10 * n))))
        m_max = min(n, max(30, int(round(0.10 * n))))
        tol = 1e-8

        # sherman5 requiere m_max = 331 y tol = 1e-9 sin precondicionador
        if mat_name == "sherman5":
            m0 = 30
            m_max = 331
            tol = 1e-9

        # Ejecutar los 4 métodos
        c_gmres, t_gmres = run_solver(bin_path, "GMRES", m0, m_max, tol=tol)
        c_lgmres, t_lgmres = run_solver(bin_path, "LGMRES", m0, m_max, tol=tol, aug_l=2)
        c_dlk, t_dlk = run_solver(bin_path, "GMRES_DMDc_LQR", m0, m_max, tol=tol)
        c_dlkl, t_dlkl = run_solver(bin_path, "LGMRES_DMDc_LQR", m0, m_max, tol=tol, aug_l=2)

        cycles_data[pid] = {
            "Problem": pid,
            "Matrix": mat_name,
            "Cycles_GMRES": c_gmres,
            "Cycles_LGMRES": c_lgmres,
            "Cycles_GMRES_DMDc_LQR": c_dlk,
            "Cycles_LGMRES_DMDc_LQR": c_dlkl
        }

        time_data[pid] = {
            "Problem": pid,
            "Matrix": mat_name,
            "avg_time_GMRES": t_gmres,
            "avg_time_LGMRES": t_lgmres,
            "avg_time_GMRES_DMDc_LQR": t_dlk,
            "avg_time_LGMRES_DMDc_LQR": t_dlkl
        }

        print(f"{pid:<3} | {mat_name:<14} | {n:<6} | {str(c_gmres):<13} | {str(c_lgmres):<14} | {str(c_dlk):<11} | {str(c_dlkl):<13} | {str(t_dlk)}")

    # Guardar cycles_experiment.csv
    with open(cycles_csv_path, 'w', newline='') as f:
        fieldnames = ["Problem", "Matrix", "Cycles_GMRES", "Cycles_LGMRES", "Cycles_GMRES_DMDc_LQR", "Cycles_LGMRES_DMDc_LQR"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for p in problems:
            pid = p["Problem"]
            if pid in cycles_data:
                writer.writerow(cycles_data[pid])

    # Guardar avg_relativa_time_experiment.csv
    with open(time_csv_path, 'w', newline='') as f:
        fieldnames = ["Problem", "Matrix", "avg_time_GMRES", "avg_time_LGMRES", "avg_time_GMRES_DMDc_LQR", "avg_time_LGMRES_DMDc_LQR"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for p in problems:
            pid = p["Problem"]
            if pid in time_data:
                writer.writerow(time_data[pid])

    print("-" * 110)
    print(f"[✓] Tablas actualizadas exitosamente:")
    print(f"    - {cycles_csv_path}")
    print(f"    - {time_csv_path}")
    print("=" * 110)

if __name__ == "__main__":
    main()
