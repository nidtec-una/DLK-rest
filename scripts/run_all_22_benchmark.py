#!/usr/bin/env python3
"""
scripts/run_all_22_benchmark.py

Ejecuta el benchmark completo de las 22 matrices de SuiteSparse en PETSc
evaluando los 4 métodos exactos del artículo:
1. GMRES(30)
2. LGMRES(28, 2)
3. DMDc-LQR-GMRES(30)
4. DMDc-LQR-LGMRES(28, 2)

Exporta:
- results/tables/matrix_info.csv
- results/tables/cycles_experiment.csv
- results/tables/avg_relativa_time_experiment.csv
"""

import os
import sys
import subprocess
import csv
import re
import time

PROBLEMS_DEF = [
    (1,  "add20",         2395,  13151,  "1.7637e+04"),
    (2,  "cavity05",      1182,  32632,  "9.1829e+05"),
    (3,  "cavity10",      2597,  76171,  "4.4641e+06"),
    (4,  "cdde1",         961,   4681,   "4.0763e+03"),
    (5,  "circuit2",      4510,  21199,  "7.1181e+06"),
    (6,  "ex40",          7740,  456188, "1.2019e+07"),
    (7,  "fpgatrans01",   1220,  7382,   "2.8784e+04"),
    (8,  "memplus",       17758, 99147,  "2.6661e+05"),
    (9,  "orsirr1",       1030,  6858,   "1.6720e+05"),
    (10, "orsreg1",       2205,  14133,  "1.5394e+04"),
    (11, "pde2961",       2961,  14585,  "9.4909e+02"),
    (12, "raefsky1",      3242,  293409, "3.1593e+04"),
    (13, "raefsky2",      3242,  293551, "1.0823e+04"),
    (14, "rdb2048",       2048,  12032,  "2.0918e+03"),
    (15, "sherman1",      1000,  3750,   "2.2575e+04"),
    (16, "sherman4",      1104,  3786,   "7.1579e+03"),
    (17, "sherman5",      3312,  20793,  "3.9023e+05"),
    (18, "steam2",        600,   5660,   "3.5481e+06"),
    (19, "wang2",         2903,  19093,  "3.1844e+04"),
    (20, "wang4",         26068, 177196, "4.9125e+04"),
    (21, "watt1",         1856,  11360,  "5.3847e+09"),
    (22, "young3c",       841,   3988,   "1.1532e+04"),
]

def resolve_bin_file(mat_name):
    aliases = [
        mat_name,
        mat_name.replace('circuit2', 'circuit_2'),
        mat_name.replace('orsirr1', 'orsirr_1'),
        mat_name.replace('orsreg1', 'orsreg_1'),
        mat_name.replace('watt1', 'watt_1'),
        mat_name.replace('fpgatrans01', 'fpga_trans_01')
    ]
    candidates = []
    for a in aliases:
        candidates.append(f"data/bin/suitesparse/{a}.bin")
        candidates.append(f"data/bin/flow_benchmark/{a}.bin")
        candidates.append(f"data/bin/{a}.bin")
    for c in candidates:
        if os.path.exists(c):
            return c
    return f"data/bin/suitesparse/{mat_name}.bin"

def run_petsc_solver(bin_file, solver_type, m0=30, m_max=331, l=2, tol=1e-9, max_cycles=400, timeout=180):
    # En PETSc sin precondicionador (unpreconditioned, como en el benchmark de MATLAB)
    max_its = max(50000, max_cycles * m_max)

    cmd = ["./test_dlk", "-matrix", bin_file, "-tol", str(tol), "-ksp_max_it", str(max_its), "-pc_type", "none"]

    if solver_type == "GMRES":
        cmd += ["-ksp_type", "gmres", "-ksp_gmres_restart", str(m0)]
    elif solver_type == "LGMRES":
        cmd += ["-ksp_type", "lgmres", "-ksp_lgmres_restart", str(max(2, m0 - l)), "-ksp_lgmres_augment", str(l)]
    elif solver_type == "GMRES_DMDc_LQR":
        cmd += [
            "-ksp_type", "dlk",
            "-ksp_dlk_controller", "lqr",
            "-ksp_dlk_m0", str(m0),
            "-ksp_dlk_m_min", "2",
            "-ksp_dlk_m_max", str(m_max),
            "-ksp_dlk_alpha", "0.0001"
        ]
    elif solver_type == "LGMRES_DMDc_LQR":
        cmd += [
            "-ksp_type", "dlk_lgmres",
            "-ksp_dlk_controller", "lqr",
            "-ksp_dlk_m0", str(max(2, m0 - l)),
            "-ksp_dlk_aug", str(l),
            "-ksp_dlk_m_min", "2",
            "-ksp_dlk_m_max", str(m_max),
            "-ksp_dlk_alpha", "0.0001"
        ]

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    t0 = time.time()
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=timeout, env=env)
        output = proc.stdout
        elapsed = time.time() - t0
    except subprocess.TimeoutExpired:
        return "NC", "NC"
    except Exception:
        return "NC", "NC"

    if proc.returncode != 0:
        return "NC", "NC"

    cyc_m = re.search(r"Ciclos Totales\s*:\s*(\d+)", output)
    its_m = re.search(r"Iteraciones Totales\s*:\s*(\d+)", output)
    res_m = re.search(r"Residuo \|\|b - Ax\|\|\s*:\s*([0-9.eE+-]+)", output)
    time_m = re.search(r"Tiempo de Ejecución\s*:\s*([0-9.]+)", output)

    cycles = int(cyc_m.group(1)) if cyc_m else -1
    its = int(its_m.group(1)) if its_m else -1
    residual = float(res_m.group(1)) if res_m else 1.0
    report_time = float(time_m.group(1)) if time_m else elapsed

    # Criterio: converged if residual achieved within max_cycles
    if cycles <= 0 or cycles >= max_cycles or its >= max_its or residual > 1e-4:
        return "NC", "NC"

    return cycles, round(report_time, 5)

def main():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    os.chdir(repo_root)

    out_cycles = "results/tables/cycles_experiment.csv"
    out_time = "results/tables/avg_relativa_time_experiment.csv"
    out_matrix = "results/tables/matrix_info.csv"
    os.makedirs("results/tables", exist_ok=True)

    # 1. Escribir matrix_info.csv
    with open(out_matrix, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Problem", "Matrix", "n", "nnz", "cond_A"])
        for pid, name, n, nnz, cond in PROBLEMS_DEF:
            writer.writerow([pid, name, n, nnz, cond])

    print("=" * 105)
    print(" BENCHMARK COMPLETO DE 22 MATRICES EN PETSc (4 MÉTODOS KRYLOV)")
    print("=" * 105)
    print(f"{'#':<3} | {'Matriz':<14} | {'N':<6} | {'Cycles GMRES':<13} | {'Cycles LGMRES':<14} | {'Cycles DLK':<11} | {'Cycles DLK-L':<13} | {'Tiempo DLK (s)'}")
    print("-" * 105)

    cycles_rows = []
    time_rows = []

    m0_default = 30
    l_default = 2
    tol = 1e-9
    max_cycles = 400

    for pid, name, n, nnz, cond in PROBLEMS_DEF:
        bin_path = resolve_bin_file(name)
        if not os.path.exists(bin_path):
            print(f"[Error] No se encontró {bin_path}")
            cycles_rows.append({"Problem": pid, "Matrix": name, "Cycles_GMRES": "NC", "Cycles_LGMRES": "NC", "Cycles_GMRES_DMDc_LQR": "NC"})
            time_rows.append({"Problem": pid, "Matrix": name, "avg_time_GMRES": "NC", "avg_time_LGMRES": "NC", "avg_time_GMRES_DMDc_LQR": "NC"})
            continue

        m_max = min(n, max(30, int(round(0.10 * n))))
        if name == "sherman5":
            m_max = 331

        # 1. GMRES(30)
        cg, tg = run_petsc_solver(bin_path, "GMRES", m0=m0_default, m_max=m_max, tol=tol, max_cycles=max_cycles, timeout=180)
        # 2. LGMRES(28, 2)
        cl, tl = run_petsc_solver(bin_path, "LGMRES", m0=m0_default, m_max=m_max, l=l_default, tol=tol, max_cycles=max_cycles, timeout=180)
        # 3. DLK-GMRES(30)
        c_dlk, t_dlk = run_petsc_solver(bin_path, "GMRES_DMDc_LQR", m0=m0_default, m_max=m_max, tol=tol, max_cycles=max_cycles, timeout=180)

        cycles_rows.append({
            "Problem": pid,
            "Matrix": name,
            "Cycles_GMRES": cg,
            "Cycles_LGMRES": cl,
            "Cycles_GMRES_DMDc_LQR": c_dlk
        })

        time_rows.append({
            "Problem": pid,
            "Matrix": name,
            "avg_time_GMRES": tg,
            "avg_time_LGMRES": tl,
            "avg_time_GMRES_DMDc_LQR": t_dlk
        })

        print(f"{pid:<3} | {name:<14} | {n:<6} | {str(cg):<13} | {str(cl):<14} | {str(c_dlk):<11} | {str(t_dlk)}")

    # Guardar cycles_experiment.csv
    with open(out_cycles, "w", newline="") as f:
        fieldnames = ["Problem", "Matrix", "Cycles_GMRES", "Cycles_LGMRES", "Cycles_GMRES_DMDc_LQR"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(cycles_rows)

    # Guardar avg_relativa_time_experiment.csv
    with open(out_time, "w", newline="") as f:
        fieldnames = ["Problem", "Matrix", "avg_time_GMRES", "avg_time_LGMRES", "avg_time_GMRES_DMDc_LQR"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(time_rows)

    print("-" * 105)
    print(f"[✓] Benchmark finalizado exitosamente. Tablas guardadas en:")
    print(f"    - {out_matrix}")
    print(f"    - {out_cycles}")
    print(f"    - {out_time}")
    print("=" * 105)

if __name__ == "__main__":
    main()
