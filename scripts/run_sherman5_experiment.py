#!/usr/bin/env python3
"""
scripts/run_sherman5_experiment.py

Ejecuta el experimento completo sobre la matriz sherman5 bajo las mismas
condiciones de publicación (sin precondicionador, tol=1e-9, OMP_NUM_THREADS=1):
1. GMRES(110)
2. LGMRES(107, 3)
3. DLK-LQR-GMRES (m0=30, m_max=331)
4. DLK-LQR-LGMRES (m0=27, aug=3, m_max=331)

Exporta historiales de ciclo a ciclo y regenera las figuras estilo SIAM en results/figures/.
"""

import os
import subprocess
import sys
import re

def main():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    os.chdir(repo_root)

    matrix_bin = "data/bin/sherman5.bin"
    if not os.path.exists(matrix_bin):
        print(f"[Error] No se encontró {matrix_bin}")
        sys.exit(1)

    os.makedirs("results/tables", exist_ok=True)
    os.makedirs("results/figures", exist_ok=True)

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"

    solvers = [
        {
            "name": "GMRES(30)",
            "cycle_file": "results/tables/cycles_sherman5_gmres30.csv",
            "cmd": [
                "./test_dlk", "-matrix", matrix_bin,
                "-ksp_type", "gmres", "-ksp_gmres_restart", "30",
                "-pc_type", "none", "-tol", "1e-9", "-ksp_max_it", "1800",
                "-cycle_file", "results/tables/cycles_sherman5_gmres30.csv"
            ]
        },
        {
            "name": "LGMRES(28, 2)",
            "cycle_file": "results/tables/cycles_sherman5_lgmres28.csv",
            "cmd": [
                "./test_dlk", "-matrix", matrix_bin,
                "-ksp_type", "lgmres", "-ksp_lgmres_restart", "28", "-ksp_lgmres_augment", "2",
                "-pc_type", "none", "-tol", "1e-9", "-ksp_max_it", "1800",
                "-cycle_file", "results/tables/cycles_sherman5_lgmres28.csv"
            ]
        },
        {
            "name": "DMDc-LQR-GMRES",
            "cycle_file": "results/tables/cycles_sherman5_dlk.csv",
            "cmd": [
                "./test_dlk", "-matrix", matrix_bin,
                "-ksp_type", "dlk", "-ksp_dlk_controller", "lqr",
                "-ksp_dlk_m0", "30", "-ksp_dlk_m_min", "2", "-ksp_dlk_m_max", "331",
                "-ksp_dlk_alpha", "0.0001",
                "-pc_type", "none", "-tol", "1e-9", "-ksp_max_it", "50000",
                "-cycle_file", "results/tables/cycles_sherman5_dlk.csv"
            ]
        },
        {
            "name": "DMDc-LQR-LGMRES",
            "cycle_file": "results/tables/cycles_sherman5_dlklgmres.csv",
            "cmd": [
                "./test_dlk", "-matrix", matrix_bin,
                "-ksp_type", "dlk_lgmres", "-ksp_dlk_controller", "lqr",
                "-ksp_dlk_m0", "28", "-ksp_dlk_aug", "2", "-ksp_dlk_m_min", "2", "-ksp_dlk_m_max", "331",
                "-ksp_dlk_alpha", "0.0001",
                "-pc_type", "none", "-tol", "1e-9", "-ksp_max_it", "50000",
                "-cycle_file", "results/tables/cycles_sherman5_dlklgmres.csv"
            ]
        }
    ]

    print("==========================================================================")
    print(" EJECUTANDO EXPERIMENTO SHERMAN5 CON AJUSTES NUMÉRICOS Y OMP_NUM_THREADS=1")
    print("==========================================================================")

    for s in solvers:
        print(f"\n---> Ejecutando {s['name']} ...")
        proc = subprocess.run(s['cmd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
        if proc.returncode != 0:
            print(f"[Error] Falló ejecución de {s['name']}:\n{proc.stderr}")
            continue

        cyc_m = re.search(r"Ciclos Totales\s*:\s*(\d+)", proc.stdout)
        its_m = re.search(r"Iteraciones Totales\s*:\s*(\d+)", proc.stdout)
        time_m = re.search(r"Tiempo de Ejecución\s*:\s*([0-9.]+)", proc.stdout)
        res_m = re.search(r"Residuo \|\|b - Ax\|\|\s*:\s*([0-9.eE+-]+)", proc.stdout)

        c_val = cyc_m.group(1) if cyc_m else "?"
        it_val = its_m.group(1) if its_m else "?"
        t_val = time_m.group(1) if time_m else "?"
        r_val = res_m.group(1) if res_m else "?"
        print(f"     [OK] Ciclos: {c_val} | Iters: {it_val} | Tiempo: {t_val}s | ||r||: {r_val}")

    print("\n==========================================================================")
    print(" GENERANDO GRÁFICAS COMPARATIVAS ESTILO SIAM")
    print("==========================================================================")

    # 1. Comparison of methods
    print("-> Generando Comparison of methods (sherman5) ...")
    subprocess.run(["python3", "scripts/plot_comparison_of_methods.py"], env=env)

    # 2. Performance with m (DLK-GMRES)
    print("-> Generando Performance with m (DLK-GMRES sherman5) ...")
    subprocess.run([
        "python3", "scripts/plot_performance_with_m.py",
        "--cycle_file", "results/tables/cycles_sherman5_dlk.csv",
        "--output_prefix", "results/figures/PerformanceofDMDcLQRGMRESmsherman5_with_m",
        "--m_max", "331", "--m_min", "2",
        "--title", "sherman5"
    ], env=env)

    # 3. Performance with m (DLK-LGMRES)
    print("-> Generando Performance with m (DLK-LGMRES sherman5) ...")
    subprocess.run([
        "python3", "scripts/plot_performance_with_m.py",
        "--cycle_file", "results/tables/cycles_sherman5_dlklgmres.csv",
        "--output_prefix", "results/figures/PerformanceofDMDcLQRLGMRESmellsherman5_with_m",
        "--m_max", "331", "--m_min", "2",
        "--title", "sherman5"
    ], env=env)

    print("\n[✓] Experimento sherman5 completado y figuras regeneradas exitosamente.")

if __name__ == '__main__':
    main()
