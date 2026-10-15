#!/usr/bin/env python3
"""
scripts/run_young3c_experiment.py

Ejecuta el experimento comparativo completo sobre la matriz young3c bajo las mismas
condiciones unificadas del benchmark del artículo:
- Sin precondicionador (-pc_type none)
- Tolerancia estricta (-tol 1e-9)
- OMP_NUM_THREADS=1
- Métodos evaluados:
  * GMRES(30)
  * LGMRES(28, 2)
  * DMDc-LQR-GMRES (m0=30, m_min=2, m_max=84, alpha=1e-4)

Posteriormente genera las figuras para publicación:
1. results/figures/Comparisonofmethodsyoung3c (PNG y EPS)
2. results/figures/PerformanceofDMDcLQRGMRESmyoung3c_with_m (PNG y EPS)
"""

import os
import subprocess
import sys

def run_command(cmd, desc):
    print(f"\n=======================================================")
    print(f"[>] {desc}")
    print(f"    Comando: {' '.join(cmd)}")
    print(f"=======================================================")
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    res = subprocess.run(cmd, env=env)
    if res.returncode != 0:
        print(f"[Error] Falló la ejecución de: {desc} (código {res.returncode})", file=sys.stderr)
        return False
    return True

def main():
    matrix_bin = "data/bin/young3c.bin"
    if not os.path.exists(matrix_bin):
        print(f"[Error] No se encontró el archivo binario: {matrix_bin}", file=sys.stderr)
        print("Ejecuta primero: python3 data/convert_mat_to_bin.py")
        sys.exit(1)

    os.makedirs("results/tables", exist_ok=True)
    os.makedirs("results/figures", exist_ok=True)

    solvers = [
        {
            "name": "GMRES(30)",
            "cycle_file": "results/tables/cycles_young3c_gmres30.csv",
            "cmd": [
                "./test_dlk", "-matrix", matrix_bin,
                "-ksp_type", "gmres", "-ksp_gmres_restart", "30",
                "-pc_type", "none", "-tol", "1e-9", "-ksp_max_it", "50000",
                "-cycle_file", "results/tables/cycles_young3c_gmres30.csv"
            ]
        },
        {
            "name": "LGMRES(28, 2)",
            "cycle_file": "results/tables/cycles_young3c_lgmres28.csv",
            "cmd": [
                "./test_dlk", "-matrix", matrix_bin,
                "-ksp_type", "lgmres", "-ksp_lgmres_restart", "28", "-ksp_lgmres_augment", "2",
                "-pc_type", "none", "-tol", "1e-9", "-ksp_max_it", "50000",
                "-cycle_file", "results/tables/cycles_young3c_lgmres28.csv"
            ]
        },
        {
            "name": "DMDc-LQR-GMRES",
            "cycle_file": "results/tables/cycles_young3c_dlk.csv",
            "cmd": [
                "./test_dlk", "-matrix", matrix_bin,
                "-ksp_type", "dlk", "-ksp_dlk_controller", "lqr",
                "-ksp_dlk_m0", "30", "-ksp_dlk_m_min", "2", "-ksp_dlk_m_max", "84",
                "-ksp_dlk_alpha", "0.0001",
                "-pc_type", "none", "-tol", "1e-9", "-ksp_max_it", "50000",
                "-cycle_file", "results/tables/cycles_young3c_dlk.csv"
            ]
        }
    ]

    print("\n" + "=" * 65)
    print(" EJECUTANDO EXPERIMENTO COMPARATIVO CON ORSIRR1")
    print("=" * 65)

    for s in solvers:
        ok = run_command(s["cmd"], f"Ejecutando {s['name']} sobre young3c")
        if not ok:
            print(f"[Aviso] {s['name']} no completó exitosamente.")

    print("\n" + "=" * 65)
    print(" REGENERANDO FIGURAS PARA PUBLICACIÓN (ORSIRR1)")
    print("=" * 65)

    # 1. Comparison of methods
    print("-> Generando Comparison of methods (young3c) ...")
    cmd_comp = [
        "python3", "scripts/plot_comparison_of_methods.py",
        "--matrix", "young3c",
        "--output_prefix", "results/figures/Comparisonofmethodsyoung3c",
        "--title", "young3c"
    ]
    subprocess.run(cmd_comp, check=True)

    # 2. Performance with m (DLK-GMRES)
    print("-> Generando Performance with m (DLK-GMRES young3c) ...")
    cmd_perf_gmres = [
        "python3", "scripts/plot_performance_with_m.py",
        "--cycle_file", "results/tables/cycles_young3c_dlk.csv",
        "--output_prefix", "results/figures/PerformanceofDMDcLQRGMRESmyoung3c_with_m",
        "--m_max", "84", "--m_min", "2",
        "--title", "young3c"
    ]
    subprocess.run(cmd_perf_gmres, check=True)

    print("\n[✓] Experimento young3c completado y figuras regeneradas exitosamente.")

if __name__ == "__main__":
    main()
