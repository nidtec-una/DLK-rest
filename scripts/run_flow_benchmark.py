#!/usr/bin/env python3
"""
scripts/run_flow_benchmark.py

Ejecuta el benchmark completo de solucionadores Krylov (GMRES, LGMRES, DMDc-LQR-GMRES,
DMDc-LQR-LGMRES) sobre el conjunto de matrices de flujo en medios fracturados
(Mixed-Dimensional Problems with High Permeability Ratios, Kf = 10^-8 a 10^8).

Exporta tablas CSV idénticas a las del benchmark de SuiteSparse:
1. results/tables/matrix_info_flow.csv
2. results/tables/cycles_flow_benchmark.csv
3. results/tables/avg_relative_time_flow_benchmark.csv

Genera y guarda progresivamente:
- Historiales de ciclo: results/tables/cycles_flow_<tag>_<solver>.csv
- Gráficos de performance en 3 paneles: results/figures/PerformanceofDMDcLQRGMRES_flow_<tag>_with_m.png / .eps
- Gráficos comparativos de barras: results/figures/bar_plot_flow_benchmark_cycles.png / .eps
"""

import os
import sys
import subprocess
import csv
import re
import time
import argparse
import numpy as np

# Configuración del conjunto de 9 matrices de permeabilidad
FLOW_PROBLEMS = [
    {
        "id": 1,
        "name": "flow_benchmark_2d_case_1_Kf1e-08",
        "tag": "Kf1e-08",
        "Kf_str": "1.0e-08",
        "regime": "Barrera casi impermeable",
        "file": "data/bin/flow_benchmark/flow_benchmark_2d_case_1_Kf1e-08.bin",
        "N": 11089,
        "nnz": 76208,
        "norm_b": 2.0e+01
    },
    {
        "id": 2,
        "name": "flow_benchmark_2d_case_1_Kf1e-06",
        "tag": "Kf1e-06",
        "Kf_str": "1.0e-06",
        "regime": "Fractura de baja conductividad",
        "file": "data/bin/flow_benchmark/flow_benchmark_2d_case_1_Kf1e-06.bin",
        "N": 11089,
        "nnz": 76208,
        "norm_b": 2.0e+01
    },
    {
        "id": 3,
        "name": "flow_benchmark_2d_case_1_Kf1e-04",
        "tag": "Kf1e-04",
        "Kf_str": "1.0e-04",
        "regime": "Fractura de resistencia moderada",
        "file": "data/bin/flow_benchmark/flow_benchmark_2d_case_1_Kf1e-04.bin",
        "N": 11089,
        "nnz": 76208,
        "norm_b": 2.0e+01
    },
    {
        "id": 4,
        "name": "flow_benchmark_2d_case_1_Kf1e-02",
        "tag": "Kf1e-02",
        "Kf_str": "1.0e-02",
        "regime": "Transición matriz-fractura",
        "file": "data/bin/flow_benchmark/flow_benchmark_2d_case_1_Kf1e-02.bin",
        "N": 11089,
        "nnz": 76208,
        "norm_b": 2.0e+01
    },
    {
        "id": 5,
        "name": "flow_benchmark_2d_case_1_Kf1e+00",
        "tag": "Kf1e+00",
        "Kf_str": "1.0e+00",
        "regime": "Permeabilidad homogénea",
        "file": "data/bin/flow_benchmark/flow_benchmark_2d_case_1_Kf1e+00.bin",
        "N": 11089,
        "nnz": 76208,
        "norm_b": 2.0e+01
    },
    {
        "id": 6,
        "name": "flow_benchmark_2d_case_1_Kf1e+02",
        "tag": "Kf1e+02",
        "Kf_str": "1.0e+02",
        "regime": "Fractura moderadamente conductora",
        "file": "data/bin/flow_benchmark/flow_benchmark_2d_case_1_Kf1e+02.bin",
        "N": 11089,
        "nnz": 76208,
        "norm_b": 2.0199e+01
    },
    {
        "id": 7,
        "name": "flow_benchmark_2d_case_1_Kf1e+04",
        "tag": "Kf1e+04",
        "Kf_str": "1.0e+04",
        "regime": "Fractura conductora",
        "file": "data/bin/flow_benchmark/flow_benchmark_2d_case_1_Kf1e+04.bin",
        "N": 11089,
        "nnz": 76208,
        "norm_b": 2.8355e+02
    },
    {
        "id": 8,
        "name": "flow_benchmark_2d_case_1_Kf1e+06",
        "tag": "Kf1e+06",
        "Kf_str": "1.0e+06",
        "regime": "Canal altamente conductor",
        "file": "data/bin/flow_benchmark/flow_benchmark_2d_case_1_Kf1e+06.bin",
        "N": 11089,
        "nnz": 76208,
        "norm_b": 2.8284e+04
    },
    {
        "id": 9,
        "name": "flow_benchmark_2d_case_1_Kf1e+08",
        "tag": "Kf1e+08",
        "Kf_str": "1.0e+08",
        "regime": "Canal de flujo hiper-conductor",
        "file": "data/bin/flow_benchmark/flow_benchmark_2d_case_1_Kf1e+08.bin",
        "N": 11089,
        "nnz": 76208,
        "norm_b": 2.8284e+06
    }
]

def run_petsc_solver(bin_file, solver_type, cycle_file=None, norm_b=1.0, m0=30, m_max=60, l=2, tol=1e-9, max_cycles=400, max_its=50000, timeout=120):
    cmd = [
        "./test_dlk",
        "-matrix", bin_file,
        "-tol", str(tol),
        "-ksp_max_it", str(max_its),
        "-pc_type", "none"
    ]

    if cycle_file:
        cmd += ["-cycle_file", cycle_file]

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
        return "NC", "NC", f"Timeout (>{timeout}s)"
    except Exception as e:
        return "NC", "NC", str(e)

    if proc.returncode != 0:
        return "NC", "NC", f"Error {proc.returncode}"

    cyc_m = re.search(r"Ciclos Totales\s*:\s*(\d+)", output)
    its_m = re.search(r"Iteraciones Totales\s*:\s*(\d+)", output)
    res_m = re.search(r"Residuo \|\|b - Ax\|\|\s*:\s*([0-9.eE+-]+)", output)
    time_m = re.search(r"Tiempo de Ejecución\s*:\s*([0-9.]+)", output)

    cycles = int(cyc_m.group(1)) if cyc_m else -1
    its = int(its_m.group(1)) if its_m else -1
    residual = float(res_m.group(1)) if res_m else 1.0
    report_time = float(time_m.group(1)) if time_m else elapsed

    # Criterio estricto de convergencia relativa: ||r_k|| / ||b|| <= tol * 50
    rel_res = residual / norm_b if norm_b > 0 else residual
    converged = (cycles > 0) and (cycles < max_cycles) and (its < max_its) and (rel_res <= tol * 50)
    if not converged:
        return "NC", "NC", f"No conv (its={its}, rel_res={rel_res:.2e})"

    return cycles, round(report_time, 5), f"Conv (its={its}, rel_res={rel_res:.2e})"

def save_csv_tables(out_matrix, out_cycles, out_time, matrix_rows, cycles_rows, time_rows):
    """Guarda inmediatamente las tablas CSV para asegurar persistencia incremental."""
    with open(out_matrix, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["Problem", "Matrix", "Tag", "Kf", "Regime", "n", "nnz", "norm_b"])
        writer.writeheader()
        writer.writerows(matrix_rows)

    with open(out_cycles, "w", newline="") as f:
        fieldnames = ["Problem", "Matrix", "Tag", "Kf", "Cycles_GMRES", "Cycles_LGMRES", "Cycles_GMRES_DMDc_LQR", "Cycles_LGMRES_DMDc_LQR"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(cycles_rows)

    with open(out_time, "w", newline="") as f:
        fieldnames = ["Problem", "Matrix", "Tag", "Kf", "avg_time_GMRES", "avg_time_LGMRES", "avg_time_GMRES_DMDc_LQR", "avg_time_LGMRES_DMDc_LQR"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(time_rows)

def generate_performance_figure(cycle_file, output_prefix, title, m_max=60, m_min=2):
    """Invoca scripts/plot_performance_with_m.py para generar figuras de 3 paneles."""
    if not os.path.exists(cycle_file):
        return False
    cmd = [
        "python3", "scripts/plot_performance_with_m.py",
        "--cycle_file", cycle_file,
        "--output_prefix", output_prefix,
        "--m_max", str(m_max),
        "--m_min", str(m_min),
        "--title", title
    ]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except Exception as e:
        print(f"      [!] Error generando figura {output_prefix}: {e}")
        return False

def generate_comparative_bar_chart(cycles_csv, time_csv):
    """Genera gráficos de barras comparativos para el benchmark de flujo."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        tags = []
        c_gmres, c_lgmres, c_dlk_gmres, c_dlk_lgmres = [], [], [], []
        t_gmres, t_lgmres, t_dlk_gmres, t_dlk_lgmres = [], [], [], []

        def parse(v):
            if v is None or v == "" or v == "NC":
                return np.nan
            return float(v)

        with open(cycles_csv, 'r') as f:
            for row in csv.DictReader(f):
                tags.append(row["Tag"])
                c_gmres.append(parse(row.get("Cycles_GMRES")))
                c_lgmres.append(parse(row.get("Cycles_LGMRES")))
                c_dlk_gmres.append(parse(row.get("Cycles_GMRES_DMDc_LQR")))
                c_dlk_lgmres.append(parse(row.get("Cycles_LGMRES_DMDc_LQR")))

        with open(time_csv, 'r') as f:
            for row in csv.DictReader(f):
                t_gmres.append(parse(row.get("avg_time_GMRES")))
                t_lgmres.append(parse(row.get("avg_time_LGMRES")))
                t_dlk_gmres.append(parse(row.get("avg_time_GMRES_DMDc_LQR")))
                t_dlk_lgmres.append(parse(row.get("avg_time_LGMRES_DMDc_LQR")))

        x = np.arange(len(tags))
        width = 0.20

        # 1. Gráfico de Ciclos
        fig, ax = plt.subplots(figsize=(10, 4.8))
        ax.bar(x - 1.5*width, c_gmres, width, label='GMRES(30)', color='#404040')
        ax.bar(x - 0.5*width, c_lgmres, width, label='LGMRES(28, 2)', color='#707070')
        ax.bar(x + 0.5*width, c_dlk_gmres, width, label='DMDc-LQR-GMRES', color='#1f77b4')
        ax.bar(x + 1.5*width, c_dlk_lgmres, width, label='DMDc-LQR-LGMRES', color='#ff7f0e')

        ax.set_ylabel('Ciclos de Reinicio (escala log)', fontsize=11)
        ax.set_title('Flow Benchmark: Número de Ciclos vs Contraste de Permeabilidad $K_f$', fontsize=12, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(tags, fontsize=10, rotation=25)
        ax.set_yscale('log')
        ax.grid(True, which='both', linestyle=':', alpha=0.6)
        ax.legend(frameon=True, fontsize=10)
        fig.tight_layout()

        fig.savefig("results/figures/bar_plot_flow_benchmark_cycles.png", dpi=300)
        fig.savefig("results/figures/bar_plot_flow_benchmark_cycles.eps", format='eps')
        plt.close(fig)

        # 2. Gráfico de Tiempos
        fig, ax = plt.subplots(figsize=(10, 4.8))
        ax.bar(x - 1.5*width, t_gmres, width, label='GMRES(30)', color='#404040')
        ax.bar(x - 0.5*width, t_lgmres, width, label='LGMRES(28, 2)', color='#707070')
        ax.bar(x + 0.5*width, t_dlk_gmres, width, label='DMDc-LQR-GMRES', color='#1f77b4')
        ax.bar(x + 1.5*width, t_dlk_lgmres, width, label='DMDc-LQR-LGMRES', color='#ff7f0e')

        ax.set_ylabel('Tiempo Computacional [s] (escala log)', fontsize=11)
        ax.set_title('Flow Benchmark: Tiempo de Cómputo vs Contraste de Permeabilidad $K_f$', fontsize=12, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(tags, fontsize=10, rotation=25)
        ax.set_yscale('log')
        ax.grid(True, which='both', linestyle=':', alpha=0.6)
        ax.legend(frameon=True, fontsize=10)
        fig.tight_layout()

        fig.savefig("results/figures/bar_plot_flow_benchmark_time.png", dpi=300)
        fig.savefig("results/figures/bar_plot_flow_benchmark_time.eps", format='eps')
        plt.close(fig)
        print("   [✓] Gráficos de barras comparativos guardados en results/figures/")
    except Exception as e:
        print(f"   [!] Error generando gráficos comparativos de barras: {e}")

def main():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    os.chdir(repo_root)

    os.makedirs("results/tables", exist_ok=True)
    os.makedirs("results/figures", exist_ok=True)

    out_matrix = "results/tables/matrix_info_flow.csv"
    out_cycles = "results/tables/cycles_flow_benchmark.csv"
    out_time = "results/tables/avg_relative_time_flow_benchmark.csv"

    parser = argparse.ArgumentParser(description="Ejecutar benchmark completo de problemas de flujo fracturado")
    parser.add_argument("--force", action="store_true", help="Forzar re-ejecución ignorando resultados previos")
    parser.add_argument("--timeout", type=int, default=120, help="Timeout por corrida en segundos (def: 120)")
    args = parser.parse_args()

    print("\n" + "=" * 115)
    print(" BENCHMARK COMPLETO: MIXED-DIMENSIONAL FLOW BENCHMARK EN MEDIOS FRACTURADOS")
    print(" (Evaluación de 4 métodos Krylov a través de 9 órdenes de magnitud de permeabilidad Kf)")
    print("=" * 115)
    print(" Condiciones de ejecución:")
    print(" - Hilos BLAS/OpenMP : OMP_NUM_THREADS = 1 (secuencial determinista)")
    print(" - Precondicionador  : Ninguno (-pc_type none, unificado con artículo)")
    print(" - Tolerancia        : tol = 1e-9 (estricta)")
    print(" - Parámetro DLK     : m0 = 30, m_min = 2, m_max = 60, alpha = 1e-4")
    print("=" * 115)
    sys.stdout.flush()

    matrix_rows = []
    cycles_rows = []
    time_rows = []

    # Cargar datos previos si existen y no es --force
    existing_cycles = {}
    existing_time = {}
    if not args.force and os.path.exists(out_cycles) and os.path.exists(out_time):
        with open(out_cycles, 'r') as f:
            for r in csv.DictReader(f):
                existing_cycles[int(r["Problem"])] = r
        with open(out_time, 'r') as f:
            for r in csv.DictReader(f):
                existing_time[int(r["Problem"])] = r

    for p in FLOW_PROBLEMS:
        pid = p["id"]
        name = p["name"]
        tag = p["tag"]
        kf_str = p["Kf_str"]
        regime = p["regime"]
        bin_file = p["file"]
        N = p["N"]
        nnz = p["nnz"]
        norm_b = p["norm_b"]

        matrix_rows.append({
            "Problem": pid,
            "Matrix": name,
            "Tag": tag,
            "Kf": kf_str,
            "Regime": regime,
            "n": N,
            "nnz": nnz,
            "norm_b": f"{norm_b:.4e}"
        })

        print(f"\n[{pid}/9] Problema: {tag} ({name})")
        print(f"      Física: {regime} | N = {N}, nnz = {nnz}, ||b|| = {norm_b:.2e}")
        sys.stdout.flush()

        if pid in existing_cycles and pid in existing_time:
            print(f"      [i] Usando resultados previos en caché (usar --force para re-ejecutar)")
            c_row = existing_cycles[pid]
            t_row = existing_time[pid]
            cycles_rows.append(c_row)
            time_rows.append(t_row)
            print(f"      GMRES(30) : {c_row['Cycles_GMRES']} ciclos | {t_row['avg_time_GMRES']} s")
            print(f"      LGMRES    : {c_row['Cycles_LGMRES']} ciclos | {t_row['avg_time_LGMRES']} s")
            print(f"      DLK-GMRES : {c_row['Cycles_GMRES_DMDc_LQR']} ciclos | {t_row['avg_time_GMRES_DMDc_LQR']} s")
            print(f"      DLK-LGMRES: {c_row['Cycles_LGMRES_DMDc_LQR']} ciclos | {t_row['avg_time_LGMRES_DMDc_LQR']} s")
            sys.stdout.flush()
            continue

        cyc_file_gmres = f"results/tables/cycles_flow_{tag}_gmres30.csv"
        cyc_file_lgmres = f"results/tables/cycles_flow_{tag}_lgmres28.csv"
        cyc_file_dlk = f"results/tables/cycles_flow_{tag}_dlk.csv"
        cyc_file_dlk_l = f"results/tables/cycles_flow_{tag}_dlk_lgmres.csv"

        # 1. GMRES(30)
        print("      -> Ejecutando GMRES(30) ... ", end="", flush=True)
        cg, tg, msg_g = run_petsc_solver(bin_file, "GMRES", cycle_file=cyc_file_gmres, norm_b=norm_b, m0=30, timeout=args.timeout)
        print(f"{cg} ciclos | {tg} s ({msg_g})")
        sys.stdout.flush()

        # 2. LGMRES(28, 2)
        print("      -> Ejecutando LGMRES(28, 2) ... ", end="", flush=True)
        cl, tl, msg_l = run_petsc_solver(bin_file, "LGMRES", cycle_file=cyc_file_lgmres, norm_b=norm_b, m0=30, l=2, timeout=args.timeout)
        print(f"{cl} ciclos | {tl} s ({msg_l})")
        sys.stdout.flush()

        # 3. DMDc-LQR-GMRES
        print("      -> Ejecutando DMDc-LQR-GMRES (m0=30, m_max=60) ... ", end="", flush=True)
        c_dlk, t_dlk, msg_dlk = run_petsc_solver(bin_file, "GMRES_DMDc_LQR", cycle_file=cyc_file_dlk, norm_b=norm_b, m0=30, m_max=60, timeout=args.timeout)
        print(f"{c_dlk} ciclos | {t_dlk} s ({msg_dlk})")
        sys.stdout.flush()

        # 4. DMDc-LQR-LGMRES
        print("      -> Ejecutando DMDc-LQR-LGMRES (m0=28, l=2, m_max=60) ... ", end="", flush=True)
        c_dlkl, t_dlkl, msg_dlkl = run_petsc_solver(bin_file, "LGMRES_DMDc_LQR", cycle_file=cyc_file_dlk_l, norm_b=norm_b, m0=30, l=2, m_max=60, timeout=args.timeout)
        print(f"{c_dlkl} ciclos | {t_dlkl} s ({msg_dlkl})")
        sys.stdout.flush()

        c_entry = {
            "Problem": pid,
            "Matrix": name,
            "Tag": tag,
            "Kf": kf_str,
            "Cycles_GMRES": cg,
            "Cycles_LGMRES": cl,
            "Cycles_GMRES_DMDc_LQR": c_dlk,
            "Cycles_LGMRES_DMDc_LQR": c_dlkl
        }
        t_entry = {
            "Problem": pid,
            "Matrix": name,
            "Tag": tag,
            "Kf": kf_str,
            "avg_time_GMRES": tg,
            "avg_time_LGMRES": tl,
            "avg_time_GMRES_DMDc_LQR": t_dlk,
            "avg_time_LGMRES_DMDc_LQR": t_dlkl
        }
        cycles_rows.append(c_entry)
        time_rows.append(t_entry)

        # Guardado incremental inmediato en disco para evitar pérdida de datos
        save_csv_tables(out_matrix, out_cycles, out_time, matrix_rows, cycles_rows, time_rows)
        print(f"      [✓] Progreso persistido incrementalmente en results/tables/")

        # Generar figura de 3 paneles para este caso
        fig_prefix = f"results/figures/PerformanceofDMDcLQRGMRES_flow_{tag}_with_m"
        title_fig = f"Flow Benchmark ({tag})"
        ok_fig = generate_performance_figure(cyc_file_dlk, fig_prefix, title_fig, m_max=60, m_min=2)
        if ok_fig:
            print(f"      [✓] Figura generada: {fig_prefix}.png / .eps")
        sys.stdout.flush()

    # Guardado final de tablas
    save_csv_tables(out_matrix, out_cycles, out_time, matrix_rows, cycles_rows, time_rows)

    # Generar gráficos comparativos de barras
    print("\n" + "=" * 115)
    print(" GENERANDO GRÁFICOS COMPARATIVOS DE BARRAS")
    print("=" * 115)
    generate_comparative_bar_chart(out_cycles, out_time)

    # Resumen en consola estilo LaTeX
    print("\n" + "=" * 115)
    print(" RESUMEN FINAL DEL BENCHMARK DE FLUJO FRACTURADO (4 MÉTODOS)")
    print("=" * 115)
    print(f"{'#':<2} | {'Tag':<10} | {'Regimen':<34} | {'GMRES(30)':<10} | {'LGMRES(28,2)':<13} | {'DLK-GMRES':<10} | {'DLK-LGMRES':<11}")
    print("-" * 115)
    for c in cycles_rows:
        pid = c["Problem"]
        tag = c["Tag"]
        regime = FLOW_PROBLEMS[pid-1]["regime"]
        print(f"{pid:<2} | {tag:<10} | {regime:<34} | {str(c['Cycles_GMRES']):<10} | {str(c['Cycles_LGMRES']):<13} | {str(c['Cycles_GMRES_DMDc_LQR']):<10} | {str(c['Cycles_LGMRES_DMDc_LQR']):<11}")
    print("=" * 115)
    print(f"[✓] Tablas guardadas en:")
    print(f"    - {out_matrix}")
    print(f"    - {out_cycles}")
    print(f"    - {out_time}")
    print(f"[✓] Todas las figuras vectoriales guardadas en results/figures/")
    print("=" * 115)

if __name__ == "__main__":
    main()
