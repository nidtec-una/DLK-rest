#!/usr/bin/env python3
"""
scripts/plot_benchmark_bars.py

Lee results/tables/benchmark_comparison.csv y genera gráficos de barras comparativos
para publicaciones (estilo SIAM / IEEE):
1. results/figures/benchmark_bars_iterations.png / .eps (Iteraciones por solver y matriz)
2. results/figures/benchmark_bars_time.png / .eps (Tiempo de CPU en segundos)
"""

import os
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_benchmark_bars(csv_path="results/tables/benchmark_comparison.csv",
                        output_dir="results/figures"):
    os.makedirs(output_dir, exist_ok=True)
    
    if not os.path.exists(csv_path):
        print(f"[Error] No se encontró el archivo {csv_path}")
        return

    # Leer datos del CSV
    data = {}
    solvers = []
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            mat = row["Matrix"]
            s_label = row["Solver_Label"]
            its = int(row["Iterations"])
            t_sec = float(row["Time_sec"])
            
            if s_label not in solvers:
                solvers.append(s_label)
                
            if mat not in data:
                data[mat] = {}
            data[mat][s_label] = {"iterations": its, "time": t_sec}

    matrices = list(data.keys())
    num_mats = len(matrices)
    num_solvers = len(solvers)
    
    if num_mats == 0 or num_solvers == 0:
        print("[Error] No hay datos suficientes en el CSV para graficar.")
        return

    # Configuración de colores y estilos profesionales (paleta editorial SIAM / IEEE)
    colors = ['#7f7f7f', '#3b528b', '#21918c', '#5ec962', '#e41a1c']
    if len(solvers) > len(colors):
        colors = plt.cm.tab10(np.linspace(0, 1, len(solvers)))

    x = np.arange(num_mats)
    bar_width = 0.8 / num_solvers

    # -------------------------------------------------------------
    # 1. Gráfico de Barras: Iteraciones Totales (Escala Logarítmica)
    # -------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8.0, 4.5))
    fig.patch.set_facecolor('white')

    for i, s_label in enumerate(solvers):
        its_vals = [max(1, data[m].get(s_label, {}).get("iterations", 0)) for m in matrices]
        pos = x - 0.4 + (i + 0.5) * bar_width
        ax.bar(pos, its_vals, width=bar_width, label=s_label, color=colors[i % len(colors)],
               edgecolor='black', linewidth=0.6, alpha=0.9)

    ax.set_ylabel("Iteraciones Totales (escala log)", fontsize=11, fontweight='semibold')
    ax.set_title("Comparación de Iteraciones de Solvers Krylov en PETSc", fontsize=12, pad=10)
    ax.set_xticks(x)
    ax.set_xticklabels(matrices, fontsize=10, rotation=25, ha='right')
    ax.set_yscale('log')
    ax.grid(axis='y', linestyle=':', alpha=0.7)
    ax.legend(frameon=True, fontsize=9, loc='upper right', ncol=min(3, num_solvers))
    ax.tick_params(direction='in', top=True, right=True, which='both')

    plt.tight_layout()
    fig_iter_png = os.path.join(output_dir, "benchmark_bars_iterations.png")
    fig_iter_eps = os.path.join(output_dir, "benchmark_bars_iterations.eps")
    plt.savefig(fig_iter_png, dpi=300, bbox_inches='tight')
    plt.savefig(fig_iter_eps, format='eps', bbox_inches='tight')
    plt.close()

    # -------------------------------------------------------------
    # 2. Gráfico de Barras: Tiempo de Ejecución (Segundos)
    # -------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8.0, 4.5))
    fig.patch.set_facecolor('white')

    for i, s_label in enumerate(solvers):
        time_vals = [max(1e-4, data[m].get(s_label, {}).get("time", 0.0)) for m in matrices]
        pos = x - 0.4 + (i + 0.5) * bar_width
        ax.bar(pos, time_vals, width=bar_width, label=s_label, color=colors[i % len(colors)],
               edgecolor='black', linewidth=0.6, alpha=0.9)

    ax.set_ylabel("Tiempo de Ejecución [s] (escala log)", fontsize=11, fontweight='semibold')
    ax.set_title("Comparación de Tiempo de Cómputo en PETSc", fontsize=12, pad=10)
    ax.set_xticks(x)
    ax.set_xticklabels(matrices, fontsize=10, rotation=25, ha='right')
    ax.set_yscale('log')
    ax.grid(axis='y', linestyle=':', alpha=0.7)
    ax.legend(frameon=True, fontsize=9, loc='upper right', ncol=min(3, num_solvers))
    ax.tick_params(direction='in', top=True, right=True, which='both')

    plt.tight_layout()
    fig_time_png = os.path.join(output_dir, "benchmark_bars_time.png")
    fig_time_eps = os.path.join(output_dir, "benchmark_bars_time.eps")
    plt.savefig(fig_time_png, dpi=300, bbox_inches='tight')
    plt.savefig(fig_time_eps, format='eps', bbox_inches='tight')
    plt.close()

    print("[✓] Gráficos de barras comparativos generados exitosamente:")
    print(f"    - Iteraciones: {fig_iter_png} y {fig_iter_eps}")
    print(f"    - Tiempo:      {fig_time_png} y {fig_time_eps}")

if __name__ == "__main__":
    plot_benchmark_bars()
