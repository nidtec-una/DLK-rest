#!/usr/bin/env python3
"""
scripts/plot_benchmark_bars.py

Genera los gráficos de barras comparativos para el artículo en LaTeX:
1. results/figures/bar_plot_Number_of_cycles.png y .eps
   (Número de ciclos para cada problema y método, a partir de results/tables/cycles_experiment.csv)
2. results/figures/bar_plot_Average_computation_relative_time.png y .eps
   (Tiempo computacional promedio relativo a GMRES, a partir de results/tables/avg_relativa_time_experiment.csv)
"""

import os
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_number_of_cycles(csv_path="results/tables/cycles_experiment.csv",
                          output_prefix="results/figures/bar_plot_Number_of_cycles"):
    if not os.path.exists(csv_path):
        print(f"[Aviso] No se encontró {csv_path}")
        return

    problems = []
    matrix_names = []
    c_gmres = []
    c_lgmres = []
    c_dlk = []
    c_dlkl = []

    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            problems.append(int(row["Problem"]))
            matrix_names.append(row["Matrix"])

            def parse_val(v):
                if v == "NC" or v == "" or v is None:
                    return np.nan
                try:
                    return float(v)
                except ValueError:
                    return np.nan

            c_gmres.append(parse_val(row.get("Cycles_GMRES")))
            c_lgmres.append(parse_val(row.get("Cycles_LGMRES")))
            c_dlk.append(parse_val(row.get("Cycles_GMRES_DMDc_LQR")))
            c_dlkl.append(parse_val(row.get("Cycles_LGMRES_DMDc_LQR")))

    n_probs = len(problems)
    x = np.arange(n_probs)
    width = 0.20

    # Tamaño estilo SIAM / IEEE: ancho ~7.0 in, alto ~3.8 in
    fig, ax = plt.subplots(figsize=(7.5, 3.8))
    fig.patch.set_facecolor('white')

    # Paleta editorial sobria para publicación
    # GMRES (gris claro), LGMRES (gris medio), DMDc-GMRES (azul petróleo), DMDc-LGMRES (verde oscuro)
    color_gmres = '#b0b0b0'
    color_lgmres = '#6e6e6e'
    color_dlk = '#1f77b4'
    color_dlkl = '#2ca02c'

    # Reemplazar NaN con 0 para plotting sin romper la escala log
    def clean_for_log(arr):
        return [val if (not np.isnan(val) and val > 0) else 0.001 for val in arr]

    b1 = ax.bar(x - 1.5 * width, clean_log(c_gmres), width, label=r'GMRES($m$)', color=color_gmres, edgecolor='black', linewidth=0.5)
    b2 = ax.bar(x - 0.5 * width, clean_log(c_lgmres), width, label=r'LGMRES($m,\ell$)', color=color_lgmres, edgecolor='black', linewidth=0.5)
    b3 = ax.bar(x + 0.5 * width, clean_log(c_dlk), width, label=r'DMDc-LQR-GMRES($m$)', color=color_dlk, edgecolor='black', linewidth=0.5)
    b4 = ax.bar(x + 1.5 * width, clean_log(c_dlkl), width, label=r'DMDc-LQR-LGMRES($m,\ell$)', color=color_dlkl, edgecolor='black', linewidth=0.5)

    ax.set_ylabel('Number of cycles (log scale)', fontsize=10, fontweight='semibold')
    ax.set_xlabel('Problem Number', fontsize=10, fontweight='semibold')
    ax.set_yscale('log')
    ax.set_ylim([0.8, 2000])

    ax.set_xticks(x)
    ax.set_xticklabels([str(p) for p in problems], fontsize=9)
    ax.grid(axis='y', which='both', linestyle=':', linewidth=0.5, alpha=0.6)
    ax.legend(frameon=True, fontsize=8.5, loc='upper right', ncol=2)
    ax.tick_params(direction='in', top=True, right=True, which='both')

    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(output_prefix)), exist_ok=True)
    png_path = f"{output_prefix}.png"
    eps_path = f"{output_prefix}.eps"
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(eps_path, format='eps', bbox_inches='tight')
    plt.close()
    print(f"[✓] Gráfico de ciclos guardado:\n    {png_path}\n    {eps_path}")

def clean_log(arr):
    return [val if (not np.isnan(val) and val > 0) else 0.001 for val in arr]

def plot_relative_time(csv_path="results/tables/avg_relativa_time_experiment.csv",
                       output_prefix="results/figures/bar_plot_Average_computation_relative_time"):
    if not os.path.exists(csv_path):
        print(f"[Aviso] No se encontró {csv_path}")
        return

    problems = []
    matrix_names = []
    t_gmres = []
    t_lgmres = []
    t_dlk = []
    t_dlkl = []

    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            problems.append(int(row["Problem"]))
            matrix_names.append(row["Matrix"])

            def parse_val(v):
                if v == "NC" or v == "" or v is None:
                    return np.nan
                try:
                    return float(v)
                except ValueError:
                    return np.nan

            t_gmres.append(parse_val(row.get("avg_time_GMRES")))
            t_lgmres.append(parse_val(row.get("avg_time_LGMRES")))
            t_dlk.append(parse_val(row.get("avg_time_GMRES_DMDc_LQR")))
            t_dlkl.append(parse_val(row.get("avg_time_LGMRES_DMDc_LQR")))

    n_probs = len(problems)
    x = np.arange(n_probs)
    width = 0.20

    # Calcular tiempos relativos respecto a GMRES(m) si GMRES convergió
    rel_gmres = []
    rel_lgmres = []
    rel_dlk = []
    rel_dlkl = []

    for i in range(n_probs):
        tg = t_gmres[i]
        tl = t_lgmres[i]
        td = t_dlk[i]
        tdl = t_dlkl[i]

        if not np.isnan(tg) and tg > 0:
            rel_gmres.append(1.0)
            rel_lgmres.append(tl / tg if not np.isnan(tl) else np.nan)
            rel_dlk.append(td / tg if not np.isnan(td) else np.nan)
            rel_dlkl.append(tdl / tg if not np.isnan(tdl) else np.nan)
        else:
            # Si GMRES no convergió, tomar como base el valor absoluto o normalizado al mínimo
            rel_gmres.append(np.nan)
            rel_lgmres.append(tl if not np.isnan(tl) else np.nan)
            rel_dlk.append(td if not np.isnan(td) else np.nan)
            rel_dlkl.append(tdl if not np.isnan(tdl) else np.nan)

    fig, ax = plt.subplots(figsize=(7.5, 3.8))
    fig.patch.set_facecolor('white')

    color_gmres = '#b0b0b0'
    color_lgmres = '#6e6e6e'
    color_dlk = '#1f77b4'
    color_dlkl = '#2ca02c'

    b1 = ax.bar(x - 1.5 * width, clean_log(rel_gmres), width, label=r'GMRES($m$)', color=color_gmres, edgecolor='black', linewidth=0.5)
    b2 = ax.bar(x - 0.5 * width, clean_log(rel_lgmres), width, label=r'LGMRES($m,\ell$)', color=color_lgmres, edgecolor='black', linewidth=0.5)
    b3 = ax.bar(x + 0.5 * width, clean_log(rel_dlk), width, label=r'DMDc-LQR-GMRES($m$)', color=color_dlk, edgecolor='black', linewidth=0.5)
    b4 = ax.bar(x + 1.5 * width, clean_log(rel_dlkl), width, label=r'DMDc-LQR-LGMRES($m,\ell$)', color=color_dlkl, edgecolor='black', linewidth=0.5)

    ax.set_ylabel('Relative Time (normalized to GMRES)', fontsize=10, fontweight='semibold')
    ax.set_xlabel('Problem Number', fontsize=10, fontweight='semibold')
    ax.set_yscale('log')
    ax.set_ylim([0.01, 50])

    ax.set_xticks(x)
    ax.set_xticklabels([str(p) for p in problems], fontsize=9)
    ax.grid(axis='y', which='both', linestyle=':', linewidth=0.5, alpha=0.6)
    ax.legend(frameon=True, fontsize=8.5, loc='upper right', ncol=2)
    ax.tick_params(direction='in', top=True, right=True, which='both')

    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(output_prefix)), exist_ok=True)
    png_path = f"{output_prefix}.png"
    eps_path = f"{output_prefix}.eps"
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(eps_path, format='eps', bbox_inches='tight')
    plt.close()
    print(f"[✓] Gráfico de tiempos relativos guardado:\n    {png_path}\n    {eps_path}")

def main():
    plot_number_of_cycles()
    plot_relative_time()

if __name__ == '__main__':
    main()
