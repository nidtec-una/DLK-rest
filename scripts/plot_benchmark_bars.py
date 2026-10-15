#!/usr/bin/env python3
"""
scripts/plot_benchmark_bars.py

Reproduce con exactitud la lógica visual y matemática de show_bar_dual_custom.m
de MATLAB para generar los gráficos de barras del artículo en LaTeX:
1. results/figures/bar_plot_Number_of_cycles.eps y .png
2. results/figures/bar_plot_Average_computation_relative_time.eps y .png

Lógica:
- Lee results/tables/cycles_experiment.csv y results/tables/avg_relativa_time_experiment.csv
- Calcula los ratios log10:
  * log10(GMRES / LGMRES)
  * log10(GMRES / DMDc-LQR-GMRES)
  * log10(GMRES / DMDc-LQR-LGMRES)
- Filtra las matrices que no convergieron en todos los métodos (valid_mask)
- Grafica 3 barras por problema activo en escala de grises [0.2, 0.5, 0.8]
"""

import os
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def parse_val(v):
    if v is None or v == "" or v == "NC" or str(v).strip().lower() == "nan":
        return np.nan
    try:
        return float(v)
    except ValueError:
        return np.nan

def generate_bar_dual_custom(csv_path, title_label, y_label, output_filename, is_cycle=True):
    if not os.path.exists(csv_path):
        print(f"[Error] No se encontró {csv_path}")
        return

    problems = []
    matrices = []
    col_gmres = []
    col_lgmres = []
    col_dlk = []

    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            problems.append(int(row["Problem"]))
            matrices.append(row["Matrix"])
            if is_cycle:
                col_gmres.append(parse_val(row.get("Cycles_GMRES")))
                col_lgmres.append(parse_val(row.get("Cycles_LGMRES")))
                col_dlk.append(parse_val(row.get("Cycles_GMRES_DMDc_LQR")))
            else:
                col_gmres.append(parse_val(row.get("avg_time_GMRES")))
                col_lgmres.append(parse_val(row.get("avg_time_LGMRES")))
                col_dlk.append(parse_val(row.get("avg_time_GMRES_DMDc_LQR")))

    n_total = len(problems)

    # 1. Calcular ratios log10(GMRES / Metodo)
    r_lgmres = []
    r_dlk = []

    for i in range(n_total):
        vg = col_gmres[i]
        vl = col_lgmres[i]
        vd = col_dlk[i]

        if not np.isnan(vg) and vg > 0 and not np.isnan(vl) and vl > 0:
            r_lgmres.append(np.log10(vg / vl))
        else:
            r_lgmres.append(np.nan)

        if not np.isnan(vg) and vg > 0 and not np.isnan(vd) and vd > 0:
            r_dlk.append(np.log10(vg / vd))
        else:
            r_dlk.append(np.nan)

    # 2. Máscara de validación (Filtrar NC/NaN para tener las 2 series válidas)
    valid_mask = []
    for i in range(n_total):
        is_valid = (not np.isnan(r_lgmres[i])) and (not np.isnan(r_dlk[i]))
        valid_mask.append(is_valid)

    active_indices = [i for i, v in enumerate(valid_mask) if v]
    active_ids = [problems[i] for i in active_indices]

    data_series = [
        [r_lgmres[i] for i in active_indices],
        [r_dlk[i] for i in active_indices]
    ]

    num_series = 2
    num_active = len(active_ids)
    if num_active == 0:
        print(f"[Aviso] No hay suficientes problemas convergentes para {title_label}")
        return

    # Proporciones SIAM exactas de show_bar_dual_custom.m: [2, 2, 16, 8.5] cm -> ~6.3 x 3.35 in
    fig, ax = plt.subplots(figsize=(6.3, 3.35))
    fig.patch.set_facecolor('white')

    # Paleta escala de grises de 2 series: [0.3 0.3 0.3; 0.7 0.7 0.7]
    colors = [(0.3, 0.3, 0.3), (0.7, 0.7, 0.7)]
    labels = [
        r'$\log(\mathrm{GMRES}(30)/\mathrm{LGMRES}(28,2))$',
        r'$\log(\mathrm{GMRES}(30)/\mathrm{DMDc\text{-}LQR\text{-}GMRES})$'
    ]

    bar_width = 0.8 / num_series
    x_pos = np.arange(num_active)

    handles = []
    for s in range(num_series):
        offset = (s - (num_series - 1) / 2.0) * bar_width
        vals = data_series[s]
        b = ax.bar(x_pos + offset, vals, width=bar_width,
                   color=colors[s], edgecolor=(0.15, 0.15, 0.15),
                   linewidth=0.5, alpha=0.9, label=labels[s])
        handles.append(b)

    # Línea base en y = 0
    ax.axhline(0, color='black', linewidth=0.6, linestyle='-')

    # Límites y padding idénticos a MATLAB
    all_vals = np.concatenate(data_series)
    max_val = np.max(all_vals)
    min_val = np.min(all_vals)
    range_val = max(1e-3, max_val - min_val)

    top_lim = max_val + 0.55 * range_val
    bottom_lim = min_val - 0.12 * range_val
    ax.set_ylim([bottom_lim, top_lim])
    ax.set_xlim([-0.6, num_active - 0.4])

    # Estilo de ejes SIAM
    ax.grid(True, which='both', axis='y', linestyle=':', linewidth=0.5, alpha=0.5)
    ax.set_axisbelow(True)
    ax.tick_params(direction='in', top=True, right=True, labelsize=8)

    ax.set_ylabel(y_label, fontsize=9)
    ax.set_xlabel('Problem', fontsize=9)
    ax.set_xticks(x_pos)
    ax.set_xticklabels([str(pid) for pid in active_ids], fontsize=8)

    # Leyenda sin borde
    ax.legend(handles=handles, labels=labels, fontsize=7.5, loc='upper right', frameon=False, ncol=1)

    plt.tight_layout()

    out_eps = f"{output_filename}.eps"
    out_png = f"{output_filename}.png"
    os.makedirs(os.path.dirname(os.path.abspath(out_eps)), exist_ok=True)

    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.savefig(out_eps, format='eps', bbox_inches='tight')
    plt.close()

    print(f"[✓] Gráfico generado exitosamente:")
    print(f"    PNG: {out_png}")
    print(f"    EPS: {out_eps}")

def main():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    os.chdir(repo_root)

    # 1. Gráfico de Ciclos (Number of cycles)
    generate_bar_dual_custom(
        csv_path="results/tables/cycles_experiment.csv",
        title_label="Number of cycles",
        y_label=r'$\log(\mathrm{GMRES}(30)/\mathrm{Krylov\ Solver})$',
        output_filename="results/figures/bar_plot_Number_of_cycles",
        is_cycle=True
    )

    # 2. Gráfico de Tiempo Relativo (Average computation relative time)
    generate_bar_dual_custom(
        csv_path="results/tables/avg_relativa_time_experiment.csv",
        title_label="Average computation relative time",
        y_label=r'$\log(\mathrm{GMRES}(30)/\mathrm{Krylov\ Solver})$',
        output_filename="results/figures/bar_plot_Average_computation_relative_time",
        is_cycle=False
    )

if __name__ == '__main__':
    main()
