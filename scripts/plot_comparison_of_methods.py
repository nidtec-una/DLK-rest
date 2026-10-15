#!/usr/bin/env python3
"""
plot_comparison_of_methods.py

Genera la gráfica de comparación de métodos de Krylov en formato estándar SIAM
para publicación (replicando src/+plots/show_results_custom.m de MATLAB):
- GMRES(110)
- LGMRES(107, 3)
- DMDc-LQR-GMRES
- DMDc-LQR-LGMRES

Exporta en formato PNG de alta resolución (300 DPI) y EPS vectorial para LaTeX.
"""

import argparse
import os
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load_cycle_data(csv_file):
    cycles = []
    rel_res = []
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            cycles.append(int(row['cycle']))
            rel_res.append(float(row['rel_residual']))
    return np.array(cycles), np.array(rel_res)

def find_csv(filename):
    candidates = [
        os.path.join("results/tables", filename),
        os.path.join("results", filename),
        filename
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    return os.path.join("results/tables", filename)

def plot_comparison(output_prefix="results/figures/Comparisonofmethodssherman5",
                    title_label="sherman5",
                    tol=1e-9):
    # Configuración de tamaño de figura SIAM: ~13 x 10 cm (5.12 x 3.94 pulgadas)
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    fig.patch.set_facecolor('white')

    series = [
        {
            'file': find_csv('cycles_sherman5_gmres110.csv'),
            'name': 'GMRES(110)',
            'color': '#808080',
            'linestyle': ':',
            'marker': 'o',
            'linewidth': 1.1,
            'markersize': 5
        },
        {
            'file': find_csv('cycles_sherman5_lgmres107.csv'),
            'name': 'LGMRES(107,3)',
            'color': '#404040',
            'linestyle': '--',
            'marker': 's',
            'linewidth': 1.2,
            'markersize': 4.5
        },
        {
            'file': find_csv('cycles_sherman5_dlk.csv'),
            'name': 'DMDc-LQR-GMRES',
            'color': '#000000',
            'linestyle': '--',
            'marker': '^',
            'linewidth': 1.4,
            'markersize': 5
        },
        {
            'file': find_csv('cycles_sherman5_dlklgmres.csv'),
            'name': 'DMDc-LQR-LGMRES',
            'color': '#1a1a1a',
            'linestyle': '-.',
            'marker': 'd',
            'linewidth': 1.3,
            'markersize': 4.5
        }
    ]

    max_x = 0
    for s in series:
        if not os.path.exists(s['file']):
            print(f"[Aviso] No se encontró {s['file']}, omitiendo curva.")
            continue
        c, r = load_cycle_data(s['file'])
        total_cycles = len(c) - 1
        if total_cycles > max_x:
            max_x = total_cycles

        # Rango de marcadores limpio (10 marcadores por curva)
        markevery = max(1, len(c) // 10)
        label_text = f"{s['name']}, cycles = {total_cycles}"

        ax.semilogy(c, r,
                    linestyle=s['linestyle'],
                    color=s['color'],
                    marker=s['marker'],
                    markersize=s['markersize'],
                    markerfacecolor='none',
                    markeredgecolor=s['color'],
                    markeredgewidth=1.1,
                    linewidth=s['linewidth'],
                    markevery=markevery,
                    label=label_text)

    ax.set_xlabel('Number of cycles ($k$)', fontsize=11)
    ax.set_ylabel(r'Residual norm $\|r_k\|$', fontsize=11)
    ax.set_ylim([tol * 0.5, 1.3])
    ax.set_xlim([0, max_x + 5])

    # Estilo de leyenda SIAM sin marco
    ax.legend(frameon=False, fontsize=9.5, loc='upper right')

    # Grid punteado y ticks hacia adentro
    ax.grid(True, which='both', linestyle=':', linewidth=0.6, alpha=0.6)
    ax.tick_params(labelsize=10, direction='in', top=True, right=True, which='both')

    plt.tight_layout()

    # Guardar en PNG y EPS
    png_path = f"{output_prefix}.png"
    eps_path = f"{output_prefix}.eps"
    os.makedirs(os.path.dirname(os.path.abspath(png_path)), exist_ok=True)

    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(eps_path, format='eps', bbox_inches='tight')
    plt.close()

    print(f"[OK] Gráfica comparativa guardada exitosamente:")
    print(f"     PNG: {png_path}")
    print(f"     EPS: {eps_path}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Graficar comparación de métodos Krylov estilo SIAM")
    parser.add_argument('--output_prefix', type=str, default="results/figures/Comparisonofmethodssherman5",
                        help="Prefijo para archivos PNG y EPS")
    parser.add_argument('--title', type=str, default="sherman5", help="Etiqueta del problema")
    parser.add_argument('--tol', type=float, default=1e-9, help="Tolerancia del residuo")
    args = parser.parse_args()

    plot_comparison(args.output_prefix, args.title, args.tol)
