#!/usr/bin/env python3
"""
plot_performance_with_m.py

Genera figuras estilo SIAM en 3 paneles:
1. Norma relativa del residuo ||r_k||/||r_0|| vs Ciclos.
2. Error de reconstrucción DMDc e = ||S_hat - S||_F / ||S||_F (barras con gradiente).
3. Dinámica del parámetro de reinicio m vs Ciclos.

Exporta tanto en formato PNG como en EPS vectorial.
"""

import argparse
import os
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

def plot_performance(csv_file, output_prefix, title_label="sherman5", m_max=331, m_min=2):
    cycles = []
    m_vals = []
    rel_res = []
    dmdc_errs = []
    
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            cycles.append(int(row['cycle']))
            m_vals.append(int(row['m']))
            rel_res.append(float(row['rel_residual']))
            dmdc_errs.append(float(row['dmdc_error']))
            
    cycles = np.array(cycles)
    m_vals = np.array(m_vals)
    rel_res = np.array(rel_res)
    dmdc_errs = np.array(dmdc_errs)
    
    num_cycles = len(cycles)
    x_end = cycles[-1]
    
    fig, axes = plt.subplots(3, 1, figsize=(7.0, 7.5), sharex=True,
                             gridspec_kw={'height_ratios': [1.2, 1.0, 1.0]})
    fig.patch.set_facecolor('white')
    
    # -------------------------------------------------------------
    # (1) Curva de Residuos
    # -------------------------------------------------------------
    ax1 = axes[0]
    ax1.semilogy(cycles, rel_res, color='#262626', linestyle='-', linewidth=1.2,
                 marker='^', markersize=4, markerfacecolor='none', markeredgecolor='#262626',
                 markevery=max(1, len(cycles) // 12), label=r'DMDc-LQR')
    ax1.set_ylabel(r'$\|r_k\|/\|r_0\|$', fontsize=11)
    ax1.set_ylim([1e-10, 1.5])
    ax1.set_xlim([0, x_end])
    ax1.grid(True, which='both', linestyle=':', linewidth=0.6, alpha=0.7)
    ax1.tick_params(labelsize=10, direction='in', top=True, right=True)
    ax1.set_title(f'Performance of DMDc-LQR ({title_label})', fontsize=11, fontweight='bold', pad=8)
    
    # -------------------------------------------------------------
    # (2) Barras de Error DMDc (con mapa de color en escala de grises SIAM)
    # -------------------------------------------------------------
    ax2 = axes[1]
    
    # Identificar snapshots donde dmdc_error > 0
    snap_indices = np.where(dmdc_errs > 0)[0]
    
    if len(snap_indices) > 0:
        valid_errs = dmdc_errs[snap_indices]
        min_err = np.min(valid_errs)
        max_err = np.max(valid_errs)
        if min_err == max_err:
            min_err = max_err * 0.1
        
        norm = mcolors.LogNorm(vmin=min_err * 0.8, vmax=max_err * 1.2)
        try:
            cmap = matplotlib.colormaps['gist_yarg']
        except AttributeError:
            cmap = cm.get_cmap('gist_yarg')
        
        # Intervalos de snapshot (cada barra cubre desde el snapshot previo hasta el actual)
        prev_idx = 0
        for s_idx in snap_indices:
            val = dmdc_errs[s_idx]
            x1 = prev_idx
            x2 = cycles[s_idx]
            width = max(1.0, float(x2 - x1))
            color = cmap(norm(val))
            
            ax2.bar(x1, val, width=width, align='edge',
                    color=color, edgecolor='#404040', linewidth=0.6, alpha=0.85)
            prev_idx = x2
            
        ax2.set_yscale('log')
        ax2.set_ylim([min_err * 0.5, max_err * 2.0])
        
        # Colorbar
        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax2, orientation='vertical', pad=0.02, fraction=0.046)
        cbar.ax.tick_params(labelsize=8)
        cbar.set_label(r'$e$', fontsize=10, rotation=0, labelpad=8)
    else:
        ax2.text(0.5, 0.5, 'No DMDc snapshots recorded', ha='center', va='center', transform=ax2.transAxes)
    
    ax2.set_ylabel(r'$e = \|\hat{S} - S\|_F / \|S\|_F$', fontsize=10)
    ax2.grid(True, which='both', linestyle=':', linewidth=0.6, alpha=0.7)
    ax2.tick_params(labelsize=10, direction='in', top=True, right=True)
    
    # -------------------------------------------------------------
    # (3) Valores dinámicos de m
    # -------------------------------------------------------------
    ax3 = axes[2]
    # Stairs plot de m vs ciclos
    ax3.step(cycles, m_vals, where='post', color='#1a1a1a', linewidth=1.3, label=r'$m_k$')
    
    # Cotas m_min y m_max
    if m_max is not None:
        ax3.axhline(m_max, color='#666666', linestyle='--', linewidth=1.0,
                    label=rf'$m_{{\max}} = {m_max}$')
        ax3.text(x_end * 0.02, m_max * 1.03, rf'$m_{{\max}} = {m_max}$',
                 fontsize=9, color='#404040', verticalalignment='bottom')
                 
    if m_min is not None and m_min > 0:
        ax3.axhline(m_min, color='#888888', linestyle=':', linewidth=0.8)
    
    y_upper = max(m_max * 1.15 if m_max else 350, float(np.max(m_vals)) * 1.15)
    ax3.set_ylim([0, y_upper])
    ax3.set_xlabel('Number of cycles ($k$)', fontsize=11)
    ax3.set_ylabel(r'$m$', fontsize=11)
    ax3.grid(True, which='both', linestyle=':', linewidth=0.6, alpha=0.7)
    ax3.tick_params(labelsize=10, direction='in', top=True, right=True)
    
    plt.tight_layout()
    
    # Asegurar directorios de salida
    png_path = f"{output_prefix}.png"
    eps_path = f"{output_prefix}.eps"
    os.makedirs(os.path.dirname(os.path.abspath(png_path)), exist_ok=True)
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(eps_path, format='eps', bbox_inches='tight')
    plt.close()
    
    print(f"[OK] Figuras guardadas exitosamente:")
    print(f"     PNG: {png_path}")
    print(f"     EPS: {eps_path}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Graficar curvas de rendimiento DLK con dinámica de m")
    parser.add_argument('--cycle_file', type=str, default="results/tables/cycles_orsirr1_dlk.csv", help="Archivo CSV con el historial de ciclos")
    parser.add_argument('--output_prefix', type=str, default="results/figures/PerformanceofDMDcLQRGMRESmorsirr1_with_m", help="Prefijo para archivos PNG y EPS")
    parser.add_argument('--title', type=str, default="orsirr1", help="Nombre o etiqueta de la matriz")
    parser.add_argument('--m_max', type=int, default=103, help="Cota máxima m_max")
    parser.add_argument('--m_min', type=int, default=2, help="Cota mínima m_min")
    args = parser.parse_args()
    
    cycle_file = args.cycle_file
    if not os.path.exists(cycle_file):
        cand = os.path.join("results/tables", os.path.basename(cycle_file))
        if os.path.exists(cand):
            cycle_file = cand
            
    plot_performance(cycle_file, args.output_prefix, args.title, args.m_max, args.m_min)
