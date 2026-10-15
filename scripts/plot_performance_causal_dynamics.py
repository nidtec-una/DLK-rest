#!/usr/bin/env python3
"""
scripts/plot_performance_causal_dynamics.py

Genera una figura explicativa y didáctica para la sección C.2.1 que ilustra de forma explícita:
1. El mecanismo de muestreo por ventanas (snapshots) en el subgráfico 1.
2. La identificación de modelos reducidos DMDc y su error de reconstrucción en el subgráfico 2.
3. La resolución del control óptimo LQR (CARE) en cada frontera k = p, 2p, ... y la actualización de m.
4. El delay causal inherente de una ventana (Delta k = p ciclos) entre los datos observados
   con m_{i-1} y la aplicación de la acción de control m_i en la etapa siguiente [k_i, k_{i+1}].
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
from matplotlib.patches import ConnectionPatch

def plot_causal_dynamics(csv_file, output_prefix, title_label="orsirr1", m_max=103, m_min=2):
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
    
    x_end = cycles[-1]
    
    # -------------------------------------------------------------
    # Extracción de etapas de control
    # -------------------------------------------------------------
    snap_indices = np.where(dmdc_errs > 0)[0]
    
    stages = []
    prev_cycle = 0
    
    for s_idx in snap_indices:
        s_cycle = cycles[s_idx]
        val = dmdc_errs[s_idx]
        
        if prev_cycle == 0:
            m_active = m_vals[0]
        else:
            idx = np.where(cycles == prev_cycle + 1)[0]
            if len(idx) > 0:
                m_active = m_vals[idx[0]]
            else:
                m_active = m_vals[s_idx]
                
        stages.append({
            'x_start': prev_cycle,
            'x_end': s_cycle,
            'm': m_active,
            'error': val
        })
        prev_cycle = s_cycle
        
    if prev_cycle < cycles[-1]:
        idx = np.where(cycles == prev_cycle + 1)[0]
        if len(idx) > 0:
            m_active = m_vals[idx[0]]
        else:
            m_active = m_vals[-1]
            
        stages.append({
            'x_start': prev_cycle,
            'x_end': cycles[-1],
            'm': m_active,
            'error': None
        })

    # Proporciones óptimas para visualización causal
    fig, axes = plt.subplots(3, 1, figsize=(8.0, 9.4), sharex=True,
                             gridspec_kw={'height_ratios': [1.25, 1.05, 1.15]})
    fig.patch.set_facecolor('white')
    
    ax1 = axes[0]
    ax2 = axes[1]
    ax3 = axes[2]
    
    # Colores temáticos suaves para las ventanas iniciales
    c_w1 = '#eaf2f8'  # Azul suave (Ventana 1)
    c_w2 = '#fef5e7'  # Naranja suave (Ventana 2)
    c_w3 = '#eafaf1'  # Verde menta suave (Ventana 3)
    c_w4 = '#f4ecf7'  # Púrpura suave (Ventana 4)
    
    # -------------------------------------------------------------
    # (1) Curva de Residuos con Sombreado de Ventanas de Muestreo
    # -------------------------------------------------------------
    ax1.axvspan(0, 5, facecolor=c_w1, alpha=0.85, zorder=0)
    ax1.axvspan(5, 10, facecolor=c_w2, alpha=0.85, zorder=0)
    ax1.axvspan(10, 15, facecolor=c_w3, alpha=0.85, zorder=0)
    ax1.axvspan(15, 20, facecolor=c_w4, alpha=0.85, zorder=0)
    
    ax1.semilogy(cycles, rel_res, color='#1f2d3d', linestyle='-', linewidth=1.3,
                 marker='^', markersize=4.5, markerfacecolor='white', markeredgecolor='#1f2d3d',
                 markevery=1, zorder=3, label=r'Relative Residual $\|r_k\|/\|r_0\|$')
                 
    ax1.set_ylabel(r'$\|r_k\|/\|r_0\|$', fontsize=11)
    ax1.set_ylim([1e-10, 1.8])
    ax1.set_title(f'Causal Feedback Dynamics and Batch Sampling Delay in DMDc-LQR ({title_label})',
                  fontsize=11.5, fontweight='bold', pad=10)
    
    # Batch annotations in Subplot 1
    ax1.text(2.5, 0.08, 'Batch 1\n($p=5$, with $m_0$)', fontsize=8, ha='center', va='center',
             bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='#7094b7', alpha=0.92))
    ax1.text(7.5, 0.008, 'Batch 2\n($p=5$, with $m_1$)', fontsize=8, ha='center', va='center',
             bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='#d9a779', alpha=0.92))
    ax1.text(12.5, 2e-5, 'Batch 3\n(with $m_2$)', fontsize=8, ha='center', va='center',
             bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='#62aa7e', alpha=0.92))
    ax1.text(17.5, 4e-8, 'Batch 4\n(with $m_3$)', fontsize=8, ha='center', va='center',
             bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='#9973b5', alpha=0.92))
    
    # -------------------------------------------------------------
    # (2) Barras de Error DMDc con Identificación de Modelos M_i
    # -------------------------------------------------------------
    valid_errs = np.array([st['error'] for st in stages if st['error'] is not None and st['error'] > 0])
    min_err = np.min(valid_errs)
    max_err = np.max(valid_errs)
    
    y_min = 10.0 ** (np.floor(np.log10(min_err)) - 0.5)
    y_max = 10.0 ** (np.ceil(np.log10(max_err)) + 0.8)
    
    norm = mcolors.LogNorm(vmin=min_err, vmax=max_err)
    cdict = {
        'red':   [(0.0, 0.82, 0.82), (1.0, 0.20, 0.20)],
        'green': [(0.0, 0.82, 0.82), (1.0, 0.20, 0.20)],
        'blue':  [(0.0, 0.82, 0.82), (1.0, 0.20, 0.20)]
    }
    cmap = mcolors.LinearSegmentedColormap('siam_gray', cdict, 256)
    
    ax2.axvspan(0, 5, facecolor=c_w1, alpha=0.5, zorder=0)
    ax2.axvspan(5, 10, facecolor=c_w2, alpha=0.5, zorder=0)
    ax2.axvspan(10, 15, facecolor=c_w3, alpha=0.5, zorder=0)
    ax2.axvspan(15, 20, facecolor=c_w4, alpha=0.5, zorder=0)
    
    for st in stages:
        if st['error'] is not None and st['error'] > 0:
            x1 = st['x_start']
            x2 = st['x_end']
            val = st['error']
            color = cmap(norm(val))
            ax2.fill([x1, x2, x2, x1], [y_min, y_min, val, val],
                     facecolor=color, edgecolor='#404040', linewidth=0.7, alpha=0.9, zorder=2)
                     
    ax2.set_yscale('log')
    ax2.set_ylim([y_min, y_max])
    ax2.set_ylabel(r'$e = \|\hat{S} - S\|_F / \|S\|_F$', fontsize=10)
    
    # Model identification labels
    ax2.text(2.5, stages[0]['error'] * 0.35, r'$\mathcal{M}_1$', fontsize=9, fontweight='bold',
             ha='center', va='center', color='#222222', zorder=4)
    ax2.text(7.5, stages[1]['error'] * 0.35, r'$\mathcal{M}_2$', fontsize=9, fontweight='bold',
             ha='center', va='center', color='#222222', zorder=4)
    ax2.text(12.5, stages[2]['error'] * 0.35, r'$\mathcal{M}_3$', fontsize=9, fontweight='bold',
             ha='center', va='center', color='#222222', zorder=4)
    ax2.text(17.5, stages[3]['error'] * 0.35, r'$\mathcal{M}_4$', fontsize=9, fontweight='bold',
             ha='center', va='center', color='#ffffff', zorder=4)
             
    # LQR decision nodes at k = 5 and 10
    ax2.plot(5, stages[0]['error'], marker='o', markersize=6.5, color='#b22222', zorder=5)
    ax2.plot(10, stages[1]['error'], marker='o', markersize=6.5, color='#b22222', zorder=5)
    
    ax2.annotate('LQR 1 at $k=5$\n' + r'$\mathcal{M}_1 \to m_1=42$', xy=(5, stages[0]['error']),
                 xytext=(2.2, 4e-11),
                 arrowprops=dict(arrowstyle='->', color='#b22222', lw=1.1, connectionstyle="arc3,rad=-0.1"),
                 fontsize=8, color='#b22222', fontweight='bold',
                 bbox=dict(boxstyle='round,pad=0.25', facecolor='#fffbfa', edgecolor='#b22222', alpha=0.95),
                 zorder=6)
                 
    ax2.annotate('LQR 2 at $k=10$\n' + r'$\mathcal{M}_2 \to m_2=90$', xy=(10, stages[1]['error']),
                 xytext=(7.2, 1.3e-10),
                 arrowprops=dict(arrowstyle='->', color='#b22222', lw=1.1, connectionstyle="arc3,rad=-0.1"),
                 fontsize=8, color='#b22222', fontweight='bold',
                 bbox=dict(boxstyle='round,pad=0.25', facecolor='#fffbfa', edgecolor='#b22222', alpha=0.95),
                 zorder=6)

    # -------------------------------------------------------------
    # (3) Valores dinámicos de m y Representación del Delay
    # -------------------------------------------------------------
    ax3.axvspan(0, 5, facecolor=c_w1, alpha=0.5, zorder=0)
    ax3.axvspan(5, 10, facecolor=c_w2, alpha=0.5, zorder=0)
    ax3.axvspan(10, 15, facecolor=c_w3, alpha=0.5, zorder=0)
    ax3.axvspan(15, 20, facecolor=c_w4, alpha=0.5, zorder=0)
    
    x_stairs = []
    y_stairs = []
    for st in stages:
        x_stairs.extend([st['x_start'], st['x_end']])
        y_stairs.extend([st['m'], st['m']])
        
    ax3.plot(x_stairs, y_stairs, color='#111111', linewidth=1.5, zorder=3, label=r'$m_k$')
    
    # Lines for m bounds
    if m_max is not None:
        ax3.axhline(m_max, color='#666666', linestyle='--', linewidth=1.0, zorder=2)
        ax3.text(x_end * 0.02, m_max * 1.03, rf'$m_{{\max}} = {m_max}$',
                 fontsize=8.5, color='#404040', verticalalignment='bottom')
    if m_min is not None and m_min > 0:
        ax3.axhline(m_min, color='#888888', linestyle=':', linewidth=0.8, zorder=2)
        
    y_upper = max(m_max * 1.22 if m_max else 350, float(np.max(m_vals)) * 1.22)
    ax3.set_ylim([0, y_upper])
    ax3.set_xlabel('Number of cycles ($k$)', fontsize=11)
    ax3.set_ylabel(r'$m$', fontsize=11)
    
    # Control action annotations in Subplot 3
    ax3.text(2.5, 16, r'$m_0 = 30$' + '\n(initial)', fontsize=8, ha='center', va='center',
             bbox=dict(boxstyle='square,pad=0.2', facecolor='white', edgecolor='#888888', alpha=0.92))
    ax3.text(7.5, 25, r'$m_1 = 42$' + '\n' + r'(from $\mathcal{M}_1$)', fontsize=8, ha='center', va='center',
             bbox=dict(boxstyle='square,pad=0.2', facecolor='#fffbfa', edgecolor='#b22222', alpha=0.92))
    ax3.text(12.5, 55, r'$m_2 = 90$' + '\n' + r'(from $\mathcal{M}_2$)', fontsize=8, ha='center', va='center',
             bbox=dict(boxstyle='square,pad=0.2', facecolor='#fffbfa', edgecolor='#b22222', alpha=0.92))
    ax3.text(17.5, 78, r'$m_3 = 103$' + '\n' + r'($m_{\max}$)', fontsize=8, ha='center', va='center',
             bbox=dict(boxstyle='square,pad=0.2', facecolor='white', edgecolor='#333333', alpha=0.92))

    # Batch Delay Indicator in Subplot 3
    ax3.annotate('', xy=(5, 55), xytext=(0, 55),
                 arrowprops=dict(arrowstyle='<->', color='#1d568c', lw=1.3))
    ax3.text(2.5, 58, 'Batch delay: ' + r'$\Delta k = p = 5$', fontsize=7.5,
             ha='center', va='bottom', color='#1d568c', fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.15', facecolor='#edf4fc', edgecolor='#1d568c', alpha=0.92))

    # Líneas verticales guía cruzando los 3 paneles en las fronteras de decisión k = 5, 10, 15, 20
    for split_k in [5, 10, 15, 20]:
        for ax in axes:
            ax.axvline(split_k, color='#b22222', linestyle='--', linewidth=0.85, alpha=0.75, zorder=1)

    # Configuración uniforme de límites horizontales y rejilla en los 3 ejes
    for ax in axes:
        ax.set_xlim([0, x_end])
        ax.grid(True, which='both', linestyle=':', linewidth=0.6, alpha=0.6)

    ax1.tick_params(labelsize=10, direction='in', top=True, right=True, labelbottom=False)
    ax2.tick_params(labelsize=10, direction='in', top=True, right=True, labelbottom=False)
    ax3.tick_params(labelsize=10, direction='in', top=True, right=True, labelbottom=True)

    # Ajuste de geometría fija
    fig.subplots_adjust(left=0.13, right=0.80, top=0.94, bottom=0.07, hspace=0.18)
    fig.align_ylabels(axes)
    
    # Colorbar desplazado a la derecha
    pos2 = ax2.get_position()
    cax = fig.add_axes([0.855, pos2.y0, 0.022, pos2.height])
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, orientation='vertical')
    cbar.ax.tick_params(labelsize=8.5)
    cbar.set_label(r'$e$', fontsize=10, rotation=0, labelpad=10)

    # Conectores causales entre subgráficos (Flechas de realimentación en lazo cerrado)
    # Flecha 1: Del fin del modelo 1 (ax2, k=5) al escalón de m_1 (ax3, k=5)
    con1 = ConnectionPatch(xyA=(5, stages[0]['error']), coordsA='data', axesA=ax2,
                           xyB=(5, 42), coordsB='data', axesB=ax3,
                           arrowstyle="-|>", color='#b22222', lw=1.6, mutation_scale=13,
                           connectionstyle="arc3,rad=-0.12", zorder=10)
    fig.add_artist(con1)
    
    # Flecha 2: Del fin del modelo 2 (ax2, k=10) al escalón de m_2 (ax3, k=10)
    con2 = ConnectionPatch(xyA=(10, stages[1]['error']), coordsA='data', axesA=ax2,
                           xyB=(10, 90), coordsB='data', axesB=ax3,
                           arrowstyle="-|>", color='#b22222', lw=1.6, mutation_scale=13,
                           connectionstyle="arc3,rad=-0.12", zorder=10)
    fig.add_artist(con2)

    # Guardar figuras
    png_path = f"{output_prefix}.png"
    eps_path = f"{output_prefix}.eps"
    os.makedirs(os.path.dirname(os.path.abspath(png_path)), exist_ok=True)
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(eps_path, format='eps', bbox_inches='tight')
    plt.close()
    
    print(f"[OK] Gráfico causal generado exitosamente:")
    print(f"     PNG: {png_path}")
    print(f"     EPS: {eps_path}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Graficar dinámica causal y delay en DMDc-LQR")
    parser.add_argument('--cycle_file', type=str, default="results/tables/cycles_orsirr1_dlk.csv")
    parser.add_argument('--output_prefix', type=str, default="results/figures/Performance_causal_dynamics_orsirr1")
    parser.add_argument('--title', type=str, default="orsirr1")
    parser.add_argument('--m_max', type=int, default=103)
    parser.add_argument('--m_min', type=int, default=2)
    args = parser.parse_args()
    
    cycle_file = args.cycle_file
    if not os.path.exists(cycle_file):
        cand = os.path.join("results/tables", os.path.basename(cycle_file))
        if os.path.exists(cand):
            cycle_file = cand
            
    plot_causal_dynamics(cycle_file, args.output_prefix, args.title, args.m_max, args.m_min)
