#!/usr/bin/env python3
"""
scripts/plot_performance_causal_dynamics.py

Generates explanatory causal dynamics figures for sections C.1.1 (young3c), C.2.1 (orsirr1),
and C.3.1 (wang2) illustrating:
1. Batch sampling mechanism ($p=5$) in Subplot 1.
2. DMDc reduced model identification and reconstruction error bars in Subplot 2.
3. Optimal control (CARE/LQR) decisions at batch boundaries $k=p, 2p, \\dots$.
4. Inherent causal delay (\\Delta k = p cycles) between observations and control actions.
5. Closed-loop feedback arrows connecting model identification to control action steps.
All annotations and terminology in English using "Batch" and "Batch delay".
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

def plot_causal_dynamics_matrix(matrix_name, csv_file, output_prefix):
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
    
    # Colores temáticos suaves para los primeros 4 lotes
    c_w1 = '#eaf2f8'  # Azul suave (Batch 1)
    c_w2 = '#fef5e7'  # Naranja suave (Batch 2)
    c_w3 = '#eafaf1'  # Verde menta suave (Batch 3)
    c_w4 = '#f4ecf7'  # Púrpura suave (Batch 4)
    
    # -------------------------------------------------------------
    # (1) Curva de Residuos con Sombreado de Lotes de Muestreo
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
    ax1.set_title(f'Causal Feedback Dynamics and Batch Sampling Delay in DMDc-LQR ({matrix_name})',
                  fontsize=11.5, fontweight='bold', pad=10)
    
    # Anotaciones de Batch en Panel 1 según matriz
    if matrix_name == 'orsirr1':
        ax1.text(2.5, 0.08, 'Batch 1\n($p=5$, with $m_0$)', fontsize=8, ha='center', va='center',
                 bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='#7094b7', alpha=0.92))
        ax1.text(7.5, 0.008, 'Batch 2\n($p=5$, with $m_1$)', fontsize=8, ha='center', va='center',
                 bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='#d9a779', alpha=0.92))
        ax1.text(12.5, 2e-5, 'Batch 3\n(with $m_2$)', fontsize=8, ha='center', va='center',
                 bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='#62aa7e', alpha=0.92))
        ax1.text(17.5, 4e-8, 'Batch 4\n(with $m_3$)', fontsize=8, ha='center', va='center',
                 bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='#9973b5', alpha=0.92))
    elif matrix_name == 'young3c':
        ax1.text(2.5, 0.20, 'Batch 1\n($p=5$)', fontsize=7.5, ha='center', va='center',
                 bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#7094b7', alpha=0.92))
        ax1.text(7.5, 2.0e-4, 'Batch 2\n($p=5$)', fontsize=7.5, ha='center', va='center',
                 bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#d9a779', alpha=0.92))
        ax1.text(12.5, 0.015, 'Batch 3\n(with $m_2$)', fontsize=7.5, ha='center', va='center',
                 bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#62aa7e', alpha=0.92))
        ax1.text(17.5, 2.0e-5, 'Batch 4\n(with $m_3$)', fontsize=7.5, ha='center', va='center',
                 bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#9973b5', alpha=0.92))
    elif matrix_name == 'wang2':
        ax1.text(2.5, 0.04, 'Batch 1\n($p=5$, with $m_0$)', fontsize=8, ha='center', va='center',
                 bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='#7094b7', alpha=0.92))
        ax1.text(7.5, 0.003, 'Batch 2\n($p=5$, with $m_1$)', fontsize=8, ha='center', va='center',
                 bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='#d9a779', alpha=0.92))
        ax1.text(12.5, 0.002, 'Batch 3\n(with $m_2$)', fontsize=8, ha='center', va='center',
                 bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='#62aa7e', alpha=0.92))
        ax1.text(17.5, 1.2e-8, 'Batch 4\n(with $m_3$)', fontsize=8, ha='center', va='center',
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
    
    # Model identification labels (primeros 4 modelos)
    ax2.text(2.5, stages[0]['error'] * 0.35, r'$\mathcal{M}_1$', fontsize=9, fontweight='bold',
             ha='center', va='center', color='#222222', zorder=4)
    ax2.text(7.5, stages[1]['error'] * 0.35, r'$\mathcal{M}_2$', fontsize=9, fontweight='bold',
             ha='center', va='center', color='#222222', zorder=4)
    ax2.text(12.5, stages[2]['error'] * 0.35, r'$\mathcal{M}_3$', fontsize=9, fontweight='bold',
             ha='center', va='center', color='#222222', zorder=4)
    ax2.text(17.5, stages[3]['error'] * 0.35, r'$\mathcal{M}_4$', fontsize=9, fontweight='bold',
             ha='center', va='center', color='#ffffff' if norm(stages[3]['error']) > 0.6 else '#222222', zorder=4)
             
    # LQR decision nodes and annotations according to matrix
    ax2.plot(5, stages[0]['error'], marker='o', markersize=6.5, color='#b22222', zorder=5)
    ax2.plot(10, stages[1]['error'], marker='o', markersize=6.5, color='#b22222', zorder=5)
    
    if matrix_name == 'orsirr1':
        m1_target = 42
        m2_target = 90
        ax2.annotate(f'LQR 1 at $k=5$\n' + rf'$\mathcal{{M}}_1 \to m_1={m1_target}$', xy=(5, stages[0]['error']),
                     xytext=(2.2, 4e-11),
                     arrowprops=dict(arrowstyle='->', color='#b22222', lw=1.1, connectionstyle="arc3,rad=-0.1"),
                     fontsize=8, color='#b22222', fontweight='bold',
                     bbox=dict(boxstyle='round,pad=0.25', facecolor='#fffbfa', edgecolor='#b22222', alpha=0.95),
                     zorder=6)
        ax2.annotate(f'LQR 2 at $k=10$\n' + rf'$\mathcal{{M}}_2 \to m_2={m2_target}$', xy=(10, stages[1]['error']),
                     xytext=(7.2, 1.3e-10),
                     arrowprops=dict(arrowstyle='->', color='#b22222', lw=1.1, connectionstyle="arc3,rad=-0.1"),
                     fontsize=8, color='#b22222', fontweight='bold',
                     bbox=dict(boxstyle='round,pad=0.25', facecolor='#fffbfa', edgecolor='#b22222', alpha=0.95),
                     zorder=6)
    elif matrix_name == 'young3c':
        m1_target = 76
        m2_target = 51
        ax2.annotate(f'LQR 1 at $k=5$\n' + rf'$\mathcal{{M}}_1 \to m_1={m1_target}$', xy=(5, stages[0]['error']),
                     xytext=(1.8, 5e-11),
                     arrowprops=dict(arrowstyle='->', color='#b22222', lw=1.1, connectionstyle="arc3,rad=-0.1"),
                     fontsize=7.5, color='#b22222', fontweight='bold',
                     bbox=dict(boxstyle='round,pad=0.2', facecolor='#fffbfa', edgecolor='#b22222', alpha=0.95),
                     zorder=6)
        ax2.annotate(f'LQR 2 at $k=10$\n' + rf'$\mathcal{{M}}_2 \to m_2={m2_target}$', xy=(10, stages[1]['error']),
                     xytext=(13.0, 1.4e-10),
                     arrowprops=dict(arrowstyle='->', color='#b22222', lw=1.1, connectionstyle="arc3,rad=-0.1"),
                     fontsize=7.5, color='#b22222', fontweight='bold',
                     bbox=dict(boxstyle='round,pad=0.2', facecolor='#fffbfa', edgecolor='#b22222', alpha=0.95),
                     zorder=6)
    elif matrix_name == 'wang2':
        m1_target = 57
        m2_target = 12
        ax2.annotate(f'LQR 1 at $k=5$\n' + rf'$\mathcal{{M}}_1 \to m_1={m1_target}$', xy=(5, stages[0]['error']),
                     xytext=(2.3, 4e-11),
                     arrowprops=dict(arrowstyle='->', color='#b22222', lw=1.1, connectionstyle="arc3,rad=-0.1"),
                     fontsize=8, color='#b22222', fontweight='bold',
                     bbox=dict(boxstyle='round,pad=0.25', facecolor='#fffbfa', edgecolor='#b22222', alpha=0.95),
                     zorder=6)
        ax2.annotate(f'LQR 2 at $k=10$\n' + rf'$\mathcal{{M}}_2 \to m_2={m2_target}$', xy=(10, stages[1]['error']),
                     xytext=(7.5, 1.2e-10),
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
    
    # Parámetros m_max y m_min según matriz
    matrix_params = {
        'orsirr1': {'m_max': 103, 'm_min': 2, 'm_pos_x': 0.85},
        'young3c': {'m_max': 84,  'm_min': 2, 'm_pos_x': 0.85},
        'wang2':   {'m_max': 290, 'm_min': 2, 'm_pos_x': 0.85}
    }
    m_max = matrix_params[matrix_name]['m_max']
    m_min = matrix_params[matrix_name]['m_min']
    m_pos_x = matrix_params[matrix_name]['m_pos_x']
    
    if m_max is not None:
        ax3.axhline(m_max, color='#666666', linestyle='--', linewidth=1.0, zorder=2)
        ax3.text(x_end * m_pos_x, m_max * 1.03, rf'$m_{{\max}} = {m_max}$',
                 fontsize=8.5, color='#404040', ha='center', verticalalignment='bottom')
    if m_min is not None and m_min > 0:
        ax3.axhline(m_min, color='#888888', linestyle=':', linewidth=0.8, zorder=2)
        
    y_upper = max(m_max * 1.22 if m_max else 350, float(np.max(m_vals)) * 1.22)
    ax3.set_ylim([0, y_upper])
    ax3.set_xlabel('Number of cycles ($k$)', fontsize=11)
    ax3.set_ylabel(r'$m$', fontsize=11)
    
    # Control action annotations in Subplot 3 & Delay arrows (sin cruzar borde izquierdo)
    if matrix_name == 'orsirr1':
        ax3.text(2.5, 16, r'$m_0 = 30$' + '\n(initial)', fontsize=8, ha='center', va='center',
                 bbox=dict(boxstyle='square,pad=0.2', facecolor='white', edgecolor='#888888', alpha=0.92))
        ax3.text(7.5, 25, r'$m_1 = 42$' + '\n' + r'(from $\mathcal{M}_1$)', fontsize=8, ha='center', va='center',
                 bbox=dict(boxstyle='square,pad=0.2', facecolor='#fffbfa', edgecolor='#b22222', alpha=0.92))
        ax3.text(12.5, 55, r'$m_2 = 90$' + '\n' + r'(from $\mathcal{M}_2$)', fontsize=8, ha='center', va='center',
                 bbox=dict(boxstyle='square,pad=0.2', facecolor='#fffbfa', edgecolor='#b22222', alpha=0.92))
        ax3.text(17.5, 78, r'$m_3 = 103$' + '\n' + r'($m_{\max}$)', fontsize=8, ha='center', va='center',
                 bbox=dict(boxstyle='square,pad=0.2', facecolor='white', edgecolor='#333333', alpha=0.92))

        # Batch Delay Indicator
        ax3.annotate('', xy=(4.8, 55), xytext=(0.2, 55),
                     arrowprops=dict(arrowstyle='<->', color='#1d568c', lw=1.3))
        ax3.text(2.5, 58, 'Batch delay\n' + r'$\Delta k = p = 5$', fontsize=7.2,
                 ha='center', va='bottom', color='#1d568c', fontweight='bold',
                 bbox=dict(boxstyle='round,pad=0.15', facecolor='#edf4fc', edgecolor='#1d568c', alpha=0.92))

    elif matrix_name == 'young3c':
        ax3.text(2.5, 14, r'$m_0 = 30$' + '\n(initial)', fontsize=7.5, ha='center', va='center',
                 bbox=dict(boxstyle='square,pad=0.15', facecolor='white', edgecolor='#888888', alpha=0.92))
        ax3.text(7.5, 62, r'$m_1 = 76$' + '\n' + r'(from $\mathcal{M}_1$)', fontsize=7.5, ha='center', va='center',
                 bbox=dict(boxstyle='square,pad=0.15', facecolor='#fffbfa', edgecolor='#b22222', alpha=0.92))
        ax3.text(12.5, 38, r'$m_2 = 51$' + '\n' + r'(from $\mathcal{M}_2$)', fontsize=7.5, ha='center', va='center',
                 bbox=dict(boxstyle='square,pad=0.15', facecolor='#fffbfa', edgecolor='#b22222', alpha=0.92))
        ax3.text(17.5, 22, r'$m_3 = 34$', fontsize=7.5, ha='center', va='center',
                 bbox=dict(boxstyle='square,pad=0.15', facecolor='white', edgecolor='#333333', alpha=0.92))

        # Batch Delay Indicator en la zona superior abierta
        ax3.annotate('', xy=(4.8, 92), xytext=(0.2, 92),
                     arrowprops=dict(arrowstyle='<->', color='#1d568c', lw=1.3))
        ax3.text(6.0, 92, r'$\leftarrow$ Batch delay: $\Delta k = p = 5$', fontsize=7.8,
                 ha='left', va='center', color='#1d568c', fontweight='bold',
                 bbox=dict(boxstyle='round,pad=0.2', facecolor='#edf4fc', edgecolor='#1d568c', alpha=0.95))

    elif matrix_name == 'wang2':
        ax3.text(2.5, 16, r'$m_0 = 30$' + '\n(initial)', fontsize=8, ha='center', va='center',
                 bbox=dict(boxstyle='square,pad=0.2', facecolor='white', edgecolor='#888888', alpha=0.92))
        ax3.text(7.5, 38, r'$m_1 = 57$' + '\n' + r'(from $\mathcal{M}_1$)', fontsize=8, ha='center', va='center',
                 bbox=dict(boxstyle='square,pad=0.2', facecolor='#fffbfa', edgecolor='#b22222', alpha=0.92))
        ax3.text(12.5, 25, r'$m_2 = 12$' + '\n' + r'(from $\mathcal{M}_2$)', fontsize=8, ha='center', va='center',
                 bbox=dict(boxstyle='square,pad=0.2', facecolor='#fffbfa', edgecolor='#b22222', alpha=0.92))
        ax3.text(17.5, 45, r'$m_3 = 63$' + '\n' + r'(from $\mathcal{M}_3$)', fontsize=8, ha='center', va='center',
                 bbox=dict(boxstyle='square,pad=0.2', facecolor='#fffbfa', edgecolor='#b22222', alpha=0.92))
        ax3.text(21.5, 42, r'$m_4 = 59$', fontsize=7.5, ha='center', va='center',
                 bbox=dict(boxstyle='square,pad=0.15', facecolor='white', edgecolor='#333333', alpha=0.92))

        # Batch Delay Indicator en espacio superior amplio
        ax3.annotate('', xy=(4.8, 120), xytext=(0.2, 120),
                     arrowprops=dict(arrowstyle='<->', color='#1d568c', lw=1.3))
        ax3.text(2.5, 126, 'Batch delay\n' + r'$\Delta k = p = 5$', fontsize=7.5,
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

    # Ajuste de geometría fija garantizando perfecta alineación horizontal
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

    # Conectores causales rectos entre subgráficos (Flechas de realimentación directas en x=5 y x=10)
    con1 = ConnectionPatch(xyA=(5, stages[0]['error']), coordsA='data', axesA=ax2,
                           xyB=(5, m1_target), coordsB='data', axesB=ax3,
                           arrowstyle="-|>", color='#b22222', lw=1.6, mutation_scale=13,
                           connectionstyle="arc3,rad=0", zorder=10)
    fig.add_artist(con1)
    
    con2 = ConnectionPatch(xyA=(10, stages[1]['error']), coordsA='data', axesA=ax2,
                           xyB=(10, m2_target), coordsB='data', axesB=ax3,
                           arrowstyle="-|>", color='#b22222', lw=1.6, mutation_scale=13,
                           connectionstyle="arc3,rad=0", zorder=10)
    fig.add_artist(con2)

    # Guardar figuras
    png_path = f"{output_prefix}.png"
    eps_path = f"{output_prefix}.eps"
    os.makedirs(os.path.dirname(os.path.abspath(png_path)), exist_ok=True)
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(eps_path, format='eps', bbox_inches='tight')
    plt.close()
    
    print(f"[OK] Gráfico causal ({matrix_name}) generado exitosamente:")
    print(f"     PNG: {png_path}")
    print(f"     EPS: {eps_path}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Graficar dinámica causal y delay en DMDc-LQR")
    parser.add_argument('--matrix', type=str, default='all', choices=['all', 'orsirr1', 'young3c', 'wang2'])
    args = parser.parse_args()
    
    matrices = ['young3c', 'orsirr1', 'wang2'] if args.matrix == 'all' else [args.matrix]
    
    for mat in matrices:
        cycle_file = f"results/tables/cycles_{mat}_dlk.csv"
        out_prefix = f"results/figures/Performance_causal_dynamics_{mat}"
        plot_causal_dynamics_matrix(mat, cycle_file, out_prefix)
