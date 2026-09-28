#!/usr/bin/env python3
import os, csv, glob, argparse
import matplotlib.pyplot as plt

def read_history_csv(fpath):
    iterations, rel_residuals = [], []
    with open(fpath, mode='r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            iterations.append(int(row['iteration']))
            rel_residuals.append(float(row.get('rel_residual', row['residual'])))
    return iterations, rel_residuals

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="*")
    parser.add_argument("-o", "--output", default="convergence_dlk.png")
    parser.add_argument("-t", "--title", default="Comparación DLK-GMRES vs GMRES Nativos")
    args = parser.parse_args()

    files = args.files or glob.glob("results/history_*.csv")
    if not files: return

    plt.figure(figsize=(10, 6), dpi=200)
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    
    styles = ['-', '--', '-.', ':']
    markers = ['o', 's', '^', 'd', 'v', 'x']

    for idx, fpath in enumerate(sorted(files)):
        label = os.path.basename(fpath).replace("history_", "").replace(".csv", "").replace("_", " ")
        try:
            x, y = read_history_csv(fpath)
            if not x: continue
            me = max(1, len(x) // 15)
            plt.semilogy(x, y, label=label, linestyle=styles[idx % len(styles)],
                         marker=markers[idx % len(markers)], markevery=me, linewidth=2, markersize=5)
        except Exception as e:
            print(f"Error {fpath}: {e}")

    plt.title(args.title, fontsize=14, fontweight='bold', pad=15)
    plt.xlabel("Iteraciones de Krylov ($k$)", fontsize=12)
    plt.ylabel(r"Residuo Relativo $\|r_k\| / \|r_0\|$", fontsize=12)
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend(fontsize=10, loc="best", framealpha=0.9)
    plt.tight_layout()
    plt.savefig(args.output, bbox_inches='tight')
    plt.close()
    print(f"[✓] Curva de convergencia guardada en: {args.output}")

if __name__ == '__main__':
    main()
