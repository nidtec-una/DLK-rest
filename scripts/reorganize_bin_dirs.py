#!/usr/bin/env python3
"""
scripts/reorganize_bin_dirs.py

Reorganiza los archivos .bin dentro de data/bin en subdirectorios temáticos:
- data/bin/suitesparse/ : Matrices de la colección SuiteSparse (22 matrices del benchmark)
- data/bin/flow_benchmark/ : Matrices de problemas de flujo en medios fracturados (9 matrices de contraste)
Crea enlaces simbólicos relativos en data/bin/ para preservar total compatibilidad con scripts existentes.
"""

import os
import shutil
import glob

def main():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    bin_dir = os.path.join(repo_root, "data", "bin")
    
    suitesparse_dir = os.path.join(bin_dir, "suitesparse")
    flow_dir = os.path.join(bin_dir, "flow_benchmark")
    
    os.makedirs(suitesparse_dir, exist_ok=True)
    os.makedirs(flow_dir, exist_ok=True)
    
    # 1. Identificar y mover flow_benchmark matrices
    flow_files = glob.glob(os.path.join(bin_dir, "flow_benchmark_2d_case_1_Kf*.bin"))
    for f in flow_files:
        base = os.path.basename(f)
        dest = os.path.join(flow_dir, base)
        if not os.path.islink(f):
            shutil.move(f, dest)
            print(f"Moved to flow_benchmark: {base}")
            
    # 2. Identificar y mover SuiteSparse matrices que estén directamente en data/bin
    suitesparse_names = [
        "add20", "cavity05", "cavity10", "cdde1", "circuit_2", "circuit2",
        "ex40", "fpga_trans_01", "fpgatrans01", "memplus", "orsirr_1", "orsirr1",
        "orsreg_1", "orsreg1", "pde2961", "raefsky1", "raefsky2", "rdb2048",
        "sherman1", "sherman4", "sherman5", "steam2", "wang2", "wang4",
        "watt_1", "watt1", "young3c"
    ]
    
    for name in suitesparse_names:
        cand = os.path.join(bin_dir, f"{name}.bin")
        dest = os.path.join(suitesparse_dir, f"{name}.bin")
        if os.path.exists(cand) and not os.path.islink(cand):
            shutil.move(cand, dest)
            print(f"Moved to suitesparse: {name}.bin")
        elif os.path.islink(cand):
            # Si era symlink a otro archivo en bin_dir, recrearlo en suitesparse
            target = os.readlink(cand)
            os.remove(cand)
            # Recrear en suitesparse si no existe
            if not os.path.exists(dest):
                os.symlink(target, dest)
                print(f"Created symlink in suitesparse: {name}.bin -> {target}")

    # 3. Crear symlinks relativos en data/bin para compatibilidad hacia atrás
    # Para suitesparse:
    for f in sorted(os.listdir(suitesparse_dir)):
        if f.endswith(".bin"):
            symlink_path = os.path.join(bin_dir, f)
            rel_target = os.path.join("suitesparse", f)
            if os.path.islink(symlink_path) or os.path.exists(symlink_path):
                os.remove(symlink_path)
            os.symlink(rel_target, symlink_path)
            print(f"Symlink: data/bin/{f} -> {rel_target}")
            
    # Para flow_benchmark:
    for f in sorted(os.listdir(flow_dir)):
        if f.endswith(".bin"):
            symlink_path = os.path.join(bin_dir, f)
            rel_target = os.path.join("flow_benchmark", f)
            if os.path.islink(symlink_path) or os.path.exists(symlink_path):
                os.remove(symlink_path)
            os.symlink(rel_target, symlink_path)
            print(f"Symlink: data/bin/{f} -> {rel_target}")
            
    print("[OK] Reorganización completada con éxito.")

if __name__ == "__main__":
    main()
