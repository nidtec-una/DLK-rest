#!/usr/bin/env python3
"""
scripts/compute_condition_numbers.py

Calcula el número de condicionamiento cond_1(A) = ||A||_1 * ||A^-1||_1
mediante el algoritmo de estimación matricial de 1-norma de Higham & Tisseur (2000)
acoplado con factorización SuperLU (splu) para las 9 matrices de flujo fracturado.
"""

import os
import sys
import time
import csv
import numpy as np
import scipy.io as sio
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu, onenormest, LinearOperator

def compute_condition_numbers(mat_dir=None):
    if mat_dir is None:
        mat_dir = "/mnt/c/Users/Carlos Benitez/Downloads/Matrices ixed-Dimensional Problems with High Permeability Ratios-20261015T190037Z-1-001/Matrices ixed-Dimensional Problems with High Permeability Ratios"

    files = [
        ("Kf1e-08", "flow_benchmark_2d_case_1_Kf1e-08.mat", "1.0e-08", "Barrera casi impermeable"),
        ("Kf1e-06", "flow_benchmark_2d_case_1_Kf1e-06.mat", "1.0e-06", "Fractura de baja conductividad"),
        ("Kf1e-04", "flow_benchmark_2d_case_1_Kf1e-04.mat", "1.0e-04", "Fractura de resistencia moderada"),
        ("Kf1e-02", "flow_benchmark_2d_case_1_Kf1e-02.mat", "1.0e-02", "Transición matriz-fractura"),
        ("Kf1e+00", "flow_benchmark_2d_case_1_Kf1e+00.mat", "1.0e+00", "Permeabilidad homogénea"),
        ("Kf1e+02", "flow_benchmark_2d_case_1_Kf1e+02.mat", "1.0e+02", "Fractura moderadamente conductora"),
        ("Kf1e+04", "flow_benchmark_2d_case_1_Kf1e+04.mat", "1.0e+04", "Fractura conductora"),
        ("Kf1e+06", "flow_benchmark_2d_case_1_Kf1e+06.mat", "1.0e+06", "Canal altamente conductor"),
        ("Kf1e+08", "flow_benchmark_2d_case_1_Kf1e+08.mat", "1.0e+08", "Canal hiper-conductor"),
    ]

    print("=" * 95)
    print(" CÁLCULO DE CONDICIONAMIENTO cond(A) PARA MATRICES DE FLUJO (HIGHAM & TISSEUR 1-NORM)")
    print("=" * 95)
    print(f"{'#':<2} | {'Tag':<10} | {'Permeabilidad Kf':<18} | {'||A||_1':<12} | {'||A^-1||_1':<12} | {'cond_1(A)':<14} | {'Tiempo'}")
    print("-" * 95)
    sys.stdout.flush()

    cond_results = {}

    for idx, (tag, fname, kf_val, regime) in enumerate(files, 1):
        fpath = os.path.join(mat_dir, fname)
        if not os.path.exists(fpath):
            print(f"{idx:<2} | {tag:<10} | No encontrado: {fname}")
            continue

        t0 = time.time()
        try:
            m = sio.loadmat(fpath)
            A = csc_matrix(m['A'].astype(np.float64))

            # 1. Norma 1 de A
            norm_A = float(onenormest(A))

            # 2. Factorización SuperLU
            lu = splu(A, diag_pivot_thresh=0.1)

            # 3. Norma 1 de A^-1 mediante operador inverso
            A_inv_op = LinearOperator(A.shape, matvec=lu.solve, rmatvec=lambda x: lu.solve(x, trans='T'))
            norm_A_inv = float(onenormest(A_inv_op))

            # 4. Condicionamiento
            cond_val = norm_A * norm_A_inv
            elapsed = time.time() - t0

            cond_str = f"{cond_val:.4e}"
            cond_results[tag] = cond_str

            print(f"{idx:<2} | {tag:<10} | {kf_val:<18} | {norm_A:<12.4e} | {norm_A_inv:<12.4e} | {cond_str:<14} | {elapsed:.2f}s")
            sys.stdout.flush()

        except Exception as e:
            print(f"{idx:<2} | {tag:<10} | Error calculando condicionamiento: {e}")
            cond_results[tag] = "N/A"
            sys.stdout.flush()

    print("=" * 95)
    return cond_results

if __name__ == "__main__":
    compute_condition_numbers()
