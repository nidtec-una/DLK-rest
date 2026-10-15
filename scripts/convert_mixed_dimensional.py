#!/usr/bin/env python3
"""
scripts/convert_mixed_dimensional.py

Convierte el conjunto de matrices de problemas dimensionales mixtos con altos contrastes de permeabilidad
(Mixed-Dimensional Problems with High Permeability Ratios) desde formato .mat (MATLAB / SciPy) a formato
binario nativo de PETSc (.bin).
"""

import os
import sys
import glob
import numpy as np
import scipy.io as sio
from scipy.sparse import csr_matrix, issparse

# Configuración de ruta para PetscBinaryIO
sys.path.append("/usr/share/petsc/3.24/lib/petsc/bin")
try:
    import PetscBinaryIO
except ImportError:
    petsc_dir = os.environ.get("PETSC_DIR", "/usr/lib/petscdir/petsc3.24/x86_64-linux-gnu-real")
    for cand in ["lib/petsc/bin", "share/petsc/bin"]:
        p = os.path.join(petsc_dir, cand)
        if p not in sys.path:
            sys.path.append(p)
    import PetscBinaryIO

def convert_matrices(source_dir, dest_dirs):
    mat_files = sorted(glob.glob(os.path.join(source_dir, "*.mat")))
    if not mat_files:
        print(f"[ERROR] No se encontraron archivos .mat en {source_dir}")
        return []

    for d in dest_dirs:
        os.makedirs(d, exist_ok=True)

    io = PetscBinaryIO.PetscBinaryIO()
    results = []

    print("=" * 115)
    print(" CONVERSIÓN A FORMATO BINARIO PETSc (.bin): MIXED-DIMENSIONAL HIGH PERMEABILITY RATIO MATRICES")
    print("=" * 115)
    print(f"{'#':<2} | {'Matriz':<35} | {'Dimensiones':<13} | {'NNZ':<8} | {'Contrast':<10} | {'||b||_2':<10} | Estado")
    print("-" * 115)

    for idx, mat_path in enumerate(mat_files, 1):
        base_name = os.path.splitext(os.path.basename(mat_path))[0]
        try:
            mat = sio.loadmat(mat_path)
            
            # Extraer A
            if 'A' not in mat:
                raise KeyError(f"Clave 'A' no encontrada en {mat_path}")
            raw_A = mat['A']
            
            if not issparse(raw_A):
                A_csr = csr_matrix(raw_A)
            elif not isinstance(raw_A, csr_matrix):
                A_csr = raw_A.tocsr()
            else:
                A_csr = raw_A
                
            M, N = A_csr.shape
            nnz = A_csr.nnz
            
            # Asegurar tipos estándar para PETSc (CSR)
            indptr = np.asarray(A_csr.indptr, dtype=np.int32)
            indices = np.asarray(A_csr.indices, dtype=np.int32)
            data = np.asarray(A_csr.data, dtype=np.float64)
            
            # Extraer b
            if 'b' in mat and mat['b'] is not None and mat['b'].size > 0:
                b = mat['b'].flatten().astype(np.float64)
            else:
                b = np.ones(N, dtype=np.float64)
                
            norm_b = np.linalg.norm(b)
            
            # Parámetros físicos adicionales si existen
            contrast_val = float(mat['contrast'].flatten()[0]) if 'contrast' in mat else float('nan')
            
            # Construir objetos binarios de PETSc
            petsc_mat = PetscBinaryIO.MatSparse(((M, N), (indptr, indices, data)))
            petsc_vec = b.view(PetscBinaryIO.Vec)
            
            # Guardar en cada directorio destino
            saved_paths = []
            for out_dir in dest_dirs:
                bin_path = os.path.join(out_dir, f"{base_name}.bin")
                io.writeBinaryFile(bin_path, [petsc_mat, petsc_vec])
                saved_paths.append(bin_path)
                
            # Verificación de integridad: leer de vuelta
            test_read = io.readBinaryFile(saved_paths[0])
            if len(test_read) != 2:
                raise ValueError("La verificación de lectura binaria no retornó [Mat, Vec]")
                
            status_str = "[OK] Verificado"
            print(f"{idx:<2} | {base_name:<35} | {f'{M} x {N}':<13} | {nnz:<8} | {contrast_val:<10.1e} | {norm_b:<10.2e} | {status_str}")
            
            results.append({
                'name': base_name,
                'M': M,
                'N': N,
                'nnz': nnz,
                'contrast': contrast_val,
                'norm_b': norm_b,
                'destinations': saved_paths
            })
            
        except Exception as e:
            print(f"{idx:<2} | {base_name:<35} | ERROR: {e}")

    print("-" * 115)
    print(f"[✓] Éxito: {len(results)} de {len(mat_files)} matrices convertidas y validadas correctamente.")
    print("=" * 115)
    return results

if __name__ == '__main__':
    source = "/mnt/c/Users/Carlos Benitez/Downloads/Matrices ixed-Dimensional Problems with High Permeability Ratios-20261015T190037Z-1-001/Matrices ixed-Dimensional Problems with High Permeability Ratios"
    
    # Destinos:
    # 1. En la misma carpeta de descargas del usuario (al lado de los .mat)
    # 2. En el repositorio DLK-rest (data/bin)
    repo_bin = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data", "bin"))
    destinations = [source, repo_bin]
    
    convert_matrices(source, destinations)
