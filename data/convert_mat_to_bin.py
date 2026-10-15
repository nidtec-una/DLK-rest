#!/usr/bin/env python3
"""
convert_all_mat_to_petsc.py

Convierte ordenadamente todos los archivos .mat (SuiteSparse / MATLAB) a formato
binario nativo de PETSc (.bin), empaquetando la matriz A y el vector del lado derecho b.
Si la matriz en el archivo .mat no contiene el vector b, se genera un vector b = ones(N, 1),
tal como se hace en la versión de MATLAB.
"""

import os
import sys
import glob
import numpy as np
import scipy.io as sio
from scipy.sparse import issparse, csr_matrix

# Configurar path para PetscBinaryIO
sys.path.append("/usr/share/petsc/3.24/lib/petsc/bin")
try:
    import PetscBinaryIO
except ImportError:
    petsc_dir = os.environ.get("PETSC_DIR", "/usr/lib/petscdir/petsc3.24/x86_64-linux-gnu-real")
    sys.path.append(os.path.join(petsc_dir, "lib/petsc/bin"))
    sys.path.append(os.path.join(petsc_dir, "share/petsc/bin"))
    import PetscBinaryIO

def inspect_and_extract(mat_file):
    mat = sio.loadmat(mat_file)
    A = None
    b = None
    source_b = "generado ones(N, 1)"
    
    # Caso 1: Estructura Problem estándar de SuiteSparse
    if 'Problem' in mat:
        p = mat['Problem'][0, 0]
        fields = p.dtype.names
        if 'A' in fields:
            A = p['A']
        if 'b' in fields:
            b_cand = p['b']
            if b_cand is not None and b_cand.size > 0:
                b = b_cand.flatten().astype(np.float64)
                source_b = "original del .mat"
    # Caso 2: Claves directas A y b
    elif 'A' in mat:
        A = mat['A']
        if 'b' in mat:
            b_cand = mat['b']
            if b_cand is not None and b_cand.size > 0:
                b = b_cand.flatten().astype(np.float64)
                source_b = "original del .mat"
                
    if A is None:
        raise ValueError(f"No se encontró la matriz A en {mat_file}")
        
    # Extraer arrays CSR asegurando tipos
    indptr = getattr(A, 'indptr', None)
    indices = getattr(A, 'indices', None)
    data = getattr(A, 'data', None)
    
    if indptr is None or indices is None or data is None:
        # Si es matriz densa u otro formato, convertir a CSR
        A_csr = csr_matrix(A)
        indptr = A_csr.indptr
        indices = A_csr.indices
        data = A_csr.data
        shape = A_csr.shape
    else:
        shape = A.shape
        
    M, N = shape
    
    # Si no tiene b, completar con vector de 1s como en MATLAB
    if b is None or len(b) != N:
        b = np.ones(N, dtype=np.float64)
        source_b = "generado ones(N, 1)"
        
    # Construir objeto MatSparse nativo de PetscBinaryIO
    # ((M, N), (indptr, indices, data))
    petsc_mat = PetscBinaryIO.MatSparse(((M, N), (np.asarray(indptr), np.asarray(indices), np.asarray(data, dtype=np.float64))))
    petsc_vec = b.astype(np.float64).view(PetscBinaryIO.Vec)
    
    return petsc_mat, petsc_vec, M, N, len(data), source_b

def main():
    io = PetscBinaryIO.PetscBinaryIO()
    
    # Determinar rutas relativas al script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.abspath(os.path.join(script_dir, ".."))
    
    search_dirs = [
        os.path.join(script_dir, "mat_collection_test_full"),
        os.path.join(script_dir, "mat_collection_test")
    ]
    
    mat_files = {}
    for d in search_dirs:
        if os.path.exists(d):
            for f in sorted(os.listdir(d)):
                if f.endswith(".mat"):
                    name = f[:-4]
                    if name not in mat_files:
                        mat_files[name] = os.path.join(d, f)
                        
    output_dir = os.path.join(script_dir, "bin")
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 105)
    print(" CONVERSIÓN MASIVA Y ORDENADA DE MATRICES (.mat -> .bin)")
    print("=" * 105)
    print(f"{'#':<3} | {'Matriz':<18} | {'Filas x Cols':<16} | {'NNZ':<10} | {'Origen Vector b':<22} | {'Archivo Destino'}")
    print("-" * 105)
    
    idx = 1
    converted_records = []
    
    for name in sorted(mat_files.keys()):
        mat_path = mat_files[name]
        try:
            petsc_mat, petsc_vec, M, N, nnz, source_b = inspect_and_extract(mat_path)
            
            # Guardar en data/bin/<name>.bin
            bin_path = os.path.join(output_dir, f"{name}.bin")
            
            io.writeBinaryFile(bin_path, [petsc_mat, petsc_vec])
            
            dim_str = f"{M} x {N}"
            rel_bin_path = os.path.relpath(bin_path, repo_root)
            print(f"{idx:<3} | {name:<18} | {dim_str:<16} | {nnz:<10} | {source_b:<22} | {rel_bin_path}")
            converted_records.append({
                'name': name,
                'M': M,
                'N': N,
                'nnz': nnz,
                'source_b': source_b,
                'bin_path': bin_path
            })
            idx += 1
        except Exception as e:
            print(f"[ERROR] Error procesando {name} ({mat_path}): {e}")
            
    print("-" * 105)
    print(f"[✓] Se convirtieron exitosamente {len(converted_records)} de {len(mat_files)} matrices a formato binario PETSc (.bin).")
    print("=" * 105)

if __name__ == "__main__":
    main()
