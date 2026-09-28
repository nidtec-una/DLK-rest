import sys
sys.path.append("/usr/share/petsc/3.24/lib/petsc/bin")
import PetscBinaryIO
import scipy.io as sio
from scipy.sparse import csr_matrix
import numpy as np

mat = sio.loadmat("data/mat_collection_test/sherman5.mat")
problem = mat['Problem'][0, 0]
A = csr_matrix(problem['A'])
b = problem['b'].flatten().astype(np.float64)

io = PetscBinaryIO.PetscBinaryIO()
# En PetscBinaryIO, los vectores deben ser de tipo PetscBinaryIO.Vec
petsc_vec = b.view(PetscBinaryIO.Vec)

io.writeBinaryFile("data/sherman5.bin", [A, petsc_vec])
print("Matriz A y vector b guardados exitosamente en data/sherman5.bin")
