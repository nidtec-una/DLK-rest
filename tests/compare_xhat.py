import scipy.io as sio
import numpy as np
import scipy.linalg as la
from scipy.sparse.linalg import gmres

mat = sio.loadmat('data/mat_collection_test_full/sherman5.mat')
A = mat['Problem']['A'][0, 0].tocsr()
b = mat['Problem']['b'][0, 0].flatten()
N = len(b)

x = np.zeros(N)
for snap in range(2):
    m = 30 if snap == 0 else 331
    MatR = [(b - A @ x) / np.linalg.norm(b - A @ x)]
    for c in range(5):
        x, _ = gmres(A, b, x0=x, restart=m, maxiter=1, atol=1e-12)
        r = b - A @ x
        MatR.append(r / np.linalg.norm(r))
    
    StateData = np.column_stack(MatR)
    X = StateData[:, :-1]
    Xp = StateData[:, 1:]
    InputData = m * np.ones((1, StateData.shape[1]))
    Up = InputData[:, :-1]
    Omega = np.vstack([X, Up])
    
    U, Sig, Vt = np.linalg.svd(Omega, full_matrices=False)
    Util = U[:, :5]
    Sigtil = np.diag(Sig[:5])
    Vtil = Vt.T[:, :5]
    
    U_s, _, _ = np.linalg.svd(StateData, full_matrices=False)
    Uhat = U_s[:, :6]
    
    U_1 = Util[:N, :]
    U_2 = Util[N:, :]
    Htil = Xp @ Vtil @ np.linalg.inv(Sigtil) @ U_1.T
    Btil = Xp @ Vtil @ np.linalg.inv(Sigtil) @ U_2.T
    Hbar = Uhat.T @ Htil @ Uhat
    Bbar = Uhat.T @ Btil
    
    xhat_act = Uhat.T @ StateData[:, -1]
    
    Co = np.column_stack([np.linalg.matrix_power(Hbar, k) @ Bbar for k in range(6)])
    r = np.linalg.matrix_rank(Co)
    U_co, S_co, _ = np.linalg.svd(Co, full_matrices=False)
    T = U_co[:, :r]
    Hc = T.T @ Hbar @ T
    Bc = T.T @ Bbar
    P = la.solve_continuous_are(Hc, Bc, np.eye(r), np.eye(Bc.shape[1]))
    Kc = np.linalg.inv(np.eye(Bc.shape[1])) @ Bc.T @ P
    
    xc_act = T.T @ xhat_act
    xc_act /= np.linalg.norm(xc_act)
    u_opt = -Kc @ xc_act
    u_act = int(np.floor(np.linalg.norm(u_opt)))
    
    print(f"Snap {snap+1} (m={m}):")
    print(f"  u_opt: {u_opt}")
    print(f"  u_act: {u_act}")
    print(f"  S_co: {S_co}")
    print(f"  rank(Co): {r}")
    if snap == 1:
        print(f"  Hbar:\n{Hbar}")
        print(f"  Bbar:\n{Bbar.flatten()}")
        print(f"  xhat_act:\n{xhat_act}")
