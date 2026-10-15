import scipy.io as sio
import numpy as np
import scipy.linalg as la
from scipy.sparse.linalg import gmres

mat = sio.loadmat("data/mat_collection_test/sherman5.mat")
A = mat['Problem']['A'][0, 0].tocsr()
b = mat['Problem']['b'][0, 0].flatten()
N = len(b)
norm_b = np.linalg.norm(b)

m_max = int(N * 0.1)
m_min = 2
m = 30
tol = 1e-9
cycle_max = 50
p = 5
q = 5
t = 6
alpha = 1e-4

x = np.zeros(N)
tot_cycles = 0
R_history = []
m_history = []
r0 = b - A @ x
R_history.append(np.linalg.norm(r0) / norm_b)
m_history.append(m)

print(f"Iniciando simulación con m0 = {m}, m_max = {m_max}, tol = {tol}")

while tot_cycles < cycle_max and R_history[-1] >= tol:
    print(f"Snapshot con m = {m} ...")
    MatR = []
    MatR.append((b - A @ x) / np.linalg.norm(b - A @ x))
    cycles_window = 0
    
    for c in range(p):
        x_new, info = gmres(A, b, x0=x, restart=m, maxiter=1, atol=1e-12)
        r = b - A @ x_new
        rel_res = np.linalg.norm(r) / norm_b
        R_history.append(rel_res)
        m_history.append(m)
        MatR.append(r / np.linalg.norm(r))
        x = x_new
        tot_cycles += 1
        cycles_window += 1
        if rel_res < tol:
            break
            
    if R_history[-1] < tol:
        break
        
    # DMDc
    StateData = np.column_stack(MatR)
    X = StateData[:, :-1]
    Xp = StateData[:, 1:]
    InputData = m * np.ones((1, StateData.shape[1]))
    Up = InputData[:, :-1]
    Omega = np.vstack([X, Up])

    U, Sig, Vt = np.linalg.svd(Omega, full_matrices=False)
    q_curr = min(q, len(Sig))
    Util = U[:, :q_curr]
    Sigtil = np.diag(Sig[:q_curr])
    Vtil = Vt.T[:, :q_curr]

    U_s, _, _ = np.linalg.svd(StateData, full_matrices=False)
    t_curr = min(t, U_s.shape[1])
    Uhat = U_s[:, :t_curr]

    U_1 = Util[:N, :]
    U_2 = Util[N:, :]

    Htil = Xp @ Vtil @ np.linalg.inv(Sigtil) @ U_1.T
    Btil = Xp @ Vtil @ np.linalg.inv(Sigtil) @ U_2.T

    Hbar = Uhat.T @ Htil @ Uhat
    Bbar = Uhat.T @ Btil

    # Trajectory of xhat_p
    Xhat_p = np.zeros((t_curr, StateData.shape[1]))
    Xhat_p[:, 0] = Uhat.T @ StateData[:, 0]
    for i in range(StateData.shape[1] - 1):
        Xhat_p[:, i + 1] = Hbar @ Xhat_p[:, i] + Bbar.flatten() * m

    xhat = Xhat_p[:, -1]

    # Controllability
    Co = np.column_stack([np.linalg.matrix_power(Hbar, k) @ Bbar for k in range(t_curr)])
    r = np.linalg.matrix_rank(Co)

    try:
        if r < t_curr:
            U_co, _, _ = np.linalg.svd(Co, full_matrices=False)
            T = U_co[:, :r]
            Hc = T.T @ Hbar @ T
            Bc = T.T @ Bbar
            Qc = np.eye(r)
            Rc = np.eye(Bc.shape[1])
            P = la.solve_continuous_are(Hc, Bc, Qc, Rc)
            Kc = np.linalg.inv(Rc) @ Bc.T @ P
            xhatc = T.T @ xhat
            xhatc = xhatc / np.linalg.norm(xhatc)
            u_opt = -Kc @ xhatc
            u = int(np.floor(np.linalg.norm(u_opt)))
        else:
            Q = np.eye(t_curr)
            R = np.array([[alpha]])
            P = la.solve_continuous_are(Hbar, Bbar, Q, R)
            K = (1.0 / alpha) * Bbar.T @ P
            xhat_norm = xhat / np.linalg.norm(xhat)
            u_opt = -K @ xhat_norm
            u = int(np.floor(np.linalg.norm(u_opt)))
    except Exception as e:
        print(f"LQR error: {e}, keeping m")
        u = m

    m_new = min(u, m_max)
    m_new = max(m_new, m_min)
    print(f"  Ciclo {tot_cycles}: rel_res = {R_history[-1]:.2e}, LQR u = {u} -> nuevo m = {m_new}")
    m = m_new

print(f"\nFinalizado en {tot_cycles} ciclos con residuo final = {R_history[-1]:.2e}")
