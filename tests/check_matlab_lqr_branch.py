import scipy.io as sio
import numpy as np
import scipy.linalg as la

mat = sio.loadmat("data/mat_collection_test/sherman5.mat")
A = mat['Problem']['A'][0, 0].tocsr()
b = mat['Problem']['b'][0, 0].flatten()
N = len(b)

from scipy.sparse.linalg import gmres

x = np.zeros(N)
MatR = []
r0 = b - A @ x
MatR.append(r0 / np.linalg.norm(r0))

p = 5
m = 30
for c in range(p):
    x_new, info = gmres(A, b, x0=x, restart=m, maxiter=1, atol=1e-12)
    r = b - A @ x_new
    MatR.append(r / np.linalg.norm(r))
    x = x_new

StateData = np.column_stack(MatR)

q = 5
t = 6
X = StateData[:, :-1]
Xp = StateData[:, 1:]
InputData = m * np.ones((1, StateData.shape[1]))
Up = InputData[:, :-1]
Omega = np.vstack([X, Up])

U, Sig, Vt = np.linalg.svd(Omega, full_matrices=False)
Util = U[:, :q]
Sigtil = np.diag(Sig[:q])
Vtil = Vt.T[:, :q]

U_s, _, _ = np.linalg.svd(StateData, full_matrices=False)
Uhat = U_s[:, :t]

U_1 = Util[:N, :]
U_2 = Util[N:, :]

Htil = Xp @ Vtil @ np.linalg.inv(Sigtil) @ U_1.T
Btil = Xp @ Vtil @ np.linalg.inv(Sigtil) @ U_2.T

Hbar = Uhat.T @ Htil @ Uhat
Bbar = Uhat.T @ Btil

# Trajectory of xhat_p
Xhat_p = np.zeros((t, p + 1))
Xhat_p[:, 0] = Uhat.T @ StateData[:, 0]
for i in range(p):
    Xhat_p[:, i + 1] = Hbar @ Xhat_p[:, i] + Bbar.flatten() * m

xhat = Xhat_p[:, -1]
print("xhat:", xhat)

Co = np.column_stack([np.linalg.matrix_power(Hbar, k) @ Bbar for k in range(t)])
r = np.linalg.matrix_rank(Co)

if r < t:
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
    print("u calculado (rama r < t):", u)
else:
    Q = np.eye(t)
    R = np.array([[1e-4]])
    P = la.solve_continuous_are(Hbar, Bbar, Q, R)
    K = (1.0 / 1e-4) * Bbar.T @ P
    xhat_norm = xhat / np.linalg.norm(xhat)
    u_opt = -K @ xhat_norm
    u = int(np.floor(np.linalg.norm(u_opt)))
    print("u calculado (rama r == t):", u)

m_max = int(N * 0.1)
m_new = min(u, m_max)
print(f"m_new limitado por m_max ({m_max}):", m_new)
