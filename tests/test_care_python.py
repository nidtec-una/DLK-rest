import numpy as np
import scipy.linalg as la

def solve_care_hamiltonian(H, B, Q, alpha):
    t = H.shape[0]
    R_inv = 1.0 / alpha
    # Hamiltonian matrix (2t x 2t)
    Z = np.zeros((2 * t, 2 * t))
    Z[:t, :t] = H
    Z[:t, t:] = -R_inv * (B @ B.T)
    Z[t:, :t] = -Q
    Z[t:, t:] = -H.T

    # Eigenvalues and eigenvectors
    w, vr = np.linalg.eig(Z)
    
    # Sort by real part, select the t eigenvalues with Re(lambda) < 0
    stable_idx = [i for i, val in enumerate(w) if val.real < -1e-12]
    if len(stable_idx) != t:
        # Sort by real part ascending and take first t
        stable_idx = np.argsort(w.real)[:t]

    # Form real invariant subspace basis
    V_real = np.zeros((2 * t, t))
    i = 0
    used = set()
    for idx in stable_idx:
        if idx in used:
            continue
        val = w[idx]
        vec = vr[:, idx]
        if abs(val.imag) < 1e-10:
            V_real[:, i] = vec.real
            used.add(idx)
            i += 1
        else:
            # find conjugate
            V_real[:, i] = vec.real
            V_real[:, i + 1] = vec.imag
            used.add(idx)
            # find matching conjugate
            for c_idx in stable_idx:
                if c_idx not in used and abs(w[c_idx] - np.conj(val)) < 1e-8:
                    used.add(c_idx)
                    break
            i += 2
        if i >= t:
            break

    U11 = V_real[:t, :t]
    U21 = V_real[t:, :t]
    
    # P = U21 @ inv(U11)
    P = np.linalg.solve(U11.T, U21.T).T
    P = 0.5 * (P + P.T)
    K = R_inv * (B.T @ P)
    return P, K

# Test on 10 random systems
for seed in range(5):
    np.random.seed(seed)
    t = 6
    H = np.random.randn(t, t) * 0.5
    B = np.random.randn(t, 1)
    Q = np.eye(t)
    alpha = 1e-4

    P_scipy = la.solve_continuous_are(H, B, Q, np.array([[alpha]]))
    K_scipy = (1.0 / alpha) * (B.T @ P_scipy)

    P_ham, K_ham = solve_care_hamiltonian(H, B, Q, alpha)

    diff_P = np.linalg.norm(P_scipy - P_ham)
    diff_K = np.linalg.norm(K_scipy - K_ham)
    print(f"Seed {seed}: diff_P = {diff_P:.2e}, diff_K = {diff_K:.2e}, Norm K = {np.linalg.norm(K_ham):.1f}")
