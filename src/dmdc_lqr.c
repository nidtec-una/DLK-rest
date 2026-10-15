#include "dmdc_lqr.h"
#include <petscblaslapack.h>
#include <math.h>

const char *const DLKControllerTypes[] = {"NONE", "PD", "LQR", "DLKControllerType", "DLK_CTRL_", NULL};

/* Ley de control Proporcional-Derivativa (PD) */
PetscErrorCode DLK_UpdateByPD(PetscReal r_k, PetscReal r_km1, PetscReal r_km2,
                              PetscInt m_last, PetscInt *m_next)
{
    PetscFunctionBeginUser;
    if (r_km1 <= 1e-14) {
        *m_next = m_last;
        PetscFunctionReturn(PETSC_SUCCESS);
    }
    PetscReal term1 = -3.0 * r_k / r_km1;
    PetscReal term2 = 9.0 * (r_k - r_km2) / (2.0 * r_km1);
    PetscInt delta_m = - (PetscInt)PetscFloorReal(term1) + (PetscInt)PetscFloorReal(term2);
    *m_next = m_last + delta_m;
    PetscFunctionReturn(PETSC_SUCCESS);
}

/* SVD / Eigendecomposition de matriz simétrica semidefinida positiva pequeña (M x M) */
static PetscErrorCode DLK_SmallSymmEig(PetscInt n, const PetscReal *A, PetscReal *V, PetscReal *sigma)
{
    PetscBLASInt bn, lda, lwork, info;
    PetscReal   *work, *w;
    PetscReal   *Acopy;

    PetscFunctionBeginUser;
    PetscCall(PetscBLASIntCast(n, &bn));
    lda = bn;
    PetscCall(PetscMalloc1(n * n, &Acopy));
    PetscCall(PetscArraycpy(Acopy, A, n * n));
    PetscCall(PetscMalloc1(n, &w));

    /* Consulta de tamaño óptimo de trabajo */
    lwork = -1;
    PetscReal work_query;
    LAPACKsyev_("V", "U", &bn, Acopy, &lda, w, &work_query, &lwork, &info);
    lwork = (PetscBLASInt)work_query;
    if (lwork < 3 * bn) lwork = 3 * bn;
    PetscCall(PetscMalloc1(lwork, &work));

    LAPACKsyev_("V", "U", &bn, Acopy, &lda, w, work, &lwork, &info);
    PetscCheck(info == 0, PETSC_COMM_SELF, PETSC_ERR_LIB, "Fallo en LAPACKsyev_ con info = %d", (int)info);

    /* LAPACK ordena autovalores de menor a mayor; los invertimos a orden descendente */
    for (PetscInt i = 0; i < n; i++) {
        PetscInt rev = n - 1 - i;
        sigma[i] = w[rev] > 0.0 ? PetscSqrtReal(w[rev]) : 0.0;
        for (PetscInt r = 0; r < n; r++) {
            V[r + i * n] = Acopy[r + rev * n];
        }
    }

    PetscCall(PetscFree(Acopy));
    PetscCall(PetscFree(w));
    PetscCall(PetscFree(work));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/* Identificación DMDc a partir de los vectores de snapshot */
PetscErrorCode DLK_DMDc_Identify(PetscInt num_snaps, Vec *snapshots, const PetscReal *inputs,
                                 PetscInt q_rank, PetscInt t_rank,
                                 PetscReal *Hbar, PetscReal *Bbar, PetscReal *xhat_end,
                                 PetscReal *dmdc_error)
{
    PetscInt  P, N_total;
    PetscReal *G, *G_omega;
    PetscReal *V_S, *sig_S, *W_t;
    PetscReal *V_om, *sig_om, *M_q;
    PetscReal *A1, *A2, *U2;

    PetscFunctionBeginUser;
    P = num_snaps - 1;
    N_total = num_snaps;
    PetscCheck(P >= 2, PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "Se requieren al menos 3 snapshots para DMDc");

    if (t_rank > N_total) t_rank = N_total;
    if (q_rank > P)       q_rank = P;

    /* 1. Construir matriz de Gram G = R^T * R (tamaño N_total x N_total) */
    PetscCall(PetscCalloc1(N_total * N_total, &G));
    for (PetscInt j = 0; j < N_total; j++) {
        for (PetscInt i = j; i < N_total; i++) {
            PetscScalar dot;
            PetscCall(VecDot(snapshots[i], snapshots[j], &dot));
            G[i + j * N_total] = PetscRealPart(dot);
            G[j + i * N_total] = PetscRealPart(dot);
        }
    }

    /* 2. SVD de R (StateData) mediante G */
    PetscCall(PetscMalloc1(N_total * N_total, &V_S));
    PetscCall(PetscMalloc1(N_total, &sig_S));
    PetscCall(DLK_SmallSymmEig(N_total, G, V_S, sig_S));

    /* W_t = V_S(:, 0:t-1) * inv(diag(sig_S(0:t-1))) (tamaño N_total x t_rank) */
    PetscCall(PetscMalloc1(N_total * t_rank, &W_t));
    for (PetscInt k = 0; k < t_rank; k++) {
        PetscReal inv_sig = sig_S[k] > 1e-12 ? 1.0 / sig_S[k] : 0.0;
        for (PetscInt i = 0; i < N_total; i++) {
            W_t[i + k * N_total] = V_S[i + k * N_total] * inv_sig;
        }
    }

    /* 3. Construir G_omega = X^T * X + U_p^T * U_p (tamaño P x P) */
    PetscCall(PetscCalloc1(P * P, &G_omega));
    for (PetscInt j = 0; j < P; j++) {
        for (PetscInt i = 0; i < P; i++) {
            G_omega[i + j * P] = G[i + j * N_total] + inputs[i] * inputs[j];
        }
    }

    /* SVD de Omega */
    PetscCall(PetscMalloc1(P * P, &V_om));
    PetscCall(PetscMalloc1(P, &sig_om));
    PetscCall(DLK_SmallSymmEig(P, G_omega, V_om, sig_om));

    /* M_q = V_om(:, 0:q-1) * inv(diag(sig_om(0:q-1))) (tamaño P x q_rank) */
    PetscCall(PetscMalloc1(P * q_rank, &M_q));
    for (PetscInt k = 0; k < q_rank; k++) {
        PetscReal inv_sig = sig_om[k] > 1e-12 ? 1.0 / sig_om[k] : 0.0;
        for (PetscInt i = 0; i < P; i++) {
            M_q[i + k * P] = V_om[i + k * P] * inv_sig;
        }
    }

    /* 4. Calcular A1 = W_t^T * G(:, 1:P) * M_q (tamaño t_rank x q_rank) */
    /* Temporal T1 = G(:, 1:P) * M_q (tamaño N_total x q_rank) */
    PetscReal *T1;
    PetscCall(PetscCalloc1(N_total * q_rank, &T1));
    for (PetscInt i = 0; i < N_total; i++) {
        for (PetscInt k = 0; k < q_rank; k++) {
            PetscReal sum = 0.0;
            for (PetscInt j = 0; j < P; j++) {
                sum += G[i + (j + 1) * N_total] * M_q[j + k * P];
            }
            T1[i + k * N_total] = sum;
        }
    }

    PetscCall(PetscCalloc1(t_rank * q_rank, &A1));
    for (PetscInt r = 0; r < t_rank; r++) {
        for (PetscInt c = 0; c < q_rank; c++) {
            PetscReal sum = 0.0;
            for (PetscInt i = 0; i < N_total; i++) {
                sum += W_t[i + r * N_total] * T1[i + c * N_total];
            }
            A1[r + c * t_rank] = sum;
        }
    }
    PetscCall(PetscFree(T1));

    /* 5. Calcular A2 = M_q^T * G(0:P-1, :) * W_t (tamaño q_rank x t_rank) */
    PetscReal *T2;
    PetscCall(PetscCalloc1(P * t_rank, &T2));
    for (PetscInt i = 0; i < P; i++) {
        for (PetscInt c = 0; c < t_rank; c++) {
            PetscReal sum = 0.0;
            for (PetscInt j = 0; j < N_total; j++) {
                sum += G[i + j * N_total] * W_t[j + c * N_total];
            }
            T2[i + c * P] = sum;
        }
    }

    PetscCall(PetscCalloc1(q_rank * t_rank, &A2));
    for (PetscInt r = 0; r < q_rank; r++) {
        for (PetscInt c = 0; c < t_rank; c++) {
            PetscReal sum = 0.0;
            for (PetscInt i = 0; i < P; i++) {
                sum += M_q[i + r * P] * T2[i + c * P];
            }
            A2[r + c * q_rank] = sum;
        }
    }
    PetscCall(PetscFree(T2));

    /* 6. Hbar = A1 * A2 (tamaño t_rank x t_rank) */
    for (PetscInt r = 0; r < t_rank; r++) {
        for (PetscInt c = 0; c < t_rank; c++) {
            PetscReal sum = 0.0;
            for (PetscInt k = 0; k < q_rank; k++) {
                sum += A1[r + k * t_rank] * A2[k + c * q_rank];
            }
            Hbar[r + c * t_rank] = sum;
        }
    }

    /* 7. U2 = U_p * M_q (tamaño 1 x q_rank) */
    PetscCall(PetscCalloc1(q_rank, &U2));
    for (PetscInt k = 0; k < q_rank; k++) {
        PetscReal sum = 0.0;
        for (PetscInt i = 0; i < P; i++) {
            sum += inputs[i] * M_q[i + k * P];
        }
        U2[k] = sum;
    }

    /* 8. Bbar = A1 * U2^T (tamaño t_rank x 1) */
    for (PetscInt r = 0; r < t_rank; r++) {
        PetscReal sum = 0.0;
        for (PetscInt k = 0; k < q_rank; k++) {
            sum += A1[r + k * t_rank] * U2[k];
        }
        Bbar[r] = sum;
    }

    /* 9. xhat_end = W_t^T * G(:, P) (tamaño t_rank x 1) */
    for (PetscInt r = 0; r < t_rank; r++) {
        PetscReal sum = 0.0;
        for (PetscInt i = 0; i < N_total; i++) {
            sum += W_t[i + r * N_total] * G[i + P * N_total];
        }
        xhat_end[r] = sum;
    }

    /* 10. Error de reconstrucción DMDc: ||StateData - X_aprox||_F / ||StateData||_F */
    if (dmdc_error) {
        PetscReal *Xhat_all;
        PetscCall(PetscCalloc1(t_rank * N_total, &Xhat_all));

        /* Xhat_all[:, 0] = W_t^T * G[:, 0] */
        for (PetscInt r = 0; r < t_rank; r++) {
            PetscReal sum = 0.0;
            for (PetscInt i = 0; i < N_total; i++) {
                sum += W_t[i + r * N_total] * G[i + 0 * N_total];
            }
            Xhat_all[r + 0 * t_rank] = sum;
        }

        /* Xhat_all[:, i+1] = Hbar * Xhat_all[:, i] + Bbar * inputs[i] */
        for (PetscInt i = 0; i < P; i++) {
            PetscReal u_in = inputs[i];
            for (PetscInt r = 0; r < t_rank; r++) {
                PetscReal sum = 0.0;
                for (PetscInt c = 0; c < t_rank; c++) {
                    sum += Hbar[r + c * t_rank] * Xhat_all[c + i * t_rank];
                }
                Xhat_all[r + (i + 1) * t_rank] = sum + Bbar[r] * u_in;
            }
        }

        /* C = W_t * Xhat_all (N_total x N_total) */
        /* E = I - C */
        PetscReal *E;
        PetscCall(PetscCalloc1(N_total * N_total, &E));
        for (PetscInt r = 0; r < N_total; r++) {
            for (PetscInt c = 0; c < N_total; c++) {
                PetscReal sum = 0.0;
                for (PetscInt k = 0; k < t_rank; k++) {
                    sum += W_t[r + k * N_total] * Xhat_all[k + c * t_rank];
                }
                E[r + c * N_total] = (r == c ? 1.0 : 0.0) - sum;
            }
        }
        PetscCall(PetscFree(Xhat_all));

        /* GE = G * E */
        PetscReal *GE;
        PetscCall(PetscCalloc1(N_total * N_total, &GE));
        for (PetscInt r = 0; r < N_total; r++) {
            for (PetscInt c = 0; c < N_total; c++) {
                PetscReal sum = 0.0;
                for (PetscInt k = 0; k < N_total; k++) {
                    sum += G[r + k * N_total] * E[k + c * N_total];
                }
                GE[r + c * N_total] = sum;
            }
        }

        /* Tr(E^T * GE) = sum_{i} E[i] * GE[i] */
        PetscReal tr_num = 0.0;
        for (PetscInt i = 0; i < N_total * N_total; i++) {
            tr_num += E[i] * GE[i];
        }

        /* Tr(G) = sum_i G[i + i*N_total] */
        PetscReal tr_den = 0.0;
        for (PetscInt i = 0; i < N_total; i++) {
            tr_den += G[i + i * N_total];
        }

        PetscCall(PetscFree(E));
        PetscCall(PetscFree(GE));

        if (tr_den > 1e-14 && tr_num >= 0.0) {
            *dmdc_error = PetscSqrtReal(tr_num) / PetscSqrtReal(tr_den);
        } else {
            *dmdc_error = 0.0;
        }
    }

    /* Liberación de memoria auxiliar */
    PetscCall(PetscFree(G));
    PetscCall(PetscFree(V_S));
    PetscCall(PetscFree(sig_S));
    PetscCall(PetscFree(W_t));
    PetscCall(PetscFree(G_omega));
    PetscCall(PetscFree(V_om));
    PetscCall(PetscFree(sig_om));
    PetscCall(PetscFree(M_q));
    PetscCall(PetscFree(A1));
    PetscCall(PetscFree(A2));
    PetscCall(PetscFree(U2));

    PetscFunctionReturn(PETSC_SUCCESS);
}

/* Función auxiliar para resolver la Ecuación Algebraica de Riccati en Tiempo Continuo (CARE)
   mediante el método del Hamiltoniano (Laub 1979), replicando lqr() de MATLAB.
   H^T P + P H - (1/R_val) P B B^T P + Q = 0, con Q = I.
*/
static PetscErrorCode DLK_SolveCARE(PetscInt d, const PetscReal *H, const PetscReal *B,
                                    PetscReal R_val, PetscReal *P_out, PetscReal *K_out)
{
    PetscBLASInt n2 = 2 * d;
    PetscBLASInt lda = n2, ldvr = n2, info, lwork;
    PetscReal    work_query;
    PetscReal   *Z, *wr, *wi, *vr, *work;
    PetscReal    R_inv = 1.0 / R_val;

    PetscFunctionBeginUser;
    PetscCall(PetscCalloc1(n2 * n2, &Z));
    PetscCall(PetscCalloc1(n2, &wr));
    PetscCall(PetscCalloc1(n2, &wi));
    PetscCall(PetscCalloc1(n2 * n2, &vr));

    /* Z = [  H         -R_inv * B * B^T ]
           [ -Q = -I    -H^T             ]
       Orden por columnas (Fortran/LAPACK)
    */
    for (PetscInt r = 0; r < d; r++) {
        for (PetscInt c = 0; c < d; c++) {
            Z[r + c * n2]             = H[r + c * d];
            Z[r + (c + d) * n2]       = -R_inv * B[r] * B[c];
            Z[(r + d) + c * n2]       = (r == c) ? -1.0 : 0.0;
            Z[(r + d) + (c + d) * n2] = -H[c + r * d];
        }
    }

    lwork = -1;
    LAPACKgeev_("N", "V", &n2, Z, &lda, wr, wi, NULL, &lda, vr, &ldvr, &work_query, &lwork, &info);
    if (info != 0) {
        PetscCall(PetscFree(Z));
        PetscCall(PetscFree(wr));
        PetscCall(PetscFree(wi));
        PetscCall(PetscFree(vr));
        PetscFunctionReturn(PETSC_ERR_LIB);
    }
    lwork = (PetscBLASInt)work_query;
    PetscCall(PetscMalloc1(lwork, &work));

    /* Reconstruir Z porque LAPACKgeev_ destruye la matriz de entrada */
    for (PetscInt r = 0; r < d; r++) {
        for (PetscInt c = 0; c < d; c++) {
            Z[r + c * n2]             = H[r + c * d];
            Z[r + (c + d) * n2]       = -R_inv * B[r] * B[c];
            Z[(r + d) + c * n2]       = (r == c) ? -1.0 : 0.0;
            Z[(r + d) + (c + d) * n2] = -H[c + r * d];
        }
    }

    LAPACKgeev_("N", "V", &n2, Z, &lda, wr, wi, NULL, &lda, vr, &ldvr, work, &lwork, &info);
    PetscCall(PetscFree(work));
    if (info != 0) {
        PetscCall(PetscFree(Z));
        PetscCall(PetscFree(wr));
        PetscCall(PetscFree(wi));
        PetscCall(PetscFree(vr));
        PetscFunctionReturn(PETSC_ERR_LIB);
    }

    /* Identificar y ordenar autovalores estables (Re(lambda) < 0) */
    typedef struct {
        PetscInt  idx;
        PetscReal re;
        PetscReal im;
    } EigPair;
    EigPair eigs[24];
    for (PetscInt i = 0; i < n2; i++) {
        eigs[i].idx = i;
        eigs[i].re  = wr[i];
        eigs[i].im  = wi[i];
    }
    for (PetscInt i = 0; i < n2 - 1; i++) {
        for (PetscInt j = i + 1; j < n2; j++) {
            if (eigs[j].re < eigs[i].re) {
                EigPair tmp = eigs[i];
                eigs[i] = eigs[j];
                eigs[j] = tmp;
            }
        }
    }

    /* Construir base real para el subespacio invariante estable (dimensión n2 x d) */
    PetscReal *V_real;
    PetscCall(PetscCalloc1(n2 * d, &V_real));
    PetscInt col = 0;
    while (col < d) {
        PetscInt k = eigs[col].idx;
        if (PetscAbsReal(wi[k]) < 1e-9) {
            for (PetscInt r = 0; r < n2; r++) V_real[r + col * n2] = vr[r + k * n2];
            col++;
        } else {
            for (PetscInt r = 0; r < n2; r++) {
                V_real[r + col * n2] = vr[r + k * n2];
                if (col + 1 < d) V_real[r + (col + 1) * n2] = vr[r + (k + 1) * n2];
            }
            col += 2;
        }
    }

    PetscReal *U11, *U21;
    PetscCall(PetscCalloc1(d * d, &U11));
    PetscCall(PetscCalloc1(d * d, &U21));
    for (PetscInt r = 0; r < d; r++) {
        for (PetscInt c = 0; c < d; c++) {
            U11[r + c * d] = V_real[r + c * n2];
            U21[r + c * d] = V_real[(r + d) + c * n2];
        }
    }

    PetscReal *A_mat, *B_mat;
    PetscBLASInt bd = (PetscBLASInt)d, binfo;
    PetscBLASInt *ipiv;
    PetscCall(PetscCalloc1(d * d, &A_mat));
    PetscCall(PetscCalloc1(d * d, &B_mat));
    PetscCall(PetscCalloc1(d, &ipiv));

    /* A_mat = U11^T, B_mat = U21^T para resolver U11^T * P = U21^T */
    for (PetscInt r = 0; r < d; r++) {
        for (PetscInt c = 0; c < d; c++) {
            A_mat[r + c * d] = U11[c + r * d];
            B_mat[r + c * d] = U21[c + r * d];
        }
    }

    LAPACKgesv_(&bd, &bd, A_mat, &bd, ipiv, B_mat, &bd, &binfo);
    if (binfo != 0) {
        PetscCall(PetscFree(Z));
        PetscCall(PetscFree(wr));
        PetscCall(PetscFree(wi));
        PetscCall(PetscFree(vr));
        PetscCall(PetscFree(V_real));
        PetscCall(PetscFree(U11));
        PetscCall(PetscFree(U21));
        PetscCall(PetscFree(A_mat));
        PetscCall(PetscFree(B_mat));
        PetscCall(PetscFree(ipiv));
        PetscFunctionReturn(PETSC_ERR_LIB);
    }

    if (P_out) {
        for (PetscInt r = 0; r < d; r++) {
            for (PetscInt c = 0; c < d; c++) {
                P_out[r + c * d] = 0.5 * (B_mat[c + r * d] + B_mat[r + c * d]);
            }
        }
    }

    if (K_out) {
        for (PetscInt c = 0; c < d; c++) {
            PetscReal sum = 0.0;
            for (PetscInt r = 0; r < d; r++) {
                sum += B[r] * B_mat[c + r * d];
            }
            K_out[c] = R_inv * sum;
        }
    }

    PetscCall(PetscFree(Z));
    PetscCall(PetscFree(wr));
    PetscCall(PetscFree(wi));
    PetscCall(PetscFree(vr));
    PetscCall(PetscFree(V_real));
    PetscCall(PetscFree(U11));
    PetscCall(PetscFree(U21));
    PetscCall(PetscFree(A_mat));
    PetscCall(PetscFree(B_mat));
    PetscCall(PetscFree(ipiv));

    PetscFunctionReturn(PETSC_SUCCESS);
}

/* Ley de control LQR basada en CARE con descomposición de controlabilidad (como en MATLAB updatebyLQR.m) */
PetscErrorCode DLK_UpdateByLQR(PetscInt t, const PetscReal *Hbar, const PetscReal *Bbar,
                              const PetscReal *xhat, PetscInt m_last, PetscReal alpha,
                              PetscInt *m_next)
{
    PetscFunctionBeginUser;
    if (alpha <= 0.0) alpha = 1e-4;

    /* 1. Construir matriz de controlabilidad Co = [B, H*B, H^2*B, ..., H^{t-1}*B] (t x t) */
    PetscReal *Co, *v_curr, *v_next;
    PetscCall(PetscCalloc1(t * t, &Co));
    PetscCall(PetscCalloc1(t, &v_curr));
    PetscCall(PetscCalloc1(t, &v_next));

    for (PetscInt i = 0; i < t; i++) v_curr[i] = Bbar[i];

    for (PetscInt k = 0; k < t; k++) {
        for (PetscInt i = 0; i < t; i++) Co[i + k * t] = v_curr[i];
        if (k + 1 < t) {
            for (PetscInt r = 0; r < t; r++) {
                PetscReal sum = 0.0;
                for (PetscInt c = 0; c < t; c++) sum += Hbar[r + c * t] * v_curr[c];
                v_next[r] = sum;
            }
            for (PetscInt i = 0; i < t; i++) v_curr[i] = v_next[i];
        }
    }
    PetscCall(PetscFree(v_curr));
    PetscCall(PetscFree(v_next));

    /* 2. SVD de Co para determinar rango de controlabilidad r */
    PetscBLASInt bt = (PetscBLASInt)t, info, lwork;
    PetscReal    work_q;
    PetscReal   *U_co, *S_co, *VT_co, *work, *Co_copy;
    PetscCall(PetscCalloc1(t * t, &U_co));
    PetscCall(PetscCalloc1(t, &S_co));
    PetscCall(PetscCalloc1(t * t, &VT_co));
    PetscCall(PetscCalloc1(t * t, &Co_copy));
    PetscCall(PetscArraycpy(Co_copy, Co, t * t));

    lwork = -1;
    LAPACKgesvd_("A", "N", &bt, &bt, Co_copy, &bt, S_co, U_co, &bt, VT_co, &bt, &work_q, &lwork, &info);
    PetscInt r_rank = t;
    if (info == 0) {
        lwork = (PetscBLASInt)work_q;
        PetscCall(PetscMalloc1(lwork, &work));
        LAPACKgesvd_("A", "N", &bt, &bt, Co_copy, &bt, S_co, U_co, &bt, VT_co, &bt, work, &lwork, &info);
        PetscCall(PetscFree(work));
        if (info == 0) {
            PetscReal tol_rank = (PetscReal)t * S_co[0] * 1e-12;
            r_rank = 0;
            for (PetscInt i = 0; i < t; i++) {
                if (S_co[i] > tol_rank) r_rank++;
            }
            if (r_rank < 1) r_rank = 1;
        }
    }
    PetscCall(PetscFree(VT_co));
    PetscCall(PetscFree(S_co));
    PetscCall(PetscFree(Co_copy));
    PetscCall(PetscFree(Co));

    PetscInt u_val = m_last;

    if (r_rank < t) {
        /* Subespacio controlable: T = U_co[:, 0:r_rank] (t x r_rank) */
        PetscInt r = r_rank;
        PetscReal *T, *Hc, *Bc, *Kc;
        PetscCall(PetscCalloc1(t * r, &T));
        PetscCall(PetscCalloc1(r * r, &Hc));
        PetscCall(PetscCalloc1(r, &Bc));
        PetscCall(PetscCalloc1(r, &Kc));

        for (PetscInt i = 0; i < t; i++) {
            for (PetscInt j = 0; j < r; j++) {
                T[i + j * t] = U_co[i + j * t];
            }
        }

        /* Hc = T^T * Hbar * T */
        PetscReal *HT;
        PetscCall(PetscCalloc1(t * r, &HT));
        for (PetscInt i = 0; i < t; i++) {
            for (PetscInt j = 0; j < r; j++) {
                PetscReal sum = 0.0;
                for (PetscInt k = 0; k < t; k++) sum += Hbar[i + k * t] * T[k + j * t];
                HT[i + j * t] = sum;
            }
        }
        for (PetscInt i = 0; i < r; i++) {
            for (PetscInt j = 0; j < r; j++) {
                PetscReal sum = 0.0;
                for (PetscInt k = 0; k < t; k++) sum += T[k + i * t] * HT[k + j * t];
                Hc[i + j * r] = sum;
            }
        }
        PetscCall(PetscFree(HT));

        /* Bc = T^T * Bbar */
        for (PetscInt i = 0; i < r; i++) {
            PetscReal sum = 0.0;
            for (PetscInt k = 0; k < t; k++) sum += T[k + i * t] * Bbar[k];
            Bc[i] = sum;
        }

        /* Resolver CARE en subespacio con Rc = alpha consistente con el sistema completo */
        PetscErrorCode err = DLK_SolveCARE(r, Hc, Bc, alpha, NULL, Kc);
        if (!err) {
            /* xhatc = T^T * xhat */
            PetscReal *xhatc;
            PetscCall(PetscCalloc1(r, &xhatc));
            PetscReal norm_xc = 0.0;
            for (PetscInt i = 0; i < r; i++) {
                PetscReal sum = 0.0;
                for (PetscInt k = 0; k < t; k++) sum += T[k + i * t] * xhat[k];
                xhatc[i] = sum;
                norm_xc += sum * sum;
            }
            norm_xc = PetscSqrtReal(norm_xc);
            if (norm_xc > 1e-12) {
                for (PetscInt i = 0; i < r; i++) xhatc[i] /= norm_xc;
            }

            PetscReal u_opt = 0.0;
            for (PetscInt i = 0; i < r; i++) u_opt += -Kc[i] * xhatc[i];
            u_val = (PetscInt)PetscFloorReal(PetscAbsReal(u_opt));
            PetscCall(PetscFree(xhatc));
        }

        PetscCall(PetscFree(T));
        PetscCall(PetscFree(Hc));
        PetscCall(PetscFree(Bc));
        PetscCall(PetscFree(Kc));
    } else {
        /* Sistema completamente controlable: resolver CARE con R = alpha */
        PetscReal *K;
        PetscCall(PetscCalloc1(t, &K));
        PetscErrorCode err = DLK_SolveCARE(t, Hbar, Bbar, alpha, NULL, K);
        if (!err) {
            PetscReal norm_x = 0.0;
            for (PetscInt i = 0; i < t; i++) norm_x += xhat[i] * xhat[i];
            norm_x = PetscSqrtReal(norm_x);
            if (norm_x < 1e-12) norm_x = 1.0;

            PetscReal u_opt = 0.0;
            for (PetscInt i = 0; i < t; i++) u_opt += -K[i] * (xhat[i] / norm_x);
            u_val = (PetscInt)PetscFloorReal(PetscAbsReal(u_opt));
        }
        PetscCall(PetscFree(K));
    }

    PetscCall(PetscFree(U_co));

    if (u_val <= 0) u_val = m_last;
    *m_next = u_val;

    PetscFunctionReturn(PETSC_SUCCESS);
}
