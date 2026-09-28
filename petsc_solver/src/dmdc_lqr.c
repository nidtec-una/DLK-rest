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
                                 PetscReal *Hbar, PetscReal *Bbar, PetscReal *xhat_end)
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

/* Ley de control LQR basada en iteración de punto fijo para DARE */
PetscErrorCode DLK_UpdateByLQR(PetscInt t, const PetscReal *Hbar, const PetscReal *Bbar,
                              const PetscReal *xhat, PetscInt m_last, PetscReal alpha,
                              PetscInt *m_next)
{
    PetscReal *P, *P_next, *Q, R_lqr;
    PetscReal *Ht_P, *temp;
    PetscInt  max_iter = 100;
    PetscReal tol = 1e-7;

    PetscFunctionBeginUser;
    if (alpha <= 0.0) alpha = 1e-4;
    R_lqr = alpha;

    PetscCall(PetscCalloc1(t * t, &P));
    PetscCall(PetscCalloc1(t * t, &P_next));
    PetscCall(PetscCalloc1(t * t, &Q));
    PetscCall(PetscCalloc1(t * t, &Ht_P));
    PetscCall(PetscCalloc1(t, &temp));

    /* Q = (1 - alpha) * I */
    for (PetscInt i = 0; i < t; i++) {
        Q[i + i * t] = 1.0 - alpha;
        P[i + i * t] = Q[i + i * t];
    }

    /* Iteración de Riccati: P_{k+1} = H^T P H - (H^T P B)(B^T P H) / (R + B^T P B) + Q */
    for (PetscInt iter = 0; iter < max_iter; iter++) {
        /* Ht_P = H^T * P */
        for (PetscInt r = 0; r < t; r++) {
            for (PetscInt c = 0; c < t; c++) {
                PetscReal sum = 0.0;
                for (PetscInt k = 0; k < t; k++) {
                    sum += Hbar[k + r * t] * P[k + c * t];
                }
                Ht_P[r + c * t] = sum;
            }
        }

        /* Ht_P_B = (H^T P) * B (vector t x 1) */
        for (PetscInt r = 0; r < t; r++) {
            PetscReal sum = 0.0;
            for (PetscInt k = 0; k < t; k++) {
                sum += Ht_P[r + k * t] * Bbar[k];
            }
            temp[r] = sum;
        }

        /* denom = R + B^T P B */
        PetscReal Bt_P_B = 0.0;
        for (PetscInt r = 0; r < t; r++) {
            for (PetscInt c = 0; c < t; c++) {
                Bt_P_B += Bbar[r] * P[r + c * t] * Bbar[c];
            }
        }
        PetscReal denom = R_lqr + Bt_P_B;
        if (PetscAbsReal(denom) < 1e-12) denom = 1e-12;

        /* P_next = H^T P H - (Ht_P_B * Ht_P_B^T) / denom + Q */
        PetscReal diff_max = 0.0;
        for (PetscInt r = 0; r < t; r++) {
            for (PetscInt c = 0; c < t; c++) {
                PetscReal Ht_P_H = 0.0;
                for (PetscInt k = 0; k < t; k++) {
                    Ht_P_H += Ht_P[r + k * t] * Hbar[k + c * t];
                }
                PetscReal val = Ht_P_H - (temp[r] * temp[c]) / denom + Q[r + c * t];
                PetscReal diff = PetscAbsReal(val - P[r + c * t]);
                if (diff > diff_max) diff_max = diff;
                P_next[r + c * t] = val;
            }
        }

        PetscCall(PetscArraycpy(P, P_next, t * t));
        if (diff_max < tol) break;
    }

    /* Ganancia K = (B^T P H) / denom */
    /* Bt_P_H = (H^T P B)^T = temp^T */
    PetscReal denom = R_lqr;
    for (PetscInt r = 0; r < t; r++) {
        for (PetscInt c = 0; c < t; c++) {
            denom += Bbar[r] * P[r + c * t] * Bbar[c];
        }
    }
    if (PetscAbsReal(denom) < 1e-12) denom = 1e-12;

    /* u_opt = -K * xhat_norm */
    PetscReal norm_x = 0.0;
    for (PetscInt i = 0; i < t; i++) norm_x += xhat[i] * xhat[i];
    norm_x = PetscSqrtReal(norm_x);
    if (norm_x < 1e-12) norm_x = 1.0;

    PetscReal u_opt = 0.0;
    for (PetscInt i = 0; i < t; i++) {
        PetscReal Ki = temp[i] / denom;
        u_opt += -Ki * (xhat[i] / norm_x);
    }

    PetscInt u_val = (PetscInt)PetscFloorReal(PetscAbsReal(u_opt));
    if (u_val <= 0) u_val = m_last;
    *m_next = u_val;

    PetscCall(PetscFree(P));
    PetscCall(PetscFree(P_next));
    PetscCall(PetscFree(Q));
    PetscCall(PetscFree(Ht_P));
    PetscCall(PetscFree(temp));

    PetscFunctionReturn(PETSC_SUCCESS);
}
