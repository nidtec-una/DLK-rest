#include <petscsys.h>
#include <petscblaslapack.h>
#include <stdio.h>
#include <math.h>

PetscErrorCode DLK_SolveCARE_Hamiltonian(PetscInt t, const PetscReal *H, const PetscReal *B,
                                         PetscReal alpha, PetscReal *P_out, PetscReal *K_out)
{
    PetscBLASInt n2 = 2 * t;
    PetscBLASInt lda = n2, ldvr = n2, info, lwork;
    PetscReal    work_query;
    PetscReal   *Z, *wr, *wi, *vr, *work;
    PetscReal    R_inv = 1.0 / alpha;

    PetscFunctionBeginUser;
    PetscCall(PetscCalloc1(n2 * n2, &Z));
    PetscCall(PetscCalloc1(n2, &wr));
    PetscCall(PetscCalloc1(n2, &wi));
    PetscCall(PetscCalloc1(n2 * n2, &vr));

    /* Construir matriz hamiltoniana Z (orden por columnas de Fortran / LAPACK):
       Z = [  H         -R_inv * B * B^T ]
           [ -Q = -I    -H^T             ]
    */
    for (PetscInt r = 0; r < t; r++) {
        for (PetscInt c = 0; c < t; c++) {
            /* Z11 = H */
            Z[r + c * n2] = H[r + c * t];
            /* Z12 = -R_inv * B * B^T */
            Z[r + (c + t) * n2] = -R_inv * B[r] * B[c];
            /* Z21 = -Q (con Q = I) */
            Z[(r + t) + c * n2] = (r == c) ? -1.0 : 0.0;
            /* Z22 = -H^T */
            Z[(r + t) + (c + t) * n2] = -H[c + r * t];
        }
    }

    /* Consulta de tamaño óptimo de trabajo para LAPACKgeev_ */
    lwork = -1;
    LAPACKgeev_("N", "V", &n2, Z, &lda, wr, wi, NULL, &lda, vr, &ldvr, &work_query, &lwork, &info);
    PetscCheck(info == 0, PETSC_COMM_SELF, PETSC_ERR_LIB, "Fallo query LAPACKgeev_ info = %d", (int)info);
    lwork = (PetscBLASInt)work_query;
    PetscCall(PetscMalloc1(lwork, &work));

    /* Reconstruir Z porque LAPACKgeev_ destruye la matriz de entrada */
    for (PetscInt r = 0; r < t; r++) {
        for (PetscInt c = 0; c < t; c++) {
            Z[r + c * n2] = H[r + c * t];
            Z[r + (c + t) * n2] = -R_inv * B[r] * B[c];
            Z[(r + t) + c * n2] = (r == c) ? -1.0 : 0.0;
            Z[(r + t) + (c + t) * n2] = -H[c + r * t];
        }
    }

    LAPACKgeev_("N", "V", &n2, Z, &lda, wr, wi, NULL, &lda, vr, &ldvr, work, &lwork, &info);
    PetscCall(PetscFree(work));
    PetscCheck(info == 0, PETSC_COMM_SELF, PETSC_ERR_LIB, "Fallo ejecución LAPACKgeev_ info = %d", (int)info);

    /* Identificar los t autovalores estables (Re(lambda) < 0) */
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
    /* Ordenar por parte real ascendente */
    for (PetscInt i = 0; i < n2 - 1; i++) {
        for (PetscInt j = i + 1; j < n2; j++) {
            if (eigs[j].re < eigs[i].re) {
                EigPair tmp = eigs[i];
                eigs[i] = eigs[j];
                eigs[j] = tmp;
            }
        }
    }

    /* Construir base real para el subespacio invariante estable (dimension n2 x t) */
    PetscReal *V_real;
    PetscCall(PetscCalloc1(n2 * t, &V_real));

    PetscInt col = 0;
    while (col < t) {
        PetscInt k = eigs[col].idx;
        PetscReal im_val = wr[k]; // check if complex
        if (PetscAbsReal(wi[k]) < 1e-9) {
            /* Autovalor real */
            for (PetscInt r = 0; r < n2; r++) {
                V_real[r + col * n2] = vr[r + k * n2];
            }
            col++;
        } else {
            /* Par conjugado: LAPACK almacena vr[:, k] = real, vr[:, k+1] = imag */
            for (PetscInt r = 0; r < n2; r++) {
                V_real[r + col * n2]       = vr[r + k * n2];
                if (col + 1 < t) {
                    V_real[r + (col + 1) * n2] = vr[r + (k + 1) * n2];
                }
            }
            col += 2;
        }
    }

    /* Extraer bloques U11 = V_real[0:t, :] y U21 = V_real[t:2t, :] */
    PetscReal *U11, *U21;
    PetscCall(PetscCalloc1(t * t, &U11));
    PetscCall(PetscCalloc1(t * t, &U21));
    for (PetscInt r = 0; r < t; r++) {
        for (PetscInt c = 0; c < t; c++) {
            U11[r + c * t] = V_real[r + c * n2];
            U21[r + c * t] = V_real[(r + t) + c * n2];
        }
    }

    /* Resolver P * U11 = U21  =>  U11^T * P^T = U21^T
       En LAPACK, con matrices columna: U11^T * P = U21^T.
       Resolvemos A * X = B donde A = U11^T, B = U21^T (ambas de t x t).
    */
    PetscReal *A_mat, *B_mat;
    PetscBLASInt bt = (PetscBLASInt)t, bnrhs = (PetscBLASInt)t, blda = bt, bldb = bt;
    PetscBLASInt *ipiv;
    PetscCall(PetscCalloc1(t * t, &A_mat));
    PetscCall(PetscCalloc1(t * t, &B_mat));
    PetscCall(PetscCalloc1(t, &ipiv));

    /* A_mat = U11^T (almacenado columna: A_mat[r + c*t] = U11[c + r*t]) */
    for (PetscInt r = 0; r < t; r++) {
        for (PetscInt c = 0; c < t; c++) {
            A_mat[r + c * t] = U11[c + r * t];
            B_mat[r + c * t] = U21[c + r * t];
        }
    }

    LAPACKgesv_(&bt, &bnrhs, A_mat, &blda, ipiv, B_mat, &bldb, &info);
    PetscCheck(info == 0, PETSC_COMM_SELF, PETSC_ERR_LIB, "Fallo en LAPACKgesv_ info = %d", (int)info);

    /* La solución X está en B_mat (X = P^T = P) */
    for (PetscInt r = 0; r < t; r++) {
        for (PetscInt c = 0; c < t; c++) {
            PetscReal val = 0.5 * (B_mat[c + r * t] + B_mat[r + c * t]);
            if (P_out) P_out[r + c * t] = val;
        }
    }

    /* Ganancia K = R_inv * B^T * P (1 x t) */
    if (K_out) {
        for (PetscInt c = 0; c < t; c++) {
            PetscReal sum = 0.0;
            for (PetscInt r = 0; r < t; r++) {
                sum += B[r] * B_mat[c + r * t];
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

int main(int argc, char **argv)
{
    PetscFunctionBeginUser;
    PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));

    PetscInt t = 6;
    PetscReal H[36] = {
        0.5, 0.1, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.4, 0.1, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.3, 0.1, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.2, 0.1, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.1, 0.1,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.05
    };
    PetscReal B[6] = {1.0, 0.5, 0.2, 0.1, 0.05, 0.01};
    PetscReal alpha = 1e-4;
    PetscReal P[36], K[6];

    PetscCall(DLK_SolveCARE_Hamiltonian(t, H, B, alpha, P, K));

    printf("K vector: [");
    for (int i = 0; i < t; i++) printf(" %1.4e", K[i]);
    printf(" ]\n");

    PetscReal normK = 0.0;
    for (int i = 0; i < t; i++) normK += K[i] * K[i];
    normK = sqrt(normK);
    printf("Norm of K: %1.4e\n", normK);

    PetscCall(PetscFinalize());
    return 0;
}
