static char help[] = "Driver de prueba y validación para el solver nativo KSPDLK.\n"
                      "Uso: ./test_dlk [-m <nx>] [-n <ny>] [-ksp_type dlk] [-ksp_dlk_controller lqr] ...\n\n";

#include <petscksp.h>
#include "petsckspdlk.h"

int main(int argc, char **argv)
{
    Mat            A;
    Vec            x, b, u;
    KSP            ksp;
    PetscInt       i, j, Ii, J, Istart, Iend, m = 60, n = 60;
    PetscReal      v, hx, hy, norm, tol = 1.e-8;
    PetscInt       its;
    PetscLogDouble t1, t2;
    KSPType        type;
    PetscBool      save_hist = PETSC_FALSE;
    char           hist_file[PETSC_MAX_PATH_LEN] = "";

    PetscFunctionBeginUser;
    PetscCall(PetscInitialize(&argc, &argv, (char *)0, help));

    /* Registrar el solver nativo KSPDLK */
    PetscCall(KSPRegister(KSPDLK, KSPCreate_DLK));

    PetscCall(PetscOptionsGetInt(NULL, NULL, "-m", &m, NULL));
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL));
    PetscCall(PetscOptionsGetString(NULL, NULL, "-history_file", hist_file, sizeof(hist_file), &save_hist));

    hx = 1.0 / (PetscReal)(m + 1);
    hy = 1.0 / (PetscReal)(n + 1);

    /* 1. Crear Matriz A (Laplaciano 2D / Convección-Difusión no simétrica de prueba) */
    PetscCall(MatCreate(PETSC_COMM_WORLD, &A));
    PetscCall(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, m * n, m * n));
    PetscCall(MatSetFromOptions(A));
    PetscCall(MatMPIAIJSetPreallocation(A, 5, NULL, 5, NULL));
    PetscCall(MatSeqAIJSetPreallocation(A, 5, NULL));

    PetscCall(MatGetOwnershipRange(A, &Istart, &Iend));
    for (Ii = Istart; Ii < Iend; Ii++) {
        i = Ii / n;
        j = Ii - i * n;

        /* Vecino inferior */
        if (i > 0) {
            J = Ii - n;
            v = -1.0 / (hy * hy);
            PetscCall(MatSetValues(A, 1, &Ii, 1, &J, &v, INSERT_VALUES));
        }
        /* Vecino superior con convección (rompe simetría para retar a GMRES) */
        if (i < m - 1) {
            J = Ii + n;
            v = -1.0 / (hy * hy) + 10.0 / (2.0 * hy);
            PetscCall(MatSetValues(A, 1, &Ii, 1, &J, &v, INSERT_VALUES));
        }
        /* Vecino izquierdo */
        if (j > 0) {
            J = Ii - 1;
            v = -1.0 / (hx * hx);
            PetscCall(MatSetValues(A, 1, &Ii, 1, &J, &v, INSERT_VALUES));
        }
        /* Vecino derecho con convección */
        if (j < n - 1) {
            J = Ii + 1;
            v = -1.0 / (hx * hx) - 10.0 / (2.0 * hx);
            PetscCall(MatSetValues(A, 1, &Ii, 1, &J, &v, INSERT_VALUES));
        }
        /* Diagonal principal */
        v = 2.0 / (hx * hx) + 2.0 / (hy * hy);
        PetscCall(MatSetValues(A, 1, &Ii, 1, &Ii, &v, INSERT_VALUES));
    }
    PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

    /* 2. Crear Vectores */
    PetscCall(MatCreateVecs(A, &x, &b));
    PetscCall(VecDuplicate(b, &u));
    PetscCall(VecSet(b, 1.0));
    PetscCall(VecSet(x, 0.0));

    /* 3. Configurar KSP */
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetOperators(ksp, A, A));
    PetscCall(KSPSetTolerances(ksp, tol, PETSC_DEFAULT, PETSC_DEFAULT, 2000));

    if (save_hist) {
        PetscCall(KSPSetResidualHistory(ksp, NULL, PETSC_DECIDE, PETSC_TRUE));
    }

    PetscCall(KSPSetFromOptions(ksp));

    /* 4. Resolver */
    PetscCall(PetscTime(&t1));
    PetscCall(KSPSolve(ksp, b, x));
    PetscCall(PetscTime(&t2));

    /* 5. Comprobación del residuo ||b - Ax|| */
    PetscCall(MatMult(A, x, u));
    PetscCall(VecAYPX(u, -1.0, b));
    PetscCall(VecNorm(u, NORM_2, &norm));

    PetscCall(KSPGetIterationNumber(ksp, &its));
    PetscCall(KSPGetType(ksp, &type));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n==================== REPORTE DE EJECUCIÓN KSP ====================\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, " Dimensión del Sistema : N = %d (%d x %d)\n", m * n, m, n));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, " Tipo de Solver (KSP)  : %s\n", type));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, " Iteraciones Totales   : %d\n", its));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, " Residuo ||b - Ax||    : %1.4e\n", (double)norm));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, " Tiempo de Ejecución   : %f segundos\n", t2 - t1));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "==================================================================\n\n"));

    /* Vista detallada del KSP */
    PetscCall(KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD));

    /* 6. Exportar historial CSV si se solicitó */
    if (save_hist) {
        PetscMPIInt rank;
        PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
        if (rank == 0) {
            const PetscReal *hist = NULL;
            PetscInt        num_hist = 0;
            FILE           *fp;

            PetscCall(KSPGetResidualHistory(ksp, &hist, &num_hist));
            fp = fopen(hist_file, "w");
            if (fp) {
                fprintf(fp, "iteration,residual,rel_residual\n");
                if (num_hist > 0 && hist) {
                    PetscReal r0 = hist[0] > 0.0 ? hist[0] : 1.0;
                    for (PetscInt k = 0; k < num_hist; k++) {
                        fprintf(fp, "%d,%1.16e,%1.16e\n", (int)k, (double)hist[k], (double)(hist[k] / r0));
                    }
                }
                fclose(fp);
                PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n Historial de residuos exportado a: %s\n", hist_file));
            }
        }
    }

    /* 7. Limpieza */
    PetscCall(KSPDestroy(&ksp));
    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&u));
    PetscCall(VecDestroy(&b));
    PetscCall(MatDestroy(&A));

    PetscCall(PetscFinalize());
    return 0;
}
