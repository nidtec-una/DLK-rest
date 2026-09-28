static char help[] = "Driver de prueba y validación para los solvers KSPDLK y KSPDLKLGMRES.\n"
                      "Uso: ./test_dlk [-matrix <archivo.bin>] [-m <nx>] [-n <ny>] [-ksp_type <tipo>] ...\n\n";

#include <petscksp.h>
#include "petsckspdlk.h"

int main(int argc, char **argv)
{
    Mat            A;
    Vec            x, b, u;
    KSP            ksp;
    PetscInt       i, j, Ii, J, Istart, Iend, m = 60, n = 60, total_N = 0;
    PetscReal      v, hx, hy, norm, tol = 1.e-8;
    PetscInt       its;
    PetscLogDouble t1, t2;
    KSPType        type;
    PetscBool      save_hist = PETSC_FALSE, has_matrix = PETSC_FALSE;
    char           hist_file[PETSC_MAX_PATH_LEN] = "";
    char           matrix_file[PETSC_MAX_PATH_LEN] = "";

    PetscFunctionBeginUser;
    PetscCall(PetscInitialize(&argc, &argv, (char *)0, help));

    /* Registrar los solvers nativos KSPDLK (GMRES) y KSPDLKLGMRES (LGMRES) */
    PetscCall(KSPRegister(KSPDLK, KSPCreate_DLK));
    PetscCall(KSPRegister(KSPDLKLGMRES, KSPCreate_DLKLGMRES));

    PetscCall(PetscOptionsGetString(NULL, NULL, "-matrix", matrix_file, sizeof(matrix_file), &has_matrix));
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-m", &m, NULL));
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL));
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-tol", &tol, NULL));
    PetscCall(PetscOptionsGetString(NULL, NULL, "-history_file", hist_file, sizeof(hist_file), &save_hist));

    if (has_matrix) {
        /* Cargar matriz A y vector b desde archivo binario de PETSc */
        PetscViewer viewer;
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Cargando matriz y lado derecho desde: %s ...\n", matrix_file));
        PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, matrix_file, FILE_MODE_READ, &viewer));
        PetscCall(MatCreate(PETSC_COMM_WORLD, &A));
        PetscCall(MatSetFromOptions(A));
        PetscCall(MatLoad(A, viewer));

        PetscCall(VecCreate(PETSC_COMM_WORLD, &b));
        PetscCall(VecSetFromOptions(b));
        PetscCall(VecLoad(b, viewer));
        PetscCall(PetscViewerDestroy(&viewer));

        PetscCall(MatGetSize(A, &total_N, NULL));
    } else {
        /* Generar problema sintético 2D Convección-Difusión */
        total_N = m * n;
        hx = 1.0 / (PetscReal)(m + 1);
        hy = 1.0 / (PetscReal)(n + 1);

        PetscCall(MatCreate(PETSC_COMM_WORLD, &A));
        PetscCall(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, total_N, total_N));
        PetscCall(MatSetFromOptions(A));
        PetscCall(MatMPIAIJSetPreallocation(A, 5, NULL, 5, NULL));
        PetscCall(MatSeqAIJSetPreallocation(A, 5, NULL));

        PetscCall(MatGetOwnershipRange(A, &Istart, &Iend));
        for (Ii = Istart; Ii < Iend; Ii++) {
            i = Ii / n;
            j = Ii - i * n;

            if (i > 0) {
                J = Ii - n;
                v = -1.0 / (hy * hy);
                PetscCall(MatSetValues(A, 1, &Ii, 1, &J, &v, INSERT_VALUES));
            }
            if (i < m - 1) {
                J = Ii + n;
                v = -1.0 / (hy * hy) + 10.0 / (2.0 * hy);
                PetscCall(MatSetValues(A, 1, &Ii, 1, &J, &v, INSERT_VALUES));
            }
            if (j > 0) {
                J = Ii - 1;
                v = -1.0 / (hx * hx);
                PetscCall(MatSetValues(A, 1, &Ii, 1, &J, &v, INSERT_VALUES));
            }
            if (j < n - 1) {
                J = Ii + 1;
                v = -1.0 / (hx * hx) - 10.0 / (2.0 * hx);
                PetscCall(MatSetValues(A, 1, &Ii, 1, &J, &v, INSERT_VALUES));
            }
            v = 2.0 / (hx * hx) + 2.0 / (hy * hy);
            PetscCall(MatSetValues(A, 1, &Ii, 1, &Ii, &v, INSERT_VALUES));
        }
        PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

        PetscCall(MatCreateVecs(A, &x, &b));
        PetscCall(VecSet(b, 1.0));
    }

    /* Crear vector de solución x (inicializado en 0) y vector auxiliar u */
    PetscCall(MatCreateVecs(A, &x, NULL));
    PetscCall(VecDuplicate(b, &u));
    PetscCall(VecSet(x, 0.0));

    /* Configurar KSP */
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetOperators(ksp, A, A));
    PetscCall(KSPSetTolerances(ksp, tol, PETSC_DEFAULT, PETSC_DEFAULT, 2000));

    if (save_hist) {
        PetscCall(KSPSetResidualHistory(ksp, NULL, PETSC_DECIDE, PETSC_TRUE));
    }

    PetscCall(KSPSetFromOptions(ksp));

    /* Resolver */
    PetscCall(PetscTime(&t1));
    PetscCall(KSPSolve(ksp, b, x));
    PetscCall(PetscTime(&t2));

    /* Comprobación del residuo explícito: ||b - Ax|| */
    PetscCall(MatMult(A, x, u));
    PetscCall(VecAYPX(u, -1.0, b));
    PetscCall(VecNorm(u, NORM_2, &norm));

    PetscCall(KSPGetIterationNumber(ksp, &its));
    PetscCall(KSPGetType(ksp, &type));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n==================== REPORTE DE EJECUCIÓN KSP ====================\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, " Dimensión del Sistema : N = %d\n", total_N));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, " Tipo de Solver (KSP)  : %s\n", type));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, " Iteraciones Totales   : %d\n", its));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, " Residuo ||b - Ax||    : %1.4e\n", (double)norm));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, " Tiempo de Ejecución   : %f segundos\n", t2 - t1));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "==================================================================\n\n"));

    /* Exportar historial CSV si se solicitó */
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
                PetscCall(PetscPrintf(PETSC_COMM_WORLD, " Historial de residuos exportado a: %s\n\n", hist_file));
            }
        }
    }

    /* Limpieza */
    PetscCall(KSPDestroy(&ksp));
    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&u));
    PetscCall(VecDestroy(&b));
    PetscCall(MatDestroy(&A));

    PetscCall(PetscFinalize());
    return 0;
}
