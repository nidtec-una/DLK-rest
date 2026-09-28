#include "kspdlkimpl.h"
#include <petscblaslapack.h>

/* Setup: Reserva memoria para la base de Arnoldi, matrices densas, snapshots y LGMRES */
static PetscErrorCode KSPSetUp_DLK(KSP ksp)
{
    KSP_DLK *dlk = (KSP_DLK *)ksp->data;

    PetscFunctionBeginUser;
    if (dlk->aug_l < 0) dlk->aug_l = 0;
    if (dlk->m_max < dlk->m_min) dlk->m_max = dlk->m_min + 10;
    if (dlk->max_k < dlk->m_max + dlk->aug_l + 2) dlk->max_k = dlk->m_max + dlk->aug_l + 2;
    dlk->m = dlk->m_initial;

    /* 1. Base ortonormal de Arnoldi V (max_k + 1 vectores) */
    PetscCall(KSPCreateVecs(ksp, dlk->max_k + 1, &dlk->V, 0, NULL));

    /* Vector de trabajo auxiliar */
    PetscCall(VecDuplicate(ksp->vec_rhs, &dlk->w));

    /* 2. Vectores de enriquecimiento para LGMRES */
    if (dlk->aug_l > 0) {
        PetscCall(KSPCreateVecs(ksp, dlk->aug_l, &dlk->aug_vecs, 0, NULL));
        PetscCall(VecDuplicate(ksp->vec_rhs, &dlk->delta_x));
        dlk->num_aug = 0;
    }

    /* 3. Matrices densas para GMRES: Hessenberg, rotaciones de Givens y lados derechos */
    PetscCall(PetscCalloc1((dlk->max_k + 1) * dlk->max_k, &dlk->hh));
    PetscCall(PetscCalloc1(dlk->max_k, &dlk->c));
    PetscCall(PetscCalloc1(dlk->max_k, &dlk->s));
    PetscCall(PetscCalloc1(dlk->max_k + 1, &dlk->rs));
    PetscCall(PetscCalloc1(dlk->max_k, &dlk->y));

    /* 4. Buffer de snapshots para DMDc */
    dlk->max_snapshots = dlk->p_samples + 5;
    PetscCall(KSPCreateVecs(ksp, dlk->max_snapshots, &dlk->snapshots, 0, NULL));
    PetscCall(PetscCalloc1(dlk->max_snapshots, &dlk->inputs));
    PetscCall(PetscCalloc1(dlk->max_snapshots, &dlk->cycle_rnorms));
    dlk->num_snapshots = 0;
    dlk->total_cycles = 0;

    PetscFunctionReturn(PETSC_SUCCESS);
}

/* Reset: Libera memoria dinámica sin destruir el contexto */
static PetscErrorCode KSPReset_DLK(KSP ksp)
{
    KSP_DLK *dlk = (KSP_DLK *)ksp->data;

    PetscFunctionBeginUser;
    if (dlk->V) {
        PetscCall(VecDestroyVecs(dlk->max_k + 1, &dlk->V));
        dlk->V = NULL;
    }
    PetscCall(VecDestroy(&dlk->w));

    if (dlk->aug_vecs) {
        PetscCall(VecDestroyVecs(dlk->aug_l, &dlk->aug_vecs));
        dlk->aug_vecs = NULL;
    }
    PetscCall(VecDestroy(&dlk->delta_x));
    dlk->num_aug = 0;

    PetscCall(PetscFree(dlk->hh));
    PetscCall(PetscFree(dlk->c));
    PetscCall(PetscFree(dlk->s));
    PetscCall(PetscFree(dlk->rs));
    PetscCall(PetscFree(dlk->y));

    if (dlk->snapshots) {
        PetscCall(VecDestroyVecs(dlk->max_snapshots, &dlk->snapshots));
        dlk->snapshots = NULL;
    }
    PetscCall(PetscFree(dlk->inputs));
    PetscCall(PetscFree(dlk->cycle_rnorms));
    dlk->num_snapshots = 0;

    PetscFunctionReturn(PETSC_SUCCESS);
}

/* Destroy: Destruye el objeto KSPDLK */
static PetscErrorCode KSPDestroy_DLK(KSP ksp)
{
    PetscFunctionBeginUser;
    PetscCall(KSPReset_DLK(ksp));
    PetscCall(PetscFree(ksp->data));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/* Solve: Bucle de ciclos de reinicio de GMRES/LGMRES con ajuste adaptativo de m */
static PetscErrorCode KSPSolve_DLK(KSP ksp)
{
    KSP_DLK     *dlk = (KSP_DLK *)ksp->data;
    Mat          A;
    Vec          b, x, w, *V;
    PetscReal    rnorm, beta;
    PetscScalar *hh, *rs, *y;
    PetscReal   *c, *s;
    PetscInt     m, i, j, k;

    PetscFunctionBeginUser;
    PetscCall(PCGetOperators(ksp->pc, &A, NULL));
    b = ksp->vec_rhs;
    x = ksp->vec_sol;
    w = dlk->w;
    V = dlk->V;
    hh = dlk->hh;
    rs = dlk->rs;
    y = dlk->y;
    c = dlk->c;
    s = dlk->s;

    dlk->m = dlk->m_initial;
    dlk->num_snapshots = 0;
    dlk->total_cycles = 0;
    dlk->num_aug = 0;
    ksp->its = 0;

    /* 1. Cálculo de residuo inicial r0 = b - A*x0 */
    if (!ksp->guess_zero) {
        PetscCall(KSP_MatMult(ksp, A, x, w));
        PetscCall(VecAYPX(w, -1.0, b)); /* w = b - A*x */
    } else {
        PetscCall(VecCopy(b, w));
    }
    /* Precondicionamiento por la izquierda: V[0] = M^{-1} * (b - A*x) */
    PetscCall(KSP_PCApply(ksp, w, V[0]));
    PetscCall(VecNorm(V[0], NORM_2, &rnorm));
    KSPCheckNorm(ksp, rnorm);
    ksp->rnorm = rnorm;

    PetscCall(KSPLogResidualHistory(ksp, rnorm));
    PetscCall(KSPMonitor(ksp, 0, rnorm));
    PetscCall((*ksp->converged)(ksp, 0, rnorm, &ksp->reason, ksp->cnvP));
    if (ksp->reason) PetscFunctionReturn(PETSC_SUCCESS);

    /* Guardar snapshot inicial */
    if (dlk->ctrl_type != DLK_CTRL_NONE) {
        PetscCall(VecCopy(V[0], dlk->snapshots[0]));
        if (rnorm > 1e-14) {
            PetscCall(VecScale(dlk->snapshots[0], 1.0 / rnorm));
        }
        dlk->cycle_rnorms[0] = rnorm;
        dlk->inputs[0] = (PetscReal)dlk->m;
        dlk->num_snapshots = 1;
    }

    /* Bucle externo de reinicio */
    while (ksp->its < ksp->max_it && !ksp->reason) {
        dlk->total_cycles++;
        m = dlk->m;
        if (m > dlk->m_max) m = dlk->m_max;
        if (m < dlk->m_min) m = dlk->m_min;

        /* Determinar subespacio de enriquecimiento para LGMRES */
        PetscInt d = (dlk->aug_l > 0) ? dlk->num_aug : 0;
        PetscInt s_dim = m + d;
        if (ksp->its + s_dim > ksp->max_it) {
            s_dim = ksp->max_it - ksp->its;
            if (s_dim < m) m = s_dim;
            d = s_dim - m;
        }

        /* Normalizar v_0 = V[0] / ||V[0]|| */
        beta = rnorm;
        if (beta > 1e-14) {
            PetscCall(VecScale(V[0], 1.0 / beta));
        }
        PetscCall(PetscArrayzero(rs, dlk->max_k + 1));
        rs[0] = beta;

        PetscInt actual_s = 0;

        /* Proceso de Arnoldi de longitud s_dim (m pasos Krylov + d pasos LGMRES) */
        for (j = 0; j < s_dim; j++) {
            actual_s = j + 1;
            ksp->its++;

            if (j < m) {
                /* Paso estándar de Krylov: w = M^{-1} * A * V[j] */
                PetscCall(KSP_MatMult(ksp, A, V[j], w));
                PetscCall(KSP_PCApply(ksp, w, V[j + 1]));
            } else {
                /* Paso de enriquecimiento LGMRES: w = M^{-1} * A * z_{j-m} */
                PetscInt z_idx = j - m;
                PetscCall(KSP_MatMult(ksp, A, dlk->aug_vecs[z_idx], w));
                PetscCall(KSP_PCApply(ksp, w, V[j + 1]));
            }

            /* Ortogonalización Gram-Schmidt Modificada */
            for (i = 0; i <= j; i++) {
                PetscScalar h_val;
                PetscCall(VecDot(V[j + 1], V[i], &h_val));
                hh[i + j * (dlk->max_k + 1)] = h_val;
                PetscCall(VecAXPY(V[j + 1], -h_val, V[i]));
            }
            PetscReal h_next;
            PetscCall(VecNorm(V[j + 1], NORM_2, &h_next));
            hh[(j + 1) + j * (dlk->max_k + 1)] = h_next;

            if (h_next > 1e-14) {
                PetscCall(VecScale(V[j + 1], 1.0 / h_next));
            }

            /* Aplicar rotaciones de Givens previas a la columna j de hh */
            for (i = 0; i < j; i++) {
                PetscScalar h0 = hh[i + j * (dlk->max_k + 1)];
                PetscScalar h1 = hh[(i + 1) + j * (dlk->max_k + 1)];
                hh[i + j * (dlk->max_k + 1)]       = c[i] * h0 + s[i] * h1;
                hh[(i + 1) + j * (dlk->max_k + 1)] = -s[i] * h0 + c[i] * h1;
            }

            /* Calcular nueva rotación de Givens para anular hh[j+1, j] */
            PetscScalar a_diag = hh[j + j * (dlk->max_k + 1)];
            PetscScalar b_sub  = hh[(j + 1) + j * (dlk->max_k + 1)];
            PetscReal   gamma  = PetscSqrtReal(PetscRealPart(a_diag * a_diag + b_sub * b_sub));
            if (gamma > 1e-14) {
                c[j] = PetscRealPart(a_diag) / gamma;
                s[j] = PetscRealPart(b_sub) / gamma;
            } else {
                c[j] = 1.0;
                s[j] = 0.0;
            }
            hh[j + j * (dlk->max_k + 1)] = gamma;
            hh[(j + 1) + j * (dlk->max_k + 1)] = 0.0;

            /* Aplicar rotación al vector del lado derecho rs */
            PetscScalar g0 = rs[j];
            rs[j]     = c[j] * g0;
            rs[j + 1] = -s[j] * g0;

            /* Estimación del residuo actual */
            rnorm = PetscAbsScalar(rs[j + 1]);
            ksp->rnorm = rnorm;

            PetscCall(KSPLogResidualHistory(ksp, rnorm));
            PetscCall(KSPMonitor(ksp, ksp->its, rnorm));
            PetscCall((*ksp->converged)(ksp, ksp->its, rnorm, &ksp->reason, ksp->cnvP));
            if (ksp->reason) break;
        }

        /* Resolver sistema triangular superior R * y = rs (tamaño actual_s) */
        for (i = actual_s - 1; i >= 0; i--) {
            PetscScalar sum = rs[i];
            for (k = i + 1; k < actual_s; k++) {
                sum -= hh[i + k * (dlk->max_k + 1)] * y[k];
            }
            PetscScalar diag = hh[i + i * (dlk->max_k + 1)];
            if (PetscAbsScalar(diag) > 1e-14) {
                y[i] = sum / diag;
            } else {
                y[i] = 0.0;
            }
        }

        /* Construir incremento de solución: delta_x = sum_{i=0}^{m-1} y_i V[i] + sum_{k=0}^{d-1} y_{m+k} z_k */
        if (dlk->aug_l > 0) {
            PetscCall(VecSet(dlk->delta_x, 0.0));
            PetscInt m_used = (actual_s < m) ? actual_s : m;
            for (i = 0; i < m_used; i++) {
                PetscCall(VecAXPY(dlk->delta_x, y[i], V[i]));
            }
            if (actual_s > m) {
                for (k = 0; k < actual_s - m; k++) {
                    PetscCall(VecAXPY(dlk->delta_x, y[m + k], dlk->aug_vecs[k]));
                }
            }
            /* Actualizar solución global: x = x + delta_x */
            PetscCall(VecAXPY(x, 1.0, dlk->delta_x));

            /* Almacenar en buffer circular/FIFO de LGMRES */
            if (dlk->num_aug < dlk->aug_l) {
                PetscCall(VecCopy(dlk->delta_x, dlk->aug_vecs[dlk->num_aug]));
                dlk->num_aug++;
            } else {
                for (PetscInt idx = 1; idx < dlk->aug_l; idx++) {
                    PetscCall(VecCopy(dlk->aug_vecs[idx], dlk->aug_vecs[idx - 1]));
                }
                PetscCall(VecCopy(dlk->delta_x, dlk->aug_vecs[dlk->aug_l - 1]));
            }
        } else {
            /* Modo estándar GMRES */
            for (i = 0; i < actual_s; i++) {
                PetscCall(VecAXPY(x, y[i], V[i]));
            }
        }

        /* Calcular nuevo residuo explícito: r = M^{-1} * (b - A*x) y preparar v_0 */
        PetscCall(KSP_MatMult(ksp, A, x, w));
        PetscCall(VecAYPX(w, -1.0, b));
        PetscCall(KSP_PCApply(ksp, w, V[0]));
        PetscCall(VecNorm(V[0], NORM_2, &rnorm));
        ksp->rnorm = rnorm;

        if (ksp->reason) break;

        /* Almacenar snapshot del ciclo */
        if (dlk->ctrl_type != DLK_CTRL_NONE) {
            PetscInt s_idx = dlk->num_snapshots;
            if (s_idx < dlk->max_snapshots) {
                PetscCall(VecCopy(V[0], dlk->snapshots[s_idx]));
                if (rnorm > 1e-14) {
                    PetscCall(VecScale(dlk->snapshots[s_idx], 1.0 / rnorm));
                }
                dlk->cycle_rnorms[s_idx] = rnorm;
                dlk->inputs[s_idx] = (PetscReal)dlk->m;
                dlk->num_snapshots++;
            }

            /* Control adaptativo para seleccionar el nuevo m */
            PetscInt m_new = dlk->m;
            if (dlk->ctrl_type == DLK_CTRL_PD) {
                if (dlk->num_snapshots >= 3) {
                    PetscInt last = dlk->num_snapshots - 1;
                    PetscCall(DLK_UpdateByPD(dlk->cycle_rnorms[last],
                                            dlk->cycle_rnorms[last - 1],
                                            dlk->cycle_rnorms[last - 2],
                                            dlk->m, &m_new));
                }
            } else if (dlk->ctrl_type == DLK_CTRL_LQR) {
                if (dlk->num_snapshots >= dlk->p_samples + 1) {
                    PetscReal *Hbar, *Bbar, *xhat_end;
                    PetscInt   t = dlk->t_rank;
                    PetscCall(PetscCalloc1(t * t, &Hbar));
                    PetscCall(PetscCalloc1(t, &Bbar));
                    PetscCall(PetscCalloc1(t, &xhat_end));

                    PetscCall(DLK_DMDc_Identify(dlk->p_samples + 1, dlk->snapshots, dlk->inputs,
                                                dlk->q_rank, dlk->t_rank,
                                                Hbar, Bbar, xhat_end));

                    PetscInt dim_sp = dlk->m;
                    PetscCall(DLK_UpdateByLQR(t, Hbar, Bbar, xhat_end, dlk->m, dlk->alpha_lqr, &dim_sp));

                    /* Si es LGMRES, se descuenta l como en LGMRES_DMDc_LQR_fixp.m de arvix */
                    if (dlk->aug_l > 0) {
                        m_new = dim_sp - dlk->aug_l;
                    } else {
                        m_new = dim_sp;
                    }

                    PetscCall(PetscFree(Hbar));
                    PetscCall(PetscFree(Bbar));
                    PetscCall(PetscFree(xhat_end));

                    /* Desplazar ventana de snapshots */
                    PetscCall(VecCopy(dlk->snapshots[dlk->num_snapshots - 1], dlk->snapshots[0]));
                    dlk->cycle_rnorms[0] = dlk->cycle_rnorms[dlk->num_snapshots - 1];
                    dlk->inputs[0] = (PetscReal)m_new;
                    dlk->num_snapshots = 1;
                }
            }

            /* Acotar m_new dentro de [m_min, m_max] */
            if (m_new < dlk->m_min) m_new = dlk->m_min;
            if (m_new > dlk->m_max) m_new = dlk->m_max;
            dlk->m = m_new;
        }
    }

    if (ksp->its >= ksp->max_it && !ksp->reason) {
        ksp->reason = KSP_DIVERGED_ITS;
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}

/* SetFromOptions: Soporte para opciones nativas desde la línea de comandos */
static PetscErrorCode KSPSetFromOptions_DLK(KSP ksp, PetscOptionItems PetscOptionsObject)
{
    KSP_DLK   *dlk = (KSP_DLK *)ksp->data;
    PetscBool  flg;

    PetscFunctionBeginUser;
    PetscOptionsHeadBegin(PetscOptionsObject, "Opciones de KSP DLK (Adaptive DMDc-LQR GMRES / LGMRES)");
    PetscCall(PetscOptionsInt("-ksp_dlk_m0", "Parámetro inicial de reinicio m0", "KSPDLKSetRestart", dlk->m_initial, &dlk->m_initial, &flg));
    if (flg) dlk->m = dlk->m_initial;
    PetscCall(PetscOptionsInt("-ksp_dlk_m_min", "Cota mínima de m", "KSPDLKSetRestartBounds", dlk->m_min, &dlk->m_min, NULL));
    PetscCall(PetscOptionsInt("-ksp_dlk_m_max", "Cota máxima de m", "KSPDLKSetRestartBounds", dlk->m_max, &dlk->m_max, NULL));
    PetscCall(PetscOptionsInt("-ksp_dlk_aug", "Número de vectores de enriquecimiento LGMRES (l)", "KSPDLKSetAugmentation", dlk->aug_l, &dlk->aug_l, NULL));
    PetscCall(PetscOptionsEnum("-ksp_dlk_controller", "Tipo de controlador adaptativo", "KSPDLKSetController", DLKControllerTypes, (PetscEnum)dlk->ctrl_type, (PetscEnum *)&dlk->ctrl_type, NULL));
    PetscCall(PetscOptionsReal("-ksp_dlk_alpha", "Penalización alpha para LQR", "KSPDLKSetAlpha", dlk->alpha_lqr, &dlk->alpha_lqr, NULL));
    PetscCall(PetscOptionsInt("-ksp_dlk_p", "Ventana de snapshots p", "KSPDLKSetRanks", dlk->p_samples, &dlk->p_samples, NULL));
    PetscCall(PetscOptionsInt("-ksp_dlk_q", "Rango q para entradas DMDc", "KSPDLKSetRanks", dlk->q_rank, &dlk->q_rank, NULL));
    PetscCall(PetscOptionsInt("-ksp_dlk_t", "Rango t para estados DMDc", "KSPDLKSetRanks", dlk->t_rank, &dlk->t_rank, NULL));
    PetscOptionsHeadEnd();

    PetscFunctionReturn(PETSC_SUCCESS);
}

/* View: Información en -ksp_view */
static PetscErrorCode KSPView_DLK(KSP ksp, PetscViewer viewer)
{
    KSP_DLK   *dlk = (KSP_DLK *)ksp->data;
    PetscBool  isascii;

    PetscFunctionBeginUser;
    PetscCall(PetscObjectTypeCompare((PetscObject)viewer, PETSCVIEWERASCII, &isascii));
    if (isascii) {
        if (dlk->aug_l > 0) {
            PetscCall(PetscViewerASCIIPrintf(viewer, "  DLK: Variante = LGMRES con l = %d vectores de enriquecimiento\n", dlk->aug_l));
        } else {
            PetscCall(PetscViewerASCIIPrintf(viewer, "  DLK: Variante = GMRES estándar\n"));
        }
        PetscCall(PetscViewerASCIIPrintf(viewer, "  DLK: Controlador = %s\n", DLKControllerTypes[dlk->ctrl_type]));
        PetscCall(PetscViewerASCIIPrintf(viewer, "  DLK: m_inicial = %d, m_min = %d, m_max = %d\n", dlk->m_initial, dlk->m_min, dlk->m_max));
        if (dlk->ctrl_type == DLK_CTRL_LQR) {
            PetscCall(PetscViewerASCIIPrintf(viewer, "  DLK: alpha_lqr = %g, p_samples = %d, q_rank = %d, t_rank = %d\n", (double)dlk->alpha_lqr, dlk->p_samples, dlk->q_rank, dlk->t_rank));
        }
    }
    PetscFunctionReturn(PETSC_SUCCESS);
}

/* Constructor KSPCreate_DLK (variante GMRES por defecto) */
PETSC_EXTERN PetscErrorCode KSPCreate_DLK(KSP ksp)
{
    KSP_DLK *dlk;

    PetscFunctionBeginUser;
    PetscCall(PetscNew(&dlk));
    ksp->data = (void *)dlk;

    /* Valores por defecto */
    dlk->m_initial     = 30;
    dlk->m             = 30;
    dlk->m_min         = 2;
    dlk->m_max         = 60;
    dlk->max_k         = 70;
    dlk->aug_l         = 0;  /* 0 = GMRES estándar */
    dlk->ctrl_type     = DLK_CTRL_LQR;
    dlk->p_samples     = 5;
    dlk->q_rank        = 5;
    dlk->t_rank        = 6;
    dlk->alpha_lqr     = 1e-4;

    PetscCall(KSPSetSupportedNorm(ksp, KSP_NORM_PRECONDITIONED, PC_LEFT, 2));
    PetscCall(KSPSetSupportedNorm(ksp, KSP_NORM_UNPRECONDITIONED, PC_LEFT, 1));

    ksp->ops->setup          = KSPSetUp_DLK;
    ksp->ops->solve          = KSPSolve_DLK;
    ksp->ops->reset          = KSPReset_DLK;
    ksp->ops->destroy        = KSPDestroy_DLK;
    ksp->ops->setfromoptions = KSPSetFromOptions_DLK;
    ksp->ops->view           = KSPView_DLK;
    ksp->ops->buildsolution  = KSPBuildSolutionDefault;
    ksp->ops->buildresidual  = KSPBuildResidualDefault;

    PetscFunctionReturn(PETSC_SUCCESS);
}

/* Constructor KSPCreate_DLKLGMRES (variante LGMRES con l=2 por defecto) */
PETSC_EXTERN PetscErrorCode KSPCreate_DLKLGMRES(KSP ksp)
{
    PetscFunctionBeginUser;
    PetscCall(KSPCreate_DLK(ksp));
    KSP_DLK *dlk = (KSP_DLK *)ksp->data;
    dlk->aug_l = 2; /* 2 vectores de enriquecimiento por defecto */
    PetscFunctionReturn(PETSC_SUCCESS);
}

/* Funciones de la API de configuración en C */
PetscErrorCode KSPDLKSetRestart(KSP ksp, PetscInt m)
{
    KSP_DLK *dlk = (KSP_DLK *)ksp->data;
    PetscFunctionBeginUser;
    dlk->m_initial = m;
    dlk->m = m;
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode KSPDLKSetRestartBounds(KSP ksp, PetscInt m_min, PetscInt m_max)
{
    KSP_DLK *dlk = (KSP_DLK *)ksp->data;
    PetscFunctionBeginUser;
    dlk->m_min = m_min;
    dlk->m_max = m_max;
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode KSPDLKSetController(KSP ksp, DLKControllerType type)
{
    KSP_DLK *dlk = (KSP_DLK *)ksp->data;
    PetscFunctionBeginUser;
    dlk->ctrl_type = type;
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode KSPDLKSetAlpha(KSP ksp, PetscReal alpha)
{
    KSP_DLK *dlk = (KSP_DLK *)ksp->data;
    PetscFunctionBeginUser;
    dlk->alpha_lqr = alpha;
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode KSPDLKSetRanks(KSP ksp, PetscInt p_samples, PetscInt q_rank, PetscInt t_rank)
{
    KSP_DLK *dlk = (KSP_DLK *)ksp->data;
    PetscFunctionBeginUser;
    dlk->p_samples = p_samples;
    dlk->q_rank    = q_rank;
    dlk->t_rank    = t_rank;
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode KSPDLKSetAugmentation(KSP ksp, PetscInt aug_l)
{
    KSP_DLK *dlk = (KSP_DLK *)ksp->data;
    PetscFunctionBeginUser;
    dlk->aug_l = aug_l;
    PetscFunctionReturn(PETSC_SUCCESS);
}
