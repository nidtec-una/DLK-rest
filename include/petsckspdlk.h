#ifndef PETSCKSPDLK_H
#define PETSCKSPDLK_H

#include <petscksp.h>
#include "dmdc_lqr.h"

#define KSPDLK       "dlk"
#define KSPDLKLGMRES "dlk_lgmres"

/* Funciones de la API pública de PETSc para configurar KSPDLK y KSPDLKLGMRES */
PETSC_EXTERN PetscErrorCode KSPCreate_DLK(KSP ksp);
PETSC_EXTERN PetscErrorCode KSPCreate_DLKLGMRES(KSP ksp);

PETSC_EXTERN PetscErrorCode KSPDLKSetRestart(KSP ksp, PetscInt m);
PETSC_EXTERN PetscErrorCode KSPDLKSetRestartBounds(KSP ksp, PetscInt m_min, PetscInt m_max);
PETSC_EXTERN PetscErrorCode KSPDLKSetController(KSP ksp, DLKControllerType type);
PETSC_EXTERN PetscErrorCode KSPDLKSetAlpha(KSP ksp, PetscReal alpha);
PETSC_EXTERN PetscErrorCode KSPDLKSetRanks(KSP ksp, PetscInt p_samples, PetscInt q_rank, PetscInt t_rank);
PETSC_EXTERN PetscErrorCode KSPDLKSetAugmentation(KSP ksp, PetscInt aug_l);

#endif
