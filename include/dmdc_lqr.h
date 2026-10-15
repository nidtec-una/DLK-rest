#ifndef DMDC_LQR_H
#define DMDC_LQR_H

#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>

typedef enum {
    DLK_CTRL_NONE = 0,
    DLK_CTRL_PD   = 1,
    DLK_CTRL_LQR  = 2
} DLKControllerType;

/* Declaración del tipo de controlador en texto para opciones de PETSc */
PETSC_EXTERN const char *const DLKControllerTypes[];

/* Identificación DMDc a partir de los vectores de snapshot */
PETSC_EXTERN PetscErrorCode DLK_DMDc_Identify(PetscInt num_snaps, Vec *snapshots, const PetscReal *inputs,
                                              PetscInt q_rank, PetscInt t_rank,
                                              PetscReal *Hbar, PetscReal *Bbar, PetscReal *xhat_end,
                                              PetscReal *dmdc_error);

/* Ley de control LQR basada en la ecuación algebraica discreta de Riccati (DARE) */
PETSC_EXTERN PetscErrorCode DLK_UpdateByLQR(PetscInt t, const PetscReal *Hbar, const PetscReal *Bbar,
                                            const PetscReal *xhat, PetscInt m_last, PetscReal alpha,
                                            PetscInt *m_next);

/* Ley de control Proporcional-Derivativa (PD) */
PETSC_EXTERN PetscErrorCode DLK_UpdateByPD(PetscReal r_k, PetscReal r_km1, PetscReal r_km2,
                                           PetscInt m_last, PetscInt *m_next);

#endif
