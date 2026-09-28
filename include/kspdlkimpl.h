#ifndef KSPDLKIMPL_H
#define KSPDLKIMPL_H

#include <petsc/private/kspimpl.h>
#include "petsckspdlk.h"

/* Estructura de contexto privado del solver KSPDLK / KSPDLKLGMRES */
typedef struct {
    /* --- Parámetros del espacio de Krylov y reinicio --- */
    PetscInt          m;              /* Parámetro de reinicio actual m_k */
    PetscInt          m_initial;      /* Parámetro inicial m0 (por defecto 30) */
    PetscInt          m_min;          /* Cota inferior de m (por defecto 2) */
    PetscInt          m_max;          /* Cota superior de m (por defecto 100) */
    PetscInt          max_k;          /* Capacidad máxima de memoria reservada para Arnoldi */

    /* --- Parámetros LGMRES (subespacio de enriquecimiento) --- */
    PetscInt          aug_l;          /* Número de vectores de enriquecimiento l (0 = GMRES, >0 = LGMRES) */
    Vec              *aug_vecs;       /* Buffer de vectores de aproximaciones previas z_0, ..., z_{l-1} */
    PetscInt          num_aug;        /* Número actual de vectores en el buffer (0 <= num_aug <= l) */
    Vec               delta_x;        /* Vector auxiliar de desplazamiento de solución */

    /* --- Parámetros de control adaptativo DMDc-LQR / PD --- */
    DLKControllerType ctrl_type;      /* DLK_CTRL_LQR, DLK_CTRL_PD, DLK_CTRL_NONE */
    PetscInt          p_samples;      /* Ventana de snapshots para identificación (def: 5) */
    PetscInt          q_rank;         /* Rango de truncamiento en espacio de entradas (def: 5) */
    PetscInt          t_rank;         /* Rango de truncamiento de estados (def: 6) */
    PetscReal         alpha_lqr;      /* Parámetro de penalización alpha en LQR (def: 1e-4) */

    /* --- Almacenamiento de Snapshots --- */
    Vec              *snapshots;      /* Array de vectores de residuos residuales */
    PetscInt          num_snapshots;  /* Número actual de snapshots recolectados en la ventana */
    PetscInt          max_snapshots;  /* Capacidad máxima del buffer de snapshots */
    PetscReal        *inputs;         /* Registro de los valores m utilizados en cada ciclo */
    PetscReal        *cycle_rnorms;   /* Registro de normas de residuos al final de cada ciclo */
    PetscInt          total_cycles;   /* Contador acumulado de ciclos de reinicio */

    /* --- Base ortonormal de Arnoldi y matrices densas de proyección --- */
    Vec              *V;              /* Base ortonormal de Arnoldi [v_0, v_1, ..., v_{max_k}] */
    PetscScalar      *hh;             /* Matriz de Hessenberg almacenada por columnas */
    PetscReal        *c;              /* Cosenos de las rotaciones de Givens */
    PetscReal        *s;              /* Senos de las rotaciones de Givens */
    PetscScalar      *rs;             /* Vector de residuo proyectado (g = beta * e_1) */
    PetscScalar      *y;              /* Vector de coeficientes solución de mínimos cuadrados */
    Vec               w;              /* Vector de trabajo auxiliar */
} KSP_DLK;

#endif
