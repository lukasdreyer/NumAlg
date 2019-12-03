#ifndef MATSHELL_H
#define MATSHELL_H
#include <petscmat.h>

PetscErrorCode mass_mult(Mat A,Vec x,Vec y);
PetscErrorCode stiffness_mult(Mat A,Vec x,Vec y);
PetscErrorCode boundary_mult(Mat A,Vec x, Vec y);

#endif