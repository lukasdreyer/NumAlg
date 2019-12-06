#ifndef MATSHELL_H
#define MATSHELL_H
#include <petscmat.h>


PetscErrorCode mass_mult(Mat A,Vec x,Vec y);
PetscErrorCode stiffness_mult(Mat A,Vec x,Vec y);
PetscErrorCode boundary_mult(Mat A,Vec x, Vec y);

int invertMat2x2(double** M,double** M_inv);
int transposeMat2x2(double** M,double** M_t);
int MatMatMult2x2(double** A,double** B,double** C);
int fill_Mat2x2(double** M,double x00,double x01,double x10,double x11);


#endif
