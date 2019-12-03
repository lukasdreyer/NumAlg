#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H

#include "griddata.h"

void indices(unsigned e,unsigned *i,unsigned *j,unsigned *s, GridData *data);
void rotate(double *x, double *y);
void project(double *x,double *y, unsigned i, unsigned j, GridData *data);
void coefficients(double x, double y, double *theta, double *c, double *r,GridData *data);
void element_transformation(unsigned e,double x, double y,double *tx,double *ty, GridData *data);
double determinant_jacobian_transformation(unsigned e,double x, double y, GridData *data);

#endif