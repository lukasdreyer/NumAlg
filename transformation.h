#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H

#include "griddata.h"

//ATTENTION: i always corresponds to y coordinate, j to x coordinate!!!
//for indices, coordinates, project


void indices(unsigned e,unsigned *i,unsigned *j,unsigned *s, GridData *data);
void coordinates_dof(unsigned i, double *x,double *y, GridData *data);
void project(double *x,double *y, unsigned i, unsigned j, GridData *data);


void rotate(double *x, double *y);
void coefficients(double x, double y, double *theta, double *c, double *r,GridData *data);


#endif
