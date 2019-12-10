#include "transformation.h"
#include <math.h>

//F_e = F_s * P_(i,j)
void indices(unsigned e,unsigned *i,unsigned *j,unsigned *s, GridData *data){
	unsigned M=data->M;
	*s=e/(M*M);
	e=e%(M*M);
	*i=e/M;
	*j=e%M;
}
/* Implementation of coefficients needed for transformation */
void coefficients(double x, double y, double *theta, double *c, double *r,GridData *data){
	*theta=M_PI/4*y;
	*c=data->L/(2*cos(*theta));
	*r = (x*(1-*c)+(1+*c))/2;
}

/* 90 degree rotation*/
void rotate(double *x, double *y){
	double temp = *x;
	*x = *y;
	*y = -temp;
}
/* Implement P_ij, scale down x,y and shift to finite element (i,j)*/
void project(double *x,double *y, unsigned i, unsigned j, GridData *data){
	*x = *x/data->M;
	*y = *y/data->M;
	*x = *x + (-1 + (2. * j + 1)/data->M);
	*y = *y + (1 - (2. * i + 1)/data->M);
}
