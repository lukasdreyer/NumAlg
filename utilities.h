#ifndef UTILITIES_H
#define UTILITIES_H
#include "griddata.h"

typedef double(*FunctionPointer)(double,double);
double radiussquared(double a, double b);
double one(double a, double b);

void print_int_mat(unsigned **mat,unsigned n,unsigned m);
void print_3tensor(double ***mat,unsigned n,unsigned m,unsigned l);
double norm_l2(FunctionPointer f,GridData *data);

double TensorVandermonde(unsigned p,unsigned p_ref,unsigned i, unsigned alpha, GridData *data);


#endif
