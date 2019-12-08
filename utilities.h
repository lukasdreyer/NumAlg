#ifndef UTILITIES_H
#define UTILITIES_H
#include "griddata.h"



typedef double(*FunctionPointer)(double,double);
double radiussquared(double a, double b);
double one(double a, double b);

void print_int_mat(unsigned **mat,unsigned n,unsigned m);
void print_double_mat(double **mat,unsigned n,unsigned m);
void print_3tensor2(double mat[2][2][2]);
void assert(int a);

#endif
