#ifndef UTILITIES_H
#define UTILITIES_H
#include "griddata.h"



typedef double(*FunctionPointer)(double,double);
double radiussquared(double a, double b);
double one(double a, double b);

void print_int_mat(unsigned **mat,unsigned n,unsigned m);
void print_3tensor(double ***mat,unsigned n,unsigned m,unsigned l);


#endif
