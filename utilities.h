#ifndef UTILITIES_H
#define UTILITIES_H
#include "griddata.h"
#include <stdio.h>
#include <stdlib.h>


void print_int_mat(unsigned **mat,unsigned n,unsigned m);
void print_double_mat(double **mat,unsigned n,unsigned m);
void print_3tensor2(double mat[2][2][2]);

void assert(int a);//memory allocation

int print_gnuplot_vec_circle(char* filename,Vec v,GridData *data);



#endif

