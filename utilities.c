#include "utilities.h"
#include "griddata.h"
#include <stdio.h>
#include "transformation.h"
#include <math.h>
#include <stdlib.h>

double radiussquared(double a, double b){
	return a*a+b*b;
}
double one(double a, double b){
	return 1;
}

void print_int_mat(unsigned **mat,unsigned n,unsigned m){
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			printf("%i ",mat[i][j]);
		}
		printf("\n");
	}
}
void print_3tensor(double ***mat,unsigned n,unsigned m,unsigned l){
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			for(int k=0;k<l;k++){
				printf("%f ",mat[i][j][k]);
			}
			printf("\n");
		}
		printf("\n");
		printf("\n");
	}
}
