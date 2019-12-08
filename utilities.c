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
	printf("\n");
}
void print_double_mat(double **mat,unsigned n,unsigned m){
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			printf("%f ",mat[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}
void print_3tensor2(double mat[2][2][2]){
	for(int i=0;i<2;i++){
		for(int j=0;j<2;j++){
			for(int k=0;k<2;k++){
				printf("%f ",mat[i][j][k]);
			}
			printf("\n");
		}
		printf("\n");
		printf("\n");
	}
}
void assert(int a){
	if(!a)printf("allocation error\n");
}
