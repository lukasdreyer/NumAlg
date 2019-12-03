//static char help[]="Exercise 1."
/*
 ============================================================================
 Name        : ex1.c
 Author      : Lukas, Hannes
 Description : Ex1
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "griddata.h"
#include "utilities.h"

int main(int argc, char **argv) {
	GridData data;

	double L = 0.8;//<sqrt(1/2)
	unsigned M = 1;//number of elements in 1D

	data.q_nodes[0]=-sqrt(1./3);
	data.q_nodes[1]=sqrt(1./3);
	data.q_weights[0]=1;
	data.q_weights[1]=1;


	printf("norm of constant one function\n");
	for(int i=0;i<7;i++){
		init_GridData(M,L,&data);
		printf("error: %.10f\n", fabs(norm_l2(one,&data)-sqrt(M_PI)));
		free_GridData(&data);
		M*=2;
	}

	M=1;
	printf("\nnorm of radius^2\n");

	for(int i=0;i<7;i++){
		init_GridData(M,L,&data);
		printf("error: %.10f\n", fabs(norm_l2(radiussquared,&data)-sqrt(M_PI/3)));
		free_GridData(&data);
		M*=2;
	}

	return 0;
}
