#include "utilities.h"
#include "griddata.h"
#include <stdio.h>
#include "transformation.h"
#include <math.h>

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


double TensorVandermonde(unsigned p,unsigned p_ref,unsigned i, unsigned alpha, GridData *data){
	if(p==p_ref)return data->derivativeVandermondeM[i][alpha];
	return data->VandermondeM[i][alpha];
}

double norm_l2(FunctionPointer f,GridData *data){
	double sum=0;
	double f_val,tx,ty;
	double det_dfe;

	for(unsigned e=0;e<data->E;e++){
		for(unsigned i=0;i<_QUADRATURE_NODES;i++){
			for(int j=0;j<_QUADRATURE_NODES;j++){
				element_transformation(e,data->q_nodes[i],data->q_nodes[j],&tx,&ty,data);
				f_val = f(tx,ty);
				det_dfe=fabs(determinant_jacobian_transformation(e,data->q_nodes[i],data->q_nodes[j],data));
				sum = sum + data->q_weights[i] * data->q_weights[j] * f_val * f_val * det_dfe;
			}
		}
	}
	return sqrt(sum);
}
