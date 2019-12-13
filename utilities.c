#include "utilities.h"
#include "griddata.h"
#include <stdio.h>
#include "transformation.h"
#include <math.h>
#include <stdlib.h>



int print_gnuplot_vec_circle(char* filename,Vec v,GridData *data){
	FILE *fp;
	fp = fopen(filename,"w");
	if(fp == NULL) {
		printf("Datei konnte NICHT geoeffnet werden.\n");
		return -1;
	}

	double x,y;
	const double *values;
	VecGetArrayRead(v,&values);
	for(unsigned i=0;i<data->global_dof;i++){
		coordinates_dof(i,&x,&y,data);
		fprintf(fp,"%f %f %f\n",x,y,values[i]);
	}
	VecRestoreArrayRead(v,&values);

	fclose(fp);
	return 0;
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
