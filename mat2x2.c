#include "utilities.h"
#include "mat2x2.h"
int invert_Mat2x2(double** M,double** M_inv){
	double det=determinant_Mat2x2(M);
	if(det==0)printf("det=0");
	M_inv[0][0]=M[1][1]/det;
	M_inv[0][1]=-M[0][1]/det;
	M_inv[1][0]=-M[1][0]/det;
	M_inv[1][1]=M[0][0]/det;
	return 0;
}
int transpose_Mat2x2(double** M,double** M_t){
	M_t[0][0]=M[0][0];
	M_t[0][1]=M[1][0];
	M_t[1][0]=M[0][1];
	M_t[1][1]=M[1][1];
	return 0;
}
int MatMult_Mat2x2(double** A,double** B,double** C){
	C[0][0]=A[0][0]*B[0][0]+A[0][1]*B[1][0];
	C[0][1]=A[0][0]*B[0][1]+A[0][1]*B[1][1];
	C[1][0]=A[1][0]*B[0][0]+A[1][1]*B[1][0];
	C[1][1]=A[1][0]*B[0][1]+A[1][1]*B[1][1];
	return 0;
}

int fill_Mat2x2(double** M,double x00,double x01,double x10,double x11){
	M[0][0]=x00;
	M[0][1]=x01;
	M[1][0]=x10;
	M[1][1]=x11;
	return 0;
}
int alloc_Mat2x2(double ***M_p){
	*M_p=malloc(sizeof(double*) * 2);
	for(unsigned i=0;i<2;i++){
		(*M_p)[i]=malloc(sizeof(double)*2);
	}
	return 0;
}
int free_Mat2x2(double **M){
	for(unsigned i=0;i<2;i++){
		free(M[i]);
	}
	free(M);
	return 0;
}
int differ_Mat2x2(double **M1,double **M2){
	for(unsigned i=0;i<1;i++){
		for(unsigned j=0;j<1;j++){
			if(M1[i][j]!=M2[i][j]){
				printf("The matrices differ! \n");
				return i*2+j;
			}
		}
	}
	return 0;
}

int print_Mat2x2(double **M){
	print_double_mat(M,2,2);
	return 0;
}
double determinant_Mat2x2(double **M){
	return M[0][0]*M[1][1]-M[1][0]*M[0][1];
}

