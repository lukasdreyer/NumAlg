#include <petscmat.h>
#include <petscsys.h>
#include "griddata.h"
#include "transformation.h"
#include "utilities.h"
#include "matshell.h"

PetscErrorCode mass_mult(Mat A,Vec x,Vec y){
	GridData *data;

	const double *inputvector;
	double *outputvector;
	double x_local[_DOF1D][_DOF1D];

	double tensorM1[_DOF1D][_QUADRATURE_NODES];
	double tensorM2[_QUADRATURE_NODES][_QUADRATURE_NODES];
	double tensorM3[_QUADRATURE_NODES][_DOF1D];

	VecGetArrayRead(x,&inputvector);
	VecGetArray(y,&outputvector);
	MatShellGetContext(A,&data);

	for(unsigned dof=0;dof<data->global_dof;dof++){
		outputvector[dof]=0;
	}

	for(unsigned e=0;e<data->E;e++){

//		printf("element: %i",e);
		for(unsigned i=0;i<_DOF1D;i++){
			for(unsigned j=0;j<_DOF1D;j++){
				x_local[i][j]=inputvector[data->FEtoDOF[e][i*_DOF1D+j]];//symmetry
			}
		}
//		printf(", 1");
		for(unsigned k = 0;k < _DOF1D;k++){
			for(unsigned beta = 0;beta < _QUADRATURE_NODES;beta++){
				tensorM1[k][beta] =0;
				for(unsigned l = 0; l < _DOF1D; l++){
					tensorM1[k][beta] += data->VandermondeM[beta][l]*x_local[k][l];
				}
			}
		}
//		printf(", 2");

		for(unsigned alpha = 0;alpha < _QUADRATURE_NODES;alpha++){
			for(unsigned beta = 0;beta < _QUADRATURE_NODES;beta++){
				tensorM2[alpha][beta] =0;
				for(unsigned k = 0; k < _DOF1D; k++){
					tensorM2[alpha][beta] += data->VandermondeM[alpha][k]*tensorM1[k][beta];
				}
			}
		}
//		printf(", 3");

		for(unsigned alpha = 0;alpha < _QUADRATURE_NODES;alpha++){
			for(unsigned j = 0;j < _DOF1D;j++){
				tensorM3[alpha][j] =0;
				for(unsigned beta = 0; beta < _QUADRATURE_NODES; beta++){
//					printf(",a: %i,b: %i ",alpha,beta);

					tensorM3[alpha][j] += data->VandermondeM[beta][j] * tensorM2[alpha][beta]
											* data->q_weights[alpha]*data->q_weights[beta]
											* data->det_DJe[alpha][beta][e];
				}
			}
		}
//		printf(", 4");

		for(unsigned i = 0;i < _DOF1D;i++){
			for(unsigned j = 0;j < _DOF1D;j++){
				for(unsigned alpha = 0; alpha < _QUADRATURE_NODES; alpha++){
					outputvector[data->FEtoDOF[e][i*_DOF1D+j]] += data->VandermondeM[alpha][i] * tensorM3[alpha][j];
//					outputvector[0] += data->VandermondeM[alpha][i] * tensorM3[alpha][j];
				}
			}
		}
//		printf(", 5\n");
	}//element loop
	VecRestoreArrayRead(x,&inputvector);
	VecRestoreArray(y,&outputvector);
	return 0;
}

PetscErrorCode stiffness_mult(Mat A,Vec x,Vec y){
	GridData *data;

	const double *inputvector;
	double *outputvector;
	double x_local[_DOF1D][_DOF1D];

	double tensor_vandermonde_value;

	double tensorA1[_DIM][_QUADRATURE_NODES][_DOF1D];
	double tensorA2[_DIM][_QUADRATURE_NODES][_QUADRATURE_NODES];
	double tensorA3[_DIM][_QUADRATURE_NODES][_QUADRATURE_NODES];
	double tensorA4[_DIM][_DOF1D][_QUADRATURE_NODES];
	double tensorA5[_DOF1D][_DOF1D][_QUADRATURE_NODES];

	VecGetArrayRead(x,&inputvector);
	VecGetArray(y,&outputvector);
	MatShellGetContext(A,&data);

	for(unsigned dof=0;dof<data->global_dof;dof++){
		outputvector[dof]=0;
	}

	unsigned s_idx,i_idx,j_idx;//e=(s,i,j), where s determines which coarse element, i is row of refined

	for(unsigned e=0;e<data->E;e++){
		indices(e,&i_idx,&j_idx,&s_idx,data);//		DPij = DP, DFs , Qs-1

		for(unsigned k=0;k<_DOF1D;k++){
			for(unsigned l=0;l<_DOF1D;l++){
				x_local[k][l]=inputvector[data->FEtoDOF[e][k*_DOF1D+l]];//symmetry
			}
		}
		for(unsigned p =0; p < _DIM ; p++){
			for(unsigned  beta=0;beta < _QUADRATURE_NODES ;  beta++){
				for(unsigned  k=0;k < _DOF1D ;  k++){
					tensorA1[p][beta][k]=0;
					for(unsigned  l=0; l< _DOF1D ;  l++){
						if(p==1) tensor_vandermonde_value = data->derivativeVandermondeM[beta][l];
						else tensor_vandermonde_value = data->VandermondeM[beta][l];
						tensorA1[p][beta][k]+=tensor_vandermonde_value*x_local[k][l];
					}
				}
			}
		}

		for(unsigned  p=0;p < _DIM ; p ++){
			for(unsigned beta =0; beta < _QUADRATURE_NODES ;  beta++){
				for(unsigned  alpha=0; alpha< _QUADRATURE_NODES ; alpha ++){
					tensorA2[p][beta][alpha]=0;
					for(unsigned  k=0;k < _DOF1D ; k ++){
						if(p==0) tensor_vandermonde_value = data->derivativeVandermondeM[alpha][k];
						else tensor_vandermonde_value = data->VandermondeM[alpha][k];
						tensorA2[p][beta][alpha] += tensor_vandermonde_value*tensorA1[p][beta][k];
					}
				}
			}
		}

		for(unsigned  q=0; q<  _DIM;  q++){
			for(unsigned  beta=0; beta< _QUADRATURE_NODES ;  beta++){
				for(unsigned  alpha=0; alpha< _QUADRATURE_NODES ;  alpha++){
					tensorA3[1][beta][alpha]=0;
					for(unsigned  p=0; p< _DOF1D ;  p++){
						tensorA3[q][beta][alpha] += tensorA2[p][beta][alpha]*data->Gepq[alpha][beta][e][p][q];
					}
				}
			}
		}


		for(unsigned  q=0; q<  _DIM;  q++){
			for(unsigned  j=0; j< _DOF1D ;  j++){
				for(unsigned  alpha=0; alpha< _QUADRATURE_NODES ;  alpha++){
					tensorA4[q][j][alpha]=0;
					for(unsigned  beta=0; beta< _DOF1D ;  beta++){
						if(q==1) tensor_vandermonde_value = data->derivativeVandermondeM[beta][j];
						else tensor_vandermonde_value = data->VandermondeM[beta][j];
						tensorA4[q][j][alpha] += tensorA3[q][beta][alpha] * tensor_vandermonde_value
										* data->det_DJe[alpha][beta][e]* data->q_weights[alpha]*data->q_weights[beta];
					}
				}
			}
		}

		for(unsigned  i=0; i< _DOF1D ;  i++){
			for(unsigned  j=0; j< _DOF1D ;  j++){
				for(unsigned  alpha=0; alpha< _QUADRATURE_NODES ;  alpha++){
					tensorA5[i][j][alpha]=0;
					for(unsigned  q=0; q< _DIM ;  q++){
						if(q==0) tensor_vandermonde_value = data->derivativeVandermondeM[alpha][i];
						else tensor_vandermonde_value = data->VandermondeM[alpha][i];
						tensorA5[i][j][alpha]+=tensor_vandermonde_value*tensorA4[q][j][alpha];
					}
				}
			}
		}

		for(unsigned i =0; i< _DOF1D ;  i++){
			for(unsigned j =0; j< _DOF1D ;  j++){
				for(unsigned  alpha=0; alpha< _QUADRATURE_NODES ;  alpha++){
					outputvector[data->FEtoDOF[e][i*_DOF1D+j]] += tensorA5[i][j][alpha];
				}
			}
		}

	}//end element loop

	VecRestoreArrayRead(x,&inputvector);
	VecRestoreArray(y,&outputvector);

	return 0;
}

PetscErrorCode boundary_mult(Mat A,Vec x, Vec y){
	GridData *data;
	const double *inputvector;
	double *outputvector;


	MatShellGetContext(A,&data);
	MatMult(data->StiffnessM,x,y);

	VecGetArrayRead(x,&inputvector);
	VecGetArray(y,&outputvector);

	for(unsigned i=0;i<data->boundary_dof;i++){
		outputvector[data->boundary_nodes[i]]=inputvector[data->boundary_nodes[i]];
	}

	VecRestoreArrayRead(x,&inputvector);
	VecRestoreArray(y,&outputvector);
	return 0;
}

int invertMat2x2(double** M,double** M_inv){
	double det=M[0][0]*M[1][1]-M[0][1]*M[1][0];
	M_inv[0][0]=M[1][1]/det;
	M_inv[0][1]=-M[0][1]/det;
	M_inv[1][0]=-M[1][0]/det;
	M_inv[1][1]=M[0][0]/det;
	return 0;
}
int transposeMat2x2(double** M,double** M_t){
	M_t[0][0]=M[0][0];
	M_t[0][1]=M[1][0];
	M_t[1][0]=M[0][1];
	M_t[1][1]=M[1][1];
	return 0;
}
int MatMatMult2x2(double** A,double** B,double** C){
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


