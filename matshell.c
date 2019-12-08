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

		for(unsigned i=0;i<_DOF1D;i++){
			for(unsigned j=0;j<_DOF1D;j++){
				x_local[i][j]=inputvector[data->FEtoDOF[e][i*_DOF1D+j]];//symmetry
			}
		}

		for(unsigned k = 0;k < _DOF1D;k++){
			for(unsigned beta = 0;beta < _QUADRATURE_NODES;beta++){
				tensorM1[k][beta] =0;
				for(unsigned l = 0; l < _DOF1D; l++){
					tensorM1[k][beta] += data->VandermondeM[beta][l]*x_local[k][l];
				}
			}
		}

		for(unsigned alpha = 0;alpha < _QUADRATURE_NODES;alpha++){
			for(unsigned beta = 0;beta < _QUADRATURE_NODES;beta++){
				tensorM2[alpha][beta] =0;
				for(unsigned k = 0; k < _DOF1D; k++){
					tensorM2[alpha][beta] += data->VandermondeM[alpha][k]*tensorM1[k][beta];
				}
			}
		}

		for(unsigned alpha = 0;alpha < _QUADRATURE_NODES;alpha++){
			for(unsigned j = 0;j < _DOF1D;j++){
				tensorM3[alpha][j] =0;
				for(unsigned beta = 0; beta < _QUADRATURE_NODES; beta++){
					tensorM3[alpha][j] += data->VandermondeM[beta][j] * tensorM2[alpha][beta]
											* data->q_weights[alpha]*data->q_weights[beta]
											* data->det_DJe[alpha][beta][e];
				}
			}
		}

		for(unsigned i = 0;i < _DOF1D;i++){
			for(unsigned j = 0;j < _DOF1D;j++){
				for(unsigned alpha = 0; alpha < _QUADRATURE_NODES; alpha++){
					outputvector[data->FEtoDOF[e][i*_DOF1D+j]] += data->VandermondeM[alpha][i] * tensorM3[alpha][j];
				}
			}
		}

	}//element loop
	VecRestoreArrayRead(x,&inputvector);
	VecRestoreArray(y,&outputvector);
	return 0;
}

PetscErrorCode stiffness_mult(Mat A,Vec x,Vec y){
	double x_local[_DOF1D][_DOF1D];
	double tensor_vandermonde_value;

	double tensorA1[_DIM][_QUADRATURE_NODES][_DOF1D];
	double tensorA2[_DIM][_QUADRATURE_NODES][_QUADRATURE_NODES];
	double tensorA3[_DIM][_QUADRATURE_NODES][_QUADRATURE_NODES];
	double tensorA4[_DIM][_DOF1D][_QUADRATURE_NODES];
	double tensorA5[_DOF1D][_DOF1D][_QUADRATURE_NODES];

	GridData *data;
	const double *inputvector;
	double *outputvector;

	MatShellGetContext(A,&data);
	VecGetArrayRead(x,&inputvector);
	VecGetArray(y,&outputvector);

	for(unsigned dof=0;dof<data->global_dof;dof++){
		outputvector[dof]=0;
	}

	unsigned s_idx,i_idx,j_idx;//e=(i,j,s), where s determines which coarse element, i is row of refined

	for(unsigned e=0;e<data->E;e++){
		indices(e,&i_idx,&j_idx,&s_idx,data);

		for(unsigned k=0;k<_DOF1D;k++){
			for(unsigned l=0;l<_DOF1D;l++){
				x_local[k][l]=inputvector[data->FEtoDOF[e][k*_DOF1D+l]];
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
		printf("e: %i\n",e);
		print_3tensor2(tensorA1);


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
//		print_3tensor2(tensorA2);

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

//		print_3tensor2(tensorA3);


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

//		print_3tensor2(tensorA4);


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

//		print_3tensor2(tensorA5);


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
	const double *inputvector;
	double *outputvector;
	double *temparray;
	double ierr=0;

	GridData *data;
	MatShellGetContext(A,&data);

	//PREprocessing
	Vec temp;
	VecDuplicate(x,&temp);
	VecCopy(x,temp);
	VecGetArray(temp,&temparray);

	for(unsigned i=0;i<data->boundary_dof;i++){
		temparray[data->boundary_nodes[i]]=0;
	}

	VecRestoreArray(temp,&temparray);

	//use complete stiffnessmatrix
	MatMult(data->StiffnessM,temp,y);

	//POSTprocessing
	VecGetArrayRead(x,&inputvector);
	VecGetArray(y,&outputvector);

	for(unsigned i=0;i<data->boundary_dof;i++){
		outputvector[data->boundary_nodes[i]]=inputvector[data->boundary_nodes[i]];
	}

	VecRestoreArrayRead(x,&inputvector);
	VecRestoreArray(y,&outputvector);
	return ierr;
}
