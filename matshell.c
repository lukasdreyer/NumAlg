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
				x_local[i][j]=inputvector[data->FEtoDOF[e][i][j]];//symmetry
			}
		}

		for(unsigned k = 0;k < _DOF1D;k++){
			for(unsigned beta = 0;beta < _QUADRATURE_NODES;beta++){
				tensorM1[k][beta] =0;
				for(unsigned l = 0; l < _DOF1D; l++){
					tensorM1[k][beta] += data->VandermondeM[0][beta][l]*x_local[k][l];
				}
			}
		}

		for(unsigned alpha = 0;alpha < _QUADRATURE_NODES;alpha++){
			for(unsigned beta = 0;beta < _QUADRATURE_NODES;beta++){
				tensorM2[alpha][beta] =0;
				for(unsigned k = 0; k < _DOF1D; k++){
					tensorM2[alpha][beta] += data->VandermondeM[0][alpha][k]*tensorM1[k][beta];
				}
			}
		}

		for(unsigned alpha = 0;alpha < _QUADRATURE_NODES;alpha++){
			for(unsigned j = 0;j < _DOF1D;j++){
				tensorM3[alpha][j] =0;
				for(unsigned beta = 0; beta < _QUADRATURE_NODES; beta++){
					tensorM3[alpha][j] += data->VandermondeM[0][beta][j] * tensorM2[alpha][beta]
											* data->W[alpha][beta][e];
				}
			}
		}

		for(unsigned i = 0;i < _DOF1D;i++){
			for(unsigned j = 0;j < _DOF1D;j++){
				for(unsigned alpha = 0; alpha < _QUADRATURE_NODES; alpha++){
					outputvector[data->FEtoDOF[e][i][j]] += data->VandermondeM[0][alpha][i] * tensorM3[alpha][j];
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
		for(unsigned i =0; i< _DOF1D ;  i++){
			for(unsigned j =0; j< _DOF1D ;  j++){
				outputvector[data->FEtoDOF[e][i][j]] = 0;
			}
		}
	}

	for(unsigned e=0;e<data->E;e++){
//only center coarse element
//	for(unsigned e=0;e<data->E/5;e++){
//only outer four coarse elements
//	for(unsigned e=data->E/5;e<data->E;e++){
		indices(e,&i_idx,&j_idx,&s_idx,data);

		for(unsigned k=0;k<_DOF1D;k++){
			for(unsigned l=0;l<_DOF1D;l++){
				x_local[k][l]=inputvector[data->FEtoDOF[e][k][l]];
			}
		}
#if 1
		for(unsigned i =0; i< _DOF1D ;  i++){
			for(unsigned j =0; j< _DOF1D ;  j++){
				for(unsigned alpha = 0; alpha < _QUADRATURE_NODES ; alpha++){
					for(unsigned beta= 0; beta < _QUADRATURE_NODES ;  beta++){
						for(unsigned p= 0;  p<_DIM  ; p++){
							for(unsigned q= 0;  q< _DIM ;  q++){
								for(unsigned k= 0;  k< _DOF1D ;  k++){
									for(unsigned l= 0; l < _DOF1D ; l++){
										double product=x_local[k][l];
										product*=data->W[alpha][beta][e];
										product*=data->VandermondeM[p==0][alpha][l];
										product*=data->VandermondeM[p==1][beta][k];
										product*=data->Gepq[alpha][beta][e][p][q];
										product*=data->VandermondeM[q==0][alpha][j];
										product*=data->VandermondeM[q==1][beta][i];

										outputvector[data->FEtoDOF[e][i][j]] += product;
									}
								}
							}
						}
					}
				}
			}
		}
#else
		double tensorA1[_DIM][_QUADRATURE_NODES][_DOF1D];
		double tensorA2[_DIM][_QUADRATURE_NODES][_QUADRATURE_NODES];
		double tensorA3[_DIM][_QUADRATURE_NODES][_QUADRATURE_NODES];
		double tensorA4[_DIM][_DOF1D][_QUADRATURE_NODES];
		double tensorA5[_DOF1D][_DOF1D][_QUADRATURE_NODES];

		for(unsigned p =0; p < _DIM ; p++){
			for(unsigned  beta=0;beta < _QUADRATURE_NODES ;  beta++){
				for(unsigned  k=0;k < _DOF1D ;  k++){
					tensorA1[p][beta][k]=0;
					for(unsigned  l=0; l< _DOF1D ;  l++){
						tensorA1[p][beta][k]+=data->VandermondeM[p==1][beta][l]*x_local[k][l];
					}
				}
			}
		}
//		printf("e: %i\n",e);
//		print_3tensor2(tensorA1);

		for(unsigned  p=0;p < _DIM ; p ++){
			for(unsigned beta =0; beta < _QUADRATURE_NODES ;  beta++){
				for(unsigned  alpha=0; alpha< _QUADRATURE_NODES ; alpha ++){
					tensorA2[p][beta][alpha]=0;
					for(unsigned  k=0;k < _DOF1D ; k ++){
						tensorA2[p][beta][alpha] += data->VandermondeM[alpha][k][p==0]*tensorA1[p][beta][k];
					}
				}
			}
		}
//		print_3tensor2(tensorA2);

		for(unsigned  q=0; q<  _DIM;  q++){
			for(unsigned  beta=0; beta< _QUADRATURE_NODES ;  beta++){
				for(unsigned  alpha=0; alpha< _QUADRATURE_NODES ;  alpha++){
					tensorA3[q][beta][alpha]=0;
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
						tensorA4[q][j][alpha] += tensorA3[q][beta][alpha] * data->VandermondeM[q==1][beta][j]
										* data->W[alpha][beta][e];
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
						tensorA5[i][j][alpha]+=data->VandermondeM[q==0][alpha][i]*tensorA4[q][j][alpha];
					}
				}
			}
		}

//		print_3tensor2(tensorA5);


		for(unsigned i =0; i< _DOF1D ;  i++){
			for(unsigned j =0; j< _DOF1D ;  j++){
				for(unsigned  alpha=0; alpha< _QUADRATURE_NODES ;  alpha++){
					outputvector[data->FEtoDOF[e][i][j]] += tensorA5[i][j][alpha];
				}
			}
		}
#endif
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


//	VecView(x,PETSC_VIEWER_STDOUT_WORLD);

	//PREprocessing
	Vec temp;
	VecDuplicate(x,&temp);
	VecCopy(x,temp);
	VecGetArray(temp,&temparray);

	for(unsigned i=0;i<data->boundary_dof;i++){
		temparray[data->boundary_nodes[i]]=0;
	}

	VecRestoreArray(temp,&temparray);
//	VecView(temp,PETSC_VIEWER_STDOUT_WORLD);


	//use complete stiffnessmatrix
	MatMult(data->StiffnessM,temp,y);

	//POSTprocessing
	VecGetArrayRead(x,&inputvector);
	VecGetArray(y,&outputvector);

	for(unsigned i=0;i<data->boundary_dof;i++){
		outputvector[data->boundary_nodes[i]] = inputvector[data->boundary_nodes[i]];
	}

	VecRestoreArrayRead(x,&inputvector);
	VecRestoreArray(y,&outputvector);
//	VecView(y,PETSC_VIEWER_STDOUT_WORLD);
	return ierr;
}
