#include <petscmat.h>
#include <petscsys.h>
#include "griddata.h"
#include "transformation.h"
#include "utilities.h"
#include "matshell.h"



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
	VecDestroy(&temp);

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


		//TODO: same as stiffnessmult, macro tensor2mult
		for(unsigned i=0;i<_DOF1D;i++){
			for(unsigned j=0;j<_DOF1D;j++){
				x_local[i][j]=inputvector[data->FEtoDOF[e][i][j]];
			}
		}

		for(unsigned k = 0;k < _DOF1D;k++){
			for(unsigned alpha = 0;alpha < _QUADRATURE_NODES;alpha++){
				tensorM1[k][alpha] =0;
				for(unsigned l = 0; l < _DOF1D; l++){
					tensorM1[k][alpha] += data->VandermondeM[0][alpha][l]*x_local[k][l];
				}
			}
		}

//		TENSOR2MULT(tensorM1, beta, alpha, 1, k, _QUADRATURE_NODES, _QUADRATURE_NODES, _DOF1D,
//				data->VandermondeM[0][beta][k]* data->W[alpha][beta][e],tensorM2)

		for(unsigned beta = 0;beta < _QUADRATURE_NODES;beta++){
			for(unsigned alpha = 0;alpha < _QUADRATURE_NODES;alpha++){
				tensorM2[beta][alpha] =0;
				for(unsigned k = 0; k < _DOF1D; k++){
					tensorM2[beta][alpha] += data->VandermondeM[0][beta][k]* data->W[alpha][beta][e]*
							tensorM1[k][alpha];
				}
			}
		}

//		TENSOR2MULT(tensorM2, beta, j, 2, alpha, _QUADRATURE_NODES, _DOF1D, _QUADRATURE_NODES,
//				data->VandermondeM[0][alpha][j],tensorM3)

		for(unsigned beta = 0;beta < _QUADRATURE_NODES;beta++){
			for(unsigned j = 0;j < _DOF1D;j++){
				tensorM3[beta][j] =0;
				for(unsigned alpha = 0; alpha < _QUADRATURE_NODES; alpha++){
					tensorM3[beta][j] += data->VandermondeM[0][alpha][j] * tensorM2[beta][alpha] ;
				}
			}
		}

		for(unsigned i = 0;i < _DOF1D;i++){
			for(unsigned j = 0;j < _DOF1D;j++){
				for(unsigned beta = 0; beta < _QUADRATURE_NODES; beta++){
					outputvector[data->FEtoDOF[e][i][j]] += data->VandermondeM[0][beta][i] * tensorM3[beta][j];
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

//only center coarse element
//	for(unsigned e=0;e<data->E/5;e++){
//only outer four coarse elements
//	for(unsigned e=data->E/5;e<data->E;e++){

	for(unsigned e=0;e<data->E;e++){
		for(unsigned k=0;k<_DOF1D;k++){
			for(unsigned l=0;l<_DOF1D;l++){
				x_local[k][l]=inputvector[data->FEtoDOF[e][k][l]];
			}
		}

		double tensorA1[_DIM][_QUADRATURE_NODES][_DOF1D];
		double tensorA2[_DIM][_QUADRATURE_NODES][_QUADRATURE_NODES];
		double tensorA3[_DIM][_QUADRATURE_NODES][_QUADRATURE_NODES];
		double tensorA4[_DIM][_DOF1D][_QUADRATURE_NODES];
		double tensorA5[_DOF1D][_DOF1D][_QUADRATURE_NODES];

		for(unsigned p =0; p < _DIM ; p++){
			for(unsigned  beta=0;beta < _QUADRATURE_NODES ;  beta++){
				for(unsigned  l=0;l < _DOF1D ;  l++){
					tensorA1[p][beta][l]=0;
					for(unsigned  k=0; k< _DOF1D ;  k++){
						tensorA1[p][beta][l]+=data->VandermondeM[p==1][beta][k]*x_local[k][l];
					}
				}
			}
		}
		TENSOR3MULT(tensorA1, p, beta, alpha, 3, l, _DIM, _QUADRATURE_NODES, _QUADRATURE_NODES, _DOF1D,
				data->VandermondeM[p==0][alpha][l],tensorA2)

		TENSOR3MULT(tensorA2, q, beta, alpha, 1, p, _DIM, _QUADRATURE_NODES, _QUADRATURE_NODES, _DIM,
				data->Gepq[alpha][beta][e][p][q]*data->W[alpha][beta][e],tensorA3)

		TENSOR3MULT(tensorA3, q, i, alpha, 2, beta, _DIM, _DOF1D, _QUADRATURE_NODES, _QUADRATURE_NODES,
				data->VandermondeM[q==1][beta][i],tensorA4)

		TENSOR3MULT(tensorA4, q, i, j, 3, alpha, _DIM, _DOF1D, _DOF1D, _QUADRATURE_NODES,
				data->VandermondeM[q==0][alpha][j],tensorA5)

		for(unsigned i =0; i< _DOF1D ;  i++){
			for(unsigned j =0; j< _DOF1D ;  j++){
				for(unsigned q=0;q<_DIM;q++){
					outputvector[data->FEtoDOF[e][i][j]] += tensorA5[q][i][j];
				}
			}
		}
	}//end element loop

	VecRestoreArrayRead(x,&inputvector);
	VecRestoreArray(y,&outputvector);

	return 0;
}

PetscErrorCode stiffness_mult_primitive(Mat A,Vec x,Vec y){
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

//only center coarse element
//	for(unsigned e=0;e<data->E/5;e++){
//only outer four coarse elements
//	for(unsigned e=data->E/5;e<data->E;e++){

	for(unsigned e=0;e<data->E;e++){
		for(unsigned k=0;k<_DOF1D;k++){
			for(unsigned l=0;l<_DOF1D;l++){
				x_local[k][l]=inputvector[data->FEtoDOF[e][k][l]];
			}
		}

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
	}//end element loop

	VecRestoreArrayRead(x,&inputvector);
	VecRestoreArray(y,&outputvector);

	return 0;
}
