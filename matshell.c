#include <petscmat.h>
#include<petscsys.h>
#include "griddata.h"
#include "transformation.h"

	/*
	for(unsigned  =0; <  ;  ++){
		for(unsigned  =0; <  ;  ++){
			for(unsigned  =0; <  ;  ++){
					for(unsigned  =0; <  ;  ++){
				}
			}
		}
	}
	*/

PetscErrorCode mass_mult(Mat A,Vec x,Vec y){
	GridData *data;
	unsigned *local_dof_indices;
	const double *inputvector;
	double *outputvector;
	double x_local[_DOF1D][_DOF1D];

	double tensorM1[_QUADRATURE_NODES][_DOF1D];
	double tensorM2[_QUADRATURE_NODES][_QUADRATURE_NODES];
	double tensorM3[_DOF1D][_QUADRATURE_NODES];

	VecGetArrayRead(x,&inputvector);
	VecGetArray(y,&outputvector);
	MatShellGetContext(A,&data);

	for(unsigned e=0;e<data->E;e++){
		local_dof_indices = data->FEtoDOF[e];
		for(unsigned i=0;i<_DOF1D;i++){
			for(unsigned j=0;j<_DOF1D;j++){
				x_local[i][j]=inputvector[local_dof_indices[i*_DOF1D+j]];//symmetry
			}
		}


		for(unsigned alpha = 0;alpha < _QUADRATURE_NODES;alpha++){
			for(unsigned j = 0;j < _DOF1D;j++){
				tensorM1[alpha][j] =0;
				for(unsigned i = 0; i < _DOF1D; i++){
					tensorM1[alpha][j] += data->VandermondeM[i][alpha]*x_local[i][j];
				}
			}
		}
		for(unsigned alpha = 0;alpha < _QUADRATURE_NODES;alpha++){
			for(unsigned beta = 0;beta < _QUADRATURE_NODES;beta++){
				tensorM2[alpha][beta] =0;
				for(unsigned j = 0; j < _DOF1D; j++){
					tensorM2[alpha][beta] += data->VandermondeM[j][beta]*tensorM1[alpha][j];
				}
			}
		}

		for(unsigned k = 0;k < _DOF1D;k++){
			for(unsigned beta = 0;beta < _QUADRATURE_NODES;beta++){
				tensorM3[k][beta] =0;
				for(unsigned alpha = 0; alpha < _QUADRATURE_NODES; alpha++){
					tensorM3[k][beta] += data->VandermondeM[k][alpha] * tensorM2[alpha][beta]
											* determinant_jacobian_transformation(e,data->q_nodes[alpha],data->q_nodes[beta],data)
											* data->q_weights[alpha]*data->q_weights[beta];
				}
			}
		}

		for(unsigned k = 0;k < _DOF1D;k++){
			for(unsigned l = 0;l < _DOF1D;l++){
				for(unsigned beta = 0; beta < _QUADRATURE_NODES; beta++){
					outputvector[local_dof_indices[k*_DOF1D+l]] +=data->VandermondeM[l][beta]*tensorM3[k][beta];
				}
			}
		}

	}

	VecRestoreArrayRead(x,&inputvector);
	VecRestoreArray(y,&outputvector);


	return 0;
}

PetscErrorCode stiffness_mult(Mat A,Vec x,Vec y){
	double tensorA1[_DIM][_QUADRATURE_NODES][_DOF1D];
	double tensorA2[_DIM][_QUADRATURE_NODES][_QUADRATURE_NODES];
	double tensorA3[_DIM][_QUADRATURE_NODES][_QUADRATURE_NODES];
	double tensorA4[_DIM][_QUADRATURE_NODES][_DOF1D];
	double tensorA5[_DIM][_DOF1D][_DOF1D];

	GridData *data;
	unsigned *local_dof_indices;
	const double *inputvector;
	double *outputvector;
	double x_local[_DOF1D][_DOF1D];

	VecGetArrayRead(x,&inputvector);
	VecGetArray(y,&outputvector);

	MatShellGetContext(A,&data);

	for(unsigned e=0;e<data->E;e++){
		local_dof_indices = data->FEtoDOF[e];
		for(unsigned i=0;i<_DOF1D;i++){
			for(unsigned j=0;j<_DOF1D;j++){
				x_local[i][j]=inputvector[local_dof_indices[i*_DOF1D+j]];//symmetry
			}
		}
		for(unsigned q =0; q < _DIM ; q++){
			for(unsigned  alpha=0;alpha < _QUADRATURE_NODES ;  alpha++){
				for(unsigned  j=0;j < _DOF1D ;  j++){
					tensorA1[q][alpha][j]=0;
					for(unsigned  i=0; i< _DOF1D ;  i++){
						tensorA1[q][alpha][j]+=TensorVandermonde(q,1,i,alpha,data)*x_local[i][j];
					}
				}
			}
		}

		for(unsigned  q=0;q < _DIM ; q ++){
			for(unsigned  alpha=0; alpha< _QUADRATURE_NODES ; alpha ++){
				for(unsigned beta =0; beta < _QUADRATURE_NODES ;  beta++){
					tensorA2[q][alpha][beta]=0;
					for(unsigned  j=0;j < _DOF1D ; j ++){
						tensorA2[q][alpha][beta] += TensorVandermonde(q,0,j,beta,data)*tensorA1[q][alpha][j];
					}
				}
			}
		}

		for(unsigned  p=0; p<  _DIM;  p++){
			for(unsigned  alpha=0; alpha< _QUADRATURE_NODES ;  alpha++){
				for(unsigned  beta=0; beta< _QUADRATURE_NODES ;  beta++){
					tensorA3[p][alpha][beta]=0;
					for(unsigned  q=0; q< _DOF1D ;  q++){
						tensorA3[p][alpha][beta] += geometryterm(p,q,alpha,beta,data)*tensorA2[q][alpha][beta];
					}
				}
			}
		}

		for(unsigned  p=0; p<  _DIM;  p++){
			for(unsigned  alpha=0; alpha< _QUADRATURE_NODES ;  alpha++){
				for(unsigned  l=0; l< _DOF1D ;  l++){
					tensorA4[p][alpha][l]=0;
					for(unsigned  beta=0; beta< _DOF1D ;  beta++){
						tensorA4[p][alpha][l] += tensorA3[p][alpha][beta] * TensorVandermonde(p,0,l,beta,data)
										* determinant_jacobian_transformation(e,data->q_nodes[alpha],data->q_nodes[beta],data)
													* data->q_weights[alpha]*data->q_weights[beta];
					}
				}
			}
		}

		for(unsigned  p=0; p< _DIM ;  p++){
			for(unsigned  k=0; k< _DOF1D ;  k++){
				for(unsigned  l=0; l< _DOF1D ;  l++){
					tensorA5[p][k][l]=0;
					for(unsigned  alpha=0; alpha< _QUADRATURE_NODES ;  alpha++){
						tensorA5[p][k][l]+=TensorVandermonde(p,1,k,alpha,data)*tensorA4[p][alpha][l];
					}
				}
			}
		}

		for(unsigned k =0; k< _DOF1D ;  k++){
			for(unsigned l =0; l< _DOF1D ;  l++){
				outputvector[local_dof_indices[k*_DOF1D+l]]=0;
				for(unsigned p =0; p< _DIM ; p ++){
					outputvector[local_dof_indices[k*_DOF1D+l]] += tensorA5[p][k][l];
				}
			}
		}

	}//end element loop

	VecRestoreArrayRead(x,&inputvector);
	VecRestoreArray(y,&outputvector);

	return 0;
}

PetscErrorCode boundary_mult(Mat A,Vec x, Vec y){
	return 0;
}
