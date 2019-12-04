#include <petscmat.h>
#include "griddata.h"
#include "transformation.h"

PetscErrorCode mass_mult(Mat A,Vec x,Vec y){
	GridData *data;
	unsigned *local_dof_indices;\
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
		for(unsigned alpha = 0;alpha < _QUADRATURE_NODES;alpha++){
			for(unsigned j = 0;j < _DOF1D;j++){
				data->tensorM1[alpha][j] =0;
				for(unsigned i = 0; i < _DOF1D; i++){
					data->tensorM1[alpha][j] += data->VandermondeM[i][alpha]*x_local[i][j];
				}
			}
		}
		for(unsigned alpha = 0;alpha < _QUADRATURE_NODES;alpha++){
			for(unsigned beta = 0;beta < _QUADRATURE_NODES;beta++){
				data->tensorM2[alpha][beta] =0;
				for(unsigned j = 0; j < _DOF1D; j++){
					data->tensorM2[alpha][beta] += data->VandermondeM[j][beta]*data->tensorM1[alpha][j];
				}
			}
		}

		for(unsigned k = 0;k < _DOF1D;k++){
			for(unsigned beta = 0;beta < _QUADRATURE_NODES;beta++){
				data->tensorM3[k][beta] =0;
				for(unsigned alpha = 0; alpha < _QUADRATURE_NODES; alpha++){
					data->tensorM3[k][beta] += data->VandermondeM[k][alpha] * data->tensorM2[alpha][beta]
											* determinant_jacobian_transformation(e,data->q_nodes[alpha],data->q_nodes[beta],data)
											* data->q_weights[alpha]*data->q_weights[beta];
				}
			}
		}

		for(unsigned k = 0;k < _DOF1D;k++){
			for(unsigned l = 0;l < _DOF1D;l++){
				for(unsigned beta = 0; beta < _QUADRATURE_NODES; beta++){
					outputvector[local_dof_indices[k*_DOF1D+l]] +=data->VandermondeM[l][beta]*data->tensorM3[k][beta];
				}
			}
		}

	}

	VecRestoreArrayRead(x,&inputvector);
	VecRestoreArray(y,&outputvector);


	return 0;
}

PetscErrorCode stiffness_mult(Mat A,Vec x,Vec y){
	return 0;
}

PetscErrorCode boundary_mult(Mat A,Vec x, Vec y){
	return 0;
}
