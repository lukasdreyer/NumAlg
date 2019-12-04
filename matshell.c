#include <petscmat.h>
#include "griddata.h"
PetscErrorCode mass_mult(Mat A,Vec x,Vec y){
	GridData *data;
	unsigned *local_dof_indices;
	MatShellGetContext(A,&data);
	for(unsigned e=0;e<data->E;e++){
		local_dof_indices = data->FEtoDOF[e];
	}

	//TODO Test
	return local_dof_indices[0];
}

PetscErrorCode stiffness_mult(Mat A,Vec x,Vec y){
	return 0;
}

PetscErrorCode boundary_mult(Mat A,Vec x, Vec y){
	return 0;
}
