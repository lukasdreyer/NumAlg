PetscErrorCode mass_mult(Mat A,Vec x,Vec y){
	GridData *data;
	double *local_dof_indices;
	MatShellGetContext(A,&data);
	for(unsigned e=0;e<data->E;e++){
		local_dof_indices = data->FEtoDOF[e];
	}

	//TODO Test
	return local_dof_indices;
}

PetscErrorCode stiffness_mult(Mat A,Vec x,Vec y){
	return 0;
}

PetscErrorCode boundary_mult(Mat A,Vec x, Vec y){
	return 0;
}
