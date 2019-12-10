static char help[]="Unittesting\n";
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "utilities.h"
#include "griddata.h"
//#include "transformation.h"


//#include "mat2x2.h"
#include "griddata.h"

#include <petscmat.h>

int main(int argc,char **args){
	PetscInt			M=1; 			/*Number of elements in one direction*/
	PetscScalar			L=0.8,norm=0;

	PetscErrorCode		ierr;
	PetscMPIInt			size;

	GridData			data;

	/*Initalize Petsc and check uniprocessor*/

	ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;

	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
	if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");

	/* Initialise Applicationdata */



	for(unsigned i=0;i<7;i++){
		ierr = init_GridData(M,L,&data);CHKERRQ(ierr);
		printf("MatMult VecDot\n");
		ierr = MatMult(data.MassM,data.f,data.b);CHKERRQ(ierr);

		ierr = VecDot(data.b,data.f,&norm);CHKERRQ(ierr);
		printf("error: %.10f,norm: %.10f,pi: %.10f, \n", fabs(norm-M_PI),norm,M_PI);
		ierr = free_GridData(&data);CHKERRQ(ierr);
		M*=2;
	}
	ierr = PetscFinalize();
	return ierr;
}
