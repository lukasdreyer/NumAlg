static char help[]="Solves -laplace(u)=f approximately on a circular-grid.\n\nOptions:\n\n"

		"\n\n\n"
		"Important petsc_options:"
		"\n\n"
		"-pc_type <jacobi>: Preconditioner (one of) none jacobi pbjacobi vpbjacobi bjacobi sor lu shell mg eisenstat ilu icc cholesky asm gasm ksp composite redundant nn mat fieldsplit galerkin exotic cp lsc redistribute svd gamg kaczmarz telescope patch tfs bddc lmvm (PCSetType)\n"
		"-ksp_type <richardson>: Krylov method (one of) cg groppcg pipecg pipecgrr pipelcg cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs pipebcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls fetidp (KSPSetType)\n"
		"-ksp_max_it <10000>: Maximum number of iterations (KSPSetTolerances)\n"
		"-ksp_rtol <1e-05>: Relative decrease in residual norm (KSPSetTolerances)\n"
		"-ksp_richardson_scale <1.>: damping factor (KSPRichardsonSetScale)\n"
		"-ksp_richardson_self_scale: <FALSE> dynamically determine optimal damping factor (KSPRichardsonSetSelfScale)\n"

		 "\n\n\n\n" ;


#include <petscksp.h>
#include "griddata.h"
#include "matshell.h"



int main(int argc,char **args){
	Vec				x;					/* approx solution, RHS, function values */

	KSP				ksp;				/* linear solver context */
	PC				pc;					/* preconditioner context */
	//TODO: add als console input
	PetscInt			M=1000; 			/*Number of elements in one direction*/
	PetscScalar			L=0.8,norm=0;

	PetscErrorCode		ierr;
	PetscInt			its;
	PetscMPIInt			size;

	GridData			data;

	/*Initalize Petsc and check uniprocessor*/

	ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;

	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
	if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");

	/* Initialise Applicationdata */


	ierr = init_GridData(M,L,&data);CHKERRQ(ierr);


	ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
	ierr = VecSetSizes(x,PETSC_DECIDE,data.global_dof);CHKERRQ(ierr);
	ierr = VecSetFromOptions(x);


//	VecView(data.b,PETSC_VIEWER_STDOUT_WORLD);
//	VecView(data.f,PETSC_VIEWER_STDOUT_WORLD);

	ierr = VecDot(data.b,data.f,&norm);CHKERRQ(ierr);
	printf("error: %.10f,norm: %.10f,pi: %.10f, \n", fabs(norm-M_PI),norm,M_PI);

	ierr = VecDestroy(&x);CHKERRQ(ierr);
	free_GridData(&data);
	ierr = PetscFinalize();
	return ierr;


	/* Create a Linear solver (Krylov space) */

	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,data.boundaryStiffnessM,data.boundaryStiffnessM);CHKERRQ(ierr);

	ierr = KSPSetType(ksp,KSPCG);CHKERRQ(ierr);

	ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
	ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);//PCAMG

	ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

	ierr = KSPSolve(ksp,data.f,x);CHKERRQ(ierr);

	ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n%D Iterations \n",its);CHKERRQ(ierr);

	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

	/* Always call PetscFinalize() before exiting a program. */

	ierr = VecDestroy(&x);CHKERRQ(ierr);
	free_GridData(&data);
	ierr = PetscFinalize();
	return ierr;
}
