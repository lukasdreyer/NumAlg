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
	Vec				x, b;			/* approx solution, RHS, function values */
	Mat			MassM;
//	StiffnessM, MassM, boundaryStiffnessM;			/* stifness and mass matrix */
//	KSP				ksp;				/* linear solver context */
//	PC				pc;					/* preconditioner context */
	PetscInt			M=32;			/*Number of elements in one direction*/
	PetscScalar			L=0.8,norm;

	PetscErrorCode			ierr;
//	PetscInt			its;
	PetscMPIInt			size;

	GridData			data;

	/*Initalize Petsc and check uniprocessor*/

	ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;

	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
	if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");

	/* Initialise Applicationdata */

	ierr = init_GridData(M,L,&data);CHKERRQ(ierr);

	/* Create Vectors x,b,f */

	ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
	ierr = VecSetSizes(x,PETSC_DECIDE,data.global_dof);CHKERRQ(ierr);
	ierr = VecSetFromOptions(x);

	ierr = VecDuplicate(x,&b);CHKERRQ(ierr);

	ierr = VecSet(x,1);CHKERRQ(ierr);
	ierr = VecSet(b,0);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(x);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(x);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(b);CHKERRQ(ierr);


//	ierr = VecDuplicate(x,&f);CHKERRQ(ierr);

	/* Create Matrices A and M */

	ierr = MatCreateShell(PETSC_COMM_WORLD,data.global_dof,data.global_dof,PETSC_DECIDE,PETSC_DECIDE,&data,&MassM);CHKERRQ(ierr);
//	ierr = MatCreateShell(PETSC_COMM_WORLD,data.global_dof,data.global_dof,PETSC_DECIDE,PETSC_DECIDE,&data,&StiffnessM);CHKERRQ(ierr);
//	ierr = MatCreateShell(PETSC_COMM_WORLD,data.global_dof,data.global_dof,PETSC_DECIDE,PETSC_DECIDE,&data,&boundaryStiffnessM);CHKERRQ(ierr);

	/*
	* Fill Matrices A and M
	*/

	ierr = MatShellSetOperation(MassM,MATOP_MULT,(void(*)(void))mass_mult);CHKERRQ(ierr);
//	ierr = MatShellSetOperation(StiffnessM,MATOP_MULT,(void(*)(void))stiffness_mult);CHKERRQ(ierr);
//	ierr = MatShellSetOperation(boundaryStiffnessM,MATOP_MULT,(void(*)(void))boundary_mult);CHKERRQ(ierr);

	/* Fill Vec f and b=M*f */
//TODO
//	ierr = assembleVecFct(f,&data);CHKERRQ(ierr);
	ierr = MatMult(MassM,x,b);CHKERRQ(ierr);
	ierr = VecDot(x,b,&norm);CHKERRQ(ierr);
	printf("error: %.10f,norm: %.10f,sqrt(pi): %.10f, \n", fabs(sqrt(norm)-sqrt(M_PI)),sqrt(norm),sqrt(M_PI));


	/* Create a Linear solver (Krylov space) */
/*
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,boundaryStiffnessM,boundaryStiffnessM);CHKERRQ(ierr);


	ierr = KSPSetType(ksp,KSPCG);CHKERRQ(ierr);

	ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
	ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);//PCAMG

	ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);


	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);


	ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);


	ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);


	ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n%D Iterations \n",its);CHKERRQ(ierr);

	ierr = display(x,b,f,A,M,&appctx);CHKERRQ(ierr);


	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

*/
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	ierr = VecDestroy(&b);CHKERRQ(ierr);
//	ierr = VecDestroy(&f);CHKERRQ(ierr);
//	ierr = MatDestroy(&StiffnessM);CHKERRQ(ierr);
	ierr = MatDestroy(&MassM);CHKERRQ(ierr);
//	ierr = MatDestroy(&boundaryStiffnessM);CHKERRQ(ierr);

	free_GridData(&data);
	/* Always call PetscFinalize() before exiting a program. */
	ierr = PetscFinalize();
	return ierr;


}
