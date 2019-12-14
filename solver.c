#include <petscksp.h>
#include "griddata.h"
#include "solver.h"

PetscErrorCode inverse_diag(PC pc,Vec x,Vec y ){
	GridData *data;
	PCShellGetContext(pc,(void**)&data);
	VecPointwiseMult(y,x,data->lumped_mass_diag);
	return 0;
}

int init_rhs(Vec *b,GridData *data){
	PetscErrorCode ierr = 0;
	Vec  g,b0;

	ierr = VecDuplicate(data->f,b);CHKERRQ(ierr);
	ierr = VecDuplicate(data->f,&g);CHKERRQ(ierr);
	ierr = VecDuplicate(data->f,&b0);CHKERRQ(ierr);


	ierr = VecSet(g,0);CHKERRQ(ierr);
	ierr = VecSetValues(g,data->boundary_dof,data->boundary_nodes,data->boundary_values,INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(g);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(g);CHKERRQ(ierr);

	ierr = MatMult(data->StiffnessM,g,*b);CHKERRQ(ierr);
	ierr = MatMult(data->MassM,data->f,b0);CHKERRQ(ierr);

	ierr = VecAYPX(*b,-1,b0);CHKERRQ(ierr);

	//Overwrite values in b on boundary by the given boundary values.
	ierr = VecSetValues(*b,data->boundary_dof,data->boundary_nodes,data->boundary_values,INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(*b);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(*b);CHKERRQ(ierr);

	ierr = VecDestroy(&g);CHKERRQ(ierr);
	ierr = VecDestroy(&b0);CHKERRQ(ierr);
	return 0;
}

int solve(Vec x,GridData *data){
	KSP				ksp;				/* linear solver context */
	PC				pc;					/* preconditioner context */

	PetscErrorCode		ierr;
	PetscInt			its;
	Vec 				b;

	init_rhs(&b,data);//allocates memory for b

	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,data->boundaryStiffnessM,data->boundaryStiffnessM);CHKERRQ(ierr);

	ierr = KSPSetType(ksp,KSPCG);CHKERRQ(ierr);

	ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);

//	ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);


	ierr = PCSetType(pc,PCSHELL);CHKERRQ(ierr);
	ierr = PCShellSetApply(pc,inverse_diag);CHKERRQ(ierr);
	ierr = PCShellSetContext(pc,data);CHKERRQ(ierr);


	ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,1000);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

	ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

	ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
//	ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n%D Iterations \n",its);CHKERRQ(ierr);

	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
	//Destroy_rhs
	ierr = VecDestroy(&b);CHKERRQ(ierr);
	return ierr;
}
