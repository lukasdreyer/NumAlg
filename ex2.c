static char help[]="Solves -laplace(u)=f approximately on a circular-grid.\n\nOptions:\n\n"
		"-l:   sidelength of square (should be less than sqrt(2)!)\n"
		"-m:   number of squares in 1 direction\n"
		"\n\n";


#include <petscksp.h>
#include "griddata.h"
#include "matshell.h"
#include "testfunctions.h"
#include "utilities.h"
#include "solver.h"



int main(int argc,char **args){
	PetscErrorCode		ierr;		//Petsc error handling

	GridData			*data;//Maybe allocate memory for struct as well
	Vec				x;			//stores solution

	/*Initalize Petsc */

	ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;

	//Check uniprocessor
	PetscMPIInt			size;
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
	if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");

	/* Initialise Applicationdata */


	ierr = init_GridData_Input(&data);CHKERRQ(ierr);

	//Get memory for solution vector
	ierr = VecDuplicate(data->f,&x);

	/*solve*/

	char filename1[] = "out_linear.txt";
	data->tfdata.A=4;
	set_boundary_initial(linearAx,zero_fct,data);
	printf("data1 set\n");

	solve(x,data);
	VecAXPY(x,-1,data->u_ana);//error
	print_gnuplot_vec_circle(filename1,x,data);

	/********************/

	char filename2[] = "out_radiussquared.txt";
	data->tfdata.A=4;
	set_boundary_initial(radiussquared,constantA,data);
	printf("data2 set\n");

	solve(x,data);
	VecAXPY(x,-1,data->u_ana);//error
	print_gnuplot_vec_circle(filename2,x,data);

	/* Always call PetscFinalize() before exiting a program. */

	ierr = VecDestroy(&x);CHKERRQ(ierr);
	free_GridData(data);
	ierr = PetscFinalize();
	return ierr;
}
