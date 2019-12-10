static char help[]="Unittesting\n";
#include "mat2x2.h"
#include "griddata.h"
#include "transformation.h"

#include <stdio.h>
#include <petscmat.h>

int main(int argc,char **args){
	double **M;
	double **M_inv;
	double **M_inv_t;
	double **M_inv_x_M_inv_t;
	double **G;

	alloc_Mat2x2(&M);
	alloc_Mat2x2(&M_inv);
	alloc_Mat2x2(&M_inv_t);
	alloc_Mat2x2(&M_inv_x_M_inv_t);
	alloc_Mat2x2(&G);

	fill_Mat2x2(M,1,2,3,4);
	invert_Mat2x2(M,M_inv);
	printf("M_inv:\n");
	print_Mat2x2(M_inv);
	transpose_Mat2x2(M_inv,M_inv_t);
	printf("M_inv_t:\n");
	print_Mat2x2(M_inv_t);
	MatMult_Mat2x2(M_inv,M_inv_t,M_inv_x_M_inv_t);
	printf("M_inv_x_M_inv_t:\n");
	print_Mat2x2(M_inv_x_M_inv_t);

	fill_Mat2x2(G,5,-3.5,-3.5,2.5);
	differ_Mat2x2(G,M_inv_x_M_inv_t);


	free_Mat2x2(M);
	free_Mat2x2(M_inv);
	free_Mat2x2(M_inv_t);
	free_Mat2x2(M_inv_x_M_inv_t);
	free_Mat2x2(G);

	PetscInt			elements1D=1; 			/*Number of elements in one direction*/
	PetscScalar			L=0.8;

	PetscErrorCode		ierr;
	PetscMPIInt			size;

	GridData			data;

	/*Initalize Petsc and check uniprocessor*/

	ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;

	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
	if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");

	/* Initialise Applicationdata */

	ierr = init_GridData(elements1D,L,&data);CHKERRQ(ierr);

	for(unsigned alpha = 0;alpha < _QUADRATURE_NODES;alpha++){
		for(unsigned beta = 0;beta < _QUADRATURE_NODES;beta++){
			for(unsigned e=0;e<data.E;e++){
				printf("e: %i, a:%i, b:%i, detMat:%f, detAna:%f, diff: %f\n",e,alpha,beta,determinant_Mat2x2(data.DJ[alpha][beta][e]),
						determinant_jacobian_transformation(e,data.q_nodes[alpha],data.q_nodes[beta],&data),
						determinant_Mat2x2(data.DJ[alpha][beta][e])-determinant_jacobian_transformation(e,data.q_nodes[alpha],data.q_nodes[beta],&data));
			}
		}
	}
	for(unsigned i=0;i<data.global_dof;i++){
		Vec x,y;
		ierr = VecCreateSeq(PETSC_COMM_WORLD,data.global_dof,&x);CHKERRQ(ierr);
		VecDuplicate(x,&y);
		VecSetValue(x,i,1,INSERT_VALUES);
		VecAssemblyBegin(x);
		VecAssemblyEnd(x);
		MatMult(data.StiffnessM,x,y);
		VecView(y,PETSC_VIEWER_STDOUT_WORLD);
		VecDestroy(&x);
		VecDestroy(&y);
	}

	free_GridData(&data);
	ierr = PetscFinalize();
	return ierr;
}
