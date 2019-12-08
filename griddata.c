
#include <math.h>
#include "griddata.h"
#include <petscmat.h>
#include "matshell.h"
#include "utilities.h"
#include "transformation.h"
#include "mat2x2.h"

int init_GridData(unsigned M,double L,GridData *data){
	PetscErrorCode ierr=0;
	data->M=M;
	data->E=5*M*M;
	data->L=L;
	data->global_dof=data->E+2*M+1; //5M^2+2M+1
	data->boundary_dof=4*M;

	printf("fill 1dim-quadrature information \n");
	//1 dimensional data
	init_quadrature(data);
	init_VandermondeM(data);
	init_derivativeVandermondeM(data);

	printf("allocate and fill boundary and connectivity information 1 \n");
	//allocate memory and fill matrices/tensors
	init_boundary_nodes(data);
	init_FEtoDOF(data);

	printf("init 2 \n");
	init_D(data);

	printf("init 3 \n");
	init_detDJe(data);

	printf("init 4 \n");
	init_Gepq(data);

	printf("init 5 \n");
	init_GridData_Mat_Vec(data);

	printf("initialisation ready \n");
	return ierr;
}
int free_GridData(GridData *data){
	PetscErrorCode ierr=0;
	free_boundary_nodes(data);
	free_FEtoDOF(data);
	free_GridData_Mat_Vec(data);
	free_D(data);
	free_detDJe(data);
	free_Gepq(data);
	return ierr;
}
int free_boundary_nodes(GridData* data){
	free(data->boundary_nodes);
	return 0;
}
int free_FEtoDOF(GridData *data){
	for(int e=0; e < data->E;e++){
		free(data->FEtoDOF[e]);
	}
	free(data->FEtoDOF);
	return 0;
}
int free_GridData_Mat_Vec(GridData *data){
	PetscErrorCode ierr;
	ierr = VecDestroy(&data->f);CHKERRQ(ierr);
	ierr = VecDestroy(&data->b);CHKERRQ(ierr);
	ierr = MatDestroy(&data->StiffnessM);CHKERRQ(ierr);
	ierr = MatDestroy(&data->MassM);CHKERRQ(ierr);
	ierr = MatDestroy(&data->boundaryStiffnessM);CHKERRQ(ierr);
	return ierr;
}
int free_D(GridData *data){
	for(unsigned alpha=0;alpha<_QUADRATURE_NODES;alpha++){
		for(unsigned beta=0;beta<_QUADRATURE_NODES;beta++){
			for(unsigned e=0;e<data->E;e++){
				free_Mat2x2(data->DJ[alpha][beta][e]);
			}
			free(data->DJ[alpha][beta]);
		}
	}
	return 0;
}
int free_detDJe(GridData *data){
	for(unsigned alpha=0;alpha<_QUADRATURE_NODES;alpha++){
		for(unsigned beta=0;beta<_QUADRATURE_NODES;beta++){
			free(data->det_DJe[alpha][beta]);
		}
	}
	return 0;
}
int free_Gepq(GridData *data){
	for(unsigned alpha=0;alpha<_QUADRATURE_NODES;alpha++){
		for(unsigned beta=0;beta<_QUADRATURE_NODES;beta++){
			for(unsigned e=0;e<data->E;e++){
				free_Mat2x2(data->Gepq[alpha][beta][e]);
			}
			free(data->Gepq[alpha][beta]);
		}
	}
	return 0;
}

int init_quadrature(GridData* data){
	if(_DOF2D!=4) return -1;
	data->q_nodes[0]=-sqrt(1./3);
	data->q_nodes[1]=sqrt(1./3);
	data->q_weights[0]=1;
	data->q_weights[1]=1;
	return 0;
}
int init_VandermondeM(GridData *data){
	data->VandermondeM[0][0]=(1+sqrt(1./3))/2;
	data->VandermondeM[0][1]=(1-sqrt(1./3))/2;
	data->VandermondeM[1][0]=(1-sqrt(1./3))/2;
	data->VandermondeM[1][1]=(1+sqrt(1./3))/2;
	return 0;
}
//TODO:Check
int init_derivativeVandermondeM(GridData *data){
	data->derivativeVandermondeM[0][0]=-.5;
	data->derivativeVandermondeM[0][1]=.5;
	data->derivativeVandermondeM[1][0]=-.5;
	data->derivativeVandermondeM[1][1]=.5;
	return 0;
}

int init_boundary_nodes(GridData *data){
	if(_DOF2D!=4) return -1;
	data->boundary_nodes = malloc(sizeof(int)*data->boundary_dof);
	data->boundary_values = malloc(sizeof(double)*data->boundary_dof);
	for (unsigned i=0;i<data->boundary_dof;i++){
		data->boundary_nodes[i]=data->M*(data->M+3+i);
		//TODO:Input boundary function
		data->boundary_values[i]=0;
	}
	return 0;
}
int init_FEtoDOF(GridData* data){
	if(_DOF2D!=4) return -1;

	unsigned M = data->M;

	data->FEtoDOF = malloc(sizeof(unsigned *) * (data->E));
	for(unsigned e=0;e<data->E;e++){
		data->FEtoDOF[e] = malloc(sizeof(unsigned) * _DOF2D);
	}

	for(unsigned i=0;i<M;i++){
		for(unsigned j=0;j<M;j++){
			//fill center square
			for(unsigned k=0;k<_DOF1D;k++){
				for(unsigned l=0;l<_DOF1D;l++){
					data->FEtoDOF[i*M+j][k*_DOF1D + l]=(i+k)*(M+1)+j +l;
				}
			}
			for(int s=1;s<=4;s++){//(M+1)^2-1
				/* 4 different quartercircle cuts. 
				fill boundary squares, starting right clockwise
				 * */
				for(unsigned k=0;k<_DOF1D;k++){
					for(unsigned l=0;l<_DOF1D;l++){
						data->FEtoDOF[s*M*M+ i*M + j][k*_DOF1D + l]= ((M+2)*M + (s-1)*M*M) +  (i+k)*M +j +l;
					}
				}
			}
			if(i==(M-1)){
				/* adjust boundary between F_1 and F_4 **/
				data->FEtoDOF[i*M+j+4*M*M][2]=(M+2)*M + j;
				data->FEtoDOF[i*M+j+4*M*M][3]=(M+2)*M + j +1;

			}
		}
		/* adjust boundary to center square	 */
		data->FEtoDOF[i*M+M*M][0]=M+i*(M+1);
		data->FEtoDOF[i*M+M*M][2]=M+(i+1)*(M+1);
		data->FEtoDOF[i*M+2*M*M][0]=M*(M+2)-i;
		data->FEtoDOF[i*M+2*M*M][2]=M*(M+2)-(i+1);
		data->FEtoDOF[i*M+3*M*M][0]=M*(M+1)- i*(M+1);
		data->FEtoDOF[i*M+3*M*M][2]=M*(M+1)- (i+1) * (M+1);
		data->FEtoDOF[i*M+4*M*M][0]=i;
		data->FEtoDOF[i*M+4*M*M][2]=i+1;
	}
	return 0;
}
int init_GridData_Mat_Vec(GridData *data){
	PetscErrorCode ierr=0;

	ierr = VecCreateSeq(PETSC_COMM_WORLD,data->global_dof,&(data->f));CHKERRQ(ierr);
	ierr = VecDuplicate(data->f,&(data->b));
	ierr = VecSet(data->f,1);CHKERRQ(ierr);

	ierr = MatCreateShell(PETSC_COMM_WORLD,data->global_dof,data->global_dof,
			PETSC_DECIDE,PETSC_DECIDE,data,&data->MassM);CHKERRQ(ierr);
	ierr = MatCreateShell(PETSC_COMM_WORLD,data->global_dof,data->global_dof,
			PETSC_DECIDE,PETSC_DECIDE,data,&data->StiffnessM);CHKERRQ(ierr);
	ierr = MatCreateShell(PETSC_COMM_WORLD,data->global_dof,data->global_dof,
			PETSC_DECIDE,PETSC_DECIDE,data,&data->boundaryStiffnessM);CHKERRQ(ierr);

	ierr = MatShellSetOperation(data->MassM,MATOP_MULT,(void(*)(void))mass_mult);CHKERRQ(ierr);
	ierr = MatShellSetOperation(data->StiffnessM,MATOP_MULT,(void(*)(void))stiffness_mult);CHKERRQ(ierr);
	ierr = MatShellSetOperation(data->boundaryStiffnessM,MATOP_MULT,(void(*)(void))boundary_mult);CHKERRQ(ierr);
	return ierr;
}
int init_D(GridData *data){
	double c,r,theta;

	unsigned M_sq=data->M*data->M;
	unsigned ij_idx;
	double **DP;
	double **DQ1;
	double **DC;
	double **DR;
	double **DA;
	double **DADP;
	double **DRDADP;

	double x,y;

	alloc_Mat2x2(&DP);
	alloc_Mat2x2(&DADP);
	alloc_Mat2x2(&DRDADP);
	alloc_Mat2x2(&DC);
	alloc_Mat2x2(&DR);
	alloc_Mat2x2(&DA);
	alloc_Mat2x2(&DQ1);

	fill_Mat2x2(DP,1./data->M,0,0,1./data->M);
	fill_Mat2x2(DQ1,0,1,-1,1);

	for(unsigned alpha=0;alpha<_QUADRATURE_NODES;alpha++){
		for(unsigned beta=0;beta<_QUADRATURE_NODES;beta++){
			data->DJ[alpha][beta]=malloc(sizeof(double**)*data->E);
			for(unsigned e=0;e<data->E;e++){
				alloc_Mat2x2(&data->DJ[alpha][beta][e]);
			}
			for(unsigned e=0;e<data->E/5;e++){//all elements in the center square
				fill_Mat2x2(data->DJ[alpha][beta][e], data->L/(2*data->M), 0 , 0, data->L/(2*data->M));
			}
			//TODO:CHECK
			for(unsigned i=0;i<data->M;i++){
				for(unsigned j=0;j<data->M;j++){
					x=data->q_nodes[alpha];
					y=data->q_nodes[beta];

					project(&x,&y,i,j,data);

					coefficients(x,y,&theta,&c,&r,data);
					fill_Mat2x2(DA,1, 0 , 0, M_PI/4);
					fill_Mat2x2(DC,cos(theta), - r*sin(theta) , sin(theta), r*cos(theta));
					fill_Mat2x2(DR, (1-c)/2 ,(1-x)*data->L*sin(theta)/(4*cos(theta)*cos(theta)) , 0, 1);

					//DJ1=DC DR DA DP
					MatMult_Mat2x2(DA,DP,DADP);
					MatMult_Mat2x2(DR,DADP,DRDADP);
					ij_idx=data->M*i+j;
					MatMult_Mat2x2(DC,DRDADP,data->DJ[alpha][beta][M_sq+ij_idx]);
					MatMult_Mat2x2(DQ1,data->DJ[alpha][beta][M_sq+ij_idx],data->DJ[alpha][beta][2*M_sq+ij_idx]);
					MatMult_Mat2x2(DQ1,data->DJ[alpha][beta][2*M_sq+ij_idx],data->DJ[alpha][beta][3*M_sq+ij_idx]);
					MatMult_Mat2x2(DQ1,data->DJ[alpha][beta][3*M_sq+ij_idx],data->DJ[alpha][beta][4*M_sq+ij_idx]);
				}
			}
		}
	}
	free_Mat2x2(DP);
	free_Mat2x2(DADP);
	free_Mat2x2(DRDADP);
	free_Mat2x2(DA);
	free_Mat2x2(DC);
	free_Mat2x2(DQ1);
	return 0;
}
int init_detDJe(GridData *data){
	for(unsigned alpha=0;alpha<_QUADRATURE_NODES;alpha++){
		for(unsigned beta=0;beta<_QUADRATURE_NODES;beta++){

			data->det_DJe[alpha][beta] = malloc(sizeof(double) * data->E);

			//TODO: Check if det F == determinant_jacobian_transformation
			for(unsigned e=0;e<data->E;e++){
				data->det_DJe[alpha][beta][e] = determinant_Mat2x2(data->DJ[alpha][beta][e]);
			}
		}
	}
	return 0;
}
int init_Gepq(GridData *data){
	for(unsigned alpha=0;alpha<_QUADRATURE_NODES;alpha++){
		for(unsigned beta=0;beta<_QUADRATURE_NODES;beta++){
			data->Gepq[alpha][beta]=malloc(sizeof(double**)*data->E);
			assert(data->Gepq[alpha][beta]!=NULL);
			for(unsigned e=0;e<data->E;e++){
				alloc_Mat2x2(&data->Gepq[alpha][beta][e]);
			}
		}
	}

	double** DJ_inv;
	double** DJ_inv_t;

	alloc_Mat2x2(&DJ_inv);
	alloc_Mat2x2(&DJ_inv_t);

	for(unsigned alpha=0;alpha<_QUADRATURE_NODES;alpha++){
		for(unsigned beta=0;beta<_QUADRATURE_NODES;beta++){
			for(unsigned e=1;e<data->E;e++){
				//TODO:CHECK

				invert_Mat2x2(data->DJ[alpha][beta][e],DJ_inv);

				transpose_Mat2x2(DJ_inv,DJ_inv_t);

				MatMult_Mat2x2(DJ_inv,DJ_inv_t,data->Gepq[alpha][beta][e]);//E x Q_N x Q_N x DIM x DIM
			}
		}
	}
	free_Mat2x2(DJ_inv);
	free_Mat2x2(DJ_inv_t);
	return 0;
}
