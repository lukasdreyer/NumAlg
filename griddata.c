
#include <math.h>
#include "griddata.h"
#include <petscmat.h>
#include "matshell.h"
#include "utilities.h"
#include "transformation.h"

int init_GridData(unsigned M,double L,GridData *data){
	PetscErrorCode ierr=0;
	data->M=M;
	data->E=5*M*M;
	data->L=L;
	data->global_dof=data->E+2*M+1; //5M^2+2M+1
	data->global_dof=4*M;


	//1 dimensional data
	init_quadrature(data);
	init_VandermondeM(data);
	init_derivativeVandermondeM(data);

	//allocate memory and fill matrices/tensors
	init_boundary_nodes(data);
	init_FEtoDOF(data);
	init_GridData_Mat_Vec(data);
	init_D(data);
	init_detDJe(data);
	init_Gepq(data);
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
			for(unsigned s=0;s<5;s++){
				for(unsigned i=0;i<_DIM;i++){
					free(data->DF[s][alpha][beta][i]);
				}
				free(data->DF[s][alpha][beta]);
			}
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
				for(unsigned p=0;p<_DIM;p++){
					free(data->Gepq[alpha][beta][e][p]);
				}
				free(data->Gepq[alpha][beta][e]);
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
	data->boundary_nodes = malloc(sizeof(unsigned)*data->boundary_dof);
	for (unsigned i=0;i<data->boundary_dof;i++){
		data->boundary_nodes[i]=data->M*(data->M+3+i);
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

	ierr = VecCreate(PETSC_COMM_WORLD,&(data->f));CHKERRQ(ierr);
	ierr = VecSetSizes(data->f,PETSC_DECIDE,data->global_dof);CHKERRQ(ierr);
	ierr = VecSetFromOptions(data->f);
	ierr = VecDuplicate(data->f,&(data->b));

	ierr = VecSet(data->f,1);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(data->f);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(data->f);CHKERRQ(ierr);

	ierr = VecSet(data->b,0);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(data->b);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(data->b);CHKERRQ(ierr);

	ierr = MatCreateShell(PETSC_COMM_WORLD,data->global_dof,data->global_dof,
			PETSC_DECIDE,PETSC_DECIDE,data,&data->MassM);CHKERRQ(ierr);
	ierr = MatCreateShell(PETSC_COMM_WORLD,data->global_dof,data->global_dof,
			PETSC_DECIDE,PETSC_DECIDE,data,&data->StiffnessM);CHKERRQ(ierr);
	ierr = MatCreateShell(PETSC_COMM_WORLD,data->global_dof,data->global_dof,
			PETSC_DECIDE,PETSC_DECIDE,data,&data->boundaryStiffnessM);CHKERRQ(ierr);

	ierr = MatShellSetOperation(data->MassM,MATOP_MULT,(void(*)(void))mass_mult);CHKERRQ(ierr);
	ierr = MatShellSetOperation(data->StiffnessM,MATOP_MULT,(void(*)(void))stiffness_mult);CHKERRQ(ierr);
	ierr = MatShellSetOperation(data->boundaryStiffnessM,MATOP_MULT,(void(*)(void))boundary_mult);CHKERRQ(ierr);

//	ierr = MatMult(data->MassM,data->f,data->b);CHKERRQ(ierr);
//	ierr = MatMult(data->StiffnessM,data->f,data->b);CHKERRQ(ierr);
	VecView(data->b,PETSC_VIEWER_STDOUT_WORLD);

	return ierr;
}
int init_D(GridData *data){
	double x00,x01,x10,x11;
	double c,r,theta,lambda;

	data->DQ=malloc(sizeof(double*)*_DIM);
	data->DP=malloc(sizeof(double*)*_DIM);
	for(unsigned i=0;i<_DIM;i++){
		data->DQ[i]=malloc(sizeof(double)*_DIM);
		data->DP[i]=malloc(sizeof(double)*_DIM);
	}

	for(unsigned alpha=0;alpha<_QUADRATURE_NODES;alpha++){
		for(unsigned beta=0;beta<_QUADRATURE_NODES;beta++){
			for(unsigned s=0;s<5;s++){
				data->DF[s][alpha][beta] = malloc(sizeof(double*)*_DIM);
				for(unsigned i=0;i<_DIM;i++){
					data->DF[s][alpha][beta][i]=malloc(sizeof(double)*_DIM);
				}
			}

			x00=data->L/2;
			x10=0;
			x01=0;
			x11=data->L/2;

			fill_Mat2x2(data->DF[0][alpha][beta], x00, x01, x10, x11);

			coefficients(data->q_nodes[alpha],data->q_nodes[beta],&theta,&c,&r,data);
			lambda = data->L*M_PI*(1-data->q_nodes[alpha])*sin(theta)/(16*cos(theta)*cos(theta));

			//TODO check F1 values
			x00=(1-c)/2*cos(theta);
			x01=lambda * cos(theta)- M_PI/4*r*sin(theta);
			x10=(1-c)/2*sin(theta);
			x11=lambda * sin(theta)- M_PI/4*r*cos(theta);

			fill_Mat2x2(data->DF[1][alpha][beta], x00, x01, x10, x11);

			MatMatMult2x2(data->DQ,data->DF[1][alpha][beta],data->DF[2][alpha][beta]);
			MatMatMult2x2(data->DQ,data->DF[2][alpha][beta],data->DF[3][alpha][beta]);
			MatMatMult2x2(data->DQ,data->DF[3][alpha][beta],data->DF[4][alpha][beta]);
		}
	}
	return 0;
}
int init_detDJe(GridData *data){
	for(unsigned alpha=0;alpha<_QUADRATURE_NODES;alpha++){
		for(unsigned beta=0;beta<_QUADRATURE_NODES;beta++){
			data->det_DJe[alpha][beta] = malloc(sizeof(double) * data->E);
			for(unsigned e=0;e<data->E;e++){
				data->det_DJe[alpha][beta][e] = determinant_jacobian_transformation
						(e,data->q_nodes[alpha],data->q_nodes[beta],data);
			}
		}
	}
	return 0;
}
int init_Gepq(GridData *data){
	for(unsigned alpha=0;alpha<_QUADRATURE_NODES;alpha++){
		for(unsigned beta=0;beta<_QUADRATURE_NODES;beta++){
			data->Gepq[alpha][beta]=malloc(sizeof(double**)*data->E);
			for(unsigned e=0;e<data->E;e++){
				data->Gepq[alpha][beta][e]=malloc(sizeof(double*)*_DIM);
				for(unsigned p=0;p<_DIM;p++){
					data->Gepq[alpha][beta][e][p]=malloc(sizeof(double)*_DIM);
				}
			}
		}
	}

	double** DFDP;
	double** DJ;
	double** DJ_inv;
	double** DJ_inv_t;

	DFDP=malloc(sizeof(double*)*_DIM);
	DJ=malloc(sizeof(double*)*_DIM);
	DJ_inv=malloc(sizeof(double*)*_DIM);
	DJ_inv_t=malloc(sizeof(double*)*_DIM);
	for(unsigned i=0;i<_DIM;i++){
		DFDP[i]=malloc(sizeof(double)*_DIM);
		DJ[i]=malloc(sizeof(double)*_DIM);
		DJ_inv[i]=malloc(sizeof(double)*_DIM);
		DJ_inv_t[i]=malloc(sizeof(double)*_DIM);
	}
	unsigned s_idx,i_idx,j_idx;

	for(unsigned alpha=0;alpha<_QUADRATURE_NODES;alpha++){
		for(unsigned beta=0;beta<_QUADRATURE_NODES;beta++){
			for(unsigned e=1;e<data->E;e++){
				indices(e,&s_idx,&i_idx,&j_idx,data);
				MatMatMult2x2(data->DF[s_idx][alpha][beta],data->DP,DFDP);
				MatMatMult2x2(data->DQ,DFDP,DJ);
				invertMat2x2(DJ,DJ_inv);
				transposeMat2x2(DJ_inv,DJ_inv_t);
				MatMatMult2x2(DJ_inv,DJ_inv_t,data->Gepq[alpha][beta][e]);//E x Q_N x Q_N x DIM x DIM
			}
		}
	}
	return 0;
}
