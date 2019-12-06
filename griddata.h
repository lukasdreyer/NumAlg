#ifndef GRIDDATA_H
#define GRIDDATA_H

#define _QUADRATURE_NODES 2
#define _DOF1D 2
#define _DOF2D _DOF1D*_DOF1D
#define _DIM 2

#include <stdlib.h>
#include <petscmat.h>
#include "matshell.h"
#include "utilities.h"

typedef struct GridData{
	double L;
	unsigned M;
	unsigned E;
	unsigned global_dof,boundary_dof;

	double q_weights[_QUADRATURE_NODES];
	double q_nodes[_QUADRATURE_NODES];

	double VandermondeM [_QUADRATURE_NODES][_DOF1D];
	double derivativeVandermondeM [_QUADRATURE_NODES][_DOF1D];

	double** DF[5][_QUADRATURE_NODES][_QUADRATURE_NODES];
	double** DP;
	double** DQ;

	unsigned *boundary_nodes;
	unsigned **FEtoDOF;

	double* det_DJe[_QUADRATURE_NODES][_QUADRATURE_NODES];//
	double*** Gepq[_QUADRATURE_NODES][_QUADRATURE_NODES];//[E][_DIM][_DIM]

	Mat	StiffnessM, MassM, boundaryStiffnessM;
	Vec b,f;

}GridData;

int init_GridData(unsigned M,double L,GridData *data);

int init_quadrature(GridData* data);
int init_VandermondeM(GridData *data);
int init_derivativeVandermondeM(GridData *data);

int init_boundary_nodes(GridData *data);
int init_FEtoDOF(GridData* data);
int init_GridData_Mat_Vec(GridData *data);
int init_D(GridData *data);
int init_detDJe(GridData *data);
int init_Gepq(GridData *data);

int free_GridData(GridData *data);

int free_boundary_nodes(GridData *data);
int free_FEtoDOF(GridData *data);
int free_GridData_Mat_Vec(GridData *data);
int free_D(GridData *data);
int free_detDJe(GridData *data);
int free_Gepq(GridData *data);

#endif
