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

	double VandermondeM [2][_QUADRATURE_NODES][_DOF1D];//first entry 0: 1D Vandermonde, first entry 1:1D derivativeVandermonde
//	double derivativeVandermondeM [_QUADRATURE_NODES][_DOF1D];

	int *boundary_nodes;
	double *boundary_values;
	unsigned ***FEtoDOF;//[E][DIM][DIM]

	double* W[_QUADRATURE_NODES][_QUADRATURE_NODES];//[E]
	double*** DJ[_QUADRATURE_NODES][_QUADRATURE_NODES];//[E][_DIM][_DIM]
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
int init_W(GridData *data);
int init_Gepq(GridData *data);

int set_boundary_values_const(GridData *data,double val);

int free_GridData(GridData *data);

int free_boundary_nodes(GridData *data);
int free_FEtoDOF(GridData *data);
int free_GridData_Mat_Vec(GridData *data);
int free_D(GridData *data);
int free_W(GridData *data);
int free_Gepq(GridData *data);

#endif
