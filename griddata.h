
/*
 ============================================================================
 Name        : griddata.h
 Author      : Lukas, Hannes
 Description : 
 ============================================================================
 */

#ifndef GRIDDATA_H
#define GRIDDATA_H

#define _QUADRATURE_NODES 2
#define _DOF1D 2
#define _DOF2D _DOF1D*_DOF1D
#define _DIM 2

#include<stdlib.h>

typedef struct GridData{
	double L;
	unsigned M;
	unsigned E;
	unsigned global_dof;
	double q_weights[_QUADRATURE_NODES];
	double q_nodes[_QUADRATURE_NODES];
	unsigned *boundary_nodes;
	unsigned **FEtoDOF;
	double VandermondeM [_DOF1D][_QUADRATURE_NODES];
	double derivativeVandermondeM [_DOF1D][_QUADRATURE_NODES];

}GridData;

int init_GridData(unsigned M,double L,GridData *data);

int init_boundary_nodes(GridData *data);
int init_FEtoDOF(GridData* data);
int init_quadrature(GridData* data);
int init_VandermondeM(GridData *data);

int free_FEtoDOF(GridData *data);
int free_GridData(GridData *data);



#endif

