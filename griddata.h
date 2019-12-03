
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
}GridData;

int set_boundary_nodes(GridData *data);
int set_FEtoDOF(GridData* data);
int set_quadrature(GridData* data);
int init_GridData(unsigned M,double L,GridData *data);
int free_FEtoDOF(GridData *data);
int free_GridData(GridData *data);





#endif

