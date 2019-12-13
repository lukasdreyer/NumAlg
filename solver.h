#ifndef SOLVER_H_
#define SOLVER_H_
#include <petscksp.h>
#include "griddata.h"


int solve(Vec x,GridData *data);
int init_rhs(Vec *b,GridData *data);

#endif /* SOLVER_H_ */
