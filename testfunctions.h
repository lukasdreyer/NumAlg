#ifndef TESTFUNCTIONS_H_
#define TESTFUNCTIONS_H_

#include <petscsystypes.h>
#include "griddata.h"

typedef double (*TestFunction) (double, double, GridData *);



double constantA(double x,double y,GridData *data);
double linearAx(double x,double y,GridData *data);
double bilinearAxpBy(double x,double y,GridData *data);
double zero_fct(double x,double y,GridData *data);

double radius(double x,double y,GridData *data);
double radiussquared(double x,double y,GridData *data);


//Implementier ableitungsfunktionen


#endif /* TESTFUNCTIONS_H_ */
