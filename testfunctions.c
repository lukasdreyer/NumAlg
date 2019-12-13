#include "testfunctions.h"
#include <math.h>
#include "griddata.h"
double constantA(double x,double y,GridData *data){
	return data->tfdata.A;
}
double linearAx(double x,double y,GridData *data){
	return data->tfdata.A*x;
}
double bilinearAxpBy(double x,double y,GridData *data){
	return data->tfdata.A*x+data->tfdata.B*y;
}

double radius(double x,double y,GridData *data){
	return sqrt(x*x+y*y);
}
double radiussquared(double x,double y,GridData *data){
	return x*x+y*y;
}

double zero_fct(double x,double y,GridData *data){
	return 0;
}
