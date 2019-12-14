#include "transformation.h"
#include <math.h>

//F_e = F_s * P_(i,j)
//i -> y coordinate
void indices(unsigned e,unsigned *i,unsigned *j,unsigned *s, GridData *data){
	unsigned M=data->M;
	*s=e/(M*M);
	e=e%(M*M);
	*i=e/M;
	*j=e%M;
}
/* Implementation of coefficients needed for transformation */
void coefficients(double x, double y, double *theta, double *c, double *r,GridData *data){
	//linear transformation of y to the angle theta
	*theta=M_PI/4*y;

	//distance to the square boundary from the center along a ray with angle theta
	*c=data->L/(2*cos(*theta));

	//linear transformation of [-1,1] -> [c,1]
	*r = (x*(1-*c)+(1+*c))/2;
}

void coordinates_dof(unsigned dof, double *x,double *y, GridData *data){
	unsigned i,j,outer_dof,number_of_rot; //i -> y coordinate, j-> x coordinate
	double theta,c,r;
	if(dof<(data->M+1)*(data->M+1)){
		//center square
		i = dof / (data->M+1);//ycooord
		j = dof % (data->M+1);//xcoord
		*x = -1 + (2. * j) / data->M;//from left to right
		*y = 1 - (2. * i) / data->M;//from up to down
		//apply F_0
		*x *= (data->L/2);
		*y *= (data->L/2);
		return;
	}else{
		outer_dof=dof- (data->M+1)*(data->M+1);//exclude inner square dofs
		number_of_rot = outer_dof / (data->M*data->M); //there are m^2 dofs attributed to each outer corse element
		j = (outer_dof % data->M) + 1; // all dof on the outer ring are shifted, since these are already counted by the center square
		i = (outer_dof/data->M)%data->M; //first operation gives the index in all outer rows, mod gives index in element
		*x = -1 + (2. * j) / data->M;//same as above
		*y = 1 - (2. * i) / data->M;
	}

	//apply F_1 to the coordinates, TODO: auslagern
	coefficients(*x,*y,&theta,&c,&r,data);
	*x = r * cos(theta);
	*y = r * sin(theta);

	//rotate to get F2,F3,F4
	while(number_of_rot){
		rotate(x,y);
		number_of_rot--;
	}

}
/* Implement P_ij, scale down x,y and shift to finite element (i,j)*/
void project(double *x,double *y, unsigned i, unsigned j, GridData *data){
	*x = *x/data->M;
	*y = *y/data->M;
	*x = *x + (-1 + (2. * j + 1)/data->M);
	*y = *y + (1 - (2. * i + 1)/data->M);
}
/* 90 degree rotation*/
void rotate(double *x, double *y){
	double temp = *x;
	*x = *y;
	*y = -temp;
}
