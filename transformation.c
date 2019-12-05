#include "transformation.h"
#include <math.h>

//F_e = F_s * P_(i,j)
void indices(unsigned e,unsigned *i,unsigned *j,unsigned *s, GridData *data){
	unsigned M=data->M;
	*s=e/(M*M);
	e=e%(M*M);
	*i=e/M;
	*j=e%M;
}
/* 90 degree rotation
 * */
void rotate(double *x, double *y){
	double temp = *x;
	*x = *y;
	*y = -temp;
}
/* Implement P_ij, scale down x,y and shift to finite element (i,j)
 * */
void project(double *x,double *y, unsigned i, unsigned j, GridData *data){

	*x = *x/data->M;
	*y = *y/data->M;
	*x = *x + (-1 + (2. * j + 1)/data->M);
	*y = *y + (1 - (2. * i + 1)/data->M);
}

/*
 * Implementation of coefficients needed for transformation
 *
 */
void coefficients(double x, double y, double *theta, double *c, double *r,GridData *data){
	*theta=M_PI/4*y;
	*c=data->L/(2*cos(*theta));
	*r = (x*(1-*c)+(1+*c))/2;
}



/* J_e = F_s (P_i,j(x,y))
 * F_s maps reference square to a portion of the unit circle
 * F_0 center
 * F_1 right
 * clockwise afterwards
 *
 * F_1 = (r cos theta, r sin theta) ^ (r, Id) ^ (Id,theta)
 */
void element_transformation(unsigned e,double x, double y,double *tx,double *ty, GridData *data){
	unsigned i,j,s;
	double c,r,theta;

	indices(e,&i,&j,&s,data);

	project(&x,&y,i,j,data);//P_i,j

	if(s==0){// s=0: scaling
		*tx= data->L/2*x;
		*ty= data->L/2*y;
		return;
	}

	coefficients(x,y,&theta,&c,&r,data);

	*tx = r*cos(theta);// polar transformation
	*ty = r*sin(theta);

//	printf("i: %i  j: %i s: %i \n",i,j,s);

	for(int k=1;k<s;k++){//rotate s-1 times
		rotate(tx,ty);
	}
}

double determinant_jacobian_transformation(unsigned e,double x, double y, GridData *data){
	unsigned i,j,s;
	double c,r,theta;

	unsigned M = data->M;
	double L = data->L;

	indices(e,&i,&j,&s,data);
	if(s==0) return (L*L)/(4.*M*M);
	project(&x,&y,i,j,data);

	coefficients(x,y,&theta,&c,&r,data);

	return r * ((1-c)/2) * (M_PI/4) / (M*M);
}

double geometryterm(unsigned p,unsigned q,unsigned alpha,unsigned beta,unsigned e,GridData *data){
	double c,r,theta,lambda,ret_value;
	unsigned i,j,s;
	double determinant;

	indices(e,&i,&j,&s,data);


	if(s==0){
		if((p+q) %2 == 0) return 4 * data->M *data->M* data->L *data->L;
		else return 0;
	}
	coefficients(data->q_nodes[alpha],data->q_nodes[beta],&theta,&c,&r,data);

	lambda = data->L * M_PI *(1-data->q_nodes[alpha])*sin(theta)/(16*cos(theta)*cos(theta));

	determinant = determinant_jacobian_transformation(e,data->q_nodes[alpha],data->q_nodes[beta],data);
	if(p==1&&q==1)ret_value = (1-c)*(1-c)/4;
	if(p==0&&q==0)ret_value = lambda * lambda + M_PI * M_PI * r * r / 16 - M_PI * r *cos(theta)*sin(theta);
	if(p*q==0 && p+q==1) ret_value = (1-c)/2*(-lambda+M_PI/2 *r*cos(theta)*sin(theta));

	return ret_value/(determinant*determinant)*data->M*data->M;
}
