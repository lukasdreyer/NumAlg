
#include <math.h>
#include "griddata.h"

int init_GridData(unsigned M,double L,GridData *data){
	data->M=M;
	data->E=5*M*M;
	data->L=L;
	data->global_dof=data->E+2*M+1; //5M^2+2M+1
	set_boundary_nodes(data);//allocates memory
	set_FEtoDOF(data);//allocates memory
	set_quadrature(data);
	return 0;
}

int free_FEtoDOF(GridData *data){
	for(int e=0; e < data->E;e++){
		free(data->FEtoDOF[e]);
	}
	free(data->FEtoDOF);
	return 0;
}

int free_GridData(GridData *data){
	free_FEtoDOF(data);
	free(data->boundary_nodes);
	return 0;
}


int set_FEtoDOF(GridData* data){
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
/*			data->FEtoDOF[i*M+j][0]=i*(M+1)+j;
			data->FEtoDOF[i*M+j][1]=i*(M+1)+j+1;
			data->FEtoDOF[i*M+j][2]=(i+1)*(M+1)+j;
			data->FEtoDOF[i*M+j][3]=(i+1)*(M+1)+j+1;
*/
			for(int s=1;s<=4;s++){//(M+1)^2-1
				/* 4 different quartercircle cuts. 
				fill boundary squares, starting right clockwise
				 * */
				for(unsigned k=0;k<_DOF1D;k++){
					for(unsigned l=0;l<_DOF1D;l++){
						data->FEtoDOF[s*M*M+ i*M + j][k*_DOF1D + l]= ((M+2)*M + (s-1)*M*M) +  (i+k)*M +j +l;
					}
				}
/*				data->FEtoDOF[i*M+j+s*M*M][0]=(M+2)*M + (s-1)*M*M + i*M + j;
				data->FEtoDOF[i*M+j+s*M*M][1]=(M+2)*M + (s-1)*M*M + i*M + j+1;
				data->FEtoDOF[i*M+j+s*M*M][2]=(M+2)*M + (s-1)*M*M + (i+1)*M + j;
				data->FEtoDOF[i*M+j+s*M*M][3]=(M+2)*M + (s-1)*M*M + (i+1)*M + j+1;
*/			}
			if(i==(M-1)){
				/* adjust boundary between F_1 and F_4
				 **/
				data->FEtoDOF[i*M+j+4*M*M][2]=(M+2)*M + j;
				data->FEtoDOF[i*M+j+4*M*M][3]=(M+2)*M + j +1;

			}
		}
		/* adjust boundary to center square
		 *
		 */
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

int set_boundary_nodes(GridData *data){
	if(_DOF2D!=4) return -1;
	data->boundary_nodes = malloc(sizeof(unsigned)*4*data->M);//there are 4 quartercircles->4m DOF on boundary
	for (unsigned i=0;i<4*data->M;i++){
		data->boundary_nodes[i]=data->M*(data->M+3+i);
	}
	return 0;
}
int set_quadrature(GridData* data){
	if(_DOF2D!=4) return -1;
	data->q_nodes[0]=-sqrt(1./3);
	data->q_nodes[1]=sqrt(1./3);
	data->q_weights[0]=1;
	data->q_weights[1]=1;
	return 0;
}