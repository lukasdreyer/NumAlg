#ifndef MATSHELL_H
#define MATSHELL_H
#include <petscmat.h>



// positions 1,2,3 !!
#define TENSOR3MULT(tensor, v1, v2, v3, pos, vtensor, dim1, dim2, dim3, dimtensor, factor,tensor_erg) \
for(unsigned v1 = 0; v1 < dim1; v1++){\
	for(unsigned v2 = 0; v2 < dim2; v2++){\
		for(unsigned v3 = 0; v3 < dim3; v3++){\
			tensor_erg[v1][v2][v3]=0;\
			switch(pos){\
				case 1:\
				for(unsigned vtensor = 0; vtensor < dimtensor;vtensor++){\
					tensor_erg[v1][v2][v3]+=tensor[vtensor][v2][v3]*(factor);\
				}\
				break;\
				case 2:\
				for(unsigned vtensor = 0; vtensor < dimtensor;vtensor++){\
					tensor_erg[v1][v2][v3]+=tensor[v1][vtensor][v3]*(factor);\
				}\
				break;\
				case 3:\
				for(unsigned vtensor = 0; vtensor < dimtensor;vtensor++){\
					tensor_erg[v1][v2][v3]+=tensor[v1][v2][vtensor]*(factor);\
				}\
				break;\
				default:\
				printf("You entered the wrong position!!\n");\
			}\
		}\
	}\
}\

#define TENSOR2MULT(tensor, v1, v2, pos, vtensor, dim1, dim2, dimtensor, factor,tensor_erg) \
for(unsigned v1 = 0; v1 < dim1; v1++){\
	for(unsigned v2 = 0; v2 < dim2; v2++){\
		tensor_erg[v1][v2]=0;\
		switch(pos){\
			case 1:\
			for(unsigned vtensor = 0; vtensor < dimtensor;vtensor++){\
				tensor_erg[v1][v2]+=tensor[vtensor][v2]*(factor);\
			}\
			break;\
			case 2:\
			for(unsigned vtensor = 0; vtensor < dimtensor;vtensor++){\
				tensor_erg[v1][v2]+=tensor[v1][vtensor]*(factor);\
			}\
			break;\
			default:\
			printf("You entered the wrong position!!\n");\
		}\
	}\
}\

PetscErrorCode mass_mult(Mat A,Vec x,Vec y);
PetscErrorCode stiffness_mult(Mat A,Vec x,Vec y);
PetscErrorCode boundary_mult(Mat A,Vec x, Vec y);


#endif
