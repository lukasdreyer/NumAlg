/*
 * mat2x2.h
 *
 *  Created on: Dec 8, 2019
 *      Author: lukas
 */

#ifndef MAT2X2_H_
#define MAT2X2_H_

int invert_Mat2x2(double** M,double** M_inv);
int transpose_Mat2x2(double** M,double** M_t);
int MatMult_Mat2x2(double** A,double** B,double** C);
int fill_Mat2x2(double** M,double x00,double x01,double x10,double x11);
int free_Mat2x2(double **M);
int alloc_Mat2x2(double ***M_p);
int differ_Mat2x2(double **M1,double **M2);
int print_Mat2x2(double **M);
double determinant_Mat2x2(double **M);

#endif /* MAT2X2_H_ */
