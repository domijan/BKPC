#ifndef _MATRICES_H_
#define _MATRICES_H_

void transpose(double *M, int Nrow, int Ncol);
void matrix_by_vector(double *A, double *x, double *b, int Nrow, int Ncol);
void matrix_by_vectorM(double *A, double *x, double *b, int Nrow, int Ncol, int Ncolred); 
void t_matrix_by_vector(double *A, double *x, double *b, int Nrow, int Ncol); 
void matrix_by_matrix(double *M, double *A, int Nrow, int Ncol) ;


void eigen(double *K, double *E, int Nrow );
void cholesky(double *V, double *L, int Nrow);
void LInv(double *L, double *Li, int Nrow);



#endif
