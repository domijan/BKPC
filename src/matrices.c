
#include "matrices.h"



#include <R.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Constant for Choleski 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#define TOL               0.0000000000000000000000001



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                      MATRICES - SIMPLE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*takes matrix M of size Nrow by Ncol and transposes it to a Ncol by Nrow matrix*/
void transpose(double *M, int Nrow, int Ncol)
{    
     int i, j;
     double *temp = (double *) malloc(sizeof(double)*Nrow*Ncol);
     
     memcpy (temp,M, Nrow*Ncol*sizeof(double));
     for (i=0; i<Ncol; i++)for (j=0; j<Nrow; j++)M[i*Nrow+j]=temp[i+Ncol*j];
     
     free(temp);
}


/*Multiplies a NrowxNcol matrix A by a Ncolx1 vector x and places this into a Nrowx1 vector b*/
void matrix_by_vector(double *A, double *x, double *b, int Nrow, int Ncol) 
{
 int j,k;
 for (k=Nrow;k--;){
	b[k]=0;
	for (j=Ncol;j--;)b[k]+=A[k*Ncol+j]*x[j];
 }
}

/*Multiplies a NrowxNcol matrix A by a Ncolx1 vector x and places this into a column  b of a matrix*/
void matrix_by_vectorM(double *A, double *x, double *b, int Nrow, int Ncol, int Ncolred) 
{
 int j,k;
 for (k=0;k<Nrow;k++){
	b[k*Ncolred]=0;
	for (j=0; j<Ncol;j++)b[k*Ncolred]+=A[k*Ncol+j]*x[j];
 }
}

/*Multiplies a transpose of a NrowxNcol matrix A by a Ncolx1 vector x and places this into a Nrowx1 vector b*/
void t_matrix_by_vector(double *A, double *x, double *b, int Nrow, int Ncol) 
{
 int j,k;
 for (k=Ncol;k--;){
	b[k]=0;
	for (j=Nrow;j--;)b[k]+=A[j*Ncol+k]*x[j];		
 }
}


/* A=M'M, where M is a Nrow x Ncol matrix and A is a Ncol x Ncol matrix*/
void matrix_by_matrix(double *M, double *A, int Nrow, int Ncol) 
{
  
	     char transa, transb;
	     double alpha, beta;
         int lda, ldb,ldc,  MM, N, K;


     

     transpose(M, Nrow, Ncol);

     transa='T';
     transb='N';
     MM=Ncol;
     N=Ncol;
     K=Nrow;
     
     lda = Nrow;
     ldb = Nrow;
     ldc= Ncol;
     alpha = 1.0;
     beta = 0.0;
      
     F77_CALL(dgemm)(&transa, &transb,&MM,&N,&K,&alpha, M,&lda,M,&ldb,&beta,A,&ldc);
     
     transpose(M, Ncol, Nrow);

}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                      MATRICES - MORE COMPLEX
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


// eigen decomposition of K and returns eigenvalues in E and eigenvectors in K
void eigen(double *K, double *E, int Nrow )
{

  int lda, n, info, lbuffer;
  char jobz, upper_triangle;
  
  upper_triangle = 'U';
  jobz = 'V';  
  n = Nrow;
  lda = Nrow;
  lbuffer = 3*Nrow-1;  
  
/* in FORTRAN/LAPACK/MATLAB world: a(i,j) = a[i+j*N] */
     double *buffer = (double *) malloc(sizeof(double)*lbuffer);    

     F77_CALL(dsyev)( &jobz, &upper_triangle, &n, K, &lda, E, buffer, &lbuffer, &info);

 //    printf("\n INFO=%d \n\n", info );  
     free(buffer);        

}



// decomposes Vinv=LL' and returns lower triangular matrix L  size Nrow x Nrow
void cholesky(double *V, double *L, int Nrow )
{

  int  i, j;
  int lda, n, info;
  char upper_triangle = 'U';
  
  n = Nrow;
  lda = Nrow;

  memcpy (L,V,Nrow*Nrow*sizeof(double));
  dpotrf_( &upper_triangle, &n, L, &lda, &info); //inverse of the inverse
  for(i = 0; i < Nrow; i++) {
    	for(j = i + 1; j < Nrow; j++) {
          L[i*Nrow + j] = 0;              
        }
  }    
}


// Inverse of lower triangular matrix L size Nrow x Nrow
void LInv(double *L, double *Li, int Nrow)
{
  double sum;
  int i,j,k, in, jn;

  for(i = 0; i < Nrow*Nrow; i++) Li[i]=L[i]; //init Li, above diag=0
  for(i = 0; i < Nrow; i++) {
	//diagonal
	in=i*Nrow;
	Li[in+i] = 1.0 / L[in+i]; 
	//below diagonal
	for(j = i + 1; j < Nrow; j++) {
		jn=j*Nrow;
         	sum = 0.0;
         	for(k = i; k < j; k++) sum -= Li[jn+k] * Li[k* Nrow+i];
         	Li[jn+i] = sum / L[jn+j];
      	}
  }
}


