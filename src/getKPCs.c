#include <R.h>

#include "matrices.h"


// R CMD SHLIB getKPCs.c matrices.c 


void getKPCs(double *Ktild, double *Kscaled, double *E, double *KPC, int *Nrow, int *Nred ) //Ktild is scaled kernel
{   
    
    	int i,j;
    	eigen(Ktild, E, *Nrow); //Ktild is now a matrix of eigenvectors    
    
    
        for(i=0;i<*Nred;i++)for(j=0;j<*Nrow;j++)Ktild[i* *Nrow+j]=Ktild[i* *Nrow+j]/sqrt(E[i]/ *Nrow);
                                              
   
       	j=0; // creates KPC
    	for(i=0;i<*Nred;i++){
                             matrix_by_vectorM(Kscaled, &Ktild[i* *Nrow], &KPC[j], *Nrow, *Nrow, *Nred); 
                             ++j;
                             }
}


