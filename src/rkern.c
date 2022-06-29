#include <R.h>

// R CMD SHLIB rkern.c  


void rkern(double *x1, double *x2, double *Kernel,  int *Nrow, int *Ncol, int *Nfeat)
{	
	int j,k, f, kn, kf, jf;
	double x_sum, sum_sum_sq=0.0;

 	for (k=0;k<*Nrow;k++){//Nrow
		kf=k * *Nfeat;
		kn=k * *Ncol;
		for (j=0;j<*Ncol;j++){//Ncol
			jf=j * *Nfeat; //added if(pcurrent->dtheta[f]>0) to not bother with theta=0
			for (f=0;f<*Nfeat; f++){
								x_sum=x1[kf+f]-x2[jf+f];
							//	 Rprintf(" %f  %i\n",x1[kf+f], x2[jf+f]);
								sum_sum_sq+=x_sum * x_sum;
								}
			Kernel[kn+j]=sum_sum_sq;
			sum_sum_sq=0.0;
		}
	}
}

