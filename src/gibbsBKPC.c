/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


R CMD SHLIB gibbsBKPC.c  matrices.c 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h> 

#include "matrices.h"



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               Define a structure of current parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/



struct cparameters {
       int      Nred, *dbeta, thin,  N,  N_ITER,  NM; 
       double   *y, g1, g2, g4, g3, *acceptz, sd, sigmasq, rt, sh1, sh2, *sum_beta_K,  *Km, *KPC, 
          *L,*Li, *V, *Vinv, *KtK,*z,*znew,*zold, *tau, *beta, *Q, Q2, *M, m, *Kz, *Vm,
          *initBeta, *initTau, *betas,  *taus, *zs, *sigmasqs;
          };



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                   FUNCTION DECLARATIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


void mainbkpc(double *design,  int *input_y, int *nred, int *thin, int *ntr,  
int *niter, int *nm, double *sd, 
double *g1, double *g2, double *g3, double *g4, double *initSigmasq, double *initBeta, double *initTau, 
double *betas, double *taus, double *zs, double *sigmasqs);

		
void allocate_parameters(struct cparameters *pcurrent);
//void free_parameters(struct cparameters *pcurrent);	
void allocate_data(struct cparameters *pcurrent);
//void free_data(struct cparameters *pcurrent);


void print_current(struct cparameters *pcurrent, int i);


void initialize_parameters(struct cparameters *pcurrent);
void initialize_sum_beta(struct cparameters *pcurrent);




void gibbs_iterations(struct cparameters *pcurrent);
void metropolis_step(struct cparameters *pcurrent, int j);

void p_z(double *log_prob, double *current_z, struct cparameters *pcurrent, int current_j);
void update_matrices(struct cparameters *pcurrent);
void update_beta(struct cparameters *pcurrent);
void update_sigmasq(struct cparameters *pcurrent);
void update_tau(struct cparameters *pcurrent);


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                      MAIN
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/



/* main loop */

void mainbkpc(double *design,  int *input_y, int *nred, int *thin, int *ntr,  int *niter,  int *nm, double *sd, 
double *g1, double *g2, double *g3, double *g4, double *initSigmasq, double *initBeta, double *initTau,  double *betas, 
double *taus, double *zs, double *sigmasqs){
 
  int i;
	
  /* Initialize structure */	
  struct cparameters current, *pcurrent;
  pcurrent=&current;
  
  /* arguments from command line*/ 
  current.Nred=*nred;
  current.thin=*thin;
  current.N=*ntr;                         
  current.N_ITER=*niter;                      	 
  current.NM=*nm;                             
  current.sd=*sd;

  current.g1=*g1;
  current.g2=*g2;
  current.g3=*g3;
  current.g4=*g4;

  current.initBeta=initBeta;
  current.initTau=initTau;
  
  current.betas=betas;
  current.taus=taus;
  current.zs=zs;
  current.sigmasqs=sigmasqs;
          
  /* Initialize parameters that never change */ 
  current.sh1=current.g1+0.5*current.N*current.NM;
  current.sh2=current.g3+0.5;

  
  
  allocate_data(pcurrent);
  for (i=0; i<current.NM*current.N;i++)current.y[i]=input_y[i];   
  for (i=0; i<current.N*current.Nred;i++)current.KPC[i]=design[i];   
  current.sigmasq= *initSigmasq;//THIS CHANGES LATER
 
  GetRNGstate();
  allocate_parameters(pcurrent);         
  gibbs_iterations(pcurrent);                                     	
 // free_parameters(pcurrent); 
  PutRNGstate();
  
//  free_data(pcurrent);          

}





/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                INITIALIZE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
 /*initialize sum_beta_K and z*/
 
 
void initialize_sum_beta(struct cparameters *pcurrent)
{
    	int m, row,rowred, j;
     
    	for (j=0; j<pcurrent->N;j++)(pcurrent->acceptz)[j]=0.0;
    	
    	for (m=0; m<pcurrent->NM; m++){
            row=m*pcurrent->N;
            rowred=m*pcurrent->Nred;

            for (j=0; j<pcurrent->N;j++)(pcurrent->z)[row+j]=((pcurrent->y)[row+j]-0.5)*6.*pcurrent->sigmasq;
            for (j=0; j<pcurrent->Nred;j++){
                                     (pcurrent->beta)[rowred+j]=(pcurrent->initBeta)[rowred+j];
        	                         (pcurrent->tau)[rowred+j]=(pcurrent->initTau)[rowred+j];	
                                     }
		   
           matrix_by_vector(&(pcurrent->Km)[row*pcurrent->Nred], &(pcurrent->beta)[rowred], &(pcurrent->sum_beta_K)[row], pcurrent->N ,pcurrent->Nred);
        }
}

 /*initialize parameters*/
 
 
void initialize_parameters(struct cparameters *pcurrent)
{ 
    int m, i;
	for(m=0;m<pcurrent->NM;m++)for(i=0;i<pcurrent->N*pcurrent->Nred;i++)(pcurrent->Km)[m*pcurrent->N*pcurrent->Nred+i]=pcurrent->KPC[i]; //copy KPCs into Km   
 	initialize_sum_beta(pcurrent);//initialize sum_beta_K and z
 	print_current(pcurrent, 0);
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                MCMC 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/



  /*Gibbs iterations*/
  
  
void gibbs_iterations(struct cparameters *pcurrent)
{
  int i,j, m;
  
 initialize_parameters(pcurrent);
 for (i=1;i<pcurrent->N_ITER;i++)
  {      
      R_CheckUserInterrupt();                 
      for (j = 0; j < pcurrent->N; j++)metropolis_step(pcurrent, j);//update z
	  update_matrices(pcurrent);//update theta		  
	  update_sigmasq(pcurrent); //update sigmasq
	  update_beta(pcurrent); //update beta 
	  for (m=0; m<pcurrent->NM; m++)matrix_by_vector(&(pcurrent->Km)[m*pcurrent->N*pcurrent->Nred], &(pcurrent->beta)[pcurrent->Nred*m], &(pcurrent->sum_beta_K)[pcurrent->N*m], pcurrent->N ,pcurrent->Nred);//update sum_beta_K: 	 
	  update_tau( pcurrent); //update tau	     
	  if( i % pcurrent->thin == 0 )print_current(pcurrent, i/pcurrent->thin);
     }	
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                               METROPOLIS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


/* function for a metropolis step */


void metropolis_step(struct cparameters *pcurrent, int j)
{
	double u,r, log_new_prob, log_old_prob;
	int m;
      	
       	for (m=0; m<pcurrent->NM; m++){
           	pcurrent->Q2=rnorm(0,1)*pcurrent->sd;
           	pcurrent->m=(pcurrent->z)[j+pcurrent->N*m];
           	(pcurrent->znew)[m]=pcurrent->Q2+pcurrent->m;
           	(pcurrent->zold)[m]=(pcurrent->z)[j+pcurrent->N*m];
           	}
        
	    p_z(&log_old_prob, pcurrent->zold, pcurrent, j);
	    p_z(&log_new_prob, pcurrent->znew, pcurrent, j); 
      
	    r=exp(log_new_prob-log_old_prob);
	    u=runif(0,1); 
	    if(u<r){
                pcurrent->acceptz[j]=pcurrent->acceptz[j]+1.;
                for(m=0; m<pcurrent->NM; m++)(pcurrent->z)[j+pcurrent->N*m]=(pcurrent->znew)[m];
                }
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                               METROP STEP: Z
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

//current_z is a vector 1x(NM-1)

void p_z(double *log_prob, double *current_z, struct cparameters *pcurrent, int j)
{
    int m, jrow;
    double pz=0., sum_exp_z=0.;
    for (m=0; m<pcurrent->NM; m++){
	    jrow=pcurrent->N*m+j;	
        pz+=current_z[m]*(pcurrent->y)[jrow]-0.5*(current_z[m]-(pcurrent->sum_beta_K)[jrow])*(current_z[m]-(pcurrent->sum_beta_K)[jrow])/pcurrent->sigmasq;
        sum_exp_z+=exp(current_z[m]);
        }
	*log_prob=(pz-log(1+sum_exp_z));  
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                               METROP STEP: THETA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

void update_matrices( struct cparameters *pcurrent)
{
	int j,m, row, rown, rownred, rowred;
	double z_z=0.0, mVm=0.0;
	for (m=0; m<pcurrent->NM;m++)
        {          	  
        row=m*pcurrent->N;
        rowred=m*pcurrent->Nred;
        rown=m*pcurrent->Nred*pcurrent->Nred;
        rownred=m*pcurrent->N*pcurrent->Nred;
        matrix_by_matrix(&(pcurrent->Km)[rownred], pcurrent->KtK, pcurrent->N ,pcurrent->Nred);
        for (j=0;j<pcurrent->Nred*pcurrent->Nred;j++)(pcurrent->Vinv)[rown+j]=(pcurrent->KtK)[j];
        for (j=pcurrent->Nred;j--;)(pcurrent->Vinv)[rown+j*pcurrent->Nred+j]+=(pcurrent->tau)[rowred+j];// Vinv=K'K+T
        t_matrix_by_vector(&(pcurrent->Km)[rownred], &(pcurrent->z)[row], &(pcurrent->Kz)[rowred], pcurrent->N ,pcurrent->Nred);

        cholesky(&(pcurrent->Vinv)[rown], &(pcurrent->L)[rown], pcurrent->Nred); //decompose Vinv=LL'
        LInv(&(pcurrent->L)[rown], &(pcurrent->Li)[rown], pcurrent->Nred); //find Li=inv(L)
	          	
        //calculte mean vector M and matrix V,  where VCOV=sigmasqV
        matrix_by_matrix(&(pcurrent->Li)[rown], &(pcurrent->V)[rown],pcurrent->Nred,pcurrent->Nred); //V=Li'Li
        matrix_by_vector(&(pcurrent->V)[rown], &(pcurrent->Kz)[rowred], &(pcurrent->M)[rowred],pcurrent->Nred,pcurrent->Nred); //M=VKo'z
        matrix_by_vector(&(pcurrent->Vinv)[rown], &(pcurrent->M)[rowred], &(pcurrent->Vm)[rowred],pcurrent->Nred,pcurrent->Nred);

        cholesky(&(pcurrent->V)[rown], &(pcurrent->L)[rown], pcurrent->Nred);
        }

	for (j=pcurrent->NM*pcurrent->N;j--;)z_z+=(pcurrent->z)[j]*(pcurrent->z)[j]; 
	for (j=pcurrent->NM*pcurrent->Nred;j--;)mVm+=(pcurrent->M)[j]*(pcurrent->Vm)[j];
		
    pcurrent->rt=pcurrent->g2+0.5*(z_z-mVm);//Rprintf(" %f %f  %f\n",pcurrent->rt, z_z, mVm);
}




/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                               GIBBS INDIVIDUAL UPDATES: BETA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

void update_beta(struct cparameters *pcurrent)//draw RANDOM betas~MVN(M,sigmasqV)
{
	int m,j;
	double sigma=sqrt(pcurrent->sigmasq);

	for (m=0; m<pcurrent->NM;m++){
		for (j=pcurrent->Nred;j--;)(pcurrent->Q)[j]=rnorm(0, 1); //simulate 1xN Q~iid N(0,1)
		matrix_by_vector(&(pcurrent->L)[m*(pcurrent->Nred)*(pcurrent->Nred)], pcurrent->Q, &(pcurrent->beta)[m*pcurrent->Nred], pcurrent->Nred, pcurrent->Nred);//beta=LQ
	}

	for (j=pcurrent->NM*(pcurrent->Nred);j--;)(pcurrent->beta)[j]=pcurrent->beta[j]*sigma+(pcurrent->M)[j];//beta=M+sigmaLiQ
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                               GIBBS INDIVIDUAL UPDATES: sigmasq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/* directly depends on updates to pcurrent->rt which happens in the theta step. 
pcurrent->rt depends on z, m=VKz and Vinv=(KK+T).
Hence depends on z, theta(K), tau. SO this has to be updated b4 tau.*/


void update_sigmasq(struct cparameters *pcurrent)
{	
pcurrent->sigmasq = pcurrent->rt/rgamma((double)pcurrent->sh1, 1);
//	Rprintf(" %f %f %f \n",pcurrent->rt,(double)pcurrent->sh1, pcurrent->sigmasq); //forcing to stay smaller than 50
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                               GIBBS INDIVIDUAL UPDATES: TAU
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

//only  directly depends on updates to beta and sigma_sq


void update_tau(struct cparameters *pcurrent)
{
     double rate;
     int j;
     for (j=pcurrent->NM*(pcurrent->Nred);j--;){
         rate=pcurrent->g4+0.5*(pcurrent->beta)[j]*(pcurrent->beta)[j]/pcurrent->sigmasq;
         (pcurrent->tau)[j]=rgamma((double)pcurrent->sh2, 1/rate);
         }
}
     


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                      MEMORY RELATED
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


void allocate_data(struct cparameters *pcurrent)
{
	  pcurrent->y = (double *) R_alloc(pcurrent->N*pcurrent->NM, sizeof(double));
      pcurrent->KPC = (double *) R_alloc(pcurrent->N*pcurrent->Nred, sizeof(double));
}

//void allocate_data(struct cparameters *pcurrent)
//{
//	  pcurrent->y = (double *) Calloc(pcurrent->N*pcurrent->NM, double);
//      pcurrent->KPC = (double *) Calloc(pcurrent->N*pcurrent->Nred, double);
//}

//void free_data(struct cparameters *pcurrent)
//{
//	  free(pcurrent->y);
//	  free(pcurrent->KPC);
//}


void allocate_parameters(struct cparameters *pcurrent)
{
      int NR=pcurrent->Nred;
	  pcurrent->tau = (double *) R_alloc(pcurrent->NM*NR,sizeof(double));
	  pcurrent->beta = (double *) R_alloc(pcurrent->NM*NR,sizeof(double));
	  pcurrent->Km = (double *) R_alloc(pcurrent->NM*pcurrent->N*NR,sizeof(double));
	  pcurrent->KtK = (double *) R_alloc(NR*NR,sizeof(double));
	  pcurrent->Q = (double *) R_alloc(NR,sizeof(double));
	  pcurrent->M = (double *) R_alloc(pcurrent->NM*NR,sizeof(double));
	  pcurrent->Kz = (double *) R_alloc(pcurrent->NM*NR,sizeof(double));
	  pcurrent->Vm = (double *) R_alloc(pcurrent->NM*NR,sizeof(double));
	  pcurrent->L = (double *) R_alloc(pcurrent->NM*NR*NR,sizeof(double));
	  pcurrent->Li = (double *) R_alloc(pcurrent->NM*NR*NR,sizeof(double));		
	  pcurrent->V = (double *) R_alloc(pcurrent->NM*NR*NR,sizeof(double));
	  pcurrent->Vinv = (double *) R_alloc(pcurrent->NM*NR*NR,sizeof(double));
	  
      pcurrent->acceptz = (double *) R_alloc(pcurrent->N,sizeof(double));
	  pcurrent->z = (double *) R_alloc(pcurrent->NM*pcurrent->N,sizeof(double));
	  pcurrent->znew = (double *) R_alloc(pcurrent->NM,sizeof(double));
	  pcurrent->zold = (double *) R_alloc(pcurrent->NM,sizeof(double));
	  pcurrent->sum_beta_K = (double *) R_alloc(pcurrent->NM*pcurrent->N,sizeof(double));

}


//void allocate_parameters(struct cparameters *pcurrent)
//{
//      int NR=pcurrent->Nred;
//	  pcurrent->tau = (double *) Calloc(pcurrent->NM*NR,double);
//	  pcurrent->beta = (double *) Calloc(pcurrent->NM*NR,double);
//	  pcurrent->Km = (double *) Calloc(pcurrent->NM*pcurrent->N*NR,double);
//	  pcurrent->KtK = (double *) Calloc(NR*NR,double);
//	  pcurrent->Q = (double *) Calloc(NR,double);
//	  pcurrent->M = (double *) Calloc(pcurrent->NM*NR,double);
//	  pcurrent->Kz = (double *) Calloc(pcurrent->NM*NR,double);
//	  pcurrent->Vm = (double *) Calloc(pcurrent->NM*NR,double);
//	  pcurrent->L = (double *) Calloc(pcurrent->NM*NR*NR,double);
//	  pcurrent->Li = (double *) Calloc(pcurrent->NM*NR*NR,double);		
//	  pcurrent->V = (double *) Calloc(pcurrent->NM*NR*NR,double);
//	  pcurrent->Vinv = (double *) Calloc(pcurrent->NM*NR*NR,double);
//	  
//      pcurrent->acceptz = (double *) Calloc(pcurrent->N,double);
//	  pcurrent->z = (double *) Calloc(pcurrent->NM*pcurrent->N,double);
//	  pcurrent->znew = (double *) Calloc(pcurrent->NM,double);
//	  pcurrent->zold = (double *) Calloc(pcurrent->NM,double);
//	  pcurrent->sum_beta_K = (double *) Calloc(pcurrent->NM*pcurrent->N,double);
//
//}
	

//
//	  
//void free_parameters(struct cparameters *pcurrent)
//{	 
//      Free(pcurrent->acceptz);	  
//	  Free(pcurrent->z);
//	  Free(pcurrent->znew);
//	  Free(pcurrent->zold);
//	  Free(pcurrent->sum_beta_K); 
//    
//	  Free(pcurrent->tau);
//	  Free(pcurrent->beta);     
//     
//      Free(pcurrent->Km);
//	  Free(pcurrent->KtK);  
//	  Free(pcurrent->Q);
//	  Free(pcurrent->M);
//	  Free(pcurrent->Kz);
//      Free(pcurrent->Vm);
//	  Free(pcurrent->L);
//	  Free(pcurrent->Li);
//	  Free(pcurrent->V);
//	  Free(pcurrent->Vinv);
//}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                WRITE FILE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


void print_current(struct cparameters *pcurrent, int i)
{ 	
    int j;
    
   	pcurrent->sigmasqs[i] = pcurrent->sigmasq;
   	for (j=pcurrent->NM*pcurrent->Nred;j--;)pcurrent->betas[i * pcurrent->NM * pcurrent->Nred + j] = pcurrent->beta[j];
   	for (j=pcurrent->NM*pcurrent->Nred;j--;)pcurrent->taus[i * pcurrent->NM * pcurrent->Nred + j] = pcurrent->tau[j];
   	for (j=pcurrent->NM*pcurrent->N;j--;)pcurrent->zs[i * pcurrent->NM * pcurrent->N + j] = pcurrent->z[j];
    
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
