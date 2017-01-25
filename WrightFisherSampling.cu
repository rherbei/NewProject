/******************************************************************************************
**Code to generate samples from a Wright-Fisher diffusion with selection
**For each set of parameters, repeatedly proposes N values until total of SAMPSIZE accepted
**Set T, theta, range of s, and range of p at start of main program
**Set name of output file at end of main program
**Output file contains SAMPSIZE draws for each combination of s and p  
**Compile with...
/usr/local/cuda/bin/nvcc -arch sm_20 -L /usr/local/cuda/lib64 -lcurand -o wf WrightFisherSampling.cu `pkg-config --cflags --libs gsl`
**Execute with...   ./wf &
******************************************************************************************/

#include <stddef.h>  // NULL, size_t 
#include <math.h> // expf 
#include <stdio.h> // printf
#include <time.h> // time 
#include <assert.h> 
#include <gsl/gsl_randist.h>

#include <curand.h> 
#include <cuda.h> 
#include <curand_kernel.h> 

#include <iostream> 
#include <fstream> 

using namespace std; 

#define N 25000
#define SAMPSIZE 10000

double phiplus(double s, double theta)
{
    double phimax;
    if (0<=0.5/s*(s-2.0*theta) && 0.5/s*(s-2.0*theta)<=1)
        phimax = 1.0/8.0*pow(s-2.0*theta,2)+theta*s/2;
    else if (0.5/s*(s-2*theta)<0)
        phimax=theta*s/2;
    else
        phimax=-theta*s/2;
    return phimax;
}

void quick_sort(double* arr,int* ord,int low,int high)

{

 int pivot,j,i,temp2;

 double temp;
 if(low<high)

 {
  
  pivot = low;

  i = low;

  j = high;


 
 while(i<j)

  {

   while((arr[i]<=arr[pivot])&&(i<high))

   {

    i++;

   }

 
  while(arr[j]>arr[pivot])

   {

    j--;

   }

 
  if(i<j)

   {
 
    temp=arr[i];
     arr[i]=arr[j];

     arr[j]=temp;

     temp2=ord[i];
     ord[i]=ord[j];

     ord[j]=temp2;

   }

  }
  temp=arr[pivot];

  arr[pivot]=arr[j];

  arr[j]=temp;

  temp2=ord[pivot];

  ord[pivot]=ord[j];

  ord[j]=temp2;
  quick_sort(arr,ord,low,j-1);

  quick_sort(arr,ord,j+1,high);

 }

}

__global__ void TimeDiffs(int *dev_csp, unsigned int *dev_poi, double *dev_parms, double *dev_tss, double *dev_tdiffs)
{
  int bid=blockIdx.x;
   int size, start, firstt;
   size=dev_poi[bid];
   firstt=dev_csp[bid];
   start=firstt+bid;
   if (size==0)
     dev_tdiffs[start]=dev_parms[2];
   else
   {
     dev_tdiffs[start]=dev_tss[firstt];
     for (int j=1; j<size; j++)
     {
       dev_tdiffs[start+j]=dev_tss[firstt+j]-dev_tss[firstt+j-1];
     }
     dev_tdiffs[start+size]=dev_parms[2]-dev_tss[firstt+size-1];
   }
}

__global__ void NormApproxParams(int *dev_csp, unsigned int *dev_poi, double *dev_tdiffs, double *dev_parms, double *dev_means, double *dev_vars)
{
  int bid=blockIdx.x;
   int size, start;
   size=dev_poi[bid]+1;
   start=dev_csp[bid]+bid;
   for (int j=0; j<size; j++)
   {
    double beta, eta;
    beta=0.5*(dev_parms[1]-1)*dev_tdiffs[start+j];
    eta=beta/(exp(beta)-1);
    dev_vars[start+j]=abs(2.0*eta/dev_tdiffs[start+j]*pow(eta+beta,2)*(1+eta/(eta+beta)-2.0*eta));
    dev_means[start+j]=2.0*eta/dev_tdiffs[start+j];
   }
}
   
__global__ void KeepChoice(int *dev_csp, unsigned int *dev_poi, double *dev_skels, double *dev_us, double *dev_parms, double *dev_psio, double *dev_endpts, int *dev_keep)
{
  int bid=blockIdx.x;
   int size, start;
   double phitilde;
   dev_keep[bid]=1;
   dev_endpts[bid]=dev_skels[dev_csp[bid+1]-1+bid];
   size=dev_poi[bid];
   start=dev_csp[bid];
   if(exp(dev_parms[0]*dev_endpts[bid]-max(0.0,dev_parms[0])) < dev_us[bid])
    dev_keep[bid]=0;
   for (int j=0; j<size-1; j++)
   {
    phitilde=dev_parms[0]/2.0*(-dev_parms[0]*pow(dev_skels[start+j],2)+dev_skels[start+j]*(dev_parms[0]-2.0*dev_parms[1])+dev_parms[1]);
    if (phitilde > dev_psio[start+j])
     dev_keep[bid]=0;
   }
}

int main( void ) 
{ 
 /*declare objects with fixed size*/
 double s, p, phi, samp[SAMPSIZE], params[3], endpts[N];
 int y, z, i, j, k, poisum, cspoi[N+1], keep[N], size, start, index, kept, loops;
 unsigned int poi[N];

 /*set up random number generators*/
 curandGenerator_t gen; 
 curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MTGP32); 
 curandSetPseudoRandomGeneratorSeed(gen,1234ULL); 
 const gsl_rng_type * rngT;
 gsl_rng_env_setup();
 rngT = gsl_rng_default;

 /*set parameters*/
 double T=0.1;
 double theta=0.00014;
 double ps[99];
 for (i=0; i<99; i++) 
 {
  ps[i]=((double)i+1.0)/99.0; 
 }
 double ss[59];
 for (i=0; i<1; i++) 
 {
  ss[i]=((double)i/2.0)-12.0; 
 }

 for (y=0; y<59; y++)
 {
  s=ss[y];
  phi=phiplus(s,theta);

  FILE *file3;    
  file3 = fopen("current_sb.txt","w");
   fprintf(file3, "%.1f\n", s);  
  fclose(file3); 

  for (z=0; z<99; z++)
  { 
   p=ps[z]; 

   FILE *file4;    
   file4 = fopen("current_pb.txt","w");
    fprintf(file4, "%.2f\n", p);  
   fclose(file4); 

   params[0]=s;
   params[1]=theta;
   params[2]=T;

   kept=0;
   loops=0;

   while (kept<SAMPSIZE)
   {

    //printf("%d\n",loops);
    //printf("%d\n",kept);
    
    /***sample from poisson point process***/
    //generate poisson sample
    unsigned int *dev_poi; 
     cudaMalloc( (void**)&dev_poi, N*sizeof(int) ); 
    curandGeneratePoisson(gen, dev_poi, N, phi*T); 
    cudaMemcpy(poi, dev_poi, N*sizeof(int), cudaMemcpyDeviceToHost);  

    //compute sum of poisson sample
    poisum=0;
    for (j=0; j<N; j++)
    {
     poisum=poisum+poi[j];
    }

    //generate uniform samples for ts and psis
    double *ts, *psis;
     ts = (double *) malloc(poisum*sizeof(double));
     psis = (double *) malloc(poisum*sizeof(double));
    double *dev_ts, *dev_psis;
     cudaMalloc( (void**)&dev_ts, poisum*sizeof(double) ); 
     cudaMalloc( (void**)&dev_psis, poisum*sizeof(double) );
    curandGenerateUniformDouble(gen, dev_ts, poisum); 
    cudaMemcpy(ts, dev_ts, poisum*sizeof(double), cudaMemcpyDeviceToHost); 
    cudaFree(dev_ts);  
    for (j=0; j<poisum; j++)
    {
     ts[j]=ts[j]*T;
    } 
    curandGenerateUniformDouble(gen, dev_psis, poisum); 
    cudaMemcpy(psis, dev_psis, poisum*sizeof(double), cudaMemcpyDeviceToHost); 
    cudaFree(dev_psis);  
    for (j=0; j<poisum; j++)
    {  
     psis[j]=psis[j]*phi; 
    }

    /***sort each of N sets of t and reorder N sets of psi accordingly***/
    //compute cumulative sum for indexing sets
    cspoi[0]=0;
    for (j=0; j<N; j++)
    {
     cspoi[j+1]=cspoi[j]+poi[j];
    }

    //use quicksort to sort ts and reorder psis 
    double *psis_ord;
     psis_ord = (double *) malloc(poisum*sizeof(double));  
    for (k=0; k<N; k++)
    { 
     size=poi[k];
     start=cspoi[k];
     if (size==1)
     {
      psis_ord[start]=psis[start];
     }
     if (size>1)
     {
      double *a;
       a = (double *) malloc(size*sizeof(double));
      int *o;
       o = (int *) malloc(size*sizeof(int));
      for (j=0; j<size; j++)
      {
       a[j]=ts[start+j];
       o[j]=j;
      }
      quick_sort(a,o,0,size-1);
      for (j=0; j<size; j++)
      {
       ts[start+j]=a[j];
       psis_ord[start+j]=psis[start+o[j]];
      }
      free(a);
      free(o);
     }
    }
    free(psis);

    /***compute differences among consecutive ts within each set***/
    double *tdiffs;
     tdiffs = (double *) malloc((N+poisum)*sizeof(double));
    double *dev_tdiffs, *dev_parms, *dev_tss;
     cudaMalloc( (void**)&dev_tdiffs, (N+poisum)*sizeof(double) );
     cudaMalloc( (void**)&dev_parms, 3*sizeof(double) );
     cudaMalloc( (void**)&dev_tss, poisum*sizeof(double) );
    int *dev_csp;
     cudaMalloc( (void**)&dev_csp, (N+1)*sizeof(int) );
    cudaMemcpy( dev_tss, ts, poisum*sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy( dev_csp, cspoi, (N+1)*sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy( dev_parms, params, 3*sizeof(double), cudaMemcpyHostToDevice );
    TimeDiffs<<<N,1>>>(dev_csp,dev_poi,dev_parms,dev_tss,dev_tdiffs);
    cudaMemcpy( tdiffs, dev_tdiffs, (N+poisum)*sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(dev_tss);
  
    /***compute parameters for normal approximation***/
    double *naMeans, *naVars;
     naMeans = (double *) malloc((N+poisum)*sizeof(double));
     naVars = (double *) malloc((N+poisum)*sizeof(double));
    double *dev_means, *dev_vars;
     cudaMalloc( (void**)&dev_means, (N+poisum)*sizeof(double) );
     cudaMalloc( (void**)&dev_vars, (N+poisum)*sizeof(double) );
    NormApproxParams<<<N,1>>>(dev_csp,dev_poi,dev_tdiffs,dev_parms,dev_means,dev_vars);
    cudaMemcpy(naMeans, dev_means, (N+poisum)*sizeof(double), cudaMemcpyDeviceToHost); 
    cudaMemcpy(naVars, dev_vars, (N+poisum)*sizeof(double), cudaMemcpyDeviceToHost);     
    cudaFree(dev_tdiffs);
    cudaFree(dev_means);
    cudaFree(dev_vars);
  
    /***generate draws from neutral wright-fisher process***/
    double *ms, *skels;
     ms = (double *) malloc((N+poisum)*sizeof(double));
     skels = (double *) malloc((N+poisum)*sizeof(double));
    long long *rms;
     rms = (long long *) malloc((N+poisum)*sizeof(long long));
    int *ls;
     ls = (int *) malloc((N+poisum)*sizeof(int));
    gsl_rng * r;
     r = gsl_rng_alloc (rngT);
    index=0;
    for (j=0; j<N+poisum; j++) 
    {
     if (tdiffs[j]>0.000001)
     {
      ms[j]=gsl_ran_gaussian(r,sqrt(naVars[j]));
      ms[j]=ms[j]+naMeans[j];
     }
     else
      ms[j]=(2.0*theta-1.0)/(exp(0.5*(2.0*theta-1.0)*0.000001)-1.0);
     rms[j]=round(ms[j]);
     if (j==cspoi[index]+index)
     {
      ls[j]=gsl_ran_binomial(r,p,rms[j]);
      index++;
     }
     else
      ls[j]=gsl_ran_binomial(r,skels[j-1],rms[j]);
     skels[j]=gsl_ran_beta(r,theta+(double)ls[j],theta+(double)rms[j]-(double)ls[j]);
    }
    gsl_rng_free(r);  
    free(naMeans);
    free(naVars);
    free(ms);
    free(rms);
    free(ls);

    /***choose whether to accept or reject each skeleton***/
    //generate random uniforms with which to make decisions
    double *dev_us;
     cudaMalloc( (void**)&dev_us, N*sizeof(double) ); 
    curandGenerateUniformDouble(gen, dev_us, N); 
  
    //identify endpoints of skeletons and which ones to accept
    int *dev_keep;
     cudaMalloc( (void**)&dev_keep, N*sizeof(int) );
    double *dev_skels, *dev_psio, *dev_endpts;
     cudaMalloc( (void**)&dev_skels, (N+poisum)*sizeof(double) );
     cudaMalloc( (void**)&dev_psio, poisum*sizeof(double) );
     cudaMalloc( (void**)&dev_endpts, N*sizeof(double) );
    cudaMemcpy( dev_skels, skels, (N+poisum)*sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy( dev_psio, psis_ord, poisum*sizeof(double), cudaMemcpyHostToDevice );
    KeepChoice<<<N,1>>>(dev_csp,dev_poi,dev_skels,dev_us,dev_parms,dev_psio,dev_endpts,dev_keep);
    cudaMemcpy(endpts, dev_endpts, N*sizeof(double), cudaMemcpyDeviceToHost); 
    cudaMemcpy(keep, dev_keep, N*sizeof(int), cudaMemcpyDeviceToHost);
    cudaFree(dev_poi);
    cudaFree(dev_csp);
    cudaFree(dev_keep);
    cudaFree(dev_skels);
    cudaFree(dev_parms);
    cudaFree(dev_psio);
    cudaFree(dev_endpts);
    cudaFree(dev_us);
    free(ts);
    free(psis_ord);
    free(tdiffs);
    free(skels);

    /***store accepted endpoints and keep track of number accepted***/
    for (j=0; j<N; j++)
    {
     if (keep[j]==1 && kept<SAMPSIZE)
     {  
      samp[kept]=endpts[j];
      kept++;
     }
    }
    loops++;
   }

   /***append final sample to output file***/
   FILE *file1;  
   file1 = fopen("final_samples_cuda.txt","a");
   for(j=0; j<SAMPSIZE; j++)
   { 
    fprintf(file1, "%.6f\n", samp[j]);  
   } 
   fclose(file1); 
  } 
 }  
 return 1; 
}
