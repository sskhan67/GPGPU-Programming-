extern "C"
 {
#include "cuda.h"
#include "cuda_runtime.h"
#include "stdlib.h"
#include "time.h"
#include "stdio.h"
#include "wrapper.h"
#include "PyC_types.h"
#include "stdlib.h"
#include <sys/time.h>
#define T 6


// Global variables
double *d_V1112,*d_Rcca1,*d_Ra2,*d_V1222,*d_Rcaa2,*d_Rc1;
double *h_H,*d_h;
double *d_H,*h_Hr,*d_Hr;





// working on it , sayed 
/*
__global__
void reduction(double *H,double *R)
{

extern __shared__ volatile double sdata[];
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(blockDim.x*2) + threadIdx.x;
    sdata[tid] = H[i] + H[i+blockDim.x];
    
    __syncthreads();
    // do reduction in shared mem
    for (unsigned int s=blockDim.x/2; s>32; s>>=1) {
        if (tid < s) sdata[tid] += sdata[tid + s];
        __syncthreads();
    }

    if (tid < 32)
    {
        sdata[tid] += sdata[tid + 32];
        sdata[tid] += sdata[tid + 16];
        sdata[tid] += sdata[tid + 8];
        sdata[tid] += sdata[tid + 4];
        sdata[tid] += sdata[tid + 2];
        sdata[tid] += sdata[tid + 1];
    }
    // write result for this block to global mem
    if (tid == 0)
	{ 
	R[blockIdx.x] = sdata[0];
	//printf(" %lf", sdata[0]);
	}
}


*/

__global__
void reduction(double *H,double *R)
{

extern __shared__  double sdata[];
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
    sdata[tid] = H[i];

    __syncthreads();

	for (unsigned int s=blockDim.x/2; s>0; s>>=1)
	{
	if (tid<s)
	{
	sdata[tid]+=sdata[tid+s];
	}
	__syncthreads();
	}
	if (tid==0){R[blockIdx.x]=sdata[0];}
}








__global__
void outerloop(double *V1112,double *Rcca1,double *Ra2,double *H,int n_orb1,int n_orb2,double *V1222,double *Rcaa2, double *Rc1,double *h)
{
 int p1=threadIdx.x+blockIdx.x*blockDim.x;
 int q1=threadIdx.y+blockIdx.y*blockDim.y;
 int r1=threadIdx.z+blockIdx.z*blockDim.z;

 double Hlocal=0;

 if (p1 <n_orb1 && q1 <n_orb1 && r1<n_orb1)
 {	for (int s2=0;  s2<n_orb2;  s2++)
        {
//upperloop
	  Hlocal += V1112[((p1*n_orb1 + q1)*n_orb1 + r1)*n_orb2 + s2] * Rcca1[(q1*n_orb1 + p1)*n_orb1 + r1] * Ra2[s2];
//middleloop
	    Hlocal+=V1222[((p1*n_orb2 + q1)*n_orb2 + r1)*n_orb2 + s2] *     Rc1[p1] * Rcaa2[(q1*n_orb2 + s2)*n_orb2 + r1];
	  
        }
//bottomloop	
 	H[(p1*n_orb1+q1)*n_orb1+r1]=2*Hlocal;
        H[(p1*n_orb1+q1)*n_orb1+r1]+=(r1==0)?(h[p1*n_orb2 + q1] * Rc1[p1] * Ra2[q1]):0;
 }	
//reduction still performed externally

}

//exploring that n_orb1 and n_orb2 will be same


//TODO: Syed performs cudaMallocs just once (optimization 1, save and push file as cutest_s1.cu). Syed performs reduction inside kernel (optimization 2, which builds on 1,cutest_s2.cu).
//TODO: Thor performs cudastreams on this version and then later includes optimizations 1 and 2 from Syed. Thor also works on his own version for comparison.


//VKP: Skeleton code

double test_wrapper(int n_orb1, int n_orb2, double* Rc1, double* Rcca1,  double* Ra2,  double* Rcaa2,  double* h,  double* V1112,  double* V1222,int freevaribales)

{	


if (freevaribales==0) // will not execute cudaFree and host free 
{

cudaError_t cudaResult;
struct timeval start,stop,gpustart,gpustop;

int N4=n_orb1*n_orb1*n_orb1*n_orb2;
int N3=n_orb1*n_orb1*n_orb1;


int i=0;
gettimeofday(&gpustart,NULL);

	
dim3 dimblock(T,T,T);
int nbpgrid=n_orb1/T;
dim3 dimgrid(nbpgrid,nbpgrid,nbpgrid);
int blocks=(n_orb1*n_orb1*n_orb1)/(T*T);
dim3 dimblockR(T*T);
dim3 dimgridR(blocks);


	// counter-> assign dynamic memory only once at count=0
static int count=0;


if(count==0)
{

cudaMalloc((void **)&d_V1112,sizeof(double)*N4);
cudaMalloc((void **)&d_Rcca1,sizeof(double)*N3);
cudaMalloc((void **)&d_Ra2,sizeof(double)*n_orb2);

cudaMalloc((void **)&d_V1222,sizeof(double)*n_orb1*n_orb2*n_orb2*n_orb2);
cudaMalloc((void **)&d_Rcaa2,sizeof(double)*n_orb2*n_orb2*n_orb2);
cudaMalloc((void **)&d_Rc1,sizeof(double)*n_orb1);


h_Hr=(double *)malloc(sizeof(double)*blocks);
h_H=(double *)malloc(sizeof(double)*n_orb1*n_orb1*n_orb1);

cudaMalloc((void **)&d_H,sizeof(double)*n_orb1*n_orb1*n_orb1);
cudaMalloc((void **)&d_h,sizeof(double)*n_orb1*n_orb2);
cudaMalloc((void **)&d_Hr,sizeof(double)*blocks);

}


gettimeofday(&start, NULL);

cudaMemcpy(d_V1112,V1112,sizeof(double)*N4,cudaMemcpyHostToDevice);
cudaMemcpy(d_Rcca1,Rcca1,sizeof(double)*N3,cudaMemcpyHostToDevice);
cudaMemcpy(d_Ra2,Ra2,sizeof(double)*n_orb2,cudaMemcpyHostToDevice);

cudaMemcpy(d_V1222,V1222,sizeof(double)*n_orb1*n_orb2*n_orb2*n_orb2,cudaMemcpyHostToDevice);
cudaMemcpy(d_Rcaa2,Rcaa2,sizeof(double)*n_orb2*n_orb2*n_orb2,cudaMemcpyHostToDevice);
cudaMemcpy(d_Rc1,Rc1,sizeof(double)*n_orb1,cudaMemcpyHostToDevice);
cudaMemcpy(d_h,h,sizeof(double)*n_orb1*n_orb2,cudaMemcpyHostToDevice);

cudaResult = cudaGetLastError();
if (cudaResult != cudaSuccess)
{
	printf("\n H2D failed...");
	printf("%s",cudaGetErrorString(cudaResult));
}


 outerloop<<<dimgrid,dimblock>>>(d_V1112,d_Rcca1,d_Ra2,d_H,n_orb1,n_orb2,d_V1222,d_Rcaa2,d_Rc1,d_h);
 

cudaResult = cudaGetLastError();
if (cudaResult != cudaSuccess)
{
	printf("\n Outerloop failed...");
}



	// working on reduction-> sayed	

//printf("\n reduction threads:%d,blocks:%d",T*T,dimgridR.x);
	reduction<<<dimgridR,dimblockR,sizeof(double)*T*T>>>(d_H,d_Hr);	


 cudaResult = cudaGetLastError();
if (cudaResult != cudaSuccess)
{
	printf("\n reduction failed...");
}

	// external reduction  
cudaMemcpyAsync(h_H,d_H,sizeof(double)*n_orb1*n_orb1*n_orb1,cudaMemcpyDeviceToHost,streams[0]);
	// GPU reduction 
cudaMemcpyAsync(h_Hr,d_Hr,sizeof(double)*blocks,cudaMemcpyDeviceToHost,streams[1]);


 cudaStreamSynchonize(streams[0]);
 cudaStreamSynchonize(streams[1]);
if (cudaResult != cudaSuccess)
{
	printf("\n D2H failed...");
}

	// gpu reduction 
double sum=0;
for(i=0;i<blocks;i++)
        sum+=h_Hr[i];

printf("Reduction  H: %lf \n", sum);


	// external reduction 

double H = 0;

  for(i=0;i<n_orb1*n_orb1*n_orb1;i++)
	H+=h_H[i];

printf("External reduction H : %lf \n",H);
printf("error: %lf \n", H-sum);

gettimeofday(&stop,NULL);
 double comptime=(double)(stop.tv_sec-start.tv_sec)*1000+(double)(stop.tv_usec-start.tv_usec)/1000;


gettimeofday(&gpustop,NULL);
double gputime=(double)(gpustop.tv_sec-gpustart.tv_sec)*1000+(double)(gpustop.tv_usec-gpustart.tv_usec)/1000;


 gettimeofday(&stop,NULL);

// double cpucomptime=(double)(stop.tv_sec-start.tv_sec)*1000+(double)(stop.tv_usec-start.tv_usec)/1000;

printf("\n GPU computation time: %lf ms GPU end to end: %lf ms",comptime,gputime);

count++;

return  H;

}


	// free dynamic memory 
else 
	{
	
	cudaFree(d_V1112);
	cudaFree(d_Rcca1);
	cudaFree(d_Ra2);
	cudaFree(d_V1222);
	cudaFree(d_Rcaa2);
	cudaFree(d_Rc1);
	cudaFree(d_H);
	cudaFree(d_h);
	cudaFree(d_Hr);
	free(h_Hr);
	free(h_H);

	return 0;

	}





}




}
