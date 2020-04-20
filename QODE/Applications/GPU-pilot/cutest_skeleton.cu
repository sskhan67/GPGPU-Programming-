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
__global__
void reduction(double *H,double *R,int numele)
{
  __shared__ double partialSum[T*T];
  int t=threadIdx.x;
  int tid=threadIdx.x+blockIdx.x*blockDim.x;
 
  if(tid<numele) {
    partialSum[t]=H[tid];
 
    for(int stride=blockDim.x/2;stride>=1;stride/=2)
    {
    __syncthreads();
    if(t<stride)
          partialSum[t]+=partialSum[t+stride];
    } 
    if(t==0)
    R[blockIdx.x]=partialSum[t];
  } 
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

double test_wrapper(int n_orb1, int n_orb2, double* Rc1, double* Rcca1,  double* Ra2,  double* Rcaa2,  double* h,  double* V1112,  double* V1222)

{	
//if(n_orb1!=18 && n_orb2!=18)
//	printf("\n Different values: %d %d",n_orb1,n_orb2);

cudaError_t cudaResult;
struct timeval start,stop,gpustart,gpustop;
//printf("\n n_orb1:%d n_orb2:%d",n_orb1,n_orb2);

double *d_V1112,*d_Rcca1,*d_Ra2,*d_V1222,*d_Rcaa2,*d_Rc1;
double *d_H,*h_H,*d_h,*h_Hr,*d_Hr;

int N4=n_orb1*n_orb1*n_orb1*n_orb2;
int N3=n_orb1*n_orb1*n_orb1;
int i=0;
double H_gpu = 0.0;
gettimeofday(&gpustart,NULL);

dim3 dimblock(T,T,T);
int nbpgrid=n_orb1/T;
dim3 dimgrid(nbpgrid,nbpgrid,nbpgrid);
 int blocks=(n_orb1*n_orb1*n_orb1)/(T*T);
 dim3 dimblockR(T*T);
 dim3 dimgridR(blocks);


//Refactor such that these cudaMallocs are done just once because norbs never change

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
 //middleloop<<<dimgrid,dimblock>>>(d_V1222,d_Rcaa2,d_Rc1,d_H,n_orb1,n_orb2);

cudaResult = cudaGetLastError();
if (cudaResult != cudaSuccess)
{
    printf("\n Outerloop failed...");
}

//printf("\n reduction threads:%d,blocks:%d",T*T,dimgridR.x);
//reduction<<<dimgridR,dimblockR>>>(d_H,d_Hr,n_orb1*n_orb1*n_orb1);
 cudaResult = cudaGetLastError();
if (cudaResult != cudaSuccess)
{
    printf("\n reduction failed...");
}
//cudaMemcpy(h_Hr,d_Hr,sizeof(double)*blocks,cudaMemcpyDeviceToHost);
 cudaMemcpy(h_H,d_H,sizeof(double)*n_orb1*n_orb1*n_orb1,cudaMemcpyDeviceToHost);
 cudaDeviceSynchronize();

if (cudaResult != cudaSuccess)
{
    printf("\n D2H failed...");
}
 

double H = 0;
//for(i=0;i<blocks;i++)
  for(i=0;i<n_orb1*n_orb1*n_orb1;i++)
    H_gpu +=h_H[i];

 
 gettimeofday(&stop,NULL);
 double comptime=(double)(stop.tv_sec-start.tv_sec)*1000+(double)(stop.tv_usec-start.tv_usec)/1000;

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

gettimeofday(&gpustop,NULL);
double gputime=(double)(gpustop.tv_sec-gpustart.tv_sec)*1000+(double)(gpustop.tv_usec-gpustart.tv_usec)/1000;

printf("\n GPU computation time: %lf ms GPU end to end: %lf ms",comptime,gputime);
// printf("\n GPU printed: %lf in %lf ms",H,comptime);
gettimeofday(&start,NULL);
 H=0;	

    for (int p1=0;  p1<n_orb1;  p1++)
        {
        for (int q1=0;  q1<n_orb1;  q1++)
            {
            for (int r1=0;  r1<n_orb1;  r1++)
                {
                for (int s2=0;  s2<n_orb2;  s2++)
                    {
                    H += V1112[((p1*n_orb1 + q1)*n_orb1 + r1)*n_orb2 + s2] * Rcca1[(q1*n_orb1 + p1)*n_orb1 + r1] * Ra2[s2];
                    }
                }
            }
        }



for (int p1=0;  p1<n_orb1;  p1++)
        {
        for (int q2=0;  q2<n_orb2;  q2++)
            {
            for (int r2=0;  r2<n_orb2;  r2++)
                {
                for (int s2=0;  s2<n_orb2;  s2++)
                    {
                    H += V1222[((p1*n_orb2 + q2)*n_orb2 + r2)*n_orb2 + s2] * Rc1[p1] * Rcaa2[(q2*n_orb2 + s2)*n_orb2 + r2];
                    }
                }
            }
        }

    H *= 2;
    for (int p1=0;  p1<n_orb1;  p1++)
        {
        for (int q2=0;  q2<n_orb2;  q2++)
            {
            H += h[p1*n_orb2 + q2] * Rc1[p1] * Ra2[q2];
            }
        }

 gettimeofday(&stop,NULL);

  double cpucomptime=(double)(stop.tv_sec-start.tv_sec)*1000+(double)(stop.tv_usec-start.tv_usec)/1000;
// printf("\n CPU got: %lf in %lf ms",H,comptime);
printf("\n CPU: %lf GPU:%lf error: :%lf CPU time: %lf GPU time: %lf",H,H_gpu,(H_gpu-H),cpucomptime,comptime);
*/
return  H;


}




}
