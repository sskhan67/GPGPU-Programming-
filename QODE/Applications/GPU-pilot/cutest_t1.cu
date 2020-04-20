//cudaMalloc done just once. Streams only for memcpys. No reduction kerne


extern "C"
 {
#include "cuda.h"
#include "cuda_runtime.h"
#include "stdlib.h"
#include "time.h"
#include "stdio.h"
#include "PyC_types.h"
#include <sys/time.h>
#define T 6
#define numStreams 7 // 7 for t1


// cleaner error handling; just wrap cuda library calls with gpuErrchk(foo());
#define gpuErr(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
} 


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

// create global streams once and re-use them in each call to test_wrapper_t1, then destroy them later (these are called in the loop function)


// t2:  break input into chunks
double test_wrapper_t2(int n_orb1, int n_orb2, double* Rc1, double* Rcca1,  double* Ra2,  double* Rcaa2,  double* h,  double* V1112,  double* V1222)
{	
double H = 0;
/*
double H_gpu = 0;
struct timeval start,stop,gpustart,gpustop;
//if(n_orb1!=18 && n_orb2!=18)
//	printf("\n Different values: %d %d",n_orb1,n_orb2);

cudaError_t cudaResult;
//printf("\n n_orb1:%d n_orb2:%d",n_orb1,n_orb2);

double *d_V1112,*d_Rcca1,*d_Ra2,*d_V1222,*d_Rcaa2,*d_Rc1;
double *d_H,*h_H,*d_h,*h_Hr,*d_Hr;

int N4=n_orb1*n_orb1*n_orb1*n_orb2;
int N3=n_orb1*n_orb1*n_orb1;
int N2 = n_orb1 * n_orb2;
int i=0;
gettimeofday(&gpustart,NULL);
dim3 dimblock(T,T,T);

int nbpgrid=n_orb1/T;
dim3 dimgrid(numStreams/nbpgrid,numStreams/nbpgrid,numStreams/nbpgrid);
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

int s4 = N4 / numStreams;
int s3 = N3 / numStreams;
int s2 = N2 / numStreams;
cudaMemcpyAsync(d_Ra2,Ra2,sizeof(double)*n_orb2,cudaMemcpyHostToDevice, streams[0]);
cudaMemcpyAsync(d_Rc1,Rc1,sizeof(double)*n_orb1,cudaMemcpyHostToDevice, streams[1]);
for(int j=0; j<numStreams; j++) {
	int o4 = j*s4;
	int o3 = j*s3;
	int o2 = j*s2;
	cudaMemcpyAsync(&d_V1112[o4], &V1112[o4], s4*sizeof(double), cudaMemcpyHostToDevice, streams[j]);
	cudaMemcpyAsync(&d_V1222[o4], &V1222[o4], s4*sizeof(double), cudaMemcpyHostToDevice, streams[j]);
	cudaMemcpyAsync(&d_Rcca1[o3], &Rcca1[o3], s3*sizeof(double), cudaMemcpyHostToDevice, streams[j]);
	cudaMemcpyAsync(&d_Rcaa2[o3] ,&Rcaa2[o3], s3*sizeof(double), cudaMemcpyHostToDevice, streams[j]);
	cudaMemcpyAsync(&d_h[o2],     &h[o2],     s2*sizeof(double), cudaMemcpyHostToDevice, streams[j]);
}
for(int j=0; j<numStreams; j++) {
	int o4 = j*s4;
	int o3 = j*s3;
	int o2 = j*s2;
	outerloop<<<
}



cudaResult = cudaGetLastError();
if (cudaResult != cudaSuccess)
{
	printf("\n H2D failed...");
	printf("%s",cudaGetErrorString(cudaResult));
}

 for(int j=0; j<numStreams; j++) {
    cudaStreamSynchronize(streams[j]);
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
 cudaMemcpyAsync(h_H,d_H,sizeof(double)*n_orb1*n_orb1*n_orb1,cudaMemcpyDeviceToHost, streams[0]);
 cudaStreamSynchronize(streams[0]); // using cudaDeviceSynchronize here is unnecessary

if (cudaResult != cudaSuccess)
{
	printf("\n D2H failed...");
}
 

//for(i=0;i<blocks;i++)
  for(i=0;i<n_orb1*n_orb1*n_orb1;i++)
	H_gpu+=h_H[i];

 
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
H = H_gpu;
printf("\n GPU computation time: %lf ms GPU end to end: %lf ms\n",comptime,gputime);
// printf("\n GPU printed: %lf in %lf ms",H,comptime);
gettimeofday(&start,NULL);
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
printf("\n CPU got: %lf in %lf ms",H,cpucomptime);
printf("\n CPU: %lf GPU:%lf error: :%lf CPU time: %lf GPU time: %lf",H,H_gpu,(H_gpu-H),cpucomptime,comptime);
*/
return  H;
}
int count =0;
void preMalloc(int n_orb1, int n_orb2)
{
	
}
// t1 and s1 : Asynchronous memcpys, 1 stream per input array, streams and mallocs done once
double test_wrapper(int n_orb1, int n_orb2, double* Rc1, double* Rcca1,  double* Ra2,  double* Rcaa2,  double* h,  double* V1112,  double* V1222, int freevariables)
{	
if(freevariables)
{
    for(int i=0; i<numStreams; i++) {
            cudaStreamDestroy(streams[i]);
    }
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
    return 0.0;
}
if(!count) {
    preMalloc(n_orb1,n_orb2);
}
double H = 0.0;
struct timeval start,stop,gpustart,gpustop;
//if(n_orb1!=18 && n_orb2!=18)
//	printf("\n Different values: %d %d",n_orb1,n_orb2);

cudaError_t cudaResult;
//printf("\n n_orb1:%d n_orb2:%d",n_orb1,n_orb2);

int i=0;
gettimeofday(&gpustart,NULL);
dim3 dimblock(T,T,T);
int N4=n_orb1*n_orb1*n_orb1*n_orb2;
int N3=n_orb1*n_orb1*n_orb1;
int nbpgrid=n_orb1/T;
dim3 dimgrid(nbpgrid,nbpgrid,nbpgrid);
dim3 dimblockR(T*T);
int blocks=(n_orb1*n_orb1*n_orb1)/(T*T);
dim3 dimgridR(blocks);


//Refactor such that these cudaMallocs are done just once because norbs never change
gettimeofday(&start, NULL);
cudaMemcpyAsync(d_V1112,V1112,sizeof(double)*N4,cudaMemcpyHostToDevice, streams[0]);
cudaMemcpyAsync(d_Rcca1,Rcca1,sizeof(double)*N3,cudaMemcpyHostToDevice, streams[1]);
cudaMemcpyAsync(d_Ra2,Ra2,sizeof(double)*n_orb2,cudaMemcpyHostToDevice, streams[2]);
cudaMemcpyAsync(d_V1222,V1222,sizeof(double)*n_orb1*n_orb2*n_orb2*n_orb2,cudaMemcpyHostToDevice, streams[3]);
cudaMemcpyAsync(d_Rcaa2,Rcaa2,sizeof(double)*n_orb2*n_orb2*n_orb2,cudaMemcpyHostToDevice, streams[4]);
cudaMemcpyAsync(d_Rc1,Rc1,sizeof(double)*n_orb1,cudaMemcpyHostToDevice, streams[5]);
cudaMemcpyAsync(d_h,h,sizeof(double)*n_orb1*n_orb2,cudaMemcpyHostToDevice, streams[6]);

cudaResult = cudaGetLastError();
if (cudaResult != cudaSuccess)
{
	printf("\n H2D failed...");
	printf("%s",cudaGetErrorString(cudaResult));
}

 for(int j=0; j<numStreams; j++) {
    gpuErr(cudaStreamSynchronize(streams[j]));
 }
 outerloop<<<dimgrid,dimblock>>>(d_V1112,d_Rcca1,d_Ra2,d_H,n_orb1,n_orb2,d_V1222,d_Rcaa2,d_Rc1,d_h);
 //middleloop<<<dimgrid,dimblock>>>(d_V1222,d_Rcaa2,d_Rc1,d_H,n_orb1,n_orb2);
 gpuErr(cudaPeekAtLastError());

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
 cudaMemcpyAsync(h_H,d_H,sizeof(double)*n_orb1*n_orb1*n_orb1,cudaMemcpyDeviceToHost, streams[0]);
 cudaDeviceSynchronize();

if (cudaResult != cudaSuccess)
{
	printf("\n D2H failed...");
}
 

//for(i=0;i<blocks;i++)
  for(i=0;i<n_orb1*n_orb1*n_orb1;i++)
	H += h_H[i];

 
 gettimeofday(&stop,NULL);

 double comptime=(double)(stop.tv_sec-start.tv_sec)*1000+(double)(stop.tv_usec-start.tv_usec)/1000;

gettimeofday(&gpustop,NULL);
double gputime=(double)(gpustop.tv_sec-gpustart.tv_sec)*1000+(double)(gpustop.tv_usec-gpustart.tv_usec)/1000;
//H = H_gpu;
//printf("\n GPU computation time: %lf ms GPU end to end: %lf ms\n",comptime,gputime);
// printf("\n GPU printed: %lf in %lf ms",H,comptime);

count++;
return H;

}



}
