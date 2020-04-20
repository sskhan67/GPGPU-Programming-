#include "PyC_types.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include "stdio.h"
#include <sys/time.h>
#define DEBUG 1
#define T 6
extern "C" {

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
    {	
        for (int s2=0;  s2<n_orb2;  s2++)
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

void dimer_1min1pls_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyFloat* sign, PyInt* n_orb1_in, PyInt* n_orb2_in, Double** Rc1,  Double** Rcca1, Double** Ra2, Double** Rcaa2, Double** h, Double** V1, Double** V2)
{
    struct timeval start,stop; 
    gettimeofday(&start,0);
    double *d_V1112,*d_Rcca1,*d_Ra2,*d_V1222,*d_Rcaa2,*d_Rc1,*d_H,*h_H,*d_h,*h_Hr,*d_Hr;

    // All of these assume that the values of n_orb1 and n_orb2 don't change
    const int n_orb1 =  n_orb1_in[0];
    const int n_orb2 =  n_orb2_in[0];
    const int blocks =  (n_orb1*n_orb1*n_orb1)/(T*T);
    const int nbpgrid=  n_orb1/T;
    const int N_V1112 = n_orb1*n_orb1*n_orb1*n_orb2;
    const int N_Rcca1 = n_orb1*n_orb1*n_orb1;
    const int N_Ra2 =   n_orb2;
    const int N_V1222 = n_orb1*n_orb2*n_orb2*n_orb2;
    const int N_Rcaa2 = n_orb2*n_orb2*n_orb2;
    const int N_Rc1 =   n_orb1;
    const int N_H =     n_orb1*n_orb1*n_orb1; // this assumes n_orb1 = n_orb2
    const int N_h =     n_orb1*n_orb2;
    const int N_Hr =    blocks;

    const dim3 dimblock(T,T,T);
    const dim3 dimgrid(nbpgrid,nbpgrid,nbpgrid);
    const dim3 dimblockR(T*T);
    const dim3 dimgridR(blocks);

    gpuErr(cudaMalloc((void **) &d_V1112,   sizeof(double)*N_V1112));
    gpuErr(cudaMalloc((void **) &d_Ra2,     sizeof(double)*N_Ra2));
    gpuErr(cudaMalloc((void **) &d_V1222,   sizeof(double)*N_V1222));
    gpuErr(cudaMalloc((void **) &d_Rcca1,   sizeof(double)*N_Rcca1));
    gpuErr(cudaMalloc((void **) &d_Rcaa2,   sizeof(double)*N_Rcaa2));
    gpuErr(cudaMalloc((void **) &d_Rc1,     sizeof(double)*N_Rc1));
    gpuErr(cudaMalloc((void **) &d_H,       sizeof(double)*N_H));
    gpuErr(cudaMalloc((void **) &d_h,       sizeof(double)*N_h));
    gpuErr(cudaMalloc((void **) &d_Hr,      sizeof(double)*N_Hr));
    h_Hr=(double *)malloc(sizeof(double)*N_Hr);
    h_H=(double *)malloc(sizeof(double)*N_H);

    for(int n=0; n<n_elem; n++) {
        if(DEBUG) printf("n = %d\n", n);
        if(DEBUG) printf("n_orb1: %d, n_orb2: %d\n", n_orb1_in[n], n_orb2_in[n]);
        int index = i[n]*dim[n]+j[n];
        double tmp = 0.0;
        gpuErr(cudaMemcpy(d_Rcca1,  Rcca1[n],   sizeof(double)*N_Rcca1, cudaMemcpyHostToDevice));
        if(DEBUG) printf("Rcca1 copied\n");
        gpuErr(cudaMemcpy(d_Rcaa2,  Rcaa2[n],   sizeof(double)*N_Rcaa2, cudaMemcpyHostToDevice));
        if(DEBUG) printf("Rcaa2 copied\n");
        gpuErr(cudaMemcpy(d_Rc1,    Rc1[n],     sizeof(double)*N_Rc1,   cudaMemcpyHostToDevice));
        if(DEBUG) printf("Rc1 copied\n");
        gpuErr(cudaMemcpy(d_Ra2,    Ra2[n],     sizeof(double)*N_Ra2,   cudaMemcpyHostToDevice));
        if(DEBUG) printf("Ra2 copied\n");
        gpuErr(cudaMemcpy(d_h,      h[n],       sizeof(double)*N_h,     cudaMemcpyHostToDevice));
        if(DEBUG) printf("h copied\n");
        gpuErr(cudaMemcpy(d_V1112,  V1[n],      sizeof(double)*N_V1112, cudaMemcpyHostToDevice));
        if(DEBUG) printf("V1112 copied\n");
        gpuErr(cudaMemcpy(d_V1222,  V2[n],      sizeof(double)*N_V1222, cudaMemcpyHostToDevice));
        if(DEBUG) printf("V1222 copied\n");

        if(DEBUG) printf("Launching outerloop\n");
        outerloop<<<dimgrid,dimblock>>>(d_V1112,d_Rcca1,d_Ra2,d_H,n_orb1,n_orb2,d_V1222,d_Rcaa2,d_Rc1,d_h);
        gpuErr(cudaPeekAtLastError());
        gpuErr(cudaMemcpy(h_H, d_H, sizeof(double)*N_H, cudaMemcpyDeviceToHost));
        gpuErr(cudaDeviceSynchronize());
        if(DEBUG) printf("Launching reduction\n");
        reduction<<<dimgridR,dimblockR,sizeof(double)*T*T>>>(d_H,d_Hr);	
        gpuErr(cudaPeekAtLastError());
        gpuErr(cudaMemcpy(h_Hr, d_Hr, sizeof(double)*N_Hr, cudaMemcpyDeviceToHost));
        gpuErr(cudaDeviceSynchronize());
        for(int k=0; k<blocks; k++) {
            tmp += h_Hr[k];
        }
        H[n][index] = tmp;
    }
    gpuErr(cudaFree(d_V1112));
    gpuErr(cudaFree(d_Rcca1));
    gpuErr(cudaFree(d_Ra2));
    gpuErr(cudaFree(d_V1222));
    gpuErr(cudaFree(d_Rcaa2));
    gpuErr(cudaFree(d_Rc1));
    gpuErr(cudaFree(d_H));
    gpuErr(cudaFree(d_h));
    gpuErr(cudaFree(d_Hr));
    free(h_Hr);
    free(h_H);

    gettimeofday(&stop,0);
    if(DEBUG) {
        double t = (double)(stop.tv_sec-start.tv_sec)*1000+(double)(stop.tv_usec-start.tv_usec)/1000;
        printf("dimer_1min1pls_loop inline cpu version finished in %f ms\n", t);
    }
}

}