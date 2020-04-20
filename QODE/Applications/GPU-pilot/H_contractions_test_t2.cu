/*   (C) Copyright 2018, 2020 Anthony D. Dutoi and Yuhong Liu
 *
 *   This file is part of Qode.
 *
 *   Qode is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Qode is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Qode.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "PyC_types.h"
#include<cuda.h>
//#include"wrapper.h"
#include "cuda_runtime.h"
#include "stdlib.h"
#include "time.h"
#include "stdio.h"
//#include "wrapper.h"
#include "PyC_types.h"
#include "stdlib.h"
#include <sys/time.h>
#include<omp.h>
#define T 6
#define TR 8
extern "C"
 {
/*
// Global variables
double *d_V1112,*d_Rcca1,*d_Ra2,*d_V1222,*d_Rcaa2,*d_Rc1;
double *h_H,*d_h;
double *d_H,*h_Hr,*d_Hr;
*/




// reduction works

__device__ void warpReduce(volatile double* sdata, int tid) {
sdata[tid] += sdata[tid + 32];
sdata[tid] += sdata[tid + 16];
sdata[tid] += sdata[tid + 8];
sdata[tid] += sdata[tid + 4];
sdata[tid] += sdata[tid + 2];
sdata[tid] += sdata[tid + 1];
}

__global__
void reduction(double *H,double *R, int numelements)
{
extern __shared__ volatile double sdata[];
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(blockDim.x*2) + threadIdx.x;
if(i<numelements)
{
    sdata[tid] = H[i] + H[i+blockDim.x];
    
    __syncthreads();
    // do reduction in shared mem
    for (unsigned int s=blockDim.x/2; s>32; s>>=1) {
        //if (tid < s)
	{
	 sdata[tid] += sdata[tid + s];
        __syncthreads();
	}
	    }

	if (tid<32) warpReduce(sdata, tid);
	

  // write result for this block to global mem
    if (tid == 0)
	{ 
	R[blockIdx.x] = sdata[0];
	
	}
}
}










/*
__global__
void reduction(double *H,double *R, int numelements)
{

extern __shared__  double sdata[];
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;

 if (i<numelements)
        {
        sdata[tid] = H[i];

    __syncthreads();
        //if (i<numelements)
//      {
        for (unsigned int s=blockDim.x/2; s>0; s>>=1)
                {
                if (tid<s)
                {
                        sdata[tid]+=sdata[tid+s];
                }
                __syncthreads();
                }
        if (tid==0)
        {
        R[blockIdx.x]=sdata[0];
        }
        }
}


*/




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

//double dimer_1min1pls(int n_orb1, int n_orb2, double* Rc1, double* Rcca1,  double* Ra2,  double* Rcaa2,  double* h,  double* V1112,  double* V1222,int freevaribales)
void dimer_1min1pls_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyFloat* sign, PyInt* n_orb1_in, PyInt* n_orb2_in, Double** Rc1,  Double** Rcca1, Double** Ra2, Double** Rcaa2, Double** h, Double** V1, Double** V2)

{	
	//printf("hello");

        //cudaMemcpy(d_h,      h[n],       sizeof(double)*N_h,     cudaMemcpyHostToDevice);
        //if(DEBUG) printf("h copied\n");
	double *d_V1112,*d_Rcca1,*d_Ra2,*d_V1222,*d_Rcaa2,*d_Rc1,*d_H,*h_H,*d_h,*h_Hr,*d_Hr;

    // All of these assume that the values of n_orb1 and n_orb2 don't change
    const int n_orb1 =  n_orb1_in[0];
    const int n_orb2 =  n_orb2_in[0];
    const int blocks =  (n_orb1*n_orb1*n_orb1)/(TR*TR);
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

    cudaMalloc((void **) &d_V1112,   sizeof(double)*N_V1112);
    cudaMalloc((void **) &d_Ra2,     sizeof(double)*N_Ra2);
    cudaMalloc((void **) &d_V1222,   sizeof(double)*N_V1222);
    cudaMalloc((void **) &d_Rcca1,   sizeof(double)*N_Rcca1);
    cudaMalloc((void **) &d_Rcaa2,   sizeof(double)*N_Rcaa2);
    cudaMalloc((void **) &d_Rc1,     sizeof(double)*N_Rc1);
    cudaMalloc((void **) &d_H,       sizeof(double)*N_H);
    cudaMalloc((void **) &d_h,       sizeof(double)*N_h);
    cudaMalloc((void **) &d_Hr,      sizeof(double)*N_Hr);
    h_Hr=(double *)malloc(sizeof(double)*N_Hr);
    h_H=(double *)malloc(sizeof(double)*N_H);
    #pragma unroll 
  /*  
	cudaMemcpy(d_Rcca1,  Rcca1,   sizeof(double)*N_Rcca1, cudaMemcpyHostToDevice);
        cudaMemcpy(d_Rcaa2,  Rcaa2,   sizeof(double)*N_Rcaa2, cudaMemcpyHostToDevice);
        cudaMemcpy(d_Rc1,    Rc1,     sizeof(double)*N_Rc1,   cudaMemcpyHostToDevice);
        cudaMemcpy(d_Ra2,    Ra2,     sizeof(double)*N_Ra2,   cudaMemcpyHostToDevice);
        cudaMemcpy(d_h,      h,       sizeof(double)*N_h,     cudaMemcpyHostToDevice);
        cudaMemcpy(d_V1112,  V1,      sizeof(double)*N_V1112, cudaMemcpyHostToDevice);
        cudaMemcpy(d_V1222,  V2,      sizeof(double)*N_V1222, cudaMemcpyHostToDevice);
*/
    for(int n=0; n<n_elem; n++) {
        int index = i[n]*dim[n]+j[n];
        double tmp = 0.0;
        
	cudaMemcpy(d_Rcca1,  Rcca1[n],   sizeof(double)*N_Rcca1, cudaMemcpyHostToDevice);
        cudaMemcpy(d_Rcaa2,  Rcaa2[n],   sizeof(double)*N_Rcaa2, cudaMemcpyHostToDevice);
        cudaMemcpy(d_Rc1,    Rc1[n],     sizeof(double)*N_Rc1,   cudaMemcpyHostToDevice);
        cudaMemcpy(d_Ra2,    Ra2[n],     sizeof(double)*N_Ra2,   cudaMemcpyHostToDevice);
        cudaMemcpy(d_h,      h[n],       sizeof(double)*N_h,     cudaMemcpyHostToDevice);
        cudaMemcpy(d_V1112,  V1[n],      sizeof(double)*N_V1112, cudaMemcpyHostToDevice);
        cudaMemcpy(d_V1222,  V2[n],      sizeof(double)*N_V1222, cudaMemcpyHostToDevice);
       

	outerloop<<<dimgrid,dimblock>>>(d_V1112,d_Rcca1,d_Ra2,d_H,n_orb1,n_orb2,d_V1222,d_Rcaa2,d_Rc1,d_h);
        
	cudaMemcpy(h_H, d_H, sizeof(double)*N_H, cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
        //f(DEBUG) printf("Launching reduction\n");
        reduction<<<dimgridR,dimblockR,sizeof(double)*TR*TR>>>(d_H,d_Hr,(n_orb1*n_orb1*n_orb1));	
        //gpuErr(cudaPeekAtLastError());
        cudaMemcpy(h_Hr, d_Hr, sizeof(double)*N_Hr, cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
        for(int k=0; k<blocks; k++) {
            tmp += h_Hr[k];
        }
        H[n][index] = tmp*sign[n];
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

//printf("hello ");



}









Double monomer(PyInt n_orb, Double* Rca, Double* Rccaa, Double* h, Double* V)
	{
	Double H = 0;
	//printf("%d N orbibtal: ",n_orb);
	PyInt p=0;
	for (p=0;  p<n_orb;  p++)
		{
		PyInt q=0;
		for (q=0;  q<n_orb;  q++)
			{
			PyInt r=0;
			for ( r=0;  r<n_orb;  r++)
				{
				PyInt s=0;
				for (s=0;  s<n_orb;  s++)
					{
					H += V[((p*n_orb + q)*n_orb + r)*n_orb + s] * Rccaa[((p*n_orb + q)*n_orb + s)*n_orb + r];
					}
				}
			}
		}
	PyInt p1=0;
	for (p1=0;  p1<n_orb;  p1++)
		{
		PyInt q1=0;
		for (q1=0;  q1<n_orb;  q1++)
			{
			H += h[p1*n_orb + q1] * Rca[p1*n_orb + q1];
			}
		}
	return H;
	}



Double dimer_2min2pls(PyInt n_orb1, PyInt n_orb2, Double* Rcc1, Double* Raa2, Double* V)
	{
	Double H = 0;
	PyInt p1=0;
	for (p1=0;  p1<n_orb1;  p1++)
		{
		PyInt q1=0;
		for (q1=0;  q1<n_orb1;  q1++)
			{
			PyInt r2=0;
			for (r2=0;  r2<n_orb2;  r2++)
				{
				PyInt s2=0;
				for ( s2=0;  s2<n_orb2;  s2++)
					{
					H += V[((p1*n_orb1 + q1)*n_orb2 + r2)*n_orb2 + s2] * Rcc1[p1*n_orb1 + q1] * Raa2[s2*n_orb2 + r2];
					}
				}
			}
		}
	return H;
 
	}
/*
static int count=0;
Double dimer_1min1pls(PyInt n_orb1, PyInt n_orb2, Double* Rc1, Double* Rcca1, Double* Ra2, Double* Rcaa2, Double* h, Double* V1112, Double* V1222, int freevariables)
	{
	test_wrapper(n_orb1, n_orb2,  Rc1, Rcca1,  Ra2, Rcaa2,  h,  V1112,  V1222,freevariables);
	count++;
//if(count>10)
	//exit(0);

	}*/






/*

Double dimer_1min1pls(PyInt n_orb1, PyInt n_orb2, Double* Rc1, Double* Rcca1, Double* Ra2, Double* Rcaa2, Double* h, Double* V1112, Double* V1222)
	{
	Double H = 0;
	//printf("%d ord1: ",n_orb1);
	//printf(" n_orb2:%d\n ",n_orb2);
	//printf(" Rc1:%f \n",*Rc1);
	//printf(" Rcca1:%f\n ",*Rcca1);
	//printf(" Ra2: %f\n ",*Ra2);
	//printf(" h:%f",*h);
	//printf(" V1112:%f",* V1112);
	PyInt p1=0;
	for ( p1=0;  p1<n_orb1;  p1++)
		{
		PyInt q1=0;
		for ( q1=0;  q1<n_orb1;  q1++)
			{
			PyInt r1=0;
			for ( r1=0;  r1<n_orb1;  r1++)
				{
				PyInt s2=0;
				for ( s2=0;  s2<n_orb2;  s2++)
					{
					H += V1112[((p1*n_orb1 + q1)*n_orb1 + r1)*n_orb2 + s2] * Rcca1[(q1*n_orb1 + p1)*n_orb1 + r1] * Ra2[s2];
					}
				}
			}
		}
	PyInt p11=0;
	for (p11=0;  p11<n_orb1;  p11++)
		{
		PyInt q2=0;
		for (q2=0;  q2<n_orb2;  q2++)
			{
			PyInt r2=0;
			for (r2=0;  r2<n_orb2;  r2++)
				{
				PyInt s2=0;
				for ( s2=0;  s2<n_orb2;  s2++)
					{
					H += V1222[((p11*n_orb2 + q2)*n_orb2 + r2)*n_orb2 + s2] * Rc1[p1] * Rcaa2[(q2*n_orb2 + s2)*n_orb2 + r2];
					}
				}
			}
		}
	H *= 2;
	PyInt p12=0;
	for (p12=0;  p12<n_orb1;  p12++)
		{
		PyInt q22=0;
		for (q22=0;  q22<n_orb2;  q22++)
			{
			H += h[p12*n_orb2 + q22] * Rc1[p12] * Ra2[q22];
			}
		}
	//printf("%f H: ",H);
	return H;
	}

*/
Double dimer_00(PyInt n_orb1, PyInt n_orb2, Double* Rca1, Double* Rca2, Double* V)
	{
	Double H = 0;
	PyInt p1=0;
	for ( p1=0;  p1<n_orb1;  p1++)
		{
		PyInt q2=0;
		for ( q2=0;  q2<n_orb2;  q2++)
			{
			PyInt r1=0;
			for ( r1=0;  r1<n_orb1;  r1++)
				{
				PyInt s2=0;
				for ( s2=0;  s2<n_orb2;  s2++)
					{
					H += V[((p1*n_orb2 + q2)*n_orb1 + r1)*n_orb2 + s2] * Rca1[p1*n_orb1 + r1] * Rca2[q2*n_orb2 + s2];
					}
				}
			}
		}
	return 4*H;
	}



Double dimer_1pls1min(PyInt n_orb1, PyInt n_orb2, Double* Ra1, Double* Rcaa1, Double* Rc2, Double* Rcca2, Double* h, Double* V2221, Double* V2111)
	{
	Double H = 0;
	PyInt p2=0;
	for (p2=0;  p2<n_orb2;  p2++)
		{
		PyInt q2=0;
		for ( q2=0;  q2<n_orb2;  q2++)
			{
			PyInt r2=0;
			for ( r2=0;  r2<n_orb2;  r2++)
				{
				PyInt s1=0;
				for ( s1=0;  s1<n_orb1;  s1++)
					{
					H += V2221[((p2*n_orb2 + q2)*n_orb2 + r2)*n_orb1 + s1] * Rcca2[(q2*n_orb2 + p2)*n_orb2 + r2] * Ra1[s1];
					}
				}
			}
		}
	PyInt p22=0;
	for (p22=0;  p22<n_orb2;  p22++)
		{
		PyInt q1=0;
		for ( q1=0;  q1<n_orb1;  q1++)
			{
			PyInt r1=0;
			for (r1=0;  r1<n_orb1;  r1++)
				{
				PyInt s2=0;
				for (s2=0;  s2<n_orb1;  s2++)
					{
					H += V2111[((p22*n_orb1 + q1)*n_orb1 + r1)*n_orb1 + s2] * Rc2[p22] * Rcaa1[(q1*n_orb1 + s2)*n_orb1 + r1];
					}
				}
			}
		}
	H *= 2;
	PyInt p3=0;
	for (p3=0;  p3<n_orb2;  p3++)
		{
		PyInt q3=0;
		for (q3=0;  q3<n_orb1;  q3++)
			{
			H += h[p3*n_orb1 + q3] * Rc2[p3] * Ra1[q3];
			}
		}
	return H;
	}



Double dimer_2pls2min(PyInt n_orb1, PyInt n_orb2, Double* Raa1, Double* Rcc2, Double* V)
	{
	Double H = 0;
	PyInt p2=0;
	for (p2=0;  p2<n_orb2;  p2++)
		{
		PyInt q2=0;
		for (q2=0;  q2<n_orb2;  q2++)
			{
			PyInt r1=0;
			for (r1=0;  r1<n_orb1;  r1++)
				{
				PyInt s1=0;
				for (s1=0;  s1<n_orb1;  s1++)
					{
					H += V[((p2*n_orb2 + q2)*n_orb1 + r1)*n_orb1 + s1] * Rcc2[p2*n_orb2 + q2] * Raa1[s1*n_orb1 + r1];
					}
				}
			}
		}
	return H;
	}

/*

void dimer_1min1pls_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyFloat* sign, PyInt* n_orb1, PyInt* n_orb2, Double** Rc1,  Double** Rcca1, Double** Ra2, Double** Rcaa2, Double** h, Double** V1, Double** V2)
{
    int n;
    for(n=0; n<n_elem; n++) {
        int index = i[n]*dim[n]+j[n];
        double tmp = 0.0;

        int p1,r1,q1,s2;
        for(p1=0; p1<n_orb1[n]; p1++) {
            for(r1=0; r1<n_orb1[n]; r1++) {
                for(q1=0; q1<n_orb1[n]; q1++) {
                    for( s2=0; s2<n_orb2[n]; s2++) {
                        tmp += V1[n][((p1*n_orb1[n] + q1)*n_orb1[n] + r1)*n_orb2[n] + s2] * Rcca1[n][(q1*n_orb1[n] + p1)*n_orb1[n] + r1] * Ra2[n][s2];
                    }
                }
            }
        }

     	int p11,q2,r2;
        for(p11=0; p11<n_orb1[n]; p11++) {
            for(q2=0; q2<n_orb2[n]; q2++) {
                for(r2=0; r2<n_orb2[n]; r2++) {
                    for(s2=0; s2<n_orb2[n]; s2++) {
                        tmp += V2[n][((p11*n_orb2[n] + q2)*n_orb2[n] + r2)*n_orb2[n] + s2] * Rc1[n][p11] * Rcaa2[n][(q2*n_orb2[n] + s2)*n_orb2[n] + r2];
                    }
                }
            }
        }

        tmp *= 2;

        int p12,q22;
        for(p12=0; p12<n_orb1[n]; p12++) {
            for(q22=0; q22<n_orb2[n]; q22++) {
                tmp += h[n][p12*n_orb2[n]+q22] * Rc1[n][p12] * Ra2[n][q22];
            }
        }
        H[n][index] = tmp*sign[n];
    }
    
    
}

*/

void monomer_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyInt* n_orb, Double** Rca,  Double** Rccaa, Double** h, Double** V)
	{
	PyInt n=0;
	for ( n=0;  n<n_elem;  n++)
		{
		PyInt index = i[n]*dim[n] + j[n];
		H[n][index] = monomer(n_orb[n], Rca[n], Rccaa[n], h[n], V[n]);
		}
	return;
	}

void dimer_2min2pls_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyInt* n_orb1, PyInt* n_orb2, Double** Rcc1, Double** Raa2, Double** V)
	{
	PyInt n=0;
	for ( n=0;  n<n_elem;  n++)
		{
		PyInt index = i[n]*dim[n] + j[n];
		H[n][index] = dimer_2min2pls(n_orb1[n], n_orb2[n], Rcc1[n], Raa2[n], V[n]);
		}
	return;
	}


/*

void dimer_1min1pls_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyFloat* sign, PyInt* n_orb1, PyInt* n_orb2, Double** Rc1,  Double** Rcca1, Double** Ra2, Double** Rcaa2, Double** h, Double** V1, Double** V2)
	{

	int freevariables=0;
	PyInt n=0;
	#pragma unroll
	for ( n=0;  n<n_elem;  n++)
		{
		PyInt index = i[n]*dim[n] + j[n];
		H[n][index] = sign[n] * dimer_1min1pls(n_orb1[n], n_orb2[n], Rc1[n], Rcca1[n], Ra2[n], Rcaa2[n], h[n], V1[n], V2[n],freevariables);
		}
 	freevariables=1;
	dimer_1min1pls(n_orb1[n], n_orb2[n], Rc1[n], Rcca1[n], Ra2[n], Rcaa2[n], h[n], V1[n], V2[n],freevariables);
	return;
	}
*/



/*
void dimer_1min1pls_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyFloat* sign, PyInt* n_orb1, PyInt* n_orb2, Double** Rc1,  Double** Rcca1, Double** Ra2, Double** Rcaa2, Double** h, Double** V1, Double** V2)
        {
	 test_wrapper( n_elem,H,i, j,dim,sign,n_orb1,n_orb2,Rc1,Rcca1,Ra2,Rcaa2, h,V1,V2);
	//return 
	}

*/
void dimer_00_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyInt* n_orb1, PyInt* n_orb2, Double** Rca1, Double** Rca2, Double** V)
	{
	PyInt n=0;
	for ( n=0;  n<n_elem;  n++)
		{
		PyInt index = i[n]*dim[n] + j[n];
		H[n][index] = dimer_00(n_orb1[n], n_orb2[n], Rca1[n], Rca2[n], V[n]);
		}
	return;
	}

void dimer_1pls1min_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyFloat* sign, PyInt* n_orb1, PyInt* n_orb2, Double** Ra1,  Double** Rcaa1, Double** Rc2, Double** Rcca2, Double** h, Double** V1, Double** V2)
	{
	PyInt n=0;
	for ( n=0;  n<n_elem;  n++)
		{
		PyInt index = i[n]*dim[n] + j[n];
		H[n][index] = sign[n] * dimer_1pls1min(n_orb1[n], n_orb2[n], Ra1[n], Rcaa1[n], Rc2[n], Rcca2[n], h[n], V1[n], V2[n]);
		}
	return;
	}

void dimer_2pls2min_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyInt* n_orb1, PyInt* n_orb2, Double** Raa1, Double** Rcc2, Double** V)
	{
	PyInt n=0;
	for ( n=0;  n<n_elem;  n++)
		{
		PyInt index = i[n]*dim[n] + j[n];
		H[n][index] = dimer_2pls2min(n_orb1[n], n_orb2[n], Raa1[n], Rcc2[n], V[n]);
		}
	return;
	}


}
