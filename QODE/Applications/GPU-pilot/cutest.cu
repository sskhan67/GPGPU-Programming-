#include "PyC_types.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include "stdio.h"
#include <sys/time.h>
#include <iostream>

#define DEBUG 1
#define T 6

extern "C" {


void dimer_1min1pls_loop(PyInt n_elem, Double** H, PyInt* i, PyInt* j, PyInt* dim, PyFloat* sign, PyInt* n_orb1, PyInt* n_orb2, Double** Rc1,  Double** Rcca1, Double** Ra2, Double** Rcaa2, Double** h, Double** V1, Double** V2)
{
    struct timeval start,stop; 
    gettimeofday(&start,0);
    for(int n=0; n<n_elem; n++) {
        int index = i[n]*dim[n]+j[n];

        // Upper loop
        for(int p1=0; p1<n_orb1[n]; p1++) {
            for(int r1=0; r1<n_orb1[n]; r1++) {
                for(int q1=0; q1<n_orb1[n]; q1++) {
                    for(int s2=0; s2<n_orb2[n]; s2++) {
                        H[n][index] += sign[n] * V1[n][((p1*n_orb1[n] + q1)*n_orb1[n] + r1)*n_orb2[n] + s2] * Rcca1[n][(q1*n_orb1[n] + p1)*n_orb1[n] + r1] * Ra2[n][s2];
                    }
                }
            }
        }

        // Middle loop
        for(int p11=0; p11<n_orb1[n]; p11++) {
            for(int q2=0; q2<n_orb2[n]; q2++) {
                for(int r2=0; r2<n_orb2[n]; r2++) {
                    for(int s2=0; s2<n_orb2[n]; s2++) {
                        H[n][index] += sign[n] * V2[n][((p11*n_orb2[n] + q2)*n_orb2[n] + r2)*n_orb2[n] + s2] * Rc1[n][p1] * Rcaa2[n][(q2*n_orb2[n] + s2)*n_orb2[n] + r2];
                    }
                }
            }
        }

        H[n][index] *= 2;

        // Bottom Loop
        for(int p12=0; p12<n_orb1[n]; p12++) {
            for(int q22=0; q22<n_orb2[n]; q22++) {
                H[n][index] += h[n][p12*n_orb2[n]+q22] * Rc1[n][p12] * Ra2[q22];
            }
        }
    }
    gettimeofday(&stop,0);
    if(DEBUG) {
        double t = (double)(stop.tv_sec-start.tv_sec)*1000+(double)(stop.tv_usec-start.tv_usec)/1000;
        printf("dimer_1min1pls_loop finished in %f ms\n", t);
    }
}


}

 /*
    //std::vector<std::thread> hostThreads;
    // lambda expression to create streams
    auto createStreams = [streams]() {
        for(int i=0; i<numStreams; i++) {
            gpuErr(cudaStreamCreate(&streams[i]));
        }
    };
   //hostThreads.push_back(std::move(std::thread(createStreams))); 
   for(std::thread& t : hostThreads) {
        if(t.joinable())
        t.join();  
    }
    */


    // not sure if cudaMalloc is threadsafe if called from multiple host threads
    // will experiment with this -- concurrent allocation would save significant time
    // generalized lambda expression to perform a cudamalloc
    /*
    auto cudaPreMalloc = [](double *arr, unsigned long size) {
        gpuErr(cudaMalloc((void **) &arr, size));
    };
    threads.push(std::thread(cudaPreMalloc, &d_V1112, sizeof(double)*N4*n_elem/numChunks));
    threads.push(std::thread(cudaPreMalloc, &d_Rcca1, sizeof(double)*N3*n_elem/numChunks));
    threads.push(std::thread(cudaPreMalloc, &d_Ra2,   sizeof(double)*n_orb2*n_elem/numChunks));
    threads.push(std::thread(cudaPreMalloc, &d_V1222, sizeof(double)*n_orb1*n_orb2*n_orb2*n_orb2*n_elem/numChunks));
    threads.push(std::thread(cudaPreMalloc, &d_Rcaa2, sizeof(double)*n_orb2*n_orb2*n_orb2*n_elem/numChunks));
    threads.push(std::thread(cudaPreMalloc, &d_Rc1,   sizeof(double)*n_orb1*n_elem/numChunks));
    threads.push(std::thread(cudaPreMalloc, &d_H,     sizeof(double)*n_orb1*n_orb1*n_orb1*n_elem/numChunks));
    threads.push(std::thread(cudaPreMalloc, &d_h,     sizeof(double)*n_orb1*n_orb2*n_elem/numChunks));
    threads.push(std::thread(cudaPreMalloc, &d_Hr,    sizeof(double)*blocks*n_elem/numChunks));
    */

