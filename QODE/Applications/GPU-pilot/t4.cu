#include <cuda.h>
#include <cuda_runtime.h>
#include "stdio.h"
#include <sys/time.h>
#define DEBUG 1

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

  __device__ int get_tid()
  {
    int blockId = blockIdx.x + blockIdx.y * gridDim.x
      + blockIdx.z * gridDim.x * gridDim.y;
    return  blockId * (blockDim.x * blockDim.y * blockDim.z)
      + (threadIdx.z * (blockDim.x * blockDim.y))
      + (threadIdx.y * blockDim.x)
      + threadIdx.x;
  }

  __global__ void dimer_1min1pls_t4(int n_elem, double *d_H, int *d_i, int *d_j, int *d_dim, float *d_sign, int *d_n_orb1, int *d_n_orb2, double *d_Rc1, double *d_Rcca1, double *d_Ra2, double *d_Rcaa2, double *d_h, double *d_V1, double *d_V2)
  {
    const int tid = get_tid();
    if(tid < n_elem) {
      const int n_orb1 = d_n_orb1[tid];
      const int n_orb2 = d_n_orb2[tid];
      double tmp = 0.0;
      int p,q,r,s;

      // In practice the n calculated here is not equal to the "actual" n here
      // since multiple kernel calls will be needed.
      // n_elem is the length

      // Top Loop
      // Use the tid to calculate variables in five dimensions
      //int n = tid;
      s = tid % n_orb2;
      r = ((tid - s) / (n_orb2)) % n_orb1;
      q = ((tid - r - s) * (n_orb2 * n_orb1)) % n_orb1;
      p = ((tid - q - r - s) / (n_orb2 * n_orb1 * n_orb1)) % n_orb1;
      tmp = d_V1[(((p * n_orb1 + q ) * n_orb1 + r) * n_orb2 + s) * n_elem + tid]
        * d_Rcca1[((q * n_orb1 + p)* n_orb1 + r) * n_elem + tid]
        * d_Ra2[s * n_elem + tid];

      // Middle Loop
      // 5d calculation
      s = tid % n_orb2;
      r = ((tid - s) / (n_orb2)) % n_orb2;
      q = ((tid - r - s) * (n_orb2 * n_orb2)) % n_orb2;
      p = ((tid - q - r - s) / (n_orb2 * n_orb2 * n_orb2)) % n_orb1;

      tmp += d_V2[(((p * n_orb2 + q ) * n_orb2 + r) * n_orb2 + s) * n_elem + tid]
        * d_Rcaa2[((q * n_orb2 + s) * n_orb2 + r) * n_elem + tid]
        * d_Rc1[p * n_elem + tid];
      tmp *= 2.0;

      // Bottom loop
      // 3d calculation
      s = tid % n_orb2;
      r = ((tid - s) / (n_orb2)) % n_orb1;
      tmp += d_h[(r * n_orb2 * s) * n_elem + tid]
        * d_Rc1[r * n_elem + tid]
        * d_Ra2[s * n_elem + tid];
      tmp *= d_sign[tid];

      // Assignment
      int index = d_i[tid] * d_dim[tid] + d_j[tid];
      d_H[index * n_elem + tid] = tmp;

      }
  }

void dimer_1min1pls_loop(int n_elem, double** H, int* i, int* j, int* dim, float* sign, int* n_orb1, int* n_orb2, double** Rc1,  double** Rcca1, double** Ra2, double** Rcaa2, double** h, double** V1, double** V2)
{
    struct timeval start,stop; 
    gettimeofday(&start,0);

    double *d_H, *d_Rc1, *d_Rcca1, *d_Ra2, *d_Rcaa2, *d_h, *d_V1, *d_V2;
    int *d_i, *d_j, *d_n_orb1, *d_n_orb2, *d_dim;
    float *d_sign;

    const int n1 = n_orb1[0];
    const int n2 = n_orb2[0];

    const int V1_len = n_elem * n2 * n1 * n1 * n1;
    const int V2_len = n_elem * n2 * n2 * n2 * n1;
    const int Rc1_len = n_elem * n1;
    const int Ra2_len = n_elem * n2;
    const int Rcca1_len = n_elem * n1 * n1 *n1;
    const int Rcaa2_len = n_elem * n2 * n2 * n2;
    const int h_len = n_elem * n1 * n2;
    const int H_len = n_elem * n1 * n1 * n1 * n2; // could be an issue if n1 != n2


    const int chunks = n1 * n2; // for maximum thread utilization this should go evenly into H_len
    const int elemPerChunk = ceil(n_elem / chunks);
    const int threadsPerChunk = ceil(H_len / chunks);
    const dim3 dimblock(4, 4, 4);
    const dim3 dimgrid(ceil(threadsPerChunk / dimblock.x), ceil(threadsPerChunk / dimblock.y), ceil(threadsPerChunk / dimblock.z));

    if(DEBUG) {
      printf("Chunks: %d\tThreads per chunk: %d\n",chunks,threadsPerChunk);
      printf("dimblock:\t<%d\t%d\t%d>\ndimgrid:\t<%d\t%d\t%d>\n",dimblock.x,dimblock.y,dimblock.z,dimgrid.x,dimgrid.y,dimgrid.z);
    }

    const int V1_size = sizeof(double) * ceil(V1_len / chunks);
    const int V2_size = sizeof(double) * ceil(V2_len / chunks);
    const int Rc1_size = sizeof(double) * ceil(Rc1_len / chunks);
    const int Ra2_size = sizeof(double) * ceil(Ra2_len / chunks);
    const int Rcca1_size = sizeof(double) * ceil(Rcca1_len / chunks);
    const int Rcaa2_size = sizeof(double) * ceil(Rcaa2_len / chunks);
    const int h_size = sizeof(double) * ceil(h_len / chunks);
    const int i_size = sizeof(int) * elemPerChunk; // for all the n-elem int arrays
    const int sign_size = sizeof(float) * elemPerChunk;
    const int H_size = sizeof(double) * threadsPerChunk;

    if(DEBUG) {
      printf("Starting device memory allocation\n");
    }

    gpuErr(cudaMalloc((void **) &d_H,       H_size));
    gpuErr(cudaMalloc((void **) &d_Rc1,     Rc1_size));
    gpuErr(cudaMalloc((void **) &d_V1,      V1_size));
    gpuErr(cudaMalloc((void **) &d_V2,      V2_size));
    gpuErr(cudaMalloc((void **) &d_Ra2,     Ra2_size));
    gpuErr(cudaMalloc((void **) &d_Rcca1,   Rcca1_size));
    gpuErr(cudaMalloc((void **) &d_Rcaa2,   Rcaa2_size));
    gpuErr(cudaMalloc((void **) &d_h,       h_size));
    gpuErr(cudaMalloc((void **) &d_i,       i_size));
    gpuErr(cudaMalloc((void **) &d_j,       i_size));
    gpuErr(cudaMalloc((void **) &d_dim,     i_size));
    gpuErr(cudaMalloc((void **) &d_n_orb1,  i_size));
    gpuErr(cudaMalloc((void **) &d_n_orb2,  i_size));
    gpuErr(cudaMalloc((void **) &d_sign,    sign_size));

    if(DEBUG) {
      printf("Finished device memory allocation\n");
    }

    for(int k=0; k<chunks; k++) {

      if(DEBUG) {
        printf("Chunk %d: starting host to device memcpys\n",k);
      }

      gpuErr(cudaMemcpy(d_Rc1, Rc1 + k * Rc1_len, Rc1_size, cudaMemcpyHostToDevice));
      gpuErr(cudaMemcpy(d_V1, V1 + k * V1_len, V1_size, cudaMemcpyHostToDevice));
      gpuErr(cudaMemcpy(d_V2, V2 + k * V2_len , V2_size, cudaMemcpyHostToDevice));
      gpuErr(cudaMemcpy(d_Ra2, Ra2 + k * Ra2_len, Ra2_size, cudaMemcpyHostToDevice));
      gpuErr(cudaMemcpy(d_Rcca1, Rcca1 + k * Rcca1_len, Rcca1_size, cudaMemcpyHostToDevice));
      gpuErr(cudaMemcpy(d_Rcaa2, Rcaa2 + k * Rcaa2_len, Rcaa2_size, cudaMemcpyHostToDevice));
      gpuErr(cudaMemcpy(d_h, h + k * h_len, h_size, cudaMemcpyHostToDevice));
      gpuErr(cudaMemcpy(d_i, i + k * elemPerChunk, i_size, cudaMemcpyHostToDevice));
      gpuErr(cudaMemcpy(d_j, j + k * elemPerChunk, i_size, cudaMemcpyHostToDevice));
      gpuErr(cudaMemcpy(d_n_orb1, n_orb1 + k * elemPerChunk, i_size, cudaMemcpyHostToDevice));
      gpuErr(cudaMemcpy(d_n_orb2, n_orb2 + k * elemPerChunk, i_size, cudaMemcpyHostToDevice));
      gpuErr(cudaMemcpy(d_sign, sign + k * elemPerChunk, sign_size, cudaMemcpyHostToDevice));

      if(DEBUG) {
        printf("Chunk %d: finished host to device memcpys\n",k);
        printf("Chunk %d: launching dimer_1min1pls_t4 kernel\n",k);
      }

      dimer_1min1pls_t4<<<dimblock,dimgrid>>>(n_elem,d_H,d_i,d_j,d_dim,d_sign,d_n_orb1,d_n_orb2,d_Rc1,d_Rcca1,d_Ra2,d_Rcaa2,d_h,d_V1,d_V2);
      gpuErr(cudaPeekAtLastError());
      gpuErr(cudaDeviceSynchronize());

      memcpy(H+k*H_len, d_H, H_size);


      if(DEBUG) {
        printf("Chunk %d: finished executing dimer_1min1pls_t4 kernel \n",k);
      }


    }

    gpuErr(cudaFree(d_Rc1));
    gpuErr(cudaFree(d_V1));
    gpuErr(cudaFree(d_V2));
    gpuErr(cudaFree(d_Ra2));
    gpuErr(cudaFree(d_Rcca1));
    gpuErr(cudaFree(d_Rcaa2));
    gpuErr(cudaFree(d_h));
    gpuErr(cudaFree(d_i));
    gpuErr(cudaFree(d_j));
    gpuErr(cudaFree(d_n_orb1));
    gpuErr(cudaFree(d_n_orb2));
    gpuErr(cudaFree(d_sign))
    gettimeofday(&stop,0);
    if(DEBUG) {
      double t = (double)(stop.tv_sec-start.tv_sec)*1000+(double)(stop.tv_usec-start.tv_usec)/1000;
      printf("dimer_1min1pls_loop_t4 finished in %f ms\n", t);
    }


}
}
