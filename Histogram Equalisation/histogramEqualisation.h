#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <cuda.h>



//histrogram


__global__
void hist1(float *image,float *hist, int width, int height)
{
        int tidx,tidy;
        tidx=blockIdx.x*blockDim.x + threadIdx.x;
        tidy=blockIdx.y*blockDim.y + threadIdx.y;

	int index=tidx*width+tidy;       

if (tidx<height && tidy<width)
{
{
                atomicAdd(&hist[(int)image[index]],1);        
               
 }      
}
}

// normalization 


__global__
void hist_f(float *image,float *cdf, int width, int height)

{
int tidx,tidy;
        tidx=blockIdx.x*blockDim.x + threadIdx.x;
        tidy=blockIdx.y*blockDim.y + threadIdx.y;
	int index=tidy+tidx*width;
     
              if (tidx<height && tidy<width)
              {
		{
       			 image[index]=(255*cdf[(int)image[index]]);
                
		}
	       }

	      
}





 
