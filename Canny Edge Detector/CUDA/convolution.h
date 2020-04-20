#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <cuda.h>

// convoultion function : core of the program, also main to source to parrallel programming, Inout: kernel, temp image, Output: Convoluated image 

__global__
void convoultion(float *image, float *output, float *g_kernel,float *d_flop1,float *d_global1, int height, int width, int k_height, int k_width)

{
	int tidx,tidy,offseti,offsetj;
	tidx=blockIdx.x*blockDim.x + threadIdx.x;
	tidy=blockIdx.y*blockDim.y + threadIdx.y;
	
	float kerw=(k_width>k_height)?k_width:k_height;

	float flop1=0, global1=0;


	// share memory 
	int tx=threadIdx.x;
	int ty=threadIdx.y;
        
	
extern __shared__ float Ashared[];
Ashared[tx*blockDim.y+ty]=image[tidx*width+tidy];
__syncthreads();
	

	if (tidx<height && tidy<width)
		{

		float sum=0;
		for (int m=0;m<kerw;m++)
		{
                         offseti = k_height>1?(-1*(k_height/2)+m):0;
                         offsetj =k_width>1? (-1*(k_width/2)+m):0;
                        if( (tx+offseti)>=0 && (tx+offseti)<blockDim.x && (ty+offsetj)>=0 && (ty+offsetj)<blockDim.y)
			{


                                sum+= Ashared[(tx+offseti)*blockDim.y+(ty+offsetj)]*g_kernel[m];
				flop1=flop1+2; // *,+ =2 floating point operation
				global1++; // g_kernel is in global memory

			}
			else if ((tidx+offseti>=0) && (tidx+offseti<height) && (tidy+offsetj)>=0 && (tidy+offsetj)<width)  
			{
			        
				sum+= image[(tidx+offseti)*width+(tidy+offsetj)]*g_kernel[m];
				
				flop1=flop1+2; // *,+ =2 floating point operation
                                global1++; // g_kernel is in global memory

			}
		
		}
                              output[(tidx*width)+tidy]=sum;
 	
			      global1++; // 1 global memory access
			      d_global1[tidx*width+tidy]=global1;
	                      d_flop1[tidx*width+tidy]=flop1;

		}






/*
	// 2. Convolution using global memory
   if (tidx<height && tidy < width
{	float sum = 0;

			// Iteration through  gaussian kernel
	for(int k=0; k<kerw;k++)
		{
		
			{
			int offseti = k_height>1?(-1*(k_height/2)+k):0; 
			int offsetj =k_width>1? (-1*(k_width/2)+k):0;
					
					
			if((tidx+offseti>=0) && (tidx+offseti<height) && (tidy+offsetj)>=0 && (tidy+offsetj)<width)
				{
				sum+= (image[(tidx+offseti)*width+tidy+offsetj])*(g_kernel[k]);
						

				}
			}
		}
		*(output+(tidx*width)+tidy)=sum;
}*/


}
