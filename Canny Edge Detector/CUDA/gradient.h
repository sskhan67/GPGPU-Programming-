#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include<cuda.h>

// Mag and phase calculation -> Input: Horizontal and Vertical Gradient, Output: Magnitude and Phase
__global__
void gradient(float *v_grad, float *h_grad, float *magnitude , float *phase,float *d_flop2, float *d_global2,  int height, int width)
{
	
	int tidx,tidy;
        tidx=blockIdx.x*blockDim.x + threadIdx.x;
        tidy=blockIdx.y*blockDim.y + threadIdx.y;

        
        // 1. Iteration through pixels
        //    if (tidx<height && tidy < width)	
	
	float flop=0, global=0; 


	if (tidx<height && tidy<width   )

	{
			magnitude[tidx*width+tidy]=sqrt( pow( v_grad[tidx*width+tidy], 2) + pow ( h_grad[tidx*width+tidy],2) );//magnitude
			phase[tidx*width+tidy]=atan2( v_grad[tidx*width+tidy],h_grad[tidx*width+tidy] );//Phase
	
			
			flop=flop+5;
			global=global+6;

			d_flop2[tidx*width+tidy]=flop;
                        d_global2[tidx*width+tidy]=global;

	}
	
	
}
