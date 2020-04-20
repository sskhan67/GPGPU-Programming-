#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include<cuda.h>
#include <thrust/device_vector.h>

#include <thrust/host_vector.h>

#include <thrust/sort.h>

#include <thrust/copy.h>


__global__

void hysteresis(float *sup, float *hyst,float t_high,float t_low,float *d_global_h, int height, int width)

{
	
	
	
		int tidx,tidy;
		tidx=threadIdx.x+blockIdx.x*blockDim.x;
                tidy=threadIdx.y+blockIdx.y*blockDim.y;
		
                hyst[tidx*width+tidy]=sup[tidx*width+tidy];
		
		float global=0;

		//thrust::device_ptr<float> thr_d(sup);

		//thrust::device_vector<float>sup(thr_d,thr_d+(height*width));

		//thrust::sort(sup.begin(),sup.end());


		//copy hyst to sorted varibale 
                //sorted[tidx*width+tidy]=hyst[tidx*width+tidy]; 
		//thrust::sort(sorted);
		//int s=round(0.95*height*width);
		//int t_high=sorted[s];
		//t_low=t_high/5;
		
		//int index = (int)((0.95)*height*width);


		//float t_high = sup[index];

		//float t_low =t_high/5;


		//printf("high:",&t_high);
                if(tidx<height && tidy <width)
			{
                            if(sup[tidx*width+tidy]>=t_high)
				{
                                                    hyst[tidx*width+tidy]=255;
						   global=global+2;
				}			
                            else if(sup[tidx*width+tidy]<=t_low)
                                                    {
						    hyst[tidx*width+tidy]=0;
						   global=global+2;
							}
                            else if(sup[tidx*width+tidy]<t_high && sup[tidx*width+tidy]>t_low)
					{
                                                    hyst[tidx*width+tidy]=125;
						   global=global+2;
					}
		

			d_global_h[tidx*width+tidy]=global;
                        }
}
