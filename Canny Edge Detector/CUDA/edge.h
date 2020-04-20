#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include<cuda.h>


__global__

void finaledge(float *sup, float *edges,float *hyst,float t_high, float t_low, float *d_global_h,  int height, int width)


{
	
		int tidx,tidy;
                tidx=threadIdx.x+blockIdx.x*blockDim.x;
                tidy=threadIdx.y+blockIdx.y*blockDim.y;

		float global=0;

                hyst[tidx*width+tidy]=sup[tidx*width+tidy];	
	

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
					//printf("global: %f",global);
					d_global_h[tidx*width+tidy]=global;
			}




		edges[tidx*width+tidy]=hyst[tidx*width+tidy];
		float global1=0;

			for (int y=-1; y<=1; y++)
				{
                                for (int x=-1; x<=1; x++)
					{

					
                                    	if(tidx+y<height && tidx+y>0 && tidy+x<width && tidy+x> 0)
						{
                                        		if (hyst[tidx*width+tidy]==255)
								{
                                            			edges[tidx*width+tidy]=255;
								global1=global1+2;
								}
                                        		else
								{
                                             			edges[tidx*width+tidy]=0;
								global1++;
								}
                                            	}
					}
				d_global_h[tidx*width+tidy]=global1;
				}


                              


	
}
