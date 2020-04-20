#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include<cuda.h>

__global__


void suppress(float *mag, float *phase, float *sup, float *d_flop3,float *d_global3, int height, int width)
{
	

	int tidx,tidy;
        tidx=blockIdx.x*blockDim.x + threadIdx.x;
        tidy=blockIdx.y*blockDim.y + threadIdx.y;
	
	// copy magnitude to supression, so that we can be supressed pixel from magnitude using several condtions
	sup[tidx*width+tidy]=mag[tidx*width+tidy];

	float flop=0,global=0;


		if (tidx<height && tidy<<width)
		{
			float theta = phase[tidx*width+tidy];
			
			global++;

			if(theta<0)
				theta+=M_PI;
			
			theta*=(180/M_PI);			
			
			if(theta<=22.5 || theta >157.5)
			{
				if (tidx-1 >= 0 && tidx+1 < height)
				 {
                                	if (mag[(tidx-1)*width+tidy]>mag[tidx*width+tidy] || mag[(tidx+1)*width+tidy]>mag[tidx*width+tidy])
						{

							sup[tidx*width+tidy]=0; 
							
							global=global+5;
							//printf("g: %f\n",global);

						}
				}
			}
			else if(theta>22.5 && theta<=67.5)
			{
				
                            if ( (tidx-1 >= 0 && tidy-1 >= 0) || ( tidx+1 < height && tidy+1 < width))
				 {
                               		if (mag[(tidx-1)*width+(tidx-1)]> mag[tidx*width+tidy] || mag[(tidx+1)*width+(tidy+1)]>mag[tidx*width+tidy])
						{
							sup[tidx*width+tidy]=0;
                                                        global=global+5;
   
 
						}
				}
			}
			else if(theta>67.5 && theta<=112.5)
			{
				if (tidx-1 >= 0 && tidy+1 < width)
				 {
                                	if (mag[tidx*width+(tidy-1)]>mag[tidx*width+tidy] || mag[tidx*width+(tidy+1)]>mag[tidx*width+tidy])
					 {
                                     		sup[tidx*width+tidy]=0;
                              			                                                        global=global+5;

					}
				}
                            
			}

			else if(theta>112.5 && theta<=157.5)
			{
			

				if ((tidy-1 >= 0 && tidx-1 >= 0 ) &&(tidx+1 < height && tidy+1 < width))
					 {
                                		if (mag[(tidx+1)*width+(tidy-1)]>mag[tidx*width+tidy] || mag[(tidx-1)*width+(tidy+1)]>mag[tidx*width+tidy])
							 {
  
                                   				sup[tidx*width+tidy]=0;

								                                                        global=global+5;

//								printf("g: %f\n",global);
								d_global3[tidx*width+tidy]= global; 
                                			}
                            		}
                        }

			
//			d_global3[tidx*width+tidy]= global;	
		}
//                        d_global3[tidx*width+tidy]= global;

}
