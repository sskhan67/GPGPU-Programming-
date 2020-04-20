#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include<cuda.h>

// The function provides the feature detection algothim from Cannny edge detector, The feature detection is implemented to put extra computational stress 

	
__global__
	

void feature_detector(float *d_corner,float *d_flop,float *d_global,int height, int width, float *vertical, float *horizon,int blocksize)
{
                // varibale declaration 
                float k = 0.04;  // constant for cornernees formula 
                int window_width = 7;// window size 
		float Ixx,Iyy,IxIy;
                int count,locatI,locatJ,countK,countM,locatI1[100000],locatJ1[100000];

		

		// Device IDs
		int locaIx = threadIdx.x;
                int locaIy = threadIdx.y;
                int globIx = threadIdx.x+blockIdx.x*blockDim.x;
                int globIy = threadIdx.y+blockIdx.y*blockDim.y;
		
		float flop=0,global=0;

		extern __shared__ float Ashared[];
                __shared__ float *Vshared, *Hshared;
                Vshared = Ashared;
                Hshared = Ashared+blocksize*blocksize;

                Vshared[locaIx*blockDim.y+locaIy] = vertical[globIx*width+globIy];
                Hshared[locaIx*blockDim.y+locaIy] = horizon[globIx*width+globIy];
                        __syncthreads();
                        Ixx = 0;
                        Iyy = 0;
                        IxIy = 0;





		for(int k = -window_width/2; k < window_width/2 ; k++)
		{
                     for(int m = -window_width/2; m < window_width/2 ; m++)
			{

                         if(locaIx+k >= 0 && locaIx+k < blockDim.x && locaIy+m >= 0 && locaIy+m < blockDim.y)
			{
                                int offseti = locaIx+k;
                                int offsetj = locaIy+m;

                                Ixx = Ixx + pow(Vshared[offseti*blockDim.y+offsetj],2);
                                Iyy = Iyy + pow(Hshared[offseti*blockDim.y+offsetj],2);
                                IxIy = IxIy + Vshared[offseti*blockDim.y+offsetj] * Hshared[offseti*blockDim.y+offsetj];
                         
				flop=flop+6;
			}

                        else if(globIx+k >= 0 && globIx+k < height && globIy+m >= 0 && globIy+m < width)
			{
                                int offseti = globIx+k;
                                int offsetj = globIy+m;

                                Ixx = Ixx + pow(vertical[offseti*width+offsetj],2);
                                Iyy = Iyy + pow(horizon[offseti*width+offsetj],2);
                                IxIy = IxIy + vertical[offseti*width+offsetj] * horizon[offseti*width+offsetj];
                		
				global=global+4;

				flop=flop+6;	
			}
                        
			}
                      }
		  __syncthreads();
                  d_corner[globIx*width+globIy]= (Ixx*Iyy) - (IxIy*IxIy) - 0.04*(Ixx+Iyy)*(Ixx+Iyy);

		  d_flop[globIx*width+globIy]=flop;

		  global=global+1;


                  d_global[globIx*width+globIy]=global;





}

__global__

void find_featureGPU(float *output, float *d_cornerness,int height, int width, int blocksize)
{
	
                int locatI,locatJ,stride,localIndex;
                int window_size = blockDim.x;
		int window_width = blockDim.y;
             	int locaIx = threadIdx.x;
		int locaIy = threadIdx.y;
		
 		int globIx = threadIdx.x+blockIdx.x*blockDim.x;
		int globIy = threadIdx.y+blockIdx.y*blockDim.y;
		
		extern __shared__ float Shared[];
		__shared__ float *AShared,*indexShared;
                AShared = Shared;
		indexShared = Shared+blocksize*blocksize;

		
		AShared[locaIx*blockDim.y+locaIy] = d_cornerness[globIx*width+globIy];
		indexShared[locaIx*blockDim.y+locaIy] = globIx*width+globIy;
		__syncthreads();
               
		for ( stride = ((window_size*window_size)/2);stride >= 1; stride/=2)
        		{
			__syncthreads();
			localIndex = locaIx*blockDim.y+locaIy;
			if(localIndex < stride)
			{
		
        		if(AShared[localIndex]<AShared[localIndex+stride])
        		{
			AShared[localIndex]=AShared[localIndex+stride];
			indexShared[localIndex]=indexShared[localIndex+stride];
				
			
			}
			}
		
		
			}

			if(locaIx == 0 && locaIy == 0)
			{
			output[globIx*width+globIy]=indexShared[0];
			
			//int a[2000];
			int a = indexShared[0];
			int localI, localJ;
			localI=a/width;
			localJ=a%width;
			//FILE *fp;
			//printf("index= %d i=%d j=%d	",a,localI,localJ);

			}
}







