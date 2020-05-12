
//author : Syeduzzaman Ikhan
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "image_template.h"
#include "histogramEqualisation.h"
#include<cuda.h>
#include "sys/time.h"

// main function starts here 
int main(int argc, char **argv)
{
		
		

		// Host variable 
		float *image,*temp,*h_hist,*h_cdf;	//  create pointer variable 
		int height, width;// image height, width, and hernel width size 
		
		// device variable
		float *d_image,*d_temp,*d_hist,*d_cdf;
		
// step 1. call function to read image and pass the funtion to the function 
                read_image_template(argv[1],&image,&width,&height);



		//Malloc for Host 
		temp = (float*)malloc(sizeof(float)*height*width);
		h_hist=(float*)malloc(sizeof(float)*256);
		h_cdf=(float*)malloc(sizeof(float)*256);
		
		

		// Malloc for DEVICE GPU
		cudaMalloc((void **)&d_image,sizeof(float)*width*height);
		cudaMalloc((void **)&d_temp,sizeof(float)*width*height);
		cudaMalloc((void **)&d_hist,sizeof(float)*256);
		cudaMalloc((void **)&d_cdf,sizeof(float)*256);

		cudaMemset(d_hist, 0, 256*sizeof(float));		
		cudaMemset(d_cdf, 0, 256*sizeof(float));



        	//copy the items from CPU to GPU
		cudaMemcpy(d_image,image,sizeof(float)*width*height,cudaMemcpyHostToDevice);
		int block_dim = 16;
		dim3 dmBlock(block_dim, block_dim, 1);
		dim3 dmGrid(ceil(height/block_dim), ceil(width/block_dim), 1);
		





		// global memory Kernel , histrogram 
		hist1<<<dmGrid,dmBlock>>>(d_image,d_hist,width,height);
		cudaDeviceSynchronize();
		
		
	//normalization 
		
	for (int i=0;i<256;i++)
		{	
			
        		
			h_hist[i]=h_hist[i]/(height*width);
	
			
		}
		
		h_cdf[0]=h_hist[0];

		//CDF calculation 

		for (int i=1;i<256;i++)
		{
        		h_cdf[i]=h_cdf[i-1]+h_hist[i];
	
		}
	
 		cudaMemcpy(d_cdf,h_cdf,sizeof(float)*256,cudaMemcpyHostToDevice);
		
		


		//Hist Normalization kernel 
		hist_f<<<dmGrid,dmBlock>>>(d_image,d_cdf,width,height);
        cudaDeviceSynchronize();






		cudaMemcpy(temp,d_image,sizeof(float)*width*height,cudaMemcpyDeviceToHost);


		// write image on disk 

		write_image_template("histrogram.pgm",temp, width, height);

		cudaFree(image);
		cudaFree(d_image);
		cudaFree(temp);
		cudaFree(d_temp);


}



