
//author : Syeduzzaman Ikhan
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "image_template.h"
#include "convolution.h"
#include"gaussian_kernel.h"
#include"gradient.h"
#include"suppress.h"
#include"hysteresis.h"
#include"edge.h"
#include"cornerness.h"
#include "sys/time.h"
#include <cuda.h>
#include<thrust/host_vector.h>
#include<thrust/device_vector.h>
#include<thrust/sort.h>


/*
Project Goal: To implement Canny Edge detection's 7 Feature detection and Parallelization using Open MP 

I/O File :  Test Image
Main function: Takes two argument as image path and sigma value=0.6 or 1 or 1.25

Setps: 

	User defined Header Files: image.h -> Image I/O handler, gaussian_kernel.h-> create gaussian kernel using user defined sigma value =0.6 or 1 or 1.25, convolution.h-> To perform
	convolution using Gaussian and Gaussian derivative kernel,  gradient.h-> perform horizontal and vertical gradient , boundary.h-> to handle boudary for convolution  
 hystersis+ suppress+edges 



*/


// main function starts here 
int main(int argc, char **argv)
{
		struct timeval start,end,start_k1,end_k1,start_k2,end_k2,start_k3,end_k3,start_k4,end_k4,start_k5,end_k5,start_c1,end_c1,start_c2,end_c2,start_c3,end_c3,start_nf,end_nf,start_f1,end_f1,start_f2,end_f2;

		
		// Host variable 
		float *image,*g_kernel, *deriv_kernel, *temp_hor, *temp_ver,  *horizontal_grad, *vertical_grad, *magnitude, *phase,*sup,*hyst,*edge,*corner,*feature;	//  create pointer variable 
		int height, width, k_width;// image height, width, and hernel width size 
		float sigma=atof(argv[2]);
		int a=round(2.5*sigma-0.5);
		k_width=2*a+1;
		
			//global memory access variable

		float *flop,*global,*flop1,*global1,*global2,*flop2,*flop3,*global4,*global_h;

		// Device variable 
float *d_image,*d_temp_hor,*d_horizontal_grad,*d_temp_ver,*d_vertical_grad,*d_magnitude,*d_phase,*d_sup,*d_hyst,*d_edge,*d_corner,*d_feature;
float *d_kernel,*d_deriv_kernel;

float *d_flop,*d_global,*d_flop1,*d_global1,*d_global2,*d_flop2,*d_flop3,*d_global3,*d_global_h;

// with file I/O time start
gettimeofday(&start,NULL);



// start file I/O time 
gettimeofday(&start_f1,NULL);

        	// step 1. call function to read image and pass the funtion to the function 
		read_image_template(argv[1],&image,&width,&height);
// end time for File I/O
gettimeofday(&end_f1,NULL);



//Malloc for Host 
temp_hor = (float*)malloc(sizeof(float)*height*width);
horizontal_grad = (float*)malloc(sizeof(float)*height*width);
temp_ver = (float*)malloc(sizeof(float)*height*width);
vertical_grad = (float*)malloc(sizeof(float)*height*width);
magnitude = (float*)malloc(sizeof(float)*height*width);
phase = (float*)malloc(sizeof(float)*height*width);
sup = (float*)malloc(sizeof(float)*height*width);
hyst = (float*)malloc(sizeof(float)*height*width);
edge = (float*)malloc(sizeof(float)*height*width);
corner = (float*)malloc(sizeof(float)*height*width);
feature = (float*)malloc(sizeof(float)*height*width);

flop = (float*)malloc(sizeof(float)*height*width);
global = (float*)malloc(sizeof(float)*height*width);

flop1 = (float*)malloc(sizeof(float)*height*width);
global1 = (float*)malloc(sizeof(float)*height*width);

flop2= (float*)malloc(sizeof(float)*height*width);
global2 = (float*)malloc(sizeof(float)*height*width);

flop3= (float*)malloc(sizeof(float)*height*width);
global4 = (float*)malloc(sizeof(float)*height*width);
global_h = (float*)malloc(sizeof(float)*height*width);

		
		// Malloc for DEVICE GPU
cudaMalloc((void **)&d_image,sizeof(float)*width*height);
cudaMalloc((void **)&d_temp_hor,sizeof(float)*width*height);
cudaMalloc((void **)&d_horizontal_grad,sizeof(float)*width*height);
cudaMalloc((void **)&d_temp_ver,sizeof(float)*width*height);
cudaMalloc((void **)&d_vertical_grad,sizeof(float)*width*height);
cudaMalloc((void **)&d_magnitude,sizeof(float)*width*height);
cudaMalloc((void **)&d_phase,sizeof(float)*width*height);
cudaMalloc((void **)&d_sup,sizeof(float)*width*height);
cudaMalloc((void **)&d_hyst,sizeof(float)*width*height);
cudaMalloc((void **)&d_edge,sizeof(float)*width*height);
cudaMalloc((void **)&d_corner,sizeof(float)*width*height);
cudaMalloc((void **)&d_feature,sizeof(float)*width*height);
//udaMalloc((void **)&d_feature,sizeof(float)*width*height);



cudaMalloc((void **)&d_flop,sizeof(float)*width*height);
cudaMalloc((void **)&d_global,sizeof(float)*width*height);
cudaMalloc((void **)&d_flop1,sizeof(float)*width*height);
cudaMalloc((void **)&d_global1,sizeof(float)*width*height);
cudaMalloc((void **)&d_flop2,sizeof(float)*width*height);
cudaMalloc((void **)&d_global2,sizeof(float)*width*height);
cudaMalloc((void **)&d_flop3,sizeof(float)*width*height);
cudaMalloc((void **)&d_global3,sizeof(float)*width*height);
cudaMalloc((void **)&d_global_h,sizeof(float)*width*height);




cudaMalloc((void **)&d_kernel,sizeof(float)*k_width);
cudaMalloc((void **)&d_deriv_kernel,sizeof(float)*k_width);

		// step 2. Create gaussian kernel and gaussian derivative kernel using sigma 
		
		gaussian(&g_kernel, &deriv_kernel, atof(argv[2]), &k_width); // atof function converts sigma string value to type double 



	// time ->  communication time 

gettimeofday(&start_c1,NULL);
	
	//copy the items from CPU to GPU
cudaMemcpy(d_image,image,sizeof(float)*width*height,cudaMemcpyHostToDevice);
cudaMemcpy(d_kernel,g_kernel,sizeof(float)*k_width,cudaMemcpyHostToDevice);
cudaMemcpy(d_deriv_kernel,deriv_kernel,sizeof(float)*k_width,cudaMemcpyHostToDevice);
cudaDeviceSynchronize();
	// end communication time 

gettimeofday(&end_c1,NULL);

float c1=((end_c1.tv_sec*1000000+end_c1.tv_usec)-(start_c1.tv_sec*1000000+start_c1.tv_usec));


//Horizontal
int block_dim=atof(argv[3]);
//int block_dim = 32;
dim3 dmBlock(block_dim, block_dim, 1);
dim3 dmGrid(ceil(height/block_dim), ceil(width/block_dim), 1);


	// kernel time start 
gettimeofday(&start_k1,NULL);

/*
convoultion<<<dmGrid,dmBlock>>>(d_image, d_temp_hor,d_kernel,height,width,k_width,1);
convoultion<<<dmGrid,dmBlock>>>(d_temp_hor,d_horizontal_grad,d_deriv_kernel,height,width,1,k_width);
cudaDeviceSynchronize();

*/
convoultion<<<dmGrid,dmBlock,sizeof(float)*block_dim*block_dim>>>(d_image, d_temp_hor,d_kernel,d_flop1,d_global1,height,width,k_width,1);	


	// flop & global memory access calculation
	cudaMemcpy(flop1,d_flop1,sizeof(float)*width*height,cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();


	cudaMemcpy(global1,d_global1,sizeof(float)*width*height,cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();

	


convoultion<<<dmGrid,dmBlock,sizeof(float)*block_dim*block_dim>>>(d_temp_hor,d_horizontal_grad,d_deriv_kernel,d_flop1,d_global1,height,width,1,k_width); 
cudaDeviceSynchronize();



	    
// Vertical
/*
convoultion<<<dmGrid,dmBlock>>>(d_image, d_temp_ver,d_kernel,height,width,k_width,1);
convoultion<<<dmGrid,dmBlock>>>(d_temp_ver,d_vertical_grad,d_deriv_kernel,height,width,1,k_width);

*/

convoultion<<<dmGrid,dmBlock,sizeof(float)*block_dim*block_dim>>>(d_image, d_temp_ver,d_kernel,d_flop1,d_global1,height,width,k_width,1);
convoultion<<<dmGrid,dmBlock,sizeof(float)*block_dim*block_dim>>>(d_temp_ver,d_vertical_grad,d_deriv_kernel,d_flop1,d_global1,height,width,1,k_width);


cudaDeviceSynchronize();
	
	//kernel time end

gettimeofday(&end_k1,NULL);

float k1=((end_k1.tv_sec*1000000+end_k1.tv_usec)-(start_k1.tv_sec*1000000+start_k1.tv_usec));
printf("\nHorizontal and Vertical Kernel time(ms): %f\n",k1/1000);

// Magnitude & Phase
 
	// kernel time start
gettimeofday(&start_k2,NULL);

gradient<<<dmGrid,dmBlock>>>(d_vertical_grad,d_horizontal_grad,d_magnitude,d_phase,d_flop2,d_global2,height,width);
cudaDeviceSynchronize();


	// flop & global memory access calculation
        cudaMemcpy(flop2,d_flop2,sizeof(float)*width*height,cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();


        cudaMemcpy(global2,d_global2,sizeof(float)*width*height,cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();


	// kernel time end
gettimeofday(&end_k2,NULL);

float k2=((end_k2.tv_sec*1000000+end_k2.tv_usec)-(start_k2.tv_sec*1000000+start_k2.tv_usec));
printf("Mag & Phase  Kernel time(ms): %f\n",k2/1000);


// Suppression

	// kernel time start

gettimeofday(&start_k3,NULL);

suppress<<<dmGrid,dmBlock>>>(d_magnitude,d_phase,d_sup,d_flop3,d_global3,height,width);
cudaDeviceSynchronize();


	
	// flop & global memory access calculation
        
	
	cudaMemcpy(flop3,d_flop3,sizeof(float)*width*height,cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();


        cudaMemcpy(global4,d_global3,sizeof(float)*width*height,cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();




	//kernel time end
gettimeofday(&end_k3,NULL);

float k3=((end_k3.tv_sec*1000000+end_k3.tv_usec)-(start_k3.tv_sec*1000000+start_k3.tv_usec));
printf("Supression Kernel time(ms): %f\n",k3/1000);



thrust::device_ptr<float> thr_d(d_sup);

thrust::device_vector<float>sup_vec(thr_d,thr_d+(height*width));

thrust::sort(sup_vec.begin(),sup_vec.end());

int index = (int)((0.95)*height*width);

float t_high = sup_vec[index];

float t_low =t_high/5;


// hysteresis and edge linking 


	//kernel time start
gettimeofday(&start_k4,NULL);

finaledge<<<dmGrid,dmBlock>>>(d_sup,d_edge,d_hyst,t_high,t_low,d_global_h,height,width);
cudaDeviceSynchronize();

	//kernel time end
gettimeofday(&end_k4,NULL);


	 // flop & global memory access calculation


        //cudaMemcpy(flop3,d_flop3,sizeof(float)*width*height,cudaMemcpyDeviceToHost);
        //cudaDeviceSynchronize();


        cudaMemcpy(global_h,d_global_h,sizeof(float)*width*height,cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();





float k4=((end_k4.tv_sec*1000000+end_k4.tv_usec)-(start_k4.tv_sec*1000000+start_k4.tv_usec));
printf("Hyster. & Final edge Kernel time(ms): %f\n",k4/1000);




	//start communication time
gettimeofday(&start_c2,NULL);


cudaMemcpy(edge, d_edge, sizeof(float)*width*height, cudaMemcpyDeviceToHost);
cudaDeviceSynchronize();

	//end communication time 
gettimeofday(&end_c2,NULL);

float  c2=((end_c2.tv_sec*1000000+end_c2.tv_usec)-(start_c2.tv_sec*1000000+start_c2.tv_usec));


// CORNERESS funtion 

	//kernel time start
gettimeofday(&start_k5,NULL);

feature_detector<<<dmGrid,dmBlock,2*sizeof(float)*block_dim*block_dim>>>(d_corner,d_global,d_flop,height, width, d_vertical_grad, d_horizontal_grad, block_dim);
//cudaMemcpy(corner,d_corner,sizeof(float)*width*height,cudaMemcpyDeviceToHost);
        
find_featureGPU<<<dmGrid,dmBlock,2*sizeof(float)*block_dim*block_dim>>>(d_feature,d_corner,height, width,block_dim);


 	//kernel time end
gettimeofday(&end_k5,NULL);

float k5=((end_k5.tv_sec*1000000+end_k5.tv_usec)-(start_k5.tv_sec*1000000+start_k5.tv_usec));
printf(" Corner & Feature  detection  Kernel time(ms): %f\n",k5/1000);



	 //start communication time
gettimeofday(&start_c3,NULL);
cudaMemcpy(feature,d_feature,sizeof(float)*width*height,cudaMemcpyDeviceToHost);
cudaDeviceSynchronize();

//end communication time 
gettimeofday(&end_c3,NULL);

float  c3=((end_c3.tv_sec*1000000+end_c3.tv_usec)-(start_c3.tv_sec*1000000+start_c3.tv_usec));


// CSV File writting for Features 
	FILE *fp ;
    	fp=fopen("corner.csv","w+");
	 fprintf(fp,"i,j");
	int LocI,LocJ;
	for (int i=0;i<height*width;i++)
	{
	if (*(feature+i)>0)
		{

		int a=*(feature+i);
		LocI=a/width;
		LocJ=a % width;
		
		fprintf(fp,"\n%d, %d\n",LocI,LocJ);
		}
	}
	// Flop and global for Feature detetector  

cudaMemcpy(flop,d_flop,sizeof(float)*width*height,cudaMemcpyDeviceToHost);
cudaDeviceSynchronize();


cudaMemcpy(global,d_global,sizeof(float)*width*height,cudaMemcpyDeviceToHost);
cudaDeviceSynchronize();


// flops

float sum_ff=0,sum_f2=0,sum_f3,sum_f4;

for(int i=0;i<height*width;i++)
{

	if(*(flop+i)>0)
	{	
	sum_ff=sum_ff+flop[i];
	sum_f2=sum_f2+flop1[i];
	sum_f3=sum_f3+flop2[i];


	}	
}

printf("------Flop Calculation-------\n");
printf("\n Horizontal Flops: %f\n",4*sum_f2);

printf("\n Mag & Phase Flops: %f\n",sum_f3);


printf("\n Feature Flopsa: %f\n",sum_ff);


printf("------Memory access  Calculation-------\n");


// global memory

float sum_gf=0,sum_g2=0,sum_g3=0,sum_g4=0,sum_g5=0;

for(int i=0;i<height*width;i++)
{

        if(*(global+i)>0)
        {
        sum_gf=sum_gf+global[i];
        
        sum_g2=sum_g2+global1[i];
        sum_g3=sum_g3+global2[i];
        //sum_g4=sum_g4+global4[i];
        //sum_g5=sum_g5+global_h[i];



        }
}

printf("Horizontal global : %f\n",4*sum_g2);
printf("Mag & Phase Global : %f\n",sum_g3);
printf("feature_global: %f\n\n",sum_gf);


printf("-------------\n");







	// File I/O time start
gettimeofday(&start_f2,NULL);



write_image_template("edge.pgm",edge, width, height);

	// file i/o time over 
gettimeofday(&end_f2,NULL);


	// end to end time
gettimeofday(&end,NULL);

float time_io=((end_f1.tv_sec*1000000+end_f1.tv_usec)-(start_f1.tv_sec*1000000+start_f1.tv_usec))+   ((end_f2.tv_sec*1000000+end_f2.tv_usec)-(start_f2.tv_sec*1000000+start_f2.tv_usec));
float time_k=(k1+k2+k3+k4)/1000;

float time_c=(c1+c2+c3)/1000;
float time_t1=((end.tv_sec*1000000+end.tv_usec)-(start.tv_sec*1000000+start.tv_usec));
printf("Image width: %d,Image Height: %d, Sigma: %f,Block size: %d, Kernel time (ms): %f, communication time (ms): %f, Parallel time (ms): %f, file i/o time: %f, end-to-end time with file i/o: %f\n",width,height,sigma,block_dim,time_k,time_c,time_k+time_c,time_io/1000,time_t1/1000);

printf("C1: %f",c1/1000);
printf("C2: %f",c2/1000);

	cudaFree(image);
 	cudaFree(d_image);
	cudaFree(d_temp_hor);
	cudaFree(d_horizontal_grad);
	cudaFree(d_temp_ver);
    	cudaFree(d_vertical_grad);
	cudaFree(d_magnitude);
	cudaFree(d_phase);	
	cudaFree(sup);
	cudaFree(hyst);
	cudaFree(edge);
	cudaFree(d_corner);
	cudaFree(d_feature);
        cudaFree(d_global);
        cudaFree(d_global1);
        cudaFree(d_global2);
        //cudaFree(d_global4);

        cudaFree(d_flop);
        cudaFree(d_flop1);
        cudaFree(d_flop2);



// Tp= communication time +kernel time 
// Ts= compution time without File I/O 

}



