#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<cuda.h>
#include<time.h>
#include<sys/time.h>
#include<math.h>
int main(int argc, char ** argv){


	float *host, *device; 
	struct timeval start_timeS,end_timeE;
	int size;
	
	//int a = atof(argv[1]);	
	int a[]={9400,26200,59000,1,2,4,8,16,32,64,128,256,512,1000,2000,4000,8000,16000,32000,64000,128000,256000,512000};
	

//        int a[]={1,2,4,8,16,32,64,128,256,512,1000,2000,4000,8000,16000,32000,64000,128000,256000,512000};

	FILE *fp ;
        fp=fopen("d2h.csv","w+");
        fprintf(fp,"KB, time, BW");
	
for (int k=0;k<20;k++)
{	
	size = a[k]*1024; 
	host = (float*)malloc(sizeof(float)*size);
	cudaMalloc((void**)&device,sizeof(float)*size);

	for (int i = 0; i<size/4;i++){
		host[i]=rand();
	}


	cudaMemcpy(device,host,sizeof(float)*size,cudaMemcpyHostToDevice);

	cudaDeviceSynchronize();



    gettimeofday(&start_timeS,NULL);
    for (int i = 0; i<30; i++){
    cudaMemcpy(host,device,sizeof(float)*size,cudaMemcpyDeviceToHost);
    }
    gettimeofday(&end_timeE,NULL);

    cudaDeviceSynchronize();


	float timeM = ((end_timeE.tv_sec*1000000+end_timeE.tv_usec)-(start_timeS.tv_sec*1000000+start_timeS.tv_usec));
	float t=timeM/(1000*30);// time in ms
	float BW=a[k]/(t/1000); // kbps	
	printf("KB: %d, time from device to host : %f ms, BW[kpbs]: %f \n",a[k],t,BW);
	fprintf(fp,"\n%d, %f, %f \n",a[k],t,BW);

}
}
