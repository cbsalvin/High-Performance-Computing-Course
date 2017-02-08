//
//  main.c
//  problem 1
//
//  Created by Alex CONG on 2/6/2016.
//  Copyright © 2016 Alex CONG. All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <curand.h>
#include <curand_kernel.h>


__global__ void mcquad(double* x, int N, int T, long int seed)
{
curandState_t state;
// generate the random seed
curand_init(blockIdx.x, 0,0, &state);
int j;
double temp;
int tid = blockIdx.x;
if (tid<T)
{
// run the Monte-Carlo
for(j=0;j<N;j++)
{
temp=curand(&state) /(float)(0x0FFFFFFFFUL);
x[tid]=x[tid]+cos(-log(temp));
}
x[tid]=x[tid]/N;
}

}
int main(int argc, char* argv[])
{
int N,T;
int i=0;


// create time function
float elapsedTime;
cudaEvent_t start,stop;
cudaEventCreate(&start);
cudaEventCreate(&stop);
cudaEventRecord(start,0);

FILE* inputfile = fopen(argv[1], "r");

fscanf(inputfile, "%d", &N); // read an integer N for example
fscanf(inputfile, "%d", &T);// read time period T
fclose(inputfile);
// initial the data
double x[T];
double* dev_x;
int dev_N;
int dev_T;

dev_N=N;
dev_T=T;

cudaMalloc((void**)&dev_x, T*sizeof(double));

for (i=0;i<T;i++)
{
x[i]=0.;
}

cudaMemcpy(dev_x, x, T*sizeof(double), cudaMemcpyHostToDevice);

mcquad<<<T,1>>>(dev_x,dev_N,dev_T,time(NULL));

cudaMemcpy(x, dev_x, T*sizeof(double), cudaMemcpyDeviceToHost);

// Time function
cudaFree(dev_x);
cudaEventRecord(stop,0);
cudaEventSynchronize(stop);
cudaEventElapsedTime(&elapsedTime,start,stop);
printf("The Elapsed Time is %f seconds \n",elapsedTime/(float)1000);

FILE* outfile=fopen(argv[2], "w");
fwrite(x, sizeof(double), T, outfile);
fclose(outfile);
return 0;
}

