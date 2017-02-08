//
//  main.c
//  problem_2
//
//  Created by Alex CONG on 5/6/2016.
//  Copyright © 2016 Alex CONG. All rights reserved.
//  Generate the matrix of U[y][x]
//  With the boundary condition u(-1,y,t)=y, u(1,y,t)=-y, u(x,-1,t)=x, u(x,1,t)=-x
//  Generate some blocks with threads to diff ut=nu*(uxx+uyy);
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cuda.h>


#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

__global__ void heat(double *u1, double *u2, double *dx, double *dt, double *nu, int *N)
{
int id=threadIdx.x+ blockIdx.x*blockDim.x;         //index in N*N array u

int id_xp=id+1;     // x+1
int id_xm=id-1;     // x-1
int id_yp=id-*N;    // y-1
int id_ym=id+*N;    // y+1

if(id<*N || id>=*N*(*N-1)) // boundary for y=-1 and y=1
{
u2[id]=u1[id];  // update the new matrix
}
else if(id%*N==0 || (id+1)%*N==0) // boundary for x=-1 and x=1
{
u2[id]=u1[id];  // update the new matrix
}
else        // update the difference condition
{

u2[id]=u1[id]+(*nu*(*dt)/(*dx*(*dx)))*(u1[id_xp]+u1[id_xm]-2*u1[id])+(*nu*(*dt)/(*dx*(*dx)))*(u1[id_yp]+u1[id_ym]-2*u1[id]);  //finite difference to calculate heat equation

}

}
int main(int argc,char * argv[]) {

// initial the cuda
cudaDeviceProp prop;
int dev;
memset(&prop, 0, sizeof(cudaDeviceProp));
prop.multiProcessorCount=13;
cudaChooseDevice(&dev, &prop);
cudaSetDevice(dev);
cudaGetDeviceProperties(&prop, dev);
int num_threads=prop.maxThreadsPerBlock; // number of threads
// create time function
cudaEvent_t start,stop;
cudaEventCreate(&start);
cudaEventCreate(&stop);
cudaEventRecord(start,0);
FILE* outfile=fopen(argv[2], "w");
// start reading input data using function fscanf here
int N;
double dt;
double T;
double nu;

FILE* inputfile = fopen(argv[1], "r");
fscanf(inputfile, "%d", &N); // read an integer N for example
fscanf(inputfile, "%lf", &dt);
fscanf(inputfile, "%lf", &T);
fscanf(inputfile, "%lf", &nu);
fclose(inputfile);

double dx,dtest; // delta_x
float elapsedTime; // time

int i,j; // counter
int nsteps; // time step

fwrite(&N, sizeof(int), 1, outfile); // output N

dtest=N;
dx=2/(dtest-1); // find the distance for two x grid points
nsteps = round(T/dt); // how many steps for the whole process

int k=0;             //k=0 is an index to mark the alternation of dev_u
int num_blocks=(N*N)/num_threads+((N*N)%num_threads?1:0);  //number of blocks


int *dev_N;
double *x=(double*)malloc(N*sizeof(double));
double *u=(double*)malloc(N*N*sizeof(double));
double *dev_u[2], *dev_dx, *dev_dt, *dev_nu;

cudaMalloc((void **)&dev_N, sizeof(int));
cudaMalloc((void **)&dev_nu, sizeof(double));
cudaMalloc((void **)&dev_dt, sizeof(double));
cudaMalloc((void **)&dev_dx, sizeof(double));
cudaMalloc((void **)&dev_u[0], N*N*sizeof(double));
cudaMalloc((void **)&dev_u[1], N*N*sizeof(double));
// initial the array x
for (i=0;i<N;i++)
{
x[i]=-1+i*dx;
}
// generate the matrix u

for(i=0;i<N;i++)
{
for(j=0;j<N;j++)
{
u[i+j*N]=-x[i]*x[j]+cos(11*M_PI*x[i]/2)*sin(8*M_PI*x[j]);
}
}

fwrite(x, sizeof(double),N,outfile);    // output x
fwrite(u, sizeof(double),N*N,outfile);  // output u

cudaMemcpy(dev_dt, &dt, sizeof(double), cudaMemcpyHostToDevice);
cudaMemcpy(dev_dx, &dx, sizeof(double), cudaMemcpyHostToDevice);
cudaMemcpy(dev_nu, &nu, sizeof(double), cudaMemcpyHostToDevice);
cudaMemcpy(dev_u[0], u, N*N*sizeof(double), cudaMemcpyHostToDevice);
cudaMemcpy(dev_u[1], u, N*N*sizeof(double), cudaMemcpyHostToDevice);
cudaMemcpy(dev_N, &N, sizeof(int), cudaMemcpyHostToDevice);

for(i=0;i<nsteps;i++)
{
heat<<<num_blocks, num_threads>>>(dev_u[k], dev_u[1-k], dev_dx, dev_dt, dev_nu, dev_N);
if((i+1)%4000==0)
{
// output the intended result
cudaMemcpy(u, dev_u[1-k], N*N*sizeof(double), cudaMemcpyDeviceToHost);
fwrite(u, sizeof(double), N*N, outfile);
}
k=1-k;
}

// output the time used.
cudaEventRecord(stop,0);
cudaEventSynchronize(stop);
cudaEventElapsedTime(&elapsedTime,start,stop);
printf("The Elapsed Time is %f seconds \n",elapsedTime/(float)1000);
// free memory
cudaFree(dev_u);
cudaFree(dev_dx);
cudaFree(dev_N);
free(x);
free(u);
fclose(outfile);
return 0;
}
