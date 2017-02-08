//
//  main.c
//  problem_2
//
//  Created by Alex CONG on 6/6/2016.
//  Copyright © 2016 Alex CONG. All rights reserved.
//  Update from problem_2 of project 3
//  Generate the matrix of U[y][x]
//  With the boundary condition u(-1,y,t)=y, u(1,y,t)=-y, u(x,-1,t)=x, u(x,1,t)=-x
//  Apply ADI method to find ut=nu*(uxx+uyy);
//  in kernel
//  {
//  explicit the localU
//  transpose them
//  Lapack to solve diagonal matrix
//  explicit the localU
//  transpose them
//  Lapack to solve diagonal matrix
//  }
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <magma.h>
#include <magma_types.h>
#include <magma_lapack.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

__device__ static int dev_N;            // global variable for dimension N in kernel
__device__ static double dev_lamda;     // lamda is coefficent nu*dt/2/(dx*dx)

const int gendimT=1024;// number of general threads number

__global__ void solve(double* u, double* unew){
__shared__ double localu[gendimT]; // local memory
int row,col; // row and column for the local u
int id_local;
int id=blockIdx.x*dev_N+threadIdx.x;
if (threadIdx.x<dev_N && blockIdx.x<dev_N)
{
id_local=threadIdx.x;
row=id/dev_N;
col=id%dev_N;
localu[id_local]=u[id];
__syncthreads();

if (row!=0 && row!=dev_N-1 && col!=0 && col!=dev_N-1)
{
unew[id]=localu[id_local]+dev_lamda*(localu[id_local+1]-2*localu[id_local]+localu[id_local-1]);
}
}
}

__global__ void transpose(double* u,double* unew){
__shared__ double localu[gendimT]; // local memory
int id=blockIdx.x*dev_N+threadIdx.x;
if (threadIdx.x<dev_N && blockIdx.x<dev_N)
{
int row=id/dev_N;
int col=id%dev_N;
localu[threadIdx.x]=u[id];
__syncthreads();
unew[col*dev_N+row]=localu[threadIdx.x];
}
}


int main(int argc,char * argv[]) {

//initial the cuda
cudaDeviceProp prop;
int dev;
memset( &prop, 0, sizeof(cudaDeviceProp));
prop.multiProcessorCount = 13;
prop.major = 3;
prop.minor = 5;
cudaChooseDevice(&dev, &prop);
cudaSetDevice(dev);

//create the event and record times
cudaEvent_t start,stop;
cudaEventCreate(&start);
cudaEventCreate(&stop);
cudaEventRecord(start,0);
float elapsedTime;

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

double dx,dtest;
double lamda;
int i,j; // counter
int nsteps; // time step

fwrite(&N, sizeof(int), 1, outfile); // output N

dtest=N;
dx=2/(dtest-1); // find the distance for two x grid points
nsteps=round(T/dt); // how many steps for the whole process
lamda=nu*dt/2/(dx*dx);  // generate the lamda

// give number to global kernal dev_N and dev_lamda
cudaMemcpyToSymbol(dev_N, &N, sizeof(int));
cudaMemcpyToSymbol(dev_lamda, &lamda, sizeof(double));

// initial the matrix and copy the matrix into kernel
double *x=(double*)malloc(N*sizeof(double));
double *u=(double*)malloc(N*N*sizeof(double));
double *dev_unew;
double *dev_uold;

cudaMalloc((void **)&dev_unew, N*N*sizeof(double));
cudaMalloc((void **)&dev_uold, N*N*sizeof(double));

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


// Initialize the MAGMA system
magma_init();
magma_int_t *piv, info;
magma_int_t dimen_N=N;
magma_int_t dimen_b=N-2;

double* A;                  // matrix Ax=b
double* dev_A;              // matrix A for kernel
piv=(magma_int_t*)malloc(dimen_N*sizeof(magma_int_t));
// assign memory for the cuda lapack
magma_dmalloc_cpu(&A,dimen_N*dimen_N);
magma_dmalloc(&dev_A,dimen_N*dimen_N);

for (i=0;i<N*N;i++)
{
A[i]=0.;
}
for(i=1; i<N-1; i++)
{
// wipe out the first row and last row
A[N*(i-1)+i]=-lamda;
A[N*(i)+i]=1+2*lamda;
A[N*(i+1)+i]=-lamda;
}
// the head and tail for the A
A[0]=1.0;
A[N*N-1]=1.0;

// Send this matrix to the device
magma_dsetmatrix(dimen_N,dimen_N,A,dimen_N,dev_A,dimen_N);
// Get the first part of solver
magma_dgetrf_gpu(dimen_N,dimen_N,dev_A,dimen_N,piv,&info);

// initial the state
cudaMemcpy(dev_uold, u, N*N*sizeof(double), cudaMemcpyHostToDevice);
cudaMemcpy(dev_unew, u, N*N*sizeof(double), cudaMemcpyHostToDevice);
// Main loop
for(i=0;i<nsteps;i++)
{
// step A
solve<<<N, N>>>(dev_uold, dev_unew);
transpose<<<N,N>>>(dev_unew,dev_uold);
// step B
magma_dgetrs_gpu(MagmaNoTrans,dimen_N,dimen_b,dev_A,dimen_N,piv,&(dev_uold[N]),dimen_N,&info);
// step C
solve<<<N, N>>>(dev_uold, dev_unew);
transpose<<<N,N>>>(dev_unew,dev_uold);
// step D
magma_dgetrs_gpu(MagmaNoTrans,dimen_N,dimen_b,dev_A,dimen_N,piv,&(dev_uold[N]),dimen_N,&info);
if((i+1)%1000==0)
{
// output the intended result
cudaMemcpy(u, dev_uold, N*N*sizeof(double), cudaMemcpyDeviceToHost);
fwrite(u, sizeof(double), N*N, outfile);
}
}
// output the time used
cudaEventRecord(stop,0);
cudaEventSynchronize(stop);
cudaEventElapsedTime(&elapsedTime,start,stop);
printf("The Elapsed Time is %f seconds \n",elapsedTime/(float)1000);

//close the output file
fclose(outfile);
//free memory
free(u);
free(x);
free(A);
magma_free(dev_A);
free(piv);
cudaFree(dev_unew);
cudaFree(dev_uold);
magma_finalize();
return 0;
}
