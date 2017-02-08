//
//  main.c
//  Problem_3
//  This program is write to do Crank-Nicolson method for PDF problem
//  Inputs are:
//  N: number of point.
//  dt: timestep for iteration.
//  T: totally time
//  The outputs are:
//  X
//  Result U when Time=0;
//  Result U when Time=0.5;
//  Result U when Time=1;
//  Result U when Time=1.5;
//  Result U when Time=2;
//  Created by Alex CONG on 13/4/2016.
//  Copyright © 2016 Alex CONG. All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <Accelerate/Accelerate.h>
//#include <lapacke.h>

double f(double in)
{
    return (in*in)/2;
}
void print_array(double *b,int n,FILE* outfile)
{
    int i;
    for (i=0;i<n;i++)
    {
        printf("%f\t",b[i]);
    }
    printf ("\t\n");
}
void burger2(double *u,int N,double dt,double T,double dx,FILE* outfile)
{
    
    
    int i,j;
    int one=1;
    int info;
    char tr='N';
    double pi=M_PI;
    double nu,lamda,timer;
    int nsteps,matrix_Di;
    double* uhat=(double*)malloc(N*sizeof(double));
    
    double* B=(double*)malloc((N-2)*sizeof(double));
    double* ld=(double*)malloc((N-2)*sizeof(double));
    double* d=(double*)malloc((N-2)*sizeof(double));
    double* ud=(double*)malloc((N-2)*sizeof(double));
    double* uud=(double*)malloc((N-2)*sizeof(double));
    int* pivot = malloc((N-2)*sizeof(int));
    
    matrix_Di=N-2;
    nu = 0.01/pi;
    nsteps = round(T/dt);
    lamda=nu*dt*0.5/(dx*dx);
    // initial all the arrays
    for (i=0;i<N;i++)
    {
        uhat[i]=0.;
    }
    for (i=0;i<N-2;i++)
    {
        B[i]=0.;
        ld[i]=0.;
        d[i]=0.;
        ud[i]=0.;
        uud[i]=0.;
        pivot[i]=0;
    }
    
    // set up the d,ld,ud,uud
    for (i=0;i<N-2;i++)
    {
       ld[i]=-lamda;
       ud[i]=-lamda;
       d[i]=1+2*lamda;
    }
    ud[N-3]=0.;
    ld[N-3]=0.;
   // upper triagnoal sovle
    dgttrf_(&matrix_Di,ld, d, ud, uud, pivot, &info);
    

    for (i=0;i<nsteps;i++)
    {
    
        
        if ((i+1)%2==0)
        {
            for (j=1;j<N-1;j++) // crank-nicolson FB
            {
                uhat[j] = u[j] - dt/dx*(f(u[j+1])-f(u[j]));
            }
            for (j=1;j<N-1;j++) // crank-nicolson FB
            {
                u[j] = 0.5*(uhat[j]+u[j]-dt/dx*(f(uhat[j])-f(uhat[j-1])));
            }
        }
        else
        {
            for (j=1;j<N-1;j++) // crank-nicolson BF
            {
                uhat[j] = u[j] - dt/dx*(f(u[j])-f(u[j-1]));
            }
            for (j=1;j<N-1;j++) // crank-nicolson BF
            {
                u[j] = 0.5*(uhat[j]+u[j]-dt/dx*(f(uhat[j+1])-f(uhat[j])));
            }
        }
        // solve the implicit matrix
    
        for (j=1;j<N-1;j++)
        {
            B[j-1]=lamda*(u[j-1]+u[j+1])+(1-2.0*lamda)*u[j];
        }
        // solve the matrix
        dgttrs_(&tr, &matrix_Di, &one, ld, d, ud, uud, pivot, B, &N, &info);
        

        for (j=1;j<N-1;j++)
        {
            u[j]=B[j-1];
        }

        
        if ((i+1)%2000==0)
        {
            timer=i+1;
            printf ("The u on time %f is:\t\n",timer/2000/2);
            print_array(u,N,outfile);
            fwrite(u, sizeof(double), N, outfile);
            printf ("\t\n");
        }
    }
    free(uhat);
    free (B);
    free (ld);
    free (ud);
    free (uud);
    free (d);
    free (pivot);
    
}
int main(int argc, char* argv[])
{
     if (argc != 3) {
     printf("Incorrect usage: only enter the input data file name\n");
     return 0;
     }
     FILE* inputfile = fopen(argv[1], "r");

    if (!inputfile) {
        printf("Unable to open input file\n");
        return 0;
    }
    // start reading input data using function fscanf here

    FILE* outfile=fopen(argv[2], "w");
   
    int N;
    fscanf(inputfile, "%d", &N); // read an integer N for example
    // read rest of parameters here
    double dt;
    fscanf(inputfile, "%lf", &dt);
    double T;
    fscanf(inputfile, "%lf", &T);
    fclose(inputfile);
    
    double pi=M_PI;
    double dx,dtest;
    double* x=(double*)malloc(N*sizeof(double));
    double* u=(double*)malloc(N*sizeof(double));
    int i;
    
    
    clock_t start = clock();
    
    // ... rest of program ...
    printf("The N is %d.\t\n",N);  // notify output N
    fwrite(&N, sizeof(int), 1, outfile); // output N
    
    dtest=N;  // choose a different N for double
    dx=2/(dtest-1);
    
    printf("The x is:\t\n");  // notify output x

    for (i=0;i<N;i++)
    {
        x[i]=-1+i*dx;
        u[i]=-sin(pi*x[i]);
        printf("%f\t",x[i]);   // output the x
    }
    fwrite(x, sizeof(double), N, outfile);
    printf ("\t\n");  // tab and enter into another line


    printf ("The u on time %f is:\t\n",0.);
    print_array(u,N,outfile);
    fwrite(u, sizeof(double), N, outfile);
    
    burger2(u,N,dt,T,dx,outfile);
    
    free (x);
    free (u);
    
    printf("Time elapsed: %g seconds\n", (float)(clock()-start)/CLOCKS_PER_SEC);
    fclose(outfile);
    return 0;
}

