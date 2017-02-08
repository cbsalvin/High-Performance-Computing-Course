//
//  main.c
//  Problem_2
//  The program is designed to solve PDF problem with MacCormack’s method
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
//  Created by Alex CONG on 12/4/2016.
//  Copyright © 2016 Alex CONG. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

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
void burger(double *uhat,double *a,int N,double dt,double T,double dx,FILE* outfile)
{
    int i,j;
    double pi=M_PI;
    double nu,timer;
    int nsteps;
    
    nu = 0.01/pi;
    nsteps = round(T/dt);
    // initial the array
    for (i=0;i<N;i++)
    {
        uhat[i]=0.;
    }

  
    for (i=0;i<nsteps;i++)
    {
        if ((i+1)%2==0) //BF
        {
           for (j=1;j<N-1;j++)
           {
               uhat[j]=a[j]+(dt/dx)*(-f(a[j+1])+nu*(a[j+1]-a[j])/dx+f(a[j])-nu*(a[j]-a[j-1])/dx);
           }
           for (j=1;j<N-1;j++)
            {
                a[j]=0.5*(uhat[j]+a[j]+dt/dx*(-f(uhat[j])+nu*(uhat[j+1]-uhat[j])/dx+f(uhat[j-1])-nu*(uhat[j]-uhat[j-1])/dx));
            }
        }
        else
        {
            for (j=1;j<N-1;j++) // CF
            {
                uhat[j]=a[j]+ (dt/dx)*(-f(a[j])+nu*(a[j+1]-a[j])/dx+f(a[j-1])-nu*(a[j]-a[j-1])/dx);
            }
            for (j=1;j<N-1;j++)
            {
                a[j]=0.5*(uhat[j]+a[j]+dt/dx*(-f(uhat[j+1])+nu*(uhat[j+1]-uhat[j])/dx+f(uhat[j])-nu*(uhat[j]-uhat[j-1])/dx));
            }
    
        }
        if ((i+1)%50000==0)// print out and write out the output
        {
            timer=i+1;
            //printf ("The u on time %f is:\t\n",timer/50000/2);
            //print_array(a,N,outfile);
            fwrite(a, sizeof(double), N, outfile);
            //printf ("\t\n");
        }
    }


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
  
    clock_t start = clock();
    double* x=(double*)malloc(N*sizeof(double));
    double* u=(double*)malloc(N*sizeof(double));
    double* uhat=(double*)malloc(N*sizeof(double));
    double pi=M_PI;
    double dx,dtest;

    int i;
    // ... rest of program ...
    printf("The N is:%d. \t\n\n",N);  // notify output N
    //fprintf(outfile,"The N is:%d. \t\n\n",N);
    fwrite(&N, sizeof(int), 1, outfile); // output N
    
    dtest=N;
    dx=2/(dtest-1);
    
    //printf("The x is:\t\n");  // notify output x
    for (i=0;i<N;i++)
    {
        x[i]=-1+i*dx;
        u[i]=-sin(pi*x[i]);
        //printf("%f\t",x[i]);   // output the x
    }
    fwrite(x, sizeof(double), N, outfile); // output N
    //printf ("\t\n");  // tab and enter into another line
 
    
    //printf ("The u on time %f is:\t\n",0.);
    fwrite(u, sizeof(double), N, outfile);
    //print_array(u,N,outfile);
    
    burger(uhat,u,N,dt,T,dx,outfile);
    // free all the malloc
    free (x);
    free (u);
    free (uhat);
    
    printf("Time elapsed: %g seconds\n", (float)(clock()-start)/CLOCKS_PER_SEC);
    fclose(outfile);
    return 0;
}

