//
//  main.c
//  problem_2
//
//  Created by Alex CONG on 8/5/2016.
//  Copyright © 2016 Alex CONG. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>


// LAPACK function declarations
// These are here because there is no lapacke.h file on Quest
void dgttrf_(int*, double*, double*, double*, double*, int*, int*);
void dgttrs_(char*, int*, int*, double*, double*, double*, double*,
             int*, double*, int*, int*);
void dgtsv_(int*, int*, double*, double*, double*, double*, int*, int*);
//


int main(int argc, const char * argv[]) {
     MPI_Init(&argc, &argv);
     
     int rank, size;
     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     MPI_Comm_size(MPI_COMM_WORLD, &size);
     
     double precision = MPI_Wtick();
     double starttime = MPI_Wtime();
  

    FILE* inputfile = fopen(argv[1], "r");
    FILE* outfile=fopen(argv[2], "w");
    
    // start reading input data using function fscanf here
    // read rest of parameters here
    int N;
    fscanf(inputfile, "%d", &N); // read an integer N for example
    // read rest of parameters here
    double dt;
    fscanf(inputfile, "%lf", &dt);
    double T;
    fscanf(inputfile, "%lf", &T);
    double nu;
    fscanf(inputfile, "%lf", &nu);
    
    fclose(inputfile);

    //double* x=(double*)malloc(N*sizeof(double));
    //double* u=(double*)malloc(N*sizeof(double));
    //double* uhat=(double*)malloc(N*sizeof(double));
    
    
    double* x=(double*)malloc(N*sizeof(double));
    double* y=(double*)malloc(N*sizeof(double));
    double* new_x=(double*)malloc(N*sizeof(double));
    double new_y;
    double pi=M_PI;
    double dx,dtest;
    
    int i,j,k; // counter
    int localN,p; // local grid points & size P
    int nsteps; // time step
    
    p=sqrt(size);// number of processors for each dimension

    fwrite(&N, sizeof(int), 1, outfile); // output N
    
    dtest=N;
    dx=2/(dtest-1); // find the distance for two x grid points
    nsteps = round(T/dt); // how many steps for the whole process
    localN=N/p+(N%p>rank?1:0)+(rank>0?1:0)+(rank<p-1?1:0);
    
    double u[2][localN][localN]; // u[New/Old][x axix][y axix]
    
    // initial the matrix u
    for (i=0;i<N;i++)
    {
        x[i]=-1+i*dx;
        for (k=0;k<N;k++)
        {
            y[k]=-1+k*dx;
        }
    }
    for (i=0;i<localN;i++)
    {
        new_x[i]=x[(rank%p)*localN]+i*dx;
        for (k=0;k<localN;k++)
        {
            new_y=y[(rank/p)*localN]+k*dx;
            u[0][i][k]=-new_x[i]*new_y+cos(11*new_x[i]*pi/2)*sin(8*new_y*pi);
            u[1][i][k]=0.;
        }
    }
    
    // Boundary condition
    // boundary on x axix
    if (rank%p==0)
    {
        for (j=0;j<localN;j++)
            u[1][0][j]=u[0][0][j];
    }
    if (rank%p==p-1)
    {
        for (j=0;j<localN;j++)
            u[1][localN-1][j]=u[0][localN-1][j];
    }
    // boundary on y axix
    if (rank/p==0)
    {
        for (j=0;j<localN;j++)
            u[1][j][0]=u[0][j][0];
    }
    if (rank/p==p-1)
    {
        for (j=0;j<localN;j++)
            u[1][j][localN-1]=u[0][j][localN-1];
    }
    // storage for tracking communication info
    MPI_Request sendLeftRequest;
    MPI_Request sendRightRequest;
    MPI_Request recvLeftRequest;
    MPI_Request recvRightRequest;
    
    MPI_Request sendUpRequest;
    MPI_Request sendRightRequest;
    MPI_Request recvLeftRequest;
    MPI_Request recvRightRequest;
    
    // main loop
    for (i=0;i<nsteps;i++)
    {
        // x axix
        // send to left
        if (rank%p>0)
        {
            for (j=0;j<localN;j++)
            {
                MPI_Isend(&(u[i%2][1][j]),1,MPI_DOUBLE,rank-1,i,MPI_COMM_WORLD,&sendLeftRequest);
                MPI_Irecv(&(u[i%2][0][j]),1,MPI_DOUBLE,rank-1,i,MPI_COMM_WORLD,&recvLeftRequest);
            }
        }
        // send to right
        if (rank%p<p-1)
        {
            for (j=0;j<localN;j++)
            {
                MPI_Isend(&(u[i%2][p-2][j]),1,MPI_DOUBLE,rank+1,i,MPI_COMM_WORLD,&sendRightRequest);
                MPI_Irecv(&(u[i%2][p-1][j]),1,MPI_DOUBLE,rank+1,i,MPI_COMM_WORLD,&recvRightRequest);
            }
        }
        // y axix
        // send to up
        if (rank/p>0)
        {
            for (j=0;j<localN;j++)
            {
                MPI_Isend(&(u[i%2][j][1]),1,MPI_DOUBLE,rank-1,i,MPI_COMM_WORLD,&sendUpRequest);
                MPI_Irecv(&(u[i%2][j][0]),1,MPI_DOUBLE,rank-1,i,MPI_COMM_WORLD,&recvUpRequest);
            }
        }
        // send to down
        if (rank/p<p-1)
        {
            for (j=0;j<localN;j++)
            {
                MPI_Isend(&(u[i%2][j][p-2]),1,MPI_DOUBLE,rank+1,i,MPI_COMM_WORLD,&sendDownRequest);
                MPI_Irecv(&(u[i%2][j][p-1]),1,MPI_DOUBLE,rank+1,i,MPI_COMM_WORLD,&recvDownRequest);
            }
        }
        
        // update the interior where nothing to do with communication
        for (j=2;j<localN-2;j++)
        {
            for (k=2;k<localN-2;k++)
            {
                u[(i+1)%2][j][k]=u[i%2][j][k]+nu*dt/dx*(dt/dx)*(u[i%2][j+1][k]-2*u[i%2][j][k]+u[i%2][j-1][k])+nu*dt/dx*(dt/dx)*(u[i%2][j][k+1]-2*u[i%2][j][k]+u[i%2][j][k-1]);
            }
            
        }
        
        // update the points next to doundary
        // x axix
        if (rank%p==0)
        {
            for (k=0;k<localN;k++)
            {
                u[(i+1)%2][1][k]=u[i%2][1][k]+nu*dt/dx*(dt/dx)*(u[i%2][2][k]-2*u[i%2][1][k]+u[i%2][0][k])+nu*dt/dx*(dt/dx)*(u[i%2][1][k+1]-2*u[i%2][1][k]+u[i%2][1][k-1]);
            }
        }
        if (rank%p==p-1)
        {
            for (k=0;k<localN;k++)
            {
                u[(i+1)%2][localN-2][k]=u[i%2][localN-2][k]+nu*dt/dx*(dt/dx)*(u[i%2][localN-1][k]-2*u[i%2][localN-2][k]+u[i%2][localN-3][k])+nu*dt/dx*(dt/dx)*(u[i%2][localN-2][k+1]-2*u[i%2][localN-2][k]+u[i%2][localN-2][k-1]);
            }
        }
        // boundary on y axix
        if (rank/p==0)
        {
            for (j=0;j<localN;j++)
            {
                k=1;
                u[(i+1)%2][1][k]=u[i%2][1][k]+nu*dt/dx*(dt/dx)*(u[i%2][2][k]-2*u[i%2][1][k]+u[i%2][0][k])+nu*dt/dx*(dt/dx)*(u[i%2][1][k+1]-2*u[i%2][1][k]+u[i%2][1][k-1]);
            }
        }
        if (rank/p==p-1)
        {
            for (j=0;j<localN;j++)
            {
                k=localN-2;
                u[(i+1)%2][1][k]=u[i%2][1][k]+nu*dt/dx*(dt/dx)*(u[i%2][2][k]-2*u[i%2][1][k]+u[i%2][0][k])+nu*dt/dx*(dt/dx)*(u[i%2][1][k+1]-2*u[i%2][1][k]+u[i%2][1][k-1]);
            }
        }
        
        // gather the result together
        if ((i+1)%4000==0)||(i==0)// print out and write out the output
        {
            if (rank==0)
            {
                double* finalu=(double*)malloc((N*N)*sizeof(double));
                double final_want[N][N];
            }
            for (k=0;k<p;k++)
            {
                for (j=1;j<localN-1;j++)
                {
                    MPI_Gather(&(u[i%2][j+k*p][1]),localN-2,MPI_DOUBLE,&(final_want[j+k*p-1][0]),localN-2,MPI_DOUBLE,p*(k+1),MPI_COMM_WORLD);
                }
            }
             MPI_Gather(final_want,N,MPI_DOUBLE,finalu,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
            fwrite(finalu, sizeof(double), N*N, outfile);
        }
        
    }
    free (finalu);
    free (x);
    free (y);
    free (new_x);

    fclose(outfile);
    
    
    
    MPI_Finalize();
    return 0;
}
