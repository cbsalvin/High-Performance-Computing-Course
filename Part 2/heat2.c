//
//  main.c
//  problem_3
//
//  Created by Alex CONG on 12/5/2016.
//  Copyright © 2016 Alex CONG. All rights reserved.
//  Generate the matrix of U[y][x]
//  With the boundary condition u(-1,y,t)=y, u(1,y,t)=-y, u(x,-1,t)=x, u(x,1,t)=-x
//  ADI Method to find ut=nu*(uxx+uyy);
//  for loop
//  {
//  Split the U into localU
//  explicit the localU
//  Gather them and transpose
//  Split the U into localU
//  Lapack to solve diagonal matrix
//  explicit the localU
//  Gather them
//  }
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


int main(int argc, char* argv[]){
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    double precision = MPI_Wtick();
    double starttime = MPI_Wtime();
    
    FILE* outfile=fopen(argv[2], "w");
    
    // start reading input data using function fscanf here
    int N;
    double dt;
    double T;
    double nu;
    
    if (rank==0)
    {
        FILE* inputfile = fopen(argv[1], "r");
        fscanf(inputfile, "%d", &N);
        fscanf(inputfile, "%lf", &dt);
        fscanf(inputfile, "%lf", &T);
        fscanf(inputfile, "%lf", &nu);
        fclose(inputfile);
    }
    MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&T,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&nu,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    double x[N];
    double y;
    double temp;
    double pi=M_PI;
    double dx,dtest,lamda;
    
    int i,j,k; // counter
    int localN,p; // local grid points & size P
    int nsteps; // time step
    
    p=size;// number of processors
    
    fwrite(&N, sizeof(int), 1, outfile); // output N
    
    dtest=N;
    dx=2/(dtest-1); // find the distance for two x grid points
    nsteps = round(T/dt); // how many steps for the whole process
    localN=N/p;
    lamda=nu*dt/2/(dx*dx);
    
    double u[N][N]; // u[x axix][y axix]
    
    double u_local[localN][N];
    double uhat_local[localN][N];

    // initial the matrix u
   
    for (j=0;j<N;j++)
        {
            x[j]=-1.0+j*dx;
        }

    for (i=0;i<N;i++)
    {
            for (k=0;k<N;k++)
            {
                y=-1+dx*k;
                u[i][k]=-x[i]*y+cos(11*x[i]*M_PI/2)*sin(8*y*M_PI);
            }
    }

  
    // write the output
    if (rank==0)
    {
        fwrite(x, sizeof(double), N, outfile);
        fwrite(u, sizeof(double), N*N, outfile);
    }
    //// %%% lapack set up starts %%%///
    // set up the lapack solution matrix.
    
    char tr='N';
    int one=1;
    double* ld=(double*)malloc(N*sizeof(double));
    double* d=(double*)malloc(N*sizeof(double));
    double* ud=(double*)malloc(N*sizeof(double));
    double* uud=(double*)malloc(N*sizeof(double));
    int* pivot = malloc(N*sizeof(int));
    int info;
    
    // initial the arrays
    for (i=0;i<N;i++)
    {
        ld[i]=0.;
        d[i]=0.;
        ud[i]=0.;
        uud[i]=0.;
        pivot[i]=0;
    }
    // set up the arrays
    for (i=0;i<N-1;i++)
    {
        ld[i]=-lamda;
        d[i]=1+2*lamda;
        ud[i]=-lamda;
    }
    ld[N-2]=0.;
    ud[0]=0.;
    d[0]=1.0;
    d[N-1]=1.0;
    // upper triagnoal sovle
    dgttrf_(&N,ld, d, ud, uud, pivot, &info);
    
    // %%% lapack set up ends %%%///
    
    
    // main loop
    for (i=0;i<nsteps;i++)
    {
        // step A
        // split the matrix into sub ones
        for (j=0;j<localN;j++)
        {
            for (k=0;k<N;k++)
            {
                u_local[j][k]=u[rank*localN+j][k];
            }
        }
        // calculate the u with time step
        for (j=0;j<localN;j++)
        {
            for (k=1;k<N-1;k++)
            {
                uhat_local[j][k]=u_local[j][k]+lamda*(u_local[j][k+1]-2*u_local[j][k]+u_local[j][k-1]);
            }
            uhat_local[j][0]=u_local[j][0];
            uhat_local[j][N-1]=u_local[j][N-1];
        }

        //step B
        // gather the result together
        MPI_Allgather(uhat_local,localN*N,MPI_DOUBLE,u,localN*N,MPI_DOUBLE,MPI_COMM_WORLD);
        // transpose the matrix
        for (j=0;j<N;j++)
        {
                for (k=0;k<N;k++)
                {
                    temp=u[j][k];
                    u[j][k]=u[k][j];
                    u[k][j]=temp;
                }
        }
        
        // split the matrix into sub ones
        for (j=0;j<localN;j++)
        {
            for (k=0;k<N;k++)
            {
                u_local[j][k]=u[rank*localN+j][k];
            }
        }
        // solve the matrix
        for (j=0;j<localN;j++)
        {
            dgttrs_(&tr, &N, &one, ld, d, ud, uud, pivot, u_local[j], &N, &info);
        }
        // step C
        // calculate the u with time step
        for (j=0;j<localN;j++)
        {
            for (k=1;k<N-1;k++)
            {
                uhat_local[j][k]=u_local[j][k]+lamda*(u_local[j][k+1]-2*u_local[j][k]+u_local[j][k-1]);
            }
            uhat_local[j][0]=u_local[j][0];
            uhat_local[j][N-1]=u_local[j][N-1];
        }
        // step D
        // gather the result together
        MPI_Allgather(uhat_local,localN*N,MPI_DOUBLE,u,localN*N,MPI_DOUBLE,MPI_COMM_WORLD);
        // transpose the matrix
        for (j=0;j<N;j++)
        {
            for (k=0;k<N;k++)
            {
                temp=u[j][k];
                u[j][k]=u[k][j];
                u[k][j]=temp;
            }
        }
        // split the matrix into sub ones
        for (j=0;j<localN;j++)
        {
            for (k=0;k<N;k++)
            {
                u_local[j][k]=u[rank*localN+j][k];
            }
        }
        // solve the matrix
        for (j=0;j<localN;j++)
        {
            dgttrs_(&tr, &N, &one, ld, d, ud, uud, pivot, u_local[j], &N, &info);
        }
        MPI_Allgather(u_local,localN*N,MPI_DOUBLE,u,localN*N,MPI_DOUBLE,MPI_COMM_WORLD);
        // print out and write out the output

        if ((i+1)%(nsteps/4)==0)
        {
            if (rank==0)
            {
                fwrite(u, sizeof(double), N*N, outfile);
            }
        }
    }
    free (d);
    free (ld);
    free (ud);
    free (uud);
    free (pivot);
    fclose(outfile);
    if (rank == 0)
    {
        double elapsed_time = MPI_Wtime()-starttime;
        printf("Execution time = %le seconds, with precision %le seconds\n", elapsed_time, precision);
    }
    MPI_Finalize();
    return 0;
}
