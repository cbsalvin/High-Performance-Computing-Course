//
//  main.c
//  Pro_1
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

double f_x(double in)
{
    return cos(in);
}
double mcquad(double a[],int n)
{
    
    double out=0.0;
    int i=0;
    for (i=0;i<n;i++)
    {
        out=out+a[i];
    }
    out=out/n;
    return out;
}
int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    double precision = MPI_Wtick();
    double starttime = MPI_Wtime();
   
   int N_read,T;
    if (rank==0)
    {
        FILE* inputfile = fopen(argv[1], "r");
        
        fscanf(inputfile, "%d", &N_read); // read an integer N for example
        fscanf(inputfile, "%d", &T);// read time period T
        fclose(inputfile);
    }
    MPI_Bcast(&N_read,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&T,1,MPI_INT,0,MPI_COMM_WORLD);
    
    // initial the data
    int i,k,N;
    //int data;
    N=N_read/size;
    double y,x[N],result_each[T],result[T];
    double estimate_m=0.0;
    clock_t start = clock();

    // choose different seed for random number generate
    srand48((unsigned)time(NULL));
    
    
    for (k=0;k<T;k++)
    {
        
        for (i=0;i<N;i++)
        {
            y=drand48();
            x[i]=f_x(-log(y));
        }
        estimate_m=mcquad(x,N);// monte-carlo
        result_each[k]=estimate_m/size;
        //printf ("Mean is %f, k is %d\n",result_each[k],k);
    }

    
    MPI_Reduce(result_each,result,T,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    if (rank==0)
    {
        FILE* outfile=fopen(argv[2], "w");
        fwrite(result, sizeof(double), T, outfile);
        fclose(outfile);
        double elapsed_time = MPI_Wtime()-starttime;
        printf("Execution time = %le seconds, with precision %le seconds\n", elapsed_time, precision);
    }

    
    MPI_Finalize();

    return 0;
}

