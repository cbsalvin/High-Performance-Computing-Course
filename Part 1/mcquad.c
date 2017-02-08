//
//  main.c
//  Problem_1
//
//  Created by Alex CONG on 11/4/2016.
//  Copyright © 2016 Alex CONG. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
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
    if (argc != 2) {
        printf("Incorrect usage: only enter the input data file name\n");
        return 0;
    }
    FILE* inputfile = fopen(argv[1], "r");
   
    if (!inputfile) {
        printf("Unable to open input file\n");
        return 0;
    }
    // start reading input data using function fscanf here
    
    int N;
    
    fscanf(inputfile, "%d", &N); // read an integer N for example
    // read rest of parameters here
    fclose(inputfile);
    // initial the data
    int i=0;
    double y,x[N],temp;
    double estimate_m, estimate_var;
    clock_t start = clock();
    
    // ... rest of program ...
    // choose different seed for random number generate
    srand48(time(NULL));
    for (i=0;i<N;i++)
    {
        y=drand48();
        x[i]=f_x(-log(y));
    }
    estimate_m=mcquad(x,N);// monte-carlo
    // get the variance
    temp=0.0;
    for (i=0;i<N;i++)
    {
        temp=temp+(x[i]-estimate_m)*(x[i]-estimate_m);
    }
    estimate_var=temp/(N-1);
    printf("Time elapsed: %g seconds\n", (float)(clock()-start)/CLOCKS_PER_SEC);
    printf ("Mean is %f\n",estimate_m);
    printf ("Variance is %f\n",estimate_var);
    return 0;
}
