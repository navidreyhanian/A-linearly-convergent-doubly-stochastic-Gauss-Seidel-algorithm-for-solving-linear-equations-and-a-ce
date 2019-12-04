#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define HEIGHT 1200
#define WIDTH  1
#define PI         3.14159265   // The value of pi
double norm(double mean, double std_dev);  // Returns a normal rv
double rand_val(int seed);                 // Jain's RNG
clock_t start, end;
double cpu_time_used;



// C program for implementation of Bubble sort
#include <stdio.h>



#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int FindInArray(double *SortedArray, int ArraySize, double query);


int main()
{
    FILE *fPointer;
    fPointer = fopen("it_ds.txt","w");
    FILE *gPointer;
    gPointer = fopen("time_ds.txt","w");


    float epsilon=pow(10,(-10));
    int m = 40;
    int n = 30;

    float alpha;
    alpha = .9*pow(n,-1); // not correct for under-determined system
    printf("alpha=%f\n",alpha);
    printf("epsilon %f\n",epsilon);
    int T = 1000000000;
    int i, j;
    double **A;
    A=malloc(sizeof(double*)*m);
    for (i = 0; i < m; i++)
    {
        A[i]=malloc(sizeof(double)*n);
    }
    double **AA;
    AA=malloc(sizeof(double*)*m*n);
    for (i = 0; i < m*n; i++)
    {
        AA[i]=malloc(sizeof(double));
    }
    FILE *myfile;
    double myvariable;

    printf("test1\n");

    myfile=fopen("Mymatrix.txt", "r");

  for(i = 0; i < HEIGHT; i++)
  {
    for (j = 0 ; j < WIDTH; j++)
    {
      fscanf(myfile,"%lf",&myvariable);
      AA[i][j]=myvariable;
    }
    printf("\n");
  }
  fclose(myfile);

  printf("test2\n");

  for(i = 0; i < HEIGHT; i++)
  {
    for (j = 0 ; j < WIDTH; j++)
    {
      printf("%.15f ",AA[i][j]);
    }
    printf("\n");
  }
  for(i = 0; i < m; i++)
  {
    for (j = 0 ; j < n; j++)
    {
      A[i][j]=AA[(i)*n+j][0];
    }
  }
    printf("test3\n");

  srand(time(0));
    double error;
    double rd1=0;
    double rd2=0;
    double x[n][1];
    double x_true[n][1];
    double b[m][1];
    int k;
    int t;
    int ell;
    double sum0=0;
    float loc;
    double sum1=0;
    double mul;
    double mul1;
    double err=0;
    double add1=0;
    int add2;
    int hh=0;
    double rh[m][1];
    float lambda=1.5;
    double nrm=0;
    int idxx = 0;
    int beta = 1;
    double **res;
    res=malloc(sizeof(double*)*beta);
    for (i = 0; i < beta; i++)
    {
        res[i]=malloc(sizeof(double));
    }
    printf("test4\n");

    int count=0;
    int in=0;
    double idx;
    int add;
    double **tau;
    tau=malloc(sizeof(double*)*m);
    for (i = 0; i < m; i++)
    {
        tau[i]=malloc(sizeof(double)*n);
    }
    int bin;
    double **P_double;
    P_double=malloc(sizeof(double*)*m);
    for (i = 0; i < m; i++)
    {
        P_double[i]=malloc(sizeof(double)*n);
    }
    double **A1;
    A1=malloc(sizeof(double*)*m);
    for (i = 0; i < m; i++)
    {
        A1[i]=malloc(sizeof(double)*n);
    }
    int ind=0;
    for ( i = 0; i < n; i++ )
    {

        for ( j = 0; j < 1; j++ )
        {
            x[i][j] = (double)rand() / (double)RAND_MAX;
        }
    }
    for ( i = 0; i < n; i++ )
    {

        for ( j = 0; j < 1; j++ )
        {
            x_true[i][j] = x[i][j];
        }
    }
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < 1; j++)
        {
            b[i][j] = 0;
            for (k = 0; k < n; k++)
                b[i][j] += A[i][k]*x[k][j];
        }
    }
    for ( i = 0; i < n; i++ )
    {

        for ( j = 0; j < 1; j++ )
        {
            x[i][j] = (double)rand() / (double)RAND_MAX;
        }
    }
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            P_double[i][j]=1;
        }
    }
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            P_double[i][j]=A[i][j]*A[i][j];
        }
    }
    add1=0;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            add1 = add1 + P_double[i][j];
        }
    }


     for (i = 0; i < m; i++)
     {
        for (j = 0; j < n; j++)
        {
            P_double[i][j]=P_double[i][j]/add1;
        }
    }



    double *SumP;  
    SumP = (double*)calloc(n*m, sizeof(double));

    double temp = 0;
    for (i = 0; i<m; i++)
    {    
        for (j = 0; j<n; j++)
        {
            SumP [j + n*i] = temp;
            temp  = temp  + P_double[i][j];  
        }
    }


/*
    int Idx,IdxI,IdxJ;

int counter = 0;
    for (i=0; i <m; i++)
    {
        for (j=0; j < n; j++)
        {
//            printf("%.3f  ", A[i][j]);
     //       printf("%.3f  ", SumP[i+m*j]);
            printf("%.3f\n", SumP[counter]);
            counter++;
        }
        printf("\n");
    }
    for (t = 0; t< 100; t++)
    {
        loc = (double)rand() / (double)RAND_MAX;
        printf("%.3f  ", loc);
        Idx = FindInArray(SumP, m*n, loc);
        IdxI = Idx/n;
        IdxJ = Idx%n;
        printf("%d %d\n", IdxI+1, IdxJ+1);
    }

    */
    int Idx;

    for ( t = 0; t < T; t++ )
    {
        for ( ell = 0; ell < n; ell++)
        {
            loc = (double)rand() / (double)RAND_MAX;
            Idx = FindInArray(SumP, m*n, loc);
            i = Idx/n;
            j = Idx%n;
            mul=0;
            //printf("i=%d,j=%d \n", i,j);
            for (k = 0; k < n; k++)
            {
                mul +=  - A[i][k]*x[k][0];
            }
            mul += b[i][0];
            //printf("mul = %f\n", mul );
            x[j][0]= x[j][0] + alpha * mul / A[i][j];

        }
        err=0;
        if (t%100 == 0)
        {
            for (i=0 ; i <m ; i++)
            {
                mul=0;
                for (k = 0; k < n; k++)
                {
                    mul +=  - A[i][k]*x[k][0];
                }
                mul += b[i][0];
                err += pow(mul,2);
            }
            error=sqrt(err);
            printf("error = %f\n", log(error) );
            //printf("t = %d\n", t);
            if (error <= epsilon){break;}
        }
    }
    fprintf(fPointer,"%d\n", t);
    printf("t = %d\n", t);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("cpu_time_used = %f\n", cpu_time_used);
    fprintf(gPointer,"%f\n", cpu_time_used);



return 0;
}


int FindInArray(double *SortedArray, int ArraySize, double query)
{
    int i = 0;
    int j = ArraySize ;
    int mid;
    while(j-i>1.1)
    {
        mid = (i+j)/2;
        if (query < SortedArray[mid])
            j = mid;
        if (query >= SortedArray[mid])
            i = mid;
    }
    return i;
}
















