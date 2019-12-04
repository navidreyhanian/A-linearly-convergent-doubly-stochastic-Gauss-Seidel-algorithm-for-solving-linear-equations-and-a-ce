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

int main()
{
	FILE *fPointer;
	fPointer = fopen("it_rk.txt","w");
	FILE *gPointer;
	gPointer = fopen("time_rk.txt","w");


	float epsilon=pow(10,(-10));
	int m = 40;
	int n = 30;
	float alpha;
	alpha = 2*pow(n,-1);
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
	myfile=fopen("Mymatrix.txt", "r");

	for(i = 0; i < HEIGHT; i++)
	{
		for (j = 0 ; j < WIDTH; j++)
		{
			fscanf(myfile,"%lf",&myvariable);
			AA[i][j]=myvariable;
		}
		//printf("\n");
	}
	fclose(myfile);

	for(i = 0; i < HEIGHT; i++)
	{
		for (j = 0 ; j < WIDTH; j++)
		{
		//	printf("%.15f ",AA[i][j]);
		}
		//printf("\n");
	}
	for(i = 0; i < m; i++)
	{
		for (j = 0 ; j < n; j++)
		{
			A[i][j]=AA[(i)*n+j][0];
		}
	}
	double **P_RK;
	P_RK=malloc(sizeof(double*)*m);
	for (i = 0; i < m; i++)
	{
		P_RK[i]=malloc(sizeof(double));
	}
	double **P_RK1;
	P_RK1=malloc(sizeof(double*)*m);
	for (i = 0; i < m; i++)
	{
		P_RK1[i]=malloc(sizeof(double));
	}
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
		for (j = 0; j < 1; j++)
		{
			P_RK[i][j]=1;
		}
	}
	for (i = 0; i < m; i++)
	{
		sum0=0;
		for (j = 0; j < n; j++)
		{
			sum0=sum0+A[i][j]*A[i][j];
		}
		P_RK1[i][0]=sum0;
	}
	sum0=0;
	for (i = 0; i < m; i++)
	{
		sum0=sum0+P_RK1[i][0];
	}
	for (i = 0; i < m; i++)
	{
		P_RK[i][0]=P_RK1[i][0]/sum0;
	}




	for ( t = 0; t < T; t++ )
	{
		loc = (double)rand() / (double)RAND_MAX;
		sum1=0;
		sum0=0;
		for (i = 0; i < m; i++)
		{
			sum1 += P_RK[i][0];
			if (loc < sum1 && sum0 <= loc)
			{
				add1=0;
				for (j = 0; j < n; j++)
				{
					add1=add1+A[i][j]*x[j][0];
				}
					//printf("ii = %d\n", ii );
                    //printf("j = %d\n", j );
				//printf("add1 = %f\n", add1 );
				for (k = 0; k < n; k++)
				{
					x[k][0]= x[k][0] + (b[i][0]-add1)*A[i][k]/P_RK1[i][0];
                        //printf("x = %f\n", x[k][0] );
				}
			}
			sum0 = sum1;
		}


		err=0;
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
		//printf("error = %f\n", log(error) );
		if (error <= epsilon){break;}

	}
	fprintf(fPointer,"%d\n", t);
	printf("t = %d\n", t);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("cpu_time_used = %f\n", cpu_time_used);
    fprintf(gPointer,"%f\n", cpu_time_used);

	return 0;
}


