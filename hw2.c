/*
	CS440 Numerical Computation

	Homework 2: Jacobi SVD

	Allen Wang

	4/8/16
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct singular_value{ // my struct 

	int index;
	double value;
} singular_value;
int comparison (const void * a, const void * b); // compare function
double sign(double val); // returns negative or positive one based on the sign of double
double g_eigen(double* a, int n, int row, int col); //generates the eigenvalue of particular index/position
singular_value* singular_value_e(double* a, int n); // struct pointer 
void jacobi(double* a, int n, double *s, double *u, double *v); // the jacobi algorithm with cosine and sine
void arrange(double* a, int n);
void rotation(double* a, int n, int p, int q, double* u); //rotation when precision is not perfect, used in jacobi

void print(double *a, int m, int n) { 
  int i, k;
  for (i = 0; i < m; i++) {
    for (k = 0; k < n; k++) {
      printf("\t%.10g\t", a[i * m + k]);
    }
    printf("\n");
  }
}
void jacobi(double* a, int n, double *s, double *u, double *v){
	int i = 0;
	int j = 0;
	double ep = 1e-15;
	double comp;
	int count = 1;

	double* T;
	T= malloc (n*n*sizeof(double));
	arrange(T, n);
	arrange(u, n);
	arrange(v, n);

	while (count != 0){
		count = 0;
		for (i= 0; i< (n-1); ++i){
			for (j=i+1; j<n; ++j){
				comp = ep * (sqrt(g_eigen(a, n, i, i) * g_eigen(a, n, j, j))); 
				if (fabs(g_eigen(a, n, i, j)) > comp){ //apply rotation
					rotation(a, n, i, j, T);
					count = count + 1;
				}
			}
		}
	}
	
	singular_value *sigma;
	sigma = singular_value_e(a, n);
	for (i = 0; i < n; ++i){
		for (j = 0; j < n; ++j){
			u[i * n + j] = a[i * n + sigma[j].index];
			v[i * n + j] = T[i * n + sigma[j].index];
		}
		s[i] = sigma[i].value; // set s = the singular value
	}	

	free(T);
	free(sigma);
}

double g_eigen(double* a, int n, int row, int col){
	int i;
	double sum = 0;
	for(i = 0; i < n; ++i){
		sum += a[i * n + row] * a[i * n + col];
	}
	return sum;
}
	
double sign(double val){
	if(val < 0){
		return -1;
	}
	return 1;
}

void arrange(double* a, int n){
	int i;
	int j;
	for (i = 0; i < n; ++i){
		for (j = 0; j < n; ++j){
			if (i != j){
				a[i * n + j] = 0.0;
			}
			else{
				a[i * n + j] = 1.0;
			}
		}
	}
}



void rotation(double *a, int n, int p, int q, double* u){
	int i;
	double a_pp, a_pq, a_qq, t, c, s, temp;
	a_pp = g_eigen(a, n, p, p);
	a_pq = g_eigen(a, n, p, q);
	a_qq = g_eigen(a, n, q, q);
	t = (a_pp - a_qq) / (2 *(a_pq));
	t = sign(t) / (fabs(t) + sqrt(1 + t * t));
	c = 1 / sqrt(1 + t * t);
	s = (c * t);
	for (i = 0; i < n; ++i){
		temp = a[i * n + p];
		a[i * n + p] = s * a[i * n + q] + c * temp;
		a[i * n + q] = c * a[i * n + q] - s * temp;

		temp = u[i * n + p];
		u[i * n + p] = s * u[i * n + q] + c * temp;
		u[i * n + q] = c * u[i * n + q] - s * temp;
	}
}

int comparison (const void * a, const void * b)
{
	if ((*(singular_value*)a).value > (*(singular_value*)b).value){
		return -1;
	}
	if ((*(singular_value*)a).value == (*(singular_value*)b).value){
		return 0;
	}
	if ((*(singular_value*)a).value <  (*(singular_value*)b).value){
		return 1;
	}
	return 0;
}

singular_value* singular_value_e(double* a, int n){
	singular_value* temp;
	temp = malloc(n * sizeof(singular_value));
	int i;
	int j;
	double sum, sing_val;
	for (j = 0; j < n; ++j){
		sum = 0;
		for (i = 0; i < n; ++i){
			sum += (a[i * n + j] * a[i * n + j]);
		}
		sing_val = sqrt(sum);
		temp[j].value = sing_val;
		temp[j].index = j;
		for (i = 0; i < n; ++i){
			a[i * n + j] /= sing_val;
		}
	}
	qsort(temp, n, sizeof(singular_value), comparison);
	return temp;
}	

int main(int argc, char **argv) { 

  int n = 4;
  int i, j;

  double a[n * n];
  double s[n];
  double u[n * n];
  double v[n * n];
  double id[n * n];
  for (i = 0; i < n; i++){
  	for (j = 0; j < n; j++){
  		if (j != i){
  			id[i * n + j] = 0;
  		}
  		else
  			id[i * n + j] = 1;
  	}
  }

  // a[0] = 15;
  // a[1] = 11;
  // a[2] = 20;
  // a[3] = 32;
  // a[4] = 42;
  // a[5] = 32;
  // a[6] = 11;
  // a[7] = 1;
  // a[8] = 5;
  // a[9] = 7;
  // a[10] = 25;
  // a[11] = 81;
  // a[12] = 1;
  // a[13] = 0;
  // a[14] = 32;
  // a[15] = 11;
//test 2
    a[0] = 1;
  a[1] = 2;
  a[2] = 9;
  a[3] = 11;
  a[4] = 8;
  a[5] = 2;
  a[6] = 1;
  a[7] = 8;
  a[8] = 3;
  a[9] = 9;
  a[10] = 0;
  a[11] = -23;
  a[12] = 14;
  a[13] = 7;
  a[14] = 5;
  a[15] = 6;

  printf("\n");
  print(a, n, n);
  printf("\n");

  jacobi(a, n, s, u, v);

  printf("s\n");
  print(s, 1, n);
  printf("\n");
  printf("u\n");
  print(u, n, n);
  printf("\n");
  printf("v\n");
  print(v, n, n);
  printf("\n");

}

