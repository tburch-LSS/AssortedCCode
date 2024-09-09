/*
  test_mat.c

  to compile:   cc -o test_mat test_mat.c linalg.c -lm
  to execute:   test_mat < test_mat_5.dat

  9/2000 tb
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "linalg.h"


int main(){
  double **a,*x,*b,*s;
  int i,j,n=5;

  /* allocate memory for arrays */
  x = (double *)malloc(n*sizeof(double));
  b = (double *)malloc(n*sizeof(double));
  a = (double **)malloc(n*sizeof(double *));
  a[0] = (double *)malloc(n*n*sizeof(double));
  /* and set up pointers for 2d array */
  for( i=1; i<n; i++ ) a[i] = a[0] + i*n;

  /* read in matrix "a" and vector "b" from stdin */
  for( i=0; i<n; i++ ){
    for( j=0, s=a[i]; j<n; j++, s++ ) scanf("%le",s);
    scanf("%le",&b[i]);
  }

  /* solve the linear system for vector "x" */
  lineq(a[0],b,x,n);

  /* print result */
  for( i=0; i<n; i++ ) printf("%g\n",x[i]);

  /* free allocated memory */
  free(x); free(b); free(a[0]); free(a);

  return 0;
}
