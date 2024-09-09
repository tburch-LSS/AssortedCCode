/*
  basis_fit.c

  general linear least-squares fit of a function:
  f(x) = Sum_(j=0-->order)[ a_j * basis(x,j) ]
  to a given set of data points (x,y,sigma_y)

  to compile:   cc -o basis_fit basis_fit.c linalg.c -lm
  to execute:   basis_fit [order] < [datafile]

  3/2000 tb
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "linalg.h"

double basis( double x, int j );


int main( int argc, char **argv ){
  int order,i,j,n;
  double **M,**Minv,*a,*b,x,y,sy,chi2,yfit;

  /* get order of fit from command line */
  argc = 2;
  sscanf(argv[1],"%d",&order);

  /* allocate memory for arrays */
  a = (double *)malloc((order+1)*sizeof(double));
  b = (double *)malloc((order+1)*sizeof(double));
  M = (double **)malloc((order+1)*sizeof(double *));
  M[0] = (double *)malloc((order+1)*(order+1)*sizeof(double));
  for( i=1; i<=order; i++ ) M[i] = M[0] + i*(order+1);
  Minv = (double **)malloc((order+1)*sizeof(double *));
  Minv[0] = (double *)malloc((order+1)*(order+1)*sizeof(double));
  for( i=1; i<=order; i++ ) Minv[i] = Minv[0] + i*(order+1);

  /* initialize array values */
  for( i=0; i<=order; i++ ){
    b[i] = 0.0;
    for( j=0; j<=order; j++ ) M[i][j] = 0.0;
  }

  /* read in data points from stdin */
  n = 0;
  while( scanf("%le%le%le",&x,&y,&sy) == 3 ){
    for( i=0; i<=order; i++ ){
      b[i] += y*basis(x,i)/(sy*sy);
      for( j=0; j<=order; j++ )
	M[i][j] += basis(x,i)*basis(x,j)/(sy*sy);
    }
    n++;
  }

  /* solve the linear set of equations for the fit parameters "a" */
  lineq(M[0],b,a,order+1);

  /* get the inverse of the matrix "M" for the errors in the parameters */
  matinv(M[0],Minv[0],order+1);

  /* write out the results */
  for( i=0; i<=order; i++ )
    printf("a[%d] = %g +- %g\n",i,a[i],sqrt(Minv[i][i]));

  /* determine the value of chi^2 */
  chi2 = 0.0;
  rewind(stdin);
  while( scanf("%le%le%le",&x,&y,&sy) == 3 ){
    for( j=0, yfit=0.0; j<=order; j++ ) yfit += a[j]*basis(x,j);
    chi2 += (y-yfit)*(y-yfit)/(sy*sy);
  }

  printf("chi^2/dof = %g\n",chi2/(n-order-1));

  free(a); free(b); free(M[0]); free(M); free(Minv[0]); free(Minv);

  return 0;
}


double basis( double x, int j ){
  return pow(x,(double)j);   /* polynomials */
}
