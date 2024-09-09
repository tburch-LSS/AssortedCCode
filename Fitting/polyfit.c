/*
  polyfit.c

  to compile: cc -o polyfit polyfit.c linalg.c -lm
  to execute: polyfit [order] < [datafile]
*/

#include <stdio.h>
#include <math.h>
#include "linalg.h"


main( int argc, char **argv ){
  int order,i,j,n;
  double **M,**Minv,*a,*b,x,y,sy,xrow,xcol,dum,chi2,yfit;

  argc = 2;
  sscanf(argv[1],"%d",&order);

  a = (double *)malloc((order+1)*sizeof(double));
  b = (double *)malloc((order+1)*sizeof(double));

  M = (double **)malloc((order+1)*sizeof(double *));
  M[0] = (double *)malloc((order+1)*(order+1)*sizeof(double));
  for( i=1; i<=order; i++ ) M[i] = M[0] + i*(order+1);

  Minv = (double **)malloc((order+1)*sizeof(double *));
  Minv[0] = (double *)malloc((order+1)*(order+1)*sizeof(double));
  for( i=1; i<=order; i++ ) Minv[i] = Minv[0] + i*(order+1);

  for( i=0; i<=order; i++ ){
    b[i] = 0.0;
    for( j=0; j<=order; j++ ) M[i][j] = 0.0;
  }

  n = 0;
  while( scanf("%le%le%le",&x,&y,&sy) == 3 ){
    for( i=0, xrow=1.0; i<=order; i++, xrow*=x ){
      b[i] += y*xrow/(sy*sy);
      for( j=0, xcol=1.0; j<=order; j++, xcol*=x )
	M[i][j] += xrow*xcol/(sy*sy);
    }
    n++;
  }

  lineq(M[0],b,a,order+1);

  matinv(M[0],Minv[0],order+1);

  for( i=0; i<=order; i++ )
    printf("a[%d] = %g +- %g\n",i,a[i],sqrt(Minv[i][i]));

  chi2 = 0.0;
  rewind(stdin);
  while( scanf("%le%le%le",&x,&y,&sy) == 3 ){
    for( j=0, yfit=0.0, xrow=1.0; j<=order; j++, xrow*=x ) yfit += a[j]*xrow;
    chi2 += (y-yfit)*(y-yfit)/(sy*sy);
  }

  printf("chi^2/dof = %g\n",chi2/(n-order-1));
}
