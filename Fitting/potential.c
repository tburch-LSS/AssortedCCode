/* for use with "fit.c".
* fit to form  V(r) = b[0] - b[1]/r + b[2]*r
*
* compile with cofit.c min.c linalg.c
*/
#include <stdio.h>
#include <math.h>
#include "min.h"

extern int npar;

void f_init(FILE *fp)
{
  fprintf(stderr,"static-quark potential fitting\n");
  npar = 3;
}

double f(double x, double *b)
/* b[0] = V_0
*  b[1] = alpha
*  b[2] = (string tension)^2
*/
{
  double z;

  z  = b[0] - b[1]/x + b[2]*x;
  return(z);
}

double df(double x, double *b, int i)
{
  double z;

  if( i==0 ) z = 1.0;
  if( i==1 ) z = -1.0/x;
  if( i==2 ) z = x;
  return(z);
}

double ddf(double x, double *b, int i, int j)
{
  return(0.0);
}
