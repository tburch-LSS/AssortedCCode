/*
  ode.c

  routines for integrating ordinary differential equations;
  requires the function 'derivs' (possibly external) which calculates the
  dydx's using the DE's

  3/2000 tb
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ode.h"

/*
void derivs( double *dydx, double x, double *y, int n ){
  double ;

  dydx[0] = ;
  dydx[1] = ;
}
*/

void rk2( double *y, double *dydx, double x, double dx, int n ){
  /* 2nd order Runge-Kutta algorithm for one step (dx) of integration */
  int i;
  double k,*yold;

  yold = (double *)malloc(n*sizeof(double));

  derivs(dydx,x,y,n);
  for( i=0; i<n; i++ ){
    k = dx*dydx[i];
    yold[i] = y[i];
    y[i] = yold[i] + k/2.0;
  }
  x += 0.5*dx;
  derivs(dydx,x,y,n);
  for( i=0; i<n; i++ ){
    k = dx*dydx[i];
    y[i] = yold[i] + k;
  }

  free(yold);
}


void rk4( double *y, double *dydx, double x, double dx, int n ){
  /* 4th order Runge-Kutta algorithm for one step (dx) of integration */
  int i;
  double k,*yold,*dy;

  yold = (double *)malloc(n*sizeof(double));
  dy = (double *)malloc(n*sizeof(double));

  derivs(dydx,x,y,n);
  for( i=0; i<n; i++ ){
    k = dx*dydx[i];
    yold[i] = y[i];
    y[i] = yold[i] + k/2.0;
    dy[i] = k/6.0;
  }
  x += 0.5*dx;
  derivs(dydx,x,y,n);
  for( i=0; i<n; i++ ){
    k = dx*dydx[i];
    y[i] = yold[i] + k/2.0;
    dy[i] += k/3.0;
  }
  derivs(dydx,x,y,n);
  for( i=0; i<n; i++ ){
    k = dx*dydx[i];
    y[i] = yold[i] + k;
    dy[i] += k/3.0;
  }
  x += 0.5*dx;
  derivs(dydx,x,y,n);
  for( i=0; i<n; i++ ){
    k = dx*dydx[i];
    dy[i] += k/6.0;
    y[i] = yold[i] + dy[i];
  }

  free(yold); free(dy);
}
