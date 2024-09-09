/*
  zode.c 

  routines for integrating complex ordinary differential equations; 
  requires the function 'derivs' (possibly external) which calculates the 
  dydx's using the DE's 

  12/2018 tb 
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "zcomplex.h"
#include "zode.h"

/*
void zderivs( dcomplex *dydx, double x, dcomplex *y, int n ){
  double ;

  dydx[0] = ;
  dydx[1] = ;
}
*/

void zrk2( dcomplex *y, dcomplex *dydx, double x, double dx, int n ){
  /* 2nd order Runge-Kutta algorithm for one step (dx) of integration */
  int i;
  dcomplex k,*yold;

  yold = (dcomplex *)malloc(n*sizeof(dcomplex));

  zderivs(dydx,x,y,n);
  for( i=0; i<n; i++ ){
    ZMULREAL(k,0.5*dx,dydx[i]);
    ZEQ(yold[i],y[i]);
    ZADD(y[i],yold[i],k);
  }
  x += 0.5*dx;
  zderivs(dydx,x,y,n);
  for( i=0; i<n; i++ ){
    ZMULREAL(k,dx,dydx[i]);    
    ZADD(y[i],yold[i],k);
  }

  free(yold);
  return;
}


void zrk4( dcomplex *y, dcomplex *dydx, double x, double dx, int n ){
  /* 4th order Runge-Kutta algorithm for one step (dx) of integration */
  int i;
  dcomplex k,*yold,*dy;

  yold = (dcomplex *)malloc(n*sizeof(dcomplex));
  dy = (dcomplex *)malloc(n*sizeof(dcomplex));

  zderivs(dydx,x,y,n);
  for( i=0; i<n; i++ ){
    ZMULREAL(k,0.5*dx,dydx[i]);
    ZEQ(yold[i],y[i]);
    ZADD(y[i],yold[i],k);
    ZMULREAL(dy[i],0.333333333333,k);
  }
  x += 0.5*dx;
  zderivs(dydx,x,y,n);
  for( i=0; i<n; i++ ){
    ZMULREAL(k,0.5*dx,dydx[i]);
    ZADD(y[i],yold[i],k);
    ZMULEQREAL(k,0.666666666667);
    ZPEQ(dy[i],k);
  }
  zderivs(dydx,x,y,n);
  for( i=0; i<n; i++ ){
    ZMULREAL(k,dx,dydx[i]);
    ZADD(y[i],yold[i],k);
    ZMULEQREAL(k,0.333333333333);
    ZPEQ(dy[i],k);
  }
  x += 0.5*dx;
  zderivs(dydx,x,y,n);
  for( i=0; i<n; i++ ){
    ZMULREAL(k,0.166666666667*dx,dydx[i]);
    ZPEQ(dy[i],k);
    ZADD(y[i],yold[i],dy[i]);
  }

  free(yold); free(dy);
  return;
}

