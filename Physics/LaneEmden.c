/*
  numerical solver for the Lane-Emden Equation 
  (1/xi^2) d_( xi^2 d_theta/d_xi )/d_xi + theta^n = 0 

  written as two first-order DEs: 
  d_theta/d_xi = phi / xi^2 
  d_phi/d_xi = -theta^n xi^2 

  tb 11/2018 
*/

#include<math.h>
#include<stdio.h>
#include "ode.h"
#include "phys_consts.h"

/* Polytropic index, n */
double poly_n;


int main( int argc, char **argv ){

  int i=0,imax;
  double dxi,xi;

  /* "thetaphi" array stores present theta value [0], present phi [1], 
     derivs theta' and phi' [2-3], and previous theta and phi [4-5] */
  double thetaphi[6],*dtpdx;
  dtpdx = thetaphi + 2;

  argc = 4;
  sscanf(argv[1],"%le",&poly_n);
  sscanf(argv[2],"%le",&dxi);
  sscanf(argv[3],"%d",&imax);

  xi = 0.001*dxi;
  thetaphi[0] = 1.0;
  thetaphi[1] = 0.0;
  printf("%g %g %g\n",xi,thetaphi[0],thetaphi[1]);

  while( i<imax ){
    thetaphi[4] = thetaphi[0];
    thetaphi[5] = thetaphi[1];
    rk4(thetaphi,dtpdx,xi,dxi,2);
    printf("%g %g %g\n",xi,thetaphi[0],thetaphi[1]);
    i++;
    xi = i*dxi;
  }

  printf("Lane-Emden Solution\n%d steps of RK4:\n",i);

}


void derivs( double *dydx, double x, double *y, int n ){
  double x2=x*x;
  dydx[0] = y[1]/x2;
  dydx[1] = -pow(y[0],poly_n)*x2;
}
