/*
  Thomas-Fermi equation for electrostatic corrections to the 
  (non-relativistic, degenerate) electron equation of state (EoS)

  phi'' = phi^3/2 / x^1/2

  with BC's: phi(0) = 1  ;  phi'(x_0) = phi(x_0) / x_0

  Then  phi'(x_0)^5/2  can be used to calculate the (cell-averaged) 
  electron pressure [see Eq2.4.32 of Shapiro&Teukolsky(1983)]

  tb 10/2010
*/

#include<math.h>
#include<stdio.h>
#include "ode.h"

#define PI 3.14159265358979

/* (m_e c^2) in ergs */
#define MeC2 8.187104385e-7
/* 1/lambda_e^3 in cm^-3 */
#define INV_Le3 1.736603307e31
/* e^2 / (4 pi eps_0) in erg*cm */
#define K_e 2.307077130e-19
/* a_0 (Bohr radius) in cm */
#define a_0 5.291772086e-9


int main( int argc, char **argv ){

  int n=0,Z;
  double dx,x,dpdx52,invx3,mu,ne,Pe;

  /* "phi" array stores present phi value [0] and its first two derivatives 
     [1-2], along with previous value [3] and previous first derivative [4] */
  double phi[5],*dphidx;
  dphidx = phi + 1;

  argc = 3;
  sscanf(argv[1],"%le",&dx);
  sscanf(argv[2],"%le",&dphidx[0]);
  /* phi'(0) > -1.5880710 (the zero-pressure, infinite-radius case) */
  sscanf(argv[3],"%d",&Z);

  mu = pow( 9.0*PI*PI/(128.0*(double)Z), 0.33333333333333 ) * a_0;

  x = 0.01*dx;
  phi[0] = 1.0;

  while( x*dphidx[0] < phi[0] ){
    phi[3] = phi[0];
    dphidx[3] = dphidx[0];
    rk4(phi,dphidx,x,dx,2);
    n++;
    x = n*dx;
  }

  /* back up and "refine" last step to find more precise BC: 
     phi'(x_0)=phi(x_0)/x_0 */
  /*
  x -= dx;
  dx *= 0.01;
  phi[0] = phi[3];
  dphidx[0] = dphidx[3];
  while( x*dphidx[0] < phi[0] ){
    rk4(phi,dphidx,x,dx,2);
    x += dx;
  }
  */

  dpdx52 = dphidx[0]*dphidx[0]*sqrt(dphidx[0]);
  invx3 = 1.0/(x*x*x);
  ne = 3.0*invx3*(double)Z / (4.0*PI*mu*mu*mu);
  Pe = 0.1*K_e*dpdx52*(double)(Z*Z) / (PI*mu*mu*mu*mu);
  printf("Thomas-Fermi Eq. Solution\n%d steps of RK4:\n",n);
  printf("x_0 = %le\nphi'(x_0) = %le\n\n",x,dphidx[0]);
  printf(" n_e ~ 1/x_0^3 = %le\n P(n_e) ~ phi'(x_0)^5/2 = %le\n\n",invx3,dpdx52);
  printf(" n_e = %le cm^-3 \t P(n_e) = %le dyn/cm^2\n\n",ne,Pe);

//  printf("%le\t%le\n",invx3,dpdx52);

}


void derivs( double *dydx, double x, double *y, int n ){
  dydx[0] = y[1];
  dydx[1] = y[0]*sqrt(y[0]/x);
}
