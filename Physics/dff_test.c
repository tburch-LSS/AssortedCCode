/*
  testing code for...
  functions and integrals relating to highly degenerate fermionic matter

  tb 10/2010
*/

#include<math.h>
#include<stdio.h>
#define MAIN 
#include "calc.h"
#include "EoS_defines.h"
#include "deg_ferm_func.h"


int main( int argc, char **argv ){
  int n,N=1.e4,Z;
  double dx,x,xh,beta,eta,I1,I2,I3,I2c,I3c,dIdeta,deta;

  argc = 7;
  sscanf(argv[1],"%le",&xP);
  sscanf(argv[2],"%le",&yP);
  sscanf(argv[3],"%le",&eP);
  sscanf(argv[4],"%le",&fP);
  sscanf(argv[5],"%le",&Beta);
  sscanf(argv[6],"%le",&Eta);

  xh = sqrt(Eta*(Eta+2.0));
  x_range_from_integrand1(Beta,Eta);

/* x_Deg=0.0; */

TP_check = 1;

  dx = (x_End-x_Deg) / (double)N;
  for( n=0; n<N; n++ ){
    x = x_Deg+n*dx;
    printf("%le\t%le\n",x,integrand(x));
  }
/*
  beta = Beta; Beta = 10000.0;
  printf("#cm 1\n");
  for( n=0; n<N; n++ ){
    x = x_Deg+n*dx;
    printf("%le\t%le\n",x,integrand(x));
  }
  Beta = beta;
*/
  fflush(stdout);

TP_check = 0;

  fprintf(stderr,"xP = %g\tyP = %g\nbeta = %g\teta = %g\nx_Deg = %g\tx_End = %g\n",
	  xP,yP,Beta,Eta,x_Deg,x_End);
  I1 = Integral_1(Beta,Eta);
  I2 = Integral_2(Beta,Eta);
  I3 = Integral_3(Beta,Eta);
  Z = 6;
  elec_stat_corr(I1,Z,&I2c,&I3c);
  fprintf(stderr,"1/(3(x_half)^3) = %g\nI1(beta,eta) = %g\n",
	  INV_PI2*xh*xh*xh/3.0,I1);
  fprintf(stderr,"phi(x_half) = %g\nI2(beta,eta) = %g\tI2_corr(Z=%2d) = %g\n",
	  INV_PI2*phi(xh),I2,Z,I2c);
  fprintf(stderr,"chi(x_half) = %g\nI3(beta,eta) = %g\tI3_corr(Z=%2d) = %g\n",
	  INV_PI2*chi(xh),I3,Z,I3c);
  fprintf(stderr," n_e = %g cm^-3\tP(n_e,Z) = %g dyn/cm^2\te(n_e,Z) = %g erg/cm^3\n",
	  I1*INV_Le3,(I2+I2c)*MeC2*INV_Le3,(I3+I3c)*MeC2*INV_Le3);

  eta = pow( 3.0*PI*PI*I1, 0.33333333333333 );
  fprintf(stderr,"Initial guess: \n eta = %g\n",eta);
  eta = eta_from_beta_I1(Beta,I1);
  fprintf(stderr,"Successful inversion? \n eta = %.14f\n",eta);
}
