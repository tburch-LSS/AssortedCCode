/*
  functions and integrals relating to relatively degenerate fermionic 
  matter

  fermion chemical potential: 
  eta = mu / m c^2

  tb 10/2010
*/

#include<math.h>
#include<stdio.h>
#include "calc.h"
#include "EoS_defines.h"
#define DEG_FERM_FUNC 
#include "deg_ferm_func.h"


double integrand( double x ){
  /*
    x = p / m c  ;  beta = m c^2 / k T  ;  eta = mu / m c^2  ; 
    y - 1 = sqrt(1 + x^2) - 1  (relativistic Kinetic Energy) ; 
    f = 1/(exp(beta(y-1-eta)) + 1)  (relativistic Fermi-Dirac distribution) ; 
    integrand = x^xP * y^yP * exp^eP(beta(y-eta)) f^fP(beta,eta,x)
  */
  double y,e,f,ii;

  y = sqrt(1.0 + x*x);
  e = exp(Beta*(y - 1.0 - Eta));
  f = 1.0 / (e + 1.0);

  if( TP_check == 1 ) ii = - pow(x,xP) * log(e * f);
  else ii = pow(x,xP) * pow(y,yP) * pow(e,eP) * pow(f,fP);

  return(ii);
}


/* for obtaining the number density */
double Integral_1( double beta, double eta ){
  double I1;

  Beta = beta; Eta = eta;
  xP = 2.0; yP = 0.0;
  eP = 0.0; fP = 1.0;

  I1 = x_Deg*x_Deg*x_Deg/3.0;
  /*  I1 += gauss_legendre(x_Deg,x_End,32); */
  I1 += romberg_fast(x_Deg,x_End,10000,2);
  I1 *= INV_PI2;
  return(I1);
}


double dintegrand1_dx_zero( double x ){
  /*
    function to zero in order to find maximum of integrand1: 
    g(beta,eta,x) = 1 + exp(beta(1+eta-y)) - beta*x^2/2y = 0 
    where  y = sqrt(1 + x^2)
  */
  double y,e,g;

  y = sqrt(1.0 + x*x);
  e = exp(Beta*(1.0 + Eta - y));
  g = 1.0 + e - 0.5*Beta*x*x/y;
  return(g);
}


void x_range_from_integrand1( double beta, double eta ){
  int n=0;
  double xmid,xlo,xhi,gmid,glo,ghi;

  Beta = beta; Eta = eta;

  /* degenerate case */
  if( beta*eta > 10.0 ){
    xmid = sqrt(eta*(eta+2.0));
    xlo = xmid - 15.0*(eta+1.0)/(beta*xmid);
    xhi = xmid + 15.0*(eta+1.0)/(beta*xmid);
    if( xlo < 0.0 ) xlo = 0.0;
    /* set upper global bound */
    x_End = xhi;
  }
  else{
    xlo = 0.0;
    xhi = 15.0*(1.0+sqrt(beta))/beta;
  }

  /* set lower global bound */
  x_Deg = xlo;

  /* bisection method to find maximum of integrand1 */
  glo = dintegrand1_dx_zero(xlo);
  ghi = dintegrand1_dx_zero(xhi);
  while( glo*ghi > 0.0 ){
    xlo = xhi;
    glo = ghi;
    xhi *= 2.0;
    ghi = dintegrand1_dx_zero(xhi);
  }
  while( (xhi-xlo) > 1.e-6 ){
    n++;
    xmid = 0.5*(xlo+xhi);
    gmid = dintegrand1_dx_zero(xmid);
    if( gmid*glo <= 0.0 ){
      xhi = xmid;
      ghi = gmid;
    }
    else{
      xlo = xmid;
      glo = gmid;
    }
  }
  xmid = 0.5*(xlo+xhi);

  /* set x(integrand1_max) */
  x_Max = xmid;

  /* if non-degenerate, reset upper global bound */
  if( beta*eta <= 10.0 ) x_End = 6.0*xmid;

  fprintf(stderr,"x(integrand1_max) = %g \t n = %d\n",xmid,n);
  fprintf(stderr,"x_Deg = %g \t x_End = %g\n",x_Deg,x_End);
}


/* for obtaining the pressure */
double Integral_2( double beta, double eta ){
  double I2;

  Beta = beta; Eta = eta;
  xP = 4.0; yP = -1.0;
  eP = 0.0; fP = 1.0;

  I2 = 3.0*phi(x_Deg);
  /*  I2 += gauss_legendre(x_Deg,x_End,32); */
  I2 += romberg_fast(x_Deg,x_End,10000,2);
  I2 *= INV_PI2 / 3.0;
  return(I2);
}


/* for obtaining the energy density */
double Integral_3( double beta, double eta ){
  double I3;

  Beta = beta; Eta = eta;
  xP = 2.0; yP = 1.0;
  eP = 0.0; fP = 1.0;

  I3 = chi(x_Deg);
  /*  I3 += gauss_legendre(x_Deg,x_End,32); */
  I3 += romberg_fast(x_Deg,x_End,10000,2);
  I3 *= INV_PI2;
  return(I3);
}


/* for obtaining the pressure in the completely degenerate case */
double phi( double x ){
  double y,ff;

  y = sqrt(1.0 + x*x);
  ff = x * y * (2.0*x*x/3.0 - 1.0);
  ff += log(x + y);
  ff *= 0.125;
  return(ff);
}


/* for obtaining the energy density in the completely degenerate case */
double chi( double x ){
  double y,gg;

  y = sqrt(1.0 + x*x);
  gg = x * y * (2.0*x*x + 1.0);
  gg -= log(x + y);
  gg *= 0.125;
  return(gg);
}


/* electrostatic correction for degenerate electron pressure and 
   energy density */
void elec_stat_corr( double I1, int Z, double *I2_corr, double *I3_corr ){
  double prefact,corr;

  /* Chandrasekhar's ee & eZ interactions */
/*  prefact = -0.3 * ( K_e / ( INV_Le3 * MeC2 ) ) * 
    pow( 1.33333333333333 * PI * (double)(Z*Z) , 0.33333333333333 );  */

  /* BCC lattice */
  prefact = -0.48141 * (K_e / ( INV_Le3 * MeC2 )) * 
      pow( (double)Z , 0.66666666666667 );

  corr = pow( I1 * INV_Le3 , 1.33333333333333 );
  *I2_corr = prefact * corr;
  *I3_corr = 3.0 * prefact * corr;
}


/* used for calculating  d^n Integral / d beta^n-m d eta^m  */
double dI_dbeta_deta( double beta, double eta, 
		      double xp, double yp, double ep, double fp ){
  double dI;

  Beta = beta; Eta = eta;
  xP = xp; yP = yp;
  eP = ep; fP = fp;

  if( eP >= fP ) fprintf(stderr,"Warning: eP >= fP in dI_dbeta_deta\n");

  dI = 0.0;
  /*  dI += gauss_legendre(x_Deg,x_End,32); */
  dI += romberg_fast(x_Deg,x_End,10000,2);
  dI *= INV_PI2;
  return(dI);
}


/* "inversion" of Fermi-Dirac equation of state: 
   from beta and I1(beta,eta) (i.e., from temperature and number density), 
   determine eta (i.e., chemical potential) */
double eta_from_beta_I1( double beta, double I1 ){
  double deta,eta,I1i,dIde;

  deta = 1.0;
  /* first guess from degenerate case:  I1 = x^3 / 3 pi^2 ;  eta ~= x  */
  eta = pow( 3.0*PI*PI*I1, 0.33333333333333 );
  x_range_from_integrand1(beta,eta);
  while( fabs(deta) > 1.e-6 ){
    I1i = Integral_1(beta,eta);
    dIde = beta * dI_dbeta_deta(beta,eta,2,0,1,2);
    deta = ( I1 - I1i ) / dIde;
    fprintf(stderr,"eta = %g  +  %g  =  %g\n",eta,deta,eta+deta);
    eta += deta;
    x_range_from_integrand1(beta,eta);
  }
  return(eta);
}


double Thermo_Potent_over_MC2( double beta, double eta ){
  double asinh,ThPot;

  TP_check = 1;
  Beta = beta; Eta = eta;
  xP = 2.0;

  asinh = log( x_Deg + sqrt(1.0 + x_Deg*x_Deg) );
//  if( x_Deg > ?? ){
//    ThPot = beta * (eta+1.0) * x_Deg*x_Deg*x_Deg / 3.0;
//    ThPot -= beta * x_Deg*x_Deg*x_Deg*x_Deg / 4.0;
//    ThPot += beta * log(2.0*x_Deg) / 8.0;
//    ThPot += romberg_fast(x_Deg,x_End,10000,2);
//  }
  ThPot = beta * (eta+1.0) * x_Deg*x_Deg*x_Deg / 3.0;
  ThPot -= beta * sinh(4.0*asinh) / 32.0;
  ThPot += beta * asinh / 8.0;
  ThPot += romberg_fast(x_Deg,x_End,10000,2);
  ThPot *= beta * 4.0 * PI * PI;

  TP_check = 0;
  return(ThPot);
}



/* needed for calc.c, but not used at the moment */
double integrand3( double x, double y, double z ){
  return(0);
}
double integrand2( double x, double y ){
  return(0);
}
double function( double x ){
  return(0);
}
