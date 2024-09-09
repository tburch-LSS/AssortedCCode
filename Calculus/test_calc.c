/*
  test_calc.c

  tests integration and differentiation routines

  to compile:   cc -o test_calc test_calc.c calc.c -lm
  to execute:   test_calc

  3/2000 tb
*/

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include "calc.h"

main(){
  double a,b,c,d,h,x,result,exact,error,pi,aa[3],bb[3];
  int n,m,k,nxyz[3];

  a = 1.e-12;
  b = 35.0;
  pi = acos(-1.0);

  exact = pi*pi*pi*pi/15.0;   /* integral of x^3/(exp(x)-1) from 0 to Inf */

  for( n=20; n<=100; n+=20 ){
    result = trapezoid(a,b,n);
    error = fabs( result - exact );
    printf("trap result (n=%d) = %24.22f\nerror = %24.22f\n",n,result,error);
  }

  k = 2;
  for( n=5; n<=25; n+=5 ){
    result = romberg(a,b,n,k);
    error = fabs( result - exact );
    printf("romberg (k=%d, n=%d) result = %24.22f\nerror = %24.22f\n",k,n,
	   result,error);
  }

  n = 32;
  result = gauss_legendre(a,b,n);
  error = fabs( result - exact );
  printf("gauss_legendre(%d) result = %24.22f\nerror = %24.22f\n",n,result,
	 error);

  a = 0.0;
  n = 18;
  result = gauss_laguerre(a,n);
  error = fabs( result - exact );
  printf("gauss_laguerre(%d) result = %24.22f\nerror = %24.22f\n\n",n,
	 result,error);

  x = 1.0;
  exact = -exp(-x);   /* derivative of exp(-x) at x=1 */

  for( h=1.0; h>0.2; h/=2.0 ){
    result = cent_diff(x,h);
    error = fabs( result - exact );
    printf("cent_diff result (h=%g) = %24.22f\nerror = %24.22f\n",h,result,
	   error);
  }

  k = 1;
  for( h=1.0; h>0.2; h/=2.0 ){
    result = richardson(x,h,k);
    error = fabs( result - exact );
    printf("richardson (k=%d) result (h=%g) = %24.22f\nerror = %24.22f\n",k,
	   h,result,error);
  }

  a = 0.0;
  b = 1.0;
  c = 1.0;
  d = 2.0;
  n = 200;
  m = 200;
  exact = 7.0*(exp(1.0)-1.0)/3.0;   /* integral of y^2*exp(x) over ([0,1],[1,2]) */
  result = trap_2d(a,b,c,d,n,m);
  error = fabs( result - exact )/exact;
  printf("\n2d trap result = %24.22f\nrel. error = %24.22f\n",result,error);

  aa[0] = 0.0;
  aa[1] = aa[2] = 1.0;
  bb[0] = 1.0;
  bb[1] = bb[2] = 2.0;
  nxyz[0] = nxyz[1] = nxyz[2] = 200;
  exact = 35.0*(exp(1.0)-1.0)/4.0;   /* integral of y^2*z^3*exp(x) over ([0,1],[1,2],[1,2]) */
  result = trap_3d(aa,bb,nxyz);
  error = fabs( result - exact )/exact;
  printf("\n3d trap result = %24.22f\nrel. error = %24.22f\n",result,error);
}


double function( double x ){
  return exp(-x);
}


double integrand( double x ){
  return x*x*x/(exp(x)-1.0);
}


double integrand2( double x, double y ){
  return y*y*exp(x);
}


double integrand3( double x, double y, double z ){
  return y*y*z*z*z*exp(x);
}
