/*
  cubic.c

  given the 4 real coefficients of a cubic polynomial, find the three 
  (in general, complex) roots
*/

#include <stdio.h>
#include <math.h>
#include "cubic.h"


double discriminant( double *C ) 
{
  double D,bc,ad;

  bc = C[2]*C[1];
  ad = C[3]*C[0];
  D = 4.0*C[2]*C[2]*C[2]*C[0] - bc*bc;
  D += 4.0*C[3]*C[1]*C[1]*C[1] - 18.0*ad*bc;
  D += 27.0*ad*ad;

  return D;
}


/* find COMPLEX roots of C[3]*x^3 + C[2]*x^2 + C[1]*x + C[0] = 0 
   and store solutions in x[i] */
void cubic_solve ( dcomplex *x, double *C ) 
{
  int i,j;
  double discrim,a,b,c,p,q,t2,r3,theta3,r,theta,u2,tmp;
  dcomplex u[3];

  /*  discrim = discriminant(C); */

  /* Cardano's method */

  /* first, make cubic coefficient = 1 */
  a = 1.0/C[3];
  b = a*C[1];
  c = a*C[0];
  a *= C[2];

  /* roots: 
     x = p/u - u - a/3
     first from principal cube root of 
     u = { q +- [ q^2 + p^3 ]^1/2 }^1/3 */
  p = b*THIRD - a*a*NINTH;
  q = 0.5*c + a*(a*a-4.5*b)*TWENSEVTH;
  t2 = q*q + p*p*p;
  if( t2 >= 0.0 ) {
    u[0].imag = 0.0;
    r3 = q + sqrt(t2);  /* + or - works here */
    if( r3 > 0.0 ) u[0].real = exp(THIRD*log(r3));
    else if( r3 < 0.0 ) u[0].real = -exp(THIRD*log(-r3));
    else {
      x[0].real = x[1].real = x[2].real = - THIRD*a;
      x[0].imag = x[1].imag = x[2].imag = 0.0;
      return;
    }
    u2 = u[0].real*u[0].real;
    tmp = p/u2;
    x[0].real = tmp*u[0].real - u[0].real - THIRD*a;
    x[0].imag = 0.0;
  }
  else {
    r3 = sqrt( q*q - t2 );
    theta3 = atan2(sqrt(-t2),q);
    r = exp(THIRD*log(r3));
    theta = THIRD*theta3;
    u[0].real = r*cos(theta);
    u[0].imag = r*sin(theta);
    u2 = r*r;
    tmp = p/u2;
    x[0].real = x[0].imag = tmp;
    x[0].real *= u[0].real; x[0].imag *= -u[0].imag;
    x[0].real -= u[0].real; x[0].imag -= u[0].imag;
    x[0].real -= THIRD*a;
  }

  /* other 2 roots: u_+- = u0 * [ -1/2 +- i sqrt(3)/2 ] */
  u[1].real = -0.5*u[0].real - HALFRADTHREE*u[0].imag;
  u[1].imag = -0.5*u[0].imag + HALFRADTHREE*u[0].real;
  u[2].real = -0.5*u[0].real + HALFRADTHREE*u[0].imag;
  u[2].imag = -0.5*u[0].imag - HALFRADTHREE*u[0].real;

  x[1].real = x[1].imag = tmp;
  x[1].real *= u[1].real; x[1].imag *= -u[1].imag;
  x[1].real -= (u[1].real + THIRD*a);
  x[1].imag -= u[1].imag;

  x[2].real = x[2].imag = tmp;
  x[2].real *= u[2].real; x[2].imag *= -u[2].imag;
  x[2].real -= (u[2].real + THIRD*a);
  x[2].imag -= u[2].imag;

  return;

} /* cubic_solve */


/* find REAL roots of C[3]*x^3 + C[2]*x^2 + C[1]*x + C[0] = 0 
   and store solutions in x[i] */
void cubic_solve_real ( double *x, double *C ) 
{
  int i,j;
  double a,b,c,p,q,t2,r3,theta3,r,theta,u2,ur,tmp,eps=1.e-6;
  dcomplex u0;

  /* Cardano's method */

  /* first, make cubic coefficient = 1 */
  a = 1.0/C[3];
  b = a*C[1];
  c = a*C[0];
  a *= C[2];

  /* roots: 
     x = p/u - u - a/3 
     first from principal cube root of 
     u = { q +- [ q^2 + p^3 ]^1/2 }^1/3 = r e^(i theta) */
  p = b*THIRD - a*a*NINTH;
  q = 0.5*c + a*(a*a-4.5*b)*TWENSEVTH;
  t2 = q*q + p*p*p;
  if( t2 >= 0.0 ) {
    u0.imag = 0.0;
    r3 = q + sqrt(t2);  /* + or - works here */
    if( r3 > 0.0 ) u0.real = exp(THIRD*log(r3));
    else if( r3 < 0.0 ) u0.real = -exp(THIRD*log(-r3));
    else {
      x[0] = x[1] = x[2] = - THIRD * a;
      return;
    }
    u2 = u0.real*u0.real;
    tmp = p/u2;
  }
  else {
    r3 = sqrt( q*q - t2 );
    theta3 = atan2(sqrt(-t2),q);
    r = exp(THIRD*log(r3));
    theta = THIRD*theta3;
    u0.real = r*cos(theta);
    u0.imag = r*sin(theta);
    u2 = r*r;
    tmp = p/u2;
  }

  if( (u0.imag != 0.0) && (fabs(tmp+1.0) > eps) ) {
    printf("Complex prinicpal root in cubic_real_solve!!!\n");
    exit(1);
  }
  x[0] = x[1] = x[2] = - THIRD * a;
  x[0] += (tmp*u0.real - u0.real);

  /* other 2 roots: u_+- = u0 * [ -1/2 +- i sqrt(3)/2 ] */
  ur = -0.5*u0.real - HALFRADTHREE*u0.imag;
  x[1] += (tmp*ur - ur);
  ur = -0.5*u0.real + HALFRADTHREE*u0.imag;
  x[2] += (tmp*ur - ur);

  return;

} /* cubic_solve_real */


/* testing code */
/*
int main() 
{
  double C[4],xr[3],dd;
  dcomplex x[3];

  printf("a x^3 + b x^2 + c x + d = 0\n");
  printf("enter a:");
  scanf("%le",&C[3]);
  printf("enter b:");
  scanf("%le",&C[2]);
  printf("enter c:");
  scanf("%le",&C[1]);
  printf("enter d:");
  scanf("%le",&C[0]);

  cubic_solve( x, C );

  printf("x0 = %g + %g i\n",x[0].real,x[0].imag);
  printf("x1 = %g + %g i\n",x[1].real,x[1].imag);
  printf("x2 = %g + %g i\n",x[2].real,x[2].imag);

  dd = discriminant(C);
  if( dd < 0.0 ) {
    printf("discriminant < 0 ; roots are real:\n");

    cubic_solve_real( xr, C );

    printf("x0 = %g\n",xr[0]);
    printf("x1 = %g\n",xr[1]);
    printf("x2 = %g\n",xr[2]);
  }

  return 0;
}
*/
