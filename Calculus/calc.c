/*
  calc.c

  integration (of 'integrand') & differentiation (of 'function') routines; 
  'integrand' and 'function' may be defined externally

  3/2000 tb 
  8/2011 (added trad_3d) tb 
*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "calc.h"


double trapezoid( double a, double b, int n ){
  /* integrates the function 'integrand(x)' using the trapezoid rule over n 
     intervals from a to b */
  double x,dx,sum;
  int i;

  dx = (b-a)/((double)n);
  for( i=1, sum=0.0; i<n; i++ ){
    x = a + dx*i;
    sum += integrand(x);
  }
  sum += 0.5*( integrand(a) + integrand(b) );
  return dx*sum;
}


double trap_2d( double a, double b, double c, double d, int nx, int ny ){
  /* integrates the function 'integrand2(x,y)' over the rectangular domain 
     ([a,b],[c,d]) using the trapezoid rule for both dimensions */
  double x,dx,y,dy,sum,xsum;
  int i,j;

  dx = (b-a)/((double)nx);
  dy = (d-c)/((double)ny);
  for( j=1, sum=0.0; j<ny; j++ ){
    y = c + dy*j;
    for( i=1, xsum=0.0; i<nx; i++ ){
      x = a + dx*i;
      xsum += integrand2(x,y);
    }
    xsum += 0.5*( integrand2(a,y) + integrand2(b,y) );
    sum += dx*xsum;
  }
  for( i=1, xsum=0.0; i<nx; i++ ){
    x = a + dx*i;
    xsum += 0.5*( integrand2(x,c) + integrand2(x,d) );
  }
  xsum += 0.25*( integrand2(a,c) + integrand2(b,c) + integrand2(a,d) +
		 integrand2(b,d) );
  sum += dx*xsum;
  return dy*sum;
}


double trap_3d( double *a, double *b, int *nxyz ){
  /* integrates the function 'integrand3(x,y,z)' over the domain 
     ([a0,b0],[a1,b1],[a2,b2]) using the trapezoid rule for all dimensions */
  double x,dx,y,dy,z,dz,sum,xsum,yxsum;
  int i,j,k;

  dx = (b[0]-a[0])/((double)nxyz[0]);
  dy = (b[1]-a[1])/((double)nxyz[1]);
  dz = (b[2]-a[2])/((double)nxyz[2]);

  for( k=1, sum=0.0; k<nxyz[2]; k++ ){
    z = a[2] + dz*k;
    for( j=1, yxsum=0.0; j<nxyz[1]; j++ ){
      y = a[1] + dy*j;
      for( i=1, xsum=0.0; i<nxyz[0]; i++ ){
	x = a[0] + dx*i;
	xsum += integrand3(x,y,z);
      }
      xsum += 0.5*( integrand3(a[0],y,z) + integrand3(b[0],y,z) );
      yxsum += dx*xsum;
    }
    for( i=1, xsum=0.0; i<nxyz[0]; i++ ){
      x = a[0] + dx*i;
      xsum += 0.5*( integrand3(x,a[1],z) + integrand3(x,b[1],z) );
    }
    xsum += 0.25*( integrand3(a[0],a[1],z) + integrand3(a[0],b[1],z) + integrand3(b[0],a[1],z) 
		   + integrand3(b[0],b[1],z) );
    yxsum += dx*xsum;
    sum += dy*yxsum;
  }

  for( j=1, yxsum=0.0; j<nxyz[1]; j++ ){
    y = a[1] + dy*j;
    for( i=1, xsum=0.0; i<nxyz[0]; i++ ){
      x = a[0] + dx*i;
      xsum += 0.5*( integrand3(x,y,a[2]) + integrand3(x,y,b[2]) );
    }
    xsum += 0.25*( integrand3(a[0],y,a[2]) + integrand3(b[0],y,a[2]) + integrand3(a[0],y,b[2]) 
		   + integrand3(b[0],y,b[2]) );
    yxsum += dx*xsum;
  }
  for( i=1, xsum=0.0; i<nxyz[0]; i++ ){
    x = a[0] + dx*i;
    xsum += 0.25*( integrand3(x,a[1],a[2]) + integrand3(x,b[1],a[2]) + integrand3(x,a[1],b[2]) 
		   + integrand3(x,b[1],b[2]) );
  }
  xsum += 0.125*( integrand3(a[0],a[1],a[2]) + integrand3(a[0],a[1],b[2]) + 
		  integrand3(a[0],b[1],a[2]) + integrand3(a[0],b[1],b[2]) + 
		  integrand3(b[0],a[1],a[2]) + integrand3(b[0],a[1],b[2]) + 
		  integrand3(b[0],b[1],a[2]) + integrand3(b[0],b[1],b[2]) );
  yxsum += dx*xsum;
  sum += dy*yxsum;

  return dz*sum;
}


double romberg( double a, double b, int n, int k ){
  /* combines trapezoid answers using n --> (2^k)*n intervals to get 
     better estimate of integral (e.g., k=1 --> Simpson's rule, k=2 --> 
     Bode's rule) */
  int i,j;
  double result[k+1][k+1],temp;

  if( k <= 0 ) fprintf(stderr,"Bad value of k in romberg\n");

  for( i=0; i<=k; i++, n*=2 ) result[i][0] = trapezoid(a,b,n);

  for( j=1, temp=4.0; j<=k; j++, temp*=4.0 ) for( i=j; i<=k; i++ )
    result[i][j] = ( temp*result[i][j-1] - result[i-1][j-1] )/(temp - 1.0);

  return result[k][k];
}


double romberg_fast( double a, double b, int n, int k ){
  /* combines trapezoidal answers using n --> (2^k)*n intervals to get 
     better estimate of integral (e.g., k=1 --> Simpson's rule, k=2 --> 
     Bode's rule); this routine minimizes the number of function calls by 
     performing the trapezoidal integration locally */
  int i,j,m,mm;
  double result[k+1][k+1],temp,x,dx;

  if( k <= 0 ) fprintf(stderr,"Bad value of k in romberg_fast\n");

  for( j=0; j<=k; j++, n*=2 ) for( i=j; i<=k; i++ )
    result[i][j] = 0.0;

  n /= 2;
  dx = (b-a)/(double)n;

  for( m=1; m<n; m++ ){
    x = a + dx*m;
    temp = integrand(x);
    result[k][0] += temp;
    for( i=k-1, mm=2; i>=0; i--, mm*=2 ) if( m%mm == 0 ) result[i][0] += temp;
  }
  temp = 0.5*integrand(a);
  for( i=0; i<=k; i++ ) result[i][0] += temp;
  temp = 0.5*integrand(b);
  for( i=0; i<=k; i++ ) result[i][0] += temp;
  for( i=k, mm=1; i>=0; i--, mm*=2 ) result[i][0] *= dx*mm;

  for( j=1, temp=4.0; j<=k; j++, temp*=4.0 ) for( i=j; i<=k; i++ )
    result[i][j] = ( temp*result[i][j-1] - result[i-1][j-1] )/(temp - 1.0);

  return result[k][k];
}


double gauss_legendre( double a, double b, int n ){
  /* integrates from a to b using n-point Gauss-Legendre integration; 
     requires 'gauss_legendre.[n]' files which contain the weights and 
     abscissas */
  double weight,abscissa,sum,x;
  char filename[20];
  FILE *wx_file;

  sprintf(filename,"gauss_legendre.%02d",n);
  wx_file = fopen(filename,"r");

  sum = 0.0;
  while( fscanf(wx_file,"%le%le",&abscissa,&weight) == 2 ){
    x = 0.5*(b-a)*(abscissa+1.0) + a;
    sum += weight*integrand(x);
  }

  fclose(wx_file);

  return 0.5*(b-a)*sum;
}


double gauss_laguerre( double a, int n ){
  /* integrates from a to Inf using n-point Gauss-Laguerre integration; 
     requires 'gauss_laguerre.[n]' files which contain the weights and 
     abscissas */
  double weight,abscissa,sum,x;
  char filename[20];
  FILE *wx_file;

  sprintf(filename,"gauss_laguerre.%02d",n);
  wx_file = fopen(filename,"r");

  sum = 0.0;
  while( fscanf(wx_file,"%le%le",&abscissa,&weight) == 2 ){
    x = abscissa + a;
    sum += weight*exp(abscissa)*integrand(x);
  }

  fclose(wx_file);

  return sum;
}


double cent_diff( double x, double h ){
  /* centered-difference form for approximating the derivative of 
     'function(x)' */
  return ( function(x+h) - function(x-h) )/(2.0*h);
}


double richardson( double x, double h, int k ){
  /* uses centered-difference scheme with intervals of h --> (2^k)*h to 
     better approximate the derivative */
  double result[k+1][k+1],temp;
  int i,j;

  if( k <= 0 ) fprintf(stderr,"Bad value for k in richardson\n");

  for( i=0; i<=k; i++, h*=2.0 ) result[i][0] = cent_diff(x,h);

  for( j=1, temp=4.0; j<=k; j++, temp*=4.0 ) for( i=0; i<=k-j; i++ )
    result[i][j] = ( temp*result[i][j-1] - result[i+1][j-1] )/(temp - 1.0);

  return result[0][k];
}

