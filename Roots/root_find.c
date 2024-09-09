/*
  root_find.c

  root-finding algorithms (for 'function')

  3/2000 tb
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "root_find.h"


double bisect( double xlow, double xhigh, double prec, int *n ){
  /* bisection method to find root of function(x) within a certain
     precision */
  double xmid,fmid,flow,fhigh;

  *n = 0;
  if( xlow > xhigh ){
    fprintf(stderr,"Error in bisect: x_low > x_high\n");
    exit(1);
  }
  flow = function(xlow);
  fhigh = function(xhigh);
  if( flow*fhigh > 0.0 ){
    fprintf(stderr,"Error in bisect: no root surrounded\n");
    exit(1);
  }
  while( (xhigh-xlow) > prec ){
    (*n)++;
    xmid = 0.5*(xlow + xhigh);
    fmid = function(xmid);
    if ( fmid*flow <= 0.0 ){
      xhigh = xmid;
      fhigh = fmid;
    }
    else{
      xlow = xmid;
      flow = fmid;
    }
  }
  return 0.5*(xlow + xhigh);
}


double newt_raph( double x, double tol, int *n ){
  /* Newton-Raphson technique to find root within a given "y-tolerance"
     [ |function(x)| <= tol ] */
  double func;

  *n = 0;
  func = function(x);
  while( fabs(func) > tol ){
    (*n)++;
    x = x - func/dfdx(x);
    func = function(x);
  }
  return x;
}
