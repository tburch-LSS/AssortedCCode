/*
  bilinear interpolation 

  tb 02/2011
*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "bilin_interp.h"


double bilin_interpol( double fracx1, double fracy1, 
		       double z11, double z12, double z21, double z22 ){
  double z;

/*
  fracx1 = (x2 - x) / (x2 - x1);
  fracx2 = (x - x1) / (x2 - x1);  = 1.0 - fracx1;
  fracy1 = (y2 - y) / (y2 - y1);
  fracy2 = (y - y1) / (y2 - y1);  = 1.0 - fracy1;
*/
  z = z11 * fracx1 * fracy1 + z12 * fracx1 * (1.0 - fracy1) + 
      z21 * (1.0 - fracx1) * fracy1 + z22 * (1.0 - fracx1) * (1.0 - fracy1);
  return(z);
}


double inv_bilin_interpol( double x1, double x2, 
			   double fracy1, double z, 
			   double z11, double z12, double z21, double z22 ){
  double f,den,num,x;

  f = z11 * fracy1;
  den = -f;
  num = f * x2;

  f = z12 * (1.0 - fracy1);
  den += -f;
  num += f * x2;

  f = z21 * fracy1;
  den += f;
  num += -f * x1;

  f = z22 * (1.0 - fracy1);
  den += f;
  num += -f * x1;

  x = ( z * (x2 - x1) - num ) / den;
  return(x);
}

