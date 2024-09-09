/*
  zroot.c 

  routines for finding complex roots of a complex function 

  11/2019 tb 
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "zcomplex.h"


dcomplex zFunction( dcomplex z ){
  dcomplex zF;

  // direct F(z): 
  
  return(zF);
}


void decompF( dcomplex *zF , dcomplex z ){
  double x,y,u,v;

  // u(x,y) & v(x,y): 
  x = z.real;
  y = z.imag;
  u = ;
  v = ;
  ZCMPLX(zF[0],u,v);
  return;
}


dcomplex zF_recurs_n( dcomplex z , int n ){
  int j;
  dcomplex zF0,zF1,zFn,ctmp;

  ZCMPLX(zF0,1,0);
  ZEQ(zF1,z);
  ZEQ(zFn,zF1);
  // direct F(z): 
  for( j=1 ; j<n ; j++ ){
    ZEQ(ctmp,zFn);
    // recursion relation: 
    
    // end recursion 
    ZEQ(zF0,zF1);
    ZEQ(zF1,ctmp);
  }
  return(zFn);
}


double sq_mag_F( dcomplex z ){
  double FstarF;
  //  double u,v;
  dcomplex zF;

  zF = zFunction(z);
  //  u = zF.real; v = zF.imag;
  //  FstarF = u*u + v*v;
  FstarF = zmag2(zF);
  return(FstarF);
}

