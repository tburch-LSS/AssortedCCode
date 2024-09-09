/*
  test_interp.c 

  tests for interpolation routines in multiD_interp.c 

  11/2018 tb 
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "multiD_interp.h"

double func( double x );


int main(){
  double x[maxD],dx[maxD],xstart[maxD],fracx[maxD],deltax[maxD],params[maxD];
  double zexact,zint,p,*zcorn,*ztab,*xtab;
  int seed,i,d,d2,dim,iz,ix,index,ntab[maxD],ntot[maxD],step[maxD];

  /* set grid here */
  for( d=0 ; d<maxD ; d++ ) ntab[d] = 11;
  for( d=0 ; d<maxD ; d++ ) dx[d] = 0.1;
  for( d=0 ; d<maxD ; d++ ) xstart[d] = 1.0;

  for( d=0 ; d<maxD ; d++ ) ntot[d] = 1;
  for( d=0 ; d<maxD ; d++ ) for( d2=d ; d2<maxD ; d2++ ) ntot[d2] *= ntab[d];
  for( d=0 ; d<maxD ; d++ ){ step[d] = 1; for( d2=0 ; d2<d ; d2++ ) step[d] *= ntab[d2];}

  for( d=0 ; d<maxD ; d++ ) params[d] = 1.0;

  printf("Enter a seed \n");
  scanf("%d",&seed);
  srand48(seed);

  for( d=0 ; d<maxD ; d++ ) x[d] = xstart[d] + drand48();

  /* set dimension of test here */
  dim = 3;

  /* allocate memory for arrays */
  zcorn = (double *)malloc(Ncorn[dim]*sizeof(double));
  ztab = (double *)malloc(ntot[dim-1]*sizeof(double));
  xtab = (double *)malloc(dim*ntot[dim-1]*sizeof(double));

  /* exact result */
  zexact = 1.0;
  for( d=0 ; d<dim ; d++ ){
    zexact *= func(x[d]);
    printf("x[%2d] = %g   ",d,x[d]);
  }
  printf("\n Exact result = %g\n",zexact);

  /* set up tables for z and x's */
  for( iz=0 , ix=0 ; iz<ntot[dim-1] ; iz++ ){
    for( d=0 ; d<dim ; d++ , ix++ ){
      /* find index vector component for current iz value */
      index = ((int) (iz / step[d])) % ntab[d];
      xtab[ix] = xstart[d] + index * dx[d];
    }
    ix -= dim;
    ztab[iz] = 1.0;
    for( d=0 ; d<dim ; d++ , ix++ ) ztab[iz] *= func(xtab[ix]);
  }

  set_corners_fracs(zcorn,fracx,ztab,ntab,x,xstart,dx,dim);

  zint = Dlin_interpol(fracx,zcorn,dim);
  printf("Linearly interpolated result = %g\n",zint);

  p = params[0];
  zint = InvFracDistHcub_interpol(fracx,zcorn,p,dim);
  printf("InvFracDistHcub result = %g\n",zint);

  set_corners_deltas(zcorn,deltax,ztab,ntab,x,xstart,dx,dim);

  zint = InvDistHcub_interpol(deltax,zcorn,dx,p,dim);
  printf("InvDistHcub result = %g\n",zint);

  zint = InvDistAll_interpol(ztab,ntab,x,xstart,dx,p,dim);
  printf("InvDistAll result = %g\n",zint);

  zint = RadBasFuncScat_interpol(ztab,xtab,ntot[dim-1],params,x,dim);
  printf("RadBasFuncScat result = %g\n",zint);

  free(ztab);
  free(xtab);
  free(zcorn);

  return(0);
}


double func( double x ){
  //  return( exp(x) );
  return( log(x) );
}

