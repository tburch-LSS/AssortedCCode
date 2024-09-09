/*
  multi-D interpolation routines 

  tb 11/2018 
*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "multiD_interp.h"


/* based upon fractional location (fracx) of abscissa vector x within local hypercuboid and 
   the corner values of z, perform multi-D linear interpolation for z */
double Dlin_interpol( double *fracx, double *zcorn, int dim ){
  int c,d,a;
  double z,fcorn,*one_min_f;

  if( dim > maxD ){
    fprintf(stderr,"Error in Dlin_interpol: maximum number of dimensions limited to %d\n",maxD);
    exit(1);
  }

  /* especially beneficial for higher-D */
  one_min_f = (double *)malloc(dim*sizeof(double));
  for( d=0 ; d<dim ; d++ ) one_min_f[d] = 1.0 - fracx[d];

  // number of corners per dimension is now statically set in Ncorn[] array 
  //  Ncorn = 2;
  //  for( d=1 ; d<dim ; d++ ) Ncorn *= 2;

  z = 0.0;
  a = Astart[dim];

  for( c=0 ; c<Ncorn[dim] ; c++ ){

    /* overall fractional contribution of current corner: 
       each corner's weight is given by the fractional volume OPPOSITE said corner */
    fcorn = 1.0;
    for( d=0 ; d<dim ; d++ , a++ ){
      if( Above[a] ) fcorn *= fracx[d];
      else fcorn *= one_min_f[d];
      //      else fcorn *= (1.0 - fracx[d]);
    }

    /* add corner's contribution to overall interpolated value */
    z += fcorn * zcorn[c];
  }

  free(one_min_f);
  return(z);
}


/* based upon fractional location (fracx) of abscissa vector x within local hypercuboid and 
   the corner values of z, perform inverse-distance weighted interpolation for z */
double InvFracDistHcub_interpol( double *fracx , double *zcorn , double p , int dim ){
  int c,d,a;
  double z,dist2,invdp,tinvdp,*fx2,*omfx2;

  if( dim > maxD ){
    fprintf(stderr,"Error in InvFracDistHcub_interpol: maximum number of dimensions limited to %d\n",maxD);
    exit(1);
  }

  /* especially beneficial for higher-D */
  fx2 = (double *)malloc(dim*sizeof(double));
  omfx2 = (double *)malloc(dim*sizeof(double));
  for( d=0 ; d<dim ; d++ ) fx2[d] = fracx[d] * fracx[d];
  for( d=0 ; d<dim ; d++ ) omfx2[d] = (1.0 - fracx[d]) * (1.0 - fracx[d]);

  z = 0.0;
  a = Astart[dim];
  tinvdp = 0.0;

  for( c=0 ; c<Ncorn[dim] ; c++ ){

    /* (fractional-) distance^2 between current corner and fracx */
    dist2 = 0.0;
    for( d=0 ; d<dim ; d++ , a++ ){
      if( Above[a] ) dist2 += omfx2[d];
      else dist2 += fx2[d];
    }

    /* current weight = 1/(distance)^p */
    if( dist2 > 0.0 ){
      invdp = pow( dist2 , -0.5*p );
      tinvdp += invdp;
    }
    else return( zcorn[c] );

    /* add corner's contribution to overall interpolated value */
    z += invdp * zcorn[c];
  }

  free(fx2);
  free(omfx2);
  return( z / tinvdp );
}


/* based upon abscissa vector x and tabular values of z, find the appropriate corner values of 
   the local z hypercuboid and the fractional location of x within this h-cuboid */
void set_corners_fracs( double *zcorn , double *fracx, 
			double *ztab , int *ntab , double *x , double *xstart , double *dx, int dim ){
  int d,d2,c0,c,nc,a,step[dim],index[dim];
  double xdiff;

  if( dim > maxD ){
    fprintf(stderr,"Error in set_corners_fracs: maximum number of dimensions limited to %d\n",maxD);
    exit(1);
  }

  /* find index vector for first corner of hypercuboid and fractional location of x within */
  for( d=0 ; d<dim ; d++ ){
    xdiff = x[d] - xstart[d];
    fracx[d] = xdiff/dx[d];
    index[d] = (int)( fracx[d] );
    fracx[d] -= (double)( index[d] );
    //    fracx[d] = (xdiff - index[d]*dx[d]) / dx[d];
    if( index[d] >= ntab[d] || index[d] < 0 ){
      fprintf(stderr,"Error in set_corners_fracs: x[%d] is out of bounds of table\n",d);
      exit(1);
    }
  }

  /* find first corner index for z table, along with size of steps needed for advance in each dimension */
  c0 = 0;
  for( d=0 ; d<dim ; d++ ){
    step[d] = 1;
    for( d2=0 ; d2<d ; d2++ ) step[d] *= ntab[d2];
    c0 += step[d] * index[d];
  }

  /* using steps, find and set corner z values */
  zcorn[0] = ztab[c0];
  a = Astart[dim] + dim;
  for( nc=1 ; nc<Ncorn[dim] ; nc++ ){
    c = c0;
    for( d=0 ; d<dim ; d++,a++ ) if( Above[a] ) c += step[d];
    zcorn[nc] = ztab[c];
  }

  return;
}


/* based upon location (deltax) of abscissa vector x within local hypercuboid and 
   the corner values of z, perform inverse-distance weighted interpolation for z */
double InvDistHcub_interpol( double *deltax , double *zcorn , double *dx , double p , int dim ){
  int c,d,a;
  double z,dist2,invdp,tinvdp,*dx2,*dmdx2;

  if( dim > maxD ){
    fprintf(stderr,"Error in InvDistHcub_interpol: maximum number of dimensions limited to %d\n",maxD);
    exit(1);
  }

  /* especially beneficial for higher-D */
  dx2 = (double *)malloc(dim*sizeof(double));
  dmdx2 = (double *)malloc(dim*sizeof(double));
  for( d=0 ; d<dim ; d++ ) dx2[d] = deltax[d] * deltax[d];
  for( d=0 ; d<dim ; d++ ) dmdx2[d] = (dx[d] - deltax[d]) * (dx[d] - deltax[d]);

  z = 0.0;
  a = Astart[dim];
  tinvdp = 0.0;

  for( c=0 ; c<Ncorn[dim] ; c++ ){

    /* distance^2 between current corner and x */
    dist2 = 0.0;
    for( d=0 ; d<dim ; d++ , a++ ){
      if( Above[a] ) dist2 += dmdx2[d];
      else dist2 += dx2[d];
    }

    /* current weight = 1/(distance)^p */
    if( dist2 > 0.0 ){
      invdp = pow( dist2 , -0.5*p );
      tinvdp += invdp;
    }
    else return( zcorn[c] );

    /* add corner's contribution to overall interpolated value */
    z += invdp * zcorn[c];
  }

  free(dx2);
  free(dmdx2);
  return( z / tinvdp );
}


/* based upon abscissa vector x and tabular values of z, find the appropriate corner values of 
   the local z hypercuboid and the location of x within this h-cuboid */
void set_corners_deltas( double *zcorn , double *deltax, 
			 double *ztab , int *ntab , double *x , double *xstart , double *dx, int dim ){
  int d,d2,c0,c,nc,a,step[dim],index[dim];
  double xdiff;

  if( dim > maxD ){
    fprintf(stderr,"Error in set_corners_deltas: maximum number of dimensions limited to %d\n",maxD);
    exit(1);
  }

  /* find index vector for first corner of hypercuboid and location of x within */
  for( d=0 ; d<dim ; d++ ){
    xdiff = x[d] - xstart[d];
    index[d] = (int)( xdiff/dx[d] );
    deltax[d] = xdiff - index[d]*dx[d];
    if( index[d] >= ntab[d] || index[d] < 0 ){
      fprintf(stderr,"Error in set_corners_deltas: x[%d] is out of bounds of table\n",d);
      exit(1);
    }
  }

  /* find first corner index for z table, along with size of steps needed for advance in each dimension */
  c0 = 0;
  for( d=0 ; d<dim ; d++ ){
    step[d] = 1;
    for( d2=0 ; d2<d ; d2++ ) step[d] *= ntab[d2];
    c0 += step[d] * index[d];
  }

  /* using steps, find and set corner z values */
  zcorn[0] = ztab[c0];
  a = Astart[dim] + dim;
  for( nc=1 ; nc<Ncorn[dim] ; nc++ ){
    c = c0;
    for( d=0 ; d<dim ; d++,a++ ) if( Above[a] ) c += step[d];
    zcorn[nc] = ztab[c];
  }

  return;
}


/* based upon location of abscissa vector x, perform inverse-distance weighted interpolation for z using ALL points in regular table */
double InvDistAll_interpol( double *ztab , int *ntab , double *x , double *xstart , double *dx , double p , int dim ){
  int iz,itot,d,d2,step[dim],index;
  double z,dist2,invdp,tinvdp,deltax,*xdiff;

  /* especially beneficial for higher-D */
  xdiff = (double *)malloc(dim*sizeof(double));
  for( d=0 ; d<dim ; d++ ) xdiff[d] = x[d] - xstart[d];

  /* find size of steps needed in ztab for advance in each dimension */
  for( d=0 , itot=1 ; d<dim ; d++ ){
    step[d] = 1;
    for( d2=0 ; d2<d ; d2++ ) step[d] *= ntab[d2];
    itot *= ntab[d];
  }

  z = 0.0;
  tinvdp = 0.0;

  /* loop over all values in z-table */
  for( iz=0 ; iz<itot ; iz++ ){

    /* distance^2 between current x-table-value and x */
    dist2 = 0.0;
    for( d=0 ; d<dim ; d++ ){
      /* find index vector component for current iz value */
      index = ((int) (iz / step[d])) % ntab[d];
      deltax = xdiff[d] - index * dx[d];
      dist2 += deltax * deltax;
    }

    /* current weight = 1/(distance)^p */
    if( dist2 > 0.0 ){
      invdp = pow( dist2 , -0.5*p );
      tinvdp += invdp;
    }
    else return( ztab[iz] );

    /* add corner's contribution to overall interpolated value */
    z += invdp * ztab[iz];
  }

  free(xdiff);
  return( z / tinvdp );
}


double RBFunction( double r , double *params ){
  double rbf;

  //  rbf = exp(-params[0] * r);

  /* inverse distance weighting with p=params[0] */
  rbf = pow( r , -params[0] );

  return( rbf );
}


/* based upon location of abscissa vector x, perform radial-basis-function weighted interpolation for z using ALL points of scattered data */
double RadBasFuncScat_interpol( double *ztab , double *xtab , int ntot , double *params , double *x , int dim ){
  int iz,ix,d;
  double z,dist2,r,rbf,totw,deltax;

  z = 0.0;
  totw = 0.0;

  /* loop over all values in z-table */
  for( iz=0 , ix=0 ; iz<ntot ; iz++ ){

    /* distance^2 between current x-table-value and x */
    dist2 = 0.0;
    for( d=0 ; d<dim ; d++ , ix++ ){
      deltax = x[d] - xtab[ix];
      dist2 += deltax * deltax;
    }
    r = sqrt( dist2 );

    /* current weight = RBFunc(r) */
    if( dist2 > 0.0 ) rbf = RBFunction( r , params );
    else return( ztab[iz] );
    totw += rbf;

    /* add corner's contribution to overall interpolated value */
    z += rbf * ztab[iz];
  }

  return( z / totw );
}

