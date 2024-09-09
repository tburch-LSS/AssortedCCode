/*
  relax_laplace_2d.c

  solves laplace's equation for the potential inside a box with constant
  boundary conditions using the relaxation method

  to compile:   gcc relax_laplace_2d.c -lm -o relax_laplace_2d.x
  to execute:   relax_laplace_2d.x > [datafile]

  11/2000 tb, altered 10/2006
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NDIV 100
#define FCSTART 30
#define FCSIZE 40
#define FCSKIP 5
#define DPR 2


double dielec( double x, double y );


int main(){
  int i,j,ndiv,check,iter_num=0;
  double **phi,**phi_new,left,right,bottom,top,tol,r,x,y,dx;
  double DE[5],DEavg;

  /* get boundary conditions, # divisions, tolerance, and relaxation
     parameter */
  fprintf(stderr,"potential(x=0,y) = ");
  scanf("%le",&left);
  fprintf(stderr,"potential(x=1,y) = ");
  scanf("%le",&right);
  fprintf(stderr,"potential(x,y=0) = ");
  scanf("%le",&bottom);
  fprintf(stderr,"potential(x,y=1) = ");
  scanf("%le",&top);

  ndiv = NDIV;
  fprintf(stderr,"using %d divisions\n",ndiv);
  /*
  fprintf(stderr,"# divisions = ");
  scanf("%d",&ndiv);
  */
  dx = 1.0/((double)ndiv);

  fprintf(stderr,"tolerance = ");
  scanf("%le",&tol);
  fprintf(stderr,"relaxation parameter (0 < r < 2) = ");
  scanf("%le",&r);

  /* allocate memory for potential (old and updated) */
  phi = (double **)malloc((ndiv+1)*sizeof(double *));
  phi[0] = (double *)malloc((ndiv+1)*(ndiv+1)*sizeof(double));
  for( i=1; i<=ndiv; i++ ) phi[i] = phi[0] + i*(ndiv+1);
  phi_new = (double **)malloc((ndiv+1)*sizeof(double *));
  phi_new[0] = (double *)malloc((ndiv+1)*(ndiv+1)*sizeof(double));
  for( i=1; i<=ndiv; i++ ) phi_new[i] = phi_new[0] + i*(ndiv+1);

  /* set the BC's and initial values */
  for( i=0; i<=ndiv; i++ )
    for( j=0; j<=ndiv; j++ ) {
      if( i == 0 ) phi[i][j] = left;
      else if( i == ndiv ) phi[i][j] = right;
      else if( j == 0 ) phi[i][j] = bottom;
      else if( j == ndiv ) phi[i][j] = top;
      else phi[i][j] = 0.0;
    }

#ifdef DIPOLE
  phi[ndiv/2][ndiv/2-DPR] = 1.0;
  phi[ndiv/2][ndiv/2+DPR] = -1.0;
#endif

#ifdef QPOLE
  phi[ndiv/2-DPR][ndiv/2-DPR] = 1.0;
  phi[ndiv/2-DPR][ndiv/2+DPR] = -1.0;
  phi[ndiv/2+DPR][ndiv/2-DPR] = -1.0;
  phi[ndiv/2+DPR][ndiv/2+DPR] = 1.0;
#endif

  /* relaxation method */
  check = 1;
  while( check > 0 ) {
    for( i=1, check=0; i<ndiv; i++ )
      for( j=1; j<ndiv; j++ ) {

	x=dx*i; y=dx*j;

	/* update potential */
#ifndef DIELEC
	phi_new[i][j] = (1.0-r)*phi[i][j] + 0.25*r*( phi[i-1][j] +
						     phi[i+1][j] +
						     phi[i][j-1] +
						     phi[i][j+1] );
#else
	DE[0] = dielec(x,y);
	DE[1] = dielec(x-dx,y);
	DE[2] = dielec(x+dx,y);
	DE[3] = dielec(x,y-dx);
	DE[4] = dielec(x,y+dx);
	DEavg = 2.0*DE[0] + 0.5*(DE[1]+DE[2]+DE[3]+DE[4]);
	phi_new[i][j] = (1.0-r)*phi[i][j] + 
	  r*(0.5/DEavg)*( (DE[0]+DE[1])*phi[i-1][j] +
			  (DE[0]+DE[2])*phi[i+1][j] +
			  (DE[0]+DE[3])*phi[i][j-1] +
			  (DE[0]+DE[4])*phi[i][j+1] );
#endif

#ifdef FCAGE
	/* check if (i,j) is part of Faraday cage */
	if ( (i==FCSTART) || (i==FCSTART+FCSIZE) ){
	  if ( (j>=FCSTART) && (j<=FCSTART+FCSIZE) && (j%FCSKIP==0) )
	    phi_new[i][j] = phi[i][j];
	}
	if ( (j==FCSTART) || (j==FCSTART+FCSIZE) ){
	  if ( (i>=FCSTART) && (i<=FCSTART+FCSIZE) && (i%FCSKIP==0) )
	    phi_new[i][j] = phi[i][j];
	}
#endif

#ifdef DIPOLE
	if ( i==ndiv/2 ){
	  if ( j==ndiv/2-DPR || j==ndiv/2+DPR ) phi_new[i][j] = phi[i][j];
	}
#endif

#ifdef QPOLE
	if ( i==ndiv/2-DPR || i==ndiv/2+DPR ){
	  if ( j==ndiv/2-DPR || j==ndiv/2+DPR ) phi_new[i][j] = phi[i][j];
	}
#endif

	/* test tolerances */
	if( fabs(phi_new[i][j] - phi[i][j]) > tol ) check++;
      }

    /* old = new */
    for( i=1; i<ndiv; i++ ) for( j=1; j<ndiv; j++ )
      phi[i][j] = phi_new[i][j];
    iter_num++;
  }

  /* output x, y, phi(x,y) */
  for( i=0; i<=ndiv; i++ ) for( j=0; j<=ndiv; j++ )
    printf( "%le %le %le\n", (double)i/(double)ndiv, (double)j/(double)ndiv,
	    phi[i][j] );

  /* output x, phi(x,const) */
  /*
  j = ndiv/2;
  for( i=0; i<=ndiv; i++ )
    printf( "%le %le\n", (double)i/(double)ndiv, phi[i][j] );
  */

  /* output y, phi(const,y) */
  /*
  i = ndiv/2;
  for( j=0; j<=ndiv; j++ )
    printf( "%le %le\n", (double)j/(double)ndiv, phi[i][j] );
  */

  fprintf(stderr,"%d iterations performed.\n",iter_num);

  return 0;
}


double dielec ( double x, double y ){
  double D;

  if ( y <= 0.5 ) D=80.0;
  else D=1.0;
  return D;
}
