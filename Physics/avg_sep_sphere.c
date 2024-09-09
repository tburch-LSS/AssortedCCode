/*
  avg_sep_sphere.c 

  average distance between two random points within the unit sphere in 1,2,3 dimensions 

  01/2019 tb 
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


int main(){
  int n,npairs,seed;
  double dist,dist2,dsig;
  double r[2],theta[2],phi[2],x[2],y[2],z[2];
  double dx,dy,dz,dr2,dr;

  /* get seed value to start random number generator */
  fprintf(stderr,"Enter a seed\n");
  scanf("%d",&seed);
  srand48(seed);

  npairs = 100000;

  // 1D 
  dist = 0.0;
  dist2 = 0.0;
  /* loop over pairs of points */
  for ( n=0; n<npairs; n++ ){
    x[0] = 2.0 * (drand48() - 0.5);
    x[1] = 2.0 * (drand48() - 0.5);
    dr = fabs(x[1] - x[0]);
    dr2 = dr * dr;
    dist += dr;
    dist2 += dr2;
  }
  dist /= npairs;
  dist2 /= npairs;
  dsig = sqrt((dist2 - dist*dist)/npairs);
  printf("1D: %le +- %le\n",dist,dsig);

  // 2D 
  dist = 0.0;
  dist2 = 0.0;
  /* loop over pairs of points */
  for ( n=0; n<npairs; n++ ){
    r[0] = drand48();
    theta[0] = 2.0 * M_PI * drand48();
    r[1] = drand48();
    theta[1] = 2.0 * M_PI * drand48();
    x[0] = r[0] * cos( theta[0] );
    y[0] = r[0] * sin( theta[0] );
    x[1] = r[1] * cos( theta[1] );
    y[1] = r[1] * sin( theta[1] );
    dx = x[1] - x[0];
    dy = y[1] - y[0];
    dr2 = dx*dx + dy*dy;
    dr = sqrt(dr2);
    dist += dr;
    dist2 += dr2;
  }
  dist /= npairs;
  dist2 /= npairs;
  dsig = sqrt((dist2 - dist*dist)/npairs);
  printf("2D: %le +- %le\n",dist,dsig);

  // 3D 
  dist = 0.0;
  dist2 = 0.0;
  /* loop over pairs of points */
  for ( n=0; n<npairs; n++ ){
    r[0] = drand48();
    theta[0] = acos(2.0*drand48() - 1.0);
    phi[0] = 2.0 * M_PI * drand48();
    r[1] = drand48();
    theta[1] = acos(2.0*drand48() - 1.0);
    phi[1] = 2.0 * M_PI * drand48();
    x[0] = r[0] * cos( theta[0] ) * cos( phi[0] );
    y[0] = r[0] * cos( theta[0] ) * sin( phi[0] );
    z[0] = r[0] * sin( theta[0] );
    x[1] = r[1] * cos( theta[1] ) * cos( phi[1] );
    y[1] = r[1] * cos( theta[1] ) * sin( phi[1] );
    z[1] = r[1] * sin( theta[1] );
    dx = x[1] - x[0];
    dy = y[1] - y[0];
    dz = z[1] - z[0];
    dr2 = dx*dx + dy*dy + dz*dz;
    dr = sqrt(dr2);
    dist += dr;
    dist2 += dr2;
  }
  dist /= npairs;
  dist2 /= npairs;
  dsig = sqrt((dist2 - dist*dist)/npairs);
  printf("3D: %le +- %le\n",dist,dsig);

  return 0;
}
