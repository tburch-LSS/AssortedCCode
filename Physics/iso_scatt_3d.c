/*
  iso_scatt_3d.c

  3D random walk with steps of exponentially distributed lengths

  to compile:   gcc iso_scatt_3d.c -lm -o iso_scatt_3d.x
  to execute:   iso_scatt_3d.x > out.dat

  10/2006 tb
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


int main(){
  int min,max,nwalks,i,j,k,seed;
  double max_rad,tot_dist;
  double x,y,z,dx,dy,dz,dr,r2,theta,phi;

  /* get seed value to start random number generator */
  fprintf(stderr,"Enter a seed\n");
  scanf("%d",&seed);
  srand48(seed);

  /* get max radius, # walks */
  fprintf(stderr,"Enter max radius\n");
  scanf("%le",&max_rad);
  fprintf(stderr,"Enter # walks\n");
  scanf("%d",&nwalks);

  /* loop over walk # */
  for ( j=0; j<nwalks; j++ ){
    x = 0.0; y = 0.0; z = 0.0;
    r2 = 0.0;
    tot_dist = 0.0;

    /* take steps */
    while ( r2 < max_rad*max_rad ){

      /* uniform distribution for d(cos(theta)) */
      theta = acos(2.0*drand48() - 1.0);

      /* uniform distribution in phi */
      phi = 2.0*M_PI*drand48();

      /* exponential distribution (m.f.p.=1) in distance to next scattering */
      dr = -log(drand48());

      dx = dr*sin(theta)*cos(phi);
      dy = dr*sin(theta)*sin(phi);
      dz = dr*cos(theta);
      tot_dist += dr;
      x += dx; y += dy; z += dz;
      r2 = x*x + y*y + z*z;
    }

    /* correct last step for better estimate of boundary crossing */
    tot_dist -= dr;
    x -= dx; y -= dy; z -= dz;
    r2 = x*x + y*y + z*z;
    dr /= 1000.0;
    dx /= 1000.0; dy /= 1000.0; dz /= 1000.0;
    while ( r2 < max_rad*max_rad ){
      tot_dist += dr;
      x += dx; y += dy; z += dz;
      r2 = x*x + y*y + z*z;
    }

    printf("%le\n",tot_dist);

  }

  return 0;
}
