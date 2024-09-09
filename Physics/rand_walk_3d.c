/*
  rand_walk.c

  3D random walk with unit steps

  to compile:   cc -o rand_walk_3d rand_walk_3d.c -lm
  to execute:   rand_walk_3d > rand_walk_3d.dat

  10/1999 tb
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


main(){
  int min,max,nsteps,nwalks,ntrials,i,j,k,seed;
  double x,y,z,r2,theta,phi,sum_r2,sum_r4,avg_r2,SUM_r2,SUM_r4,AVG_r2,
    std_err;

  /* get seed value to start random number generator */
  fprintf(stderr,"Enter a seed\n");
  scanf("%d",&seed);
  srand48(seed);

  /* get min and max # steps, # walks, and # trials */
  fprintf(stderr,"Enter min # steps\n");
  scanf("%d",&min);
  fprintf(stderr,"Enter max # steps\n");
  scanf("%d",&max);
  fprintf(stderr,"Enter # walks\n");
  scanf("%d",&nwalks);
  fprintf(stderr,"Enter # trials\n");
  scanf("%d",&ntrials);

  /* loop over # steps */
  for ( nsteps=min; nsteps<=max; nsteps++ ){
    SUM_r2 = 0.0;
    SUM_r4 = 0.0;

    /* loop over trial # */
    for ( i=0; i<ntrials; i++ ){
      sum_r2 = 0.0;
      sum_r4 = 0.0;

      /* loop over walk # */
      for ( j=0; j<nwalks; j++ ){
	x = 0.0;
	y = 0.0;
	z = 0.0;

	/* take steps */
	for ( k=0; k<nsteps; k++ ){
	  theta = acos(2.0*drand48() - 1.0);
	  phi = 2.0*M_PI*drand48();
	  x += sin(theta)*cos(phi);
	  y += sin(theta)*sin(phi);
	  z += cos(theta);
	}

	/* calculate sums for R^2 */
	r2 = x*x + y*y + z*z;
	sum_r2 += r2;
	sum_r4 += r2*r2;
      }

      /* calculate sums for <R^2> and std error */
      avg_r2 = sum_r2/(double)nwalks;
      SUM_r2 += avg_r2;
      SUM_r4 += avg_r2*avg_r2;
    }

    AVG_r2 = SUM_r2/(double)ntrials;
    std_err = sqrt((SUM_r4 - ((double)ntrials)*AVG_r2*AVG_r2)/
		   (double)(ntrials*ntrials - ntrials));
    printf("%d %le %le\n",nsteps,AVG_r2,std_err);
  }
}
