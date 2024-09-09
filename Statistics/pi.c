/*
  pi.c

  compute pi using monte carlo

  to compile:   cc -o pi pi.c -lm
  to execute:   pi

  10/2000 tb
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double f( double x );


int main(){
  double x,y,pi,avg_pi,std_err,var,error;
  int i,j,seed,npoints,ntrials,nbelow;

  /* get seed value for random number generator */
  fprintf(stderr,"Enter a seed\n");
  scanf("%d",&seed);
  srand48(seed);

  /* get # points and trials */
  fprintf(stderr,"Enter # points\n");
  scanf("%d",&npoints);
  fprintf(stderr,"Enter # trials\n");
  scanf("%d",&ntrials);

  avg_pi = 0.0;
  var = 0.0;

  /* loop over trial # */
  for( i=0; i<ntrials; i++ ){

    /* loop over points */
    for( j=0, nbelow=0; j<npoints; j++ ){

      /* generate random point in the box (0-1,0-1) */
      x = drand48();
      y = drand48();

      /* is the point within the unit circle? */
      if( y < f(x) ) nbelow++;
    }

    /* calculate pi estimate for this trial */
    pi = 4.0*(double)nbelow/(double)npoints;

    avg_pi += pi/(double)ntrials;
    var += pi*pi/(double)(ntrials - 1);
  }
  var -= avg_pi*avg_pi*(double)ntrials/(double)(ntrials - 1);

  /* standard error of the mean */
  std_err = sqrt(var/(double)ntrials);

  /* actual deviation from pi */
  error = fabs(M_PI - avg_pi);

  printf("pi = %g +- %g\nactual error = %g\n",avg_pi,std_err,error);

  return 0;
}


double f( double x ){
  return sqrt( 1.0 - x*x );
}
