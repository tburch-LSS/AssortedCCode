/*
  binn.c

  bin a one-column list of numbers into a set number of bins between
  minimum and maximum values

  to compile:   cc -o binn binn.c -lm
  to execute:   binn [min] [max] [nbins] < [datafile]

  11/2000 tb
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


int main( int argc, char **argv ){
  int i,n,nbins,*bin_count;
  double min,max,x,dx;

  /* get minimum and maximum values, and # bins, from command line */
  argc = 4;
  sscanf(argv[1],"%le",&min);
  sscanf(argv[2],"%le",&max);
  sscanf(argv[3],"%d",&nbins);

  /* allocate memory for bin counts and initialize */
  bin_count = (int *)malloc(nbins*sizeof(int));
  for( i=0; i<nbins; i++ ) bin_count[i] = 0;

  /* bin width */
  dx = (max - min)/(double)nbins;

  /* read in data and bin */
  n = 0;
  while( scanf("%le",&x) == 1 )
    if( x >= min && x < max) {
      i = (int)(nbins*(x - min)/(max - min));
      bin_count[i]++;
      n++;
    }

  /* output avg bin values, bin counts and associated errors (poisson
     statistics) */
  for( i=0; i<nbins; i++ ) printf( "%le %d %le\n", min+((double)i+0.5)*dx,
				   bin_count[i], sqrt((double)bin_count[i]) );

  fprintf(stderr,"%d numbers binned between %g and %g\n",n,min,max);

  return 0;
}
