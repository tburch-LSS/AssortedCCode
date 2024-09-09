/* *** binning.c *** */
/* to compile: "cc -o binning binning.c -lm" */
/* to run: "binning [#bins] < [datafile (single column)]" */
/* output (3 columns): bin-center, frequency, Poissonian error (sqrt(freq)) */

#include <math.h>
#include <stdio.h>

#define Mostmeas 1000000
#define Mostbins 50000

main(int argc,char **argv) {
  int i,N,Nbins,bin,freq[Mostbins];
  double *x,low,high,dx;

  argc = 2;
  sscanf(argv[1],"%d",&Nbins);
  x = (double *)malloc(Mostmeas*sizeof(double));

  N = 0;
  while ( scanf("%le",&x[N]) == 1 ) {
    if ( N == 0 ) {
      low = x[N];
      high = x[N];}
    else if ( x[N] < low ) low = x[N];
    else if ( x[N] > high ) high = x[N];
    N += 1;
  }

  dx = (high-low)/Nbins;
  for ( i=0; i<N; i++ ) {
    bin = (x[i]-low)/dx;
    if ( bin == Nbins ) bin -= 1;
    freq[bin] += 1;
  }

  for ( bin=0; bin<Nbins; bin++ ) printf("%le %i %le\n",(low+(bin+0.5)*dx),
					 freq[bin],sqrt(1.0*freq[bin]));
}
