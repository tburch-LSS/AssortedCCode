#include <math.h>
#include <stdio.h>

#define Mostmeas 1000000

main(int argc, char **argv) {
  int t,tau,N,narg;
  double *x,xavg,xtxT,xt2,autocorr;

  x = (double *)malloc(Mostmeas*sizeof(double));
  argc=2;
  narg=sscanf(argv[1],"%d",&tau);
/*  printf("%i arguments scanned from command line\n",narg); */

  N = 0;
  xavg=0.0;
  while ( scanf("%le",&x[N]) == 1 ){
    xavg += x[N];
    N += 1;
  }
  xavg /= N;

  xtxT=0.0;
  xt2=0.0;
  for(t=0;t<(N-tau);t++) {
    xtxT += (x[t] - xavg)*(x[t+tau] - xavg);
    xt2 += (x[t] - xavg)*(x[t] - xavg);
  }
  autocorr = xtxT/xt2;
  printf("autocorr(%i) = %e\n",tau,autocorr);
}
