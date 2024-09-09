/* *** error.c *** */
/* to compile: "cc -o error error.c -lm" */
/* to run: error [blocksize] < [datafile (single column)]" */
/* output: #measurements, blocksize, #blocks, average +- std. error */

#include <math.h>
#include <stdio.h>

#define Mostblocks 100000

main( int argc, char **argv ) {
  int nmeas,blocksize,nblocks,iblock;
  double *avg,grandavg,temp,stderror;

  avg = (double *)malloc(Mostblocks*sizeof(double));

  argc = 2;
  sscanf(argv[1],"%d",&blocksize);

  for ( iblock=0; iblock<Mostblocks; iblock++ ) avg[iblock] = 0.0;
  grandavg = 0.0;
  nmeas = nblocks = 0;
  while ( scanf("%le",&temp) == 1 ) {
    nmeas += 1;
    avg[nblocks] += temp;
    if ( nmeas%blocksize == 0 ) {
      avg[nblocks] /= blocksize;
      grandavg += avg[nblocks];
      nblocks += 1;
    }
  }
  grandavg /= nblocks;

  stderror = 0.0;
  for ( iblock=0; iblock<nblocks; iblock++ ) stderror += (avg[iblock] -
	                                grandavg)*(avg[iblock]-grandavg);
  stderror = sqrt(stderror/(1.0*(nblocks*nblocks-nblocks)));

  printf("\nNmeas = %d\nBlocksize = %d\nNblocks = %d\nAvg = %e +- %e\n",
	 nmeas,blocksize,nblocks,grandavg,stderror);
}
