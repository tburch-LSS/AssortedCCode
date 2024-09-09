#include <math.h>
#include <stdio.h>

void srand48( int seed );
double drand48( void );

main( int argc, char **argv ){
  int seed,n_samp,n_data,i,j,select;

  argc = 4;
  sscanf(argv[1],"%d",&seed);
  sscanf(argv[2],"%d",&n_samp);
  sscanf(argv[3],"%d",&n_data);
  srand48(seed);

  for ( i=0; i<n_samp; i++ )
    for ( j=0; j<n_data; j++ ) {
      select = n_data*drand48()+1.0;
      printf("%d %d %d\n",i,j,select);
    }
}
