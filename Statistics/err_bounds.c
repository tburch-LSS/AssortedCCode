#include <math.h>
#include <stdio.h>

main( int argc, char **argv ){
  int i,n=0,area;
  double avg,step,data[10000],x,min=1000.0,max=0.0,low,high;

  argc = 3;
  sscanf(argv[1],"%le",&avg);
  sscanf(argv[2],"%le",&step);

  while ( scanf("%le",&data[n]) == 1 ) {
    if ( data[n] < min ) min = data[n];
    if ( data[n] > max ) max = data[n];
    n++;
  }

  area = 0;
  x = min;
  while ( area < 0.15866*n ) {
    for ( i=0; i<n; i++ ) if ( data[i] <= x && data[i] > (x-step) ) area += 1;
    x += step;
  }
  low = x;

  area = 0;
  x = min;
  while ( area < 0.84134*n ) {
    for ( i=0; i<n; i++ ) if ( data[i] <= x && data[i] > (x-step) ) area += 1;
    x += step;
  }
  high = x;

  area = 0;
  x = max;
  while ( area < 0.84134*n ) {
    for ( i=0; i<n; i++ ) if ( data[i] >= x && data[i] < (x+step) ) area += 1;
    x -= step;
  }
  low += x;
  low *= 0.5;

  area = 0;
  x = max;
  while ( area < 0.15866*n ) {
    for ( i=0; i<n; i++ ) if ( data[i] >= x && data[i] < (x+step) ) area += 1;
    x -= step;
  }
  high += x;
  high *= 0.5;

  printf("Avg = %le + %le - %le\n",avg,(high-avg),(avg-low));
}
