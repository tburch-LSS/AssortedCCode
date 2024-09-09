/*
  square.c

  average distance between two randomly chosen points in the unit square
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


int main(){
  double x1,x2,y1,y2,rsq,avg_r,avg_rsq,sum,sum2;
  int seed,i;

  /* get seed value to start random number generator */
  fprintf(stderr,"Enter a seed\n");
  scanf("%d",&seed);
  srand48(seed);

  for( i=0, sum=sum2=0.0; i<1000000; i++ ){
    x1 = drand48();
    x2 = drand48();
    y1 = drand48();
    y2 = drand48();

    rsq = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
    sum += sqrt( rsq );
    sum2 += rsq;
  }

  avg_r = sum/((double)(i-1));
  avg_rsq = sum2/((double)(i-1));
  printf("<r>           = %g\n",avg_r);
  printf("<r^2>         = %g\n",avg_rsq);
  printf("<r^2> - <r>^2 = %g\n",avg_rsq-avg_r*avg_r);
}
