/*
  logistic.c

  generates points (mu,x_fp) of the logistic map

  3/2000 tb
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define N_PTS 50
#define TOL 1.e-4


main(){
  double x,x_fp[N_PTS],mu;
  int i,j,check;

  for( mu=0.25; mu<=1.0; mu+=0.0002 ){
    for( i=0, x=0.1; i<300; i++ ) x = 4.0*mu*x*(1.0 - x);
    for( i=0; i<N_PTS; i++ ){
      x = 4.0*mu*x*(1.0 - x);
      x_fp[i] = x;
      for( j=0, check=0; j<i; j++ ) if( fabs(x-x_fp[j]) < TOL ) check++;
      if( check == 0 ) printf("%e %e\n",mu,x);
    }
  }
}
