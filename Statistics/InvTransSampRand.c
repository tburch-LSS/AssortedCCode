/*
  InvTransSampRand.c 

  generate random numbers with a specific distribution given the 
  Quantile Function (or inverse Cumulative Distribution Function) 

  12/2018 tb 
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "InvTransSampRand.h"

int seed;
bool G2bool;
double G2xval,G2yval;


double RandWithDist( int type , double *params ){
  double x,y;

  x = drand48();
  switch( type ){
  case 0: 
    // exponential dist: params = ( alpha )  ;  exp(-alpha*x) 
    y = -log( x ) / params[0];
    break;
  case 1: 
    // Breit-Wigner dist: params = ( x0 , gamma )  ;  f(x;x0,gamma) = (1/pi) * gamma /[(x-x0)^2 + gamma^2] 
    y = params[0] + params[1] * tan(M_PI*(x - 0.5));
    break;
  case 2: 
    // Gaussian dist: params = ( sigma , mean )  ;  Box-Muller method (requires two x's, generates two y's) 
    if( G2bool ){
      y = G2yval;
      G2xval = x;
      G2bool = 0;
    }
    else{
      y = params[0] * sqrt(-2.0 * log(x)) * cos(2.0 * M_PI * G2xval) + params[1];
      G2yval = params[0] * sqrt(-2.0 * log(x)) * sin(2.0 * M_PI * G2xval) + params[1];
      G2bool = 1;
    }
    break;
  default: 
    fprintf(stderr,"Error in RandWithDist: type is out of bounds\n");
  }
  return( y );
}


void InitRNG(){
  srand48(seed);
}


void InitRNGPrompt(){
  /* get a seed value to start the random number generator */
  fprintf(stderr,"Enter a seed\n");
  scanf("%d",&seed);
  srand48(seed);
}


/* testing code */
int main(){
  int i,n,type;
  double x,y,params[3],*p;

  /* initialize the random number generator */
  InitRNGPrompt();

  fprintf(stderr,"Enter choice of distribution (0:exp , 1:BreitWigner , 2:Gauss)\n");
  scanf("%d",&type);
  if( type == 2 ){ G2bool = 0; G2xval = drand48(); }

  fprintf(stderr,"Enter values for params\n");
  p = params;
  scanf("%le",p);
  p++;
  scanf("%le",p);

  /* how many random numbers? */
  fprintf(stderr,"Enter # of random #'s to be made\n");
  scanf("%d",&n);

  /* create and output n random numbers with a given distribution */
  for( i=0; i<n; i++ ){
    y = RandWithDist(type,params);
    printf("%le\n",y);
    //    printf("%le %le %le\n",G2xval,G2yval,y);
  }

  return(0);
}
