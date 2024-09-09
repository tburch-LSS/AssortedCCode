/*
  test_root.c

  tests root finding algorithms

  to compile:   cc -o test_root test_root.c root_find.c -lm
  to execute:   test_root

  3/2000 tb
*/



#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include "root_find.h"


main(){
  double x,a,b,p;
  int n;

  a = -0.1;
  b = 0.9;
  p = 1.e-8;
  x = bisect(a,b,p,&n);
  printf("bisect root = %e in %d iters.\n",x,n);

  x = -0.5;
  x = newt_raph(x,p,&n);
  printf("NR root = %e in %d iters.\n",x,n);
}


double function( double x ){
  return sin(x);
}


double dfdx( double x ){
  return cos(x);
}
