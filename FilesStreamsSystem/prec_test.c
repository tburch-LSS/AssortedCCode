/*
  prec_test.c

  tests machine precision of single- and double-precision floating point
  numbers by adding (base)^(power) to 1 and writing out the result to see
  if it makes a difference ('power' should be < 0)

  to compile:   cc -o prec_test prec_test.c -lm
  to execute:   prec_test

  3/2000 tb
*/

#include<math.h>
#include<stdio.h>

main(){
  double d_one,base,power;
  float one;

  d_one = 1.0;
  fprintf(stderr,"Enter a base\n");
  scanf("%le",&base);
  fprintf(stderr,"Enter a power ( < 0 )\n");
  scanf("%le",&power);

  d_one += pow(base,power);
  one = (float)d_one;

  printf("single-precision 1.0+base^power = %24.22f\n",one);
  printf("double-precision 1.0+base^power = %24.22f\n",d_one);
}
