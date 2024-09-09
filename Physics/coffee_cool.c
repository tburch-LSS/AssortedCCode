/*
  coffee_cool.c

  to compile:   cc -o coffee_cool coffee_cool.c ode.c -lm
  to execute:   coffee_cool > coffee_cool.dat

  9/2000 tb
*/

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include "ode.h"

double r,Tr;


int main(){
  double t,dt,tf,T0,T,dTdt,T_exact,error;
  int i,nt;

  /* get parameters and initial values */
  fprintf(stderr,"Enter a cooling parameter (min^-1)\n");
  scanf("%le",&r);
  fprintf(stderr,"Enter room temperature (C)\n");
  scanf("%le",&Tr);
  fprintf(stderr,"Enter an initial temperature (C)\n");
  scanf("%le",&T0);
  T = T0;
  fprintf(stderr,"Enter a final time (min)\n");
  scanf("%le",&tf);
  fprintf(stderr,"Enter a time step (min)\n");
  scanf("%le",&dt);

  /* calculate number of time steps (if tf/dt is not an integer, O(dt)
     error will appear in the final result) */
  nt = (int)(tf/dt + 0.5);

  /* loop over time steps */
  for( i=0; i<=nt; i++ ){

    /* calculate time */
    t = dt*i;

    /* print results to stdout */
    printf("%le %le\n",t,T);

    /* update temperature */
    if( i != nt ) rk2(&T,&dTdt,t,dt,1);

  }

  T_exact = Tr + (T0 - Tr)*exp(-r*tf);
  error = fabs(T_exact - T);
  fprintf(stderr,"T_exact = %24.22f\n",T_exact);
  fprintf(stderr,"T_numer = %24.22f\n",T);
  fprintf(stderr,"error   = %24.22f\n",error);

  return 0;
}


/* coffee-cooling differential equation */
void derivs( double *dTdt, double t, double *T, int n ){
  dTdt[0] = -r*(T[0]-Tr);
}
