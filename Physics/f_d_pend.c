/*
  f_d_pend.c

  simulation of a forced, damped pendulum

  to compile:   cc -o f_d_pend f_d_pend.c ode.c -lm
  to execute:   f_d_pend > f_d_pend.dat

  output is 3 columns:
  t  theta  omega

  to make a two-column file to plot, e.g., t vs. theta:
  awk '{print $1,$2}' f_d_pend.dat > t_vs_theta.dat

  3/2000 tb
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ode.h"

/* driving frequency & amplitude, acceleration due to gravity, length,
   and damping constant */
double w,a,g,l,b;


main(){
  double t,*theta,*omega,dt,t_final,pi;
  int i,nt;

  /* allocate adjacent memory locations for theta, omega=dtheta/dt, and
     domega/dt */
  theta = (double *)malloc(3*sizeof(double));
  omega = theta + 1;

  pi = acos(-1.0);

  /* get parameters and initial conditions */
  g = 980.0;   /* cm/s^2 */
  fprintf(stderr,"Enter a final time (s)\n");
  scanf("%le",&t_final);
  fprintf(stderr,"Enter a time step (s)\n");
  scanf("%le",&dt);
  fprintf(stderr,"Enter an initial angle (rad)\n");
  scanf("%le",theta);
  fprintf(stderr,"Enter an initial angular velocity (rad/s)\n");
  scanf("%le",omega);
  fprintf(stderr,"Enter a length (cm)\n");
  scanf("%le",&l);
  fprintf(stderr,"Enter a driving frequency (1/s)\n");
  scanf("%le",&w);
  fprintf(stderr,"Enter a driving amplitude (cm)\n");
  scanf("%le",&a);
  fprintf(stderr,"Enter a damping constant (1/(cm s))\n");
  scanf("%le",&b);

  /* calculate number of time steps needed */
  nt = (int)(t_final/dt + 0.5);

  /* loop over time */
  for( i=0; i<=nt; i++ ){

    /* keep theta between +pi and -pi */
    while( *theta < -pi ) *theta += 2.0*pi;
    while( *theta > pi ) *theta -= 2.0*pi;

    /* calculate time */
    t = i*dt;

    /* write out t, theta and omega */
    printf("%e %e %e\n",t,theta[0],omega[0]);

    /* update theta & omega */
    rk4(theta,omega,t,dt,2);

  }

}


void derivs( double *dydt, double t, double *y, int n ){

  /*  dtheta/dt = omega  */
  dydt[0] = y[1];

  /*  domega/dt = -(g/l)sin(theta) - b*omega + (aw^2/l)sin(wt)sin(theta)  */
  dydt[1] = -(g/l)*sin(y[0]) - b*y[1] + (a*w*w/l)*sin(w*t)*sin(y[0]);
}
