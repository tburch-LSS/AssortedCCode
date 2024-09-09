/*
  orbit.c

  computes the orbit of a mass moving in the presence of a central potential:
  ( e.g., U/m = -GM/r )

  to compile:   cc -o orbit orbit.c ode.c -lm
  to execute:   orbit > orbit.dat

  output is 7 columns of numbers:
  t  x  y  K/m  U/m  E/m  L/m

  to make a 2-column data file from this: e.g., to make an (x,y) plot use:
  awk '{print $2,$3}' orbit.dat > xy.dat

  9/2000 tb
*/

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include "ode.h"

double GM,L_over_m;


main(){
  double *r,*theta,*vr,*dtheta_dt,*dvr_dt,dt,t,tmax,x,y,K,U,E,Lz,pi;
  int i,j,nt;

  pi = acos(-1.0);
  GM = 4.0*pi*pi;   /* in AU & yr */

  /* allocate adjacent memory for r, theta, v_r=dr/dt, dtheta/dt, dv_r/dt */
  r = (double *)malloc(5*sizeof(double));
  theta = r + 1;
  vr = theta + 1;
  dtheta_dt = vr + 1;
  dvr_dt = dtheta_dt + 1;

  /* get total time, time step, and initial conditions */
  fprintf(stderr,"Enter t_max (in yr)\n");
  scanf("%le",&tmax);
  fprintf(stderr,"Enter dt (in yr)\n");
  scanf("%le",&dt);
  fprintf(stderr,"Enter r0 (in AU)\n",i);
  scanf("%le",r);
  fprintf(stderr,"Enter r0*(dtheta/dt)0 (in AU/yr)\n",i);
  scanf("%le",&L_over_m);
  *dtheta_dt = L_over_m/(*r);
  L_over_m *= (*r);
  *theta = 0.0;
  *vr = 0.0;
  t = 0.0;

  /* calculate initial acceleration */
  derivs(vr,t,r,3);

  /* calculate # of time intervals */
  nt = (int)(tmax/dt + 0.5);

  /* loop over time steps */
  for( j=0; j<=nt; j++ ){

    /* calculate the time */
    t = j*dt;

    /* calculate x, y, U/m, K/m, E/m, & L/m and print output */
    x = (*r)*cos(*theta);
    y = (*r)*sin(*theta);
    U = -GM/(*r);
    K = 0.5*((*vr)*(*vr) + (*r)*(*r)*(*dtheta_dt)*(*dtheta_dt));
    E = K + U;
    Lz = (*r)*(*r)*(*dtheta_dt);
    printf("%le %le %le %le %le %le %le\n",t,x,y,K,U,E,Lz);

    /* update the position and velocity */
    rk2(r,vr,t,dt,3);

  }

  free(r);
}


/* differential equations */
void derivs( double *dydt, double t, double *y, int n ){

  /*   dr/dt = v_r   */
  dydt[0] = y[2];

  /*   dtheta/dt = (L/m)/r^2   */
  dydt[1] = L_over_m/(y[0]*y[0]);

  /*   dv_r/dt = (1/r^3)(L/m)^2 - GM/r^2   */
  dydt[2] = L_over_m*L_over_m/(y[0]*y[0]*y[0]) - GM/(y[0]*y[0]);
}
