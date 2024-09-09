/* 
   PIDtest.c 

   test simulation for a digital PID controller (time domain) 

   needs to be compiled with PID.c and ode.c 

   12/2018 tb 
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MAIN 
#include "PID.h"
#include "ode.h"

// Control Variable, u(t), made available here for Process DEs 
double uCV;


/* set-point function (r) for which the Process Variable (y) should aim */
double SetPoint( double t ){
  double sp;

  // unit step at t = 0 
  //  if( t < 0 ) sp = 0.0;
  //  else sp = 1.0;

  // unit step at t = dt_PID 
  if( t < p_inp.dt_PID ) sp = 0.0;
  else sp = 1.0;

  return( sp );
}


int main(){
  int i,imax,n=0;
  double t,dt,yy[6],*dydt;
  dPID_info pid;

  // "yy" array contains the Process Variable, y, possibly related (through 1st-order DEs) functions, 
  //  and their time derivs (in "dydt" part) 
  dydt = yy + 1;
  //  dydt = yy + 2;

  GetParams();
  Init_dPID( &pid );

  // initialize the Process 
  yy[0] = p_inp.y0;
  dydt[0] = p_inp.dydt0;

  // time step for simulation of Process (dt << dt_PID) 
  dt = p_inp.dt_PID / (double)(p_inp.Rdt);

  imax = (int)( p_inp.tmax / dt );
  for( i=0 ; i<imax ; i++ ){
    t = i*dt;

    // when another dt_PID has passed, update the PID controller 
    if( (n * p_inp.dt_PID - t) < 0.5*dt ){
      n++;
      // update set-point if necessary 
      pid.rr = SetPoint( t );
      // measurement of process 
      pid.yy = yy[0];
      // PID updates the control variable 
      dPID_update( &pid );
      uCV = pid.uu;
    }

    // output  t, e, u, x, y 
    printf("%le %le %le %le %le\n",t,pid.ee[0],pid.uu,yy[0],pid.yy);

    // update the Process 
    rk4( yy , dydt , t , dt , 2 );
  }

  return(0);
}


/* Process DEs */
void derivs( double *dydt , double t , double *y , int nn ){
  /* Control Variable, uCV, and parameters, p_inp.param[0..2], should be used here: */
  //  2nd-order process with u(t) as forcing term: ( dydt = yy + 1 should stand above ) 
  dydt[0] = y[1];
  dydt[1] = p_inp.param[0] * uCV - p_inp.param[1] * y[1] - p_inp.param[2] * y[0];
  //  resonant (sinusoidal) process: ( dydt = yy + 1 should stand above ) 
  //  dydt[0] = y[1];
  //  dydt[1] = p_inp.param[0] * uCV - p_inp.param[1] * y[0];
}

