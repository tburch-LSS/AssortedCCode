/* 
   PID.c 

   routines for simulating a PID controller (time domain) 

   12/2018 tb 
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "PID.h"


void GetParams(){
  fprintf(stderr,"Enter a PID time step and max time\n");
  scanf("%le",&(p_inp.dt_PID));
  scanf("%le",&(p_inp.tmax));
  fprintf(stderr,"Enter an integer ratio for PID time step / Process-simulation time step\n");
  scanf("%d",&(p_inp.Rdt));
  fprintf(stderr,"Enter Process init conditions, y0 and dydt0\n");
  scanf("%le",&(p_inp.y0));
  scanf("%le",&(p_inp.dydt0));
  fprintf(stderr,"Enter Kp, Ki, Kd, Kt\n");
  scanf("%le",&(p_inp.K[0]));
  scanf("%le",&(p_inp.K[1]));
  scanf("%le",&(p_inp.K[2]));
  scanf("%le",&(p_inp.K[3]));
  fprintf(stderr,"Enter initial Integral value\n");
  scanf("%le",&(p_inp.initI));
  fprintf(stderr,"Enter additional params 0...5\n");
  scanf("%le",&(p_inp.param[0]));
  scanf("%le",&(p_inp.param[1]));
  scanf("%le",&(p_inp.param[2]));
  scanf("%le",&(p_inp.param[3]));
  scanf("%le",&(p_inp.param[4]));
  scanf("%le",&(p_inp.param[5]));
  return;
}


/* set-point function (r) for which the Process Variable (y) should aim */
/* 
double SetPoint( double t ){
  double sp;
  sp = 1;
  return( sp );
} 
*/


void Init_dPID( dPID_info *pid ){
  (*pid).Kp = p_inp.K[0];
  (*pid).Ki = p_inp.K[1] * p_inp.dt_PID;
  (*pid).Kd = p_inp.K[2] / p_inp.dt_PID;
  (*pid).Kt = p_inp.K[3];
  (*pid).ee[0] = 0.0;
  (*pid).Prop = (*pid).Kp * (*pid).ee[0];
  (*pid).Integ[0] = p_inp.initI * (*pid).Ki;
  (*pid).Deriv[0] = 0.0;
  (*pid).vv[0] = (*pid).Prop + (*pid).Integ[0] + (*pid).Deriv[0];
  (*pid).vv[1] = 0.0;
  (*pid).yy = p_inp.y0;
  (*pid).ff = 0.0;
  (*pid).rr = SetPoint(0);
  return;
}


void dPID_update( dPID_info *pid ){
  // Error, e(t): 
  (*pid).ee[1] = (*pid).ee[0];
  (*pid).ee[0] = (*pid).rr - (*pid).yy;

  // Proportional branch: 
  (*pid).Prop = (*pid).Kp * (*pid).ee[0];

  // Derivative branch: 
  (*pid).Deriv[1] = (*pid).Deriv[0];
  // backward Euler 
  (*pid).Deriv[0] = (*pid).Kd * ((*pid).ee[0] - (*pid).ee[1]);

  // Integral branch: 
  (*pid).Integ[1] = (*pid).Integ[0];
#ifdef TRAPEZOID
  // trapezoid 
  (*pid).Integ[0] += (*pid).Ki * 0.5 * ((*pid).ee[0] + (*pid).ee[1]);
#else
  // simple rectangles,  e_old * dt 
  (*pid).Integ[0] += (*pid).Ki * (*pid).ee[1];
#endif

  // Initial PID output: v 
  (*pid).vv[1] = (*pid).vv[0];
  (*pid).vv[0] = (*pid).Prop + (*pid).Integ[0] + (*pid).Deriv[0];

  // Anti-windup: w-v 
  if( (*pid).Kt != 0 ){
    dAntiWU( pid );

    // Final PID output: control variable, u 
    (*pid).uu = (*pid).Prop + (*pid).Integ[0] + (*pid).Deriv[0] + (*pid).ff;
  }
  else (*pid).uu = (*pid).vv[0] + (*pid).ff;

  return;
}


void dAntiWU( dPID_info *pid ){
#ifdef ACTMODEL
  // Actuator model, w(v_old) 
  (*pid).ww = dActModel( pid );
#endif

  // ... or just external Tracking signal, w 

  // ... then 
  (*pid).Integ[0] += (*pid).Kt * ((*pid).ww - (*pid).vv[1]);
  return;
}


double dActModel( dPID_info *pid ){
  /* parameters p_inp.param[3..5] should be used here */
  double w;

  //  w = somefunc ( (*pid).vv[1] );
  w = (*pid).vv[1];
  return( w );
}


void Init_fPID( fPID_info *pid ){
  (*pid).Kp = (float)(p_inp.K[0]);
  (*pid).Ki = (float)(p_inp.K[1] * p_inp.dt_PID);
  (*pid).Kd = (float)(p_inp.K[2] / p_inp.dt_PID);
  (*pid).Kt = (float)(p_inp.K[3]);
  (*pid).ee[0] = 0.0;
  (*pid).Prop = (*pid).Kp * (*pid).ee[0];
  (*pid).Integ[0] = (float)(p_inp.initI) * (*pid).Ki;
  (*pid).Deriv[0] = 0.0;
  (*pid).vv[0] = (*pid).Prop + (*pid).Integ[0] + (*pid).Deriv[0];
  (*pid).vv[1] = 0.0;
  (*pid).yy = (float)(p_inp.y0);
  (*pid).ff = 0.0;
  (*pid).rr = (float)(SetPoint(0));
  return;
}


void fPID_update( fPID_info *pid ){
  // Error, e(t): 
  (*pid).ee[1] = (*pid).ee[0];
  (*pid).ee[0] = (*pid).rr - (*pid).yy;

  // Proportional branch: 
  (*pid).Prop = (*pid).Kp * (*pid).ee[0];

  // Derivative branch: 
  (*pid).Deriv[1] = (*pid).Deriv[0];
  // backward Euler 
  (*pid).Deriv[0] = (*pid).Kd * ((*pid).ee[0] - (*pid).ee[1]);

  // Integral branch: 
  (*pid).Integ[1] = (*pid).Integ[0];
#ifdef TRAPEZOID
  // trapezoid 
  (*pid).Integ[0] += (*pid).Ki * 0.5 * ((*pid).ee[0] + (*pid).ee[1]);
#else
  // simple rectangles,  e_old * dt 
  (*pid).Integ[0] += (*pid).Ki * (*pid).ee[1];
#endif

  // Initial PID output: v 
  (*pid).vv[1] = (*pid).vv[0];
  (*pid).vv[0] = (*pid).Prop + (*pid).Integ[0] + (*pid).Deriv[0];

  // Anti-windup: w-v 
  if( (*pid).Kt != 0 ){
    fAntiWU( pid );

    // Final PID output: control variable, u 
    (*pid).uu = (*pid).Prop + (*pid).Integ[0] + (*pid).Deriv[0] + (*pid).ff;
  }
  else (*pid).uu = (*pid).vv[0] + (*pid).ff;

  return;
}


void fAntiWU( fPID_info *pid ){
#ifdef ACTMODEL
  // Actuator model, w(v_old) 
  (*pid).ww = dActModel( pid );
#endif

  // ... or just external Tracking signal, w 

  // ... then 
  (*pid).Integ[0] += (*pid).Kt * ((*pid).ww - (*pid).vv[1]);
  return;
}


float fActModel( fPID_info *pid ){
  /* parameters p_inp.param[3..5] should be used here */
  float w;

  //  w = somefunc ( (*pid).vv[1] );
  w = (*pid).vv[1];
  return( w );
}

