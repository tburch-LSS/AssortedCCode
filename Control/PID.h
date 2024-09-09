
#ifdef MAIN
#define EXTERN 
#else
#define EXTERN extern
#endif


// main structure for current status of PID controller 
typedef struct{

  // PID parameters (1st 3 for P,I,D branches; last for anti-windup: w-v in I branch) 
  // Ki and Kd will absorb the time step (Ki --> dt*Ki and Kd --> Kd/dt) 
  double Kp,Ki,Kd,Kt;

  // Propotional (P), Integral (I), and Derivative (D) terms  [*to,of e(t))] 
  double Prop,Integ[2],Deriv[2];

  // current [0] and previous [1] errors, e(t) 
  double ee[2];

  // inputs: set-point (r), process variable (y), feedforward (f), and tracking signal (w) 
  double rr,yy,ff,ww;

  // outputs: control variable (u) and initial PID output (v) 
  double uu,vv[2];

} dPID_info;


// same, in single precision 
typedef struct{
  float Kp,Ki,Kd,Kt;
  float Prop,Integ[2],Deriv[2];
  float ee[2];
  float rr,yy,ff,ww;
  float uu,vv[2];
} fPID_info;


typedef struct{
  // time step of PID and max time of sim 
  double dt_PID,tmax;
  // ratio of dt_PID / dt 
  int Rdt;
  // initial conditions for Process 
  double y0,dydt0;
  // initial K parameters of PID 
  double K[4];
  // initial value of Integral term 
  double initI;
  // additional parameters; e.g., for Process DEs (0-2) or Actuator Model (3-5) 
  double param[6];
} initParams;

EXTERN initParams p_inp;

// EXTERN dPID_info dpid;
// EXTERN fPID_info fpid;

void GetParams();
double SetPoint( double t );

void Init_dPID( dPID_info *pid );
void dPID_update( dPID_info *pid );
void dAntiWU( dPID_info *pid );
double dActModel( dPID_info *pid );

void Init_fPID( fPID_info *pid );
void fPID_update( fPID_info *pid );
void fAntiWU( fPID_info *pid );
float fActModel( fPID_info *pid );
