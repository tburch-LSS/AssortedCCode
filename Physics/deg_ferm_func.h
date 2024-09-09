
#ifdef MAIN
#define EXTERN 
#else
#define EXTERN extern
#endif

EXTERN int TP_check;
EXTERN double xP,yP,eP,fP,Beta,Eta,x_Max,x_Deg,x_End;

double Integral_1( double beta, double eta );
double dintegrand1_dx_zero( double x );
void x_range_from_integrand1( double beta, double eta );
double Integral_2( double beta, double eta );
double Integral_3( double beta, double eta );
double phi( double x );
double chi( double x );
void elec_stat_corr( double I1, int Z, double *I2_corr, double *I3_corr );
double dI_dbeta_deta( double beta, double eta, 
		      double xp, double yp, double ep, double fp );
double eta_from_beta_I1( double beta, double I1 );
double Thermo_Potent_over_MC2( double beta, double eta );
