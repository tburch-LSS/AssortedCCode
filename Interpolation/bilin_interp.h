
#ifdef MAIN
#define EXTERN 
#else
#define EXTERN extern
#endif

double bilin_interpol( double fracx1, double fracy1, 
		       double z11, double z12, double z21, double z22 );
double inv_bilin_interpol( double x1, double x2, 
			   double fracy1, double z, 
			   double z11, double z12, double z21, double z22 );
