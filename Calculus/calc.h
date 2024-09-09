
double function( double x );
double integrand( double x );
double integrand2( double x, double y );
double integrand3( double x, double y, double z );
double trapezoid( double a, double b, int n );
double trap_2d( double a, double b, double c, double d, int nx, int ny );
double trap_3d( double *a, double *b, int *nxyz );
double romberg( double a, double b, int n, int k );
double romberg_fast( double a, double b, int n, int k );
double gauss_legendre( double a, double b, int n );
double gauss_laguerre( double a, int n );
double cent_diff( double x, double h );
double richardson( double x, double h, int k );
