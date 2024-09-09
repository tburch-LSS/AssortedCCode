
int RKWrkspcSize( int n , int ord );
void derivs( double *dydx, double x, double *y, int n );
void LinDecomp( double *uu, double **AA, double x, int n );
void DerivsJacob( double *ff, double **dfdy, double x, double *y, int n );
void rk2( double *y, double *wrkspc, double x, double dx, int n );
void rk4( double *y, double *wrkspc, double x, double dx, int n );
void TrapLinSolve( double *y, double *wrkspc, double x, double dx, int n );
void TrapAlgSolve( double *y, double *wrkspc, double x, double dx, int n );
void DiffAlgSolve( double *y, double *wrkspc, double x, double dx, int n , int ord );
