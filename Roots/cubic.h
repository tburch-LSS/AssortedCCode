#include "zcomplex.h"

#define THIRD 0.3333333333333
#define NINTH 0.1111111111111
#define TWENSEVTH 0.03703703703704
#define HALFRADTHREE 0.866025404

double discriminant( double *C );
void cubic_solve ( dcomplex *x, double *C );
void cubic_solve_real ( double *x, double *C );
