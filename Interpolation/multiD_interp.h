
// only up to 5D now... 
static int maxD=5;
static int Ncorn[6]={0,2,4,8,16,32};
static int Astart[6]={0,1,3,11,35,99};

// each block assoc. w/ a h-cuboid corner; within each block are the above-the-internal-point (x) truth values in each dimension 
// first index fastest running 
static bool Above[259]={0, 
			0,         1, 
			0,0,       1,0,       0,1,       1,1, 
			0,0,0,     1,0,0,     0,1,0,     1,1,0,     0,0,1,     1,0,1,     0,1,1,     1,1,1, 
			0,0,0,0,   1,0,0,0,   0,1,0,0,   1,1,0,0,   0,0,1,0,   1,0,1,0,   0,1,1,0,   1,1,1,0, 
			0,0,0,1,   1,0,0,1,   0,1,0,1,   1,1,0,1,   0,0,1,1,   1,0,1,1,   0,1,1,1,   1,1,1,1, 
			0,0,0,0,0, 1,0,0,0,0, 0,1,0,0,0, 1,1,0,0,0, 0,0,1,0,0, 1,0,1,0,0, 0,1,1,0,0, 1,1,1,0,0, 
			0,0,0,1,0, 1,0,0,1,0, 0,1,0,1,0, 1,1,0,1,0, 0,0,1,1,0, 1,0,1,1,0, 0,1,1,1,0, 1,1,1,1,0, 
			0,0,0,0,1, 1,0,0,0,1, 0,1,0,0,1, 1,1,0,0,1, 0,0,1,0,1, 1,0,1,0,1, 0,1,1,0,1, 1,1,1,0,1, 
			0,0,0,1,1, 1,0,0,1,1, 0,1,0,1,1, 1,1,0,1,1, 0,0,1,1,1, 1,0,1,1,1, 0,1,1,1,1, 1,1,1,1,1  };

// each block assoc. w/ a h-cuboid corner; within each block are the above-the-internal-point (x) truth values in each dimension 
// first index slowest running? (probably not) 
/* 
static bool Above2[259]={0, 
			 0,         1, 
			 0,0,       0,1,       1,0,       1,1, 
			 0,0,0,     0,0,1,     0,1,0,     0,1,1,     1,0,0,     1,0,1,     1,1,0,     1,1,1, 
			 0,0,0,0,   0,0,0,1,   0,0,1,0,   0,0,1,1,   0,1,0,0,   0,1,0,1,   0,1,1,0,   0,1,1,1, 
			 1,0,0,0,   1,0,0,1,   1,0,1,0,   1,0,1,1,   1,1,0,0,   1,1,0,1,   1,1,1,0,   1,1,1,1, 
			 0,0,0,0,0, 0,0,0,0,1, 0,0,0,1,0, 0,0,0,1,1, 0,0,1,0,0, 0,0,1,0,1, 0,0,1,1,0, 0,0,1,1,1, 
			 0,1,0,0,0, 0,1,0,0,1, 0,1,0,1,0, 0,1,0,1,1, 0,1,1,0,0, 0,1,1,0,1, 0,1,1,1,0, 0,1,1,1,1, 
			 1,0,0,0,0, 1,0,0,0,1, 1,0,0,1,0, 1,0,0,1,1, 1,0,1,0,0, 1,0,1,0,1, 1,0,1,1,0, 1,0,1,1,1, 
			 1,1,0,0,0, 1,1,0,0,1, 1,1,0,1,0, 1,1,0,1,1, 1,1,1,0,0, 1,1,1,0,1, 1,1,1,1,0, 1,1,1,1,1  }; 
*/

void set_corners_fracs( double *zcorn , double *fracx, 
			double *ztab , int *ntab , double *x , double *xstart , double *dx, int dim );
double Dlin_interpol( double *fracx , double *zcorn , int dim );
double InvFracDistHcub_interpol( double *fracx , double *zcorn , double p , int dim );

void set_corners_deltas( double *zcorn , double *deltax, 
			 double *ztab , int *ntab , double *x , double *xstart , double *dx, int dim );
double InvDistHcub_interpol( double *deltax , double *zcorn , double *dx , double p , int dim );

double InvDistAll_interpol( double *ztab , int *ntab , double *x , double *xstart , double *dx , double p , int dim );

double RBFunction( double r , double *params );
double RadBasFuncScat_interpol( double *ztab , double *xtab , int ntot , double *params , double *x , int dim );
