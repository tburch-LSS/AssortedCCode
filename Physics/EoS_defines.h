
#define PI 3.14159265358979
#define INV_PI2 0.101321183642

#define MeV_per_ERG 6.241509647e5

/* (m_e c^2) in ergs */
#define MeC2 8.187104385e-7
/* (m_p c^2) in ergs */
#define MpC2 1.503277360e-3
/* (m_n c^2) in ergs */
#define MnC2 1.505349508e-3
/* (m_u c^2) in ergs (^12_C scale) */
#define MuC2 1.492417830e-3
/* (m_n - m_p - m_e)c^2 in ergs */
#define DELTA 1.253437562e-6

/* 1/lambda_e^3 in cm^-3 */
#define INV_Le3 1.736603307e31
/* 1/lambda_p^3 in cm^-3 */
#define INV_Lp3 1.075045860e41
/* 1/lambda_n^3 in cm^-3 */
#define INV_Ln3 1.079497585e41

/* c in cm/s */
#define c_LIGHT 2.99792458e10
/* e^2 / (4 pi eps_0) in erg*cm */
#define K_e 2.307077130e-19
/* k_B in erg/K */
#define k_B 1.3806504e-16
/* sigma_B in erg/s/cm^2/K^4 */
#define sigma_B 5.670400e-5
/* a_0 (Bohr radius) in cm */
#define a_0 5.291772086e-9

/* least-squares-fit coefficients for the semi-empirical mass formula 
   (in ergs) */
#define a_VOLU 2.531e-5
#define a_SURF 2.932e-5
#define a_COUL 1.144e-6
#define a_ASYM 3.717e-5
#define a_PAIR 1.923e-5

/* max Z,N,A values and the length of the nuclei mass table */
#define ZMAX 118
#define NMAX 176
#define AMAX 293
#define TBL_LENGTH 3178
