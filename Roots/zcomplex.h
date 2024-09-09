
typedef struct { double real; double imag; } dcomplex;
typedef struct { float real; float imag; } fcomplex;

#define zmag(z) ( sqrt( (z).real*(z).real + (z).imag*(z).imag ) )
#define zarg(z) ( atan2( (z).imag , (z).real ) )
#define zmag2(z) ( (z).real*(z).real + (z).imag*(z).imag )

#define ZEQ(z,c) { (z).real = (c).real; (z).imag = (c).imag; }
#define ZEQA(z,c) { (z).real = (c).real; (z).imag = -(c).imag; }
#define ZCMPLX(z,x,y) { (z).real = (x); (z).imag = (y); }
#define ZADD(c,a,b) { (c).real = (a).real + (b).real;	\
                      (c).imag = (a).imag + (b).imag; }
#define ZSUB(c,a,b) { (c).real = (a).real - (b).real;	\
                      (c).imag = (a).imag - (b).imag; }
#define ZPEQ(a,b) { (a).real += (b).real; (a).imag += (b).imag; }
#define ZMEQ(a,b) { (a).real -= (b).real; (a).imag -= (b).imag; }
#define ZMUL(c,a,b) { (c).real = (a).real*(b).real - (a).imag*(b).imag; \
                      (c).imag = (a).real*(b).imag + (a).imag*(b).real; }
#define ZDIV(c,a,b) { double t = (b).real*(b).real + (b).imag*(b).imag; \
                      (c).real = ((a).real*(b).real + (a).imag*(b).imag)/t; \
                      (c).imag = ((a).imag*(b).real - (a).real*(b).imag)/t; }
#define ZDIVREAL(z,r,c) { double t = (c).real*(c).real + (c).imag*(c).imag; \
                      (z).real = (r)*(c).real / t; \
                      (z).imag = -(r)*(c).imag / t; }
#define ZMULREAL(z,r,c) { (z).real = (r)*(c).real; (z).imag = (r)*(c).imag; }
#define ZMULIMAG(z,i,c) { (z).real = -(i)*(c).imag; (z).imag = (i)*(c).real; }
#define ZADDREAL(z,r,c) { (z).real = (r)+(c).real; (z).imag = (c).imag; }
#define ZADDIMAG(z,i,c) { (z).real = (c).real; (z).imag = (i)+(c).imag; }
#define ZMULEQREAL(z,r) { (z).real *= (r); (z).imag *= (r); }
#define ZPEQREAL(z,r) { (z).real += (r); }
#define ZPEQIMAG(z,i) { (z).imag += (i); }

