/*
  hvy_q_pot2.c

  to compile:
  cc -o hvy_q_pot2 hvy_q_pot2.c bisect.c rk4.c -lm -DBISECT
*/

#include<math.h>
#include<stdio.h>

void derivs( double *dydx, double x, double *y, int n );
double bisect( double elo, double ehi, double tol );
void rk4( double *y, double *dydx, double x, double dx, int n );
double potential( double x );
double energy,*wf,*dwfdx,xmax,dx,alpha,mass,kk,bb,cc,dd;
int l,nsteps;


main(){
  int i;
  double x,elo,ehi,e,prec,e_conv,x_conv,integral,wf_old,psi2,dp4,p4,*wave,
    *wave1,*wave2,K,Vc,Vs,V,rdVdr;
  char e_units[10],x_units[10];
  FILE *file1,*file2;

  wf = (double *)malloc(4*sizeof(double));
  dwfdx = wf + 2;

  file1 = fopen("wf.dat","w");
  file2 = fopen("dwf2.dat","w");

  /* -hydrogen (test):
     energy in eV, e_conv=hbar^2/(2*m_e*a0)=13.6eV, x in a0, alpha=1.999
     -heavy-quark potentials:
     mass & energies in units of sqrt(string-tension)=K^0.5=440MeV
     radial distance in units of K^-0.5=0.448fm
     beta=7.40: alpha=0.39
     beta=7.75: alpha=0.32 */
  /*
  alpha = 0.32;
  mass = 2100.0;
  l = 0;
  */
  sprintf(e_units,"MeV");
  sprintf(x_units,"fm");
  e_conv = 440.0;
  x_conv = 0.448;

#ifdef BISECT
  fprintf(stderr,"Enter E_low in %s\n",e_units);
  scanf("%le",&elo);
  fprintf(stderr,"Enter E_high in %s\n",e_units);
  scanf("%le",&ehi);
  elo /= e_conv;
  ehi /= e_conv;
  fprintf(stderr,"Enter WF_precision\n");
  scanf("%le",&prec);
#else
  fprintf(stderr,"Enter E in %s\n",e_units);
  scanf("%le",&energy);
  energy /= e_conv;
#endif

  /*
  fprintf(stderr,"Enter alpha\n");
  scanf("%le",&alpha);
  */

  fprintf(stderr,"Enter kappa'/kappa\n");
  scanf("%le",&kk);
  fprintf(stderr,"Enter bb/kappa^3/2\n");
  scanf("%le",&bb);
  fprintf(stderr,"Enter cc/kappa\n");
  scanf("%le",&cc);
  fprintf(stderr,"Enter dd/kappa^1/2\n");
  scanf("%le",&dd);

  fprintf(stderr,"Enter reduced mass in %s\n",e_units);
  scanf("%le",&mass);
  mass /= e_conv;
  fprintf(stderr,"Enter l (0 = S-wave, 1 = P-wave, etc.)\n");
  scanf("%d",&l);
  fprintf(stderr,"Enter x_max in %s\n",x_units);
  scanf("%le",&xmax);
  xmax /= x_conv;
  fprintf(stderr,"Enter Nsteps\n");
  scanf("%d",&nsteps);

  wave = (double *)malloc((nsteps+1)*sizeof(double));
  wave1 = (double *)malloc((nsteps+1)*sizeof(double));
  wave2 = (double *)malloc((nsteps+1)*sizeof(double));

  dx = xmax/(double)nsteps;

#ifdef BISECT
  e = bisect(elo,ehi,prec);
#endif

  wave[0] = 0.0;
  wave1[0] = 1.0;
  wf[0] = dx;
  wf[1] = 1.0;
  integral=0.5*dx*dx*dx;
  wf_old = wf[0];
  for( i=1; i<nsteps; i++ ){
    x = (double)i*dx;
    derivs(dwfdx,x,wf,2);
    wave[i] = wf[0];
    wave1[i] = wf[1];
    wave2[i] = dwfdx[1];
    rk4(wf,dwfdx,x,dx,2);
    integral += 0.5*dx*(wf_old*wf_old+wf[0]*wf[0]);
    wf_old = wf[0];
  }
  wave[i] = wf[0];
  wave1[i] = wf[1];
  wave2[i] = dwfdx[1];
  wave2[0] = 2.0*wave2[1] - wave2[2];
  psi2 = 1.0/(4.0*acos(-1.0)*integral*x_conv*x_conv*x_conv);

  fprintf(file1,"%s alpha=%le ; m=%le %s ; l=%d ; kappa^1/2 = %le %s\n","%",
	  alpha,mass*e_conv,e_units,l,e_conv,e_units);

  fprintf(stderr,"wavefunction(%le %s) = %le\n",xmax*x_conv,x_units,wf[0]);
  fprintf(stderr,"energy = %le %s\n",energy*e_conv,e_units);
  fprintf(file1,"%s energy = %le %s\n","%",energy*e_conv,e_units);
  fprintf(stderr,"psi2(0) = %le %s^-3 (only valid for ground state)\n",
	  psi2,x_units);
  fprintf(file1,"%s psi2(0) = %le %s^-3 (only valid for ground state)\n",
	  "%",psi2,x_units);
  fprintf(stderr,"psi(0) = %le %s^-1.5 (only valid for ground state)\n",
	  sqrt(psi2),x_units);
  fprintf(file1,"%s psi(0) = %le %s^-1.5 (only valid for ground state)\n",
	  "%",sqrt(psi2),x_units);

  i = nsteps/1000;
  x = (double)i*dx;
  p4 = 0.5*x*wave[i]*(-wave2[i-2] + 16.0*wave2[i-1] - 30.0*wave2[i]
		      + 16.0*wave2[i+1] - wave2[i+2])/(12.0*dx*dx);
  p4 -= ((double)(l*(l+1)))*wave[i]*wave2[i]/x;
  p4 += 2.0*((double)(l*(l+1)))*wave[i]*wave1[i]/(x*x);
  p4 += ((double)(l*(l+1)))*(((double)(l*(l+1))-6.0))*wave[i]*wave[i]/(x*x*x);
  K = 0.5*x*wave[i]*wave2[i];
  K -= 0.5*((double)(l*(l+1)))*wave[i]*wave[i]/x;
  V = 0.5*x*potential(x)*wave[i]*wave[i];
  /*
  Vc = 0.5*(-alpha)*wave[i]*wave[i];
  Vs = 0.5*(x*x)*wave[i]*wave[i];
  */
  /*
  Vs = 0.5*x*sqrt(x*x+bb*x+cc)*wave[i]*wave[i];
  rdVdr = 0.5*x*(x*x+0.5*bb*x)*wave[i]*wave[i]/sqrt(x*x+bb*x+cc);
  */
  for( ; i<nsteps; i++ ){
    x = (double)i*dx;
    dp4 = wave[i]*(-wave2[i-2] + 16.0*wave2[i-1] - 30.0*wave2[i]
		   + 16.0*wave2[i+1] - wave2[i+2])/(12.0*dx);
    dp4 -= 2.0*((double)(l*(l+1)))*dx*wave[i]*wave2[i]/(x*x);
    dp4 += 4.0*((double)(l*(l+1)))*dx*wave[i]*wave1[i]/(x*x*x);
    dp4 += ((double)(l*(l+1)))*(((double)(l*(l+1))-6.0))
      *dx*wave[i]*wave[i]/(x*x*x*x);
    p4 += dp4;
    K += wave[i]*wave2[i]*dx;
    K -= ((double)(l*(l+1)))*dx*wave[i]*wave[i]/(x*x);
    V += dx*potential(x)*wave[i]*wave[i];
    /*
    Vc += dx*(-alpha/x)*wave[i]*wave[i];
    Vs += dx*x*wave[i]*wave[i];
    */
    /*
    Vs += dx*sqrt(x*x+bb*x+cc)*wave[i]*wave[i];
    rdVdr += dx*(x*x+0.5*bb*x)*wave[i]*wave[i]/sqrt(x*x+bb*x+cc);
    */
    /*
    printf("%le %le\n",x*x_conv,K);
    */
  }
  p4 *= -e_conv/(8.0*mass*mass*mass*integral);
  K *= -e_conv/(2.0*mass*integral);

  rdVdr = Vs - Vc;

  V *= e_conv/integral;
  Vc *= e_conv/integral;
  Vs *= e_conv/integral;
  rdVdr *= e_conv/integral;
  fprintf(stderr,"<psi|-p^4/8m^3|psi> = %le %s\n",p4,e_units);
  fprintf(file1,"%s <psi|-p^4/8m^3|psi> = %le %s\n","%",p4,e_units);
  fprintf(stderr,"<psi|p^2/2m|psi> = %le %s\n",K,e_units);
  fprintf(file1,"%s <psi|p^2/2m|psi> = %le %s\n","%",K,e_units);
  fprintf(stderr,"<psi|V(r)|psi> = %le %s\n",V,e_units);
  fprintf(file1,"%s <psi|V(r)|psi> = %le %s\n","%",V,e_units);
  /*
  fprintf(stderr,"<psi|-alpha/r|psi> = %le %s\n",Vc,e_units);
  fprintf(file1,"%s <psi|-alpha/r|psi> = %le %s\n","%",Vc,e_units);
  fprintf(stderr,"<psi|kappa*r|psi> = %le %s\n",Vs,e_units);
  fprintf(file1,"%s <psi|kappa*r|psi> = %le %s\n","%",Vs,e_units);
  */
  /*
  fprintf(stderr,"<psi|sqrt(kappa^2*r^2 + bb*r + cc)|psi> = %le %s\n",Vs,
	  e_units);
  fprintf(file1,"%s <psi|sqrt(kappa^2*r^2 + bb*r + cc)|psi> = %le %s\n","%",
	  Vs,e_units);
  */
  /*
  fprintf(stderr,"<psi|r(dV/dr)/2|psi> = %le %s\n",0.5*rdVdr,e_units);
  fprintf(file1,"%s <psi|r(dV/dr)/2|psi> = %le %s\n","%",0.5*rdVdr,e_units);
  */

  for( i=0; i<nsteps; i++ ){
    x = (double)i*dx;
    fprintf(file1,"%le %le\n",x*x_conv,wave[i]/sqrt(integral));
    fprintf(file2,"%le %le\n",x*x_conv,wave2[i]/sqrt(integral));
  }

  free(wf);
  free(wave); free(wave1); free(wave2);
  fclose(file1); fclose(file2);
}


/* f(E) for bisection method */
double f( double e ){
  double x;
  int i;

  energy = e;
  wf[0] = dx;
  wf[1] = 1.0;

  for( i=1; i<nsteps; i++ ){
    x = (double)i*dx;
    rk4(wf,dwfdx,x,dx,2);
  }
  fprintf(stderr,"wavefunction(end) = %le ; energy = %.13f\n",wf[0],energy);

  return wf[0];
}


/*** DE's used by rk4 ***/
void derivs( double *dydx, double x, double *y, int n ){
  dydx[0] = y[1];
  dydx[1] = ((double)(l*(l+1))/(x*x) + 2.0*mass*(potential(x) - energy))*y[0];
}


double potential( double x ){
  /* singlet */
  /*
  if( x > 0.223361 )
    return(x - 0.3058/x);
  else
    return(0.223361 - 0.3058/0.223361 - 1.71399 + 1.71399*x/0.223361);
  */
  /*
  if( x > 0.315595 )
    return(x - 0.2967/x);
  else
    return(0.315595 - 0.2967/0.315595 - 1.3788 + 1.3788*x/0.315595);
  */
  /* hybrid (octet) */
  return(sqrt( kk*kk*x*x + bb*x + cc ) + dd);
}
