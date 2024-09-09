/*
  ode2.c

  routines for integrating ordinary differential equations;
  requires the function 'derivs' or 'LinDecomp' (possibly external) 
  which calculates the dydx's using the DE's

  3/2000 ode.c , tb
  1/2022 altered for more general implicit RK methods and for external workspace, tb
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "linalg.h"
#include "ode2.h"

/*
void derivs( double *dydx, double x, double *y, int n ){
  double ;

  dydx[0] = ;
  dydx[1] = ;
}

void LinDecomp( double *uu, double **AA, double x, int n ){
  // dy_i/dx = uu_i(x) + AA_ij(x) y_j 
  double ;

  uu[0] = ; uu[1] = ;
  AA[0][0] = ; AA[0][1] = ;
}

void DerivsJacob( double *ff, double **dfdy, double x, double *y, int n ){
  // dy_i/dx = f_i(x,y) ~ f_i(x,y0) + df_i/dy_j (y_j-y0_j) + ... 
  double ;

  ff[0] = ; ff[1] = ;
  dfdy[0][0] = ; dfdy[0][1] = ;
}
*/


/* determine necessary workspace size */
int RKWrkspcSize( int n , int ord ){
  return((ord+1)*(n+1)*n);
}


void rk2( double *y, double *wrkspc, double x, double dx, int n ){
  /* 2nd order Runge-Kutta algorithm for one step (dx) of integration */
  int i;
  double k,*dydx,*yold;

  /* workspace pointers */
  dydx = wrkspc; yold = dydx + n;

  derivs(dydx,x,y,n);
  for( i=0; i<n; i++ ){
    k = dx*dydx[i];
    yold[i] = y[i];
    y[i] = yold[i] + 0.5*k;
  }
  x += 0.5*dx;
  derivs(dydx,x,y,n);
  for( i=0; i<n; i++ ){
    k = dx*dydx[i];
    y[i] = yold[i] + k;
  }
  return;
}


void rk4( double *y, double *wrkspc, double x, double dx, int n ){
  /* 4th order Runge-Kutta algorithm for one step (dx) of integration */
  int i;
  double k,*dydx,*yold,*dy;

  /* workspace pointers */
  dydx = wrkspc; yold = dydx + n; dy = yold + n;

  derivs(dydx,x,y,n);
  for( i=0; i<n; i++ ){
    k = dx*dydx[i];
    yold[i] = y[i];
    y[i] = yold[i] + 0.5*k;
    dy[i] = 0.16666666666667*k;
  }
  x += 0.5*dx;
  derivs(dydx,x,y,n);
  for( i=0; i<n; i++ ){
    k = dx*dydx[i];
    y[i] = yold[i] + 0.5*k;
    dy[i] += 0.33333333333333*k;
  }
  derivs(dydx,x,y,n);
  for( i=0; i<n; i++ ){
    k = dx*dydx[i];
    y[i] = yold[i] + k;
    dy[i] += 0.33333333333333*k;
  }
  x += 0.5*dx;
  derivs(dydx,x,y,n);
  for( i=0; i<n; i++ ){
    k = dx*dydx[i];
    dy[i] += 0.16666666666667*k;
    y[i] = yold[i] + dy[i];
  }
  return;
}


void TrapLinSolve( double *y, double *wrkspc, double x, double dx, int n ){
  /* Implicit Trapezoidal algorithm for one step (dx) of integration where the terms and cross-terms 
     involving the dependent variables are all Linear: 
      dy_i/dx = u_i(x) + A_ij(x) y_j = f_i(x,y) 
       ==> 
      [I - h/2 A(x+h)]_ij y_j(x+h) = y_i(x) + h/2 [u_i(x) + u_i(x+h) + A_ik(x) y_k(x)]  ...or... 
      [I - h/2 dfdy(x+h,0)]_ij y_j(x+h) = y_i(x) + h/2 [f_i(x,y(x)) + f_i(x+h,0)]  */
  int i,j,ij;
  double *dydx,*uu,*vv,*BB,**AA;
  
  /* workspace pointers */
  dydx = wrkspc; uu = dydx + n; vv = uu + n; BB = vv + n;
  AA = (double **)malloc(n*sizeof(double *));
  AA[0] = BB + (n*n); for(i=1;i<n;i++) AA[i] = AA[0] + (i*n);

  /*  u_i(x) + A_ik(x) y_k(x)  */
  derivs( dydx, x, y, n );
  x += dx;
  /*  u_i(x+h) ; A_ik(x+h)  */
  LinDecomp( uu, AA, x, n );
  for( i=0,ij=0 ; i<n ; i++ ){
    for( j=0 ; j<n ; j++,ij++ ){
      /*  B_ij(x+h) = [I - h/2 A(x+h)]_ij  */
      BB[ij] = -0.5*dx*AA[i][j];
      if(i==j) BB[ij] += 1.0;
    }
    /*  v_i(x,x+h) = y_i(x) + h/2 [u_i(x) + A_ij(x) y_j(x) + u_i(x+h)]  */
    vv[i] = y[i] + 0.5*dx*(dydx[i] + uu[i]);
  }
  /*  y_i(x+h) = B^-1_ij(x+h) v_j(x,x+h)  */
  lineq( BB, vv, y, n );
  
  free(AA);
  return;
}


void TrapAlgSolve( double *y, double *wrkspc, double x, double dx, int n ){
  /* Implicit Trapezoidal algorithm for one step (dx) of integration: 
      dy_i/dx = f_i(x,y) 
       ==> 
      y_i(x+h) = y_i(x) + h/2 [f_i(x,y(x)) + f_i(x+h,y(x+h))]  */
  int i,j,ij;
  double *dydx,*uu,*vv,*BB,**AA;
  
  /* workspace pointers */
  dydx = wrkspc; uu = dydx + n; vv = uu + n; BB = vv + n;
  AA = (double **)malloc(n*sizeof(double *));
  AA[0] = BB + (n*n); for(i=1;i<n;i++) AA[i] = AA[0] + (i*n);
  
  
  
  free(AA);
  return;
}


void DiffAlgSolve( double *y, double *wrkspc, double x, double dx, int n , int ord ){
  /* Implicit RK algorithm for one step (dx) of integration: 
      dy_I/dx = f_I(x,y) 
       ==> 
      y_I(x+h) = y_I(x) + h Sum_j b_j k_I,j
      
                  /                                         \ 
      k_I,j = f_I | x + c_j h , y_J(x) + h Sum_l A_jl k_J,l | 
                  \                                         / 
      
      where A_jl , b_j , c_j are predetermined constants depending 
       on the order (Butcher table)  */
  int i,j,ij;
  double *ff;
  double **dfdy,**kk;
  
  /* workspace pointers */
  ff = wrkspc;
  dfdy = (double **)malloc(n*sizeof(double *));
  dfdy[0] = ff + n; for(i=1;i<n;i++) dfdy[i] = dfdy[0] + (i*n);
  kk = (double **)malloc(n*sizeof(double *));
  kk[0] = dfdy + (n*n); for(i=1;i<n;i++) kk[i] = kk[0] + (i*ord);
  
  DerivsJacob( ff, dfdy, x, y, n );
  
  
  free(dfdy);
  return;
}

