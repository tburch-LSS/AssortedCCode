/*
  3d.c
*/

#include<math.h>
#include<stdio.h>
#include<stdlib.h>


main(){
  double ***matrix3,****matrix4,*s;
  int i,j,k,l,N=2,M=3,P=5,Q=7;

  matrix3 = (double ***)malloc(N*sizeof(double **));
  matrix3[0] = (double **)malloc(N*M*sizeof(double *));
  matrix3[0][0] = (double *)malloc(N*M*P*sizeof(double));

  for( i=1; i<N; i++ ) matrix3[i] = matrix3[0] + i*M;
  for( i=0; i<N; i++ ) for( j=0; j<M; j++ )
    matrix3[i][j] = matrix3[0][0] + i*M*P + j*P;

  for( i=0; i<N; i++ ) for( j=0; j<M; j++ )
    for( k=0, s=matrix3[i][j]; k<P; k++, s++ ) *s = 0.5*j;

  for( i=0; i<N; i++ ) for( j=0; j<M; j++ ) for( k=0; k<P; k++ )
    printf("m3[%2d][%2d][%2d] = %g\n",i,j,k,matrix3[i][j][k]);


  matrix4 = (double ****)malloc(N*sizeof(double ***));
  matrix4[0] = (double ***)malloc(N*M*sizeof(double **));
  matrix4[0][0] = (double **)malloc(N*M*P*sizeof(double *));
  matrix4[0][0][0] = (double *)malloc(N*M*P*Q*sizeof(double));

  for( i=1; i<N; i++ ) matrix4[i] = matrix4[0] + i*M;
  for( i=0; i<N; i++ ) for( j=0; j<M; j++ )
    matrix4[i][j] = matrix4[0][0] + i*M*P + j*P;
  for( i=0; i<N; i++ ) for( j=0; j<M; j++ ) for( k=0; k<P; k++ )
    matrix4[i][j][k] = matrix4[0][0][0] + i*M*P*Q + j*P*Q + k*Q;

  for( i=0; i<N; i++ ) for( j=0; j<M; j++ ) for( k=0; k<P; k++ )
    for( l=0, s=matrix4[i][j][k]; l<Q; l++, s++ ) *s = 0.5*k;

  for( i=0; i<N; i++ ) for( j=0; j<M; j++ ) for( k=0; k<P; k++ )
    for( l=0; l<Q; l++ )
      printf("m4[%2d][%2d][%2d][%2d] = %g\n",i,j,k,l,matrix4[i][j][k][l]);
}
