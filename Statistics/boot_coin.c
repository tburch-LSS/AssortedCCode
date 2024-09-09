#include <math.h>
#include <stdio.h>

#define N 100
#define p 0.05
#define Nboot 1000
#define BOOT

void srand48(int seed);
double drand48(void);

main(){
  int seed,i,j,k,dist[N];
  double avg[Nboot],Avg,Low,High,fract_area,x;

  fprintf(stderr,"Enter seed\n");
  scanf("%d",&seed);
  srand48(seed);

#ifdef BOOT
  for ( i=0; i<N; i++ ) {
    if ( drand48() < p ) dist[i]=1;
    else dist[i]=0;
  }

  Avg=0.0;
  for ( j=0; j<Nboot; j++ ) {
    avg[j]=0.0;
    for ( i=0; i<N; i++ ) {
      k=N*drand48();
      avg[j]+=1.0*dist[k];
    }
    avg[j]/=1.0*N;
    printf("%le\n",avg[j]);
    Avg+=avg[j];
  }
  Avg/=1.0*Nboot;

  /**
  fract_area=0.0;
  x=0.0;
  while ( fract_area < 0.1585 ) {
    x+=1.0/N;
    for ( j=0; j<Nboot; j++ ) if ( avg[j] < x ) fract_area+=1.0/Nboot;
  }
  Low=x;

  fract_area=0.0;
  x=1.0;
  while ( fract_area < 0.1585 ) {
    x-=1.0/N;
    for ( j=0; j<Nboot; j++ ) if ( avg[j] > x ) fract_area+=1.0/Nboot;
  }
  High=x;

  Low=Avg-Low;
  High=High-Avg;
  printf("Avg = %le + %le - %le\n",Avg,High,Low);
  **/
#endif

#ifdef ACTUAL
  for ( j=0; j<Nboot; j++ ) {
    avg[j]=0.0;
    for ( i=0; i<N; i++ ) if ( drand48() < p ) avg[j]+=1.0;
    avg[j]/=1.0*N;
    printf("%le\n",avg[j]);
  }
#endif
}
