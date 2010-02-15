#include <iostream>
#include <math.h>
#include <time.h>
#include "Rcpp.h"
#include "misc.h"
#include <vector>
//Implementations of bootstrap algorithms for simultaneous confidence intervals
//Confidence intervals based on Student's t-test in bifactorial designs
RcppExport SEXP kritstudent2(SEXP Yr,SEXP nr,SEXP parNr){
  int a,b,i,k=0,m; SEXP rl=0; RcppResultSet rs; RcppVector<int> parN(parNr),n(nr);
  int nsim=parN(0),A=parN(1),B=parN(2);
  int ** anfang; anfang=new int *[A+1];
  for(i=0;i<A+1;++i){anfang[i]=new int[B+1];}
  RcppVector<double> Y(Yr),maxi(nsim),mini(nsim); 
  vector<double> X(0),Z(0),tst(A*B);
  double ** mx; mx=new double *[A+1]; for(i=0;i<A+1;++i){mx[i]=new double[B+1];}
  double ** vx; vx=new double *[A+1]; for(i=0;i<A+1;++i){vx[i]=new double[B+1];}
  for(a=0;a<=A;++a){for(b=0;b<=B;++b){
    anfang[a][b]=0; for(i=0;i<=c2(a,b,B)-1;++i){anfang[a][b]+=n(i);}
  }}
  srand((unsigned) time(NULL));
  while(k<=nsim-1){
    for(a=0;a<=A;++a){for(b=0;b<=B;++b){
      X.resize(0);
      for(m=0;m<=n(c2(a,b,B))-1;++m){X.push_back(Y(anfang[a][b]+m));}
      Z=bootstrap(X); mx[a][b]=mean(Z); vx[a][b]=var(Z);
    }}
    tst.resize(0);
    for(a=0;a<=A;++a){for(b=0;b<=B;++b){
      tst.push_back((mx[a][b]-mx[a][0])/sqrt((vx[a][b]/n(c2(a,b,B)))+(vx[a][0]/n(c2(a,0,B)))));
      tst.push_back((mx[a][b]-mx[0][b])/sqrt((vx[a][b]/n(c2(a,b,B)))+(vx[0][b]/n(c2(0,b,B)))));
    }}
    mini(k)=minimalwert(tst); maxi(k)=maximalwert(tst); ++k;
  }
  rs.add("mini",mini); rs.add("maxi",maxi); return rs.getReturnList();
}
//Confidence intervals based on Student's t-test in trifactorial designs
RcppExport SEXP kritstudent3(SEXP Yr,SEXP nr,SEXP parNr){
  int a,b,c,i,j,k=0,m; SEXP rl=0; RcppResultSet rs; RcppVector<int> parN(parNr),n(nr);
  int nsim=parN(0),A=parN(1),B=parN(2),C=parN(3);
  int *** anfang; anfang=new int **[A+1];
  for(i=0;i<A+1;++i){anfang[i]=new int *[B+1]; for(j=0;j<B+1;++j){anfang[i][j]=new int[C+1];}}
  RcppVector<double> Y(Yr),maxi(nsim),mini(nsim); 
  vector<double> X(0),Z(0),tst(A*B);
  double *** mx; mx=new double **[A+1];
  for(i=0;i<A+1;++i){mx[i]=new double *[B+1]; for(j=0;j<B+1;++j){mx[i][j]=new double[C+1];}}
  double *** vx; vx=new double **[A+1];
  for(i=0;i<A+1;++i){vx[i]=new double *[B+1]; for(j=0;j<B+1;++j){vx[i][j]=new double[C+1];}}
  for(a=0;a<=A;++a){for(b=0;b<=B;++b){for(c=0;c<=C;++c){
    anfang[a][b][c]=0; for(i=0;i<=c3(a,b,c,B,C)-1;++i){anfang[a][b][c]+=n(i);}
  }}}
  srand((unsigned) time(NULL));
  while(k<=nsim-1){
    for(a=0;a<=A;++a){for(b=0;b<=B;++b){for(c=0;c<=C;++c){
      X.resize(0);
      for(m=0;m<=n(c3(a,b,c,B,C))-1;++m){X.push_back(Y(anfang[a][b][c]+m));}
      Z=bootstrap(X); mx[a][b][c]=mean(Z); vx[a][b][c]=var(Z);
    }}}
    tst.resize(0);
    for(a=0;a<=A;++a){for(b=0;b<=B;++b){for(c=0;c<=C;++c){
      tst.push_back((mx[a][b][c]-mx[a][b][0])/sqrt((vx[a][b][c]/n(c3(a,b,c,B,C)))+(vx[a][b][0]/n(c3(a,b,0,B,C)))));
      tst.push_back((mx[a][b][c]-mx[a][0][c])/sqrt((vx[a][b][c]/n(c3(a,b,c,B,C)))+(vx[a][0][c]/n(c3(a,0,c,B,C)))));
      tst.push_back((mx[a][b][c]-mx[0][b][c])/sqrt((vx[a][b][c]/n(c3(a,b,c,B,C)))+(vx[0][b][c]/n(c3(0,b,c,B,C)))));
    }}}
    mini(k)=minimalwert(tst); maxi(k)=maximalwert(tst); ++k;
  }
  rs.add("mini",mini); rs.add("maxi",maxi); return rs.getReturnList();
}
//Confidence intervals based on a Z statistic in bifactorial designs
RcppExport SEXP kritbinomial2(SEXP pr,SEXP nr,SEXP parNr){
  int a,b,i,j,k=0,m; SEXP rl=0; RcppResultSet rs; RcppVector<int> parN(parNr),n(nr);
  int nsim=parN(0),A=parN(1),B=parN(2);
  RcppVector<double> p(pr),maxi(nsim),mini(nsim); 
  vector<int> X(0); vector<double> zst(A*B); double px[A+1][B+1],vx[A+1][B+1],pm;
  while(k<=nsim-1){
    zst.resize(0);
    for(a=1;a<=A;++a){for(b=1;b<=B;++b){
      pm=(p(c2(a,b,B))+p(c2(a,0,B)))/2;
      X=binomial_sample(n(c2(a,b,B)),pm); px[a][b]=mean(X); vx[a][b]=px[a][b]*(1-px[a][b]);
      X=binomial_sample(n(c2(a,0,B)),pm); px[a][0]=mean(X); vx[a][0]=px[a][0]*(1-px[a][0]);
      pm=(p(c2(a,b,B))+p(c2(0,b,B)))/2;
      X=binomial_sample(n(c2(a,b,B)),pm); px[a][b]=mean(X); vx[a][b]=px[a][b]*(1-px[a][b]);
      X=binomial_sample(n(c2(0,b,B)),pm); px[0][b]=mean(X); vx[0][b]=px[0][b]*(1-px[0][b]);
      zst.push_back((px[a][b]-px[a][0])/sqrt((vx[a][b]/n(c2(a,b,B)))+(vx[a][0]/n(c2(a,0,B)))));
      zst.push_back((px[a][b]-px[0][b])/sqrt((vx[a][b]/n(c2(a,b,B)))+(vx[0][b]/n(c2(0,b,B)))));
    }}
    mini(k)=minimalwert(zst); maxi(k)=maximalwert(zst); ++k;
  }
  rs.add("mini",mini); rs.add("maxi",maxi); return rs.getReturnList();
}
//Confidence intervals based on a Z statistic in trifactorial designs
RcppExport SEXP kritbinomial3(SEXP pr,SEXP nr,SEXP parNr){
  int a,b,c,i,j,k=0,m; RcppResultSet rs; RcppVector<int> parN(parNr),n(nr);
  int nsim=parN(0),A=parN(1),B=parN(2),C=parN(3);
  RcppVector<double> p(pr),maxi(nsim),mini(nsim); 
  vector<int> X(0); vector<double> zst(A*B*C); 
  double *** px; px=new double **[A+1];
  for(i=0;i<A+1;++i){px[i]=new double *[B+1]; for(j=0;j<B+1;++j){px[i][j]=new double[C+1];}}
  double *** vx; vx=new double **[A+1];
  for(i=0;i<A+1;++i){vx[i]=new double *[B+1]; for(j=0;j<B+1;++j){vx[i][j]=new double[C+1];}}
  double pm;
  while(k<=nsim-1){
    zst.resize(0);
    for(a=1;a<=A;++a){for(b=1;b<=B;++b){for(c=1;c<=C;++c){
      if(p(c3(a,b,0,B,C))>=max(p(c3(a,0,c,B,C)),p(c3(0,b,c,B,C)))){
	pm=(p(c3(a,b,c,B,C))+p(c3(a,b,0,B,C)))/2;
	X=binomial_sample(n(c3(a,b,c,B,C)),pm); px[a][b][c]=mean(X); vx[a][b][c]=px[a][b][c]*(1-px[a][b][c]);
	X=binomial_sample(n(c3(a,b,0,B,C)),pm); px[a][b][0]=mean(X); vx[a][b][0]=px[a][b][0]*(1-px[a][b][0]);
	zst.push_back((px[a][b][c]-px[a][b][0])/sqrt((vx[a][b][c]/n(c3(a,b,c,B,C)))+(vx[a][b][0]/n(c3(a,b,0,B,C)))));
	pm=(p(c3(a,b,c,B,C))+p(c3(a,0,c,B,C)))/2;
	X=binomial_sample(n(c3(a,b,c,B,C)),pm); px[a][b][c]=mean(X); vx[a][b][c]=px[a][b][c]*(1-px[a][b][c]);
	X=binomial_sample(n(c3(a,0,c,B,C)),pm); px[a][0][c]=mean(X); vx[a][0][c]=px[a][0][c]*(1-px[a][0][c]);
	zst.push_back((px[a][b][c]-px[a][0][c])/sqrt((vx[a][b][c]/n(c3(a,b,c,B,C)))+(vx[a][0][c]/n(c3(a,0,c,B,C)))));
	pm=(p(c3(a,b,c,B,C))+p(c3(0,b,c,B,C)))/2;
	X=binomial_sample(n(c3(a,b,c,B,C)),pm); px[a][b][c]=mean(X); vx[a][b][c]=px[a][b][c]*(1-px[a][b][c]);
	X=binomial_sample(n(c3(0,b,c,B,C)),pm); px[0][b][c]=mean(X); vx[0][b][c]=px[0][b][c]*(1-px[0][b][c]);
	zst.push_back((px[a][b][c]-px[0][b][c])/sqrt((vx[a][b][c]/n(c3(a,b,c,B,C)))+(vx[0][b][c]/n(c3(0,b,c,B,C)))));
      }
    }}}
    mini(k)=minimalwert(zst); maxi(k)=maximalwert(zst); ++k;
  }
  rs.add("mini",mini); rs.add("maxi",maxi); return rs.getReturnList();
}
