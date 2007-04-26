#include <iostream>
#include <math.h>
#include <time.h>
#include "Rcpp.hpp"
#include "misc.hpp"
#include <vector>
//Implementations of bootstrap algorithms for the AVE-test (Hung, 2000)
//AVE-test based on Student's t-test in bifactorial designs
RcppExport SEXP avestudent2(SEXP Yr,SEXP nr,SEXP sonstNr,SEXP sonstRr){
  RcppResultSet rs; RcppVector<int> sonstN(sonstNr); RcppMatrix<int> n(nr);
  int a,b,i,j,k=1,l,m,nsim=sonstN(0),A=sonstN(1),B=sonstN(2),anfang[A+1][B+1],count=0;
  RcppVector<double> sonstR(sonstRr),Y(Yr); vector<double> X(0),Z(0),tminst(A*B);
  double taveSc=sonstR(0),simerror=sonstR(1),mx[A+1][B+1],vx[A+1][B+1],oldmean[A+1][B+1],delta[A+1][B+1];
  for(a=0;a<=A;++a){for(b=0;b<=B;++b){
    anfang[a][b]=0; 
    for(i=0;i<=(a-1);++i){for(j=0;j<=B;++j){anfang[a][b]+=n(i,j);}}
    for(j=0;j<=b;++j){anfang[a][b]+=n(a,j);}
    anfang[a][b]-=n(a,b);
  }}
  for(a=0;a<=A;++a){for(b=0;b<=B;++b){
    for(m=0;m<=n(a,b)-1;++m){X.push_back(Y(anfang[a][b]+m));}
    oldmean[a][b]=mean(X); X.resize(0);
  }}
  while(k<=nsim){
    for(a=1;a<=A;++a){for(b=1;b<=B;++b){
      for(m=0;m<=n(a,b)-1;++m){X.push_back(Y(anfang[a][b]+m)-oldmean[a][b]);}
      Z=bootstrap(X); mx[a][b]=mean(Z); vx[a][b]=var(Z); X.resize(0);    
      if(oldmean[a][0]-oldmean[0][b]>=0){
	for(m=0;m<=n(a,0)-1;++m){X.push_back(Y(anfang[a][0]+m)-oldmean[a][0]);}
	tminst[(a-1)*B+b-1]=(mx[a][b]-mx[a][0])/sqrt((vx[a][b]/n(a,b))+(vx[a][0]/n(a,0)));
      }
      else{
	for(m=0;m<=n(0,b)-1;++m){X.push_back(Y(anfang[0][b]+m)-oldmean[0][b]);}
	tminst[(a-1)*B+b-1]=(mx[a][b]-mx[0][b])/sqrt((vx[a][b]/n(a,b))+(vx[0][b]/n(0,b)));
      }
    }}
    if((mean(tminst)/sqrt(var(tminst)/(A*B)))>=taveSc){count+=1;}
    if(k==900){if(simerror!=9){
      double p=double(count)/double(k);
      nsim=max(nsim,gaussklammer(double(((1-p)*p))/double(simerror*simerror))+100);
    }}
    ++k;
  }
  rs.add("count",count); rs.add("nsim",nsim); return rs.getReturnList();
}
//AVE-test based on Student's t-test in trifactorial designs
RcppExport SEXP avestudent3(SEXP Yr,SEXP nr,SEXP sonstNr,SEXP sonstRr){
  RcppResultSet rs; RcppVector<int> sonstN(sonstNr), n(nr); RcppVector<double> sonstR(sonstRr),Y(Yr);
  int a,b,c,i,j,k=1,l,m,nsim=sonstN(0),A=sonstN(1),B=sonstN(2),C=sonstN(3),anfang[A+1][B+1][C+1];
  double taveSc=sonstR(0),simerror=sonstR(1),mx[A+1][B+1][C+1],vx[A+1][B+1][C+1],count=0,oldmean[A+1][B+1][C+1];
  vector<double> X(0),Z(0),tminst(A*B*C);
  for(a=0;a<=A;++a){for(b=0;b<=B;++b){for(c=0;c<=C;++c){
    anfang[a][b][c]=0; 
    for(i=0;i<=c3(a,b,c,B,C)-1;++i){anfang[a][b][c]+=n(i);}
  }}}
  while(k<=nsim){
    for(a=0;a<=A;++a){for(b=0;b<=B;++b){for(c=0;c<=C;++c){
      X.resize(0);
      for(m=0;m<=n((a*(B+1)*(C+1))+(b*(C+1))+c)-1;++m){X.push_back(Y(anfang[a][b][c]+m));}
      oldmean[a][b][c]=mean(X);
    }}}
    for(a=1;a<=A;++a){for(b=1;b<=B;++b){for(c=0;c<=C;++c){
      X.resize(0);
      for(m=0;m<=n((a*(B+1)*(C+1))+(b*(C+1))+c)-1;++m){X.push_back(Y(anfang[a][b][c]+m)-oldmean[a][b][c]);}
      Z=bootstrap(X); mx[a][b][c]=mean(Z); vx[a][b][c]=var(Z); X.resize(0);
      if(oldmean[0][b][c]>=max(oldmean[a][0][c],oldmean[a][b][0])){
	for(m=0;m<=n((b*(C+1))+c)-1;++m){X.push_back(Y(anfang[0][b][c]+m)-oldmean[0][b][c]);}
        Z=bootstrap(X); mx[0][b][c]=mean(Z); vx[0][b][c]=var(Z);
	tminst[(a-1)*B+b-1]=(mx[a][b][c]-mx[0][b][c])/sqrt((vx[a][b][c]/n((a*(B+1)*(C+1))+(b*(C+1))+c))+(vx[0][b][c]/n((b*(C+1))+c)));
      }
      if(oldmean[a][0][c]>=max(oldmean[0][b][c],oldmean[a][b][0])){
	for(m=0;m<=n((a*(B+1)*(C+1))+c)-1;++m){X.push_back(Y(anfang[a][0][c]+m)-oldmean[a][0][c]);}
        Z=bootstrap(X); mx[a][0][c]=mean(Z); vx[a][0][c]=var(Z);
	tminst[(a-1)*B+b-1]=(mx[a][b][c]-mx[a][0][c])/sqrt((vx[a][b][c]/n((a*(B+1)*(C+1))+(b*(C+1))+c))+(vx[a][0][c]/n((a*(B+1)*(C+1))+c)));
      }
      if(oldmean[a][b][0]>=max(oldmean[a][0][c],oldmean[0][b][c])){
	for(m=0;m<=n((a*(B+1)*(C+1))+(b*(C+1)))-1;++m){X.push_back(Y(anfang[a][b][0]+m)-oldmean[a][b][0]);}
        Z=bootstrap(X); mx[a][b][0]=mean(Z); vx[a][b][0]=var(Z);
	tminst[(a-1)*B+b-1]=(mx[a][b][c]-mx[a][b][0])/sqrt((vx[a][b][c]/n((a*(B+1)*(C+1))+(b*(C+1))+c))+(vx[a][b][0]/n((a*(B+1)*(C+1))+(b*(C+1)))));
      }
    }}}
    if((mean(tminst)/sqrt(var(tminst)/(A*B*C)))>=taveSc){count+=1;}
    if(k==900){if(simerror!=9){
      double p=double(count)/double(k);
      nsim=max(nsim,gaussklammer(double(((1-p)*p))/double(simerror*simerror))+100);
    }}
    ++k;
  }
  rs.add("count",count); rs.add("nsim",nsim); return rs.getReturnList();
}
//AVE-test based on a Z statistic in bifactorial designs
RcppExport SEXP avebinomial2(SEXP nr,SEXP pr,SEXP sonstNr,SEXP sonstRr){
  RcppResultSet rs; RcppVector<int> sonstN(sonstNr); RcppMatrix<int> n(nr);
  RcppVector<double> sonstR(sonstRr); RcppMatrix<double> p(pr); 
  int a,b,i,j,k=1,l,m,nsim=sonstN(0),A=sonstN(1),B=sonstN(2);
  double zave=sonstR(0),simerror=sonstR(1),count=0,p1,p2;
  vector<double> zminst(A*B); vector<int> Z1(0),Z2(0);
  while(k<=nsim){
    zminst.resize(0);
    for(a=1;a<=A;++a){for(b=1;b<=B;++b){
      if(p(a,0)>=p(0,b)){
	Z1=binomial_sample(n(a,b),p(a,0)); p1=mean(Z1);
	Z2=binomial_sample(n(a,0),p(a,0)); p2=mean(Z2);
	zminst.push_back((p1-p2)/sqrt(((p1*(1-p1))/n(a,b))+((p2*(1-p2))/n(a,0))));
      }
      else{
	Z1=binomial_sample(n(a,b),p(0,b)); p1=mean(Z1);
	Z2=binomial_sample(n(0,b),p(0,b)); p2=mean(Z2);
	zminst.push_back((p1-p2)/sqrt(((p1*(1-p1))/n(a,b))+((p2*(1-p2))/n(0,b))));
      }
    }}
    if((mean(zminst)/sqrt(var(zminst)/(A*B)))>=zave){count+=1;}
    if(simerror!=9){if(k==900){
      double p=double(count)/double(k);
      nsim=max(nsim,gaussklammer(double(((1-p)*p))/double(simerror*simerror))+100);
    }}
    ++k;
  }
  rs.add("count",count); rs.add("nsim",nsim); return rs.getReturnList();
}
//AVE-test based on a Z statistic in trifactorial designs
RcppExport SEXP avebinomial3(SEXP nr,SEXP pr,SEXP sonstNr,SEXP sonstRr){
  RcppResultSet rs; RcppVector<int> sonstN(sonstNr), n(nr); RcppVector<double> sonstR(sonstRr), p(pr); 
  int a,b,c,i,j,k=1,l,m,nsim=sonstN(0),A=sonstN(1),B=sonstN(2),C=sonstN(3);
  double zave=sonstR(0),simerror=sonstR(1),count=0,p1,p2,pi;
  vector<double> zminst(A*B*C); vector<int> Z1(0),Z2(0);
  while(k<=nsim){
    for(a=1;a<=A;++a){for(b=1;b<=B;++b){for(c=1;c<=C;++c){
      if(p(c3(a,b,0,B,C))>=max(p(c3(a,0,c,B,C)),p(c3(0,b,c,B,C)))){
	pi=(p(c3(a,b,c,B,C))+p(c3(a,b,0,B,C)))/2;
	Z1=binomial_sample(n(c3(a,b,c,B,C)),pi); p1=mean(Z1);
	Z2=binomial_sample(n(c3(a,b,0,B,C)),pi); p2=mean(Z2);
	zminst[c30(a,b,c,B,C)]=(p1-p2)/sqrt(((p1*(1-p1))/n(c3(a,b,c,B,C)))+((p2*(1-p2))/n(c3(a,b,0,B,C))));
      }
      if(p(c3(a,0,c,B,C))>=max(p(c3(a,b,0,B,C)),p(c3(0,b,c,B,C)))){
	pi=(p(c3(a,b,c,B,C))+p(c3(a,0,c,B,C)))/2;
	Z1=binomial_sample(n(c3(a,b,c,B,C)),pi); p1=mean(Z1);
	Z2=binomial_sample(n(c3(a,0,c,B,C)),pi); p2=mean(Z2);
	zminst[c30(a,b,c,B,C)]=(p1-p2)/sqrt(((p1*(1-p1))/n(c3(a,b,c,B,C)))+((p2*(1-p2))/n(c3(a,0,c,B,C))));
      }
      if(p(c3(0,b,c,B,C))>=max(p(c3(a,b,0,B,C)),p(c3(a,0,c,B,C)))){
	pi=(p(c3(a,b,c,B,C))+p(c3(0,b,c,B,C)))/2;
	Z1=binomial_sample(n(c3(a,b,c,B,C)),pi); p1=mean(Z1);
	Z2=binomial_sample(n(c3(0,b,c,B,C)),pi); p2=mean(Z2);
	zminst[c30(a,b,c,B,C)]=(p1-p2)/sqrt(((p1*(1-p1))/n(c3(a,b,c,B,C)))+((p2*(1-p2))/n(c3(0,b,c,B,C))));
      }
    }}}
    if((mean(zminst)/sqrt(var(zminst)/(A*B*C)))>=zave){count+=1;}
    if(simerror!=9){if(k==900){
      double p=double(count)/double(k);
      nsim=max(nsim,gaussklammer(double(((1-p)*p))/double(simerror*simerror))+100);
    }}
    ++k;
  }
  rs.add("count",count); rs.add("nsim",nsim); return rs.getReturnList();
}
