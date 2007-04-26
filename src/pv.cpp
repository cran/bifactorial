#include <iostream>
#include <math.h>
#include <time.h>
#include "Rcpp.hpp"
#include "misc.hpp"
#include <vector>
//Implementations of bootstrap algorithms for the min-test (Laska and Meisner, 1989)
//min-test based on Student's t-test in bifactorial designs
RcppExport SEXP student2(SEXP Yr,SEXP nr,SEXP tminr,SEXP parNr,SEXP parRr){
  SEXP rl=0; srand((unsigned) time(NULL)); char* exceptionMesg=NULL;
  RcppVector<int> parN(parNr); RcppVector<double> parR(parRr),Y(Yr);
  int a,b,i,j,k=1,m,nsim=parN(0),A=parN(1),B=parN(2),anfang[A+1][B+1];
  double simerror=parR(0),mx[A+1][B+1];
  RcppVector<int> count(A*B); RcppMatrix<double> tmin(tminr); RcppMatrix<int> n(nr);
  RcppResultSet rs; vector<double> tminst(A*B),X1(0),X2(0),Z1(0),Z2(0);
  for(a=0;a<=A;++a){for(b=0;b<=B;++b){
    anfang[a][b]=0; 
    for(i=0;i<=(a-1);++i){for(j=0;j<=B;++j){anfang[a][b]+=n(i,j);}}
    for(j=0;j<=b;++j){anfang[a][b]+=n(a,j);}
    anfang[a][b]-=n(a,b);
  }}
  for(i=0;i<=(A*B)-1;++i){count(i)=0;}
  for(a=0;a<=A;++a){for(b=0;b<=B;++b){
    X1.resize(0); for(m=0;m<=n(a,b)-1;++m){X1.push_back(Y(anfang[a][b]+m));} 
    mx[a][b]=mean(X1);
  }}
  while(k<=nsim){
    for(a=1;a<=A;++a){for(b=1;b<=B;++b){
      X1.resize(0); for(m=0;m<=n(a,b)-1;++m){X1.push_back(Y(anfang[a][b]+m)-mx[a][b]);} 
      Z1=bootstrap(X1); X2.resize(0);
      if(mx[a][0]-mx[0][b]>=0){
	for(m=0;m<=n(a,0)-1;++m){X2.push_back(Y(anfang[a][0]+m)-mx[a][0]);}
	Z2=bootstrap(X2);
	tminst[((a-1)*B)+b-1]=(mean(Z1)-mean(Z2))/sqrt((var(Z1)/n(a,b))+(var(Z2)/n(a,0)));
      }
      else{
	for(m=0;m<=n(0,b)-1;++m){X2.push_back(Y(anfang[0][b]+m)-mx[0][b]);}
	Z2=bootstrap(X2);
	tminst[((a-1)*B)+b-1]=(mean(Z1)-mean(Z2))/sqrt((var(Z1)/n(a,b))+(var(Z2)/n(0,b)));
      }
    }}
    for(a=1;a<=A;++a){for(b=1;b<=B;++b){
      if(maximalwert(tminst)>=tmin(a-1,b-1)){count(((a-1)*B)+b-1)+=1;}
    }}
    nsim=newnsim(count,k,nsim,simerror,A*B); ++k;
  }
  rs.add("count",count); rs.add("nsim",nsim); rl=rs.getReturnList(); 
  return rl;
}
//min-test based on Student's t-test in trifactorial designs
RcppExport SEXP student3(SEXP Yr,SEXP nr,SEXP tminr,SEXP parNr,SEXP parRr){
  SEXP rl=0; srand((unsigned) time(NULL)); char* exceptionMesg=NULL;
  RcppVector<int> parN(parNr); RcppVector<double> parR(parRr),Y(Yr);
  int a,b,c,i,j,k=1,m,nsim=parN(0),A=parN(1),B=parN(2),C=parN(3),anfang[A+1][B+1][C+1];
  double simerror=parR(0),mx[A+1][B+1][C+1];
  RcppVector<int> count(A*B*C); RcppVector<double> tmin(tminr); RcppVector<int> n(nr);
  RcppResultSet rs; vector<double> tminst(A*B*C),X1(0),X2(0),Z1(0),Z2(0);
  for(a=0;a<=A;++a){for(b=0;b<=B;++b){for(c=0;c<=C;++c){
    anfang[a][b][c]=0;
    for(i=0;i<=c3(a,b,c,B,C)-1;++i){anfang[a][b][c]+=n(i);}
  }}}
  for(i=0;i<=(A*B*C)-1;++i){count(i)=0;}
  for(a=0;a<=A;++a){for(b=0;b<=B;++b){for(c=0;c<=C;++c){
    X1.resize(0); for(m=0;m<=n(c3(a,b,c,B,C))-1;++m){X1.push_back(Y(anfang[a][b][c]+m));} 
    mx[a][b][c]=mean(X1);
  }}}
  while(k<=nsim){
    for(a=1;a<=A;++a){for(b=1;b<=B;++b){for(c=1;c<=C;++c){
      X1.resize(0);
      for(m=0;m<=n(c3(a,b,c,B,C))-1;++m){X1.push_back(Y(anfang[a][b][c]+m)-mx[a][b][c]);} 
      Z1=bootstrap(X1); X2.resize(0);
      if(mx[a][b][0]>=max(mx[0][b][c],mx[a][0][c])){
	for(m=0;m<=n(c3(a,b,0,B,C))-1;++m){X2.push_back(Y(anfang[a][b][0]+m)-mx[a][b][0]);}
	Z2=bootstrap(X2);
       tminst[c30(a,b,c,B,C)]=(mean(Z1)-mean(Z2))/sqrt((var(Z1)/n(c3(a,b,c,B,C)))+(var(Z2)/n(c3(a,b,0,B,C))));
      }
      if(mx[a][0][c]>=max(mx[0][b][c],mx[a][b][0])){
	for(m=0;m<=n(c3(a,b,c,B,C))-1;++m){X2.push_back(Y(anfang[a][0][c]+m)-mx[a][0][c]);}
	Z2=bootstrap(X2);
       tminst[c30(a,b,c,B,C)]=(mean(Z1)-mean(Z2))/sqrt((var(Z1)/n(c3(a,b,c,B,C)))+(var(Z2)/n(c3(a,0,c,B,C))));
      }
      if(mx[0][b][c]>=max(mx[a][b][0],mx[a][0][c])){
	for(m=0;m<=n(c3(0,b,c,B,C))-1;++m){X2.push_back(Y(anfang[0][b][c]+m)-mx[0][b][c]);}
	Z2=bootstrap(X2);
       tminst[c30(a,b,c,B,C)]=(mean(Z1)-mean(Z2))/sqrt((var(Z1)/n(c3(a,b,c,B,C)))+(var(Z2)/n(c3(0,b,c,B,C))));
      }
    }}}
    for(a=1;a<=A;++a){for(b=1;b<=B;++b){for(c=1;c<=C;++c){
      if(maximalwert(tminst)>=tmin(c30(a,b,c,B,C))){count(c30(a,b,c,B,C))+=1;}
    }}}
    nsim=newnsim(count,k,nsim,simerror,A*B*C); ++k;
  }
  rs.add("count",count); rs.add("nsim",nsim); rl=rs.getReturnList(); return rl;
}
//min-test based on a Z statistic in bifactorial designs
RcppExport SEXP binomial2(SEXP nr,SEXP pr,SEXP zminr,SEXP parNr,SEXP parRr){
  SEXP rl=0; srand((unsigned) time(NULL)); char* exceptionMesg=NULL;
  RcppVector<int> parN(parNr); RcppVector<double> parR(parRr); 
  RcppMatrix<int> n(nr); RcppMatrix<double> p(pr);
  int a,b,i,k=1,l,m,lmax=0,nsim=parN(0),A=parN(1),B=parN(2),anfang[A+1][B+1];
  double p1,p2,pi,simerror=parR(0);
  RcppVector<int> count(A*B); RcppMatrix<double> zmin(zminr); RcppResultSet rs;
  vector<double> zminst(A*B); vector<int> Z1(0),Z2(0);
  for(i=0;i<=(A*B)-1;++i){count(i)=0;}
  while(k<=nsim){
    for(a=1;a<=A;++a){for(b=1;b<=B;++b){
      Z1=binomial_sample(n(a,b),p(a,b)); p1=summe(Z1)/n(a,b);
      if(p(a,0)>=p(0,b)){
	pi=(p(a,b)+p(a,0))/2;
	Z1=binomial_sample(n(a,b),pi); p1=summe(Z1)/n(a,b);
	Z2=binomial_sample(n(a,0),pi); p2=summe(Z2)/n(a,0);
	zminst[((a-1)*B)+b-1]=(p1-p2)/sqrt(((p1*(1-p1))/n(a,b))+((p2*(1-p2))/n(a,0)));
      }
      else{
	pi=(p(a,b)+p(0,b))/2;
	Z1=binomial_sample(n(a,b),pi); p1=summe(Z1)/n(a,b);
	Z2=binomial_sample(n(0,b),pi); p2=summe(Z2)/n(0,b);
	zminst[((a-1)*B)+b-1]=(p1-p2)/sqrt(((p1*(1-p1))/n(a,b))+((p2*(1-p2))/n(0,b)));
      }
    }}
    for(a=1;a<=A;++a){for(b=1;b<=B;++b){
      if(maximalwert(zminst)>=zmin(a-1,b-1)){count(((a-1)*B)+b-1)+=1;}
    }}
    nsim=newnsim(count,k,nsim,simerror,A*B); ++k;
  }
  rs.add("count",count); rs.add("nsim",nsim); rl=rs.getReturnList(); return rl;
}
//min-test based on a Z statistic in trifactorial designs
RcppExport SEXP binomial3(SEXP nr,SEXP pr,SEXP zminr,SEXP parNr,SEXP parRr){
  SEXP rl=0; srand((unsigned) time(NULL)); char* exceptionMesg=NULL;
  RcppVector<int> parN(parNr); RcppVector<double> parR(parRr);
  int a,b,c,i,j,k=1,m,nsim=parN(0),A=parN(1),B=parN(2),C=parN(3);
  double p1,p2,pi,simerror=parR(0);
  RcppVector<double> p(pr),zmin(zminr); 
  RcppVector<int> count(A*B*C),n(nr);
  RcppResultSet rs; vector<double> zminst(A*B*C);
  vector<int> X1(0),X2(0),Z1(0),Z2(0);
  for(i=0;i<=(A*B*C)-1;++i){count(i)=0;}
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
  for(a=1;a<=A;++a){for(b=1;b<=B;++b){for(c=1;c<=C;++c){
    if(maximalwert(zminst)>=zmin(c30(a,b,c,B,C))){count(c30(a,b,c,B,C))+=1;}
  }}}
  nsim=newnsim(count,k,nsim,simerror,A*B*C); ++k;
  }
  rs.add("count",count); rs.add("nsim",nsim); rl=rs.getReturnList(); return rl;
}
