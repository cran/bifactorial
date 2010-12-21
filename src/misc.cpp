#include <math.h>
#include <stdlib.h>
#include <vector> 
using namespace std;
#include <iostream>
#include <valarray> 
using std::valarray;
#include <time.h>
#include <Rcpp.h>
//Determine largest integer smaller than x///////////////////////////
int gaussklammer(double x){
  int k=0;
  while(k<=x){++k;}
  return k-1;
}
//Draw a binomial distributed sample/////////////////////////////////
vector<int> binomial_sample(int n, double p){
  double * U; U=new double[n];
  vector<int> Y(n);
  for(int m=0;m<=n-1;++m){
    U[m]=(1+(rand()%10000))*.0001;
    while(U[m]>=1){U[m]=(1+(rand()%10000))*.0001;}
    if(U[m]<=p){Y[m]=1;}
    else{Y[m]=0;}
  }
  return Y;
}
//Draw a bootstrap sample of X///////////////////////////////////////
vector<double> bootstrap(vector<double>& X){
  int m, n=X.size();
  double * U; U=new double[n];
  int * index; index=new int[n];
  vector<double> Y(n);
  for(m=0;m<=n-1;++m){
    U[m]=(1+(rand()%10000))*.0001;
    while(U[m]>=1){U[m]=(1+(rand()%10000))*.0001;}
    index[m]=gaussklammer(n*U[m]);
    Y[m]=X[index[m]];
  }
  return Y;
}
//Calculate mean of metric data X////////////////////////////////////
double mean(vector<double>& X){
  int i, n=X.size();
  double summe=0;
  for(i=1;i<=n;++i){
    summe+=X[i-1];
  }
  return summe/n;
}
//Calculate mean of binary data X////////////////////////////////////
double mean(vector<int>& X){
  int i, n=X.size();
  int summe=0;
  for(i=1;i<=n;++i){
    summe+=X[i-1];
    //cout<<X[i-1]<<endl;
  }
  //    cout<<summe<<endl;
  return double(summe)/double(n);
}
//Calculate variance of metric data X////////////////////////////////
double var(vector<double>& X){
  int i, n=X.size();
  double varianz=0;
  for(i=0;i<=n-1;++i){
    varianz+=(X[i]-mean(X))*(X[i]-mean(X));
  }
  return varianz/(n-1);
}
//Calculate sum of binary data X/////////////////////////////////////
double summe(vector<int>& X){
  int i,n=X.size();
  double summe=0;
  for(i=0;i<=n-1;++i){summe+=X[i];}
  return summe;
}
//Determine the minimum of x and y///////////////////////////////////
double min(double x, double y){
  if(x<=y){ return x; }
  return y;
}
//Determine the maximum of x and y///////////////////////////////////
double max(double x, double y){
  if(x>=y){ return x; }
  return y;
}
//Determine the minimum of the data X////////////////////////////////
double minimalwert(vector<double>& X){
  int n=X.size();
  double minimum=8000;
  for(int i=0;i<=n-1;++i){
    if(X[i]<minimum){minimum=X[i];}
  }
  return minimum;
}
//Determine the maximum of the data X////////////////////////////////
double maximalwert(vector<double>& X){
  int n=X.size();
  double maximum=-100000;
  for(int i=0;i<=n-1;++i){
    if(X[i]>maximum){maximum=X[i];}
  }
  return maximum;
}
//Determine absolute value of x//////////////////////////////////////
double betrag(double x){
  if(x>=0){return x;}
  else{return -x;}
}
//Calculate number of simulations according to prespecified standard error
int newnsim(Rcpp::IntegerVector& count,int k, int nsim,double simerror,int D){
  vector<double> abst(D);
  int l,lmax=0;
  if(simerror!=9 & k==900){
    for(l=0;l<=D-1;++l){abst[l]=betrag(double(count(l))-double(450));}
    for(l=0;l<=D-1;++l){if(abst[l]==minimalwert(abst)){lmax=l;}}//Der Fehler wird für p=.5 maximal
    double p=double(count(lmax))/double(k);
    nsim=max(nsim,gaussklammer(double(((1-p)*p))/double(simerror*simerror))+100);
  }
  return nsim;
}
//Determine position in a vector for k=2/////////////////////////////
int c2(int a,int b,int B){
  return (a*(B+1))+b;
}
//Determine position in a reduced vector for k=2/////////////////////
int c20(int a,int b,int B){
  return ((a-1)*B)+b-1;
}

//Determine position in a vector for k=3/////////////////////////////
int c3(int a,int b,int c,int B,int C){
  return (a*(B+1)*(C+1))+(b*(C+1))+c;
}

//Determine position in a reduced vector for k=3/////////////////////
int c30(int a, int b,int c,int B,int C){
  return ((a-1)*B*C)+((b-1)*C)+c-1;
}
