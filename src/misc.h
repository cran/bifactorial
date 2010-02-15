#include <vector>
using namespace std;
double mean(vector<double>& X);
double mean(vector<int>& X);
int gaussklammer(double x);
vector<double> bootstrap(vector<double>& X);
vector<int> binomial_sample(int n, double p);
double var(vector<double>& X);
double summe(vector<int>& X);
double min(double x, double y);
double max(double x, double y);
double maximalwert(vector<double>& X);
double minimalwert(vector<double>& X);
double betrag(double x);
int newnsim(RcppVector<int>& count, int k, int nsim,double simerror,int D);
int c2(int a,int b,int B);
int c20(int a,int b,int B);
int c3(int a,int b,int c,int B,int C);
int c30(int a,int b,int c,int B,int C);
