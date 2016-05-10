#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double Kcpp(double u){  
  // double res = (std::abs(u)<=1)/2.0;
  double res = exp(-pow(u,2)/2)/sqrt(2*M_PI);
  return res;
}

// [[Rcpp::export]]
double PKcpp(NumericVector x, NumericVector xi, NumericVector h){
  int p = xi.size();
  double Kvec = 1.0;
  for(int i=0;i<p;i++){
    Kvec=Kvec*Kcpp((xi[i]-x[i])/h[i])/h[i];
    if(Kvec==0){
      break;
    }
  }
  return Kvec;
}

// [[Rcpp::export]]
double kdendcpp(NumericVector x, NumericMatrix X, NumericVector h){
  int n=X.nrow(); // sample size
  int p=X.ncol(); // dimension
  double Kx=0;
  NumericVector xi;
  for(int i=0;i<n;i++){
    xi=X(i,_);
    Kx=Kx+PKcpp(x,xi,h);
  }
  return Kx/n;
}



