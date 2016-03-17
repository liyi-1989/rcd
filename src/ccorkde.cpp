#include <Rcpp.h>
#include <cmath>
// #define M 200
using namespace Rcpp;

// [[Rcpp::export]]

double ccorecpp(NumericVector u, NumericVector v, NumericVector bw, int M){
  // Calculate Ccor using positive part, |1-c|_+
  double n=u.size();
  //double h=0.25*pow(n,-1.0/4.0);
  double h=bw[0];
  double A[M][M];
  int i(0),j(0),k(0),ul(0),uu(0),vl(0),vu(0),ir(0),ic(0);  
  double s(0);
  
  for(k=0;k<n;k++){
    ul=std::max(std::floor((u[k]-h)*M),0.0);
    uu=std::min(std::floor((u[k]+h)*M),(double)M-1);
    vl=std::max(std::floor((v[k]-h)*M),0.0);
    vu=std::min(std::floor((v[k]+h)*M),(double)M-1);
    
    for(ir=ul;ir<=uu;ir++){
      for(ic=vl;ic<=vu;ic++){
        A[ir][ic]++;
      }
    }
    
  }
    
  for(i=0;i<M;i++){
    for(j=0;j<M;j++){
      A[i][j]=std::max(A[i][j]/(n*h*h*4)-1,0.0);
      s+=A[i][j];
    }
  }
  

  return(s/(M*M));
}

// [[Rcpp::export]]

double ccorecppabs(NumericVector u, NumericVector v, NumericVector bw, int M){
  double n=u.size();
  double h=bw[0];
  //double h=0.25*pow(n,-1.0/4.0);
  double A[M][M];
  int i(0),j(0),k(0),ul(0),uu(0),vl(0),vu(0),ir(0),ic(0);  
  double s(0);
  
  for(k=0;k<n;k++){
    ul=std::max(std::floor((u[k]-h)*M),0.0);
    uu=std::min(std::floor((u[k]+h)*M),(double)M-1);
    vl=std::max(std::floor((v[k]-h)*M),0.0);
    vu=std::min(std::floor((v[k]+h)*M),(double)M-1);
    
    for(ir=ul;ir<=uu;ir++){
      for(ic=vl;ic<=vu;ic++){
        A[ir][ic]++;
      }
    }
    
  }
    
  for(i=0;i<M;i++){
    for(j=0;j<M;j++){
      A[i][j]=std::abs(1-A[i][j]/(n*h*h*4));
      s+=A[i][j];
    }
  }
  

  return(s/(2*M*M));
}
