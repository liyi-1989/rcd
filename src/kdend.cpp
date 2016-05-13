#include <Rcpp.h>
#include <cmath>
// #include <vector>
using namespace Rcpp;


double Kcpp(double u){  
  double res = (std::abs(u)<=1)/2.0;
  // double res = exp(-pow(u,2)/2)/sqrt(2*M_PI);
  return res;
}



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
  int n=X.nrow(); // sample size; int p=X.ncol(); // dimension
  double Kx=0;
  NumericVector xi;
  for(int i=0;i<n;i++){
    xi=X(i,_);
    Kx=Kx+PKcpp(x,xi,h);
  }
  return Kx/n;
}



// double kdendcpp2(NumericVector x, NumericMatrix X, NumericVector h){
//   int n=X.nrow(),p=X.ncol(); // sample size & dimension
//   double Kx=0,Kvec = 1.0;
//   // NumericVector xi;
//   for(int i=0;i<n;i++){
//     //xi=X(i,_);
//     // =============== PKcpp(Product cpp) =============== 
//     Kvec = 1.0;
//     for(int j=0;j<p;j++){
//       // Kvec=Kvec*Kcpp((xi[j]-x[j])/h[j])/h[j]; (xi[j]-x[j])/h[j]
//       Kvec=Kvec*(std::abs((X(i,j)-x[j])/h[j])<=1)/(2*h[j]); // using square kernel
//       if(Kvec==0){
//         break;
//       }
//     }
//     // =============== End PKcpp(return the product Kvec) ===============
//     Kx=Kx+Kvec;
//   }
//   return Kx/n;
// }
// 

// NumericVector kdendveccpp2(NumericMatrix X, NumericVector h){
//   int n=X.nrow(); // sample size
//   NumericVector cn=X(_,1);
//   for(int i=0;i<n;i++){
//     cn(i)=kdendcpp2(X(i,_),X,h);
//   }
//   return cn; 
// }


// [[Rcpp::export]]
NumericVector kdendveccpp(NumericMatrix X, NumericVector h){
  int n=X.nrow(),p=X.ncol(); // sample size & dimension
  double Kx=0,Kvec = 1.0;
  NumericVector cn=X(_,1);
  for(int irow=0;irow<n;irow++){
    Kx=0;
    // %%%%%%%%%%%%%%% Density for ith obs %%%%%%%%%%%%%%%
    for(int i=0;i<n;i++){
      //xi=X(i,_);
      // ============== PKcpp(Product cpp) =============== 
      Kvec = 1.0;
      for(int j=0;j<p;j++){
        // Kvec=Kvec*Kcpp((xi[j]-x[j])/h[j])/h[j]; (xi[j]-x[j])/h[j]
        Kvec=Kvec*((std::abs((X(i,j)-X(irow,j))/h[j])<=1)/(2.0*h[j])); // using square kernel
        if(Kvec==0){
          break;
        }
      }
      // ====== End PKcpp(return the product Kvec) =======
      Kx=Kx+Kvec;
    }
    cn(irow)=Kx/n;
    // %%%%%%%%%%%%% End density for ith obs %%%%%%%%%%%%%
  }
  return cn; 
}















