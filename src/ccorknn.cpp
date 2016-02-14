#include <Rcpp.h>
using namespace Rcpp;
#include <cmath>
#include <algorithm>
/* 
S1. generic function for divergence
*/
template <typename InputIterator1, typename InputIterator2>
inline double divergence(InputIterator1 begin1, InputIterator1 end1, InputIterator2 begin2) {    
   double rval = 0; // value to return
   double rele = 0;     
   InputIterator1 it1 = begin1; // set iterators to beginning of ranges
   InputIterator2 it2 = begin2;      
   while (it1 != end1) { // for each input item         
      double d1 = *it1++; // take the value and increment the iterator
      double d2 = *it2++;
      //rval += std::pow(d1-d2,2);
      rele = std::abs(d1-d2);
      rval = (rval>rele)?rval:rele;
   }
   //return std::sqrt(rval);  
   return rval;
}


/* 
S2. serial version of distance matrix
*/

// [[Rcpp::export]]
NumericMatrix rcpp_distance(NumericMatrix mat) {  
   NumericMatrix rmat(mat.nrow(), mat.nrow()); // allocate the matrix we will return   
   for (int i = 0; i < rmat.nrow(); i++) {
      for (int j = 0; j < i; j++) {
         NumericMatrix::Row row1 = mat.row(i); // rows we will operate on
         NumericMatrix::Row row2 = mat.row(j);
         // write to output matrix
         rmat(i,j) = divergence(row1.begin(), row1.end(), row2.begin());
         rmat(j,i) = rmat(i,j);
      }
   } 
   return rmat;
}

/* 
S3. Parallel version of distance matrix
*/

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

struct MyDistance : public Worker {
   
   // input matrix to read from
   const RMatrix<double> mat;
   
   // output matrix to write to
   RMatrix<double> rmat;
   
   // initialize from Rcpp input and output matrixes (the RMatrix class
   // can be automatically converted to from the Rcpp matrix type)
   MyDistance(const NumericMatrix mat, NumericMatrix rmat)
      : mat(mat), rmat(rmat) {}
   
   // function call operator that work for the specified range (begin/end)
   void operator()(std::size_t begin, std::size_t end) {
      for (std::size_t i = begin; i < end; i++) {
         for (std::size_t j = 0; j < i; j++) {
            
            // rows we will operate on
            RMatrix<double>::Row row1 = mat.row(i);
            RMatrix<double>::Row row2 = mat.row(j);

            rmat(i,j) = divergence(row1.begin(), row1.end(), row2.begin());
            rmat(j,i) = rmat(i,j);
         
         }
      }
   }
};


// [[Rcpp::export]]
NumericMatrix rcpp_parallel_distance(NumericMatrix mat) {
  
   // allocate the matrix we will return
   NumericMatrix rmat(mat.nrow(), mat.nrow());

   // create the worker
   MyDistance myDistance(mat, rmat);
     
   // call it with parallelFor
   parallelFor(0, mat.nrow(), myDistance);

   return rmat;
}