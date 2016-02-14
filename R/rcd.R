#' Calculate robust copula dependence
#' 
#' This is the function that used to calculate the robust copula dependence (RCD) 
#' between two random variables \code{x} and \code{y}. Note that the length of 
#' \code{x} and \code{y} should be the same. 
#' 
#'  @param x The sampled data, a vector, for the first random variable
#'  @param y The sampled data, a vector, for the second random variable  
#'  @param method The method used in density estimation, i.e. knn or kde.
#'  @param bw The bandwidth need to be used in the density estimation. This parameter could be missing and a deault value will be applied. 
#'  @param k The parameter K in KNN that need to be used in the density estimation. This parameter could be missing and a deault value will be applied. 
#'  @param typ The ways to calculate the density. "0" means the density is calculated by using the result in the paper Kraskov(2003). "1" means the density is calculated by using the paper Loftsgaarden and Quesenberry (1965).
#'  @param cpp The ways to calculate the distance matrix. "0" means using R distance function, "1" means using serial version of Rcpp, "2" means using parallel version of Rcpp.
#'  @return The RCD of \code{x} and \code{y}
#'  @examples
#'  n <- 1000
#'  x <- runif(n)
#'  y <- x + 0.5 * runif(n)
#'  rcd(x,y,method="kde")
#'  @export
rcd=function(x,y,method="knn",k,bw,typ=0,cpp=2){  
  
  if(method=="knn"){
    x=as.matrix(x)
    y=as.matrix(y)
    n=dim(x)[1]
    if(missing(k)){
      k=(1/4)*n^(4/(dim(x)[2]+dim(y)[2]+6))
      k=round(max(1,k))
    }
    score=rcdknn(x,y,k,typ=typ,cpp=cpp)
    return(score)
  }else if(method=="kde"){
    n=length(x)
    if(missing(bw)){
      bw=0.25*n^(-1/4)
    }
    score=rcdkde(x,y,bw)
    return(score)
  }else{
    stop("method is either kde or knn!")
  }
   
}