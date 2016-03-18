#' Calculate robust copula dependence
#' 
#' This is the function that used to calculate the robust copula dependence (RCD) 
#' between two random variables \code{x} and \code{y}. Note that the length of 
#' \code{x} and \code{y} should be the same. 
#' 
#'  @param x The vector for the first random variable
#'  @param y The vector for the first random variable 
#'  @param method The method used in density estimation, i.e. knn or kde.
#'  @param bandwidth The bandwidth used in the density estimation. This parameter could be missing and a default value will be applied. 
#'  @param M Number of steps used in numerical integration. This parameter could be missing and a default value will be applied. 
#'  @param k The parameter K in KNN density estimation. This parameter could be missing and a default value will be applied. 
#'  @param cpp The ways to calculate the distance matrix. "none" means using R distance function, "serial" means using serial version of Rcpp, "parallel" means using parallel version of Rcpp.
#'  @return The RCD of \code{x} and \code{y}
#'  @examples
#'  n <- 1000
#'  x <- runif(n)
#'  y <- x^2 + 2*runif(n)
#'  rcd(x,y,method="kde")
#'  @export
rcd=function(x,y,method="knn",k,bandwidth,M=200,cpp="parallel"){  
  
  if(missing(y)){
    if(method=="kde"){
      warning("Single multivariate variable, using knn method.")
    }
    n=dim(x)[1]
    if(missing(k)){
      k=(1/4)*n^(4/(dim(x)[2]+6))
      k=round(max(1,k))
    }
    score=rcdknn_single(x,k,cpp=cpp)
    return(ifelse(score>1,1,ifelse(score<0,0,score)))
  }else{
    x=as.matrix(x)
    y=as.matrix(y)
    
    if(dim(x)[1]!=dim(y)[1]){
      stop("x and y must have the same length!")
    }
    
    if((dim(x)[2]>1 | dim(y)[2]>1) & method=="kde"){
      method="knn"
      warning("Multidimensional case, using method = knn.")
    }
    
    if(method=="knn"){
      n=dim(x)[1]
      if(missing(k)){
        k=(1/4)*n^(4/(dim(x)[2]+dim(y)[2]+6))
        k=round(max(1,k))
      }
      score=rcdknn(x,y,k,cpp=cpp)
      return(ifelse(score>1,1,ifelse(score<0,0,score)))
    }else if(method=="kde"){
      n=length(x)
      if(missing(bandwidth)){
        bandwidth=0.25*n^(-1/4)
      }
      score=rcdkde(x,y,bandwidth,M=M)
      return(ifelse(score>1,1,ifelse(score<0,0,score)))
    }else{
      stop("method is either kde or knn!")
    }
  }
  
 
}