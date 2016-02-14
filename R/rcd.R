#' Calculate robust copula dependence
#' 
#' This is the function that used to calculate the robust copula dependence (RCD) 
#' between two random variables \code{x} and \code{y}. Note that the length of 
#' \code{x} and \code{y} should be the same. 
#' 
#'  @param x The sampled data, a vector, for the first random variable
#'  @param y The sampled data, a vector, for the second random variable  
#'  @param k The bandwidth need to be used in the density estimation. This parameter could be missing and a deault value will be applied. 
#'  @param bw The ways to calculate the density
#'  @param method Using cpp or not
#'  @return The RCD of \code{x} and \code{y}
#'  @examples
#'  n <- 1000
#'  x <- runif(n)
#'  y <- x + runif(n)
#'  rcdknn(x,y)
#'  @export
rcd=function(x,y,method,k,bw){  
  
  if(method=="knn"){
    n=length(x)
    if(missing(bw)){
      k=0.25*sqrt(n)
    }
    score=rcdknn(x,y,cc=1/4,typ=0,cpp=2)
    return(score)
  }else if(method=="kde"){
    n=length(x)
    if(missing(bw)){
      bw=0.25*n^(-1/4)
    }
    score=rcdkde(x,y,bw,cpp=1)
    return(score)
  }else{
    stop("method is either kde or knn!")
  }
   
}