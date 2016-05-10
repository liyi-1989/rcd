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

rcd=function(x,y,type="internal",integral="ecdf",density="kde",k,bandwidth,cpp="parallel",S=F){
  if(missing(x)){
    # Case 1: X missing, impossible
    stop("X is missing!")
  }else if(missing(y)){
    # Case 2: single variable X (no Y) -> must be internal version
    if(type=="external") warning("Input is only X, switching to internal version.")
    if(density=="kde"){
      # Case 2.1: internal with kde
      return(rcd.internal.kde(x,integral=integral,bandwidth=bandwidth,cpp=cpp,S=S))
    }else if(density=="knn"){
      # Case 2.2: internal with knn
      return(rcd.internal.knn(x,k=k,cpp=cpp,S=S))
    }else{
      stop("Density must be kde or knn!")
    }
    
  }else{
    # Case 3: two variables X and Y -> must be external version
    if(type=="internal") warning("Input is X and Y, switching to external version.")
    if(density=="kde"){
      # Case 3.1: external with kde
      if((ncol(x)==1)&(ncol(y)==1)){
        # Case 3.1.1: external with kde but 2-dim -> using(coincide with) internal version
        return(rcd.internal.kde(cbind(x,y),integral=integral,bandwidth=bandwidth,cpp=cpp,S=S))
      }else{
        # Case 3.1.2: external with kde > 2-dim -> true multivariate version
        return(rcd.external.kde(x,y,bandwidth=bandwidth,cpp=cpp,S=S))
      }
    }else if(density=="knn"){
      # Case 3.2: external with knn
      return(rcd.external.knn(x,y,k=k,cpp=cpp,S=S))
    }else{
      stop("Density must be kde or knn!")
    }
  }
}

# 
# 
# rcd_ori=function(x,y,method="knn",k,bandwidth,M=200,cpp="parallel"){  
#   
#   if(missing(y)){
#     if(method=="kde"){
#       warning("Single multivariate variable, using knn method.")
#     }
#     n=dim(x)[1]
#     if(missing(k)){
#       k=(1/4)*n^(4/(dim(x)[2]+6))
#       k=round(max(1,k))
#     }
#     score=rcdknn_single(x,k,cpp=cpp)
#     return(ifelse(score>1,1,ifelse(score<0,0,score)))
#   }else{
#     x=as.matrix(x)
#     y=as.matrix(y)
#     
#     if(dim(x)[1]!=dim(y)[1]){
#       stop("x and y must have the same length!")
#     }
#     
#     if((dim(x)[2]>1 | dim(y)[2]>1) & method=="kde"){
#       method="knn"
#       warning("Multidimensional case, using method = knn.")
#     }
#     
#     if(method=="knn"){
#       n=dim(x)[1]
#       if(missing(k)){
#         k=(1/4)*n^(4/(dim(x)[2]+dim(y)[2]+6))
#         k=round(max(1,k))
#       }
#       score=rcdknn(x,y,k,cpp=cpp)
#       return(ifelse(score>1,1,ifelse(score<0,0,score)))
#     }else if(method=="kde"){
#       n=length(x)
#       if(missing(bandwidth)){
#         bandwidth=0.25*n^(-1/4)
#       }
#       score=rcdkde(x,y,bandwidth,M=M)
#       return(ifelse(score>1,1,ifelse(score<0,0,score)))
#     }else{
#       stop("method is either kde or knn!")
#     }
#   }
#   
#  
# }