#' @useDynLib rcd
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs

ccorercpp=function(x,y,bw,M){
  n=length(x)
  u=rank(x)/(n+1)
  v=rank(y)/(n+1)
  return(ccorecpp(u,v,bw,M))
}

minfc.p=function(n,bw,M){
  k=ceiling(2*bw*(n+1))
  x=floor(n/k)
  xin=1:n
  yin=NULL
  for(i in 1:k){
    yin=cbind(yin,t(i+k*(0:x)))
  }
  return(ccorercpp(xin,yin[yin<=n],bw=bw,M=M))
}

#' Calculate robust copula dependence with the KDE estimator
#' 
#' This is the function used to calculate the robust copula dependence (RCD) 
#' between two random variables \code{x} and \code{y} with the KDE estimator. Note that the length of 
#' \code{x} and \code{y} should be the same. 
#' 
#'  @param x The vector for the first random variable
#'  @param y The vector for the second random variable  
#'  @param bandwidth The bandwidth used in the density estimation. This parameter could be missing and a default value will be applied. 
#'  @param M Number of steps used in numerical integration. This parameter could be missing and a default value will be applied.
#'  @return The RCD of \code{x} and \code{y}
#'  @examples
#'  n <- 1000
#'  x <- runif(n)
#'  y <- x^2 + 2*runif(n)
#'  rcdkde(x,y)
#'  @export
rcdkde=function(x,y,bandwidth,M){ 
  n=length(x)
  if(n!=length(y)){
    stop("x and y must have the same length!")
  }
  if(missing(bandwidth)){
    bandwidth=0.25*n^(-1/4)
  }
  if(missing(M)){
    M=200
  }
  maxc=ccorercpp(1:n,1:n,bw=bandwidth,M=M)
  minc=minfc.p(n,bandwidth,M) 
  score=(ccorercpp(x,y,bw=bandwidth,M=M)-minc)/(maxc-minc)
  return(score)
}
