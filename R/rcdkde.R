#' @useDynLib rcd
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs

ccorercpp=function(x,y,bw,cpp=1){
  n=length(x)
  u=rank(x)/(n+1)
  v=rank(y)/(n+1)
  
  
  if(cpp==1){
    return(ccorecpp(u,v,bw))
  }else{
    h=l=bw
    m=200
    A=matrix(0,m,m) #calculate the density on a m by m grid
    pos = seq(m)/(m+1)
    for(k in 1:n){
      ul=max(1,floor(m*(u[k]-h)))
      uu=min(m,floor(m*(u[k]+h)))
      vl=max(1,floor(m*(v[k]-l)))
      vu=min(m,floor(m*(v[k]+l)))    
      A[ul:uu,vl:vu]=A[ul:uu,vl:vu]+1  
      
    }
    A=A/(n*h*l*4)
    return(mean(pmax(1-A,0)))
  }
  
  
}

minfc.p=function(n,bw,cpp=1){
  h=bw
  k=ceiling(2*h*(n+1))
  x=floor(n/k)
  xin=1:n
  yin=NULL
  for(i in 1:k){
    yin=cbind(yin,t(i+k*(0:x)))
  }
  return(ccorercpp(xin,yin[yin<=n],bw=bw,cpp=cpp))
}

#' Calculate robust copula dependence
#' 
#' This is the function that used to calculate the robust copula dependence (RCD) 
#' between two random variables \code{x} and \code{y}. Note that the length of 
#' \code{x} and \code{y} should be the same. 
#' 
#'  @param x The sampled data, a vector, for the first random variable
#'  @param y The sampled data, a vector, for the second random variable  
#'  @param bw The bandwidth need to be used in the density estimation. This parameter could be missing and a deault value will be applied. 
#'  @param cpp use cpp or not
#'  @return The RCD of \code{x} and \code{y}
#'  @examples
#'  n <- 1000
#'  x <- runif(n)
#'  y <- x + runif(n)
#'  rcdkde(x,y)
#'  @export
rcdkde=function(x,y,bw,cpp=1){
  n=length(x)
  if(missing(bw)){
    bw=0.25*n^(-1/4)
  }
  maxc=ccorercpp(1:n,1:n,bw=bw,cpp=cpp)
  minc=minfc.p(n,bw) 
  score=(ccorercpp(x,y,bw=bw,cpp=cpp)-minc)/(maxc-minc)
  return(score)
}
