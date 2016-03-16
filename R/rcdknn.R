#-----------------------------
# Ccor function with no scale
# k: k in knn
# typ=0: using the density(i) estimation from the paper(log(c) and then e^log(c))
# typ=1: using the density(i) estimation directly
#-----------------------------
ccorknne=function(x,y,k,cpp="parallel"){
  x=as.matrix(x)
  y=as.matrix(y)
  n=dim(x)[1]
  u=apply(x,2,rank)/(n+1)
  v=apply(y,2,rank)/(n+1)
  
#   k=cc*n^(4/(dim(x)[2]+dim(y)[2]+6))
#   k=round(max(1,k))
  
  if(cpp=="serial"){
    D=rcpp_distance(cbind(u,v))
  }else if(cpp=="parallel"){
    D=rcpp_parallel_distance(cbind(u,v))
  }else{
    df=data.frame(u,v) #data points in each row
    D=as.matrix(dist(df,method="maximum"))# distance matrix with maximum norm
  }
  
  e=rep(0,n)#e: epsilon(i) = 2*distance to the k nearest neighbour(knn)
  for(i in 1:n){
    e[i]=2*sort(D[i,])[k+1] # sort the distance matrix(col&row eighter is OK) and take the knn
  }
  d=dim(x)[2]+dim(y)[2]#d:dimension of (u,v)/Cd: vol of ball in dim d
  Cd=1#pi/4#(pi^(d/2))/(gamma(1+d/2)*(2^d))
# typ The ways to calculate the density. "0" means the density is calculated by using the result in the paper Kraskov(2003). "1" means the density is calculated by using the paper Loftsgaarden and Quesenberry (1965).
  cpl=exp(digamma(k)-digamma(n)-(d)*log(e)-log(Cd)) # copula density by using the paper
#   if(typ==0){
#     cpl=exp(digamma(k)-digamma(n)-(d)*log(e)-log(Cd)) # copula density by using the paper
#   }else if(typ==1){
#     cpl=k/(n*Cd*e^d) # copula density by calculating directly
#   }
  intg=1-1/cpl #ccor=0.5*\int |1-1/cpl|*cpl*dudv = AVG|1-1/cpl|/2 = AVG(|intg|)/2
  vp=mean(ifelse(intg>0,intg,0))#mean(max(intg,0)) # ccor = AVG(intg)_+
  return(vp)
}

#==========================================================================
minfc_knn=function(n,K,cpp=cpp){
  h=0.25*n^(-1/4)
  k=ceiling(2*h*(n+1))
  x=floor(n/k)
  xin=1:n
  yin=NULL
  for(i in 1:k){
    yin=cbind(yin,t(i+k*(0:x)))
  }
  yin=yin[yin<=n]
  r=ccorknne(xin,yin,K,cpp=cpp)
  #plot(xin,yin)
  return(r)
}

#' Calculate robust copula dependence with the KNN estimator
#' 
#' This is the function used to calculate the robust copula dependence (RCD) 
#' between two random variables \code{x} and \code{y} with the KNN estimator. Note that the length of 
#' \code{x} and \code{y} should be the same. 
#' 
#'  @param x The vector for the first random variable
#'  @param y The vector for the second random variable  
#'  @param k The parameter K in KNN that need to be used in the density estimation. This parameter could be missing and a deault value will be applied. 
#'  @param cpp The ways to calculate the distance matrix. "none" means using R distance function, "serial" means using serial version of Rcpp, "parallel" means using parallel version of Rcpp. This parameter could be missing and a default value will be applied.
#'  @return The RCD of \code{x} and \code{y}
#'  @examples
#'  n <- 1000
#'  x <- runif(n)
#'  y <- x^2 + 2*runif(n)
#'  rcdknn(x,y)
#'  @export
rcdknn=function(x,y,k,cpp="parallel"){  
  x=as.matrix(x)
  y=as.matrix(y)
  n=dim(x)[1]
  if(missing(k)){
    k=(1/4)*n^(4/(dim(x)[2]+dim(y)[2]+6))
    k=round(max(1,k))
  }
  maxc=ccorknne(1:n,1:n,k,cpp=cpp)
  minc=minfc_knn(n,k,cpp=cpp) 
  score=(ccorknne(x,y,k,cpp=cpp)-minc)/(maxc-minc)
  return(score)
}
