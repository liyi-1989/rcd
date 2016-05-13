#' @useDynLib rcd
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs



rcd.external.knn=function(x,y,k,cpp="parallel",S=F,verbose=F){
  # 1. if scaled, use scaled version
  if(S){return(rcd.external.knn.scale(x,y,k=k,cpp=cpp))}
  
  if((ncol(x)==1)&(ncol(y)==1)){
    return(rcd.internal.knn(cbind(x,y),k=k,cpp=cpp,S=S,verbose=verbose))
  }
  
  if(verbose){print("Using external knn ...")}
  
  if(ncol(x)<ncol(y)){ # if dimension of y is larger, swap x and y.
    tmp=x;x=y;y=tmp;rm(tmp)
  }
  
  n=nrow(x); dx=ncol(x); dy=ncol(y); 
  U=apply(x,2,rank)/(n+1); V=apply(y,2,rank)/(n+1);
  
  nus=rep(0,n) # nu(f(xi))
  cu=cv=cuv=rep(0,n)
  
  if(dy==1){ # Case 1: dim(x)>dim(y)=1, special case, no need to estimate y
    
    if(missing(k)){
      kx=0.25*n^(2/(2+dx))
      kxy=0.25*n^(2/(2+dx+dy))
      k=c(round(max(1,kx)),round(max(1,kxy)))
    }
    
    cu=dist2density(U,k[1],cpp)
    cuv=dist2density(cbind(U,V),k[2],cpp)

    nus=pmax(1-cu/cuv,0)#0.5*abs(cu/cuv-1)
    return(mean(nus))
    
  }else{ # Case 2: dim(x)>dim(y)>1, general case
    
    if(missing(k)){
      kx=0.25*n^(2/(2+dx))
      ky=0.25*n^(2/(2+dy))
      kxy=0.25*n^(2/(2+dx+dy))
      k=c(round(max(1,kx)),round(max(1,ky)),round(max(1,kxy)))
    }
    
    cu=dist2density(U,k[1],cpp)
    cv=dist2density(V,k[2],cpp)
    cuv=dist2density(cbind(U,V),k[3],cpp)

    nus=pmax(1-cu*cv/cuv,0)#0.5*abs(cu*cv/cuv-1)
    return(mean(nus))
  }
}


rcd.external.knn.scale=function(X,y,k=k,cpp="parallel",verbose=F){
  if(verbose){print("Using scaled external knn ...")}
  n=nrow(X);dx=ncol(X);dy=ncol(y)
  bw=0.25*n^(-1/(2+dx))
  kk=ceiling(2*bw*(n+1));x=floor(n/kk);xin=1:n;yin=NULL
  for(i in 1:kk){
    yin=cbind(yin,t(i+kk*(0:x)))
  }
  
  maxc=rcd.external.knn(matrix(rep(xin,dx),n,dx), matrix(rep(xin,dy),n,dy),k,cpp=cpp,S=F,verbose=F)
  minc=rcd.external.knn(as.matrix(xin),as.matrix(yin[yin<=n]),k,cpp=cpp,S=F,verbose=F)
  score=rcd.external.knn(X,y,k,cpp=cpp,S=F,verbose=F)
  r=(score-minc)/(maxc-minc)
  return(min(max(r,0),1))
}


dist2density=function(X,k,cpp){
  if(cpp=="serial"){
    D=rcpp_distance(X)
  }else if(cpp=="parallel"){
    D=rcpp_parallel_distance(X)
  }else{
    D=as.matrix(dist(data.frame(X),method="maximum"))# distance matrix with maximum norm
  }
  e=2*apply(D,1,sort)[k+1,]
  n=nrow(X);d=ncol(X) # d:dimension of (u,v)/Cd: vol of ball in dim d
  Cd=1#pi/4#(pi^(d/2))/(gamma(1+d/2)*(2^d))
  cpl=exp(digamma(k)-digamma(n)-(d)*log(e)-log(Cd)) # copula density by using the paper
  return(cpl)
}
