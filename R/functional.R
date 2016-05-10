library(Rcpp)
source('kdend.R')
sourceCpp("kdend.cpp")

# 1. Application of kdend: T(f), eg. rcd
# Need to estimate integral of nu(f(x))f(x)dx, and average nu(f(xi))

# nu(x)
nu=function(x){
  # res=0.5*abs(1-1/x)
  if(x==0){
    res=0
  }else{
    res=max(1-1/x,0)
    #res=log(x)
  }
  return(res)
}

# KDE based functional estimator: Multivariate rcd(X)
fkde=function(X,h,type="normal",cpp=T){
  X=as.matrix(X)
  n=dim(X)[1] # n: sample size
  if(missing(h)){
    h=rep(0.25*n^(-1/4),dim(X)[2])
  }
  X=apply(X,2,rank)/(n+1) # Copula transform
  nus=rep(0,n) # nu(f(xi))
  for(i in 1:n){
    if(cpp){
      fxi=kdendcpp(X[i,],X[-i,],h) 
    }else{
      fxi=kdend(X[i,],X[-i,],h,type) 
    }
    nus[i]=nu(fxi)
  }
  return(mean(nus))
}

# 2. Application of kdend: T(f,g,h), eg. conditional rcd
# rcd(X,Y|Z)
# Need to estimate mean(0.5*|c(u,w)c(v,w)/c(u,v,w)-1|)
# KDE based functional estimator: Conditional rcd(X,Y|Z)
cfkde=function(X,Y,Z,hx,hy,hz,type="normal",cpp=T){
  X=as.matrix(X)
  Y=as.matrix(Y)
  n=dim(X)[1] # n: sample size
  if(missing(hx)){
    hx=rep(0.25*n^(-1/4),dim(X)[2])
  }
  if(missing(hy)){
    hy=rep(0.25*n^(-1/4),dim(Y)[2])
  }
  
  if(missing(Z)){
    Z=0
    return(fkde(cbind(X,Y),c(hx,hy),type,cpp))
  }
  
  Z=as.matrix(Z)
  if(missing(hz)){
    hz=rep(0.25*n^(-1/4),dim(Z)[2])
  }
 
  
  U=apply(X,2,rank)/(n+1) # Copula transform
  V=apply(Y,2,rank)/(n+1) # Copula transform
  W=apply(Z,2,rank)/(n+1) # Copula transform
  nus=rep(0,n) # nu(f(xi))
  cuvwis=rep(0,n)
  for(i in 1:n){
    if(cpp){
      cuwi=kdendcpp(cbind(U,W)[i,],cbind(U,W)[-i,],c(hx,hz))
      cvwi=kdendcpp(cbind(V,W)[i,],cbind(V,W)[-i,],c(hy,hz))
      cuvwi=kdendcpp(cbind(U,V,W)[i,],cbind(U,V,W)[-i,],c(hx,hy,hz))
      cuvwis[i]=cuvwi
    }else{
      cuwi=kdend(cbind(U,W)[i,],cbind(U,W),c(hx,hz),type)
      cvwi=kdend(cbind(V,W)[i,],cbind(V,W),c(hy,hz),type)
      cuvwi=kdend(cbind(U,V,W)[i,],cbind(U,V,W),c(hx,hy,hz),type)
    }

    nus[i]=0.5*abs(cuwi*cvwi/cuvwi-1)
  }
  return(mean(nus))
}

###############################################
# Mutual Information Multivariate
# 3. Application of kdend: T(f), eg. rcd
# Need to estimate integral of nu(f(x))f(x)dx, and average nu(f(xi))

# nu(x)
numi=function(x){
  # res=0.5*abs(1-1/x)
  #res=max(1-1/x,0)
  res=log(x+.Machine$double.eps)
  return(res)
}

# KDE based functional estimator: Multivariate rcd(X)
fkdemi=function(X,h,type="normal",cpp=T){
  X=as.matrix(X)
  n=dim(X)[1] # n: sample size
  if(missing(h)){
    h=rep(0.25*n^(-1/4),dim(X)[2])
  }
  X=apply(X,2,rank)/(n+1) # Copula transform
  nus=rep(0,n) # nu(f(xi))
  for(i in 1:n){
    if(cpp){
      fxi=kdendcpp(X[i,],X[-i,],h) 
    }else{
      fxi=kdend(X[i,],X[-i,],h,type) 
    }
    nus[i]=numi(fxi)
  }
  return(mean(nus))
}
