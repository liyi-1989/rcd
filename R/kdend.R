#' @useDynLib rcd
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs


# N-dim KDE estimation: used in ecdf+kde

# Kernel function in 1-dim
K=function(u,type="normal"){
  switch(type,
    rec={
      res=(abs(u)<=1)/2
    },
    tri={
      res=(1-abs(u))*(abs(u)<=1)
    },
    para={
      res=0.75*(1-u^2)*(abs(u)<=1)
    },
    normal={
      res=exp(-u^2)/sqrt(2*pi)
    },
    stop("No kernel in this type!")
    )
  return(res)
}

# Product Kernel for p-dim point
# x: Need to estimate density at x, f(x), where dim(x)=p
# xi: ith sample point. If sample is x1,..,xn, then this is xi
# h: bandwidth
PK=function(x,xi,h,type="normal"){
  p=length(xi) #p: dimension
  Kvec=rep(0,p) # simple kernel in each dimension
  for(i in 1:p){
    Kvec[i]=K((xi[i]-x[i])/h[i],type)/h[i]
    if(Kvec[i]==0){
      break
    }
  }
  res=prod(Kvec)
  return(res)
}

# N-dim KDE
# x: need to estimate density at point x, i.e. f(x)
# X: data matrix
# h: bandwidth 
kdend2=function(x,X,h,type="normal"){
  d=dim(X); n=d[1]; p=d[2] # n: sample size; p: dimension
  Kx=rep(0,n) 
  for(i in 1:n){
    #============ Product Kernel ============
    Kvec=rep(0,p) # simple kernel in each dimension
    for(j in 1:p){
      Kvec[j]=K((X[i,j]-x[j])/h[j],type)/h[j]
      if(Kvec[j]==0){
        break
      }
    }
    Kx[i]=prod(Kvec)#PK(x,X[i,],h,type)
    # ============ End of Product Kernel ============
  }
  return(mean(Kx))
}

kdendvec=function(X,h,type="rec"){
  d=dim(X); n=d[1]; p=d[2] # n: sample size; p: dimension
  cn=rep(0,n)
  for(irow in 1:n){
    Kx=rep(0,n)
    for(i in 1:n){
      #============ Product Kernel ============
      Kvec=rep(0,p) # simple kernel in each dimension
      for(j in 1:p){
        Kvec[j]=K((X[i,j]-X[irow,j])/h[j],type)/h[j]
        if(Kvec[j]==0){
          break
        }
      }
      Kx[i]=prod(Kvec)#PK(x,X[i,],h,type)
      # ============ End of Product Kernel ============
    }
    cn[irow]=mean(Kx)
  }
  return(cn)
}

# N-dim KDE
# x: need to estimate density at point x, i.e. f(x)
# X: data matrix
# h: bandwidth 
kdend=function(x,X,h,type="normal"){
  d=dim(X); n=d[1]; p=d[2] # n: sample size; p: dimension
  Kx=rep(0,n) 
  for(i in 1:n){
    Kx[i]=PK(x,X[i,],h,type)
  }
  return(mean(Kx))
}