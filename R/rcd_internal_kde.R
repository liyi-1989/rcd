rcd.internal.kde=function(x,integral="ecdf",bandwidth,cpp=T,S=F,verbose=F){
  # 1. if scaled, use scaled version
  if(S){return(rcd.internal.kde.scale(x,integral=integral,bandwidth=bandwidth,cpp=cpp,verbose=verbose))}
  if(verbose){print("Using internal kde ...")}
  n=nrow(x)
  d=ncol(x)
  X=apply(x,2,rank)/(n+1)
  if(d==1){stop("X is only one dimensional!")}
  if(missing(bandwidth)){
    bandwidth=rep(0.25*n^(-1/(2+d)),d)
  }
  
  # 2. if not scaled 
  if((d>2)&(integral=="quad")){ # 2.0 Just checking quad is only for bivariate case 
    integral=="ecdf"
    warning("quad is only for bivariate x, switching to ecdf.")
  }
  
  if(integral=="quad"){ # 2.1 bivariate with quad 
    if(cpp!=FALSE){ # 2.1.1 bivariate with quad in C++  #########CALC HERE#########
      return(ccorecpp(X[,1],X[,2],bandwidth))
    }else{ # 2.1.2 bivariate with quad in R             #########CALC HERE#########
      # R code for bivariate kde with quad
      m=200
      u=X[,1]; v=X[,2]
      h=bandwidth[1]; l=bandwidth[2]
      A=matrix(0,m,m) #calculate the density on a m by m grid
      #pos = seq(m)/(m+1) pos = (2*seq(m)+1)/(2*m)
      for(k in 1:n){
        ul=max(1,floor(m*(u[k]-h)))
        uu=min(m,floor(m*(u[k]+h)))
        vl=max(1,floor(m*(v[k]-l)))
        vu=min(m,floor(m*(v[k]+l)))    
        A[ul:uu,vl:vu]=A[ul:uu,vl:vu]+1  
        
      }
      A=A/(n*h*l*4)
      return(mean(pmax(A-1,0)))
    }
  }else if(integral=="ecdf"){ # 2.2 multivariate with ecdf 
    X=as.matrix(X) # nus=rep(0,n) # nu(f(xi))
    if(cpp!=FALSE){                                     #########CALC HERE#########
      fxi=kdendveccpp(X,bandwidth)
    }else{                                              #########CALC HERE#########
      fxi=kdendvec(X,bandwidth,"rec")
    }
    nus=pmax(1-1/fxi,0) #nus=(fxi!=0)*pmax(1-1/fxi,0)
    # fxi=kdendcpp(X[i,],X,bandwidth);kdend(X[i,],X,bandwidth,"rec");nus[i]=ifelse(fxi==0,0,max(1-1/fxi,0)) #nus[i]=nu(fxi)
    return(mean(nus))
  }else{
   stop("The integration method is either ecdf or quad!") 
  }
}


rcd.internal.kde.scale=function(X,integral="ecdf",bandwidth,cpp=T,verbose=F){
  if(verbose){print("Using scaled internal kde ...")}
  n=nrow(X);dx=ncol(X)
  if(missing(bandwidth)){
    bw=0.25*n^(-1/(2+dx))
  }else{
    bw=bandwidth[1]
  }
  
  k=ceiling(2*bw*(n+1));x=floor(n/k);xin=1:n;yin=NULL
  for(i in 1:k){
    yin=cbind(yin,t(i+k*(0:x)))
  }
  
  maxc=rcd.internal.kde(matrix(rep(xin,dx),n,dx),integral=integral,bandwidth=bandwidth,cpp=cpp,S=F,verbose=F)
  minc=rcd.internal.kde(as.matrix(cbind(xin,yin[yin<=n])),integral=integral,bandwidth=bandwidth,cpp=cpp,S=F,verbose=F)
  score=rcd.internal.kde(X,integral=integral,bandwidth=bandwidth,cpp=cpp,S=F,verbose=F)
  r=(score-minc)/(maxc-minc)
  return(min(max(r,0),1))#return(ifelse(r>0,r,0))
}








