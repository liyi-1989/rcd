rcd.external.kde=function(x,y,bandwidth,cpp=T,S=F){
  # 1. if scaled, use scaled version
  if(S){return(rcd.external.kde.scale(x,y,bandwidth=bandwidth,cpp=cpp))}
  print("Using external kde ...")
  if(ncol(x)<ncol(y)){ # if dimension of y is larger, swap x and y.
    tmp=x;x=y;y=tmp;rm(tmp)
  }
  
  n=nrow(x); dx=ncol(x); dy=ncol(y); 
  U=apply(x,2,rank)/(n+1); V=apply(y,2,rank)/(n+1);

  nus=rep(0,n) # nu(f(xi))
  cu=cv=cuv=rep(0,n)
  
  if(dy==1){ # Case 1: dim(x)>dim(y)=1, special case, no need to estimate y
    
    if(missing(bandwidth)){
      hx=rep(0.25*n^(-1/(2+dx)),dx)
      hxy=rep(0.25*n^(-1/(2+dx+dy)),dx+dy)
      bandwidth=c(hx,hxy)
    }
    
    if(cpp!=FALSE){
      for(i in 1:n){
        cu[i]=kdendcpp(U[i,],U,bandwidth[1:dx])
        cuv[i]=kdendcpp(cbind(U,V)[i,],cbind(U,V),bandwidth[(dx+1):(dx+dx+dy)])
      }
    }else{
      for(i in 1:n){
        cu[i]=kdend(U[i,],U,bandwidth[1:dx],"rec")
        cuv[i]=kdend(cbind(U,V)[i,],cbind(U,V),bandwidth[(dx+dy+1):(dx+dy+dx+dy)],"rec")
      }
    }
    nus=0.5*abs(cu/cuv-1)
    return(mean(nus))
    
  }else{ # Case 2: dim(x)>dim(y)>1, general case
    
    if(missing(bandwidth)){
      hx=rep(0.25*n^(-1/(2+dx)),dx)
      hy=rep(0.25*n^(-1/(2+dy)),dy)
      hxy=rep(0.25*n^(-1/(2+dx+dy)),dx+dy)
      bandwidth=c(hx,hy,hxy)
    }
    
    if(cpp!=FALSE){
      for(i in 1:n){
        cu[i]=kdendcpp(U[i,],U,bandwidth[1:dx])
        cv[i]=kdendcpp(V[i,],V,bandwidth[(dx+1):(dx+dy)])
        cuv[i]=kdendcpp(cbind(U,V)[i,],cbind(U,V),bandwidth[(dx+dy+1):(dx+dy+dx+dy)])
      }
    }else{
      for(i in 1:n){
        cu[i]=kdend(U[i,],U,bandwidth[1:dx],"rec")
        cv[i]=kdend(V[i,],V,bandwidth[(dx+1):(dx+dy)],"rec")
        cuv[i]=kdend(cbind(U,V)[i,],cbind(U,V),bandwidth[(dx+dy+1):(dx+dy+dx+dy)],"rec")
      }
    }
    nus=0.5*abs(cu*cv/cuv-1)
    return(mean(nus))
  }
  
}


rcd.external.kde.scale=function(X,y,bandwidth,cpp=T){
  print("Using scaled external kde ...")
  n=nrow(X);dx=ncol(X);dy=ncol(y)
  if(missing(bandwidth)){
    bw=0.25*n^(-1/(2+dx))
  }else{
    bw=bandwidth[1]  
  }
  
  k=ceiling(2*bw*(n+1));x=floor(n/k);xin=1:n;yin=NULL
  for(i in 1:k){
    yin=cbind(yin,t(i+k*(0:x)))
  }

  maxc=rcd.external.kde(matrix(rep(xin,dx),n,dx),matrix(rep(xin,dy),n,dy),integral=integral,bandwidth=bandwidth,cpp=cpp,S=F)
  minc=rcd.internal.kde(as.matrix(cbind(xin,yin[yin<=n])),integral=integral,bandwidth=bandwidth,cpp=cpp,S=F)
  score=rcd.external.kde(X,y,integral=integral,bandwidth=bandwidth,cpp=cpp,S=F)
  r=(score-minc)/(maxc-minc)
  return(ifelse(r>0,r,score))
}






