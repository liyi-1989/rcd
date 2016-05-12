crcd.knn=function(x,y,z,k,cpp="parallel",verbose=F){
  if(verbose){print("Using knn for conditional rcd...")}
  n=nrow(x); dx=ncol(x); dy=ncol(y); dz=ncol(z);
  U=apply(x,2,rank)/(n+1); V=apply(y,2,rank)/(n+1); W=apply(z,2,rank)/(n+1);
  
  nus=rep(0,n)
  cuw=cvw=cuvw=rep(0,n)
  
  if(missing(k)){
    kuw=0.25*n^(2/(2+dx+dz))
    kvw=0.25*n^(2/(2+dy+dz))
    kuvw=0.25*n^(2/(2+dx+dy+dz))
    k=c(round(max(1,kuw)),round(max(1,kvw)),round(max(1,kuvw)))
  }
  
  cuw=dist2density(cbind(U,W),k[1],cpp)
  cvw=dist2density(cbind(V,W),k[2],cpp)
  cuvw=dist2density(cbind(U,V,W),k[3],cpp)
  
  nus=pmax(1-cuw*cvw/cuvw,0)
  return(mean(nus))
}