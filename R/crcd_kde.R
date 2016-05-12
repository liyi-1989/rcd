crcd.kde=function(x,y,z,bandwidth,cpp=T,verbose=F){
  if(verbose){print("Using kde for conditional rcd...")}
  n=nrow(x); dx=ncol(x); dy=ncol(y); dz=ncol(z);
  U=apply(x,2,rank)/(n+1); V=apply(y,2,rank)/(n+1); W=apply(z,2,rank)/(n+1);
  
  nus=rep(0,n)
  cuw=cvw=cuvw=rep(0,n)
  
  if(missing(bandwidth)){
    huw=rep(0.25*n^(-1/(2+dx+dz)),dx+dz)
    hvw=rep(0.25*n^(-1/(2+dy+dz)),dy+dz)
    huvw=rep(0.25*n^(-1/(2+dx+dy+dz)),dx+dy+dz)
    bandwidth=c(huw,hvw,huvw)
  }else{
    huw=bandwidth[1:(dx+dz)]
    hvw=bandwidth[(dx+dz+1):(dx+dz+dy+dz)]
    huvw=bandwidth[(dx+dz+dy+dz+1):(dx+dz+dy+dz+dx+dy+dz)]
  }
  
  if(cpp!=FALSE){
    cuw=kdendveccpp(cbind(U,W),huw)
    cvw=kdendveccpp(cbind(V,W),hvw)
    cuvw=kdendveccpp(cbind(U,V,W),huvw)
  }else{
    cuw=kdendvec(cbind(U,W),huw)
    cvw=kdendvec(cbind(V,W),hvw)
    cuvw=kdendvec(cbind(U,V,W),huvw)
  }

  nus=pmax(1-cuw*cvw/cuvw,0)
  return(mean(nus))
  
  
}