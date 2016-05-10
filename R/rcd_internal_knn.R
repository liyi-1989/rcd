rcd.internal.knn=function(x,k,cpp="parallel",S=F){
  # 1. if scaled, use scaled version
  if(S){return(rcd.internal.knn.scale(x,k=k,cpp=cpp))}
  n=nrow(x)
  d=ncol(x)
  X=apply(x,2,rank)/(n+1)
  if(d==1){stop("X is only one dimensional!")}
  if(missing(k)){
    k=0.25*n^(4/(d+6))
    k=round(max(1,k))
  }
  
  if(cpp=="serial"){
    D=rcpp_distance(X)
  }else if(cpp=="parallel"){
    D=rcpp_parallel_distance(X)
  }else{
    #df=data.frame(X) #data points in each row
    D=as.matrix(dist(data.frame(X),method="maximum"))# distance matrix with maximum norm
  }
  
  # e=rep(0,n)#e: epsilon(i) = 2*distance to the k nearest neighbour(knn)
  # for(i in 1:n){
  #   e[i]=2*sort(D[i,])[k+1] # sort the distance matrix(col&row eighter is OK) and take the knn
  # }
  
  e=2*apply(D,1,sort)[,k+1]
  
  #d=dim(x)[2] # d:dimension of (u,v)/Cd: vol of ball in dim d
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


rcd.internal.knn.scale=function(x,k,cpp="parallel"){
  n=nrow(x);dx=ncol(x)
  bw=0.25*n^(-1/(2+dx))
  kk=ceiling(2*bw*(n+1));x=floor(n/kk);xin=1:n;yin=NULL
  for(i in 1:kk){
    yin=cbind(yin,t(i+kk*(0:x)))
  }
  
  maxc=rcd.internal.knn(matrix(rep(xin,dx),n,dx),k,cpp=cpp,S=F)
  minc=rcd.internal.knn(cbind(xin,yin),k,cpp=cpp,S=F)
  score=rcd.internal.knn(x,k,cpp=cpp,S=F)
  r=(score-minc)/(maxc-minc)
  return(ifelse(r>0,r,score))
}










