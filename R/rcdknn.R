#-----------------------------
# Ccor function with no scale
# k: k in knn
# typ=0: using the density(i) estimation from the paper(log(c) and then e^log(c))
# typ=1: using the density(i) estimation directly
#-----------------------------
ccorknne=function(x,y,cc=1/16,typ=0,cpp=2){
  x=as.matrix(x)
  y=as.matrix(y)
  n=dim(x)[1]
  u=apply(x,2,rank)/(n+1)
  v=apply(y,2,rank)/(n+1)
  
  k=cc*n^(4/(dim(x)[2]+dim(y)[2]+6))
  k=round(max(1,k))
  
  if(cpp==1){
    D=rcpp_distance(cbind(u,v))
  }else if(cpp==2){
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
  if(typ==0){
    cpl=exp(digamma(k)-digamma(n)-(d)*log(e)-log(Cd)) # copula density by using the paper
  }else if(typ==1){
    cpl=k/(n*Cd*e^d) # copula density by calculating directly
  }
  intg=1-1/cpl #ccor=0.5*\int |1-1/cpl|*cpl*dudv = AVG|1-1/cpl|/2 = AVG(|intg|)/2
  #vpm=mean(abs(intg))/2 # ccor = AVG(|intg|)/2
  vp=mean(ifelse(intg>0,intg,0))#mean(max(intg,0)) # ccor = AVG(intg)_+
  #vm=mean(ifelse(intg<0,-intg,0))#mean(-min(intg,0)) # ccor = AVG(intg)_-
  #np=sum(intg>0) # number of sample points acctually for (intg)_+
  #nm=sum(intg<0) # number of sample points acctually for (intg)_-
  #best=vpm # ccor is equal to AVG(|intg|)/2
  #Pdiff=abs(np-nm)/(np+nm) # The difference of number of sample points in (intg)_+ and (intg)_-
  
  #   if(Pdiff>0.9){ # If differs too much, modify ccor by using the one with larger (sub)sample size.
  #     if(np>nm){
  #       best=vp # ccor = AVG(intg)_+
  #     }else{
  #       best=vm # ccor = AVG(intg)_-
  #     }
  #   }
  # return: 1. modified ccor 2. AVG(|intg|)/2 3. AVG(intg)_+ 4. AVG(intg)_- 
  # 5. # of sample in (intg)_+ 6. # of sample in (intg)_-  
  #return(list("ccore"=best,"v"=vpm,"vp"=vp,"vm"=vm,"np"=np,"nm"=nm))
  return(vp)
}

#==========================================================================
minfc_knn=function(n,cc=cc,typ=typ,cpp=cpp){
  h=0.25*n^(-1/4)
  k=ceiling(2*h*(n+1))
  x=floor(n/k)
  xin=1:n
  yin=NULL
  for(i in 1:k){
    yin=cbind(yin,t(i+k*(0:x)))
  }
  yin=yin[yin<=n]
  r=ccorknne(xin,yin,cc=cc,typ=typ,cpp=cpp)
  #plot(xin,yin)
  return(r)
}

#' Calculate robust copula dependence
#' 
#' This is the function that used to calculate the robust copula dependence (RCD) 
#' between two random variables \code{x} and \code{y}. Note that the length of 
#' \code{x} and \code{y} should be the same. 
#' 
#'  @param x The sampled data, a vector, for the first random variable
#'  @param y The sampled data, a vector, for the second random variable  
#'  @param cc The bandwidth need to be used in the density estimation. This parameter could be missing and a deault value will be applied. 
#'  @param typ The ways to calculate the density
#'  @param cpp Using cpp or not
#'  @return The RCD of \code{x} and \code{y}
#'  @examples
#'  n <- 1000
#'  x <- runif(n)
#'  y <- x + runif(n)
#'  rcdknn(x,y)
#'  @export
rcdknn=function(x,y,cc=1/4,typ=0,cpp=2){  
  x=as.matrix(x)
  y=as.matrix(y)
  n=dim(x)[1]
  k=cc*n^(4/(dim(x)[2]+dim(y)[2]+6))
  k=round(max(1,k))
  maxc=ccorknne(1:n,1:n,cc=cc,typ=typ,cpp=cpp)
  minc=minfc_knn(n,cc=cc,typ=typ,cpp=cpp) 
  score=(ccorknne(x,y,cc=cc,typ=typ,cpp=cpp)-minc)/(maxc-minc)
  return(score)
}
