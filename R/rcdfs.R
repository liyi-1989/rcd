plot.rcdfs=function(fit){
  p=dim(fit$X)[2]; id=fit$id
  barplot(order(id),1:p,width = rep(1,p),space=rep(0,p),ylim=c(0,p+0.5),col="yellow",border = F,
          main="Feature Selection Results",xlab="Feature (column) Number",ylab="Rank")
  lines(1:p-0.5,order(id),col="blue",type="o",pch=16)
}

graph=function(fit,type="dot"){
  M=t(apply(as.matrix(fit$value), MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))[-dim(fit$X)[2],]
  X=melt(M)
  p=nrow(M)
  Max=apply(fit$value,1,which.max)[1:p]
  
  X[,"na"]=X[,2]==paste0("X",rep(Max,time=p+1))
  names(X)=c("Var1","Var2","value","na")
  p1 = ggplot(data=X, aes(x=Var2, y=Var1, fill=value)) +
    geom_raster() +
    scale_fill_gradient2(low="blue", high="yellow", na.value="black", name="") +
    geom_point(colour = "blue",aes(size=ifelse(na, "dot", "no_dot"))) +
    scale_size_manual(values=c(dot=6, no_dot=NA), guide="none") +
    labs(title="Feature Selection Results",x="Feature (column) Number",y="Iteration")+
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 
  #============================================================
  frames = X[X$na, c("Var1", "Var2")]
  frames[,2]=as.numeric(gsub("[^0-9]", "", frames[,2]))
  
  p2 = ggplot(data=X) +
    geom_raster(aes(x=Var2, y=Var1, fill=value)) +
    scale_fill_gradient2(low="blue", high="yellow", na.value="black", name="") +
    geom_rect(data=frames, size=1, fill=NA, colour="blue",
              aes(xmin=Var2 - 0.5, xmax=Var2 + 0.5, ymin=Var1 - 0.5, ymax=Var1 + 0.5)) +
    labs(title="Feature Selection Results",x="Feature (column) Number",y="Iteration")+
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 
  
  if(type=="dot"){
    return(p1)
  }else if(type=="rec"){
    return(p2)
  }else{
    stop("Type must be either dot or rec!")
  }
  
}

# library(gplots)
# heatmap.2(t(scale(t(as.matrix(fit$value))))[-dim(fit$X)[2],],
#           dendrogram = "none",key = T, Rowv=FALSE, Colv=FALSE,trace='none',
#           col=colorRampPalette(c("black", "red"))(n = 1000))
# 
# heatmap.2(t(apply(as.matrix(fit$value), MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))[-dim(fit$X)[2],],
#           dendrogram = "none",key = T, Rowv=FALSE, Colv=FALSE,trace='none',col=colorRampPalette(c("yellow", "orange"))(n = 1000),
#           xlab="Features",ylab="Iteration")
# heatmap(t(apply(as.matrix(fit$value), MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))[-dim(fit$X)[2],],
#         Rowv=NA, Colv=NA,
#         col=colorRampPalette(c("yellow", "orange"))(n = 1000),
#         xlab="Features",ylab="Iteration")
# 
# 
# points(2,3,pch=19,cex=1,col="blue")





summary.rcdfs=function(fit){
  cat("******* Summary of the Feature Selection Result *******\n")
  cat("* Number of observations: ", dim(fit$X)[1],"\n")
  cat("* Number of features:     ", dim(fit$X)[2],"\n")
  cat("* Ranked features:        ", fit$id,"\n")
}


