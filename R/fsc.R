fsc=function(X,y){
  n=dim(X)[1]; p=dim(X)[2]
  vm=data.frame(matrix(0,p,p)) # save conditional dependence in ith feature, and jth avaiable feature
  selected_index=rep(NA,p) # save order of the selected index (some permutation of 1:p)
  # 1. loop: no condition, and save p candidates value to the first row of vm
  for(i in 1:p){
    vm[1,i]=rcd(cbind(y,X[,i]))
  }
  selected_index[1]=which.max(vm[1,]) # first selected feature(index) set to first value
  selected=X[,selected_index[1]] # "selected": selected features used to condition on
  # 2. loop: from the 2 to p-1, all are condition version. The last(p) is what's left, no need to compute.
  for(i in 2:(p-1)){ # selection of the ith related feature
    # print(i)
    selected_max_value=0 # max conditional value of the rest p-i+1 candidate features
    restfeature=sort(setdiff(1:p,selected_index))
    restfeature=restfeature[!is.na(restfeature)]
    for(j in restfeature){# loop all the rest (unselected) features(index ordered)
      # print(paste("j=",j))
      vm[i,j]=crcd(y,X[,j],selected)# calculate the candidate feature j (for the ith selection)
      if(selected_max_value<vm[i,j]){
        selected_max_value=vm[i,j] # update max value
        selected_index[i]=j # update max index 
        # Since max_value is always initialized to 0, there must be one j is the max.
      }
    }
    selected=cbind(selected,X[,selected_index[i]]) # update the selected features(data) used for condition
  }
  selected_index[p]=setdiff(1:p,selected_index)
  fit=list(id=selected_index,value=vm,X=X)
  class(fit)="rcdfs"
  return(fit)
}
