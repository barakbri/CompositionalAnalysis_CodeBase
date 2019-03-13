wilcoxon_taxa_wise = function(X,y,normalize = F,normalize.P = 1.0,paired = F){
  pvals = rep(NA,ncol(X))
  if(normalize){
    normalizer.vec = rep(1,nrow(X))
    q= floor(normalize.P * ncol(X))
    for(i in 1:nrow(X)){
      sorted = sort(as.numeric(X[i,1:ncol(X)]))
      normalizer.vec[i] = sum(sorted[1:q]) + 1
      X[i,1:ncol(X)] = X[i,1:ncol(X)] /normalizer.vec[i]
    }
  }
  for(i in 1:ncol(X)){
      pvals[i] = wilcox.test(X[y==0,i],X[y==1,i],exact = F,paired = paired)$p.value
    if(is.nan(pvals[i])){pvals[i] = 1}
  }
  ret = list()
  ret$p.values = pvals
  return(ret)
}

#temp = Scenario_simulation_advanced_null_10()  
#temp2 = wilcoxon_taxa_wise(temp$X, as.numeric(as.factor(temp$Y))-1)
#temp3 = wilcoxon_taxa_wise(temp$X, as.numeric(as.factor(temp$Y))-1,normalize = T)
#temp2$p.values
#temp3$p.values

#X = temp$X
#y = as.numeric(as.factor(temp$Y))-1
