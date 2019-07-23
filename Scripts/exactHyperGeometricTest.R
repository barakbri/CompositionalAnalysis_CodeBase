# Implementation of a fisher exact test for counts data.
# WARNING: this is not a valid procedure - the function is given here only for the comparison done in the manuscript
# See manuscript for additional details
# X is a matrix of counts, samples being rows.
# Y is a vector of 0\1's - the group lableing
# SelectedReferences is a vector of indices of reference taxa
# the function returns an array of P-values, NA for reference taxa.
exactHyperGeometricTest = function(X,Y,SelectedReferences){
  X_Y0 = apply(X[which(Y==0),],2,sum) # sum of counts, for group Y=0, for each taxon
  X_Y1 = apply(X[which(Y==1),],2,sum) # sum of counts, for group Y=1, for each taxon
  reference_counts_Y0 = sum(X_Y0[SelectedReferences]) #total number of counts in reference taxa by group
  reference_counts_Y1 = sum(X_Y1[SelectedReferences])
  
  #compute P-values , for fisher exact tests:
  p.values = rep(NA,ncol(X))
  for(j in 1:ncol(X)){
    if(j %in% SelectedReferences)
      next
    q0 = X_Y0[j]
    q1 = X_Y1[j]
 

    X_contingency = matrix(NA,nrow = 2,ncol = 2)
    X_contingency[1,1] = q0
    X_contingency[1,2] = reference_counts_Y0
    X_contingency[2,1] = q1
    X_contingency[2,2] = reference_counts_Y1
    p.values[j] = fisher.test(X_contingency)$p.value
    
    
  }
    
  return(p.values)
}