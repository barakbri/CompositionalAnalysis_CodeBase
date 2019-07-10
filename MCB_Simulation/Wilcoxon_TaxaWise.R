#' Function for performing Wilcoxon rank sum tests over a dataset
#' The function performes Wilcoxon rank sum tests over the marginal distributions of taxa, and returns a list with a single item 'p.values' - the array of P-values
#' WARNING: this is not a valid procedure - see manuscript for additional details
#' @param X matrix of counts, rows are samples
#' @param y a vector of 0\1's - group labeling
#' @param normalize should normalization be performed
#' @param normalize.P the ratio of sorted counts used for normalization. 1 is equivlant to TSS normalization.
#' @param nr.perms number of permutations for computing P-values
#'
#' @return
#' @export
#'
#' @examples
wilcoxon_taxa_wise = function(X,y,normalize = F,normalize.P = 1.0,nr.perms = 1/(0.05/ncol(X))){
  pvals = rep(NA,ncol(X))
  #handle normalization, if required
  if(normalize){
    normalizer.vec = rep(1,nrow(X))
    q= floor(normalize.P * ncol(X))
    for(i in 1:nrow(X)){
      sorted = sort(as.numeric(X[i,1:ncol(X)]))
      normalizer.vec[i] = sum(sorted[1:q]) + 1
      X[i,1:ncol(X)] = X[i,1:ncol(X)] /normalizer.vec[i]
    }
  }
  #compute P-values by exact permutations
  for(i in 1:ncol(X)){
    pvals[i] = subzero::PermSumCount.test(X = X[,i],Y = y,B = nr.perms,DoWald = as.integer(0))$p.value
    if(is.nan(pvals[i])){pvals[i] = 1}
  }
  ret = list()
  ret$p.values = pvals
  return(ret)
}
