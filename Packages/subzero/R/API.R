

#' Permutation based wilcoxon rank sum test
#' Function for performing the Wilcoxon rank sum test, with permutations
#' @param X a vector of counts
#' @param Y a vector of 0,1's - group labeling
#' @param B number of permutations
#' @param DoWald Should an early stopping rule be used
#' @param return_perms should the permutation P-values be given
#' @param disable.ties.correction should the correction for ties be disabled 
#'
#' @return a list with the test statistic, P.value (2-sided), number of permutations and the asymptotic P.value. Item perms will be included if permutations were requested
#' @export
#'
#' @examples
PermSumCount.test = function(X,Y,B = as.integer(c(40000)),DoWald = as.integer(c(1)),return_perms = F,disable.ties.correction = F){
  ranked_X = rank(X,ties.method = 'average') # rank, ties get average score - this is known to give optimal powe in monotone scenarios
  perms_val = as.integer(c(0))
  if(return_perms)
    perms_val = as.integer(c(1))
  res = subzero::rcpp_Wilcoxon_PermTest(ranked_X,Y,B,DoWald,perms_val) # Run C level code - this is extremely fast...
  #compute correction for ties - 
  NR_Y0 = sum(Y == 0)
  N = length(Y)
  NR_Y1 = N - NR_Y0
  
  E_H0 = sum(ranked_X) * NR_Y1 / N
  
  unique_value_counts = table(ranked_X)
  which_unique_value_counts_are_ties = which(unique_value_counts>1)
  unique_value_counts_only_ties_correction = 0
  if(length(which_unique_value_counts_are_ties)>0 & !disable.ties.correction){
    t_r = unique_value_counts[which_unique_value_counts_are_ties]
    unique_value_counts_only_ties_correction  = sum(t_r * (t_r^2 - 1))
  }
  
  V_H0 = NR_Y0 * NR_Y1 * (N+1) / 12 - NR_Y0 * NR_Y1 * unique_value_counts_only_ties_correction / (12 * N * (N - 1))
  Z_score = (res[[1]] - E_H0)/sqrt(V_H0)
  pval.by.normal = 2*(1-pnorm(abs(Z_score)))
  
  ret = list(stat = res[[1]],p.value = res[[2]],b = res[[3]],pval.by.normal=pval.by.normal)
  if(return_perms)
    ret$perms = res[[7]]
  return(ret)
}


#' Function for performing the Wilcoxon rank test, over multiple rarefactions of the data.
#' Pvalue is computed by averaging over the P.values of multiple rarefactions.
#' WARNING: This function was used to test the this approach is invalid - do not use in practice.
#' See paper for further details
#' @param List_of_X list of vectors, corresponding to multiple rarefactions
#' @param Y vector of 0/1's - group labeling
#' @param B number of permutations
#' @param is.Asymptotic logical - should use asymptotic P-values
#' @param disable.ties.correction logical - should the correction for ties be disabled 
#'
#' @return list with a single item - the average P.value across different rarefactions
#' @export
#'
#' @examples
PermSumCount.multipleX.test = function(List_of_X,Y,B = as.integer(c(1000)),is.Asymptotic = F,disable.ties.correction = F){
  p = length(List_of_X)
  pval.vec = rep(NA,p)
  if(is.Asymptotic)
    B = 0
  
  for(i in 1:p){
    res = subzero::PermSumCount.test(List_of_X[[i]],Y,B,DoWald = as.integer(c(0)),return_perms = F,disable.ties.correction = disable.ties.correction)
    if(is.Asymptotic)
      pval.vec[i] = res$pval.by.normal
    else
      pval.vec[i] = res$p.value
  }
  
  p.value = mean(pval.vec)
  return(list(p.value=p.value))
}
