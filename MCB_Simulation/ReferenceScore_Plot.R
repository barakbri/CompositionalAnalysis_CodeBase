

#' The function is used to produce a plot of the median SD scores.
#' Note that this function has been replaced by a better function in the DACOMP function, given here only for reproducibility
#' @param ref_select refernce selection object
#' @param label label for chart
#' @param quantiles_of_scores_to_plot pecentiles for vertical lines
#' @param breaks_param number of bins
#'
#' @return
#' @export
#'
#' @examples
plot_ref_select_scores = function(ref_select,label,quantiles_of_scores_to_plot = c(0.5,0.7,0.9), breaks_param = 30){
  Target_MinAbundance_values = ref_select$target_abundance
  hist(ref_select$scores,breaks = breaks_param,main = label,xlab = "Score metric for reference selection")
  sorted_scores = sort(ref_select$scores)
  threshold_ind = 1
  for(threshold_ind in 1:length(Target_MinAbundance_values)){
    abline(v = sorted_scores[ref_select$arg_obj[threshold_ind]],col = threshold_ind+1,lwd = 3)
  }  
  qunatiles_in_scores = quantile(ref_select$scores,probs = quantiles_of_scores_to_plot)
  for(i in 1:length(quantiles_of_scores_to_plot)){
    abline(v = qunatiles_in_scores,col = 1,lty = 2,lwd = 3)
  }
  
}




