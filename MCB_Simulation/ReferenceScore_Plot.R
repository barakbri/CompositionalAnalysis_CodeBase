
plot_ref_select_scores = function(ref_select,label,quantiles_of_scores_to_plot = c(0.5,0.7,0.9), breaks_param = 30){
  Target_MinAbundance_values = ref_select$target_abundance
  hist(ref_select$scores,breaks = breaks_param,main = label,xlab = "Score metric for reference selection")
  sorted_scores = sort(ref_select$scores)
  threshold_ind = 1
  for(threshold_ind in 1:length(Target_MinAbundance_values)){
    abline(v = sorted_scores[ref_select$arg_obj[threshold_ind]],col = threshold_ind+1,lwd = 3)
    text(x = sorted_scores[ref_select$arg_obj[threshold_ind]]+0.1,y=3+1*threshold_ind,labels = Target_MinAbundance_values[threshold_ind],col = threshold_ind+1)
  }  
  qunatiles_in_scores = quantile(ref_select$scores,probs = quantiles_of_scores_to_plot)
  for(i in 1:length(quantiles_of_scores_to_plot)){
    abline(v = qunatiles_in_scores,col = 1,lty = 2,lwd = 3)
    text(x = qunatiles_in_scores[i] + 0.1, y=0+i,labels = quantiles_of_scores_to_plot[i],col = 1)
  }
  
}




