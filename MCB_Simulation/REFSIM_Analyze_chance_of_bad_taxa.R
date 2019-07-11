#This script is used to produce the table found in appendix B.1, detailing the number of diff.abun. taxa that have entered the reference set by mistake
compute_by_method = 'DACOMP,Wilcoxon,rarefaction,S = 1.3' # we show it for one method - since all of the S=1.3 variants, used the same reference set
nr.reps.per.case = 96 #must be consistent with the definition in REFSIM_Single_Scenario.R

#cases to iterate over
cases_arr = c(1:27)
diag_mat = matrix(NA,nrow = length(cases_arr),ncol = 3) # record mean number of diff abundand taxa, with sd

#for each scenario:
for(q in 1:length(cases_arr)){
  #load simulation results:
  filename_results = paste0("../../Results//REFSIM_results_",cases_arr[q],".RData")
  load(filename_results)
  #collect results across cores
  case_res = NULL
  for(i in 1:length(res)){
    case_res = rbind(case_res,res[[i]]$results[which(res[[i]]$results$methodlabel == compute_by_method),])
  }
  
  
  case_res$nr_bad_reference = as.numeric(case_res$nr_bad_reference)
  #compute mean and Sd, and record results
  mean_nr_bad_ref = mean(case_res$nr_bad_reference)
  sd_nr_bad_ref = sd(case_res$nr_bad_reference)
  
  
  diag_mat[q,2:3] = c(mean_nr_bad_ref,sd_nr_bad_ref)
}

#save to file, and print out as latex table:
diag_dt = data.frame(Names = rep(NA,length(cases_arr)), BAD_REF = rep(NA,length(cases_arr)))
for(q in 1:length(cases_arr)){
  diag_dt$Names[q] = cases_arr[q]
  diag_dt$BAD_REF[q] = paste0(round(diag_mat[q,2],2),' (',round(diag_mat[q,3]/sqrt(nr.reps.per.case),2),')')
}
diag_dt
library(xtable)
print(xtable(diag_dt), include.rownames=FALSE)
