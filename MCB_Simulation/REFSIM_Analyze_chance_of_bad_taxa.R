
cases_arr = c(30:40,66:75,76:79)
diag_mat = matrix(NA,nrow = length(cases_arr),ncol = 7)
for(q in 1:length(cases_arr)){
  filename_results = paste0("C:/MCB2/Results_Previous/REFSIM_results_",cases_arr[q],".RData")
  load(filename_results)
  case_res = NULL
  for(i in 1:length(res)){
    case_res = rbind(case_res,res[[i]]$results[which(res[[i]]$results$methodlabel == 'S = 1.3'),])
  }
  
  
  case_res$nr_bad_reference = as.numeric(case_res$nr_bad_reference)
  case_res$reference_significant = as.numeric(case_res$reference_significant)
  
  mean_nr_bad_ref = mean(case_res$nr_bad_reference)
  sd_nr_bad_ref = sd(case_res$nr_bad_reference)
  
  mean_nr_ref_sig = mean(case_res$reference_significant)
  sd_nr_ref_sig = sd(case_res$reference_significant) / sqrt(nrow(case_res))
  
  mean_nr_ref_sig_cond = NA#mean(case_res$reference_significant[case_res$nr_bad_reference>0])
  sd_nr_ref_sig_cond = NA#sd(case_res$reference_significant[case_res$nr_bad_reference>0]) / sqrt(sum(case_res$nr_bad_reference>0))
  diag_mat[q,2:7] = c(mean_nr_bad_ref,sd_nr_bad_ref,mean_nr_ref_sig,sd_nr_ref_sig,mean_nr_ref_sig_cond,sd_nr_ref_sig_cond)
}

diag_dt = data.frame(Names = rep(NA,length(cases_arr)), BAD_REF = rep(NA,length(cases_arr)),DIAG_SIG = rep(NA,length(cases_arr)))
for(q in 1:length(cases_arr)){
  diag_dt$Names[q] = cases_arr[q]
  diag_dt$BAD_REF[q] = paste0(round(diag_mat[q,2],2),' (',round(diag_mat[q,3],2),')')
  diag_dt$DIAG_SIG[q] = paste0(round(diag_mat[q,4],2),'(',round(diag_mat[q,5],2),')')
}
diag_dt
library(xtable)
print(xtable(diag_dt), include.rownames=FALSE)
