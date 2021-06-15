#This script is used to produce the table found in SM S4.4 detailing the number of diff.abun. taxa that have entered the reference set by mistake
load('../../Results/dacomp_results_ref_by_counts.rdata')
names(dacomp_results)

dacomp_results$nr_diff_abundant_in_ref_avg_and_se =  paste0(
  round(dacomp_results[,c(19)],2),
  '(', round(dacomp_results[,c(20)],2),')'
)

library(xtable)
print(xtable(
  dacomp_results[,'nr_diff_abundant_in_ref_avg_and_se',drop=F]
  ), include.rownames=T)
max(dacomp_results[,c(20)])


# Results for contamination study in 4.2
load('../../Results/contamination_sim_mode_1_reads_thres_3000.Rdata')
combined_res = as.data.frame(combined_res)
print(xtable(combined_res[,c(1,6,5,2,10,9)]), include.rownames=T)
load('../../Results/contamination_sim_mode_1_reads_thres_3000_sd.Rdata')
combined_res_sd = as.data.frame(combined_res_sd)
(as.numeric(round(apply(combined_res_sd,2,max)/sqrt(100),2)))[c(1,6,5,2,10,9)]

load('../../Results/contamination_sim_mode_1_reads_thres_7000.Rdata')
combined_res = as.data.frame(combined_res)
print(xtable(combined_res[,c(1,6,5,2,10,9)]), include.rownames=T)
load('../../Results/contamination_sim_mode_1_reads_thres_7000_sd.Rdata')
combined_res_sd = as.data.frame(combined_res_sd)
(as.numeric(round(apply(combined_res_sd,2,max)/sqrt(100),2)))[c(1,6,5,2,10,9)]

