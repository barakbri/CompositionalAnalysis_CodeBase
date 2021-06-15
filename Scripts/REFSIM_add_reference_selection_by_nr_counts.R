#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Script is used to run the ZINB-Wave method on the simulation settings discussed in the paper.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#install package if needed


#load packages and settigns index
library(dacomp)
source(paste0('REFSIM_GenerateSettings_Index.R')) # scenario library
source('Validation_and_reselection_procedure.R')
#Parameters
PARAM_SCENARIOS = 1:47 # Scenarios to be run
PARAM_REPS = 100 #number of repetitions per scenario
PARAM_Q_LEVEL = 0.1 #FDR level


# Function for evaluating FDR and average number of true positives for a scenario.
evaluate_scenario_Power_and_FDR_for_dacomp_reference_by_counts = function(scenario_id = 1,
                                                                          Q_LEVEL = 0.1,
                                                                          nr_reps = 50){
  
  set.seed(1)
  TP_vec_Wilcoxon = rep(NA,nr_reps)
  FDR_vec_Wilcoxon = rep(NA,nr_reps)
  rejected_vec_Wilcoxon = rep(NA,nr_reps)
  TP_vec_Wilcoxon_Ratio = rep(NA,nr_reps)
  FDR_vec_Wilcoxon_Ratio = rep(NA,nr_reps)
  rejected_vec_Wilcoxon_Ratio = rep(NA,nr_reps)
  TP_vec_Welch = rep(NA,nr_reps)
  FDR_vec_Welch = rep(NA,nr_reps)
  rejected_vec_Welch = rep(NA,nr_reps)
  nr_diff_abundant_in_ref = rep(NA,nr_reps)
  nr_diff_abundant_in_ref_before_clean = rep(NA,nr_reps)
  
  for(r in 1:nr_reps){
    print(paste0('Doing data generation ',r,'/',nr_reps))
    #generate data
    data_generation = REFSIM_generate_setting_wrapper(setting_parameters = REFSIM_SETTINGS_LIST[[scenario_id]])
    
    select_references = dacomp.select_references(X = data_generation$X,
                                                 median_SD_threshold = 0.0001,
                                                 minimal_TA = 100,
                                                 maximal_TA = 500,
                                                 Pseudo_Count_used = 1)
    
    select_references = select_references$selected_references
    
    
    nr_diff_abundant_in_ref[r] = sum(select_references %in% data_generation$select_diff_abundant )
    
    res_dacomp_Wilcoxon = dacomp.test(X = data_generation$X,y = data_generation$Y,ind_reference_taxa = select_references,test = DACOMP.TEST.NAME.WILCOXON,disable_DSFDR = T,compute_ratio_normalization = T,q = PARAM_Q_LEVEL,nr_perm = 1/(Q_LEVEL/(ncol(data_generation$X))))
    res_dacomp_Welch = dacomp.test(X = data_generation$X,y = data_generation$Y,ind_reference_taxa = select_references,test = DACOMP.TEST.NAME.WELCH_LOGSCALE,disable_DSFDR = T,q = PARAM_Q_LEVEL,nr_perm = 1/(Q_LEVEL/(ncol(data_generation$X))))
    
    #check the rejected and collect statistics
    rejected = which(p.adjust(res_dacomp_Wilcoxon$p.values.test,method = 'BH')<=Q_LEVEL)
    rejected_vec_Wilcoxon[r] = length(rejected)
    TP = sum(rejected %in% data_generation$select_diff_abundant)
    FP = length(rejected) - TP
    FDR = FP/max(length(rejected),1)
    TP_vec_Wilcoxon[r] = TP
    FDR_vec_Wilcoxon[r] = FDR
    
    rejected = which(p.adjust(res_dacomp_Wilcoxon$p.values.test.ratio.normalization,method = 'BH')<=Q_LEVEL)
    rejected_vec_Wilcoxon_Ratio[r] = length(rejected)
    TP = sum(rejected %in% data_generation$select_diff_abundant)
    FP = length(rejected) - TP
    FDR = FP/max(length(rejected),1)
    TP_vec_Wilcoxon_Ratio[r] = TP
    FDR_vec_Wilcoxon_Ratio[r] = FDR
    
    rejected = which(p.adjust(res_dacomp_Welch$p.values.test,method = 'BH')<=Q_LEVEL)
    rejected_vec_Welch[r] = length(rejected)
    TP = sum(rejected %in% data_generation$select_diff_abundant)
    FP = length(rejected) - TP
    FDR = FP/max(length(rejected),1)
    TP_vec_Welch[r] = TP
    FDR_vec_Welch[r] = FDR
    
  } # done loop over cases
  
  summary_statistics = c(
          TP_avg_Wilcoxon = mean(TP_vec_Wilcoxon),
          TP_se_Wilcoxon = sd(TP_vec_Wilcoxon)/sqrt(nr_reps),
          FDR_avg_Wilcoxon = mean(FDR_vec_Wilcoxon),
          FDR_se_Wilcoxon = sd(FDR_vec_Wilcoxon)/sqrt(nr_reps),
          rejected_avg_Wilcoxon = mean(rejected_vec_Wilcoxon),
          rejected_se_Wilcoxon = sd(rejected_vec_Wilcoxon)/sqrt(nr_reps),
          
          TP_avg_Wilcoxon_Ratio = mean(TP_vec_Wilcoxon_Ratio),
          TP_se_Wilcoxon_Ratio = sd(TP_vec_Wilcoxon_Ratio)/sqrt(nr_reps),
          FDR_avg_Wilcoxon_Ratio = mean(FDR_vec_Wilcoxon_Ratio),
          FDR_se_Wilcoxon_Ratio = sd(FDR_vec_Wilcoxon_Ratio)/sqrt(nr_reps),
          rejected_avg_Wilcoxon_Ratio = mean(rejected_vec_Wilcoxon_Ratio),
          rejected_se_Wilcoxon_Ratio = sd(rejected_vec_Wilcoxon_Ratio)/sqrt(nr_reps),
          
          TP_avg_Welch = mean(TP_vec_Welch),
          TP_se_Welch = sd(TP_vec_Welch)/sqrt(nr_reps),
          FDR_avg_Welch = mean(FDR_vec_Welch),
          FDR_se_Welch = sd(FDR_vec_Welch)/sqrt(nr_reps),
          rejected_avg_Welch = mean(rejected_vec_Welch),
          rejected_se_Welch = sd(rejected_vec_Welch)/sqrt(nr_reps),
          
          nr_diff_abundant_in_ref_avg = mean(nr_diff_abundant_in_ref),
          nr_diff_abundant_in_ref_se = sd(nr_diff_abundant_in_ref)/sqrt(nr_reps)
          )
  # original_vectors = list(TP_vec = TP_vec,FDR_vec = FDR_vec)
  # ret = list(summary_statistics = summary_statistics,original_vectors = original_vectors)
  return(summary_statistics)
}

#temp = evaluate_scenario_Power_and_FDR_for_dacomp_reference_by_counts(7,0.1,3)
#temp

# Parallelize over cases
  
  library(doParallel)
  library(doRNG)
  cl = NULL
  rm(cl)
  cl <- makeCluster(parallel::detectCores() - 1)
  registerDoParallel(cl)
  
  
  print(paste0('Starting sim'))
  print(Sys.time())
  Start = Sys.time()
  
  #length(REFSIM_SETTINGS_LIST)
  dacomp_parallel_res <- foreach(i=PARAM_SCENARIOS, .options.RNG=1,.combine = 'rbind') %dorng% {
    
    library(dacomp)
    
    ret = list()
    tryCatch(expr = {ret = evaluate_scenario_Power_and_FDR_for_dacomp_reference_by_counts(scenario_id = i,
                                                                        nr_reps = PARAM_REPS,
                                                                        Q_LEVEL = PARAM_Q_LEVEL)})
    return(ret)
  }
  
  stopCluster(cl)
  
  print(paste0('Ending sim'))
  print(Sys.time())
  print('Duration:')
  print(Sys.time() - Start)

#put meaning full row names, by scenarios
rownames(dacomp_parallel_res) = unlist(lapply(REFSIM_SETTINGS_LIST,function(x){x$label}))[PARAM_SCENARIOS] #paste0('Scenario ',PARAM_SCENARIOS)

dacomp_results = as.data.frame(dacomp_parallel_res)
save(dacomp_results,file = '../../Results/dacomp_results_ref_by_counts.rdata')
