

source(paste0('REFSIM_GenerateSettings_Index.R')) # scenario library
library(dacomp)
library(doRNG)
library(doParallel)
scenario_names = unlist(lapply(REFSIM_SETTINGS_LIST,function(x){x$label}))
scenarios_to_run = c(1:11)

set.seed(1)
res_list = list()


Reps = 100
mode = 1
lambda_threshold = 7000
nr.cores = 11
PARAM_NR_PERM = 10000

scenario_worker = function(scenario_id){
  source(paste0('REFSIM_GenerateSettings_Index.R')) # scenario library
  library(dacomp)
  library(doRNG)
  library(doParallel)
  
  ret_obj = NULL
  
  for(b in 1:Reps){
    set.seed(b)

    
    library(dacomp)
    data_generation = REFSIM_generate_setting_wrapper(setting_parameters = REFSIM_SETTINGS_LIST[[scenario_id]])
    ref_obj = dacomp:::parallel_reference_select(X = data_generation$X,
                                                 median_SD_threshold = 0.0001,
                                                 minimal_TA = lambda_threshold,
                                                 maximal_TA = 100000,
                                                 Pseudo_Count_used = 1,verbose = F,run_in_parallel = F)
    
    if(mode == 2 & length(data_generation$select_diff_abundant)>1){
      S_crit_with_contamination = quantile(ref_obj$scores[data_generation$select_diff_abundant],probs = 0.5)
      ref_obj = dacomp:::parallel_reference_select(X = data_generation$X,
                                                   median_SD_threshold = S_crit_with_contamination,
                                                   minimal_TA = 10,
                                                   maximal_TA = 100000,
                                                   Pseudo_Count_used = 1,verbose = F,
                                                   Previous_Reference_Selection_Object = ref_obj,
                                                   run_in_parallel = F)
    }
  
    references_to_use_by_lambda_threshold = ref_obj$selected_references  
    
    
    Nr.Contaminated = sum(data_generation$select_diff_abundant %in% references_to_use_by_lambda_threshold)
    
    cleaned_reference = dacomp.validate_references(X = data_generation$X,
                                                   Y = data_generation$Y,
                                                   ref_obj = ref_obj,
                                                   test = DACOMP.TEST.NAME.WILCOXON,
                                                   Q_validation = 0.1,
                                                   Minimal_Counts_in_ref_threshold = 10,
                                                   Reduction_Factor = 0.5,
                                                   NR_perm = PARAM_NR_PERM,
                                                   Verbose = T,
                                                   disable_DSFDR = T)
    
    Nr.Contaminated.after.cleaning = sum(data_generation$select_diff_abundant %in% cleaned_reference)
    if(length(references_to_use_by_lambda_threshold)< ncol(data_generation$X)){
      testing_res = dacomp::dacomp.test(X = data_generation$X,
                                        y = data_generation$Y,
                                        ind_reference_taxa = references_to_use_by_lambda_threshold,
                                        test = DACOMP.TEST.NAME.WILCOXON,q = 0.1,disable_DSFDR = T,
                                        nr_perm  = PARAM_NR_PERM,verbose  = T)
      rejected = which(p.adjust(testing_res$p.values.test,method = 'BH') <= 0.1)
      R = length(rejected)
      V = 0
      tp = 0
      if(length(data_generation$select_diff_abundant)>0)
        tp = sum(rejected %in% data_generation$select_diff_abundant)
      if(R>0){
        V = R - tp
        R.FDR = R
      }else{
        R.FDR = 1
      }
      
    }else{
      R=NA
      V=NA
      tp=NA
      R.FDR = NA
    }
    
    if(length(cleaned_reference)<=ncol(data_generation$X)){
      testing_res_cleaned = dacomp::dacomp.test(X = data_generation$X,
                                                y = data_generation$Y,
                                                ind_reference_taxa = cleaned_reference,
                                                test = DACOMP.TEST.NAME.WILCOXON,q = 0.1,
                                                nr_perm  = PARAM_NR_PERM, verbose  = T,disable_DSFDR = T)
      rejected = which(p.adjust(testing_res_cleaned$p.values.test,method = 'BH') <= 0.1)
      R_cleaned = length(rejected)
      V_cleaned = 0
      tp_cleaned = 0
      if(length(data_generation$select_diff_abundant)>0)
        tp_cleaned = sum(rejected %in% data_generation$select_diff_abundant)
      if(R_cleaned>0){
        V_cleaned = R_cleaned - tp_cleaned
        R.FDR_cleaned = R_cleaned
      }else{
        R.FDR_cleaned = 1
      }
      
    }else{
      R_cleaned = NA
      V_cleaned = NA
      tp_cleaned = NA
      R.FDR_cleaned = NA
    }
    temp = (c(Nr.Contaminated = Nr.Contaminated,
              Nr.Contaminated.after.cleaning = Nr.Contaminated.after.cleaning,
              V = V,
              R=R,
              TP = R-V,
              FDR = V/R.FDR, 
              V_cleaned = V_cleaned,
              R_cleaned = R_cleaned,
              TP_cleaned = R_cleaned - V_cleaned,
              FDR_cleaned = V_cleaned/R.FDR_cleaned))
    
    if(b==1){
      ret_obj = temp
    }else{
      ret_obj = rbind(ret_obj,temp)
    }
  }
  return(ret_obj)
}

#temp = scenario_worker(3)

cl <- makeCluster(nr.cores)
registerDoParallel(cl)


res_list <- foreach(scenario_id=scenarios_to_run, .options.RNG=1234) %dorng% { #scenarios_to_run
  return(scenario_worker(scenario_id))
}


stopCluster(cl)


res_list_2 = res_list[which(unlist(lapply(res_list,function(x){!is.null(x)})))]

save(res_list,file = paste0('../../Results/contamination_sim_mode_',mode,'_reads_thres_',lambda_threshold,'_res_list.Rdata'))
save(res_list_2,file = paste0('../../Results/contamination_sim_mode_',mode,'_reads_thres_',lambda_threshold,'_res_list_2.Rdata'))


combined_res = matrix(data = unlist(lapply(res_list_2,FUN = function(x){apply(x,2,mean,na.rm=T)})),ncol = 10,byrow = T)
combined_res_sd = matrix(data = unlist(lapply(res_list_2,FUN = function(x){apply(x,2,sd,na.rm=T)})),ncol = 10,byrow = T)
colnames(combined_res) = colnames(res_list_2[[1]])
row.names(combined_res) = scenario_names[scenarios_to_run]
colnames(combined_res_sd) = colnames(res_list_2[[1]])
row.names(combined_res_sd) = scenario_names[scenarios_to_run]

save(combined_res,file = paste0('../../Results/contamination_sim_mode_',mode,'_reads_thres_',lambda_threshold,'.Rdata'))
save(combined_res_sd,file = paste0('../../Results/contamination_sim_mode_',mode,'_reads_thres_',lambda_threshold,'_sd.Rdata'))
#View(combined_res)
#View(combined_res_sd)
