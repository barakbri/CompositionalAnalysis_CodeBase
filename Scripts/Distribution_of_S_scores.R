
source(paste0('REFSIM_GenerateSettings_Index.R')) # scenario library
library(dacomp)
library(doRNG)
library(doParallel)
scenario_names = unlist(lapply(REFSIM_SETTINGS_LIST,function(x){x$label}))
scenarios_to_run = c(1:47)

set.seed(1)
data_gen_reps = 100
nr.cores = 47



run_S_sim_for_scenario=function(scenario_id){
  S_for_null = NULL
  S_for_diff_abundant = NULL
  min_S_for_diff_abundant = NULL
  min_lambda_for_diff_abundant = NULL
  for(sim_rep in 1:data_gen_reps){
    print(paste0('Scenario:',scenario_id,', Doing sim rep: ',sim_rep,'/',data_gen_reps))
    
    data_generation = REFSIM_generate_setting_wrapper(setting_parameters = REFSIM_SETTINGS_LIST[[scenario_id]])
    ref_obj = dacomp::dacomp.select_references(X = data_generation$X,
                                               median_SD_threshold = 1.3,
                                               minimal_TA = 10,maximal_TA = 200,
                                               Pseudo_Count_used = 1,verbose = F,
                                               run_in_parallel = F,Nr.Cores = nr.cores)
    if(length(data_generation$select_diff_abundant)>1){
      diff_abun = data_generation$select_diff_abundant
      min_S_for_diff_abundant = c(min(ref_obj$scores[diff_abun]),min_S_for_diff_abundant)
      S_for_null = c(ref_obj$scores[-diff_abun],S_for_null)
      S_for_diff_abundant = c(ref_obj$scores[diff_abun],S_for_diff_abundant)
      taxa_at_contamination_point = which(sort(ref_obj$scores) <= min(ref_obj$scores[diff_abun]))
      min_lambda_for_diff_abundant = c(max(ref_obj$min_abundance_over_the_sorted[taxa_at_contamination_point]),min_lambda_for_diff_abundant)
    }else{
      S_for_null = c(ref_obj$scores,S_for_null)
    }
  }
  return(list(S_for_null = S_for_null,
              S_for_diff_abundant = S_for_diff_abundant,
              min_S_for_diff_abundant = min_S_for_diff_abundant,
              min_lambda_for_diff_abundant = min_lambda_for_diff_abundant))
}

cl <- makeCluster(nr.cores)
registerDoParallel(cl)

parallel_res <- foreach(scenario_id=scenarios_to_run, .options.RNG=1234) %dorng% {
    library(dacomp)
    return(run_S_sim_for_scenario(scenario_id))
}
# stop cluster
stopCluster(cl)

save(parallel_res,file = '../../Results/S_sim_results.rdata')

results_matrix = matrix(NA,nrow = length(scenarios_to_run),ncol = 2)
for(s in 1:length(parallel_res)){
  if(is.null( parallel_res[[s]]$min_S_for_diff_abundant))
    next
  results_matrix[s,1] = paste0(round(mean(parallel_res[[s]]$min_S_for_diff_abundant),2),'(',round(sd(parallel_res[[s]]$min_S_for_diff_abundant)/sqrt(data_gen_reps),2),')')
  results_matrix[s,2] = paste0(round(mean(parallel_res[[s]]$min_lambda_for_diff_abundant),0),'(',round(sd(parallel_res[[s]]$min_lambda_for_diff_abundant)/sqrt(data_gen_reps),0),')')
}

rownames(results_matrix) = scenario_names
colnames(results_matrix) = c('minSDiffAbundant','minCountsDiffAbundant')

library(xtable)
sink(file = '../../Results/S_simulations_results.txt')
xtable(results_matrix)
sink()

extract_scores_table = function(Settings_indices,labels = Settings_indices){
  dt_gg = NULL
  
  for(s in 1:length(Settings_indices)){
    temp = rbind(
      data.frame(Scenario = labels[s],Taxon = 'Diff. Abundant',
                 Scores = parallel_res[[Settings_indices[s]]]$S_for_diff_abundant),
      data.frame(Scenario = labels[s],Taxon = 'Non-Diff. Abundant',
                 Scores = parallel_res[[Settings_indices[s]]]$S_for_null)
    )
    if(is.null(dt_gg)){
      dt_gg = temp
    }else{
      dt_gg = rbind(dt_gg,temp)
    }
  }
  return(dt_gg)
}


library(ggplot2)
labels = c("m1=10, Delta_ml = 0.5",  "m1=100, Delta_ml = 0.5",
           "m1=10, Delta_ml = 1.0",  "m1=100, Delta_ml = 1.0",
           "m1=10, Delta_ml = 1.5",  "m1=100, Delta_ml = 1.5",
           "m1=10, Delta_ml = 2.0",  "m1=100, Delta_ml = 2.0",
           "m1=10, Delta_ml = 2.5",  "m1=100, Delta_ml = 2.5",
           "m1=10, Delta_ml = 3.0",  "m1=100, Delta_ml = 3.0")
dt_gg = extract_scores_table(Settings_indices = c(2:11,26,27),labels = labels)
p = ggplot(dt_gg)+geom_histogram(aes(x = Scores,fill = Taxon),bins = 100,alpha = 0.5,position = "identity")+
  facet_wrap(~Scenario,ncol = 6)+theme_bw()+scale_y_sqrt()+theme(strip.text.x = element_text(size = 8))
ggsave(filename = '../../Results/Sj_distribution_CD.pdf',plot = p,width = 10,height = 4)

labels = c("p_A = 0.9, m1 = 120",
           "p_A = 0.8, m1 = 120",
           "p_A = 0.7, m1 = 120",
           "p_A = 0.6, m1 = 120",
           "p_A = 0.5, m1 = 120",
           "p_A = 0.9, m1 = 60",
           "p_A = 0.8, m1 = 60",
           "p_A = 0.7, m1 = 60",
           "p_A = 0.6, m1 = 60",
           "p_A = 0.5, m1 = 60")
dt_gg = extract_scores_table(Settings_indices = c(12:21),labels = labels)
p = ggplot(dt_gg)+geom_histogram(aes(x = Scores,fill = Taxon),bins = 100,alpha = 0.5,position = "identity")+
  facet_wrap(~Scenario,nrow = 2)+theme_bw()+scale_y_sqrt()+theme(strip.text.x = element_text(size = 8))
ggsave(filename = '../../Results/Sj_distribution_signal_in_rare.pdf',plot = p,width = 10,height = 4)

labels = c("15 samples in each group",
           "20 samples in each group",
           "25 samples in each group",
           "30 samples in each group")
dt_gg = extract_scores_table(Settings_indices = c(22:25),labels = labels)
p = ggplot(dt_gg)+geom_histogram(aes(x = Scores,fill = Taxon),bins = 100,alpha = 0.5,position = "identity")+
  facet_wrap(~Scenario)+theme_bw()+scale_y_sqrt()+theme(strip.text.x = element_text(size = 8))
ggsave(filename = '../../Results/Sj_distribution_no_compositionality.pdf',plot = p,width = 10,height = 4)


labels = c("m1=10, Delta_ml = 1.0,\n group S oversampled",  "m1=10, Delta_ml = 1.0,\n group S undersampled",
           "m1=100, Delta_ml = 1.0,\n group S oversampled",  "m1=100, Delta_ml = 1.0,\n group S undersampled",
           "m1=10, Delta_ml = 2.0,\n group S oversampled",  "m1=10, Delta_ml = 2.0,\n group S undersampled",
           "m1=100, Delta_ml = 2.0,\n group S oversampled",  "m1=100, Delta_ml = 2.0,\n group S undersampled",
           "m1=10, Delta_ml = 3.0,\n group S oversampled",  "m1=10, Delta_ml = 3.0,\n group S undersampled",
           "m1=100, Delta_ml = 3.0,\n group S oversampled",  "m1=100, Delta_ml = 3.0,\n group S undersampled")
dt_gg = extract_scores_table(Settings_indices = c(30:41),labels = labels)
p = ggplot(dt_gg)+geom_histogram(aes(x = Scores,fill = Taxon),bins = 100,alpha = 0.5,position = "identity")+
  facet_wrap(~Scenario,nrow = 2)+theme_bw()+scale_y_sqrt()+theme(strip.text.x = element_text(size = 7))
ggsave(filename = '../../Results/Sj_distribution_confounded_seq_depth.pdf',plot = p,width = 10,height = 4)

labels = c("m1=10,\n 10% increased variance",  "m1=100, 10% increased variance",
           "m1=10,\n 20% increased variance",  "m1=100, 20% increased variance",
           "m1=10,\n 30% increased variance",  "m1=100, 30% increased variance")
           
dt_gg = extract_scores_table(Settings_indices = c(42:47),labels = labels)
p = ggplot(dt_gg)+geom_histogram(aes(x = Scores,fill = Taxon),bins = 100,alpha = 0.5,position = "identity")+
  facet_wrap(~Scenario,nrow = 2)+theme_bw()+scale_y_sqrt()+theme(strip.text.x = element_text(size = 8))
ggsave(filename = '../../Results/Sj_distribution_NB_seq_depth.pdf',plot = p,width = 10,height = 4)