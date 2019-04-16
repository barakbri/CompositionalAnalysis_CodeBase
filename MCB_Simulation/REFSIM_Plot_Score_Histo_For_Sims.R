
#Modes
CORRECTION_TYPE_SUBZERO = 'BH'
CORRECTION_TYPE_WILCOXON = 'BH'

library(subzero)
source('ReferenceScore_Plot.R')
source(paste0('REFSIM_GenerateSettings_Index.R'))
RESULTS_DIR = paste0("../../Results/")

ref_selection_results = data.frame(label = NA,nr_taxa = NA,nr_selected = NA,actual_abundance = NA)
selected_obj_list = list()
set.seed(1)
SCENARIO_ID = 30

DO_PLOTS = T
DO_PLOTS_TO_FILE = T



for(SCENARIO_ID in 1:25){
  print(paste0('Doing scenario ',SCENARIO_ID))
  if(DO_PLOTS_TO_FILE)
    pdf(file =  paste0(RESULTS_DIR,'/RefSelection_Scores_Analysis_SCENARIO_',SCENARIO_ID,'.pdf'),width = 8,height = 5)
  current_setting_generator = REFSIM_SETTINGS_LIST[[SCENARIO_ID]]  
  data = REFSIM_generate_setting_wrapper(current_setting_generator)  
  
  ref_select = select.references.Median.SD.Threshold(X = data$X,median_SD_threshold = 1.3,minimal_TA = 10,maximal_TA = 200)
  
  Selected_References = ref_select$selected_references
  ActualMinAbundance = NULL
  if(length(Selected_References)>1){
    ActualMinAbundance = min(apply(data$X[,Selected_References],1,sum)  )
  }else{
    ActualMinAbundance = min(data$X[,Selected_References])
  }
  selected_obj_list[[SCENARIO_ID]] = ref_select
  ref_selection_results[SCENARIO_ID,] = c(current_setting_generator$label,ncol(data$X),length(Selected_References),ActualMinAbundance)
  
  if(DO_PLOTS){
    plot_ref_select_scores(ref_select = ref_select,
                           label = "")#paste0('Case: ',SCENARIO_ID," : ",current_setting_generator$label)
    # hist(ref_select$scores,breaks = 30,main = paste0('Case: ',SCENARIO_ID," : ",current_setting_generator$label))
    # sorted_scores = sort(ref_select$scores)
    # threshold_ind = 1
    # for(threshold_ind in 1:length(Target_MinAbundance_values)){
    #   abline(v = sorted_scores[ref_select$arg_obj[threshold_ind]],col = threshold_ind+1,lwd = 3)
    #   text(x = sorted_scores[ref_select$arg_obj[threshold_ind]]+0.1,y=10,labels = Target_MinAbundance_values[threshold_ind],col = threshold_ind+1)
    # }  
  }
  if(DO_PLOTS_TO_FILE)
    dev.off()
  
}


write.csv(ref_selection_results,file = paste0(RESULTS_DIR,'/Ref_Selction_Compare_results.csv'))


set.seed(1)
SCENARIO_REPS = 10
m1_scores_list = list()
SCENARIO_TO_COMPUTE_vec = c(1:25)
vec_of_scores = NULL
ref_size_list = list()
vec_of_ref_size = NULL
for(SCENARIO_ID in SCENARIO_TO_COMPUTE_vec){
  print(paste0('Doing scenario ',SCENARIO_ID))
  current_setting_generator = REFSIM_SETTINGS_LIST[[SCENARIO_ID]]  
  vec_of_scores = NULL
  vec_of_ref_size = NULL
  for(r in 1:SCENARIO_REPS){
    print(paste0('r: ',r))
    data = REFSIM_generate_setting_wrapper(current_setting_generator)  
    ref_select = select.references.Median.SD.Threshold(X = data$X,median_SD_threshold = 1.3,minimal_TA = 10,maximal_TA = 200)
    vec_of_scores = c(vec_of_scores,ref_select$scores[data$select_diff_abundant])
    vec_of_ref_size = c(vec_of_ref_size,length(ref_select$selected_references))
  }
  m1_scores_list[[SCENARIO_ID]] = vec_of_scores
  ref_size_list[[SCENARIO_ID]] = vec_of_ref_size
  print(paste0('Mean ref size: ',mean(vec_of_ref_size),' SD: ',sd(vec_of_ref_size)))
}

MIN_REFSCORE_H1 = SCENARIO_TO_COMPUTE_vec
for(i in 1:length(MIN_REFSCORE_H1)){
  MIN_REFSCORE_H1[i] = min(m1_scores_list[[SCENARIO_TO_COMPUTE_vec[i]]])
}

dt_refscore_m1_sim = data.frame(SCENARIO_TO_COMPUTE_vec,MIN_REFSCORE_H1)
dt_refscore_m1_sim
write.csv(dt_refscore_m1_sim,paste0(RESULTS_DIR,'/REFSCORE_H1.csv'),row.names = F)
