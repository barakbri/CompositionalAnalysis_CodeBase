#This script is used for generating plots of the reference selection scores to file.
# The plots shown in appendix B are generated using this script

# The script partially loads the simulation framework given in REFSIM_single_scenario.R

library(subzero)
source('ReferenceScore_Plot.R')
source(paste0('REFSIM_GenerateSettings_Index.R'))
RESULTS_DIR = paste0("../../Results/")

# Will be used to store the number of taxa in each scenario, reference size and its abundance (minium across subjects)
ref_selection_results = data.frame(label = NA,nr_taxa = NA,nr_selected = NA,actual_abundance = NA)
selected_obj_list = list()
set.seed(1)
SCENARIO_ID = 30 # Will be updated based on the current scenario

#Should do plots, and should output be to file
DO_PLOTS = T
DO_PLOTS_TO_FILE = T

#Iterate over scenarios
for(SCENARIO_ID in 1:25){
  print(paste0('Doing scenario ',SCENARIO_ID))
  if(DO_PLOTS_TO_FILE)
    pdf(file =  paste0(RESULTS_DIR,'/RefSelection_Scores_Analysis_SCENARIO_',SCENARIO_ID,'.pdf'),width = 8,height = 5)
  #Generate data
  current_setting_generator = REFSIM_SETTINGS_LIST[[SCENARIO_ID]]  
  data = REFSIM_generate_setting_wrapper(current_setting_generator)  
  
  #select references
  ref_select = dacomp::dacomp.select_references(X = data$X,median_SD_threshold = 1.3,minimal_TA = 10,maximal_TA = 200)
  
  Selected_References = ref_select$selected_references
  ActualMinAbundance = NULL
  
  #compute minimal abundance in reference taxa across samples
  if(length(Selected_References)>1){
    ActualMinAbundance = min(apply(data$X[,Selected_References],1,sum)  )
  }else{
    ActualMinAbundance = min(data$X[,Selected_References])
  }
  #record results
  selected_obj_list[[SCENARIO_ID]] = ref_select
  ref_selection_results[SCENARIO_ID,] = c(current_setting_generator$label,ncol(data$X),length(Selected_References),ActualMinAbundance)
  
  if(DO_PLOTS){
    plot_ref_select_scores(ref_select = ref_select,
                           label = "")
  }
  if(DO_PLOTS_TO_FILE)
    dev.off()
  
}


write.csv(ref_selection_results,file = paste0(RESULTS_DIR,'/Ref_Selction_Compare_results.csv'))

# This section of the is used to estimate the distribution of the references scores for the differentially 
# abundant taxa across the different scenarios and size of reference sets.
# This was used to better understand how to set the medianSD threshold.
set.seed(1)
SCENARIO_REPS = 10
m1_scores_list = list() # distribution of differentially abundant scores across scenarios
SCENARIO_TO_COMPUTE_vec = c(1:25)
vec_of_scores = NULL
ref_size_list = list() # sizes of reference sets across scenarios
vec_of_ref_size = NULL
for(SCENARIO_ID in SCENARIO_TO_COMPUTE_vec){ #iterate over scenarios
  print(paste0('Doing scenario ',SCENARIO_ID))
  current_setting_generator = REFSIM_SETTINGS_LIST[[SCENARIO_ID]]  
  vec_of_scores = NULL
  vec_of_ref_size = NULL
  for(r in 1:SCENARIO_REPS){ #several repetitions for each scenario
    print(paste0('r: ',r))
    #generate data and select references:
    data = REFSIM_generate_setting_wrapper(current_setting_generator)
    ref_select = select.references.Median.SD.Threshold(X = data$X,median_SD_threshold = 1.3,minimal_TA = 10,maximal_TA = 200)
    vec_of_scores = c(vec_of_scores,ref_select$scores[data$select_diff_abundant])
    vec_of_ref_size = c(vec_of_ref_size,length(ref_select$selected_references))
  }
  m1_scores_list[[SCENARIO_ID]] = vec_of_scores
  ref_size_list[[SCENARIO_ID]] = vec_of_ref_size
  print(paste0('Mean ref size: ',mean(vec_of_ref_size),' SD: ',sd(vec_of_ref_size)))
}

#record results
MIN_REFSCORE_H1 = SCENARIO_TO_COMPUTE_vec
for(i in 1:length(MIN_REFSCORE_H1)){
  MIN_REFSCORE_H1[i] = min(m1_scores_list[[SCENARIO_TO_COMPUTE_vec[i]]])
}

dt_refscore_m1_sim = data.frame(SCENARIO_TO_COMPUTE_vec,MIN_REFSCORE_H1)
dt_refscore_m1_sim
write.csv(dt_refscore_m1_sim,paste0(RESULTS_DIR,'/REFSCORE_H1.csv'),row.names = F)
