# This script is used to compile the aggregated results (across reps of data generations)
# for the different scenarios to to file, giving the mean (power, FDR) and sd (power, FDR) of the different scenarios

BATCH_VEC = c(1:47) #scenarios to run over
source(paste0('REFSIM_GenerateSettings_Index.R')) #scenario library - scenario names will be used from here
RESULTS_DIR = paste0("../../Results/") # where results will be outputed to 

#auxiliarry function for the data file with aggregeted reusults. added an option for entering a suffix (for SD file)
REFSIM_aggregated_results_file = function(RESULTS_DIR,SCENARIO_ID,suffix = ""){ #suffix "sd_" for sd
  return(paste0(RESULTS_DIR,'/REFSIM_aggregated_results_',suffix,SCENARIO_ID,'.RData'))
}

#main function for combining results, called twice at the end of file
combine_to_file = function(file,suffix){
  #read first file
  filename = REFSIM_aggregated_results_file(RESULTS_DIR,BATCH_VEC[1],suffix = suffix)
  load(file = filename)
  if(suffix == '')
    results_all_settings = aggregated_results
  if(suffix == 'sd_'){
    temp = aggregated_results_sd
    temp$setting_id = 1
    results_all_settings = temp
  }
    
  #aggregate results over all files
  for(s in (2):length(BATCH_VEC)){
    
    filename = REFSIM_aggregated_results_file(RESULTS_DIR,BATCH_VEC[s],suffix = suffix)
    load(file = filename)
    if(suffix == '')
      results_all_settings = rbind(results_all_settings, aggregated_results)
    if(suffix == 'sd_'){
      temp = aggregated_results_sd
      temp$setting_id = s
      results_all_settings = rbind(results_all_settings, temp,suffix = suffix)
    }
      
    
  }
  
  #add setting names
  settings_names = rep(NA,length(REFSIM_SETTINGS_LIST))
  for(i in 1:length(settings_names))
    settings_names[i] = REFSIM_SETTINGS_LIST[[i]]$label
  results_all_settings$setting_name = settings_names[results_all_settings$setting_id]
  
  #save to file
  write.csv(results_all_settings,file = paste0(RESULTS_DIR,file))
  
}

#combine point estiamte and SDs across scenarios
combine_to_file('REFSIM_Combined_Results.csv','')
combine_to_file('REFSIM_Combined_Results_sd.csv','sd_')

