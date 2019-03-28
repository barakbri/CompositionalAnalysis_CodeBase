BATCH_VEC = c(1:25)
source(paste0('REFSIM_GenerateSettings_Index.R'))
RESULTS_DIR = paste0("../../Results/")
REFSIM_aggregated_results_file = function(RESULTS_DIR,SCENARIO_ID,suffix = ""){ #suffix "sd_" for sd
  return(paste0(RESULTS_DIR,'/REFSIM_aggregated_results_',suffix,SCENARIO_ID,'.RData'))
}


combine_to_file = function(file,suffix){
  filename = REFSIM_aggregated_results_file(RESULTS_DIR,BATCH_VEC[1],suffix = suffix)
  load(file = filename)
  if(suffix == '')
    results_all_settings = aggregated_results
  if(suffix == 'sd_'){
    temp = aggregated_results_sd
    temp$setting_id = 1
    results_all_settings = temp
  }
    
  
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
  
  settings_names = rep(NA,length(REFSIM_SETTINGS_LIST))
  for(i in 1:length(settings_names))
    settings_names[i] = REFSIM_SETTINGS_LIST[[i]]$label
  results_all_settings$setting_name = settings_names[results_all_settings$setting_id]
  
  
  write.csv(results_all_settings,file = paste0(RESULTS_DIR,file))
  
}

combine_to_file('REFSIM_Combined_Results.csv','')
combine_to_file('REFSIM_Combined_Results_sd.csv','sd_')

