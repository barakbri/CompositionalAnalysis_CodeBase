BATCH_VEC = c(30:40,66:75,76:79)



MAINDIR = '~/'
RESULTS_DIR = paste0(MAINDIR,"/Results/")
REFSIM_aggregated_results_file = function(RESULTS_DIR,SCENARIO_ID){
  return(paste0(RESULTS_DIR,'/REFSIM_aggregated_results_',SCENARIO_ID,'.RData'))
}

ROW_ORDER = c(2:3)#c(1:4)

filename = REFSIM_aggregated_results_file(RESULTS_DIR,BATCH_VEC[1])
load(file = filename)
results_all_settings = aggregated_results#[ROW_ORDER,]

for(s in (2):length(BATCH_VEC)){
  
  filename = REFSIM_aggregated_results_file(RESULTS_DIR,BATCH_VEC[s])
  load(file = filename)
  results_all_settings = rbind(results_all_settings, aggregated_results)#[ROW_ORDER,])
}

settings_names = rep(NA,length(REFSIM_SETTINGS_LIST))
for(i in 1:length(settings_names))
  settings_names[i] = REFSIM_SETTINGS_LIST[[i]]$label
results_all_settings$setting_name = settings_names[results_all_settings$setting_id]


write.csv(results_all_settings,file = paste0(RESULTS_DIR,'REFSIM_Combined_Results.csv'))
