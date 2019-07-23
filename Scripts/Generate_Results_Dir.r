#This script is used for creating the results directory

RESULTS_DIR = paste0("../../Results/")
if(!dir.exists(RESULTS_DIR)){
  dir.create(RESULTS_DIR)
  print('Results directory created successfully!')
}