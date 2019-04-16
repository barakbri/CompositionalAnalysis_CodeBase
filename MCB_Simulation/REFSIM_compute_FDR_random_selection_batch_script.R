
SCENARIO_VEC = c(2,3,10,11,12,16,17,21) #c(2,3,10,11,12,16,17,21,25)


for(s in SCENARIO_VEC){
  SCENARIO_ID = s
  source(file = paste0('REFSIM_compute_FDR_random_selection.R'), echo = F)
}
