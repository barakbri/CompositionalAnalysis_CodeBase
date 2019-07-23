#This script is used for running REFSIM_SingleScenario.R across multiple scenarios
# See readme on how this file fits in the overall pipeline for getting results for the paper

SCENARIO_VEC = c(1:27)

for(s in SCENARIO_VEC){
  SCENARIO_ID = s
  source(file = paste0('REFSIM_SingleScenario.R'), echo = F)
}
