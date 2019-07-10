SCENARIO_VEC = c(26:29)#c(1:25)


for(s in SCENARIO_VEC){
  SCENARIO_ID = s
  source(file = paste0('REFSIM_SingleScenario.R'), echo = F)
}
