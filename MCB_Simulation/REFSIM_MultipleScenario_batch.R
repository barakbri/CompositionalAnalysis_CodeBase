#BATCH_START_SCENARIO = 30#53#1#
#BATCH_END_SCENARIO = 40#59#34#

SCENARIO_VEC = c(1:11)#c(1:25)#c(1:25)


for(s in SCENARIO_VEC){
  SCENARIO_ID = s
  source(file = paste0('REFSIM_SingleScenario.R'), echo = F)
}
