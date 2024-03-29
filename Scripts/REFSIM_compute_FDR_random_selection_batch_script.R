
# This is the batch script for running `REFSIM_compute_FDR_random_selection.R` with various parameters (scenarios and activation mode)
# Simulation output is text file, with the results used for constructing tables in appendix B.
# See `REFSIM_compute_FDR_random_selection.R` for additional details.

#Run simulations for estimating FDR for DACOMP when using naive reference selection methods. discussed in appendix B.
SCENARIO_VEC = c(3,10,16,25)

MODE_COMPUTE_GLOBAL_NULL = F
MODE_COMPUTE_RANDOM_SELECT_FDR = T

for(s in SCENARIO_VEC){
  SCENARIO_ID = s
  source(file = paste0('REFSIM_compute_FDR_random_selection.R'), echo = F)
}
