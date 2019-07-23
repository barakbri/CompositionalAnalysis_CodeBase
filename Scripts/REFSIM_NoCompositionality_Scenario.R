# Script is used for the scenarios with no compositionality
#- some taxa go up, some go down - for null taxa, the marginal distribution of counts stays the same

library(vegan)

#scenario name
REFSIM_NOCOMP_TYPE_SCENARIO_DEF = 'REFSIM_NOCOMP_TYPE_SCENARIO_DEF'

#generate scenario definition
REFSIM_generate_NOCOMP_TYPE_Scenario = function(label = "MISSING_LABEL",
                               m_diff_abundant = 50, # number of diff abundant taxa
                               effect_relative_to_total_sample  = 0.05, #effect size, relative to the total sample
                               n0 = 20, # sample sizes
                               n1 = 20,
                               global_NULL = F, #is it a global null scenario (with no differentially abundant taxa)
                               mean_vec = c(rep(200,50),rep(25,150),rep(5,800)), #vector of means
                               size = 5, #parameter for negative binomial variance
                               chance_for_structural_zero = 0.0, #change for a structural zero, a coin flip to turn a count to a zero for a sample, at random
                               effect_multiplier = 3, #effect size
                               select_diff_abundant = c(1) #select taxa to be differentially abundant
                               ){
  m = length(mean_vec)
  
  
  
  #USed only in const siZe variant
  index_in_selected_for_reducing = sample(x = 1:length(select_diff_abundant),size = ceiling(length(select_diff_abundant)/2),replace = F)
  
  
  ret = list()
  ret$label = label
  ret$m = m
  ret$m_diff_abundant = m_diff_abundant
  ret$effect_relative_to_total_sample = effect_relative_to_total_sample
  ret$select_diff_abundant = select_diff_abundant
  ret$n0 = n0
  ret$n1 = n1
  ret$global_NULL = global_NULL
  ret$lambda_vec = mean_vec
  ret$chance_for_structural_zero = chance_for_structural_zero
  ret$effect_multiplier = effect_multiplier
  ret$size = size
  class(ret) = REFSIM_NOCOMP_TYPE_SCENARIO_DEF
  return(ret)
}


# Function for generating data, for the no compositionality scenarios, based on a scenario definition
REFSIM_generate_data_for_NOCOMP_Type_Scenario = function(setting_def){
    
  if(class(setting_def) != REFSIM_NOCOMP_TYPE_SCENARIO_DEF){
    stop("Invalid class object used for generation")
  }
  
  #retreive definitions from definition
  label = setting_def$label
  m = setting_def$m
  m_diff_abundant = setting_def$m_diff_abundant
  effect_relative_to_total_sample = setting_def$effect_relative_to_total_sample
  select_diff_abundant = setting_def$select_diff_abundant
  
  n0 = setting_def$n0
  n1 = setting_def$n1
  global_NULL = setting_def$global_NULL
  lambda_vec = setting_def$lambda_vec
  chance_for_structural_zero = setting_def$chance_for_structural_zero
  effect_multiplier = setting_def$effect_multiplier
  size = setting_def$size
  
  #allocate space for returned values
  X_unsampled = matrix(NA,nrow = setting_def$n0 + setting_def$n1,ncol = setting_def$m)
  X = matrix(NA,nrow = setting_def$n0 + setting_def$n1,ncol = setting_def$m)
  Y = c(rep(0,setting_def$n0),rep(1,setting_def$n1))  
  
  # for non global null scenarios, set the vector of multipliers for taxa that go up or down
  multiplier_vec = rep(1,m)
  if(!global_NULL){
    increase_vec = select_diff_abundant[seq(1,to = length(select_diff_abundant),by = 2)]
    decrease_vec = select_diff_abundant[seq(2,to = length(select_diff_abundant),by = 2)]
    multiplier_vec[increase_vec] = 1 + effect_multiplier
    multiplier_vec[decrease_vec] = 1 - effect_multiplier
  }
  
  # generate data, sample by sample
  for(i in 1:nrow(X)){
    zeroizer_vec = rbinom(m,1,chance_for_structural_zero)
    if(Y[i]==0){
      temp_row = rnbinom(m,size = size,mu = lambda_vec) 
    }else{#1
      temp_row = rnbinom(m,size = size,mu = multiplier_vec*lambda_vec) #for the second sample, insert effects
    }
    temp_row = (1-zeroizer_vec) * temp_row
    
    X[i,] = temp_row
  }
  
  #return results
  ret = list()
  ret$X = X
  ret$Y = Y
  ret$select_diff_abundant = select_diff_abundant
  ret$n0 = n0
  ret$n1 = n1
  return(ret)
}

