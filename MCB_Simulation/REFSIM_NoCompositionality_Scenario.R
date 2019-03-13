library(vegan)


REFSIM_NOCOMP_TYPE_SCENARIO_DEF = 'REFSIM_NOCOMP_TYPE_SCENARIO_DEF'

#generate frequencies
REFSIM_generate_NOCOMP_TYPE_Scenario = function(label = "MISSING_LABEL",
                               m_diff_abundant = 50,
                               effect_relative_to_total_sample  = 0.05,
                               n0 = 20,
                               n1 = 20,
                               global_NULL = F,
                               mean_vec = c(rep(200,50),rep(25,150),rep(5,800)),
                               size = 5,
                               chance_for_structural_zero = 0.3,
                               effect_multiplier = 3,
                               select_diff_abundant = c(1)
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


REFSIM_generate_data_for_NOCOMP_Type_Scenario = function(setting_def){
    
  if(class(setting_def) != REFSIM_NOCOMP_TYPE_SCENARIO_DEF){
    stop("Invalid class object used for generation")
  }
  
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
  
  X_unsampled = matrix(NA,nrow = setting_def$n0 + setting_def$n1,ncol = setting_def$m)
  X = matrix(NA,nrow = setting_def$n0 + setting_def$n1,ncol = setting_def$m)
  Y = c(rep(0,setting_def$n0),rep(1,setting_def$n1))  
  
  multiplier_vec = rep(1,m)
  if(!global_NULL)
    multiplier_vec[select_diff_abundant] = effect_multiplier
  for(i in 1:nrow(X)){
    zeroizer_vec = rbinom(m,1,chance_for_structural_zero)
    if(Y[i]==0){
      temp_row = rnbinom(m,size = size,mu = lambda_vec)#rpois(m,lambda_vec)  
    }else{#1
      temp_row = rnbinom(m,size = size,mu = multiplier_vec*lambda_vec) #rpois(m,multiplier_vec*lambda_vec)  
      
    }
    temp_row = (1-zeroizer_vec) * temp_row
    
    X[i,] = temp_row
  }
  
  
  ret = list()
  ret$X = X
  ret$Y = Y
  ret$select_diff_abundant = select_diff_abundant
  ret$n0 = n0
  ret$n1 = n1
  return(ret)
}

