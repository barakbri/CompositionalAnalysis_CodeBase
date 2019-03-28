library(vegan)
load(paste0('qPCR_data.RData')) #=> qPCR_data

REFSIM_QPCR_TYPE_SCENARIO_DEF = 'REFSIM_QPCR_TYPE_SCENARIO_DEF'

select_diff_abundant_10 = sample(1:(ncol(qPCR_data$counts_matrix)),size = 10,replace = F)
select_diff_abundant_scalar_10 = sample(c(0.5,1,1.5),size = 10,replace = T)

select_diff_abundant_100 = sample(1:(ncol(qPCR_data$counts_matrix)),size = 100,replace = F)
select_diff_abundant_scalar_100 = sample(c(0.5,1,1.5),size = 100,replace = T)

#generate frequencies
REFSIM_generate_qPCR_TYPE_Scenario = function(label = "MISSING_LABEL",
                               m_diff_abundant = 10,
                               effect_relative_to_total_sample  = 0.05,
                               n0 = 200,
                               n1 = 200,
                               poisson_mean_reads = 100000,
                               global_NULL = F,
                               Const_Size_Variant = F,
                               Const_Size_Effect_Multiplier = 2,
                               Percent_to_add = 0.5
                               ){
  m = ncol(qPCR_data$counts_matrix)
  select_diff_abundant = 0
  select_diff_abundant_scalar = 1
  if(!global_NULL){
    if(m_diff_abundant == 10){
      select_diff_abundant = select_diff_abundant_10
      select_diff_abundant_scalar = select_diff_abundant_scalar_10
    }else if (m_diff_abundant == 100){
      select_diff_abundant = select_diff_abundant_100
      select_diff_abundant_scalar = select_diff_abundant_scalar_100
    }else{
      select_diff_abundant = select_diff_abundant_100
      select_diff_abundant_scalar = select_diff_abundant_scalar_100
    }
    
  }
  
  
  
  #USed only in const siZe variant
  index_in_selected_for_reducing = sample(x = 1:length(select_diff_abundant),size = ceiling(length(select_diff_abundant)/2),replace = F)
  
  
  ret = list()
  ret$label = label
  ret$m = m
  ret$m_diff_abundant = m_diff_abundant
  ret$effect_relative_to_total_sample = effect_relative_to_total_sample
  ret$select_diff_abundant = select_diff_abundant
  ret$select_diff_abundant_scalar = select_diff_abundant_scalar
  ret$n0 = n0
  ret$n1 = n1
  ret$poisson_mean_reads = poisson_mean_reads
  ret$global_NULL = global_NULL
  ret$Const_Size_Variant = Const_Size_Variant
  ret$index_in_selected_for_reducing = index_in_selected_for_reducing
  ret$Const_Size_Effect_Multiplier = Const_Size_Effect_Multiplier
  ret$Percent_to_add = Percent_to_add
  class(ret) = REFSIM_QPCR_TYPE_SCENARIO_DEF
  return(ret)
}


REFSIM_generate_data_for_qPCR_Type_Scenario = function(setting_def){
    
  if(class(setting_def) != REFSIM_QPCR_TYPE_SCENARIO_DEF){
    stop("Invalid class object used for generation")
  }
  
  OTU_Counts = qPCR_data$counts_matrix
  qPCR_Counts = qPCR_data$qPCR_counts
  
  label = setting_def$label
  m = setting_def$m
  m_diff_abundant = setting_def$m_diff_abundant
  effect_relative_to_total_sample = setting_def$effect_relative_to_total_sample
  select_diff_abundant = setting_def$select_diff_abundant
  taxon_type = setting_def$taxon_type
  n0 = setting_def$n0
  n1 = setting_def$n1
  poisson_mean_reads = setting_def$poisson_mean_reads
  global_NULL = setting_def$global_NULL
  Const_Size_Variant = setting_def$Const_Size_Variant
  Const_Size_Effect_Multiplier = setting_def$Const_Size_Effect_Multiplier
  index_in_selected_for_reducing = setting_def$index_in_selected_for_reducing
  select_diff_abundant_scalar = setting_def$select_diff_abundant_scalar
  Percent_to_add = setting_def$Percent_to_add
  
  X_unsampled = matrix(NA,nrow = setting_def$n0 + setting_def$n1,ncol = setting_def$m)
  X = matrix(NA,nrow = setting_def$n0 + setting_def$n1,ncol = setting_def$m)
  Y = c(rep(0,setting_def$n0),rep(1,setting_def$n1))  
  Total_Original_Counts = rep(NA,length(Y))
  
  for(i in 1:nrow(X)){
    
    sampled_ind = sample(1:nrow(OTU_Counts),size = 1)
    temp_row = as.numeric(OTU_Counts[sampled_ind,])
    
    
    qPCR_counts_in_sample = qPCR_Counts[sampled_ind]
    X_unsampled[i,] = qPCR_counts_in_sample * temp_row/sum(temp_row)
    
    
    if(Y[i] == 1){
      if(!Const_Size_Variant){
        for(s in 1:length(select_diff_abundant)){
          lambda_for_effect = effect_relative_to_total_sample * qPCR_counts_in_sample* select_diff_abundant_scalar[s]*rbinom(1,1,Percent_to_add)/length(select_diff_abundant_scalar)
          X_unsampled[i,select_diff_abundant[s]] = X_unsampled[i,select_diff_abundant[s]] + round(rnorm(1,lambda_for_effect,sqrt(lambda_for_effect)))
        }  
      }else{
        
        index_for_reducing = select_diff_abundant[index_in_selected_for_reducing]
        index_for_adding = select_diff_abundant[-index_in_selected_for_reducing]
        
        total_size_of_reduced = sum(X_unsampled[i,index_for_reducing])
        X_unsampled[i,index_for_reducing] = X_unsampled[i,index_for_reducing]/Const_Size_Effect_Multiplier
        value_to_add = total_size_of_reduced - sum(X_unsampled[i,index_for_reducing])
        X_unsampled[i,index_for_adding] = X_unsampled[i,index_for_adding] + ceiling(value_to_add/length(index_for_adding))
      }
      
    }
    
    Total_Original_Counts[i] = sum(X_unsampled[i,])
    counts_to_sample = rpois(1,poisson_mean_reads)
    
    prob_vector = X_unsampled[i,]/sum(X_unsampled[i,])
    X[i,] = rmultinom(1,counts_to_sample , prob = prob_vector)
    
  }
  
  
  ret = list()
  ret$X = X
  ret$Y = Y
  ret$Total_Original_Counts = Total_Original_Counts
  ret$X_unsampled = X_unsampled
  ret$select_diff_abundant = select_diff_abundant
  ret$taxon_type = taxon_type
  ret$n0 = n0
  ret$n1 = n1
  return(ret)
}

