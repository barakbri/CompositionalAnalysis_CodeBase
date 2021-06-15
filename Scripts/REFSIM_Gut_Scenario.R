#This file is used for building generator (definition objects) for the scenarios based on resampling from real data.


library(vegan)
load(paste0('Gut_Flow_data.RData')) #=> Gut_Flow_data
#Scenario object name
REFSIM_GUT_TYPE_SCENARIO_DEF = 'REFSIM_GUT_TYPE_SCENARIO_DEF'

#select differentially abundant taxa and effect sizes, for simulations.
select_diff_abundant_10 = sample(1:(ncol(Gut_Flow_data$counts_matrix)),size = 10,replace = F)
select_diff_abundant_scalar_10 = sample(c(0.5,1,1.5),size = 10,replace = T)

select_diff_abundant_100 = sample(1:(ncol(Gut_Flow_data$counts_matrix)),size = 100,replace = F)
select_diff_abundant_scalar_100 = sample(c(0.5,1,1.5),size = 100,replace = T)

#generate frequencies
REFSIM_generate_Gut_TYPE_Scenario = function(label = "MISSING_LABEL", #label for scenario
                               m_diff_abundant = 10,# number of diff. abundent taxa
                               effect_relative_to_total_sample  = 0.05, #effect size, over all taxa, in ratio to microbial load
                               n0 = 200, #number of samples, in each group
                               n1 = 200,
                               poisson_mean_reads = 100000, # mean number of reads, for observations
                               global_NULL = F,#is it a global null scenario
                               Const_Size_Variant = F, #for a variant with some taxa going up and some going down, but for null taxa, marginal distribtuions are the same
                               Const_Size_Effect_Multiplier = 2, # Parameters for Const size variant - not used in simulations 
                               Percent_to_add = 0.5,
                               Reads_Multiplier_Group_1 = 1, # for creating a confounder effect between the number of reads and the group labeling
                               NB_dispersion_extra_variance_as_part_mu = NULL
                               ){
  #Pack scenario definitions into a list object
  m = ncol(Gut_Flow_data$counts_matrix)
  select_diff_abundant = 0
  select_diff_abundant_scalar = 1
  if(!global_NULL){
    if(m_diff_abundant == 10){
      select_diff_abundant = select_diff_abundant_10
      select_diff_abundant_scalar = select_diff_abundant_scalar_10
    }else if (m_diff_abundant == 100){
      select_diff_abundant = select_diff_abundant_100
      select_diff_abundant_scalar = select_diff_abundant_scalar_100
    }else{# we don
      select_diff_abundant = select_diff_abundant_100
      select_diff_abundant_scalar = select_diff_abundant_scalar_100
    }
    
  }
  
  
  
  #Used only in const size variant
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
  ret$Reads_Multiplier_Group_1 = Reads_Multiplier_Group_1
  ret$NB_dispersion_extra_variance_as_part_mu = NB_dispersion_extra_variance_as_part_mu
  class(ret) = REFSIM_GUT_TYPE_SCENARIO_DEF
  return(ret)
}

#function for generating data, based on gut scenario definition
REFSIM_generate_data_for_Gut_Type_Scenario = function(setting_def){
    
  if(class(setting_def) != REFSIM_GUT_TYPE_SCENARIO_DEF){
    stop("Invalid class object used for generation")
  }
  #get data from stored object
  OTU_Counts = Gut_Flow_data$counts_matrix
  Flow_counts = Gut_Flow_data$Flow_counts
  # get parameters from scenario definition
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
  Reads_Multiplier_Group_1 = setting_def$Reads_Multiplier_Group_1
  NB_dispersion_extra_variance_as_part_mu = setting_def$NB_dispersion_extra_variance_as_part_mu
  
  # allocate memory for returned types
  X_unsampled = matrix(NA,nrow = setting_def$n0 + setting_def$n1,ncol = setting_def$m)
  X = matrix(NA,nrow = setting_def$n0 + setting_def$n1,ncol = setting_def$m)
  Y = c(rep(0,setting_def$n0),rep(1,setting_def$n1))  
  Total_Original_Counts = rep(NA,length(Y))
  
  #generate, sample by sample
  for(i in 1:nrow(X)){
    #take a base line of counts from the real data, used for modeling overdispersion in counts
    sampled_ind = sample(1:nrow(OTU_Counts),size = 1)
    temp_row = as.numeric(OTU_Counts[sampled_ind,])
    
    
    Flow_counts_in_sample = Flow_counts[sampled_ind]
    #generate unknown total abundance
    X_unsampled[i,] = Flow_counts_in_sample * temp_row/sum(temp_row) 
    
    #if in group 1 , insert signal
    if(Y[i] == 1){
      if(!Const_Size_Variant){
        # for each differentially abundant taxon, increase by an effect size defined in the paper
        for(s in 1:length(select_diff_abundant)){
          lambda_for_effect = effect_relative_to_total_sample * Flow_counts_in_sample* select_diff_abundant_scalar[s]*rbinom(1,1,Percent_to_add)/length(select_diff_abundant_scalar)
          X_unsampled[i,select_diff_abundant[s]] = X_unsampled[i,select_diff_abundant[s]] + round(rnorm(1,lambda_for_effect,sqrt(lambda_for_effect)))
        }  
      }else{
        # Const size variant, some go up, some go down - according to scenario definition
        index_for_reducing = select_diff_abundant[index_in_selected_for_reducing]
        index_for_adding = select_diff_abundant[-index_in_selected_for_reducing]
        
        total_size_of_reduced = sum(X_unsampled[i,index_for_reducing])
        X_unsampled[i,index_for_reducing] = X_unsampled[i,index_for_reducing]/Const_Size_Effect_Multiplier
        value_to_add = total_size_of_reduced - sum(X_unsampled[i,index_for_reducing])
        X_unsampled[i,index_for_adding] = X_unsampled[i,index_for_adding] + ceiling(value_to_add/length(index_for_adding))
      }
      
    }
    
    Total_Original_Counts[i] = sum(X_unsampled[i,])
    
    Confounder_Multiplier = (1*(Y[i]==0)) + Reads_Multiplier_Group_1 * (Y[i]==1)
    if(is.null(NB_dispersion_extra_variance_as_part_mu)){
      counts_to_sample = rpois(1,poisson_mean_reads * Confounder_Multiplier) #pick the observed sampling depth      
    }else{
      mu_to_use = poisson_mean_reads * Confounder_Multiplier
      variance_to_add =   mu_to_use*NB_dispersion_extra_variance_as_part_mu
      size_to_use =    mu_to_use^2/variance_to_add
      counts_to_sample = rnbinom(n = 1,size = size_to_use,mu = mu_to_use) #pick the observed sampling depth            
    }
    
    
    prob_vector = X_unsampled[i,]/sum(X_unsampled[i,])
    X[i,] = rmultinom(1,counts_to_sample , prob = prob_vector) # Do actual sampling
    
  }
  
  #Return results
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

