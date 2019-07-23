# Script is used for generating data from a multinomial setting - no overdispersion
# For this scenario, some of the taxa are of high abundance - but they are not Diff.Abundant.
# Only the low abundance taxa may be dif.abundent

#SCENARIO NAME
REFSIM_NO_OVERDISPERSION_SCENARIO_DEF = 'REFSIM_NO_OVERDISPERSION_SCENARIO_DEF'

#function for generating scenario definition
REFSIM_generate_NO_OVERDISPERSION_TYPE_Scenario = function(label = "MISSING_LABEL",
                                               m = 300, #nr taxa
                                               m_high= 17, #nr taxa of high abundance
                                               m1 = 120, #nr taxa that are diff.abun.
                                               signal_proportion = 0.35,# proportion of load that belongs to the disease - see manuscript for the generating model
                                               p_high = 0.95, # proportion of microbial load held by a the high abundance taxa
                                               N.reads = 2500, # number of reads for microbial sampling
                                               n_1 = 50,# sample sizes
                                               n_2 = 50
){
  
  ret = list()
  ret$label = label
  ret$m = m
  ret$m_high = m_high
  ret$m1 = m1
  ret$signal_proportion = signal_proportion
  ret$p_high = p_high
  ret$N.reads = N.reads
  ret$n_1 = n_1
  ret$n_2 = n_2
  class(ret) = REFSIM_NO_OVERDISPERSION_SCENARIO_DEF
  return(ret)
}

# function for generating data based on scenario definitions
REFSIM_generate_data_for_NO_OVERDISPERSION_Type_Scenario = function(setting_def){
  if(class(setting_def) != REFSIM_NO_OVERDISPERSION_SCENARIO_DEF){
    stop("Invalid class object used for generation")
  }
  #extract parameters from definitions object
  m = setting_def$m
  m_high = setting_def$m_high
  m1 = setting_def$m1
  signal_proportion = setting_def$signal_proportion
  p_high = setting_def$p_high
  N.reads = setting_def$N.reads
  n_1 = setting_def$n_1
  n_2 = setting_def$n_2
  
  #generate probability vectors
  prob1 = c(rep(p_high/m_high,m_high),rep((1-p_high)/(m-m_high),m-m_high))
  prob1 = prob1/sum(prob1)
  
  #generate signal vector
  signal_vec = c(rep(0,m_high),rep(1,m1),rep(0,m-m_high-m1))
  signal_vec = signal_vec / sum(signal_vec)
  prob2 = (1-signal_proportion)* prob1 + signal_proportion*signal_vec
  # perform multinomial samping
  X_Y0 = t(rmultinom(n = n_1,size = N.reads,prob = prob1))
  X_Y1 = t(rmultinom(n = n_2,size = N.reads,prob = prob2))
  
  #return results
  X = rbind(X_Y0,X_Y1)
  Y = c(rep(0,n_1),rep(1,n_2))
  
  
  
  ret = list()
  ret$X = X
  ret$Y = Y
  ret$select_diff_abundant = which(signal_vec != 0)
  ret$prob1 = prob1
  ret$prob2 = prob2
  return(ret)
}

