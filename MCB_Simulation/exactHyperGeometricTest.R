exactHyperGeometricTest = function(X,Y,SelectedReferences){
  X_Y0 = apply(X[which(Y==0),],2,sum)
  X_Y1 = apply(X[which(Y==1),],2,sum)
  reference_counts_Y0 = sum(X_Y0[SelectedReferences])
  reference_counts_Y1 = sum(X_Y1[SelectedReferences])
  p.values = rep(NA,ncol(X))
  for(j in 1:ncol(X)){
    if(j %in% SelectedReferences)
      next
    q0 = X_Y0[j]
    q1 = X_Y1[j]
    #k = q0 + q1
    #m = q0 + reference_counts_Y0
    #n = q1 + reference_counts_Y1

    X_contingency = matrix(NA,nrow = 2,ncol = 2)
    X_contingency[1,1] = q0
    X_contingency[1,2] = reference_counts_Y0
    X_contingency[2,1] = q1
    X_contingency[2,2] = reference_counts_Y1
    p.values[j] = fisher.test(X_contingency)$p.value

    #p.values[j] =prop.test(c(q0,q1),c(q0+reference_counts_Y0,q1+reference_counts_Y1))$p.value
    
    
  }
    
  return(p.values)
}

if(F){
  source('REFSIM_ANCOM_SCENARIOS_2.R')
  set.seed(1)
  #setting_def = REFSIM_generate_ANCOM_TYPE_2_Scenario(m_high = 30,p_high = 0.9,n_1 = 20,n_2 = 20)
  #setting_def = REFSIM_generate_ANCOM_TYPE_2_Scenario(m_high = 30,p_high = 0.8,n_1 = 20,n_2 = 20)
  setting_def = REFSIM_generate_ANCOM_TYPE_2_Scenario(m_high = 30,p_high = 0.7,n_1 = 20,n_2 = 20)
  
  data = REFSIM_generate_data_for_ANCOM_Type_2_Scenario(setting_def)
  
  set.seed(1)
  MAINDIR = 'E:/MCB2/'
  source('E:/MCB2/MCB2/MCB_Simulation/REFSIM_qPCR_SCENARIO.R', echo=TRUE)
  scenario = REFSIM_generate_qPCR_TYPE_Scenario(label = 'TestScenario', 
                                                m_diff_abundant = 20,effect_relative_to_total_sample = 2.0,
                                                n0 = 30,n1 = 30,
                                                poisson_mean_reads = 22000)
  
  data = REFSIM_generate_data_for_qPCR_Type_Scenario(scenario)
  
  
  image(t(log(data$X[,data$select_diff_abundant]+1))) 
  image(t(log(data$X+1))) 
  plot(log(data$X[,data$select_diff_abundant[3]]+1))
  image(log10(t(data$X)+1))
  dim(data$X)
  
  
  source('SelectReferences_MedianSD_Threshold.R')
  library(subzero)
  ref_obj = select.references.Median.SD.Threshold(data$X,median_SD_threshold = 1.2)  
  ref_obj$selected_references
  min(apply(data$X[,ref_obj$selected_references,drop = F],1,sum))
  which(ref_obj$selected_references %in% data$select_diff_abundant)
  
  pvals = exactHyperGeometricTest(data$X,data$Y,ref_obj$selected_references)
  hist(pvals)
  rejected = which(p.adjust(pvals,method = 'BH')<=0.1)
  data$select_diff_abundant
  
  tp = sum(data$select_diff_abundant %in% rejected)
  fp = length(rejected) - tp
  tp  
  fp
}