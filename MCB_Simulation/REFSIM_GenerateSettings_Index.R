
source(paste0('REFSIM_NO_OVERDISPERSION_SCENARIOS.R'))
source(paste0('REFSIM_qPCR_SCENARIO.R'))
source(paste0('REFSIM_NoCompositionality_Scenario.R'))

REFSIM_generate_setting_wrapper = function(setting_parameters){
  
  if(class(setting_parameters) == REFSIM_QPCR_TYPE_SCENARIO_DEF){
    return(REFSIM_generate_data_for_qPCR_Type_Scenario(setting_parameters))
    
  }else if(class(setting_parameters) == REFSIM_NO_OVERDISPERSION_SCENARIO_DEF){
    return(REFSIM_generate_data_for_NO_OVERDISPERSION_Type_Scenario(setting_parameters))
    
  }else if(class(setting_parameters) == REFSIM_NOCOMP_TYPE_SCENARIO_DEF){
    return(REFSIM_generate_data_for_NOCOMP_Type_Scenario(setting_parameters))
    
  }else{
    stop('Error: Unidentified setting class')
  }
}

set.seed(1)
EFFECT_SCALAR = 0.1
n_qpcr = 20
qPCR_MEDIAN_READS = 22449

REFSIM_SETTINGS_LIST = list()


# 1 - 11 # data resample

REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_qPCR_TYPE_Scenario(label = 'qPCR, GlobalNull, 20:20', m_diff_abundant = 0,global_NULL = T, effect_relative_to_total_sample = 0.00, n0 = n_qpcr,n1 = n_qpcr, poisson_mean_reads = qPCR_MEDIAN_READS)
REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_qPCR_TYPE_Scenario(label = 'qPCR, m1 = 10, 20:20, effect = 0.5', m_diff_abundant = 10,global_NULL = F, effect_relative_to_total_sample = 0.5, n0 = n_qpcr,n1 = n_qpcr, poisson_mean_reads = qPCR_MEDIAN_READS) 
REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_qPCR_TYPE_Scenario(label = 'qPCR, m1 = 100, 20:20, effect = 0.5', m_diff_abundant = 100,global_NULL = F, effect_relative_to_total_sample = 0.5, n0 = n_qpcr,n1 = n_qpcr, poisson_mean_reads = qPCR_MEDIAN_READS) 
REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_qPCR_TYPE_Scenario(label = 'qPCR, m1 = 10, 20:20, effect = 1.0', m_diff_abundant = 10,global_NULL = F, effect_relative_to_total_sample = 1.0, n0 = n_qpcr,n1 = n_qpcr, poisson_mean_reads = qPCR_MEDIAN_READS) 
REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_qPCR_TYPE_Scenario(label = 'qPCR, m1 = 100, 20:20, effect = 1.0', m_diff_abundant = 100,global_NULL = F, effect_relative_to_total_sample = 1.0, n0 = n_qpcr,n1 = n_qpcr, poisson_mean_reads = qPCR_MEDIAN_READS) 
REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_qPCR_TYPE_Scenario(label = 'qPCR, m1 = 10, 20:20, effect = 1.5', m_diff_abundant = 10,global_NULL = F, effect_relative_to_total_sample = 1.5, n0 = n_qpcr,n1 = n_qpcr, poisson_mean_reads = qPCR_MEDIAN_READS) 
REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_qPCR_TYPE_Scenario(label = 'qPCR, m1 = 100, 20:20, effect = 1.5', m_diff_abundant = 100,global_NULL = F, effect_relative_to_total_sample = 1.5, n0 = n_qpcr,n1 = n_qpcr, poisson_mean_reads = qPCR_MEDIAN_READS)
REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_qPCR_TYPE_Scenario(label = 'qPCR, m1 = 10, 20:20, effect = 2.0', m_diff_abundant = 10,global_NULL = F, effect_relative_to_total_sample = 2.0, n0 = n_qpcr,n1 = n_qpcr, poisson_mean_reads = qPCR_MEDIAN_READS) 
REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_qPCR_TYPE_Scenario(label = 'qPCR, m1 = 100, 20:20, effect = 2.0', m_diff_abundant = 100,global_NULL = F, effect_relative_to_total_sample = 2.0, n0 = n_qpcr,n1 = n_qpcr, poisson_mean_reads = qPCR_MEDIAN_READS) 
REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_qPCR_TYPE_Scenario(label = 'qPCR, m1 = 10, 20:20, effect = 2.5', m_diff_abundant = 10,global_NULL = F, effect_relative_to_total_sample = 2.5, n0 = n_qpcr,n1 = n_qpcr, poisson_mean_reads = qPCR_MEDIAN_READS) 
REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_qPCR_TYPE_Scenario(label = 'qPCR, m1 = 100, 20:20, effect = 2.5', m_diff_abundant = 100,global_NULL = F, effect_relative_to_total_sample = 2.5, n0 = n_qpcr,n1 = n_qpcr, poisson_mean_reads = qPCR_MEDIAN_READS) 

# 12 - 21 #NO OD
REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_NO_OVERDISPERSION_TYPE_Scenario(label = "Signal in low, m_high = 30, p_high = 0.9, n_1=n_2=20,m1 = 120",m_high = 30,p_high = 0.9,n_1 = 20,n_2 = 20)
REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_NO_OVERDISPERSION_TYPE_Scenario(label = "Signal in low, m_high = 30, p_high = 0.8, n_1=n_2=20,m1 = 120",m_high = 30,p_high = 0.8,n_1 = 20,n_2 = 20)
REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_NO_OVERDISPERSION_TYPE_Scenario(label = "Signal in low, m_high = 30, p_high = 0.7, n_1=n_2=20,m1 = 120",m_high = 30,p_high = 0.7,n_1 = 20,n_2 = 20)
REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_NO_OVERDISPERSION_TYPE_Scenario(label = "Signal in low, m_high = 30, p_high = 0.6, n_1=n_2=20,m1 = 120",m_high = 30,p_high = 0.6,n_1 = 20,n_2 = 20)
REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_NO_OVERDISPERSION_TYPE_Scenario(label = "Signal in low, m_high = 30, p_high = 0.5, n_1=n_2=20,m1 = 120",m_high = 30,p_high = 0.5,n_1 = 20,n_2 = 20)
REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_NO_OVERDISPERSION_TYPE_Scenario(label = "Signal in low, m_high = 30, p_high = 0.9, n_1=n_2=20,m1 = 60",m_high = 30,m1 = 60,p_high = 0.9,n_1 = 20,n_2 = 20)
REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_NO_OVERDISPERSION_TYPE_Scenario(label = "Signal in low, m_high = 30, p_high = 0.8, n_1=n_2=20,m1 = 60",m_high = 30,m1 = 60,p_high = 0.8,n_1 = 20,n_2 = 20)
REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_NO_OVERDISPERSION_TYPE_Scenario(label = "Signal in low, m_high = 30, p_high = 0.7, n_1=n_2=20,m1 = 60",m_high = 30,m1 = 60,p_high = 0.7,n_1 = 20,n_2 = 20)
REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_NO_OVERDISPERSION_TYPE_Scenario(label = "Signal in low, m_high = 30, p_high = 0.6, n_1=n_2=20,m1 = 60",m_high = 30,m1 = 60,p_high = 0.6,n_1 = 20,n_2 = 20)
REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_NO_OVERDISPERSION_TYPE_Scenario(label = "Signal in low, m_high = 30, p_high = 0.5, n_1=n_2=20,m1 = 60",m_high = 30,m1 = 60,p_high = 0.5,n_1 = 20,n_2 = 20)

# 22 -25, no compositionality
select_diff_abundant = c(1:10,51:60,201:230)

REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_NOCOMP_TYPE_Scenario(label = "Scenario with no compositionality, 15:15 ",m_diff_abundant=100, mean_vec = 1*c(rep(200,50),rep(20,150),rep(1,800)),effect_multiplier = 5,n0 = 15,n1=15,select_diff_abundant = select_diff_abundant)
REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_NOCOMP_TYPE_Scenario(label = "Scenario with no compositionality, 20:20 ",m_diff_abundant=100, mean_vec = 1*c(rep(200,50),rep(20,150),rep(1,800)),effect_multiplier = 5,n0 = 20,n1=20,select_diff_abundant = select_diff_abundant)
REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_NOCOMP_TYPE_Scenario(label = "Scenario with no compositionality, 25:25 ",m_diff_abundant=100, mean_vec = 1*c(rep(200,50),rep(20,150),rep(1,800)),effect_multiplier = 5,n0 = 25,n1=25,select_diff_abundant = select_diff_abundant)
REFSIM_SETTINGS_LIST[[length(REFSIM_SETTINGS_LIST) + 1]] = REFSIM_generate_NOCOMP_TYPE_Scenario(label = "Scenario with no compositionality, 30:30 ",m_diff_abundant=100, mean_vec = 1*c(rep(200,50),rep(20,150),rep(1,800)),effect_multiplier = 5,n0 = 30,n1=30,select_diff_abundant = select_diff_abundant)