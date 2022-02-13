# This script is used for running two simulations found in the appendices:
# 1) Analyzing the FDR of DACOMP when reference taxa are picked at random (or the most abundant taxa are selected as reference)
# 2) Analyzing the RVP procedure described in appendix B - 
#the Reference Validation Procedure is used for detecting if a diff. abundant taxon has entered the reference set

#%%%%%%%%%%%%%%%%%%%%
#OUTPUT format:
# This scripts outputs text files with:
# -> The T1E estimates of different versions of the RVP, if MODE_COMPUTE_GLOBAL_NULL is TRUE.
#    Output is done to files named REFSIM_RVP_aggregated_results_SCENARIO_ID.txt with SCENARIO_ID being the ID of the scenario.
# -> FDR estimate and SE of FDR estimate, if MODE_COMPUTE_RANDOM_SELECT_FDR is true.
#    Output is done to files named 'REFSIM_FDR_BAD_SELECTION_aggregated_results_SCENARIO_ID.txt' with SCENARIO_ID being the ID of the scenario.
#%%%%%%%%%%%%%%%%%%%%

#Two operation modes: - mutually exclusive !. The next three lines should be commented out and controlled form the batch script
#MODE_COMPUTE_GLOBAL_NULL = T 
#MODE_COMPUTE_RANDOM_SELECT_FDR = F
#SCENARIO_ID = 21 # see batch file for computed values


# 2021/05/31 note that MODE_COMPUTE_GLOBAL_NULL is currently deprecated, and not supported anymore
if(MODE_COMPUTE_GLOBAL_NULL & MODE_COMPUTE_RANDOM_SELECT_FDR){
  stop('MODE_COMPUTE_GLOBAL_NULL and MODE_COMPUTE_RANDOM_SELECT_FDR are mutually exclusive!')
}

NR.WORKERS = 200 # number of workers - for reproducible results - do not change
if(MODE_COMPUTE_GLOBAL_NULL){
  NR.WORKERS = 1000
}
NR.CORES = 96
MAINDIR = './'
Q_LEVEL = 0.1

print(paste0('Start Time:'))
print(Sys.time())

#Modes
CORRECTION_TYPE_DACOMP = 'BH'

MODE_RUN_SIMULATION = T #used to run simulations on multiple cores
if(SCENARIO_ID %in% c(22:25)) #These scenarios do not support RVP, since there is a chance for a single taxon being selected as a refernece set
  MODE_COMPUTE_GLOBAL_NULL = F
MODE_PROCESS_RESULTS = T
#used for debugging:
DEBUG = F #for limited number of reps
DEBUG.SINGLE.CORE = F #for running on single core

GlobalNull.Test.Alpha = 0.1 #T1E for RVP testing
ref_select_method_Median_SD_Threshold = 1.3 #median SD threshold for reference selection RVP
NR.Perms.GN = 200 

#load packages and simulation framework:
library(dacomp)
library(subzero)
library(foreach)
library(doParallel)
library(doRNG)
library(ancom.R)
library(phyloseq)
library(metagenomeSeq)
source(paste0('REFSIM_GenerateSettings_Index.R'))
source('GlobalNull_Over_Reference.R')

#auxilary functions - used for generating names for saving files
REFSIM_results_file = function(RESULTS_DIR,SCENARIO_ID){
  return(paste0(RESULTS_DIR,'/REFSIM_FDR_BAD_SELECTION_results_',SCENARIO_ID,'.RData'))
}

REFSIM_aggregated_results_file = function(RESULTS_DIR,SCENARIO_ID){
  return(paste0(RESULTS_DIR,'/REFSIM_FDR_BAD_SELECTION_aggregated_results_',SCENARIO_ID,'.txt'))
}

REFSIM_RVP_results_file = function(RESULTS_DIR,SCENARIO_ID){
  return(paste0(RESULTS_DIR,'/REFSIM_RVP_results_',SCENARIO_ID,'.RData'))
}

REFSIM_RVP_aggregated_results_file = function(RESULTS_DIR,SCENARIO_ID){
  return(paste0(RESULTS_DIR,'/REFSIM_RVP_aggregated_results_',SCENARIO_ID,'.txt'))
}

RESULTS_DIR = paste0("../../Results/")

RNG_SEED = 1
REPS_PER_SETTING = 1*NR.WORKERS # each worker will run three random generations of the data (in either mdoes)

if(DEBUG){
  NR.WORKERS = 7
  REPS_PER_SETTING = 1*NR.WORKERS
}

NR_REPS_PER_WORKER = ceiling(REPS_PER_SETTING/NR.WORKERS)

# function for scenario based worker
Worker_Function = function(core_nr){
  
  library(dacomp)
  
  TP_random_vec = rep(NA,NR_REPS_PER_WORKER) #TP and FDR for two NAIVE reference selection methods
  FDR_random_vec = rep(NA,NR_REPS_PER_WORKER)
  TP_most_abundant_vec = rep(NA,NR_REPS_PER_WORKER)
  FDR_most_abundant_vec = rep(NA,NR_REPS_PER_WORKER)
  GN_NULL_LIST = list() #results for the RVP procedure
  CONTAINS_DIFF_ABUNDANT = rep(NA,NR_REPS_PER_WORKER)# will indicate if reference tested by RVP has differentially abundant taxa in it.
  #These realizations will be removed when computing T1E of the RVP
  
  for(b in 1:NR_REPS_PER_WORKER){ # for each data realization
    
    if(DEBUG)
      cat(paste0('Case ',b,' out of ',NR_REPS_PER_WORKER,'\n\r'))
    current_setting_generator = REFSIM_SETTINGS_LIST[[SCENARIO_ID]]
    data = REFSIM_generate_setting_wrapper(current_setting_generator) #generate data
    
    #handle naive selection of taxa
    if(MODE_COMPUTE_RANDOM_SELECT_FDR){
      
      Selected_References_random = sample(1:ncol(data$X),size = 50,replace = F) #select at random
      abundance_rank = (ncol(data$X) + 1 - rank(apply(data$X,2,sum),ties.method = 'random'))
      Selected_References_most_abundant = which(abundance_rank<=50) #select the 50 most abundant taxa, by total number of counts in taxa
      
      #test, with reference taxa selected at random
      res = dacomp::dacomp.test(X = data$X,y = data$Y,
                          ind_reference_taxa = Selected_References_random,
                          verbose = T,
                          test = DACOMP.TEST.NAME.WILCOXON,
                          nr_perm = 1/(Q_LEVEL/(ncol(data$X))),
                          q = Q_LEVEL,disable_DSFDR = F)
      
      #compute FDR and TP, for selection at random
      bad_select = length(which(Selected_References_random %in% data$select_diff_abundant))
      
      pvals = res$p.values.test
      pvals[Selected_References_random] = NA 
      rejected = which(p.adjust(pvals,method = CORRECTION_TYPE_DACOMP) < Q_LEVEL) #pvals
      true_positive = length(which(rejected %in% data$select_diff_abundant))
      false_positive = length(rejected) - true_positive
      FDR = 0
      if(length(rejected)>0)
        FDR = false_positive/length(rejected)
      
      TP_random_vec[b] = true_positive
      FDR_random_vec[b] = FDR
      
      # test with the most abundant taxa selected as reference
      res = dacomp::dacomp.test(X = data$X,y = data$Y,
                                ind_reference_taxa = Selected_References_most_abundant,
                                verbose = T,
                                test = DACOMP.TEST.NAME.WILCOXON,
                                nr_perm = 1/(Q_LEVEL/(ncol(data$X))),
                                q = Q_LEVEL,disable_DSFDR = F)
      
      #compute TP and FDR for picking the most abundant taxa as reference:
      bad_select = length(which(Selected_References_most_abundant %in% data$select_diff_abundant))
      
      pvals = res$p.values.test
      pvals[Selected_References_most_abundant] = NA 
      rejected = which(p.adjust(pvals,method = CORRECTION_TYPE_DACOMP) < Q_LEVEL) #pvals
      true_positive = length(which(rejected %in% data$select_diff_abundant))
      false_positive = length(rejected) - true_positive
      FDR = 0
      if(length(rejected)>0)
        FDR = false_positive/length(rejected)
      
      TP_most_abundant_vec[b] = true_positive
      FDR_most_abundant_vec[b] = FDR
    }
    
    #handle running the RVP, if needed by the defined mode:
  } # 
  
  return_object = list(TP_random_vec = TP_random_vec,
                       FDR_random_vec = FDR_random_vec,
                       TP_most_abundant_vec = TP_most_abundant_vec,
                       FDR_most_abundant_vec = FDR_most_abundant_vec,
                       GN_NULL_LIST = GN_NULL_LIST,
                       CONTAINS_DIFF_ABUNDANT = CONTAINS_DIFF_ABUNDANT)
  return(return_object)
} # End of worker function

#handle running of simulation on multiple cores:
if(MODE_RUN_SIMULATION){
  cat(paste0('Running simulation for scenario ID ',SCENARIO_ID,'\n\r'))
  Start = Sys.time()
  rm(res)
  res=NULL
  if(!DEBUG.SINGLE.CORE){
    cat(paste0(' Using multiple cores \n\r'))
    cl = makeCluster(NR.CORES)
    registerDoParallel(cl)
    res = foreach(core = 1:NR.WORKERS, .options.RNG = RNG_SEED) %dorng% {Worker_Function(core)} #.options.RNG=c(1:NR.WORKERS) #%dorng%
    stopCluster(cl)
    cat(paste0(' Run done \n\r'))
  }else{
    cat(paste0(' Using single core \n\r'))
    res = list()
    for(core in 1:NR.WORKERS){
      set.seed(core)
      print(paste0('doing core ',core))
      print(Sys.time())
      res[[core]] = Worker_Function(core)
      print(paste0('done core ',core))
    }
  }
  
  End = Sys.time()
  cat(paste0(' Total time for setting: \n\r'))
  print(End-Start)
  if(MODE_COMPUTE_RANDOM_SELECT_FDR){
    filename = REFSIM_results_file(RESULTS_DIR,SCENARIO_ID)
  }else{# we are in  MODE_COMPUTE_GLOBAL_NULL
    filename = REFSIM_RVP_results_file(RESULTS_DIR,SCENARIO_ID)
  } 
  save(res,file = filename)
}

#process results
if(MODE_PROCESS_RESULTS){
  cat(paste0('Combining results for scenario ',SCENARIO_ID,' \n\r'))
  #the filename changes based on the operation mode
  if(MODE_COMPUTE_RANDOM_SELECT_FDR){
    filename = REFSIM_results_file(RESULTS_DIR,SCENARIO_ID)
  }else{# we are in  MODE_COMPUTE_GLOBAL_NULL
    filename = REFSIM_RVP_results_file(RESULTS_DIR,SCENARIO_ID)
  } 
  load(file = filename) #=> res
  
  #aggregate the results from different cores
  #by mode:
  if(MODE_COMPUTE_RANDOM_SELECT_FDR){ #get FDR of the naive reference selection methods:
    combined_results_FDR_random = res[[1]]$FDR_random_vec
    combined_results_FDR_most_abundant = res[[1]]$FDR_most_abundant_vec
    if(length(res)>1){
      for(i in 2:length(res)){
        combined_results_FDR_random = c(combined_results_FDR_random,res[[i]]$FDR_random_vec)
        combined_results_FDR_most_abundant = c(combined_results_FDR_most_abundant,res[[i]]$FDR_most_abundant_vec)
      }    
    }  
  }
  
  #save to file
  if(MODE_COMPUTE_RANDOM_SELECT_FDR){
    filename = REFSIM_aggregated_results_file(RESULTS_DIR,SCENARIO_ID)
  }else{# we are in  MODE_COMPUTE_GLOBAL_NULL
    filename = REFSIM_RVP_aggregated_results_file(RESULTS_DIR,SCENARIO_ID)
  } 
  
  #Write results to file:
  sink(filename)
  if(MODE_COMPUTE_RANDOM_SELECT_FDR){
    print('FDR - random')
    print(mean(combined_results_FDR_random))
    print('FDR - random - se')
    print(sd(combined_results_FDR_random)/sqrt(length(combined_results_FDR_random)))
    print('FDR - most abundant')
    print(mean(combined_results_FDR_most_abundant))
    print('FDR -  most abundant- se')
    print(sd(combined_results_FDR_most_abundant)/sqrt(length(combined_results_FDR_most_abundant)))
    
  }
  sink()
  
  cat(paste0('Done scenario ',SCENARIO_ID,' \n\r'))
}
