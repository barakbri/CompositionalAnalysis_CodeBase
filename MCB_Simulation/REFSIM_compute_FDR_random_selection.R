NR.WORKERS = 8 
SCENARIO_ID = 10 #3,10,16,25 #need to compute for all of these
MAINDIR = './'
Q_LEVEL = 0.1

print(paste0('Start Time:'))
print(Sys.time())

#Modes
CORRECTION_TYPE_SUBZERO = 'BH'
CORRECTION_TYPE_WILCOXON = 'BH'

MODE_RUN_SIMULATION = T
MODE_PROCESS_RESULTS = T
DEBUG = F
DEBUG.SINGLE.CORE = F

DIAG_PLOTS_SELECTION = F

library(subzero)
library(foreach)
library(doParallel)
library(doRNG)
library(ancom.R)
library(phyloseq)
library(metagenomeSeq)
source(paste0('REFSIM_GenerateSettings_Index.R'))

#auxilary functions
REFSIM_results_file = function(RESULTS_DIR,SCENARIO_ID){
  return(paste0(RESULTS_DIR,'/REFSIM_FDR_BAD_SELECTION_results_',SCENARIO_ID,'.RData'))
}

REFSIM_aggregated_results_file = function(RESULTS_DIR,SCENARIO_ID){
  return(paste0(RESULTS_DIR,'/REFSIM_FDR_BAD_SELECTION_aggregated_results_',SCENARIO_ID,'.txt'))
}


RESULTS_DIR = paste0("../../Results/")

RNG_SEED = 1
REPS_PER_SETTING = 9*NR.WORKERS 

if(DEBUG){
  
  NR.WORKERS = 7
  REPS_PER_SETTING = 1*NR.WORKERS
  
}
NR_REPS_PER_WORKER = ceiling(REPS_PER_SETTING/NR.WORKERS)




# function for scenario based worker

Worker_Function = function(core_nr){
  library(subzero)

  TP_random_vec = rep(NA,NR_REPS_PER_WORKER)
  FDR_random_vec = rep(NA,NR_REPS_PER_WORKER)
  TP_most_abundant_vec = rep(NA,NR_REPS_PER_WORKER)
  FDR_most_abundant_vec = rep(NA,NR_REPS_PER_WORKER)
  
  for(b in 1:NR_REPS_PER_WORKER){
    
    if(DEBUG)
      cat(paste0('Case ',b,' out of ',NR_REPS_PER_WORKER,'\n\r'))
    current_setting_generator = REFSIM_SETTINGS_LIST[[SCENARIO_ID]]
    data = REFSIM_generate_setting_wrapper(current_setting_generator)
    
    
    Selected_References_random = sample(1:ncol(data$X),size = 50,replace = F)
    abundance_rank = (ncol(data$X) + 1 - rank(apply(data$X,2,sum),ties.method = 'random'))
    Selected_References_most_abundant = which(abundance_rank<=50)
    
    res = subzero.dfdr.test(X = data$X, y = data$Y,
                            nr_reference_taxa = Selected_References_random, verbose = T, nr_perm = 1/(Q_LEVEL/(ncol(data$X))),
                            q = Q_LEVEL,disable_DS.FDR = F)
    
    bad_select = length(which(Selected_References_random %in% data$select_diff_abundant))
    
    pvals = res$p.values
    pvals[Selected_References_random] = NA 
    rejected = which(p.adjust(pvals,method = CORRECTION_TYPE_SUBZERO) < Q_LEVEL) #pvals
    true_positive = length(which(rejected %in% data$select_diff_abundant))
    false_positive = length(rejected) - true_positive
    FDR = 0
    if(length(rejected)>0)
      FDR = false_positive/length(rejected)
    
    TP_random_vec[b] = true_positive
    FDR_random_vec[b] = FDR
    
    res = subzero.dfdr.test(X = data$X, y = data$Y,
                            nr_reference_taxa = Selected_References_most_abundant, verbose = T, nr_perm = 1/(Q_LEVEL/(ncol(data$X))),
                            q = Q_LEVEL,disable_DS.FDR = F)
    
    bad_select = length(which(Selected_References_most_abundant %in% data$select_diff_abundant))
    
    pvals = res$p.values
    pvals[Selected_References_most_abundant] = NA 
    rejected = which(p.adjust(pvals,method = CORRECTION_TYPE_SUBZERO) < Q_LEVEL) #pvals
    true_positive = length(which(rejected %in% data$select_diff_abundant))
    false_positive = length(rejected) - true_positive
    FDR = 0
    if(length(rejected)>0)
      FDR = false_positive/length(rejected)
    
    TP_most_abundant_vec[b] = true_positive
    FDR_most_abundant_vec[b] = FDR
    
  }  
  return_object = list(TP_random_vec = TP_random_vec,
                       FDR_random_vec = FDR_random_vec,
                       TP_most_abundant_vec = TP_most_abundant_vec,
                       FDR_most_abundant_vec = FDR_most_abundant_vec)
  
  return(return_object)
} # End of worker function



if(MODE_RUN_SIMULATION){
  cat(paste0('Running simulation for scenario ID ',SCENARIO_ID,'\n\r'))
  Start = Sys.time()
  res=NULL
  rm(res)
  res=NULL
  if(!DEBUG.SINGLE.CORE){
    cat(paste0(' Using multiple cores \n\r'))
    cl = makeCluster(NR.WORKERS)
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
  
  filename = REFSIM_results_file(RESULTS_DIR,SCENARIO_ID)
  save(res,file = filename)
}

if(MODE_PROCESS_RESULTS){
  cat(paste0('Combining results for scenario ',SCENARIO_ID,' \n\r'))
  filename = REFSIM_results_file(RESULTS_DIR,SCENARIO_ID)
  load(file = filename) #=> res
  
  combined_results_FDR_random = res[[1]]$FDR_random_vec
  combined_results_FDR_most_abundant = res[[1]]$FDR_most_abundant_vec
  if(length(res)>1){
    for(i in 2:length(res)){
      combined_results_FDR_random = c(combined_results_FDR_random,res[[i]]$FDR_random_vec)
      combined_results_FDR_most_abundant = c(combined_results_FDR_most_abundant,res[[i]]$FDR_most_abundant_vec)
    }    
  }
  
  filename = REFSIM_aggregated_results_file(RESULTS_DIR,SCENARIO_ID)
  sink(filename)
  print('FDR - random')
  print(mean(combined_results_FDR_random))
  print('FDR - random - sd')
  print(sd(combined_results_FDR_random)/sqrt(length(combined_results_FDR_random)))
  print('FDR - most abundant')
  print(mean(combined_results_FDR_most_abundant))
  print('FDR -  most abundant- sd')
  print(sd(combined_results_FDR_most_abundant)/sqrt(length(combined_results_FDR_most_abundant)))
  sink()
  
  cat(paste0('Done scenario ',SCENARIO_ID,' \n\r'))
}
