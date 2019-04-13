NR.WORKERS = 7
SCENARIO_ID = 21 #3,10,16,25 #need to compute for all of these
MAINDIR = './'
Q_LEVEL = 0.1

print(paste0('Start Time:'))
print(Sys.time())

#Modes
CORRECTION_TYPE_SUBZERO = 'BH'
CORRECTION_TYPE_WILCOXON = 'BH'

MODE_RUN_SIMULATION = T
MODE_COMPUTE_GLOBAL_NULL = T
MODE_COMPUTE_RANDOM_SELECT_FDR = F
MODE_PROCESS_RESULTS = T
DEBUG = F
DEBUG.SINGLE.CORE = F

GlobalNull.Test.Alpha = 0.1
ref_select_method_Median_SD_Threshold = 1.3
NR.Perms.GN = 200

library(subzero)
library(foreach)
library(doParallel)
library(doRNG)
library(ancom.R)
library(phyloseq)
library(metagenomeSeq)
source(paste0('REFSIM_GenerateSettings_Index.R'))
source('GlobalNull_Over_Reference.R')

#auxilary functions
REFSIM_results_file = function(RESULTS_DIR,SCENARIO_ID){
  return(paste0(RESULTS_DIR,'/REFSIM_FDR_BAD_SELECTION_results_',SCENARIO_ID,'.RData'))
}

REFSIM_aggregated_results_file = function(RESULTS_DIR,SCENARIO_ID){
  return(paste0(RESULTS_DIR,'/REFSIM_FDR_BAD_SELECTION_aggregated_results_',SCENARIO_ID,'.txt'))
}


RESULTS_DIR = paste0("../../Results/")

RNG_SEED = 1
REPS_PER_SETTING = 30*NR.WORKERS 

if(DEBUG){
  
  NR.WORKERS = 7
  REPS_PER_SETTING = 30*NR.WORKERS
  
}
NR_REPS_PER_WORKER = ceiling(REPS_PER_SETTING/NR.WORKERS)




# function for scenario based worker

Worker_Function = function(core_nr){
  library(subzero)

  TP_random_vec = rep(NA,NR_REPS_PER_WORKER)
  FDR_random_vec = rep(NA,NR_REPS_PER_WORKER)
  TP_most_abundant_vec = rep(NA,NR_REPS_PER_WORKER)
  FDR_most_abundant_vec = rep(NA,NR_REPS_PER_WORKER)
  GN_NULL_LIST = list()
  CONTAINS_DIFF_ABUNDANT = rep(NA,NR_REPS_PER_WORKER)
  for(b in 1:NR_REPS_PER_WORKER){
    
    if(DEBUG)
      cat(paste0('Case ',b,' out of ',NR_REPS_PER_WORKER,'\n\r'))
    current_setting_generator = REFSIM_SETTINGS_LIST[[SCENARIO_ID]]
    data = REFSIM_generate_setting_wrapper(current_setting_generator)
    
    if(MODE_COMPUTE_RANDOM_SELECT_FDR){
      
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
    
    if(MODE_COMPUTE_GLOBAL_NULL){
      ref_select = select.references.Median.SD.Threshold(X = data$X,
                                                         median_SD_threshold = ref_select_method_Median_SD_Threshold,
                                                         minimal_TA = 10,maximal_TA = 200,Psuedo_Count_used = 1,factor_by_Median_Score = F)
      Selected_References = ref_select$selected_references
      X_ref = data$X[,Selected_References,drop = F]
      print(dim(X_ref))
      # if(nrow(X_ref) ==1){
      #   aaa =1
      #   aaa = aaa +1
      # }
      CONTAINS_DIFF_ABUNDANT[b] = length(which(Selected_References %in% data$select_diff_abundant))
      gn_res = compute_GlobalNull_Reference(X_ref,Y = data$Y,nr.perm = NR.Perms.GN)
      #Rejections.counter.global.null = Rejections.counter.global.null + (unlist(gn_res)<=GlobalNull.Test.Alpha)
      GN_NULL_LIST[[b]] = gn_res
    }
  }  
  return_object = list(TP_random_vec = TP_random_vec,
                       FDR_random_vec = FDR_random_vec,
                       TP_most_abundant_vec = TP_most_abundant_vec,
                       FDR_most_abundant_vec = FDR_most_abundant_vec,
                       GN_NULL_LIST = GN_NULL_LIST,
                       CONTAINS_DIFF_ABUNDANT = CONTAINS_DIFF_ABUNDANT)
  return(return_object)
} # End of worker function

#stop('Stopping!')

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
  if(MODE_COMPUTE_RANDOM_SELECT_FDR){
    combined_results_FDR_random = res[[1]]$FDR_random_vec
    combined_results_FDR_most_abundant = res[[1]]$FDR_most_abundant_vec
    if(length(res)>1){
      for(i in 2:length(res)){
        combined_results_FDR_random = c(combined_results_FDR_random,res[[i]]$FDR_random_vec)
        combined_results_FDR_most_abundant = c(combined_results_FDR_most_abundant,res[[i]]$FDR_most_abundant_vec)
      }    
    }  
  }
  if(MODE_COMPUTE_GLOBAL_NULL){
    flatten_GN_NULL_LIST = function(l,contains_diff_abun){
      temp = 1*(unlist(l[[1]]) < GlobalNull.Test.Alpha)
      if(length(l)>1){
        for(i in 2:length(l)){
          temp = temp +1*(unlist(l[[i]]) < GlobalNull.Test.Alpha)  
        }
      }
      return(temp)
    }
    nr_reject_GN = flatten_GN_NULL_LIST(res[[1]]$GN_NULL_LIST,res[[1]]$CONTAINS_DIFF_ABUNDANT)
    if(length(res)>1){
      for(i in 2:length(res)){
        nr_reject_GN = nr_reject_GN + flatten_GN_NULL_LIST(res[[i]]$GN_NULL_LIST,res[[i]]$CONTAINS_DIFF_ABUNDANT)
      }    
    }  
    nr_reject_GN = nr_reject_GN/(NR_REPS_PER_WORKER*NR.WORKERS)
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
  
  print('nr_reject_GN')
  print(nr_reject_GN)
  print('Error in T1E estimation;2*SE')
  print(round(1.96*sqrt(GlobalNull.Test.Alpha*(1-GlobalNull.Test.Alpha)/(NR_REPS_PER_WORKER*NR.WORKERS)),3))
  sink()
  
  cat(paste0('Done scenario ',SCENARIO_ID,' \n\r'))
}
