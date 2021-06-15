#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Script is used to run the ZINB-Wave method on the simulation settings discussed in the paper.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#install package if needed


#load packages and settigns index
library(zinbwave)
library(phyloseq)
library(DESeq2)
source(paste0('REFSIM_GenerateSettings_Index.R')) # scenario library

#Parameters
PARAM_SCENARIOS = 1:47 # Scenarios to be run
PARAM_REPS = 100 #number of repetitions per scenario
PARAM_Q_LEVEL = 0.1 #FDR level

source('zinbwave_imported_functions.R')

# Function for evaluating FDR and average number of true positives for a scenario.
evaluate_scenario_Power_and_FDR_for_ZINBWAVE = function(scenario_id = 1,Q_LEVEL = 0.1,nr_reps = 50){
  #scenario_id = 29
  #Q_LEVEL = 0.1
  #nr_reps = 50
  set.seed(1)
  TP_vec = rep(NA,nr_reps)
  FDR_vec = rep(NA,nr_reps)
  rejected_vec = rep(NA,nr_reps)
  for(r in 1:nr_reps){
    print(paste0('Doing data generation ',r,'/',nr_reps))
    #generate data
    data_generation = REFSIM_generate_setting_wrapper(setting_parameters = REFSIM_SETTINGS_LIST[[scenario_id]])
    
    #remove taxa with zero counts in all samples
    ind_to_remove_from_original =which(apply(data_generation$X,2,sum)==0)
    if(length(ind_to_remove_from_original)>0){
      diff_abundant_as_binary = rep(0,ncol(data_generation$X))
      diff_abundant_as_binary[data_generation$select_diff_abundant] = 1
      
      data_generation$X = data_generation$X[,-ind_to_remove_from_original]
      data_generation$X_unsampled = data_generation$X_unsampled[,-ind_to_remove_from_original]
      diff_abundant_as_binary = diff_abundant_as_binary[-ind_to_remove_from_original]
      data_generation$select_diff_abundant = which(diff_abundant_as_binary == 1)  
    }
    
    
    # Pack as phyloseq object
    physeq = phyloseq(otu_table(t(data_generation$X),taxa_are_rows = TRUE),
                      sample_data(
                        data.frame(
                          grp = c(rep("Y0",sum(data_generation$Y==0)),rep("Y1",sum(data_generation$Y==1))),
                          col.names = paste0('Sample_',1:(nrow(data_generation$X))) #row.names = paste0('Taxon_',1:(ncol(data_generation$X))
                        )
                      )
    )
    
    #used for debugging
    #dim(otu_table(physeq))
    #otu_table(physeq)
    #sample_data(physeq)
    
    #convert to DESEQ2 object and run methods as the eval_functions.R script in in https://github.com/mcalgaro93/sc2meta
    physeq <- normDESeq2(physeq = physeq) 
    
    epsilon = 1e10
    zinbmodel <- zinbFit(Y = physeq@otu_table@.Data, 
                         X = model.matrix(~ physeq@sam_data$grp), K = 0,
                         epsilon = epsilon, commondispersion = TRUE, verbose = FALSE, BPPARAM = SerialParam())
    
    weights <- computeExactWeights(model = zinbmodel,x = physeq@otu_table@.Data)
    colnames(weights) <- colnames(physeq@otu_table)
    rownames(weights) <- rownames(physeq@otu_table)
    
    
    DESeq2_poscounts_zinbwave <- negBinTestDESeq2_zinbweights(physeq, normFacts = "poscounts",weights = weights)
    
    #check the rejected and collect statistics
    rejected = which(p.adjust(DESeq2_poscounts_zinbwave$pValMat[,1],method = 'BH')<=Q_LEVEL)
    
    
    rejected_vec[r] = length(rejected)
    TP = sum(rejected %in% data_generation$select_diff_abundant)
    FP = length(rejected) - TP
    FDR = FP/max(length(rejected),1)
    
    TP_vec[r] = TP
    FDR_vec[r] = FDR
  } # done loop over cases
  
  summary_statistics = c(TP_avg = mean(TP_vec),
          TP_se = sd(TP_vec)/sqrt(nr_reps),
          FDR_avg = mean(FDR_vec),
          FDR_se = sd(FDR_vec)/sqrt(nr_reps),
          rejected_avg = mean(rejected_vec),
          rejected_se = sd(rejected_vec)/sqrt(nr_reps))
  # original_vectors = list(TP_vec = TP_vec,FDR_vec = FDR_vec)
  # ret = list(summary_statistics = summary_statistics,original_vectors = original_vectors)
  return(summary_statistics)
}


# Parallelize over cases
  
  library(doParallel)
  library(doRNG)
  cl = NULL
  rm(cl)
  cl <- makeCluster(parallel::detectCores() - 1)
  registerDoParallel(cl)
  
  
  print(paste0('Starting sim'))
  print(Sys.time())
  Start = Sys.time()
  
  #length(REFSIM_SETTINGS_LIST)
  zinbwave_deseq2_parallel_res <- foreach(i=PARAM_SCENARIOS, .options.RNG=1,.combine = 'rbind') %dorng% {
    
    library(zinbwave)
    library(phyloseq)
    library(DESeq2)
    
    ret = list()
    tryCatch(expr = {ret = evaluate_scenario_Power_and_FDR_for_ZINBWAVE(scenario_id = i,
                                                                        nr_reps = PARAM_REPS,
                                                                        Q_LEVEL = PARAM_Q_LEVEL)})
    return(ret)
  }
  
  stopCluster(cl)
  
  print(paste0('Ending sim'))
  print(Sys.time())
  print('Duration:')
  print(Sys.time() - Start)

#put meaning full row names, by scenarios
rownames(zinbwave_deseq2_parallel_res) = unlist(lapply(REFSIM_SETTINGS_LIST,function(x){x$label}))[PARAM_SCENARIOS] #paste0('Scenario ',PARAM_SCENARIOS)

zinbwave_deseq2_results = as.data.frame(zinbwave_deseq2_parallel_res)
# zinbwave_deseq2_results$TP_sd = zinbwave_deseq2_results$TP_sd/sqrt(PARAM_REPS)
# zinbwave_deseq2_results$FDR_sd = zinbwave_deseq2_results$FDR_sd/sqrt(PARAM_REPS)
# zinbwave_deseq2_results$rejected_sd = zinbwave_deseq2_results$rejected_sd/sqrt(PARAM_REPS)
# names(zinbwave_deseq2_results)[c(2,4,6)] = c("TP_se","FDR_se","rejected_se")

# max(zinbwave_deseq2_results$FDR_se)
# View(zinbwave_deseq2_results)
# 
# sim_res_zinbwave = list()
# sim_res_zinbwave$zinbwave_deseq2_parallel_res = zinbwave_deseq2_parallel_res
# sim_res_zinbwave$zinbwave_deseq2_results = zinbwave_deseq2_results
save(zinbwave_deseq2_results,file = '../../Results/zinbwave_deseq2_results.rdata')
