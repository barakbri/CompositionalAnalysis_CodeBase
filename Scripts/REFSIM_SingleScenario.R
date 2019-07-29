# This is the main simulation script, used for reproducing the results for the simulation in
# section 4 of the manuscript, along with the results in appendicies A, 
# and the data used for appendix B.1 (results for appendix B.1 are compiled via 'REFSIM_Analyze_chance_of_bad_taxa')
# see the read me file on how this script is fitted in the complete pipeline for reproducing results

# This script is for running a single scenario at a time. The following line is commented out,
# since this file is called by REFSIM_MultipleScenario_batch.R with multiple parameters.
# To run this file for a single scenario, comment out the following line.

#SCENARIO_ID = 10

#This script was run on a machine with 96 cores. You can change this number
# if you want to run on a different machine

NR.WORKERS = 100
NR.CORES = 96

MAINDIR = './'


print(paste0('Start Time:'))
print(Sys.time())

#Modes and parameters

#Parameters for multiplicity correction
Q_LEVEL = 0.1
CORRECTION_TYPE_MULTIPLICITY = 'BH'


MODE_RUN_SIMULATION = T # should the simulation be run
MODE_PROCESS_RESULTS = T #should results from the simulation be proccessed (single core - compute Power , FDR, etc.. across sample repetitions)
DEBUG = F # used for printing messages, running low number of reps
DEBUG.SINGLE.CORE = F #run in single core...
DISABLE_ANCOM = F #ANCOM takes time to run. for debug purposes - you can turn off ANCOM

#Load libraries
library(dacomp)
library(foreach)
library(doParallel)
library(doRNG)
library(ancom.R)
library(phyloseq)
library(metagenomeSeq)
library(ALDEx2)
library(Wrench)
library(DESeq2)

#load additional files used
source(paste0('Wilcoxon_TaxaWise.R')) #competitor tests
source(paste0('exactHyperGeometricTest.R'))
source(paste0('REFSIM_GenerateSettings_Index.R')) # scenario library

RESULTS_DIR = paste0("../../Results/")

RNG_SEED = 1
REPS_PER_SETTING = 1*NR.WORKERS 

if(DEBUG){
  NR.WORKERS = 7
  REPS_PER_SETTING = 1*NR.WORKERS
}
# scenarios 22,25 are memory intensive. ALDEx2 could not run 96 instances in parallel, on the C5 amazon machine I used, so i use only half of the cores (so I wouldnt run out of RAM)
if(SCENARIO_ID %in% c(22:25) & !DEBUG){
  NR.CORES = 50
}

NR_REPS_PER_WORKER = ceiling(REPS_PER_SETTING/NR.WORKERS)

#labels for the different methods
METHOD_LABEL_ANCOM = "ANCOM"
METHOD_LABEL_HG = "HG"
METHOD_LABEL_WILCOXON_NAIVE   = "WILCOXON_NAIVE"
METHOD_LABEL_WILCOXON_FLOW    = "WILCOXON_FLOW"
METHOD_LABEL_WILCOXON_PERCENT = "WILCOXON_PERCENT"
METHOD_LABEL_WILCOXON_PAULSON = "WILCOXON_PAULSON"
METHOD_LABEL_ALDEx2_Welch     = "ALDEx2_Welch"
METHOD_LABEL_ALDEx2_Wilcoxon  = "ALDEx2_Wilcoxon"
METHOD_LABEL_Wrench  = "Wrench"

#for the DACOMP method, we will have prefixes
PREFIX_DACOMP_WILCOXON = 'DACOMP,Wilcoxon,'
PREFIX_DACOMP_WELCH = 'DACOMP,Welch,'

PREFIX_DACOMP_rarefaction = 'rarefaction,'
PREFIX_DACOMP_division = 'division,'

PS_used = 1 #psuedo count used

SET_OF_SCENARIOS_ANCOM_ALL_FEATURES = -1 # can be use for running ANCOM with other than the default parameters, see in code below. Currently, not used

#Define the median threshold scores used, along with the oracle methods used. 
df_selection_method = data.frame(MethodName = NA, is_Oracle = NA, Median_SD_Threshold = NA, All_Oracle = NA, RandomSize = NA, stringsAsFactors = F)
df_selection_method[1,]                            = c('S = 1.3'          ,F, 1.3,F,NA)

if(!DEBUG){
 df_selection_method[1,]                            = c('S = 1.1'          ,F, 1.1, F,NA)
 df_selection_method[nrow(df_selection_method) +1,] = c('S = 1.2'          ,F, 1.2, F,NA)
 df_selection_method[nrow(df_selection_method) +1,] = c('S = 1.3'          ,F, 1.3, F,NA)
 df_selection_method[nrow(df_selection_method) +1,] = c('S = 1.3, Oracle'  ,T, 1.3, F,NA)
 df_selection_method[nrow(df_selection_method) +1,] = c('S = 1.4'          ,F, 1.4, F,NA)
 df_selection_method[nrow(df_selection_method) +1,] = c('S = 1.4, Oracle'  ,T, 1.4, F,NA)
}

#Number of methods run
NR_METHODS = 4*(nrow(df_selection_method)) +11 # 9+2 for ANCOMX3 and other competitor
NR_METHODS_ROWS = 4*(nrow(df_selection_method)) +11

#auxilary functions for file names
REFSIM_results_file = function(RESULTS_DIR,SCENARIO_ID){ #to store results
  return(paste0(RESULTS_DIR,'/REFSIM_results_',SCENARIO_ID,'.RData'))
}

REFSIM_warnings_file = function(RESULTS_DIR,SCENARIO_ID){ #to give warnings
  return(paste0(RESULTS_DIR,'/REFSIM_warnings_',SCENARIO_ID,'.txt'))
}

REFSIM_aggregated_results_file = function(RESULTS_DIR,SCENARIO_ID){ # aggregated results, average power, FDR etc...
  return(paste0(RESULTS_DIR,'/REFSIM_aggregated_results_',SCENARIO_ID,'.RData'))
}

REFSIM_aggregated_results_file_sd = function(RESULTS_DIR,SCENARIO_ID){ # sd of aggregates
  return(paste0(RESULTS_DIR,'/REFSIM_aggregated_results_sd_',SCENARIO_ID,'.RData'))
}

# function for scenario based worker. this is the main function being run

Worker_Function = function(core_nr){
  library(dacomp) #load the framework
  library(ancom.R)
  library(phyloseq)
  library(metagenomeSeq)
  library(ALDEx2)
  library(Wrench)
  library(DESeq2)
  
  source(paste0('Wilcoxon_TaxaWise.R'))
  source(paste0('exactHyperGeometricTest.R'))
  
  #this is the function from the flow cytometry paper. See comment below - they had a round error I had to fix
  rarefy_even_sampling_depth <- function(cnv_corrected_abundance_table, cell_counts_table) 
  {
    try(if(all(row.names(cnv_corrected_abundance_table) == row.names(cell_counts_table))==FALSE) stop("Cnv_corrected_abundance_table and cell_counts_table do not have the same sample-names, Please check!"))
    cnv_corrected_abundance_table = ceiling(cnv_corrected_abundance_table) # data values are rounded up in order to make use of integer values during the calculations
    cell_counts_table = t(cell_counts_table[order(row.names(cnv_corrected_abundance_table)),]) # make sure the order of the samples is the same in both files  
    sample_sizes = rowSums(cnv_corrected_abundance_table) # sample size of each sample (total nr of reads)
    sampling_depths = sample_sizes / cell_counts_table # sampling depth of each sample (total nr of reads divided by the cell count)
    minimum_sampling_depth = min(sampling_depths) # minimum of all sampling depths
    rarefy_to = cell_counts_table * minimum_sampling_depth # nr of reads to rarefy in each sample in order to get to an even sampling depth over all samples
    cnv_corrected_abundance_table_phyloseq = otu_table(cnv_corrected_abundance_table, taxa_are_rows = FALSE) # convert to phyloseq otutable
    rarefied_matrix=matrix(nrow = nrow(cnv_corrected_abundance_table_phyloseq), ncol = ncol(cnv_corrected_abundance_table_phyloseq), dimnames = list(rownames(cnv_corrected_abundance_table_phyloseq), colnames(cnv_corrected_abundance_table_phyloseq)))
    for (i in 1:nrow(cnv_corrected_abundance_table_phyloseq))
    {
      # the rarefaction depth is taken to be one counts less (a single sequence) from the one they set.
      # It seems they have a rounding error and you can get a generated dataset where they try to rarefy
      # one of the samples to a larger depth than it originally has (at most one read bigger, due to rounding)
      # The "minus one" fixes that
      x <- rarefy_even_depth(cnv_corrected_abundance_table_phyloseq[i,], sample.size = rarefy_to[i]-1, rngseed = FALSE, replace = FALSE, trimOTUs = F, verbose = FALSE)
      rarefied_matrix[i,] = x
    }
    normalised_rarefied_matrix = rarefied_matrix/rowSums(rarefied_matrix)
    QMP = normalised_rarefied_matrix*cell_counts_table[1,]
    return(QMP)
  }
  
  #data structure for results
  current_row = 1
  RESULTS_LENGTH = NR_REPS_PER_WORKER*NR_METHODS_ROWS
  results = data.frame(setting_id = rep(NA,RESULTS_LENGTH),
                       methodlabel = rep(NA,RESULTS_LENGTH),
                       rejected = rep(NA,RESULTS_LENGTH),
                       tp = rep(NA,RESULTS_LENGTH),
                       fp = rep(NA,RESULTS_LENGTH),
                       nr_bad_reference = rep(NA,RESULTS_LENGTH) #relevant only for DACOMP
                       )
  
  for(b in 1:NR_REPS_PER_WORKER){ #for each repetition, in the current scenario:
    
    if(DEBUG)
      cat(paste0('Case ',b,' out of ',NR_REPS_PER_WORKER,'\n\r'))
    #generate data:
    current_setting_generator = REFSIM_SETTINGS_LIST[[SCENARIO_ID]]
    data = REFSIM_generate_setting_wrapper(current_setting_generator)
    
    #run our methods, with both Welch and Wilcoxon tests, and all reference selections methods:
    for(m in 1: nrow(df_selection_method)){ 
      if(DEBUG)
        cat(paste0('Method ',m,' out of ',NR_METHODS,'\n\r'))
      
      # get reference selection method details:
      
      ref_select_method_Label             = df_selection_method$MethodName[m]
      ref_select_method_Median_SD_Threshold  = as.numeric(df_selection_method$Median_SD_Threshold[m])
      ref_select_method_is_Oracle         = as.logical(df_selection_method$is_Oracle[m])
      ref_select_method_All_Oracle        = as.logical(df_selection_method$All_Oracle[m])
      ref_select_method_RandomSize        = as.numeric(df_selection_method$RandomSize[m])
      
      # either select references, or for an oracle method, select with oracle knowledge
      Real_Nulls = (1:(dim(data$X)[2]))
      if(length(data$select_diff_abundant) >0 & length(which(data$select_diff_abundant == 0))==0)
        Real_Nulls = (1:(dim(data$X)[2]))[-data$select_diff_abundant]
      
      Selected_References = NULL
      if(!is.na(ref_select_method_RandomSize)){
        Selected_References = sample(1:(dim(data$X)[2]),size = ref_select_method_RandomSize,replace = F)
      }else if (ref_select_method_All_Oracle){
        Selected_References = Real_Nulls
      }else{
        to_select_from = (1:(dim(data$X)[2]))
        if(ref_select_method_is_Oracle)
          to_select_from = Real_Nulls
        
        ref_select = dacomp::dacomp.select_references(X = data$X,
                                                      median_SD_threshold = ref_select_method_Median_SD_Threshold,
                                                      minimal_TA = 10,
                                                      maximal_TA = 200,
                                                      Pseudo_Count_used = PS_used,
                                                      verbose = DEBUG.SINGLE.CORE,
                                                      select_from = to_select_from)
        
        Selected_References = ref_select$selected_references
                                                      
      }
      
      # run both wilcoxon and welch tests
      res_wilcoxon = dacomp.test(X = data$X,
                                 y = data$Y,
                                 ind_reference_taxa = Selected_References,
                                 test = DACOMP.TEST.NAME.WILCOXON,
                                 q = Q_LEVEL,
                                 nr_perm = 1/(Q_LEVEL/(ncol(data$X))),compute_ratio_normalization = T,verbose = DEBUG.SINGLE.CORE)
      
      res_welch = dacomp.test(X = data$X,
                                 y = data$Y,
                                 ind_reference_taxa = Selected_References,
                                 test = DACOMP.TEST.NAME.WELCH_LOGSCALE,
                                 q = Q_LEVEL,
                                 nr_perm = 1/(Q_LEVEL/(ncol(data$X))),compute_ratio_normalization = T,verbose = DEBUG.SINGLE.CORE)
      
      #number of diff.abun. taxa that have entered the reference set by error
      bad_select = length(which(ref_select$selected_references %in% data$select_diff_abundant))
      
      #auxiliary function for wrapping the results
      construct_results_row = function(p.values.test,ref_select_method_Label){
        rejected = which(p.adjust(p.values.test,method = CORRECTION_TYPE_MULTIPLICITY) <= Q_LEVEL) #pvals
        true_positive = length(which(rejected %in% data$select_diff_abundant))
        false_positive = length(rejected) - true_positive
        ret = c(SCENARIO_ID,
                ref_select_method_Label,
                length(rejected),
                true_positive,
                false_positive,
                bad_select)
        return(ret)
      }
      
      #docuemnt our results for wilcoxon/Welch, and for subsampling/ normalization by ratio.
      results[current_row,] = construct_results_row(res_wilcoxon$p.values.test ,
                                                    paste0(PREFIX_DACOMP_WILCOXON,PREFIX_DACOMP_rarefaction,ref_select_method_Label)
                                                    ); current_row = current_row + 1
      
      results[current_row,] = construct_results_row(res_wilcoxon$p.values.test.ratio.normalization ,
                                                    paste0(PREFIX_DACOMP_WILCOXON,PREFIX_DACOMP_division,ref_select_method_Label)
                                                    ); current_row = current_row + 1
      
      results[current_row,] = construct_results_row(res_welch$p.values.test ,
                                                    paste0(PREFIX_DACOMP_WELCH,PREFIX_DACOMP_rarefaction,ref_select_method_Label)
                                                    ); current_row = current_row + 1
      
      results[current_row,] = construct_results_row(res_welch$p.values.test.ratio.normalization ,
                                                    paste0(PREFIX_DACOMP_WELCH,PREFIX_DACOMP_division,ref_select_method_Label)
                                                    ); current_row = current_row + 1
      

      #run variant with hypergeometric test. Only valid if there is now over dispersion of counts      
      res_HG = exactHyperGeometricTest(data$X,data$Y,Selected_References)
      rejected = which(p.adjust(res_HG,method = CORRECTION_TYPE_MULTIPLICITY) <= Q_LEVEL) #pvals
      true_positive = length(which(rejected %in% data$select_diff_abundant))
      false_positive = length(rejected) - true_positive
      results[current_row,] = c(SCENARIO_ID,paste0(METHOD_LABEL_HG,',', ref_select_method_Label),length(rejected),true_positive,
                                false_positive, bad_select) ; current_row = current_row + 1
      
      
    }# end of loop over methods
    
    #run marginal wilcoxon tests
    if(DEBUG)
      cat(paste0('Computing Wilcoxon \n\r'))
    
    #without normalization
    res = wilcoxon_taxa_wise(data$X,data$Y)
    rejected = which(p.adjust(res$p.values,method = CORRECTION_TYPE_MULTIPLICITY) <= Q_LEVEL)
    true_positive = length(which(rejected %in% data$select_diff_abundant))
    false_positive = length(rejected) - true_positive
    results[current_row,] = c(SCENARIO_ID,METHOD_LABEL_WILCOXON_NAIVE,length(rejected),true_positive,false_positive,0) ; current_row = current_row + 1
    
    #TSS normalization
    res = wilcoxon_taxa_wise(data$X,data$Y,normalize = T,normalize.P = 1)
    rejected = which(p.adjust(res$p.values,method = CORRECTION_TYPE_MULTIPLICITY) <= Q_LEVEL)
    true_positive = length(which(rejected %in% data$select_diff_abundant))
    false_positive = length(rejected) - true_positive
    results[current_row,] = c(SCENARIO_ID,METHOD_LABEL_WILCOXON_PERCENT,length(rejected),true_positive,false_positive,0) ; current_row = current_row + 1
    
    #CSS normalization
    CSS_normalized = t(cumNormMat(t(data$X)))
    res = wilcoxon_taxa_wise(CSS_normalized,data$Y,normalize = F)
    rejected = which(p.adjust(res$p.values,method = CORRECTION_TYPE_MULTIPLICITY) <= Q_LEVEL)
    true_positive = length(which(rejected %in% data$select_diff_abundant))
    false_positive = length(rejected) - true_positive
    results[current_row,] = c(SCENARIO_ID,METHOD_LABEL_WILCOXON_PAULSON,length(rejected),true_positive,false_positive,0) ; current_row = current_row + 1
    
    # for scenarios with simulated flow cytometry- run van de putte method
    if(class(current_setting_generator) %in% c(REFSIM_GUT_TYPE_SCENARIO_DEF) & "Total_Original_Counts" %in% names(data)){
      
      
      if(DEBUG)
        cat(paste0('Computing Wilcoxon over Flow-Cytometry adjusted values \n\r'))
      
      X_corrected = data$X
      row.names(X_corrected) = 1:(nrow(X_corrected))
      Average_Cell_Count_matrix = matrix(data$Total_Original_Counts,ncol = 1)
      rownames(Average_Cell_Count_matrix) = rownames(X_corrected)
      X_corrected = rarefy_even_sampling_depth(X_corrected, Average_Cell_Count_matrix)
      
      res = wilcoxon_taxa_wise(X_corrected,data$Y)
      rejected = which(p.adjust(res$p.values,method = CORRECTION_TYPE_MULTIPLICITY) <= Q_LEVEL)
      true_positive = length(which(rejected %in% data$select_diff_abundant))
      false_positive = length(rejected) - true_positive
      results[current_row,] = c(SCENARIO_ID,METHOD_LABEL_WILCOXON_FLOW,length(rejected),true_positive,false_positive,0) ; current_row = current_row + 1
    }
    
    #run ANCOM, if required    
    if(DEBUG)
      cat(paste0('Computing ANCOM \n\r'))
    
    mat = (as.data.frame(data$X))
    mat[,ncol(mat) + 1] =  as.factor(data$Y)
    names(mat) = c(1:ncol(data$X),'Y')
    if(!DISABLE_ANCOM){
      for(ai in 1:3){ #iterate over ANCOM parameters
        #run ANCOM with the default parameter, or if it is a scenario that requries checking all configurations
        if(ai == 3 | SCENARIO_ID %in% SET_OF_SCENARIOS_ANCOM_ALL_FEATURES){ 
          ancom.out <- ANCOM(mat, sig = Q_LEVEL, multcorr = ai, repeated=FALSE)
          ANCOM_Rejected = ancom.out$detected
          for(i in 1:length(ANCOM_Rejected)){
            ANCOM_Rejected[i] = substr(ANCOM_Rejected[i],start = 2,stop = nchar(ANCOM_Rejected[i]))
          }
          ANCOM_Rejected = as.numeric(ANCOM_Rejected)
          rejected = ((ANCOM_Rejected))
          true_positive = length(which(rejected %in% data$select_diff_abundant))
          false_positive = length(rejected) - true_positive
          current_label = METHOD_LABEL_ANCOM
          if(ai<=2)
            current_label = paste0(METHOD_LABEL_ANCOM,'_',ai)
          results[current_row,] = c(SCENARIO_ID,current_label,length(rejected),true_positive,false_positive,0) ; current_row = current_row + 1
        }
      }
    }else{
      results[current_row,] = c(SCENARIO_ID,METHOD_LABEL_ANCOM,-1,-1,-1,-1) ; current_row = current_row + 1  
    }
    
    #ALDEx2 (note that we alter the data object here, we remove taxa which are empty, and keep track of how the inices chagned for computing Power and FDR):
    indicies_not_empty = which(apply(data$X>0,2,sum)>0) # to keep
    m_original = ncol(data$X)
    flags_H1 = rep(0,m_original)
    flags_H1[data$select_diff_abundant] = 1
    flags_H1 = flags_H1[indicies_not_empty]
    data$X = data$X[,indicies_not_empty] # keep only relevant columns
    data$select_diff_abundant = which(flags_H1==1) #the vector of indices of diff.abun. taxa, under the new indexing
    
    #run iqlr normalizaed aldex 2
    aldex.res.iqlr <- aldex(t(data$X), as.character(data$Y), mc.samples=128, denom="iqlr",
                            test="t", effect=FALSE,verbose = DEBUG.SINGLE.CORE)
    
    # record both wilcoxon and welch variants
    rejected.iqlr.wi = p.adjust(aldex.res.iqlr$wi.ep, method = CORRECTION_TYPE_MULTIPLICITY) <= Q_LEVEL
    rejected.iqlr.we = p.adjust(aldex.res.iqlr$we.ep, method = CORRECTION_TYPE_MULTIPLICITY) <= Q_LEVEL
    
    ALDEx2.Nr.Rejected.Wi = sum(rejected.iqlr.wi)
    ALDEx2.Nr.Rejected.We = sum(rejected.iqlr.we)
    ALDEx2.Nr.TP.Wi = sum(which(rejected.iqlr.wi) %in% data$select_diff_abundant)
    ALDEx2.Nr.TP.We = sum(which(rejected.iqlr.we) %in% data$select_diff_abundant)
    ALDEx2.Nr.FP.Wi = ALDEx2.Nr.Rejected.Wi - ALDEx2.Nr.TP.Wi
    ALDEx2.Nr.FP.We = ALDEx2.Nr.Rejected.We - ALDEx2.Nr.TP.We
    
    results[current_row,] = c(SCENARIO_ID,METHOD_LABEL_ALDEx2_Wilcoxon, ALDEx2.Nr.Rejected.Wi, ALDEx2.Nr.TP.Wi ,ALDEx2.Nr.FP.Wi,0) ; current_row = current_row + 1
    results[current_row,] = c(SCENARIO_ID,METHOD_LABEL_ALDEx2_Welch, ALDEx2.Nr.Rejected.We, ALDEx2.Nr.TP.We ,ALDEx2.Nr.FP.We,0) ; current_row = current_row + 1
    
    #Wrench - note that we removed otus with zeroes before
    W <- wrench( t(data$X), condition=data$Y  ) #normalize. Note that this function throws off warnings if the scenario is too sparse..
    compositionalFactors <- W$ccf
    normalizationFactors <- W$nf
    
    library(DESeq2) #test via deseq 2.
    deseq.obj <- DESeq2::DESeqDataSetFromMatrix(countData = t(data$X),
                                                DataFrame(Y = factor(data$Y)),
                                                ~ Y )
    sizeFactors(deseq.obj) <- normalizationFactors
    deseq2_res = DESeq2::DESeq(deseq.obj)
    deseq2_res2 = results(deseq2_res)
    
    
    wrench_rejected = which(p.adjust(deseq2_res2$pvalue,method = CORRECTION_TYPE_MULTIPLICITY) <= Q_LEVEL)
    wrench_nr_rejected = length(wrench_rejected)
    wrench_nr_tp = sum(wrench_rejected%in% data$select_diff_abundant)
    wrench_nr_fp = wrench_nr_rejected - wrench_nr_tp
    results[current_row,] = c(SCENARIO_ID,METHOD_LABEL_Wrench, wrench_nr_rejected, wrench_nr_tp ,wrench_nr_fp,0) ;
    current_row = current_row + 1
    
  }# end of reps per worker  
  return_object = list()
  return_object$results = results

  return(return_object)
} # End of worker function

#Code snippet for testing the worker function
# if(F){
#   results_demo = Worker_Function(1)
#   results_demo$results
#   
#   hist(results_demo$lambda_selected_list[[1]])
# }

#if needed to run simulation, handle:
if(MODE_RUN_SIMULATION){
  cat(paste0('Running simulation for scenario ID ',SCENARIO_ID,'\n\r'))
  Start = Sys.time()
  res=NULL
  rm(res)
  res=NULL
  if(!DEBUG.SINGLE.CORE){ #run on multiple cores
    cat(paste0(' Using multiple cores \n\r'))
    cl = makeCluster(NR.CORES)
    registerDoParallel(cl)
    res = foreach(core = 1:NR.WORKERS, .options.RNG = RNG_SEED) %dorng% {Worker_Function(core)} #.options.RNG=c(1:NR.WORKERS) #%dorng%
    stopCluster(cl)
    cat(paste0(' Run done \n\r'))
  }else{ #run on single core, can see warnings and debug 
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
  #save results from this scenario (with all repetitions) to file
  filename = REFSIM_results_file(RESULTS_DIR,SCENARIO_ID)
  save(res,file = filename)
}

sink(file = REFSIM_warnings_file(RESULTS_DIR,SCENARIO_ID))
warnings()
sink()

#process results if needed
if(MODE_PROCESS_RESULTS){
  #combine results across all cores
  cat(paste0('Combining results for scenario ',SCENARIO_ID,' \n\r'))
  filename = REFSIM_results_file(RESULTS_DIR,SCENARIO_ID)
  load(file = filename) #=> res
  combined_results = res[[1]]$results
  if(length(res)>1){
    for(i in 2:length(res)){
      combined_results = rbind(combined_results,res[[i]]$results)
    }    
  }
  combined_results$rejected = as.numeric(combined_results$rejected)
  combined_results$tp = as.numeric(combined_results$tp)
  combined_results$fp = as.numeric(combined_results$fp)
  combined_results$setting_id = as.numeric(combined_results$setting_id)
  combined_results$nr_bad_reference = as.numeric(combined_results$nr_bad_reference)
  
  # compute mean number of rejections, fdr across scenarios. Also compute the SD of population values (not SEs, need to divide by sqrt(reps) for final answer)
  rejected_for_fdr = combined_results$rejected
  if(length(which(rejected_for_fdr==0))>0)
    rejected_for_fdr[rejected_for_fdr==0] = 1
  combined_results$fdr = combined_results$fp / rejected_for_fdr
  
  
  aggregated_results = aggregate(x= combined_results[,-2],
                       by = list(methodlabel = as.factor(combined_results$methodlabel)),
                       FUN = mean)
  filename = REFSIM_aggregated_results_file(RESULTS_DIR,SCENARIO_ID)
  save(aggregated_results,file = filename)
  
  aggregated_results_sd = aggregate(x= combined_results[,-2],
                                 by = list(methodlabel = as.factor(combined_results$methodlabel)),
                                 FUN = sd)
  filename_SD = REFSIM_aggregated_results_file_sd(RESULTS_DIR,SCENARIO_ID)
  save(aggregated_results_sd,file = filename_SD)
  
  cat(paste0('Done scenario ',SCENARIO_ID,' \n\r'))
}
