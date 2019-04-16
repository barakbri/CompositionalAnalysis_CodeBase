NR.WORKERS = 72 #63#72-1#63
#SCENARIO_ID = 25
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
DISABLE_ANCOM = F
DIAG_PLOTS_SELECTION = F

library(subzero)
library(foreach)
library(doParallel)
library(doRNG)
library(ancom.R)
library(phyloseq)
library(metagenomeSeq)

source(paste0('Wilcoxon_TaxaWise.R'))
source(paste0('exactHyperGeometricTest.R'))
source(paste0('REFSIM_GenerateSettings_Index.R'))

RESULTS_DIR = paste0("../../Results/")

RNG_SEED = 1
REPS_PER_SETTING = 1*NR.WORKERS 

if(DEBUG){
  
  NR.WORKERS = 7
  REPS_PER_SETTING = 1*NR.WORKERS
  
}

NR_REPS_PER_WORKER = ceiling(REPS_PER_SETTING/NR.WORKERS)
METHOD_LABEL_ANCOM = "ANCOM"
METHOD_LABEL_HG = "HG"
METHOD_LABEL_WILCOXON_NAIVE   = "WILCOXON_NAIVE"
METHOD_LABEL_WILCOXON_qPCR    = "WILCOXON_qPCR"
METHOD_LABEL_WILCOXON_PERCENT = "WILCOXON_PERCENT"
METHOD_LABEL_WILCOXON_PAULSON = "WILCOXON_PAULSON"

PS_used = 1
Use_Median = F
SET_OF_SCENARIOS_ANCOM_ALL_FEATURES = -1

df_selection_method = data.frame(MethodName = NA, is_Oracle = NA, Median_SD_Threshold = NA, All_Oracle = NA, RandomSize = NA, stringsAsFactors = F)
df_selection_method[1,]                            = c('S = 1.3'          ,F, 1.3,F,NA)

if(!DEBUG){
 #df_selection_method[1,] = c('S = 1.3'          ,F, 1.3, F,NA)
 #df_selection_method[nrow(df_selection_method) +1,] = c('S = 1.4, Oracle'    ,T, 1.4, F,NA)
 
 df_selection_method[1,]                            = c('S = 1.1'          ,F, 1.1, F,NA)
 df_selection_method[nrow(df_selection_method) +1,] = c('S = 1.2'          ,F, 1.2, F,NA)
 df_selection_method[nrow(df_selection_method) +1,] = c('S = 1.3'          ,F, 1.3, F,NA)
 df_selection_method[nrow(df_selection_method) +1,] = c('S = 1.4'          ,F, 1.4, F,NA)
 df_selection_method[nrow(df_selection_method) +1,] = c('S = 1.4, Oracle'    ,T, 1.4, F,NA)
}

NR_METHODS = 3*(nrow(df_selection_method)) +7 # 7 for ANCOMX3 and WILCOXON + percent + Paulson + qPCR,  for dfdr with "mean log +1"
NR_METHODS_ROWS = 3*(nrow(df_selection_method)) +7 

#auxilary functions
REFSIM_results_file = function(RESULTS_DIR,SCENARIO_ID){
  return(paste0(RESULTS_DIR,'/REFSIM_results_',SCENARIO_ID,'.RData'))
}

REFSIM_aggregated_results_file = function(RESULTS_DIR,SCENARIO_ID){
  return(paste0(RESULTS_DIR,'/REFSIM_aggregated_results_',SCENARIO_ID,'.RData'))
}

REFSIM_aggregated_results_file_sd = function(RESULTS_DIR,SCENARIO_ID){
  return(paste0(RESULTS_DIR,'/REFSIM_aggregated_results_sd_',SCENARIO_ID,'.RData'))
}

#QPCR testing function



# function for scenario based worker

Worker_Function = function(core_nr){
  library(subzero)
  library(ancom.R)
  library(phyloseq)
  library(metagenomeSeq)
  
  source(paste0('Wilcoxon_TaxaWise.R'))
  source(paste0('exactHyperGeometricTest.R'))
  #source(paste0('../Packages/QMP-master/QMP-master/QMP.R'))
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
      #need to document the -1 here
      x <- rarefy_even_depth(cnv_corrected_abundance_table_phyloseq[i,], sample.size = rarefy_to[i]-1, rngseed = FALSE, replace = FALSE, trimOTUs = F, verbose = FALSE)
      rarefied_matrix[i,] = x
    }
    normalised_rarefied_matrix = rarefied_matrix/rowSums(rarefied_matrix)
    QMP = normalised_rarefied_matrix*cell_counts_table[1,]
    return(QMP)
  }
  
  current_row = 1
  RESULTS_LENGTH = NR_REPS_PER_WORKER*NR_METHODS_ROWS
  results = data.frame(setting_id = rep(NA,RESULTS_LENGTH),
                       methodlabel = rep(NA,RESULTS_LENGTH),
                       rejected = rep(NA,RESULTS_LENGTH),
                       tp = rep(NA,RESULTS_LENGTH),
                       fp = rep(NA,RESULTS_LENGTH),
                       nr_bad_reference = rep(NA,RESULTS_LENGTH)
                       )
  
  for(b in 1:NR_REPS_PER_WORKER){
    
    if(DEBUG)
      cat(paste0('Case ',b,' out of ',NR_REPS_PER_WORKER,'\n\r'))
    current_setting_generator = REFSIM_SETTINGS_LIST[[SCENARIO_ID]]
    data = REFSIM_generate_setting_wrapper(current_setting_generator)
    
    for(m in 1: nrow(df_selection_method)){
      if(DEBUG)
        cat(paste0('Method ',m,' out of ',NR_METHODS,'\n\r'))
      
      
      #MethodName = NA, is_Oracle = NA, target_prevalence = NA, Min_Prev = NA, All_Oracle = NA, RandomSize = NA
      ref_select_method_Label             = df_selection_method$MethodName[m]
      ref_select_method_Median_SD_Threshold  = as.numeric(df_selection_method$Median_SD_Threshold[m])
      ref_select_method_is_Oracle         = as.logical(df_selection_method$is_Oracle[m])
      ref_select_method_All_Oracle        = as.logical(df_selection_method$All_Oracle[m])
      ref_select_method_RandomSize        = as.numeric(df_selection_method$RandomSize[m])
      
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
        
        ref_select = select.references.Median.SD.Threshold(X = data$X,
                                                           median_SD_threshold = ref_select_method_Median_SD_Threshold,
                                                           minimal_TA = 10,maximal_TA = 200,select_from = to_select_from,Psuedo_Count_used = PS_used,factor_by_Median_Score = Use_Median)
        
        Selected_References = ref_select$selected_references
                                                      
      }
      
      res = subzero.dfdr.test(X = data$X, y = data$Y,
                                           nr_reference_taxa = Selected_References, verbose = T, nr_perm = 1/(Q_LEVEL/(ncol(data$X))),
                              q = Q_LEVEL,disable_DS.FDR = F)
      
      if(DIAG_PLOTS_SELECTION & !ref_select_method_is_Oracle & !ref_select_method_All_Oracle){
        ind_plot = function(res,data,ref_select,filename_1,filename_2){
          pval_H1_test = res$p.values.test
          pval_H1_ref = res$p.values.ref
          pval_H0_test = res$p.values.test
          pval_H0_ref = res$p.values.ref
          
          pval_H1_test[-data$select_diff_abundant] = NA
          pval_H1_ref[-data$select_diff_abundant] = NA
          pval_H0_test[data$select_diff_abundant] = NA
          pval_H0_ref[data$select_diff_abundant] = NA
          png(filename = filename_1,width = 800,height = 800)
          plot(pval_H0_test,ref_select$scores,xlab = "Pvalue of null hyptheses outside of reference set",ylab = 'reference score',main = 'reference score by Pvalue for null hypotheses tested')  
          dev.off()
          png(filename = filename_2,width = 800,height = 800)
          plot(pval_H0_ref,ref_select$scores,xlab = "Pvalue of null hyptheses inside reference set",ylab = 'reference score',main = 'reference score by Pvalue for null hypotheses not tested')  
          dev.off()
        }
        plots_dir = paste0(RESULTS_DIR,'diagnostic_plots/')
        if(!dir.exists(plots_dir))
          dir.create(plots_dir)
        ind_plot(res,data,ref_select,
                 paste0(plots_dir,'diag_Scenario_',SCENARIO_ID,'_core_',core_nr,'_rep_',b,'_method_',m,'_test.png'),
                 paste0(plots_dir,'diag_Scenario_',SCENARIO_ID,'_core_',core_nr,'_rep_',b,'_method_',m,'_ref.png'))  
      }
      
      
      bad_select = length(which(ref_select$selected_references %in% data$select_diff_abundant))
      
      pvals = res$p.values
      pvals[Selected_References] = NA 
      pvals_ref = res$p.values
      pvals_ref[-Selected_References] = NA 
      rejected = which(p.adjust(pvals,method = CORRECTION_TYPE_SUBZERO) < Q_LEVEL) #pvals
      
      true_positive = length(which(rejected %in% data$select_diff_abundant))
      false_positive = length(rejected) - true_positive
      
      
      
      results[current_row,] = c(SCENARIO_ID,ref_select_method_Label,length(rejected),true_positive,
                                false_positive, bad_select) ; current_row = current_row + 1
      
      res_HG = exactHyperGeometricTest(data$X,data$Y,Selected_References)
      rejected = which(p.adjust(res_HG,method = CORRECTION_TYPE_SUBZERO) < Q_LEVEL) #pvals
      true_positive = length(which(rejected %in% data$select_diff_abundant))
      false_positive = length(rejected) - true_positive
      results[current_row,] = c(SCENARIO_ID,paste0(METHOD_LABEL_HG,',', ref_select_method_Label),length(rejected),true_positive,
                                false_positive, bad_select) ; current_row = current_row + 1
      
      
    }# end of loop over methods
    
    if(DEBUG)
      cat(paste0('Computing Wilcoxon \n\r'))
    
    res = wilcoxon_taxa_wise(data$X,data$Y)
    rejected = which(p.adjust(res$p.values,method = CORRECTION_TYPE_WILCOXON) < Q_LEVEL)
    true_positive = length(which(rejected %in% data$select_diff_abundant))
    false_positive = length(rejected) - true_positive
    results[current_row,] = c(SCENARIO_ID,METHOD_LABEL_WILCOXON_NAIVE,length(rejected),true_positive,false_positive,0) ; current_row = current_row + 1
    
    res = wilcoxon_taxa_wise(data$X,data$Y,normalize = T,normalize.P = 1)
    rejected = which(p.adjust(res$p.values,method = CORRECTION_TYPE_WILCOXON) < Q_LEVEL)
    true_positive = length(which(rejected %in% data$select_diff_abundant))
    false_positive = length(rejected) - true_positive
    results[current_row,] = c(SCENARIO_ID,METHOD_LABEL_WILCOXON_PERCENT,length(rejected),true_positive,false_positive,0) ; current_row = current_row + 1
    
    CSS_normalized = t(cumNormMat(t(data$X)))
    res = wilcoxon_taxa_wise(CSS_normalized,data$Y,normalize = F)
    rejected = which(p.adjust(res$p.values,method = CORRECTION_TYPE_WILCOXON) < Q_LEVEL)
    true_positive = length(which(rejected %in% data$select_diff_abundant))
    false_positive = length(rejected) - true_positive
    results[current_row,] = c(SCENARIO_ID,METHOD_LABEL_WILCOXON_PAULSON,length(rejected),true_positive,false_positive,0) ; current_row = current_row + 1
    
    
    if(class(current_setting_generator) %in% c(REFSIM_QPCR_TYPE_SCENARIO_DEF) & "Total_Original_Counts" %in% names(data)){
      
      
      if(DEBUG)
        cat(paste0('Computing Wilcoxon over qPCR adjusted values \n\r'))
      
      X_corrected = data$X
      row.names(X_corrected) = 1:(nrow(X_corrected))
      Average_Cell_Count_matrix = matrix(data$Total_Original_Counts,ncol = 1)
      rownames(Average_Cell_Count_matrix) = rownames(X_corrected)
      X_corrected = rarefy_even_sampling_depth(X_corrected, Average_Cell_Count_matrix)
      
      res = wilcoxon_taxa_wise(X_corrected,data$Y)
      rejected = which(p.adjust(res$p.values,method = CORRECTION_TYPE_WILCOXON) < Q_LEVEL)
      true_positive = length(which(rejected %in% data$select_diff_abundant))
      false_positive = length(rejected) - true_positive
      results[current_row,] = c(SCENARIO_ID,METHOD_LABEL_WILCOXON_qPCR,length(rejected),true_positive,false_positive,0) ; current_row = current_row + 1
      
    }
    
        
    if(DEBUG)
      cat(paste0('Computing ANCOM \n\r'))
    
    mat = (as.data.frame(data$X))
    mat[,ncol(mat) + 1] =  as.factor(data$Y)
    names(mat) = c(1:ncol(data$X),'Y')
    if(!DISABLE_ANCOM){
      for(ai in 1:3){
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
    
  }  
  return_object = list()
  return_object$results = results

  return(return_object)
} # End of worker function

#Tester for worker function
if(F){
  results_demo = Worker_Function(1)
  results_demo$results
  
  hist(results_demo$lambda_selected_list[[1]])
}

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
