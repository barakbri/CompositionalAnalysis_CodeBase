NR.WORKERS = 96 #63#72-1#63
#SCENARIO_ID = 11
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

source(paste0('Wilcoxon_TaxaWise.R'))
source(paste0('exactHyperGeometricTest.R'))
source(paste0('REFSIM_GenerateSettings_Index.R'))

RESULTS_DIR = paste0("../../Results/")

RNG_SEED = 1
REPS_PER_SETTING = 1*NR.WORKERS 

if(DEBUG){
  NR.WORKERS = 1
  REPS_PER_SETTING = 1*NR.WORKERS
}

if(SCENARIO_ID %in% c(22:25) & !DEBUG){
  NR.WORKERS = 96/2
  REPS_PER_SETTING = 2*NR.WORKERS
}

NR_REPS_PER_WORKER = ceiling(REPS_PER_SETTING/NR.WORKERS)
METHOD_LABEL_ANCOM = "ANCOM"
METHOD_LABEL_HG = "HG"
METHOD_LABEL_WILCOXON_NAIVE   = "WILCOXON_NAIVE"
METHOD_LABEL_WILCOXON_qPCR    = "WILCOXON_qPCR"
METHOD_LABEL_WILCOXON_PERCENT = "WILCOXON_PERCENT"
METHOD_LABEL_WILCOXON_PAULSON = "WILCOXON_PAULSON"
METHOD_LABEL_ALDEx2_Welch     = "ALDEx2_Welch"
METHOD_LABEL_ALDEx2_Wilcoxon  = "ALDEx2_Wilcoxon"
METHOD_LABEL_Wrench  = "Wrench"

PREFIX_DACOMP_WILCOXON = 'DACOMP,Wilcoxon,'
PREFIX_DACOMP_WELCH = 'DACOMP,Welch,'

PREFIX_DACOMP_rarefaction = 'rarefaction,'
PREFIX_DACOMP_division = 'division,'

PS_used = 1
Use_Median = F
SET_OF_SCENARIOS_ANCOM_ALL_FEATURES = -1

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

NR_METHODS = 4*(nrow(df_selection_method)) +11 # 9+2 for ANCOMX3 and other competitor
NR_METHODS_ROWS = 4*(nrow(df_selection_method)) +11

#auxilary functions
REFSIM_results_file = function(RESULTS_DIR,SCENARIO_ID){
  return(paste0(RESULTS_DIR,'/REFSIM_results_',SCENARIO_ID,'.RData'))
}

REFSIM_warnings_file = function(RESULTS_DIR,SCENARIO_ID){
  return(paste0(RESULTS_DIR,'/REFSIM_warnings_',SCENARIO_ID,'.txt'))
}

REFSIM_aggregated_results_file = function(RESULTS_DIR,SCENARIO_ID){
  return(paste0(RESULTS_DIR,'/REFSIM_aggregated_results_',SCENARIO_ID,'.RData'))
}

REFSIM_aggregated_results_file_sd = function(RESULTS_DIR,SCENARIO_ID){
  return(paste0(RESULTS_DIR,'/REFSIM_aggregated_results_sd_',SCENARIO_ID,'.RData'))
}

# function for scenario based worker

Worker_Function = function(core_nr){
  library(dacomp)
  library(ancom.R)
  library(phyloseq)
  library(metagenomeSeq)
  library(ALDEx2)
  library(Wrench)
  library(DESeq2)
  
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
        
        ref_select = dacomp::dacomp.select_references(X = data$X,
                                                      median_SD_threshold = ref_select_method_Median_SD_Threshold,
                                                      minimal_TA = 10,
                                                      maximal_TA = 200,
                                                      Pseudo_Count_used = PS_used,
                                                      verbose = DEBUG.SINGLE.CORE,
                                                      select_from = to_select_from)
        
        Selected_References = ref_select$selected_references
                                                      
      }
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
      
      bad_select = length(which(ref_select$selected_references %in% data$select_diff_abundant))
      
      construct_results_row = function(p.values.test,ref_select_method_Label){
        rejected = which(p.adjust(p.values.test,method = CORRECTION_TYPE_SUBZERO) < Q_LEVEL) #pvals
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
      
      results[current_row,] = construct_results_row(res_wilcoxon$p.values.test.ratio ,
                                                    paste0(PREFIX_DACOMP_WILCOXON,PREFIX_DACOMP_rarefaction,ref_select_method_Label)
                                                    ); current_row = current_row + 1
      
      results[current_row,] = construct_results_row(res_wilcoxon$p.values.test.ratio.normalization ,
                                                    paste0(PREFIX_DACOMP_WILCOXON,PREFIX_DACOMP_division,ref_select_method_Label)
                                                    ); current_row = current_row + 1
      
      results[current_row,] = construct_results_row(res_welch$p.values.test.ratio ,
                                                    paste0(PREFIX_DACOMP_WELCH,PREFIX_DACOMP_rarefaction,ref_select_method_Label)
                                                    ); current_row = current_row + 1
      
      results[current_row,] = construct_results_row(res_welch$p.values.test.ratio.normalization ,
                                                    paste0(PREFIX_DACOMP_WELCH,PREFIX_DACOMP_division,ref_select_method_Label)
                                                    ); current_row = current_row + 1
      
      # results[current_row,] = c(SCENARIO_ID,ref_select_method_Label,length(rejected),true_positive,
      #                           false_positive, bad_select) ; current_row = current_row + 1
      # 
      
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
    
    #ALDEx2 (note that we alter the data object here):
    indicies_not_empty = which(apply(data$X>0,2,sum)>0)
    m_original = ncol(data$X)
    flags_H1 = rep(0,m_original)
    flags_H1[data$select_diff_abundant] = 1
    flags_H1 = flags_H1[indicies_not_empty]
    data$X = data$X[,indicies_not_empty]
    data$select_diff_abundant = which(flags_H1==1)
    
    aldex.res.iqlr <- aldex(t(data$X), as.character(data$Y), mc.samples=128, denom="iqlr",
                            test="t", effect=FALSE,verbose = DEBUG.SINGLE.CORE)
    
    rejected.iqlr.wi = p.adjust(aldex.res.iqlr$wi.ep, method = CORRECTION_TYPE_SUBZERO) <= Q_LEVEL
    rejected.iqlr.we = p.adjust(aldex.res.iqlr$we.ep, method = CORRECTION_TYPE_SUBZERO) <= Q_LEVEL
    
    ALDEx2.Nr.Rejected.Wi = sum(rejected.iqlr.wi)
    ALDEx2.Nr.Rejected.We = sum(rejected.iqlr.we)
    ALDEx2.Nr.TP.Wi = sum(which(rejected.iqlr.wi) %in% data$select_diff_abundant)
    ALDEx2.Nr.TP.We = sum(which(rejected.iqlr.we) %in% data$select_diff_abundant)
    ALDEx2.Nr.FP.Wi = ALDEx2.Nr.Rejected.Wi - ALDEx2.Nr.TP.Wi
    ALDEx2.Nr.FP.We = ALDEx2.Nr.Rejected.We - ALDEx2.Nr.TP.We
    
    results[current_row,] = c(SCENARIO_ID,METHOD_LABEL_ALDEx2_Wilcoxon, ALDEx2.Nr.Rejected.Wi, ALDEx2.Nr.TP.Wi ,ALDEx2.Nr.FP.Wi,0) ; current_row = current_row + 1
    results[current_row,] = c(SCENARIO_ID,METHOD_LABEL_ALDEx2_Welch, ALDEx2.Nr.Rejected.We, ALDEx2.Nr.TP.We ,ALDEx2.Nr.FP.We,0) ; current_row = current_row + 1
    
    #Wrench - note that we removed otus with zeroes before
    W <- wrench( t(data$X), condition=data$Y  )
    compositionalFactors <- W$ccf
    normalizationFactors <- W$nf
    
    library(DESeq2)
    deseq.obj <- DESeq2::DESeqDataSetFromMatrix(countData = t(data$X),
                                                DataFrame(Y = factor(data$Y)),
                                                ~ Y )
    sizeFactors(deseq.obj) <- normalizationFactors
    deseq2_res = DESeq2::DESeq(deseq.obj)
    deseq2_res2 = results(deseq2_res)
    
    
    wrench_rejected = which(p.adjust(deseq2_res2$pvalue,method = CORRECTION_TYPE_SUBZERO) <= Q_LEVEL)
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

sink(file = REFSIM_warnings_file(RESULTS_DIR,SCENARIO_ID))
warnings()
sink()

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
