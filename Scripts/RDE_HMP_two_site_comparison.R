# This script is used for analyzing differential abundance between two body sites, in the HMP data example.
# This script is called repeatedly from RDE_HMP_main.R
# The parameters for the analyzed body sites and if to load the data from disk are commented out,
# since they are manager from the batch script. To run this script for a single comparison of body sites,
# uncomment Y0, Y1 and LOAD_DATA


# List of body sites, by groups:
# G1: "Stool"

# G2:  "Saliva"                       "Tongue_dorsum"  
# "Hard_palate"                  "Buccal_mucosa"                "Attached_Keratinized_gingiva"
# "Palatine_Tonsils"             "Throat"                       "Supragingival_plaque"        
# "Subgingival_plaque"           

# G3: "Right_Antecubital_fossa" "Left_Retroauricular_crease"   "Right_Retroauricular_crease"
#"Left_Antecubital_fossa" "Anterior_nares"              

# G4:   "Vaginal_introitus" "Posterior_fornix" "Mid_vagina"           

       
#Y0 = 'Supragingival_plaque' #'Vaginal_introitus'#  'Supragingival_plaque'
#Y1 = 'Palatine_Tonsils' #'Mid_vagina'#  'Palatine_Tonsils'
#LOAD_DATA = T

DISABLE_ANCOM = F

HMP_RESULTS_DIR = "../../Results/HMP/"
if(!dir.exists(HMP_RESULTS_DIR)){
  dir.create(HMP_RESULTS_DIR)
}
# This will eventually be saved to file:
results_to_save = list()

#used to extract body site name from the data file
util_name_from_bodysite = function(X){
  return(strsplit(X,'UBERON:')[[1]][2])
}

####
#Section I: Parameter Declaration:
##########################################################

Q = 0.1 #for BH

MINIMAL_RATIO_OTHER = 0.025  #under this level of prevalence - taxa will be removed from the study
AGG_LEVEL = 6 #Genus
GLOBAL_RAREFACTION_QUANTILE = 0.1 # under this number - measurements will be discarded.
# parameters for reference selection:
MIN_TA = 10
MAX_TA = 200
MED_SD_THRES = 1.3

NR_PERMUTATIONS = 30000
VERBOSE_MSGS = T

results_to_save$parameters_list = list(Q = Q, MINIMAL_RATIO_OTHER = MINIMAL_RATIO_OTHER,
                                       AGG_LEVEL = AGG_LEVEL, GLOBAL_RAREFACTION_QUANTILE = GLOBAL_RAREFACTION_QUANTILE,
                                       MIN_TA = MIN_TA,MAX_TA = MAX_TA,
                                       MED_SD_THRES = MED_SD_THRES)
results_to_save$body_sites = list(Y0 = Y0,Y1 = Y1)
results_to_save$computation_parameters = list(NR_PERMUTATIONS = NR_PERMUTATIONS)


####
# Section II: Load data and parse it, Select body sites
##########################################################

if(VERBOSE_MSGS){
  cat(paste0('Loading data, body sites ',Y0,' and ',Y1,'\n\r'))
}

source('ReferenceScore_Plot.R')

library(dacomp)
library(ancom.R)
library(Matrix)
set.seed(1)

AGGREGATE_TO_GENUS = T
if(LOAD_DATA){ # load data - this is called only once if using the batch script
  counts_matrix = read.csv(file = '../HMP/OTU_Counts.csv')
  sample_data = read.csv(file = '../HMP/HMP_sample_data.csv')
  genus_type = as.character(counts_matrix[,ncol(counts_matrix)]) # last col contains genus level data
  counts_as_matrix_original = t(counts_matrix[,-ncol(counts_matrix)])
  counts_as_matrix_original = counts_as_matrix_original[-1,]
  sample_data$X.SampleID = paste0('X',sample_data$X.SampleID) # moving names to the same format across the two files
  counts_as_matrix_original = (Matrix::as.matrix(counts_as_matrix_original))  
  counts_as_matrix_original_cast = matrix(as.numeric(counts_as_matrix_original),ncol = ncol(counts_as_matrix_original),nrow=nrow(counts_as_matrix_original))
  rownames(counts_as_matrix_original_cast) = rownames(counts_as_matrix_original)
  counts_as_matrix_original = counts_as_matrix_original_cast
  
  # aggregate to genus level if needed
  if(AGGREGATE_TO_GENUS){
    unique_genera = unique(genus_type)
    list_OTUS_for_genera = list()
    counts_as_matrix_aggregated = matrix(nrow = nrow(counts_as_matrix_original),ncol = length(unique_genera))
    for(i in 1:length(unique_genera)){
     list_OTUS_for_genera[[i]] = which(genus_type == unique_genera[i])
     if(length(list_OTUS_for_genera[[i]])>1){
       counts_as_matrix_aggregated[,i] = apply(counts_as_matrix_original[, list_OTUS_for_genera[[i]]],1,sum)
     }else{
       counts_as_matrix_aggregated[,i] = counts_as_matrix_original[, list_OTUS_for_genera[[i]]]
     }
    }
    genus_type = unique_genera
    row.names(counts_as_matrix_aggregated) = row.names(counts_as_matrix_original)
    counts_as_matrix_original = counts_as_matrix_aggregated
  }
  
}
#from this point on we do not edit counts_as_matrix_original and sample_data - the will be used by different comparison

#keep only body sites in the current pair.
dt_ind_to_keep = which(sample_data$HMPbodysubsite %in% c(Y0,Y1))
dt = sample_data[dt_ind_to_keep,]
dt$HMPbodysubsite = factor(dt$HMPbodysubsite)
dt$HMPbodysubsite = as.numeric(dt$HMPbodysubsite) - 1

#filter the X_matrix as well:
X_ind_to_keep = which(rownames(counts_as_matrix_original) %in% dt$X.SampleID)
counts_as_matrix = counts_as_matrix_original[X_ind_to_keep,]
#function returns the permutation that will order to_order by order_by
order_one_by_another = function(to_order, order_by){
  ret = rep(NA,length(order_by))
  for(i in 1:length(order_by)){
    ret[i] = which(to_order == order_by[i])
  }
  return(ret)
}

X = counts_as_matrix
Y = as.factor(dt$HMPbodysubsite)
subject_id = dt$RSID

ind_to_keep = which(dt$X.SampleID %in% rownames(X))
ordering_perm = order_one_by_another(as.character(dt$X.SampleID), as.character(rownames(X)))

#check ordering is valid
if(sum(is.na(ordering_perm)) != 0){
  stop('Ordering Perm invalid')
}

if(!all.equal(rownames(X),as.character(dt$X.SampleID[ordering_perm]) )){
  stop("ERROR ordering rows of X and Y")
}


Y = Y[ordering_perm]
Y = as.numeric(Y) - 1 # convert to 0,1
subject_id = subject_id[ordering_perm]

#select only one sample, from each subject (one of the body sites)
selected_unique = !duplicated(subject_id)
if(unique(table(subject_id[selected_unique])) != 1)
  stop("ERROR in selecting unique from subjects")

X = X[selected_unique, ]
Y = Y[selected_unique]
subject_id = subject_id[selected_unique]

####
#Section III: Run pipeline over subjects
##########################################################


if(VERBOSE_MSGS){
  cat(paste0('Running pipeline\n\r'))
}
reads_per_subject = apply(X,1,sum)


#remove all rare taxa ("other")

ind_to_other = which(as.numeric(apply(1*(X >0),2,sum)) < MINIMAL_RATIO_OTHER*nrow(X))
other_counts = apply(X[,ind_to_other],1,sum) #sometimes usefull to look at this,to see how much load there is at all the small taxa
X_2 = X[,-c(ind_to_other)]
genus_type_aggregated = as.character(genus_type[-ind_to_other])

#Remove observations with a low number of counts:
to_keep = which( reads_per_subject >= quantile(reads_per_subject,probs = GLOBAL_RAREFACTION_QUANTILE))
X2_Original = X_2
Y_Original = Y

X_2 = X_2[to_keep,]
Y = Y[to_keep]

results_to_save[["dim(X)"]] = dim(X)
results_to_save[["dim(X_2)"]] = dim(X_2)

#document taxa sums
taxa_sums = apply(X_2,2,sum)
results_to_save$taxa_sums = taxa_sums
results_to_save$empty_taxa = which(taxa_sums == 0)

#select references
selected_references_obj = dacomp::dacomp.select_references(X = X_2,median_SD_threshold = MED_SD_THRES,minimal_TA = MIN_TA,maximal_TA = MAX_TA,select_from = 1:(ncol(X_2)),verbose = F,Pseudo_Count_used  = 1)

#plot histogram for selected references
ref_analysis_save_file = paste0(HMP_RESULTS_DIR,"Ref_selection_analysis_",(Y0),"_",Y1,'.pdf')
pdf(file = ref_analysis_save_file,width = 8,height = 5)
plot_ref_select_scores(selected_references_obj,label = "") #paste0('Reference Selection for methods: ',(Y0),"_",Y1)
dev.off()

results_to_save$selected_references_obj = selected_references_obj

selected_references = selected_references_obj$selected_references
if(!is.null(selected_references) & length(selected_references) <= dim(X_2)[2] - 2 ){ #if we have a valid set of refernces, perform analysis:
  ref_values = NULL
  if(length(selected_references)>1){
    ref_values = sort(as.numeric(apply(X_2[,selected_references],MARGIN = 1,FUN = sum )))
  }else{
    ref_values = sort(as.numeric(X_2[,selected_references]))
  }
  
  results_to_save$ref_values = ref_values
  results_to_save$ref_values_sorted = sort(ref_values)
  
  count_rejections_for_DACOMP = function(res_x,alpha = 0.05){
    return(length(which(p.adjust(res_x$pvalue,method = 'BH')<=alpha)))
  }
  
  
  set.seed(1)
  flag = T
  if(VERBOSE_MSGS)
    print(paste0('Computing dacomp'))
  #res_dsfdr is a misleading name, because we will use P-values eventully, so that different method are comperable
  res_dsfdr = dacomp::dacomp.test(X = X_2, y = Y, test = DACOMP.TEST.NAME.WILCOXON, ind_reference_taxa = as.numeric(selected_references), nr_perm = 30000, q = Q, compute_ratio_normalization = T)
  res_dsfdr_Welch = dacomp::dacomp.test(X = X_2, y = Y, test = DACOMP.TEST.NAME.WELCH_LOGSCALE, ind_reference_taxa = as.numeric(selected_references), nr_perm = 30000, q = Q, compute_ratio_normalization = T)
  
  
  results_to_save$res_dsfdr = res_dsfdr
  results_to_save$res_dsfdr_Welch = res_dsfdr_Welch
  
  rejected = res_dsfdr$rejected
  rejected_by_pval = which(p.adjust(res_dsfdr$p.values.test,method = 'BH') <= Q)
  rejected_by_pval_division = which(p.adjust(res_dsfdr$p.values.test.ratio.normalization,method = 'BH') <= Q)
  
  rejected_by_pval_Welch = which(p.adjust(res_dsfdr$p.values.test,method = 'BH') <= Q)
  rejected_by_pval_division_Welch = which(p.adjust(res_dsfdr$p.values.test.ratio.normalization,method = 'BH') <= Q)
  
  results_to_save$rejected = rejected
  results_to_save$rejected_by_pval = rejected_by_pval
  results_to_save$rejected_by_pval_division = rejected_by_pval_division
  results_to_save$rejected_by_pval_Welch = rejected_by_pval_Welch
  results_to_save$rejected_by_pval_division_Welch = rejected_by_pval_division_Welch

}


####
#Section IV: Compare to ANCOM
##########################################################

if(VERBOSE_MSGS){
  cat(paste0('Running ANCOM \n\r'))
}

X_ANCOM = X2_Original
library(ancom.R)
otu_dt = as.matrix(cbind(X_ANCOM, as.numeric(Y_Original) - 1 ))
colnames(otu_dt) = as.character(1:(ncol(otu_dt)))
otu_dt = as.data.frame(otu_dt)

ANCOM_res = NULL
ANCOM_Detected = integer(0)
if(!DISABLE_ANCOM){
  ANCOM_res = ANCOM(otu_dt, sig = Q, multcorr = 3)
  
  #plots ecdf for ancom W statistics
  ANCOM_GRAPH_FILENAME = paste0(HMP_RESULTS_DIR,"",strsplit(Y0,'UBERON:')[[1]][2],"_",strsplit(Y1,'UBERON:')[[1]][2],"_ANCOM_GRAPH",'.png')
  png(filename = ANCOM_GRAPH_FILENAME)
  plot(ecdf(ANCOM_res$W),xlim = c(0,dim(X_ANCOM)[2]),ylim = c(0,1),main="ECDF of # of Wilcoxon Rejection")
  dev.off()

  ANCOM_Detected = which(paste0('X',colnames(otu_dt)) %in% ANCOM_res$detected)  
}

#auxiliray function for computing intersection witn ANCOM
intersect_with_ANCOM = function(res_dsfdr,ANCOM_Detected_Original,alpha = Q){
  ret = list()
  pvalues_for_test = res_dsfdr$p.values
  pvalues_for_test[selected_references] = NA
  rejected = which(p.adjust(res_dsfdr$p.values.test,method = 'BH')<=alpha)
  ret$shared_with_ANCOM = intersect(rejected,ANCOM_Detected_Original)
  ret$unique = setdiff(rejected,ret$shared_with_ANCOM)
  ret$unique_ANCOM = setdiff(ANCOM_Detected_Original,ret$shared_with_ANCOM)
  return(ret)
}


int_obj = intersect_with_ANCOM(res_dsfdr,ANCOM_Detected)
results_to_save$intersect_with_ANCOM = int_obj
results_to_save$ANCOM_res = ANCOM_res


pvalues_for_test = res_dsfdr$p.values
pvalues_for_test[selected_references] = NA
rej = which(p.adjust(pvalues_for_test,method = 'BH')<=Q)
non_rej = which(p.adjust(pvalues_for_test,method = 'BH')>Q)

#the length of these vectors are printed out in the output, but not displayed in the paper. can be used to discuss agreement between the different methods
shared_rejections_ind = rej[which(rej %in% ANCOM_Detected)] #between ANCOM and DACOMP
non_shared_rejections_ind = non_rej[which(non_rej %in% ANCOM_Detected)] # detected by ANCOM alone
ANCOM_rejections_in_reference = ANCOM_Detected[which(ANCOM_Detected %in% selected_references)] #detected by ANCOM, but was in our reference set...

#investigate abundance and prevalnce by shared or unique discoveries
abundance_shared_discoveries = median(as.numeric(apply(X_2[,shared_rejections_ind,drop = F],2,mean)))
abundance_ANCOM_UNIQUE_discoveries = median(as.numeric(apply(X_2[,non_shared_rejections_ind,drop = F],2,mean)))
abundance_ANCOM_in_ref = median((as.numeric(apply(X_2[,ANCOM_rejections_in_reference,drop = F],2,mean))))
prevalence_shared_discoveries = median(as.numeric(apply(X_2[,shared_rejections_ind,drop = F],2,function(x){mean(x>0)})))
prevalence_ANCOM_UNIQUE_discoveries = median(as.numeric(apply(X_2[,non_shared_rejections_ind,drop = F],2,function(x){mean(x>0)})))
prevalence_ANCOM_in_ref = median((as.numeric(apply(X_2[,ANCOM_rejections_in_reference,drop = F],2,function(x){mean(x>0)}))))


results_to_save$abundance_shared_discoveries = abundance_shared_discoveries
results_to_save$abundance_ANCOM_UNIQUE_discoveries = abundance_ANCOM_UNIQUE_discoveries
results_to_save$abundance_ANCOM_in_ref = abundance_ANCOM_in_ref
results_to_save$prevalence_shared_discoveries = prevalence_shared_discoveries
results_to_save$prevalence_ANCOM_UNIQUE_discoveries = prevalence_ANCOM_UNIQUE_discoveries
results_to_save$prevalence_ANCOM_in_ref = prevalence_ANCOM_in_ref

#compare to subset - run the method on 2/3 of the observations - check how power is affected

subset_ratio = 2/3
subset_size = ceiling(nrow(X_2)*subset_ratio)
selected_for_subset = sample(1:nrow(X_2),size = subset_size)
X_2_subset = X_2[selected_for_subset,]
Y_subset = Y[selected_for_subset]

selected_references_obj_subset = dacomp.select_references(X = X_2_subset,median_SD_threshold = MED_SD_THRES,minimal_TA = MIN_TA,maximal_TA = MAX_TA,select_from = 1:(ncol(X_2_subset)))


res_dsfdr_subset = dacomp.test(X = X_2_subset,
                               y = Y_subset,
                               ind_reference_taxa = as.numeric(selected_references_obj_subset$selected_references),
                               test = DACOMP.TEST.NAME.WILCOXON,
                               q = Q,
                              nr_perm = max(ceiling(1/(Q/(dim(X_2)[2]))),1000),
                              verbose = F)

pvalues_subset = res_dsfdr_subset$p.values
pvalues_subset[selected_references_obj_subset$selected_references] = NA
nr_rejections_subset = length(which(p.adjust(pvalues_subset,method = 'BH')<=Q))
results_to_save$res_dsfdr_subset = res_dsfdr_subset
results_to_save$nr_rejections_subset = nr_rejections_subset

####
#Section V: Wilcoxon
##########################################################

if(VERBOSE_MSGS){
  cat(paste0('Running Wilcoxon (naive)...\n\r'))
}

source( 'Wilcoxon_TaxaWise.R' )
Wilcoxon_res           = wilcoxon_taxa_wise(X_2,Y)
Wilcoxon_res_normalize = wilcoxon_taxa_wise(X_2,Y,normalize = T)
X_CSS = t(metagenomeSeq::cumNormMat(t(X_2)))
Wilcoxon_res_normalize_CSS = wilcoxon_taxa_wise(X_CSS,Y,normalize = F,normalize.P = 1)

results_to_save$Wilcoxon_res = Wilcoxon_res
results_to_save$Wilcoxon_res_normalize = Wilcoxon_res_normalize
results_to_save$Wilcoxon_res_normalize_CSS = Wilcoxon_res_normalize_CSS

results_to_save$Nr_Wilcoxon_Rejections = length( which(p.adjust(Wilcoxon_res$p.values,method = 'BH') <= Q) )
results_to_save$Nr_Wilcoxon_Normalize_Rejections = length( which(p.adjust(Wilcoxon_res_normalize$p.values,method = 'BH') <= Q) )
results_to_save$Nr_Wilcoxon_Normalize_Rejections_CSS = length( which(p.adjust(Wilcoxon_res_normalize_CSS$p.values,method = 'BH') <= Q) )


####
#Section VI: ALDEx2 and Wrench
##########################################################
library(ALDEx2)

#ALDEx2 
X_2_copy = X_2
indicies_not_empty = which(apply(X_2>0,2,sum)>0)
m_original = ncol(X_2)
original_ind = 1:m_original

X_2_copy = X_2_copy[,indicies_not_empty]
original_ind = original_ind[indicies_not_empty]

aldex.res.iqlr <- aldex(t(X_2_copy), as.character(Y), mc.samples=128, denom="iqlr",
                        test="t", effect=FALSE,verbose = T)

rejected.iqlr.wi = original_ind[which(p.adjust(aldex.res.iqlr$wi.ep, method = 'BH') <= Q)]
rejected.iqlr.we = original_ind[which(p.adjust(aldex.res.iqlr$we.ep, method = 'BH') <= Q)]

results_to_save$rejected.iqlr.wi = rejected.iqlr.wi
results_to_save$rejected.iqlr.we = rejected.iqlr.we

library(Wrench)
#Wrench - note that we removed otus with zeroes before
W <- wrench( t(X_2_copy), condition=Y )
compositionalFactors <- W$ccf
normalizationFactors <- W$nf

library(DESeq2)
deseq.obj <- DESeq2::DESeqDataSetFromMatrix(countData = t(X_2_copy),
                                            DataFrame(Y = factor(Y)),
                                            ~ Y )
sizeFactors(deseq.obj) <- normalizationFactors
deseq2_res = DESeq2::DESeq(deseq.obj)
deseq2_res2 = results(deseq2_res)


wrench_rejected = which(p.adjust(deseq2_res2$pvalue,method = 'BH') <= Q)
wrench_rejected = original_ind[wrench_rejected]

results_to_save$wrench_rejected = wrench_rejected

####
#Section VII: Save results
##########################################################


results_save_file = paste0(HMP_RESULTS_DIR,"RESULTS_FILE_",(Y0),"_",Y1,'.Rdata')
save(results_to_save,file = results_save_file)


if(VERBOSE_MSGS){
  cat(paste0('Done body sites ',Y0,' and ',Y1,'\n\r'))
}

