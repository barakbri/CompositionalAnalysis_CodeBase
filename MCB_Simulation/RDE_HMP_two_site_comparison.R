# Analyze HMP Data
# This file analyzes a single pair of body sites


#Take two sites
# G1: "Stool"
# G2:  "Saliva"                       "Tongue_dorsum"  
# "Hard_palate"                  "Buccal_mucosa"                "Attached_Keratinized_gingiva"
# "Palatine_Tonsils"             "Throat"                       "Supragingival_plaque"        
# "Subgingival_plaque"           

# G3: "Right_Antecubital_fossa" "Left_Retroauricular_crease"   "Right_Retroauricular_crease"  "Left_Antecubital_fossa" 
# "Anterior_nares"              
# G4:   "Vaginal_introitus" "Posterior_fornix" "Mid_vagina"           

           
#if you run from meta, you need to comment these
#Y0 = 'Supragingival_plaque' #'Vaginal_introitus'#  'Supragingival_plaque'
#Y1 = 'Palatine_Tonsils' #'Mid_vagina'#  'Palatine_Tonsils'
#LOAD_DATA = T

DISABLE_ANCOM = F

HMP_RESULTS_DIR = "../../Results/HMP/"
# This will eventually be saved to file:
results_to_save = list()

util_name_from_bodysite = function(X){
  return(strsplit(X,'UBERON:')[[1]][2])
}

####
#Section I: Parameter Declaration:
##########################################################

Q = 0.1

MINIMAL_RATIO_OTHER = 0.025
AGG_LEVEL = 6 #Genus
GLOBAL_RAREFACTION_QUANTILE = 0.1 # under this number - measurements will be discarded.
MIN_TA = 10
MAX_TA = 200
ADJUSTMENT = 0
MIN_PREV = 0.00
MED_SD_THRES = 1.3
SENSITIVITY_ANALYSIS_BY_BH = T

NORMAL_APPROX = T
NR_PERMUTATIONS = 30000
TEST_USED = 'Wilcoxon.MultipleX'#'Wilcoxon.MultipleX'  TwoPartWilcoxon.MultipleX
NR_MULTIPLE_RAREFACTION = 1
VERBOSE_MSGS = T
LAMBDA_ITERATION_MULTIPLIER = 0.1
LAMBDA_MULTIPLIER_MINIMAL = 0.1#0.01

results_to_save$parameters_list = list(Q = Q, MINIMAL_RATIO_OTHER = MINIMAL_RATIO_OTHER,
                                       AGG_LEVEL = AGG_LEVEL, GLOBAL_RAREFACTION_QUANTILE = GLOBAL_RAREFACTION_QUANTILE,
                                       MIN_TA = MIN_TA,MAX_TA = MAX_TA, MIN_PREV = MIN_PREV,
                                       ADJUSTMENT = ADJUSTMENT,MED_SD_THRES = MED_SD_THRES, SENSITIVITY_ANALYSIS_BY_BH = SENSITIVITY_ANALYSIS_BY_BH)
results_to_save$body_sites = list(Y0 = Y0,Y1 = Y1)
results_to_save$computation_parameters = list(NORMAL_APPROX = NORMAL_APPROX,
                                              NR_PERMUTATIONS = NR_PERMUTATIONS,
                                              TEST_USED = TEST_USED,
                                              NR_MULTIPLE_RAREFACTION = NR_MULTIPLE_RAREFACTION)


####
# Section II: Load data and parse it, Select body sites
##########################################################

if(VERBOSE_MSGS){
  cat(paste0('Loading data, body sites ',Y0,' and ',Y1,'\n\r'))
}

#source('C:/MCB2/MCB2/MCB_Simulation/Aggregate_by_Level.R')
#source('C:/MCB2/MCB2/MCB_Simulation/REFSIM_Selection_Methods.R')
#source(paste0('C:/MCB2/','MCB2/MCB_Simulation/SelectReferencesByRatios_and_Abundance.R')) #for plotting
#source('C:/MCB2/MCB2/MCB_Simulation/SelectReferences_MedianSD_Threshold.R')

source('ReferenceScore_Plot.R')

library(subzero)
library(ancom.R)
library(Matrix)
set.seed(1)

AGGREGATE_TO_GENUS = T
if(LOAD_DATA){
  counts_matrix = read.csv(file = '../../HMP2/OTU_Counts.csv')
  sample_data = read.csv(file = '../../HMP2/HMP_sample_data.csv')
  genus_type = as.character(counts_matrix[,ncol(counts_matrix)])
  counts_as_matrix_original = t(counts_matrix[,-ncol(counts_matrix)])
  counts_as_matrix_original = counts_as_matrix_original[-1,]
  sample_data$X.SampleID = paste0('X',sample_data$X.SampleID)
  counts_as_matrix_original = (Matrix::as.matrix(counts_as_matrix_original))  
  counts_as_matrix_original_cast = matrix(as.numeric(counts_as_matrix_original),ncol = ncol(counts_as_matrix_original),nrow=nrow(counts_as_matrix_original))
  rownames(counts_as_matrix_original_cast) = rownames(counts_as_matrix_original)
  counts_as_matrix_original = counts_as_matrix_original_cast
  
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



#as.character(unique(sample_data$HMPbodysubsite))

dt_ind_to_keep = which(sample_data$HMPbodysubsite %in% c(Y0,Y1))
dt = sample_data[dt_ind_to_keep,]
dt$HMPbodysubsite = factor(dt$HMPbodysubsite)
dt$HMPbodysubsite = as.numeric(dt$HMPbodysubsite) - 1

#filter the X_matrix as well:
X_ind_to_keep = which(rownames(counts_as_matrix_original) %in% dt$X.SampleID)
counts_as_matrix = counts_as_matrix_original[X_ind_to_keep,]
#function returns the permutation that will order to order by orer by
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
if(sum(is.na(ordering_perm)) != 0){
  stop('Ordering Perm invalid')
}

#head(rownames(X))
#head(as.character(dt$X.SampleID  [ordering_perm]) )
if(!all.equal(rownames(X),as.character(dt$X.SampleID[ordering_perm]) )){
  stop("ERROR ordering rows of X and Y")
}


Y = Y[ordering_perm]
Y = as.numeric(Y) - 1
subject_id = subject_id[ordering_perm]

selected_unique = !duplicated(subject_id)
if(unique(table(subject_id[selected_unique])) != 1)
  stop("ERROR in selecting unique from subjects")

X = X[selected_unique, ]
Y = Y[selected_unique]
subject_id = subject_id[selected_unique]

#X_cast = matrix(as.numeric(X),ncol = ncol(X),nrow=nrow(X))
#all.equal(X_cast[,2],as.numeric(X[,2])) #check
#X = X_cast

####
#Section III: Run pipeline over subjects
##########################################################


if(VERBOSE_MSGS){
  cat(paste0('Running pipeline\n\r'))
}
reads_per_subject = apply(X,1,sum)


#Aggregate to genus and compute other
#aggregated_data = aggregate_counts_by_level(counts_as_matrix = X,x_metadata = x_metadata,L = AGG_LEVEL,VERBOSE = F)
#no need to aggregatre here, data is already at genus level

ind_to_other = which(as.numeric(apply(1*(X >0),2,sum)) < MINIMAL_RATIO_OTHER*nrow(X))
other_counts = apply(X[,ind_to_other],1,sum)
X_2 = X[,-c(ind_to_other)]##cbind(X[,-c(ind_to_other)] , other_counts)
genus_type_aggregated = as.character(genus_type[-ind_to_other])#c(as.character(genus_type[-ind_to_other]),'Other')

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


# selected_references_obj_for_plot = select.references.by.pair.ratios.and.abundance(X = X_2,
#                                                                          target_abundance =  c(2,5,10,20,50),
#                                                                          adjustment = ADJUSTMENT,
#                                                                          MIN_PREV = MIN_PREV,verbose = F)


selected_references_obj = select.references.Median.SD.Threshold(X = X_2,median_SD_threshold = MED_SD_THRES,minimal_TA = MIN_TA,
                                                                maximal_TA = MAX_TA,select_from = 1:(ncol(X_2)),verbose = F)

ref_analysis_save_file = paste0(HMP_RESULTS_DIR,"Ref_selection_analysis_",(Y0),"_",Y1,'.pdf')
pdf(file = ref_analysis_save_file,width = 8,height = 5)
plot_ref_select_scores(selected_references_obj,label = "") #paste0('Reference Selection for methods: ',(Y0),"_",Y1)
dev.off()

results_to_save$selected_references_obj = selected_references_obj
#mean(selected_references_obj$mean_prevalence >= MIN_PREV)
selected_references = selected_references_obj$selected_references
if(!is.null(selected_references) & length(selected_references) <= dim(X_2)[2] - 2 ){
  ref_values = NULL
  if(length(selected_references)>1){
    ref_values = sort(as.numeric(apply(X_2[,selected_references],MARGIN = 1,FUN = sum )))
  }else{
    ref_values = sort(as.numeric(X_2[,selected_references]))
  }
  
  results_to_save$ref_values = ref_values
  results_to_save$ref_values_sorted = sort(ref_values)
  
  count_rejections_for_subzero = function(res_x,alpha = 0.05){
    return(length(which(p.adjust(res_x$pvalue,method = 'BH')<=alpha)))
  }
  
  
  set.seed(1)
  flag = T
  if(VERBOSE_MSGS)
    print(paste0('Computing DS-FDR'))
  
  res_dsfdr = subzero.dfdr.test(X = X_2,Y,as.numeric(selected_references),
                                nr_perm = ceiling(1/(Q/(dim(X_2)[2]))),verbose = F,q = Q)
  results_to_save$res_dsfdr = res_dsfdr
  
  rejected = res_dsfdr$rejected
  rejected_by_pval = which(p.adjust(res_dsfdr$p.values.test,method = 'BH')<=Q)
  
  results_to_save$rejected = rejected
  results_to_save$rejected_by_pval = rejected_by_pval
  
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
  
  # ANCOM_GRAPH_FILENAME = paste0(HMP_RESULTS_DIR,"",strsplit(Y0,'UBERON:')[[1]][2],"_",strsplit(Y1,'UBERON:')[[1]][2],"_ANCOM_GRAPH",'.png')
  # png(filename = ANCOM_GRAPH_FILENAME)
  # plot(ecdf(ANCOM_res$W),xlim = c(0,dim(X_ANCOM)[2]),ylim = c(0,1),main="ECDF of # of Wilcoxon Rejection")
  # dev.off()
  # 
  ANCOM_Detected = which(paste0('X',colnames(otu_dt)) %in% ANCOM_res$detected)  
}

#length(ANCOM_Detected)

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

discovered_subzero = c(int_obj$shared_with_ANCOM,int_obj$unique)

Prevalence_Y0 = as.numeric( apply(X_2[Y == 0, ], 2, function(x){mean(x>0)}) )
Prevalence_Y1 = as.numeric( apply(X_2[Y == 1, ], 2, function(x){mean(x>0)}) )

if(length(ANCOM_Detected)>0){
  Prevalence_ANCOM_Y0 = as.numeric( apply(X_2[Y == 0, ANCOM_Detected,drop = F], 2, function(x){mean(x>0)}) )
  Prevalence_ANCOM_Y1 = as.numeric( apply(X_2[Y == 1, ANCOM_Detected,drop = F], 2, function(x){mean(x>0)}) )  
}

if(length(discovered_subzero)>0){
  Prevalence_SUBZERO_Y0 = as.numeric( apply(X_2[Y == 0, discovered_subzero,drop = F], 2, function(x){mean(x>0)}) )
  Prevalence_SUBZERO_Y1 = as.numeric( apply(X_2[Y == 1, discovered_subzero,drop = F], 2, function(x){mean(x>0)}) )  
}



#sort(pmin(Prevalence_ANCOM_Y1,Prevalence_ANCOM_Y0))
#sort(pmin(Prevalence_Y0,Prevalence_Y1))


ecdf_prev_total_x = sort(pmin(Prevalence_Y0,Prevalence_Y1))
if(length(ANCOM_Detected)>0)
  ecdf_prev_ANCOM_x = sort(sort(pmin(Prevalence_ANCOM_Y1,Prevalence_ANCOM_Y0)))
if(length(discovered_subzero)>0)
  ecdf_prev_SUBZERO_x = sort(sort(pmin(Prevalence_SUBZERO_Y0,Prevalence_SUBZERO_Y1)))

ecdf_prev_total_x_max = sort(pmax(Prevalence_Y0,Prevalence_Y1))
if(length(ANCOM_Detected)>0)
  ecdf_prev_ANCOM_x_max = sort(sort(pmax(Prevalence_ANCOM_Y1,Prevalence_ANCOM_Y0)))
if(length(discovered_subzero)>0)
  ecdf_prev_SUBZERO_x_max = sort(sort(pmax(Prevalence_SUBZERO_Y0,Prevalence_SUBZERO_Y1)))

# ecdf_prev_total_x = sort(abs(Prevalence_Y0-Prevalence_Y1))
# ecdf_prev_ANCOM_x = sort(abs(Prevalence_ANCOM_Y1-Prevalence_ANCOM_Y0))
# ecdf_prev_SUBZERO_x = sort(abs(Prevalence_SUBZERO_Y0-Prevalence_SUBZERO_Y1))
# 
# 
ecdf_prev_total = ecdf(ecdf_prev_total_x)
if(length(ANCOM_Detected)>0)
  ecdf_prev_ANCOM = ecdf(ecdf_prev_ANCOM_x)
if(length(discovered_subzero)>0)
  ecdf_prev_SUBZERO = ecdf(ecdf_prev_SUBZERO_x)

ecdf_prev_total_max = ecdf(ecdf_prev_total_x_max)
if(length(ANCOM_Detected)>0)
  ecdf_prev_ANCOM_max = ecdf(ecdf_prev_ANCOM_x_max)
if(length(discovered_subzero)>0)
  ecdf_prev_SUBZERO_max = ecdf(ecdf_prev_SUBZERO_x_max)

ecdf_prev_total_y = ecdf_prev_total( ecdf_prev_total_x )
if(length(ANCOM_Detected)>0)
  ecdf_prev_ANCOM_y = ecdf_prev_ANCOM( ecdf_prev_ANCOM_x )
if(length(discovered_subzero)>0)
  ecdf_prev_SUBZERO_y = ecdf_prev_SUBZERO( ecdf_prev_SUBZERO_x )

ecdf_prev_total_y_max = ecdf_prev_total_max( ecdf_prev_total_x_max )
if(length(ANCOM_Detected)>0)
  ecdf_prev_ANCOM_y_max = ecdf_prev_ANCOM_max( ecdf_prev_ANCOM_x_max )
if(length(discovered_subzero)>0)
  ecdf_prev_SUBZERO_y_max = ecdf_prev_SUBZERO_max( ecdf_prev_SUBZERO_x_max )



png(paste0(HMP_RESULTS_DIR,'\\Rejections_Compare_',Y0,'_',Y1,'.png'),width = 800,height = 600)

par(mfrow = c(1,2))
plot(ecdf_prev_total_x, ecdf_prev_total_y,type = 'l', xlab = 'ECDF of Minor prevalence of taxa, by subgroups',ylab = 'F')
if(length(ANCOM_Detected)>0)
  lines(ecdf_prev_ANCOM_x, length(ANCOM_Detected)/ncol(X_2)*ecdf_prev_ANCOM_y,type = 'l',col = 'red')
if(length(discovered_subzero)>0)
  lines(ecdf_prev_SUBZERO_x, length(discovered_subzero)/ncol(X_2)* ecdf_prev_SUBZERO_y,type = 'l',col = 'blue')
legend(0.5, 0.5, legend=c("All Taxa", "ANCOM discoveries","SubZero discoveries"),
       col=c("black","red", "blue"), lty=1, cex=0.8)

plot(ecdf_prev_total_x_max, ecdf_prev_total_y_max,type = 'l', xlab = 'ECDF of Major prevalence of taxa, by subgroups',ylab = 'F',xlim = c(0,1))
if(length(ANCOM_Detected)>0)
  lines(ecdf_prev_ANCOM_x_max, length(ANCOM_Detected)/ncol(X_2)*ecdf_prev_ANCOM_y_max,type = 'l',col = 'red')
if(length(discovered_subzero)>0)
  lines(ecdf_prev_SUBZERO_x_max, length(discovered_subzero)/ncol(X_2)*ecdf_prev_SUBZERO_y_max,type = 'l',col = 'blue')
legend(0.0, 1.0, legend=c("All Taxa", "ANCOM discoveries","SubZero discoveries"),
       col=c("black","red", "blue"), lty=1, cex=0.8)

par(mfrow = c(1,1))
dev.off()


pvalues_for_test = res_dsfdr$p.values
pvalues_for_test[selected_references] = NA
rej = which(p.adjust(pvalues_for_test,method = 'BH')<=Q)
non_rej = which(p.adjust(pvalues_for_test,method = 'BH')>Q)

shared_rejections_ind = rej[which(rej %in% ANCOM_Detected)]
non_shared_rejections_ind = non_rej[which(non_rej %in% ANCOM_Detected)]
ANCOM_rejections_in_reference = ANCOM_Detected[which(ANCOM_Detected %in% selected_references)]

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

#compare to subset

subset_ratio = 2/3
subset_size = ceiling(nrow(X_2)*subset_ratio)
selected_for_subset = sample(1:nrow(X_2),size = subset_size)
X_2_subset = X_2[selected_for_subset,]
Y_subset = Y[selected_for_subset]

selected_references_obj_subset = select.references.Median.SD.Threshold(X = X_2_subset,median_SD_threshold = MED_SD_THRES,minimal_TA = MIN_TA,maximal_TA = MAX_TA,select_from = 1:(ncol(X_2_subset)-1))

res_dsfdr_subset = subzero.dfdr.test(X = X_2_subset,Y_subset,as.numeric(selected_references_obj_subset$selected_references),q = Q,
                              nr_perm = ceiling(1/(Q/(dim(X_2)[2]))),verbose = F)
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
#Section VI: Save results
##########################################################


results_save_file = paste0(HMP_RESULTS_DIR,"RESULTS_FILE_",(Y0),"_",Y1,'.Rdata')
save(results_to_save,file = results_save_file)


if(VERBOSE_MSGS){
  cat(paste0('Done body sites ',Y0,' and ',Y1,'\n\r'))
  # pvals = res_dsfdr$p.values
  # pvals_test = pvals; pvals_test[selected_references] = NA
  # pvals_ref = pvals; pvals_ref[-selected_references] = NA
  # rej_lambda1 = length(results_to_save$res_dsfdr_lambda_1_0_rejected)
  # nr_rej =length(which(p.adjust(pvals_test,method = 'BH')<=Q))
  # nr_rej_ref =length(which(p.adjust(pvals_ref,method = 'BH')<=Q))
  # cat(paste0('Rej: ',rej_lambda1 ,', Rej_sens:',nr_rej,' , Rej-ref',nr_rej_ref,'\n\r'))
}

