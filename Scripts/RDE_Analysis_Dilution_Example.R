#This script analyzes the data from the experiment of Stammler et al. (2016)
# "Adjusting  microbiome  profiles  for  differences  in  microbial  load  by spike-in bacteria."

#Required FDR level
Q_LVL = 0.1

#%%%%%%%%%%%%%%%%%
# Step I: load taxonomy data, reads, microbial loads, and experiment design matrix
#%%%%%%%%%%%%%%%%%

dt_meta = read.csv(file = '../Spike_in_Experiment/40168_2016_175_MOESM7_ESM.txt',sep = '\t',stringsAsFactors = F)
dt_loads = read.csv(file = '../Spike_in_Experiment/Loads.csv',stringsAsFactors = F)
dt_reads = read.csv(file = '../Spike_in_Experiment/40168_2016_175_MOESM9_ESM.txt',sep = '\t',stringsAsFactors = F)
dt_design = ((read.csv('../Spike_in_Experiment/Spike_in_Design.csv',sep = '\t')))
dt_design_sample = as.numeric(substr(colnames(dt_design)[-1],start = 8,nchar(colnames(dt_design)[-1])))
dt_design_S1 = (dt_design)[1,-1] #These are the designed loads
dt_design_S2 = (dt_design)[2,-1]
dt_design_S3 = (dt_design)[3,-1]

#These are the diluation factors, by sample
data_dilution_sample_id = c(65,66,67,68,69,70,
                       71,72,73,74,75,76,
                       77,78,79,80,81,82,
                       83,84,85,86,87,88,
                       89,90,91,92,93,94,
                       95,96,97,98,99,100)


data_dilution_factor = 1/rep(c(1,2.15,3.75,6.53,11.37,19.82),6)

#%%%%%%%%%%%%%%%%%
# Step II: order samples, from the different tables. Take the measured load as the trait variable:
#%%%%%%%%%%%%%%%%%

#length(dt_meta$Description) #these are the given samples
#length(dt_loads$Sample)

reads_OTU_ID = dt_reads[,1]
reads_counts = t(dt_reads[,2:37])

reads_taxonomy = dt_reads[,38]
row.names(reads_counts)

#which(row.names(reads_counts) %in%dt_meta$SampleID) #check that we have all the rows.

ind_in_loads = rep(NA,nrow(reads_counts))
index_in_other_system = rep(NA,nrow(reads_counts))
row_in_other_system = rep(NA,nrow(reads_counts))
load = rep(NA,nrow(reads_counts))
S1 = rep(NA,nrow(reads_counts))  #These will hold the absolute abundance (up to some coefficient) of the spikein, by spike in type
S2 = rep(NA,nrow(reads_counts)) 
S3 = rep(NA,nrow(reads_counts)) 
dilution_factor = rep(NA,nrow(reads_counts)) #This will hold the dilution factor, for the original fecal sample, by sample

#Note that two of the 38 samples are used as controls, for the spikein itself, and for an original sample with no spike-in, and no dilution.
#Therefore we have 36 samples in one table, and 38 in another.

for(i in 1:nrow(reads_counts)){
  ind_in_loads[i] = which(dt_meta$SampleID == row.names(reads_counts)[i])  
  index_in_other_system[i] = dt_meta$Description[ind_in_loads[i]]
  ind_temp = which(dt_design_sample == index_in_other_system[i])
   
  if(length(ind_temp)==0)
    next
  row_in_other_system[i] = ind_temp
  
  # A second option, not taken would have been to take the design matrix as the trait, rather than the actual spike-in measured
  # (The intended amount to spike-in, rather than the actual amount spiked in).
  # S1[i] = dt_design_S1[row_in_other_system[i]]
  # S2[i] = dt_design_S2[row_in_other_system[i]]
  # S3[i] = dt_design_S3[row_in_other_system[i]]

  row_in_other_system[i] = which(dt_loads$Sample == index_in_other_system[i])
  load[i] = dt_loads$Total.16S.rDNA.copies[row_in_other_system[i]]
  S1[i] = dt_loads$A..acidiphilus[row_in_other_system[i]]
  S2[i] = dt_loads$S..ruber[row_in_other_system[i]]
  S3[i] = dt_loads$R..radiobacter[row_in_other_system[i]]
  
  row_in_dilution = which(data_dilution_sample_id == index_in_other_system[i])
  dilution_factor[i] = data_dilution_factor[row_in_dilution]
  
}
S1 = as.numeric(S1)
S2 = as.numeric(S2)
S3 = as.numeric(S3)

#%%%%%%%%%%%%%%%%%
# Step III: remove obs with missing 16S
#%%%%%%%%%%%%%%%%%

ind_to_remove = which(is.na(S1) | is.na(S2) | is.na(S3)) #note that two of the samples did not have 16S reads.
#Therefore, they are removed from the data, despite haveing trait (spikein level) values.

S1 = S1[-ind_to_remove]
S2 = S2[-ind_to_remove]
S3 = S3[-ind_to_remove]
reads_counts = reads_counts[-ind_to_remove,]
load = load[-ind_to_remove]
dilution_factor = dilution_factor[-ind_to_remove]

#%%%%%%%%%%%%%%%%%
# Step IV: remove taxa which are all zeroes
#%%%%%%%%%%%%%%%%%

ind_to_keep = which(apply(reads_counts,2,sum)>0) 
reads_counts = reads_counts[,ind_to_keep]
reads_taxonomy = reads_taxonomy[ind_to_keep]
reads_OTU_ID = reads_OTU_ID[ind_to_keep]

#%%%%%%%%%%%%%%%%%
# Step V: locate the spikeins
#%%%%%%%%%%%%%%%%%

spikeins_id = c(grep('acidiphilus',reads_taxonomy),grep('Rhizobium',reads_taxonomy),grep('ruber',reads_taxonomy))

reads_taxonomy[spikeins_id]# these are the spikein

tot_non_spike_in = apply(reads_counts[,-spikeins_id],1,sum) # totals of non spikein

#remove samples where the spikein has taken over the sample.
#Note: actually we get to keep all samples...
taxa_where_spike_in_is_not_all_sample = which(tot_non_spike_in>100) 

S1 = S1[taxa_where_spike_in_is_not_all_sample]
S2 = S2[taxa_where_spike_in_is_not_all_sample]
S3 = S3[taxa_where_spike_in_is_not_all_sample]
reads_counts = reads_counts[taxa_where_spike_in_is_not_all_sample,]
load = load[taxa_where_spike_in_is_not_all_sample]
dilution_factor = dilution_factor[taxa_where_spike_in_is_not_all_sample]


#%%%%%%%%%%%%%%%%%
# Step VI: Preparations for testing: load packages, find references, 
# create TSS,CSS and CLR normalizations of the data.
#%%%%%%%%%%%%%%%%%


# create a matrix of all outcomes, so we iterate over them, one by one
Outcomes = list(S1 = S1, S3 = S3,dilution_factor = dilution_factor)


library(dacomp) #load packages
library(ALDEx2)

#find references for dacomp
ref_obj = dacomp::dacomp.select_references(X = reads_counts,
                                           median_SD_threshold = 0.0001,
                                           minimal_TA = 100,
                                           maximal_TA = 200,
                                           verbose = T)
#length(ref_obj$selected_references) = 951
#ref_obj$selected_MinAbundance = 107
cleaned_reference_by_3 = dacomp.validate_references(X = reads_counts,
                                               Y = Outcomes[[3]],
                                               ref_obj = ref_obj,
                                               test = DACOMP.TEST.NAME.SPEARMAN,
                                               Q_validation = 0.1,
                                               Minimal_Counts_in_ref_threshold = 10,
                                               Reduction_Factor = 0.9,
                                               Verbose = T,
                                               disable_DSFDR = T)

cleaned_reference_by_2 = dacomp.validate_references(X = reads_counts,
                                                    Y = Outcomes[[2]],
                                                    ref_obj = ref_obj,
                                                    test = DACOMP.TEST.NAME.SPEARMAN,
                                                    Q_validation = 0.1,
                                                    Minimal_Counts_in_ref_threshold = 10,
                                                    Reduction_Factor = 0.9,
                                                    Verbose = T,
                                                    disable_DSFDR = T)

cleaned_reference_by_1 = dacomp.validate_references(X = reads_counts,
                                                    Y = Outcomes[[1]],
                                                    ref_obj = ref_obj,
                                                    test = DACOMP.TEST.NAME.SPEARMAN,
                                                    Q_validation = 0.1,
                                                    Minimal_Counts_in_ref_threshold = 10,
                                                    Reduction_Factor = 0.9,
                                                    Verbose = T,
                                                    disable_DSFDR = T)

dacomp::dacomp.plot_reference_scores(ref_obj)
ref_obj$selected_MinAbundance

if(!all.equal(
length(ref_obj$selected_references),
length(cleaned_reference_by_3),
length(cleaned_reference_by_2),
length(cleaned_reference_by_1)
)){
  stop('Stopping: validation procedure found discrepencies between cleanned references')
}

#create matrix for CSS normalized data
X_CSS = t(metagenomeSeq::cumNormMat(t(reads_counts)))

#create matrix for TSS normalized data
X_TSS = reads_counts
for(i in 1:nrow(X_TSS)){
  X_TSS[i,] = X_TSS[i,] / sum(X_TSS[i,] )
}

# A wrapper function, for a marginal spearman test, to be used with CSS, TSS normalizations.
compute_spearman_marginal = function(X,Y,nr_perm = 10^4){
  Y_matrix = matrix(NA, ncol = nr_perm+1, nrow = nrow(X))
  Y_matrix[,1] = rank(Y,ties.method = 'average')
  Y_matrix[,1] = Y_matrix[,1] - mean(Y_matrix[,1])
  for( i in 2:ncol(Y_matrix)){
    Y_matrix[,i] = sample(Y_matrix[,1])
  }
  p.value.vec = rep(NA,ncol(X))
  for(i in 1:ncol(X)){
    ranked_X = rank(X[,i], ties.method = 'average') #break ties by average rank
    
    ranked_X = ranked_X - mean(ranked_X)
    current_stats = (dacomp:::rcpp_Spearman_PermTest_Given_Permutations(ranked_X,Y_matrix)[[1]]) #compute statistic and a sample of test statistic given from the null hypothesis
    current_stats = current_stats^2
    p.value.vec[i] = mean(current_stats>=current_stats[1])
    
    #Note that the asymp. P-value can also be found under:
    #p.value.vec[i] = cor.test(X[,i],Y,method = 'spearman')$p.value 
  }
  return(p.value.vec)
}

#create CLR transformations of the data:

#install.packages('compositions')
#install.packages('coda')
library(compositions)
clr_X_1 = compositions::clr(reads_counts+1)
clr_X_0_5 = compositions::clr(reads_counts+0.5)

#%%%%%%%%%%%%%%%%%
# Step VII: Run all methods, one-by-one, on each of the outputs.
#%%%%%%%%%%%%%%%%%

#these will hold the results objects.
res_DACOMP_objects = list()
res_ALDEx2_objects = list()
res_CSS_objects = list()
res_TSS_objects = list()
res_CLR_1_objects = list()
res_CLR_0_5_objects = list()

# Matrices for false and true positives. A true discovery is considered a discovery which is one of the spikeins.
# A false discovery is considered a discovery which is not one of the original spikeins.
# The term 'false' may be misleading, since the ground truth is not exacly known - we can still have contaminations toghether with out spikeins,
# or so the data suggests...

TP_Table = matrix(NA,nrow = length(Outcomes),ncol = 6)
FP_Table = matrix(NA,nrow = length(Outcomes),ncol = 6)
colnames(TP_Table) = c('DACOMP','ALDEx2','CSS','TSS','CLR_1','CLR_0_5')
colnames(FP_Table) = colnames(TP_Table)
rownames(TP_Table) = c('S1','S3','DilutionFactor')
rownames(FP_Table) = rownames(TP_Table)

#Iterate over outcomes:
for(outcome_id in 1:length(Outcomes)){
  print(paste0('Computing outcome id:',outcome_id))
  outcome = Outcomes[[outcome_id]]
  #and do testing, by method:
  print('DACOMP')
  set.seed(1)
  DACOMP_test = dacomp::dacomp.test(reads_counts,y = outcome,ind_reference_taxa = ref_obj,test = DACOMP.TEST.NAME.SPEARMAN,verbose = F,nr_perm = 10^4)
  res_DACOMP_objects[[outcome_id]] = which(p.adjust(DACOMP_test$p.values.test,method = 'BH')<=Q_LVL)
  
  TP_Table[outcome_id,1] = sum(res_DACOMP_objects[[outcome_id]] %in% spikeins_id)
  FP_Table[outcome_id,1] = length(res_DACOMP_objects[[outcome_id]]) - TP_Table[outcome_id,1]
  
  print('ALDEx2')
  dt_glm = data.frame(outcome=outcome)
  model_matrix_glm = model.matrix(~.,data = dt_glm)
  aldex_model_1 = aldex(reads = t(reads_counts),conditions = model_matrix_glm,test = 'glm',verbose = F)
  res_ALDEx2_objects[[outcome_id]] = which(p.adjust(aldex_model_1$model.outcome.Pr...t..,method = 'BH')<=Q_LVL)

  TP_Table[outcome_id,2] = sum(res_ALDEx2_objects[[outcome_id]] %in% spikeins_id)
  FP_Table[outcome_id,2] = length(res_ALDEx2_objects[[outcome_id]]) - TP_Table[outcome_id,2]

  
  print('CSS')
  pval_CSS = compute_spearman_marginal(X_CSS,outcome)
  res_CSS_objects[[outcome_id]] = which(p.adjust(pval_CSS,method = 'BH')<=Q_LVL)
  
  TP_Table[outcome_id,3] = sum(res_CSS_objects[[outcome_id]] %in% spikeins_id)
  FP_Table[outcome_id,3] = length(res_CSS_objects[[outcome_id]]) - TP_Table[outcome_id,3]
  
  
  print('TSS')
  pval_TSS = compute_spearman_marginal(X_TSS,outcome)
  res_TSS_objects[[outcome_id]] = which(p.adjust(pval_TSS,method = 'BH')<=Q_LVL)
  
  TP_Table[outcome_id,4] = sum(res_TSS_objects[[outcome_id]] %in% spikeins_id)
  FP_Table[outcome_id,4] = length(res_TSS_objects[[outcome_id]]) - TP_Table[outcome_id,4]
  
  print('CLR-1')
  pval_clr_1 = compute_spearman_marginal(clr_X_0_5,outcome)
  res_CLR_1_objects[[outcome_id]] = which(p.adjust(pval_clr_1,method = 'BH')<=Q_LVL)
  
  TP_Table[outcome_id,5] = sum(res_CLR_1_objects[[outcome_id]] %in% spikeins_id)
  FP_Table[outcome_id,5] = length(res_CLR_1_objects[[outcome_id]]) - TP_Table[outcome_id,5]
  
  
  print('CLR-0.5')
  pval_clr_0_5 = compute_spearman_marginal(clr_X_0_5,outcome)
  res_CLR_0_5_objects[[outcome_id]] = which(p.adjust(pval_clr_0_5,method = 'BH')<=Q_LVL)
  
  TP_Table[outcome_id,6] = sum(res_CLR_0_5_objects[[outcome_id]] %in% spikeins_id)
  FP_Table[outcome_id,6] = length(res_CLR_0_5_objects[[outcome_id]]) - TP_Table[outcome_id,6]
  
}


#%%%%%%%%%%%%%%%%%
# Step VIII: print tables for results, to be used at the paper
#%%%%%%%%%%%%%%%%%

library(xtable)
xtable(TP_Table)

xtable(FP_Table)

sink(file = '../../Results/Dilution_Results.txt')
xtable(TP_Table)

xtable(FP_Table)
sink()
stop('DONE')

#%%%%%%%%%%%%%%%%%
# Step IX: in - depth look at a taxon which is probably a contamination, see paper for exaplanation and details.
#%%%%%%%%%%%%%%%%%


cor.test(reads_counts[,spikeins_id[1]],reads_counts[,1561],method = 'spearman')
cor.test(reads_counts[,spikeins_id[2]],reads_counts[,1561],method = 'spearman')
cor.test(reads_counts[,spikeins_id[3]],reads_counts[,1561],method = 'spearman')

plot(reads_counts[,spikeins_id[2]],reads_counts[,1561])
plot(reads_counts[,spikeins_id[1]],reads_counts[,1561])

reads_taxonomy[1561]
