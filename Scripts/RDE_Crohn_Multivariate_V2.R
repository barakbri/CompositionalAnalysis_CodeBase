# This script runs the multivariate tests, over genera of taxa,
#for the Crohn's Disease dataset example.

# Function for Mahalanobis distance 
# This function was adapted from the SM of "Design of observational studies" by Paul Rosenbaum.
# The original code for this function was adapted from http://www-stat.wharton.upenn.edu/~rosenbap/software.html
# Under match functions for design of observational studies.

smahal = function (Xmat) 
{
  
  Xmat <- as.matrix(Xmat)
  n <- dim(Xmat)[1]
  rownames(Xmat) <- 1:n
  k <- dim(Xmat)[2]
  for (j in 1:k) Xmat[, j] <- rank(Xmat[, j])
  cv <- cov(Xmat)
  vuntied <- var(1:n)
  rat <- sqrt(vuntied/diag(cv))
  if(length(rat) > 1){
    diag_rat = diag(rat)
  }else{
    diag_rat = matrix(rat,nrow = 1,ncol = 1)
  }
  
  cv <- diag_rat %*% cv %*% diag_rat
  out <- matrix(NA, n, n)
  rownames(out) <- rownames(Xmat)
  colnames(out) <- rownames(Xmat)
  library(MASS)
  icov <- ginv(cv)
  for (i in 1:n) out[i, ] <- mahalanobis(Xmat, Xmat[i, ], icov, 
                                         inverted = T)
  out
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%
#Step I: Load the data
#%%%%%%%%%%%%%%%%%%%%%%%%%%
source('./RDE_Crohn_Load_Dataset.R')
OTHER = 2
prevalence_matrix = 1*(X>0)
to_remove = which(as.numeric(apply(prevalence_matrix, 2, sum))<OTHER)
X = X[,-to_remove]

load('../../Results/Gut_temp_v2_file_2.RData') # This is for S_crit = 1.3

#%%%%%%%%%%%%%%%%%%%%%%%%%%
#Step II: train a taxonomy classifier, and run classify the different sOTUs
#%%%%%%%%%%%%%%%%%%%%%%%%%%

#BiocManager::install('dada2')
library(dada2); packageVersion("dada2")
library(vegan)
NR_PERMS_MULTIVARIATE = 5000
Q_LVL = 0.1

set.seed(1) # Initialize random number generator for reproducibility
seqs = colnames(X)

taxa <- assignTaxonomy(seqs,
                       "../Crohn/gg_13_8_train_set_97.fa.gz",
                       multithread=T,
                       verbose = T,
                       minBoot = 80)

#%%%%%%%%%%%%%%%%%%%%%%%%%%
#Step III: check which genera have been identified in the data, and have more than 1 sOTU in them
#%%%%%%%%%%%%%%%%%%%%%%%%%%

genera_labels = unname(taxa)[,6]
genera_labels_to_test = (which(table(genera_labels)>1))
genera_labels_to_test = genera_labels_to_test[!(names(genera_labels_to_test) %in% c('g__'))]

genera_labels_to_test = names(genera_labels_to_test)

#number of sOTUs in genera for testing
#length(which((genera_labels %in% genera_labels_to_test)))


#%%%%%%%%%%%%%%%%%%%%%%%%%%
#Step IV: select references, and test the different genera
#%%%%%%%%%%%%%%%%%%%%%%%%%%

#select references
reference_obj = temp_obj$current_selected_ref_obj #Selected in the univariate example
#dacomp::dacomp.plot_reference_scores(reference_obj)

pval_rarefied             = rep(1,length(genera_labels_to_test))
pval_rarefied_amalgamated = rep(1,length(genera_labels_to_test))
pval_ratio                = rep(1,length(genera_labels_to_test))
pval_ratio_amalgamated    = rep(1,length(genera_labels_to_test))

set.seed(1)
for(genera_id_to_test in 1:length(genera_labels_to_test)){
  #genera_id_to_test = 2
  print(paste0('Testing genera ',genera_id_to_test,'/',length(genera_labels_to_test),': ', genera_labels_to_test[genera_id_to_test]))
  
  ind_reference = reference_obj$selected_references
  ind_genera = which(genera_labels == genera_labels_to_test[genera_id_to_test])
  ind_reference = ind_reference[!(ind_reference %in% ind_genera)]
  total_counts_in_reference = apply(X[,ind_reference],1,sum)
  
  X_genera_to_rarefy             = cbind(X[,ind_genera],total_counts_in_reference)
  X_genera_to_rarefy_amalgamated = cbind(apply(X[,ind_genera],1,sum),total_counts_in_reference)
  lambda_j = min(apply(X_genera_to_rarefy,1,sum))
  
  X_genera_rarefied = vegan::rrarefy(X_genera_to_rarefy,sample = lambda_j)
  X_genera_rarefied_amalgamated = vegan::rrarefy(X_genera_to_rarefy_amalgamated,sample = lambda_j)
  
  X_genera_TSS = X_genera_to_rarefy
  X_genera_TSS_amalgamated = X_genera_to_rarefy_amalgamated
  for(i in 1:nrow(X_genera_TSS)){
    X_genera_TSS[i,] = X_genera_TSS[i,]/sum(X_genera_TSS[i,])
    X_genera_TSS_amalgamated[i,] = X_genera_TSS_amalgamated[i,]/sum(X_genera_TSS_amalgamated[i,])
  }
  
  
  X_genera_rarefied = X_genera_rarefied[,-ncol(X_genera_to_rarefy),drop = F]
  X_genera_TSS = X_genera_TSS[,-ncol(X_genera_to_rarefy),drop = F]
  
  ind_to_keep = which(apply(X_genera_rarefied, 2, sum)>0)
  X_genera_rarefied = X_genera_rarefied[,ind_to_keep,drop = F]
  
  #X_genera_rarefied_amalgamated = X_genera_rarefied_amalgamated[,-ncol(X_genera_to_rarefy),drop = F]
  #X_genera_TSS_amalgamated = X_genera_TSS_amalgamated[,-ncol(X_genera_to_rarefy),drop = F]
  
  #Dist_L1_rarefied = dist(X_genera_rarefied,method = 'manhattan')
  #Dist_L1 = dist(X_genera_TSS,method = 'manhattan')
  if(!is.null(X_genera_rarefied) & ncol(X_genera_rarefied) > 0 ){
    Dist_L1_rarefied = smahal(X_genera_rarefied) 
  }
  Dist_L1 = smahal(X_genera_TSS)
  
  #Dist_L1_amalgamated_rarefied = dist(X_genera_rarefied_amalgamated,method = 'manhattan')
  #Dist_L1_amalgamated = dist(X_genera_TSS_amalgamated,method = 'manhattan')
  
  if(!is.null(X_genera_rarefied) & ncol(X_genera_rarefied) > 0){
    res_permanova = vegan::adonis(Dist_L1_rarefied~Y,permutations = NR_PERMS_MULTIVARIATE)
    pval_rarefied[genera_id_to_test] = res_permanova$aov.tab[1,6]  
  }
  
  res_permanova = vegan::adonis(Dist_L1~Y,permutations = NR_PERMS_MULTIVARIATE)
  pval_ratio[genera_id_to_test] = res_permanova$aov.tab[1,6]

  #res_permanova = vegan::adonis(Dist_L1_amalgamated_rarefied~Y,permutations = NR_PERMS_MULTIVARIATE)
  #pval_rarefied_amalgamated[genera_id_to_test] = res_permanova$aov.tab[1,6]

  #res_permanova = vegan::adonis(Dist_L1_amalgamated~Y,permutations = NR_PERMS_MULTIVARIATE)
  #pval_ratio_amalgamated[genera_id_to_test] = res_permanova$aov.tab[1,6]
  


  res_wilcoxon = subzero::PermSumCount.test(X = X_genera_rarefied_amalgamated[,1],Y = Y,B = NR_PERMS_MULTIVARIATE,DoWald = as.integer(0),return_perms = F,disable.ties.correction = F)
  pval_rarefied_amalgamated[genera_id_to_test] = res_wilcoxon$p.value

  res_wilcoxon = subzero::PermSumCount.test(X = X_genera_TSS_amalgamated[,1],Y = Y,B = NR_PERMS_MULTIVARIATE,DoWald = as.integer(0),return_perms = F,disable.ties.correction = F)
  pval_ratio_amalgamated[genera_id_to_test] = res_wilcoxon$p.value
}

pval_rarefied[is.na(pval_rarefied)] = 1
pval_rarefied_amalgamated[is.na(pval_rarefied_amalgamated)] = 1

adj_pval_rarefied = p.adjust(pval_rarefied,method = 'BH')
adj_pval_rarefied_amalgamated = p.adjust(pval_rarefied_amalgamated,method = 'BH')
adj_pval_ratio = p.adjust(pval_ratio,method = 'BH')
adj_pval_ratio_amalgamated = p.adjust(pval_ratio_amalgamated,method = 'BH')


#%%%%%%%%%%%%%%%%%%%%%%%%%%
#Step V: Analyze results:
#%%%%%%%%%%%%%%%%%%%%%%%%%%


# a list, containing the vecotrs of indices of discoveries, by method
disc_list = list(
  which(adj_pval_rarefied<=Q_LVL),
  which(adj_pval_ratio<=Q_LVL),
  which(adj_pval_rarefied_amalgamated<=Q_LVL),
  which(adj_pval_ratio_amalgamated<=Q_LVL)
)

#compute a matrix of shared discoveries. Diagonal entries are the number of discoveries of each method
method_names_rows = c('Multi.','Ratio, Multi.','Uni.','Ratio, Uni.')
method_names_col = method_names_rows#c('DACOMP, Mult.','DACOMP-Ratio, Mult.','DACOMP, Uni.','DACOMP-Ratio, Uni.')
shared_disc_mat = matrix(NA,nrow = length(method_names_rows),ncol = length(method_names_rows))
rownames(shared_disc_mat) = method_names_rows
colnames(shared_disc_mat) = method_names_col
for(i in 1:length(method_names_rows)){
  for(j in i:length(method_names_rows)){
    shared_disc_mat[i,j] = sum(disc_list[[i]] %in% disc_list[[j]])
  }
}

library(xtable)
xtable(shared_disc_mat)


