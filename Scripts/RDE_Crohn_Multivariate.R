# This script runs the multivariate tests, over genera of taxa,
#for the Crohn's Disease dataset example.

smahal = function (X) 
{
  X <- as.matrix(X)
  n <- dim(X)[1]
  rownames(X) <- 1:n
  k <- dim(X)[2]
  for (j in 1:k) X[, j] <- rank(X[, j])
  cv <- cov(X)
  vuntied <- var(1:n)
  rat <- sqrt(vuntied/diag(cv))
  cv <- diag(rat) %*% cv %*% diag(rat)
  out <- matrix(NA, n, n)
  rownames(out) <- rownames(X)
  colnames(out) <- rownames(X)
  library(MASS)
  icov <- ginv(cv)
  for (i in 1:m) out[i, ] <- mahalanobis(X, X[i, ], icov, 
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

pval_rarefied             = rep(NA,length(genera_labels_to_test))
pval_rarefied_amalgamated = rep(NA,length(genera_labels_to_test))
pval_ratio                = rep(NA,length(genera_labels_to_test))
pval_ratio_amalgamated    = rep(NA,length(genera_labels_to_test))

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
  
  #X_genera_rarefied = X_genera_rarefied[,-ncol(X_genera_to_rarefy),drop = F]
  #X_genera_TSS = X_genera_TSS[,-ncol(X_genera_to_rarefy),drop = F]
  X_genera_rarefied_amalgamated = X_genera_rarefied_amalgamated[,-ncol(X_genera_to_rarefy),drop = F]
  X_genera_TSS_amalgamated = X_genera_TSS_amalgamated[,-ncol(X_genera_to_rarefy),drop = F]
  
  Dist_L1_rarefied = dist(X_genera_rarefied,method = 'manhattan')
  Dist_L1 = dist(X_genera_TSS,method = 'manhattan')
  Dist_L1_amalgamated_rarefied = dist(X_genera_rarefied_amalgamated,method = 'manhattan')
  Dist_L1_amalgamated = dist(X_genera_TSS_amalgamated,method = 'manhattan')
  
  
  res_permanova = vegan::adonis(Dist_L1_rarefied~Y,permutations = NR_PERMS_MULTIVARIATE)
  pval_rarefied[genera_id_to_test] = res_permanova$aov.tab[1,6]
  
  res_permanova = vegan::adonis(Dist_L1~Y,permutations = NR_PERMS_MULTIVARIATE)
  pval_ratio[genera_id_to_test] = res_permanova$aov.tab[1,6]

  res_permanova = vegan::adonis(Dist_L1_amalgamated_rarefied~Y,permutations = NR_PERMS_MULTIVARIATE)
  pval_rarefied_amalgamated[genera_id_to_test] = res_permanova$aov.tab[1,6]

  res_permanova = vegan::adonis(Dist_L1_amalgamated~Y,permutations = NR_PERMS_MULTIVARIATE)
  pval_ratio_amalgamated[genera_id_to_test] = res_permanova$aov.tab[1,6]
  
  #
  # 
  # res_wilcoxon = subzero::PermSumCount.test(X = X_genera_rarefied_amalgamated[,1],Y = Y,B = NR_PERMS_MULTIVARIATE,DoWald = as.integer(0),return_perms = F,disable.ties.correction = F)
  # pval_rarefied_amalgamated[genera_id_to_test] = res_wilcoxon$p.value
  # 
  # res_wilcoxon = subzero::PermSumCount.test(X = X_genera_TSS_amalgamated[,1],Y = Y,B = NR_PERMS_MULTIVARIATE,DoWald = as.integer(0),return_perms = F,disable.ties.correction = F)
  # pval_ratio_amalgamated[genera_id_to_test] = res_wilcoxon$p.value
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

#plot the rejections by test
plot(-log(adj_pval_rarefied),-log(adj_pval_rarefied_amalgamated))
abline(h = -log(Q_LVL),col = 'red')
abline(a = 0,b=1,col = 'red')
abline(v = -log(Q_LVL),col = 'red')

sum(adj_pval_rarefied<=Q_LVL) #rejections for rarefied multivariate
sum(adj_pval_rarefied_amalgamated<=Q_LVL) #rarefied univariate
sum(adj_pval_rarefied<=Q_LVL & adj_pval_rarefied_amalgamated >Q_LVL) # multivariate discoveries that are unique
sum(adj_pval_rarefied<=adj_pval_rarefied_amalgamated &
      adj_pval_rarefied_amalgamated <=Q_LVL) #pvalues lower in multivariate than in univariate

#discoveries for the ratio test
sum(adj_pval_ratio<=Q_LVL)
sum(adj_pval_ratio_amalgamated<=Q_LVL)
sum(adj_pval_ratio<=Q_LVL & adj_pval_ratio_amalgamated >Q_LVL)

#sizes of genera, in terms of OTUs
gen_size = rep(NA,length(genera_labels_to_test))
for(genera_id_to_test in 1:length(genera_labels_to_test))
  gen_size[genera_id_to_test] = length(which(genera_labels == genera_labels_to_test[genera_id_to_test]))

median(gen_size) 
mean(gen_size)
