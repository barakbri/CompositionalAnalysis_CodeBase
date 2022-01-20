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

load('../../Results/Gut_temp_v2_file_dacomp.RData')

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
#Step III: run msa
#%%%%%%%%%%%%%%%%%%%%%%%%%%

#try unifrac
# using https://f1000research.com/articles/5-1492/v1
library(msa)
seqs <- colnames(otu_table)
names(seqs) <- seqs # This propagates to the tip labels of the tree
mult <- msa(seqs, method="ClustalW", type="dna", order="input")
#install.packages('phangorn')
library("phangorn")
phang.align <- as.phyDat(mult, type="DNA", names=getSequence(seqtab))
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
## negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
#detach("package:phangorn", unload=TRUE)

library(phyloseq)

#%%%%%%%%%%%%%%%%%%%%%%%%%%
#Step IV: check which genera have been identified in the data, and have more than 1 sOTU in them
#%%%%%%%%%%%%%%%%%%%%%%%%%%
Agg_level = 6
genera_labels = unname(taxa)[,Agg_level]
genera_labels_to_test = (which(table(genera_labels)>1))
genera_labels_to_test = genera_labels_to_test[!(names(genera_labels_to_test) %in% c('g__'))]
genera_labels_to_test = names(genera_labels_to_test)

#number of sOTUs in genera for testing
#length(which((genera_labels %in% genera_labels_to_test)))
#temp = sort(table(genera_labels[genera_labels %in% genera_labels_to_test]))
#median(temp)

#%%%%%%%%%%%%%%%%%%%%%%%%%%
#Step V: select references, and test the different genera
#%%%%%%%%%%%%%%%%%%%%%%%%%%

#select references
reference_taxa = temp_obj$Selected_references #Selected in the univariate example
#dacomp::dacomp.plot_reference_scores(reference_obj)

pval_rarefied             = rep(1,length(genera_labels_to_test))
pval_rarefied_amalgamated = rep(1,length(genera_labels_to_test))
pval_ratio                = rep(1,length(genera_labels_to_test))
pval_ratio_amalgamated    = rep(1,length(genera_labels_to_test))
Nr_Samples_with_zero_in_genus = rep(1,length(genera_labels_to_test))

Subvectors_rarefied = list()
Subvectors_TSS = list()
DIST_rarefied = list()
DIST_TSS = list()

set.seed(1)
for(genera_id_to_test in 1:length(genera_labels_to_test)){
  #genera_id_to_test = 2
  print(paste0('Testing genera ',genera_id_to_test,'/',length(genera_labels_to_test),': ', genera_labels_to_test[genera_id_to_test]))
  
  ind_reference =reference_taxa
  ind_genera = which(genera_labels == genera_labels_to_test[genera_id_to_test])
  ind_reference = ind_reference[!(ind_reference %in% ind_genera)]
  total_counts_in_reference = apply(X[,ind_reference],1,sum)
  Nr_Samples_with_zero_in_genus[genera_id_to_test] = sum(apply(X[,ind_genera],1,sum)==0)
  
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
  
  Subvectors_rarefied[[genera_id_to_test]] = X_genera_rarefied
  Subvectors_TSS[[genera_id_to_test]] = X_genera_TSS
  
  X_genera_rarefied = X_genera_rarefied[,-ncol(X_genera_to_rarefy),drop = F]
  X_genera_TSS = X_genera_TSS[,-ncol(X_genera_to_rarefy),drop = F]
  
  ind_to_keep = which(apply(X_genera_rarefied, 2, sum)>0)
  X_genera_rarefied = X_genera_rarefied[,ind_to_keep,drop = F]
  
  if(!is.null(X_genera_rarefied) & ncol(X_genera_rarefied) > 0 ){
    Dist_L1_rarefied = smahal(X_genera_rarefied) 
  }
  Dist_L1 = smahal(X_genera_TSS)
  DIST_rarefied[[genera_id_to_test]] = Dist_L1_rarefied
  DIST_TSS[[genera_id_to_test]] = Dist_L1
  
  if(!is.null(X_genera_rarefied) & ncol(X_genera_rarefied) > 0){
    res_permanova = vegan::adonis(Dist_L1_rarefied~Y,permutations = NR_PERMS_MULTIVARIATE)
    pval_rarefied[genera_id_to_test] = res_permanova$aov.tab[1,6]  
  }
  
  res_permanova = vegan::adonis(Dist_L1~Y,permutations = NR_PERMS_MULTIVARIATE)
  pval_ratio[genera_id_to_test] = res_permanova$aov.tab[1,6]
  
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
#Step VI: Analyze results:
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

sink(file =  paste0('../../Results/Crohn_MV_res_',Agg_level,'.txt'))
print(shared_disc_mat)
sink()

#summary statistics for the number of zeros, for section 2
length(Nr_Samples_with_zero_in_genus)
sum(Nr_Samples_with_zero_in_genus==0)
sum(Nr_Samples_with_zero_in_genus>=10)
sum(Nr_Samples_with_zero_in_genus>=30)
sum(Nr_Samples_with_zero_in_genus>=1)
output_discoveries = data.frame(Genera_Names = genera_labels_to_test,
                                Discovered_by_Multivariate = rep(F,length(genera_labels_to_test)),
                                Discovered_by_Multivariate_Ratio = rep(F,length(genera_labels_to_test)),
                                Discovered_by_Univariate = rep(F,length(adj_pval_rarefied_amalgamated)),
                                Discovered_by_Univariate_Ratio = rep(F,length(genera_labels_to_test)))
output_discoveries$Discovered_by_Multivariate[which(adj_pval_rarefied<=Q_LVL)] = T
output_discoveries$Discovered_by_Multivariate_Ratio[which(adj_pval_ratio<=Q_LVL)] = T
output_discoveries$Discovered_by_Univariate[which(adj_pval_rarefied_amalgamated<=Q_LVL)] = T
output_discoveries$Discovered_by_Univariate_Ratio[which(adj_pval_ratio_amalgamated<=Q_LVL)] = T

#write.csv(output_discoveries,file = '../../Results/multivariate_discoveries.csv',quote = F,row.names = F)
write.csv(output_discoveries,file = paste0('../../Results/multivariate_discoveries_',Agg_level,'.csv'),quote = F,row.names = F)


rarefaction_Mult_alone = which(output_discoveries$Discovered_by_Multivariate &
                                 !output_discoveries$Discovered_by_Univariate)
rarefaction_Uni_alone = which(!output_discoveries$Discovered_by_Multivariate &
                                output_discoveries$Discovered_by_Univariate)
rarefaction_Uni_and_Mult = which(output_discoveries$Discovered_by_Multivariate &
                                   output_discoveries$Discovered_by_Univariate)
rarefaction_None = which(!output_discoveries$Discovered_by_Multivariate &
                           !output_discoveries$Discovered_by_Univariate)



#%%%%%%%%%%%%%%%%%%%%%%%%%%
#Step VII: create plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%
set.seed(1)
source('./RDE_Crohn_Multivariate_PCOA_functions.R')

library(ape)

prevalence = apply(X>0,2,mean)
output_univariate_discoveries = read.csv(file = '../../Results/univariate_discoveries.csv')

discovered_by_at_least_one = which(output_univariate_discoveries$ANCOM|
output_univariate_discoveries$WIL.CSS|
output_univariate_discoveries$DACOMP|
output_univariate_discoveries$DACOMP.ratio|
output_univariate_discoveries$ALDEx2.t|
output_univariate_discoveries$WIL.FLOW)


TAXON_PREVALENCE_RANK_TO_PICK = 1

X_rarefied = vegan::rrarefy(X,sample = min(apply(X,1,sum)))  

res_unifrac = run_configuration_with_raw(method = 'unifrac',remove_uni_discoveries = F,modes = 1:4)
res_Rmahal = run_configuration_with_raw(method = 'Rmahal',remove_uni_discoveries = F,modes = 1:4)
res_BC = run_configuration_with_raw(method = 'BC',remove_uni_discoveries = F,modes = 1:4)

res_unifrac_filtered = run_configuration_with_raw(method = 'unifrac',remove_uni_discoveries = T,modes = 1:5)
res_Rmahal_filtered = run_configuration_with_raw(method = 'Rmahal',remove_uni_discoveries = T,modes = 1:5)
res_BC_filtered = run_configuration_with_raw(method = 'BC',remove_uni_discoveries = T,modes = 1:5)

res_ratio = run_configuration_with_raw(method = 'ratio',remove_uni_discoveries = F,modes = 1:5,disable_numbers_in_ratio_plot = T)

res_L1 = run_configuration_with_raw(method = 'L1',remove_uni_discoveries = F,modes = 1:4,do_not_include_ref_in_unnormalized = T)

#plots for each group by BC
perm_for_draw = sample(1:length(Y))
PCH_vec = rep('+',length(Y))
PCH_vec[Y==1] = 'x'
ASP = NA
pdf(file = '../../Results/PCoA_mult_alone.pdf',height = 2/3*7)
ind_to_take =rarefaction_Mult_alone
par(mfrow = c(2,3))
for(i in 1:length(ind_to_take)){
  if(!is.null(res_BC[[ind_to_take[i]]]$X_1)){
    plot(res_BC[[ind_to_take[i]]]$X_1[perm_for_draw], 
         res_BC[[ind_to_take[i]]]$Y_1[perm_for_draw], 
         col = Y[perm_for_draw]+1,asp=ASP, pch=PCH_vec[perm_for_draw],cex.axis = 0.8, cex = 0.8, xlab = 'PCOA axis 1', ylab = 'PCOA axis 2',main = paste0(genera_labels_to_test[ind_to_take[i]]) )
    abline(h=0,col='blue',lty=2)
    abline(v=0,col='blue',lty=2)
  }
}
dev.off()
par(mfrow = c(1,1))


pdf(file = '../../Results/PCoA_uni_alone.pdf')
ind_to_take =rarefaction_Uni_alone
par(mfrow = c(3,3))
for(i in 1:length(ind_to_take)){
  if(!is.null(res_BC[[ind_to_take[i]]]$X_1)){
    plot(res_BC[[ind_to_take[i]]]$X_1[perm_for_draw], 
         res_BC[[ind_to_take[i]]]$Y_1[perm_for_draw], 
         col = Y[perm_for_draw]+1,asp=ASP, pch=PCH_vec[perm_for_draw],cex.axis = 0.8, cex = 0.8, xlab = 'PCOA axis 1', ylab = 'PCOA axis 2',main = paste0(genera_labels_to_test[ind_to_take[i]]) )
    abline(h=0,col='blue',lty=2)
    abline(v=0,col='blue',lty=2)
  }
}
dev.off()
par(mfrow = c(1,1))

pdf(file = '../../Results/PCoA_uni_and_mult.pdf',height = 12,width = 12)
ind_to_take =rarefaction_Uni_and_Mult
par(mfrow = c(4,4))
for(i in 1:length(ind_to_take)){
  if(!is.null(res_BC[[ind_to_take[i]]]$X_1)){
    plot(res_BC[[ind_to_take[i]]]$X_1[perm_for_draw], 
         res_BC[[ind_to_take[i]]]$Y_1[perm_for_draw], 
         col = Y[perm_for_draw]+1,asp=ASP, pch=PCH_vec[perm_for_draw],cex.axis = 0.8, cex = 0.8, xlab = 'PCOA axis 1', ylab = 'PCOA axis 2',main = paste0(genera_labels_to_test[ind_to_take[i]]) )
    abline(h=0,col='blue',lty=2)
    abline(v=0,col='blue',lty=2)
  }
}
dev.off()
par(mfrow = c(1,1))



#Masked signal - g__Oscillospira, g__Dorea, g__Clostridium - both correct and incorrect analysis


pdf(file = '../../Results/PCoA_masked_signal.pdf',width = 7,height = 7)
ind_to_take = which(genera_labels_to_test %in% c('g__Oscillospira', 'g__Clostridium'))
par(mfrow = c(length(ind_to_take),2))
for(i in 1:length(ind_to_take)){
  plot(res_BC[[ind_to_take[i]]]$X_1, 
       res_BC[[ind_to_take[i]]]$Y_1, 
       col = Y+1,asp=ASP, pch=PCH_vec,cex.axis = 0.8, cex = 0.8, xlab = 'PCOA axis 1', ylab = 'PCOA axis 2',main = paste0(genera_labels_to_test[ind_to_take[i]],', correct analysis: ') )
  abline(h=0,col='blue',lty=2)
  abline(v=0,col='blue',lty=2)
  plot(res_L1[[ind_to_take[i]]]$X_2, 
       res_L1[[ind_to_take[i]]]$Y_2, 
       col = Y+1,asp=ASP, pch=PCH_vec,cex.axis = 0.8, cex = 0.8, xlab = 'PCOA axis 1', ylab = 'PCOA axis 2',main = 'naive PCoA:' )
  abline(h=0,col='blue',lty=2)
  abline(v=0,col='blue',lty=2)
}
dev.off()
par(mfrow = c(1,1))


#None with fake signals: g__[Eubacterium], g__[Prevotella]
pdf(file = '../../Results/PCoA_fake_signal.pdf',height = 3.5)
ind_to_take = which(genera_labels_to_test %in% c('g__[Eubacterium]'))
par(mfrow = c(length(ind_to_take),2))
for(i in 1:length(ind_to_take)){
  plot(res_Rmahal[[ind_to_take[i]]]$X_1, 
       res_Rmahal[[ind_to_take[i]]]$Y_1, 
       col = Y+1, pch=PCH_vec,cex.axis = 0.8,asp=ASP, cex = 0.8, xlab = 'PCOA axis 1', ylab = 'PCOA axis 2',main = paste0(genera_labels_to_test[ind_to_take[i]],',\n\r correct analysis:') )
  abline(h=0,col='blue',lty=2)
  abline(v=0,col='blue',lty=2)
  plot(res_Rmahal[[ind_to_take[i]]]$X_2, 
       res_Rmahal[[ind_to_take[i]]]$Y_2, 
       col = Y+1, pch=PCH_vec,cex.axis = 0.8,asp=ASP, cex = 0.8, xlab = 'PCOA axis 1', ylab = 'PCOA axis 2',main = 'naive PCoA:' )
  abline(h=0,col='blue',lty=2)
  abline(v=0,col='blue',lty=2)
}
dev.off()
par(mfrow = c(1,1))



#Removal of discoveries - still see signal.
#g__[Ruminococcus],g__Blautia,g__Ruminococcus
ind_to_take = sort(unique(c(rarefaction_Mult_alone,rarefaction_Uni_alone,rarefaction_Uni_and_Mult)))
ind_to_take = c(12,34,43)#c(3,12,15,22,25,29,31,32,34,43,52)
pdf(file = '../../Results/PCoA_after_ASV_discoveries_filtered.pdf',height = 12*3/4,width = 12)
par(mfrow = c(1,3))
for(i in 1:length(ind_to_take)){
  if(!is.null(res_BC_filtered[[ind_to_take[i]]]$X_1)){
    plot(res_BC_filtered[[ind_to_take[i]]]$X_1, 
         res_BC_filtered[[ind_to_take[i]]]$Y_1, 
         col = Y+1,asp=ASP,cex.axis = 0.8, pch=PCH_vec, cex = 0.8, xlab = 'PCOA axis 1', ylab = 'PCOA axis 2',main = paste0(genera_labels_to_test[ind_to_take[i]]) )
    abline(h=0,col='blue',lty=2)
    abline(v=0,col='blue',lty=2)
  }
}
dev.off()
par(mfrow = c(1,1))

pdf(file = paste0('../../Results/PCoA_after_ASV_discoveries_filtered_comparison.pdf'),height = 8*2/3,width = 8)
par(mfcol = c(2,3))
subfigure_labels = c('A','B','C','D','E','F')
subfigure_ind = 1
for(i in  1:length(ind_to_take)){
  if(!is.null(res_BC_filtered[[ind_to_take[i]]]$X_1)){
    plot(res_BC_filtered[[ind_to_take[i]]]$X_1, 
         res_BC_filtered[[ind_to_take[i]]]$Y_1, 
         col = Y+1,asp=ASP,cex.axis = 0.8, pch=PCH_vec, cex = 0.8, xlab = 'PCOA axis 1', ylab = 'PCOA axis 2',main = paste0(subfigure_labels[subfigure_ind],') ',genera_labels_to_test[ind_to_take[i]],'\n\r filtered ASV discoveries') )
    abline(h=0,col='blue',lty=2)
    abline(v=0,col='blue',lty=2)
    subfigure_ind = subfigure_ind +1
    plot(res_BC[[ind_to_take[i]]]$X_1, 
         res_BC[[ind_to_take[i]]]$Y_1, 
         col = Y+1,asp=ASP,cex.axis = 0.8, pch=PCH_vec, cex = 0.8, xlab = 'PCOA axis 1', ylab = 'PCOA axis 2',main = paste0(subfigure_labels[subfigure_ind],') ',genera_labels_to_test[ind_to_take[i]],'\n\r all ASVs') )
    abline(h=0,col='blue',lty=2)
    abline(v=0,col='blue',lty=2)
    subfigure_ind = subfigure_ind +1
  }
}
dev.off()
par(mfrow = c(1,1))

shuffling_perm = sample(1:length(Y))
pdf(file = paste0('../../Results/PCoA_after_ASV_discoveries_filtered_comparison_permuted.pdf'),height = 8*2/3,width = 8)
par(mfcol = c(2,3))
for(i in  1:length(ind_to_take)){
  if(!is.null(res_BC_filtered[[ind_to_take[i]]]$X_1)){
    plot(res_BC_filtered[[ind_to_take[i]]]$X_1, 
         res_BC_filtered[[ind_to_take[i]]]$Y_1, 
         col = (Y)[shuffling_perm]+1,asp=ASP,cex.axis = 0.8, pch=PCH_vec[shuffling_perm], cex = 0.8, xlab = 'PCOA axis 1', ylab = 'PCOA axis 2',main = paste0(genera_labels_to_test[ind_to_take[i]],'\n\r filtered ASV discoveries\n\r labels shuffled') )
    abline(h=0,col='blue',lty=2)
    abline(v=0,col='blue',lty=2)
    plot(res_BC[[ind_to_take[i]]]$X_1, 
         res_BC[[ind_to_take[i]]]$Y_1, 
         col = (Y)[shuffling_perm]+1,asp=ASP,cex.axis = 0.8, pch=PCH_vec[shuffling_perm], cex = 0.8, xlab = 'PCOA axis 1', ylab = 'PCOA axis 2',main = paste0(genera_labels_to_test[ind_to_take[i]],'\n\r all ASVs\n\r labels shuffled') )
    abline(h=0,col='blue',lty=2)
    abline(v=0,col='blue',lty=2)
  }
}
dev.off()
par(mfrow = c(1,1))

pdf(file = '../../Results/PCoA_other_metrics.pdf',height = 5,width = 10)
par(mfrow=c(1,2))
plot(res_unifrac[[34]]$X_1, 
     res_unifrac[[34]]$Y_1, 
     col = Y+1,asp=ASP, pch=PCH_vec,cex.axis = 0.8, cex = 0.8, xlab = 'PCOA axis 1', ylab = 'PCOA axis 2',main = paste0('A) \"within-genus\" PCOA for \n',genera_labels_to_test[34],'\n using the Unifrac metric' ))
abline(h=0,col='blue',lty=2)
abline(v=0,col='blue',lty=2)


plot(res_Rmahal[[27]]$X_1, 
     res_Rmahal[[27]]$Y_1, 
     col = Y+1,asp=ASP, pch=PCH_vec,cex.axis = 0.8, cex = 0.8, xlab = 'PCOA axis 1', ylab = 'PCOA axis 2',main = paste0('B) \"within-genus\" PCOA for \n',genera_labels_to_test[27],'\n using the robust Mahalanobis metric') )
abline(h=0,col='blue',lty=2)
abline(v=0,col='blue',lty=2)
dev.off()
par(mfrow=c(1,1))

col_vector_GrayScale = as.character(Y+1)
col_vector_GrayScale[col_vector_GrayScale == '1'] = 'dimgray'
col_vector_GrayScale[col_vector_GrayScale == '2'] = 'black'

pdf(file = '../../Results/PCoA_examples_for_paper.pdf',height = 6,width = 6)
par(mfrow=c(2,2))
plot(res_BC[[3]]$X_1, 
     res_BC[[3]]$Y_1, 
     col = col_vector_GrayScale,asp=ASP, pch=PCH_vec, cex = 0.8, xlab = 'PCOA axis 1',
     ylab = 'PCOA axis 2',main = paste0('A) ',genera_labels_to_test[3]),cex.axis = 0.8)
abline(h=0,col='gray',lty=2)
abline(v=0,col='gray',lty=2)

plot(res_BC[[15]]$X_1, 
     res_BC[[15]]$Y_1, 
     col = col_vector_GrayScale,asp=ASP, pch=PCH_vec, cex = 0.8, xlab = 'PCOA axis 1',
     ylab = 'PCOA axis 2',main = paste0('B) ',genera_labels_to_test[15]),cex.axis = 0.8)
abline(h=0,col='gray',lty=2)
abline(v=0,col='gray',lty=2)

boxplot(res_ratio[[3]],main = paste0('C) ',genera_labels_to_test[3]),ylim = c(0,1),cex.axis = 0.8)

boxplot(res_ratio[[15]],main = paste0('D) ',genera_labels_to_test[15]),ylim = c(0,1),cex.axis = 0.8)
dev.off()


#%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step VIII : permanova with rarefaction over all taxa and then testing
#%%%%%%%%%%%%%%%%%%%%%%%%%%
PV_naive_sum = rep(NA,length(genera_labels_to_test))
PV_naive_mult = rep(NA,length(genera_labels_to_test))
set.seed(1)
for(genera_id_to_test in 1:length(genera_labels_to_test)){
  #genera_id_to_test = 2
  print(paste0('Testing genera ',genera_id_to_test,'/',length(genera_labels_to_test),': ', genera_labels_to_test[genera_id_to_test]))
  
  ind_reference =reference_taxa
  ind_genera = which(genera_labels == genera_labels_to_test[genera_id_to_test])
  ind_reference = ind_reference[!(ind_reference %in% ind_genera)]
  total_counts_in_reference = apply(X_rarefied[,ind_reference],1,sum)
  
  X_sum = apply(X_rarefied[,ind_genera],1,sum)
  X_genera = cbind(X_rarefied[,ind_genera],total_counts_in_reference)
  Dist_L1 = smahal(X_genera)
  res_permanova = vegan::adonis(Dist_L1~Y,permutations = NR_PERMS_MULTIVARIATE)
  PV_naive_mult[genera_id_to_test] = res_permanova$aov.tab[1,6]  
  
  res_wilcoxon = subzero::PermSumCount.test(X = X_sum,Y = Y,B = NR_PERMS_MULTIVARIATE,DoWald = as.integer(0),return_perms = F,disable.ties.correction = F)
  PV_naive_sum[genera_id_to_test] = res_wilcoxon$p.value
  
}

PV_naive_mult_adjusted = p.adjust(PV_naive_mult,method = 'BH')
PV_naive_sum_adjusted = p.adjust(PV_naive_sum,method = 'BH')
#Compare P-values
genera_labels_to_test[c(22,29,43)]
PV_naive_mult_adjusted[c(22,29,43)]
adj_pval_rarefied[c(22,29,43)]
PV_naive_sum_adjusted[c(22,29,43)]
adj_pval_rarefied_amalgamated[c(22,29,43)]
