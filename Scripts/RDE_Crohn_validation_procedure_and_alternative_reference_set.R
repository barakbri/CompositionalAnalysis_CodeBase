###
# Run validation for references
###
temp_ref = ref_obj
sort(temp_ref$scores)[length(Selected_references)] #1.163254 
temp_ref$min_abundance_over_the_sorted[length(Selected_references) + (-1:1)]
(sort(temp_ref$scores))[length(Selected_references) + (-1:1)]
library(dacomp)
set.seed(1)
pvals_validation_ref = dacomp:::leave_one_out_validation(X = X,Y = Y,
                                                ref_obj_for_validation = Selected_references,
                                                Verbose = T,
                                                NR_PERMS = 1/(0.1/ncol(X)))

head(sort(p.adjust(pvals_validation_ref$pval.vec,method = 'BH')))

###
# Run DACOMP with other references
###

X_TSS = X
for(i in 1:nrow(X_TSS)){
  X_TSS[i,] = X_TSS[i,]/sum(X_TSS[i,])
}

plot(rank(apply(X_TSS,2,mean)),res_Wilcoxon_corrected$p.values)

other_references = which(apply(X_TSS>0,2,mean)>=0.1 & res_Wilcoxon_corrected$p.values>=0.3)
length(other_references)

hist(apply(X[,other_references],1,sum))
min(apply(X[,other_references],1,sum))



# run DACOMP tests
set.seed(1)
library(dacomp)
res_perm_Wilcoxon_other_ref = dacomp.test(X = X,y = Y,ind_reference_taxa = other_references,
                                verbose = T,
                                q = Q_LEVEL,
                                compute_ratio_normalization = T,
                                test = DACOMP.TEST.NAME.WILCOXON,disable_DSFDR = T,
                                nr_perm = 1000)

disc_other_ref = which(res_perm_Wilcoxon_other_ref$p.values.test.adjusted<=Q_LEVEL)
disc_DACOMP = which(p.adjust(res_Wilcoxon_list[[1]]$p.values.test,method = 'BH')<=Q_LEVEL)
length(disc_other_ref)
length(disc_DACOMP)
sum(disc_other_ref %in% disc_DACOMP)
sum(disc_other_ref %in% disc_Wilcoxon_corrected)


sum((disc_other_ref %in% ind_to_prevalent))
sum((disc_DACOMP %in% ind_to_prevalent))
sum(!(disc_other_ref %in% ind_to_prevalent))
sum(!(disc_DACOMP %in% ind_to_prevalent))

sum(!(disc_Wilcoxon_corrected %in% disc_other_ref) & !(disc_Wilcoxon_corrected %in% ind_to_prevalent))

disc_other_ref = which(res_perm_Wilcoxon_other_ref$p.values.test.adjusted.ratio.normalization<=Q_LEVEL)
disc_DACOMP = which(p.adjust(res_Wilcoxon_list[[1]]$p.values.test.ratio.normalization,method = 'BH')<=Q_LEVEL)
length(disc_other_ref)
length(disc_DACOMP)
sum(disc_other_ref %in% disc_DACOMP)
sum(disc_other_ref %in% disc_Wilcoxon_corrected)

sum((disc_other_ref %in% ind_to_prevalent))
sum((disc_DACOMP %in% ind_to_prevalent))
sum(!(disc_other_ref %in% ind_to_prevalent))
sum(!(disc_DACOMP %in% ind_to_prevalent))

sum(!(disc_Wilcoxon_corrected %in% disc_other_ref) & !(disc_Wilcoxon_corrected %in% ind_to_prevalent))
sum(!(which(!(1:ncol(X) %in% ind_to_prevalent)) %in% other_references))
