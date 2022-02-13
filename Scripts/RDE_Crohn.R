# Script for processing the crohn dataset, from Van De Putte et al. (2017)
# The data analyzed is the data from the third analysis done in the paper.
source('./RDE_Crohn_Load_Dataset.R')

# median number of counts, for OTUs - sequencing depths are similar across groups
median(sort(as.numeric(apply(X,1,sum))))
#[1] 20437
mean(sort(as.numeric(apply(X,1,sum))))
#[1] 22886

#prevalence of taxa across subjects
prevalence_matrix = 1*(X>0)
prev_taxa = as.numeric(apply(prevalence_matrix,2,mean))
hist(prev_taxa)

X_as_percentage = X
for(i in 1:nrow(X_as_percentage)){
  X_as_percentage[i,] = (X_as_percentage[i,])/sum(X_as_percentage[i,])
}
mean_percent_taxa  = as.numeric(apply(X_as_percentage,2,mean))


#correct by volume counts, according to the method of Van De Putte
library(phyloseq)
source('./../Packages/QMP-master/QMP-master/QMP.R')

rownames(X)
Average_Cell_Count_matrix = matrix(Average_Cell_Count,ncol = 1)
rownames(Average_Cell_Count_matrix) = rownames(X)
X_corrected = rarefy_even_sampling_depth(X,Average_Cell_Count_matrix)

# remove taxa in less than 2 people. Note that this is done after correction because the Van De Putte correction uses the original samples library size
OTHER = 2
to_remove = which(as.numeric(apply(prevalence_matrix, 2, sum))<OTHER)
X = X[,-to_remove]
X_corrected = X_corrected[,-to_remove]

# number of permutations for computing P-values 
NR.PERMS = 100000

#Wilcoxon over the corrected (by flow cytometry)
source('Wilcoxon_TaxaWise.R')
res_Wilcoxon_corrected = wilcoxon_taxa_wise(X_corrected,y = Y,nr.perms = NR.PERMS)
hist(res_Wilcoxon_corrected$p.values)
total_counts_in_taxa = as.numeric(apply(prevalence_matrix,2,sum))

disc_Wilcoxon_corrected = which(p.adjust(res_Wilcoxon_corrected$p.values,method = 'BH')<=Q_LEVEL)
length(disc_Wilcoxon_corrected) 


#Run ancom, and define a function for comparing intersection of results with ANCOM
library(ancom.R)
ANCOM_otu_dat = X
ANCOM_otu_dat = cbind(ANCOM_otu_dat,Y)
ANCOM_otu_dat = data.frame(ANCOM_otu_dat)
ANCOM_default_res = ancom.R::ANCOM(OTUdat = ANCOM_otu_dat,multcorr = 3,sig = Q_LEVEL)

intersect_with_ANCOM = function(ANCOM_res,dacomp_rejections){
  print(paste0('NUMBER OF ANCOM REJ:',length(ANCOM_res$detected)))
  disc_ANCOM =which(colnames(ANCOM_otu_dat) %in% ANCOM_res$detected)
  print(paste0('ANCOM in WILCOX corrected:',sum(disc_ANCOM %in% disc_Wilcoxon_corrected)))
  print(paste0('DACOMP in ANCOM corrected:',sum(dacomp_rejections %in% disc_ANCOM)))
  
  possibly_novel_disc = dacomp_rejections[which(!(dacomp_rejections %in% disc_ANCOM))]
  print(paste0('Verified, Unique:',sum(possibly_novel_disc %in% disc_Wilcoxon_corrected),'/',length(possibly_novel_disc)))
  ANCOM_TO_INVESTIGATE = disc_ANCOM[which(!(disc_ANCOM %in% dacomp_rejections ))]
  return(ANCOM_TO_INVESTIGATE)
}

# Run DACOMP across different median SD thresholds (defined in median_SD_thres_Vec)
set.seed(1)
library(subzero)
library(dacomp)
PS_value = 1 # Pseudocount value used

set.seed(1)

#stored file path
filepath = '../../Results/Gut_temp_v2_file_dacomp.RData'

# select references

print(paste0('Selecting references'))
ref_obj = dacomp.select_references(X = X,median_SD_threshold = 1,maximal_TA = 200,Pseudo_Count_used = PS_value,  verbose = F)
nr_taxa_to_select_as_ref = sum(ref_obj$min_abundance_over_the_sorted<=100)+1
Scrit = sort(ref_obj$scores)[nr_taxa_to_select_as_ref] #1.163254 
Selected_references = which(ref_obj$scores<=Scrit)

library(latex2exp)
pdf('../../Results/Gut_Scrit_select.pdf',width = 5,height = 5/3,pointsize = 8)
par(mfrow = c(1,3))
plot(sort(ref_obj$scores),ref_obj$min_abundance_over_the_sorted,xlab = TeX("S_{crit}"),ylab = "Min. counts in reference",pch=20,cex=0.3,main = 'A')
abline(v = Scrit,col = 'azure4',lty=2)
plot(sort(ref_obj$scores),ref_obj$min_abundance_over_the_sorted,xlim = c(1,1.5),xlab = TeX("S_{crit}"),ylab = "Min. counts in reference",ylim = c(0,300),pch=20,cex=0.3,main = 'B')
abline(v = Scrit,col = 'azure4',lty=2)
hist(ref_obj$scores,main = 'C',breaks = 50,xlab = 'Reference selection scores')
dev.off()
par(mfrow = c(1,1))

pdf('../../Results/Gut_Scrit_select_for_paper.pdf',height = 3.5,pointsize = 8)
par(mfrow = c(1,2))
plot(sort(ref_obj$scores),ref_obj$min_abundance_over_the_sorted,xlab = TeX("S_{crit}"),ylab = "Min. counts in reference",pch=20,cex=0.1,main = 'A')
abline(v = Scrit,col = 'azure4',lty=2)
plot(sort(ref_obj$scores),ref_obj$min_abundance_over_the_sorted,xlim = c(1,1.5),xlab = TeX("S_{crit}"),ylab = "Min. counts in reference",ylim = c(0,300),pch=20,cex=0.1,main = 'B')
abline(v = Scrit,col = 'azure4',lty=2)
dev.off()
par(mfrow = c(1,1))

pdf('../../Results/Gut_Scrit_dist.pdf',height = 3.5)
hist(ref_obj$scores,breaks = 50,xlab = 'Reference selection scores')
dev.off()

# run DACOMP tests
res_perm_Wilcoxon = dacomp.test(X = X,y = Y,ind_reference_taxa = Selected_references,
                                verbose = T,
                                q = Q_LEVEL,
                                compute_ratio_normalization = T,
                                test = DACOMP.TEST.NAME.WILCOXON,disable_DSFDR = T,
                                nr_perm = NR.PERMS)

res_perm_Welch= dacomp.test(X = X,y = Y,ind_reference_taxa = Selected_references,
                            verbose = T,
                            q = Q_LEVEL,
                            compute_ratio_normalization = T,
                            test = DACOMP.TEST.NAME.WELCH_LOGSCALE,disable_DSFDR = T,
                            nr_perm = NR.PERMS)

temp_obj = list(ref_obj = ref_obj,
                res_perm_Wilcoxon = res_perm_Wilcoxon,
                res_perm_Welch = res_perm_Welch,
                Scrit = Scrit,
                Selected_references = Selected_references)
save(temp_obj,file = filepath)


#W-CSS
X_CSS = t(metagenomeSeq::cumNormMat(t(X)))
res_Wilcoxon_Paulson = wilcoxon_taxa_wise(X_CSS,Y,normalize = F,normalize.P = 1,nr.perms = NR.PERMS)

#look at results:
hist(res_Wilcoxon_Paulson$p.values) #dist of P-values
disc_Wilcoxon_Paulson = which(p.adjust(res_Wilcoxon_Paulson$p.values,method = 'BH')<=Q_LEVEL) #which taxa are discovered as diff. abun.
length(disc_Wilcoxon_Paulson) #how many discoveries

#W-TSS, with similar look at results:
res_Wilcoxon_Percent = wilcoxon_taxa_wise(X_corrected,y = Y,normalize = T,normalize.P = 1,nr.perms = NR.PERMS)
hist(res_Wilcoxon_Percent$p.values)
disc_Wilcoxon_percent = which(p.adjust(res_Wilcoxon_Percent$p.values,method = 'BH')<=Q_LEVEL)
length(disc_Wilcoxon_percent) 

# check we have no zeros so aldex 2 can run on original indices
if(any(as.numeric(apply((X > 0),2,sum))==0)){
  stop('zeros in some columns, cannot run aldex 2!')
}
  

library(ALDEx2)
aldex.res.iqlr <- aldex(t(X), as.character(Y), mc.samples=128, denom="iqlr",
                        test="t", effect=FALSE,verbose = T)
aldex.res.zero <- aldex(t(X), as.character(Y), mc.samples=128, denom="zero",
                        test="t", effect=FALSE,verbose = T)
#run Wrench

library(Wrench)
W <- wrench( t(X), condition=Y  )
compositionalFactors <- W$ccf
normalizationFactors <- W$nf

library(DESeq2)
deseq_counts = t(X)
deseq_cols_data = DataFrame(Y = factor(Y))
deseq.obj <- DESeq2::DESeqDataSetFromMatrix(countData = deseq_counts,
                                            deseq_cols_data,
                                            ~ Y )
sizeFactors(deseq.obj) <- normalizationFactors
deseq2_res = DESeq2::DESeq(deseq.obj)
deseq2_res2 = results(deseq2_res)


#run ZINB-WAVE with deseq2
set.seed(1)
library(zinbwave)
library(phyloseq)
library(DESeq2)
source('zinbwave_imported_functions.R')


# Pack as phyloseq object
reads_counts_no_names = X
rownames(reads_counts_no_names) = NULL
physeq = phyloseq(otu_table(t(reads_counts_no_names),taxa_are_rows = TRUE),
                  sample_data(
                    data.frame(
                      grp = factor(Y),
                      col.names = rownames(X) 
                    )
                  )
)

#convert to DESEQ2 object and run methods as the eval_functions.R script in in https://github.com/mcalgaro93/sc2meta
physeq <- normDESeq2(physeq = physeq) 

epsilon = 1e10
zinbmodel <- zinbFit(Y = physeq@otu_table@.Data, 
                     X = model.matrix(~ physeq@sam_data$grp), K = 0,
                     epsilon = epsilon, commondispersion = TRUE, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())

weights <- computeExactWeights(model = zinbmodel,x = physeq@otu_table@.Data)
colnames(weights) <- colnames(physeq@otu_table)
rownames(weights) <- rownames(physeq@otu_table)

library(BiocParallel)
DESeq2_poscounts_zinbwave <- negBinTestDESeq2_zinbweights(physeq, normFacts = "poscounts",weights = weights)

#check the rejected and collect statistics
zinb_wave_deseq2_rejected = which(p.adjust(DESeq2_poscounts_zinbwave$pValMat[,1],method = 'BH')<=Q_LEVEL)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Combine results
#%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Combine to get results, by parameter value - compare different parameter values to different methods
res_matrix = data.frame(NA,row.names = 'Results')
res_matrix$nr_rejected_W = rep(NA,nrow(res_matrix)) #number of rejections - DACOMP with Wilcoxon (rarefaction- the default)
res_matrix$nr_rejected_W_dsfdr = rep(NA,nrow(res_matrix)) #number of rejections - DACOMP with Wilcoxon - DSFDR adjustment for multiplicity
res_matrix$nr_rejected_W_ratio = rep(NA,nrow(res_matrix)) #number of rejections - DACOMP with Wilcoxon - normalization by ratio
res_matrix$nr_rejected_t = rep(NA,nrow(res_matrix)) #number of rejections for DACOMP with a welch t-test
res_matrix$nr_rejected_t_ratio = rep(NA,nrow(res_matrix)) #number of rejections for DACOMP with a welch t-test and normalization by ratio
res_matrix$median_lambda = rep(NA,nrow(res_matrix)) #median lambda across tests
res_matrix$mean_lambda = rep(NA,nrow(res_matrix)) #mean lambda across tests
res_matrix$ANCOM_3_intersect = rep(NA,nrow(res_matrix)) #intersection of rejections with ANCOM (default parameter for ANCOM)
res_matrix$cell_intersect = rep(NA,nrow(res_matrix)) #intersection with W-FLOW
res_matrix$Wrench_intersect = rep(NA,nrow(res_matrix)) #intersection with Wrench
res_matrix$ALDEx2_W_intersect = rep(NA,nrow(res_matrix)) #intersection with ALDEx2 - two variantes
res_matrix$ALDEx2_t_intersect = rep(NA,nrow(res_matrix)) 
res_matrix$ref_size = rep(NA,nrow(res_matrix)) #reference set size (in number of taxa) for parameter configuration
res_matrix$ZINB_WAVE_DESEQ2_INTERSECT = rep(NA,nrow(res_matrix)) #intersection with ZINB-Wave+DESEQ2

res_Wilcoxon_list = list(res_perm_Wilcoxon)
res_Welch_list = list(res_perm_Welch)
ref_list = list(ref_obj)
# fill out data structure
for(i in 1:nrow(res_matrix)){
  res_matrix$nr_rejected_W[i] = length(which(p.adjust(res_Wilcoxon_list[[i]]$p.values.test,method = 'BH')<=Q_LEVEL))
  res_matrix$nr_rejected_W_dsfdr[i] = length(res_Wilcoxon_list[[i]]$dsfdr_rejected)
  res_matrix$median_lambda[i] = median(res_Wilcoxon_list[[i]]$lambda,na.rm = T)
  res_matrix$mean_lambda[i] = mean(res_Wilcoxon_list[[i]]$lambda,na.rm = T)
  res_matrix$nr_rejected_W_ratio[i] = length(which(p.adjust(res_Wilcoxon_list[[i]]$p.values.test.ratio.normalization,method = 'BH')<=Q_LEVEL))

  res_matrix$nr_rejected_t[i] = length(which(p.adjust(res_Welch_list[[i]]$p.values.test ,method = 'BH')<=Q_LEVEL))
  res_matrix$nr_rejected_t_ratio[i] = length(which(p.adjust(res_Welch_list[[i]]$p.values.test.ratio.normalization,method = 'BH')<=Q_LEVEL))  
  
  WCOMP_rejections = which(p.adjust(res_Wilcoxon_list[[i]]$p.values.test,method = 'BH')<=Q_LEVEL)
  disc_ANCOM_3 =which(colnames(ANCOM_otu_dat) %in% ANCOM_default_res$detected)
  res_matrix$ANCOM_3_intersect[i] = sum(WCOMP_rejections %in% disc_ANCOM_3)
  res_matrix$cell_intersect[i] = sum(WCOMP_rejections %in% disc_Wilcoxon_corrected)
  res_matrix$ref_size[i] = length(ref_list[[i]]$selected_references) 
  
  res_matrix$Wrench_intersect[i] = sum(WCOMP_rejections%in% which(deseq2_res2$padj<=Q_LEVEL))
  res_matrix$ALDEx2_W_intersect[i] = sum(WCOMP_rejections%in% which(aldex.res.iqlr$wi.eBH<=Q_LEVEL))
  res_matrix$ALDEx2_t_intersect[i] = sum(WCOMP_rejections%in% which(aldex.res.iqlr$we.eBH<=Q_LEVEL))
  res_matrix$ZINB_WAVE_DESEQ2_INTERSECT[i] = sum(WCOMP_rejections%in% zinb_wave_deseq2_rejected)
}

#save results
write.csv(res_matrix,file = '../../Results/Gut_Data_Results.csv',row.names =F)



# a list, containing the vecotrs of indices of discoveries, by method
disc_list = list(
disc_vec_ANCOM = which(colnames(ANCOM_otu_dat)%in% ANCOM_default_res$detected),
disc_Wilcoxon_corrected = disc_Wilcoxon_corrected ,
disc_Wilcoxon_Paulson = disc_Wilcoxon_Paulson,
disc_Wilcoxon_percent = disc_Wilcoxon_percent,
disc_W_COMP = which(p.adjust(res_Wilcoxon_list[[1]]$p.values.test,method = 'BH')<=Q_LEVEL),
disc_ALDEx2_iqlr_Wi = which(aldex.res.iqlr$wi.eBH <= Q_LEVEL),
disc_ALDEx2_iqlr_We = which(aldex.res.iqlr$we.eBH <= Q_LEVEL),
disc_Wrench = which(deseq2_res2$padj <= Q_LEVEL),
disc_W_COMP_Ratio = which(p.adjust(res_Wilcoxon_list[[1]]$p.values.test.ratio.normalization,method = 'BH')<=Q_LEVEL),
disc_ZINBWAVE_DESEQ2 = zinb_wave_deseq2_rejected
)

#compute a matrix of shared discoveries. Diagonal entries are the number of discoveries of each method
method_names = c('ANCOM','WIL-FLOW','WIL-CSS','WIL-TSS','DACOMP','ALDEx2-W','ALDEx2-t','Wrench','DACOMP-ratio','ZINBWAVE-DESEQ2')
shared_disc_mat = matrix(NA,nrow = length(method_names),ncol = length(method_names))
rownames(shared_disc_mat) = method_names
colnames(shared_disc_mat) = method_names
for(i in 1:length(method_names)){
  for(j in i:length(method_names)){
    shared_disc_mat[i,j] = sum(disc_list[[i]] %in% disc_list[[j]])
  }
}

#write to file
write.csv(shared_disc_mat,file = '../../Results/gut_cell_shared_disc_mat.csv')
focus_ind = c(1,2,3,7,5,9,10)
write.csv(shared_disc_mat[focus_ind,focus_ind],file = '../../Results/gut_cell_shared_disc_mat_for_paper.csv')


#Plot, compare large taxa, small taxa, and different subsets of methods
ind_to_prevalent = which(apply(X,2,mean)>=10)


disc_list_plot = disc_list[1:8]
names(disc_list_plot) = c('ANCOM','W-FLOW','W-CSS','W-TSS','DACOMP','ALDEx2-W','ALDEx2-t','Wrench')
disc_list_plot_prevalent = disc_list_plot
disc_list_plot_rare = disc_list_plot

for(i in 1:length(disc_list_plot_prevalent)){
  disc_list_plot_prevalent[[i]] = intersect(disc_list_plot_prevalent[[i]],ind_to_prevalent) 
  disc_list_plot_rare[[i]] = setdiff(disc_list_plot_rare[[i]],ind_to_prevalent)
}

temp = disc_list_plot[1:3]

library(gplots)
pdf(file = '../../Results/Only_Compositional.pdf',height = 6,width = 7)
par(mfrow=c(1,1))
venn(disc_list_plot[-c(2,3,4)])
par(mfrow=c(1,1))
dev.off()


pdf(file = '../../Results/Crohn_shared.pdf',height = 6,width = 14)
par(mfrow=c(1,2))
venn(disc_list_plot_prevalent[c(1,2,3,5,7)],simplify =T,small = 0.9)
venn(disc_list_plot_rare[c(1,2,3,5,7)],simplify = T,small = 0.9)
par(mfrow=c(1,1))
dev.off()

pdf(file = '../../Results/Crohn_shared_prevalent.pdf',height = 6,width = 7)
par(mar=c(0,0,0,0))
venn(disc_list_plot_prevalent[c(1,2,3,5,7)],small = 0.85)
dev.off()

pdf(file = '../../Results/Crohn_shared_rare.pdf',height = 6,width = 7)
par(mar=c(0,0,0,0))
venn(disc_list_plot_rare[c(1,2,3,5,7)],small = 0.85)
dev.off()

pdf(file = '../../Results/One_From_Each_FrameWork.pdf',height = 6,width = 7)
par(mfrow=c(1,1))
venn(disc_list_plot[-c(3,4,7)])
par(mfrow=c(1,1))
dev.off()

is_rejected = matrix(0, nrow = length(names(disc_list_plot)), ncol= ncol(X))
colnames(is_rejected) = 1:ncol(X)
rownames(is_rejected) = names(disc_list_plot)

for(i in 1:nrow(is_rejected)){
  is_rejected[i,disc_list_plot[[i]]] = 1  
}

is_rejected = is_rejected[,(which(apply(is_rejected,2,sum)>0))]
ord = order(apply(is_rejected,2,sum),decreasing = T)
is_rejected = is_rejected[,ord]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#output discoveries to file
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

dt_taxonomy = (unname(taxa))
colnames(dt_taxonomy) = paste0('L',1:7)
discoveries_matrix_for_output = matrix(F,nrow = nrow(dt_taxonomy),ncol = length(method_names))
colnames(discoveries_matrix_for_output) = method_names
for(i in 1:length(method_names)){
  discoveries_matrix_for_output[disc_list[[i]],i] = T
}

output_univariate_discoveries = cbind(seqs,dt_taxonomy,discoveries_matrix_for_output)
output_univariate_discoveries = as.data.frame(output_univariate_discoveries)

effect_string = rep(NA,nrow(output_univariate_discoveries))
effect_direction_string = rep(NA,nrow(output_univariate_discoveries))
sum_ref = apply(X[,Selected_references],1,sum)
for(taxon_id in 1:nrow(output_univariate_discoveries)){
  #taxon_id = 1
  sum_ref_for_compute = sum_ref
  if(taxon_id %in% Selected_references){
    sum_ref_for_compute = sum_ref_for_compute - X[,taxon_id]
  }
  ratio = X[,taxon_id]/(X[,taxon_id]+sum_ref_for_compute)
  ratio_CD = mean(ratio[Y==1])
  ratio_H = mean(ratio[Y==0])
  if(ratio_CD>=ratio_H){
    effect_direction_string[taxon_id] = 'CD>=H'
    effect_string[taxon_id] = paste0('CD>=H;mean ratio CD : ',(ratio_CD),' ;mean ratio H : ',(ratio_H))
  }else{
    effect_direction_string[taxon_id] = 'H>CD'
    effect_string[taxon_id] = paste0('H>CD;mean ratio CD : ',(ratio_CD),' ;mean ratio H : ',(ratio_H))
  }
}
output_univariate_discoveries$effect_direction = effect_direction_string
output_univariate_discoveries$effect = effect_string
write.csv(output_univariate_discoveries,file = '../../Results/univariate_discoveries.csv',quote = F,row.names = F)
