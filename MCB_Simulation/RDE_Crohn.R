# Script for processing the crohn dataset, from Van De Putte et al. (2017)
# The data analyzed is the data from the third analysis done in the paper.

memory.limit(25000) # scale memory
library(biomformat)
Q_LEVEL = 0.1 #BH threshold level

file_path <- "./gut_otu_table.RData" #load counts data
load(file_path)
Average_Cell_Count = read.csv(file = 'FlowCytometryCounts.csv')$Average_Cell_Count #load flow cytometry counts
Y = rep(0,95) #group labeling is such that the first 29 subjects are CD, the rest are healthy.
Y[1:29] = 1

#show a box plot for the distribution of average cell counts across subjects
boxplot(Average_Cell_Count[Y==1],Average_Cell_Count[Y==0],names = c('Sick','Healthy'),main='Boxplot for cell counts in Healthy/Sick')
mean(Average_Cell_Count[Y==0])/(mean(Average_Cell_Count[Y==1])) #ratio between mean number of reads - discussed in paper
#2.508608
median(Average_Cell_Count[Y==0])/(median(Average_Cell_Count[Y==1])) # ratio between median number of reads
#3.086826

# order OTU table (from deblur) by original sample names
otu_table = t(otu_table)
reordering_permutation = order(rownames(otu_table))
otu_table = otu_table[reordering_permutation,]
rownames(otu_table) #it is now ordered

dim(otu_table)
X = otu_table

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
median_SD_thres_Vec = seq(1.2,1.4,0.1) 
PS_value = 1 # Pseudocount value used
parameter_matrix = expand.grid(median_SD_thres_Vec = median_SD_thres_Vec)

current_selected_ref_obj = NULL

ref_list = list() #used to store referecne objects
res_Wilcoxon_list = list() #used to store results with Wilcoxon rank sum tests
res_Welch_list = list() # used to store results with Welch t-tests
set.seed(1)

#stored file path
filepath_i = function(i){return(paste0('../../Results/Gut_temp_v2_file_',i,'.RData'))}

#iterate over parameters
for(i in 1:nrow(parameter_matrix)){
  set.seed(1) # for reproducability of permutations
  current_thres = parameter_matrix$median_SD_thres_Vec[i]
  
  # select references
  print(paste0('Selecting references , thres = ',current_thres))
  current_selected_ref_obj = dacomp.select_references(X = X,median_SD_threshold = current_thres,maximal_TA = 200,Pseudo_Count_used = PS_value,  verbose = F)
  

  ref_list[[i]] =   current_selected_ref_obj
  
  # run DACOMP tests
  res_perm_Wilcoxon = dacomp.test(X = X,y = Y,ind_reference_taxa = current_selected_ref_obj$selected_references,
                                  verbose = F,
                                  q = Q_LEVEL,
                                  compute_ratio_normalization = T,
                                  test = DACOMP.TEST.NAME.WILCOXON,
                                                 nr_perm = NR.PERMS)
  
  res_perm_Welch= dacomp.test(X = X,y = Y,ind_reference_taxa = current_selected_ref_obj$selected_references,
                                  verbose = F,
                                  q = Q_LEVEL,
                                  compute_ratio_normalization = T,
                                  test = DACOMP.TEST.NAME.WELCH_LOGSCALE,
                                  nr_perm = NR.PERMS)
  
  res_Wilcoxon_list[[i]] = res_perm_Wilcoxon
  res_Welch_list[[i]] = res_perm_Welch
  temp_obj = list(current_selected_ref_obj = current_selected_ref_obj,
                  res_perm_Wilcoxon = res_perm_Wilcoxon,
                  res_perm_Welch = res_perm_Welch)
  current_temp_filepath_i = filepath_i(i)
  save(temp_obj,file = current_temp_filepath_i)
}


set.seed(1)
#load
ref_list = list()
res_Wilcoxon_list = list()
res_Welch_list = list()
for(i in 1:nrow(parameter_matrix)){
  print(paste0('Loading case ',i))
  load(filepath_i(i))
  ref_list[[i]] = temp_obj$current_selected_ref_obj
  res_Wilcoxon_list[[i]] = temp_obj$res_perm_Wilcoxon
  res_Welch_list[[i]] = temp_obj$res_perm_Welch
}


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
deseq.obj <- DESeq2::DESeqDataSetFromMatrix(countData = t(X),
                                            DataFrame(Y = factor(Y)),
                                            ~ Y )
sizeFactors(deseq.obj) <- normalizationFactors
deseq2_res = DESeq2::DESeq(deseq.obj)
deseq2_res2 = results(deseq2_res)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Combine results
#%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Combine to get results, by parameter value - compare different parameter values to different methods
res_matrix = parameter_matrix # each row is a parameter configuration
res_matrix = as.data.frame(res_matrix)
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
}

#save results
write.csv(res_matrix,file = '../../Results/Gut_Data_Results.csv',row.names =F)




disc_list = list(
disc_vec_ANCOM = which(colnames(ANCOM_otu_dat)%in% ANCOM_default_res$detected),
disc_Wilcoxon_corrected = disc_Wilcoxon_corrected ,
disc_Wilcoxon_Paulson = disc_Wilcoxon_Paulson,
disc_Wilcoxon_percent = disc_Wilcoxon_percent,
disc_W_COMP = which(p.adjust(res_Wilcoxon_list[[which(median_SD_thres_Vec == 1.3)]]$p.values.test,method = 'BH')<=Q_LEVEL),
disc_ALDEx2_iqlr_Wi = which(aldex.res.iqlr$wi.eBH <= Q_LEVEL),
disc_ALDEx2_iqlr_We = which(aldex.res.iqlr$we.eBH <= Q_LEVEL),
disc_Wrench = which(deseq2_res2$padj <= Q_LEVEL)
)

method_names = c('ANCOM','WIL-FLOW','WIL-CSS','WIL-TSS','DACOMP','ALDEx2-W','ALDEx2-t','Wrench')
shared_disc_mat = matrix(NA,nrow = length(method_names),ncol = length(method_names))
rownames(shared_disc_mat) = method_names
colnames(shared_disc_mat) = method_names
for(i in 1:length(method_names)){
  for(j in i:length(method_names)){
    shared_disc_mat[i,j] = sum(disc_list[[i]] %in% disc_list[[j]])
  }
}

write.csv(shared_disc_mat,file = '../../Results/gut_cell_shared_disc_mat.csv')


#Plot, compare large taxa

ind_to_prevalent = which(apply(X,2,mean)>=10)


disc_list_plot = disc_list
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
venn(disc_list_plot_prevalent[c(1,2,3,5,7)])
venn(disc_list_plot_rare[c(1,2,3,5,7)])
par(mfrow=c(1,1))
dev.off()

pdf(file = '../../Results/Crohn_shared_prevalent.pdf',height = 6,width = 7)
venn(disc_list_plot_prevalent[c(1,2,3,5,7)])
dev.off()

pdf(file = '../../Results/Crohn_shared_rare.pdf',height = 6,width = 7)
venn(disc_list_plot_rare[c(1,2,3,5,7)])
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

dim(is_rejected)
#image(t(is_rejected))

library(reshape2)
melted_is_rejected = melt(data.frame(Method=rownames(is_rejected), is_rejected), id.vars="Method")


library(yarrr)
pallete = yarrr::piratepal("basel",
                 plot.result = FALSE,
                 trans = 0)[-c(8)]          # Slightly transparent

op <- par(bg = "#f4f5f7")

pdf(file = '../../Results/Shared_Plot.pdf',height = 6,width = 7)
plot(c(-80, ncol(is_rejected)), c(0, 8), type = "n", xlab = "", ylab = "",
     main = "")
reorder_vec = c(3,1,4,2,5,6,7,8)
for(i in 1:8){
  p = which(is_rejected[reorder_vec[i],]==1)
  rect(xleft = p-1,xright = p,ybottom = i-1,ytop = i,col = pallete[i],border = NA)  #rainbow(n=8)
  text(x=-50,y= i-0.5,labels = names(disc_list_plot)[reorder_vec[i]])
}
dev.off()


