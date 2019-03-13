# Script 

library(biomformat)
file_path <- "C:/MCB2/MCB2/MCB_Simulation/gut_otu_table.RData"
load(file_path)
Q_LEVEL = 0.1

#dim(otu_table)

Average_Cell_Count = c(5314887324,
61042558217,
8409441603,
32494562535,
18764726648,
12110411381,
14660465393,
94291923126,
11438369157,
7208134178,
17961581253,
112399396469,
48073771982,
48871825000,
57857074420,
16078707428,
49526146305,
15501173550,
25077773966,
83992769634,
17572221569,
112523183751,
81790788530,
37664569474,
60352099591,
98691161158,
68314356782,
12194333170,
57346163417,
99183906832,
217221266674,
151329659639,
90299128192,
95809600876,
132035763643,
39495621712,
119271846132,
138185802083,
121428187249,
145939517130,
224754815252,
72756226469,
35927546063,
107628776723,
130776743130,
41841212649,
22764340926,
172987119956,
46425487958,
153062883269,
95590999034,
208495566551,
53070915155,
18359008555,
70503715565,
58500810719,
89307497088,
147599455753,
19187237081,
217276666820,
253525709584,
25395316913,
22507873775,
142169041331,
116275992152,
142346591541,
161477761409,
5853245978,
47942791702,
28769973004,
168930811359,
126572467016,
117499264965,
116251948461,
177107799526,
74197893718,
48330742911,
203676977127,
187211301217,
174670957469,
75065628771,
54160688666,
165005851861,
215345198361,
150051113739,
191391240184,
20159773225,
38956568917,
123582072411,
98056602368,
184826494881,
54837809766,
70590688175,
132689255189,
98343121358)

Y = rep(0,95)
Y[1:29] = 1

boxplot(Average_Cell_Count[Y==1],Average_Cell_Count[Y==0],names = c('Sick','Healthy'),main='Boxplot for qPCR counts in Healthy/Sick')
mean(Average_Cell_Count[Y==0])/(mean(Average_Cell_Count[Y==1]))
#2.508608
median(Average_Cell_Count[Y==0])/(median(Average_Cell_Count[Y==1]))
#3.086826

otu_table = t(otu_table)
reordering_permutation = order(rownames(otu_table))
otu_table = otu_table[reordering_permutation,]
rownames(otu_table) #it is now ordered

dim(otu_table)
X = otu_table


median(sort(as.numeric(apply(X,1,sum))))
#[1] 20437
mean(sort(as.numeric(apply(X,1,sum))))
#[1] 22886

prevalence_matrix = 1*(X>0)
prev_taxa = as.numeric(apply(prevalence_matrix,2,mean))
hist(prev_taxa)

X_as_percentage = X
for(i in 1:nrow(X_as_percentage)){
  X_as_percentage[i,] = (X_as_percentage[i,])/sum(X_as_percentage[i,])
}
mean_percent_taxa  = as.numeric(apply(X_as_percentage,2,mean))
#correct by volume counts

library(phyloseq)

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
    x <- rarefy_even_depth(cnv_corrected_abundance_table_phyloseq[i,], sample.size = rarefy_to[i], rngseed = 711, replace = FALSE, trimOTUs = F, verbose = FALSE)
    rarefied_matrix[i,] = x
  }
  normalised_rarefied_matrix = rarefied_matrix/rowSums(rarefied_matrix)
  QMP = normalised_rarefied_matrix*cell_counts_table[1,]
  return(QMP)
}



# Example
#a = matrix( c(4,4,2,1,8,5,2,0,3,5,3,1,10,8,3,0,0,6,4,3), nrow=5, ncol=4, byrow = TRUE, dimnames = list(c("Sample A", "Sample B", "Sample C", "Sample D", "Sample E"),c("taxa1", "taxa2", "taxa3", "taxa4"))) # my cnv_corrected_abundance_table
#b = matrix(c(10,20,34,21,12), nrow=5, ncol=1, byrow = TRUE, dimnames = list(c("Sample A", "Sample B", "Sample C", "Sample D", "Sample E"),c("#")))*100000 # my cell_counts_table
#rarefy_even_sampling_depth(a,b)

rownames(X)
Average_Cell_Count_matrix = matrix(Average_Cell_Count,ncol = 1)
rownames(Average_Cell_Count_matrix) = rownames(X)
X_corrected = rarefy_even_sampling_depth(X,Average_Cell_Count_matrix)

# X_corrected = X
# for(i in 1:nrow(X_corrected)){
#   X_corrected[i,] = X_corrected[i,] * Average_Cell_Count[i]
# }
# image(t(log10(X_corrected+1)))

#remove taxa in less than 4 people
OTHER = 2 #4
to_remove = which(as.numeric(apply(prevalence_matrix, 2, sum))<OTHER)
X = X[,-to_remove]
X_corrected = X_corrected[,-to_remove]





#Wilcoxon over the corrected
source('C:/MCB2/MCB2/MCB_Simulation/Wilcoxon_TaxaWise.R')
res_Wilcoxon_corrected = wilcoxon_taxa_wise(X_corrected,y = Y)
hist(res_Wilcoxon_corrected$p.values)
total_counts_in_taxa = as.numeric(apply(prevalence_matrix,2,sum))

disc_Wilcoxon_corrected = which(p.adjust(res_Wilcoxon_corrected$p.values,method = 'BH')<=Q_LEVEL)
length(disc_Wilcoxon_corrected) #211 



library(ancom.R)
ANCOM_otu_dat = X
ANCOM_otu_dat = cbind(ANCOM_otu_dat,Y)
ANCOM_otu_dat = data.frame(ANCOM_otu_dat)
ANCOM_default_res = ancom.R::ANCOM(OTUdat = ANCOM_otu_dat,multcorr = 3,sig = Q_LEVEL)
#ANCOM_multcorr_2_res = ancom.R::ANCOM(OTUdat = ANCOM_otu_dat,multcorr = 2,sig = Q_LEVEL)
#ANCOM_multcorr_1_res = ancom.R::ANCOM(OTUdat = ANCOM_otu_dat,multcorr = 1,sig = Q_LEVEL)

intersect_with_ANCOM = function(ANCOM_res,subzero_rejections){
  print(paste0('NUMBER OF ANCOM REJ:',length(ANCOM_res$detected)))
  disc_ANCOM =which(colnames(ANCOM_otu_dat) %in% ANCOM_res$detected)
  print(paste0('ANCOM in WILCOX corrected:',sum(disc_ANCOM %in% disc_Wilcoxon_corrected)))
  print(paste0('SUBZERO in ANCOM corrected:',sum(subzero_rejections %in% disc_ANCOM)))
  
  possibly_novel_disc = subzero_rejections[which(!(subzero_rejections %in% disc_ANCOM))]
  print(paste0('Verified, Unique:',sum(possibly_novel_disc %in% disc_Wilcoxon_corrected),'/',length(possibly_novel_disc)))
  ANCOM_TO_INVESTIGATE = disc_ANCOM[which(!(disc_ANCOM %in% subzero_rejections ))]
  return(ANCOM_TO_INVESTIGATE)
}

#SAVE POINT
#load("E:/MCB2/Results/gut_V2_savepoint.RData")

#Checking subzero/thac0
library(subzero)
median_SD_thres_Vec = seq(1.0,1.4,0.1)#seq(1.0,1.4,0.05)
lambda_multiplier_Vec = c(1.0)
parameter_matrix = expand.grid(lambda_multiplier_Vec = lambda_multiplier_Vec ,
                               median_SD_thres_Vec = median_SD_thres_Vec)
source('C:/MCB2/MCB2/MCB_Simulation/SelectReferences_MedianSD_Threshold.R')

current_selected_ref_obj = NULL
last_computed_thres = NULL
ref_list = list()
res_Wilcoxon_list = list()
set.seed(1)

filepath_i = function(i){return(paste0('C:/MCB2/Results/Gut_temp_v2_file_',i,'.RData'))}

for(i in 1:nrow(parameter_matrix)){
  current_thres = parameter_matrix$median_SD_thres_Vec[i]
  current_multiplier = parameter_matrix$lambda_multiplier_Vec[i]
  
  NEED_TO_SELECT = T
  if(!is.null(last_computed_thres))
    if(last_computed_thres == current_thres)
      NEED_TO_SELECT = F
  
  if(NEED_TO_SELECT){
    print(paste0('Selecting references , thres = ',current_thres))
    current_selected_ref_obj = select.references.Median.SD.Threshold(X,median_SD_threshold = current_thres,maximal_TA = 200,  verbose = F)
    last_computed_thres = current_thres
  }
  ref_list[[i]] =   current_selected_ref_obj
  
  
  print(paste0('Running subzero , mult = ',current_multiplier))
  res_perm_Wilcoxon = subzero::subzero.dfdr.test(X = X,y = Y,
                                                 nr_reference_taxa = current_selected_ref_obj$selected_references,verbose = F,q = Q_LEVEL,
                                                 nr_rarefactions_multiple_X = 1,
                                                 nr_perm = 30000,lambda_multiplier = current_multiplier)
  res_Wilcoxon_list[[i]] = res_perm_Wilcoxon
  temp_obj = list(current_selected_ref_obj = current_selected_ref_obj,
                  res_perm_Wilcoxon = res_perm_Wilcoxon)
  current_temp_filepath_i = filepath_i(i)
  save(temp_obj,file = current_temp_filepath_i)
}
#load
ref_list = list()
res_Wilcoxon_list = list()
for(i in 1:nrow(parameter_matrix)){
  print(paste0('Loading case ',i))
  load(filepath_i(i))
  ref_list[[i]] = temp_obj$current_selected_ref_obj
  res_Wilcoxon_list[[i]] = temp_obj$res_perm_Wilcoxon
}

#Combine to get results
res_matrix = parameter_matrix
res_matrix = as.data.frame(res_matrix)
res_matrix$nr_rejected = rep(NA,nrow(res_matrix))
res_matrix$nr_rejected_dsfdr = rep(NA,nrow(res_matrix))
res_matrix$nr_rejected_ref = rep(NA,nrow(res_matrix))
res_matrix$nr_rejected_ref_dsfdr = rep(NA,nrow(res_matrix))
res_matrix$median_lambda = rep(NA,nrow(res_matrix))
res_matrix$mean_lambda = rep(NA,nrow(res_matrix))
res_matrix$ANCOM_3_intersect = rep(NA,nrow(res_matrix))
res_matrix$ANCOM_2_intersect = rep(NA,nrow(res_matrix))
res_matrix$ANCOM_1_intersect = rep(NA,nrow(res_matrix))
res_matrix$qPCR_intersect = rep(NA,nrow(res_matrix))
res_matrix$ref_size = rep(NA,nrow(res_matrix))


for(i in 1:nrow(res_matrix)){
  
  res_matrix$nr_rejected[i] = length(which(p.adjust(res_Wilcoxon_list[[i]]$p.values.test,method = 'BH')<=Q_LEVEL))
  res_matrix$nr_rejected_dsfdr[i] = length(res_Wilcoxon_list[[i]]$rejected)
  res_matrix$nr_rejected_ref[i] = length(which(p.adjust(res_Wilcoxon_list[[i]]$rejected_ref,method = 'BH')<=Q_LEVEL))
  res_matrix$nr_rejected_ref_dsfdr[i] = length(res_Wilcoxon_list[[i]]$rejected_ref)
  res_matrix$median_lambda[i] = median(res_Wilcoxon_list[[i]]$min_value_array[-ref_list[[i]]$selected_references])
  res_matrix$mean_lambda[i] = mean(res_Wilcoxon_list[[i]]$min_value_array[-ref_list[[i]]$selected_references])
  
  subzero_rejections = which(p.adjust(res_Wilcoxon_list[[i]]$p.values.test,method = 'BH')<=Q_LEVEL)
  disc_ANCOM_3 =which(colnames(ANCOM_otu_dat) %in% ANCOM_default_res$detected)
  disc_ANCOM_2 =-1#which(colnames(ANCOM_otu_dat) %in% ANCOM_multcorr_2_res$detected)
  disc_ANCOM_1 =-1#which(colnames(ANCOM_otu_dat) %in% ANCOM_multcorr_1_res$detected)
  res_matrix$ANCOM_3_intersect[i] = sum(subzero_rejections %in% disc_ANCOM_3)
  res_matrix$ANCOM_2_intersect[i] = -1#sum(subzero_rejections %in% disc_ANCOM_2)
  res_matrix$ANCOM_1_intersect[i] = -1#sum(subzero_rejections %in% disc_ANCOM_1)
  res_matrix$qPCR_intersect[i] = sum(subzero_rejections %in% disc_Wilcoxon_corrected)
  res_matrix$ref_size[i] = length(ref_list[[i]]$selected_references) 
}


write.csv(res_matrix,file = 'C:/MCB2/Results/Gut_Data_Results.csv',row.names =F)
#View(res_matrix)


length(ANCOM_default_res$detected)
#length(ANCOM_multcorr_2_res$detected)
#length(ANCOM_multcorr_1_res$detected)
length(disc_Wilcoxon_corrected)

X_CSS = t(metagenomeSeq::cumNormMat(t(X)))
res_Wilcoxon_Paulson = wilcoxon_taxa_wise(X_CSS,Y,normalize = F,normalize.P = 0.75)
hist(res_Wilcoxon_Paulson$p.values)
disc_Wilcoxon_Paulson = which(p.adjust(res_Wilcoxon_Paulson$p.values,method = 'BH')<=Q_LEVEL)
length(disc_Wilcoxon_Paulson)
sum(disc_Wilcoxon_Paulson %in% disc_Wilcoxon_corrected)

res_Wilcoxon_Percent = wilcoxon_taxa_wise(X_corrected,y = Y,normalize = T,normalize.P = 1)
hist(res_Wilcoxon_Percent$p.values)


disc_Wilcoxon_percent = which(p.adjust(res_Wilcoxon_Percent$p.values,method = 'BH')<=Q_LEVEL)
length(disc_Wilcoxon_percent) #211 

disc_list = list(
disc_vec_ANCOM = which(colnames(ANCOM_otu_dat)%in% ANCOM_default_res$detected),
disc_Wilcoxon_corrected = disc_Wilcoxon_corrected ,
disc_Wilcoxon_Paulson = disc_Wilcoxon_Paulson,
disc_Wilcoxon_percent = disc_Wilcoxon_percent,
RAR = which(p.adjust(res_Wilcoxon_list[[4]]$p.values.test,method = 'BH')<=Q_LEVEL)
)

method_names = c('ANCOM','WIL-FLOW','WIL-CSS','WIL-TSS','RAR')
shared_disc_mat = matrix(NA,nrow = length(method_names),ncol = length(method_names))
rownames(shared_disc_mat) = method_names
colnames(shared_disc_mat) = method_names
for(i in 1:length(method_names)){
  for(j in i:length(method_names)){
    shared_disc_mat[i,j] = sum(disc_list[[i]] %in% disc_list[[j]])
  }
}

write.csv(shared_disc_mat,file = 'C:/MCB2/Results/gut_qPCR_shared_disc_mat.csv')

# 
# > length(ANCOM_default_res$detected)
# [1] 216
# > length(ANCOM_multcorr_2_res$detected)
# [1] 198
# > length(ANCOM_multcorr_1_res$detected)
# [1] 121
# > length(disc_Wilcoxon_corrected)
# [1] 211
# > res_Wilcoxon_Paulson = wilcoxon_taxa_wise(X,Y,normalize = T,normalize.P = 0.75)
# > hist(res_Wilcoxon_Paulson$p.values)
# > disc_Wilcoxon_Paulson = which(p.adjust(res_Wilcoxon_Paulson$p.values,method = 'BH')<=Q_LEVEL)
# > length(disc_Wilcoxon_Paulson)
# [1] 283
# > sum(disc_Wilcoxon_Paulson %in% disc_Wilcoxon_corrected)
# [1] 171
# > 
# 
# > dim(X)
# [1]   95 1569
# > table(Y)
# Y
# 0  1 
# 66 29


stop('STOPPED!!')

source('C:/MCB2/MCB2/MCB_Simulation/SelectReferences_MedianSD_Threshold.R')
  current_selected_ref_obj = select.references.Median.SD.Threshold(X,
                                                                   median_SD_threshold = 1.2,
                                                                   maximal_TA = 200,  verbose = F)
  

# 
# set.seed(1)
# res_perm_Deconv = subzero::subzero.dfdr.test(X = X,y = Y,
#                                                nr_reference_taxa = current_selected_ref_obj$selected_references,verbose = T,
#                                                q = Q_LEVEL,
#                                                nr_rarefactions_multiple_X = 1,test = 'Log.Avg.Diff', # 'Deconv.CVM.2',
#                                                nr_perm = 10000,lambda_multiplier = 1)
# 

set.seed(1)
res_perm_Deconv_LAD = subzero::subzero.dfdr.test(X = X,y = Y,
                                             nr_reference_taxa = current_selected_ref_obj$selected_references,verbose = T,
                                             q = Q_LEVEL,
                                             nr_rarefactions_multiple_X = 1,test = 'Log.Avg.Diff',
                                             nr_perm = 30000,lambda_multiplier = 1)

length(res_perm_Deconv_LAD$rejected)

set.seed(1)
res_perm_Deconv_AD = subzero::subzero.dfdr.test(X = X,y = Y,
                                                 nr_reference_taxa = current_selected_ref_obj$selected_references,verbose = T,
                                                 q = Q_LEVEL,
                                                 nr_rarefactions_multiple_X = 1,test = 'Avg.Diff',
                                                 nr_perm = 30000,lambda_multiplier = 1)

length(res_perm_Deconv_AD$rejected)

res_perm_Deconv_KS$rejected
hist(res_perm_Deconv_KS$p.values.test)
sum(res_perm_Deconv_KS$p.values.test<= 0.1,na.rm = T)
length(which(p.adjust(res_perm_Deconv_KS$p.values.test,method = 'BH')<=Q_LEVEL))

res_perm_Deconv$rejected
hist(res_perm_Deconv$p.values.test)
length(which(p.adjust(res_perm_Deconv$p.values.test,method = 'BH')<=Q_LEVEL))


ind = 25
n_vec = apply(X,1,sum)
temp = deconv_stat(Counts = X[,ind],Labels = Y,TotalCounts = n_vec)
plot(temp$G_X)
points(temp$G_Y,col = 'red')
temp$stat.ks

temp2 = deconv_stat(Counts = X[,ind],Labels = sample(Y),TotalCounts = n_vec)
plot(temp2$G_X)
points(temp2$G_Y,col = 'red')
temp2$stat.ks


library(lattice)
lattice.options(default.theme = standard.theme(color = FALSE))
X_plot = X
rownames(X_plot) = NULL
colnames(X_plot) = NULL
lattice::levelplot(log10(t(X_plot)+1),aspect = 1,xlab = "OTUs",ylab = "Samples", col.regions = gray(0:100/100))
