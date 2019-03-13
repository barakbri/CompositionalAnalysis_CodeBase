#graphs for paper

library(ggplot2)

dt = read.csv('C:/MCB2/Results/REFSIM_Combined_Results.csv')

#table 1
dt_1 = dt[dt$setting_id %in% c(30:40),]
dt_1$m1 = rep(NA,nrow(dt_1))
dt_1$effect = rep(NA,nrow(dt_1))

dt_1$m1[dt_1$setting_id %in% c(31,33,35,37,39) ] = 10
dt_1$m1[dt_1$setting_id%in% (c(31,33,35,37,39)+1) ] = 20
dt_1$m1[dt_1$setting_id == 30] = 0

dt_1$effect[dt_1$setting_id == 30] = 0
dt_1$effect[dt_1$setting_id %in% c(31,32)] = 0.5
dt_1$effect[dt_1$setting_id %in% c(33,34)] = 1.0
dt_1$effect[dt_1$setting_id %in% c(35,36)] = 1.5
dt_1$effect[dt_1$setting_id %in% c(37,38)] = 2.0
dt_1$effect[dt_1$setting_id %in% c(39,40)] = 2.5

#install.packages('reshape2')
library(reshape2)
names(dt_1)
method_names = as.character(unique(dt_1$methodlabel))
dt_tp = matrix(NA,nrow = length(30:40),ncol  = length(method_names))
dt_rejected = matrix(NA,nrow = length(30:40),ncol  = length(method_names))
dt_FDR = matrix(NA,nrow = length(30:40),ncol  = length(method_names))
for(i in 1:length(method_names)){
  for(j in 1:11){
    tp_row = which(dt_1$setting_id == j+29 & dt_1$methodlabel == method_names[i])
    dt_tp[j,i] = dt_1$tp[tp_row[1]]
    dt_rejected[j,i] = dt_1$rejected[tp_row[1]]
    dt_FDR[j,i] = dt_1$fdr[tp_row[1]]
  }
}

#View(dt_1)

#ggplot(dt_1)+geom_point(aes(x = effect,y = tp,color = m1,shape = methodlabel))
dt_combined = round(dt_FDR,2)
# for(i in 1:ncol(dt_combined)){
#   dt_combined[,i] = as.character(round(as.numeric(dt_combined[,i]),0))
# }
# for(i in 1:nrow(dt_combined)){
#   for(j in 1:ncol(dt_combined)){
#     dt_combined[i,j] = paste0(dt_combined[i,j],', ',round(dt_FDR[i,j],2))
#   }
# }

colnames(dt_combined) = method_names
#dt_combined = dt_combined[,-7]
dt_combined = as.data.frame(dt_combined)
dt_combined$m1 =c(0,rep(c(10,20),5))
dt_combined$effect =c(0,rep(0.5,2),rep(1.0,2),rep(1.5,2),rep(2.0,2),rep(2.5,2))

library(xtable)
names(dt_combined)
cols_to_keep = c(16,17,1,13,14,15,9,6)
dt_combined = dt_combined[,cols_to_keep]
names(dt_combined)
xtable::xtable(dt_combined)


print(xtable(dt_combined), include.rownames=FALSE)


#table 2

dt_2 = dt[dt$setting_id %in% c(66:75),]
dt_2$m1 = rep(NA,nrow(dt_2))
dt_2$p_high = rep(NA,nrow(dt_2))

dt_2$m1[dt_2$setting_id %in% c(66:70) ] = 120
dt_2$m1[dt_2$setting_id %in% c(71:75)] = 60



dt_2$p_high[dt_2$setting_id %in% c(66,71)] = 0.9
dt_2$p_high[dt_2$setting_id %in% c(67,72)] = 0.8
dt_2$p_high[dt_2$setting_id %in% c(68,73)] = 0.7
dt_2$p_high[dt_2$setting_id %in% c(69,74)] = 0.6
dt_2$p_high[dt_2$setting_id %in% c(70,75)] = 0.5

#install.packages('reshape2')
library(reshape2)
names(dt_2)
method_names = as.character(unique(dt_2$methodlabel))
dt_tp = matrix(NA,nrow = length(66:75),ncol  = length(method_names))
dt_rejected = matrix(NA,nrow = length(66:75),ncol  = length(method_names))
dt_FDR = matrix(NA,nrow = length(66:75),ncol  = length(method_names))
for(i in 1:length(method_names)){
  for(j in 1:10){
    tp_row = which(dt_2$setting_id == j+65 & dt_2$methodlabel == method_names[i])
    dt_tp[j,i] = dt_2$tp[tp_row[1]]
    dt_rejected[j,i] = dt_2$rejected[tp_row[1]]
    dt_FDR[j,i] = dt_2$fdr[tp_row[1]]
  }
}

#View(dt_1)

#ggplot(dt_1)+geom_point(aes(x = effect,y = tp,color = m1,shape = methodlabel))
dt_combined = dt_tp
for(i in 1:ncol(dt_combined)){
  dt_combined[,i] = as.character(round(as.numeric(dt_combined[,i]),0))
}
for(i in 1:nrow(dt_combined)){
  for(j in 1:ncol(dt_combined)){
    dt_combined[i,j] = paste0(dt_combined[i,j],', ',round(dt_FDR[i,j],2))
  }
}

colnames(dt_combined) = method_names
dt_combined = dt_combined[,-7]
dt_combined = as.data.frame(dt_combined)
dt_combined$m1 =c(rep(120,5),rep(60,5))
dt_combined$p_high =rep(c(0.9,0.8,0.7,0.6,0.5),2)

names(dt_combined)
cols_to_keep = c(14,15,1,12,13,8,6)
dt_combined = dt_combined[,cols_to_keep]
names(dt_combined)
xtable::xtable(dt_combined)


print(xtable(dt_combined), include.rownames=FALSE)
#####

#table 3

dt_3 = dt[dt$setting_id %in% c(76:79),]
#dt_3$n = c('15:15','20:20','25:25','30:30')



#install.packages('reshape2')
library(reshape2)

method_names = as.character(unique(dt_3$methodlabel))
dt_tp = matrix(NA,nrow = length(76:79),ncol  = length(method_names))
dt_rejected = matrix(NA,nrow = length(76:79),ncol  = length(method_names))
dt_FDR = matrix(NA,nrow = length(76:79),ncol  = length(method_names))
for(i in 1:length(method_names)){
  for(j in 1:4){
    tp_row = which(dt_3$setting_id == j+75 & dt_3$methodlabel == method_names[i])
    dt_tp[j,i] = dt_3$tp[tp_row[1]]
    dt_rejected[j,i] = dt_3$rejected[tp_row[1]]
    dt_FDR[j,i] = dt_3$fdr[tp_row[1]]
  }
}

#View(dt_1)

dt_tp = cbind(c('15:15','20:20','25:25','30:30'),round(dt_tp,2))
dt_FDR = cbind(c('15:15','20:20','25:25','30:30'),round(dt_FDR,2))
colnames(dt_tp) = c('n',method_names)
colnames(dt_FDR) = c('n',method_names)

cols_to_keep = c(1,14,15,2,10,7)
dt_tp = dt_tp[,cols_to_keep]
dt_FDR = dt_FDR[,cols_to_keep]

print(xtable(dt_tp[,-c(3,6)]), include.rownames=FALSE)
print(xtable(dt_FDR), include.rownames=FALSE)





