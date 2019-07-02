#graphs for paper
library(ggplot2)

dt = read.csv('../../Results/REFSIM_Combined_Results.csv')
dt_sd = read.csv('../../Results/REFSIM_Combined_Results_sd.csv')

methods_to_keep = c('ALDEx2_Welch',
'ALDEx2_Wilcoxon',
'ANCOM',
'DACOMP,Welch,division,S = 1.3',
'DACOMP,Welch,rarefaction,S = 1.3',
'DACOMP,Wilcoxon,division,S = 1.3',
'DACOMP,Wilcoxon,rarefaction,S = 1.3',
'HG,S = 1.3, Oracle',
'WILCOXON_PAULSON',
'WILCOXON_PERCENT',
'WILCOXON_qPCR',
'Wrench')

dt = dt[dt$methodlabel %in% methods_to_keep,]
dt_sd = dt_sd[dt$methodlabel %in% methods_to_keep,]


#table 1
dt_1 = dt[dt$setting_id %in% c(1:11),]
dt_1$m1 = rep(NA,nrow(dt_1))
dt_1$effect = rep(NA,nrow(dt_1))

dt_1$m1[dt_1$setting_id %in% (c(1,3,5,7,9))+1 ] = 10
dt_1$m1[dt_1$setting_id%in% (c(1,3,5,7,9))+2 ] = 100
dt_1$m1[dt_1$setting_id == 30] = 0

dt_1$effect[dt_1$setting_id == 1] = 0
dt_1$effect[dt_1$setting_id %in% c(2,3)] = 0.5
dt_1$effect[dt_1$setting_id %in% c(4,5)] = 1.0
dt_1$effect[dt_1$setting_id %in% c(6,7)] = 1.5
dt_1$effect[dt_1$setting_id %in% c(8,9)] = 2.0
dt_1$effect[dt_1$setting_id %in% c(10,11)] = 2.5

#install.packages('reshape2')
library(reshape2)
names(dt_1)
method_names = as.character(unique(dt_1$methodlabel))
dt_tp = matrix(NA,nrow = length(1:11),ncol  = length(method_names))
dt_rejected = matrix(NA,nrow = length(1:11),ncol  = length(method_names))
dt_FDR = matrix(NA,nrow = length(1:11),ncol  = length(method_names))
for(i in 1:length(method_names)){
  for(j in 1:11){
    tp_row = which(dt_1$setting_id == j & dt_1$methodlabel == method_names[i])
    dt_tp[j,i] = dt_1$tp[tp_row[1]]
    dt_rejected[j,i] = dt_1$rejected[tp_row[1]]
    dt_FDR[j,i] = dt_1$fdr[tp_row[1]]
  }
}



dt_FDR = round(dt_FDR,2)
dt_tp = round(dt_tp,2)

colnames(dt_FDR) = method_names
colnames(dt_tp) = method_names

dt_tp = as.data.frame(dt_tp)
dt_tp$m1 =c(0,rep(c(10,100),5))
dt_tp$effect =c(0,rep(0.5,2),rep(1.0,2),rep(1.5,2),rep(2.0,2),rep(2.5,2))

dt_FDR = as.data.frame(dt_FDR)
dt_FDR$m1 =c(0,rep(c(10,100),5))
dt_FDR$effect =c(0,rep(0.5,2),rep(1.0,2),rep(1.5,2),rep(2.0,2),rep(2.5,2))

library(xtable)
names(dt_tp)
cols_reorder = 1:ncol(dt_tp)
dt_tp = dt_tp[,cols_reorder]
names(dt_tp)

dt_FDR= dt_FDR[,cols_reorder]
names(dt_FDR)


print(xtable(dt_tp), include.rownames=FALSE)

print(xtable(dt_FDR), include.rownames=FALSE)



Compute_SE_for_Scenarios = function(scenarios_vec,B=96,col = 'fdr'){
  return(max(as.numeric(as.character(dt_sd[as.character(dt_sd$setting_id) %in% as.character(c(scenarios_vec)),col])))/sqrt(B))
}

Compute_SE_for_Scenarios(1:11)
Compute_SE_for_Scenarios(1:11,col = 'tp')

library(latex2exp)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


new_col_names =  c('ALDEx2-t','ALDEx2-W','ANCOM','DACOMP-t-ratio','DACOMP-t-subsample','DACOMP-ratio','DACOMP-subsample','HG','W-CSS','W-TSS','W-FLOW','Wrench','m1','effect')
new_col_names_no_flow = new_col_names[-c(11,13,14)]
dt_FDR_Plot = dt_FDR

names(dt_FDR_Plot) =new_col_names
write.csv(dt_FDR_Plot,file = '../../Results/sim1_fdr.csv')
dt_FDR_Plot = dt_FDR_Plot[, !(colnames(dt_FDR_Plot)%in% c('WCOMP-ratio','WCOMP-t-ratio'))]

library(yarrr)
pallete = yarrr::piratepal("xmen",
                           plot.result = FALSE)#[c(1,2,3,4,10)]          
pallete = substr(pallete,1,7)

methods_main_body = c('ALDEx2-t','ANCOM','DACOMP-subsample','HG','W-CSS','W-FLOW')
methods_appendix = c('ALDEx2-W','DACOMP-t-subsample','W-TSS','Wrench')
methods_list = list(methods_main_body = methods_main_body, methods_appendix = methods_appendix)

for(graph_id in c(1,2)){
  dt_FDR_Plot_melt = reshape2::melt(dt_FDR_Plot,id.vars = c('m1','effect'), value.name = 'FDR')
  dt_FDR_Plot_melt = dt_FDR_Plot_melt[dt_FDR_Plot_melt$variable  %in% methods_list[[graph_id]],]
  dt_FDR_Plot_melt_GLOBAL_NULL_m1_10 = dt_FDR_Plot_melt[which(dt_FDR_Plot_melt$m1==0),]
  dt_FDR_Plot_melt_GLOBAL_NULL_m1_100 = dt_FDR_Plot_melt[which(dt_FDR_Plot_melt$m1==0),]
  dt_FDR_Plot_melt_GLOBAL_NULL_m1_10$m1 = 10
  dt_FDR_Plot_melt_GLOBAL_NULL_m1_100$m1 = 100
  dt_FDR_Plot_melt = rbind(dt_FDR_Plot_melt,dt_FDR_Plot_melt_GLOBAL_NULL_m1_10,dt_FDR_Plot_melt_GLOBAL_NULL_m1_100)
  names(dt_FDR_Plot_melt)[3] = 'Method'
  names(dt_FDR_Plot_melt)[2] = 'Effect'
  dt_FDR_Plot_melt = dt_FDR_Plot_melt[-which(dt_FDR_Plot_melt$m1==0),]
  dt_FDR_Plot_melt$m1 = factor(dt_FDR_Plot_melt$m1,levels = c(10,100),labels = c('m1 = 10','m1 = 100'))
  p_1_FDR = ggplot(dt_FDR_Plot_melt,aes(x=Effect,y = FDR,color = Method))+geom_line(lwd = 0.75) + geom_point(size = 0.4) +
    facet_wrap(m1~.) + theme_bw()+geom_hline(yintercept = 0.1,color = 'black',lty = 1,alpha=0.2)+xlab(TeX("$\\lambda_{effect}$")) +
    scale_color_manual(values=as.character(pallete[c(1,2,3,4,5,6)]))  #gg_color_hue(10)[1:10]
  ggsave(p_1_FDR,filename = paste0('../../Results/sim_p1_FDR',graph_id,'.pdf'),width = 7,height = 3)
}



dt_TP_Plot = dt_tp
names(dt_TP_Plot) =new_col_names
write.csv(dt_TP_Plot,file = '../../Results/sim1_tp.csv')
dt_TP_Plot = dt_TP_Plot[, !(colnames(dt_TP_Plot)%in% c('WCOMP-ratio','WCOMP-t-ratio'))]

for(graph_id in c(1,2)){
  dt_TP_Plot_melt = reshape2::melt(dt_TP_Plot,id.vars = c('m1','effect'),value.name = 'TP')  
  dt_TP_Plot_melt = dt_TP_Plot_melt[dt_TP_Plot_melt$variable   %in% methods_list[[graph_id]],]
  dt_TP_Plot_melt = dt_TP_Plot_melt[-which(dt_TP_Plot_melt$m1 == 0),]
  names(dt_TP_Plot_melt)[3] = 'Method'
  names(dt_TP_Plot_melt)[2] = 'Effect'
  dt_TP_Plot_melt$m1 = factor(dt_TP_Plot_melt$m1,levels = c(10,100),labels = c('m1 = 10','m1 = 100'))
  p_1_TP = ggplot(dt_TP_Plot_melt,aes(x=Effect,y = TP,color = Method))+geom_line(lwd = 0.75)+geom_point(size = 0.4)+
    facet_wrap(m1~.,scales = "free") + theme_bw() +xlab(TeX("$\\lambda_{effect}$"))+
    scale_color_manual(values=as.character(pallete[c(1,2,3,4,5,6)]))
  ggsave(p_1_TP,filename = paste0('../../Results/sim_p1_TP',graph_id,'.pdf'),width = 7,height = 3)
}


#table 2

dt_2 = dt[dt$setting_id %in% c(12:21),]
dt_2$m1 = rep(NA,nrow(dt_2))
dt_2$p_high = rep(NA,nrow(dt_2))

dt_2$m1[dt_2$setting_id %in% c(17:21) ] = 120
dt_2$m1[dt_2$setting_id %in% c(12:16) ] = 60


dt_2$p_high[dt_2$setting_id %in% c(12,17)] = 0.9
dt_2$p_high[dt_2$setting_id %in% c(13,18)] = 0.8
dt_2$p_high[dt_2$setting_id %in% c(14,19)] = 0.7
dt_2$p_high[dt_2$setting_id %in% c(15,20)] = 0.6
dt_2$p_high[dt_2$setting_id %in% c(16,21)] = 0.5

#install.packages('reshape2')
library(reshape2)
names(dt_2)
method_names = as.character(unique(dt_2$methodlabel))
dt_tp = matrix(NA,nrow = length(12:21),ncol  = length(method_names))
dt_rejected = matrix(NA,nrow = length(12:21),ncol  = length(method_names))
dt_FDR = matrix(NA,nrow = length(12:21),ncol  = length(method_names))
for(i in 1:length(method_names)){
  for(j in 1:10){
    tp_row = which(dt_2$setting_id == j+11 & dt_2$methodlabel == method_names[i])
    dt_tp[j,i] = dt_2$tp[tp_row[1]]
    dt_rejected[j,i] = dt_2$rejected[tp_row[1]]
    dt_FDR[j,i] = dt_2$fdr[tp_row[1]]
  }
}


dt_tp
dt_FDR

for(i in 1:ncol(dt_tp)){
  dt_tp[,i] = as.character(round(as.numeric(dt_tp[,i]),0))
}
for(i in 1:ncol(dt_tp)){
  dt_FDR[,i] = as.character(round(as.numeric(dt_FDR[,i]),2))
}

print_by_Xtable = function(dt_combined,cols_to_remove){
  colnames(dt_combined) = new_col_names_no_flow
  dt_combined = dt_combined
  dt_combined = as.data.frame(dt_combined)
  dt_combined$m1 =c(rep(120,5),rep(60,5))
  dt_combined$p_high =rep(c(0.9,0.8,0.7,0.6,0.5),2)
  print(xtable(dt_combined[,-cols_to_remove]), include.rownames=FALSE)
  return(dt_combined)
}

cols_remove_main = c(2,4,5,6,10,11)
cols_remove_appendix = c(1,3,7,8,9)


write.csv(print_by_Xtable(dt_FDR,cols_remove_main),file = '../../Results/sim2_FDR.csv')
write.csv(print_by_Xtable(dt_FDR,cols_remove_appendix),file = '../../Results/sim2_FDR.csv')
write.csv(print_by_Xtable(dt_tp,cols_remove_main),file = '../../Results/sim2_tp.csv')
write.csv(print_by_Xtable(dt_tp,cols_remove_appendix),file = '../../Results/sim2_tp.csv')

Compute_SE_for_Scenarios(12:21)
Compute_SE_for_Scenarios(12:21,col = 'tp')


#####

#table 3

dt_3 = dt[dt$setting_id %in% c(22:25),]




#install.packages('reshape2')
library(reshape2)

method_names = as.character(unique(dt_3$methodlabel))
dt_tp = matrix(NA,nrow = length(22:25),ncol  = length(method_names))
dt_rejected = matrix(NA,nrow = length(22:25),ncol  = length(method_names))
dt_FDR = matrix(NA,nrow = length(22:25),ncol  = length(method_names))
for(i in 1:length(method_names)){
  for(j in 1:4){
    tp_row = which(dt_3$setting_id == j+21 & dt_3$methodlabel == method_names[i])
    dt_tp[j,i] = dt_3$tp[tp_row[1]]
    dt_rejected[j,i] = dt_3$rejected[tp_row[1]]
    dt_FDR[j,i] = dt_3$fdr[tp_row[1]]
  }
}


dt_tp = cbind(c('15:15','20:20','25:25','30:30'),round(dt_tp,2))
dt_FDR = cbind(c('15:15','20:20','25:25','30:30'),round(dt_FDR,2))
colnames(dt_tp) = c('n',new_col_names_no_flow)
colnames(dt_FDR) = c('n',new_col_names_no_flow)


write.csv(dt_FDR,row.names = F,file = '../../Results/sim3_FDR.csv')
write.csv(dt_tp,row.names = F,file = '../../Results/sim3_TP.csv')

print(xtable(dt_FDR[,-c(3,5,6,7,11,12)]), include.rownames=FALSE)#main body
print(xtable(dt_FDR[,-c(2,4,8,9)]), include.rownames=FALSE)
print(xtable(dt_tp[,-c(3,5,6,7,11,12,9)]), include.rownames=FALSE) # -9 = remove HG
print(xtable(dt_tp[,-c(2,4,8,9)]), include.rownames=FALSE) # -9 = remove HG

Compute_SE_for_Scenarios(22:25)
Compute_SE_for_Scenarios(22:25,col = 'tp')




