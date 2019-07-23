#This script produces the graphs and tables shown in section 4, and appendix A.

#NOTE: NEED TO REVISE THIS FILE IN TERM OF INDICES!!

#graphs for paper
library(ggplot2)

NR.REPS.FOR.SE = 96 # this is the number of data generations used for converting SDs to SEs (for the point estiamtes of power and FDR)

#load point estimates and results
dt = read.csv('../../Results/REFSIM_Combined_Results.csv')
dt_sd = read.csv('../../Results/REFSIM_Combined_Results_sd.csv')

# We show our methods for a specific S, results for other S were similar.
methods_to_keep = c('ALDEx2_Welch',
'ALDEx2_Wilcoxon',
'ANCOM',
'DACOMP,Welch,rarefaction,S = 1.3',
'DACOMP,Wilcoxon,division,S = 1.3',
'DACOMP,Wilcoxon,rarefaction,S = 1.3',
'HG,S = 1.3, Oracle',
'WILCOXON_PAULSON',
'WILCOXON_PERCENT',
'WILCOXON_FLOW',
'Wrench')

dt = dt[dt$methodlabel %in% methods_to_keep,]
dt_sd = dt_sd[dt$methodlabel %in% methods_to_keep,]


# Results for subsection 4.1 - simulations based on sampling from a dataset
#parse table for plotting:
dt_1 = dt[dt$setting_id %in% c(1:11,26:27),] #crop only the scenarios we want - note that 26,27 were added in a later version of the manuscript, hence the gap in ordering
dt_1$m1 = rep(NA,nrow(dt_1)) #create columns for nr of diff abun taxa, and effect size
dt_1$effect = rep(NA,nrow(dt_1))
# place values in new columns - this is consistent with the defs in REFSIM_GenerateSettins_Index.R
dt_1$m1[dt_1$setting_id %in% (c(1,3,5,7,9))+1 ] = 10
dt_1$m1[dt_1$setting_id%in% (c(1,3,5,7,9))+2 ] = 100
dt_1$m1[dt_1$setting_id %in% c(26) ] = 10
dt_1$m1[dt_1$setting_id %in% c(27) ] = 100
dt_1$m1[dt_1$setting_id%in% (c(1,3,5,7,9))+2 ] = 100

dt_1$effect[dt_1$setting_id == 1] = 0
dt_1$effect[dt_1$setting_id %in% c(2,3)] = 0.5
dt_1$effect[dt_1$setting_id %in% c(4,5)] = 1.0
dt_1$effect[dt_1$setting_id %in% c(6,7)] = 1.5
dt_1$effect[dt_1$setting_id %in% c(8,9)] = 2.0
dt_1$effect[dt_1$setting_id %in% c(10,11)] = 2.5
dt_1$effect[dt_1$setting_id %in% c(26,27)] = 3.0

#install.packages('reshape2')
library(reshape2)
# create table for number of TP and FR, by method and scenario
names(dt_1)
settings_index = c(1:11,26:27)
method_names = as.character(unique(dt_1$methodlabel))
dt_tp = matrix(NA,nrow = length(settings_index),ncol  = length(method_names))
dt_rejected = matrix(NA,nrow = length(settings_index),ncol  = length(method_names))
dt_FDR = matrix(NA,nrow = length(settings_index),ncol  = length(method_names))
for(i in 1:length(method_names)){
  for(j in 1:length(settings_index)){
    tp_row = which(dt_1$setting_id == settings_index[j] & dt_1$methodlabel == method_names[i])
    dt_tp[j,i] = dt_1$tp[tp_row[1]]
    dt_rejected[j,i] = dt_1$rejected[tp_row[1]]
    dt_FDR[j,i] = dt_1$fdr[tp_row[1]]
  }
}


#additional columsn and method names
dt_FDR = round(dt_FDR,2)
dt_tp = round(dt_tp,2)

colnames(dt_FDR) = method_names
colnames(dt_tp) = method_names

dt_tp = as.data.frame(dt_tp)
dt_tp$m1 =c(0,rep(c(10,100),6))
dt_tp$effect =c(0,rep(0.5,2),rep(1.0,2),rep(1.5,2),rep(2.0,2),rep(2.5,2),rep(3.0,2))

dt_FDR = as.data.frame(dt_FDR)
dt_FDR$m1 =c(0,rep(c(10,100),6))
dt_FDR$effect =c(0,rep(0.5,2),rep(1.0,2),rep(1.5,2),rep(2.0,2),rep(2.5,2),rep(3.0,2))

#generate latex tables
library(xtable)
names(dt_tp)
cols_reorder = 1:ncol(dt_tp)
dt_tp = dt_tp[,cols_reorder]
names(dt_tp)

dt_FDR= dt_FDR[,cols_reorder]
names(dt_FDR)


print(xtable(dt_tp), include.rownames=FALSE)

print(xtable(dt_FDR), include.rownames=FALSE)


# compute SE for FDR and TP, this is the maximum standard error
Compute_SE_for_Scenarios = function(scenarios_vec,B=NR.REPS.FOR.SE,col = 'fdr'){
  return(max(as.numeric(as.character(dt_sd[as.character(dt_sd$setting_id) %in% as.character(c(scenarios_vec)),col])))/sqrt(B))
}

se_fdr = Compute_SE_for_Scenarios(1:11)
se_fdr
Compute_SE_for_Scenarios(1:11,col = 'tp')

library(latex2exp)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# The next section of the script produces the graphs for power and FDR, given in subsection 4.1 and appendix A.

new_col_names =  c('ALDEx2-t','ALDEx2-W','ANCOM','DACOMP-t','DACOMP-ratio','DACOMP','HG','W-CSS','W-TSS','W-FLOW','Wrench','m1','effect')
new_col_names_no_flow = new_col_names[-c(11,13,14)]
dt_FDR_Plot = dt_FDR

names(dt_FDR_Plot) =new_col_names
write.csv(dt_FDR_Plot,file = '../../Results/sim1_fdr.csv')
dt_FDR_Plot = dt_FDR_Plot[, !(colnames(dt_FDR_Plot)%in% c('WCOMP-ratio','WCOMP-t-ratio'))]

library(yarrr)
pallete = yarrr::piratepal("xmen",
                           plot.result = FALSE)#[c(1,2,3,4,10)]          
pallete = substr(pallete,1,7)

methods_main_body = c('ALDEx2-t','ANCOM','DACOMP-ratio','DACOMP','HG','W-CSS','W-FLOW')
methods_appendix = c('ALDEx2-W','DACOMP-t','W-TSS','Wrench')
methods_list = list(methods_main_body = methods_main_body, methods_appendix = methods_appendix)

# two graphs for FDR- main body and appendix
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
    facet_wrap(m1~.) + theme_bw()+geom_hline(yintercept = 0.1,color = 'black',lty = 1,alpha=0.2)+xlab(TeX("$\\lambda_{effect}$")) + xlim(c(0,3))+
    scale_color_manual(values=as.character(pallete[c(1,2,3,4,5,6,7)])) +  geom_hline(yintercept = 0.1+2*se_fdr ,color = 'black',lty = 2,alpha=0.2) #gg_color_hue(10)[1:10]
  ggsave(p_1_FDR,filename = paste0('../../Results/sim_p1_FDR',graph_id,'.pdf'),width = 7,height = 3)
}

# two graphs for Power - main body and appendix
dt_TP_Plot = dt_tp
names(dt_TP_Plot) =new_col_names
write.csv(dt_TP_Plot,file = '../../Results/sim1_tp.csv')
dt_TP_Plot = dt_TP_Plot[, !(colnames(dt_TP_Plot)%in% c('HG'))]

for(graph_id in c(1,2)){
  dt_TP_Plot_melt = reshape2::melt(dt_TP_Plot,id.vars = c('m1','effect'),value.name = 'TP')  
  dt_TP_Plot_melt = dt_TP_Plot_melt[dt_TP_Plot_melt$variable   %in% methods_list[[graph_id]],]
  dt_TP_Plot_melt = dt_TP_Plot_melt[-which(dt_TP_Plot_melt$m1 == 0),]
  names(dt_TP_Plot_melt)[3] = 'Method'
  names(dt_TP_Plot_melt)[2] = 'Effect'
  dt_TP_Plot_melt$m1 = factor(dt_TP_Plot_melt$m1,levels = c(10,100),labels = c('m1 = 10','m1 = 100'))
  p_1_TP = ggplot(dt_TP_Plot_melt,aes(x=Effect,y = TP,color = Method))+geom_line(lwd = 0.75)+geom_point(size = 0.4)+
    facet_wrap(m1~.,scales = "free") + theme_bw() +xlab(TeX("$\\lambda_{effect}$"))+xlim(c(0.5,3))+
    scale_color_manual(values=as.character(pallete[c(1,2,3,4,6,7,5)]))
  ggsave(p_1_TP,filename = paste0('../../Results/sim_p1_TP',graph_id,'.pdf'),width = 7,height = 3)
}


# This section is used for producing the power and FDR results, for the simulations
# in subsection 4.2 - settings with no over dispersion

dt_2 = dt[dt$setting_id %in% c(12:21),] #crop scenario results from the main table

# put columns for additional parameters, describing each scenario
dt_2$m1 = rep(NA,nrow(dt_2)) 
dt_2$p_high = rep(NA,nrow(dt_2))

dt_2$m1[dt_2$setting_id %in% c(17:21) ] = 120
dt_2$m1[dt_2$setting_id %in% c(12:16) ] = 60


dt_2$p_high[dt_2$setting_id %in% c(12,17)] = 0.9
dt_2$p_high[dt_2$setting_id %in% c(13,18)] = 0.8
dt_2$p_high[dt_2$setting_id %in% c(14,19)] = 0.7
dt_2$p_high[dt_2$setting_id %in% c(15,20)] = 0.6
dt_2$p_high[dt_2$setting_id %in% c(16,21)] = 0.5

#compile tables: 


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


colnames(dt_tp) = method_names
colnames(dt_FDR) = method_names

#round results, 
for(i in 1:ncol(dt_tp)){
  dt_tp[,i] = as.character(round(as.numeric(dt_tp[,i]),0))
}
for(i in 1:ncol(dt_tp)){
  dt_FDR[,i] = as.character(round(as.numeric(dt_FDR[,i]),2))
}

#aux function used for printing both to latex and CSV
print_by_Xtable = function(dt_combined,cols_to_remove){
  dt_combined = dt_combined
  dt_combined = as.data.frame(dt_combined)
  m1 =c(rep(120,5),rep(60,5))
  p_high =rep(c(0.9,0.8,0.7,0.6,0.5),2)
  dt_combined =cbind(m1,p_high,dt_combined[,-cols_to_remove])
  print(xtable(dt_combined), include.rownames=FALSE)
  return(dt_combined)
}

# we now produce two sets of tables, for the appendix and the main test, for FDR and power.
cols_remove_main = c(2,4,9,10) #cols will not be displayed in main body
cols_remove_appendix = c(1,3,5,6,7,8) #cols will not be displayed in appendix

# print to file AND show latex
write.csv(print_by_Xtable(dt_FDR,cols_remove_main),file = '../../Results/sim2_FDR.csv')
write.csv(print_by_Xtable(dt_FDR,cols_remove_appendix),file = '../../Results/sim2_FDR.csv')
write.csv(print_by_Xtable(dt_tp,cols_remove_main),file = '../../Results/sim2_tp.csv')
write.csv(print_by_Xtable(dt_tp,cols_remove_appendix),file = '../../Results/sim2_tp.csv')

# compute max standard errors (across scenarios) for this FDR and TP
Compute_SE_for_Scenarios(12:21)
Compute_SE_for_Scenarios(12:21,col = 'tp')


#####

# This section of the script is used for generating the Power and FDR tables shown in subsection 4.1 (and the correspong part of appendix A)

dt_3 = dt[dt$setting_id %in% c(22:25),] #crop scenarios

#install.packages('reshape2')
library(reshape2)

#aggregate power and FDR to a single table
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

# add sample sizes
dt_tp = cbind(c('15:15','20:20','25:25','30:30'),round(dt_tp,2))
dt_FDR = cbind(c('15:15','20:20','25:25','30:30'),round(dt_FDR,2))
colnames(dt_tp) = c('n',method_names)
colnames(dt_FDR) = c('n',method_names)

#write results to file and latex table
write.csv(dt_FDR,row.names = F,file = '../../Results/sim3_FDR.csv')
write.csv(dt_tp,row.names = F,file = '../../Results/sim3_TP.csv')

print(xtable(dt_FDR[,-c(3,5,10,11)]), include.rownames=FALSE)#main body
print(xtable(dt_FDR[,c(1,3,5,10,11)]), include.rownames=FALSE)
print(xtable(dt_tp[,-c(3,5,10,11,8)]), include.rownames=FALSE) #  remove HG - was not valid at any scenario!
print(xtable(dt_tp[,c(1,3,5,10,11)]), include.rownames=FALSE) #  remove HG

#compute maximum SE for power and FDR (across scenarios)
Compute_SE_for_Scenarios(22:25)
Compute_SE_for_Scenarios(22:25,col = 'tp')




