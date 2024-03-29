#Version 5 add a gray scale version of the graph with results in the main text

#graphs for paper
library(ggplot2)

#Used to filter HG from the power results of graphs for 4.2.
FILTER_HG_FROM_POWER_GUT_SIM = FALSE

NR.REPS.FOR.SE = 100 # this is the number of data generations used for converting SDs to SEs (for the point estiamtes of power and FDR)

#load point estimates and results
dt = read.csv('../../Results/REFSIM_Combined_Results.csv',stringsAsFactors = F)
dt_sd = read.csv('../../Results/REFSIM_Combined_Results_sd.csv',stringsAsFactors = F)

#Addition to version 3 - add in the ZINB-WAVE results
load('../../Results/zinbwave_deseq2_results.rdata')
for(i in 1:nrow(zinbwave_deseq2_results)){
  setting_name = row.names(zinbwave_deseq2_results)[i] #as.character(dt[dt$setting_id==i,"setting_name"][1])
  row_to_add_dt = c(nrow(dt)+1,
                    'ZINB-WAVE',
                    i,
                    zinbwave_deseq2_results$rejected_avg[i],
                    zinbwave_deseq2_results$TP_avg[i],
                    NA,NA,
                    zinbwave_deseq2_results$FDR_avg[i],
                    setting_name)
  row_to_add_dt_SD = c(nrow(dt_sd)+1,
                    'ZINB-WAVE',
                    i,
                    zinbwave_deseq2_results$rejected_se[i]*sqrt(NR.REPS.FOR.SE),
                    zinbwave_deseq2_results$TP_se[i]*sqrt(NR.REPS.FOR.SE),
                    NA,NA,
                    zinbwave_deseq2_results$FDR_se[i]*sqrt(NR.REPS.FOR.SE),
                    setting_name)
  dt = rbind(dt,row_to_add_dt)
  dt_sd = rbind(dt_sd,row_to_add_dt_SD)
}


#Addition to version 4 - add in the DACOMP variants with automatic tuning of Scrit
load('../../Results/dacomp_results_ref_by_counts.rdata')

for(i in 1:nrow(dacomp_results)){
  setting_name = row.names(dacomp_results)[i]
  
  #DACOMP with Wilcoxon
  row_to_add_dt = c(nrow(dt)+1,
                    'DACOMP,Wilcoxon',
                    i,
                    dacomp_results$rejected_avg_Wilcoxon[i],
                    dacomp_results$TP_avg_Wilcoxon[i],
                    NA,NA,
                    dacomp_results$FDR_avg_Wilcoxon[i],
                    setting_name)
  
  row_to_add_dt_SD = c(nrow(dt_sd)+1,
                       'DACOMP,Wilcoxon',
                       i,
                       dacomp_results$rejected_se_Wilcoxon[i]*sqrt(NR.REPS.FOR.SE),
                       dacomp_results$TP_se_Wilcoxon[i]*sqrt(NR.REPS.FOR.SE),
                       NA,NA,
                       dacomp_results$FDR_se_Wilcoxon[i]*sqrt(NR.REPS.FOR.SE),
                       setting_name)
  
  dt = rbind(dt,row_to_add_dt)
  dt_sd = rbind(dt_sd,row_to_add_dt_SD)
  
  #DACOMP-ratio with Wilcoxon
  row_to_add_dt = c(nrow(dt)+1,
                    'DACOMP-Ratio,Wilcoxon',
                    i,
                    dacomp_results$rejected_avg_Wilcoxon_Ratio[i],
                    dacomp_results$TP_avg_Wilcoxon_Ratio[i],
                    NA,NA,
                    dacomp_results$FDR_avg_Wilcoxon_Ratio[i],
                    setting_name)
  row_to_add_dt_SD = c(nrow(dt_sd)+1,
                       'DACOMP-Ratio,Wilcoxon',
                       i,
                       dacomp_results$rejected_se_Wilcoxon_Ratio[i]*sqrt(NR.REPS.FOR.SE),
                       dacomp_results$TP_se_Wilcoxon_Ratio[i]*sqrt(NR.REPS.FOR.SE),
                       NA,NA,
                       dacomp_results$FDR_se_Wilcoxon_Ratio[i]*sqrt(NR.REPS.FOR.SE),
                       setting_name)
  dt = rbind(dt,row_to_add_dt)
  dt_sd = rbind(dt_sd,row_to_add_dt_SD)
  
  #DACOMP with Welch
  row_to_add_dt = c(nrow(dt)+1,
                    'DACOMP,Welch',
                    i,
                    dacomp_results$rejected_avg_Welch[i],
                    dacomp_results$TP_avg_Welch[i],
                    NA,NA,
                    dacomp_results$FDR_avg_Welch[i],
                    setting_name)
  
  row_to_add_dt_SD = c(nrow(dt_sd)+1,
                       'DACOMP,Welch',
                       i,
                       dacomp_results$rejected_se_Welch[i]*sqrt(NR.REPS.FOR.SE),
                       dacomp_results$TP_se_Welch[i]*sqrt(NR.REPS.FOR.SE),
                       NA,NA,
                       dacomp_results$FDR_se_Welch[i]*sqrt(NR.REPS.FOR.SE),
                       setting_name)
  
  dt = rbind(dt,row_to_add_dt)
  dt_sd = rbind(dt_sd,row_to_add_dt_SD)
}

# We show our methods for a specific S, results for other S were similar.
methods_to_keep = c('ALDEx2_Welch',
'ALDEx2_Wilcoxon',
'ANCOM',
'DACOMP,Wilcoxon',
'DACOMP-Ratio,Wilcoxon',
'DACOMP,Welch',
'HG,S = 1.3, Oracle',
'WILCOXON_PAULSON',
'WILCOXON_PERCENT',
'WILCOXON_FLOW',
'Wrench',
'ZINB-WAVE')

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
dt_FDR = matrix(round(as.numeric(dt_FDR),2),ncol =ncol(dt_FDR))
dt_tp = matrix(round(as.numeric(dt_tp),2),ncol =ncol(dt_tp))

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
max(Compute_SE_for_Scenarios(c(2,4,6,8,10,26),col = 'tp')/10,
Compute_SE_for_Scenarios(c(2,4,6,8,10,26)+1,col = 'tp')/100)
library(latex2exp)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# The next section of the script produces the graphs for power and FDR, given in subsection 4.1 and appendix A.

new_col_names =  c('ALDEx2-t','ALDEx2-W','ANCOM','HG','W-FLOW','W-CSS','W-TSS','Wrench','ZINB-WAVE','DACOMP','DACOMP-ratio','DACOMP-t','m1','effect')
#new_col_names_no_flow = new_col_names[-c(11,13,14)]
dt_FDR_Plot = dt_FDR

#Note - check via: cbind(names(dt_FDR_Plot),new_col_names)
names(dt_FDR_Plot) =new_col_names
write.csv(dt_FDR_Plot,file = '../../Results/sim1_fdr.csv')
#dt_FDR_Plot = dt_FDR_Plot[, !(colnames(dt_FDR_Plot)%in% c('WCOMP-ratio','WCOMP-t-ratio'))]

library(yarrr)
pallete = yarrr::piratepal("xmen",
                           plot.result = FALSE)#[c(1,2,3,4,10)]          
pallete = substr(pallete,1,7)

methods_main_body = c('ALDEx2-t','ANCOM','DACOMP-ratio','DACOMP','HG','W-CSS','W-FLOW')
methods_appendix = c('ALDEx2-W','DACOMP-t','W-TSS','Wrench','ZINB-WAVE')
methods_list = list(methods_main_body = methods_main_body, methods_appendix = methods_appendix)

Melted_list_FDR = list()
Melted_list_TP = list()
Melted_list_Power = list()
plots_list_FDR = list()
plots_list_TP = list()
plots_list_Power = list()

plots_list_FDR_vertical = list()
plots_list_Power_vertical = list()
plots_list_Power_vertical_with_HG = list()
BASE_FONT_SIZE_VERTICAL = 15
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
  Melted_list_FDR[[graph_id]] = dt_FDR_Plot_melt
  p_1_FDR = ggplot(dt_FDR_Plot_melt,aes(x=Effect,y = FDR,color = Method))+geom_line(lwd = 0.75) + geom_point(size = 0.4) +
    facet_wrap(m1~.) + theme_bw()+geom_hline(yintercept = 0.1,color = 'black',lty = 1,alpha=0.2,lwd = 1.25)+xlab(TeX("$\\\\Delta_{ml}$")) + xlim(c(0,3))+
    scale_color_manual(values=as.character(pallete[c(1,2,3,4,5,6,7)]))
  plots_list_FDR[[graph_id]] = p_1_FDR
  
  p_1_vertical = ggplot(dt_FDR_Plot_melt,aes(x=Effect,y = FDR,color = Method))+geom_line(lwd = 0.75) + geom_point(size = 0.4) +
    facet_grid(m1~.) + theme_bw(base_size = BASE_FONT_SIZE_VERTICAL)+geom_hline(yintercept = 0.1,color = 'black',lty = 1,alpha=0.2,lwd = 1.25)+xlab(TeX("$\\\\Delta_{ml}$")) + xlim(c(0,3))+
    scale_color_manual(values=as.character(pallete[c(1,2,3,4,5,6,7)]))
  
  plots_list_FDR_vertical[[graph_id]] = p_1_vertical
  ggsave(p_1_FDR,filename = paste0('../../Results/sim_p1_FDR',graph_id,'.pdf'),width = 7,height = 3)
}

# two graphs for Power - main body and appendix
dt_TP_Plot = dt_tp
names(dt_TP_Plot) =new_col_names
write.csv(dt_TP_Plot,file = '../../Results/sim1_tp.csv')
if(FILTER_HG_FROM_POWER_GUT_SIM)
  dt_TP_Plot = dt_TP_Plot[, !(colnames(dt_TP_Plot)%in% c('HG'))]

for(graph_id in c(1,2)){
  dt_TP_Plot_melt = reshape2::melt(dt_TP_Plot,id.vars = c('m1','effect'),value.name = 'TP')  
  dt_TP_Plot_melt = dt_TP_Plot_melt[dt_TP_Plot_melt$variable   %in% methods_list[[graph_id]],]
  dt_TP_Plot_melt = dt_TP_Plot_melt[-which(dt_TP_Plot_melt$m1 == 0),]
  names(dt_TP_Plot_melt)[3] = 'Method'
  names(dt_TP_Plot_melt)[2] = 'Effect'
  dt_TP_Plot_melt$m1 = factor(dt_TP_Plot_melt$m1,levels = c(10,100),labels = c('m1 = 10','m1 = 100'))
  Melted_list_TP[[graph_id]] = dt_TP_Plot_melt
  p_1_TP = ggplot(dt_TP_Plot_melt,aes(x=Effect,y = TP,color = Method))+geom_line(lwd = 0.75)+geom_point(size = 0.4)+
    facet_wrap(m1~.,scales = "free") + theme_bw() +xlab(TeX("$\\Delta_{ml}$"))+xlim(c(0.5,3))+
    scale_color_manual(values=as.character(pallete[c(1,2,3,4,5,6,7)]))
  plots_list_TP[[graph_id]] = p_1_TP
  ggsave(p_1_TP,filename = paste0('../../Results/sim_p1_TP',graph_id,'.pdf'),width = 7,height = 3)
  
  dt_TP_Plot_melt$Power = dt_TP_Plot_melt$TP/(10^as.numeric(dt_TP_Plot_melt$m1))
  Melted_list_Power[[graph_id]] = dt_TP_Plot_melt$Power
  p_1_Power = ggplot(dt_TP_Plot_melt,aes(x=Effect,y = Power,color = Method))+geom_line(lwd = 0.75)+geom_point(size = 0.4)+
    facet_wrap(m1~.) + theme_bw() +xlab(TeX("$\\\\Delta_{ml}$"))+xlim(c(0.5,3))+
    scale_color_manual(values=as.character(pallete[c(1,2,3,4,5,6,7)]))
  
  plots_list_Power[[graph_id]] = p_1_Power
  
  p_1_Power_vertical = ggplot(dt_TP_Plot_melt,aes(x=Effect,y = Power,color = Method))+geom_line(lwd = 0.75)+geom_point(size = 0.4)+
    facet_grid(m1~.) + theme_bw(base_size = BASE_FONT_SIZE_VERTICAL) +xlab(TeX("$\\\\Delta_{ml}$"))+xlim(c(0.5,3))+
    scale_color_manual(values=as.character(pallete[c(1,2,3,4,5,6,7)]))
  plots_list_Power_vertical[[graph_id]] = p_1_Power_vertical
}

#Generate graphs with no margin
pvt_TP_m1_10 = Melted_list_TP[[1]][Melted_list_TP[[1]]$m1 == "m1 = 10",]
pvt_TP_m1_100 = Melted_list_TP[[1]][Melted_list_TP[[1]]$m1 == "m1 = 100",]
pvt_FDR_m1_10 = Melted_list_FDR[[1]][Melted_list_FDR[[1]]$m1 == "m1 = 10",]
pvt_FDR_m1_100 = Melted_list_FDR[[1]][Melted_list_FDR[[1]]$m1 == "m1 = 100",]
names(pvt_TP_m1_10)[4] = "Value"
names(pvt_TP_m1_100)[4] = "Value"
names(pvt_FDR_m1_10)[4] = "Value"
names(pvt_FDR_m1_100)[4] = "Value"

pvt_TP_m1_10$Type = "Power"
pvt_TP_m1_100$Type = "Power"
pvt_FDR_m1_10$Type = "FDR"
pvt_FDR_m1_100$Type = "FDR"

pvt_TP_m1_10$Value = pvt_TP_m1_10$Value/10
pvt_TP_m1_100$Value = pvt_TP_m1_100$Value/100

combined_melted = rbind(pvt_TP_m1_10,pvt_TP_m1_100,pvt_FDR_m1_10,pvt_FDR_m1_100)

pvt_FDR_horizontal_line = pvt_FDR_m1_10[pvt_FDR_m1_10$Method == "ANCOM",]
pvt_FDR_horizontal_line_copy = pvt_FDR_horizontal_line
pvt_FDR_horizontal_line_copy$m1 = "m1 = 100"
pvt_FDR_horizontal_line = rbind(pvt_FDR_horizontal_line,pvt_FDR_horizontal_line_copy)

pvt_FDR_horizontal_line$Method = as.character(pvt_FDR_horizontal_line$Method)
pvt_FDR_horizontal_line$Method = rep("FDR = 0.1",nrow(pvt_FDR_horizontal_line))
pvt_FDR_horizontal_line$Value= 0.1

pvt_FDR_horizontal_line$alpha = 0.2
pvt_FDR_horizontal_line$lwd = 1.25
pvt_FDR_horizontal_line$lty = 2

combined_melted$alpha = 0
combined_melted$lwd = 0.75
combined_melted$lty = 1

combined_melted_plus_horizontal = combined_melted
levels_to_keep = levels(combined_melted_plus_horizontal$Method)
combined_melted_plus_horizontal$Method = as.character(combined_melted_plus_horizontal$Method)


#combined_melted_plus_horizontal = rbind(combined_melted_plus_horizontal,pvt_FDR_horizontal_line)
combined_melted_plus_horizontal$Method = factor(x = combined_melted_plus_horizontal$Method,levels = c('ALDEx2-t','ANCOM','DACOMP-ratio','DACOMP','HG','W-FLOW','W-CSS')) #,'FDR = 0.1'

pallete_combined = c(pallete,"grey" = "grey")

# p_1_Combined = ggplot(combined_melted,aes(x=Effect,y = Value,color = Method))+geom_line(lwd = 0.75)+geom_point(size = 0.4)+
#   facet_grid(Type~m1,scales = "free_y") + theme_bw(base_size = BASE_FONT_SIZE_VERTICAL) +xlab(TeX("$\\\\Delta_{ml}$"))+xlim(c(0,3))+
#   scale_color_manual(values=as.character(pallete[c(1,2,3,4,6,7,5)])) +geom_hline(yintercept = 0.1,color = 'black',lty = 1,alpha=0.2,lwd = 1.25)
combined_melted_plus_horizontal$lty = as.factor(combined_melted_plus_horizontal$lty)
p_1_Combined_h_line = ggplot(combined_melted_plus_horizontal,aes(x=Effect,y = Value,color = Method,lty = lty))+geom_line(lwd=0.75)+geom_point(size = 0.4)+
  facet_grid(Type~m1,scales = "free_y") + theme_bw(base_size = BASE_FONT_SIZE_VERTICAL) +xlab(TeX("$\\\\Delta_{ml}$"))+xlim(c(0,3))+
  scale_color_manual(values=as.character(pallete_combined[c(1,2,3,4,6,7,5,9)]))+ylab("")  + guides(linetype = FALSE)
p_1_Combined_h_line+theme(panel.spacing = unit(0, "lines"))

ggsave(p_1_Combined_h_line+theme(panel.spacing = unit(0, "lines")),filename = paste0('../../Results/AOAS_combined.pdf'),width = 7*1.5,height = 4.5)
#This is a revised graph for the AOAS version

# gray scale version for AOAS
override.shape <- c(4,3,15,16,4,15,16)
override.lty <- c(1,1,1,1,2,2,2)
col_vec = pallete_combined[c(1,2,3,4,6,7,5,9)[-8]]

FDR_line_dummy_dt = data.frame(value = c(0.1,0.1,NA,NA),
                               m1 = c('m1 = 10','m1 = 100','m1 = 10','m1 = 100'),
                               Type = c('FDR','FDR','Power','Power'))

p_1_Combined_h_line_gs = ggplot(combined_melted_plus_horizontal,
                                aes(x=Effect,
                                    y = Value,
                                    color = Method,
                                    shape = Method,
                                    lty = Method))+
  geom_line(lwd=0.4)+geom_point(size = 1.7)+
  facet_grid(Type~m1,scales = "free_y") +
  theme_bw(base_size = BASE_FONT_SIZE_VERTICAL) +
  xlab(TeX("$\\\\Delta_{ml}$"))+
  xlim(c(0,3))+
  scale_color_manual(values=as.character(col_vec))+
  scale_shape_manual(values = override.shape)+
  scale_linetype_manual(values = override.lty)+
  scale_alpha_manual(values = c(rep(0.3,7),1))+
  ylab("")+ 
  theme(panel.spacing = unit(0, "lines")) +
  #scale_linetype(guide='none') +
  guides(color = guide_legend(override.aes = list(shape = c(override.shape),
                                                  linetype = c(override.lty),
                                                  color = as.character(col_vec))))+
  geom_hline(mapping = aes(yintercept = value),data = FDR_line_dummy_dt,
             col='black',lty=2)
p_1_Combined_h_line_gs
ggsave(p_1_Combined_h_line_gs,filename = paste0('../../Results/AOAS_combined_greyscale.pdf'),width = 7*1.5,height = 4.5)


library(cowplot)
temp = plot_grid(plots_list_FDR[[1]]+xlab(""), NULL, plots_list_Power[[1]]+ theme(legend.position = "none",plot.margin=unit(c(5.5, 112, 5.5, 5.5), "points")) ,
          nrow = 3,ncol = 1,rel_heights = c(0.5, -0.04, 0.5),rel_widths = c(1.0, 1.0, 0.5))
ggsave(temp,filename = paste0('../../Results/AOAS_combined_attempt_1.pdf'),width = 7*1.5,height = 4.5)
temp = plot_grid(plots_list_Power_vertical[[1]]+ theme(legend.position = "none"), plots_list_FDR_vertical[[1]] ,
          nrow = 1,ncol =2,rel_widths = c(0.41, 0.59)) #c(0.4305, 0.5695)
ggsave(temp,filename = paste0('../../Results/AOAS_combined_attempt_2.pdf'),width = 7*1.5,height = 4.5)


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
  if(!is.null(cols_to_remove)){
    dt_combined =cbind(m1,p_high,dt_combined[,-cols_to_remove])  
  }else{
    dt_combined =cbind(m1,p_high,dt_combined)
  }
  print(xtable(t(dt_combined)), include.rownames=TRUE)
  return(dt_combined)
}

# we now produce two sets of tables, for the appendix and the main test, for FDR and power.
cols_remove_main = c(2,4,9,10,11) #cols will not be displayed in main body
cols_remove_appendix = c(1,3,5,6,7,8) #cols will not be displayed in appendix

# print to file AND show latex
#write.csv(print_by_Xtable(dt_FDR,cols_remove_main),file = '../../Results/sim2_FDR.csv')
#write.csv(print_by_Xtable(dt_FDR,cols_remove_appendix),file = '../../Results/sim2_FDR.csv')
#write.csv(print_by_Xtable(dt_tp,cols_remove_main),file = '../../Results/sim2_tp.csv')
#write.csv(print_by_Xtable(dt_tp,cols_remove_appendix),file = '../../Results/sim2_tp.csv')
write.csv(print_by_Xtable(dt_FDR,NULL),file = '../../Results/sim2_FDR.csv')
write.csv(print_by_Xtable(dt_tp,NULL),file = '../../Results/sim2_tp.csv')

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
dt_tp = cbind(c('15:15','20:20','25:25','30:30'),matrix(round(as.numeric(dt_tp),2),ncol =ncol(dt_tp)))
dt_FDR = cbind(c('15:15','20:20','25:25','30:30'),matrix(round(as.numeric(dt_FDR),2),ncol =ncol(dt_FDR)))
colnames(dt_tp) = c('n',method_names)
colnames(dt_FDR) = c('n',method_names)

#write results to file and latex table
write.csv(dt_FDR,row.names = F,file = '../../Results/sim3_FDR.csv')
write.csv(dt_tp,row.names = F,file = '../../Results/sim3_TP.csv')

#print(xtable(dt_FDR[,-c(3,5,10,11,12)]), include.rownames=FALSE)#main body
#print(xtable(dt_FDR[,c(1,3,5,10,11,12)]), include.rownames=FALSE)
#print(xtable(dt_tp[,-c(3,5,10,11,8,12)]), include.rownames=FALSE) #  remove HG - was not valid at any scenario!
#print(xtable(dt_tp[,c(1,3,5,10,11,12)]), include.rownames=FALSE) #  remove HG


print(xtable(t(dt_FDR)), include.rownames=T)
print(xtable(t(dt_tp[,-5])), include.rownames=T)

#compute maximum SE for power and FDR (across scenarios)
Compute_SE_for_Scenarios(22:25)
Compute_SE_for_Scenarios(22:25,col = 'tp')

#Cases with confounder

dt_4 = dt[dt$setting_id %in% c(28,30,32,34,36,38,40),] #crop scenario results from the main table
dt_5 = dt[dt$setting_id %in% (c(28,30,32,34,36,38,40)+1),] #crop scenario results from the main table


dt_4$m1 = rep(NA,nrow(dt_4))
dt_5$m1 = rep(NA,nrow(dt_5))
dt_4$lambda = rep(NA,nrow(dt_4))
dt_5$lambda = rep(NA,nrow(dt_5))

dt_4$m1[dt_4$setting_id == 28] = 10
dt_4$m1[dt_4$setting_id == 30] = 10
dt_4$m1[dt_4$setting_id == 32] = 100
dt_4$m1[dt_4$setting_id == 34] = 10
dt_4$m1[dt_4$setting_id == 36] = 100
dt_4$m1[dt_4$setting_id == 38] = 10
dt_4$m1[dt_4$setting_id == 40] = 100

dt_5$m1[dt_5$setting_id == 28+1] = 10
dt_5$m1[dt_5$setting_id == 30+1] = 10
dt_5$m1[dt_5$setting_id == 32+1] = 100
dt_5$m1[dt_5$setting_id == 34+1] = 10
dt_5$m1[dt_5$setting_id == 36+1] = 100
dt_5$m1[dt_5$setting_id == 38+1] = 10
dt_5$m1[dt_5$setting_id == 40+1] = 100

dt_4$lambda[dt_4$setting_id == 28] = 0
dt_4$lambda[dt_4$setting_id == 30] = 1
dt_4$lambda[dt_4$setting_id == 32] = 1
dt_4$lambda[dt_4$setting_id == 34] = 2
dt_4$lambda[dt_4$setting_id == 36] = 2
dt_4$lambda[dt_4$setting_id == 38] = 3
dt_4$lambda[dt_4$setting_id == 40] = 3


dt_5$lambda[dt_5$setting_id == 28+1] = 0
dt_5$lambda[dt_5$setting_id == 30+1] = 1
dt_5$lambda[dt_5$setting_id == 32+1] = 1
dt_5$lambda[dt_5$setting_id == 34+1] = 2
dt_5$lambda[dt_5$setting_id == 36+1] = 2
dt_5$lambda[dt_5$setting_id == 38+1] = 3
dt_5$lambda[dt_5$setting_id == 40+1] = 3

crop_global_null = dt_4[dt_4$lambda==0,]
crop_global_null$m1 = 100
dt_4 = rbind(dt_4,crop_global_null)
dt_4$Artifical_Disease_Group = "Group S Oversampled,"
crop_global_null = dt_5[dt_5$lambda==0,]
crop_global_null$m1 = 100
dt_5 = rbind(dt_5,crop_global_null)
dt_5$Artifical_Disease_Group = "Group S Undersampled,"

dt_confounder = rbind(dt_4,dt_5)
dt_confounder$methodlabel = as.character(dt_confounder$methodlabel)
colnames(dt_confounder)[colnames(dt_confounder) == 'methodlabel'] = 'Method'
replace_vec = unique(dt_confounder$Method)
replace_to = c('ALDEx2-t','ALDEx2-W','ANCOM','HG','W-FLOW','W-CSS','W-TSS','Wrench','ZINB-WAVE','DACOMP','DACOMP-ratio','DACOMP-t')
for(i in 1:length(replace_vec)){
  dt_confounder$Method[dt_confounder$Method == replace_vec[i]] = replace_to[i]
}
dt_confounder$m1 = factor(dt_confounder$m1,levels = c(10,100),labels = c('m1 = 10','m1 = 100'))
dt_confounder_main_method = dt_confounder[dt_confounder$Method %in% methods_main_body,]
dt_confounder_main_method$Method = factor(dt_confounder_main_method$Method,levels = methods_main_body)
dt_confounder_appendix_method = dt_confounder[dt_confounder$Method %in% methods_appendix,]
dt_confounder_appendix_method$Method  = factor(dt_confounder_appendix_method$Method,levels = methods_appendix)

library(ggplot2)
library(yarrr)
pallete = yarrr::piratepal("xmen",
                           plot.result = FALSE)#[c(1,2,3,4,10)]          
pallete = substr(pallete,1,7)

for(i in 1:2){
  #i=1
  dt_list = list(dt_confounder_main_method,dt_confounder_appendix_method)
  dt_current = dt_list[[i]]
  dt_current$fdr = as.numeric(dt_current$fdr)
  dt_current$tp = as.numeric(dt_current$tp)
  dt_current$Method = as.character(dt_current$Method)
  
  pal_vec = c(1,2,3,4,5,7,6)
  if(i == 2)
    pal_vec = c(1,2,3,4,7,6)
  
  p_fdr = ggplot(dt_current,aes(x = lambda,y=fdr,color = Method))+
    geom_line(lwd = 0.75) + geom_point(size = 0.4)+
    facet_wrap(Artifical_Disease_Group~m1) +
    theme_bw()+geom_hline(yintercept = 0.1,color = 'black',lty = 1,alpha=0.2,lwd = 1.25)+
    xlab(TeX("$\\\\Delta_{ml}$"))+scale_color_manual(values=as.character(pallete[pal_vec]))+ylab('FDR')
  
  
  
  p_tp = ggplot(dt_current[dt_current$Method != 'HG',],aes(x = lambda,y=tp,color = Method))+
    geom_line(lwd = 0.75) + geom_point(size = 0.4)+
    facet_wrap(Artifical_Disease_Group~m1) +
    theme_bw()+geom_hline(yintercept = 0.1,color = 'black',lty = 1,alpha=0.2,lwd = 1.25)+
    xlab(TeX("$\\\\Delta_{ml}$"))+scale_color_manual(values=as.character(pallete[c(1,2,3,4,7,6)]))+ylab('Diff. Abun. Taxa Discovered')
  
  
  
  p_tp_m1_100 = ggplot(dt_current[dt_current$m1 == 'm1 = 100' & dt_current$Method != 'HG',],aes(x = lambda,y=tp,color = Method))+
    geom_line(lwd = 0.75) + geom_point(size = 0.4)+
    facet_wrap(Artifical_Disease_Group~m1) +
    theme_bw()+geom_hline(yintercept = 0.1,color = 'black',lty = 1,alpha=0.2,lwd = 1.25)+
    xlab(TeX("$\\\\Delta_{ml}$"))+scale_color_manual(values=as.character(pallete[c(1,2,3,4,7,6)]))+ylab('Diff. Abun. Taxa Discovered')
  ret = list()
  
  ggsave(p_fdr,filename = paste0('../../Results/sim_p4_FDR_',i,'.pdf'),width = 7,height = 5)
  ggsave(p_tp,filename = paste0('../../Results/sim_p4_TP_',i,'.pdf'),width = 7,height = 5)
  ggsave(p_tp_m1_100,filename = paste0('../../Results/sim_p4_TP_m1_100_',i,'.pdf'),width = 5,height = 5)
  
  
}

se_fdr = Compute_SE_for_Scenarios(28:41)
se_fdr
Compute_SE_for_Scenarios(28:41,col = 'tp')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Settings with negative binomial sequencing depth
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dt_6 = dt[dt$setting_id %in% c(42:47),] #crop scenario results from the main table

# put columns for additional parameters, describing each scenario
dt_6$m1 = rep(NA,nrow(dt_6)) 
dt_6$Percent_Increase = rep(NA,nrow(dt_6))

dt_6$m1[dt_6$setting_id %in% c(42,44,46) ] = 10
dt_6$m1[dt_6$setting_id %in% c(43,45,47) ] = 100


dt_6$Percent_Increase[dt_6$setting_id %in% c(42:43)] = 0.1
dt_6$Percent_Increase[dt_6$setting_id %in% c(44:45)] = 0.2
dt_6$Percent_Increase[dt_6$setting_id %in% c(46:47)] = 0.3

#compile tables: 


library(reshape2)
names(dt_6)
method_names = as.character(unique(dt_6$methodlabel))
dt_tp = matrix(NA,nrow = length(42:47),ncol  = length(method_names))
dt_rejected = matrix(NA,nrow = length(42:47),ncol  = length(method_names))
dt_FDR = matrix(NA,nrow = length(42:47),ncol  = length(method_names))
for(i in 1:length(method_names)){
  for(j in 1:6){
    tp_row = which(dt_6$setting_id == j+41 & dt_6$methodlabel == method_names[i])
    dt_tp[j,i] = dt_6$tp[tp_row[1]]
    dt_rejected[j,i] = dt_6$rejected[tp_row[1]]
    dt_FDR[j,i] = dt_6$fdr[tp_row[1]]
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
  m1 =c(10,10,10,100,100,100)
  percent_increase =c(0.1,0.2,0.3,0.1,0.2,0.3)
  if(!is.null(cols_to_remove)){
    dt_combined =cbind(m1,percent_increase,dt_combined[,-cols_to_remove])  
  }else{
    dt_combined =cbind(m1,percent_increase,dt_combined)
  }
  print(xtable(t(dt_combined)), include.rownames=TRUE)
  return(dt_combined)
}

write.csv(print_by_Xtable(dt_FDR[c(1,3,5,2,4,6),],NULL),file = '../../Results/sim6_FDR.csv')
write.csv(print_by_Xtable(dt_tp[c(1,3,5,2,4,6),],NULL),file = '../../Results/sim6_tp.csv')

# compute max standard errors (across scenarios) for this FDR and TP
Compute_SE_for_Scenarios(42:47)
Compute_SE_for_Scenarios(42:47,col = 'tp')


