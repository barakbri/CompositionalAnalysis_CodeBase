library(cowplot)
library(ggplot2)

MAINDIR = 'E:/MCB2/'
RESULTS_DIR = paste0(MAINDIR,"/Results/")
REFSIM_aggregated_results_file = function(RESULTS_DIR,SCENARIO_ID){
  return(paste0(RESULTS_DIR,'/REFSIM_aggregated_results_',SCENARIO_ID,'.RData'))
}

SCENARIO_TO_PLOT = 36
for(SCENARIO_TO_PLOT in 18:68){ #1:59
  
  filename = REFSIM_aggregated_results_file(RESULTS_DIR,SCENARIO_TO_PLOT)
  load(file = filename)
  
  Rows_subzero = (grep('TargetAbun',aggregated_results$methodlabel ))
  Row_ANCOM    = which(aggregated_results$methodlabel == 'ANCOM')#(grep('ANCOM',aggregated_results$methodlabel ))
  Rows_RANDOM    = (grep('Random',aggregated_results$methodlabel ))
  Rows_Wilcoxon    = which(aggregated_results$methodlabel == 'WILCOXON')
  Rows_Wilcoxon_PERCENT    = which(aggregated_results$methodlabel == 'WILCOXON_PERCENT')#
  Rows_Wilcoxon_PAULSON    = which(aggregated_results$methodlabel == 'WILCOXON_PAULSON')#
  Rows_Wilcoxon_qPCR    = which(aggregated_results$methodlabel == 'WILCOXON_qPCR')#
  
  Rows_Oracle_All  = (grep('Oracle, All',aggregated_results$methodlabel ))
  MODE_SHOW_ANCOM = T
  
  Rows_subzero_wilcoxon = Rows_subzero[13:24]
  Rows_subzero_LAD = Rows_subzero[1:12]
  
  
  #rows_to_take = c(10:12,14:16,5,7,9,6,8,10,2)
  reorder = c(10,11,12,1,2,3,4,6,8,5,7,9)
  dt_plot = aggregated_results[Rows_subzero_wilcoxon[reorder],]
  dt_plot = rbind(dt_plot,aggregated_results[Rows_Oracle_All[2],])
  dt_plot_LAD = aggregated_results[Rows_subzero_LAD[reorder],]
  dt_plot_LAD = rbind(dt_plot_LAD,aggregated_results[Rows_Oracle_All[2],])
  
  PS_VEC = c(rep(c('+1','Abundance','Prev'),4),'')
  
  Target_Abundance_vec = c(rep(5,3),rep(10,3),rep(25,6),'')
  Oracle_vec = c(rep(F,9),rep(T,3),T)
  
  #TwoPart_vec = c(rep(c(F,T),12),F,T)
  Labels = c('TA = 5',rep('',2),
             'TA = 10',rep('',2),
             'TA = 25',rep('',2),
             'Oracle, TA = 25,',rep('',2),
             'Oracle, All')
  
  x_axis = 1:(length(Labels))
  
  dt_plot$x_axis  = as.factor(x_axis)
  dt_plot$PS_VEC = as.factor(PS_VEC)
  dt_plot$Target_Abundance_vec = as.factor(Target_Abundance_vec)
  dt_plot$Oracle_vec = as.factor(Oracle_vec)
  
  dt_plot_LAD$x_axis  = as.factor(x_axis)
  dt_plot_LAD$PS_VEC = as.factor(PS_VEC)
  dt_plot_LAD$Target_Abundance_vec = as.factor(Target_Abundance_vec)
  dt_plot_LAD$Oracle_vec = as.factor(Oracle_vec)
  
  #dt_plot$Labels = Labels
  #dt_plot$TwoPart_vec = TwoPart_vec
  Rejections = ggplot(data = dt_plot) + geom_point(aes(x = x_axis,y = tp,color = PS_VEC)) + scale_x_discrete(labels = Labels) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle('# True Positive') + xlab('Ref. selection methods') + ylab('')
  BadReference = ggplot(data = dt_plot) + geom_point(aes(x = x_axis,y = nr_bad_reference,color = PS_VEC)) + scale_x_discrete(labels = Labels) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle('# H1 in reference set') + xlab('Ref. selection methods')+ ylab('')
  References_Significant = ggplot(data = dt_plot) + geom_point(aes(x = x_axis,y = reference_significant,color = PS_VEC)) + scale_x_discrete(labels = Labels) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle('# References Significant') + xlab('Ref. selection methods')+ ylab('')
  FDR = ggplot(data = dt_plot) + geom_point(aes(x = x_axis,y = fdr,color = PS_VEC)) + scale_x_discrete(labels = Labels) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle('FDR') + xlab('Ref. selection methods')+ ylab('')
  
  Rejections_LAD = ggplot(data = dt_plot_LAD) + geom_point(aes(x = x_axis,y = tp,color = PS_VEC)) + scale_x_discrete(labels = Labels) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle('# True Positive') + xlab('Ref. selection methods') + ylab('')
  BadReference_LAD = ggplot(data = dt_plot_LAD) + geom_point(aes(x = x_axis,y = nr_bad_reference,color = PS_VEC)) + scale_x_discrete(labels = Labels) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle('# H1 in reference set') + xlab('Ref. selection methods')+ ylab('')
  References_Significant_LAD = ggplot(data = dt_plot_LAD) + geom_point(aes(x = x_axis,y = reference_significant,color = PS_VEC)) + scale_x_discrete(labels = Labels) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle('# References Significant') + xlab('Ref. selection methods')+ ylab('')
  FDR_LAD = ggplot(data = dt_plot_LAD) + geom_point(aes(x = x_axis,y = fdr,color = PS_VEC)) + scale_x_discrete(labels = Labels) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle('FDR') + xlab('Ref. selection methods')+ ylab('')
  
  
  #
  FDR_LAD = FDR + geom_hline(aes(yintercept=0.1),color = 'gray') + geom_text(aes(10,0.05,label = "nominal 0.1", vjust = - 2),color = 'gray')
  FDR_LAD = FDR + geom_hline(aes(yintercept=0.1),color = 'gray') + geom_text(aes(10,0.05,label = "nominal 0.1", vjust = - 2),color = 'gray')
  
  title_size = 10
  xlab_size = 10
  
  Rejections = Rejections + theme(plot.title = element_text(size=title_size),axis.title.x = element_text(size = xlab_size))
  BadReference = BadReference + theme(plot.title = element_text(size=title_size),axis.title.x = element_text(size = xlab_size))
  References_Significant = References_Significant + theme(plot.title = element_text(size=title_size),axis.title.x = element_text(size = xlab_size))
  FDR = FDR + theme(plot.title = element_text(size=title_size),axis.title.x = element_text(size = xlab_size))
  
  
  Rejections_LAD = Rejections_LAD + theme(plot.title = element_text(size=title_size),axis.title.x = element_text(size = xlab_size))
  BadReference_LAD = BadReference_LAD + theme(plot.title = element_text(size=title_size),axis.title.x = element_text(size = xlab_size))
  References_Significant_LAD = References_Significant_LAD + theme(plot.title = element_text(size=title_size),axis.title.x = element_text(size = xlab_size))
  FDR_LAD = FDR_LAD + theme(plot.title = element_text(size=title_size),axis.title.x = element_text(size = xlab_size))
  
  
  Value_Wilcoxon_Rejections = aggregated_results$tp[Rows_Wilcoxon_PAULSON]
  Rejections = Rejections + geom_hline(aes(yintercept=Value_Wilcoxon_Rejections)) + geom_text(aes(7,Value_Wilcoxon_Rejections,label = "Wilcoxon, Paulson", vjust = 2))
  Rejections_LAD = Rejections_LAD + geom_hline(aes(yintercept=Value_Wilcoxon_Rejections)) + geom_text(aes(7,Value_Wilcoxon_Rejections,label = "Wilcoxon, Paulson", vjust = 2))
  Value_Wilcoxon_FDR = aggregated_results$fdr[Rows_Wilcoxon_PAULSON]
  FDR_LAD = FDR_LAD + geom_hline(aes(yintercept=Value_Wilcoxon_FDR)) + geom_text(aes(7,Value_Wilcoxon_FDR,label = "Wilcoxon, Paulson", vjust = -4))
  
  
  
  
  if(length(Rows_Wilcoxon_qPCR)>1){
    Value_Wilcoxon_qPCR_Rejections = aggregated_results$tp[Rows_Wilcoxon_qPCR]
    Rejections = Rejections + geom_hline(aes(yintercept=Value_Wilcoxon_qPCR_Rejections),color = 'red') + geom_text(aes(10,Value_Wilcoxon_qPCR_Rejections,label = "Wilcoxon, qPCR Corrected", vjust = 4),color = 'red')
    Rejections_LAD = Rejections_LAD + geom_hline(aes(yintercept=Value_Wilcoxon_qPCR_Rejections),color = 'red') + geom_text(aes(10,Value_Wilcoxon_qPCR_Rejections,label = "Wilcoxon, qPCR Corrected", vjust = 4),color = 'red')
    Value_Wilcoxon_qPCR_FDR = aggregated_results$fdr[Rows_Wilcoxon_qPCR]
    FDR = FDR + geom_hline(aes(yintercept=Value_Wilcoxon_qPCR_FDR),color = 'red') + geom_text(aes(10,Value_Wilcoxon_qPCR_FDR,label = "Wilcoxon, qPCR Corrected", vjust = - 2),color = 'red')
    FDR_LAD = FDR_LAD + geom_hline(aes(yintercept=Value_Wilcoxon_qPCR_FDR),color = 'red') + geom_text(aes(10,Value_Wilcoxon_qPCR_FDR,label = "Wilcoxon, qPCR Corrected", vjust = - 2),color = 'red')
  }
  
  if(MODE_SHOW_ANCOM){
    Value_ANCOM_Rejections = aggregated_results$tp[Row_ANCOM]
    Rejections = Rejections + geom_hline(aes(yintercept=Value_ANCOM_Rejections),color = 'green') + geom_text(aes(10,Value_ANCOM_Rejections,label = "ANCOM", vjust = 4),color = 'green')
    Rejections_LAD = Rejections_LAD + geom_hline(aes(yintercept=Value_ANCOM_Rejections),color = 'green') + geom_text(aes(10,Value_ANCOM_Rejections,label = "ANCOM", vjust = 4),color = 'green')
    Value_ANCOM_FDR = aggregated_results$fdr[Row_ANCOM]
    FDR = FDR + geom_hline(aes(yintercept=Value_ANCOM_FDR),color = 'green') + geom_text(aes(10,Value_ANCOM_FDR,label = "ANCOM", vjust = - 2),color = 'green')    
    FDR_LAD = FDR_LAD + geom_hline(aes(yintercept=Value_ANCOM_FDR),color = 'green') + geom_text(aes(10,Value_ANCOM_FDR,label = "ANCOM", vjust = - 2),color = 'green')    
  }
  
  
  
  plot_grid(Rejections, FDR,References_Significant,BadReference,nrow = 2,ncol = 2) + draw_label(REFSIM_SETTINGS_LIST[[SCENARIO_TO_PLOT]]$label,vjust = -47, fontface='bold',size = 10)
  
  ggsave(paste0('E:/MCB2/Results/REFSIM_plot_scenario_',SCENARIO_TO_PLOT,'.png'),width = 10,height = 10)
  
  plot_grid(Rejections_LAD, FDR_LAD,References_Significant_LAD,BadReference_LAD,nrow = 2,ncol = 2) + draw_label(REFSIM_SETTINGS_LIST[[SCENARIO_TO_PLOT]]$label,vjust = -47, fontface='bold',size = 10)
  
  ggsave(paste0('E:/MCB2/Results/REFSIM_plot_scenario_',SCENARIO_TO_PLOT,'_LAD.png'),width = 10,height = 10)
  
}

