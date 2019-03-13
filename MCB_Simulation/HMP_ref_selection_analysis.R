#HMP ref selection analysis

dt_list = list()
dt_list[[1]] = read.csv('E:/MCB2/HMP2/Results_Thres_1_0_Select_by_BH/results.csv')
dt_list[[2]] = read.csv('E:/MCB2/HMP2/Results_Thres_1_1_Select_by_BH/results.csv')
dt_list[[3]] = read.csv('E:/MCB2/HMP2/Results_Thres_1_2_Select_by_BH/results.csv')
dt_list[[4]] = read.csv('E:/MCB2/HMP2/Results_Thres_1_3_Select_by_BH/results.csv')
dt_list[[5]] = read.csv('E:/MCB2/HMP2/Results_Thres_1_4_Select_by_BH/results.csv')
# 
# ind_rows_to_take = 1:36
# for(i in 1:length(dt_list)){
#   dt_list[[i]] = dt_list[[i]][ind_rows_to_take,]
# }
ind_rows_to_take = 1:nrow(dt_list[[1]])
x_vec = c(1,1.1,1.2,1.3,1.4)
#View(dt_list[[1]])
pdf('E:/MCB2/HMP2/ref_selection_thres_comp.pdf',width = 16,height = 12)
#comp_to_look_at = 4
for(comp_to_look_at in ind_rows_to_take){
  print(comp_to_look_at)
  N = length(dt_list)
  ref_size_selected = rep(NA,N)
  taxa_for_test = rep(NA,N)
  ds_fdr_total_rejections = rep(NA,N)
  ds_fdr_shared = rep(NA,N)
  median_lambda = rep(NA,N)
  Mult1_Rej_in_Ref = rep(NA,N)
  lambda_multiplier = rep(NA,N)
  rejected_in_reference = rep(NA,N)
  rejections = rep(NA,N)
  rejections_shared = rep(NA,N)
  Mult1_Rej = rep(NA,N)
  for(i in 1:N){
    
    taxa_for_test[i] = dt_list[[i]]$taxa_for_test[comp_to_look_at]
    ds_fdr_total_rejections[i] = dt_list[[i]]$ds_fdr_total_rejections[comp_to_look_at]
    ds_fdr_shared[i] = dt_list[[i]]$ds_fdr_shared[comp_to_look_at]
    rejections[i] = dt_list[[i]]$rejections[comp_to_look_at]
    rejections_shared[i] = dt_list[[i]]$shared[comp_to_look_at]
    ref_size_selected[i] = dt_list[[i]]$ref_size_selected[comp_to_look_at]
    rejected_in_reference[i] = dt_list[[i]]$rejected_in_reference[comp_to_look_at]
    Mult1_Rej_in_Ref[i] = dt_list[[i]]$Mult1_Rej_in_Ref[comp_to_look_at]
    
    median_lambda[i] = dt_list[[i]]$median_lambda[comp_to_look_at]
    
    lambda_multiplier[i] = dt_list[[i]]$lambda_multiplier[comp_to_look_at]
    Mult1_Rej[i] = dt_list[[i]]$Mult1_Rejections[comp_to_look_at]
  }
  
  ylim_current = max(max(ref_size_selected),max(taxa_for_test)) +50
  ylim_current_ref = ylim_current
  if(comp_to_look_at>36){
    ylim_current = 30
  }
  PCH =18
  CEX = 1.5
  par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
  if(comp_to_look_at<=36){
    plot(x_vec,taxa_for_test,main = 'Taxa-tested, and discovered',pch=PCH,cex=  CEX,ylim = c(0,ylim_current))
    points(x_vec,ds_fdr_total_rejections,col = 'red',pch=PCH,cex=  CEX,ylim = c(0,ylim_current))
    points(x_vec,rejections,col = 'red',pch=1,cex=  CEX,ylim = c(0,ylim_current))
    points(x_vec,ds_fdr_shared,col = 'blue',pch=PCH,cex=  CEX,ylim = c(0,ylim_current))
    legend(1.0, ylim_current, legend=c("Tested", "DS-FDR","BH",'DS-FDR Shared'),
           col=c("black", "red",'blue','red'),pch=c(18,18,1,18) , cex=0.8)  
  }else{
    plot(x_vec,taxa_for_test,main = 'Taxa-tested, and discovered',pch=PCH,cex=  CEX,ylim = c(0,ylim_current))
    points(x_vec,rejections,col = 'blue',pch=1,cex=  CEX,ylim = c(0,ylim_current))
    points(x_vec,rejections_shared,col = 'blue',pch=2,cex=  CEX,ylim = c(0,ylim_current))
    legend(1.0, ylim_current, legend=c("Tested","BH",'BH Shared'),
           col=c("black",'blue','blue'),pch=c(18,1,2) , cex=0.8)  
  }
  
  
  plot(x_vec,ref_size_selected,main = 'Ref. and rejected in Ref',pch=PCH,cex=  CEX,ylim = c(0,ylim_current_ref))
  points(x_vec,rejected_in_reference,col = 'red',pch=PCH,cex=  CEX,ylim = c(0,ylim_current))
  points(x_vec,Mult1_Rej_in_Ref,col = 'blue',pch=PCH,cex=  CEX,ylim = c(0,ylim_current))
  points(x_vec,Mult1_Rej,col = 'green',pch=PCH,cex=  CEX,ylim = c(0,ylim_current))
  legend(1.0, ylim_current_ref, legend=c("Ref. Size", "Rejected in ref","Mult1.rej in ref",'Mult1.rej'),
         col=c("black", "red",'blue','green'),pch=c(18,18,18,18) , cex=0.8)
  
  plot(x_vec,median_lambda,main = 'median lambda')
  if(sum(is.na(lambda_multiplier))==0)
    plot(x_vec,lambda_multiplier,main = 'lambda multiplier')
  
  
  Sites_title = paste0(dt_list[[1]]$Site1[comp_to_look_at],' vs. ',dt_list[[1]]$Site2[comp_to_look_at])
  mtext(Sites_title, outer = TRUE, cex = 1.5)
  
  par(mfrow=c(1,1))
}


dev.off()