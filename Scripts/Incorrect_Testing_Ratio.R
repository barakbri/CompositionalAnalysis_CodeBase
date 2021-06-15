# This script is used to reproduce Figure 1 on subsection 2.1
# This script compares several approaches for tests of differential abundance.

#Output file:

set.seed(1)
library(dacomp)
library(dirmult)

generate_data = function(n_0 = 50,
                         n_1 = 50,
                         N = 10000,
                         P_0 = c(90,0.1,10),
                         w = 0.75){
  P_1 = (1-w)*P_0 + (w)*c(1,0,0)
  X_0 = matrix(NA,nrow = n_0,ncol = 3)
  X_1 = matrix(NA,nrow = n_1,ncol = 3)
  for(i in 1:n_0){
    P_Sample = rdirichlet(n=1, alpha=P_0)
    X_0[i,] = rmultinom(n = 1,size = N,P_Sample)[,1]
  }
  for(i in 1:n_1){
    P_Sample = rdirichlet(n=1, alpha=P_1)
    X_1[i,] = rmultinom(n = 1,size = N,P_Sample)[,1]
  }
  Y = c(rep(0,n_0),rep(1,n_1))
  return(list(X = rbind(X_0,X_1),Y = Y))
}


nr.reps = 10000
alpha  = 0.1
W_grid = c(0,1/4,2/4,3/4)

T1E_Wilcoxon = rep(NA,length(W_grid))
T1E_Welch = rep(NA,length(W_grid))


for(w_ind in 1:length(W_grid)){
  print(paste0('W: ',W_grid[w_ind]))
  
  Pval_ratio_Wilcoxon = rep(NA,nr.reps)
  Pval_ratio_Welch = rep(NA,nr.reps)

  for(b in 1:nr.reps){
    if(b%%100==1){
      print(b) 
    }
    
    data = generate_data(w = W_grid[w_ind])
    dacomp_res = dacomp.test(X = data$X,y = data$Y,ind_reference_taxa = 3,
                             test = DACOMP.TEST.NAME.WILCOXON,disable_DSFDR = T,
                             nr_perm = 1000,compute_ratio_normalization = T)
    Pval_ratio_Wilcoxon[b] = dacomp_res$p.values.test.ratio.normalization[2]
    dacomp_res = dacomp.test(X = data$X,y = data$Y,ind_reference_taxa = 3,
                             test = DACOMP.TEST.NAME.WELCH,disable_DSFDR = T,
                             nr_perm = 1000,compute_ratio_normalization = T)
    Pval_ratio_Welch[b] = dacomp_res$p.values.test.ratio.normalization[2]  
  }
  if(any(is.na(Pval_ratio_Welch))){
   stop('Zero rarefaction depth!') 
  }
  T1E_Wilcoxon[w_ind] = mean(Pval_ratio_Wilcoxon<=alpha,na.rm = T)
  T1E_Welch[w_ind] = mean(Pval_ratio_Welch<=alpha,na.rm = T)
  
}


#plot(ecdf(Pval_ratio_Welch))
#abline(0,1,col='red')

res = data.frame(W = W_grid,T1E_Wilcoxon = T1E_Wilcoxon,T1E_Welch = T1E_Welch)
res

save(res,file = '../../Results/Ratio_Inflation_Results_temp.Rdata')

set.seed(1)
data_for_plot_w_0 = generate_data(n_0 = 100000,n_1 = 100000,w = 0, P_0 = c(90,0.1,10))
data_for_plot_w_0_5 = generate_data(n_0 = 100000,n_1 = 100000,w = 0.5, P_0 = c(90,0.1,10))
data_for_plot_w_0_75 = generate_data(n_0 = 100000,n_1 = 100000,w = 0.75, P_0 = c(90,0.1,10))
ind = 1:100000+ 100000
dt_comb_1 = data.frame(ratio =  data_for_plot_w_0$X[ind,2]/
                         (data_for_plot_w_0$X[ind,2]+data_for_plot_w_0$X[ind,3]),
                       w = rep(0,100000))
dt_comb_2 = data.frame(ratio =  data_for_plot_w_0_5$X[ind,2]/
                         (data_for_plot_w_0_5$X[ind,2]+data_for_plot_w_0_5$X[ind,3]),
                       w = rep(0.5,100000))
dt_comb_3 = data.frame(ratio =  data_for_plot_w_0_75$X[ind,2]/
                         (data_for_plot_w_0_75$X[ind,2]+data_for_plot_w_0_75$X[ind,3]),
                       w = rep(0.75,100000))

ecdf_1 = ecdf(dt_comb_1$ratio)
ecdf_2 = ecdf(dt_comb_2$ratio)
ecdf_3 = ecdf(dt_comb_3$ratio)

library(latex2exp)
x_points = seq(0,1,0.0001)
pdf(file = '../../Results/DACOMP_ratio_approximation_fail.pdf',width = 7,height = 3)
par(mar=c(4.3, 3.75, 0.1, 0)+0.1)
plot(x_points,ecdf_1(x_points),ylim = c(0.6,1),lwd = 2,xlim = c(0,0.1),type = 'l',lty = 1,xlab = TeX('$ \\textbf{X}(\\textbf{e}_2) / (\\textbf{X}(\\textbf{e}_2)+ \\textbf{X}(\\textbf{e}_3) )$'),ylab = 'CDF')
lines(x_points,ecdf_2(x_points),type = 'l',lty = 2,lwd = 2)
lines(x_points,ecdf_3(x_points),type = 'l',lty = 3,lwd = 2)
legend(x = 0.06,y = 0.85,legend = c('w=0','w=0.5','w=0.75'),lty = c(1,2,3),lwd = 2)
dev.off()

