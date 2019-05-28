
pdf(file = '../../Results/Incorrect_Tests.pdf', width = 7*1.5-1, height = 7/3*1.5 )
par(mfrow = c(1,3))
set.seed(1)
library(subzero)
effect_grid = c(0.25, 0.33, 0.5)

for(e_i in 1:3){
  reps = 10000
  n_X = 50
  n_Y = 50
  N = 5000
  
  p_vector_X = c( 1- 6/N, 1/N, 5/N )
  effect = effect_grid[e_i]
  p_vector_Y = (1-effect) * p_vector_X + effect * c(1,0,0) 
  
  pval.vec.PC = rep(NA,reps)
  pval.vec.NewDef = rep(NA,reps)
  
  for(r in 1:reps){
    counts_X = rmultinom(n_X,N,prob = p_vector_X)
    counts_Y = rmultinom(n_Y,N,prob = p_vector_Y)
    
    values_X = counts_X[2,]/max(counts_X[2,]+counts_X[3,],1)
    values_Y = counts_Y[2,]/max(counts_Y[2,]+counts_Y[3,],1)
    
    values_pseudocount_X = (counts_X[2,]+1)/(counts_X[2,]+counts_X[3,]+2)
    values_pseudocount_Y = (counts_Y[2,]+1)/(counts_Y[2,]+counts_Y[3,]+2)
    
    X_wilcox_PC = c(values_pseudocount_X, values_pseudocount_Y)
    Y_wilcox = c(rep(0,n_X),rep(1,n_Y))
    X_wilcox_NewDef = c(values_X, values_Y)
    
    res_PC = subzero::PermSumCount.test(X = X_wilcox_PC, Y = Y_wilcox, B = 1000)
    pval.vec.PC[r] = res_PC$p.value  
    res_NewDef = subzero::PermSumCount.test(X = X_wilcox_NewDef, Y = Y_wilcox, B = 1000)
    pval.vec.NewDef[r] = res_NewDef$p.value
  }
  
  X = seq(0,1,0.001)
  ecdf_newdef = ecdf(pval.vec.NewDef)
  points_F_NewDef = (ecdf_newdef)(X)
  ecdf_PC = ecdf(pval.vec.PC)
  points_F_PC = (ecdf_PC)(X)
  
  plot(X, points_F_PC, type = 'l',main = paste0('w = ',effect),xlab = 'P-value',ylab = 'CDF',asp = 1)
  lines(X, points_F_NewDef, type = 'l',col = 'black',lty = 2)
  lines(c(0,1),c(0,1),lty = 1,col = 'gray')
  print(c(paste0('T1E, "+1" = ',round(mean(pval.vec.PC<=0.1),2)),
                                        paste0('T1E, \"0/0 -> 0\" = ',round(mean(pval.vec.NewDef<=0.1),2))
                                        ))
  #legend(0.25, 0.18, legend = c(paste0('T1E, "+1" = ',round(mean(pval.vec.PC<=0.1),2)),
  #                              paste0('T1E, \"0/0 -> 0\" = ',round(mean(pval.vec.NewDef<=0.1),2))
  #                              ), text.col = c('black','black'),lty = c(1,2))
  
}


dev.off()
par(mfrow = c(1,1))
