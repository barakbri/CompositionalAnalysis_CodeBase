# This script is used to reproduce Figure 1 on subsection 2.1
# This script compares several approaches for tests of differential abundance.

#Output file:
pdf(file = '../../Results/Incorrect_Tests.pdf', width = 7*1.5-1, height = 7/3*1.5 )
par(mfrow = c(1,3))
set.seed(1)
library(subzero)
#These are the effects for which we will compare, see paper for definition of the effect size
effect_grid = c(0.25, 0.33, 0.5) 

for(e_i in 1:3){ # For each effect, run sim:
  reps = 10000 #number of repetitions
  n_X = 50 #Sample sizes
  n_Y = 50
  N = 5000 #Sequencing depthj
  
  p_vector_X = c( 1- 6/N, 1/N, 5/N ) # baseline proportions
  effect = effect_grid[e_i]
  p_vector_Y = (1-effect) * p_vector_X + effect * c(1,0,0) #insert effect to first taxon alone
  #vector of P-values for each test variant
  pval.vec.PC = rep(NA,reps)
  pval.vec.NewDef = rep(NA,reps)
  #over repetitions:
  for(r in 1:reps){
    counts_X = rmultinom(n_X,N,prob = p_vector_X) #sample counts
    counts_Y = rmultinom(n_Y,N,prob = p_vector_Y)
    
    values_X = counts_X[2,]/max(counts_X[2,]+counts_X[3,],1) #values for test, for both samples, variant 1
    values_Y = counts_Y[2,]/max(counts_Y[2,]+counts_Y[3,],1)
    
    values_pseudocount_X = (counts_X[2,]+1)/(counts_X[2,]+counts_X[3,]+2)# values for test, variant 2
    values_pseudocount_Y = (counts_Y[2,]+1)/(counts_Y[2,]+counts_Y[3,]+2)
    
    #Run tests and collect results:
    X_wilcox_PC = c(values_pseudocount_X, values_pseudocount_Y)
    Y_wilcox = c(rep(0,n_X),rep(1,n_Y))
    X_wilcox_NewDef = c(values_X, values_Y)
    
    res_PC = subzero::PermSumCount.test(X = X_wilcox_PC, Y = Y_wilcox, B = 1000)
    pval.vec.PC[r] = res_PC$p.value  
    res_NewDef = subzero::PermSumCount.test(X = X_wilcox_NewDef, Y = Y_wilcox, B = 1000)
    pval.vec.NewDef[r] = res_NewDef$p.value
  }
  #plot graph:
  X = seq(0,1,0.001)
  ecdf_newdef = ecdf(pval.vec.NewDef)
  points_F_NewDef = (ecdf_newdef)(X)
  ecdf_PC = ecdf(pval.vec.PC)
  points_F_PC = (ecdf_PC)(X)
  
  plot(X, points_F_PC, type = 'l',main = paste0('w = ',effect),xlab = 'P-value',ylab = 'CDF',asp = 1)
  lines(X, points_F_NewDef, type = 'l',col = 'black',lty = 2)
  lines(c(0,1),c(0,1),lty = 1,col = 'gray')
  
}

#Close file:
dev.off()
par(mfrow = c(1,1))

effect_grid = c(0.25,0.33,0.5)
for(e_i in 1:length(effect_grid)){
  N = 5000 #Sequencing depth
  p_vector_X = c( 1- 6/N, 1/N, 5/N ) # baseline proportions
  effect = effect_grid[e_i]
  p_vector_Y = (1-effect) * p_vector_X + effect * c(1,0,0)
  print(paste0('effect: ',effect))
  print(paste0('Prob for zero in Y=0: ',round(dbinom(x = 0,size = N,prob = sum(p_vector_X[2:3])),3)))
  print(paste0('Prob for zero in Y=1: ',round(dbinom(x = 0,size = N,prob = sum(p_vector_Y[2:3])),3)))
}
