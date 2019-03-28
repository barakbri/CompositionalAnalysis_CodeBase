#install.packages('MCMCpack')

set.seed(1)
library(MCMCpack)
library(subzero)
alpha_vec = rep(1,200)
m = length(alpha_vec)
n = 50
prop_mat = rdirichlet(n,alpha = alpha_vec)
N = 5000

B_selected_original = list()
B_selected_permuted = list()
b_max = 1000 #per_core
j_to_flip = 1
Bscript = 1:m
Bscript_without_j = setdiff(Bscript,j_to_flip)
NR_cores = 8

generate_B = function(to_permute =F){
  B_selected = list()
  for(b in 1:b_max){
    print(paste0('Doing iter ',b))
    X = prop_mat
    if(to_permute){
      #print(1)
      B_total_proportion = apply(X[,Bscript],1,sum)
      #print(2)
      B_total_proportion_without_j = apply(X[,Bscript_without_j],1,sum)
      #print(3)
      rho = 1 - B_total_proportion_without_j/B_total_proportion
      #print(4)
      random_perm = sample(1:n)
      rho = rho[random_perm]
      
      temp = X[,Bscript_without_j]
      for(j in 1:nrow(temp)){
        temp[j,] = temp[j,]/sum(temp[j,])
      }
      print(5)
      X[,Bscript_without_j] = B_total_proportion*(1-rho)* temp
      X[,j_to_flip] = B_total_proportion*(rho)
    }
    
    for(i in 1:nrow(X)){
      X[i,] = rmultinom(1,size = N,prob = X[i,])
    }
    
    selected_ref = subzero::select.references.Median.SD.Threshold(X)
    B_selected[[b]]=selected_ref$selected_references
    
  }
  return(B_selected)
  
}

library(doRNG)
library(doParallel)
cl <- makeCluster(NR_cores)
registerDoParallel(cl)


B_selected_original = foreach(i=1:NR_cores, .options.RNG=1234,.combine = 'c') %dorng% {generate_B(to_permute = F) }

B_selected_permuted = foreach(i=1:NR_cores, .options.RNG=1234,.combine = 'c') %dorng% {generate_B(to_permute = T) }

stopCluster(cl)


#install.packages('hash')

aux_hash = function(vec){
  return(digest::digest(sort(vec)))
}

hash_original = lapply(B_selected_original,aux_hash)
hash_permuted = lapply(B_selected_permuted,aux_hash)


dt_o = as.data.frame(t(t(table(unlist(hash_original)))))
dt_p = as.data.frame(t(t(table(unlist(hash_permuted)))))

dt_merge = merge(x = dt_o, y = dt_p, by = "Var1", all = TRUE)[,c(1,3,5)]
names(dt_merge) = c('B-hashcode','Freq.Original','Freq.Permuted')

#dt_merge

dt_merge$Freq.Original[which(is.na(dt_merge$Freq.Original))] = 0
dt_merge$Freq.Permuted[which(is.na(dt_merge$Freq.Permuted))] = 0

plot(sqrt(dt_merge$Freq.Original),sqrt(dt_merge$Freq.Permuted),cex = 0.8,pch=20,xlab = 'Sqrt(Frequency of B, in original R)',ylab = 'Sqrt(Frequency of B, in permuted R)',asp=1)
abline(0,1,col='red')
plot((dt_merge$Freq.Original),(dt_merge$Freq.Permuted),cex = 0.8,pch=20,xlab = '(Frequency of B, in original R)',ylab = '(Frequency of B, in permuted R)',asp=1)
abline(0,1,col='red')
#(dt_merge$Freq.Original - dt_merge$Freq.Permuted)/sqrt


