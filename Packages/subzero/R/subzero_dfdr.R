disable.ties.correction = F

Compute.ChiSq.test = function(X_matrix,Y_matrix,lambda,statistic = 'Wilcoxon'){
  nr_subsamples = ncol(X_matrix)
  nr_bootstraps  = ncol(Y_matrix)
  stats = rep(0,nr_bootstraps)
  current_stats = rep(0,nr_bootstraps)
  n = nrow(X_matrix)
  for(s in 1:nr_subsamples){
    ranked_X = rank(X_matrix[,s], ties.method = 'average')
    if(statistic == 'Wilcoxon'){
      current_stats = (subzero::rcpp_Wilcoxon_PermTest_Given_Permutations(ranked_X,Y_matrix)[[1]])
      current_stats = current_stats - sum(Y_matrix[,1])/nrow(Y_matrix) * sum(ranked_X)
      
      NR_Y0 = sum(Y_matrix[,1] == 0)
      N = nrow(Y_matrix)
      NR_Y1 = N - NR_Y0
      E_H0 = sum(ranked_X) * NR_Y1 / N
      
      unique_value_counts = table(ranked_X)
      which_unique_value_counts_are_ties = which(unique_value_counts>1)
      unique_value_counts_only_ties_correction = 0
      if(length(which_unique_value_counts_are_ties)>0 & !disable.ties.correction){
        t_r = unique_value_counts[which_unique_value_counts_are_ties]
        unique_value_counts_only_ties_correction  = sum(t_r * (t_r^2 - 1))
      }
      
      V_H0 = NR_Y0 * NR_Y1 * (N+1) / 12 - NR_Y0 * NR_Y1 * unique_value_counts_only_ties_correction / (12 * N * (N - 1))
      current_stats = current_stats/sqrt(V_H0) #convert to Z score
    }
    if(statistic == 'Avg.Diff'){
      current_stats = rep(0,nr_bootstraps)
      for(b in 1:ncol(Y_matrix)){
        current_stats[b] = mean(X_matrix[Y_matrix[,b]==0,s]) - mean(X_matrix[Y_matrix[,b]==1,s])
      }
      #current_stats = current_stats^2
    } 
    if(statistic == 'Log.Avg.Diff'){
      current_stats = rep(0,nr_bootstraps)
      for(b in 1:ncol(Y_matrix)){
        current_stats[b] = mean(log10(X_matrix[Y_matrix[,b]==0,s]+1)) - mean(log10(X_matrix[Y_matrix[,b]==1,s]+1))
      }
      #current_stats = current_stats^2
    }   
    if(statistic == 'TwoPartWilcoxon')
      current_stats = (subzero::rcpp_TwoPartTest_Given_Permutations(ranked_X,Y_matrix)[[1]])
    if(statistic == 'SignedWilcoxon'){
      
      for(j in 1:nr_bootstraps){
        g1 = X_matrix[Y_matrix[1:(n/2),j],s]
        g2 = X_matrix[Y_matrix[(n/2 +1):(n),j],s]
        add_stat = SignedRankWilcoxon.statistic(g1,g2)  #(qnorm(wilcox.test(g1, g2, paired = TRUE,correct = F,exact = F,alternative = 'greater')$p.value))^2
        current_stats[j] = current_stats[j] + add_stat
      }
      
    }
    if(statistic == 'KW'){
      for(j in 1:nr_bootstraps){
        add_stat = (kruskal.test(X_matrix[,s], as.factor(Y_matrix[,j]))$statistic)
        if(is.nan(add_stat))
          add_stat = 0
        current_stats[j] = current_stats[j] + add_stat
      }
    }
    if(statistic %in% c('Deconv.KS','Deconv.CVM')){
      # if(verbose)
      #   print(paste0('Deconv on taxa ',i))
      for(j in 1:nr_bootstraps){
        temp = deconv_stat(Counts = X_matrix[,s],Labels = Y_matrix[,j],TotalCounts = rep(lambda,nrow(X_matrix)))
        if(statistic == 'Deconv.KS')
          add_stat = temp$stat.ks
        if(statistic == 'Deconv.CVM')
          add_stat = temp$stat.cvm
        
        if(is.nan(add_stat))
          add_stat = 0
        current_stats[j] = current_stats[j] + add_stat
      }
    }
    stats = stats + (current_stats)
  }
  stat = stats/nr_subsamples
  if(statistic %in% c('Wilcoxon','Avg.Diff','Log.Avg.Diff','SignedWilcoxon')){
    stats = stats^2
  }
  return(stats)
}

#subzero test function
subzero.dfdr.test = function(X,y,nr_reference_taxa,
                   verbose = T,
                   z = NULL,
                   test = 'Wilcoxon',
                   nr_perm = 200,q=0.05, disable_DS.FDR = F
                   ){
  nr_rarefactions_multiple_X = 1
  lambda_multiplier = 1
  if(nr_rarefactions_multiple_X != 1){
    stop("Multiple rarefactions currently not supported")
  }
  p = ncol(X)
  n = nrow(X)
  min_value_array = rep(NA,p)
  rarefied_data_length_array = rep(NA,p)
  pval_res = rep(NA,p)
  pval_res.CVM = rep(NA,p)
  taxa_nr_res = rep(NA,p)
  reference_values = rep(NA,n)
  lambda_selected = rep(NA,p)
  remaining_Y = list()
  # Compute reference values
  if(length(nr_reference_taxa)>1){
    for(i in 1:n){
      reference_values[i] = sum(X[i,nr_reference_taxa])
    }
  }else{
    reference_values = X[,nr_reference_taxa] 
  }
  
  stats_matrix = matrix(NA, ncol = p, nrow = nr_perm+1)
  rarefaction_matrix = matrix(NA,nrow = n,ncol = nr_rarefactions_multiple_X)
  
  #We compute the permutation matrix
  Y_matrix = matrix(NA, ncol = nr_perm+1, nrow = n)
  if(test == 'SignedWilcoxon'){
    Y_matrix[,1] = 1:n
    for( i in 2:ncol(Y_matrix)){
      #generate a permutation over signed subjects
      ind = 1:n
      for(t in 1:(n/2)){
        if(runif(1)<=0.5){
          ind[t] = t +(n/2)
          ind[t +(n/2)] = t
        }
      }
      Y_matrix[,i] = ind
    }
  }else{
    Y_matrix[,1] = y
    for( i in 2:ncol(Y_matrix)){
      Y_matrix[,i] = sample(y)
    }
  }
  
  
  
  
  
  for(i in 1:p){
    
    if(verbose)
      if(i%% ceiling(p/100) == 1)
        print(i)
    
    nom = X[,i]
    dnom = reference_values
    
    #
    if(i %in% nr_reference_taxa){
      nom = X[,i]
      dnom = reference_values - nom 
    }
      
    
    #else, we compute the test
    
    
    #choose for rarefaction
    total_reads_per_subject = nom+dnom
    
    min_value = ceiling( lambda_multiplier * min(total_reads_per_subject) )
    min_value_array[i] = min_value
    to_keep = which(total_reads_per_subject >= min_value)
    
    z_keep = z
    
    #enforcing pairs to have both items
    if(!is.null(z)){
      z_keep = z[to_keep]
      counts_for_subject = table(z_keep)
      z_to_keep = as.numeric(which(counts_for_subject == 2))
      z_values_to_keep = as.numeric(names(counts_for_subject)[z_to_keep])
      to_keep = which(z %in% z_values_to_keep) # we already filtered once by minimum value
      z_keep = z[to_keep]
    }
    
    nom_keep = nom[to_keep]
    dnom_keep = dnom[to_keep]
    y_keep = y[to_keep]
    remaining_Y[[i]] = y_keep
    rarefied_data_length_array[i] = length(nom_keep)
    
    check_1_minimum_size_group = T #true for check passed
    
    if(length(y_keep)<1){
      check_1_minimum_size_group = F
    }
    if(min(table(y_keep))<1 | (length(table(y_keep))!=length(table(y)))){check_1_minimum_size_group = F}          
    if(!check_1_minimum_size_group){next}
    
    # perform the actual subsample,
    nom_keep_original = nom_keep
    dnom_keep_original = dnom_keep
    
    if(check_1_minimum_size_group){
      
      #perform subsample
      for(s in 1:nr_rarefactions_multiple_X){
        temp_subsampled =  rhyper(n, nom_keep_original, dnom_keep_original,min_value)  
        rarefaction_matrix[,s] = temp_subsampled 
      }
      stats_matrix[,i] = Compute.ChiSq.test(rarefaction_matrix,Y_matrix,min_value,statistic = test)
      
      
      
    }
  }
  
  stats = stats_matrix
  p.values = rep(NA,ncol(stats))
  for(i in 1:ncol(stats)){
    p.values[i] = mean(stats[,i]>=stats[1,i])
  }
  
  if(!disable_DS.FDR){
    C_test = dfdr_find_thresholds(stats[,-nr_reference_taxa,drop=F],q)  
  }
  
  
  
  p.values.test = p.values; p.values.test[nr_reference_taxa] = NA
  p.values.ref = p.values; p.values.ref[-nr_reference_taxa] = NA
  
  ret = list()
  ret$rarefied_data_length_array = rarefied_data_length_array
  ret$min_value_array = min_value_array
  ret$lambda_selected = lambda_selected
  ret$remaining_Y = remaining_Y
  ret$stats_matrix = stats_matrix
  
  ret$p.values = p.values
  ret$p.values.test = p.values.test
  ret$p.values.ref = p.values.ref
  
  if(!disable_DS.FDR){
    ret$rejected = which(p.values.test<=C_test)
    ret$C_test = C_test  
  }
  return(ret)
}

dfdr_find_thresholds = function(stats,q=0.05){
  #convert to pvalues
  for(j in 1:ncol(stats)){
    stats[,j] = (nrow(stats) +1+1 - rank(stats[,j],ties.method = 'min'))/(nrow(stats)+1) # one +1 is for permutations, one +1 is because rank is 1 based
  }
  
  C_possible_values = sort(unique(stats[1,]))
  C_possible_values = C_possible_values[C_possible_values<=q]
  FDR_hat_vec = rep(NA,length(C_possible_values))
  for(c in 1:length(C_possible_values)){
    C = C_possible_values[c]
    Vhat = sum(stats<=C)/nrow(stats)
    Rhat = sum(stats[1,]<=C)
    FDR_hat = Vhat/Rhat
    FDR_hat_vec[c] = FDR_hat
  }
  
  selected_c_ind = (which(FDR_hat_vec<= q))
  selected_c = -1 #no rejections are allowed, this is in P-value "units"
  if(length(selected_c_ind)>0){
    selected_c = max(C_possible_values[selected_c_ind])
  }
  return(selected_c)
}
