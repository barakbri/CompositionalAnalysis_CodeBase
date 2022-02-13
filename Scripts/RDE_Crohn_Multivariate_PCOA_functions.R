get_PCoA_for_single_genus = function(plotted_genus_id = 1,
                                     method = 'unifrac',
                                     remove_uni_discoveries = F,
                                     do_not_normalize = F,
                                     ind_reference = reference_taxa,
                                     col_of_ones_order = TAXON_PREVALENCE_RANK_TO_PICK,
                                     return_D_instead = F,
                                     do_not_include_ref_in_unnormalized = F){
  
  
  taxonomy_ind_for_genera = which(genera_labels == genera_labels_to_test[plotted_genus_id])
  taxonomy_ind_for_genera_before_filter = taxonomy_ind_for_genera
  if(remove_uni_discoveries){
    taxonomy_ind_for_genera = taxonomy_ind_for_genera[!(taxonomy_ind_for_genera %in% discovered_by_at_least_one)]
  }
  if(method == 'unifrac'){
    if(do_not_include_ref_in_unnormalized){
      stop(' Unifrac and do_not_include_ref_in_unnormalized mutually exclusive')
    }
    ind_to_take_as_reference   = ind_reference
    ind_to_take_as_reference = ind_to_take_as_reference[!(ind_to_take_as_reference %in% taxonomy_ind_for_genera_before_filter)]
    ranks_of_prevalence = rank(prevalence[ind_to_take_as_reference],ties.method = 'first')
    ind_to_take_as_reference = ind_to_take_as_reference[which(ranks_of_prevalence == length(ranks_of_prevalence)+1-col_of_ones_order)]
    ind_reference = ind_reference[ !(ind_reference %in% taxonomy_ind_for_genera)]
    
    combined_ind = c(taxonomy_ind_for_genera,ind_reference)
    if(!do_not_normalize){
      X_current = (X[,combined_ind])
    }else{
      X_current = (X_rarefied[,combined_ind])
    }
   
    
    if(!do_not_normalize){
      if(!is.null(Subvectors_rarefied[[plotted_genus_id]])  & !remove_uni_discoveries){
        X_current = Subvectors_rarefied[[plotted_genus_id]]
      }else{
        lambda = min(apply(X_current,1,sum))
        X_current = vegan::rrarefy(X_current,sample = lambda)  
      }
    }    
    
    X_current = X_current[,1:length(taxonomy_ind_for_genera),drop=F]
    X_current = cbind(X_current,X[,ind_to_take_as_reference,drop=F])
    X_current[,ncol(X_current)] = 1
    
    rows_permutation = sample(1:nrow(X_current))
    X_current = t(X_current)
    #rownames(X_current)
    colnames(X_current) <- paste0("Sample", 1:ncol(X_current))
    
    taxmat = unname(taxa)[c(taxonomy_ind_for_genera,ind_to_take_as_reference),]
    rownames(taxmat) <- rownames(X_current)
    colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    
    OTU = otu_table(X_current, taxa_are_rows = TRUE)
    TAX = tax_table(taxmat)
    
    sampledata = sample_data(data.frame(
      CD = Y,
      row.names=colnames(OTU),
      stringsAsFactors=FALSE
    ))
    
    physeq = phyloseq(OTU, TAX,sampledata,phy_tree(fitGTR$tree))
    
    d = UniFrac(physeq,weighted = F)
    
  }
  
  if(method %in%  c('Rmahal','BC','L1')){
    ind_reference = ind_reference[!(ind_reference %in% taxonomy_ind_for_genera_before_filter)]
    if(!do_not_normalize){
      total_counts_in_reference = apply(X[,ind_reference],1,sum)
      X_current   = cbind(X[,taxonomy_ind_for_genera],total_counts_in_reference)
      lambda_j = min(apply(X_genera_to_rarefy,1,sum))
      X_current = vegan::rrarefy(X_current,sample = lambda_j)  
    }else{
      if(do_not_include_ref_in_unnormalized){
        X_current   = cbind(X_rarefied[,taxonomy_ind_for_genera])
      }else{
        total_counts_in_reference = apply(X_rarefied[,ind_reference],1,sum)
        X_current   = cbind(X_rarefied[,taxonomy_ind_for_genera],total_counts_in_reference)
      }
        
    }
    if(method == 'Rmahal'){
      X_current = X_current[,-ncol(X_current),drop=F]
      
      if(!do_not_normalize & !is.null( DIST_rarefied[[plotted_genus_id]] ) &  !remove_uni_discoveries){ #best to use previous rarefaction if possible
        d = DIST_rarefied[[plotted_genus_id]] 
      }else{
        d = smahal(X_current)
      }
    }
    if(method == 'BC'){
      if(!do_not_normalize & !remove_uni_discoveries)  #best to use previous rarefaction if possible
        X_current = Subvectors_rarefied[[plotted_genus_id]]
      d = vegan::vegdist(X_current)
    }
    if(method == 'L1'){
      if(!do_not_normalize & !remove_uni_discoveries)  #best to use previous rarefaction if possible
        X_current = Subvectors_rarefied[[plotted_genus_id]]
      d = dist(X_current,method = 'manhattan')
    }
  }
  if(return_D_instead){
    return(d)
  }
  PCoA_obj = pcoa(d)
  return(PCoA_obj)
}

library(phyloseq)


run_configuration_with_raw = function(method = 'unifrac',
                             remove_uni_discoveries = F,plot_with_ranks = T,
                             disable_numbers_in_ratio_plot = F,modes = c(1:5),do_not_include_ref_in_unnormalized = F){
  if(!(method %in% c('unifrac','Rmahal', 'BC' ,'ratio','L1'))){
    stop("method must be in 'unifrac','Rmahal', 'BC' ,'ratio','L1'")
  }
  ret = list()
  for(mode in c(1:5)){
    if(!(mode %in% modes)){
      next
    }
    disc = 'MISSING'
    if(mode==1){
      id_list = rarefaction_Mult_alone
      disc = 'Mult'
      
    }
    if(mode==2){
      id_list = rarefaction_Uni_alone
      disc = 'Uni'
      
    }
    if(mode==3){
      id_list = rarefaction_Uni_and_Mult
      disc = 'Both'
      
    }
    if(mode==4){
      id_list = rarefaction_None
      disc = 'None'
      
    }
    if(mode==5){
      id_list = sort(c(rarefaction_Mult_alone,rarefaction_Uni_alone,rarefaction_Uni_and_Mult))
      disc = 'Discoveries'
    }
    unnormalized_label = '_Combined'
    excluded_label = ''
    if(remove_uni_discoveries)
      excluded_label = '_filtered'
    ranks_label = ''
    pdf(file = paste0('../../Results/Crohn_PCoA_',method,'_',disc,'_',Agg_level,unnormalized_label,excluded_label,ranks_label,'.pdf'),width = 6,height = 6)
    
    if(method %in% c('BC','Rmahal','unifrac','L1')){
      PCOA_list = list()
      if(plot_with_ranks)
        par(mfrow=c(2,2))
      else
        par(mfrow(c(1,2)))
      for(i in 1:length(id_list)){
        try({
          D= get_PCoA_for_single_genus(plotted_genus_id = id_list[[i]],
                                       method = method,
                                       remove_uni_discoveries = remove_uni_discoveries,
                                       do_not_normalize = F,return_D_instead = T,do_not_include_ref_in_unnormalized = do_not_include_ref_in_unnormalized);
          D_raw= get_PCoA_for_single_genus(plotted_genus_id = id_list[[i]],
                                       method = method,
                                       remove_uni_discoveries = remove_uni_discoveries,
                                       do_not_normalize = T,return_D_instead = T,do_not_include_ref_in_unnormalized = do_not_include_ref_in_unnormalized);
          perm_vec = 1:nrow(X)#sample(1:nrow(X))
          PCoA_obj =  pcoa(as.matrix(D)[perm_vec,perm_vec])
          PCoA_obj_raw = pcoa(as.matrix(D_raw)[perm_vec,perm_vec])
          X_1 = PCoA_obj$vectors[,1]
          Y_1 = PCoA_obj$vectors[,2]
          X_2 = PCoA_obj_raw$vectors[,1]
          Y_2 = PCoA_obj_raw$vectors[,2]
          
          #jitter
          X_1 = jitter(X_1,amount = mean(dist(X_1))/20,factor = 5)
          Y_1 = jitter(Y_1,amount = mean(dist(Y_1))/20,factor = 5)
          X_2 = jitter(Y_1,amount = mean(dist(Y_1))/20,factor = 5)
          Y_2 = jitter(Y_2,amount = mean(dist(Y_2))/20,factor = 5)
          
          X_1_r = rank(X_1,ties.method = 'average')
          Y_1_r = rank(Y_1,ties.method = 'average')
          X_2_r = rank(X_2,ties.method = 'average')
          Y_2_r = rank(Y_2,ties.method = 'average')
      
          
          plot(X_1, 
               Y_1, 
               col = Y[perm_vec]+1, pch='+', cex = 0.8, xlab = 'PCOA axis 1', ylab = 'PCOA axis 2',main = paste0(id_list[[i]],')',genera_labels_to_test[id_list[[i]]]) )
          abline(h=0,col='blue',lty=2)
          abline(v=0,col='blue',lty=2)
          plot(X_2, 
               Y_2, 
               col = Y[perm_vec]+1, pch='+', cex = 0.8, xlab = 'PCOA axis 1', ylab = 'PCOA axis 2',main = paste0('unnormalized:') )
          abline(h=0,col='blue',lty=2)
          abline(v=0,col='blue',lty=2)
          
          if(plot_with_ranks){
            plot(X_1_r, 
                 Y_1_r, 
                 col = Y[perm_vec]+1, pch='+', cex = 0.8, xlab = 'ranks of PCOA axis 1', ylab = 'ranks of PCOA axis 2',main = 'ranked:' )
            
            plot(X_2_r, 
                 Y_2_r, 
                 col = Y[perm_vec]+1, pch='+', cex = 0.8, xlab = 'ranks of PCOA axis 1', ylab = 'ranks of PCOA axis 2',main = paste0(' unnormalized,ranked:') )
            
          }
         
          current_PCoA_obj = list(D = D, D_raw = D_raw,
                                  X_1 = X_1,Y_1 = Y_1,X_2 = X_2,Y_2 = Y_2,
                                  X_1_r = X_1_r,Y_1_r = Y_1_r,X_2_r = X_2_r,Y_2_r = Y_2_r,
                                  perm_vec = perm_vec, PCoA_obj = PCoA_obj, PCoA_obj_raw = PCoA_obj_raw)
          ret[[id_list[[i]]]] = current_PCoA_obj
        })
      }
      
      par(mfcol=c(1,1))
    }
    
    if(method %in% c('ratio')){
      par(mfcol=c(3,2))
      for(i in 1:length(id_list)){
        X_to_plot = Subvectors_TSS[[id_list[[i]]]]
        X_to_plot = X_to_plot[,-ncol(X_to_plot),drop=F]
        relative_sum_of_group_vs_ref = apply(X_to_plot,1,sum)
        l = list(CD = relative_sum_of_group_vs_ref[Y==1],
                 H = relative_sum_of_group_vs_ref[Y==0])
        if(disable_numbers_in_ratio_plot){
          current_main = paste0(genera_labels_to_test[id_list[[i]]])
        }else{
          current_main = paste0(id_list[[i]],')',genera_labels_to_test[id_list[[i]]])
        }
        boxplot(l,main = current_main)
        ret[[id_list[[i]]]] = l
      }
      
      par(mfcol=c(1,1))
    }
    #par(mfrow=c(1,1))
    
    dev.off()
  }
  return(ret)
}