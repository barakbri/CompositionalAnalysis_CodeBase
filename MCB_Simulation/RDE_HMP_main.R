
Valid_Body_Sites = c(
  "Stool",
  
  "Saliva",
  "Tongue_dorsum",
  "Hard_palate",
  "Buccal_mucosa",
  "Attached_Keratinized_gingiva",
  "Palatine_Tonsils",
  "Throat",
  "Supragingival_plaque",
  "Subgingival_plaque",
  
  "Right_Antecubital_fossa", #elbow pit
  "Left_Retroauricular_crease", #Ear
  "Right_Retroauricular_crease", #Ear
  "Left_Antecubital_fossa", #elbow pit
  "Anterior_nares", # nose- nostrils
  
  "Vaginal_introitus",
  "Posterior_fornix",
  "Mid_vagina"
  
)


region = c(1,rep(2,9),rep(3,5),rep(4,3))

# length(Valid_Body_Sites); choose(length(Valid_Body_Sites),2)
# choose(9,2) + choose(5,2) + choose(3,2)

MODE_RUN = T
MODE_ANALYZE = T
START_FROM = 1
START_FROM_SITE_2 = 2
DATA_LOADED = F
DS.FDR.LIMIT = 10
DS.FORCE = T

if(MODE_RUN){
  for(site1 in START_FROM:(length(Valid_Body_Sites)-1)){
    min_val_site2 = (site1+1)
    if(site1 == START_FROM)
      min_val_site2 = START_FROM_SITE_2
    for(site2 in (min_val_site2):length(Valid_Body_Sites)){
      if(region[site1]!=region[site2])
        next
      
      Start_time  = Sys.time()
      print(paste0('########################################'))
      print(paste0('Meta running script for ',site1,'-',site2,', ',Valid_Body_Sites[site1]," - ",Valid_Body_Sites[site2]))
      print(paste0('########################################'))
      
      
      Y0 = Valid_Body_Sites[site1]
      Y1 = Valid_Body_Sites[site2]
      if(DATA_LOADED){
        LOAD_DATA = F
      }else{
        LOAD_DATA = T
      }
      
      source('C:/MCB2/MCB2/MCB_Simulation/HMP_RDE_V12.R', echo=FALSE)
      
      DATA_LOADED = T
      
      End_time  = Sys.time()
      cat(paste0('Total time for two-site comparison: '))
      print(End_time - Start_time)
    }
  }
}

# 
# util_name_from_bodysite = function(X){
#   return(strsplit(X,'UBERON:')[[1]][2])
# }


if(MODE_ANALYZE){
  
  row_pointer = 1
  nr_rows_in_results = choose(length(Valid_Body_Sites),2)
  dt_results = data.frame(Same_Region = rep(NA,nr_rows_in_results),
                          Site1 = rep(NA,nr_rows_in_results), Site2 = rep(NA,nr_rows_in_results),
                              nr_subjects = rep(NA,nr_rows_in_results), taxa_for_test = rep(NA,nr_rows_in_results),
                              Nr_Wilcoxon_Rejections = rep(NA,nr_rows_in_results),
                              Nr_Wilcoxon_Rejections_Normalized = rep(NA,nr_rows_in_results),
                              Nr_Wilcoxon_Rejections_Normalized_0_75 = rep(NA,nr_rows_in_results),
                              ANCOM_rejections= rep(NA,nr_rows_in_results), 
                              ref_size_selected = rep(NA,nr_rows_in_results),
                              rejections = rep(NA,nr_rows_in_results),
                              shared = rep(NA,nr_rows_in_results),
                              rejected_in_reference = rep(NA,nr_rows_in_results),
                              median_lambda = rep(NA,nr_rows_in_results),
                              ds_fdr_total_rejections = rep(NA,nr_rows_in_results),
                              ds_fdr_shared = rep(NA,nr_rows_in_results),
                              ds_fdr_unique = rep(NA,nr_rows_in_results),
                              ds_fdr_ref_rejected = rep(NA,nr_rows_in_results),
                              lambda_multiplier = rep(NA,nr_rows_in_results),
                              subset_disc = rep(NA,nr_rows_in_results),
                              Abundance_shared = rep(NA,nr_rows_in_results),
                              Abundance_ANCOM_Unique = rep(NA,nr_rows_in_results),
                              Abundance_ANCOM_in_ref = rep(NA,nr_rows_in_results),
                          
                              Prevalence_shared = rep(NA,nr_rows_in_results),
                              Prevalence_ANCOM_Unique = rep(NA,nr_rows_in_results),
                              Prevalence_ANCOM_in_ref = rep(NA,nr_rows_in_results),
                          
                              Mult1_Rejections = rep(NA,nr_rows_in_results),
                              Mult1_Shared = rep(NA,nr_rows_in_results),
                              Mult1_Rej_in_Ref = rep(NA,nr_rows_in_results)
                              )
  
  for(site1 in 1:(length(Valid_Body_Sites)-1)){
    for(site2 in (site1+1):length(Valid_Body_Sites)){
      if(region[site1]!=region[site2])
        next
      
      Y0 = Valid_Body_Sites[site1]
      Y1 = Valid_Body_Sites[site2]
      results_save_file = paste0(HMP_RESULTS_DIR,"RESULTS_FILE_",(Y0),"_",(Y1),'.Rdata')    
      load(file = results_save_file)
      ref_size_selected = length(results_to_save$selected_references_obj$selected_references) #results_to_save$selected_references_obj$arg_obj
      #util_name_from_bodysite(results_to_save$body_sites$Y0)
      #util_name_from_bodysite(results_to_save$body_sites$Y1)
      #results_to_save$`dim(X)`[2] # OTUs before formatting
      #results_to_save$`dim(X)`[1] # Nr Subjects
      #results_to_save$`dim(X_2)`[2] #NR OTUS after formatting
      
      nr_vanilla_shared       = 0
      nr_vanilla_unique       = 0
      nr_vanilla_unique_ANCOM = 0
      median_lambda = NA
      
      rejected_dsfdr_shared = 0
      rejected_dsfdr_unique = 0
      rejected_dsfdr_ref = 0
      
      nr_rejections_in_ref = 0
      if(!is.na(ref_size_selected)){
        nr_vanilla_shared       = length(results_to_save$intersect_with_ANCOM$shared_with_ANCOM)
        nr_vanilla_unique       = length(results_to_save$intersect_with_ANCOM$unique)
        nr_vanilla_unique_ANCOM = length(results_to_save$intersect_with_ANCOM$unique_ANCOM)  
        nr_rejections_in_ref =  length(which(p.adjust(results_to_save$res_dsfdr$p.values.ref,method = 'BH')<=Q))
        
        #mult1_intersect_obj = intersect_with_ANCOM(results_to_save$res_dsfdr_lambda_1_0,results_to_save$ANCOM_res)
        mult1_rejected = which(p.adjust(results_to_save$res_dsfdr_lambda_1_0$p.values.test,method = 'BH')<=Q)
        if(length(mult1_rejected)>DS.FDR.LIMIT | DS.FORCE){
          mult1_rejected = results_to_save$res_dsfdr_lambda_1_0$rejected
        }
        
        ANCOM_DISC = c(results_to_save$intersect_with_ANCOM$shared_with_ANCOM,results_to_save$intersect_with_ANCOM$unique_ANCOM)
        mult1_nr_vanilla_shared       = sum(mult1_rejected %in% ANCOM_DISC)
        mult1_nr_rejections_in_ref = length(which(p.adjust(results_to_save$res_dsfdr_lambda_1_0$p.values.ref,method='BH')<=Q))#length(which(p.adjust(results_to_save$res_dsfdr_lambda_1_0$pvalue_ref,method = 'BH')<=Q))
        #which(p.adjust(results_to_save$res_dsfdr_lambda_1_0$p.values.ref,method='BH')<=0.05)
        
        ANCOM_rej = union(results_to_save$intersect_with_ANCOM$shared_with_ANCOM,
                          results_to_save$intersect_with_ANCOM$unique_ANCOM)
        
        median_lambda = median(results_to_save$res_dsfdr$min_value_array)
        if(is.null(median_lambda))
          median_lambda = NA
        rejected_dsfdr_shared = length(intersect(results_to_save$res_dsfdr$rejected,ANCOM_rej))
        rejected_dsfdr_unique = length(results_to_save$res_dsfdr$rejected) - rejected_dsfdr_shared 
        rejected_dsfdr_ref = length(which(p.adjust(results_to_save$res_dsfdr$p.values.ref,method='BH')<=Q))#length(results_to_save$res_dsfdr$rejected_ref)
        #results_to_save$res_dsfdr$p.values
        
      }
      
      
      nr_Wilcoxon_rejections = results_to_save$Nr_Wilcoxon_Rejections
      nr_Wilcoxon_rejections_norm = results_to_save$Nr_Wilcoxon_Normalize_Rejections
      nr_Wilcoxon_rejections_norm_0_75 = results_to_save$Nr_Wilcoxon_Normalize_0_75_Rejections
      
      same_region = 0
      if(region[site1] == region[site2])
        same_region = 1
      multiplier = results_to_save$lambda_multiplier
      if(is.null(multiplier))
        multiplier = NA
      rejections_subset = results_to_save$nr_rejections_subset
      if(is.null(rejections_subset))
        rejections_subset = NA
      row_to_insert = c(same_region,(Y0),(Y1),
                        results_to_save$`dim(X)`[1],results_to_save$`dim(X_2)`[2],
                        nr_Wilcoxon_rejections,
                        nr_Wilcoxon_rejections_norm,
                        nr_Wilcoxon_rejections_norm_0_75,
                        nr_vanilla_shared + nr_vanilla_unique_ANCOM,
                        ref_size_selected,
                        nr_vanilla_shared+nr_vanilla_unique, nr_vanilla_shared,
                        nr_rejections_in_ref,
                        median_lambda,
                        rejected_dsfdr_shared + rejected_dsfdr_unique,
                        rejected_dsfdr_shared,
                        rejected_dsfdr_unique,
                        rejected_dsfdr_ref,
                        multiplier,
                        rejections_subset,
                        results_to_save$abundance_shared_discoveries,
                        results_to_save$abundance_ANCOM_UNIQUE_discoveries,
                        results_to_save$abundance_ANCOM_in_ref,
                        
                        results_to_save$prevalence_shared_discoveries,
                        results_to_save$prevalence_ANCOM_UNIQUE_discoveries,
                        results_to_save$prevalence_ANCOM_in_ref,
                        
                        length(mult1_rejected),
                        mult1_nr_vanilla_shared,
                        mult1_nr_rejections_in_ref
                        )
      print(paste0('Inserting row ',row_pointer,' with length : ',length(row_to_insert)))
      dt_results[row_pointer,] = row_to_insert
      row_pointer = row_pointer + 1
      
    }
  }
  
  
  write.csv(dt_results[1:(row_pointer-1),],file = 'C:/MCB2/HMP2/Results/results.csv',quote = F,row.names = F)
  
  
}
