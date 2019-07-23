# This script is used to analyze the HMP data set.
# script iterativaly calls RDE_HMP_two_site_comparison.R with different paris of body sites (from the same region)
# based on the partitioning found below.
# The script makes sure to load the data (which takes time) only for the first pairwise comparison.
# The second part of the scripts load the results from all analyzes and outputs the results


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

# code region for each body type
region = c(1,rep(2,9),rep(3,5),rep(4,3))


MODE_RUN = T #run pairwise comparisons
MODE_ANALYZE = T # should results be aggregated to a single file
START_FROM = 1 #set the starting points for the double loops
START_FROM_SITE_2 = 2
DATA_LOADED = F #has the data been loaded yet

# Run all pairwise comparisons
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
      if(DATA_LOADED){ # no need to load the data after the first time
        LOAD_DATA = F
      }else{
        LOAD_DATA = T
      }
      
      source('RDE_HMP_two_site_comparison.R', echo=FALSE)
      
      DATA_LOADED = T # data has been loaded, keep track for the next iteration
      
      End_time  = Sys.time()
      cat(paste0('Total time for two-site comparison: '))
      print(End_time - Start_time)
    }
  }
}

# aggregate all results to a single file
if(MODE_ANALYZE){
  
  row_pointer = 1
  nr_rows_in_results = choose(length(Valid_Body_Sites),2)
  #prepare data structure for results
  dt_results = data.frame(Same_Region = rep(NA,nr_rows_in_results), #indicator of two sites are in the same region
                          Site1 = rep(NA,nr_rows_in_results), Site2 = rep(NA,nr_rows_in_results), #names of sites
                              nr_subjects = rep(NA,nr_rows_in_results), taxa_for_test = rep(NA,nr_rows_in_results), #number of subjects and taxa for testing
                              Nr_Wilcoxon_Rejections = rep(NA,nr_rows_in_results), #number of rejection for unnormalized wilcoxon
                              Nr_Wilcoxon_Rejections_Normalized = rep(NA,nr_rows_in_results), #W-TSS
                              Nr_Wilcoxon_Rejections_Normalized_CSS = rep(NA,nr_rows_in_results), #W-CSS
                              ANCOM_rejections= rep(NA,nr_rows_in_results),  #number of rejections for ANCOM
                              ref_size_selected = rep(NA,nr_rows_in_results), #size of reference set selected, in number of taxa
                              rejections = rep(NA,nr_rows_in_results), #number of rejections by dacomp, with wilcoxon and subsampling
                              shared = rep(NA,nr_rows_in_results), #number of rejections DACOMP shared with ANCOM
                              median_lambda = rep(NA,nr_rows_in_results), #median lambda, across tested taxa
                              subset_disc = rep(NA,nr_rows_in_results), #number of discoveries for DACOMP, when using only 2/3 of the data (picked at random)
                              rejections_division = rep(NA,nr_rows_in_results), #number of discoveries, for DACOMP with Wilcoxon and normalization by division
                              rejections_Welch = rep(NA,nr_rows_in_results),  #number of discoveries, for DACOMP with Welch t test
                              rejections_division_Welch = rep(NA,nr_rows_in_results),#number of discoveries, for DACOMP with Welch t test and normalization by division
                              ALDEx2_Wilcoxon = rep(NA,nr_rows_in_results), #number of disceveries for ALDEx2 - two variants
                              ALDEx2_Welch = rep(NA,nr_rows_in_results),
                              Wrench = rep(NA,nr_rows_in_results), #number of discoveries for Wrench, with DESEQ2 tests
                              shared_ALDEx2_Wilcoxon = rep(NA,nr_rows_in_results), #number of discoveries shared between DACOMP (wilcoxon, rarefaction) and ALDEx2/ Wrench
                              shared_ALDEx2_Welch = rep(NA,nr_rows_in_results),
                              shared_Wrench = rep(NA,nr_rows_in_results)
                              )
  # iterate over pairs of body sites
  for(site1 in 1:(length(Valid_Body_Sites)-1)){
    for(site2 in (site1+1):length(Valid_Body_Sites)){
      if(region[site1]!=region[site2])
        next
      #load results
      Y0 = Valid_Body_Sites[site1]
      Y1 = Valid_Body_Sites[site2]
      results_save_file = paste0(HMP_RESULTS_DIR,"RESULTS_FILE_",(Y0),"_",(Y1),'.Rdata')    
      load(file = results_save_file)
      ref_size_selected = length(results_to_save$selected_references_obj$selected_references)
      
      #compute entries of the above table
      nr_vanilla_shared       = 0
      nr_vanilla_unique       = 0
      nr_vanilla_unique_ANCOM = 0
      median_lambda = NA
      if(!is.na(ref_size_selected)){
        nr_vanilla_shared       = length(results_to_save$intersect_with_ANCOM$shared_with_ANCOM)
        nr_vanilla_unique       = length(results_to_save$intersect_with_ANCOM$unique)
        nr_vanilla_unique_ANCOM = length(results_to_save$intersect_with_ANCOM$unique_ANCOM)  
        
        
        median_lambda = median(results_to_save$res_dsfdr$min_value_array)
        if(is.null(median_lambda))
          median_lambda = NA
        
      }
      
      nr_Wilcoxon_rejections = results_to_save$Nr_Wilcoxon_Rejections
      nr_Wilcoxon_rejections_norm = results_to_save$Nr_Wilcoxon_Normalize_Rejections
      nr_Wilcoxon_rejections_norm_CSS = results_to_save$Nr_Wilcoxon_Normalize_Rejections_CSS
      if(is.null(nr_Wilcoxon_rejections_norm_CSS)) #DEBUG
        nr_Wilcoxon_rejections_norm_CSS=NA 
      same_region = 0
      if(region[site1] == region[site2])
        same_region = 1
     
      #collect the row to insert
      row_to_insert = c(same_region,(Y0),(Y1),
                        results_to_save$`dim(X)`[1],results_to_save$`dim(X_2)`[2],
                        nr_Wilcoxon_rejections,
                        nr_Wilcoxon_rejections_norm,
                        nr_Wilcoxon_rejections_norm_CSS,
                        nr_vanilla_shared + nr_vanilla_unique_ANCOM,
                        ref_size_selected,
                        nr_vanilla_shared+nr_vanilla_unique, nr_vanilla_shared,
                        median_lambda,
                        results_to_save$nr_rejections_subset,
                        length(results_to_save$rejected_by_pval_division),
                        length(results_to_save$rejected_by_pval_Welch),
                        length(results_to_save$rejected_by_pval_division_Welch),
                        length(results_to_save$rejected.iqlr.wi),
                        length(results_to_save$rejected.iqlr.we),
                        length(results_to_save$wrench_rejected),
                        length(which(results_to_save$rejected.iqlr.wi %in% results_to_save$rejected_by_pval)),
                        length(which(results_to_save$rejected.iqlr.we %in% results_to_save$rejected_by_pval)),
                        length(which(results_to_save$wrench_rejected %in% results_to_save$rejected_by_pval))
                        )
      
      print(paste0('Inserting row ',row_pointer,' with length : ',length(row_to_insert)))
      dt_results[row_pointer,] = row_to_insert
      row_pointer = row_pointer + 1
    }
  }
  #write to file
  write.csv(dt_results[1:(row_pointer-1),],file = paste0('../../Results/HMP/results.csv'),quote = F,row.names = F)
  
}
