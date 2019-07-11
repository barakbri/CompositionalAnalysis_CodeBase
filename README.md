# "Testing for differential abundance in compositional counts data, with application to microbiome studies".

## Introduction

## Installation

## Scripts to run for reproducing results

### Simulation
   * REFSIM_MultipleScenario_batch.R
   * REFSIM_SingleScenario.R
   * REFSIM_Analyze_Results.R
   * REFSIM_Analyze_chance_of_bad_taxa.R
   * REFSIM_Analyze_Results_Tables.R
   
### Data analysis for the crohn study
   * FlowCytometryCounts.csv
   * gut_otu_table.RData
   * RDE_Crohn.R
   * RDE_Generate_Latex_tables.R
   
   
### Analyzing data from the Human Microbiome Project
   * "RDE_HMP_main.R" 
   * "RDE_HMP_two_site_comparison.R" 
   
### Simulation results - in appendicies
   * `APPENDIX_SIM_ANCOM_NULL_CASES.R`
   * GlobalNull_Over_Reference.R
   
   * REFSIM_compute_FDR_random_selection.R
   * REFSIM_compute_FDR_random_selection_batch_script.R
   
## Description of additional files:

### Scripts for producing graphs
* Graph_section_2_1_incorrect_testing.R
* PLOT_NO_REMOVE_OBS.R
* ReferenceScore_Plot.R
* REFSIM_Plot_Score_Histo_For_Sims.R

### Data preprocessing
* gut_otu_table_sim.RData
* Gut_Flow_data.RData
* REFSIM_preprocess_Crohn_data.R

### Scripts used for competitor methods - Wilcoxon, and Fisher exact test
* exactHyperGeometricTest.R
* Wilcoxon_TaxaWise.R

### Scripts used for data generation in simulation
* REFSIM_GenerateSettings_Index.R
* REFSIM_Gut_Scenario.R
* REFSIM_NO_OVERDISPERSION_SCENARIOS.R
* REFSIM_NoCompositionality_Scenario.R



