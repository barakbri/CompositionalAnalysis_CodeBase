# "Testing for differential abundance in compositional counts data, with application to microbiome studies".

## Introduction
This is the repository for reproducing the simulation results, real data analyses and graphs found in the manuscript with the above name [arXiv link](https://arxiv.org/pdf/1904.08937.pdf). The manuscript describes the DACOMP method ([GitHub link](https://github.com/barakbri/dacomp)), a non parametric approach for testing differential abundance in compositional counts data. The method aims for detecting differentially abundant taxa in 16S microbiome counts data.

The structure of the repository is as follows. The 'packages' subfolder contains the installation files for packages used for running simulations and data analysis. These packages need to be installed prior to running simulations and data analyses. The directory contains a list of the statistical methods used from each package, when reproducing results for the paper. See the Installtion section below for further details. The 'Scripts' subfolder contains the R scripts for simulation and real data analyses. The 'Crohn' folder contains the  Crohn Disease dataset analyzed in Section 5 of the manuscript. The 'HMP' folder contains the  metadata and OTU counts for the samples measured under the Human Microbiome Project. A results folder will be constructed in parallel to the repository directory. This directory will contain all results from real data analyses and simulations.

## Installation

1. Begin by downloading the repository and placing it in a directory. 
2. Open RStudio/R and install the packages found under the 'packages' subdirectory.
3. Unzip the file HMP/OTU_Counts.zip to its current directory. This file contains the OTU counts in CSV format (obtained from according to the instructions found in the manuscript from the HMP site).
4. Open the project file "MCB_Simulation.Rproj" found under scripts in RStudio. All scripts will be run under this project.
5. Run the script "Generate_Results_Dir.r". This script will generate the results directory, where all results will be placed.


## Scripts to run for reproducing results

The next subsection goes over the scripts that needed to be run to reproduce results from the paper. Scripts that are run by the user will have "RUN TO..." (capitalized) under their description. In Addition, scripts that are run below each subsubsection will have a number (**#1,#2**) at the start of their description to indicate the running order. Scripts called by other scripts,e.g. library scripts, will not have "RUN TO..." at the start of their description.

The scripts make use of two RData objects of sOTU reads obtained from the deblurr package. The deblurr package is not part of the R codebase, and was called manually using the CLI supplied by the package. See the manuscript on instructions on how and where the data for "deblurring" was obtained.
    
### Simulation
   * `REFSIM_MultipleScenario_batch.R` - **(#1)** RUN TO perform the simulations desribed under Section 4 of the   manuscript, and appendix A (additional simulations). The script performs multiple calls to `REFSIM_SingleScenario.R`. Note that this script took ~2 days to run on a 96 core Amazon C5 machine. It is not feasible to run this script on a PC.
   * `REFSIM_SingleScenario.R` -  The script for performing simulations, for a single scenario found in Section 4, for all DACOMP variants and competitors. Saves output to the results directory.
   * `REFSIM_Analyze_Results.R`- **(#2)** RUN TO compile results from all scenarios, after running  `REFSIM_MultipleScenario_batch.R`.
   * `REFSIM_additional_competitors.R` - **(#3)** An additional script used for running more competitor method, for late edition to the manuscript. Currently runs ZINB-WAVE combined with Deseq2. The script calls functions from `zinbwave_imported_functions` (taken originally from https://github.com/mcalgaro93/sc2meta , see additional notes below and in file.)
   * `REFSIM_add_reference_selection_by_nr_counts.R` - **(#4)** An additional script used for running DACOMP with references such that $S_{crit}$ is set to have the lowest value, so that at least 100 counts are available in the reference taxa across all samples. Also tracks for the number of differentially abundant taxa entering the reference set. 
   * `validation_procedure_sim.R` - **(#5)** Script used for running the sim found in SM S6.2, checking the risk of differentially abundant taxa enter the reference set, whether the diagnostic check can detect the differentially abundant taxa, and whether testing can be done on a shrunken dataset without the differentially abundant taxa.
   * `REFSIM_Analyze_chance_of_bad_taxa.R` - **(#6)** RUN TO generate the table found in SM S6.4, describing the mean number of differentially abundant taxa that have erronuosly entered the reference set. Script is run after 1-5. Also - generates the results table for SM S6.2
   * `REFSIM_Analyze_Results_Tables.R` - **(#7)** RUN TO compile graphs and tables for Section 6 and the Supplementary Material, after running steps 1-4.

### Analyzing data from the Human Microbiome Project
   * `RDE_HMP_main.R` - **(#8)** RUN TO analyze results for the HMP dataset. This script calls `RDE_HMP_two_site_comparison.R` for each pair of body sites, and then outputs the results to a subdirectory under the results directory. See `RDE_Generate_Latex_tables.R` for processing outputs of this data example to a lateX table.
   * `RDE_HMP_two_site_comparison.R` - Analyze differential abundance of taxa between a pair of body sites. This script is called repeatedly by `RDE_HMP_main.R`, each time with a different pair of body sites as parameters.
      
### Data analysis for the Crohn Disease study
   * `RDE_Crohn.R / RDE_Crohn_Multivariate.R` - **(#9)** RUN TO analyze the crohn disease dataset. the two scripts perform the univariate tests for differential abundance, along with the multivariate tests for differential abundance, testing if sOTUs in each genera maintained their respective ratios in the presence of CD.
   * `RDE_Generate_Latex_tables.R` - **(#10)** RUN TO compile the tables generated by `RDE_HMP_main.R` and `RDE_Crohn.R` to lateX format.
   * `RDE_Analysis_Dilution_Example.R` - **(#11)** RUN TO analyze the dilution experiment example discussed in the Supplementary Material.
   * `RDE_Crohn_Multivariate_PCOA_functions.R` - A function script called by `RDE_Crohn_Multivariate.R` in order to produce "within-genus" PCoA plots for the data. 
   * `RDE_Crohn_validation_procedure_and_alternative_reference_set.R` - **(#12)** Additional script used for (a) checking the validity of the reference script, and (b) running DACOMP with an alternative reference set, chosen also using flow-cytometric measurements.

### Simulation results - in appendicies
   * `APPENDIX_SIM_ANCOM_NULL_CASES.R` - **(#13)** RUN TO generate the results found in appendix C - the simulation study for FWER of ANCOM under the global null.
   
   * `REFSIM_compute_FDR_random_selection_batch_script.R` - **(#14)** RUN TO obtain the results found in SM S4.3 (FDR of DACOMP with naive reference selection methods) Results are outputd to text files, see additional details in the script `REFSIM_compute_FDR_random_selection.R`.
   * `REFSIM_compute_FDR_random_selection.R` - A script called by `REFSIM_compute_FDR_random_selection_batch_script.R`, for either estimating FDR of naive reference selection methods or estimationing the T1E of the RVP, each time for a single scenario.

## Description of additional files:

### Scripts for producing graphs
* `Graph_section_2_1_incorrect_testing.R` - **(#15)** RUN TO generate the graph for the CDF of the test statistics commonly used for testing differential abundance, found under subsection 2.1 of the manuscript.
* `PLOT_NO_REMOVE_OBS.R`  - **(#16)** RUN to generate the graph found under Section 3, showing why removing samples based on the number of counts in a subvector is invalid.
* `ReferenceScore_Plot.R` - A function used for plotting the distribution of reference scores (as defined in subsection 3.1). The DACOMP package has a fancier version of this function, the plain version was used when producing the manuscript.
* `Distribution_of_S_scores.R` -  **(#17)** Used for plotting the histogram of reference scores for both diff-abundant and non-diff-abundant taxa. The plots are featured in SM S6.1.
* `Incorrect_Testing_Ratio.R` -  **(#18)** Used for Example 2, showing DACOMP-ratio may not be exactly valid.

### Data preprocessing
* `gut_otu_table_sim.RData` - A table of deblurred sOTU counts for the healthy cohort from the Crohn Disease study, used for simulations.
* `REFSIM_preprocess_Crohn_data.R` - A script for reading the table of deblurred sOTU reads (`gut_otu_table_sim.RData`) and combining it with the flow cytometry measurements of the healthy cohort to create a single object loaded for the data based simulations of Subsection 4.1 of the manuscript.
* `Gut_Flow_data.RData` - A single object for containing both the table of sOTU reads and flow cytometru measurements used for simulations.


### Scripts used for competitor methods - Wilcoxon, and Fisher exact test
* `exactHyperGeometricTest.R` - An implementation of the HG test described in Section 4.
* `Wilcoxon_TaxaWise.R` - A wrapper for using Wilcoxon rank sum tests (from the subzero auxiliary package) on each taxon in a given data set (marginal tests for equality of distributions).
* `zinbwave_imported_functions` - A script with additional functions used for running ZINB-WAVE+DESEQ2 as in Calgaro et al. 2020. The functions were taken from the GitHub repository https://github.com/mcalgaro93/sc2meta showing how to use Deseq2 + Zinb-Wave for Microbiome data. See additional comments in file.

### Scripts used for data generation in simulation
* `REFSIM_GenerateSettings_Index.R` - A single script generating the index for all simulatin scenarios and loading all additional scripts used for simulations.
* `REFSIM_Gut_Scenario.R` - Script for generating datasets for the scenarios of 4.1 - scenarios based on sampling from real data.
* `REFSIM_NO_OVERDISPERSION_SCENARIOS.R` - Scripts for generating simulated datasets for the scenarios of 4.2 - multinomial cases where the most abundant taxa are not differentially abundant.
* `REFSIM_NoCompositionality_Scenario.R` - Scripts for generating datasets with no compositionality - the scenarios of 4.3.



