# Packages used for simulation:

This directory stores the packages used for simulation:

* ancom_software  - a zip file of ancom.R and ancom.R_1.1-3.tar - two version of ANCOM for windows and linux.
* ALDEx2_1.16.0.tar - is used for the ALDEX variants
* dacomp_1.1.tar -  upgraded from 1.0 - no contains the welch test and normalization by ratio as well.
* metagenomeSeq_1.24.1.tar - used for W-CSS
* Wrench_1.2.0.tar and DESeq2_1.24.0.tar - used for running Wrench
* QMP-master (also unzipped to a directory) - the Flow cytometry method of Van de Putte et al. (2017) - given using their research paper.
* subzero (directory and  tar ball) - a small RCPP package used for a permutation based implementation of the Wilcoxon rank sum test. Used by the W-TSS,W-CSS and W-FLOW in simulation. called by Wilcoxon_TaxaWise.R.
