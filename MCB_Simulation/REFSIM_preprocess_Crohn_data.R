# This script is used for preprocessing sOTU counts for the crohn dataset into a data structure used for simulation
# This is data taken from the healthy group, used for the analysis in the third section of van de putte et al.

#Load data
library(biomformat)
file_path <- "./gut_otu_table_sim.RData"
load(file_path)

#A vector of cell counts, this is taken from the supplementry data for the paper).
Average_Cell_Count = c(99183906832,
                       217221266674,
                       151329659639,
                       90299128192,
                       95809600876,
                       132035763643,
                       39495621712,
                       119271846132,
                       138185802083,
                       121428187249,
                       145939517130,
                       224754815252,
                       72756226469,
                       35927546063,
                       107628776723,
                       130776743130,
                       41841212649,
                       22764340926,
                       172987119956,
                       46425487958,
                       153062883269,
                       95590999034,
                       208495566551,
                       53070915155,
                       18359008555,
                       70503715565,
                       58500810719,
                       89307497088,
                       147599455753,
                       19187237081,
                       217276666820,
                       253525709584,
                       25395316913,
                       22507873775,
                       142169041331,
                       116275992152,
                       142346591541,
                       161477761409,
                       5853245978,
                       47942791702,
                       28769973004,
                       168930811359,
                       126572467016,
                       117499264965,
                       116251948461,
                       177107799526,
                       74197893718,
                       48330742911,
                       203676977127,
                       187211301217,
                       174670957469,
                       75065628771,
                       54160688666,
                       165005851861,
                       215345198361,
                       150051113739,
                       191391240184,
                       20159773225,
                       38956568917,
                       123582072411,
                       98056602368,
                       184826494881,
                       54837809766,
                       70590688175,
                       132689255189,
                       98343121358)
  
# flip the OTU table, reorder rows
otu_table = t(otu_table)
reordering_permutation = order(rownames(otu_table))
otu_table = otu_table[reordering_permutation,]
rownames(otu_table) #it is now ordered

X = otu_table
Average_Cell_Count = Average_Cell_Count

# We keep only taxa that appeared in at least 4 subject, 4 being a threshold chosen for the simulation
prevalence_matrix = 1*(X>0)
total_counts_in_taxa = as.numeric(apply(prevalence_matrix,2,sum))
threshold = 4
to_keep = which(total_counts_in_taxa >= threshold) 
X2 = X[,to_keep]

# Pack as an object and save
Gut_Flow_data = list()
Gut_Flow_data$counts_matrix = X2
Gut_Flow_data$Flow_counts = Average_Cell_Count

save(Gut_Flow_data,file = 'Gut_Flow_data.RData')

