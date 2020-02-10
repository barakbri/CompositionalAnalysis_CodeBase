
memory.limit(25000) # scale memory
library(biomformat)
Q_LEVEL = 0.1 #BH threshold level

file_path <- "../Crohn/gut_otu_table.RData" #load counts data
load(file_path)
Average_Cell_Count = read.csv(file = '../Crohn/FlowCytometryCounts.csv')$Average_Cell_Count #load flow cytometry counts
Y = rep(0,95) #group labeling is such that the first 29 subjects are CD, the rest are healthy.
Y[1:29] = 1

#show a box plot for the distribution of average cell counts across subjects
boxplot(Average_Cell_Count[Y==1],Average_Cell_Count[Y==0],names = c('Sick','Healthy'),main='Boxplot for cell counts in Healthy/Sick')
mean(Average_Cell_Count[Y==0])/(mean(Average_Cell_Count[Y==1])) #ratio between mean number of reads - discussed in paper
#2.508608
median(Average_Cell_Count[Y==0])/(median(Average_Cell_Count[Y==1])) # ratio between median number of reads
#3.086826

# order OTU table (from deblur) by original sample names
otu_table = t(otu_table)
reordering_permutation = order(rownames(otu_table))
otu_table = otu_table[reordering_permutation,]
rownames(otu_table) #it is now ordered

dim(otu_table)
X = otu_table