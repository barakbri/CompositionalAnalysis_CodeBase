library(xtable)

dt = read.csv('C:/MCB2/Results/gut_qPCR_shared_disc_mat.csv')
xtable::xtable(dt)


dt = read.csv('C:/MCB2/Results/Gut_Data_Results.csv')
names(dt)[c(2,4,9,12,13)]
xtable::xtable(dt[,c(2,4,9,12,13)])



dt = read.csv('C:/MCB2/HMP2/Results/results.csv')
names(dt)

Assumption_Test_Failed = 1*(dt$Mult1_Rej_in_Ref>0)
Assumption_Test_Failed[Assumption_Test_Failed == 1] = '*'
Assumption_Test_Failed[Assumption_Test_Failed == '0'] = ''

dt_export = data.frame(Site1 = dt$Site1,
Site2 = dt$Site2,
NR.Taxa = dt$taxa_for_test,
ANCOM = dt$ANCOM_rejections,
DISC_Wilcoxon_CSS = dt$Nr_Wilcoxon_Rejections_Normalized_0_75,
DISC_Wilcoxon_TSS = dt$Nr_Wilcoxon_Rejections_Normalized,
Discoveries = paste0(dt$Mult1_Rejections,Assumption_Test_Failed),
Shared = dt$Mult1_Shared,
ReferenceSize = dt$ref_size_selected)


  library(xtable)
print(xtable(dt_export), include.rownames=FALSE)

