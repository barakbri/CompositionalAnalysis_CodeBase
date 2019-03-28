library(xtable)

dt = read.csv('../../Results/gut_qPCR_shared_disc_mat.csv')
xtable::xtable(dt)


dt = read.csv('../../Results/Gut_Data_Results.csv')
names(dt)[c(2,3,7,8,9)]
xtable::xtable(dt[,c(2,3,7,8,9)])



dt = read.csv('../../Results/HMP/results.csv')
names(dt)


dt_export = data.frame(Site1 = dt$Site1,
Site2 = dt$Site2,
NR.Taxa = dt$taxa_for_test,
ANCOM = dt$ANCOM_rejections,
DISC_Wilcoxon_CSS = dt$Nr_Wilcoxon_Rejections_Normalized_CSS,
DISC_Wilcoxon_TSS = dt$Nr_Wilcoxon_Rejections_Normalized,
Discoveries = dt$rejections,
Shared = dt$shared,
ReferenceSize = dt$ref_size_selected)


library(xtable)
print(xtable(dt_export), include.rownames=FALSE)

