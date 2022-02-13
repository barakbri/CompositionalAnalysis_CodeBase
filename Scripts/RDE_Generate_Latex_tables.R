#Script is used for generating latex format tables, for the result of data analysis

library(xtable)

# tables for the gut CD data set
dt = read.csv('../../Results/gut_cell_shared_disc_mat_for_paper.csv')
row.names(dt) = NULL
xtable::xtable(dt[-7,-8])

# table for the HMP study:
dt = read.csv('../../Results/HMP/results.csv')
names(dt)


dt_export = data.frame(Site1 = dt$Site1,
Site2 = dt$Site2,
NR.Taxa = dt$taxa_for_test,
ANCOM = paste0(dt$ANCOM_rejections,' (',dt$shared,')'),
DISC_Wilcoxon_CSS = paste0(dt$Nr_Wilcoxon_Rejections_Normalized_CSS),
DISC_ALDEx2t = paste0(dt$ALDEx2_Welch,' (',dt$shared_ALDEx2_Welch,')'),
Discoveries = dt$rejections,
ReferenceSize = dt$ref_size_selected)

library(xtable)
print(xtable(dt_export), include.rownames=FALSE)

dt_export_sensitivity  = data.frame(Site1 = dt$Site1,
                       Site2 = dt$Site2,
                       contamination = dt$nr_contamination_ref,
                       nr_removed_from_ref = dt$nr_removed_from_ref,
                       contamination_in_reduced = dt$nr_contamination_ref_cleaned,
                       nr_discoveries_in_original = dt$rejections,
                       nr_discoveries_in_reduced = paste0(dt$nr_rejections_cleaned,'(',dt$nr_rejections_cleaned_shared,')')
                       )
dt_export_sensitivity$contamination = ifelse(dt_export_sensitivity$contamination>0,'Y',no = 'N')
dt_export_sensitivity$contamination_in_reduced = ifelse(dt_export_sensitivity$contamination_in_reduced>0,'Y',no = 'N')
ind_to_dash = which((dt_export_sensitivity$contamination == 'Y' &
                      dt_export_sensitivity$nr_removed_from_ref==0) | dt_export_sensitivity$contamination == 'N')
dt_export_sensitivity$contamination_in_reduced[ind_to_dash] = '-'
dt_export_sensitivity$nr_discoveries_in_reduced[ind_to_dash] = '-'
print(xtable(dt_export_sensitivity), include.rownames=FALSE)

