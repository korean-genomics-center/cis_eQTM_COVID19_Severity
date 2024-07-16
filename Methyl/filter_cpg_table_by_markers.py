#%%
import pandas as pd

path_cpg_tot = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/MethylCpGTable/Infectomics.Copy_From_HjRyu/MethylCpGTable.Control.Mild.Case.Severe.Filtered.DMP.Hyper_Hypo.Sev_vs_Mild.Visit1.LOO_common.tsv"
path_save_cpg_filtered = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/MethylCpGTable/Infectomics.Copy_From_HjRyu/MethylCpGTable.Control.Mild.Case.Severe.Filtered.DMP.Hyper_Hypo.Sev_vs_Mild.Visit1.LOO_common.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_6.tsv"
path_dmp = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/Methyl_RNA_Relation/Methyl_RNA_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Mild_Sev_Visit1.LOO_common.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_6.DMP_DEG_LOO_annot.20240402.tsv"

cols_check = ["chr", "start", "end"]

table_cpg_tot = pd.read_csv(path_cpg_tot, sep = '\t')
table_dmp = pd.read_csv(path_dmp, sep = '\t')

# featname_include = table_dmp[cols_check].apply(lambda row : '_'.join(list(map(str, row))), axis = 1).to_list()
featname_include = table_dmp["Methyl"].to_list()

table_cpg_tot.index = table_cpg_tot[cols_check].apply(lambda row : '_'.join(list(map(str, row))), axis = 1)

table_cpg_tot = table_cpg_tot[table_cpg_tot.index.isin(featname_include)]

table_cpg_tot.to_csv(path_save_cpg_filtered, sep = '\t', index = False)
# %%
