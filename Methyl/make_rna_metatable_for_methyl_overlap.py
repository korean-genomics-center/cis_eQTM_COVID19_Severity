#%%
import pandas as pd
import os


path_rna_meta = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Resources/Data/COVID19_master_table_20231007.txt"
path_methyl_cpg = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231211/combined_all_samples_meta_Rds/MethylKitTable/MethylCpGTable.Control.Mild.Case.Severe.tsv"

col_methyl_nonsample = ["chr", "start", "end"]

with open(path_methyl_cpg, 'r') as fr:
    header = fr.readline().strip('\n').split('\t')
    list_methyl_samples = list(filter(lambda x : x not in col_methyl_nonsample, header))
    
# %%
table_rna_meta = pd.read_csv(path_rna_meta, sep = '\t')
# %%
table_rna_meta_with_methyl = table_rna_meta[table_rna_meta["Project_ID_Alias"].isin(list_methyl_samples)]
# %%
table_rna_meta_with_methyl = table_rna_meta_with_methyl.drop_duplicates(subset = ["Project_ID_Alias"])

table_rna_meta_with_methyl = table_rna_meta_with_methyl[table_rna_meta_with_methyl["Sample_Trait"] == "ViralInfection"]

samples_exclude = ["C19-C014-V1", "C19-C014-V4"]

table_rna_meta_with_methyl = table_rna_meta_with_methyl[list(map(lambda x : x not in samples_exclude, table_rna_meta_with_methyl["Project_ID_Alias"]))]

# %%
table_rna_meta_with_methyl.to_csv("/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Resources/MetaTable/COVID19_master_table_20231007.Methyl_Overlap.20240402.txt", sep = '\t', index = False)
# %%
table_rna_meta_with_methyl[table_rna_meta_with_methyl["Sample_Sex"] == "M"].to_csv("/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Resources/MetaTable/COVID19_master_table_20231007.Methyl_Overlap.20240402.Sex_M.txt", sep = '\t', index = False)
table_rna_meta_with_methyl[table_rna_meta_with_methyl["Sample_Sex"] == "F"].to_csv("/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Resources/MetaTable/COVID19_master_table_20231007.Methyl_Overlap.20240402.Sex_F.txt", sep = '\t', index = False)
# %%
from collections import Counter
Counter(table_rna_meta_with_methyl["Visit_order"])
# %%
table_meta_with_severity = pd.read_csv("/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/9_clinical/Infectomics_Severity_Information_Methyl_20240102.tsv", sep = '\t')
# %%
dict_sample_to_sev = dict(zip(table_meta_with_severity["Sample_ID"], table_meta_with_severity["Severity"]))
table_rna_meta_with_methyl["Severity"] = table_rna_meta_with_methyl["Project_ID_Alias"].apply(dict_sample_to_sev.__getitem__)

#%%
table_rna_meta_with_methyl.to_csv("/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Resources/MetaTable/COVID19_master_table_20231007.Methyl_Overlap.with_Severity.20240402.txt", sep = '\t', index = False)
# %%
