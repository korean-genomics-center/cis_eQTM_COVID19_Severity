# %%
import os
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu, kruskal, ranksums

# %%
list_exclude_checking = ["Number", "KU10K_ID", "Birth-date", "PCI-date", "last-FU_date", "Death_date", "AMI_date", "Any_repeat_revascularization-date", "종류", "기타", "AMI_여부_from_의사"]
list_additional_categorical = ["HTN", "DM", "Current-smok", "FHX", "prev-MI", "Prev-PCI", "Prev-CABG", "ECG", "Disch-aspirin", "Disch-P2Y12inh", "Disch-statin", "Disch-Bblock", "Death", "AMI", "Any_repeat_revascularization", "TypeI_MI", "Stent_PCI_여부"]

#%%
file_crf = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/9_clinical/Infectomics_Severity_Information_Methyl_20240102.tsv"
file_dmp = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/DMPDEG/first/cpgbeta_overlap.tsv"
table_crf = pd.read_csv(file_crf, sep="\t")
table_cpg = pd.read_csv(file_dmp, sep="\t")

# %%
def calcuate_ratio(a, b):
    if pd.isna(a) or pd.isna(b):
        return pd.NA
    else:
        return a/b

table_crf["Neutrophil_to_Lymphocyte_Ratio"] = table_crf.apply(lambda x : calcuate_ratio(x["neutrophil"], x["lymphocyte"]), axis = 1)

dict_column_type = dict()
for col in table_crf.columns:
    if col in list_exclude_checking:
        continue
    if table_crf[col].dtypes == object:
        dict_column_type[col] = "Cat"
    elif table_crf[col].dtypes == float:
        if col in list_additional_categorical:
            dict_column_type[col] = "Cat"
        else:
            dict_column_type[col] = "Seq"
    else:
        raise Exception(col)
dict_column_type["Neutrophil_to_Lymphocyte_Ratio"] = "Seq"

# %%
list_checking_only = ["Age","Sex","BMI","HTN","DM","Current-smok","T-chol","LDL-chol","HDL-chol","TG","WBC","RBC","platelets","neutrophil","lymphocyte","monocyte","eosinophils","Neutrophil_to_Lymphocyte_Ratio","Basophil","CPK_peak","CKMB_peak","TnI_peak","EF","Base_TIMI","Number_diseased_vessel"]
table_clinical_significance = pd.DataFrame(columns = ["chr", "start", "end"] + list_checking_only)
list_sample = table_crf["KU10K_ID"].to_list()
for ind, row in table_cpg.iterrows():
    dict_row = row.to_dict()
    list_sample_check = list(set(list_sample) & set(dict_row.keys()))
    table_clinical_significance.loc[ind, "chr"] = dict_row["chr"]
    table_clinical_significance.loc[ind, "start"] = dict_row["start"]
    table_clinical_significance.loc[ind, "end"] = dict_row["end"]
    
    list_methylation = list(map(dict_row.__getitem__, list_sample_check))
    for col in list_checking_only:
        dtype = dict_column_type[col]
        dict_sample_to_clinical = dict(zip(table_crf["KU10K_ID"], table_crf[col]))
        list_clinical = list(map(dict_sample_to_clinical.__getitem__, list_sample_check))
        
        zip_methyl_clinical = list(zip(list_methylation, list_clinical))
        zip_methyl_clinical_nonna = list(filter(lambda x : pd.notna(x[1]), zip_methyl_clinical))
        
        list_methyl_check = list(map(lambda x : x[0], zip_methyl_clinical_nonna))
        list_clinical_check = list(map(lambda x : x[1], zip_methyl_clinical_nonna))
        
        if dtype == "Cat":
            list_cat_type = list(set(list_clinical_check))
            list_matrix = list()
            for cat_type in list_cat_type:
                zip_methyl_clinical_cat = list(filter(lambda x : x[1] == cat_type, zip_methyl_clinical_nonna))
                list_methyl_this_cat = list(map(lambda x : x[0], zip_methyl_clinical_cat))
                list_matrix.append(list_methyl_this_cat)
            if len(list_cat_type) == 2:
                _, pval = ranksums(list_matrix[0], list_matrix[1])
                stat = abs(np.mean(list_matrix[0]) - np.mean(list_matrix[1]))
            elif len(list_cat_type) == 1:
                pval = 1
            else:
                _, pval = kruskal(*list_matrix)
        else:
            stat, pval = spearmanr(list_methyl_check, list_clinical_check)
        table_clinical_significance.loc[ind, col] = pval
        table_clinical_significance.loc[ind, f"{col}_Stat"] = stat