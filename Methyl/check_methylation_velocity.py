#%%
import pandas as pd

path_crf = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Resources/Data/infectomics_CRF_20230410_edit.xlsx"
path_meta = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Resources/MetaTable/metatable_combined_all_231114.tsv.cpg_table_file_path.tsv"
path_cpg = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/MethylCpGTable/Infectomics.Copy_From_HjRyu/MethylCpGTable.Control.Mild.Case.Severe.tsv"

#%%
table_crf = pd.read_excel(path_crf, engine = "openpyxl", skiprows = 1)
# %%
from datetime import datetime

table_crf["Visit date_datetime"] = table_crf["Visit date"].apply(lambda x : x.date())
# %%
dict_sample_to_date = dict(zip(table_crf["Subject NO.(고유번호)"], table_crf["Visit date_datetime"]))
# %%
table_cpg = pd.read_csv(path_cpg, sep = '\t')
table_meta = pd.read_csv(path_meta, sep = '\t')
# %%
dict_subject_to_samples = dict()
for subjid in table_meta["Subject_ID"].dropna().unique():
    list_samples_subj = table_meta[table_meta["Subject_ID"] == subjid]["Sample_ID"].to_list()
    assert len(list_samples_subj) == 2, list_samples_subj
    dict_subject_to_samples[subjid] = list_samples_subj
# %%
dict_subject_to_date_interval = dict()
for subj, list_sample in dict_subject_to_samples.items():
    sample1, sample2 = list_sample
    date1, date2 = dict_sample_to_date[sample1], dict_sample_to_date[sample2]
    dict_subject_to_date_interval[subj] = (date2 - date1).days
#%%
dict_subject_to_diffs = dict()
for subj, list_sample in dict_subject_to_samples.items():
    sample1, sample2 = list_sample
    diff = table_cpg[sample2] - table_cpg[sample1]
    dict_subject_to_diffs[subj] = diff.to_numpy()
# %%
dict_sample_to_severity = dict(zip(table_crf["Subject NO.(고유번호)"], table_crf["중증도분류"]))
dict_subject_to_severity = dict()
for subj, list_sample in dict_subject_to_samples.items():
    sample1, sample2 = list_sample
    sev1, sev2 = dict_sample_to_severity[sample1], dict_sample_to_severity[sample2]
    assert sev1 == sev2, subj
    dict_subject_to_severity[subj] = sev1
# %%
from matplotlib import pyplot as plt
import numpy as np
##################################
fig, ax = plt.subplots(2, 1, figsize = (12, 7))
axes = ax.flatten()

dict_subject_to_meandiff = dict()
for subj, list_diff in dict_subject_to_diffs.items():
    dict_subject_to_meandiff[subj] = np.mean(list_diff)
sorted_subj_by_meandiff = sorted(dict_subject_to_meandiff.keys(), key = dict_subject_to_meandiff.__getitem__, reverse = True)
axes[0].boxplot(list(map(dict_subject_to_diffs.__getitem__, sorted_subj_by_meandiff)))
axes[0].set_xticks([])
axes[0].set_ylabel("Methylation_Difference")
axes[1].scatter(sorted_subj_by_meandiff, list(map(dict_subject_to_severity.__getitem__, sorted_subj_by_meandiff)))
axes[1].set_xticks([])
axes[1].set_xlim([-0.5, len(dict_subject_to_severity)-0.5])
axes[1].set_ylabel("Severity")
plt.show()
# %%
##################################
list_subj = list(dict_subject_to_meandiff.keys())
list_mean_diff = list(map(dict_subject_to_meandiff.__getitem__, list_subj))
list_severity = list(map(dict_subject_to_severity.__getitem__, list_subj))

import seaborn as sns

sns.regplot(x = list_severity, y = list_mean_diff)
plt.xlim(0.5, 4.5)
plt.xticks([1,2,3,4])
plt.xlabel("Severity")
plt.ylabel("Mean Methylation Difference\nAcross Whole Genome (%)")
plt.show()
# %%
##################################
list_subj = list(dict_subject_to_meandiff.keys())
list_mean_diff_rate = list(map(lambda subj : dict_subject_to_meandiff[subj] / dict_subject_to_date_interval[subj], list_subj))
list_severity = list(map(dict_subject_to_severity.__getitem__, list_subj))

import seaborn as sns

sns.regplot(x = list_severity, y = list_mean_diff_rate)
plt.xlim(0.5, 4.5)
plt.xticks([1,2,3,4])
plt.xlabel("Severity")
plt.ylabel("Mean Methylation Difference per Day\nAcross Whole Genome (%/day)")
plt.show()
# %%
##################################
dict_subject_to_absmeandiff = dict()
for subj, list_diff in dict_subject_to_diffs.items():
    dict_subject_to_absmeandiff[subj] = np.mean(abs(list_diff))

list_subj = list(dict_subject_to_absmeandiff.keys())
list_absmean_diff_rate = list(map(lambda subj : dict_subject_to_absmeandiff[subj] / dict_subject_to_date_interval[subj], list_subj))
list_severity = list(map(dict_subject_to_severity.__getitem__, list_subj))

import seaborn as sns

sns.regplot(x = list_severity, y = list_absmean_diff_rate)
plt.xlim(0.5, 4.5)
plt.xticks([1,2,3,4])
plt.xlabel("Severity")
plt.ylabel("Mean |Methylation Difference| per Day\nAcross Whole Genome (%/day)")
plt.show()
# %%
from scipy.stats import spearmanr

spearmanr(list_severity, list_absmean_diff_rate)
# %%
###############################
list_cols_crf_check = ["WBC","RBC","Hb","Hct","PLT","Neutrophil","Lymphocytes","Monocytes","Eosinophils","Basophils","Calcium ","Phosphorus","Glucose","BUN","Creatinine","Sodium","Potassium","Chloride","Uric acid","Total cholesterol","LDL-cholesterol","HDL-cholesterol","Triglyceride","Total protein","Albumin","AST(SGOT)","ALT(SGPT)","Alkaline Phosphatase","r-GTP","Total bilirubin ","L D H","CPK","CK-MB","Troponin I","D-dimer","FDP","Anti-thrombin III","PT","INR","Activated partial\nthromboplastin time","Protein C Activity","Protein S Activity","Lupus-like anticoagulant Screening","Anti beta 2-GP 1 lgG","Anti beta 2-GP 1 lgG","Anti beta 2-GP 1 lgM","Anti beta 2-GP 1 lgM","Anti Cardiolipin Ab lgG","Anti Cardiolipin Ab lgG","Anti Cardiolipin Ab lgM","Anti Cardiolipin Ab lgM","Anti-Phospholipid lgM","Anti-Phospholipid lgM","Anti-Phospholipid lgG","Anti-Phospholipid lgG","Ferritin","CRP","ESR","IL-6","HbA1C","TSH","T3","Free T4"]
# %%
dict_subject_to_absmeandiff_per_day = dict()
for subj, meanabsdiff in dict_subject_to_absmeandiff.items():
    dict_subject_to_absmeandiff_per_day[subj] = meanabsdiff/dict_subject_to_date_interval[subj]
#%%
list_p = list()

list_sample_visit1 = list(map(lambda x : dict_subject_to_samples[x][0], list_subj))
for col_check in list_cols_crf_check:
    dict_sample_to_val_check = dict(zip(table_crf["Subject NO.(고유번호)"], table_crf[col_check]))
    sample_visit1_to_val_check_pair = list(zip(list_sample_visit1, list(map(dict_sample_to_val_check.__getitem__, list_sample_visit1))))
    sample_visit1_to_val_check_pair_dropna = list(filter(lambda pair : pd.notna(pair[1]), sample_visit1_to_val_check_pair))
    
    subjs_nonna = list(map(lambda x : x[0], sample_visit1_to_val_check_pair_dropna))
    values_nonna = list(map(lambda x : x[1], sample_visit1_to_val_check_pair_dropna))
    absmeandiff_per_day_nonna = list(map(lambda sample : dict_subject_to_absmeandiff_per_day['-'.join(sample.split('-')[:-1])], subjs_nonna))
    
    corr, p_spearman = spearmanr(values_nonna, absmeandiff_per_day_nonna)
    list_p.append(p_spearman)
    if p_spearman * len(list_cols_crf_check) < 0.05:
        sns.regplot(x = values_nonna, y = absmeandiff_per_day_nonna)
        plt.xlabel(col_check)
        plt.ylabel("Mean |Methylation Difference| per Day\nAcross Whole Genome (%/day)")
        plt.show()
        print(corr, p_spearman)
# %%
from statsmodels.stats.multitest import fdrcorrection

col_to_pval_pair = list(zip(list_cols_crf_check, list_p))

col_to_pval_pair_nonna = list(filter(lambda pair : pd.notna(pair[1]), col_to_pval_pair))

_, fdr = fdrcorrection(list(map(lambda x : x[1], col_to_pval_pair_nonna)))

col_to_fdr_pair_nonna = list(zip(list(map(lambda x : x[0], col_to_pval_pair_nonna)), fdr))
# %%
