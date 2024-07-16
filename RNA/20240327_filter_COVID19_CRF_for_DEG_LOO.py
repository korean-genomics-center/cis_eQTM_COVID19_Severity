# %%
import os
import pandas as pd
import glob


# %%
path_crf = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Resources/Data/COVID19_master_table_added_CRF_20231007.txt"
path_exp = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/3_rsem/Confirmed_Matched_with_Methyl_FU_loss"
dir_out = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Resources/Data"
# %%
list_file_exp = os.listdir(path_exp)
df_crf = pd.read_csv(path_crf, sep="\t")
df_crf_exp = df_crf[df_crf["Project_ID_Alias"].isin(list_file_exp)]
list_severe_subjects = df_crf_exp[df_crf_exp["Severity_Binary"]=="Severe"]["Project_ID"].unique()
for sev_subj in list_severe_subjects:
    filename = os.path.join(dir_out, f"COVID19_master_table_removed_CRF_FU_loss_removed_loo_{sev_subj}.txt")
    df_crf_exp_loo = df_crf_exp.copy()
    df_crf_exp_loo_removed_one_sample = df_crf_exp_loo[df_crf_exp_loo["Project_ID"] != sev_subj]
    df_crf_exp_loo_removed_one_sample.to_csv(filename, sep="\t", index=False)
    list_severe_subjects_loo = df_crf_exp_loo_removed_one_sample[df_crf_exp_loo_removed_one_sample["Severity_Binary"]=="Severe"]["Project_ID"].unique()
    print(list_severe_subjects_loo)
    print(len(list_severe_subjects_loo))
# %%

# %%
