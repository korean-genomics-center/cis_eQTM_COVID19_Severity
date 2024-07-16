#%%
import os, math, random
from copy import deepcopy
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold

#%%
# Config ############################################################
path_meta = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Resources/MetaTable/metatable_combined_all_231114.tsv.cpg_table_file_path.relink_CpG_Genecode42.remove_highNA_Control.tsv"
path_meta_sev = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/9_clinical/Infectomics_Severity_Information_Methyl_20240102.tsv"

path_save_fold_splitted_total = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Resources/MetaTable/metatable_combined_all_231114.tsv.cpg_table_file_path.relink_CpG_Genecode42.remove_highNA_Control.Case_Sev234.Control_Sev1.CV_Splitted.filter.Methyl_Visit_Order_First.tsv"

dir_save = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Resources/MetaTable/Cross_Validation_Split/Case_Sev234.Control_Sev1.CV_4Fold"
filename_format = "metatable_combined_all_231114.remove_highNA_Control.Case_Sev{case_name}.Control_Sev{control_name}.CV_4.Fold_{num_fold}.filter.{extrafilter}.tsv"

col_sampleid = "Sample_ID"

col_group = "Severity"
groupnames_case = [2,3,4]
groupnames_control = [1]

dict_filter = {
    "Methyl_Visit_Order" : ["First"]
}

col_assign_fold = "Fold"
random_state = 0
num_fold = 4
cv_stratify = "Severity"
#######################################################################

#%%
def filter_table_by_dict_filters(table, dict_filter):
    table_filt = table.copy()
    for col, values in dict_filter.items():
        table_filt = table_filt[table_filt[col].isin(values)]
    return table_filt

def assign_control_case_to_table(table, col_group, values_control, values_case):
    dict_value_to_group = dict()
    for val_control in values_control:
        dict_value_to_group[val_control] = "Control"
    for val_case in values_case:
        dict_value_to_group[val_case] = "Case"
    table[f"{col_group}_DMP_Group"] = table[col_group].apply(dict_value_to_group.__getitem__)
    table = table[table[f"{col_group}_DMP_Group"] != None]
    return table

def random_assign_fold_number(table, num_fold_split, random_state, col_add_fold_info):
    num_data = table.shape[0]
    list_fold_pool = deepcopy(list(range(num_fold_split)) * math.ceil(num_data / num_fold_split))
    print(list_fold_pool)
    list_fold_pool = list_fold_pool[:num_data]
    
    np.random.seed(random_state)
    list_shuffled_fold = np.random.permutation(list_fold_pool)
    print(list_shuffled_fold)
    table[col_add_fold_info] = list_shuffled_fold
    return table

#%%
table_meta = pd.read_csv(path_meta, sep = '\t')
table_meta_sev = pd.read_csv(path_meta_sev, sep = '\t')

#%%
table_meta_filtered = filter_table_by_dict_filters(table_meta, dict_filter)
dict_sample_to_sev = dict(zip(table_meta_sev["Sample_ID"], table_meta_sev["Severity"]))
table_meta_filtered["Severity"] = table_meta_filtered["Sample_ID"].apply(dict_sample_to_sev.__getitem__)
#%%
table_meta_filtered_group = assign_control_case_to_table(table_meta_filtered, col_group, groupnames_control, groupnames_case)
#%%
table_meta_filtered_group_control = table_meta_filtered_group[table_meta_filtered_group[f"{col_group}_DMP_Group"] == "Control"]
table_meta_filtered_group_case = table_meta_filtered_group[table_meta_filtered_group[f"{col_group}_DMP_Group"] == "Case"]
#%%
def split_data_kfold_stratified(table, col_sampleid, col_stratify, num_fold, random_state, colname_fold):
    table = table.reset_index(drop = True)
    kfold_splitter = StratifiedKFold(n_splits=num_fold, shuffle=True, random_state=random_state)

    for i, (train_index, test_index) in enumerate(kfold_splitter.split(table[col_sampleid].to_list(), table[col_stratify].to_list())):
        table.loc[test_index, colname_fold] = i
    return table

table_meta_filtered_group_control_fold = split_data_kfold_stratified(table_meta_filtered_group_control, "Sample_ID", col_group, num_fold, random_state, "Fold")
table_meta_filtered_group_case_fold = split_data_kfold_stratified(table_meta_filtered_group_case, "Sample_ID", col_group, num_fold, random_state, "Fold")

#%%
table_meta_filtered_group_control_fold["Severity_Weight"] = 1
table_meta_filtered_group_case_fold["Severity_Weight"] = table_meta_filtered_group_case_fold["Severity"] - 1
# %%
table_tot_fold = pd.concat([table_meta_filtered_group_control_fold, table_meta_filtered_group_case_fold]).reset_index(drop = True)

table_tot_fold.to_csv(path_save_fold_splitted_total, sep = '\t', index = False)

os.makedirs(dir_save, exist_ok=True)
list_name_filtered = list()
for filt_key, filt_values in dict_filter.items():
    list_name_filtered.append(f"{filt_key}_{'_'.join(filt_values)}")
name_filtered = '.'.join(list_name_filtered)
case_name = ''.join(list(map(str, groupnames_case)))
control_name = ''.join(list(map(str, groupnames_control)))

for fold_ind in range(num_fold):
    filename_save = filename_format.format(case_name = case_name, control_name = control_name, num_fold = fold_ind, extrafilter = name_filtered)
    path_save = os.path.join(dir_save, filename_save)
    table_fold = table_tot_fold[table_tot_fold["Fold"] != fold_ind]
    
    table_fold.to_csv(path_save, sep = '\t', index = False)
# %%
