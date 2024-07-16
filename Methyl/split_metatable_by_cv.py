#%%
import os, math, random
from copy import deepcopy
import pandas as pd
import numpy as np

#%%
# Config ############################################################
path_meta = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Resources/MetaTable/metatable_combined_all_231114.tsv.cpg_table_file_path.relink_CpG_Genecode42.remove_highNA_Control.tsv"

path_save_fold_splitted_total = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Resources/MetaTable/metatable_combined_all_231114.tsv.cpg_table_file_path.relink_CpG_Genecode42.remove_highNA_Control.Case_Severe.Control_Mild.CV_Splitted.filter.Methyl_Visit_Order_First.tsv"

dir_save = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Resources/MetaTable/Cross_Validation_Split/Case_Severe.Control_Mild.CV_Case_9.CV_Control_1"
filename_format = "metatable_combined_all_231114.remove_highNA_Control.Case_{case_name}.Control_{control_name}.CV_Case_{num_cv_case}.Fold_{case_fold}.CV_Control_{num_cv_control}.Fold_{control_fold}.filter.{extrafilter}.tsv"

col_sampleid = "Sample_ID"

col_group = "Severity_group"
groupnames_case = ["Severe"]
groupnames_control = ["Mild"]

dict_filter = {
    "Methyl_Visit_Order" : ["First"]
}

col_assign_fold = "CV_Fold"
random_state = 0
num_fold_split_control = 1
num_fold_split_case = 9
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
    
def save_table_by_cv_fold_combination(table_control, table_case, names_control, names_case, num_fold_control, num_fold_case, dict_filter, col_fold, dir_save, filename_format):
    os.makedirs(dir_save, exist_ok=True)
    list_name_filtered = list()
    for filt_key, filt_values in dict_filter.items():
        list_name_filtered.append(f"{filt_key}_{'_'.join(filt_values)}")
    name_filtered = '.'.join(list_name_filtered)
    for ind_fold_exclude_control in range(num_fold_control):
        for ind_fold_exclude_case in range(num_fold_case):
            if num_fold_control > 1:
                table_control_exc_fold = table_control[table_control[col_fold] != ind_fold_exclude_control]
            else:
                table_control_exc_fold = table_control
            if num_fold_case > 1:
                table_case_exc_fold = table_case[table_case[col_fold] != ind_fold_exclude_case]
            else:
                table_case_exc_fold = table_case
            table_all_fold_combination = pd.concat([table_control_exc_fold, table_case_exc_fold])
            print(table_all_fold_combination.shape)
            filename_save = filename_format.format(case_name = '_'.join(names_case), control_name = '_'.join(names_control), num_cv_case = num_fold_case, num_cv_control = num_fold_control, case_fold = ind_fold_exclude_case, control_fold = ind_fold_exclude_control, extrafilter = name_filtered)
            path_save = os.path.join(dir_save, filename_save)
            table_all_fold_combination.to_csv(path_save, sep = '\t', index = False)


#%%
if __name__ == "__main__":
    table_meta = pd.read_csv(path_meta, sep = '\t')

    table_meta_filtered = filter_table_by_dict_filters(table_meta, dict_filter)
    table_meta_filtered_group = assign_control_case_to_table(table_meta_filtered, col_group, groupnames_control, groupnames_case)
    table_meta_filtered_group_control = table_meta_filtered_group[table_meta_filtered_group[f"{col_group}_DMP_Group"] == "Control"]
    table_meta_filtered_group_case = table_meta_filtered_group[table_meta_filtered_group[f"{col_group}_DMP_Group"] == "Case"]
    
    table_meta_filtered_group_control_fold = random_assign_fold_number(table_meta_filtered_group_control, num_fold_split_control, random_state, "Fold")
    table_meta_filtered_group_case_fold = random_assign_fold_number(table_meta_filtered_group_case, num_fold_split_case, random_state, "Fold")
    save_table_by_cv_fold_combination(table_meta_filtered_group_control_fold, table_meta_filtered_group_case_fold, groupnames_control, groupnames_case, num_fold_split_control, num_fold_split_case, dict_filter, "Fold", dir_save, filename_format)

    pd.concat([table_meta_filtered_group_control_fold, table_meta_filtered_group_case_fold]).to_csv(path_save_fold_splitted_total, sep = '\t', index = False)
# %%
