#%%
import pandas as pd

path_table = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Resources/MetaTable/metatable_combined_all_231114.tsv.cpg_table_file_path.remove_00446.tsv"
path_save = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Resources/MetaTable/metatable_combined_all_231114.tsv.cpg_table_file_path.relink_CpG_Genecode42.remove_00446.tsv"

path_cpg_formats = [
    "/BiO/Research/Project2/Infectomics_COVID-19_Host/Results/Infectomics_COVID-19_Methyl/Include_Second_Year_Collection/MethylCpGMin/HealthyControl/{0}.pair_merged.methyl_cpg_min.tsv",
    "/BiO/Research/Project2/Infectomics_COVID-19_Host/Results/Infectomics_COVID-19_Methyl/Include_Second_Year_Collection/MethylCpGMin/{0}/{0}.pair_merged.methyl_cpg_min.tsv"
]

#%%
table_meta = pd.read_csv(path_table, sep = '\t')
# %%
import os

def find_available_file_path(list_path_format, val):
    list_available = list()
    for path_format in list_path_format:
        path = path_format.format(val)
        if os.path.exists(path):
            list_available.append(path)
    if len(list_available) != 1:
        raise Exception(val)
    return list_available[0]

table_meta["Path_to_CpG_Table"] = table_meta["Sample_ID"].apply(lambda x : find_available_file_path(path_cpg_formats, x))
# %%
table_meta.to_csv(path_save, sep = '\t', index = False)
# %%
