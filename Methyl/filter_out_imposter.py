#%%
import pandas as pd

#%%
path_meta = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Resources/MetaTable/metatable_combined_all_231114.tsv.cpg_table_file_path.relink_CpG_Genecode42.remove_00446.tsv"
path_save = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Resources/MetaTable/metatable_combined_all_231114.tsv.cpg_table_file_path.relink_CpG_Genecode42.remove_highNA_Control.tsv"

dict_imposter = {
    "U10K-00446" : "Low Breadth of Coverage",
    "KU10K-05913" : "High NA Value",
    "KU10K-05954" : "High NA Value",
    "KU10K-05575" : "High NA Value",
    "KU10K-04492" : "High NA Value",
    "KU10K-05044" : "High NA Value"
}

table = pd.read_csv(path_meta, sep = '\t')
table_remove_imposter = table[list(map(lambda x : dict_imposter.get(x) == None, table["Sample_ID"]))]
table_remove_imposter.to_csv(path_save, sep = '\t', index = False)
# %%
