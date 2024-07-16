#%%
import pandas as pd
import os


col_diff = "wmeth.diff"
list_cutoff = [1, 2, 3, 4, 5, 10]
num_fold = 4

path_dmp_format = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/DMPExtract/Infectomics.Case_Sev234.Control_Sev1.CV_4.filter.Methyl_Visit_Order_First.Weight_Severity/DMPExtract.Infectomics.Case_Sev234.Control_Sev1.CV_4.Fold_{{fold}}.Weight_Severity.20240314.tsv.FDR_05.{col_diff}_abs_{{cutoff}}.tsv"

dict_marker = dict()

def check_direction_consistency(list_values):
    is_neg = list(map(lambda x : x < 0, list_values))
    if sum(is_neg) == 0 or sum(is_neg) == len(list_values):
        return True
    else:
        return False

for cutoff in list_cutoff:
    dict_marker[cutoff] = dict()
    dict_posname_to_methdiff = dict()
    for fold in range(num_fold):
        path_dmp = path_dmp_format.format(fold = fold, cutoff = cutoff)
        table_dmp = pd.read_csv(path_dmp, sep = '\t')
        table_dmp["PosName"] = table_dmp.apply(lambda x : f"{x['chr']}:{x['start']}", axis = 1)
        for posname, methdiff in zip(table_dmp["PosName"], table_dmp["meth.diff"]):
            if dict_posname_to_methdiff.get(posname) == None:
                dict_posname_to_methdiff[posname] = list()
            dict_posname_to_methdiff[posname].append(methdiff)
    list_posname_consistent = list(filter(lambda posname: check_direction_consistency(dict_posname_to_methdiff[posname]), dict_posname_to_methdiff.keys()))
    list_posname_inconsistent = list(set(dict_posname_to_methdiff.keys()) - set(list_posname_consistent))
    for pos_inconsistent in list_posname_inconsistent:
        dict_posname_to_methdiff.pop(pos_inconsistent)
    dict_marker[cutoff] = dict_posname_to_methdiff
# %%
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from collections import Counter

plt.rcParams["font.size"] = 16

for cutoff in list_cutoff:
    list_num_sig_per_cv = list(map(len, dict_marker[cutoff].values()))
    count_consistent_sig = Counter(list_num_sig_per_cv)
    plt.bar(range(1, num_fold+1), list(map(count_consistent_sig.__getitem__, range(1, num_fold+1))), color = "gray")
    plt.xlabel("Number of Overlaps")
    plt.ylabel("Number of Biomarkers")
    plt.title(f"Methdiff Cutoff : {cutoff}%")
    plt.yscale("log")
    plt.show()
    # break
# %%
