# %%
import os
import glob
import numpy as np
import pandas as pd
import pickle
import gzip
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from sklearn.metrics import roc_curve

from sklearn.metrics import roc_auc_score
# %%
def main(path_sev_info, dir_methylcpgmin, list_drop_samples, infilename, outfilename):
    list_files_methylcpgmin = get_list_files_methylcpgmin(dir_methylcpgmin)
    list_methyl_sample = get_list_methyl_sample(list_files_methylcpgmin)
    list_methyl_sample = list(filter(lambda x: "C19-C" in x, list_methyl_sample))
    list_methyl_sample = list(filter(lambda x: x not in list_drop_samples , list_methyl_sample))
    list_methyl_sample = list(filter(lambda x: "L1" not in x, list_methyl_sample))
    # dict_severity_group = get_dict_severity_group(path_sev_info, list_methyl_sample)
    # dict_severity_group = remove_severity_group_no_content(dict_severity_group)
    dict_cpgmin_all_samples = load_pickle(infilename)     
    dict_cpgmin_all_samples_marker_fixed = dict()
    for sampleid, dict_cpgs in dict_cpgmin_all_samples.items():
        if sampleid in list_methyl_sample:
            dict_cpgs_marker_fixed = dict()
            for marker, freqC in dict_cpgs.items():
                fixed_marker = ":".join(marker)
                dict_cpgs_marker_fixed[fixed_marker] = freqC
            
            dict_cpgmin_all_samples_marker_fixed[sampleid] = dict_cpgs_marker_fixed

    df_cpgmin_all_sample = pd.DataFrame(dict_cpgmin_all_samples_marker_fixed).astype(float)
    df_all_sample_cpgmin = df_cpgmin_all_sample.T.reset_index(drop=False).rename(columns={"index": "Sample_ID"})

    # # imputation by median if needed (deprecated)
    # df_all_sample_cpgmin_imputed_na = df_all_sample_cpgmin.copy()
    # list_markers = list(df_all_sample_cpgmin_imputed_na.columns[1:])
    # for marker in list_markers:
    #     median = df_all_sample_cpgmin_imputed_na[marker].median()
    #     df_all_sample_cpgmin_imputed_na[marker] = df_all_sample_cpgmin_imputed_na[marker].fillna(median)

    df_sev = pd.read_csv(path_sev_info, sep="\t")
    df_merged = pd.merge(df_all_sample_cpgmin, df_sev, on="Sample_ID", how="inner")
    df_merged = df_merged.set_index("Sample_ID")

    # list_columns = list(df_merged.columns)
    # list_markers = list(filter(lambda x: x.startswith("chr"), list_columns))
    # mean_beta_sev_per_marker_first = df_merged.groupby("Severity_visit")[list_markers[0]].apply(np.mean)
    # index_col = mean_beta_sev_per_marker_first.index
    # df_mean_beta_sev_per_marker = pd.DataFrame(index=index_col)

    # for marker in list_markers:
    #     mean_beta_sev_per_marker = df_merged.groupby("Severity_visit")[marker].apply(np.mean)
    #     array_mean_beta = (mean_beta_sev_per_marker.to_numpy())
    #     df_mean_beta_sev_per_marker[marker] = array_mean_beta

    direction = "hyper" 
    annot_methyl = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231020/annotation/methyl_{direction}_severe_mild_annotatr.tsv"
    df_annot_methyl = pd.read_csv(annot_methyl, sep="\t")
    df_annot_methyl["marker"] = df_annot_methyl["seqnames"].astype(str) + ":" + df_annot_methyl["start"].apply(lambda x: int(x)-1).astype(str)

    dict_marker_symbol = dict()
    for marker, symbol in zip(df_annot_methyl["marker"], df_annot_methyl["annot.symbol"]):
        if dict_marker_symbol.get(marker, None) == None:
            dict_marker_symbol[marker] = list()
        dict_marker_symbol[marker].append(symbol)

    for key, list_val in dict_marker_symbol.items():
        list_unique_val = list(set(list_val))
        list_unique_val.remove(np.nan)
        if len(list_unique_val) > 1:
            unique_val = list_unique_val[-1]
        else:
            unique_val = "-".join(list_unique_val)
        if unique_val == "":
            unique_val = "None"
        dict_marker_symbol[key] = unique_val

    df_all_sample_cpgmin_indexed = df_annot_methyl.set_index("Sample_ID")
    df_all_sample_cpgmin_indexed = df_all_sample_cpgmin_indexed.rename(columns=dict_marker_symbol)
    hm = sns.clustermap(df_all_sample_cpgmin_indexed.T, xticklabels=True, yticklabels=True)
    labels = hm.ax_heatmap.get_xticklabels()

    list_sampleid = list()
    for label in labels:
        sampleid = label.get_text()
        list_sampleid.append(sampleid)

    df_merged["Sample_Age_Group"] = df_merged["Sample_Age"].apply(lambda x: "Old" if int(x) > 50 else "Young")
    list_sev = list(map(lambda x: df_merged.loc[x, "Severity"], list_sampleid))
    list_sex = list(map(lambda x: df_merged.loc[x, "Sample_Sex"], list_sampleid))
    list_age_group = list(map(lambda x: df_merged.loc[x, "Sample_Age_Group"], list_sampleid))
    dict_color_sev_visit = { 1: "g", 2: "yellow", 3: "orange", 4: "red"}
    dict_color_sex = {1:"b", 2:"r"}
    dict_color_age = {"Old":"grey", "Young": "k"}
    list_sev_color = list(map(dict_color_sev_visit.__getitem__, list_sev))
    list_sex_color = list(map(dict_color_sex.__getitem__, list_sex))
    list_age_color = list(map(dict_color_age.__getitem__, list_age_group))
    color_df = pd.DataFrame({"Sex": list_sex_color, "Age": list_age_color, "Severity": list_sev_color}, index=list_sampleid)

    hm = sns.clustermap(df_all_sample_cpgmin_indexed.T, col_colors=color_df, xticklabels=True, yticklabels=False, figsize=(20,5))
    hm.ax_row_dendrogram.set_visible(False)
    plt.show()
    plt.close()

# %%
def get_list_files_methylcpgmin(dir_methylcpgmin):
    list_files_methylcpgmin = glob.glob(f"{dir_methylcpgmin}/**/*pair_merged.methyl_cpg_min.tsv", recursive=True)

    return list_files_methylcpgmin

def get_list_methyl_sample(list_files_methylcpgmin):
    list_methyl_sample = list()
    for file_methylcpgmin in list_files_methylcpgmin:
        dir_methylcpgmin = os.path.basename(os.path.dirname(file_methylcpgmin))
        if dir_methylcpgmin == "HealthyControl":
            name_sample = os.path.basename(file_methylcpgmin).split(".")[0]
        else:
            name_sample = os.path.basename(dir_methylcpgmin)
        list_methyl_sample.append(name_sample)

    return list_methyl_sample

def get_dict_sample_severity(path_sev_info, list_methyl_sample):
    dict_sample_severity = dict()
    with open (path_sev_info, mode="r") as fr:
        _skiprow = fr.readline()
        for line in fr:
            record = line.rstrip("\n").split("\t")
            sampleid = record[0]
            severity_visit = record[-1]
            if sampleid in list_methyl_sample:
                dict_sample_severity[sampleid] = severity_visit

    return dict_sample_severity

# def get_dict_severity_group(path_sev_info, list_methyl_sample):
#     dict_severity_group = dict()
#     with open(path_sev_info, mode="r") as fr:
#         _skiprow = fr.readline()
#         for line in fr:
#             record = line.rstrip("\n").split("\t")
#             sampleid = record[0]
#             severity_visit = record[-1]

#             if dict_severity_group.get(severity_visit) == None:
#                 dict_severity_group[severity_visit] = list()
            
#             for methyl_sample in list_methyl_sample:
#                 if sampleid == methyl_sample:
#                     dict_severity_group[severity_visit].append(sampleid)

#     return dict_severity_group

# def remove_severity_group_no_content(dict_severity_group):
#     list_severity_group_no_content = list()
#     for key, list_val in dict_severity_group.items():
#         if list_val == list():
#             list_severity_group_no_content.append(key)
    
#     for sev_grp in list_severity_group_no_content:
#         dict_severity_group.pop(sev_grp)
            
#     return dict_severity_group    

def load_pickle(loadfilename):
    with gzip.open(loadfilename,'rb') as fr:
        data = pickle.load(fr)
    
    return data

def get_dict_mean_beta_all_samples(dict_cpgmin_all_samples):
    dict_meanbeta_all_samples = dict()
    for sampleid, dict_cpgmin in dict_cpgmin_all_samples.items():
        list_beta = list(map(float, dict_cpgmin.values()))
        mean_beta = np.mean(list_beta)
        dict_meanbeta_all_samples[sampleid] = mean_beta

    return dict_meanbeta_all_samples
# %%
if __name__ == "__main__":
    visit = "first"
    direction = "hyper"
    path_sev_info = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/9_clinical/Infectomics_Severity_Information_Methyl_20231106.tsv"
    dir_methylcpgmin = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/MethylCpGMin"
    list_drop_samples = ["C19-C045-V2",
                        "C19-C045-V3",
                        "C19-C047-V2",
                        "C19-C047-V3",
                        "C19-C050-V2",
                        "C19-C050-V3",
                        "C19-C051-V2",
                        "C19-C051-V3",
                        "C19-C052-V2",
                        "C19-C052-V3",
                        "C19-C053-V2",
                        "C19-C053-V3",
                        "C19-C055-V2",
                        "C19-C055-V3",
                        "C19-C056-V2",
                        'C19-C056-V3',
                        'C19-C060-V2',
                        'C19-C060-V3',
                        'C19-C061-V2',
                        'C19-C061-V3']
    infilename = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{visit}/{direction}/dictionary_marker_freqC_all_samples_20240102.pk.gz"
    outfilename = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{visit}/{direction}/mean_beta_by_severiy_table_20240102.tsv"
# %%
    main(path_sev_info, dir_methylcpgmin, list_drop_samples, infilename, outfilename)

# %%
annot_methyl = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231020/annotation/methyl_hyper_severe_mild_annotatr.tsv"
df_annot_methyl = pd.read_csv(annot_methyl, sep="\t")
df_annot_methyl["marker"] = df_annot_methyl["seqnames"].astype(str) + ":" + df_annot_methyl["start"].apply(lambda x: int(x)-1).astype(str)

dict_marker_symbol = dict()
for marker, symbol in zip(df_annot_methyl["marker"], df_annot_methyl["annot.symbol"]):
    if dict_marker_symbol.get(marker, None) == None:
        dict_marker_symbol[marker] = list()
    dict_marker_symbol[marker].append(symbol)

for key, list_val in dict_marker_symbol.items():
    list_unique_val = list(set(list_val))
    list_unique_val.remove(np.nan)
    if len(list_unique_val) > 1:
        unique_val = list_unique_val[-1]
    else:
        unique_val = "-".join(list_unique_val)
    if unique_val == "":
        unique_val = "None"
    dict_marker_symbol[key] = unique_val