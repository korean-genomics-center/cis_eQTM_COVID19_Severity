# %%
import glob
import gzip
import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (accuracy_score, classification_report,
                             confusion_matrix, roc_auc_score, roc_curve)
from sklearn.model_selection import train_test_split


# %%
def main(path_sev_info, dir_methylcpgmin, list_drop_samples, infilename):
    list_files_methylcpgmin = get_list_files_methylcpgmin(dir_methylcpgmin)
    list_methyl_sample = get_list_methyl_sample(list_files_methylcpgmin)
    list_methyl_sample = list(filter(lambda x: "C19-C" in x, list_methyl_sample))
    list_methyl_sample = list(filter(lambda x: x not in list_drop_samples , list_methyl_sample))
    list_methyl_sample = list(filter(lambda x: "L1" not in x, list_methyl_sample))
    dict_cpgmin_all_samples = load_pickle(infilename)     
    dict_cpgmin_all_samples_marker_fixed = dict()
    for sampleid, dict_cpgs in dict_cpgmin_all_samples.items():
        if sampleid in list_methyl_sample:
            dict_cpgs_marker_fixed = dict()
            for marker, freqC in dict_cpgs.items():
                fixed_marker = ":".join(marker)
                dict_cpgs_marker_fixed[fixed_marker] = freqC
            
            dict_cpgmin_all_samples_marker_fixed[sampleid] = dict_cpgs_marker_fixed

    df_cpgmin_all_sample = pd.DataFrame(dict_cpgmin_all_samples_marker_fixed)
    df_all_sample_cpgmin = df_cpgmin_all_sample.T.reset_index(drop=False).rename(columns={"index": "Sample_ID"})
    df_sev = pd.read_csv(path_sev_info, sep="\t")
    df_merged = pd.merge(df_all_sample_cpgmin, df_sev, on="Sample_ID", how="inner")
    df_merged = df_merged.set_index("Sample_ID")

    #####
    #SEX
    df_merged_per_subj = df_merged[df_merged.index.str.endswith("V1")]
    df_grp_sev_sex = df_merged_per_subj.groupby("Severity_group")["Sample_Sex"].apply(np.array)

    from collections import Counter
    dict_sex = {1:"M", 2:"F"}
    list_sex_num_mild = list(dict(Counter(df_grp_sev_sex["Mild"])).keys())
    list_sex_mild = list(map(lambda x: dict_sex[x], list_sex_num_mild))

    list_cnt_mild = list(dict(Counter(df_grp_sev_sex["Mild"])).values())
    list_cnt_sev = list(dict(Counter(df_grp_sev_sex["Severe"])).values())

    fig, ax = plt.subplots(figsize=(5,5),layout="constrained")

    N = 2
    ind = np.arange(N)
    width = 0.3

    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#D68C78', edgecolor='k',label='Severe'), Patch(facecolor='#ECECEC', edgecolor='k',label='Mild')]
    ax.bar(ind, list_cnt_mild, width, label="Mild", color="#aaaaaa")
    ax.bar(ind+width, list_cnt_sev, width, label="Severe", color="#D68C78")
    for i in ax.containers:
        ax.bar_label(i,)
    plt.xticks(ind+width/2, list_sex_mild)
    # plt.title("Mild (N=41)")
    plt.xlabel("Sex", fontdict={"fontsize":12})
    plt.ylabel("No. Patients", fontdict={"fontsize":12})
    plt.title("Mild(N=41) vs Severe(N=10)")
    plt.legend(handles=legend_elements, loc="upper left")
    plt.show()
    plt.savefig("/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/9_clinical/BaselineCharacteristics/sex_distributions.png", dpi=600)
    plt.close()

    #AGE
    df_grp_sev_age = df_merged_per_subj.groupby("Severity_group")["Sample_Age"].apply(np.array)
    plt.hist(df_grp_sev_age["Mild"], bins=10, color="#aaaaaa", alpha=0.5)
    plt.vlines(x=np.mean(df_grp_sev_age["Mild"]), ymin=0, ymax=7, color="k")


    plt.hist(df_grp_sev_age["Severe"], bins=10, color="#D68C78", alpha=0.5)
    plt.vlines(x=np.mean(df_grp_sev_age["Severe"]), ymin=0, ymax=7, color="k")

    plt.text(x=np.mean(df_grp_sev_age["Mild"])-17, y=6, s=f"Mean(Mild)={round(np.mean(df_grp_sev_age['Mild']), 3)}", color="k")
    plt.text(x=np.mean(df_grp_sev_age["Severe"])+1, y=2, s=f"Mean(Severe)={round(np.mean(df_grp_sev_age['Severe']), 3)}", color="k")

    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#D68C78', edgecolor='k',label='Severe'), Patch(facecolor='#ECECEC', edgecolor='k',label='Mild')]

    plt.legend(handles=legend_elements, loc='upper right')

    plt.xlabel("Age")
    plt.ylabel("Count")
    plt.title("Mild (N=41) vs Severe (N=10)")
    plt.show()
    plt.savefig("/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/9_clinical/BaselineCharacteristics/age_distributions.png", dpi=600)
    plt.close()
    ####

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

if __name__ == "__main__":
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
    infilename = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{direction}/dictionary_marker_freqC_all_samples_20231101.pk.gz"
# %%
    main(path_sev_info, dir_methylcpgmin, list_drop_samples, infilename)