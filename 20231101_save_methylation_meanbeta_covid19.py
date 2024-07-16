# %%
import glob
import gzip
import os
import pickle

import numpy as np
import pandas as pd


# %%
def main(path_marker_info, path_sev_info, dir_methylcpgmin, infilename, outfilename):
    list_files_methylcpgmin = get_list_files_methylcpgmin(dir_methylcpgmin)
    list_methyl_sample = get_list_methyl_sample(list_files_methylcpgmin)
    dict_cpgmin_all_samples = load_pickle(infilename)   
    list_markers = get_marker_info(path_marker_info)  
    dict_meanbeta_all_samples = get_dict_mean_beta_all_samples(dict_cpgmin_all_samples, list_markers)
    dict_sample_severity = get_dict_sample_severity(path_sev_info, list_methyl_sample)
    save_meanbeta_cpgmin(dict_sample_severity, dict_meanbeta_all_samples, outfilename)

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

def load_pickle(loadfilename):
    with gzip.open(loadfilename,'rb') as fr:
        data = pickle.load(fr)
    
    return data

def get_marker_info(path_marker_info, marker_col="DMP_Markers"):
    df_marker = pd.read_csv(path_marker_info, sep="\t")
    list_dmp = df_marker[marker_col].to_list()
    list_markers= list()
    for dmp in list_dmp:
        dmpsplit = dmp.split("/")
        list_markers.extend(dmpsplit)
    
    return list_markers

def get_dict_mean_beta_all_samples(dict_cpgmin_all_samples, list_markers):
    dict_meanbeta_all_samples = dict()
    for sampleid, dict_cpgmin in dict_cpgmin_all_samples.items():
        if list_markers == []:
            list_beta = list(map(float, dict_cpgmin.values()))
        else:
            list_beta = list()
            for key, val in dict_cpgmin.items():
                pos = key[0] + ":" + key[1]
                if pos in list_markers:
                    list_beta.append(val)
        list_beta = list(map(float, list_beta))
        print(len(list_beta))
        mean_beta = np.mean(list_beta)
        dict_meanbeta_all_samples[sampleid] = mean_beta

    return dict_meanbeta_all_samples

def save_meanbeta_cpgmin(dict_sample_severity, dict_meanbeta_all_samples, outfilename):
    with open(outfilename, mode="w") as fw:
        list_columns = ["Sample_ID", "Severity", "Visit", "Severity_Visit", "MeanBeta"]
        header = "\t".join(list_columns)
        fw.write(header + "\n")
        for sampleid, severity_visit in dict_sample_severity.items():
            severity = severity_visit.split("_")[0]
            visit = severity_visit.split("_")[1]
            meanbeta = str(dict_meanbeta_all_samples[sampleid])
            fw.write("\t".join([sampleid, severity, visit, severity_visit, meanbeta]) + "\n")
# %%
if __name__ == "__main__":
    direction = "hyper"
    visit = "first"
    path_marker_info = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/DMPDEG/{visit}/table_deg_dmp_overlap_abslog2fc_1.3_qval_0.05_20240229_sorted.tsv"
    path_sev_info = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/9_clinical/Infectomics_Severity_Information_Methyl_20240102.tsv"
    dir_methylcpgmin = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/MethylCpGMin"
    infilename = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{visit}/{direction}/dictionary_marker_freqC_all_samples_20240220.pk.gz"
    outfilename = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{visit}/{direction}/mean_beta_by_severity_table_dmpdeg_20240229.tsv"
    # outfilename = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{visit}/{direction}/mean_beta_by_severity_scale_in_detail_table_20240102.tsv"
    main(path_marker_info, path_sev_info, dir_methylcpgmin, infilename, outfilename)
# %%
