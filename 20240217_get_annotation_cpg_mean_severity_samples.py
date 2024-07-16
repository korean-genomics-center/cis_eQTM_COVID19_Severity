# %%
import glob
import gzip
import os
import pickle

import numpy as np
import pandas as pd


# %%
def main(path_sev_info, dir_methylcpgmin, methylannot, num_threads, infilename, outfilename, final_output):
    list_files_methylcpgmin = get_list_files_methylcpgmin(dir_methylcpgmin)
    list_methyl_sample = get_list_methyl_sample(list_files_methylcpgmin)
    dict_cpgmin_all_samples = load_pickle(infilename)
    dict_sample_severity = get_dict_sample_severity(path_sev_info, list_methyl_sample)
    dict_marker_genename = annotate_gene_name_marker(methylannot)
        
    with open(outfilename, mode="w") as fw:
        fw.write("\t".join(["pos", "gene_name", "cpgmin", "severity"]) + "\n")
        for sample_id, dict_marker in dict_cpgmin_all_samples.items():
            try:
                severity_visit = dict_sample_severity[sample_id]
            except:
                print(f"Skipping {sample_id}... No value")
            # severity_only = severity_visit.split("_")[0]
            for marker_name, cpg in dict_marker.items():
                pos = marker_name[0] + ":" + marker_name[1]
                genename = dict_marker_genename[pos]
                fw.write("\t".join([pos, genename, cpg, severity_visit]) + "\n")
            
    list_dict_cpg_mild_first = list()
    list_dict_cpg_mild_last = list()
    list_dict_cpg_severe_first = list()
    list_dict_cpg_severe_last = list()
    df_pos_cpg_sev = pd.read_csv(outfilename, sep="\t")
    df_pos_sev_grby_mean_beta = df_pos_cpg_sev.groupby(["pos", "gene_name", "severity"])["cpgmin"].apply(np.mean).reset_index(drop=False)
    for idx, row in df_pos_sev_grby_mean_beta.iterrows():
        dict_row = dict(row)
        if dict_row["severity"] == "Mild_First":
            list_dict_cpg_mild_first.append(dict_row)
        elif dict_row["severity"] == "Mild_Last":
            list_dict_cpg_mild_last.append(dict_row)
        elif dict_row["severity"] == "Severe_First":
            list_dict_cpg_severe_first.append(dict_row)
        elif dict_row["severity"] == "Severe_Last":
            list_dict_cpg_severe_last.append(dict_row)
    
    with open(final_output, mode="w") as fw:
        fw.write("\t".join(["pos", "gene_name", "mild_first_mean_beta", "mild_last_mean_beta", "severe_first_mean_beta", "severe_last_mean_beta", "diff_first_mean_beta", "diff_last_mean_beta", "delta_diff_mean_beta"]) + "\n")
        list_pos = list(map(lambda x: x["pos"], list_dict_cpg_mild_first))
        for pos in list_pos:
            for dict_cpg_mild_first, dict_cpg_mild_last, dict_cpg_severe_first, dict_cpg_severe_last in zip(list_dict_cpg_mild_first, list_dict_cpg_mild_last, list_dict_cpg_severe_first, list_dict_cpg_severe_last):
                if dict_cpg_mild_first["pos"] == pos and dict_cpg_severe_first["pos"] == pos:
                    fin_pos = dict_cpg_mild_first["pos"]
                    fin_name = dict_cpg_mild_first["gene_name"]
                    fin_mild_first_mean_cpg = str(dict_cpg_mild_first["cpgmin"])
                    fin_mild_last_mean_cpg = str(dict_cpg_mild_last["cpgmin"])
                    fin_sev_first_mean_cpg = str(dict_cpg_severe_first["cpgmin"])
                    fin_sev_last_mean_cpg = str(dict_cpg_severe_last["cpgmin"])
                    fin_diff_first_mean_cpg = float(fin_sev_first_mean_cpg) - float(fin_mild_first_mean_cpg)
                    fin_diff_last_mean_cpg = float(fin_sev_last_mean_cpg) - float(fin_mild_last_mean_cpg)
                    fin_diff_diff_last_first_mean_cpg = float(fin_diff_last_mean_cpg) - float(fin_diff_first_mean_cpg)
                    fw.write("\t".join([fin_pos, fin_name, fin_mild_first_mean_cpg, fin_mild_last_mean_cpg, fin_sev_first_mean_cpg, fin_sev_last_mean_cpg, fin_diff_first_mean_cpg, fin_diff_last_mean_cpg, fin_diff_diff_last_first_mean_cpg]) + "\n")

    df_fin = pd.read_csv(final_output, sep="\t")
    df_fin["abs_diff_first_mean_beta"] = df_fin["diff_first_mean_beta"].apply(lambda x: abs(float(x)))
    df_fin["abs_diff_last_mean_beta"] = df_fin["diff_last_mean_beta"].apply(lambda x: abs(float(x)))
    df_fin_diff_first_sorted = df_fin.sort_values(by=["abs_diff_first_mean_beta"], ascending=False)
    df_fin_diff_last_sorted = df_fin.sort_values(by=["abs_diff_last_mean_beta"], ascending=False)
    
    list_gene_first = list(filter(lambda x: x != "None", df_fin_diff_first_sorted["gene_name"].to_list()))[:30]
    list_gene_last = list(filter(lambda x: x != "None", df_fin_diff_last_sorted["gene_name"].to_list()))[:30]
    intsc = list(set(list_gene_first).intersection(set(list_gene_last)))
    for gene in intsc:
        print(gene)
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

def get_dict_mean_beta_all_samples(dict_cpgmin_all_samples):
    dict_meanbeta_all_samples = dict()
    for sampleid, dict_cpgmin in dict_cpgmin_all_samples.items():
        list_beta = list(map(float, dict_cpgmin.values()))
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
            meanbeta = dict_meanbeta_all_samples[sampleid]
            fw.write("\t".join([sampleid, severity, visit, severity_visit, meanbeta]) + "\n")
            
def annotate_gene_name_marker(annot_methyl):
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
    
    return dict_marker_symbol

# %%
if __name__ == "__main__":
    direction = "hyper"
    visit = "first"
    path_sev_info = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/9_clinical/Infectomics_Severity_Information_Methyl_20240102.tsv"
    dir_methylcpgmin = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/MethylCpGMin"
    methylannot = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231211/severe_mild_{visit}Visit/annotation/methyl_{direction}_severe_mild_{visit}Visit_annotatr.tsv"
    infilename = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{visit}/{direction}/dictionary_marker_freqC_all_samples_20240220.pk.gz"
    outfilename = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{visit}/{direction}/all_beta_by_severity_table_20240220.tsv"
    final_output = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{visit}/{direction}/diff_beta_each_marker_by_severity_table_20240220.tsv"
    # outfilename = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{visit}/{direction}/mean_beta_by_severity_scale_in_detail_table_20240102.tsv"
    # main(path_sev_info, dir_methylcpgmin, methylannot, num_threads, infilename, outfilename, final_output)
# %%
