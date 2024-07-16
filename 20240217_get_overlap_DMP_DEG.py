# %%
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd


# %%
def main(dir_dmp, list_file_dmp, dir_deg, list_file_deg, gene_chr_info, fcthres, qvalthres, dir_output, outfilename):
    dict_ensid_chr = get_dict_ensid_chr_gene_info(gene_chr_info)
    path_output = os.path.join(dir_output, outfilename)

    with open(path_output, mode="w") as fw:
        list_header = ["DMP_direction", "DEG_visit_number", "Gene_Symbol", "DMP_Markers", "Log2FC", "qValue"]
        fw.write("\t".join(list_header) + "\n")
        dict_dmp = get_dictionary_annotation_dmp(dir_dmp, list_file_dmp)
        dict_deg = get_dictionary_stat_deg(dir_deg, list_file_deg, dict_ensid_chr, fcthres=fcthres, qvalthres=qvalthres)
        for direction, dict_dmp_marker in dict_dmp.items():
            for visit, dict_deg_stat in dict_deg.items():
                for dmp_gene_symbol, set_markers in dict_dmp_marker.items():
                    for deg_gene, deg_stat in dict_deg_stat.items():
                        list_markers = "/".join(list(set_markers))
                        chr_dmp = list_markers.split("/")[0].split(":")[0]
                        chr_deg = deg_stat["chr"]
                        deg_gene_symbol = "_".join(deg_gene.split("_")[1:])
                        if (dmp_gene_symbol == deg_gene_symbol) and (chr_dmp == chr_deg):
                            fc = deg_stat["log2fc"]
                            qval = deg_stat["qval"]
                            list_content = [direction, visit, deg_gene, list_markers, fc, qval]
                            fw.write("\t".join(list_content) + "\n")

    df_overlap = sort_by_overlap_deg_fc(dir_output, outfilename)
    print(df_overlap)
# %%
def get_dict_ensid_chr_gene_info(gene_chr_info):
    dict_ensid_chr = dict()
    with open(gene_chr_info, mode="r") as fr:
        list_header = fr.readline().rstrip("\n").split("\t")
        idx_ensid = list_header.index("ensembl_gene_id")
        idx_loc = list_header.index("location")
        for line in fr:
            record = line.rstrip("\n").split("\t")
            ensid = record[idx_ensid]
            loc = record[idx_loc]
            chr = "chr" + re.split("q|p", loc)[0]
            dict_ensid_chr[ensid] = chr
    
    return dict_ensid_chr

def get_list_dmp(path_file_dmp):
    with open(path_file_dmp, mode="r") as fr:
        list_records = fr.readlines()
        list_record_strip = list(map(lambda x: x.rstrip(), list_records))
        
    return list_record_strip

def get_dictionary_annotation_dmp(dir_dmp, list_file_dmp):
    dict_pos_gene_all = dict()
    list_path_file_dmp = list(map(lambda x: os.path.join(dir_dmp, x), list_file_dmp))
    for path_file_dmp in list_path_file_dmp:
        # dmp_direction = os.path.basename(path_file_dmp).split("_")[1]
        dmp_direction = os.path.basename(path_file_dmp).split("_")[-3]
        table_dmp = pd.read_csv(path_file_dmp, sep="\t")
        # table_dmp["pos"] = table_dmp["seqnames"].astype(str) + ":" + table_dmp["start"].apply(lambda x: int(x)-1).astype(str)
        table_dmp["pos"] = table_dmp["Methyl"].apply(lambda x: ":".join(x.split("_")[:2]))
        list_pos_dmp = list(map(str, table_dmp["pos"].to_list()))
        list_annot_symbol = list(map(str, table_dmp["RNA"].to_list()))

        dict_pos_gene = defaultdict(set)
        for pos, gene in zip(list_pos_dmp, list_annot_symbol):
            if gene == "nan":
                continue
            else:
                dict_pos_gene[gene].add(pos)
        
        dict_pos_gene_all[dmp_direction] = dict(dict_pos_gene)
    
    return dict_pos_gene_all

def get_dictionary_stat_deg(dir_deg, list_file_deg, dict_ensid_chr, fcthres=1.5, qvalthres=0.05):
    list_path_file_deg = list(map(lambda x: os.path.join(dir_deg, x), list_file_deg))
    
    dict_stat_deg_all_visit = dict()
    for path_file_deg in list_path_file_deg:
        deg_visit_num = os.path.basename(path_file_deg).split("_")[0]
        df_deg = pd.read_csv(path_file_deg, sep="\t")
        df_deg_fc_cut = df_deg[abs(df_deg["log2FoldChange"]) > fcthres]
        df_deg_pval_cut = df_deg_fc_cut[df_deg_fc_cut["padj"] < qvalthres]

        dict_stat_deg_each_visit = dict()
        for _, rows in df_deg_pval_cut.iterrows():
            dict_rows = dict(rows)
            gene_id = dict_rows["ID"]
            ens_id = gene_id.split(".")[0]
            if ens_id in dict_ensid_chr.keys():
                chr = dict_ensid_chr[ens_id]
            else:
                chr = None
            fc = str(dict_rows["log2FoldChange"])
            qval = str(dict_rows["padj"])
            dict_stat_deg_each_visit[gene_id] = {"chr": chr, "log2fc": fc, "qval": qval}

        dict_stat_deg_all_visit[deg_visit_num] = dict_stat_deg_each_visit

    return dict_stat_deg_all_visit

def sort_by_overlap_deg_fc(dir_output, outfilename):
    path_output = os.path.join(dir_output, outfilename)
    df_out = pd.read_csv(path_output, sep="\t")
    df_out["absLog2FC"] = df_out["Log2FC"].apply(abs)
    df_out_sorted = df_out.sort_values(by=["DMP_direction", "DEG_visit_number", "absLog2FC"], ascending=False)
    df_out_sorted = df_out_sorted.drop(columns=["absLog2FC"])
    df_out_sorted = df_out_sorted.reset_index(drop=True)
    path_sorted_output = path_output.split(".tsv")[0] + "_sorted.tsv"
    df_out_sorted.to_csv(path_sorted_output, sep="\t", index=False)

    return df_out_sorted
# %%
if __name__ == "__main__":
    visit = "first"
    list_file_deg = ["Visit1_Severe__Visit1_Mild_20240327.tsv", "Visit2_Severe__Visit2_Mild_20240327.tsv"]
    dir_deg = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/5_deg"
    list_file_dmp = [f"Methyl_RNA_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Mild_Sev_Visit1.LOO_common.Filtered.rho_75.methyl_direction_annot.hyper.20240326.tsv", "Methyl_RNA_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Mild_Sev_Visit1.LOO_common.Filtered.rho_75.methyl_direction_annot.hypo.20240326.tsv"]
    dir_dmp = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/Methyl_RNA_Relation"
    gene_chr_info = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Resources/Data/hgnc_gene_info.txt"
    dir_output = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/DMPDEG/{visit}"
    os.makedirs(dir_output, exist_ok=True)
    fcthres = 1.3
    qvalthres = 0.05
    outfilename = f"table_deg_dmp_overlap_abslog2fc_{fcthres}_qval_{qvalthres}_20240327.tsv"
    main(dir_dmp, list_file_dmp, dir_deg, list_file_deg, gene_chr_info, fcthres, qvalthres, dir_output, outfilename)
    
# %%
# if __name__ == "__main__":
#     visit = "first"
#     list_file_deg = ["Visit1_Severe__Visit1_Mild_20240326.tsv", "Visit2_Severe__Visit2_Mild_20240326.tsv"]
#     dir_deg = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/5_deg"
#     list_file_dmp = [f"methyl_hyper_severe_mild_{visit}Visit_annotatr.tsv", f"methyl_hypo_severe_mild_{visit}Visit_annotatr.tsv"]
#     dir_dmp = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231211/severe_mild_{visit}Visit/annotation"
#     gene_chr_info = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Resources/Data/hgnc_gene_info.txt"
#     dir_output = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/DMPDEG/{visit}"
#     os.makedirs(dir_output, exist_ok=True)
#     fcthres = 1.0
#     qvalthres = 0.05
#     outfilename = f"table_deg_dmp_overlap_abslog2fc_{fcthres}_qval_{qvalthres}_20240326.tsv"
#     main(dir_dmp, list_file_dmp, dir_deg, list_file_deg, gene_chr_info, fcthres, qvalthres, dir_output, outfilename)
# %%
