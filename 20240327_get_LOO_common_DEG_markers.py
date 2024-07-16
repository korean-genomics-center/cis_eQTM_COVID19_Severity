# %% 
import os
import glob
import pandas as pd
import numpy as np
from collections import defaultdict
from itertools import chain
# %%
dir_deg = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/5_deg"
fc_thres = 1.3
padj_thres = 0.05
list_num_overlap_folds = list(range(1, 9, 1))

df_dmp = pd.read_csv("/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/Methyl_RNA_Relation/Methyl_RNA_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Mild_Sev_Visit1.LOO_common.Filtered.rho_75.methyl_direction_annot.20240326.tsv", sep = '\t')
# %%
path_deg = os.path.join(dir_deg, f"Visit*_Severe__Visit*_Mild_removed_loo*.tsv")
list_deg_results = glob.glob(path_deg)
list_visits = sorted(list(set(list(map(lambda x: os.path.basename(x).split("_")[0], list_deg_results)))))
dict_visit_deg_files = defaultdict(list)
for visit in list_visits:
    for deg_res in list_deg_results:
        if visit in deg_res:
            dict_visit_deg_files[visit].append(deg_res)

dict_visit_fold_list_common = dict()
dict_visit_list_fcs = dict()
for visit, list_deg_files_per_visit in dict(dict_visit_deg_files).items():
    dict_gene_list_fc_each_fold = defaultdict(list)
    for deg_file_per_visit in list_deg_files_per_visit:
        df_deg_per_visit = pd.read_csv(deg_file_per_visit, sep="\t")
        df_deg_per_visit["abslog2FoldChange"] = df_deg_per_visit["log2FoldChange"].apply(abs)
        df_deg_per_visit_sig = df_deg_per_visit[np.logical_and(df_deg_per_visit["abslog2FoldChange"] > fc_thres, df_deg_per_visit["padj"] < padj_thres)]
        for idx, row in df_deg_per_visit_sig.iterrows():
            dict_row = dict(row)
            gene = dict_row["ID"]
            fc = dict_row["log2FoldChange"]
            dict_gene_list_fc_each_fold[gene].append(fc)

    dict_overlap_genes_per_fold = dict()
    dict_visit_list_fcs[visit] = dict()
    for num_overlap_folds in list_num_overlap_folds:
        dict_gene_list_fc_all_folds = dict()
        for geneid, list_fc in dict_gene_list_fc_each_fold.items():
            cond_direction = (len(set(np.sign(list_fc))) == 1)
            cond_all_fold = (len(list_fc) == num_overlap_folds)
            if np.logical_and(cond_all_fold, cond_direction):
                dict_gene_list_fc_all_folds[geneid] = list_fc
        
        dict_visit_list_fcs[visit][num_overlap_folds] = dict_gene_list_fc_all_folds
        df_common_loo_deg_each_visit = pd.DataFrame.from_dict(dict_gene_list_fc_all_folds, orient="index")
        df_common_loo_deg_each_visit_reidx = df_common_loo_deg_each_visit.reset_index(drop=False)
        # df_common_loo_deg_each_visit_reidx.columns = ["GeneID"] + list(map(lambda x: f"log2FC(fold{x})", list(range(1, num_overlap_folds+1, 1))))
        path_common_loo_deg = os.path.join(dir_deg, f"LOO_Common_DEG_markers_visit_{visit}_overlap{num_overlap_folds}folds_20240327.tsv")
        # df_common_loo_deg_each_visit_reidx.to_csv(path_common_loo_deg, sep="\t", index=False)

        list_deg_common = df_common_loo_deg_each_visit_reidx["index"].to_list()
        dict_overlap_genes_per_fold[num_overlap_folds] = list_deg_common
    
    dict_visit_fold_list_common[visit] = dict_overlap_genes_per_fold

    list_fold = range(1, 9, 1)
    fold_cutoff = 7
    list_folds_over_cutoff = list(filter(lambda x : x >= fold_cutoff, list_fold))

    list_genes_over_cutoff = list(chain(*list(map(lambda x : dict_visit_fold_list_common["Visit1"].get(x, list()), list_folds_over_cutoff))))

    df_dmp_over_cutoff_with_RNA = df_dmp[df_dmp["RNA"].isin(list_genes_over_cutoff)]

    print(len(list_genes_over_cutoff), df_dmp_over_cutoff_with_RNA.shape)

    for x in df_dmp_over_cutoff_with_RNA["RNA"].unique():
        print(x.split("_")[-1])
# %%
for visit in dict_visit_fold_list_common.keys():
    from matplotlib import pyplot as plt
    plt.figure()
    plt.bar(range(1, 9, 1), list(map(lambda x : len(dict_visit_fold_list_common[visit].get(x, list())), range(1, 9, 1))))
    plt.title(visit)
    plt.xticks(range(1, 9, 1))
    plt.show()
    plt.close()

# %%
