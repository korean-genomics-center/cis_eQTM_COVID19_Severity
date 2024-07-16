#%%
import os
import pandas as pd
import numpy as np
from itertools import chain
from collections import Counter
from glob import glob

##########################
col_p = "corr_fdr"
padj_cutoff = 0.05
fc_cutoff = 1.3
corr_cutoffs = np.linspace(0.5, 0.75, 6)
loo_overlap_cutoffs = range(1, 8+1, 1)


path_deg_loo_format = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/5_deg/Visit4_Severe__Visit4_Mild_removed_loo_*_20240327.tsv"
path_dmp_correlated_rna = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/Methyl_RNA_Relation/Methyl_RNA_Correlation.Infectomics.Visit3_4.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Mild_Sev_VisitLast.LOO_common.cis.distance_cf.1000000.20240326.tsv"

list_files_loo = glob(path_deg_loo_format)

#########################
#%%
# RNA
dict_genes_deg_significant_of_all_loo = dict()
dict_genes_deg_significant_of_all_loo_sample_wo = dict()
for deg_file in list_files_loo:
    table_deg_tmp = pd.read_csv(deg_file, sep = '\t')
    table_deg_tmp_sig = table_deg_tmp[np.logical_and(table_deg_tmp["log2FoldChange"].abs() > fc_cutoff, table_deg_tmp["padj"] < padj_cutoff)]
    for geneid, fc in zip(table_deg_tmp_sig["ID"].to_list(), table_deg_tmp_sig["log2FoldChange"].to_list()):
        if dict_genes_deg_significant_of_all_loo.get(geneid) == None:
            dict_genes_deg_significant_of_all_loo[geneid] = list()
            dict_genes_deg_significant_of_all_loo_sample_wo[geneid] = list()
        dict_genes_deg_significant_of_all_loo[geneid].append(fc)
        dict_genes_deg_significant_of_all_loo_sample_wo[geneid].append(deg_file.split('_')[-2])

for geneid in list(dict_genes_deg_significant_of_all_loo.keys()):
    directions = list(map(lambda x : x < 0, dict_genes_deg_significant_of_all_loo[geneid]))
    is_consistent = sum(directions) == 0 or sum(directions) == len(directions)
    if not is_consistent:
        print(geneid)
        dict_genes_deg_significant_of_all_loo.pop(geneid)
    
count_deg_gene = dict(zip(dict_genes_deg_significant_of_all_loo.keys(), list(map(lambda g : len(dict_genes_deg_significant_of_all_loo[g]), dict_genes_deg_significant_of_all_loo.keys()))))

dict_genes_over_loo_overlap_threshold = dict()
for overlap_cutoff in loo_overlap_cutoffs:
    dict_genes_over_loo_overlap_threshold[overlap_cutoff] = list(filter(lambda gene : count_deg_gene[gene] >= overlap_cutoff, count_deg_gene.keys()))
    
    
dict_genes_equal_loo_overlap_threshold = dict()
for overlap_cutoff in loo_overlap_cutoffs:
    dict_genes_equal_loo_overlap_threshold[overlap_cutoff] = list(filter(lambda gene : count_deg_gene[gene] == overlap_cutoff, count_deg_gene.keys()))

#%%
# Methyl ~ RNA
table_dmp_correlated_rna = pd.read_csv(path_dmp_correlated_rna, sep = '\t')
table_dmp_correlated_rna_fdr_sig = table_dmp_correlated_rna[table_dmp_correlated_rna[col_p] < padj_cutoff]
#%%
dict_genes_over_corr_threshold = dict()
for corr_cutoff in corr_cutoffs:
    table_dmp_correlated_rna_fdr_sig_corr_sig = table_dmp_correlated_rna_fdr_sig[table_dmp_correlated_rna_fdr_sig["corr_rho"].abs() > corr_cutoff]
    dict_genes_over_corr_threshold[corr_cutoff] = list(table_dmp_correlated_rna_fdr_sig_corr_sig["RNA"].unique())
    
# %%
# LOO Stats
import seaborn as sns
from matplotlib import pyplot as plt

plt.rcParams["font.size"] = 13
plt.bar(loo_overlap_cutoffs, list(map(lambda cf : len(dict_genes_equal_loo_overlap_threshold[cf]), loo_overlap_cutoffs)), color = "gray")
for loo_cf in loo_overlap_cutoffs:
    cnt_genes = len(dict_genes_equal_loo_overlap_threshold[loo_cf])
    plt.text(loo_cf, cnt_genes, cnt_genes, ha = "center", va = "bottom")
plt.xlabel("Number of Overlaps on DEG-LOO\n(LOO on Severe, 8 samples)")
plt.ylabel("Number of Genes")
plt.tight_layout()
#%%
total_loo_samples = list(map(lambda x : x.split('_')[-2], list_files_loo))
table_loo_sample_specific = pd.DataFrame(columns = total_loo_samples, index = list(dict_genes_deg_significant_of_all_loo_sample_wo.keys())).fillna(0)
for gene, list_sample in dict_genes_deg_significant_of_all_loo_sample_wo.items():
    table_loo_sample_specific.loc[gene, list_sample] = 1
table_loo_sample_specific["discovery"] = table_loo_sample_specific.apply(sum, axis = 1)
table_loo_sample_specific = table_loo_sample_specific.sort_values(by = "discovery", ascending=False)
dict_gene_first_cnt = dict()
for loo_cnt in sorted(table_loo_sample_specific["discovery"].unique()):
    ind_first_loo_cnt = list(table_loo_sample_specific["discovery"]).index(loo_cnt)
    dict_gene_first_cnt[loo_cnt] = list(table_loo_sample_specific.index)[ind_first_loo_cnt]
#%%
table_loo_sample_specific = table_loo_sample_specific.drop(columns = ["discovery"])
#%%
plt.figure(figsize = (5, 12))
sns.heatmap(table_loo_sample_specific,cbar_kws = {"label":"Whether marker discovered (True/False)"})
plt.yticks([])
for loo_cnt, gene_first in dict_gene_first_cnt.items():
    if loo_cnt == len(total_loo_samples):continue
    plt.axhline(list(table_loo_sample_specific.index).index(gene_first), linestyle = "--", linewidth = 1, color = 'r')
plt.xlabel("Excluded Sample")
plt.ylabel("DEG significant genes (Sorted by Number of LOO-Overlaps)")
#%%
# Correlation Stats
list_num_dmps_with_correlation_between = list()
for ind, corr_cf in enumerate(corr_cutoffs):
    num_dmps_over_corr_cf = len(dict_genes_over_corr_threshold[corr_cf])
    if ind == len(corr_cutoffs)-1:
        num_dmps_over_next_corr_cf = 0
    else:
        num_dmps_over_next_corr_cf = len(dict_genes_over_corr_threshold[corr_cutoffs[ind+1]])
    list_num_dmps_with_correlation_between.append(num_dmps_over_corr_cf - num_dmps_over_next_corr_cf)
plt.bar(range(1, len(list_num_dmps_with_correlation_between)+1, 1), list_num_dmps_with_correlation_between, color = "gray", width = 1)
for ind, corr_cf in enumerate(corr_cutoffs, start = 1):
    plt.text(ind, list_num_dmps_with_correlation_between[ind-1], list_num_dmps_with_correlation_between[ind-1], ha = "center", va = "bottom")
for ind in np.arange(1.5, len(list_num_dmps_with_correlation_between), 1):
    plt.axvline(ind, linestyle = "--", linewidth = 1, color = 'r')
plt.xticks(np.linspace(0.5, len(list_num_dmps_with_correlation_between)+0.5, len(list_num_dmps_with_correlation_between)+1), list(corr_cutoffs) + ["~"])
plt.xlabel("Cutoff Range of DMS~RNA Correlation")
plt.ylabel("Number of Genes")
plt.show()
#%%
# Gene Count Heatmap
table_gene_overlap_counter_per_thresholds = pd.DataFrame(columns = list(loo_overlap_cutoffs), index = reversed(list(corr_cutoffs)))

for loo_cf in loo_overlap_cutoffs:
    for corr_cf in corr_cutoffs:
        list_genes_over_loo_cf = dict_genes_over_loo_overlap_threshold[loo_cf]
        list_genes_over_corr_cf = dict_genes_over_corr_threshold[corr_cf]
        list_genes_overlap_between_two = list(set(list_genes_over_corr_cf) & set(list_genes_over_loo_cf))
        table_gene_overlap_counter_per_thresholds.loc[corr_cf, loo_cf] = len(list_genes_overlap_between_two)
        
sns.heatmap(data = table_gene_overlap_counter_per_thresholds.astype(float), annot = True, fmt = ".0f")
plt.xlabel("Number of DEG overlaps during 8 time of LOO")
plt.ylabel("Correlation cutoffs between\nLOO-common DMSs and Gene expression")
plt.title("Number of Genes overlap\n(DMSs-Gene Correlation & LOO DEG)")
plt.show()
# %%
# Sig Relation
import math

corr_cutoff = 0.5
loo_cnt_cutoff = 6

table_rna_methylcorr_sig = table_dmp_correlated_rna_fdr_sig[table_dmp_correlated_rna_fdr_sig["corr_rho"].abs() > corr_cutoff]
table_rna_methylcorr_sig = table_rna_methylcorr_sig[table_rna_methylcorr_sig["RNA"].isin(set(dict_genes_over_corr_threshold[corr_cutoff]) & set(dict_genes_over_loo_overlap_threshold[loo_cnt_cutoff]))]
table_rna_methylcorr_sig["-log10(p)"] = table_rna_methylcorr_sig["corr_p"].apply(lambda x : -math.log10(x))

# table_rna_methylcorr_sig["Mean_"]

sns.scatterplot(data = table_rna_methylcorr_sig, x = "corr_rho", y = "-log10(p)", s = 9, color = "k")
plt.xlabel("Spearman Rho")
plt.axhline(table_rna_methylcorr_sig["-log10(p)"].min(), linestyle = "--", linewidth = 1, color = 'gray')
# plt.text(0, table_rna_methylcorr_sig["-log10(p)"].min(), f'Minimum p-value : {format(table_rna_methylcorr_sig["corr_p"].max(), ".2E")}', ha = "center", va = "bottom")
plt.show()

# table_rna_methylcorr_sig.to_csv(f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/Methyl_RNA_Relation/Methyl_RNA_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Mild_Sev_Visit1.LOO_common.cis.distance_cf.1000000.Overlap_DEG.corr_{corr_cutoff}.log2fc_{fc_cutoff}.loo_{loo_cnt_cutoff}.20240326.tsv", sep = '\t', index = False)
# %%
