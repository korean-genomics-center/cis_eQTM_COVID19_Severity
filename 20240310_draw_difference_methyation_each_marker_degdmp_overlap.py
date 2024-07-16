# %%
import glob
import gzip
import os
import pickle
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import levene, mannwhitneyu, ranksums, ttest_ind
from statannot import add_stat_annotation
from statsmodels.stats.multitest import fdrcorrection

# %%
path_sev_info = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/9_clinical/Infectomics_Severity_Information_Methyl_20240102.tsv"
dir_methylcpgmin = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/MethylCpGMin"
infilenamehyper = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/first/hyper/dictionary_marker_freqC_all_samples_20240102.pk.gz"
infilenamehypo = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/first/hypo/dictionary_marker_freqC_all_samples_20240102.pk.gz"
file_deg_dmp_overlap = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/DMPDEG/first/table_deg_dmp_overlap_abslog2fc_1.3_qval_0.05_20240229_sorted.tsv"
list_drop_samples = [
                    "C19-C009-V1",
                    "C19-C011-V1",
                    "C19-C016-V1",
                    "C19-C021-V1",
                    "C19-C022-V1",
                    "C19-C045-V2",
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
outdir = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/InfectomicsPaper1/boxplot_per_marker"
os.makedirs(outdir, exist_ok=True)

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

def load_pickle(loadfilename):
    with gzip.open(loadfilename,'rb') as fr:
        data = pickle.load(fr)
    
    return data

def get_dict_deg_dmp_overlap_markers(file_deg_dmp_overlap, col_marker="DMP_Markers", col_gene = "Gene_Symbol", col_direction="DMP_direction", direction="hyper"):
    df_deg_dmp_overlap = pd.read_csv(file_deg_dmp_overlap, sep="\t")
    df_deg_dmp_overlap_direction = df_deg_dmp_overlap[df_deg_dmp_overlap[col_direction] == direction]
    list_marker = df_deg_dmp_overlap_direction[col_marker].to_list()
    list_genesym = df_deg_dmp_overlap_direction[col_gene].to_list()
    dict_all_markers = dict()
    for marker, genesym in zip(list_marker, list_genesym):
        all_markers = marker.split("/")
        for splited_marker in all_markers:
            dict_all_markers[splited_marker] = genesym

    return dict_all_markers

def get_dataframe_methylation_beta_samples(infilename, list_selected_samples):
    dict_cpgmin_all_samples = load_pickle(infilename)     
    dict_cpgmin_all_samples_marker_fixed = dict()
    for sampleid, dict_cpgs in dict_cpgmin_all_samples.items():
        dict_cpgs_marker_fixed = dict()
        for marker, freqC in dict_cpgs.items():
            fixed_marker = ":".join(marker)
            dict_cpgs_marker_fixed[fixed_marker] = freqC
        dict_cpgmin_all_samples_marker_fixed[sampleid] = dict_cpgs_marker_fixed

    df_beta = pd.DataFrame.from_dict(dict_cpgmin_all_samples_marker_fixed).astype(float)
    df_beta_transposed = df_beta.T.reset_index(drop=False).rename(columns={"index": "Sample_ID"})
    df_sev = pd.read_csv(path_sev_info, sep="\t")
    df_beta_sev = pd.merge(df_beta_transposed, df_sev, how="inner", on="Sample_ID")
    df_beta_sev["Visit"] = df_beta_sev["Visit"].fillna("Healthy")
    df_beta_sev["Visit_order"] = df_beta_sev["Visit_order"].fillna("Visit0")
    df_beta_sev["Severity_group"] = df_beta_sev["Severity_group"].fillna("Healthy")
    df_beta_sev["Subject_ID"] = df_beta_sev["Subject_ID"].fillna(df_beta_sev["Sample_ID"])
    df_beta_sev = df_beta_sev[df_beta_sev["Visit_order"]!="Visit6"]
    df_beta_set_idx = df_beta_sev.set_index("Sample_ID")
    # df_beta_set_idx_select = df_beta_set_idx.loc[list_selected_samples, :]
    df_beta_set_idx_select = df_beta_set_idx[df_beta_set_idx.index.isin(list_selected_samples)]
    
    return df_beta_set_idx_select

def meandiff(a, b):
    deltamean = np.mean(b) - np.mean(a)

    return deltamean

def get_dict_combinatorial_beta_sev_grp_compared(dict_beta_sev_grp, list_combination_sev_grps_filtered):
    from collections import defaultdict
    dict_beta_combination_groups = defaultdict(list)
    for combination_groups in list_combination_sev_grps_filtered:
        for grp in combination_groups:
            beta = dict_beta_sev_grp[grp]
            dict_beta_combination_groups[combination_groups].append(beta)
    
    return dict_beta_combination_groups

def get_stat_results_marker_beta_across_groups(outdir, dict_marker_beta_group):
    path_stattest = os.path.join(outdir, f"stattest_{marker}_{gensym}.tsv")
    dict_group_stats = dict()
    for groups, betas in dict_marker_beta_group.items():
        deltamean = meandiff(*betas)
        # lstat, lpval = levene(*betas, center="mean")
        # if lpval < 0.05:
        #     stat, pval = ttest_ind(*betas, equal_var=False)
        # else:
        #     stat, pval = ttest_ind(*betas, equal_var=True)
        stat, pval = ranksums(*betas)
        dict_group_stats[groups] = {"Delta": deltamean, "Statistics": stat, "P-value": pval}
    result_test = pd.DataFrame.from_dict(dict_group_stats, orient="index").reset_index(drop=False)
    result_test = result_test.rename(columns={"level_0": "Group1", "level_1": "Group2"})
    list_sigtrue = fdrcorrection(result_test["P-value"])[0]
    list_padj = fdrcorrection(result_test["P-value"])[1]
    result_test["Padj"] = list_padj
    result_test["TestSignificant"] = list_sigtrue
    result_test.to_csv(path_stattest, sep="\t", index=False)

    return result_test

def draw_boxplot_beta_each_marker(outdir, marker, gensym, df_beta_marker, result_test_marker, list_severity):
    path_boxplot = os.path.join(outdir, f"boxplot_{marker}_{gensym}.png")
    ## for sample count 
    # dict_cnt = df_beta_marker.groupby(["Severity_visit"]).apply(len).to_dict()
    
    # list_severity_cnt = list()
    # for sev in list_severity:
    #     if sev in dict_cnt.keys():
    #         cnt = dict_cnt[sev]
    #         if sev.split("_")[0] == "Healthy":
    #             sev = sev.split("_")[0] + "\n" + "Controls"
    #         else:
    #             sev = sev.split("_")[0] + "\n" + "Group"
    #         sev_nline = '\n'.join(sev.split('_'))
    #         list_severity_cnt.append(sev_nline)
    xticklabels = list(map(lambda x: "\n".join(x.split("_")), list_severity))
    xticklabels = list(map(lambda x: x.replace("Convalescent", "Conval"), xticklabels))
    
    title = f"{marker} ({'_'.join(gensym.split('_')[1:])})"
            
    palette = dict()
    for sev in list_severity:
        if "Mild" in sev:
            palette[sev] = "forestgreen"
        elif "Severe" in sev:
            palette[sev] = "firebrick"
        else:
            palette[sev] = "royalblue"

    plt.figure(figsize=(5,5))
    flierprops = dict(marker='o', markerfacecolor='None', markersize=5, markeredgecolor='black')
    ax = sns.boxplot(data=df_beta_marker, 
                     x="Severity_visit", 
                     y=marker, 
                     order=list_severity, 
                     flierprops=flierprops, 
                     palette=palette)
    pairs = [(i[1]["Group1"], i[1]["Group2"]) for i in result_test_marker.iterrows()]   
    pvalues = [i[1]["P-value"] for i in result_test_marker.iterrows()]
    add_stat_annotation(ax, data=df_beta_marker, x="Severity_visit", y=marker, order=list_severity,
                    box_pairs=pairs,
                    perform_stat_test=False, pvalues=pvalues,
                    test=None, text_format='star', loc='outside', verbose=2)
    
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_xlabel("Severity", fontsize=0)
    ax.set_ylabel("Proportion of methylated CpGs ($\\beta$-value)", fontsize=12)
    ax.set_xticks(ticks=list(range(len(list_severity))), labels=xticklabels)
    ax.plot([0.7,2.3],[-.15,-.15], color="k", transform=ax.get_xaxis_transform(), linewidth=2.5, clip_on=False)
    ax.plot([2.7,4.3],[-.15,-.15], color="k", transform=ax.get_xaxis_transform(), linewidth=2.5, clip_on=False)
    plt.text(x=0.34, y=-0.05, s="Acute Phase", weight="bold", transform=plt.gcf().transFigure)
    plt.text(x=0.63, y=-0.05, s="Recovery Phase", weight="bold", transform=plt.gcf().transFigure)
    plt.savefig(path_boxplot, dpi=600, bbox_inches="tight")
    plt.ylim(-10, 100)
    # plt.title(title)
    plt.show()
    plt.close()

# %%
list_files_methylcpgmin = get_list_files_methylcpgmin(dir_methylcpgmin)
list_methyl_sample = get_list_methyl_sample(list_files_methylcpgmin)
list_methyl_sample_filtered = list(filter(lambda x: not str(x).startswith("C19-R"), list_methyl_sample))
df_beta_all_hyper = get_dataframe_methylation_beta_samples(infilenamehyper, list_methyl_sample_filtered)
df_beta_all_hypo = get_dataframe_methylation_beta_samples(infilenamehypo, list_methyl_sample_filtered)
list_metainfo = list(df_beta_all_hyper.iloc[:, -8:].columns)
dict_markers_hyper = get_dict_deg_dmp_overlap_markers(file_deg_dmp_overlap, direction="hyper")
dict_markers_hypo = get_dict_deg_dmp_overlap_markers(file_deg_dmp_overlap, direction="hypo")
df_beta_selec_hyper = df_beta_all_hyper[list(dict_markers_hyper.keys()) + list_metainfo]
df_beta_selec_hypo = df_beta_all_hypo[list(dict_markers_hypo.keys()) + list_metainfo]

dict_order_visit = {"Control": 0, "First": 0.1, "Last": 0.2}
dict_order_group = {"Healthy": 0, "Mild": 0.1, "Severe": 0.2}
list_severity = list(df_beta_all_hyper["Severity_visit"].unique())
list_severity_sorted = sorted(list_severity, key=lambda x: (dict_order_visit[x.split("_")[-1]], dict_order_group[x.split("_")[0]]))
list_comb_sev_grp = list(combinations(list_severity_sorted, 2))
list_combination_sev_grps_filtered = list(filter(lambda x: (x[0].split("_")[1] == x[1].split("_")[1] or (x[0] == "Healthy_Control")), list_comb_sev_grp))

for marker, gensym in dict_markers_hyper.items():
    df_beta_selec_hyper_copy = df_beta_selec_hyper.copy()
    dict_beta_hyper_each_marker = df_beta_selec_hyper_copy.groupby("Severity_visit")[marker].apply(np.array).to_dict()
    dict_combinatorial_beta_hyper_each_marker = get_dict_combinatorial_beta_sev_grp_compared(dict_beta_hyper_each_marker, list_combination_sev_grps_filtered)
    outdirhyper = os.path.join(outdir,"hyper")
    os.makedirs(outdirhyper, exist_ok=True)
    res_stat_hyper = get_stat_results_marker_beta_across_groups(outdirhyper, dict_combinatorial_beta_hyper_each_marker)
    draw_boxplot_beta_each_marker(outdirhyper, marker, gensym, df_beta_selec_hyper_copy, res_stat_hyper, list_severity_sorted)

for marker, gensym in dict_markers_hypo.items():
    df_beta_selec_hypo_copy = df_beta_selec_hypo.copy()
    dict_beta_hypo_each_marker = df_beta_selec_hypo_copy.groupby("Severity_visit")[marker].apply(np.array).to_dict()
    dict_combinatorial_beta_hypo_each_marker = get_dict_combinatorial_beta_sev_grp_compared(dict_beta_hypo_each_marker, list_combination_sev_grps_filtered)
    outdirhypo = os.path.join(outdir,"hypo")
    os.makedirs(outdirhypo, exist_ok=True)
    res_stat_hypo = get_stat_results_marker_beta_across_groups(outdirhypo, dict_combinatorial_beta_hypo_each_marker)
    draw_boxplot_beta_each_marker(outdirhypo, marker, gensym, df_beta_selec_hypo_copy, res_stat_hypo, list_severity_sorted)

# %%
