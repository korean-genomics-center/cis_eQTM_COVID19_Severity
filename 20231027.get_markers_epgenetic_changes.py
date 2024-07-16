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
from scipy.stats import mannwhitneyu
from statannot import add_stat_annotation

# %%
# path_methyl = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Results/Infectomics_COVID-19_Methyl/MethylKitTable/methylkittable.by_hjryu.20231027/Methylseq_allsamples.20231027.methyl_hypo_severe_mild.tsv"

# df_beta = pd.read_csv(path_methyl, sep="\t")
# %%
direction = "hyper"
path_sev_info = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/9_clinical/Infectomics_Severity_Information_Methyl_20231106.tsv"
dir_methylcpgmin = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/MethylCpGMin"
infilename = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/first/{direction}/dictionary_marker_freqC_all_samples_20231211.pk.gz"
annot_methyl = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231211/severe_mild_firstVisit/annotation/methyl_{direction}_severe_mild_firstVisit_annotatr.tsv"
list_drop_samples = ['C19-C056-V3',
                    'C19-C053-V3',
                    'C19-C050-V3',
                    'C19-C051-V3',
                    'C19-C060-V3',
                    'C19-C055-V3',
                    'C19-C061-V3',
                    'C19-C045-V3',
                    'C19-C052-V3',
                    'C19-C047-V3']
outdir = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/first/{direction}/boxplot_20231205"
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

def get_dict_mean_beta_all_samples(dict_cpgmin_all_samples):
    dict_meanbeta_all_samples = dict()
    for sampleid, dict_cpgmin in dict_cpgmin_all_samples.items():
        list_beta = list(map(float, dict_cpgmin.values()))
        mean_beta = np.mean(list_beta)
        dict_meanbeta_all_samples[sampleid] = mean_beta

    return dict_meanbeta_all_samples

def meandiff(a, b):
    deltamean = np.mean(b) - np.mean(a)

    return deltamean

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

def draw_boxplot_per_marker(df_beta, result_t, list_severity, dict_marker_symbol, outdir):
    list_severity_cnt = list()
    for sev in list_severity:
        if sev in dict_cnt.keys():
            cnt = dict_cnt[sev]
            sev = sev.split("_")[0][0] + "." + sev.split("_")[1]
            sev_nline = '\n'.join(sev.split('_'))
            sev_cnt = f"{sev_nline}\n({cnt})"
            list_severity_cnt.append(sev_cnt)

    palette = dict()
    for sev in list_severity:
        if "Mild" in sev:
            palette[sev] = "#ECECEC"
        elif "Severe" in sev:
            palette[sev] = "#D68C78"
        else:
            palette[sev] = "#75A7C3"

    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#D68C78', edgecolor='k',label='Severe'),
                    Patch(facecolor='#ECECEC', edgecolor='k',label='Mild'),
                    Patch(facecolor='#75A7C3', edgecolor='k',label='Healthy')]

    plt.figure(figsize=(10,5))
    flierprops = dict(marker='o', markerfacecolor='None', markersize=5, markeredgecolor='black')
    # meanprops = dict(marker='s', markerfacecolor='white', markersize=7, markeredgecolor='black')
    ax = sns.boxplot(data=df_beta, x="Severity_visit", y=marker, order=list_severity, flierprops=flierprops, palette=palette)
    plt.xticks(list(range(len(list_severity_cnt))), list_severity_cnt)

    pairs = [(i[1]["Group1"], i[1]["Group2"]) for i in result_t.iterrows()]   
    pvalues = [i[1]["Padj"] for i in result_t.iterrows()]
    add_stat_annotation(ax, data=df_beta, x="Severity_visit", y=marker, order=list_severity,
                    box_pairs=pairs,
                    perform_stat_test=False, pvalues=pvalues,
                    test=None, text_format='star', loc='inside', verbose=2)

    list_xticks = list()
    for x in ax.get_xticks():
        if x%2 == 0:
            list_xticks.append(x)
    [ax.axvline(x+.5,color='k',lw=1, linestyle="dashed") for x in list_xticks] 
    
    ax.legend(handles=legend_elements, loc="upper right")
    
    plt.xlabel("Severity", fontsize=11)
    markername = f"{marker}({dict_marker_symbol[marker]})"
    plt.ylabel(markername, fontsize=11)
    plt.tight_layout()
    plt.show()
    path_boxplot = os.path.join(outdir, f"boxplot_severity_{'_'.join(marker.split(':'))}_20231205.png")
    plt.savefig(path_boxplot, dpi=300)
    plt.close()


#%%
list_files_methylcpgmin = get_list_files_methylcpgmin(dir_methylcpgmin)
list_methyl_sample = get_list_methyl_sample(list_files_methylcpgmin)
# dict_severity_group = get_dict_severity_group(path_sev_info, list_methyl_sample)
# dict_severity_group = remove_severity_group_no_content(dict_severity_group)
dict_cpgmin_all_samples = load_pickle(infilename)     
dict_cpgmin_all_samples_marker_fixed = dict()
for sampleid, dict_cpgs in dict_cpgmin_all_samples.items():
    dict_cpgs_marker_fixed = dict()
    for marker, freqC in dict_cpgs.items():
        fixed_marker = ":".join(marker)
        dict_cpgs_marker_fixed[fixed_marker] = freqC
    
    dict_cpgmin_all_samples_marker_fixed[sampleid] = dict_cpgs_marker_fixed

df_beta = pd.DataFrame(dict_cpgmin_all_samples_marker_fixed).astype(float).reset_index(drop=False).rename(columns={"index": "Posname"})
df_beta_transposed = df_beta.set_index("Posname").T
df_beta_transposed_reidx = df_beta_transposed.reset_index(drop=False).rename(columns={"index": "Sample_ID"})
df_sev_sample = pd.read_csv(path_sev_info, sep="\t") 
df_merged = pd.merge(df_beta_transposed_reidx, df_sev_sample, how="inner", on="Sample_ID")
df_merged = df_merged[~df_merged["Sample_ID"].isin(list_drop_samples)]
df_merged["Severity_group"] = df_merged["Severity_group"].fillna("Healthy")
df_merged["Visit"] = df_merged["Visit"].fillna("Control")

dict_order_visit = {"Control": 0, "First": 1, "Last": 2, "Convalescent": 3, "LongCOVID": 4}
dict_order_sev = {"Healthy": 0, "Mild": 0.1, "Severe": 0.2}

dict_cnt = dict()
dict_grpby_sev_visit = df_merged.groupby(["Severity_group", "Visit"]).apply(len).to_dict()
for key, cnt in dict_grpby_sev_visit.items():
    sev_visit = "_".join(key)
    dict_cnt[sev_visit] = cnt

list_severity = sorted(list(df_merged["Severity_visit"].unique()), key=lambda x:dict_order_visit[x.split("_")[-1]]+dict_order_sev[x.split("_")[0]])

list_severity_cnt = list()
for sev in list_severity:
    if sev in dict_cnt.keys():
        cnt = dict_cnt[sev]
        sev = sev.split("_")[0][0] + "." + sev.split("_")[1]
        sev_nline = '\n'.join(sev.split('_'))
        sev_cnt = f"{sev_nline}\n({cnt})"
        list_severity_cnt.append(sev_cnt)

list_comb_sev_visit = list(combinations(list_severity, 2))
list_comb_sev_visit_rmv_diff_visit = list(filter(lambda x: (x[0].split("_")[1] == x[1].split("_")[1]) or (x[0] == "Healthy_Control" and "First" in x[1]), list_comb_sev_visit))

dict_marker_symbol = annotate_gene_name_marker(annot_methyl)

list_epi_change_markers_sev = list()
list_markers = list(filter(lambda x: str(x).startswith("chr"), list(df_merged.columns)))
for marker in list_markers:
    df_merged_copy = df_merged.copy()
    dict_sev_visit = df_merged_copy.groupby("Severity_visit")[marker].apply(np.array).to_dict()

    from collections import defaultdict
    dict_marker_comb_visits = defaultdict(list)
    for comb_visits in list_comb_sev_visit_rmv_diff_visit:
        for visit in comb_visits:
            meanbeta = dict_sev_visit[visit]
            dict_marker_comb_visits[comb_visits].append(meanbeta)

    from scipy.stats import ttest_ind
    from statsmodels.stats.multitest import fdrcorrection

    dict_visit_stats = dict()
    for visits, meanbetas in dict_marker_comb_visits.items():
        deltamean = meandiff(*meanbetas)
        tstat, pval = ttest_ind(*meanbetas, equal_var=False)
        dict_visit_stats[visits] = {"Delta": deltamean, "Statistics": tstat, "P-value": pval}
    result_t = pd.DataFrame.from_dict(dict_visit_stats, orient="index").reset_index(drop=False)
    result_t = result_t.rename(columns={"level_0": "Group1", "level_1": "Group2"})
    list_sigtrue = fdrcorrection(result_t["P-value"])[0]
    list_padj = fdrcorrection(result_t["P-value"])[1]
    result_t["Padj"] = list_padj
    result_t["TestSignificant"] = list_sigtrue
    result_t["Marker"] = marker
    result_t["Gene"] = dict_marker_symbol[marker]
    list_epi_change_markers_sev.append(result_t)

df_markers_diff_beta = pd.concat(list_epi_change_markers_sev, axis=0)
df_markers_diff_beta_sigonly = df_markers_diff_beta[df_markers_diff_beta["TestSignificant"] == True]
df_markers_diff_beta_sigonly_conval = df_markers_diff_beta_sigonly[df_markers_diff_beta_sigonly["Group1"].str.contains("Convalescent")]
df_markers_diff_beta_sigonly_long = df_markers_diff_beta_sigonly[df_markers_diff_beta_sigonly["Group1"].str.contains("Mild_LongCOVID")]
# list_marker_epi_change_after = list(df_markers_diff_beta_sigonly_conval["Marker"].unique()) + list(df_markers_diff_beta_sigonly_long["Marker"].unique())
list_marker_sig_diff = list(df_markers_diff_beta_sigonly["Marker"].unique())

for marker in list_marker_sig_diff:
    df_merged_copy = df_merged.copy()
    dict_sev_visit = df_merged_copy.groupby("Severity_visit")[marker].apply(np.array).to_dict()

    dict_marker_comb_visits = defaultdict(list)
    for comb_visits in list_comb_sev_visit_rmv_diff_visit:
        for visit in comb_visits:
            meanbeta = dict_sev_visit[visit]
            dict_marker_comb_visits[comb_visits].append(meanbeta)

    dict_visit_stats = dict()
    for visits, meanbetas in dict_marker_comb_visits.items():
        deltamean = meandiff(*meanbetas)
        tstat, pval = ttest_ind(*meanbetas, equal_var=False)
        dict_visit_stats[visits] = {"Delta": deltamean, "Statistics": tstat, "P-value": pval}
    result_t = pd.DataFrame.from_dict(dict_visit_stats, orient="index").reset_index(drop=False)
    result_t = result_t.rename(columns={"level_0": "Group1", "level_1": "Group2"})
    list_sigtrue = fdrcorrection(result_t["P-value"])[0]
    list_padj = fdrcorrection(result_t["P-value"])[1]
    result_t["Padj"] = list_padj
    result_t["TestSignificant"] = list_sigtrue
    result_t["Marker"] = marker
    result_t["Gene"] = dict_marker_symbol[marker]
    
    draw_boxplot_per_marker(df_merged, result_t, list_severity, dict_marker_symbol, outdir)

# #%%
# df_sev_sample = pd.read_csv(path_sev_info, sep="\t") 
# list_colheader = list(df_beta.columns)
# list_samples = list(filter(lambda x: x.startswith("C19"), list_colheader))
# df_sev_sample = df_sev_sample[df_sev_sample["Sample_ID"].isin(list_samples)]

# dict_sev_sample = dict()
# for sev, id in zip(df_sev_sample["Severity_visit"], df_sev_sample["Sample_ID"]):
#     if dict_sev_sample.get(sev, None) == None:
#         dict_sev_sample[sev] = list()
#     dict_sev_sample[sev].append(id)

# # list_keys = list(dict_sev_sample.keys())
# # list_second_keys = list(filter(lambda x: "Second" in x, list_keys))
# # for key in list_second_keys:
# #     dict_sev_sample.pop(key)

# dict_cpgs = dict()
# for _, row in df_beta.iterrows():
#     dict_row = row.to_dict()
#     posname = dict_row["Posname"]
#     dict_cpgs[posname] = dict()
#     for key, samples in dict_sev_sample.items():
#         dict_cpgs[posname][key] = list(map(dict_row.__getitem__, samples))

# first_pos = list(dict_cpgs.keys())[0]
# list_comb = list(combinations(dict_cpgs[first_pos].keys(), 2))
# list_comb = list(filter(lambda x: x[0].split("_")[-1] == x[1].split("_")[-1], list_comb))
# list_visitpair = list(map(lambda x : ':'.join(sorted(x)), list_comb))
# table_comb = pd.DataFrame(columns = ["Posname"] + list_visitpair)

# from scipy.stats import ranksums

# for ind, row in df_beta.iterrows():
#     dict_row = row.to_dict()
#     dict_tmp = dict()
#     for key, samples in dict_sev_sample.items():
#         dict_tmp[key] = list(map(dict_row.__getitem__, samples))
#     posname = row["Posname"]

#     table_comb.loc[ind, "Posname"] = posname
#     for g1, g2 in list_comb:
#         vals1 = dict_tmp[g1]
#         vals2 = dict_tmp[g2]
#         median1 = np.median(vals1)
#         median2 = np.median(vals2)
#         # deltamedian = median1 - median2
#         combname = ':'.join(sorted([g1, g2]))
#         _, pval = ranksums(vals1, vals2)
#         table_comb.loc[ind, combname] = pval

# from statsmodels.stats.multitest import fdrcorrection

# table_comb_before_correction = table_comb.copy()
# table_comb_after_correction = table_comb.copy()
# for pair in list_visitpair:
#     pval = table_comb_before_correction[pair]
#     _, padj = fdrcorrection(pval)
#     table_comb_after_correction[pair] = padj

# table_comb_sig = table_comb_after_correction[np.logical_and(table_comb["Mild_First:Severe_First"] < 0.05, table_comb["Mild_Last:Severe_Last"]< 0.05)]
# table_comb_sig = table_comb_sig.set_index("Posname")
# table_comb_sig_logged = table_comb_sig.applymap(lambda x: 1-float(x))
# table_comb_sig_logged["order"] = 1-(table_comb_sig_logged["Mild_First:Severe_First"].astype(float) * table_comb_sig_logged["Mild_Last:Severe_Last"].astype(float))
# table_comb_sig_sorted = table_comb_sig_logged.sort_values(by=["order"], ascending=True)
# table_comb_sig_sorted_reidx = table_comb_sig_sorted.reset_index(drop=False)

# # table_comb_sig_logged = table_comb_sig.applymap(lambda x: -math.log10(float(x)))
# # table_comb_sig_logged["order"] = table_comb_sig_logged["Mild_First:Severe_First"].astype(float) + table_comb_sig_logged["Mild_Last:Severe_Last"].astype(float)
# # table_comb_sig_sorted = table_comb_sig_logged.sort_values(by=["order"], ascending=False)


# # %%
# annot_methyl = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231120_withoutOV/annotation/methyl_hypo_severe_mild_annotatr.tsv"
# df_annot_methyl = pd.read_csv(annot_methyl, sep="\t")
# df_annot_methyl["marker"] = df_annot_methyl["seqnames"].astype(str) + ":" + df_annot_methyl["start"].apply(lambda x: int(x)-1).astype(str)

# dict_marker_symbol = dict()
# for marker, symbol in zip(df_annot_methyl["marker"], df_annot_methyl["annot.symbol"]):
#     if dict_marker_symbol.get(marker, None) == None:
#         dict_marker_symbol[marker] = list()
#     dict_marker_symbol[marker].append(symbol)

# for key, list_val in dict_marker_symbol.items():
#     list_unique_val = list(set(list_val))
#     list_unique_val.remove(np.nan)
#     if len(list_unique_val) > 1:
#         unique_val = list_unique_val[-1]
#     else:
#         unique_val = "-".join(list_unique_val)
#     if unique_val == "":
#         unique_val = "None"
#     dict_marker_symbol[key] = unique_val

# table_comb_sig_sorted_reidx["GeneSymbol"] = table_comb_sig_sorted_reidx["Posname"].apply(lambda x: dict_marker_symbol.get(x, "None"))
# table_comb_sig_sorted_reidx = table_comb_sig_sorted_reidx.drop(columns=["Posname", "order"])
# table_comb_sig_sorted_reidx = table_comb_sig_sorted_reidx.set_index("GeneSymbol", drop=True)

# # %%
# import math

# import matplotlib.pyplot as plt
# import seaborn as sns

# table_comb_sig_sorted_reidx_top = table_comb_sig_sorted_reidx.head(50)
# plt.figure(figsize=(10, 10))
# ax = sns.heatmap(data=table_comb_sig_sorted_reidx_top, xticklabels=True, yticklabels=True, cmap="Greens_r")
# ax.set_xticklabels(table_comb_sig_sorted_reidx_top.columns, fontsize=13)
# ax.set_yticklabels(table_comb_sig_sorted_reidx_top.index, fontsize=13)
# plt.ylabel("Delayed Methylation Markers", fontsize=15)
# plt.xlabel("Visit Pairs", fontsize=15)
# plt.show()
# %%
