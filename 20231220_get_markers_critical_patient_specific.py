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
from scipy.stats import f_oneway, kruskal, mannwhitneyu, ttest_ind
from statannot import add_stat_annotation
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import fdrcorrection

# %%
direction = "hypo"
path_sev_info = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/9_clinical/Infectomics_Severity_Information_Methyl_20231106.tsv"
dir_methylcpgmin = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/MethylCpGMin"
infilename = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/first/{direction}/dictionary_marker_freqC_all_samples_20231211.pk.gz"
annot_methyl = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231211/severe_mild_firstVisit/methyl_{direction}_severe_mild_commonMarkers_LOO_markerlist.tsv"
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
outdir = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/first/{direction}/boxplot_20231221"
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

# %%
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

list_markers = list(filter(lambda x: str(x).startswith("chr"), list(df_beta_set_idx.columns)))

import numpy as np
from joblib import Parallel, delayed
# %%
from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd


def get_critical_tukey_results_between_sev_mild_visit(df_beta_set_idx, marker, visit="First"):
    df_beta_set_idx_copy = df_beta_set_idx.copy()
    df_beta_set_idx_copy["Severity_Visit"] = df_beta_set_idx_copy["Severity"].astype(str) + "_" + df_beta_set_idx_copy["Visit"] 
    beta_marker_by_sev = df_beta_set_idx_copy.groupby(["Severity_Visit"])[marker].apply(np.array)
    beta_marker_by_sev_reset_idx = beta_marker_by_sev.reset_index(drop=False)
    list_severity_visit = beta_marker_by_sev_reset_idx["Severity_Visit"]
    list_select_bool = list(map(lambda x: x.split("_")[-1] in ["Healthy", visit], list_severity_visit))
    list_severity_visit_select = [x for x, flagx in zip(list_severity_visit, list_select_bool) if flagx]
    list_beta_vals = beta_marker_by_sev_reset_idx[marker].to_list()
    list_beta_vals_select = [y for y, flagy in zip(list_beta_vals, list_select_bool) if flagy]
    
    stat, pval = f_oneway(*list_beta_vals_select)
    critical_tukey = None
    critical_marker = None
    if pval < 0.05:
        df_beta_set_idx_copy_first = df_beta_set_idx_copy[df_beta_set_idx_copy["Visit"].isin([visit, "Healthy"])]
        posthoc = pairwise_tukeyhsd(df_beta_set_idx_copy_first[marker], df_beta_set_idx_copy_first["Severity_Visit"], alpha=0.05)
        result_tukey = pd.DataFrame(data=posthoc._results_table.data[1:], columns=posthoc._results_table.data[0])
        bool_select = np.logical_and(result_tukey["group1"] == f"3_{visit}", result_tukey["group2"] == f"4_{visit}")
        result_tukey_select = result_tukey[bool_select]

        if result_tukey_select["reject"].values[0] == True:
            result_tukey["marker"] = marker
            critical_marker = marker
            critical_tukey = result_tukey
    
    return (critical_marker, critical_tukey)

# %%
with Parallel(n_jobs=150) as parallel:
    list_results = parallel(delayed(get_critical_tukey_results_between_sev_mild_visit)(df_beta_set_idx, marker, visit="Last") for marker in list_markers)
list_results_critical = list(filter(lambda x : type(x[0]) != type(None), list_results))

list_critical_markers = list(map(lambda x : x[0], list_results_critical))
list_critical_tukey = list(map(lambda x : x[1], list_results_critical))

result_tukey_all = pd.concat(list_critical_tukey, axis=0)
result_tukey_all.to_csv(os.path.join(outdir, "tukey_hsd_results_critcal_markers_20231221.tsv"), sep="\t", index=False)

# %%
dict_order_visit = {"Healthy": 0, "First": 0.1, "Last": 0.2, "Convalescent": 0.3}

dict_cnt = dict()
dict_grpby_sev_visit = df_beta_set_idx.groupby(["Severity", "Visit"]).apply(len).to_dict()
for key, cnt in dict_grpby_sev_visit.items():
    sev_visit = str(key[0]) + "_" + str(key[1])
    dict_cnt[sev_visit] = cnt

list_severity = sorted(list(df_beta_set_idx["Severity"].unique()))
list_visit = sorted(list(df_beta_set_idx["Visit"].unique()), key=lambda x: dict_order_visit[x])

list_severity_visit = ["0_Healthy"]
for sev in list_severity[1:]:
    for visit in list_visit[1:]:
        sev_visit = f"{sev}_{visit}"
        list_severity_visit.append(sev_visit)

list_severity_cnt = list()
for sev in list_severity_visit:
    if sev in dict_cnt.keys():
        cnt = dict_cnt[sev]
        sev = sev.split("_")[0][0] + "." + sev.split("_")[1]
        sev_nline = '\n'.join(sev.split('_'))
        sev_cnt = f"{sev_nline}\n({cnt})"
        list_severity_cnt.append(sev_cnt)

palette = dict()
for sev_visit in list_severity_visit:
    sev = int(sev_visit.split("_")[0])
    if 0 == sev:
        palette[sev_visit] = "#75A7C3"
    elif sev <= 2:
        palette[sev_visit] = "#aaaaaa"
    elif sev == 3:
        palette[sev_visit] = "#D68C78"
    else:
        palette[sev_visit] = '#9e3737'

from matplotlib.patches import Patch

legend_elements = [Patch(facecolor='#9e3737', edgecolor='k',label='Critical'),
                Patch(facecolor='#D68C78', edgecolor='k',label='Severe'),
                Patch(facecolor='#aaaaaa', edgecolor='k',label='Mild'),
                Patch(facecolor='#75A7C3', edgecolor='k',label='Healthy')]

for critical_marker, result_tukey in zip(list_critical_markers , list_critical_tukey):
    df_beta_set_idx_copy = df_beta_set_idx.copy()
    df_beta_set_idx_copy_marker = df_beta_set_idx_copy[[critical_marker, "Severity", "Visit"]]
    df_beta_set_idx_copy_marker["Severity_Visit"] = df_beta_set_idx_copy_marker["Severity"].astype(str) + "_" + df_beta_set_idx_copy_marker["Visit"] 
    plt.figure(figsize=(20,5))
    flierprops = dict(marker='o', markerfacecolor='None', markersize=5, markeredgecolor='black')
    ax = sns.boxplot(data=df_beta_set_idx_copy_marker, x="Severity_Visit", y=critical_marker, order=list_severity_visit, flierprops=flierprops, palette=palette)
    plt.xticks(list(range(len(list_severity_cnt))), list_severity_cnt)
    pairs = [(i[1]["group1"], i[1]["group2"]) for i in result_tukey.iterrows()]
    pvalues = [i[1]["p-adj"] for i in result_tukey.iterrows()]
    add_stat_annotation(ax, data=df_beta_set_idx_copy_marker, x="Severity_Visit", y=critical_marker, order=list_severity_visit,
                    box_pairs=pairs,
                    perform_stat_test=False, pvalues=pvalues,
                    test=None, text_format='star', loc='inside', verbose=2)

    list_xticks = list()
    for x in ax.get_xticks():
        if x%3 == 0:
            list_xticks.append(x)
    [ax.axvline(x+.5,color='k',lw=1, linestyle="dashed") for x in list_xticks] 

    # if direction == "hyper":
    #     plt.ylim(40, 100)
    # else:
    #     plt.ylim(0, 100)
    plt.legend()
    ax.set_xlabel("Severity", fontsize=11)
    ax.set_ylabel("MeanBeta", fontsize=11)

    if direction == "hyper":
        ax.legend(handles=legend_elements, loc="upper right")
    else:
        ax.legend(handles=legend_elements, loc="lower right")


    plt.tight_layout()
    # plt.savefig(path_boxplot, dpi=600)
    plt.show()
    plt.close()

# %%
direction = "hypo"
annot_methyl = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231211/severe_mild_firstVisit/annotation/methyl_{direction}_severe_mild_firstVisit_annotatr.tsv"
dict_marker_symbol = annotate_gene_name_marker(annot_methyl)
