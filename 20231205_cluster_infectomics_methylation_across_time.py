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
outdir = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/first/{direction}/boxplot_20231211"
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
df_beta_sev["Visit_order"] = df_beta_sev["Visit_order"].fillna("Visit0")
df_beta_sev["Severity_group"] = df_beta_sev["Severity_group"].fillna("Healthy")
df_beta_sev["Subject_ID"] = df_beta_sev["Subject_ID"].fillna(df_beta_sev["Sample_ID"])
df_beta_sev = df_beta_sev[df_beta_sev["Visit_order"]!="Visit6"]
df_beta_set_idx = df_beta_sev.set_index("Sample_ID")
series_grpby_beta_sev = df_beta_set_idx.groupby(["Severity_group", "Subject_ID", "Visit_order"]).apply(np.array)
df_grpby_beta_sev = series_grpby_beta_sev.reset_index(drop=False).rename(columns={0:"beta"})
df_grpby_beta_sev["beta"] = df_grpby_beta_sev["beta"].apply(lambda x: x[0][0:-8])

list_markers = list(filter(lambda x: str(x).startswith("chr"), list(df_beta_set_idx.columns)))

# %%
from collections import defaultdict

dict_subj_beta = defaultdict(list)
prev = list()
for ind, row in df_grpby_beta_sev.iterrows():
    dict_row = row.to_dict()
    subjid = dict_row["Subject_ID"]
    beta = dict_row["beta"]

    if subjid not in prev:
        dict_subj_beta[subjid].append(beta)
        prev.append(subjid)
    else:
        dict_subj_beta[subjid].append(beta)
    
    prev = list()

dict_subj_diff_arr = dict()
list_subj_id_no_paired_samples = list()
for key, list_value in dict_subj_beta.items():
    num_vals = len(list_value)
    if num_vals > 1:
        arr1 = np.array(list_value[0], dtype=np.float64)
        arr2 = np.array(list_value[1], dtype=np.float64)
        diff_arr = (arr2 - arr1)
        dict_subj_diff_arr[key] = diff_arr
    else:
        print(f"{key} only has only single sample")
        list_subj_id_no_paired_samples.append(key)

df_subj_sev = df_sev[["Subject_ID", "Severity_group"]].drop_duplicates(keep="first")
df_subj_diff = pd.DataFrame.from_dict(dict_subj_diff_arr, orient="index", columns=list_markers).reset_index(drop=False).rename(columns={"index": "Subject_ID"})
df_subj_diff_sev = pd.merge(df_subj_diff, df_subj_sev, how="inner", on="Subject_ID")

alist = list()
for marker in list_markers:
    df_subj_diff_sev_copy = df_subj_diff_sev.copy()
    subj_diff_sev = df_subj_diff_sev_copy.groupby("Severity_group")[marker].apply(np.mean)
    arr = subj_diff_sev.values
    alist.append(arr)
mat = np.vstack(alist)

df_sev_marker = pd.DataFrame(mat, index=list_markers, columns=["Mild", "Severe"])

# %%
# num_thres = 10
# dict_clust = dict()
# for ind, row in df_sev_marker.iterrows():
#     dict_row = row.to_dict()
#     mild_val = dict_row["Mild"]
#     sev_val = dict_row["Severe"]
#     if (mild_val == num_thres) or (sev_val == num_thres):
#         dict_clust[ind] = "center"
#     elif (mild_val > num_thres) and (sev_val > num_thres):
#         dict_clust[ind] = "top_right"
#     elif (mild_val > num_thres) and (sev_val < -num_thres):
#         dict_clust[ind] = "bottom_right"
#     elif (mild_val < -num_thres) and (sev_val > num_thres):
#         dict_clust[ind] = "top_left"
#     elif (mild_val < -num_thres) and (sev_val < -num_thres):
#         dict_clust[ind] = "bottom_left"
#     else:
#         print(f"{ind} not included")

num_thres = 10
center_thres = 5
dict_clust = dict()
for ind, row in df_sev_marker.iterrows():
    dict_row = row.to_dict()
    mild_val = dict_row["Mild"]
    sev_val = dict_row["Severe"]
    if (-num_thres < mild_val < num_thres) and (-num_thres < sev_val < num_thres):
        dict_clust[ind] = "center"
    elif (mild_val > num_thres) and (sev_val > num_thres):
        dict_clust[ind] = "top_right"
    elif (mild_val > num_thres) and (sev_val < -num_thres):
        dict_clust[ind] = "bottom_right"
    elif (mild_val < -num_thres) and (sev_val > num_thres):
        dict_clust[ind] = "top_left"
    elif (mild_val < -num_thres) and (sev_val < -num_thres):
        dict_clust[ind] = "bottom_left"
    elif (mild_val < -num_thres)  and (-num_thres < sev_val < num_thres):
        dict_clust[ind] = "middle_left"
    elif (mild_val > num_thres)  and (-num_thres < sev_val < num_thres):
        dict_clust[ind] = "middle_right"
    elif (-num_thres < mild_val < num_thres)  and (sev_val < -num_thres):
        dict_clust[ind] = "middle_bottom"
    elif (-num_thres < mild_val < num_thres)  and (sev_val > num_thres):
        dict_clust[ind] = "middle_top"
    else:
        print(f"{ind} not included")

    if (-center_thres < mild_val < center_thres) and (-center_thres < sev_val < center_thres):
        dict_clust[ind] = "center"

df_sev_marker_reindx = df_sev_marker.reset_index(drop=False).rename(columns={"index": "Marker"})
df_sev_marker_reindx["Cluster"] = df_sev_marker_reindx["Marker"].apply(lambda x: dict_clust.get(x, "None"))
df_sev_marker_reindx.to_csv("/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/marker_cluster.txt", sep="\t", index=True)

# %%
# import numpy as np
# from sklearn.cluster import AgglomerativeClustering
# dict_marker_cluster = dict()
# cluster = AgglomerativeClustering(n_clusters=5, affinity='euclidean', linkage='average')
# list_clust = cluster.fit_predict(mat_sign)

# %%
plt.figure(figsize=(5, 5))
sns.scatterplot(data=df_sev_marker_reindx, x="Mild", y="Severe", hue="Cluster", palette="Spectral", zorder=3)
# plt.axhline(y=0, linestyle="dashed", color="k", zorder=5)
# plt.axvline(x=0, linestyle="dashed", color="k", zorder=5)
# plt.hlines(y=0, xmin=-20, xmax=20, linestyles="dashed", colors="k", zorder=5)
# plt.vlines(x=0, ymin=-20, ymax=20, linestyles="dashed", colors="k", zorder=5)

# 배수로 반올림/내림
import math
xlim_left = min(-15, math.floor(df_sev_marker_reindx["Mild"].min()/5)*5)
xlim_right = max(15, math.ceil(df_sev_marker_reindx["Mild"].max()/5)*5)
ylim_bottom = min(-15, math.floor(df_sev_marker_reindx["Severe"].min()/5)*5)
ylim_top = max(15, math.ceil(df_sev_marker_reindx["Severe"].max()/5)*5)
plt.xticks(list(range(xlim_left, xlim_right+1, 5)))
plt.yticks(list(range(ylim_bottom, ylim_top+1, 5)))
plt.xlim(xlim_left, xlim_right)
plt.ylim(ylim_bottom, ylim_top)
plt.legend(loc="upper right", bbox_to_anchor=(1.5, 1.0))
plt.xlabel("Mild (Last Beta - First Beta)")
plt.ylabel("Severe (Last Beta - First Beta)")
plt.grid(visible=True,axis='both', linestyle = "--",zorder=0)
plt.show()
plt.close()

# %%
df_sev_visit_beta = pd.DataFrame(df_grpby_beta_sev["beta"].to_list(), index=series_grpby_beta_sev.index, columns=list_markers).reset_index(drop=False)
df_sev_visit_beta_melted = pd.melt(df_sev_visit_beta, ["Severity_group", "Subject_ID", "Visit_order"]).rename(columns={"variable": "Marker", "value": "Beta"})
dict_rename_visit = {"Visit0": "Healthy", "Visit1": "First", "Visit3": "Last", "Visit4": "Last", "Visit5": "Convalescent"}
df_sev_visit_beta_melted["Visit_group"] = df_sev_visit_beta_melted["Visit_order"].apply(lambda x: dict_rename_visit[x])

# %% [remove unpaired subjects from the beta matrix]
df_sev_visit_beta_melted_filtered_unpaired = df_sev_visit_beta_melted[~np.logical_and(df_sev_visit_beta_melted["Visit_order"] == "Visit1", df_sev_visit_beta_melted["Subject_ID"].isin(list_subj_id_no_paired_samples))]

# %%
def get_palette(df):
    palette = dict()
    for sev in df["Severity_group"].unique():
        if "Mild" in sev:
            palette[sev] = "#aaaaaa"
        elif "Severe" in sev:
            palette[sev] = "#D68C78"
        else:
            palette[sev] = "#75A7C3"
            
    return palette

# %%
def center_methylation_value_into_healthy_controls(df_clustered_beta_visit, col_sev, col_marker, col_beta, hc_name = "Healthy"):
    col_centered_beta = f"{col_beta}_Centered"
    df_centered = pd.DataFrame(columns = list(df_clustered_beta_visit)+[col_centered_beta])
    for marker in df_clustered_beta_visit[col_marker].unique():
        df_clustered_beta_visit_single_marker = df_clustered_beta_visit[df_clustered_beta_visit[col_marker] == marker].copy()
        list_beta_hc = df_clustered_beta_visit_single_marker[df_clustered_beta_visit_single_marker[col_sev] == hc_name][col_beta].to_list()
        mean_beta_hc = np.mean(list_beta_hc)
        df_clustered_beta_visit_single_marker[col_centered_beta] = df_clustered_beta_visit_single_marker[col_beta] - mean_beta_hc
        df_centered = pd.concat([df_centered, df_clustered_beta_visit_single_marker])
    
    return df_centered

# %%
dict_top10_cluster = dict()

series_cluster_markers = df_sev_marker_reindx.groupby("Cluster")["Marker"].apply(list)
cluster_names = list(series_cluster_markers.index)
for cluster_name, cluster_markers in zip(cluster_names, series_cluster_markers):
    df_sev_visit_beta_melted_copy = df_sev_visit_beta_melted_filtered_unpaired.copy()
    df_clustered_beta_visit = df_sev_visit_beta_melted_copy[df_sev_visit_beta_melted_copy["Marker"].isin(cluster_markers)]

    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#D68C78', edgecolor='k',label='Severe'), Patch(facecolor='#aaaaaa', edgecolor='k',label='Mild'), Patch(facecolor='#75A7C3', edgecolor='k',label='Healthy')]
    
    df_clustered_beta_visit_centered = center_methylation_value_into_healthy_controls(df_clustered_beta_visit, "Severity_group", "Marker", "Beta")
    palette = get_palette(df_clustered_beta_visit_centered)

    df_clustered_beta_visit_centered_grpby_beta = (df_clustered_beta_visit_centered.groupby(["Marker", "Severity_group", "Visit_group"])["Beta_Centered"].apply(np.mean).reset_index(drop=False))
    df_clustered_beta_visit_centered_grpby_beta_sev_first = df_clustered_beta_visit_centered_grpby_beta[np.logical_and(df_clustered_beta_visit_centered_grpby_beta["Visit_group"]=="First", df_clustered_beta_visit_centered_grpby_beta["Severity_group"]=="Severe")]

    df_clustered_beta_visit_centered_grpby_beta_sev_first_sorted = df_clustered_beta_visit_centered_grpby_beta_sev_first.sort_values(by=["Beta_Centered"], ascending=False)

    dict_top10_cluster[cluster_name] = df_clustered_beta_visit_centered_grpby_beta_sev_first_sorted.head(10)["Marker"].to_list()

# %%
    # fig, ax = plt.subplots(1,1)

    # ax.set_title(cluster_name)
    # for marker_info in df_clustered_beta_visit_centered.groupby("Marker"):
    #     marker_name = marker_info[0]
    #     marker_df = marker_info[1]
    #     pointplot_ax = sns.pointplot(data=marker_df, x="Visit_group", y="Beta_Centered", hue="Severity_group", palette=palette, capsize=.15, errwidth=0.5, legend=False, ax=ax, zorder = 3)
    #     plt.setp(pointplot_ax.collections, alpha=.4, zorder = 3)
    #     plt.setp(pointplot_ax.lines, alpha=.4, zorder = 3)
    #     ax.get_legend().remove()

    # darker_palette = {
    #     "Healthy" : "#347091",
    #     "Mild" : "#3d3d3d",
    #     "Severe" : "#ad482d"
    # }
    # pointplot_ax = sns.pointplot(data=df_clustered_beta_visit_centered, x="Visit_group", y="Beta_Centered", hue="Severity_group", palette=darker_palette, capsize=.15, errwidth=0.5, legend=False, zorder=999, ax = ax, edgecolor = 'k', edgewidth = 5)
    # plt.setp(pointplot_ax.collections, zorder = 999)
    # plt.setp(pointplot_ax.lines, zorder = 999)

    # ax.grid(visible=True, axis="y", linestyle="--", linewidth=1, zorder = 200)
    # ax.legend(handles=legend_elements, loc="upper right")
    # plt.show()
    # plt.close()

# %%
direction = "hypo"
annot_methyl = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231211/severe_mild_firstVisit/annotation/methyl_{direction}_severe_mild_firstVisit_annotatr.tsv"
dict_marker_symbol = annotate_gene_name_marker(annot_methyl)
list_marker = df_sev_marker_reindx["Marker"].to_list()
list_symbol = list(map(lambda x: dict_marker_symbol[x], list_marker))
df_sev_marker_reindx["Symbol"] = list_symbol
series_cluster_symbol = df_sev_marker_reindx.groupby("Cluster")["Symbol"].apply(list)
list_cluster = list(series_cluster_symbol.index)
for cluster, list_gene in zip(list_cluster, series_cluster_symbol):
    print(cluster)
    print("----")
    for gene in list(set(list_gene)):
        print(gene)
# %%
direction = "hypo"
annot_methyl = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231211/severe_mild_firstVisit/annotation/methyl_{direction}_severe_mild_firstVisit_annotatr.tsv"
dict_marker_symbol = annotate_gene_name_marker(annot_methyl)
list_marker = df_sev_marker_reindx["Marker"].to_list()
for category, list_markers in dict_top10_cluster.items():
    list_symbols = list(map(lambda x: dict_marker_symbol[x], list_markers))
    print(category, list_symbols)
# %%
