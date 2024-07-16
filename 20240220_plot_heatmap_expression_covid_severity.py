# %%
import glob
import gzip
import os
import pickle
from collections import defaultdict
from itertools import chain

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import math
from sklearn.preprocessing import StandardScaler


# %%
def main(pathsevinfo, fileexp, filetargetgene, outfig, list_drop_samples):
    list_rna_sample = get_list_rna_sample(fileexp)
    dict_sample_severity = get_dict_sample_severity(pathsevinfo, list_rna_sample)
    dict_sample_severity_sorted = sort_dict_sample_severity_visit(dict_sample_severity)
    df_dmp_deg_sorted = get_dmp_deg_stats(filetargetgene, colmethdiff="dmp_median_methdiff", colfc="deg_median_log2FC")
    dict_marker_gene = get_dict_target_marker_gene(df_dmp_deg_sorted)
    dict_marker_direction = get_dict_target_marker_direction(df_dmp_deg_sorted)
    dict_marker_hyper = {k : v for k, v in dict_marker_direction.items() if v == "hyper"}
    list_marker_hyper = [v for k, v in dict_marker_gene.items() if k in dict_marker_hyper]
    list_unique_marker_hyper = list(set(list(chain(*list_marker_hyper))))
    dict_marker_hypo = {k : v for k, v in dict_marker_direction.items() if v == "hypo"}
    list_marker_hypo = [v for k, v in dict_marker_gene.items() if k in dict_marker_hypo]
    list_unique_marker_hypo = list(set(list((chain(*list_marker_hypo)))))
    df_marker_sample_hyper = pd.read_csv(fileexp, sep="\t").loc[list_unique_marker_hyper]
    df_marker_sample_hypo = pd.read_csv(fileexp, sep="\t").loc[list_unique_marker_hypo]
    df_marker_sample_all = pd.concat([df_marker_sample_hypo, df_marker_sample_hyper], axis=0).reset_index(drop=False).rename(columns={"index": "ID"})

    df_marker_sample_hyper["Cluster"] = get_list_cluster_markers(df_marker_sample_hyper, n_clust=1)
    df_marker_sample_hyper_cluster = df_marker_sample_hyper.sort_values(by=["Cluster"], ascending=True)
    df_marker_sample_hyper_cluster = df_marker_sample_hyper_cluster.drop(columns=["Cluster"])
    df_marker_sample_hypo["Cluster"] = get_list_cluster_markers(df_marker_sample_hypo, n_clust=3)
    df_marker_sample_hypo_cluster = df_marker_sample_hypo.sort_values(by=["Cluster"], ascending=True)
    df_marker_sample_hypo_cluster = df_marker_sample_hypo_cluster.drop(columns=["Cluster"])
    df_marker_sample_gene_cluster = pd.concat([df_marker_sample_hypo_cluster, df_marker_sample_hyper_cluster], axis=0)
    list_sample_sort = dict_sample_severity_sorted.keys()
    df_marker_sample_gene_cluster_sorted = df_marker_sample_gene_cluster[list_sample_sort]
    df_marker_sample_gene_cluster_sorted = df_marker_sample_gene_cluster_sorted.applymap(lambda x: math.log2(x+1))
    df_marker_sample_gene_cluster_sorted_reset_idx = df_marker_sample_gene_cluster_sorted.reset_index(drop=False)
    df_marker_sample_gene_cluster_sorted_reset_idx.iloc[9, 0] += "_1"
    df_marker_sample_gene_cluster_sorted_reset_idx.iloc[-1, 0] += "_2"
    df_marker_sample_gene_cluster_sorted_plot = df_marker_sample_gene_cluster_sorted_reset_idx.set_index("index")
    df_marker_sample_gene_cluster_sorted.to_csv("/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/DMPDEG/first/gene_expression_matrix.tsv", sep="\t", index=False)

    df_meta = pd.read_csv(pathsevinfo, sep="\t").set_index("Sample_ID")
    df_meta_filtered = df_meta.loc[list(df_marker_sample_all.columns)[1:]].reset_index(drop=False)
    list_sampleid = list(df_meta_filtered.groupby("Sample_ID")["Severity_group"].apply(list).index)
    list_sev = df_meta_filtered.groupby("Sample_ID")["Severity_group"].apply(list).str[0].to_list()
    
    dict_color_sev = {"Mild": "forestgreen", "Severe": "firebrick"}
    dict_color_vis = {"First": "pink", "Last": "powderblue"}
    dict_color_direction = {"hypo": "grey", "hyper": "black"}
    list_sev_color = list(map(dict_color_sev.__getitem__, list_sev))
    # list_vis_color = list(map(dict_color_vis.__getitem__, list_vis))
    list_directions = ["hypo"]*17 + ["hyper"]*2
    list_direction_color = list(map(dict_color_direction.__getitem__, list_directions))
    # color_df = pd.DataFrame({"Infection Phase": list_vis_color, "Severity Group": lis·t_sev_color}, index=list_sampleid)
    color_df = pd.DataFrame({"Severity Group": list_sev_color}, index=list_sampleid)
    row_df = pd.DataFrame({"Direction": list_direction_color}, index=df_marker_sample_all["ID"].to_list())
    row_df_reset_idx = row_df.reset_index(drop=False)
    row_df_reset_idx.iloc[6, 0] += "_1"
    row_df_reset_idx.iloc[-1, 0] += "_2"
    row_df_new = row_df_reset_idx.set_index("index")

    g = sns.clustermap(df_marker_sample_gene_cluster_sorted_plot,
                       row_colors=row_df_new,
                       col_colors=color_df,
                       colors_ratio=0.04,
                       cmap="viridis", 
                       row_cluster=False, 
                       col_cluster=False,
                       xticklabels=True, 
                       yticklabels=True, 
                       square=True,
                       figsize=(10, 5), 
                       cbar_kws={"ticks":list(range(0, 110, 25)), "orientation":'horizontal'})

    x0, _y0, _w, _h = g.cbar_pos
    g.ax_cbar.set_position([x0+0.7, -0.1, g.ax_row_dendrogram.get_position().width-0.1, 0.02])
    g.ax_cbar.set_title('$\\beta$-value', fontdict={"fontweight":"bold","fontsize":16})
    g.ax_cbar.tick_params(axis='x', length=0.5)
    g.ax_heatmap.set_xlabel("COVID-19 Patients", size=20)
    g.ax_heatmap.set_ylabel("", size=20)
    g.ax_heatmap.axvline(x=0, color="white", linewidth=5, linestyle="solid")
    g.ax_heatmap.axvline(x=35, color="white", linewidth=5, linestyle="solid")
    # g.ax_heatmap.axvline(x=46, color="white", linewidth=3, linestyle="solid")
    # g.ax_heatmap.axvline(x=83, color="white", linewidth=3, linestyle="solid")
    g.ax_heatmap.axhline(y=0, color="white", linewidth=5, linestyle="solid")
    g.ax_heatmap.axhline(y=17, color="white", linewidth=5, linestyle="solid")

    # from matplotlib.patches import Patch
    # legend_elements1 = [Patch(facecolor='pink', edgecolor='k',label='Acute'),
    #                     Patch(facecolor='powderblue', edgecolor='k',label='Recovery')]
    # legend_elements2 = [Patch(facecolor='firebrick', edgecolor='k',label='Severe'),
    #                     Patch(facecolor='forestgreen', edgecolor='k',label='Mild')]
    # legend_elements3 = [Patch(facecolor='grey', edgecolor='k',label='Hypomethylation'),
    #                     Patch(facecolor='black', edgecolor='k',label='Hypermethylation')]

    # legend1 = plt.legend(handles=legend_elements1, loc="center left", bbox_to_anchor=(-7.0, -0.1), frameon=False, fontsize=14)
    # legend1.set_title("Infection Phase", prop={"weight":"bold", "size":16})  
    # legend2 = plt.legend(handles=legend_elements2, loc="center left", bbox_to_anchor=(-5.0, -0.1), frameon=False, fontsize=14)
    # legend2.set_title("Severity Group", prop={"weight":"bold", "size":16})  
    # legend3 = plt.legend(handles=legend_elements3, loc="center left", bbox_to_anchor=(-3.0, -0.1), frameon=False, fontsize=14)
    # legend3.set_title("Expression", prop={"weight":"bold", "size":16})  
    # legend3.get_frame().set_alpha(None)
    # plt.gca().add_artist(legend1)
    # plt.gca().add_artist(legend2)

    # plt.gca().axes.xaxis.set_ticklabels([])
    # plt.rcParams['xtick.bottom'] = False
    # plt.rcParams['xtick.labelbottom'] = False
    plt.tight_layout()
    # plt.savefig(outfig, dpi=600, bbox_inches="tight")
    plt.show()
    plt.close()
    
    # ------

    # list_target_gene_sort = df_dmp_deg_sorted["Gene_Symbol"].to_list()
    # list_idx_gene_sort = list(range(len(list_target_gene_sort)))
    # dict_gene_sort = dict(zip(list_target_gene_sort, list_idx_gene_sort))
    
    # df_exp = pd.read_csv(fileexp, sep="\t")
    # df_exp["C19-C014-V1"] = np.nan
    # df_exp_idxed = df_exp.set_index("ID")
    # list_gene_id = df_exp_idxed.index
    # list_genesymbol = list(map(lambda x: "_".join(x.split("_")[1:]), list_gene_id))
    # list_all_samples = df_exp_idxed.columns
    # mtx_exp_scaled = StandardScaler().fit_transform(df_exp_idxed)
    # df_exp_scaled = pd.DataFrame(mtx_exp_scaled, index=list_gene_id, columns=list_all_samples)
    # list_expsamples = list(df_exp_scaled.columns)
    # list_colsamples = list(filter(lambda x: x in list_expsamples, list_sample_sort))
    # df_exp_selec = df_exp_scaled[list_colsamples]
    # df_exp_selec["gene_symbol"] = list_genesymbol
    # df_exp_selec_filt = df_exp_selec[df_exp_selec["gene_symbol"].isin(list_target_gene_sort)]
    # df_exp_selec_filt.index = list(map(lambda x: "_".join(x.split("_")[1:]) + "(" + x.split("_")[0] + ")", df_exp_selec_filt.index))
    # df_exp_selec_filt["Order"] = df_exp_selec_filt["gene_symbol"].map(dict_gene_sort)
    # df_exp_selec_filt_order = df_exp_selec_filt.sort_values(by=["Order"], ascending=True)
    # df_exp_selec_filt_order_rmv_dup = df_exp_selec_filt_order[~df_exp_selec_filt_order["gene_symbol"].duplicated(keep="first")]
    # df_exp_selec_filt_drop_genesym = df_exp_selec_filt_order_rmv_dup.drop(columns=["gene_symbol", "Order"])

    # sns.clustermap(df_exp_selec_filt_drop_genesym, cmap="rocket_r", row_cluster=False, col_cluster=False, yticklabels=True, xticklabels=True, linewidths=.5, square=True, robust=True, mask=df_exp_selec_filt_drop_genesym.isnull(), figsize=(20, 8))
    # plt.show()
    # plt.close()


# %%
def get_list_rna_sample(fileexp):
    df = pd.read_csv(fileexp, sep="\t")
    list_rna_samples = list(df.columns)

    return list_rna_samples

def get_dict_sample_severity(pathsevinfo, list_methyl_sample):
    dict_sample_severity = dict()
    with open (pathsevinfo, mode="r") as fr:
        list_header = fr.readline().rstrip("\n").split("\t")
        idx_id = list_header.index("Sample_ID")
        idx_sev_visit = list_header.index("Severity_visit")
        for line in fr:
            record = line.rstrip("\n").split("\t")
            sampleid = record[idx_id]
            sev_vis = record[idx_sev_visit]
            if sampleid in list_methyl_sample:
                dict_sample_severity[sampleid] = sev_vis

    return dict_sample_severity  

def sort_dict_sample_severity_visit(dict_sample_severity):
    dict_sample_severity_sort = dict(sorted(dict_sample_severity.items(), key=lambda x: (x[1].split("_")[1], x[1].split("_")[0])))
    
    return dict_sample_severity_sort

def get_dmp_deg_stats(filetargetgene, colmethdiff="dmp_median_methdiff", colfc="deg_median_log2FC"):
    colfcabs = "abs" + colfc
    df_dmp_deg = pd.read_csv(filetargetgene, sep="\t")
    df_dmp_deg[colfcabs] = df_dmp_deg[colfc].apply(abs)
    df_dmp_deg_sorted = df_dmp_deg.sort_values(by=[colmethdiff, colfcabs], ascending=False)
    
    return df_dmp_deg_sorted

def get_dict_target_marker_gene(df_dmp_deg_sorted, colmeth="Methyl", colrna="RNA"):
    dict_marker_gene = defaultdict(set)
    for _, rows in df_dmp_deg_sorted.iterrows():
        dict_rows = dict(rows)
        gene_sym = dict_rows[colrna]
        dmp_markers = ":".join(dict_rows[colmeth].split("_")[:-1])
        dict_marker_gene[dmp_markers].add(gene_sym)
    
    return dict(dict_marker_gene)

def get_dict_target_marker_direction(df_dmp_deg_sorted, colmeth="Methyl", colmethdiff="dmp_median_methdiff"):
    dict_marker_gene = dict()
    for _, rows in df_dmp_deg_sorted.iterrows():
        dict_rows = dict(rows)
        direction = dict_rows[colmethdiff]
        if np.sign(direction) > 0:
            direction = "hyper"
        else:
            direction = "hypo"
        dmp_markers = ":".join(dict_rows[colmeth].split("_")[:-1])
        dict_marker_gene[dmp_markers] = direction
    
    return dict(dict_marker_gene)

def get_samplewise_genewise_cpg(outfilename, dict_markers, dict_sample_severity, dict_marker_gene):
    with open(outfilename, mode='w') as fw:
        list_header = ["ID", "Severity", "Marker", "GeneSymbol", "CpGbeta"]
        fw.write("\t".join(list_header) + "\n")
        for sample_id, dict_marker in dict_markers.items():
            if sample_id in dict_sample_severity.keys():
                severity_visit = dict_sample_severity[sample_id]
            else:
                continue
            for marker_name, cpg in dict_marker.items():
                pos = marker_name[0] + ":" + marker_name[1]
                if pos in dict_marker_gene.keys():
                    list_gene = dict_marker_gene[pos]
                    for gene in list_gene:
                        list_content = [sample_id, severity_visit, pos, gene, cpg]
                        fw.write("\t".join(list_content) + "\n")

def modify_samplewise_genewise_cpg_table(outfile):
    df_cpg_sample_all = pd.read_csv(outfile, sep="\t")
    df_cpg_sample_all["GeneSymbol_plot"] = df_cpg_sample_all["GeneSymbol"].apply(lambda x: "_".join(x.split("_")[1:]))
    df_cpg_sample_all["Marker_plot"] = df_cpg_sample_all["GeneSymbol_plot"] + " (" + df_cpg_sample_all["Marker"] + ")"
    df_cpg_sample_all["Severity_Group"] = df_cpg_sample_all["Severity"].str.split("_").str[0]
    df_cpg_sample_all["Visit"] = df_cpg_sample_all["Severity"].str.split("_").str[1]

    return df_cpg_sample_all


def get_list_cluster_markers(df_marker_sample_gene_cluster, n_clust):
    import numpy as np
    from sklearn.cluster import AgglomerativeClustering
    data = df_marker_sample_gene_cluster.values
    cluster = AgglomerativeClustering(n_clusters=n_clust, affinity='euclidean', linkage='average')
    list_clust = cluster.fit_predict(data)

    return list_clust
# %%
if __name__ == "__main__":
    visit = "first"
    pathsevinfo = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/9_clinical/Infectomics_Severity_Information_Methyl_20240102.tsv"
    filetargetgene = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/Methyl_RNA_Relation/Methyl_RNA_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Mild_Sev_Visit1.LOO_common.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_6.DMP_DEG_LOO_annot.20240402.tsv"
    fileexp = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/5_deg/Visit1_Severe__Visit1_Mild_20240327.tsv.normcount" 
    outfig = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/InfectomicsPaper1/methylation_heatmap_severity_and_infectionphase.png"
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
                        'C19-C061-V3'
                        ]
        
    main(pathsevinfo, fileexp, filetargetgene, outfig, list_drop_samples)
# %%
