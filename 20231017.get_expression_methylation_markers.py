# %%
import math
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu, ttest_ind
from statannot import add_stat_annotation
from statsmodels.stats.multitest import fdrcorrection


# %%
def main(path_exp, path_methyl_marker, path_meta, list_drop_samples, list_annotation, outdir):
    df_TPM = read_TPM_matrix(path_exp, list_drop_samples)

    list_filter_genes = list()
    if list_filter_genes == list():
        list_methyl_markers = get_list_methylation_markers(path_methyl_marker)
        df_TPM_gene_filtered = filter_list_genes(df_TPM, list_methyl_markers, mode="genesymbol")
        df_TPM_transposed = transpose_TPM_matrix(df_TPM_gene_filtered)
        df_TPM_split_phase = split_sample_by_phase(df_TPM_transposed)
        df_TPM_meta = merge_TPM_matrix_meta(df_TPM_split_phase, path_meta)
        df_TPM_meta = attach_recovered_severity(df_TPM_meta)
        list_select_gene = df_TPM_gene_filtered["ID"].to_list()
        list_filter_genes.extend(list_select_gene)
        df_gene_visit_sorted = make_dataframe_stat_test(df_TPM_meta, list_select_gene, outdir)
        # df_gene_visit_sorted_sig_only = df_gene_visit_sorted[df_gene_visit_sorted["testsignificant"]==True]
        

    df_TPM_gene_filtered = filter_list_genes(df_TPM, list_filter_genes, mode="geneonly")
    df_TPM_transposed = transpose_TPM_matrix(df_TPM_gene_filtered)
    df_TPM_split_phase = split_sample_by_phase(df_TPM_transposed)
    df_TPM_meta = merge_TPM_matrix_meta(df_TPM_split_phase, path_meta)
    df_TPM_meta = attach_recovered_severity(df_TPM_meta)
    list_xticklabels = get_xticklabels(df_TPM_meta, colsample="ID", colsev="Severity_group", coltime="Visit")
    list_visit = list(map(lambda x:x.split("\n")[0], list_xticklabels))
    draw_mean_difference(df_TPM_meta, df_gene_visit_sorted, list_select_gene, list_annotation, outdir)

######
    # df_TPM_meta_visit2 = df_TPM_meta[np.logical_or(df_TPM_meta["Visit_order"] == "Visit2", df_TPM_meta["Visit_order"] == "Recover")]
    # for gene in list_filter_genes:
    #     df_sev_mean_gene_exp = df_TPM_meta_visit2.groupby("Severity_group")[gene].apply(np.mean).reset_index(drop=False).rename(columns={"index": gene})
    #     mean_severe = df_sev_mean_gene_exp[df_sev_mean_gene_exp["Severity_group"]=="Severe"][gene].to_numpy()[0]
    #     mean_mild = df_sev_mean_gene_exp[df_sev_mean_gene_exp["Severity_group"]=="Mild"][gene].to_numpy()[0]
    #     diff_mean = (mean_severe - mean_mild)
    #     if diff_mean > 0:
    #         print(gene, diff_mean)
# #########
#     sig_visit_cnt = df_gene_visit_sorted_sig_only.groupby("Visit").apply(len).reset_index(drop=False)
#     sig_visit_cnt.columns=["Visit_order", "Count"]
#     ax = sns.barplot(data=sig_visit_cnt, x="Visit_order", y="Count", order=list_visit)
#     for i in ax.containers:
#         ax.bar_label(i,)
#     plt.yscale("log")
#     plt.title("Number of Visits with Sigificant Difference in mean Expression \n between Severe and Mild Patients for Methylation Markers")

#     plt.show()
# #########
#     list_filter_genesymbols = list(map(lambda x: "_".join(x.split("_")[1:]), list_filter_genes))
#     print(list_filter_genesymbols)

# %%
def read_TPM_matrix(path_exp, list_drop_samples):
    df_TPM = pd.read_csv(path_exp, sep="\t")
    df_TPM  = df_TPM.drop(columns=list_drop_samples)

    return df_TPM

def select_TPM_matrix(df_TPM, select_pattern, delim_id="-", namepos=0):
    list_colname = list(df_TPM.columns)
    list_colname_filt = list(filter(lambda x: x.split(delim_id)[namepos][0] != select_pattern if len(x.split("-")) > 1 else x, list_colname))
    df_TPM = df_TPM.loc[:, list_colname_filt]

    return df_TPM

def get_list_methylation_markers(path_methyl_marker):
    list_methyl_markers = list()
    with open(path_methyl_marker, mode="r") as fr:
        for line in fr:
            record = line.rstrip("\n").split("\n")
            list_methyl_markers.extend(record)
    
    return list_methyl_markers

def filter_list_genes(df_TPM, list_gene, mode, colgene="ID"):
    if len(list_gene) > 0:
        if mode == "genesymbol":
            df_TPM["tmp"] = df_TPM[colgene].apply(lambda x: "_".join(x.split("_")[1:]))
            df_TPM_gene_filtered = df_TPM[df_TPM["tmp"].isin(list_gene)].drop(columns=["tmp"])
        if mode == "geneonly":
            df_TPM_gene_filtered = df_TPM[df_TPM[colgene].isin(list_gene)].drop(columns=["tmp"])

    else:
        df_TPM_gene_filtered = df_TPM
        
    return df_TPM_gene_filtered

def transpose_TPM_matrix(df_TPM, colsample="ID"):
    list_geneid = df_TPM[colsample].to_list()
    df_TPM_transposed = df_TPM.T.iloc[1:, :]
    df_TPM_transposed.columns = list_geneid
    df_TPM_transposed = df_TPM_transposed.reset_index(drop=False).rename(columns={"index": colsample})
    
    return df_TPM_transposed

def get_visit_order(sampleid, dict_manual={"0": "V0", "R": "V5", "L": "V6"}):
    visit = sampleid.split("-")[-1]
    for key, val in dict_manual.items():
        if sampleid.split("-")[1][0] == key:
            visit = val
        if sampleid.split("-")[-1][0] == key:
            visit = val
    visit_order = int(visit[1:])

    return visit_order

def split_sample_by_phase(df_TPM_gene_filtered, colsample="ID", colTPM="TPM", coltime="Time"):
    df_split_phase = df_TPM_gene_filtered.copy()
    df_split_phase[coltime] = df_split_phase[colsample].apply(get_visit_order)
    df_split_phase = df_split_phase.sort_values(by=coltime, ascending=True)
    df_split_phase[coltime] = df_split_phase[coltime].astype(str)

    return df_split_phase
    

def merge_TPM_matrix_meta(df_TPM_split_phase, path_meta, colsample="ID"):
    df_meta = pd.read_csv(path_meta, sep="\t")
    df_meta_rename_colsample = df_meta.rename(columns={"Sample_ID": colsample})
    df_meta_rename_colsample_set_idx = df_meta_rename_colsample.set_index(colsample)
    # SAMPLES WITH NO RNA DATA
    df_meta_rename_colsample_drop_no_exp = df_meta_rename_colsample_set_idx.drop(index=["C19-C014-V1", "C19-C059-V2"])
    df_meta_rename_colsample_drop_no_exp = df_meta_rename_colsample_drop_no_exp.reset_index(drop=False).rename(columns={"index": colsample})
    df_merged = pd.merge(df_TPM_split_phase, df_meta_rename_colsample_drop_no_exp, how="outer", on=colsample)
    df_merged = df_merged[df_merged["ID"].str.split("-").str[-1].str.len() < 3]
    # df_merged["Severity_visit"] = df_merged["Severity_visit"].fillna("Healthy_Control")
    # df_merged["Severity_group"] = df_merged["Severity_group"].fillna("Control")
    # df_merged["Visit"] = df_merged["Visit"].fillna("Healthy")
    
    return df_merged

def attach_recovered_severity(df_TPM_meta):
    df_TPM_meta["Severity"] = df_TPM_meta["Severity"].apply(lambda x: 0 if x =="-" else x)
    df_TPM_meta[df_TPM_meta["ID"].str.contains("C19-R")]["Severity_group"] = "Convalescent"
    df_TPM_meta[df_TPM_meta["ID"].str.contains("L1")]["Severity_group"] = "LongCOVID"

    return df_TPM_meta

def meandiff(a, b):
    deltamean = np.mean(b) - np.mean(a)

    return deltamean

def make_dataframe_stat_test(df_TPM_meta, list_select_gene, outdir):
    dict_gene_visit_pval = dict()
    for gene in list_select_gene:
        df_group_visit = df_TPM_meta.groupby(["Visit", "Severity_group"])[gene].apply(np.array).reset_index(drop=False)

        dict_visit_pval = dict()
        for ind1 in list(df_group_visit.index):
            if ind1 < len(df_group_visit)-1:
                ind2 = int(ind1) + 1
            comp1 = df_group_visit.loc[ind1, :] 
            comp2 = df_group_visit.loc[ind2, :]
            visit_num1 = comp1["Visit"]
            visit_num2 = comp2["Visit"]
            if visit_num1 == visit_num2:
                list_case = comp1[gene]
                list_case = list(map(float, list_case))
                list_control = comp2[gene]
                list_control = list(map(float, list_control))

                if len(set(list_control)) > 1:
                    deltamean = meandiff(list_case, list_control)
                    stat, pval = ttest_ind(list_case, list_control, equal_var=False)
                    stat_res = {"delta": deltamean, "stat": stat, "pval": pval}
                    dict_visit_pval[visit_num1] = stat_res

        if dict_gene_visit_pval.get(gene) == None:
            dict_gene_visit_pval[gene] = dict()
        
        dict_gene_visit_pval[gene].update(dict_visit_pval)

    df_gene_visit = pd.concat({k: pd.DataFrame(v).T for k, v in dict_gene_visit_pval.items()}, axis=0)
    df_gene_visit = df_gene_visit.reset_index(drop=False).rename(columns={"level_0": "ID", "level_1": "Visit"})

    df_gene_visit_sorted = df_gene_visit.sort_values(by=["pval"], ascending=True)
    list_pval = df_gene_visit_sorted["pval"].to_list()
    list_sig = fdrcorrection(list_pval)[0]
    list_padj = fdrcorrection(list_pval)[1]
    df_gene_visit_sorted["padj"] = list_padj
    df_gene_visit_sorted["testsignificant"] = list_sig
    df_gene_visit_sorted.to_csv(os.path.join(outdir, "stattest_expressed_methylation_markers.tsv"), sep="\t", index=False)

    return df_gene_visit_sorted

def get_xticklabels(df_TPM_meta, colsample="ID", colsev="Severity_group", coltime="Visit"):
    df_num = df_TPM_meta.groupby([colsev, coltime])[colsample].unique().apply(len).reset_index(drop=False)
    list_visit = df_TPM_meta[coltime].unique()

    list_xticklabels = list()
    for visit in list_visit:
        list_label_visit = list()
        for ind, row in df_num.iterrows():
            if str(row[coltime]) == str(visit):
                list_label_elem = row.to_list()
                label_visit = list_label_elem[0] + f"\n(N={list_label_elem[-1]})"
                list_label_visit.append(label_visit)
        list_xticklabels.append(list_label_visit)
    list_xticklabels = list(map(lambda x: "\n".join(x), list_xticklabels))

    list_new_xticklabels = list(map(lambda x: "\n\n".join(x), list(zip(list_visit, list_xticklabels))))
    list_new_xticklabels = list(map(lambda x: f"{x}\n", list_new_xticklabels))

    return list_new_xticklabels

   
def draw_mean_difference(df_TPM_meta, df_gene_visit_sorted, list_select_gene, list_annotation, outdir):
    list_xticklabels = get_xticklabels(df_TPM_meta)
    list_visit = list(map(lambda x:x.split("\n")[0], list_xticklabels))
    dict_sort = dict(zip(list_visit, list(range(len(list_visit)))))

    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#D68C78', edgecolor='k',label='Severe'), Patch(facecolor='#ECECEC', edgecolor='k',label='Mild')]


    dict_severity_visit_order = dict(zip(['Mild_First', 'Severe_First', 'Mild_Second',
       'Severe_Second', 'Mild_Third', 'Severe_Third', 'Mild_Fourth',
       'Severe_Fourth', 'Mild_Convalescent', 'Severe_Convalescent',
       'Mild_LongCOVID', 'Severe_LongCOVID'], range(12)))

    df_TPM_meta_order = df_TPM_meta.copy()
    df_TPM_meta_order["X_order"] = df_TPM_meta_order["Severity_visit"].apply(dict_severity_visit_order.__getitem__)
    df_TPM_meta_order["Severity_visit"] = df_TPM_meta_order.apply(lambda x : x["Severity_visit"] if pd.notna(x["Severity_visit"]) else x["X_order"], axis = 1)
    df_TPM_meta_order = df_TPM_meta_order.sort_values(by = ["X_order"])

    palette = dict()
    for sev in df_TPM_meta["Severity_visit"].unique():
        if "Mild" in sev:
            palette[sev] = "#ECECEC"
        elif "Severe" in sev:
            palette[sev] = "#D68C78"
        else:
            palette[sev] = "#75A7C3"

    palette_flexible = dict()
    for key in df_TPM_meta_order["Severity_visit"].unique():
        if palette.get(key) == None:
            palette_flexible[key] = 'k'
        else:
            palette_flexible[key] = palette[key]

    from itertools import combinations
    fig, axes = plt.subplots(ncols = math.ceil(math.sqrt(len(list_select_gene))), nrows = math.ceil(math.sqrt((len(list_select_gene)))), figsize=(30, 20))
    ax = axes.flatten()
    for ind, gene in enumerate(list_select_gene):
        plot = sns.boxplot(data=df_TPM_meta_order, x="Severity_visit", y=gene, palette=palette_flexible, ax=ax[ind])
        list_xticks = plot.get_xticklabels()
        list_hide = list(map(str,[0.5, 2.5, 4.5, 6.5, 8.5, 10.5]))
        dict_tickindex_to_tick = dict(zip(list(range(0, len(df_TPM_meta_order["Severity_visit"].unique()))), list(map(lambda x : x.get_text() if x.get_text() not in list_hide else "", list_xticks))))
        tickindex_to_show = list(filter(lambda x : dict_tickindex_to_tick[x] != "", dict_tickindex_to_tick.keys()))
        tickvalue_to_show = list(filter(lambda x : x != "", dict_tickindex_to_tick.values()))


        if "_".join(gene.split("_")[1:]) in list_annotation:
            ax[ind].set_title(gene, fontsize = 11, color="r")
        else:
            ax[ind].set_title(gene, fontsize = 11, color="k")

        list_categories = list(df_TPM_meta.groupby(["Severity_visit"])[gene].apply(np.array).index)
        list_combination = list(combinations(list_categories, 2))
        list_combination_same_visit = list(filter(lambda x: x[0].split("_")[-1] == x[-1].split("_")[-1], list_combination))
        list_comb_same_visit_sorted = sorted(list_combination_same_visit, key=lambda x: dict_sort[x[0].split("_")[-1]])
        # pairs = list_comb_same_visit_sorted

        list_visit = list(map(lambda x:x.split("\n")[0], list_xticklabels))
        df_gene_visit_select = df_gene_visit_sorted[df_gene_visit_sorted["ID"]==gene]
        dict_visit_pval = dict(zip(df_gene_visit_select["Visit"], df_gene_visit_select["padj"]))
        dict_visit_pval_sorted = dict(sorted(dict_visit_pval.items(), key=lambda x: dict_sort[x[0]]))
        list_visit_pval_sorted = list(dict_visit_pval_sorted.values())
        # pvalues = list_visit_pval_sorted

        # add_stat_annotation(ax[ind], data=df_TPM_meta, x="Severity_visit", y=gene,
        #                     box_pairs=pairs,perform_stat_test=False, pvalues=pvalues,test=None, text_format='star', loc='inside', verbose=2)

        ax[ind].set_xlabel("TimePoints", fontsize=12)
        ax[ind].set_xticks(tickindex_to_show, tickvalue_to_show, fontsize=10)
        ax[ind].set_ylabel("TPM", fontsize=12)

    outfilename = os.path.join(outdir, "boxplot_expression_methylation_markers.png")
    plt.tight_layout()
    plt.savefig(outfilename, dpi=600)
    plt.show()
    plt.close()

# %%
if __name__ == "__main__":
    direction = "hyper"
    visit = "first"
    path_exp = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/4_expmtx/ConfirmedRecovered/expression_matrix_genes.results_TPM.tsv"
    # path_methyl_marker = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{visit}/{direction}/target_genelist.tsv"
    path_methyl_marker = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{visit}/{direction}/ranomd_genelist.tsv"
    path_meta = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/9_clinical/Infectomics_Severity_Information_RNA_20240102.tsv"
    list_drop_samples = []
    list_annotation = []
    list_annotation = []
    list_drop_samples = ["C19-C045-V2",
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

    outdir = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/{visit}/{direction}"
    os.makedirs(outdir, exist_ok=True)

# %%
    main(path_exp, path_methyl_marker, path_meta, list_drop_samples, list_annotation, outdir)
# %%
