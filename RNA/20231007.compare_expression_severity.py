 # %%
import math
import os
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.preprocessing import StandardScaler


# %%
def main(path_exp, path_deg, path_meta, outdir, list_drop_samples, list_upreg_genes, list_downreg_genes, log2fc_thres, padj_thres, colsample, colTPM, colsev, coltime, nareplace):
    df_TPM = read_TPM_matrix(path_exp, list_drop_samples)
    df_TPM = select_TPM_matrix(df_TPM, select_pattern="L", delim_id="-", namepos=2)
    df_summary = make_summary_table_deg(df_TPM, path_deg, log2fc_thres, padj_thres, outdir)
    df_stats = report_statistics(df_summary, outdir)
    df_merged = merge_metadata(df_summary, path_meta)
    df_merged.loc[(df_merged["Severity_Binary"].isna()) & df_merged["ID"].str.contains("C19-R"), "Severity_Binary"] = 'Convalescent'
    df_merged.loc[(df_merged["Severity_Binary"].isna()) & df_merged["ID"].str.contains("U10K"), "Severity_Binary"] = 'Healthy'
    df_merged_downreg = separate_dataframe_by_direction(df_merged, "downreg")
    df_merged_upreg = separate_dataframe_by_direction(df_merged, "upreg")

    ax = sns.catplot(kind="point", data=df_merged_downreg, x="Time", y="NormTPM", hue="Severity_Binary", errorbar=('ci', 95), order=["0", "1", "2", "3", "4", "5"], palette={"Healthy":"#2E8B57" ,"Severe":"#D68C78", "Mild":"#aaaaaa", "Convalescent":"#75A7C3"}, capsize=.15, errwidth=0.5, legend=False, height=5, aspect=7/5)
    list_xticklabels = get_xticklabels(df_merged_downreg)[::-1]
    ax.set_xticklabels(list_xticklabels, fontsize=8)
    figpathpoint = os.path.join(outdir, "severity_NormTPM_downreg_pointplot.png")
    plt.legend(loc = "upper left", bbox_to_anchor = (1.05, 0.6))
    plt.xlabel("Visit", fontsize=11)
    plt.ylabel("Gene-wise Normalized TPM", fontsize=11)
    plt.tight_layout()
    plt.show()
    plt.savefig(figpathpoint, dpi=300)
    plt.close()

    #     ax = sns.catplot(kind="point", data=df_merged_upreg, x="Time", y="NormTPM", hue="Severity_Binary", errorbar=('ci', 95), order=["1", "2", "3", "4", "5"], capsize=.15, errwidth=0.5, legend=False, height=8.27, aspect=11.7/8.27)
    #     list_xticklabels = get_xticklabels(df_merged_upreg)[::-1]
    #     ax.set_xticklabels(list_xticklabels, fontsize=8)
    #     figpathpoint = os.path.join(outdir, "severity_NormTPM_upreg_pointplot.png")
    #     plt.legend(loc = "upper left", bbox_to_anchor = (1.05, 0.6))
    #     plt.tight_layout()
    #     plt.savefig(figpathpoint, dpi=300)
    #     plt.close()

    #     ax = sns.catplot(kind="box", data=df_merged_downreg, x="Time", y="NormTPM", hue="Severity_Binary", order=["1", "2", "3", "4", "5"], legend=False, height=8.27, aspect=11.7/8.27)
    #     list_xticklabels = get_xticklabels(df_merged_downreg)[::-1]
    #     ax.set_xticklabels(list_xticklabels, fontsize=8)
    #     figpathbox = os.path.join(outdir, "severity_NormTPM_downreg_boxplot.png")
    #     plt.legend(loc = "upper left", bbox_to_anchor = (1.05, 0.6))
    #     plt.ylim(top=30)
    #     plt.tight_layout()
    #     plt.savefig(figpathbox, dpi=300)
    #     plt.close()

    #     ax = sns.catplot(kind="box", data=df_merged_upreg, x="Time", y="NormTPM", hue="Severity_Binary", order=["1", "2", "3", "4", "5"], legend=False, height=8.27, aspect=11.7/8.27)
    #     list_xticklabels = get_xticklabels(df_merged_upreg)[::-1]
    #     ax.set_xticklabels(list_xticklabels, fontsize=8)
    #     figpathbox = os.path.join(outdir, "severity_NormTPM_upreg_boxplot.png")
    #     plt.legend(loc = "upper left", bbox_to_anchor = (1.05, 0.6))
    #     plt.ylim(top=30, bottom=0)
    #     plt.tight_layout()
    #     plt.savefig(figpathbox, dpi=300)
    #     plt.close()


        # plot_gene_severity(dict_deg, df_merged_downreg, path_analysis, tag="downreg")
        # plot_gene_severity(dict_deg, df_merged_upreg, path_analysis, tag="upreg")
        # plot_total_tpm_severity(df_merged_downreg, path_analysis, tag="downreg", plottype="pointplot", logscale=True)
        # plot_total_tpm_severity(df_merged_upreg, path_analysis, tag="upreg", plottype="pointplot", logscale=True)
        # plot_total_tpm_severity(df_merged_downreg, path_analysis, tag="downreg", plottype="boxplot", logscale=True)
        # plot_total_tpm_severity(df_merged_upreg, path_analysis, tag="upreg", plottype="boxplot", logscale=True)
        # plot_trend(df_merged_downreg, dict_deg, path_analysis, tag="downreg")
        # plot_trend(df_merged_upreg, dict_deg, path_analysis, tag="upreg")
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

def get_dict_degs(path_deg, log2fc_thres, padj_thres):
    dict_deg = dict()
    with open(path_deg, mode="r") as fr:
        header = fr.readline().rstrip("\n").split("\t")
        idx_gene = header.index("ID")
        idx_log2fc = header.index("log2FoldChange")  
        idx_padj = header.index("padj") 
        for line in fr:
            record = line.rstrip("\n").split("\t")
            gene = record[idx_gene]
            log2fc = record[idx_log2fc]
            padj = record[idx_padj]

            if (abs(float(log2fc)) > log2fc_thres and str(padj) != "NA" and float(padj) < float(padj_thres)):
                dict_deg[gene] = dict()
                dict_deg[gene] = {"log2fc":log2fc, "padj":padj}
    
    return dict_deg

def select_gene_of_interest(df_TPM, target_gene, mode, norm=True, colgene="ID", colTPM="TPM"):
    if mode == "geneonly":
        df_gene = df_TPM[df_TPM[colgene] == target_gene]
        df_selected = df_gene.T.iloc[1:, :].reset_index(drop=False)
        df_selected.columns = [colgene, colTPM]
    if mode == "genesymbol":
        df_TPM["tmp"] = df_TPM[colgene].apply(lambda x: "_".join(x.split("_")[1:]))
        df_gene = df_TPM[df_TPM["tmp"] == target_gene].drop(columns=["tmp"])
        df_selected = df_gene.T.iloc[1:, :].reset_index(drop=False)
        df_selected.columns = [colgene, colTPM]
    if mode == "noversion":
        df_TPM["tmp"] = df_TPM[colgene].apply(lambda x: x.split(".")[0] + "_" + "_".join(x.split("_")[1:]))
        df_gene = df_TPM[df_TPM["tmp"] == target_gene].drop(columns=["tmp"])
        df_selected = df_gene.T.iloc[1:, :].reset_index(drop=False)
        df_selected.columns = [colgene, colTPM]
    if norm:
        df_selected["MeanTPMGeneWise"] = df_selected[colTPM].mean()
        df_selected["NormTPM"] = df_selected[colTPM]/df_selected["MeanTPMGeneWise"]

    return df_selected

def get_visit_order(sampleid, namepos=1, dict_manual={"R": "V5", "0": "V0"}):
    visit = sampleid.split("-")[-1]
    for key, val in dict_manual.items():
        if str(sampleid.split("-")[namepos][0])== key:
            visit = val
    visit_order = int(visit[1:])

    return visit_order

def split_sample_by_phase(df_selected, colsample="ID", colTPM="TPM", coltime="Time"):
    df_split_phase = df_selected.copy()
    df_split_phase[coltime] = df_split_phase[colsample].apply(get_visit_order)
    df_split_phase = df_split_phase.sort_values(by=coltime, ascending=False)
    df_split_phase[coltime] = df_split_phase[coltime].astype(str)
    df_split_phase[colsample] = df_split_phase[colsample].apply(lambda x:"-".join(x.split("-")[:-1]))

    return df_split_phase

def add_gene_info_column_dataframe(df, gene, deg_info, colgene="Gene", colgeneinfo="Gene_info", coldirec="Direction"):
    df_modif = df.copy()
    log2fc = round(float(deg_info["log2fc"]), 3)
    padj = round(float(deg_info["padj"]), 3)
    gene_info = f"{gene}\nlog2FC = {log2fc}\nPadj = {padj}"
    df_modif[colgene] = gene
    df_modif[colgeneinfo] = gene_info
    if float(log2fc) < 0:
        direction = "negative"
    else:
        direction = "positive"
    df_modif[coldirec] = direction

    return df_modif

def make_summary_table_deg(df_TPM, path_deg, log2fc_thres, padj_thres, outdir):
    list_df_split_phase = list()
    dict_deg = get_dict_degs(path_deg, log2fc_thres, padj_thres)   
    for gene, deg_info in dict_deg.items():
        df_selected = select_gene_of_interest(df_TPM, gene, mode="geneonly")
        df_split_phase = split_sample_by_phase(df_selected)
        df_split_phase_modified = add_gene_info_column_dataframe(df_split_phase, gene, deg_info)
        list_df_split_phase.append(df_split_phase_modified)
    df_summary = pd.concat(list_df_split_phase, axis=0)
    df_summary.to_csv(os.path.join(outdir, "deg_results.txt"), sep="\t", index=False)
    
    return df_summary

def report_statistics(df_summary, outdir):
    downreg_genes = list(df_summary[df_summary["Direction"]=="negative"]["Gene"].unique())
    upreg_genes = list(df_summary[df_summary["Direction"]=="positive"]["Gene"].unique())
    dict_stat = {"downreg":{"number":len(downreg_genes), "genes":downreg_genes}, "upreg":{"number":len(upreg_genes), "genes":upreg_genes}}
    df_stat = pd.DataFrame.from_dict(dict_stat, orient='index')

    df_stat.to_csv(os.path.join(outdir, "deg_stats.txt"), sep="\t")

    return df_stat

def merge_metadata(df_summary, path_meta, colsample="ID", colsev="Severity_Binary"):
    df_meta = pd.read_csv(path_meta, sep="\t")
    df_meta = df_meta.rename(columns={"SampleID": colsample})
    df_merged = pd.merge(df_summary, df_meta, how="outer", on=colsample)

    return df_merged

def separate_dataframe_by_direction(df_merged, tag):
    if tag == "downreg":
        df_merged_reg = df_merged[df_merged["Direction"]=="negative"]
    if tag == "upreg":
        df_merged_reg = df_merged[df_merged["Direction"]=="positive"]

    return df_merged_reg

def get_xticklabels(df_merged, colsample="ID", colsev="Severity_Binary", coltime="Time"):
    df_num = df_merged.groupby([colsev, coltime])[colsample].unique().apply(len).reset_index(drop=False)
    df_num = df_num.sort_values(by=[coltime], ascending=True).reset_index(drop=True)
    list_visit = df_merged[coltime].unique()

    from itertools import zip_longest
    list_xticklabels = list()
    for visit in list_visit:
        list_label_visit = list()
        for x in df_num.iterrows():
            if str(x[1][coltime]) == str(visit):
                list_label_elem = x[1].to_list()
                label_visit = list_label_elem[0] + f"\n(N={list_label_elem[-1]})"
                list_label_visit.append(label_visit)
        list_xticklabels.append(list_label_visit)
    list_xticklabels = list(map(lambda x: "\n".join(x), list_xticklabels))
    list_xticklabels = list(map(lambda x: "\n".join(x), list(zip_longest(["Healthy", "First", "Second", "Third", "Fourth", "Recovered"][::-1], list_xticklabels))))

    return list_xticklabels

def plot_gene_severity(dict_deg, df_merged, path_analysis, tag, hue_list = "Set2", colTPM="TPM", colsev="Severity_Binary", coltime="Time", colgene="Gene"):
    if tag == "downreg":
        list_deg = list(filter(lambda x : float(dict_deg[x]["log2fc"])<0, dict_deg.keys()))
        sorted_list_deg = sorted(list_deg, key = lambda x : float(dict_deg[x]["log2fc"]), reverse = True)
    elif tag == "upreg":
        list_deg = list(filter(lambda x : float(dict_deg[x]["log2fc"])>0, dict_deg.keys()))
        sorted_list_deg = sorted(list_deg, key = lambda x : float(dict_deg[x]["log2fc"]), reverse = True)
    else:
        raise Exception

    num_col = int(math.sqrt(len(dict_deg))*0.7)

    fig, ax = plt.subplots(ncols = num_col, nrows = math.ceil(len(sorted_list_deg)/num_col))
    axes = ax.flatten()

    severity_order = sorted(df_merged[colsev].unique())

    for ind, geneid in enumerate(sorted_list_deg):
        geneid = geneid.split('_')[-1]
        table_part = df_merged[df_merged[colgene] == geneid].copy()
        sns.boxplot(data = table_part, x = coltime, y = colTPM, hue = colsev, hue_order=severity_order, palette=hue_list[ind])
        axes[ind].set_title(geneid, fontsize = 10)
        list_xticklabels = get_xticklabels(df_merged)
        axes[ind].set_xlabel("TimePoints", fontsize=10)
        axes[ind].set_xticklabels(list_xticklabels, fontsize=10)
        axes[ind].set_ylim(bottom = 0)
        if ind == num_col-1:
            axes[ind].get_legend().set_visible(True)
            axes[ind].legend(loc = "upper left", bbox_to_anchor = (1.05, 0.95))
        else:
            axes[ind].get_legend().set_visible(False)

    plt.tight_layout()
    plt.savefig(os.path.join(path_analysis, f"severity_TPM_each_gene_{tag}.png"), dpi=300)
    plt.show()
    plt.close()

def plot_total_tpm_severity(df_merged, path_analysis, tag, plottype, logscale=True, colsample="ID", colsev="Severity_Binary", colTPM="TPM", coltime="Time"):
    plt.figure(figsize=(12,5))
    df_mean_tpm_severity = df_merged.groupby([colsev, coltime, colsample])[colTPM].mean().reset_index(drop=False)
    df_std_tpm_severity = df_merged.groupby([colsev, coltime, colsample])[colTPM].std().reset_index(drop=False)
    if plottype == "boxplot":
        ax = sns.boxplot(data=df_mean_tpm_severity, x=coltime, y=colTPM, hue=colsev)
    if plottype == "pointplot":
        ax = sns.pointplot(data=df_mean_tpm_severity, x=coltime, y=colTPM, hue=colsev, dodge=0.5, err_wdith=0)
    list_xticklabels = get_xticklabels(df_merged)
    ax.set_xticklabels(list_xticklabels, fontsize=8)
    ax.set_xlabel("TimePoints", fontsize=12)
    ax.set_ylabel("SumTPM", fontsize=12)
    if logscale:
        plt.yscale("log")
        ax.set_ylabel("$log10$(SumTPM)", fontsize=12)
    plt.legend(loc = "upper left", bbox_to_anchor = (1.05, 0.95))
    plt.tight_layout()
    plt.savefig(os.path.join(path_analysis, f"severity_SumTPM_across_time_{tag}_{plottype}.png"), dpi=300)
    plt.show()
    plt.close()

def plot_trend(df, dict_deg, path_analysis, tag="downreg", split_by="Severity_Binary", colgene="Gene", colTPM="TPM", coltime="Time", coltrial="Trial", coldirec="Direction"):
    # https://github.com/webermarcolivier/statannot/issues/27
    # https://github.com/mwaskom/seaborn/issues/1166
    order = df[coltime].unique()
    hue_order = df[split_by].unique()
    palette = "Set2"
    figure = sns.FacetGrid(data=df, col=colgene, col_wrap=int(math.sqrt(len(dict_deg))), hue=split_by, sharex=False, sharey=False, legend_out=True)
    figure.map(sns.violinplot, coltime, colTPM, split_by, palette=palette, hue_order=hue_order, order=order)
    labels = hue_order
    colors = sns.color_palette(palette).as_hex()[:len(labels)]
    handles = [patches.Patch(color=col, label=lab) for col, lab in zip(colors, labels)]
    plt.legend(handles=handles, title=colsev, loc="center left", bbox_to_anchor=(1, 0.75))
    plt.savefig(os.path.join(path_analysis, f"default_split_by_{split_by}_{tag}.png"), dpi=300)
    plt.show()
    plt.close()

    # try:
    #     figure = sns.FacetGrid(data=df, col=colgene, hue=coltrial, col_wrap=int(math.sqrt(len(dict_deg))), sharex=False, sharey=False, legend_out=True)
    #     figure.map(sns.lineplot, coltime, colTPM)
    #     plt.savefig(os.path.join(path_analysis, "div_by_trial.png"), dpi=300)
    #     plt.close()
    # except:
    #     pass

    # try:
    #     figure = sns.FacetGrid(data=df, col=colgene, hue=coldirec, col_wrap=int(math.sqrt(len(dict_deg))), sharex=False, sharey=False, legend_out=True)
    #     figure.map(sns.violinplot, coltime, colTPM)
    #     figure.map(sns.lineplot, coltime, colTPM)
    #     plt.savefig(os.path.join(path_analysis, "div_by_direction.png"), dpi=300)
    #     plt.close()
    # except:
    #     pass
# %%     
if __name__ == "__main__":
    log2fc_thres = 1.5
    padj_thres = 0.05

    colsample = "ID"
    colTPM = "TPM"
    colsev = "Severity_Binary"
    coltime = "Time"
    nareplace = "Convalescent"

    path_exp = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/4_expmtx/ConfirmedRecoveredVaccine/expression_matrix_genes.results_TPM.tsv"
    path_meta = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Resources/Data/Infectomics_CRF_severity_only_20230410.tsv"

    list_upreg_genes = []
    list_downreg_genes = []
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
    list_compare = ["Visit2_Severe__Visit2_Mild"]
    # list_compare = ["Visit1_Severe__Visit1_Mild", "Visit2_Severe__Visit2_Mild", "Visit3_Severe__Visit3_Mild", "Visit4_Severe__Visit4_Mild"]
    for compare in list_compare:
        path_deg = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/5_deg/{compare}_20231007.tsv"
        outdir = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/7_compare_severity/{colsev}_{str(log2fc_thres).replace('.','')}_{str(padj_thres).replace('.','')}/{compare}"
        os.makedirs(outdir, exist_ok=True)
# %%
        main(path_exp, path_deg, path_meta, outdir, list_drop_samples, list_upreg_genes, list_downreg_genes, log2fc_thres, padj_thres, colsample, colTPM, colsev, coltime, nareplace)


# %%

direction = "hypo"
path_marker = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231020/methyl_{direction}_severe_mild_annot_genelist.tsv"
list_marker = list()
with open(path_marker, mode="r") as fr:
    _skiprow = fr.readline()
    for line in fr:
        record = line.rstrip("\n").split("\t")
        marker = record[0]
        list_marker.append(marker)

df_TPM = read_TPM_matrix(path_exp, list_drop_samples)
df_TPM = select_TPM_matrix(df_TPM, select_pattern="L", delim_id="-", namepos=2)
df_TPM_gene_filtered = filter_list_genes(df_TPM, list_marker, mode="genesymbol")
df_TPM_gene_filtered_drop_dup = remove_duplicate_genes_not_max_exp(df_TPM_gene_filtered)
df_summary = make_summary_table_markers(df_TPM_gene_filtered_drop_dup)
df_merged = merge_metadata(df_summary, path_meta)
df_merged.loc[(df_merged["Severity_Binary"].isna()) & df_merged["ID"].str.contains("C19-R"), "Severity_Binary"] = 'Convalescent'
df_merged.loc[(df_merged["Severity_Binary"].isna()) & df_merged["ID"].str.contains("U10K"), "Severity_Binary"] = 'Healthy'
ax = sns.catplot(kind="point", data=df_merged, x="Time", y="NormTPM", hue="Severity_Binary", errorbar=('ci', 95), order=["0", "1", "2", "3", "4", "5"], palette={"Healthy":"#2E8B57" ,"Severe":"#D68C78", "Mild":"#aaaaaa", "Convalescent":"#75A7C3"}, capsize=.15, errwidth=0.5, legend=False, height=5, aspect=7/5)
list_xticklabels = get_xticklabels(df_merged)[::-1]
ax.set_xticklabels(list_xticklabels, fontsize=8)
figpathpoint = os.path.join(outdir, "severity_NormTPM_downreg_pointplot.png")
plt.legend(loc = "upper left", bbox_to_anchor = (1.05, 0.6))
plt.xlabel("Visit", fontsize=11)
plt.ylabel("Gene-wise Normalized TPM", fontsize=11)
plt.tight_layout()
plt.show()
plt.savefig(figpathpoint, dpi=300)
plt.close()

# %%
def filter_list_genes(df_TPM, list_gene, mode, colgene="ID"):
    if len(list_gene) > 0:
        if mode == "genesymbol":
            df_TPM["tmp"] = df_TPM[colgene].apply(lambda x: "_".join(x.split("_")[1:]))
            df_TPM_gene_filtered = df_TPM[df_TPM["tmp"].isin(list_gene)].drop(columns=["tmp"])

        if mode == "geneonly":
            df_TPM_gene_filtered = df_TPM[df_TPM[colgene].isin(list_gene)]

        if mode == "noversion":
            df_TPM["tmp"] = df_TPM[colgene].apply(lambda x: x.split(".")[0] + "_" + "_".join(x.split("_")[1:]))
            df_TPM_gene_filtered = df_TPM[df_TPM["tmp"].isin(list_gene)].drop(columns=["tmp"])
        
    else:
        df_TPM_gene_filtered = df_TPM

    return df_TPM_gene_filtered

def remove_duplicate_genes_not_max_exp(df_TPM_gene_filtered):
    df_TPM_gene_filtered["Meanexp"] = df_TPM_gene_filtered.mean(axis=1)
    df_TPM_gene_filtered["Genesymbol"] = df_TPM_gene_filtered["ID"].apply(lambda x: "_".join(x.split("_")[1:]))
    df_TPM_gene_filtered = df_TPM_gene_filtered.reset_index(drop=True)
    dict_gene_mean_exp = dict()
    for ind, (gene, meanexp) in enumerate(zip(df_TPM_gene_filtered["Genesymbol"], df_TPM_gene_filtered["Meanexp"])):
        if dict_gene_mean_exp.get(gene) == None:
            dict_gene_mean_exp[gene] = list()
        dict_gene_mean_exp[gene].append({ind: meanexp})

    list_drop_indices = list()
    for key, values in dict_gene_mean_exp.items():
        if len(values) > 1:
            sorted_values = sorted(values, key=lambda x: list(x.values())[0], reverse=True)
            drop_indices  = list(map(lambda x: list(x.keys())[0], sorted_values[1:]))
            list_drop_indices.extend(drop_indices)

    list_drop_indices = list(map(int, list_drop_indices))
    df_TPM_gene_filtered_drop_dup = df_TPM_gene_filtered.drop(index=list_drop_indices)

    df_TPM_gene_filtered_drop_dup = df_TPM_gene_filtered_drop_dup.drop(columns="Genesymbol")

    df_TPM_gene_filtered_drop_dup = df_TPM_gene_filtered_drop_dup.reset_index(drop=True)

    return df_TPM_gene_filtered_drop_dup


def make_summary_table_markers(df_TPM_gene_filtered_drop_dup):
    list_df = list()
    gene_id = df_TPM_gene_filtered_drop_dup["ID"]
    for ind, gene in enumerate(gene_id):
        series_selected = df_TPM_gene_filtered_drop_dup.copy().iloc[ind, 1:].astype(float)
        rowmean = df_TPM_gene_filtered_drop_dup.copy().iloc[ind, -1]
        series_norm = series_selected/rowmean
        series_norm_mean_exp_remv = series_norm.iloc[:-1]
        df_split_phase = pd.DataFrame(series_norm_mean_exp_remv).reset_index(drop=False)
        df_split_phase.columns = ["ID", "TPM"]
        df_split_phase[coltime] = df_split_phase[colsample].apply(get_visit_order)
        df_split_phase = df_split_phase.sort_values(by=coltime, ascending=False)
        df_split_phase[coltime] = df_split_phase[coltime].astype(str)
        df_split_phase[colsample] = df_split_phase[colsample].apply(lambda x:"-".join(x.split("-")[:-1]))
        df_split_phase["Gene"] = gene
        list_df.append(df_split_phase)
    df_summary = pd.concat(list_df, axis=0)

    return df_summary


def merge_metadata(df_summary, path_meta, colsample="ID", colsev="Severity_Binary"):
    df_meta = pd.read_csv(path_meta, sep="\t")
    df_meta = df_meta.rename(columns={"SampleID": colsample})
    df_merged = pd.merge(df_summary, df_meta, how="outer", on=colsample)

    return df_merged


