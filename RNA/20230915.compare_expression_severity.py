 # %%
import copy
import math
import os
import sys
from functools import reduce

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


# %%
def main(path_exp, path_deg, path_crf, path_analysis, list_drop_samples, list_upreg_genes, list_downreg_genes, trial_flag, dict_trial, log2fc_thres, padj_thres):
    dict_deg = get_dict_degs(path_deg, log2fc_thres, padj_thres)   

    df_TPM = read_TPM_matrix(path_exp, list_drop_samples)

    list_colname = list(df_TPM.columns)
    list_colname_filt = list(filter(lambda x: x.split("-")[1][0] != "R" if len(x.split("-")) > 1 else x, list_colname))
    df_TPM = df_TPM.loc[:, list_colname_filt]

    list_df_split_phase = list()
    for gene, deg_info in dict_deg.items():
        if len(gene.split(".")) > 1:
            gene = gene.split(".")[0] + "_"  + gene.split("_")[-1]
        df_selected = select_gene_of_interest(df_TPM, gene)
        try:
            df_split_phase = split_sample_by_phase(df_selected)
        except:
            print(f"{gene} not found!")
        genesym = get_gene_info(gene, deg_info)
        df_split_phase_modified = add_gene_column_dataframe(df_split_phase, genesym)
        df_split_phase_modified = add_effect_direction_column_dataframe(df_split_phase_modified, deg_info)
        if trial_flag == True:
            df_split_phase_modified = add_trial_column_dataframe(df_split_phase_modified, dict_trial)
        list_df_split_phase.append(df_split_phase_modified)

    df_all_sample_split_phase = pd.concat(list_df_split_phase, axis=0)
    df_save = save_dataframe(df_all_sample_split_phase, path_analysis)
    report_statistics(df_save, path_analysis)

    path_deg_res = os.path.join(path_analysis, "deg_results.txt")
    df_analysis = pd.read_csv(path_deg_res, sep="\t")
    df_crf = pd.read_csv(path_crf, sep="\t")
    df_merged = pd.merge(df_analysis, df_crf, how="outer", on="SampleID") 
    df_merged["Severity"] = df_merged["Severity"].fillna("Healthy")
    df_merged_downreg = df_merged[df_merged["Direction"]=="negative"]
    if len(list_downreg_genes) > 0:
        df_merged_downreg = df_merged_downreg[df_merged_downreg["Gene"].isin(list_downreg_genes)]
    df_merged_upreg = df_merged[df_merged["Direction"]=="positive"]
    if len(list_upreg_genes) > 0:
        df_merged_upreg = df_merged_upreg[df_merged_upreg["Gene"].isin(list_upreg_genes)]
    
    plot_(dict_deg, df_merged_downreg, path_analysis, "downreg")
    plot_(dict_deg, df_merged_upreg, path_analysis, "upreg")
    plot_total_line_(df_merged_downreg, path_analysis, "downreg")
    plot_total_line_(df_merged_upreg, path_analysis, "upreg")
    plot_total_box_(df_merged_downreg, path_analysis, "downreg")
    plot_total_box_(df_merged_upreg, path_analysis, "upreg")
    # plot_trend(df_merged_downreg, dict_deg, path_analysis, tag="downreg")
    # plot_trend(df_merged_upreg, dict_deg, path_analysis, tag="upreg")


def read_TPM_matrix(path_exp, list_drop_samples):
    df_TPM = pd.read_csv(path_exp, sep="\t")
    df_TPM  = df_TPM.drop(columns=list_drop_samples)

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

def select_gene_of_interest(df_TPM, target_gene, colgene="ID"):
    df_selected = df_TPM[df_TPM[colgene] == target_gene]
    df_selected = df_selected.set_index(colgene)
    df_selected.index = df_selected.index.str.split("_").str[1]
    df_selected = df_selected.T

    return df_selected

def split_sample_by_phase(df_selected, colsample="SampleID", colTPM="TPM", coltime="Time"):
    df_selected_reidx = df_selected.reset_index(drop=False)
    df_selected_reidx.columns = [colsample, colTPM]
    df_split_phase = df_selected_reidx.copy()
    df_split_phase[coltime] = df_split_phase[colsample].apply(lambda x: x.replace("L1", "V5").replace("A1", "V7"))
    df_split_phase[coltime] = df_split_phase[coltime].apply(lambda x: x.split("-")[-1])
    df_split_phase[coltime] = df_split_phase[coltime].apply(lambda x: int(x[1:]))
    df_split_phase[colsample] = df_split_phase[colsample].apply(lambda x:x.split("-")[0] + "-" + x.split("-")[1])
    df_split_phase = df_split_phase.sort_values(by=coltime, ascending=True)
    df_split_phase[coltime] = df_split_phase[coltime].astype(str)
    # df_split_phase = df_split_phase.drop(columns=["tmp"])

    return df_split_phase

def get_gene_info(gene, deg_info):
    genesymbol = gene.split("_")[-1]
    log2fc = round(float(deg_info["log2fc"]), 3)
    padj = round(float(deg_info["padj"]), 3)
    gene_info = f"{genesymbol}\nlog2FC = {log2fc}\nPadj = {padj}"

    return gene_info

def add_gene_column_dataframe(df, genesym, colgene="Gene"):
    df[colgene] = genesym

    return df

def add_effect_direction_column_dataframe(df, deg_info, coleffect="Effect", coldirec="Direction"):
    log2fc = float(deg_info["log2fc"])
    df[coleffect] = log2fc
    if float(log2fc) < 0:
        direction = "negative"
    else:
        direction = "positive"
    df[coldirec] = direction
    
    return df

def add_trial_column_dataframe(df, dict_trial, coltime="Time", coltrial="Trial"):
    list_trials = list()
    for time in df[coltime].to_list():
        if isinstance(time, int):
           time = str(time) 
        for trial, list_times in dict_trial.items():
            if time in list_times:
                list_trials.append(trial)

    df[coltrial] = list_trials

    return df

def save_dataframe(df, path_analysis):
    df_save = copy.deepcopy(df)
    df_save["Gene"] = df_save["Gene"].apply(lambda x: x.split("\n")[0])
    df_save.to_csv(os.path.join(path_analysis, "deg_results.txt"), sep="\t", index=False)
    
    return df_save

def report_statistics(df_save, path_analysis):
    downreg_genes = list(df_save[df_save["Direction"]=="negative"]["Gene"].unique())
    upreg_genes = list(df_save[df_save["Direction"]=="positive"]["Gene"].unique())
    dict_stat = {"downreg":{"number":len(downreg_genes), "genes":downreg_genes}, "upreg":{"number":len(upreg_genes), "genes":upreg_genes}}
    df_stat = pd.DataFrame.from_dict(dict_stat, orient='index')

    df_stat.to_csv(os.path.join(path_analysis, "deg_stats.txt"), sep="\t")

    return df_stat

# def get_pairs(df, x, y):
#     from itertools import product
#     list_x = list(set(df[x].to_list()))
#     list_y = list(set(df[y].to_list()))
#     list_pair = list(product(list_x, list_y))
#     list_chunk_pairs = [list_pair[x: x+len(list_y)] for x in range(0, len(list_pair), len(list_y))]

#     return list_chunk_pairs

# def stat_annotation(ax, **kwargs):
#     pairs = get_pairs(kwargs["data"], kwargs["x"], kwargs["hue"])
#     annotator = Annotator(ax, pairs, **kwargs)
#     annotator.configure(test='Mann-Whitney', text_format='star', loc='outside')
#     annotator.apply_and_annotate()

def plot_(dict_deg, df_merged, path_analysis, tag):
    list_deg_down = list(filter(lambda x : float(dict_deg[x]["log2fc"])<0, dict_deg.keys()))
    list_deg_up = list(filter(lambda x : float(dict_deg[x]["log2fc"])>0, dict_deg.keys()))
    sorted_deg_down = sorted(list_deg_down, key = lambda x : float(dict_deg[x]["log2fc"]))
    sorted_deg_up = sorted(list_deg_up, key = lambda x : float(dict_deg[x]["log2fc"]), reverse = True)

    import math
    num_col = int(math.sqrt(len(dict_deg))*0.7)

    fig, ax = plt.subplots(ncols = num_col, nrows = math.ceil(len(sorted_deg_up)/num_col), figsize = (15, 40))
    axes = ax.flatten()

    severity_order = sorted(df_merged["Severity"].unique())
    hue_list = ["#CD313A", "#0047A0"]

    for ind, geneid in enumerate(sorted_deg_up):
        geneid = geneid.split('_')[-1]
        table_part = df_merged[df_merged["Gene"] == geneid].copy()
        sns.violinplot(data = table_part, x = "Time", y = "TPM", hue = "Severity", hue_order=severity_order, palette=hue_list, ax = axes[ind], legend = False)
        axes[ind].set_title(geneid, fontsize = 10)
        axes[ind].set_ylim(bottom = 0)
        if ind == num_col-1:
            axes[ind].get_legend().set_visible(True)
            axes[ind].legend(loc = "upper left", bbox_to_anchor = (1.05, 0.95))
        else:
            axes[ind].get_legend().set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(path_analysis, f"severity_TPM_each_gene_across_time_{tag}.png"), dpi=300)
    plt.close()

def plot_total_line_(df_merged, path_analysis, tag):
    df_mean_tpm_severity = df_merged.groupby(["Severity", "Time", "SampleID"])["TPM"].sum().reset_index(drop=False)
    df_num = df_mean_tpm_severity.groupby(["Severity", "Time"])["SampleID"].count().reset_index(drop=False)
    list_cnt = df_num["SampleID"].to_list()
    ax = sns.pointplot(data=df_mean_tpm_severity, x="Time", y="TPM", hue="Severity", errwidth=0.5)
    plt.yscale("log")
    xticklabel = [f"1\n(1):{list_cnt[0]}\n(2):{list_cnt[5]}\n(3):{list_cnt[10]}\n(4):{list_cnt[14]}",
                f"2\n(1):{list_cnt[1]}\n(2):{list_cnt[6]}\n(3):{list_cnt[11]}\n(4):{list_cnt[15]}",
                f"3\n(1):{list_cnt[2]}\n(2):{list_cnt[7]}\n(3):{list_cnt[12]}\n(4):{list_cnt[16]}",
                f"4\n(1):{list_cnt[3]}\n(2):{list_cnt[8]}\n(3):{list_cnt[13]}\n(4):{list_cnt[17]}",
                f"L\n(1):{list_cnt[4]}\n(2):{list_cnt[9]}\n(3):0\n(4):0",
                f"HC\nHealthy:{list_cnt[18]}"]
    ax.set_xticklabels(xticklabel, fontsize=7)
    ax.set_xlabel("TimePoints", fontsize=12)
    ax.set_ylabel("$log10$(MeanTPM)", fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(path_analysis, f"severity_mean_TPM_across_time_line_{tag}.png"), dpi=300)
    plt.close()

def plot_total_box_(df_merged, path_analysis, tag):
    df_mean_tpm_severity = df_merged.groupby(["Severity", "Time", "SampleID"])["TPM"].sum().reset_index(drop=False)
    df_num = df_mean_tpm_severity.groupby(["Severity", "Time"])["SampleID"].count().reset_index(drop=False)
    list_cnt = df_num["SampleID"].to_list()
    ax = sns.boxplot(data=df_mean_tpm_severity, x="Time", y="TPM", hue="Severity")
    plt.yscale("log")
    xticklabel = [f"1\n(1):{list_cnt[0]}\n(2):{list_cnt[5]}\n(3):{list_cnt[10]}\n(4):{list_cnt[14]}",
                f"2\n(1):{list_cnt[1]}\n(2):{list_cnt[6]}\n(3):{list_cnt[11]}\n(4):{list_cnt[15]}",
                f"3\n(1):{list_cnt[2]}\n(2):{list_cnt[7]}\n(3):{list_cnt[12]}\n(4):{list_cnt[16]}",
                f"4\n(1):{list_cnt[3]}\n(2):{list_cnt[8]}\n(3):{list_cnt[13]}\n(4):{list_cnt[17]}",
                f"L\n(1):{list_cnt[4]}\n(2):{list_cnt[9]}\n(3):0\n(4):0",
                f"HC\nHealthy:{list_cnt[18]}"]
    ax.set_xticklabels(xticklabel, fontsize=10)
    ax.set_xlabel("TimePoints", fontsize=12)
    ax.set_ylabel("$log10$(MeanTPM)", fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(path_analysis, f"severity_mean_TPM_across_time_{tag}.png"), dpi=300)
    plt.close()

# def plot_trend(df, dict_deg, path_analysis, colgene="Gene", colTPM="TPM", coltime="Time", coltrial="Trial", coldirec="Direction", split_by="Severity", tag="downreg"):
#     # https://github.com/webermarcolivier/statannot/issues/27
#     # https://github.com/mwaskom/seaborn/issues/1166
#     order = df[coltime].unique()
#     hue_order = df[split_by].unique()
#     palette = "Set2"
#     for ind, geneid in enumerate(sorted_deg_up):
#         geneid = geneid.split('_')[-1]
#         table_part = df_merged[df_merged["Gene"] == geneid].copy()
#         sns.violinplot(data = table_part, x = "Time", y = "TPM", hue = "Severity", hue_order=severity_order, palette=hue_list, ax = axes[ind], legend = False)
#         axes[ind].set_title(geneid, fontsize = 10)
#         axes[ind].set_ylim(bottom = 0)
#         if ind == num_col-1:
#             axes[ind].get_legend().set_visible(True)
#             axes[ind].legend(loc = "upper left", bbox_to_anchor = (1.05, 0.95))
#         else:
#             axes[ind].get_legend().set_visible(False)
#     plt.tight_layout()
#     plt.savefig(os.path.join(path_analysis, f"severity_TPM_each_gene_across_time_{tag}.png"), dpi=300)
#     plt.close()
    
    # figure = sns.FacetGrid(data=df, col=colgene, col_wrap=int(math.sqrt(len(dict_deg))), hue=split_by, sharex=False, sharey=False, legend_out=True)
    # figure.map(sns.violinplot, coltime, colTPM, split_by, palette=palette, hue_order=hue_order, order=order)
    # labels = hue_order
    # colors = sns.color_palette(palette).as_hex()[:len(labels)]
    # handles = [patches.Patch(color=col, label=lab) for col, lab in zip(colors, labels)]
    # plt.legend(handles=handles, title="Severity", loc="center left", bbox_to_anchor=(1, 0.75))
    # plt.savefig(os.path.join(path_analysis, f"default_split_by_{split_by}_{tag}.png"), dpi=300)
    # plt.close()

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
    # list_compare = ["Visit1_Severe__Visit4_Severe", "Visit1_Mild__Visit4_Mild"]
    list_compare = ["Visit1_Severe__Visit1_Mild", "Visit2_Severe__Visit2_Mild", "Visit3_Severe__Visit3_Mild", "Visit4_Severe__Visit4_Mild"]
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
    list_upreg_genes = []
    list_downreg_genes = []
    dict_trial = dict()
    trial_flag = False
    log2fc_thres = 2.0
    padj_thres = 0.05
    for compare in list_compare:
        path_exp = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/4_expmtx/Infectomics/expression_matrix_genes.results_TPM.tsv"
        path_deg = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/5_deg/{compare}_20231007.tsv"
        path_crf = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Resources/Data/COVID19_master_table_added_CRF_20231007.txt"
        path_analysis = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/7_compare_severity/severity_20_005/{compare}"
        os.makedirs(path_analysis, exist_ok=True)

        main(path_exp, path_deg, path_crf, path_analysis, list_drop_samples, list_upreg_genes, list_downreg_genes, trial_flag, dict_trial, log2fc_thres, padj_thres)
