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
def main(path_exp, path_deg, path_crf, path_analysis, list_drop_samples, trial_flag, dict_trial, log2fc_thres, padj_thres):
    dict_deg = get_dict_degs(path_deg, log2fc_thres, padj_thres)   

    df_TPM = read_TPM_matrix(path_exp, list_drop_samples)

    list_df_split_phase = list()
    for gene, deg_info in dict_deg.items():
        df_selected = select_gene_of_interest(df_TPM, gene)
        df_split_phase = split_sample_by_phase(df_selected)
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
    df_merged = pd.merge(df_analysis, df_crf, how="inner", on="SampleID") 

    plot_trend(df_merged, dict_deg, path_analysis)


def read_TPM_matrix(path_exp, list_drop_samples):
    df_TPM = pd.read_csv(path_exp, sep="\t")
    if len(list_drop_samples) > 0:
        df_TPM = df_TPM.drop(columns=list_drop_samples)
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
    df_split_phase = df_selected.reset_index()
    df_split_phase.columns = [colsample, colTPM]
    df_split_phase[coltime] = df_split_phase[colsample].apply(lambda x:x.split("-")[-1])
    df_split_phase[coltime] = df_split_phase[coltime].apply(lambda x: x.replace("L1", "V5"))
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

    df_stat.to_csv(os.path.join(path_analysis, "deg_stats.txt"), sep="\t", index=False)

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

def plot_trend(df, dict_deg, path_analysis, colgene="Gene", colTPM="TPM", coltime="Time", coltrial="Trial", coldirec="Direction", split_by="Severity_Binary"):
    hue_order=df[split_by].unique()
    palette="Set1"
    # https://github.com/webermarcolivier/statannot/issues/27
    figure = sns.FacetGrid(data=df, col=colgene, col_wrap=int(math.sqrt(len(dict_deg))), hue=split_by, sharex=False, sharey=False, legend_out=True).add_legend()
    figure.map(sns.violinplot, coltime, colTPM, split_by, hue_order=hue_order, order=df[coltime].unique(), palette=palette)
    # figure.map(sns.lineplot, coltime, colTPM)
    labels = hue_order
    colors = sns.color_palette(palette).as_hex()[:len(labels)]
    handles = [patches.Patch(color=col, label=lab) for col, lab in zip(colors, labels)]
    plt.legend(handles=handles, title="Severity", loc="center left", bbox_to_anchor=(1, 0.75))
    plt.savefig(os.path.join(path_analysis, f"default_split_by_{split_by}.png"), dpi=300)
    plt.close()

    try:
        figure = sns.FacetGrid(data=df, col=colgene, hue=coltrial, col_wrap=int(math.sqrt(len(dict_deg))), sharex=False, sharey=False, legend_out=True)
        figure.map(sns.lineplot, coltime, colTPM)
        plt.savefig(os.path.join(path_analysis, "div_by_trial.png"), dpi=300)
        plt.close()
    except:
        pass

    try:
        figure = sns.FacetGrid(data=df, col=colgene, hue=coldirec, col_wrap=int(math.sqrt(len(dict_deg))), sharex=False, sharey=False, legend_out=True)
        figure.map(sns.violinplot, coltime, colTPM)
        figure.map(sns.lineplot, coltime, colTPM)
        plt.savefig(os.path.join(path_analysis, "div_by_direction.png"), dpi=300)
        plt.close()
    except:
        pass
# %%     
if __name__ == "__main__":
    list_visit_compare = ["Visit2__Visit1", "Visit3__Visit1", "Visit4__Visit1", "Visit3__Visit2", "Visit4__Visit2", "Visit4__Visit3"]
    for visit_compare in list_visit_compare:
        path_exp = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/4_expmtx/Confirmed/expression_matrix_genes.results_TPM.tsv"
        path_deg = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/5_deg/{visit_compare}_DEG_20230925.tsv"
        path_crf = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Resources/Data/Infectomics_CRF_20230410.tsv"
        path_analysis = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/6_time_series/paired_20_005/{visit_compare}"
        os.makedirs(path_analysis, exist_ok=True)
        list_drop_samples = []
        dict_trial = dict()
        trial_flag = False
        log2fc_thres = 2.0
        padj_thres = 0.05
        main(path_exp, path_deg, path_crf, path_analysis, list_drop_samples, trial_flag, dict_trial, log2fc_thres, padj_thres)

# %%
# if __name__ == "__main__":
#     list_visit_compare = ["Visit2__Visit1", "Visit4__Visit1", "Visit5__Visit1", "Visit6__Visit1", "Visit5__Visit4", "Visit6__Visit4"]
#     for visit_compare in list_visit_compare:
#         path_exp = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Results/20230813/4_expmtx/expression_matrix_genes.results_TPM.tsv"
#         path_deg = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Results/20230902/5_deg/unpaired/{visit_compare}_DEG_20230902.tsv"
#         path_analysis = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Results/20230902/6_time_series/unpaired_05_01/{visit_compare}"
#         os.makedirs(path_analysis, exist_ok=True)
#         list_drop_samples = ["KU10K-00922-V11", "KU10K-00922-V12"]
#         trial_flag = True
#         dict_trial = {"trial1":["1", "2", "3"], "trial2":["4", "5", "6"], "trial3":["7", "8", "9","10"]}
#         log2fc_thres = 0.5
#         padj_thres = 0.1

#         main(path_exp, path_deg, path_analysis, list_drop_samples, trial_flag, dict_trial, log2fc_thres, padj_thres)

# %%
# if __name__ == "__main__":
#     list_visit_compare = ["Visit7__Visit4", "Visit9__Visit4", "Visit10__Visit4", "Visit10__Visit1"]
#     for visit_compare in list_visit_compare:
#         path_exp = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Results/20230813/4_expmtx/expression_matrix_genes.results_TPM.tsv"
#         path_deg = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Results/20230904/5_deg/unpaired_V789/{visit_compare}_DEG_20230902.tsv"
#         path_analysis = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Results/20230904/6_time_series/unpaired_V789/{visit_compare}"
#         os.makedirs(path_analysis, exist_ok=True)
#         list_drop_samples = ["KU10K-00922-V11", "KU10K-00922-V12"]
#         trial_flag = True
#         dict_trial = {"trial1":["1", "2", "3"], "trial2":["4", "5", "6"], "trial3":["7", "8", "9","10"]}
#         log2fc_thres = 0.5
#         padj_thres = 0.1
#         main(path_exp, path_deg, path_analysis, list_drop_samples, trial_flag, dict_trial, log2fc_thres, padj_thres)
# %%
# def get_target_files_in_directory(path_target,pattern):
#     list_target_files = []
#     list_files = os.listdir(path_target)
#     for file in list_files:
#         path_target_file = os.path.join(path_target,file)
#         if os.path.isdir(path_target_file):
#             list_traget_files = os.listdir(path_target_file)
#             for files in list_traget_files:
#                 if pattern in files:
#                     path_target_file = os.path.join(path_target_file,files)
#                     list_target_files.append(path_target_file)
#         else:
#             list_target_files.append(path_target_file)
#     return list_target_files

# def merge_all_expression_values(list_rsem_expr_file,which_measure="TPM"):
#     list_df_expr = []
#     for file_expr in list_rsem_expr_file:
#         sample_id = get_sample_id(file_expr)
#         df_expr = pd.read_csv(file_expr,sep="\t",usecols=["gene_id",which_measure])
#         df_expr.columns = ["gene_id",sample_id]
#         df_expr["gene_id"] = df_expr["gene_id"].apply(lambda x:x.split("_")[-1])
#         list_df_expr.append(df_expr)
#     df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['gene_id'],how='outer'), list_df_expr)
#     return df_merged

# %%
# %%
# list_deg_down = list(filter(lambda x : float(dict_deg[x]["log2fc"])<0, dict_deg.keys()))
# list_deg_up = list(filter(lambda x : float(dict_deg[x]["log2fc"])>0, dict_deg.keys()))
# sorted_deg_down = sorted(list_deg_down, key = lambda x : float(dict_deg[x]["log2fc"]))
# sorted_deg_up = sorted(list_deg_up, key = lambda x : float(dict_deg[x]["log2fc"]), reverse = True)

# import math

# num_col = int(math.sqrt(len(dict_deg))*0.7)

# fig, ax = plt.subplots(ncols = num_col, nrows = math.ceil(len(sorted_deg_up)/num_col), figsize = (15, 40))
# axes = ax.flatten()

# severity_order = sorted(df_merged["Severity_Binary"].unique())
# hue_list = ["#CD313A", "#0047A0"]

# for ind, geneid in enumerate(sorted_deg_up):
#     geneid = geneid.split('_')[-1]
#     table_part = df_merged[df_merged["Gene"] == geneid].copy()
#     sns.violinplot(data = table_part, x = "Time", y = "TPM", hue = "Severity_Binary", hue_order=severity_order, palette=hue_list, ax = axes[ind], legend = False)
#     axes[ind].set_title(geneid, fontsize = 10)
#     axes[ind].set_ylim(bottom = 0)
#     if ind == num_col-1:
#         axes[ind].get_legend().set_visible(True)
#         axes[ind].legend(loc = "upper left", bbox_to_anchor = (1.05, 0.95))
#     else:
#         axes[ind].get_legend().set_visible(False)
# plt.tight_layout()
# plt.show()