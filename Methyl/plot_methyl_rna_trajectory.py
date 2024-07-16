#%%
import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import gridspec
import seaborn as sns
from matplotlib.cm import ScalarMappable
from matplotlib import patches as mpatches
from scipy.stats import ttest_rel, ranksums, wilcoxon
from statannotations.Annotator import Annotator

import warnings
warnings.filterwarnings("ignore")

import math, os

SAMPLE_CONFIG = {
    "col_id" : "Project_ID_Alias",
    "col_subj" : "Project_ID",
    "col_severity" : "Severity",
    "col_trait" : "Visit_order",
    "trait_info" : {
        "Longitudinal" : {
            "First" : ["Visit1"],
            "Last" : ["Visit3", "Visit4"]
        },
        "Single" : ["Recover"]
    }
}
PLOT_CONFIG = {
    "Longitudinal": {
        1 : dict(facecolor = 'cornflowerblue', headwidth = 15, width = 3, alpha = 1),
        2 : dict(facecolor = 'cornflowerblue', headwidth = 15, width = 3, alpha = 1),
        3 : dict(facecolor = 'firebrick', headwidth = 15, width = 3, alpha = 1),
        4 : dict(facecolor = 'firebrick', headwidth = 15, width = 3, alpha = 1)
    },
    "Single": {
        1 : 'dimgray',
        2 : 'dimgray',
        3 : 'dimgray',
        4 : 'dimgray'
    }
}
SEV_CONFIG = {
    1 : "Mild",
    2 : "Mild",
    3 : "Severe",
    4 : "Severe"
}

def main(path_methyl_rna_corr, path_cpg, path_exp, path_meta, path_save_png_format, path_save_visit_relation_pval):
    table_methyl_rna_corr = pd.read_csv(path_methyl_rna_corr, sep = '\t')
    table_cpg = pd.read_csv(path_cpg, sep = '\t')
    table_exp = pd.read_csv(path_exp, sep = '\t')
    table_meta = pd.read_csv(path_meta, sep = '\t')
    
    dict_sample_info = organize_sample_from_meta_table(table_meta)
    dict_sample_to_sev = get_sample_to_severity_relation(table_meta)
    dict_sample_to_cpg = get_sample_to_feat(table_cpg, ["chr", "start", "end"])
    dict_sample_to_exp = get_sample_to_feat(table_exp, ["index"])
    
    table_visit_compare = plot_rna_methyl_trajectory(table_methyl_rna_corr, dict_sample_info, dict_sample_to_cpg, dict_sample_to_exp, dict_sample_to_sev, path_save_png_format)
    table_visit_compare.to_csv(path_save_visit_relation_pval, sep = '\t', index = False)
    return table_visit_compare
    
def organize_sample_from_meta_table(table_meta, n_fu = 2):
    dict_sample_info = dict()
    
    from itertools import chain
    trait_values_longitudinal = list(chain(*list(SAMPLE_CONFIG["trait_info"]["Longitudinal"].values())))
    trait_values_single = SAMPLE_CONFIG["trait_info"]["Single"]
    table_meta_longitudinal = table_meta[table_meta[SAMPLE_CONFIG["col_trait"]].isin(trait_values_longitudinal)]
    table_meta_single = table_meta[table_meta[SAMPLE_CONFIG["col_trait"]].isin(trait_values_single)]
    
    dict_sample_info["Single"] = table_meta_single[SAMPLE_CONFIG["col_id"]].to_list()
    dict_sample_info["Longitudinal"] = list()
    
    for subjid in table_meta_longitudinal[SAMPLE_CONFIG["col_subj"]].unique():
        table_meta_single_subj = table_meta_longitudinal[table_meta_longitudinal[SAMPLE_CONFIG["col_subj"]] == subjid]
        rows_first = table_meta_single_subj[table_meta_single_subj[SAMPLE_CONFIG["col_trait"]].isin(SAMPLE_CONFIG["trait_info"]["Longitudinal"]["First"])]
        rows_last = table_meta_single_subj[table_meta_single_subj[SAMPLE_CONFIG["col_trait"]].isin(SAMPLE_CONFIG["trait_info"]["Longitudinal"]["Last"])]
        
        assert rows_first.shape[0] == 1 and rows_last.shape[0] == 1, table_meta_single_subj
        
        id_first = rows_first[SAMPLE_CONFIG["col_id"]].to_list()[0]
        id_last = rows_last[SAMPLE_CONFIG["col_id"]].to_list()[0]
        
        dict_sample_info["Longitudinal"].append((id_first, id_last))
    return dict_sample_info
        
def get_sample_to_severity_relation(table_meta):
    dict_sample_to_sev = dict(zip(table_meta[SAMPLE_CONFIG["col_id"]], table_meta[SAMPLE_CONFIG["col_severity"]]))
    return dict_sample_to_sev
        
def get_sample_to_feat(table, cols_feat):
    list_sample = list(filter(lambda x : x not in cols_feat, table.columns))
    table["index"] = table.index
    table["featname"] = table[cols_feat].apply(lambda row : '_'.join(list(map(str, row))), axis = 1)
    dict_sample_to_feat = dict()
    for sample in list_sample:
        dict_single_feat = dict(zip(table["featname"], table[sample]))
        dict_sample_to_feat[sample] = dict_single_feat
    return dict_sample_to_feat
        
def plot_rna_methyl_trajectory(table_methyl_rna_corr, dict_sample_info, dict_sample_to_cpg, dict_sample_to_exp, dict_sample_to_sev, path_save_png_format = None):
    from itertools import chain
    
    table_visit_compare = pd.DataFrame(columns = ["Methyl", "RNA", "p_Mild_RNA_T", "p_Mild_RNA_W", "p_Severe_RNA_T", "p_Severe_RNA_W", "p_Mild_Methyl_T", "p_Mild_Methyl_W", "p_Severe_Methyl_T", "p_Severe_Methyl_W", "p_Acute_Methyl", "p_Acute_RNA", "p_recover_Methyl", "p_recover_RNA"])
    for ind, corr_row in table_methyl_rna_corr.iterrows():
        dict_corr = corr_row.to_dict()
        fig = plt.figure(figsize = (11, 9))
        plt.rcParams["font.size"] = 13
        
        gsfig = gridspec.GridSpec(
            100, 123,
            left = 0, right = 1, bottom = 0,
            top = 1, wspace = 1, hspace = 1
        )
        gs1 = gsfig[0:70, 0:25]
        ax1 = fig.add_subplot(gs1)
        
        gs2 = gsfig[0:70, 33:103]
        ax2 = fig.add_subplot(gs2)
        
        gs3 = gsfig[75:100, 33:103]
        ax3 = fig.add_subplot(gs3)
        
        gs4 = gsfig[20:40, 103:123]
        ax4 = fig.add_subplot(gs4)
        # plot Single-Sample
        list_single_sample = dict_sample_info["Single"]
        try:
            list_single_exp = list(map(lambda sample : dict_sample_to_exp[sample][dict_corr["RNA"]], list_single_sample))
            list_single_cpg = list(map(lambda sample : dict_sample_to_cpg[sample][dict_corr["Methyl"]], list_single_sample))
            list_single_color = list(map(lambda sample : PLOT_CONFIG["Single"][dict_sample_to_sev[sample]], list_single_sample))
        except Exception as e:
            print(e)
            continue
        ax2.scatter(list_single_cpg, list_single_exp, c = list_single_color)
        # plot Longitudinal-Sample
        for sample_first, sample_last in dict_sample_info["Longitudinal"]:
            exp_first, exp_last = dict_sample_to_exp[sample_first][dict_corr["RNA"]], dict_sample_to_exp[sample_last][dict_corr["RNA"]]
            cpg_first, cpg_last = dict_sample_to_cpg[sample_first][dict_corr["Methyl"]], dict_sample_to_cpg[sample_last][dict_corr["Methyl"]]
            arrow_config = PLOT_CONFIG["Longitudinal"][dict_sample_to_sev[sample_first]]
            ax2.annotate("", xy = (cpg_last, exp_last), xytext = (cpg_first, exp_first), arrowprops = arrow_config, zorder = 5 if SEV_CONFIG[dict_sample_to_sev[sample_first]] == "Mild" else 10)
        
        list_longitudinal_sample_first = list(map(lambda x : x[0], dict_sample_info["Longitudinal"]))
        list_longitudinal_sample_first_severe = list(filter(lambda x : SEV_CONFIG[dict_sample_to_sev[x]] == "Severe", list_longitudinal_sample_first))
        list_longitudinal_sample_first_mild = list(filter(lambda x : SEV_CONFIG[dict_sample_to_sev[x]] == "Mild", list_longitudinal_sample_first))
        list_longitudinal_sample_last = list(map(lambda x : x[1], dict_sample_info["Longitudinal"]))
        list_longitudinal_sample_last_severe = list(filter(lambda x : SEV_CONFIG[dict_sample_to_sev[x]] == "Severe", list_longitudinal_sample_last))
        list_longitudinal_sample_last_mild = list(filter(lambda x : SEV_CONFIG[dict_sample_to_sev[x]] == "Mild", list_longitudinal_sample_last))
        list_first_mild_exp = list(map(lambda sample : dict_sample_to_exp[sample][dict_corr["RNA"]], list_longitudinal_sample_first_mild))
        list_first_mild_cpg = list(map(lambda sample : dict_sample_to_cpg[sample][dict_corr["Methyl"]], list_longitudinal_sample_first_mild))
        list_last_mild_exp = list(map(lambda sample : dict_sample_to_exp[sample][dict_corr["RNA"]], list_longitudinal_sample_last_mild))
        list_last_mild_cpg = list(map(lambda sample : dict_sample_to_cpg[sample][dict_corr["Methyl"]], list_longitudinal_sample_last_mild))
        list_first_severe_exp = list(map(lambda sample : dict_sample_to_exp[sample][dict_corr["RNA"]], list_longitudinal_sample_first_severe))
        list_first_severe_cpg = list(map(lambda sample : dict_sample_to_cpg[sample][dict_corr["Methyl"]], list_longitudinal_sample_first_severe))
        list_last_severe_exp = list(map(lambda sample : dict_sample_to_exp[sample][dict_corr["RNA"]], list_longitudinal_sample_last_severe))
        list_last_severe_cpg = list(map(lambda sample : dict_sample_to_cpg[sample][dict_corr["Methyl"]], list_longitudinal_sample_last_severe))
        
        dict_sevname_to_level = dict(zip(SEV_CONFIG.values(), SEV_CONFIG.keys()))
        dict_palette = {
            "Conv." : PLOT_CONFIG["Single"][dict_sevname_to_level["Mild"]],
            "Mild First" : PLOT_CONFIG["Longitudinal"][dict_sevname_to_level["Mild"]]["facecolor"],
            "Mild Last" : PLOT_CONFIG["Longitudinal"][dict_sevname_to_level["Mild"]]["facecolor"],
            "Sev. First" : PLOT_CONFIG["Longitudinal"][dict_sevname_to_level["Severe"]]["facecolor"],
            "Sev. Last" : PLOT_CONFIG["Longitudinal"][dict_sevname_to_level["Severe"]]["facecolor"]
        }
        boxplot_order = ["Conv.", "Mild First", "Mild Last", "Sev. First", "Sev. Last"]
        # plot rna boxplot
        df_rna = pd.DataFrame({
            "Group":["Conv."] * len(list_single_exp) + ["Mild First"] * len(list_first_mild_exp) + ["Mild Last"] * len(list_last_mild_exp) + ["Sev. First"] * len(list_first_severe_exp) + ["Sev. Last"] * len(list_last_severe_exp),
            "Exp" : list_single_exp + list_first_mild_exp + list_last_mild_exp + list_first_severe_exp + list_last_severe_exp
        })
        sns.boxplot(data = df_rna, x = "Group", y = "Exp", palette = dict_palette, ax = ax1, width = 0.55, order = boxplot_order)
        ax1.set_xticks([0,1,2,3,4], ["Conv.", "Mild\nFirst", "Mild\nLast", "Sev.\nFirst", "Sev.\nLast"])   
        for interval in [0.5, 2.5]:
            ax1.axvline(interval, linestyle = '--', linewidth = 2, c = 'gray')
        _, pval_mild_exp_t = ttest_rel(list_first_mild_exp, list_last_mild_exp)
        _, pval_severe_exp_t = ttest_rel(list_first_severe_exp, list_last_severe_exp)
        _, pval_mild_exp_w = wilcoxon(list_first_mild_exp, list_last_mild_exp)
        _, pval_severe_exp_w = wilcoxon(list_first_severe_exp, list_last_severe_exp)
        _, pval_acute_exp = ranksums(list_first_mild_exp, list_first_severe_exp)
        _, pval_recover_exp = ranksums(list_last_mild_exp, list_last_severe_exp)
        annot = Annotator(ax1, pairs = [("Mild First", "Mild Last"), ("Sev. First", "Sev. Last"), ("Mild First", "Sev. First"), ("Mild Last", "Sev. Last")], data = df_rna, x = "Group", y = "Exp", order = boxplot_order)
        annot.configure(test = None, loc = "outside", verbose = 0)
        annot.set_pvalues([pval_mild_exp_t, pval_severe_exp_t, pval_acute_exp, pval_recover_exp])
        annot.annotate()        
        
        # plot methyl boxplot
        df_methyl = pd.DataFrame({
            "Group":["Conv."] * len(list_single_cpg) + ["Mild First"] * len(list_first_mild_cpg) + ["Mild Last"] * len(list_last_mild_cpg) + ["Sev. First"] * len(list_first_severe_cpg) + ["Sev. Last"] * len(list_last_severe_cpg),
            "CpG" : list_single_cpg + list_first_mild_cpg + list_last_mild_cpg + list_first_severe_cpg + list_last_severe_cpg
        })
        sns.boxplot(data = df_methyl, y = "Group", x = "CpG", palette = dict_palette, ax = ax3, width = 0.55, order = boxplot_order)
        ax3.set_yticks([0,1,2,3,4], list(["Conv.", "Mild First", "Mild Last", "Sev. First", "Sev. Last"]))
        for interval in [0.5, 2.5]:
            ax3.axhline(interval, linestyle = '--', linewidth = 2, c = 'gray')
        _, pval_mild_cpg_t = ttest_rel(list_first_mild_cpg, list_last_mild_cpg)
        _, pval_severe_cpg_t = ttest_rel(list_first_severe_cpg, list_last_severe_cpg)
        _, pval_mild_cpg_w = wilcoxon(list_first_mild_cpg, list_last_mild_cpg)
        _, pval_severe_cpg_w = wilcoxon(list_first_severe_cpg, list_last_severe_cpg)
        _, pval_acute_cpg = ranksums(list_first_mild_cpg, list_first_severe_cpg)
        _, pval_recover_cpg = ranksums(list_last_mild_cpg, list_last_severe_cpg)
        annot = Annotator(ax3, pairs = [("Mild First", "Mild Last"), ("Sev. First", "Sev. Last"), ("Mild First", "Sev. First"), ("Mild Last", "Sev. Last")], data = df_methyl, y = "Group", x = "CpG", order = boxplot_order, orient = 'h')
        annot.configure(test = None, loc = "outside", verbose = 0)
        annot.set_pvalues([pval_mild_cpg_t, pval_severe_cpg_t, pval_acute_cpg, pval_recover_cpg])
        annot.annotate()        
        
        table_visit_compare.loc[ind, :] = [dict_corr["Methyl"], dict_corr["RNA"], pval_mild_exp_t, pval_mild_exp_w, pval_severe_exp_t, pval_severe_exp_w, pval_mild_cpg_t, pval_mild_cpg_w, pval_severe_cpg_t, pval_severe_cpg_w, pval_acute_cpg, pval_acute_exp, pval_recover_cpg, pval_recover_exp]
        
        # make customized legend
        levels_severe = list(filter(lambda lev : SEV_CONFIG[lev] == "Severe", SEV_CONFIG.keys()))
        levels_mild = list(filter(lambda lev : SEV_CONFIG[lev] == "Mild", SEV_CONFIG.keys()))
        color_severe_arrow = set(list(map(lambda lev : PLOT_CONFIG["Longitudinal"][lev]["facecolor"], levels_severe)))
        color_mild_arrow = set(list(map(lambda lev : PLOT_CONFIG["Longitudinal"][lev]["facecolor"], levels_mild)))
        assert len(color_severe_arrow) == 1 and len(color_mild_arrow) == 1
        color_severe_arrow = list(color_severe_arrow)[0]
        color_mild_arrow = list(color_mild_arrow)[0]
        color_convalescent = set(PLOT_CONFIG["Single"].values())
        assert len(color_convalescent) == 1
        color_convalescent = list(color_convalescent)[0]
        handles = [
            mpatches.Patch(color = color_convalescent, label = "Convalescent"),
            mpatches.Patch(color = color_mild_arrow, label = "Mild\n(Confirmed)"),
            mpatches.Patch(color = color_severe_arrow, label = "Severe\n(Confirmed)")
        ]
        ax4.legend(handles = handles, loc = "center left", bbox_to_anchor = (0, 0.5))
        ax4.axis("off")
        
        list_all_sample = list(chain(*dict_sample_info["Longitudinal"])) + dict_sample_info["Single"]
        list_all_cpg = list(map(lambda sample : dict_sample_to_cpg[sample][dict_corr["Methyl"]], list_all_sample))
        list_all_exp = list(map(lambda sample : dict_sample_to_exp[sample][dict_corr["RNA"]], list_all_sample))
        cpg_minmax_scale = max(list_all_cpg) - min(list_all_cpg)
        exp_minmax_scale = max(list_all_exp) - min(list_all_exp)
        ax2.set_xlim([min(list_all_cpg) - cpg_minmax_scale*0.05, max(list_all_cpg) + cpg_minmax_scale*0.05])
        ax2.set_ylim([min(list_all_exp) - exp_minmax_scale*0.05, max(list_all_exp) + exp_minmax_scale*0.05])
        ax1.set_ylim([min(list_all_exp) - exp_minmax_scale*0.05, max(list_all_exp) + exp_minmax_scale*0.05])
        ax3.set_xlim([min(list_all_cpg) - cpg_minmax_scale*0.05, max(list_all_cpg) + cpg_minmax_scale*0.05])
        
        ax3.set_xlabel("Methylation (%)")
        ax3.set_ylabel("")
        ax1.set_xlabel("")
        ax1.set_ylabel("Expression Level (DESeq2 Normalized)")
        ax2.set_title(f"{dict_corr['Methyl']} ~ {dict_corr['RNA']}")
        ax1.yaxis.set_tick_params(labelsize = plt.rcParams["font.size"] - 1.5)
        ax2.xaxis.set_tick_params(labelsize = plt.rcParams["font.size"] - 1.5)
        ax2.yaxis.set_tick_params(labelsize = plt.rcParams["font.size"] - 1.5)
        ax3.xaxis.set_tick_params(labelsize = plt.rcParams["font.size"] - 1.5)
        gsfig.tight_layout(fig)
        # plt.show()
        plt.savefig(path_save_png_format.format(Methyl = dict_corr["Methyl"], RNA = dict_corr["RNA"]), bbox_inches = "tight")
        plt.close()
        # raise Exception
    return table_visit_compare
    
#%%
if __name__ == "__main__":
    path_methyl_rna_corr = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/Methyl_RNA_Relation/Methyl_RNA_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Mild_Sev_Visit1.LOO_common.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_6.DMP_DEG_LOO_annot.20240402.tsv"
    path_cpg = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/MethylCpGTable/Infectomics.Copy_From_HjRyu/MethylCpGTable.Control.Mild.Case.Severe.Filtered.DMP.Hyper_Hypo.Sev_vs_Mild.Visit1.LOO_common.tsv"
    path_exp = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/DEGExtract/DEGExtract.RNA_samples_with_Methyl.DEG_by_Sex.Control_M.Case_F.Cov_Age.mincount_1.20240402.tsv.normcount"
    path_meta = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Resources/MetaTable/COVID19_master_table_20231007.Methyl_Overlap.with_Severity.20240402.txt"
    
    path_save_png_format = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/Methyl_RNA_Relation/Trajactory_Plot/Trajactory_plot.Methyl_LOO_common.cis.distance_cf.1000000.Overlap_DEG_LOO.corr_0.5.log2fc_1.3.loo_6.{Methyl}.{RNA}.png"
    path_save_visit_relation_pval = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/Methyl_RNA_Relation/Trajactory_Plot/Trajactory_comparison.Methyl_LOO_common.cis.distance_cf.1000000.Overlap_DEG_LOO.corr_0.5.log2fc_1.3.loo_6.pvalue_between_visit.tsv"
    table_visit_compare = main(path_methyl_rna_corr, path_cpg, path_exp, path_meta, path_save_png_format, path_save_visit_relation_pval)
# %%
