#%%
import pandas as pd
import numpy as np
from glob import glob
import seaborn as sns
from matplotlib import pyplot as plt


CONFIG_DMP = {
    "featname" : ["chr", "start", "end"],
    "fdr" : "qvalue",
    "magnitude" : "meth.diff",
    "fdr_cf" : 0.05,
    "magnitude_cf" : 10
}

CONFIG_DEG = {
    "featname" : ["ID"],
    "fdr" : "padj",
    "magnitude" : "log2FoldChange",
    "fdr_cf" : 0.05,
    "magnitude_cf" : 1.3
}

N_LOO = 8

def main(path_dmp_loo_glob, path_deg_loo_glob, path_dmp_correlation, path_save_median_values):
    dict_dmp_loo_stat = read_loo_stat(path_dmp_loo_glob, CONFIG_DMP)
    dict_deg_loo_stat = read_loo_stat(path_deg_loo_glob, CONFIG_DEG)
    
    table_corr = pd.read_csv(path_dmp_correlation, sep = '\t')
    table_corr_loo_stat_added = add_median_loo_statistic_to_correlation_table(table_corr, dict_dmp_loo_stat, "Methyl", "dmp_median_fdr", "dmp_median_methdiff")
    table_corr_loo_stat_added = add_median_loo_statistic_to_correlation_table(table_corr_loo_stat_added, dict_deg_loo_stat, "RNA", "deg_median_fdr", "deg_median_log2FC")
    
    table_corr_loo_stat_added.to_csv(path_save_median_values, sep = '\t', index = False)
    visualize_dmp_deg_correlation(table_corr_loo_stat_added)

def read_loo_stat(path_loo_glob, config):
    list_path_loo = glob(path_loo_glob)
    dict_loo_stat = dict()
    for path_file in list_path_loo:
        table = pd.read_csv(path_file, sep = '\t')
        table["featname"] = table[config["featname"]].apply(lambda row : '_'.join(list(map(str, row.to_list()))), axis = 1)
        for _, row in table.iterrows():
            featname = row["featname"]
            fdr = row[config["fdr"]]
            mag = row[config["magnitude"]]
            if dict_loo_stat.get(featname) == None:
                dict_loo_stat[featname] = {
                    "fdr" : list(),
                    "magnitude" : list()
                }
            dict_loo_stat[featname]["fdr"].append(fdr)
            dict_loo_stat[featname]["magnitude"].append(mag)
    return dict_loo_stat

def add_median_loo_statistic_to_correlation_table(table, dict_stat, col_name, col_add_fdr, col_add_magnitude, check_n_overlap = 8):
    for name in table[col_name].unique():
        assert len(dict_stat[name]["fdr"]) >= check_n_overlap, f"{name} : {len(dict_stat[name]['fdr'])}"
        median_fdr = np.median(dict_stat[name]["fdr"])
        median_mag = np.median(dict_stat[name]["magnitude"])
        table.loc[table[col_name] == name, col_add_fdr] = median_fdr
        table.loc[table[col_name] == name, col_add_magnitude] = median_mag
    return table
        
def visualize_dmp_deg_correlation(table):
    sns.scatterplot(data = table, x = "dmp_median_methdiff", y = "deg_median_log2FC", hue = "corr_rho", palette = "coolwarm", vmin = -0.7, vmax = 0.7, hue_norm=(-0.7, 0.7))
    plt.xlabel("Median Methylation Difference (%)")
    plt.ylabel("Median Log2 Fold Change")
    plt.legend(title = "Spearman Rho")
    plt.axhline(1.3, linestyle = "--", linewidth = "2", c = 'moccasin')
    plt.axhline(-1.3, linestyle = "--", linewidth = "2", c = 'moccasin')
    plt.axvline(10, linestyle = "--", linewidth = "2", c = 'paleturquoise')
    plt.axvline(-10, linestyle = "--", linewidth = "2", c = 'paleturquoise')
    plt.show()

#%%
if __name__ == "__main__":
    path_dmp_loo_glob = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231211/severe_mild_firstVisit/DMPExtract/Methylation_DMP_Extract.Control_Mild.Case_Severe_WO_*.filtered.fdr_05.methdiff_abs10.tsv"
    path_deg_loo_glob = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/5_deg/Visit1_Severe__Visit1_Mild_removed_loo_*_20240327.tsv"
    path_dmp_correlation = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/Methyl_RNA_Relation/Methyl_RNA_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Mild_Sev_Visit1.LOO_common.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_6.20240326.tsv"
    
    path_save_median_values = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/Methyl_RNA_Relation/Methyl_RNA_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Mild_Sev_Visit1.LOO_common.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_6.DMP_DEG_LOO_annot.20240402.tsv"
    
    main(path_dmp_loo_glob, path_deg_loo_glob, path_dmp_correlation, path_save_median_values)
# %%
