# %%
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import levene, mannwhitneyu, ranksums, ttest_ind
from statannot import add_stat_annotation
from statsmodels.stats.multitest import fdrcorrection


# %%
def main(path_mean_beta,path_test_result, path_boxplot):

    df_mean_beta = pd.read_csv(path_mean_beta, sep="\t")
    df_mean_beta_filt = df_mean_beta[df_mean_beta["Visit"] != "LongCOVID"]
    df_mean_beta_filt = df_mean_beta_filt[df_mean_beta_filt["Visit"] != "Convalescent"]
    
    # dict_order_visit = {"Control": 0, "First": 1, "Last": 2, "Convalescent": 3, "LongCOVID": 4}
    dict_order_visit = {"Control": 0, "First": 1, "Last": 2}
    dict_order_sev = {"Healthy": 0, "Mild": 0.1, "Severe": 0.2}

    dict_cnt = dict()
    dict_grpby_sev_visit = df_mean_beta_filt.groupby(["Severity", "Visit"]).apply(len).to_dict()
    for key, cnt in dict_grpby_sev_visit.items():
        sev_visit = "_".join(key)
        dict_cnt[sev_visit] = cnt

    list_severity = sorted(list(df_mean_beta_filt["Severity_Visit"].unique()), key=lambda x:dict_order_visit[x.split("_")[-1]]+dict_order_sev[x.split("_")[0]])
    
    list_severity_cnt = list()
    for sev in list_severity:
        if sev in dict_cnt.keys():
            cnt = dict_cnt[sev]
            if sev.split("_")[0] == "Healthy":
                sev = sev.split("_")[0] + "\n" + "Controls"
            else:
                sev = sev.split("_")[0] + "\n" + "Group"
            sev_nline = '\n'.join(sev.split('_'))
            list_severity_cnt.append(sev_nline)

    list_comb_sev_visit = list(combinations(list_severity, 2))
    # list_comb_sev_visit_rmv_diff_visit = list(filter(lambda x: (x[0].split("_")[1] == x[1].split("_")[1]) or (x[0] == "Healthy_Control" and "First" in x[1]), list_comb_sev_visit))
    list_comb_sev_visit_rmv_diff_visit = list(filter(lambda x: (x[0].split("_")[1] == x[1].split("_")[1] or (x[0] == "Healthy_Control")), list_comb_sev_visit))
    dict_sev_visit = df_mean_beta_filt.groupby("Severity_Visit")["MeanBeta"].apply(np.array).to_dict()

    from collections import defaultdict
    dict_meanbeta_comb_visits = defaultdict(list)
    for comb_visits in list_comb_sev_visit_rmv_diff_visit:
        for visit in comb_visits:
            meanbeta = dict_sev_visit[visit]
            dict_meanbeta_comb_visits[comb_visits].append(meanbeta)

    def meandiff(a, b):
        deltamean = np.mean(b) - np.mean(a)

        return deltamean

    dict_visit_stats = dict()
    for visits, meanbetas in dict_meanbeta_comb_visits.items():
        deltamean = meandiff(*meanbetas)
        ### ranksum test
        stat, pval = ranksums(*meanbetas)
        ## t-test
        # lstat, lpval = levene(*meanbetas, center="mean")
        # if lpval < 0.05:
        #     stat, pval = ttest_ind(*meanbetas, equal_var=False)
        # else:
        #     stat, pval = ttest_ind(*meanbetas, equal_var=True)
        dict_visit_stats[visits] = {"Delta": deltamean, "Statistics": stat, "P-value": pval}
    result_test = pd.DataFrame.from_dict(dict_visit_stats, orient="index").reset_index(drop=False)
    result_test = result_test.rename(columns={"level_0": "Group1", "level_1": "Group2"})
    list_sigtrue = fdrcorrection(result_test["P-value"])[0]
    list_padj = fdrcorrection(result_test["P-value"])[1]
    result_test["Padj"] = list_padj
    result_test["TestSignificant"] = list_sigtrue
    print(result_test)
    result_test.to_csv(path_test_result, sep="\t", index=False)

    palette = dict()
    for sev in list_severity:
        if "Mild" in sev:
            palette[sev] = "forestgreen"
        elif "Severe" in sev:
            palette[sev] = "firebrick"
        else:
            palette[sev] = "royalblue"

    plt.figure(figsize=(5,5))
    flierprops = dict(marker='o', markerfacecolor='None', markersize=5, markeredgecolor='black')
    ax = sns.boxplot(data=df_mean_beta, 
                     x="Severity_Visit", 
                     y="MeanBeta", 
                     order=list_severity, 
                     flierprops=flierprops, 
                     palette=palette)
    pairs = [(i[1]["Group1"], i[1]["Group2"]) for i in result_test.iterrows()]   
    pvalues = [i[1]["P-value"] for i in result_test.iterrows()]
    add_stat_annotation(ax, data=df_mean_beta, x="Severity_Visit", y="MeanBeta", order=list_severity,
                    box_pairs=pairs,
                    perform_stat_test=False, pvalues=pvalues,
                    test=None, text_format='star', loc='inside', verbose=2)

    if direction == "hyper":
        plt.ylim(0, 100)
    else:
        plt.ylim(0, 100)
    
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_xlabel("Severity", fontsize=0)
    ax.set_ylabel("Mean ratio of methylated CpGs ($\\beta$-value)", fontsize=12)
    ax.set_xticks(ticks=list(range(len(list_severity_cnt))), labels=list_severity_cnt)
    ax.plot([0.7,2.3],[-.15,-.15], color="k", transform=ax.get_xaxis_transform(), linewidth=2.5, clip_on=False)
    ax.plot([2.6,4.2],[-.15,-.15], color="k", transform=ax.get_xaxis_transform(), linewidth=2.5, clip_on=False)
    plt.text(x=0.34, y=-0.05, s="Acute Phase", weight="bold", transform=plt.gcf().transFigure)
    plt.text(x=0.61, y=-0.05, s="Recovery Phase", weight="bold", transform=plt.gcf().transFigure)
    plt.savefig(path_boxplot, dpi=600, bbox_inches="tight")
    plt.show()
    plt.close()

# %%
if __name__ == "__main__":
    list_direction = ["hyper", "hypo"]
    list_visit = ["first"]
    # list_visit = ["first", "last"]
    for direction in list_direction:
        for visit in list_visit:
            # path_mean_beta = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{visit}/{direction}/mean_beta_by_severity_scale_in_detail_table_20231220.tsv"
            path_mean_beta = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{visit}/{direction}/mean_beta_by_severity_table_dmpdeg_20240229.tsv"
            path_test_result = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{visit}/{direction}/stat_test_meanbeta_severity_table_dmpdeg_20240229.tsv"
            # path_boxplot = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{visit}/{direction}/boxplot_meanbeta_severity_20240229.png"
            path_boxplot = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/InfectomicsPaper1/boxplot_meanbeta_severity.png"
            main(path_mean_beta, path_test_result, path_boxplot)
    
# %%
