# %%
import glob
import gzip
import os
import pickle
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
from statsmodels.stats.multitest import fdrcorrection

# %%
def filter_df_methyl_by_sample_id(df_methyl, startswith, endswith, colsample="ID"):
    df_methyl = df_methyl[np.logical_and(df_methyl[colsample].str.startswith(startswith), df_methyl[colsample].str.endswith(endswith))]

    return df_methyl

def get_dict_sample_age(df_sev_info):
    dict_sample_age = dict()
    for _, rows in df_sev_info.iterrows():
        dict_rows = dict(rows)
        sampleid = dict_rows["Sample_ID"]
        age = dict_rows["Sample_Age"]
        dict_sample_age[sampleid] = age
    
    return dict_sample_age

# %%
def main(pathsevinfo, outfilehyper, outfilehypo, list_drop_samples, outfilename):
    df_sev_info = pd.read_csv(pathsevinfo, sep="\t")
    df_methyl = df_sev_info[~df_sev_info["Sample_ID"].isin(list_drop_samples)]
    dict_sample_age = get_dict_sample_age(df_methyl)
    df_dmp_deg_hyper = pd.read_csv(outfilehyper, sep="\t")
    df_dmp_deg_hypo = pd.read_csv(outfilehypo, sep="\t")
    df_dmp_deg_all = pd.concat([df_dmp_deg_hyper, df_dmp_deg_hypo], axis=0)
    df_dmp_deg_all["Sample_Age"] = df_dmp_deg_all["ID"].apply(lambda x: dict_sample_age[x])
    list_visit = ["R1"]
    list_phase = ["Convalescent"]
    for visit, phase in zip(list_visit, list_phase):
        df_dmp_deg_all_copy = df_dmp_deg_all.copy()
        df_dmp_deg_all_filt = filter_df_methyl_by_sample_id(df_dmp_deg_all_copy, "C19", visit)
        list_array_age = df_dmp_deg_all_filt.groupby(["Marker", "GeneSymbol"])["Sample_Age"].apply(np.array).to_list()
        list_array_cpg = df_dmp_deg_all_filt.groupby(["Marker", "GeneSymbol"])["CpGbeta"].apply(np.array).to_list()
        list_markers = df_dmp_deg_all_filt.groupby(["Marker", "GeneSymbol"])["CpGbeta"].apply(np.array).index.to_list()
        # list_corr = list()
        # list_pval = list()
        # for marker, array_cpg, array_age in zip(list_markers, list_array_cpg, list_array_age):
        #     corr, pval = spearmanr(array_age, array_cpg)
        #     list_corr.append(corr)
        #     list_pval.append(pval)
        # list_rej, list_padj = fdrcorrection(list_pval, alpha=0.05)
        # print(list_rej)
        # print(list_padj)
        # print(list(zip(list_markers, list_padj)))
        fig, axes = plt.subplots(nrows=int(np.sqrt(len(list_markers)))+1, ncols=int(np.sqrt(len(list_markers)))+1, figsize=(30,30), sharex=True, sharey=True)
        for marker, array_cpg, array_age, ax in zip(list_markers, list_array_cpg, list_array_age, axes.flatten()):
            corr, pval = spearmanr(array_age, array_cpg)
            if pval < 0.05:
                ax.scatter(array_age, array_cpg, c="red")
            else:
                ax.scatter(array_age, array_cpg, c="black")
            titlename = f"{marker[1].split('_')[-1]} ({marker[0]})\ncorr:{round(corr, 2)} pval:{round(pval, 2)}"
            if pval < 0.05:
                ax.set_title(titlename, fontdict={"fontsize":16, "fontweight":"bold"})
            else:
                ax.set_title(titlename, fontdict={"fontsize":16})
        plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
        fig.supxlabel("Sample Age (Yrs)", fontsize=30)
        fig.supylabel("$\\beta$-value (%)", fontsize=30)
        fig.tight_layout()
        plt.show()
        plt.close()
# %%
if __name__ == "__main__":
    visit = "first"
    pathsevinfo = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/9_clinical/Infectomics_Severity_Information_Methyl_20240102.tsv"
    filetargetgene = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/DMPDEG/{visit}/table_deg_dmp_overlap_abslog2fc_1.3_qval_0.05_sorted.tsv"
    outfilehyper = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{visit}/hyper/samplewise_genewise_cpg_dmpdeg_overlap_20240220.tsv"
    outfilehypo = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{visit}/hypo/samplewise_genewise_cpg_dmpdeg_overlap_20240220.tsv"
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
    outfilename = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes"
    main(pathsevinfo, outfilehyper, outfilehypo, list_drop_samples, outfilename)
# %%
# Severity vs Age
# df_sev_info_filt = df_sev_info[df_sev_info["Severity"]!=0]
# sns.violinplot(data=df_sev_info_filt, y="Sample_Age", x="Severity")

# %%
