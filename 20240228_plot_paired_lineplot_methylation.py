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
from scipy.stats import ttest_rel
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

def make_dataframe_paired_CpG(dict_phase_marker_methyl, list_sample, marker, columns=["Sample_ID", "Phase", "CpG"]):
    ind = 0
    df_paired_marker = pd.DataFrame(columns=columns)
    for sample in list_sample:
        sample_acute = sample
        acute_val = dict_phase_marker_methyl["Acute"][marker].get(sample_acute)
        sample_recover = sample.replace("-V1", "-V4")
        rec_val = dict_phase_marker_methyl["Recovery"][marker].get(sample_recover)
        if rec_val == None:
            sample_recover = sample.replace("-V1", "-V3")
            rec_val = dict_phase_marker_methyl["Recovery"][marker].get(sample_recover)
        if acute_val == None or rec_val == None:
            continue
    
        df_paired_marker.loc[ind, :] = [sample, 0.25, acute_val]
        ind += 1
        df_paired_marker.loc[ind, :] = [sample, 0.75, rec_val]
        ind += 1

    return df_paired_marker


def draw_paired_lineplot(df_paired_marker, stat_test, color, ax, draw_ytick=True):
    df_paired_marker[["Phase", "CpG"]] = df_paired_marker[["Phase", "CpG"]].astype(float)
    sns.scatterplot(data = df_paired_marker, x = "Phase", y = "CpG", edgecolor='k', color=color, ax = ax, legend = False, zorder = 999)
    sns.lineplot(data = df_paired_marker, x = "Phase", y = "CpG", linestyle="dotted", linewidth=1.5, palette=["grey"]*len(df_paired_marker["Sample_ID"].unique()), hue="Sample_ID", ax = ax, legend = False)
    sns.regplot(data = df_paired_marker, x = "Phase", y = "CpG", color="k", scatter=False, ci=False, scatter_kws={"zorder":100, "alpha":1}, ax = ax)

    pairs = [0.25, 0.75]
    pvalue = round(float(stat_test.loc["P-value"]), 3)
    star = "ns"
    if 0.01 <= pvalue < 0.05:
        star = "*"
    elif 0.005 <= pvalue < 0.01:
        star = "**"
    elif 0.001 <= pvalue < 0.005:
        star = "***"
    elif pvalue < 0.001:
        star = "****"
    x1, x2 = pairs 
    y, h, color = 100+5, 2, 'k'

    ax.plot([x1, x2], [y, y], lw=1.5, c=color)
    ax.text(x=(x1+x2)*.5, y=y+h, s=star, fontdict={'fontsize':10}, ha='center', va='bottom', color=color)    
    
    ax.set_xlabel("")
    ax.set_ylabel("Proportion of methylated CpGs \n ($\\beta$-value)")
    ax.set_xticks([0.25, 0.75], ["Acute", "Recovery"])
    ax.set_xlim(0, 1)
    ax.set_ylim(-10, 110)
    yticks = list(range(0, 110, 20))
    ax.set_yticks(yticks)
    ax.grid(axis="y")
    ax.spines[['top', 'right']].set_visible(False)
    if not draw_ytick:
        ax.spines['left'].set_visible(False)
        # ax.set_yticks(yticks, [""]*len(yticks))
        ax.tick_params(left=False, labelleft=False)
        ax.set_ylabel("")

def meandiff(a, b):
    deltamean = np.mean(b) - np.mean(a)

    return deltamean

def conduct_paired_t_test(df_paired_CpG):
    array_cpg_first = list(df_paired_CpG.groupby("Phase")["CpG"].apply(np.array))[0]
    array_cpg_last = list(df_paired_CpG.groupby("Phase")["CpG"].apply(np.array))[1]
    deltamean = meandiff(array_cpg_first, array_cpg_last)
    tstat, pval = ttest_rel(array_cpg_first, array_cpg_last)
    dict_stat = {"Delta": deltamean, "Statistics": tstat, "P-value": pval}
    result_test = pd.DataFrame.from_dict(dict_stat, orient="index")

    return result_test


# %%
def main(pathsevinfo, filehyper, filehypo, list_drop_samples, outdir):
    df_sev_info = pd.read_csv(pathsevinfo, sep="\t")
    df_methyl = df_sev_info[~df_sev_info["Sample_ID"].isin(list_drop_samples)]
    dict_sample_age = get_dict_sample_age(df_methyl)
    df_dmp_deg_hyper = pd.read_csv(filehyper, sep="\t")
    df_dmp_deg_hypo = pd.read_csv(filehypo, sep="\t")
    df_dmp_deg_all = pd.concat([df_dmp_deg_hyper, df_dmp_deg_hypo], axis=0)
    df_dmp_deg_all_filtered = df_dmp_deg_all[~df_dmp_deg_all["ID"].isin(list_drop_samples)]
    df_dmp_deg_all["Sample_Age"] = df_dmp_deg_all_filtered["ID"].apply(lambda x: dict_sample_age[x])
    list_visit = ["V1", ("V3", "V4")]
    list_phase = ["Acute", "Recovery"]
    
    dict_phase_marker_methyl = dict()
    for visit, phase in zip(list_visit, list_phase):
        df_dmp_deg_all_copy = df_dmp_deg_all.copy()
        df_dmp_deg_all_filt = filter_df_methyl_by_sample_id(df_dmp_deg_all_copy, "C19", visit)
        list_samples = df_dmp_deg_all_filt.groupby(["Marker", "GeneSymbol"])["ID"].apply(np.array).to_list()
        list_array_cpg = df_dmp_deg_all_filt.groupby(["Marker", "GeneSymbol"])["CpGbeta"].apply(np.array).to_list()
        list_markers = df_dmp_deg_all_filt.groupby(["Marker", "GeneSymbol"])["CpGbeta"].apply(np.array).index.to_list()
        dict_marker_methyl = dict()
        for sample, marker, array_cpg in zip(list_samples, list_markers, list_array_cpg):
            dict_marker_methyl[marker] = dict(zip(sample, array_cpg))
            
        dict_phase_marker_methyl[phase] = dict_marker_methyl
    
    list_marker = list(dict_phase_marker_methyl["Acute"].keys())
    list_sample = list(list(dict_phase_marker_methyl["Acute"].values())[0].keys())
    dict_sample_severity = dict(zip(df_sev_info["Sample_ID"], df_sev_info["Severity_group"]))

    for marker in list_marker:
        list_sample_mild = list(filter(lambda x: dict_sample_severity[x] == "Mild", list_sample))
        list_sample_severe = list(filter(lambda x: dict_sample_severity[x] == "Severe", list_sample))

        df_paired_CpG_mild = make_dataframe_paired_CpG(dict_phase_marker_methyl, list_sample_mild, marker)
        df_paired_CpG_severe = make_dataframe_paired_CpG(dict_phase_marker_methyl, list_sample_severe, marker)
        
        result_test_mild = conduct_paired_t_test(df_paired_CpG_mild)
        result_test_severe = conduct_paired_t_test(df_paired_CpG_severe)

        plt.rcParams["font.size"] = 14
        fig, ax = plt.subplots(1, 2, figsize = (5, 4), layout='constrained')
        draw_paired_lineplot(df_paired_CpG_mild, result_test_mild, color="forestgreen", ax=ax[0])
        draw_paired_lineplot(df_paired_CpG_severe, result_test_severe, color="firebrick", ax=ax[1], draw_ytick=False)
        fig.suptitle(marker[0] + " " + f"({'_'.join(marker[-1].split('_')[1:])})")
        fig.supxlabel("Infection Phase", y=-0.05)
        plt.subplots_adjust(wspace=0)

        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor='forestgreen',label='Mild'),
                            Patch(facecolor='firebrick',label='Severe')]
        legend = plt.legend(handles=legend_elements, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=plt.rcParams["font.size"]-5)
        legend.set_title("Severity", prop={"weight":"bold", "size":plt.rcParams["font.size"]-3})  

        outpirpairedplot = os.path.join(outdir, "pairedplot")
        os.makedirs(outpirpairedplot, exist_ok=True)
        outfilename = os.path.join(outpirpairedplot, f"pairedplot_{marker[0] + '_' + '_'.join(marker[-1].split('_')[1:])}.png")
        plt.savefig(outfilename, dpi=600, bbox_inches="tight")
        plt.show()
        plt.close()

# %%
if __name__ == "__main__":
    visit = "first"
    pathsevinfo = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/9_clinical/Infectomics_Severity_Information_Methyl_20240102.tsv"
    filehyper = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{visit}/hyper/samplewise_genewise_cpg_dmpdeg_overlap_20240220.tsv"
    filehypo = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{visit}/hypo/samplewise_genewise_cpg_dmpdeg_overlap_20240220.tsv"
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
    outdir = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/InfectomicsPaper1"
    main(pathsevinfo, filehyper, filehypo, list_drop_samples, outdir)
 # %%
