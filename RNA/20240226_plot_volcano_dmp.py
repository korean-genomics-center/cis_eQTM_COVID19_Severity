# Courtesy of Yoonsung

# %%
import math
import os

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


# %%
def get_marker_gene(df_dmp_deg):
    dict_marker_gene = dict()
    for _, rows in df_dmp_deg.iterrows():
        dict_rows = dict(rows)
        dmp_markers = dict_rows["PosName"]
        gene_symbol = dict_rows["Gene_Symbol"].split("_")[-1]
        for dmp_marker in dmp_markers:
            dict_marker_gene[dmp_marker] = gene_symbol

    return dict_marker_gene

# %%
def main(path_all_dmp, path_marker_loo, path_dmp_deg, path_volcano):
    
    plt.rcParams["font.size"] = 24

    df_dmp = pd.read_csv(path_all_dmp, sep="\t")
    df_dmp_loo = pd.read_csv(path_marker_loo, sep = '\t')
    df_dmp_deg = pd.read_csv(path_dmp_deg, sep='\t').rename(columns={"DMP_Markers": "PosName"})
    df_dmp_deg["PosName"] = df_dmp_deg["PosName"].str.split("/")
    df_dmp["PosName"] = df_dmp.apply(lambda x : f"{x['chr']}:{x['start']}", axis = 1)
    df_dmp_loo["PosName"] = df_dmp_loo.apply(lambda x : f"{x['chr']}:{x['start']}", axis = 1)

    list_posname_df_dmp_loo = df_dmp_loo["PosName"].to_list()
    list_list_posname_df_dmp_deg_per_deg = df_dmp_deg["PosName"].to_list()
    list_posname_df_dmp_deg = list()
    for list_dmp_deg_per_deg in list_list_posname_df_dmp_deg_per_deg:
        list_posname_df_dmp_deg.extend(list_dmp_deg_per_deg)

    df_dmp["MarkerType"] = "All"
    df_dmp.loc[df_dmp["PosName"].isin(list_posname_df_dmp_loo), "MarkerType"] = "Common"
    df_dmp.loc[df_dmp["PosName"].isin(list_posname_df_dmp_deg), "MarkerType"] = "DEG"

    df_dmp["minuslog10qvalue"] = df_dmp["qvalue"].apply(lambda x:-math.log10(x))

    dict_style = {
        "size" : {"All":7, "Common":14, "DEG":105},
        "color": {"All":(229/255,219/255,193/255), "Common":(254/255,207/255,141/255), "DEG":(196/255,22/255,28/255)},
        "marker":{"All":'o', "Common":'o', "DEG":'*'}
    }

    ylim  = [0.80, 1.0]
    ylim2 = [0.0, 0.80]
    ylimratio = (ylim[1]-ylim[0])/(ylim2[1]-ylim2[0]+ylim[1]-ylim[0])
    ylim2ratio = (ylim2[1]-ylim2[0])/(ylim2[1]-ylim2[0]+ylim[1]-ylim[0])
    gs = gridspec.GridSpec(2, 1, height_ratios=[ylimratio, ylim2ratio])
    fig = plt.figure(figsize=(8,10))
    ax_top = fig.add_subplot(gs[0])
    ax_bottom = fig.add_subplot(gs[1])

    sns.scatterplot(data=df_dmp[np.logical_and(df_dmp["MarkerType"] == "DEG", df_dmp["minuslog10qvalue"] > 7.5)],
                    x="meth.diff", 
                    y="minuslog10qvalue",
                    color="white", 
                    edgecolor="k", 
                    hue = "MarkerType",
                    palette = dict_style["color"],
                    size = "MarkerType",
                    sizes = dict_style["size"],
                    style = "MarkerType",
                    markers = dict_style["marker"],
                    ax=ax_top,
                    zorder=100)
    sns.scatterplot(data=df_dmp[np.logical_and(df_dmp["MarkerType"]  == "Common", df_dmp["minuslog10qvalue"] > 7.5)], 
                    x="meth.diff", 
                    y="minuslog10qvalue",
                    color="white",
                    edgecolor="k", 
                    hue = "MarkerType",
                    palette = dict_style["color"],
                    size = "MarkerType",
                    sizes = dict_style["size"],
                    style = "MarkerType",
                    markers = dict_style["marker"],
                    ax=ax_top,
                    zorder=10)
    sns.scatterplot(data=df_dmp[np.logical_and(df_dmp["MarkerType"]  == "DEG", df_dmp["minuslog10qvalue"] < 7.5)],
                    x="meth.diff", 
                    y="minuslog10qvalue",
                    color="white", 
                    edgecolor="k", 
                    hue = "MarkerType",
                    palette = dict_style["color"],
                    size = "MarkerType",
                    sizes = dict_style["size"],
                    style = "MarkerType",
                    markers = dict_style["marker"],
                    ax=ax_bottom,
                    zorder=100)
    sns.scatterplot(data=df_dmp[np.logical_and(df_dmp["MarkerType"]  == "Common", df_dmp["minuslog10qvalue"] < 7.5)], 
                    x="meth.diff", 
                    y="minuslog10qvalue",
                    color="white",
                    edgecolor="k", 
                    hue = "MarkerType",
                    palette = dict_style["color"],
                    size = "MarkerType",
                    sizes = dict_style["size"],
                    style = "MarkerType",
                    markers = dict_style["marker"],
                    ax=ax_bottom,
                    zorder=10)
    sns.scatterplot(data=df_dmp[np.logical_and(df_dmp["MarkerType"]  == "All", df_dmp["minuslog10qvalue"] < 7.5)], 
                    x="meth.diff", 
                    y="minuslog10qvalue",
                    color="white",
                    edgecolor="k", 
                    hue = "MarkerType",
                    palette = dict_style["color"],
                    size = "MarkerType",
                    sizes = dict_style["size"],
                    style = "MarkerType",
                    markers = dict_style["marker"],
                    ax=ax_bottom,
                    zorder=-999)
    
    ax_top.set_xticks([])
    ax_top.set_ylabel("")
    ax_bottom.set_ylabel("")
    ax_top.set_xlabel("")
    ax_bottom.set_xlabel("")
    ax_top.set_xlim(-45, 45)
    ax_bottom.set_xlim(-45,45)
    ax_top.set_ylim(7.5, 60)
    ax_bottom.set_ylim(-0.5,7.5)
    ax_top.get_legend().remove()

    dict_gene_loc = {
        "xy" : {
            "ANK1" : (-10, 0),
            "HPD" : (-7, 0),
            "SMIM24" : (-10, -3),
            "GALNT14" : (-10, 3),
            "RAB13" : (10, -3),
            "SGIP1": (10, -3),
            "SEC14L5": (10, -3),
            "S100P": (-40,-3),
            "SMIM24": (-45,-1),
            "TESC": (-35, -3),
            "PTPRF": (-40, 3),
            "RCVRN": (-40, 3),
            "LCN2": (-30, -1),
            "FOXO3": (40, -5),
            "CHIT1": (-35, 5),
            "HDGF": (35, 3),
            "SRXN1": (10, -7),
        },
        "ha" : {
            "ANK1" : "right",
            "HPD" : "right",
            "SMIM24" : "right",
            "GALNT14" : "right",
            "S100P": "left",
            "SMIM24": "left",
            "TESC" : "left",
            "PTPRF" : "left",
            "RCVRN" : "left",
            "LCN2": "left",
            "FOXO3": "right",
            "CHIT1": "left",
            "HDGF": "right",
            "SRNX1": "right"
        }
    }

    dict_marker_gene = get_marker_gene(df_dmp_deg)
    df_dmp_selected_only_deg = df_dmp[df_dmp["MarkerType"] == "DEG"]
    df_dmp_selected_only_deg["GeneSymbol"] = df_dmp_selected_only_deg["PosName"].apply(lambda x: dict_marker_gene[x] if str(x) in dict_marker_gene.keys() else x)
    for i, rows in df_dmp_selected_only_deg.iterrows():
        dict_rows = dict(rows)
        gene = str(dict_rows["GeneSymbol"])
        x = float(dict_rows["meth.diff"])
        y = float(dict_rows["minuslog10qvalue"])
        ax_top.annotate(gene, xy=(x,y), xycoords='data', xytext=dict_gene_loc["xy"].get(gene, (10,0)), textcoords='offset points', color='k', fontsize=plt.rcParams["font.size"]-15, weight="bold", zorder=200, ha = dict_gene_loc["ha"].get(gene, "left"))
        ax_bottom.annotate(gene, xy=(x,y), xycoords='data', xytext=dict_gene_loc["xy"].get(gene, (10,0)), textcoords='offset points', color='k', fontsize=plt.rcParams["font.size"]-15, weight="bold", zorder=200, ha = dict_gene_loc["ha"].get(gene, "left"))

    ax_top.axvline(x=-10, linestyle="dashed", color="grey", zorder=0)
    ax_top.axvline(x=10, linestyle="dashed", color="grey", zorder=0)
    ax_bottom.axvline(x=-10, linestyle="dashed", color="grey", zorder=0)
    ax_bottom.axvline(x=10, linestyle="dashed", color="grey", zorder=0)
    ax_bottom.axhline(y=-math.log10(0.05), linestyle="dashed", color="grey", zorder=0)
    fig.supxlabel("Mean Methylation Difference (%)", fontsize=plt.rcParams["font.size"]-3, x=0.55, y=0.05)
    fig.supylabel("-log10(FDR)", fontsize=plt.rcParams["font.size"]-3, x=0.1, y=0.5)
    legend = plt.legend(loc = "lower right", title = "Marker Type", fontsize=plt.rcParams["font.size"]-9)
    legend.get_title().set_fontsize(fontsize=plt.rcParams["font.size"]-7)
    legend.get_title().set_fontweight(weight="bold")
    legend.get_texts()[0].set_text("LOOC DMP & DEG")
    legend.get_texts()[1].set_text("LOOC DMP")
    legend.get_texts()[2].set_text("NS DMP")
    fig.tight_layout()
    plt.savefig(path_volcano, dpi=600, bbox_inches="tight")
    plt.show()
    plt.close()

# %%
if __name__ == "__main__":
    visit = "first"
    path_all_dmp = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_allSamples/DMPExtract/Methylation_DMP_Extract.Control_Mild.Case_Severe.tsv"
    path_marker_loo = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231211/severe_mild_firstVisit/methyl_hyperhypo_severe_mild_commonMarkers_LOO.tsv"
    path_dmp_deg = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/DMPDEG/first/table_deg_dmp_overlap_abslog2fc_1.3_qval_0.05_sorted.tsv"
    path_volcano = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/InfectomicsPaper1/volcanoplot_methyl_dmpdeg_overlap.png"
    main(path_all_dmp, path_marker_loo, path_dmp_deg, path_volcano)
# %%
