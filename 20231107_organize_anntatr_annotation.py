#courtsey of yoonsung

#%%
import pandas as pd
import os
visit = "first"
direction = "hyper"
path_annotatr_result = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231211/severe_mild_{visit}Visit/annotation/methyl_{direction}_severe_mild_{visit}Visit_annotatr.tsv"
table_annot = pd.read_csv(path_annotatr_result, sep = '\t')
table_annot["position"] = table_annot.apply(lambda x : f"{x['seqnames']}_{x['start']}_{x['end']}", axis = 1)

# %%
def get_cpg_info(table):
    table_cpg = table[list(map(lambda x : "cpg" in x, table["annot.type"]))]
    if table_cpg.shape[0] == 0:
        cpg_info = "unrecog"
    else:
        cpg_info = table_cpg["annot.type"].to_list()[0].split('_')[-1]
    if cpg_info == "inter":
        cpg_info = "OpenSea"
    elif cpg_info == "shores":
        cpg_info = "CpG Shore"
    elif cpg_info == "shelves":
        cpg_info = "CpG Shelf"
    elif cpg_info == "islands":
        cpg_info = "CpG Island"
    elif cpg_info == "unrecog":
        cpg_info = "Unknown"
    else:
        raise Exception(f"Undefined CpG type : {cpg_info}")
    return cpg_info

def is_enhancer(table):
    table_enhanc = table[table["annot.type"] == "hg38_enhancers_fantom"]
    is_enhanc = table_enhanc.shape[0] > 0
    return is_enhanc

#%%

def organize_table_by_chr(table_chr):
    table_organized = pd.DataFrame(columns = ["chr", "start", "end", "CpGs", "Gene_Symbol", "Gene_Attrib", "Is_Enhancer", "PosName"])

    name_conv = {
        "1to5kb" : "Far_Promoter",
        "promoter" : "Promoter",
        "3UTR" : "3'UTR",
        "5UTR" : "5'UTR",
        "exon" : "Exon",
        "intron" : "Intron",
        "intronexonboundary" : "IntronExonBoundary"
    }
    list_chr = list()
    list_start = list()
    list_end = list()
    list_cpg = list()
    list_symb = list()
    list_attr = list()
    list_enh = list()
    list_pos = list()

    ind_table_org = 0
    for pos in table_chr["position"].unique():
        chrname, start, end = pos.split('_')
        table_part = table_chr.loc[pos, :].copy()
        if isinstance(table_part, pd.DataFrame):
            table_part = table_part.reset_index(drop = True)
        else:
            table_part = pd.DataFrame(table_part).T.reset_index(drop = True)
        cpg_type = get_cpg_info(table_part)
        is_enhanc = is_enhancer(table_part)
        table_gene = table_part.dropna(subset = ["annot.gene_id", "annot.tx_id"], how = "all")
        list_annot_id = list(map(lambda x : x.split(':')[0], table_gene["annot.id"]))
        list_gene_symbol = table_gene["annot.symbol"].to_list()
        list_gene_tx = table_gene["annot.tx_id"].to_list()
        list_gene_name = list(map(lambda x, y: x if not pd.isna(x) else y, list_gene_symbol, list_gene_tx))
        list_gene_annot = list(map(lambda x,y : f"{x}::{name_conv[y]}", list_gene_name, list_annot_id))
        if len(list_gene_annot) == 0:
            gene_symbol = "None"
            annot_id = "Intergenic"
            list_chr.append(chrname)
            list_start.append(start)
            list_end.append(end)
            list_cpg.append(cpg_type)
            list_symb.append(gene_symbol)
            list_attr.append(annot_id)
            list_enh.append(is_enhanc)
            list_pos.append(pos)
            # list_row = [chrname, start, end, cpg_type, gene_symbol, annot_id, is_enhanc, pos]
            # table_organized.loc[ind_table_org, :] = list_row
            # ind_table_org += 1
        else:
            for gene_annot in set(list_gene_annot):
                gene_symbol, annot_id = gene_annot.split("::")
                
                list_chr.append(chrname)
                list_start.append(start)
                list_end.append(end)
                list_cpg.append(cpg_type)
                list_symb.append(gene_symbol)
                list_attr.append(annot_id)
                list_enh.append(is_enhanc)
                list_pos.append(pos)
                # list_row = [chrname, start, end, cpg_type, gene_symbol, annot_id, is_enhanc, pos]
                # table_organized.loc[ind_table_org, :] = list_row
                # ind_table_org += 1
    table_organized["chr"] = list_chr
    table_organized["start"] = list_start
    table_organized["end"] = list_end
    table_organized["CpGs"] = list_cpg
    table_organized["Gene_Symbol"] = list_symb
    table_organized["Gene_Attrib"] = list_attr
    table_organized["Is_Enhancer"] = list_enh
    table_organized["PosName"] = list_pos

    # table_organized.to_csv(f"/BiO/Research/Project2/CardiomicsMethylome/Cardiomics_10_fold_cross_validation_DMP.20230302/Results/MethylCpGTable/Cardiomics.Raw/Annotation/Cardiomics.MinPerGroup.NULL.CovCutoff.10.control_case.pca_filtered.case_filtered.extra_control_filtered.20230227.annotation.Annoatr.organized.{chrname}.tsv", sep = '\t', index = False)
    return table_organized

#%%
from multiprocessing import Pool

table_annot = table_annot.set_index("position", drop = False)

list_chr = table_annot["seqnames"].unique()
list_table_chr = list(map(lambda x : table_annot[table_annot["seqnames"] == x], list_chr))
#%%
with Pool(processes=len(list_chr)) as p:
    list_table_organized = p.map(organize_table_by_chr, list_table_chr)
#%%
df_tot = pd.concat(list_table_organized)
df_tot = df_tot.set_index("PosName", drop = False)
# list_chr = list(range(1,23)) + ['X']
# path_format = "/BiO/Research/Project2/CardiomicsMethylome/Cardiomics_10_fold_cross_validation_DMP.20230302/Results/MethylCpGTable/Cardiomics.Raw/Annotation/Cardiomics.MinPerGroup.NULL.CovCutoff.10.control_case.pca_filtered.case_filtered.extra_control_filtered.20230227.annotation.Annoatr.organized.chr{chrname}.tsv"
# df_tot = pd.DataFrame()
# for chrname in list_chr:
#     path_table = path_format.format(chrname = chrname)
#     table = pd.read_csv(path_table, sep = '\t')
#     df_tot = pd.concat([df_tot, table])

# %%
df_tot_uniq = df_tot.drop_duplicates(subset = ["chr", "start"])
df_tot_uniq["start"] = df_tot_uniq["start"].astype(int)
df_tot_uniq = df_tot_uniq.sort_values(by = ["chr", "start"])

def get_neighbor_methyl(table_part):
    list_neighbor = [[table_part.iloc[0, :]["PosName"]]]
    for _, row in table_part.iloc[1:, :].iterrows():
        bef_pos = int(list_neighbor[-1][-1].split('_')[-1])
        row_pos = int(row["start"])
        if abs(row_pos - bef_pos) <= 100:
            list_neighbor[-1].append(row["PosName"])
        else:
            list_neighbor.append([row["PosName"]])
    return list_neighbor

from joblib import Parallel, delayed

with Parallel(n_jobs = len(list_chr)) as parallel:
    list_results_per_chr = parallel(delayed(get_neighbor_methyl)(df_tot_uniq[df_tot_uniq["chr"] == chrname]) for chrname in list_chr)

from itertools import chain
list_neighbors = list(chain(*list_results_per_chr))

from collections import Counter
neighbor_size = list(map(lambda x : len(x), list_neighbors))
count_neighbor_size = Counter(neighbor_size)

# %%
import matplotlib.pyplot as plt
# plt.rcParams["font.size"] = 16
# plt.rcParams["font.family"] = "serif"
plt.bar(range(1, max(count_neighbor_size.keys())+1), list(map(lambda x : count_neighbor_size[x], range(1, max(count_neighbor_size.keys())+1))), color = "gray", zorder = 3)
plt.xticks(range(1, max(count_neighbor_size.keys())+1))
plt.yscale("log")
plt.ylim(0, 1000)
plt.grid(axis = "y", linestyle = "--", zorder = 1)
plt.ylabel("Number of Biomarkers")
plt.xlabel("Size of Neighbor CpG Biomarkers")
print('\n'.join(df_tot_uniq.loc[list(chain(*list(filter(lambda x : len(x) > 2, list_neighbors))))]["Gene_Symbol"].unique()))

# %%
df_tot_target = df_tot_uniq.loc[list(chain(*list(filter(lambda x : len(x) > 2, list_neighbors))))]

# %%
import numpy as np
def get_gene_symbol_from_enst_id(df_enst):
    dict_enst = dict()
    for transcript_id, gene_name in zip(df_enst["ensembl_transcript_id_version"], df_enst["external_gene_name"]):
            dict_enst[transcript_id] = gene_name

    return dict_enst

# %%
print('\n'.join(list(filter(lambda x : x.startswith("ENST"), df_tot_target["Gene_Symbol"].unique()))))
# %%
enst_hyper = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231211/severe_mild_firstVisit/annotation/methyl_hyper_severe_mild_firstVisit_annotatr_neighbor_methylation_enst_conv.tsv"
df_enst_hyper = pd.read_csv(enst_hyper, sep='\t')
dict_enst_hyper = get_gene_symbol_from_enst_id(df_enst_hyper)
df_tot_target["Gene_Symbol"] = df_tot_target["Gene_Symbol"].apply(lambda x: dict_enst_hyper[x] if str(x).startswith("ENST") else x)

# %%
enst_hypo = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231211/severe_mild_firstVisit/annotation/methyl_hypo_severe_mild_firstVisit_annotatr_neighbor_methylation_enst_conv.tsv"
df_enst_hypo = pd.read_csv(enst_hypo, sep='\t')
dict_enst_hypo = get_gene_symbol_from_enst_id(df_enst_hypo)
df_tot_target["Gene_Symbol"] = df_tot_target["Gene_Symbol"].apply(lambda x: dict_enst_hypo[x] if str(x).startswith("ENST") else x)

#%%
def get_cpg_attrib_heatmap(table_org):
    table_heatmap = pd.DataFrame(columns = ["OpenSea", "CpG Shelf", "CpG Shore", "CpG Island", "Unknown"], index = ["Enhancer", "Promoter", "Far_Promoter", "5'UTR", "Intron", "Exon", "IntronExonBoundary", "3'UTR", "Intergenic", "Total"])
    table_heatmap = table_heatmap.fillna(0)
    for pos in set(table_org.index):
        table_part = table_org.loc[pos, :].copy()
        if not isinstance(table_part, pd.DataFrame):
            table_part = pd.DataFrame(table_part).T
        table_part = table_part.reset_index(drop = True)
        row = table_part.iloc[0, :]
        is_enhanc = row["Is_Enhancer"]
        cpg_type = row["CpGs"]
        if is_enhanc:
            table_heatmap.loc["Enhancer", cpg_type] += 1
        for attr in table_part["Gene_Attrib"].unique():
            table_heatmap.loc[attr, cpg_type] += 1
        table_heatmap.loc["Total", cpg_type] += 1
    return table_heatmap

# %%
table_heatmap = get_cpg_attrib_heatmap(df_tot_target)
#%%
table_heatmap["Total"] = table_heatmap.sum(axis = 1)
table_heatmap_perc = table_heatmap/table_heatmap.loc["Total", "Total"] * 100
# table_heatmap_perc["Total"] = table_heatmap_perc.sum(axis = 1)
# %%
import seaborn as sns
from matplotlib import pyplot as plt
outdir = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/{visit}/{direction}"

# ax = sns.heatmap(data = table_heatmap_perc[list(filter(lambda x : x != "Total", table_heatmap_perc.columns))], annot = True, fmt = ".2f", cmap = "Oranges", cbar = False, vmin=0, vmax = 100, linewidths = 0.1)
ax = sns.heatmap(data = table_heatmap_perc[list(filter(lambda x : x != "Total", table_heatmap_perc.columns))], annot = True, fmt = ".2f", cmap = "rocket_r", cbar = False, vmin=0, vmax = 100, linewidths = 0.1)
for t in ax.texts:
    t.set_text(t.get_text() + " %")
plt.tight_layout()
plt.savefig(os.path.join(outdir, "annotation_heatmap.png"), dpi=600)
plt.show()
plt.close()

plt.figure(figsize=(3, 5))
plt.barh(width = table_heatmap_perc["Total"], y = table_heatmap_perc.index, color = "darkgray")
plt.gca().invert_yaxis()
plt.yticks([])
plt.xlabel("Percentage of CpG site involved (%)")
plt.tight_layout()
plt.savefig(os.path.join(outdir, "annotation_distibution.png"), dpi=600)
plt.show()
plt.close()

# %%
for symbol in df_tot_target["Gene_Symbol"].unique():
    if type(symbol) != float:
        print(symbol)

# %%
################################
## survey promoter methylation
################################
df_promot = df_tot[df_tot["Gene_Attrib"].isin(["Promoter", "Far_Promoter"])]
print('\n'.join(list(filter(lambda x : not x.startswith("ENST"), df_promot["Gene_Symbol"].unique()))))
enst_hypo = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231211/severe_mild_firstVisit/annotation/methyl_hypo_severe_mild_firstVisit_annotatr_promoter_only_enst_conv.tsv"
df_enst_hypo = pd.read_csv(enst_hypo, sep='\t')
print("\n".join(df_enst_hypo.loc[:, "external_gene_name"].dropna().unique()))

# %%
