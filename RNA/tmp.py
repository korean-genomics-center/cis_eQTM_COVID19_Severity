#%%
import pandas as pd
import os
from multiprocessing import Pool

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

# %%
def get_annot_result(path_annotatr_result):
    table_annot = pd.read_csv(path_annotatr_result, sep = '\t')
    table_annot["position"] = table_annot.apply(lambda x : f"{x['seqnames']}_{int(x['start'])-1}_{int(x['end'])-1}", axis = 1)
    table_annot = table_annot.set_index("position", drop = False)
    list_chr = table_annot["seqnames"].unique()
    list_table_chr = list(map(lambda x : table_annot[table_annot["seqnames"] == x], list_chr))
    with Pool(processes=len(list_chr)) as p:
        list_table_organized = p.map(organize_table_by_chr, list_table_chr)
    df_tot = pd.concat(list_table_organized)

    return df_tot

visit = "first"
path_annotatr_result_hyper = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231211/severe_mild_{visit}Visit/annotation/methyl_hyper_severe_mild_{visit}Visit_annotatr.tsv"
path_annotatr_result_hypo = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231211/severe_mild_{visit}Visit/annotation/methyl_hypo_severe_mild_{visit}Visit_annotatr.tsv"
df_annot_hyper = get_annot_result(path_annotatr_result_hyper)
df_annot_hypo = get_annot_result(path_annotatr_result_hypo)
df_annot = pd.concat([df_annot_hyper, df_annot_hypo], axis=0)

# %%
file_overlap_dmp_deg = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/DMPDEG/first/table_deg_dmp_overlap_abslog2fc_1.3_qval_0.05_20240229_sorted.tsv"
df_overlap = pd.read_csv(file_overlap_dmp_deg,sep="\t")
list_overlap_dmp = df_overlap["DMP_Markers"].to_list()
list_overlap_gene = list(map(lambda x: "_".join(x.split("_")[1:]), df_overlap["Gene_Symbol"].to_list()))

list_dmp_deg_overlap_dmp_all = list()
for overlap_dmp in list_overlap_dmp:
    marker = overlap_dmp.split("/")
    list_dmp_deg_overlap_dmp_all.extend(marker)

list_dmp_deg_overlap_gene_all = list()
for overlap_gene in list_overlap_gene:
    gene = overlap_gene.split("/")
    list_dmp_deg_overlap_gene_all.extend(gene)

list_dmp_deg_overlap_all_posname = list(map(lambda x: x.split(":")[0] + "_" + x.split(":")[1] + "_" + x.split(":")[1], list_dmp_deg_overlap_dmp_all))
df_annot_degdmp_overlap_marker = df_annot[df_annot["PosName"].isin(list_dmp_deg_overlap_all_posname)]
df_annot_degdmp_overlap_gene = df_annot_degdmp_overlap_marker[df_annot_degdmp_overlap_marker["Gene_Symbol"].isin(list_dmp_deg_overlap_gene_all)]
df_annot_degdmp_attrib = df_annot_degdmp_overlap_gene.groupby(["PosName", "Gene_Symbol"])["Gene_Attrib"].apply(lambda x: "/".join(list(set(x)))).reset_index(drop=False)
df_annot_degdmp_cpgs = df_annot_degdmp_overlap_gene.groupby(["PosName", "Gene_Symbol"])["CpGs"].apply(lambda x: "/".join(list(set(x)))).reset_index(drop=False)
list_cpgs = df_annot_degdmp_cpgs["CpGs"].to_list()
df_annot_degdmp_attrib["CpGs"] = list_cpgs
df_annot_degdmp_attrib["PosName"] = df_annot_degdmp_attrib["PosName"].apply(lambda x: x.split("_")[0] + ":" + x.split("_")[1])
df_annot_degdmp_attrib.columns = ["Markers", "Gene_Symbol", "Location", "CpG_Site"]
df_annot_degdmp_attrib.to_csv("/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/DMPDEG/first/table_deg_dmp_overlap_abslog2fc_1.3_qval_0.05_20240229_sorted_annotated.tsv", sep="\t", index=False)
# %%
