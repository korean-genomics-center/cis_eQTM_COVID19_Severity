# %%
import glob
import math
import os
import warnings

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu, pearsonr
from sklearn.linear_model import LinearRegression
from statsmodels.stats.multitest import fdrcorrection

warnings.filterwarnings("ignore")

# %%
def main(path_exp, path_meta, path_crf, outdir, list_drop_samples, list_gene, colsample, colgene, colclin, corr_thres, pval_thres):
    df_TPM = read_TPM_matrix(path_exp, list_drop_samples)
    df_TPM = select_TPM_matrix(df_TPM, select_pattern="U", delim_id="-", namepos=0)
    df_TPM = select_TPM_matrix(df_TPM, select_pattern="K", delim_id="-", namepos=0)
    df_TPM = select_TPM_matrix(df_TPM, select_pattern="L", delim_id="-", namepos=2)
    df_TPM_gene_filtered = filter_list_genes(df_TPM, list_gene, colgene)
    df_TPM_transposed = transpose_TPM_matrix(df_TPM_gene_filtered, colsample)
    df_TPM_meta = merge_TPM_matrix_meta(df_TPM_transposed, path_meta, colsample)
    df_master_table = merge_TPM_matrix_CRF(df_TPM_meta, path_crf, colsample)
    find_genes_associated_crf(df_master_table, df_TPM_transposed, outdir, corr_thres, pval_thres, colgene, colclin)
    
    

def read_TPM_matrix(path_exp, list_drop_samples):
    df_TPM = pd.read_csv(path_exp, sep="\t")
    df_TPM  = df_TPM.drop(columns=list_drop_samples)

    return df_TPM

def select_TPM_matrix(df_TPM, select_pattern, delim_id="-", namepos=0):
    list_colname = list(df_TPM.columns)
    list_colname_filt = list(filter(lambda x: x.split(delim_id)[namepos][0] != select_pattern if len(x.split("-")) > 1 else x, list_colname))
    df_TPM = df_TPM.loc[:, list_colname_filt]

    return df_TPM

def filter_list_genes(df_TPM, list_gene, colgene="ID"):
    if len(list_gene) > 0:
        df_TPM["tmp"] = df_TPM[colgene].apply(lambda x: "_".join(x.split("_")[1:]))
        df_TPM_gene_filtered = df_TPM[df_TPM["tmp"].isin(list_gene)]
        df_TPM_gene_filtered = df_TPM_gene_filtered.drop(columns=["tmp"])
    else:
        df_TPM_gene_filtered = df_TPM
        
    return df_TPM_gene_filtered

def transpose_TPM_matrix(df_TPM, colsample):
    list_geneid = df_TPM[colsample].to_list()
    df_TPM_transposed = df_TPM.T.iloc[1:, :]
    df_TPM_transposed.columns = list_geneid
    df_TPM_transposed = df_TPM_transposed.reset_index(drop=False).rename(columns={"index": colsample})
    
    return df_TPM_transposed

def merge_TPM_matrix_meta(df_TPM, path_meta, colsample):
    df_meta = pd.read_csv(path_meta, sep="\t")
    df_meta = df_meta.rename(columns={"Project_ID_Alias": colsample})
    df_merged = pd.merge(df_TPM, df_meta, how="inner", on=colsample)
    
    return df_merged

def remove_useless_spaces_crf(list_crf_names):
    list_new_crf_names = list(map(lambda x: x.replace("\n", "").replace("\t", "_").replace("\xa0", "_").replace("\xa1", "_").replace("\xa2", "_").replace(" ", "_"), list_crf_names))

    return list_new_crf_names

def merge_TPM_matrix_CRF(df_TPM, path_crf, colsample, num=2):
    df_crf_raw = pd.read_excel(path_crf, engine="openpyxl", sheet_name="21-22등록 대상자_modify", skiprows=1)
    df_crf_select = df_crf_raw.iloc[:, num:]
    df_crf = df_crf_select.rename(columns={"Subject NO.(고유번호)": colsample})
    list_colname = list(df_crf.columns)
    list_new_colname = remove_useless_spaces_crf(list_colname)
    df_crf.columns = list_new_colname
    df_merged = pd.merge(df_TPM, df_crf, how="inner", on=colsample) 
    
    return df_merged

def filter_list_elem_float_only(df_master_table, list_all):
    list_filtered = list(filter(lambda x: isinstance(df_master_table[x].to_list()[0], float) and df_master_table[x].to_list()[0]!=np.nan, list_all))

    return list_filtered

def get_dataframe_gene_exp_clinical_value(df_master_table, gene, value):
    # select only a specific gene and clinical value
    df_gene_value = df_master_table[[gene, value]]
    # remove non-numeric characters within a float value
    df_gene_value = df_gene_value.replace(r'[^0-9]+', '', regex=True)
    # remove empty characters with na
    df_gene_value = df_gene_value.replace(r'^\s*$', np.nan, regex=True)
    # drop na row/sample-wise
    df_gene_value = df_gene_value.dropna(axis=0)

    return df_gene_value

def find_genes_associated_crf(df_master_table, df_TPM_transposed, outdir, corr_thres, pval_thres, colgene, colclin):
    gene_table = df_master_table.iloc[:, 1:len(df_TPM_transposed.columns)]
    value_table = df_master_table.iloc[:, len(df_TPM_transposed.columns):]
    list_gene_comp = list(gene_table.columns)
    list_value_comp = list(value_table.columns)
    list_value_comp_filt = filter_list_elem_float_only(df_master_table, list_value_comp)
    
    dict_gene_value_mann = dict()
    dict_gene_value_corr = dict()
    for gene in list_gene_comp:

        dict_value_mann = dict()
        dict_value_corr = dict()
        for value in list_value_comp_filt:
            df_gene_value = get_dataframe_gene_exp_clinical_value(df_master_table, gene, value)
            array_exp = df_gene_value[gene].astype(float).values
            array_value = df_gene_value[value].astype(float).values

            # remove any array(s) with only one value
            if len(np.unique(array_value)) <= 1:
                continue

            # compare two mannwhitney
            if len(np.unique(array_value)) == 2:
                list_list_cat = list()
                list_list_exp = list()
                grby_val = df_gene_value.groupby(value)[gene].apply(list)
                grpby_cnt = df_gene_value.groupby(value)[gene].apply(len)
                grpby_idx = list(grby_val.index)

                for cat, cnt in zip(grpby_idx, grpby_cnt):
                    list_list_cat.append(f"{cat}\nN={cnt}")

                for val in grby_val:
                    list_list_exp.append(val)

                list_control = list_list_exp[0]
                list_case = list_list_exp[1]

                if len(set(list_control)) > 1:
                    stat, pval = mannwhitneyu(list_case, list_control)
                    stat_res = {"stat": stat, "pval": pval, "categories": list_list_cat, "expression": list_list_exp}
                    dict_value_mann[value] = stat_res
                            
            # compare multiple correlation
            if len(np.unique(array_value)) > 2:
                # get array of expression and values
                corr, pval = pearsonr(array_value, array_exp)
                corr_res = {"corr": corr, "pval": pval, "clinical": array_value, "expression": array_exp}
                dict_value_corr[value] = corr_res

        # get gene-clinical compare 
        if dict_gene_value_mann.get(gene) == None:
            dict_gene_value_mann[gene] = dict()
        dict_gene_value_mann[gene].update(dict_value_mann)

        # get gene-clinical correlation 
        if dict_gene_value_corr.get(gene) == None:
            dict_gene_value_corr[gene] = dict()
        dict_gene_value_corr[gene].update(dict_value_corr)

    df_stats = pd.concat({k: pd.DataFrame(v).T for k, v in dict_gene_value_mann.items()}, axis=0)
    df_stats = df_stats.reset_index(drop=False).rename(columns={"level_0": colgene, "level_1": colclin})
    df_stats_sorted = df_stats.sort_values(by=["pval"], ascending=True)
    df_stats_sorted.to_csv(os.path.join(outdir, "compare_statistics.tsv"), sep="\t", index=False)

    df_corr = pd.concat({k: pd.DataFrame(v).T for k, v in dict_gene_value_corr.items()}, axis=0)
    df_corr = df_corr.reset_index(drop=False).rename(columns={"level_0": colgene, "level_1": colclin})
    df_corr["abs_corr"] = df_corr["corr"].apply(abs)
    df_corr_sorted = df_corr.sort_values(by=["abs_corr"], ascending=False)
    df_corr_sorted.to_csv(os.path.join(outdir, "correlation_statistics.tsv"), sep="\t", index=False)

    return df_stats_sorted, df_corr_sorted

def set_korean_font(fontpath = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Resources/Data/NanumGothicBold.ttf"):
    fontname = os.path.basename(fontpath).replace(".ttf", "")
    fonttype = fm.FontEntry(fname=fontpath, name=fontname)
    fm.fontManager.ttflist.insert(0, fonttype)
    plt.rcParams.update({"font.size": 10, "font.family":fontname})

def plot_cat(value, gene, list_list_cat, list_list_exp, stat, pval, outdir):
    plt.boxplot(list_list_exp)
    plt.xticks(ticks=list(range(1, len(list_list_cat)+1, 1)), labels=list_list_cat)
    plt.xlabel(value)
    plt.ylabel(gene)
    # plt.annotate(f"$u$:{round(float(stat),3)}\n$p$:{round(float(pval),3)}", xy=(0.10, 0.70), xycoords='axes fraction')
    outfigpath = os.path.join(outdir, value)
    os.makedirs(outfigpath, exist_ok=True)
    set_korean_font()
    plt.savefig(f"{outfigpath}/{gene}_{value}_barplot.png", dpi=100)
    plt.show()
    plt.close()

def plot_num(value, gene, array_value, array_exp, corr, pval, outdir):
    x_reshape = array_value.reshape(-1,1)
    y = array_exp
    model = LinearRegression().fit(x_reshape, y)
    array_exp_pred = model.predict(x_reshape)
    rsq = model.score(x_reshape, y)
    plt.scatter(array_value, array_exp)
    plt.plot(array_value, array_exp_pred)
    plt.xlabel(value)
    plt.ylabel(gene)
    plt.annotate(f"$R^2$:{round(float(rsq),3)}\n$r$:{round(float(corr),3)}\n$p$:{round(float(pval),3)}", xy=(0.10, 0.70), xycoords='axes fraction')
    outfigpath = os.path.join(outdir, value)
    os.makedirs(outfigpath, exist_ok=True)
    set_korean_font()
    plt.savefig(f"{outfigpath}/{gene}_{value}_linearplot.png", dpi=100)
    plt.show()
    plt.close()

# %%
if __name__ == "__main__":
    colsample = "ID"
    colgene = "ID"
    colclin = "Clinical"
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
    path_exp = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/4_expmtx/Infectomics/expression_matrix_genes.results_TPM.tsv"
    path_crf = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Resources/Data/infectomics_CRF_20230410_edit.xlsx"
    path_meta = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Resources/Data/COVID19_master_table_added_CRF_20231007.txt"
    outdir = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/9_clinical"
    list_gene = ["PRRC2A"]
    # list_gene = ['ENSG00000284820', 'ENSG00000290242', 'ENSG00000254680', 'RPIA', 'GABARAPL2', 'STIM2-AS1', 'METTL7B', 'MARCHF8', 'RNF123', 'GPX1', 'NCOA4', 'SMIM24', 'ANKH', 'GATA1', 'ST6GALNAC4', 'TNXB', 'CLRN1-AS1', 'NEU1', 'AOAH-IT1', 'SAMD14', 'ENSG00000257607', 'LINC02207', 'FKBP8', 'ENSG00000288524', 'ABALON', 'NFIX', 'RNF10', 'ENSG00000287166', 'VNN1', 'FAM104A', 'TNS1', 'MBNL3', 'ACHE', 'CREG1', 'YPEL4', 'PHOSPHO1', 'FAXDC2', 'FBXO7', 'MAOB', 'BBOF1', 'RBM38', 'TMC5', 'FOXO4', 'ARL4A', 'ENSG00000275719', 'TPGS2', 'PIP5K1B', 'IGHV3-13', 'LINC02967', 'HAGH', 'MSH5-SAPCD1', 'H1-2', 'ENSG00000268119', 'SLCO2B1', 'ENSG00000259033', 'DMTN', 'YBX3P1', 'ENSG00000239920', 'SLC6A8', 'TFDP1', 'RIOK3', 'ENSG00000289514', 'EYA2', 'ANK1', 'ENSG00000282024', 'YBX3', 'ENSG00000287642', 'C4B', 'LINC01036', 'TSPAN7', 'CLIC2', 'ADIPOR1', 'KANK2', 'RNF11', 'C17orf99', 'SIAH2', 'ENSG00000253986', 'TBCEL-TECTA', 'ACKR1', 'ENSG00000289377', 'LRRC2', 'ENSG00000291105', 'ENSG00000267952', 'SLC6A9', 'SHISA7', 'TENT5C', 'VSTM2L', 'MYO18B', 'BCL2L1', 'OR2W3', 'CR1L', 'KEL', 'ENSG00000267780', 'ENSG00000287632', 'ENSG00000246203', 'ENSG00000239219', 'DNAJC6', 'PLEK2', 'KDM7A-DT', 'TSPO2', 'EPB42', 'MFSD2B', 'GLRX5', 'GMPR', 'ENSG00000235105', 'PAQR9', 'SLC4A1', 'TMCC2', 'NF1P8', 'ENSG00000265401', 'TRHDE-AS1', 'TGM2', 'IGHV1-3', 'OSBP2', 'GYPE', 'C4A', 'MXI1', 'KRT1', 'GFUS', 'SOX6', 'ABCG2', 'ZDHHC19', 'ITLN1', 'FECH', 'IFI27', 'SNCA', 'SGIP1', 'ENSG00000287860', 'LIFR', 'RHCE', 'C9orf153', 'WFDC1', 'CTNNAL1', 'ABCC13', 'OR10Z1', 'HEPACAM2', 'ENSG00000289933', 'ENSG00000224091', 'PAGE2B', 'ENSG00000264066', 'SELENBP1', 'ALAS2', 'ART4', 'RELN', 'HBG2', 'HEMGN', 'THEM5', 'ENSG00000290318', 'ENSG00000285498', 'GYPB', 'ENSG00000255477', 'ENSG00000287510', 'ENSG00000285117', 'BPGM', 'ESRG', 'RHAG', 'TMEM158', 'ABCB5', 'FAM83A', 'SLC6A19', 'HLA-C', 'RING1', 'UICLM', 'CCL2', 'IGHV7-4-1', 'ENSG00000282678', 'RPL19P16', 'ENSG00000274422', 'ENSG00000289278', 'SIGLEC8', 'PRSS33', 'MOAP1', 'IGHV3-11', 'CTSG', 'LILRB5', 'MAOB', 'SAMD14', 'C19orf33', 'GJB6', 'APOE', 'KRT23', 'SULT1A4', 'ENSG00000232220', 'DEFA3', 'ENSG00000265401', 'ENSG00000288853', 'CA1', 'CD177', 'ENSG00000290242', 'ORM2', 'OLAH', 'KCNMA1', 'WFDC1', 'DEFA1', 'ENSG00000287510','ENSG00000230699', 'HOXC5', 'LY6G5B', 'MSH5', 'HLA-E', 'ENSG00000285314', 'NHIP', 'PRPF31-AS1', 'ROBO1', 'KRT5', 'KLRC2', 'ALPK2', 'DMRTC1', 'ENSG00000282484', 'RASGRF2-AS1', 'MYZAP', 'MICU3', 'IGKV2-29', 'ENSG00000287277', 'NAT8L', 'IGHV1-2','IGHV3-74', 'MSH5-SAPCD1', 'ENSG00000288681', 'HLA-DQB1', 'HLA-H', 'MTND1P23', 'ATF6B','CSNK2B', 'ZNF780B', 'IGHV3-7']
    # list_gene = ['HLA-C', 'RING1', 'UICLM', 'CCL2', 'IGHV7-4-1', 'ENSG00000282678', 'RPL19P16', 'ENSG00000274422', 'ENSG00000289278', 'SIGLEC8', 'PRSS33', 'MOAP1', 'IGHV3-11', 'ENSG00000284820', 'ENSG00000290242', 'ENSG00000254680', 'RPIA', 'GABARAPL2', 'STIM2-AS1', 'METTL7B', 'MARCHF8', 'RNF123', 'GPX1', 'NCOA4', 'SMIM24', 'ANKH', 'GATA1', 'ST6GALNAC4', 'TNXB', 'CLRN1-AS1', 'NEU1', 'AOAH-IT1', 'SAMD14', 'ENSG00000257607', 'LINC02207', 'FKBP8', 'ENSG00000288524', 'ABALON', 'NFIX', 'RNF10', 'ENSG00000287166', 'VNN1', 'FAM104A', 'TNS1', 'MBNL3', 'ACHE', 'CREG1', 'YPEL4', 'PHOSPHO1', 'FAXDC2', 'FBXO7', 'MAOB', 'BBOF1', 'RBM38', 'TMC5', 'FOXO4', 'ARL4A', 'ENSG00000275719', 'TPGS2', 'PIP5K1B', 'IGHV3-13', 'LINC02967', 'HAGH', 'MSH5-SAPCD1', 'H1-2', 'ENSG00000268119', 'SLCO2B1', 'ENSG00000259033', 'DMTN', 'YBX3P1', 'ENSG00000239920', 'SLC6A8', 'TFDP1', 'RIOK3', 'ENSG00000289514', 'EYA2', 'ANK1', 'ENSG00000282024', 'YBX3', 'ENSG00000287642', 'C4B', 'LINC01036', 'TSPAN7', 'CLIC2', 'ADIPOR1', 'KANK2', 'RNF11', 'C17orf99', 'SIAH2', 'ENSG00000253986', 'TBCEL-TECTA', 'ACKR1', 'ENSG00000289377', 'LRRC2', 'ENSG00000291105', 'ENSG00000267952', 'SLC6A9', 'SHISA7', 'TENT5C', 'VSTM2L', 'MYO18B', 'BCL2L1', 'OR2W3', 'CR1L', 'KEL', 'ENSG00000267780', 'ENSG00000287632', 'ENSG00000246203', 'ENSG00000239219', 'DNAJC6', 'PLEK2', 'KDM7A-DT', 'TSPO2', 'EPB42', 'MFSD2B', 'GLRX5', 'GMPR', 'ENSG00000235105', 'PAQR9', 'SLC4A1', 'TMCC2', 'NF1P8', 'ENSG00000265401', 'TRHDE-AS1', 'TGM2', 'IGHV1-3', 'OSBP2', 'GYPE', 'C4A', 'MXI1', 'KRT1', 'GFUS', 'SOX6', 'ABCG2', 'ZDHHC19', 'ITLN1', 'FECH', 'IFI27', 'SNCA', 'SGIP1', 'ENSG00000287860', 'LIFR', 'RHCE', 'C9orf153', 'WFDC1', 'CTNNAL1', 'ABCC13', 'OR10Z1', 'HEPACAM2', 'ENSG00000289933', 'ENSG00000224091', 'PAGE2B', 'ENSG00000264066', 'SELENBP1', 'ALAS2', 'ART4', 'RELN', 'HBG2', 'HEMGN', 'THEM5', 'ENSG00000290318', 'ENSG00000285498', 'GYPB', 'ENSG00000255477', 'ENSG00000287510', 'ENSG00000285117', 'BPGM', 'ESRG', 'RHAG', 'TMEM158', 'ABCB5', 'FAM83A', 'SLC6A19']
    corr_thres=0.5
    pval_thres=0.05
    # main(path_exp, path_meta, path_crf, outdir, list_drop_samples, list_gene, colsample, colgene, colclin, corr_thres, pval_thres)

# %%