# %%
import math
import os

import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from sklearn.preprocessing import StandardScaler
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import cm
from matplotlib.colors import rgb2hex
from matplotlib.patches import Ellipse
from sklearn.decomposition import PCA

# %%
def main(expmtx, metadata, list_gene, list_drop_samples, colsample, n_components, outdir, pcx_num, pcy_num, dict_order):
    pcx = f"PC{pcx_num}"
    pcy = f"PC{pcy_num}"

    df_TPM = read_TPM_matrix(expmtx, list_drop_samples)
    df_TPM_gene_selected = filter_list_genes(df_TPM, list_gene , mode="gene")
    # list_colgenes = list(df_TPM_gene_selected["ID"])
    df_TPM_transposed = transpose_TPM_matrix(df_TPM_gene_selected, colsample)
    scaler = StandardScaler().fit(df_TPM_transposed)
    tpm_scaled = scaler.transform(df_TPM_transposed)
    pca = PCA(n_components=n_components).fit(tpm_scaled)
    pca_tpm = pca.transform(tpm_scaled)

    # # Draw scree plot
    # plt.plot(np.cumsum(pca.explained_variance_ratio_))
    # plt.xlabel('number of components')filter_list_genes
    # plt.ylabel('cumulative explained variance')
    # plt.show()
    # plt.close()

    list_PC = list(map(lambda x: f"PC{x}", list(range(1, n_components+1, 1))))
    list_samples = list(df_TPM_transposed.index)
    df_pca_tpm = pd.DataFrame(pca_tpm, columns=list_PC, index=list_samples).reset_index(drop=False).rename(columns={"index":colsample})
    df_meta = read_metadata(metadata, colsample)
    df_meta = df_meta[df_meta["Sample_Age"] != "-"]
    df_meta["Sample_Age_Group"] = df_meta["Sample_Age"].apply(lambda x: "Old" if int(x) > 60 else "Young")
    df_meta["Severity_Visit"] = df_meta["Severity_Binary"] + "_" + df_meta["Visit_order"]
    list_meta = list(df_meta.columns)[2:]
    df_pca_meta = pd.merge(df_pca_tpm, df_meta, how="inner", on=colsample)
    df_pca_meta_filt = df_pca_meta[df_pca_meta["Visit_order"]!="Control"]
    df_pca_meta_filt = df_pca_meta_filt[df_pca_meta_filt["Visit_order"]!="LongCOVID"]
    df_pca_meta_filt["Severity"] = df_pca_meta_filt["Severity"].apply(lambda x: "Convalescent" if x == "-" else x)
    df_pca_meta_filt["Severity_Binary"] = df_pca_meta_filt["Severity_Binary"].apply(lambda x: "Convalescent" if x == "-" else x)
    df_pca_meta_filt["Severity_Visit"] = df_pca_meta_filt["Severity_Binary"] + "_" + df_pca_meta_filt["Visit_order"]
    
    # # Extract genes involved in making components
    # # number of components
    # n_pcs= pca.components_.shape[0]
    # # get the index of the most important feature on EACH component
    # most_important = [np.abs(pca.components_[i]).argsort()[::-1][:5] for i in range(n_pcs)]
    # initial_feature_names = list_colgenes
    # # get the names
    # most_important_names = [initial_feature_names[most_important[i]] for i in range(n_pcs)]
    # dic = {'PC{}'.format(i): most_important_names[i] for i in range(n_pcs)}
    # # build the dataframe
    # df = pd.DataFrame(dic.items())

    # df_pca_meta_filt = df_pca_meta_filt[~df_pca_meta_filt["ID"].isin(list_outlier)]
    outfilename = os.path.join(outdir, "pca_meta_dataframe.tsv")
    df_pca_meta_filt.to_csv(outfilename, sep="\t", index=False)

    fig, axes = plt.subplots(nrows=int(math.sqrt(len(list_meta))+1), ncols=int(math.sqrt(len(list_meta))+1), figsize=(30, 30))
    for meta, ax in zip(list_meta, axes.flatten()):
        if meta == "Sample_Age":
            continue

        val_groupby_meta_pcx = df_pca_meta_filt.groupby(meta)[pcx].apply(np.array).tolist()
        val_groupby_meta_pcy = df_pca_meta_filt.groupby(meta)[pcy].apply(np.array).tolist()
        cmap = cm.get_cmap('Spectral', len(val_groupby_meta_pcx))
        list_colors = [rgb2hex(cmap(ind)) for ind in range(len(val_groupby_meta_pcx))]
        list_legends = list(df_pca_meta_filt.groupby(meta)[pcy].apply(np.array).index)

        for x, y, color, legend in zip(val_groupby_meta_pcx, val_groupby_meta_pcy, list_colors, list_legends):
            # ax.scatter(x, y, c=color, s=0, label=list_legends)
            confidence_ellipse(x, y, ax, n_std=3.0, edgecolor=color, facecolor='none')
            ax.scatter(np.mean(x), np.mean(y), c=color, s=20, label=legend)
            ax.legend(loc='best')  
            ax.set_title(meta, fontsize=20)
            ax.set_ylabel(f"{pcy}\nVariance_Explained: {round(pca.explained_variance_ratio_[pcy_num-1]*100, 1)}%", fontsize=20)
            ax.set_xlabel(f"{pcx}\nVariance_Explained: {round(pca.explained_variance_ratio_[pcx_num-1]*100, 1)}%", fontsize=20)

        if meta in dict_order.keys():
            list_meta_order = dict_order.get(meta)
            for meta_order in list_meta_order:
                zip_draw_data = zip(val_groupby_meta_pcx, val_groupby_meta_pcy, list_colors, list_legends)
                sorted_zip_draw_data = sorted(zip_draw_data, key = lambda x : meta_order.get(x[3], math.inf))
                for ind_1 in range(len(meta_order)-1):
                    ind_2 = ind_1 + 1
                    if meta_order.get(sorted_zip_draw_data[ind_2][3]) == None: break

                    prior_x, prior_y, _, _ = sorted_zip_draw_data[ind_1]
                    post_x, post_y, _, _ = sorted_zip_draw_data[ind_2]

                    mean_prior_x = np.mean(prior_x)
                    mean_prior_y = np.mean(prior_y)
                    mean_post_x = np.mean(post_x)
                    mean_post_y = np.mean(post_y)

                    prop = dict(arrowstyle="-|>,head_width=0.1,head_length=0.5", connectionstyle="arc3,rad=.5", shrinkA=0, shrinkB=0, color="k")
                    ax.annotate("", xytext=(mean_prior_x,mean_prior_y), xy=((mean_post_x),(mean_post_y)), arrowprops=prop)

    fig.tight_layout()
    outfigname = os.path.join(outdir, "pca.png")
    plt.savefig(outfigname, dpi=300)
    plt.show()
    plt.close()

    # x = df_pca_meta_filt[pcx].to_list()
    # y = df_pca_meta_filt[pcy].to_list()
    # top_x = sorted(df_pca_meta_filt[pcx].to_list())[-5:]
    # top_y = sorted(df_pca_meta_filt[pcy].to_list())[-5:]

    # list_outlier = list()
    # for sample, pc2 in zip(df_pca_meta_filt["ID"], df_pca_meta_filt[pcx]):
    #     if pc2 in top_x:
    #         list_outlier.append(sample)
    # for sample, pc1 in zip(df_pca_meta_filt["ID"], df_pca_meta_filt[pcy]):
    #     if pc1 in top_y:
    #         list_outlier.append(sample)
    # list_outlier = sorted(list_outlier)

    # annot_x = df_pca_meta_filt[df_pca_meta_filt["ID"].isin(list_outlier)][pcx].to_list()
    # annot_y = df_pca_meta_filt[df_pca_meta_filt["ID"].isin(list_outlier)][pcy].to_list()

    # plt.ylabel(fpcynVariance_Explained: {round(pca.explained_variance_ratio_[0]*100,1)}%")
    # plt.xlabel(fpcxnVariance_Explained: {round(pca.explained_variance_ratio_[1]*100,1)}%")
    # plt.title(meta)
    # for i in range(len(list_outlier)):
    #     plt.annotate(list_outlier[i], (annot_x[i], annot_y[i] + 0.2))
    # plt.legend(loc = "upper left", bbox_to_anchor = (1.05, 0.95))
    # plt.show()
    # plt.close()

def read_TPM_matrix(expmtx, list_drop_samples):
    df_TPM = pd.read_csv(expmtx, sep="\t")
    df_TPM = df_TPM.drop(columns=list_drop_samples)

    return df_TPM

def filter_list_genes(df_TPM, list_gene, mode, colgene="ID"):
    if len(list_gene) > 0:
        if mode == "genesymbol":
            df_TPM["tmp"] = df_TPM[colgene].apply(lambda x: "_".join(x.split("_")[1:]))
            df_TPM_gene_filtered = df_TPM[df_TPM["tmp"].isin(list_gene)].drop(columns=["tmp"])
        if mode == "geneonly":
            df_TPM_gene_filtered = df_TPM[df_TPM[colgene].isin(list_gene)].drop(columns=["tmp"])

    else:
        df_TPM_gene_filtered = df_TPM
        
    return df_TPM_gene_filtered

def transpose_TPM_matrix(df_TPM, colsample):
    list_geneid = df_TPM[colsample].to_list()
    df_TPM_transposed = df_TPM.T.iloc[1:, :]
    df_TPM_transposed.columns = list_geneid
    
    return df_TPM_transposed

def read_metadata(metadata, colsample):
    df_meta = pd.read_csv(metadata, sep="\t")
    df_meta = df_meta.rename(columns={"Project_ID_Alias": colsample})

    return df_meta

def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensional dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the standard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the standard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)
# %%
if __name__ == "__main__":
    expmtx = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/4_expmtx/Infectomics/expression_matrix_genes.results_TPM.tsv"
    metadata = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Resources/Data/COVID19_master_table_added_CRF_20231007.txt"
    # list_gene = ['ENSG00000284820', 'ENSG00000290242', 'ENSG00000254680', 'RPIA', 'GABARAPL2', 'STIM2-AS1', 'METTL7B', 'MARCHF8', 'RNF123', 'GPX1', 'NCOA4', 'SMIM24', 'ANKH', 'GATA1', 'ST6GALNAC4', 'TNXB', 'CLRN1-AS1', 'NEU1', 'AOAH-IT1', 'SAMD14', 'ENSG00000257607', 'LINC02207', 'FKBP8', 'ENSG00000288524', 'ABALON', 'NFIX', 'RNF10', 'ENSG00000287166', 'VNN1', 'FAM104A', 'TNS1', 'MBNL3', 'ACHE', 'CREG1', 'YPEL4', 'PHOSPHO1', 'FAXDC2', 'FBXO7', 'MAOB', 'BBOF1', 'RBM38', 'TMC5', 'FOXO4', 'ARL4A', 'ENSG00000275719', 'TPGS2', 'PIP5K1B', 'IGHV3-13', 'LINC02967', 'HAGH', 'MSH5-SAPCD1', 'H1-2', 'ENSG00000268119', 'SLCO2B1', 'ENSG00000259033', 'DMTN', 'YBX3P1', 'ENSG00000239920', 'SLC6A8', 'TFDP1', 'RIOK3', 'ENSG00000289514', 'EYA2', 'ANK1', 'ENSG00000282024', 'YBX3', 'ENSG00000287642', 'C4B', 'LINC01036', 'TSPAN7', 'CLIC2', 'ADIPOR1', 'KANK2', 'RNF11', 'C17orf99', 'SIAH2', 'ENSG00000253986', 'TBCEL-TECTA', 'ACKR1', 'ENSG00000289377', 'LRRC2', 'ENSG00000291105', 'ENSG00000267952', 'SLC6A9', 'SHISA7', 'TENT5C', 'VSTM2L', 'MYO18B', 'BCL2L1', 'OR2W3', 'CR1L', 'KEL', 'ENSG00000267780', 'ENSG00000287632', 'ENSG00000246203', 'ENSG00000239219', 'DNAJC6', 'PLEK2', 'KDM7A-DT', 'TSPO2', 'EPB42', 'MFSD2B', 'GLRX5', 'GMPR', 'ENSG00000235105', 'PAQR9', 'SLC4A1', 'TMCC2', 'NF1P8', 'ENSG00000265401', 'TRHDE-AS1', 'TGM2', 'IGHV1-3', 'OSBP2', 'GYPE', 'C4A', 'MXI1', 'KRT1', 'GFUS', 'SOX6', 'ABCG2', 'ZDHHC19', 'ITLN1', 'FECH', 'IFI27', 'SNCA', 'SGIP1', 'ENSG00000287860', 'LIFR', 'RHCE', 'C9orf153', 'WFDC1', 'CTNNAL1', 'ABCC13', 'OR10Z1', 'HEPACAM2', 'ENSG00000289933', 'ENSG00000224091', 'PAGE2B', 'ENSG00000264066', 'SELENBP1', 'ALAS2', 'ART4', 'RELN', 'HBG2', 'HEMGN', 'THEM5', 'ENSG00000290318', 'ENSG00000285498', 'GYPB', 'ENSG00000255477', 'ENSG00000287510', 'ENSG00000285117', 'BPGM', 'ESRG', 'RHAG', 'TMEM158', 'ABCB5', 'FAM83A', 'SLC6A19', 'HLA-C', 'RING1', 'UICLM', 'CCL2', 'IGHV7-4-1', 'ENSG00000282678', 'RPL19P16', 'ENSG00000274422', 'ENSG00000289278', 'SIGLEC8', 'PRSS33', 'MOAP1', 'IGHV3-11', 'CTSG', 'LILRB5', 'MAOB', 'SAMD14', 'C19orf33', 'GJB6', 'APOE', 'KRT23', 'SULT1A4', 'ENSG00000232220', 'DEFA3', 'ENSG00000265401', 'ENSG00000288853', 'CA1', 'CD177', 'ENSG00000290242', 'ORM2', 'OLAH', 'KCNMA1', 'WFDC1', 'DEFA1', 'ENSG00000287510','ENSG00000230699', 'HOXC5', 'LY6G5B', 'MSH5', 'HLA-E', 'ENSG00000285314', 'NHIP', 'PRPF31-AS1', 'ROBO1', 'KRT5', 'KLRC2', 'ALPK2', 'DMRTC1', 'ENSG00000282484', 'RASGRF2-AS1', 'MYZAP', 'MICU3', 'IGKV2-29', 'ENSG00000287277', 'NAT8L', 'IGHV1-2','IGHV3-74', 'MSH5-SAPCD1', 'ENSG00000288681', 'HLA-DQB1', 'HLA-H', 'MTND1P23', 'ATF6B','CSNK2B', 'ZNF780B', 'IGHV3-7']
    # list_gene = ['HLA-C', 'RING1', 'UICLM', 'CCL2', 'IGHV7-4-1', 'ENSG00000282678', 'RPL19P16', 'ENSG00000274422', 'ENSG00000289278', 'SIGLEC8', 'PRSS33', 'MOAP1', 'IGHV3-11', 'ENSG00000284820', 'ENSG00000290242', 'ENSG00000254680', 'RPIA', 'GABARAPL2', 'STIM2-AS1', 'METTL7B', 'MARCHF8', 'RNF123', 'GPX1', 'NCOA4', 'SMIM24', 'ANKH', 'GATA1', 'ST6GALNAC4', 'TNXB', 'CLRN1-AS1', 'NEU1', 'AOAH-IT1', 'SAMD14', 'ENSG00000257607', 'LINC02207', 'FKBP8', 'ENSG00000288524', 'ABALON', 'NFIX', 'RNF10', 'ENSG00000287166', 'VNN1', 'FAM104A', 'TNS1', 'MBNL3', 'ACHE', 'CREG1', 'YPEL4', 'PHOSPHO1', 'FAXDC2', 'FBXO7', 'MAOB', 'BBOF1', 'RBM38', 'TMC5', 'FOXO4', 'ARL4A', 'ENSG00000275719', 'TPGS2', 'PIP5K1B', 'IGHV3-13', 'LINC02967', 'HAGH', 'MSH5-SAPCD1', 'H1-2', 'ENSG00000268119', 'SLCO2B1', 'ENSG00000259033', 'DMTN', 'YBX3P1', 'ENSG00000239920', 'SLC6A8', 'TFDP1', 'RIOK3', 'ENSG00000289514', 'EYA2', 'ANK1', 'ENSG00000282024', 'YBX3', 'ENSG00000287642', 'C4B', 'LINC01036', 'TSPAN7', 'CLIC2', 'ADIPOR1', 'KANK2', 'RNF11', 'C17orf99', 'SIAH2', 'ENSG00000253986', 'TBCEL-TECTA', 'ACKR1', 'ENSG00000289377', 'LRRC2', 'ENSG00000291105', 'ENSG00000267952', 'SLC6A9', 'SHISA7', 'TENT5C', 'VSTM2L', 'MYO18B', 'BCL2L1', 'OR2W3', 'CR1L', 'KEL', 'ENSG00000267780', 'ENSG00000287632', 'ENSG00000246203', 'ENSG00000239219', 'DNAJC6', 'PLEK2', 'KDM7A-DT', 'TSPO2', 'EPB42', 'MFSD2B', 'GLRX5', 'GMPR', 'ENSG00000235105', 'PAQR9', 'SLC4A1', 'TMCC2', 'NF1P8', 'ENSG00000265401', 'TRHDE-AS1', 'TGM2', 'IGHV1-3', 'OSBP2', 'GYPE', 'C4A', 'MXI1', 'KRT1', 'GFUS', 'SOX6', 'ABCG2', 'ZDHHC19', 'ITLN1', 'FECH', 'IFI27', 'SNCA', 'SGIP1', 'ENSG00000287860', 'LIFR', 'RHCE', 'C9orf153', 'WFDC1', 'CTNNAL1', 'ABCC13', 'OR10Z1', 'HEPACAM2', 'ENSG00000289933', 'ENSG00000224091', 'PAGE2B', 'ENSG00000264066', 'SELENBP1', 'ALAS2', 'ART4', 'RELN', 'HBG2', 'HEMGN', 'THEM5', 'ENSG00000290318', 'ENSG00000285498', 'GYPB', 'ENSG00000255477', 'ENSG00000287510', 'ENSG00000285117', 'BPGM', 'ESRG', 'RHAG', 'TMEM158', 'ABCB5', 'FAM83A', 'SLC6A19']
    list_drop_samples = []
    # list_drop_samples = ["C19-C045-V2",
    #                     "C19-C045-V3",
    #                     "C19-C047-V2",
    #                     "C19-C047-V3",
    #                     "C19-C050-V2",
    #                     "C19-C050-V3",
    #                     "C19-C051-V2",
    #                     "C19-C051-V3",
    #                     "C19-C052-V2",
    #                     "C19-C052-V3",
    #                     "C19-C053-V2",
    #                     "C19-C053-V3",
    #                     "C19-C055-V2",
    #                     "C19-C055-V3",
    #                     "C19-C056-V2",
    #                     'C19-C056-V3',
    #                     'C19-C060-V2',
    #                     'C19-C060-V3',
    #                     'C19-C061-V2',
    #                     'C19-C061-V3']
    colsample = "ID"
    n_components = 20
    outdir = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/8_pca"
    os.makedirs(outdir, exist_ok=True)
    pcx_num = 2
    pcy_num = 1 
    dict_order = {"Visit_order":[{"Visit1":0, "Visit2":1, "Visit3":2, "Visit4":3, "Recover":4}], "Severity":[{"1.0":0, "2.0":1, "3.0":2, "4.0":3, "Convalescent":-1}], "Severity_Visit":[{"Mild_Visit1":0, "Mild_Visit2":1, "Mild_Visit3":2, "Mild_Visit4":3, "Convalescent_Recover":4}, {"Severe_Visit1":0, "Severe_Visit2":1, "Severe_Visit3":2, "Severe_Visit4":3, "Convalescent_Recover":4}]}

# %%    
    main(expmtx, metadata, list_gene, list_drop_samples, colsample, n_components, outdir, pcx_num, pcy_num, dict_order)