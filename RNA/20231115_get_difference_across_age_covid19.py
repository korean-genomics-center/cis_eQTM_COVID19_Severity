# %%
import os
import glob
import numpy as np
import pandas as pd
import pickle
import gzip
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
# %%
def main(path_sev_info, dir_methylcpgmin, list_drop_samples, infilename, outfilename):
    list_files_methylcpgmin = get_list_files_methylcpgmin(dir_methylcpgmin)
    list_methyl_sample = get_list_methyl_sample(list_files_methylcpgmin)
    list_methyl_sample = list(filter(lambda x: "C19-C" in x, list_methyl_sample))
    list_methyl_sample = list(filter(lambda x: x not in list_drop_samples , list_methyl_sample))
    list_methyl_sample = list(filter(lambda x: "L1" not in x, list_methyl_sample))
    # dict_severity_group = get_dict_severity_group(path_sev_info, list_methyl_sample)
    # dict_severity_group = remove_severity_group_no_content(dict_severity_group)
    dict_cpgmin_all_samples = load_pickle(infilename)     
    dict_cpgmin_all_samples_marker_fixed = dict()
    for sampleid, dict_cpgs in dict_cpgmin_all_samples.items():
        if sampleid in list_methyl_sample:
            dict_cpgs_marker_fixed = dict()
            for marker, freqC in dict_cpgs.items():
                fixed_marker = ":".join(marker)
                dict_cpgs_marker_fixed[fixed_marker] = freqC
            
            dict_cpgmin_all_samples_marker_fixed[sampleid] = dict_cpgs_marker_fixed

    df_cpgmin_all_sample = pd.DataFrame(dict_cpgmin_all_samples_marker_fixed).astype(float)
    df_all_sample_cpgmin = df_cpgmin_all_sample.T.reset_index(drop=False).rename(columns={"index": "Sample_ID"})

    # imputation by median if needed (deprecated)
    df_all_sample_cpgmin_imputed_na = df_all_sample_cpgmin.copy()
    list_markers = list(df_all_sample_cpgmin_imputed_na.columns[1:])
    for marker in list_markers:
        median = df_all_sample_cpgmin_imputed_na[marker].median()
        df_all_sample_cpgmin_imputed_na[marker] = df_all_sample_cpgmin_imputed_na[marker].fillna(median)
    df_all_sample_cpgmin = df_all_sample_cpgmin_imputed_na.copy()

    df_sev = pd.read_csv(path_sev_info, sep="\t")
    df_merged = pd.merge(df_all_sample_cpgmin, df_sev, on="Sample_ID", how="inner")
    df_merged = df_merged.set_index("Sample_ID")

    list_markers = list(filter(lambda x: x.startswith("chr"), list(df_merged.columns)))

    direction = "hypo" 
    path_exp = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/4_expmtx/ConfirmedRecovered/expression_matrix_genes.results_TPM.tsv"
    df_TPM = read_TPM_matrix(path_exp, list_drop_samples)
    # df_TPM = select_TPM_matrix(df_TPM, select_pattern="L", delim_id="-", namepos=2)
    path_meta = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/9_clinical/Infectomics_Severity_Information_RNA_20231106.tsv"
    path_methyl_marker = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231020/methyl_{direction}_severe_mild_annot_genelist.tsv"
    outdir = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/{direction}"
    list_filter_genes = list()
    if list_filter_genes == list():
        list_methyl_markers = get_list_methylation_markers(path_methyl_marker)
        df_TPM_gene_filtered = filter_list_genes(df_TPM, list_methyl_markers, mode="genesymbol")
        df_TPM_transposed = transpose_TPM_matrix(df_TPM_gene_filtered)
        df_TPM_split_phase = split_sample_by_phase(df_TPM_transposed)
        df_TPM_meta = merge_TPM_matrix_meta(df_TPM_split_phase, path_meta)
        df_TPM_meta = attach_recovered_severity(df_TPM_meta)
        list_select_gene = df_TPM_gene_filtered["ID"].to_list()
        df_gene_visit_sorted = make_dataframe_stat_test(df_TPM_meta, list_select_gene, outdir)
        df_gene_visit_sorted_sig_only = df_gene_visit_sorted[df_gene_visit_sorted["testsignificant"]==True]
        # df_gene_visit_sorted_sig_only = df_gene_visit_sorted_sig_only.head(25)
        list_filter_genes.extend(df_gene_visit_sorted_sig_only["ID"].to_list())

    list_filter_genesymbols = list(map(lambda x: "_".join(x.split("_")[1:]), list_filter_genes))
    dict_marker_genesymbol_exp = dict()
    for marker in list_markers:
        genesymbol = dict_marker_symbol[marker]
        if genesymbol in list_filter_genesymbols:
            dict_marker_genesymbol_exp[marker] = genesymbol

    list_exp_altered_markers = list(dict_marker_genesymbol_exp.keys())

    df_all_sample_cpgmin_indexed = df_all_sample_cpgmin.set_index("Sample_ID")
    df_all_sample_cpgmin_exp_changed = df_all_sample_cpgmin_indexed[list_exp_altered_markers].reset_index(drop=False)
    df_all_sample_cpgmin = df_all_sample_cpgmin_exp_changed

    from scipy.stats import pearsonr
    from sklearn.linear_model import LinearRegression
    import numpy as np
    import matplotlib.pyplot as plt


    direction = "hypo" 
    annot_methyl = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231020/annotation/methyl_{direction}_severe_mild_annotatr.tsv"
    df_annot_methyl = pd.read_csv(annot_methyl, sep="\t")
    df_annot_methyl["marker"] = df_annot_methyl["seqnames"].astype(str) + ":" + df_annot_methyl["start"].apply(lambda x: int(x)-1).astype(str)

    dict_marker_symbol = dict()
    for marker, symbol in zip(df_annot_methyl["marker"], df_annot_methyl["annot.symbol"]):
        if dict_marker_symbol.get(marker, None) == None:
            dict_marker_symbol[marker] = list()
        dict_marker_symbol[marker].append(symbol)

    for key, list_val in dict_marker_symbol.items():
        list_unique_val = list(set(list_val))
        list_unique_val.remove(np.nan)
        if len(list_unique_val) > 1:
            unique_val = list_unique_val[-1]
        else:
            unique_val = "-".join(list_unique_val)
        if unique_val == "":
            unique_val = "None"
        dict_marker_symbol[key] = unique_val

    df_severe_sorted_by_age = df_merged[df_merged["Severity_visit"] == "Severe_First"]

    list_age_marker = list()
    for marker in list_markers:
        genesymbol = dict_marker_symbol[marker]
        meanbeta = df_severe_sorted_by_age[marker].to_numpy()
        age = df_severe_sorted_by_age["Sample_Age"].to_numpy()
        corr, pval = pearsonr(meanbeta, age)
        if pval < 0.05:
            X = age
            y = meanbeta
            list_age_marker.append(genesymbol)
            m, b = np.polyfit(X, y, 1)
            plt.scatter(X, y, color="k")
            plt.plot(X, m*X+b, color="r")
            plt.title(f"{marker} ({genesymbol})")
            plt.xlabel("Sample_Age")
            plt.ylabel("Mean_Beta")
            plt.show()
            plt.close()

    dict_marker_age = dict()
    for marker in list_markers:
        genesymbol = dict_marker_symbol[marker]
        meanbeta = df_severe_sorted_by_age[marker].to_numpy()
        age = df_severe_sorted_by_age["Sample_Age"].to_numpy()
        corr, pval = pearsonr(meanbeta, age)
        if pval < 0.05:
            dict_marker_age[genesymbol] = {"corr": corr, "pval": pval}
    
    df = pd.DataFrame.from_dict(dict_marker_age, orient="index")
    df["abs_corr"] = df["corr"].apply(abs)
    df_corr_sorted = df.sort_values(by=["abs_corr"], ascending=False)





def stratify_split_train_test_data(X: pd.DataFrame, y: pd.DataFrame, stratify="", train_ratio=0.70, random_state=1) -> pd.DataFrame:
    X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=stratify, test_size= 1 - train_ratio, random_state=random_state)

    return X_train, X_test, y_train, y_test


def _reset_indices(df: pd.DataFrame) -> pd.DataFrame:
    df_reidx = df.reset_index(drop=True)

    return df_reidx

def draw_histogram_train_test(y_train, y_test, **kwargs):
    fig, axes = plt.subplots(1,2, sharey=True)
    axes[0].hist(y_train, label="Train", color="skyblue", **kwargs)
    axes[1].hist(y_test, label="Test", color="orange", **kwargs)
    axes[0].legend(loc="best")
    axes[1].legend(loc="best")
    fig.supylabel("Count", fontsize=12)
    fig.supxlabel("Severity_group", fontsize=12)
    fig.suptitle("Histogram of Stratified Split Between Train and Test by Severity_group")
    plt.show()
    plt.close()

def get_array_array_x(X):
    array_array_x = np.array(X.T.iloc[1:, :].values).astype(np.float64)

    return array_array_x

def get_array_y(y):
    array_array_y = np.array(y)

    return array_array_y


def train_logistic_regression_model(X_train, y_train) -> dict:
    model = LogisticRegression(solver='liblinear').fit(X_train, y_train)

    return model

def test_logistic_regression_model(model: object, X_test, y_test) -> dict:
    slope = model.coef_[0]
    score = model.score(X_test, y_test)
    pred = model.predict(X_test)

    return slope, score, pred


def get_list_files_methylcpgmin(dir_methylcpgmin):
    list_files_methylcpgmin = glob.glob(f"{dir_methylcpgmin}/**/*pair_merged.methyl_cpg_min.tsv", recursive=True)

    return list_files_methylcpgmin

def get_list_methyl_sample(list_files_methylcpgmin):
    list_methyl_sample = list()
    for file_methylcpgmin in list_files_methylcpgmin:
        dir_methylcpgmin = os.path.basename(os.path.dirname(file_methylcpgmin))
        if dir_methylcpgmin == "HealthyControl":
            name_sample = os.path.basename(file_methylcpgmin).split(".")[0]
        else:
            name_sample = os.path.basename(dir_methylcpgmin)
        list_methyl_sample.append(name_sample)

    return list_methyl_sample

def get_dict_sample_severity(path_sev_info, list_methyl_sample):
    dict_sample_severity = dict()
    with open (path_sev_info, mode="r") as fr:
        _skiprow = fr.readline()
        for line in fr:
            record = line.rstrip("\n").split("\t")
            sampleid = record[0]
            severity_visit = record[-1]
            if sampleid in list_methyl_sample:
                dict_sample_severity[sampleid] = severity_visit

    return dict_sample_severity

# def get_dict_severity_group(path_sev_info, list_methyl_sample):
#     dict_severity_group = dict()
#     with open(path_sev_info, mode="r") as fr:
#         _skiprow = fr.readline()
#         for line in fr:
#             record = line.rstrip("\n").split("\t")
#             sampleid = record[0]
#             severity_visit = record[-1]

#             if dict_severity_group.get(severity_visit) == None:
#                 dict_severity_group[severity_visit] = list()
            
#             for methyl_sample in list_methyl_sample:
#                 if sampleid == methyl_sample:
#                     dict_severity_group[severity_visit].append(sampleid)

#     return dict_severity_group

# def remove_severity_group_no_content(dict_severity_group):
#     list_severity_group_no_content = list()
#     for key, list_val in dict_severity_group.items():
#         if list_val == list():
#             list_severity_group_no_content.append(key)
    
#     for sev_grp in list_severity_group_no_content:
#         dict_severity_group.pop(sev_grp)
            
#     return dict_severity_group    

def load_pickle(loadfilename):
    with gzip.open(loadfilename,'rb') as fr:
        data = pickle.load(fr)
    
    return data

def get_dict_mean_beta_all_samples(dict_cpgmin_all_samples):
    dict_meanbeta_all_samples = dict()
    for sampleid, dict_cpgmin in dict_cpgmin_all_samples.items():
        list_beta = list(map(float, dict_cpgmin.values()))
        mean_beta = np.mean(list_beta)
        dict_meanbeta_all_samples[sampleid] = mean_beta

    return dict_meanbeta_all_samples
# %%
if __name__ == "__main__":
    direction = "hypo"
    path_sev_info = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/9_clinical/Infectomics_Severity_Information_Methyl_20231106.tsv"
    dir_methylcpgmin = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/MethylCpGMin"
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
    infilename = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{direction}/dictionary_marker_freqC_all_samples_20231101.pk.gz"
    outfilename = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{direction}/mean_beta_by_severiy_table_20231101.tsv"
# %%
    main(path_sev_info, dir_methylcpgmin, list_drop_samples, infilename, outfilename)