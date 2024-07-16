# %%
import os
import glob
import numpy as np
import pandas as pd
import pickle
import gzip
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection
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
    
#################
    path_exp = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/4_expmtx/ConfirmedRecovered/expression_matrix_genes.results_TPM.tsv"
    df_TPM = read_TPM_matrix(path_exp, list_drop_samples)
    df_TPM = select_TPM_matrix(df_TPM, select_pattern="L", delim_id="-", namepos=2)
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
        list_filter_genes.extend(df_gene_visit_sorted_sig_only["ID"].to_list())
    
    df_TPM_gene_filtered = filter_list_genes(df_TPM, list_filter_genes, mode="geneonly")
    df_TPM_transposed = transpose_TPM_matrix(df_TPM_gene_filtered)
    df_TPM_split_phase = split_sample_by_phase(df_TPM_transposed)
    df_TPM_meta = merge_TPM_matrix_meta(df_TPM_split_phase, path_meta)
    df_TPM_meta = attach_recovered_severity(df_TPM_meta)

    list_markers = list(filter(lambda x: x.startswith("ENSG"), list(df_TPM_meta.columns)))

    df_TPM_meta_selected = df_TPM_meta[df_TPM_meta["ID"].str.contains("C19-C") & df_TPM_meta["ID"].str.contains("V1")]

    X_ = df_TPM_meta_selected[list_markers]
    y = df_TPM_meta_selected["Severity_group"]
    X_train, X_test, y_train, y_test = stratify_split_train_test_data(X_, y, stratify=y, train_ratio=0.70, random_state=1)
    draw_histogram_train_test(y_train, y_test, bins=2, alpha=0.5)

    model = train_logistic_regression_model(X_train.values, y_train.values)
    slope, score, y_pred = test_logistic_regression_model(model, X_test.values, y_test.values)

    # from sklearn.ensemble import RandomForestClassifier
    # clf = RandomForestClassifier(n_estimators=50, criterion='entropy',max_depth=5, max_features='sqrt', bootstrap=True,  oob_score=True, random_state=100)

    # model = clf.fit(X_train, y_train)
    # y_pred = model.predict(X_test)

    print('Training set score: {:.4f}'.format(model.score(X_train, y_train)))
    print('Test set score: {:.4f}'.format(model.score(X_test, y_test)))

    from sklearn.metrics import confusion_matrix

    cm = confusion_matrix(y_test, y_pred)

    print('Confusion matrix\n\n', cm, end = "\n\n")

    TP = cm[0,0]
    TN = cm[1,1]
    FP = cm[0,1]
    FN = cm[1,0]

    print('\nTrue Positives(TP) = ', TP)

    print('\nTrue Negatives(TN) = ', TN)

    print('\nFalse Positives(FP) = ', FP)

    print('\nFalse Negatives(FN) = ', FN)

    false_positive_rate = FP / float(FP + TN)
    print('False Positive Rate : {0:0.4f}'.format(false_positive_rate))

    specificity = TN / (TN + FP)
    print('Specificity : {0:0.4f}'.format(specificity))

    cm_matrix = pd.DataFrame(data=cm, columns=['Actual Positive:1', 'Actual Negative:0'], 
                                 index=['Predict Positive:1', 'Predict Negative:0'])

    sns.heatmap(cm_matrix, annot=True, fmt='d', cmap='YlGnBu')

    from sklearn.metrics import classification_report
    
    print(classification_report(y_test, y_pred))


    from sklearn.metrics import roc_curve

    y_pred1 = model.predict_proba(X_test)[:,1]

    fpr, tpr, thresholds = roc_curve(y_test, y_pred1, pos_label = 'Severe')

    plt.figure(figsize=(6,4))

    plt.plot(fpr, tpr, linewidth=2)

    plt.plot([0,1], [0,1], 'k--' )

    plt.rcParams['font.size'] = 12

    plt.title('ROC curve')

    plt.xlabel('False Positive Rate (1 - Specificity)')

    plt.ylabel('True Positive Rate (Sensitivity)')

    plt.show()


    from sklearn.metrics import roc_auc_score

    ROC_AUC = roc_auc_score(y_test, y_pred1)

    print('ROC AUC : {:.4f}'.format(ROC_AUC))


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


def read_TPM_matrix(path_exp, list_drop_samples):
    df_TPM = pd.read_csv(path_exp, sep="\t")
    df_TPM  = df_TPM.drop(columns=list_drop_samples)

    return df_TPM

def select_TPM_matrix(df_TPM, select_pattern, delim_id="-", namepos=0):
    list_colname = list(df_TPM.columns)
    list_colname_filt = list(filter(lambda x: x.split(delim_id)[namepos][0] != select_pattern if len(x.split("-")) > 1 else x, list_colname))
    df_TPM = df_TPM.loc[:, list_colname_filt]

    return df_TPM

def get_list_methylation_markers(path_methyl_marker):
    list_methyl_markers = list()
    with open(path_methyl_marker, mode="r") as fr:
        for line in fr:
            record = line.rstrip("\n").split("\n")
            list_methyl_markers.extend(record)
    
    return list_methyl_markers

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

def transpose_TPM_matrix(df_TPM, colsample="ID"):
    list_geneid = df_TPM[colsample].to_list()
    df_TPM_transposed = df_TPM.T.iloc[1:, :]
    df_TPM_transposed.columns = list_geneid
    df_TPM_transposed = df_TPM_transposed.reset_index(drop=False).rename(columns={"index": colsample})
    
    return df_TPM_transposed

def get_visit_order(sampleid, dict_manual={"R": "V5", "L": "V6"}):
    visit = sampleid.split("-")[-1]
    for key, val in dict_manual.items():
        if sampleid.split("-")[1][0] == key:
            visit = val
        if sampleid.split("-")[-1][0] == key:
            visit = val
    visit_order = int(visit[1:])

    return visit_order

def split_sample_by_phase(df_TPM_gene_filtered, colsample="ID", colTPM="TPM", coltime="Time"):
    df_split_phase = df_TPM_gene_filtered.copy()
    df_split_phase[coltime] = df_split_phase[colsample].apply(get_visit_order)
    df_split_phase = df_split_phase.sort_values(by=coltime, ascending=True)
    df_split_phase[coltime] = df_split_phase[coltime].astype(str)

    return df_split_phase
    

def merge_TPM_matrix_meta(df_TPM, path_meta, colsample="ID"):
    df_meta = pd.read_csv(path_meta, sep="\t")
    df_meta = df_meta.rename(columns={"Sample_ID": colsample})
    df_merged = pd.merge(df_TPM, df_meta, how="inner", on=colsample)
    
    return df_merged

def attach_recovered_severity(df_TPM_meta):
    df_TPM_meta["Severity"] = df_TPM_meta["Severity"].apply(lambda x: 0.0 if x =="-" else x)
    df_TPM_meta[df_TPM_meta["ID"].str.contains("C19-R")]["Severity_group"] = "Convalescent"
    df_TPM_meta[df_TPM_meta["ID"].str.contains("L1")]["Severity_group"] = "LongCOVID"

    return df_TPM_meta

def mediandiff(a, b):
    deltamedian = np.median(b) - np.median(a)

    return deltamedian

def make_dataframe_stat_test(df_TPM_meta, list_select_gene, outdir):
    dict_gene_visit_pval = dict()
    for gene in list_select_gene:
        df_group_visit = df_TPM_meta.groupby(["Visit", "Severity_group"])[gene].apply(np.array).reset_index(drop=False)

        dict_visit_pval = dict()
        for ind1 in list(df_group_visit.index):
            if ind1 < len(df_group_visit)-1:
                ind2 = int(ind1) + 1
            comp1 = df_group_visit.loc[ind1, :] 
            comp2 = df_group_visit.loc[ind2, :]
            visit_num1 = comp1["Visit"]
            visit_num2 = comp2["Visit"]
            if visit_num1 == visit_num2:
                list_case = comp1[gene]
                list_case = list(map(float, list_case))
                list_control = comp2[gene]
                list_control = list(map(float, list_control))

                if len(set(list_control)) > 1:
                    deltamedian = mediandiff(list_case, list_control)
                    stat, pval = mannwhitneyu(list_case, list_control)
                    stat_res = {"delta": deltamedian, "stat": stat, "pval": pval}
                    dict_visit_pval[visit_num1] = stat_res

        if dict_gene_visit_pval.get(gene) == None:
            dict_gene_visit_pval[gene] = dict()
        
        dict_gene_visit_pval[gene].update(dict_visit_pval)

    df_gene_visit = pd.concat({k: pd.DataFrame(v).T for k, v in dict_gene_visit_pval.items()}, axis=0)
    df_gene_visit = df_gene_visit.reset_index(drop=False).rename(columns={"level_0": "ID", "level_1": "Visit"})

    df_gene_visit_sorted = df_gene_visit.sort_values(by=["pval"], ascending=True)
    list_pval = df_gene_visit_sorted["pval"].to_list()
    list_sig = fdrcorrection(list_pval)[0]
    list_padj = fdrcorrection(list_pval)[1]
    df_gene_visit_sorted["padj"] = list_padj
    df_gene_visit_sorted["testsignificant"] = list_sig
    df_gene_visit_sorted.to_csv(os.path.join(outdir, "stattest_expressed_methylation_markers.tsv"), sep="\t", index=False)

    return df_gene_visit_sorted
# %%
if __name__ == "__main__":
    direction = "hyper"
    visit = "first"
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
    infilename = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{visit}/{direction}/dictionary_marker_freqC_all_samples_20231106.pk.gz"
    outfilename = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/Epigenetic_changes/{visit}/{direction}/mean_beta_by_severiy_table_20231106.tsv"
# %%
    main(path_sev_info, dir_methylcpgmin, list_drop_samples, infilename, outfilename)

# # %%
# annot_methyl = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231020/annotation/methyl_hyper_severe_mild_annotatr.tsv"
# df_annot_methyl = pd.read_csv(annot_methyl, sep="\t")
# df_annot_methyl["marker"] = df_annot_methyl["seqnames"].astype(str) + ":" + df_annot_methyl["start"].apply(lambda x: int(x)-1).astype(str)

# dict_marker_symbol = dict()
# for marker, symbol in zip(df_annot_methyl["marker"], df_annot_methyl["annot.symbol"]):
#     if dict_marker_symbol.get(marker, None) == None:
#         dict_marker_symbol[marker] = list()
#     dict_marker_symbol[marker].append(symbol)

# for key, list_val in dict_marker_symbol.items():
#     list_unique_val = list(set(list_val))
#     list_unique_val.remove(np.nan)
#     if len(list_unique_val) > 1:
#         unique_val = list_unique_val[-1]
#     else:
#         unique_val = "-".join(list_unique_val)
#     if unique_val == "":
#         unique_val = "None"
#     dict_marker_symbol[key] = unique_val

# ##################
#     # imputation by median if needed (deprecated)
#     df_all_sample_cpgmin_imputed_na = df_all_sample_cpgmin.copy()
#     list_markers = list(df_all_sample_cpgmin_imputed_na.columns[1:])
#     for marker in list_markers:
#         median = df_all_sample_cpgmin_imputed_na[marker].median()
#         df_all_sample_cpgmin_imputed_na[marker] = df_all_sample_cpgmin_imputed_na[marker].fillna(median)
# ##################
    # df_sev = pd.read_csv(path_sev_info, sep="\t")
    # df_merged = pd.merge(df_all_sample_cpgmin, df_sev, on="Sample_ID", how="inner")
    # df_merged = df_merged.set_index("Sample_ID")
# ##################
    # df_merged = df_merged[df_merged.index.str.contains("-V4")]
    # list_markers = list(filter(lambda x: x.startswith("chr"), list(df_merged.columns)))
#     list_filter_genesymbols = list(map(lambda x: "_".join(x.split("_")[1:]), list_filter_genes))
#     dict_marker_genesymbol_exp = dict()
#     for marker in list_markers:
#         genesymbol = dict_marker_symbol[marker]
#         if genesymbol in list_filter_genesymbols:
#             dict_marker_genesymbol_exp[marker] = genesymbol
#     list_exp_altered_markers = list(dict_marker_genesymbol_exp.keys())
#     df_all_sample_cpgmin_indexed = df_all_sample_cpgmin.set_index("Sample_ID")
#     df_all_sample_cpgmin_exp_changed = df_all_sample_cpgmin_indexed[list_exp_altered_markers].reset_index(drop=False)
#     df_all_sample_cpgmin = df_all_sample_cpgmin_exp_changed
# ######################
