
# %%
import math
import os
import glob

# from scipy.stats import ttest_ind
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from functools import partial
import warnings
warnings.filterwarnings("ignore")


# %% [clinicals]
def main(path_excel, dir_methylcpgmin, list_drop_samples, outfilename, mode):
    df_crf = read_crf_file(path_excel, list_drop_samples)
    dict_rename = {"Subject_NO.(고유번호)":"Sample_ID", "성별":"Sample_Sex", "만나이":"Sample_Age", "중증도분류":"Severity"}
    list_columns_needed = [x for x in df_crf.columns if x in dict_rename.keys()]
    df_crf_columns_needed = df_crf[list_columns_needed]
    df_sev_info = df_crf_columns_needed.rename(columns=dict_rename)
    df_sev_info["Visit_order"] = df_sev_info["Sample_ID"].apply(get_visit_order)
    df_sev_info["Visit"] = df_sev_info["Visit_order"].apply(partial(get_visit_info, mode=mode))
    df_sev_info["Subject_ID"] = df_sev_info["Sample_ID"].apply(get_subject_info)
    df_sev_info["Sample_Sex"] = df_sev_info["Sample_Sex"].apply(get_sex_info)
    df_sev_info["Severity_group"] = df_sev_info["Severity"].apply(get_sev_group)
    df_sev_info["Severity_visit"] = df_sev_info["Severity_group"] + "_" + df_sev_info["Visit"]

    list_clincal_samples = df_sev_info["Sample_ID"].to_list()
    list_clinical_subjects = df_sev_info["Subject_ID"].to_list()
    list_files_methylcpgmin = get_list_files_methylcpgmin(dir_methylcpgmin)
    list_name_sample = get_list_name_sample(list_files_methylcpgmin)
    list_name_subject = list(map(lambda x: "-".join(x.split("-")[:2]), list_name_sample))
    
    list_missing_recovered_subjid = list()
    for methyl_subject, methyl_sample in zip(list_name_subject, list_name_sample):
        if methyl_sample not in list_clincal_samples and methyl_subject in list_clinical_subjects:
            list_missing_recovered_subjid.append(methyl_subject)
    
    df_sev_info_missing_recovered = df_sev_info[df_sev_info["Subject_ID"].isin(list_missing_recovered_subjid)]
    df_sev_info_missing_recovered["Sample_ID"] = df_sev_info_missing_recovered["Sample_ID"].apply(lambda x: x.replace("L", "V"))
    df_sev_info_missing_recovered["Severity_visit"] = df_sev_info_missing_recovered["Severity_visit"].apply(lambda x: x.replace("LongCOVID", "Convalescent"))
    df_sev_info_missing_recovered_added = pd.concat([df_sev_info, df_sev_info_missing_recovered], axis=0)

    list_healthy_control = list()
    for methyl_sample in list_name_sample:
        if methyl_sample not in list_clincal_samples and methyl_sample.startswith("KU10K") and ("V" not in methyl_sample) and ("CT" not in methyl_sample) and ("CR" not in methyl_sample):
            list_healthy_control.append(methyl_sample)
    
    df_healthy_control = pd.DataFrame(columns=df_sev_info.columns)
    df_healthy_control["Sample_ID"] = list_healthy_control
    df_healthy_control["Severity"] = "0"
    df_healthy_control["Severity_visit"] = "Healthy_Control"
    df_sev_info_added_healthy_control = pd.concat([df_sev_info_missing_recovered_added, df_healthy_control], axis=0)

    df_sev_info_final = df_sev_info_added_healthy_control.sort_values(by=["Sample_ID"], ascending=True)
    df_sev_info_final.to_csv(outfilename, sep="\t", index=False)


def get_subject_info(sampleid):
    subjid = "-".join(sampleid.split("-")[:2])
    
    return subjid
    
def get_visit_order(sampleid):
    if sampleid.startswith("C19-C"):
        visit_info = sampleid.split("-")[-1]
        visit_info = visit_info.replace("V", "Visit")
    if sampleid.startswith("C19-R"):
        visit_info = "Visit5"
    if sampleid.endswith("L1"):
        visit_info = "Visit6"
    
    visit_info = str(visit_info)

    return visit_info

def get_visit_info(visit_order, mode=""):
    if visit_order == "Visit1":
        visit_info = "First"
    if visit_order == "Visit2":
        visit_info = "Second"
    if visit_order == "Visit3":
        if mode == "Methyl":
            visit_info = "Last"
        else:
            visit_info = "Third"
    if visit_order == "Visit4":
        if mode == "Methyl":
             visit_info = "Last"
        else:
            visit_info = "Fourth"
    if visit_order == "Visit5":
        visit_info = "Convalescent"
    if visit_order == "Visit6":
        visit_info = "LongCOVID"

    return visit_info

def get_sex_info(sex):
    if sex == "남자":
        sex = 1
    if sex == "여자":
        sex = 2
    
    return sex

def get_sev_group(sev):
    if sev == 1 or sev == 2:
        sev_group = "Mild"
    else:
        sev_group = "Severe"

    return sev_group


def read_excel(*args, **kwargs):
    df_excel = pd.read_excel(*args, **kwargs)

    return df_excel

def correct_header_dataframe(df_crf_raw):
    list_columns = list(df_crf_raw.columns)
    list_new_columns = list()
    for colname in list_columns:
        new_colname = colname.replace("\n", "_").replace(" ", "_").replace("\t", "_").replace("3-1. ", "").replace("3-2. ", "").replace("3-3. ", "").rstrip("_")
        list_new_columns.append(new_colname)
    
    df_crf_raw.columns = list_new_columns

    return df_crf_raw

def filter_dataframe(df_crf, list_drop_samples):
    if len(list_drop_samples) > 0:
        df_crf_filtered = df_crf[~df_crf.iloc[:, 0].isin(list_drop_samples)]
    else:
        df_crf_filtered = df_crf
        
    return df_crf_filtered

def select_dataframe(df_crf, num=2):
    df_crf_select = df_crf.iloc[:, num:]

    return df_crf_select

def read_crf_file(path_excel, list_drop_samples):
    df_crf = pd.read_excel(path_excel, engine="openpyxl", sheet_name="21-22등록 대상자_modify", skiprows=1)
    df_crf = correct_header_dataframe(df_crf)
    df_crf = select_dataframe(df_crf, num=2)
    df_crf = filter_dataframe(df_crf, list_drop_samples)

    return df_crf

def get_list_files_methylcpgmin(dir_methylcpgmin):
    list_files_methylcpgmin = glob.glob(f"{dir_methylcpgmin}/**/*pair_merged.methyl_cpg_min.tsv", recursive=True)

    return list_files_methylcpgmin

def get_list_name_sample(list_files_methylcpgmin):
    list_name_sample = list()
    for file_methylcpgmin in list_files_methylcpgmin:
        dir_methylcpgmin = os.path.basename(os.path.dirname(file_methylcpgmin))
        if dir_methylcpgmin == "HealthyControl":
            name_sample = os.path.basename(file_methylcpgmin).split(".")[0]
        else:
            name_sample = os.path.basename(dir_methylcpgmin)
        list_name_sample.append(name_sample)

    return list_name_sample

# %%
if __name__ == "__main__":
    mode = "Methyl"
    path_excel = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Resources/Data/infectomics_CRF_20230410_edit.xlsx" 
    dir_methylcpgmin = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/MethylCpGMin"
    # list_drop_samples = []
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
    outfilename = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/9_clinical/Infectomics_Severity_Information_{mode}_20240102.tsv"
    main(path_excel, dir_methylcpgmin, list_drop_samples, outfilename, mode)
# %%
#%%
import pandas as pd

table_crf = pd.read_csv("/BiO/Research/Project2/CardiomicsMethylome/Cardi_Stress_Methyl_vs_Clinical.20230517/Results/Clinical_Outlier_Removal/AMI_Clinical_information_20230531_re.fix_column_name.NA_processed.fix_abnormal_values.fix_date.methylation_sequencing.outlier_IQR_filtered.tsv", sep = '\t')
# %%
table_dmp = pd.read_csv("/BiO/Research/Project2/CardiomicsMethylome/Cardiomics_10_fold_cross_validation_DMP.20230302/Results/DMPExtract/Cardiomics.Final_Train_Set.10Fold.OV_MN.20230822/Methylation_DMP_Extract.Control_Control.Case_Case.seed0.kmeans.train_set.OV_MN.Test_Chisq.Cov_Sex_Age.10_fold_overlap.tsv", sep = '\t')
# %%
table_cpg = pd.read_csv("/BiO/Research/Project2/CardiomicsMethylome/Cardiomics_10_fold_cross_validation_DMP.20230302/Results/MethylCpGTable/Cardiomics.10_fold_overlap.Train_set.OV_MN/Cardiomics.MinPerGroup.NULL.CovCutoff.10.control_case.pca_filtered.case_filtered.extra_control_filtered.filter_dmp.Train_set.OV_MN.10_fold_overlap.tsv", sep = '\t')
table_cpg["posname"] = table_cpg.apply(lambda x : f"{x['chr']}_{x['start']}_{x['end']}", axis = 1)
#%%
####################################
table_cpg_filter = pd.read_csv("/BiO/Research/Project2/CardiomicsMethylome/Cardiomics_10_fold_cross_validation_DMP.20230302/Results/DMP_CoxRegression/Death/Filtered/DMPExtract.Cox_Regression.total_samples.sample_filtered.pca_filtered.case_filtered.extra_control_filtered.seed0.kmeans.Train_set_only.Cov_Sex_Age.20231208.limit_365day.filtered.pval_atleast_single_05.meandiff_atleast_single_10.filtered.pval_05.dmp.tsv", sep = '\t')
# table_cpg = table_cpg[table_cpg["posname"].isin(table_cpg_filter["position"].to_list())]
####################################
# %%
list_exclude_checking = ["Number", "KU10K_ID", "Birth-date", "PCI-date", "last-FU_date", "Death_date", "AMI_date", "Any_repeat_revascularization-date", "종류", "기타", "AMI_여부_from_의사"]
list_additional_categorical = ["HTN", "DM", "Current-smok", "FHX", "prev-MI", "Prev-PCI", "Prev-CABG", "ECG", "Disch-aspirin", "Disch-P2Y12inh", "Disch-statin", "Disch-Bblock", "Death", "AMI", "Any_repeat_revascularization", "TypeI_MI", "Stent_PCI_여부"]

def calcuate_ratio(a, b):
    if pd.isna(a) or pd.isna(b):
        return pd.NA
    else:
        return a/b

table_crf["Neutrophil_to_Lymphocyte_Ratio"] = table_crf.apply(lambda x : calcuate_ratio(x["neutrophil"], x["lymphocyte"]), axis = 1)

dict_column_type = dict()
for col in table_crf.columns:
    if col in list_exclude_checking:
        continue
    if table_crf[col].dtypes == object:
        dict_column_type[col] = "Cat"
    elif table_crf[col].dtypes == float:
        if col in list_additional_categorical:
            dict_column_type[col] = "Cat"
        else:
            dict_column_type[col] = "Seq"
    else:
        raise Exception(col)
dict_column_type["Neutrophil_to_Lymphocyte_Ratio"] = "Seq"
# %%
############################################
list_checking_only = ["Age","Sex","BMI","HTN","DM","Current-smok","T-chol","LDL-chol","HDL-chol","TG","WBC","RBC","platelets","neutrophil","lymphocyte","monocyte","eosinophils","Neutrophil_to_Lymphocyte_Ratio","Basophil","CPK_peak","CKMB_peak","TnI_peak","EF","Base_TIMI","Number_diseased_vessel"]

############################################
#%%
import numpy as np
table_clinical_significance = pd.DataFrame(columns = ["chr", "start", "end"] + list_checking_only)

from scipy.stats import spearmanr, mannwhitneyu, kruskal, ranksums

list_sample = table_crf["KU10K_ID"].to_list()
for ind, row in table_cpg.iterrows():
    dict_row = row.to_dict()
    list_sample_check = list(set(list_sample) & set(dict_row.keys()))
    table_clinical_significance.loc[ind, "chr"] = dict_row["chr"]
    table_clinical_significance.loc[ind, "start"] = dict_row["start"]
    table_clinical_significance.loc[ind, "end"] = dict_row["end"]
    
    list_methylation = list(map(dict_row.__getitem__, list_sample_check))
    for col in list_checking_only:
        dtype = dict_column_type[col]
        dict_sample_to_clinical = dict(zip(table_crf["KU10K_ID"], table_crf[col]))
        list_clinical = list(map(dict_sample_to_clinical.__getitem__, list_sample_check))
        
        zip_methyl_clinical = list(zip(list_methylation, list_clinical))
        zip_methyl_clinical_nonna = list(filter(lambda x : pd.notna(x[1]), zip_methyl_clinical))
        
        list_methyl_check = list(map(lambda x : x[0], zip_methyl_clinical_nonna))
        list_clinical_check = list(map(lambda x : x[1], zip_methyl_clinical_nonna))
        
        if dtype == "Cat":
            list_cat_type = list(set(list_clinical_check))
            list_matrix = list()
            for cat_type in list_cat_type:
                zip_methyl_clinical_cat = list(filter(lambda x : x[1] == cat_type, zip_methyl_clinical_nonna))
                list_methyl_this_cat = list(map(lambda x : x[0], zip_methyl_clinical_cat))
                list_matrix.append(list_methyl_this_cat)
            if len(list_cat_type) == 2:
                _, pval = ranksums(list_matrix[0], list_matrix[1])
                stat = abs(np.mean(list_matrix[0]) - np.mean(list_matrix[1]))
            elif len(list_cat_type) == 1:
                pval = 1
            else:
                _, pval = kruskal(*list_matrix)
        else:
            stat, pval = spearmanr(list_methyl_check, list_clinical_check)
        table_clinical_significance.loc[ind, col] = pval
        table_clinical_significance.loc[ind, f"{col}_Stat"] = stat
# %%
from statsmodels.stats.multitest import fdrcorrection

for col in list_checking_only:
    _, fdr = fdrcorrection(table_clinical_significance[col])
    table_clinical_significance[f"{col}_FDR"] = fdr
# %%
list_column_to_priority = list()
for col in list_checking_only:
    num_significant = sum(table_clinical_significance[col] < 0.05)
    priority_adjust = 0 if dict_column_type[col] == "Seq" else -500
    priority = num_significant + priority_adjust
    list_column_to_priority.append((col, priority))

list_column_ordered_to_plot = list(map(lambda x : x[0], sorted(list_column_to_priority, key = lambda x : x[1], reverse = True)))
# %%
table_clinical_significance["posname"] = table_clinical_significance.apply(lambda x : f"{x['chr']}:{x['start']}", axis = 1)
table_clinical_significance["posname_"] = table_clinical_significance.apply(lambda x : f"{x['chr']}_{x['start']}_{x['end']}", axis = 1)
# %%
table_dmp_group = pd.read_csv("/BiO/Research/Project2/CardiomicsMethylome/Cardiomics_10_fold_cross_validation_DMP.20230302/Results/DMPExtract/Cardiomics.Final_Train_Set.10Fold.OV_MN.20230822/Methylation_DMP_Extract.Control_Control.Case_Case.seed0.kmeans.train_set.OV_MN.Test_Chisq.Cov_Sex_Age.10_fold_overlap.mean_methdiff.group_marker_only.group_interval_100bp.tsv", sep = '\t')
table_dmp_direction = pd.read_csv("/BiO/Research/Project2/CardiomicsMethylome/Cardiomics_10_fold_cross_validation_DMP.20230302/Results/DMPExtract/Cardiomics.Final_Train_Set.10Fold.OV_MN.20230822/Methylation_DMP_Extract.Control_Control.Case_Case.seed0.kmeans.train_set.OV_MN.Test_Chisq.Cov_Sex_Age.10_fold_overlap.mean_methdiff.tsv", sep = '\t')
# %%
table_dmp_group["posname"] = table_dmp_group.apply(lambda x : f"{x['chr']}_{x['start']}_{x['end']}", axis = 1)
table_dmp_direction["posname"] = table_dmp_direction.apply(lambda x : f"{x['chr']}_{x['start']}_{x['end']}", axis = 1)
dict_posname_to_direction = dict(zip(table_dmp_direction["posname"], list(map(lambda x : int(x > 0),  table_dmp_direction["mean_methdiff"]))))
table_cpg_posonly = table_cpg[["chr", "start", "end", "posname"]].copy()
table_cpg_posonly["Non_Neighbor"] = table_cpg_posonly["posname"].apply(lambda x : int(x not in table_dmp_group["posname"].to_list()))
table_cpg_posonly["Is_Hyper"] = table_cpg_posonly["posname"].apply(dict_posname_to_direction.__getitem__)
dict_sexchr_to_int = {
    "X":"23",
    "Y":"24",
    "M":"25"
}
table_cpg_posonly["chr_int"] = table_cpg_posonly["chr"].apply(lambda x : int(dict_sexchr_to_int.get(x.replace("chr", ''), x.replace("chr", ''))))
table_cpg_posonly = table_cpg_posonly.sort_values(by = ["Is_Hyper", "chr_int", "start"])
# %%
table_cpg_posonly = table_cpg_posonly.reset_index(drop = True)

dict_posname_to_order = dict(zip(table_cpg_posonly["posname"], table_cpg_posonly.index))
# %%
table_clinical_significance = table_clinical_significance.set_index("posname_", drop = False)
# %%
list_column_ordered_to_plot_stat = list(map(lambda x : f"{x}_Stat", list_column_ordered_to_plot))
list_column_ordered_to_plot_fdr = list(map(lambda x : f"{x}_FDR", list_column_ordered_to_plot))
table_clinical_significance_plot = table_clinical_significance.loc[table_cpg_posonly["posname"].to_list(), list_column_ordered_to_plot_stat]
table_clinical_significance_plot = table_clinical_significance_plot.rename(columns = dict(zip(list_column_ordered_to_plot_stat, list_column_ordered_to_plot)))
# %%
table_clinical_significance_pval = table_clinical_significance.loc[table_cpg_posonly["posname"].to_list(), list_column_ordered_to_plot_fdr]
table_clinical_significance_pval = table_clinical_significance_pval.rename(columns = dict(zip(list_column_ordered_to_plot_fdr, list_column_ordered_to_plot)))
table_clinical_significance_pval = table_clinical_significance_pval.applymap(lambda x : "*" if x < 0.05 else "")
# %%
table_clinical_significance_plot_T = table_clinical_significance_plot.T
table_clinical_significance_pval_T = table_clinical_significance_pval.T
# %%
table_clinical_significance_plot_T = table_clinical_significance_plot_T.rename(columns = dict(zip(table_clinical_significance["posname_"], table_clinical_significance["posname"])))
list_sequential_traits = list(filter(lambda x : dict_column_type[x] == "Seq", table_clinical_significance_plot_T.index))
list_categorical_traits = list(filter(lambda x : dict_column_type[x] == "Cat", table_clinical_significance_plot_T.index))
#%%
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import gridspec
import seaborn as sns
from matplotlib.cm import ScalarMappable

plt.rcParams['font.family'] = 'DeJavu Serif'
plt.rcParams['font.serif'] = ['Arial']
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams['font.size'] = 24

fig = plt.figure(figsize=(35, 12))
row = 100
col = 100
fontsize_label = 40
gsfig = gridspec.GridSpec(
    row, col, 
    left=0, right=1, bottom=0,
    top=1, wspace=1, hspace=1)

gs1 = gsfig[0:80, 0:95] 
ax1 = fig.add_subplot(gs1)

gs2 = gsfig[80:95, 0:95]
ax2 = fig.add_subplot(gs2)

gs3 = gsfig[96:100, 0:95] 
ax3 = fig.add_subplot(gs3)

gs4 = gsfig[0:40, 95:97]
ax4 = fig.add_subplot(gs4)

gs5 = gsfig[80:95, 95:97]
ax5 = fig.add_subplot(gs5)

sns.heatmap(data = table_clinical_significance_plot_T.loc[list_sequential_traits, :], annot = table_clinical_significance_pval_T.loc[list_sequential_traits, :].to_numpy(), cmap = "coolwarm", vmin = -0.5, vmax = 0.5, fmt = '', annot_kws = {"color":"k", "fontsize":12}, xticklabels = False, ax = ax1, cbar = True, cbar_ax = ax4, cbar_kws={"label" : "Spearman rho", "shrink":0.4})
ax1.set_xlabel("")

sns.heatmap(data = table_clinical_significance_plot_T.loc[list_categorical_traits, :], annot = table_clinical_significance_pval_T.loc[list_categorical_traits, :].to_numpy(), cmap = "Reds", vmin = 0, vmax = 10, fmt = '', annot_kws = {"color":"k", "fontsize":12}, xticklabels = False, ax = ax2, cbar = True, cbar_ax = ax5, cbar_kws={"label" : "|Mean Methylation\nDifference| (%)"})
ax2.set_xlabel("")

# table_cpg_isneighbor = pd.DataFrame(columns = table_clinical_significance_plot_T.columns)
# table_cpg_posonly["posname_"] = table_cpg_posonly.apply(lambda x : f"{x['chr']}:{x['start']}", axis = 1)
# dict_posname_to_neighbor = dict(zip(table_cpg_posonly["posname_"], table_cpg_posonly["Non_Neighbor"]))
# table_cpg_isneighbor.loc["Neighbor_Methylation", :] = list(map(lambda x : dict_posname_to_neighbor[x], list(table_cpg_isneighbor.columns)))
# norm = mpl.colors.Normalize(0,1)
# colors = ['k', 'w']
# blackwhite_cmap = mpl.colors.LinearSegmentedColormap.from_list("", colors)
# table_cpg_isneighbor = table_cpg_isneighbor.astype(int)

# g = sns.heatmap(data = table_cpg_isneighbor, cmap = blackwhite_cmap, xticklabels = False, ax = ax3, cbar = False)
# g.set_yticklabels(g.get_yticklabels(), rotation = 0)
# ax3.set_xlabel("")

table_cpg_ishyper = pd.DataFrame(columns = table_clinical_significance_plot_T.columns)
table_cpg_posonly["posname_"] = table_cpg_posonly.apply(lambda x : f"{x['chr']}:{x['start']}", axis = 1)
dict_posname_to_hyper = dict(zip(table_cpg_posonly["posname_"], table_cpg_posonly["Is_Hyper"]))
table_cpg_ishyper.loc["Direction", :] = list(map(lambda x : dict_posname_to_hyper[x], list(table_cpg_ishyper.columns)))
table_cpg_ishyper = table_cpg_ishyper.astype(int)

hypo_annot_position = sum(table_cpg_posonly["Is_Hyper"] == 0) / 2
hyper_annot_position = sum(table_cpg_posonly["Is_Hyper"] == 0) + sum(table_cpg_posonly["Is_Hyper"] == 1) / 2

g = sns.heatmap(data = table_cpg_ishyper, cmap = "coolwarm", xticklabels = False, ax = ax3, cbar = False)
g.set_yticklabels(g.get_yticklabels(), rotation = 0)
ax3.annotate("Hypo", xy = (hypo_annot_position, 0.5), xytext = (hypo_annot_position, 0.5), ha = "center", va = "center", color = "white", fontweight = "bold")
ax3.annotate("Hyper", xy = (hyper_annot_position, 0.5), xytext = (hyper_annot_position, 0.5), ha = "center", va = "center", color = "white", fontweight = "bold")
ax3.set_xlabel("")

gsfig.tight_layout(fig)
# %%
# ★ 𝜌
