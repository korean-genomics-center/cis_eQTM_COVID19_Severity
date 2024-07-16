
# %%
import math
import os
import warnings
from collections import Counter

# from scipy.stats import ttest_ind
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scikit_posthocs import posthoc_dunn, posthoc_tukey
from scipy.stats import f_oneway, kruskal, mannwhitneyu, sem, t
from statannot import add_stat_annotation
from statsmodels.stats.multitest import fdrcorrection
from collections import defaultdict
warnings.filterwarnings("ignore")

# %% [clinicals]
def main(path_excel, target, list_drop_samples, fontpath, dict_palette, dict_idx_rename, outdir):
    # set_korean_font(fontpath)
    
    df_crf_raw = pd.read_excel(path_excel, engine="openpyxl", sheet_name="21-22등록 대상자_modify", skiprows=1)
    df_crf_raw = correct_header_dataframe(df_crf_raw)
    df_crf_raw_selec = select_dataframe(df_crf_raw, num=2)
    df_crf_raw_filt = filter_dataframe(df_crf_raw_selec, list_drop_samples)
    df_crf_raw_add_sev = define_severity(df_crf_raw_filt)
    dict_sample_phase = get_dict_sample_visit(df_crf_raw_add_sev)
    list_sample_pattern = get_sample_pattern_need_visit(dict_sample_phase)
    df_crf = extract_dataframe(df_crf_raw_add_sev, list_sample_pattern)
    list_df_compare = get_clinical_assoc_severity(df_crf)

    dict_value_pval_two_groups = defaultdict(dict)
    dict_value_pval_multi_groups = defaultdict(dict)
    for df_compare in list_df_compare:
        value = list(df_compare.columns)[-1]
        series_sev_val = df_compare.groupby(target)[value].apply(np.array)
        list_groups = list(series_sev_val.index)
        list_array_val = list(map(lambda x: series_sev_val[x], list_groups))
        list_cnt = list(map(lambda x: len(series_sev_val[x]), list_groups))
        if (len(list_cnt) < 2):
            continue

        if (len(set(list_array_val[0])) < 3 and max(list_array_val[0]) < 3):
            continue

        if (len(set(list_array_val[0])) < 3):
            continue

        if (len(list_cnt) == 2 and min(list_cnt) > 2):
            df_compare_copy = df_compare.copy()
            df_sev_val = df_compare_copy.groupby(target)[value].apply(np.array)
            control = df_sev_val.index[0]
            list_control = df_sev_val[control]
            case = df_sev_val.index[1]
            list_case = df_sev_val[case]
            try:
                u_stat, pval = mannwhitneyu(list_case, list_control)
                dict_stat_res = {"compar": list(df_sev_val.index), "stat": u_stat, "pval": pval}
                dict_value_pval_two_groups[value].update(dict_stat_res)

                colname = list(df_compare.columns)[-1]
                if colname not in dict_idx_rename.keys():
                    new_colname = colname
                else:
                    new_colname = dict_idx_rename[colname]

                plt.figure(figsize=(3,5))
                flierprops = dict(marker='o', markerfacecolor='None', markersize=5, markeredgecolor='black')
                order = sorted(list(df_compare[target].unique()), reverse=True)
                palette = list(map(lambda x: dict_palette[x], order))

                ax = sns.boxplot(data=df_compare_copy, x=target, y=value, palette=palette, order=order, flierprops=flierprops)
                ax.set_xlabel(target, fontsize=12)
                ax.set_ylabel(new_colname, fontsize=12)
                order = dict_stat_res["compar"]
                pairs = [(case, control)]
                pvalues = [pval]
                add_stat_annotation(ax, data=df_compare, x=target, y=value, order=order,
                    box_pairs=pairs,
                    perform_stat_test=False, pvalues=pvalues,
                    test=None, text_format=f'star', loc='inside', verbose=2)
            
                list_target = df_compare[target].to_list()
                count_target = Counter(list_target)
                new_xticklabels = list()
                new_legend = list()
                for ind_xtick, count_pair in enumerate(count_target.items(), start = 1):
                    key, cnt = count_pair
                    shortkey = key[0].upper()
                    longkey = key[0].upper() + key[1:].lower()
                    xticklabel = f"{shortkey}\n(N={cnt})"
                    new_xticklabels.append(xticklabel)
                plt.xticks(ticks=list(range(0, len(new_xticklabels))), labels=new_xticklabels, fontsize=10)
                plt.tight_layout()
                plt.savefig(f"{outdir}/boxplot_{target}_vs_{new_colname}.png", dpi=300)
                plt.show()
                plt.close()

            except Exception as e:
                print(e)
                    
        if len(list_cnt) > 2 and min(list_cnt) > 2:
            df_compare_copy = df_compare.copy()

            value = list(df_compare.columns)[-1]
            df_dunn = posthoc_dunn(df_compare, val_col=value, group_col=target, p_adjust="fdr_bh")
            remove = np.tril(np.ones(df_dunn.shape), k=0).astype("bool")
            df_dunn[remove] = np.nan
            result_dunn = df_dunn.melt(ignore_index=False).reset_index().dropna()
            result_dunn.columns = ["index", "variable", "Post_hoc_P"]
            result_dunn["category"] = value

            df_sev_val = df_compare_copy.groupby(target)[value].apply(np.array)
            one = df_sev_val.index[0]
            list_one = df_sev_val[one]
            two = df_sev_val.index[1]
            list_two = df_sev_val[two]
            three = df_sev_val.index[2]
            list_three = df_sev_val[three]

            try:
                h_stat, pval = kruskal(list_one, list_two, list_three)
                dict_stat_res = {"compar": list(df_sev_val.index), "stat": h_stat, "pval": pval}
                dict_value_pval_multi_groups[value].update(dict_stat_res)

                colname = list(df_compare.columns)[-1]
                if colname not in dict_idx_rename.keys():
                    new_colname = colname
                else:
                    new_colname = dict_idx_rename[colname]

                plt.figure(figsize=(3,5))
                flierprops = dict(marker='o', markerfacecolor='None', markersize=5, markeredgecolor='black')
                order = sorted(list(df_compare[target].unique()), reverse=True)
                palette = list(map(lambda x: dict_palette[x], order))
                
                ax = sns.boxplot(data=df_compare_copy, x=target, y=value, palette=palette, order=order, flierprops=flierprops)
                ax.set_xlabel(target, fontsize=12)
                ax.set_ylabel(new_colname, fontsize=12)
                order = dict_stat_res["compar"]
                pairs = [(i[1]["index"], i[1]["variable"]) for i in result_dunn.iterrows()]   
                pvalues = [i[1]["Post_hoc_P"] for i in result_dunn.iterrows()]
                add_stat_annotation(ax, data=df_compare, x=target, y=value, order=order,
                    box_pairs=pairs,
                    perform_stat_test=False, pvalues=pvalues,
                    test=None, text_format=f'star', loc='inside', verbose=2)
                
                list_target = df_compare[target].to_list()
                count_target = Counter(list_target)
                new_xticklabels = list()
                new_legend = list()
                for ind_xtick, count_pair in enumerate(count_target.items(), start = 1):
                    key, cnt = count_pair
                    shortkey = key[0].upper()
                    longkey = key[0].upper() + key[1:].lower()
                    xticklabel = f"{shortkey}\n(N={cnt})"
                    new_xticklabels.append(xticklabel)
                plt.xticks(ticks=list(range(0, len(new_xticklabels))), labels=new_xticklabels, fontsize=10)
                plt.tight_layout()
                plt.savefig(f"{outdir}/boxplot_{target}_vs_{new_colname}.png", dpi=300)
                plt.show()
                plt.close()

            except Exception as e:
                print(e)      

    df_value_pval_two_groups = pd.DataFrame.from_dict(dict_value_pval_two_groups, orient="index")
    list_pval_two_groups = df_value_pval_two_groups['pval'].to_list()
    list_sig_two_groups, list_qval_two_groups = fdrcorrection(list_pval_two_groups)
    df_value_pval_two_groups["qval"] = list_qval_two_groups
    df_value_pval_two_groups["sig"] = list_sig_two_groups

    df_value_pval_multi_groups = pd.DataFrame.from_dict(dict_value_pval_multi_groups, orient="index")
    list_pval_multi_groups = df_value_pval_multi_groups['pval'].to_list()
    list_sig_multi_groups, list_qval_multi_groups = fdrcorrection(list_pval_multi_groups)
    df_value_pval_multi_groups["qval"] = list_qval_multi_groups
    df_value_pval_multi_groups["sig"] = list_sig_multi_groups

    df_value_pval_merged = pd.concat([df_value_pval_two_groups, df_value_pval_multi_groups])
    df_value_pval_merged_sorted = df_value_pval_merged.sort_values(by="qval", ascending=True)
    df_value_pval_merged_sorted_idx_rename = df_value_pval_merged_sorted.rename(index=dict_idx_rename)
    df_value_pval_merged_sorted_reidx = df_value_pval_merged_sorted_idx_rename.reset_index(drop=False).rename(columns={"index": "category"})

    df_value_pval_merged_sorted_reidx.to_csv("/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Resources/Data/ClinicalInformation_stat_test_20231121.tsv", sep="\t", index=False)


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

def change_severity_class_to_group(sev_class):
    sev_class = int(sev_class)
    if sev_class == 1 or sev_class == 2:
        sev_group = "mild"

    else:
        sev_group = "severe"
    
    return sev_group

def define_severity(df_crf):
    df_crf[target] = df_crf["중증도분류"].apply(change_severity_class_to_group)
    for sampleid in df_crf.iloc[:, 0]:
        if str(sampleid).split("-")[1][0] == "R":
            df_crf_indexed = df_crf.set_index(df_crf.columns[0])
            df_crf_indexed.loc[sampleid, target] = "Convalescent"
            df_crf = df_crf_indexed.reset_index(drop=False)

    return df_crf

def get_dict_sample_visit(df_crf):
    dict_sample_phase = dict()
    list_sampleid = df_crf.iloc[:, 0].to_list()
    for sampleid in list_sampleid:
        subjectid = "-".join(sampleid.split("-")[:2])
        phase = sampleid.split("-")[-1]
        if dict_sample_phase.get(subjectid) == None:
            dict_sample_phase[subjectid] = list()
        dict_sample_phase[subjectid].append(phase)

    return dict_sample_phase

def get_sample_pattern_need_visit(dict_sample_phase):
    list_sample_pattern = list()
    for subj, list_visit in dict_sample_phase.items():
        sample_pattern = subj + "-" + list_visit[0]
        list_sample_pattern.append(sample_pattern)

    return list_sample_pattern

def extract_dataframe(df_crf, list_sample_pattern):
    df_crf_filt = df_crf[df_crf.iloc[:, 0].isin(list_sample_pattern)]
    
    return df_crf_filt

def get_clinical_assoc_severity(df_crf):

    list_df_compare = list()
    list_columns = list(df_crf.columns)
    df_crf_selec = df_crf.copy()
    for colname in list_columns:
        list_crf_values = df_crf_selec[colname].to_list()
        list_crf_values_nona = list(filter(lambda x : not isinstance(x, float) or not math.isnan(x), list_crf_values))
        list_crf_values_nona = list(map(lambda x: x.replace("PY", "") if "PY" in str(x) else x, list_crf_values_nona))
        if list_crf_values_nona == list():
            continue

        is_numeric = isinstance(list_crf_values_nona[0], int) or isinstance(list_crf_values_nona[0], float)
        if is_numeric:
            try:
                df_crf_selec[colname] = df_crf_selec[colname].apply(float)
                df_compare = df_crf_selec[[target, colname]].dropna()
                list_df_compare.append(df_compare)
            except Exception as e:
                print(e)
    
    return list_df_compare

def set_korean_font(fontpath):
    fontname = os.path.basename(fontpath).replace(".ttf", "")
    fonttype = fm.FontEntry(fname=fontpath, name=fontname)
    fm.fontManager.ttflist.insert(0, fonttype)
    plt.rcParams.update({"font.size": 10, "font.family":fontname})

def get_sample_size(df_compare, target):
    value = list(df_compare.columns)[-1]
    sample_size = df_compare.groupby(target)[value].apply(len)

    return sample_size

def calc_ci(x):
    ci = t.interval(0.95, len(x)-1, loc=np.mean(x), scale=sem(x))
    ub = ci[1]
    mean_ci = f"{round(np.mean(x),3)} ± {round(ub - np.mean(x),3)}"
    
    return mean_ci

def get_confidence_interval(df_compare, target, ci=0.95):
    value = list(df_compare.columns)[-1]
    compare_mean_ci = df_compare.groupby(target)[value].apply(calc_ci)

    return compare_mean_ci

def do_mwu(df_compare, target):
    value = list(df_compare.columns)[-1]
    array_severe =  df_compare[df_compare[target].isin(["severe"])][value].to_numpy()
    array_mild = df_compare[df_compare[target].isin(["mild"])][value].to_numpy()
    u_stat, mwu_p = mannwhitneyu(array_severe, array_mild)
    result_mwu = pd.DataFrame({"Category": [value], "U-stat": [u_stat] , "MWU_P": [mwu_p]})

    return result_mwu

def do_kruskal_wallis(df_compare, target):
    value = list(df_compare.columns)[-1]
    array_severe = df_compare[df_compare[target].isin(["severe"])][value].to_numpy()
    array_mild = df_compare[df_compare[target].isin(["mild"])][value].to_numpy()
    array_recov = df_compare[df_compare[target].isin(["Convalescent"])][value].to_numpy()
    h_stat, kruskal_p = kruskal(array_severe, array_mild, array_recov)
    result_kruskal = pd.DataFrame({"Category": [value], "H-stat": [h_stat] , "Kruskal_P": [kruskal_p]})

    return result_kruskal

def do_posthoc_dunn(df_compare, target):
    # https://blog.4dcu.be/programming/2021/12/30/Posthoc-Statannotations.html
    value = list(df_compare.columns)[-1]
    df_dunn = posthoc_dunn(df_compare, val_col=value, group_col=target, p_adjust="fdr_bh")
    remove = np.tril(np.ones(df_dunn.shape), k=0).astype("bool")
    df_dunn[remove] = np.nan
    result_dunn = df_dunn.melt(ignore_index=False).reset_index().dropna()
    result_dunn.columns = ["index", "variable", "Post_hoc_P"]
    result_dunn["category"] = value

    return result_dunn

def draw_histogram(list_samples, target, colname, tag, outdir):
    outfigdir = os.path.join(outdir, colname)
    os.makedirs(outfigdir, exist_ok=True)
    sns.histplot(list_samples)
    plt.title(f"{target}_{'_'.join(sorted(list(set(list(map(str, [tag]))))))}")
    plt.xlabel("value")
    plt.ylabel("count")
    plt.tight_layout()
    plt.savefig(f"{outfigdir}/histplot_{target}_{'_'.join(sorted(list(set(list(map(str, [tag]))))))}.png", dpi=300)
    # plt.show()
    plt.close()

def draw_boxplot(df_compare, target, result_pair, dict_idx_rename, outdir):
    colname = list(df_compare.columns)[-1]
    if colname not in dict_idx_rename.keys():
        new_colname = colname
    else:
        new_colname = dict_idx_rename[colname]
    
    plt.figure(figsize=(3,5))
    flierprops = dict(marker='o', markerfacecolor='None', markersize=5,  markeredgecolor='black')
    if len(df_compare[target].unique()) == 2:
        order = sorted(list(df_compare[target].unique()), reverse=True)
        dict_palette = {"severe":"#D68C78", "mild":"#aaaaaa", "Convalescent":"#75A7C3"}
        palette = list(map(lambda x: dict_palette[x], order))

    elif len(df_compare[target].unique()) > 2:
        order = ["severe", "mild", "Convalescent"]
        palette = {"severe":"#D68C78", "mild":"#aaaaaa", "Convalescent":"#75A7C3"}

    ax = sns.boxplot(data=df_compare, x=target, y=colname, palette=palette, order=order, flierprops=flierprops)
    ax.set_xlabel(target,fontsize=12)
    ax.set_ylabel(new_colname,fontsize=12)

    if len(df_compare[target].unique()) == 2:
        pairs = [("severe", "mild")]   
        pvalues = [i[1]["MWU_P"] for i in result_pair.iterrows()]

    elif len(df_compare[target].unique()) > 2:
        pairs = [(i[1]["index"], i[1]["variable"]) for i in result_pair.iterrows()]   
        pvalues = [i[1]["Post_hoc_P"] for i in result_pair.iterrows()]

    add_stat_annotation(ax, data=df_compare, x=target, y=colname, order=sorted(order),
                    box_pairs=pairs,
                    perform_stat_test=False, pvalues=pvalues,
                    test=None, text_format=f'star', loc='inside', verbose=2)

    list_target = df_compare[target].to_list()
    count_target = Counter(list_target)
    new_xticklabels = list()
    new_legend = list()
    for ind_xtick, count_pair in enumerate(count_target.items(), start = 1):
        key, cnt = count_pair
        shortkey = key[0].upper()
        longkey = key[0].upper() + key[1:].lower()
        xticklabel = f"{shortkey}\n(N={cnt})"
        new_xticklabels.append(xticklabel)
    plt.xticks(ticks=list(range(0, len(new_xticklabels))), labels=new_xticklabels, fontsize=10)
    plt.ylim(0, 100)
    plt.tight_layout()
    plt.savefig(f"{outdir}/boxplot_{target}_vs_{colname}.png", dpi=300)
    plt.show()
    plt.close()
        



# %%
if __name__ == "__main__":
    path_excel = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Resources/Data/infectomics_CRF_20230410_edit.xlsx"
    target = "Severity"
    case_tag = "Severe"
    control_tag = "Mild"
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
    fontpath = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Resources/Data/NanumGothicBold.ttf"
    outdir = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/9_clinical/figure"
    dict_palette = {"severe":"#D68C78", "mild":"#aaaaaa", "Convalescent":"#75A7C3"}
    dict_idx_rename = {"흉부_X선_검사_실시_여부": "ChestXrayExamination", "심전도_검사_검사_실시_여부": "ECGExamination", "심전도_검사검사결과": "ECGResults", "코로나_백신_유0_0._무__1._유": "VaccinationStatus", "CT_검사_실시_여부": "CTExamination", "2._현재_담배를_피우고_계십니까?": "SmokingStatus", "3-2._하루_평균_흡연_개피_수_(개피)": "DailyTobaccoConsumption", "후유증유무": "", "3-1._흡연_시작_연령(만_나이)": "AgeStartSmoking", "진단시_CT(PCR)__N": "CT(PCR)_N", "진단시_CT(PCR)__R": "CT(PCR)_R", "진단시_CT(PCR)_값____E": "CT(PCR)_E", "진단시_CT(PCR)_값___ORF1ab": "CT(PCR)_ORF1ab", "호흡곤란정도": "SelfReportedDyspnea", "3-3._흡연_총_기간_(년)": "SmokingYears", "1._살아오는_동안_총_5갑_이상의_흡연": "EverSmoked5PacksTobacco", "Anti Cardiolipin Ab lgG.1": "Anti beta 2-GP 1 IgG level", "Anti beta 2-GP 1 lgM.1": "Anti beta 2-GP 1 IgM level", "Anti Cardiolipin Ab lgG.1": "Anti Cardiolipin Ab lgG level", "Anti-Phospholipid lgG.1": "Anti-Phospholipid lgG level", "Anti Cardiolipin Ab lgM.1": "Anti Cardiolipin Ab lgM level"}
    os.makedirs(outdir, exist_ok=True)

# %%
    main(path_excel, target, list_drop_samples, fontpath, dict_palette, dict_idx_rename, outdir)


# %%
# if __name__ == "__main__":
#     path_excel = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Resources/Data/infectomics_CRF_20230410_edit.xlsx"
#     target = "중증도분류"
#     case_tag = ["3", "4"]
#     control_tag = ["1", "2"]
#     fontpath = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Resources/Data/NanumGothicBold.ttf"
#     outdir = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Results/Infectomics_COVID-19_CRF"
#     main(path_excel, target, case_tag, control_tag, fontpath, outdir)


# %% [severity]
# def main(crf_excel, crf_txt):
#     df_crf = pd.read_excel(crf_excel, engine="openpyxl",sheet_name="중증도 그룹")
#     df_crf = df_crf.iloc[:, :2]
#     df_crf = df_crf.rename(columns={"ID": "SampleID", "Sevierity_group": target})
#     df_crf["SampleID"] = df_crf["SampleID"].apply(lambda x: "-".join(x.split("-")[:-1]))
#     df_crf["Severity_Binary"] = df_crf[target].apply(lambda x: "Severe" if int(x) > 2 else "Mild")
#     df_crf.to_csv(crf_txt, sep="\t", index=False)


# if __name__ == "__main__":
#     crf_excel = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Resources/Data/Infectomics_CRF_20230915.xlsx"
#     crf_txt = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Resources/Data/Infectomics_CRF_202300915.tsv"
#     main(crf_excel, crf_txt)

# %% [deprecated]
# def main(path_excel, target, case_tag, control_tag, list_drop_samples, fontpath, outdir):
#     set_korean_font(fontpath)
    
#     df_crf_raw = pd.read_excel(path_excel, engine="openpyxl", sheet_name="21-22등록 대상자_modify", skiprows=1)
#     df_crf_raw = correct_header_dataframe(df_crf_raw)
#     df_crf_raw_selec = select_dataframe(df_crf_raw, num=2)
#     df_crf_raw_filt = filter_dataframe(df_crf_raw_selec, list_drop_samples)
#     df_crf_raw_add_sev = define_severity(df_crf_raw_filt)
#     dict_sample_phase = get_dict_sample_visit(df_crf_raw_add_sev)
#     list_sample_pattern = get_sample_pattern_need_visit(dict_sample_phase)
#     df_crf = extract_dataframe(df_crf_raw_add_sev, list_sample_pattern)
#     list_df_compare = get_clinical_assoc_severity(df_crf)

#     list_sample_size = list()
#     list_mean_ci = list()
#     list_pair = list()
#     list_kruskal = list()
#     list_dunn = list()
#     for df_compare in list_df_compare:
#         df_compare_copy = df_compare.copy()
#         value = list(df_compare_copy.columns)[-1]
#         try:
#             sample_size = get_sample_size(df_compare_copy, target)
#             list_sample_size.append(sample_size)

#             mean_ci = get_confidence_interval(df_compare_copy, target, ci=0.95)
#             list_mean_ci.append(mean_ci)

#             if len(df_compare_copy[target].unique()) == 2:
#                 result_pair = do_mwu(df_compare_copy, target)
#                 list_pair.append(result_pair)
#                 # draw_boxplot(df_compare_copy, target, result_pair, outdir)

#             elif len(df_compare_copy[target].unique()) > 2:
#                 result_kruskal = do_kruskal_wallis(df_compare_copy, target)
#                 list_kruskal.append(result_kruskal)

#                 result_dunn = do_posthoc_dunn(df_compare_copy, target)
#                 list_dunn.append(result_dunn)
#                 # draw_boxplot(df_compare_copy, target, result_dunn, outdir)
            
#             else:
#                 print(value)
#                 break

#         except Exception as e:
#             print(value, e)
    
#     df_n = pd.concat(list_sample_size, axis=1)
#     df_n_dropna = df_n.dropna(axis=1)
#     df_n_T = df_n_dropna.T.reset_index(drop=False).rename(columns={"index": "Category"})    

#     df_ci = pd.concat(list_mean_ci, axis=1)
#     df_ci_dropna = df_ci.dropna(axis=1)
#     df_ci_T = df_ci_dropna.T.reset_index(drop=False).rename(columns={"index": "Category"})
    
#     overall_kruskal = pd.concat(list_kruskal)
#     overall_kruskal_dropna = overall_kruskal.dropna()
#     overall_kruskal_dropna["Significant"] = fdrcorrection(overall_kruskal_dropna["Kruskal_P"])[0]
#     overall_kruskal_dropna["Padj"] = fdrcorrection(overall_kruskal_dropna["Kruskal_P"])[1]
#     overall_kruskal_sig = overall_kruskal_dropna[overall_kruskal_dropna["Significant"]==True]

#     overall_dunn = pd.concat(list_dunn)
#     overall_dunn_dropna = overall_dunn.dropna()
#     overall_dunn_dropna["Significant"] = fdrcorrection(overall_dunn_dropna["Post_hoc_P"])[0]
#     overall_dunn_dropna["Padj"] = fdrcorrection(overall_dunn_dropna["Post_hoc_P"])[1]
#     overall_dunn_sig = overall_dunn_dropna[overall_dunn_dropna["Significant"]==True]  

