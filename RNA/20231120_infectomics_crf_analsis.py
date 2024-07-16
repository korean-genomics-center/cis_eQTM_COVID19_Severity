# %%
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import t, ttest_ind, mannwhitneyu, kruskal
from statsmodels.stats.multitest import fdrcorrection
from matplotlib import pyplot as plt
import math

# %%
def str_to_num(x):
    try:
        x = float(x)

    except Exception: 
        pass
    
    return x

def conduct_stat_test(crf, column, alpha):
    group_yes_dex = crf[crf["DEXA"]==True][column].to_list()
    group_yes_dex = [yes for yes in group_yes_dex if str(yes) != 'nan']
    group_no_dex = crf[crf["DEXA"]==False][column].to_list()
    group_no_dex = [no for no in group_no_dex if str(no) != 'nan']
    stat, pval = ttest_ind(group_yes_dex, group_no_dex, equal_var=False)
    if pval < alpha:
        sig = 1
    else:
        sig = 0

    return (round(stat,2), round(pval,2), sig)

def draw_boxplot(crf, column):
    sns.boxplot(data=crf, x="DEXA", y=column)
    plt.show()
    # fig = p.get_figure()
    # fig.savefig(f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Results/Infectomics_COVID-19_CRF/swarm_{crf_col}.png")
    plt.close()

# %%
path_crf_excel = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Resources/Data/infectomics_CRF_20230410_edit.xlsx"

crf = pd.read_excel(path_crf_excel, engine="openpyxl", header=1)
dex_filt = (crf.loc[:,"치료제 및 주사제 투약"].str.lower().str.contains("dexa"))
crf["DEXA"] = dex_filt
list_yes_dex_samples = crf[crf["DEXA"]==True].loc[:, "Subject NO.(고유번호)"]

dict_t_result = {}
list_crf_col = list(crf.columns)
for crf_col in list_crf_col:
    list_values = crf[crf_col].to_list()
    list_values_nona = list(filter(lambda x : not isinstance(x, float) or not math.isnan(x), list_values))
    if len(list_values_nona) == 0: continue
    if type(list_values_nona[0]) != float: continue
    is_numeric = isinstance(str_to_num(list_values_nona[0]), float)
    if is_numeric:
        crf[crf_col] = crf[crf_col].dropna()
        crf[crf_col] = crf[crf_col].apply(str)
        crf[crf_col] = crf[crf_col].apply(lambda x: x.replace("..", ".").replace("R ", ""))
        crf[crf_col] = crf[crf_col].apply(lambda x: float(x.lstrip("'")))
        result = conduct_stat_test(crf, crf_col, 0.05)
        dict_t_result[crf_col] = result
df_t = pd.DataFrame.from_dict(dict_t_result, orient="index")
df_t = df_t.reset_index()
df_t.columns = ["categories","t-statistcs","pvalue","significant(yes/no)"]

df_t_sig = df_t[df_t["significant(yes/no)"]==1]
list_selec_col = df_t_sig["categories"].to_list()
for column in list_selec_col:
    draw_boxplot(crf, column)
# df_t.to_csv("/BiO/Research/Project2/Infectomics_COVID-19_Host/Results/Infectomics_COVID-19_CRF/t_test_results.tsv", sep="\t", index=False)
# %%
df_dex_sig_var = pd.read_csv("/BiO/Research/Project2/Infectomics_COVID-19_Host/Results/Infectomics_COVID-19_CRF/t_test_results.tsv", sep="\t")
df_dex_sig_var[df_dex_sig_var["significant(yes/no)"]==1]["categories"].tolist()
# %%
