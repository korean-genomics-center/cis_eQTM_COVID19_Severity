import os
import pandas as pd
import numpy as np
from scipy.stats import f_oneway, kruskal
from scikit_posthocs import posthoc_tukey, posthoc_dunn
from statannot import add_stat_annotation
import matplotlib.pyplot as plt
import seaborn as sns

def main(path_analysis, gene):
    df_analysis = read_input_files(path_analysis, gene)
    result_kruskal = do_kruskal_wallis(df_analysis, gene)
    print(result_kruskal)
    result_dunn = do_posthoc_dunn(df_analysis, gene)
    print(result_dunn)
    draw_boxplot(df_analysis, result_dunn, gene)
    result_anova = do_one_way_anova(df_analysis, gene)
    print(result_anova)
    result_tukey = do_posthoc_tukey_hsd(df_analysis, gene)
    print(result_tukey)

def read_input_files(path_analysis, gene):
    

def do_kruskal_wallis(df_analysis, gene):
    df_analysis = df_analysis.pivot(index='sample', columns='severity',values=gene)
    asymptomatic = [x for x in df_analysis.asymptomatic.to_list() if str(x) != "nan"]
    mild = [x for x in df_analysis.mild.to_list() if str(x) != "nan"]
    moderate =  [x for x in df_analysis.moderate.to_list() if str(x) != "nan"]
    severe =  [x for x in df_analysis.severe.to_list() if str(x) != "nan"]
    h_stat, kruskal_p = kruskal(asymptomatic, mild, moderate, severe)
    result_kruskal = pd.DataFrame({"H-stat": [h_stat] , "Kruskal_P": [kruskal_p]})
    return result_kruskal

def do_one_way_anova(df_analysis, gene):
    df_analysis = df_analysis.pivot(index='sample', columns='severity',values=gene)
    asymptomatic = [x for x in df_analysis.asymptomatic.to_list() if str(x) != "nan"]
    mild = [x for x in df_analysis.mild.to_list() if str(x) != "nan"]
    moderate =  [x for x in df_analysis.moderate.to_list() if str(x) != "nan"]
    severe =  [x for x in df_analysis.severe.to_list() if str(x) != "nan"]
    f_stat, anova_p = f_oneway(asymptomatic, mild, moderate, severe)
    result_anova = pd.DataFrame({"F-stat": [f_stat] , "ANOVA_P": [anova_p]})
    return result_anova

def do_posthoc_dunn(df_analysis, gene):
    df_dunn = posthoc_dunn(df_analysis, val_col=gene, group_col="severity", p_adjust="fdr_bh")
    remove = np.tril(np.ones(df_dunn.shape), k=0).astype("bool")
    df_dunn[remove] = np.nan
    result_dunn = df_dunn.melt(ignore_index=False).reset_index().dropna()
    # result_dunn = result_dunn[result_dunn.value < 0.05]
    result_dunn.columns = ["group1", "group2", "Post_hoc_P"]
    return result_dunn

def do_posthoc_tukey_hsd(df_analysis, gene):
    df_tukey = posthoc_tukey(df_analysis, val_col=gene, group_col="severity")
    remove = np.tril(np.ones(df_tukey.shape), k=0).astype("bool")
    df_tukey[remove] = np.nan
    result_tukey = df_tukey.melt(ignore_index=False).reset_index().dropna()
    # result_tukey = result_tukey[result_tukey.value < 0.05]
    result_tukey.columns = ["group1", "group2", "Post_hoc_P"]
    return result_tukey

def draw_boxplot(df_analysis, result_posthoc, gene):
    order = ["asymptomatic", "mild", "moderate", "severe"]
    ax = sns.boxplot(data=df_analysis, x="severity", y=gene, order=order)
    n_asym = df_analysis[df_analysis["severity"] == 1].shape[0]
    n_mild = df_analysis[df_analysis["severity"] == 2].shape[0]
    n_mod = df_analysis[df_analysis["severity"] == 3].shape[0]
    n_sev = df_analysis[df_analysis["severity"] == 4].shape[0]
    xticklabel = [f"asymptomatic \n (N={n_asym})", f"mild \n  (N={n_mild})", f"moderate \n  (N={n_mod})", f"severe \n  (N={n_sev})"]
    ax.set_xticklabels(xticklabel)
    pairs = [(i[1]["group1"], i[1]["group2"]) for i in result_posthoc.iterrows()]   
    pvalues = [i[1]["Post_hoc_P"] for i in result_posthoc.iterrows()]
    add_stat_annotation(ax, data=df_analysis, x="severity", y=gene, order=order,
                    box_pairs=pairs,
                    perform_stat_test=False, pvalues=pvalues,
                    test=None, text_format=f'star', loc='inside', verbose=2)
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"./{gene}_expression_by_severity_boxplot.png")

if __name__ == "__main__":
    path_analysis = ""
    gene = "LY6E"
    main(path_analysis, gene)