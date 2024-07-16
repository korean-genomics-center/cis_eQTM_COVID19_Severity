#%%
import pandas as pd

list_cutoff = [1, 2, 3, 4, 5, 10]
folds = range(4)
col_check_diff = ["meth.diff", "wmeth.diff"]

path_file_format = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/DMPExtract/Infectomics.Case_Sev234.Control_Sev1.CV_4.filter.Methyl_Visit_Order_First.Weight_Severity/DMPExtract.Infectomics.Case_Sev234.Control_Sev1.CV_4.Fold_{fold}.Weight_Severity.20240314.tsv.FDR_05.tsv"
path_save_format = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/DMPExtract/Infectomics.Case_Sev234.Control_Sev1.CV_4.filter.Methyl_Visit_Order_First.Weight_Severity/DMPExtract.Infectomics.Case_Sev234.Control_Sev1.CV_4.Fold_{fold}.Weight_Severity.20240314.tsv.FDR_05.{col_diff}_abs_{cutoff}.tsv"

for fold in folds:
    table_dmp = pd.read_csv(path_file_format.format(fold = fold), sep = '\t')
    for cutoff in list_cutoff:
        for col_diff in col_check_diff:
            table_dmp_sig = table_dmp[list(map(lambda x : abs(x) > cutoff, table_dmp[col_diff]))]
            table_dmp_sig.to_csv(path_save_format.format(fold = fold, col_diff = col_diff, cutoff = cutoff), sep = '\t', index = False)
# %%
