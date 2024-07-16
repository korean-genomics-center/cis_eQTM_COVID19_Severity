# %%
import os
import pandas as pd

path_recovered1 = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/4_expmtx/Recovered1/expression_matrix_genes.results_TPM.tsv"
path_recovered2 = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/4_expmtx/Recovered2/expression_matrix_genes.results_TPM.tsv"

df_recov1 = pd.read_csv(path_recovered1, sep="\t")
df_recov2 = pd.read_csv(path_recovered2, sep="\t")

df_recovall = pd.merge(df_recov1, df_recov2, how="inner", on="ID")
df_recovall = df_recovall.set_index("ID")
df_recovall = df_recovall[sorted(list(df_recovall.columns))]
df_recovall.to_csv("/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/4_expmtx/RecoveredAll/expression_matrix_genes.results_TPM.tsv", sep="\t", index=True)

# %%

path_confirmed = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/4_expmtx/Confirmed/expression_matrix_genes.results_TPM.tsv"
path_recoveredall = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/4_expmtx/RecoveredAll/expression_matrix_genes.results_TPM.tsv"

df_conf = pd.read_csv(path_confirmed, sep="\t")
df_recov = pd.read_csv(path_recoveredall, sep="\t")

df_confrecov = pd.merge(df_conf, df_recov, how="inner", on="ID")
df_confrecov = df_confrecov.set_index("ID")
df_confrecov = df_confrecov[sorted(list(df_confrecov.columns))]
df_confrecov.to_csv("/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/4_expmtx/ConfirmedRecovered/expression_matrix_genes.results_TPM.tsv", sep="\t", index=True)