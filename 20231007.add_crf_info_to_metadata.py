import os

import pandas as pd


def main(path_master, path_crf):    
    df_master_table = pd.read_csv(path_master, sep="\t")
    df_crf = pd.read_csv(path_crf, sep="\t")
    df_crf = df_crf.rename(columns={"SampleID": "Project_ID"})
    df_crf = df_crf[["Project_ID", "Severity"]]
    df_master_table_added_crf = pd.merge(df_master_table, df_crf, how="outer", on="Project_ID")
    df_master_table_added_crf = df_master_table_added_crf.fillna("-")
    df_master_table_added_crf.to_csv(outfilepath, sep="\t", index=False)

if __name__ == "__main__":
    path_master = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Resources/Data/COVID19_master_table_20231007.txt"
    path_crf = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/9_clinical/Infectomics_Severity_Information_20231106.tsv"
    outfilepath = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Resources/Data/COVID19_master_table_added_CRF_20231007.txt"
    main(path_master, path_crf)