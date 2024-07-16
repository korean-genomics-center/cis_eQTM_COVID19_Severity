# %%
import os

import numpy as np
import pandas as pd


def make_master_table_metadata(metadata, outdir, outfilename, dict_sampledir, sample_column, list_target_columns, list_new_colnames):
    print("start...")
    df_metadata = read_metadata(metadata)
    df_metadata = df_metadata[list_target_columns]
    df_metadata.columns = list_new_colnames
    df_metadata_add_ox = add_analysis_ox_column(df_metadata, dict_sampledir, sample_column)
    df_metadata_add_ox_vis_order = parse_visit_order(df_metadata_add_ox, sample_column)
    df_master_table = fillna_df_master_table(df_metadata_add_ox_vis_order)
    outfilepath = os.path.join(outdir, outfilename)
    df_master_table.to_csv(outfilepath, sep="\t", index=False)
    print("end...")

def read_metadata(metadata):
    df = pd.read_excel(metadata, skiprows=2, engine="openpyxl", sheet_name="RNA-seq")
    df = df.iloc[:, 2:30]
    df_dropna = df.dropna(axis = 0, how = 'all')
    df_metadata = df_dropna.astype(str)
    
    return df_metadata

def __get_sample_list_from_project_directory(path):
    try:
        list_file = os.listdir(path)
    except FileNotFoundError:
        return list()
    return list_file

def __check_is_sample_in_master_table(master_table, list_sample, sample_column):
    list_sample_in_master = master_table.loc[:, sample_column].to_list()
    list_bool_is_sample_in = list()
    for sample in list_sample:
        answer_bool = sample in list_sample_in_master
        list_bool_is_sample_in.append(answer_bool)
    return list_bool_is_sample_in
    
def __divide_list_by_bool(list_data, list_bool):
    list_sample_in = list()
    list_sample_not_in = list()
    for idx, data_bool in enumerate(list_bool):
        if data_bool:
            list_sample_in.append(list_data[idx])
        else:
            list_sample_not_in.append(list_data[idx])
    return list_sample_in, list_sample_not_in

def __get_binary_list_is_in(df_master_table, sample_column, list_data_compare_bool):
    list_compared_binary = list()
    list_data_original = df_master_table[sample_column].to_list()
    for data_original in list_data_original:
        if data_original in list_data_compare_bool:
            list_compared_binary.append("1")
        else:
            list_compared_binary.append("0")
    return list_compared_binary

def __add_not_in_sample_column(sample_column, sample_not_in, analysis_type, df_master_table, dict_sampledir):
    list_project_name = list(dict_sampledir.keys())
    list_sample_column = list(df_master_table.columns)
    list_not_in_sample_column = ["-"] * len(list_sample_column)
    add_series = pd.Series(index =list_sample_column, data=list_not_in_sample_column)
    add_series[sample_column] = sample_not_in
    for project_type in list_project_name:
        if project_type in list_sample_column:
            add_series[project_type] = 0
    add_series[analysis_type] = 1
    df_master_table = df_master_table.append(add_series, ignore_index = True)
    return df_master_table

def add_analysis_ox_column(df_master_table, dict_sampledir, sample_column):
    list_data_analysis = list(dict_sampledir.keys())
    for analysis_type in list_data_analysis:
        path_analysis = dict_sampledir[analysis_type]
        list_sample_from_project_dir = __get_sample_list_from_project_directory(path_analysis)
        list_bool_is_sample_in = __check_is_sample_in_master_table(df_master_table, list_sample_from_project_dir, sample_column)
        list_sample_in, list_sample_not_in = __divide_list_by_bool(list_sample_from_project_dir, list_bool_is_sample_in)
        list_col_add = __get_binary_list_is_in(df_master_table, sample_column, list_sample_in)
        df_master_table[analysis_type] = list_col_add
        for sample_not_in in list_sample_not_in:
            df_master_table = __add_not_in_sample_column(sample_column, sample_not_in, analysis_type, df_master_table, dict_sampledir)
    return df_master_table

def __get_visit(data):
    import re
    sample_parser = re.compile("-[A-Z][0-9]+")
    list_find_visit_order = sample_parser.findall(data)
    
    visit = '-'
    if len(list_find_visit_order) > 0:
        tmp_visit = list_find_visit_order[-1][1:]
        visit = tmp_visit.replace("V", "Visit").replace("A1", "Control").replace("L1", "LongCOVID")
    
    if data.startswith("C19-R"):
        visit = "Recover"

    return visit

def parse_data(df_master, column, func):
    list_data = df_master[column].to_list()
    list_parsed = list(map(func, list_data))
    return list_parsed

def parse_visit_order(df_master, sample_column):
    df_master["Visit_order"] = parse_data(df_master, sample_column, __get_visit)
    return df_master

def fillna_df_master_table(df):
    df = df.fillna("-")
    return df

if __name__ == "__main__":
    metadata = "/BiO/Research/RNASeqReference/Vaccination/Resources/Data/KOGIC_Metadata_RNA_20230902_version1.xlsx"
    outdir = "/BiO/Research/RNASeqReference/Vaccination/Resources/Data"
    outfilename = "COVID19_master_table_20231007.txt"
    dict_sampledir = dict()
    dict_sampledir["Vaccination"] = "/BiO/Research/RNASeqReference/Vaccination/Results/20230813/3_rsem"
    dict_sampledir["HealthyControl"] = "/BiO/Research/RNASeqReference/HealthyControl/3_rsem"
    dict_sampledir["Confirmed"] = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/kyungwhan1998/COVID19Infected/Results/3_rsem/Confirmed"
    dict_sampledir["Recovered"] = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/kyungwhan1998/COVID19Infected/Results/3_rsem/Recovered"
    list_target_columns = ["Project-ID", "Project-ID-Alias", "Sample_Age", "Sample_Sex", "Project_Name", "Sample_Trait", "Sequencing_Type_Platform"]
    list_new_colnames = ["Project_ID", "Project_ID_Alias", "Sample_Age", "Sample_Sex", "Project_Name", "Sample_Trait", "Sequencing_Type_Platform"]
    sample_column = "Project_ID_Alias"
    make_master_table_metadata(metadata, outdir, outfilename, dict_sampledir, sample_column, list_target_columns, list_new_colnames)