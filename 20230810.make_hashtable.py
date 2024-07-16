import pandas as pd
import os
from collections import Counter

def main(input_format, path_id_table_xlsx, path_id_table_txt, path_hashtable, sample_type, list_outheader):
    if input_format == "xlsx":
        make_txt_from_excel(path_id_table_xlsx, sample_type, path_id_table_txt)
    assert os.path.exists(path_id_table_txt), "Error: input file does not exist!"
    dict_id = make_dictionary_id(path_id_table_txt)
    dict_expt = make_dictionary_expt(dict_id)
    make_hashtable(dict_expt, path_hashtable, list_outheader)


def make_txt_from_excel(path_id_table_xlsx, sample_type, path_id_table_txt):
    df_id = pd.read_excel(path_id_table_xlsx, engine="openpyxl")
    df_id = df_id[df_id["Sample type"] == sample_type]
    df_id.to_csv(path_id_table_txt, sep="\t", index=False)

def make_dictionary_id(path_id_table_txt):
    dict_id = dict()
    with open (path_id_table_txt, mode="r") as fr:
        skiprow = fr.readline()
        for line in fr:
            record = line.rstrip("\n").split("\t")
            sample_id = record[-2]
            flowcell_id = record[-4]
            adapter = record[-1]
            flowcell_id = flowcell_id + "_L01_" + adapter
            dict_id[flowcell_id] = sample_id

    return dict_id


def make_dictionary_expt(dict_id):
    dict_expt = dict()
    list_check_dup = []
    for flowcell_id, sample_id in dict_id.items():
        
        cnt = 1
        if sample_id in list_check_dup:
            cnt +=1 
            for i in range(1, cnt+1, 1):
                dict_expt[flowcell_id] = f"{sample_id}_L0{i}"
            continue
        
        else:
            dict_expt[flowcell_id] = f"{sample_id}_L0{cnt}"

        list_check_dup.append(sample_id)
    
    # for sample, cnt in Counter(list_check_dup).items():
    #     if cnt > 1:
    #         print(sample)

    return dict_expt


def make_hashtable(dict_id, path_hashtable, list_outheader):
    if len(list_outheader) == 0:
        list_outheader = ["source", "destination"]
    with open (path_hashtable, mode="w") as fw:
        header = "\t".join(list_outheader) + "\n"
        fw.write(header)
        for key, value in dict_id.items():
            fw.write(key + "\t" + value + "\n")
    

if __name__ == "__main__":
    input_format = "txt"
    path_id_table_xlsx = None
    path_id_table_txt = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/kyungwhan1998/Vaccination/Resources/u10k_00017_sample_vacc_metadata.txt"
    path_hashtable = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/kyungwhan1998/Vaccination/Resources/hashtable_u10k_00017_sample_vacc_metadata.txt"
    sample_type = "RNA"
    list_outheader = []
    main(input_format, path_id_table_xlsx, path_id_table_txt, path_hashtable, sample_type, list_outheader)