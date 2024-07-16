import os
import glob
import sys
import shutil
import subprocess


def main(dir_input, dir_output, path_hashtable, pattern):
    check_if_dir_input_dir_output_same(dir_input, dir_output)
    list_src_files = get_list_src_files(dir_input, pattern)
    hashtable = get_hashtable(path_hashtable)
    dict_src_dest = get_dict_pair_src_dest_files(dir_output, list_src_files, hashtable)
    make_symlinks(dict_src_dest)

def check_if_dir_input_dir_output_same(dir_input, dir_output):
    if dir_input == dir_output:
        print("Critical Error: DO NOT set the input and output directories the same!")
        raise Exception
    else:
        print("Appropriate paths! Proceeding...")

def get_list_src_files(dir_input, pattern):
    list_src_files = glob.glob(f'{dir_input}/**/*{pattern}', recursive=True)

    return list_src_files

def get_hashtable(path_hashtable):
    hashtable = {}
    with open (path_hashtable, mode="r") as fr:
        _skiprow = fr.readline()
        for line in fr:
            record = line.rstrip().split("\t")
            ori_id = record[0]
            new_id = record[1]

            hashtable[ori_id] = new_id
    
    return hashtable

def __make_dest_dir(dest_dir):
    os.makedirs(dest_dir, exist_ok=True)

def __get_filename_prefix(filename):
    prefix = "_".join(os.path.basename(filename).split("_")[:-1])

    return prefix

def __get_filename_suffix(filename):
    suffix = os.path.basename(filename).split("_")[-1]

    return suffix

def get_dict_pair_src_dest_files(dir_output, list_src_files, hashtable):
    dict_src_dest = dict()
    list_prefix = list(map(__get_filename_prefix, list_src_files))
    list_suffix = list(map(__get_filename_suffix, list_src_files))
    for src_prefix, src_suffix, src_file in zip(list_prefix, list_suffix, list_src_files):
        if src_prefix in hashtable.keys():
            dest_prefix = hashtable[src_prefix]
            dest_suffix = f"_R{src_suffix}"
            new_subject_id = dest_prefix.split("_")[0]
            dest_dir = os.path.join(dir_output, new_subject_id)
            __make_dest_dir(dest_dir)
            dest_file = dest_prefix + dest_suffix
            dest_file = os.path.join(dest_dir, dest_file)
            dict_src_dest[src_file] = dest_file

    return dict_src_dest
    
def make_symlinks(dict_src_dest):
    for src, des in dict_src_dest.items():
        if src == des:
            raise Exception
        else:
            if not os.path.exists(src):
                raise FileNotFoundError
            else:
                if os.path.exists(des):
                    raise Exception
                else:
                    subprocess.run(f"ln -s {src} {des}", shell=True)

if __name__ == "__main__":
    dir_input = "/BiO/Store/RawData_BMS_2022/RawData_BMS_2022"
    dir_output = "/BiO/Store/Vaccination"
    os.makedirs(dir_output, exist_ok=True)
    path_hashtable = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/kyungwhan1998/Vaccination/Resources/hashtable_u10k_00017_sample_vacc_metadata.txt"
    pattern = ".fq.gz"
    main(dir_input, dir_output, path_hashtable, pattern)