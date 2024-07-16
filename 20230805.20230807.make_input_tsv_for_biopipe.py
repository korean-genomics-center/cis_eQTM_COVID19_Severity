import os
import glob
import pandas as pd
from functools import reduce


def main(input_dir, 
         output_dir, 
         results_dir, 
         list_samples, 
         include,
         exclude, 
         run_date):

    check_if_input_dir_exist(input_dir)
    make_directories(results_dir, output_dir)
    make_input = MakeBioPipeInput(input_dir, 
                                  output_dir, 
                                  results_dir, 
                                  list_samples,
                                  include,
                                  exclude,
                                  run_date)
    make_input.provide_full_path_info()


class MakeBioPipeInput:

    def __init__(self, 
                input_dir, 
                output_dir, 
                results_dir, 
                list_samples, 
                include,
                exclude, 
                run_date):

        self.input_dir = input_dir
        self.output_dir = output_dir
        self.results_dir = results_dir
        self.list_samples = list_samples
        self.include = include
        self.exclude = exclude
        self.run_date = run_date
    
    def provide_full_path_info(self):
        print("start...")
        self.__get_fastp_input()
        self.__make_tsv_fastp_input()
        self.__get_star_input()
        self.__make_tsv_star_input()
        print("end...")
    
    def __get_fastp_input(self):
        path_input = self.input_dir
        path_result = self.results_dir
        path_fastp_result = f"{path_result}/1_fastp"
        list_files_abspath = self.search_files_path_source(path_input, "gz")

        if self.include != None:
            list_files_abspath = list(filter(lambda x : self.include in x, list_files_abspath))
        elif self.exclude != None:
            list_files_abspath = list(filter(lambda x : self.exclude not in x, list_files_abspath))

        list_fastp_input_read1_abspath = [file_abspath for file_abspath in list_files_abspath if "_R1" in file_abspath]
        list_fastp_input_read2_abspath = [file_abspath for file_abspath in list_files_abspath if "_R2" in file_abspath]
        list_fastp_input_read1 = [os.path.basename(file_abspath) for file_abspath in list_fastp_input_read1_abspath] 
        list_fastp_input_read2 = [os.path.basename(file_abspath) for file_abspath in list_fastp_input_read2_abspath]

        list_sample_id = ["_".join(file.split("_")[:-1]) for file in list_fastp_input_read1]

        list_fastp_output_read1 = [self.add_string_element_list(elem, sep=".", index=1, string="trimmed") for elem in list_fastp_input_read1]
        list_fastp_output_read2 = [self.add_string_element_list(elem, sep=".", index=1, string="trimmed") for elem in list_fastp_input_read2]
        list_fastp_output_read1_abspath = [f"{path_fastp_result}/{sample_id}/{file}" for sample_id, file in zip(list_sample_id, list_fastp_output_read1)]
        list_fastp_output_read2_abspath = [f"{path_fastp_result}/{sample_id}/{file}" for sample_id, file in zip(list_sample_id, list_fastp_output_read2)]
        list_fastp_report_name_abspath = [fastp_output_abspath.split(".")[0][:-3] for fastp_output_abspath in list_fastp_output_read1_abspath]
        list_fastp_html_report_name_abspath = [f"{fastp_report_name_abspath}.html" for fastp_report_name_abspath in list_fastp_report_name_abspath]
        list_fastp_json_report_name_abspath = [f"{fastp_report_name_abspath}.json" for fastp_report_name_abspath in list_fastp_report_name_abspath]

        self.path_fastp_result = path_fastp_result
        self.list_sample_id = list_sample_id
        self.list_fastp_input_read1_abspath = list_fastp_input_read1_abspath
        self.list_fastp_input_read2_abspath = list_fastp_input_read2_abspath
        self.list_fastp_output_read1_abspath = list_fastp_output_read1_abspath
        self.list_fastp_output_read2_abspath = list_fastp_output_read2_abspath
        self.list_fastp_html_report_name_abspath = list_fastp_html_report_name_abspath
        self.list_fastp_json_report_name_abspath = list_fastp_json_report_name_abspath

    def __make_tsv_fastp_input(self):
        path_output = self.output_dir
        path_output_tsv_fastp_input = f"{path_output}/fastp_input_for_biopipe_{self.run_date}.tsv"
        dict_fastp_input = {"SampleID":self.list_sample_id,
                         "FastpInputRead1":self.list_fastp_input_read1_abspath,
                         "FastpInputRead2":self.list_fastp_input_read2_abspath,
                         "FastpOutputRead1":self.list_fastp_output_read1_abspath,
                         "FastpOutputRead2":self.list_fastp_output_read2_abspath, 
                         "FastpHtmlReportName":self.list_fastp_html_report_name_abspath,
                         "FastpJsonReportName":self.list_fastp_json_report_name_abspath}

        df_fastp_input = pd.DataFrame.from_dict(dict_fastp_input, orient="columns")
        df_fastp_input = df_fastp_input.sort_values(by=["SampleID"])
        if self.list_samples != []:
            df_fastp_input = df_fastp_input.set_index("SampleID")
            df_fastp_input = df_fastp_input.loc[self.list_samples,:]
        df_fastp_input.to_csv(path_output_tsv_fastp_input, sep="\t", index=False)
        return df_fastp_input

    def __get_star_input(self):
        path_input = self.input_dir
        list_files_abspath = self.search_files_path_source(path_input, "gz")

        if self.include != None:
            list_files_abspath = list(filter(lambda x : self.include in x, list_files_abspath))
        elif self.exclude != None:
            list_files_abspath = list(filter(lambda x : self.exclude not in x, list_files_abspath))


        list_sample_id_abspath = [os.path.basename(file) for file in list_files_abspath]
        list_sample_id = ["_".join(file.split("_")[:-2]) for file in list_sample_id_abspath]
        
        list_files = [os.path.basename(file_abspath) for file_abspath in list_files_abspath]
        list_star_input = [self.add_string_element_list(elem, sep=".", index=1, string="trimmed") for elem in list_files]
        list_star_input = [f"{self.path_fastp_result}/{sample_id}/{star_input}" for sample_id, star_input in zip(list_sample_id, list_star_input)]
        list_sample_id_star_input_pair = self.__pair_two_lists(list_sample_id, list_star_input)
        dict_sample_id_star_input_pair= self.__groupby_tuple(list_sample_id_star_input_pair)
        self.dict_sample_id_star_input_pair = dict_sample_id_star_input_pair


    def __make_tsv_star_input(self):
        path_result = self.results_dir
        path_star_result = f"{path_result}/2_star"
        path_output = self.output_dir
        path_output_tsv_star = f"{path_output}/star_input_for_biopipe_{self.run_date}.tsv"
        dict_star_input_read1 = {}
        dict_star_input_read2 = {}
        df_star_output_prefix = {}
        for k, v in self.dict_sample_id_star_input_pair.items():
            list_r1 = []
            list_r2 = []
            for elem in v:
                if "_R1" in elem:
                    list_r1.append(elem)
                elif "_R2" in elem:
                    list_r2.append(elem)
                else:
                    exit(f"Error: Check your fastq files. {elem} has no R1 or R2 signatures for generating STAR input.")
            str_r1_cs = ",".join(list_r1)
            str_r2_cs = ",".join(list_r2)
            dict_star_input_read1[k] = str_r1_cs
            dict_star_input_read2[k] = str_r2_cs
            df_star_output_prefix[k] = f"{path_star_result}/{k}/{k}_"
        df_star_input_read1 = pd.DataFrame(dict_star_input_read1.items(), columns = ["SampleID", "STARInputRead1"])
        df_star_input_read2 = pd.DataFrame(dict_star_input_read2.items(), columns = ["SampleID", "STARInputRead2"])
        df_star_output_prefix = pd.DataFrame(df_star_output_prefix.items(), columns = ["SampleID", "STAROutputPrefix"])
        list_dfs = [df_star_input_read1, df_star_input_read2, df_star_output_prefix]
        df_star_input_total = reduce(lambda left, right: pd.merge(left, right, how="inner", on="SampleID"), list_dfs)
        df_star_input_total = df_star_input_total.sort_values(by=["SampleID"])
        if self.list_samples != []:
            df_star_input_total = df_star_input_total.set_index("SampleID")
            df_star_input_total = df_star_input_total.loc[self.list_samples,:]
        df_star_input_total.to_csv(path_output_tsv_star, sep="\t", index=False)
        return df_star_input_total
    
    @staticmethod
    def search_files_path_source(path_input, exclude):
        list_src_files = []
        for src in glob.glob(f'{path_input}/**/*{exclude}', recursive=True):

            list_src_files.append(src)
        return list_src_files

    @staticmethod
    def add_string_element_list(elem, sep, index, string):
        list_elem = elem.split(sep)
        list_elem.insert(index, string)
        elem = sep.join(list_elem)
        return elem

    @staticmethod
    def __pair_two_lists(list1, list2):   
        list_zipped_tuples = list(zip(list1, list2))
        return list_zipped_tuples
    
    @staticmethod
    def __groupby_tuple(list_zipped_tuples):
        dict_groupby = {}
        for k, v in list_zipped_tuples:
            if k not in dict_groupby:
                dict_groupby[k] = [v]
            else:
                dict_groupby[k].append(v)
        return dict_groupby

def get_path_both_input_output(input_dir, output_dir):
    path_input = os.path.abspath(input_dir)
    path_output = os.path.abspath(output_dir)
    return path_input, path_output

def get_group_name_samples(path):
    group_name = os.path.basename(path)
    return group_name

def get_files_under_directory(path):
    list_files = os.listdir(path)
    return list_files

def get_abs_path_files(path, file):
    file_abs_path = os.path.join(path,file)
    return file_abs_path

def get_unique_list(target_list):
    unique_list = list(set(target_list))
    return unique_list

def get_sorted_list(target_list):
    sorted_list = sorted(target_list)
    return sorted_list

def make_directories(results_dir, output_dir):
    os.makedirs(results_dir,exist_ok=True)
    os.makedirs(output_dir,exist_ok=True)

def check_if_input_dir_exist(input_dir):
    if os.path.exists(input_dir):
        pass
    else:
        print(f"'{input_dir}' does not exist!")
        raise FileNotFoundError
            

if __name__ == "__main__":
    run_date = "20230805"
    input_dir = "/BiO/Store/Vaccination"
    output_dir = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/kyungwhan1998/Vaccination/Analysis/BioPipeConfig/BioPipeInput"
    results_dir = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/kyungwhan1998/Vaccination/Results/{run_date}"
    list_samples = []
    include = ""
    exclude = ""
    main(input_dir, output_dir, results_dir, list_samples, include, exclude, run_date)
    
