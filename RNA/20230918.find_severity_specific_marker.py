# %%
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles


# %%
def main(severe_vs_mild, severe_only):
    list_severe_mild_deg = read_deg_file(severe_vs_mild)
    dict_severe_mild_deg = get_dictionary_diffreg(list_severe_mild_deg)
    list_severe_only_deg = read_deg_file(severe_only)
    dict_severe_only_deg = get_dictionary_diffreg(list_severe_only_deg)

    set_severe_mild_downreg = set(dict_severe_mild_deg["downreg"])
    set_severe_mild_upreg = set(dict_severe_mild_deg["upreg"])
    set_severe_only_downreg = set(dict_severe_only_deg["downreg"])
    set_severe_only_upreg = set(dict_severe_only_deg["upreg"])

    set_severe_overlap_downreg = set_severe_mild_downreg.intersection(set_severe_only_downreg) 
    set_severe_overlap_upreg = set_severe_mild_upreg.intersection(set_severe_only_upreg) 

    print(set_severe_overlap_upreg)
    print(set_severe_overlap_downreg)

    # list_mild_deg = read_deg_file(mild_only)
    # dict_mild_diffreg = get_dictionary_diffreg(list_mild_deg)
    # list_severe_deg = read_deg_file(severe_only)
    # dict_severe_diffreg = get_dictionary_diffreg(list_severe_deg)

    # set_severe_downreg = set(dict_severe_diffreg["downreg"])
    # set_severe_upreg = set(dict_severe_diffreg["upreg"])
    # set_mild_downreg = set(dict_mild_diffreg["downreg"])
    # set_mild_upreg = set(dict_mild_diffreg["upreg"])

    # set_severe_specific = set_severe_downreg.difference(set_mild_downreg) 
    # set_severe_specific = set_severe_upreg.difference(set_mild_upreg) 

def read_deg_file(filepath):
    list_deg = list()
    with open(filepath, mode='r') as fr:
        for line in fr:
            record = line.rstrip("\n").split("\t")
            list_deg.append(record)
    
    return list_deg
            
def get_dictionary_diffreg(list_deg):
    dict_diffreg = dict()
    for ind, deg_info in enumerate(list_deg):
        if ind == 0 :
            continue 

        deg_genes = deg_info[-1][1:-1].split(",")
        deg_genes = list(map(lambda x: x.lstrip(" ")[1:-1], deg_genes))

        if dict_diffreg.get(deg_info[0]) == None:
            dict_diffreg[deg_info[0]] = list()
        dict_diffreg[deg_info[0]].extend(deg_genes)

    return dict_diffreg
# %%
if __name__ == "__main__":
    # mild_only = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/6_time_series/severity_30_005/Visit1_Mild__Visit4_Mild/deg_stats.txt"
    # severe_only = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/6_time_series/severity_30_005/Visit1_Severe__Visit4_Severe/deg_stats.txt"
    severe_vs_mild = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/6_time_series/severity_15_005/Visit1_Severe__Visit1_Mild/deg_stats.txt"
    severe_only = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/7_compare_severity/severity_15_005/Visit1_Severe__Visit4_Severe/deg_stats.txt" 

    main(severe_vs_mild, severe_only)