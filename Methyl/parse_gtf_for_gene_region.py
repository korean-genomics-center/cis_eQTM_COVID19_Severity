#%%
import pandas as pd
import re


def main(path_gtf, path_save, cols = ["chr", "src", "feature", "start", "end", "score", "strand", "frame", "attrib"], cols_save = ["chr", "start", "end", "strand"], attrib_save = ["gene_id", "gene_name", "gene_type"]):
    gtf_gene_line_generator = get_gene_line_from_gtf_generator(path_gtf, cols)
    cols_save_tot = cols_save + attrib_save
    with open(path_save, 'w') as fw:
        fw.write('\t'.join(cols_save_tot) + '\n')
        for line in gtf_gene_line_generator:
            list_values_line = line.strip('\n').split('\t')
            dict_values_line = dict(zip(cols, list_values_line))
            dict_attrib = parse_gtf_attribute(dict_values_line["attrib"])
            dict_values_line.update(dict_attrib)
            list_values_write = list(map(dict_values_line.__getitem__, cols_save_tot))
            fw.write('\t'.join(list_values_write) + '\n')

def get_gene_line_from_gtf_generator(path_gtf, cols):
    where_feat = cols.index("feature")
    with open(path_gtf, 'r') as fr:
        for line in fr:
            list_values = line.strip('\n').split('\t')
            val_type = list_values[min(len(list_values)-1, where_feat)]
            if val_type == "gene":
                yield line

def parse_gtf_attribute(attrib):
    dict_attrib = dict()
    re_attrib_parsing = re.compile('(.*) "?(.*)"?')
    for attrib_pair in attrib.split(';'):
        attrib_pair = attrib_pair.strip(' ')
        if len(attrib_pair) == 0: continue
        key, val = re_attrib_parsing.match(attrib_pair).groups()
        dict_attrib[key] = val.strip('"')
    return dict_attrib
    

#%%
if __name__ == "__main__":
    path_gtf = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/References/Genome/Gencode/v42/gencode.v42.chr_patch_hapl_scaff.annotation.gtf"
    
    path_save = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/GTF_Gene_Parsing/GTF.genecode.v42.chr_patch_hapl_scaff.gene_only.tsv"
    
    main(path_gtf, path_save)
# %%
