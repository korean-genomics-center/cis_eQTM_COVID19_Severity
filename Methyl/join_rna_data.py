#%%
import pandas as pd
from itertools import chain


def main(path_meta, col_sample, path_rsem_format, path_deseq_norm, path_save):
    table_meta = pd.read_csv(path_meta, sep = '\t')
    list_samples = table_meta[col_sample].to_list()
    
    table_deseq_norm = pd.read_csv(path_deseq_norm, sep = '\t')
    list_expression_dictionary = list(map(lambda sample : get_expression_values_for_single_sample(path_rsem_format, table_deseq_norm, sample), list_samples))
    list_expression_table_per_sample = list(map(lambda sample, dict_exp : organize_expression_dictionaries_as_dataframe(sample, dict_exp), list_samples, list_expression_dictionary))
    table_tot_expression = pd.concat(list_expression_table_per_sample)
    table_tot_expression.to_csv(path_save, sep = '\t', index = False)
    
    
def get_expression_values_for_single_sample(path_rsem_format, table_deseq_norm, sampleid):
    path_rsem = path_rsem_format.format(sampleid)
    dict_exp = read_rsem_results(path_rsem)
    dict_exp["DESeq_Norm"] = table_deseq_norm[sampleid].to_dict()
    return dict_exp
    
def read_rsem_results(path_rsem, col_gene = "gene_id", col_exps = ["expected_count", "TPM"]):
    table_rsem = pd.read_csv(path_rsem, sep = '\t')
    dict_exp = dict()
    for col in col_exps:
        dict_exp[col] = dict(zip(table_rsem[col_gene], table_rsem[col]))
    return dict_exp

def organize_expression_dictionaries_as_dataframe(sampleid, dict_expression):
    cols_exp = list(dict_expression.keys())
    table_exp = pd.DataFrame(columns = ["Sample_ID", "Gene_ID"] + cols_exp)
    set_geneid = set(list(chain(*list(map(lambda col_exp : dict_expression[col_exp].keys(), cols_exp)))))
    geneid_sortedname = sorted(set_geneid)
    table_exp["Gene_ID"] = geneid_sortedname
    table_exp["Sample_ID"] = sampleid
    for col_exp in cols_exp:
        table_exp[col_exp] = table_exp["Gene_ID"].apply(lambda gene : dict_expression[col_exp].get(gene, pd.NA))
    print(sampleid)
    
    return table_exp

#%%
if __name__ == "__main__":
    path_meta = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Resources/MetaTable/COVID19_master_table_20231007.Methyl_Overlap.with_Severity.20240402.txt"
    col_sample = "Project_ID_Alias"
    path_rsem_format = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/3_rsem/ConfirmedRecovered/{0}/{0}.genes.results"
    path_deseq_norm = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/DEGExtract/DEGExtract.RNA_samples_with_Methyl.DEG_by_Sex.Control_M.Case_F.Cov_Age.mincount_1.20240402.tsv.normcount"
    
    path_save = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/RNAExpressionTable/RNAExpression.COVID19.RNA_samples_with_Methyl.tsv"
    
    main(path_meta, col_sample, path_rsem_format, path_deseq_norm, path_save)
#%%