#%%
import pandas as pd


def main(path_meta, col_sample, path_rsem_format, col_rsem_genename, path_filter, col_filter_genename, path_save_format):
    table_meta = pd.read_csv(path_meta, sep = '\t')
    list_sample = table_meta[col_sample].to_list()
    
    table_filter = pd.read_csv(path_filter, sep = '\t')
    list_genes_filter = table_filter[col_filter_genename].to_list()
    
    for sample in list_sample:
        path_rsem = path_rsem_format.format(sample)
        table_rsem = pd.read_csv(path_rsem, sep = '\t')
        table_rsem_filtered = table_rsem[table_rsem[col_rsem_genename].isin(list_genes_filter)]
        table_rsem_filtered.to_csv(path_save_format.format(sample), sep = '\t', index = False)



#%%
if __name__ == "__main__":
    path_meta = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Resources/MetaTable/COVID19_master_table_20231007.Methyl_Overlap.with_Severity.20240402.txt"
    col_sample = "Project_ID_Alias"
    path_rsem_format = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/3_rsem/ConfirmedRecovered/{0}/{0}.genes.results"
    col_rsem_genename = "gene_id"
    
    path_filter = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/Methyl_RNA_Relation/Methyl_RNA_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Mild_Sev_Visit1.LOO_common.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_6.DMP_DEG_LOO_annot.20240402.tsv"
    col_filter_genename = "RNA"
    
    path_save_format = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/RSEM_Filtered/Genes_Mild_Sev_Visit1.18_genes/{0}.genes.results.filtered"
    
    main(path_meta, col_sample, path_rsem_format, col_rsem_genename, path_filter, col_filter_genename, path_save_format)