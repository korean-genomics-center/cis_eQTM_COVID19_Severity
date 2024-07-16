#%%
import pandas as pd


def main(path_rna, col_rna_geneid, path_filter, col_filter_geneid, path_save):
    table_rna = pd.read_csv(path_rna, sep = '\t')
    table_filter = pd.read_csv(path_filter, sep = '\t')
    
    list_genes_filter = table_filter[col_filter_geneid].to_list()
    table_rna_filtered = table_rna[table_rna[col_rna_geneid].isin(list_genes_filter)]
    
    table_rna_filtered.to_csv(path_save, sep = '\t', index = False)
    
    
if __name__ == "__main__":
    path_rna = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/RNAExpressionTable/RNAExpression.COVID19.RNA_samples_with_Methyl.tsv"
    col_rna_geneid = "Gene_ID"
    path_filter = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/Methyl_RNA_Relation/Methyl_RNA_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Mild_Sev_Visit1.LOO_common.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_6.DMP_DEG_LOO_annot.20240402.tsv"
    col_filter_geneid = "RNA"
    path_save = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/RNAExpressionTable/RNAExpression.COVID19.RNA_samples_with_Methyl.filter_Genes.Mild_Sev_Visit1.LOO_common.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_6.tsv"
    
    main(path_rna, col_rna_geneid, path_filter, col_filter_geneid, path_save)
# %%
