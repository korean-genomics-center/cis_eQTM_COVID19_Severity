library(biomaRt)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

transcript_ids <- c(
'ENST00000492017.3',
'ENST00000565570.1',
'ENST00000590595.1',
'ENST00000623521.1',
'ENST00000664895.1',
'ENST00000513652.1',
'ENST00000518675.1'

)

res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'ensembl_gene_id', 
                            'external_transcript_name',
                            'external_gene_name'),
             filters = 'ensembl_transcript_id_version', 
             values = transcript_ids,
             mart = mart)

write.table(res, "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP/disovery_markers/marker_231211/severe_mild_firstVisit/annotation/methyl_hypo_severe_mild_firstVisit_annotatr_neighbor_methylation_enst_conv.tsv", sep = "\t", quote = FALSE, row.names=FALSE)