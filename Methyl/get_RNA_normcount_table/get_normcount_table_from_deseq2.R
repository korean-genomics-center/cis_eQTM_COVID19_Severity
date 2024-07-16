{   #parse arguments
    library(argparse)
    parse_arguments <- function(){
        parser <- ArgumentParser()
        parser$add_argument("--subject_info_table_1", required=TRUE)
        parser$add_argument("--subject_info_table_2", required=TRUE)
        parser$add_argument("--id_column", required=TRUE)
        parser$add_argument("--filein", required=TRUE)
        parser$add_argument("--fileout", required=TRUE)
        parser$add_argument("--design", required=TRUE)
        parser$add_argument("--compare_var", required=TRUE)
        parser$add_argument("--meta_columns", required=TRUE)
        parser$add_argument("--control_group", required=TRUE)
        parser$add_argument("--case_group", required=TRUE)
        parser$add_argument("--mincount", required=TRUE) 
        parser$add_argument("--num_cores", required=TRUE)
        parser$add_argument("--sep", required=TRUE)
        inputargs <- parser$parse_args()
        return (inputargs)
    }

}

{   # get input data
    library("stringr")
    library("tximport")
    library("jsonlite")

    inputargs <- parse_arguments()

    # get sample metadata
    samplemetadf1 <- read.csv(inputargs$subject_info_table_1, sep="\t", colClasses="character")
    samplemetadf2 <- read.csv(inputargs$subject_info_table_2, sep="\t", colClasses="character")
    all_samples <- rbind(samplemetadf1, samplemetadf2)

    # get gene expression matrix
    template <- inputargs$filein

    filepathfun <- function(samplename){
        return(str_replace_all(template, "\\[id_column\\]",samplename))
    }

    files_in <- unlist(lapply(all_samples[[inputargs$id_column]], filepathfun))
    txi_rsem <- tximport(files_in, type="rsem", txIn=FALSE, txOut=FALSE)
    counts <- as.data.frame(txi_rsem$counts)
    colnames(counts) <- all_samples[[inputargs$id_column]]
    metainfo <- all_samples[, unlist(str_split(inputargs$meta_columns, pattern=inputargs$sep))]
}

{   # make DESeq2 object
    library(DESeq2)
    dds <- DESeqDataSetFromMatrix(countData=round(counts), colData=metainfo, design=as.formula(inputargs$design))
}

{   # filter out by mincounts
    if("mincount" %in% names(inputargs)){
        keep <- rowSums(counts(dds)) >= inputargs$mincount
        ddsfiltered <- dds[keep,]
    }
}

{   # relevel to clearly define case and control groups
    condition <- inputargs$compare_var
    ddsrelevel <- ddsfiltered
    ddsrelevel[[condition]] <- relevel(factor(ddsrelevel[[condition]]), ref=inputargs$control_group)
}

{   # run DESeq2 normalization and test
        # estimate size factor by all genes (by default) to set minimum expression as 0 plus alpha - i.e. 샘플별 상대적 전체 생산량 (기본 발현값이 0이 아니라 그 이상) 
        # This function estimates the size factors using the "median ratio method" described by Equation 5 in Anders and Huber (2010).
    ddssize <- estimateSizeFactors(ddsrelevel)
        # estimate dispersion or alpha (a form of measuring spread of data; 0.01 dispersion = 10% variation around expected mean; var = mu + alpha*mu^sq;) 
        # dispersion depends on size factor and covariate
    ddsdisp <- estimateDispersions(ddssize)
        # DEseq2 model uses a negative binomial generalized linear models 
        # K[i,j] ~ NegBin(mu[i,j], alpha[i])
        # Counts (K) for gene i, sample j are modelled using a negative binomial distribution with fitted mean mu[i,j] and a gene-specific dispersion parameter alpha[i]
        # Hypothesis test for DEG: Wald test - we use the estimated SE of a log2FC to test if it is equal to zero 
    ddsnorm <- nbinomWaldTest(ddsdisp, maxit=1000)
}

{   # save normalized count data (i.e. normalized = raw_count/sizefactor) N.B. dds object only has the raw count in dds$count throughout the analysis thus far 
    normcount <- counts(ddsnorm, normalized=TRUE)
    write.table(as.data.frame(normcount), file=paste0(inputargs$fileout,".normcount"), row.names=TRUE, sep="\t", quote=FALSE)
}

{   # save vst (variance-stabilization transformation) object used as an input for downstream visualization during EDA (QC)
    vsd <- vst(ddsnorm, blind=FALSE, nsub = nrow(ddsnorm))
    saveRDS(vsd, paste0(inputargs$fileout,".vstobj"))
}