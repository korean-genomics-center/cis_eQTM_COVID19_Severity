# code to test the leave one out validation 
# for detail, check
# http://in.genomelab.org/Stressomics:Validataion:Leave_One_out
# The code is written in R
# #%% indicate code block. 
# creating conda environment
#   conda install -c conda-forge r-base=4.1
#   conda install -c conda-forge r-essentials
#   conda install -c conda-forge r-xml
{
    #Install required packages and load packages
    print("Load Required Packages")
    installed_packages = installed.packages()[,"Package"]
    required_packages = c("BiocManager", "readr", "stringr", "data.table", "argparse", "ggfortify", "ggrepel", "ggplot2", "RColorBrewer")
    Bioc_required_packages = c("methylKit", "BiocParallel")
    for (req_package in required_packages){
        # if(! req_package %in% installed_packages){
        #     print(paste0("Installing ", req_package))
        #     suppressMessages(install.packages(req_package, repos="https://cran.biodisk.org/", quiet = TRUE))
        # }
        suppressMessages(library(req_package, character.only=TRUE))
        print(paste0("Load ", req_package))
    }

    #Update packages
    # print("Updating Packages")
    # suppressMessages(update.packages(ask=FALSE, checkBuilt=TRUE, repos="https://cran.biodisk.org/"))
    # suppressMessages(BiocManager::install(version = '3.14', ask = FALSE))

    print("Load BiocManager Packages")
    for (req_bioc_package in Bioc_required_packages){
        # if(! req_bioc_package %in% installed_packages){
        #     print(paste0("Installing ", req_bioc_package))
        #     BiocManager::install(req_bioc_package, version = '3.14', ask = FALSE, update = TRUE)
        # }
        suppressMessages(library(req_bioc_package, character.only=TRUE))
        print(paste0("Load ", req_bioc_package))
    }
}

### parse arguments
{
    parser <- ArgumentParser()
    parser$add_argument("--tsv", help = "Path to tsv file contains id and file path")
    parser$add_argument("--RDS", help = "Path to RDS compressed MethylBase file (Pre-united methyl CpG Table)")
    parser$add_argument("--col_id", help = "Column to sample ID")
    parser$add_argument("--col_path", help = "Column to methylation file path")
    parser$add_argument("--cov_low_count_cutoff", help="Cutoff of Lower Count Coverage [default : 10]", required = FALSE, type = 'integer', default = 10)
    parser$add_argument("--output", help="Output path of result")
    parser$add_argument("--threads", help="The number of threads", default=1, type='integer')
    args <- parser$parse_args();
}

### set arguments value
{
    numcore <- as.integer(args$threads);
    register(MulticoreParam(numcore));
}

read_tsv_custom <- function(path){
    table <- read_tsv(path, show_col_types = FALSE);
    return(table);
}
read_methyl_CpG_tables <- function(total_paths, total_ids, treatment, args){
    #%%
    # Read in methylation files required. 
    # myobj has the full information used below.  
    # The data are initially filtered out by the coverage
    # and united into one file.
    methyl_obj <- methRead(as.list(total_paths), sample.id =as.list(total_ids), assembly = "hg38",treatment = treatment, context = "CpG", pipeline=list(fraction=FALSE, chr.col=1, start.col=2, end.col=2, coverage.col=4, strand.col=3, freqC.col=5))
    methyl_obj_coverage_filtered <- filterByCoverage(methyl_obj, lo.count = args$cov_low_count_cutoff, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)
    methyl_obj <- NULL
    methyl_united <- unite(methyl_obj_coverage_filtered, mc.cores = args$threads)
    methyl_obj_coverage_filtered <- NULL
    return(methyl_united);
}
read_methyl_DB_table <- function(args, total_ids, treatment){
    print("Reading RDS formatted methylation table...")
    methyl_united <- readRDS(args$RDS)
    return(methyl_united);
}
{
    if(is.null(args$RDS)){
        config_table <- read_tsv_custom(args$tsv)
        total_ids <- c(config_table[[args$col_id]])
        total_paths <- c(config_table[[args$col_path]])
        treatment <- c(rep(0, nrow(config_table)))
        methyl_united <- read_methyl_CpG_tables(total_paths, total_ids, treatment, args)
    }
    else{
        methyl_united <- read_methyl_DB_table(args, total_ids, treatment)
    }
}
get_position_identifier <- function(chr, start, end){
    identifier <- paste(chr, start, end, sep = '_')
    return(identifier);
}
{
    mat <- getData(methyl_united)
    # mat$position <- do.call(function(chr, start, end, ...) get_position_identifier(chr, start, end), mat)
    # position <- mat[, "position"]
    position <- mat[, c("chr", "start", "end")]
    methyl_met <-  mat[, methyl_united@numCs.index]/
      (mat[,methyl_united@numCs.index] + mat[,methyl_united@numTs.index] ) * 100
    names(methyl_met) <- methyl_united@sample.ids
    methyl_met <- cbind(position, methyl_met)
    write.table(methyl_met, args$output, row.names=FALSE, sep="\t", quote=FALSE)
}