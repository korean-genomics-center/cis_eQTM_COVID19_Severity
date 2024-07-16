{   
    library(ggplot2)
    library(ggrepel)
    library(readr)
    library(stringr)
}
{
    pcaData <- readr::read_tsv("/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/8_pca/pca_meta_dataframe.tsv")
}
{
    percentVar1 <- round(100*attr(pcaData, "percentVar"))
}
# {
#     draw_PC_biplot_with_lineplots <- function(pcaData, meta, outdir){
#         options(ggrepel.max.overlaps = Inf)
#         g <- ggplot(pcaData, aes(x=PC1, y=PC2, color=meta) )+
#             geom_point(size=1.5) +
#             # geom_path(aes(group=pcaData$ID),
#             #           size = 0.1, 
#             #           color="black",
#             #           arrow=arrow(type="open",
#             #                       length=unit(0.05, "inches"))) +
#             geom_text_repel(label=pcaData$ID, size=1.5, fontface="bold") + 
#             xlab(paste0("PC1 (",as.character(pcaData$PC1), "% variance)")) +
#             ylab(paste0("PC2 (",as.character(pcaData$PC2), "% variance)")) +
#             stat_ellipse(type="norm", linetype=2) +
#             scale_color_brewer(palette="Set2")+
#             # scale_color_grey() +
#             coord_cartesian(clip="off") +
#             ggtitle(title) +
#             theme_bw(base_size=10) +
#             theme(
#                 plot.title = element_text(size=10, face="bold", hjust = 0.5),
#                 axis.ticks = element_line(size=0.5),
#                 axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
#                 axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
#                 legend.position = "right"
#             )
#         ggsave(g, filename=paste0("PCA_",meta,".png"), device="png", path=outdir, width=8, height=8, dpi=300)
#     }
# }
{
    
}
{
    draw_pca_plot <- function(result_pr, table, args, PC, col_mark, draw_label=TRUE, draw_cluster=TRUE ){
        if(is.null(args$PC)){
            PC <- c(1, 2)
        }
        else{
            PC <- args$PC
        }
        qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
        color_palette <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
        unique_mark <- unique(table[[args$col_mark]])
        count_treatment <- length(unique_mark)
        color_palette_use <- color_palette[1:count_treatment]
        colordict <- make_color_dict(color_palette_use, unique_mark)
        for(ind_pc1 in 1:(length(PC)-1)){
            for(ind_pc2 in (ind_pc1+1):length(PC)){
                PC1 <- PC[ind_pc1]
                PC2 <- PC[ind_pc2]
                label <- args$draw_label
                cluster <- args$draw_cluster
                list_total_mean_se <- get_mean_se_of_each_mark(result_pr, table, args, PC1, PC2)
                output_prefix_split <- strsplit(args$output_prefix, split="/")[[1]]
                output_dir <- paste(c(output_prefix_split[1:length(output_prefix_split)-1]), collapse = '/')
                filename_prefix <- output_prefix_split[[length(output_prefix_split)]]
                filename <- paste(filename_prefix, "pca_plot", "PC", PC1, "PC", PC2, "png", sep = '.')
                fig <- autoplot(result_pr, data = table, x=PC1, y=PC2, scale=0, colour = args$col_mark, label = label, label.repel = label, frame = cluster, frame.type = "norm")+
                    ggtitle(filename)+
                    scale_color_manual(values = colordict)+
                    scale_fill_manual(values =  colordict)
                if(args$draw_center){
                    total_mark <- table[[args$col_mark]]
                    uniq_mark <- unique(total_mark)
                    for(ind_mark in 1:length(uniq_mark)){
                        mean_se_of_mark <- list_total_mean_se[[ind_mark]]
                        mark <- uniq_mark[ind_mark]
                        color_mark <- colordict[mark]
                        cross_v_x <- mean_se_of_mark$PC1$mean
                        cross_h_y <- mean_se_of_mark$PC2$mean
                        cross_v_y1 <- cross_h_y - mean_se_of_mark$PC2$se
                        cross_v_y2 <- cross_h_y + mean_se_of_mark$PC2$se
                        cross_h_x1 <- cross_v_x - mean_se_of_mark$PC1$se
                        cross_h_x2 <- cross_v_x + mean_se_of_mark$PC1$se
                        fig <- fig + geom_segment(aes_string(x=cross_v_x, y=cross_v_y1, xend=cross_v_x, yend=cross_v_y2),
                                            color=color_mark, size=1)+
                            geom_segment(aes_string(x=cross_h_x1, y=cross_h_y, xend=cross_h_x2, yend=cross_h_y),
                                            color=color_mark, size=1)
                    }
                }
                path_save <- paste(output_dir, filename, sep="/")
                ggsave(file=path_save, plot=fig)
            }
        }
    }
}
    draw_scree_plot <- function(result_pr, args){
        output_prefix_split <- strsplit(args$output_prefix, split="/")[[1]]
        output_dir <- paste(c(output_prefix_split[1:length(output_prefix_split)-1]), collapse = '/')
        filename_prefix <- output_prefix_split[[length(output_prefix_split)]]
        filename <- paste(filename_prefix, "scree_plot", "png", sep = '.')

        var_explained <- result_pr$sdev^2*100/sum(result_pr$sdev^2)
        var_df <- data.frame(PC = 1:length(var_explained), var_exp = var_explained)
        var_df$cumsum_var_exp <- cumsum(var_df$var_exp)

        if(! is.null(args$num_screeplot)){
            if(args$num_screeplot < length(var_explained)){
                var_df <- var_df[1:args$num_screeplot,]
            }
        }
        
        ggplot(var_df, aes(x=PC))+
            geom_col(aes(y=var_exp, fill="Proportion"))+
            scale_fill_manual(values="#777777")+
            geom_line(aes(y=cumsum_var_exp, color="Cumulative"))+
            scale_color_manual(values="black")+
            geom_point(aes(y=cumsum_var_exp))+
            scale_x_continuous(breaks=1:nrow(var_df))+
            theme(legend.key=element_blank(), legend.title=element_blank(), legend.box="horizontal",legend.position = "bottom")+
            ylim(0,100)+
            xlab("Principal Component")+
            ylab("Percentage Explained Variance [%]")+
            ggtitle(filename)
        path_save <- paste(output_dir, filename, sep="/")
        ggsave(file=path_save)
}
