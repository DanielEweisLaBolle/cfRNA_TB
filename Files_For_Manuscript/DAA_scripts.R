
##------------------------------------
# Make heatmap
##------------------------------------

make_heatmap <- function(res,
                        counts,
                        samples,
                        SAMPLE_ID_VAR,
                        MYCOLORS,
                        VARS,
                        SIG_THRESH,
                        FC_THRESH
                        ){


    sig_genes <- data.frame(res) %>% 
        filter(padj < SIG_THRESH & abs(log2FoldChange) > FC_THRESH) %>% 
        # filter(gene_type == "protein_coding") %>% 
        rownames()

    length(sig_genes)

    ##---------------------
    # Prepare count matrix and normalize
    counts <- DESeq2::varianceStabilizingTransformation(as.matrix(counts)) %>% data.frame() 
    mat <- counts[sig_genes,]

    ##---------------------
    # Prepare meta data and colors for annotations
    anno = samples[match(colnames(mat),samples[[SAMPLE_ID_VAR]]),] %>% 
        dplyr::select(all_of(VARS))
    rownames(anno) = colnames(mat)

    mycolors <- MYCOLORS

    color = colorRampPalette(c("blue","yellow"))(50)


    breaksList = seq(-2, 2, length = 51)

    mat <- t(scale(t(mat)))

    ##---------------------
    # Plot

    heatmap_plt <- pheatmap(mat,
                          # Colors
                 col=color,
                 breaks=breaksList, 
                 annotation_col=anno,
                 annotation_colors=mycolors,
                #  na_col = "#FFFFFF",

                 # Fonts
                 show_colnames=F,
                 show_rownames=F,
                 fontsize=12,
                 fontsize_col=3,
                 annotation_names_col=F,
                 annotation_names_row=F,

                 # Clustering

                 clustering_distance_cols="correlation",
                 clustering_distance_rows="correlation", 
                 treeheight_row=0,
                treeheight_col= 15,

                 # Misc.
                 border_color=NA,
                legend=FALSE,
                annotation_legend=FALSE)
    
    return(heatmap_plt)
    }






########################################################
########################################################
########################################################
#### OLD 
########################################################
########################################################

##------------------------------------
# Load metadata and counts
##------------------------------------
load_data <- function(meta_data_file,
                      count_matrix_file,
                      SAMPLE_ID_VAR,
                      COMP_VAR,
                      GROUPS,
                      gene_name_key,
                     subset=NA){


    ##------------------------------------
    # Load sample metadata file
    samples <- read.delim(meta_data_file) %>%
        filter(passQC == TRUE) %>% 
        filter(.data[[COMP_VAR]] %in% all_of(GROUPS)) %>%
        filter(!is.na(.data[[SAMPLE_ID_VAR]])) %>%
        mutate(raw_id = sample_id) %>%
        mutate(sample_id = ifelse(!grepl("cfrna",sample_id),paste0("X",sample_id),sample_id))%>%
        mutate(sample_id = gsub("-",".",sample_id))

    samples$expGroup <- samples[[COMP_VAR]]
    rownames(samples) <- samples$sample_id
    
    if (!is.na(subset)){samples <- samples %>% filter(raw_id %in% subset)}
    
    print(table(samples$expGroup))      
    # head(samples)

    ##------------------------------------
    # Load count data
    sample_ids = unique(samples[[SAMPLE_ID_VAR]])                                                                    # Get sample ids that pass qc
    counts = read.delim(count_matrix_file,row.names=1)
    # rownames(counts) <- counts$geneID
    
    # return(list("samples"=samples,"counts"=counts))
    
    counts = counts[,sample_ids]

    ##------------------------------------
    # Remove ChrX, ChrY, ChrM, and RB genes
    gene.list <- read.delim(gene_name_key,col.names = c("type,","ENSMBL","gene_symbol"))
    gene.ids <- gsub("\\..*","",rownames(counts))
    exclude.idx <- gene.ids %in% gene.list[,2]
    counts = counts[!exclude.idx,]     
    # head(counts)

    return(list("samples"=samples,"counts"=counts))
}

##------------------------------------
# RUN DESeq
##------------------------------------

run_DESeq <- function(counts,
                      samples,
                      DESIGN,
                      gene_name_key,
                      GROUPS){
    ##------------------------------------
    # Contstruct DESeq Data Set
    dds <- DESeqDataSetFromMatrix(counts,
                                    colData = samples,
                                    design = formula(DESIGN))


    ##------------------------------------
    # Add Gene metadata
    annotation = fread(file=gene_name_key)
    annotation <- annotation[match(rownames(dds), annotation$gene_id),]
    all(rownames(dds) == annotation$ftcount_id)
    mcols(dds) <- cbind(mcols(dds), annotation)


    ##------------------------------------
    # Re-factor
    dds$Category <- factor(dds$Category, levels = GROUPS)

    ##------------------------------------
    # Pre-filter
    # keep <- rowSums(counts(dds)) >= 10
    # dds <- dds[keep,]

    ##------------------------------------
    # DAA
    dds <- DESeq(dds)
    
    return(dds)
    }


make_complex_heatmap <- function(res,
                            counts,
                            samples,
                            SAMPLE_ID_VAR,
                            SIG_THRESH = 0.05,
                            FC_THRESH = 1){
    
    sig_genes <- data.frame(res) %>% 
        filter(padj < SIG_THRESH & abs(log2FoldChange) > FC_THRESH) %>% 
        filter(gene_type == "protein_coding") %>% 
        rownames()

    length(sig_genes)

    ##---------------------
    # Prepare count matrix and normalize
    cpm_counts <- data.frame(edgeR::cpm(counts))
    mat <- cpm_counts %>% filter(row.names(cpm_counts) %in% all_of(sig_genes))

    ##---------------------
    # Prepare meta data and colors for annotations
    anno =  samples[match(colnames(mat),samples[[SAMPLE_ID_VAR]]),] %>% 
        dplyr::select(Diagnosis,Category)
    
    rownames(anno) = colnames(mat)

    mycolors <- list(
        Diagnosis = c("KD" = "red4",
                     "MIS-C" = "blue",
                     "FC" = "green"),
        Category = c("ALT normal" = "#2CDB43",
                     "ALT high" = "#21A433",
                     "Cluster 1" = "#9C2A1F", 
                     "Cluster 4" = "#F44130",
                     'lowest EF<55' = "#1F6A9C",
                     'lowest 60>EF≥55' = "#2988C8",
                     'lowest EF≥60' = "#33ADff",
                     'lowest EF≥55' = "#FFFFFF")
    )
    
    column_ha = HeatmapAnnotation (
        # Cateogry = anno$Category,
        # Diagnosis = anno$Diagnosis,
        df =anno,
        col = mycolors,
        annotation_legend_param = list(
            Category = list(
                        title = "Category"),
            Diagnosis = list(
                        title = "Diagnosis",
                        nrow = 2)
                
        )
        )

    
    # col_fun = colorRamp2(c(-2,2), hcl_palette = "Reds", reverse = TRUE)
    # col_fun = colorRamp2(c(-2,2), hcl_palette = "Blue-Red 3")
    # col_fun = colorRamp2(c(-2,2), hcl_palette = "RdBu",reverse=TRUE)
    
    # col_fun = colorRamp2(c(-2,2), c("blue", "yellow"))
    col_fun = colorRamp2(c(-2,0,2), c("blue", "black", "yellow"))

    mat <- t(scale(t(mat)))
    
    lgd = Legend


    htmap <- Heatmap(mat,
           name = "Z-Score (CPM)",
           col = col_fun,
            
           show_row_dend = FALSE,
           column_dend_height = unit(2, "cm"),
            clustering_distance_columns = "pearson",
            
            show_row_names = FALSE,
            show_column_names = FALSE,
            
            top_annotation = column_ha,
            heatmap_legend_param = list(
                legend_direction = "vertical", 
                legend_width = unit(2, "cm")), 
                )
    
#     htmap <- draw(htmap, heatmap_legend_side="right", 
#             annotation_legend_side="bottom",
#            legend_grouping = "original")
    
    return(htmap)
    
    }