get_enrichment <- function(Seurat_data, selected_cells, background_cells = NULL, directory = 'Enrichment_Results', short = 'Group1_vs_Group2'){
    
    library(viper)
    library(limma)
    if (!dir.exists(directory)){
        dir.create(directory)
    }
    if (is.element('enrichment', colnames(Seurat_data@meta.data))){
        Seurat_data$enrichment <- NULL
    }
    enrichment_ident <- rep(0, ncol(Seurat_data))
    names(enrichment_ident) <- colnames(Seurat_data)
    enrichment_ident[selected_cells] <- 1
    Seurat_data <- AddMetaData(object = Seurat_data, metadata = enrichment_ident, col.name = 'enrichment')
    Seurat_data <- SetIdent(object = Seurat_data, value = 'enrichment')
    All_genes <- rownames(Seurat_data)
    non_ribosome_genes <- setdiff(All_genes, All_genes[grep('^RPL|^RPS', All_genes)])
    
    if (is.null(background_cells)){
        DEG_result <- FindMarkers(object = Seurat_data, ident.1 = 1, logfc.threshold = log(1), min.pct = 0, features = non_ribosome_genes)
    }else{
        enrichment_ident[background_cells] <- 2
        DEG_result <- FindMarkers(object = Seurat_data, ident.1 = 1, logfc.threshold = log(1), min.pct = 0, features = non_ribosome_genes, ident.1 = 2)
    }
    DEG_result$Gene <- rownames(DEG_result)
    DEG_result <- DEG_result[, c(6,2,3,4,1,5)]
    write.table(DEG_result, file = paste0(directory, '/', short, '_DEG_infos.txt'), sep = '\t', quote = F, row.names = F)
    
    rds_paths <- '/Users/sund/Desktop/scissor/data/offical_rds'
    signature <- DEG_result[,'avg_logFC'] * log10(1/(DEG_result[,'p_val_adj'] + 1e-300))
    gene_sets <- dir(rds_paths)
    names(signature) <- DEG_result$Gene
    for (i in 1:length(gene_sets)){
        hg_v6  <- readRDS(paste(rds_paths, gene_sets[i], sep = '/'))
        index  <- ids2indices(hg_v6, DEG_result$Gene)
        camera <- cameraPR(signature, index)
        if(!is.null(camera)){
            write.csv(camera, file = paste0(directory, '/', short, '_Camera_', gsub('.v6.2.symbols.gmt.rds', '', gene_sets[i]), '.csv'), quote = F)
        }
    }
    print('Enrichment analysis finished!')
    
    pdf(paste0(directory, '/MARINA_', short, '.pdf'))
    regulons <- c('regulon', 'kinase')
    for(choose in regulons){
        if (choose == 'regulon') load('/Users/sund/Desktop/scissor/data/regulon/post-processing-regulon-TF.RData')
        if (choose == 'kinase')  load('/Users/sund/Desktop/scissor/data/regulon/post-processing-regulon-kinase-level1.RData')
        mrs <- msviper(signature, regul, minsize = 5, verbose = FALSE)
                plot(mrs, cex = 1.2, 30)
        write.table(summary(mrs,length(mrs$es$nes)), file = paste0(directory, '/MARINA', short, '_.txt'), sep = '\t', quote = F, row.names = F)
    }
    dev.off()
    print('Regulon analysis finished!')
}
