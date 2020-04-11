#' Scissor: Single-Cell Identification of Subsets with bulk RNA-Seq phenotype cORrelation
#'
#' Scissor is a novel proposed network-based model, which can identify the cooperative
#' biomarkers for heterogeneous complex disease diagnoses.
#'
#' Scissor is a novel proposed network-based model, which can identify the cooperative
#' biomarkers for heterogeneous complex disease diagnoses. Scissor uses the Seurat package for the
#' preprocessings of single-cell RNA-seq data.
#'
#' @param alpha The parameter balances the effect of l1 norm and network-based penalties. The range of alpha is in \code{[0,1]}.
#' A larger alpha will lay more emphasis on l1 norm, i.e. Scissor will select fewer cells with more sparsity.
#' @param family Type of outcome. Can be binomial for discrete phenotype or cox for clinical outcome.
#' @param Given_Inputs Whether a saved regression input is given for tuning parameter. See tutorial for more details.
#' @param bulk_dataset The bulk RNA-seq expression data for Scissor. The rows are genes and the columns are samples.
#' @param sc_dataset The single-cell RNA-seq expression data for Scissor. The rows are genes and the columns are samples.
#' @param survival_info A three-column format of clinical data used for cox regression. The first column contains the
#' sample ID (observations). The second column is the survival time and the third column stands for the event
#' (e.g.recurrence of cancer or death), which is indicated by the variable "Status" with 0 (no event, censored) and 1 (event).
#' @param phenotype The discrete vector for logistic regression.
#' @param cor_type The correlation method in constructing the correlation matrix between bulk samples and single cells.
#' @param counts Whether the input single-cell data is count quantification. If yes, an additional normalization step is
#' applied in preprocessing steps.
#' @param save_file The file name for saving. Scissor can save the preprocessed regression inputs into a R data, which
#' can help users to directly tune model parameters. See tutorial for more details.
#'
#' @return This function will return a list with the following components:
#'   \item{Coefs}{The regression coefficient for each cell.}
#'   \item{Cell_positive}{The cell IDs that have positive coefficient.}
#'   \item{Cell_negative}{The cell IDs that have negative coefficient.}
#'   \item{Seurat_data}{The preprocessed Seurat object that contains the single cell information.}
#'
#' @references Duanchen Sun and Zheng Xia (2020): Scissor: Phenotype guided single-cell
#' subset identification from bulk RNA-seq.
#'
#' @import Seurat Matrix preprocessCore
#'
#' @export
Scissor <- function(alpha, family = c('binomial', 'cox'), Given_Inputs = NULL, bulk_dataset, sc_dataset,
                    survival_info = NULL, phenotype = NULL, cor_type = c('pearson', 'spearman'), counts = TRUE,
                    save_file = 'Regression_Inputs.RData'){
    library(Seurat)
    library(Matrix)
    library(preprocessCore)

    if (is.null(Given_Inputs)){
        data <- CreateSeuratObject(counts = sc_dataset, project = 'Scissor_Single_Cell', min.cells = 400, min.features = 0)
        if(counts) data <- NormalizeData(object = data, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = F)
        data <- FindVariableFeatures(object = data, selection.method = 'vst')
        data <- ScaleData(object = data, verbose = F)
        data <- RunPCA(object = data, features = VariableFeatures(data), verbose = F)
        data <- FindNeighbors(object = data, dims = 1:10, verbose = F)
        data <- FindClusters(object = data, resolution = 0.6, verbose = F)
        data <- RunTSNE(object = data, dims = 1:10)
        data <- RunUMAP(object = data, dims = 1:10, verbose = F) # umap.method = 'umap-learn', metrix = 'correlation'

        network <- as.matrix(data@graphs$RNA_snn)
        diag(network) <- 0
        network[which(network != 0)] <- 1

        common <- intersect(rownames(bulk_dataset), rownames(data))
        if (length(common) == 0) {
            stop("There is no common genes between bulk and single-cell RNA-seq datasets.")
        }
        dataset0 <- cbind(bulk_dataset[common,], as.matrix(data@assays$RNA@data[common,]))  # Dataset before quantile normalization.
        dataset1 <- normalize.quantiles(dataset0)                                           # Dataset after  quantile normalization.
        rownames(dataset1) <- rownames(dataset0)
        colnames(dataset1) <- colnames(dataset0)

        Expression_bulk <- dataset1[,1:ncol(bulk_dataset)]
        Expression_pbmc <- dataset1[,(ncol(bulk_dataset) + 1):ncol(dataset1)]
        X <- cor(Expression_bulk, Expression_pbmc, method = cor_type)

        if (is.null(phenotype)){
            Y <- cbind(time = survival_info[,2], status = survival_info[,3])
            save(X, Y, data, network, survival_info, file = save_file)
        }else{
            Y <- phenotype
            save(X, Y, data, network, file = save_file)
        }
    }else{
        load(Given_Inputs)
    }

    set.seed(123)
    fit0 <- APML1(X, Y, family = family, penalty = 'Net', alpha = alpha, Omega = network, nlambda = 100, nfolds = 10)
    fit1 <- APML1(X, Y, family = family, penalty = 'Net', alpha = alpha, Omega = network, lambda = fit0$lambda.min)
    Coefs <- as.numeric(fit1$Beta)

    Cell1 <- colnames(X)[which(Coefs > 0)]
    Cell2 <- colnames(X)[which(Coefs < 0)]
    scissor_ident <- rep(0, ncol(data))
    names(scissor_ident) <- colnames(data)
    scissor_ident[Cell1] <- 1
    scissor_ident[Cell2] <- 2
    data <- AddMetaData(object = data, metadata = scissor_ident, col.name = "scissor")
    print(sprintf('Scissor output information: Positive cells: %d. Negative cells: %d.', length(Cell1), length(Cell2)))

    return(list(Coefs = Coefs,
                Cell_positive = Cell1,
                Cell_negative = Cell2,
                Seurat_data = data))
}

