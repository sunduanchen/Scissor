#' Preprocess the single-cell raw data using functions in the \code{Seurat} package
#'
#' This function provide a simplified-version of Seurat analysis pipeline for single-cell RNA-seq data. It contains the following steps in the pipeline:
#' \itemize{
#'    \item Create a \code{Seurat} object from raw data.
#'    \item Normalize the count data present in a given assay.
#'    \item Identify the variable features.
#'    \item Scales and centers features in the dataset.
#'    \item Run a PCA dimensionality reduction.
#'    \item Constructs a Shared Nearest Neighbor (SNN) Graph for a given dataset.
#'    \item Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm.
#'    \item Run t-distributed Stochastic Neighbor Embedding (t-SNE) dimensionality reduction on selected features.
#'    \item Runs the Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique.
#' }
#'
#' @param counts A \code{matrix}-like object with unnormalized data with cells as columns and features as rows.
#' @param project Project name for the \code{Seurat} object.
#' @param min.cells Include features detected in at least this many cells. Will subset the counts matrix as well.
#' To reintroduce excluded features, create a new object with a lower cutoff.
#' @param min.features Include cells where at least this many features are detected.
#' @param normalization.method Method for normalization.
#'   \itemize{
#'   \item LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor.
#'   This is then natural-log transformed using log1p.
#'   \item CLR: Applies a centered log ratio transformation.
#'   \item RC: Relative counts. Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor.
#'   No log-transformation is applied. For counts per million (CPM) set \code{scale.factor = 1e6}.
#' }
#' @param scale.factor Sets the scale factor for cell-level normalization.
#' @param selection.method How to choose top variable features. Choose one of :
#'   \itemize{
#'   \item vst: First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess).
#'   Then standardizes the feature values using the observed mean and expected variance (given by the fitted line).
#'   Feature variance is then calculated on the standardized values after clipping to a maximum (see clip.max parameter).
#'   \item mean.var.plot (mvp): First, uses a function to calculate average expression (mean.function) and dispersion (dispersion.function)
#'   for each feature. Next, divides features into num.bin (deafult 20) bins based on their average expression, and calculates
#'   z-scores for dispersion within each bin. The purpose of this is to identify variable features while controlling for the strong
#'   relationship between variability and average expression.
#'   \item dispersion (disp): selects the genes with the highest dispersion values
#'   }
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' @param dims_Neighbors Dimensions of reduction to use as input.
#' @param dims_TSNE Which dimensions to use as input features for t-SNE.
#' @param dims_UMAP Which dimensions to use as input features for UMAP.
#' @param verbose Print output.
#'
#' @return A \code{Seurat} object containing t-SNE and UMAP representations.
#' @import Seurat
#' @export
Seurat_preprocessing <- function(counts, project = "Scissor_Single_Cell", min.cells = 400, min.features = 0,
                                 normalization.method = "LogNormalize", scale.factor = 10000,
                                 selection.method = "vst", resolution = 0.6,
                                 dims_Neighbors = 1:10, dims_TSNE = 1:10, dims_UMAP = 1:10,
                                 verbose = TRUE){
    library(Seurat)
    data <- CreateSeuratObject(counts = counts, project = project, min.cells = min.cells, min.features = min.features)
    data <- NormalizeData(object = data, normalization.method = normalization.method, scale.factor = scale.factor, verbose = verbose)
    data <- FindVariableFeatures(object = data, selection.method = selection.method, verbose = verbose)
    data <- ScaleData(object = data, verbose = verbose)
    data <- RunPCA(object = data, features = VariableFeatures(data), verbose = verbose)
    data <- FindNeighbors(object = data, dims = dims_Neighbors, verbose = verbose)
    data <- FindClusters( object = data, resolution = resolution, verbose = verbose)
    data <- RunTSNE(object = data, dims = dims_TSNE)
    data <- RunUMAP(object = data, dims = dims_UMAP, verbose = verbose)

    return(data)
}

