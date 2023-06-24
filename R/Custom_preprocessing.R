#' Add custom training datasets to perform prediction
#'
#' @description \code{Custom_preprocessing} Preprocessing for custom training datasets.
#'
#' @param counts The counts matrix of scRNA, the rowname is GeneSymbol, the colname is cell names.
#' @param annotation The annotation of scRNA with two columns, first column is cell name(matched with counts),
#' second column is cell type (tumor or non-tumor).
#' @importFrom Seurat CreateSeuratObject NormalizeData GetAssayData
#'
#' @return A dataframe
#' @export
#'
#' @examples


Custom_preprocessing <- function(counts, annotation, project = "PanCanCell", min.cells = 400, min.features = 0,
                                 normalization.method = "LogNormalize", scale.factor = 10000,
                                 verbose = TRUE){
  library(Seurat)
  data("feature",package = "PreCanCell") # loading feature datasets
  data <- CreateSeuratObject(counts = counts, project = project, min.cells = min.cells, min.features = min.features)
  data <- NormalizeData(object = data, normalization.method = normalization.method, scale.factor = scale.factor, verbose = verbose)
  data <- GetAssayData(data, slot="data")
  data <- data[rownames(data) %in% feature$Symbol,]
  data <- data[,which(Matrix::colSums(data)!=0)]

  if (dim(data)[1]!=259)
    stop("There are not enough feature genes")
  
  data <- apply(t(as.data.frame(data)), 2, function(x) {
    (x - min(x)) / (max(x) - min(x))
  })
  data <- as.data.frame(CancerLabel=annotation[,2],data)
  return(data)
}
