#' Data Preprocessing
#'
#' @description \code{PreCanCell_data} This function selects matched genes and scales gene expression values (maximum=1, minimum=0).
#'
#' @details The "input" refers to a normalized expression data with gene symbol.
#' Please note that \code{feature} must be present in row of \code{input} and
#' should not contain any duplicated feature names.
#' @param input The normalized expression dataframe, with rows being genes, and columns being cells.
#' @return A m*n matrix contains the expression values of n matched genes in m cells.
#' @export
#'
#' @examples
#' path <- system.file("extdata", "example.txt", package = "PreCanCell", mustWork = TRUE)
#' input <- read.table(path, stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, sep = "\t", quote = "", row.names = 1)
#' testdata <- PreCanCell_data(input)


PreCanCell_data <- function(input) {

  ## Check arguments
  if (missing(input) || class(input) != "data.frame"|| dim(input)[2] < 5)
    stop("'input' data is missing or incorrect")

  ## Filter and min-max scaling
  if (TRUE %in% is.na(match(feature[,"Symbol"], rownames(input)))) {
    stop("Predictor variables with missing values are presented in the current input data")
  }
  input <- input[match(feature[,"Symbol"], rownames(input)),, drop = FALSE]

  input <- apply(t(input), 2, function(x) {
    (x - min(x)) / (max(x) - min(x))
  })
  input[is.nan(input)] <- 0

  return(input)

}
