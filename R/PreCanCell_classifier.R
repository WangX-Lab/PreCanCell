#' Predicting Cancer and Non-Cancer Cells from Single-cell Transcriptomes
#'
#' @description \code{PreCanCell_classifier} Identify malignant and non-malignant cells across cancer types.
#'
#' @param testdata Cell type to be predicted.
#' (Note: Gene expression values of the samples need to be scaled to the range [0,1] by the function \link{PreCanCell_data} first)
#' @param cores Number of cores for parallel computing.
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom stats predict
#'
#' @return A dataframe with 2 columns:
#' \item{Sample}{Cell id}
#' \item{pred_labels}{Prediction results}
#'
#' @export
#'
#' @examples
#' path <- system.file("extdata", "example.txt", package = "PreCanCell", mustWork = TRUE)
#' input <- read.table(path, stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, sep = "\t", quote = "", row.names = 1)
#' testdata <- PreCanCell_data(input)
#' results <- PreCanCell_classifier(testdata, 2)


PreCanCell_classifier <- function(testdata, cores) {

  ## Check arguments
  if (missing(testdata) || !class(testdata) %in% c("matrix", "data.frame"))
    stop("'testdata' is missing or incorrect")

  if (TRUE %in% is.na(feature[,"Symbol"] %in% colnames(testdata))) {
    stop("Predictor variables are missing or incorrect")
  }

  ## Predict malignant and non-malignant cells
  set.seed(1234)

  cl <- makeCluster(cores)
  registerDoParallel(cl)
  res <- suppressWarnings(
    foreach(i = seq_len(length(ModelList)), .combine = "cbind", .packages = "caret", .export = c("ModelList","testdata"), .errorhandling = "remove") %dopar% {
      model <- ModelList[[i]]
      pred <- predict(model, newdata = testdata)
    })
  stopCluster(cl)

  resC <- apply(res, 1, function(x){
    names(which.max(table(x)))
  })

  result <- data.frame(Sample = rownames(testdata), pred_labels = resC, row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)

  result$pred_labels <- factor(result$pred_labels,levels = c("1","2"),labels = c("cancer","non_cancer"))

  return(result)

}

