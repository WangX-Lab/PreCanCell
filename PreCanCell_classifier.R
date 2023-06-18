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
#' @importFrom foreach %dopar%
#' @importFrom fastknn fastknn
#'
#' @return A dataframe with 2 columns:
#' \item{Sample}{Cell id}
#' \item{pred_labels}{Prediction results}
#'
#' @export
#'
#' @examples
#' path <- system.file("extdata", "example.txt", package = "PreCanCell", mustWork = TRUE)
#' input <- read.table(path, stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
#' testdata <- PreCanCell_data(input)
#' results <- PreCanCell_classifier(testdata, 2)


PreCanCell_classifier <- function(testdata, cores) {

  ## Check arguments
  if (missing(testdata) || !class(testdata) %in% c("matrix", "data.frame") || dim(testdata)[2] < 100)
    stop("'testdata' is missing or incorrect")

  ## Predict malignant and non-malignant cells
  set.seed(1234)

  cl <- makeCluster(cores) ## Threading enables steps to run in parallel
  registerDoParallel(cl)
  res <- suppressWarnings(
    foreach(i = seq_len(length(TrainDataList)), .combine = "cbind", .packages = "fastknn", .export = c("TrainDataList","testdata"), .errorhandling = "stop") %dopar% {
      TrainData <- TrainDataList[[i]]
      TrainData <- TrainData[,c(1, match(colnames(testdata), colnames(TrainData)))]
      TrainData[-1] <- apply(TrainData[-1], 2, function(x) {
        (x - min(x)) / (max(x) - min(x))
      })
      pred_knn <- fastknn(data.matrix(TrainData[,-1, drop = FALSE]),as.factor(TrainData[,1]), data.matrix(testdata), k = 5)
      label <- pred_knn$class
    })
  stopCluster(cl)

  pred <- apply(res, 1, function(x) c({sum(x == "1")/5}, {sum(x == "2")/5}, {names(which.max(table(x)))}))
  pred <- t(pred)

  result <- data.frame(Sample = rownames(testdata),
                       freq_cancer = pred[,1],
                       freq_non_cancer = pred[,2],
                       pred_label = factor(pred[,3], levels = c("1","2"),labels = c("cancer","non_cancer")),
                       row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)

  return(result)

}


