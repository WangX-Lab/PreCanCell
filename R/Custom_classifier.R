#' Add custom training datasets to perform prediction
#'
#' @description \code{Custom_classifier} Preprocessing for custom training datasets.
#'
#' @param testdata Cell type to be predicted.
#' (Note: Gene expression values of the samples need to be scaled to the range [0,1] by the function \link{PreCanCell_data} first)
#' @param CustomDataList A list of training datasets(added to TrainDataList).
#' @param cores Number of cores for parallel computing.
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom fastknn fastknn
#'
#' @return A dataframe with 4 columns:
#' \item{Sample}{Cell id}
#' \item{freq_cancer}{the ratio of the number of classifiers predicting that the cell is a cancer cell to 5 classifiers}
#' \item{freq_non_cancer}{the ratio of the number of classifiers predicting that the cell is a non-cancer cell to 5 classifiers}
#' \item{pred_labels}{Prediction results}
#'
#' @export
#'
#' @examples


Custom_classifier <- function(testdata, cores, CustomDataList) {

  CustomDataList <- c(CustomDataList,TrainDataList) # Add to TrainDataList

  ## Check arguments
  if (missing(testdata) || !class(testdata) %in% c("matrix", "data.frame") || dim(testdata)[2] < 100)
    stop("'testdata' is missing or incorrect")

  if (length(CustomDataList)%%2 == 0)
    stop("The number of TrainingDataList must be odd")

  ## Predict malignant and non-malignant cells
  set.seed(1234)

  cl <- makeCluster(cores) ## Threading enables steps to run in parallel
  registerDoParallel(cl)
  res <- suppressWarnings(
    foreach(i = seq_len(length(CustomDataList)), .combine = "cbind", .packages = "fastknn", .export = "testdata", .errorhandling = "stop") %dopar% {
      CustomData <- CustomDataList[[i]]
      CustomData <- CustomData[,c(1, match(colnames(testdata), colnames(CustomData)))]
      pred_knn <- fastknn(data.matrix(CustomData[,-1, drop = FALSE]),as.factor(CustomData[,1]), data.matrix(testdata), k = 5)
      label <- pred_knn$class
    })
  stopCluster(cl)

  if (length(CustomDataList)==1){
    result <- data.frame(
      Sample = rownames(testdata),
      pred_label = res,
      row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE
    )
  }
  else {
    pred <- apply(res, 1, function(x) c({sum(x == "1")/length(CustomDataList)}, {sum(x == "2")/length(CustomDataList)}, {names(which.max(table(x)))}))
    pred <- t(pred)
    result <- data.frame(Sample = rownames(testdata),
                         freq_cancer = pred[,1],
                         freq_non_cancer = pred[,2],
                         pred_label = factor(pred[,3], levels = c("1","2"),labels = c("cancer","non_cancer")),
                         row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
  }

  return(result)

}

