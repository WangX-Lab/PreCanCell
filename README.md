# PreCanCell
A simple and effective ensemble learning algorithm for predicting cancer and non-cancer cells from single-cell transcriptomes.

# Description
PreCanCell first identified the differentially expressed genes (DEGs) between malignant and non-malignant cells commonly in five common cancer-associated single-cell transcriptome datasets. With each of the five datasets as the training set and the DEGs as the features, a single cell is classified as malignant or non-malignant by *k*-NN (*k* = 5). Finally, the single cell is classified by the majority vote of the five *k*-NN classification results.

# Details
+ The function `PreCanCell_data()` is used to data preprocessing. Its input should be normalized expression matrix with rownames being genes and colnames being cells.
+ The function `PreCanCell_classifier()` is used to identify malignant and non-malignant cells from single-cell transcriptomes, containing 2 parameters: testdata and cores.
  + "testdata" is a output matrix of the function `PreCanCell_data()`.
  + "cores" is the number of threads.

# Installation
Users can install the released version of PreCanCell with:
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("WangX-Lab/PreCanCell")
```

# Examples
```
## Data preprocessing (select matched genes and [0,1]-scaled gene expression values) --------------
library(PreCanCell)
path <- system.file("extdata", "example.txt", package = "PreCanCell", mustWork = TRUE)
input <- read.table(path, stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
testdata <- PreCanCell_data(input)
```

```
## Prediction of malignant and non-malignant cells -----------------------------------------------
library(PreCanCell)
results <- PreCanCell_classifier(testdata, 2)
```

# Contact
E-mail any questions to Xiaosheng Wang (xiaosheng.wang@cpu.edu.cn)
