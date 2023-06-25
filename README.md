# PreCanCell

A simple and effective ensemble learning algorithm for predicting cancer and non-cancer cells from single-cell transcriptomes.

<img width="1022" alt="image" src="https://github.com/WangX-Lab/PreCanCell/master/docs/PreCanCell_overview.png">

&nbsp;
&nbsp;
# Description

PreCanCell first identified the differentially expressed genes (DEGs) between malignant and non-malignant cells commonly in five common cancer-associated single-cell transcriptome datasets. With each of the five datasets as the training set and the DEGs as the features, a single cell is classified as malignant or non-malignant by *k*-NN (*k* = 5). Finally, the single cell is classified by the majority vote of the five *k*-NN classification results.

&nbsp;
&nbsp;
# Details

+ The function `PreCanCell_data()` is used to data preprocessing. Its input should be normalized expression matrix with rownames being genes and colnames being cells. The input data can be any library-depth normalization (e.g. TPM, CPM), but not log-transformed.
+ The function `PreCanCell_classifier()` is used to identify malignant and non-malignant cells from single-cell transcriptomes, containing 2 parameters: testdata and cores.
  + "testdata" is a output matrix of the function `PreCanCell_data()`.
  + "cores" is the number of threads.

&nbsp;
&nbsp;
# Installation

- The `Seurat` package (*version* >= 4.3.0) is used for data preprocessing.

- To install `PreCanCell` , first install the `fastknn` package, which can be installed as follows:
&nbsp;
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
    
devtools::install_github("WangX-Lab/PreCanCell")
```

- Finally, users can install the released version of `PreCanCell` with:
&nbsp;
```
devtools::install_github("WangX-Lab/PreCanCell")
```

&nbsp;
&nbsp;
# Examples

**Prepare data:**
&nbsp;
```
library(PreCanCell)
path <- system.file("extdata", "example.txt", package = "PreCanCell", mustWork = TRUE)
input <- read.table(path, stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
input[1:5,1:5]
```

|         |  0   |  1   |  2   |  3   |  4   |
| :-----: | :--: | :--: | :--: | :--: | :--: |
| SERINC2 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
|  PTPRF  | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| S100A1  | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
|  EFNA1  | 0.00 | 0.00 | 0.00 | 6.33 | 0.00 |
|  SOX4   | 2.59 | 0.00 | 0.00 | 5.23 | 0.00 |

**Data preprocessing (select matched genes and [0,1]-scaled gene expression values):**
&nbsp;
```
testdata <- PreCanCell_data(input)
testdata[1:5,1:5]
```

|      | SERINC2 | PTPRF | S100A16 | EFNA1 | SOX4 |
| :--: | :-----: | :---: | :-----: | :---: | :--: |
|  0   |  0.00   | 0.00  |  0.00   | 0.00  | 0.46 |
|  1   |  0.00   | 0.00  |  0.00   | 0.00  | 0.00 |
|  2   |  0.00   | 0.00  |  0.00   | 0.00  | 0.00 |
|  3   |  0.00   | 0.00  |  0.00   | 1.00  | 0.94 |
|  4   |  0.00   | 0.00  |  0.00   | 0.00  | 0.00 |


**Prediction of malignant and non-malignant cells:**
`cores` represents the number of cores to use for parallel execution. 
&nbsp;
```
results <- PreCanCell_classifier(testdata, 2)
head(results)
```

| Sample | freq_cancer | freq_non_cancer | pred_label |
| :----: | :---------: | :-------------: | :--------: |
|   0    |      1      |        0        |   cancer   |
|   1    |      0      |        1        | non_cancer |
|   2    |     0.8     |       0.2       |   cancer   |
|   3    |      1      |        0        |   cancer   |
|   4    |      0      |        1        | non_cancer |
|   5    |      0      |        1        | non_cancer |


# Vignettes
[Predicting malignant and non-malignant cells from single-cell transcriptomes](https://wangx-lab.github.io/PreCanCell/docs/index.html))

&nbsp;
&nbsp;
# Contact

E-mail any questions to Xiaosheng Wang (xiaosheng.wang@cpu.edu.cn)
