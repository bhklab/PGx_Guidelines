

# Evaluating feature selection on cross domain prediction


## Requirements
- R
- Required R packages: glmnet, caret, doMC, PharmacoGx, SummarizedExperiment, mRMRe

## Running
Download the psets from https://orcestra.ca and and specify the path to the psets in line 6 of `LoadData.R`

For l1000 results, set the fs.method to `l1000` in `cross_domain.R`, and for mRMR set it to `mRMR`.

If you move the the l1000.geneList.txt file, you need to provide the absolute path in line 81 of `LoadData.R`.

Finally, you can run the script by the following commands within the R enviroment. 
```
source("cross_domain.R")
```
Alternatively, you can run it from your terminal by the following command. 
```
Rscript cross_domain.R > output.txt
```
The results will be printed to the `output.txt` file. 

This program was successfully tested using the following OS and packages:
Linux Mint 19.2 Cinnamon x_64
R version 4.0.3 (2020-10-10)
glmnet 4.0-2
caret 6.0-86
doMC 1.3.6
PharmacoGx 2.0.3
SummarizedExperiment 1.18.1
mRMRe 2.1.0
