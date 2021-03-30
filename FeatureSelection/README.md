

# Evaluating feature selection on cross domain prediction


## Requirements
- R
- Required R packages: glmnet, caret, doMC, PharmacoGx, SummarizedExperiment, mRMRe

## Running
Download the psets from https://orcestra.ca and and specify the path to the psets in line 6 of `LoadData.R`.<br/>
For l1000 results, set the fs.method to `l1000` in `cross_domain.R`, and for mRMR set it to `mRMR`.<br/>
If you move the the l1000.geneList.txt file, you need to provide the absolute path in line 81 of `LoadData.R`.<br/>
Finally, you can run the script by the following commands within the R enviroment. <br/>
```
source("cross_domain.R")
```
Alternatively, you can run it from your terminal by the following command. <br/>
```
Rscript cross_domain.R > output.txt
```
The results will be printed to the `output.txt` file.

This program was successfully tested using the following OS and packages:<br/>
<i>Linux Mint 19.2 Cinnamon x_64<br/>
R version 4.0.3 (2020-10-10)<br/>
glmnet 4.0-2<br/>
caret 6.0-86<br/>
doMC 1.3.6<br/>
PharmacoGx 2.0.3<br/>
SummarizedExperiment 1.18.1<br/>
mRMRe 2.1.0<br/>
</i>
