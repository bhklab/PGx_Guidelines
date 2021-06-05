#############################################
# Author: Soheil Jahangiri-Tazehkand - BHKLab
#############################################
rm(list = ls())

no.cores <- 6 # Number of CPU cores to use.

source("LoadData.R")

set.seed(42)

#### Feature Selection method can be "mRMR" or "Landmark"
# mRMR: 200 high scoring protein coding genes are selected. You can change the number high scoring genes in the code.
# Landmark: l1000 landmark genes are selected
fs.method <- "mRMR" 

#### Response can be "AAC" or "IC50"
# AAC: Area above the viability-concentration curve will be used as the summary of drug response
# IC50: Log2 of the half maximal inhibitory concentration (IC50) will be used as the summary of drug response
response <- "IC50" 


############################################################################
if (fs.method == "Landmark"){
  subset_genes <- geneId2Name$Name %in% l1000$V1
}else if (fs.method == "mRMR"){
  subset_genes <- geneId2Name$ID %in% protein.coding.genes$ID
}

ctrp_exprs <- ctrp_exprs[, subset_genes]
gdsc_exprs <- gdsc_exprs[, subset_genes]
gcsi_exprs <- gcsi_exprs[, subset_genes]

#############################################################################

if (response == "AAC"){
  ctrp.resp_mat <- ctrp_sens_aac
  gdsc.resp_mat <- gdsc_sens_aac
  gcsi.resp_mat <- gcsi_sens_aac
}else if (response == "IC50"){
  ctrp.resp_mat <- log2(ctrp_sens_ic50)
  ctrp.resp_mat[is.infinite(ctrp.resp_mat)] <- NA
  gdsc.resp_mat <- log2(gdsc_sens_ic50)
  gdsc.resp_mat[is.infinite(gdsc.resp_mat)] <- NA
  gcsi.resp_mat <- log2(gcsi_sens_ic50)
  gcsi.resp_mat[is.infinite(gcsi.resp_mat)] <- NA
}

#############################################################################
resp_mat <- ctrp.resp_mat
exprs_mat <- ctrp_exprs
#############################################################################
#### Removing haematopoietic cell lines and those with unmapped tissue types
ccle.tissueTypes <- ccle.tissueTypes[-which(ccle.tissueTypes$tissueid %in% c("Lymphoid", "Myeloid", "Other")), c(1, 2), drop = FALSE]
non_haema <- intersect(rownames(ccle.tissueTypes), rownames(exprs_mat))

exprs_mat <- exprs_mat[non_haema, ]
resp_mat <- resp_mat[non_haema, ]

#####################  ADD TISSUE TYPE AS FEATURE TO CTRP #####################
###############################################################################
# Only keep necessary tissue types
ccle.tissueTypes.subset <- ccle.tissueTypes[rownames(exprs_mat), ]
# Mapping categorical tissue id to numerical tissue id
for (tissueType in unique(ccle.tissueTypes.subset$tissueid)){
  ccle.tissueTypes.subset[tissueType] <- ifelse(ccle.tissueTypes.subset$tissueid == tissueType, 1, 0)
}
# Bind tissue type as a feature to expression matrix
exprs_mat <- cbind(as.matrix(ccle.tissueTypes.subset[rownames(exprs_mat), -c(1, 2)]), exprs_mat)
###################################################################################################
#####################  SUBSET FEATURES AND ADD TISSUE TYPE AS FEATURE TO GDSC #####################
# Keep GDSC cell lines with tissue ids in CTRP
gdsc.keep <- gdsc.tissueTypes[gdsc.tissueTypes$tissueid %in% unique(ccle.tissueTypes.subset$tissueid), ]$cellid
gdsc.keep <- intersect(gdsc.keep, rownames(gdsc_exprs))
gdsc_exprs <- gdsc_exprs[gdsc.keep, ] # features are either PC or l1000
gdsc.tissueTypes.subset <- gdsc.tissueTypes[gdsc.keep, ]

# Assign tissue id numbers to GDSC tissue types
for (tissueType in unique(ccle.tissueTypes.subset$tissueid)){
  gdsc.tissueTypes.subset[tissueType] <- ifelse(gdsc.tissueTypes.subset$tissueid == tissueType, 1, 0)
}
# Add tissue id number as a feature to expression matrix
gdsc_exprs <- cbind(as.matrix(gdsc.tissueTypes.subset[rownames(gdsc_exprs), -c(1, 2)]), gdsc_exprs)
###################################################################################################
#####################  SUBSET FEATURES AND ADD TISSUE TYPE AS FEATURE TO gCSI #####################
# Keep gCSI cell lines with tissue ids in CTRP
gcsi.keep <- gcsi.tissueTypes[gcsi.tissueTypes$tissueid %in% unique(ccle.tissueTypes.subset$tissueid), ]$cellid
gcsi.keep <- intersect(gcsi.keep, rownames(gcsi_exprs))
gcsi_exprs <- gcsi_exprs[gcsi.keep, ] # features are either PC or l1000
gcsi.tissueTypes.subset <- gcsi.tissueTypes[gcsi.keep, ]

# Assign tissue id numbers to gCSI tissue types
for (tissueType in unique(ccle.tissueTypes.subset$tissueid)){
  gcsi.tissueTypes.subset[tissueType] <- ifelse(gcsi.tissueTypes.subset$tissueid == tissueType, 1, 0)
}

# Add tissue id number as a feature to expression matrix
gcsi_exprs <- cbind(as.matrix(gcsi.tissueTypes.subset[rownames(gcsi_exprs), -c(1, 2)]), gcsi_exprs)
###############################################################################

predictions.gdsc <- list()
predictions.gcsi <- list()

for (drug.id in 1:ncol(resp_mat)){
  drugName <- colnames(resp_mat)[drug.id]
  print(drugName)
  cellLines <- complete.cases(resp_mat[, drugName])
  sensVec <- resp_mat[cellLines, drugName]
  exprsMat.cc <- exprs_mat[cellLines, ]
  
  #### Scale the training expression matrix
  preproc.model <- suppressWarnings(caret::preProcess(exprsMat.cc, method = c("center", "scale")))
  train.Exp.trans <- predict(preproc.model, exprsMat.cc)
  
  # Normalizing test matrices by mean and variance of the training data .
  gdsc_exprs.trans <- predict(preproc.model, gdsc_exprs)
  gcsi_exprs.trans <- predict(preproc.model, gcsi_exprs)
  
  if (fs.method == "Landmark"){ # if feature selection is set to Landmark, we have already subseted to l1000 genes above
    features <- colnames(train.Exp.trans)
  }else if (fs.method  == "mRMR"){
    dataMat <- mRMRe::mRMR.data(data=as.data.frame(cbind(sensVec, train.Exp.trans), stringAsFactor=FALSE))
    features <- mRMRe::mRMR.ensemble(data = dataMat, target_indices = 1, solution_count = 1, feature_count = 200)
    features <- features@feature_names[unlist(features@filters)]
  }

  # Subset train matrix to selected features
  train.Exp.trans <- train.Exp.trans[, features]
  
  # Subset test matrices to selected features
  gdsc_exprs.trans <- gdsc_exprs.trans[, features]
  gcsi_exprs.trans <- gcsi_exprs.trans[, features]
  
  predictions.gdsc[[drugName]] <- vector(mode = "numeric", length = nrow(gdsc_exprs.trans))
  predictions.gcsi[[drugName]] <- vector(mode = "numeric", length = nrow(gcsi_exprs.trans))
  
  # Training and Prediction
  doMC::registerDoMC(cores = no.cores)
  fit <- glmnet::cv.glmnet(x = train.Exp.trans, y = sensVec, alpha = 0, parallel = TRUE)
  
  predictions.gdsc[[drugName]] <- predict(fit, s = fit$lambda.min, newx = gdsc_exprs.trans)
  predictions.gcsi[[drugName]] <- predict(fit, s = fit$lambda.min, newx = gcsi_exprs.trans)
}

###############################################################################################################
##### Evaluate the predictions using Pearson and Spearman correlations, and RMSE
###############################################################################################################
RMSE  <-  function(P, L){
  sqrt(mean((P - L)^2))
}

drugs <- c("Bortezomib", "Entinostat", "Sirolimus", "Docetaxel", "Gemcitabine", "Crizotinib", "Lapatinib", "Vorinostat", "Erlotinib", "Paclitaxel", "Pictilisib")

results <- data.frame(Drug = character(length = 0),
                      Response = character(length = 0),
                      Method = character(length = 0),
                      Evaluation = character(length = 0),
                      GDSCv2 = numeric(length = 0),
                      gCSI = numeric(length = 0))

for (drugName in drugs){
  print(drugName)
  cellLines <- complete.cases(resp_mat[, drugName])
  sensVec <- resp_mat[cellLines, drugName]
  # exprsMat.cc <- exprs_mat[cellLines, ]
  # ctrp.cellLines <- rownames(exprs_mat[cellLines, ])
  
  gdsc.prediction.cellLines <- rownames(as.matrix(predictions.gdsc[[drugName]]))
  gcsi.prediction.cellLines <- rownames(as.matrix(predictions.gcsi[[drugName]]))
  
  # Getting the ground truth vectors for the drug and the cell lines we have response for
  gt.resp.gdsc <- gdsc.resp_mat[gdsc.prediction.cellLines, drugName]
  gt.resp.gdsc <- gt.resp.gdsc[!is.na(gt.resp.gdsc)]
  gt.resp.gcsi <- gcsi.resp_mat[gcsi.prediction.cellLines, drugName]
  gt.resp.gcsi <- gt.resp.gcsi[!is.na(gt.resp.gcsi)]
  
  gdsc.cellLines <- names(gt.resp.gdsc)
  gcsi.cellLines <- names(gt.resp.gcsi)
  
  res.all_gdsc.pcc <- cor(gt.resp.gdsc, as.matrix(predictions.gdsc[[drugName]])[gdsc.cellLines, 1], method = "pearson")
  res.all_gdsc.sp <- cor(gt.resp.gdsc, as.matrix(predictions.gdsc[[drugName]])[gdsc.cellLines, 1], method = "spearman")
  res.all_gdsc.rmse <- RMSE(gt.resp.gdsc, as.matrix(predictions.gdsc[[drugName]])[gdsc.cellLines, 1])
  res.all_gcsi.pcc <- cor(gt.resp.gcsi, as.matrix(predictions.gcsi[[drugName]])[gcsi.cellLines, 1], method = "pearson")
  res.all_gcsi.sp <- cor(gt.resp.gcsi, as.matrix(predictions.gcsi[[drugName]])[gcsi.cellLines, 1], method = "spearman")
  res.all_gcsi.rmse <- RMSE(gt.resp.gcsi, as.matrix(predictions.gcsi[[drugName]])[gcsi.cellLines, 1])
  
  results <- rbind(results, data.frame(Drug = drugName, Response = response, Method = fs.method, Evaluation = "Pearson", GDSCv2 = round(res.all_gdsc.pcc, 2), gCSI = round(res.all_gcsi.pcc, 2)))
  results <- rbind(results, data.frame(Drug = drugName, Response = response, Method = fs.method, Evaluation = "Spearman", GDSCv2 = round(res.all_gdsc.sp, 2), gCSI = round(res.all_gcsi.sp, 2)))
  results <- rbind(results, data.frame(Drug = drugName, Response = response, Method = fs.method, Evaluation = "RMSE", GDSCv2 = round(res.all_gdsc.rmse, 2), gCSI = round(res.all_gcsi.rmse, 2)))
  
}
