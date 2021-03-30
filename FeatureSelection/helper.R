#############################################
# Author: Soheil Jahangiri-Tazehkand - BHKLab
#############################################

# NOTE: pSets are downloaded from https://orcestra.ca
pathToPSets <- "../Data/psets/"
ccle <- readRDS(paste0(pathToPSets, "CCLE.rds"))
ctrp <- readRDS(paste0(pathToPSets, "CTRPv2.rds"))
gdsc <- readRDS(paste0(pathToPSets, "GDSC2.rds"))
gcsi <- readRDS(paste0(pathToPSets, "gCSI.rds"))

# =============================================================
# Intersecting to common drugs
# =============================================================
intersection <- PharmacoGx::intersectPSet(pSets = list(ctrp, gcsi, gdsc), intersectOn = ("drugs"))
ctrp <- intersection[["CTRPv2"]]
gdsc <- intersection[["GDSC_v2"]]
gcsi <- intersection[["gCSI"]]
rm(intersection)

# =============================================================
# Extracting sensitivity matrices (AAC and IC50) from pSets
# =============================================================
ctrp_sens_aac <- PharmacoGx::summarizeSensitivityProfiles(pSet = ctrp, sensitivity.measure = "aac_recomputed", summary.stat = "median")
gdsc_sens_aac <- PharmacoGx::summarizeSensitivityProfiles(pSet = gdsc, sensitivity.measure = "aac_recomputed", summary.stat = "median")
gcsi_sens_aac <- PharmacoGx::summarizeSensitivityProfiles(pSet = gcsi, sensitivity.measure = "aac_recomputed", summary.stat = "median")

# =============================================================
# Extracting gene expression profiles (TPM) from each pSet. 
# =============================================================
# For CTRP, we use CCLE to extract the gene expression profiles. 
ctrp_exprs <- PharmacoGx::summarizeMolecularProfiles(pSet = ccle, mDataType = "Kallisto_0.46.1.rnaseq", cell.lines = colnames(ctrp_sens_aac))
ctrp_exprs <- SummarizedExperiment::assay(ctrp_exprs)

gdsc_exprs <- PharmacoGx::summarizeMolecularProfiles(pSet = gdsc, mDataType = "Kallisto_0.46.1.rnaseq", cell.lines = colnames(gdsc_sens_aac))
gdsc_exprs <- SummarizedExperiment::assay(gdsc_exprs)

gcsi_exprs <- PharmacoGx::summarizeMolecularProfiles(pSet = gcsi, mDataType = "Kallisto_0.46.1.rnaseq", cell.lines = colnames(gcsi_sens_aac))
gcsi_exprs <- SummarizedExperiment::assay(gcsi_exprs)

# =============================================================
# Removing cell lines with no RNA expression in CTRP and 
# subseting the sensitivity matrices to cell lines having RNAseq
# =============================================================
ctrp_exprs <- ctrp_exprs[, !is.na(colSums(ctrp_exprs))]
ctrp_sens_aac <- ctrp_sens_aac[, colnames(ctrp_exprs)]

gdsc_exprs <- gdsc_exprs[, !is.na(colSums(gdsc_exprs))]
gdsc_sens_aac <- gdsc_sens_aac[, colnames(gdsc_exprs)]

gcsi_exprs <- gcsi_exprs[, !is.na(colSums(gcsi_exprs))]
gcsi_sens_aac <- gcsi_sens_aac[, colnames(gcsi_exprs)]

# =============================================================
# Transposing matrices to have cell lines on rows and features 
# and drugs on columns
# =============================================================
ctrp_exprs <- t(ctrp_exprs)
ctrp_sens_aac <- t(ctrp_sens_aac)

gcsi_exprs <- t(gcsi_exprs)
gcsi_sens_aac <- t(gcsi_sens_aac)

gdsc_exprs <- t(gdsc_exprs)
gdsc_sens_aac <- t(gdsc_sens_aac)

# =============================================================
# Reading cell line tissue type info from psets
# =============================================================
ft.info <- PharmacoGx::featureInfo(ccle, mDataType = "Kallisto_0.46.1.rnaseq")
geneId2Name <- data.frame(ID = ft.info$gene_id, Name = ft.info$gene_name, row.names = ft.info$gene_id)

ccle.tissueTypes <- PharmacoGx::cellInfo(ccle)[, c("cellid", "tissueid")]
gdsc.tissueTypes <- PharmacoGx::cellInfo(gdsc)[, c("cellid", "tissueid")]
gcsi.tissueTypes <- PharmacoGx::cellInfo(gcsi)[, c("cellid", "tissueid")]


ft.info.pc <- ft.info[ft.info$gene_type == "protein_coding", ]
protein.coding.genes <- data.frame(ID = ft.info.pc$gene_id, Name = ft.info.pc$gene_name, row.names = ft.info.pc$gene_id)

l1000 <- read.table(paste0("Hossein_Project/finalCodes/", "l1000.geneList.txt"))

rm(ctrp, ccle, gdsc, gcsi, ft.info, ft.info.pc)
gc(full = TRUE)


