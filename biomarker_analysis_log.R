# biomarker analysis - excluding non-solid tissues and log2 on ic50
library(PharmacoGx); library(SummarizedExperiment); library(knitr); library(readxl)

# loading PSets
CTRPv2 <- readRDS("../code/PSets/CTRPv2.rds")
GDSC2 <- readRDS("../code/PSets/GDSC2.rds")
CCLE <- readRDS("../code/PSets/CCLE.rds")
gCSI <- readRDS("../code/PSets/gCSI.rds")
GDSC1 <- readRDS("../code/PSets/GDSC1.rds")

# intersect PSets by common drugs
common <- PharmacoGx::intersectPSet(list('gCSI' = gCSI, 'GDSC2' = GDSC2, 'CTRPv2' = CTRPv2, 'GDSCv1' = GDSC1), intersectOn = c("drugs"))
ctrp <- common[["CTRPv2"]] # no molprof info available
gdsc2 <- common[["GDSC2"]] 
gcsi <- common[["gCSI"]]
gdsc1 <- common[["GDSCv1"]]
rm(common)

common_drugs <- rownames(ctrp@drug)

#use uniprot id from biomarker list to get ensembl genes
gene_map <- read_excel("genes_annotations.xls")
genes <- na.omit(gene_map$ensemblvar)

# mapping between CCLE and CTRPv2 cell lines
ccle_celllines <- CCLE@molecularProfiles[["Kallisto_0.46.1.rnaseq"]]@colData@listData[["Cell_Line"]]
ctrp_celllines <- ctrp@cell$cellid

ccle_molprof <- as.data.frame(CCLE@molecularProfiles[["Kallisto_0.46.1.rnaseq"]]@colData@listData)
common_celllines <- subset(ccle_molprof, ctrp_celllines %in% ccle_celllines) # rnaseq info for cell lines tested in CTRPv2
rownames(common_celllines) <- common_celllines$cellid

ctrp@molecularProfiles[["Kallisto_0.46.1.rnaseq"]] <- CCLE@molecularProfiles[["Kallisto_0.46.1.rnaseq"]]
ctrp_non_haema <- ctrp@cell[-which(ctrp@cell$tissueid == "Lymphoid" | ctrp@cell$tissueid == "Myeloid" | ctrp@cell$tissueid == "Other" ),]

# drug sensitivity signatures
#CTRP
#log
ctrp@sensitivity[["profiles"]][["ic50_recomputed"]] <- log2(ctrp@sensitivity[["profiles"]][["ic50_recomputed"]])
ctrp.aac.sigs <- PharmacoGx::drugSensitivitySig(ctrp, "Kallisto_0.46.1.rnaseq", sensitivity.measure = "aac_recomputed", drugs=common_drugs, cells = rownames(ctrp_non_haema), features = genes)
ctrp.ic50.sigs <- PharmacoGx::drugSensitivitySig(ctrp, "Kallisto_0.46.1.rnaseq", sensitivity.measure = "ic50_recomputed", drugs=common_drugs, cells = rownames(ctrp_non_haema), features = genes)

#GDSC
#log
gdsc2_non_haema <- GDSC2@cell[-which(GDSC2@cell$tissueid == "Lymphoid" | GDSC2@cell$tissueid == "Myeloid" | GDSC2@cell$tissueid == "Other"),]

GDSC2@sensitivity[["profiles"]][["ic50_recomputed"]] <- log2(GDSC2@sensitivity[["profiles"]][["ic50_recomputed"]])
gdsc2.aac.sigs <- PharmacoGx::drugSensitivitySig(GDSC2, "Kallisto_0.46.1.rnaseq", sensitivity.measure = "aac_recomputed", drugs= common_drugs, cells = rownames(gdsc2_non_haema), features = genes)
gdsc2.ic50.sigs <- PharmacoGx::drugSensitivitySig(GDSC2, "Kallisto_0.46.1.rnaseq", sensitivity.measure = "ic50_recomputed", drugs= common_drugs, cells = rownames(gdsc2_non_haema),  features = genes)

#gCSI
#log of ic50_recomputed
gcsi_non_haema <- gCSI@cell[-which(gCSI@cell$tissueid == "Lymphoid" | gCSI@cell$tissueid == "Myeloid" | gCSI@cell$tissueid == "Other"),]
gCSI@sensitivity[["profiles"]][["ic50_recomputed"]] <- log2(gCSI@sensitivity[["profiles"]][["ic50_recomputed"]])
gcsi.aac.sigs <- PharmacoGx::drugSensitivitySig(gCSI, "Kallisto_0.46.1.rnaseq", sensitivity.measure = "aac_recomputed", drugs= common_drugs, cells = rownames(gcsi_non_haema), features = genes)
gcsi.ic50.sigs <- PharmacoGx::drugSensitivitySig(gCSI, "Kallisto_0.46.1.rnaseq", sensitivity.measure = "ic50_recomputed", drugs= common_drugs, cells = rownames(gcsi_non_haema), features = genes)

# plots 
ctrp.ic50.sigs <- as.data.frame(ctrp.ic50.sigs@.Data)
gdsc2.ic50.sigs <- as.data.frame(gdsc2.ic50.sigs@.Data)
gcsi.ic50.sigs <- as.data.frame(gcsi.ic50.sigs@.Data)

ctrp.aac.sigs <- as.data.frame(ctrp.aac.sigs@.Data)
gdsc2.aac.sigs <- as.data.frame(gdsc2.aac.sigs@.Data)
gcsi.aac.sigs <- as.data.frame(gcsi.aac.sigs@.Data)

# Lapatinib
lap = as.data.frame(gene_map[which(gene_map$compound == "Lapatinib"), ])
lap$ctrp.ic50 <- NA; lap$gdsc2.ic50 <- NA; lap$gcsi.ic50 <- NA
lap$ctrp.aac <- NA; lap$gdsc2.aac <- NA; lap$gcsi.aac <- NA

y <- "Lapatinib.estimate"
for (row in rownames(lap)) {
  x <- lap[row, "ensemblvar"]
  
  lap[row, "ctrp.ic50"] <- ctrp.ic50.sigs[x, y]
  lap[row, "gdsc2.ic50"] <- gdsc2.ic50.sigs[x, y]
  lap[row, "gcsi.ic50"] <- gcsi.ic50.sigs[x, y]
  lap[row, "ctrp.aac"] <- ctrp.aac.sigs[x, y]
  lap[row, "gdsc2.aac"] <- gdsc2.aac.sigs[x, y]
  lap[row, "gcsi.aac"] <- gcsi.aac.sigs[x, y]
}

#Docetaxel
doc = as.data.frame(gene_map[which(gene_map$compound == "Docetaxel"), ])
doc$ctrp.ic50 <- NA; doc$gdsc2.ic50 <- NA; doc$gcsi.ic50 <- NA
doc$ctrp.aac <- NA; doc$gdsc2.aac <- NA; doc$gcsi.aac <- NA

y <- "Docetaxel.estimate"
for (row in rownames(doc)) {
  x <- doc[row, "ensemblvar"]
  doc[row, "ctrp.ic50"] <- ctrp.ic50.sigs[x, y]
  doc[row, "gdsc2.ic50"] <- gdsc2.ic50.sigs[x, y]
  doc[row, "gcsi.ic50"] <- gcsi.ic50.sigs[x, y]
  doc[row, "ctrp.aac"] <- ctrp.aac.sigs[x, y]
  doc[row, "gdsc2.aac"] <- gdsc2.aac.sigs[x, y]
  doc[row, "gcsi.aac"] <- gcsi.aac.sigs[x, y]
}

#Pictilisib
pic = as.data.frame(gene_map[which(gene_map$compound == "Pictilisib"), ])
pic$ctrp.ic50 <- NA; pic$gdsc2.ic50 <- NA; pic$gcsi.ic50 <- NA
pic$ctrp.aac <- NA; pic$gdsc2.aac <- NA; pic$gcsi.aac <- NA
y <- "Pictilisib.estimate"
for (row in rownames(pic)) {
  x <- pic[row, "ensemblvar"]
  pic[row, "ctrp.ic50"] <- ctrp.ic50.sigs[x, y]
  pic[row, "gdsc2.ic50"] <- gdsc2.ic50.sigs[x, y]
  pic[row, "gcsi.ic50"] <- gcsi.ic50.sigs[x, y]
  pic[row, "ctrp.aac"] <- ctrp.aac.sigs[x, y]
  pic[row, "gdsc2.aac"] <- gdsc2.aac.sigs[x, y]
  pic[row, "gcsi.aac"] <- gcsi.aac.sigs[x, y]
}

#Gemcitabine
gem = as.data.frame(gene_map[which(gene_map$compound == "Gemcitabine"), ])
gem$ctrp.ic50 <- NA; gem$gdsc2.ic50 <- NA; gem$gcsi.ic50 <- NA
gem$ctrp.aac <- NA; gem$gdsc2.aac <- NA; gem$gcsi.aac <- NA
y <- "Gemcitabine.estimate"
for (row in rownames(gem)) {
  x <- gem[row, "ensemblvar"]
  gem[row, "ctrp.ic50"] <- ctrp.ic50.sigs[x, y]
  gem[row, "gdsc2.ic50"] <- gdsc2.ic50.sigs[x, y]
  gem[row, "gcsi.ic50"] <- gcsi.ic50.sigs[x, y]
  gem[row, "ctrp.aac"] <- ctrp.aac.sigs[x, y]
  gem[row, "gdsc2.aac"] <- gdsc2.aac.sigs[x, y]
  gem[row, "gcsi.aac"] <- gcsi.aac.sigs[x, y]
}

#Vorinostat
vor = as.data.frame(gene_map[which(gene_map$compound == "Vorinostat"), ])
vor$ctrp.ic50 <- NA; vor$gdsc2.ic50 <- NA; vor$gcsi.ic50 <- NA
vor$ctrp.aac <- NA; vor$gdsc2.aac <- NA; vor$gcsi.aac <- NA
y <- "Vorinostat.estimate"
for (row in rownames(vor)) {
  x <- vor[row, "ensemblvar"]
  vor[row, "ctrp.ic50"] <- ctrp.ic50.sigs[x, y]
  vor[row, "gdsc2.ic50"] <- gdsc2.ic50.sigs[x, y]
  vor[row, "gcsi.ic50"] <- gcsi.ic50.sigs[x, y]
  vor[row, "ctrp.aac"] <- ctrp.aac.sigs[x, y]
  vor[row, "gdsc2.aac"] <- gdsc2.aac.sigs[x, y]
  vor[row, "gcsi.aac"] <- gcsi.aac.sigs[x, y]
  }

#Paclitaxel
pac = as.data.frame(gene_map[which(gene_map$compound == "Paclitaxel"), ])
pac$ctrp.ic50 <- NA; pac$gdsc2.ic50 <- NA; pac$gcsi.ic50 <- NA
pac$ctrp.aac <- NA; pac$gdsc2.aac <- NA; pac$gcsi.aac <- NA
y <- "Paclitaxel.estimate"
for (row in rownames(pac)) {
  x <- pac[row, "ensemblvar"]
  pac[row, "ctrp.ic50"] <- ctrp.ic50.sigs[x, y]
  pac[row, "gdsc2.ic50"] <- gdsc2.ic50.sigs[x, y]
  pac[row, "gcsi.ic50"] <- gcsi.ic50.sigs[x, y]
  pac[row, "ctrp.aac"] <- ctrp.aac.sigs[x, y]
  pac[row, "gdsc2.aac"] <- gdsc2.aac.sigs[x, y]
  pac[row, "gcsi.aac"] <- gcsi.aac.sigs[x, y]
}

#Crizotinib
cri = as.data.frame(gene_map[which(gene_map$compound == "Crizotinib"), ])
cri$ctrp.ic50 <- NA; cri$gdsc2.ic50 <- NA; cri$gcsi.ic50 <- NA
cri$ctrp.aac <- NA; cri$gdsc2.aac <- NA; cri$gcsi.aac <- NA
y <- "Crizotinib.estimate"
for (row in rownames(cri)) {
  x <- cri[row, "ensemblvar"]
  cri[row, "ctrp.ic50"] <- ctrp.ic50.sigs[x, y]
  cri[row, "gdsc2.ic50"] <- gdsc2.ic50.sigs[x, y]
  cri[row, "gcsi.ic50"] <- gcsi.ic50.sigs[x, y]
  cri[row, "ctrp.aac"] <- ctrp.aac.sigs[x, y]
  cri[row, "gdsc2.aac"] <- gdsc2.aac.sigs[x, y]
  cri[row, "gcsi.aac"] <- gcsi.aac.sigs[x, y]
}

#Erlotinib
erl = as.data.frame(gene_map[which(gene_map$compound == "Erlotinib"), ])
erl$ctrp.ic50 <- NA; erl$gdsc2.ic50 <- NA; erl$gcsi.ic50 <- NA
erl$ctrp.aac <- NA; erl$gdsc2.aac <- NA; erl$gcsi.aac <- NA
y <- "Erlotinib.estimate"
for (row in rownames(erl)) {
  x <- erl[row, "ensemblvar"]
  erl[row, "ctrp.ic50"] <- ctrp.ic50.sigs[x, y]
  erl[row, "gdsc2.ic50"] <- gdsc2.ic50.sigs[x, y]
  erl[row, "gcsi.ic50"] <- gcsi.ic50.sigs[x, y]
  erl[row, "ctrp.aac"] <- ctrp.aac.sigs[x, y]
  erl[row, "gdsc2.aac"] <- gdsc2.aac.sigs[x, y]
  erl[row, "gcsi.aac"] <- gcsi.aac.sigs[x, y]
}

biomarker_log <- do.call("rbind", list(lap, doc, pic, gem, vor, pac, cri, erl))
WriteXLS::WriteXLS(biomarker_log, "../results/biomarker_log.xls")
