library(SummarizedExperiment)
library(PharmacoGx)
library(readr)

tissueType_encoding <- read_csv("tissueType_encoding.csv")
tissueType_encoding<-tissueType_encoding[!duplicated(tissueType_encoding$tissueid),3:28]

memory.limit(size=12000)

drugs <- c("Bortezomib", "Entinostat", "Sirolimus","Docetaxel","Gemcitabine", "Crizotinib",  "Lapatinib","Vorinostat","Erlotinib","Paclitaxel","Pictilisib")

CCLE <- readRDS("~/PSet/CCLE.rds")
features <- featureInfo(CCLE, mDataType = "Kallisto_0.46.1.rnaseq")
features <- features[which(features$gene_type=="protein_coding"),]
protein_genes <- rownames(features) 
ctrp.exprs <- summarizeMolecularProfiles(pSet = CCLE, mDataType = "Kallisto_0.46.1.rnaseq")
ctrp.exprs.protein_genes <- ctrp.exprs@assays@data$expr[protein_genes,]
ctrp.cell.info <- PharmacoGx::cellInfo(CCLE)
ttypes <- merge(ctrp.cell.info, tissueType_encoding, by = "tissueid")
comm <- intersect(ttypes$cellid, colnames(ctrp.exprs.protein_genes))
ctrp.exprs.protein_genes<-rbind(ctrp.exprs.protein_genes[,comm], t(ttypes[,24:48]))
ctrp.cell.info$Tumor <- ifelse(ctrp.cell.info$tissueid %in% c("Lymphoid", "Myeloid", "Other"), yes = 1, no = 0)
rm(CCLE)
rm(ctrp.exprs)


gCSI <- readRDS("~/PSet/gCSI.rds")
gcsi.exprs <- summarizeMolecularProfiles(pSet = gCSI, mDataType = "Kallisto_0.46.1.rnaseq")
gcsi.exprs.protein_genes <- gcsi.exprs@assays@data$expr[protein_genes,]
gcsi.aac <- PharmacoGx::summarizeSensitivityProfiles(
  pSet=gCSI,
  drugs = drugs,
  sensitivity.measure='aac_recomputed', 
  fill.missing = TRUE,
  summary.stat="median",
  verbose=FALSE)
gcsi.cell.info <- PharmacoGx::cellInfo(gCSI)
ttypes <- merge(gcsi.cell.info,tissueType_encoding, by = "tissueid")
comm <- intersect(ttypes$cellid,colnames(gcsi.exprs.protein_genes))
gcsi.exprs.protein_genes<-rbind(gcsi.exprs.protein_genes[,comm], t(ttypes[,39:63]))
gcsi.cell.info$Tumor <- ifelse(gcsi.cell.info$tissueid %in% c("Lymphoid", "Myeloid", "Other"), yes = 1, no = 0)
rm(gCSI)
rm(gcsi.exprs)

GDSCv2 <- readRDS("~/PSet/GDSC2.rds")
gdsc.exprs <- summarizeMolecularProfiles(pSet = GDSCv2, mDataType = "Kallisto_0.46.1.rnaseq")
gdsc.exprs.protein_genes <- gdsc.exprs@assays@data$expr[protein_genes,]
gdsc.aac <- PharmacoGx::summarizeSensitivityProfiles(
  pSet=GDSCv2,
  drugs = drugs,
  sensitivity.measure='aac_recomputed', 
  fill.missing = TRUE,
  summary.stat="median",
  verbose=FALSE)
gdsc.cell.info <- PharmacoGx::cellInfo(GDSCv2)
ttypes <- merge(gdsc.cell.info,tissueType_encoding, by = "tissueid")
comm <- intersect(ttypes$cellid,colnames(gdsc.exprs.protein_genes))
gdsc.exprs.protein_genes<-rbind(gdsc.exprs.protein_genes[,comm], t(ttypes[,22:46]))
gdsc.cell.info$Tumor <- ifelse(gdsc.cell.info$tissueid %in% c("Lymphoid", "Myeloid", "Other"), yes = 1, no = 0)
rm(GDSCv2)
rm(gdsc.exprs)

CTRPv2 <- readRDS("~/PSet/CTRPv2.rds")
ctrp.aac <- PharmacoGx::summarizeSensitivityProfiles(
  pSet=CTRPv2,
  drugs = drugs,
  sensitivity.measure='aac_recomputed', 
  fill.missing = TRUE,
  summary.stat="median",
  verbose=FALSE)
rm(CTRPv2)

GDSC1 <- readRDS("~/PSet/GDSC1.rds")
gdscv1.exprs <- summarizeMolecularProfiles(pSet = GDSC1, mDataType = "Kallisto_0.46.1.rnaseq")
gdscv1.exprs.protein_genes <- gdscv1.exprs@assays@data$expr[protein_genes,]
gdscv1.aac <- PharmacoGx::summarizeSensitivityProfiles(
  pSet=GDSC1,
  drugs = drugs,
  sensitivity.measure='aac_recomputed', 
  fill.missing = TRUE,
  summary.stat="median",
  verbose=FALSE)
gdscv1.cell.info <- PharmacoGx::cellInfo(GDSC1)
comm <- intersect(ttypes$cellid,colnames(gdscv1.exprs.protein_genes))
gdscv1.exprs.protein_genes<-rbind(gdscv1.exprs.protein_genes[,comm], t(ttypes[,22:46]))
gdscv1.cell.info$Tumor <- ifelse(gdscv1.cell.info$tissueid %in% c("Lymphoid", "Myeloid", "Other"), yes = 1, no = 0)
rm(gdscv1.exprs)
rm(GDSC1)


ctrp.exprs.protein_genes <- ctrp.exprs.protein_genes[ , !is.na(colSums(ctrp.exprs.protein_genes))]
gdsc.exprs.protein_genes <- gdsc.exprs.protein_genes[ , !is.na(colSums(gdsc.exprs.protein_genes))]
gdscv1.exprs.protein_genes <- gdscv1.exprs.protein_genes[ , !is.na(colSums(gdscv1.exprs.protein_genes))]
gcsi.exprs.protein_genes <- gcsi.exprs.protein_genes[ , !is.na(colSums(gcsi.exprs.protein_genes))]

ctrp.exprs.protein_genes <- na.omit(ctrp.exprs.protein_genes)
gdsc.exprs.protein_genes <- na.omit(gdsc.exprs.protein_genes)
gdscv1.exprs.protein_genes <- na.omit(gdscv1.exprs.protein_genes)
gcsi.exprs.protein_genes <- na.omit(gcsi.exprs.protein_genes)

write.table(gdsc.exprs.protein_genes, file = "GDSCv2.exprsALL.tsv", row.names=TRUE, sep="\t", col.names = TRUE)
write.table(gdscv1.exprs.protein_genes, file = "GDSCv1.exprsALL.tsv", row.names=TRUE, sep="\t", col.names = TRUE)
write.table(gcsi.exprs.protein_genes, file = "gCSI.exprsALL.tsv", row.names=TRUE, sep="\t", col.names = TRUE)
write.table(ctrp.exprs.protein_genes, file = "CTRP.exprsALL.tsv", row.names=TRUE, sep="\t", col.names = TRUE)

write.table(gdsc.aac, file = "GDSCv2.aacALL.tsv", row.names=TRUE, sep="\t", col.names = TRUE)
write.table(gdscv1.aac, file = "GDSCv1.aacALL.tsv", row.names=TRUE, sep="\t", col.names = TRUE)
write.table(gcsi.aac, file = "gCSI.aacALL.tsv", row.names=TRUE, sep="\t", col.names = TRUE)
write.table(ctrp.aac, file = "CTRP.aacALL.tsv", row.names=TRUE, sep="\t", col.names = TRUE)

write.table(gdsc.cell.info, file = "GDSCv2.infoALL.tsv", row.names=TRUE, sep="\t", col.names = TRUE)
write.table(gdscv1.cell.info, file = "GDSCv1.infoALL.tsv", row.names=TRUE, sep="\t", col.names = TRUE)
write.table(gcsi.cell.info, file = "gCSI.infoALL.tsv", row.names=TRUE, sep="\t", col.names = TRUE)
write.table(ctrp.cell.info, file = "CTRP.infoALL.tsv", row.names=TRUE, sep="\t", col.names = TRUE)

