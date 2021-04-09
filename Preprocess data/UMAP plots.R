library(umap)
library(ggplot2)
library(SummarizedExperiment)
library(PharmacoGx)

memory.limit(size=12000)

drugs <- c("Bortezomib", "Entinostat", "Sirolimus","Docetaxel","Gemcitabine", "Crizotinib",  "Lapatinib","Vorinostat","Erlotinib","Paclitaxel","Pictilisib")

gCSI <- readRDS("~/PSet/gCSI.rds")
features <- featureInfo(gCSI, mDataType = "Kallisto_0.46.1.rnaseq")
features <- features[which(features$gene_type=="protein_coding"),]
protein_genes <- rownames(features) 
gcsi.exprs <- summarizeMolecularProfiles(pSet = gCSI, mDataType = "Kallisto_0.46.1.rnaseq")
gcsi.exprs.protein_genes <- gcsi.exprs@assays@data$expr[protein_genes,]
# gcsi.aac <- PharmacoGx::summarizeSensitivityProfiles(
#   pSet=gCSI,
#   drugs = drugs,
#   sensitivity.measure='aac_recomputed', 
#   fill.missing = TRUE,
#   summary.stat="median",
#   verbose=FALSE)
# gcsi.aac.comp <- gcsi.aac[ , !is.na(colSums(gcsi.aac))]
# gcsi.aac.comp <- gcsi.aac.comp[ ,-which(gcsi.cell.info$tissueid == "Other")]
# gcsi.cell.info = gcsi.cell.info[-which(gcsi.cell.info$tissueid == "Other"),]
# gcsi.cell.info$tissue.type <- 0
# gcsi.cell.info$tissue.type[which(gcsi.cell.info$tissueid == "Lymphoid" | gcsi.cell.info$tissueid == "Myeloid")] <- "non-solid"
# gcsi.cell.info$tissue.type[-which(gcsi.cell.info$tissueid == "Lymphoid" | gcsi.cell.info$tissueid == "Myeloid")] <- "solid"
# gcsi.umap.acc <- umap::umap(t(gcsi.aac.comp))
# gcsi.umap.acc.df <- as.data.frame(gcsi.umap.acc$layout)
# 
# tiff('gCSIumapresponse.png', units="in", width=10, height=7, res=600, compression = 'lzw')
# par(cex.axis=0.9)
# Tumor <- gcsi.cell.info[rownames(gcsi.umap.acc.df), c("tissue.type")]
# ggplot(gcsi.umap.acc.df, aes(V1, V2, color = Tumor))+ 
#   labs(x ="Umap 1", y = "Umap 2")+
#   geom_point(cex = 1.9)
# dev.off()

gcsi.exprs.protein_genes <- gcsi.exprs.protein_genes[ , !is.na(colSums(gcsi.exprs.protein_genes))]
gcsi.exprs.protein_genes <- na.omit(gcsi.exprs.protein_genes)
gcsi.cell.info <- PharmacoGx::cellInfo(gCSI)
gcsi.cell.info$tissue.type <- 0
gcsi.cell.info$tissue.type[which(gcsi.cell.info$tissueid == "Lymphoid" | gcsi.cell.info$tissueid == "Myeloid" | gcsi.cell.info$tissueid == "Other")] <- "non-solid"
gcsi.cell.info$tissue.type[-which(gcsi.cell.info$tissueid == "Lymphoid" | gcsi.cell.info$tissueid == "Myeloid" | gcsi.cell.info$tissueid == "Other")] <- "solid"

idx.gcsi <- intersect(row.names(gcsi.cell.info),colnames(gcsi.exprs.protein_genes))
gcsi.cell.info <- gcsi.cell.info[idx.gcsi,]
gcsi.exprs.protein_genes <- gcsi.exprs.protein_genes[,idx.gcsi]

gcsi.umap.exprs <- umap::umap(t(gcsi.exprs.protein_genes))
gcsi.umap.exprs.df <- as.data.frame(gcsi.umap.exprs$layout)

png('gCSIumapexprsv2.png', units="in", width=10, height=7, res=600)
par(cex.axis=0.9)
Tumor <- gcsi.cell.info[rownames(gcsi.umap.exprs.df), c("tissue.type")]
ggplot(gcsi.umap.exprs.df, aes(V1, V2, color = Tumor))+ 
  labs(x ="Umap 1", y = "Umap 2")+
  geom_point(cex = 1.9)#+geom_text(aes(label=ifelse(Tumor == "solid", row.names(gcsi.umap.exprs.df)," "),hjust=0, vjust=0))
dev.off()
rm(gCSI)
rm(gcsi.exprs)

CCLE <- readRDS("~/PSet/CCLE.rds")
ctrp.exprs <- summarizeMolecularProfiles(pSet = CCLE, mDataType = "Kallisto_0.46.1.rnaseq")
ctrp.exprs.protein_genes <- ctrp.exprs@assays@data$expr[protein_genes,]
ctrp.exprs.protein_genes <- ctrp.exprs.protein_genes[ , !is.na(colSums(ctrp.exprs.protein_genes))]
ctrp.exprs.protein_genes <- na.omit(ctrp.exprs.protein_genes)

ctrp.cell.info <- PharmacoGx::cellInfo(CCLE)
ctrp.cell.info$tissue.type <- 0
ctrp.cell.info$tissue.type[which(ctrp.cell.info$tissueid == "Lymphoid" | ctrp.cell.info$tissueid == "Myeloid" | ctrp.cell.info$tissueid == "Other")] <- "non-solid"
ctrp.cell.info$tissue.type[-which(ctrp.cell.info$tissueid == "Lymphoid" | ctrp.cell.info$tissueid == "Myeloid" | ctrp.cell.info$tissueid == "Other")] <- "solid"

idx.ctrp <- intersect(row.names(ctrp.cell.info),colnames(ctrp.exprs.protein_genes))
ctrp.cell.info <- ctrp.cell.info[idx.ctrp,]
ctrp.exprs.protein_genes <- ctrp.exprs.protein_genes[,idx.ctrp]

ctrp.umap.exprs <- umap::umap(t(ctrp.exprs.protein_genes))
ctrp.umap.exprs.df <- as.data.frame(ctrp.umap.exprs$layout)

png('CTRPv2umapexprsv2.png', units="in", width=10, height=7, res=600)
par(cex.axis=0.9)
Tumor <- ctrp.cell.info[rownames(ctrp.umap.exprs.df), c("tissue.type")]
ggplot(ctrp.umap.exprs.df, aes(V1, V2, color = Tumor))+ 
  labs(x ="Umap 1", y = "Umap 2")+
  geom_point(cex = 1.9)#+geom_text(aes(label=ifelse(Tumor == "solid", row.names(ctrp.umap.exprs.df)," "),hjust=0, vjust=0))
dev.off()
rm(CCLE)
rm(ctrp.exprs)

# CTRPv2 <- readRDS("~/PSet/CTRPv2.rds")
# ctrp.aac <- PharmacoGx::summarizeSensitivityProfiles(
#   pSet=CTRPv2,
#   drugs = drugs,
#   sensitivity.measure='aac_recomputed', 
#   fill.missing = TRUE,
#   summary.stat="median",
#   verbose=FALSE)
# ctrp.aac.comp <- ctrp.aac[ , !is.na(colSums(ctrp.aac))]
# ctrp.cell.info <- PharmacoGx::cellInfo(CTRPv2)
# ctrp.cell.info$tissue.type <- 0
# ctrp.cell.info$tissue.type[which(ctrp.cell.info$tissueid == "Lymphoid" | ctrp.cell.info$tissueid == "Myeloid" )] <- "non-solid"
# ctrp.cell.info$tissue.type[-which(ctrp.cell.info$tissueid == "Lymphoid" | ctrp.cell.info$tissueid == "Myeloid")] <- "solid"
# ctrp.umap.acc <- umap::umap(t(ctrp.aac.comp))
# ctrp.umap.acc.df <- as.data.frame(ctrp.umap.acc$layout)
# 
# tiff('CTRPv2umapresponse.png', units="in", width=10, height=7, res=600, compression = 'lzw')
# par(cex.axis=0.9)
# Tumor <- ctrp.cell.info[rownames(ctrp.umap.acc.df), c("tissue.type")]
# ggplot(ctrp.umap.acc.df, aes(V1, V2, color = Tumor))+ 
#   labs(x ="Umap 1", y = "Umap 2")+
#   geom_point(cex = 1.9)
# dev.off()
# rm(CTRPv2)

GDSCv2 <- readRDS("~/PSet/GDSC2.rds")
gdsc.exprs <- summarizeMolecularProfiles(pSet = GDSCv2, mDataType = "Kallisto_0.46.1.rnaseq")
gdsc.exprs.protein_genes <- gdsc.exprs@assays@data$expr[protein_genes,]
# gdsc.aac <- PharmacoGx::summarizeSensitivityProfiles(
#   pSet=GDSCv2,
#   drugs = drugs,
#   sensitivity.measure='aac_recomputed', 
#   fill.missing = TRUE,
#   summary.stat="median",
#   verbose=FALSE)
# gdsc.aac.comp <- gdsc.aac[ , !is.na(colSums(gdsc.aac))]
# gdsc.cell.info <- PharmacoGx::cellInfo(GDSCv2)
# gdsc.cell.info$tissue.type <- 0
# gdsc.cell.info$tissue.type[which(gdsc.cell.info$tissueid == "Lymphoid" | gdsc.cell.info$tissueid == "Myeloid" | gdsc.cell.info$tissueid == "Other")] <- "non-solid"
# gdsc.cell.info$tissue.type[-which(gdsc.cell.info$tissueid == "Lymphoid" | gdsc.cell.info$tissueid == "Myeloid" | gdsc.cell.info$tissueid == "Other")] <- "solid"
# gdsc.umap.acc <- umap::umap(t(gdsc.aac.comp))
# gdsc.umap.acc.df <- as.data.frame(gdsc.umap.acc$layout)
# 
# tiff('GDSCv2umapresponse.png', units="in", width=10, height=7, res=600, compression = 'lzw')
# par(cex.axis=0.9)
# Tumor <- gdsc.cell.info[rownames(gdsc.umap.acc.df), c("tissue.type")]
# ggplot(gdsc.umap.acc.df, aes(V1, V2, color = Tumor))+ 
#   labs(x ="Umap 1", y = "Umap 2")+
#   geom_point(cex = 1.9)
# dev.off()

gdsc.exprs.protein_genes <- gdsc.exprs.protein_genes[ , !is.na(colSums(gdsc.exprs.protein_genes))]
gdsc.exprs.protein_genes <- na.omit(gdsc.exprs.protein_genes)
gdsc.cell.info <- PharmacoGx::cellInfo(GDSCv2)
gdsc.cell.info$tissue.type <- 0
gdsc.cell.info$tissue.type[which(gdsc.cell.info$tissueid == "Lymphoid" | gdsc.cell.info$tissueid == "Myeloid" | gdsc.cell.info$tissueid == "Other")] <- "non-solid"
gdsc.cell.info$tissue.type[-which(gdsc.cell.info$tissueid == "Lymphoid" | gdsc.cell.info$tissueid == "Myeloid" | gdsc.cell.info$tissueid == "Other")] <- "solid"

idx.gdsc <- intersect(row.names(gdsc.cell.info),colnames(gdsc.exprs.protein_genes))
gdsc.cell.info <- gdsc.cell.info[idx.gdsc,]
gdsc.exprs.protein_genes <- gdsc.exprs.protein_genes[,idx.gdsc]

gdsc.umap.exprs <- umap::umap(t(gdsc.exprs.protein_genes))
gdsc.umap.exprs.df <- as.data.frame(gdsc.umap.exprs$layout)

png('GDSCv2umapexprsv2.png', units="in", width=10, height=7, res=600)
par(cex.axis=0.9)
Tumor <- gdsc.cell.info[rownames(gdsc.umap.exprs.df), c("tissue.type")]
ggplot(gdsc.umap.exprs.df, aes(V1, V2, color = Tumor))+ 
  labs(x ="Umap 1", y = "Umap 2")+
  geom_point(cex = 1.9)#+geom_text(aes(label=ifelse(Tumor == "solid", row.names(gdsc.umap.exprs.df)," "),hjust=0, vjust=0))
dev.off()
rm(GDSCv2)
rm(gdsc.exprs)


