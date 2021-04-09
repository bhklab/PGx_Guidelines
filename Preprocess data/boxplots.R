library(ggplot2)
library(reshape2)
library(ggpubr)
library(SummarizedExperiment)
library(PharmacoGx)

memory.limit(size=12000)

drugs <- c("Bortezomib", "Entinostat", "Sirolimus","Docetaxel","Gemcitabine", "Crizotinib",  "Lapatinib","Vorinostat","Erlotinib","Paclitaxel","Pictilisib")

gCSI <- readRDS("~/PSet/gCSI.rds")
gcsi.aac <- PharmacoGx::summarizeSensitivityProfiles(
  pSet=gCSI,
  drugs = drugs,
  sensitivity.measure='aac_recomputed', 
  fill.missing = TRUE,
  summary.stat="median",
  verbose=FALSE)
gcsi.cell.info <- PharmacoGx::cellInfo(gCSI)
df <- reshape2::melt(t(gcsi.aac[, row.names(gcsi.cell.info)])) 
df <- df[!is.na(df$value), ]
df$tissueid <- gcsi.cell.info[df$Var1, ]$tissueid
png('gcsiboxplotsv4.png', units="in", width=10, height=7, res=600)
df$Tumor <- ifelse(df$tissueid %in% c(c("Lymphoid", "Myeloid", "Other")), yes = "Non-solid", no = "Solid")
ggplot(df, aes(x = Var2, y = value, fill = Tumor)) + geom_boxplot() + ylab("AAC") + xlab("Drugs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  stat_compare_means(aes(group = Tumor), label = "p.signif")
dev.off()

CTRPv2 <- readRDS("~/PSet/CTRPv2.rds")
ctrp.aac <- PharmacoGx::summarizeSensitivityProfiles(
  pSet=CTRPv2,
  drugs = drugs,
  sensitivity.measure='aac_recomputed', 
  fill.missing = TRUE,
  summary.stat="median",
  verbose=FALSE)
ctrp.cell.info <- PharmacoGx::cellInfo(CTRPv2)
df <- reshape2::melt(t(ctrp.aac[,row.names(ctrp.cell.info)])) 
df <- df[!is.na(df$value), ]
df$tissueid <- ctrp.cell.info[df$Var1, ]$tissueid
df$Tumor <- ifelse(df$tissueid %in% c(c("Lymphoid", "Myeloid", "Other")), yes = "Non-solid", no = "Solid")
png('ctrpboxplotsv4.png', units="in", width=10, height=7, res=600)
ggplot(df, aes(x = Var2, y = value, fill = Tumor)) + geom_boxplot() + ylab("AAC") + xlab("Drugs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  stat_compare_means(aes(group = Tumor), label = "p.signif")
dev.off()

GDSCv2 <- readRDS("~/PSet/GDSC2.rds")
gdsc.aac <- PharmacoGx::summarizeSensitivityProfiles(
  pSet=GDSCv2,
  drugs = drugs,
  sensitivity.measure='aac_recomputed', 
  fill.missing = TRUE,
  summary.stat="median",
  verbose=FALSE)
gdsc.cell.info <- PharmacoGx::cellInfo(GDSCv2)
df <- reshape2::melt(t(gdsc.aac[, row.names(gdsc.cell.info)])) 
df <- df[!is.na(df$value), ]
df$tissueid <- gdsc.cell.info[df$Var1, ]$tissueid
df$Tumor <- ifelse(df$tissueid %in% c(c("Lymphoid", "Myeloid", "Other")), yes = "Non-solid", no = "Solid")
png('gdscboxplotsv4.png', units="in", width=10, height=7, res=600)
ggplot(df, aes(x = Var2, y = value, fill = Tumor)) + geom_boxplot() + ylab("AAC") + xlab("Drugs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  stat_compare_means(aes(group = Tumor), label = "p.signif")
dev.off()


gCSI <- readRDS("~/PSet/gCSI.rds")
for (i in 1:length(gCSI@sensitivity$profiles$ic50_recomputed)) {
  t1 <- gCSI@sensitivity$profiles$ic50_recomputed[i]
  t2 <- gCSI@sensitivity$info$chosen.max.range[i]
  t3 <- gCSI@sensitivity$info$chosen.min.range[i]
  if (!is.na(t1)) {
    if (t1>t2) {
      gCSI@sensitivity$profiles$ic50_recomputed[i] <- t2
    }
    if (t1<t3){
      gCSI@sensitivity$profiles$ic50_recomputed[i] <- t3
    }
  }
}  
gcsi.ic50 <- PharmacoGx::summarizeSensitivityProfiles(
  pSet=gCSI,
  drugs = drugs,
  sensitivity.measure='ic50_recomputed', 
  fill.missing = TRUE,
  summary.stat="median",
  verbose=FALSE)
gcsi.ic50 = log2(gcsi.ic50)

gcsi.cell.info <- PharmacoGx::cellInfo(gCSI)
df <- reshape2::melt(t(gcsi.ic50[, row.names(gcsi.cell.info)])) 
df <- df[!is.na(df$value), ]
df$tissueid <- gcsi.cell.info[df$Var1, ]$tissueid
png('gcsiboxplotsv3ic50.png', units="in", width=10, height=7, res=600, compression = 'lzw')
df$Tumor <- ifelse(df$tissueid %in% c(c("Lymphoid", "Myeloid", "Other")), yes = "Non-solid", no = "Solid")
ggplot(df, aes(x = Var2, y = value, fill = Tumor)) + geom_boxplot() + ylab("ic50") + xlab("Drugs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  stat_compare_means(aes(group = Tumor), label = "p.signif")
dev.off()

CTRPv2 <- readRDS("~/PSet/CTRPv2.rds")
for (i in 1:length(CTRPv2@sensitivity$profiles$ic50_recomputed)) {
  t1 <- CTRPv2@sensitivity$profiles$ic50_recomputed[i]
  t2 <- CTRPv2@sensitivity$info$chosen.max.range[i]
  t3 <- CTRPv2@sensitivity$info$chosen.min.range[i]
  if (!is.na(t1)) {
    if (t1>t2) {
      CTRPv2@sensitivity$profiles$ic50_recomputed[i] <- t2
    }
    if (t1<t3){
      CTRPv2@sensitivity$profiles$ic50_recomputed[i] <- t3
    }
  }
}
ctrp.ic50 <- PharmacoGx::summarizeSensitivityProfiles(
  pSet=CTRPv2,
  drugs = drugs,
  sensitivity.measure='ic50_recomputed', 
  fill.missing = TRUE,
  summary.stat="median",
  verbose=FALSE)
ctrp.ic50 = log2(ctrp.ic50)

ctrp.cell.info <- PharmacoGx::cellInfo(CTRPv2)
df <- reshape2::melt(t(ctrp.ic50[,row.names(ctrp.cell.info)])) 
df <- df[!is.na(df$value), ]
df$tissueid <- ctrp.cell.info[df$Var1, ]$tissueid
df$Tumor <- ifelse(df$tissueid %in% c(c("Lymphoid", "Myeloid", "Other")), yes = "Non-solid", no = "Solid")
png('ctrpboxplotsv3ic50.png', units="in", width=10, height=7, res=600, compression = 'lzw')
ggplot(df, aes(x = Var2, y = value, fill = Tumor)) + geom_boxplot() + ylab("ic50") + xlab("Drugs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  stat_compare_means(aes(group = Tumor), label = "p.signif")
dev.off()

GDSCv2 <- readRDS("~/PSet/GDSC2.rds")
for (i in 1:length(GDSCv2@sensitivity$profiles$ic50_recomputed)) {
  t1 <- GDSCv2@sensitivity$profiles$ic50_recomputed[i]
  t2 <- GDSCv2@sensitivity$info$MAX.CONC[i]
  t3 <- GDSCv2@sensitivity$info$MIN.CONC[i]
  if (!is.na(t1)) {
    if (t1>t2) {
      GDSCv2@sensitivity$profiles$ic50_recomputed[i] <- t2
    }
    if (t1<t3){
      GDSCv2@sensitivity$profiles$ic50_recomputed[i] <- t3
    }
  }
}
gdsc.ic50 <- PharmacoGx::summarizeSensitivityProfiles(
  pSet=GDSCv2,
  drugs = drugs,
  sensitivity.measure='ic50_recomputed', 
  fill.missing = TRUE,
  summary.stat="median",
  verbose=FALSE)
gdsc.ic50 = log2(gdsc.ic50)

gdsc.cell.info <- PharmacoGx::cellInfo(GDSCv2)
df <- reshape2::melt(t(gdsc.ic50[, row.names(gdsc.cell.info)])) 
df <- df[!is.na(df$value), ]
df$tissueid <- gdsc.cell.info[df$Var1, ]$tissueid
df$Tumor <- ifelse(df$tissueid %in% c("Lymphoid", "Myeloid", "Other"), yes = "Non-solid", no = "Solid")
png('gdscboxplotsv3ic50.png', units="in", width=10, height=7, res=600, compression = 'lzw')
ggplot(df, aes(x = Var2, y = value, fill = Tumor)) + geom_boxplot() + ylab("ic50") + xlab("Drugs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  stat_compare_means(aes(group = Tumor), label = "p.signif")
dev.off()

