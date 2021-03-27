library(SummarizedExperiment)
library(PharmacoGx)

memory.limit(size=12000)

drugs <- c("Bortezomib", "Entinostat", "Sirolimus","Docetaxel","Gemcitabine", "Crizotinib",  "Lapatinib","Vorinostat","Erlotinib","Paclitaxel","Pictilisib")

gCSI <- readRDS("~/PSet/gCSI.rds")
non_haema2 <- gCSI@cell[-which(gCSI@cell$tissueid == "Lymphoid" |
                                 gCSI@cell$tissueid == "Myeloid" | 
                                 gCSI@cell$tissueid == "Other"),]

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
  cell.lines = row.names(non_haema2),
  sensitivity.measure='ic50_recomputed', 
  fill.missing = TRUE,
  summary.stat="median",
  verbose=FALSE)
rm(gCSI)
log.gcsi.ic50 = log2(gcsi.ic50)

GDSCv2 <- readRDS("~/PSet/GDSC2.rds")
non_haema1 <- GDSCv2@cell[-which(GDSCv2@cell$tissueid == "Lymphoid"|
                                   GDSCv2@cell$tissueid == "Myeloid"|
                                   GDSCv2@cell$tissueid == "Other"),]
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
  #cell.lines =  row.names(non_haema1),
  sensitivity.measure='ic50_recomputed', 
  fill.missing = TRUE,
  summary.stat="median",
  verbose=FALSE)
rm(GDSCv2)
log.gdsc.ic50 = log2(gdsc.ic50)

CTRPv2 <- readRDS("~/PSet/CTRPv2.rds")
non_haema4 <- CTRPv2@cell[-which(CTRPv2@cell$tissueid == "Lymphoid" |
                                   CTRPv2@cell$tissueid == "Myeloid" | 
                                   CTRPv2@cell$tissueid == "Other"),]
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
  #cell.lines = row.names(non_haema4),
  sensitivity.measure='ic50_recomputed', 
  fill.missing = TRUE,
  summary.stat="median",
  verbose=FALSE)
rm(CTRPv2)
log.ctrp.ic50 = log2(ctrp.ic50)

GDSC1 <- readRDS("~/PSet/GDSC1.rds")
non_haema5 <- GDSC1@cell[-which(GDSC1@cell$tissueid == "Lymphoid" |
                                  GDSC1@cell$tissueid == "Myeloid"|
                                  GDSC1@cell$tissueid == "Other"),]
for (i in 1:length(GDSC1@sensitivity$profiles$ic50_recomputed)) {
  t1 <- GDSC1@sensitivity$profiles$ic50_recomputed[i]
  t2 <- GDSC1@sensitivity$info$chosen.max.range[i]
  t3 <- GDSC1@sensitivity$info$chosen.min.range[i]
  if (!is.na(t1)) {
    if (t1>t2) {
      GDSC1@sensitivity$profiles$ic50_recomputed[i] <- t2
    }
    if (t1<t3){
      GDSC1@sensitivity$profiles$ic50_recomputed[i] <- t3
    }
  }
}

gdscv1.ic50 <- PharmacoGx::summarizeSensitivityProfiles(
  pSet=GDSC1,
  drugs = drugs,
  #cell.lines = row.names(non_haema5),
  sensitivity.measure='ic50_recomputed', 
  fill.missing = TRUE,
  summary.stat="median",
  verbose=FALSE)
rm(GDSC1)

log.gdscv1.ic50 = log2(gdscv1.ic50)

write.table(log.ctrp.ic50, file = "CTRP.logIC50.tsv", row.names=TRUE, sep="\t", col.names = TRUE)
write.table(log.gdsc.ic50, file = "GDSC.logIC50.tsv", row.names=TRUE, sep="\t", col.names = TRUE)
write.table(log.gcsi.ic50, file = "gCSI.logIC50.tsv", row.names=TRUE, sep="\t", col.names = TRUE)
write.table(log.gdscv1.ic50, file = "GDSCv1.logIC50.tsv", row.names=TRUE, sep="\t", col.names = TRUE)

write.table(log.ctrp.ic50, file = "CTRP.logIC50All.tsv", row.names=TRUE, sep="\t", col.names = TRUE)
write.table(log.gdsc.ic50, file = "GDSC.logIC50All.tsv", row.names=TRUE, sep="\t", col.names = TRUE)
write.table(log.gcsi.ic50, file = "gCSI.logIC50All.tsv", row.names=TRUE, sep="\t", col.names = TRUE)
write.table(log.gdscv1.ic50, file = "GDSCv1.logIC50All.tsv", row.names=TRUE, sep="\t", col.names = TRUE)
