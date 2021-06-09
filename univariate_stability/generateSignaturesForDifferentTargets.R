library(SummarizedExperiment)
library(PharmacoGx)
# memory.limit(size=12000)
drugs <- c("Bortezomib", "Entinostat", "Sirolimus","Docetaxel","Gemcitabine", "Crizotinib",  "Lapatinib","Vorinostat","Erlotinib","Paclitaxel","Pictilisib")
gCSI <- readRDS("PSet/gCSI.rds")

system.time(gCSI.sig.aac <- drugSensitivitySig(gCSI,mDataType="Kallisto_0.46.1.rnaseq", 
  sensitivity.measure="aac_recomputed", nthread=23, drugs=drugs) )

gCSI@sensitivity$profiles$ic50_trunc<- gCSI@sensitivity$profiles$ic50_recomputed
for (i in 1:length(gCSI@sensitivity$profiles$ic50_recomputed)) {
  t1 <- gCSI@sensitivity$profiles$ic50_recomputed[i]
  t2 <- gCSI@sensitivity$info$chosen.max.range[i]
  t3 <- gCSI@sensitivity$info$chosen.min.range[i]
  if (!is.na(t1)) {
    if (t1>t2) {
      gCSI@sensitivity$profiles$ic50_trunc[i] <- t2
    }
    if (t1<t3){
      gCSI@sensitivity$profiles$ic50_trunc[i] <- t3
    }
  }
}  

gCSI@sensitivity$profiles$ic50_log <- log10(gCSI@sensitivity$profiles$ic50_recomputed)

gCSI@sensitivity$profiles$ic50_logtunc <- log10(gCSI@sensitivity$profiles$ic50_tunc)

system.time(gCSI.sig.ic50_log <- drugSensitivitySig(gCSI,mDataType="Kallisto_0.46.1.rnaseq", 
  sensitivity.measure="ic50_log", nthread=23, drugs=drugs) )

system.time(gCSI.sig.ic50_logtunc <- drugSensitivitySig(gCSI,mDataType="Kallisto_0.46.1.rnaseq", 
  sensitivity.measure="ic50_logtunc", nthread=23, drugs=drugs) )

system.time(gCSI.sig.ic50_recomputed <- drugSensitivitySig(gCSI,mDataType="Kallisto_0.46.1.rnaseq", 
  sensitivity.measure="ic50_recomputed", nthread=23, drugs=drugs) )

system.time(gCSI.sig.ic50_trunc <- drugSensitivitySig(gCSI,mDataType="Kallisto_0.46.1.rnaseq", 
  sensitivity.measure="ic50_trunc", nthread=23, drugs=drugs) )



save(gCSI.sig.ic50_log,
 gCSI.sig.ic50_recomputed,
  gCSI.sig.aac,
  gCSI.sig.ic50_logtunc,
  gCSI.sig.ic50_trunc,
   file='gCSI.sigs.RData')
rm(list=ls())


drugs <- c("Bortezomib", "Entinostat", "Sirolimus","Docetaxel","Gemcitabine", "Crizotinib",  "Lapatinib","Vorinostat","Erlotinib","Paclitaxel","Pictilisib")


GDSCv2 <- readRDS("PSet/GDSC2.rds")

GDSCv2@sensitivity$profiles$ic50_tunc<- GDSCv2@sensitivity$profiles$ic50_recomputed
for (i in 1:length(GDSCv2@sensitivity$profiles$ic50_recomputed)) {
  t1 <- GDSCv2@sensitivity$profiles$ic50_recomputed[i]
  t2 <- GDSCv2@sensitivity$info$MAX.CONC[i]
  t3 <- GDSCv2@sensitivity$info$MIN.CONC[i]
  if (!is.na(t1)) {
    if (t1>t2) {
      GDSCv2@sensitivity$profiles$ic50_trunc[i] <- t2
    }
    if (t1<t3){
      GDSCv2@sensitivity$profiles$ic50_trunc[i] <- t3
    }
  }
}


GDSCv2@sensitivity$profiles$ic50_log <- log10(GDSCv2@sensitivity$profiles$ic50_recomputed)

GDSCv2@sensitivity$profiles$ic50_logtunc <- log10(GDSCv2@sensitivity$profiles$ic50_tunc)


system.time(GDSCv2.sig.aac <- drugSensitivitySig(GDSCv2,mDataType="Kallisto_0.46.1.rnaseq", 
  sensitivity.measure="aac_recomputed", nthread=23, drugs=drugs) )


system.time(GDSCv2.sig.ic50_log <- drugSensitivitySig(GDSCv2,mDataType="Kallisto_0.46.1.rnaseq", 
  sensitivity.measure="ic50_log", nthread=23, drugs=drugs) )

system.time(GDSCv2.sig.ic50_logtunc <- drugSensitivitySig(GDSCv2,mDataType="Kallisto_0.46.1.rnaseq", 
  sensitivity.measure="ic50_logtunc", nthread=23, drugs=drugs) )

system.time(GDSCv2.sig.ic50_recomputed <- drugSensitivitySig(GDSCv2,mDataType="Kallisto_0.46.1.rnaseq", 
  sensitivity.measure="ic50_recomputed", nthread=23, drugs=drugs) )


system.time(GDSCv2.sig.ic50_trunc <- drugSensitivitySig(GDSCv2,mDataType="Kallisto_0.46.1.rnaseq", 
  sensitivity.measure="ic50_trunc", nthread=23, drugs=drugs) )



save(GDSCv2.sig.ic50_log,
 GDSCv2.sig.ic50_recomputed,
  GDSCv2.sig.aac,
  GDSCv2.sig.ic50_logtunc,
  GDSCv2.sig.ic50_trunc,
   file='GDSCv2.sigs.RData')
rm(list=ls())


drugs <- c("Bortezomib", "Entinostat", "Sirolimus","Docetaxel","Gemcitabine", "Crizotinib",  "Lapatinib","Vorinostat","Erlotinib","Paclitaxel","Pictilisib")




CTRPv2 <- readRDS("PSet/CTRPv2.rds")

CCLE <- readRDS("PSet/CCLE.rds")

source("mergePSets.R")

CTRPv2 <- mergePSets(CCLE, CTRPv2)


CTRPv2@sensitivity$profiles$ic50_tunc<- CTRPv2@sensitivity$profiles$ic50_recomputed
for (i in 1:length(CTRPv2@sensitivity$profiles$ic50_recomputed)) {
  t1 <- CTRPv2@sensitivity$profiles$ic50_recomputed[i]
  t2 <- CTRPv2@sensitivity$info$chosen.max.range[i]
  t3 <- CTRPv2@sensitivity$info$chosen.min.range[i]
  if (!is.na(t1)) {
    if (t1>t2) {
      CTRPv2@sensitivity$profiles$ic50_trunc[i] <- t2
    }
    if (t1<t3){
      CTRPv2@sensitivity$profiles$ic50_trunc[i] <- t3
    }
  }
}


CTRPv2@sensitivity$profiles$ic50_log <- log10(CTRPv2@sensitivity$profiles$ic50_recomputed)

CTRPv2@sensitivity$profiles$ic50_logtunc <- log10(CTRPv2@sensitivity$profiles$ic50_tunc)


system.time(CTRPv2.sig.aac <- drugSensitivitySig(CTRPv2,mDataType="Kallisto_0.46.1.rnaseq", 
  sensitivity.measure="aac_recomputed", nthread=23, drugs=drugs) )


system.time(CTRPv2.sig.ic50_log <- drugSensitivitySig(CTRPv2,mDataType="Kallisto_0.46.1.rnaseq", 
  sensitivity.measure="ic50_log", nthread=23, drugs=drugs) )

system.time(CTRPv2.sig.ic50_logtunc <- drugSensitivitySig(CTRPv2,mDataType="Kallisto_0.46.1.rnaseq", 
  sensitivity.measure="ic50_logtunc", nthread=23, drugs=drugs) )

system.time(CTRPv2.sig.ic50_recomputed <- drugSensitivitySig(CTRPv2,mDataType="Kallisto_0.46.1.rnaseq", 
  sensitivity.measure="ic50_recomputed", nthread=23, drugs=drugs) )

system.time(CTRPv2.sig.ic50_trunc <- drugSensitivitySig(CTRPv2,mDataType="Kallisto_0.46.1.rnaseq", 
  sensitivity.measure="ic50_trunc", nthread=23, drugs=drugs) )


save(CTRPv2.sig.ic50_log,
 CTRPv2.sig.ic50_recomputed,
  CTRPv2.sig.aac,
  CTRPv2.sig.ic50_logtunc,
  CTRPv2.sig.ic50_trunc,
   file='CTRPv2.sigs.RData')


table(apply(GDSC1@sensitivity$raw[,,1], 1, function(x) return(sum(!is.na(x))) ))
table(apply(GDSCv2@sensitivity$raw[,,1], 1, function(x) return(sum(!is.na(x))) ))
table(apply(gCSI@sensitivity$raw[,,1], 1, function(x) return(sum(!is.na(x))) ))
table(apply(CTRPv2@sensitivity$raw[,,1], 1, function(x) return(sum(!is.na(x))) ))


