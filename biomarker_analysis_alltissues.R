# biomarker analysis - all tissues included
library(PharmacoGx); library(SummarizedExperiment); library(knitr); library(readxl); library(readr)

# loading PSets
CTRPv2 <- readRDS("../code/PSets/CTRPv2.rds")
GDSC2 <- readRDS("../code/PSets/GDSC2.rds")
CCLE <- readRDS("../code/PSets/CCLE.rds")
gCSI <- readRDS("../code/PSets/gCSI.rds")
GDSC1 <- readRDS("../code/PSets/GDSC1.rds")

# intersect PSets by common drugs
common <- PharmacoGx::intersectPSet(list('gCSI' = gCSI, 'GDSCv2' = GDSC2, 'CTRPv2' = CTRPv2, 'GDSCv1' = GDSC1), intersectOn = c("drugs"))
ctrp <- common[["CTRPv2"]] # no molprof info available
gdsc2 <- common[["GDSCv2"]] 
gcsi <- common[["gCSI"]]
gdsc1 <- common[["GDSCv1"]]
rm(common)

common_drugs <- rownames(ctrp@drug)

#use uniprot id from biomarker list to get ensembl genes
uniprot <- readxl::read_excel("uniprot-id-list.xlsx")
genes <- as.data.frame(uniprot$`Ensembl Gene`, row.names = uniprot$`Ensembl Gene`)
colnames(genes) <- "genes"

genes["ENSG00000138413", ] <- "ENSG00000138413.13"
genes["ENSG00000169083", ] <- "ENSG00000169083.17"
genes["ENSG00000146648", ] <- "ENSG00000146648.18"
genes["ENSG00000140443", ] <- "ENSG00000140443.15"
genes["ENSG00000139618", ] <- "ENSG00000139618.15"
genes["ENSG00000136244", ] <- "ENSG00000136244.12"
genes["ENSG00000171791", ] <- "ENSG00000171791.13"
genes["ENSG00000017427", ] <- "ENSG00000017427.16"
genes["ENSG00000118705", ] <- "ENSG00000118705.17"
genes["ENSG00000125968", ] <- "ENSG00000125968.9"
genes["ENSG00000121879", ] <- "ENSG00000121879.5"
genes["ENSG00000176890", ] <- "ENSG00000176890.16"
genes["ENSG00000156136", ] <- "ENSG00000156136.10"
genes["ENSG00000101966", ] <- "ENSG00000101966.12"
genes["ENSG00000167325", ] <- "ENSG00000167325.15"
genes["ENSG00000181143", ] <- "ENSG00000181143.15"
genes["ENSG00000112759", ] <- "ENSG00000112759.19"
genes["ENSG00000204673", ] <- "ENSG00000204673.10"
genes["ENSG00000155066", ] <- "ENSG00000155066.16"
genes["ENSG00000133703", ] <- "ENSG00000133703.12"
genes["ENSG00000002330", ] <- "ENSG00000002330.13"
genes["ENSG00000119318", ] <- "ENSG00000119318.13"
genes["ENSG00000168610", ] <- "ENSG00000168610.14"
genes["ENSG00000111087", ] <- "ENSG00000111087.10"
genes["ENSG00000258947", ] <- "ENSG00000258947.7"
genes["ENSG00000167085", ] <- "ENSG00000167085.11"
genes["ENSG00000153944", ] <- "ENSG00000153944.11"
genes["ENSG00000171094", ] <- "ENSG00000171094.18"
genes["ENSG00000047936", ] <- "ENSG00000047936.10"
genes["ENSG00000141736", ] <- "ENSG00000141736.13"
genes["ENSG00000105976", ] <- "ENSG00000105976.15"
genes["ENSG00000065361", ] <- "ENSG00000065361.15"
genes["ENSG00000140009", ] <- "ENSG00000140009.18"
genes["ENSG00000171862", ] <- "ENSG00000171862.11"
genes["ENSG00000069869", ] <- "ENSG00000069869.16"
genes["ENSG00000101856", ] <- "ENSG00000101856.10"
genes["ENSG00000141646", ] <- "ENSG00000141646.14"
genes["ENSG00000047230", ] <- "ENSG00000047230.15"
genes["ENSG00000141510", ] <- "ENSG00000141510.17"
genes["ENSG00000170017", ] <- "ENSG00000170017.12"
genes["ENSG00000124762", ] <- "ENSG00000124762.13"
genes["ENSG00000044574", ] <- "ENSG00000044574.8"
genes["ENSG00000114423", ] <- "ENSG00000114423.22"
genes["ENSG00000097007", ] <- "ENSG00000097007.19"
genes["ENSG00000197122", ] <- "ENSG00000197122.11"
genes["ENSG00000182054", ] <- "ENSG00000182054.9"
genes["ENSG00000162733", ] <- "ENSG00000162733.19"
genes["ENSG00000157404", ] <- "ENSG00000157404.15"
genes["ENSG00000134853", ] <- "ENSG00000134853.12"
genes["ENSG00000119535", ] <- "ENSG00000119535.17"
genes["ENSG00000110395", ] <- "ENSG00000110395.7"
genes["ENSG00000213281", ] <- "ENSG00000213281.5"
genes <- genes$genes

# mapping between CCLE and CTRPv2 cell lines
ccle_celllines <- CCLE@molecularProfiles[["Kallisto_0.46.1.rnaseq"]]@colData@listData[["Cell_Line"]]
ctrp_celllines <- ctrp@cell$cellid

ccle_molprof <- as.data.frame(CCLE@molecularProfiles[["Kallisto_0.46.1.rnaseq"]]@colData@listData)
common_celllines <- subset(ccle_molprof, ctrp_celllines %in% ccle_celllines) # rnaseq info for cell lines tested in CTRPv2
rownames(common_celllines) <- common_celllines$cellid

ctrp@molecularProfiles[["Kallisto_0.46.1.rnaseq"]] <- CCLE@molecularProfiles[["Kallisto_0.46.1.rnaseq"]]

ctrp.aac.sigs <- PharmacoGx::drugSensitivitySig(ctrp, "Kallisto_0.46.1.rnaseq", sensitivity.measure = "aac_recomputed", drugs=common_drugs, features = genes)
ctrp.ic50.sigs <- PharmacoGx::drugSensitivitySig(ctrp, "Kallisto_0.46.1.rnaseq", sensitivity.measure = "ic50_recomputed", drugs=common_drugs, features = genes)


#GDSC2
gdsc2.aac.sigs <- PharmacoGx::drugSensitivitySig(GDSC2, "Kallisto_0.46.1.rnaseq", sensitivity.measure = "aac_recomputed", drugs= common_drugs, features = genes)
gdsc2.ic50.sigs <- PharmacoGx::drugSensitivitySig(GDSC2, "Kallisto_0.46.1.rnaseq", sensitivity.measure = "ic50_recomputed", drugs= common_drugs,  features = genes)

#gCSI
gcsi.aac.sigs <- PharmacoGx::drugSensitivitySig(gCSI, "Kallisto_0.46.1.rnaseq", sensitivity.measure = "aac_recomputed", drugs= common_drugs, features = genes)
gcsi.ic50.sigs <- PharmacoGx::drugSensitivitySig(gCSI, "Kallisto_0.46.1.rnaseq", sensitivity.measure = "ic50_recomputed", drugs= common_drugs, features = genes)

Top_biomarkers_Farnoosh_Sheet1 <- read_csv("Top biomarkers_Farnoosh - Sheet1.csv")
gene_map = as.data.frame(Top_biomarkers_Farnoosh_Sheet1[, 1:3])
gene_map[31, "Uniprot ID"] <- "P08581" #missing value in csv
gene_map$ensembl <- NA
uniprot_id_list <- read_excel("uniprot-id-list.xlsx")
uniprot_id_list <- as.data.frame(uniprot_id_list)
rownames(uniprot_id_list) <- uniprot_id_list$Entry
for (row in rownames(gene_map)) {
  x <- gene_map[row, "Uniprot ID"]
  gene_map[row, "ensembl"] <- uniprot_id_list[x, "Ensembl Gene"]
}
gene_map[63, "ensembl"] <- 	"ENSG00000157764"
gene_map[66, "ensembl"] <- 	"ENSG00000157764"
gene_map$ensemblvar <- NA
gene_map[which(gene_map$ensembl == "ENSG00000138413"), "ensemblvar"] <- "ENSG00000138413.13"
gene_map[which(gene_map$ensembl == "ENSG00000169083"),"ensemblvar" ] <- "ENSG00000169083.17"
gene_map[which(gene_map$ensembl == "ENSG00000146648"), "ensemblvar" ] <- "ENSG00000146648.18"
gene_map[which(gene_map$ensembl == "ENSG00000140443"), "ensemblvar" ] <- "ENSG00000140443.15"
gene_map[which(gene_map$ensembl == "ENSG00000139618"), "ensemblvar" ] <- "ENSG00000139618.15"
gene_map[which(gene_map$ensembl == "ENSG00000136244"), "ensemblvar" ] <- "ENSG00000136244.12"
gene_map[which(gene_map$ensembl == "ENSG00000171791"), "ensemblvar" ] <- "ENSG00000171791.13"
gene_map[which(gene_map$ensembl == "ENSG00000017427"), "ensemblvar" ] <- "ENSG00000017427.16"
gene_map[which(gene_map$ensembl == "ENSG00000118705"), "ensemblvar" ] <- "ENSG00000118705.17"
gene_map[which(gene_map$ensembl == "ENSG00000125968"), "ensemblvar" ] <- "ENSG00000125968.9"
gene_map[which(gene_map$ensembl == "ENSG00000121879"), "ensemblvar" ] <- "ENSG00000121879.5"
gene_map[which(gene_map$ensembl == "ENSG00000176890"), "ensemblvar" ] <- "ENSG00000176890.16"
gene_map[which(gene_map$ensembl == "ENSG00000156136"), "ensemblvar" ] <- "ENSG00000156136.10"
gene_map[which(gene_map$ensembl == "ENSG00000101966"), "ensemblvar" ] <- "ENSG00000101966.12"
gene_map[which(gene_map$ensembl == "ENSG00000167325"), "ensemblvar" ] <- "ENSG00000167325.15"
gene_map[which(gene_map$ensembl == "ENSG00000181143"), "ensemblvar" ] <- "ENSG00000181143.15"
gene_map[which(gene_map$ensembl == "ENSG00000112759"), "ensemblvar" ] <- "ENSG00000112759.19"
gene_map[which(gene_map$ensembl == "ENSG00000204673"), "ensemblvar" ] <- "ENSG00000204673.10"
gene_map[which(gene_map$ensembl == "ENSG00000155066"), "ensemblvar" ] <- "ENSG00000155066.16"
gene_map[which(gene_map$ensembl == "ENSG00000133703"), "ensemblvar" ] <- "ENSG00000133703.12"
gene_map[which(gene_map$ensembl == "ENSG00000002330"), "ensemblvar" ] <- "ENSG00000002330.13"
gene_map[which(gene_map$ensembl == "ENSG00000119318"), "ensemblvar" ] <- "ENSG00000119318.13"
gene_map[which(gene_map$ensembl == "ENSG00000168610"), "ensemblvar" ] <- "ENSG00000168610.14"
gene_map[which(gene_map$ensembl == "ENSG00000111087"), "ensemblvar" ] <- "ENSG00000111087.10"
gene_map[which(gene_map$ensembl == "ENSG00000258947"), "ensemblvar" ] <- "ENSG00000258947.7"
gene_map[which(gene_map$ensembl == "ENSG00000167085"), "ensemblvar" ] <- "ENSG00000167085.11"
gene_map[which(gene_map$ensembl == "ENSG00000153944"), "ensemblvar" ] <- "ENSG00000153944.11"
gene_map[which(gene_map$ensembl == "ENSG00000171094"), "ensemblvar" ] <- "ENSG00000171094.18"
gene_map[which(gene_map$ensembl == "ENSG00000047936"), "ensemblvar" ] <- "ENSG00000047936.10"
gene_map[which(gene_map$ensembl == "ENSG00000141736"), "ensemblvar" ] <- "ENSG00000141736.13"
gene_map[which(gene_map$ensembl == "ENSG00000105976"), "ensemblvar" ] <- "ENSG00000105976.15"
gene_map[which(gene_map$ensembl == "ENSG00000065361"), "ensemblvar" ] <- "ENSG00000065361.15"
gene_map[which(gene_map$ensembl == "ENSG00000140009"), "ensemblvar" ] <- "ENSG00000140009.18"
gene_map[which(gene_map$ensembl == "ENSG00000171862"), "ensemblvar" ] <- "ENSG00000171862.11"
gene_map[which(gene_map$ensembl == "ENSG00000069869"), "ensemblvar" ] <- "ENSG00000069869.16"
gene_map[which(gene_map$ensembl == "ENSG00000101856"), "ensemblvar" ] <- "ENSG00000101856.10"
gene_map[which(gene_map$ensembl == "ENSG00000141646"), "ensemblvar" ] <- "ENSG00000141646.14"
gene_map[which(gene_map$ensembl == "ENSG00000047230"), "ensemblvar" ] <- "ENSG00000047230.15"
gene_map[which(gene_map$ensembl == "ENSG00000141510"), "ensemblvar" ] <- "ENSG00000141510.17"
gene_map[which(gene_map$ensembl == "ENSG00000170017"), "ensemblvar" ] <- "ENSG00000170017.12"
gene_map[which(gene_map$ensembl == "ENSG00000124762"), "ensemblvar" ] <- "ENSG00000124762.13"
gene_map[which(gene_map$ensembl == "ENSG00000044574"), "ensemblvar" ] <- "ENSG00000044574.8"
gene_map[which(gene_map$ensembl == "ENSG00000114423"), "ensemblvar" ] <- "ENSG00000114423.22"
gene_map[which(gene_map$ensembl == "ENSG00000097007"), "ensemblvar" ] <- "ENSG00000097007.19"
gene_map[which(gene_map$ensembl == "ENSG00000197122"), "ensemblvar" ] <- "ENSG00000197122.11"
gene_map[which(gene_map$ensembl == "ENSG00000182054"), "ensemblvar" ] <- "ENSG00000182054.9"
gene_map[which(gene_map$ensembl == "ENSG00000162733"), "ensemblvar" ] <- "ENSG00000162733.19"
gene_map[which(gene_map$ensembl == "ENSG00000157404"), "ensemblvar" ] <- "ENSG00000157404.15"
gene_map[which(gene_map$ensembl == "ENSG00000134853"), "ensemblvar" ] <- "ENSG00000134853.12"
gene_map[which(gene_map$ensembl == "ENSG00000119535"), "ensemblvar" ] <- "ENSG00000119535.17"
gene_map[which(gene_map$ensembl == "ENSG00000110395"), "ensemblvar" ] <- "ENSG00000110395.7"
gene_map[which(gene_map$ensembl == "ENSG00000213281"), "ensemblvar" ] <- "ENSG00000213281.5"
gene_map[which(gene_map$ensembl == "ENSG00000157764"), "ensemblvar" ] <- "ENSG00000157764.13"

WriteXLS::WriteXLS(as.data.frame(gene_map), "genes_annotations.xls")

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

for (row in rownames(lap)) {
  x <- lap[row, "ensemblvar"]
  y <- "Lapatinib.estimate"
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

biomarker_alltissues <- do.call("rbind", list(lap, doc, pic, gem, vor, pac, cri, erl))
WriteXLS::WriteXLS(biomarker_alltissues, "../results/biomarker_alltissues.xls")
