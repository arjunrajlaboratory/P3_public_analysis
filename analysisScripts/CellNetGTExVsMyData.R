#!/usr/bin/env Rscript
# don't run lines 3-13 if running from R workspace
args = commandArgs(trailingOnly = T)

# check if there are exactly 3 arguments. If not, return an error
if (length(args) < 3 | length(args) > 3) {
  stop("Exactly three arguments must be supplied: projectDir, procDataSubdir and graphSubdir.", call.=FALSE)
} 
if (length(args)==3){
  projectDir = args[1]
  procDataSubdir = args[2]
  graphSubdir = args[3]
}

library(tidyverse)
library(DESeq2)
# library(SummarizedExperiment)
library(ggrepel)

# projectDir = "~/Dropbox (RajLab)/Projects/cellid/"
setwd(projectDir)
procDataDir = paste0(projectDir, procDataSubdir)
graphDir = paste0(projectDir, graphSubdir)

source("analysisScripts/RNAtagUtilities.R")
tfGeneIDfile = "annotations/TF_gene_ids.csv"
geneIDfile = "annotations/hg19gene_idToGeneSymbol.tsv"

tfTab = as_tibble(read.csv(tfGeneIDfile))
colnames(tfTab) = "gene_id"
geneIDtoGeneSymbol <- as_tibble(read.csv(geneIDfile, sep='\t'))

markerFile = "miscellaneous/markerGeneList.csv"

markerTab = as_tibble(read.csv(markerFile, header = T)) %>%
  left_join(geneIDtoGeneSymbol)

## Load data
# # CellNet
# tempSEpath = paste0(projectDir, "procDataScripted/CellNet/2017Apr05TrainingSet/TPMs")
# cellNetTPMsSE = loadHDF5SummarizedExperiment(tempSEpath)

# GTEx - not yet
tempSEpath = paste0(procDataDir, "/GTEx/V7/TPMs/V7_TPMs.rds")
GTExTPMsSE = readRDS(tempSEpath)

# # my data
# tempSEpath = paste0(projectDir, "processedData/RNAtagSeq/fibsAndiCardsTest4and5")
# tagTPMsSE = loadHDF5SummarizedExperiment(tempSEpath)


## Filter and format data
# CellNet
# sampleData = as.data.frame(colData(cellNetTPMsSE))
# sampleData$combID = rownames(sampleData)
# rownames(sampleData) = NULL
# 
# geneData = as.data.frame(rowData(cellNetTPMsSE)) #already includes gene_id column
# rownames(geneData) = geneData$GeneSymbol
# 
# datWide = as.data.frame(assay(cellNetTPMsSE))
# nSamp = dim(datWide)[2]
# 
# datWide$GeneSymbol = rownames(datWide)
# rownames(datWide) = NULL
# datTall = as_tibble(datWide) %>%
#   gather(combID, TPM, 1:nSamp)
# 
# cellNetTall = left_join(datTall, as_tibble(sampleData), by = "combID") %>%
#   left_join(as_tibble(geneData), by = "GeneSymbol")
# rm(datTall); rm(datWide)
# 
# cellNetSumStats = cellNetTall %>%
#   group_by(GeneSymbol, description1) %>%
#   summarise(meanTPM = mean(TPM), sdTPM = sd(TPM), CV = sdTPM/meanTPM, nSamp = length(TPM), nZero = sum(TPM == 0)) %>%
#   mutate(dataSource = "CellNetTraining")



# GTEx compact
sampleData = as.data.frame(colData(GTExTPMsSE))
sampleData$combID = rownames(sampleData)
rownames(sampleData) = NULL

geneData = as.data.frame(rowData(GTExTPMsSE)) #already includes gene_id column
rownames(geneData) = geneData$gene_id

dat = as.data.frame(assay(GTExTPMsSE))
nSamp = dim(dat)[2]

SMTSlist = unique(sampleData$SMTS)
SMTSDlist = unique(sampleData$SMTSD)

GTExSmtsSumStats = list()
n = 1
for (tis in SMTSlist) {
  
  cat(paste0("Working on ", tis, " (", as.character(n), "/", as.character(length(SMTSlist)), ")...\n"))
  tisSampDat = sampleData %>%
    filter(SMTS == tis)
  tempNsamp = dim(tisSampDat)[1]
  tempDat = dat[,(colnames(dat) %in% tisSampDat$combID)] # only columns for this tissue
  
  gene_ids = rownames(tempDat)
  tempMean = rowMeans(tempDat)
  tempSD = apply(tempDat,1, sd, na.rm = TRUE)
  temp0s = rowSums(tempDat == 0)
  
  tempStats = as_tibble(data.frame(
    gene_id = gene_ids,
    meanTPM = tempMean,
    sdTPM = tempSD,
    nSamp = tempNsamp,
    nZero = temp0s,
    SMTS = tis,
    dataSource = "GTEx-SMTS"
  )) %>%
    mutate(CV = sdTPM/meanTPM)
  
  if (is.null(dim(GTExSmtsSumStats))) {
    GTExSmtsSumStats = tempStats
  } else {
    GTExSmtsSumStats = bind_rows(GTExSmtsSumStats, tempStats)
  }
  
  cat(paste0("Completed ", tis, ".\n"))
  n = n+1
}
write.table(GTExSmtsSumStats, file = paste0(procDataDir,"/GTEx/V7/SMTS_summaryStats.txt"), sep = "\t", quote = F, row.names = F)

GTExSmtsdSumStats = list()
n = 1
for (tis in SMTSDlist) {
  
  cat(paste0("Working on ", tis, " (", as.character(n), "/", as.character(length(SMTSDlist)), ")...\n"))
  tisSampDat = sampleData %>%
    filter(SMTSD == tis)
  tempNsamp = dim(tisSampDat)[1]
  tempDat = dat[,(colnames(dat) %in% tisSampDat$combID)] # only columns for this tissue
  
  gene_ids = rownames(tempDat)
  tempMean = rowMeans(tempDat)
  tempSD = apply(tempDat,1, sd, na.rm = TRUE)
  temp0s = rowSums(tempDat == 0)
  
  tempStats = as_tibble(data.frame(
    gene_id = gene_ids,
    meanTPM = tempMean,
    sdTPM = tempSD,
    nSamp = tempNsamp,
    nZero = temp0s,
    SMTSD = tis,
    dataSource = "GTEx-SMTSD"
  )) %>%
    mutate(CV = sdTPM/meanTPM)
  
  if (is.null(dim(GTExSmtsdSumStats))) {
    GTExSmtsdSumStats = tempStats
  } else {
    GTExSmtsdSumStats = bind_rows(GTExSmtsdSumStats, tempStats)
  }
  
  cat(paste0("Completed ", tis, ".\n"))
  n = n+1
}
write.table(GTExSmtsdSumStats, file = paste0(procDataDir,"/GTEx/V7/SMTSD_summaryStats.txt"), sep = "\t", quote = F, row.names = F)



# # my data
# tagTall = SummarizedExperimentToTallTibble(tagTPMsSE)
# 
# tagSumStats = tagTall %>%
#   filter(Product.Name == "control", cellType != "HeLa") %>%
#   group_by(gene_id, GeneSymbol, cellType) %>%
#   summarise(meanTPM = mean(TPM), sdTPM = sd(TPM), CV = sdTPM/meanTPM, nSamp = length(TPM), nZero = sum(TPM == 0)) %>%
#   mutate(dataSource = "RNAtagSeq_IAM")
# 
# read.table("procDataScripted/GTEx/V7/SMTSD_summaryStats.txt", sep = "\t", stringsAsFactors = F)
# 
# 
# 
# cellNetSumStatsMarkers = left_join(markerTab, cellNetSumStats, by = "GeneSymbol")
# tagSumStatsMarkers = left_join(markerTab, tagSumStats, by = c("GeneSymbol", "gene_id"))
# gtexSmtsSumStatsMarkers = left_join(markerTab, GTExSmtsSumStats, by = "gene_id")
# 
# allSumMarks = bind_rows(cellNetSumStatsMarkers, tagSumStatsMarkers) %>%
#   mutate(cellTypeJoined = ifelse(dataSource == "RNAtagSeq_IAM", cellType, description1)) %>%
#   filter(meanTPM > 5) %>%
#   filter(cellTypeMark == "iCard") %>%
#   unique()
# 
# allSumMarks = bind_rows(gtexSmtsSumStatsMarkers, tagSumStatsMarkers) %>%
#   mutate(cellTypeJoined = ifelse(dataSource == "RNAtagSeq_IAM", cellType, SMTS)) %>%
#   filter(meanTPM > 5) %>%
#   filter(cellTypeMark == "iCard") %>%
#   unique()
# 
# 
# markerDots = ggplot() +
#   geom_point(data = allSumMarks, aes(log(meanTPM), CV, color = dataSource), size = 2) +
#   geom_text_repel(data = allSumMarks, aes(log(meanTPM), CV, label = cellTypeJoined)) + 
#   facet_wrap(~GeneSymbol) + 
#   theme_bw() +
#   ggtitle("iCard Marker Genes: Mean TPM vs CV(TPM)\nRed Dots: CellNet Training data (several cell types)\nBlue dots: cellid RNAtag-seq CONTROL iCards and fibroblasts")
# plot(markerDots)
# ggsave(markerDots, file = "graphs/CellNet/CellNetVsTests45controls_iCardMarkers_meanVsCV.pdf", width = 16, height = 16)
# 
# 
# markerDots = ggplot() +
#   geom_point(data = allSumMarks, aes(log(meanTPM), CV, color = dataSource), size = 2) +
#   geom_text_repel(data = allSumMarks, aes(log(meanTPM), CV, label = cellTypeJoined)) + 
#   facet_wrap(~GeneSymbol) + 
#   theme_bw() +
#   ggtitle("iCard Marker Genes: Mean TPM vs CV(TPM)\nRed Dots: GTEx v7 data (several tissue types)\nBlue dots: cellid RNAtag-seq CONTROL iCards and fibroblasts")
# plot(markerDots)
# ggsave(markerDots, file = "graphs/GTEx/GTExSMTSvsTests45controls_iCardMarkers_meanVsCV.pdf", width = 30, height = 16)

