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
library(fuzzyjoin)
# library(GenomicFeatures)

# projectDir = "/Users/iamellis/Dropbox (RajLab)/Shared_IanM/cellid_201807_onward/"
setwd(projectDir)
procDataDir = paste0(projectDir, procDataSubdir)
graphDir = paste0(projectDir, graphSubdir)


# ## specify necessary files
# # take out "tinySet." when ready to do full dataset
# # GTExFile = "extractedData/RNAseq/GTEx/allCountDataV6p/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.tinySet.gct"
# GTExFile = "extractedData/RNAseq/GTEx/allCountDataV6p/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct"
# GTExSampleDataFile1 = "extractedData/RNAseq/GTEx/allCountDataV6p/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt"
# GTExSampleDataFile2 = "extractedData/RNAseq/GTEx/allCountDataV6p/GTEx_Data_V6_Annotations_SampleAttributesDS.txt"
# 
# # GENCODEv19file = "extractedData/RNAseq/GTEx/allCountDataV6p/gencode.v19.annotation.gtf"
# geneIDlengthsFile = "annotations/hg19_lengths_per_gene_id.tsv"
# geneIDtoGeneSymbolFile = "annotations/hg19gene_idToGeneSymbol.tsv"
# 
# 
# ## process sample data
# GTExSampleData1 = as_tibble(read.table(GTExSampleDataFile1, sep = "\t", header= T))
# GTExSampleData1$SUBJID = as.character(GTExSampleData1$SUBJID)
# GTExSampleData = as_tibble(read.table(GTExSampleDataFile2, sep = "\t", header= T)) 
# GTExSampleData$SAMPID = as.character(GTExSampleData$SAMPID)
# GTExSampleData = GTExSampleData %>%
#   regex_left_join(GTExSampleData1, by = c(SAMPID = "SUBJID"))
# # GTExSampleData = as.data.frame(GTExSampleData)
# # rownames(GTExSampleData) = gsub("-",".", GTExSampleData$SAMPID)
# 
# 
# ## get gene_id attributes
# geneIDtoGeneSymbol = as_tibble(read.table(geneIDtoGeneSymbolFile, sep = "\t", header = T))
# geneIDlengths = as_tibble(read.table(file = geneIDlengthsFile, sep = '\t', header = T))
# 
# 
# ## load GTEx counts
# gtexCounts = as_tibble(read.table(file = GTExFile, skip = 2, sep = "\t", header = T)) %>%
#   separate(Name, into = c("gene_id", "ver"), by = "\\.") %>%
#   dplyr::select(-ver)
# 
# 
# ## make geneData (rowData in summarizedExperiment)
# gtexGeneData = gtexCounts %>%
#   dplyr::select(gene_id, Description) %>%
#   inner_join(geneIDlengths, by = "gene_id")
# colnames(gtexGeneData) = c("gene_id", "GeneSymbol", "length")
# gtexGeneData = as.data.frame(gtexGeneData)
# rownames(gtexGeneData) = gtexGeneData$gene_id
# 
# 
# # txdb <- makeTxDbFromGFF(GENCODEv19file, format="gtf")
# # lengthsPergeneidGTEx <- sum(width(IRanges::reduce(exonsBy(txdb, by = "gene")))) # a named int(), in bases
# # lengthtbl<-as_tibble(list(
# #   gene_id = names(lengthsPergeneidGTEx),
# #   length = lengthsPergeneidGTEx
# # ))
# 
# 
# ## Prespecify sample filters to be used for downstream analysis
# # Sample with metadata (7790 in SampleAttributes.txt file, but 8555 sample columns in count matrix...)
# # Annotated "SMTSD" : a more refined sample source location ID
# # RIN >= 7 : RNA integrity
# # SMMAPRT >= 0.5 : mapping rate >= 50%
# # SMMPPD >= 1e6 : 1MM or more mapped reads
# ## with these filters, there are 3623 samples left
# 
# GTExSampleDataFilt = GTExSampleData %>%
#   filter(!is.na(SMTSD),
#          SMRIN >= 7,
#          SMMAPRT >= 0.5,
#          SMMPPD >= 1e6)
# GTExSampleDataFilt = as.data.frame(GTExSampleDataFilt)
# rownames(GTExSampleDataFilt) = gsub("-",".", as.character(GTExSampleDataFilt$SAMPID))
# 
# 
# ## filter count matrix to include only pre-filtered samples
# gtexCounts = as.data.frame(gtexCounts)
# rownames(gtexCounts) = gtexCounts$gene_id
# gtexCounts$gene_id = NULL
# gtexCounts$Description = NULL
# 
# gtexCountsFilt = gtexCounts[,rownames(GTExSampleDataFilt)[rownames(GTExSampleDataFilt) %in% colnames(gtexCounts)]]
# rm(gtexCounts)
# 
# 
# ## make count matrix tall
# nSamp = dim(gtexCountsFilt)[2]
# gtexCountsFilt$gene_id = rownames(gtexCountsFilt)
# rownames(gtexCountsFilt) = NULL
# 
# gtexCountsFilt = as_tibble(gtexCountsFilt) %>%
#   gather(SAMPID, counts, 1:nSamp)
# 
# gtexTPM = gtexCountsFilt %>%
#   inner_join(geneIDlengths, by = "gene_id") %>%
#   group_by(SAMPID) %>%
#   mutate(rpk = 1000*counts/length, rpkScalePerMillion = sum(rpk)/1000000, TPM = rpk/rpkScalePerMillion) %>%
#   dplyr::select(-c(rpk, rpkScalePerMillion))
# rm(gtexCountsFilt)
# 
# GTExTPMfilt = gtexTPM %>%
#   ungroup() %>%
#   dplyr::select(gene_id, SAMPID, TPM) %>%
#   spread(SAMPID, TPM)
# GTExTPMfilt = as.data.frame(GTExTPMfilt)
# rownames(GTExTPMfilt) = GTExTPMfilt$gene_id
# GTExTPMfilt$gene_id = NULL
# 
# GTExSampleDataFilt2 = GTExSampleDataFilt[rownames(GTExSampleDataFilt)[rownames(GTExSampleDataFilt) %in% colnames(GTExTPMfilt)],]
# GTExGeneData2 = gtexGeneData[rownames(GTExTPMfilt),]
# 
# tempSE = SummarizedExperiment(assays = list(TPMs= as.matrix(GTExTPMfilt)), colData = GTExSampleDataFilt2, rowData = GTExGeneData2)
# tempSEpath = "procDataScripted/GTEx/V6p/TPMs"
# 
# saveHDF5SummarizedExperiment(tempSE, dir = tempSEpath, verbose = T, replace = T)
# 
# # lengthtblGTEx = lengthtbl %>%
# #   separate(gene_id, into = c("geneid", "ver"), by = "\\.")
# # colnames(lengthtblGTEx) = c("gene_id", "ver", "GENCODElength")
# # lengthtblGTEx2 = lengthtblGTEx %>%
# #   inner_join(geneIDlengths, by = "gene_id")
# # lengthTblNoGTEx = geneIDlengths %>%
# #   anti_join(lengthtblGTEx, by = "gene_id")
# # lengthTblGTExNoInhouse = lengthtblGTEx %>%
# #   anti_join(geneIDlengths, by = "gene_id")
# 
# iCardNoGTExGenes = iCardTestATPMse[rowData(iCardTestATPMse)$gene_id %in% lengthTblNoGTEx$gene_id,]
# iCardInGTExGenes = iCardTestATPMse[rowData(iCardTestATPMse)$gene_id %in% lengthtblGTEx2$gene_id,]



# filter and save SummarizedExperiment for V7 GTEx data
GTExFile = "extractedData/RNAseq/GTEx/TPM_V7/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct"
GTExSampleDataFile1 = "extractedData/RNAseq/GTEx/TPM_V7/GTEx_v7_Annotations_SubjectPhenotypesDS.txt"
GTExSampleDataFile2 = "extractedData/RNAseq/GTEx/TPM_V7/GTEx_v7_Annotations_SampleAttributesDS.txt"

# GENCODEv19file = "extractedData/RNAseq/GTEx/allCountDataV6p/gencode.v19.annotation.gtf"
geneIDlengthsFile = "annotations/hg19_lengths_per_gene_id.tsv"
geneIDtoGeneSymbolFile = "annotations/hg19gene_idToGeneSymbol.tsv"


## process sample data
GTExSampleData1 = as_tibble(read.table(GTExSampleDataFile1, sep = "\t", header= T, quote = NULL, fill = T, stringsAsFactors = F))
GTExSampleData = as_tibble(read.table(GTExSampleDataFile2, sep = "\t", header= T, quote = "", stringsAsFactors = F)) 
cat('Regex_left_joining...\n')
GTExSampleData = GTExSampleData %>%
  regex_left_join(GTExSampleData1, by = c(SAMPID = "SUBJID"))
cat('Regex_left_join done.\n')
## get gene_id attributes
geneIDtoGeneSymbol = as_tibble(read.table(geneIDtoGeneSymbolFile, sep = "\t", header = T))
geneIDlengths = as_tibble(read.table(file = geneIDlengthsFile, sep = '\t', header = T))

cat('Loading GTEx TPMs...\n')
## load GTEx counts
gtexTPMs = as_tibble(read.table(file = GTExFile, skip = 2, sep = "\t", header = T, stringsAsFactors = F)) %>%
  tidyr::separate(Name, into = c("gene_id", "ver"), sep = "\\.") %>%
  dplyr::select(-ver)
cat('Loaded. Creating geneData...\n')

## make geneData (rowData in summarizedExperiment)
gtexGeneData = gtexTPMs %>%
  dplyr::select(gene_id, Description) %>%
  left_join(geneIDlengths, by = "gene_id")
colnames(gtexGeneData) = c("gene_id", "GeneSymbol", "length")
gtexGeneData = as.data.frame(gtexGeneData)
rownames(gtexGeneData) = gtexGeneData$gene_id
cat('Done. Filtering samples...\n')

## Prespecify sample filters to be used for downstream analysis
# Samples with metadata (15598 in SampleAttributes.txt file, but 11690 sample columns in count matrix...)
# Annotated "SMTSD" : a more refined sample source location ID
# RIN >= 7 : RNA integrity
# SMMAPRT >= 0.5 : mapping rate >= 50%
# SMMPPD >= 1e6 : 1MM or more mapped reads
## with these filters, there are 3623 samples left

GTExSampleDataFilt = GTExSampleData %>%
  filter(!is.na(SMTSD),
         SMRIN >= 7,
         SMMAPRT >= 0.5,
         SMMPPD >= 1e6)
GTExSampleDataFilt = as.data.frame(GTExSampleDataFilt)
rownames(GTExSampleDataFilt) = gsub("-",".", as.character(GTExSampleDataFilt$SAMPID))


## filter count matrix to include only pre-filtered samples
gtexTPMs = as.data.frame(gtexTPMs)
rownames(gtexTPMs) = gtexTPMs$gene_id
gtexTPMs$gene_id = NULL
gtexTPMs$Description = NULL

gtexTPMs = gtexTPMs[,rownames(GTExSampleDataFilt)[rownames(GTExSampleDataFilt) %in% colnames(gtexTPMs)]]

cat('Done. Saving...\n')

# ## make count matrix tall
# nSamp = dim(gtexTPMs)[2]
# gtexTPMs$gene_id = rownames(gtexTPMs)
# rownames(gtexCountsFilt) = NULL
# 
# gtexCountsFilt = as_tibble(gtexCountsFilt) %>%
#   gather(SAMPID, counts, 1:nSamp)
# 
# gtexTPM = gtexCountsFilt %>%
#   inner_join(geneIDlengths, by = "gene_id") %>%
#   group_by(SAMPID) %>%
#   mutate(rpk = 1000*counts/length, rpkScalePerMillion = sum(rpk)/1000000, TPM = rpk/rpkScalePerMillion) %>%
#   dplyr::select(-c(rpk, rpkScalePerMillion))
# rm(gtexCountsFilt)
# 
# GTExTPMfilt = gtexTPM %>%
#   ungroup() %>%
#   dplyr::select(gene_id, SAMPID, TPM) %>%
#   spread(SAMPID, TPM)
# GTExTPMfilt = as.data.frame(GTExTPMfilt)
# rownames(GTExTPMfilt) = GTExTPMfilt$gene_id
# GTExTPMfilt$gene_id = NULL

GTExSampleDataFilt2 = GTExSampleDataFilt[rownames(GTExSampleDataFilt)[rownames(GTExSampleDataFilt) %in% colnames(gtexTPMs)],]
GTExGeneData2 = gtexGeneData[rownames(gtexTPMs),]

tempSE = SummarizedExperiment(assays = list(TPMs= as.matrix(gtexTPMs)), colData = GTExSampleDataFilt2, rowData = GTExGeneData2)
tempSEpath = paste0(procDataDir, "/GTEx/V7/TPMs")

if(!dir.exists(paste0(procDataDir, "/GTEx"))){
  dir.create(paste0(procDataDir, "/GTEx"))
}
if(!dir.exists(paste0(procDataDir, "/GTEx/V7"))){
  dir.create(paste0(procDataDir, "/GTEx/V7"))
}
if(!dir.exists(paste0(procDataDir, "/GTEx/V7/TPMs"))){
  dir.create(paste0(procDataDir, "/GTEx/V7/TPMs"))
}

# saveHDF5SummarizedExperiment(tempSE, dir = tempSEpath, verbose = T, replace = T)
saveRDS(tempSE, paste0(tempSEpath, '/V7_TPMs.rds'))
cat('Done.\n')