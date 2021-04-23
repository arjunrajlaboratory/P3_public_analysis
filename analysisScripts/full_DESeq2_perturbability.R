#!/usr/bin/env Rscript
# don't run lines 3-13 if running from R workspace
args = commandArgs(trailingOnly = T)

# check if there are exactly 3 arguments. If not, return an error
if (length(args) < 4 | length(args) > 4) {
  stop("Exactly four arguments must be supplied: projectDir, procDataSubdir, graphSubdir, and runDEseq.", call.=FALSE)
} 
if (length(args)==4){
  projectDir = args[1]
  procDataSubdir = args[2]
  graphSubdir = args[3]
  runDEseq = as.logical(args[4]) # TRUE OR FALSE
}

# Main script for analyzing all 6 plates' worth of RNAtag-seq data
## Basic clustering/PCA on TPM data
## Running DESeq2
### Comparing gene-level results per-condition
## Calculating varaibility stats
### Plotting markers
### Doing KEGG/GO on high/low pert sets

library(DESeq2)
# library(HDF5Array)
# library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(magrittr)
library(e1071)

# if running manually in Rstudio, start here and use:
# projectDir = '~/Dropbox (RajLab)/Shared_IanM/cellid_201807_onward/'
# procDataSubdir = 'procDataScripted'
# graphSubdir = 'graphs'

procDataDir = paste0(projectDir, procDataSubdir)
graphDir = paste0(projectDir, graphSubdir)

if(!dir.exists(graphDir)){
  dir.create(graphDir)
}
setwd(projectDir)
cat("Starting full_DESeq2_perturbability.R...\n")
source(paste0(projectDir, "analysisScripts/RNAtagUtilities.R"))
cat(paste0("Working in ", getwd(), "\n"))

tfGeneIDfile = "annotations/TF_gene_ids.csv"
tfTab = as_tibble(read.csv(tfGeneIDfile, stringsAsFactors = F))
colnames(tfTab) = "gene_id"

geneIDfile = "annotations/hg19gene_idToGeneSymbol.tsv"
geneIDtoGeneSymbol <- as_tibble(read.csv(geneIDfile, sep='\t', stringsAsFactors = F))

markerFile = "miscellaneous/markerGeneList.csv"
markerTab = as_tibble(read.csv(markerFile, header = T, stringsAsFactors = F)) %>%
  left_join(geneIDtoGeneSymbol, by = "GeneSymbol")

barrierFile = "miscellaneous/fibro_barriers_GeneSymbols.txt"
barrierTab = as_tibble(read.table(barrierFile, header = T, stringsAsFactors = F)) %>%
  left_join(geneIDtoGeneSymbol, by = "GeneSymbol")

zhouOlson_activators = as_tibble(read.table("miscellaneous/ZhouOlson2017_TableS1_activators_reorg.txt", sep = "\t", header = T, stringsAsFactors = F))
zhouOlson_inhibitors = as_tibble(read.table("miscellaneous/ZhouOlson2017_TableS1_inhibitors_reorg.txt", sep = "\t", header = T, stringsAsFactors = F))
## Load data

allpathc = paste0(procDataDir, "/allExperiments/mappedCounts")
# allsec = loadHDF5SummarizedExperiment(allpathc)
allsec = readRDS(paste0(allpathc, '/se_mappedCounts.rds'))

allpatht = paste0(procDataDir, "/allExperiments/allTPMs")
# allset = loadHDF5SummarizedExperiment(allpatht)
allset = readRDS(paste0(allpatht, '/se_allTPMs.rds'))

# filtTPMs_old <- readRDS('~/Dropbox (RajLab)/Projects/cellid/procDataScripted/allExperiments/allTPMs_geneFilt_readFilt_manualFilt/se.rds')
## Get list of non-low-expressing genes
## i.e., genes whose average TPM in controls (per cell type) passes a threshold value. By default, accept meanTPM > 5.
tpm_filter = 5

allset_controls = allset[,grepl("DMSO", as.character(colData(allset)$condition)) &
                           (colData(allset)$cellType != "HeLa" | grepl("iCard", as.character(colData(allset)$experiment)))]

iCard_controls = allset_controls[,grepl("iCard", as.character(colData(allset_controls)$experiment))]
fibro_controls = allset_controls[,grepl("RNAtagTest", as.character(colData(allset_controls)$experiment))]

iCard_control_meanTPMs = rowMeans(assay(iCard_controls))
fibro_control_meanTPMs = rowMeans(assay(fibro_controls))

iCard_control_meanTPMs_filt = as_tibble(data.frame(
  meanTPM = iCard_control_meanTPMs[iCard_control_meanTPMs > tpm_filter],
  gene_id = names(iCard_control_meanTPMs)[iCard_control_meanTPMs > tpm_filter]
))
iCard_control_meanTPMs_all = as_tibble(data.frame(
  meanTPM = iCard_control_meanTPMs,
  gene_id = names(iCard_control_meanTPMs)
)) %>% mutate(cellType = "iCard")

fibro_control_meanTPMs_filt = as_tibble(data.frame(
  meanTPM = fibro_control_meanTPMs[fibro_control_meanTPMs > tpm_filter],
  gene_id = names(fibro_control_meanTPMs)[fibro_control_meanTPMs > tpm_filter]
))
fibro_control_meanTPMs_all = as_tibble(data.frame(
  meanTPM = fibro_control_meanTPMs,
  gene_id = names(fibro_control_meanTPMs)
)) %>% mutate(cellType = "GM00942")

controls_meanTPMs <- bind_rows(iCard_control_meanTPMs_all, fibro_control_meanTPMs_all)

  ## Filter allTPMs gene-wise and with minimum totalReads required (by default, at least 700,000)
read_filter = 700000

# fibroblasts
fibro_tpms = allset[,colData(allset)$cellType == "GM00942"]
fibro_tpms_geneFilt_readFilt = fibro_tpms[rowData(fibro_tpms)$gene_id %in% fibro_control_meanTPMs_filt$gene_id, 
                                          colData(fibro_tpms)$totalCounts > read_filter]

# iCards
iCard_tpms = allset[,grepl("iCard", as.character(colData(allset)$experiment))]
iCard_tpms_geneFilt_readFilt = iCard_tpms[rowData(iCard_tpms)$gene_id %in% iCard_control_meanTPMs_filt$gene_id, 
                                          colData(iCard_tpms)$totalCounts > read_filter]


dim(fibro_tpms_geneFilt_readFilt)
dim(iCard_tpms_geneFilt_readFilt)



## Clustering/PCA on TPMs
fibro_forClust = assay(fibro_tpms_geneFilt_readFilt)
dst = dist(t(fibro_forClust))
tree = hclust(dst, method = "average")
# plot(tree, main = "Hierarchical Clustering of fibroblast wells, 700k reads minimum")

iCard_forClust = assay(iCard_tpms_geneFilt_readFilt)
dst = dist(t(iCard_forClust))
tree = hclust(dst, method = "average")
# plot(tree, main = "Hierarchical Clustering of iCard wells, 700k reads minimum")


# consider genes that meet meanTPM threshold in both cell types
all_control_meanTPMs_filt = inner_join(iCard_control_meanTPMs_filt, fibro_control_meanTPMs_filt, by = "gene_id")

all_tpms_geneFilt_readFilt = allset[rowData(allset)$gene_id %in% all_control_meanTPMs_filt$gene_id, 
                                    colData(allset)$totalCounts > read_filter]

dst = dist(t(assay(all_tpms_geneFilt_readFilt)))
tree = hclust(dst, method = "average")
# plot(tree, main = "Hierarchical Clustering of all wells, 700k reads minimum")

sampleTypes <- as_tibble(colData(all_tpms_geneFilt_readFilt)) %>%
  dplyr::select(experiment, poolID, tagID, cellType)
  
if (!dir.exists(paste0(graphDir, '/exploratoryAnalysis'))){
  dir.create(paste0(graphDir, '/exploratoryAnalysis'))
}
tpmPCA<-prcomp(t(assay(all_tpms_geneFilt_readFilt)), scale = T)
tpmPCAforPlot<-as.data.frame(tpmPCA$x)
tpmPCAforPlot$combID <- rownames(tpmPCA$x)
rownames(tpmPCAforPlot) <- NULL
tpmPCAforPlot<-as_tibble(tpmPCAforPlot) %>%
  separate(combID, into = c("experiment","pool", "tag"), sep = "\\_")
tpmPCAsum<-summary(tpmPCA)

tpmPCAplot12 <- ggplot(tpmPCAforPlot, aes(PC1, PC2)) +
  geom_point(aes(color = experiment)) +
  theme_bw() + theme(axis.title = element_text(size = rel(1.5))) +
  xlab(paste0("PC1: ", as.character(tpmPCAsum$importance[2,1]), " variance explained")) +
  ylab(paste0("PC2: ", as.character(tpmPCAsum$importance[2,2]), " variance explained")) +
  ggtitle("PCA on filtered TPM matrix of all 6 plates: first two PCs\nMin TPM = 5")
# plot(tpmPCAplot12)
ggsave(tpmPCAplot12, file = paste0(graphDir, "/exploratoryAnalysis/fibsAndiCards_geneFilt_readFilt_manualFilt_PCA_TPMs.pdf"), width = 8, height = 8, useDingbats = F)


## Filter only to include drugs with average manualInclude >= 0.8 (out of 1)

all_sampleData = colData(allset)
rownames(all_sampleData) = NULL
all_sampleData_summary = as_tibble(all_sampleData) %>%
  # group_by(cellType, condition, Product.Name, Product.Name, Product.Name_drug1, Product.Name_drug2) %>%
  group_by(cellType, condition, Product.Name, Product.Name_drug1, Product.Name_drug2) %>%
  filter(!is.na(manualInclude_ManualFileTitle), totalCounts > read_filter) %>%
  summarise(includeSum = mean(as.numeric(manualInclude_ManualFileTitle)),
            numImages = length(cellType),
            meanReads = mean(totalCounts))

as_tibble(all_sampleData) %>%
  # group_by(cellType, condition, Product.Name, Product.Name, Product.Name_drug1, Product.Name_drug2) %>%
  group_by(cellType, condition, Product.Name, Product.Name_drug1, Product.Name_drug2) %>%
  filter(!is.na(manualInclude_ManualFileTitle), totalCounts > read_filter, as.numeric(manualInclude_ManualFileTitle) < 0.8)

badDrugs = all_sampleData_summary %>%
  filter(includeSum < 0.8)

# check clustering for other filterable drugs
# look at d19 - FLI-06; d40 - AEE788 (NVP-AEE788); d57 - AZ191; and pair_7 - GNF-2 + NVP-AEW541 
# if other drugs cluster with them maybe pull them (some drugs had 1 or 0 images)
iCard_forClust_drugNamed = iCard_forClust
colnames(iCard_forClust_drugNamed) = paste0(colData(allset)[colnames(iCard_forClust),"condition"],"_",
                                            colData(allset)[colnames(iCard_forClust),"experiment"],"_",
                                            colData(allset)[colnames(iCard_forClust),"rowLet"],
                                            as.character(colData(allset)[colnames(iCard_forClust),"colNum"]))
dst = dist(t(iCard_forClust_drugNamed))
tree = hclust(dst, method = "average")
# plot(tree, main = "Hierarchical Clustering of iCard wells, 700k reads minimum", cex=0.7)


# weirdly seeing one control well clustering with FLI-06 (iCardTest6set10x A10). 
# Maybe this well got plated weirdly (e.g., bad geltrex or something). 
# Since it's an edge well I don't have an image; will exclude out of an abundance of caution.

# Let's check to make sure other controls don't cluster too closely to it.
# check iCard control wells only
iCard_forClust_drugNamed_controls = iCard_forClust_drugNamed[,grepl("DMSO" ,colnames(iCard_forClust_drugNamed))]
dst = dist(t(iCard_forClust_drugNamed_controls))
tree = hclust(dst, method = "average")
# plot(tree, main = "Hierarchical Clustering of iCard control wells, 700k reads minimum", cex=0.7)


## not seeing anything too crazy. Will exclude just iCard samples from iCardTest6set10x_A10, pair_7, AEE788 (NVP-AEE788), AZ191, and FLI-06.
# fibroblast samples from those conditions either died in read-filter step or weren't manually excluded.

all_tpms_geneFilt_readFilt_manualFilt = all_tpms_geneFilt_readFilt[,colData(all_tpms_geneFilt_readFilt)$cellType == "GM00942" |
                                                                     (colData(all_tpms_geneFilt_readFilt)$condition != "d19" & 
                                                                        colData(all_tpms_geneFilt_readFilt)$condition != "d40" &
                                                                        colData(all_tpms_geneFilt_readFilt)$condition != "d57" &
                                                                        colData(all_tpms_geneFilt_readFilt)$condition != "pair_7" &
                                                                        paste0(colData(all_tpms_geneFilt_readFilt)$experiment,
                                                                               colData(all_tpms_geneFilt_readFilt)$rowLet,
                                                                               as.character(colData(all_tpms_geneFilt_readFilt)$colNum)) != "iCardTest6set10xA10")]

if(!dir.exists(paste0(procDataDir, "/allExperiments"))){
  dir.create(paste0(procDataDir, "/allExperiments"))
}
if(!dir.exists(paste0(procDataDir, "/allExperiments/allTPMs_geneFilt_readFilt_manualFilt"))){
  dir.create(paste0(procDataDir, "/allExperiments/allTPMs_geneFilt_readFilt_manualFilt"))
}

# saveHDF5SummarizedExperiment(all_tpms_geneFilt_readFilt_manualFilt, dir = paste0(procDataDir, "/allExperiments/allTPMs_geneFilt_readFilt_manualFilt"))
saveRDS(all_tpms_geneFilt_readFilt_manualFilt, paste0(procDataDir, "/allExperiments/allTPMs_geneFilt_readFilt_manualFilt/se_allTPMs_geneFilt_readFilt_manualFilt.rds"))
# saveRDS(all_tpms_geneFilt_readFilt_manualFilt_rerun, paste0(procDataDir, "/allExperiments/allTPMs_geneFilt_readFilt_manualFilt/se_allTPMs_geneFilt_readFilt_manualFilt_rerun.rds"))

dst = dist(t(assay(all_tpms_geneFilt_readFilt_manualFilt)))
tree = hclust(dst, method = "average")
# plot(tree, main = "Hierarchical Clustering of all wells, 700k reads minimum")

if (!dir.exists(paste0(graphDir, '/exploratoryAnalysis'))){
  dir.create(paste0(graphDir, '/exploratoryAnalysis'))
}
tpmPCA<-prcomp(t(assay(all_tpms_geneFilt_readFilt_manualFilt)), scale = T)
tpmPCAforPlot<-as.data.frame(tpmPCA$x)
tpmPCAforPlot$combID <- rownames(tpmPCA$x)

sampleDatPCA <- colData(all_tpms_geneFilt_readFilt_manualFilt)
condTblPCA <- tibble(
  combID = rownames(sampleDatPCA),
  condition = sampleDatPCA$condition,
  cellType = sampleDatPCA$cellType
)

tpmPCAforPlot %<>% inner_join(condTblPCA)%>%
  separate(combID, into = c("experiment","pool", "tag"), sep = "\\_")
rownames(tpmPCAforPlot) <- NULL
tpmPCAsum<-summary(tpmPCA)

tpmPCAplot12_byCondition <- ggplot(tpmPCAforPlot, aes(PC1, PC2)) +
  geom_point(aes(color = condition)) +
  theme_bw() + theme(axis.title = element_text(size = rel(1.5))) +
  xlab(paste0("PC1: ", as.character(tpmPCAsum$importance[2,1]), " variance explained")) +
  ylab(paste0("PC2: ", as.character(tpmPCAsum$importance[2,2]), " variance explained")) +
  ggtitle("PCA on filtered TPM matrix of all quality controlled samples:\nfirst two PCs\nMin TPM = 5")
# plot(tpmPCAplot12)
ggsave(tpmPCAplot12_byCondition, file = paste0(graphDir, "/exploratoryAnalysis/fibsAndiCards_geneFilt_readFilt_manualFilt_PCA_TPMs_byCondition.pdf"), width = 8, height = 8, useDingbats = F)


tpmPCAplot12_byCondition_byType <- ggplot(tpmPCAforPlot, aes(PC1, PC2)) +
  geom_point(aes(color = condition, shape = cellType)) +
  theme_bw() + theme(axis.title = element_text(size = rel(1.5))) +
  xlab(paste0("PC1: ", as.character(tpmPCAsum$importance[2,1]), " variance explained")) +
  ylab(paste0("PC2: ", as.character(tpmPCAsum$importance[2,2]), " variance explained")) +
  ggtitle("PCA on filtered TPM matrix of all quality controlled samples:\nfirst two PCs\nMin TPM = 5")
# plot(tpmPCAplot12)
ggsave(tpmPCAplot12_byCondition_byType, file = paste0(graphDir, "/exploratoryAnalysis/fibsAndiCards_geneFilt_readFilt_manualFilt_PCA_TPMs_byCondition_byType.pdf"), width = 8, height = 8, useDingbats = F)

tpmPCAplot12_byType <- ggplot(tpmPCAforPlot, aes(PC1, PC2)) +
  geom_point(aes(color = cellType)) +
  theme_bw() + theme(axis.title = element_text(size = rel(1.5))) +
  xlab(paste0("PC1: ", as.character(tpmPCAsum$importance[2,1]), " variance explained")) +
  ylab(paste0("PC2: ", as.character(tpmPCAsum$importance[2,2]), " variance explained")) +
  ggtitle("PCA on filtered TPM matrix of all quality controlled samples:\nfirst two PCs\nMin TPM = 5")
# plot(tpmPCAplot12)
ggsave(tpmPCAplot12_byType, file = paste0(graphDir, "/exploratoryAnalysis/fibsAndiCards_geneFilt_readFilt_manualFilt_PCA_TPMs_byCellType.pdf"), width = 8, height = 8, useDingbats = F)

## calculate transcripts per sample
# i.e., sum(all TPMs)*millions of mapped reads) = sum(all TPMs) * (totalCounts * (1-unmappedPerc) / 1e6)

tx_per_sample = colSums(assay(all_tpms_geneFilt_readFilt_manualFilt))*
                          colData(all_tpms_geneFilt_readFilt_manualFilt)$totalCounts*
                          (1-colData(all_tpms_geneFilt_readFilt_manualFilt)$unmappedPercentageOfTotal)/1e6

if(!dir.exists(paste0(graphDir, "/sampleQC/"))){
  dir.create(paste0(graphDir, "/sampleQC/"))
}
if(!dir.exists(paste0(graphDir, "/sampleQC/RNAtagSeq"))){
  dir.create(paste0(graphDir, "/sampleQC/RNAtagSeq"))
}

pdf(paste0(graphDir, "/sampleQC/RNAtagSeq/transcripts_per_sample_filtered.pdf"), width = 8, height = 8)
hist(log10(tx_per_sample), xlim = c(1, 8), xlab = "log10(Transcripts per sample)", main = "Estimated transcripts captured per sample")
dev.off();

pdf(paste0(graphDir, "/sampleQC/RNAtagSeq/counts_per_sample_filtered.pdf"), width = 8, height = 8)
hist(log10(colData(all_tpms_geneFilt_readFilt_manualFilt)$totalCounts), xlim = c(1, 8), xlab = "log10(Mapped reads per sample)", main = "Mapped reads per sample")
dev.off();
  
# by cell type  
iCard_tpms_geneFilt_readFilt_manualFilt = iCard_tpms_geneFilt_readFilt[,colData(iCard_tpms_geneFilt_readFilt)$condition != "d19" & 
                                                                         colData(iCard_tpms_geneFilt_readFilt)$condition != "d40" &
                                                                         colData(iCard_tpms_geneFilt_readFilt)$condition != "d57" &
                                                                         colData(iCard_tpms_geneFilt_readFilt)$condition != "pair_7" &
                                                                         paste0(colData(iCard_tpms_geneFilt_readFilt)$experiment,
                                                                                colData(iCard_tpms_geneFilt_readFilt)$rowLet,
                                                                                as.character(colData(iCard_tpms_geneFilt_readFilt)$colNum)) != "iCardTest6set10xA10"]

# saveHDF5SummarizedExperiment(iCard_tpms_geneFilt_readFilt_manualFilt, dir = paste0(procDataDir, "/allExperiments/iCards_geneFilt_readFilt_manualFilt"))
# saveRDS(iCard_tpms_geneFilt_readFilt_manualFilt, paste0(procDataDir, "/allExperiments/iCard_tpms_geneFilt_readFilt_manualFilt/se_iCard_tpms_geneFilt_readFilt_manualFilt.rds"))

## DESeq2 ####
# basic iCard controls vs fibroblast controls
controls_counts_readFilt_manualFilt <- allsec[,colData(allsec)$cellType != "HeLa" &
                                                colData(allsec)$condition %in% c("DMSOonly-control", "DMSOHeLa-HeLa")]
rowTemp = rowData(controls_counts_readFilt_manualFilt)
rownames(rowTemp) = rowTemp$gene_id

controls_counts_readFilt_manualFilt = SummarizedExperiment(assays = list(counts = as.matrix(assay(controls_counts_readFilt_manualFilt))), 
                                                           colData = colData(controls_counts_readFilt_manualFilt), 
                                                           rowData = rowTemp)

if(!dir.exists(paste0(graphDir, "/differentialExpression"))){
  dir.create(paste0(graphDir, "/differentialExpression"))
}
if(!dir.exists(paste0(graphDir, "/differentialExpression/allSamples"))){
  dir.create(paste0(graphDir, "/differentialExpression/allSamples"))
}

if(!runDEseq){
  print('Skipping DESeq2 step, using pre-run results.')
}
if(runDEseq){
  ddsSE <- DESeqDataSet(controls_counts_readFilt_manualFilt, design = ~ cellType)
  ddsSE <- ddsSE[ rowSums(counts(ddsSE)) > 1, ] #some super-basic gene_id filtering
  gc()
  ddsSE <- DESeq(ddsSE)
  results <- results(ddsSE, contrast = c("cellType", "iCard", "GM00942"))
  results$gene_id <- rownames(results)
  rownames(results) <- NULL
  write.table(results, paste0(graphDir, "/differentialExpression/allSamples/controls_both_DESeq2results.txt"), sep = "\t", quote = F, row.names = F)
  
  
  controls_both_DESeq2 <- read.table(paste0(graphDir, "/differentialExpression/allSamples/controls_both_DESeq2results.txt"), sep = "\t", stringsAsFactors = F, header = T)
  
  
  all_tpms_geneFilt_readFilt_manualFilt <- readRDS(paste0(procDataDir, "/allExperiments/allTPMs_geneFilt_readFilt_manualFilt/se_allTPMs_geneFilt_readFilt_manualFilt.rds"))
  # against all controls, per cell type
  #fibroblasts
  fibro_counts_readFilt_manualFilt = allsec[,colnames(assay(all_tpms_geneFilt_readFilt_manualFilt))]
  fibro_counts_readFilt_manualFilt = fibro_counts_readFilt_manualFilt[,colData(fibro_counts_readFilt_manualFilt)$cellType == "GM00942"]
  rowTemp = rowData(fibro_counts_readFilt_manualFilt)
  rownames(rowTemp) = rowTemp$gene_id
  
  fibro_counts_readFilt_manualFilt = SummarizedExperiment(assays = list(counts = as.matrix(assay(fibro_counts_readFilt_manualFilt))), 
                                                          colData = colData(fibro_counts_readFilt_manualFilt), 
                                                          rowData = rowTemp)
  if(!dir.exists(paste0(graphDir, "/differentialExpression/allSamples/fibroblasts"))){
    dir.create(paste0(graphDir, "/differentialExpression/allSamples/fibroblasts"))
  }
  if(!dir.exists(paste0(graphDir, "/differentialExpression/allSamples/iCards"))){
    dir.create(paste0(graphDir, "/differentialExpression/allSamples/iCards"))
  }
  print('Starting DESeq2 for fibroblast samples...\n')
  ddsSE <- DESeqDataSet(fibro_counts_readFilt_manualFilt, design = ~ condition)
  ddsSE <- ddsSE[ rowSums(counts(ddsSE)) > 1, ] #some super-basic gene_id filtering
  gc()
  ddsSE <- DESeq(ddsSE)
  saveRDS(ddsSE, paste0(graphDir, "/differentialExpression/allSamples/fibroblasts/fibro_dds_se.rds"))
  
  fgraphDir = paste0(graphDir, "/differentialExpression/allSamples/fibroblasts")
  allPFdrugsfile = paste0(fgraphDir, "/fibro_readFilt_manualFilt_MAplots.pdf")
  
  ds = exp(seq(0.1, 8, 0.1))
  
  print('Extracting and plotting DESeq2 results for fibroblasts...')
  pdf(allPFdrugsfile, width = 32, height = 36)
  par(mfrow = c(8,9))
  deTab = list()
  allResults = list()
  for (cond in unique(colData(ddsSE)$condition)[unique(colData(ddsSE)$condition) != "DMSOonly-control"]) {
    
    # overall demonstration of DESeq2 LFC detection for genes of different expression levels
    # per drug
    # 1. MA plots: mean expression vs LFC colored by significance at 3 FDR-adjusted p-val cutoffs
    # 2. calculation of nDE
    # 3. loess of positive lfc detected against mean and negative lfc against mean
    
    resultsTemp = results(ddsSE, contrast = c("condition", as.character(cond), "DMSOonly-control"))
    plotMA(resultsTemp, ylim=c(-2,2), main = paste0(cond, "_fibro_vsAllControls"))
    
    # nCtl = sum(colData(ddsSE)$condition == "DMSOonly-control")
    # nDrug = sum(colData(ddsSE)$condition == as.character(cond))
    # 
    nDE = sum(resultsTemp$padj[!is.na(resultsTemp$padj)] < 0.1)
    
    # # get CVs for all genes 
    # cvctl = sd(ddsSE[,colData(ddsSE)$condition == "DMSOonly-control"])/mean(ddsSE[,colData(ddsSE)$condition == "DMSOonly-control"])
    # 
    # cvs = left_join(cvctl, cvdrug)
    # 
    # for (dept in ds) {
    #   
    #   ps <- rnapower(depth = ds, cv = cvs[dept,], effect = 2, alpha = 0.1, n = c(nCtl, nDrug))
    #   
    #   
    #   
    # }
    deTemp = as.data.frame(list(
      condition = cond,
      nDE = nDE
    ))
    
    if(is.null(dim(deTab))){
      deTab = as_tibble(deTemp)
    } else {
      deTab = bind_rows(deTab, deTemp)
    }
    
    resultsTemp$gene_id = rownames(resultsTemp)
    rownames(resultsTemp) = NULL
    resultsTemp$Product.Name = cond
    
    if(is.null(dim(allResults))){
      allResults = as_tibble(resultsTemp)
    } else {
      allResults = bind_rows(allResults, as_tibble(resultsTemp))
    } 
    
    print(paste0("Finished ", cond, ". nDE = ", as.character(nDE)))
  }
  dev.off()
  write.table(allResults, file = paste0(graphDir, "/differentialExpression/allSamples/fibroblasts/fibro_readFilt_manualFilt_DESeqResults.txt"), sep = "\t", quote = F, row.names = F)
  write.table(deTab, file = paste0(graphDir, "/differentialExpression/allSamples/fibroblasts/fibro_readFilt_manualFilt_nDEtab.txt"), sep = "\t", row.names = F, quote = F)
  
  print('Fibroblasts done. Starting iCard samples...')
  #iCards
  iCard_counts_readFilt_manualFilt = allsec[,colnames(assay(all_tpms_geneFilt_readFilt_manualFilt))]
  iCard_counts_readFilt_manualFilt = iCard_counts_readFilt_manualFilt[,colData(iCard_counts_readFilt_manualFilt)$cellType == "iCard"]
  rowTemp = rowData(iCard_counts_readFilt_manualFilt)
  rownames(rowTemp) = rowTemp$gene_id
  
  iCard_counts_readFilt_manualFilt = SummarizedExperiment(assays = list(counts = as.matrix(assay(iCard_counts_readFilt_manualFilt))), 
                                                          colData = colData(iCard_counts_readFilt_manualFilt), 
                                                          rowData = rowTemp)
  
  ddsSE <- DESeqDataSet(iCard_counts_readFilt_manualFilt, design = ~ condition)
  ddsSE <- ddsSE[ rowSums(counts(ddsSE)) > 1, ] #some super-basic gene_id filtering
  gc()
  ddsSE <- DESeq(ddsSE)
  saveRDS(ddsSE, paste0(graphDir, "/differentialExpression/allSamples/iCards/iCards_dds_se.rds"))
  # graphDir = paste0(graphDir, "differentialExpression/allSamples/iCards/")
  
  allPFdrugsfile = paste0(graphDir, "/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_MAplots.pdf")
  
  print('Extracting and plotting DESeq2 results for iCards...')
  pdf(allPFdrugsfile, width = 32, height = 36)
  par(mfrow = c(8,9))
  deTab = list()
  allResults = list()
  for (cond in unique(colData(ddsSE)$condition)[unique(colData(ddsSE)$condition) != "DMSOonly-control"]) {
    
    resultsTemp = results(ddsSE, contrast = c("condition", as.character(cond), "DMSOonly-control"))
    plotMA(resultsTemp, ylim=c(-2,2), main = paste0(cond, "_iCard_vsAllControls"))
    
    nDE = sum(resultsTemp$padj[!is.na(resultsTemp$padj)] < 0.1)
    
    deTemp = as.data.frame(list(
      condition = cond,
      nDE = nDE
    ))
    
    if(is.null(dim(deTab))){
      deTab = as_tibble(deTemp)
    } else {
      deTab = bind_rows(deTab, deTemp)
    }
    
    resultsTemp$gene_id = rownames(resultsTemp)
    rownames(resultsTemp) = NULL
    resultsTemp$Product.Name = cond
    
    if(is.null(dim(allResults))){
      allResults = as_tibble(resultsTemp)
    } else {
      allResults = bind_rows(allResults, as_tibble(resultsTemp))
    } 
    
    print(paste0("Finished ", cond, ". nDE = ", as.character(nDE)))
  }
  dev.off()
  write.table(allResults, file = paste0(graphDir, "/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_DESeqResults.txt"), sep = "\t", quote = F, row.names = F)
  write.table(deTab, file = paste0(graphDir, "/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_nDEtab.txt"), sep = "\t", row.names = F, quote = F)
  
  print('iCards done. Starting downsampling analysis...')
  ## downsample counts and calculate nDE per drug
  
  downsamp_list <- c(0.5, 0.7, 0.9)
  seed_list <- c(2380, 53250)
  
  for (seedID in seed_list) {
    # seedID <- 2380
    set.seed(seedID)
    cat(paste0('Using random seed = ', as.character(seedID), '.\n'))
    for (fracReads in downsamp_list) {
      
      cat(paste0('Downsampling counts to ', as.character(fracReads), ' of total.\n'))
      tempCounts <- downsample_counts(as.matrix(assay(iCard_counts_readFilt_manualFilt)), fracReads) 
      
      iCard_counts_readFilt_manualFilt_temp = SummarizedExperiment(assays = list(counts = tempCounts), 
                                                                   colData = colData(iCard_counts_readFilt_manualFilt), 
                                                                   rowData = rowTemp)
      
      ddsSE <- DESeqDataSet(iCard_counts_readFilt_manualFilt_temp, design = ~ condition)
      ddsSE <- ddsSE[ rowSums(counts(ddsSE)) > 1, ] #some super-basic gene_id filtering
      gc()
      ddsSE <- DESeq(ddsSE)
      
      # graphDir = paste0(graphDir, "differentialExpression/allSamples/iCards/")
      if(!dir.exists(paste0(graphDir, "/differentialExpression/allSamples/iCards/downsampled"))){
        dir.create(paste0(graphDir, "/differentialExpression/allSamples/iCards/downsampled"))
      }
      
      sePath = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled')
      saveRDS(iCard_counts_readFilt_manualFilt_temp, paste0(sePath, '/se_mappedCounts_',as.character(fracReads*100), "percent_seed_", as.character(seedID),'.rds'))
      
      allPFdrugsfile = paste0(graphDir, "/differentialExpression/allSamples/iCards/downsampled/iCard_readFilt_manualFilt_MAplots_", as.character(fracReads*100), "percent_seed_", as.character(seedID),".pdf")
      
      pdf(allPFdrugsfile, width = 32, height = 36)
      par(mfrow = c(8,9))
      deTab = list()
      allResults = list()
      for (cond in unique(colData(ddsSE)$condition)[unique(colData(ddsSE)$condition) != "DMSOonly-control"]) {
        
        resultsTemp = results(ddsSE, contrast = c("condition", as.character(cond), "DMSOonly-control"))
        plotMA(resultsTemp, ylim=c(-2,2), main = paste0(cond, "_iCard_vsAllControls"))
        
        nDE = sum(resultsTemp$padj[!is.na(resultsTemp$padj)] < 0.1)
        
        deTemp = as.data.frame(list(
          condition = cond,
          nDE = nDE
        ))
        
        if(is.null(dim(deTab))){
          deTab = as_tibble(deTemp)
        } else {
          deTab = bind_rows(deTab, deTemp)
        }
        
        resultsTemp$gene_id = rownames(resultsTemp)
        rownames(resultsTemp) = NULL
        resultsTemp$Product.Name = cond
        
        if(is.null(dim(allResults))){
          allResults = as_tibble(resultsTemp)
        } else {
          allResults = bind_rows(allResults, as_tibble(resultsTemp))
        } 
        
        print(paste0("Finished ", cond, ". nDE = ", as.character(nDE)))
      }
      dev.off()
      write.table(allResults, file = paste0(graphDir, "/differentialExpression/allSamples/iCards/downsampled/iCard_readFilt_manualFilt_DESeqResults_", as.character(fracReads*100), "percent_seed_",as.character(seedID),".txt"), sep = "\t", quote = F, row.names = F)
      write.table(deTab, file = paste0(graphDir, "/differentialExpression/allSamples/iCards/downsampled/iCard_readFilt_manualFilt_nDEtab_", as.character(fracReads*100), "percent_seed_",as.character(seedID), ".txt"), sep = "\t", row.names = F, quote = F)
      
      
    }
  }
  nDE_ds_iCard <- list()
  for (fracReads in downsamp_list) {
    
    nDE_tab <- as_tibble(read.table(paste0(graphDir, "/differentialExpression/allSamples/iCards/downsampled/iCard_readFilt_manualFilt_nDEtab_", as.character(fracReads*100), "percent.txt"), stringsAsFactors = F, header = T, sep = '\t')) %>%
      mutate(fracReads = fracReads)
    
    if (is.null(dim(nDE_ds_iCard))) {
      nDE_ds_iCard <- nDE_tab
    } else {
      nDE_ds_iCard %<>% bind_rows(nDE_tab)
    }
    
  }
  # ggplot(nDE_ds_iCard, aes(fracReads, nDE, group = condition)) +
  #   geom_line()
  print('Done with downsampling.')
} # end first runDEseq
controls_both_DESeq2 <- read.table(paste0(graphDir, "/differentialExpression/allSamples/controls_both_DESeq2results.txt"), sep = "\t", stringsAsFactors = F, header = T)
 
# per-plate DESeq2 calls
if(!runDEseq){print('Skipping DESeq2 step, using pre-run results.')}
if(runDEseq){
  print('Running per-plate DESeq...')
  if(!dir.exists(paste0(graphDir, "/differentialExpression/allSamples/iCards"))){
    dir.create(paste0(graphDir, "/differentialExpression/allSamples/iCards"))
  }
  allPFdrugsfile = paste0(graphDir, "/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_MAplots_perPlate.pdf")
  
  pdf(allPFdrugsfile, width = 32, height = 36)
  par(mfrow = c(8,9))
  deTab = list()
  allResults = list()
  
  plates = unique(colData(iCard_counts_readFilt_manualFilt)$experiment)
  
  for (plate in plates) {
    
    cat('Filtering to plate ', plate, '\n')
    ddsSEp <- DESeqDataSet(iCard_counts_readFilt_manualFilt[,colData(iCard_counts_readFilt_manualFilt)$experiment == plate], design = ~ condition)
    ddsSEp <- ddsSEp[ rowSums(counts(ddsSEp)) > 1, ] #some super-basic gene_id filtering
    gc()
    cat('Doing plate-wise Diff exp calculation\n')
    ddsSEp <- DESeq(ddsSEp)
    
    for (cond in unique(colData(ddsSEp)$condition)[unique(colData(ddsSEp)$condition) != "DMSOonly-control"]) {
      
      resultsTemp = results(ddsSEp, contrast = c("condition", as.character(cond), "DMSOonly-control"))
      plotMA(resultsTemp, ylim=c(-2,2), main = paste0(cond, "_iCard_vsAllControls"))
      
      nDE = sum(resultsTemp$padj[!is.na(resultsTemp$padj)] < 0.1)
      
      deTemp = as.data.frame(list(
        condition = cond,
        plateID = plate,
        nDE = nDE
      ))
      
      if(is.null(dim(deTab))){
        deTab = as_tibble(deTemp)
      } else {
        deTab = bind_rows(deTab, deTemp)
      }
      
      resultsTemp$gene_id = rownames(resultsTemp)
      rownames(resultsTemp) = NULL
      resultsTemp$Product.Name = cond
      resultsTemp$plateID = plate
      
      if(is.null(dim(allResults))){
        allResults = as_tibble(resultsTemp)
      } else {
        allResults = bind_rows(allResults, as_tibble(resultsTemp))
      } 
      
      print(paste0("Finished ", cond, ". nDE = ", as.character(nDE)))
    }
  }
  dev.off()
  write.table(allResults, file = paste0(graphDir, "/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_DESeqResults_perPlate.txt"), sep = "\t", quote = F, row.names = F)
  write.table(deTab, file = paste0(graphDir, "/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_nDEtab_perPlate.txt"), sep = "\t", row.names = F, quote = F)
} # end runDEseq

print('Starting Effect Size filter analysis...')
## DESeq2 frequency ####
## DESeq2-based perturbability stats = nDESeq2conditionsAll, nDESeq2conditionsUp, nDESeq2conditionsDown
deResultsiCard_temp = as_tibble(read.table(paste0(graphDir, "/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_DESeqResults.txt"), header = T, sep = "\t", stringsAsFactors = F)) %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.1) %>%
  mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) %>%
  dplyr::select(gene_id, Product.Name, direction, log2FoldChange) %>%
  mutate(cellType = "iCard") %>%
  group_by(cellType, gene_id) %>%
  summarise(nDESeq2conditionsAll = length(direction), 
            nDESeq2conditionsUp = sum(direction == "up"), 
            nDESeq2conditionsDown = sum(direction == "down"),
            nDESeq2conditionsAll_0_5 = sum(abs(log2FoldChange) >= 0.5),
            nDESeq2conditionsUp_0_5 = sum(abs(log2FoldChange) >= 0.5 & direction == "up"),
            nDESeq2conditionsDown_0_5 = sum(abs(log2FoldChange) >= 0.5 & direction == "down"),
            nDESeq2conditionsAll_0_75 = sum(abs(log2FoldChange) >= 0.75),
            nDESeq2conditionsUp_0_75 = sum(abs(log2FoldChange) >= 0.75 & direction == "up"),
            nDESeq2conditionsDown_0_75 = sum(abs(log2FoldChange) >= 0.75 & direction == "down"))

deResultsfibro_temp = as_tibble(read.table(paste0(graphDir, "/differentialExpression/allSamples/fibroblasts/fibro_readFilt_manualFilt_DESeqResults.txt"),header = T, sep = "\t", stringsAsFactors = F)) %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.1) %>%
  mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) %>%
  dplyr::select(gene_id, Product.Name, direction, log2FoldChange) %>%
  mutate(cellType = "GM00942") %>%
  group_by(cellType, gene_id) %>%
  summarise(nDESeq2conditionsAll = length(direction), 
            nDESeq2conditionsUp = sum(direction == "up"), 
            nDESeq2conditionsDown = sum(direction == "down"),
            nDESeq2conditionsAll_0_5 = sum(abs(log2FoldChange) >= 0.5),
            nDESeq2conditionsUp_0_5 = sum(abs(log2FoldChange) >= 0.5 & direction == "up"),
            nDESeq2conditionsDown_0_5 = sum(abs(log2FoldChange) >= 0.5 & direction == "down"),
            nDESeq2conditionsAll_0_75 = sum(abs(log2FoldChange) >= 0.75),
            nDESeq2conditionsUp_0_75 = sum(abs(log2FoldChange) >= 0.75 & direction == "up"),
            nDESeq2conditionsDown_0_75 = sum(abs(log2FoldChange) >= 0.75 & direction == "down"))

DESeq2_perturbability = bind_rows(deResultsiCard_temp, deResultsfibro_temp)

## DESeq2 effectSize ####
## Calculate average up or down effect sizes of perturbations to each gene
## For conditions in which perturbed
deResultsiCard_temp2 = as_tibble(read.table(paste0(graphDir, "/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_DESeqResults.txt"), header = T, sep = "\t", stringsAsFactors = F)) %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.1)

deResultsiCard_effectSize = deResultsiCard_temp2 %>%
  mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) %>%
  group_by(gene_id, direction) %>%
  summarise(meanLFC = mean(log2FoldChange),
            sdLFC = sd(log2FoldChange),
            nLFC = length(log2FoldChange),
            semLFC = sdLFC/sqrt(nLFC)) %>%
  gather('variable', 'value', meanLFC:semLFC) %>%
  unite(temp, direction, variable) %>%
  spread(temp, value) %>%
  mutate(cellType = "iCard")

deResultsfibro_temp2 = as_tibble(read.table(paste0(graphDir, "/differentialExpression/allSamples/fibroblasts/fibro_readFilt_manualFilt_DESeqResults.txt"), header = T, sep = "\t", stringsAsFactors = F)) %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.1)

deResultsfibro_effectSize = deResultsfibro_temp2 %>%
  mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) %>%
  group_by(gene_id, direction) %>%
  summarise(meanLFC = mean(log2FoldChange),
            sdLFC = sd(log2FoldChange),
            nLFC = length(log2FoldChange),
            semLFC = sdLFC/sqrt(nLFC)) %>%
  gather('variable', 'value', meanLFC:semLFC) %>%
  unite(temp, direction, variable) %>%
  spread(temp, value) %>%
  mutate(cellType = "GM00942")


DESeq2_perturbability_effectsize <- bind_rows(deResultsiCard_effectSize, deResultsfibro_effectSize)
DESeq2_perturbability_effectsize[is.na(DESeq2_perturbability_effectsize)] <- 0

DESeq2_perturbability_controls_byCellTypeTemp <- as_tibble(controls_both_DESeq2) %>%
  dplyr::select(log2FoldChange, padj, gene_id) %>%
  dplyr::rename(LFC_iCardOverFib = log2FoldChange,
                padj_iCardOverFib = padj)

DESeq2_perturbability %<>% 
  left_join(DESeq2_perturbability_effectsize, by = c("cellType", "gene_id")) %>%
  left_join(DESeq2_perturbability_controls_byCellTypeTemp, by = "gene_id")


## Bootstrap fracUP
# use deResultsiCard_temp
# iCard_DEres <- as_tibble(read.table(paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_DESeqResults.txt'), sep = '\t', header = T, stringsAsFactors = F))
# use deResultsfibro_temp
# fibro_DEres <- as_tibble(read.table(paste0(graphDir, '/differentialExpression/allSamples/fibroblasts/fibro_readFilt_manualFilt_DESeqResults.txt'), sep = '\t', header = T, stringsAsFactors = F))

# iCards

if(runDEseq){
  print('Bootstrapping DE condition counts per gene...')
  # filter by: padj < 0.1
  iCard_DEres_sig <- as_tibble(read.table(paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_DESeqResults.txt'), sep = '\t', header = T, stringsAsFactors = F)) %>%
    filter(padj < 0.1) %>%
    mutate(directionEffect = ifelse(log2FoldChange < 0, 'DOWN', 'UP'))
  
  iCard_drugs <- unique(iCard_DEres_sig$Product.Name)
  nDrugs <- length(iCard_drugs)
  
  iCard_DEres_sig_wide <- iCard_DEres_sig %>%
    dplyr::select(gene_id, Product.Name, directionEffect) %>%
    spread(Product.Name, directionEffect)
  
  iCard_DEres_sig_wide_mat <- as.data.frame(iCard_DEres_sig_wide)
  rownames(iCard_DEres_sig_wide_mat) <- iCard_DEres_sig_wide_mat$gene_id
  iCard_DEres_sig_wide_mat$gene_id <- NULL
  iCard_DEres_sig_wide_mat[is.na(iCard_DEres_sig_wide_mat)] <- 'NoSigDiff'
  
  nIter <- 1000
  set.seed(97603)
  
  iCard_boot_up <- matrix(0, nrow(iCard_DEres_sig_wide_mat), nIter)
  iCard_boot_down <- matrix(0, nrow(iCard_DEres_sig_wide_mat), nIter)
  iCard_boot_noch <- matrix(0, nrow(iCard_DEres_sig_wide_mat), nIter)
  
  rownames(iCard_boot_up) <- rownames(iCard_DEres_sig_wide_mat)
  rownames(iCard_boot_down) <- rownames(iCard_DEres_sig_wide_mat)
  rownames(iCard_boot_noch) <- rownames(iCard_DEres_sig_wide_mat)
  for (i in 1:nIter){
    
    # sample with replacement from drugs list
    dlist <- sample(iCard_drugs, nDrugs, replace = T)
    
    # get up/down/none results (including replicated results)
    temp_wide_mat <- iCard_DEres_sig_wide_mat[,dlist]
    
    # calculate nUp, nDown, nAll, fracUp
    temp_up <- rowSums(temp_wide_mat == 'UP')
    temp_down <- rowSums(temp_wide_mat == 'DOWN')
    temp_noch <- rowSums(temp_wide_mat == 'NoSigDiff')
    
    # store values in separate tables per gene
    iCard_boot_up[names(temp_up),i] <- as.numeric(temp_up)
    iCard_boot_down[names(temp_down),i] <- as.numeric(temp_down)
    iCard_boot_noch[names(temp_noch),i] <- as.numeric(temp_noch)
    
    
  }
  
  write.csv(iCard_boot_up, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/bootstrap_DESeq2_up.csv'), quote = F, row.names = T)
  write.csv(iCard_boot_down, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/bootstrap_DESeq2_down.csv'), quote = F, row.names = T)
  write.csv(iCard_boot_noch, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/bootstrap_DESeq2_noch.csv'), quote = F, row.names = T)
  print('Done with bootstrapping.')
} # end runDEseq

##### Calculate distribution-based perturbability stats - TO INCLUDE?
## then merge with DESeq2 summary
print('Calculating moment-based perturbability values. Individual sample-based analysis first...')
# on individual samples
all_counts_readFilt_manualFilt = allsec[,colnames(assay(all_tpms_geneFilt_readFilt_manualFilt))]
all_counts_readFilt_manualFilt_noHeLa = all_counts_readFilt_manualFilt[,colData(all_counts_readFilt_manualFilt)$cellType != "HeLa"]

coldat_for_tall <- as.data.frame(colData(all_counts_readFilt_manualFilt_noHeLa)) %>%
  dplyr::select(condition, cellType)

all_counts_readFilt_manualFilt_noHeLa <- SummarizedExperiment(assays = list(counts = assay(all_counts_readFilt_manualFilt_noHeLa)), 
                                                              colData = coldat_for_tall, 
                                                              rowData = rowData(all_counts_readFilt_manualFilt_noHeLa))

all_rpm_readFilt_manualFilt_tall = SummarizedExperimentToTallTibble(all_counts_readFilt_manualFilt_noHeLa) %>%
  mutate(isControl = ifelse(grepl("DMSO", condition), "control", "perturbed")) %>%
  mutate(counts = TPM) %>%
  dplyr::select(-TPM) %>%
  group_by(combID) %>%
  mutate(totalMappedCounts = sum(counts), RPM = 1e6*counts/totalMappedCounts) %>%
  dplyr::select(combID, condition, cellType, isControl, gene_id, RPM)
gc()

control_meanRPM_meanTPM = all_rpm_readFilt_manualFilt_tall %>%
  filter(isControl == "control") %>%
  group_by(cellType, gene_id) %>%
  summarise(meanRPM = mean(RPM)) %>%
  left_join(controls_meanTPMs, by = c("cellType", "gene_id"))

all_rpm_readFilt_manualFilt_variability_metrics = calculate_perturbability(all_rpm_readFilt_manualFilt_tall, control_meanRPM_meanTPM) %>%
  left_join(DESeq2_perturbability, by = c("cellType", "gene_id")) 
all_rpm_readFilt_manualFilt_variability_metrics[is.na(all_rpm_readFilt_manualFilt_variability_metrics)] <- 0 # for genes with no DEconds detected

print('Done with individual samples. Working on grouping by drug...')
# on grouped samples
# first, group samples
perturbed_grouped_meanRPM_readFilt_manualFilt_tall = all_rpm_readFilt_manualFilt_tall %>%
  ungroup() %>% filter(isControl == "perturbed") %>%
  group_by(cellType, isControl, condition, gene_id) %>%
  summarise(meanRPM = mean(RPM))

## Group controls into pseudogroups
controlsPerPlate = all_rpm_readFilt_manualFilt_tall %>% ungroup() %>%
  separate(combID, into = c("experiment", "poolID", "tagID"), sep = "\\_", remove = F) %>%
  filter(isControl == "control") %>%
  dplyr::select(cellType, experiment, combID) %>%
  unique() %>%
  group_by(cellType, experiment) %>%
  summarise(nCtls = length(cellType))

samplesPerConditionPerturbed = all_rpm_readFilt_manualFilt_tall %>% ungroup() %>%
  separate(combID, into = c("experiment", "poolID", "tagID"), sep = "\\_", remove = F) %>%
  filter(isControl == "perturbed") %>%
  dplyr::select(cellType, experiment, condition, combID) %>%
  unique() %>%
  group_by(cellType, experiment, condition) %>%
  summarise(nSamps = length(cellType))

controlCombIDs = all_rpm_readFilt_manualFilt_tall %>% ungroup() %>%
  separate(combID, into = c("experiment", "poolID", "tagID"), sep = "\\_", remove = F) %>%
  filter(isControl == "control") %>%
  dplyr::select(cellType, experiment, combID) %>%
  unique()

controlGroupList = list()
controlTab = list()
allControlPseudoGroups = list()
controlGroupIDs = list() #
nIter = 100

set.seed(2054) # 2054 for arranging all nIter pseudogrouping draws
# for control pseudogroup sizes, sample from the distribution of filter-passing drug group sizes for each plate.

# pick combID-pseudogroupID assignment draws
for (i in 1:nrow(controlsPerPlate)){
  
  cells = controlsPerPlate$cellType[i]
  plate = controlsPerPlate$experiment[i]
  nCtls = controlsPerPlate$nCtls[i]
  
  controlGroups = integer()
  controlIDs = character()
  g=0
  while(sum(controlGroups) < nCtls){
    g = g+1
    thisSamp = sample(samplesPerConditionPerturbed$nSamps[samplesPerConditionPerturbed$experiment == plate], 1)
    # thisSamp = sample(iCardGroupCounts, 1)
    controlGroups = c(controlGroups, thisSamp)
    controlIDs = c(controlIDs, rep(paste0(plate, "_control", as.character(g)), times = thisSamp))
    
  }
  
  controlCombIDsThisPlate = controlCombIDs$combID[controlCombIDs$experiment == plate]
  
  # loop here for bootstrapping
  reorderedControlIDs = sample(controlIDs, nCtls)
  controlGroupIDs[[plate]] = matrix("empty", nCtls, nIter)
  rownames(controlGroupIDs[[plate]]) = controlCombIDsThisPlate
  for (j in 1:nIter) {
    
    controlGroupIDs[[plate]][,j] = sample(controlIDs, nCtls)[1:nCtls]
    
  }
  
}

# use above control group assignments to calculate meanRPM values for each gene over pseudogroups
# pick one draw at random and proceed with analysis based on that
# bootstrapping to come later

control_rpm_readFilt_manualFilt_tall = all_rpm_readFilt_manualFilt_tall %>% ungroup() %>%
  filter(isControl == "control") 
seedID = 4930
set.seed(seedID) # 4930 for pseudogroup selection
i = sample(1:nIter, 1)

ids = c(rownames(controlGroupIDs[[1]]), 
        rownames(controlGroupIDs[[2]]),
        rownames(controlGroupIDs[[3]]),
        rownames(controlGroupIDs[[4]]),
        rownames(controlGroupIDs[[5]]),
        rownames(controlGroupIDs[[6]]))
groupTab = as_tibble(data.frame(
  combID = as.character(ids),
  controlGroup = as.character(c(controlGroupIDs[[1]][,i], 
                                controlGroupIDs[[2]][,i],
                                controlGroupIDs[[3]][,i],
                                controlGroupIDs[[4]][,i],
                                controlGroupIDs[[5]][,i],
                                controlGroupIDs[[6]][,i]))))

control_grouped_meanRPM_readFilt_manualFilt_tall = all_rpm_readFilt_manualFilt_tall %>%
  ungroup() %>% 
  inner_join(groupTab, by = "combID") %>%
  mutate(condition = controlGroup) %>%
  group_by(cellType, isControl, condition, gene_id) %>%
  summarise(meanRPM = mean(RPM))

all_grouped_meanRPM_readFilt_manualFilt_tall = bind_rows(perturbed_grouped_meanRPM_readFilt_manualFilt_tall, control_grouped_meanRPM_readFilt_manualFilt_tall)
gc()


all_groupedMeanRPM_readFilt_manualFilt_variability_metrics = calculate_perturbability(all_grouped_meanRPM_readFilt_manualFilt_tall %>% mutate(RPM = meanRPM) %>% dplyr::select(-meanRPM), 
                                                                                      control_meanRPM_meanTPM) %>%
  left_join(DESeq2_perturbability, by = c("cellType", "gene_id")) 
gc()
write.table(all_rpm_readFilt_manualFilt_variability_metrics, 
            file = paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_variability_metrics.txt"),
            sep = "\t", quote = F, row.names = F)

write.table(all_groupedMeanRPM_readFilt_manualFilt_variability_metrics, 
            file = paste0(procDataDir, "/allExperiments/all_groupedMeanRPM_readFilt_manualFilt_variability_metrics_seed", as.character(seedID), ".txt"),
            sep = "\t", quote = F, row.names = F)
print('Done with grouped analysis. Now meanRPM sliding window normalization...')

### Normalized perturbability measures to expression
## 1. sliding window along expression of all genes passing tpm_filter
## 2., sliding window along expression of TF genes, and compare only to nearest N genes in expression level
## Start with N = 40 (20 higher, 20 lower) "windowRadius = 20"
## normalize to control RPM
cat("Working on pertuability normalization with the default windowRadius = 20 and gene-level filter of meanTPM > 10.\n")
tpm_filter = 5 # originally did 10 to save time.
windowRadius = 20
# rank by expression
# all genes passing TPM filter
all_rpm_readFilt_manualFilt_tpmFilt_variability_RPMrank = all_rpm_readFilt_manualFilt_variability_metrics %>% 
  filter(meanTPM > tpm_filter) %>% 
  group_by(cellType) %>%
  mutate(expRank = rank(meanRPM))

all_groupedMeanRPM_readFilt_manualFilt_tpmFilt_variability_RPMrank = all_groupedMeanRPM_readFilt_manualFilt_variability_metrics %>% 
  filter(meanTPM > tpm_filter) %>% 
  group_by(cellType) %>%
  mutate(expRank = rank(meanRPM))

# only TFs
all_rpm_readFilt_manualFilt_tfOnly_variability_tfOnlyRPMrank = all_rpm_readFilt_manualFilt_variability_metrics %>% 
  inner_join(tfTab, by = "gene_id") %>%
  group_by(cellType) %>%
  mutate(expRank = rank(meanRPM))


# only TFs passing TPM filter
all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_tfOnlyRPMrank = all_rpm_readFilt_manualFilt_variability_metrics %>% 
  filter(meanTPM > tpm_filter) %>% 
  inner_join(tfTab, by = "gene_id") %>%
  group_by(cellType) %>%
  mutate(expRank = rank(meanRPM))

all_groupedMeanRPM_readFilt_manualFilt_tfOnly_tpmFilt_variability_tfOnlyRPMrank = all_groupedMeanRPM_readFilt_manualFilt_variability_metrics %>% 
  filter(meanTPM > tpm_filter) %>% 
  inner_join(tfTab, by = "gene_id") %>%
  group_by(cellType) %>%
  mutate(expRank = rank(meanRPM))

all_rpm_readFilt_manualFilt_tfOnly_variability_tfOnlyRPMrank[is.na(all_rpm_readFilt_manualFilt_tfOnly_variability_tfOnlyRPMrank)] <- 0
all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_tfOnlyRPMrank[is.na(all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_tfOnlyRPMrank)] <- 0

all_rpm_readFilt_manualFilt_tfOnly_variability_RPMnormMetrics = normalize_perturbability(all_rpm_readFilt_manualFilt_tfOnly_variability_tfOnlyRPMrank)
all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics = normalize_perturbability(all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_tfOnlyRPMrank)
all_rpm_readFilt_manualFilt_tpmFilt_variability_RPMnormMetrics = normalize_perturbability(all_rpm_readFilt_manualFilt_tpmFilt_variability_RPMrank)

all_groupedMeanRPM_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics = normalize_perturbability(all_groupedMeanRPM_readFilt_manualFilt_tfOnly_tpmFilt_variability_tfOnlyRPMrank)
all_groupedMeanRPM_readFilt_manualFilt_tpmFilt_variability_RPMnormMetrics = normalize_perturbability(all_groupedMeanRPM_readFilt_manualFilt_tpmFilt_variability_RPMrank)
gc()
write.table(all_rpm_readFilt_manualFilt_tfOnly_variability_RPMnormMetrics,
            file = paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_tfOnly_variability_RPMnormMetrics_window", as.character(windowRadius),  ".txt"),
            sep = "\t", quote = F, row.names = F)
write.table(all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics,
            file = paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_window", as.character(windowRadius),  ".txt"),
            sep = "\t", quote = F, row.names = F)
write.table(all_rpm_readFilt_manualFilt_tpmFilt_variability_RPMnormMetrics,
            file = paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_tpmFilt_variability_RPMnormMetrics_window", as.character(windowRadius),  ".txt"),
            sep = "\t", quote = F, row.names = F)
write.table(all_groupedMeanRPM_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics,
            file = paste0(procDataDir, "/allExperiments/all_groupedMeanRPM_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_window", as.character(windowRadius), "_seed", as.character(seedID),  ".txt"),
            sep = "\t", quote = F, row.names = F)
write.table(all_groupedMeanRPM_readFilt_manualFilt_tpmFilt_variability_RPMnormMetrics,
            file = paste0(procDataDir, "/allExperiments/all_groupedMeanRPM_readFilt_manualFilt_tpmFilt_variability_RPMnormMetrics_window", as.character(windowRadius), "_seed", as.character(seedID),  ".txt"),
            sep = "\t", quote = F, row.names = F)
print('Done with sliding window normalization.')
gc()
print('Analyzing cumulative DE genes per cell type...')
## compare DESeq2 results within each drug or drug-pair
indiv_perturbability = as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_variability_metrics.txt"), header = T, stringsAsFactors = F, sep = "\t"))
deResultsiCard = as_tibble(read.table(paste0(graphDir, "/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_DESeqResults.txt"), header = T, sep = "\t", stringsAsFactors = F)) %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.1) %>%
  dplyr::select(gene_id, Product.Name) %>%
  mutate(cellType = "iCard")

deResultsfibro = as_tibble(read.table(paste0(graphDir, "/differentialExpression/allSamples/fibroblasts/fibro_readFilt_manualFilt_DESeqResults.txt"),header = T, sep = "\t", stringsAsFactors = F)) %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.1) %>%
  dplyr::select(gene_id, Product.Name) %>%
  mutate(cellType = "fibroblast")

onlyiCardPerDrug = anti_join(deResultsiCard, deResultsfibro, by = c("Product.Name", "gene_id")) %>%
  group_by(cellType, Product.Name) %>%
  summarise(nUni = length(gene_id))
onlyfibroPerDrug = anti_join(deResultsfibro, deResultsiCard, by = c("Product.Name", "gene_id")) %>%
  group_by(cellType, Product.Name) %>%
  summarise(nUni = length(gene_id))

sumDEiCard = deResultsiCard %>%
  group_by(cellType, Product.Name) %>%
  summarise(nDE = length(gene_id)) %>%
  left_join(onlyiCardPerDrug, by = c("cellType", "Product.Name")) %>%
  filter(!is.na(nUni)) %>%
  mutate(nOverlap = nDE - nUni) %>%
  dplyr::select(-nDE) %>%
  gather("geneMembership", "nGenes", 3:4)

sumDEfibro = deResultsfibro %>%
  group_by(cellType, Product.Name) %>%
  summarise(nDE = length(gene_id)) %>%
  left_join(onlyfibroPerDrug, by = c("cellType", "Product.Name")) %>%
  filter(!is.na(nUni)) %>%
  mutate(nOverlap = nDE - nUni) %>%
  dplyr::select(-nDE) %>%
  gather("geneMembership", "nGenes", 3:4)

b1 = ggplot(bind_rows(sumDEiCard, sumDEfibro) %>% filter(Product.Name != 'DMSOHeLa-HeLa'), aes(x = Product.Name, y = nGenes, fill = geneMembership)) +
  geom_bar(stat = "identity") +
  facet_grid(cellType~.) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90))
ggsave(b1, file = paste0(graphDir, "/differentialExpression/allSamples/cellTypeOverlapBars.pdf"))


nDEsumTab = bind_rows(deResultsfibro, deResultsiCard) %>%
  group_by(cellType, Product.Name) %>%
  summarise(nDE = length(gene_id))

iCardConds = nDEsumTab$Product.Name[nDEsumTab$cellType == "iCard"]
fibroConds = nDEsumTab$Product.Name[nDEsumTab$cellType == "fibroblast"]

iCardConds_notInFibro = iCardConds[!(iCardConds %in% fibroConds)]
fibroConds_notIniCard = fibroConds[!(fibroConds %in% iCardConds)]

zeroTab = as_tibble(data.frame(
  cellType = c(rep("iCard", times = length(iCardConds_notInFibro)), rep("fibroblast", times = length(fibroConds_notIniCard))),
  Product.Name = c(iCardConds_notInFibro, fibroConds_notIniCard),
  nDE = 0
))

nDEsumTab = bind_rows(nDEsumTab, zeroTab) %>% 
  separate(Product.Name, into = c('d','num'), sep = 'd') %>% 
  mutate(Product.Name = ifelse(d == '', ifelse(as.numeric(num) < 10, paste0('d0',num),paste0('d',num)),d)) %>% # fix ordering
  dplyr::select(cellType, Product.Name, nDE) %>%
  separate(Product.Name, into = c('pair','num'), sep = 'pair_') %>% 
  mutate(Product.Name = ifelse(pair == '', ifelse(as.numeric(num) < 10, paste0('pair_0',num),paste0('pair_',num)),pair)) %>% # fix ordering
  dplyr::select(cellType, Product.Name, nDE)

b2 = ggplot(nDEsumTab %>% filter(Product.Name != 'DMSOHeLa-HeLa'), aes(x = Product.Name, y = nDE)) +
  geom_bar(stat = "identity") +
  facet_grid(cellType~.) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90))
ggsave(b2, file = paste0(graphDir, "/differentialExpression/allSamples/nDE_bars_perCondition.pdf"), width = 12, height = 7)

b2log = ggplot(nDEsumTab %>% filter(Product.Name != 'DMSOHeLa-HeLa'), aes(x = Product.Name, y = log10(nDE))) +
  geom_bar(stat = "identity") +
  facet_grid(cellType~.) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90))
ggsave(b2log, file = paste0(graphDir, "/differentialExpression/allSamples/log10nDE_bars_perCondition.pdf"), width = 12, height = 7)


# cumulativeDE <- ggplot(inner_join(deResultsiCard, indiv_perturbability, by = c('gene_id', 'cellType') %>%
#                                     filter(meanRPM > 20) %>%
#                                     group_by))

deResultsfibro_minRPM20 <- deResultsfibro %>%
  inner_join(indiv_perturbability %>% filter(cellType == 'GM00942') %>% dplyr::select(gene_id, meanRPM), by = "gene_id") %>%
  filter(meanRPM >= 20)

deResultsiCard_minRPM20 <- deResultsiCard %>%
  inner_join(indiv_perturbability %>% filter(cellType == 'iCard') %>% dplyr::select(gene_id, meanRPM), by = "gene_id") %>%
  filter(meanRPM >= 20)

set.seed(97153)
fibconds = unique(deResultsfibro_minRPM20$Product.Name)
cardconds = unique(deResultsiCard_minRPM20$Product.Name)
fibcondsR = fibconds[sample(1:length(fibconds), length(fibconds))]
cardcondsR = cardconds[sample(1:length(cardconds), length(cardconds))]

fibUniqueResults = checkUniqueDESeq(deResultsfibro_minRPM20, fibcondsR)

cardUniqueResults = checkUniqueDESeq(deResultsiCard_minRPM20, cardcondsR)

fibUniqueResults_TF = checkUniqueDESeqTFs(deResultsfibro_minRPM20, fibcondsR, tfTab)

cardUniqueResults_TF = checkUniqueDESeqTFs(deResultsiCard_minRPM20, cardcondsR, tfTab)


nGenesExp <- indiv_perturbability %>% 
  dplyr::select(gene_id, meanRPM, cellType) %>%
  unique() %>%
  group_by(cellType) %>%
  summarise(nGE20rpm = sum(meanRPM >= 20)) %>%
  mutate(cellType = case_when(
    cellType == 'GM00942' ~ 'fibroblast',
    cellType == 'iCard' ~ 'iCard'
  ))

cumulative_uniques <- bind_rows(fibUniqueResults$cumulativeUniqueDEgenes %>% mutate(cellType = 'fibroblast'),
                                cardUniqueResults$cumulativeUniqueDEgenes %>% mutate(cellType = 'iCard')) %>%
  inner_join(nGenesExp, by = 'cellType')

write.csv(cumulative_uniques, file = paste0(graphDir, "/differentialExpression/allSamples/cumulative_unique_DEgenes_rpm20.csv"), quote = F)

cumulative_unique_plot <- ggplot() +
  geom_point(data = cumulative_uniques, aes(order, totalNew)) +
  geom_line(data = cumulative_uniques, aes(order, totalNew)) +
  geom_hline(data = cumulative_uniques, aes(yintercept = nGE20rpm), color = 'blue') +
  geom_text_repel(data = nGenesExp, aes(x = 20, y = nGE20rpm), label = '# Genes > 20 RPM', color = 'blue') +
  facet_wrap(~cellType, scales = 'free_y') +
  xlab('Number of conditions considered') +
  ylab('Cumulative differentially expressed\ngenes detected') +
  theme_bw()
ggsave(cumulative_unique_plot, file = paste0(graphDir, "/differentialExpression/allSamples/cumulative_unique_DEgenes_rpm20.pdf"), width = 7.5, height = 3, useDingbats = F)

nGenesExp_TF <- indiv_perturbability %>% ungroup() %>%
  dplyr::select(gene_id, meanRPM, cellType) %>%
  inner_join(tfTab, by = 'gene_id') %>%
  unique() %>%
  group_by(cellType) %>%
  summarise(nGE20rpm = sum(meanRPM >= 20)) %>%
  mutate(cellType = case_when(
    cellType == 'GM00942' ~ 'fibroblast',
    cellType == 'iCard' ~ 'iCard'
  ))


cumulative_uniques_TF <- bind_rows(fibUniqueResults_TF$cumulativeUniqueDEgenes %>% mutate(cellType = 'fibroblast'),
                                cardUniqueResults_TF$cumulativeUniqueDEgenes %>% mutate(cellType = 'iCard')) %>%
  inner_join(nGenesExp, by = 'cellType')

write.csv(cumulative_uniques_TF, file = paste0(graphDir, "/differentialExpression/allSamples/cumulative_unique_DEgenesTFonly_rpm20.csv"), quote = F)

cumulative_unique_plot_TF <- ggplot() +
  geom_point(data = cumulative_uniques_TF, aes(order, totalNew)) +
  geom_line(data = cumulative_uniques_TF, aes(order, totalNew)) +
  geom_hline(data = nGenesExp_TF, aes(yintercept = nGE20rpm), color = 'blue') +
  geom_text_repel(data = nGenesExp_TF, aes(x = 20, y = nGE20rpm), label = '# TF genes > 20 RPM', color = 'blue') +
  facet_wrap(~cellType, scales = 'free_y') +
  xlab('Number of conditions considered') +
  ylab('Cumulative differentially expressed\ntranscription factor genes detected') +
  theme_bw()
ggsave(cumulative_unique_plot_TF, file = paste0(graphDir, "/differentialExpression/allSamples/cumulative_unique_DEgenes_TFonly_rpm20.pdf"), width = 7.5, height = 3, useDingbats = F)

print('Done with cumulative DE analysis.')

print('Analyzing downsampled results...')
downsamp_list <- c(0.5, 0.7, 0.9) 
seed_list <- c(2380, 53250)
all_sig_res_ds <- list()
for (fracReads in downsamp_list) {
  
  sig_genes <- as_tibble(read.table(paste0(graphDir, "/differentialExpression/allSamples/iCards/downsampled/iCard_readFilt_manualFilt_DESeqResults_", as.character(fracReads*100), "percent.txt"), stringsAsFactors = F, header = T, sep = '\t')) %>%
    mutate(fracReads = fracReads) %>%
    filter(padj <= 0.1)
  
  if(is.null(dim(all_sig_res_ds))) {
    all_sig_res_ds <- sig_genes
  } else {
    all_sig_res_ds %<>% bind_rows(sig_genes)
  }
  
}

all_sig_res_full <- as_tibble(read.table(paste0(graphDir, "/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_DESeqResults.txt"), stringsAsFactors = F, header = T, sep = '\t')) %>%
  mutate(fracReads = 1) %>%
  filter(padj <= 0.1)

all_res_full <- as_tibble(read.table(paste0(graphDir, "/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_DESeqResults.txt"), stringsAsFactors = F, header = T, sep = '\t')) %>%
  mutate(fracReads = 1) 

all_sig_res_ds %<>%
  bind_rows(all_sig_res_full)

temp_means_tpm10 <- as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_variability_metrics.txt"), header = T, stringsAsFactors = F, sep = "\t")) %>% 
  filter(meanTPM > 10, cellType == "iCard") %>%
  dplyr::select(gene_id, meanTPM)

all_sig_res_temp_tpm10 <- all_sig_res_full %>%
  inner_join(temp_means_tpm10) %>%
  group_by(Product.Name) %>%
  summarise(total_nDE = length(Product.Name))

temp_means_rpm20 <- as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_variability_metrics.txt"), header = T, stringsAsFactors = F, sep = "\t")) %>% 
  filter(meanRPM > 20, cellType == "iCard") %>%
  dplyr::select(gene_id, meanRPM)

all_sig_res_temp_rpm20 <- all_sig_res_full %>%
  inner_join(temp_means_rpm20) %>%
  group_by(Product.Name) %>%
  summarise(total_nDE = length(Product.Name))

all_sig_res_ds_all_res_full_temp_rpm20 <- all_res_full %>%
  bind_rows(all_sig_res_ds) %>%
  inner_join(temp_means_rpm20) 

all_sig_res_ds_sum_tpm10 <- all_sig_res_ds %>%
  inner_join(temp_means_tpm10) %>%
  group_by(fracReads, Product.Name) %>%
  summarise(nDE = length(fracReads)) %>%
  inner_join(all_sig_res_temp_tpm10, by = 'Product.Name') %>%
  mutate(frac_nDE = nDE/total_nDE)

# ggplot(all_sig_res_ds_sum_tpm10, aes(fracReads, nDE, group = Product.Name)) +
#   geom_line()

# ggplot(all_sig_res_ds_sum_tpm10 %>% filter(total_nDE > 100), aes(fracReads, frac_nDE, group = Product.Name)) +
#   geom_line() +
#   ylim(c(0,1.2)) +
#   xlim(c(0,1)) +
#   geom_abline(slope = 1, intercept = 0, color = 'blue')

all_sig_res_ds_sum_rpm20 <- all_sig_res_ds %>%
  inner_join(temp_means_rpm20) %>%
  group_by(fracReads, Product.Name) %>%
  summarise(nDE = length(fracReads)) %>%
  inner_join(all_sig_res_temp_rpm20, by = 'Product.Name') %>%
  mutate(frac_nDE = nDE/total_nDE)

# ggplot(all_sig_res_ds_sum_rpm20, aes(fracReads, nDE, group = Product.Name)) +
#   geom_line()
# ggplot(all_sig_res_ds_sum_rpm20 %>% filter(total_nDE > 100), aes(fracReads, frac_nDE, group = Product.Name)) +
#   geom_line() +
#   ylim(c(0,1.2)) +
#   xlim(c(0,1)) +
#   geom_abline(slope = 1, intercept = 0, color = 'blue')

## mutli-seed downsampling analysis
all_sig_res_ds <- list()
for (seedID in seed_list){
  for (fracReads in downsamp_list) {
    
    sig_genes <- as_tibble(read.table(paste0(graphDir, "/differentialExpression/allSamples/iCards/downsampled/iCard_readFilt_manualFilt_DESeqResults_", as.character(fracReads*100), "percent_seed_",as.character(seedID),".txt"), stringsAsFactors = F, header = T, sep = '\t')) %>%
      mutate(fracReads = fracReads,
             seed = seedID) %>%
      filter(padj <= 0.1)
    
    if(is.null(dim(all_sig_res_ds))) {
      all_sig_res_ds <- sig_genes
    } else {
      all_sig_res_ds %<>% bind_rows(sig_genes)
    }
    
  }
}
all_sig_res_full <- as_tibble(read.table(paste0(graphDir, "/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_DESeqResults.txt"), stringsAsFactors = F, header = T, sep = '\t')) %>%
  mutate(fracReads = 1) %>%
  filter(padj <= 0.1)

for (seedID in seed_list) {
  
  all_sig_res_ds %<>%
    bind_rows(all_sig_res_full %>% mutate(seed = seedID))
  
}

all_res_full <- as_tibble(read.table(paste0(graphDir, "/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_DESeqResults.txt"), stringsAsFactors = F, header = T, sep = '\t')) %>%
  mutate(fracReads = 1) 

# all_sig_res_ds %<>%
#   bind_rows(all_sig_res_full)

temp_means_tpm10 <- as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_variability_metrics.txt"), header = T, stringsAsFactors = F, sep = "\t")) %>% 
  filter(meanTPM > 10, cellType == "iCard") %>%
  dplyr::select(gene_id, meanTPM)

all_sig_res_temp_tpm10 <- all_sig_res_full %>%
  inner_join(temp_means_tpm10) %>%
  group_by(Product.Name) %>%
  summarise(total_nDE = length(Product.Name))

temp_means_rpm20 <- as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_variability_metrics.txt"), header = T, stringsAsFactors = F, sep = "\t")) %>% 
  filter(meanRPM > 20, cellType == "iCard") %>%
  dplyr::select(gene_id, meanRPM)

temp_means_rpm10 <- as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_variability_metrics.txt"), header = T, stringsAsFactors = F, sep = "\t")) %>% 
  filter(meanRPM > 10, cellType == "iCard") %>%
  dplyr::select(gene_id, meanRPM)

temp_means_rpm5 <- as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_variability_metrics.txt"), header = T, stringsAsFactors = F, sep = "\t")) %>% 
  filter(meanRPM > 5, cellType == "iCard") %>%
  dplyr::select(gene_id, meanRPM)

temp_means_nofilt <- as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_variability_metrics.txt"), header = T, stringsAsFactors = F, sep = "\t")) %>% 
  filter(cellType == "iCard") %>%
  dplyr::select(gene_id, meanRPM)

all_sig_res_temp_rpm20 <- all_sig_res_full %>%
  inner_join(temp_means_rpm20) %>%
  group_by(Product.Name) %>%
  summarise(total_nDE = length(Product.Name))
all_sig_res_temp_rpm10 <- all_sig_res_full %>%
  inner_join(temp_means_rpm10) %>%
  group_by(Product.Name) %>%
  summarise(total_nDE = length(Product.Name))
all_sig_res_temp_rpm5 <- all_sig_res_full %>%
  inner_join(temp_means_rpm5) %>%
  group_by(Product.Name) %>%
  summarise(total_nDE = length(Product.Name))
all_sig_res_temp_rpmNone <- all_sig_res_full %>%
  inner_join(temp_means_nofilt) %>%
  group_by(Product.Name) %>%
  summarise(total_nDE = length(Product.Name))

all_sig_res_ds_sum_rpm20 <- all_sig_res_ds %>%
  inner_join(temp_means_rpm20) %>%
  group_by(fracReads, seed, Product.Name) %>%
  summarise(nDE = length(fracReads)) %>%
  inner_join(all_sig_res_temp_rpm20, by = 'Product.Name') %>%
  mutate(frac_nDE = nDE/total_nDE,
         rpmFilter = 20)

all_sig_res_ds_sum_rpm10 <- all_sig_res_ds %>%
  inner_join(temp_means_rpm10) %>%
  group_by(fracReads, seed, Product.Name) %>%
  summarise(nDE = length(fracReads)) %>%
  inner_join(all_sig_res_temp_rpm10, by = 'Product.Name') %>%
  mutate(frac_nDE = nDE/total_nDE,
         rpmFilter = 10)

all_sig_res_ds_sum_rpm5 <- all_sig_res_ds %>%
  inner_join(temp_means_rpm5) %>%
  group_by(fracReads, seed, Product.Name) %>%
  summarise(nDE = length(fracReads)) %>%
  inner_join(all_sig_res_temp_rpm5, by = 'Product.Name') %>%
  mutate(frac_nDE = nDE/total_nDE,
         rpmFilter = 5)

all_sig_res_ds_sum_rpmNone <- all_sig_res_ds %>%
  inner_join(temp_means_nofilt) %>%
  group_by(fracReads, seed, Product.Name) %>%
  summarise(nDE = length(fracReads)) %>%
  inner_join(all_sig_res_temp_rpmNone, by = 'Product.Name') %>%
  mutate(frac_nDE = nDE/total_nDE,
         rpmFilter = 0)

all_sig_res_ds_rpmFilters <- bind_rows(all_sig_res_ds_sum_rpm20, all_sig_res_ds_sum_rpm10) %>%
  bind_rows(all_sig_res_ds_sum_rpm5) %>%
  bind_rows(all_sig_res_ds_sum_rpmNone)

all_sig_res_ds_avg_rpm20 <- all_sig_res_ds_sum_rpm20 %>%
  filter(total_nDE > 100) %>%
  group_by(seed, fracReads, rpmFilter) %>%
  summarise(mean_frac_nDE = mean(frac_nDE),
            median_frac_nDE = median(frac_nDE),
            sem_frac_nDE = sd(frac_nDE)/sqrt(length(frac_nDE)))

all_sig_res_ds_avg_rpmFilters <- all_sig_res_ds_rpmFilters %>%
  filter(total_nDE > 100) %>%
  group_by(seed, fracReads, rpmFilter) %>%
  summarise(mean_frac_nDE = mean(frac_nDE),
            median_frac_nDE = median(frac_nDE),
            sem_frac_nDE = sd(frac_nDE)/sqrt(length(frac_nDE)))

iCard_ds_rpm20_sum <- ggplot() +
  geom_line(data = all_sig_res_ds_sum_rpm20 %>% filter(total_nDE > 100), aes(fracReads, frac_nDE, group = Product.Name)) +
  geom_point(data = all_sig_res_ds_avg_rpm20, aes(fracReads, mean_frac_nDE), color = 'green', size = 3) +
  geom_errorbar(data = all_sig_res_ds_avg_rpm20, aes(fracReads, ymin = mean_frac_nDE - sem_frac_nDE, ymax = mean_frac_nDE + sem_frac_nDE), color = 'green', width = 0.07) +
  facet_wrap(~seed) +
  ylim(c(0,1.2)) +
  xlim(c(0,1)) +
  xlab('fraction of reads downsampled') +
  ylab('fraction of total diff. expr. genes detected') +
  ggtitle('Downsampling of iCard RNAtag-seq\nTwo random seeds\nOnly genes with mean RPM > 20, Padj cutoff = 0.1') +
  geom_abline(slope = 1, intercept = 0, color = 'blue') +
  theme_bw()
ggsave(iCard_ds_rpm20_sum, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/iCard_downsampled_frac_nDE_twoSeeds.pdf'), width = 10, height = 5, useDingbats = F)

iCard_ds_rpm20_sum_median_2380only <- ggplot() +
  geom_line(data = all_sig_res_ds_sum_rpm20 %>% filter(total_nDE > 100, seed == 2380), aes(fracReads, frac_nDE, group = Product.Name)) +
  geom_line(data = all_sig_res_ds_avg_rpm20 %>% filter(seed == 2380), aes(fracReads, median_frac_nDE), color = 'green', size = 2) +
  # facet_wrap(~seed) +
  ylim(c(0,1.2)) +
  xlim(c(0,1)) +
  xlab('fraction of reads downsampled') +
  ylab('fraction of total diff. expr. genes detected') +
  ggtitle('Downsampling of iCard RNAtag-seq\nOnly genes with mean RPM > 20, Padj cutoff = 0.1,\nGreen = median') +
  geom_abline(slope = 1, intercept = 0, color = 'blue') +
  theme_bw()
ggsave(iCard_ds_rpm20_sum_median_2380only, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/iCard_downsampled_frac_nDE_twoSeeds_median_2380only.pdf'), width = 5, height = 5, useDingbats = F)

iCard_ds_rpm20_sum_median_53250only <- ggplot() +
  geom_line(data = all_sig_res_ds_sum_rpm20 %>% filter(total_nDE > 100, seed == 53250), aes(fracReads, frac_nDE, group = Product.Name)) +
  geom_line(data = all_sig_res_ds_avg_rpm20 %>% filter(seed == 53250), aes(fracReads, median_frac_nDE), color = 'green', size = 2) +
  # facet_wrap(~seed) +
  ylim(c(0,1.2)) +
  xlim(c(0,1)) +
  xlab('fraction of reads downsampled') +
  ylab('fraction of total diff. expr. genes detected') +
  ggtitle('Downsampling of iCard RNAtag-seq\nOnly genes with mean RPM > 20, Padj cutoff = 0.1,\nGreen = median') +
  geom_abline(slope = 1, intercept = 0, color = 'blue') +
  theme_bw()
ggsave(iCard_ds_rpm20_sum_median_53250only, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/iCard_downsampled_frac_nDE_twoSeeds_median_53250only.pdf'), width = 5, height = 5, useDingbats = F)

all_sig_res_ds_avg_rpm20_fracChanges <- all_sig_res_ds_avg_rpm20 %>%
  group_by(seed) %>%
  dplyr::select(median_frac_nDE, seed, fracReads) %>%
  spread(fracReads, median_frac_nDE) %>%
  mutate(perc0_50_gain = `0.5`,
         perc50_70_gain = `0.7` - `0.5`,
         perc70_90_gain = `0.9` - `0.7`,
         perc90_100_gain = `1` - `0.9`) %>%
  dplyr::select(-c(`1` , `0.9`, `0.7` , `0.5`)) %>%
  gather('read_change','delta_median_frac_nDE',perc0_50_gain:perc90_100_gain) %>%
  mutate(delta_median_frac_nDE_per10pc = case_when(
    read_change == 'perc0_50_gain' ~ delta_median_frac_nDE/5,
    read_change == 'perc50_70_gain' ~ delta_median_frac_nDE/2,
    read_change == 'perc70_90_gain' ~ delta_median_frac_nDE/2,
    read_change == 'perc90_100_gain' ~ delta_median_frac_nDE
  ),
  frac_reads = case_when(
    read_change == 'perc0_50_gain' ~ 0.5,
    read_change == 'perc50_70_gain' ~ 0.7,
    read_change == 'perc70_90_gain' ~ 0.9,
    read_change == 'perc90_100_gain' ~ 1
  ))

all_sig_res_ds_avg_rpm20_fracChanges_meanDelta <- all_sig_res_ds_avg_rpm20 %>%
  group_by(seed) %>%
  dplyr::select(median_frac_nDE, seed, fracReads) %>%
  spread(fracReads, median_frac_nDE) %>%
  mutate(perc0_50_gain = `0.5`,
         perc50_70_gain = `0.7` - `0.5`,
         perc70_90_gain = `0.9` - `0.7`,
         perc90_100_gain = `1` - `0.9`) %>%
  dplyr::select(-c(`1` , `0.9`, `0.7` , `0.5`)) %>%
  gather('read_change','delta_median_frac_nDE',perc0_50_gain:perc90_100_gain) %>%
  mutate(delta_median_frac_nDE_per10pc = case_when(
    read_change == 'perc0_50_gain' ~ delta_median_frac_nDE/5,
    read_change == 'perc50_70_gain' ~ delta_median_frac_nDE/2,
    read_change == 'perc70_90_gain' ~ delta_median_frac_nDE/2,
    read_change == 'perc90_100_gain' ~ delta_median_frac_nDE
  ),
  frac_reads = case_when(
    read_change == 'perc0_50_gain' ~ 0.5,
    read_change == 'perc50_70_gain' ~ 0.7,
    read_change == 'perc70_90_gain' ~ 0.9,
    read_change == 'perc90_100_gain' ~ 1
  )) %>%
  group_by(read_change) %>%
  summarise(mean_delta = mean(delta_median_frac_nDE),
            mean_delta_per10pc = mean(delta_median_frac_nDE_per10pc),
            frac_reads = mean(frac_reads))


iCard_ds_rpm20_deltaFrac_nDE_2seeds <- ggplot() +
  geom_point(data = all_sig_res_ds_avg_rpm20_fracChanges, aes(read_change, delta_median_frac_nDE)) +
  geom_bar(data = all_sig_res_ds_avg_rpm20_fracChanges_meanDelta, aes(read_change, mean_delta), stat = 'identity', alpha = 0.4) +
  theme_bw() +
  ylab('Fractional change in total DE genes\nBy change in downsampled reads') +
  ggtitle('Fraction of detected DE genes downsampling analysis\nMin RPM = 20, Adj. p-val <= 0.1\nTwo random seeds tested')
ggsave(iCard_ds_rpm20_deltaFrac_nDE_2seeds, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/iCard_ds_rpm20_deltaFrac_nDE_2seeds.pdf'), width = 6, height = 6, useDingbats = F)

iCard_ds_rpm20_deltaFracPer10pc_nDE_2seeds <- ggplot() +
  geom_point(data = all_sig_res_ds_avg_rpm20_fracChanges, aes(read_change, delta_median_frac_nDE_per10pc)) +
  geom_bar(data = all_sig_res_ds_avg_rpm20_fracChanges_meanDelta, aes(read_change, mean_delta_per10pc), stat = 'identity', alpha = 0.4) +
  theme_bw() +
  ylab('Fractional change in total DE genes\nNormalized to 10% increment\nBy change in downsampled reads') +
  ggtitle('Fraction of detected DE genes downsampling analysis\nNormalized to 10% increment, Min RPM = 20, Adj. p-val <= 0.1\nTwo random seeds tested')
ggsave(iCard_ds_rpm20_deltaFracPer10pc_nDE_2seeds, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/iCard_ds_rpm20_deltaFracPer10pc_nDE_2seeds.pdf'), width = 6, height = 6, useDingbats = F)


iCard_ds_rpm20_deltaFrac_vsfracReads_nDE_2seeds <- ggplot() +
  geom_point(data = all_sig_res_ds_avg_rpm20_fracChanges, aes(frac_reads, delta_median_frac_nDE)) +
  geom_line(data = all_sig_res_ds_avg_rpm20_fracChanges_meanDelta, aes(frac_reads, mean_delta)) +
  theme_bw() +
  ylab('Fractional change in total DE genes\nBy change in downsampled reads') +
  xlab('Fraction of downsampled reads') +
  ggtitle('Fraction of detected DE genes downsampling analysis\nMin RPM = 20, Adj. p-val <= 0.1\nTwo random seeds tested')
ggsave(iCard_ds_rpm20_deltaFrac_vsfracReads_nDE_2seeds, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/iCard_ds_rpm20_deltaFrac_vsfracReads_nDE_2seeds.pdf'), width = 6, height = 6, useDingbats = F)

iCard_ds_rpm20_deltaFracPer10pc_vsfracReads_nDE_2seeds <- ggplot() +
  geom_point(data = all_sig_res_ds_avg_rpm20_fracChanges, aes(frac_reads, delta_median_frac_nDE_per10pc)) +
  geom_line(data = all_sig_res_ds_avg_rpm20_fracChanges_meanDelta, aes(frac_reads, mean_delta_per10pc)) +
  ylim(c(0,0.15)) +
  theme_bw() +
  ylab('Fractional change in total DE genes\nNormalized to 10% increment\nBy change in downsampled reads') +
  xlab('Fraction of downsampled reads') +
  ggtitle('Fraction of detected DE genes downsampling analysis\nNormalized to 10% increment, Min RPM = 20, Adj. p-val <= 0.1\nTwo random seeds tested')
ggsave(iCard_ds_rpm20_deltaFracPer10pc_vsfracReads_nDE_2seeds, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/iCard_ds_rpm20_deltaFracPer10pc_vsfracReads_nDE_2seeds.pdf'), width = 6, height = 6, useDingbats = F)

write.csv(all_sig_res_ds_avg_rpm20_fracChanges, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/all_sig_res_ds_avg_rpm20_fracChanges.csv'))
write.csv(all_sig_res_ds_avg_rpm20_fracChanges_meanDelta, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/all_sig_res_ds_avg_rpm20_fracChanges_meanDelta.csv'))


iCard_ds_rpm20_sum_median <- ggplot() +
  geom_line(data = all_sig_res_ds_sum_rpm20 %>% filter(total_nDE > 100), aes(fracReads, frac_nDE, group = Product.Name)) +
  geom_line(data = all_sig_res_ds_avg_rpm20, aes(fracReads, median_frac_nDE), color = 'green', size = 2) +
  facet_wrap(~seed) +
  ylim(c(0,1.2)) +
  xlim(c(0,1)) +
  xlab('fraction of reads downsampled') +
  ylab('fraction of total diff. expr. genes detected') +
  ggtitle('Downsampling of iCard RNAtag-seq\nTwo random seeds\nOnly genes with mean RPM > 20, Padj cutoff = 0.1, Green = median') +
  geom_abline(slope = 1, intercept = 0, color = 'blue') +
  theme_bw()
ggsave(iCard_ds_rpm20_sum_median, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/iCard_downsampled_frac_nDE_twoSeeds_median.pdf'), width = 10, height = 5, useDingbats = F)


iCard_ds_rpmFilters_sum <- ggplot() +
  geom_line(data = all_sig_res_ds_rpmFilters %>% filter(total_nDE > 100), aes(fracReads, frac_nDE, group = Product.Name)) +
  geom_point(data = all_sig_res_ds_avg_rpmFilters, aes(fracReads, mean_frac_nDE), color = 'green', size = 3) +
  geom_errorbar(data = all_sig_res_ds_avg_rpmFilters, aes(fracReads, ymin = mean_frac_nDE - sem_frac_nDE, ymax = mean_frac_nDE + sem_frac_nDE), color = 'green', width = 0.07) +
  facet_grid(seed~rpmFilter) +
  ylim(c(0,1.2)) +
  xlim(c(0,1)) +
  xlab('fraction of reads downsampled') +
  ylab('fraction of total diff. expr. genes detected') +
  ggtitle('Downsampling of iCard RNAtag-seq\nTwo random seeds\nFour minimum RPM filters, Padj cutoff = 0.1, green = mean') +
  geom_abline(slope = 1, intercept = 0, color = 'blue') +
  theme_bw()
ggsave(iCard_ds_rpmFilters_sum, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/iCard_downsampled_frac_nDE_twoSeeds_rpmFilters.pdf'), width = 20, height = 10, useDingbats = F)

iCard_ds_rpmFilters_sum_median <- ggplot() +
  geom_line(data = all_sig_res_ds_rpmFilters %>% filter(total_nDE > 100), aes(fracReads, frac_nDE, group = Product.Name)) +
  geom_line(data = all_sig_res_ds_avg_rpmFilters, aes(fracReads, median_frac_nDE), color = 'green', size = 3) +
  facet_grid(seed~rpmFilter) +
  ylim(c(0,1.2)) +
  xlim(c(0,1)) +
  xlab('fraction of reads downsampled') +
  ylab('fraction of total diff. expr. genes detected') +
  ggtitle('Downsampling of iCard RNAtag-seq\nTwo random seeds\nFour minimum RPM filters, Padj cutoff = 0.1, green = median') +
  geom_abline(slope = 1, intercept = 0, color = 'blue') +
  theme_bw()
ggsave(iCard_ds_rpmFilters_sum_median, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/iCard_downsampled_frac_nDE_twoSeeds_median_rpmFilters.pdf'), width = 20, height = 10, useDingbats = F)

iCard_ds_rpmFilters_sum_smoothed <- ggplot() +
  geom_line(data = all_sig_res_ds_rpmFilters %>% filter(total_nDE > 100), aes(fracReads, frac_nDE, group = Product.Name)) +
  geom_smooth(data = all_sig_res_ds_rpmFilters %>% filter(total_nDE > 100), aes(fracReads, frac_nDE), color = 'green') +
  facet_grid(seed~rpmFilter) +
  ylim(c(0,1.2)) +
  xlim(c(0,1)) +
  xlab('fraction of reads downsampled') +
  ylab('fraction of total diff. expr. genes detected') +
  ggtitle('Downsampling of iCard RNAtag-seq\nTwo random seeds\nFour minimum RPM filters, Padj cutoff = 0.1, green = smoothed by LOESS') +
  geom_abline(slope = 1, intercept = 0, color = 'blue') +
  theme_bw()
ggsave(iCard_ds_rpmFilters_sum_smoothed, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/iCard_downsampled_frac_nDE_twoSeeds_smoothed_rpmFilters.pdf'), width = 20, height = 10, useDingbats = F)

# also do abs(LFC) cutoff tuning
lfc_list = c(0.5, 0.75, 1)
all_sig_res_lfcFilters_rpmFilters <- all_sig_res_ds_rpmFilters %>%
  mutate(lfc_cutoff = 0)
for (lfc_temp in lfc_list) {
  
  cat(paste0('Downsampling analysis of LFC cutoff of ', as.character(lfc_temp), '\n'))
  
  all_sig_res_ds_lfc0_5 <- all_sig_res_ds %>%
    filter(abs(log2FoldChange) >= lfc_temp) %>%
    mutate(lfc_cutoff = lfc_temp)
  
  all_sig_res_lfc0_5_temp_rpm20 <- all_sig_res_ds_lfc0_5 %>%
    filter(fracReads == 1, seed == seed_list[1]) %>%
    inner_join(temp_means_rpm20) %>%
    group_by(Product.Name) %>%
    summarise(total_nDE = length(Product.Name))
  all_sig_res_lfc0_5_temp_rpm10 <- all_sig_res_ds_lfc0_5 %>%
    filter(fracReads == 1, seed == seed_list[1]) %>%
    inner_join(temp_means_rpm10) %>%
    group_by(Product.Name) %>%
    summarise(total_nDE = length(Product.Name))
  all_sig_res_lfc0_5_temp_rpm5 <- all_sig_res_ds_lfc0_5 %>%
    filter(fracReads == 1, seed == seed_list[1]) %>%
    inner_join(temp_means_rpm5) %>%
    group_by(Product.Name) %>%
    summarise(total_nDE = length(Product.Name))
  all_sig_res_lfc0_5_temp_rpmNone <- all_sig_res_ds_lfc0_5 %>%
    filter(fracReads == 1, seed == seed_list[1]) %>%
    inner_join(temp_means_nofilt) %>%
    group_by(Product.Name) %>%
    summarise(total_nDE = length(Product.Name))
  
  all_sig_res_lfc0_5_ds_sum_rpm20 <- all_sig_res_ds_lfc0_5 %>%
    inner_join(temp_means_rpm20) %>%
    group_by(fracReads, seed, Product.Name) %>%
    summarise(nDE = length(fracReads)) %>%
    inner_join(all_sig_res_lfc0_5_temp_rpm20, by = 'Product.Name') %>%
    mutate(frac_nDE = nDE/total_nDE,
           rpmFilter = 20,
           lfc_cutoff = lfc_temp)
  
  all_sig_res_lfc0_5_ds_sum_rpm10 <- all_sig_res_ds_lfc0_5 %>%
    inner_join(temp_means_rpm10) %>%
    group_by(fracReads, seed, Product.Name) %>%
    summarise(nDE = length(fracReads)) %>%
    inner_join(all_sig_res_lfc0_5_temp_rpm10, by = 'Product.Name') %>%
    mutate(frac_nDE = nDE/total_nDE,
           rpmFilter = 10,
           lfc_cutoff = lfc_temp)
  
  all_sig_res_lfc0_5_ds_sum_rpm5 <- all_sig_res_ds_lfc0_5 %>%
    inner_join(temp_means_rpm5) %>%
    group_by(fracReads, seed, Product.Name) %>%
    summarise(nDE = length(fracReads)) %>%
    inner_join(all_sig_res_lfc0_5_temp_rpm5, by = 'Product.Name') %>%
    mutate(frac_nDE = nDE/total_nDE,
           rpmFilter = 5,
           lfc_cutoff = lfc_temp)
  
  all_sig_res_lfc0_5_ds_sum_rpmNone <- all_sig_res_ds_lfc0_5 %>%
    inner_join(temp_means_nofilt) %>%
    group_by(fracReads, seed, Product.Name) %>%
    summarise(nDE = length(fracReads)) %>%
    inner_join(all_sig_res_lfc0_5_temp_rpmNone, by = 'Product.Name') %>%
    mutate(frac_nDE = nDE/total_nDE,
           rpmFilter = 0,
           lfc_cutoff = lfc_temp)
  
  all_sig_res_lfc0_5_ds_rpmFilters <- bind_rows(all_sig_res_lfc0_5_ds_sum_rpm20, all_sig_res_lfc0_5_ds_sum_rpm10) %>%
    bind_rows(all_sig_res_lfc0_5_ds_sum_rpm5) %>%
    bind_rows(all_sig_res_lfc0_5_ds_sum_rpmNone)
  
  all_sig_res_lfcFilters_rpmFilters %<>% bind_rows(all_sig_res_lfc0_5_ds_rpmFilters)
  
}

iCard_ds_rpmFilters_sum_lfc0_5_smoothed <- ggplot() +
  geom_line(data = all_sig_res_lfcFilters_rpmFilters %>% filter(total_nDE > 100, lfc_cutoff == 0.5), aes(fracReads, frac_nDE, group = Product.Name)) +
  geom_smooth(data = all_sig_res_lfcFilters_rpmFilters %>% filter(total_nDE > 100, lfc_cutoff == 0.5), aes(fracReads, frac_nDE), color = 'green') +
  facet_grid(seed~rpmFilter) +
  ylim(c(0,1.2)) +
  xlim(c(0,1)) +
  xlab('fraction of reads downsampled') +
  ylab('fraction of total diff. expr. genes detected') +
  ggtitle('Downsampling of iCard RNAtag-seq\nTwo random seeds\nFour minimum RPM filters, Padj cutoff = 0.1, LFC cutoff = 0.5, green = smoothed by LOESS') +
  geom_abline(slope = 1, intercept = 0, color = 'blue') +
  theme_bw()
ggsave(iCard_ds_rpmFilters_sum_lfc0_5_smoothed, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/iCard_downsampled_frac_nDE_twoSeeds_smoothed_lfc0_5_rpmFilters.pdf'), width = 20, height = 10, useDingbats = F)

iCard_ds_rpmFilters_sum_lfc0_median <- ggplot() +
  geom_line(data = all_sig_res_lfcFilters_rpmFilters %>% filter(total_nDE > 100, lfc_cutoff == 0), aes(fracReads, frac_nDE, group = Product.Name)) +
  geom_line(data = all_sig_res_lfcFilters_rpmFilters %>% filter(lfc_cutoff == 0) %>% group_by(seed, rpmFilter, fracReads, lfc_cutoff) %>% summarise(medianfrac = median(frac_nDE)), 
            aes(fracReads, medianfrac), color = 'green', size = 2) +
  facet_grid(seed~rpmFilter) +
  ylim(c(0,1.2)) +
  xlim(c(0,1)) +
  xlab('fraction of reads downsampled') +
  ylab('fraction of total diff. expr. genes detected') +
  ggtitle('Downsampling of iCard RNAtag-seq\nTwo random seeds\nFour minimum RPM filters, Padj cutoff = 0.1, LFC cutoff = None, green = median') +
  geom_abline(slope = 1, intercept = 0, color = 'blue') +
  theme_bw()
ggsave(iCard_ds_rpmFilters_sum_lfc0_median, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/iCard_downsampled_frac_nDE_twoSeeds_median_lfc0_rpmFilters.pdf'), width = 20, height = 10, useDingbats = F)


iCard_ds_rpmFilters_sum_lfc0_5_median <- ggplot() +
  geom_line(data = all_sig_res_lfcFilters_rpmFilters %>% filter(total_nDE > 100, lfc_cutoff == 0.5), aes(fracReads, frac_nDE, group = Product.Name)) +
  geom_line(data = all_sig_res_lfcFilters_rpmFilters %>% filter(lfc_cutoff == 0.5) %>% group_by(seed, rpmFilter, fracReads, lfc_cutoff) %>% summarise(medianfrac = median(frac_nDE)), 
            aes(fracReads, medianfrac), color = 'green', size = 2) +
  facet_grid(seed~rpmFilter) +
  ylim(c(0,1.2)) +
  xlim(c(0,1)) +
  xlab('fraction of reads downsampled') +
  ylab('fraction of total diff. expr. genes detected') +
  ggtitle('Downsampling of iCard RNAtag-seq\nTwo random seeds\nFour minimum RPM filters, Padj cutoff = 0.1, LFC cutoff = 0.5, green = median') +
  geom_abline(slope = 1, intercept = 0, color = 'blue') +
  theme_bw()
ggsave(iCard_ds_rpmFilters_sum_lfc0_5_median, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/iCard_downsampled_frac_nDE_twoSeeds_median_lfc0_5_rpmFilters.pdf'), width = 20, height = 10, useDingbats = F)

iCard_ds_rpmFilters_sum_lfc0_75_median <- ggplot() +
  geom_line(data = all_sig_res_lfcFilters_rpmFilters %>% filter(total_nDE > 100, lfc_cutoff == 0.75), aes(fracReads, frac_nDE, group = Product.Name)) +
  geom_line(data = all_sig_res_lfcFilters_rpmFilters %>% filter(lfc_cutoff == 0.75) %>% group_by(seed, rpmFilter, fracReads, lfc_cutoff) %>% summarise(medianfrac = median(frac_nDE)), 
            aes(fracReads, medianfrac), color = 'green', size = 2) +
  facet_grid(seed~rpmFilter) +
  ylim(c(0,1.2)) +
  xlim(c(0,1)) +
  xlab('fraction of reads downsampled') +
  ylab('fraction of total diff. expr. genes detected') +
  ggtitle('Downsampling of iCard RNAtag-seq\nTwo random seeds\nFour minimum RPM filters, Padj cutoff = 0.1, LFC cutoff = 0.75, green = median') +
  geom_abline(slope = 1, intercept = 0, color = 'blue') +
  theme_bw()
ggsave(iCard_ds_rpmFilters_sum_lfc0_75_median, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/iCard_downsampled_frac_nDE_twoSeeds_median_lfc0_75_rpmFilters.pdf'), width = 20, height = 10, useDingbats = F)

iCard_ds_rpmFilters_sum_lfc1_median <- ggplot() +
  geom_line(data = all_sig_res_lfcFilters_rpmFilters %>% filter(total_nDE > 100, lfc_cutoff == 1), aes(fracReads, frac_nDE, group = Product.Name)) +
  geom_line(data = all_sig_res_lfcFilters_rpmFilters %>% filter(lfc_cutoff == 1) %>% group_by(seed, rpmFilter, fracReads, lfc_cutoff) %>% summarise(medianfrac = median(frac_nDE)), 
            aes(fracReads, medianfrac), color = 'green', size = 2) +
  facet_grid(seed~rpmFilter) +
  ylim(c(0,1.2)) +
  xlim(c(0,1)) +
  xlab('fraction of reads downsampled') +
  ylab('fraction of total diff. expr. genes detected') +
  ggtitle('Downsampling of iCard RNAtag-seq\nTwo random seeds\nFour minimum RPM filters, Padj cutoff = 0.1, LFC cutoff = 1, green = median') +
  geom_abline(slope = 1, intercept = 0, color = 'blue') +
  theme_bw()
ggsave(iCard_ds_rpmFilters_sum_lfc1_median, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/iCard_downsampled_frac_nDE_twoSeeds_median_lfc1_rpmFilters.pdf'), width = 20, height = 10, useDingbats = F)


# compare across LFC cutoffs
iCard_ds_lfcFilters_sum_rpm0_median <- ggplot() +
  geom_line(data = all_sig_res_lfcFilters_rpmFilters %>% filter(total_nDE > 100, rpmFilter == 0), aes(fracReads, frac_nDE, group = Product.Name)) +
  geom_line(data = all_sig_res_lfcFilters_rpmFilters %>% filter(rpmFilter == 0) %>% group_by(seed, rpmFilter, fracReads, lfc_cutoff) %>% summarise(medianfrac = median(frac_nDE)), 
            aes(fracReads, medianfrac), color = 'green', size = 2) +
  facet_grid(seed~lfc_cutoff) +
  ylim(c(0,1.2)) +
  xlim(c(0,1)) +
  xlab('fraction of reads downsampled') +
  ylab('fraction of total diff. expr. genes detected') +
  ggtitle('Downsampling of iCard RNAtag-seq\nTwo random seeds\nFour minimum LFC filters, Padj cutoff = 0.1, RPM cutoff = 0, green = median') +
  geom_abline(slope = 1, intercept = 0, color = 'blue') +
  theme_bw()
ggsave(iCard_ds_lfcFilters_sum_rpm0_median, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/iCard_downsampled_frac_nDE_twoSeeds_median_rpm0_lfcFilters.pdf'), width = 20, height = 10, useDingbats = F)

iCard_ds_lfcFilters_sum_rpm5_median <- ggplot() +
  geom_line(data = all_sig_res_lfcFilters_rpmFilters %>% filter(total_nDE > 100, rpmFilter == 5), aes(fracReads, frac_nDE, group = Product.Name)) +
  geom_line(data = all_sig_res_lfcFilters_rpmFilters %>% filter(rpmFilter == 5) %>% group_by(seed, rpmFilter, fracReads, lfc_cutoff) %>% summarise(medianfrac = median(frac_nDE)), 
            aes(fracReads, medianfrac), color = 'green', size = 2) +
  facet_grid(seed~lfc_cutoff) +
  ylim(c(0,1.2)) +
  xlim(c(0,1)) +
  xlab('fraction of reads downsampled') +
  ylab('fraction of total diff. expr. genes detected') +
  ggtitle('Downsampling of iCard RNAtag-seq\nTwo random seeds\nFour minimum LFC filters, Padj cutoff = 0.1, RPM cutoff = 5, green = median') +
  geom_abline(slope = 1, intercept = 0, color = 'blue') +
  theme_bw()
ggsave(iCard_ds_lfcFilters_sum_rpm5_median, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/iCard_downsampled_frac_nDE_twoSeeds_median_rpm5_lfcFilters.pdf'), width = 20, height = 10, useDingbats = F)

iCard_ds_lfcFilters_sum_rpm10_median <- ggplot() +
  geom_line(data = all_sig_res_lfcFilters_rpmFilters %>% filter(total_nDE > 100, rpmFilter == 10), aes(fracReads, frac_nDE, group = Product.Name)) +
  geom_line(data = all_sig_res_lfcFilters_rpmFilters %>% filter(rpmFilter == 10) %>% group_by(seed, rpmFilter, fracReads, lfc_cutoff) %>% summarise(medianfrac = median(frac_nDE)), 
            aes(fracReads, medianfrac), color = 'green', size = 2) +
  facet_grid(seed~lfc_cutoff) +
  ylim(c(0,1.2)) +
  xlim(c(0,1)) +
  xlab('fraction of reads downsampled') +
  ylab('fraction of total diff. expr. genes detected') +
  ggtitle('Downsampling of iCard RNAtag-seq\nTwo random seeds\nFour minimum LFC filters, Padj cutoff = 0.1, RPM cutoff = 10, green = median') +
  geom_abline(slope = 1, intercept = 0, color = 'blue') +
  theme_bw()
ggsave(iCard_ds_lfcFilters_sum_rpm10_median, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/iCard_downsampled_frac_nDE_twoSeeds_median_rpm10_lfcFilters.pdf'), width = 20, height = 10, useDingbats = F)


iCard_ds_lfcFilters_sum_rpm20_median <- ggplot() +
  geom_line(data = all_sig_res_lfcFilters_rpmFilters %>% filter(total_nDE > 100, rpmFilter == 20), aes(fracReads, frac_nDE, group = Product.Name)) +
  geom_line(data = all_sig_res_lfcFilters_rpmFilters %>% filter(rpmFilter == 20) %>% group_by(seed, rpmFilter, fracReads, lfc_cutoff) %>% summarise(medianfrac = median(frac_nDE)), 
            aes(fracReads, medianfrac), color = 'green', size = 2) +
  facet_grid(seed~lfc_cutoff) +
  ylim(c(0,1.2)) +
  xlim(c(0,1)) +
  xlab('fraction of reads downsampled') +
  ylab('fraction of total diff. expr. genes detected') +
  ggtitle('Downsampling of iCard RNAtag-seq\nTwo random seeds\nFour minimum LFC filters, Padj cutoff = 0.1, RPM cutoff = 20, green = median') +
  geom_abline(slope = 1, intercept = 0, color = 'blue') +
  theme_bw()
ggsave(iCard_ds_lfcFilters_sum_rpm20_median, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/iCard_downsampled_frac_nDE_twoSeeds_median_rpm20_lfcFilters.pdf'), width = 20, height = 10, useDingbats = F)


# compare 1 drug across 2 seeds, when frac_nDE goes above 1
# just to show that random downsampling can shuffle genes differently, leading to semi-unintuitive increase in nDE with fewer total reads.
# pick d29 given a LFC cutoff of 0.5, just for example. For seed 2380 it does go above 1 at fracReads 0.7, but not for seed 53250
iCard_d29_lfc0_5_rpm20_sum <- all_sig_res_lfcFilters_rpmFilters %>% filter(rpmFilter == 20, lfc_cutoff == 0.5, Product.Name == 'd29')

iCard_d29_lfc0_5_rpm20_sig_res_DEseq <- all_sig_res_ds %>% filter(Product.Name == 'd29') %>% inner_join(temp_means_rpm20) %>% filter(meanRPM > 20)

iCard_d29_lfc0_5_rpm20_sig_res_DEseq_wide <- iCard_d29_lfc0_5_rpm20_sig_res_DEseq %>%
  dplyr::select(seed, gene_id, meanRPM, log2FoldChange, fracReads) %>%
  group_by(seed, gene_id, meanRPM) %>%
  spread(fracReads, log2FoldChange)

# consider ENSG00000008988
all_res_ds <- list()
for (seedID in seed_list){
  for (fracReads in downsamp_list) {
    
    all_genes <- as_tibble(read.table(paste0(graphDir, "/differentialExpression/allSamples/iCards/downsampled/iCard_readFilt_manualFilt_DESeqResults_", as.character(fracReads*100), "percent_seed_",as.character(seedID),".txt"), stringsAsFactors = F, header = T, sep = '\t')) %>%
      mutate(fracReads = fracReads,
             seed = seedID) 
    
    if(is.null(dim(all_res_ds))) {
      all_res_ds <- all_genes
    } else {
      all_res_ds %<>% bind_rows(all_genes)
    }
    
  }
}
all_res_full <- as_tibble(read.table(paste0(graphDir, "/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_DESeqResults.txt"), stringsAsFactors = F, header = T, sep = '\t')) %>%
  mutate(fracReads = 1)

all_res_ds %<>%
  bind_rows(all_res_full %>% mutate(seed = seed_list[1])) %>%
  bind_rows(all_res_full %>% mutate(seed = seed_list[2]))

all_res_ds_ENSG00000008988_d29 <- all_res_ds %>%
  filter(Product.Name == 'd29', gene_id == 'ENSG00000008988') %>%
  dplyr::select(seed, gene_id, Product.Name, fracReads, baseMean, log2FoldChange, padj) %>%
  arrange(seed)

write.table(all_res_ds_ENSG00000008988_d29, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/example_downsample_padjChange_ENSG00000008988_d29.txt'), quote = F, row.names = F, sep = '\t')

# downsampled reads for this gene in this drug
# switched to TPM
inds_ENSG00000008988_d29 <- colData(allset)$cellType == 'iCard' & colData(allset)$condition %in% c('d29', 'DMSOonly-control')
temp_tpms_ds_ENSG00000008988_d29_ctl <- as.data.frame(assay(allset[,inds_ENSG00000008988_d29]))
temp_tpms_ds_ENSG00000008988_d29_ctl$gene_id <- rownames(temp_tpms_ds_ENSG00000008988_d29_ctl)
rownames(temp_tpms_ds_ENSG00000008988_d29_ctl) <- NULL

tpms_ds_ENSG00000008988_d29_ctl <- temp_tpms_ds_ENSG00000008988_d29_ctl %>% as_tibble() %>%
  gather('sampleID', 'tpm',1:(ncol(temp_tpms_ds_ENSG00000008988_d29_ctl)-1)) %>%
  mutate(fracReads = 1,
         seed = seed_list[1]) %>%
  filter(gene_id == 'ENSG00000008988')
tpms_ds_ENSG00000008988_d29_ctl %<>% bind_rows(temp_tpms_ds_ENSG00000008988_d29_ctl %>% as_tibble() %>%
                                                 gather('sampleID', 'tpm',1:(ncol(temp_tpms_ds_ENSG00000008988_d29_ctl)-1)) %>%
                                                 mutate(fracReads = 1,
                                                        seed = seed_list[2])%>%
                                                 filter(gene_id == 'ENSG00000008988'))

sampleInfo <- as.data.frame(colData(allset))
sampleInfo$sampleID <- rownames(sampleInfo)
rownames(sampleInfo) <- NULL
for (seedID in seed_list) {
  cat(paste0('Working on seed ', as.character(seedID),'\n'))
  for (fracReads in downsamp_list) {
    cat(paste0('Working on fracReads ', as.character(fracReads), '\n'))
    sePath = paste0(graphDir, 'differentialExpression/allSamples/iCards/downsampled')
    iCard_counts_readFilt_manualFilt_temp = readRDS(paste0(sePath, '/se_mappedCounts_',as.character(fracReads*100), "percent_seed_", as.character(seedID),'.rds'))
    geneLengths <- as.data.frame(rowData(iCard_counts_readFilt_manualFilt_temp)) %>% dplyr::select(-GeneSymbol)
    rownames(geneLengths) <- NULL
    inds_ENSG00000008988_d29 <- colData(iCard_counts_readFilt_manualFilt_temp)$cellType == 'iCard' & colData(iCard_counts_readFilt_manualFilt_temp)$condition %in% c('d29', 'DMSOonly-control')
    
    temp_counts_ds_ENSG00000008988_d29_ctl <- as.data.frame(assay(iCard_counts_readFilt_manualFilt_temp[,inds_ENSG00000008988_d29]))
    temp_counts_ds_ENSG00000008988_d29_ctl$gene_id <- rownames(temp_counts_ds_ENSG00000008988_d29_ctl)
    rownames(temp_counts_ds_ENSG00000008988_d29_ctl) <- NULL
    temp_counts_ds_ENSG00000008988_d29_ctl %<>% as_tibble() %>%
      gather('combID', 'counts',1:(ncol(temp_counts_ds_ENSG00000008988_d29_ctl)-1)) %>%
      separate(combID, into = c('experiment', 'sampleID', 'tagID'), sep = '\\_')
    
    temp_tpms <- generateTPMfromCounts(temp_counts_ds_ENSG00000008988_d29_ctl, as_tibble(geneLengths)) %>%
      mutate(fracReads = fracReads,
             seed = seedID) %>%
      unite(sampleID, experiment:tagID, sep = '_') %>%
      dplyr::select(-c(length, counts))
    
    tpms_ds_ENSG00000008988_d29_ctl %<>% bind_rows(temp_tpms)
    
  }
}
tpms_ds_ENSG00000008988_d29_ctl %<>%
  filter(gene_id == 'ENSG00000008988') %>%
  inner_join(sampleInfo %>% dplyr::select(sampleID, condition) %>% dplyr::rename(Product.Name = condition))

tpms_ds_ENSG00000008988_d29_colorByPlate <- ggplot() +
  geom_histogram(data = tpms_ds_ENSG00000008988_d29_ctl %>% filter(Product.Name == 'DMSOonly-control') %>% separate(sampleID, into = c('plate', 'pool', 'tag')), aes(tpm, fill = plate), binwidth = 30) +
  geom_vline(data = tpms_ds_ENSG00000008988_d29_ctl %>% filter(Product.Name == 'd29'), aes(xintercept = tpm), color = 'orange') +
  geom_text(data = all_res_ds_ENSG00000008988_d29, aes(400, 15, label = paste0('LFC ', as.character(round(log2FoldChange, 3)), '\nPadj ', as.character(round(padj,3))))) +
  facet_grid(fracReads ~ seed) +
  ggtitle('ENSG00000008988 TPM values in d29\nControl values in histogram\nd29 values in orange lines\nDownsampled to 4 fractions with 2 seeds')
ggsave(tpms_ds_ENSG00000008988_d29_colorByPlate, filename = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/ENSG00000008988_d29_TPMs_histogram_colorByPlate.pdf'), width = 8, height = 8, useDingbats = F)
tpms_ds_ENSG00000008988_d29 <- ggplot() +
  geom_histogram(data = tpms_ds_ENSG00000008988_d29_ctl %>% filter(Product.Name == 'DMSOonly-control') %>% separate(sampleID, into = c('plate', 'pool', 'tag')), aes(tpm), binwidth = 30) +
  geom_vline(data = tpms_ds_ENSG00000008988_d29_ctl %>% filter(Product.Name == 'd29'), aes(xintercept = tpm), color = 'orange') +
  geom_text(data = all_res_ds_ENSG00000008988_d29, aes(400, 15, label = paste0('LFC ', as.character(round(log2FoldChange, 3)), '\nPadj ', as.character(round(padj,3))))) +
  facet_grid(fracReads ~ seed) +
  ggtitle('ENSG00000008988 TPM values in d29\nControl values in histogram\nd29 values in orange lines\nDownsampled to 4 fractions with 2 seeds')
ggsave(tpms_ds_ENSG00000008988_d29, filename = paste0(graphDir, '/differentialExpression/allSamples/iCards/downsampled/ENSG00000008988_d29_TPMs_histogram.pdf'), width = 8, height = 8, useDingbats = F)

print('Done with downsampled differential expression analysis.')
print('Starting representative gene-level (FOXJ3, MEF2C) perturbability analysis...') ### UP TO HERE
## consider FOXJ3 and MEF2C for main figures, d21, d44, d53
genes <- c('FOXJ3', 'MEF2C')
genes_tbl <- geneIDtoGeneSymbol %>% filter(GeneSymbol %in% genes)
drugs <- c('DMSOonly-control', 'd21', 'd44', 'd26')
inds_temp <- colData(all_tpms_geneFilt_readFilt_manualFilt)$cellType == 'iCard' & colData(all_tpms_geneFilt_readFilt_manualFilt)$condition %in% drugs
tpms_temp <- as.data.frame(assay(all_tpms_geneFilt_readFilt_manualFilt[,inds_temp]))
tpms_temp$gene_id <- rownames(tpms_temp)
rownames(tpms_temp) <- NULL

sampleData_temp <- colData(all_tpms_geneFilt_readFilt_manualFilt)[inds_temp, c('Product.Name', 'condition')]
sampleData_temp$sampleID <- rownames(sampleData_temp)
rownames(sampleData_temp) <- NULL
sampleData_temp <- as_tibble(sampleData_temp) 

temp_DE_res <- all_res_full %>%
  inner_join(genes_tbl, by = 'gene_id') %>%
  filter(Product.Name %in% drugs) %>%
  mutate(Product_col = Product.Name)

tpms_temp_tall <- tpms_temp %>% as_tibble() %>%
  inner_join(genes_tbl) %>%
  gather('sampleID', 'tpm', 1:(ncol(tpms_temp) - 1)) %>%
  inner_join(sampleData_temp, by = 'sampleID')

tpms_temp_ctl <- list()
for(di in drugs) {
  
  if(di != 'DMSOonly-control'){
    
    tpms_temp_di <- tpms_temp_tall %>%
      filter(condition == 'DMSOonly-control') %>%
      mutate(Product_col = di,
             isPerturbed = F)
    
    if(is.null(dim(tpms_temp_ctl))){
      tpms_temp_ctl <- tpms_temp_di
    } else {
      tpms_temp_ctl %<>% bind_rows(tpms_temp_di)
    }
  }
}

tpms_temp_drugs <- tpms_temp_tall %>%
  filter(condition != 'DMSOonly-control') %>%
  mutate(Product_col = condition,
         isPerturbed = T)

tpms_temp_bound <- bind_rows(tpms_temp_drugs, tpms_temp_ctl)

tpms_temp_sums_ctl <- tpms_temp_ctl %>% group_by(GeneSymbol, Product_col) %>% summarise(meanTPM = mean(tpm), sdTPM = sd(tpm)) %>% mutate(isPerturbed = F)
tpms_temp_sums_drugs <- tpms_temp_drugs %>% group_by(GeneSymbol, Product_col) %>% summarise(meanTPM = mean(tpm)) %>% mutate(isPerturbed = T)

tpms_temp_sums_bound <- bind_rows(tpms_temp_sums_ctl, tpms_temp_sums_drugs)


tpms_FOXJ3MEF2C_someDrugs <- ggplot() +
  geom_histogram(data = tpms_temp_ctl, aes(tpm), binwidth = 20) +
  geom_vline(data = tpms_temp_drugs, aes(xintercept = tpm), color = 'orange') +
  geom_text(data = temp_DE_res, aes(250, 9, label = paste0('LFC ', as.character(round(log2FoldChange, 3)), '\nPadj ', as.character(round(padj,3))))) +
  facet_grid(GeneSymbol ~ Product_col) +
  ggtitle('MEF2C and FOXJ3 values in 3 drugs') +
  theme_bw()
ggsave(tpms_FOXJ3MEF2C_someDrugs, filename = paste0(graphDir, '/differentialExpression/allSamples/iCards/tpms_FOXJ3MEF2C_someDrugs.pdf'), width = 12, height = 8, useDingbats = F)

tpms_FOXJ3MEF2C_someDrugs_bees <- ggplot() +
  geom_jitter(data = tpms_temp_ctl, aes(Product_col, tpm), height = 0, width = 0.1) +
  geom_point(data = tpms_temp_drugs, aes(Product_col, tpm), color = 'orange', size = 2) +
  geom_bar(data = tpms_temp_sums_ctl, aes(Product_col, meanTPM), alpha = 0, color = 'black', width = 0.4, stat = 'identity') +
  geom_bar(data = tpms_temp_sums_drugs, aes(Product_col, meanTPM), alpha = 0, color = 'orange', width = 0.4, stat = 'identity') +
  geom_text(data = temp_DE_res, aes(Product_col, 320, label = paste0('LFC ', as.character(round(log2FoldChange, 3)), '\nPadj ', as.character(round(padj,3))))) +
  facet_grid(GeneSymbol ~ .) +
  ggtitle('MEF2C and FOXJ3 values in 3 drugs') +
  theme_bw() +
  ylim(c(0,350))
ggsave(tpms_FOXJ3MEF2C_someDrugs_bees, filename = paste0(graphDir, '/differentialExpression/allSamples/iCards/tpms_FOXJ3MEF2C_someDrugs_bees.pdf'), width = 12, height = 8, useDingbats = F)

set.seed(42639)
tpms_FOXJ3MEF2C_someDrugs_bees2 <- ggplot() +
  geom_point(data = tpms_temp_bound, aes(Product_col, tpm, color = isPerturbed, group = as.character(isPerturbed)), position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0, dodge.width = 0.35)) +
  geom_errorbar(data = tpms_temp_sums_bound, aes(Product_col, ymin = meanTPM, ymax = meanTPM, color = as.character(isPerturbed)), width = 0.4, position = position_dodge2(width = 0.9), alpha = 0.8, stat = 'identity') +
  scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'orange')) +
  geom_text(data = temp_DE_res, aes(Product_col, 320, label = paste0('LFC ', as.character(round(log2FoldChange, 3)), '\nPadj ', as.character(round(padj,3))))) +
  facet_grid(GeneSymbol ~ .) +
  ggtitle('MEF2C and FOXJ3 values in 3 drugs') +
  theme_bw() +
  ylim(c(0,350))
ggsave(tpms_FOXJ3MEF2C_someDrugs_bees2, filename = paste0(graphDir, '/differentialExpression/allSamples/iCards/tpms_FOXJ3MEF2C_someDrugs_bees2.pdf'), width = 12, height = 8, useDingbats = F)


# repeat for DESeq2 normalized counts
dds_iCard <- readRDS(paste0(graphDir, '/differentialExpression/allSamples/iCards/iCards_dds_se.rds'))

drugs2_tbl <- all_sig_res_full %>% inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>% filter(GeneSymbol == 'MEF2C') 
drugs2 <- c(drugs2_tbl$Product.Name, 'DMSOonly-control')

inds_temp_desnc <- colData(dds_iCard)$cellType == 'iCard' & colData(dds_iCard)$condition %in% drugs
desnc_temp <- as.data.frame(counts(dds_iCard, normalized = T)[,inds_temp_desnc])
desnc_temp$gene_id <- rownames(desnc_temp)
rownames(desnc_temp) <- NULL

sampleData_temp_desnc <- colData(dds_iCard)[inds_temp_desnc, c('Product.Name', 'condition')]
sampleData_temp_desnc$sampleID <- rownames(sampleData_temp_desnc)
rownames(sampleData_temp_desnc) <- NULL
sampleData_temp_desnc <- as_tibble(sampleData_temp_desnc) 

desnc_temp_tall <- desnc_temp %>% as_tibble() %>%
  inner_join(genes_tbl) %>%
  gather('sampleID', 'norm_counts', 1:(ncol(desnc_temp) - 1)) %>%
  inner_join(sampleData_temp_desnc, by = 'sampleID')

desnc_temp_ctl <- list()
for(di in drugs) {
  
  if(di != 'DMSOonly-control'){
    
    desnc_temp_di <- desnc_temp_tall %>%
      filter(condition == 'DMSOonly-control') %>%
      mutate(Product_col = di,
             isPerturbed = F)
    
    if(is.null(dim(desnc_temp_ctl))){
      desnc_temp_ctl <- desnc_temp_di
    } else {
      desnc_temp_ctl %<>% bind_rows(desnc_temp_di)
    }
  }
}

desnc_temp_drugs <- desnc_temp_tall %>%
  filter(condition != 'DMSOonly-control') %>%
  mutate(Product_col = condition,
         isPerturbed = T)

desnc_temp_bound <- bind_rows(desnc_temp_drugs, desnc_temp_ctl)

desnc_temp_sums_ctl <- desnc_temp_ctl %>% group_by(GeneSymbol, Product_col) %>% summarise(meanNC = mean(norm_counts), sdNC = sd(norm_counts)) %>% mutate(isPerturbed = F)
desnc_temp_sums_drugs <- desnc_temp_drugs %>% group_by(GeneSymbol, Product_col) %>% summarise(meanNC = mean(norm_counts)) %>% mutate(isPerturbed = T)

desnc_temp_sums_bound <- bind_rows(desnc_temp_sums_ctl, desnc_temp_sums_drugs)

set.seed(42639)
DESeq2NC_FOXJ3MEF2C_someDrugs_bees2 <- ggplot() +
  geom_point(data = desnc_temp_bound, aes(Product_col, norm_counts, color = isPerturbed, group = as.character(isPerturbed)), position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0, dodge.width = 0.35)) +
  geom_errorbar(data = desnc_temp_sums_bound, aes(Product_col, ymin = meanNC, ymax = meanNC, color = as.character(isPerturbed)), width = 0.4, position = position_dodge2(width = 0.9), alpha = 0.8, stat = 'identity') +
  scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'orange')) +
  geom_text(data = temp_DE_res, aes(Product_col, 320, label = paste0('LFC ', as.character(round(log2FoldChange, 3)), '\nPadj ', as.character(round(padj,3))))) +
  facet_grid(GeneSymbol ~ .) +
  ggtitle('MEF2C and FOXJ3 DESeq2 normalized counts values in 3 drugs') +
  theme_bw() +
  ylim(c(0,350))
ggsave(DESeq2NC_FOXJ3MEF2C_someDrugs_bees2, filename = paste0(graphDir, '/differentialExpression/allSamples/iCards/DESeq2NC_FOXJ3MEF2C_someDrugs_bees2.pdf'), width = 12, height = 8, useDingbats = F)

set.seed(42639)
DESeq2NC_FOXJ3MEF2C_someDrugs_bees2_forfig <- ggplot() +
  geom_point(data = desnc_temp_bound, aes(Product_col, norm_counts, color = isPerturbed, group = as.character(isPerturbed)), position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0, dodge.width = 0.35)) +
  geom_errorbar(data = desnc_temp_sums_bound, aes(Product_col, ymin = meanNC, ymax = meanNC, color = as.character(isPerturbed)), width = 0.4, position = position_dodge2(width = 0.9), alpha = 0.8, stat = 'identity') +
  scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'orange')) +
  # geom_text(data = temp_DE_res, aes(Product_col, 320, label = paste0('LFC ', as.character(round(log2FoldChange, 3)), '\nPadj ', as.character(round(padj,3))))) +
  facet_grid(GeneSymbol ~ .) +
  # ggtitle('MEF2C and FOXJ3 DESeq2 normalized counts values in 3 drugs') +
  theme_bw() +
  theme(legend.position = 'none') +
  ylim(c(0,200))
ggsave(DESeq2NC_FOXJ3MEF2C_someDrugs_bees2_forfig, filename = paste0(graphDir, '/differentialExpression/allSamples/iCards/DESeq2NC_FOXJ3MEF2C_someDrugs_bees2_forfig.pdf'), width = 3, height = 2, useDingbats = F)

desnc_temp_bound_noDodge <- desnc_temp_bound %>%
  filter(isPerturbed == T) %>%
  dplyr::select(-condition) %>%
  bind_rows(desnc_temp_bound %>% filter(isPerturbed == F) %>% dplyr::select(-condition) %>% mutate(Product_col = 'control', Product.Name = 'DMSO only'))
desnc_temp_sums_bound_noDodge <- desnc_temp_sums_bound %>%
  filter(isPerturbed == T) %>%
  bind_rows(desnc_temp_sums_bound %>% filter(isPerturbed == F) %>% mutate(Product_col = 'control') %>% unique())
  
set.seed(42639)
DESeq2NC_FOXJ3MEF2C_someDrugs_bees3_forfig <- ggplot() +
  geom_boxplot(data = desnc_temp_bound_noDodge %>% filter(isPerturbed == F), aes(Product_col, norm_counts), outlier.shape = NA) +
  geom_jitter(data = desnc_temp_bound_noDodge, aes(Product_col, norm_counts, color = isPerturbed, alpha = isPerturbed), width = 0.1, height = 0) +
  geom_errorbar(data = desnc_temp_sums_bound_noDodge %>% filter(isPerturbed == T), aes(Product_col, ymin = meanNC, ymax = meanNC, color = as.character(isPerturbed)), width = 0.7, alpha = 0.8, stat = 'identity') +
  scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'chocolate2')) +
  scale_alpha_manual(values = c('FALSE' = 0.2, 'TRUE' = 1)) +
  facet_grid(GeneSymbol ~ .) +
  # ggtitle('MEF2C and FOXJ3 DESeq2 normalized counts values in 3 drugs') +
  theme_classic() +
  theme(legend.position = 'none') +
  ylim(c(0,200))
ggsave(DESeq2NC_FOXJ3MEF2C_someDrugs_bees3_forfig, filename = paste0(graphDir, '/differentialExpression/allSamples/iCards/DESeq2NC_FOXJ3MEF2C_someDrugs_bees3_forfig.pdf'), width = 3, height = 2, useDingbats = F)


# more drugs
inds_temp_desnc2 <- colData(dds_iCard)$cellType == 'iCard' & colData(dds_iCard)$condition %in% drugs2
desnc_temp2 <- as.data.frame(counts(dds_iCard, normalized = T)[,inds_temp_desnc2])
desnc_temp2$gene_id <- rownames(desnc_temp2)
rownames(desnc_temp2) <- NULL

sampleData_temp_desnc2 <- colData(dds_iCard)[inds_temp_desnc2, c('Product.Name', 'condition')]
sampleData_temp_desnc2$sampleID <- rownames(sampleData_temp_desnc2)
rownames(sampleData_temp_desnc2) <- NULL
sampleData_temp_desnc2 <- as_tibble(sampleData_temp_desnc2) 

desnc_temp_tall2 <- desnc_temp2 %>% as_tibble() %>%
  inner_join(genes_tbl) %>%
  gather('sampleID', 'norm_counts', 1:(ncol(desnc_temp2) - 1)) %>%
  inner_join(sampleData_temp_desnc2, by = 'sampleID')

desnc_temp_ctl2 <- list()
for(di in drugs2) {
  
  if(di != 'DMSOonly-control'){
    
    desnc_temp_di <- desnc_temp_tall2 %>%
      filter(condition == 'DMSOonly-control') %>%
      mutate(Product_col = di,
             isPerturbed = F)
    
    if(is.null(dim(desnc_temp_ctl2))){
      desnc_temp_ctl2 <- desnc_temp_di
    } else {
      desnc_temp_ctl2 %<>% bind_rows(desnc_temp_di)
    }
  }
}

desnc_temp_drugs2 <- desnc_temp_tall2 %>%
  filter(condition != 'DMSOonly-control') %>%
  mutate(Product_col = condition,
         isPerturbed = T)

desnc_temp_bound2 <- bind_rows(desnc_temp_drugs2, desnc_temp_ctl2)

desnc_temp_sums_ctl2 <- desnc_temp_ctl2 %>% group_by(GeneSymbol, Product_col) %>% summarise(meanNC = mean(norm_counts), sdNC = sd(norm_counts)) %>% mutate(isPerturbed = F)
desnc_temp_sums_drugs2 <- desnc_temp_drugs2 %>% group_by(GeneSymbol, Product_col) %>% summarise(meanNC = mean(norm_counts)) %>% mutate(isPerturbed = T)

desnc_temp_sums_bound2 <- bind_rows(desnc_temp_sums_ctl2, desnc_temp_sums_drugs2)

temp_DE_res2 <- all_res_full %>%
  inner_join(genes_tbl, by = 'gene_id') %>%
  filter(Product.Name %in% drugs2) %>%
  mutate(Product_col = Product.Name)


set.seed(42639)
DESeq2NC_FOXJ3MEF2C_moreDrugs_bees2 <- ggplot() +
  geom_point(data = desnc_temp_bound2, aes(Product_col, norm_counts, color = isPerturbed, group = as.character(isPerturbed)), position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.35)) +
  geom_errorbar(data = desnc_temp_sums_bound2, aes(Product_col, ymin = meanNC, ymax = meanNC, color = as.character(isPerturbed)), width = 0.8, position = position_dodge2(width = 0.9), alpha = 0.8, stat = 'identity') +
  scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'orange')) +
  geom_text(data = temp_DE_res2, aes(Product_col, 320, label = paste0('LFC ', as.character(round(log2FoldChange, 3)), '\nPadj ', as.character(round(padj,3))))) +
  facet_grid(GeneSymbol ~ .) +
  ggtitle('MEF2C and FOXJ3 DESeq2 normalized counts values in 11 drugs') +
  theme_bw() +
  ylim(c(0,350))
ggsave(DESeq2NC_FOXJ3MEF2C_moreDrugs_bees2, filename = paste0(graphDir, '/differentialExpression/allSamples/iCards/DESeq2NC_FOXJ3MEF2C_moreDrugs_bees2.pdf'), width = 12, height = 8, useDingbats = F)

set.seed(42639)
DESeq2NC_FOXJ3MEF2C_moreDrugs_bees2_forfig_wLabs <- ggplot() +
  geom_point(data = desnc_temp_bound2, aes(Product_col, norm_counts, color = isPerturbed, group = as.character(isPerturbed)), position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.35), size = 0.5) +
  geom_errorbar(data = desnc_temp_sums_bound2, aes(Product_col, ymin = meanNC, ymax = meanNC, color = as.character(isPerturbed)), width = 0.8, position = position_dodge2(width = 0.9), alpha = 0.8, stat = 'identity') +
  scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'orange')) +
  # geom_text(data = temp_DE_res, aes(Product_col, 320, label = paste0('LFC ', as.character(round(log2FoldChange, 3)), '\nPadj ', as.character(round(padj,3))))) +
  facet_grid(GeneSymbol ~ .) +
  # ggtitle('MEF2C and FOXJ3 DESeq2 normalized counts values in 3 drugs') +
  theme_bw() +
  theme(legend.position = 'none') +
  ylim(c(0,200))
ggsave(DESeq2NC_FOXJ3MEF2C_moreDrugs_bees2_forfig_wLabs, filename = paste0(graphDir, '/differentialExpression/allSamples/iCards/DESeq2NC_FOXJ3MEF2C_moreDrugs_bees2_forfig_wLabs.pdf'), width = 6, height = 4, useDingbats = F)

set.seed(42639)
DESeq2NC_FOXJ3MEF2C_moreDrugs_bees2_forfig_noLabs <- ggplot() +
  geom_point(data = desnc_temp_bound2, aes(Product_col, norm_counts, color = isPerturbed, group = as.character(isPerturbed)), position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.35), size = 0.1) +
  geom_errorbar(data = desnc_temp_sums_bound2, aes(Product_col, ymin = meanNC, ymax = meanNC, color = as.character(isPerturbed)), width = 0.8, position = position_dodge2(width = 0.9), alpha = 0.8, stat = 'identity') +
  scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'orange')) +
  # geom_text(data = temp_DE_res, aes(Product_col, 320, label = paste0('LFC ', as.character(round(log2FoldChange, 3)), '\nPadj ', as.character(round(padj,3))))) +
  facet_grid(GeneSymbol ~ ., ) +
  # ggtitle('MEF2C and FOXJ3 DESeq2 normalized counts values in 3 drugs') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  ylim(c(0,200))
ggsave(DESeq2NC_FOXJ3MEF2C_moreDrugs_bees2_forfig_noLabs, filename = paste0(graphDir, '/differentialExpression/allSamples/iCards/DESeq2NC_FOXJ3MEF2C_moreDrugs_bees2_forfig_noLabs.pdf'), width = 2.5, height = 1.2, useDingbats = F)

desnc_temp_bound_noDodge2 <- desnc_temp_bound2 %>%
  filter(isPerturbed == T) %>%
  dplyr::select(-condition) %>%
  bind_rows(desnc_temp_bound %>% filter(isPerturbed == F) %>% dplyr::select(-condition) %>% mutate(Product_col = 'control', Product.Name = 'DMSO only'))
desnc_temp_sums_bound_noDodge2 <- desnc_temp_sums_bound2 %>%
  filter(isPerturbed == T) %>%
  bind_rows(desnc_temp_sums_bound %>% filter(isPerturbed == F) %>% mutate(Product_col = 'control') %>% unique())

set.seed(42639)
DESeq2NC_FOXJ3MEF2C_moreDrugs_bees3_forfig <- ggplot() +
  geom_boxplot(data = desnc_temp_bound_noDodge2 %>% filter(isPerturbed == F), aes(Product_col, norm_counts), outlier.shape = NA) +
  geom_jitter(data = desnc_temp_bound_noDodge2, aes(Product_col, norm_counts, color = isPerturbed, alpha = isPerturbed), width = 0.1, height = 0) +
  geom_errorbar(data = desnc_temp_sums_bound_noDodge2 %>% filter(isPerturbed == T), aes(Product_col, ymin = meanNC, ymax = meanNC, color = as.character(isPerturbed)), width = 0.7, alpha = 0.8, stat = 'identity') +
  scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'chocolate2')) +
  scale_alpha_manual(values = c('FALSE' = 0.2, 'TRUE' = 1)) +
  facet_grid(GeneSymbol ~ .) +
  # ggtitle('MEF2C and FOXJ3 DESeq2 normalized counts values in 3 drugs') +
  theme_classic() +
  theme(legend.position = 'none') +
  ylim(c(0,200))
ggsave(DESeq2NC_FOXJ3MEF2C_moreDrugs_bees3_forfig, filename = paste0(graphDir, '/differentialExpression/allSamples/iCards/DESeq2NC_FOXJ3MEF2C_moreDrugs_bees3_forfig.pdf'), width = 3, height = 2, useDingbats = F)

set.seed(42639)
DESeq2NC_FOXJ3MEF2C_moreDrugs_bees3_forfig_orange_noLabs <- ggplot() +
  geom_boxplot(data = desnc_temp_bound_noDodge2 %>% filter(isPerturbed == F), aes(Product_col, norm_counts), outlier.shape = NA) +
  geom_jitter(data = desnc_temp_bound_noDodge2, aes(Product_col, norm_counts, color = isPerturbed, alpha = isPerturbed), width = 0.1, height = 0) +
  geom_errorbar(data = desnc_temp_sums_bound_noDodge2 %>% filter(isPerturbed == T), aes(Product_col, ymin = meanNC, ymax = meanNC, color = as.character(isPerturbed)), width = 0.7, alpha = 0.8, stat = 'identity') +
  scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'chocolate2')) +
  scale_alpha_manual(values = c('FALSE' = 0.2, 'TRUE' = 1)) +
  facet_grid(GeneSymbol ~ .) +
  # ggtitle('MEF2C and FOXJ3 DESeq2 normalized counts values in 3 drugs') +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  ylim(c(0,200))
ggsave(DESeq2NC_FOXJ3MEF2C_moreDrugs_bees3_forfig_orange_noLabs, filename = paste0(graphDir, '/differentialExpression/allSamples/iCards/DESeq2NC_FOXJ3MEF2C_moreDrugs_bees3_forfig_orange_noLabs.pdf'), width = 6, height = 3, useDingbats = F)

set.seed(42639)
DESeq2NC_FOXJ3MEF2C_moreDrugs_bees3_forfig_magenta_noLabs <- ggplot() +
  geom_boxplot(data = desnc_temp_bound_noDodge2 %>% filter(isPerturbed == F), aes(Product_col, norm_counts), outlier.shape = NA) +
  geom_jitter(data = desnc_temp_bound_noDodge2, aes(Product_col, norm_counts, color = isPerturbed, alpha = isPerturbed), width = 0.1, height = 0) +
  geom_errorbar(data = desnc_temp_sums_bound_noDodge2 %>% filter(isPerturbed == T), aes(Product_col, ymin = meanNC, ymax = meanNC, color = as.character(isPerturbed)), width = 0.7, alpha = 0.8, stat = 'identity') +
  scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'magenta')) +
  scale_alpha_manual(values = c('FALSE' = 0.2, 'TRUE' = 1)) +
  facet_grid(GeneSymbol ~ .) +
  # ggtitle('MEF2C and FOXJ3 DESeq2 normalized counts values in 3 drugs') +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  ylim(c(0,200))
ggsave(DESeq2NC_FOXJ3MEF2C_moreDrugs_bees3_forfig_magenta_noLabs, filename = paste0(graphDir, '/differentialExpression/allSamples/iCards/DESeq2NC_FOXJ3MEF2C_moreDrugs_bees3_forfig_magenta_noLabs.pdf'), width = 6, height = 3, useDingbats = F)

set.seed(42639)
DESeq2NC_FOXJ3MEF2C_moreDrugs_bees3_forfig_grey_noLabs <- ggplot() +
  geom_boxplot(data = desnc_temp_bound_noDodge2 %>% filter(isPerturbed == F), aes(Product_col, norm_counts), outlier.shape = NA) +
  geom_jitter(data = desnc_temp_bound_noDodge2, aes(Product_col, norm_counts, color = isPerturbed, alpha = isPerturbed), width = 0.1, height = 0) +
  geom_errorbar(data = desnc_temp_sums_bound_noDodge2 %>% filter(isPerturbed == T), aes(Product_col, ymin = meanNC, ymax = meanNC, color = as.character(isPerturbed)), width = 0.7, alpha = 0.8, stat = 'identity') +
  scale_color_manual(values = c('FALSE' = 'grey30', 'TRUE' = 'black')) +
  scale_alpha_manual(values = c('FALSE' = 0.2, 'TRUE' = 1)) +
  facet_grid(GeneSymbol ~ .) +
  # ggtitle('MEF2C and FOXJ3 DESeq2 normalized counts values in 3 drugs') +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  ylim(c(0,200))
ggsave(DESeq2NC_FOXJ3MEF2C_moreDrugs_bees3_forfig_grey_noLabs, filename = paste0(graphDir, '/differentialExpression/allSamples/iCards/DESeq2NC_FOXJ3MEF2C_moreDrugs_bees3_forfig_grey_noLabs.pdf'), width = 6, height = 3, useDingbats = F)


# #### GO analysis moved to separate script
# library(clusterProfiler)
# library(org.Hs.eg.db)
# 
# iCard_highUpreg_indivSamps <- indiv_perturbability %>%
#   filter(cellType == 'iCard') %>%
#   filter(nDESeq2conditionsUp >= 4, log(meanRPM) > 2.5) %>%
#   dplyr::select(cellType, gene_id, nDESeq2conditionsUp)
# 
# iCard_lowUpreg_indivSamps <- indiv_perturbability %>%
#   filter(cellType == 'iCard') %>%
#   filter(nDESeq2conditionsUp <= 1, log(meanRPM) > 2.5) %>%
#   dplyr::select(cellType, gene_id, nDESeq2conditionsUp)
# 
# iCard_all_indivSamps <- indiv_perturbability %>%
#   filter(cellType == 'iCard') %>%
#   filter(log(meanRPM) > 2.5) %>%
#   dplyr::select(cellType, gene_id, nDESeq2conditionsUp)
# 
# iCard_highUpreg_indivSamps_df = bitr(iCard_highUpreg_indivSamps$gene_id,
#                                              fromType = "ENSEMBL",
#                                              toType = c("SYMBOL", "ENTREZID"),
#                                              OrgDb = org.Hs.eg.db)
# 
# iCard_lowUpreg_indivSamps_df = bitr(iCard_lowUpreg_indivSamps$gene_id,
#                                      fromType = "ENSEMBL",
#                                      toType = c("SYMBOL", "ENTREZID"),
#                                      OrgDb = org.Hs.eg.db)
# 
# 
# iCard_all_indivSamps_df = bitr(iCard_all_indivSamps$gene_id,
#                                fromType = "ENSEMBL",
#                                toType = c("SYMBOL", "ENTREZID"),
#                                OrgDb = org.Hs.eg.db)
# 
# iCard_go_overrep_tfs <- enrichGO(gene = iCard_highUpreg_indivSamps_df$ENTREZID,
#                                  universe      = iCard_all_indivSamps_df$ENTREZID,
#                                  OrgDb         = org.Hs.eg.db,
#                                  # keytype       = "ENSEMBL",
#                                  ont           = "CC",
#                                  pAdjustMethod = "BH",
#                                  # pvalueCutoff  = 1,
#                                  # qvalueCutoff  = 1,
#                                  readable      = TRUE)
# head(summary(iCard_go_overrep_tfs))
# dotplot(iCard_go_overrep_tfs, showCategory=10, by = "p.adjust")
# 
# upreg_forplot <- summary(iCard_go_overrep_tfs)[rank(summary(iCard_go_overrep_tfs)$qvalue),]
# rownames(upreg_forplot) <- NULL
# upreg_forplot2 <- upreg_forplot %>%
#   dplyr::select(-geneID) %>%
#   unique()
# upreg_forplot2$plotOrder <- nrow(upreg_forplot2):1
# upreg_forplot2$Description <- factor(upreg_forplot2$Description)
# upreg_forplot2$Description <- factor(upreg_forplot2$Description, levels = upreg_forplot2$Description[upreg_forplot2$plotOrder])
# 
# upreg_simp <- simplify(iCard_go_overrep_tfs)
# 
# iCard_highUpreg_GO_CC_dots_top10 <- ggplot(upreg_forplot2[1:10,], aes(Description, -log10(qvalue))) +
#   geom_point() +
#   coord_flip() +
#   theme_bw() +
#   ggtitle('Over-represented GO-CC\nin frequently up-regulated\nexpressed genes')
# ggsave(iCard_highUpreg_GO_CC_dots_top10, file = '~/Dropbox (RajLab)/Shared_IanM/cellid_201807_onward/iCard_highUpreg_GO_CC_dots_top10.pdf')
# 
# 
# iCard_go_low_overrep_tfs <- enrichGO(gene = iCard_lowUpreg_indivSamps_df$ENTREZID,
#                                  universe      = iCard_all_indivSamps_df$ENTREZID,
#                                  OrgDb         = org.Hs.eg.db,
#                                  # keytype       = "ENSEMBL",
#                                  ont           = "CC",
#                                  pAdjustMethod = "BH",
#                                  # pvalueCutoff  = 1,
#                                  # qvalueCutoff  = 1,
#                                  readable      = TRUE)
# head(summary(iCard_go_low_overrep_tfs))
# 
# notupreg_forplot <- summary(iCard_go_low_overrep_tfs)
# rownames(notupreg_forplot) <- NULL
# notupreg_forplot2 <- notupreg_forplot %>%
#   dplyr::select(-geneID) %>%
#   unique()
# notupreg_forplot2$plotOrder <- nrow(notupreg_forplot2):1
# notupreg_forplot2$Description <- factor(notupreg_forplot2$Description)
# notupreg_forplot2$Description <- factor(notupreg_forplot2$Description, levels = notupreg_forplot2$Description[notupreg_forplot2$plotOrder])
# 
# iCard_notUpreg_GO_CC_dots_top10 <- ggplot(notupreg_forplot2, aes(Description, -log10(qvalue))) +
#   geom_point() +
#   coord_flip() +
#   theme_bw()+
#   ggtitle('Over-represented GO-CC\nin infrequently up-regulated\nexpressed genes')
# ggsave(iCard_notUpreg_GO_CC_dots_top10, file = '~/Dropbox (RajLab)/Shared_IanM/cellid_201807_onward/iCard_notUpreg_GO_CC_dots_top10.pdf')
# 
