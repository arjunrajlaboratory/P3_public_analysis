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

# Main script for analyzing all 6 plates' worth of RNAtag-seq data
## Basic clustering/PCA on TPM data
## Running DESeq2
### Comparing gene-level results per-condition
## Calculating varaibility stats
### Plotting markers
### Doing KEGG/GO on high/low pert sets

library(metaRNASeq)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(magrittr)

# if running manually in Rstudio, start here and use:
# projectDir = '~/Dropbox (RajLab)/Shared_IanM/cellid_201807_onward/'
# procDataSubdir = 'procDataScripted_test429'
# graphSubdir = 'graphs'

procDataDir = paste0(projectDir, procDataSubdir)
graphDir = paste0(projectDir, graphSubdir)

setwd(projectDir)
cat("Starting DE_metaanalysis.R...\n")
source(paste0(projectDir, "analysisScripts/RNAtagUtilities.R"))
cat(paste0("Working in ", getwd(), "\n"))

tfGeneIDfile = "annotations/TF_gene_ids.csv"
tfTab = as_tibble(read.csv(tfGeneIDfile, stringsAsFactors = F))
colnames(tfTab) = "gene_id"
# 
# newTFfile = "annotations/2018_HTFreview.xls"
# newTFtab = as_tibble(read_xls(newTFfile, col_names = T))

geneIDfile = "annotations/hg19gene_idToGeneSymbol.tsv"
geneIDtoGeneSymbol <- as_tibble(read.csv(geneIDfile, sep='\t', stringsAsFactors = F))

markerFile = "miscellaneous/markerGeneList.csv"
markerTab = as_tibble(read.csv(markerFile, header = T, stringsAsFactors = F)) %>%
  left_join(geneIDtoGeneSymbol, by = "GeneSymbol")
tempMarks <- geneIDtoGeneSymbol %>%
  filter(GeneSymbol %in% c('ESRRG', 'ZFPM2', 'MESP1', 'MYOCD')) %>%
  mutate(cellTypeMark = 'iCard')
markerTab2 <- markerTab %>%
  bind_rows(tempMarks)

factors_7F <- c('GATA4', 'MEF2C', 'TBX5', 'MESP1', 'ESRRG', 'MYOCD', 'ZFPM2')

barrierFile = "miscellaneous/fibro_barriers_GeneSymbols.txt"
barrierTab = as_tibble(read.table(barrierFile, header = T, stringsAsFactors = F)) %>%
  left_join(geneIDtoGeneSymbol, by = "GeneSymbol")

targs16 <- c('SKIL', 'YBX3', 'TSC22D1', 'CERS2', 'KLF13', 'TBX3', 'ID1', 'ATOH8', 'ZNF652', 'NFATC4', 'ZBTB38', 'LARP1', 'CEBPB', 'ID3', 'PRRX2', 'RUNX1')
targs8 <- c('CERS2', 'KLF13', 'ATOH8', 'ZNF652', 'ZBTB38', 'CEBPB', 'PRRX2', 'RUNX1')


# load DE and mean expression data
iCard_DE_results_all <- as_tibble(read.table(paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_DESeqResults.txt'), header = T, stringsAsFactors = F)) %>%
  mutate(directionEffect = ifelse(ifelse(padj >= 0.1, 'NONSIG', ifelse(log2FoldChange < 0, 'DOWN', 'UP'))))

indiv_perturbability = as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_variability_metrics.txt"), header = T, stringsAsFactors = F, sep = "\t"))
indiv_perturbability[is.na(indiv_perturbability)] = 0

iCard_meanRPMs <- indiv_perturbability %>%
  filter(cellType == 'iCard') %>%
  dplyr::select(gene_id, meanRPM)

# reformat
iCard_DE_wide_pval <- iCard_DE_results_all %>%
  dplyr::select(gene_id, pvalue, Product.Name) %>%
  dplyr::filter(Product.Name != 'DMSOHeLa-HeLa') %>%
  pivot_wider(names_from = Product.Name, values_from = pvalue)
iCard_DE_wide_pval <- as.data.frame(iCard_DE_wide_pval)
rownames(iCard_DE_wide_pval) <- iCard_DE_wide_pval$gene_id
iCard_DE_wide_pval %<>% dplyr::select(-gene_id)

# perform Fisher's method for combining p-values for each gene across all perturbation conditions tested, set FDR = 0.05
fishcomb_iCard <- fishercomb(iCard_DE_wide_pval, BHth = 0.05)

# log10 transform (+ 1e-16, which is 1 order of magnitude smaller than the smallest non-0 p-value)
pseudo = 1e-16
fishcomb_iCard_tbl <- tibble(
  gene_id = rownames(iCard_DE_wide_pval),
  fishcomb_pvalue = fishcomb_iCard$rawpval,
  fishcomb_BHadjpvalue = fishcomb_iCard$adjpval,
  cellType = 'iCard'
) %>%
  mutate(log10_fishcomb_BHadjpvalue = log10(fishcomb_BHadjpvalue + pseudo))

fishcomb_iCard_tbl_TF <- fishcomb_iCard_tbl %>%
  inner_join(geneIDtoGeneSymbol) %>%
  inner_join(tfTab) 

fishcomb_iCard_tbl_marks <- fishcomb_iCard_tbl %>% inner_join(markerTab2) %>% filter(cellTypeMark == 'iCard') %>% unique()
set.seed(432784) # add jitter to linesegs to enable visualization
fishcomb_iCard_tbl_marks$log10_fishcomb_BHadjpvalue_jit <- fishcomb_iCard_tbl_marks$log10_fishcomb_BHadjpvalue + runif(min = -0.25, max = 0.25, n = length(fishcomb_iCard_tbl_marks$fishcomb_BHadjpvalue))
fishcomb_iCard_tbl_marks %<>%
  dplyr::rename(log10_fishcomb_BHadjpvalue_orig = log10_fishcomb_BHadjpvalue,
                log10_fishcomb_BHadjpvalue = log10_fishcomb_BHadjpvalue_jit)

iCard_allTFs_metaanalysis_BHcorrected_FishersMethod_markers <- ggplot() + 
  geom_histogram(data = fishcomb_iCard_tbl_TF, aes(log10_fishcomb_BHadjpvalue)) +
  geom_linerange(data = fishcomb_iCard_tbl_marks, aes(x = log10_fishcomb_BHadjpvalue, ymin = 0, ymax = 100), color = 'red') +
  geom_text_repel(data = fishcomb_iCard_tbl_marks, aes(x = log10_fishcomb_BHadjpvalue, y = 100, label = GeneSymbol), color = 'red') +
  xlab("log10(combined p-value + 1e-16)\nFisher's method, Benjamini-Hochberg corrected") +
  ylab("Number of transcription factor genes") +
  theme_classic()

iCard_TFsMinRPM20_metaanalysis_BHcorrected_FishersMethod_markers <- ggplot() + 
  geom_histogram(data = fishcomb_iCard_tbl_TF %>% inner_join(iCard_meanRPMs) %>% filter(meanRPM > 20), aes(log10_fishcomb_BHadjpvalue)) +
  geom_linerange(data = fishcomb_iCard_tbl_marks %>% inner_join(iCard_meanRPMs) %>% filter(meanRPM > 20), aes(x = log10_fishcomb_BHadjpvalue, ymin = 0, ymax = 100), color = 'red') +
  geom_text_repel(data = fishcomb_iCard_tbl_marks %>% inner_join(iCard_meanRPMs) %>% filter(meanRPM > 20), aes(x = log10_fishcomb_BHadjpvalue, y = 100, label = GeneSymbol), color = 'red') +
  xlab("log10(combined p-value + 1e-16)\nFisher's method, Benjamini-Hochberg corrected") +
  ylab("Number of transcription factor genes") +
  ggtitle("Statistical meta-analysis of all iPSC-CM perturbations\nAll transcription factors > 20 RPM") +
  theme_classic()
ggsave(iCard_TFsMinRPM20_metaanalysis_BHcorrected_FishersMethod_markers, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_TFsMinRPM20_metaanalysis_BHcorrected_FishersMethod_markers_withLabs.pdf'))

iCard_TFsMinRPM20_metaanalysis_BHcorrected_FishersMethod_markers_noLabs <- ggplot() + 
  geom_histogram(data = fishcomb_iCard_tbl_TF %>% inner_join(iCard_meanRPMs) %>% filter(meanRPM > 20), aes(log10_fishcomb_BHadjpvalue)) +
  geom_linerange(data = fishcomb_iCard_tbl_marks %>% inner_join(iCard_meanRPMs) %>% filter(meanRPM > 20), aes(x = log10_fishcomb_BHadjpvalue, ymin = 0, ymax = 100), color = 'red') +
  geom_text_repel(data = fishcomb_iCard_tbl_marks %>% inner_join(iCard_meanRPMs) %>% filter(meanRPM > 20), aes(x = log10_fishcomb_BHadjpvalue, y = 100, label = GeneSymbol), color = 'red') +
  # xlab("log10(combined p-value + 1e-16)\nFisher's method, Benjamini-Hochberg corrected") +
  # ylab("Number of transcription factor genes") +
  # ggtitle("Statistical meta-analysis of all iPSC-CM perturbations\nAll transcription factors > 20 RPM") +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.title = element_blank())
ggsave(iCard_TFsMinRPM20_metaanalysis_BHcorrected_FishersMethod_markers_noLabs, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_TFsMinRPM20_metaanalysis_BHcorrected_FishersMethod_markers_noLabs.pdf'), width = 5, height = 4)
