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

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(readxl)
library(magrittr)
library(e1071)
library(matrixStats)

# if running manually in Rstudio, start here and use:
# projectDir = '~/Dropbox (RajLab)/Shared_IanM/cellid_201807_onward/'
# procDataSubdir = 'procDataScripted_test429'
# graphSubdir = 'graphs'

procDataDir = paste0(projectDir, procDataSubdir)
graphDir = paste0(projectDir, graphSubdir)

setwd(projectDir)
cat("Starting graphing_perturbability.R...\n")
source(paste0(projectDir, "analysisScripts/RNAtagUtilities.R"))
cat(paste0("Working in ", getwd(), "\n"))

tfGeneIDfile = "annotations/TF_gene_ids.csv"
tfTab = as_tibble(read.csv(tfGeneIDfile, stringsAsFactors = F))
colnames(tfTab) = "gene_id"

newTFfile = "annotations/2018_HTFreview.xls"
newTFtab = as_tibble(read_xls(newTFfile, col_names = T))

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


zhouOlson_activators = as_tibble(read.table("miscellaneous/ZhouOlson2017_TableS1_activators_reorg.txt", sep = "\t", header = T, stringsAsFactors = F))
zhouOlson_inhibitors = as_tibble(read.table("miscellaneous/ZhouOlson2017_TableS1_inhibitors_reorg.txt", sep = "\t", header = T, stringsAsFactors = F))
zhouOlson_allORFs = as_tibble(read.table("miscellaneous/ZhouOlson2017_TableS1_allORFs_reorg.txt", sep = "\t", header = T, stringsAsFactors = F))

## Load perturbability data
# for normalized values, working with default windowRadius = 20
indiv_perturbability = as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_variability_metrics.txt"), header = T, stringsAsFactors = F, sep = "\t"))
indiv_perturbability_RPMnorm20 = as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_tpmFilt_variability_RPMnormMetrics_window20.txt"), header = T, stringsAsFactors = F, sep = "\t"))
indiv_perturbability_tfOnlyRPMnorm20 = as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_tfOnly_variability_RPMnormMetrics_window20.txt"), header = T, stringsAsFactors = F, sep = "\t"))
# indiv_perturbability_newTfOnlyRPMnorm20 = as_tibble(read.table(paste0(procDataDir, "/allExperiments/newTFdef/all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_window20.txt"), header = T, stringsAsFactors = F, sep = "\t")) %>% left_join(geneIDtoGeneSymbol, by = "gene_id")

indiv_perturbability[is.na(indiv_perturbability)] = 0
indiv_perturbability_RPMnorm20[is.na(indiv_perturbability_RPMnorm20)] = 0
indiv_perturbability_tfOnlyRPMnorm20[is.na(indiv_perturbability_tfOnlyRPMnorm20)] = 0

# for grouped values, working with default seed = 4930
# grouped_perturbability = as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_groupedMeanRPM_readFilt_manualFilt_variability_metrics_seed4930.txt"), header = T, stringsAsFactors = F, sep = "\t"))
# grouped_perturbability_RPMnorm20 = as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_groupedMeanRPM_readFilt_manualFilt_tpmFilt_variability_RPMnormMetrics_window20_seed4930.txt"), header = T, stringsAsFactors = F, sep = "\t"))
# grouped_perturbability_tfOnlyRPMnorm20 = as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_groupedMeanRPM_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_window20_seed4930.txt"), header = T, stringsAsFactors = F, sep = "\t"))

iCard_DE_results_sig <- as_tibble(read.table(paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_DESeqResults.txt'), header = T, stringsAsFactors = F)) %>%
  mutate(directionEffect = ifelse(ifelse(padj >= 0.1, 'NONSIG', ifelse(log2FoldChange < 0, 'DOWN', 'UP'))))

iCard_DE_results_all <- as_tibble(read.table(paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_DESeqResults.txt'), header = T, stringsAsFactors = F))

# ggplot() + 
#   geom_density_2d(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(meanTPM > 10, cellType == "iCard"),
#              aes(log(meanTPM), deltaKurtosis), alpha = 0.2) +
#   geom_point(data = indiv_perturbability %>% 
#                inner_join(tfTab, by = "gene_id") %>%
#                inner_join(geneIDtoGeneSymbol,by = "gene_id") %>% 
#                inner_join(markerTab, by = "GeneSymbol") %>% 
#                filter(meanTPM > 10, cellType == "iCard", cellTypeMark == "iCard"),
#              aes(log(meanTPM), deltaKurtosis), color = "red") +
#   geom_text(data = indiv_perturbability %>% 
#                     inner_join(tfTab, by = "gene_id") %>%
#                     inner_join(geneIDtoGeneSymbol,by = "gene_id") %>% 
#                     inner_join(markerTab, by = "GeneSymbol") %>% 
#                     filter(meanTPM > 10, cellType == "iCard", cellTypeMark == "iCard") %>% unique(),
#                   aes(log(meanTPM), deltaKurtosis, label = GeneSymbol), color = "red", vjust = 0, hjust = 0)


# ggplot(indiv_perturbability %>% filter(meanTPM > 10), aes(log(meanRPM), controlKurtosis)) +
#   geom_point() +
#   facet_grid(~cellType) +
#   ylim(c(-5,30))
# ggplot(indiv_perturbability %>% filter(meanTPM > 10), aes(log(meanRPM), perturbedKurtosis)) +
#   geom_point() +
#   facet_grid(~cellType)+
#   ylim(c(0,30))
# 

## perturbability measures vs control expression

# ggplot()+
#   geom_density2d(data = indiv_perturbability_tfOnlyRPMnorm20 %>% filter(cellType == "iCard"), aes(log(meanTPM), normDeltaSkewness)) +
#   geom_point(data = indiv_perturbability_tfOnlyRPMnorm20 %>% filter(cellType == "iCard") %>%  inner_join(markerTab, by = "gene_id"), 
#              aes(log(meanTPM), normDeltaSkewness), color = "red") +
#   geom_text_repel(data = indiv_perturbability_tfOnlyRPMnorm20 %>% filter(cellType == "iCard") %>% inner_join(markerTab, by = "gene_id"), 
#                   aes(log(meanTPM), normDeltaSkewness, label = GeneSymbol), color = "red")
# 
# ggplot()+
#   geom_density2d(data = indiv_perturbability_newTfOnlyRPMnorm20 %>% filter(cellType == "iCard"), aes(log(meanTPM), normDeltaSkewness)) +
#   geom_point(data = indiv_perturbability_newTfOnlyRPMnorm20 %>% filter(cellType == "iCard") %>%  inner_join(markerTab, by = "gene_id"), 
#              aes(log(meanTPM), normDeltaSkewness), color = "red") +
#   geom_text_repel(data = indiv_perturbability_newTfOnlyRPMnorm20 %>% filter(cellType == "iCard") %>% inner_join(markerTab, by = "gene_id"), 
#                   aes(log(meanTPM), normDeltaSkewness, label = GeneSymbol), color = "red")

if(!dir.exists(paste0(graphDir, "/perturbability/"))){
  dir.create(paste0(graphDir, "/perturbability/"))
}
iCard_markers_tfOnly_normDeltaKurtosis <- ggplot()+
  geom_density2d(data = indiv_perturbability_tfOnlyRPMnorm20 %>% filter(cellType == "iCard"), aes(log(meanTPM), normDeltaKurtosis)) +
  geom_point(data = indiv_perturbability_tfOnlyRPMnorm20 %>% filter(cellType == "iCard") %>%  inner_join(markerTab, by = "gene_id"), 
             aes(log(meanTPM), normDeltaKurtosis), color = "red") +
  geom_text_repel(data = indiv_perturbability_tfOnlyRPMnorm20 %>% filter(cellType == "iCard") %>% inner_join(markerTab, by = "gene_id") %>% unique(), 
                  aes(log(meanTPM), normDeltaKurtosis, label = GeneSymbol), color = "red", size = 8) +
  theme_bw() +
  ggtitle("TF genes only perturbability in iCards")
ggsave(iCard_markers_tfOnly_normDeltaKurtosis, file = paste0(graphDir, "/perturbability/iCard_markers_tfOnly_normDeltaKurtosis.pdf"), width = 8, height = 8)

iCard_markers_tfOnly_deltaKurtosis <- ggplot()+
  geom_density2d(data = indiv_perturbability_tfOnlyRPMnorm20 %>% filter(cellType == "iCard"), aes(log(meanTPM), deltaKurtosis)) +
  geom_point(data = indiv_perturbability_tfOnlyRPMnorm20 %>% filter(cellType == "iCard") %>%  inner_join(markerTab, by = "gene_id"), 
             aes(log(meanTPM), deltaKurtosis), color = "red") +
  geom_text(data = indiv_perturbability_tfOnlyRPMnorm20 %>% filter(cellType == "iCard") %>% inner_join(markerTab, by = "gene_id") %>% unique(), 
                  aes(log(meanTPM), deltaKurtosis, label = GeneSymbol), color = "red", size = 5, vjust = 0, hjust = 0) +
  theme_bw() +
  ggtitle("TF genes only perturbability in iCards")
ggsave(iCard_markers_tfOnly_deltaKurtosis, file = paste0(graphDir, "/perturbability/iCard_markers_tfOnly_deltaKurtosis.pdf"), width = 8, height = 8, useDingbats = F)


iCard_markers_tpmFilt_normDeltaKurtosis <- ggplot()+
  geom_density2d(data = indiv_perturbability_RPMnorm20 %>% filter(cellType == "iCard"), aes(log(meanTPM), normDeltaKurtosis)) +
  geom_point(data = indiv_perturbability_RPMnorm20 %>% filter(cellType == "iCard") %>%  inner_join(markerTab, by = "gene_id") %>% unique(), 
             aes(log(meanTPM), normDeltaKurtosis), color = "red") +
  geom_text_repel(data = indiv_perturbability_RPMnorm20 %>% filter(cellType == "iCard") %>% inner_join(markerTab, by = "gene_id") %>% unique(), 
                  aes(log(meanTPM), normDeltaKurtosis, label = GeneSymbol), color = "red") +
  theme_bw() +
  ggtitle("All genes passing TPM > 10 perturbability in iCards")
ggsave(iCard_markers_tpmFilt_normDeltaKurtosis, file = paste0(graphDir, "/perturbability/iCard_markers_tpmFilt_normDeltaKurtosis.pdf"), width = 8, height = 8)


iCard_markers_tfOnly_high_NormDeltaKurtosis_vs_upregfreq <- ggplot() +
  geom_point(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
               filter(cellType == "iCard", log(meanTPM) > 4, normDeltaKurtosis > 0.7) %>% inner_join(geneIDtoGeneSymbol, by = "gene_id"),
             aes(normDeltaKurtosis, nDESeq2conditionsUp/nDESeq2conditionsAll)) +
  geom_text_repel(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
               filter(cellType == "iCard", log(meanTPM) > 4, normDeltaKurtosis > 0.7) %>% anti_join(markerTab, by = "gene_id") %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>% unique(),
             aes(normDeltaKurtosis, nDESeq2conditionsUp/nDESeq2conditionsAll, label = GeneSymbol)) +
  geom_point(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
               filter(cellType == "iCard", log(meanTPM) > 4, normDeltaKurtosis > 0.7) %>% inner_join(markerTab, by = "gene_id"),
             aes(normDeltaKurtosis, nDESeq2conditionsUp/nDESeq2conditionsAll), color = "red") +
  geom_text_repel(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
               filter(cellType == "iCard", log(meanTPM) > 4, normDeltaKurtosis > 0.7) %>% inner_join(markerTab, by = "gene_id") %>% unique(),
             aes(normDeltaKurtosis, nDESeq2conditionsUp/nDESeq2conditionsAll, label = GeneSymbol), color = "red") +
  theme_classic() + 
  ggtitle("TFs highly perturbable in iCards\nwith up-regulatory frequency\nln(meanTPM) > 4 only") +
  ylab("Fraction of DiffExp conditions in which UP-regulated") +
  xlab("Perturbability (normalized change in kurtosis)")
# plot(iCard_markers_tfOnly_normDeltaKurtosis_vs_upregfreq)
ggsave(iCard_markers_tfOnly_high_NormDeltaKurtosis_vs_upregfreq, file = paste0(graphDir, "/perturbability/iCard_markers_tfOnly_high_NormDeltaKurtosis_vs_upregfreq.pdf"), width = 10, height = 10, useDingbats = F)  

iCard_noMarkers_tfOnly_high_NormDeltaKurtosis_vs_upregfreq <- ggplot() +
  geom_point(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
               filter(cellType == "iCard", log(meanTPM) > 4, normDeltaKurtosis > 0.7) %>% inner_join(geneIDtoGeneSymbol, by = "gene_id"),
             aes(normDeltaKurtosis, nDESeq2conditionsUp/nDESeq2conditionsAll)) +
  geom_text_repel(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
                    filter(cellType == "iCard", log(meanTPM) > 4, normDeltaKurtosis > 0.7) %>% anti_join(markerTab, by = "gene_id") %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>% unique(),
                  aes(normDeltaKurtosis, nDESeq2conditionsUp/nDESeq2conditionsAll, label = GeneSymbol)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle("Gene group codename: black whale") +
  ylab("Fraction of DiffExp conditions in which UP-regulated") +
  xlab("Perturbability (normalized change in kurtosis)")
# plot(iCard_markers_tfOnly_normDeltaKurtosis_vs_upregfreq)
ggsave(iCard_noMarkers_tfOnly_high_NormDeltaKurtosis_vs_upregfreq, file = paste0(graphDir, "/perturbability/iCard_noMarkers_tfOnly_BLACKWHALE_NormDeltaKurtosis_vs_upregfreq.pdf"), width = 10, height = 10, useDingbats = F)  

iCard_noMarkers_tfOnly_mid_NormDeltaKurtosis_vs_upregfreq <- ggplot() +
  geom_point(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
               filter(cellType == "iCard", log(meanTPM) > 4, normDeltaKurtosis <= 0.7, normDeltaKurtosis > 0.3) %>% inner_join(geneIDtoGeneSymbol, by = "gene_id"),
             aes(normDeltaKurtosis, nDESeq2conditionsUp/nDESeq2conditionsAll)) +
  geom_text_repel(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
                    filter(cellType == "iCard", log(meanTPM) > 4, normDeltaKurtosis <= 0.7, normDeltaKurtosis > 0.3) %>% anti_join(markerTab, by = "gene_id") %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>% unique(),
                  aes(normDeltaKurtosis, nDESeq2conditionsUp/nDESeq2conditionsAll, label = GeneSymbol)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle("Gene group codename: scarlet cow") +
  ylab("Fraction of DiffExp conditions in which UP-regulated") +
  xlab("Perturbability (normalized change in kurtosis)")
# plot(iCard_markers_tfOnly_normDeltaKurtosis_vs_upregfreq)
ggsave(iCard_noMarkers_tfOnly_mid_NormDeltaKurtosis_vs_upregfreq, file = paste0(graphDir, "/perturbability/iCard_noMarkers_tfOnly_SCARLETCOW_NormDeltaKurtosis_vs_upregfreq.pdf"), width = 10, height = 10, useDingbats = F)  

iCard_noMarkers_tfOnly_low_NormDeltaKurtosis_vs_upregfreq <- ggplot() +
  geom_point(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
               filter(cellType == "iCard", log(meanTPM) > 4, normDeltaKurtosis < 0.3) %>% inner_join(geneIDtoGeneSymbol, by = "gene_id"),
             aes(normDeltaKurtosis, nDESeq2conditionsUp/nDESeq2conditionsAll)) +
  geom_text_repel(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
                    filter(cellType == "iCard", log(meanTPM) > 4, normDeltaKurtosis < 0.3) %>% anti_join(markerTab, by = "gene_id") %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>% unique(),
                  aes(normDeltaKurtosis, nDESeq2conditionsUp/nDESeq2conditionsAll, label = GeneSymbol)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle("Gene group codename: blue falcon") +
  ylab("Fraction of DiffExp conditions in which UP-regulated") +
  xlab("Perturbability (normalized change in kurtosis)")
# plot(iCard_markers_tfOnly_normDeltaKurtosis_vs_upregfreq)
ggsave(iCard_noMarkers_tfOnly_low_NormDeltaKurtosis_vs_upregfreq, file = paste0(graphDir, "/perturbability/iCard_noMarkers_tfOnly_BLUEFALCON_NormDeltaKurtosis_vs_upregfreq.pdf"), width = 10, height = 10, useDingbats = F)  




#norm_DESeq2conditionsUp
set.seed(3561)
iCard_norm_nDESeqConditionsUpvslogmeanRPM <- ggplot() +
  geom_jitter(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
                filter(cellType == "iCard"),
              aes(log(meanRPM), norm_nDESeq2conditionsUp), alpha = 0.3, height = 0.01) +
  geom_point(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
               filter(cellType == "iCard") %>% inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == 'iCard'),
             aes(log(meanRPM), norm_nDESeq2conditionsUp), color = 'red') +
  geom_text_repel(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
               filter(cellType == "iCard") %>% inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == 'iCard') %>% unique(),
             aes(log(meanRPM), norm_nDESeq2conditionsUp, label = GeneSymbol), color = 'red') +
  theme_bw() +
  xlim(c(0,7.6)) +
  ggtitle('Perturbability normalized to average expression level\nnDESeq2conditionsUp in iCards, markers in red') +
  xlab('Log(mean RPM in controls)') +
  ylab('Normalized # conditions UP')
ggsave(iCard_norm_nDESeqConditionsUpvslogmeanRPM, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_tfOnly_normUPvsRPM.pdf'), width = 7, height = 7)

set.seed(3561)
fibro_norm_nDESeqConditionsUpvslogmeanRPM_barriersOnly <- ggplot() +
  geom_jitter(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
                filter(cellType == "GM00942"),
              aes(log(meanRPM), norm_nDESeq2conditionsUp), alpha = 0.3, height = 0.01) +
  geom_point(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
               filter(cellType == "GM00942") %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>% filter(GeneSymbol %in% targs8),
             aes(log(meanRPM), norm_nDESeq2conditionsUp), color = 'darkolivegreen4') +
  geom_text_repel(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
                    filter(cellType == "GM00942") %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>% filter(GeneSymbol %in% targs8) %>% unique(),
                  aes(log(meanRPM), norm_nDESeq2conditionsUp, label = GeneSymbol), color = 'darkolivegreen4') +
  theme_bw() +
  xlim(c(0,7.6)) +
  ggtitle('Perturbability normalized to average expression level\nnDESeq2conditionsUp in fibroblasts, barriers in green') +
  xlab('Log(mean RPM in controls)') +
  ylab('Normalized # conditions UP')
ggsave(fibro_norm_nDESeqConditionsUpvslogmeanRPM_barriersOnly, file = paste0(graphDir, '/differentialExpression/allSamples/fibroblasts/fibro_tfOnly_normUPvsRPM_barriersOnly.pdf'), width = 7, height = 7)

set.seed(3561)
fibro_norm_nDESeqConditionsUpvslogmeanRPM_barriersAndNonbarriers <- ggplot() +
  geom_jitter(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
                filter(cellType == "GM00942"),
              aes(log(meanRPM), norm_nDESeq2conditionsUp), alpha = 0.3, height = 0.01) +
  geom_point(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
               filter(cellType == "GM00942") %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>% filter(GeneSymbol %in% targs8),
             aes(log(meanRPM), norm_nDESeq2conditionsUp), color = 'darkolivegreen4') +
  geom_text_repel(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
                    filter(cellType == "GM00942") %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>% filter(GeneSymbol %in% targs8) %>% unique(),
                  aes(log(meanRPM), norm_nDESeq2conditionsUp, label = GeneSymbol), color = 'darkolivegreen4') +
  geom_point(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
               filter(cellType == "GM00942") %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>% filter(GeneSymbol %in% targs16, !(GeneSymbol %in% targs8)),
             aes(log(meanRPM), norm_nDESeq2conditionsUp), color = 'darkolivegreen1') +
  geom_text_repel(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
                    filter(cellType == "GM00942") %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>% filter(GeneSymbol %in% targs16, !(GeneSymbol %in% targs8)) %>% unique(),
                  aes(log(meanRPM), norm_nDESeq2conditionsUp, label = GeneSymbol), color = 'darkolivegreen1') +
  theme_bw() +
  xlim(c(0,7.6)) +
  ggtitle('Perturbability normalized to average expression level\nnDESeq2conditionsUp in fibroblasts\nbarriers in dark green, no-effect in light green') +
  xlab('Log(mean RPM in controls)') +
  ylab('Normalized # conditions UP')
ggsave(fibro_norm_nDESeqConditionsUpvslogmeanRPM_barriersAndNonbarriers, file = paste0(graphDir, '/differentialExpression/allSamples/fibroblasts/fibro_tfOnly_normUPvsRPM_barriersAndNonbarriers.pdf'), width = 7, height = 7)


# fibroblast
fibro_barriers_tfOnly_normDeltaKurtosis <- ggplot()+
  geom_density2d(data = indiv_perturbability_tfOnlyRPMnorm20 %>% filter(cellType == "GM00942"), aes(log(meanTPM), normDeltaKurtosis)) +
  geom_point(data = indiv_perturbability_tfOnlyRPMnorm20 %>% filter(cellType == "GM00942") %>%  inner_join(barrierTab, by = "gene_id"), 
             aes(log(meanTPM), normDeltaKurtosis), color = "orange") +
  geom_text_repel(data = indiv_perturbability_tfOnlyRPMnorm20 %>% filter(cellType == "GM00942") %>% inner_join(barrierTab, by = "gene_id") %>% unique(), 
                  aes(log(meanTPM), normDeltaKurtosis, label = GeneSymbol), color = "orange", size = 8) +
  theme_bw() +
  ggtitle("TF genes only perturbability in fibroblasts")
ggsave(fibro_barriers_tfOnly_normDeltaKurtosis, file = paste0(graphDir, "/perturbability/fibro_barriers_tfOnly_normDeltaKurtosis.pdf"), width = 8, height = 8)


fibro_noMarkers_tfOnly_high_NormDeltaKurtosis_vs_upregfreq <- ggplot() +
  geom_point(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
               filter(cellType == "GM00942", log(meanTPM) > 4, normDeltaKurtosis > 0.7) %>% inner_join(geneIDtoGeneSymbol, by = "gene_id"),
             aes(normDeltaKurtosis, nDESeq2conditionsUp/nDESeq2conditionsAll)) +
  geom_text_repel(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
                    filter(cellType == "GM00942", log(meanTPM) > 4, normDeltaKurtosis > 0.7) %>% anti_join(barrierTab, by = "gene_id") %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>% unique(),
                  aes(normDeltaKurtosis, nDESeq2conditionsUp/nDESeq2conditionsAll, label = GeneSymbol)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle("Gene group codename: purple porcupine") +
  ylab("Fraction of DiffExp conditions in which UP-regulated") +
  xlab("Perturbability (normalized change in kurtosis)")
# plot(iCard_markers_tfOnly_normDeltaKurtosis_vs_upregfreq)
ggsave(fibro_noMarkers_tfOnly_high_NormDeltaKurtosis_vs_upregfreq, file = paste0(graphDir, "/perturbability/fibro_noMarkers_tfOnly_PURPLEPORCUPINE_NormDeltaKurtosis_vs_upregfreq.pdf"), width = 10, height = 10, useDingbats = F)  

fibro_noMarkers_tfOnly_mid_NormDeltaKurtosis_vs_upregfreq <- ggplot() +
  geom_point(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
               filter(cellType == "GM00942", log(meanTPM) > 4, normDeltaKurtosis <= 0.7, normDeltaKurtosis > 0.3) %>% inner_join(geneIDtoGeneSymbol, by = "gene_id"),
             aes(normDeltaKurtosis, nDESeq2conditionsUp/nDESeq2conditionsAll)) +
  geom_text_repel(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
                    filter(cellType == "GM00942", log(meanTPM) > 4, normDeltaKurtosis <= 0.7, normDeltaKurtosis > 0.3) %>% anti_join(barrierTab, by = "gene_id") %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>% unique(),
                  aes(normDeltaKurtosis, nDESeq2conditionsUp/nDESeq2conditionsAll, label = GeneSymbol)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle("Gene group codename: green gerbil") +
  ylab("Fraction of DiffExp conditions in which UP-regulated") +
  xlab("Perturbability (normalized change in kurtosis)")
# plot(iCard_markers_tfOnly_normDeltaKurtosis_vs_upregfreq)
ggsave(fibro_noMarkers_tfOnly_mid_NormDeltaKurtosis_vs_upregfreq, file = paste0(graphDir, "/perturbability/fibro_noMarkers_tfOnly_GREENGERBIL_NormDeltaKurtosis_vs_upregfreq.pdf"), width = 10, height = 10, useDingbats = F)  

fibro_noMarkers_tfOnly_low_NormDeltaKurtosis_vs_upregfreq <- ggplot() +
  geom_point(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
               filter(cellType == "GM00942", log(meanTPM) > 4, normDeltaKurtosis < 0.3) %>% inner_join(geneIDtoGeneSymbol, by = "gene_id"),
             aes(normDeltaKurtosis, nDESeq2conditionsUp/nDESeq2conditionsAll)) +
  geom_text_repel(data = indiv_perturbability_tfOnlyRPMnorm20 %>% 
                    filter(cellType == "GM00942", log(meanTPM) > 4, normDeltaKurtosis < 0.3) %>% anti_join(barrierTab, by = "gene_id") %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>% unique(),
                  aes(normDeltaKurtosis, nDESeq2conditionsUp/nDESeq2conditionsAll, label = GeneSymbol)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle("Gene group codename: snow leopard") +
  ylab("Fraction of DiffExp conditions in which UP-regulated") +
  xlab("Perturbability (normalized change in kurtosis)")
# plot(iCard_markers_tfOnly_normDeltaKurtosis_vs_upregfreq)
ggsave(fibro_noMarkers_tfOnly_low_NormDeltaKurtosis_vs_upregfreq, file = paste0(graphDir, "/perturbability/fibro_noMarkers_tfOnly_SNOWLEOPARD_NormDeltaKurtosis_vs_upregfreq.pdf"), width = 10, height = 10, useDingbats = F)  



## DESeq2 up and down ####
# replace NAs with 0s, since left_join of DESeq results inappropriately made NAs

# indiv_perturbability[is.na(indiv_perturbability)] = 0

iCard_mean_hist2_horiz <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20), aes(nDESeq2conditionsUp),
                 binwidth = 1) +
  geom_vline(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(xintercept = meanTPM), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(x = meanTPM, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() + 
  xlim(c(-0.5, 12.5)) +
  xlab("mean TPM in controls")
ggsave(iCard_mean_hist2_horiz, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_mean_histogram_withMarkers.pdf'), width = 8, height = 4)


set.seed(3561)
iCard_markers_tfOnly_nUpDown <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard", meanTPM > 10) %>% anti_join(markerTab, by = "gene_id"), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), color = "blue", alpha = 0.3) +
  geom_point(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanTPM > 10) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanTPM > 10) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp, label = GeneSymbol), color = "red", size = 7) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  ylim(c(-0.5, 12.5)) +
  xlab("Number of conditions DOWN") + ylab("Number of conditions UP")
plot(iCard_markers_tfOnly_nUpDown)
ggsave(iCard_markers_tfOnly_nUpDown, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_markers_tfOnly_nUpDown.pdf'), width = 7, height = 7, useDingbats = F)

iCard_up_hist <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanTPM > 10), aes(nDESeq2conditionsUp),
                 binwidth = 1) +
  geom_vline(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanTPM > 10) %>% unique(), aes(xintercept = nDESeq2conditionsUp), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanTPM > 10) %>% unique(), aes(x = nDESeq2conditionsUp, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() + 
  xlim(c(-0.5, 12.5)) +
  xlab("Number of conditions UP") +
  coord_flip()

iCard_down_hist <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanTPM > 10), aes(nDESeq2conditionsDown),
                 binwidth = 1) +
  geom_vline() +
  geom_vline(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanTPM > 10) %>% unique(), aes(xintercept = nDESeq2conditionsDown), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanTPM > 10) %>% unique(), aes(x = nDESeq2conditionsDown, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() +
  xlab("Number of conditions DOWN")

empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank())

pdf(paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_markers_tfOnly_nUpDown_withMarginals.pdf'), width = 15, height = 10, useDingbats = F)
grid.arrange(iCard_down_hist, empty, iCard_markers_tfOnly_nUpDown, iCard_up_hist, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()

# grey points

set.seed(3561)
iCard_grey_markers_tfOnly_nUpDown <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard", meanTPM > 10) %>% anti_join(markerTab, by = "gene_id"), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), alpha = 0.3) +
  geom_point(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanTPM > 10) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanTPM > 10) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp, label = GeneSymbol), color = "red", size = 6) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  ylim(c(-0.5, 12.5)) +
  xlab("Number of conditions DOWN") + ylab("Number of conditions UP")
plot(iCard_grey_markers_tfOnly_nUpDown)
ggsave(iCard_grey_markers_tfOnly_nUpDown, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_grey_markers_tfOnly_nUpDown.pdf'), width = 7, height = 7, useDingbats = F)

pdf(paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_grey_markers_tfOnly_nUpDown_withMarginals.pdf'), width = 9, height = 6, useDingbats = F)
grid.arrange(iCard_down_hist, empty, iCard_grey_markers_tfOnly_nUpDown, iCard_up_hist, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()


iCard_up_hist2 <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanTPM > 10), aes(nDESeq2conditionsUp),
                 binwidth = 1) +
  geom_vline(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanTPM > 10) %>% unique(), aes(xintercept = nDESeq2conditionsUp), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanTPM > 10) %>% unique(), aes(x = nDESeq2conditionsUp, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() + 
  xlim(c(-0.5, 12.5)) +
  xlab("Number of conditions UP") +
  coord_flip()

iCard_up_hist2_horiz <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanTPM > 10), aes(nDESeq2conditionsUp),
                 binwidth = 1) +
  geom_vline(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanTPM > 10) %>% unique(), aes(xintercept = nDESeq2conditionsUp), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanTPM > 10) %>% unique(), aes(x = nDESeq2conditionsUp, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() + 
  xlim(c(-0.5, 12.5)) +
  xlab("Number of conditions UP")
ggsave(iCard_up_hist2_horiz, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_nDEup_histogram_withMarkers.pdf'), width = 8, height = 4)


iCard_down_hist2 <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanTPM > 10), aes(nDESeq2conditionsDown),
                 binwidth = 1) +
  geom_vline() +
  geom_vline(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanTPM > 10) %>% unique(), aes(xintercept = nDESeq2conditionsDown), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanTPM > 10) %>% unique(), aes(x = nDESeq2conditionsDown, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() +
  xlab("Number of conditions DOWN")

set.seed(3561)
iCard_grey_markers2_tfOnly_nUpDown <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard", meanTPM > 10) %>% anti_join(markerTab2, by = "gene_id"), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), alpha = 0.3) +
  geom_point(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanTPM > 10) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanTPM > 10) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp, label = GeneSymbol), color = "red", size = 6) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  ylim(c(-0.5, 12.5)) +
  xlab("Number of conditions DOWN") + ylab("Number of conditions UP")
plot(iCard_grey_markers2_tfOnly_nUpDown)
pdf(paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_grey_markers2_tfOnly_nUpDown_withMarginals.pdf'), width = 9, height = 6, useDingbats = F)
grid.arrange(iCard_down_hist2, empty, iCard_grey_markers2_tfOnly_nUpDown, iCard_up_hist2, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()

# grey points with 3P
set.seed(3561)
iCard_grey_markers_tfOnly_nUpDown_w3P <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard", meanTPM > 10) %>% anti_join(markerTab, by = "gene_id"), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), alpha = 0.3) +
  geom_point(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanTPM > 10) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
                    filter(cellType == "iCard", GeneSymbol %in% c('SP3', 'ZBTB10', 'ZBTB44'), meanTPM > 10) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp, label = GeneSymbol), color = "blue", size = 6) +
  geom_point(data = indiv_perturbability %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
               filter(cellType == "iCard", GeneSymbol %in% c('SP3', 'ZBTB10', 'ZBTB44'), meanTPM > 10) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), color = "blue") +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  ylim(c(-0.5, 12.5)) +
  xlab("Number of conditions DOWN") + ylab("Number of conditions UP")
# plot(iCard_grey_markers_tfOnly_nUpDown)
ggsave(iCard_grey_markers_tfOnly_nUpDown, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_grey_markers_tfOnly_nUpDown.pdf'), width = 7, height = 7, useDingbats = F)

pdf(paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_grey_markers_tfOnly_nUpDown_withMarginals_w3P.pdf'), width = 9, height = 6, useDingbats = F)
grid.arrange(iCard_down_hist, empty, iCard_grey_markers_tfOnly_nUpDown_w3P, iCard_up_hist, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()


# grey points with 3P and markers
set.seed(3561)
iCard_grey_markers_tfOnly_nUpDown_w3P_wMarkers <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard", meanTPM > 10) %>% anti_join(markerTab, by = "gene_id"), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), alpha = 0.3) +
  geom_point(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanTPM > 10) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), size = 2, color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanTPM > 10) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp, label = GeneSymbol), color = "red", size = 6) +
  geom_text_repel(data = indiv_perturbability %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
                    filter(cellType == "iCard", GeneSymbol %in% c('SP3', 'ZBTB10', 'ZBTB44'), meanTPM > 10) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp, label = GeneSymbol), color = "blue", size = 6) +
  geom_point(data = indiv_perturbability %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
               filter(cellType == "iCard", GeneSymbol %in% c('SP3', 'ZBTB10', 'ZBTB44'), meanTPM > 10) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), size = 2, color = "blue") +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  ylim(c(-0.5, 12.5)) +
  xlab("Number of conditions DOWN") + ylab("Number of conditions UP")
# plot(iCard_grey_markers_tfOnly_nUpDown_w3P_wMarkers)
ggsave(iCard_grey_markers_tfOnly_nUpDown_w3P_wMarkers, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_grey_markers_tfOnly_nUpDown_w3P_wMarkers.pdf'), width = 7, height = 7, useDingbats = F)

pdf(paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_grey_markers_tfOnly_nUpDown_withMarginals_w3P_wMarkers.pdf'), width = 9, height = 6, useDingbats = F)
grid.arrange(iCard_down_hist, empty, iCard_grey_markers_tfOnly_nUpDown_w3P_wMarkers, iCard_up_hist, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()

# grey points with 7UP and markers
f7UP <- c('SP3', 'ZBTB10', 'ZBTB44', 'SP3', 'SSH2', 'ZNF770', 'ZFP91')
set.seed(3561)
iCard_grey_markers_tfOnly_nUpDown_w7UP_wMarkers <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% inner_join(geneIDtoGeneSymbol) %>%
                filter(cellType == "iCard", meanRPM > 20, !(GeneSymbol %in% f7UP)) %>% anti_join(markerTab2, by = "gene_id"), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), alpha = 0.3) +
  geom_point(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), size = 2, color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp, label = GeneSymbol), color = "red", size = 6) +
  geom_text_repel(data = indiv_perturbability %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
                    filter(cellType == "iCard", GeneSymbol %in% f7UP) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp, label = GeneSymbol), color = "magenta", size = 6) +
  geom_point(data = indiv_perturbability %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
               filter(cellType == "iCard", GeneSymbol %in% f7UP) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), size = 2, color = "magenta") +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  ylim(c(-0.5, 12.5)) +
  xlab("Number of conditions DOWN") + ylab("Number of conditions UP")
# plot(iCard_grey_markers_tfOnly_nUpDown_w3P_wMarkers)
ggsave(iCard_grey_markers_tfOnly_nUpDown_w7UP_wMarkers, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_grey_markers_tfOnly_nUpDown_w7UP_wMarkers.pdf'), width = 12, height = 7, useDingbats = F)

# pdf(paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_grey_markers_tfOnly_nUpDown_withMarginals_w3P_wMarkers.pdf'), width = 9, height = 6, useDingbats = F)
# grid.arrange(iCard_down_hist, empty, iCard_grey_markers_tfOnly_nUpDown_w7UP_wMarkers, iCard_up_hist, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
# dev.off()


# filter to predicted LOESS >= 1, corresponding to mean RPM = e^3 = 20
set.seed(3561)
iCard_grey_markers_tfOnly_nUpDown_w3P_wMarkers_RPMfilt <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard", meanRPM > 20) %>% anti_join(markerTab, by = "gene_id"), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), alpha = 0.3) +
  geom_point(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), size = 2, color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp, label = GeneSymbol), color = "red", size = 6) +
  geom_text_repel(data = indiv_perturbability %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
                    filter(cellType == "iCard", GeneSymbol %in% c('SP3', 'ZBTB10', 'ZBTB44'), meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp, label = GeneSymbol), color = "blue", size = 6) +
  geom_point(data = indiv_perturbability %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
               filter(cellType == "iCard", GeneSymbol %in% c('SP3', 'ZBTB10', 'ZBTB44'), meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), size = 2, color = "blue") +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  ylim(c(-0.5, 12.5)) +
  xlab("Number of conditions DOWN") + ylab("Number of conditions UP")
# plot(iCard_grey_markers_tfOnly_nUpDown_w3P_wMarkers_RPMfilt)

iCard_grey_markers_tfOnly_nUpDown_wMarkers_RPMfilt <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard", meanRPM > 20) %>% anti_join(markerTab2, by = "gene_id"), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), alpha = 0.3) +
  geom_point(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), size = 2, color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp, label = GeneSymbol), color = "red", size = 6) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  ylim(c(-0.5, 12.5)) +
  xlab("Number of conditions DOWN") + ylab("Number of conditions UP")
# plot(iCard_grey_markers_tfOnly_nUpDown_wMarkers_RPMfilt)

iCard_up_hist2_RPMfilt <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20), aes(nDESeq2conditionsUp),
                 binwidth = 1) +
  geom_vline(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(xintercept = nDESeq2conditionsUp), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(x = nDESeq2conditionsUp, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() + 
  xlim(c(-0.5, 12.5)) +
  xlab("Number of conditions UP") +
  coord_flip()

iCard_down_hist2_RPMfilt <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20), aes(nDESeq2conditionsDown),
                 binwidth = 1) +
  geom_vline() +
  geom_vline(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(xintercept = nDESeq2conditionsDown), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(x = nDESeq2conditionsDown, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() +
  xlab("Number of conditions DOWN")

pdf(paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_grey_markers_tfOnly_RPMfilt_nUpDown_withMarginals_w3P_wMarkers.pdf'), width = 9, height = 6, useDingbats = F)
grid.arrange(iCard_down_hist2_RPMfilt, empty, iCard_grey_markers_tfOnly_nUpDown_w3P_wMarkers_RPMfilt, iCard_up_hist2_RPMfilt, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()

pdf(paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_grey_markers_tfOnly_RPMfilt_nUpDown_withMarginals_wMarkers.pdf'), width = 9, height = 6, useDingbats = F)
grid.arrange(iCard_down_hist2_RPMfilt, empty, iCard_grey_markers_tfOnly_nUpDown_wMarkers_RPMfilt, iCard_up_hist2_RPMfilt, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()

iCard_grey_markers_tfOnly_nUpDown_wMarkerDots_RPMfilt_ZOpoints <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard", meanRPM > 20) %>% anti_join(markerTab2, by = "gene_id"), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), alpha = 0.3) +
  geom_point(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), size = 2, color = "red") +
  # geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
  #                   filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp, label = GeneSymbol), color = "red", size = 6) +
  geom_point(data = zhouOlson_activators %>% inner_join(geneIDtoGeneSymbol) %>% inner_join(indiv_perturbability, by = 'gene_id') %>%
               filter(cellType == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), size = 2, color = "black") +
  geom_text_repel(data = zhouOlson_activators %>% inner_join(geneIDtoGeneSymbol) %>% inner_join(indiv_perturbability, by = 'gene_id') %>%
                    filter(cellType == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp, label = GeneSymbol), color = "black", size = 6) +
  geom_point(data = zhouOlson_inhibitors %>% head(48) %>% inner_join(geneIDtoGeneSymbol) %>% inner_join(indiv_perturbability, by = 'gene_id') %>%
               filter(cellType == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), size = 2, color = "yellow") +
  geom_text_repel(data = zhouOlson_inhibitors %>% head(48) %>% inner_join(geneIDtoGeneSymbol) %>% inner_join(indiv_perturbability, by = 'gene_id') %>%
                    filter(cellType == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp, label = GeneSymbol), color = "yellow", size = 6) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  ylim(c(-0.5, 12.5)) +
  xlab("Number of conditions DOWN") + ylab("Number of conditions UP")

iCard_grey_markers_tfOnly_nUpDown_wMarkerDots_RPMfilt_ZOcontours <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard", meanRPM > 20) %>% anti_join(markerTab2, by = "gene_id"), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), alpha = 0.3) +
  geom_point(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), size = 2, color = "red") +
  # geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
  #                   filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp, label = GeneSymbol), color = "red", size = 6) +
  geom_density_2d(data = zhouOlson_activators %>% inner_join(geneIDtoGeneSymbol) %>% inner_join(indiv_perturbability, by = 'gene_id') %>%
               filter(cellType == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), size = 2, color = "black") +
  geom_density_2d(data = zhouOlson_inhibitors %>% head(48) %>% inner_join(geneIDtoGeneSymbol) %>% inner_join(indiv_perturbability, by = 'gene_id') %>%
               filter(cellType == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), size = 2, color = "yellow") +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  ylim(c(-0.5, 12.5)) +
  xlab("Number of conditions DOWN") + ylab("Number of conditions UP")


## LFC filter 
# no filter
iCard_up_hist_RPMfilt_lfc0 <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20), aes(nDESeq2conditionsUp),
                 binwidth = 1) +
  geom_vline(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(xintercept = nDESeq2conditionsUp), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(x = nDESeq2conditionsUp, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() + 
  xlim(c(-0.5, 12.5)) +
  xlab("Number of conditions UP") +
  coord_flip()

iCard_down_hist_RPMfilt_lfc0 <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20), aes(nDESeq2conditionsDown),
                 binwidth = 1) +
  geom_vline() +
  geom_vline(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(xintercept = nDESeq2conditionsDown), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(x = nDESeq2conditionsDown, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() +
  xlab("Number of conditions DOWN")

set.seed(3561)
iCard_grey_markers_tfOnly_nUpDown_lfc0 <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard", meanRPM > 20) %>% anti_join(markerTab2, by = "gene_id"), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), alpha = 0.3) +
  geom_point(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown, nDESeq2conditionsUp, label = GeneSymbol), color = "red", size = 6) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  ylim(c(-0.5, 12.5)) +
  xlab("Number of conditions DOWN") + ylab("Number of conditions UP") +
  ggtitle('Minimum |LFC| = 0')
# plot(iCard_grey_markers_tfOnly_nUpDown_lfc0)
ggsave(iCard_grey_markers_tfOnly_nUpDown_lfc0, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_grey_markers_tfOnly_nUpDown_lfc0_5.pdf'), width = 7, height = 7, useDingbats = F)

pdf(paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_grey_markers_tfOnly_nUpDown_withMarginals_lfc0.pdf'), width = 9, height = 6, useDingbats = F)
grid.arrange(iCard_down_hist_RPMfilt_lfc0, empty, iCard_grey_markers_tfOnly_nUpDown_lfc0, iCard_up_hist_RPMfilt_lfc0, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()


iCard_up_hist_RPMfilt_lfc0_5 <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20), aes(nDESeq2conditionsUp_0_5),
                 binwidth = 1) +
  geom_vline(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(xintercept = nDESeq2conditionsUp_0_5), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(x = nDESeq2conditionsUp_0_5, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() + 
  xlim(c(-0.5, 12.5)) +
  xlab("Number of conditions UP") +
  coord_flip()

iCard_down_hist_RPMfilt_lfc0_5 <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20), aes(nDESeq2conditionsDown_0_5),
                 binwidth = 1) +
  geom_vline() +
  geom_vline(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(xintercept = nDESeq2conditionsDown_0_5), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(x = nDESeq2conditionsDown_0_5, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() +
  xlab("Number of conditions DOWN")

set.seed(3561)
iCard_grey_markers_tfOnly_nUpDown_lfc0_5 <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard", meanRPM > 20) %>% anti_join(markerTab2, by = "gene_id"), aes(nDESeq2conditionsDown_0_5, nDESeq2conditionsUp_0_5), alpha = 0.3) +
  geom_point(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown_0_5, nDESeq2conditionsUp_0_5), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown_0_5, nDESeq2conditionsUp_0_5, label = GeneSymbol), color = "red", size = 6) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  ylim(c(-0.5, 12.5)) +
  xlab("Number of conditions DOWN") + ylab("Number of conditions UP") +
  ggtitle('Minimum |LFC| = 0.5')
# plot(iCard_grey_markers_tfOnly_nUpDown_lfc0_5)
ggsave(iCard_grey_markers_tfOnly_nUpDown_lfc0_5, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_grey_markers_tfOnly_nUpDown_lfc0_5.pdf'), width = 7, height = 7, useDingbats = F)

pdf(paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_grey_markers_tfOnly_nUpDown_withMarginals_lfc0_5.pdf'), width = 9, height = 6, useDingbats = F)
grid.arrange(iCard_down_hist_RPMfilt_lfc0_5, empty, iCard_grey_markers_tfOnly_nUpDown_lfc0_5, iCard_up_hist_RPMfilt_lfc0_5, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()


# LFC > 0.75
iCard_up_hist_RPMfilt_lfc0_75 <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20), aes(nDESeq2conditionsUp_0_75),
                 binwidth = 1) +
  geom_vline(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(xintercept = nDESeq2conditionsUp_0_75), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(x = nDESeq2conditionsUp_0_75, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() + 
  xlim(c(-0.5, 12.5)) +
  xlab("Number of conditions UP") +
  coord_flip()

iCard_down_hist_RPMfilt_lfc0_75 <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20), aes(nDESeq2conditionsDown_0_75),
                 binwidth = 1) +
  geom_vline() +
  geom_vline(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(xintercept = nDESeq2conditionsDown_0_75), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(x = nDESeq2conditionsDown_0_75, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() +
  xlab("Number of conditions DOWN")

set.seed(3561)
iCard_grey_markers_tfOnly_nUpDown_lfc0_75 <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard", meanRPM > 20) %>% anti_join(markerTab2, by = "gene_id"), aes(nDESeq2conditionsDown_0_75, nDESeq2conditionsUp_0_75), alpha = 0.3) +
  geom_point(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown_0_75, nDESeq2conditionsUp_0_75), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDown_0_75, nDESeq2conditionsUp_0_75, label = GeneSymbol), color = "red", size = 6) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  ylim(c(-0.5, 12.5)) +
  xlab("Number of conditions DOWN") + ylab("Number of conditions UP") +
  ggtitle('Minimum |LFC| = 0.75')
# plot(iCard_grey_markers_tfOnly_nUpDown_lfc0_75)
ggsave(iCard_grey_markers_tfOnly_nUpDown_lfc0_75, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_grey_markers_tfOnly_nUpDown_lfc0_75.pdf'), width = 7, height = 7, useDingbats = F)

pdf(paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_grey_markers_tfOnly_nUpDown_withMarginals_lfc0_75.pdf'), width = 9, height = 6, useDingbats = F)
grid.arrange(iCard_down_hist_RPMfilt_lfc0_75, empty, iCard_grey_markers_tfOnly_nUpDown_lfc0_75, iCard_up_hist_RPMfilt_lfc0_75, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()


# per-plate results
iCard_DEres_perPlate <- as_tibble(read.table(paste0(graphDir, "/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_DESeqResults_perPlate.txt"), sep = "\t", stringsAsFactors = F, header = T))
iCard_upDown_perPlate <- iCard_DEres_perPlate %>%
  filter(padj <= 0.1) %>%
  group_by(gene_id) %>%
  summarise(nDESeq2conditionsUpPP = sum(log2FoldChange > 0),
            nDESeq2conditionsDownPP = sum(log2FoldChange < 0),
            nDESeq2conditionsAllPP = nDESeq2conditionsUpPP + nDESeq2conditionsDownPP)

iCard_up_hist_RPMfilt_PP <- ggplot() +
  geom_histogram(data = iCard_upDown_perPlate %>% inner_join(indiv_perturbability, by = 'gene_id') %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20), aes(nDESeq2conditionsUpPP),
                 binwidth = 1) +
  geom_vline(data = iCard_upDown_perPlate %>% inner_join(indiv_perturbability, by = 'gene_id') %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(xintercept = nDESeq2conditionsUpPP), color = "red") +
  geom_text_repel(data = iCard_upDown_perPlate %>% inner_join(indiv_perturbability, by = 'gene_id') %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(x = nDESeq2conditionsUpPP, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() + 
  xlim(c(-0.5, 12.5)) +
  xlab("Number of conditions UP") +
  coord_flip()

iCard_down_hist_RPMfilt_PP <- ggplot() +
  geom_histogram(data = iCard_upDown_perPlate %>% inner_join(indiv_perturbability, by = 'gene_id') %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20), aes(nDESeq2conditionsDownPP),
                 binwidth = 1) +
  geom_vline() +
  geom_vline(data = iCard_upDown_perPlate %>% inner_join(indiv_perturbability, by = 'gene_id') %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(xintercept = nDESeq2conditionsDownPP), color = "red") +
  geom_text_repel(data = iCard_upDown_perPlate %>% inner_join(indiv_perturbability, by = 'gene_id') %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(x = nDESeq2conditionsDownPP, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() +
  xlab("Number of conditions DOWN")

set.seed(3561)
iCard_grey_markers_tfOnly_nUpDown_perPlate <- ggplot() +
  geom_jitter(data = iCard_upDown_perPlate %>% inner_join(tfTab, by = "gene_id") %>% inner_join(indiv_perturbability, by = 'gene_id') %>%
                filter(cellType == "iCard", meanRPM > 20) %>% anti_join(markerTab2, by = "gene_id"), aes(nDESeq2conditionsDownPP, nDESeq2conditionsUpPP), alpha = 0.3) +
  geom_point(data = iCard_upDown_perPlate %>% inner_join(markerTab2, by = "gene_id") %>% inner_join(indiv_perturbability, by = 'gene_id') %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDownPP, nDESeq2conditionsUpPP), color = "red") +
  geom_text_repel(data = iCard_upDown_perPlate %>% inner_join(markerTab2, by = "gene_id") %>% inner_join(indiv_perturbability, by = 'gene_id') %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(nDESeq2conditionsDownPP, nDESeq2conditionsUpPP, label = GeneSymbol), color = "red", size = 6) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  ylim(c(-0.5, 12.5)) +
  xlab("Number of conditions DOWN") + ylab("Number of conditions UP") +
  ggtitle('Per-plate analysis (fewer controls per contrast)')
pdf(paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_grey_markers_tfOnly_nUpDown_perPlate.pdf'), width = 9, height = 6, useDingbats = F)
grid.arrange(iCard_down_hist_RPMfilt_PP, empty, iCard_grey_markers_tfOnly_nUpDown_perPlate, iCard_up_hist_RPMfilt_PP, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()

## fibroblasts ####

fibro_markers_tfOnly_nUpDown <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "GM00942", meanTPM > 10) %>% anti_join(markerTab, by = "gene_id"), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), color = "blue", alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  ylim(c(-0.5, 12.5)) + xlim(c(-0.5, 20.5)) +
  xlab("Number of conditions DOWN") + ylab("Number of conditions UP")
# plot(iCard_markers_tfOnly_nUpDown)
# ggsave(iCard_markers_tfOnly_nUpDown, file = paste0(graphDir, '/differentialExpression/allSamples/fibroblasts/fibro_tfOnly_nUpDown.pdf'), width = 7, height = 7, useDingbats = F)

pdf(paste0(graphDir, '/differentialExpression/allSamples/fibroblasts/fibro_tfOnly_nUpDown.pdf'), width = 15, height = 10, useDingbats = F)
grid.arrange(empty, empty, fibro_markers_tfOnly_nUpDown, empty, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()

# grey
fibro_grey_tfOnly_nUpDown <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "GM00942", meanTPM > 10) %>% anti_join(markerTab, by = "gene_id"), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), color = "black", alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  ylim(c(-0.5, 12.5)) + xlim(c(-0.5, 20.5)) +
  xlab("Number of conditions DOWN") + ylab("Number of conditions UP")
# plot(fibro_grey_tfOnly_nUpDown)
# ggsave(iCard_markers_tfOnly_nUpDown, file = paste0(graphDir, '/differentialExpression/allSamples/fibroblasts/fibro_tfOnly_nUpDown.pdf'), width = 7, height = 7, useDingbats = F)

pdf(paste0(graphDir, '/differentialExpression/allSamples/fibroblasts/fibro_grey_tfOnly_nUpDown.pdf'), width = 9, height = 6, useDingbats = F)
grid.arrange(empty, empty, fibro_grey_tfOnly_nUpDown, empty, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()


fibro_grey_tfOnly_nUpDown_w4targ <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "GM00942", meanTPM > 10) %>% anti_join(markerTab, by = "gene_id"), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), color = "black", alpha = 0.3) +
  geom_point(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% inner_join(geneIDtoGeneSymbol) %>%
               filter(cellType == "GM00942", meanTPM > 10, GeneSymbol %in% c('CEBPB', 'ID3', 'PRRX2', 'RUNX1')) %>% anti_join(markerTab, by = "gene_id"), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), color = "blue") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% inner_join(geneIDtoGeneSymbol) %>%
                    filter(cellType == "GM00942", meanTPM > 10, GeneSymbol %in% c('CEBPB', 'ID3', 'PRRX2', 'RUNX1')) %>% anti_join(markerTab, by = "gene_id"), aes(nDESeq2conditionsDown, nDESeq2conditionsUp, label = GeneSymbol), color = "blue", size = 6) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  # ylim(c(-0.5, 12.5)) + xlim(c(-0.5, 20.5)) +
  xlab("Number of conditions DOWN") + ylab("Number of conditions UP")
# plot(fibro_grey_tfOnly_nUpDown_w4targ)
# ggsave(iCard_markers_tfOnly_nUpDown, file = paste0(graphDir, '/differentialExpression/allSamples/fibroblasts/fibro_tfOnly_nUpDown.pdf'), width = 7, height = 7, useDingbats = F)

targs16 <- c('SKIL', 'YBX3', 'TSC22D1', 'CERS2', 'KLF13', 'TBX3', 'ID1', 'ATOH8', 'ZNF652', 'NFATC4', 'ZBTB38', 'LARP1', 'CEBPB', 'ID3', 'PRRX2', 'RUNX1')
targs8 <- c('CERS2', 'KLF13', 'ATOH8', 'ZNF652', 'ZBTB38', 'CEBPB', 'PRRX2', 'RUNX1')

fibro_grey_tfOnly_nUpDown_w16targ <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "GM00942", meanTPM > 10) %>% anti_join(markerTab, by = "gene_id"), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), color = "black", alpha = 0.3) +
  geom_point(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% inner_join(geneIDtoGeneSymbol) %>%
               filter(cellType == "GM00942", meanTPM > 10, GeneSymbol %in% targs16) %>% anti_join(markerTab, by = "gene_id"), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), color = "blue") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% inner_join(geneIDtoGeneSymbol) %>%
                    filter(cellType == "GM00942", meanTPM > 10, GeneSymbol %in% targs16) %>% anti_join(markerTab, by = "gene_id"), aes(nDESeq2conditionsDown, nDESeq2conditionsUp, label = GeneSymbol), color = "blue", size = 6) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  # ylim(c(-0.5, 12.5)) + xlim(c(-0.5, 20.5)) +
  xlab("Number of conditions DOWN") + ylab("Number of conditions UP")
# plot(fibro_grey_tfOnly_nUpDown_w16targ)

fibro_grey_tfOnly_nUpDown_w8targ <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "GM00942", meanTPM > 10) %>% anti_join(markerTab, by = "gene_id"), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), color = "black", alpha = 0.3) +
  geom_point(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% inner_join(geneIDtoGeneSymbol) %>%
               filter(cellType == "GM00942", meanTPM > 10, GeneSymbol %in% targs8) %>% anti_join(markerTab, by = "gene_id"), aes(nDESeq2conditionsDown, nDESeq2conditionsUp), color = "blue") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% inner_join(geneIDtoGeneSymbol) %>%
                    filter(cellType == "GM00942", meanTPM > 10, GeneSymbol %in% targs8) %>% anti_join(markerTab, by = "gene_id"), aes(nDESeq2conditionsDown, nDESeq2conditionsUp, label = GeneSymbol), color = "blue", size = 6) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  # ylim(c(-0.5, 12.5)) + xlim(c(-0.5, 20.5)) +
  xlab("Number of conditions\nin which DOWN-regulated") + ylab("Number of conditions\nin which UP-regulated")
# plot(fibro_grey_tfOnly_nUpDown_w8targ)


pdf(paste0(graphDir, '/differentialExpression/allSamples/fibroblasts/fibro_grey_tfOnly_nUpDown_w8targ.pdf'), width = 9, height = 6, useDingbats = F)
grid.arrange(empty, empty, fibro_grey_tfOnly_nUpDown_w8targ, empty, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()


set.seed(96438)
fibro_TFs_nUpVsExp_noLOESS_vertJit_wTested <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% 
                inner_join(tfTab, by = "gene_id") %>%
                inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
                filter(cellType == "GM00942",
                       !(GeneSymbol %in% targs16)) %>% 
                mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
              aes(log(meanRPM), nDESeq2conditionsUp), width = 0, height = 0.2, alpha = 0.5) +
  geom_point(data = indiv_perturbability %>% 
               inner_join(tfTab, by = "gene_id") %>%
               inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
               filter(cellType == "GM00942", GeneSymbol %in% targs16) %>%
               mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
             aes(log(meanRPM), nDESeq2conditionsUp), color = "forestgreen") +
  geom_text_repel(data = indiv_perturbability %>% 
                    inner_join(tfTab, by = "gene_id") %>%
                    inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
                    filter(cellType == "GM00942", GeneSymbol %in% targs16) %>%
                    mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
                  aes(log(meanRPM), nDESeq2conditionsUp, label = GeneSymbol), color = "forestgreen", nudge_y = 0.1, size = 6) +
  xlim(c(0,8)) +
  theme_bw() +
  theme(axis.text = element_text(size = rel(2)),
        axis.title = element_blank())
# plot(iCard_TFs_nUpVsExp_noLOESS_vertJit_wMarkers)
setwd(projectDir)
ggsave(fibro_TFs_nUpVsExp_noLOESS_vertJit_wTested, file = paste0(graphDir, '/differentialExpression/allSamples/fibroblasts/fibro_TFs_nUpVsExp_noLOESS_vertJit_wTested.pdf'), width = 9, height = 7, useDingbats = F)



## upregulation stats ####
set.seed(3561)
iCard_upreg_effectSize_nUp_logMeanPCH <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                mutate(up_meanLFC = ifelse(is.na(up_meanLFC), 0, up_meanLFC)) %>%
                filter(cellType == "iCard", meanTPM > 50) %>% anti_join(markerTab, by = "gene_id"), aes(nDESeq2conditionsUp, up_meanLFC, size = log(meanTPM)), color = "blue", alpha = 0.3, height = 0) +
  geom_point(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
               mutate(up_meanLFC = ifelse(is.na(up_meanLFC), 0, up_meanLFC)) %>%
               filter(cellType == "iCard", cellTypeMark == "iCard") %>% unique(), aes(nDESeq2conditionsUp, up_meanLFC, size = log(meanTPM)), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
                    mutate(up_meanLFC = ifelse(is.na(up_meanLFC), 0, up_meanLFC)) %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard") %>% unique(), aes(nDESeq2conditionsUp, up_meanLFC, size = log(meanTPM), label = GeneSymbol), color = "red") +
  theme_bw() +
  xlab("Number of conditions UP") + ylab("Mean log2FC when UP") + ggtitle("TF genes in iCards, meanTPM > 50\nUp-regulation stats")
# plot(iCard_upreg_effectSize_nUp_logMeanPCH)
ggsave(iCard_upreg_effectSize_nUp_logMeanPCH, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_upreg_effectSize_nUp_logMeanPCH.pdf'), width = 7, height = 7, useDingbats = F)

set.seed(3561)
iCard_upreg_effectSize_nUp_logMeanPCH_withNames <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                mutate(up_meanLFC = ifelse(is.na(up_meanLFC), 0, up_meanLFC)) %>%
                filter(cellType == "iCard", meanTPM > 50) %>% anti_join(markerTab, by = "gene_id"), aes(nDESeq2conditionsUp, up_meanLFC, size = log(meanTPM)), color = "blue", alpha = 0.3, height = 0) +
  geom_text_repel(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
                    mutate(up_meanLFC = ifelse(is.na(up_meanLFC), 0, up_meanLFC)) %>% 
                    filter(cellType == "iCard", meanTPM > 50, nDESeq2conditionsUp >= 6) %>% anti_join(markerTab, by = "gene_id"), aes(nDESeq2conditionsUp, up_meanLFC, size = log(meanTPM), label = GeneSymbol), color = "blue", alpha = 0.6) +
  geom_point(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
               mutate(up_meanLFC = ifelse(is.na(up_meanLFC), 0, up_meanLFC)) %>%
               filter(cellType == "iCard", cellTypeMark == "iCard") %>% unique(), aes(nDESeq2conditionsUp, up_meanLFC, size = log(meanTPM)), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
                    mutate(up_meanLFC = ifelse(is.na(up_meanLFC), 0, up_meanLFC)) %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard") %>% unique(), aes(nDESeq2conditionsUp, up_meanLFC, size = log(meanTPM), label = GeneSymbol), color = "red") +
  theme_bw() +
  xlab("Number of conditions UP") + ylab("Mean log2FC when UP") + ggtitle("TF genes in iCards, meanTPM > 50\nUp-regulation stats")
# plot(iCard_upreg_effectSize_nUp_logMeanPCH_withNames)


set.seed(3561)
iCard_upreg_effectSize_nUp_logMeanPCH_withNames <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                mutate(up_meanLFC = ifelse(is.na(up_meanLFC), 0, up_meanLFC)) %>%
                filter(cellType == "iCard", meanTPM > 50) %>% anti_join(markerTab, by = "gene_id"), aes(nDESeq2conditionsUp, up_meanLFC, size = log(meanTPM)), color = "blue", alpha = 0.3, height = 0) +
  geom_text_repel(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
                    mutate(up_meanLFC = ifelse(is.na(up_meanLFC), 0, up_meanLFC)) %>% 
                    filter(cellType == "iCard", meanTPM > 50, nDESeq2conditionsUp >= 6) %>% anti_join(markerTab, by = "gene_id"), aes(nDESeq2conditionsUp, up_meanLFC, size = log(meanTPM), label = GeneSymbol), color = "blue", alpha = 0.6) +
  geom_point(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
               mutate(up_meanLFC = ifelse(is.na(up_meanLFC), 0, up_meanLFC)) %>%
               filter(cellType == "iCard", cellTypeMark == "iCard") %>% unique(), aes(nDESeq2conditionsUp, up_meanLFC, size = log(meanTPM)), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
                    mutate(up_meanLFC = ifelse(is.na(up_meanLFC), 0, up_meanLFC)) %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard") %>% unique(), aes(nDESeq2conditionsUp, up_meanLFC, size = log(meanTPM), label = GeneSymbol), color = "red") +
  theme_bw() +
  xlab("Number of conditions UP") + ylab("Mean log2FC when UP") + ggtitle("TF genes in iCards, meanTPM > 50\nUp-regulation stats")
# plot(iCard_upreg_effectSize_nUp_logMeanPCH_withNames)


iCard_upreg_effectSize_freqUp_logMeanPCH <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard", meanTPM > 50) %>% anti_join(markerTab, by = "gene_id"), aes(up_nLFC/(up_nLFC+down_nLFC), up_meanLFC, size = log(meanTPM)), color = "blue", alpha = 0.3, height = 0) +
  geom_point(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard") %>% unique(), aes(up_nLFC/(up_nLFC+down_nLFC), up_meanLFC, size = log(meanTPM)), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard") %>% unique(), aes(up_nLFC/(up_nLFC+down_nLFC), up_meanLFC, size = log(meanTPM), label = GeneSymbol), color = "red") +
  theme_bw() +
  xlab("Number of conditions UP") + ylab("Mean log2FC when UP")
# plot(iCard_upreg_effectSize_freqUp_logMeanPCH)

set.seed(8375)
iCard_upreg_effectSize_nUp_grey <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard", meanRPM > 20) %>% anti_join(markerTab, by = "gene_id"), aes(up_nLFC, up_meanLFC), alpha = 0.3, height = 0) +
  geom_point(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard") %>% unique(), aes(up_nLFC, up_meanLFC), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard") %>% unique(), aes(up_nLFC, up_meanLFC, label = GeneSymbol), color = "red", size = 6) +
  theme_bw() +
  theme(axis.text = element_text(size = rel(2)),
        axis.title = element_text(size = rel(2))) +
  xlab("Number of conditions UP") + ylab("Mean log2FC when UP")
# plot(iCard_upreg_effectSize_nUp_grey)
ggsave(iCard_upreg_effectSize_nUp_grey, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_upreg_effectSize_nUp.pdf'), width = 7, height = 7, useDingbats = F)


iCard_upreg_effectSize_nUp_logMeanPCH_contour <- ggplot() +
  geom_density_2d(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                    mutate(up_meanLFC = ifelse(is.na(up_meanLFC), 0, up_meanLFC)) %>%
                    filter(cellType == "iCard", meanTPM > 50) %>% anti_join(markerTab, by = "gene_id"), aes(up_nLFC, up_meanLFC, size = log(meanTPM)), color = "blue", alpha = 0.3) +
  geom_point(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
               mutate(up_meanLFC = ifelse(is.na(up_meanLFC), 0, up_meanLFC)) %>%
               filter(cellType == "iCard", cellTypeMark == "iCard") %>% unique(), aes(up_nLFC, up_meanLFC, size = log(meanTPM)), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
                    mutate(up_meanLFC = ifelse(is.na(up_meanLFC), 0, up_meanLFC)) %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard") %>% unique(), aes(up_nLFC, up_meanLFC, size = log(meanTPM), label = GeneSymbol), color = "red") +
  theme_bw() +
  xlab("Number of conditions UP") + ylab("Mean log2FC when UP") + ggtitle("TF genes in iCards, meanTPM > 50\nUp-regulation stats")
# plot(iCard_upreg_effectSize_nUp_logMeanPCH_contour)
ggsave(iCard_upreg_effectSize_nUp_logMeanPCH_contour, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard__upreg_effectSize_nUp_logMeanPCH_contour.pdf'), width = 7, height = 7, useDingbats = F)


iCard_upreg_effectSize_freqUp_logMeanPCH_contour <- ggplot() +
  geom_density_2d(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard", meanTPM > 50) %>% anti_join(markerTab, by = "gene_id"), aes(up_nLFC/(up_nLFC+down_nLFC), up_meanLFC, size = log(meanTPM)), color = "blue", alpha = 0.3) +
  geom_point(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard") %>% unique(), aes(up_nLFC/(up_nLFC+down_nLFC), up_meanLFC, size = log(meanTPM)), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard") %>% unique(), aes(up_nLFC/(up_nLFC+down_nLFC), up_meanLFC, size = log(meanTPM), label = GeneSymbol), color = "red") +
  theme_bw() +
  xlab("Fraction of conditions UP") + ylab("Mean log2FC when UP")
# plot(iCard_upreg_effectSize_freqUp_logMeanPCH_contour)


iCard_TFs_highNUP_highMean_highFracUP <- indiv_perturbability %>% 
  inner_join(tfTab, by = "gene_id") %>%
  mutate(up_meanLFC = ifelse(is.na(up_meanLFC), 0, up_meanLFC)) %>%
  filter(cellType == "iCard", meanTPM > 50, nDESeq2conditionsUp >= 4) %>% 
  anti_join(markerTab, by = "gene_id") %>% 
  inner_join(geneIDtoGeneSymbol, by = "gene_id") %>% 
  mutate(fracDESeq2up = nDESeq2conditionsUp / nDESeq2conditionsAll) %>% 
  dplyr::select(GeneSymbol, meanTPM, nDESeq2conditionsUp, fracDESeq2up) %>%
  arrange(-nDESeq2conditionsUp, -fracDESeq2up) %>%
  mutate(inOlsonScreen = ifelse(GeneSymbol %in% zhouOlson_allORFs$GeneSymbol, 'tested in mice', 'untested'))

iCard_TFs_highNUP_highMean_highFracUP_top30 = iCard_TFs_highNUP_highMean_highFracUP[1:30,] %>%
  dplyr::select(GeneSymbol) %>% inner_join(geneIDtoGeneSymbol, by = "GeneSymbol")
markerGeneSymbols = markerTab %>% filter(cellTypeMark == "iCard") %>% dplyr::select(GeneSymbol, gene_id)

write.table(bind_rows(iCard_TFs_highNUP_highMean_highFracUP_top30, markerGeneSymbols), file = paste0(graphDir, "/differentialExpression/allSamples/iCards/iCard_markers_AndCandidatesTop30_GeneSymbolsAndIDs.txt"), row.names = F, quote = F)

iCard_GTEx_indiv_nDESeq2ConditionsUp_cdf <- ggplot() + 
  stat_ecdf(data = indiv_perturbability %>%
              inner_join(tfTab, by = 'gene_id') %>%
              filter(cellType == "iCard", meanRPM > 20), 
            aes(nDESeq2conditionsUp), alpha = 0.2) +
  theme_bw() + 
  xlab('nDESeq2conditionsUp') +
  ylab('cumulative distribution function') +
  ggtitle("Responsiveness to perturbation (up-regulation)\nCumulative distribution function\nAll TF genes >20 RPM in control iCards")
ggsave(iCard_GTEx_indiv_nDESeq2ConditionsUp_cdf, file = paste0(graphDir, "/differentialExpression/allSamples/iCards/iCard_GTEx_indiv_nDESeq2ConditionsUp_cdf.pdf"), width = 8, height = 8, useDingbats = F)


# up alone and down alone vs expression ####
iCard_TFs_fracUpVsExp <- ggplot() +
  geom_point(data = indiv_perturbability %>% 
               inner_join(tfTab, by = "gene_id") %>%
               filter(cellType == "iCard", nDESeq2conditionsAll > 0) %>% 
               anti_join(markerTab, by = "gene_id") %>%
               mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
                      aes(log(meanRPM), fracUp)) +
  geom_smooth(data = indiv_perturbability %>% 
                inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard", nDESeq2conditionsAll > 0) %>% 
                anti_join(markerTab, by = "gene_id") %>%
                mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
              aes(log(meanRPM), fracUp)) +
  geom_point(data = indiv_perturbability %>% 
               inner_join(tfTab, by = "gene_id") %>%
               filter(cellType == "iCard", nDESeq2conditionsAll > 0) %>% 
               inner_join(markerTab, by = "gene_id") %>%
               mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
             aes(log(meanRPM), fracUp), color = "red")

iCard_TFs_nUpVsExp_wLOESS_wMarkers <- ggplot() +
  geom_point(data = indiv_perturbability %>% 
               inner_join(tfTab, by = "gene_id") %>%
               inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
               filter(cellType == "iCard") %>% 
               anti_join(markerTab, by = "gene_id") %>%
               filter(!(GeneSymbol %in% c('ZFPM2', 'ESRRG'))) %>%
               mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
             aes(log(meanRPM), nDESeq2conditionsUp)) +
  geom_smooth(data = indiv_perturbability %>% 
                inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard") %>% 
                mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
              aes(log(meanRPM), nDESeq2conditionsUp),
              method = 'loess', color = 'black') +
  geom_point(data = indiv_perturbability %>% 
               inner_join(tfTab, by = "gene_id") %>%
               inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
               filter(cellType == "iCard", gene_id %in% markerTab$gene_id) %>%
               mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
             aes(log(meanRPM), nDESeq2conditionsUp), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% 
                    inner_join(tfTab, by = "gene_id") %>%
                    inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
                    filter(cellType == "iCard", gene_id %in% markerTab$gene_id) %>%
                    mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
             aes(log(meanRPM), nDESeq2conditionsUp, label = GeneSymbol), color = "red", nudge_y = 0.1, size = 6) +
  xlim(c(0,8)) +
  theme_bw() +
  theme(axis.text = element_text(size = rel(2)),
        axis.title = element_blank())
# plot(iCard_TFs_nUpVsExp_wLOESS_wMarkers)
ggsave(iCard_TFs_nUpVsExp_wLOESS_wMarkers, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_TFs_nUpVsExp_wLOESS_wMarkers.pdf'), width = 9, height = 7, useDingbats = F)

set.seed(96438)
iCard_TFs_nUpVsExp_noLOESS_vertJit_wMarkers <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% 
               inner_join(tfTab, by = "gene_id") %>%
               inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
               filter(cellType == "iCard") %>% 
               anti_join(markerTab2, by = "gene_id") %>%
               filter(!(GeneSymbol %in% c('ZFPM2', 'ESRRG'))) %>%
               mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
             aes(log(meanRPM), nDESeq2conditionsUp), width = 0, height = 0.2, alpha = 0.5) +
  geom_point(data = indiv_perturbability %>% 
               inner_join(tfTab, by = "gene_id") %>%
               inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
               filter(cellType == "iCard", gene_id %in% markerTab2$gene_id) %>%
               mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
             aes(log(meanRPM), nDESeq2conditionsUp), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% 
                    inner_join(tfTab, by = "gene_id") %>%
                    inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
                    filter(cellType == "iCard", gene_id %in% markerTab2$gene_id) %>%
                    mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
                  aes(log(meanRPM), nDESeq2conditionsUp, label = GeneSymbol), color = "red", nudge_y = 0.1, size = 6) +
  xlim(c(0,8)) +
  theme_bw() +
  theme(axis.text = element_text(size = rel(2)),
        axis.title = element_blank())
# plot(iCard_TFs_nUpVsExp_noLOESS_vertJit_wMarkers)
setwd(projectDir)
ggsave(iCard_TFs_nUpVsExp_noLOESS_vertJit_wMarkers, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_TFs_nUpVsExp_noLOESS_vertJit_wMarkers.pdf'), width = 9, height = 7, useDingbats = F)


iCard_TFs_nUpVsExp_wLOESS_wMarkers_w3P <- ggplot() +
  geom_point(data = indiv_perturbability %>% 
               inner_join(tfTab, by = "gene_id") %>%
               inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
               filter(cellType == "iCard") %>% 
               anti_join(markerTab, by = "gene_id") %>%
               filter(!(GeneSymbol %in% c('ZFPM2', 'ESRRG'))) %>%
               mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
             aes(log(meanRPM), nDESeq2conditionsUp)) +
  geom_smooth(data = indiv_perturbability %>% 
                inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard") %>% 
                mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
              aes(log(meanRPM), nDESeq2conditionsUp),
              method = 'loess', color = 'black') +
  geom_point(data = indiv_perturbability %>% 
               inner_join(tfTab, by = "gene_id") %>%
               inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
               filter(cellType == "iCard", gene_id %in% markerTab$gene_id) %>%
               mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
             aes(log(meanRPM), nDESeq2conditionsUp), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% 
                    inner_join(tfTab, by = "gene_id") %>%
                    inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
                    filter(cellType == "iCard", gene_id %in% markerTab$gene_id) %>%
                    mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
                  aes(log(meanRPM), nDESeq2conditionsUp, label = GeneSymbol), color = "red", nudge_y = 0.1, size = 6) +
  geom_point(data = indiv_perturbability %>% 
               inner_join(tfTab, by = "gene_id") %>%
               inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
               filter(cellType == "iCard", GeneSymbol %in% c('SP3', 'ZBTB10', 'ZBTB44')) %>%
               mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
             aes(log(meanRPM), nDESeq2conditionsUp), color = "blue") +
  geom_text_repel(data = indiv_perturbability %>% 
                    inner_join(tfTab, by = "gene_id") %>%
                    inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
                    filter(cellType == "iCard", GeneSymbol %in% c('SP3', 'ZBTB10', 'ZBTB44')) %>%
                    mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
                  aes(log(meanRPM), nDESeq2conditionsUp, label = GeneSymbol), color = "blue", nudge_y = 0.1, size = 6) +
  xlim(c(0,8)) +
  theme_bw() +
  theme(axis.text = element_text(size = rel(2)),
        axis.title = element_blank())
# plot(iCard_TFs_nUpVsExp_wLOESS_wMarkers_w3P)
ggsave(iCard_TFs_nUpVsExp_wLOESS_wMarkers_w3P, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nUpVsExp_wLOESS_wMarkers_w3P.pdf'), width = 9, height = 7, useDingbats = F)



iCard_TFs_nDownVsExp_wLOESS_wMarkers <- ggplot() +
  geom_point(data = indiv_perturbability %>% 
               inner_join(tfTab, by = "gene_id") %>%
               inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
               filter(cellType == "iCard") %>% 
               anti_join(markerTab, by = "gene_id") %>%
               mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
             aes(log(meanRPM), nDESeq2conditionsDown)) +
  geom_smooth(data = indiv_perturbability %>% 
                inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard") %>% 
                mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
              aes(log(meanRPM), nDESeq2conditionsDown),
              method = 'loess', color = 'black') +
  geom_point(data = indiv_perturbability %>% 
               inner_join(tfTab, by = "gene_id") %>%
               inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
               filter(cellType == "iCard", gene_id %in% markerTab$gene_id) %>%
               mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
             aes(log(meanRPM), nDESeq2conditionsDown), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% 
                    inner_join(tfTab, by = "gene_id") %>%
                    inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
                    filter(cellType == "iCard", gene_id %in% markerTab$gene_id) %>%
                    mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
                  aes(log(meanRPM), nDESeq2conditionsDown, label = GeneSymbol), color = "red", nudge_y = -0.1, size = 6) +
  xlim(c(0,8)) +
  theme_bw() +
  theme(axis.text = element_text(size = rel(2)),
        axis.title = element_blank())
# plot(iCard_TFs_nDownVsExp_wLOESS_wMarkers)
ggsave(iCard_TFs_nDownVsExp_wLOESS_wMarkers, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nDownVsExp_wLOESS_wMarkers.pdf'), width = 9, height = 7, useDingbats = F)


# with bootstraps
iCard_boot_up <- read.csv(paste0(graphDir, '/differentialExpression/allSamples/iCards/bootstrap_DESeq2_up.csv'), header = T, stringsAsFactors = F) 
rownames(iCard_boot_up) = iCard_boot_up$X
iCard_boot_up$X <-NULL
iCard_boot_down <- read.csv(paste0(graphDir, '/differentialExpression/allSamples/iCards/bootstrap_DESeq2_down.csv'), header = T, stringsAsFactors = F)
iCard_boot_down$X <-NULL
iCard_boot_noch <- read.csv(paste0(graphDir, '/differentialExpression/allSamples/iCards/bootstrap_DESeq2_noch.csv'), header = T, stringsAsFactors = F)
iCard_boot_noch$X <-NULL


iCard_boot_sds <- as_tibble(data.frame( 
  gene_id = rownames(iCard_boot_up),
  nDEupSD = rowSds(as.matrix(iCard_boot_up)),
  nDEdownSD = rowSds(as.matrix(iCard_boot_down)),
  nDEallSD = rowSds(as.matrix(iCard_boot_down)+as.matrix(iCard_boot_up)),
  nDEnochSD = rowSds(as.matrix(iCard_boot_noch))))

iCard_TFs_nUpVsExp_wLOESS_wMarkers_wMarkerBoot <- ggplot() +
  geom_point(data = indiv_perturbability %>% 
               inner_join(tfTab, by = "gene_id") %>%
               inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
               filter(cellType == "iCard") %>% 
               anti_join(markerTab, by = "gene_id") %>%
               mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
             aes(log(meanRPM), nDESeq2conditionsUp)) +
  # geom_errorbar(data = indiv_perturbability %>% 
  #                 inner_join(iCard_boot_sds, by = 'gene_id') %>%
  #                 inner_join(tfTab, by = "gene_id") %>%
  #                 inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
  #                 filter(cellType == "iCard") %>% 
  #                 anti_join(markerTab, by = "gene_id") %>%
  #                 mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
  #               aes(x = log(meanRPM), ymin = nDESeq2conditionsUp - nDEupSD, ymax = nDESeq2conditionsUp + nDEupSD)) +
  geom_smooth(data = indiv_perturbability %>% 
                inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard") %>% 
                mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
              aes(log(meanRPM), nDESeq2conditionsUp),
              method = 'loess', color = 'black') +
  geom_point(data = indiv_perturbability %>% 
               inner_join(tfTab, by = "gene_id") %>%
               inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
               filter(cellType == "iCard", gene_id %in% markerTab$gene_id) %>%
               mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
             aes(log(meanRPM), nDESeq2conditionsUp), color = "red") +
  geom_errorbar(data = indiv_perturbability %>% 
                  inner_join(iCard_boot_sds, by = 'gene_id') %>%
                  inner_join(tfTab, by = "gene_id") %>%
                  inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
                  filter(cellType == "iCard", gene_id %in% markerTab$gene_id) %>%
                  mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
                aes(x = log(meanRPM), ymin = nDESeq2conditionsUp - nDEupSD, ymax = nDESeq2conditionsUp + nDEupSD), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% 
                    inner_join(tfTab, by = "gene_id") %>%
                    inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
                    filter(cellType == "iCard", gene_id %in% markerTab$gene_id) %>%
                    mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
                  aes(log(meanRPM), nDESeq2conditionsUp, label = GeneSymbol), color = "red", nudge_y = -0.1, size = 6) +
  xlim(c(0,8)) +
  theme_bw() +
  theme(axis.text = element_text(size = rel(2)),
        axis.title = element_blank())
# plot(iCard_TFs_nUpVsExp_wLOESS_wMarkers_wMarkerBoot)
ggsave(iCard_TFs_nUpVsExp_wLOESS_wMarkers_wMarkerBoot, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nUpVsExp_wLOESS_wMarkers_wMarkerBoot.pdf'), width = 9, height = 7, useDingbats = F)

iCard_TFs_nUpVsExp_wLOESS_wMarkers_wMarkerBoot_wLabels <- ggplot() +
  geom_point(data = indiv_perturbability %>% 
               inner_join(tfTab, by = "gene_id") %>%
               inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
               filter(cellType == "iCard") %>% 
               anti_join(markerTab, by = "gene_id") %>%
               mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
             aes(log(meanRPM), nDESeq2conditionsUp)) +
  # geom_errorbar(data = indiv_perturbability %>% 
  #                 inner_join(iCard_boot_sds, by = 'gene_id') %>%
  #                 inner_join(tfTab, by = "gene_id") %>%
  #                 inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
  #                 filter(cellType == "iCard") %>% 
  #                 anti_join(markerTab, by = "gene_id") %>%
  #                 mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
  #               aes(x = log(meanRPM), ymin = nDESeq2conditionsUp - nDEupSD, ymax = nDESeq2conditionsUp + nDEupSD)) +
  geom_smooth(data = indiv_perturbability %>% 
                inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard") %>% 
                mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
              aes(log(meanRPM), nDESeq2conditionsUp),
              method = 'loess', color = 'black') +
  geom_point(data = indiv_perturbability %>% 
               inner_join(tfTab, by = "gene_id") %>%
               inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
               filter(cellType == "iCard", gene_id %in% markerTab$gene_id) %>%
               mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
             aes(log(meanRPM), nDESeq2conditionsUp), color = "red") +
  geom_errorbar(data = indiv_perturbability %>% 
                  inner_join(iCard_boot_sds, by = 'gene_id') %>%
                  inner_join(tfTab, by = "gene_id") %>%
                  inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
                  filter(cellType == "iCard", gene_id %in% markerTab$gene_id) %>%
                  mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
                aes(x = log(meanRPM), ymin = nDESeq2conditionsUp - nDEupSD, ymax = nDESeq2conditionsUp + nDEupSD), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% 
                    inner_join(tfTab, by = "gene_id") %>%
                    inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
                    filter(cellType == "iCard", gene_id %in% markerTab$gene_id) %>%
                    mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
                  aes(log(meanRPM), nDESeq2conditionsUp, label = GeneSymbol), color = "red", nudge_y = -0.1, size = 6) +
  xlim(c(0,8)) +
  theme_bw() +
  theme(axis.text = element_text(size = rel(2)))
# plot(iCard_TFs_nUpVsExp_wLOESS_wMarkers_wMarkerBoot_wLabels)
ggsave(iCard_TFs_nUpVsExp_wLOESS_wMarkers_wMarkerBoot_wLabels, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nUpVsExp_wLOESS_wMarkers_wMarkerBoot_wLabels.pdf'), width = 9, height = 7, useDingbats = F)


iCard_TFs_nDownVsExp_wLOESS_wMarkers_wMarkerBoot <- ggplot() +
  geom_point(data = indiv_perturbability %>% 
               inner_join(tfTab, by = "gene_id") %>%
               inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
               filter(cellType == "iCard") %>% 
               anti_join(markerTab, by = "gene_id") %>%
               mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
             aes(log(meanRPM), nDESeq2conditionsDown)) +
  # geom_errorbar(data = indiv_perturbability %>% 
  #                 inner_join(iCard_boot_sds, by = 'gene_id') %>%
  #                 inner_join(tfTab, by = "gene_id") %>%
  #                 inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
  #                 filter(cellType == "iCard") %>% 
  #                 anti_join(markerTab, by = "gene_id") %>%
  #                 mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
  #               aes(x = log(meanRPM), ymin = nDESeq2conditionsDown - nDEdownSD, ymax = nDESeq2conditionsDown + nDEdownSD)) +
  geom_smooth(data = indiv_perturbability %>% 
                inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard") %>% 
                mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
              aes(log(meanRPM), nDESeq2conditionsDown),
              method = 'loess', color = 'black') +
  geom_point(data = indiv_perturbability %>% 
               inner_join(tfTab, by = "gene_id") %>%
               inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
               filter(cellType == "iCard", gene_id %in% markerTab$gene_id) %>%
               mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
             aes(log(meanRPM), nDESeq2conditionsDown), color = "red") +
  geom_errorbar(data = indiv_perturbability %>% 
                  inner_join(iCard_boot_sds, by = 'gene_id') %>%
                  inner_join(tfTab, by = "gene_id") %>%
                  inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
                  filter(cellType == "iCard", gene_id %in% markerTab$gene_id) %>%
                  mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
                aes(x = log(meanRPM), ymin = nDESeq2conditionsDown - nDEdownSD, ymax = nDESeq2conditionsDown + nDEdownSD), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% 
                    inner_join(tfTab, by = "gene_id") %>%
                    inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
                    filter(cellType == "iCard", gene_id %in% markerTab$gene_id) %>%
                    mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
                  aes(log(meanRPM), nDESeq2conditionsDown, label = GeneSymbol), color = "red", nudge_y = -0.1, size = 6) +
  xlim(c(0,8)) +
  theme_bw() +
  theme(axis.text = element_text(size = rel(2)),
        axis.title = element_blank())
plot(iCard_TFs_nDownVsExp_wLOESS_wMarkers_wMarkerBoot)
ggsave(iCard_TFs_nDownVsExp_wLOESS_wMarkers_wMarkerBoot, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nDownVsExp_wLOESS_wMarkers_wMarkerBoot.pdf'), width = 9, height = 7, useDingbats = F)


iCard_TFs_nUpOrDownVsExp_wLOESS_wMarkers_wMarkerBoot <- ggplot() +
  geom_point(data = indiv_perturbability %>% 
               inner_join(tfTab, by = "gene_id") %>%
               inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
               filter(cellType == "iCard") %>% 
               anti_join(markerTab, by = "gene_id") %>%
               mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
             aes(log(meanRPM), nDESeq2conditionsAll)) +
  # geom_errorbar(data = indiv_perturbability %>% 
  #                 inner_join(iCard_boot_sds, by = 'gene_id') %>%
  #                 inner_join(tfTab, by = "gene_id") %>%
  #                 inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
  #                 filter(cellType == "iCard") %>% 
  #                 anti_join(markerTab, by = "gene_id") %>%
  #                 mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
  #               aes(x = log(meanRPM), ymin = nDESeq2conditionsDown - nDEdownSD, ymax = nDESeq2conditionsDown + nDEdownSD)) +
  geom_smooth(data = indiv_perturbability %>% 
                inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard") %>% 
                mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
              aes(log(meanRPM), nDESeq2conditionsAll),
              method = 'loess', color = 'black') +
  geom_point(data = indiv_perturbability %>% 
               inner_join(tfTab, by = "gene_id") %>%
               inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
               filter(cellType == "iCard", gene_id %in% markerTab$gene_id) %>%
               mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
             aes(log(meanRPM), nDESeq2conditionsAll), color = "red") +
  geom_errorbar(data = indiv_perturbability %>% 
                  inner_join(iCard_boot_sds, by = 'gene_id') %>%
                  inner_join(tfTab, by = "gene_id") %>%
                  inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
                  filter(cellType == "iCard", gene_id %in% markerTab$gene_id) %>%
                  mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
                aes(x = log(meanRPM), ymin = nDESeq2conditionsAll - nDEallSD, ymax = nDESeq2conditionsAll + nDEallSD), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% 
                    inner_join(tfTab, by = "gene_id") %>%
                    inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
                    filter(cellType == "iCard", gene_id %in% markerTab$gene_id) %>%
                    mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
                  aes(log(meanRPM), nDESeq2conditionsAll, label = GeneSymbol), color = "red", nudge_y = -0.1, size = 6) +
  xlim(c(0,8)) +
  theme_bw() +
  theme(axis.text = element_text(size = rel(2)),
        axis.title = element_blank())
# plot(iCard_TFs_nUpOrDownVsExp_wLOESS_wMarkers_wMarkerBoot)
ggsave(iCard_TFs_nUpOrDownVsExp_wLOESS_wMarkers_wMarkerBoot, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nUpOrDownVsExp_wLOESS_wMarkers_wMarkerBoot.pdf'), width = 9, height = 7, useDingbats = F)

iCard_TFs_nUpOrDownVsExp_wLOESS_wMarkers_wMarkerBoot_wLabels <- ggplot() +
  geom_point(data = indiv_perturbability %>% 
               inner_join(tfTab, by = "gene_id") %>%
               inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
               filter(cellType == "iCard") %>% 
               anti_join(markerTab, by = "gene_id") %>%
               mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
             aes(log(meanRPM), nDESeq2conditionsAll)) +
  # geom_errorbar(data = indiv_perturbability %>% 
  #                 inner_join(iCard_boot_sds, by = 'gene_id') %>%
  #                 inner_join(tfTab, by = "gene_id") %>%
  #                 inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
  #                 filter(cellType == "iCard") %>% 
  #                 anti_join(markerTab, by = "gene_id") %>%
  #                 mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
  #               aes(x = log(meanRPM), ymin = nDESeq2conditionsDown - nDEdownSD, ymax = nDESeq2conditionsDown + nDEdownSD)) +
  geom_smooth(data = indiv_perturbability %>% 
                inner_join(tfTab, by = "gene_id") %>%
                filter(cellType == "iCard") %>% 
                mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
              aes(log(meanRPM), nDESeq2conditionsAll),
              method = 'loess', color = 'black') +
  geom_point(data = indiv_perturbability %>% 
               inner_join(tfTab, by = "gene_id") %>%
               inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
               filter(cellType == "iCard", gene_id %in% markerTab$gene_id) %>%
               mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
             aes(log(meanRPM), nDESeq2conditionsAll), color = "red") +
  geom_errorbar(data = indiv_perturbability %>% 
                  inner_join(iCard_boot_sds, by = 'gene_id') %>%
                  inner_join(tfTab, by = "gene_id") %>%
                  inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
                  filter(cellType == "iCard", gene_id %in% markerTab$gene_id) %>%
                  mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
                aes(x = log(meanRPM), ymin = nDESeq2conditionsAll - nDEallSD, ymax = nDESeq2conditionsAll + nDEallSD), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% 
                    inner_join(tfTab, by = "gene_id") %>%
                    inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
                    filter(cellType == "iCard", gene_id %in% markerTab$gene_id) %>%
                    mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
                  aes(log(meanRPM), nDESeq2conditionsAll, label = GeneSymbol), color = "red", nudge_y = -0.1, size = 6) +
  xlim(c(0,8)) +
  theme_bw() +
  theme(axis.text = element_text(size = rel(2)))
# plot(iCard_TFs_nUpOrDownVsExp_wLOESS_wMarkers_wMarkerBoot_wLabels)
ggsave(iCard_TFs_nUpOrDownVsExp_wLOESS_wMarkers_wMarkerBoot_wLabels, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nUpOrDownVsExp_wLOESS_wMarkers_wMarkerBoot_wLabels.pdf'), width = 9, height = 7, useDingbats = F)



# 1D and L1 and avg expression histograms
iCard_TFs_nUp_hist_wMarkers <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20), aes(nDESeq2conditionsUp),
                 binwidth = 1) +
  geom_vline(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(xintercept = nDESeq2conditionsUp), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(x = nDESeq2conditionsUp, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() + 
  xlim(c(-0.5, 12.5)) +
  xlab("Number of conditions UP") +
  ggtitle('Number of drugs in which UP\nAll TFs with mean RPM >20 in iCard controls\niCard marker genes in red')
ggsave(iCard_TFs_nUp_hist_wMarkers, file =  paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_TFs_nUp_hist_wMarkers.pdf'), width = 9, height = 4, useDingbats = F)

iCard_TFs_nDown_hist_wMarkers <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20), aes(nDESeq2conditionsDown),
                 binwidth = 1) +
  geom_vline(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(xintercept = nDESeq2conditionsDown), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(x = nDESeq2conditionsDown, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() + 
  xlim(c(-0.5, 12.5)) +
  xlab("Number of conditions DOWN")+
  ggtitle('Number of drugs in which DOWN\nAll TFs with mean RPM >20 in iCard controls\niCard marker genes in red')
ggsave(iCard_TFs_nDown_hist_wMarkers, file =  paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_TFs_nDown_hist_wMarkers.pdf'), width = 9, height = 4, useDingbats = F)


iCard_TFs_nUpAndDown_hist_wMarkers <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20), aes(nDESeq2conditionsAll),
                 binwidth = 1) +
  geom_vline(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(xintercept = nDESeq2conditionsAll), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(x = nDESeq2conditionsAll, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() + 
  xlim(c(-0.5, 12.5)) +
  xlab("Number of conditions UP+DOWN") +
  ggtitle('Number of drugs in which UP OR DOWN\nAll TFs with mean RPM >20 in iCard controls\niCard marker genes in red')
ggsave(iCard_TFs_nUpAndDown_hist_wMarkers, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nUpAndDown_hist_wMarkers.pdf'), width = 9, height = 4, useDingbats = F)

iCard_TFs_meanTPM_hist_wMarkers <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20), aes(meanTPM)) +
  geom_vline(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(xintercept = meanTPM), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(x = meanTPM, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() + 
  # xlim(c(-0.5, 12.5)) +
  xlab("mean TPM in controls")
ggsave(iCard_TFs_meanTPM_hist_wMarkers, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_meanTPM_hist_wMarkers.pdf'), width = 9, height = 4, useDingbats = F)

iCard_TFs_logMeanTPM_hist_wMarkers <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20), aes(log(meanTPM))) +
  geom_vline(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(xintercept = log(meanTPM)), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique(), aes(x = log(meanTPM), y = 200, label = GeneSymbol), color = "red") +
  theme_bw() + 
  # xlim(c(-0.5, 12.5)) +
  xlab("log(meanTPM) in controls") +
  ggtitle('Log(meanTPM) of all TFs with mean RPM >20 in iCard controls\niCard marker genes in red')
ggsave(iCard_TFs_logMeanTPM_hist_wMarkers, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_logMeanTPM_hist_wMarkers.pdf'), width = 9, height = 4, useDingbats = F)

# add jitter
set.seed(26054)

iCard_markers_updown_jitter_min20rpm <- indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
  filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique() %>%
  dplyr::select(gene_id, meanRPM, meanTPM, GeneSymbol, nDESeq2conditionsAll, nDESeq2conditionsUp, nDESeq2conditionsDown)
iCard_markers_updown_jitter_min20rpm$nDESeq2conditionsAll_jit = iCard_markers_updown_jitter_min20rpm$nDESeq2conditionsAll + runif(nrow(iCard_markers_updown_jitter_min20rpm), min = -0.2, max = 0.2)
iCard_markers_updown_jitter_min20rpm$nDESeq2conditionsUp_jit = iCard_markers_updown_jitter_min20rpm$nDESeq2conditionsUp + runif(nrow(iCard_markers_updown_jitter_min20rpm), min = -0.2, max = 0.2)
iCard_markers_updown_jitter_min20rpm$nDESeq2conditionsDown_jit = iCard_markers_updown_jitter_min20rpm$nDESeq2conditionsDown + runif(nrow(iCard_markers_updown_jitter_min20rpm), min = -0.2, max = 0.2)

set.seed(79406)
iCard_otherTFs_updown_jitter_min20rpm <- indiv_perturbability %>% anti_join(markerTab2, by = "gene_id") %>%
  filter(cellType == "iCard", meanRPM > 20) %>% unique() %>%
  dplyr::select(gene_id, meanRPM, meanTPM, gene_id, nDESeq2conditionsAll, nDESeq2conditionsUp, nDESeq2conditionsDown)
iCard_otherTFs_updown_jitter_min20rpm$nDESeq2conditionsAll_jit = iCard_otherTFs_updown_jitter_min20rpm$nDESeq2conditionsAll + runif(nrow(iCard_otherTFs_updown_jitter_min20rpm), min = -0.4, max = 0.4)
iCard_otherTFs_updown_jitter_min20rpm$nDESeq2conditionsUp_jit = iCard_otherTFs_updown_jitter_min20rpm$nDESeq2conditionsUp + runif(nrow(iCard_otherTFs_updown_jitter_min20rpm), min = -0.4, max = 0.4)
iCard_otherTFs_updown_jitter_min20rpm$nDESeq2conditionsDown_jit = iCard_otherTFs_updown_jitter_min20rpm$nDESeq2conditionsDown + runif(nrow(iCard_otherTFs_updown_jitter_min20rpm), min = -0.4, max = 0.4)


set.seed(26053)
fibro_markers_updown_jitter_min20rpm <- indiv_perturbability %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
  filter(cellType == "GM00942", meanRPM > 20, GeneSymbol %in% targs16) %>% unique() %>% mutate(isBarrier = GeneSymbol %in% targs8) %>%
  dplyr::select(gene_id, meanRPM, meanTPM, GeneSymbol, isBarrier, nDESeq2conditionsAll, nDESeq2conditionsUp, nDESeq2conditionsDown)
fibro_markers_updown_jitter_min20rpm$nDESeq2conditionsAll_jit = fibro_markers_updown_jitter_min20rpm$nDESeq2conditionsAll + runif(nrow(fibro_markers_updown_jitter_min20rpm), min = -0.2, max = 0.2)
fibro_markers_updown_jitter_min20rpm$nDESeq2conditionsUp_jit = fibro_markers_updown_jitter_min20rpm$nDESeq2conditionsUp + runif(nrow(fibro_markers_updown_jitter_min20rpm), min = -0.2, max = 0.2)
fibro_markers_updown_jitter_min20rpm$nDESeq2conditionsDown_jit = fibro_markers_updown_jitter_min20rpm$nDESeq2conditionsDown + runif(nrow(fibro_markers_updown_jitter_min20rpm), min = -0.2, max = 0.2)

set.seed(79405)
fibro_otherTFs_updown_jitter_min20rpm <- indiv_perturbability %>% anti_join(markerTab2, by = "gene_id") %>%
  filter(cellType == "GM00942", meanRPM > 20) %>% unique() %>%
  dplyr::select(gene_id, meanRPM, meanTPM, nDESeq2conditionsAll, nDESeq2conditionsUp, nDESeq2conditionsDown)
fibro_otherTFs_updown_jitter_min20rpm$nDESeq2conditionsAll_jit = fibro_otherTFs_updown_jitter_min20rpm$nDESeq2conditionsAll + runif(nrow(fibro_otherTFs_updown_jitter_min20rpm), min = -0.4, max = 0.4)
fibro_otherTFs_updown_jitter_min20rpm$nDESeq2conditionsUp_jit = fibro_otherTFs_updown_jitter_min20rpm$nDESeq2conditionsUp + runif(nrow(fibro_otherTFs_updown_jitter_min20rpm), min = -0.4, max = 0.4)
fibro_otherTFs_updown_jitter_min20rpm$nDESeq2conditionsDown_jit = fibro_otherTFs_updown_jitter_min20rpm$nDESeq2conditionsDown + runif(nrow(fibro_otherTFs_updown_jitter_min20rpm), min = -0.4, max = 0.4)


iCard_TFs_nUp_hist_wMarkers_wJitter <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20), aes(nDESeq2conditionsUp),
                 binwidth = 1) +
  geom_vline(data = iCard_markers_updown_jitter_min20rpm, aes(xintercept = nDESeq2conditionsUp_jit), color = "red") +
  geom_text_repel(data = iCard_markers_updown_jitter_min20rpm, aes(x = nDESeq2conditionsUp_jit, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() + 
  xlim(c(-0.5, 12.5)) +
  xlab("Number of conditions UP") +
  ggtitle('Number of drugs in which UP\nAll TFs with mean RPM >20 in iCard controls\niCard marker genes in red')
ggsave(iCard_TFs_nUp_hist_wMarkers_wJitter, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nUp_hist_wMarkers_wJitter.pdf'), width = 9, height = 4, useDingbats = F)

iCard_TFs_nDown_hist_wMarkers_wJitter <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20), aes(nDESeq2conditionsDown),
                 binwidth = 1) +
  geom_vline(data = iCard_markers_updown_jitter_min20rpm, aes(xintercept = nDESeq2conditionsDown_jit), color = "red") +
  geom_text_repel(data = iCard_markers_updown_jitter_min20rpm, aes(x = nDESeq2conditionsDown_jit, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() + 
  xlim(c(-0.5, 12.5)) +
  xlab("Number of conditions DOWN")+
  ggtitle('Number of drugs in which DOWN\nAll TFs with mean RPM >20 in iCard controls\niCard marker genes in red')
ggsave(iCard_TFs_nDown_hist_wMarkers_wJitter, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nDown_hist_wMarkers_wJitter.pdf'), width = 9, height = 4, useDingbats = F)


iCard_TFs_nUpAndDown_hist_wMarkers_wJitter <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20), aes(nDESeq2conditionsAll),
                 binwidth = 1) +
  geom_vline(data = iCard_markers_updown_jitter_min20rpm, aes(xintercept = nDESeq2conditionsAll_jit), color = "red") +
  geom_text_repel(data = iCard_markers_updown_jitter_min20rpm, aes(x = nDESeq2conditionsAll_jit, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() + 
  xlim(c(-0.5, 12.5)) +
  xlab("Number of conditions UP+DOWN") +
  ggtitle('Number of drugs in which UP OR DOWN\nAll TFs with mean RPM >20 in iCard controls\niCard marker genes in red')
ggsave(iCard_TFs_nUpAndDown_hist_wMarkers_wJitter, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nUpAndDown_hist_wMarkers_wJitter.pdf'), width = 9, height = 4, useDingbats = F)

# also add rugs
set.seed(26054)
iCard_TFs_nUp_hist_wMarkers_wRugs_wJitter <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20) %>% mutate(graphType = 'histogram'), aes(nDESeq2conditionsUp),
                 binwidth = 1) +
  geom_vline(data = iCard_markers_updown_jitter_min20rpm %>% mutate(graphType = 'histogram'), aes(xintercept = nDESeq2conditionsUp_jit), color = "red") +
  geom_text_repel(data = iCard_markers_updown_jitter_min20rpm %>% mutate(graphType = 'histogram'), aes(x = nDESeq2conditionsUp_jit, y = 200, label = GeneSymbol), color = "red") +
  geom_linerange(data = iCard_otherTFs_updown_jitter_min20rpm %>% inner_join(tfTab) %>% mutate(graphType = 'rug'), aes(nDESeq2conditionsUp_jit, ymin = 0, ymax = 1), alpha = 0.2, color = 'black') +
  geom_linerange(data = iCard_markers_updown_jitter_min20rpm %>% mutate(graphType = 'rug'), aes(nDESeq2conditionsUp_jit, ymin = 1, ymax = 2), color = "red") +
  geom_text_repel(data = iCard_markers_updown_jitter_min20rpm %>% mutate(graphType = 'rug'), aes(x = nDESeq2conditionsUp_jit, y = 1.75, label = GeneSymbol), color = "red") +
  theme_bw() + 
  facet_grid(graphType ~ ., scales = 'free_y') +
  xlim(c(-0.5, 12.5)) +
  xlab("Number of conditions UP") +
  ggtitle('Number of drugs in which UP\nAll TFs with mean RPM >20 in iCard controls\niCard marker genes in red')
ggsave(iCard_TFs_nUp_hist_wMarkers_wRugs_wJitter, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nUp_hist_wMarkers_wRugs_wJitter.pdf'), width = 9, height = 4, useDingbats = F)

iCard_TFs_nDown_hist_wMarkers_wRugs_wJitter <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20) %>% mutate(graphType = 'histogram'), aes(nDESeq2conditionsDown),
                 binwidth = 1) +
  geom_vline(data = iCard_markers_updown_jitter_min20rpm %>% mutate(graphType = 'histogram'), aes(xintercept = nDESeq2conditionsDown_jit), color = "red") +
  geom_text_repel(data = iCard_markers_updown_jitter_min20rpm %>% mutate(graphType = 'histogram'), aes(x = nDESeq2conditionsDown_jit, y = 200, label = GeneSymbol), color = "red") +
  geom_linerange(data = iCard_otherTFs_updown_jitter_min20rpm %>% inner_join(tfTab) %>% mutate(graphType = 'rug'), aes(nDESeq2conditionsDown_jit, ymin = 0, ymax = 1), alpha = 0.2, color = 'black') +
  geom_linerange(data = iCard_markers_updown_jitter_min20rpm %>% mutate(graphType = 'rug'), aes(nDESeq2conditionsDown_jit, ymin = 1, ymax = 2), color = "red") +
  geom_text_repel(data = iCard_markers_updown_jitter_min20rpm %>% mutate(graphType = 'rug'), aes(x = nDESeq2conditionsDown_jit, y = 1.75, label = GeneSymbol), color = "red") +
  theme_bw() + 
  facet_grid(graphType ~ ., scales = 'free_y') +
  xlim(c(-0.5, 12.5)) +
  xlab("Number of conditions DOWN")+
  ggtitle('Number of drugs in which DOWN\nAll TFs with mean RPM >20 in iCard controls\niCard marker genes in red')
ggsave(iCard_TFs_nDown_hist_wMarkers_wRugs_wJitter, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nDown_hist_wMarkers_wRugs_wJitter.pdf'), width = 9, height = 4, useDingbats = F)


iCard_TFs_nUpAndDown_hist_wMarkers_wRugs_wJitter <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20) %>% mutate(graphType = 'histogram'), aes(nDESeq2conditionsAll),
                 binwidth = 1) +
  geom_vline(data = iCard_markers_updown_jitter_min20rpm %>% mutate(graphType = 'histogram'), aes(xintercept = nDESeq2conditionsAll_jit), color = "red") +
  geom_text_repel(data = iCard_markers_updown_jitter_min20rpm %>% mutate(graphType = 'histogram'), aes(x = nDESeq2conditionsAll_jit, y = 200, label = GeneSymbol), color = "red") +
  geom_linerange(data = iCard_otherTFs_updown_jitter_min20rpm %>% inner_join(tfTab) %>% mutate(graphType = 'rug'), aes(nDESeq2conditionsAll_jit, ymin = 0, ymax = 1), alpha = 0.2, color = 'black') +
  geom_linerange(data = iCard_markers_updown_jitter_min20rpm %>% mutate(graphType = 'rug'), aes(nDESeq2conditionsAll_jit, ymin = 1, ymax = 2), color = "red") +
  geom_text_repel(data = iCard_markers_updown_jitter_min20rpm %>% mutate(graphType = 'rug'), aes(x = nDESeq2conditionsAll_jit, y = 1.75, label = GeneSymbol), color = "red") +
  facet_grid(graphType ~ ., scales = 'free_y') +
  theme_bw() + 
  xlim(c(-0.5, 12.5)) +
  xlab("Number of conditions UP+DOWN") +
  ggtitle('Number of drugs in which UP OR DOWN\nAll TFs with mean RPM >20 in iCard controls\niCard marker genes in red')
ggsave(iCard_TFs_nUpAndDown_hist_wMarkers_wRugs_wJitter, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nUpAndDown_hist_wMarkers_wRugs_wJitter.pdf'), width = 9, height = 4, useDingbats = F)

iCard_TFs_logMeanTPM_hist_wMarkers_wRugs <- ggplot() +
  geom_histogram(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20) %>% mutate(graphType = 'histogram'), aes(log(meanTPM))) +
  geom_vline(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique() %>% mutate(graphType = 'histogram'), aes(xintercept = log(meanTPM)), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique() %>% mutate(graphType = 'histogram'), aes(x = log(meanTPM), y = 200, label = GeneSymbol), color = "red") +
  geom_linerange(data = indiv_perturbability %>% inner_join(tfTab, by = "gene_id") %>% filter(cellType == "iCard", meanRPM > 20) %>% mutate(graphType = 'rug'), aes(log(meanTPM), ymin = 0, ymax = 1)) +
  geom_linerange(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
               filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique() %>% mutate(graphType = 'rug'), aes(log(meanTPM), ymin = 1, ymax = 2), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% inner_join(markerTab2, by = "gene_id") %>%
                    filter(cellType == "iCard", cellTypeMark == "iCard", meanRPM > 20) %>% unique() %>% mutate(graphType = 'rug'), aes(x = log(meanTPM), y = 1.75, label = GeneSymbol), color = "red") +
  theme_bw() + 
  facet_grid(graphType ~ ., scales = 'free_y') +
  # xlim(c(-0.5, 12.5)) +
  xlab("log(meanTPM) in controls") +
  ggtitle('Log(meanTPM) of all TFs with mean RPM >20 in iCard controls\niCard marker genes in red')
ggsave(iCard_TFs_logMeanTPM_hist_wMarkers_wRugs, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_logMeanTPM_hist_wMarkers_wRugs.pdf'), width = 9, height = 4, useDingbats = F)


# histograms faceted
indiv_perturbability[is.na(indiv_perturbability)] <- 0
pert_for_hists <- indiv_perturbability %>% 
  inner_join(tfTab, by = "gene_id") %>% 
  filter(cellType == "iCard", meanRPM > 20) %>% 
  dplyr::select(gene_id, nDESeq2conditionsAll, nDESeq2conditionsUp) %>%
  gather('Measure', 'nDrugs', 2:3)

markers_for_hists <- iCard_markers_updown_jitter_min20rpm %>%
  dplyr::select(GeneSymbol, nDESeq2conditionsAll_jit, nDESeq2conditionsUp_jit) %>%
  dplyr::rename(nDESeq2conditionsAll = nDESeq2conditionsAll_jit,
                nDESeq2conditionsUp = nDESeq2conditionsUp_jit) %>%
  gather('Measure', 'nDrugs', 2:3)


# mean_spec_for_hists <- indiv_perturbability %>% 
#   inner_join(tfTab, by = "gene_id") %>% 
#   filter(cellType == "iCard", meanRPM > 20) %>% 
#   dplyr::select(gene_id, meanTPM) %>%
#   gather('Measure', 'nDrugs', 2:3)

set.seed(26054)
iCard_TFs_nAll_nUp_hists_wMarkers_wJitter <- ggplot() +
  geom_histogram(data = pert_for_hists, aes(nDrugs),
                 binwidth = 1) +
  geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() + 
  facet_grid(. ~ Measure, scales = 'free_x') +
  xlab("Perturbability measure") +
  ggtitle('Number of drugs in which perturbed (any) or perturbed-up\nAll TFs with mean RPM >20 in iCard controls\niCard marker genes in red')
ggsave(iCard_TFs_nAll_nUp_hists_wMarkers_wJitter, file = paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nAll_nUp_hists_wMarkers_wJitter.pdf'), width = 3.5, height = 1.5, useDingbats = F)

set.seed(26054)
iCard_TFs_nAll_nUp_hists_wMarkers_wJitter_noLabs <- ggplot() +
  geom_histogram(data = pert_for_hists, aes(nDrugs),
                 binwidth = 1) +
  geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() + 
  facet_grid(. ~ Measure, scales = 'free_x') +
  theme(axis.text = element_blank(),
        axis.title = element_blank())
ggsave(iCard_TFs_nAll_nUp_hists_wMarkers_wJitter_noLabs, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nAll_nUp_hists_wMarkers_wJitter_noLabs.pdf'), width = 3.5, height = 1.5, useDingbats = F)

set.seed(26054)
iCard_TFs_nAll_nUp_hists_wMarkerHists_wJitter_noLabs <- ggplot() +
  geom_histogram(data = pert_for_hists, aes(nDrugs),
                 binwidth = 1) +
  geom_histogram(data = markers_for_hists, aes(nDrugs), fill = 'red',
                 binwidth = 1) +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() + 
  facet_grid(. ~ Measure, scales = 'free_x') +
  theme(axis.text = element_blank(),
        axis.title = element_blank())
ggsave(iCard_TFs_nAll_nUp_hists_wMarkerHists_wJitter_noLabs, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nAll_nUp_hists_wMarkerHists_wJitter_noLabs.pdf'), width = 3.5, height = 1.5, useDingbats = F)

set.seed(26054)
iCard_TFs_nAll_nUp_hists_wMarkerDens_wJitter_noLabs <- ggplot() +
  geom_density(data = pert_for_hists, aes(nDrugs), fill = 'grey', alpha = 0.5) +
  geom_density(data = markers_for_hists, aes(nDrugs), fill = 'red', alpha = 0.5) +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_bw() + 
  facet_grid(. ~ Measure, scales = 'free_x') +
  theme(axis.text = element_blank(),
        axis.title = element_blank())
ggsave(iCard_TFs_nAll_nUp_hists_wMarkerDens_wJitter_noLabs, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nAll_nUp_hists_wMarkerDens_wJitter_noLabs.pdf'), width = 3.5, height = 1.5, useDingbats = F)

set.seed(26054)
iCard_TFs_nAll_nUp_hists_wMarkerLinesegs_wJitter_noLabs <- ggplot() +
  geom_histogram(data = pert_for_hists, aes(nDrugs), fill = 'grey', alpha = 0.9,
                 binwidth = 1) +
  geom_linerange(data = markers_for_hists, aes(nDrugs, ymin = 0, ymax = 25), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 0, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  # ylim(c(-100, 200)) +
  facet_grid(. ~ Measure, scales = 'free_x') +
  theme(axis.text = element_blank(),
        axis.title = element_blank())
ggsave(iCard_TFs_nAll_nUp_hists_wMarkerLinesegs_wJitter_noLabs, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nAll_nUp_hists_wMarkerLinesegs_wJitter_noLabs.pdf'), width = 4, height = 1, useDingbats = F)

set.seed(26054)
iCard_TFs_nAll_nUp_hists_wMarkerLinesegs_wJitter_noLabs_thin <- ggplot() +
  geom_histogram(data = pert_for_hists, aes(nDrugs), fill = 'grey', alpha = 0.9,
                 binwidth = 1) +
  geom_linerange(data = markers_for_hists, aes(nDrugs, ymin = 0, ymax = 40), size = 0.2, color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 0, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  # ylim(c(-100, 200)) +
  facet_grid(. ~ Measure, scales = 'free_x') +
  theme(axis.text = element_blank(),
        axis.title = element_blank())
ggsave(iCard_TFs_nAll_nUp_hists_wMarkerLinesegs_wJitter_noLabs_thin, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nAll_nUp_hists_wMarkerLinesegs_wJitter_noLabs_thin.pdf'), width = 4, height = 1, useDingbats = F)

set.seed(26054)
iCard_TFs_nAll_nUp_hists_w7FLinesegs_wJitter_noLabs_thin <- ggplot() +
  geom_histogram(data = pert_for_hists, aes(nDrugs), fill = 'grey', alpha = 0.9,
                 binwidth = 1) +
  geom_linerange(data = markers_for_hists %>% filter(GeneSymbol %in% factors_7F), aes(nDrugs, ymin = 0, ymax = 40), size = 0.2, color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 0, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  # ylim(c(-100, 200)) +
  facet_grid(. ~ Measure, scales = 'free_x') +
  theme(axis.text = element_blank(),
        axis.title = element_blank())
ggsave(iCard_TFs_nAll_nUp_hists_w7FLinesegs_wJitter_noLabs_thin, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nAll_nUp_hists_w7FLinesegs_wJitter_noLabs_thin.pdf'), width = 4, height = 1, useDingbats = F)

# fracUp hist and scatter
pert_for_hists2 <- indiv_perturbability %>% 
  inner_join(tfTab, by = "gene_id") %>% 
  filter(cellType == "iCard", meanRPM > 20) %>% 
  dplyr::select(gene_id, nDESeq2conditionsAll, nDESeq2conditionsUp) %>%
  mutate(fracConditionsUp = nDESeq2conditionsUp/nDESeq2conditionsAll) %>%
  gather('Measure', 'nDrugs', 2:4)

markers_for_hists2 <- iCard_markers_updown_jitter_min20rpm %>%
  dplyr::select(GeneSymbol, nDESeq2conditionsAll_jit, nDESeq2conditionsUp_jit) %>%
  dplyr::rename(nDESeq2conditionsAll = nDESeq2conditionsAll_jit,
                nDESeq2conditionsUp = nDESeq2conditionsUp_jit) %>%
  mutate(fracConditionsUp = nDESeq2conditionsUp/nDESeq2conditionsAll) %>%
  gather('Measure', 'nDrugs', 2:4)

set.seed(26054)
iCard_TFs_nAll_nUp_fracUp_hists_wMarkerLinesegs_wJitter_noLabs_thin <- ggplot() +
  geom_histogram(data = pert_for_hists2, aes(nDrugs), fill = 'grey', alpha = 0.9,
                 binwidth = 1) +
  geom_linerange(data = markers_for_hists2, aes(nDrugs, ymin = 0, ymax = 40), size = 0.2, color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 0, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  # ylim(c(-100, 200)) +
  facet_grid(. ~ Measure, scales = 'free_x') +
  theme(axis.text = element_blank(),
        axis.title = element_blank())
ggsave(iCard_TFs_nAll_nUp_fracUp_hists_wMarkerLinesegs_wJitter_noLabs_thin, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nAll_nUp_fracUp_hists_wMarkerLinesegs_wJitter_noLabs_thin.pdf'), width = 4, height = 1, useDingbats = F)


fibro_pert_for_hists <- indiv_perturbability %>% 
  inner_join(tfTab, by = "gene_id") %>% 
  filter(cellType == "GM00942", meanRPM > 20) %>% 
  dplyr::select(gene_id, nDESeq2conditionsAll, nDESeq2conditionsUp) %>%
  gather('Measure', 'nDrugs', 2:3)

fibro_markers_for_hists <- fibro_markers_updown_jitter_min20rpm %>%
  dplyr::select(GeneSymbol, isBarrier, nDESeq2conditionsAll_jit, nDESeq2conditionsUp_jit) %>%
  dplyr::rename(nDESeq2conditionsAll = nDESeq2conditionsAll_jit,
                nDESeq2conditionsUp = nDESeq2conditionsUp_jit) %>%
  gather('Measure', 'nDrugs', 3:4)
set.seed(26054)
fibro_TFs_nAll_nUp_hists_wBarrierLinesegs_wJitter_noLabs_thin <- ggplot() +
  geom_histogram(data = fibro_pert_for_hists, aes(nDrugs), fill = 'grey', alpha = 0.9,
                 binwidth = 1) +
  geom_linerange(data = fibro_markers_for_hists %>% filter(isBarrier == T), aes(nDrugs, ymin = 0, ymax = 40), size = 0.2, color = "green") +
  # geom_linerange(data = fibro_markers_for_hists %>% filter(isBarrier == F), aes(nDrugs, ymin = 0, ymax = 40), size = 0.2, color = "darkolivegreen2") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 0, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  # ylim(c(-100, 200)) +
  facet_grid(. ~ Measure, scales = 'free_x') +
  theme(axis.text = element_blank(),
        axis.title = element_blank())
ggsave(fibro_TFs_nAll_nUp_hists_wBarrierLinesegs_wJitter_noLabs_thin, file =  paste0(graphDir,'/differentialExpression/allSamples/fibroblasts/fibro_TFs_nAll_nUp_hists_wBarrierLinesegs_wJitter_noLabs_thin.pdf'), width = 4, height = 1, useDingbats = F)

fibro_TFs_nAll_nUp_hists_wBarrierLinesegs_wJitter_thin <- ggplot() +
  geom_histogram(data = fibro_pert_for_hists, aes(nDrugs), fill = 'grey', alpha = 0.9,
                 binwidth = 1) +
  geom_linerange(data = fibro_markers_for_hists %>% filter(isBarrier == T), aes(nDrugs, ymin = 0, ymax = 40), size = 0.2, color = "green") +
  # geom_linerange(data = fibro_markers_for_hists %>% filter(isBarrier == F), aes(nDrugs, ymin = 0, ymax = 40), size = 0.2, color = "darkolivegreen1") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 0, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  # ylim(c(-100, 200)) +
  facet_grid(. ~ Measure, scales = 'free_x') 

ggsave(fibro_TFs_nAll_nUp_hists_wBarrierLinesegs_wJitter_thin, file =  paste0(graphDir,'/differentialExpression/allSamples/fibroblasts/fibro_TFs_nAll_nUp_hists_wBarrierLinesegs_wJitter_thin.pdf'), width = 12, height = 3, useDingbats = F)

# wide histograms
iCard_avg_pert_for_wides_expressedTFs <- indiv_perturbability %>%
  filter(cellType == 'iCard', meanRPM > 20) %>%
  inner_join(tfTab, by = 'gene_id') %>%
  dplyr::select(gene_id, meanTPM, nDESeq2conditionsAll, nDESeq2conditionsUp, nDESeq2conditionsDown)
  
iCard_avg_pert_for_wides_markers <- iCard_markers_updown_jitter_min20rpm %>%
  inner_join(markerTab2, by = c('gene_id', 'GeneSymbol')) %>%
  filter(cellTypeMark == 'iCard') %>% unique() %>%
  dplyr::select(gene_id, GeneSymbol, meanTPM, nDESeq2conditionsAll_jit, nDESeq2conditionsUp_jit, nDESeq2conditionsDown_jit) %>%
  dplyr::rename(nDESeq2conditionsAll = nDESeq2conditionsAll_jit,
                nDESeq2conditionsUp = nDESeq2conditionsUp_jit,
                nDESeq2conditionsDown = nDESeq2conditionsDown_jit)

iCard_TFs_nAll_hist_wMarkerLinesegs_wJitter_noLabs_thin_wide <- ggplot() +
  geom_histogram(data = iCard_avg_pert_for_wides_expressedTFs, aes(nDESeq2conditionsAll), fill = 'grey', alpha = 0.9,
                 binwidth = 1) +
  geom_linerange(data = iCard_avg_pert_for_wides_markers, aes(nDESeq2conditionsAll, ymin = 0, ymax = 40), size = 0.2, color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 0, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  # ylim(c(-100, 200)) +
  # facet_grid(. ~ Measure, scales = 'free_x') +
  theme(axis.text = element_blank(),
        axis.title = element_blank())
ggsave(iCard_TFs_nAll_hist_wMarkerLinesegs_wJitter_noLabs_thin_wide, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nAll_hist_wMarkerLinesegs_wJitter_noLabs_thin_wide.pdf'), width = 3.2, height = 0.68, useDingbats = F)

iCard_TFs_nAll_hist_wMarkerLinesegs_wJitter_wLabs_thin_wide <- ggplot() +
  geom_histogram(data = iCard_avg_pert_for_wides_expressedTFs, aes(nDESeq2conditionsAll), fill = 'grey', alpha = 0.9,
                 binwidth = 1) +
  geom_linerange(data = iCard_avg_pert_for_wides_markers, aes(nDESeq2conditionsAll, ymin = 0, ymax = 40), size = 0.2, color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 0, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() #+ 
  # ylim(c(-100, 200)) +
  # facet_grid(. ~ Measure, scales = 'free_x') +
  # theme(axis.text = element_blank(),
  #       axis.title = element_blank())
ggsave(iCard_TFs_nAll_hist_wMarkerLinesegs_wJitter_wLabs_thin_wide, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nAll_hist_wMarkerLinesegs_wJitter_wLabs_thin_wide.pdf'), width = 3.2, height = 0.68, useDingbats = F)


iCard_TFs_nUp_hist_wMarkerLinesegs_wJitter_wLabs_thin_wide <- ggplot() +
  geom_histogram(data = iCard_avg_pert_for_wides_expressedTFs, aes(nDESeq2conditionsUp), fill = 'grey', alpha = 0.9,
                 binwidth = 1) +
  geom_linerange(data = iCard_avg_pert_for_wides_markers, aes(nDESeq2conditionsUp, ymin = 0, ymax = 40), size = 0.2, color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 0, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  # ylim(c(-100, 200)) +
  # facet_grid(. ~ Measure, scales = 'free_x') +
  theme(axis.text = element_blank(),
        axis.title = element_blank())
ggsave(iCard_TFs_nUp_hist_wMarkerLinesegs_wJitter_wLabs_thin_wide, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nUp_hist_wMarkerLinesegs_wJitter_wLabs_thin_wide.pdf'), width = 3.2, height = 0.68, useDingbats = F)

iCard_TFs_nDown_hist_wMarkerLinesegs_wJitter_noLabs_thin_wide <- ggplot() +
  geom_histogram(data = iCard_avg_pert_for_wides_expressedTFs, aes(nDESeq2conditionsDown), fill = 'grey', alpha = 0.9,
                 binwidth = 1) +
  geom_linerange(data = iCard_avg_pert_for_wides_markers, aes(nDESeq2conditionsDown, ymin = 0, ymax = 40), size = 0.2, color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 0, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  # ylim(c(-100, 200)) +
  # facet_grid(. ~ Measure, scales = 'free_x') +
  theme(axis.text = element_blank(),
        axis.title = element_blank())
ggsave(iCard_TFs_nDown_hist_wMarkerLinesegs_wJitter_noLabs_thin_wide, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nDown_hist_wMarkerLinesegs_wJitter_noLabs_thin_wide.pdf'), width = 3.2, height = 0.68, useDingbats = F)

iCard_TFs_nDown_hist_wMarkerLinesegs_wJitter_wLabs_thin_wide <- ggplot() +
  geom_histogram(data = iCard_avg_pert_for_wides_expressedTFs, aes(nDESeq2conditionsDown), fill = 'grey', alpha = 0.9,
                 binwidth = 1) +
  geom_linerange(data = iCard_avg_pert_for_wides_markers, aes(nDESeq2conditionsDown, ymin = 0, ymax = 40), size = 0.2, color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 0, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() #+ 
  # ylim(c(-100, 200)) +
  # facet_grid(. ~ Measure, scales = 'free_x') +
  # theme(axis.text = element_blank(),
  #       axis.title = element_blank())
ggsave(iCard_TFs_nDown_hist_wMarkerLinesegs_wJitter_wLabs_thin_wide, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_nDown_hist_wMarkerLinesegs_wJitter_wLabs_thin_wide.pdf'), width = 3.2, height = 1.68, useDingbats = F)


iCard_TFs_meanTPM_hist_wMarkerLinesegs_wJitter_noLabs_thin_wide <- ggplot() +
  geom_histogram(data = iCard_avg_pert_for_wides_expressedTFs, aes(log(meanTPM)), fill = 'grey', alpha = 0.9) +
  geom_linerange(data = iCard_avg_pert_for_wides_markers, aes(log(meanTPM), ymin = 0, ymax = 30), size = 0.2, color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 0, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  xlim(c(0,8)) +
  # ylim(c(-100, 200)) +
  # facet_grid(. ~ Measure, scales = 'free_x') +
  theme(axis.text = element_blank(),
        axis.title = element_blank())
ggsave(iCard_TFs_meanTPM_hist_wMarkerLinesegs_wJitter_noLabs_thin_wide, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_meanTPM_hist_wMarkerLinesegs_wJitter_noLabs_thin_wide.pdf'), width = 3.2, height = 0.68, useDingbats = F)

iCard_TFs_meanTPM_hist_wMarkerLinesegs_wJitter_wLabs_thin_wide <- ggplot() +
  geom_histogram(data = iCard_avg_pert_for_wides_expressedTFs, aes(log(meanTPM)), fill = 'grey', alpha = 0.9) +
  geom_linerange(data = iCard_avg_pert_for_wides_markers, aes(log(meanTPM), ymin = 0, ymax = 30), size = 0.2, color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 0, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  xlim(c(0,8)) #+
  # ylim(c(-100, 200)) +
  # facet_grid(. ~ Measure, scales = 'free_x') +
  # theme(axis.text = element_blank(),
  #       axis.title = element_blank())
ggsave(iCard_TFs_meanTPM_hist_wMarkerLinesegs_wJitter_wLabs_thin_wide, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_meanTPM_hist_wMarkerLinesegs_wJitter_wLabs_thin_wide.pdf'), width = 3.2, height = 1.68, useDingbats = F)

iCard_TFs_meanTPM_hist_wMarkerLinesegs_wGeneNames_wJitter_wLabs_thin_wide <- ggplot() +
  geom_histogram(data = iCard_avg_pert_for_wides_expressedTFs, aes(log(meanTPM)), fill = 'grey', alpha = 0.9) +
  geom_linerange(data = iCard_avg_pert_for_wides_markers, aes(log(meanTPM), ymin = 0, ymax = 30), size = 0.2, color = "red") +
  geom_text_repel(data = iCard_avg_pert_for_wides_markers, aes(x = log(meanTPM), y = 30, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  xlim(c(0,8)) #+
# ylim(c(-100, 200)) +
# facet_grid(. ~ Measure, scales = 'free_x') +
# theme(axis.text = element_blank(),
#       axis.title = element_blank())
ggsave(iCard_TFs_meanTPM_hist_wMarkerLinesegs_wGeneNames_wJitter_wLabs_thin_wide, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_TFs_meanTPM_hist_wMarkerLinesegs_wGeneNames_wJitter_wLabs_thin_wide.pdf'), width = 3.2, height = 1.68, useDingbats = F)


## Effect sizes for marker genes
iCard_meanRPMs <- indiv_perturbability %>%
  filter(cellType == 'iCard') %>%
  dplyr::select(gene_id, meanRPM)
iCard_DE_marker_effectsizes <- ggplot() +
  geom_boxplot(data = iCard_DE_results_all %>% filter(padj < 0.1) %>% inner_join(markerTab2, by = 'gene_id') %>% filter(cellTypeMark == 'iCard') %>% unique() %>% inner_join(iCard_meanRPMs, by = 'gene_id') %>% filter(meanRPM > 20),
               aes(GeneSymbol, abs(log2FoldChange))) +
  geom_point(data = iCard_DE_results_all %>% filter(padj < 0.1) %>% inner_join(markerTab2, by = 'gene_id') %>% filter(cellTypeMark == 'iCard') %>% unique() %>% inner_join(iCard_meanRPMs, by = 'gene_id') %>% filter(meanRPM > 20),
             aes(GeneSymbol, abs(log2FoldChange))) +
  geom_hline(aes(yintercept = 0.5), linetype = 'dashed') +
  ylim(c(0,2.5)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, color = 'red'))
ggsave(iCard_DE_marker_effectsizes, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_DE_marker_effectsizes.pdf'), width = 4, height = 4, useDingbats = F)

iCard_DE_sum_effectSize_sigOnly <- iCard_DE_results_all %>% 
  filter(padj < 0.1) %>% 
  inner_join(iCard_meanRPMs, by = 'gene_id') %>%
  group_by(gene_id, meanRPM) %>%
  summarise(meanAbsEffSize = mean(abs(log2FoldChange)))

set.seed(827)
iCard_DE_allTF_effectSizevsmeanRPM_sigOnly <- ggplot() +
  geom_point(data = iCard_DE_sum_effectSize_sigOnly %>% inner_join(tfTab) %>% filter(meanRPM > 20),
             aes(log(meanRPM), meanAbsEffSize)) +
  geom_text_repel(data = iCard_DE_sum_effectSize_sigOnly %>% inner_join(markerTab2, by = 'gene_id') %>% filter(cellTypeMark == 'iCard') %>% unique() %>% filter(meanRPM > 20),
               aes(log(meanRPM), meanAbsEffSize, label = GeneSymbol), color = 'red') +
  geom_point(data = iCard_DE_sum_effectSize_sigOnly %>% inner_join(markerTab2, by = 'gene_id') %>% filter(cellTypeMark == 'iCard') %>% unique() %>% filter(meanRPM > 20),
             aes(log(meanRPM), meanAbsEffSize), color = 'red') +
  # geom_hline(aes(yintercept = 0.5), linetype = 'dashed') +
  # ylim(c(0,2.5)) +
  theme_classic() +
  theme(axis.text = element_text(color = 'black'))
ggsave(iCard_DE_allTF_effectSizevsmeanRPM_sigOnly, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_DE_allTF_effectSizevsmeanRPM_sigOnly.pdf'), width = 4, height = 4, useDingbats = F)

iCard_DE_sum_effectSize_allConds <- iCard_DE_results_all %>% 
  inner_join(iCard_meanRPMs, by = 'gene_id') %>%
  group_by(gene_id, meanRPM) %>%
  summarise(meanAbsEffSize = mean(abs(log2FoldChange)))

set.seed(827)
iCard_DE_allTF_effectSizevsmeanRPM_allConds <- ggplot() +
  geom_point(data = iCard_DE_sum_effectSize_allConds %>% inner_join(tfTab) %>% filter(meanRPM > 20),
             aes(log(meanRPM), meanAbsEffSize), alpha = 0.3) +
  geom_text_repel(data = iCard_DE_sum_effectSize_allConds %>% inner_join(markerTab2, by = 'gene_id') %>% filter(cellTypeMark == 'iCard') %>% unique() %>% filter(meanRPM > 20),
                  aes(log(meanRPM), meanAbsEffSize, label = GeneSymbol), color = 'red') +
  geom_point(data = iCard_DE_sum_effectSize_allConds %>% inner_join(markerTab2, by = 'gene_id') %>% filter(cellTypeMark == 'iCard') %>% unique() %>% filter(meanRPM > 20),
             aes(log(meanRPM), meanAbsEffSize), color = 'red') +
  # geom_hline(aes(yintercept = 0.5), linetype = 'dashed') +
  ylim(c(0,1.6)) +
  theme_classic() +
  theme(axis.text = element_text(color = 'black'))
ggsave(iCard_DE_allTF_effectSizevsmeanRPM_allConds, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_DE_allTF_effectSizevsmeanRPM_allConds.pdf'), width = 4, height = 4, useDingbats = F)

iCard_DE_sum_effectSizeSigned_allConds <- iCard_DE_results_all %>% 
  inner_join(iCard_meanRPMs, by = 'gene_id') %>%
  group_by(gene_id, meanRPM) %>%
  summarise(meanEffSize = mean(log2FoldChange))

set.seed(827)
iCard_DE_allTF_signedEffectSizevsmeanRPM_allConds <- ggplot() +
  geom_point(data = iCard_DE_sum_effectSizeSigned_allConds %>% inner_join(tfTab) %>% filter(meanRPM > 20),
             aes(log(meanRPM), meanEffSize), alpha = 0.3) +
  geom_text_repel(data = iCard_DE_sum_effectSizeSigned_allConds %>% inner_join(markerTab2, by = 'gene_id') %>% filter(cellTypeMark == 'iCard') %>% unique() %>% filter(meanRPM > 20),
                  aes(log(meanRPM), meanEffSize, label = GeneSymbol), color = 'red') +
  geom_point(data = iCard_DE_sum_effectSizeSigned_allConds %>% inner_join(markerTab2, by = 'gene_id') %>% filter(cellTypeMark == 'iCard') %>% unique() %>% filter(meanRPM > 20),
             aes(log(meanRPM), meanEffSize), color = 'red') +
  # geom_hline(aes(yintercept = 0.5), linetype = 'dashed') +
  # ylim(c(0,1.6)) +
  theme_classic() +
  theme(axis.text = element_text(color = 'black'))
ggsave(iCard_DE_allTF_signedEffectSizevsmeanRPM_allConds, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_DE_allTF_signedEffectSizevsmeanRPM_allConds.pdf'), width = 4, height = 4, useDingbats = F)



## Frac up vs nUp
set.seed(96438)
iCard_TFs_nAllVsfracUp_noLOESS_vertJit_wMarkers <- ggplot() +
  geom_jitter(data = indiv_perturbability %>% 
                inner_join(tfTab, by = "gene_id") %>%
                inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
                filter(cellType == "iCard") %>% 
                anti_join(markerTab2, by = "gene_id") %>%
                filter(!(GeneSymbol %in% c('ZFPM2', 'ESRRG')), meanRPM > 20) %>%
                mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
              aes(nDESeq2conditionsAll, fracUp), width = 0.2, height = 0.05, alpha = 0.5) +
  geom_point(data = indiv_perturbability %>% 
               # inner_join(tfTab, by = "gene_id") %>%
               inner_join(markerTab2, by = 'gene_id') %>%
               filter(cellType == "iCard", cellTypeMark == 'iCard', meanRPM > 20) %>%
               mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
             aes(nDESeq2conditionsAll, fracUp), color = "red") +
  geom_text_repel(data = indiv_perturbability %>% 
                    # inner_join(tfTab, by = "gene_id") %>%
                    inner_join(markerTab2, by = 'gene_id') %>%
                    filter(cellType == "iCard", cellTypeMark == 'iCard', meanRPM > 20) %>% unique() %>%
                    mutate(fracUp = ifelse(nDESeq2conditionsAll > 0, nDESeq2conditionsUp/nDESeq2conditionsAll, 0)), 
                  aes(nDESeq2conditionsAll, fracUp, label = GeneSymbol), color = "red", nudge_y = 0.05, size = 6) +
  # xlim(c(0,8)) +
  theme_bw() +
  theme(axis.text = element_text(size = rel(2)),
        axis.title = element_blank())
# plot(iCard_TFs_nUpVsExp_noLOESS_vertJit_wMarkers)
setwd(projectDir)
ggsave(iCard_TFs_nAllVsfracUp_noLOESS_vertJit_wMarkers, file = paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_TFs_nAllVsfracUp_noLOESS_vertJit_wMarkers.pdf'), width = 9, height = 7, useDingbats = F)





## zhouOlson perturbability ####
# zhouOlson_indiv_perturbability_RPMnorm20 <- indiv_perturbability_RPMnorm20 %>%
#   filter(cellType == "iCard") %>%
#   inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
#   inner_join(zhouOlson_allORFs, by = "GeneSymbol") %>%
#   mutate(effectOnXdiff = ifelse(GeneSymbol %in% zhouOlson_activators$GeneSymbol, "activator",
#                                 ifelse(GeneSymbol %in% zhouOlson_inhibitors$GeneSymbol, "inhibitor", "no_effect")))
# zhouOlson_indiv_perturbability_tfOnlyRPMnorm20 <- indiv_perturbability_tfOnlyRPMnorm20 %>%
#   filter(cellType == "iCard") %>%
#   inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
#   inner_join(zhouOlson_allORFs, by = "GeneSymbol") %>%
#   mutate(effectOnXdiff = ifelse(GeneSymbol %in% zhouOlson_activators$GeneSymbol, "activator",
#                                 ifelse(GeneSymbol %in% zhouOlson_inhibitors$GeneSymbol, "inhibitor", "no_effect")))
# 
# 
# 
# zhouOlson_tfOnly_normDeltaKurtosis <- ggplot()+
#   geom_density2d(data = zhouOlson_indiv_perturbability_tfOnlyRPMnorm20, aes(log(meanTPM), normDeltaKurtosis)) +
#   theme_bw() +
#   ggtitle("perturbability of Zhou, Olson ORFs\niCard data, only TF genes")
# ggsave(zhouOlson_tfOnly_normDeltaKurtosis, file = paste0(graphDir, "/perturbability/zhouOlson_tfOnly_normDeltaKurtosis.pdf"), width = 8, height = 8)
# 
# zhouOlson_tfOnly_normDeltaKurtosis_withInhibitors <- ggplot()+
#   geom_density2d(data = zhouOlson_indiv_perturbability_tfOnlyRPMnorm20, aes(log(meanTPM), normDeltaKurtosis)) +
#   geom_point(data = zhouOlson_indiv_perturbability_tfOnlyRPMnorm20 %>% inner_join(zhouOlson_inhibitors, by = "GeneSymbol"), 
#              aes(log(meanTPM), normDeltaKurtosis), color = "darkolivegreen") +
#   geom_text_repel(data = zhouOlson_indiv_perturbability_tfOnlyRPMnorm20 %>% inner_join(zhouOlson_inhibitors, by = "GeneSymbol"), 
#                   aes(log(meanTPM), normDeltaKurtosis, label = GeneSymbol), color = "darkolivegreen") +
#   theme_bw() +
#   ggtitle("perturbability of Zhou, Olson ORFs\niCard data, only TF genes")
# ggsave(zhouOlson_tfOnly_normDeltaKurtosis_withInhibitors, file = paste0(graphDir, "/perturbability/zhouOlson_tfOnly_normDeltaKurtosis_withInhibitors.pdf"), width = 8, height = 8, useDingbats = F)
# 
# zhouOlson_tfOnly_normDeltaKurtosis_withActivatorsAndInhibitors <- ggplot()+
#   geom_density2d(data = zhouOlson_indiv_perturbability_tfOnlyRPMnorm20, aes(log(meanTPM), normDeltaKurtosis)) +
#   geom_point(data = zhouOlson_indiv_perturbability_tfOnlyRPMnorm20 %>% inner_join(zhouOlson_activators, by = "GeneSymbol"), 
#              aes(log(meanTPM), normDeltaKurtosis), color = "red") +
#   geom_text_repel(data = zhouOlson_indiv_perturbability_tfOnlyRPMnorm20 %>% inner_join(zhouOlson_activators, by = "GeneSymbol"), 
#                   aes(log(meanTPM), normDeltaKurtosis, label = GeneSymbol), color = "red", size = 8) +
#   geom_point(data = zhouOlson_indiv_perturbability_tfOnlyRPMnorm20 %>% inner_join(zhouOlson_inhibitors, by = "GeneSymbol"), 
#              aes(log(meanTPM), normDeltaKurtosis), color = "darkolivegreen") +
#     theme_bw() +
#   ggtitle("perturbability of Zhou, Olson ORFs\niCard data, only TF genes")
# ggsave(zhouOlson_tfOnly_normDeltaKurtosis_withActivatorsAndInhibitors, file = paste0(graphDir, "/perturbability/zhouOlson_tfOnly_normDeltaKurtosis_withActivatorsAndInhibitors.pdf"), width = 8, height = 8, useDingbats = F)
# 
# 
# zhouOlson_boxplots_normDeltaKurtosis = ggplot(zhouOlson_indiv_perturbability_RPMnorm20, aes(effectOnXdiff, normDeltaKurtosis)) +
#   geom_boxplot() +
#   ylab("Normalized change in kurtosis")
# plot(zhouOlson_boxplots_normDeltaKurtosis)
# 
# zhouOlson_boxplots_normDeltaSkewness = ggplot(zhouOlson_indiv_perturbability_RPMnorm20, aes(effectOnXdiff, normDeltaSkewness)) +
#   geom_boxplot() +
#   ylab("Normalized change in skewness")
# plot(zhouOlson_boxplots_normDeltaSkewness)
# 
# zhouOlson_tfOnly_boxplots_normDeltaKurtosis = ggplot(zhouOlson_indiv_perturbability_tfOnlyRPMnorm20, aes(effectOnXdiff, normDeltaKurtosis)) +
#   geom_boxplot() +
#   ylab("Normalized change in kurtosis")
# plot(zhouOlson_tfOnly_boxplots_normDeltaKurtosis)
# 
# zhouOlson_tfOnly_boxplots_normDeltaSkewness = ggplot(zhouOlson_indiv_perturbability_tfOnlyRPMnorm20, aes(effectOnXdiff, normDeltaSkewness)) +
#   geom_boxplot() +
#   ylab("Normalized change in skewness")
# plot(zhouOlson_tfOnly_boxplots_normDeltaSkewness)

