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

# graphing network perturbability

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(readxl)
library(magrittr)
library(e1071)
library(tibble)

# # if running manually in Rstudio, start here and use:
# projectDir = '~/Dropbox (RajLab)/Projects/cellid/'
# procDataSubdir = 'procDataScripted'
# graphSubdir = 'graphs'

procDataDir = paste0(projectDir, procDataSubdir)
graphDir = paste0(projectDir, graphSubdir)

setwd(projectDir)
cat("Starting graphing_regulons.R...\n")
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

## load Marbach 2016 heart (general) network
heartNet = as_tibble(read.table("extractedData/Networks/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/32_high-level_networks/21_heart.txt", sep = "\t", header = F, stringsAsFactors = F))
colnames(heartNet) = c("TF", "Target", "edgeWt")
heartTFs = unique(heartNet$TF)

## load Marbach 2016 integumental (general) network - inlcudes all skin fibroblast nets
integNet = as_tibble(read.table(gzfile("extractedData/Networks/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/32_high-level_networks/09_connective_tissue_integumental_cells.txt.gz"), sep = "\t", header = F, stringsAsFactors = F))
colnames(integNet) = c("TF", "Target", "edgeWt")
integTFs = unique(integNet$TF)

## load Marbach 2016 "normal skin fibroblast" individual network
nlSkinNet = as_tibble(read.table(gzfile("extractedData/Networks/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/394_individual_networks/fibroblast_-_skin_normal.txt.gz"), sep = "\t", header = F, stringsAsFactors = F))
colnames(nlSkinNet) = c("TF", "Target", "edgeWt")
nlSkinTFs = unique(nlSkinNet$TF)

## load Marbach 2016 "dermal fibroblast" nets
dermNet = as_tibble(read.table(gzfile("extractedData/Networks/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/394_individual_networks/fibroblast_-_dermal.txt.gz"), sep = "\t", header = F, stringsAsFactors = F))
colnames(dermNet) = c("TF", "Target", "edgeWt")
dermTFs = unique(dermNet$TF)




## Load perturbability data
# for normalized values, working with default windowRadius = 20
indiv_perturbability = as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_variability_metrics.txt"), header = T, stringsAsFactors = F, sep = "\t")) %>%
  inner_join(geneIDtoGeneSymbol, by = "gene_id")
indiv_perturbability_RPMnorm20 = as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_tpmFilt_variability_RPMnormMetrics_window20.txt"), header = T, stringsAsFactors = F, sep = "\t")) %>%
  inner_join(geneIDtoGeneSymbol, by = "gene_id")
indiv_perturbability[is.na(indiv_perturbability)] = 0
indiv_perturbability_RPMnorm20[is.na(indiv_perturbability_RPMnorm20)] = 0

## fraction of expressed targets in heart dysregulated
iCard_heartNet_regulonStats = list()
iCard_heartNet_regulonStats_minRPM20 = list()
n = 1
for (tf in heartTFs) {
  
  tfMeanTPM = indiv_perturbability %>%
    filter(cellType == "iCard", GeneSymbol == tf) %>%
    dplyr::select(GeneSymbol, meanTPM, meanRPM) %>%
    dplyr::rename(TF = GeneSymbol,
                  tf.meanTPM = meanTPM,
                  tf.meanRPM = meanRPM)
  
  tempNet = heartNet %>%
    filter(TF == tf) %>%
    mutate(GeneSymbol = Target) %>%
    dplyr::select(GeneSymbol, edgeWt)
  
  temp_heartNet_regulonStats = indiv_perturbability %>%
    filter(cellType == "iCard", meanTPM > 10) %>%
    inner_join(tempNet, by = "GeneSymbol") %>%
    summarise(edgeWt_NDEup = sum(nDESeq2conditionsUp * edgeWt, na.rm = T)/sum(edgeWt),
              edgeWt_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_NDEup = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM), na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM) * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt)) %>%
    mutate(TF = tf) %>%
    left_join(tfMeanTPM, by = "TF")
  
  temp_heartNet_regulonStats_minRPM20 = indiv_perturbability %>%
    filter(cellType == "iCard", meanRPM > 20) %>%
    inner_join(tempNet, by = "GeneSymbol") %>%
    summarise(edgeWt_NDEup = sum(nDESeq2conditionsUp * edgeWt, na.rm = T)/sum(edgeWt),
              edgeWt_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_NDEup = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM), na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM) * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt)) %>%
    mutate(TF = tf) %>%
    left_join(tfMeanTPM, by = "TF")
  
  if(is.null(dim(iCard_heartNet_regulonStats))){
    iCard_heartNet_regulonStats = temp_heartNet_regulonStats
    iCard_heartNet_regulonStats_minRPM20 = temp_heartNet_regulonStats_minRPM20
  } else {
    iCard_heartNet_regulonStats = bind_rows(iCard_heartNet_regulonStats, temp_heartNet_regulonStats)
    iCard_heartNet_regulonStats_minRPM20 = bind_rows(iCard_heartNet_regulonStats_minRPM20, temp_heartNet_regulonStats_minRPM20)
  }
  
  if(n %% 25 == 0){cat(paste0("Completed ", as.character(n), "/", as.character(length(heartTFs)), " TF target sets...\n"))}
  
  n = n+1
  
}

iCard_heartNet_regulonStats[is.na(iCard_heartNet_regulonStats)]<-0
iCard_heartNet_regulonStats_minRPM20[is.na(iCard_heartNet_regulonStats_minRPM20)]<-0

iCard_heartNet_regulonStats_tall <- iCard_heartNet_regulonStats %>%
  group_by(TF, tf.meanTPM) %>%
  gather("regulonStat", "value", 1:(ncol(iCard_heartNet_regulonStats) -3))

iCard_heartNet_regulonStats_tall_minRPM20 <- iCard_heartNet_regulonStats_minRPM20 %>%
  group_by(TF, tf.meanTPM) %>%
  gather("regulonStat", "value", 1:(ncol(iCard_heartNet_regulonStats) -3))


iCard_heartNet_regulonStats = ggplot() +
  facet_wrap(~regulonStat, scales = "free") +
  geom_point(data = iCard_heartNet_regulonStats_tall %>% filter(tf.meanRPM > 20), aes(tf.meanTPM, value)) +
  geom_point(data = iCard_heartNet_regulonStats_tall %>% dplyr::rename(GeneSymbol = TF) %>% inner_join(markerTab2, by = "GeneSymbol") %>% filter(cellTypeMark == "iCard", tf.meanTPM > 5) %>% unique(), 
             aes(tf.meanTPM, value), color = "red") +
  geom_text_repel(data = iCard_heartNet_regulonStats_tall %>% dplyr::rename(GeneSymbol = TF) %>% inner_join(markerTab2, by = "GeneSymbol") %>% filter(cellTypeMark == "iCard", tf.meanTPM > 5) %>% unique(), 
                  aes(tf.meanTPM, value, label = GeneSymbol), color = "red") +
  geom_rug(data = iCard_heartNet_regulonStats_tall %>% filter(tf.meanTPM > 20), aes(tf.meanTPM, value),
           sides = "l") +
  geom_rug(data = iCard_heartNet_regulonStats_tall %>% dplyr::rename(GeneSymbol = TF) %>% inner_join(markerTab2, by = "GeneSymbol") %>% filter(cellTypeMark == "iCard", tf.meanTPM > 5) %>% unique(), 
           aes(tf.meanTPM, value), color = "red", sides = "r")
# plot(iCard_heartNet_regulonStats)

iCard_heartNet_regulonStats_minRPM20 = ggplot() +
  facet_wrap(~regulonStat, scales = "free") +
  geom_point(data = iCard_heartNet_regulonStats_tall_minRPM20 %>% filter(tf.meanRPM > 20), aes(tf.meanTPM, value)) +
  geom_point(data = iCard_heartNet_regulonStats_tall_minRPM20 %>% dplyr::rename(GeneSymbol = TF) %>% inner_join(markerTab2, by = "GeneSymbol") %>% filter(cellTypeMark == "iCard", tf.meanTPM > 5) %>% unique(), 
             aes(tf.meanTPM, value), color = "red") +
  geom_text_repel(data = iCard_heartNet_regulonStats_tall_minRPM20 %>% dplyr::rename(GeneSymbol = TF) %>% inner_join(markerTab2, by = "GeneSymbol") %>% filter(cellTypeMark == "iCard", tf.meanTPM > 5) %>% unique(), 
                  aes(tf.meanTPM, value, label = GeneSymbol), color = "red") +
  geom_rug(data = iCard_heartNet_regulonStats_tall_minRPM20 %>% filter(tf.meanRPM > 20), aes(tf.meanTPM, value),
           sides = "l") +
  geom_rug(data = iCard_heartNet_regulonStats_tall_minRPM20 %>% dplyr::rename(GeneSymbol = TF) %>% inner_join(markerTab, by = "GeneSymbol") %>% filter(cellTypeMark == "iCard", tf.meanTPM > 5) %>% unique(), 
           aes(tf.meanTPM, value), color = "red", sides = "r")
# plot(iCard_heartNet_regulonStats_minRPM20)

iCard_heartNet_regulonStats_minRPM20_histograms_noLabs <- ggplot() +
  geom_histogram(data = iCard_heartNet_regulonStats_tall_minRPM20 %>% filter(tf.meanRPM > 20), aes(value), fill = 'grey', alpha = 0.9) +
  geom_linerange(data = iCard_heartNet_regulonStats_tall_minRPM20 %>% dplyr::rename(GeneSymbol = TF) %>% inner_join(markerTab2, by = "GeneSymbol") %>% filter(cellTypeMark == "iCard", tf.meanTPM > 5) %>% unique(), 
                 aes(value, ymin = 0, ymax = 10), size = 0.2, color = "red") +
  geom_text_repel(data = iCard_heartNet_regulonStats_tall_minRPM20 %>% dplyr::rename(GeneSymbol = TF) %>% inner_join(markerTab2, by = "GeneSymbol") %>% filter(cellTypeMark == "iCard", tf.meanTPM > 5) %>% unique(), 
                  aes(x = value, y = 10, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  # ylim(c(-100, 200)) +
  facet_grid(. ~ regulonStat, scales = 'free_x') +
  theme(axis.text = element_blank(),
        axis.title = element_blank())

iCard_heartNet_regulon_nDEupfracUpEdgeWt_minRPM20_histograms_noLabs <- ggplot() +
  geom_histogram(data = iCard_heartNet_regulonStats_tall_minRPM20 %>% filter(tf.meanRPM > 20, regulonStat == 'edgeWt_NDEup_fracUp'), aes(value), fill = 'grey', alpha = 0.9) +
  geom_linerange(data = iCard_heartNet_regulonStats_tall_minRPM20 %>% dplyr::rename(GeneSymbol = TF) %>% inner_join(markerTab2, by = "GeneSymbol") %>% filter(cellTypeMark == "iCard", tf.meanTPM > 5, regulonStat == 'edgeWt_NDEup_fracUp') %>% unique(), 
                 aes(value, ymin = 0, ymax = 10), size = 0.2, color = "red") +
  geom_text_repel(data = iCard_heartNet_regulonStats_tall_minRPM20 %>% dplyr::rename(GeneSymbol = TF) %>% inner_join(markerTab2, by = "GeneSymbol") %>% filter(cellTypeMark == "iCard", tf.meanTPM > 5, regulonStat == 'edgeWt_NDEup_fracUp') %>% unique(), 
                  aes(x = value, y = 10, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  # ylim(c(-100, 200)) +
  facet_grid(. ~ regulonStat, scales = 'free_x') +
  theme(axis.text = element_blank(),
        axis.title = element_blank())
ggsave(iCard_heartNet_regulon_nDEupfracUpEdgeWt_minRPM20_histograms_noLabs, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_heartNet_regulon_nDEupfracUpEdgeWt_minRPM20_histograms_noLabs.pdf'), width = 4, height = 2, useDingbats = F, units = 'in')

iCard_heartNet_regulon_nDEupfracUpEdgeWt_minRPM20_histograms_wLabs <- ggplot() +
  geom_histogram(data = iCard_heartNet_regulonStats_tall_minRPM20 %>% filter(tf.meanRPM > 20, regulonStat == 'edgeWt_NDEup_fracUp'), aes(value), fill = 'grey', alpha = 0.9) +
  geom_linerange(data = iCard_heartNet_regulonStats_tall_minRPM20 %>% dplyr::rename(GeneSymbol = TF) %>% inner_join(markerTab2, by = "GeneSymbol") %>% filter(cellTypeMark == "iCard", tf.meanTPM > 5, regulonStat == 'edgeWt_NDEup_fracUp') %>% unique(), 
                 aes(value, ymin = 0, ymax = 10), size = 0.2, color = "red") +
  geom_text_repel(data = iCard_heartNet_regulonStats_tall_minRPM20 %>% dplyr::rename(GeneSymbol = TF) %>% inner_join(markerTab2, by = "GeneSymbol") %>% filter(cellTypeMark == "iCard", tf.meanTPM > 5, regulonStat == 'edgeWt_NDEup_fracUp') %>% unique(), 
                  aes(x = value, y = 10, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  # ylim(c(-100, 200)) +
  facet_grid(. ~ regulonStat, scales = 'free_x')
ggsave(iCard_heartNet_regulon_nDEupfracUpEdgeWt_minRPM20_histograms_wLabs, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_heartNet_regulon_nDEupfracUpEdgeWt_minRPM20_histograms_wLabs.pdf'), width = 4, height = 2.5, useDingbats = F)


## Scatter vs gene-level
iCtemp_reg <- iCard_heartNet_regulonStats_tall_minRPM20 %>% filter(tf.meanRPM > 20, regulonStat == 'edgeWt_NDEup_fracUp') %>% dplyr::rename(Regulon_responsiveness = value)
iCtemp_reg_gene <- indiv_perturbability %>%
  filter(cellType == 'iCard') %>%
  # inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
  inner_join(iCtemp_reg %>% dplyr::rename(GeneSymbol = TF))

set.seed(37392)
iCard_heartNet_regulon_nDEupfracUpEdgeWt_vsGeneLevelUpreg_scatter <- ggplot() +
  geom_jitter(data = iCtemp_reg_gene %>% anti_join(markerTab2 %>% filter(cellTypeMark == 'iCard'), by = c('gene_id', 'GeneSymbol')), aes(nDESeq2conditionsUp, Regulon_responsiveness), width = 0.2, height = 0.02) +
  geom_point(data = iCtemp_reg_gene %>% inner_join(markerTab2 %>% filter(cellTypeMark == 'iCard'), by = c('gene_id', 'GeneSymbol')), aes(nDESeq2conditionsUp, Regulon_responsiveness), color = 'red') +
  geom_text_repel(data = iCtemp_reg_gene %>% inner_join(markerTab2 %>% filter(cellTypeMark == 'iCard'), by = c('gene_id', 'GeneSymbol')) %>% unique(), aes(nDESeq2conditionsUp, Regulon_responsiveness, label = GeneSymbol), color = 'red') +
  theme_classic()
ggsave(iCard_heartNet_regulon_nDEupfracUpEdgeWt_vsGeneLevelUpreg_scatter, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_heartNet_regulon_nDEupfracUpEdgeWt_vsGeneLevelUpreg_scatter.pdf'), width = 3, height = 2.5, useDingbats = F)


set.seed(37392)
iCard_heartNet_regulon_nDEupfracUpEdgeWt_vsGeneLevelUpreg_scatter <- ggplot() +
  geom_jitter(data = iCtemp_reg_gene %>% anti_join(markerTab2 %>% filter(cellTypeMark == 'iCard'), by = c('gene_id', 'GeneSymbol')), aes(nDESeq2conditionsUp, Regulon_responsiveness), width = 0.2, height = 0.02) +
  geom_point(data = iCtemp_reg_gene %>% inner_join(markerTab2 %>% filter(cellTypeMark == 'iCard'), by = c('gene_id', 'GeneSymbol')), aes(nDESeq2conditionsUp, Regulon_responsiveness), color = 'red') +
  geom_text_repel(data = iCtemp_reg_gene %>% inner_join(markerTab2 %>% filter(cellTypeMark == 'iCard'), by = c('gene_id', 'GeneSymbol')) %>% unique(), aes(nDESeq2conditionsUp, Regulon_responsiveness, label = GeneSymbol), color = 'red') +
  theme_classic()
ggsave(iCard_heartNet_regulon_nDEupfracUpEdgeWt_vsGeneLevelUpreg_scatter, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_heartNet_regulon_nDEupfracUpEdgeWt_vsGeneLevelUpreg_scatter.pdf'), width = 3, height = 2.5, useDingbats = F)


iChNcor <- cor.test(iCtemp_reg_gene$nDESeq2conditionsUp, iCtemp_reg_gene$Regulon_responsiveness)
iChNcor_tbl <- tibble(
  R = iChNcor$estimate,
  p.val = iChNcor$p.value
)
write.table(iChNcor_tbl, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_heartNet_regulon_nUP_pearsoncor.txt'), quote = F, row.names = F, sep = '\t')

### normalized responsiveness

iCtemp_reg <- iCard_heartNet_regulonStats_tall_minRPM20 %>% filter(tf.meanRPM > 20, regulonStat == 'edgeWt_NDEup_fracUp') %>% dplyr::rename(Regulon_responsiveness = value)
iCtemp_reg_gene_RPMnorm20 <- indiv_perturbability_RPMnorm20 %>%
  filter(cellType == 'iCard') %>%
  # inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
  inner_join(iCtemp_reg %>% dplyr::rename(GeneSymbol = TF))

set.seed(37392)
iCard_heartNet_regulon_normnDEupfracUpEdgeWt_vsGeneLevelUpreg_scatter_RPMnorm20 <- ggplot() +
  geom_jitter(data = iCtemp_reg_gene_RPMnorm20 %>% anti_join(markerTab2 %>% filter(cellTypeMark == 'iCard'), by = c('gene_id', 'GeneSymbol')), aes(norm_nDESeq2conditionsUp, Regulon_responsiveness), width = 0.02, height = 0.02) +
  geom_point(data = iCtemp_reg_gene_RPMnorm20 %>% inner_join(markerTab2 %>% filter(cellTypeMark == 'iCard'), by = c('gene_id', 'GeneSymbol')), aes(norm_nDESeq2conditionsUp, Regulon_responsiveness), color = 'red') +
  geom_text_repel(data = iCtemp_reg_gene_RPMnorm20 %>% inner_join(markerTab2 %>% filter(cellTypeMark == 'iCard'), by = c('gene_id', 'GeneSymbol')) %>% unique(), aes(norm_nDESeq2conditionsUp, Regulon_responsiveness, label = GeneSymbol), color = 'red') +
  theme_classic()

# normalized responsiveness of all target genes
iCard_heartNet_regulonStats_targRPMnorm20 = list()
iCard_heartNet_regulonStats_minRPM20_targRPMnorm20 = list()
n = 1
for (tf in heartTFs) {
  
  tfMeanTPM = indiv_perturbability_RPMnorm20 %>%
    filter(cellType == "iCard", GeneSymbol == tf) %>%
    dplyr::select(GeneSymbol, meanTPM, meanRPM) %>%
    dplyr::rename(TF = GeneSymbol,
                  tf.meanTPM = meanTPM,
                  tf.meanRPM = meanRPM)
  
  tempNet = heartNet %>%
    filter(TF == tf) %>%
    mutate(GeneSymbol = Target) %>%
    dplyr::select(GeneSymbol, edgeWt)
  
  temp_heartNet_regulonStats = indiv_perturbability_RPMnorm20 %>%
    filter(cellType == "iCard") %>%
    inner_join(tempNet, by = "GeneSymbol") %>%
    summarise(edgeWt_NDEup = sum(nDESeq2conditionsUp * edgeWt, na.rm = T)/sum(edgeWt),
              edgeWt_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_NDEup = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM), na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM) * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt),
              edgeWt_normNDEup = sum(norm_nDESeq2conditionsUp * edgeWt, na.rm = T)/sum(edgeWt),
              edgeWt_normNDEup_fracUp = sum(norm_nDESeq2conditionsUp * edgeWt * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_normNDEup = sum(norm_nDESeq2conditionsUp * edgeWt * log(meanTPM), na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_normNDEup_fracUp = sum(norm_nDESeq2conditionsUp * edgeWt * log(meanTPM) * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt)) %>%
    mutate(TF = tf) %>%
    left_join(tfMeanTPM, by = "TF")
  
  temp_heartNet_regulonStats_minRPM20 = indiv_perturbability_RPMnorm20 %>%
    filter(cellType == "iCard", meanRPM > 20) %>%
    inner_join(tempNet, by = "GeneSymbol") %>%
    summarise(edgeWt_NDEup = sum(nDESeq2conditionsUp * edgeWt, na.rm = T)/sum(edgeWt),
              edgeWt_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_NDEup = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM), na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM) * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt),
              edgeWt_normNDEup = sum(norm_nDESeq2conditionsUp * edgeWt, na.rm = T)/sum(edgeWt),
              edgeWt_normNDEup_fracUp = sum(norm_nDESeq2conditionsUp * edgeWt * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_normNDEup = sum(norm_nDESeq2conditionsUp * edgeWt * log(meanTPM), na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_normNDEup_fracUp = sum(norm_nDESeq2conditionsUp * edgeWt * log(meanTPM) * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt)) %>%
    mutate(TF = tf) %>%
    left_join(tfMeanTPM, by = "TF")
  
  if(is.null(dim(iCard_heartNet_regulonStats_targRPMnorm20))){
    iCard_heartNet_regulonStats_targRPMnorm20 = temp_heartNet_regulonStats
    iCard_heartNet_regulonStats_minRPM20_targRPMnorm20 = temp_heartNet_regulonStats_minRPM20
  } else {
    iCard_heartNet_regulonStats_targRPMnorm20 = bind_rows(iCard_heartNet_regulonStats_targRPMnorm20, temp_heartNet_regulonStats)
    iCard_heartNet_regulonStats_minRPM20_targRPMnorm20 = bind_rows(iCard_heartNet_regulonStats_minRPM20_targRPMnorm20, temp_heartNet_regulonStats_minRPM20)
  }
  
  if(n %% 25 == 0){cat(paste0("Completed ", as.character(n), "/", as.character(length(heartTFs)), " TF target sets...\n"))}
  
  n = n+1
  
}


iCard_heartNet_regulonStats_targRPMnorm20[is.na(iCard_heartNet_regulonStats_targRPMnorm20)]<-0
iCard_heartNet_regulonStats_minRPM20_targRPMnorm20[is.na(iCard_heartNet_regulonStats_minRPM20_targRPMnorm20)]<-0

iCard_heartNet_regulonStats_tall_targRPMnorm20 <- iCard_heartNet_regulonStats_targRPMnorm20 %>%
  group_by(TF, tf.meanTPM) %>%
  gather("regulonStat", "value", 1:(ncol(iCard_heartNet_regulonStats_targRPMnorm20) -3))

iCard_heartNet_regulonStats_tall_minRPM20_targRPMnorm20 <- iCard_heartNet_regulonStats_minRPM20_targRPMnorm20 %>%
  group_by(TF, tf.meanTPM) %>%
  gather("regulonStat", "value", 1:(ncol(iCard_heartNet_regulonStats_minRPM20_targRPMnorm20) -3))

iCtemp_reg_targRPMnorm20 <- iCard_heartNet_regulonStats_tall_minRPM20_targRPMnorm20 %>% filter(tf.meanRPM > 20, regulonStat == 'edgeWt_normNDEup_fracUp') %>% dplyr::rename(Norm_Regulon_responsiveness = value)
iCtemp_reg_gene_RPMnorm20_targRPMnorm20 <- indiv_perturbability_RPMnorm20 %>%
  filter(cellType == 'iCard') %>%
  # inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
  inner_join(iCtemp_reg_targRPMnorm20 %>% dplyr::rename(GeneSymbol = TF))

set.seed(37392)
iCard_heartNet_regulon_normnDEupfracUpEdgeWt_vsGeneLevelnormUpreg_scatter_RPMnorm20_targRPMnorm20 <- ggplot() +
  geom_jitter(data = iCtemp_reg_gene_RPMnorm20_targRPMnorm20 %>% anti_join(markerTab2 %>% filter(cellTypeMark == 'iCard'), by = c('gene_id', 'GeneSymbol')), aes(norm_nDESeq2conditionsUp, Norm_Regulon_responsiveness), width = 0.02, height = 0) +
  geom_point(data = iCtemp_reg_gene_RPMnorm20_targRPMnorm20 %>% inner_join(markerTab2 %>% filter(cellTypeMark == 'iCard'), by = c('gene_id', 'GeneSymbol')), aes(norm_nDESeq2conditionsUp, Norm_Regulon_responsiveness), color = 'red') +
  geom_text_repel(data = iCtemp_reg_gene_RPMnorm20_targRPMnorm20 %>% inner_join(markerTab2 %>% filter(cellTypeMark == 'iCard'), by = c('gene_id', 'GeneSymbol')) %>% unique(), aes(norm_nDESeq2conditionsUp, Norm_Regulon_responsiveness, label = GeneSymbol), color = 'red', size = 2) +
  theme_classic()
ggsave(iCard_heartNet_regulon_normnDEupfracUpEdgeWt_vsGeneLevelnormUpreg_scatter_RPMnorm20_targRPMnorm20, file =  paste0(graphDir,'/differentialExpression/allSamples/iCards/iCard_heartNet_regulon_normnDEupfracUpEdgeWt_vsGeneLevelnormUpreg_scatter_RPMnorm20_targRPMnorm20.pdf'), width = 3, height = 2.5, useDingbats = F)

iCard_heartNet_regulon_norm_nDEupfracUpEdgeWt_minRPM20_histograms_targRPMnorm20_noLabs <- ggplot() +
  geom_histogram(data = iCard_heartNet_regulonStats_minRPM20_targRPMnorm20 %>% filter(tf.meanRPM > 20), aes(edgeWt_normNDEup_fracUp), fill = 'grey', alpha = 0.9) +
  geom_linerange(data = iCard_heartNet_regulonStats_minRPM20_targRPMnorm20 %>% dplyr::rename(GeneSymbol = TF) %>% inner_join(markerTab2, by = "GeneSymbol") %>% filter(cellTypeMark == "iCard", tf.meanTPM > 5) %>% unique(), 
                 aes(edgeWt_normNDEup_fracUp, ymin = 0, ymax = 10), size = 0.2, color = "red") +
  geom_text_repel(data = iCard_heartNet_regulonStats_minRPM20_targRPMnorm20 %>% dplyr::rename(GeneSymbol = TF) %>% inner_join(markerTab2, by = "GeneSymbol") %>% filter(cellTypeMark == "iCard", tf.meanTPM > 5) %>% unique(), 
                  aes(x = edgeWt_normNDEup_fracUp, y = 10, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  # ylim(c(-100, 200)) +
  # facet_grid(. ~ regulonStat, scales = 'free_x') +
  theme(axis.text = element_blank(),
        axis.title = element_blank())

iCard_heartNet_regulon_norm_nDEupfracUpEdgeWt_minRPM20_histograms_targRPMnorm20_wLabs <- ggplot() +
  geom_histogram(data = iCard_heartNet_regulonStats_minRPM20_targRPMnorm20 %>% filter(tf.meanRPM > 20), aes(edgeWt_normNDEup_fracUp), fill = 'grey', alpha = 0.9) +
  geom_linerange(data = iCard_heartNet_regulonStats_minRPM20_targRPMnorm20 %>% dplyr::rename(GeneSymbol = TF) %>% inner_join(markerTab2, by = "GeneSymbol") %>% filter(cellTypeMark == "iCard", tf.meanTPM > 5) %>% unique(), 
                 aes(edgeWt_normNDEup_fracUp, ymin = 0, ymax = 10), size = 0.2, color = "red") +
  geom_text_repel(data = iCard_heartNet_regulonStats_minRPM20_targRPMnorm20 %>% dplyr::rename(GeneSymbol = TF) %>% inner_join(markerTab2, by = "GeneSymbol") %>% filter(cellTypeMark == "iCard", tf.meanTPM > 5) %>% unique(), 
                  aes(x = edgeWt_normNDEup_fracUp, y = 10, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() 




#### Fibroblast
## use integument network

## fraction of expressed targets in integument dysregulated
fibro_integNet_regulonStats = list()
fibro_integNet_regulonStats_minRPM20 = list()
n = 1
for (tf in integTFs) {
  
  tfMeanTPM = indiv_perturbability %>%
    filter(cellType == "GM00942", GeneSymbol == tf) %>%
    dplyr::select(GeneSymbol, meanTPM, meanRPM) %>%
    dplyr::rename(TF = GeneSymbol,
                  tf.meanTPM = meanTPM,
                  tf.meanRPM = meanRPM)
  
  tempNet = integNet %>%
    filter(TF == tf) %>%
    mutate(GeneSymbol = Target) %>%
    dplyr::select(GeneSymbol, edgeWt)
  
  temp_heartNet_regulonStats = indiv_perturbability %>%
    filter(cellType == "GM00942", meanTPM > 10) %>%
    inner_join(tempNet, by = "GeneSymbol") %>%
    summarise(edgeWt_NDEup = sum(nDESeq2conditionsUp * edgeWt, na.rm = T)/sum(edgeWt),
              edgeWt_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_NDEup = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM), na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM) * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt)) %>%
    mutate(TF = tf) %>%
    left_join(tfMeanTPM, by = "TF")
  
  temp_heartNet_regulonStats_minRPM20 = indiv_perturbability %>%
    filter(cellType == "GM00942", meanRPM > 20) %>%
    inner_join(tempNet, by = "GeneSymbol") %>%
    summarise(edgeWt_NDEup = sum(nDESeq2conditionsUp * edgeWt, na.rm = T)/sum(edgeWt),
              edgeWt_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_NDEup = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM), na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM) * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt)) %>%
    mutate(TF = tf) %>%
    left_join(tfMeanTPM, by = "TF")
  
  if(is.null(dim(fibro_integNet_regulonStats))){
    fibro_integNet_regulonStats = temp_heartNet_regulonStats
    fibro_integNet_regulonStats_minRPM20 = temp_heartNet_regulonStats_minRPM20
  } else {
    fibro_integNet_regulonStats = bind_rows(fibro_integNet_regulonStats, temp_heartNet_regulonStats)
    fibro_integNet_regulonStats_minRPM20 = bind_rows(fibro_integNet_regulonStats_minRPM20, temp_heartNet_regulonStats_minRPM20)
  }
  
  if(n %% 25 == 0){cat(paste0("Completed ", as.character(n), "/", as.character(length(integTFs)), " TF target sets...\n"))}
  
  n = n+1
  
}


fibro_integNet_regulonStats[is.na(fibro_integNet_regulonStats)]<-0
fibro_integNet_regulonStats_minRPM20[is.na(fibro_integNet_regulonStats_minRPM20)]<-0

fibro_integNet_regulonStats_tall <- fibro_integNet_regulonStats %>%
  group_by(TF, tf.meanTPM) %>%
  gather("regulonStat", "value", 1:(ncol(fibro_integNet_regulonStats) -3))

fibro_integNet_regulonStats_tall_minRPM20 <- fibro_integNet_regulonStats_minRPM20 %>%
  group_by(TF, tf.meanTPM) %>%
  gather("regulonStat", "value", 1:(ncol(fibro_integNet_regulonStats) -3))

fibintegtemp_reg <- fibro_integNet_regulonStats_tall_minRPM20 %>% filter(tf.meanRPM > 20, regulonStat == 'edgeWt_NDEup_fracUp') %>% dplyr::rename(Regulon_responsiveness_integ = value)
fibinteg_reg_gene_RPMnorm20 <- indiv_perturbability_RPMnorm20 %>%
  filter(cellType == 'GM00942') %>%
  # inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
  inner_join(fibintegtemp_reg %>% dplyr::rename(GeneSymbol = TF))

## fraction of expressed targets in normalSkin dysregulated
fibro_nlSkinNet_regulonStats = list()
fibro_nlSkinNet_regulonStats_minRPM20 = list()
n = 1
for (tf in nlSkinTFs) {
  
  tfMeanTPM = indiv_perturbability %>%
    filter(cellType == "GM00942", GeneSymbol == tf) %>%
    dplyr::select(GeneSymbol, meanTPM, meanRPM) %>%
    dplyr::rename(TF = GeneSymbol,
                  tf.meanTPM = meanTPM,
                  tf.meanRPM = meanRPM)
  
  tempNet = nlSkinNet %>%
    filter(TF == tf) %>%
    mutate(GeneSymbol = Target) %>%
    dplyr::select(GeneSymbol, edgeWt)
  
  temp_heartNet_regulonStats = indiv_perturbability %>%
    filter(cellType == "GM00942", meanTPM > 10) %>%
    inner_join(tempNet, by = "GeneSymbol") %>%
    summarise(edgeWt_NDEup = sum(nDESeq2conditionsUp * edgeWt, na.rm = T)/sum(edgeWt),
              edgeWt_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_NDEup = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM), na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM) * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt)) %>%
    mutate(TF = tf) %>%
    left_join(tfMeanTPM, by = "TF")
  
  temp_heartNet_regulonStats_minRPM20 = indiv_perturbability %>%
    filter(cellType == "GM00942", meanRPM > 20) %>%
    inner_join(tempNet, by = "GeneSymbol") %>%
    summarise(edgeWt_NDEup = sum(nDESeq2conditionsUp * edgeWt, na.rm = T)/sum(edgeWt),
              edgeWt_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_NDEup = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM), na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM) * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt)) %>%
    mutate(TF = tf) %>%
    left_join(tfMeanTPM, by = "TF")
  
  if(is.null(dim(fibro_nlSkinNet_regulonStats))){
    fibro_nlSkinNet_regulonStats = temp_heartNet_regulonStats
    fibro_nlSkinNet_regulonStats_minRPM20 = temp_heartNet_regulonStats_minRPM20
  } else {
    fibro_nlSkinNet_regulonStats = bind_rows(fibro_nlSkinNet_regulonStats, temp_heartNet_regulonStats)
    fibro_nlSkinNet_regulonStats_minRPM20 = bind_rows(fibro_nlSkinNet_regulonStats_minRPM20, temp_heartNet_regulonStats_minRPM20)
  }
  
  if(n %% 25 == 0){cat(paste0("Completed ", as.character(n), "/", as.character(length(nlSkinTFs)), " TF target sets...\n"))}
  
  n = n+1
  
}


fibro_nlSkinNet_regulonStats[is.na(fibro_nlSkinNet_regulonStats)]<-0
fibro_nlSkinNet_regulonStats_minRPM20[is.na(fibro_nlSkinNet_regulonStats_minRPM20)]<-0

fibro_nlSkinNet_regulonStats_tall <- fibro_nlSkinNet_regulonStats %>%
  group_by(TF, tf.meanTPM) %>%
  gather("regulonStat", "value", 1:(ncol(fibro_nlSkinNet_regulonStats) -3))

fibro_nlSkinNet_regulonStats_tall_minRPM20 <- fibro_nlSkinNet_regulonStats_minRPM20 %>%
  group_by(TF, tf.meanTPM) %>%
  gather("regulonStat", "value", 1:(ncol(fibro_nlSkinNet_regulonStats) -3))

fibnlSkintemp_reg <- fibro_nlSkinNet_regulonStats_tall_minRPM20 %>% filter(tf.meanRPM > 20, regulonStat == 'edgeWt_NDEup_fracUp') %>% dplyr::rename(Regulon_responsiveness_nlSkin = value)
fibnlSkin_reg_gene_RPMnorm20 <- indiv_perturbability_RPMnorm20 %>%
  filter(cellType == 'GM00942') %>%
  # inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
  inner_join(fibnlSkintemp_reg %>% dplyr::rename(GeneSymbol = TF))


## fraction of expressed targets in dermal fibroblast dysregulated
fibro_dermNet_regulonStats = list()
fibro_dermNet_regulonStats_minRPM20 = list()
n = 1
for (tf in dermTFs) {
  
  tfMeanTPM = indiv_perturbability %>%
    filter(cellType == "GM00942", GeneSymbol == tf) %>%
    dplyr::select(GeneSymbol, meanTPM, meanRPM) %>%
    dplyr::rename(TF = GeneSymbol,
                  tf.meanTPM = meanTPM,
                  tf.meanRPM = meanRPM)
  
  tempNet = dermNet %>%
    filter(TF == tf) %>%
    mutate(GeneSymbol = Target) %>%
    dplyr::select(GeneSymbol, edgeWt)
  
  temp_heartNet_regulonStats = indiv_perturbability %>%
    filter(cellType == "GM00942", meanTPM > 10) %>%
    inner_join(tempNet, by = "GeneSymbol") %>%
    summarise(edgeWt_NDEup = sum(nDESeq2conditionsUp * edgeWt, na.rm = T)/sum(edgeWt),
              edgeWt_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_NDEup = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM), na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM) * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt)) %>%
    mutate(TF = tf) %>%
    left_join(tfMeanTPM, by = "TF")
  
  temp_heartNet_regulonStats_minRPM20 = indiv_perturbability %>%
    filter(cellType == "GM00942", meanRPM > 20) %>%
    inner_join(tempNet, by = "GeneSymbol") %>%
    summarise(edgeWt_NDEup = sum(nDESeq2conditionsUp * edgeWt, na.rm = T)/sum(edgeWt),
              edgeWt_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_NDEup = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM), na.rm = T)/sum(edgeWt),
              edgeWt_meanTarg_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM) * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt)) %>%
    mutate(TF = tf) %>%
    left_join(tfMeanTPM, by = "TF")
  
  if(is.null(dim(fibro_dermNet_regulonStats))){
    fibro_dermNet_regulonStats = temp_heartNet_regulonStats
    fibro_dermNet_regulonStats_minRPM20 = temp_heartNet_regulonStats_minRPM20
  } else {
    fibro_dermNet_regulonStats = bind_rows(fibro_dermNet_regulonStats, temp_heartNet_regulonStats)
    fibro_dermNet_regulonStats_minRPM20 = bind_rows(fibro_dermNet_regulonStats_minRPM20, temp_heartNet_regulonStats_minRPM20)
  }
  
  if(n %% 25 == 0){cat(paste0("Completed ", as.character(n), "/", as.character(length(dermTFs)), " TF target sets...\n"))}
  
  n = n+1
  
}


fibro_dermNet_regulonStats[is.na(fibro_dermNet_regulonStats)]<-0
fibro_dermNet_regulonStats_minRPM20[is.na(fibro_dermNet_regulonStats_minRPM20)]<-0

fibro_dermNet_regulonStats_tall <- fibro_dermNet_regulonStats %>%
  group_by(TF, tf.meanTPM) %>%
  gather("regulonStat", "value", 1:(ncol(fibro_dermNet_regulonStats) -3))

fibro_dermNet_regulonStats_tall_minRPM20 <- fibro_dermNet_regulonStats_minRPM20 %>%
  group_by(TF, tf.meanTPM) %>%
  gather("regulonStat", "value", 1:(ncol(fibro_dermNet_regulonStats) -3))

fibdermtemp_reg <- fibro_dermNet_regulonStats_tall_minRPM20 %>% filter(tf.meanRPM > 20, regulonStat == 'edgeWt_NDEup_fracUp') %>% dplyr::rename(Regulon_responsiveness_derm = value)
fibderm_reg_gene_RPMnorm20 <- indiv_perturbability_RPMnorm20 %>%
  filter(cellType == 'GM00942') %>%
  # inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
  inner_join(fibdermtemp_reg %>% dplyr::rename(GeneSymbol = TF))

set.seed(37392)
fibro_dermNet_regulon_nDEupfracUpEdgeWt_vsGeneLevelUpreg_scatter_minRPM20 <- ggplot() +
  geom_jitter(data = fibderm_reg_gene_RPMnorm20 %>% filter(!(GeneSymbol %in% targs16)), aes(nDESeq2conditionsUp, Regulon_responsiveness_derm), width = 0.2, height = 0) +
  geom_point(data = fibderm_reg_gene_RPMnorm20 %>% filter(GeneSymbol %in% targs8), aes(nDESeq2conditionsUp, Regulon_responsiveness_derm), color = 'darkolivegreen3') +
  geom_text_repel(data = fibderm_reg_gene_RPMnorm20 %>% filter(GeneSymbol %in% targs8), aes(nDESeq2conditionsUp, Regulon_responsiveness_derm, label = GeneSymbol), color = 'darkolivegreen3') +
  geom_point(data = fibderm_reg_gene_RPMnorm20 %>% filter(GeneSymbol %in% targs16 & !(GeneSymbol %in% targs8)), aes(nDESeq2conditionsUp, Regulon_responsiveness_derm), color = 'blue') +
  geom_text_repel(data = fibderm_reg_gene_RPMnorm20 %>% filter(GeneSymbol %in% targs16 & !(GeneSymbol %in% targs8)), aes(nDESeq2conditionsUp, Regulon_responsiveness_derm, label = GeneSymbol), color = 'blue') +
  theme_classic() +
  theme(axis.text = element_text(color = 'black'))
ggsave(fibro_dermNet_regulon_nDEupfracUpEdgeWt_vsGeneLevelUpreg_scatter_minRPM20, file =  paste0(graphDir,'/differentialExpression/allSamples/fibroblasts/fibro_dermNet_regulon_nDEupfracUpEdgeWt_vsGeneLevelUpreg_scatter_minRPM20.pdf'), width = 3, height = 2.5, useDingbats = F)

set.seed(37351)
fibro_nlSkinNet_regulon_nDEupfracUpEdgeWt_vsGeneLevelUpreg_scatter_minRPM20 <- ggplot() +
  geom_jitter(data = fibnlSkin_reg_gene_RPMnorm20 %>% filter(!(GeneSymbol %in% targs16)), aes(nDESeq2conditionsUp, Regulon_responsiveness_nlSkin), width = 0.2, height = 0) +
  geom_point(data = fibnlSkin_reg_gene_RPMnorm20 %>% filter(GeneSymbol %in% targs8), aes(nDESeq2conditionsUp, Regulon_responsiveness_nlSkin), color = 'darkolivegreen3') +
  geom_text_repel(data = fibnlSkin_reg_gene_RPMnorm20 %>% filter(GeneSymbol %in% targs8), aes(nDESeq2conditionsUp, Regulon_responsiveness_nlSkin, label = GeneSymbol), color = 'darkolivegreen3') +
  geom_point(data = fibnlSkin_reg_gene_RPMnorm20 %>% filter(GeneSymbol %in% targs16 & !(GeneSymbol %in% targs8)), aes(nDESeq2conditionsUp, Regulon_responsiveness_nlSkin), color = 'blue') +
  geom_text_repel(data = fibnlSkin_reg_gene_RPMnorm20 %>% filter(GeneSymbol %in% targs16 & !(GeneSymbol %in% targs8)), aes(nDESeq2conditionsUp, Regulon_responsiveness_nlSkin, label = GeneSymbol), color = 'blue') +
  theme_classic() +
  theme(axis.text = element_text(color = 'black'))
ggsave(fibro_nlSkinNet_regulon_nDEupfracUpEdgeWt_vsGeneLevelUpreg_scatter_minRPM20, file =  paste0(graphDir,'/differentialExpression/allSamples/fibroblasts/fibro_nlSkinNet_regulon_nDEupfracUpEdgeWt_vsGeneLevelUpreg_scatter_minRPM20.pdf'), width = 3, height = 2.5, useDingbats = F)

set.seed(39392)
fibro_integNet_regulon_nDEupfracUpEdgeWt_vsGeneLevelUpreg_scatter_minRPM20 <- ggplot() +
  geom_jitter(data = fibinteg_reg_gene_RPMnorm20 %>% filter(!(GeneSymbol %in% targs16)), aes(nDESeq2conditionsUp, Regulon_responsiveness_integ), width = 0.2, height = 0) +
  geom_point(data = fibinteg_reg_gene_RPMnorm20 %>% filter(GeneSymbol %in% targs8), aes(nDESeq2conditionsUp, Regulon_responsiveness_integ), color = 'darkolivegreen3') +
  geom_text_repel(data = fibinteg_reg_gene_RPMnorm20 %>% filter(GeneSymbol %in% targs8), aes(nDESeq2conditionsUp, Regulon_responsiveness_integ, label = GeneSymbol), color = 'darkolivegreen3') +
  geom_point(data = fibinteg_reg_gene_RPMnorm20 %>% filter(GeneSymbol %in% targs16 & !(GeneSymbol %in% targs8)), aes(nDESeq2conditionsUp, Regulon_responsiveness_integ), color = 'blue') +
  geom_text_repel(data = fibinteg_reg_gene_RPMnorm20 %>% filter(GeneSymbol %in% targs16 & !(GeneSymbol %in% targs8)), aes(nDESeq2conditionsUp, Regulon_responsiveness_integ, label = GeneSymbol), color = 'blue') +
  theme_classic() +
  theme(axis.text = element_text(color = 'black'))
ggsave(fibro_integNet_regulon_nDEupfracUpEdgeWt_vsGeneLevelUpreg_scatter_minRPM20, file =  paste0(graphDir,'/differentialExpression/allSamples/fibroblasts/fibro_integNet_regulon_nDEupfracUpEdgeWt_vsGeneLevelUpreg_scatter_minRPM20.pdf'), width = 3, height = 2.5, useDingbats = F)


fibro_dermNet_regulon_nDEupfracUpEdgeWt_minRPM20_histograms_noLabs <- ggplot() +
  geom_histogram(data = fibro_dermNet_regulonStats_tall_minRPM20 %>% filter(tf.meanRPM > 20, regulonStat == 'edgeWt_NDEup_fracUp'), aes(value), fill = 'grey', alpha = 0.9) +
  geom_linerange(data = fibro_dermNet_regulonStats_tall_minRPM20 %>% dplyr::rename(GeneSymbol = TF) %>% filter(GeneSymbol %in% targs16, regulonStat == 'edgeWt_NDEup_fracUp') %>% unique(), 
                 aes(value, ymin = 0, ymax = 10), size = 0.2, color = "darkolivegreen3") +
  geom_text_repel(data = fibro_dermNet_regulonStats_tall_minRPM20 %>% dplyr::rename(GeneSymbol = TF) %>% filter(GeneSymbol %in% targs16, regulonStat == 'edgeWt_NDEup_fracUp') %>% unique(), 
                  aes(x = value, y = 10, label = GeneSymbol), color = "darkolivegreen3") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  # ylim(c(-100, 200)) +
  facet_grid(. ~ regulonStat, scales = 'free_x') +
  theme(axis.text = element_blank(),
        axis.title = element_blank())
ggsave(fibro_dermNet_regulon_nDEupfracUpEdgeWt_minRPM20_histograms_noLabs, file =  paste0(graphDir,'/differentialExpression/allSamples/fibroblasts/fibro_dermNet_regulon_nDEupfracUpEdgeWt_minRPM20_histograms_noLabs.pdf'), width = 4, height = 2, useDingbats = F, units = 'in')

fibro_dermNet_regulon_nDEupfracUpEdgeWt_minRPM20_histograms_wLabs <- ggplot() +
  geom_histogram(data = fibro_dermNet_regulonStats_tall_minRPM20 %>% filter(tf.meanRPM > 20, regulonStat == 'edgeWt_NDEup_fracUp'), aes(value), fill = 'grey', alpha = 0.9) +
  geom_linerange(data = fibro_dermNet_regulonStats_tall_minRPM20 %>% dplyr::rename(GeneSymbol = TF) %>% filter(GeneSymbol %in% targs16, regulonStat == 'edgeWt_NDEup_fracUp') %>% unique(), 
                 aes(value, ymin = 0, ymax = 10), size = 0.2, color = "darkolivegreen3") +
  geom_text_repel(data = fibro_dermNet_regulonStats_tall_minRPM20 %>% dplyr::rename(GeneSymbol = TF) %>% filter(GeneSymbol %in% targs16, regulonStat == 'edgeWt_NDEup_fracUp') %>% unique(), 
                  aes(x = value, y = 10, label = GeneSymbol), color = "darkolivegreen3") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  # ylim(c(-100, 200)) +
  facet_grid(. ~ regulonStat, scales = 'free_x') +
  theme(axis.text = element_text(color = 'black'))
ggsave(fibro_dermNet_regulon_nDEupfracUpEdgeWt_minRPM20_histograms_wLabs, file =  paste0(graphDir,'/differentialExpression/allSamples/fibroblasts/fibro_dermNet_regulon_nDEupfracUpEdgeWt_minRPM20_histograms_wLabs.pdf'), width = 4, height = 2.5, useDingbats = F)


# # normalized responsiveness of all target genes
# iCard_heartNet_regulonStats_targRPMnorm20 = list()
# iCard_heartNet_regulonStats_minRPM20_targRPMnorm20 = list()
# n = 1
# for (tf in heartTFs) {
#   
#   tfMeanTPM = indiv_perturbability_RPMnorm20 %>%
#     filter(cellType == "iCard", GeneSymbol == tf) %>%
#     dplyr::select(GeneSymbol, meanTPM, meanRPM) %>%
#     dplyr::rename(TF = GeneSymbol,
#                   tf.meanTPM = meanTPM,
#                   tf.meanRPM = meanRPM)
#   
#   tempNet = heartNet %>%
#     filter(TF == tf) %>%
#     mutate(GeneSymbol = Target) %>%
#     dplyr::select(GeneSymbol, edgeWt)
#   
#   temp_heartNet_regulonStats = indiv_perturbability_RPMnorm20 %>%
#     filter(cellType == "iCard") %>%
#     inner_join(tempNet, by = "GeneSymbol") %>%
#     summarise(edgeWt_NDEup = sum(nDESeq2conditionsUp * edgeWt, na.rm = T)/sum(edgeWt),
#               edgeWt_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt),
#               edgeWt_meanTarg_NDEup = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM), na.rm = T)/sum(edgeWt),
#               edgeWt_meanTarg_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM) * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt),
#               edgeWt_normNDEup = sum(norm_nDESeq2conditionsUp * edgeWt, na.rm = T)/sum(edgeWt),
#               edgeWt_normNDEup_fracUp = sum(norm_nDESeq2conditionsUp * edgeWt * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt),
#               edgeWt_meanTarg_normNDEup = sum(norm_nDESeq2conditionsUp * edgeWt * log(meanTPM), na.rm = T)/sum(edgeWt),
#               edgeWt_meanTarg_normNDEup_fracUp = sum(norm_nDESeq2conditionsUp * edgeWt * log(meanTPM) * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt)) %>%
#     mutate(TF = tf) %>%
#     left_join(tfMeanTPM, by = "TF")
#   
#   temp_heartNet_regulonStats_minRPM20 = indiv_perturbability_RPMnorm20 %>%
#     filter(cellType == "iCard", meanRPM > 20) %>%
#     inner_join(tempNet, by = "GeneSymbol") %>%
#     summarise(edgeWt_NDEup = sum(nDESeq2conditionsUp * edgeWt, na.rm = T)/sum(edgeWt),
#               edgeWt_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt),
#               edgeWt_meanTarg_NDEup = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM), na.rm = T)/sum(edgeWt),
#               edgeWt_meanTarg_NDEup_fracUp = sum(nDESeq2conditionsUp * edgeWt * log(meanTPM) * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt),
#               edgeWt_normNDEup = sum(norm_nDESeq2conditionsUp * edgeWt, na.rm = T)/sum(edgeWt),
#               edgeWt_normNDEup_fracUp = sum(norm_nDESeq2conditionsUp * edgeWt * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt),
#               edgeWt_meanTarg_normNDEup = sum(norm_nDESeq2conditionsUp * edgeWt * log(meanTPM), na.rm = T)/sum(edgeWt),
#               edgeWt_meanTarg_normNDEup_fracUp = sum(norm_nDESeq2conditionsUp * edgeWt * log(meanTPM) * nDESeq2conditionsUp/nDESeq2conditionsAll, na.rm = T)/sum(edgeWt)) %>%
#     mutate(TF = tf) %>%
#     left_join(tfMeanTPM, by = "TF")
#   
#   if(is.null(dim(iCard_heartNet_regulonStats_targRPMnorm20))){
#     iCard_heartNet_regulonStats_targRPMnorm20 = temp_heartNet_regulonStats
#     iCard_heartNet_regulonStats_minRPM20_targRPMnorm20 = temp_heartNet_regulonStats_minRPM20
#   } else {
#     iCard_heartNet_regulonStats_targRPMnorm20 = bind_rows(iCard_heartNet_regulonStats_targRPMnorm20, temp_heartNet_regulonStats)
#     iCard_heartNet_regulonStats_minRPM20_targRPMnorm20 = bind_rows(iCard_heartNet_regulonStats_minRPM20_targRPMnorm20, temp_heartNet_regulonStats_minRPM20)
#   }
#   
#   if(n %% 25 == 0){cat(paste0("Completed ", as.character(n), "/", as.character(length(heartTFs)), " TF target sets...\n"))}
#   
#   n = n+1
#   
# }
# 
# 
# iCard_heartNet_regulonStats_targRPMnorm20[is.na(iCard_heartNet_regulonStats_targRPMnorm20)]<-0
# iCard_heartNet_regulonStats_minRPM20_targRPMnorm20[is.na(iCard_heartNet_regulonStats_minRPM20_targRPMnorm20)]<-0
# 
# iCard_heartNet_regulonStats_tall_targRPMnorm20 <- iCard_heartNet_regulonStats_targRPMnorm20 %>%
#   group_by(TF, tf.meanTPM) %>%
#   gather("regulonStat", "value", 1:(ncol(iCard_heartNet_regulonStats_targRPMnorm20) -3))
# 
# iCard_heartNet_regulonStats_tall_minRPM20_targRPMnorm20 <- iCard_heartNet_regulonStats_minRPM20_targRPMnorm20 %>%
#   group_by(TF, tf.meanTPM) %>%
#   gather("regulonStat", "value", 1:(ncol(iCard_heartNet_regulonStats_minRPM20_targRPMnorm20) -3))
# 
# iCtemp_reg_targRPMnorm20 <- iCard_heartNet_regulonStats_tall_minRPM20_targRPMnorm20 %>% filter(tf.meanRPM > 20, regulonStat == 'edgeWt_normNDEup_fracUp') %>% dplyr::rename(Norm_Regulon_responsiveness = value)
# iCtemp_reg_gene_RPMnorm20_targRPMnorm20 <- indiv_perturbability_RPMnorm20 %>%
#   filter(cellType == 'iCard') %>%
#   # inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
#   inner_join(iCtemp_reg_targRPMnorm20 %>% dplyr::rename(GeneSymbol = TF))

