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

## Compare tissue-specificity of mean gene expression in GTEx data with gene perturbability in cellid

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(magrittr)

# if running manually in Rstudio, start here and use:
# projectDir = '~/Dropbox (RajLab)/Shared_IanM/cellid_201807_onward/'
# procDataSubdir = 'procDataScripted'
# graphSubdir = 'graphs'

procDataDir = paste0(projectDir, procDataSubdir)
graphDir = paste0(projectDir, graphSubdir)

setwd(projectDir)
cat("Starting GTEx_fibsAndiCards_all6_comparePerturbabilityAndSpecificity.R...\n")
source(paste0(projectDir, "analysisScripts/RNAtagUtilities.R"))
cat(paste0("Working in ", getwd(), "\n"))

tfGeneIDfile = "annotations/TF_gene_ids.csv"
tfTab = as_tibble(read.csv(tfGeneIDfile))
colnames(tfTab) = "gene_id"

geneIDfile = "annotations/hg19gene_idToGeneSymbol.tsv"
geneIDtoGeneSymbol <- as_tibble(read.csv(geneIDfile, sep='\t', stringsAsFactors = F))

markerFile = "miscellaneous/markerGeneList.csv"
markerTab = as_tibble(read.csv(markerFile, header = T, stringsAsFactors = F)) %>%
  left_join(geneIDtoGeneSymbol, by = "GeneSymbol")

barrierFile = "miscellaneous/fibro_barriers_GeneSymbols.txt"
barrierTab = as_tibble(read.table(barrierFile, sep = "\t", header = T, stringsAsFactors = F)) %>%
  left_join(geneIDtoGeneSymbol, by = "GeneSymbol")

hiFT_factorSumFile <- 'miscellaneous/hiFT_factorSummary - Sheet1.tsv'
hiFT_factorSum <- as_tibble(read.table(hiFT_factorSumFile, sep = "\t", header = T, stringsAsFactors = F)) %>%
  left_join(geneIDtoGeneSymbol, by = "GeneSymbol")

zhouOlson_activators = as_tibble(read.table("miscellaneous/ZhouOlson2017_TableS1_activators_reorg.txt", sep = "\t", header = T, stringsAsFactors = F))
zhouOlson_inhibitors = as_tibble(read.table("miscellaneous/ZhouOlson2017_TableS1_inhibitors_reorg.txt", sep = "\t", header = T, stringsAsFactors = F))
zhouOlson_allORFs = as_tibble(read.table("miscellaneous/ZhouOlson2017_TableS1_allORFs_reorg.txt",  sep = "\t", header = T, stringsAsFactors = F))


# GTEx tissue-specificity scores
gtex_jssp_smts <- as_tibble(read.table(paste0(procDataSubdir, "/GTEx/SMTS_JSsp_TissueSpecificity.txt"), sep = "\t", header = T, stringsAsFactors = F)) %>%
  dplyr::select(gene_id, Heart, Skin, Max.Specificity.GTEx) %>%
  dplyr::rename(Max.Specificity.GTEx.SMTS = Max.Specificity.GTEx)
# gtex_jssp_smtsd <- as_tibble(read.table(paste0(procDataSubdir, "/GTEx/SMTSD_JSsp_TissueSpecificity.txt"), sep = "\t", header = T, stringsAsFactors = F)) %>%
#   dplyr::select(gene_id, Heart...Left.Ventricle, Heart...Atrial.Appendage, Cells...Transformed.fibroblasts, Max.Specificity.GTEx.SMTSD)
gtex_jssp_smts_942_noSkin <- as_tibble(read.table(paste0(procDataSubdir, "/GTEx/SMTS_JSsp_TissueSpecificity_942_noSkin.txt"), sep = "\t", header = T, stringsAsFactors = F)) %>%
  dplyr::select(gene_id, Heart, GM00942, Max.Specificity.GTEx) %>%
  dplyr::rename(Max.Specificity.GTEx.SMTS_942_noSkin = Max.Specificity.GTEx,
         Heart_942_noSkin = Heart,
         GM00942_noSkin = GM00942) 

gtex_jssp_smts_iCard_noHeart <- as_tibble(read.table(paste0(procDataSubdir, "/GTEx/SMTS_JSsp_TissueSpecificity_iCard_noHeart.txt"), sep = "\t", header = T, stringsAsFactors = F)) %>%
  dplyr::select(gene_id, iCard, Skin, Max.Specificity.GTEx) %>%
  dplyr::rename(Max.Specificity.GTEx.SMTS_iCard_noHeart = Max.Specificity.GTEx,
         Skin_iCard_noHeart = Skin,
         iCard_iCard_noHeart = iCard)

gtex_jssp_smts_iCard_fibs <- as_tibble(read.table(paste0(procDataSubdir, "/GTEx/SMTS_JSsp_TissueSpecificity_iCard_fibs.txt"), sep = "\t", header = T, stringsAsFactors = F)) %>%
  dplyr::select(gene_id, iCard, GM00942, Skin, Heart, Max.Specificity.GTEx) %>%
  dplyr::rename(Max.Specificity.GTEx.SMTS_iCard_fibs = Max.Specificity.GTEx,
                GM00942_iCard_fibs = GM00942,
                iCard_iCard_fibs = iCard,
                Skin_iCard_fibs = Skin,
                Heart_iCard_fibs = Heart) 

gtex_jssp_smts_iCard_fibs_noHeart_noSkin <- as_tibble(read.table(paste0(procDataSubdir, "/GTEx/SMTS_JSsp_TissueSpecificity_iCard_fibs_noHeart_noSkin.txt"), sep = "\t", header = T, stringsAsFactors = F)) %>%
  dplyr::select(gene_id, iCard, GM00942, Max.Specificity.GTEx) %>%
  dplyr::rename(Max.Specificity.GTEx.SMTS_iCard_fibs_noHeart_noSkin = Max.Specificity.GTEx,
                GM00942_iCard_fibs_noHeart_noSkin = GM00942,
                iCard_iCard_fibs_noHeart_noSkin = iCard) 


gtex_jssp = inner_join(gtex_jssp_smts, gtex_jssp_smts_942_noSkin, by = "gene_id") %>%
  inner_join(gtex_jssp_smts_iCard_noHeart, by = "gene_id") %>%
  inner_join(gtex_jssp_smts_iCard_fibs, by = 'gene_id') %>%
  inner_join(gtex_jssp_smts_iCard_fibs_noHeart_noSkin, by = 'gene_id')

# perturbability calculated on individual samples
all_rpm_readFilt_manualFilt_variability_metrics <- as_tibble(read.table(
  paste0(procDataSubdir, "/allExperiments/all_rpm_readFilt_manualFilt_variability_metrics.txt"),
  sep = "\t", stringsAsFactors = F, header = T))

# all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics <- as_tibble(read.table(
#   paste0(procDataSubdir, "/allExperiments/all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_window20.txt"),
#   sep = "\t", stringsAsFactors = F, header = T))


# # perturbability calculated on grouped samples
# # change grouping seed if you modified the set seed in fibsAndiCards_all6_DEandPert.R
# all_groupedMeanRPM_readFilt_manualFilt_variability_metrics <- as_tibble(read.table(
#   paste0(procDataSubdir, "/allExperiments/all_grouped_meanRPM_readFilt_manualFilt_variability_metrics_seed4930.txt"),
#   sep = "\t", stringsAsFactors = F, header = T))
# 
# all_groupedMeanRPM_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics <- as_tibble(read.table(
#   paste0(procDataSubdir, "/allExperiments/all_grouped_meanRPM_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_window20_seed4930.txt"),
#   sep = "\t", stringsAsFactors = F, header = T))


## bind tables
all_rpm_readFilt_manualFilt_variability_metrics_GTEx <- all_rpm_readFilt_manualFilt_variability_metrics %>%
  inner_join(gtex_jssp, by = "gene_id")
all_rpm_readFilt_manualFilt_variability_metrics_GTEx[is.na(all_rpm_readFilt_manualFilt_variability_metrics_GTEx)] <- 0

set.seed(7538)
# all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx <- all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics %>%
#   inner_join(gtex_jssp, by = "gene_id")
# 
# all_groupedMeanRPM_readFilt_manualFilt_variability_metrics_GTEx <- all_groupedMeanRPM_readFilt_manualFilt_variability_metrics %>%
#   inner_join(gtex_jssp, by = "gene_id")
# 
# all_groupedMeanRPM_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx <- all_groupedMeanRPM_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics %>%
#   inner_join(gtex_jssp, by = "gene_id")


## plot tissue-specificity against perturbability
if(!dir.exists(paste0(graphSubdir, "/perturbability/iCards"))){
  dir.create(paste0(graphSubdir, "/perturbability/iCards"))
}
if(!dir.exists(paste0(graphSubdir, "/perturbability/fibroblasts"))){
  dir.create(paste0(graphSubdir, "/perturbability/fibroblasts"))
}
if(!dir.exists(paste0(graphSubdir, "/GTEx"))){
  dir.create(paste0(graphSubdir, "/GTEx"))
}
# iCards vs Heart(SMTS)
# iCard_GTEx_indiv_normDeltaKurtosis_vs_JSspHeart <- ggplot() + 
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "iCard"), 
#              aes(normDeltaKurtosis, Heart, size = log(controlMeanTPM)), alpha = 0.2, color = "blue") +
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "iCard") %>%
#                inner_join(markerTab, by = "gene_id") %>% filter(cellTypeMark == "iCard"), 
#              aes(normDeltaKurtosis, Heart, size = log(controlMeanTPM)), color = "red") +
#   geom_text_repel(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                     filter(cellType == "iCard") %>%
#                     inner_join(markerTab, by = "gene_id") %>% filter(cellTypeMark == "iCard") %>% unique(), 
#                   aes(normDeltaKurtosis, Heart, label = GeneSymbol), color = "red") +
#   geom_text_repel(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                     filter(cellType == "iCard") %>%
#                     anti_join(markerTab, by = "gene_id") %>%
#                     inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
#                     filter(rank(Heart)/length(Heart) > 0.98), 
#                   aes(normDeltaKurtosis, Heart, label = GeneSymbol), color = "black") +
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "iCard") %>%
#                anti_join(markerTab, by = "gene_id") %>%
#                filter(rank(Heart)/length(Heart) > 0.98), 
#              aes(normDeltaKurtosis, Heart, size = log(controlMeanTPM)), color = "black") +
#   theme_bw() + 
#   ylim(c(0,0.5)) +
#   ggtitle("TF gene perturbability vs. tissue-specificity to GTEx heart samples\nPerturbability = RPM-normalized delta(Kurtosis), individual samples, window radius = 20\nPoint size = control expression level\nRed = marker TF genes\nBlack = top 2% most tissue-specific TF genes\nBlue = all expressed TF genes (>10 TPM)")
# ggsave(iCard_GTEx_indiv_normDeltaKurtosis_vs_JSspHeart, file = paste0(graphSubdir, "/perturbability/iCards/iCard_GTEx_indiv_normDeltaKurtosis_vs_JSspHeart.pdf"), width = 8, height = 8, useDingbats = F)
# 
# iCard_GTEx_indiv_normDeltaMM_vs_JSspHeart <- ggplot() + 
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "iCard"), 
#              aes(normDeltaMM, Heart, size = log(controlMeanTPM)), alpha = 0.2, color = "blue") +
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "iCard") %>%
#                inner_join(markerTab, by = "gene_id") %>% filter(cellTypeMark == "iCard"), 
#              aes(normDeltaMM, Heart, size = log(controlMeanTPM)), color = "red") +
#   geom_text_repel(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                     filter(cellType == "iCard") %>%
#                     inner_join(markerTab, by = "gene_id") %>% filter(cellTypeMark == "iCard") %>% unique(), 
#                   aes(normDeltaMM, Heart, label = GeneSymbol), color = "red") +
#   geom_text_repel(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                     filter(cellType == "iCard") %>%
#                     anti_join(markerTab, by = "gene_id") %>%
#                     inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
#                     filter(rank(Heart)/length(Heart) > 0.98), 
#                   aes(normDeltaMM, Heart, label = GeneSymbol), color = "black") +
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "iCard") %>%
#                anti_join(markerTab, by = "gene_id") %>%
#                filter(rank(Heart)/length(Heart) > 0.98), 
#              aes(normDeltaMM, Heart, size = log(controlMeanTPM)), color = "black") +
#   theme_bw() + 
#   ylim(c(0,0.5)) +
#   ggtitle("TF gene perturbability vs. tissue-specificity to GTEx heart samples\nPerturbability = RPM-normalized MaxMin extremal frac., individual samples, window radius = 20\nPoint size = control expression level\nRed = marker TF genes\nBlack = top 2% most tissue-specific TF genes\nBlue = all expressed TF genes (>10 TPM)")
# ggsave(iCard_GTEx_indiv_normDeltaMM_vs_JSspHeart, file = paste0(graphSubdir, "/perturbability/iCards/iCard_GTEx_indiv_normDeltaMM_vs_JSspHeart.pdf"), width = 8, height = 8, useDingbats = F)
# 
# 
# # iCards vs Heart - Left Ventricle (SMTSD)
# iCard_GTEx_indiv_normDeltaKurtosis_vs_JSspHeart...Left.Ventricle <- ggplot() + 
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "iCard"), 
#              aes(normDeltaKurtosis, Heart...Left.Ventricle, size = log(controlMeanTPM)), alpha = 0.2, color = "blue") +
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "iCard") %>%
#                inner_join(markerTab, by = "gene_id") %>% filter(cellTypeMark == "iCard"), 
#              aes(normDeltaKurtosis, Heart...Left.Ventricle, size = log(controlMeanTPM)), color = "red") +
#   geom_text_repel(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                     filter(cellType == "iCard") %>%
#                     inner_join(markerTab, by = "gene_id") %>% filter(cellTypeMark == "iCard") %>% unique(), 
#                   aes(normDeltaKurtosis, Heart...Left.Ventricle, label = GeneSymbol), color = "red") +
#   geom_text_repel(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                     filter(cellType == "iCard") %>%
#                     anti_join(markerTab, by = "gene_id") %>%
#                     inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
#                     filter(rank(Heart...Left.Ventricle)/length(Heart...Left.Ventricle) > 0.98), 
#                   aes(normDeltaKurtosis, Heart...Left.Ventricle, label = GeneSymbol), color = "black") +
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "iCard") %>%
#                anti_join(markerTab, by = "gene_id") %>%
#                filter(rank(Heart...Left.Ventricle)/length(Heart...Left.Ventricle) > 0.98), 
#              aes(normDeltaKurtosis, Heart...Left.Ventricle, size = log(controlMeanTPM)), color = "black") +
#   theme_bw() + 
#   ylim(c(0,0.5)) +
#   ggtitle("TF gene perturbability vs. tissue-specificity to GTEx Heart - Left Ventricle samples\nPerturbability = RPM-normalized delta(Kurtosis), individual samples, window radius = 20\nPoint size = control expression level\nRed = marker TF genes\nBlack = top 2% most tissue-specific TF genes\nBlue = all expressed TF genes (>10 TPM)")
# ggsave(iCard_GTEx_indiv_normDeltaKurtosis_vs_JSspHeart...Left.Ventricle, file = paste0(graphSubdir, "/perturbability/iCards/iCard_GTEx_indiv_normDeltaKurtosis_vs_JSspHeartLV.pdf"), width = 8, height = 8, useDingbats = F)
# 
# iCard_GTEx_indiv_normDeltaMM_vs_JSspHeart...Left.Ventricle <- ggplot() + 
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "iCard"), 
#              aes(normDeltaMM, Heart...Left.Ventricle, size = log(controlMeanTPM)), alpha = 0.2, color = "blue") +
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "iCard") %>%
#                inner_join(markerTab, by = "gene_id") %>% filter(cellTypeMark == "iCard"), 
#              aes(normDeltaMM, Heart...Left.Ventricle, size = log(controlMeanTPM)), color = "red") +
#   geom_text_repel(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                     filter(cellType == "iCard") %>%
#                     inner_join(markerTab, by = "gene_id") %>% filter(cellTypeMark == "iCard") %>% unique(), 
#                   aes(normDeltaMM, Heart...Left.Ventricle, label = GeneSymbol), color = "red") +
#   geom_text_repel(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                     filter(cellType == "iCard") %>%
#                     anti_join(markerTab, by = "gene_id") %>%
#                     inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
#                     filter(rank(Heart...Left.Ventricle)/length(Heart...Left.Ventricle) > 0.98), 
#                   aes(normDeltaMM, Heart...Left.Ventricle, label = GeneSymbol), color = "black") +
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "iCard") %>%
#                anti_join(markerTab, by = "gene_id") %>%
#                filter(rank(Heart...Left.Ventricle)/length(Heart...Left.Ventricle) > 0.98), 
#              aes(normDeltaMM, Heart...Left.Ventricle, size = log(controlMeanTPM)), color = "black") +
#   theme_bw() + 
#   ylim(c(0,0.5)) +
#   ggtitle("TF gene perturbability vs. tissue-specificity to GTEx Heart - Left Ventricle samples\nPerturbability = RPM-normalized MaxMin extremal frac., individual samples, window radius = 20\nPoint size = control expression level\nRed = marker TF genes\nBlack = top 2% most tissue-specific TF genes\nBlue = all expressed TF genes (>10 TPM)")
# ggsave(iCard_GTEx_indiv_normDeltaMM_vs_JSspHeart...Left.Ventricle, file = paste0(graphSubdir, "/perturbability/iCards/iCard_GTEx_indiv_normDeltaMM_vs_JSspHeartLV.pdf"), width = 8, height = 8, useDingbats = F)
# 
# 
# # iCards vs Heart - Atrial Appendage (SMTSD)
# iCard_GTEx_indiv_normDeltaKurtosis_vs_JSspHeart...Atrial.Appendage <- ggplot() + 
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "iCard"), 
#              aes(normDeltaKurtosis, Heart...Atrial.Appendage, size = log(controlMeanTPM)), alpha = 0.2, color = "blue") +
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "iCard") %>%
#                inner_join(markerTab, by = "gene_id") %>% filter(cellTypeMark == "iCard"), 
#              aes(normDeltaKurtosis, Heart...Atrial.Appendage, size = log(controlMeanTPM)), color = "red") +
#   geom_text_repel(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                     filter(cellType == "iCard") %>%
#                     inner_join(markerTab, by = "gene_id") %>% filter(cellTypeMark == "iCard") %>% unique(), 
#                   aes(normDeltaKurtosis, Heart...Atrial.Appendage, label = GeneSymbol), color = "red") +
#   geom_text_repel(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                     filter(cellType == "iCard") %>%
#                     anti_join(markerTab, by = "gene_id") %>%
#                     inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
#                     filter(rank(Heart...Atrial.Appendage)/length(Heart...Atrial.Appendage) > 0.98), 
#                   aes(normDeltaKurtosis, Heart...Atrial.Appendage, label = GeneSymbol), color = "black") +
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "iCard") %>%
#                anti_join(markerTab, by = "gene_id") %>%
#                filter(rank(Heart...Atrial.Appendage)/length(Heart...Atrial.Appendage) > 0.98), 
#              aes(normDeltaKurtosis, Heart...Atrial.Appendage, size = log(controlMeanTPM)), color = "black") +
#   theme_bw() + 
#   ylim(c(0,0.5)) +
#   ggtitle("TF gene perturbability vs. tissue-specificity to GTEx Heart - Atrial Appendage samples\nPerturbability = RPM-normalized delta(Kurtosis), individual samples, window radius = 20\nPoint size = control expression level\nRed = marker TF genes\nBlack = top 2% most tissue-specific TF genes\nBlue = all expressed TF genes (>10 TPM)")
# ggsave(iCard_GTEx_indiv_normDeltaKurtosis_vs_JSspHeart...Atrial.Appendage, file = paste0(graphSubdir, "/perturbability/iCards/iCard_GTEx_indiv_normDeltaKurtosis_vs_JSspHeartAA.pdf"), width = 8, height = 8, useDingbats = F)
# 
# iCard_GTEx_indiv_normDeltaMM_vs_JSspHeart...Atrial.Appendage <- ggplot() + 
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "iCard"), 
#              aes(normDeltaMM, Heart...Atrial.Appendage, size = log(controlMeanTPM)), alpha = 0.2, color = "blue") +
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "iCard") %>%
#                inner_join(markerTab, by = "gene_id") %>% filter(cellTypeMark == "iCard"), 
#              aes(normDeltaMM, Heart...Atrial.Appendage, size = log(controlMeanTPM)), color = "red") +
#   geom_text_repel(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                     filter(cellType == "iCard") %>%
#                     inner_join(markerTab, by = "gene_id") %>% filter(cellTypeMark == "iCard") %>% unique(), 
#                   aes(normDeltaMM, Heart...Atrial.Appendage, label = GeneSymbol), color = "red") +
#   geom_text_repel(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                     filter(cellType == "iCard") %>%
#                     anti_join(markerTab, by = "gene_id") %>%
#                     inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
#                     filter(rank(Heart...Atrial.Appendage)/length(Heart...Atrial.Appendage) > 0.98), 
#                   aes(normDeltaMM, Heart...Atrial.Appendage, label = GeneSymbol), color = "black") +
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "iCard") %>%
#                anti_join(markerTab, by = "gene_id") %>%
#                filter(rank(Heart...Atrial.Appendage)/length(Heart...Atrial.Appendage) > 0.98), 
#              aes(normDeltaMM, Heart...Atrial.Appendage, size = log(controlMeanTPM)), color = "black") +
#   theme_bw() + 
#   ylim(c(0,0.5)) +
#   ggtitle("TF gene perturbability vs. tissue-specificity to GTEx Heart - Atrial Appendage samples\nPerturbability = RPM-normalized MaxMin extremal frac., individual samples, window radius = 20\nPoint size = control expression level\nRed = marker TF genes\nBlack = top 2% most tissue-specific TF genes\nBlue = all expressed TF genes (>10 TPM)")
# ggsave(iCard_GTEx_indiv_normDeltaMM_vs_JSspHeart...Atrial.Appendage, file = paste0(graphSubdir, "/perturbability/iCards/iCard_GTEx_indiv_normDeltaMM_vs_JSspHeartAA.pdf"), width = 8, height = 8, useDingbats = F)

tempMarks <- geneIDtoGeneSymbol %>%
  filter(GeneSymbol %in% c('ESRRG', 'ZFPM2', 'MESP1', 'MYOCD')) %>%
  mutate(cellTypeMark = 'iCard')
markerTab2 <- markerTab %>%
  bind_rows(tempMarks)

# nDE vs JSspHeart
iCard_GTEx_indiv_nDEAll_vs_JSspHeart <- ggplot() + 
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               inner_join(tfTab, by = 'gene_id') %>%
               filter(cellType == "iCard", meanRPM > 20), 
             aes(nDESeq2conditionsAll, Heart), alpha = 0.2) +
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "iCard", meanRPM > 20) %>%
               inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == "iCard"), 
             aes(nDESeq2conditionsAll, Heart), color = "red") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    inner_join(tfTab, by = 'gene_id') %>%
                    filter(cellType == "iCard", meanRPM > 20) %>%
                    inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == "iCard") %>% unique(), 
                  aes(nDESeq2conditionsAll, Heart, label = GeneSymbol), color = "red") +
  # geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #                   inner_join(tfTab, by = 'gene_id') %>%
  #                   filter(cellType == "iCard", meanRPM > 20) %>%
  #                   anti_join(markerTab, by = "gene_id") %>%
  #                   inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
  #                   filter(rank(Heart)/length(Heart) > 0.98), 
  #                 aes(nDESeq2conditionsAll, Heart, label = GeneSymbol), color = "black") +
  # geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #              inner_join(tfTab, by = 'gene_id') %>%
  #              filter(cellType == "iCard", meanRPM > 20) %>%
  #              anti_join(markerTab, by = "gene_id") %>%
  #              filter(rank(Heart)/length(Heart) > 0.98), 
  #            aes(nDESeq2conditionsAll, Heart), color = "black") +
  theme_bw() + 
  ylim(c(0,0.5)) +
  ggtitle("TF gene perturbability vs. tissue-specificity to GTEx heart samples\nPerturbability = RPM-normalized delta(Kurtosis), individual samples, window radius = 20\nPoint size = control expression level\nRed = marker TF genes\nBlack = top 2% most tissue-specific TF genes\nBlue = all expressed TF genes (>10 TPM)")
ggsave(iCard_GTEx_indiv_nDEAll_vs_JSspHeart, file = paste0(graphSubdir, "/GTEx/iCard_GTEx_indiv_nDEAll_vs_JSspHeart.pdf"), width = 8, height = 8, useDingbats = F)

# nDE vs JSspHeart - without iCards
set.seed(83274)
iCard_GTEx_indiv_nDEAll_vs_JSspHeart_flip <- ggplot() + 
  geom_jitter(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                anti_join(markerTab2, by = 'gene_id') %>%
               inner_join(tfTab, by = 'gene_id') %>%
               filter(cellType == "iCard", meanRPM > 20), 
             aes(Heart, nDESeq2conditionsAll), alpha = 0.2) +
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "iCard", meanRPM > 20) %>%
               inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == "iCard"), 
             aes(Heart, nDESeq2conditionsAll), color = "red") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    inner_join(tfTab, by = 'gene_id') %>%
                    filter(cellType == "iCard", meanRPM > 20) %>%
                    inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == "iCard") %>% unique(), 
                  aes(Heart, nDESeq2conditionsAll, label = GeneSymbol), color = "red") +
  # geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #                   inner_join(tfTab, by = 'gene_id') %>%
  #                   filter(cellType == "iCard", meanRPM > 20) %>%
  #                   anti_join(markerTab, by = "gene_id") %>%
  #                   inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
  #                   filter(rank(Heart)/length(Heart) > 0.98), 
  #                 aes(nDESeq2conditionsAll, Heart, label = GeneSymbol), color = "black") +
  # geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #              inner_join(tfTab, by = 'gene_id') %>%
  #              filter(cellType == "iCard", meanRPM > 20) %>%
  #              anti_join(markerTab, by = "gene_id") %>%
#              filter(rank(Heart)/length(Heart) > 0.98), 
#            aes(nDESeq2conditionsAll, Heart), color = "black") +
theme_bw() + 
  # ylim(c(0,0.5)) +
  ggtitle("TF gene perturbability vs. tissue-specificity to GTEx heart samples\nPerturbability = n drugs in which diff expr")
ggsave(iCard_GTEx_indiv_nDEAll_vs_JSspHeart_flip, file = paste0(graphSubdir, "/GTEx/iCard_GTEx_indiv_nDEAll_vs_JSspHeart_flip.pdf"), width = 8, height = 8, useDingbats = F)



set.seed(83274)
iCard_GTEx_indiv_nDEup_vs_JSspHeart <- ggplot() + 
  geom_jitter(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               inner_join(tfTab, by = 'gene_id') %>%
                anti_join(markerTab2, by = 'gene_id') %>%
               filter(cellType == "iCard", meanRPM > 20), 
             aes(nDESeq2conditionsUp, Heart), alpha = 0.2) +
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "iCard", meanRPM > 20) %>%
               inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == "iCard"), 
             aes(nDESeq2conditionsUp, Heart), color = "red") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "iCard", meanRPM > 20) %>%
                    inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == "iCard") %>% unique(), 
                  aes(nDESeq2conditionsUp, Heart, label = GeneSymbol), color = "red") +
  # geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #                   inner_join(tfTab, by = 'gene_id') %>%
  #                   filter(cellType == "iCard", meanRPM > 20) %>%
  #                   anti_join(markerTab, by = "gene_id") %>%
  #                   inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
  #                   filter(rank(Heart)/length(Heart) > 0.98), 
  #                 aes(nDESeq2conditionsAll, Heart, label = GeneSymbol), color = "black") +
  # geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #              inner_join(tfTab, by = 'gene_id') %>%
  #              filter(cellType == "iCard", meanRPM > 20) %>%
  #              anti_join(markerTab, by = "gene_id") %>%
#              filter(rank(Heart)/length(Heart) > 0.98), 
#            aes(nDESeq2conditionsAll, Heart), color = "black") +
theme_bw() + 
  ylim(c(0,0.5)) +
  ylab('Heart specificity score (JS divergence)') +
  ggtitle("TF gene perturbability vs. tissue-specificity to GTEx heart samples\nPerturbability = n drugs in which UP\n")
ggsave(iCard_GTEx_indiv_nDEup_vs_JSspHeart, file = paste0(graphSubdir, "/GTEx/iCard_GTEx_indiv_nDEup_vs_JSspHeart.pdf"), width = 8, height = 8, useDingbats = F)

set.seed(83274)
iCard_GTEx_indiv_nDEup_vs_JSspHeart <- ggplot() + 
  geom_jitter(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                inner_join(tfTab, by = 'gene_id') %>%
                anti_join(markerTab2, by = 'gene_id') %>%
                filter(cellType == "iCard", meanRPM > 20), 
              aes(nDESeq2conditionsUp, Heart), alpha = 0.2) +
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "iCard", meanRPM > 20) %>%
               inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == "iCard"), 
             aes(nDESeq2conditionsUp, Heart), color = "red") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "iCard", meanRPM > 20) %>%
                    inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == "iCard") %>% unique(), 
                  aes(nDESeq2conditionsUp, Heart, label = GeneSymbol), color = "red") +
  # geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #                   inner_join(tfTab, by = 'gene_id') %>%
  #                   filter(cellType == "iCard", meanRPM > 20) %>%
  #                   anti_join(markerTab, by = "gene_id") %>%
  #                   inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
  #                   filter(rank(Heart)/length(Heart) > 0.98), 
  #                 aes(nDESeq2conditionsAll, Heart, label = GeneSymbol), color = "black") +
  # geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #              inner_join(tfTab, by = 'gene_id') %>%
  #              filter(cellType == "iCard", meanRPM > 20) %>%
  #              anti_join(markerTab, by = "gene_id") %>%
#              filter(rank(Heart)/length(Heart) > 0.98), 
#            aes(nDESeq2conditionsAll, Heart), color = "black") +
theme_bw() + 
  ylim(c(0,0.5)) +
  ylab('Heart specificity score (JS divergence)') +
  ggtitle("TF gene perturbability vs. tissue-specificity to GTEx heart samples\nPerturbability = n drugs in which UP\n")
ggsave(iCard_GTEx_indiv_nDEup_vs_JSspHeart, file = paste0(graphSubdir, "/GTEx/iCard_GTEx_indiv_nDEup_vs_JSspHeart.pdf"), width = 8, height = 8, useDingbats = F)


iCard_GTEx_indiv_nDEup_vs_JSspHeart_flip <- ggplot() + 
  geom_jitter(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                inner_join(tfTab, by = 'gene_id') %>%
                anti_join(markerTab2, by = 'gene_id') %>%
                filter(cellType == "iCard", meanRPM > 20), 
              aes(Heart, nDESeq2conditionsUp), alpha = 0.2) +
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "iCard", meanRPM > 20) %>%
               inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == "iCard"), 
             aes(Heart, nDESeq2conditionsUp), color = "red") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "iCard", meanRPM > 20) %>%
                    inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == "iCard") %>% unique(), 
                  aes(Heart, nDESeq2conditionsUp, label = GeneSymbol), color = "red") +
  # geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #                   inner_join(tfTab, by = 'gene_id') %>%
  #                   filter(cellType == "iCard", meanRPM > 20) %>%
  #                   anti_join(markerTab, by = "gene_id") %>%
  #                   inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
  #                   filter(rank(Heart)/length(Heart) > 0.98), 
  #                 aes(nDESeq2conditionsAll, Heart, label = GeneSymbol), color = "black") +
  # geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #              inner_join(tfTab, by = 'gene_id') %>%
  #              filter(cellType == "iCard", meanRPM > 20) %>%
  #              anti_join(markerTab, by = "gene_id") %>%
#              filter(rank(Heart)/length(Heart) > 0.98), 
#            aes(nDESeq2conditionsAll, Heart), color = "black") +
theme_bw() + 
  # xlim(c(0,0.5)) +
  xlab('Heart specificity score (JS divergence)') +
  ggtitle("TF gene perturbability vs. tissue-specificity to GTEx heart samples\nPerturbability = n drugs in which UP\n")
ggsave(iCard_GTEx_indiv_nDEup_vs_JSspHeart_flip, file = paste0(graphSubdir, "/GTEx/iCard_GTEx_indiv_nDEup_vs_JSspHeart_flip.pdf"), width = 8, height = 8, useDingbats = F)


set.seed(9823)
iCard_GTEx_indiv_nDEup_vs_JSspHeart_withfibsAndiCard_flip <- ggplot() + 
  geom_jitter(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                inner_join(tfTab, by = 'gene_id') %>%
                anti_join(markerTab2, by = 'gene_id') %>%
                filter(cellType == "iCard", meanRPM > 20), 
              aes(Heart_iCard_fibs, nDESeq2conditionsUp), alpha = 0.2) +
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "iCard", meanRPM > 20) %>%
               inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == "iCard"), 
             aes(Heart_iCard_fibs, nDESeq2conditionsUp), color = "red") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "iCard", meanRPM > 20) %>%
                    inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == "iCard") %>% unique(), 
                  aes(Heart_iCard_fibs, nDESeq2conditionsUp, label = GeneSymbol), color = "red") +
  # geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #                   inner_join(tfTab, by = 'gene_id') %>%
  #                   filter(cellType == "iCard", meanRPM > 20) %>%
  #                   anti_join(markerTab, by = "gene_id") %>%
  #                   inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
  #                   filter(rank(Heart)/length(Heart) > 0.98), 
  #                 aes(nDESeq2conditionsAll, Heart, label = GeneSymbol), color = "black") +
  # geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #              inner_join(tfTab, by = 'gene_id') %>%
  #              filter(cellType == "iCard", meanRPM > 20) %>%
  #              anti_join(markerTab, by = "gene_id") %>%
#              filter(rank(Heart)/length(Heart) > 0.98), 
#            aes(nDESeq2conditionsAll, Heart), color = "black") +
theme_bw() + 
  # xlim(c(0,0.5)) +
  xlab('Heart specificity score (JS divergence)') +
  ggtitle("TF gene perturbability vs. tissue-specificity to GTEx heart samples (including GM942 fibs and iCards in dataset)\nPerturbability = n drugs in which UP\n")
ggsave(iCard_GTEx_indiv_nDEup_vs_JSspHeart_withfibsAndiCard_flip, file = paste0(graphSubdir, "/GTEx/iCard_GTEx_indiv_nDEup_vs_JSspHeart_withfibsAndiCard_flip.pdf"), width = 8, height = 8, useDingbats = F)



set.seed(9823)
iCard_GTEx_indiv_nDEup_vs_JSspiCard_withSkinAndHeart_flip <- ggplot() + 
  geom_jitter(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                inner_join(tfTab, by = 'gene_id') %>%
                anti_join(markerTab2, by = 'gene_id') %>%
                filter(cellType == "iCard", meanRPM > 20), 
              aes(iCard_iCard_fibs, nDESeq2conditionsUp), alpha = 0.2) +
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "iCard", meanRPM > 20) %>%
               inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == "iCard"), 
             aes(iCard_iCard_fibs, nDESeq2conditionsUp), color = "red") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "iCard", meanRPM > 20) %>%
                    inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == "iCard") %>% unique(), 
                  aes(iCard_iCard_fibs, nDESeq2conditionsUp, label = GeneSymbol), color = "red") +
  # geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #                   inner_join(tfTab, by = 'gene_id') %>%
  #                   filter(cellType == "iCard", meanRPM > 20) %>%
  #                   anti_join(markerTab, by = "gene_id") %>%
  #                   inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
  #                   filter(rank(Heart)/length(Heart) > 0.98), 
  #                 aes(nDESeq2conditionsAll, Heart, label = GeneSymbol), color = "black") +
  # geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #              inner_join(tfTab, by = 'gene_id') %>%
  #              filter(cellType == "iCard", meanRPM > 20) %>%
  #              anti_join(markerTab, by = "gene_id") %>%
#              filter(rank(Heart)/length(Heart) > 0.98), 
#            aes(nDESeq2conditionsAll, Heart), color = "black") +
theme_bw() + 
  # xlim(c(0,0.5)) +
  xlab('Heart specificity score (JS divergence)') +
  ggtitle("TF gene perturbability vs. tissue-specificity to iCard samples (including GTEx Heart and Skin in dataset)\nPerturbability = n drugs in which UP\n")
ggsave(iCard_GTEx_indiv_nDEup_vs_JSspiCard_withSkinAndHeart_flip, file = paste0(graphSubdir, "/GTEx/iCard_GTEx_indiv_nDEup_vs_JSspiCard_withSkinAndHeart_flip.pdf"), width = 8, height = 8, useDingbats = F)

iCard_GTEx_indiv_JSspiCard_noSkinAndHeart_cdf <- ggplot() + 
  stat_ecdf(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
              inner_join(tfTab, by = 'gene_id') %>%
              filter(cellType == "iCard", meanRPM > 20), 
            aes(iCard_iCard_fibs_noHeart_noSkin), alpha = 0.2) +
  theme_bw() + 
  xlab('iCard specificity score (JS divergence)') +
  ggtitle("Tissue-specificity to iCard samples (excluding GTEx Heart and Skin in dataset)\nCumulative distribution function\nAll TF genes >20 RPM in control iCards")
ggsave(iCard_GTEx_indiv_JSspiCard_noSkinAndHeart_cdf, file = paste0(graphSubdir, "/GTEx/iCard_GTEx_indiv_JSspiCard_noSkinAndHeart_cdf.pdf"), width = 8, height = 8, useDingbats = F)


fibro_GTEx_indiv_JSspiCard_noSkinAndHeart_cdf <- ggplot() + 
  stat_ecdf(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
              inner_join(tfTab, by = 'gene_id') %>%
              filter(cellType == "GM00942", meanRPM > 20), 
            aes(GM00942_iCard_fibs_noHeart_noSkin), alpha = 0.2) +
  theme_bw() + 
  xlab('Fibroblast specificity score (JS divergence)') +
  ggtitle("Tissue-specificity to GM00942 fibroblast samples (excluding GTEx Heart and Skin in dataset)\nCumulative distribution function\nAll TF genes >20 RPM in control iCards")
ggsave(fibro_GTEx_indiv_JSspiCard_noSkinAndHeart_cdf, file = paste0(graphSubdir, "/GTEx/fibro_GTEx_indiv_JSspiCard_noSkinAndHeart_cdf.pdf"), width = 8, height = 8, useDingbats = F)


set.seed(9823)
iCard_GTEx_indiv_nDEup_vs_JSspiCard_noSkinAndHeart_flip <- ggplot() + 
  geom_jitter(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                inner_join(tfTab, by = 'gene_id') %>%
                anti_join(markerTab2, by = 'gene_id') %>%
                filter(cellType == "iCard", meanRPM > 20), 
              aes(iCard_iCard_fibs_noHeart_noSkin, nDESeq2conditionsUp), alpha = 0.2) +
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "iCard", meanRPM > 20) %>%
               inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == "iCard"), 
             aes(iCard_iCard_fibs_noHeart_noSkin, nDESeq2conditionsUp), color = "red") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "iCard", meanRPM > 20) %>%
                    inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == "iCard") %>% unique(), 
                  aes(iCard_iCard_fibs_noHeart_noSkin, nDESeq2conditionsUp, label = GeneSymbol), color = "red") +
  # geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #                   inner_join(tfTab, by = 'gene_id') %>%
  #                   filter(cellType == "iCard", meanRPM > 20) %>%
  #                   anti_join(markerTab, by = "gene_id") %>%
  #                   inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
  #                   filter(rank(Heart)/length(Heart) > 0.98), 
  #                 aes(nDESeq2conditionsAll, Heart, label = GeneSymbol), color = "black") +
  # geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #              inner_join(tfTab, by = 'gene_id') %>%
  #              filter(cellType == "iCard", meanRPM > 20) %>%
  #              anti_join(markerTab, by = "gene_id") %>%
#              filter(rank(Heart)/length(Heart) > 0.98), 
#            aes(nDESeq2conditionsAll, Heart), color = "black") +
theme_bw() + 
  # xlim(c(0,0.5)) +
  xlab('Heart specificity score (JS divergence)') +
  ggtitle("TF gene perturbability vs. tissue-specificity to iCard samples (excluding GTEx Heart and Skin in dataset)\nPerturbability = n drugs in which UP\n")
ggsave(iCard_GTEx_indiv_nDEup_vs_JSspiCard_noSkinAndHeart_flip, file = paste0(graphSubdir, "/GTEx/iCard_GTEx_indiv_nDEup_vs_JSspiCard_noSkinAndHeart_flip.pdf"), width = 8, height = 8, useDingbats = F)



YangXiSpec_heart <- as_tibble(readxl::read_xlsx(paste0('miscellaneous/Yang_Xi_bioRxiv2018/311563-16.xlsx'))) %>%
  separate(Name, into = c('gene_id', 'ver'), sep = '\\.') %>%
  filter(tissue == 'Heart - Left Ventricle') %>%
  dplyr::rename(GeneSymbol = symbol)

iCard_GTEx_indiv_nDEup_vs_JSspHeart_flip_YX <- ggplot() + 
  geom_jitter(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                inner_join(tfTab, by = 'gene_id') %>%
                anti_join(markerTab2, by = 'gene_id') %>%
                filter(cellType == "iCard", meanRPM > 20), 
              aes(Heart, nDESeq2conditionsUp), alpha = 0.2) +
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "iCard", meanRPM > 20) %>%
               inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == "iCard"), 
             aes(Heart, nDESeq2conditionsUp), color = "red") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "iCard", meanRPM > 20) %>%
                    inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == "iCard") %>% unique(), 
                  aes(Heart, nDESeq2conditionsUp, label = GeneSymbol), color = "red") +
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "iCard", meanRPM > 20) %>%
               inner_join(YangXiSpec_heart, by = "gene_id")  %>% inner_join(tfTab, by = 'gene_id'), 
             aes(Heart, nDESeq2conditionsUp), color = "blue") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "iCard", meanRPM > 20) %>%
                    inner_join(YangXiSpec_heart, by = "gene_id") %>% unique() %>% inner_join(tfTab, by = 'gene_id'), 
                  aes(Heart, nDESeq2conditionsUp, label = GeneSymbol), color = "blue") +
  # geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #                   inner_join(tfTab, by = 'gene_id') %>%
  #                   filter(cellType == "iCard", meanRPM > 20) %>%
  #                   anti_join(markerTab, by = "gene_id") %>%
  #                   inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
  #                   filter(rank(Heart)/length(Heart) > 0.98), 
  #                 aes(nDESeq2conditionsAll, Heart, label = GeneSymbol), color = "black") +
  # geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #              inner_join(tfTab, by = 'gene_id') %>%
  #              filter(cellType == "iCard", meanRPM > 20) %>%
  #              anti_join(markerTab, by = "gene_id") %>%
#              filter(rank(Heart)/length(Heart) > 0.98), 
#            aes(nDESeq2conditionsAll, Heart), color = "black") +
theme_bw() + 
  # xlim(c(0,0.5)) +
  xlab('Heart specificity score (JS divergence)') +
  ggtitle("TF gene perturbability vs. tissue-specificity to GTEx heart samples\nPerturbability = n drugs in which UP\n")
ggsave(iCard_GTEx_indiv_nDEup_vs_JSspHeart_flip_YX, file = paste0(graphSubdir, "/GTEx/iCard_GTEx_indiv_nDEup_vs_JSspHeart_flip_YX.pdf"), width = 8, height = 8, useDingbats = F)


## targets for KD
set.seed(9823)
KDtargs <- c('HMGB2', 'MEF2C', 'NKX2-5', 'HAND1', 'SOX11', 'IRX4', 'SKIDA1', 'COPS2', 'ZFPM2')

iCard_GTEx_indiv_nDEup_vs_JSspiCard_noSkinAndHeart_flip_KDtargets <- ggplot() + 
  geom_jitter(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                inner_join(tfTab, by = 'gene_id') %>%
                inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
                filter(cellType == "iCard", meanRPM > 20, !(GeneSymbol %in% KDtargs)), 
              aes(iCard_iCard_fibs_noHeart_noSkin, nDESeq2conditionsUp), alpha = 0.2) +
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
               filter(cellType == "iCard", GeneSymbol %in% KDtargs), 
             aes(iCard_iCard_fibs_noHeart_noSkin, nDESeq2conditionsUp), color = "forestgreen") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
                    filter(cellType == "iCard", GeneSymbol %in% KDtargs), 
                  aes(iCard_iCard_fibs_noHeart_noSkin, nDESeq2conditionsUp, label = GeneSymbol), color = "forestgreen") +
  # geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #                   inner_join(tfTab, by = 'gene_id') %>%
  #                   filter(cellType == "iCard", meanRPM > 20) %>%
  #                   anti_join(markerTab, by = "gene_id") %>%
  #                   inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
  #                   filter(rank(Heart)/length(Heart) > 0.98), 
  #                 aes(nDESeq2conditionsAll, Heart, label = GeneSymbol), color = "black") +
  # geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #              inner_join(tfTab, by = 'gene_id') %>%
  #              filter(cellType == "iCard", meanRPM > 20) %>%
  #              anti_join(markerTab, by = "gene_id") %>%
#              filter(rank(Heart)/length(Heart) > 0.98), 
#            aes(nDESeq2conditionsAll, Heart), color = "black") +
theme_bw() + 
  # xlim(c(0,0.5)) +
  xlab('Heart specificity score (JS divergence)') +
  ggtitle("TF gene perturbability vs. tissue-specificity to iCard samples (excluding GTEx Heart and Skin in dataset)\nPerturbability = n drugs in which UP\n")
ggsave(iCard_GTEx_indiv_nDEup_vs_JSspiCard_noSkinAndHeart_flip_KDtargets, file = paste0(graphSubdir, "/GTEx/iCard_GTEx_indiv_nDEup_vs_JSspiCard_noSkinAndHeart_flip_KDtargets.pdf"), width = 8, height = 8, useDingbats = F)


# iCard_GTEx_indiv_nDEup_vs_JSspHeartLV_flip_YX <- ggplot() + 
#   geom_jitter(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
#                 inner_join(tfTab, by = 'gene_id') %>%
#                 anti_join(markerTab2, by = 'gene_id') %>%
#                 filter(cellType == "iCard", meanRPM > 20), 
#               aes(Heart...Left.Ventricle, nDESeq2conditionsUp), alpha = 0.2) +
#   geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
#                filter(cellType == "iCard", meanRPM > 20) %>%
#                inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == "iCard"), 
#              aes(Heart...Left.Ventricle, nDESeq2conditionsUp), color = "red") +
#   geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
#                     filter(cellType == "iCard", meanRPM > 20) %>%
#                     inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == "iCard") %>% unique(), 
#                   aes(Heart...Left.Ventricle, nDESeq2conditionsUp, label = GeneSymbol), color = "red") +
#   geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
#                filter(cellType == "iCard", meanRPM > 20) %>%
#                inner_join(YangXiSpec_heart, by = "gene_id")  %>% inner_join(tfTab, by = 'gene_id'), 
#              aes(Heart...Left.Ventricle, nDESeq2conditionsUp), color = "blue") +
#   geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
#                     filter(cellType == "iCard", meanRPM > 20) %>%
#                     inner_join(YangXiSpec_heart, by = "gene_id") %>% unique() %>% inner_join(tfTab, by = 'gene_id'), 
#                   aes(Heart...Left.Ventricle, nDESeq2conditionsUp, label = GeneSymbol), color = "blue") +
#   # geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
#   #                   inner_join(tfTab, by = 'gene_id') %>%
#   #                   filter(cellType == "iCard", meanRPM > 20) %>%
#   #                   anti_join(markerTab, by = "gene_id") %>%
#   #                   inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
#   #                   filter(rank(Heart)/length(Heart) > 0.98), 
#   #                 aes(nDESeq2conditionsAll, Heart, label = GeneSymbol), color = "black") +
#   # geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
#   #              inner_join(tfTab, by = 'gene_id') %>%
#   #              filter(cellType == "iCard", meanRPM > 20) %>%
#   #              anti_join(markerTab, by = "gene_id") %>%
# #              filter(rank(Heart)/length(Heart) > 0.98), 
# #            aes(nDESeq2conditionsAll, Heart), color = "black") +
# theme_bw() + 
#   # xlim(c(0,0.5)) +
#   xlab('Heart - Left Ventricle specificity score (JS divergence)') +
#   ggtitle("TF gene perturbability vs. tissue-specificity to GTEx heart LV samples\nPerturbability = n drugs in which UP\n")
# ggsave(iCard_GTEx_indiv_nDEup_vs_JSspHeartLV_flip_YX, file = paste0(graphSubdir, "/GTEx/iCard_GTEx_indiv_nDEup_vs_JSspHeartLV_flip_YX.pdf"), width = 8, height = 8)
# 
# 
# iCard_GTEx_indiv_JSspHeartLV_hist_YX <- ggplot() +
#   geom_histogram(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
#                    inner_join(tfTab, by = 'gene_id') %>%
#                    filter(cellType == "iCard", meanRPM > 20),
#                  aes(Heart...Left.Ventricle), binwidth = 0.01) +
#   geom_vline(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
#                filter(cellType == "iCard", meanRPM > 20) %>%
#                inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == "iCard"), 
#              aes(xintercept = Heart...Left.Ventricle), color = 'red')
# plot(iCard_GTEx_indiv_JSspHeartLV_hist_YX)

# GM00942 fibroblasts vs. GTEx
# fibroblasts vs Skin(SMTS)
# fibro_GTEx_indiv_normDeltaKurtosis_vs_JSspSkin <- ggplot() + 
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "GM00942"), 
#              aes(normDeltaKurtosis, Skin, size = log(controlMeanTPM)), alpha = 0.2, color = "blue") +
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "GM00942") %>%
#                inner_join(barrierTab, by = "gene_id"), 
#              aes(normDeltaKurtosis, Skin, size = log(controlMeanTPM)), color = "darkolivegreen3") +
#   geom_text_repel(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                     filter(cellType == "GM00942") %>%
#                     inner_join(barrierTab, by = "gene_id"), 
#                   aes(normDeltaKurtosis, Skin, label = GeneSymbol), color = "darkolivegreen3") +
#   geom_text_repel(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                     filter(cellType == "GM00942") %>%
#                     anti_join(barrierTab, by = "gene_id") %>%
#                     inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
#                     filter(rank(Skin)/length(Skin) > 0.98), 
#                   aes(normDeltaKurtosis, Skin, label = GeneSymbol), color = "black") +
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "GM00942") %>%
#                anti_join(barrierTab, by = "gene_id") %>%
#                filter(rank(Skin)/length(Skin) > 0.98), 
#              aes(normDeltaKurtosis, Skin, size = log(controlMeanTPM)), color = "black") +
#   theme_bw() + 
#   ylim(c(0,0.5)) +
#   ggtitle("TF gene perturbability in fibroblasts vs. tissue-specificity to GTEx Skin samples\nPerturbability = RPM-normalized delta(Kurtosis), individual samples, window radius = 20\nPoint size = control expression level\nRed = marker TF genes\nBlack = top 2% most tissue-specific TF genes\nBlue = all expressed TF genes (>10 TPM)")
# ggsave(fibro_GTEx_indiv_normDeltaKurtosis_vs_JSspSkin, file = paste0(graphSubdir, "/perturbability/fibroblasts/fibro_GTEx_indiv_normDeltaKurtosis_vs_JSspSkin.pdf"), width = 8, height = 8, useDingbats = F)
# 
# fibro_GTEx_indiv_normDeltaMM_vs_JSspSkin <- ggplot() + 
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "GM00942"), 
#              aes(normDeltaMM, Skin, size = log(controlMeanTPM)), alpha = 0.2, color = "blue") +
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "GM00942") %>%
#                inner_join(barrierTab, by = "gene_id"), 
#              aes(normDeltaMM, Skin, size = log(controlMeanTPM)), color = "darkolivegreen3") +
#   geom_text_repel(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                     filter(cellType == "GM00942") %>%
#                     inner_join(barrierTab, by = "gene_id"), 
#                   aes(normDeltaMM, Skin, label = GeneSymbol), color = "darkolivegreen3") +
#   geom_text_repel(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                     filter(cellType == "GM00942") %>%
#                     anti_join(barrierTab, by = "gene_id") %>%
#                     inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
#                     filter(rank(Skin)/length(Skin) > 0.98), 
#                   aes(normDeltaMM, Skin, label = GeneSymbol), color = "black") +
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "GM00942") %>%
#                anti_join(barrierTab, by = "gene_id") %>%
#                filter(rank(Skin)/length(Skin) > 0.98), 
#              aes(normDeltaMM, Skin, size = log(controlMeanTPM)), color = "black") +
#   theme_bw() + 
#   ylim(c(0,0.5)) +
#   ggtitle("TF gene perturbability in fibroblasts vs. tissue-specificity to GTEx Skin samples\nPerturbability = RPM-normalized MaxMin extremal frac., individual samples, window radius = 20\nPoint size = control expression level\nRed = marker TF genes\nBlack = top 2% most tissue-specific TF genes\nBlue = all expressed TF genes (>10 TPM)")
# ggsave(fibro_GTEx_indiv_normDeltaMM_vs_JSspSkin, file = paste0(graphSubdir, "/perturbability/fibroblasts/fibro_GTEx_indiv_normDeltaMM_vs_JSspSkin.pdf"), width = 8, height = 8, useDingbats = F)
# 
# 
# # fibroblasts vs Cells - Transformed fibroblasts (SMTSD)
# fibro_GTEx_indiv_normDeltaKurtosis_vs_JSspCells...Transformed.fibroblasts <- ggplot() + 
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "GM00942"), 
#              aes(normDeltaKurtosis, Cells...Transformed.fibroblasts, size = log(controlMeanTPM)), alpha = 0.2, color = "blue") +
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "GM00942") %>%
#                inner_join(barrierTab, by = "gene_id"), 
#              aes(normDeltaKurtosis, Cells...Transformed.fibroblasts, size = log(controlMeanTPM)), color = "darkolivegreen3") +
#   geom_text_repel(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                     filter(cellType == "GM00942") %>%
#                     inner_join(barrierTab, by = "gene_id"), 
#                   aes(normDeltaKurtosis, Cells...Transformed.fibroblasts, label = GeneSymbol), color = "darkolivegreen3") +
#   geom_text_repel(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                     filter(cellType == "GM00942") %>%
#                     anti_join(barrierTab, by = "gene_id") %>%
#                     inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
#                     filter(rank(Cells...Transformed.fibroblasts)/length(Cells...Transformed.fibroblasts) > 0.98), 
#                   aes(normDeltaKurtosis, Cells...Transformed.fibroblasts, label = GeneSymbol), color = "black") +
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "GM00942") %>%
#                anti_join(barrierTab, by = "gene_id") %>%
#                filter(rank(Cells...Transformed.fibroblasts)/length(Cells...Transformed.fibroblasts) > 0.98), 
#              aes(normDeltaKurtosis, Cells...Transformed.fibroblasts, size = log(controlMeanTPM)), color = "black") +
#   theme_bw() + 
#   ylim(c(0,0.5)) +
#   ggtitle("TF gene perturbability in fibroblasts vs. tissue-specificity to GTEx Cells - Transformed fibroblasts samples\nPerturbability = RPM-normalized delta(Kurtosis), individual samples, window radius = 20\nPoint size = control expression level\nRed = marker TF genes\nBlack = top 2% most tissue-specific TF genes\nBlue = all expressed TF genes (>10 TPM)")
# ggsave(fibro_GTEx_indiv_normDeltaKurtosis_vs_JSspCells...Transformed.fibroblasts, file = paste0(graphSubdir, "/perturbability/fibroblasts/fibro_GTEx_indiv_normDeltaKurtosis_vs_JSspTxfmFibs.pdf"), width = 8, height = 8, useDingbats = F)
# 
# fibro_GTEx_indiv_normDeltaMM_vs_JSspCells...Transformed.fibroblasts <- ggplot() + 
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "GM00942"), 
#              aes(normDeltaMM, Cells...Transformed.fibroblasts, size = log(controlMeanTPM)), alpha = 0.2, color = "blue") +
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "GM00942") %>%
#                inner_join(barrierTab, by = "gene_id"), 
#              aes(normDeltaMM, Cells...Transformed.fibroblasts, size = log(controlMeanTPM)), color = "darkolivegreen3") +
#   geom_text_repel(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                     filter(cellType == "GM00942") %>%
#                     inner_join(barrierTab, by = "gene_id"), 
#                   aes(normDeltaMM, Cells...Transformed.fibroblasts, label = GeneSymbol), color = "darkolivegreen3") +
#   geom_text_repel(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                     filter(cellType == "GM00942") %>%
#                     anti_join(barrierTab, by = "gene_id") %>%
#                     inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
#                     filter(rank(Cells...Transformed.fibroblasts)/length(Cells...Transformed.fibroblasts) > 0.98), 
#                   aes(normDeltaMM, Cells...Transformed.fibroblasts, label = GeneSymbol), color = "black") +
#   geom_point(data = all_rpm_readFilt_manualFilt_tfOnly_tpmFilt_variability_RPMnormMetrics_GTEx %>%
#                filter(cellType == "GM00942") %>%
#                anti_join(barrierTab, by = "gene_id") %>%
#                filter(rank(Cells...Transformed.fibroblasts)/length(Cells...Transformed.fibroblasts) > 0.98), 
#              aes(normDeltaMM, Cells...Transformed.fibroblasts, size = log(controlMeanTPM)), color = "black") +
#   theme_bw() + 
#   ylim(c(0,0.5)) +
#   ggtitle("TF gene perturbability in fibroblasts vs. tissue-specificity to GTEx Skin samples\nPerturbability = RPM-normalized MaxMin extremal frac., individual samples, window radius = 20\nPoint size = control expression level\nRed = marker TF genes\nBlack = top 2% most tissue-specific TF genes\nBlue = all expressed TF genes (>10 TPM)")
# ggsave(fibro_GTEx_indiv_normDeltaMM_vs_JSspCells...Transformed.fibroblasts, file = paste0(graphSubdir, "/perturbability/fibroblasts/fibro_GTEx_indiv_normDeltaMM_vs_JSspTxfmFibs.pdf"), width = 8, height = 8, useDingbats = F)

# nDE comparisons
set.seed(4286)
fibro_GTEx_indiv_nDEall_vs_JSsp_GM00942 <- ggplot() + 
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "GM00942",
                      meanRPM >= 20) %>%
               inner_join(tfTab, by = 'gene_id'), 
             aes(nDESeq2conditionsAll, GM00942), alpha = 0.2) +
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "GM00942") %>%
               inner_join(hiFT_factorSum, by = "gene_id"), 
             aes(nDESeq2conditionsAll, GM00942), color = "darkolivegreen3") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "GM00942") %>%
                    inner_join(hiFT_factorSum, by = "gene_id"), 
                  aes(nDESeq2conditionsAll, GM00942, label = GeneSymbol), color = "darkolivegreen3") +
  theme_bw() + 
  ylim(c(0,0.5)) +
  ggtitle('Fibroblast TF perturbability vs GM942 average specificity')
# plot(fibro_GTEx_indiv_nDEall_vs_JSsp_GM00942)

set.seed(4286)
fibro_GTEx_indiv_nDEup_vs_JSsp_GM00942 <- ggplot() + 
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "GM00942",
                      meanRPM >= 20) %>%
               inner_join(tfTab, by = 'gene_id'), 
             aes(nDESeq2conditionsUp, GM00942_iCard_fibs_noHeart_noSkin), alpha = 0.2) +
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "GM00942") %>%
               inner_join(hiFT_factorSum, by = "gene_id"), 
             aes(nDESeq2conditionsUp, GM00942_iCard_fibs_noHeart_noSkin), color = "forestgreen") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "GM00942") %>%
                    inner_join(hiFT_factorSum, by = "gene_id"), 
                  aes(nDESeq2conditionsUp, GM00942_iCard_fibs_noHeart_noSkin, label = GeneSymbol), color = "forestgreen") +
  theme_bw() + 
  ylim(c(0,0.5)) +
  ggtitle('Fibroblast TF perturbability vs GM942 average specificity')
# plot(fibro_GTEx_indiv_nDEup_vs_JSsp_GM00942)

set.seed(4286)
fibro_GTEx_indiv_nDEup_vs_JSsp_GM00942_flip <- ggplot() + 
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "GM00942",
                      meanRPM >= 20) %>%
               inner_join(tfTab, by = 'gene_id'), 
             aes(GM00942_iCard_fibs_noHeart_noSkin, nDESeq2conditionsUp), alpha = 0.2) +
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "GM00942") %>%
               inner_join(hiFT_factorSum, by = "gene_id"), 
             aes(GM00942_iCard_fibs_noHeart_noSkin, nDESeq2conditionsUp), color = "forestgreen") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "GM00942") %>%
                    inner_join(hiFT_factorSum, by = "gene_id"), 
                  aes(GM00942_iCard_fibs_noHeart_noSkin, nDESeq2conditionsUp, label = GeneSymbol), color = "forestgreen") +
  theme_bw() + 
  # ylim(c(0,0.5)) +
  ggtitle('Fibroblast TF perturbability vs GM942 average specificity\nExcludes GTEx Heart and Skin\nAll tested factors')
# plot(fibro_GTEx_indiv_nDEup_vs_JSsp_GM00942_flip)
ggsave(fibro_GTEx_indiv_nDEup_vs_JSsp_GM00942_flip, file = paste0(graphSubdir, "/GTEx/fibro_GTEx_indiv_nDEup_vs_JSsp_GM00942_noSkinHeart_flip.pdf"), width = 8, height = 8, useDingbats = F)


set.seed(4286)
fibro_GTEx_indiv_nDEup_vs_JSsp_GM00942_flip_onlyBarriers <- ggplot() + 
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "GM00942",
                      meanRPM >= 20) %>%
               inner_join(tfTab, by = 'gene_id'), 
             aes(GM00942_iCard_fibs, nDESeq2conditionsUp), alpha = 0.2) +
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "GM00942") %>%
               inner_join(hiFT_factorSum, by = "gene_id") %>%
               filter(isBarrier == 1), 
             aes(GM00942_iCard_fibs, nDESeq2conditionsUp), color = "forestgreen") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "GM00942") %>%
                    inner_join(hiFT_factorSum, by = "gene_id") %>%
                    filter(isBarrier == 1), 
                  aes(GM00942_iCard_fibs, nDESeq2conditionsUp, label = GeneSymbol), color = "forestgreen") +
  theme_bw() + 
  # ylim(c(0,0.5)) +
  ggtitle('Fibroblast TF perturbability vs GM942 average specificity\nIncludes GTEx Heart and Skin\nOnly barriers')
# plot(fibro_GTEx_indiv_nDEup_vs_JSsp_GM00942_flip_onlyBarriers)
ggsave(fibro_GTEx_indiv_nDEup_vs_JSsp_GM00942_flip_onlyBarriers, file = paste0(graphSubdir, "/GTEx/fibro_GTEx_indiv_nDEup_vs_JSsp_GM00942_withSkinHeart_flip_onlyBarriers.pdf"), width = 8, height = 8, useDingbats = F)


set.seed(4286)
fibro_GTEx_indiv_nDEup_vs_JSsp_GM00942_noSkinAndHeart_flip_onlyBarriers <- ggplot() + 
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "GM00942",
                      meanRPM >= 20) %>%
               inner_join(tfTab, by = 'gene_id'), 
             aes(GM00942_iCard_fibs_noHeart_noSkin, nDESeq2conditionsUp), alpha = 0.2) +
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "GM00942") %>%
               inner_join(hiFT_factorSum, by = "gene_id") %>%
               filter(isBarrier == 1), 
             aes(GM00942_iCard_fibs_noHeart_noSkin, nDESeq2conditionsUp), color = "forestgreen") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "GM00942") %>%
                    inner_join(hiFT_factorSum, by = "gene_id") %>%
                    filter(isBarrier == 1), 
                  aes(GM00942_iCard_fibs_noHeart_noSkin, nDESeq2conditionsUp, label = GeneSymbol), color = "forestgreen") +
  theme_bw() + 
  # ylim(c(0,0.5)) +
  ggtitle('Fibroblast TF perturbability vs GM942 average specificity\nExcludes GTEx Heart and Skin\nOnly barriers')
# plot(fibro_GTEx_indiv_nDEup_vs_JSsp_GM00942_noSkinAndHeart_flip_onlyBarriers)
ggsave(fibro_GTEx_indiv_nDEup_vs_JSsp_GM00942_noSkinAndHeart_flip_onlyBarriers, file = paste0(graphSubdir, "/GTEx/fibro_GTEx_indiv_nDEup_vs_JSsp_GM00942_noSkinHeart_flip_onlyBarriers.pdf"), width = 8, height = 8, useDingbats = F)

set.seed(4286)
fibro_GTEx_indiv_nDEup_vs_JSsp_Skin_withFibsAndiCard_flip_onlyBarriers <- ggplot() + 
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "GM00942",
                      meanRPM >= 20) %>%
               inner_join(tfTab, by = 'gene_id'), 
             aes(Skin_iCard_fibs, nDESeq2conditionsUp), alpha = 0.2) +
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "GM00942") %>%
               inner_join(hiFT_factorSum, by = "gene_id") %>%
               filter(isBarrier == 1), 
             aes(Skin_iCard_fibs, nDESeq2conditionsUp), color = "forestgreen") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "GM00942") %>%
                    inner_join(hiFT_factorSum, by = "gene_id") %>%
                    filter(isBarrier == 1), 
                  aes(Skin_iCard_fibs, nDESeq2conditionsUp, label = GeneSymbol), color = "forestgreen") +
  theme_bw() + 
  # ylim(c(0,0.5)) +
  ggtitle('Fibroblast TF perturbability vs Skin average specificity\nIncludes GM942 fibroblasts and iCards\nOnly barriers')
# plot(fibro_GTEx_indiv_nDEup_vs_JSsp_Skin_withFibsAndiCard_flip_onlyBarriers)
ggsave(fibro_GTEx_indiv_nDEup_vs_JSsp_Skin_withFibsAndiCard_flip_onlyBarriers, file = paste0(graphSubdir, "/GTEx/fibro_GTEx_indiv_nDEup_vs_JSsp_Skin_withFibsAndiCard_flip_onlyBarriers.pdf"), width = 8, height = 8, useDingbats = F)

set.seed(4286)
fibro_GTEx_indiv_nDEup_vs_JSsp_Skin_noFibsAndiCard_flip_onlyBarriers <- ggplot() + 
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "GM00942",
                      meanRPM >= 20) %>%
               inner_join(tfTab, by = 'gene_id'), 
             aes(Skin, nDESeq2conditionsUp), alpha = 0.2) +
  geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "GM00942") %>%
               inner_join(hiFT_factorSum, by = "gene_id") %>%
               filter(isBarrier == 1), 
             aes(Skin, nDESeq2conditionsUp), color = "forestgreen") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "GM00942") %>%
                    inner_join(hiFT_factorSum, by = "gene_id") %>%
                    filter(isBarrier == 1), 
                  aes(Skin, nDESeq2conditionsUp, label = GeneSymbol), color = "forestgreen") +
  theme_bw() + 
  # ylim(c(0,0.5)) +
  ggtitle('Fibroblast TF perturbability vs Skin average specificity\nExcludes GM942 fibroblasts and iCards\nOnly barriers')
# plot(fibro_GTEx_indiv_nDEup_vs_JSsp_Skin_noFibsAndiCard_flip_onlyBarriers)
ggsave(fibro_GTEx_indiv_nDEup_vs_JSsp_Skin_noFibsAndiCard_flip_onlyBarriers, file = paste0(graphSubdir, "/GTEx/fibro_GTEx_indiv_nDEup_vs_JSsp_Skin_noFibsAndiCard_flip_onlyBarriers.pdf"), width = 8, height = 8, useDingbats = F)


fibro_TFs_JSsp942_hist_wMarkers_wRugs <- ggplot() +
  geom_histogram(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "GM00942",
                          meanRPM >= 20) %>%
                   inner_join(tfTab, by = 'gene_id') %>% mutate(graphType = 'histogram'), aes(GM00942), bins = 60) +
  geom_vline(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "GM00942") %>%
               inner_join(hiFT_factorSum, by = "gene_id") %>% mutate(graphType = 'histogram'), aes(xintercept = GM00942), color = "forestgreen") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "GM00942") %>%
                    inner_join(hiFT_factorSum, by = "gene_id") %>% mutate(graphType = 'histogram'), aes(x = GM00942, y = 200, label = GeneSymbol), color = "forestgreen") +
  geom_linerange(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "GM00942",
                          meanRPM >= 20) %>%
                   inner_join(tfTab, by = 'gene_id') %>% mutate(graphType = 'rug'), aes(GM00942, ymin = 0, ymax = 1)) +
  geom_linerange(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "GM00942") %>%
                   inner_join(hiFT_factorSum, by = "gene_id") %>% unique() %>% mutate(graphType = 'rug'), aes(GM00942, ymin = 1, ymax = 2), color = "forestgreen") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "GM00942") %>%
                    inner_join(hiFT_factorSum, by = "gene_id") %>% unique() %>% mutate(graphType = 'rug'), aes(x = GM00942, y = 1.75, label = GeneSymbol), color = "forestgreen") +
  theme_bw() + 
  facet_grid(graphType ~ ., scales = 'free_y') +
  # xlim(c(-0.5, 12.5)) +
  xlab("JSsp_GM00942 in controls") +
  ggtitle('JSsp_GM00942 of all TFs with mean RPM >20 in iCard controls\nall tested fibroblast marker genes in green')
# plot(fibro_TFs_JSsp942_hist_wMarkers_wRugs)


fibro_TFs_JSsp942_hist_wMarkers_wRugs_onlyBarriers_noSkin <- ggplot() +
  geom_histogram(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "GM00942",
                          meanRPM >= 20) %>%
                   inner_join(tfTab, by = 'gene_id') %>% mutate(graphType = 'histogram'), aes(GM00942_noSkin), bins = 60) +
  geom_vline(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "GM00942") %>%
               inner_join(hiFT_factorSum, by = "gene_id") %>% filter(isBarrier == 1) %>% mutate(graphType = 'histogram'), aes(xintercept = GM00942_noSkin), color = "forestgreen") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "GM00942") %>%
                    inner_join(hiFT_factorSum, by = "gene_id") %>% filter(isBarrier == 1) %>% mutate(graphType = 'histogram'), aes(x = GM00942_noSkin, y = 200, label = GeneSymbol), color = "forestgreen") +
  geom_linerange(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "GM00942",
                          meanRPM >= 20) %>%
                   inner_join(tfTab, by = 'gene_id') %>% mutate(graphType = 'rug'), aes(GM00942_noSkin, ymin = 0, ymax = 1)) +
  geom_linerange(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "GM00942") %>%
                   inner_join(hiFT_factorSum, by = "gene_id") %>% unique() %>% filter(isBarrier == 1) %>% mutate(graphType = 'rug'), aes(GM00942_noSkin, ymin = 1, ymax = 2), color = "forestgreen") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "GM00942") %>%
                    inner_join(hiFT_factorSum, by = "gene_id") %>% unique() %>% filter(isBarrier == 1) %>% mutate(graphType = 'rug'), aes(x = GM00942_noSkin, y = 1.75, label = GeneSymbol), color = "forestgreen") +
  theme_bw() + 
  facet_grid(graphType ~ ., scales = 'free_y') +
  # xlim(c(-0.5, 12.5)) +
  xlab("JSsp_GM00942 in controls") +
  ggtitle('JSsp_GM00942 of all TFs with mean RPM >20 in fibroblast controls\nfibroblast markers that are barriers to reprogramming in green\nScores calculated without GTEx Skin')
# plot(fibro_TFs_JSsp942_hist_wMarkers_wRugs_onlyBarriers_noSkin)

fibro_TFs_JSsp942_hist_wMarkers_wRugs_onlyBarriers_withSkinHeart <- ggplot() +
  geom_histogram(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "GM00942",
                          meanRPM >= 20) %>%
                   inner_join(tfTab, by = 'gene_id') %>% mutate(graphType = 'histogram'), aes(GM00942_iCard_fibs), bins = 60) +
  geom_vline(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "GM00942") %>%
               inner_join(hiFT_factorSum, by = "gene_id") %>% filter(isBarrier == 1) %>% mutate(graphType = 'histogram'), aes(xintercept = GM00942_iCard_fibs), color = "forestgreen") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "GM00942") %>%
                    inner_join(hiFT_factorSum, by = "gene_id") %>% filter(isBarrier == 1) %>% mutate(graphType = 'histogram'), aes(x = GM00942_iCard_fibs, y = 200, label = GeneSymbol), color = "forestgreen") +
  geom_linerange(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "GM00942",
                          meanRPM >= 20) %>%
                   inner_join(tfTab, by = 'gene_id') %>% mutate(graphType = 'rug'), aes(GM00942_iCard_fibs, ymin = 0, ymax = 1)) +
  geom_linerange(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "GM00942") %>%
                   inner_join(hiFT_factorSum, by = "gene_id") %>% unique() %>% filter(isBarrier == 1) %>% mutate(graphType = 'rug'), aes(GM00942_iCard_fibs, ymin = 1, ymax = 2), color = "forestgreen") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "GM00942") %>%
                    inner_join(hiFT_factorSum, by = "gene_id") %>% unique() %>% filter(isBarrier == 1) %>% mutate(graphType = 'rug'), aes(x = GM00942_iCard_fibs, y = 1.75, label = GeneSymbol), color = "forestgreen") +
  theme_bw() + 
  facet_grid(graphType ~ ., scales = 'free_y') +
  # xlim(c(-0.5, 12.5)) +
  xlab("JSsp_GM00942 in controls") +
  ggtitle('TF specificity score in: GM00942 fibroblasts\nDataset includes iCard-942, GTEx-Heart, and GTEx-Skin\nAll markers against all TFs with mean RPM >20 in iCard controls, barriers in green')
# plot(fibro_TFs_JSsp942_hist_wMarkers_wRugs_onlyBarriers_withSkinHeart)
ggsave(fibro_TFs_JSsp942_hist_wMarkers_wRugs_onlyBarriers_withSkinHeart, file = paste0(graphSubdir, "/GTEx/fibro_TFs_JSsp942_hist_wMarkers_wRugs_onlyBarriers_withSkinHeart.pdf"), width = 8, height = 8)

fibro_TFs_JSsp942_hist_wMarkers_wRugs_onlyBarriers_noSkin <- ggplot() +
  geom_histogram(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "GM00942",
                          meanRPM >= 20) %>%
                   inner_join(tfTab, by = 'gene_id') %>% mutate(graphType = 'histogram'), aes(GM00942_noSkin), bins = 60) +
  geom_vline(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "GM00942") %>%
               inner_join(hiFT_factorSum, by = "gene_id") %>% filter(isBarrier == 1) %>% mutate(graphType = 'histogram'), aes(xintercept = GM00942_noSkin), color = "forestgreen") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "GM00942") %>%
                    inner_join(hiFT_factorSum, by = "gene_id") %>% filter(isBarrier == 1) %>% mutate(graphType = 'histogram'), aes(x = GM00942_noSkin, y = 200, label = GeneSymbol), color = "forestgreen") +
  geom_linerange(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "GM00942",
                          meanRPM >= 20) %>%
                   inner_join(tfTab, by = 'gene_id') %>% mutate(graphType = 'rug'), aes(GM00942_noSkin, ymin = 0, ymax = 1)) +
  geom_linerange(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "GM00942") %>%
                   inner_join(hiFT_factorSum, by = "gene_id") %>% unique() %>% filter(isBarrier == 1) %>% mutate(graphType = 'rug'), aes(GM00942_noSkin, ymin = 1, ymax = 2), color = "forestgreen") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "GM00942") %>%
                    inner_join(hiFT_factorSum, by = "gene_id") %>% unique() %>% filter(isBarrier == 1) %>% mutate(graphType = 'rug'), aes(x = GM00942_noSkin, y = 1.75, label = GeneSymbol), color = "forestgreen") +
  theme_bw() + 
  facet_grid(graphType ~ ., scales = 'free_y') +
  # xlim(c(-0.5, 12.5)) +
  xlab("JSsp_GM00942 in controls") +
  ggtitle('TF specificity score in: GM00942 fibroblasts\nDataset excludes iCard-942 and GTEx-Skin\nAll markers against all TFs with mean RPM >20 in iCard controls, barriers in green')
# plot(fibro_TFs_JSsp942_hist_wMarkers_wRugs_onlyBarriers_noSkin)
ggsave(fibro_TFs_JSsp942_hist_wMarkers_wRugs_onlyBarriers_noSkin, file = paste0(graphSubdir, "/GTEx/fibro_TFs_JSsp942_hist_wMarkers_wRugs_onlyBarriers_noSkin.pdf"), width = 8, height = 8)


fibro_TFs_JSspSkin_hist_wMarkers_wRugs_onlyBarriers_iCardfibs <- ggplot() +
  geom_histogram(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "GM00942",
                          meanRPM >= 20) %>%
                   inner_join(tfTab, by = 'gene_id') %>% mutate(graphType = 'histogram'), aes(Skin_iCard_fibs), bins = 60) +
  geom_vline(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "GM00942") %>%
               inner_join(hiFT_factorSum, by = "gene_id") %>% filter(isBarrier == 1) %>% mutate(graphType = 'histogram'), aes(xintercept = Skin_iCard_fibs), color = "forestgreen") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "GM00942") %>%
                    inner_join(hiFT_factorSum, by = "gene_id") %>% filter(isBarrier == 1) %>% mutate(graphType = 'histogram'), aes(x = Skin_iCard_fibs, y = 200, label = GeneSymbol), color = "forestgreen") +
  geom_linerange(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "GM00942",
                          meanRPM >= 20) %>%
                   inner_join(tfTab, by = 'gene_id') %>% mutate(graphType = 'rug'), aes(Skin_iCard_fibs, ymin = 0, ymax = 1)) +
  geom_linerange(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "GM00942") %>%
                   inner_join(hiFT_factorSum, by = "gene_id") %>% unique() %>% filter(isBarrier == 1) %>% mutate(graphType = 'rug'), aes(Skin_iCard_fibs, ymin = 1, ymax = 2), color = "forestgreen") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "GM00942") %>%
                    inner_join(hiFT_factorSum, by = "gene_id") %>% unique() %>% filter(isBarrier == 1) %>% mutate(graphType = 'rug'), aes(x = Skin_iCard_fibs, y = 1.75, label = GeneSymbol), color = "forestgreen") +
  theme_bw() + 
  facet_grid(graphType ~ ., scales = 'free_y') +
  # xlim(c(-0.5, 12.5)) +
  xlab("JSsp_Skin in controls") +
  ggtitle('TF specificity score in: GTEx-Skin\nDataset includes iCard-942 and GM00942 fibroblasts\nAll markers against all TFs with mean RPM >20 in iCard controls, barriers in green')
# plot(fibro_TFs_JSspSkin_hist_wMarkers_wRugs_onlyBarriers_iCardfibs)
ggsave(fibro_TFs_JSspSkin_hist_wMarkers_wRugs_onlyBarriers_iCardfibs, file = paste0(graphSubdir, "/GTEx/fibro_TFs_JSspSkin_hist_wMarkers_wRugs_onlyBarriers_iCardfibs.pdf"), width = 8, height = 8)

fibro_TFs_JSspSkin_hist_wMarkers_wRugs_onlyBarriers_GTExOnly <- ggplot() +
  geom_histogram(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "GM00942",
                          meanRPM >= 20) %>%
                   inner_join(tfTab, by = 'gene_id') %>% mutate(graphType = 'histogram'), aes(Skin), bins = 60) +
  geom_vline(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "GM00942") %>%
               inner_join(hiFT_factorSum, by = "gene_id") %>% filter(isBarrier == 1) %>% mutate(graphType = 'histogram'), aes(xintercept = Skin), color = "forestgreen") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "GM00942") %>%
                    inner_join(hiFT_factorSum, by = "gene_id") %>% filter(isBarrier == 1) %>% mutate(graphType = 'histogram'), aes(x = Skin, y = 200, label = GeneSymbol), color = "forestgreen") +
  geom_linerange(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "GM00942",
                          meanRPM >= 20) %>%
                   inner_join(tfTab, by = 'gene_id') %>% mutate(graphType = 'rug'), aes(Skin, ymin = 0, ymax = 1)) +
  geom_linerange(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "GM00942") %>%
                   inner_join(hiFT_factorSum, by = "gene_id") %>% unique() %>% filter(isBarrier == 1) %>% mutate(graphType = 'rug'), aes(Skin, ymin = 1, ymax = 2), color = "forestgreen") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "GM00942") %>%
                    inner_join(hiFT_factorSum, by = "gene_id") %>% unique() %>% filter(isBarrier == 1) %>% mutate(graphType = 'rug'), aes(x = Skin, y = 1.75, label = GeneSymbol), color = "forestgreen") +
  theme_bw() + 
  facet_grid(graphType ~ ., scales = 'free_y') +
  # xlim(c(-0.5, 12.5)) +
  xlab("JSsp_Skin in controls") +
  ggtitle('TF specificity score in: GTEx-Skin\nDataset excludes iCard-942 and GM00942 fibroblasts\nAll markers against all TFs with mean RPM >20 in iCard controls, barriers in green')
# plot(fibro_TFs_JSspSkin_hist_wMarkers_wRugs_onlyBarriers_GTExOnly)
ggsave(fibro_TFs_JSspSkin_hist_wMarkers_wRugs_onlyBarriers_GTExOnly, file = paste0(graphSubdir, "/GTEx/fibro_TFs_JSspSkin_hist_wMarkers_wRugs_onlyBarriers_GTExOnly.pdf"), width = 8, height = 8)


fibro_TFs_nDEup_hist_wMarkers_wRugs_onlyBarriers <- ggplot() +
  geom_histogram(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "GM00942",
                          meanRPM >= 20) %>%
                   inner_join(tfTab, by = 'gene_id') %>% mutate(graphType = 'histogram'), aes(nDESeq2conditionsUp)) +
  geom_vline(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "GM00942") %>%
               inner_join(hiFT_factorSum, by = "gene_id") %>% filter(isBarrier == 1) %>% mutate(graphType = 'histogram'), aes(xintercept = nDESeq2conditionsUp), color = "forestgreen") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "GM00942") %>%
                    inner_join(hiFT_factorSum, by = "gene_id") %>% filter(isBarrier == 1) %>% mutate(graphType = 'histogram'), aes(x = nDESeq2conditionsUp, y = 200, label = GeneSymbol), color = "forestgreen") +
  geom_linerange(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "GM00942",
                          meanRPM >= 20) %>%
                   inner_join(tfTab, by = 'gene_id') %>% mutate(graphType = 'rug'), aes(nDESeq2conditionsUp, ymin = 0, ymax = 1)) +
  geom_linerange(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "GM00942") %>%
                   inner_join(hiFT_factorSum, by = "gene_id") %>% unique() %>% filter(isBarrier == 1) %>% mutate(graphType = 'rug'), aes(nDESeq2conditionsUp, ymin = 1, ymax = 2), color = "forestgreen") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "GM00942") %>%
                    inner_join(hiFT_factorSum, by = "gene_id") %>% unique() %>% filter(isBarrier == 1) %>% mutate(graphType = 'rug'), aes(x = nDESeq2conditionsUp, y = 1.75, label = GeneSymbol), color = "forestgreen") +
  theme_bw() + 
  facet_grid(graphType ~ ., scales = 'free_y') +
  # xlim(c(-0.5, 12.5)) +
  xlab("Num. conditions in which UP") +
  ggtitle('nDESeq2conditionsUp of all TFs with mean RPM >20 in iCard controls\nfibroblast markers that are barriers to reprogramming in green')
# plot(fibro_TFs_nDEup_hist_wMarkers_wRugs_onlyBarriers)
ggsave(fibro_TFs_nDEup_hist_wMarkers_wRugs_onlyBarriers, file = paste0(graphSubdir, "/GTEx/fibro_TFs_nDEup_hist_wMarkers_wRugs_onlyBarriers.pdf"), width = 8, height = 8)


## iCards
iCard_TFs_JSspiCard_hist_wMarkers_wRugs_onlyBarriers_withSkinHeart <- ggplot() +
  geom_histogram(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "iCard",
                          meanRPM >= 20) %>%
                   inner_join(tfTab, by = 'gene_id') %>% mutate(graphType = 'histogram'), aes(iCard_iCard_fibs), bins = 60) +
  geom_vline(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "iCard") %>%
               inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == 'iCard') %>% unique() %>% mutate(graphType = 'histogram'), aes(xintercept = iCard_iCard_fibs), color = "red") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "iCard") %>%
                    inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == 'iCard') %>% unique() %>% mutate(graphType = 'histogram'), aes(x = iCard_iCard_fibs, y = 200, label = GeneSymbol), color = "red") +
  geom_linerange(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "iCard",
                          meanRPM >= 20) %>%
                   inner_join(tfTab, by = 'gene_id') %>% mutate(graphType = 'rug'), aes(iCard_iCard_fibs, ymin = 0, ymax = 1)) +
  geom_linerange(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "iCard") %>%
                   inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == 'iCard') %>% unique() %>% mutate(graphType = 'rug'), aes(iCard_iCard_fibs, ymin = 1, ymax = 2), color = "red") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "iCard") %>%
                    inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == 'iCard') %>% unique() %>% mutate(graphType = 'rug'), aes(x = iCard_iCard_fibs, y = 1.75, label = GeneSymbol), color = "red") +
  theme_bw() + 
  facet_grid(graphType ~ ., scales = 'free_y') +
  # xlim(c(-0.5, 12.5)) +
  xlab("JSsp_iCard in controls") +
  ggtitle('TF specificity score in: iCards\nDataset includes GM00942 fibroblasts, GTEx-Heart, and GTEx-Skin\nAll markers against all TFs with mean RPM >20 in iCard controls, markers in red')
# plot(iCard_TFs_JSspiCard_hist_wMarkers_wRugs_onlyBarriers_withSkinHeart)
ggsave(iCard_TFs_JSspiCard_hist_wMarkers_wRugs_onlyBarriers_withSkinHeart, file = paste0(graphSubdir, "/GTEx/iCard_TFs_JSspiCard_hist_wMarkers_wRugs_onlyBarriers_withSkinHeart.pdf"), width = 8, height = 8)


iCard_TFs_JSspHeart_hist_wMarkers_wRugs_onlyBarriers_with942iCards <- ggplot() +
  geom_histogram(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "iCard",
                          meanRPM >= 20) %>%
                   inner_join(tfTab, by = 'gene_id') %>% mutate(graphType = 'histogram'), aes(Heart_iCard_fibs), bins = 60) +
  geom_vline(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "iCard") %>%
               inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == 'iCard') %>% unique() %>% mutate(graphType = 'histogram'), aes(xintercept = Heart_iCard_fibs), color = "red") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "iCard") %>%
                    inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == 'iCard') %>% unique() %>% mutate(graphType = 'histogram'), aes(x = Heart_iCard_fibs, y = 200, label = GeneSymbol), color = "red") +
  geom_linerange(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "iCard",
                          meanRPM >= 20) %>%
                   inner_join(tfTab, by = 'gene_id') %>% mutate(graphType = 'rug'), aes(Heart_iCard_fibs, ymin = 0, ymax = 1)) +
  geom_linerange(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "iCard") %>%
                   inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == 'iCard') %>% unique() %>% mutate(graphType = 'rug'), aes(Heart_iCard_fibs, ymin = 1, ymax = 2), color = "red") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "iCard") %>%
                    inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == 'iCard') %>% unique() %>% mutate(graphType = 'rug'), aes(x = Heart_iCard_fibs, y = 1.75, label = GeneSymbol), color = "red") +
  theme_bw() + 
  facet_grid(graphType ~ ., scales = 'free_y') +
  # xlim(c(-0.5, 12.5)) +
  xlab("JSsp_Heart in controls") +
  ggtitle('TF specificity score in: GTEx-Heart\nDataset includes GM00942 fibroblasts, iCard-942, and GTEx-Skin\nAll markers against all TFs with mean RPM >20 in iCard controls, markers in red')
# plot(iCard_TFs_JSspHeart_hist_wMarkers_wRugs_onlyBarriers_with942iCards)
ggsave(iCard_TFs_JSspHeart_hist_wMarkers_wRugs_onlyBarriers_with942iCards, file = paste0(graphSubdir, "/GTEx/iCard_TFs_JSspHeart_hist_wMarkers_wRugs_onlyBarriers_with942iCards.pdf"), width = 8, height = 8)

iCard_TFs_JSspHeart_hist_wMarkers_wRugs_onlyBarriers_justGTEx <- ggplot() +
  geom_histogram(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "iCard",
                          meanRPM >= 20) %>%
                   inner_join(tfTab, by = 'gene_id') %>% mutate(graphType = 'histogram'), aes(Heart), bins = 60) +
  geom_vline(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "iCard") %>%
               inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == 'iCard') %>% unique() %>% mutate(graphType = 'histogram'), aes(xintercept = Heart), color = "red") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "iCard") %>%
                    inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == 'iCard') %>% unique() %>% mutate(graphType = 'histogram'), aes(x = Heart, y = 200, label = GeneSymbol), color = "red") +
  geom_linerange(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "iCard",
                          meanRPM >= 20) %>%
                   inner_join(tfTab, by = 'gene_id') %>% mutate(graphType = 'rug'), aes(Heart, ymin = 0, ymax = 1)) +
  geom_linerange(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "iCard") %>%
                   inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == 'iCard') %>% unique() %>% mutate(graphType = 'rug'), aes(Heart, ymin = 1, ymax = 2), color = "red") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "iCard") %>%
                    inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == 'iCard') %>% unique() %>% mutate(graphType = 'rug'), aes(x = Heart, y = 1.75, label = GeneSymbol), color = "red") +
  theme_bw() + 
  facet_grid(graphType ~ ., scales = 'free_y') +
  # xlim(c(-0.5, 12.5)) +
  xlab("JSsp_Heart in controls") +
  ggtitle('TF specificity score in: GTEx-Heart\nDataset excludes GM00942 fibroblasts and iCard-942\nAll markers against all TFs with mean RPM >20 in iCard controls, markers in red')
# plot(iCard_TFs_JSspHeart_hist_wMarkers_wRugs_onlyBarriers_justGTEx)
ggsave(iCard_TFs_JSspHeart_hist_wMarkers_wRugs_onlyBarriers_justGTEx, file = paste0(graphSubdir, "/GTEx/iCard_TFs_JSspHeart_hist_wMarkers_wRugs_onlyBarriers_justGTEx.pdf"), width = 8, height = 8)

iCard_TFs_JSspiCard_hist_wMarkers_wRugs_onlyBarriers_noHeart <- ggplot() +
  geom_histogram(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "iCard",
                          meanRPM >= 20) %>%
                   inner_join(tfTab, by = 'gene_id') %>% mutate(graphType = 'histogram'), aes(iCard_iCard_noHeart), bins = 60) +
  geom_vline(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
               filter(cellType == "iCard") %>%
               inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == 'iCard') %>% unique() %>% mutate(graphType = 'histogram'), aes(xintercept = iCard_iCard_noHeart), color = "red") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "iCard") %>%
                    inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == 'iCard') %>% unique() %>% mutate(graphType = 'histogram'), aes(x = iCard_iCard_noHeart, y = 200, label = GeneSymbol), color = "red") +
  geom_linerange(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "iCard",
                          meanRPM >= 20) %>%
                   inner_join(tfTab, by = 'gene_id') %>% mutate(graphType = 'rug'), aes(iCard_iCard_noHeart, ymin = 0, ymax = 1)) +
  geom_linerange(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                   filter(cellType == "iCard") %>%
                   inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == 'iCard') %>% unique() %>% mutate(graphType = 'rug'), aes(iCard_iCard_noHeart, ymin = 1, ymax = 2), color = "red") +
  geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
                    filter(cellType == "iCard") %>%
                    inner_join(markerTab2, by = "gene_id") %>% filter(cellTypeMark == 'iCard') %>% unique() %>% mutate(graphType = 'rug'), aes(x = iCard_iCard_noHeart, y = 1.75, label = GeneSymbol), color = "red") +
  theme_bw() + 
  facet_grid(graphType ~ ., scales = 'free_y') +
  # xlim(c(-0.5, 12.5)) +
  xlab("JSsp_Heart in controls") +
  ggtitle('TF specificity score in: iCard\nDataset excludes GM00942 fibroblasts and GTEx-Heart\nAll markers against all TFs with mean RPM >20 in iCard controls, markers in red')
# plot(iCard_TFs_JSspiCard_hist_wMarkers_wRugs_onlyBarriers_noHeart)
ggsave(iCard_TFs_JSspiCard_hist_wMarkers_wRugs_onlyBarriers_noHeart, file = paste0(graphSubdir, "/GTEx/iCard_TFs_JSspiCard_hist_wMarkers_wRugs_onlyBarriers_noHeart.pdf"), width = 8, height = 8)


# for fig 1

mean_spec_for_hist_heart <- all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  mutate(logMeanTPM = log(meanTPM)) %>%
  filter(cellType == 'iCard', meanRPM > 20) %>%
  inner_join(tfTab, by = 'gene_id') %>%
  dplyr::select(gene_id, logMeanTPM, Heart) %>%
  gather('Measure', 'value', 2:3)

mean_spec_for_hist <- all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  mutate(logMeanTPM = log(meanTPM)) %>%
  filter(cellType == 'iCard', meanRPM > 20) %>%
  inner_join(tfTab, by = 'gene_id') %>%
  dplyr::select(gene_id, logMeanTPM, iCard_iCard_fibs_noHeart_noSkin) %>%
  gather('Measure', 'value', 2:3)
mean_spec_for_hist$Measure <- factor(mean_spec_for_hist$Measure)
mean_spec_for_hist$Measure <- factor(mean_spec_for_hist$Measure, levels = c('logMeanTPM', 'iCard_iCard_fibs_noHeart_noSkin'))


iCard_TFs_logMeanTPM_JSspHeart_hists_wMarkers_thin <- ggplot() +
  geom_histogram(data = mean_spec_for_hist_heart, aes(value), fill = 'grey', alpha = 0.9) +
  geom_linerange(data = mean_spec_for_hist_heart %>% inner_join(markerTab2, by = 'gene_id') %>% filter(cellTypeMark == 'iCard'), aes(value, ymin = 0, ymax = 40), size = 0.2, color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 0, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  # ylim(c(-100, 200)) +
  facet_grid(. ~ Measure, scales = 'free_x') #+
  # theme(axis.text = element_blank(),
  #       axis.title = element_blank())
ggsave(iCard_TFs_logMeanTPM_JSspHeart_hists_wMarkers_thin, file = paste0(graphSubdir,'/differentialExpression/allSamples/iCards/iCard_TFs_logMeanTPM_JSspHeart_hists_wMarkers_thin.pdf'), width = 12, height = 3, useDingbats = F)

iCard_TFs_logMeanTPM_JSspiCard_noHeartSkin_hists_wMarkers_thin <- ggplot() +
  geom_histogram(data = mean_spec_for_hist, aes(value), fill = 'grey', alpha = 0.9) +
  geom_linerange(data = mean_spec_for_hist %>% inner_join(markerTab2, by = 'gene_id') %>% filter(cellTypeMark == 'iCard'), aes(value, ymin = 0, ymax = 40), size = 0.2, color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 0, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  # ylim(c(-100, 200)) +
  facet_grid(. ~ Measure, scales = 'free_x') #+
# theme(axis.text = element_blank(),
#       axis.title = element_blank())
ggsave(iCard_TFs_logMeanTPM_JSspiCard_noHeartSkin_hists_wMarkers_thin, file = paste0(graphSubdir,'/differentialExpression/allSamples/iCards/iCard_TFs_logMeanTPM_JSspiCard_noHeartSkin_hists_wMarkers_thin.pdf'), width = 12, height = 3, useDingbats = F)
  
iCard_TFs_logMeanTPM_JSspiCard_noHeartSkin_hists_wMarkers_thin_nolabs <- ggplot() +
  geom_histogram(data = mean_spec_for_hist, aes(value), fill = 'grey', alpha = 0.9) +
  geom_linerange(data = mean_spec_for_hist %>% inner_join(markerTab2, by = 'gene_id') %>% filter(cellTypeMark == 'iCard'), aes(value, ymin = 0, ymax = 40), size = 0.2, color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 0, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  # ylim(c(-100, 200)) +
  facet_grid(. ~ Measure, scales = 'free_x') +
theme(axis.text = element_blank(),
      axis.title = element_blank())
ggsave(iCard_TFs_logMeanTPM_JSspiCard_noHeartSkin_hists_wMarkers_thin_nolabs, file = paste0(graphSubdir,'/differentialExpression/allSamples/iCards/iCard_TFs_logMeanTPM_JSspiCard_noHeartSkin_hists_wMarkers_thin_nolabs.pdf'), width = 4, height = 1, useDingbats = F)

iCard_TFs_logMeanTPM_JSspiCard_noHeartSkin_hists_wMarkers_wGeneNames_thin_nolabs <- ggplot() +
  geom_histogram(data = mean_spec_for_hist, aes(value), fill = 'grey', alpha = 0.9) +
  geom_linerange(data = mean_spec_for_hist %>% inner_join(markerTab2, by = 'gene_id') %>% filter(cellTypeMark == 'iCard'), aes(value, ymin = 0, ymax = 40), size = 0.2, color = "red") +
  geom_text_repel(data = mean_spec_for_hist %>% inner_join(markerTab2, by = 'gene_id') %>% filter(cellTypeMark == 'iCard'), aes(x = value, y = 40, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  # ylim(c(-100, 200)) +
  facet_grid(. ~ Measure, scales = 'free_x') +
  theme(axis.text = element_blank(),
        axis.title = element_blank())
ggsave(iCard_TFs_logMeanTPM_JSspiCard_noHeartSkin_hists_wMarkers_wGeneNames_thin_nolabs, file = paste0(graphSubdir,'/differentialExpression/allSamples/iCards/iCard_TFs_logMeanTPM_JSspiCard_noHeartSkin_hists_wMarkers_wGeneNames_thin_nolabs.pdf'), width = 4, height = 1, useDingbats = F)


iCard_TFs_JSiCard_hist_wMarkerLinesegs_wJitter_wLabs_thin_wide <- ggplot() +
  geom_histogram(data = mean_spec_for_hist %>% filter(Measure == 'iCard_iCard_fibs_noHeart_noSkin'), aes(value), fill = 'grey', alpha = 0.9) +
  geom_linerange(data = mean_spec_for_hist %>% filter(Measure == 'iCard_iCard_fibs_noHeart_noSkin') %>% inner_join(markerTab2, by = 'gene_id') %>% filter(cellTypeMark == 'iCard') %>% unique(), aes(value, ymin = 0, ymax = 40), size = 0.2, color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 0, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() #+ 
  # ylim(c(-100, 200)) +
  # facet_grid(. ~ Measure, scales = 'free_x') +
  # theme(axis.text = element_blank(),
  #       axis.title = element_blank())
ggsave(iCard_TFs_JSiCard_hist_wMarkerLinesegs_wJitter_wLabs_thin_wide, file = paste0(graphSubdir,'/differentialExpression/allSamples/iCards/iCard_TFs_JSiCard_hist_wMarkerLinesegs_wJitter_wLabs_thin_wide.pdf'), width = 3.2, height = 1.68, useDingbats = F)

iCard_TFs_JSiCard_hist_wMarkerLinesegs_wJitter_noLabs_thin_wide <- ggplot() +
  geom_histogram(data = mean_spec_for_hist %>% filter(Measure == 'iCard_iCard_fibs_noHeart_noSkin'), aes(value), fill = 'grey', alpha = 0.9) +
  geom_linerange(data = mean_spec_for_hist %>% filter(Measure == 'iCard_iCard_fibs_noHeart_noSkin') %>% inner_join(markerTab2, by = 'gene_id') %>% filter(cellTypeMark == 'iCard') %>% unique(), aes(value, ymin = 0, ymax = 40), size = 0.2, color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 0, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  # ylim(c(-100, 200)) +
  # facet_grid(. ~ Measure, scales = 'free_x') +
  theme(axis.text = element_blank(),
        axis.title = element_blank())
ggsave(iCard_TFs_JSiCard_hist_wMarkerLinesegs_wJitter_noLabs_thin_wide, file = paste0(graphSubdir,'/differentialExpression/allSamples/iCards/iCard_TFs_JSiCard_hist_wMarkerLinesegs_wJitter_noLabs_thin_wide.pdf'), width = 3.2, height = 0.68, useDingbats = F)

