#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = T)

# check if there are exactly 3 arguments. If not, return an error
if (length(args) < 3 | length(args) > 3) {
  stop("Exactly three arguments must be supplied: projectDir, procDataSubdir, and graphSubdir.", call.=FALSE)
} 
if (length(args)==3){
  projectDir = args[1]
  procDataSubdir = args[2]
  graphSubdir = args[3]
}

# heatmap draft

library(DESeq2)
library(pheatmap)
library(tidyr)
library(dplyr)
library(magrittr)

# if running manually in Rstudio, start here and use:
# projectDir = '~/Dropbox (RajLab)/Shared_IanM/cellid_201807_onward/'
# procDataSubdir = 'procDataScripted'
# graphSubdir = 'graphs'

procDataDir = paste0(projectDir, procDataSubdir)
graphDir = paste0(projectDir, graphSubdir)

setwd(projectDir)
# load TPM, DE, and perturbability data
all_tpms_geneFilt_readFilt_manualFilt <- readRDS(paste0(procDataDir, "/allExperiments/allTPMs_geneFilt_readFilt_manualFilt/se_allTPMs_geneFilt_readFilt_manualFilt.rds"))
# all_tpms_geneFilt_readFilt_manualFilt <- readRDS(paste0(procDataDir, "/allExperiments/allTPMs_geneFilt_readFilt_manualFilt/se_allTPMs_geneFilt_readFilt_manualFilt_rerun.rds"))

tfGeneIDfile = "annotations/TF_gene_ids.csv"
tfTab = as_tibble(read.csv(tfGeneIDfile, stringsAsFactors = F))
colnames(tfTab) = "gene_id"

iCard_DE_results_sig <- as_tibble(read.table(paste0(graphDir, '/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_DESeqResults.txt'), header = T, stringsAsFactors = F)) %>%
  filter(padj <= 0.1)

fibro_DE_results_sig <- as_tibble(read.table(paste0(graphDir, '/differentialExpression/allSamples/fibroblasts/fibro_readFilt_manualFilt_DESeqResults.txt'), header = T, stringsAsFactors = F)) %>%
  filter(padj <= 0.1)

DE_genes <- unique(c(iCard_DE_results_sig$gene_id, fibro_DE_results_sig$gene_id))
indiv_perturbability = as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_variability_metrics.txt"), header = T, stringsAsFactors = F, sep = "\t"))

# filter to DE genes only with mean control TPM > 5 (or minimum cutoff in full_DESeq2_perturbability.R)
all_tpms_geneFilt_readFilt_manualFilt_DE <- all_tpms_geneFilt_readFilt_manualFilt[rownames(assay(all_tpms_geneFilt_readFilt_manualFilt)) %in% DE_genes,]


#### All genes (not just TFs)
## Different variability cutoffs for inclusion in heatmap (memory limits) in 

# Top 500 fano factors of raw TPM values
geneFanos = rowVars(as.matrix(assay(all_tpms_geneFilt_readFilt_manualFilt_DE)))/rowMeans(as.matrix(assay(all_tpms_geneFilt_readFilt_manualFilt_DE)))

all_tpms_geneFilt_readFilt_manualFilt_DE_top500fano = all_tpms_geneFilt_readFilt_manualFilt_DE[order(-geneFanos)[1:500],]

sampleData = colData(all_tpms_geneFilt_readFilt_manualFilt_DE_top500fano)
sampleDataP = sampleData[,c("cellType", "condition", "experiment")]

heatmapDir <- paste0(graphDir, '/heatmaps')
if (!dir.exists(heatmapDir)) {
  dir.create(heatmapDir)
}

# pdf(file = paste0(heatmapDir, "/fibsAndiCards_top500fano_minTPM5_anyDEgene_TPMval.pdf"), width = 8, height = 20)
# pheatmap(
#   mat = t(scale(t(as.matrix(assay(all_tpms_geneFilt_readFilt_manualFilt_DE_top500fano))))),
#   border_color = NA,
#   show_rownames = F,
#   show_colnames = F,
#   annotation_col = as.data.frame(sampleDataP),
#   drop_levels = T
# )
# dev.off()
# 
# # Top 500 fano factors of log2(TPM) values
# log2GeneFanos <- rowVars(log2(as.matrix(assay(all_tpms_geneFilt_readFilt_manualFilt_DE))))/rowMeans(log2(as.matrix(assay(all_tpms_geneFilt_readFilt_manualFilt_DE))))
# all_tpms_geneFilt_readFilt_manualFilt_DE_top500log2fano = all_tpms_geneFilt_readFilt_manualFilt_DE[order(-log2GeneFanos)[1:500],]
# 
# sampleData = colData(all_tpms_geneFilt_readFilt_manualFilt_DE_top500log2fano)
# sampleDataP = sampleData[,c("cellType", "condition", "experiment")]
# 
# pdf(file = paste0(heatmapDir, "/fibsAndiCards_top500fano_minTPM5_anyDEgene_log2TPMval.pdf"), width = 8, height = 20)
# pheatmap(
#   mat = t(scale(t(log2(as.matrix(assay(all_tpms_geneFilt_readFilt_manualFilt_DE_top500log2fano)))))),
#   border_color = NA,
#   show_rownames = F,
#   show_colnames = F,
#   annotation_col = as.data.frame(sampleDataP),
#   drop_levels = T
# )
# dev.off()
# 
# # Top 500 CV of log2(TPM) values
# log2GeneCVs <- rowSds(log2(as.matrix(assay(all_tpms_geneFilt_readFilt_manualFilt_DE))))/rowMeans(log2(as.matrix(assay(all_tpms_geneFilt_readFilt_manualFilt_DE))))
# all_tpms_geneFilt_readFilt_manualFilt_DE_top500log2CV = all_tpms_geneFilt_readFilt_manualFilt_DE[order(-log2GeneCVs)[1:500],]
# 
# sampleData = colData(all_tpms_geneFilt_readFilt_manualFilt_DE_top500log2CV)
# sampleDataP = sampleData[,c("cellType", "condition", "experiment")]
# 
# pdf(file = paste0(heatmapDir, "/fibsAndiCards_top500CV_minTPM5_anyDEgene_log2TPMval.pdf"), width = 8, height = 20)
# pheatmap(
#   mat = t(scale(t(log2(as.matrix(assay(all_tpms_geneFilt_readFilt_manualFilt_DE_top500log2CV)))))),
#   border_color = NA,
#   show_rownames = F,
#   show_colnames = F,
#   annotation_col = as.data.frame(sampleDataP),
#   drop_levels = T
# )
# dev.off()

## TFs only
all_tpms_geneFilt_readFilt_manualFilt_DE_TFsOnly <- all_tpms_geneFilt_readFilt_manualFilt_DE[rownames(assay(all_tpms_geneFilt_readFilt_manualFilt_DE)) %in% tfTab$gene_id,]
log2GeneCVs_TFonly <- rowSds(log2(as.matrix(assay(all_tpms_geneFilt_readFilt_manualFilt_DE))))/rowMeans(log2(as.matrix(assay(all_tpms_geneFilt_readFilt_manualFilt_DE))))

# Try all TFs (only 784 passing TPM > 5 cutoff and DE)
sampleData = colData(all_tpms_geneFilt_readFilt_manualFilt_DE_TFsOnly)
sampleDataP = sampleData[,c("cellType", "condition", "experiment")]

# scaled log2 TPM values
# pdf(file = paste0(heatmapDir, "/fibsAndiCards_TFonly_minTPM5_DE_log2TPMval.pdf"), width = 8, height = 25)
# pheatmap(
#   mat = t(scale(t(log2(as.matrix(assay(all_tpms_geneFilt_readFilt_manualFilt_DE_TFsOnly)) + 0.01)))),
#   border_color = NA,
#   show_rownames = F,
#   show_colnames = F,
#   annotation_col = as.data.frame(sampleDataP),
#   drop_levels = T
# )
# dev.off()

# TFs with meanRPM > 20 in both cell types (570 total)
TFs_min20_bothTypes <- indiv_perturbability %>% 
  inner_join(tfTab) %>%
  dplyr::select(gene_id, cellType, meanRPM) %>%
  spread(cellType, meanRPM) %>%
  filter(iCard >= 20 & GM00942 >= 20)
all_tpms_geneFilt_readFilt_manualFilt_DE_TFsOnly_minRPM20both <- all_tpms_geneFilt_readFilt_manualFilt_DE_TFsOnly[rownames(assay(all_tpms_geneFilt_readFilt_manualFilt_DE_TFsOnly)) %in% TFs_min20_bothTypes$gene_id,]
sampleData = colData(all_tpms_geneFilt_readFilt_manualFilt_DE_TFsOnly_minRPM20both)
sampleDataP = sampleData[,c("cellType", "condition", "experiment")]

# log2 TPM values
# pdf(file = paste0(heatmapDir, "/fibsAndiCards_TFonly_minTPM5_DE_log2TPMval_minRPM20both.pdf"), width = 8, height = 25)
# pheatmap(
#   mat = t(scale(t(log2(as.matrix(assay(all_tpms_geneFilt_readFilt_manualFilt_DE_TFsOnly_minRPM20both)) + 0.01)))),
#   border_color = NA,
#   show_rownames = F,
#   show_colnames = F,
#   annotation_col = as.data.frame(sampleDataP),
#   drop_levels = T
# )
# dev.off()

# log2 TPM values
# pdf(file = paste0(heatmapDir, "/fibsAndiCards_TFonly_minTPM5_DE_log2TPMval_unscaled_minRPM20both.pdf"), width = 8, height = 25)
# pheatmap(
#   mat = log2(as.matrix(assay(all_tpms_geneFilt_readFilt_manualFilt_DE_TFsOnly_minRPM20both)) + 0.01),
#   border_color = NA,
#   show_rownames = F,
#   show_colnames = F,
#   annotation_col = as.data.frame(sampleDataP),
#   drop_levels = T
# )
# dev.off()

sampleDataP = sampleData[,c("cellType", "condition")]
pdf(file = paste0(heatmapDir, "/fibsAndiCards_TFonly_minTPM5_DE_log2TPMval_unscaled_minRPM20both_cellTypeDrug_rerun.pdf"), width = 8, height = 25)
pheatmap(
  mat = log2(as.matrix(assay(all_tpms_geneFilt_readFilt_manualFilt_DE_TFsOnly_minRPM20both)) + 0.01),
  border_color = NA,
  show_rownames = F,
  show_colnames = F,
  annotation_col = as.data.frame(sampleDataP),
  drop_levels = T
)
dev.off()

sampleDataP = data.frame(cellType = sampleData[,c("cellType")])
rownames(sampleDataP) = rownames(sampleData)
pdf(file = paste0(heatmapDir, "/fibsAndiCards_TFonly_minTPM5_DE_log2TPMval_unscaled_minRPM20both_cellType_rerun.pdf"), width = 8, height = 25)
pheatmap(
  mat = log2(as.matrix(assay(all_tpms_geneFilt_readFilt_manualFilt_DE_TFsOnly_minRPM20both)) + 0.01),
  border_color = NA,
  show_rownames = F,
  show_colnames = F,
  annotation_col = as.data.frame(sampleDataP),
  drop_levels = T
)
dev.off()

# # TFs with meanRPM > 20 in EITHER cell types (870 total)
# TFs_min20_eitherType <- indiv_perturbability %>% 
#   inner_join(tfTab) %>%
#   dplyr::select(gene_id, cellType, meanRPM) %>%
#   spread(cellType, meanRPM) %>%
#   filter(iCard >= 20 | GM00942 >= 20)
# all_tpms_geneFilt_readFilt_manualFilt_DE_TFsOnly_minRPM20either <- all_tpms_geneFilt_readFilt_manualFilt_DE_TFsOnly[rownames(assay(all_tpms_geneFilt_readFilt_manualFilt_DE_TFsOnly)) %in% TFs_min20_eitherType$gene_id,]
# sampleData = colData(all_tpms_geneFilt_readFilt_manualFilt_DE_TFsOnly_minRPM20either)
# sampleDataP = sampleData[,c("cellType", "condition", "experiment")]
# 
# # scaled log2 TPM values
# pdf(file = paste0(heatmapDir, "/fibsAndiCards_TFonly_minTPM5_DE_log2TPMval_scaled_minRPM20eitherType.pdf"), width = 8, height = 25)
# pheatmap(
#   mat = t(scale(t(log2(as.matrix(assay(all_tpms_geneFilt_readFilt_manualFilt_DE_TFsOnly_minRPM20either)) + 0.01)))),
#   border_color = NA,
#   show_rownames = F,
#   show_colnames = F,
#   annotation_col = as.data.frame(sampleDataP),
#   drop_levels = T
# )
# dev.off()
# 
# sampleDataP = sampleData[,c("cellType", "condition")]
# # scaled log2 TPM values
# pdf(file = paste0(heatmapDir, "/fibsAndiCards_TFonly_minTPM5_DE_log2TPMval_scaled_minRPM20eitherType_cellTypeDrug.pdf"), width = 8, height = 25)
# pheatmap(
#   mat = t(scale(t(log2(as.matrix(assay(all_tpms_geneFilt_readFilt_manualFilt_DE_TFsOnly_minRPM20either)) + 0.01)))),
#   border_color = NA,
#   show_rownames = F,
#   show_colnames = F,
#   annotation_col = as.data.frame(sampleDataP),
#   drop_levels = T
# )
# dev.off()
# 
# # log2 TPM values
# pdf(file = paste0(heatmapDir, "/fibsAndiCards_TFonly_minTPM5_DE_log2TPMval_unscaled_minRPM20eitherType.pdf"), width = 8, height = 25)
# pheatmap(
#   mat = log2(as.matrix(assay(all_tpms_geneFilt_readFilt_manualFilt_DE_TFsOnly_minRPM20either)) + 0.01),
#   border_color = NA,
#   show_rownames = F,
#   show_colnames = F,
#   annotation_col = as.data.frame(sampleDataP),
#   drop_levels = T
# )
# dev.off()
# 
# ### pick genes that are most variable within iCards and most variable within fibroblasts, then plot together
# iCard_tpms_geneFilt_readFilt_manualFilt_DE <- all_tpms_geneFilt_readFilt_manualFilt_DE[,colData(all_tpms_geneFilt_readFilt_manualFilt_DE)$cellType == 'iCard']
# log2GeneCVs_iCard <- rowSds(log2(as.matrix(assay(iCard_tpms_geneFilt_readFilt_manualFilt_DE))))/rowMeans(log2(as.matrix(assay(iCard_tpms_geneFilt_readFilt_manualFilt_DE))))
# 
# fibro_tpms_geneFilt_readFilt_manualFilt_DE <- all_tpms_geneFilt_readFilt_manualFilt_DE[,colData(all_tpms_geneFilt_readFilt_manualFilt_DE)$cellType == 'GM00942']
# log2GeneCVs_fibro <- rowSds(log2(as.matrix(assay(fibro_tpms_geneFilt_readFilt_manualFilt_DE))))/rowMeans(log2(as.matrix(assay(fibro_tpms_geneFilt_readFilt_manualFilt_DE))))
# 
# 
# #### iCards only
# iCard_tpms_geneFilt_readFilt_manualFilt_DE <- all_tpms_geneFilt_readFilt_manualFilt_DE[,colData(all_tpms_geneFilt_readFilt_manualFilt_DE)$cellType == 'iCard']
# log2GeneCVs_iCard <- rowSds(log2(as.matrix(assay(iCard_tpms_geneFilt_readFilt_manualFilt_DE))))/rowMeans(log2(as.matrix(assay(iCard_tpms_geneFilt_readFilt_manualFilt_DE))))
# 
# iCard_tpms_geneFilt_readFilt_manualFilt_DE_top500log2CV = iCard_tpms_geneFilt_readFilt_manualFilt_DE[order(-log2GeneCVs)[1:500],]
# 
# sampleData = colData(iCard_tpms_geneFilt_readFilt_manualFilt_DE_top500log2CV)
# sampleDataP = sampleData[,c("cellType", "condition", "experiment")]
# 
# pdf(file = paste0(heatmapDir, "/iCards_top500CV_minTPM5_anyDEgene_log2TPMval.pdf"), width = 8, height = 20)
# pheatmap(
#   mat = t(scale(t(log2(as.matrix(assay(iCard_tpms_geneFilt_readFilt_manualFilt_DE_top500log2CV)))))),
#   border_color = NA,
#   show_rownames = F,
#   show_colnames = F,
#   annotation_col = as.data.frame(sampleDataP),
#   drop_levels = T
# )
# dev.off()
# 
# # min RPM = 20
# iCard_minRPM20 <- indiv_perturbability %>%
#   dplyr::select(gene_id, cellType, meanRPM) %>%
#   spread(cellType, meanRPM) %>%
#   filter(iCard >= 20)
# fibro_minRPM20 <- indiv_perturbability %>%
#   dplyr::select(gene_id, cellType, meanRPM) %>%
#   spread(cellType, meanRPM) %>%
#   filter(GM00942 >= 20)
# 
# log2GeneCVs_iCard <- rowSds(log2(as.matrix(assay(iCard_tpms_geneFilt_readFilt_manualFilt_DE)) + 0.01))/rowMeans(log2(as.matrix(assay(iCard_tpms_geneFilt_readFilt_manualFilt_DE)) + 0.01))
# names(log2GeneCVs_iCard) <- rownames(assay(iCard_tpms_geneFilt_readFilt_manualFilt_DE))
# 
# log2GeneCVs_fibro <- rowSds(log2(as.matrix(assay(fibro_tpms_geneFilt_readFilt_manualFilt_DE)) + 0.01))/rowMeans(log2(as.matrix(assay(fibro_tpms_geneFilt_readFilt_manualFilt_DE)) + 0.01))
# names(log2GeneCVs_fibro) <- rownames(assay(fibro_tpms_geneFilt_readFilt_manualFilt_DE))
# 
# log2GeneCVs_both <- rowSds(log2(as.matrix(assay(all_tpms_geneFilt_readFilt_manualFilt_DE)) + 0.01))/rowMeans(log2(as.matrix(assay(all_tpms_geneFilt_readFilt_manualFilt_DE)) + 0.01))
# names(log2GeneCVs_both) <- rownames(assay(all_tpms_geneFilt_readFilt_manualFilt_DE))
# 
# top500log2CV_iCard <- names(log2GeneCVs_iCard)[order(-log2GeneCVs_iCard)[1:500]]
# top500log2CV_fibro <- names(log2GeneCVs_fibro)[order(-log2GeneCVs_fibro)[1:500]]
# top500log2CV_both <- names(log2GeneCVs_both)[order(-log2GeneCVs_both)[1:500]]
# top500log2CV_each <- unique(c(top500log2CV_iCard, top500log2CV_fibro, top500log2CV_both))
# 
# all_tpms_geneFilt_readFilt_manualFilt_DE_topsEach <- all_tpms_geneFilt_readFilt_manualFilt_DE[rownames(assay(all_tpms_geneFilt_readFilt_manualFilt_DE)) %in% top500log2CV_each,]
# sampleData = colData(all_tpms_geneFilt_readFilt_manualFilt_DE_topsEach)
# sampleDataP = sampleData[,c("cellType", "condition", "experiment")]
# 
# # log2 TPM values
# pdf(file = paste0(heatmapDir, "/fibsAndiCards_TFonly_minTPM5_DE_log2TPMval_toplog2CVsEach.pdf"), width = 8, height = 25)
# pheatmap(
#   mat = t(scale(t(log2(as.matrix(assay(all_tpms_geneFilt_readFilt_manualFilt_DE_topsEach)) + 0.01)))),
#   border_color = NA,
#   show_rownames = F,
#   show_colnames = F,
#   annotation_col = as.data.frame(sampleDataP),
#   drop_levels = T
# )
# dev.off()
# 
# 
# pdf(file = paste0(heatmapDir, "/iCards_top500CV_minTPM5_anyDEgene_log2TPMval.pdf"), width = 8, height = 20)
# pheatmap(
#   mat = t(scale(t(log2(as.matrix(assay(iCard_tpms_geneFilt_readFilt_manualFilt_DE_top500log2CV_minRPM20)))))),
#   border_color = NA,
#   show_rownames = F,
#   show_colnames = F,
#   annotation_col = as.data.frame(sampleDataP),
#   drop_levels = T
# )
# dev.off()
# 
# 
# ## Use DESeq2 counts
# # from http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#heatmap-of-the-count-matrix
# dds_fibro <- readRDS(paste0(graphDir, '/differentialExpression/allSamples/fibroblasts/fibro_dds_se.rds'))
# dds_iCard <- readRDS(paste0(graphDir, '/differentialExpression/allSamples/iCards/iCards_dds_se.rds'))
# 
# #vst
# dds_fibro_vst <- vst(dds_fibro)
# dds_iCard_vst <- vst(dds_iCard)
# 
# dds_fibro_vst <- vst(dds_fibro)
# dds_iCard_vst <- vst(dds_iCard)
# 
# # filter by minimum RPM and TPM < 10000
# 
# iCard_meanRPM <- indiv_perturbability %>% filter(cellType == 'iCard', meanRPM >= 20, meanTPM < 10000) %>% dplyr::select(cellType, gene_id, meanRPM, meanTPM, nDESeq2conditionsAll, nDESeq2conditionsUp, nDESeq2conditionsDown)
# iCard_meanRPM[is.na(iCard_meanRPM)] <- 0
# 
# fibro_meanRPM <- indiv_perturbability %>% filter(cellType == 'GM00942', meanRPM >= 20, meanRPM < 5000) %>% dplyr::select(cellType, gene_id, meanRPM, meanTPM, nDESeq2conditionsAll, nDESeq2conditionsUp, nDESeq2conditionsDown)
# fibro_meanRPM[is.na(fibro_meanRPM)] <- 0
# 
# fibro_meanRPM_tf <- inner_join(fibro_meanRPM, tfTab)
# 
# dds_fibro_vst_min20rpm <- assay(dds_fibro_vst)[fibro_meanRPM$gene_id,]
# 
# meanSdPlot(assay(dds_fibro_vst)[fibro_meanRPM$gene_id,])
# 
# # fibroblast
# 
# select_fibro_vst_rpmFilt_mean <- order(rowMeans(dds_fibro_vst_min20rpm),
#                          decreasing=TRUE)[1:500]
# 
# pdf(file = paste0(heatmapDir, "/fibro_DESeq2_vst_top500fibroMean_scaled.pdf"), width = 8, height = 20)
# pheatmap(
#   mat = as.matrix(dds_fibro_vst_min20rpm[select_fibro_vst_rpmFilt_mean,]),
#   scale = 'row',
#   border_color = NA,
#   show_rownames = F,
#   show_colnames = F,
#   annotation_col = as.data.frame(colData(dds_fibro_vst)[,c("experiment","cellType","condition")]),
#   drop_levels = T
# )
# dev.off()
# 
# select_fibro_vst_rpmFilt_CV <- order(rowSds(dds_fibro_vst_min20rpm)/rowMeans(dds_fibro_vst_min20rpm),
#                                        decreasing=TRUE)[1:500]
# 
# pdf(file = paste0(heatmapDir, "/fibro_DESeq2_vst_top500fibroCV_scaled.pdf"), width = 8, height = 20)
# pheatmap(
#   mat = as.matrix(dds_fibro_vst_min20rpm[select_fibro_vst_rpmFilt_CV,]),
#   scale = 'row',
#   border_color = NA,
#   show_rownames = F,
#   show_colnames = F,
#   annotation_col = as.data.frame(colData(dds_fibro_vst)[,c("experiment","cellType","condition")]),
#   drop_levels = T
# )
# dev.off()
# 
# pdf(file = paste0(heatmapDir, "/fibro_DESeq2_vst_top100fibroCV_scaled.pdf"), width = 8, height = 20)
# pheatmap(
#   mat = as.matrix(dds_fibro_vst_min20rpm[select_fibro_vst_rpmFilt_CV[1:100],]),
#   scale = 'row',
#   border_color = NA,
#   show_rownames = F,
#   show_colnames = F,
#   annotation_col = as.data.frame(colData(dds_fibro_vst)[,c("experiment","cellType","condition")]),
#   drop_levels = T
# )
# dev.off()
# 
# 
# pdf(file = paste0(heatmapDir, "/fibro_DESeq2_vst_meanRPMfilt_TFs_scaled.pdf"), width = 8, height = 20)
# pheatmap(
#   mat = as.matrix(dds_fibro_vst_min20rpm[fibro_meanRPM_tf$gene_id,]),
#   scale = 'row',
#   border_color = NA,
#   show_rownames = F,
#   show_colnames = F,
#   annotation_col = as.data.frame(colData(dds_fibro_vst)[,c("experiment","cellType","condition")]),
#   drop_levels = T
# )
# dev.off()
# 
# pdf(file = paste0(heatmapDir, "/fibro_DESeq2_vst_meanRPMfilt_TFs_unscaled.pdf"), width = 8, height = 20)
# pheatmap(
#   mat = as.matrix(dds_fibro_vst_min20rpm[fibro_meanRPM_tf$gene_id,]),
#   # scale = 'row',
#   border_color = NA,
#   show_rownames = F,
#   show_colnames = F,
#   annotation_col = as.data.frame(colData(dds_fibro_vst)[,c("experiment","cellType","condition")]),
#   drop_levels = T
# )
# dev.off()
# 
# pdf(file = paste0(heatmapDir, "/fibro_DESeq2_vst_meanRPMfilt_TFs_unscaled_nosuperhigh.pdf"), width = 8, height = 20)
# pheatmap(
#   mat = as.matrix(dds_fibro_vst_min20rpm[fibro_meanRPM_tf$gene_id,]),
#   # scale = 'row',
#   border_color = NA,
#   show_rownames = F,
#   show_colnames = F,
#   annotation_col = as.data.frame(colData(dds_fibro_vst)[,c("experiment","cellType","condition")]),
#   drop_levels = T
# )
# dev.off()
# 
# 
# 
# # iCard
# select_iCard_mean <- order(rowMeans(counts(dds_iCard,normalized=TRUE)),
#                 decreasing=TRUE)[1:500]
# 
# select_iCard_CV <- order(rowSds(counts(dds_iCard,normalized=TRUE))/rowMeans(counts(dds_iCard,normalized=TRUE)),
#                       decreasing=TRUE)[1:500]
# 
# 
# pdf(file = paste0(heatmapDir, "/iCards_DESeq2_vst_top500iCardCV_unscaled.pdf"), width = 8, height = 20)
# pheatmap(
#   mat = as.matrix(assay(dds_iCard_vst)[select_iCard_CV,]),
#   border_color = NA,
#   show_rownames = F,
#   show_colnames = F,
#   annotation_col = as.data.frame(colData(dds_iCard_vst)[,c("experiment","cellType","condition")]),
#   drop_levels = T
# )
# dev.off()
# 
# 
# df <- as.data.frame(colData(dds)[,c("condition","type")])
# pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#          cluster_cols=T, annotation_col=df)
# 
# 
# # combine iCard and fibroblast
# 
# genes_in_both <- rownames(assay(dds_fibro_vst))[rownames(assay(dds_fibro_vst)) %in% rownames(assay(dds_iCard_vst))]
# 
# vst_for_plot <- cbind(assay(dds_fibro_vst)[genes_in_both,], assay(dds_iCard_vst)[genes_in_both,])
# 
# sampleDat_vst <- rbind(colData(dds_fibro_vst), colData(dds_iCard_vst))
# sampleDat_vstP = sampleDat_vst[,c("cellType", "condition", "experiment")]
# 
# top500vstCV_total <- rowSds(vst_for_plot)/rowMeans(vst_for_plot)
# names(top500vstCV_total) <- rownames(vst_for_plot)
# 
# # rank(top500vstCV_total)
# 
# pdf(file = paste0(heatmapDir, "/iCards_DESeq2_vst_top500totalCV.pdf"), width = 8, height = 20)
# pheatmap(
#   mat = as.matrix(t(scale(t(vst_for_plot[rank(top500vstCV_total) <= 500,])))),
#   border_color = NA,
#   show_rownames = F,
#   show_colnames = F,
#   annotation_col = as.data.frame(sampleDat_vstP),
#   drop_levels = T
# )
# dev.off()
# 
# # rank(top500rowMean)
# 
# pdf(file = paste0(heatmapDir, "/iCards_DESeq2_vst_top500totalMean.pdf"), width = 8, height = 20)
# pheatmap(
#   mat = as.matrix(t(scale(t(vst_for_plot[rank(rowMeans(vst_for_plot)) <= 500,])))),
#   border_color = NA,
#   show_rownames = F,
#   show_colnames = F,
#   annotation_col = as.data.frame(sampleDat_vstP),
#   drop_levels = T
# )
# dev.off()


# #plot
# select <- order(rowMeans(counts(dds,normalized=TRUE)),
#                 decreasing=TRUE)[1:20]
# df <- as.data.frame(colData(dds)[,c("condition","type")])
# pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#          cluster_cols=T, annotation_col=df)

