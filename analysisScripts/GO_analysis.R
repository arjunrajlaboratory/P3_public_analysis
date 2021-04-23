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

# projectDir = '~/Dropbox (RajLab)/Shared_IanM/cellid_201807_onward/'
# procDataSubdir = 'procDataScripted/'
# graphSubdir = 'graphs'

procDataDir = paste0(projectDir, procDataSubdir)
graphDir = paste0(projectDir, graphSubdir)

library(clusterProfiler)
library(org.Hs.eg.db)

if(!dir.exists(paste0(graphDir, '/GSEA/'))){
  dir.create(paste0(graphDir, '/GSEA/'))
}

tfGeneIDfile = "annotations/TF_gene_ids.csv"
tfTab = as_tibble(read.csv(paste0(projectDir, tfGeneIDfile), stringsAsFactors = F))
colnames(tfTab) = "gene_id"

indiv_perturbability = as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_variability_metrics.txt"), header = T, stringsAsFactors = F, sep = "\t"))
indiv_perturbability[is.na(indiv_perturbability)] <- 0

iCard_highFreqUP_indivSamps <- indiv_perturbability %>%
  filter(cellType == "iCard", nDESeq2conditionsUp >= 4, meanRPM > 20) %>% unique()

iCard_highFreqDown_indivSamps <- indiv_perturbability %>%
  filter(cellType == "iCard", nDESeq2conditionsDown >= 4, meanRPM > 20) %>% unique()

iCard_lowFreq_indivSamps <- indiv_perturbability %>%
  filter(cellType == "iCard", nDESeq2conditionsAll <= 1, meanRPM > 20) %>% unique()

iCard_all_indivSamps <- indiv_perturbability %>%
  filter(cellType == "iCard", meanRPM > 20) %>% unique()

iCard_highFreqUP_df = bitr(iCard_highFreqUP_indivSamps$gene_id,
                                             fromType = "ENSEMBL",
                                             toType = c("SYMBOL", "ENTREZID"),
                                             OrgDb = org.Hs.eg.db)

iCard_highFreqDOWN_df = bitr(iCard_highFreqDown_indivSamps$gene_id,
                           fromType = "ENSEMBL",
                           toType = c("SYMBOL", "ENTREZID"),
                           OrgDb = org.Hs.eg.db)

iCard_lowFreq_df = bitr(iCard_lowFreq_indivSamps$gene_id,
                             fromType = "ENSEMBL",
                             toType = c("SYMBOL", "ENTREZID"),
                             OrgDb = org.Hs.eg.db)

iCard_all_indivSamps_df = bitr(iCard_all_indivSamps$gene_id,
                               fromType = "ENSEMBL",
                               toType = c("SYMBOL", "ENTREZID"),
                               OrgDb = org.Hs.eg.db)

iCard_go_cc_overrep_tfs <- enrichGO(gene = iCard_highFreqUP_df$ENTREZID,
                                 universe      = iCard_all_indivSamps_df$ENTREZID,
                                 OrgDb         = org.Hs.eg.db,
                                 # keytype       = "ENSEMBL",
                                 ont           = "CC",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 1,
                                 qvalueCutoff  = 1,
                                 readable      = TRUE)
# head(summary(iCard_go_cc_overrep_tfs))
# dotplot(iCard_go_cc_overrep_tfs, showCategory=25)

iCard_go_cc_overrep_tfs_sig <- enrichGO(gene = iCard_highFreqUP_df$ENTREZID,
                                    universe      = iCard_all_indivSamps_df$ENTREZID,
                                    OrgDb         = org.Hs.eg.db,
                                    # keytype       = "ENSEMBL",
                                    ont           = "CC",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = 0.01,
                                    qvalueCutoff  = 0.05,
                                    readable      = TRUE)
# head(summary(iCard_go_cc_overrep_tfs_sig))

pdf(paste0(graphDir, '/GSEA/iCard_highUp_allGenes_DotPlot.pdf'), width = 9, height = 12, useDingbats = F)
dotplot(iCard_go_cc_overrep_tfs_sig, showCategory=25)
dev.off()

iCard_go_bp_overrep_tfs_sig <- enrichGO(gene = iCard_highFreqUP_df$ENTREZID,
                                        universe      = iCard_all_indivSamps_df$ENTREZID,
                                        OrgDb         = org.Hs.eg.db,
                                        # keytype       = "ENSEMBL",
                                        ont           = "BP",
                                        pAdjustMethod = "BH",
                                        pvalueCutoff  = 0.01,
                                        qvalueCutoff  = 0.05,
                                        readable      = TRUE)
# head(summary(iCard_go_bp_overrep_tfs_sig))
pdf(paste0(graphDir, '/GSEA/iCard_highUp_allGenes_DotPlot_GO-BP.pdf'), width = 9, height = 12, useDingbats = F)
dotplot(iCard_go_bp_overrep_tfs_sig, showCategory=25)
dev.off()

pdf(paste0(graphDir, '/GSEA/iCard_highUp_allGenes_goPlot_GO-BP.pdf'), width = 25, height = 12, useDingbats = F)
goplot(iCard_go_bp_overrep_tfs_sig, showCategory=25)
dev.off()

pdf(paste0(graphDir, '/GSEA/iCard_highUp_allGenes_DotPlot_GO-CC.pdf'), width = 9, height = 12, useDingbats = F)
dotplot(iCard_go_cc_overrep_tfs_sig, showCategory=25)
dev.off()

pdf(paste0(graphDir, '/GSEA/iCard_highUp_allGenes_GOplot.pdf'), width = 9, height = 9, useDingbats = F)
goplot(iCard_go_cc_overrep_tfs_sig) + ggtitle('GO Cellular Component overenrichement analysis\nFrequently up-regulated genes')
dev.off()

iCard_go_cc_overrepDown_tfs_sig <- enrichGO(gene = iCard_highFreqDOWN_df$ENTREZID,
                                        universe      = iCard_all_indivSamps_df$ENTREZID,
                                        OrgDb         = org.Hs.eg.db,
                                        # keytype       = "ENSEMBL",
                                        ont           = "CC",
                                        pAdjustMethod = "BH",
                                        pvalueCutoff  = 0.01,
                                        qvalueCutoff  = 0.05,
                                        readable      = TRUE)
# head(summary(iCard_go_cc_overrepDown_tfs_sig))
# dotplot(iCard_go_cc_overrepDown_tfs_sig, showCategory=25)

iCard_go_cc_overrepLowFreq_tfs_relax <- enrichGO(gene = iCard_lowFreq_df$ENTREZID,
                                               universe      = iCard_all_indivSamps_df$ENTREZID,
                                               OrgDb         = org.Hs.eg.db,
                                               # keytype       = "ENSEMBL",
                                               ont           = "CC",
                                               pAdjustMethod = "BH",
                                               pvalueCutoff  = 0.1,
                                               qvalueCutoff  = 0.1,
                                               readable      = TRUE)
# head(summary(iCard_go_cc_overrepLowFreq_tfs_relax))
# dotplot(iCard_go_cc_overrepLowFreq_tfs_relax, showCategory=25)

iCard_go_cc_overrepLowFreq_tfs_sig <- enrichGO(gene = iCard_lowFreq_df$ENTREZID,
                                            universe      = iCard_all_indivSamps_df$ENTREZID,
                                            OrgDb         = org.Hs.eg.db,
                                            # keytype       = "ENSEMBL",
                                            ont           = "CC",
                                            pAdjustMethod = "BH",
                                            pvalueCutoff  = 0.01,
                                            qvalueCutoff  = 0.05,
                                            readable      = TRUE)
# head(summary(iCard_go_cc_overrepLowFreq_tfs_sig))
# dotplot(iCard_go_cc_overrepLowFreq_tfs_sig, showCategory=25)
pdf(paste0(graphDir, '/GSEA/iCard_lowFreq_allGenes_GOplot.pdf'), width = 9, height = 9, useDingbats = F)
goplot(iCard_go_cc_overrepLowFreq_tfs_relax) + ggtitle('GO Cellular Component overenrichement analysis\nNon-perturbable genes (<=1 DE)')
dev.off()


# fibroblasts

fibro_highFreqUP_indivSamps <- indiv_perturbability %>%
  filter(cellType == "GM00942", nDESeq2conditionsUp >= 4, meanRPM > 20) %>% unique()

fibro_highFreqDown_indivSamps <- indiv_perturbability %>%
  filter(cellType == "GM00942", nDESeq2conditionsDown >= 4, meanRPM > 20) %>% unique()

fibro_highFreqAll_indivSamps <- indiv_perturbability %>%
  filter(cellType == "GM00942", nDESeq2conditionsAll >= 6, meanRPM > 20) %>% unique()

fibro_lowFreq_indivSamps <- indiv_perturbability %>%
  filter(cellType == "GM00942", nDESeq2conditionsAll <= 1, meanRPM > 20) %>% unique()

fibro_all_indivSamps <- indiv_perturbability %>%
  filter(cellType == "GM00942", meanRPM > 20) %>% unique()

fibro_highFreqUP_df = bitr(fibro_highFreqUP_indivSamps$gene_id,
                           fromType = "ENSEMBL",
                           toType = c("SYMBOL", "ENTREZID"),
                           OrgDb = org.Hs.eg.db)

fibro_highFreqDOWN_df = bitr(fibro_highFreqDown_indivSamps$gene_id,
                             fromType = "ENSEMBL",
                             toType = c("SYMBOL", "ENTREZID"),
                             OrgDb = org.Hs.eg.db)

fibro_highFreqALL_df = bitr(fibro_highFreqAll_indivSamps$gene_id,
                             fromType = "ENSEMBL",
                             toType = c("SYMBOL", "ENTREZID"),
                             OrgDb = org.Hs.eg.db)

fibro_lowFreq_df = bitr(fibro_lowFreq_indivSamps$gene_id,
                        fromType = "ENSEMBL",
                        toType = c("SYMBOL", "ENTREZID"),
                        OrgDb = org.Hs.eg.db)

fibro_all_indivSamps_df = bitr(fibro_all_indivSamps$gene_id,
                               fromType = "ENSEMBL",
                               toType = c("SYMBOL", "ENTREZID"),
                               OrgDb = org.Hs.eg.db)


fibro_go_cc_overrepLowFreq_tfs_relax <- enrichGO(gene = fibro_lowFreq_df$ENTREZID,
                                                 universe      = fibro_all_indivSamps_df$ENTREZID,
                                                 OrgDb         = org.Hs.eg.db,
                                                 # keytype       = "ENSEMBL",
                                                 ont           = "CC",
                                                 pAdjustMethod = "BH",
                                                 pvalueCutoff  = 0.1,
                                                 qvalueCutoff  = 0.1,
                                                 readable      = TRUE)
# head(summary(fibro_go_cc_overrepLowFreq_tfs_relax))

pdf(paste0(graphDir, '/GSEA/fibro_lowFreq_allGenes_Dotplot.pdf'), width = 9, height = 9, useDingbats = F)
dotplot(fibro_go_cc_overrepLowFreq_tfs_relax, showCategory=25)
dev.off()
pdf(paste0(graphDir, '/GSEA/fibro_lowFreq_allGenes_GOplot.pdf'), width = 9, height = 9, useDingbats = F)
goplot(fibro_go_cc_overrepLowFreq_tfs_relax)
dev.off()

fibro_go_cc_overrep_tfs_sig <- enrichGO(gene = fibro_highFreqUP_df$ENTREZID,
                                        universe      = fibro_all_indivSamps_df$ENTREZID,
                                        OrgDb         = org.Hs.eg.db,
                                        # keytype       = "ENSEMBL",
                                        ont           = "CC",
                                        pAdjustMethod = "BH",
                                        pvalueCutoff  = 0.01,
                                        qvalueCutoff  = 0.05,
                                        readable      = TRUE)
# head(summary(fibro_go_cc_overrep_tfs_sig))
pdf(paste0(graphDir, '/GSEA/fibro_highFreqUP_allGenes_Dotplot.pdf'), width = 9, height = 9, useDingbats = F)
dotplot(fibro_go_cc_overrep_tfs_sig, showCategory=25)
dev.off()
pdf(paste0(graphDir, '/GSEA/fibro_highFreqUP_allGenes_GOplot.pdf'), width = 9, height = 9, useDingbats = F)
goplot(fibro_go_cc_overrep_tfs_sig)
dev.off()

fibro_go_cc_overrepDEall_tfs_sig <- enrichGO(gene = fibro_highFreqALL_df$ENTREZID,
                                        universe      = fibro_all_indivSamps_df$ENTREZID,
                                        OrgDb         = org.Hs.eg.db,
                                        # keytype       = "ENSEMBL",
                                        ont           = "CC",
                                        pAdjustMethod = "BH",
                                        pvalueCutoff  = 0.01,
                                        qvalueCutoff  = 0.05,
                                        readable      = TRUE)
# head(summary(fibro_go_cc_overrepDEall_tfs_sig))
# dotplot(fibro_go_cc_overrepDEall_tfs_sig, showCategory=25)

fibro_go_bp_overrepDEall_tfs_sig <- enrichGO(gene = fibro_highFreqALL_df$ENTREZID,
                                             universe      = fibro_all_indivSamps_df$ENTREZID,
                                             OrgDb         = org.Hs.eg.db,
                                             # keytype       = "ENSEMBL",
                                             ont           = "BP",
                                             pAdjustMethod = "BH",
                                             pvalueCutoff  = 0.01,
                                             qvalueCutoff  = 0.05,
                                             readable      = TRUE)
# head(summary(fibro_go_bp_overrepDEall_tfs_sig))
pdf(paste0(graphDir, '/GSEA/fibro_highFreqALL_allGenes_Dotplot.pdf'), width = 9, height = 9, useDingbats = F)
dotplot(fibro_go_bp_overrepDEall_tfs_sig, showCategory=25)
dev.off()
pdf(paste0(graphDir, '/GSEA/fibro_highFreqALL_allGenes_GOplot.pdf'), width = 9, height = 9, useDingbats = F)
goplot(fibro_go_bp_overrepDEall_tfs_sig)
dev.off()


fibro_go_bp_overrepLowFreq_tfs_sig <- enrichGO(gene = fibro_lowFreq_df$ENTREZID,
                                             universe      = fibro_all_indivSamps_df$ENTREZID,
                                             OrgDb         = org.Hs.eg.db,
                                             # keytype       = "ENSEMBL",
                                             ont           = "BP",
                                             pAdjustMethod = "BH",
                                             pvalueCutoff  = 0.01,
                                             qvalueCutoff  = 0.05,
                                             readable      = TRUE)
# head(summary(fibro_go_bp_overrepLowFreq_tfs_sig))
# dotplot(fibro_go_bp_overrepLowFreq_tfs_sig, showCategory=25)
# goplot(fibro_go_bp_overrepLowFreq_tfs_sig)
