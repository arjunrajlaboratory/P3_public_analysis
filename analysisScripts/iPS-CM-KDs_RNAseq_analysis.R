#!/usr/bin/env Rscript
# don't run lines 3-13 if running from R workspace
args = commandArgs(trailingOnly = T)

# check if there are exactly 4 arguments. If not, return an error
if (length(args) < 4 | length(args) > 4) {
  stop("Exactly four arguments must be supplied: projectDir, procDataSubdir, graphSubdir, and runDEseq.", call.=FALSE)
} 
if (length(args)==4){
  projectDir = args[1]
  procDataSubdir = args[2]
  graphSubdir = args[3]
  runDEseq = as.logical(args[4]) # TRUE OR FALSE
}

# projectDir = '~/Dropbox (RajLab)/Shared_IanM/cellid_201807_onward/'
# procDataSubdir = 'procDataScripted_test429'
# graphSubdir = 'graphs'

procDataDir = paste0(projectDir, procDataSubdir)
graphDir = paste0(projectDir, graphSubdir)

library(tidyverse)
library(magrittr)
library(ggrepel)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)

setwd(projectDir)
tfGeneIDfile = "annotations/TF_gene_ids.csv"
tfTab = as_tibble(read.csv(tfGeneIDfile))
colnames(tfTab) = "gene_id"

geneIDlengthsFile = "annotations/hg19_lengths_per_gene_id.tsv"
geneIDfile = "annotations/hg19gene_idToGeneSymbol.tsv"
geneIDtoGeneSymbol <- as_tibble(read.csv(geneIDfile, sep='\t', stringsAsFactors = F))
geneLengths <- as_tibble(read.table(geneIDlengthsFile, sep='\t', stringsAsFactors = F, header = T))

counts_tab <- as_tibble(read.table(paste0(projectDir, 'extractedData/RNAseq/iPS-CM-KDs/iPSCMKDs1_meltedData.tsv'), sep = '\t', header = T, stringsAsFactors = F)) %>%
  separate(sampleID, into = c('target', 'col1', 'col2', 'col3', 'rest'), sep = "-", remove = F) %>%
  mutate(target = ifelse(target == 'NKX2', 'NKX2_5', target),
         repID = ifelse(target == 'NKX2_5', col3, col2)) %>%
  dplyr::select(-c(col1,col2,col3))

# filter out unmapped reads
counts_tab %<>%
  filter(gene_id != "__alignment_not_unique" & gene_id != "__ambiguous" & gene_id != "__no_feature" & gene_id != "__not_aligned" & gene_id != "__too_low_aQual")

# join with lengths per gene_id
counts_tab %<>% inner_join(geneLengths, by = 'gene_id')

tpms_tab <- counts_tab %>%
  group_by(target, repID) %>%
  mutate(rpk = 1000*counts/length, rpkScalePerMillion = sum(rpk)/1000000, tpm = rpk/rpkScalePerMillion) %>%
  dplyr::select(-c(rpk, rpkScalePerMillion))

avg_tpms_tab <- tpms_tab %>%
  group_by(target, gene_id) %>%
  summarise(mean_TPM = mean(tpm),
            nSamp = length(tpm),
            sem_TPM = sd(tpm)/sqrt(nSamp),
            mean_counts = mean(counts),
            sem_counts = sd(counts)/sqrt(nSamp),
            length = unique(length))

nosi_avg_tpms <- avg_tpms_tab %>%
  filter(target == 'NOSI') %>%
  dplyr::rename(nosi_mean_TPM = mean_TPM,
                nosi_sem_TPM = sem_TPM)


## DESeq2
#reformat
library(DESeq2)

count_mat <- counts_tab %>%
  dplyr::select(sampleID, gene_id, counts) %>%
  spread(sampleID, counts) %>%
  filter(!(gene_id %in% c("__alignment_not_unique", "__ambiguous", "__no_feature", "__not_aligned", "__too_low_aQual")))
count_mat <- as.data.frame(count_mat)
rownames(count_mat) <- count_mat$gene_id
count_mat$gene_id <- NULL
count_mat <- as.matrix(count_mat)

col_dat <- counts_tab %>%
  dplyr::select(experiment, sampleID, target, repID) %>%
  unique() 
col_dat <- as.data.frame(col_dat)
rownames(col_dat) <- col_dat$sampleID
col_dat$target<-factor(col_dat$target)
col_dat$target<-factor(col_dat$target, levels = c('SCR', levels(col_dat$target)[1:7], levels(col_dat$target)[9:11]))

row_dat <- left_join(geneLengths, geneIDtoGeneSymbol, by = 'gene_id')
row_dat <- as.data.frame(row_dat)
rownames(row_dat) <- row_dat$gene_id

KD_SE <- SummarizedExperiment(assays = list(counts= as.matrix(count_mat)), colData = col_dat, rowData = row_dat)

if(!runDEseq){
  print('Skipping DESeq2 step, using pre-run results.')
  allResults <- read.table(paste0(graphDir, '/iPS-CM-KDs/differentialExpression/DE_results.txt'), header = T, stringsAsFactors = F)
  ddsSE <- DESeqDataSet(KD_SE, design = ~ target)
  keep <- rowSums(counts(ddsSE)) >= 10
  ddsSE <- ddsSE[keep,]
  }
if(runDEseq){
  ddsSE <- DESeqDataSet(KD_SE, design = ~ target)
  
  keep <- rowSums(counts(ddsSE)) >= 10
  ddsSE <- ddsSE[keep,]
  
  dds <- DESeq(ddsSE)
  
  if (!dir.exists(paste0(graphDir, '/iPS-CM-KDs/'))) {
    dir.create(paste0(graphDir, '/iPS-CM-KDs/'))
  }
  if (!dir.exists(paste0(graphDir, '/iPS-CM-KDs/differentialExpression'))) {
    dir.create(paste0(graphDir, '/iPS-CM-KDs/differentialExpression'))
  }
  iPSCMKDs_MAfile = paste0(graphDir, '/iPS-CM-KDs/differentialExpression/MA_plots_vsSCR.pdf')
  pdf(iPSCMKDs_MAfile, width = 16, height = 12)
  par(mfrow = c(3,4))
  deTab = list()
  allResults = list()
  for (cond in unique(colData(ddsSE)$target)[unique(colData(ddsSE)$target) != "SCR"]) {
    
    resultsTemp = results(dds, contrast = c("target", as.character(cond), "SCR"))
    plotMA(resultsTemp, ylim=c(-2,2), main = paste0(cond, "_vs_SCRcontrol"))
    
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
    resultsTemp$target = cond
    
    if(is.null(dim(allResults))){
      allResults = as_tibble(resultsTemp)
    } else {
      allResults = bind_rows(allResults, as_tibble(resultsTemp))
    } 
    
    print(paste0("Finished ", cond, ". nDE = ", as.character(nDE)))
  }
  dev.off()
  write.table(allResults,file = paste0(graphDir, '/iPS-CM-KDs/differentialExpression/DE_results.txt'), sep = '\t', quote = F, row.names = F)
}
targets <- unique(counts_tab$target)[!(unique(counts_tab$target) %in% c('SCR', 'NOSI'))]
targetsGS <- targets
targetsGS[targetsGS == 'NKX2_5'] <- 'NKX2-5'

DESeq2forDisplay <- allResults %>%
  inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
  filter(target %in% targets,
         GeneSymbol %in% targetsGS) %>%
  mutate(percKD = 100*(1-2^log2FoldChange),
         percKDupper = 100*(1-2^(log2FoldChange - lfcSE)),
         percKDlower = 100*(1-2^(log2FoldChange + lfcSE)))

DE_KD_alltargets_perKD <- ggplot() +
  geom_bar(data = DESeq2forDisplay, aes(GeneSymbol, -percKD), stat = 'identity', alpha = 0.5, width = 0.8) +
  geom_errorbar(data = DESeq2forDisplay, aes(GeneSymbol, ymax = -percKDupper, ymin = -percKDlower), width = 0.25) +
  facet_grid(target ~ .) +
  theme_classic() +
  xlab('Gene Symbol (expression)') +
  ylab('% knockdown (by DESeq2)') +
  ggtitle('Gene expression after target knockdown\nFacet labels = siRNA target\nKnockdown percentage vs. scrambled controls')
ggsave(DE_KD_alltargets_perKD, file = paste0(graphDir, '/iPS-CM-KDs/differentialExpression/PercentKD_vsSCR.pdf'), width = 6, height = 10)  

if(!dir.exists(paste0(graphDir, '/iPS-CM-KDs/GOanalysis/'))) {
  dir.create(paste0(graphDir, '/iPS-CM-KDs/GOanalysis/'))
}

targs = unique(as.character(colData(ddsSE)$target))

lfcthresh = 0.5
pthresh = 0.05

up_GO_res <- list()
down_GO_res <- list()
# simcut = 0.5
for (targ in targs[targs != 'SCR']){
  
  cat(paste0('Working on ', targ, '\n'))
  
  temp_UP <- allResults %>%
    filter(target == targ,
           log2FoldChange >= lfcthresh,
           padj < pthresh)
  
  temp_DOWN <- allResults %>%
    filter(target == targ,
           log2FoldChange <= -lfcthresh,
           padj < pthresh)
  
  temp_up_df = bitr(temp_UP$gene_id,
                     fromType = "ENSEMBL",
                     toType = c("SYMBOL", "ENTREZID"),
                     OrgDb = org.Hs.eg.db)
  
  temp_down_df = bitr(temp_DOWN$gene_id,
                     fromType = "ENSEMBL",
                     toType = c("SYMBOL", "ENTREZID"),
                     OrgDb = org.Hs.eg.db)
  
  up_GO_res[[targ]] <- simplify(enrichGO(gene = temp_up_df$ENTREZID,
                                         OrgDb         = org.Hs.eg.db,
                                         # keytype       = "ENSEMBL",
                                         ont           = "CC",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05,
                                         qvalueCutoff  = 0.05,
                                         readable      = TRUE))#, cutoff = simcut)
  down_GO_res[[targ]] <- simplify(enrichGO(gene = temp_down_df$ENTREZID,
                                           OrgDb         = org.Hs.eg.db,
                                           # keytype       = "ENSEMBL",
                                           ont           = "CC",
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05,
                                           qvalueCutoff  = 0.05,
                                           readable      = TRUE))
}

lv = 4#GO filter level
for (targ in targs[targs != 'SCR']){
  
  up_res <- gofilter(up_GO_res[[targ]], level = lv)@result
  height_up = length(unique(up_res$Description))
  up_res$Description <- factor(up_res$Description)
  up_res$Description <- factor(up_res$Description, levels = up_res$Description[order(-up_res$p.adjust)])
  
  up_pval_plot <- ggplot(up_res, aes(Description, -log10(p.adjust))) +
    geom_bar(stat = 'identity') +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
    coord_flip() +
    theme_classic() +
    theme(axis.text = element_text(size = 6))
  
  dn_res <- gofilter(down_GO_res[[targ]], level = lv)@result
  height_dn = length(unique(dn_res$Description))
  dn_res$Description <- factor(dn_res$Description)
  dn_res$Description <- factor(dn_res$Description, levels = dn_res$Description[order(-dn_res$p.adjust)])
  
  dn_pval_plot <- ggplot(dn_res, aes(Description, -log10(p.adjust))) +
    geom_bar(stat = 'identity') +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
    coord_flip() +
    theme_classic()+
    theme(axis.text = element_text(size = 6))
  
  ggsave(up_pval_plot, file = paste0(graphDir, '/iPS-CM-KDs/GOanalysis/', targ, '_upRegulated_GO_pval.pdf'), units='in', height = 0.5+0.1*height_up)#, width = 1.75)#, height = 0.5+0.1*height_up)
  # ggsave(dn_pval_plot, file = paste0(graphDir, '/iPS-CM-KDs/GOanalysis/', targ, '_downRegulated_GO_pval.pdf'), units='in', height = 0.5+0.1*height_dn)#, width = 1.75)#, height = 0.5+0.1*height_dn)
  write.table(up_res, file = paste0(graphDir, '/iPS-CM-KDs/GOanalysis/', targ, '_upRegulated_GO_results.txt'), quote = F, sep = '\t')
  write.table(dn_res, file = paste0(graphDir, '/iPS-CM-KDs/GOanalysis/', targ, '_downRegulated_GO_results.txt'), quote = F, sep = '\t')
  
}

# heatmap view

lv = 4 #GO filter level
all_up_res<-list()
all_dn_res<-list()
up_terms <- tibble(
  GO = character(),
  Description = character()
)
dn_terms <- tibble(
  GO = character(),
  Description = character()
)
up_nores <- character()
dn_nores <- character()
KDtargs <- c('HMGB2', 'MEF2C', 'NKX2_5', 'HAND1', 'SOX11', 'IRX4', 'SKIDA1', 'ZFPM2', 'COPS2')
for (targ in KDtargs){
  
  up_res <- gofilter(up_GO_res[[targ]], level = lv)@result %>%
    mutate(target = targ) %>% filter(p.adjust < 0.05)
  rownames(up_res) <- NULL
  up_terms %<>% bind_rows(up_res %>% dplyr::select(ID, Description)) %>% unique() 
  
  if(length(up_res$Description) == 0){
    up_nores <- c(up_nores, targ)
  }
  
  dn_res <- gofilter(down_GO_res[[targ]], level = lv)@result %>%
    mutate(target = targ) #%>% filter(p.adjust < 0.05)
  rownames(dn_res) <- NULL
  dn_terms %<>% bind_rows(dn_res %>% dplyr::select(ID, Description)) %>% unique()
  
  if(length(dn_res$Description) == 0){
    dn_nores <- c(dn_nores, targ)
  }
  
  if (is.null(dim(all_up_res))){
    all_up_res <- as_tibble(up_res)
    all_dn_res <- as_tibble(dn_res)
  } else {
    all_up_res %<>% bind_rows(as_tibble(up_res))
    all_dn_res %<>% bind_rows(as_tibble(dn_res))
  }
  
}
all_up_res_plot <- all_up_res %>%
  dplyr::select(target, Description, p.adjust, qvalue) %>%
  bind_rows(tibble(target = KDtargs,
                   Description = 'GO-CC Term',
                   p.adjust = NA,
                   qvalue = NA)) %>%
  mutate(classOrder = case_when( #manual ordering
    grepl('sarc', Description) ~ 1,
    grepl('synap', Description) ~ 4,
    grepl('neur', Description) ~ 4,
    grepl('axon', Description) ~ 4,
    grepl('collag', Description) ~ 3,
    Description == 'GO-CC Term' ~ 5,
    TRUE ~ 2
  )) %>%
  arrange(-classOrder)

all_dn_res_plot <- all_dn_res %>%
  dplyr::select(target, Description, p.adjust, qvalue) %>%
  bind_rows(tibble(target = KDtargs,
                   Description = 'GO-CC Term',
                   p.adjust = NA,
                   qvalue = NA)) %>%
  mutate(classOrder = case_when( #manual ordering
    grepl('sarc', Description) ~ 1,
    grepl('contract', Description) ~ 1,
    grepl('fibril', Description) ~ 1,
    grepl('band', Description) ~ 1,
    grepl('intercalated disc', Description) ~ 1,
    grepl('synap', Description) ~ 4,
    grepl('neur', Description) ~ 4,
    grepl('axon', Description) ~ 4,
    grepl('collag', Description) ~ 3,
    grepl('matrix', Description) ~ 3,
    Description == 'GO-CC Term' ~ 5,
    TRUE ~ 2
  )) %>%
  arrange(-classOrder)

all_up_res_plot$target <- factor(all_up_res_plot$target)
all_up_res_plot$target <- factor(all_up_res_plot$target, levels = KDtargs)
all_up_res_plot$Description <- factor(all_up_res_plot$Description)
all_up_res_plot$Description <- factor(all_up_res_plot$Description, levels = rev(unique(all_up_res_plot$Description)))



all_dn_res_plot$target <- factor(all_dn_res_plot$target)
all_dn_res_plot$target <- factor(all_dn_res_plot$target, levels = KDtargs)
all_dn_res_plot$Description <- factor(all_dn_res_plot$Description)
all_dn_res_plot$Description <- factor(all_dn_res_plot$Description, levels = rev(unique(all_dn_res_plot$Description)))



all_up_res_plot_heatmap <- ggplot() +
  geom_tile(data = all_up_res_plot, aes(x = target, y = Description, fill = -log10(p.adjust))) +
  scale_fill_distiller(na.value = 'white', direction = 1) +
  scale_x_discrete(position = "top") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_blank(),
        legend.position = 'top')
ggsave(all_up_res_plot_heatmap, file = paste0(graphDir, '/iPS-CM-KDs/GOanalysis/heatmap_upRegulated_GO_pval2.pdf'), units = 'in', width = 4, height = 1+0.12*length(unique(all_up_res_plot$Description)))

# all_dn_res_plot_heatmap <- ggplot() +
#   geom_tile(data = all_dn_res_plot, aes(x = target, y = Description, fill = -log10(p.adjust))) +
#   scale_fill_distiller(na.value = 'white', direction = 1) +
#   scale_x_discrete(position = "top") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, size = 6),
#         axis.text.y = element_text(size = 6),
#         axis.title = element_blank(),
#         legend.position = 'top')
# ggsave(all_dn_res_plot_heatmap, file = paste0(graphDir, '/iPS-CM-KDs/GOanalysis/heatmap_downRegulated_GO_pval.pdf'), units = 'in', width = 3.5, height = 1+0.16*length(unique(all_up_res_plot$Description)))


