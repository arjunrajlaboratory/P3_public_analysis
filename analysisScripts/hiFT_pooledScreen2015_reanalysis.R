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

procDataDir = paste0(projectDir, procDataSubdir)
graphDir = paste0(projectDir, graphSubdir)

if(!dir.exists(paste0(graphDir, '/hiFT/'))){
  dir.create(paste0(graphDir, '/hiFT/'))
}
if(!dir.exists(paste0(graphDir, '/hiFT/hiFT_pool_2015'))){
  dir.create(paste0(graphDir, '/hiFT/hiFT_pool_2015'))
}

library(tidyverse)
library(magrittr)
library(ggrepel)

setwd(projectDir)
tfGeneIDfile = "annotations/TF_gene_ids.csv"
tfTab = as_tibble(read.csv(tfGeneIDfile))
colnames(tfTab) = "gene_id"

geneIDfile = "annotations/hg19gene_idToGeneSymbol.tsv"
geneIDtoGeneSymbol <- as_tibble(read.csv(geneIDfile, sep='\t', stringsAsFactors = F))

## Reanalysis of Cacchiarelli, et al. pooled shRNA screen
# See processedData/hiFT/README.sh for data download information
# the dataset is a table with columns "ID"         "Hairpin"    "Gene"       "Sample_90"  "Sample_91"  "Sample_106" "Sample_104"

hits_in_paper <- c('EZH1', 'KDM1A', 'KTI12', 'LBR', 'NAP1L3', 'RSF1', 'SHPRH')
hits_in_pertub <- c('ZNF652', 'ZBTB38', 'ATOH8', 'CERS2', 'KLF13', 'CEBPB', 'RUNX1', 'PRRX2')

pooled_sh_hiFT <- as_tibble(read.table(paste0(projectDir, 'extractedData/hiFT/hiFT_pool_2015/hiFT_pooled_shRNA_counts.txt'), sep = '\t', header = T, stringsAsFactors = F))
indiv_perturbability = as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_variability_metrics.txt"), header = T, stringsAsFactors = F, sep = "\t"))
indiv_perturbability_RPMnorm20 = as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_tpmFilt_variability_RPMnormMetrics_window20.txt"), header = T, stringsAsFactors = F, sep = "\t"))
indiv_perturbability[is.na(indiv_perturbability)] <- 0
indiv_perturbability_RPMnorm20[is.na(indiv_perturbability_RPMnorm20)] <- 0

# based on replicate sample IDs in GSE62777 sample list description fields
pseudo = 0.0001
sum_pool_persh <- pooled_sh_hiFT %>% 
  group_by(Gene, ID, Hairpin) %>%
  summarise(mean_pre = (Sample_90 + Sample_91)/2,
            mean_TRA = (Sample_106 + Sample_104)/2) %>%
  filter(mean_pre > 0) %>%
  mutate(TRA_over_pre = (mean_TRA + pseudo) / (mean_pre + pseudo)) %>%
  arrange(-TRA_over_pre)

TRA95 <- quantile(sum_pool_persh$TRA_over_pre, probs = 0.95)
TRA90 <- quantile(sum_pool_persh$TRA_over_pre, probs = 0.90)

sum_pool_perGene <- sum_pool_persh %>%
  group_by(Gene) %>%
  summarise(mean_TRA_over_pre = mean(TRA_over_pre),
            frac_top05 = sum(TRA_over_pre > TRA95)/length(TRA_over_pre),
            frac_top10 = sum(TRA_over_pre > TRA90)/length(TRA_over_pre),
            num_non0_sh = length(TRA_over_pre)) %>%
  arrange(-mean_TRA_over_pre)

# based on replicate sample IDs in GSE62777 sample list description fields
sum_pool_persh_per_rep <- pooled_sh_hiFT %>% 
  group_by(Gene, ID, Hairpin) %>%
  summarise(TRA_over_pre_1 = (Sample_106 + pseudo)/(Sample_90 + pseudo),
            TRA_over_pre_2 = (Sample_104 + pseudo)/(Sample_91 + pseudo)) %>%
  mutate(mean_TRA_over_pre = (TRA_over_pre_1 + TRA_over_pre_2)/2) %>%
  arrange(-mean_TRA_over_pre)

TRA95 <- quantile(sum_pool_persh_per_rep$mean_TRA_over_pre, probs = 0.95)
TRA90 <- quantile(sum_pool_persh_per_rep$mean_TRA_over_pre, probs = 0.90)
TRA50 <- median(sum_pool_persh_per_rep$mean_TRA_over_pre)

sum_pool_perGene <- sum_pool_persh_per_rep %>%
  group_by(Gene) %>%
  summarise(mean_mean_TRA_over_pre = mean(mean_TRA_over_pre),
            median_mean_TRA_over_pre = median(mean_TRA_over_pre),
            frac_top05 = sum(mean_TRA_over_pre > TRA95)/length(mean_TRA_over_pre),
            frac_top10 = sum(mean_TRA_over_pre > TRA90)/length(mean_TRA_over_pre),
            frac_top50 = sum(mean_TRA_over_pre > TRA50)/length(mean_TRA_over_pre),
            num_non0_sh = length(mean_TRA_over_pre)) %>%
  arrange(-mean_mean_TRA_over_pre)

# filter 0-count pre-pool sh
sum_pool_persh_per_rep <- pooled_sh_hiFT %>% 
  group_by(Gene, ID, Hairpin) %>%
  summarise(TRA_over_pre_1 = ifelse(Sample_90 > 0, (Sample_106 + pseudo)/(Sample_90 + pseudo), NA),
            TRA_over_pre_2 = ifelse(Sample_91 > 0, (Sample_104 + pseudo)/(Sample_91 + pseudo), NA)) %>%
  mutate(mean_TRA_over_pre = mean(c(TRA_over_pre_1, TRA_over_pre_2), na.rm = T)) %>%
  arrange(-mean_TRA_over_pre)

TRA95 <- quantile(sum_pool_persh_per_rep$mean_TRA_over_pre, probs = 0.95, na.rm = T)
TRA90 <- quantile(sum_pool_persh_per_rep$mean_TRA_over_pre, probs = 0.90, na.rm = T)
TRA50 <- median(sum_pool_persh_per_rep$mean_TRA_over_pre, na.rm = T)
TRA10 <-quantile(sum_pool_persh_per_rep$mean_TRA_over_pre, probs = 0.10, na.rm = T)

sum_pool_perGene <- sum_pool_persh_per_rep %>%
  filter(!is.na(mean_TRA_over_pre)) %>%
  group_by(Gene) %>%
  summarise(mean_mean_TRA_over_pre = mean(mean_TRA_over_pre),
            median_mean_TRA_over_pre = median(mean_TRA_over_pre),
            max_mean_TRA_over_pre = max(mean_TRA_over_pre),
            frac_top05 = sum(mean_TRA_over_pre > TRA95)/length(mean_TRA_over_pre),
            frac_top10 = sum(mean_TRA_over_pre > TRA90)/length(mean_TRA_over_pre),
            frac_top50 = sum(mean_TRA_over_pre > TRA50)/length(mean_TRA_over_pre),
            frac_bot10 = sum(mean_TRA_over_pre < TRA10)/length(mean_TRA_over_pre),
            num_non0_sh = length(mean_TRA_over_pre)) %>%
  arrange(-mean_mean_TRA_over_pre)

# filter 0-count pre-pool sh and post-pool sh
sum_pool_persh_per_rep <- pooled_sh_hiFT %>% 
  group_by(Gene, ID, Hairpin) %>%
  summarise(TRA_over_pre_1 = ifelse(Sample_90 > 0 & Sample_106 > 0, (Sample_106 + pseudo)/(Sample_90 + pseudo), NA),
            TRA_over_pre_2 = ifelse(Sample_91 > 0 & Sample_104 > 0, (Sample_104 + pseudo)/(Sample_91 + pseudo), NA)) %>%
  mutate(mean_TRA_over_pre = mean(c(TRA_over_pre_1, TRA_over_pre_2), na.rm = T),
         max_TRA_over_pre = max(c(TRA_over_pre_1, TRA_over_pre_2), na.rm = T)) %>%
  arrange(-mean_TRA_over_pre)

TRA95 <- quantile(sum_pool_persh_per_rep$mean_TRA_over_pre, probs = 0.95, na.rm = T)
TRA90 <- quantile(sum_pool_persh_per_rep$mean_TRA_over_pre, probs = 0.90, na.rm = T)
TRA50 <- median(sum_pool_persh_per_rep$mean_TRA_over_pre, na.rm = T)
TRA10 <-quantile(sum_pool_persh_per_rep$mean_TRA_over_pre, probs = 0.10, na.rm = T)

TRA95x <- quantile(sum_pool_persh_per_rep$max_TRA_over_pre, probs = 0.95, na.rm = T)
TRA90x <- quantile(sum_pool_persh_per_rep$max_TRA_over_pre, probs = 0.90, na.rm = T)
TRA50x <- median(sum_pool_persh_per_rep$max_TRA_over_pre, na.rm = T)
TRA10x <-quantile(sum_pool_persh_per_rep$max_TRA_over_pre, probs = 0.10, na.rm = T)


sum_pool_perGene <- sum_pool_persh_per_rep %>%
  filter(!is.na(mean_TRA_over_pre)) %>%
  group_by(Gene) %>%
  summarise(mean_mean_TRA_over_pre = mean(mean_TRA_over_pre),
            median_mean_TRA_over_pre = median(mean_TRA_over_pre),
            mean_max_TRA_over_pre = mean(max_TRA_over_pre),
            max_mean_TRA_over_pre = max(mean_TRA_over_pre),
            median_max_TRA_over_pre = median(max_TRA_over_pre),
            frac_top05 = sum(mean_TRA_over_pre > TRA95)/length(mean_TRA_over_pre),
            frac_top10 = sum(mean_TRA_over_pre > TRA90)/length(mean_TRA_over_pre),
            frac_top50 = sum(mean_TRA_over_pre > TRA50)/length(mean_TRA_over_pre),
            frac_bot10 = sum(mean_TRA_over_pre < TRA10)/length(mean_TRA_over_pre),
            num_non0_sh = length(mean_TRA_over_pre)) %>%
  arrange(-mean_max_TRA_over_pre)

# without pseudo
sum_pool_persh_per_rep_nopseudo <- pooled_sh_hiFT %>% 
  group_by(Gene, ID, Hairpin) %>%
  summarise(TRA_over_pre_1 = ifelse(Sample_90 > 0, (Sample_106 )/(Sample_90 ), NA),
            TRA_over_pre_2 = ifelse(Sample_91 > 0, (Sample_104 )/(Sample_91 ), NA)) %>%
  mutate(mean_TRA_over_pre = mean(c(TRA_over_pre_1, TRA_over_pre_2), na.rm = T),
         max_TRA_over_pre = max(c(TRA_over_pre_1, TRA_over_pre_2), na.rm = T)) %>%
  arrange(-mean_TRA_over_pre)

# per Davide, just the ratio of normalized reads. 
# The results of the two replicates are simply averaged. 
# We did not aggregate the hairpins targeting the same gene, 
# but in figure 6a we simply ranked the individual shRNA result in the screen 
# and in the individual validation. 
# We ranked the averages and we considered for follow up in 6B only those genes 
# showing at least two effective shRNAs in the top 10%
pool_persh_per_rep_nopseudo <- pooled_sh_hiFT %>% 
  group_by(Gene, ID, Hairpin) %>%
  summarise(TRA_over_pre_1 = ifelse(Sample_90 > 0, Sample_106/Sample_90, 0),
            TRA_over_pre_2 = ifelse(Sample_91 > 0, Sample_104/Sample_91, 0)) %>%
  ungroup() %>%
  mutate(geom_mean_FC = sqrt(TRA_over_pre_1 * TRA_over_pre_2),
         rank_1 = rank(TRA_over_pre_1),
         rank_2 = rank(TRA_over_pre_2),
         avg_rank = (rank_1 + rank_2)/2,
         top10ile = quantile(avg_rank, probs = 0.9))

pool_persh_per_rep_maxPerGene <- pooled_sh_hiFT %>% 
  group_by(Gene, ID, Hairpin) %>%
  summarise(TRA_over_pre_1 = ifelse(Sample_90 > 0, Sample_106/Sample_90, 0),
            TRA_over_pre_2 = ifelse(Sample_91 > 0, Sample_104/Sample_91, 0)) %>%
  ungroup() %>%
  mutate(geom_mean_FC = sqrt(TRA_over_pre_1 * TRA_over_pre_2),
         rank_1 = rank(TRA_over_pre_1),
         rank_2 = rank(TRA_over_pre_2),
         avg_rank = (rank_1 + rank_2)/2,
         top10ile = quantile(avg_rank, probs = 0.9)) %>%
  group_by(Gene) %>%
  filter(avg_rank == max(avg_rank)) %>%
  dplyr::select(Gene, avg_rank, geom_mean_FC) %>%
  unique()


indiv_perturbability_pool <- indiv_perturbability %>%
  filter(cellType == 'GM00942') %>%
  inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
  inner_join(pool_persh_per_rep_maxPerGene %>% dplyr::rename(GeneSymbol = Gene), by = 'GeneSymbol')

indiv_perturbability_RPMnorm20_pool <- indiv_perturbability_RPMnorm20 %>%
  filter(cellType == 'GM00942') %>%
  inner_join(geneIDtoGeneSymbol, by = 'gene_id') %>%
  inner_join(pool_persh_per_rep_maxPerGene %>% dplyr::rename(GeneSymbol = Gene), by = 'GeneSymbol') 


fibro_nUP_avgPoolRank_withDCbarriers <- ggplot() + 
  geom_point(data = indiv_perturbability_pool, 
             aes(nDESeq2conditionsUp, avg_rank), alpha = 0.4) +
  geom_point(data = indiv_perturbability_pool %>% filter(GeneSymbol %in% hits_in_paper),
             aes(nDESeq2conditionsUp, avg_rank), color = 'darkolivegreen3') +
  geom_text_repel(data = indiv_perturbability_pool %>% filter(GeneSymbol %in% hits_in_paper),
                  aes(nDESeq2conditionsUp, avg_rank, label = GeneSymbol), color = 'darkolivegreen3') +
  xlab('# conditions in which UP\nin fibroblasts') +
  ylab('Average ranked enrichment\ninhiF-T pooled screen') +
  theme_bw() +
  xlim(c(0,10))

fibro_nUP_avgPoolRank_withP3barriers <- ggplot() + 
  geom_point(data = indiv_perturbability_pool, 
             aes(nDESeq2conditionsUp, avg_rank), alpha = 0.4) +
  geom_point(data = indiv_perturbability_pool %>% filter(GeneSymbol %in% hits_in_pertub),
             aes(nDESeq2conditionsUp, avg_rank), color = 'darkolivegreen3') +
  geom_text_repel(data = indiv_perturbability_pool %>% filter(GeneSymbol %in% hits_in_pertub),
                  aes(nDESeq2conditionsUp, avg_rank, label = GeneSymbol), color = 'darkolivegreen3') +
  xlab('# conditions in which UP\nin fibroblasts') +
  ylab('Average ranked enrichment\nin hiF-T pooled screen\nMaximum individual shRNA per gene') +
  theme_bw() +
  xlim(c(0,10))
  
set.seed(29347)
fibro_nUP_avgPoolRank_withP3barriers <- ggplot() + 
  geom_jitter(data = indiv_perturbability_pool, 
             aes(avg_rank, nDESeq2conditionsUp), alpha = 0.4, width = 20, height = 0) +
  ylab('# conditions in which UP\nin fibroblasts') +
  xlab('Average ranked enrichment\nin hiF-T pooled screen\nMaximum individual shRNA per gene') +
  theme_bw() 


fibro_nUPnorm_avgPoolRank_withDCbarriers <- ggplot() + 
  geom_point(data = indiv_perturbability_RPMnorm20_pool %>% filter(meanRPM > 20), 
             aes(norm_nDESeq2conditionsUp, avg_rank), alpha = 0.4) +
  geom_point(data = indiv_perturbability_RPMnorm20_pool %>% filter(GeneSymbol %in% hits_in_paper),
             aes(norm_nDESeq2conditionsUp, avg_rank), color = 'darkolivegreen3') +
  geom_text_repel(data = indiv_perturbability_RPMnorm20_pool %>% filter(GeneSymbol %in% hits_in_paper),
                  aes(norm_nDESeq2conditionsUp, avg_rank, label = GeneSymbol), color = 'darkolivegreen3') +
  xlab('# conditions in which UP\nin fibroblasts\nNormalized to meanRPM') +
  ylab('Average ranked enrichment\ninhiF-T pooled screen') +
  theme_bw()

pseudo = 1e-4
fibro_nUPnorm_avgPoollogFC_withDCbarriers <- ggplot() + 
  geom_point(data = indiv_perturbability_RPMnorm20_pool %>% filter(meanRPM > 20), 
             aes(norm_nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), alpha = 0.4) +
  geom_point(data = indiv_perturbability_RPMnorm20_pool %>% filter(GeneSymbol %in% hits_in_paper),
             aes(norm_nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), color = 'darkolivegreen3') +
  geom_text_repel(data = indiv_perturbability_RPMnorm20_pool %>% filter(GeneSymbol %in% hits_in_paper),
                  aes(norm_nDESeq2conditionsUp, log10(geom_mean_FC + pseudo), label = GeneSymbol), color = 'darkolivegreen3') +
  xlab('# conditions in which UP\nin fibroblasts\nNormalized to meanRPM') +
  ylab('Average enrichment (log10(sort/start))\ninhiF-T pooled screen') +
  theme_bw() +
  ylim(c(-4.2, 2.2))

fibro_nUPnorm_avgPoollogFC_tfOnly_withDCbarriers <- ggplot() + 
  geom_point(data = indiv_perturbability_RPMnorm20_pool %>% inner_join(tfTab) %>% filter(meanRPM > 20, !(GeneSymbol %in% hits_in_paper)), 
             aes(norm_nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), alpha = 0.4) +
  geom_point(data = indiv_perturbability_RPMnorm20_pool %>% inner_join(tfTab) %>% filter(GeneSymbol %in% hits_in_paper),
             aes(norm_nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), color = 'darkolivegreen3') +
  geom_text_repel(data = indiv_perturbability_RPMnorm20_pool %>% inner_join(tfTab) %>% filter(GeneSymbol %in% hits_in_paper),
                  aes(norm_nDESeq2conditionsUp, log10(geom_mean_FC + pseudo), label = GeneSymbol), color = 'darkolivegreen3') +
  xlab('# conditions in which UP\nin fibroblasts\nNormalized to meanRPM') +
  ylab('Average enrichment (log10(sort/start))\ninhiF-T pooled screen') +
  theme_bw() +
  ylim(c(-4.2, 2.2))

fibro_nUP_avgPoollogFC_withDCbarriers <- ggplot() + 
  geom_point(data = indiv_perturbability_pool %>% filter(meanRPM > 20, !(GeneSymbol %in% hits_in_paper)), 
             aes(nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), alpha = 0.4) +
  geom_point(data = indiv_perturbability_pool %>% filter(GeneSymbol %in% hits_in_paper),
             aes(nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), color = 'darkolivegreen3') +
  geom_text_repel(data = indiv_perturbability_pool %>% filter(GeneSymbol %in% hits_in_paper),
                  aes(nDESeq2conditionsUp, log10(geom_mean_FC + pseudo), label = GeneSymbol), color = 'darkolivegreen3') +
  xlab('# conditions in which UP\nin fibroblasts') +
  ylab('Average enrichment (log10(sort/start))\ninhiF-T pooled screen') +
  theme_bw() +
  ylim(c(-4.2, 2.2))

fibro_nUP_avgPoollogFC_withDCbarriers_boxplots <- ggplot() + 
  geom_violin(data = indiv_perturbability_pool %>% filter(meanRPM > 20), 
             aes(nDESeq2conditionsUp, log10(geom_mean_FC + pseudo), group = nDESeq2conditionsUp), alpha = 0.4) +
  # geom_point(data = indiv_perturbability_pool %>% filter(GeneSymbol %in% hits_in_paper),
  #            aes(nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), color = 'darkolivegreen3') +
  # geom_text_repel(data = indiv_perturbability_pool %>% filter(GeneSymbol %in% hits_in_paper),
  #                 aes(nDESeq2conditionsUp, log10(geom_mean_FC + pseudo), label = GeneSymbol), color = 'darkolivegreen3') +
  xlab('# conditions in which UP\nin fibroblasts') +
  ylab('Average enrichment (log10(sort/start))\ninhiF-T pooled screen') +
  theme_bw() +
  ylim(c(-4.2, 2.2))

set.seed(349978)
fibro_nUPnorm_avgPoollogFC_withDCbarriers_jitter <- ggplot() + 
  geom_jitter(data = indiv_perturbability_RPMnorm20_pool %>% filter(meanRPM > 20, !(GeneSymbol %in% hits_in_paper)), 
             aes(norm_nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), alpha = 0.4, width = 0.005, height = 0.05) +
  geom_point(data = indiv_perturbability_RPMnorm20_pool %>% filter(GeneSymbol %in% hits_in_paper),
             aes(norm_nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), color = 'darkolivegreen3') +
  geom_text_repel(data = indiv_perturbability_RPMnorm20_pool %>% filter(GeneSymbol %in% hits_in_paper),
                  aes(norm_nDESeq2conditionsUp, log10(geom_mean_FC + pseudo), label = GeneSymbol), color = 'darkolivegreen3') +
  xlab('# conditions in which UP\nin fibroblasts\nNormalized to meanRPM') +
  ylab('Average enrichment (log10(sort/start))\ninhiF-T pooled screen') +
  theme_bw() +
  ylim(c(-4.2, 2.2))
ggsave(fibro_nUPnorm_avgPoollogFC_withDCbarriers_jitter, file = paste0(graphDir, '/hiFT/hiFT_pool_2015/fibro_nUPnorm_avgPoollogFC_withDCbarriers_jitter.pdf'), width = 7, height = 7, useDingbats = F)  

set.seed(349978)
fibro_nUPnorm_avgPoollogFC_tfOnly_withDCbarriers_jitter <- ggplot() + 
  geom_jitter(data = indiv_perturbability_RPMnorm20_pool %>% inner_join(tfTab) %>% filter(meanRPM > 20, !(GeneSymbol %in% hits_in_paper)), 
              aes(norm_nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), alpha = 0.4, width = 0.005, height = 0.05) +
  geom_point(data = indiv_perturbability_RPMnorm20_pool %>% inner_join(tfTab) %>% filter(GeneSymbol %in% hits_in_paper),
             aes(norm_nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), color = 'darkolivegreen3') +
  geom_text_repel(data = indiv_perturbability_RPMnorm20_pool %>% inner_join(tfTab) %>% filter(GeneSymbol %in% hits_in_paper),
                  aes(norm_nDESeq2conditionsUp, log10(geom_mean_FC + pseudo), label = GeneSymbol), color = 'darkolivegreen3') +
  xlab('# conditions in which UP\nin fibroblasts\nNormalized to meanRPM') +
  ylab('Average enrichment (log10(sort/start))\ninhiF-T pooled screen') +
  theme_bw() +
  ylim(c(-4.2, 2.2))
ggsave(fibro_nUPnorm_avgPoollogFC_tfOnly_withDCbarriers_jitter, file = paste0(graphDir, '/hiFT/hiFT_pool_2015/fibro_nUPnorm_avgPoollogFC_tfOnly_withDCbarriers_jitter.pdf'), width = 7, height = 7, useDingbats = F)  



set.seed(349978)
fibro_nUPnorm_avgPoollogFC_onlyDCbarriers <- ggplot() + 
  # geom_jitter(data = indiv_perturbability_RPMnorm20_pool %>% filter(meanRPM > 20, !(GeneSymbol %in% hits_in_paper)), 
  #             aes(norm_nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), alpha = 0.4, width = 0.005, height = 0.05) +
  geom_point(data = indiv_perturbability_RPMnorm20_pool %>% filter(GeneSymbol %in% hits_in_paper),
             aes(norm_nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), color = 'darkolivegreen3') +
  geom_text_repel(data = indiv_perturbability_RPMnorm20_pool %>% filter(GeneSymbol %in% hits_in_paper),
                  aes(norm_nDESeq2conditionsUp, log10(geom_mean_FC + pseudo), label = GeneSymbol), color = 'darkolivegreen3') +
  xlab('# conditions in which UP\nin fibroblasts\nNormalized to meanRPM') +
  ylab('Average enrichment (log10(sort/start))\ninhiF-T pooled screen') +
  theme_bw() 
ggsave(fibro_nUPnorm_avgPoollogFC_onlyDCbarriers, file = paste0(graphDir, '/hiFT/hiFT_pool_2015/fibro_nUPnorm_avgPoollogFC_onlyDCbarriers.pdf'), width = 3, height = 3, useDingbats = F)  

# set.seed(349978)
fibro_nUPnorm_avgPoollogFC_withDCbarriers_contour <- ggplot() + 
  # geom_jitter(data = indiv_perturbability_RPMnorm20_pool %>% filter(meanRPM > 20), 
  #             aes(norm_nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), alpha = 0.4, height = 0.1, width = 0.05) +
  geom_density2d(data = indiv_perturbability_RPMnorm20_pool %>% filter(meanRPM > 20, !(GeneSymbol %in% hits_in_paper)), 
                 aes(norm_nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), alpha = 0.4, color = 'black') +
  geom_point(data = indiv_perturbability_RPMnorm20_pool %>% filter(GeneSymbol %in% hits_in_paper),
             aes(norm_nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), color = 'darkolivegreen3') +
  geom_text_repel(data = indiv_perturbability_RPMnorm20_pool %>% filter(GeneSymbol %in% hits_in_paper),
                  aes(norm_nDESeq2conditionsUp, log10(geom_mean_FC + pseudo), label = GeneSymbol), color = 'darkolivegreen3') +
  xlab('# conditions in which UP\nin fibroblasts\nNormalized to meanRPM') +
  ylab('Average enrichment (log10(sort/start))\ninhiF-T pooled screen') +
  theme_bw() +
  ylim(c(-4.2, 2.2))
ggsave(fibro_nUPnorm_avgPoollogFC_withDCbarriers_contour, file = paste0(graphDir, '/hiFT/hiFT_pool_2015/fibro_nUPnorm_avgPoollogFC_withDCbarriers_contour.pdf'), width = 7, height = 7, useDingbats = F)  


set.seed(349978)
fibro_nUP_avgPoollogFC_withDCbarriers_jitter <- ggplot() + 
  geom_jitter(data = indiv_perturbability_pool %>% filter(meanRPM > 20, !(GeneSymbol %in% hits_in_paper)), 
             aes(nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), alpha = 0.4, height = 0.1, width = 0.3) +
  geom_point(data = indiv_perturbability_pool %>% filter(GeneSymbol %in% hits_in_paper),
             aes(nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), color = 'darkolivegreen3') +
  geom_text_repel(data = indiv_perturbability_pool %>% filter(GeneSymbol %in% hits_in_paper),
                  aes(nDESeq2conditionsUp, log10(geom_mean_FC + pseudo), label = GeneSymbol), color = 'darkolivegreen3') +
  xlab('# conditions in which UP\nin fibroblasts') +
  ylab('Average enrichment (log10(sort/start))\ninhiF-T pooled screen') +
  theme_bw() +
  ylim(c(-4.2, 2.2))
ggsave(fibro_nUP_avgPoollogFC_withDCbarriers_jitter, file = paste0(graphDir, '/hiFT/hiFT_pool_2015/fibro_nUP_avgPoollogFC_withDCbarriers_jitter.pdf'), width = 10, height = 10, useDingbats = F)  


set.seed(349978)
fibro_nUP_avgPoollogFC_tfOnly_withDCbarriers_jitter <- ggplot() + 
  geom_jitter(data = indiv_perturbability_pool %>% inner_join(tfTab) %>% filter(meanRPM > 20, !(GeneSymbol %in% hits_in_paper)), 
              aes(nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), alpha = 0.4, height = 0.1, width = 0.3) +
  geom_point(data = indiv_perturbability_pool %>% inner_join(tfTab) %>% filter(GeneSymbol %in% hits_in_paper),
             aes(nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), color = 'darkolivegreen3') +
  geom_text_repel(data = indiv_perturbability_pool %>% inner_join(tfTab) %>% filter(GeneSymbol %in% hits_in_paper),
                  aes(nDESeq2conditionsUp, log10(geom_mean_FC + pseudo), label = GeneSymbol), color = 'darkolivegreen3') +
  xlab('# conditions in which UP\nin fibroblasts') +
  ylab('Average enrichment (log10(sort/start))\ninhiF-T pooled screen') +
  theme_bw() +
  ylim(c(-4.2, 2.2))
ggsave(fibro_nUP_avgPoollogFC_tfOnly_withDCbarriers_jitter, file = paste0(graphDir, '/hiFT/hiFT_pool_2015/fibro_nUP_avgPoollogFC_tfOnly_withDCbarriers_jitter.pdf'), width = 10, height = 10, useDingbats = F)  


fibro_nUP_avgPoollogFC_onlyDCbarriers <- ggplot() + 
  # geom_jitter(data = indiv_perturbability_pool %>% filter(meanRPM > 20, !(GeneSymbol %in% hits_in_paper)), 
  #             aes(nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), alpha = 0.4, height = 0.1, width = 0.3) +
  geom_point(data = indiv_perturbability_pool %>% filter(GeneSymbol %in% hits_in_paper),
             aes(nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), color = 'darkolivegreen3') +
  geom_text_repel(data = indiv_perturbability_pool %>% filter(GeneSymbol %in% hits_in_paper),
                  aes(nDESeq2conditionsUp, log10(geom_mean_FC + pseudo), label = GeneSymbol), color = 'darkolivegreen3') +
  xlab('# conditions in which UP\nin fibroblasts') +
  ylab('Average enrichment (log10(sort/start))\ninhiF-T pooled screen') +
  theme_bw() 
ggsave(fibro_nUP_avgPoollogFC_onlyDCbarriers, file = paste0(graphDir, '/hiFT/hiFT_pool_2015/fibro_nUP_avgPoollogFC_onlyDCbarriers.pdf'), width = 3, height = 3, useDingbats = F)  

fibro_nUP_avgPoollogFC_onlyDCbarriers <- ggplot() + 
  # geom_jitter(data = indiv_perturbability_pool %>% filter(meanRPM > 20, !(GeneSymbol %in% hits_in_paper)), 
  #             aes(nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), alpha = 0.4, height = 0.1, width = 0.3) +
  geom_point(data = indiv_perturbability_pool %>% filter(GeneSymbol %in% hits_in_paper),
             aes(nDESeq2conditionsUp, log10(geom_mean_FC + pseudo)), color = 'darkolivegreen3') +
  geom_text_repel(data = indiv_perturbability_pool %>% filter(GeneSymbol %in% hits_in_paper),
                  aes(nDESeq2conditionsUp, log10(geom_mean_FC + pseudo), label = GeneSymbol), color = 'darkolivegreen3') +
  xlab('# conditions in which UP\nin fibroblasts') +
  ylab('Average enrichment (log10(sort/start))\ninhiF-T pooled screen') +
  theme_bw() 
ggsave(fibro_nUP_avgPoollogFC_onlyDCbarriers, file = paste0(graphDir, '/hiFT/hiFT_pool_2015/fibro_nUP_avgPoollogFC_onlyDCbarriers.pdf'), width = 3, height = 3, useDingbats = F)  



ct <- cor.test(log10(indiv_perturbability_pool$geom_mean_FC+1e-4), indiv_perturbability_pool$nDESeq2conditionsUp, method = 's')

indiv_perturbability_pool_noDO <- indiv_perturbability_pool %>%
  filter(geom_mean_FC > 0)
indiv_perturbability_pool_noDO_tfs <- indiv_perturbability_pool %>%
  filter(geom_mean_FC > 0) %>%
  inner_join(tfTab)
indiv_perturbability_RPMnorm20_pool_noDO <- indiv_perturbability_RPMnorm20_pool %>%
  filter(geom_mean_FC > 0)
hits_norm <- indiv_perturbability_RPMnorm20_pool %>% filter(GeneSymbol %in% hits_in_paper)
hits_res <- indiv_perturbability_pool %>% filter(GeneSymbol %in% hits_in_paper)

ct_noDO <- cor.test(log10(indiv_perturbability_pool_noDO$geom_mean_FC), indiv_perturbability_pool_noDO$nDESeq2conditionsUp, method = 's')
ct_noDO_norm <- cor.test(log10(indiv_perturbability_RPMnorm20_pool_noDO$geom_mean_FC), indiv_perturbability_RPMnorm20_pool_noDO$norm_nDESeq2conditionsUp, method = 's')
ct_hits_norm <- cor.test(log10(hits_norm$geom_mean_FC), hits_norm$norm_nDESeq2conditionsUp, method = 's')
ct_hits <- cor.test(log10(hits_res$geom_mean_FC), hits_res$nDESeq2conditionsUp, method = 's')

# write.table(ct, file = paste0(graphDir, '/hiFT/hiFT_pool_2015/pool_FCvsnUP_cortest.txt'), quote = F)

# set.seed(26053)
# fibro_markers_updown_jitter_min20rpm <- indiv_perturbability %>% inner_join(geneIDtoGeneSymbol, by = "gene_id") %>%
#   filter(cellType == "GM00942", meanRPM > 20, GeneSymbol %in% targs16) %>% unique() %>% mutate(isBarrier = GeneSymbol %in% targs8) %>%
#   dplyr::select(gene_id, meanRPM, meanTPM, GeneSymbol, isBarrier, nDESeq2conditionsAll, nDESeq2conditionsUp, nDESeq2conditionsDown)
# fibro_markers_updown_jitter_min20rpm$nDESeq2conditionsAll_jit = fibro_markers_updown_jitter_min20rpm$nDESeq2conditionsAll + runif(nrow(fibro_markers_updown_jitter_min20rpm), min = -0.2, max = 0.2)
# fibro_markers_updown_jitter_min20rpm$nDESeq2conditionsUp_jit = fibro_markers_updown_jitter_min20rpm$nDESeq2conditionsUp + runif(nrow(fibro_markers_updown_jitter_min20rpm), min = -0.2, max = 0.2)
# fibro_markers_updown_jitter_min20rpm$nDESeq2conditionsDown_jit = fibro_markers_updown_jitter_min20rpm$nDESeq2conditionsDown + runif(nrow(fibro_markers_updown_jitter_min20rpm), min = -0.2, max = 0.2)

fibro_pert_for_hists <- indiv_perturbability_pool %>% 
  # inner_join(tfTab, by = "gene_id") %>% 
  filter(cellType == "GM00942", meanRPM > 20) %>% 
  dplyr::select(gene_id, nDESeq2conditionsAll, nDESeq2conditionsUp) %>%
  gather('Measure', 'nDrugs', 2:3)

# fibro_markers_for_hists <- fibro_markers_updown_jitter_min20rpm %>%
#   dplyr::select(GeneSymbol, isBarrier, nDESeq2conditionsAll_jit, nDESeq2conditionsUp_jit) %>%
#   dplyr::rename(nDESeq2conditionsAll = nDESeq2conditionsAll_jit,
#                 nDESeq2conditionsUp = nDESeq2conditionsUp_jit) %>%
#   gather('Measure', 'nDrugs', 3:4)
# set.seed(26054)
fibro_epiRegs_nAll_nUp_hists_noLabs_thin <- ggplot() +
  geom_histogram(data = fibro_pert_for_hists, aes(nDrugs), fill = 'grey', alpha = 0.9,
                 binwidth = 1) +
  # geom_linerange(data = fibro_markers_for_hists %>% filter(isBarrier == T), aes(nDrugs, ymin = 0, ymax = 40), size = 0.2, color = "green") +
  # geom_linerange(data = fibro_markers_for_hists %>% filter(isBarrier == F), aes(nDrugs, ymin = 0, ymax = 40), size = 0.2, color = "darkolivegreen2") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 0, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  # ylim(c(-100, 200)) +
  facet_grid(. ~ Measure, scales = 'free_x') +
  theme(axis.text = element_blank(),
        axis.title = element_blank())
ggsave(fibro_epiRegs_nAll_nUp_hists_noLabs_thin, file =  paste0(graphDir,'/differentialExpression/allSamples/fibroblasts/fibro_epiRegs_nAll_nUp_hists_noLabs_thin.pdf'), width = 4, height = 1, useDingbats = F)

fibro_epiRegs_nAll_nUp_hists_wLabs_thin <- ggplot() +
  geom_histogram(data = fibro_pert_for_hists, aes(nDrugs), fill = 'grey', alpha = 0.9,
                 binwidth = 1) +
  # geom_linerange(data = fibro_markers_for_hists %>% filter(isBarrier == T), aes(nDrugs, ymin = 0, ymax = 40), size = 0.2, color = "green") +
  # geom_linerange(data = fibro_markers_for_hists %>% filter(isBarrier == F), aes(nDrugs, ymin = 0, ymax = 40), size = 0.2, color = "darkolivegreen2") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 0, label = GeneSymbol), color = "red") +
  # geom_vline(data = markers_for_hists, aes(xintercept = nDrugs), color = "red") +
  # geom_text_repel(data = markers_for_hists, aes(x = nDrugs, y = 200, label = GeneSymbol), color = "red") +
  theme_classic() + 
  # ylim(c(-100, 200)) +
  facet_grid(. ~ Measure, scales = 'free_x') #+
  # theme(axis.text = element_blank(),
  #       axis.title = element_blank())
ggsave(fibro_epiRegs_nAll_nUp_hists_wLabs_thin, file =  paste0(graphDir,'/differentialExpression/allSamples/fibroblasts/fibro_epiRegs_nAll_nUp_hists_wLabs_thin.pdf'), width = 4, height = 1, useDingbats = F)



