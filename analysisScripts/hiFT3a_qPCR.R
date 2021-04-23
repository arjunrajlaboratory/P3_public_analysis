#!/usr/bin/env Rscript
# don't run lines 3-13 if running from R workspace
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

library(tidyverse)
library(magrittr)

# projectDir = '~/Dropbox (RajLab)/Shared_IanM/cellid_201807_onward/'
# procDataSubdir = 'procDataScripted/'
# graphSubdir = 'graphs'

procDataDir = paste0(projectDir, 'extractedData/')
graphDir = paste0(projectDir, graphSubdir)

if (!dir.exists(paste0(graphDir, '/hiFT'))){
  dir.create(paste0(graphDir, '/hiFT'))
}
if (!dir.exists(paste0(graphDir, '/hiFT/IAMhiFT3a'))){
  dir.create(paste0(graphDir, '/hiFT/IAMhiFT3a'))
}

hits <- c('control', 'CEBPB', 'RUNX1', 'PRRX2')

# first, the colony counts
hiFT3a_colonies <- as_tibble(read.table(paste0(procDataDir, 'hiFT/IAMhiFT3a_ First pass set - candidate shRNAs with controls - 3 sh per target - Sheet1.tsv'), header = T, stringsAsFactors = F, skip = 110, sep = '\t')) %>%
  filter(!(shID %in% c('RUNX1sh1', 'SHC003'))) 

hiFT3a_colonies_sum <- hiFT3a_colonies %>%
  group_by(target, shIndex) %>%
  summarise(meanCount = mean(colonies),
            semCount = sd(colonies)/sqrt(length(colonies)),
            errbMax = meanCount + semCount, 
            errbMin = ifelse(meanCount - semCount > 0, meanCount - semCount, 0))

hiFT3a_colonies$target <- factor(hiFT3a_colonies$target)
hiFT3a_colonies$target <- factor(hiFT3a_colonies$target, levels = c('control', 'CEBPB', 'PRRX2', 'RUNX1', 'ID3', 'ALL'))

hiFT3a_colonies_sum$target <- factor(hiFT3a_colonies_sum$target)
hiFT3a_colonies_sum$target <- factor(hiFT3a_colonies_sum$target, levels = c('control', 'CEBPB', 'PRRX2', 'RUNX1', 'ID3', 'ALL'))


controlMeans <- hiFT3a_colonies %>%
  filter(target == 'control') %>%
  summarise(meanCtlCount = mean(colonies),
            semCtlCount = sd(colonies)/sqrt(length(colonies)))

colony_dotsNbars_forSlide <- ggplot() + 
  geom_point(data = hiFT3a_colonies %>% filter(target != 'PO'), aes(shIndex, colonies), position = position_jitter(w = 0.15, h = 0), alpha = 0.2) +
  geom_bar(data = hiFT3a_colonies_sum %>% filter(target != 'PO'), aes(shIndex, meanCount), alpha = 0.4, width = 0.7, stat = 'identity') +
  geom_errorbar(data = hiFT3a_colonies_sum %>% filter(target != 'PO'), aes(x = shIndex, ymin = errbMin, ymax = errbMax), width = 0.25, alpha = 0.5) +
  facet_grid(~target) +
  geom_hline(yintercept = controlMeans$meanCtlCount, color = 'blue', alpha = 0.8) +
  theme_classic() +
  ylab('Alkaline Phosphatase-positive colonies per well') +
  xlab('shRNA ID for target') +
  ggtitle('hiF-T iPSC colony counts after perturbable factor knockdown') +
  theme(axis.text.y = element_text(size = rel(2)),
        axis.title = element_blank(),
        plot.title = element_blank(),
        strip.text = element_text(size = rel(1.6)))
ggsave(colony_dotsNbars_forSlide, file = paste0(graphDir, '/hiFT/IAMhiFT3a/colony_counts_dotsNbars_forSlide.pdf'), width = 10, height = 5, useDingbats = F)  


# first, the colony counts
hiFT3a_colonies_hits <- hiFT3a_colonies %>%
  filter(target %in% hits,
         shID != 'SCH003') %>%
  mutate(shIndex = ifelse(shID == 'SHC007', 3, shIndex))

hiFT3a_colonies_hits_sum <- hiFT3a_colonies_hits %>%
  group_by(target, shIndex) %>%
  summarise(meanCount = mean(colonies),
            semCount = sd(colonies)/sqrt(length(colonies)),
            errbMax = meanCount + semCount, 
            errbMin = ifelse(meanCount - semCount > 0, meanCount - semCount, 0))

hiFT3a_colonies_hits$target <- factor(hiFT3a_colonies_hits$target)
hiFT3a_colonies_hits$target <- factor(hiFT3a_colonies_hits$target, levels = c('control', 'CEBPB', 'RUNX1', 'PRRX2'))

hiFT3a_colonies_hits_sum$target <- factor(hiFT3a_colonies_hits_sum$target)
hiFT3a_colonies_hits_sum$target <- factor(hiFT3a_colonies_hits_sum$target, levels = c('control', 'CEBPB', 'RUNX1', 'PRRX2'))


controlMeans <- hiFT3a_colonies_hits %>%
  filter(target == 'control') %>%
  summarise(meanCtlCount = mean(colonies),
            semCtlCount = sd(colonies)/sqrt(length(colonies)))

colony_dotsNbars_forSlide_hits <- ggplot() + 
  geom_point(data = hiFT3a_colonies_hits %>% filter(target != 'PO'), aes(as.character(shIndex), colonies), position = position_jitter(w = 0.15, h = 0), alpha = 0.2) +
  geom_bar(data = hiFT3a_colonies_hits_sum %>% filter(target != 'PO'), aes(as.character(shIndex), meanCount), alpha = 0.4, width = 0.7, stat = 'identity') +
  geom_errorbar(data = hiFT3a_colonies_hits_sum %>% filter(target != 'PO'), aes(x = as.character(shIndex), ymin = errbMin, ymax = errbMax), width = 0.25, alpha = 0.5) +
  facet_grid(~target) +
  geom_hline(yintercept = controlMeans$meanCtlCount, color = 'blue', alpha = 0.8) +
  theme_classic() +
  ylab('Alkaline Phosphatase-positive colonies per well') +
  xlab('shRNA ID for target') +
  ggtitle('hiF-T iPSC colony counts after perturbable factor knockdown') +
  theme(axis.text.y = element_text(size = rel(2)),
        axis.title = element_blank(),
        plot.title = element_blank(),
        strip.text = element_text(size = rel(1.6)),
        panel.spacing = unit(0.3, "lines"))
ggsave(colony_dotsNbars_forSlide_hits, file = paste0(graphDir, '/hiFT/IAMhiFT3a/colony_counts_hits_dotsNbars_forSlide.pdf'), width = 7, height = 5, useDingbats = F)  

#fisher's exact

targets = unique(as.character(hiFT3a_colonies$target))[!unique(as.character(hiFT3a_colonies$target)) %in% c('PO', 'control')]

ctlDat <- hiFT3a_colonies %>%
  filter(shID %in% c('SHC001', 'SHC002', 'SHC007'))

FET_pvals_3a <- list()
for (targ in targets){
  
  tempDat <- hiFT3a_colonies %>%
    filter(target == targ)
  
  tempFET <- fisher.test(matrix(c(sum(tempDat$colonies), length(tempDat$colonies)*1e4 - sum(tempDat$colonies),  
                                  sum(ctlDat$colonies), length(ctlDat$colonies)*1e4 - sum(ctlDat$colonies)), ncol = 2), alternative = 'g')
  
  if(is.null(dim(FET_pvals_3a))){
    FET_pvals_3a = tibble(target = targ,
                          FET_pval = tempFET$p.value)
  } else {
    FET_pvals_3a %<>% bind_rows(tibble(target = targ, FET_pval = tempFET$p.value))
  }
  
}

write.csv(FET_pvals_3a, file = paste0(graphDir, '/hiFT/IAMhiFT3a/FET_pvals_3a.csv'))

# 
# colonies <- colonies[,1:4]
# colonies$Target <- factor(colonies$Target)
# colonies$Target <- factor(colonies$Target, levels = c('SHC002', 'CEBPB', 'PRRX2', 'RUNX1', 'ID3', 'ALL'))
# colonyPlot <- ggplot(colonies, aes(as.character(shID), colonies)) +
#   facet_grid(.~Target) +
#   geom_jitter() +
#   theme_bw() +
#   theme(axis.text = element_text(size = rel(2)),
#     axis.title = element_blank(),
#     strip.text = element_text(size = rel(1.5)))
# plot(colonyPlot)
# ggsave(colonyPlot, file = '~/Dropbox (RajLab)/Projects/cellid/graphs/hiFT/IAMhiFT3a/colonyPlot.pdf', width = 9, height = 4)
# 
# coloniesSum <- colonies %>%
#   mutate(targSH = paste0(as.character(Target),'sh',as.character(shID))) %>%
#   filter(Target != "ALL") %>%
#   group_by(targSH) %>%
#   summarise(meanColCt = mean(colonies),
#             sdColCt = sd(colonies))
# coloniesSum$targSH <- factor(coloniesSum$targSH)
# coloniesSum$targSH <- factor(coloniesSum$targSH, levels = levels(coloniesSum$targSH)[c(10,1:9)])
# 
# gaps <- as_tibble(data.frame(
#   targSH = c('gap1', 'gap2', 'gap3', 'gap4'),
#   meanColCt = c(0,0,0,0),
#   sdColCt = c(0,0,0,0)
# ))
# 
# coloniesSum2 <- bind_rows(coloniesSum, gaps)
# 
# coloniesSum2$targSH <- factor(coloniesSum2$targSH)
# coloniesSum2$targSH <- factor(coloniesSum2$targSH, levels = levels(coloniesSum2$targSH)[c(14,3,1:2,4,7:9,5,10:11,6,12:13)])
# 
# 
# colonyPlot2 <- ggplot() +
#   geom_bar(data = coloniesSum2, stat = "identity", aes(targSH, meanColCt)) +
#   geom_errorbar(data = coloniesSum2, aes(targSH, ymin = ifelse(meanColCt - sdColCt > 0, meanColCt - sdColCt, 0), ymax = meanColCt + sdColCt), width = 0.25) + 
#   theme_classic() +
#   theme(axis.text = element_text(size = rel(2)),
#         axis.title = element_blank(),
#         strip.text = element_text(size = rel(1.5)))
# plot(colonyPlot2)
# ggsave(colonyPlot2, file = '~/Dropbox (RajLab)/Projects/cellid/graphs/hiFT/IAMhiFT3a/colonyPlot_bars.pdf', width = 9, height = 4)
# 
# 
# 
# hiFT3a_qPCRplate1_layout <- as_tibble(read.csv('~/Downloads/IAMhiFT3-qPCR1-layout - Sheet1 (1).csv', header = T, stringsAsFactors = F)) %>%
#   mutate(Well = paste0(Row, Col))
# hiFT3a_qPCRplate1_results <- as_tibble(read.csv('~/Downloads/IAMhiFT3a_qPCRplate1_results_someCols - Sheet1.csv', header = T, stringsAsFactors = F))
# 
# hiFT3a_qPCRplate2_layout <- as_tibble(read.csv('~/Downloads/IAMhiFT3a-qPCRplate2 20180830 - Sheet1.csv', header = T, stringsAsFactors = F, skip = 22)) %>%
#   mutate(Well = paste0(Row, Col))
# hiFT3a_qPCRplate2_results <- as_tibble(read.csv('~/Downloads/IAMhiFT3a_qPCRplate2_results_someCols - Sheet1.csv', header = T, stringsAsFactors = F))
# 
# 
# 
# gapdhCt <- inner_join(hiFT3a_qPCRplate1_layout, hiFT3a_qPCRplate1_results, by = 'Well') %>%
#   filter(!(Well %in% c('F10', 'G10', 'G11', 'G12', 'H9', 'H10', 'H11', 'H12')),
#          primer == 'GAPDH') %>%
#   group_by(sample, primer) %>%
#   summarise(meanGAPDHct = mean(as.numeric(Ct)))
# 
# CtPlot <- ggplot(inner_join(hiFT3a_qPCRplate1_layout, hiFT3a_qPCRplate1_results, by = 'Well') %>%
#                       filter(!(Well %in% c('F10', 'F12', 'G10', 'G11', 'G12', 'H9', 'H10', 'H11', 'H12')), primer != 'empty') %>%
#                       left_join(gapdhCt %>% dplyr::select(-primer), by = c('sample')), aes(sample, as.numeric(Ct))) +
#   facet_grid(primer~., scales = 'free') +
#   geom_point()
# plot(CtPlot)
# 
# gapdhCt2 <- inner_join(hiFT3a_qPCRplate2_layout, hiFT3a_qPCRplate2_results, by = 'Well') %>%
#   filter(primer == 'GAPDH') %>%
#   group_by(sample, primer) %>%
#   summarise(meanGAPDHct = mean(as.numeric(Ct)))
# 
# CtPlot2 <- ggplot(inner_join(hiFT3a_qPCRplate2_layout, hiFT3a_qPCRplate2_results, by = 'Well') %>%
#                    filter(primer != 'empty') %>%
#                    left_join(gapdhCt %>% dplyr::select(-primer), by = c('sample')), aes(sample, as.numeric(Ct))) +
#   facet_grid(primer~., scales = 'free') +
#   geom_point()
# plot(CtPlot2)
# 
# obs_dCt <- inner_join(hiFT3a_qPCRplate1_layout, hiFT3a_qPCRplate1_results, by = 'Well') %>%
#   filter(!(Well %in% c('F10', 'G10', 'G11', 'G12', 'H9', 'H10', 'H11', 'H12'))) %>%
#   left_join(gapdhCt %>% dplyr::select(-primer), by = c('sample')) %>%
#   mutate(dCt = as.numeric(Ct) - meanGAPDHct)  %>%
#   dplyr::select(sample, primer, repID, dCt) %>%
#   group_by(sample, repID) %>%
#   filter(!is.na(dCt)) %>%
#   spread(primer, dCt)
# 
# controls_CEBPB_dCt <- obs_dCt %>%
#   filter(sample %in% c('SHC002', 'SHC007')) %>%
#   dplyr::select(-c('GAPDH', 'ID3')) %>%
#   group_by(sample) %>%
#   summarise(meanCEBPBdCt = mean(CEBPB)) %>%
#   spread(sample, meanCEBPBdCt) %>%
#   dplyr::rename(SHC002_CEBPB = SHC002,
#                 SHC007_CEBPB = SHC007) %>%
#   mutate(plateID = 1)
# 
# controls_ID3_dCt <- obs_dCt %>%
#   filter(sample %in% c('SHC002', 'SHC007')) %>%
#   dplyr::select(-c('GAPDH', 'CEBPB')) %>%
#   group_by(sample) %>%
#   summarise(meanID3dCt = mean(ID3)) %>%
#   spread(sample, meanID3dCt) %>%
#   dplyr::rename(SHC002_ID3 = SHC002,
#                 SHC007_ID3 = SHC007) %>%
#   mutate(plateID = 1)
# 
# obs_ddCt <- obs_dCt %>%
#   mutate(plateID = 1) %>%
#   left_join(inner_join(controls_CEBPB_dCt, controls_ID3_dCt, by = 'plateID'), by = 'plateID') %>%
#   mutate(ddCt_CEBPB = CEBPB - SHC002_CEBPB,
#          ddCt_ID3 = ID3 - SHC002_ID3,
#          percent_reduction_CEBPB = 1 - 2^(-ddCt_CEBPB),
#          percent_reduction_ID3 = 1 - 2^(-ddCt_ID3),
#          inferred_percKD_infectedCells_CEBPB = percent_reduction_CEBPB/0.8,
#          inferred_percKD_infectedCells_ID3 = percent_reduction_ID3/0.8)
# 
# write.csv(obs_ddCt, file = 'graphs/hiFT/IAMhiFT3a/qPCRplate1_ddCt_perRep.csv', quote = F)
# 
# obs_mean_ddCt <- obs_dCt %>%
#   mutate(plateID = 1) %>%
#   left_join(inner_join(controls_CEBPB_dCt, controls_ID3_dCt, by = 'plateID'), by = 'plateID') %>%
#   mutate(ddCt_CEBPB = CEBPB - SHC002_CEBPB,
#          ddCt_ID3 = ID3 - SHC002_ID3) %>%
#   group_by(sample) %>%
#   summarise(mean_ddCt_CEBPB = mean(ddCt_CEBPB),
#             mean_ddCt_ID3 = mean(ddCt_ID3)) %>%
#   mutate(percent_reduction_CEBPB = 1 - 2^(-mean_ddCt_CEBPB),
#          percent_reduction_ID3 = 1 - 2^(-mean_ddCt_ID3),
#          inferred_percKD_infectedCells_CEBPB = percent_reduction_CEBPB/0.8,
#          inferred_percKD_infectedCells_ID3 = percent_reduction_ID3/0.8)
#          

            