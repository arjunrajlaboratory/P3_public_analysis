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
library(magrittr)

# projectDir = '~/Dropbox (RajLab)/Shared_IanM/cellid_201807_onward/'
# procDataSubdir = 'procDataScripted/'
# graphSubdir = 'graphs'

procDataDir = paste0(projectDir, 'extractedData/')
graphDir = paste0(projectDir, graphSubdir)

if (!dir.exists(paste0(graphDir, '/hiFT'))){
  dir.create(paste0(graphDir, '/hiFT'))
}
if (!dir.exists(paste0(graphDir, '/hiFT/IAMhiFT4a'))){
  dir.create(paste0(graphDir, '/hiFT/IAMhiFT4a'))
}

hits <- c('control', 'ZNF652', 'ZBTB38', 'ATOH8', 'CERS2', 'KLF13')

hiFT4a_colonies <- as_tibble(read.table(paste0(procDataDir, 'hiFT/IAMhiFT4a_ second set of shRNAs - Sheet1.tsv'), header = T, stringsAsFactors = F, skip = 69, sep = '\t'))
hiFT4a_colonies <- hiFT4a_colonies[,1:5]

hiFT4a_colonies_hits <- hiFT4a_colonies %>% filter(target %in% hits, shID != 'SHC003') %>% 
  mutate(shIndex = ifelse(shID == 'SHC007', 3, shIndex))

hiFT4a_colonies_hits_sum <- hiFT4a_colonies_hits %>%
  group_by(target, shIndex) %>%
  summarise(meanCount = mean(colonies),
            semCount = sd(colonies)/sqrt(length(colonies)),
            errbMax = meanCount + semCount, 
            errbMin = ifelse(meanCount - semCount > 0, meanCount - semCount, 0))

hiFT4a_colonies_hits$target <- factor(hiFT4a_colonies_hits$target)
hiFT4a_colonies_hits$target <- factor(hiFT4a_colonies_hits$target, levels = c('control', 'ZNF652', 'ZBTB38', 'ATOH8', 'CERS2', 'KLF13', 'CEBPB', 'RUNX1', 'PRRX2'))

hiFT4a_colonies_hits_sum$target <- factor(hiFT4a_colonies_hits_sum$target)
hiFT4a_colonies_hits_sum$target <- factor(hiFT4a_colonies_hits_sum$target, levels = c('control', 'ZNF652', 'ZBTB38', 'ATOH8', 'CERS2', 'KLF13', 'CEBPB', 'RUNX1', 'PRRX2'))

controlMeans <- hiFT4a_colonies %>%
  filter(target == 'control') %>%
  summarise(meanCtlCount = mean(colonies),
            semCtlCount = sd(colonies)/sqrt(length(colonies)))

colony_dotsNbars <- ggplot() + 
  geom_point(data = hiFT4a_colonies_hits %>% filter(target != 'PO'), aes(shIndex, colonies), position = position_jitter(w = 0.15, h = 0), alpha = 0.2) +
  geom_bar(data = hiFT4a_colonies_hits_sum %>% filter(target != 'PO'), aes(shIndex, meanCount), alpha = 0.4, width = 0.7, stat = 'identity') +
  geom_errorbar(data = hiFT4a_colonies_hits_sum %>% filter(target != 'PO'), aes(x = shIndex, ymin = errbMin, ymax = errbMax), width = 0.25, alpha = 0.5) +
  facet_grid(~target) +
  geom_hline(yintercept = controlMeans$meanCtlCount, color = 'blue', alpha = 0.4) +
  theme_classic() +
  ylab('Alkaline Phosphatase-positive colonies per well') +
  xlab('shRNA ID for target') +
  ggtitle('hiF-T iPSC colony counts after perturbable factor knockdown') +
  theme(axis.text.y = element_text(size = rel(2)),
        axis.title = element_text(size = rel(1.5)),
        panel.spacing = unit(0.3, "lines"))
ggsave(colony_dotsNbars, file = paste0(graphDir, '/hiFT/IAMhiFT4a/colony_counts_dotsNbars_forSlide1.pdf'), width = 10, height = 7, useDingbats = F)  

colony_dotsNbars_forSlide <- ggplot() + 
  geom_point(data = hiFT4a_colonies_hits %>% filter(target != 'PO'), aes(as.character(shIndex), colonies), position = position_jitter(w = 0.15, h = 0), alpha = 0.2) +
  geom_bar(data = hiFT4a_colonies_hits_sum %>% filter(target != 'PO'), aes(as.character(shIndex), meanCount), alpha = 0.4, width = 0.7, stat = 'identity') +
  geom_errorbar(data = hiFT4a_colonies_hits_sum %>% filter(target != 'PO'), aes(x = as.character(shIndex), ymin = errbMin, ymax = errbMax), width = 0.25, alpha = 0.5) +
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
ggsave(colony_dotsNbars_forSlide, file = paste0(graphDir, '/hiFT/IAMhiFT4a/colony_counts_dotsNbars_forSlide.pdf'), width = 10, height = 5, useDingbats = F)  

#fisher's exact

targets = unique(as.character(hiFT4a_colonies$target))[!unique(as.character(hiFT4a_colonies$target)) %in% c('PO', 'control')]

ctlDat <- hiFT4a_colonies %>%
  filter(shID %in% c('SHC001', 'SHC002', 'SHC007'))

FET_pvals_4a <- list()
for (targ in targets){
  
  tempDat <- hiFT4a_colonies %>%
    filter(target == targ)
  
  tempFET <- fisher.test(matrix(c(sum(tempDat$colonies), length(tempDat$colonies)*1e4 - sum(tempDat$colonies),  
                                  sum(ctlDat$colonies), length(ctlDat$colonies)*1e4 - sum(ctlDat$colonies)), ncol = 2), alternative = 'g')
  
  if(is.null(dim(FET_pvals_4a))){
    FET_pvals_4a = tibble(target = targ,
                          FET_pval = tempFET$p.value)
  } else {
    FET_pvals_4a %<>% bind_rows(tibble(target = targ, FET_pval = tempFET$p.value))
  }
  
}

write.csv(FET_pvals_4a, file = paste0(graphDir, '/hiFT/IAMhiFT4a/FET_pvals_4a.csv'))

# 
# # qPCR ####
# 
# genes <- c('SKIL', 'YBX3', 'TSC22D1', 'ID1', 'KLF13', 'CERS2', 'LARP1', 'TBX3', 'ZBTB38', 'ATOH8', 'ZNF652', 'NFATC4')
# plates <- 1:6
# exp <- 'IAMhiFT4a'
# 
# all_results<-list()
# for (plate in plates) {
#   
#   genesInPlate = genes[c(2*(plate-1) + 1, 2*(plate-1) + 2)]
#   
#   results <- as_tibble(read.csv(paste0('~/Google Drive/Cell ID/Reprogramming designs/hiF-T/qPCR results/', exp, '-qPCRplate', as.character(plate), '_somecols - Sheet1.csv'), header = T, stringsAsFactors = F))
#   
#   if (plate != 4) {
#   samples <- as_tibble(read.csv('~/Downloads/IAMhiFT4a-qPCRplate1 20181018 - Sheet1.csv', skip = 1, stringsAsFactors = F, header = T, nrows = 8)) %>%
#     dplyr::rename(row = cDNA) %>%
#     gather(col, sample, 2:13) %>%
#     mutate(col = substring(col, 2)) %>%
#     mutate(target = ifelse(sample %in% c('SKILsh1', 'SKILsh2', 'SKILsh3'), 'SKIL', 
#                            ifelse(sample %in% c('YBX3sh2', 'YBX3sh3', 'YBX3sh5'), 'YBX3',
#                                   ifelse(sample %in% c('SHC001', 'SHC002', 'SHC003', 'SHC007'), 'control', 'empty'))),
#            shID = ifelse(sample %in% c('SKILsh1', 'YBX3sh2', 'SHC001'), 1,
#                          ifelse(sample %in% c('SKILsh2', 'YBX3sh3', 'SHC002'), 2,
#                                 ifelse(sample %in% c('SKILsh3', 'YBX3sh5', 'SHC003'), 3,
#                                        ifelse(sample == 'SHC007', 4, NA)))))
#   
#   primers <-  as_tibble(read.csv('~/Downloads/IAMhiFT4a-qPCRplate1 20181018 - Sheet1.csv', skip = 11, stringsAsFactors = F, header = T, nrows = 8)) %>%
#     dplyr::rename(row = Primers) %>%
#     gather(col, primer, 2:13) %>%
#     mutate(col = substring(col, 2))
#   
#   } else {
#     samples <- as_tibble(read.csv('~/Downloads/IAMhiFT4a-qPCRplate1 20181018 ALT - Sheet1.csv', skip = 1, stringsAsFactors = F, header = T, nrows = 8)) %>%
#       dplyr::rename(row = cDNA) %>%
#       gather(col, sample, 2:13) %>%
#       mutate(col = substring(col, 2)) %>%
#       mutate(target = ifelse(sample %in% c('SKILsh1', 'SKILsh2', 'SKILsh3'), 'SKIL', 
#                              ifelse(sample %in% c('YBX3sh2', 'YBX3sh3', 'YBX3sh5'), 'YBX3',
#                                     ifelse(sample %in% c('SHC001', 'SHC002', 'SHC003', 'SHC007'), 'control', 'empty'))),
#              shID = ifelse(sample %in% c('SKILsh1', 'YBX3sh2', 'SHC001'), 1,
#                            ifelse(sample %in% c('SKILsh2', 'YBX3sh3', 'SHC002'), 2,
#                                   ifelse(sample %in% c('SKILsh3', 'YBX3sh5', 'SHC003'), 3,
#                                          ifelse(sample == 'SHC007', 4, NA)))))
#     
#     primers <-  as_tibble(read.csv('~/Downloads/IAMhiFT4a-qPCRplate1 20181018 ALT - Sheet1.csv', skip = 11, stringsAsFactors = F, header = T, nrows = 8)) %>%
#       dplyr::rename(row = Primers) %>%
#       gather(col, primer, 2:13) %>%
#       mutate(col = substring(col, 2))
#     
#   }
#   
#   
#   
#   layout_results <- inner_join(samples, primers, by = c('row', 'col')) %>%
#     mutate(target = ifelse(target == 'SKIL', genesInPlate[1],
#                            ifelse(target == 'YBX3', genesInPlate[2], 'control')),
#            primers = ifelse(primer == 'SKIL', genesInPlate[1],
#                             ifelse(primer == 'YBX3', genesInPlate[2], primer))) %>%
#     mutate(Well = paste0(row, col)) %>%
#     inner_join(results, by = 'Well') %>%
#     dplyr::select(-c(sample, Sample.Name, Task)) %>%
#     mutate(plateID = plate)
#   
#   
#   if (is.null(dim(all_results))) {
#     all_results <- layout_results
#   } else {
#     all_results %<>% bind_rows(layout_results)
#   }
#   
# }
# 
# CTC <- all_results %>%
#   filter(target == 'control', !is.na(shID), shID != 3) %>%
#   group_by(plateID, primers, shID) %>%
#   summarise(ControlMeanCt = mean(as.numeric(Ct), na.rm = T))
# 
# CTE <- all_results %>%
#   filter(target != 'control', !is.na(shID)) %>%
#   group_by(plateID, target, primers, shID) %>%
#   summarise(ControlMeanCt = ifelse(sum(is.na(as.numeric(Ct))) > 1, 40, mean(as.numeric(Ct), na.rm = T)))
# 
# HE <- CTE %>%
#   filter(primers == 'GAPDH') %>%
#   dplyr::rename(HE = ControlMeanCt)
# TE <- CTE %>%
#   filter(primers != 'GAPDH') %>%
#   dplyr::rename(TE = ControlMeanCt)
# dCTE <- inner_join(HE, TE, by = c('plateID', 'target', 'shID')) %>%
#   mutate(dCTE = TE-HE)
# 
# HC <- CTC %>%
#   filter(primers == 'GAPDH') %>%
#   dplyr::rename(HC = ControlMeanCt)
# TC <- CTC %>%
#   filter(primers != 'GAPDH') %>%
#   dplyr::rename(TC = ControlMeanCt)
# dCTC <- inner_join(HC, TC, by = c('plateID', 'shID')) %>%
#   mutate(dCTC = TC-HC)
# 
# dCTCavg <- dCTC %>%
#   group_by(plateID, primers.y) %>%
#   summarise(dCTCavg = mean(dCTC))
# 
# ddCT <- left_join(dCTE, dCTCavg, by = c('plateID', 'primers.y')) %>%
#   mutate(ddCT = dCTE - dCTCavg,
#          FC = 2^(-ddCT),
#          percKD = 1-FC) %>%
#   filter(dCTE > 0)
# 
# ddCT$target <- factor(ddCT$target)
# ddCT$target <- factor(ddCT$target, levels = levels(hiFT4a_colonies_sum$target)[levels(hiFT4a_colonies_sum$target) != 'PO'])
# 
# 
# KDbarplot1 <- ggplot() +
#   geom_bar(data = ddCT, aes(shID, percKD), stat = 'identity') +
#   facet_grid(.~target) +
#   ylim(c(0,1)) +
#   theme_classic() +
#   ylab('Average % KD (1 - 2^(-ddCt))') +
#   xlab('shRNA ID for target') +
#   ggtitle('Estimated KD efficiency of shRNAs') +
#   theme(axis.text.y = element_text(size = rel(2)),
#         axis.title = element_blank(),
#         plot.title = element_blank(),
#         strip.text = element_text(size = rel(1.6)))
# ggsave(KDbarplot1, file = '~/Dropbox (RajLab)/Projects/cellid/graphs/hiFT/IAMhiFT4a/KDbarplot_forSlide.pdf', width = 10, height = 5)  
# 
# 
# ## merge the two
# 
# colony_percKD_dotsNbars_forSlide <- ggplot() + 
#   geom_point(data = hiFT4a_colonies %>% filter(target != 'PO') %>% mutate(measurement = 'colonies'), aes(shIndex, colonies), position = position_jitter(w = 0.15, h = 0), alpha = 0.2) +
#   geom_bar(data = hiFT4a_colonies_sum %>% filter(target != 'PO') %>% mutate(measurement = 'colonies'), aes(shIndex, meanCount), alpha = 0.4, width = 0.7, stat = 'identity') +
#   geom_errorbar(data = hiFT4a_colonies_sum %>% filter(target != 'PO') %>% mutate(measurement = 'colonies'), aes(x = shIndex, ymin = errbMin, ymax = errbMax), width = 0.25, alpha = 0.5) +
#   facet_grid(measurement~target, scales = 'free_y') +
#   geom_bar(data = ddCT %>% dplyr::rename(shIndex = shID) %>% mutate(measurement = 'qPCR'), aes(shIndex, percKD), stat = 'identity') +
#   geom_hline(yintercept = controlMeans$meanCtlCount, color = 'blue', alpha = 0.8) +
#   theme_classic() +
#   ylab('Alkaline Phosphatase-positive colonies per well') +
#   xlab('shRNA ID for target') +
#   ggtitle('hiF-T iPSC colony counts after perturbable factor knockdown') +
#   theme(axis.text.y = element_text(size = rel(2)),
#         axis.title = element_blank(),
#         plot.title = element_blank())
# ggsave(colony_percKD_dotsNbars_forSlide, file = '~/Dropbox (RajLab)/Projects/cellid/graphs/hiFT/IAMhiFT4a/colony_percKD_counts_dotsNbars_forSlide.pdf', width = 10, height = 5)  
#  
