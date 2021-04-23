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
# procDataSubdir = 'extractedData/'
# graphSubdir = 'graphs'

procDataDir = paste0(projectDir, 'extractedData/')
graphDir = paste0(projectDir, graphSubdir)

if (!dir.exists(paste0(graphDir, '/hiFT'))){
  dir.create(paste0(graphDir, '/hiFT'))
}
if (!dir.exists(paste0(graphDir, '/hiFT/IAMhiFT3c4c'))){
  dir.create(paste0(graphDir, '/hiFT/IAMhiFT3c4c'))
}

# first, the colony counts
hiFT3c4c_colonies <- as_tibble(read.table(paste0(procDataDir, 'hiFT/IAMhiFT3c4c_ first and second set of shRNAs third rep - Sheet1.tsv'), header = T, stringsAsFactors = F, skip = 103, sep = '\t'))
hiFT3c4c_colonies <- hiFT3c4c_colonies[,1:5]
hits <- c('control', 'ZNF652', 'ZBTB38', 'ATOH8', 'CERS2', 'KLF13', 'CEBPB', 'RUNX1', 'PRRX2')

hiFT3c4c_colonies_hits <- hiFT3c4c_colonies %>% filter(target %in% hits, shID != 'SHC003') %>% 
  mutate(shIndex = ifelse(shID == 'SHC007', 3, shIndex))

hiFT3c4c_colonies_hits_sum <- hiFT3c4c_colonies_hits %>%
  group_by(target, shIndex) %>%
  summarise(meanCount = mean(colonies),
            semCount = sd(colonies)/sqrt(length(colonies)),
            errbMax = meanCount + semCount, 
            errbMin = ifelse(meanCount - semCount > 0, meanCount - semCount, 0))

hiFT3c4c_colonies_hits$target <- factor(hiFT3c4c_colonies_hits$target)
hiFT3c4c_colonies_hits$target <- factor(hiFT3c4c_colonies_hits$target, levels = c('control', 'ZNF652', 'ZBTB38', 'ATOH8', 'CERS2', 'KLF13', 'CEBPB', 'RUNX1', 'PRRX2'))

hiFT3c4c_colonies_hits_sum$target <- factor(hiFT3c4c_colonies_hits_sum$target)
hiFT3c4c_colonies_hits_sum$target <- factor(hiFT3c4c_colonies_hits_sum$target, levels = c('control', 'ZNF652', 'ZBTB38', 'ATOH8', 'CERS2', 'KLF13', 'CEBPB', 'RUNX1', 'PRRX2'))


controlMeans <- hiFT3c4c_colonies_hits %>%
  filter(target == 'control') %>%
  summarise(meanCtlCount = mean(colonies),
            semCtlCount = sd(colonies)/sqrt(length(colonies)))

colony_dotsNbars_forSlide <- ggplot() + 
  geom_point(data = hiFT3c4c_colonies_hits %>% filter(target != 'PO'), aes(as.character(shIndex), colonies), position = position_jitter(w = 0.15, h = 0), alpha = 0.2) +
  geom_bar(data = hiFT3c4c_colonies_hits_sum %>% filter(target != 'PO'), aes(as.character(shIndex), meanCount), alpha = 0.4, width = 0.7, stat = 'identity') +
  geom_errorbar(data = hiFT3c4c_colonies_hits_sum %>% filter(target != 'PO'), aes(x = as.character(shIndex), ymin = errbMin, ymax = errbMax), width = 0.25, alpha = 0.5) +
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
ggsave(colony_dotsNbars_forSlide, file = paste0(graphDir, '/hiFT/IAMhiFT3c4c/colony_counts_dotsNbars_forSlide.pdf'), width = 15, height = 5, useDingbats = F)  

#fisher's exact

targets = unique(as.character(hiFT3c4c_colonies_hits$target))[!unique(as.character(hiFT3c4c_colonies_hits$target)) %in% c('PO', 'control')]

ctlDat <- hiFT3c4c_colonies_hits %>%
  filter(shID %in% c('SHC001', 'SHC002', 'SHC007'))

FET_pvals_3c4c <- list()
for (targ in targets){
  
  tempDat <- hiFT3c4c_colonies_hits %>%
    filter(target == targ)
  
  tempFET <- fisher.test(matrix(c(sum(tempDat$colonies), length(tempDat$colonies)*1e4 - sum(tempDat$colonies),  
                                  sum(ctlDat$colonies), length(ctlDat$colonies)*1e4 - sum(ctlDat$colonies)), ncol = 2), alternative = 'g')
  
  if(is.null(dim(FET_pvals_3c4c))){
    FET_pvals_3c4c = tibble(target = targ,
                          FET_pval = tempFET$p.value)
  } else {
    FET_pvals_3c4c %<>% bind_rows(tibble(target = targ, FET_pval = tempFET$p.value))
  }
  
}

write.csv(FET_pvals_3c4c, file = paste0(graphDir, '/hiFT/IAMhiFT3c4c/FET_pvals_3c4c.csv'))


