#!/usr/bin/env Rscript

# hiFT4b.R loads, analyzes, and plots hiFT4b colony counts and qPCR data.
# 	quality filters:

# don't run lines 22-32 if running from R workspace
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

# if running manually in Rstudio, start here and use:
# projectDir = '~/Dropbox (RajLab)/Shared_IanM/cellid_201807_onward/'
# procDataSubdir = 'extractedData'
# graphSubdir = 'graphs'

procDataDir = paste0(projectDir, 'extractedData')
graphDir = paste0(projectDir, graphSubdir)

if (!dir.exists(paste0(graphDir, '/hiFT'))){
  dir.create(paste0(graphDir, '/hiFT'))
}
if (!dir.exists(paste0(graphDir, '/hiFT/IAMhiFT4b'))){
  dir.create(paste0(graphDir, '/hiFT/IAMhiFT4b'))
}

library(tidyverse)
library(magrittr)

hiFT4b_colonies <- as_tibble(read.table(paste0(procDataDir, '/hiFT/IAMhiFT4b_ second set of shRNAs second rep - Sheet1.tsv'), sep = '\t', header = T, skip = 71))
hiFT4b_colonies <- hiFT4b_colonies[,1:5]

hiFT4b_colonies_sum <- hiFT4b_colonies %>%
  group_by(target, shIndex) %>%
  summarise(meanCount = mean(colonies),
            semCount = sd(colonies)/sqrt(length(colonies)),
            errbMax = meanCount + semCount, 
            errbMin = ifelse(meanCount - semCount > 0, meanCount - semCount, 0))

hiFT4b_colonies_sum_maxes <- hiFT4b_colonies_sum %>%
  group_by(target) %>%
  summarise(maxMeanCount = max(meanCount)) %>%
  arrange(-maxMeanCount) %>%
  filter(target != "control")

hiFT4b_colonies$target <- factor(hiFT4b_colonies$target, levels = c('control', as.character(hiFT4b_colonies_sum_maxes$target)))
hiFT4b_colonies_sum$target <- factor(hiFT4b_colonies_sum$target, levels = c('control', as.character(hiFT4b_colonies_sum_maxes$target)))

controlMeans <- hiFT4b_colonies %>%
  filter(target == 'control') %>%
  summarise(meanCtlCount = mean(colonies),
            semCtlCount = sd(colonies)/sqrt(length(colonies)))

controlMeans_127 <- hiFT4b_colonies %>%
  filter(target == 'control', shID != 'SHC003') %>%
  summarise(meanCtlCount = mean(colonies),
            semCtlCount = sd(colonies)/sqrt(length(colonies)))

hiFT4b_colonies_127 <- hiFT4b_colonies %>%
  filter(shID != 'SHC003') %>%
  mutate(shIndex = ifelse(shID == 'SHC007', 3, shIndex))

hiFT4b_colonies_sum_127 <- hiFT4b_colonies_127 %>%
  group_by(target, shIndex) %>%
  summarise(meanCount = mean(colonies),
            semCount = sd(colonies)/sqrt(length(colonies)),
            errbMax = meanCount + semCount, 
            errbMin = ifelse(meanCount - semCount > 0, meanCount - semCount, 0))

hiFT4b_colonies_sum_maxes_127 <- hiFT4b_colonies_sum_127 %>%
  group_by(target) %>%
  summarise(maxMeanCount = max(meanCount)) %>%
  arrange(-maxMeanCount) %>%
  filter(target != "control")


colony_dotsNbars <- ggplot() + 
  geom_point(data = hiFT4b_colonies %>% filter(target != 'PO'), aes(shIndex, colonies), position = position_jitter(w = 0.15, h = 0), alpha = 0.2) +
  geom_bar(data = hiFT4b_colonies_sum %>% filter(target != 'PO'), aes(shIndex, meanCount), alpha = 0.4, width = 0.7, stat = 'identity') +
  geom_errorbar(data = hiFT4b_colonies_sum %>% filter(target != 'PO'), aes(x = shIndex, ymin = errbMin, ymax = errbMax), width = 0.25, alpha = 0.5) +
  facet_grid(~target) +
  geom_hline(yintercept = controlMeans$meanCtlCount, color = 'blue', alpha = 0.4) +
  theme_classic() +
  ylab('Alkaline Phosphatase-positive colonies per well') +
  xlab('shRNA ID for target') +
  ggtitle('hiF-T iPSC colony counts after perturbable factor knockdown') +
  theme(axis.text.y = element_text(size = rel(2)),
        axis.title = element_text(size = rel(1.5)))
ggsave(colony_dotsNbars, file = paste0(graphDir, '/hiFT/IAMhiFT4b/colony_counts_dotsNbars.pdf'), width = 20, height = 15)  


colony_dotsNbars_forfig <- ggplot() + 
  geom_point(data = hiFT4b_colonies %>% filter(target != 'PO'), aes(shIndex, colonies), position = position_jitter(w = 0.15, h = 0), alpha = 0.2) +
  geom_bar(data = hiFT4b_colonies_sum %>% filter(target != 'PO'), aes(shIndex, meanCount), alpha = 0.4, width = 0.7, stat = 'identity') +
  geom_errorbar(data = hiFT4b_colonies_sum %>% filter(target != 'PO'), aes(x = shIndex, ymin = errbMin, ymax = errbMax), width = 0.25, alpha = 0.5) +
  facet_grid(~target) +
  geom_hline(yintercept = controlMeans$meanCtlCount, color = 'blue', alpha = 0.4) +
  theme_classic() +
  ylab('Alkaline Phosphatase-positive colonies per well') +
  xlab('shRNA ID for target') +
  ggtitle('hiF-T iPSC colony counts after perturbable factor knockdown') +
  theme(axis.text.y = element_text(size = rel(2)),
        axis.title = element_text(size = rel(1.5)))
ggsave(colony_dotsNbars_forfig, file = paste0(graphDir, '/hiFT/IAMhiFT4b/colony_counts_dotsNbars_forfig.pdf'), width = 4, height = 2)  

## t-test comparing each gene against SHC001, 002, 007, the real negative controls

targets = unique(as.character(hiFT4b_colonies$target))[!unique(as.character(hiFT4b_colonies$target)) %in% c('PO', 'control')]

ctlDat <- hiFT4b_colonies %>%
  filter(shID %in% c('SHC001', 'SHC002', 'SHC007'))
for (targ in targets){
  
  tempDat <- hiFT4b_colonies %>%
    filter(target == targ)
  
  tempT <- t.test(tempDat$colonies, ctlDat$colonies, alternative = 'greater')
  
  
}

#fisher's exact
FET_pvals <- list()
for (targ in targets){
  
  tempDat <- hiFT4b_colonies %>%
    filter(target == targ)
  
  tempFET <- fisher.test(matrix(c(sum(tempDat$colonies), length(tempDat$colonies)*1e4 - sum(tempDat$colonies),  
                                sum(ctlDat$colonies), length(ctlDat$colonies)*1e4 - sum(ctlDat$colonies)), ncol = 2), alternative = 'g')
  
  if(is.null(dim(FET_pvals))){
    FET_pvals = tibble(target = targ,
                       FET_pval = tempFET$p.value)
  } else {
    FET_pvals %<>% bind_rows(tibble(target = targ, FET_pval = tempFET$p.value))
  }
  
}

write.csv(FET_pvals, file = paste0(graphDir, '/hiFT/IAMhiFT4b/FET_pvals_4b.csv'))

# qPCR ####

genes <- c('NFATC4', 'YBX3', 'TSC22D1', 'ID1', 'KLF13', 'CERS2', 'LARP1', 'TBX3', 'ZBTB38', 'ATOH8', 'ZNF652', 'NFATC4')
plateFiles <- list.files(paste0(procDataDir, '/hiFT/'))[grepl('hiFT4b_plate', list.files(paste0(procDataDir, '/hiFT/')))]
layoutFiles <- list.files(paste0(procDataDir, '/hiFT/'))[grepl('IAMhiFT4b-qPCR', list.files(paste0(procDataDir, '/hiFT/')))]
conditionMap <- as_tibble(read.table(paste0(procDataDir, '/hiFT/IAMhiFT4b_conditionMap - Sheet1.tsv'), sep = '\t', header = T, stringsAsFactors = F))
targetFileMaster <- as_tibble(read.table(paste0(procDataDir, '/hiFT/IAMhiFT4b qPCR filelist master - Sheet1.tsv'), sep = '\t', header = T, stringsAsFactors = F))
exp <- 'IAMhiFT4b'

# utility function
calculate_knockdown_qPCR <- function(ct_file, layout_file, targ, exludeFlag) { # exclude flag is relevant only to plate 7, with a few edge wells that had a bad seal
  
  ctvals <- as_tibble(read.csv(paste0(procDataDir, '/hiFT/', ct_file),skip = 27, header = T, stringsAsFactors = F)) %>%
    dplyr::select(Well, Ct)
  cDNA <- as_tibble(read.table(paste0(procDataDir, '/hiFT/', layout_file),skip = 1, nrows = 10, sep = '\t', header = T, stringsAsFactors = F)) %>%
    filter(cDNA %in% LETTERS[1:8]) %>%
    dplyr::rename(rowlet = cDNA) %>%
    gather(xcolnum,cDNA,2:13) %>%
    mutate(colnum = sub('.', '', xcolnum)) %>%
    unite(Well, c(rowlet, colnum), sep = '') %>%
    dplyr::select(-xcolnum)
  primers <- as_tibble(read.table(paste0(procDataDir, '/hiFT/', layout_file),skip = 11, nrows = 10, sep = '\t', header = T, stringsAsFactors = F)) %>%
    filter(Primers %in% LETTERS[1:8]) %>%
    dplyr::rename(rowlet = Primers) %>%
    gather(xcolnum,Primers,2:13) %>%
    mutate(colnum = sub('.', '', xcolnum)) %>%
    unite(Well, c(rowlet, colnum), sep = '') %>%
    dplyr::select(-xcolnum)
  results <- inner_join(ctvals, cDNA) %>% inner_join(primers) %>% inner_join(conditionMap, by = 'cDNA') %>%
    filter(Primers %in% c('GAPDH', targ), target %in% c('control', targ))
  
  shcGAPDH <- results %>%
    filter(cDNA == 'SHC002', Primers == 'GAPDH') %>%
    group_by(target, Primers, shID) %>%
    summarise(ControlMeanCt = mean(as.numeric(Ct), na.rm = T),
              ControlSEMCt = sd(as.numeric(Ct), na.rm = T)/sqrt(length(as.numeric(Ct))))
  shcTARG <- results %>%
    filter(cDNA == 'SHC002', Primers == targ) %>%
    group_by(target, Primers, shID) %>%
    summarise(meanCt = mean(as.numeric(Ct), na.rm = T),
              SEMCt = sd(as.numeric(Ct), na.rm = T)/sqrt(length(as.numeric(Ct))))
  
  targGAPDH <- results %>%
    filter(target == targ, Primers == 'GAPDH') %>%
    group_by(target, Primers, shID) %>%
    summarise(ControlMeanCt = mean(as.numeric(Ct), na.rm = T),
              ControlSEMCt = sd(as.numeric(Ct), na.rm = T)/sqrt(length(as.numeric(Ct))))
  if(excludeFlag){
    targGAPDH <- results %>%
      filter(target == targ, Primers == 'GAPDH', !(Well %in% c('F1', 'G1', 'H1'))) %>%
      group_by(target, Primers, shID) %>%
      summarise(ControlMeanCt = mean(as.numeric(Ct), na.rm = T),
                ControlSEMCt = sd(as.numeric(Ct), na.rm = T)/sqrt(length(as.numeric(Ct))))
  }
  targTARG <- results %>%
    filter(target == targ, Primers == targ) %>%
    group_by(target, Primers, shID) %>%
    summarise(meanCt = mean(as.numeric(Ct), na.rm = T),
              SEMCt = sd(as.numeric(Ct), na.rm = T)/sqrt(length(as.numeric(Ct))))
  
  targTARG_perRep <- results %>%
    filter(target == targ, Primers == targ)
  
  shc_dCt <- inner_join(shcGAPDH, shcTARG, by = c('target', 'shID')) %>% mutate(SHC_dCt = meanCt - ControlMeanCt,
                                                                                  SHC_sem_dCt = SEMCt + ControlSEMCt) %>% ungroup() %>%
    mutate(plateID = 1) %>%
    dplyr::select(plateID, SHC_dCt, SHC_sem_dCt)
  
  targ_dCt_ddCt_percKD <- inner_join(targGAPDH, targTARG, by = c('target', 'shID')) %>% mutate(dCt = meanCt - ControlMeanCt,
                                                                                                 sem_dCt = SEMCt + ControlSEMCt) %>% mutate(plateID = 1) %>%
    inner_join(shc_dCt, by = 'plateID') %>%
    mutate(mean_ddCt = dCt - SHC_dCt,
           sem_ddCt = sqrt(SHC_sem_dCt^2 + sem_dCt^2),
           mean_fracKD = 2^(-mean_ddCt) - 1,
           mean_fracKD_lower = 2^(-(mean_ddCt+sem_ddCt)) - 1,
           mean_fracKD_upper = 2^(-(mean_ddCt-sem_ddCt)) - 1) %>% ungroup() %>%
    dplyr::select(target, shID, mean_ddCt, mean_fracKD, mean_fracKD_lower, mean_fracKD_upper)
  
  targ_dCt_ddCt_percKD_perRep <- inner_join(targTARG_perRep, targGAPDH, by = c('target', 'shID')) %>% mutate(dCt = as.numeric(Ct) - ControlMeanCt,
                                                                                                               sem_dCt = ControlSEMCt) %>% mutate(plateID = 1) %>%
    inner_join(shc_dCt, by = 'plateID') %>%
    mutate(ddCt = dCt - SHC_dCt,
           fracKD = 2^(-ddCt) - 1) %>% ungroup() %>%
    dplyr::select(target, shID, ddCt, fracKD)
  
  reslist <- list(
    avgRes = targ_dCt_ddCt_percKD,
    repRes = targ_dCt_ddCt_percKD_perRep
  )
  return(reslist)
}

qPCR_avg_results <- list()
qPCR_rep_results <- list()
for (i in 1:12) {
  
  targ = targetFileMaster$target[i]
  ct_file = targetFileMaster$ct_file[i]
  layout_file = targetFileMaster$layout_file[i]
  excludeFlag = targetFileMaster$excludeFlag[i]
  
  qres = calculate_knockdown_qPCR(ct_file, layout_file, targ, excludeFlag)
  
  temp_avg = qres[['avgRes']]
  temp_rep = qres[['repRes']]
  
  if(is.null(dim(qPCR_avg_results))) {
    qPCR_avg_results <- temp_avg
    qPCR_rep_results <- temp_rep
  } else {
    qPCR_avg_results %<>% bind_rows(temp_avg)
    qPCR_rep_results %<>% bind_rows(temp_rep)
  }
  
}


####
qPCR_avg_results$target <- factor(qPCR_avg_results$target)
qPCR_avg_results$target <- factor(qPCR_avg_results$target, levels = levels(hiFT4b_colonies_sum$target)[levels(hiFT4b_colonies_sum$target) != 'control'])

qPCR_rep_results$target <- factor(qPCR_rep_results$target)
qPCR_rep_results$target <- factor(qPCR_rep_results$target, levels = levels(hiFT4b_colonies_sum$target)[levels(hiFT4b_colonies_sum$target) != 'control'])

colony_percKD_dotsNbars_both_bothErr_forSlide <- ggplot() + 
  geom_point(data = hiFT4b_colonies %>% filter(target != 'PO') %>% mutate(measurement = 'colonies'), aes(shIndex, colonies), position = position_jitter(w = 0.15, h = 0), alpha = 0.2) +
  geom_bar(data = hiFT4b_colonies_sum %>% filter(target != 'PO') %>% mutate(measurement = 'colonies'), aes(shIndex, meanCount), alpha = 0.4, width = 0.7, stat = 'identity') +
  geom_errorbar(data = hiFT4b_colonies_sum %>% filter(target != 'PO') %>% mutate(measurement = 'colonies'), aes(x = shIndex, ymin = errbMin, ymax = errbMax), width = 0.25, alpha = 0.5) +
  facet_grid(measurement~target, scales = 'free_y') +
  geom_bar(data = qPCR_avg_results %>% dplyr::rename(shIndex = shID) %>% mutate(measurement = 'qPCR'), aes(shIndex, mean_fracKD), stat = 'identity', alpha = 0.4) +
  geom_errorbar(data = qPCR_avg_results %>% dplyr::rename(shIndex = shID) %>% mutate(measurement = 'qPCR'), aes(x = shIndex, ymin = mean_fracKD_lower, ymax = mean_fracKD_upper), width = 0.25, alpha = 0.5) +
  geom_point(data = qPCR_rep_results %>% dplyr::rename(shIndex = shID) %>% mutate(measurement = 'qPCR'), aes(shIndex, fracKD), position = position_jitter(w = 0.15, h = 0), alpha = 0.2) +
  geom_hline(data = controlMeans %>% mutate(measurement = 'colonies'), aes(yintercept = meanCtlCount), color = 'blue', alpha = 0.8) +
  theme_classic() +
  ylab('Alkaline Phosphatase-positive colonies per well') +
  xlab('shRNA ID for target') +
  ggtitle('hiF-T iPSC colony counts after perturbable factor knockdown') +
  theme(axis.text.y = element_text(size = rel(2)),
        axis.title = element_blank(),
        plot.title = element_blank())
ggsave(colony_percKD_dotsNbars_both_bothErr_forSlide, file = paste0(graphDir, '/hiFT/IAMhiFT4b/colony_percKD_dotsNbars_both_bothErr_forSlide.pdf'), width = 10, height = 5)  

hits <- c('control', 'ZNF652', 'ZBTB38', 'ATOH8', 'CERS2', 'KLF13')
colony_percKD_dotsNbars_both_bothErr_hitsOnly <- ggplot() + 
  geom_point(data = hiFT4b_colonies %>% filter(target != 'PO') %>% mutate(measurement = 'colonies') %>% filter(target %in% hits), aes(shIndex, colonies), position = position_jitter(w = 0.15, h = 0), alpha = 0.2) +
  geom_bar(data = hiFT4b_colonies_sum %>% filter(target != 'PO') %>% mutate(measurement = 'colonies') %>% filter(target %in% hits), aes(shIndex, meanCount), alpha = 0.4, width = 0.7, stat = 'identity') +
  geom_errorbar(data = hiFT4b_colonies_sum %>% filter(target != 'PO') %>% mutate(measurement = 'colonies') %>% filter(target %in% hits), aes(x = shIndex, ymin = errbMin, ymax = errbMax), width = 0.25, alpha = 0.5) +
  facet_grid(measurement~target, scales = 'free_y') +
  geom_bar(data = qPCR_avg_results %>% dplyr::rename(shIndex = shID) %>% mutate(measurement = 'qPCR') %>% filter(target %in% hits), aes(shIndex, mean_fracKD), stat = 'identity', alpha = 0.4) +
  geom_errorbar(data = qPCR_avg_results %>% dplyr::rename(shIndex = shID) %>% mutate(measurement = 'qPCR') %>% filter(target %in% hits), aes(x = shIndex, ymin = mean_fracKD_lower, ymax = mean_fracKD_upper), width = 0.25, alpha = 0.5) +
  geom_point(data = qPCR_rep_results %>% dplyr::rename(shIndex = shID) %>% mutate(measurement = 'qPCR') %>% filter(target %in% hits), aes(shIndex, fracKD), position = position_jitter(w = 0.15, h = 0), alpha = 0.2) +
  geom_hline(data = controlMeans %>% mutate(measurement = 'colonies'), aes(yintercept = meanCtlCount), color = 'blue', alpha = 0.8) +
  theme_classic() +
  ylab('Alkaline Phosphatase-positive colonies per well') +
  xlab('shRNA ID for target') +
  ggtitle('hiF-T iPSC colony counts after perturbable factor knockdown') +
  theme(axis.text.y = element_text(size = rel(2)),
        axis.title = element_blank(),
        plot.title = element_blank())
ggsave(colony_percKD_dotsNbars_both_bothErr_hitsOnly, file = paste0(graphDir, '/hiFT/IAMhiFT4b/colony_percKD_dotsNbars_both_bothErr_hitsOnly.pdf'), width = 10, height = 5)  


colony_dotsNbars_hitsOnly <- ggplot() + 
  geom_point(data = hiFT4b_colonies %>% filter(target != 'PO') %>% mutate(measurement = 'colonies') %>% filter(target %in% hits), aes(shIndex, colonies), position = position_jitter(w = 0.15, h = 0), alpha = 0.2) +
  geom_bar(data = hiFT4b_colonies_sum %>% filter(target != 'PO') %>% mutate(measurement = 'colonies') %>% filter(target %in% hits), aes(shIndex, meanCount), alpha = 0.4, width = 0.7, stat = 'identity') +
  geom_errorbar(data = hiFT4b_colonies_sum %>% filter(target != 'PO') %>% mutate(measurement = 'colonies') %>% filter(target %in% hits), aes(x = shIndex, ymin = errbMin, ymax = errbMax), width = 0.25, alpha = 0.5) +
  facet_grid(.~target) +
  # geom_bar(data = qPCR_avg_results %>% dplyr::rename(shIndex = shID) %>% mutate(measurement = 'qPCR') %>% filter(target %in% hits), aes(shIndex, mean_fracKD), stat = 'identity', alpha = 0.4) +
  # geom_errorbar(data = qPCR_avg_results %>% dplyr::rename(shIndex = shID) %>% mutate(measurement = 'qPCR') %>% filter(target %in% hits), aes(x = shIndex, ymin = mean_fracKD_lower, ymax = mean_fracKD_upper), width = 0.25, alpha = 0.5) +
  # geom_point(data = qPCR_rep_results %>% dplyr::rename(shIndex = shID) %>% mutate(measurement = 'qPCR') %>% filter(target %in% hits), aes(shIndex, fracKD), position = position_jitter(w = 0.15, h = 0), alpha = 0.2) +
  geom_hline(data = controlMeans %>% mutate(measurement = 'colonies'), aes(yintercept = meanCtlCount), color = 'blue', alpha = 0.8) +
  theme_classic() +
  ylab('Alkaline Phosphatase-positive colonies per well') +
  xlab('shRNA ID for target') +
  ggtitle('hiF-T iPSC colony counts after perturbable factor knockdown') +
  theme(axis.text.y = element_text(size = rel(2)),
        axis.title = element_blank(),
        plot.title = element_blank())
ggsave(colony_dotsNbars_hitsOnly, file = paste0(graphDir, '/hiFT/IAMhiFT4b/colony_dotsNbars_hitsOnly.pdf'), width = 10, height = 5)  

colony_dotsNbars_hitsOnly_forfig <- ggplot() + 
  geom_point(data = hiFT4b_colonies_127 %>% filter(target != 'PO') %>% mutate(measurement = 'colonies') %>% filter(target %in% hits), aes(as.character(shIndex), colonies), position = position_jitter(w = 0.15, h = 0), alpha = 0.2) +
  geom_bar(data = hiFT4b_colonies_sum_127 %>% filter(target != 'PO') %>% mutate(measurement = 'colonies') %>% filter(target %in% hits), aes(as.character(shIndex), meanCount), alpha = 0.4, width = 0.7, stat = 'identity') +
  geom_errorbar(data = hiFT4b_colonies_sum_127%>% filter(target != 'PO') %>% mutate(measurement = 'colonies') %>% filter(target %in% hits), aes(x = as.character(shIndex), ymin = errbMin, ymax = errbMax), width = 0.25, alpha = 0.5) +
  facet_grid(.~target) +
  # geom_bar(data = qPCR_avg_results %>% dplyr::rename(shIndex = shID) %>% mutate(measurement = 'qPCR') %>% filter(target %in% hits), aes(shIndex, mean_fracKD), stat = 'identity', alpha = 0.4) +
  # geom_errorbar(data = qPCR_avg_results %>% dplyr::rename(shIndex = shID) %>% mutate(measurement = 'qPCR') %>% filter(target %in% hits), aes(x = shIndex, ymin = mean_fracKD_lower, ymax = mean_fracKD_upper), width = 0.25, alpha = 0.5) +
  # geom_point(data = qPCR_rep_results %>% dplyr::rename(shIndex = shID) %>% mutate(measurement = 'qPCR') %>% filter(target %in% hits), aes(shIndex, fracKD), position = position_jitter(w = 0.15, h = 0), alpha = 0.2) +
  geom_hline(data = controlMeans_127 %>% mutate(measurement = 'colonies'), aes(yintercept = meanCtlCount), color = 'blue', alpha = 0.8) +
  theme_classic() +
  ylab('Alkaline Phosphatase-positive colonies per well') +
  xlab('shRNA ID for target') +
  ggtitle('hiF-T iPSC colony counts after perturbable factor knockdown') +
  theme(axis.text.y = element_text(size = rel(2)),
        axis.title = element_blank(),
        plot.title = element_blank(),
        panel.spacing = unit(0.3, "lines"))
ggsave(colony_dotsNbars_hitsOnly_forfig, file = paste0(graphDir, '/hiFT/IAMhiFT4b/colony_dotsNbars_hitsOnly_forfig.pdf'), width = 4, height = 2)  

colony_percKD_dotsNbars_both_bothErr_hitsOnly_forfig <- ggplot() + 
  geom_point(data = hiFT4b_colonies_127 %>% filter(target != 'PO') %>% mutate(measurement = 'colonies') %>% filter(target %in% hits), aes(shIndex, colonies), position = position_jitter(w = 0.15, h = 0), alpha = 0.2) +
  geom_bar(data = hiFT4b_colonies_sum_127 %>% filter(target != 'PO') %>% mutate(measurement = 'colonies') %>% filter(target %in% hits), aes(shIndex, meanCount), alpha = 0.4, width = 0.7, stat = 'identity') +
  geom_errorbar(data = hiFT4b_colonies_sum_127 %>% filter(target != 'PO') %>% mutate(measurement = 'colonies') %>% filter(target %in% hits), aes(x = shIndex, ymin = errbMin, ymax = errbMax), width = 0.25, alpha = 0.5) +
  facet_grid(measurement~target, scales = 'free_y') +
  geom_bar(data = qPCR_avg_results %>% dplyr::rename(shIndex = shID) %>% mutate(measurement = 'qPCR') %>% filter(target %in% hits), aes(shIndex, mean_fracKD), stat = 'identity', alpha = 0.4) +
  geom_errorbar(data = qPCR_avg_results %>% dplyr::rename(shIndex = shID) %>% mutate(measurement = 'qPCR') %>% filter(target %in% hits), aes(x = shIndex, ymin = mean_fracKD_lower, ymax = mean_fracKD_upper), width = 0.25, alpha = 0.5) +
  geom_point(data = qPCR_rep_results %>% dplyr::rename(shIndex = shID) %>% mutate(measurement = 'qPCR') %>% filter(target %in% hits), aes(shIndex, fracKD), position = position_jitter(w = 0.15, h = 0), alpha = 0.2) +
  geom_hline(data = controlMeans_127 %>% mutate(measurement = 'colonies'), aes(yintercept = meanCtlCount), color = 'blue', alpha = 0.8) +
  theme_classic() +
  ylab('Alkaline Phosphatase-positive colonies per well') +
  xlab('shRNA ID for target') +
  ggtitle('hiF-T iPSC colony counts after perturbable factor knockdown') +
  theme(axis.text.y = element_text(size = rel(2)),
        axis.title = element_blank(),
        plot.title = element_blank())
ggsave(colony_percKD_dotsNbars_both_bothErr_hitsOnly_forfig, file = paste0(graphDir, '/hiFT/IAMhiFT4b/colony_percKD_dotsNbars_both_bothErr_hitsOnly_forfig.pdf'), width = 10, height = 5)  


# 
# ctvals <- as_tibble(read.csv(paste0(procDataDir, '/hiFT/hiFT4b_plate1_ZNF652_NFATC4_20190728.csv'),skip = 27, header = T, stringsAsFactors = F)) %>%
#   dplyr::select(Well, Ct)
# cDNA <- as_tibble(read.table(paste0(procDataDir, '/hiFT/IAMhiFT4b-qPCRplate1 20190728 - Sheet1.tsv'),skip = 1, nrows = 10, sep = '\t', header = T, stringsAsFactors = F)) %>%
#   filter(cDNA %in% LETTERS[1:8]) %>%
#   dplyr::rename(rowlet = cDNA) %>%
#   gather(xcolnum,cDNA,2:13) %>%
#   mutate(colnum = sub('.', '', xcolnum)) %>%
#   unite(Well, c(rowlet, colnum), sep = '') %>%
#   dplyr::select(-xcolnum)
# primers <- as_tibble(read.table(paste0(procDataDir, '/hiFT/IAMhiFT4b-qPCRplate1 20190728 - Sheet1.tsv'),skip = 11, nrows = 10, sep = '\t', header = T, stringsAsFactors = F)) %>%
#   filter(Primers %in% LETTERS[1:8]) %>%
#   dplyr::rename(rowlet = Primers) %>%
#   gather(xcolnum,Primers,2:13) %>%
#   mutate(colnum = sub('.', '', xcolnum)) %>%
#   unite(Well, c(rowlet, colnum), sep = '') %>%
#   dplyr::select(-xcolnum)
# results <- inner_join(ctvals, cDNA) %>% inner_join(primers) %>% inner_join(conditionMap, by = 'cDNA') %>%
#   filter(Primers %in% c('GAPDH', 'NFATC4'), target %in% c('control', 'NFATC4'))
# 
# shcGAPDH <- results %>%
#   filter(cDNA == 'SHC002', Primers == 'GAPDH') %>%
#   group_by(target, Primers, shID) %>%
#   summarise(ControlMeanCt = mean(as.numeric(Ct), na.rm = T),
#             ControlSEMCt = sd(as.numeric(Ct), na.rm = T)/sqrt(length(as.numeric(Ct))))
# shcNFATC4 <- results %>%
#   filter(cDNA == 'SHC002', Primers == 'NFATC4') %>%
#   group_by(target, Primers, shID) %>%
#   summarise(meanCt = mean(as.numeric(Ct), na.rm = T),
#             SEMCt = sd(as.numeric(Ct), na.rm = T)/sqrt(length(as.numeric(Ct))))
# 
# nfatGAPDH <- results %>%
#   filter(target == 'NFATC4', Primers == 'GAPDH') %>%
#   group_by(target, Primers, shID) %>%
#   summarise(ControlMeanCt = mean(as.numeric(Ct), na.rm = T),
#             ControlSEMCt = sd(as.numeric(Ct), na.rm = T)/sqrt(length(as.numeric(Ct))))
# nfatNFATC4 <- results %>%
#   filter(target == 'NFATC4', Primers == 'NFATC4') %>%
#   group_by(target, Primers, shID) %>%
#   summarise(meanCt = mean(as.numeric(Ct), na.rm = T),
#             SEMCt = sd(as.numeric(Ct), na.rm = T)/sqrt(length(as.numeric(Ct))))
# 
# nfatNFATC4_perRep <- results %>%
#   filter(target == 'NFATC4', Primers == 'NFATC4')
# 
# shc_dCt <- inner_join(shcGAPDH, shcNFATC4, by = c('target', 'shID')) %>% mutate(SHC_dCt = meanCt - ControlMeanCt,
#                                                                                 SHC_sem_dCt = SEMCt + ControlSEMCt) %>% ungroup() %>%
#   mutate(plateID = 1) %>%
#   dplyr::select(plateID, SHC_dCt, SHC_sem_dCt)
# 
# nfat_dCt_ddCt_percKD <- inner_join(nfatGAPDH, nfatNFATC4, by = c('target', 'shID')) %>% mutate(dCt = meanCt - ControlMeanCt,
#                                                                                                sem_dCt = SEMCt + ControlSEMCt) %>% mutate(plateID = 1) %>%
#   inner_join(shc_dCt, by = 'plateID') %>%
#   mutate(mean_ddCt = dCt - SHC_dCt,
#          sem_ddCt = sqrt(SHC_sem_dCt^2 + sem_dCt^2),
#          mean_fracKD = 2^(-mean_ddCt) - 1,
#          mean_fracKD_lower = 2^(-(mean_ddCt+sem_ddCt)) - 1,
#          mean_fracKD_upper = 2^(-(mean_ddCt-sem_ddCt)) - 1) %>% ungroup() %>%
#   dplyr::select(target, shID, mean_ddCt, mean_fracKD, mean_fracKD_lower, mean_fracKD_upper)
# 
# nfat_dCt_ddCt_percKD_perRep <- inner_join(nfatNFATC4_perRep, nfatGAPDH, by = c('target', 'shID')) %>% mutate(dCt = as.numeric(Ct) - ControlMeanCt,
#                                                                                                              sem_dCt = ControlSEMCt) %>% mutate(plateID = 1) %>%
#   inner_join(shc_dCt, by = 'plateID') %>%
#   mutate(ddCt = dCt - SHC_dCt,
#          fracKD = 2^(-ddCt) - 1) %>% ungroup() %>%
#   dplyr::select(target, shID, ddCt, fracKD)
# 
# hiFT4b_KD_results_avg <- nfat_dCt_ddCt_percKD
# hiFT4b_KD_results_perRep <- nfat_dCt_ddCt_percKD_perRep
# 
# # YBX3 and SKIL
# ctvals <- as_tibble(read.csv(paste0(procDataDir, '/hiFT/hiFT4b_plate6_SKIL_YBX3_20190802.csv'),skip = 27, header = T, stringsAsFactors = F)) %>%
#   dplyr::select(Well, Ct)
# cDNA <- as_tibble(read.table(paste0(procDataDir, '/hiFT/IAMhiFT4b-qPCRplate6 20190802 - Sheet1.tsv'),skip = 1, nrows = 10, sep = '\t', header = T, stringsAsFactors = F)) %>%
#   filter(cDNA %in% LETTERS[1:8]) %>%
#   dplyr::rename(rowlet = cDNA) %>%
#   gather(xcolnum,cDNA,2:13) %>%
#   mutate(colnum = sub('.', '', xcolnum)) %>%
#   unite(Well, c(rowlet, colnum), sep = '') %>%
#   dplyr::select(-xcolnum)
# primers <- as_tibble(read.table(paste0(procDataDir, '/hiFT/IAMhiFT4b-qPCRplate6 20190802 - Sheet1.tsv'),skip = 11, nrows = 10, sep = '\t', header = T, stringsAsFactors = F)) %>%
#   filter(Primers %in% LETTERS[1:8]) %>%
#   dplyr::rename(rowlet = Primers) %>%
#   gather(xcolnum,Primers,2:13) %>%
#   mutate(colnum = sub('.', '', xcolnum)) %>%
#   unite(Well, c(rowlet, colnum), sep = '') %>%
#   dplyr::select(-xcolnum)
# results <- inner_join(ctvals, cDNA) %>% inner_join(primers) %>% inner_join(conditionMap, by = 'cDNA') %>%
#   filter(Primers %in% c('GAPDH', 'SKIL', 'YBX3'), target %in% c('control', 'SKIL', 'YBX3'))
# 
# shcGAPDH <- results %>%
#   filter(cDNA == 'SHC002', Primers == 'GAPDH') %>%
#   group_by(target, Primers, shID) %>%
#   summarise(ControlMeanCt = mean(as.numeric(Ct), na.rm = T),
#             ControlSEMCt = sd(as.numeric(Ct), na.rm = T)/sqrt(length(as.numeric(Ct))))
# shcSKIL <- results %>%
#   filter(cDNA == 'SHC002', Primers == 'SKIL') %>%
#   group_by(target, Primers, shID) %>%
#   summarise(meanCt = mean(as.numeric(Ct), na.rm = T),
#             SEMCt = sd(as.numeric(Ct), na.rm = T)/sqrt(length(as.numeric(Ct))))
# shcYBX3 <- results %>%
#   filter(cDNA == 'SHC002', Primers == 'YBX3') %>%
#   group_by(target, Primers, shID) %>%
#   summarise(meanCt = mean(as.numeric(Ct), na.rm = T),
#             SEMCt = sd(as.numeric(Ct), na.rm = T)/sqrt(length(as.numeric(Ct))))
# 
# skilGAPDH <- results %>%
#   filter(target == 'SKIL', Primers == 'GAPDH') %>%
#   group_by(target, Primers, shID) %>%
#   summarise(ControlMeanCt = mean(as.numeric(Ct), na.rm = T),
#             ControlSEMCt = sd(as.numeric(Ct), na.rm = T)/sqrt(length(as.numeric(Ct))))
# skilSKIL <- results %>%
#   filter(target == 'SKIL', Primers == 'SKIL') %>%
#   group_by(target, Primers, shID) %>%
#   summarise(meanCt = mean(as.numeric(Ct), na.rm = T),
#             SEMCt = sd(as.numeric(Ct), na.rm = T)/sqrt(length(as.numeric(Ct))))
# 
# skilSKIL_perRep <- results %>%
#   filter(target == 'SKIL', Primers == 'SKIL')
# 
# 
# ybx3GAPDH <- results %>%
#   filter(target == 'YBX3', Primers == 'GAPDH') %>%
#   group_by(target, Primers, shID) %>%
#   summarise(ControlMeanCt = mean(as.numeric(Ct), na.rm = T),
#             ControlSEMCt = sd(as.numeric(Ct), na.rm = T)/sqrt(length(as.numeric(Ct))))
# ybx3YBX3 <- results %>%
#   filter(target == 'YBX3', Primers == 'YBX3') %>%
#   group_by(target, Primers, shID) %>%
#   summarise(meanCt = mean(as.numeric(Ct), na.rm = T),
#             SEMCt = sd(as.numeric(Ct), na.rm = T)/sqrt(length(as.numeric(Ct))))
# 
# ybx3YBX3_perRep <- results %>%
#   filter(target == 'YBX3', Primers == 'YBX3')
# 
# shc_dCt <- inner_join(shcGAPDH, shcSKIL, by = c('target', 'shID')) %>% mutate(SHC_dCt = meanCt - ControlMeanCt,
#                                                                                 SHC_sem_dCt = SEMCt + ControlSEMCt) %>% ungroup() %>%
#   mutate(plateID = 1) %>%
#   dplyr::select(plateID, SHC_dCt, SHC_sem_dCt)
# 
# skil_dCt_ddCt_percKD <- inner_join(skilGAPDH, skilSKIL, by = c('target', 'shID')) %>% mutate(dCt = meanCt - ControlMeanCt,
#                                                                                                sem_dCt = SEMCt + ControlSEMCt) %>% mutate(plateID = 6) %>%
#   inner_join(shc_dCt, by = 'plateID') %>%
#   mutate(mean_ddCt = dCt - SHC_dCt,
#          sem_ddCt = sqrt(SHC_sem_dCt^2 + sem_dCt^2),
#          mean_fracKD = 2^(-mean_ddCt) - 1,
#          mean_fracKD_lower = 2^(-(mean_ddCt+sem_ddCt)) - 1,
#          mean_fracKD_upper = 2^(-(mean_ddCt-sem_ddCt)) - 1) %>% ungroup() %>%
#   dplyr::select(target, shID, mean_ddCt, mean_fracKD, mean_fracKD_lower, mean_fracKD_upper)
# 
# skil_dCt_ddCt_percKD_perRep <- inner_join(skilSKIL_perRep, skilGAPDH, by = c('target', 'shID')) %>% mutate(dCt = as.numeric(Ct) - ControlMeanCt,
#                                                                                                              sem_dCt = ControlSEMCt) %>% mutate(plateID = 6) %>%
#   inner_join(shc_dCt, by = 'plateID') %>%
#   mutate(ddCt = dCt - SHC_dCt,
#          fracKD = 2^(-ddCt) - 1) %>% ungroup() %>%
#   dplyr::select(target, shID, ddCt, fracKD)
# 
# ybx3_dCt_ddCt_percKD <- inner_join(ybx3GAPDH, ybx3YBX3, by = c('target', 'shID')) %>% mutate(dCt = meanCt - ControlMeanCt,
#                                                                                                sem_dCt = SEMCt + ControlSEMCt) %>% mutate(plateID = 6) %>%
#   inner_join(shc_dCt, by = 'plateID') %>%
#   mutate(mean_ddCt = dCt - SHC_dCt,
#          sem_ddCt = sqrt(SHC_sem_dCt^2 + sem_dCt^2),
#          mean_fracKD = 2^(-mean_ddCt) - 1,
#          mean_fracKD_lower = 2^(-(mean_ddCt+sem_ddCt)) - 1,
#          mean_fracKD_upper = 2^(-(mean_ddCt-sem_ddCt)) - 1) %>% ungroup() %>%
#   dplyr::select(target, shID, mean_ddCt, mean_fracKD, mean_fracKD_lower, mean_fracKD_upper)
# 
# ybx3_dCt_ddCt_percKD_perRep <- inner_join(ybx3YBX3_perRep, ybx3GAPDH, by = c('target', 'shID')) %>% mutate(dCt = as.numeric(Ct) - ControlMeanCt,
#                                                                                                              sem_dCt = ControlSEMCt) %>% mutate(plateID = 6) %>%
#   inner_join(shc_dCt, by = 'plateID') %>%
#   mutate(ddCt = dCt - SHC_dCt,
#          fracKD = 2^(-ddCt) - 1) %>% ungroup() %>%
#   dplyr::select(target, shID, ddCt, fracKD)
# 






