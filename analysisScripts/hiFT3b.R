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

library(tidyverse)
library(magrittr)

# projectDir = '~/Dropbox (RajLab)/Shared_IanM/cellid_201807_onward/'
# procDataSubdir = 'extractedData'
# graphSubdir = 'graphs'

procDataDir = paste0(projectDir, 'extractedData/')
graphDir = paste0(projectDir, graphSubdir)

if (!dir.exists(paste0(graphDir, '/hiFT'))){
  dir.create(paste0(graphDir, '/hiFT'))
}
if (!dir.exists(paste0(graphDir, '/hiFT/IAMhiFT3b'))){
  dir.create(paste0(graphDir, '/hiFT/IAMhiFT3b'))
}

# # first, the colony counts
hiFT3b_colonies <- as_tibble(read.csv(paste0(procDataDir, 'hiFT/IAMhiFT3b_ First pass set - candidate shRNAs with controls - 3 sh per target - Sheet1.csv'), header = T, stringsAsFactors = F, skip = 110))

hiFT3b_colonies_sum <- hiFT3b_colonies %>%
  group_by(target, shIndex) %>%
  summarise(meanCount = mean(colonies),
            semCount = sd(colonies)/sqrt(length(colonies)),
            errbMax = meanCount + semCount,
            errbMin = ifelse(meanCount - semCount > 0, meanCount - semCount, 0))

hiFT3b_colonies$target <- factor(hiFT3b_colonies$target)
hiFT3b_colonies$target <- factor(hiFT3b_colonies$target, levels = c('control', 'CEBPB', 'PRRX2', 'RUNX1', 'ID3'))

hiFT3b_colonies_sum$target <- factor(hiFT3b_colonies_sum$target)
hiFT3b_colonies_sum$target <- factor(hiFT3b_colonies_sum$target, levels = c('control', 'CEBPB', 'PRRX2', 'RUNX1', 'ID3'))


controlMeans <- hiFT3b_colonies %>%
  filter(target == 'control') %>%
  summarise(meanCtlCount = mean(colonies),
            semCtlCount = sd(colonies)/sqrt(length(colonies)))

# colony_dotsNbars_forSlide <- ggplot() + 
#   geom_point(data = hiFT3b_colonies %>% filter(target != 'PO'), aes(shIndex, colonies), position = position_jitter(w = 0.15, h = 0), alpha = 0.2) +
#   geom_bar(data = hiFT3b_colonies_sum %>% filter(target != 'PO'), aes(shIndex, meanCount), alpha = 0.4, width = 0.7, stat = 'identity') +
#   geom_errorbar(data = hiFT3b_colonies_sum %>% filter(target != 'PO'), aes(x = shIndex, ymin = errbMin, ymax = errbMax), width = 0.25, alpha = 0.5) +
#   facet_grid(~target) +
#   geom_hline(yintercept = controlMeans$meanCtlCount, color = 'blue', alpha = 0.8) +
#   theme_classic() +
#   ylab('Alkaline Phosphatase-positive colonies per well') +
#   xlab('shRNA ID for target') +
#   ggtitle('hiF-T iPSC colony counts after perturbable factor knockdown') +
#   theme(axis.text.y = element_text(size = rel(2)),
#         axis.title = element_blank(),
#         plot.title = element_blank(),
#         strip.text = element_text(size = rel(1.6)))
# ggsave(colony_dotsNbars_forSlide, file = paste0(graphDir,'/hiFT/IAMhiFT3b/colony_counts_dotsNbars_forSlide.pdf'), width = 10, height = 5)  
# 
# 
# first, the colony counts
hiFT3b_colonies_127 <- as_tibble(read.csv(paste0(procDataDir, 'hiFT/IAMhiFT3b_ First pass set - candidate shRNAs with controls - 3 sh per target - Sheet1.csv'), header = T, stringsAsFactors = F, skip = 110))%>%
  filter(shID != 'SCH003',
         target != 'ID3') %>%
  mutate(shIndex = ifelse(shID == 'SHC007', 3, shIndex))

hiFT3b_colonies_sum_127 <- hiFT3b_colonies_127 %>%
  group_by(target, shIndex) %>%
  summarise(meanCount = mean(colonies),
            semCount = sd(colonies)/sqrt(length(colonies)),
            errbMax = meanCount + semCount,
            errbMin = ifelse(meanCount - semCount > 0, meanCount - semCount, 0))

hiFT3b_colonies_127$target <- factor(hiFT3b_colonies_127$target)
hiFT3b_colonies_127$target <- factor(hiFT3b_colonies_127$target, levels = c('control', 'CEBPB', 'RUNX1', 'PRRX2'))

hiFT3b_colonies_sum_127$target <- factor(hiFT3b_colonies_sum_127$target)
hiFT3b_colonies_sum_127$target <- factor(hiFT3b_colonies_sum_127$target, levels = c('control', 'CEBPB', 'RUNX1', 'PRRX2'))


controlMeans_127 <- hiFT3b_colonies_127 %>%
  filter(target == 'control') %>%
  summarise(meanCtlCount = mean(colonies),
            semCtlCount = sd(colonies)/sqrt(length(colonies)))

colony_dotsNbars_forSlide_127 <- ggplot() +
  geom_point(data = hiFT3b_colonies_127 %>% filter(target != 'PO'), aes(as.character(shIndex), colonies), position = position_jitter(w = 0.15, h = 0), alpha = 0.2) +
  geom_bar(data = hiFT3b_colonies_sum_127 %>% filter(target != 'PO'), aes(as.character(shIndex), meanCount), alpha = 0.4, width = 0.7, stat = 'identity') +
  geom_errorbar(data = hiFT3b_colonies_sum_127 %>% filter(target != 'PO'), aes(x = as.character(shIndex), ymin = errbMin, ymax = errbMax), width = 0.25, alpha = 0.5) +
  facet_grid(~target) +
  geom_hline(yintercept = controlMeans_127$meanCtlCount, color = 'blue', alpha = 0.8) +
  theme_classic() +
  ylab('Alkaline Phosphatase-positive colonies per well') +
  xlab('shRNA ID for target') +
  ggtitle('hiF-T iPSC colony counts after perturbable factor knockdown') +
  theme(axis.text.y = element_text(size = rel(2)),
        axis.title = element_blank(),
        plot.title = element_blank(),
        # strip.text = element_text(size = rel(1.6)),
        panel.spacing = unit(0.3, "lines"))
ggsave(colony_dotsNbars_forSlide_127, file = paste0(graphDir, '/hiFT/IAMhiFT3b/colony_counts_dotsNbars_127.pdf'), width = 8/3, height = 2)
# 





# # qPCR ####
# 
# genes <- c('CEBPB', 'ID3', 'PRRX2', 'RUNX1')
# plates <- c(1,2)
# exp <- 'IAMhiFT3b'
# 
# all_results<-list()
# for (plate in plates) {
#   
#   genesInPlate = genes[c(2*(plate-1) + 1, 2*(plate-1) + 2)]
#   
#   results <- as_tibble(read.csv(paste0(procDataDir, exp, '-qPCRplate', as.character(plate), '_somecols - Sheet1.csv'), header = T, stringsAsFactors = F))
#   
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
#   filter(target == 'control', !is.na(shID)) %>%
#   group_by(plateID, primers, shID) %>%
#   summarise(ControlMeanCt = mean(as.numeric(Ct), na.rm = T))
#   
# CTE <- all_results %>%
#   filter(target != 'control', !is.na(shID)) %>%
#   group_by(plateID, target, primers, shID) %>%
#   summarise(ControlMeanCt = mean(as.numeric(Ct), na.rm = T))
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
#   # correct for shift caused by having only 2 PRRX2 samples
#   ungroup() %>%
#   mutate(target = ifelse(paste0(target, shID) == 'PRRX23', 'RUNX1', target),
#          shID = ifelse(paste0(target, primers.y) == 'RUNX1PRRX2', 1, 
#                        ifelse(paste0(target,shID,primers.y) == 'RUNX11RUNX1', 2,
#                               ifelse(paste0(target,shID,primers.y) == 'RUNX12RUNX1', 3, shID)))) %>%
#   filter(dCTE > 0)
# 
# ddCT$target <- factor(ddCT$target)
# ddCT$target <- factor(ddCT$target, levels = c('CEBPB', 'PRRX2', 'RUNX1', 'ID3'))
# 
# 
# qPCR_avg_results <- list()
# qPCR_rep_results <- list()
# for (i in 1:12) {
#   
#   targ = targetFileMaster$target[i]
#   ct_file = targetFileMaster$ct_file[i]
#   layout_file = targetFileMaster$layout_file[i]
#   excludeFlag = targetFileMaster$excludeFlag[i]
#   
#   qres = calculate_knockdown_qPCR(ct_file, layout_file, targ, excludeFlag)
#   
#   temp_avg = qres[['avgRes']]
#   temp_rep = qres[['repRes']]
#   
#   if(is.null(dim(qPCR_avg_results))) {
#     qPCR_avg_results <- temp_avg
#     qPCR_rep_results <- temp_rep
#   } else {
#     qPCR_avg_results %<>% bind_rows(temp_avg)
#     qPCR_rep_results %<>% bind_rows(temp_rep)
#   }
#   
# }

### updated format
genes <- c('CEBPB', 'ID3', 'PRRX2', 'RUNX1')
plateFiles <- list.files(paste0(procDataDir, '/hiFT/'))[grepl('hiFT3b_q', list.files(paste0(procDataDir, '/hiFT/')))]
layoutFiles <- list.files(paste0(procDataDir, '/hiFT/'))[grepl('IAMhiFT3b-qPCR', list.files(paste0(procDataDir, '/hiFT/')))]
conditionMap <- as_tibble(read.table(paste0(procDataDir, '/hiFT/IAMhiFT3b_conditionMap - Sheet1.tsv'), sep = '\t', header = T, stringsAsFactors = F))
targetFileMaster <- as_tibble(read.table(paste0(procDataDir, '/hiFT/IAMhiFT3b qPCR filelist master - Sheet1.tsv'), sep = '\t', header = T, stringsAsFactors = F))
exp <- 'IAMhiFT3b'


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
for (i in 1:4) {
  
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
qPCR_avg_results$target <- factor(qPCR_avg_results$target, levels = levels(hiFT3b_colonies_sum$target)[levels(hiFT3b_colonies_sum$target) != 'control'])

qPCR_rep_results$target <- factor(qPCR_rep_results$target)
qPCR_rep_results$target <- factor(qPCR_rep_results$target, levels = levels(hiFT3b_colonies_sum$target)[levels(hiFT3b_colonies_sum$target) != 'control'])


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
# ggsave(KDbarplot1, file = '~/Dropbox (RajLab)/Projects/cellid/graphs/hiFT/IAMhiFT3b/KDbarplot_forSlide.pdf', width = 10, height = 5)  
# 

## merge the two

# colony_percKD_dotsNbars_forSlide <- ggplot() + 
#   geom_point(data = hiFT3b_colonies %>% filter(target != 'PO') %>% mutate(measurement = 'colonies'), aes(shIndex, colonies), position = position_jitter(w = 0.15, h = 0), alpha = 0.2) +
#   geom_bar(data = hiFT3b_colonies_sum %>% filter(target != 'PO') %>% mutate(measurement = 'colonies'), aes(shIndex, meanCount), alpha = 0.4, width = 0.7, stat = 'identity') +
#   geom_errorbar(data = hiFT3b_colonies_sum %>% filter(target != 'PO') %>% mutate(measurement = 'colonies'), aes(x = shIndex, ymin = -errbMin, ymax = errbMax), width = 0.25, alpha = 0.5) +
#   facet_grid(measurement~target, scales = 'free_y') +
#   geom_bar(data = ddCT %>% dplyr::rename(shIndex = shID) %>% mutate(measurement = 'qPCR'), aes(shIndex, percKD), stat = 'identity') +
#   geom_hline(yintercept = controlMeans$meanCtlCount, color = 'blue', alpha = 0.8) +
#   theme_classic() +
#   ylab('Alkaline Phosphatase-positive colonies per well') +
#   xlab('shRNA ID for target') +
#   ggtitle('hiF-T iPSC colony counts after perturbable factor knockdown') +
#   theme(axis.text.y = element_text(size = rel(2)),
#         axis.title = element_blank(),
#         plot.title = element_blank(),
#         strip.text = element_text(size = rel(1.6)))
# ggsave(colony_percKD_dotsNbars_forSlide, file = paste0(projectDir, 'graphs/hiFT/IAMhiFT3b/colony_percKD_counts_dotsNbars_forSlide.pdf'), width = 10, height = 5)  

hits <- c('control', 'CEBPB', 'RUNX1', 'PRRX2')
colony_percKD_dotsNbars_both_bothErr_hitsOnly <- ggplot() + 
  geom_point(data = hiFT3b_colonies_127 %>% filter(target != 'PO') %>% mutate(measurement = 'colonies') %>% filter(target %in% hits), aes(shIndex, colonies), position = position_jitter(w = 0.15, h = 0), alpha = 0.2) +
  geom_bar(data = hiFT3b_colonies_sum_127 %>% filter(target != 'PO') %>% mutate(measurement = 'colonies') %>% filter(target %in% hits), aes(shIndex, meanCount), alpha = 0.4, width = 0.7, stat = 'identity') +
  geom_errorbar(data = hiFT3b_colonies_sum_127 %>% filter(target != 'PO') %>% mutate(measurement = 'colonies') %>% filter(target %in% hits), aes(x = shIndex, ymin = errbMin, ymax = errbMax), width = 0.25, alpha = 0.5) +
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
ggsave(colony_percKD_dotsNbars_both_bothErr_hitsOnly, file = paste0(graphDir, '/hiFT/IAMhiFT3b/colony_percKD_dotsNbars_both_bothErr_hitsOnly.pdf'), width = 5, height = 5)  

#fisher's exact

targets = unique(as.character(hiFT3b_colonies$target))[!unique(as.character(hiFT3b_colonies$target)) %in% c('PO', 'control')]

ctlDat <- hiFT3b_colonies %>%
  filter(shID %in% c('SHC001', 'SHC002', 'SHC007'))

FET_pvals_3b <- list()
for (targ in targets){
  
  tempDat <- hiFT3b_colonies %>%
    filter(target == targ)
  
  tempFET <- fisher.test(matrix(c(sum(tempDat$colonies), length(tempDat$colonies)*1e4 - sum(tempDat$colonies),  
                                  sum(ctlDat$colonies), length(ctlDat$colonies)*1e4 - sum(ctlDat$colonies)), ncol = 2), alternative = 'g')
  
  if(is.null(dim(FET_pvals_3b))){
    FET_pvals_3b = tibble(target = targ,
                          FET_pval = tempFET$p.value)
  } else {
    FET_pvals_3b %<>% bind_rows(tibble(target = targ, FET_pval = tempFET$p.value))
  }
  
}

write.csv(FET_pvals_3b, file = paste0(graphDir, '/hiFT/IAMhiFT3b/FET_pvals_3b.csv'))


