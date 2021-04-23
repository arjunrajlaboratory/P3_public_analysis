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

# projectDir = '~/Dropbox (RajLab)/Shared_IanM/cellid_201807_onward/'
# procDataSubdir = 'processedImages'
# graphSubdir = 'graphs'

setwd(projectDir)
library(tidyverse)
library(magrittr)
library(yaml)

procDataDir = paste0(projectDir, procDataSubdir)
graphDir = paste0(projectDir, graphSubdir)

procDataSubdir2 <- paste0(procDataSubdir, '/transdiff/')

experiments = c('20190629_HCF30-subset/')

# stackSpotCounts <- list()
scanSpotCounts <- list() 
for (exp in experiments) {
  
  readme <- read_yaml(paste0(projectDir, procDataSubdir2, exp, '/readMe.yaml'))
  
  wellMap <- tibble(
    wellName = strsplit(readme[["wellMap"]][['wellName']], ',')[[1]],
    condition = strsplit(readme[["wellMap"]][['condition']], ',')[[1]]
  )
  
  expt = readme$ID
  cells = readme$sample$cellLine
  
  for (wind in 1:length(wellMap$wellName)) {
    wellN <- wellMap$wellName[wind]
    cond <- wellMap$condition[wind]
    
    # stackSpotCounts_temp <- as_tibble(read.csv(paste0(projectDir, procDataSubdir2, exp, wellN, '_stacks/counts.csv'), header = T, stringsAsFactors = F)) %>%
    #   dplyr::rename(objectNumber = objNum, arrayNumber = objArrayNum) %>% dplyr::select(-isGood) %>%
    #   inner_join(as_tibble(read.csv(paste0(projectDir, procDataSubdir2, exp, wellN, '_stacks/cellMaskSize.csv'), header = T, stringsAsFactors = F)), by = c('arrayNumber', 'objectNumber')) %>%
    #   mutate(condition = cond,
    #          wellName = wellN,
    #          cellLine = cells,
    #          experiment = expt) %>% filter(isGood == 1)
    
    scanSpotCounts_temp <- as_tibble(read.csv(paste0(projectDir, procDataSubdir2, exp, wellN, '_scan/centroidsAndTheirSpots.csv'), header = T, stringsAsFactors = F)) %>%
      mutate(wellName = wellN) %>%
      inner_join(wellMap)
    
    if(is.null(dim(scanSpotCounts))) {
      # stackSpotCounts <- stackSpotCounts_temp
      scanSpotCounts <- scanSpotCounts_temp
    } else {
      # stackSpotCounts %<>% bind_rows(stackSpotCounts_temp)
      scanSpotCounts %<>% bind_rows(scanSpotCounts_temp)
    }
    
  }
}
# stackSpotCounts$condition <- factor(stackSpotCounts$condition)
# stackSpotCounts$condition <- factor(stackSpotCounts$condition, levels = levels(stackSpotCounts$condition)[c(4,2,1,3)])

scanSpotCounts$condition <- factor(scanSpotCounts$condition)
scanSpotCounts$condition <- factor(scanSpotCounts$condition, levels = levels(scanSpotCounts$condition)[c(8,2,3,5,4,1,7,6)])

fracsMarkerPos <- scanSpotCounts %>%
  filter(wellName %in% c('w1', 'w3', 'w4', 'w5')) %>%
  group_by(wellName, condition) %>%
  summarise(fracNPPAonlyPos = sum(alexa >= 15 & cy < 20)/length(alexa),
            fracTNNT2onlyPos = sum(alexa < 15 & cy >= 20)/length(alexa),
            fracDoublePos = sum(alexa >= 15 & cy >= 20)/length(alexa),
            fracNegative = sum(alexa < 15 & cy < 20)/length(alexa),
            nCells = length(alexa))

sz=0.1
set.seed(7435)
scanScatters1245 <- ggplot() +
  geom_point(data = scanSpotCounts %>% filter(wellName %in% c('w1', 'w3', 'w4', 'w5')), aes(alexa, cy),
             size = 0.1) +
  geom_text(data = fracsMarkerPos %>%
              gather('quadrantName', 'fraction',fracNPPAonlyPos:fracNegative) %>%
              mutate(alexa = ifelse(quadrantName %in% c('fracNPPAonlyPos', 'fracDoublePos'), 1000, -10),
                     cy = ifelse(quadrantName %in% c('fracDoublePos', 'fracTNNT2onlyPos'), 360, -10)), 
            aes(alexa, cy, label = as.character(round(fraction, 3)))) +
  geom_hline(data = scanSpotCounts %>% filter(wellName %in% c('w1', 'w3', 'w4', 'w5')), aes(yintercept = 20)) + 
  geom_vline(data = scanSpotCounts %>% filter(wellName %in% c('w1', 'w3', 'w4', 'w5')), aes(xintercept = 15)) +
  facet_grid(experiment~condition) +
  theme_classic() +
  xlim(c(-50,1250)) +
  ylim(c(-20,400))
ggsave(scanScatters1245, file = paste0(graphDir, '/transdiff/HCF_169-30/immHCF-30_scan_counts_empty_7F_14_7UP_scatterFull.pdf'), width = 15, height = 4, useDingbats = F)  

scanScatters1245z <- ggplot() +
  geom_point(data = scanSpotCounts %>% filter(wellName %in% c('w1', 'w3', 'w4', 'w5')), aes(alexa, cy),
             size = sz) +
  geom_text(data = fracsMarkerPos %>%
              gather('quadrantName', 'fraction',fracNPPAonlyPos:fracNegative) %>%
              mutate(alexa = ifelse(quadrantName %in% c('fracNPPAonlyPos', 'fracDoublePos'), 1000, -10),
                     cy = ifelse(quadrantName %in% c('fracDoublePos', 'fracTNNT2onlyPos'), 360, -10)), 
            aes(alexa, cy, label = as.character(round(fraction, 3)))) +
  geom_hline(data = scanSpotCounts %>% filter(wellName %in% c('w1', 'w3', 'w4', 'w5')), aes(yintercept = 20)) + 
  geom_vline(data = scanSpotCounts %>% filter(wellName %in% c('w1', 'w3', 'w4', 'w5')), aes(xintercept = 15)) +
  facet_grid(experiment~condition) +
  theme_classic() +
  xlim(c(-20,150)) +
  ylim(c(-10,400))
ggsave(scanScatters1245z, file = paste0(graphDir, '/transdiff/HCF_169-30/immHCF-30_scan_counts_empty_7F_14_7UP_scatterZoom.pdf'), width = 7.5, height = 2, useDingbats = F)  

write.csv(fracsMarkerPos, file = paste0(graphDir, '/transdiff/HCF_169-30/immHCF-30_scan_fractionsMarkerPos_empty_7F_14_7UP.csv'))
  
