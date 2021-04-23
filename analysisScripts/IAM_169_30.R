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

experiments = c('20190701_IAM11169-30-subset/')

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
scanSpotCounts$condition <- factor(scanSpotCounts$condition, levels = levels(scanSpotCounts$condition)[c(8,2,5,4,1,7,6,3)])

set.seed(7435)
scanScatters1245 <- ggplot(scanSpotCounts %>% filter(wellName %in% c('w1', 'w2', 'w4', 'w5')), aes(alexa, cy)) +
  geom_point() +
  facet_grid(experiment~condition) +
  theme_bw()
ggsave(scanScatters1245, file = paste0(graphDir, '/transdiff/HCF_169-30/GM11169-30_scan_counts_empty_7F_14_7UP_scatterFull.pdf'), width = 10, height = 10, useDingbats = F)  

scanScatters1245z <- ggplot(scanSpotCounts %>% filter(wellName %in% c('w1', 'w2', 'w4', 'w5')), aes(alexa, cy)) +
  geom_point() +
  facet_grid(experiment~condition) +
  theme_bw() +
  xlim(c(0,150))
ggsave(scanScatters1245z, file = paste0(graphDir, '/transdiff/HCF_169-30/GM11169-30_scan_counts_empty_7F_14_7UP_scatterZoom.pdf'), width = 10, height = 10, useDingbats = F)  


