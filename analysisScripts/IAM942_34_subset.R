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

procDataSubdir2 <- paste0(procDataSubdir, '/transdiff/')
graphDir = paste0(projectDir, graphSubdir)

setwd(projectDir)
library(tidyverse)
library(magrittr)
library(yaml)
### IAM942-34-subset

if(!dir.exists(paste0(graphDir, '/transdiff'))){
  dir.create(paste0(graphDir, '/transdiff'))
}

if(!dir.exists(paste0(graphDir, '/transdiff/IAM942-34-subset'))){
  dir.create(paste0(graphDir, '/transdiff/IAM942-34-subset'))
}


experimentDir <- '20190731_IAM942-34-subset/'

readme <- read_yaml(paste0(projectDir, procDataSubdir2, experimentDir, 'readMe.yaml'))

wellMap <- tibble(
  wellName = strsplit(readme[["wellMap"]][['wellName']], ',')[[1]],
  condition = strsplit(readme[["wellMap"]][['condition']], ',')[[1]]
)

stackSpotCounts <- list()
scanSpotCounts <- list()
for (wind in 1:4) {
  wellN <- wellMap$wellName[wind]
  cond <- wellMap$condition[wind]
  
  stackSpotCounts_temp <- as_tibble(read.csv(paste0(projectDir, procDataSubdir2, experimentDir, wellN, '_stacks/counts.csv'), header = T, stringsAsFactors = F)) %>%
    dplyr::rename(objectNumber = objNum, arrayNumber = objArrayNum) %>% dplyr::select(-isGood) %>%
    inner_join(as_tibble(read.csv(paste0(projectDir, procDataSubdir2, experimentDir, wellN, '_stacks/cellMaskSize.csv'), header = T, stringsAsFactors = F)), by = c('arrayNumber', 'objectNumber')) %>%
    mutate(condition = cond,
           wellName = wellN,
           cellLine = 'GM00942',
           experiment = 'IAM942-34-subset') %>% filter(isGood == 1)
  
  scanSpotCounts_temp <- as_tibble(read.csv(paste0(projectDir, procDataSubdir2, experimentDir, wellN, '_scan/centroidsAndTheirSpots.csv'), header = T, stringsAsFactors = F)) %>%
    mutate(wellName = wellN) %>%
    inner_join(wellMap)
  
  if(is.null(dim(stackSpotCounts_temp))) {
    stackSpotCounts <- stackSpotCounts_temp
    scanSpotCounts <- scanSpotCounts_temp
  } else {
    stackSpotCounts %<>% bind_rows(stackSpotCounts_temp)
    scanSpotCounts %<>% bind_rows(scanSpotCounts_temp)
  }
  
}

stackSpotCounts$condition <- factor(stackSpotCounts$condition)
stackSpotCounts$condition <- factor(stackSpotCounts$condition, levels = levels(stackSpotCounts$condition)[c(4,2,1,3)])

scanSpotCounts$condition <- factor(scanSpotCounts$condition)
scanSpotCounts$condition <- factor(scanSpotCounts$condition, levels = levels(scanSpotCounts$condition)[c(4,2,1,3)])


# scatters
set.seed(51243)
stackScatter12 = ggplot(stackSpotCounts %>% filter(condition %in% c('MXs-empty', '7F')), aes(cy.RNACounts, alexa.RNACounts)) + 
  geom_jitter(width = 10, height = 10) +
  facet_wrap(~condition) +
  theme_bw() +
  xlab('TNNT2 smFISH counts per cell') + ylab('NPPA smFISH counts per cell') +
  ggtitle('Transdifferentiation of GM00942 dermal fibroblasts')
ggsave(stackScatter12, file = paste0(graphDir, '/transdiff/IAM942-34-subset/stacks_counts_scatter_empty_7F.pdf'), width = 10, height = 5, useDingbats = F)  

stackScatter12_pixel = ggplot(stackSpotCounts %>% filter(condition %in% c('MXs-empty', '7F')), aes(cy.RNACounts/nPixels, alexa.RNACounts/nPixels)) + 
  geom_jitter() +
  facet_wrap(~condition) +
  theme_bw() +
  xlab('TNNT2 smFISH counts per pixel per cell') + ylab('NPPA smFISH counts per pixel per cell') +
  ggtitle('Transdifferentiation of GM00942 dermal fibroblasts')
ggsave(stackScatter12_pixel, file = paste0(graphDir, '/transdiff/IAM942-34-subset/stacks_countsPerArea_scatter_empty_7F.pdf'), width = 10, height = 5, useDingbats = F)  

stackScatter12_gapdh = ggplot(stackSpotCounts %>% filter(condition %in% c('MXs-empty', '7F')), aes(cy.RNACounts/tmr.RNACounts, alexa.RNACounts/tmr.RNACounts)) + 
  geom_jitter() +
  facet_wrap(~condition) +
  theme_bw() +
  xlab('TNNT2 smFISH counts normalized to GAPDH count per cell') + ylab('NPPA smFISH counts normalized to GAPDH count per cell') +
  ggtitle('Transdifferentiation of GM00942 dermal fibroblasts')
ggsave(stackScatter12_gapdh, file = paste0(graphDir, '/transdiff/IAM942-34-subset/stacks_countsPerGAPDH_scatter_empty_7F.pdf'), width = 10, height = 5, useDingbats = F)  

set.seed(96127)
scanScatterw12<- ggplot(scanSpotCounts %>% filter(wellName %in% c('well1', 'well2')), aes(alexa, cy)) +
  geom_jitter() +
  facet_wrap(~condition) +
  theme_bw()
ggsave(scanScatterw12, file = paste0(graphDir, '/transdiff/IAM942-34-subset/scan_counts_scatter_empty_7F.pdf'), width = 10, height = 5, useDingbats = F)  

fracsMarkerPos <- scanSpotCounts %>%
  filter(wellName %in% c('well1', 'well2')) %>%
  group_by(wellName, condition) %>%
  summarise(fracNPPAonlyPos = sum(alexa >= 25 & cy < 20)/length(alexa),
            fracTNNT2onlyPos = sum(alexa < 25 & cy >= 20)/length(alexa),
            fracDoublePos = sum(alexa >= 25 & cy >= 20)/length(alexa),
            fracNegative = sum(alexa < 25 & cy < 20)/length(alexa),
            nCells = length(alexa))

write.csv(fracsMarkerPos, file = paste0(graphDir, '/transdiff/IAM942-34-subset/scan_markerFracs_empty_7F.csv') )

# thumbnails
scanThumbnail12_tnnt2<- ggplot(scanSpotCounts %>% 
                           filter(wellName %in% c('well1', 'well2')), aes(xPos, yPos, color = cy)) +
  geom_point()+
  facet_wrap(~condition) +
  scale_color_distiller(palette = 'Spectral', values = c(0,0.4,0.5,0.7,1)) +
  theme(panel.background = element_rect(fill = 'black'),
          panel.grid = element_line(colour = NA))
ggsave(scanThumbnail12_tnnt2, file = paste0(graphDir, '/transdiff/IAM942-34-subset/scan_thumbnail_TNNT2_empty_7F.pdf'), width = 16, height = 8, useDingbats = F)  

scanThumbnail12_nppa<- ggplot(scanSpotCounts %>% 
                                 filter(wellName %in% c('well1', 'well2')), aes(xPos, yPos, color = alexa)) +
  geom_point()+
  facet_wrap(~condition) +
  scale_color_distiller(palette = 'Spectral', values = c(0,0.4,0.5,0.7,1)) +
  theme(panel.background = element_rect(fill = 'black'),
        panel.grid = element_line(colour = NA))
ggsave(scanThumbnail12_nppa, file = paste0(graphDir, '/transdiff/IAM942-34-subset/scan_thumbnail_NPPA_empty_7F.pdf'), width = 16, height = 8, useDingbats = F)  

scanThumbnail12_gapdh<- ggplot(scanSpotCounts %>% 
                                filter(wellName %in% c('well1', 'well2')), aes(xPos, yPos, color = tmr)) +
  geom_point()+
  facet_wrap(~condition) +
  scale_color_distiller(palette = 'Spectral', values = c(0,0.1,0.2,0.3,1)) +
  theme(panel.background = element_rect(fill = 'black'),
        panel.grid = element_line(colour = NA))
ggsave(scanThumbnail12_gapdh, file = paste0(graphDir, '/transdiff/IAM942-34-subset/scan_thumbnail_GAPDH_empty_7F.pdf'), width = 16, height = 8, useDingbats = F)  

