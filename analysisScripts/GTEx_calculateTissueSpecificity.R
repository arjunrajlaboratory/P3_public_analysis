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

## Calculate Tissue-specificity of mean gene expression in GTEx data

library(dplyr)
library(tidyr)
library(tibble)
library(magrittr)

# if running manually in Rstudio, start here and use:
# projectDir = '~/Dropbox (RajLab)/Shared_IanM/cellid_201807_onward/'
# procDataSubdir = 'procDataScripted'
# graphSubdir = 'graphs'

procDataDir = paste0(projectDir, procDataSubdir)
graphDir = paste0(projectDir, graphSubdir)

# setwd(projectDir)
# cat("Starting fibsAndiCards_all6_DEandPert.R...\n")
# source(paste0(projectDir, "analysisScripts/RNAtagUtilities.R"))
# cat(paste0("Working in ", getwd(), "\n"))
# 
# tfGeneIDfile = "annotations/TF_gene_ids.csv"
# tfTab = as_tibble(read.csv(tfGeneIDfile))
# colnames(tfTab) = "gene_id"
# 
# geneIDfile = "annotations/hg19gene_idToGeneSymbol.tsv"
# geneIDtoGeneSymbol <- as_tibble(read.csv(geneIDfile, sep='\t', stringsAsFactors = F))
# 
# markerFile = "miscellaneous/markerGeneList.csv"
# markerTab = as_tibble(read.csv(markerFile, header = T, stringsAsFactors = F)) %>%
#   left_join(geneIDtoGeneSymbol)
# 
# zhouOlson_activators = as_tibble(read.table("miscellaneous/ZhouOlson2017_TableS1_activators_reorg.txt", sep = "\t", header = T, stringsAsFactors = F))
# zhouOlson_inhibitors = as_tibble(read.table("miscellaneous/ZhouOlson2017_TableS1_inhibitors_reorg.txt", sep = "\t", header = T, stringsAsFactors = F))
# 

gtex_smts_tpms_meanAndCV = as_tibble(read.table(paste0(procDataDir, "/GTEx/V7/SMTS_summaryStats.txt"), stringsAsFactors = F, header = T, sep = "\t"))
# gtex_smtsd_tpms_meanAndCV = as_tibble(read.table(paste0("procDataScripted/GTEx/V7/SMTSD_summaryStats.txt"), stringsAsFactors = F, header = T, sep = "\t"))
indiv_perturbability = as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_variability_metrics.txt"), header = T, stringsAsFactors = F, sep = "\t"))


## tissue-specificity score of http://genesdev.cshlp.org/content/25/18/1915.full#sec-16
## briefly, compare each gene mean vector against vector e^t where e_i^t = {1, i == t; 0, otherwise}
## JS_sp (gene|tissue) = 1 - JS_dist(gene mean vector, e^tissue)
## JS_dist (g, e) = sqrt(JS(g, e))
## JS(g, e) = H((g + e)/2) - (H(g) + H(e))/2
## H(g) = entropy = -sum(plog1/p); by convention -0log0 = 0
## p = log2(mean + 1)^tissue / sum(log2(mean + 1)^tissues)

## For SMTS sample breakdown
H <- function(v) { # from https://stackoverflow.com/questions/11226627/jensen-shannon-divergence-in-r
  v <- v[v > 0]
  return(sum(-v * log(v)))
}

gtex_smts_meanTPMprime_perTissue = gtex_smts_tpms_meanAndCV %>%
  dplyr::select(gene_id, SMTS, meanTPM) %>%
  group_by(gene_id) %>%
  mutate(geneDenom = sum(log2(meanTPM + 1)),
         meanTPMprime = log2(meanTPM + 1)/geneDenom,
         check = sum(meanTPMprime)) %>%
  dplyr::select(gene_id, SMTS, meanTPMprime) %>%
  spread(SMTS, meanTPMprime)
gtex_smts_meanTPMprime_perTissue = data.frame(gtex_smts_meanTPMprime_perTissue)
rownames(gtex_smts_meanTPMprime_perTissue) = gtex_smts_meanTPMprime_perTissue$gene_id
gtex_smts_meanTPMprime_perTissue$gene_id = NULL

testgtex = gtex_smts_meanTPMprime_perTissue[1:4,]
rAdip = matrix(0, dim(testgtex)[1], dim(testgtex)[2])
rAdip[,1] = 1 

nGene = dim(gtex_smts_meanTPMprime_perTissue)[1]
nTis = dim(gtex_smts_meanTPMprime_perTissue)[2]
gtex_jssp_smts = list()
for (tisCol in 1:nTis){
  
  et = matrix(0, nGene, nTis)
  et[,tisCol] = 1 
  jssp = 1 - sqrt(apply((gtex_smts_meanTPMprime_perTissue + et)/2, 1, function(x) H(x)) - apply(gtex_smts_meanTPMprime_perTissue, 1, function(x) H(x))/2) # note exclusion of H(et) in right-hand side because for unit vector H(et) = 0
  
  jsspTemp = as_tibble(data.frame(
    gene_id = names(jssp),
    value = as.numeric(jssp)
  ))
  colnames(jsspTemp)[2] = colnames(gtex_smts_meanTPMprime_perTissue)[tisCol]
  
  if(is.null(dim(gtex_jssp_smts))) {
    gtex_jssp_smts = jsspTemp 
  } else {
    gtex_jssp_smts = left_join(gtex_jssp_smts, jsspTemp, by = "gene_id")
  }
  
  cat(paste0("Done with ", colnames(gtex_smts_meanTPMprime_perTissue)[tisCol], "\n"))
  
}

gtex_jssp_smts %<>%
  gather(SMTS, Specificity, -gene_id) %>%
  group_by(gene_id) %>%
  mutate(Max.Specificity.GTEx = max(Specificity)) %>%
  group_by(gene_id, Max.Specificity.GTEx) %>%
  spread(SMTS, Specificity)

write.table(gtex_jssp_smts, file = paste0(procDataDir, "/GTEx/SMTS_JSsp_TissueSpecificity.txt"), sep = "\t", quote = F, row.names = F)


# ## For SMTSD sample breakdown
# gtex_smtsd_meanTPMprime_perTissue = gtex_smtsd_tpms_meanAndCV %>%
#   dplyr::select(gene_id, SMTSD, meanTPM) %>%
#   group_by(gene_id) %>%
#   mutate(geneDenom = sum(log2(meanTPM + 1)),
#          meanTPMprime = log2(meanTPM + 1)/geneDenom,
#          check = sum(meanTPMprime)) %>%
#   dplyr::select(gene_id, SMTSD, meanTPMprime) %>%
#   spread(SMTSD, meanTPMprime)
# gtex_smtsd_meanTPMprime_perTissue = data.frame(gtex_smtsd_meanTPMprime_perTissue)
# rownames(gtex_smtsd_meanTPMprime_perTissue) = gtex_smtsd_meanTPMprime_perTissue$gene_id
# gtex_smtsd_meanTPMprime_perTissue$gene_id = NULL
# 
# nGene = dim(gtex_smtsd_meanTPMprime_perTissue)[1]
# nTis = dim(gtex_smtsd_meanTPMprime_perTissue)[2]
# gtex_jssp_smtsd = list()
# for (tisCol in 1:nTis){
#   
#   et = matrix(0, nGene, nTis)
#   et[,tisCol] = 1 
#   jssp = 1 - sqrt(apply((gtex_smtsd_meanTPMprime_perTissue + et)/2, 1, function(x) H(x)) - apply(gtex_smtsd_meanTPMprime_perTissue, 1, function(x) H(x))/2)
#   
#   jsspTemp = as_tibble(data.frame(
#     gene_id = names(jssp),
#     value = as.numeric(jssp)
#   ))
#   colnames(jsspTemp)[2] = colnames(gtex_smtsd_meanTPMprime_perTissue)[tisCol]
#   
#   if(is.null(dim(gtex_jssp_smtsd))) {
#     gtex_jssp_smtsd = jsspTemp 
#   } else {
#     gtex_jssp_smtsd = left_join(gtex_jssp_smtsd, jsspTemp, by = "gene_id")
#   }
#   
#   cat(paste0("Done with ", colnames(gtex_smtsd_meanTPMprime_perTissue)[tisCol], "\n"))
#   
# }
# 
# gtex_jssp_smtsd %<>%
#   gather(SMTSD, Specificity, -gene_id) %>%
#   group_by(gene_id) %>%
#   mutate(Max.Specificity.GTEx.SMTSD = max(Specificity)) %>%
#   group_by(gene_id, Max.Specificity.GTEx.SMTSD) %>%
#   spread(SMTSD, Specificity)
# 
# write.table(gtex_jssp_smtsd, file = paste0(projectDir, "procDataScripted/GTEx/SMTSD_JSsp_TissueSpecificity.txt"), sep = "\t", quote = F, row.names = F)


#### Add my fibroblast data and take out Skin
# SMTS

# gtex_smts_meanTPMprime_perTissue_942 <- gtex_smts_meanTPMprime_perTissue %>%
#   inner_join(indiv_perturbability %>% filter(cellType == 'GM00942') %>% dplyr::select(gene_id, meanTPM), by = 'gene_id') %>%
#   dplyr::select(-Skin)

gtex_smts_meanTPMprime_perTissue_942_noSkin = gtex_smts_tpms_meanAndCV %>%
  filter(SMTS != 'Skin',
         gene_id %in% unique(indiv_perturbability$gene_id)) %>%
  dplyr::select(gene_id, SMTS, meanTPM) %>%
  bind_rows(indiv_perturbability %>% 
              filter(cellType == 'GM00942', 
                     gene_id %in% unique(gtex_smts_tpms_meanAndCV$gene_id)) %>% 
              dplyr::select(cellType, gene_id, meanTPM) %>%
              dplyr::rename(SMTS = cellType)) %>%
  group_by(gene_id) %>%
  mutate(geneDenom = sum(log2(meanTPM + 1)),
         meanTPMprime = log2(meanTPM + 1)/geneDenom,
         check = sum(meanTPMprime)) %>%
  dplyr::select(gene_id, SMTS, meanTPMprime) %>% 
  ungroup() %>% 
  group_by(gene_id) %>%
  spread(SMTS, meanTPMprime)

gtex_smts_meanTPMprime_perTissue_942_noSkin = data.frame(gtex_smts_meanTPMprime_perTissue_942_noSkin)
rownames(gtex_smts_meanTPMprime_perTissue_942_noSkin) = gtex_smts_meanTPMprime_perTissue_942_noSkin$gene_id
gtex_smts_meanTPMprime_perTissue_942_noSkin$gene_id = NULL

gtex_smts_meanTPMprime_perTissue_942_noSkin = gtex_smts_meanTPMprime_perTissue_942_noSkin[rowMeans(is.na(gtex_smts_meanTPMprime_perTissue_942_noSkin)) < 1,] # filter out rows with no expression

testgtex = gtex_smts_meanTPMprime_perTissue_942_noSkin[1:4,]
rAdip = matrix(0, dim(testgtex)[1], dim(testgtex)[2])
rAdip[,1] = 1 

nGene = dim(gtex_smts_meanTPMprime_perTissue_942_noSkin)[1]
nTis = dim(gtex_smts_meanTPMprime_perTissue_942_noSkin)[2]
gtex_jssp_smts_942_noSkin = list()
for (tisCol in 1:nTis){
  
  et = matrix(0, nGene, nTis)
  et[,tisCol] = 1 
  jssp = 1 - sqrt(apply((gtex_smts_meanTPMprime_perTissue_942_noSkin + et)/2, 1, function(x) H(x)) - apply(gtex_smts_meanTPMprime_perTissue_942_noSkin, 1, function(x) H(x))/2)
  
  jsspTemp = as_tibble(data.frame(
    gene_id = names(jssp),
    value = as.numeric(jssp)
  ))
  colnames(jsspTemp)[2] = colnames(gtex_smts_meanTPMprime_perTissue_942_noSkin)[tisCol]
  
  if(is.null(dim(gtex_jssp_smts_942_noSkin))) {
    gtex_jssp_smts_942_noSkin = jsspTemp 
  } else {
    gtex_jssp_smts_942_noSkin = left_join(gtex_jssp_smts_942_noSkin, jsspTemp, by = "gene_id")
  }
  
  cat(paste0("Done with ", colnames(gtex_smts_meanTPMprime_perTissue_942_noSkin)[tisCol], "\n"))
  
}


gtex_jssp_smts_942_noSkin %<>%
  gather(SMTS, Specificity, -gene_id) %>%
  group_by(gene_id) %>%
  mutate(Max.Specificity.GTEx = max(Specificity)) %>%
  group_by(gene_id, Max.Specificity.GTEx) %>%
  spread(SMTS, Specificity)

write.table(gtex_jssp_smts_942_noSkin, file = paste0(procDataDir, "/GTEx/SMTS_JSsp_TissueSpecificity_942_noSkin.txt"), sep = "\t", quote = F, row.names = F)


## Add my iCards and take away heart

gtex_smts_meanTPMprime_perTissue_iCard_noHeart = gtex_smts_tpms_meanAndCV %>%
  filter(SMTS != 'Heart',
         gene_id %in% unique(indiv_perturbability$gene_id)) %>%
  dplyr::select(gene_id, SMTS, meanTPM) %>%
  bind_rows(indiv_perturbability %>% 
              filter(cellType == 'iCard', 
                     gene_id %in% unique(gtex_smts_tpms_meanAndCV$gene_id)) %>% 
              dplyr::select(cellType, gene_id, meanTPM) %>%
              dplyr::rename(SMTS = cellType)) %>%
  group_by(gene_id) %>%
  mutate(geneDenom = sum(log2(meanTPM + 1)),
         meanTPMprime = log2(meanTPM + 1)/geneDenom,
         check = sum(meanTPMprime)) %>%
  dplyr::select(gene_id, SMTS, meanTPMprime) %>% 
  ungroup() %>% 
  group_by(gene_id) %>%
  spread(SMTS, meanTPMprime)

gtex_smts_meanTPMprime_perTissue_iCard_noHeart = data.frame(gtex_smts_meanTPMprime_perTissue_iCard_noHeart)
rownames(gtex_smts_meanTPMprime_perTissue_iCard_noHeart) = gtex_smts_meanTPMprime_perTissue_iCard_noHeart$gene_id
gtex_smts_meanTPMprime_perTissue_iCard_noHeart$gene_id = NULL

gtex_smts_meanTPMprime_perTissue_iCard_noHeart = gtex_smts_meanTPMprime_perTissue_iCard_noHeart[rowMeans(is.na(gtex_smts_meanTPMprime_perTissue_iCard_noHeart)) < 1,] # filter out rows with no expression

testgtex = gtex_smts_meanTPMprime_perTissue_iCard_noHeart[1:4,]
rAdip = matrix(0, dim(testgtex)[1], dim(testgtex)[2])
rAdip[,1] = 1 

nGene = dim(gtex_smts_meanTPMprime_perTissue_iCard_noHeart)[1]
nTis = dim(gtex_smts_meanTPMprime_perTissue_iCard_noHeart)[2]
gtex_jssp_smts_iCard_noHeart = list()
for (tisCol in 1:nTis){
  
  et = matrix(0, nGene, nTis)
  et[,tisCol] = 1 
  jssp = 1 - sqrt(apply((gtex_smts_meanTPMprime_perTissue_iCard_noHeart + et)/2, 1, function(x) H(x)) - apply(gtex_smts_meanTPMprime_perTissue_iCard_noHeart, 1, function(x) H(x))/2)
  
  jsspTemp = as_tibble(data.frame(
    gene_id = names(jssp),
    value = as.numeric(jssp)
  ))
  colnames(jsspTemp)[2] = colnames(gtex_smts_meanTPMprime_perTissue_iCard_noHeart)[tisCol]
  
  if(is.null(dim(gtex_jssp_smts_iCard_noHeart))) {
    gtex_jssp_smts_iCard_noHeart = jsspTemp 
  } else {
    gtex_jssp_smts_iCard_noHeart = left_join(gtex_jssp_smts_iCard_noHeart, jsspTemp, by = "gene_id")
  }
  
  cat(paste0("Done with ", colnames(gtex_smts_meanTPMprime_perTissue_iCard_noHeart)[tisCol], "\n"))
  
}


gtex_jssp_smts_iCard_noHeart %<>%
  gather(SMTS, Specificity, -gene_id) %>%
  group_by(gene_id) %>%
  mutate(Max.Specificity.GTEx = max(Specificity)) %>%
  group_by(gene_id, Max.Specificity.GTEx) %>%
  spread(SMTS, Specificity)

write.table(gtex_jssp_smts_iCard_noHeart, file = paste0(procDataDir, "/GTEx/SMTS_JSsp_TissueSpecificity_iCard_noHeart.txt"), sep = "\t", quote = F, row.names = F)


## Add my fibs and iCards and don't take away skin or heart

gtex_smts_meanTPMprime_perTissue_iCard_fibs = gtex_smts_tpms_meanAndCV %>%
  filter(gene_id %in% unique(indiv_perturbability$gene_id)) %>%
  dplyr::select(gene_id, SMTS, meanTPM) %>%
  bind_rows(indiv_perturbability %>% 
              filter(cellType == 'iCard' | cellType == 'GM00942', 
                     gene_id %in% unique(gtex_smts_tpms_meanAndCV$gene_id)) %>% 
              dplyr::select(cellType, gene_id, meanTPM) %>%
              dplyr::rename(SMTS = cellType)) %>%
  group_by(gene_id) %>%
  mutate(geneDenom = sum(log2(meanTPM + 1)),
         meanTPMprime = log2(meanTPM + 1)/geneDenom,
         check = sum(meanTPMprime)) %>%
  dplyr::select(gene_id, SMTS, meanTPMprime) %>% 
  ungroup() %>% 
  group_by(gene_id) %>%
  spread(SMTS, meanTPMprime)

gtex_smts_meanTPMprime_perTissue_iCard_fibs = data.frame(gtex_smts_meanTPMprime_perTissue_iCard_fibs)
rownames(gtex_smts_meanTPMprime_perTissue_iCard_fibs) = gtex_smts_meanTPMprime_perTissue_iCard_fibs$gene_id
gtex_smts_meanTPMprime_perTissue_iCard_fibs$gene_id = NULL

gtex_smts_meanTPMprime_perTissue_iCard_fibs = gtex_smts_meanTPMprime_perTissue_iCard_fibs[rowMeans(is.na(gtex_smts_meanTPMprime_perTissue_iCard_fibs)) < 1,] # filter out rows with no expression

testgtex = gtex_smts_meanTPMprime_perTissue_iCard_fibs[1:4,]
rAdip = matrix(0, dim(testgtex)[1], dim(testgtex)[2])
rAdip[,1] = 1 

nGene = dim(gtex_smts_meanTPMprime_perTissue_iCard_fibs)[1]
nTis = dim(gtex_smts_meanTPMprime_perTissue_iCard_fibs)[2]
gtex_jssp_smts_iCard_fibs = list()
for (tisCol in 1:nTis){
  
  et = matrix(0, nGene, nTis)
  et[,tisCol] = 1 
  jssp = 1 - sqrt(apply((gtex_smts_meanTPMprime_perTissue_iCard_fibs + et)/2, 1, function(x) H(x)) - apply(gtex_smts_meanTPMprime_perTissue_iCard_fibs, 1, function(x) H(x))/2)
  
  jsspTemp = as_tibble(data.frame(
    gene_id = names(jssp),
    value = as.numeric(jssp)
  ))
  colnames(jsspTemp)[2] = colnames(gtex_smts_meanTPMprime_perTissue_iCard_fibs)[tisCol]
  
  if(is.null(dim(gtex_jssp_smts_iCard_fibs))) {
    gtex_jssp_smts_iCard_fibs = jsspTemp 
  } else {
    gtex_jssp_smts_iCard_fibs = left_join(gtex_jssp_smts_iCard_fibs, jsspTemp, by = "gene_id")
  }
  
  cat(paste0("Done with ", colnames(gtex_smts_meanTPMprime_perTissue_iCard_fibs)[tisCol], "\n"))
  
}


gtex_jssp_smts_iCard_fibs %<>%
  gather(SMTS, Specificity, -gene_id) %>%
  group_by(gene_id) %>%
  mutate(Max.Specificity.GTEx = max(Specificity)) %>%
  group_by(gene_id, Max.Specificity.GTEx) %>%
  spread(SMTS, Specificity)

write.table(gtex_jssp_smts_iCard_fibs, file = paste0(procDataDir, "/GTEx/SMTS_JSsp_TissueSpecificity_iCard_fibs.txt"), sep = "\t", quote = F, row.names = F)


## Add my fibs and iCards and do take away skin and heart

gtex_smts_meanTPMprime_perTissue_iCard_fibs_noHeart_noSkin = gtex_smts_tpms_meanAndCV %>%
  filter(gene_id %in% unique(indiv_perturbability$gene_id),
         SMTS != 'Heart', SMTS != 'Skin') %>%
  dplyr::select(gene_id, SMTS, meanTPM) %>%
  bind_rows(indiv_perturbability %>% 
              filter(cellType == 'iCard' | cellType == 'GM00942', 
                     gene_id %in% unique(gtex_smts_tpms_meanAndCV$gene_id)) %>% 
              dplyr::select(cellType, gene_id, meanTPM) %>%
              dplyr::rename(SMTS = cellType)) %>%
  group_by(gene_id) %>%
  mutate(geneDenom = sum(log2(meanTPM + 1)),
         meanTPMprime = log2(meanTPM + 1)/geneDenom,
         check = sum(meanTPMprime)) %>%
  dplyr::select(gene_id, SMTS, meanTPMprime) %>% 
  ungroup() %>% 
  group_by(gene_id) %>%
  spread(SMTS, meanTPMprime)

gtex_smts_meanTPMprime_perTissue_iCard_fibs_noHeart_noSkin = data.frame(gtex_smts_meanTPMprime_perTissue_iCard_fibs_noHeart_noSkin)
rownames(gtex_smts_meanTPMprime_perTissue_iCard_fibs_noHeart_noSkin) = gtex_smts_meanTPMprime_perTissue_iCard_fibs_noHeart_noSkin$gene_id
gtex_smts_meanTPMprime_perTissue_iCard_fibs_noHeart_noSkin$gene_id = NULL

gtex_smts_meanTPMprime_perTissue_iCard_fibs_noHeart_noSkin = gtex_smts_meanTPMprime_perTissue_iCard_fibs_noHeart_noSkin[rowMeans(is.na(gtex_smts_meanTPMprime_perTissue_iCard_fibs_noHeart_noSkin)) < 1,] # filter out rows with no expression

testgtex = gtex_smts_meanTPMprime_perTissue_iCard_fibs_noHeart_noSkin[1:4,]
rAdip = matrix(0, dim(testgtex)[1], dim(testgtex)[2])
rAdip[,1] = 1 

nGene = dim(gtex_smts_meanTPMprime_perTissue_iCard_fibs_noHeart_noSkin)[1]
nTis = dim(gtex_smts_meanTPMprime_perTissue_iCard_fibs_noHeart_noSkin)[2]
gtex_jssp_smts_iCard_fibs_noHeart_noSkin = list()
for (tisCol in 1:nTis){
  
  et = matrix(0, nGene, nTis)
  et[,tisCol] = 1 
  jssp = 1 - sqrt(apply((gtex_smts_meanTPMprime_perTissue_iCard_fibs_noHeart_noSkin + et)/2, 1, function(x) H(x)) - apply(gtex_smts_meanTPMprime_perTissue_iCard_fibs_noHeart_noSkin, 1, function(x) H(x))/2)
  
  jsspTemp = as_tibble(data.frame(
    gene_id = names(jssp),
    value = as.numeric(jssp)
  ))
  colnames(jsspTemp)[2] = colnames(gtex_smts_meanTPMprime_perTissue_iCard_fibs_noHeart_noSkin)[tisCol]
  
  if(is.null(dim(gtex_jssp_smts_iCard_fibs_noHeart_noSkin))) {
    gtex_jssp_smts_iCard_fibs_noHeart_noSkin = jsspTemp 
  } else {
    gtex_jssp_smts_iCard_fibs_noHeart_noSkin = left_join(gtex_jssp_smts_iCard_fibs_noHeart_noSkin, jsspTemp, by = "gene_id")
  }
  
  cat(paste0("Done with ", colnames(gtex_smts_meanTPMprime_perTissue_iCard_fibs_noHeart_noSkin)[tisCol], "\n"))
  
}


gtex_jssp_smts_iCard_fibs_noHeart_noSkin %<>%
  gather(SMTS, Specificity, -gene_id) %>%
  group_by(gene_id) %>%
  mutate(Max.Specificity.GTEx = max(Specificity)) %>%
  group_by(gene_id, Max.Specificity.GTEx) %>%
  spread(SMTS, Specificity)

write.table(gtex_jssp_smts_iCard_fibs_noHeart_noSkin, file = paste0(procDataDir, "/GTEx/SMTS_JSsp_TissueSpecificity_iCard_fibs_noHeart_noSkin.txt"), sep = "\t", quote = F, row.names = F)


## Diff exp-based method for tissue specificity.
#  As in this paper, , calculate lfc of each gene in each tissue against average of all other tissues in the dataset.

# tistypes = unique(gtex_smts_tpms_meanAndCV$SMTS)
# pseudo = 1
# gtex_smts_lfc <- list()
# for (ref_tis in tistypes) {
#   
#   DEfc_spec_temp <- gtex_smts_tpms_meanAndCV %>%
#     dplyr::select(gene_id, SMTS, meanTPM) %>%
#     mutate(ref_tissue = ref_tis,
#            isRef = ref_tissue == SMTS) %>%
#     group_by(gene_id, isRef, ref_tissue) %>%
#     summarise(meanMeanTPM = mean(meanTPM),
#               meanMeanTPM_plusPseudo = meanMeanTPM + pseudo) %>%
#     dplyr::select(-meanMeanTPM) %>%
#     ungroup() %>%
#     spread(isRef, meanMeanTPM_plusPseudo) %>%
#     mutate(lfc_ref = log2(`TRUE`/`FALSE`),
#            ref_meanTPM = `TRUE`) %>%
#     dplyr::select(gene_id, ref_tissue, ref_meanTPM, lfc_ref)
#   
#   if(is.null(dim(gtex_smts_lfc))) {
#     gtex_smts_lfc <- DEfc_spec_temp
#   } else {
#     gtex_smts_lfc %<>% bind_rows(DEfc_spec_temp)
#   }
#   
#   cat(paste0("Done with ", ref_tis, "\n"))
#   
# }
# 
# 
# # with iCards and fibroblasts (including skin and heart)
# 
# gtex_smts_tpms_means_fibs_iCards <- gtex_smts_tpms_meanAndCV %>%
#   dplyr::select(gene_id, SMTS, meanTPM) %>%
#   bind_rows(indiv_perturbability %>%
#               dplyr::select(gene_id, cellType, meanTPM) %>%
#               dplyr::rename(SMTS = cellType))
# 
# tistypes = unique(gtex_smts_tpms_means_fibs_iCards$SMTS)
# pseudo = 1
# gtex_smts_lfc_fibs_iCards <- list()
# for (ref_tis in tistypes) {
#   
#   DEfc_spec_temp <- gtex_smts_tpms_means_fibs_iCards %>%
#     dplyr::select(gene_id, SMTS, meanTPM) %>%
#     mutate(ref_tissue = ref_tis,
#            isRef = ref_tissue == SMTS) %>%
#     group_by(gene_id, isRef, ref_tissue) %>%
#     summarise(meanMeanTPM = mean(meanTPM),
#               meanMeanTPM_plusPseudo = meanMeanTPM + pseudo) %>%
#     dplyr::select(-meanMeanTPM) %>%
#     ungroup() %>%
#     spread(isRef, meanMeanTPM_plusPseudo) %>%
#     mutate(lfc_ref = log2(`TRUE`/`FALSE`),
#            ref_meanTPM = `TRUE`) %>%
#     dplyr::select(gene_id, ref_tissue, ref_meanTPM, lfc_ref)
#   
#   if(is.null(dim(gtex_smts_lfc_fibs_iCards))) {
#     gtex_smts_lfc_fibs_iCards <- DEfc_spec_temp
#   } else {
#     gtex_smts_lfc_fibs_iCards %<>% bind_rows(DEfc_spec_temp)
#   }
#   
#   cat(paste0("Done with ", ref_tis, "\n"))
#   
# }




