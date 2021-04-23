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
# procDataSubdir = 'procDataScripted_test429'
# graphSubdir = 'graphs'

procDataDir = paste0(projectDir, procDataSubdir)
graphDir = paste0(projectDir, graphSubdir)

if(!dir.exists(paste0(graphDir, '/hiFT/'))){
  dir.create(paste0(graphDir, '/hiFT/'))
}
# if(!dir.exists(paste0(graphDir, '/hiFT/hiFT_pool_2015'))){
#   dir.create(paste0(graphDir, '/hiFT/hiFT_pool_2015'))
# }

library(tidyverse)
library(magrittr)
library(ggrepel)
library(readxl)

setwd(projectDir)
tfGeneIDfile = "annotations/TF_gene_ids.csv"
tfTab = as_tibble(read.csv(tfGeneIDfile))
colnames(tfTab) = "gene_id"

geneIDfile = "annotations/hg19gene_idToGeneSymbol.tsv"
geneIDtoGeneSymbol <- as_tibble(read.csv(geneIDfile, sep='\t', stringsAsFactors = F))

## Reanalysis of Borkent et al. MEF screen
# See processedData/hiFT/Borkent/README.sh for data download information
# the dataset is a table with columns 

# hits_in_paper <- c('EZH1', 'KDM1A', 'KTI12', 'LBR', 'NAP1L3', 'RSF1', 'SHPRH')
hits_in_pertub <- c('ZNF652', 'ZBTB38', 'ATOH8', 'CERS2', 'KLF13', 'CEBPB', 'RUNX1', 'PRRX2')

pooled_sh_borkent <- as_tibble(read_xlsx(paste0(projectDir, 'extractedData/hiFT/borkent_2016/mmc2.xlsx'))) %>%
  separate(Gene_Symbol, into = 'first_Gene_Symbol', sep = ',', remove = F)
indiv_perturbability = as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_variability_metrics.txt"), header = T, stringsAsFactors = F, sep = "\t"))
indiv_perturbability_RPMnorm20 = as_tibble(read.table(paste0(procDataDir, "/allExperiments/all_rpm_readFilt_manualFilt_tpmFilt_variability_RPMnormMetrics_window20.txt"), header = T, stringsAsFactors = F, sep = "\t"))
indiv_perturbability[is.na(indiv_perturbability)] <- 0
indiv_perturbability_RPMnorm20[is.na(indiv_perturbability_RPMnorm20)] <- 0

# GTEx tissue-specificity scores
gtex_jssp_smts <- as_tibble(read.table(paste0(procDataSubdir, "/GTEx/SMTS_JSsp_TissueSpecificity.txt"), sep = "\t", header = T, stringsAsFactors = F)) %>%
  dplyr::select(gene_id, Heart, Skin, Max.Specificity.GTEx) %>%
  dplyr::rename(Max.Specificity.GTEx.SMTS = Max.Specificity.GTEx)
# gtex_jssp_smtsd <- as_tibble(read.table(paste0(procDataSubdir, "/GTEx/SMTSD_JSsp_TissueSpecificity.txt"), sep = "\t", header = T, stringsAsFactors = F)) %>%
#   dplyr::select(gene_id, Heart...Left.Ventricle, Heart...Atrial.Appendage, Cells...Transformed.fibroblasts, Max.Specificity.GTEx.SMTSD)
gtex_jssp_smts_942_noSkin <- as_tibble(read.table(paste0(procDataSubdir, "/GTEx/SMTS_JSsp_TissueSpecificity_942_noSkin.txt"), sep = "\t", header = T, stringsAsFactors = F)) %>%
  dplyr::select(gene_id, Heart, GM00942, Max.Specificity.GTEx) %>%
  dplyr::rename(Max.Specificity.GTEx.SMTS_942_noSkin = Max.Specificity.GTEx,
                Heart_942_noSkin = Heart,
                GM00942_noSkin = GM00942) 

gtex_jssp_smts_iCard_noHeart <- as_tibble(read.table(paste0(procDataSubdir, "/GTEx/SMTS_JSsp_TissueSpecificity_iCard_noHeart.txt"), sep = "\t", header = T, stringsAsFactors = F)) %>%
  dplyr::select(gene_id, iCard, Skin, Max.Specificity.GTEx) %>%
  dplyr::rename(Max.Specificity.GTEx.SMTS_iCard_noHeart = Max.Specificity.GTEx,
                Skin_iCard_noHeart = Skin,
                iCard_iCard_noHeart = iCard)

gtex_jssp_smts_iCard_fibs <- as_tibble(read.table(paste0(procDataSubdir, "/GTEx/SMTS_JSsp_TissueSpecificity_iCard_fibs.txt"), sep = "\t", header = T, stringsAsFactors = F)) %>%
  dplyr::select(gene_id, iCard, GM00942, Skin, Heart, Max.Specificity.GTEx) %>%
  dplyr::rename(Max.Specificity.GTEx.SMTS_iCard_fibs = Max.Specificity.GTEx,
                GM00942_iCard_fibs = GM00942,
                iCard_iCard_fibs = iCard,
                Skin_iCard_fibs = Skin,
                Heart_iCard_fibs = Heart) 

gtex_jssp_smts_iCard_fibs_noHeart_noSkin <- as_tibble(read.table(paste0(procDataSubdir, "/GTEx/SMTS_JSsp_TissueSpecificity_iCard_fibs_noHeart_noSkin.txt"), sep = "\t", header = T, stringsAsFactors = F)) %>%
  dplyr::select(gene_id, iCard, GM00942, Max.Specificity.GTEx) %>%
  dplyr::rename(Max.Specificity.GTEx.SMTS_iCard_fibs_noHeart_noSkin = Max.Specificity.GTEx,
                GM00942_iCard_fibs_noHeart_noSkin = GM00942,
                iCard_iCard_fibs_noHeart_noSkin = iCard) 


gtex_jssp = inner_join(gtex_jssp_smts, gtex_jssp_smts_942_noSkin, by = "gene_id") %>%
  inner_join(gtex_jssp_smts_iCard_noHeart, by = "gene_id") %>%
  inner_join(gtex_jssp_smts_iCard_fibs, by = 'gene_id') %>%
  inner_join(gtex_jssp_smts_iCard_fibs_noHeart_noSkin, by = 'gene_id')

# perturbability calculated on individual samples
all_rpm_readFilt_manualFilt_variability_metrics <- as_tibble(read.table(
  paste0(procDataSubdir, "/allExperiments/all_rpm_readFilt_manualFilt_variability_metrics.txt"),
  sep = "\t", stringsAsFactors = F, header = T))

## bind tables
all_rpm_readFilt_manualFilt_variability_metrics_GTEx <- all_rpm_readFilt_manualFilt_variability_metrics %>%
  inner_join(gtex_jssp, by = "gene_id")
all_rpm_readFilt_manualFilt_variability_metrics_GTEx[is.na(all_rpm_readFilt_manualFilt_variability_metrics_GTEx)] <- 0


tfTab_GS <- inner_join(tfTab, geneIDtoGeneSymbol)

library(biomaRt)
ensembl <- useMart('ensembl')
ensembl.human <- useEnsembl(biomart="ensembl", dataset='hsapiens_gene_ensembl', version = 102)
ensembl.mouse <- useEnsembl(biomart="ensembl", dataset='mmusculus_gene_ensembl', version = 102)

HsMm_tfOrthologs<-getLDS(attributes = c('hgnc_symbol'),
                         filters = 'hgnc_symbol', values = tfTab_GS$GeneSymbol, mart = ensembl.human,
                         attributesL = c('mgi_symbol') , martL = ensembl.mouse) 
HsMm_tfOrthologs_forJoin <- HsMm_tfOrthologs %>%
  as_tibble() %>%
  dplyr::rename(GeneSymbol = HGNC.symbol,
                first_Gene_Symbol = MGI.symbol)

pooled_sh_borkent_sumPerGene <- inner_join(pooled_sh_borkent, HsMm_tfOrthologs_forJoin) %>%
  group_by(first_Gene_Symbol) %>%
  filter(Overall_Max == max(Overall_Max))

pooled_sh_borkent_sumPerGene_multiSHenhance <- inner_join(pooled_sh_borkent, HsMm_tfOrthologs_forJoin) %>%
  group_by(first_Gene_Symbol) %>%
  filter(sum(Overall_Max > 0) > 1,
         Overall_Max == max(Overall_Max))

pooled_sh_borkent_tf_responsiveness_all <- pooled_sh_borkent_sumPerGene %>%
  inner_join(all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>% filter(cellType == 'GM00942') %>% inner_join(geneIDtoGeneSymbol), by = 'GeneSymbol')

pooled_sh_borkent_tf_responsiveness_multiSHenhance <- pooled_sh_borkent_sumPerGene_multiSHenhance %>%
  inner_join(all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>% filter(cellType == 'GM00942') %>% inner_join(geneIDtoGeneSymbol), by = 'GeneSymbol')

if(!dir.exists(paste0(graphDir, "/hiFT"))){
  dir.create(paste0(graphDir, "/hiFT"))
}
if(!dir.exists(paste0(graphDir, "/hiFT/borkent_2016"))){
  dir.create(paste0(graphDir, "/hiFT/borkent_2016"))
}

set.seed(38923)
borkent_tfOnly_responsiveUpVsUnresponsiveUp <- ggplot() +
  geom_boxplot(data = pooled_sh_borkent_tf_responsiveness_all %>%
               filter(meanTPM > 50, 
                      nDESeq2conditionsUp <= 1 | nDESeq2conditionsUp >= 5) %>%
               mutate(isResponsive = case_when(
                 nDESeq2conditionsUp <= 1 ~ F,
                 nDESeq2conditionsUp >= 5 ~ T
               )),
             aes(isResponsive, Overall_Max), alpha = 0.2) +
  geom_jitter(data = pooled_sh_borkent_tf_responsiveness_all %>%
                 filter(meanTPM > 50, 
                        nDESeq2conditionsUp <= 1 | nDESeq2conditionsUp >= 5) %>%
                 mutate(isResponsive = case_when(
                   nDESeq2conditionsUp <= 1 ~ F,
                   nDESeq2conditionsUp >= 5 ~ T
                 )),
               aes(isResponsive, Overall_Max), height = 0, width = 0.3, alpha = 0.2) +
  xlab("Human TF gene responsive?") +
  ylab("Maximum log2(fold change)\nafter knockdown of mouse ortholog") +
  theme_classic()
ggsave(borkent_tfOnly_responsiveUpVsUnresponsiveUp, file = paste0(graphDir, "/hiFT/borkent_2016/borkent_tfOnly_responsiveUpVsUnresponsiveUp.pdf"), width = 4, height = 8, useDingbats = F)  

borkent_tfOnly_responsiveUpVsDownVsUnresponsive <- ggplot() +
  geom_boxplot(data = pooled_sh_borkent_tf_responsiveness_all %>%
                 filter(meanTPM > 50, 
                        nDESeq2conditionsAll <= 1 | nDESeq2conditionsUp >= 5 | nDESeq2conditionsDown >= 5) %>%
                 mutate(isResponsive = case_when(
                   nDESeq2conditionsAll <= 1 ~ 'Unresponsive',
                   nDESeq2conditionsUp >= 5 & nDESeq2conditionsUp/nDESeq2conditionsAll >= 0.6 ~ 'Responsive up',
                   nDESeq2conditionsDown >= 5 & nDESeq2conditionsUp/nDESeq2conditionsAll <= 0.4 ~ 'Responsive down',
                   nDESeq2conditionsUp >= 5 & nDESeq2conditionsUp/nDESeq2conditionsAll > 0.4 &nDESeq2conditionsUp/nDESeq2conditionsAll < 0.6   ~ 'Responsive both'
                 )),
               aes(isResponsive, Overall_Max), alpha = 0.2) +
  geom_jitter(data = pooled_sh_borkent_tf_responsiveness_all %>%
                 filter(meanTPM > 50, 
                        nDESeq2conditionsUp <= 1 | nDESeq2conditionsUp >= 5 | nDESeq2conditionsDown >= 5) %>%
                 mutate(isResponsive = case_when(
                   nDESeq2conditionsAll <= 1 ~ 'Unresponsive',
                   nDESeq2conditionsUp >= 5 & nDESeq2conditionsUp/nDESeq2conditionsAll >= 0.6 ~ 'Responsive up',
                   nDESeq2conditionsDown >= 5 & nDESeq2conditionsDown/nDESeq2conditionsAll >= 0.6 ~ 'Responsive down',
                   nDESeq2conditionsUp >= 5 & nDESeq2conditionsUp/nDESeq2conditionsAll > 0.4 &nDESeq2conditionsUp/nDESeq2conditionsAll < 0.6   ~ 'Responsive both'
                 )),
               aes(isResponsive, Overall_Max), height = 0, width = 0.2, alpha = 0.2) +
  xlab("Human TF gene responsive?") +
  ylab("Maximum log2(fold change)\nafter knockdown of mouse ortholog") +
  theme_classic()
ggsave(borkent_tfOnly_responsiveUpVsDownVsUnresponsive, file = paste0(graphDir, "/hiFT/borkent_2016/borkent_tfOnly_responsiveUpVsDownVsUnresponsive.pdf"), width = 4, height = 8, useDingbats = F)  

specFiltHi = 0.24
specFiltLo = 0.24
borkent_tfOnly_responsiveUpVsSpecific <- ggplot() +
  geom_boxplot(data = pooled_sh_borkent_tf_responsiveness_all %>%
                 filter(meanTPM > 50, 
                        nDESeq2conditionsUp <= 1 | nDESeq2conditionsUp >= 5) %>%
                 mutate(isResponsive = case_when(
                   nDESeq2conditionsUp <= 1 & iCard_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Unresponsive\nNon-specific',
                   nDESeq2conditionsUp >= 5 & iCard_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Responsive\nNon-specific',
                   nDESeq2conditionsUp <= 1 & iCard_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Unresponsive\niPS-CM-specific',
                   nDESeq2conditionsUp >= 5 & iCard_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Responsive\niPS-CM-pecific'
                 )),
               aes(isResponsive, Overall_Max), alpha = 0.2) +
  geom_jitter(data = pooled_sh_borkent_tf_responsiveness_all %>%
                filter(meanTPM > 50, 
                       nDESeq2conditionsUp <= 1 | nDESeq2conditionsUp >= 5 | nDESeq2conditionsDown >= 5) %>%
                mutate(isResponsive = case_when(
                  nDESeq2conditionsUp <= 1 & iCard_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Unresponsive\nNon-specific',
                  nDESeq2conditionsUp >= 5 & iCard_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Responsive\nNon-specific',
                  nDESeq2conditionsUp <= 1 & iCard_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Unresponsive\niPS-CM-specific',
                  nDESeq2conditionsUp >= 5 & iCard_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Responsive\niPS-CM-pecific'
                )),
              aes(isResponsive, Overall_Max), height = 0, width = 0.2, alpha = 0.2) +
  xlab("Human TF gene responsive and specific?") +
  ylab("Maximum log2(fold change)\nafter knockdown of mouse ortholog") +
  theme_classic()
# ggsave(borkent_tfOnly_responsiveUpVsSpecific, file = paste0(graphDir, "/hiFT/borkent_2016/borkent_tfOnly_responsiveUpVsSpecific.pdf"), width = 4, height = 8, useDingbats = F)  

enrFilt = 2
borkent_tfOnly_responsiveUpVsSpecific_minTPM20 <- ggplot() +
  geom_boxplot(data = pooled_sh_borkent_tf_responsiveness_all %>%
                 filter(meanTPM > 20, 
                        nDESeq2conditionsUp <= 1 | nDESeq2conditionsUp >= 5) %>%
                 mutate(isResponsive = case_when(
                   nDESeq2conditionsUp <= 1 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Unresponsive\nNon-specific',
                   nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Responsive\nNon-specific',
                   nDESeq2conditionsUp <= 1 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Unresponsive\nfibroblast-specific',
                   nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Responsive\nfibroblast-pecific'
                 )),
               aes(isResponsive, Overall_Max), alpha = 0.2) +
  geom_jitter(data = pooled_sh_borkent_tf_responsiveness_all %>%
                filter(meanTPM > 20, 
                       nDESeq2conditionsUp <= 1 | nDESeq2conditionsUp >= 5 | nDESeq2conditionsDown >= 5) %>%
                mutate(isResponsive = case_when(
                  nDESeq2conditionsUp <= 1 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Unresponsive\nNon-specific',
                  nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Responsive\nNon-specific',
                  nDESeq2conditionsUp <= 1 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Unresponsive\nfibroblast-specific',
                  nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Responsive\nfibroblast-pecific'
                )),
              aes(isResponsive, Overall_Max), height = 0, width = 0.2, alpha = 0.2) +
  geom_text(data = pooled_sh_borkent_tf_responsiveness_all %>%
              filter(meanTPM > 20, 
                     nDESeq2conditionsUp <= 1 | nDESeq2conditionsUp >= 5 | nDESeq2conditionsDown >= 5) %>%
              mutate(isResponsive = case_when(
                nDESeq2conditionsUp <= 1 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Unresponsive\nNon-specific',
                nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Responsive\nNon-specific',
                nDESeq2conditionsUp <= 1 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Unresponsive\nfibroblast-specific',
                nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Responsive\nfibroblast-pecific'
              )) %>% ungroup() %>%
              group_by(isResponsive) %>%
              summarise(fracEnriched = sum(Overall_Max >= enrFilt)/length(Overall_Max)) %>%
              mutate(Overall_Max = 7),
            aes(isResponsive, Overall_Max, label = paste0(as.character(round(100*fracEnriched, 3)), '%'))) +
  geom_hline(data = NULL, aes(yintercept = enrFilt), linetype = "dashed", color = 'blue') +
  xlab("Human TF gene responsive and specific?") +
  ylab("Maximum log2(fold change)\nafter knockdown of mouse ortholog") +
  theme_classic()
# ggsave(borkent_tfOnly_responsiveUpVsSpecific_minTPM20, file = paste0(graphDir, "/hiFT/borkent_2016/borkent_tfOnly_responsiveUpVsSpecific_minTPM20.pdf"), width = 3.5, height = 3, useDingbats = F)  


enrFilt = 2
borkent_tfOnly_responsiveUpVsSpecific_minTPM20_UnrLE4 <- ggplot() +
  geom_boxplot(data = pooled_sh_borkent_tf_responsiveness_all %>%
                 filter(meanTPM > 20, 
                        nDESeq2conditionsUp <= 4 | nDESeq2conditionsUp >= 5) %>%
                 mutate(isResponsive = case_when(
                   nDESeq2conditionsUp <= 4 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Unresponsive\nNon-specific',
                   nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Responsive\nNon-specific',
                   nDESeq2conditionsUp <= 4 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Unresponsive\nfibroblast-specific',
                   nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Responsive\nfibroblast-pecific'
                 )),
               aes(isResponsive, Overall_Max), alpha = 0.2) +
  geom_jitter(data = pooled_sh_borkent_tf_responsiveness_all %>%
                filter(meanTPM > 20, 
                       nDESeq2conditionsUp <= 4 | nDESeq2conditionsUp >= 5) %>%
                mutate(isResponsive = case_when(
                  nDESeq2conditionsUp <= 4 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Unresponsive\nNon-specific',
                  nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Responsive\nNon-specific',
                  nDESeq2conditionsUp <= 4 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Unresponsive\nfibroblast-specific',
                  nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Responsive\nfibroblast-pecific'
                )),
              aes(isResponsive, Overall_Max), height = 0, width = 0.2, alpha = 0.2) +
  geom_text(data = pooled_sh_borkent_tf_responsiveness_all %>%
              filter(meanTPM > 20, 
                     nDESeq2conditionsUp <= 4 | nDESeq2conditionsUp >= 5) %>%
              mutate(isResponsive = case_when(
                nDESeq2conditionsUp <= 4 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Unresponsive\nNon-specific',
                nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Responsive\nNon-specific',
                nDESeq2conditionsUp <= 4 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Unresponsive\nfibroblast-specific',
                nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Responsive\nfibroblast-pecific'
              )) %>% ungroup() %>%
              group_by(isResponsive) %>%
              summarise(fracEnriched = sum(Overall_Max >= enrFilt)/length(Overall_Max)) %>%
              mutate(Overall_Max = 7),
            aes(isResponsive, Overall_Max, label = paste0(as.character(round(100*fracEnriched, 3)), '%'))) +
  geom_hline(data = NULL, aes(yintercept = enrFilt), linetype = "dashed", color = 'blue') +
  xlab("Human TF gene responsive and specific?") +
  ylab("Maximum log2(fold change)\nafter knockdown of mouse ortholog") +
  theme_classic()
# ggsave(borkent_tfOnly_responsiveUpVsSpecific_minTPM20_UnrLE4, file = paste0(graphDir, "/hiFT/borkent_2016/borkent_tfOnly_responsiveUpVsSpecific_minTPM20.pdf"), width = 3.5, height = 3, useDingbats = F)  

set.seed(4286)
fibro_GTEx_indiv_nDEup_vs_JSsp_GM00942_noSkinAndHeart_flip_onlyBorkent <- ggplot() + 
  geom_point(data = pooled_sh_borkent_tf_responsiveness_all %>%
               filter(meanTPM > 20, 
                      nDESeq2conditionsUp <= 4 | nDESeq2conditionsUp >= 5) , 
             aes(GM00942_iCard_fibs_noHeart_noSkin, nDESeq2conditionsUp), alpha = 0.2) +
  # geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #              filter(cellType == "GM00942") %>%
  #              inner_join(hiFT_factorSum, by = "gene_id") %>%
  #              filter(isBarrier == 1), 
  #            aes(GM00942_iCard_fibs_noHeart_noSkin, nDESeq2conditionsUp), color = "forestgreen") +
  # geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #                   filter(cellType == "GM00942") %>%
  #                   inner_join(hiFT_factorSum, by = "gene_id") %>%
  #                   filter(isBarrier == 1), 
  #                 aes(GM00942_iCard_fibs_noHeart_noSkin, nDESeq2conditionsUp, label = GeneSymbol), color = "forestgreen") +
  theme_bw() + 
  # ylim(c(0,0.5)) +
  ggtitle('Fibroblast TF perturbability vs GM942 average specificity\nExcludes GTEx Heart and Skin\nOnly factors in Borkent screen')
# plot(fibro_GTEx_indiv_nDEup_vs_JSsp_GM00942_noSkinAndHeart_flip_onlyBarriers)
# ggsave(fibro_GTEx_indiv_nDEup_vs_JSsp_GM00942_noSkinAndHeart_flip_onlyBorkent, file = paste0(graphSubdir, "/hiFT/borkent_2016/fibro_GTEx_indiv_nDEup_vs_JSsp_GM00942_noSkinHeart_flip_onlyBorkent.pdf"), width = 8, height = 8, useDingbats = F)



borkent_tfOnly_responsiveUpVsSpecific_minRPM20 <- ggplot() +
  geom_boxplot(data = pooled_sh_borkent_tf_responsiveness_all %>%
                 filter(meanRPM > 20, 
                        nDESeq2conditionsUp <= 1 | nDESeq2conditionsUp >= 5) %>%
                 mutate(isResponsive = case_when(
                   nDESeq2conditionsUp <= 1 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Unresponsive\nNon-specific',
                   nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Responsive\nNon-specific',
                   nDESeq2conditionsUp <= 1 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Unresponsive\nfibroblast-specific',
                   nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Responsive\nfibroblast-pecific'
                 )),
               aes(isResponsive, Overall_Max), alpha = 0.2) +
  geom_jitter(data = pooled_sh_borkent_tf_responsiveness_all %>%
                filter(meanRPM > 20, 
                       nDESeq2conditionsUp <= 1 | nDESeq2conditionsUp >= 5 | nDESeq2conditionsDown >= 5) %>%
                mutate(isResponsive = case_when(
                  nDESeq2conditionsUp <= 1 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Unresponsive\nNon-specific',
                  nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Responsive\nNon-specific',
                  nDESeq2conditionsUp <= 1 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Unresponsive\nfibroblast-specific',
                  nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Responsive\nfibroblast-pecific'
                )),
              aes(isResponsive, Overall_Max), height = 0, width = 0.2, alpha = 0.2) +
  geom_text(data = pooled_sh_borkent_tf_responsiveness_all %>%
              filter(meanRPM > 20, 
                     nDESeq2conditionsUp <= 1 | nDESeq2conditionsUp >= 5 | nDESeq2conditionsDown >= 5) %>%
              mutate(isResponsive = case_when(
                nDESeq2conditionsUp <= 1 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Unresponsive\nNon-specific',
                nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Responsive\nNon-specific',
                nDESeq2conditionsUp <= 1 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Unresponsive\nfibroblast-specific',
                nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Responsive\nfibroblast-pecific'
              )) %>% ungroup() %>%
              group_by(isResponsive) %>%
              summarise(fracEnriched = sum(Overall_Max >= enrFilt)/length(Overall_Max)) %>%
              mutate(Overall_Max = 7),
            aes(isResponsive, Overall_Max, label = paste0(as.character(round(100*fracEnriched, 3)), '%'))) +
  geom_hline(data = NULL, aes(yintercept = enrFilt), linetype = "dashed", color = 'blue') +
  xlab("Human TF gene responsive and specific?") +
  ylab("Maximum log2(fold change)\nafter knockdown of mouse ortholog") +
  theme_classic()

enrFilt = 2.5
borkent_tfOnly_responsiveUpVsSpecificVsUnresponsive_minRPM20 <- ggplot() +
  geom_boxplot(data = pooled_sh_borkent_tf_responsiveness_all %>%
                 filter(meanRPM > 20, 
                        nDESeq2conditionsUp <= 1 | nDESeq2conditionsUp >= 5) %>%
                 mutate(isResponsive = case_when(
                   nDESeq2conditionsUp <= 2 & iCard_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Unresponsive\nNon-specific',
                   nDESeq2conditionsUp >= 5 & iCard_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Responsive\nNon-specific',
                   nDESeq2conditionsUp <= 2 & iCard_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Unresponsive\niPS-CM-specific',
                   nDESeq2conditionsUp >= 5 & iCard_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Responsive\niPS-CM-pecific'
                 )),
               aes(isResponsive, Overall_Max), alpha = 0.2) +
  geom_jitter(data = pooled_sh_borkent_tf_responsiveness_all %>%
                filter(meanRPM > 20, 
                       nDESeq2conditionsUp <= 1 | nDESeq2conditionsUp >= 5) %>%
                mutate(isResponsive = case_when(
                  nDESeq2conditionsUp <= 2 & iCard_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Unresponsive\nNon-specific',
                  nDESeq2conditionsUp >= 5 & iCard_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Responsive\nNon-specific',
                  nDESeq2conditionsUp <= 2 & iCard_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Unresponsive\niPS-CM-specific',
                  nDESeq2conditionsUp >= 5 & iCard_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Responsive\niPS-CM-pecific'
                )),
              aes(isResponsive, Overall_Max), height = 0, width = 0.2, alpha = 0.2) +
  geom_text(data = pooled_sh_borkent_tf_responsiveness_all %>%
              filter(meanRPM > 20, 
                     nDESeq2conditionsUp <= 1 | nDESeq2conditionsUp >= 5) %>%
              mutate(isResponsive = case_when(
                nDESeq2conditionsUp <= 2 & iCard_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Unresponsive\nNon-specific',
                nDESeq2conditionsUp >= 5 & iCard_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Responsive\nNon-specific',
                nDESeq2conditionsUp <= 2 & iCard_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Unresponsive\niPS-CM-specific',
                nDESeq2conditionsUp >= 5 & iCard_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Responsive\niPS-CM-pecific'
              )) %>% ungroup() %>%
              group_by(isResponsive) %>%
              summarise(fracEnriched = sum(Overall_Max >= enrFilt)/length(Overall_Max)) %>%
              mutate(Overall_Max = 7),
            aes(isResponsive, Overall_Max, label = paste0(as.character(round(100*fracEnriched, 3)), '%'))) +
  geom_hline(data = NULL, aes(yintercept = enrFilt), linetype = "dashed", color = 'blue') +
  xlab("Human TF gene responsive and specific?") +
  ylab("Maximum log2(fold change)\nafter knockdown of mouse ortholog") +
  theme_classic()


### Sliding threshold analysis

minUpResp = 5
enrFilt = 2

filts = seq(0.1,3,0.1)
slidingRes = list()
for (fil in filts) {
  enrFiltT = fil
  tempResSlide <- pooled_sh_borkent_tf_responsiveness_all %>%
  dplyr::select(GeneSymbol, gene_id, meanTPM, meanRPM, nDESeq2conditionsUp, GM00942_iCard_fibs_noHeart_noSkin, Overall_Max) %>%
  filter(meanTPM > 30) %>%
  mutate(isResponsive = case_when(
    nDESeq2conditionsUp < minUpResp ~ 'Unresponsive',
    nDESeq2conditionsUp >= minUpResp ~ 'Responsive'
  )) %>%
  ungroup() %>%
  group_by(isResponsive) %>%
  summarise(fracEnriched = sum(Overall_Max >= enrFiltT)/length(Overall_Max)) %>%
  spread(isResponsive, fracEnriched) %>%
  mutate(absChange = Responsive - Unresponsive,
         foldChange = Responsive / Unresponsive,
         fracChange = absChange / Unresponsive,
         enrFilt = enrFiltT)
 if(is.null(dim(slidingRes))){
   slidingRes <- tempResSlide
 } else {
   slidingRes %<>% bind_rows(tempResSlide)
 }
}

slidingRes_diffs <- slidingRes %>%
  group_by(enrFilt,Responsive,Unresponsive) %>%
  gather(metric, value, 3:5)

borkent_fracEnriched_perFilter <- ggplot() +
  geom_point(data = slidingRes, aes(enrFilt, Responsive), color = 'orange') +
  geom_point(data = slidingRes, aes(enrFilt, Unresponsive)) +
  geom_vline(data = NULL, aes(xintercept = enrFilt), linetype = 'dashed', color = 'blue') +
  ylab('Fraction enriched') +
  theme_classic() +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), limits = c(0,0.5))
ggsave(borkent_fracEnriched_perFilter, file = paste0(graphDir, "/hiFT/borkent_2016/borkent_fracEnriched_perFilter_minTPM30.pdf"), width = 3, height = 1.5, useDingbats = F)  

borkent_foldChangeEnriched_perFilter <- ggplot() +
  geom_point(data = slidingRes, aes(enrFilt, foldChange), shape = 4) +
  geom_vline(data = NULL, aes(xintercept = enrFilt), linetype = 'dashed', color = 'blue') +
  ylab('Responsive/Unresponsive') +
  theme_classic() +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), limits = c(0,5.2))
ggsave(borkent_foldChangeEnriched_perFilter, file = paste0(graphDir, "/hiFT/borkent_2016/borkent_foldChangeEnriched_perFilter_minTPM30.pdf"), width = 3, height = 1.5, useDingbats = F)  

borkent_absChangeEnriched_perFilter <- ggplot() +
  geom_point(data = slidingRes, aes(enrFilt, absChange), shape = 3) +
  geom_vline(data = NULL, aes(xintercept = enrFilt), linetype = 'dashed', color = 'blue') +
  ylab('Responsive - Unresponsive') +
  theme_classic() +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))
ggsave(borkent_absChangeEnriched_perFilter, file = paste0(graphDir, "/hiFT/borkent_2016/borkent_absChangeEnriched_perFilter_minTPM30.pdf"), width = 3, height = 1.5, useDingbats = F)  

# ggplot(slidingRes, aes(enrFilt, foldChange)) +
#   geom_point() +
#   ylim(c(0,5))
# 
# ggplot(slidingRes, aes(enrFilt, absChange)) +
#   geom_point()
#   
enrFilt = 2
set.seed(17182)
borkent_tfOnly_responsiveUpVsSpecific_minTPM30_UnrLE4 <- ggplot() +
  geom_boxplot(data = pooled_sh_borkent_tf_responsiveness_all %>%
                 filter(meanTPM > 30, 
                        nDESeq2conditionsUp <= 4 | nDESeq2conditionsUp >= 5) %>%
                 mutate(isResponsive = case_when(
                   nDESeq2conditionsUp <= 4 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'F_Unresponsive\nA_Non-specific',
                   nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Responsive\nA_Non-specific',
                   nDESeq2conditionsUp <= 4 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'F_Unresponsive\nfibroblast-specific',
                   nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Responsive\nfibroblast-pecific'
                 )),
               aes(isResponsive, Overall_Max), alpha = 0.2, outlier.shape = NA) +
  geom_jitter(data = pooled_sh_borkent_tf_responsiveness_all %>%
                filter(meanTPM > 30,
                       nDESeq2conditionsUp <= 4 | nDESeq2conditionsUp >= 5) %>%
                mutate(isResponsive = case_when(
                  nDESeq2conditionsUp <= 4 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'F_Unresponsive\nA_Non-specific',
                  nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Responsive\nA_Non-specific',
                  nDESeq2conditionsUp <= 4 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'F_Unresponsive\nfibroblast-specific',
                  nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Responsive\nfibroblast-pecific'
                )),
              aes(isResponsive, Overall_Max), height = 0, width = 0.2, alpha = 0.2) +
  geom_text(data = pooled_sh_borkent_tf_responsiveness_all %>%
              filter(meanTPM > 30, 
                     nDESeq2conditionsUp <= 4 | nDESeq2conditionsUp >= 5) %>%
              mutate(isResponsive = case_when(
                nDESeq2conditionsUp <= 4 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'F_Unresponsive\nA_Non-specific',
                nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin <= specFiltLo ~ 'Responsive\nA_Non-specific',
                nDESeq2conditionsUp <= 4 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'F_Unresponsive\nfibroblast-specific',
                nDESeq2conditionsUp >= 5 & GM00942_iCard_fibs_noHeart_noSkin >= specFiltHi ~ 'Responsive\nfibroblast-pecific'
              )) %>% ungroup() %>%
              group_by(isResponsive) %>%
              summarise(fracEnriched = sum(Overall_Max >= enrFilt)/length(Overall_Max)) %>%
              mutate(Overall_Max = 7),
            aes(isResponsive, Overall_Max, label = paste0(as.character(round(100*fracEnriched, 3)), '%'))) +
  geom_hline(data = NULL, aes(yintercept = enrFilt), linetype = "dashed", color = 'blue') +
  xlab("Human TF gene responsive and specific?") +
  ylab("Maximum log2(fold change)\nafter knockdown of mouse ortholog") +
  theme_classic()
ggsave(borkent_tfOnly_responsiveUpVsSpecific_minTPM30_UnrLE4, file = paste0(graphDir, "/hiFT/borkent_2016/borkent_tfOnly_responsiveUpVsSpecific_minTPM30.pdf"), width = 3.5, height = 3, useDingbats = F)  

borkent_responsiveness_noSpec <- pooled_sh_borkent_tf_responsiveness_all %>%
  filter(meanTPM > 30) %>%
  mutate(isResponsive = nDESeq2conditionsUp >= 5)

set.seed(17812)
borkent_tfOnly_responsiveUp_minTPM30_UnrLE4 <- ggplot() +
  geom_boxplot(data = borkent_responsiveness_noSpec,
               aes(isResponsive, Overall_Max), alpha = 0.2, outlier.shape = NA) +
  geom_jitter(data = borkent_responsiveness_noSpec,
              aes(isResponsive, Overall_Max), height = 0, width = 0.2, alpha = 0.2) +
  geom_text(data = borkent_responsiveness_noSpec %>% ungroup() %>%
              group_by(isResponsive) %>%
              summarise(fracEnriched = sum(Overall_Max >= enrFilt)/length(Overall_Max)) %>%
              mutate(Overall_Max = 7),
            aes(isResponsive, Overall_Max, label = paste0(as.character(round(100*fracEnriched, 3)), '%'))) +
  geom_hline(data = NULL, aes(yintercept = enrFilt), linetype = "dashed", color = 'blue') +
  xlab("Human TF gene responsive?") +
  ylab("Maximum log2(fold change)\nafter knockdown of mouse ortholog") +
  theme_classic()
ggsave(borkent_tfOnly_responsiveUp_minTPM30_UnrLE4, file = paste0(graphDir, "/hiFT/borkent_2016/borkent_tfOnly_responsiveUp_minTPM30.pdf"), width = 2.5, height = 2, useDingbats = F)  


set.seed(4286)
fibro_GTEx_indiv_nDEup_vs_JSsp_GM00942_noSkinAndHeart_flip_onlyBorkent_minTPM30 <- ggplot() + 
  geom_jitter(data = pooled_sh_borkent_tf_responsiveness_all %>%
               filter(meanTPM > 30, 
                      nDESeq2conditionsUp <= 4 | nDESeq2conditionsUp >= 5) , 
             aes(GM00942_iCard_fibs_noHeart_noSkin, nDESeq2conditionsUp), alpha = 0.2, height = 0.1) +
  # geom_point(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #              filter(cellType == "GM00942") %>%
  #              inner_join(hiFT_factorSum, by = "gene_id") %>%
  #              filter(isBarrier == 1), 
  #            aes(GM00942_iCard_fibs_noHeart_noSkin, nDESeq2conditionsUp), color = "forestgreen") +
  # geom_text_repel(data = all_rpm_readFilt_manualFilt_variability_metrics_GTEx %>%
  #                   filter(cellType == "GM00942") %>%
  #                   inner_join(hiFT_factorSum, by = "gene_id") %>%
  #                   filter(isBarrier == 1), 
  #                 aes(GM00942_iCard_fibs_noHeart_noSkin, nDESeq2conditionsUp, label = GeneSymbol), color = "forestgreen") +
  theme_bw() + 
  # ylim(c(0,0.5)) +
  ggtitle('Fibroblast TF perturbability vs GM942 average specificity\nExcludes GTEx Heart and Skin\nOnly transcription factors in Borkent screen')
# plot(fibro_GTEx_indiv_nDEup_vs_JSsp_GM00942_noSkinAndHeart_flip_onlyBarriers)
ggsave(fibro_GTEx_indiv_nDEup_vs_JSsp_GM00942_noSkinAndHeart_flip_onlyBorkent_minTPM30, file = paste0(graphSubdir, "/hiFT/borkent_2016/fibro_GTEx_indiv_nDEup_vs_JSsp_GM00942_noSkinHeart_flip_onlyBorkent_minTPM30.pdf"), width = 8, height = 8, useDingbats = F)

