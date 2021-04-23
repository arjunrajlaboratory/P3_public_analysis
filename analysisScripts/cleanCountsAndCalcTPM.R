#!/usr/bin/env Rscript

# cleanCountsAndCalcTPM.R loads all raw count data, quality-filters it, and transforms to TPM.
# 	quality filters:
# 		- no sample with less than 60% mapping
# 		- no gene with median TPM per 96-well plate <= 0.25
# 	generates: 
#		- allMappedCounts: summarizedExperiment object saved in procDataDir/RNAtagSeq/allMappedCounts/
#				genes x samples table of HTseq-mapped counts over all samples in all RNAtag-seq runs
#				includes all genes
#				includes all samples
#		- allQcFilteredTPMs: summarizedExperiment object saved in procDataDir/RNAtagSeq/allQcFilteredTPMs/
#				genes x samples table of TPM values of QC-filtered genes and samples over all RNAtag-seq runs
#   - for each 96-well plate, $PLATENAME/mappedCounts: summarizedExperiment object saved in procDataDir/RNAtagSeq/$PLATENAME/mappedCounts/
#				genes x samples table of HTseq-mapped counts
#				includes all genes
#				includes all samples
#   - for each 96-well plate, $PLATENAME/qcFilteredTPMs: summarizedExperiment object saved in procDataDir/RNAtagSeq/$PLATENAME/qcFilteredTPMs/
#				genes x samples table of TPM values of QC-filtered genes and samples

# don't run lines 22-32 if running from R workspace
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

# if running manually in Rstudio, start here and use:
# projectDir = '~/Dropbox (RajLab)/Shared_IanM/cellid_201807_onward/'
# procDataSubdir = 'procDataScripted'
# graphSubdir = 'graphsScripted'

procDataDir = paste0(projectDir, procDataSubdir)
graphDir = paste0(projectDir, graphSubdir)

setwd(projectDir)
cat("Starting cleanCountsAndCalcTPM.R...\n")
cat(paste0("Working in ", getwd(), "\n"))


# load relevant libraries and source custom functions
# comment out if running analyzeAll.sh instead of just this script (analyzeAll.sh loads all libraries)
library(DESeq2)
# library(HDF5Array)
library(tidyverse)
# library(rtracklayer)
source(paste0(projectDir, "analysisScripts/RNAtagUtilities.R"))

setwd(projectDir)

# load generic metadata
cat("Loading generic metadata.\n")
geneIDlengthsFile = "annotations/hg19_lengths_per_gene_id.tsv"                  # gene lengths
geneIDfile = "annotations/hg19gene_idToGeneSymbol.tsv"                          # gene symbols
tfGeneIDfile = "annotations/TF_gene_ids.csv"                                    # gene IDs corresponding to transcription factors
drugIDfile = "metadata/workingSuppTable_drugIDsWithTargetsAndAttributes.txt"    # drug IDs and corresponding reference information, e.g., drug name, main targets, medianIC50

geneIDtoGeneSymbol <- as_tibble(read.csv(geneIDfile, sep='\t', stringsAsFactors = F))
geneIDlengths <- as_tibble(read.csv(geneIDlengthsFile, sep='\t', stringsAsFactors = F))
drugIDtable <- as_tibble(read.csv(drugIDfile, sep="\t", stringsAsFactors = F)) %>%
  mutate(drug = drugID)
cat("Generic metadata loaded.\n")

cat("Assembling full geneData (rowData) table...\n")
geneData = full_join(geneIDlengths, geneIDtoGeneSymbol)
geneData = as.data.frame(geneData)
rownames(geneData)<-geneData$gene_id
cat("geneData assembled.\n")

# for single-drug samples
# load count data
runList = c(#"RNAtagTest2"
             #"RNAtagTest3"
             #"iCardTestA"
              "RNAtagTest4",
              "RNAtagTest5",
              "iCardTest4set100x",
              "iCardTest5set100x"
            )

countDir = "extractedData/RNAtagSeq/"

countFileName = "htseq_meltedData.tsv"

cat("Looping through the following runs:\n")
cat(runList)
counts = list()
for (run in runList){
  
  cat(paste0("Working on ", run, "...\n"))
  
  cat("Loading run-specific metadata...\n")
  # load per-run metadata
  
  wellLayoutFile = paste0("metadata/", run, "/wellIDpoolTagMapping.csv")      # columns: rowLet, colNum, poolID, TagID 
  wellContentFile = paste0("metadata/", run, "/conditionsPositions_tall.csv") # Sample data per well: drug, dose, if it's HeLa
  wellManualInspectionFile = paste0("metadata/", run, "/transImageNotes.txt") # Manual inspection of interior wells
  # note: cell type is inferred from plate ID; iCard/GM00942 grown on separate plates.
  
  # check for custom layout, otherwise assume a default
  if (file.exists(wellLayoutFile)){
    wellLayoutTab <- as_tibble(read.csv(wellLayoutFile, stringsAsFactors = F))
  } else {
    wellLayoutTab = as_tibble(data.frame(
      rowLet = rep(LETTERS[1:8], times = 12),
      colNum = rep(1:12, each = 8),
      poolID = paste0(run, "-pool", as.character(rep(1:3,each = 32))),
      tagID = rep(c(paste0("Tag0", as.character(1:9)), paste0("Tag", as.character(10:32))), times = 3)
    ))
  }
  #edit the following out when updated...
  wellLayoutTab2 = wellLayoutTab %>% 
    separate(poolID, into = c("run", "pID"), sep = "\\-", remove = T) %>% 
    mutate(runID = run) %>% 
    unite(poolID, c(runID, pID), sep = "-") %>% 
    dplyr::select(-run) 
  
  wellContentTab <- as_tibble(read.csv(wellContentFile, stringsAsFactors = F)) %>%
    separate(platePosition, into = c("rowLet", "colNum"), sep = "\\_") %>%
    mutate(colNum = as.integer(colNum))
  
  wellManualTab <- as_tibble(read.table(wellManualInspectionFile, sep = "\t", header = T, stringsAsFactors = F)) %>%
    separate(platePosition, into = c("rowLet", "colNum"), sep = "\\_") %>%
    mutate(colNum = as.integer(colNum))
  colnames(wellManualTab) = c("imageID", "rowLet", "colNum", paste0(colnames(wellManualTab)[4:ncol(wellManualTab)], "_", "ManualFileTitle"))
  
  wellPositions = inner_join(wellLayoutTab2, wellContentTab, by = c("rowLet", "colNum")) %>%
    mutate(experiment = run) %>%
    unite(combID, c(experiment, poolID, tagID), sep = "_", remove = T) %>%
    full_join(wellManualTab, by = c("rowLet", "colNum"))
  
  cat("Metadata loaded. Loading count data...\n")
  # load run count data
  tempFile = paste0(countDir, run, '/', countFileName)
  rawCounts = as_tibble(read.table(tempFile, header = T))
  
  cat("Count data loaded. Checking mapping stats...\n")
  # check mapping stats
  unmappingResults = unmappedRNAtagReadsPerSample(rawCounts)
  unmapDat = unmappingResults$unmappingSum %>%
    filter(readStatus == "unmappedPerc") %>%
    mutate(unmappedPercentageOfTotal = percentage) %>%
    dplyr::select(-c(readStatus, percentage))
  
  cat("Mapping stats done. Checking counts per sample...\n")
  # check reads per sample
  countsPerSample = rawCounts %>%
    group_by(experiment, sampleID, tagID) %>%
    summarise(totalCounts = sum(counts)) %>% ungroup() %>%
    unite(combID, experiment:tagID, sep = "_", remove = T)
  
  cat("Counts per sample done. Assembing sampleData (colData)...\n")
  
  tempDrugIDfile = paste0("metadata/", run, "/", run, "_drugs_stockPositions.tsv")
  # drugIDtable <- as_tibble(read.csv(tempDrugIDfile, sep="\t")) %>%
    # mutate(drug = drugIDforPlating)
  
  
  sampleData = inner_join(countsPerSample, wellPositions, by = "combID") %>%
    mutate(foldMedianIC50 = "100x", 
           drug = ifelse(grepl("DMSO", as.character(condition)), "DMSO", 
                         ifelse(grepl("PBS", condition), "PBS", as.character(condition))),
           cellType = ifelse(grepl("iCard", combID), "iCard", 
                             ifelse(grepl("HeLa", as.character(condition)), "HeLa", "GM00942")),
           finalDMSOconcentration.percentByVolume = 0.28,
           treatmentDuration.hours = 96) %>%
    left_join(drugIDtable, by = "drug") %>%
    mutate(finalDrugConcentration.nM = ifelse(is.na(medianIC50..nM.), 0, medianIC50..nM. * 100)) %>%
    separate(combID, into = c("experiment","poolID", "tagID"), sep = "\\_", remove = F) %>%
    inner_join(unmapDat, by = c("experiment","poolID", "tagID")) %>%
    dplyr::select(-c(drugID_ManualFileTitle,
                     Class_ManualFileTitle,
                     Targets_ManualFileTitle, 
                     medianIC50..nM._ManualFileTitle,
                     InfoFromSelleckchem_ManualFileTitle))
  sampleData$Product.Name = as.character(sampleData$Product.Name)
  sampleData$Product.Name[is.na(sampleData$Product.Name)] = "control"
  sampleData = as.data.frame(sampleData)
  rownames(sampleData) = sampleData$combID
  sampleData$combID <- NULL
  
  cat("sampleData assembled. Saving run-specific mappedCounts summarizedExperiment object...\n")
  mappedCounts = rawCounts %>%
    filter(gene_id != "__alignment_not_unique" & gene_id != "__ambiguous" & gene_id != "__no_feature" & gene_id != "__not_aligned" & gene_id != "__too_low_aQual") %>%
    filter(tagID != "Unassigned") %>%
    unite(combID, experiment:tagID, sep = "_", remove = T) %>%
    spread(combID, counts)
  mappedCounts<-as.data.frame(mappedCounts)
  rownames(mappedCounts) = mappedCounts$gene_id
  mappedCounts$gene_id<-NULL
  
  seTemp = SummarizedExperiment(assays = list(counts = as.matrix(mappedCounts)), colData = sampleData, rowData = geneData)
  
  sePath = paste0(procDataDir, '/', run, '/mappedCounts')
  
  if(!dir.exists(paste0(procDataDir, '/', run))){
    dir.create(paste0(procDataDir, '/', run))
  }
  if(!dir.exists(sePath)){
    dir.create(sePath)
  }
  
  # saveHDF5SummarizedExperiment(seTemp, dir = sePath, verbose = T, replace = T) # deprecated since summarizedExperiment update
  saveRDS(seTemp, paste0(sePath, '/se_mappedCounts.rds'))
  
  cat("Run-specific mappedCounts summarizedExperiment saved. Filtering out genes that account for a median of >1% of all reads over all samples and saving as mappedCountsForTPM summarized experiment.\n")
  
  badGenes = rawCounts %>%
    filter(gene_id != "__alignment_not_unique" & gene_id != "__ambiguous" & gene_id != "__no_feature" & gene_id != "__not_aligned" & gene_id != "__too_low_aQual") %>%
    filter(tagID != "Unassigned") %>%
    unite(combID, experiment:tagID, sep = "_", remove = T) %>%
    group_by(combID) %>%
    mutate(totalCounts = sum(counts), fracCounts = counts/totalCounts) %>%
    group_by(gene_id) %>%
    summarise(medianFrac = median(fracCounts), meanFrac = mean(fracCounts), maxFrac = max(fracCounts)) %>%
    filter(medianFrac > 0.01) %>%
    left_join(geneIDtoGeneSymbol, by = "gene_id") %>%
    arrange(-medianFrac)
  
  cat(paste0("Found ", as.character(nrow(badGenes)), " genes with counts accounting for a median of >1% of reads over all samples. Filtering them out and saving their info in: ", run, "/highReadCounts.txt\n"))
  
  write.table(badGenes, file = paste0(procDataDir, '/', run, '/highReadCounts.txt'), row.names = F, quote = F)
  
  mappedCountsForTPM = rawCounts %>%
    filter(gene_id != "__alignment_not_unique" & gene_id != "__ambiguous" & gene_id != "__no_feature" & gene_id != "__not_aligned" & gene_id != "__too_low_aQual") %>%
    filter(tagID != "Unassigned") %>%
    filter(!(gene_id %in% badGenes$gene_id)) %>%
    unite(combID, experiment:tagID, sep = "_", remove = T) %>%
    spread(combID, counts)
  mappedCountsForTPM<-as.data.frame(mappedCountsForTPM)
  rownames(mappedCountsForTPM) = mappedCountsForTPM$gene_id
  mappedCountsForTPM$gene_id<-NULL
  
  geneDataForTPM  = geneData[rownames(mappedCountsForTPM),]
  seTemp = SummarizedExperiment(assays = list(counts = as.matrix(mappedCountsForTPM)), colData = sampleData, rowData = geneDataForTPM)
  
  sePath = paste0(procDataDir, '/', run, '/mappedCountsWithoutHighFracGenes')
  # saveHDF5SummarizedExperiment(seTemp, dir = sePath, verbose = T, replace = T)
  if(!dir.exists(sePath)){
    dir.create(sePath)
  }
  saveRDS(seTemp, paste0(sePath, '/se_mappedCountsWithoutHighFracGenes.rds'))
  
  cat("Run-specific mappedCountsWithoutHighFracGenes summarizedExperiment saved. Calculating TPMs (on all mappedCounts) for all genes and all samples...\n")
  
  tpms = generateTPMfromCounts(rawCounts, geneIDlengths)
  tpmSpread = tpms %>%
    filter(tagID != "Unassigned") %>%
    unite(combID, experiment:tagID, sep = "_", remove = T) %>%
    dplyr::select(combID, gene_id, tpm) %>%
    spread(combID, tpm)
  tpmSpread <- as.data.frame(tpmSpread)
  rownames(tpmSpread) = tpmSpread$gene_id
  tpmSpread$gene_id<-NULL
  
  cat("TPMs calculated for all genes in all samples. Saving run-specific TPM summarizedExperiment object...\n")
  
  seTemp = SummarizedExperiment(assays = list(TPMs= as.matrix(tpmSpread)), colData = sampleData, rowData = geneData)
  sePath = paste0(procDataDir, '/', run, '/allTPMs')
  # saveHDF5SummarizedExperiment(seTemp, dir = sePath, verbose = T, replace = T)
  if(!dir.exists(sePath)){
    dir.create(sePath)
  }
  saveRDS(seTemp, paste0(sePath, '/se_allTPMs.rds'))
  cat(paste0("Run-specific TPM SE object saved. Done with ", run, ".\n"))
  
#   print("Run-specific all-TPM summarizedExperiment object saved. Filtering out samples with <50% mapping...")
#   unmapDatFilt = unmapDat %>%
# 	  mutate(sampleID = poolID) %>%
# 	  filter(unmappedPercentageOfTotal <= 0.5)
#   rawCountsFilt = rawCounts %>%
#   	inner_join(unmapDatFilt, by = c("experiment","sampleID", "tagID"))
#   sampleDataFilt = sampleData %>% 
#   	inner_join(unmapDatFilt, by = c("experiment","poolID", "tagID", "unmappedPercentageOfTotal")) %>%
#     mutate(combID = paste0(experiment, '_', poolID, '_', tagID))
#   sampleDataFilt = as.data.frame(sampleDataFilt)
#   rownames(sampleDataFilt) <- sampleDataFilt$combID
#   sampleDataFilt$combID <- NULL
#   
#   print("Calculating TPM and filtering out genes with median(TPM) < 0.25 over filtered samples in this run...")
#   tpms = generateTPMfromCounts(rawCountsFilt, geneIDlengths)
#   tpmSpread = tpms %>%
#     filter(tagID != "Unassigned") %>%
#     unite(combID, experiment:tagID, sep = "_", remove = T) %>%
#     dplyr::select(combID, gene_id, tpm) %>%
#     group_by(gene_id) %>%
#     mutate(medianTPM = median(tpm)) %>%
#     filter(medianTPM >= 0.25) %>%
#     dplyr::select(-medianTPM) %>%
#     spread(combID, tpm)
#   
#   filtGenes = as_tibble(as.data.frame(list(
#   	gene_id = tpmSpread$gene_id)))
#   geneDataFilt = inner_join(geneData, filtGenes, by = "gene_id")
#   geneDataFilt = as.data.frame(geneDataFilt)
#   rownames(geneDataFilt)<-geneDataFilt$gene_id
#   
#   tpmSpread <- as.data.frame(tpmSpread)
#   rownames(tpmSpread) = tpmSpread$gene_id
#   tpmSpread$gene_id<-NULL
#   print("TPMs calculated for filtered samples. Saving run-specific TPM summarizedExperiment object...")
# 
#   seTemp = SummarizedExperiment(assays = list(filteredTPMs= as.matrix(tpmSpread)), colData = sampleDataFilt, rowData = geneDataFilt)
#   sePath = paste0(procDataDir, '/', run, '/filteredTPMs')
#   saveHDF5SummarizedExperiment(seTemp, dir = sePath, verbose = T, replace = T)
#   
#   print(paste0("Run-specific TPM summarizedExperiment saved. Appending ", run, " mappedCounts and TPMs to full tables..."))
  
  
  
}


##################
# for paired-drug samples
##################
# load count data
runList = c("RNAtagTest6",
  "iCardTest6set10x"
)

countDir = "extractedData/RNAtagSeq/"

countFileName = "htseq_meltedData.tsv"

cat("Looping through the following runs:\n")
cat(runList)
counts = list()
for (run in runList){
  
  cat(paste0("Working on ", run, "...\n"))
  
  cat("Loading run-specific metadata...\n")
  # load per-run metadata
  
  wellLayoutFile = paste0("metadata/", run, "/wellIDpoolTagMapping.csv")      # columns: rowLet, colNum, poolID, TagID 
  wellContentFile = paste0("metadata/", run, "/conditionsPositions_tall.csv") # Sample data per well: drug, dose, if it's HeLa
  wellManualInspectionFile = paste0("metadata/", run, "/transImageNotes.txt") # Manual inspection of interior wells
  pairInfoFile = "metadata/RNAtagTest6/20171205_drugPairsWithPairID_withStockPositions.txt"
  # note: cell type is inferred from plate ID; iCard/GM00942 grown on separate plates.
  
  # check for custom layout, otherwise assume a default
  if (file.exists(wellLayoutFile)){
    wellLayoutTab <- as_tibble(read.csv(wellLayoutFile, stringsAsFactors = F))
  } else {
    wellLayoutTab = as_tibble(data.frame(
      rowLet = as.character(rep(LETTERS[1:8], times = 12)),
      colNum = rep(1:12, each = 8),
      poolID = paste0(run, "-pool", as.character(rep(1:3,each = 32))),
      tagID = rep(c(paste0("Tag0", as.character(1:9)), paste0("Tag", as.character(10:32))), times = 3)
    ))
  }
  #edit the following out when updated...
  wellLayoutTab2 = wellLayoutTab %>% 
    separate(poolID, into = c("run", "pID"), sep = "\\-", remove = T) %>% 
    mutate(runID = run) %>% 
    unite(poolID, c(runID, pID), sep = "-") %>% 
    dplyr::select(-run) 
  
  wellContentTab <- as_tibble(read.csv(wellContentFile, stringsAsFactors = F)) %>%
    separate(platePosition, into = c("rowLet", "colNum"), sep = "\\_") %>%
    mutate(colNum = as.integer(colNum))
  
  wellManualTab <- as_tibble(read.table(wellManualInspectionFile, sep = "\t", header = T, stringsAsFactors = F)) %>%
    separate(platePosition, into = c("rowLet", "colNum"), sep = "\\_") %>%
    mutate(colNum = as.integer(colNum))
  colnames(wellManualTab) = c("imageID", "rowLet", "colNum", paste0(colnames(wellManualTab)[4:ncol(wellManualTab)], "_", "ManualFileTitle"))
  
  wellPositions = inner_join(wellLayoutTab2, wellContentTab, by = c("rowLet", "colNum")) %>%
    mutate(experiment = run) %>%
    unite(combID, c(experiment, poolID, tagID), sep = "_", remove = T) %>%
    full_join(wellManualTab, by = c("rowLet", "colNum"))
  
  cat("Metadata loaded. Loading count data...\n")
  # load run count data
  tempFile = paste0(countDir, run, '/', countFileName)
  rawCounts = as_tibble(read.table(tempFile, header = T))
  
  cat("Count data loaded. Checking mapping stats...\n")
  # check mapping stats
  unmappingResults = unmappedRNAtagReadsPerSample(rawCounts)
  unmapDat = unmappingResults$unmappingSum %>%
    filter(readStatus == "unmappedPerc") %>%
    mutate(unmappedPercentageOfTotal = percentage) %>%
    dplyr::select(-c(readStatus, percentage))
  
  cat("Mapping stats done. Checking counts per sample...\n")
  # check reads per sample
  countsPerSample = rawCounts %>%
    group_by(experiment, sampleID, tagID) %>%
    summarise(totalCounts = sum(counts)) %>% ungroup() %>%
    unite(combID, experiment:tagID, sep = "_", remove = T)
  
  cat("Counts per sample done. Assembing sampleData (colData)...\n")
  
  pairInfo = as_tibble(read.table(pairInfoFile, sep = "\t", header = T, stringsAsFactors = F))
  
  sampleData = inner_join(countsPerSample, wellPositions, by = "combID") %>%
    mutate(foldMedianIC50 = "10x", 
           pairID = ifelse(grepl("DMSO", as.character(condition)), "DMSO", 
                           ifelse(grepl("PBS", condition), "PBS", as.character(condition))),
           cellType = ifelse(grepl("iCard", combID), "iCard", "GM00942"),
           finalDMSOconcentration.percentByVolume = 0.28,
           treatmentDuration.hours = 96) %>%
    left_join(pairInfo, by = "pairID") %>%
    separate(combID, into = c("experiment","poolID", "tagID"), sep = "\\_", remove = F) %>%
    inner_join(unmapDat, by = c("experiment","poolID", "tagID"))
  
  sampleData1 = sampleData %>%
    dplyr::select(combID, drug1) %>%
    mutate(Product.Name = drug1) %>%
    left_join(drugIDtable, by = "Product.Name") %>%
    mutate(finalDrugConcentration.nM = ifelse(is.na(medianIC50..nM.), 0, medianIC50..nM. * 100))
  colnames(sampleData1) = c("combID", "drug1", paste0(colnames(sampleData1)[3:ncol(sampleData1)], "_drug1"))
  
  sampleData2 = sampleData %>%
    dplyr::select(combID, drug2) %>%
    mutate(Product.Name = drug2) %>%
    left_join(drugIDtable, by = "Product.Name") %>%
    mutate(finalDrugConcentration.nM = ifelse(is.na(medianIC50..nM.), 0, medianIC50..nM. * 100))
  colnames(sampleData2) = c("combID", "drug2", paste0(colnames(sampleData2)[3:ncol(sampleData2)], "_drug2"))
  
  sampleData = left_join(sampleData, sampleData1, by = c("combID", "drug1")) %>%
    left_join(sampleData2, by = c("combID", "drug2"))
 
  sampleData = as.data.frame(sampleData)
  rownames(sampleData) = sampleData$combID
  sampleData$combID <- NULL
  
  cat("sampleData assembled. Saving run-specific mappedCounts summarizedExperiment object...\n")
  mappedCounts = rawCounts %>%
    filter(gene_id != "__alignment_not_unique" & gene_id != "__ambiguous" & gene_id != "__no_feature" & gene_id != "__not_aligned" & gene_id != "__too_low_aQual") %>%
    filter(tagID != "Unassigned") %>%
    unite(combID, experiment:tagID, sep = "_", remove = T) %>%
    spread(combID, counts)
  mappedCounts<-as.data.frame(mappedCounts)
  rownames(mappedCounts) = mappedCounts$gene_id
  mappedCounts$gene_id<-NULL
  
  seTemp = SummarizedExperiment(assays = list(counts = as.matrix(mappedCounts)), colData = sampleData, rowData = geneData)
  if(!dir.exists(paste0(procDataDir, '/', run))){
    dir.create(paste0(procDataDir, '/', run))
  }
  sePath = paste0(procDataDir, '/', run, '/mappedCounts')
  if(!dir.exists(sePath)){
    dir.create(sePath)
  }
  # saveHDF5SummarizedExperiment(seTemp, dir = sePath, verbose = T, replace = T)
  saveRDS(seTemp, paste0(sePath, '/se_mappedCounts.rds'))
  
  cat("Run-specific mappedCounts summarizedExperiment saved. Filtering out genes that account for a median of >1% of all reads over all samples and saving as mappedCountsForTPM summarized experiment.\n")
  
  badGenes = rawCounts %>%
    filter(gene_id != "__alignment_not_unique" & gene_id != "__ambiguous" & gene_id != "__no_feature" & gene_id != "__not_aligned" & gene_id != "__too_low_aQual") %>%
    filter(tagID != "Unassigned") %>%
    unite(combID, experiment:tagID, sep = "_", remove = T) %>%
    group_by(combID) %>%
    mutate(totalCounts = sum(counts), fracCounts = counts/totalCounts) %>%
    group_by(gene_id) %>%
    summarise(medianFrac = median(fracCounts), meanFrac = mean(fracCounts), maxFrac = max(fracCounts)) %>%
    filter(medianFrac > 0.01) %>%
    left_join(geneIDtoGeneSymbol, by = "gene_id") %>%
    arrange(-medianFrac)
  
  cat(paste0("Found ", as.character(nrow(badGenes)), " genes with counts accounting for a median of >1% of reads over all samples. Filtering them out and saving their info in: ", run, "/highReadCounts.txt\n"))
  
  write.table(badGenes, file = paste0(procDataDir, '/', run, '/highReadCounts.txt'), row.names = F, quote = F)
  
  mappedCountsForTPM = rawCounts %>%
    filter(gene_id != "__alignment_not_unique" & gene_id != "__ambiguous" & gene_id != "__no_feature" & gene_id != "__not_aligned" & gene_id != "__too_low_aQual") %>%
    filter(tagID != "Unassigned") %>%
    filter(!(gene_id %in% badGenes$gene_id)) %>%
    unite(combID, experiment:tagID, sep = "_", remove = T) %>%
    spread(combID, counts)
  mappedCountsForTPM<-as.data.frame(mappedCountsForTPM)
  rownames(mappedCountsForTPM) = mappedCountsForTPM$gene_id
  mappedCountsForTPM$gene_id<-NULL
  
  geneDataForTPM  = geneData[rownames(mappedCountsForTPM),]
  seTemp = SummarizedExperiment(assays = list(counts = as.matrix(mappedCountsForTPM)), colData = sampleData, rowData = geneDataForTPM)
  
  sePath = paste0(procDataDir, '/', run, '/mappedCountsWithoutHighFracGenes')
  if(!dir.exists(sePath)){
    dir.create(sePath)
  }
  # saveHDF5SummarizedExperiment(seTemp, dir = sePath, verbose = T, replace = T)
  saveRDS(seTemp, paste0(sePath, '/se_mappedCountsWithoutHighFracGenes.rds'))
  cat("Run-specific mappedCountsWithoutHighFracGenes summarizedExperiment saved. Calculating TPMs for all genes and all samples...\n")
  
  tpms = generateTPMfromCounts(rawCounts, geneIDlengths)
  tpmSpread = tpms %>%
    filter(tagID != "Unassigned") %>%
    unite(combID, experiment:tagID, sep = "_", remove = T) %>%
    dplyr::select(combID, gene_id, tpm) %>%
    spread(combID, tpm)
  tpmSpread <- as.data.frame(tpmSpread)
  rownames(tpmSpread) = tpmSpread$gene_id
  tpmSpread$gene_id<-NULL
  
  cat("TPMs calculated for all genes in all samples. Saving run-specific TPM summarizedExperiment object...\n")
  
  seTemp = SummarizedExperiment(assays = list(TPMs= as.matrix(tpmSpread)), colData = sampleData, rowData = geneData)
  sePath = paste0(procDataDir, '/', run, '/allTPMs')
  if(!dir.exists(sePath)){
    dir.create(sePath)
  }
  # saveHDF5SummarizedExperiment(seTemp, dir = sePath, verbose = T, replace = T)
  saveRDS(seTemp, paste0(sePath, '/se_allTPMs.rds'))
  
  #   print("Run-specific all-TPM summarizedExperiment object saved. Filtering out samples with <50% mapping...")
  #   unmapDatFilt = unmapDat %>%
  # 	  mutate(sampleID = poolID) %>%
  # 	  filter(unmappedPercentageOfTotal <= 0.5)
  #   rawCountsFilt = rawCounts %>%
  #   	inner_join(unmapDatFilt, by = c("experiment","sampleID", "tagID"))
  #   sampleDataFilt = sampleData %>% 
  #   	inner_join(unmapDatFilt, by = c("experiment","poolID", "tagID", "unmappedPercentageOfTotal")) %>%
  #     mutate(combID = paste0(experiment, '_', poolID, '_', tagID))
  #   sampleDataFilt = as.data.frame(sampleDataFilt)
  #   rownames(sampleDataFilt) <- sampleDataFilt$combID
  #   sampleDataFilt$combID <- NULL
  #   
  #   print("Calculating TPM and filtering out genes with median(TPM) < 0.25 over filtered samples in this run...")
  #   tpms = generateTPMfromCounts(rawCountsFilt, geneIDlengths)
  #   tpmSpread = tpms %>%
  #     filter(tagID != "Unassigned") %>%
  #     unite(combID, experiment:tagID, sep = "_", remove = T) %>%
  #     dplyr::select(combID, gene_id, tpm) %>%
  #     group_by(gene_id) %>%
  #     mutate(medianTPM = median(tpm)) %>%
  #     filter(medianTPM >= 0.25) %>%
  #     dplyr::select(-medianTPM) %>%
  #     spread(combID, tpm)
  #   
  #   filtGenes = as_tibble(as.data.frame(list(
  #   	gene_id = tpmSpread$gene_id)))
  #   geneDataFilt = inner_join(geneData, filtGenes, by = "gene_id")
  #   geneDataFilt = as.data.frame(geneDataFilt)
  #   rownames(geneDataFilt)<-geneDataFilt$gene_id
  #   
  #   tpmSpread <- as.data.frame(tpmSpread)
  #   rownames(tpmSpread) = tpmSpread$gene_id
  #   tpmSpread$gene_id<-NULL
  #   print("TPMs calculated for filtered samples. Saving run-specific TPM summarizedExperiment object...")
  # 
  #   seTemp = SummarizedExperiment(assays = list(filteredTPMs= as.matrix(tpmSpread)), colData = sampleDataFilt, rowData = geneDataFilt)
  #   sePath = paste0(procDataDir, '/', run, '/filteredTPMs')
  #   saveHDF5SummarizedExperiment(seTemp, dir = sePath, verbose = T, replace = T)
  #   
  #   print(paste0("Run-specific TPM summarizedExperiment saved. Appending ", run, " mappedCounts and TPMs to full tables..."))
  
  
  
}


########### 
# Store one SE for all runs considered
###########

runList = c("RNAtagTest4",
            "RNAtagTest5",
            "RNAtagTest6",
            "iCardTest4set100x",
            "iCardTest5set100x",
            "iCardTest6set10x")

amcall = list()
cdmcall = list()
rdmcall = list()
atall = list()
cdtall = list()
rdtall = list()
for (run in runList){
  
  cat(paste0("Working on ", run, "...\n"))
  
  sePath = paste0(procDataDir, '/', run, '/mappedCounts/se_mappedCounts.rds')
  # mcse = loadHDF5SummarizedExperiment(sePath)
  mcse = readRDS(sePath)
  amc = assay(mcse)
  cdmc = colData(mcse)
  rdmc = rowData(mcse)
  
  sePath = paste0(procDataDir, '/', run, '/allTPMs/se_allTPMs.rds')
  # atse = loadHDF5SummarizedExperiment(sePath)
  atse = readRDS(sePath)
  aat = assay(atse)
  cdt = colData(atse)
  rdt = rowData(atse)
  
  if(is.null(dim(amcall))){
    amcall = amc
    
    cdmcT = cdmc
    cdmcT$combID = rownames(cdmcT)
    rownames(cdmcT) = NULL
    cdmcT = as_tibble(cdmcT)
    cdmcall = cdmcT
    
    rdmcall = rdmc
    
    atall = aat
    
    cdtT = cdt
    cdtT$combID = rownames(cdtT)
    rownames(cdtT) = NULL
    cdtT = as_tibble(cdtT)
    cdtall = cdtT
    
    rdtall = rdt
  } else {
    amcall = cbind(amcall, amc)
    
    cdmcT = cdmc
    cdmcT$combID = rownames(cdmcT)
    rownames(cdmcT) = NULL
    cdmcT = as_tibble(cdmcT)
    cdmcall = bind_rows(cdmcall, cdmcT)
    
    if(sum(rdmc$gene_id != rdmcall$gene_id) > 0){
      cat(paste0("Warning: ", run, " has different mappedCounts genedata.\n"))
    }
    
    atall = cbind(atall, aat)
    
    cdtT = cdt
    cdtT$combID = rownames(cdtT)
    rownames(cdtT) = NULL
    cdtT = as_tibble(cdtT)
    cdtall = bind_rows(cdtall, cdtT)
    
    if(sum(rdt$gene_id != rdtall$gene_id) > 0){
      cat(paste0("Warning: ", run, " has different allTPMs genedata.\n"))
    }
  }
  
  cat(paste0("Finished binding ", run, "\n"))
  
}

cdmcall = as.data.frame(cdmcall)
rownames(cdmcall) = cdmcall$combID

cdtall = as.data.frame(cdtall)
rownames(cdtall) = cdtall$combID

if(!dir.exists(paste0(procDataDir, '/allExperiments'))){
  dir.create(paste0(procDataDir, '/allExperiments'))
}

rownames(rdmcall) = rdmcall$gene_id
seTemp = SummarizedExperiment(assays = list(counts = as.matrix(amcall)), colData = cdmcall, rowData = rdmcall)
sePath = paste0(procDataDir, '/allExperiments/mappedCounts')
if(!dir.exists(sePath)){
  dir.create(sePath)
}
# saveHDF5SummarizedExperiment(seTemp, dir = sePath, verbose = T, replace = T)
saveRDS(seTemp, paste0(sePath, '/se_mappedCounts.rds'))

rownames(rdtall) = rdtall$gene_id
seTemp = SummarizedExperiment(assays = list(counts = as.matrix(atall)), colData = cdtall, rowData = rdtall)
sePath = paste0(procDataDir, '/allExperiments/allTPMs')
if(!dir.exists(sePath)){
  dir.create(sePath)
}
# saveHDF5SummarizedExperiment(seTemp, dir = sePath, verbose = T, replace = T)
saveRDS(seTemp, paste0(sePath, '/se_allTPMs.rds'))


