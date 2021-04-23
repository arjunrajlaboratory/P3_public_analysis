#!/usr/bin/env Rscript

# exploratoryDataAnalysis.R performs analyses toward visualizing the sources of greatest variability in the data.
# Using TPM data (min-expression filtered genes):
# Sample-to-sample variability:
#   - PCA on all samples per experiment
#   - PCA on quality-filtered samples per experiment
#   - PCA on all quality-filtered samples in all experiments
#   - hierarchical clustering of all quality-filtered samples over all genes passing QC
#   - hierarchical clustering of all quality-filtered samples over all transcription factors
# Gene-to-gene variability
#   - Within controls: CV x mean expression plots


# don't run lines 16-26 if running from R workspace
args = commandArgs(trailingOnly = T)

# check if there are exactly 3 arguments. If not, return an error
if (length(args) < 3 | length(args) > 3) {
  stop("Exactly three arguments must be supplied: projectDir, procDataSubdir and graphSubdir.", call.=FALSE)
} 
if (length(args)==2){
  projectDir = args[1]
  procDataSubdir = args[2]
  graphSubdir = args[3]
}

# if running manually in Rstudio, start here and use:
projectDir = '~/Dropbox (RajLab)/Projects/cellid/'
procDataSubdir = 'procDataScripted'
graphSubdir = 'graphsScripted'

procDataDir = paste0(projectDir, procDataSubdir)
graphDir = paste0(projectDir, graphSubdir)

setwd(projectDir)
print("Starting exploratoryDataAnalysis.R...")
print(paste0("Working in ", getwd()))


# load relevant libraries and source custom functions
# comment out if running analyzeAll.sh instead of just this script (analyzeAll.sh loads all libraries)
library(SummarizedExperiment)
library(tidyverse)
library(RColorBrewer)
source(paste0(projectDir, "analysisScripts/RNAtagUtilities.R"))

tpmPseudo = 0.00001

# for one run
tpmData = loadHDF5SummarizedExperiment(dir="procDataScripted/RNAtagTest3/filteredTPMs")

# for all runs
tpmData = loadHDF5SummarizedExperiment(dir="procDataScripted/RNAtagSeq/allRuns/filteredTPMs")

sampleData = colData(tpmData)


tpmFiltPseudoSpreadMat = as.matrix(assays(tpmData)[[1]][,sampleData$totalCounts>=1000000]) + tpmPseudo
sampleDataFilt = sampleData[sampleData$totalCounts>=1000000,]


## Hierarchical clustering
dst = dist(t(tpmFiltPseudoSpreadMat));
tree = hclust(dst, method = "average")
plot(tree)

tpmRename = t(tpmFiltPseudoSpreadMat)
rownames(tpmRename) = paste0(sampleDataFilt$Product.Name, "_", sampleDataFilt$concentration)
dst2 = dist(tpmRename)
tree2 = hclust(dst2, method = "average")
plot(tree2)

# cluster the controls
tpmFiltPseudoSpreadMatControls = t(as.matrix(assays(tpmData)[[1]][,sampleData$Product.Name == "control" & sampleData$totalCounts>=1000000]) + tpmPseudo)
sampleDataFiltControls = sampleData[sampleData$Product.Name == "control" & sampleData$totalCounts>=1000000,]
rownames(tpmFiltPseudoSpreadMatControls) = paste0("Beat-", sampleDataFiltControls$anyBeating_ManualFileTitle, "_", 
                                                  "Bleb-", sampleDataFiltControls$anyBlebbing_ManualFileTitle, "_",
                                                  sampleDataFiltControls$rowLet, sampleDataFiltControls$colNum)
dst3 = dist(tpmFiltPseudoSpreadMatControls)
tree3 = hclust(dst3, method = "average")
plot(tree3, main = "Hierarchical Clustering of iCard control wells, 1MM reads minimum")


# PCA
tpmPCA<-prcomp(t(tpmFiltPseudoSpreadMat), scale = T)
tpmPCAforPlot<-as.data.frame(tpmPCA$x)
tpmPCAforPlot$combID <- rownames(tpmPCA$x)
rownames(tpmPCAforPlot) <- NULL
tpmPCAforPlot<-as_tibble(tpmPCAforPlot) %>%
  separate(combID, into = c("experiment","poolID", "tagID"), sep = "\\_")
tpmPCAsum<-summary(tpmPCA)

tpmPCAforPlotAnnotated = tpmPCAforPlot %>%
  inner_join(as_tibble(sampleData), by = c("experiment","poolID", "tagID"))


tpmPCAplot12 <- ggplot(tpmPCAforPlotAnnotated, aes(PC1, PC2)) +
  geom_point(aes(shape = as.factor(paste0(anyBeating_ManualFileTitle, "_")),
                 color = as.factor(paste0(Product.Name, "_")))) +
  xlab(paste0("PC1: ", as.character(tpmPCAsum$importance[2,1]), " variance explained")) +
  ylab(paste0("PC2: ", as.character(tpmPCAsum$importance[2,2]), " variance explained")) +
  ggtitle("PCA on scaled TPM matrix: first two PCs\nonly including genes with medianTPM >= 0.25")
plot(tpmPCAplot12)

# controls PCA
tpmPCA<-prcomp(tpmFiltPseudoSpreadMatControls, scale = T)
tpmPCAforPlot<-as.data.frame(tpmPCA$x)
tpmPCAforPlot$combID <- rownames(tpmPCA$x)
rownames(tpmPCAforPlot) <- NULL
tpmPCAforPlot<-as_tibble(tpmPCAforPlot) %>%
  separate(combID, into = c("experiment","poolID", "tagID"), sep = "\\_")
tpmPCAsum<-summary(tpmPCA)

# tpmPCAforPlotAnnotated = tpmPCAforPlot %>%
  # inner_join(sampleDataFiltControls, by = c("experiment","poolID", "tagID"))


tpmPCAplot12 <- ggplot(tpmPCAforPlot, aes(PC1, PC2)) +
  geom_point(aes(shape = as.factor(paste0(experiment, "_")),
                 color = as.factor(paste0(poolID, "_")))) +
  xlab(paste0("PC1: ", as.character(tpmPCAsum$importance[2,1]), " variance explained")) +
  ylab(paste0("PC2: ", as.character(tpmPCAsum$importance[2,2]), " variance explained")) +
  ggtitle("PCA on scaled TPM matrix of Controls: first two PCs\nonly including genes with medianTPM >= 0.25")
plot(tpmPCAplot12)
