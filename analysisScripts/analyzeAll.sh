#!/bin/bash

# analyzeAll.sh performs all pre-figure analysis steps on RNAtag-seq count data. Run from 
# the command line in your top-level cellid project directory. Be sure to edit the following
# variables below this paragraph:
#	projectDir = top-level directory for cellid
# 	analysisSubdir = subdirectory basename of projectDir, in which you ALREADY HAVE analysis 
#		scripts run within the following analyzeAll.sh
# 	metadataSubdir = subdirectory basename of projectDir, in which you ALREADY HAVE files 
# 		containing the metadata needed for analyzeAll.sh Rscripts
# 	procDataSubdir = subdirectory basename of projectDir, in which processed count data will be saved
# 	graphSubdir = subdirectory basename of projectDir, in which QC graphs will be saved
#	logSubdir = subdirectory basename of projectDir, in which logs from this scripts will
#		be saved. One log file, analyzeAll/$DATE_$TIME.log per execution of this script.
#
# Example usage:
# 	bash analysisScripts/analyzeAll.sh

projectDir='~/Dropbox\ \(RajLab\)/Shared_IanM/cellid_201807_onward/'
analysisSubdir='analysisScripts/'
metadataSubdir='metadata/'
procDataSubdir='procDataScripted_test421/'
procImgSubdir='processedImages/'
graphSubdir='graphsScripted_test421/'
logSubdir='logsForScriptedAnalysis_test421/'
runDEseq='FALSE' # switch to TRUE if you would like to rerun all DESeq2 steps - takes several days on a modern desktop due to hundreds of samples in 120+ groups, hundreds of comparisons, and downsampling analysis.

# make relevant subdirectories if they don't exist
if [ ! -d $procDataSubdir ]; then
	echo "Making $procDataSubdir."
	mkdir $procDataSubdir
fi

if [ ! -d $graphSubdir ]; then
	echo "Making $graphSubdir."
	mkdir $graphSubdir
fi

if [ ! -d $logSubdir ]; then
	echo "Making $logSubdir."
	mkdir $logSubdir
	echo "Making analyzeAll subdirectory in $logSubdir"
	mkdir $logSubdir/analyzeAll
fi

# make the journal
JOURNAL=$logSubdir/analyzeAll/$(date +%Y-%m-%d_%H-%M).log
echo "Starting analyzeAll.sh..." >> $JOURNAL
date >> $JOURNAL



## Load counts, do QC, and transform to TPM
# cleanCountsAndCalcTPM.R loads all raw count data, quality-filters it, and transforms to TPM.
# 	quality filters:
# 		- no sample with less than 60% mapping
# 		- no gene with median TPM per plate <= 0.25
# 	generates: 
#		- allMappedCounts: summarizedExperiment object in procDataDir/allMappedCounts/
#				genes x samples table of HTseq-mapped counts
#				includes all genes
#				includes all samples
#		- allQcFilteredTPMs: summarizedExperiment object in procDataDir/allQcFilteredTPMs/
#				genes x samples table of TPM values of QC-filtered genes and samples

ccmd="Rscript $analysisSubdir/cleanCountsAndCalcTPM.R $projectDir $procDataSubdir $graphSubdir"
echo "Starting..." >> $JOURNAL
date >> $JOURNAL
echo "$ccmd" >> $JOURNAL
eval "$ccmd"
date >> $JOURNAL
echo "Finished." # working 4/22/21

# ## Move previously run DESeq2 results to graphSubdir if runDEseq=F
# if [ $runDEseq=="FALSE" ]; then
	mvcmd="bash $analysisSubdir/move_DESeq2_results.sh $projectDir $graphSubdir $JOURNAL"
	echo "runDEseq is FALSE, transferring already-run DESeq2 results. Starting..." >> $JOURNAL
	date >> $JOURNAL
	echo "$mvcmd" >> $JOURNAL
	eval "$mvcmd" # working 4/22/21
	date >> $JOURNAL
	echo "Finished."
# fi

# Check simple differential expression between each perturbation/drug and controls
# and calculate perturbability statistics.
# full_DESeq2_perturbability.R uses DESeq2 to check for differential expression of all genes in
# all samples passing sample-level QC. It also calculates perturbability statistics on all genes
# in all samples passing sample-level QC.
dpmd="Rscript $analysisSubdir/full_DESeq2_perturbability.R $projectDir $procDataSubdir $graphSubdir $runDEseq"
echo "Starting..." >> $JOURNAL
date >> $JOURNAL
echo "$dpmd" >> $JOURNAL
eval "$dpmd"
date >> $JOURNAL
echo "Finished." # seems to be working 4/22/21. Decimal rounding errors introduced (3-4 places out), but they seem to be floating point issues.

## Graph perturbability statistics
# graphing_perturbability.R plots perturbability statistics in exploratory analysis and to
# compare with other gene-level features.
gpmd="Rscript $analysisSubdir/graphing_perturbability.R $projectDir $procDataSubdir $graphSubdir"
echo "Starting..." >> $JOURNAL
date >> $JOURNAL
echo "$gpmd" >> $JOURNAL
eval "$gpmd"
date >> $JOURNAL
echo "Finished." # working 4/22/21

## Caclulate tissue-specificity scores for all genes in filtered GTEx data
# GTEx_calculateTissueSpecificity.R loads all GTEx TPM data and calculates JSD-based
# tissue-specificity of mean gene expression scores, as explained in 
# http://genesdev.cshlp.org/content/25/18/1915.full#sec-16
# 	generates: 
#		- SMTS_JSsp_TissueSpecificity*.txt: tab-delimited tables of JSsp values for each
#				gene in each SMTS-designated tissue type.
gpcmd="Rscript $analysisSubdir/GTEx_cleanCountsAndCalcTPM.R $projectDir $procDataSubdir $graphSubdir"
gtcmd="Rscript $analysisSubdir/CellNetGTExVsMyData.R $projectDir $procDataSubdir $graphSubdir"
gscmd="Rscript $analysisSubdir/GTEx_calculateTissueSpecificity.R $projectDir $procDataSubdir $graphSubdir"
ggcmd="Rscript $analysisSubdir/GTEx_fibsAndiCards_all6_comparePerturbabilityAndSpecificity.R $projectDir $procDataSubdir $graphSubdir"
echo "Starting..." >> $JOURNAL
date >> $JOURNAL
echo "$gpcmd" >> $JOURNAL
eval "$gpcmd" # seems to be working 4/22/21
date >> $JOURNAL
echo "$gtcmd" >> $JOURNAL
eval "$gtcmd" # seems to be working 4/22/21
date >> $JOURNAL
echo "$gscmd" >> $JOURNAL
eval "$gscmd" # seems to be working 4/22/21
date >> $JOURNAL
echo "$ggcmd" >> $JOURNAL
eval "$ggcmd" # seems to be working 4/22/21
date >> $JOURNAL
echo "Finished."

## Perform GSEA on perturbability-based gene sets using GO terms
# GO_analysis.R loads perturbability data and performs gene set overenrichment analysis
# 	generates:
# 		- 
gocmd="Rscript $analysisSubdir/GO_analysis.R $projectDir $procDataSubdir $graphSubdir"
echo "Starting..." >> $JOURNAL
date >> $JOURNAL
echo "$gocmd" >> $JOURNAL
eval "$gocmd" # seems to be working
date >> $JOURNAL
echo "Finished."

## Plot heatmaps from P3 data
# plotting_heatmaps.R generates heatmaps for exploratory analysis of RNAtag-seq results from P3.
hmcmd="Rscript $analysisSubdir/plotting_heatmaps.R $projectDir $procDataSubdir $graphSubdir"
echo "Starting..." >> $JOURNAL
date >> $JOURNAL
echo "$hmcmd" >> $JOURNAL
eval "$hmcmd" # seems to be working
date >> $JOURNAL
echo "Finished."

## Plot Borkent mouse reprogramming screen
# borkent_reanalysis.R summarizes the MEF reprogramming screen in Borkent et al., 2016
# and compares the results to orthologous gene responsiveness in this study. Generates plots.
brcmd="Rscript $analysisSubdir/borkent_reanalysis.R $projectDir $procDataSubdir $graphSubdir"
echo "Starting..." >> $JOURNAL
date >> $JOURNAL
echo "$brcmd" >> $JOURNAL
eval "$brcmd" # seems to be working
date >> $JOURNAL
echo "Finished."

## Plot regulon responsiveness results
# graphing_regulons.R calculates regulon-level responsiveness scores for each transcription factor
# by taking a weighted average over all its targets of nDESeq2conditionsUp*fraction of conditions
# in which up-regulated*edge weight. Then plots against gene-level responsiveness
regmd="Rscript $analysisSubdir/graphing_regulons.R $projectDir $procDataSubdir $graphSubdir"
echo "Starting..." >> $JOURNAL
date >> $JOURNAL
echo "$regmd" >> $JOURNAL
eval "$regmd" # seems to be working
date >> $JOURNAL
echo "Finished."

## Plot Differential expression meta-analysis
# DE_metaanalysis.R performs meta-analysis by Fisher's method on DESeq2-generated adj. p-values,
# and then plots the results
fshmd="Rscript $analysisSubdir/DE_metaanalysis.R $projectDir $procDataSubdir $graphSubdir"
echo "Starting..." >> $JOURNAL
date >> $JOURNAL
echo "$fshmd" >> $JOURNAL
eval "$fshmd" # seems to be working
date >> $JOURNAL
echo "Finished."

# ## Move previously run DESeq2 iPSC-CM knockdown results to graphSubdir if runDEseq=F
# if [ $runDEseq=="FALSE" ]; then
	mvcmd="bash $analysisSubdir/move_DESeq2_results_iPSCMKDs.sh $projectDir $graphSubdir $JOURNAL"
	echo "runDEseq is FALSE, transferring already-run DESeq2 results. Starting..." >> $JOURNAL
	date >> $JOURNAL
	echo "$mvcmd" >> $JOURNAL
	eval "$mvcmd" # working 4/22/21
	date >> $JOURNAL
	echo "Finished."
# fi

## Plot results of iPSC-CM knockdown experiments
# iPS-CM-KDs_RNAseq_analysis.R performed differential expression analysis with DESeq2 (or skips it
# if previously done), performed GO analysis on the results, and plots heatmaps of p-values for each
# enriched term
ikdmd="Rscript $analysisSubdir/iPS-CM-KDs_RNAseq_analysis.R $projectDir $procDataSubdir $graphSubdir $runDEseq"
echo "Starting..." >> $JOURNAL
date >> $JOURNAL
echo "$ikdmd" >> $JOURNAL
eval "$ikdmd" # seems to be working
date >> $JOURNAL
echo "Finished."

## Plot results of CellNet analysis
# CellNet_setup.R plots the output of CellNet analysis of a subset of P3 samples
cnemd="Rscript $analysisSubdir/CellNet_setup.R $projectDir $procDataSubdir $graphSubdir"
echo "Starting..." >> $JOURNAL
date >> $JOURNAL
echo "$cnemd" >> $JOURNAL
eval "$cnemd" # seems to be working
date >> $JOURNAL
echo "Finished."

## Plot results from fibroblast-cardiomyocyte transdifferentiation experiments
# scan results
xd1cmd="Rscript $analysisSubdir/IAM942_34_subset.R $projectDir $procImgSubdir $graphSubdir"
xd2cmd="Rscript $analysisSubdir/IAMimmHCF_169_30_subset.R $projectDir $procImgSubdir $graphSubdir"
xd3cmd="Rscript $analysisSubdir/IAMimmHCF30.R $projectDir $procImgSubdir $graphSubdir"
echo "Starting..." >> $JOURNAL
date >> $JOURNAL
echo "$xd1cmd" >> $JOURNAL
eval "$xd1cmd" # seems to be working
date >> $JOURNAL
echo "$xd2cmd" >> $JOURNAL
eval "$xd2cmd" # seems to be working
date >> $JOURNAL
echo "$xd3cmd" >> $JOURNAL
eval "$xd3cmd" # seems to be working
date >> $JOURNAL
echo "Finished."

# Plot results from hiF-T reprogramming experiments
hi1cmd="Rscript $analysisSubdir/hiFT3a_qPCR.R $projectDir $procDataSubdir $graphSubdir"
hi2cmd="Rscript $analysisSubdir/hiFT4a.R $projectDir $procDataSubdir $graphSubdir"
hi3cmd="Rscript $analysisSubdir/hiFT3b.R $projectDir $procDataSubdir $graphSubdir"
hi4cmd="Rscript $analysisSubdir/hiFT4b.R $projectDir $procDataSubdir $graphSubdir"
hi5cmd="Rscript $analysisSubdir/hiFT3c4c.R $projectDir $procDataSubdir $graphSubdir"
echo "Starting..." >> $JOURNAL
date >> $JOURNAL
echo "$hi1cmd" >> $JOURNAL
eval "$hi1cmd" # seems to be working
date >> $JOURNAL
echo "$hi2cmd" >> $JOURNAL
eval "$hi2cmd" # seems to be working
date >> $JOURNAL
echo "$hi3cmd" >> $JOURNAL
eval "$hi3cmd" # seems to be working
date >> $JOURNAL
echo "$hi4cmd" >> $JOURNAL
eval "$hi4cmd" # seems to be working
date >> $JOURNAL
echo "$hi5cmd" >> $JOURNAL
eval "$hi5cmd" # seems to be working
date >> $JOURNAL
echo "Finished."

echo "Finished  analyzeAll.sh..." >> $JOURNAL
date >> $JOURNAL