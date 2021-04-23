# P3_public_analysis

This repo contains scripts needed to produce all graphs and images contained
in figures for the paper:

Mellis IA, et al. Perturbation panel profiling identifies transcription factors 
that enhance directed changes of cell identity. bioRxiv (2020); Revised version resubmitted to Cell Systems in April 2021.

For any questions, please contact IAM (ian.mellis at gmail) and the corresponding
authors for this paper, AR (arjunrajlab at gmail) and RJ (jainr at pennmedicine.upenn.edu)

In order to reproduce all graphs and images in this paper:

0. Download and install dependencies for this code:
- MATLAB_2017a or MATLAB_2019a
- rajlabimagetools (https://github.com/arjunrajlaboratory/rajlabimagetools changeset a2c6ac5)
- R v3.6.1
- R packages e1071_1.7-3, gridExtra_2.3, DESeq2_1.24.0, SummarizedExperiment_1.14.1, DelayedArray_0.10.0, BiocParallel_1.18.1, matrixStats_0.57.0, GenomicRanges_1.36.1, GenomeInfoDb_1.20.0, metaRNASeq_1.0.3, org.Hs.eg.db_3.8.2, AnnotationDbi_1.46.1, IRanges_2.18.3, S4Vectors_0.22.1, Biobase_2.44.0, BiocGenerics_0.30.0, clusterProfiler_3.12.0, readxl_1.3.1, ggrepel_0.8.2, magrittr_1.5, forcats_0.5.0, stringr_1.4.0, dplyr_1.0.2, purrr_0.3.4, readr_1.4.0, tidyr_1.1.2, tibble_3.0.3, ggplot2_3.3.2, tidyverse_1.3.0, yaml_2.2.1, and their associated dependencies
- See analysisScripts/README.txt for more details regarding R packages and script usage

1. Download all extractedData (6.3GB of intermediate-processed sequencing, network, and reprogramming experiment data) files in the structure provided for this repository from Dropbox: https://www.dropbox.com/sh/vr99twboxppn3c1/AADuO_a1O_afqZBUxl2v2LrXa?dl=0.
- your top-level project directory should have subdirectories: analysisScripts, annotations, extractedData, imagesForPaper, metadata, miscellaneous, rawData (step 2). Upon completion of step 7 of this pipeline, you will also have subdirectories processedImages, procDataScripted, graphsScripted, and logsForScripts, with names you can specify per the readme in analysisScripts.

2. Download all raw imaging data (2.2 TB) in the file structure provided from:
https://www.dropbox.com/sh/2ny7k6c4zy6zdsh/AADLNoom3YOw0Ps0B8w509R7a?dl=0

3. Edit imagesForPaper/assemblyScripts/assembleAllImages.m variables to reflect your local file paths for this project repository, raw imaging data, and your desired project
subfolder in which processed images and summary data will be stored. The default is ‘processedImages/‘, but if you would like you can specify a different folder to recreate all images and extracted data tables to compare. 

4. Edit analysisScripts/analyzeAll.sh variables to reflect your local file
paths for this project repository, processed imaging data (from 3), and your 
desired project subfolders in which summary statistics, tables, and graphs will 
be stored. The default variable names are: 
projectDir='~/Dropbox\ \(RajLab\)/Shared_IanM/cellid_201807_onward/'
analysisSubdir='analysisScripts/'
metadataSubdir='metadata/'
procDataSubdir='procDataScripted_test421/‘
procImgSubdir='processedImages/'
graphSubdir='graphsScripted_test421/‘
logSubdir='logsForScriptedAnalysis_test421/‘
runDEseq='FALSE' # switch to TRUE if you would like to rerun all DESeq2 steps - takes several days on a modern desktop due to hundreds of samples in 120+ groups, hundreds of comparisons, and downsampling analysis.

If you would like you can specify different names for procDataSubdir, procImgSubdir, and graphSubdir to generate a duplicate set of graphs and tables.

5. Open MATLAB, add rajlabimagetools with subfolders to your MATLAB path (remove hidden files), and navigate to this project repo folder. Add extractionScripts_images/ to your MATLAB path, navigate to extractionScripts_images/, and run:
>> assembleAllImages

6. Open Terminal, navigate to this project repo folder, and run:
$ bash analysisScripts/analyzeAll.sh
NOTE: this last step will require at least 12 GB of RAM to run in a reasonable
amount of time, and will take several hours on a recent model Macbook Pro or
iMac desktop.

7. See graph_figure_map_CellSystemsResub_20210423.txt for a list of which file corresponds to which figure panel