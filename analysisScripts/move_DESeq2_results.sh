#!/bin/bash

# move_DESeq2_results.sh moves (copies) previously run move_DESeq2_results with bootstraps into the target
# graphs directory. Expects to be running in projectDir. 

projectDir=$1
graphSubdir=$2
JOURNAL=$3

if [ ! -d $graphSubdir ]; then
    mkdir $graphSubdir
fi
if [ ! -d $graphSubdir/differentialExpression ]; then
    mkdir $graphSubdir/differentialExpression
fi
if [ ! -d $graphSubdir/differentialExpression/allSamples ]; then
    mkdir $graphSubdir/differentialExpression/allSamples
fi
if [ ! -d $graphSubdir/differentialExpression/allSamples/iCards ]; then
    mkdir $graphSubdir/differentialExpression/allSamples/iCards
fi
if [ ! -d $graphSubdir/differentialExpression/allSamples/iCards/downsampled ]; then
    mkdir $graphSubdir/differentialExpression/allSamples/iCards/downsampled
fi
if [ ! -d $graphSubdir/differentialExpression/allSamples/fibroblasts ]; then
    mkdir $graphSubdir/differentialExpression/allSamples/fibroblasts
fi

cmdi1="cp graphs/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_DESeqResults_perPlate.txt $graphSubdir/differentialExpression/allSamples/iCards/"
cmdi2="cp graphs/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_DESeqResults.txt $graphSubdir/differentialExpression/allSamples/iCards/"
cmdi3="cp graphs/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_nDEtab_perPlate.txt $graphSubdir/differentialExpression/allSamples/iCards/"
cmdi4="cp graphs/differentialExpression/allSamples/iCards/iCard_readFilt_manualFilt_nDEtab.txt $graphSubdir/differentialExpression/allSamples/iCards/"
cmdi5="cp graphs/differentialExpression/allSamples/iCards/iCards_dds_se.rds $graphSubdir/differentialExpression/allSamples/iCards/"
cmdi6="cp -r graphs/differentialExpression/allSamples/iCards/downsampled $graphSubdir/differentialExpression/allSamples/iCards/"
cmdi7="cp graphs/differentialExpression/allSamples/iCards/bootstrap_DESeq2_up.csv $graphSubdir/differentialExpression/allSamples/iCards/"
cmdi8="cp graphs/differentialExpression/allSamples/iCards/bootstrap_DESeq2_down.csv $graphSubdir/differentialExpression/allSamples/iCards/"
cmdi9="cp graphs/differentialExpression/allSamples/iCards/bootstrap_DESeq2_noch.csv $graphSubdir/differentialExpression/allSamples/iCards/"

cmdf1="cp graphs/differentialExpression/allSamples/fibroblasts/fibro_readFilt_manualFilt_DESeqResults.txt $graphSubdir/differentialExpression/allSamples/fibroblasts/"
cmdf2="cp graphs/differentialExpression/allSamples/fibroblasts/fibro_readFilt_manualFilt_nDEtab.txt $graphSubdir/differentialExpression/allSamples/fibroblasts/"
cmdf3="cp graphs/differentialExpression/allSamples/fibroblasts/fibro_dds_se.rds $graphSubdir/differentialExpression/allSamples/fibroblasts/"

cmdc1="cp graphs/differentialExpression/allSamples/controls_both_DESeq2results.txt $graphSubdir/differentialExpression/allSamples"

echo "Starting..." >> $JOURNAL
date >> $JOURNAL
echo "$cmdi1" >> $JOURNAL
eval "$cmdi1"
date >> $JOURNAL
echo "$cmdi2" >> $JOURNAL
eval "$cmdi2"
date >> $JOURNAL
echo "$cmdi3" >> $JOURNAL
eval "$cmdi3"
date >> $JOURNAL
echo "$cmdi4" >> $JOURNAL
eval "$cmdi4"
date >> $JOURNAL
echo "$cmdi5" >> $JOURNAL
eval "$cmdi5"
date >> $JOURNAL
echo "$cmdi6" >> $JOURNAL
eval "$cmdi6"
date >> $JOURNAL
echo "$cmdi7" >> $JOURNAL
eval "$cmdi7"
date >> $JOURNAL
echo "$cmdi8" >> $JOURNAL
eval "$cmdi8"
date >> $JOURNAL
echo "$cmdi9" >> $JOURNAL
eval "$cmdi9"
date >> $JOURNAL
echo "$cmdf1" >> $JOURNAL
eval "$cmdf1"
date >> $JOURNAL
echo "$cmdf2" >> $JOURNAL
eval "$cmdf2"
date >> $JOURNAL
echo "$cmdf3" >> $JOURNAL
eval "$cmdf3"
date >> $JOURNAL
echo "$cmdc1" >> $JOURNAL
eval "$cmdc1"
date >> $JOURNAL
echo "Finished."