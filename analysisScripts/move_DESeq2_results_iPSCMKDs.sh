#!/bin/bash

# move_DESeq2_results_iPSCMKDs.sh moves (copies) previously run move_DESeq2_results from iPS-CM-KDs into the target
# graphs directory. Expects to be running in projectDir. 

projectDir=$1
graphSubdir=$2
JOURNAL=$3

if [ ! -d $graphSubdir ]; then
    mkdir $graphSubdir
fi
if [ ! -d $graphSubdir/iPS-CM-KDs ]; then
    mkdir $graphSubdir/iPS-CM-KDs
fi
if [ ! -d $graphSubdir/iPS-CM-KDs/differentialExpression ]; then
    mkdir $graphSubdir/iPS-CM-KDs/differentialExpression
fi


cmdi1="cp graphs/iPS-CM-KDs/differentialExpression/DE_results.txt $graphSubdir/iPS-CM-KDs/differentialExpression/"

echo "Starting..." >> $JOURNAL
date >> $JOURNAL
echo "$cmdi1" >> $JOURNAL
eval "$cmdi1"
date >> $JOURNAL

echo "Finished."