#!/bin/bash

vcfFile=$1
group1=$2
group2=$3
groupFileDir=$4
outDir=$5

echo "running vcftools weir-fst-pop analysis"
echo "start time `/bin/date`"

vcftools \
--vcf $vcfFile \
--weir-fst-pop $groupFileDir/$group1.txt \
--weir-fst-pop $groupFileDir/$group2.txt \
--out $outDir/pairwFst.$group1.$group2

echo "finished vcftools diversity analysis"
echo "end time `/bin/date`"
