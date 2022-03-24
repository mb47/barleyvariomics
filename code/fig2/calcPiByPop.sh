#!/bin/bash

#SBATCH --array=1-3

#This script calculates pi for each of the three populations separately -- cultivars, landraces, wild barleys. 

#the BCF file containing our SNPs
bcfFile=variomeSNPs_final815samples_filteredSNPs.bcf

#a directory that contains the files with the sample IDs, one file per group
groupFileDir=./groupLists

#we have the following populations:
populations[1]=cultivars
populations[2]=wildBarleys
populations[3]=landraces

#pick the appropriate one for this run of the script
population=${populations[SLURM_ARRAY_TASK_ID]}
#this is the specific file with the sample IDs for this population
popFile=$groupFileDir/$population.txt

echo "running vcftools window-pi analysis for population $population"
echo "start time `/bin/date`"

vcftools \
--bcf $bcfFile \
--window-pi 10000 \
--keep $popFile \
--out pi_10k_windows.$population


echo "finished vcftools analysis"
echo "end time `/bin/date`"