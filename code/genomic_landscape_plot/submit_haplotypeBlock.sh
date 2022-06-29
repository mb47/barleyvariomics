#!/bin/bash
#$ -cwd
#$ -j yes
#$ -V
#$ -t 1-3

## Convert vcf to plink format for 3 origins

vcfFile=$1
outputDir=$2
groupDir=groupList
group[1]=cultivars
group[2]=landraces
group[3]=spontaneums
groups=${group[SGE_TASK_ID]}

conda activate delAllel

mkdir $outputDir/$groups
vcftools --vcf $vcfFile --keep $groupDir/$groups.txt --recode --out $outputDir/$groups/strictFiltersNoMAF_SNPs.$groups
plink --allow-extra-chr --vcf $outputDir/$groups/strictFiltersNoMAF_SNPs.$groups.recode.vcf --out $outputDir/$groups/plink.$groups

# plink parameter explain
# Estimate haplotype block and generate .block file
# --allow-extra-chr: for correct reading of the chromosome name
# no-pheno-req: is ommiting the need of phenotype data
# --recode beagle: to beagle format for further haplotype characterization
# --blocks-max-kb: expand the window size to calculate the pairwise LD for all marker pairs

## Calculate haplotype block

plink --allow-extra-chr --bfile $outputDir/$groups/plink.$groups --recode beagle --blocks no-pheno-req --blocks-max-kb 800000 --out $outputDir/$groups/haploblock.$groups

## Summarise results

./blockSummary.sh $outputDir/$groups/haploblock.$groups.blocks.det > $outputDir/$groups/haploblock.$groups.blocks.summary

conda deactivate
conda activate R_env

## Plot results

./formatDetFile.sh $outputDir/$groups/haploblock.$groups.blocks.det $outputDir/$groups
Rscript plotHaploBlock.r $outputDir/$groups/formatted.block.det $groups 5000

conda deactivate
