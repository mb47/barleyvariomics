#!/bin/bash
#$ -cwd
#$ -j yes
#$ -V
# -pe smp 8
#$ -m bea
#$ -t 1-3
#$ -q ln.q

group[1]=cultivars
group[2]=landraces
group[3]=spontaneums
groups=${group[SGE_TASK_ID]}

conda activate delAllel

#===============================
# thin down vcf file
#===============================

vcfFile=variomeSNPs_final815samples_filteredSNPs.vcf
vcftools --vcf $vcfFile --recode --recode-INFO-all --thin 10000 --out variomeSNPs_thinned10000

#===============================
# convert vcf to plink
#===============================

mkdir $groups

groupDir=groupList
input=variomeSNPs_thinned10000.recode.vcf

vcftools --vcf $input --keep $groupDir/$groups.txt --recode --out $groups/variomeSNPs_thinned10000.$groups
plink --allow-extra-chr --vcf $groups/variomeSNPs_thinned10000.$groups.recode.vcf --out $groups/plink.$groups

#===============================
# calculation
#===============================

plink	--allow-extra-chr \
	--bfile $groups/plink.$groups \
	--r2 \
	--ld-window-r2 0.05 \
	--ld-window 100000 \
	--ld-window-kb 15000 \
	--out $groups/r2.$groups
