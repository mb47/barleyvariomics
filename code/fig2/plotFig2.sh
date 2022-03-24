#!/bin/bash

#SBATCH --mem=5G
#SBATCH -o slurm-%x_%A.out   

datadir=`cwd`

conda activate R_env

coorFile=coordinate.txt
piFile=../pi_10k_windows_combined.txt
fstFile1=$datadir/Fst/pairwFst.cultivars.spontaneums.weir.fst
fstFile2=$datadir/Fst/pairwFst.landraces.spontaneums.weir.fst
fstFile3=$datadir/Fst/pairwFst.cultivars.landraces.weir.fst
blockFile=$datadir/haplotype_blocks/combinedBlockData.txt
geneFile=domestication_genes.bed
outName=figure2_3Mar2022.png


Rscript plotAll.r $coorFile $piFile $fstFile1 $fstFile2 $fstFile3 $blockFile $geneFile $outName
