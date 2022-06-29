#$ -cwd
#$ -j yes
#$ -V

conda activate R_env

coorFile=coordinate.txt
piFile=piBySite/combinedPiData.txt
fstFile1=Fst/pairwFst.cultivars.spontaneums.weir.fst
fstFile2=Fst/pairwFst.landraces.spontaneums.weir.fst
fstFile3=Fst/pairwFst.cultivars.landraces.weir.fst
blockFile=haplotype_blocks/combinedBlockData.txt
tajFile=TajimasD/combinedTajimaDData.cleaned.txt
muFile=selectiveSweep/combined_report.txt
geneFile=domestication_genes.bed
outName=figure2_new.png
#outName=Mu_plot.png

Rscript plotAll.r $coorFile $piFile $fstFile1 $fstFile2 $fstFile3 $blockFile $tajFile $muFile $geneFile $outName
#Rscript plotMu.r $coorFile $piFile $fstFile1 $fstFile2 $fstFile3 $blockFile $tajFile $muFile $geneFile $outName

rm Rplots.pdf
