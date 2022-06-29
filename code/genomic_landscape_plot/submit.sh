#$ -cwd
#$ -j yes
#$ -V

coorFile=coordinate.txt
haploFile=gene_haplotype_count.txt
muFile=selectiveSweep/combined_report.txt
snpFile=delSNPs.bed
geneFile=candidateGenes.bed
tajFile=TajD_toPlot.txt
outName=figure4_a.present.png

conda activate R_env
Rscript plot_figure4.r $coorFile $haploFile $muFile $snpFile $geneFile $tajFile $outName
