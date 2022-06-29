#$ -cwd
#$ -j yes
#$ -V

conda activate R_env

#head -n 1000000000 lowMAF.toPlot.tsv > dummy.txt
Rscript plotGQ.r dummy.txt              lowMafGQ.subset1t.log.png
#Rscript plotGQ.r lowMAF.toPlot.tsv	lowMafGQ.log.png
Rscript plotGQ.r otherMarkers.toPlot.tsv otherMarkerGQ.log.png
