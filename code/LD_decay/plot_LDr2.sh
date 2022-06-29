#$ -cwd
#$ -j yes
#$ -V


## Combine data in to one file
#echo -e R2'\t'DIST'\t'GROUP > combinedR2Data.txt
#awk -v OFS="\t" '{print $7,$8,"\"#4575b4\""}' <(tail -n +2 cultivars/formatted.r2.cultivars.ld) >> combinedR2Data.txt
#awk -v OFS="\t" '{print $7,$8,"\"#fc8d59\""}' <(tail -n +2 landraces/formatted.r2.landraces.ld) >> combinedR2Data.txt
#awk -v OFS="\t" '{print $7,$8,"\"#d73027\""}' <(tail -n +2 spontaneums/formatted.r2.spontaneums.ld) >> combinedR2Data.txt
 
## Plot
#conda activate R_env
Rscript plotDecayLD.r combinedR2Data.txt combinedR2plot_new.png
