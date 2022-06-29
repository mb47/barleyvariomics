#$ -cwd
#$ -j yes
#$ -V
#$ -q ln.q

snpEffDir=snpEff4.4/
vcfFile=$1
outDir=$2
outHtml=$outDir/$3.html
outFile=$outDir/$3.bcf

conda activate delAllel

echo "reannotate old BCF file"
bcftools view $vcfFile | \
java -jar $snpEffDir/SnpSift.jar rmInfo - EFF | \
java -jar $snpEffDir/snpEff.jar ann -stats $outHtml -lof Morex2017WithBart1_0.stringtie.transdecoder.homology -fastaProt - | \
bcftools view -O b -o $outFile

echo "index new BCF file"
bcftools index $outFile

echo "extract zone3 nonsyn SNPs"
bcfFile=$outFile
bedFile=zone3coordinates.bed
bcftools view --regions-file $bedFile $bcfFile | \
java -jar $snpEffDir/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') | (ANN[*].EFFECT has 'start_lost') | (ANN[*].EFFECT has 'stop_lost') | (ANN[*].EFFECT has 'stop_gained')" > $outDir/zone3.nonsynonymousSnps.vcf

echo "summarise result by genes"
cat $outDir/zone3.nonsynonymousSnps.vcf | java -jar groupbygene.jar > $outDir/nonsynonymousSnps_summary.tsv

echo "pipeline complete"
conda deactivate

