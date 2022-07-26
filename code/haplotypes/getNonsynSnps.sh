#$ -cwd
#$ -j yes
#$ -V
#$ -q ln.q

snpEffDir=snpEff4.4/
bcfFile=snpEff.815lines.full.bcf
bedFile=zone3coordinates.bed
datasets=variomeIndels

conda activate delAllel

echo "extract nonsyn variants in zone3"
bcftools view --regions-file $bedFile $bcfFile | \
java -jar $snpEffDir/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') | (ANN[*].EFFECT has 'start_lost') | (ANN[*].EFFECT has 'stop_lost') | (ANN[*].EFFECT has 'stop_gained')" > zone3.nonsynonymousSnps.vcf
