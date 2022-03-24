
#$ -q ln.q

bcfInputFile=variomeSNPs_final815samples_filteredSNPs.bcf

bcfOutputFile=variomeSNPs_final815samples_filteredSNPs.reannotated.bcf

snpEffDir=/mnt/shared/scratch/synology/mb40521/apps/snpeff/snpEff_4_3

echo "annotate BCF file"
bcftools view $bcfInputFile | \
java -Xmx20g -jar $snpEffDir/SnpSift.jar rmInfo - EFF | \
java -Xmx20g -jar $snpEffDir/snpEff.jar -lof Morex2017WithBart1_0 - | \
bcftools view -O b -o $bcfOutputFile

echo "pipeline complete"
