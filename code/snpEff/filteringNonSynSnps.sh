#$ -cwd
#$ -j yes
#$ -V
#$ -q ln.q

conda activate delAllel

## I/O
bcfFile=variome.815Lines.snpeff.bcf
snpEffDir=snpEff4.4/

## Command
## Keeping the HIGH and MODERATE impact variants
echo "extract nonsyn SNPs location"
bcftools view $bcfFile | \
java -jar $snpEffDir/SnpSift.jar filter "(ANN[*].IMPACT has 'HIGH') | (ANN[*].IMPACT has 'MODERATE')" | \
java -jar $snpEffDir/SnpSift.jar extractFields - CHROM POS POS "ANN[*].EFFECT" \
> nonSynSnpsLocation.txt
echo "done"
