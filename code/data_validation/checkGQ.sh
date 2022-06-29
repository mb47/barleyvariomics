#$ -cwd
#$ -j yes
#$ -V
#$ -q ln.q

## filtering and keeping only those with MAF < 0.05 (i.e., < 40.75 --> 41 individual) 
conda activate delAllel
vcfFile=variomeSNPs_final815samples_filteredSNPs.bcf
bcftools view -i 'COUNT(GT="AA")<41 || COUNT(GT="RR")<41' $vcfFile > lowMAF.vcf
bcftools view -i 'COUNT(GT="AA")>=41 && COUNT(GT="RR")>=41' $vcfFile > othersMarkers.vcf

## Extract GQ format info for low MAF markers
vcftools --vcf lowMAF.vcf --extract-FORMAT-info GQ --out lowMAF.GQ
vcftools --vcf othersMarkers.vcf --extract-FORMAT-info GQ --out othersMarkers.GQ

## check stats of MAF filtered vcf
bcftools stats othersMarkers.vcf > othersMarkers.stats

## check stats of low MAF vcf
bcftools stats lowMAF.vcf > lowMAF.stats
