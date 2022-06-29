vcfFile=$1
population=$2
groupFileDir=$3
popFile=$groupFileDir/$population.txt
outDir=$4

echo "running vcftools freq analysis"
echo "vcfFile = $vcfFile"
echo "start time `/bin/date`"

#vcftools \
#--vcf $vcfFile \
#--freq \
#--keep $popFile \
#--out alleleFreqs.$population


#echo "running vcftools site-pi analysis"
#echo "start time `/bin/date`"

vcftools \
--vcf $vcfFile \
--site-pi \
--keep $popFile \
--out $outDir/piBySite.$population


#echo "finished vcftools diversity analysis"
#echo "end time `/bin/date`"
