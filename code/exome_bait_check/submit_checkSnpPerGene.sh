#$ -cwd
#$ -j yes
#$ -V

inFile=geneAndRegions.txt
bcfFile=filteredSnpsAndIndel.vcf
outFile=group3mSnpsIndelsByGene.tsv

echo "GENE	NO.VAR	NO.noALT	NO.SNP	NO.MNP	NO.INDEL	NO.OTHER	NO.MULTIALLEL	NO.MULTIALLELSNP" > $outFile

conda activate delAllel
while IFS= read -r line; do
    gene=$(echo $line | awk '{print $1}')
    region=$(echo $line | awk '{print $2}')
	bcftools view --regions $region $bcfFile | bcftools stats > tempFile
	noVar=$(grep "number of records:"	tempFile | awk '{print $6}')
	noALT=$(grep "number of no-ALTs:"	tempFile | awk '{print $6}')
	noSNP=$(grep "number of SNPs:"		tempFile | awk '{print $6}')
	noMNP=$(grep "number of MNPs:"		tempFile | awk '{print $6}')
	noIND=$(grep "number of indels:"		tempFile | awk '{print $6}')
	noOTH=$(grep "number of others:"		tempFile | awk '{print $6}')
	noMAS=$(grep "number of multiallelic sites:"		tempFile | awk '{print $7}')
	noMASNP=$(grep "number of multiallelic SNP sites:"		tempFile | awk '{print $8}')
	echo -e "$gene\t$noVar\t$noALT\t$noSNP\t$noMNP\t$noIND\t$noOTH\t$noMAS\t$noMASNP" >> $outFile
done<$inFile

rm tempFile

conda deactivate

## bcftools stats example output
#SN      0       number of samples:      879
#SN      0       number of records:      530
#SN      0       number of no-ALTs:      0
#SN      0       number of SNPs: 466
#SN      0       number of MNPs: 0
#SN      0       number of indels:       86
#SN      0       number of others:       35
#SN      0       number of multiallelic sites:   86
#SN      0       number of multiallelic SNP sites:       30
