#script for extracting and plotting QUAL values from the exome capture BCF file
#split by on- and off-target variants 

#the BCF file we want to analyse
bcfFile=variomeSNPs_final815samples_filteredSNPs.bcf

#this BED file contains the on-target regions (BLAST hit regions for capture probes)
bedFileExcapTargetRegions=BLASTbarleyExomeSeqsVsFinalPseudomoleculesSep2015.topHSPsOnly.wholeChrCoords.bed

#prefix for all output files
prefix=variomeSNPs_strictFiltersNoMAF_final
#a BED file with all the variant positions
bedFileAllSNPs=$prefix.bed

#get processing

echo "make the BED file"  #chr, pos (0-based), end pos (1-based), id
bcftools query \
-f'%CHROM\t%POS0\t%END\n' \
$bcfFile \
> $bedFileAllSNPs

echo "compute overlap with exome capture probes"
bedtools intersect \
-wa \
-sorted \
-a $bedFileAllSNPs \
-b $bedFileExcapTargetRegions \
-u \
> $prefix.onTargetSNPs.bed

echo "number of on-target SNPs:"
wc -l $prefix.onTargetSNPs.bed

echo "subtract on target SNPs from whole file to get off target SNPs"
bedtools subtract \
-sorted \
-a $bedFileAllSNPs \
-b $prefix.onTargetSNPs.bed \
> $prefix.offTargetSNPs.bed

echo "number of off-target SNPs:"
wc -l $prefix.offTargetSNPs.bed

echo "extract VCF file for on target SNPs"
bcftools view \
--regions-file $prefix.onTargetSNPs.bed \
-O v \
-o $prefix.onTargetSNPs.vcf \
$bcfFile

echo "extract VCF file for off target SNPs"
bcftools view \
--regions-file $prefix.offTargetSNPs.bed \
-O v \
-o $prefix.offTargetSNPs.vcf \
$bcfFile

echo "extract QUAL values, on target SNPs"
cut -f 6 $prefix.onTargetSNPs.vcf \
| grep -v "^#" \
| tail -n +2 \
> $prefix.onTargetSNPs.QUAL.txt

echo "extract QUAL values, off target SNPs"
cut -f 6 $prefix.offTargetSNPs.vcf \
| grep -v "^#" \
| tail -n +2 \
> $prefix.offTargetSNPs.QUAL.txt

echo "format the QUAL data for plotting"
#first prefix each file with a category label
sed 's/^/onTarget\t/g' $prefix.onTargetSNPs.QUAL.txt > tmp1.txt
sed 's/^/offTarget\t/g' $prefix.offTargetSNPs.QUAL.txt > tmp2.txt
#then concatenate these
cat tmp1.txt tmp2.txt > allSNPs.QUAL.txt
#clean up
rm tmp1.txt
rm tmp2.txt

echo "plot the QUAL data"
Rscript \
violinPlot.r \
allSNPs.QUAL.txt \
allSNPs.QUAL.png \
"SNP quality score"

echo "extract read depth data for onTarget variants"
vcftools \
--vcf variomeSNPs_strictFiltersNoMAF_final.onTargetSNPs.vcf \
--get-INFO DP \
--out variomeSNPs_strictFiltersNoMAF_final.onTargetSNPs

echo "extract read depth data for offTarget variants"
vcftools \
--vcf variomeSNPs_strictFiltersNoMAF_final.offTargetSNPs.vcf \
--get-INFO DP \
--out variomeSNPs_strictFiltersNoMAF_final.offTargetSNPs

#output files from this are suffixed with .INFO
#data we want is in col 5 with header (discard the latter)
#on target first
cut -f 5 variomeSNPs_strictFiltersNoMAF_final.onTargetSNPs.INFO \
| tail -n +2 \
| sed 's/^/onTarget\t/g' \
> tmpOn.txt
#off target
cut -f 5 variomeSNPs_strictFiltersNoMAF_final.offTargetSNPs.INFO \
| tail -n +2 \
| sed 's/^/offTarget\t/g' \
> tmpOff.txt
#concat the files
cat tmpOn.txt tmpOff.txt > allSNPs.DP.txt
#clean up
rm tmpOn.txt
rm tmpOff.txt

echo "plot the DP data"
Rscript \
violinPlot.r \
allSNPs.DP.txt \
allSNPs.DP.png \
"Read depth"

echo "done"


