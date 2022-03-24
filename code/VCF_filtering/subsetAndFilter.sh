#Script for subsetting original BCF file to the final 815 samples and filtering it

#file with names of samples to include
sampleList=815accessions.txt

#this BCF file contains the unfiltered variants for the set of 879 original samples
bcfFile=/mnt/shared/projects/barley/201710_barleyVariome/rerunSep2018/jointGenotyper/jointGenotyperSNPs_allSamples.bcf

#prefix for output files
prefix=variomeSNPs_final815samples

echo "extract non-monomorphic variants for the subset of 815 samples"
bcftools view \
-o $prefix.bcf \
--output-type b \
-e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' \
--samples-file $sampleList \
$bcfFile

#the previous step has produced a BCF file
#we need VCF format to do the filtering though
echo "convert BCF file to VCF"
bcftools view \
-o $prefix.vcf \
--output-type v \
$prefix.bcf

#now apply the filtering
echo "filter the VCF file"
java \
utils.snps.vcf.parsing.gatk.ConvertVCF_GATK_ToSpreadsheetFormat \
$prefix.vcf \
GATK_VCF_filter_config.txt 

#the filtering has produced a VCF file
#we need BCF format of this for a number of purposes
echo "convert filtered VCF file to BCF"
bcftools view \
-o $prefix"_filteredSNPs.bcf" \
--output-type b \
$prefix"_filteredSNPs.vcf"

echo "index the final BCF file"
bcftools \
index \
$prefix"_filteredSNPs.bcf"

#generate stats
bcftools \
stats \
$prefix"_filteredSNPs.bcf"

echo "workflow complete"

