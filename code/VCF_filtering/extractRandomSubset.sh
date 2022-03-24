#Script for extraction of random subset of SNPs from main set

#parameters
dataDir=./analysis
baseName=jointGenotyperSNPs_allSamples
inputVCF=$dataDir/$baseName"_wholeChrCords_annotated_filteredSNPs.vcf"

refSeq=150831_barley_pseudomolecules.fasta

prefix=SNPs_allSamples_filtered_subsample5pct
outputVCF=$prefix.vcf
fraction=0.05

#path to the GATK jar file
gatkJar=./apps/GATK/3.4.0/GenomeAnalysisTK.jar

#use the GATK for this
java -jar $gatkJar \
-T SelectVariants \
-R $refSeq \
-V $inputVCF \
-o $outputVCF \
--select_random_fraction $fraction \
--nonDeterministicRandomSeed
