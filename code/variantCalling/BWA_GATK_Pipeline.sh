#!/bin/bash


##########################################################################################################
#This is a pipeline script for aligning PE Illumina reads from a single sample 
#to a FASTA reference sequence and then calling variants on it using the GATK.
#It follows the GATK best practice in broad terms and assumes there is no existing
#truth data available in terms of known SNPs (as is the case for most non-model organisms). It
#creates its own truth data in an initial run of the HaplotypeCaller and uses these as the basis for 
#the base quality score recalibration. 
#
#This pipeline does not include the production of final SNPs in a population of samples.
#For this, several runs of this pipeline on different samples are necessary, followed by the joint genotyper being run over the set of GVCF files produced. 

#Tools required: BWA, GATK, vcflib, samtools, bamtools
#This script will put the required versions of these on the PATH, before any other existing versions. 

#Micha Bayer, The James Hutton Institute
##########################################################################################################

##########################################################################################################
#configure the CLASSPATH and PATH so all the right versions of the required tools/libraries are available
##########################################################################################################
export PATH=:/mnt/apps/bamtools/bamtools/bin:/mnt/apps/bwa/bwa-0.7.10:/mnt/apps/samtools/0.1.19:/mnt/apps/vcflib/20140627/bin/:/mnt/apps/ngsQA/utils/scripts:$PATH

export CLASSPATH=/mnt/apps/ngsQA/utils/lib/utils.jar:$CLASSPATH

#path to the GATK jar file
gatkJar=/mnt/apps/gatk/3.4.0/GenomeAnalysisTK.jar

##########################################################################################################
#variables -- mostly from command line
##########################################################################################################
#check for the correct number of args
if [ $# -ne 6 ]
then
    echo "Error in $0 - Invalid Argument Count"
    echo "Syntax: $0  <R1 FASTQ file, gzipped> <R2 FASTQ file, gzipped> <sample name> <full path to reference sequence, FASTA> <min. alignment score> <numThreads>"
    exit
fi

#R1 FASTQ file, gzipped
R1=$1
#R2 FASTQ file, gzipped
R2=$2
#sample name -- can be different from file name
sampleName=$3
#full path to reference sequence, FASTA
refseq=$4
#min. alignment score
minAlignmentScore=$5
#number of threads
numThreads=$6


##########################################################################################################
#unzip the FASTQ files
##########################################################################################################
echo -e "\n========================================="
echo "processing sample $sampleName"
echo "start time `/bin/date`"
echo -e "========================================="
echo "unzipping R1 file"
baseNameR1=`basename $R1 .gz`
pigz -p $numThreads -cd $R1 > $baseNameR1

echo "unzipping R2 file"
baseNameR2=`basename $R2 .gz`
pigz -p $numThreads -cd $R2 > $baseNameR2


##########################################################################################################
#map the reads with BWA-MEM
##########################################################################################################
echo -e "\n========================================="
echo "mapping with BWA"
echo "start time `/bin/date`"
echo -e "========================================="
bwa mem $refseq $baseNameR1 $baseNameR2 -t $numThreads -R "@RG\tID:$sampleName\tSM:$sampleName\tPL:Illumina" 2> $sampleName.bwa.log > $sampleName.sam

echo "convert SAM to BAM, filter and sort"
samtools view -F 4 -S -b -h -u $sampleName.sam |  \
bamtools-2.2.3 filter -tag "AS:>=$minAlignmentScore" | \
samtools sort -l 9 - $sampleName.sorted  

#get rid of the SAM file immmediately
rm $sampleName.sam

echo "done mapping"


##########################################################################################################
#duplicate removal
##########################################################################################################
echo -e "\n========================================="
echo "removing duplicates"
echo "start time `/bin/date`"
echo -e "========================================="
samtools rmdup $sampleName.sorted.bam $sampleName.rmduped.bam 2> $sampleName.rmdup.log
samtools index $sampleName.rmduped.bam
echo "done with rmdup"


##########################################################################################################
#indel realignment
##########################################################################################################
echo -e "\n========================================="
echo "Indel realignment step 1 - generating target interval list"
echo "start time `/bin/date`"
echo -e "========================================="
java -jar $gatkJar -T RealignerTargetCreator -R $refseq -I $sampleName.rmduped.bam -o $sampleName.target_intervals.list -nt $numThreads

echo -e "\n========================================="
echo "Indel realignment step 2 - realigning"
echo "start time `/bin/date`"
echo -e "========================================="
java -jar $gatkJar -T IndelRealigner -R $refseq -I $sampleName.rmduped.bam -targetIntervals $sampleName.target_intervals.list -o $sampleName.realigned.bam
#index the newly created realigned BAM file
samtools index $sampleName.realigned.bam


##########################################################################################################
#initial run of the HaplotypeCaller -- this generates the first set of variants we will use as truth data for the base quality score recalibration
##########################################################################################################
echo -e "\n========================================="
echo "Running the HaplotypeCaller - first run"
echo "start time `/bin/date`"
echo -e "========================================="
java -jar $gatkJar -T HaplotypeCaller -R $refseq -I $sampleName.realigned.bam -o $sampleName.initial.variants.vcf -dontUseSoftClippedBases -nct $numThreads


##########################################################################################################
#filter the initial variants to remove poor quality calls with a QUAL of
#less than 20 - this will become our truth set for the BQSR step
##########################################################################################################
echo "filtering variants for >QUAL20"
vcffilter -f "QUAL > 20" $sampleName.initial.variants.vcf > $sampleName.initial.variants.filteredQ20.vcf


##########################################################################################################
#base quality score recalibration
##########################################################################################################
#produce the recalibration table
echo -e "\n========================================="
echo "producing the recalibration table"
echo "start time `/bin/date`"
echo -e "========================================="
java -jar $gatkJar -T BaseRecalibrator -R $refseq -I $sampleName.realigned.bam -knownSites $sampleName.initial.variants.filteredQ20.vcf -o $sampleName.recal.table -nct $numThreads

#having produced the table, now recalibrate the data
echo -e "\n========================================="
echo "recalibrating the data"
echo "start time `/bin/date`"
echo -e "========================================="
java -jar $gatkJar -T PrintReads -R $refseq -I $sampleName.realigned.bam -BQSR $sampleName.recal.table -o $sampleName.recalibrated.bam -nct $numThreads


##########################################################################################################
#second run of the HaplotypeCaller -- this generates the final set of variants we will use as input for the joint genotyper
##########################################################################################################
echo -e "\n========================================="
echo "Running the HaplotypeCaller - second run"
echo "start time `/bin/date`"
echo -e "========================================="
java -jar $gatkJar -T HaplotypeCaller -R $refseq -I $sampleName.recalibrated.bam -o $sampleName.final.variants.g.vcf -ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -dontUseSoftClippedBases -nct $numThreads


##########################################################################################################
#clean up any redundant files
##########################################################################################################
echo -e "\ncleaning up redundant files"
rm $sampleName.sorted.bam 
rm $sampleName.rmduped.bam
rm $sampleName.rmduped.bam.bai
rm $sampleName.realigned.bai
rm $sampleName.realigned.bam
rm $sampleName.realigned.bam.bai
#also delete the uncompressed FASTQ files 
#these were created by reading the compressed files only (rather than overwriting them)
rm $baseNameR1
rm $baseNameR2
#various other files that are now redundant
rm $sampleName.initial.variants.vcf
rm $sampleName.initial.variants.filteredQ20.vcf
rm $sampleName.recal.table
rm $sampleName.target_intervals.list
rm $sampleName.bwa.log

##########################################################################################################
#done
##########################################################################################################
echo -e "\n========================================="
echo "PIPELINE COMPLETE"
echo "END TIME `/bin/date`"
echo -e "========================================="

