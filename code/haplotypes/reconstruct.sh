#!/bin/sh

#============================================================================================================================
# USAGE
# ./reconstruct.sh [gene/transcript_name] [accession list] [reference genome] [gene feature] [variants] [output directory]
#
#   accession list: accession in single lines
#   reference genome: fasta file
#   gene feature: gff3 file with exon feature
#   variants: in bcf or vcf.gz (bgzipped)
#             by default, the ALT allele is used when encountering heterozygous variants (bcftools consensus -H A option)
#
#============================================================================================================================
# List of input files
# i is the gene/transcript name
i=$1
accessionList=$2
refGenome=$3
gffFile=$4
bcfFile=$5
outputDir=$6
#gffFile=177genes.longestTranscriptOnly.gff3
#accessionList=23SampleList.txt
#bcfFile=variomeSNPs_strictFiltersNoMAF_final.bcf
#bcfFile=strictFiltersNoMAF_SNPs.indels.vcf.gz
#outputDir=177GenesSequence

# Get sequence coordinate into bed format
grep $i $gffFile | awk -v OFS='\t' '$3 == "exon" {print $1,$4-1,$5}' > $i.temp.bed
# Loop through accession list and retain transcript variations across accessions
if [ $(cat $i.temp.bed | wc -l) -gt 0 ]; then
	echo "analysing $i"
	# Get the reference sequence into fasta format
	cat $refGenome | seqkit subseq --bed $i.temp.bed --quiet | sed 's/\.//g' | sed 's/://g' | sed 's/_/:/g' > $i.ref.temp.fa
	# Build consensus sequence all accessions
	for acc in $(cat $accessionList)
	do
		# Bcftools consensus
		cat $i.ref.temp.fa | bcftools consensus -H A $bcfFile --sample $acc | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' > $i.transcript.temp.fa
		# Gigu script convert sequence format to single fasta entry
		strand=$(grep $i $gffFile | head -n 1 | awk '{print $7}')
		python3.7 gigu_script.py $i.transcript.temp.fa $i $acc $strand >> $outputDir/$i.combined.fasta
	done
	rm $i.ref.temp.fa
	rm $i.transcript.temp.fa
else
	echo "no gff record found for: $i in $gffFile"
fi
rm $i.temp.bed
