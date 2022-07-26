#!/bin/bash

#geneName is i
i=$1
workingDir=$2

# Identify transcript ORF
TransDecoder.LongOrfs -S -t $workingDir/$i.combined.fasta
cp $i.combined.fasta.transdecoder_dir/longest_orfs.pep $workingDir/$i.longest_orfs.pep
cp $i.combined.fasta.transdecoder_dir/longest_orfs.cds $workingDir/$i.longest_orfs.cds
rm -r $i.combined.fasta.transdecoder_dir*
rm pipeliner*
# Grep all the p1 sequences
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $workingDir/$i.longest_orfs.pep > $workingDir/$i.longest_orfs.single.pep
grep -A 1 -w "p1" $workingDir/$i.longest_orfs.single.pep > $workingDir/$i.longest_orfs.p1.pep
rm $workingDir/$i.longest_orfs.single.pep
# Align pep sequences
mafft --inputorder --anysymbol --auto $workingDir/$i.longest_orfs.p1.pep > $workingDir/$i.longest_orfs.p1.aligned.pep
