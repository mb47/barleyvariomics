#!/bin/bash

fastaFile=$1
geneName=$2

# remove the -- separator and convert multiline fasta to single line fasta
grep -v "^--" $fastaFile | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' > $geneName.temp.fa
totalCount=$(grep -v "^>" $geneName.temp.fa | sort | uniq | wc -l)

# grep fasta with different origins
groupDir=groupList
group[1]=cultivars
group[2]=landraces
group[3]=spontaneums

for i in {1..3}
do
	groupFile=$groupDir/${group[i]}.txt
	for accession in $(cat $groupFile)
	do
		grep -A 1 "$accession.p1" $geneName.temp.fa >> $geneName.directory/$geneName.${group[i]}.fa
	done
done

# count number of uniq aa haplotype
cultivarCount=$(grep -v "^>" $geneName.directory/$geneName.cultivars.fa | sort | uniq | wc -l)
landraceCount=$(grep -v "^>" $geneName.directory/$geneName.landraces.fa | sort | uniq | wc -l)
spontCount=$(grep -v "^>" $geneName.directory/$geneName.spontaneums.fa | sort | uniq | wc -l)

echo -e "#cultivars	landraces	spontaneums	total" > $geneName.directory/haplotypeCount.txt
echo -e $cultivarCount'\t'$landraceCount'\t'$spontCount'\t'$totalCount >> $geneName.directory/haplotypeCount.txt

echo "finished analysing $geneName"

# remove the one line temp file
rm $geneName.temp.fa
