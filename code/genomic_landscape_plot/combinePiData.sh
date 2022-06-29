#!/bin/bash

groups[1]=cultivars
groups[2]=spontaneums
groups[3]=landraces	

dataDir=$1

#extract the column with the pi data from each group's file
for i in {1..3}
do
	group=${groups[i]}
	piFile=$dataDir/piBySite.$group.sites.pi
	echo "processing file $piFile"
	#prepend a header with the group name
	echo "$group" > $group.tmp
	#add the pi data to it
	cut -f3 $piFile | tail -n +2 >> $group.tmp
done

#extract the first two columns with CHROM and POS from one the files
cut -f 1,2 $piFile > chromPos.tmp

#now paste all the group files together
paste *.tmp > combinedPiData.txt

echo "done"
