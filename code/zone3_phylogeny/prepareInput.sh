#!/bin/bash

csvFile=treeCluster_allChrom_blank.csv
coordFile=../../../script/misc/country_coor_editted.tsv

chrom[1]=chr1H
chrom[2]=chr2H
chrom[3]=chr3H
chrom[4]=chr4H
chrom[5]=chr5H
chrom[6]=chr6H
chrom[7]=chr7H
group[1]=Cultivar
group[2]=Landrace
group[3]=Spontaneum

mkdir dataframe

for i in {1..7}
do
	chrom=${chrom[i]}
	for a in {1..3}
	do
		group=${group[a]}
		./getDataFrame.sh $csvFile $coordFile $chrom $group 
	done
done

mv df_*.tsv dataframe/
