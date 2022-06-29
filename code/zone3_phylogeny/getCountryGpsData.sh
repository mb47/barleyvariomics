#!/bin/bash

gpsFile=country_coor.tsv
passport=treeCluster_allChrom_blank.csv

## working dir
#mkdir country_gps

##Get unique list of country in passport file
tail -n +2 $passport | awk -F',' '{print $15}' | sort | uniq > country_gps/listofcountry

## Prepare header
echo -e Country"\t"Lat"\t"Long > country_gps/interFile.tsv

## Search for overlap
IFS=$'\n'
for i in $(cat country_gps/listofcountry)
do
	searchResult=$(grep -w "$i" $gpsFile | wc -l)
	if [ $searchResult == 0 ]
	then
		echo $i >> country_gps/interFile.tsv
		echo "Couldn't find $i"
	else
		awk -F"\t" -v country=$i -v OFS="\t" '$4==country {print $4,$2,$3}' $gpsFile >> country_gps/interFile.tsv
	fi
done


## note
# manually added those countries that can't be found and mv the files as script/misc/country_coord.tsv

