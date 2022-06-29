#!/bin/bash
#### This script takes the raw treeCluster grouping information and count/transform it to the input format for scatterpie plot in R
#### I/O
csvFile=$1
gpsFile=$2
chrom=$3
origin=$4
echo "----Analysing $chrom"
# Determine the column number based on chromosome
if   [ "$chrom" == "chr1H" ];then column=3
elif [ "$chrom" == "chr2H" ];then column=4
elif [ "$chrom" == "chr3H" ];then column=5
elif [ "$chrom" == "chr4H" ];then column=6
elif [ "$chrom" == "chr5H" ];then column=7
elif [ "$chrom" == "chr6H" ];then column=8
elif [ "$chrom" == "chr7H" ];then column=9
else echo "**Chromosome doesn't exist"
fi

#### Set IFS (separator) so the country names with spaces can be read correctly
IFS=$'\n'

#### Get temp Country list
# Set IFS (separator) so the country names with spaces can be read correctly
echo "--------Retrieving list of countries"
IFS=$'\n'
echo "Country" > countryList.temp
awk -v origin=$origin -F',' '$12==origin {print $15}' $csvFile | sort | uniq >> countryList.temp
#### Get temp Group list
echo "--------Retrieving list of groups"
awk -v origin=$origin -F',' -v col=$column '$12==origin {print $col}' $csvFile | sort | uniq > tempGroupList

#### Get country GPS
echo "--------Retrieving countries' GPS"
echo -e Lat"\t"Long > coord.temp
for i in $(tail -n +2 countryList.temp)
do
	awk -F'\t' -v OFS="\t" -v country=$i '$1==country {print $2,$3}' $gpsFile >> coord.temp
done

#### Loop and get barley count data
for i in $(cat tempGroupList)
do
	echo "----------Writing group$i for $chrom"
	echo Group$i > group$i.temp
	for j in $(tail -n +2 countryList.temp)
	do
		awk -F',' -v country=$j -v origin=$origin -v col=$column '$12==origin && $15==country {print $col}' $csvFile |\
			grep -w $i | wc -l >> group$i.temp
	done
done

#### Loop and get total count data as radius
echo "--------Writing total count by country"
echo Radius > radius.temp
for i in $(tail -n+2 countryList.temp)
do
	awk -F',' -v country=$i -v origin=$origin '$15==country && $12==origin' $csvFile | wc -l >> radius.temp
done

#### Paste results into data frame
echo "--------Combining all data"
paste countryList.temp coord.temp radius.temp group*.temp > df_$chrom.$origin.tsv

#### Remove all intermediate temp files
rm countryList.temp
rm tempGroupList
rm coord.temp
rm radius.temp
rm group*.temp
