#$ -cwd
#$ -j yes
#$ -V
#$ -q ln.q
#$ -t 1-3

vcfFile=$1
tempOutDir=$2
scriptDir=script
groupDir=fullAccessionInformation/groupList
group[1]=cultivars
group[2]=landraces
group[3]=spontaneums
groups=${group[SGE_TASK_ID]}

## Calculate site pi

# conda activate delAllel
# echo
# echo "Analysis results in $tempOutDir"

# echo
# echo "Calculating site pi"
# $scriptDir/compDivStatsByPop.sh $vcfFile $groups $groupDir $tempOutDir

# echo
# echo "Combining site pi results"
# $scriptDir/combinePiData.sh $tempOutDir

# conda deactivate
# conda activate R_env
# echo
# echo "plotting site pi results"
# piFile=$tempOutDir/combinedPiData.txt
# piPng=$tempOutDir/combinedPiData.scaled.png
# Rscript $scriptDir/plotCombinedData.r $piFile $piPng 10000 2 "" "Nucleotide Diversity"


#echo "Calculate pi by gene"

subset[1]=cultivars
subset[2]=landraces
subset[3]=spontaneums
subsets=${subset[SGE_TASK_ID]}

# $scriptDir/getSnpList.sh misc/filtered.gene.gtf $tempOutDir/piBySite.$subsets.sites.pi $tempOutDir/piByGene.$subsets

echo "Combine pi by gene data"

groups[1]=cultivars
groups[2]=spontaneums
groups[3]=landraces	

dataDir=$tempOutDir

#extract the column with the pi data from each group's file
for i in {1..3}
do
	group=${groups[i]}
	piFile=$dataDir/piByGene.$group.txt
	echo "processing file $piFile"
	#prepend a header with the group name
	echo "$group" > $tempOutDir/$group.tmp
	#add the pi data to it
	cut -f4 $piFile >> $tempOutDir/$group.tmp
done

# Append header to chromPos.tmp
echo -e CHROM'\t'GENE'\t'POS > $tempOutDir/chromPos.tmp

#extract the first two columns with CHROM and POS from one the files
cut -f 1,2,3 $piFile >> $tempOutDir/chromPos.tmp

#now paste all the group files together
paste $tempOutDir/*.tmp > $tempOutDir/combinedPiData.txt
rm $tempOutDir/*.tmp

echo "done"
