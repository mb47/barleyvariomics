#Script for computing Fst, pairwise comparisons between the three populations

#the VCF file containing our SNPs
vcfFile=/mnt/shared/projects/barley/201710_barleyVariome/rerunSep2018/subsets/strictFiltersNoMAF/final/variomeSNPs_strictFiltersNoMAF_final.vcf


#we need the following comparisons:
#wild versus cultivars, wild versus landraces and landraces versus cultivars
comparisons[1]=wildBarleys.cultivars
comparisons[2]=wildBarleys.landraces
comparisons[3]=landraces.cultivars
#the files containing the line names are labelled as follows:
# cultivars.txt
# landraces.txt
# wildBarleys.txt

groupFileDir=/mnt/shared/projects/barley/201710_barleyVariome/rerunSep2018/analysis/groupLists/

#iterate over the comparisons
for j in {1..3}
do
	comparison=${comparisons[j]}
	group1=`echo $comparison | cut -f1 -d "."`
	group2=`echo $comparison | cut -f2 -d "."`
	echo -e "\nsubmitting job with comparison $comparison"
	echo "groups: $group1 $group2"		
	qsub compFstPairwise.sh $vcfFile $group1 $group2 $groupFileDir
done

echo "all jobs submitted"

