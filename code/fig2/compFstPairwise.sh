#Script for computing Fst, pairwise comparisons between two populations at a time

vcfFile=$1
group1=$2
group2=$3
groupFileDir=$4


echo "running vcftools weir-fst-pop analysis"
echo "start time `/bin/date`"

vcftools \
--vcf $vcfFile \
--weir-fst-pop $groupFileDir/$group1.txt \
--weir-fst-pop $groupFileDir/$group2.txt \
--out pairwFst.$group1.$group2

echo "finished vcftools diversity analysis"
echo "end time `/bin/date`"