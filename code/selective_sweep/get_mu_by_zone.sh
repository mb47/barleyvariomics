
#$ -cwd
#$ -j yes
#$ -V

chrom[1]=chr1H
chrom[2]=chr2H
chrom[3]=chr3H
chrom[4]=chr4H
chrom[5]=chr5H
chrom[6]=chr6H
chrom[7]=chr7H
group[1]=cultivars
group[2]=landraces
group[3]=spontaneums
zone[1]="zone_1"
zone[2]="zone_2"
zone[3]="zone_3"

# Iterate through populations
for i in {1..3}
do
	# Iterate through zones
	for j in {1..3}
		do
			input=mu_to_filter.txt
			grep ${zone[j]} coordinate.bed > temp_zone.bed
			bedtools intersect -a $input -b temp_zone.bed > ${group[i]}.${zone[j]}.mu.tsv
			awk -v OFS="\t" -v zone=${zone[j]} '{print $0,zone}' ${group[i]}.${zone[j]}.mu.tsv > ${group[i]}.${zone[j]}.final
			rm temp_zone.bed
			rm ${group[i]}.${zone[j]}.mu.tsv
		done
done

#echo -e zone"\t"blockLength"\t"group > final.tsv
cat *.final >> final.tsv

rm *.final
