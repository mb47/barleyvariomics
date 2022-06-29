#$ -cwd
#$ -q ln.q
#$ -j yes
#$ -V
#$ -t 1-3

## I/O
bcfFile=variomeSNPs_final815samples_filteredSNPs.bcf

chrom[1]=chr1H
chrom[2]=chr2H
chrom[3]=chr3H
chrom[4]=chr4H
chrom[5]=chr5H
chrom[6]=chr6H
chrom[7]=chr7H

groupDir=groupList
group[1]=cultivars
group[2]=landraces
group[3]=spontaneums
groups=${group[SGE_TASK_ID]}


## Analysis
conda activate delAllel
for i in {4..5};do
	## Prepare vcf file of the designated chromosome
	bcftools view --min-af 0.01:minor --regions ${chrom[i]} $bcfFile > ${chrom[i]}.maf0.01.vcf
	## Calculate mu statistic
	RAiSD -n $groups.${chrom[i]}.maf0.01 -I ${chrom[i]}.maf0.01.vcf -S $groupDir/$groups.txt -M 1 -B 800000000 1 -f -R
	wait
	## Remove the vcf file
	rm ${chrom[i]}.vcf 
done
conda deactivate
