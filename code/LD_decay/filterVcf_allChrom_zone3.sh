#$ -cwd
#$ -j yes
#$ -V
#$ -pe smp 8
#$ -t 1-7
#$ -q ln.q

bcfFile=jointGenotyperSNPs_allSamples_wholeChrCords_annotated_filteredSNPs.bcf
sampleFile=815accessions.txt

chrom[1]=chr1H
chrom[2]=chr2H
chrom[3]=chr3H
chrom[4]=chr4H
chrom[5]=chr5H
chrom[6]=chr6H
chrom[7]=chr7H
chroms=${chrom[SGE_TASK_ID]}
region[1]=chr1H:151200000-303200000
region[2]=chr2H:256000000-401600000
region[3]=chr3H:224000000-380000000
region[4]=chr4H:163200000-406500000
region[5]=chr5H:168800000-301600000
region[6]=chr6H:199200000-316000000
region[7]=chr7H:254400000-384000000
regions=${region[SGE_TASK_ID]}

conda activate delAllel

mkdir $chroms

bcftools view --threads $NSLOTS --regions $regions --samples-file $sampleFile --output-type v --output-file $chroms/$chroms.zone3.vcf $bcfFile
#bcftools stats --threads $NSLOTS chr4H.zone3.vcf

# convert to flapjack
java utils.snps.vcf.ConvertVCFToFlapjackFormat $chroms/$chroms.zone3.vcf $chroms/flapjack.$chroms.zone3.dat $chroms/flapjack.$chroms.zone3.map

conda deactivate
