#$ -cwd
#$ -j yes
#$ -V
#$ -q ln.q
#$ -t 1-220
conda activate delAllel

accessionFile=815accessions.txt
refGenome=150831_barley_pseudomolecules.fasta
vcfFile=variomeSNPs_strictFiltersNoMAF_final.815lines.bcf
gffFile=Morex2017_BaRT.1_StringTie_Transdecoder_Homology.genome.gff3

for i in $(cat affectedGeneList.$SGE_TASK_ID.txt)
do
	# create temp gff file
	grep -w "$i" $gffFile | grep ".1.p1" > $i.temp.gff
	# create output directory
	mkdir $i.directory
	# reconstruct gene sequence
	./reconstruct.sh $i $accessionFile $refGenome $i.temp.gff $vcfFile $i.directory
	./transcriptFinder.sh $i $i.directory
	# remove temp gff file
	rm $i.temp.gff
done
