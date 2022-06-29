#$ -cwd
#$ -j yes
#$ -V
#$ -pe smp 32

conda activate delAllel

blastn \
	-db BaRT.1.0_.fa \
	-query WholeExome.fasta \
	-task megablast \
	-outfmt 6 \
	-evalue 1e-25 \
	-max_target_seqs 1 \
	-num_threads 32 \
	-out WholeExome2BaRT.1.0_.fa.outfmt6.evalue1e-25.txt

conda deactivate
