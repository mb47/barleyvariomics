
#$ -pe smp 32
#$ -N GATK_JG 
#$ -q ln.q@n13-16-256-kif

refSeq=/mnt/shared/projects/barley/201409_MTP/pseudomolecules_Aug2015/final/pseudomolecules/150831_barley_pseudomolecules_parts.fasta

gatkJar=/mnt/shared/scratch/synology/mb40521/apps/GATK/3.4.0/GenomeAnalysisTK.jar

output=jointGenotyperSNPs_allSamples.vcf 

java -Xmx245g \
-jar $gatkJar \
-T GenotypeGVCFs \
-R $refSeq \
-o $output \
-nt $NSLOTS \
-V 1.g.vcf \
-V 2.g.vcf \
-V 3.g.vcf \
-V 4.g.vcf \
-V 5.g.vcf \
-V 6.g.vcf \
-V 7.g.vcf \
-V 8.g.vcf \
-V 9.g.vcf \
-V 10.g.vcf \
-V 11.g.vcf \
-V 12.g.vcf \
-V 13.g.vcf \
-V 14.g.vcf \
-V 15.g.vcf \
-V 16.g.vcf \
-V 17.g.vcf \
-V 18.g.vcf \
-V 19.g.vcf \
-V 20.g.vcf \
-V 21.g.vcf \
-V 22.g.vcf \
-V 23.g.vcf \
-V 24.g.vcf \
-V 25.g.vcf \
-V 26.g.vcf \
-V 27.g.vcf \
-V 28.g.vcf \
-V 29.g.vcf \
-V 30.g.vcf \
-V 31.g.vcf \
-V 32.g.vcf \
-V 33.g.vcf \
-V 34.g.vcf \
-V 35.g.vcf \
-V 36.g.vcf \
-V 37.g.vcf \
-V 38.g.vcf \
-V 39.g.vcf \
-V 40.g.vcf \
-V 41.g.vcf \
-V 42.g.vcf \
-V 43.g.vcf \
-V 44.g.vcf 

