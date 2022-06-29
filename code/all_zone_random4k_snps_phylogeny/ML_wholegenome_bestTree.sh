#$ -cwd
#$ -j yes
#$ -V
#$ -q ln.q
#$ -pe smp 16

conda activate delAllel

raxml-ng --msa variomeSNPs.random4k.min4.phy \
--model GTR+G4 \
--prefix ./wholeGenome.bestTree.GTR+G4 \
--threads $NSLOTS
