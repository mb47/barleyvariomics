#$ -cwd
#$ -j yes
#$ -V
#$ -q ln.q
# -pe smp 8
# -m bea

input=variomeSNPs.random4k.vcf
outName=random10k.plink
conda activate delAllel

plink --allow-extra-chr --vcf $input --out $outName

conda deactivate
