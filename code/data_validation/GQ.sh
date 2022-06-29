#$ -cwd
#$ -j yes
#$ -V
#$ -q ln.q

## Prepare input file
grep -v "##" lowMAF.vcf | cut -f10-825 > lowMAF.table.tsv
grep -v "##" othersMarkers.vcf | cut -f10-825 > otherMarkers.table.tsv

## Read through lowMAF marker file and separate lowMAF individuals and other individuals
echo "TYPE	GQ" > lowMAF.toPlot.tsv
python sortGqByGeno.py lowMAF.table.tsv >> lowMAF.toPlot.tsv

## Read through other markers file
echo "TYPE      GQ" > otherMarkers.toPlot.tsv
python sortGqByGeno.py otherMarkers.table.tsv >> otherMarkers.toPlot.tsv

