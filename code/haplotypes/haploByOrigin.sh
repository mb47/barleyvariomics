#$ -cwd
#$ -j yes
#$ -V
#$ -q ln.q


for i in $(cat geneList2)
do
	./countHaplotype.sh $i.directory/$i.longest_orfs.p1.pep $i
done


