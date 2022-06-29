#$ -cwd
#$ -j yes
#$ -V
#$ -q ln.q
#$ -t 1-12

popSize[1]=3
popSize[2]=4
popSize[3]=5
popSize[4]=6
popSize[5]=7
popSize[6]=8
popSize[7]=9
popSize[8]=10
popSize[9]=11
popSize[10]=12
popSize[11]=13
popSize[12]=14
popSizes=${popSize[SGE_TASK_ID]}

input=random10k.plink
#mkdir structures

conda activate delAllel

structure.py -K $popSizes --input=$input --output=structures/output_simple.$popSizes

conda deactivate
