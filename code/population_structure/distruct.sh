#$ -cwd
#$ -j yes
#$ -V
#$ -t 4-8

conda activate delAllel

#mkdir graph

N=$SGE_TASK_ID

distruct.py -K $N --input=structures/output_simple.$N --popfile=popList4.txt --output=graph/K-$N.distruct.geo.png

conda deactivate
