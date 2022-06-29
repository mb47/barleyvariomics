#$ -cwd
#$ -j yes
#$ -V


conda activate delAllel

chooseK.py --input=structures/output_simple

conda deactivate
